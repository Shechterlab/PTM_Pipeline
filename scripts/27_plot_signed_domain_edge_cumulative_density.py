from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from common import BREWER_COLORS, apply_paper_style, save_figure, style_axis


def signed_distance(row: pd.Series) -> float | None:
    distance = row.get("nearest_domain_edge_distance")
    side = str(row.get("nearest_domain_side"))
    if pd.isna(distance) or side == "inside":
        return None
    distance = float(distance)
    if side == "N_terminal_to_domain":
        return -distance
    if side == "C_terminal_to_domain":
        return distance
    return None


def prepare_sites(site_context_path: Path, length_null_sites_path: Path) -> pd.DataFrame:
    length_sites = pd.read_csv(length_null_sites_path, sep="\t")
    context = pd.read_csv(site_context_path, sep="\t", low_memory=False)
    context = context[
        [
            "site_key",
            "canonical_UniProtAC",
            "corrected_position",
            "nearest_domain_edge_distance",
            "nearest_domain_side",
            "domain_context_class",
            "disorder_context_class",
        ]
    ].copy()
    context["position"] = pd.to_numeric(context["corrected_position"], errors="coerce").astype("Int64")
    context = context.drop(columns=["corrected_position"]).drop_duplicates(["canonical_UniProtAC", "position"])
    sites = length_sites[["canonical_UniProtAC", "position", "site_key", "outside_domain_idr_segment_id"]].merge(
        context,
        on=["canonical_UniProtAC", "position", "site_key"],
        how="left",
    )
    sites["signed_domain_edge_distance"] = sites.apply(signed_distance, axis=1)
    return sites[sites["signed_domain_edge_distance"].notna()].copy()


def annotate_domain_side(frame: pd.DataFrame, interpro_path: Path) -> pd.DataFrame:
    intervals = pd.read_csv(interpro_path, sep="\t", low_memory=False)
    intervals = intervals[intervals["interpro_type"].astype(str).isin({"domain", "repeat"})].copy()
    intervals["fragment_start"] = pd.to_numeric(intervals["fragment_start"], errors="coerce")
    intervals["fragment_end"] = pd.to_numeric(intervals["fragment_end"], errors="coerce")
    intervals = intervals[intervals["fragment_start"].notna() & intervals["fragment_end"].notna()].copy()
    intervals["fragment_start"] = intervals["fragment_start"].astype(int)
    intervals["fragment_end"] = intervals["fragment_end"].astype(int)
    intervals_by_acc = {
        accession: group.sort_values(["fragment_start", "fragment_end"]).reset_index(drop=True)
        for accession, group in intervals.groupby("canonical_UniProtAC")
    }
    out = frame.copy()
    out["nearest_domain_side"] = pd.NA
    for accession, group in out.groupby("canonical_UniProtAC", sort=False):
        domains = intervals_by_acc.get(accession)
        if domains is None or domains.empty:
            continue
        starts = domains["fragment_start"].to_numpy()
        ends = domains["fragment_end"].to_numpy()
        for idx, position in pd.to_numeric(group["position"], errors="coerce").items():
            if pd.isna(position):
                continue
            position = int(position)
            best_distance = None
            best_side = None
            for start, end in zip(starts, ends):
                if start <= position <= end:
                    distance = 0
                    side = "inside"
                elif position < start:
                    distance = start - position
                    side = "N_terminal_to_domain"
                else:
                    distance = position - end
                    side = "C_terminal_to_domain"
                if best_distance is None or distance < best_distance:
                    best_distance = distance
                    best_side = side
            out.at[idx, "nearest_domain_side"] = best_side
    return out


def prepare_background(background_path: Path, length_null_sites_path: Path, interpro_path: Path) -> pd.DataFrame:
    observed_segments = set(pd.read_csv(length_null_sites_path, sep="\t")["outside_domain_idr_segment_id"].dropna().astype(str))
    background = pd.read_csv(background_path, sep="\t", low_memory=False)
    # Reuse segment assignments from the length-null site table's intervals by assigning through the written segment file.
    # Background rows already have domain distance and side but not the outside-domain IDR segment ID.
    segment_sites = pd.read_csv(length_null_sites_path, sep="\t")
    segments = (
        segment_sites[
            ["canonical_UniProtAC", "outside_domain_idr_segment_id", "outside_domain_idr_start", "outside_domain_idr_end"]
        ]
        .drop_duplicates()
        .dropna()
    )
    frames = []
    for accession, group in background.groupby("canonical_UniProtAC", sort=False):
        segs = segments[segments["canonical_UniProtAC"] == accession]
        if segs.empty:
            continue
        positions = pd.to_numeric(group["position"], errors="coerce")
        for seg in segs.itertuples(index=False):
            if str(seg.outside_domain_idr_segment_id) not in observed_segments:
                continue
            mask = positions.between(int(seg.outside_domain_idr_start), int(seg.outside_domain_idr_end), inclusive="both")
            subset = group.loc[mask].copy()
            if subset.empty:
                continue
            subset["outside_domain_idr_segment_id"] = seg.outside_domain_idr_segment_id
            frames.append(subset)
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    out = out[
        out["domain_context_class"].isin(["boundary", "inter_domain_linker", "distal"])
        & out["nearest_domain_edge_distance"].notna()
        & out["outside_domain_idr_segment_id"].notna()
    ].copy()
    out = annotate_domain_side(out, interpro_path)
    out["signed_domain_edge_distance"] = out.apply(signed_distance, axis=1)
    return out[out["signed_domain_edge_distance"].notna()].copy()


def cumulative_side(frame: pd.DataFrame, side: str, thresholds: list[int]) -> pd.DataFrame:
    if side == "N":
        distances = -frame.loc[frame["signed_domain_edge_distance"] < 0, "signed_domain_edge_distance"].astype(float)
    else:
        distances = frame.loc[frame["signed_domain_edge_distance"] > 0, "signed_domain_edge_distance"].astype(float)
    total = len(distances)
    rows = []
    for threshold in thresholds:
        count = int((distances <= threshold).sum())
        rows.append(
            {
                "side": side,
                "distance_aa": threshold,
                "count_within_distance": count,
                "total_on_side": total,
                "fraction_within_distance": count / total if total else 0.0,
            }
        )
    return pd.DataFrame(rows)


def weighted_expected_side(observed: pd.DataFrame, universe: pd.DataFrame, side: str, thresholds: list[int]) -> pd.DataFrame:
    if side == "N":
        observed_side = observed[observed["signed_domain_edge_distance"] < 0].copy()
        universe_side = universe[universe["signed_domain_edge_distance"] < 0].copy()
        universe_side["abs_distance"] = -universe_side["signed_domain_edge_distance"].astype(float)
    else:
        observed_side = observed[observed["signed_domain_edge_distance"] > 0].copy()
        universe_side = universe[universe["signed_domain_edge_distance"] > 0].copy()
        universe_side["abs_distance"] = universe_side["signed_domain_edge_distance"].astype(float)

    segment_groups = {
        segment_id: group["abs_distance"].to_numpy(dtype=float)
        for segment_id, group in universe_side.groupby("outside_domain_idr_segment_id", sort=False)
    }
    observed_segment_ids = observed_side["outside_domain_idr_segment_id"].astype(str).tolist()
    rows = []
    for threshold in thresholds:
        expected = 0.0
        usable = 0
        for segment_id in observed_segment_ids:
            distances = segment_groups.get(segment_id)
            if distances is None or len(distances) == 0:
                continue
            expected += float((distances <= threshold).sum()) / float(len(distances))
            usable += 1
        rows.append(
            {
                "side": side,
                "distance_aa": threshold,
                "expected_count_within_distance": expected,
                "expected_total_site_weighted": usable,
                "expected_fraction_within_distance": expected / usable if usable else 0.0,
            }
        )
    return pd.DataFrame(rows)


def build_curve(observed: pd.DataFrame, background: pd.DataFrame, max_distance: int) -> pd.DataFrame:
    thresholds = list(range(1, max_distance + 1))
    obs = pd.concat([cumulative_side(observed, "N", thresholds), cumulative_side(observed, "C", thresholds)], ignore_index=True)
    universe_cols = ["outside_domain_idr_segment_id", "signed_domain_edge_distance"]
    universe = pd.concat([background[universe_cols], observed[universe_cols]], ignore_index=True).dropna()
    expected = pd.concat(
        [
            weighted_expected_side(observed, universe, "N", thresholds),
            weighted_expected_side(observed, universe, "C", thresholds),
        ],
        ignore_index=True,
    )
    out = obs.merge(expected, on=["side", "distance_aa"], how="left")
    out["obs_minus_expected_fraction"] = out["fraction_within_distance"] - out["expected_fraction_within_distance"]
    return out


def save_plot(fig, prefix: Path) -> None:
    save_figure(fig, prefix)
    fig.savefig(prefix.with_suffix(".svg"), bbox_inches="tight", facecolor="white")


def plot_curve(curve: pd.DataFrame, outdir: Path, max_distance: int) -> None:
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(7.7, 6.2), gridspec_kw={"height_ratios": [2.1, 1.0], "hspace": 0.34})
    palette = {"N": "#CBB2B9", "C": "#4D4D4D"}
    labels = {"N": "N-terminal side", "C": "C-terminal side"}
    for side in ["N", "C"]:
        side_df = curve[curve["side"] == side]
        ax.plot(side_df["distance_aa"], 100 * side_df["fraction_within_distance"], color=palette[side], linewidth=2.2, label=f"Observed, {labels[side]}")
        ax.plot(side_df["distance_aa"], 100 * side_df["expected_fraction_within_distance"], color=palette[side], linewidth=1.5, linestyle="--", alpha=0.75, label=f"Length null, {labels[side]}")
        ax2.plot(side_df["distance_aa"], 100 * side_df["obs_minus_expected_fraction"], color=palette[side], linewidth=2.0, label=labels[side])

    style_axis(ax)
    ax.set_ylabel("% within distance")
    ax.set_title("Cumulative methylarginine density from nearest structured-domain edge")
    ax.legend(frameon=False, loc="lower right", fontsize=8.2, ncol=2)
    ax.set_xlim(0, max_distance)
    ax.set_ylim(0, 100)

    style_axis(ax2, zero="y")
    ax2.axhline(0, color=BREWER_COLORS["dark_gray"], linewidth=0.8)
    ax2.set_xlabel("Absolute distance from nearest domain edge (aa)")
    ax2.set_ylabel("Observed - null (%)")
    ax2.set_xlim(0, max_distance)

    fig.text(
        0.99,
        0.01,
        "Site-weighted length null uses arginines from the same outside-domain IDR segment and same N/C side as each observed methyl-site.",
        ha="right",
        va="bottom",
        fontsize=7.8,
        color=BREWER_COLORS["dark_gray"],
    )
    fig.subplots_adjust(left=0.11, right=0.98, top=0.9, bottom=0.14)
    save_plot(fig, outdir / "arg_methyl_signed_domain_edge_cumulative_density_length_null")


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot bidirectional cumulative density from structured-domain edges.")
    parser.add_argument("--site-context", default="PTM_results/domain_context/sites_with_domain_context.tsv")
    parser.add_argument("--length-null-sites", default="PTM_results/domain_context/domain_edge_idr_length_null_sites.tsv")
    parser.add_argument("--background", default="PTM_results/domain_context/same_protein_nonmethyl_arginine_domain_background.tsv")
    parser.add_argument("--interpro", default="PTM_data/context/interpro_human_reviewed_domain_like_intervals.tsv")
    parser.add_argument("--outdir", default="PTM_results/domain_context")
    parser.add_argument("--max-distance", type=int, default=200)
    args = parser.parse_args()

    observed = prepare_sites(Path(args.site_context), Path(args.length_null_sites))
    background = prepare_background(Path(args.background), Path(args.length_null_sites), Path(args.interpro))
    curve = build_curve(observed, background, args.max_distance)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    curve.to_csv(outdir / "arg_methyl_signed_domain_edge_cumulative_density_length_null.tsv", sep="\t", index=False)
    qc = pd.DataFrame(
        [
            {
                "observed_sites": len(observed),
                "background_arginines": len(background),
                "observed_n_terminal_sites": int((observed["signed_domain_edge_distance"] < 0).sum()),
                "observed_c_terminal_sites": int((observed["signed_domain_edge_distance"] > 0).sum()),
                "background_nonmethyl_n_terminal_arginines": int((background["signed_domain_edge_distance"] < 0).sum()),
                "background_nonmethyl_c_terminal_arginines": int((background["signed_domain_edge_distance"] > 0).sum()),
            }
        ]
    )
    qc.to_csv(outdir / "arg_methyl_signed_domain_edge_cumulative_density_length_null_qc.tsv", sep="\t", index=False)
    apply_paper_style()
    plot_curve(curve, outdir, args.max_distance)
    print(outdir / "arg_methyl_signed_domain_edge_cumulative_density_length_null.svg")
    print(outdir / "arg_methyl_signed_domain_edge_cumulative_density_length_null.tsv")


if __name__ == "__main__":
    main()

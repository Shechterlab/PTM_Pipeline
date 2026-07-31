from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import BREWER_COLORS, add_q_values, apply_paper_style, canonical_site_table, parse_fasta, save_figure, style_axis


def residue_positions(sequence: str, residue: str = "R") -> np.ndarray:
    return np.array([idx + 1 for idx, aa in enumerate(sequence) if aa == residue], dtype=int)


def validate_sites(sites: pd.DataFrame, sequences: dict[str, str]) -> pd.DataFrame:
    out = canonical_site_table(sites, accession_col="canonical_UniProtAC", position_col="corrected_position", residue_col="residue")
    out = out[out["residue"] == "R"].copy()
    out = out[out["canonical_UniProtAC"].isin(sequences)].copy()
    valid = []
    for row in out.itertuples(index=False):
        seq = sequences.get(str(row.canonical_UniProtAC), "")
        pos = int(row.corrected_position)
        valid.append(bool(seq) and 1 <= pos <= len(seq) and seq[pos - 1] == "R")
    out = out[valid].drop_duplicates(["canonical_UniProtAC", "corrected_position"]).copy()
    return out


def neighbor_counts_for_thresholds(positions: np.ndarray, thresholds: list[int]) -> np.ndarray:
    positions = np.sort(np.asarray(positions, dtype=int))
    n = len(positions)
    counts = np.zeros((n, len(thresholds)), dtype=np.int16)
    if n <= 1:
        return counts
    for t_idx, threshold in enumerate(thresholds):
        left = 0
        right = 0
        for idx, pos in enumerate(positions):
            while positions[left] < pos - threshold:
                left += 1
            while right + 1 < n and positions[right + 1] <= pos + threshold:
                right += 1
            counts[idx, t_idx] = right - left
    return counts


def occupancy_surface(positions_by_acc: dict[str, np.ndarray], thresholds: list[int], max_occupancy: int) -> pd.DataFrame:
    total_sites = sum(len(positions) for positions in positions_by_acc.values())
    hits = np.zeros((max_occupancy, len(thresholds)), dtype=float)
    if total_sites == 0:
        return pd.DataFrame()
    for positions in positions_by_acc.values():
        counts = neighbor_counts_for_thresholds(positions, thresholds)
        if counts.size == 0:
            continue
        for occ in range(1, max_occupancy + 1):
            hits[occ - 1, :] += (counts >= occ).sum(axis=0)
    rows = []
    for occ in range(1, max_occupancy + 1):
        for t_idx, threshold in enumerate(thresholds):
            rows.append(
                {
                    "occupancy_k": occ,
                    "distance_aa": threshold,
                    "sites_with_at_least_k_neighbors": int(hits[occ - 1, t_idx]),
                    "sites_total": int(total_sites),
                    "observed_fraction": hits[occ - 1, t_idx] / total_sites,
                }
            )
    return pd.DataFrame(rows)


def null_surfaces(
    arg_positions_by_acc: dict[str, np.ndarray],
    methyl_counts_by_acc: dict[str, int],
    thresholds: list[int],
    max_occupancy: int,
    permutations: int,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    accessions = [acc for acc, count in methyl_counts_by_acc.items() if count > 0 and len(arg_positions_by_acc.get(acc, [])) >= count]
    for permutation in range(1, permutations + 1):
        sampled_by_acc = {}
        for acc in accessions:
            candidates = arg_positions_by_acc[acc]
            sampled_by_acc[acc] = np.sort(rng.choice(candidates, size=methyl_counts_by_acc[acc], replace=False))
        surface = occupancy_surface(sampled_by_acc, thresholds, max_occupancy)
        surface["permutation"] = permutation
        rows.append(surface)
    return pd.concat(rows, ignore_index=True)


def summarize_with_null(observed: pd.DataFrame, nulls: pd.DataFrame) -> pd.DataFrame:
    null_summary = (
        nulls.groupby(["occupancy_k", "distance_aa"], as_index=False)
        .agg(
            null_mean_fraction=("observed_fraction", "mean"),
            null_q025_fraction=("observed_fraction", lambda x: float(pd.Series(x).quantile(0.025))),
            null_q975_fraction=("observed_fraction", lambda x: float(pd.Series(x).quantile(0.975))),
        )
    )
    out = observed.merge(null_summary, on=["occupancy_k", "distance_aa"], how="left")
    p_values = []
    for row in out.itertuples(index=False):
        null_group = nulls[(nulls["occupancy_k"] == row.occupancy_k) & (nulls["distance_aa"] == row.distance_aa)]
        null_values = null_group["observed_fraction"].to_numpy(dtype=float)
        p_values.append((np.count_nonzero(null_values >= row.observed_fraction) + 1) / (len(null_values) + 1))
    out["empirical_p_greater"] = p_values
    ratio = out["observed_fraction"] / out["null_mean_fraction"].replace(0, np.nan)
    out["log2_observed_over_null"] = pd.Series(np.nan, index=out.index, dtype=float)
    valid = ratio.notna() & (ratio > 0)
    out.loc[valid, "log2_observed_over_null"] = np.log2(ratio.loc[valid])
    out["delta_fraction_vs_null"] = out["observed_fraction"] - out["null_mean_fraction"]
    return add_q_values(out, p_col="empirical_p_greater", out_col="q_value")


def save_plot(fig, prefix: Path) -> None:
    save_figure(fig, prefix)
    fig.savefig(prefix.with_suffix(".svg"), bbox_inches="tight", facecolor="white")


def plot_heatmaps(summary: pd.DataFrame, outdir: Path, q_threshold: float) -> None:
    distance_order = sorted(summary["distance_aa"].unique())
    occupancy_order = sorted(summary["occupancy_k"].unique(), reverse=True)

    log2_matrix = (
        summary.pivot_table(index="occupancy_k", columns="distance_aa", values="log2_observed_over_null", observed=True)
        .reindex(index=occupancy_order, columns=distance_order)
    )
    obs_matrix = (
        summary.pivot_table(index="occupancy_k", columns="distance_aa", values="observed_fraction", observed=True)
        .reindex(index=occupancy_order, columns=distance_order)
    )
    q_matrix = (
        summary.pivot_table(index="occupancy_k", columns="distance_aa", values="q_value", observed=True)
        .reindex(index=occupancy_order, columns=distance_order)
    )

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.4, 5.3), gridspec_kw={"width_ratios": [1.05, 1.0], "wspace": 0.25})
    values = log2_matrix.to_numpy(dtype=float)
    vmax = max(0.5, float(np.nanmax(np.abs(values))) if np.isfinite(values).any() else 0.5)
    im = ax.imshow(values, aspect="auto", cmap="RdGy_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(distance_order)), [str(x) for x in distance_order], rotation=35)
    ax.set_yticks(range(len(occupancy_order)), [f">={x}" for x in occupancy_order])
    ax.set_xlabel("Distance cutoff X (aa)")
    ax.set_ylabel("Other methyl-Arg sites within X aa")
    ax.set_title("Pairwise occupancy enrichment")
    for y in range(q_matrix.shape[0]):
        for x in range(q_matrix.shape[1]):
            q_value = q_matrix.iloc[y, x]
            if pd.notna(q_value) and q_value <= q_threshold:
                ax.text(x, y, "*", ha="center", va="center", fontsize=8.5, color="#111111")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label("log2(observed / arginine null)")

    obs_values = 100.0 * obs_matrix.to_numpy(dtype=float)
    im2 = ax2.imshow(obs_values, aspect="auto", cmap="Greys", vmin=0, vmax=max(5.0, float(np.nanmax(obs_values))))
    ax2.set_xticks(range(len(distance_order)), [str(x) for x in distance_order], rotation=35)
    ax2.set_yticks(range(len(occupancy_order)), [f">={x}" for x in occupancy_order])
    ax2.set_xlabel("Distance cutoff X (aa)")
    ax2.set_ylabel("")
    ax2.set_title("Observed cumulative occupancy")
    cbar2 = fig.colorbar(im2, ax=ax2, fraction=0.046, pad=0.02)
    cbar2.set_label("% of methyl-sites")

    fig.suptitle("Methylarginine pairwise occupancy and cumulative local density", y=0.98, fontsize=13)
    fig.text(
        0.99,
        0.01,
        f"Asterisks mark BH q<={q_threshold}. Null samples methyl-site counts from arginines in the same proteins.",
        ha="right",
        va="bottom",
        fontsize=7.8,
        color=BREWER_COLORS["dark_gray"],
    )
    fig.subplots_adjust(left=0.08, right=0.96, top=0.84, bottom=0.17)
    save_plot(fig, outdir / "arg_methyl_pairwise_occupancy_heatmap")


def main() -> None:
    parser = argparse.ArgumentParser(description="Create 2D pairwise methylarginine occupancy heatmaps.")
    parser.add_argument("--sites", default="PTM_results/domain_context/sites_with_domain_context.tsv")
    parser.add_argument("--fasta", default="PTM_data/context/uniprot_human_reviewed_canonical.fasta")
    parser.add_argument("--outdir", default="PTM_results/methyl_arg_clustering")
    parser.add_argument("--distances", nargs="*", type=int, default=[1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40, 50, 75, 100])
    parser.add_argument("--max-occupancy", type=int, default=10)
    parser.add_argument("--permutations", type=int, default=250)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--q-threshold", type=float, default=0.05)
    args = parser.parse_args()

    sites = pd.read_csv(args.sites, sep="\t", low_memory=False)
    sequences = parse_fasta(args.fasta)
    sites = validate_sites(sites, sequences)

    observed_by_acc = {
        accession: np.sort(group["corrected_position"].to_numpy(dtype=int))
        for accession, group in sites.groupby("canonical_UniProtAC", sort=False)
    }
    methyl_counts_by_acc = {acc: len(pos) for acc, pos in observed_by_acc.items()}
    arg_positions_by_acc = {acc: residue_positions(sequences[acc], "R") for acc in observed_by_acc if acc in sequences}

    distances = sorted(set(args.distances))
    observed = occupancy_surface(observed_by_acc, distances, args.max_occupancy)
    nulls = null_surfaces(arg_positions_by_acc, methyl_counts_by_acc, distances, args.max_occupancy, args.permutations, args.seed)
    summary = summarize_with_null(observed, nulls)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(outdir / "arg_methyl_pairwise_occupancy_surface.tsv", sep="\t", index=False)
    nulls.to_csv(outdir / "arg_methyl_pairwise_occupancy_null_permutations.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "sites_analyzed": int(sum(len(pos) for pos in observed_by_acc.values())),
                "proteins_analyzed": int(len(observed_by_acc)),
                "permutations": int(args.permutations),
                "distances": ",".join(map(str, distances)),
                "max_occupancy": int(args.max_occupancy),
            }
        ]
    ).to_csv(outdir / "arg_methyl_pairwise_occupancy_qc.tsv", sep="\t", index=False)

    apply_paper_style()
    plot_heatmaps(summary, outdir, args.q_threshold)
    print(outdir / "arg_methyl_pairwise_occupancy_heatmap.svg")
    print(outdir / "arg_methyl_pairwise_occupancy_surface.tsv")


if __name__ == "__main__":
    main()

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from common import BREWER_COLORS, apply_paper_style, save_figure, style_axis


ARCHITECTURE_ORDER = [
    "broad across IDR",
    "domain-proximal clustered",
    "domain-distal clustered",
]

ARCHITECTURE_COLORS = {
    "broad across IDR": "#CBB2B9",
    "domain-proximal clustered": BREWER_COLORS["mid_gray"],
    "domain-distal clustered": BREWER_COLORS["dark_gray"],
}

SIGNED_BINS = [-10**9, -200, -100, -40, -20, 0, 20, 40, 100, 200, 10**9]
SIGNED_LABELS = [
    "N >200",
    "N 101-200",
    "N 41-100",
    "N 21-40",
    "N <=20",
    "C <=20",
    "C 21-40",
    "C 41-100",
    "C 101-200",
    "C >200",
]


def parse_positions(value: str) -> list[int]:
    out: list[int] = []
    for token in str(value).replace(";", ",").split(","):
        token = token.strip()
        if token:
            out.append(int(token))
    return out


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


def side_label(side: str) -> str:
    if side == "N_terminal_to_domain":
        return "N-terminal side of domain"
    if side == "C_terminal_to_domain":
        return "C-terminal side of domain"
    return str(side)


def build_cluster_site_table(cluster_path: Path, profiles_path: Path, site_context_path: Path, cluster_threshold: int) -> pd.DataFrame:
    clusters = pd.read_csv(cluster_path, sep="\t")
    clusters = clusters[clusters["max_sites_in_20aa_span_within_bin"] >= cluster_threshold].copy()
    profiles = pd.read_csv(profiles_path, sep="\t")[
        ["canonical_UniProtAC", "gene", "architecture_class"]
    ].drop_duplicates("canonical_UniProtAC")
    site_context = pd.read_csv(site_context_path, sep="\t", low_memory=False)
    site_context = site_context[
        [
            "canonical_UniProtAC",
            "corrected_position",
            "site_key",
            "nearest_domain_edge_distance",
            "nearest_domain_side",
            "nearest_domain_class",
            "nearest_interpro_name",
        ]
    ].copy()
    site_context["position"] = pd.to_numeric(site_context["corrected_position"], errors="coerce").astype("Int64")
    site_context = site_context.drop(columns=["corrected_position"]).drop_duplicates(["canonical_UniProtAC", "position"])

    rows = []
    for cluster in clusters.itertuples(index=False):
        for position in parse_positions(cluster.best_window_positions):
            rows.append(
                {
                    "canonical_UniProtAC": cluster.canonical_UniProtAC,
                    "gene": cluster.gene,
                    "domain_edge_bin_compact": cluster.domain_edge_bin_compact,
                    "cluster_window_start": int(cluster.best_window_start),
                    "cluster_window_end": int(cluster.best_window_end),
                    "max_sites_in_20aa_span_within_bin": int(cluster.max_sites_in_20aa_span_within_bin),
                    "position": int(position),
                }
            )
    cluster_sites = pd.DataFrame(rows).drop_duplicates(
        ["canonical_UniProtAC", "position", "domain_edge_bin_compact", "cluster_window_start", "cluster_window_end"]
    )
    cluster_sites = cluster_sites.merge(profiles, on=["canonical_UniProtAC", "gene"], how="left")
    cluster_sites = cluster_sites.merge(site_context, on=["canonical_UniProtAC", "position"], how="left")
    cluster_sites["signed_domain_edge_distance"] = cluster_sites.apply(signed_distance, axis=1)
    cluster_sites = cluster_sites[cluster_sites["signed_domain_edge_distance"].notna()].copy()
    cluster_sites["domain_side"] = cluster_sites["nearest_domain_side"].map(side_label)
    cluster_sites["signed_distance_bin"] = pd.cut(
        cluster_sites["signed_domain_edge_distance"],
        bins=SIGNED_BINS,
        labels=SIGNED_LABELS,
        include_lowest=True,
        right=True,
    )
    cluster_sites["abs_distance_class"] = cluster_sites["nearest_domain_edge_distance"].map(
        lambda value: "<=20 aa" if float(value) <= 20 else ">20 aa"
    )
    return cluster_sites


def summarize(cluster_sites: pd.DataFrame, profiles: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    architecture_counts = (
        profiles.groupby("architecture_class", observed=True)
        .agg(n_proteins=("canonical_UniProtAC", "nunique"))
        .reindex(ARCHITECTURE_ORDER)
        .reset_index()
    )
    side_summary = (
        cluster_sites.groupby(["architecture_class", "domain_side"], observed=True)
        .agg(
            clustered_window_sites=("site_key", "nunique"),
            proteins=("canonical_UniProtAC", "nunique"),
            median_abs_distance=("nearest_domain_edge_distance", "median"),
        )
        .reset_index()
        .sort_values(["architecture_class", "domain_side"])
    )
    distance_summary = (
        cluster_sites.groupby(["architecture_class", "signed_distance_bin"], observed=True)
        .agg(clustered_window_sites=("site_key", "nunique"), proteins=("canonical_UniProtAC", "nunique"))
        .reset_index()
    )
    return architecture_counts, side_summary, distance_summary


def save_plot(fig, prefix: Path) -> None:
    save_figure(fig, prefix)
    fig.savefig(prefix.with_suffix(".svg"), bbox_inches="tight", facecolor="white")


def plot_histogram(distance_summary: pd.DataFrame, outdir: Path, cluster_threshold: int) -> None:
    pivot = distance_summary.pivot_table(
        index="signed_distance_bin",
        columns="architecture_class",
        values="clustered_window_sites",
        fill_value=0,
        observed=True,
    ).reindex(SIGNED_LABELS, fill_value=0)
    pivot = pivot.reindex(columns=ARCHITECTURE_ORDER, fill_value=0)

    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    bottom = pd.Series(0, index=pivot.index, dtype=float)
    for cls in ARCHITECTURE_ORDER:
        ax.bar(
            pivot.index,
            pivot[cls],
            bottom=bottom,
            color=ARCHITECTURE_COLORS[cls],
            edgecolor="white",
            linewidth=0.8,
            label=cls,
        )
        bottom += pivot[cls]
    style_axis(ax)
    ax.set_ylabel("Cluster-window methyl-sites")
    ax.set_xlabel("Signed distance from nearest domain edge")
    ax.set_title(f"Clustered methylarginine sites map to both domain flanks ({cluster_threshold}+ sites / 20 aa)")
    ax.tick_params(axis="x", rotation=35)
    ax.legend(frameon=False, loc="upper right", fontsize=8.3)
    ax.axvline(4.5, color=BREWER_COLORS["dark_gray"], linewidth=0.8)
    ax.text(2.0, ax.get_ylim()[1] * 0.94, "N-terminal side", ha="center", va="top", fontsize=9, color=BREWER_COLORS["dark_gray"])
    ax.text(7.0, ax.get_ylim()[1] * 0.94, "C-terminal side", ha="center", va="top", fontsize=9, color=BREWER_COLORS["dark_gray"])
    fig.text(
        0.99,
        0.01,
        "Negative bins are N-terminal to the nearest domain; positive bins are C-terminal to the nearest domain.",
        ha="right",
        va="bottom",
        fontsize=7.8,
        color=BREWER_COLORS["dark_gray"],
    )
    fig.subplots_adjust(left=0.09, right=0.98, top=0.88, bottom=0.28)
    save_plot(fig, outdir / f"domain_cluster_signed_distance_histogram_ge{cluster_threshold}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot signed domain-edge distances for clustered methylarginine sites.")
    parser.add_argument("--clusters", default="PTM_results/domain_context/domain_edge_class_20aa_cluster_examples.tsv")
    parser.add_argument("--profiles", default="")
    parser.add_argument("--site-context", default="PTM_results/domain_context/sites_with_domain_context.tsv")
    parser.add_argument("--outdir", default="PTM_results/domain_context")
    parser.add_argument("--cluster-threshold", type=int, default=3)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    profiles_path = Path(args.profiles) if args.profiles else outdir / f"domain_cluster_architecture_classes_ge{args.cluster_threshold}.tsv"
    cluster_sites = build_cluster_site_table(Path(args.clusters), profiles_path, Path(args.site_context), args.cluster_threshold)
    profiles = pd.read_csv(profiles_path, sep="\t")
    architecture_counts, side_summary, distance_summary = summarize(cluster_sites, profiles)

    suffix = f"ge{args.cluster_threshold}"
    cluster_sites.to_csv(outdir / f"domain_cluster_signed_distance_sites_{suffix}.tsv", sep="\t", index=False)
    architecture_counts.to_csv(outdir / f"domain_cluster_signed_distance_architecture_counts_{suffix}.tsv", sep="\t", index=False)
    side_summary.to_csv(outdir / f"domain_cluster_signed_distance_side_summary_{suffix}.tsv", sep="\t", index=False)
    distance_summary.to_csv(outdir / f"domain_cluster_signed_distance_histogram_counts_{suffix}.tsv", sep="\t", index=False)

    apply_paper_style()
    plot_histogram(distance_summary, outdir, args.cluster_threshold)
    print(outdir / f"domain_cluster_signed_distance_side_summary_{suffix}.tsv")
    print(outdir / f"domain_cluster_signed_distance_histogram_ge{args.cluster_threshold}.svg")


if __name__ == "__main__":
    main()

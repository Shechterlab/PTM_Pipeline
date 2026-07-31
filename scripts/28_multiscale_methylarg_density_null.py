from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import BREWER_COLORS, add_q_values, apply_paper_style, parse_fasta, save_figure, style_axis


def residue_positions(sequence: str, residue: str = "R") -> np.ndarray:
    return np.array([idx + 1 for idx, aa in enumerate(sequence) if aa == residue], dtype=int)


def max_sites_in_window(positions: np.ndarray, window: int) -> tuple[int, int, int]:
    positions = np.sort(np.asarray(positions, dtype=int))
    if len(positions) == 0:
        return 0, 0, 0
    best = 1
    best_left = int(positions[0])
    best_right = int(positions[0])
    left = 0
    for right in range(len(positions)):
        while positions[right] - positions[left] > window:
            left += 1
        count = right - left + 1
        if count > best:
            best = count
            best_left = int(positions[left])
            best_right = int(positions[right])
    return best, best_left, best_right


def window_site_count(positions: np.ndarray, start: int, end: int) -> int:
    positions = np.asarray(positions, dtype=int)
    return int(((positions >= start) & (positions <= end)).sum())


def empirical_scan(
    methyl_positions: np.ndarray,
    arg_positions: np.ndarray,
    windows: list[int],
    permutations: int,
    rng: np.random.Generator,
) -> list[dict]:
    methyl_positions = np.sort(np.asarray(methyl_positions, dtype=int))
    arg_positions = np.sort(np.asarray(arg_positions, dtype=int))
    methyl_count = len(methyl_positions)
    arg_count = len(arg_positions)
    if methyl_count == 0 or arg_count == 0 or methyl_count > arg_count:
        return []

    observed = {}
    null_values = {window: [] for window in windows}
    for window in windows:
        max_count, start, end = max_sites_in_window(methyl_positions, window)
        observed[window] = (max_count, start, end)

    for _ in range(permutations):
        sampled = np.sort(rng.choice(arg_positions, size=methyl_count, replace=False))
        for window in windows:
            null_values[window].append(max_sites_in_window(sampled, window)[0])

    rows = []
    for window in windows:
        obs_count, start, end = observed[window]
        null = np.asarray(null_values[window], dtype=float)
        null_mean = float(null.mean()) if len(null) else math.nan
        p_greater = (int((null >= obs_count).sum()) + 1) / (len(null) + 1) if len(null) else math.nan
        arg_count_in_window = window_site_count(arg_positions, start, end)
        methyl_positions_in_window = methyl_positions[(methyl_positions >= start) & (methyl_positions <= end)]
        rows.append(
            {
                "window_aa": int(window),
                "methyl_site_count": int(methyl_count),
                "protein_arg_count": int(arg_count),
                "observed_max_methyl_sites": int(obs_count),
                "null_mean_max_sites": null_mean,
                "null_q025_max_sites": float(np.quantile(null, 0.025)) if len(null) else math.nan,
                "null_q975_max_sites": float(np.quantile(null, 0.975)) if len(null) else math.nan,
                "observed_over_null_mean": obs_count / null_mean if null_mean > 0 else math.nan,
                "log2_observed_over_null_mean": math.log2(obs_count / null_mean) if null_mean > 0 and obs_count > 0 else math.nan,
                "empirical_p_greater": p_greater,
                "best_window_start": int(start),
                "best_window_end": int(end),
                "best_window_arg_count": int(arg_count_in_window),
                "best_window_methyl_fraction_of_args": obs_count / arg_count_in_window if arg_count_in_window else math.nan,
                "best_window_methyl_positions": ",".join(map(str, methyl_positions_in_window.tolist())),
            }
        )
    return rows


def plot_heatmap(scan: pd.DataFrame, outdir: Path, top_n: int, q_threshold: float) -> None:
    sig = scan[scan["q_value"] <= q_threshold].copy()
    if sig.empty:
        plot_df = scan.sort_values(["max_log2_enrichment_for_protein", "methyl_site_count"], ascending=False).head(top_n)
    else:
        keep = (
            sig.groupby("canonical_UniProtAC")
            .agg(max_log2=("log2_observed_over_null_mean", "max"), methyl_site_count=("methyl_site_count", "max"))
            .sort_values(["max_log2", "methyl_site_count"], ascending=False)
            .head(top_n)
            .index
        )
        plot_df = scan[scan["canonical_UniProtAC"].isin(keep)].copy()

    protein_order = (
        plot_df.groupby(["canonical_UniProtAC", "gene"], dropna=False)
        .agg(max_log2=("log2_observed_over_null_mean", "max"), methyl_site_count=("methyl_site_count", "max"))
        .sort_values(["max_log2", "methyl_site_count"], ascending=True)
        .reset_index()
    )
    labels = [
        str(row.gene) if pd.notna(row.gene) and str(row.gene) else str(row.canonical_UniProtAC)
        for row in protein_order.itertuples(index=False)
    ]
    matrix = (
        plot_df.pivot_table(
            index="canonical_UniProtAC",
            columns="window_aa",
            values="log2_observed_over_null_mean",
            fill_value=np.nan,
            observed=True,
        )
        .reindex(protein_order["canonical_UniProtAC"])
        .sort_index(axis=1)
    )
    q_matrix = (
        plot_df.pivot_table(index="canonical_UniProtAC", columns="window_aa", values="q_value", fill_value=np.nan, observed=True)
        .reindex(matrix.index)
        .reindex(columns=matrix.columns)
    )

    fig_h = max(4.6, 0.28 * len(matrix) + 1.3)
    fig, ax = plt.subplots(figsize=(7.6, fig_h))
    values = matrix.to_numpy(dtype=float)
    vmax = max(0.7, float(np.nanmax(np.abs(values))) if np.isfinite(values).any() else 0.7)
    im = ax.imshow(values, aspect="auto", cmap="RdGy_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(matrix.columns)), [str(int(col)) for col in matrix.columns])
    ax.set_yticks(range(len(matrix.index)), labels)
    ax.set_xlabel("Rolling window size (aa)")
    ax.set_ylabel("Protein")
    ax.set_title("Multiscale methylarginine cluster density vs arginine null")
    for y in range(q_matrix.shape[0]):
        for x in range(q_matrix.shape[1]):
            q_value = q_matrix.iloc[y, x]
            if pd.notna(q_value) and q_value <= q_threshold:
                ax.text(x, y, "*", ha="center", va="center", color="#111111", fontsize=9.5)
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("log2(observed max density / arginine-null mean)")
    fig.text(
        0.99,
        0.01,
        f"Asterisks mark BH q<={q_threshold}. Null samples the same number of methyl-sites from arginines in the same protein.",
        ha="right",
        va="bottom",
        fontsize=7.8,
        color=BREWER_COLORS["dark_gray"],
    )
    fig.subplots_adjust(left=0.19, right=0.92, top=0.92, bottom=0.1)
    save_figure(fig, outdir / "arg_methyl_multiscale_density_heatmap")
    fig.savefig(outdir / "arg_methyl_multiscale_density_heatmap.svg", bbox_inches="tight", facecolor="white")


def plot_window_summary(scan: pd.DataFrame, outdir: Path, q_threshold: float) -> None:
    summary = (
        scan.assign(significant=scan["q_value"] <= q_threshold)
        .groupby("window_aa")
        .agg(
            significant_proteins=("significant", "sum"),
            tested_proteins=("canonical_UniProtAC", "nunique"),
            median_log2_enrichment=("log2_observed_over_null_mean", "median"),
            q10_log2_enrichment=("log2_observed_over_null_mean", lambda x: float(pd.Series(x).quantile(0.1))),
            q90_log2_enrichment=("log2_observed_over_null_mean", lambda x: float(pd.Series(x).quantile(0.9))),
        )
        .reset_index()
    )
    summary.to_csv(outdir / "arg_methyl_multiscale_density_window_summary.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(summary["window_aa"], summary["significant_proteins"], color="#CBB2B9", linewidth=2.2, marker="o")
    style_axis(ax, grid_axis="y")
    ax.set_xlabel("Rolling window size (aa)")
    ax.set_ylabel(f"Proteins with q<={q_threshold}")
    ax.set_title("Significant methylarginine clustering across window sizes")
    for row in summary.itertuples(index=False):
        ax.text(row.window_aa, row.significant_proteins + 1.0, str(int(row.significant_proteins)), ha="center", va="bottom", fontsize=8.2)
    fig.subplots_adjust(left=0.12, right=0.98, top=0.87, bottom=0.16)
    save_figure(fig, outdir / "arg_methyl_multiscale_density_significant_counts")
    fig.savefig(outdir / "arg_methyl_multiscale_density_significant_counts.svg", bbox_inches="tight", facecolor="white")


def main() -> None:
    parser = argparse.ArgumentParser(description="Multiscale methylarginine density scan against a same-protein arginine-position null.")
    parser.add_argument("--sites", default="PTM_results/domain_context/sites_with_domain_context.tsv")
    parser.add_argument("--fasta", default="PTM_data/context/uniprot_human_reviewed_canonical.fasta")
    parser.add_argument("--outdir", default="PTM_results/methyl_arg_clustering")
    parser.add_argument("--windows", nargs="*", type=int, default=[5, 10, 15, 20, 30, 50, 75, 100])
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--min-sites", type=int, default=3)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--top-n", type=int, default=35)
    parser.add_argument("--q-threshold", type=float, default=0.05)
    args = parser.parse_args()

    sites = pd.read_csv(args.sites, sep="\t", low_memory=False)
    sequences = parse_fasta(args.fasta)
    sites = sites[sites["canonical_UniProtAC"].isin(sequences)].copy()
    sites["corrected_position"] = pd.to_numeric(sites["corrected_position"], errors="coerce")
    sites = sites[sites["corrected_position"].notna()].copy()
    sites["corrected_position"] = sites["corrected_position"].astype(int)
    sites = sites.drop_duplicates(["canonical_UniProtAC", "corrected_position"])

    rng = np.random.default_rng(args.seed)
    rows = []
    skipped = 0
    for accession, group in sites.groupby("canonical_UniProtAC", sort=False):
        sequence = sequences.get(str(accession), "")
        arg_positions = residue_positions(sequence, "R")
        methyl_positions = np.sort(group["corrected_position"].to_numpy(dtype=int))
        methyl_positions = methyl_positions[np.isin(methyl_positions, arg_positions)]
        if len(methyl_positions) < args.min_sites:
            skipped += 1
            continue
        gene_values = group["substrate_genename"].dropna().astype(str)
        gene = gene_values.iloc[0] if not gene_values.empty else ""
        for row in empirical_scan(methyl_positions, arg_positions, sorted(set(args.windows)), args.permutations, rng):
            row["canonical_UniProtAC"] = accession
            row["gene"] = gene
            rows.append(row)

    scan = pd.DataFrame(rows)
    scan = add_q_values(scan, p_col="empirical_p_greater", out_col="q_value")
    max_by_protein = scan.groupby("canonical_UniProtAC")["log2_observed_over_null_mean"].max().rename("max_log2_enrichment_for_protein")
    scan = scan.merge(max_by_protein, on="canonical_UniProtAC", how="left")
    scan = scan.sort_values(["q_value", "log2_observed_over_null_mean", "methyl_site_count"], ascending=[True, False, False])

    significant_regions = scan[scan["q_value"] <= args.q_threshold].copy()
    best_regions = (
        scan.sort_values(["q_value", "log2_observed_over_null_mean"], ascending=[True, False])
        .drop_duplicates("canonical_UniProtAC")
        .sort_values(["q_value", "log2_observed_over_null_mean"], ascending=[True, False])
    )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    scan.to_csv(outdir / "arg_methyl_multiscale_density_scan.tsv", sep="\t", index=False)
    significant_regions.to_csv(outdir / "arg_methyl_multiscale_density_significant_regions.tsv", sep="\t", index=False)
    best_regions.to_csv(outdir / "arg_methyl_multiscale_density_best_region_per_protein.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "proteins_tested": int(scan["canonical_UniProtAC"].nunique()) if not scan.empty else 0,
                "proteins_skipped_below_min_sites": int(skipped),
                "protein_window_tests": int(len(scan)),
                "permutations": int(args.permutations),
                "min_sites": int(args.min_sites),
                "windows": ",".join(map(str, sorted(set(args.windows)))),
                "significant_protein_window_tests": int(len(significant_regions)),
                "significant_proteins": int(significant_regions["canonical_UniProtAC"].nunique()) if not significant_regions.empty else 0,
            }
        ]
    ).to_csv(outdir / "arg_methyl_multiscale_density_qc.tsv", sep="\t", index=False)

    apply_paper_style()
    if not scan.empty:
        plot_heatmap(scan, outdir, top_n=args.top_n, q_threshold=args.q_threshold)
        plot_window_summary(scan, outdir, q_threshold=args.q_threshold)

    print(outdir / "arg_methyl_multiscale_density_scan.tsv")
    print(outdir / "arg_methyl_multiscale_density_heatmap.svg")


if __name__ == "__main__":
    main()

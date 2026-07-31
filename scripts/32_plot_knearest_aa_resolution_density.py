from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from scipy.ndimage import gaussian_filter

from common import BREWER_COLORS, add_q_values, canonical_site_table, parse_fasta, save_figure


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
    return out[valid].drop_duplicates(["canonical_UniProtAC", "corrected_position"]).copy()


def first_second_neighbor_distances(positions: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    positions = np.sort(np.asarray(positions, dtype=int))
    n = len(positions)
    d1 = np.full(n, np.inf, dtype=float)
    d2 = np.full(n, np.inf, dtype=float)
    if n <= 1:
        return d1, d2
    for idx, pos in enumerate(positions):
        distances = np.abs(positions - pos)
        distances = np.sort(distances[distances > 0])
        if len(distances) >= 1:
            d1[idx] = distances[0]
        if len(distances) >= 2:
            d2[idx] = distances[1]
    return d1, d2


def feature_table(positions_by_acc: dict[str, np.ndarray], label: str) -> pd.DataFrame:
    rows = []
    for accession, positions in positions_by_acc.items():
        d1, d2 = first_second_neighbor_distances(positions)
        for pos, first, second in zip(np.sort(positions), d1, d2):
            rows.append(
                {
                    "canonical_UniProtAC": accession,
                    "position": int(pos),
                    "nearest_methyl_arg_distance": first,
                    "second_nearest_methyl_arg_distance": second,
                    "set": label,
                }
            )
    return pd.DataFrame(rows)


def exact_count_matrix(features: pd.DataFrame, max_distance: int) -> np.ndarray:
    matrix = np.zeros((max_distance, max_distance), dtype=float)
    subset = features[
        features["nearest_methyl_arg_distance"].between(1, max_distance, inclusive="both")
        & features["second_nearest_methyl_arg_distance"].between(1, max_distance, inclusive="both")
    ].copy()
    if subset.empty:
        return matrix
    x = subset["nearest_methyl_arg_distance"].astype(int).to_numpy()
    y = subset["second_nearest_methyl_arg_distance"].astype(int).to_numpy()
    for first, second in zip(x, y):
        matrix[second - 1, first - 1] += 1
    return matrix


def null_count_matrix(
    arg_positions_by_acc: dict[str, np.ndarray],
    methyl_counts_by_acc: dict[str, int],
    max_distance: int,
    permutations: int,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    accessions = [
        acc for acc, count in methyl_counts_by_acc.items()
        if count > 0 and len(arg_positions_by_acc.get(acc, [])) >= count
    ]
    expected = np.zeros((max_distance, max_distance), dtype=float)
    for _ in range(permutations):
        sampled = {
            acc: np.sort(rng.choice(arg_positions_by_acc[acc], size=methyl_counts_by_acc[acc], replace=False))
            for acc in accessions
        }
        expected += exact_count_matrix(feature_table(sampled, "null"), max_distance)
    return expected / permutations


def save_plot(fig, prefix: Path) -> None:
    save_figure(fig, prefix)
    fig.savefig(prefix.with_suffix(".svg"), bbox_inches="tight", facecolor="white")


def smooth_with_nan_weights(matrix: np.ndarray, sigma: float) -> np.ndarray:
    values = np.where(np.isfinite(matrix), matrix, 0.0)
    weights = np.isfinite(matrix).astype(float)
    smooth_values = gaussian_filter(values, sigma=sigma, mode="constant")
    smooth_weights = gaussian_filter(weights, sigma=sigma, mode="constant")
    out = np.full_like(matrix, np.nan, dtype=float)
    np.divide(smooth_values, smooth_weights, out=out, where=smooth_weights > 0)
    return out


def matrix_to_table(observed: np.ndarray, expected_per_permutation: np.ndarray, permutations: int, observed_total: int) -> pd.DataFrame:
    expected_total = expected_per_permutation * permutations
    rows = []
    max_distance = observed.shape[0]
    for second_idx in range(max_distance):
        for first_idx in range(max_distance):
            observed_count = float(observed[second_idx, first_idx])
            null_count = float(expected_total[second_idx, first_idx])
            table = [
                [int(round(observed_count)), int(round(observed_total - observed_count))],
                [int(round(null_count)), int(round(observed_total * permutations - null_count))],
            ]
            odds_ratio, p_value = fisher_exact(table, alternative="two-sided")
            if odds_ratio == 0:
                log2_or = -np.inf
            elif np.isinf(odds_ratio):
                log2_or = np.inf
            else:
                log2_or = np.log2(odds_ratio)
            rows.append(
                {
                    "nearest_distance_aa": first_idx + 1,
                    "second_nearest_distance_aa": second_idx + 1,
                    "observed_count": observed_count,
                    "observed_total_sites": observed_total,
                    "null_count_across_permutations": null_count,
                    "null_total_sites_across_permutations": observed_total * permutations,
                    "expected_count_per_permutation": expected_per_permutation[second_idx, first_idx],
                    "odds_ratio": odds_ratio,
                    "log2_odds_ratio": log2_or,
                    "fisher_p_value": p_value,
                }
            )
    out = add_q_values(pd.DataFrame(rows), p_col="fisher_p_value", out_col="q_value")
    out["log2_odds_ratio"] = out["log2_odds_ratio"].replace([np.inf, -np.inf], np.nan)
    return out


def plot_density(
    observed: np.ndarray,
    stats_table: pd.DataFrame,
    features: pd.DataFrame,
    outdir: Path,
    max_distance: int,
    display_max_distance: int,
    smooth_sigma: float,
) -> None:
    display_n = min(max_distance, display_max_distance)
    observed_display = observed[:display_n, :display_n]
    triangle_mask = np.fromfunction(lambda row, col: row <= col, (display_n, display_n), dtype=int)
    # Plot as x=second-nearest, y=nearest. This rotates the original matrix and exposes
    # the valid nearest <= second-nearest triangle rather than implying a full square.
    observed_rotated = observed_display.T
    smoothed = gaussian_filter(observed_rotated, sigma=smooth_sigma, mode="constant")
    smoothed = np.where(triangle_mask, smoothed, np.nan)
    log2_or = (
        stats_table.pivot(index="second_nearest_distance_aa", columns="nearest_distance_aa", values="log2_odds_ratio")
        .reindex(index=range(1, max_distance + 1), columns=range(1, max_distance + 1))
        .iloc[:display_n, :display_n]
        .to_numpy(dtype=float)
    )
    log2_or = np.where(np.isfinite(log2_or), np.clip(log2_or, -4.0, 4.0), np.nan)
    smoothed_log2_or = smooth_with_nan_weights(log2_or.T, sigma=smooth_sigma)
    smoothed_log2_or = np.where(triangle_mask, smoothed_log2_or, np.nan)
    q_values = (
        stats_table.pivot(index="second_nearest_distance_aa", columns="nearest_distance_aa", values="q_value")
        .reindex(index=range(1, max_distance + 1), columns=range(1, max_distance + 1))
        .iloc[:display_n, :display_n]
        .to_numpy(dtype=float)
    )

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.0, 5.7), gridspec_kw={"width_ratios": [1.08, 1.0], "wspace": 0.28})
    positive_density = smoothed[smoothed > 0]
    vmax_density = max(1.0, float(np.nanpercentile(positive_density, 99.3)) if len(positive_density) else 1.0)
    density_cmap = plt.get_cmap("magma").copy()
    density_cmap.set_bad(color="white")
    im = ax.imshow(
        smoothed,
        origin="lower",
        aspect="equal",
        cmap=density_cmap,
        vmin=0,
        vmax=vmax_density,
        extent=[0.5, display_n + 0.5, 0.5, display_n + 0.5],
        interpolation="bilinear",
    )
    ax.set_title("Observed short-range 2D density")
    ax.set_xlabel("Second-nearest methyl-Arg distance (aa)")
    ax.set_ylabel("Nearest methyl-Arg distance (aa)")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label("Smoothed site density")

    finite_or = smoothed_log2_or[np.isfinite(smoothed_log2_or)]
    vmax = max(0.8, float(np.nanpercentile(np.abs(finite_or), 98)) if len(finite_or) else 0.8)
    or_cmap = plt.get_cmap("RdGy_r").copy()
    or_cmap.set_bad(color="white")
    im2 = ax2.imshow(
        smoothed_log2_or,
        origin="lower",
        aspect="equal",
        cmap=or_cmap,
        vmin=-vmax,
        vmax=vmax,
        extent=[0.5, display_n + 0.5, 0.5, display_n + 0.5],
        interpolation="bilinear",
    )
    ax2.set_title("Smoothed Fisher odds ratio")
    ax2.set_xlabel("Second-nearest methyl-Arg distance (aa)")
    ax2.set_ylabel("")
    cbar2 = fig.colorbar(im2, ax=ax2, fraction=0.046, pad=0.02)
    cbar2.set_label("log2 odds ratio")

    for axis in (ax, ax2):
        axis.set_xlim(0.5, display_n + 0.5)
        axis.set_ylim(0.5, display_n + 0.5)
        ticks = [tick for tick in [1, 5, 10, 20, 30] if tick <= display_n]
        axis.set_xticks(ticks)
        axis.set_yticks(ticks)
        axis.plot([0.5, display_n + 0.5], [0.5, display_n + 0.5], color="#5F6670", linewidth=0.8, alpha=0.35)
        axis.fill_between([0.5, display_n + 0.5], [0.5, display_n + 0.5], [display_n + 0.5, display_n + 0.5], color="white", zorder=3)
        for spine in axis.spines.values():
            spine.set_visible(False)

    significant = np.where((q_values <= 0.05) & (observed_display >= 3))
    if len(significant[0]):
        ax2.scatter(significant[0] + 1, significant[1] + 1, s=5, color="#111111", alpha=0.45, linewidths=0, zorder=4)

    in_window = features[
        features["nearest_methyl_arg_distance"].between(1, display_n, inclusive="both")
        & features["second_nearest_methyl_arg_distance"].between(1, display_n, inclusive="both")
    ]
    within10 = features[
        features["nearest_methyl_arg_distance"].between(1, 10, inclusive="both")
        & features["second_nearest_methyl_arg_distance"].between(1, 10, inclusive="both")
    ]
    fig.suptitle("Amino-acid-resolution methylarginine nearest-neighbor density", y=0.98, fontsize=13)
    fig.text(
        0.99,
        0.01,
        f"Smoothed from exact 1 aa bins, shown to {display_n} aa. White triangle is impossible because nearest distance cannot exceed second-nearest distance; black points mark q<=0.05 and exact-bin count >=3.",
        ha="right",
        va="bottom",
        fontsize=7.8,
        color=BREWER_COLORS["dark_gray"],
    )
    fig.subplots_adjust(left=0.08, right=0.96, top=0.86, bottom=0.13)
    save_plot(fig, outdir / "arg_methyl_knearest_aa_resolution_density")


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot amino-acid-resolution 2D density of first and second nearest methylarginine distances.")
    parser.add_argument("--sites", default="PTM_results/domain_context/sites_with_domain_context.tsv")
    parser.add_argument("--fasta", default="PTM_data/context/uniprot_human_reviewed_canonical.fasta")
    parser.add_argument("--outdir", default="PTM_results/methyl_arg_clustering")
    parser.add_argument("--max-distance", type=int, default=100)
    parser.add_argument("--display-max-distance", type=int, default=30)
    parser.add_argument("--smooth-sigma", type=float, default=1.2)
    parser.add_argument("--permutations", type=int, default=250)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    sites = pd.read_csv(args.sites, sep="\t", low_memory=False)
    sequences = parse_fasta(args.fasta)
    sites = validate_sites(sites, sequences)
    observed_by_acc = {
        accession: np.sort(group["corrected_position"].to_numpy(dtype=int))
        for accession, group in sites.groupby("canonical_UniProtAC", sort=False)
    }
    arg_positions_by_acc = {acc: residue_positions(sequences[acc], "R") for acc in observed_by_acc if acc in sequences}
    methyl_counts_by_acc = {acc: len(pos) for acc, pos in observed_by_acc.items()}
    features = feature_table(observed_by_acc, "observed")
    observed = exact_count_matrix(features, args.max_distance)
    expected = null_count_matrix(arg_positions_by_acc, methyl_counts_by_acc, args.max_distance, args.permutations, args.seed)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    features.to_csv(outdir / "arg_methyl_knearest_aa_resolution_site_features.tsv", sep="\t", index=False)
    stats_table = matrix_to_table(observed, expected, args.permutations, observed_total=len(features))
    stats_table.to_csv(outdir / "arg_methyl_knearest_aa_resolution_density.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "sites_analyzed": int(len(features)),
                "proteins_analyzed": int(len(observed_by_acc)),
                "max_distance": int(args.max_distance),
                "permutations": int(args.permutations),
                "sites_with_two_neighbors_within_max_distance": int(observed.sum()),
            }
        ]
    ).to_csv(outdir / "arg_methyl_knearest_aa_resolution_density_qc.tsv", sep="\t", index=False)
    plot_density(observed, stats_table, features, outdir, args.max_distance, args.display_max_distance, args.smooth_sigma)
    print(outdir / "arg_methyl_knearest_aa_resolution_density.svg")
    print(outdir / "arg_methyl_knearest_aa_resolution_density.tsv")


if __name__ == "__main__":
    main()

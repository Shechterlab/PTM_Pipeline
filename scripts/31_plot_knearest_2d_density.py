from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import BREWER_COLORS, canonical_site_table, parse_fasta, save_figure, style_axis


BIN_EDGES = np.array([0, 1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40, 60, 80, 100, np.inf], dtype=float)
BIN_LABELS = ["1", "2", "3", "4", "5", "6", "7-8", "9-10", "11-15", "16-20", "21-30", "31-40", "41-60", "61-80", "81-100", ">100/none"]


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
        candidates = []
        for neighbor_idx in (idx - 2, idx - 1, idx + 1, idx + 2):
            if 0 <= neighbor_idx < n:
                candidates.append(abs(int(positions[neighbor_idx]) - int(pos)))
        candidates = sorted(candidates)
        if candidates:
            d1[idx] = candidates[0]
        if len(candidates) >= 2:
            d2[idx] = candidates[1]
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


def hist2d_density(features: pd.DataFrame) -> np.ndarray:
    x = features["nearest_methyl_arg_distance"].replace(np.inf, 101).to_numpy(dtype=float)
    y = features["second_nearest_methyl_arg_distance"].replace(np.inf, 101).to_numpy(dtype=float)
    counts, _, _ = np.histogram2d(x, y, bins=[BIN_EDGES, BIN_EDGES])
    total = counts.sum()
    return counts / total if total else counts


def null_density(
    arg_positions_by_acc: dict[str, np.ndarray],
    methyl_counts_by_acc: dict[str, int],
    permutations: int,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    accs = [acc for acc, count in methyl_counts_by_acc.items() if count > 0 and len(arg_positions_by_acc.get(acc, [])) >= count]
    acc_densities = []
    for _ in range(permutations):
        sampled = {
            acc: np.sort(rng.choice(arg_positions_by_acc[acc], size=methyl_counts_by_acc[acc], replace=False))
            for acc in accs
        }
        acc_densities.append(hist2d_density(feature_table(sampled, "null")))
    return np.mean(np.stack(acc_densities), axis=0)


def save_plot(fig, prefix: Path) -> None:
    save_figure(fig, prefix)
    fig.savefig(prefix.with_suffix(".svg"), bbox_inches="tight", facecolor="white")


def matrix_to_table(observed: np.ndarray, expected: np.ndarray) -> pd.DataFrame:
    rows = []
    ratio = np.divide(observed, expected, out=np.full_like(observed, np.nan), where=expected > 0)
    log2_ratio = np.full_like(observed, np.nan, dtype=float)
    valid = ratio > 0
    log2_ratio[valid] = np.log2(ratio[valid])
    for i, x_label in enumerate(BIN_LABELS):
        for j, y_label in enumerate(BIN_LABELS):
            rows.append(
                {
                    "nearest_distance_bin": x_label,
                    "second_nearest_distance_bin": y_label,
                    "observed_probability": observed[i, j],
                    "arginine_null_probability": expected[i, j],
                    "log2_observed_over_null": log2_ratio[i, j],
                }
            )
    return pd.DataFrame(rows)


def plot_density(observed: np.ndarray, expected: np.ndarray, outdir: Path) -> None:
    table = matrix_to_table(observed, expected)
    log2_matrix = table.pivot(index="second_nearest_distance_bin", columns="nearest_distance_bin", values="log2_observed_over_null").reindex(index=BIN_LABELS[::-1], columns=BIN_LABELS)
    obs_matrix = table.pivot(index="second_nearest_distance_bin", columns="nearest_distance_bin", values="observed_probability").reindex(index=BIN_LABELS[::-1], columns=BIN_LABELS)

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.6, 5.6), gridspec_kw={"width_ratios": [1.0, 1.0], "wspace": 0.28})

    obs_pct = 100.0 * obs_matrix.to_numpy(dtype=float)
    im = ax.imshow(obs_pct, aspect="auto", cmap="magma", vmin=0)
    ax.set_xticks(range(len(BIN_LABELS)), BIN_LABELS, rotation=40, ha="right")
    ax.set_yticks(range(len(BIN_LABELS)), BIN_LABELS[::-1])
    ax.set_xlabel("Nearest methyl-Arg distance (aa)")
    ax.set_ylabel("Second-nearest methyl-Arg distance (aa)")
    ax.set_title("Observed 2D probability density")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label("% of methyl-sites")

    values = log2_matrix.to_numpy(dtype=float)
    vmax = max(0.8, float(np.nanmax(np.abs(values))) if np.isfinite(values).any() else 0.8)
    im2 = ax2.imshow(values, aspect="auto", cmap="RdGy_r", vmin=-vmax, vmax=vmax)
    ax2.set_xticks(range(len(BIN_LABELS)), BIN_LABELS, rotation=40, ha="right")
    ax2.set_yticks(range(len(BIN_LABELS)), BIN_LABELS[::-1])
    ax2.set_xlabel("Nearest methyl-Arg distance (aa)")
    ax2.set_ylabel("")
    ax2.set_title("Enrichment over same-protein arginine null")
    cbar2 = fig.colorbar(im2, ax=ax2, fraction=0.046, pad=0.02)
    cbar2.set_label("log2(observed / null)")

    for axis in (ax, ax2):
        for spine in axis.spines.values():
            spine.set_visible(False)
        axis.tick_params(length=0)

    fig.suptitle("Methylarginine nearest-neighbor density landscape", y=0.98, fontsize=13)
    fig.text(
        0.99,
        0.01,
        "Each methyl-site contributes one point: distance to nearest and second-nearest methyl-Arg in the same protein.",
        ha="right",
        va="bottom",
        fontsize=7.8,
        color=BREWER_COLORS["dark_gray"],
    )
    fig.subplots_adjust(left=0.08, right=0.97, top=0.86, bottom=0.22)
    save_plot(fig, outdir / "arg_methyl_knearest_2d_probability_density")
    table.to_csv(outdir / "arg_methyl_knearest_2d_probability_density.tsv", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot 2D methylarginine probability density from first and second nearest-neighbor distances.")
    parser.add_argument("--sites", default="PTM_results/domain_context/sites_with_domain_context.tsv")
    parser.add_argument("--fasta", default="PTM_data/context/uniprot_human_reviewed_canonical.fasta")
    parser.add_argument("--outdir", default="PTM_results/methyl_arg_clustering")
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

    observed_features = feature_table(observed_by_acc, "observed")
    observed_density = hist2d_density(observed_features)
    expected_density = null_density(arg_positions_by_acc, methyl_counts_by_acc, args.permutations, args.seed)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    observed_features.to_csv(outdir / "arg_methyl_knearest_site_features.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "sites_analyzed": int(len(observed_features)),
                "proteins_analyzed": int(len(observed_by_acc)),
                "permutations": int(args.permutations),
                "x_axis": "nearest methylarginine distance",
                "y_axis": "second-nearest methylarginine distance",
            }
        ]
    ).to_csv(outdir / "arg_methyl_knearest_2d_probability_density_qc.tsv", sep="\t", index=False)

    plot_density(observed_density, expected_density, outdir)
    print(outdir / "arg_methyl_knearest_2d_probability_density.svg")
    print(outdir / "arg_methyl_knearest_2d_probability_density.tsv")


if __name__ == "__main__":
    main()

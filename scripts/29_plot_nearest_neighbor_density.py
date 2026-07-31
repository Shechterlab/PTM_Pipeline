from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import BREWER_COLORS, apply_paper_style, save_figure, style_axis


BIN_EDGES = [0, 1, 2, 3, 5, 10, 20, 40, 100, np.inf]
BIN_LABELS = ["1", "2", "3", "4-5", "6-10", "11-20", "21-40", "41-100", ">100/no neighbor"]


def nearest_histogram(per_site: pd.DataFrame) -> pd.DataFrame:
    distances = per_site["nearest_methyl_arg_distance"].copy()
    distances = distances.fillna(np.inf)
    bins = pd.cut(distances, bins=BIN_EDGES, labels=BIN_LABELS, include_lowest=True, right=True)
    counts = bins.value_counts().reindex(BIN_LABELS, fill_value=0)
    out = counts.rename_axis("nearest_neighbor_distance_bin").reset_index(name="site_count")
    out["site_fraction"] = out["site_count"] / out["site_count"].sum()
    return out


def null_histogram_from_curve(curve: pd.DataFrame) -> pd.DataFrame:
    lookup = curve.set_index("distance_aa")["null_mean_prob_neighbor_within_x"].to_dict()
    rows = []
    previous = 0.0
    for label, low, high in zip(BIN_LABELS, BIN_EDGES[:-1], BIN_EDGES[1:]):
        if np.isinf(high):
            fraction = 1.0 - float(lookup.get(int(low), previous))
        else:
            hi = float(lookup.get(int(high), previous))
            lo = float(lookup.get(int(low), 0.0)) if low > 0 else 0.0
            fraction = hi - lo
            previous = hi
        rows.append({"nearest_neighbor_distance_bin": label, "null_mean_site_fraction": fraction})
    return pd.DataFrame(rows)


def save_plot(fig, prefix: Path) -> None:
    save_figure(fig, prefix)
    fig.savefig(prefix.with_suffix(".svg"), bbox_inches="tight", facecolor="white")


def plot_density(hist: pd.DataFrame, curve: pd.DataFrame, outdir: Path) -> None:
    plot_df = hist.copy()
    null_hist = null_histogram_from_curve(curve)
    plot_df = plot_df.merge(null_hist, on="nearest_neighbor_distance_bin", how="left")
    x = np.arange(len(plot_df))

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(7.6, 6.3), gridspec_kw={"height_ratios": [1.35, 1.65], "hspace": 0.38})
    width = 0.38
    ax.bar(
        x - width / 2,
        100 * plot_df["site_fraction"],
        width=width,
        color="#CBB2B9",
        edgecolor="white",
        linewidth=0.8,
        label="Observed",
    )
    ax.bar(
        x + width / 2,
        100 * plot_df["null_mean_site_fraction"],
        width=width,
        color=BREWER_COLORS["light_gray"],
        edgecolor="white",
        linewidth=0.8,
        label="Arginine null",
    )
    style_axis(ax, grid_axis="y")
    ax.set_xticks(x, plot_df["nearest_neighbor_distance_bin"])
    ax.set_ylabel("% of methyl-sites")
    ax.set_title("Nearest methylarginine spacing across the methylome")
    ax.legend(frameon=False, loc="upper right")

    curve40 = curve[curve["distance_aa"] <= 40].copy()
    ax2.fill_between(
        curve40["distance_aa"],
        100 * curve40["null_q025_prob_neighbor_within_x"],
        100 * curve40["null_q975_prob_neighbor_within_x"],
        color=BREWER_COLORS["light_gray"],
        alpha=0.28,
        linewidth=0,
        label="Null 95% interval",
    )
    ax2.plot(
        curve40["distance_aa"],
        100 * curve40["null_mean_prob_neighbor_within_x"],
        color=BREWER_COLORS["mid_gray"],
        linewidth=1.8,
        linestyle="--",
        label="Arginine null mean",
    )
    ax2.plot(
        curve40["distance_aa"],
        100 * curve40["prob_neighbor_within_x"],
        color="#8E747C",
        linewidth=2.4,
        label="Observed",
    )
    style_axis(ax2, grid_axis="y")
    ax2.set_xlim(1, 40)
    ax2.set_ylim(0, max(62, 100 * curve40["prob_neighbor_within_x"].max() + 4))
    ax2.set_xlabel("Distance to nearest methylarginine in same protein (aa)")
    ax2.set_ylabel("% with another methyl-Arg within X aa")
    ax2.legend(frameon=False, loc="lower right", fontsize=8.4)

    fig.text(
        0.99,
        0.01,
        "Null samples the same number of methyl-sites from arginines in the same proteins, matched by IDR/non-IDR context where available.",
        ha="right",
        va="bottom",
        fontsize=7.8,
        color=BREWER_COLORS["dark_gray"],
    )
    fig.subplots_adjust(left=0.1, right=0.98, top=0.92, bottom=0.12)
    save_plot(fig, outdir / "arg_methyl_nearest_neighbor_density")

    plot_df.to_csv(outdir / "arg_methyl_nearest_neighbor_distance_histogram.tsv", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot nearest-neighbor methylarginine distance density.")
    parser.add_argument("--per-site", default="PTM_results/methyl_arg_clustering/per_site_nearest_methyl_arg_distance.tsv")
    parser.add_argument("--curve", default="PTM_results/methyl_arg_clustering/nearest_neighbor_curve_vs_null.tsv")
    parser.add_argument("--outdir", default="PTM_results/methyl_arg_clustering")
    args = parser.parse_args()

    per_site = pd.read_csv(args.per_site, sep="\t")
    curve = pd.read_csv(args.curve, sep="\t")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    hist = nearest_histogram(per_site)
    apply_paper_style()
    plot_density(hist, curve, outdir)
    print(outdir / "arg_methyl_nearest_neighbor_density.svg")
    print(outdir / "arg_methyl_nearest_neighbor_distance_histogram.tsv")


if __name__ == "__main__":
    main()

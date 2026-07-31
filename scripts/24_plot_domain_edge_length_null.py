from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import BREWER_COLORS, apply_paper_style, save_figure, style_axis


BIN_ORDER = ["Domain edge <=20 aa", "Domain edge 21-40 aa", "Domain-distal >40 aa"]
DISPLAY_LABELS = {
    "Domain edge <=20 aa": "<=20 aa",
    "Domain edge 21-40 aa": "21-40 aa",
    "Domain-distal >40 aa": ">40 aa",
}


def save_plot(fig, path_prefix: Path) -> None:
    save_figure(fig, path_prefix)
    fig.savefig(path_prefix.with_suffix(".svg"), bbox_inches="tight", facecolor="white")


def plot_overall(length_null: pd.DataFrame, qc: pd.DataFrame, outdir: Path) -> None:
    plot_df = length_null.set_index("category").loc[BIN_ORDER].reset_index()
    labels = [DISPLAY_LABELS[item] for item in plot_df["category"]]
    x = np.arange(len(plot_df))

    observed_pct = 100.0 * plot_df["observed_fraction"].to_numpy(dtype=float)
    expected_pct = 100.0 * plot_df["expected_fraction_length_conditioned"].to_numpy(dtype=float)
    log2_values = plot_df["log2_enrichment_vs_length_null"].to_numpy(dtype=float)

    fig, (ax, ax2) = plt.subplots(
        2,
        1,
        figsize=(6.9, 5.9),
        gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.34},
    )

    bar_color = "#CBB2B9"
    marker_color = BREWER_COLORS["dark_gray"]
    ax.bar(x, observed_pct, width=0.58, color=bar_color, edgecolor="white", linewidth=0.8, label="Observed methyl-sites")
    ax.scatter(x, expected_pct, marker="D", s=42, color=marker_color, zorder=3, label="Expected from IDR length")
    for xpos, obs, exp, count in zip(x, observed_pct, expected_pct, plot_df["observed_count"]):
        high = max(obs, exp)
        ax.text(xpos, high + 2.0, f"{int(count)}", ha="center", va="bottom", fontsize=8.6, color=BREWER_COLORS["dark_gray"])

    style_axis(ax)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("% of methyl-sites")
    ax.set_title("Domain-edge proximity after conditioning on IDR segment length")
    ax.legend(frameon=False, loc="upper right", fontsize=8.6)
    ax.set_ylim(0, max(observed_pct.max(), expected_pct.max()) * 1.22)

    signed_colors = [BREWER_COLORS["mid_gray"] if value >= 0 else BREWER_COLORS["light_gray"] for value in log2_values]
    ax2.bar(x, log2_values, width=0.58, color=signed_colors, edgecolor="white", linewidth=0.8)
    style_axis(ax2, zero="y")
    ax2.axhline(0, color=BREWER_COLORS["dark_gray"], linewidth=0.8)
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels)
    ax2.set_ylabel("log2(obs/exp)")
    ax2.set_ylim(-0.18, 0.18)
    for xpos, value, q_value in zip(x, log2_values, plot_df["q_value"]):
        va = "bottom" if value >= 0 else "top"
        offset = 0.018 if value >= 0 else -0.018
        ax2.text(xpos, value + offset, f"q={q_value:.2g}", ha="center", va=va, fontsize=8.2, color=BREWER_COLORS["dark_gray"])

    row = qc.iloc[0]
    fig.text(
        0.99,
        0.01,
        (
            f"Disordered outside-domain methyl-sites only (n={int(row['methylated_disordered_outside_domain_sites'])}); "
            f"expected fractions from arginines in the same outside-domain IDR segments. "
            f"Median site-weighted IDR segment length={row['median_segment_length_methyl_sites_site_weighted']:.0f} aa."
        ),
        ha="right",
        va="bottom",
        fontsize=7.8,
        color=BREWER_COLORS["dark_gray"],
    )
    fig.subplots_adjust(left=0.11, right=0.98, top=0.9, bottom=0.14, hspace=0.42)
    save_plot(fig, outdir / "arg_methyl_domain_edge_idr_length_null")


def plot_by_segment_length(length_strata: pd.DataFrame, outdir: Path) -> None:
    length_order = ["<=40 aa", "41-80 aa", "81-160 aa", ">160 aa"]
    plot_df = length_strata[length_strata["category"].isin(BIN_ORDER)].copy()
    plot_df["outside_domain_idr_length_bin"] = pd.Categorical(
        plot_df["outside_domain_idr_length_bin"],
        categories=length_order,
        ordered=True,
    )
    plot_df["category"] = pd.Categorical(plot_df["category"], categories=BIN_ORDER, ordered=True)
    plot_df = plot_df.sort_values(["outside_domain_idr_length_bin", "category"])

    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    offsets = np.array([-0.23, 0.0, 0.23])
    palette = {
        "Domain edge <=20 aa": "#CBB2B9",
        "Domain edge 21-40 aa": BREWER_COLORS["mid_gray"],
        "Domain-distal >40 aa": BREWER_COLORS["light_gray"],
    }
    base = np.arange(len(length_order))
    for idx, category in enumerate(BIN_ORDER):
        subset = plot_df[plot_df["category"] == category].set_index("outside_domain_idr_length_bin").loc[length_order]
        xpos = base + offsets[idx]
        values = subset["log2_enrichment_vs_length_null"].to_numpy(dtype=float)
        ax.bar(
            xpos,
            values,
            width=0.21,
            color=palette[category],
            edgecolor="white",
            linewidth=0.8,
            label=DISPLAY_LABELS[category],
        )

    style_axis(ax, zero="y")
    ax.axhline(0, color=BREWER_COLORS["dark_gray"], linewidth=0.8)
    ax.set_xticks(base)
    ax.set_xticklabels(length_order)
    ax.set_ylabel("log2(obs/exp)")
    ax.set_xlabel("Outside-domain IDR segment length")
    ax.set_title("Length-null residuals by IDR segment length")
    ax.legend(frameon=False, loc="upper right", fontsize=8.4, ncol=3)
    ax.set_ylim(-0.62, 0.18)
    fig.text(
        0.99,
        0.01,
        "Expected values are calculated within the same outside-domain IDR segments as the observed methyl-sites.",
        ha="right",
        va="bottom",
        fontsize=7.8,
        color=BREWER_COLORS["dark_gray"],
    )
    fig.subplots_adjust(left=0.1, right=0.98, top=0.88, bottom=0.17)
    save_plot(fig, outdir / "arg_methyl_domain_edge_idr_length_null_by_segment_length")


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot domain-edge proximity after conditioning on IDR segment length.")
    parser.add_argument("--length-null", default="PTM_results/domain_context/domain_edge_bin_enrichment_disordered_only_compact_idr_length_null.tsv")
    parser.add_argument(
        "--length-strata",
        default="PTM_results/domain_context/domain_edge_bin_enrichment_disordered_only_compact_idr_length_null_by_segment_length.tsv",
    )
    parser.add_argument("--qc", default="PTM_results/domain_context/domain_edge_idr_length_null_qc.tsv")
    parser.add_argument("--outdir", default="PTM_results/domain_context")
    args = parser.parse_args()

    apply_paper_style()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    length_null = pd.read_csv(args.length_null, sep="\t")
    length_strata = pd.read_csv(args.length_strata, sep="\t")
    qc = pd.read_csv(args.qc, sep="\t")
    plot_overall(length_null, qc, outdir)
    plot_by_segment_length(length_strata, outdir)


if __name__ == "__main__":
    main()

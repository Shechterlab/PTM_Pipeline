from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from common import BREWER_COLORS, apply_paper_style, save_figure, style_axis


PROXIMAL_BIN = "domain-proximal <=20 aa"
DISTAL_BIN = "domain-distal >20 aa"
ARCHITECTURE_ORDER = [
    "broad across IDR",
    "domain-proximal clustered",
    "domain-distal clustered",
]


def collapse_bin(label: str) -> str:
    if str(label) == "Domain edge <=20 aa":
        return PROXIMAL_BIN
    return DISTAL_BIN


def classify_architecture(classes: set[str]) -> str:
    if PROXIMAL_BIN in classes and DISTAL_BIN in classes:
        return "broad across IDR"
    if PROXIMAL_BIN in classes:
        return "domain-proximal clustered"
    return "domain-distal clustered"


def best_window(positions: list[int], window: int = 20) -> tuple[int, int, int, str]:
    positions = sorted(set(int(pos) for pos in positions))
    if not positions:
        return 0, 0, 0, ""
    best_positions = [positions[0]]
    for idx, start in enumerate(positions):
        end = start + window - 1
        current = [pos for pos in positions[idx:] if pos <= end]
        if len(current) > len(best_positions):
            best_positions = current
    return (
        len(best_positions),
        int(best_positions[0]),
        int(best_positions[-1]),
        ",".join(str(pos) for pos in best_positions),
    )


def build_cluster_examples(site_context_path: Path, outdir: Path) -> Path:
    sites = pd.read_csv(site_context_path, sep="\t", low_memory=False)
    required = {
        "canonical_UniProtAC",
        "corrected_position",
        "domain_edge_bin_compact",
        "nearest_domain_edge_distance",
        "nearest_domain_class",
    }
    missing = required - set(sites.columns)
    if missing:
        raise ValueError(f"Cannot build cluster examples; missing columns: {', '.join(sorted(missing))}")

    sites = sites[sites["domain_edge_bin_compact"].notna()].copy()
    if "disorder_context_class" in sites.columns:
        sites = sites[sites["disorder_context_class"] == "disordered"].copy()
    sites["position"] = pd.to_numeric(sites["corrected_position"], errors="coerce")
    sites = sites[sites["position"].notna()].copy()
    sites["position"] = sites["position"].astype(int)

    rows = []
    for (accession, edge_bin), group in sites.groupby(["canonical_UniProtAC", "domain_edge_bin_compact"], dropna=False):
        positions = sorted(set(group["position"].astype(int)))
        if len(positions) < 2:
            continue
        max_count, start, end, best_positions = best_window(positions, window=20)
        rows.append(
            {
                "canonical_UniProtAC": accession,
                "gene": str(group.get("substrate_genename", pd.Series([""])).dropna().astype(str).iloc[0]) if "substrate_genename" in group else "",
                "domain_edge_bin_compact": edge_bin,
                "total_sites_in_bin": len(positions),
                "max_sites_in_20aa_span_within_bin": max_count,
                "best_window_start": start,
                "best_window_end": end,
                "best_window_positions": best_positions,
                "median_domain_edge_distance_in_bin": pd.to_numeric(group["nearest_domain_edge_distance"], errors="coerce").median(),
                "median_idr_segment_length": pd.NA,
                "nearest_domain_classes": " ; ".join(sorted(set(group["nearest_domain_class"].dropna().astype(str)))),
            }
        )
    out = pd.DataFrame(rows)
    out_path = outdir / "domain_edge_class_20aa_cluster_examples.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    return out_path


def save_plot(fig, prefix: Path) -> None:
    save_figure(fig, prefix)
    fig.savefig(prefix.with_suffix(".svg"), bbox_inches="tight", facecolor="white")


def summarize(input_path: Path, outdir: Path, cluster_threshold: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    clusters = pd.read_csv(input_path, sep="\t")
    clusters = clusters[clusters["max_sites_in_20aa_span_within_bin"] >= cluster_threshold].copy()
    clusters["collapsed_domain_adjacency_class"] = clusters["domain_edge_bin_compact"].map(collapse_bin)

    protein_profiles = (
        clusters.groupby("canonical_UniProtAC")
        .agg(
            gene=("gene", "first"),
            architecture_class=("collapsed_domain_adjacency_class", lambda x: classify_architecture(set(x))),
            clustered_adjacency_classes=("collapsed_domain_adjacency_class", lambda x: " | ".join(sorted(set(x)))),
            original_edge_bins=("domain_edge_bin_compact", lambda x: " | ".join(sorted(set(x)))),
            max_sites_in_20aa_span=("max_sites_in_20aa_span_within_bin", "max"),
            total_sites_in_clustered_bins=("total_sites_in_bin", "sum"),
            best_cluster_positions=("best_window_positions", lambda x: " ; ".join(map(str, x))),
            median_domain_edge_distance=("median_domain_edge_distance_in_bin", "median"),
            median_idr_segment_length=("median_idr_segment_length", "median"),
            nearest_domain_classes=("nearest_domain_classes", lambda x: " ; ".join(sorted(set(map(str, x))))),
        )
        .reset_index()
    )
    protein_profiles["architecture_class"] = pd.Categorical(
        protein_profiles["architecture_class"], categories=ARCHITECTURE_ORDER, ordered=True
    )
    protein_profiles = protein_profiles.sort_values(
        ["architecture_class", "max_sites_in_20aa_span", "total_sites_in_clustered_bins"],
        ascending=[True, False, False],
    )

    architecture_summary = (
        protein_profiles.groupby("architecture_class", observed=True)
        .agg(
            n_proteins=("canonical_UniProtAC", "nunique"),
            median_max_sites_in_20aa_span=("max_sites_in_20aa_span", "median"),
            median_total_sites_in_clustered_bins=("total_sites_in_clustered_bins", "median"),
            median_domain_edge_distance=("median_domain_edge_distance", "median"),
            median_idr_segment_length=("median_idr_segment_length", "median"),
        )
        .reindex(ARCHITECTURE_ORDER)
        .reset_index()
    )

    expanded_rows = []
    for row in protein_profiles.itertuples(index=False):
        classes = str(row.nearest_domain_classes).replace(" ; ", ";").split(";")
        for domain_class in sorted({item.strip() for item in classes if item.strip() and item.strip() != "nan"}):
            expanded_rows.append(
                {
                    "canonical_UniProtAC": row.canonical_UniProtAC,
                    "gene": row.gene,
                    "architecture_class": row.architecture_class,
                    "nearest_domain_class": domain_class,
                }
            )
    expanded = pd.DataFrame(expanded_rows)
    domain_class_summary = (
        expanded.groupby(["architecture_class", "nearest_domain_class"], observed=True)
        .agg(n_proteins=("canonical_UniProtAC", "nunique"))
        .reset_index()
        .sort_values(["architecture_class", "n_proteins"], ascending=[True, False])
    )

    suffix = f"ge{cluster_threshold}"
    protein_profiles.to_csv(outdir / f"domain_cluster_architecture_classes_{suffix}.tsv", sep="\t", index=False)
    architecture_summary.to_csv(outdir / f"domain_cluster_architecture_summary_{suffix}.tsv", sep="\t", index=False)
    domain_class_summary.to_csv(outdir / f"domain_cluster_architecture_domain_class_summary_{suffix}.tsv", sep="\t", index=False)
    return protein_profiles, architecture_summary, domain_class_summary


def plot_summary(summary: pd.DataFrame, domain_class_summary: pd.DataFrame, outdir: Path, cluster_threshold: int) -> None:
    apply_paper_style()
    plot_df = summary.dropna(subset=["n_proteins"]).copy()
    plot_df["architecture_class"] = pd.Categorical(plot_df["architecture_class"], categories=ARCHITECTURE_ORDER, ordered=True)
    plot_df = plot_df.sort_values("architecture_class")

    colors = {
        "broad across IDR": "#CBB2B9",
        "domain-proximal clustered": BREWER_COLORS["mid_gray"],
        "domain-distal clustered": BREWER_COLORS["dark_gray"],
    }

    fig, ax = plt.subplots(figsize=(7.1, 4.2))
    ax.bar(
        plot_df["architecture_class"],
        plot_df["n_proteins"],
        color=[colors[item] for item in plot_df["architecture_class"].astype(str)],
        edgecolor="white",
        linewidth=0.8,
    )
    style_axis(ax)
    ax.set_ylabel("Proteins with clustered methyl-sites")
    ax.set_xlabel("")
    ax.set_title(f"Clustered methylarginine architectures ({cluster_threshold}+ sites / 20 aa)")
    ax.tick_params(axis="x", rotation=18)
    for xpos, value in enumerate(plot_df["n_proteins"]):
        ax.text(xpos, value + 3, str(int(value)), ha="center", va="bottom", fontsize=9, color=BREWER_COLORS["dark_gray"])
    fig.text(
        0.99,
        0.01,
        "Broad = at least one proximal and one distal clustered region; proximal <=20 aa from nearest domain edge; distal >20 aa.",
        ha="right",
        va="bottom",
        fontsize=7.8,
        color=BREWER_COLORS["dark_gray"],
    )
    fig.subplots_adjust(left=0.11, right=0.98, top=0.86, bottom=0.25)
    suffix = f"ge{cluster_threshold}"
    save_plot(fig, outdir / f"domain_cluster_architecture_summary_{suffix}")

    top_classes = (
        domain_class_summary.groupby("nearest_domain_class")["n_proteins"].sum().sort_values(ascending=False).head(8).index
    )
    heat = domain_class_summary[domain_class_summary["nearest_domain_class"].isin(top_classes)].pivot_table(
        index="nearest_domain_class",
        columns="architecture_class",
        values="n_proteins",
        fill_value=0,
        observed=True,
    )
    heat = heat.reindex(columns=ARCHITECTURE_ORDER, fill_value=0)
    heat = heat.loc[heat.sum(axis=1).sort_values(ascending=True).index]

    fig2, ax2 = plt.subplots(figsize=(7.4, 4.8))
    left = pd.Series(0, index=heat.index, dtype=float)
    for cls in ARCHITECTURE_ORDER:
        ax2.barh(heat.index, heat[cls], left=left, color=colors[cls], edgecolor="white", linewidth=0.8, label=cls)
        left += heat[cls]
    style_axis(ax2)
    ax2.set_xlabel("Clustered proteins")
    ax2.set_ylabel("Nearest domain class")
    ax2.set_title("Protein classes carrying clustered methylarginine architectures")
    ax2.legend(frameon=False, loc="lower right", fontsize=8)
    fig2.subplots_adjust(left=0.28, right=0.98, top=0.88, bottom=0.13)
    save_plot(fig2, outdir / f"domain_cluster_architecture_domain_classes_{suffix}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Classify clustered methylarginine regions into broad, domain-proximal, and domain-distal architectures.")
    parser.add_argument("--clusters", default="PTM_results/domain_context/domain_edge_class_20aa_cluster_examples.tsv")
    parser.add_argument("--site-context", default="")
    parser.add_argument("--outdir", default="PTM_results/domain_context")
    parser.add_argument("--cluster-threshold", type=int, default=3)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cluster_path = Path(args.clusters)
    if not cluster_path.exists():
        if not args.site_context:
            raise FileNotFoundError(f"{cluster_path} not found; pass --site-context to derive cluster examples")
        cluster_path = build_cluster_examples(Path(args.site_context), outdir)
    profiles, summary, domain_summary = summarize(cluster_path, outdir, args.cluster_threshold)
    plot_summary(summary, domain_summary, outdir, args.cluster_threshold)
    print(outdir / f"domain_cluster_architecture_classes_ge{args.cluster_threshold}.tsv")
    print(outdir / f"domain_cluster_architecture_summary_ge{args.cluster_threshold}.tsv")
    print(outdir / f"domain_cluster_architecture_summary_ge{args.cluster_threshold}.svg")


if __name__ == "__main__":
    main()

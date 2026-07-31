from __future__ import annotations

import argparse
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import apply_paper_style, save_figure, save_table, style_axis


ONTOLOGY_COLORS = {
    'BP': '#739b7d',
    'CC': '#8da0cb',
    'MF': '#c7a15a',
}


def wrap_label(text: str, width: int = 38) -> str:
    return '\n'.join(textwrap.wrap(str(text), width=width)) if text else ''


def plot_one(path: Path, top_n: int, padj_cutoff: float) -> None:
    df = pd.read_csv(path, sep='\t', low_memory=False)
    if df.empty or 'ONTOLOGY' not in df.columns or 'p.adjust' not in df.columns:
        return

    df = df[df['p.adjust'].notna()].copy()
    df = df[df['p.adjust'] <= padj_cutoff].copy()
    if df.empty:
        return

    plot_df = df.sort_values(['p.adjust', 'Count'], ascending=[True, False]).head(top_n).copy()
    plot_df['neg_log10_padj'] = -np.log10(plot_df['p.adjust'].clip(lower=1e-300))
    plot_df = plot_df.sort_values(['neg_log10_padj', 'Count'], ascending=[True, False]).copy()
    plot_df['plot_label'] = [wrap_label(desc, width=34) for desc in plot_df['Description'].astype(str)]
    plot_df['plot_label'] = pd.Categorical(plot_df['plot_label'], categories=plot_df['plot_label'][::-1], ordered=True)

    fig_height = max(4.8, 0.58 * len(plot_df) + 1.8)
    fig, ax = plt.subplots(figsize=(6.3, fig_height))
    ax.scatter(
        plot_df['neg_log10_padj'],
        plot_df['plot_label'],
        s=34 + plot_df['Count'].astype(float) * 1.2,
        c=[ONTOLOGY_COLORS.get(v, '#6b6b6b') for v in plot_df['ONTOLOGY']],
        alpha=0.95,
        edgecolors='white',
        linewidths=0.65,
    )
    style_axis(ax)
    ax.set_xlabel('-log10(BH adjusted p)', fontsize=11.5, fontfamily='Arial')
    ax.set_ylabel('')
    title = path.name.replace('_go_combined.tsv', '').replace('_', ' ')
    ax.set_title(f'clusterProfiler GO enrichment: {title}', fontsize=13.0, fontfamily='Arial')
    ax.tick_params(axis='both', labelsize=10.8)
    xmax = float(plot_df['neg_log10_padj'].max()) if len(plot_df) else 1.0
    ax.set_xlim(0, xmax * 1.05)

    legend_handles = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markeredgecolor='white', markersize=8, label=ont)
        for ont, color in ONTOLOGY_COLORS.items()
        if ont in set(plot_df['ONTOLOGY'].astype(str))
    ]
    if legend_handles:
        ax.legend(handles=legend_handles, title='Ontology', loc='lower right', frameon=False)

    fig.tight_layout()
    out_prefix = path.with_name(path.name.replace('_go_combined.tsv', '_go_dotplot'))
    save_figure(fig, out_prefix)
    save_table(plot_df.drop(columns=['plot_label']), path.with_name(path.name.replace('_go_combined.tsv', '_go_combined_simplified.tsv')))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--indir', required=True)
    ap.add_argument('--top-n', type=int, default=10)
    ap.add_argument('--padj-cutoff', type=float, default=0.05)
    args = ap.parse_args()

    apply_paper_style()
    indir = Path(args.indir)
    for path in sorted(indir.glob('*_go_combined.tsv')):
        plot_one(path, top_n=args.top_n, padj_cutoff=args.padj_cutoff)


if __name__ == '__main__':
    main()

#!/usr/bin/env python
"""Figures for the methylarginine disease-gene tests (scripts 36-38).

1. disease_ontology_bar   : top 10 Disease Ontology terms (-log10 p), motor neuron disease and ALS highlighted.
2. kegg_simple_dotplot    : simplified KEGG over-representation plot (fold enrichment, one colour, ALS highlighted).
3. als_vs_tdp43_forest    : adjusted odds ratios, inherited ALS genes vs TDP-43 target sets.
Each figure is written as SVG, PDF and 300-dpi PNG, with the exact plotting table as TSV.
"""
from __future__ import annotations

import argparse
import os
import textwrap
from pathlib import Path

import matplotlib
import matplotlib.ticker
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = Path(os.environ.get('PTM_RESULTS_ROOT', ROOT / 'PTM_results'))

NEUTRAL = '#687388'
HIGHLIGHT = '#B5332A'
INK = '#1F2937'
MUTED = '#6B7280'
GRID = '#E5E7EB'
TITLE = '#1F3F7A'
HIGHLIGHT_DO = {'motor neuron disease', 'amyotrophic lateral sclerosis'}


def style() -> None:
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.linewidth': 0.8,
        'axes.edgecolor': MUTED, 'xtick.color': MUTED, 'ytick.color': INK,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    })


def save(fig, prefix: Path) -> None:
    for ext in ('svg', 'pdf', 'png'):
        fig.savefig(prefix.with_suffix(f'.{ext}'), dpi=300, facecolor='white')
    plt.close(fig)


def short(label: str, width: int = 34) -> str:
    label = label[0].upper() + label[1:]
    return '\n'.join(textwrap.wrap(label, width))


def disease_ontology_bar(d: Path, genes: pd.DataFrame, out: Path) -> None:
    do = pd.read_csv(d / 'dose_DO_ms_detectable_bg_full.tsv', sep='\t')
    n_terms = len(do)
    top = do.head(10).copy()
    extra = do[do['Description'].isin(HIGHLIGHT_DO) & (do['rank'] > 10)]
    plot = pd.concat([top, extra]).reset_index(drop=True)
    plot['neglog10p'] = -np.log10(plot['pvalue'])
    plot['k'] = plot['GeneRatio'].str.split('/').str[0].astype(int)
    plot['term_size'] = plot['BgRatio'].str.split('/').str[0].astype(int)
    plot['highlight'] = plot['Description'].isin(HIGHLIGHT_DO)
    plot.to_csv(out / 'fig_disease_ontology_bar.tsv', sep='\t', index=False,
                columns=['rank', 'ID', 'Description', 'k', 'term_size', 'GeneRatio', 'BgRatio',
                         'FoldEnrichment', 'pvalue', 'p.adjust', 'neglog10p', 'highlight', 'geneSymbol'])

    mnd = do[do['Description'] == 'motor neuron disease'].iloc[0]
    als = do[do['Description'] == 'amyotrophic lateral sclerosis'].iloc[0]
    mnd_genes = genes[genes['symbol'].isin(mnd['geneSymbol'].split('/'))].sort_values(
        ['arg_methyl_sites', 'symbol'], ascending=[False, True])
    shown = list(mnd_genes['symbol'].head(16))

    fig = plt.figure(figsize=(11, 6.2))
    fig.text(0.04, 0.94, f"Motor neuron disease ranks {int(mnd['rank'])} of {n_terms} disease terms",
             fontsize=19, fontweight='bold', color=TITLE)
    fig.text(0.04, 0.885, f"Disease Ontology enrichment (DOSE) of {genes['methyl_any'].sum():,} methylarginine "
             f"genes against an MS-detectable background. Top 10 terms, plus ALS.", fontsize=10.5, color=MUTED)
    ax = fig.add_axes([0.30, 0.20, 0.36, 0.62])
    y = np.arange(len(plot))[::-1].astype(float)
    if len(extra):
        y[10:] -= 0.6  # visual gap before rows outside the top 10
    colors = [HIGHLIGHT if h else NEUTRAL for h in plot['highlight']]
    ax.barh(y, plot['neglog10p'], color=colors, height=0.68)
    ax.set_yticks(y)
    ax.set_yticklabels([short(s) for s in plot['Description']], fontsize=9)
    for yi, r in zip(y, plot.itertuples()):
        if r.highlight:
            ax.text(r.neglog10p + 0.08, yi, f'{r.k} of {r.term_size} genes', va='center', fontsize=9,
                    fontweight='bold', color=INK)
        if r.rank > 10:
            ax.text(-0.05, yi + 0.55, f'rank {r.rank}', ha='right', va='center', fontsize=7.5,
                    color=MUTED, transform=ax.get_yaxis_transform(), style='italic')
    xmax = float(np.ceil(plot['neglog10p'].max() + 1.2))
    ax.set_xlim(0, xmax)
    ax.set_xlabel('Enrichment, −log10 p', color=MUTED)
    ax.xaxis.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis='y', length=0)

    fig.text(0.70, 0.79, 'Methylated motor neuron\ndisease genes', fontsize=11.5, fontweight='bold',
             color=HIGHLIGHT, va='top')
    lines = textwrap.wrap('  '.join(shown), 30)
    fig.text(0.70, 0.70, '\n'.join(lines), fontsize=10.5, fontweight='bold', color=INK, va='top', linespacing=1.35)
    fig.text(0.70, 0.70 - 0.043 * len(lines) - 0.012, f'+ {len(mnd_genes) - len(shown)} more (ordered by site count)',
             fontsize=9.5, color=MUTED, va='top')
    bg_n = int(mnd['BgRatio'].split('/')[1])
    fig.text(0.70, 0.46, textwrap.fill(
        f"Fold enrichment {mnd['FoldEnrichment']:.2f}. Background: {bg_n:,} DO-annotated genes "
        f"detected by MS in human cell lines (PaxDb 5.0) or carrying a methylarginine site, "
        f"so only proteins that could have been seen are counted.", 46), fontsize=8.8, color=MUTED, va='top')

    fig.text(0.04, 0.085,
             f"Motor neuron disease: rank {int(mnd['rank'])} of {n_terms} (p = {mnd['pvalue']:.1g}; FDR {mnd['p.adjust']:.2f}). "
             f"ALS: rank {int(als['rank'])} (p = {als['pvalue']:.1g}; FDR {als['p.adjust']:.2f}). "
             f"Top terms are developmental and intellectual-disability terms.", fontsize=8.8, color=MUTED)
    fig.text(0.04, 0.04, textwrap.fill(
        'DOSE 4.6 enrichDO (BH). Methylarginine: PTM_pipeline site union (Maron 2021 Table S5, PRoMetheusDB, '
        'curated PTM databases; any evidence). Background: PaxDb 5.0 cell-line proteomes (Geiger 2012; Nagaraj 2011; '
        'Kulak 2014; Müller 2020).', 170), fontsize=7.8, color=MUTED, style='italic')
    save(fig, out / 'fig_disease_ontology_bar')


def kegg_simple(d: Path, out: Path, n: int = 15) -> None:
    kg = pd.read_csv(d / 'dose_KEGG_ms_detectable_bg_full.tsv', sep='\t').head(n).copy()
    kg['k'] = kg['GeneRatio'].str.split('/').str[0].astype(int)
    kg['highlight'] = kg['Description'].str.contains('lateral sclerosis|neurodegeneration', case=False)
    kg.to_csv(out / 'fig_kegg_simple_dotplot.tsv', sep='\t', index=False,
              columns=['rank', 'ID', 'Description', 'k', 'GeneRatio', 'BgRatio', 'FoldEnrichment', 'pvalue', 'p.adjust', 'highlight'])

    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    fig.subplots_adjust(left=0.42, right=0.78, top=0.86, bottom=0.15)
    y = np.arange(len(kg))[::-1]
    colors = [HIGHLIGHT if h else NEUTRAL for h in kg['highlight']]
    ax.hlines(y, 1, kg['FoldEnrichment'], color=GRID, lw=1.6, zorder=1)
    ax.scatter(kg['FoldEnrichment'], y, s=46, color=colors, zorder=2, edgecolor='white', linewidth=1)
    ax.axvline(1, color=MUTED, lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(kg['Description'], fontsize=8.6)
    for lab, h in zip(ax.get_yticklabels(), kg['highlight']):
        if h:
            lab.set_color(HIGHLIGHT)
            lab.set_fontweight('bold')
    xmax = kg['FoldEnrichment'].max() * 1.08
    ax.set_xlim(0.9, xmax)
    for yi, k, fdr in zip(y, kg['k'], kg['p.adjust']):
        ax.annotate(f'{k}', (1.03, yi), xycoords=('axes fraction', 'data'), va='center', fontsize=8, color=INK)
        ax.annotate(f'{fdr:.0e}'.replace('e-0', 'e-'), (1.22, yi), xycoords=('axes fraction', 'data'),
                    va='center', fontsize=8, color=INK)
    ax.annotate('Genes', (1.03, len(kg) - 0.2), xycoords=('axes fraction', 'data'), fontsize=8, color=MUTED, fontweight='bold')
    ax.annotate('FDR', (1.22, len(kg) - 0.2), xycoords=('axes fraction', 'data'), fontsize=8, color=MUTED, fontweight='bold')
    ax.set_xlabel('Fold enrichment', color=MUTED)
    ax.xaxis.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis='y', length=0)
    fig.text(0.02, 0.95, 'KEGG pathways enriched among methylarginine proteins', fontsize=12, fontweight='bold', color=TITLE)
    fig.text(0.02, 0.905, f'Top {n} of {len(pd.read_csv(d / "dose_KEGG_ms_detectable_bg_full.tsv", sep=chr(9)))} pathways '
             'by p value; MS-detectable background', fontsize=8.8, color=MUTED)
    fig.text(0.02, 0.025, 'clusterProfiler enrichKEGG, BH-adjusted. Same gene set and background as the Disease Ontology panel.',
             fontsize=7.4, color=MUTED, style='italic')
    save(fig, out / 'fig_kegg_simple_dotplot')


def forest(d: Path, out: Path) -> None:
    t = pd.read_csv(d / 'tdp43_compendium_methylarginine.tsv', sep='\t')
    t = t[(t['background'] == 'MS-detectable') & (t['outcome'] == 'methyl_any')].copy()
    order = ['Inherited ALS genes (ClinGen ALS panel)',
             'All compendium host genes (any TDP-43-associated junction)',
             'TDP-43-sensitive cryptic events', 'NMD-rescued or dark cryptic events',
             'Down after TDP-43 knockdown (>=2 datasets)', 'Pathology-discriminating events (AUC >= 0.8)']
    t = t.set_index('gene_set').loc[order].reset_index()
    t.to_csv(out / 'fig_als_vs_tdp43_forest.tsv', sep='\t', index=False)

    fig, ax = plt.subplots(figsize=(7.4, 3.4))
    fig.subplots_adjust(left=0.40, right=0.74, top=0.78, bottom=0.19)
    y = np.arange(len(t))[::-1].astype(float)
    y[1:] -= 0.5
    colors = [HIGHLIGHT] + [NEUTRAL] * (len(t) - 1)
    ax.hlines(y, t['adj_OR_lo95'], t['adj_OR_hi95'], color=colors, lw=1.8)
    ax.scatter(t['adj_OR'], y, color=colors, s=44, zorder=3, edgecolor='white', linewidth=1)
    ax.axvline(1, color=MUTED, lw=0.8, ls='--')
    ax.set_xscale('log')
    ax.set_xticks([0.5, 1, 2, 4, 8])
    ax.set_xticklabels(['0.5', '1', '2', '4', '8'])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlim(0.4, 9)
    labels = ['Inherited ALS genes (ClinGen)', 'All TDP-43 compendium genes', 'TDP-43-sensitive cryptic events',
              'NMD-rescued / dark cryptic events', 'Down after TDP-43 knockdown', 'Pathology-discriminating (AUC ≥ 0.8)']
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.6)
    ax.get_yticklabels()[0].set_color(HIGHLIGHT)
    ax.get_yticklabels()[0].set_fontweight('bold')
    for yi, r in zip(y, t.itertuples()):
        ax.annotate(f'{r.n_methyl}/{r.n_set}', (1.03, yi), xycoords=('axes fraction', 'data'), va='center', fontsize=8)
        ax.annotate(f'{r.adj_p:.2g}', (1.28, yi), xycoords=('axes fraction', 'data'), va='center', fontsize=8)
    ax.annotate('Methyl', (1.03, y[0] + 0.6), xycoords=('axes fraction', 'data'), fontsize=8, color=MUTED, fontweight='bold')
    ax.annotate('p', (1.28, y[0] + 0.6), xycoords=('axes fraction', 'data'), fontsize=8, color=MUTED, fontweight='bold')
    ax.set_xlabel('Adjusted odds ratio of carrying methylarginine (95% CI)', color=MUTED)
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis='y', length=0)
    fig.text(0.02, 0.93, 'ALS genes carry methylarginine; most TDP-43 target sets do not', fontsize=12,
             fontweight='bold', color=TITLE)
    fig.text(0.02, 0.865, 'Logistic regression adjusted for RNA binding, protein length and abundance; MS-detectable background',
             fontsize=8.3, color=MUTED)
    fig.text(0.02, 0.03, 'TDP-43 sets: external cryptic-junction compendium (Shechterlab/cryptic_exon). Unadjusted for multiple sets.',
             fontsize=7.2, color=MUTED, style='italic')
    save(fig, out / 'fig_als_vs_tdp43_forest')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--dir', type=Path, default=RESULTS / 'disease_gene_tests')
    args = ap.parse_args()
    style()
    genes = pd.read_csv(args.dir / 'gene_table.tsv', sep='\t')
    genes['methyl_any'] = genes['methyl_any'].astype(str).eq('True')
    out = args.dir / 'figures'
    out.mkdir(exist_ok=True)
    disease_ontology_bar(args.dir, genes, out)
    kegg_simple(args.dir, out)
    forest(args.dir, out)
    print('wrote', sorted(p.name for p in out.iterdir()))


if __name__ == '__main__':
    main()

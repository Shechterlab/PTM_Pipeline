#!/usr/bin/env python
"""16:9 slide figures for an ALS talk, built from the outputs of scripts 36-40.

1. als_kegg_slide        : KEGG enrichment of arginine-methylated proteins (2022 analysis re-run on the union,
                           KEGG default background), top pathways by p value, ALS/neurodegeneration highlighted.
2. als_clingen_slide     : fraction of ALS-causing proteins that are arginine-methylated vs all disease proteins
                           and all human proteins (ClinGen), with the ranked panel strip.
3. als_protein_map_slide : methylarginine sites on FET proteins, TDP-43, hnRNPA1/A2B1 and MATR3, drawn on UniProt
                           domain annotations.
4. als_do_slide          : Disease Ontology (DOSE) top terms, motor neuron disease and ALS highlighted.
Each is written as SVG, PDF and 300-dpi PNG with its plotting table.
"""
from __future__ import annotations

import json
import os
import textwrap
from pathlib import Path

import matplotlib
import matplotlib.ticker
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get('PTM_DATA_ROOT', ROOT / 'PTM_data'))
RESULTS = Path(os.environ.get('PTM_RESULTS_ROOT', ROOT / 'PTM_results'))
D = RESULTS / 'disease_gene_tests'
OUT = D / 'als_talk_figures'

NEUTRAL = '#687388'
LIGHT = '#C9CED8'
HIGHLIGHT = '#B5332A'
HIGHLIGHT_SOFT = '#E9B8B3'
INK = '#1F2937'
MUTED = '#6B7280'
GRID = '#E5E7EB'
TITLE = '#1F3F7A'
DOMAIN = '#3F5E8C'
LCD = '#EFE6D2'
SLIDE = (13.33, 7.5)

ALS_KEGG = ('Amyotrophic lateral sclerosis', 'Pathways of neurodegeneration - multiple diseases')
PROTEINS = [  # (gene, accession, label)
    ('FUS', 'P35637', 'FUS'), ('EWSR1', 'Q01844', 'EWSR1'), ('TAF15', 'Q92804', 'TAF15'),
    ('TARDBP', 'Q13148', 'TDP-43'), ('HNRNPA1', 'P09651', 'hnRNPA1'), ('HNRNPA2B1', 'P22626', 'hnRNPA2B1'),
    ('MATR3', 'P43243', 'MATR3'),
]


def style() -> None:
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 12, 'axes.linewidth': 0.9,
        'axes.edgecolor': MUTED, 'xtick.color': MUTED, 'ytick.color': INK,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    })


def slide(title: str, subtitle: str):
    fig = plt.figure(figsize=SLIDE)
    fig.text(0.045, 0.915, title, fontsize=23, fontweight='bold', color=TITLE)
    fig.text(0.045, 0.862, subtitle, fontsize=13, color=MUTED)
    return fig


def footer(fig, text: str) -> None:
    fig.text(0.045, 0.035, textwrap.fill(text, 175), fontsize=9.5, color=MUTED, style='italic')


def save(fig, name: str) -> None:
    for ext in ('svg', 'pdf', 'png'):
        fig.savefig(OUT / f'{name}.{ext}', dpi=300, facecolor='white')
    plt.close(fig)


def clean(ax) -> None:
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis='y', length=0)


def kegg_slide(n: int = 12) -> None:
    kg = pd.read_csv(D / 'kegg_2022_rerun' / 'enrichKEGG_ptm_pipeline_union.tsv', sep='\t')
    kg = kg.sort_values('pvalue').reset_index(drop=True)
    kg['rank'] = np.arange(1, len(kg) + 1)
    top = kg.head(n).copy()
    top['k'] = top['GeneRatio'].str.split('/').str[0].astype(int)
    top['neglog10fdr'] = -np.log10(top['p.adjust'])
    top['highlight'] = top['Description'].isin(ALS_KEGG)
    top.to_csv(OUT / 'als_kegg_slide.tsv', sep='\t', index=False)
    als = kg[kg['Description'] == ALS_KEGG[0]].iloc[0]
    n_query = int(als['GeneRatio'].split('/')[1])

    fig = slide(f"ALS is the #{int(als['rank'])} KEGG pathway for arginine-methylated proteins",
                f"KEGG enrichment of {n_query:,} KEGG-annotated methylated proteins. Top {n} of {len(kg)} significant pathways.")
    ax = fig.add_axes([0.33, 0.14, 0.40, 0.66])
    y = np.arange(len(top))[::-1]
    ax.barh(y, top['neglog10fdr'], color=[HIGHLIGHT if h else NEUTRAL for h in top['highlight']], height=0.66)
    ax.set_yticks(y)
    ax.set_yticklabels(top['Description'].str.replace(' - multiple diseases', ''), fontsize=12.5)
    for lab, h in zip(ax.get_yticklabels(), top['highlight']):
        if h:
            lab.set_color(HIGHLIGHT)
            lab.set_fontweight('bold')
    for yi, r in zip(y, top.itertuples()):
        ax.text(r.neglog10fdr + 0.2, yi, f'{r.k}', va='center', fontsize=11,
                color=INK, fontweight='bold' if r.highlight else 'normal')
    ax.set_xlim(0, top['neglog10fdr'].max() * 1.12)
    ax.set_xlabel('Enrichment, −log10 FDR   (numbers: methylated proteins in pathway)', color=MUTED)
    ax.xaxis.grid(True, color=GRID)
    ax.set_axisbelow(True)
    clean(ax)

    genes = set(als['geneID'].split('/'))
    showcase = [g for g in ['TARDBP', 'FUS', 'TAF15', 'EWSR1', 'HNRNPA1', 'HNRNPA2B1', 'MATR3', 'ATXN2', 'VCP',
                            'SOD1', 'UBQLN2', 'SQSTM1', 'G3BP1', 'NEFH', 'PRPH', 'DCTN1', 'SETX', 'TUBA4A'] if g in genes]
    fig.text(0.765, 0.78, 'Methylated proteins in\nthe KEGG ALS pathway', fontsize=14, fontweight='bold',
             color=HIGHLIGHT, va='top')
    lines = textwrap.wrap('  '.join(showcase), 22)
    fig.text(0.765, 0.67, '\n'.join(lines), fontsize=13, fontweight='bold', color=INK, va='top', linespacing=1.4)
    fig.text(0.765, 0.67 - 0.038 * len(lines) - 0.01, f"+ {int(als['Count']) - len(showcase)} more",
             fontsize=12, color=MUTED, va='top')
    fig.text(0.765, 0.30, textwrap.fill(
        f"{int(als['Count'])} of {int(als['BgRatio'].split('/')[0])} pathway proteins are methylated "
        f"({als['FoldEnrichment']:.1f}-fold over expected; FDR {als['p.adjust']:.0e}).", 34),
        fontsize=12, color=INK, va='top')
    footer(fig, 'clusterProfiler enrichKEGG (organism hsa, default KEGG background, BH), as in the 2022 analysis of Maron et al. '
                '2021 Table S5, re-run on the PTM_pipeline methylarginine union (Table S5, PRoMetheusDB, curated PTM databases). '
                'Against a stricter MS-detectable background the ALS pathway ranks 11 (FDR 0.001).')
    save(fig, 'als_kegg_slide')


def clingen_slide() -> None:
    panels = pd.read_csv(D / 'clingen_panel_methylarginine.tsv', sep='\t')
    genes = pd.read_csv(D / 'gene_table.tsv', sep='\t')
    genes['methyl_any'] = genes['methyl_any'].astype(str).eq('True')
    als = panels[panels['panel'].str.startswith('Amyotrophic')].iloc[0]
    # Fraction across all curated disease proteins (each gene once), from the regression background.
    cur = pd.read_csv(D / 'clingen_curated_genes.tsv', sep='\t')
    disease_frac = cur['methyl_any'].mean()
    proteome_frac = genes['methyl_any'].mean()
    bars = pd.DataFrame({
        'group': ['ALS-causing proteins', 'All disease-causing proteins', 'All human proteins'],
        'fraction': [als['frac_methyl'], disease_frac, proteome_frac],
        'n': [f"{int(als['n_methyl'])} of {int(als['n_set'])}", f"{int(cur['methyl_any'].sum()):,} of {len(cur):,}",
              f"{int(genes['methyl_any'].sum()):,} of {len(genes):,}"],
    })
    bars.to_csv(OUT / 'als_clingen_slide.tsv', sep='\t', index=False)

    fig = slide(f"{als['frac_methyl']:.0%} of ALS-causing proteins are arginine-methylated",
                f"Most of any disease group: ALS ranks 1st of {len(panels)} ClinGen curation panels.")
    ax = fig.add_axes([0.06, 0.24, 0.46, 0.54])
    x = np.arange(3)
    ax.bar(x, bars['fraction'] * 100, color=[HIGHLIGHT, NEUTRAL, LIGHT], width=0.62)
    for xi, r in zip(x, bars.itertuples()):
        ax.text(xi, r.fraction * 100 + 2, f'{r.fraction:.0%}', ha='center', fontsize=22, fontweight='bold',
                color=HIGHLIGHT if xi == 0 else INK)
    ax.set_xticks(x)
    ax.set_xticklabels([f'{lab}\n({n_})' for lab, n_ in zip(['ALS-causing proteins', 'All disease proteins',
                                                              'All human proteins'], bars['n'])], fontsize=12.5)
    ax.tick_params(axis='x', length=0, pad=8)
    ax.set_ylim(0, 85)
    ax.set_yticks([])
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)

    ax2 = fig.add_axes([0.62, 0.24, 0.32, 0.54])
    p = panels.sort_values('frac_methyl', ascending=True).reset_index(drop=True)
    isals = p['panel'].str.startswith('Amyotrophic')
    ax2.barh(np.arange(len(p)), p['frac_methyl'] * 100, color=np.where(isals, HIGHLIGHT, LIGHT), height=0.8)
    ax2.axvline(disease_frac * 100, color=NEUTRAL, lw=1, ls='--')
    ax2.text(disease_frac * 100 + 1, -2.3, 'average disease protein', fontsize=10, color=NEUTRAL)
    i = int(np.flatnonzero(isals)[0])
    ax2.text(p.loc[i, 'frac_methyl'] * 100 + 1.5, i, 'ALS', va='center', fontsize=12,
             color=HIGHLIGHT, fontweight='bold')
    ax2.set_yticks([])
    ax2.set_xlim(0, 100)
    ax2.set_ylim(-3, len(p))
    ax2.set_xlabel('% of panel proteins methylated', color=MUTED)
    ax2.set_title(f'All {len(p)} ClinGen disease panels', fontsize=13, color=INK, loc='left')
    clean(ax2)
    footer(fig, f"ClinGen gene–disease validity (Definitive/Strong/Moderate; panels ≥10 proteins). Methylated = any arginine-methylation "
                f"site in the PTM_pipeline union. Adjusted for RNA binding, length and abundance, ALS odds ratio "
                f"{als['adj_OR']:.1f} (95% CI {als['adj_OR_lo95']:.1f}–{als['adj_OR_hi95']:.1f}). Methylated ALS proteins: "
                f"{als['methyl_genes'].replace(',', ', ')}.")
    save(fig, 'als_clingen_slide')


def uniprot_features(acc: str) -> tuple[int, list[tuple[str, int, int, str]]]:
    d = json.load(open(DATA / 'context' / 'uniprot_features' / f'{acc}.json'))
    feats = [(f['type'], f['location']['start']['value'], f['location']['end']['value'], f.get('description', ''))
             for f in d.get('features', [])]
    return d['sequence']['length'], feats


def protein_map_slide() -> None:
    union = pd.read_csv(RESULTS / 'integrated' / 'human_arg_methyl_union_dedup_by_site.tsv', sep='\t', low_memory=False,
                        usecols=['canonical_UniProtAC', 'corrected_position', 'source_family_count_exact', 'reference_any'])
    rows = []
    fig = slide('ALS proteins carry methylarginine in RGG / low-complexity regions',
                'Arginine-methylation sites on FET proteins, TDP-43, hnRNPA1/A2B1 and MATR3.')
    maxlen = 850
    ax = fig.add_axes([0.14, 0.16, 0.80, 0.61])
    n = len(PROTEINS)
    for i, (gene, acc, label) in enumerate(PROTEINS):
        yy = n - 1 - i
        length, feats = uniprot_features(acc)
        sites = union[union['canonical_UniProtAC'] == acc].copy()
        pos = sorted(sites['corrected_position'].astype(int))
        refs = sites['reference_any'].fillna('').str.split(';').explode().str.strip()
        n_studies = refs[refs != ''].nunique()
        for p_ in pos:
            rows.append({'protein': label, 'uniprot': acc, 'position': p_})
        ax.add_patch(Rectangle((0, yy - 0.06), length, 0.12, color='#B8BEC9', zorder=1))
        for ftype, s, e, desc in feats:
            if ftype == 'Compositional bias' and ('Gly' in desc or 'Low complexity' in desc) or \
                    (ftype == 'Region' and ('RGG' in desc or 'D-R-' in desc)):
                ax.add_patch(Rectangle((s, yy - 0.17), e - s, 0.34, color=LCD, zorder=1.5, lw=0))
        for ftype, s, e, desc in feats:
            if ftype in ('Domain', 'Zinc finger') and desc not in ('IQ',):
                ax.add_patch(FancyBboxPatch((s, yy - 0.2), e - s, 0.4, boxstyle='round,pad=0,rounding_size=6',
                                            color=DOMAIN, zorder=2, lw=0, mutation_aspect=0.02))
                short = 'RRM' if desc.startswith('RRM') else ('ZnF' if ftype == 'Zinc finger' else desc)
                ax.text((s + e) / 2, yy, short, ha='center', va='center', fontsize=8.5, color='white',
                        fontweight='bold', zorder=3)
        for p_ in pos:
            ax.plot([p_, p_], [yy + 0.2, yy + 0.36], color=HIGHLIGHT, lw=1.1, zorder=3)
            ax.scatter([p_], [yy + 0.38], s=14, color=HIGHLIGHT, zorder=4, lw=0)
        ax.text(-20, yy, label, ha='right', va='center', fontsize=15, fontweight='bold', color=INK)
        ax.text(length + 15, yy, f'{len(pos)} sites', va='center', fontsize=11.5, color=HIGHLIGHT, fontweight='bold')
        if gene == 'TARDBP':
            for p_, ha, dx in zip(pos, ['right', 'left'], [-3, 3]):
                ax.annotate(f'R{p_}', (p_, yy + 0.42), xytext=(dx, 4), textcoords='offset points', ha=ha,
                            fontsize=10, color=HIGHLIGHT, fontweight='bold')
    ax.set_xlim(0, maxlen + 90)
    ax.set_ylim(-0.6, n - 0.35)
    ax.set_yticks([])
    ax.set_xlabel('Residue', color=MUTED)
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)
    # legend
    lx, ly = 0.60, 0.80
    fig.patches.append(Rectangle((lx, ly), 0.018, 0.022, transform=fig.transFigure, color=DOMAIN))
    fig.text(lx + 0.024, ly + 0.004, 'Folded domain (RRM, ZnF)', fontsize=10.5, color=INK)
    fig.patches.append(Rectangle((lx + 0.20, ly), 0.018, 0.022, transform=fig.transFigure, color=LCD))
    fig.text(lx + 0.224, ly + 0.004, 'Gly-rich / RGG / LCD', fontsize=10.5, color=INK)
    pd.DataFrame(rows).to_csv(OUT / 'als_protein_map_slide.tsv', sep='\t', index=False)
    footer(fig, 'Sites: PTM_pipeline methylarginine union, canonical UniProt positions. Domains and Gly-rich / low-complexity '
                'regions: UniProt annotations. TDP-43 R293 (C-terminal Gly-rich LCD) is reported as monomethylated in 7 '
                'proteomic studies; R275 in 3.')
    save(fig, 'als_protein_map_slide')


def do_slide() -> None:
    do = pd.read_csv(D / 'dose_DO_reviewed_bg_full.tsv', sep='\t')
    top = do.head(10)
    extra = do[do['Description'].isin(['motor neuron disease', 'amyotrophic lateral sclerosis']) & (do['rank'] > 10)]
    plot = pd.concat([top, extra]).reset_index(drop=True)
    plot['k'] = plot['GeneRatio'].str.split('/').str[0].astype(int)
    plot['size'] = plot['BgRatio'].str.split('/').str[0].astype(int)
    plot['neglog10fdr'] = -np.log10(plot['p.adjust'])
    plot['highlight'] = plot['Description'].isin(['motor neuron disease', 'amyotrophic lateral sclerosis'])
    plot.to_csv(OUT / 'als_do_slide.tsv', sep='\t', index=False)
    mnd = do[do['Description'] == 'motor neuron disease'].iloc[0]
    n_query = int(mnd['GeneRatio'].split('/')[1])

    fig = slide(f"Motor neuron disease ranks #{int(mnd['rank'])} of {len(do):,} disease terms",
                f"Disease Ontology enrichment of {n_query:,} arginine-methylated proteins (DO-annotated). "
                f"Top 10 of {len(do):,} terms, plus ALS.")
    ax = fig.add_axes([0.33, 0.14, 0.42, 0.66])
    y = np.arange(len(plot))[::-1].astype(float)
    if len(extra):
        y[10:] -= 0.6
    ax.barh(y, plot['neglog10fdr'], color=[HIGHLIGHT if h else NEUTRAL for h in plot['highlight']], height=0.66)
    ax.set_yticks(y)
    ax.set_yticklabels(['\n'.join(textwrap.wrap(s[0].upper() + s[1:], 32)) for s in plot['Description']], fontsize=12)
    for lab, h in zip(ax.get_yticklabels(), plot['highlight']):
        if h:
            lab.set_color(HIGHLIGHT)
            lab.set_fontweight('bold')
    for yi, r in zip(y, plot.itertuples()):
        if r.highlight:
            ax.text(r.neglog10fdr + 0.15, yi, f'{r.k} of {r.size}', va='center', fontsize=12, fontweight='bold', color=INK)
        if r.rank > 10:
            ax.text(-0.01, yi + 0.55, f'rank {r.rank}', ha='right', va='center', fontsize=10, color=MUTED,
                    style='italic', transform=ax.get_yaxis_transform())
    ax.set_xlim(0, plot['neglog10fdr'].max() * 1.15)
    ax.set_xlabel('Enrichment, −log10 FDR', color=MUTED)
    ax.xaxis.grid(True, color=GRID)
    ax.set_axisbelow(True)
    clean(ax)
    fig.text(0.785, 0.78, 'Methylated motor\nneuron disease\nproteins', fontsize=14, fontweight='bold', color=HIGHLIGHT, va='top')
    shown = ['FUS', 'TAF15', 'EWSR1', 'TARDBP', 'HNRNPA1', 'MATR3', 'VCP', 'SOD1', 'ANXA11', 'NEFM', 'DYNC1H1', 'BAG3', 'HSPB1', 'UBA1']
    shown = [g for g in shown if g in set(mnd['geneSymbol'].split('/'))]
    lines = textwrap.wrap('  '.join(shown), 19)
    fig.text(0.785, 0.67, '\n'.join(lines), fontsize=13, fontweight='bold', color=INK, va='top', linespacing=1.4)
    fig.text(0.785, 0.67 - 0.038 * len(lines) - 0.01, f"+ {int(mnd['Count']) - len(shown)} more", fontsize=12, color=MUTED, va='top')
    footer(fig, f"DOSE enrichDO, reviewed human proteome background (BH). Motor neuron disease FDR {mnd['p.adjust']:.0e}; "
                f"against a stricter MS-detectable background it ranks 8 (FDR 0.06). Counts are DO-annotated proteins.")
    save(fig, 'als_do_slide')


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    style()
    kegg_slide()
    clingen_slide()
    protein_map_slide()
    do_slide()
    print('wrote', sorted(p.name for p in OUT.iterdir() if p.suffix == '.png'))


if __name__ == '__main__':
    main()

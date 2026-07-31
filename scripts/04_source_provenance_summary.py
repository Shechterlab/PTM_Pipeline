from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from common import BREWER_COLORS, annotate_barh, apply_paper_style, save_figure, save_table, style_axis


def protein_label(row: pd.Series) -> str:
    gene = str(row.get('substrate_genename', '') or '').strip()
    acc = str(row.get('substrate_UniProtAC', '') or '').strip()
    suffix = ' *' if bool(row.get('accession_is_isoform', False)) else ''
    if gene:
        return f'{gene} ({acc}){suffix}'
    return f'{acc}{suffix}'


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--integrated-sites', required=True)
    ap.add_argument('--outdir', default='results/qc/source_provenance')
    args = ap.parse_args()

    df = pd.read_csv(args.integrated_sites, sep='\t', low_memory=False)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    apply_paper_style()

    support_dist = (
        df['source_family_count_exact']
        .fillna(1)
        .astype(int)
        .value_counts()
        .rename_axis('source_family_count_exact')
        .reset_index(name='site_count')
        .sort_values('source_family_count_exact')
    )
    save_table(support_dist, outdir / 'source_family_support_distribution.tsv')
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    ax.bar(
        support_dist['source_family_count_exact'].astype(str),
        support_dist['site_count'],
        color=BREWER_COLORS['mid_gray'],
        edgecolor='white',
        linewidth=0.8,
    )
    style_axis(ax)
    ax.set_ylabel('Deduplicated sites')
    ax.set_xlabel('Exact source-family support count')
    ax.set_title('Integrated arg-methylome source-support distribution')
    fig.tight_layout()
    save_figure(fig, outdir / 'arg_methyl_source_support_distribution')

    support = pd.DataFrame({
        'support_type': ['Exact >=2', 'Fuzzy >=2'],
        'site_count': [int(df['support_exact_ge_2'].sum()), int(df['support_any_ge_2'].sum())],
    })
    save_table(support, outdir / 'support_exact_vs_fuzzy_summary.tsv')
    fig2, ax2 = plt.subplots(figsize=(5.6, 4.8))
    ax2.bar(support['support_type'], support['site_count'], color=BREWER_COLORS['mid_gray'], edgecolor='white', linewidth=0.8)
    style_axis(ax2)
    ax2.set_ylabel('Deduplicated sites')
    ax2.set_title('Multi-source support recovery with fuzzy matching')
    for i, value in enumerate(support['site_count']):
        ax2.text(i, value, str(value), ha='center', va='bottom', fontsize=9, color=BREWER_COLORS['dark_gray'])
    fig2.tight_layout()
    save_figure(fig2, outdir / 'arg_methyl_support_exact_vs_fuzzy')

    protein_counts = df.groupby(
        ['substrate_UniProtAC', 'canonical_UniProtAC', 'substrate_genename', 'accession_is_isoform', 'needs_sequence_remap', 'canonical_counting_status'],
        as_index=False,
    ).size().rename(columns={'size': 'site_count_accession_level'})
    protein_counts['label'] = protein_counts.apply(protein_label, axis=1)
    protein_counts = protein_counts.sort_values(
        ['site_count_accession_level', 'substrate_genename', 'substrate_UniProtAC'],
        ascending=[False, True, True],
    )
    save_table(protein_counts, outdir / 'protein_site_counts_accession_level.tsv')

    top = protein_counts.head(20).sort_values('site_count_accession_level')
    save_table(top, outdir / 'top_proteins_by_site_count_accession_level.tsv')
    fig3, ax3 = plt.subplots(figsize=(10.2, 7.0))
    ax3.barh(top['label'], top['site_count_accession_level'], color=BREWER_COLORS['mid_gray'], edgecolor='white', linewidth=0.8)
    style_axis(ax3)
    ax3.set_xlabel('Unique accession-level sites')
    ax3.set_ylabel('Protein')
    ax3.set_title('Top proteins by accession-level arg-methyl site count')
    ax3.set_xlim(0, float(top['site_count_accession_level'].max()) * 1.18)
    annotate_barh(
        ax3,
        top['label'],
        top['site_count_accession_level'],
        [str(v) for v in top['site_count_accession_level']],
        fontsize=8.2,
    )
    fig3.tight_layout()
    save_figure(fig3, outdir / 'arg_methyl_top_proteins_accession_level')


if __name__ == '__main__':
    main()

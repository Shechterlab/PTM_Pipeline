from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from common import save_figure, save_table


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
    ap.add_argument('--outdir', default='results/example_figures')
    args = ap.parse_args()

    df = pd.read_csv(args.integrated_sites, sep='	', low_memory=False)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    support = pd.DataFrame({
        'support_type': ['Exact >=2', 'Fuzzy >=2'],
        'site_count': [int(df['support_exact_ge_2'].sum()), int(df['support_any_ge_2'].sum())],
    })
    save_table(support, outdir / 'support_exact_vs_fuzzy_summary.tsv')
    fig2, ax2 = plt.subplots(figsize=(5.6, 4.8))
    ax2.bar(support['support_type'], support['site_count'])
    ax2.set_ylabel('Deduplicated sites')
    ax2.set_title('Multi-source support recovery with fuzzy matching')
    for i, value in enumerate(support['site_count']):
        ax2.text(i, value, str(value), ha='center', va='bottom')
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
    ax3.barh(top['label'], top['site_count_accession_level'])
    ax3.set_xlabel('Unique accession-level sites')
    ax3.set_ylabel('Protein')
    ax3.set_title('Top proteins by accession-level arg-methyl site count')
    fig3.tight_layout()
    save_figure(fig3, outdir / 'arg_methyl_top_proteins_accession_level')


if __name__ == '__main__':
    main()

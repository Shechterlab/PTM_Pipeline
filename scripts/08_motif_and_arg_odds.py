from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from common import (
    canonical_site_table,
    format_p_value,
    log2_odds_ratio,
    parse_fasta,
    residue_background_from_sequences,
    save_figure,
    save_table,
    sequence_window,
)


HYDROPHOBIC = set('AILMFWVY')


def residue_at(window: str, relative_position: int) -> str:
    center = len(window) // 2
    index = center + relative_position
    if index < 0 or index >= len(window):
        return '_'
    return window[index]


def family_flags(window: str) -> dict[str, bool]:
    neighborhood = [residue_at(window, i) for i in range(-5, 6) if i != 0]
    near3 = [residue_at(window, i) for i in range(-3, 4) if i != 0]
    return {
        'motif_RG': residue_at(window, 1) == 'G',
        'motif_RGG': residue_at(window, 1) == 'G' and residue_at(window, 2) == 'G',
        'motif_GRG': residue_at(window, -1) == 'G' and residue_at(window, 1) == 'G',
        'motif_GAR': residue_at(window, -2) == 'G' and residue_at(window, -1) == 'A',
        'motif_RXR': residue_at(window, -2) == 'R' or residue_at(window, 2) == 'R',
        'motif_CARM1_proline_rich': 'P' in neighborhood,
        'motif_CARM1_hydrophobic_proline': ('P' in neighborhood) and any(aa in HYDROPHOBIC for aa in near3),
        'motif_PRMT5_acidic_DR_like': residue_at(window, -1) in {'D', 'E'} or residue_at(window, 1) in {'D', 'E'},
        'motif_PRMT5_acidic_within2': any(residue_at(window, i) in {'D', 'E'} for i in [-2, -1, 1, 2]),
    }


def odds_ratio(a: int, b: int, c: int, d: int) -> float:
    return ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))


def motif_enrichment(target: pd.DataFrame, background: pd.DataFrame, motif_columns: list[str], background_label: str) -> pd.DataFrame:
    rows = []
    for column in motif_columns:
        a = int(target[column].sum())
        b = int(len(target) - a)
        c = int(background[column].sum())
        d = int(len(background) - c)
        rows.append({
            'motif_family': column.replace('motif_', ''),
            'target_count': a,
            'target_fraction': a / len(target) if len(target) else np.nan,
            'background_count': c,
            'background_fraction': c / len(background) if len(background) else np.nan,
            'odds_ratio': odds_ratio(a, b, c, d),
            'log2_odds_ratio': log2_odds_ratio(odds_ratio(a, b, c, d)),
            'p_value': float(fisher_exact([[a, b], [c, d]], alternative='two-sided').pvalue),
            'background_model': background_label,
        })
    out = pd.DataFrame(rows)
    return out.sort_values(['odds_ratio', 'target_count'], ascending=[False, False])


def positional_enrichment(target_windows: pd.Series, background_windows: pd.Series) -> pd.DataFrame:
    rows = []
    if target_windows.empty or background_windows.empty:
        return pd.DataFrame()
    flank = len(str(target_windows.iloc[0])) // 2
    amino_acids = sorted(set(''.join(target_windows.astype(str).tolist() + background_windows.astype(str).tolist())) - {'_'})
    for rel in range(-flank, flank + 1):
        if rel == 0:
            continue
        target_chars = target_windows.astype(str).map(lambda x: residue_at(x, rel))
        background_chars = background_windows.astype(str).map(lambda x: residue_at(x, rel))
        for aa in amino_acids:
            a = int((target_chars == aa).sum())
            b = int(len(target_chars) - a)
            c = int((background_chars == aa).sum())
            d = int(len(background_chars) - c)
            rows.append({
                'relative_position': rel,
                'amino_acid': aa,
                'target_count': a,
                'background_count': c,
                'odds_ratio': odds_ratio(a, b, c, d),
                'log2_odds_ratio': math.log2(odds_ratio(a, b, c, d)),
                'p_value': float(fisher_exact([[a, b], [c, d]], alternative='two-sided').pvalue),
            })
    out = pd.DataFrame(rows)
    return out


def top_protein_label(row: pd.Series) -> str:
    gene = str(row.get('substrate_genename', '') or '').strip()
    accession = str(row.get('canonical_UniProtAC', '') or '').strip()
    if gene:
        return f'{gene} ({accession})'
    return accession


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--integrated-sites', required=True)
    ap.add_argument('--canonical-fasta', required=True)
    ap.add_argument('--outdir', default='results/motif_arg_odds')
    ap.add_argument('--window-flank', type=int, default=7)
    ap.add_argument('--window-flank-wide', type=int, default=15)
    ap.add_argument('--min-arg-count', type=int, default=10)
    ap.add_argument('--min-site-count', type=int, default=2)
    args = ap.parse_args()

    sites = pd.read_csv(args.integrated_sites, sep='\t', low_memory=False)
    sites = canonical_site_table(sites, accession_col='canonical_UniProtAC', position_col='corrected_position')
    sequences = parse_fasta(args.canonical_fasta)

    sites = sites[sites['canonical_UniProtAC'].isin(sequences)].copy()
    sites['window_15'] = sites.apply(lambda row: sequence_window(sequences.get(row['canonical_UniProtAC'], ''), int(row['corrected_position']), args.window_flank), axis=1)
    sites['window_31'] = sites.apply(lambda row: sequence_window(sequences.get(row['canonical_UniProtAC'], ''), int(row['corrected_position']), args.window_flank_wide), axis=1)
    motif_columns = sorted(family_flags(sites['window_15'].iloc[0]).keys()) if not sites.empty else []
    motif_flags_df = pd.DataFrame([family_flags(window) for window in sites['window_15']])
    sites = pd.concat([sites.reset_index(drop=True), motif_flags_df], axis=1)

    methylated_proteins = sorted(set(sites['canonical_UniProtAC'].dropna().astype(str)))
    all_arg_proteome = residue_background_from_sequences(sequences, residue='R', flank=args.window_flank)
    all_arg_proteome = all_arg_proteome.rename(columns={f'window_{2 * args.window_flank + 1}': 'window_15'})
    all_arg_methyl_proteins = all_arg_proteome[all_arg_proteome['canonical_UniProtAC'].isin(methylated_proteins)].copy()
    methyl_site_keys = set(sites['site_key'].astype(str))
    background_nonmethyl = all_arg_methyl_proteins[~all_arg_methyl_proteins['site_key'].isin(methyl_site_keys)].copy()
    if not background_nonmethyl.empty:
        bg_flags = pd.DataFrame([family_flags(window) for window in background_nonmethyl['window_15']])
        background_nonmethyl = pd.concat([background_nonmethyl.reset_index(drop=True), bg_flags], axis=1)
    if not all_arg_proteome.empty:
        proteome_flags = pd.DataFrame([family_flags(window) for window in all_arg_proteome['window_15']])
        all_arg_proteome = pd.concat([all_arg_proteome.reset_index(drop=True), proteome_flags], axis=1)

    motif_vs_methyl_proteins = motif_enrichment(sites, background_nonmethyl, motif_columns, 'arginines_in_methylated_proteins')
    motif_vs_proteome = motif_enrichment(sites, all_arg_proteome, motif_columns, 'all_proteome_arginines')
    positional = positional_enrichment(sites['window_15'], background_nonmethyl['window_15'])

    arg_counts = residue_background_from_sequences(sequences, residue='R')[['canonical_UniProtAC', 'site_key']].groupby('canonical_UniProtAC').size().rename('arginine_count').reset_index()
    protein_counts = sites.groupby(['canonical_UniProtAC', 'substrate_genename'], as_index=False).size().rename(columns={'size': 'methyl_site_count'})
    protein_counts = protein_counts.merge(arg_counts, how='left', on='canonical_UniProtAC')
    total_methyl = int(protein_counts['methyl_site_count'].sum())
    total_arg = int(arg_counts['arginine_count'].sum())
    protein_counts['other_methyl'] = total_methyl - protein_counts['methyl_site_count']
    protein_counts['other_arg'] = total_arg - protein_counts['arginine_count']
    protein_counts['odds_ratio'] = protein_counts.apply(
        lambda row: odds_ratio(
            int(row['methyl_site_count']),
            max(int(row['arginine_count']) - int(row['methyl_site_count']), 0),
            int(row['other_methyl']),
            max(int(row['other_arg']) - int(row['other_methyl']), 0),
        ),
        axis=1,
    )
    protein_counts['log2_odds_ratio'] = protein_counts['odds_ratio'].map(log2_odds_ratio)
    protein_counts['p_value'] = protein_counts.apply(
        lambda row: float(
            fisher_exact(
                [
                    [int(row['methyl_site_count']), max(int(row['arginine_count']) - int(row['methyl_site_count']), 0)],
                    [int(row['other_methyl']), max(int(row['other_arg']) - int(row['other_methyl']), 0)],
                ],
                alternative='two-sided',
            ).pvalue
        ),
        axis=1,
    )
    protein_counts = protein_counts[(protein_counts['arginine_count'] >= args.min_arg_count) & (protein_counts['methyl_site_count'] >= args.min_site_count)].copy()
    protein_counts['label'] = protein_counts.apply(top_protein_label, axis=1)
    protein_counts = protein_counts.sort_values(['odds_ratio', 'methyl_site_count'], ascending=[False, False])

    context_rows = []
    for column in motif_columns:
        a = int(sites[column].sum())
        b = int(len(sites) - a)
        c = int(background_nonmethyl[column].sum())
        d = int(len(background_nonmethyl) - c)
        context_rows.append({
            'context': column.replace('motif_', ''),
            'site_count': a,
            'odds_ratio': odds_ratio(a, b, c, d),
            'log2_odds_ratio': log2_odds_ratio(odds_ratio(a, b, c, d)),
            'p_value': float(fisher_exact([[a, b], [c, d]], alternative='two-sided').pvalue),
            'background_model': 'arginines_in_methylated_proteins',
        })
    per_proteome_context = pd.DataFrame(context_rows).sort_values('odds_ratio', ascending=False)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(sites, outdir / 'sequence_windows.tsv')
    save_table(motif_vs_methyl_proteins, outdir / 'motif_family_enrichment_vs_methylated_protein_arginines.tsv')
    save_table(motif_vs_proteome, outdir / 'motif_family_enrichment_vs_all_proteome_arginines.tsv')
    save_table(positional, outdir / 'positional_amino_acid_enrichment.tsv')
    save_table(protein_counts, outdir / 'per_protein_arg_odds.tsv')
    save_table(per_proteome_context, outdir / 'per_proteome_arg_context_odds.tsv')

    plot_motif = motif_vs_methyl_proteins.head(12).sort_values('odds_ratio')
    fig, ax = plt.subplots(figsize=(9.0, 6.0))
    ax.barh(plot_motif['motif_family'], plot_motif['log2_odds_ratio'])
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlabel('log2(OR) vs arginines in methylated proteins')
    ax.set_ylabel('Motif family')
    ax.set_title('Methylarginine motif-family enrichment')
    for y, v, n, p in zip(plot_motif['motif_family'], plot_motif['log2_odds_ratio'], plot_motif['target_count'], plot_motif['p_value']):
        ax.text(v, y, f'  n={n}, p={format_p_value(p)}', va='center', ha='left' if v >= 0 else 'right', fontsize=8.5)
    fig.tight_layout()
    save_figure(fig, outdir / 'arg_methyl_motif_family_enrichment')

    heatmap_df = positional.copy()
    selected_aa = (
        heatmap_df.groupby('amino_acid')['target_count'].sum()
        .sort_values(ascending=False)
        .head(12)
        .index.tolist()
    )
    if selected_aa:
        matrix = heatmap_df[heatmap_df['amino_acid'].isin(selected_aa)].pivot(index='amino_acid', columns='relative_position', values='log2_odds_ratio').fillna(0.0)
        fig2, ax2 = plt.subplots(figsize=(10.5, 5.5))
        im = ax2.imshow(matrix.to_numpy(), aspect='auto', cmap='coolwarm', vmin=-3, vmax=3)
        ax2.set_yticks(range(len(matrix.index)))
        ax2.set_yticklabels(matrix.index)
        ax2.set_xticks(range(len(matrix.columns)))
        ax2.set_xticklabels(matrix.columns)
        ax2.set_xlabel('Position relative to methyl-Arg')
        ax2.set_title('Positional amino-acid enrichment around methylarginines')
        fig2.colorbar(im, ax=ax2, label='log2(OR)')
        fig2.tight_layout()
        save_figure(fig2, outdir / 'arg_methyl_positional_enrichment_heatmap')

    top_or = protein_counts.head(20).sort_values('odds_ratio')
    fig3, ax3 = plt.subplots(figsize=(10.0, 7.0))
    ax3.barh(top_or['label'], top_or['log2_odds_ratio'])
    ax3.axvline(0, color='black', linewidth=0.8)
    ax3.set_xlabel('log2 Arg-normalized odds ratio')
    ax3.set_ylabel('Protein')
    ax3.set_title('Top proteins by methylarginine Arg odds ratio')
    for y, v, n, p in zip(top_or['label'], top_or['log2_odds_ratio'], top_or['methyl_site_count'], top_or['p_value']):
        ax3.text(v, y, f'  n={n}, p={format_p_value(p)}', va='center', ha='left' if v >= 0 else 'right', fontsize=8)
    fig3.tight_layout()
    save_figure(fig3, outdir / 'top_proteins_by_arg_odds_ratio')

    fig4, ax4 = plt.subplots(figsize=(7.0, 5.2))
    ax4.scatter(protein_counts['methyl_site_count'], np.log2(protein_counts['odds_ratio']), alpha=0.55)
    ax4.set_xlabel('Methylarginine site count')
    ax4.set_ylabel('log2 Arg odds ratio')
    ax4.set_title('Protein-level methylarginine burden and Arg enrichment')
    fig4.tight_layout()
    save_figure(fig4, outdir / 'protein_site_count_vs_arg_odds_ratio')

    print({
        'methyl_sites': len(sites),
        'background_arginines_same_proteins': len(background_nonmethyl),
        'proteins_with_arg_odds': len(protein_counts),
        'output': str(outdir),
    })


if __name__ == '__main__':
    main()

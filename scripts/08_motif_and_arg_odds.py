from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest, fisher_exact

from common import (
    BREWER_COLORS,
    add_q_values,
    annotate_barh,
    apply_paper_style,
    balanced_category_subset,
    canonical_site_table,
    format_p_value,
    log2_odds_ratio,
    parse_fasta,
    residue_background_from_sequences,
    save_figure,
    save_table,
    sequence_window,
    set_symmetric_xlim,
    signed_bar_colors,
    style_axis,
)


HYDROPHOBIC = set('AILMFWVY')
HEATMAP_ORDER = ['G', 'A', 'P', 'D', 'E', 'R', 'Y', 'F', 'W', 'L', 'I', 'V', 'M', 'C', 'K', 'Q', 'N', 'S', 'T', 'H']
FOCUS_HEATMAP_ORDER = ['G', 'A', 'P', 'D', 'E', 'R']


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
        # A centered methyl-Arg can participate in RXR in either orientation:
        # R-X-[methyl-R] or [methyl-R]-X-R.
        'motif_RXR': residue_at(window, -2) == 'R' or residue_at(window, 2) == 'R',
        'motif_CARM1_proline_rich': 'P' in neighborhood,
        'motif_CARM1_hydrophobic_proline': ('P' in neighborhood) and any(aa in HYDROPHOBIC for aa in near3),
        'motif_DR': residue_at(window, -1) == 'D',
        'motif_RD': residue_at(window, 1) == 'D',
        'motif_ER': residue_at(window, -1) == 'E',
        'motif_RE': residue_at(window, 1) == 'E',
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
    out = add_q_values(pd.DataFrame(rows))
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
    return add_q_values(out)


def top_protein_label(row: pd.Series) -> str:
    gene = str(row.get('substrate_genename', '') or '').strip()
    accession = str(row.get('canonical_UniProtAC', '') or '').strip()
    if gene:
        return f'{gene} ({accession})'
    return accession


def selected_heatmap_rows(heatmap_df: pd.DataFrame, limit: int = 20) -> list[str]:
    observed = set(heatmap_df['amino_acid'].astype(str))
    ordered = [aa for aa in HEATMAP_ORDER if aa in observed]
    if len(ordered) >= limit:
        return ordered[:limit]
    remaining = (
        heatmap_df[~heatmap_df['amino_acid'].isin(ordered)]
        .assign(abs_effect=lambda x: x['log2_odds_ratio'].abs())
        .groupby('amino_acid')['abs_effect']
        .max()
        .sort_values(ascending=False)
        .index.tolist()
    )
    return (ordered + remaining)[:limit]


def build_heatmap_matrices(heatmap_df: pd.DataFrame, selected_aa: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    matrix = heatmap_df[heatmap_df['amino_acid'].isin(selected_aa)].pivot(index='amino_acid', columns='relative_position', values='log2_odds_ratio').fillna(0.0)
    matrix = matrix.reindex(index=[aa for aa in selected_aa if aa in matrix.index])
    p_matrix = heatmap_df[heatmap_df['amino_acid'].isin(selected_aa)].pivot(index='amino_acid', columns='relative_position', values='p_value')
    p_matrix = p_matrix.reindex(index=matrix.index, columns=matrix.columns)
    p_floor = np.nextafter(0.0, 1.0)
    neglog10_p = -np.log10(p_matrix.fillna(1.0).clip(lower=p_floor)).clip(upper=50)
    return matrix, p_matrix, neglog10_p


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--integrated-sites', required=True)
    ap.add_argument('--canonical-fasta', required=True)
    ap.add_argument('--outdir', default='results/motif_arg_odds')
    ap.add_argument('--window-flank', type=int, default=7)
    ap.add_argument('--window-flank-wide', type=int, default=15)
    ap.add_argument('--min-arg-count', type=int, default=10)
    ap.add_argument('--min-site-count', type=int, default=2)
    ap.add_argument('--shrinkage-prior-strength', type=float, default=200.0)
    args = ap.parse_args()

    sites = pd.read_csv(args.integrated_sites, sep='\t', low_memory=False)
    sites = canonical_site_table(sites, accession_col='canonical_UniProtAC', position_col='corrected_position')
    input_site_count = len(sites)
    sequences = parse_fasta(args.canonical_fasta)
    apply_paper_style()

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
    protein_counts = add_q_values(protein_counts)
    global_rate = total_methyl / total_arg if total_arg else np.nan
    alpha0 = global_rate * args.shrinkage_prior_strength if pd.notna(global_rate) else np.nan
    beta0 = (1 - global_rate) * args.shrinkage_prior_strength if pd.notna(global_rate) else np.nan
    protein_counts['expected_global_rate'] = global_rate
    protein_counts['posterior_rate'] = (
        protein_counts['methyl_site_count'] + alpha0
    ) / (
        protein_counts['arginine_count'] + alpha0 + beta0
    )
    protein_counts['shrinkage_log2_enrichment'] = protein_counts['posterior_rate'].map(
        lambda x: log2_odds_ratio(x / global_rate) if pd.notna(x) and global_rate > 0 else np.nan
    )
    protein_counts['binom_p_value'] = protein_counts.apply(
        lambda row: float(binomtest(int(row['methyl_site_count']), int(row['arginine_count']), global_rate, alternative='two-sided').pvalue)
        if pd.notna(global_rate) and int(row['arginine_count']) > 0 else np.nan,
        axis=1,
    )
    protein_counts = add_q_values(protein_counts, p_col='binom_p_value', out_col='binom_q_value')
    protein_counts['protein_length'] = protein_counts['canonical_UniProtAC'].map(lambda acc: len(sequences.get(acc, '')))
    protein_counts = protein_counts[(protein_counts['arginine_count'] >= args.min_arg_count) & (protein_counts['methyl_site_count'] >= args.min_site_count)].copy()
    protein_counts['label'] = protein_counts.apply(top_protein_label, axis=1)
    protein_counts = protein_counts.sort_values(['shrinkage_log2_enrichment', 'methyl_site_count'], ascending=[False, False])

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(sites, outdir / 'sequence_windows.tsv')
    save_table(motif_vs_methyl_proteins, outdir / 'motif_family_enrichment_vs_methylated_protein_arginines.tsv')
    save_table(motif_vs_proteome, outdir / 'motif_family_enrichment_vs_all_proteome_arginines.tsv')
    save_table(positional, outdir / 'positional_amino_acid_enrichment.tsv')
    save_table(protein_counts, outdir / 'per_protein_arg_odds.tsv')
    save_table(protein_counts, outdir / 'per_protein_shrinkage_prioritization.tsv')
    save_table(
        pd.DataFrame([{
            'integrated_input_sites': input_site_count,
            'analyzed_methyl_sites_with_sequence_window': len(sites),
            'excluded_sites_without_usable_sequence': input_site_count - len(sites),
            'background_arginines_same_proteins': len(background_nonmethyl),
            'all_proteome_arginines': len(all_arg_proteome),
            'proteins_with_arg_odds': len(protein_counts),
        }]),
        outdir / 'motif_analysis_qc.tsv',
    )

    plot_motif = motif_vs_methyl_proteins.sort_values('log2_odds_ratio')
    fig, ax = plt.subplots(figsize=(9.0, 6.0))
    ax.barh(
        plot_motif['motif_family'],
        plot_motif['log2_odds_ratio'],
        color=signed_bar_colors(plot_motif['log2_odds_ratio'], positive=BREWER_COLORS['green'], negative=BREWER_COLORS['purple']),
        edgecolor='white',
        linewidth=0.8,
    )
    style_axis(ax, zero='x')
    ax.set_xlabel('log2(OR) vs arginines in methylated proteins')
    ax.set_ylabel('Motif family')
    ax.set_title('Methylarginine motif-family enrichment')
    set_symmetric_xlim(ax, plot_motif['log2_odds_ratio'], annotation_pad_ratio=0.4, center_on_zero=False)
    annotate_barh(
        ax,
        plot_motif['motif_family'],
        plot_motif['log2_odds_ratio'],
        [f'n={n}, p={format_p_value(p)}, q={format_p_value(q)}' for n, p, q in zip(plot_motif['target_count'], plot_motif['p_value'], plot_motif['q_value'])],
        fontsize=8.0,
    )
    fig.tight_layout()
    save_figure(fig, outdir / 'arg_methyl_motif_family_enrichment')

    heatmap_df = positional.copy()
    selected_aa = selected_heatmap_rows(heatmap_df, limit=20)
    if selected_aa:
        matrix, p_matrix, neglog10_p = build_heatmap_matrices(heatmap_df, selected_aa)

        fig2, (ax2, ax2b) = plt.subplots(
            ncols=2,
            figsize=(15.0, 7.8),
            gridspec_kw={'width_ratios': [1.0, 1.0]},
        )
        im = ax2.imshow(matrix.to_numpy(), aspect='auto', cmap='RdBu_r', vmin=-3, vmax=3)
        ax2.set_yticks(range(len(matrix.index)))
        ax2.set_yticklabels(matrix.index)
        ax2.set_xticks(range(len(matrix.columns)))
        ax2.set_xticklabels(matrix.columns)
        ax2.set_xlabel('Position relative to methyl-Arg')
        ax2.set_title('Positional amino-acid enrichment around methylarginines')
        for row_idx, aa in enumerate(matrix.index):
            for col_idx, rel in enumerate(matrix.columns):
                pval = p_matrix.loc[aa, rel]
                if pd.isna(pval):
                    continue
                if pval < 1e-20:
                    mark = '***'
                elif pval < 1e-5:
                    mark = '**'
                elif pval < 0.05:
                    mark = '*'
                else:
                    mark = ''
                if mark:
                    ax2.text(col_idx, row_idx, mark, ha='center', va='center', fontsize=7, color='black')
        fig2.colorbar(im, ax=ax2, label='log2(OR)')

        im_p = ax2b.imshow(neglog10_p.to_numpy(), aspect='auto', cmap='YlGnBu', vmin=0, vmax=50)
        ax2b.set_yticks(range(len(neglog10_p.index)))
        ax2b.set_yticklabels(neglog10_p.index)
        ax2b.set_xticks(range(len(neglog10_p.columns)))
        ax2b.set_xticklabels(neglog10_p.columns)
        ax2b.set_xlabel('Position relative to methyl-Arg')
        ax2b.set_title('Positional Fisher exact significance')
        fig2.colorbar(im_p, ax=ax2b, label='-log10(p), capped at 50')
        fig2.tight_layout()
        save_figure(fig2, outdir / 'arg_methyl_positional_enrichment_heatmap')

        focus_aa = [aa for aa in FOCUS_HEATMAP_ORDER if aa in set(heatmap_df['amino_acid'].astype(str))]
        if focus_aa:
            focus_matrix, focus_p, _ = build_heatmap_matrices(heatmap_df, focus_aa)
            fig2_focus, ax2_focus = plt.subplots(figsize=(8.2, 3.8))
            im_focus = ax2_focus.imshow(focus_matrix.to_numpy(), aspect='auto', cmap='RdBu_r', vmin=-3, vmax=3)
            ax2_focus.set_yticks(range(len(focus_matrix.index)))
            ax2_focus.set_yticklabels(focus_matrix.index)
            ax2_focus.set_xticks(range(len(focus_matrix.columns)))
            ax2_focus.set_xticklabels(focus_matrix.columns)
            ax2_focus.set_xlabel('Position relative to methyl-Arg')
            ax2_focus.set_title('Focused positional enrichment around methylarginines')
            for row_idx, aa in enumerate(focus_matrix.index):
                for col_idx, rel in enumerate(focus_matrix.columns):
                    pval = focus_p.loc[aa, rel]
                    if pd.isna(pval):
                        continue
                    if pval < 1e-20:
                        mark = '***'
                    elif pval < 1e-5:
                        mark = '**'
                    elif pval < 0.05:
                        mark = '*'
                    else:
                        mark = ''
                    if mark:
                        ax2_focus.text(col_idx, row_idx, mark, ha='center', va='center', fontsize=7, color='black')
            fig2_focus.colorbar(im_focus, ax=ax2_focus, label='log2(OR)')
            fig2_focus.tight_layout()
            save_figure(fig2_focus, outdir / 'arg_methyl_positional_enrichment_focus')

        acidic_focus = positional[
            positional['amino_acid'].isin(['D', 'E']) & positional['relative_position'].isin([-2, -1, 1, 2])
        ].copy()
        acidic_focus['position_label'] = acidic_focus['relative_position'].map({-2: '-2', -1: '-1', 1: '+1', 2: '+2'})
        save_table(acidic_focus.sort_values(['amino_acid', 'relative_position']), outdir / 'acidic_context_focus.tsv')
        fig2b, ax2c = plt.subplots(figsize=(6.4, 4.6))
        width = 0.38
        positions = np.arange(4)
        d_rows = acidic_focus[acidic_focus['amino_acid'] == 'D'].sort_values('relative_position')
        e_rows = acidic_focus[acidic_focus['amino_acid'] == 'E'].sort_values('relative_position')
        ax2c.bar(positions - width / 2, d_rows['log2_odds_ratio'], width=width, color=BREWER_COLORS['orange'], label='Asp (D)', edgecolor='white', linewidth=0.8)
        ax2c.bar(positions + width / 2, e_rows['log2_odds_ratio'], width=width, color=BREWER_COLORS['blue'], label='Glu (E)', edgecolor='white', linewidth=0.8)
        style_axis(ax2c, zero='y', grid_axis='y')
        ax2c.set_xticks(positions)
        ax2c.set_xticklabels(['-2', '-1', '+1', '+2'])
        ax2c.set_xlabel('Position relative to methyl-Arg')
        ax2c.set_ylabel('log2(OR) vs same-protein non-methyl Arg')
        ax2c.set_title('Localized acidic asymmetry near methylarginines')
        ax2c.legend(frameon=False, loc='upper right')
        for shift, rows in [(-width / 2, d_rows), (width / 2, e_rows)]:
            for idx, row in enumerate(rows.itertuples(index=False)):
                ax2c.text(
                    positions[idx] + shift,
                    float(row.log2_odds_ratio),
                    f"p={format_p_value(row.p_value)}",
                    ha='center',
                    va='bottom' if float(row.log2_odds_ratio) >= 0 else 'top',
                    fontsize=7.2,
                    color=BREWER_COLORS['dark_gray'],
                )
        fig2b.tight_layout()
        save_figure(fig2b, outdir / 'arg_methyl_acidic_asymmetry')

    top_shrunk = protein_counts.head(20).sort_values('shrinkage_log2_enrichment')
    fig3, ax3 = plt.subplots(figsize=(10.0, 7.0))
    ax3.barh(
        top_shrunk['label'],
        top_shrunk['shrinkage_log2_enrichment'],
        color=signed_bar_colors(top_shrunk['shrinkage_log2_enrichment'], positive=BREWER_COLORS['green'], negative=BREWER_COLORS['purple']),
        edgecolor='white',
        linewidth=0.8,
    )
    style_axis(ax3, zero='x')
    ax3.set_xlabel('log2 shrinkage-adjusted enrichment')
    ax3.set_ylabel('Protein')
    ax3.set_title('Top proteins by shrinkage-adjusted methylarginine enrichment')
    set_symmetric_xlim(ax3, top_shrunk['shrinkage_log2_enrichment'], annotation_pad_ratio=0.42, center_on_zero=False)
    annotate_barh(
        ax3,
        top_shrunk['label'],
        top_shrunk['shrinkage_log2_enrichment'],
        [f'sites={n}, q={format_p_value(q)}' for n, q in zip(top_shrunk['methyl_site_count'], top_shrunk['binom_q_value'])],
        fontsize=7.6,
    )
    fig3.tight_layout()
    save_figure(fig3, outdir / 'top_proteins_by_shrinkage_score')

    fig4, ax4 = plt.subplots(figsize=(7.0, 5.2))
    ax4.scatter(
        protein_counts['methyl_site_count'],
        protein_counts['shrinkage_log2_enrichment'],
        alpha=0.6,
        s=46,
        color=BREWER_COLORS['blue'],
        edgecolor='white',
        linewidth=0.4,
    )
    style_axis(ax4, grid_axis='both')
    ax4.set_xlabel('Methylarginine site count')
    ax4.set_ylabel('log2 shrinkage-adjusted enrichment')
    ax4.set_title('Protein-level methylarginine burden and shrinkage score')
    label_points = protein_counts.nlargest(10, 'shrinkage_log2_enrichment')
    for row in label_points.itertuples(index=False):
        ax4.text(
            float(row.methyl_site_count) + 0.45,
            float(row.shrinkage_log2_enrichment) + 0.03,
            str(row.label).split(' (')[0],
            fontsize=7.3,
            color=BREWER_COLORS['dark_gray'],
        )
    fig4.tight_layout()
    save_figure(fig4, outdir / 'protein_site_count_vs_shrinkage_score')

    top_or = protein_counts.head(20).sort_values('odds_ratio')
    fig5, ax5 = plt.subplots(figsize=(10.0, 7.0))
    ax5.barh(
        top_or['label'],
        top_or['log2_odds_ratio'],
        color=signed_bar_colors(top_or['log2_odds_ratio'], positive=BREWER_COLORS['green'], negative=BREWER_COLORS['purple']),
        edgecolor='white',
        linewidth=0.8,
    )
    style_axis(ax5, zero='x')
    ax5.set_xlabel('log2 Arg-normalized odds ratio')
    ax5.set_ylabel('Protein')
    ax5.set_title('Top proteins by raw Arg odds ratio')
    set_symmetric_xlim(ax5, top_or['log2_odds_ratio'], annotation_pad_ratio=0.42, center_on_zero=False)
    annotate_barh(
        ax5,
        top_or['label'],
        top_or['log2_odds_ratio'],
        [f'n={n}, p={format_p_value(p)}, q={format_p_value(q)}' for n, p, q in zip(top_or['methyl_site_count'], top_or['p_value'], top_or['q_value'])],
        fontsize=7.6,
    )
    fig5.tight_layout()
    save_figure(fig5, outdir / 'top_proteins_by_arg_odds_ratio')

    fig6, ax6 = plt.subplots(figsize=(7.0, 5.2))
    ax6.scatter(
        protein_counts['methyl_site_count'],
        np.log2(protein_counts['odds_ratio']),
        alpha=0.6,
        s=46,
        color=BREWER_COLORS['purple'],
        edgecolor='white',
        linewidth=0.4,
    )
    style_axis(ax6, grid_axis='both')
    ax6.set_xlabel('Methylarginine site count')
    ax6.set_ylabel('log2 Arg odds ratio')
    ax6.set_title('Protein-level methylarginine burden and Arg enrichment')
    fig6.tight_layout()
    save_figure(fig6, outdir / 'protein_site_count_vs_arg_odds_ratio')

    print({
        'methyl_sites': len(sites),
        'background_arginines_same_proteins': len(background_nonmethyl),
        'proteins_with_arg_odds': len(protein_counts),
        'output': str(outdir),
    })


if __name__ == '__main__':
    main()

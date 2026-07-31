from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from common import (
    BREWER_COLORS,
    add_q_values,
    annotate_barh,
    apply_paper_style,
    canonical_site_table,
    count_over_total_labels,
    fisher_like_enrichment,
    log2_odds_ratio,
    parse_fasta,
    q_threshold_note,
    residue_background_from_sequences,
    save_figure,
    save_table,
    set_symmetric_xlim,
    signed_bar_colors,
    style_axis,
)


def normalize_disorder_intervals(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    rename_map = {}
    for source, target in [
        ('accession', 'canonical_UniProtAC'),
        ('Entry', 'canonical_UniProtAC'),
        ('start', 'fragment_start'),
        ('end', 'fragment_end'),
    ]:
        if source in out.columns and target not in out.columns:
            rename_map[source] = target
    out = out.rename(columns=rename_map)
    required = ['canonical_UniProtAC', 'fragment_start', 'fragment_end']
    missing = [col for col in required if col not in out.columns]
    if missing:
        raise ValueError(f'missing disorder interval columns: {missing}')
    out['canonical_UniProtAC'] = out['canonical_UniProtAC'].astype(str).str.strip()
    out['fragment_start'] = pd.to_numeric(out['fragment_start'], errors='coerce')
    out['fragment_end'] = pd.to_numeric(out['fragment_end'], errors='coerce')
    out = out[out['fragment_start'].notna() & out['fragment_end'].notna()].copy()
    out['fragment_start'] = out['fragment_start'].astype(int)
    out['fragment_end'] = out['fragment_end'].astype(int)
    return out.sort_values(['canonical_UniProtAC', 'fragment_start', 'fragment_end']).reset_index(drop=True)


def classify_positions(intervals: pd.DataFrame, positions: np.ndarray, boundary_window: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = len(positions)
    classes = np.full(n, 'no_disorder_annotation', dtype=object)
    distances = np.full(n, np.nan, dtype=float)
    has_annotation = np.zeros(n, dtype=bool)
    if intervals.empty or n == 0:
        return classes, distances, has_annotation

    starts = intervals['fragment_start'].to_numpy(dtype=int)
    ends = intervals['fragment_end'].to_numpy(dtype=int)
    min_distance = np.full(n, np.iinfo(np.int32).max, dtype=int)
    in_interval = np.zeros(n, dtype=bool)

    for start, end in zip(starts, ends):
        overlap = (positions >= start) & (positions <= end)
        in_interval |= overlap
        dist = np.where(
            positions < start,
            start - positions,
            np.where(positions > end, positions - end, np.minimum(positions - start, end - positions)),
        )
        min_distance = np.minimum(min_distance, dist)

    classes[in_interval] = 'disordered'
    distances[in_interval] = min_distance[in_interval]
    boundary = (~in_interval) & (min_distance <= boundary_window)
    classes[boundary] = 'disorder_boundary'
    distances[boundary] = min_distance[boundary]
    ordered = (~in_interval) & (min_distance > boundary_window)
    classes[ordered] = 'ordered'
    distances[ordered] = min_distance[ordered]
    # If any intervals exist for this protein, all queried positions are classifiable
    # as disordered, boundary-adjacent, or ordered relative to that annotation layer.
    has_annotation[:] = True
    return classes, distances, has_annotation


def idr_edge_bin(distance) -> str | pd.NA:
    if pd.isna(distance):
        return pd.NA
    distance = float(distance)
    if distance <= 5:
        return 'IDR edge <=5 aa'
    if distance <= 10:
        return 'IDR edge 6-10 aa'
    if distance <= 20:
        return 'IDR edge 11-20 aa'
    if distance <= 40:
        return 'IDR edge 21-40 aa'
    return 'Distal IDR >40 aa'


def cumulative_edge_curve(target_dist: pd.Series, background_dist: pd.Series, max_distance: int = 40) -> pd.DataFrame:
    target = target_dist.dropna().astype(float)
    background = background_dist.dropna().astype(float)
    rows = []
    for threshold in range(1, max_distance + 1):
        a = int((target <= threshold).sum())
        b = int(len(target) - a)
        c = int((background <= threshold).sum())
        d = int(len(background) - c)
        odds = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)) if (b + 0.5) * (c + 0.5) > 0 else np.nan
        rows.append({
            'distance_aa': threshold,
            'target_within_distance': a,
            'target_total': len(target),
            'background_within_distance': c,
            'background_total': len(background),
            'target_fraction_within_distance': a / len(target) if len(target) else np.nan,
            'background_fraction_within_distance': c / len(background) if len(background) else np.nan,
            'odds_ratio': odds,
            'log2_odds_ratio': log2_odds_ratio(odds) if pd.notna(odds) else np.nan,
            'p_value': float(fisher_exact([[a, b], [c, d]], alternative='two-sided').pvalue),
        })
    return add_q_values(pd.DataFrame(rows))


def rolling_edge_curve(
    target_dist: pd.Series,
    background_dist: pd.Series,
    window_size: int = 5,
    max_distance: int = 40,
) -> pd.DataFrame:
    target = target_dist.dropna().astype(float)
    background = background_dist.dropna().astype(float)
    rows = []
    for start in range(1, max_distance - window_size + 2):
        end = start + window_size - 1
        a = int(((target >= start) & (target <= end)).sum())
        b = int(len(target) - a)
        c = int(((background >= start) & (background <= end)).sum())
        d = int(len(background) - c)
        odds = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)) if (b + 0.5) * (c + 0.5) > 0 else np.nan
        rows.append({
            'window_start_aa': start,
            'window_end_aa': end,
            'window_center_aa': (start + end) / 2.0,
            'window_label': f'{start}-{end} aa',
            'target_in_window': a,
            'target_total': len(target),
            'background_in_window': c,
            'background_total': len(background),
            'target_fraction_in_window': a / len(target) if len(target) else np.nan,
            'background_fraction_in_window': c / len(background) if len(background) else np.nan,
            'odds_ratio': odds,
            'log2_odds_ratio': log2_odds_ratio(odds) if pd.notna(odds) else np.nan,
            'p_value': float(fisher_exact([[a, b], [c, d]], alternative='two-sided').pvalue),
        })
    return add_q_values(pd.DataFrame(rows))


def disorder_label(label: str) -> str:
    return {
        'disordered': 'Disordered',
        'disorder_boundary': 'Outside-IDR edge (<=20 aa)',
        'ordered': 'Ordered',
        'no_disorder_annotation': 'No disorder annotation',
    }.get(str(label), str(label))


def annotate_frame(df: pd.DataFrame, intervals_by_acc: dict[str, pd.DataFrame], accession_col: str, position_col: str, boundary_window: int) -> pd.DataFrame:
    out = df.copy()
    out[accession_col] = out[accession_col].astype(str).str.strip()
    out[position_col] = pd.to_numeric(out[position_col], errors='coerce')
    out['disorder_context_class'] = 'no_position'
    out['nearest_disorder_edge_distance'] = pd.NA
    out['has_disorder_annotation'] = False

    valid_mask = (~out[accession_col].isin(['', 'nan', 'None'])) & out[position_col].notna()
    if valid_mask.any():
        valid = out.loc[valid_mask, [accession_col, position_col]].copy()
        valid[position_col] = valid[position_col].astype(int)
        for accession, group in valid.groupby(accession_col, sort=False):
            positions = group[position_col].to_numpy(dtype=int)
            intervals = intervals_by_acc.get(accession, pd.DataFrame())
            classes, distances, has_annotation = classify_positions(intervals, positions, boundary_window)
            out.loc[group.index, 'disorder_context_class'] = classes
            out.loc[group.index, 'nearest_disorder_edge_distance'] = distances
            out.loc[group.index, 'has_disorder_annotation'] = has_annotation & (classes != 'no_disorder_annotation')

    out['disorder_context_boundary_window'] = boundary_window
    out['disorder_context_layer'] = 'disorder_context'
    out['idr_binary_class'] = pd.NA
    out.loc[out['disorder_context_class'] == 'disordered', 'idr_binary_class'] = 'IDR'
    out.loc[out['disorder_context_class'].isin(['ordered', 'disorder_boundary']), 'idr_binary_class'] = 'non_IDR'
    out['idr_proximal_binary_class'] = pd.NA
    out.loc[out['disorder_context_class'].isin(['disordered', 'disorder_boundary']), 'idr_proximal_binary_class'] = 'IDR_proximal'
    out.loc[out['disorder_context_class'] == 'ordered', 'idr_proximal_binary_class'] = 'ordered_distal'
    out['idr_edge_bin'] = pd.NA
    disordered_mask = out['disorder_context_class'] == 'disordered'
    out.loc[disordered_mask, 'idr_edge_bin'] = out.loc[disordered_mask, 'nearest_disorder_edge_distance'].map(idr_edge_bin)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--sites', required=True)
    ap.add_argument('--disorder-intervals', required=True)
    ap.add_argument('--canonical-fasta', default='')
    ap.add_argument('--outdir', default='results/disorder_context')
    ap.add_argument('--accession-col', default='canonical_UniProtAC')
    ap.add_argument('--position-col', default='corrected_position')
    ap.add_argument('--boundary-window', type=int, default=20)
    args = ap.parse_args()

    sites = pd.read_csv(args.sites, sep='\t', low_memory=False)
    sites = canonical_site_table(sites, accession_col=args.accession_col, position_col=args.position_col)
    intervals = normalize_disorder_intervals(pd.read_csv(args.disorder_intervals, sep='\t', low_memory=False))
    apply_paper_style()
    intervals_by_acc = {
        accession: group.reset_index(drop=True)
        for accession, group in intervals.groupby('canonical_UniProtAC')
    }

    annotated_sites = annotate_frame(
        sites,
        intervals_by_acc=intervals_by_acc,
        accession_col='canonical_UniProtAC',
        position_col='corrected_position',
        boundary_window=args.boundary_window,
    )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(annotated_sites, outdir / 'sites_with_disorder_context.tsv')

    summary = annotated_sites['disorder_context_class'].value_counts(dropna=False).rename_axis('disorder_context_class').reset_index(name='site_count')
    summary['boundary_window'] = args.boundary_window
    save_table(summary, outdir / 'disorder_context_summary.tsv')

    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    plot_df = summary[summary['disorder_context_class'] != 'no_position'].copy()
    bar_colors = [
        {
            'disordered': BREWER_COLORS['dark_gray'],
            'disorder_boundary': BREWER_COLORS['mid_gray'],
            'ordered': BREWER_COLORS['light_gray'],
            'no_disorder_annotation': BREWER_COLORS['mid_gray'],
        }.get(cat, BREWER_COLORS['mid_gray'])
        for cat in plot_df['disorder_context_class']
    ]
    display_labels = plot_df['disorder_context_class'].map(disorder_label)
    ax.bar(display_labels, plot_df['site_count'], color=bar_colors, edgecolor='white', linewidth=0.8)
    style_axis(ax)
    ax.set_ylabel('Methylarginine sites')
    ax.set_title('Methylarginine disorder-context distribution')
    ax.tick_params(axis='x', rotation=20)
    fig.tight_layout()
    save_figure(fig, outdir / 'arg_methyl_disorder_context_distribution')

    if args.canonical_fasta:
        sequences = parse_fasta(args.canonical_fasta)
        proteome_arg = residue_background_from_sequences(sequences, residue='R')
        proteome_arg = proteome_arg.rename(columns={'position': 'corrected_position'})
        proteome_arg = annotate_frame(
            proteome_arg,
            intervals_by_acc=intervals_by_acc,
            accession_col='canonical_UniProtAC',
            position_col='corrected_position',
            boundary_window=args.boundary_window,
        )
        methyl_site_keys = set(annotated_sites['site_key'].astype(str))
        proteome_arg = proteome_arg[~proteome_arg['site_key'].isin(methyl_site_keys)].copy()
        methylated_proteins = sorted(set(annotated_sites['canonical_UniProtAC'].dropna().astype(str)))
        methyl_protein_arg = proteome_arg[proteome_arg['canonical_UniProtAC'].isin(methylated_proteins)].copy()
        enrich_all = fisher_like_enrichment(annotated_sites['disorder_context_class'], proteome_arg['disorder_context_class'])
        enrich_same = fisher_like_enrichment(annotated_sites['disorder_context_class'], methyl_protein_arg['disorder_context_class'])
        binary_all = fisher_like_enrichment(
            annotated_sites[annotated_sites['idr_binary_class'].notna()]['idr_binary_class'],
            proteome_arg[proteome_arg['idr_binary_class'].notna()]['idr_binary_class'],
        )
        binary_same = fisher_like_enrichment(
            annotated_sites[annotated_sites['idr_binary_class'].notna()]['idr_binary_class'],
            methyl_protein_arg[methyl_protein_arg['idr_binary_class'].notna()]['idr_binary_class'],
        )
        proximal_binary_all = fisher_like_enrichment(
            annotated_sites[annotated_sites['idr_proximal_binary_class'].notna()]['idr_proximal_binary_class'],
            proteome_arg[proteome_arg['idr_proximal_binary_class'].notna()]['idr_proximal_binary_class'],
        )
        proximal_binary_same = fisher_like_enrichment(
            annotated_sites[annotated_sites['idr_proximal_binary_class'].notna()]['idr_proximal_binary_class'],
            methyl_protein_arg[methyl_protein_arg['idr_proximal_binary_class'].notna()]['idr_proximal_binary_class'],
        )
        save_table(enrich_all, outdir / 'idr_enrichment_vs_all_arginines.tsv')
        save_table(enrich_same, outdir / 'idr_enrichment_vs_arginines_in_methylated_proteins.tsv')
        save_table(binary_all, outdir / 'idr_binary_enrichment_vs_all_arginines.tsv')
        save_table(binary_same, outdir / 'idr_binary_enrichment_vs_arginines_in_methylated_proteins.tsv')
        save_table(proximal_binary_all, outdir / 'idr_proximal_binary_enrichment_vs_all_arginines.tsv')
        save_table(proximal_binary_same, outdir / 'idr_proximal_binary_enrichment_vs_arginines_in_methylated_proteins.tsv')
        summary_rows = []
        for label, series in [
            ('strict_IDR', annotated_sites['idr_binary_class']),
            ('IDR_proximal', annotated_sites['idr_proximal_binary_class']),
        ]:
            valid = series.dropna()
            positive = 'IDR' if label == 'strict_IDR' else 'IDR_proximal'
            summary_rows.append({
                'metric': label,
                'classified_site_count': int(len(valid)),
                'positive_site_count': int((valid == positive).sum()),
                'positive_fraction': float((valid == positive).mean()) if len(valid) else np.nan,
            })
        save_table(pd.DataFrame(summary_rows), outdir / 'idr_binary_summary.tsv')

        edge_bin_all = fisher_like_enrichment(
            annotated_sites[annotated_sites['idr_edge_bin'].notna()]['idr_edge_bin'],
            proteome_arg[proteome_arg['idr_edge_bin'].notna()]['idr_edge_bin'],
        )
        edge_bin_same = fisher_like_enrichment(
            annotated_sites[annotated_sites['idr_edge_bin'].notna()]['idr_edge_bin'],
            methyl_protein_arg[methyl_protein_arg['idr_edge_bin'].notna()]['idr_edge_bin'],
        )
        edge_curve_same = cumulative_edge_curve(
            annotated_sites.loc[annotated_sites['idr_edge_bin'].notna(), 'nearest_disorder_edge_distance'],
            methyl_protein_arg.loc[methyl_protein_arg['idr_edge_bin'].notna(), 'nearest_disorder_edge_distance'],
            max_distance=40,
        )
        edge_roll_same = rolling_edge_curve(
            annotated_sites.loc[annotated_sites['idr_edge_bin'].notna(), 'nearest_disorder_edge_distance'],
            methyl_protein_arg.loc[methyl_protein_arg['idr_edge_bin'].notna(), 'nearest_disorder_edge_distance'],
            window_size=5,
            max_distance=40,
        )
        save_table(edge_bin_all, outdir / 'idr_edge_bin_enrichment_vs_all_disordered_arginines.tsv')
        save_table(edge_bin_same, outdir / 'idr_edge_bin_enrichment_vs_disordered_arginines_in_methylated_proteins.tsv')
        save_table(edge_curve_same, outdir / 'idr_edge_cumulative_enrichment.tsv')
        save_table(edge_roll_same, outdir / 'idr_edge_rolling_enrichment.tsv')

        compare = enrich_all[enrich_all['category'].isin(['disordered', 'disorder_boundary', 'ordered'])].sort_values('odds_ratio').copy()
        compare['display_category'] = compare['category'].map(disorder_label)
        fig2, ax2 = plt.subplots(figsize=(7.6, 4.9))
        colors = [
            {
                'disordered': BREWER_COLORS['dark_gray'],
                'disorder_boundary': BREWER_COLORS['mid_gray'],
                'ordered': BREWER_COLORS['light_gray'],
            }.get(cat, BREWER_COLORS['mid_gray'])
            for cat in compare['category']
        ]
        ax2.barh(compare['display_category'], compare['log2_odds_ratio'], color=colors, edgecolor='white', linewidth=0.8)
        style_axis(ax2, zero='x')
        ax2.set_xlabel('log2(OR) vs all proteome arginines')
        ax2.set_ylabel('Disorder class')
        ax2.set_title('Global methylarginine disorder-context enrichment')
        set_symmetric_xlim(ax2, compare['log2_odds_ratio'], annotation_pad_ratio=0.34, center_on_zero=False)
        annotate_barh(
            ax2,
            compare['display_category'],
            compare['log2_odds_ratio'],
            count_over_total_labels(compare),
            fontsize=9.0,
        )
        fig2.text(
            0.99,
            0.01,
            f'{q_threshold_note(compare["q_value"])}. Global class-level ORs; within-IDR edge bins are shown separately and use only disordered residues.',
            ha='right',
            va='bottom',
            fontsize=8,
            color=BREWER_COLORS['dark_gray'],
        )
        fig2.tight_layout()
        save_figure(fig2, outdir / 'arg_methyl_idr_context_enrichment')

        binary_plot = binary_all[binary_all['category'] == 'IDR'].copy()
        if not binary_plot.empty:
            fig3, ax3 = plt.subplots(figsize=(4.8, 4.8))
            ax3.bar(['IDR'], binary_plot['log2_odds_ratio'], color=BREWER_COLORS['mid_gray'], edgecolor='white', linewidth=0.8)
            style_axis(ax3, zero='y')
            ax3.set_ylabel('log2(OR) vs all proteome arginines')
            ax3.set_title('Methylarginine IDR enrichment')
            value = float(binary_plot['log2_odds_ratio'].iloc[0])
            qval = float(binary_plot['q_value'].iloc[0])
            label = count_over_total_labels(binary_plot)[0]
            ax3.text(0, value, label, ha='center', va='bottom' if value >= 0 else 'top', fontsize=9.0)
            fig3.text(0.99, 0.01, q_threshold_note([qval]), ha='right', va='bottom', fontsize=8, color=BREWER_COLORS['dark_gray'])
            fig3.tight_layout()
            save_figure(fig3, outdir / 'arg_methyl_idr_binary_enrichment')

        proximal_plot = proximal_binary_all[proximal_binary_all['category'] == 'IDR_proximal'].copy()
        if not proximal_plot.empty:
            fig3b, ax3b = plt.subplots(figsize=(5.8, 4.8))
            display_label = 'Disordered or outside-IDR edge (<=20 aa)'
            ax3b.bar([display_label], proximal_plot['log2_odds_ratio'], color=BREWER_COLORS['mid_gray'], edgecolor='white', linewidth=0.8)
            style_axis(ax3b, zero='y')
            ax3b.set_ylabel('log2(OR) vs all proteome arginines')
            ax3b.set_title('Methylarginine enrichment in IDR-proximal sequence')
            value = float(proximal_plot['log2_odds_ratio'].iloc[0])
            qval = float(proximal_plot['q_value'].iloc[0])
            label = count_over_total_labels(proximal_plot)[0]
            ax3b.text(0, value, label, ha='center', va='bottom' if value >= 0 else 'top', fontsize=9.0)
            ax3b.tick_params(axis='x', rotation=15)
            fig3b.text(
                0.99,
                0.01,
                f'{q_threshold_note([qval])}. IDR-proximal = inside an IDR or <=20 aa outside an IDR edge.',
                ha='right',
                va='bottom',
                fontsize=8,
                color=BREWER_COLORS['dark_gray'],
            )
            fig3b.tight_layout()
            save_figure(fig3b, outdir / 'arg_methyl_idr_proximal_binary_enrichment')

        edge_plot = edge_bin_same.sort_values('log2_odds_ratio')
        edge_palette = {
            'IDR edge <=5 aa': BREWER_COLORS['light_gray'],
            'IDR edge 6-10 aa': '#c9c9c9',
            'IDR edge 11-20 aa': BREWER_COLORS['mid_gray'],
            'IDR edge 21-40 aa': '#8d8d8d',
            'Distal IDR >40 aa': BREWER_COLORS['dark_gray'],
        }
        fig4, ax4 = plt.subplots(figsize=(8.3, 5.2))
        ax4.barh(
            edge_plot['category'],
            edge_plot['log2_odds_ratio'],
            color=[edge_palette.get(cat, BREWER_COLORS['mid_gray']) for cat in edge_plot['category']],
            edgecolor='white',
            linewidth=0.8,
        )
        style_axis(ax4, zero='x')
        ax4.set_xlabel('log2(OR) vs disordered arginines in methylated proteins')
        ax4.set_ylabel('Within-IDR edge-distance bin')
        ax4.set_title('Methylarginine localization within IDRs')
        set_symmetric_xlim(ax4, edge_plot['log2_odds_ratio'], annotation_pad_ratio=0.32, center_on_zero=False)
        annotate_barh(
            ax4,
            edge_plot['category'],
            edge_plot['log2_odds_ratio'],
            count_over_total_labels(edge_plot),
            fontsize=8.8,
        )
        fig4.text(
            0.99,
            0.01,
            q_threshold_note(edge_plot['q_value']),
            ha='right',
            va='bottom',
            fontsize=8,
            color=BREWER_COLORS['dark_gray'],
        )
        fig4.tight_layout()
        save_figure(fig4, outdir / 'arg_methyl_idr_edge_bin_enrichment')

        fig5, ax5 = plt.subplots(figsize=(8.8, 4.8))
        ax5.plot(edge_curve_same['distance_aa'], edge_curve_same['log2_odds_ratio'], color=BREWER_COLORS['dark_gray'], linewidth=2.2)
        sig = edge_curve_same[edge_curve_same['q_value'] <= 0.05]
        if not sig.empty:
            ax5.scatter(sig['distance_aa'], sig['log2_odds_ratio'], color='#111111', s=18, zorder=3)
        style_axis(ax5, zero='y', grid_axis='both')
        ax5.set_xlabel('Within-IDR distance from nearest edge (<= X aa)')
        ax5.set_ylabel('log2(OR) vs disordered arginines in methylated proteins')
        ax5.set_title('Cumulative enrichment within IDRs by edge distance')
        ax5.set_xlim(1, 40)
        for checkpoint in [5, 10, 20, 40]:
            row = edge_curve_same.loc[edge_curve_same['distance_aa'] == checkpoint]
            if row.empty:
                continue
            row = row.iloc[0]
            ax5.axvline(checkpoint, color=BREWER_COLORS['light_gray'], linestyle='--', linewidth=0.8)
            ax5.text(checkpoint, row['log2_odds_ratio'], f' {checkpoint}', fontsize=8, va='bottom', color=BREWER_COLORS['dark_gray'])
        fig5.text(
            0.99,
            0.01,
            'Disordered sites only; cumulative curve shown for positions within 40 aa of an IDR edge.',
            ha='right',
            va='bottom',
            fontsize=8,
            color=BREWER_COLORS['dark_gray'],
        )
        fig5.tight_layout()
        save_figure(fig5, outdir / 'arg_methyl_idr_edge_cumulative_enrichment')

        fig6, ax6 = plt.subplots(figsize=(8.8, 4.8))
        ax6.plot(edge_roll_same['window_center_aa'], edge_roll_same['log2_odds_ratio'], color=BREWER_COLORS['dark_gray'], linewidth=2.2)
        sig_roll = edge_roll_same[edge_roll_same['q_value'] <= 0.05]
        if not sig_roll.empty:
            ax6.scatter(sig_roll['window_center_aa'], sig_roll['log2_odds_ratio'], color='#111111', s=18, zorder=3)
        style_axis(ax6, zero='y', grid_axis='both')
        ax6.set_xlabel('Within-IDR distance from nearest edge (rolling 5-aa window)')
        ax6.set_ylabel('log2(OR) vs disordered arginines in methylated proteins')
        ax6.set_title('Rolling enrichment within IDRs by edge distance')
        ax6.set_xlim(3, 38)
        for checkpoint in [5, 10, 20, 40]:
            ax6.axvline(checkpoint, color=BREWER_COLORS['light_gray'], linestyle='--', linewidth=0.8)
        fig6.text(
            0.99,
            0.01,
            'Disordered sites only; each point is a 5-aa rolling window and points mark q<=0.05.',
            ha='right',
            va='bottom',
            fontsize=8,
            color=BREWER_COLORS['dark_gray'],
        )
        fig6.tight_layout()
        save_figure(fig6, outdir / 'arg_methyl_idr_edge_rolling_enrichment')

    print({
        'annotated_sites': len(annotated_sites),
        'sites_with_disorder_annotation': int(annotated_sites['has_disorder_annotation'].fillna(False).sum()),
        'output': str(outdir / 'sites_with_disorder_context.tsv'),
    })


if __name__ == '__main__':
    main()

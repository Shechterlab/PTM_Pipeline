from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from common import (
    BREWER_COLORS,
    add_q_values,
    annotate_barh,
    apply_paper_style,
    balanced_category_subset,
    compact_q_label,
    count_over_total_labels,
    q_threshold_note,
    save_figure,
    set_symmetric_xlim,
    signed_bar_colors,
    style_axis,
)


def domain_context_label(label: str) -> str:
    return {
        'in_domain': 'In domain/repeat',
        'boundary': 'Near domain edge (<=20 aa)',
        'inter_domain_linker': 'Inter-domain linker',
        'distal': 'Domain-distal',
        'no_domain_annotation': 'No annotated domain/repeat',
    }.get(str(label), str(label))


def domain_class_label(label: str) -> str:
    return {
        'No domain annotation': 'No annotated domain/repeat',
    }.get(str(label), str(label))


def format_domain_panels(results_root: Path) -> None:
    outdir = results_root / 'domain_context'
    apply_paper_style()

    context = pd.read_csv(outdir / 'domain_context_enrichment.tsv', sep='\t')
    plot_context = context[context['target_count'] >= 20].sort_values('log2_odds_ratio').copy()
    plot_context['display_category'] = plot_context['category'].map(domain_context_label)
    context_palette = {
        'in_domain': BREWER_COLORS['dark_gray'],
        'boundary': BREWER_COLORS['mid_gray'],
        'inter_domain_linker': '#8d8d8d',
        'distal': '#bdbdbd',
        'no_domain_annotation': BREWER_COLORS['light_gray'],
    }
    fig, ax = plt.subplots(figsize=(9.2, 5.8))
    ax.barh(
        plot_context['display_category'],
        plot_context['log2_odds_ratio'],
        color=[context_palette.get(cat, BREWER_COLORS['mid_gray']) for cat in plot_context['category']],
        edgecolor='white',
        linewidth=0.8,
    )
    style_axis(ax, zero='x')
    ax.set_xlabel('log2(OR) vs arginines in methylated proteins')
    ax.set_ylabel('Domain context')
    ax.set_title('Methylarginine enrichment by domain context')
    set_symmetric_xlim(ax, plot_context['log2_odds_ratio'], annotation_pad_ratio=0.38, center_on_zero=False)
    annotate_barh(
        ax,
        plot_context['display_category'],
        plot_context['log2_odds_ratio'],
        count_over_total_labels(plot_context),
        fontsize=8.9,
    )
    fig.text(
        0.99,
        0.01,
        f'{q_threshold_note(plot_context["q_value"])}. Site-level counts using strict InterPro domain/repeat intervals.',
        ha='right',
        va='bottom',
        fontsize=8,
        color=BREWER_COLORS['dark_gray'],
    )
    fig.tight_layout()
    save_figure(fig, outdir / 'arg_methyl_domain_context_enrichment')

    domain_class = pd.read_csv(outdir / 'domain_class_enrichment.tsv', sep='\t')
    plot_domain = balanced_category_subset(domain_class[domain_class['q_value'] <= 0.05].copy(), min_count=20, top_positive=7, top_negative=5)
    plot_domain = plot_domain[~plot_domain['category'].isin({'Other domain', 'No domain annotation'})].copy()
    plot_domain['display_category'] = plot_domain['category'].map(domain_class_label)
    colors = signed_bar_colors(plot_domain['log2_odds_ratio'])
    fig2, ax2 = plt.subplots(figsize=(9.4, 6.2))
    ax2.barh(plot_domain['display_category'], plot_domain['log2_odds_ratio'], color=colors, edgecolor='white', linewidth=0.8)
    style_axis(ax2, zero='x')
    ax2.set_xlabel('log2(OR) vs arginines in methylated proteins')
    ax2.set_ylabel('Nearest domain class')
    ax2.set_title('Nearest domain-class enrichment around methylarginines')
    set_symmetric_xlim(ax2, plot_domain['log2_odds_ratio'], annotation_pad_ratio=0.42, center_on_zero=False)
    annotate_barh(
        ax2,
        plot_domain['display_category'],
        plot_domain['log2_odds_ratio'],
        count_over_total_labels(plot_domain),
        fontsize=8.6,
    )
    fig2.text(
        0.99,
        0.01,
        f'{q_threshold_note(plot_domain["q_value"])}. Site-level counts; selected enriched/depleted nearest-domain classes.',
        ha='right',
        va='bottom',
        fontsize=8,
        color=BREWER_COLORS['dark_gray'],
    )
    fig2.tight_layout()
    save_figure(fig2, outdir / 'arg_methyl_nearest_domain_class_enrichment')

    fraction = pd.read_csv(outdir / 'domain_edge_bin_enrichment_disordered_only_compact_fraction.tsv', sep='\t')
    fraction['display_category'] = fraction['category']
    compact_palette = {
        'Domain edge <=20 aa': BREWER_COLORS['light_gray'],
        'Domain edge 21-40 aa': BREWER_COLORS['mid_gray'],
        'Domain-distal >40 aa': BREWER_COLORS['dark_gray'],
    }
    fig3, ax3 = plt.subplots(figsize=(7.8, 4.6))
    ax3.barh(
        fraction['display_category'],
        100.0 * fraction['methylated_fraction_in_bin'],
        color=[compact_palette.get(cat, BREWER_COLORS['mid_gray']) for cat in fraction['category']],
        edgecolor='white',
        linewidth=0.8,
    )
    style_axis(ax3)
    ax3.set_xlabel('% of arginines methylated in bin')
    ax3.set_ylabel('Outside-domain distance from nearest edge')
    ax3.set_title('Methylarginine frequency in disordered outside-domain bins')
    annotate_barh(
        ax3,
        fraction['display_category'],
        100.0 * fraction['methylated_fraction_in_bin'],
        [
            f"{n}/{t} ({100.0 * frac:.1f}%)"
            for n, t, frac in zip(
                fraction['target_count'],
                fraction['total_arginines_in_bin'],
                fraction['methylated_fraction_in_bin'],
            )
        ],
        fontsize=8.8,
    )
    fig3.text(
        0.99,
        0.01,
        'Disordered outside-domain arginines only; denominator = methylated + non-methyl arginines in each bin.',
        ha='right',
        va='bottom',
        fontsize=8,
        color=BREWER_COLORS['dark_gray'],
    )
    fig3.tight_layout()
    save_figure(fig3, outdir / 'arg_methyl_domain_edge_bin_frequency_disordered_only_compact')


def format_cross_ptm_panel(results_root: Path) -> None:
    outdir = results_root / 'ptm_compare'
    apply_paper_style()

    binary = pd.read_csv(outdir / 'ptm_idr_binary_comparison.tsv', sep='\t')
    counts = pd.read_csv(outdir / 'ptm_idr_binary_counts.tsv', sep='\t')
    plot_df = binary[binary['category'] == 'IDR'].copy().merge(counts, on='ptm_group', how='left')
    plot_df['observed_idr_fraction_pct'] = 100.0 * plot_df['idr_site_count'] / plot_df['classified_site_count']
    plot_df['background_idr_fraction_pct'] = 100.0 * plot_df['background_idr_site_count'] / plot_df['background_classified_site_count']
    plot_df = plot_df.sort_values('log2_odds_ratio')

    fig, ax = plt.subplots(figsize=(9.2, 5.6))
    colors = signed_bar_colors(
        plot_df['log2_odds_ratio'],
        positive=BREWER_COLORS['dark_gray'],
        negative=BREWER_COLORS['light_gray'],
    )
    ax.barh(plot_df['ptm_group'], plot_df['log2_odds_ratio'], color=colors, edgecolor='white', linewidth=0.8)
    style_axis(ax, zero='x')
    ax.set_xlabel('log2(OR) for IDR vs target-residue proteome background')
    ax.set_ylabel('PTM class')
    ax.set_title('Cross-PTM IDR enrichment')
    set_symmetric_xlim(ax, plot_df['log2_odds_ratio'], annotation_pad_ratio=0.34, center_on_zero=False)
    annotate_barh(
        ax,
        plot_df['ptm_group'],
        plot_df['log2_odds_ratio'],
        [f'{n}/{t}' for n, t in zip(plot_df['idr_site_count'], plot_df['classified_site_count'])],
        fontsize=8.7,
    )
    fig.text(0.99, 0.01, q_threshold_note(plot_df['q_value']), ha='right', va='bottom', fontsize=8, color=BREWER_COLORS['dark_gray'])
    fig.tight_layout()
    save_figure(fig, outdir / 'cross_ptm_idr_binary_comparison')

    fig2, ax2 = plt.subplots(figsize=(9.2, 5.6))
    ax2.barh(
        plot_df['ptm_group'],
        plot_df['observed_idr_fraction_pct'],
        color=BREWER_COLORS['mid_gray'],
        edgecolor='white',
        linewidth=0.8,
    )
    ax2.scatter(
        plot_df['background_idr_fraction_pct'],
        plot_df['ptm_group'],
        color=BREWER_COLORS['dark_gray'],
        marker='D',
        s=28,
        zorder=3,
        label='Target-residue proteome background',
    )
    style_axis(ax2)
    ax2.set_xlabel('% of classified sites in IDRs')
    ax2.set_ylabel('PTM class')
    ax2.set_title('Cross-PTM IDR site burden')
    ax2.legend(loc='lower right', frameon=False)
    annotate_barh(
        ax2,
        plot_df['ptm_group'],
        plot_df['observed_idr_fraction_pct'],
        [f'{n}/{t}' for n, t in zip(plot_df['idr_site_count'], plot_df['classified_site_count'])],
        fontsize=8.0,
    )
    fig2.tight_layout()
    save_figure(fig2, outdir / 'cross_ptm_idr_site_fraction')


def format_neighbor_summary(results_root: Path) -> None:
    outdir = results_root / 'methyl_arg_clustering'
    path = outdir / 'nearest_neighbor_distance_bin_summary.tsv'
    if not path.exists():
        return

    bins = pd.read_csv(path, sep='\t')
    if 'empirical_q_greater' not in bins.columns and 'empirical_p_greater' in bins.columns:
        bins = add_q_values(bins, p_col='empirical_p_greater', out_col='empirical_q_greater')

    plot_bins = bins[bins['empirical_q_greater'] <= 0.05].copy()
    if plot_bins.empty:
        return

    apply_paper_style()
    plot_bins = plot_bins.sort_values('ratio_vs_null', ascending=True).reset_index(drop=True)
    bar_colors = [BREWER_COLORS['teal'], BREWER_COLORS['blue']][-len(plot_bins):]
    fig, ax = plt.subplots(figsize=(5.2, 3.7))
    ax.barh(
        plot_bins['distance_bin'],
        plot_bins['ratio_vs_null'],
        color=bar_colors,
        edgecolor='white',
        linewidth=0.8,
    )
    style_axis(ax, grid_axis='x')
    ax.set_xlabel('Observed/null nearest-neighbor probability')
    ax.set_ylabel('Nearest methyl-Arg distance bin')
    ax.set_title('Short-range methyl-Arg clustering')
    for idx, row in plot_bins.iterrows():
        compare_text = f'{100.0 * float(row["observed_probability"]):.1f}% vs {100.0 * float(row["null_mean_probability"]):.1f}%'
        ax.text(
            float(row['ratio_vs_null']) + 0.05,
            idx,
            f'{compare_text}\n{compact_q_label(row["empirical_q_greater"])}',
            ha='left',
            va='center',
            fontsize=8.7,
            color=BREWER_COLORS['dark_gray'],
        )
    ax.set_xlim(0, max(2.6, float(plot_bins['ratio_vs_null'].max()) + 0.65))
    fig.text(
        0.99,
        0.01,
        'Only bins with BH q<0.05 are shown. Short-range clustering summary; full curve and null tables are provided separately.',
        ha='right',
        va='bottom',
        fontsize=8,
        color=BREWER_COLORS['dark_gray'],
    )
    fig.tight_layout()
    save_figure(fig, outdir / 'arg_methyl_neighbor_clustering_summary')

    checkpoint_path = outdir / 'nearest_neighbor_checkpoint_summary.tsv'
    if checkpoint_path.exists():
        checkpoints = pd.read_csv(checkpoint_path, sep='\t')
    else:
        checkpoint_curve = pd.read_csv(outdir / 'nearest_neighbor_curve_checkpoints.tsv', sep='\t')
        checkpoints = checkpoint_curve[checkpoint_curve['distance_aa'].isin([5, 10, 20])].copy()
        checkpoints = checkpoints.rename(columns={
            'prob_neighbor_within_x': 'observed_probability',
            'null_mean_prob_neighbor_within_x': 'null_mean_probability',
            'ratio_prob_neighbor_within_x_vs_null': 'ratio_vs_null',
            'empirical_p_greater_prob_neighbor_within_x': 'empirical_p_greater',
        })
        checkpoints['label'] = checkpoints['distance_aa'].map(lambda x: f'<={int(x)} aa')
        checkpoints = checkpoints[['distance_aa', 'label', 'sites_total', 'sites_with_neighbor_within_x', 'observed_probability', 'null_mean_probability', 'ratio_vs_null', 'empirical_p_greater']]
        checkpoints = add_q_values(checkpoints, p_col='empirical_p_greater', out_col='empirical_q_greater')
        checkpoints.to_csv(checkpoint_path, sep='\t', index=False)

    if not checkpoints.empty:
        fig2, ax2 = plt.subplots(figsize=(5.9, 4.1))
        x = pd.Series(range(len(checkpoints)), dtype=float).to_numpy()
        width = 0.34
        observed_pct = 100.0 * checkpoints['observed_probability'].to_numpy(dtype=float)
        null_pct = 100.0 * checkpoints['null_mean_probability'].to_numpy(dtype=float)
        ax2.bar(x - width / 2, observed_pct, width=width, color=BREWER_COLORS['blue'], edgecolor='white', linewidth=0.8, label='Observed')
        ax2.bar(x + width / 2, null_pct, width=width, color=BREWER_COLORS['light_gray'], edgecolor='white', linewidth=0.8, label='Within-protein null')
        style_axis(ax2, grid_axis='y')
        ax2.set_xticks(x, checkpoints['label'])
        ax2.set_xlabel('Nearest methyl-Arg distance checkpoint')
        ax2.set_ylabel('% of methyl sites with another methyl-Arg nearby')
        ax2.set_title('Methyl-Arg sites often occur near another methyl-Arg')
        ax2.legend(loc='upper right')
        ymax = max(float(observed_pct.max()), float(null_pct.max()))
        ax2.set_ylim(0, ymax + 12)
        for xpos, obs, qv in zip(x - width / 2, observed_pct, checkpoints['empirical_q_greater']):
            ax2.text(xpos, obs + 1.6, compact_q_label(qv), ha='center', va='bottom', fontsize=8.6, color=BREWER_COLORS['dark_gray'])
        fig2.text(
            0.99,
            0.01,
            'Cumulative checkpoint view using the same within-protein permutation null as the full spacing curve.',
            ha='right',
            va='bottom',
            fontsize=8,
            color=BREWER_COLORS['dark_gray'],
        )
        fig2.tight_layout()
        save_figure(fig2, outdir / 'arg_methyl_neighbor_checkpoint_summary')


def format_neighbor_probability(results_root: Path) -> None:
    outdir = results_root / 'methyl_arg_clustering'
    stats_path = outdir / 'nearest_neighbor_probability_curve_stats.tsv'
    path = stats_path if stats_path.exists() else (outdir / 'nearest_neighbor_curve_vs_null.tsv')
    if not path.exists():
        return

    curve = pd.read_csv(path, sep='\t')
    curve40 = curve[curve['distance_aa'] <= 40].copy()
    if curve40.empty:
        return

    apply_paper_style()
    fig, ax = plt.subplots(figsize=(3.8, 4.0))
    observed_col = 'observed_probability' if 'observed_probability' in curve40.columns else 'prob_neighbor_within_x'
    null_mean_col = 'null_mean_probability' if 'null_mean_probability' in curve40.columns else 'null_mean_prob_neighbor_within_x'
    null_lo_col = 'null_q025_probability' if 'null_q025_probability' in curve40.columns else 'null_q025_prob_neighbor_within_x'
    null_hi_col = 'null_q975_probability' if 'null_q975_probability' in curve40.columns else 'null_q975_prob_neighbor_within_x'
    observed_pct = 100.0 * curve40[observed_col].to_numpy(dtype=float)
    null_pct = 100.0 * curve40[null_mean_col].to_numpy(dtype=float)
    null_lo_pct = 100.0 * curve40[null_lo_col].to_numpy(dtype=float)
    null_hi_pct = 100.0 * curve40[null_hi_col].to_numpy(dtype=float)
    distance = curve40['distance_aa'].to_numpy(dtype=float)

    ax.fill_between(distance, null_lo_pct, null_hi_pct, alpha=0.14, color=BREWER_COLORS['light_gray'])
    ax.plot(distance, null_pct, linewidth=2.0, linestyle='--', color=BREWER_COLORS['mid_gray'])
    ax.plot(distance, observed_pct, linewidth=2.4, color=BREWER_COLORS['blue'])
    style_axis(ax, grid_axis='y')
    ax.set_xlim(1, 40)
    ax.set_xticks([10, 20, 30, 40])
    ax.set_xlabel('Within X aa')
    ax.set_ylabel('% of methyl-Arg sites with another nearby')
    ax.set_title('Nearest methyl-Arg neighbor probability')

    end_x = float(distance[-1])
    ax.text(end_x + 0.5, float(observed_pct[-1]), 'Observed', color=BREWER_COLORS['blue'], fontsize=9.0, va='center', ha='left')
    ax.text(end_x + 0.5, float(null_pct[-1]), 'Within-protein null', color=BREWER_COLORS['dark_gray'], fontsize=9.0, va='center', ha='left')

    ax.set_ylim(0, max(62.0, float(observed_pct.max()) + 4.0))
    fig.tight_layout()
    save_figure(fig, outdir / 'nearest_neighbor_probability_curve')


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--results-root', default='PTM_results')
    args = ap.parse_args()

    results_root = Path(args.results_root)
    format_domain_panels(results_root)
    format_cross_ptm_panel(results_root)
    format_neighbor_probability(results_root)
    format_neighbor_summary(results_root)


if __name__ == '__main__':
    main()

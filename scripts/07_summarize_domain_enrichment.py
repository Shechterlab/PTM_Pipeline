from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError

from common import (
    classify_by_patterns,
    fisher_like_enrichment,
    format_p_value,
    load_json,
    parse_fasta,
    residue_background_from_sequences,
    save_figure,
    save_table,
)


TYPE_PRIORITY = {
    'domain': 0,
    'repeat': 1,
    'homologous_superfamily': 2,
    'family': 3,
}


EMPTY_INTERVAL_COLUMNS = [
    'canonical_UniProtAC', 'interpro_accession', 'interpro_name', 'interpro_type',
    'member_database_accessions', 'fragment_start', 'fragment_end',
]


def interval_priority(interpro_type: str) -> int:
    return TYPE_PRIORITY.get(str(interpro_type), 9)


def load_intervals(path: str) -> pd.DataFrame:
    try:
        intervals = pd.read_csv(path, sep='\t', low_memory=False)
    except EmptyDataError:
        intervals = pd.DataFrame(columns=EMPTY_INTERVAL_COLUMNS)
    if intervals.empty:
        return pd.DataFrame(columns=EMPTY_INTERVAL_COLUMNS)
    for col in ['interpro_accession', 'interpro_name', 'interpro_type', 'member_database_accessions']:
        if col not in intervals.columns:
            intervals[col] = ''
    required = ['canonical_UniProtAC', 'fragment_start', 'fragment_end']
    if any(col not in intervals.columns for col in required):
        return pd.DataFrame(columns=EMPTY_INTERVAL_COLUMNS)
    intervals['fragment_start'] = pd.to_numeric(intervals['fragment_start'], errors='coerce')
    intervals['fragment_end'] = pd.to_numeric(intervals['fragment_end'], errors='coerce')
    intervals = intervals[intervals['fragment_start'].notna() & intervals['fragment_end'].notna()].copy()
    if intervals.empty:
        return pd.DataFrame(columns=EMPTY_INTERVAL_COLUMNS)
    intervals['fragment_start'] = intervals['fragment_start'].astype(int)
    intervals['fragment_end'] = intervals['fragment_end'].astype(int)
    return intervals


def edge_distance(position: int, start: int, end: int) -> tuple[int, str]:
    if start <= position <= end:
        return 0, 'inside'
    if position < start:
        return start - position, 'N_terminal_to_domain'
    return position - end, 'C_terminal_to_domain'


def choose_interval(intervals: pd.DataFrame, position: int) -> tuple[pd.Series | None, int]:
    if intervals.empty:
        return None, -1
    overlaps = intervals[(intervals['fragment_start'] <= position) & (intervals['fragment_end'] >= position)].copy()
    if not overlaps.empty:
        overlaps['span'] = overlaps['fragment_end'] - overlaps['fragment_start']
        overlaps['priority'] = overlaps['interpro_type'].map(interval_priority)
        overlaps = overlaps.sort_values(['priority', 'span', 'fragment_start', 'interpro_accession'])
        return overlaps.iloc[0], 0
    candidates = intervals.copy()
    candidates['distance'] = candidates.apply(lambda r: edge_distance(position, int(r['fragment_start']), int(r['fragment_end']))[0], axis=1)
    candidates['span'] = candidates['fragment_end'] - candidates['fragment_start']
    candidates['priority'] = candidates['interpro_type'].map(interval_priority)
    candidates = candidates.sort_values(['distance', 'priority', 'span', 'fragment_start', 'interpro_accession'])
    best = candidates.iloc[0]
    return best, int(best['distance'])


def is_inter_domain_linker(intervals: pd.DataFrame, position: int) -> bool:
    if len(intervals) < 2:
        return False
    ordered = intervals[['fragment_start', 'fragment_end']].sort_values(['fragment_start', 'fragment_end']).drop_duplicates()
    previous_end = None
    for row in ordered.itertuples(index=False):
        start = int(row.fragment_start)
        end = int(row.fragment_end)
        if previous_end is not None and previous_end < position < start:
            return True
        previous_end = max(previous_end or end, end)
    return False


def classify_residue(intervals_by_acc: dict[str, pd.DataFrame], accession: str, position: int, boundary_window: int, ontology: dict) -> dict:
    intervals = intervals_by_acc.get(accession, pd.DataFrame())
    if intervals.empty:
        return {
            'domain_context_class': 'no_domain_annotation',
            'nearest_domain_class': 'No domain annotation',
            'nearest_domain_edge_distance': pd.NA,
            'is_within_5aa_domain_edge': False,
            'is_within_10aa_domain_edge': False,
            'is_within_20aa_domain_edge': False,
        }
    best, distance = choose_interval(intervals, position)
    if best is None:
        return {
            'domain_context_class': 'no_domain_annotation',
            'nearest_domain_class': 'No domain annotation',
            'nearest_domain_edge_distance': pd.NA,
            'is_within_5aa_domain_edge': False,
            'is_within_10aa_domain_edge': False,
            'is_within_20aa_domain_edge': False,
        }
    if distance == 0:
        context_class = 'in_domain'
    elif distance <= boundary_window:
        context_class = 'boundary'
    elif is_inter_domain_linker(intervals, position):
        context_class = 'inter_domain_linker'
    else:
        context_class = 'distal'
    label_text = ' '.join([
        str(best.get('interpro_name', '')),
        str(best.get('interpro_type', '')),
        str(best.get('member_database_accessions', '')),
    ])
    return {
        'domain_context_class': context_class,
        'nearest_domain_class': classify_by_patterns(label_text, ontology['domain_classes'], default='Other domain'),
        'nearest_domain_edge_distance': distance,
        'is_within_5aa_domain_edge': bool(distance >= 0 and distance <= 5),
        'is_within_10aa_domain_edge': bool(distance >= 0 and distance <= 10),
        'is_within_20aa_domain_edge': bool(distance >= 0 and distance <= 20),
    }


def annotate_background(background: pd.DataFrame, intervals_by_acc: dict[str, pd.DataFrame], boundary_window: int, ontology: dict) -> pd.DataFrame:
    if background.empty:
        return background.copy()
    annotations = background.apply(
        lambda row: pd.Series(
            classify_residue(
                intervals_by_acc=intervals_by_acc,
                accession=str(row['canonical_UniProtAC']),
                position=int(row['position']),
                boundary_window=boundary_window,
                ontology=ontology,
            )
        ),
        axis=1,
    )
    return pd.concat([background.reset_index(drop=True), annotations], axis=1)


def matched_background_weights(background: pd.DataFrame, observed: pd.DataFrame) -> pd.DataFrame:
    bg_counts = background.groupby('canonical_UniProtAC').size().rename('background_arg_count')
    obs_counts = observed.groupby('canonical_UniProtAC').size().rename('methyl_site_count')
    weights = bg_counts.to_frame().merge(obs_counts, how='left', left_index=True, right_index=True).fillna(0)
    weights['matched_weight'] = weights.apply(
        lambda row: float(row['methyl_site_count']) / float(row['background_arg_count']) if row['background_arg_count'] > 0 else 0.0,
        axis=1,
    )
    out = background.merge(weights[['matched_weight']], how='left', left_on='canonical_UniProtAC', right_index=True)
    out['matched_weight'] = out['matched_weight'].fillna(0.0)
    return out


def weighted_category_summary(observed: pd.DataFrame, background: pd.DataFrame, category_col: str) -> pd.DataFrame:
    obs_counts = observed[category_col].fillna('NA').value_counts().rename_axis('category').reset_index(name='observed_count')
    bg = background.copy()
    bg[category_col] = bg[category_col].fillna('NA')
    bg_counts = bg.groupby(category_col)['matched_weight'].sum().rename_axis('category').reset_index(name='matched_control_weight')
    merged = obs_counts.merge(bg_counts, how='outer', on='category').fillna(0)
    merged['observed_fraction'] = merged['observed_count'] / max(float(len(observed)), 1.0)
    total_bg = max(float(merged['matched_control_weight'].sum()), 1.0)
    merged['matched_control_fraction'] = merged['matched_control_weight'] / total_bg
    merged['fraction_ratio_observed_vs_matched'] = merged['observed_fraction'] / merged['matched_control_fraction'].replace(0, np.nan)
    return merged.sort_values(['observed_count', 'fraction_ratio_observed_vs_matched'], ascending=[False, False])


def cumulative_edge_table(observed: pd.DataFrame, background: pd.DataFrame, cutoffs: list[int]) -> pd.DataFrame:
    rows = []
    obs_dist = pd.to_numeric(observed['nearest_domain_edge_distance'], errors='coerce')
    bg_dist = pd.to_numeric(background['nearest_domain_edge_distance'], errors='coerce')
    total_bg = max(float(background['matched_weight'].sum()), 1.0)
    for cutoff in cutoffs:
        obs_frac = float(((obs_dist >= 0) & (obs_dist <= cutoff)).mean()) if len(obs_dist) else np.nan
        bg_mask = (bg_dist >= 0) & (bg_dist <= cutoff)
        bg_weight = float(background.loc[bg_mask.fillna(False), 'matched_weight'].sum()) if len(background) else np.nan
        bg_frac = bg_weight / total_bg if np.isfinite(bg_weight) else np.nan
        rows.append({
            'distance_cutoff_aa': cutoff,
            'observed_fraction_within_cutoff': obs_frac,
            'matched_control_fraction_within_cutoff': bg_frac,
            'fraction_ratio_observed_vs_matched': obs_frac / bg_frac if pd.notna(bg_frac) and bg_frac != 0 else np.nan,
        })
    return pd.DataFrame(rows)


def edge_distance_distribution(observed: pd.DataFrame, background: pd.DataFrame, max_distance: int = 100) -> pd.DataFrame:
    obs_dist = pd.to_numeric(observed['nearest_domain_edge_distance'], errors='coerce')
    bg_dist = pd.to_numeric(background['nearest_domain_edge_distance'], errors='coerce')
    bg_total = max(float(background['matched_weight'].sum()), 1.0)
    rows = []
    for dist in range(0, max_distance + 1):
        rows.append({
            'distance_aa': dist,
            'observed_fraction_at_distance': float((obs_dist == dist).mean()) if len(obs_dist) else np.nan,
            'matched_control_fraction_at_distance': float(background.loc[(bg_dist == dist).fillna(False), 'matched_weight'].sum()) / bg_total if len(background) else np.nan,
        })
    out = pd.DataFrame(rows)
    out['observed_over_matched'] = out['observed_fraction_at_distance'] / out['matched_control_fraction_at_distance'].replace(0, np.nan)
    return out


def plot_bar_enrichment(df: pd.DataFrame, outpath: Path) -> None:
    plot_df = df[df['category'].isin(['in_domain', 'boundary', 'inter_domain_linker', 'distal'])].sort_values('log2_odds_ratio')
    if plot_df.empty:
        return
    fig, ax = plt.subplots(figsize=(8.3, 5.0))
    ax.barh(plot_df['category'], plot_df['log2_odds_ratio'])
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlabel('log2(OR) vs non-methyl arginines in same proteins')
    ax.set_ylabel('Domain context')
    ax.set_title('Methylarginine enrichment by domain context')
    for y, v, n, q in zip(plot_df['category'], plot_df['log2_odds_ratio'], plot_df['target_count'], plot_df['q_value_bh']):
        ax.text(v, y, f'  n={n}, q={format_p_value(q)}', va='center', ha='left' if v >= 0 else 'right', fontsize=9)
    fig.tight_layout()
    save_figure(fig, outpath)


def plot_domain_class(df: pd.DataFrame, outpath: Path) -> None:
    plot_df = df[(df['target_count'] >= 15) & (df['category'] != 'No domain annotation')].sort_values(['q_value_bh', 'odds_ratio', 'target_count'], ascending=[True, False, False]).head(12)
    if plot_df.empty:
        return
    plot_df = plot_df.sort_values('log2_odds_ratio')
    fig, ax = plt.subplots(figsize=(9.7, max(5.4, 0.42 * len(plot_df) + 1.7)))
    ax.barh(plot_df['category'], plot_df['log2_odds_ratio'])
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlabel('log2(OR) vs non-methyl arginines in same proteins')
    ax.set_ylabel('Nearest domain class')
    ax.set_title('Nearest domain-class enrichment around methylarginines')
    for y, v, n, q in zip(plot_df['category'], plot_df['log2_odds_ratio'], plot_df['target_count'], plot_df['q_value_bh']):
        ax.text(v, y, f'  n={n}, q={format_p_value(q)}', va='center', ha='left' if v >= 0 else 'right', fontsize=8.5)
    fig.tight_layout()
    save_figure(fig, outpath)


def plot_edge_curve(edge_table: pd.DataFrame, outpath: Path) -> None:
    if edge_table.empty:
        return
    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    ax.plot(edge_table['distance_cutoff_aa'], edge_table['observed_fraction_within_cutoff'], linewidth=2.2, label='Methyl-Arg sites')
    ax.plot(edge_table['distance_cutoff_aa'], edge_table['matched_control_fraction_within_cutoff'], linewidth=2.0, linestyle='--', label='Matched non-methyl Arg controls')
    ax.set_xlabel('Distance to nearest domain edge cutoff (aa)')
    ax.set_ylabel('Fraction of sites within cutoff')
    ax.set_title('Methylarginines accumulate near domain boundaries')
    ax.legend()
    fig.tight_layout()
    save_figure(fig, outpath)


def plot_matched_ratio(summary: pd.DataFrame, outpath: Path, title: str, xlabel: str) -> None:
    plot_df = summary.copy()
    plot_df = plot_df[~plot_df['category'].isin(['NA', 'No domain annotation'])]
    if plot_df.empty:
        return
    plot_df = plot_df.sort_values('fraction_ratio_observed_vs_matched')
    fig, ax = plt.subplots(figsize=(8.0, max(4.8, 0.38 * len(plot_df) + 1.7)))
    ax.barh(plot_df['category'], plot_df['fraction_ratio_observed_vs_matched'])
    ax.axvline(1.0, color='black', linewidth=0.8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Category')
    ax.set_title(title)
    fig.tight_layout()
    save_figure(fig, outpath)


def plot_edge_distance_distribution(edge_dist: pd.DataFrame, outpath: Path) -> None:
    if edge_dist.empty:
        return
    fig, ax = plt.subplots(figsize=(8.8, 5.1))
    ax.plot(edge_dist['distance_aa'], edge_dist['observed_fraction_at_distance'], linewidth=2.0, label='Methyl-Arg sites')
    ax.plot(edge_dist['distance_aa'], edge_dist['matched_control_fraction_at_distance'], linewidth=2.0, linestyle='--', label='Matched non-methyl Arg controls')
    ax.set_xlabel('Exact distance to nearest domain edge (aa)')
    ax.set_ylabel('Fraction of sites')
    ax.set_title('Domain-edge distance distribution for methylarginines')
    ax.legend()
    fig.tight_layout()
    save_figure(fig, outpath)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--annotated-sites', required=True)
    ap.add_argument('--interpro-intervals', required=True)
    ap.add_argument('--canonical-fasta', required=True)
    ap.add_argument('--ontology', default='config/domain_class_ontology.json')
    ap.add_argument('--boundary-window', type=int, default=20)
    ap.add_argument('--outdir', default='results/domain_context')
    args = ap.parse_args()

    annotated_sites = pd.read_csv(args.annotated_sites, sep='\t', low_memory=False)
    intervals = load_intervals(args.interpro_intervals)
    intervals_by_acc = {
        accession: group.sort_values(['fragment_start', 'fragment_end', 'interpro_accession']).reset_index(drop=True)
        for accession, group in intervals.groupby('canonical_UniProtAC')
    } if not intervals.empty else {}
    ontology = load_json(args.ontology)

    if 'nearest_domain_class' not in annotated_sites.columns:
        text_cols = [c for c in ['nearest_interpro_name', 'nearest_interpro_type', 'nearest_member_database_accessions'] if c in annotated_sites.columns]
        if text_cols:
            class_text = annotated_sites.fillna('').astype(str)[text_cols].agg(' '.join, axis=1)
            annotated_sites['nearest_domain_class'] = class_text.map(lambda x: classify_by_patterns(x, ontology['domain_classes'], default='Other domain'))
        else:
            annotated_sites['nearest_domain_class'] = 'Other domain'

    sequences = parse_fasta(args.canonical_fasta)
    proteome_arg = residue_background_from_sequences(sequences, residue='R')
    methylated_proteins = sorted(set(annotated_sites['canonical_UniProtAC'].dropna().astype(str)))
    methyl_protein_arg = proteome_arg[proteome_arg['canonical_UniProtAC'].isin(methylated_proteins)].copy()
    methyl_site_keys = set(annotated_sites['site_key'].astype(str)) if 'site_key' in annotated_sites.columns else set(annotated_sites['canonical_UniProtAC'].astype(str) + ':' + annotated_sites['corrected_position'].astype(str))
    if 'site_key' not in methyl_protein_arg.columns:
        methyl_protein_arg['site_key'] = methyl_protein_arg['canonical_UniProtAC'].astype(str) + ':' + methyl_protein_arg['site'].astype(str)
    nonmethyl_arg = methyl_protein_arg[~methyl_protein_arg['site_key'].isin(methyl_site_keys)].copy()
    nonmethyl_arg = annotate_background(nonmethyl_arg, intervals_by_acc, args.boundary_window, ontology)
    weighted_background = matched_background_weights(nonmethyl_arg, annotated_sites)

    context_enrichment = fisher_like_enrichment(annotated_sites['domain_context_class'], nonmethyl_arg['domain_context_class']) if not nonmethyl_arg.empty else pd.DataFrame()
    domain_class_enrichment = fisher_like_enrichment(annotated_sites['nearest_domain_class'], nonmethyl_arg['nearest_domain_class']) if not nonmethyl_arg.empty else pd.DataFrame()
    edge_summary = annotated_sites[['nearest_domain_edge_distance', 'domain_context_class']].copy()
    matched_context = weighted_category_summary(annotated_sites, weighted_background, 'domain_context_class') if not weighted_background.empty else pd.DataFrame()
    matched_domain_class = weighted_category_summary(annotated_sites, weighted_background, 'nearest_domain_class') if not weighted_background.empty else pd.DataFrame()
    edge_curve = cumulative_edge_table(annotated_sites, weighted_background, cutoffs=[0, 1, 2, 3, 5, 10, 15, 20, 30, 50, 75, 100]) if not weighted_background.empty else pd.DataFrame()
    edge_dist = edge_distance_distribution(annotated_sites, weighted_background, max_distance=100) if not weighted_background.empty else pd.DataFrame()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(nonmethyl_arg, outdir / 'nonmethyl_arginines_with_domain_context.tsv')
    save_table(weighted_background, outdir / 'matched_nonmethyl_arginines_with_weights.tsv')
    save_table(context_enrichment, outdir / 'domain_context_enrichment.tsv')
    save_table(domain_class_enrichment, outdir / 'domain_class_enrichment.tsv')
    save_table(edge_summary, outdir / 'domain_edge_distance_summary.tsv')
    save_table(matched_context, outdir / 'domain_context_matched_control_summary.tsv')
    save_table(matched_domain_class, outdir / 'domain_class_matched_control_summary.tsv')
    save_table(edge_curve, outdir / 'domain_edge_distance_cumulative_vs_matched_control.tsv')
    save_table(edge_dist, outdir / 'domain_edge_distance_distribution_vs_matched_control.tsv')

<<<<<<< HEAD
    plot_bar_enrichment(context_enrichment, outdir / 'arg_methyl_domain_context_enrichment')
    plot_domain_class(domain_class_enrichment, outdir / 'arg_methyl_nearest_domain_class_enrichment')
    plot_edge_curve(edge_curve, outdir / 'arg_methyl_domain_edge_proximity_curve')
    plot_matched_ratio(matched_context, outdir / 'arg_methyl_domain_context_matched_control_shift', 'Matched-control domain-context shift', 'Observed / matched-control fraction')
    plot_matched_ratio(matched_domain_class, outdir / 'arg_methyl_domain_class_matched_control_shift', 'Matched-control nearest-domain-class shift', 'Observed / matched-control fraction')
    plot_edge_distance_distribution(edge_dist, outdir / 'arg_methyl_domain_edge_distance_distribution')
=======
    plot_context = context_enrichment[context_enrichment['category'].isin(['in_domain', 'boundary', 'inter_domain_linker', 'distal'])].sort_values('odds_ratio')
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    ax.barh(plot_context['category'], plot_context['log2_odds_ratio'])
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlabel('log2(OR) vs arginines in methylated proteins')
    ax.set_ylabel('Domain context')
    ax.set_title('Methylarginine enrichment by domain context')
    for y, v, n, q in zip(plot_context['category'], plot_context['log2_odds_ratio'], plot_context['target_count'], plot_context['q_value_bh']):
        ax.text(v, y, f'  n={n}, q={format_p_value(q)}', va='center', ha='left' if v >= 0 else 'right', fontsize=9)
    fig.tight_layout()
    save_figure(fig, outdir / 'arg_methyl_domain_context_enrichment')

    plot_domain = domain_class_enrichment[domain_class_enrichment['target_count'] >= 20].head(12).sort_values('odds_ratio')
    fig2, ax2 = plt.subplots(figsize=(9.0, 6.0))
    ax2.barh(plot_domain['category'], plot_domain['log2_odds_ratio'])
    ax2.axvline(0, color='black', linewidth=0.8)
    ax2.set_xlabel('log2(OR) vs arginines in methylated proteins')
    ax2.set_ylabel('Nearest domain class')
    ax2.set_title('Nearest domain-class enrichment around methylarginines')
    for y, v, n, q in zip(plot_domain['category'], plot_domain['log2_odds_ratio'], plot_domain['target_count'], plot_domain['q_value_bh']):
        ax2.text(v, y, f'  n={n}, q={format_p_value(q)}', va='center', ha='left' if v >= 0 else 'right', fontsize=9)
    fig2.tight_layout()
    save_figure(fig2, outdir / 'arg_methyl_nearest_domain_class_enrichment')
>>>>>>> 2c1a8ff7d8603628918f8ffc773cdcc98c903bff

    print({
        'annotated_sites': len(annotated_sites),
        'background_arginines_same_proteins': len(nonmethyl_arg),
        'matched_control_weight_total': float(weighted_background['matched_weight'].sum()) if not weighted_background.empty else 0.0,
        'outputs': str(outdir),
    })


if __name__ == '__main__':
    main()

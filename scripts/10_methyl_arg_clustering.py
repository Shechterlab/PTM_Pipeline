from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import (
    BREWER_COLORS,
    add_q_values,
    apply_paper_style,
    canonical_site_table,
    format_p_value,
    parse_fasta,
    save_figure,
    save_table,
    style_axis,
)
from clustering_utils import (
    add_empirical_p_values,
    aggregate_null_curve,
    checkpoint_summary,
    histograms_from_positions,
    merged_intervals_by_accession,
    normalize_disorder_intervals,
    permutation_summary,
    permutation_summary_grouped,
    position_is_in_idr,
    protein_density_summary,
    residue_positions,
    summarize_curve_from_histograms,
)


def validate_sites(sites: pd.DataFrame, sequences: dict[str, str]) -> pd.DataFrame:
    out = canonical_site_table(sites, accession_col='canonical_UniProtAC', position_col='corrected_position', residue_col='residue')
    out = out[out['residue'] == 'R'].copy()
    out = out[out['canonical_UniProtAC'].isin(sequences)].copy()

    def matches_sequence(row: pd.Series) -> bool:
        seq = sequences.get(str(row['canonical_UniProtAC']), '')
        pos = int(row['corrected_position'])
        return bool(seq) and 1 <= pos <= len(seq) and seq[pos - 1] == 'R'

    out = out[out.apply(matches_sequence, axis=1)].copy()
    return out.drop_duplicates(['canonical_UniProtAC', 'corrected_position']).reset_index(drop=True)


def annotate_idr_binary(df: pd.DataFrame, intervals_by_acc: dict[str, pd.DataFrame]) -> pd.DataFrame:
    out = df.copy()
    out['idr_binary_class'] = 'non_IDR'
    for accession, idx in out.groupby('canonical_UniProtAC').groups.items():
        positions = out.loc[idx, 'corrected_position'].to_numpy(dtype=int)
        flags = position_is_in_idr(positions, intervals_by_acc.get(accession, pd.DataFrame(columns=['fragment_start', 'fragment_end'])))
        out.loc[idx, 'idr_binary_class'] = np.where(flags, 'IDR', 'non_IDR')
    return out


def distance_bin_summary(curve: pd.DataFrame, permutation_df: pd.DataFrame, thresholds: list[int]) -> pd.DataFrame:
    thresholds = sorted(set(int(x) for x in thresholds if int(x) > 0))
    rows = []
    previous = 0
    for threshold in thresholds:
        obs_hi = float(curve.loc[curve['distance_aa'] == threshold, 'prob_neighbor_within_x'].iloc[0])
        obs_lo = float(curve.loc[curve['distance_aa'] == previous, 'prob_neighbor_within_x'].iloc[0]) if previous > 0 else 0.0
        obs_bin = obs_hi - obs_lo

        null_group = permutation_df[permutation_df['distance_aa'].isin(([previous] if previous > 0 else []) + [threshold])].copy()
        if previous > 0:
            wide = null_group.pivot(index='permutation', columns='distance_aa', values='prob_neighbor_within_x').fillna(0.0)
            null_bin = wide[threshold] - wide[previous]
        else:
            wide = null_group.pivot(index='permutation', columns='distance_aa', values='prob_neighbor_within_x').fillna(0.0)
            null_bin = wide[threshold]
        label = f'<={threshold} aa' if previous == 0 else f'{previous + 1}-{threshold} aa'
        p_tail = (np.count_nonzero(null_bin.to_numpy(dtype=float) >= obs_bin) + 1) / (len(null_bin) + 1)
        rows.append({
            'distance_bin': label,
            'bin_start_exclusive': previous,
            'bin_end_inclusive': threshold,
            'observed_probability': obs_bin,
            'null_mean_probability': float(null_bin.mean()),
            'null_q025_probability': float(null_bin.quantile(0.025)),
            'null_q975_probability': float(null_bin.quantile(0.975)),
            'ratio_vs_null': obs_bin / float(null_bin.mean()) if float(null_bin.mean()) > 0 else np.nan,
            'log2_ratio_vs_null': np.log2(obs_bin / float(null_bin.mean())) if float(null_bin.mean()) > 0 and obs_bin > 0 else np.nan,
            'empirical_p_greater': p_tail,
        })
        previous = threshold

    obs_prev = float(curve.loc[curve['distance_aa'] == thresholds[-1], 'prob_neighbor_within_x'].iloc[0])
    obs_tail = 1.0 - obs_prev
    tail_group = permutation_df[permutation_df['distance_aa'] == thresholds[-1]].copy()
    null_tail = 1.0 - tail_group['prob_neighbor_within_x'].astype(float)
    p_tail = (np.count_nonzero(null_tail.to_numpy(dtype=float) >= obs_tail) + 1) / (len(null_tail) + 1)
    rows.append({
        'distance_bin': f'>{thresholds[-1]} aa',
        'bin_start_exclusive': thresholds[-1],
        'bin_end_inclusive': pd.NA,
        'observed_probability': obs_tail,
        'null_mean_probability': float(null_tail.mean()),
        'null_q025_probability': float(null_tail.quantile(0.025)),
        'null_q975_probability': float(null_tail.quantile(0.975)),
        'ratio_vs_null': obs_tail / float(null_tail.mean()) if float(null_tail.mean()) > 0 else np.nan,
        'log2_ratio_vs_null': np.log2(obs_tail / float(null_tail.mean())) if float(null_tail.mean()) > 0 and obs_tail > 0 else np.nan,
        'empirical_p_greater': p_tail,
    })
    return pd.DataFrame(rows)


def rolling_probability_ratio(curve: pd.DataFrame, permutation_df: pd.DataFrame, window_size: int = 5, max_distance: int = 40) -> pd.DataFrame:
    rows = []
    for start in range(1, max_distance - window_size + 2):
        end = start + window_size - 1
        obs_hi = float(curve.loc[curve['distance_aa'] == end, 'prob_neighbor_within_x'].iloc[0])
        obs_lo = float(curve.loc[curve['distance_aa'] == start - 1, 'prob_neighbor_within_x'].iloc[0]) if start > 1 else 0.0
        obs_window = obs_hi - obs_lo

        needed = [end] + ([start - 1] if start > 1 else [])
        null_group = permutation_df[permutation_df['distance_aa'].isin(needed)].copy()
        wide = null_group.pivot(index='permutation', columns='distance_aa', values='prob_neighbor_within_x').fillna(0.0)
        null_window = wide[end] - wide[start - 1] if start > 1 else wide[end]
        null_mean = float(null_window.mean())
        p_tail = (np.count_nonzero(null_window.to_numpy(dtype=float) >= obs_window) + 1) / (len(null_window) + 1)

        rows.append({
            'window_start_aa': start,
            'window_end_aa': end,
            'window_center_aa': (start + end) / 2.0,
            'window_label': f'{start}-{end} aa',
            'observed_probability': obs_window,
            'null_mean_probability': null_mean,
            'null_q025_probability': float(null_window.quantile(0.025)),
            'null_q975_probability': float(null_window.quantile(0.975)),
            'ratio_vs_null': obs_window / null_mean if null_mean > 0 else np.nan,
            'log2_ratio_vs_null': np.log2(obs_window / null_mean) if null_mean > 0 and obs_window > 0 else np.nan,
            'empirical_p_greater': p_tail,
        })
    return add_q_values(pd.DataFrame(rows), p_col='empirical_p_greater', out_col='empirical_q_greater')


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--integrated-sites', required=True)
    ap.add_argument('--canonical-fasta', required=True)
    ap.add_argument('--outdir', default='results/methyl_arg_clustering')
    ap.add_argument('--max-distance', type=int, default=100)
    ap.add_argument('--permutations', type=int, default=100000)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--progress-every', type=int, default=5000)
    ap.add_argument('--checkpoints', nargs='*', type=int, default=[3, 5, 10, 20, 50])
    ap.add_argument('--summary-checkpoints', nargs='*', type=int, default=[5, 10, 20, 40])
    ap.add_argument('--disorder-intervals', default='')
    ap.add_argument('--null-match-disorder', action='store_true')
    ap.add_argument('--restrict-context', choices=['all', 'IDR', 'non_IDR'], default='all')
    args = ap.parse_args()

    sites = pd.read_csv(args.integrated_sites, sep='\t', low_memory=False)
    input_site_count = len(sites)
    sequences = parse_fasta(args.canonical_fasta)
    sites = validate_sites(sites, sequences=sequences)
    apply_paper_style()

    intervals_by_acc = None
    if args.disorder_intervals:
        disorder = normalize_disorder_intervals(pd.read_csv(args.disorder_intervals, sep='\t', low_memory=False))
        intervals_by_acc = merged_intervals_by_accession(disorder)
        sites = annotate_idr_binary(sites, intervals_by_acc)
        if args.restrict_context != 'all':
            sites = sites[sites['idr_binary_class'] == args.restrict_context].copy()
    elif args.null_match_disorder or args.restrict_context != 'all':
        raise ValueError('--disorder-intervals is required for --null-match-disorder or --restrict-context')

    methylated_proteins = sorted(set(sites['canonical_UniProtAC'].astype(str)))
    arg_positions_by_acc = {}
    arg_positions_by_key = {}
    for accession in methylated_proteins:
        if accession not in sequences:
            continue
        positions = residue_positions(sequences[accession], {'R'})
        if intervals_by_acc is None:
            if len(positions):
                arg_positions_by_acc[accession] = positions
            continue
        flags = position_is_in_idr(positions, intervals_by_acc.get(accession, pd.DataFrame(columns=['fragment_start', 'fragment_end'])))
        if args.restrict_context == 'IDR':
            subset = positions[flags]
            if len(subset):
                arg_positions_by_acc[accession] = subset
                arg_positions_by_key[(accession, 'IDR')] = subset
        elif args.restrict_context == 'non_IDR':
            subset = positions[~flags]
            if len(subset):
                arg_positions_by_acc[accession] = subset
                arg_positions_by_key[(accession, 'non_IDR')] = subset
        else:
            if len(positions):
                arg_positions_by_acc[accession] = positions
            idr = positions[flags]
            non_idr = positions[~flags]
            if len(idr):
                arg_positions_by_key[(accession, 'IDR')] = idr
            if len(non_idr):
                arg_positions_by_key[(accession, 'non_IDR')] = non_idr
    methyl_counts_by_acc = sites.groupby('canonical_UniProtAC').size().to_dict()
    methyl_counts_by_key = {
        (str(accession), str(label)): int(count)
        for (accession, label), count in sites.groupby(['canonical_UniProtAC', 'idr_binary_class']).size().to_dict().items()
    } if intervals_by_acc is not None else {}
    observed_positions_by_acc = {
        accession: np.sort(group['corrected_position'].to_numpy(dtype=int))
        for accession, group in sites.groupby('canonical_UniProtAC')
    }

    observed_nearest_hist, observed_ordered_pair_hist, per_site = histograms_from_positions(
        observed_positions_by_acc,
        max_distance=args.max_distance,
        position_col='corrected_position',
        neighbor_distance_col='nearest_methyl_arg_distance',
        neighbor_flag_col='has_neighbor_within_max_distance',
    )
    observed_curve = summarize_curve_from_histograms(
        nearest_hist=observed_nearest_hist,
        ordered_pair_hist=observed_ordered_pair_hist,
        total_sites=len(sites),
        max_distance=args.max_distance,
    )

    if args.null_match_disorder and intervals_by_acc is not None:
        permutations = permutation_summary_grouped(
            candidate_positions_by_key=arg_positions_by_key,
            site_counts_by_key=methyl_counts_by_key,
            max_distance=args.max_distance,
            permutations=args.permutations,
            seed=args.seed,
            progress_every=args.progress_every,
        )
    else:
        permutations = permutation_summary(
            candidate_positions_by_acc=arg_positions_by_acc,
            site_counts_by_acc=methyl_counts_by_acc,
            max_distance=args.max_distance,
            permutations=args.permutations,
            seed=args.seed,
            progress_every=args.progress_every,
        )
    null_curve = aggregate_null_curve(permutations)
    curve = observed_curve.merge(null_curve, how='left', on='distance_aa')
    curve = add_empirical_p_values(curve, permutation_df=permutations)
    curve['delta_prob_neighbor_within_x_vs_null'] = curve['prob_neighbor_within_x'] - curve['null_mean_prob_neighbor_within_x']
    curve['ratio_prob_neighbor_within_x_vs_null'] = curve['prob_neighbor_within_x'] / curve['null_mean_prob_neighbor_within_x'].replace(0, np.nan)
    curve['delta_other_methyl_args_within_x_vs_null'] = curve['mean_other_methyl_args_within_x'] - curve['null_mean_other_methyl_args_within_x']
    curve['ratio_other_methyl_args_within_x_vs_null'] = curve['mean_other_methyl_args_within_x'] / curve['null_mean_other_methyl_args_within_x'].replace(0, np.nan)
    curve['restrict_context'] = args.restrict_context
    curve['null_match_disorder'] = bool(args.null_match_disorder and intervals_by_acc is not None)

    protein_summary = protein_density_summary(observed_positions_by_acc, windows=args.checkpoints)
    rename_map = {'site_count': 'methyl_site_count'}
    for window in args.checkpoints:
        rename_map[f'max_sites_in_{window}aa_window'] = f'max_methyl_sites_in_{window}aa_window'
    protein_summary = protein_summary.rename(columns=rename_map)
    if 'substrate_genename' in sites.columns:
        gene_map = (
            sites[['canonical_UniProtAC', 'substrate_genename']]
            .drop_duplicates('canonical_UniProtAC')
            .rename(columns={'substrate_genename': 'gene'})
        )
        protein_summary = protein_summary.merge(gene_map, how='left', on='canonical_UniProtAC')

    per_site = per_site.merge(
        sites[['canonical_UniProtAC', 'corrected_position', 'substrate_genename', 'source_family_count_exact', 'source_family_count_fuzzy']].drop_duplicates(
            ['canonical_UniProtAC', 'corrected_position']
        ),
        how='left',
        on=['canonical_UniProtAC', 'corrected_position'],
    )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(per_site, outdir / 'per_site_nearest_methyl_arg_distance.tsv')
    save_table(curve, outdir / 'nearest_neighbor_curve_vs_null.tsv')
    checkpoint_df = checkpoint_summary(curve, args.checkpoints)
    concise_df = checkpoint_summary(curve, args.summary_checkpoints).copy()
    concise_df['log2_ratio_prob_neighbor_within_x_vs_null'] = np.log2(concise_df['ratio_prob_neighbor_within_x_vs_null'])
    concise_df['log2_ratio_other_methyl_args_within_x_vs_null'] = np.log2(concise_df['ratio_other_methyl_args_within_x_vs_null'])
    distance_bin_df = distance_bin_summary(curve, permutations, args.summary_checkpoints)
    rolling_df = rolling_probability_ratio(curve, permutations, window_size=5, max_distance=min(40, args.max_distance))
    save_table(checkpoint_df, outdir / 'nearest_neighbor_curve_checkpoints.tsv')
    save_table(concise_df, outdir / 'nearest_neighbor_summary_cumulative.tsv')
    save_table(concise_df, outdir / 'nearest_neighbor_summary_5_10_20.tsv')
    save_table(distance_bin_df, outdir / 'nearest_neighbor_distance_bin_summary.tsv')
    save_table(rolling_df, outdir / 'nearest_neighbor_rolling_enrichment.tsv')
    save_table(permutations, outdir / 'nearest_neighbor_curve_null_permutations.tsv')
    save_table(protein_summary, outdir / 'protein_local_density_summary.tsv')
    save_table(
        pd.DataFrame([{
            'integrated_input_sites': input_site_count,
            'sites_analyzed_after_sequence_validation': len(sites),
            'excluded_sites_without_usable_sequence': input_site_count - len(sites),
            'proteins_analyzed': len(observed_positions_by_acc),
            'permutations': args.permutations,
            'max_distance': args.max_distance,
            'restrict_context': args.restrict_context,
            'null_match_disorder': bool(args.null_match_disorder and intervals_by_acc is not None),
        }]),
        outdir / 'clustering_qc.tsv',
    )

    fig, ax = plt.subplots(figsize=(9.0, 5.5))
    ax.plot(curve['distance_aa'], curve['prob_neighbor_within_x'], label='Observed', linewidth=2.2, color=BREWER_COLORS['orange'])
    ax.plot(curve['distance_aa'], curve['null_mean_prob_neighbor_within_x'], label='Null mean', linewidth=1.8, linestyle='--', color=BREWER_COLORS['blue'])
    ax.fill_between(
        curve['distance_aa'],
        curve['null_q025_prob_neighbor_within_x'],
        curve['null_q975_prob_neighbor_within_x'],
        alpha=0.18,
        color=BREWER_COLORS['blue'],
        label='Null 95% interval',
    )
    style_axis(ax, grid_axis='both')
    ax.set_xlabel('Distance X from methyl-Arg (aa)')
    ax.set_ylabel('P(nearest methyl-Arg <= X)')
    ax.set_title('Nearest methylarginine neighbor probability curve')
    ax.legend()
    fig.tight_layout()
    save_figure(fig, outdir / 'nearest_neighbor_probability_curve')

    curve40 = curve[curve['distance_aa'] <= 40].copy()
    fig1b, ax1b = plt.subplots(figsize=(8.8, 4.8))
    ax1b.plot(
        curve40['distance_aa'],
        np.log2(curve40['ratio_prob_neighbor_within_x_vs_null'].replace(0, np.nan)),
        linewidth=2.2,
        color=BREWER_COLORS['orange'],
    )
    sig40 = curve40[curve40['empirical_p_greater_prob_neighbor_within_x'] <= 0.05]
    if not sig40.empty:
        ax1b.scatter(
            sig40['distance_aa'],
            np.log2(sig40['ratio_prob_neighbor_within_x_vs_null'].replace(0, np.nan)),
            color=BREWER_COLORS['purple'],
            s=18,
            zorder=3,
        )
    style_axis(ax1b, zero='y', grid_axis='both')
    ax1b.set_xlabel('Distance from methyl-Arg (<= X aa)')
    ax1b.set_ylabel('log2(observed/null nearest-neighbor probability)')
    ax1b.set_title('Cumulative methyl-Arg neighbor enrichment')
    ax1b.set_xlim(1, 40)
    for checkpoint in [5, 10, 20, 40]:
        row = curve40.loc[curve40['distance_aa'] == checkpoint]
        if row.empty:
            continue
        row = row.iloc[0]
        ax1b.axvline(checkpoint, color=BREWER_COLORS['light_gray'], linestyle='--', linewidth=0.8)
        ax1b.text(checkpoint, np.log2(row['ratio_prob_neighbor_within_x_vs_null']), f' {checkpoint}', fontsize=8, va='bottom', color=BREWER_COLORS['dark_gray'])
    fig1b.text(
        0.99,
        0.01,
        'Points mark empirical p<=0.05 from the within-protein permutation null.',
        ha='right',
        va='bottom',
        fontsize=8,
        color=BREWER_COLORS['dark_gray'],
    )
    fig1b.tight_layout()
    save_figure(fig1b, outdir / 'arg_methyl_neighbor_cumulative_enrichment')

    fig1c, ax1c = plt.subplots(figsize=(8.8, 4.8))
    ax1c.plot(rolling_df['window_center_aa'], rolling_df['log2_ratio_vs_null'], linewidth=2.2, color=BREWER_COLORS['green'])
    sig_roll = rolling_df[rolling_df['empirical_q_greater'] <= 0.05]
    if not sig_roll.empty:
        ax1c.scatter(sig_roll['window_center_aa'], sig_roll['log2_ratio_vs_null'], color=BREWER_COLORS['purple'], s=18, zorder=3)
    style_axis(ax1c, zero='y', grid_axis='both')
    ax1c.set_xlabel('Nearest methyl-Arg distance (rolling 5-aa window)')
    ax1c.set_ylabel('log2(observed/null nearest-neighbor probability)')
    ax1c.set_title('Rolling methyl-Arg neighbor enrichment')
    ax1c.set_xlim(3, 38)
    for checkpoint in [5, 10, 20, 40]:
        ax1c.axvline(checkpoint, color=BREWER_COLORS['light_gray'], linestyle='--', linewidth=0.8)
    fig1c.text(
        0.99,
        0.01,
        'Each point is a 5-aa rolling window; points mark empirical q<=0.05.',
        ha='right',
        va='bottom',
        fontsize=8,
        color=BREWER_COLORS['dark_gray'],
    )
    fig1c.tight_layout()
    save_figure(fig1c, outdir / 'arg_methyl_neighbor_rolling_enrichment')

    fig2, ax2 = plt.subplots(figsize=(9.0, 5.5))
    ax2.plot(curve['distance_aa'], curve['mean_other_methyl_args_within_x'], label='Observed', linewidth=2.2, color=BREWER_COLORS['green'])
    ax2.plot(curve['distance_aa'], curve['null_mean_other_methyl_args_within_x'], label='Null mean', linewidth=1.8, linestyle='--', color=BREWER_COLORS['purple'])
    ax2.fill_between(
        curve['distance_aa'],
        curve['null_q025_other_methyl_args_within_x'],
        curve['null_q975_other_methyl_args_within_x'],
        alpha=0.18,
        color=BREWER_COLORS['purple'],
        label='Null 95% interval',
    )
    style_axis(ax2, grid_axis='both')
    ax2.set_xlabel('Distance X from methyl-Arg (aa)')
    ax2.set_ylabel('Mean other methyl-Args within X aa')
    ax2.set_title('Local methylarginine density curve')
    ax2.legend()
    fig2.tight_layout()
    save_figure(fig2, outdir / 'local_methyl_arg_density_curve')

    top_dense = protein_summary.head(20).copy()
    if not top_dense.empty:
        label_col = 'gene' if 'gene' in top_dense.columns else 'canonical_UniProtAC'
        labels = top_dense[label_col].fillna(top_dense['canonical_UniProtAC']).astype(str)
        window = max(args.checkpoints)
        density_col = f'max_methyl_sites_in_{window}aa_window'
        fig3, ax3 = plt.subplots(figsize=(10.0, 7.0))
        top_dense = top_dense.sort_values([density_col, 'methyl_site_count'], ascending=[True, True])
        ax3.barh(labels.loc[top_dense.index], top_dense[density_col], color=BREWER_COLORS['teal'], edgecolor='white', linewidth=0.8)
        style_axis(ax3, grid_axis='x')
        ax3.set_xlabel(f'Max methyl sites in {window} aa window')
        ax3.set_ylabel('Protein')
        ax3.set_title('Top proteins by local methylarginine density')
        fig3.tight_layout()
        save_figure(fig3, outdir / 'top_proteins_by_local_methyl_arg_density')

    if not distance_bin_df.empty:
        fig4, ax4 = plt.subplots(figsize=(6.8, 4.8))
        ax4.bar(
            distance_bin_df['distance_bin'],
            distance_bin_df['log2_ratio_vs_null'],
            color=BREWER_COLORS['orange'],
            edgecolor='white',
            linewidth=0.8,
        )
        style_axis(ax4, zero='y', grid_axis='y')
        ax4.set_xlabel('Nearest methyl-Arg distance bin')
        ax4.set_ylabel('log2(observed/null nearest-neighbor probability)')
        ax4.set_title('Nearest-neighbor probability enrichment summary')
        for idx, row in distance_bin_df.reset_index(drop=True).iterrows():
            p_text = f"p={format_p_value(row['empirical_p_greater'])}" if float(row['empirical_p_greater']) <= 0.05 else 'n.s.'
            ax4.text(
                idx,
                float(row['log2_ratio_vs_null']),
                f"obs={row['observed_probability']:.3f}\nnull={row['null_mean_probability']:.3f}\n{p_text}",
                ha='center',
                va='bottom' if float(row['log2_ratio_vs_null']) >= 0 else 'top',
                fontsize=7.4,
                color=BREWER_COLORS['dark_gray'],
            )
        fig4.text(
            0.99,
            0.01,
            'Not an odds ratio: each bar is the observed nearest-neighbor probability in the bin divided by the within-protein permutation null mean.',
            ha='right',
            va='bottom',
            fontsize=8,
            color=BREWER_COLORS['dark_gray'],
        )
        fig4.tight_layout()
        save_figure(fig4, outdir / 'arg_methyl_neighbor_clustering_summary')

    print({
        'sites_analyzed': len(sites),
        'proteins_analyzed': len(observed_positions_by_acc),
        'permutations': args.permutations,
        'max_distance': args.max_distance,
        'output': str(outdir),
    })


if __name__ == '__main__':
    main()

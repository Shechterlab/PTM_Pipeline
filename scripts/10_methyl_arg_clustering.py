from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import canonical_site_table, parse_fasta, save_figure, save_table
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--integrated-sites', required=True)
    ap.add_argument('--canonical-fasta', required=True)
    ap.add_argument('--outdir', default='results/methyl_arg_clustering')
    ap.add_argument('--max-distance', type=int, default=100)
    ap.add_argument('--permutations', type=int, default=250)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--progress-every', type=int, default=25)
    ap.add_argument('--checkpoints', nargs='*', type=int, default=[3, 5, 10, 20, 50])
    ap.add_argument('--disorder-intervals', default='')
    ap.add_argument('--null-match-disorder', action='store_true')
    ap.add_argument('--restrict-context', choices=['all', 'IDR', 'non_IDR'], default='all')
    args = ap.parse_args()

    sites = pd.read_csv(args.integrated_sites, sep='\t', low_memory=False)
    sequences = parse_fasta(args.canonical_fasta)
    sites = validate_sites(sites, sequences=sequences)

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
    save_table(checkpoint_summary(curve, args.checkpoints), outdir / 'nearest_neighbor_curve_checkpoints.tsv')
    save_table(permutations, outdir / 'nearest_neighbor_curve_null_permutations.tsv')
    save_table(protein_summary, outdir / 'protein_local_density_summary.tsv')

    fig, ax = plt.subplots(figsize=(9.0, 5.5))
    ax.plot(curve['distance_aa'], curve['prob_neighbor_within_x'], label='Observed', linewidth=2.2)
    ax.plot(curve['distance_aa'], curve['null_mean_prob_neighbor_within_x'], label='Null mean', linewidth=1.8, linestyle='--')
    ax.fill_between(
        curve['distance_aa'],
        curve['null_q025_prob_neighbor_within_x'],
        curve['null_q975_prob_neighbor_within_x'],
        alpha=0.2,
        label='Null 95% interval',
    )
    ax.set_xlabel('Distance X from methyl-Arg (aa)')
    ax.set_ylabel('P(nearest methyl-Arg <= X)')
    ax.set_title('Nearest methylarginine neighbor probability curve')
    ax.legend()
    fig.tight_layout()
    save_figure(fig, outdir / 'nearest_neighbor_probability_curve')

    fig2, ax2 = plt.subplots(figsize=(9.0, 5.5))
    ax2.plot(curve['distance_aa'], curve['mean_other_methyl_args_within_x'], label='Observed', linewidth=2.2)
    ax2.plot(curve['distance_aa'], curve['null_mean_other_methyl_args_within_x'], label='Null mean', linewidth=1.8, linestyle='--')
    ax2.fill_between(
        curve['distance_aa'],
        curve['null_q025_other_methyl_args_within_x'],
        curve['null_q975_other_methyl_args_within_x'],
        alpha=0.2,
        label='Null 95% interval',
    )
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
        ax3.barh(labels.loc[top_dense.index], top_dense[density_col])
        ax3.set_xlabel(f'Max methyl sites in {window} aa window')
        ax3.set_ylabel('Protein')
        ax3.set_title('Top proteins by local methylarginine density')
        fig3.tight_layout()
        save_figure(fig3, outdir / 'top_proteins_by_local_methyl_arg_density')

    print({
        'sites_analyzed': len(sites),
        'proteins_analyzed': len(observed_positions_by_acc),
        'permutations': args.permutations,
        'max_distance': args.max_distance,
        'output': str(outdir),
    })


if __name__ == '__main__':
    main()

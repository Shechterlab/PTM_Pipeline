from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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
    residue_positions,
    summarize_curve_from_histograms,
)
from common import canonicalize_uniprot_accession, parse_fasta, save_figure, save_table


DEFAULT_PTM_GROUPS = [
    'Arg methylation',
    'Phosphorylation',
    'Lys acylation',
    'Lys methylation',
    'Ubiquitin/SUMO',
]

PTM_GROUP_TARGET_RESIDUES = {
    'Arg methylation': {'R'},
    'Phosphorylation': {'S', 'T', 'Y'},
    'Lys acylation': {'K'},
    'Lys methylation': {'K'},
    'Ubiquitin/SUMO': {'K'},
}


def normalize_master_sites(master: pd.DataFrame, sequences: dict[str, str], ptm_groups: list[str]) -> pd.DataFrame:
    out = master[master['ptm_group'].isin(ptm_groups)].copy()
    out['canonical_UniProtAC'] = out['substrate_UniProtAC'].astype(str).map(canonicalize_uniprot_accession)
    out['position'] = pd.to_numeric(out['position'], errors='coerce')
    out = out[out['position'].notna()].copy()
    out['position'] = out['position'].astype(int)
    out['residue'] = out['residue'].astype(str).str.upper().str[:1]
    out = out[
        out.apply(
            lambda row: row['residue'] in PTM_GROUP_TARGET_RESIDUES.get(str(row['ptm_group']), {str(row['residue'])}),
            axis=1,
        )
    ].copy()

    def valid_row(row: pd.Series) -> bool:
        seq = sequences.get(str(row['canonical_UniProtAC']), '')
        pos = int(row['position'])
        return bool(seq) and 1 <= pos <= len(seq) and seq[pos - 1] == row['residue']

    out = out[out.apply(valid_row, axis=1)].copy()
    out['site'] = out['residue'] + out['position'].astype(str)
    out['site_key'] = out['canonical_UniProtAC'] + ':' + out['site']
    return out.drop_duplicates(['ptm_group', 'canonical_UniProtAC', 'residue', 'position']).copy()


def normalize_arg_methyl_sites(df: pd.DataFrame, sequences: dict[str, str]) -> pd.DataFrame:
    out = df.copy()
    out['canonical_UniProtAC'] = out['canonical_UniProtAC'].astype(str).map(canonicalize_uniprot_accession)
    out['position'] = pd.to_numeric(out['corrected_position'], errors='coerce')
    out = out[out['position'].notna()].copy()
    out['position'] = out['position'].astype(int)
    out['residue'] = 'R'
    out['ptm_group'] = 'Arg methylation'

    def valid_row(row: pd.Series) -> bool:
        seq = sequences.get(str(row['canonical_UniProtAC']), '')
        pos = int(row['position'])
        return bool(seq) and 1 <= pos <= len(seq) and seq[pos - 1] == 'R'

    out = out[out.apply(valid_row, axis=1)].copy()
    out['site'] = out['residue'] + out['position'].astype(str)
    out['site_key'] = out['canonical_UniProtAC'] + ':' + out['site']
    return out.drop_duplicates(['ptm_group', 'canonical_UniProtAC', 'residue', 'position']).copy()


def plot_curve_family(curves: pd.DataFrame, value_col: str, ylabel: str, title: str, outpath: Path) -> None:
    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    for ptm_group, group in curves.groupby('ptm_group'):
        group = group.sort_values('distance_aa')
        ax.plot(group['distance_aa'], group[value_col], label=ptm_group, linewidth=2)
    ax.set_xlabel('Distance X from modified residue (aa)')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    save_figure(fig, outpath)


def annotate_idr_binary(df: pd.DataFrame, intervals_by_acc: dict[str, pd.DataFrame]) -> pd.DataFrame:
    out = df.copy()
    out['idr_binary_class'] = 'non_IDR'
    for accession, idx in out.groupby('canonical_UniProtAC').groups.items():
        positions = out.loc[idx, 'position'].to_numpy(dtype=int)
        flags = position_is_in_idr(positions, intervals_by_acc.get(accession, pd.DataFrame(columns=['fragment_start', 'fragment_end'])))
        out.loc[idx, 'idr_binary_class'] = np.where(flags, 'IDR', 'non_IDR')
    return out


def build_candidate_positions(
    proteins: list[str],
    target_residues: set[str],
    sequences: dict[str, str],
    intervals_by_acc: dict[str, pd.DataFrame] | None,
    restrict_context: str,
) -> tuple[dict[str, np.ndarray], dict[tuple[str, str], np.ndarray]]:
    flat = {}
    grouped = {}
    for accession in proteins:
        if accession not in sequences:
            continue
        positions = residue_positions(sequences[accession], target_residues)
        if len(positions) == 0:
            continue
        if intervals_by_acc is None:
            flat[accession] = positions
            grouped[(accession, 'all')] = positions
            continue
        flags = position_is_in_idr(positions, intervals_by_acc.get(accession, pd.DataFrame(columns=['fragment_start', 'fragment_end'])))
        if restrict_context == 'IDR':
            subset = positions[flags]
            if len(subset):
                flat[accession] = subset
                grouped[(accession, 'IDR')] = subset
        elif restrict_context == 'non_IDR':
            subset = positions[~flags]
            if len(subset):
                flat[accession] = subset
                grouped[(accession, 'non_IDR')] = subset
        else:
            flat[accession] = positions
            idr = positions[flags]
            non_idr = positions[~flags]
            if len(idr):
                grouped[(accession, 'IDR')] = idr
            if len(non_idr):
                grouped[(accession, 'non_IDR')] = non_idr
    return flat, grouped


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--base-master', required=True)
    ap.add_argument('--integrated-arg-sites', required=True)
    ap.add_argument('--canonical-fasta', required=True)
    ap.add_argument('--outdir', default='results/ptm_clustering')
    ap.add_argument('--ptm-groups', nargs='*', default=list(DEFAULT_PTM_GROUPS))
    ap.add_argument('--max-distance', type=int, default=100)
    ap.add_argument('--permutations', type=int, default=100000)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--progress-every', type=int, default=5000)
    ap.add_argument('--checkpoints', nargs='*', type=int, default=[3, 5, 10, 20, 50])
    ap.add_argument('--disorder-intervals', default='')
    ap.add_argument('--null-match-disorder', action='store_true')
    ap.add_argument('--restrict-context', choices=['all', 'IDR', 'non_IDR'], default='all')
    args = ap.parse_args()

    sequences = parse_fasta(args.canonical_fasta)
    master = pd.read_csv(args.base_master, sep='\t', low_memory=False)
    arg_sites = pd.read_csv(args.integrated_arg_sites, sep='\t', low_memory=False)

    ptm_sites = normalize_master_sites(master, sequences=sequences, ptm_groups=[group for group in args.ptm_groups if group != 'Arg methylation'])
    arg_sites = normalize_arg_methyl_sites(arg_sites, sequences=sequences)
    ptm_sites = pd.concat([ptm_sites, arg_sites], ignore_index=True, sort=False)
    ptm_sites = ptm_sites[ptm_sites['ptm_group'].isin(args.ptm_groups)].copy()

    intervals_by_acc = None
    if args.disorder_intervals:
        disorder = normalize_disorder_intervals(pd.read_csv(args.disorder_intervals, sep='\t', low_memory=False))
        intervals_by_acc = merged_intervals_by_accession(disorder)
        ptm_sites = annotate_idr_binary(ptm_sites, intervals_by_acc)
        if args.restrict_context != 'all':
            ptm_sites = ptm_sites[ptm_sites['idr_binary_class'] == args.restrict_context].copy()
    elif args.null_match_disorder or args.restrict_context != 'all':
        raise ValueError('--disorder-intervals is required for --null-match-disorder or --restrict-context')

    curve_rows = []
    permutation_rows = []
    checkpoint_rows = []
    summary_rows = []

    for ptm_group, group in ptm_sites.groupby('ptm_group'):
        target_residues = PTM_GROUP_TARGET_RESIDUES.get(ptm_group)
        if not target_residues:
            continue

        proteins = sorted(set(group['canonical_UniProtAC'].astype(str)))
        candidate_positions_by_acc, candidate_positions_by_key = build_candidate_positions(
            proteins=proteins,
            target_residues=target_residues,
            sequences=sequences,
            intervals_by_acc=intervals_by_acc,
            restrict_context=args.restrict_context,
        )
        observed_positions_by_acc = {
            accession: np.sort(acc_group['position'].to_numpy(dtype=int))
            for accession, acc_group in group.groupby('canonical_UniProtAC')
        }
        site_counts_by_acc = group.groupby('canonical_UniProtAC').size().to_dict()
        site_counts_by_key = {
            (str(accession), str(label)): int(count)
            for (accession, label), count in group.groupby(['canonical_UniProtAC', 'idr_binary_class']).size().to_dict().items()
        } if intervals_by_acc is not None else {}

        nearest_hist, ordered_pair_hist, _ = histograms_from_positions(observed_positions_by_acc, max_distance=args.max_distance)
        observed_curve = summarize_curve_from_histograms(
            nearest_hist=nearest_hist,
            ordered_pair_hist=ordered_pair_hist,
            total_sites=len(group),
            max_distance=args.max_distance,
        )
        if args.null_match_disorder and intervals_by_acc is not None:
            null_permutations = permutation_summary_grouped(
                candidate_positions_by_key=candidate_positions_by_key,
                site_counts_by_key=site_counts_by_key,
                max_distance=args.max_distance,
                permutations=args.permutations,
                seed=args.seed,
                progress_every=args.progress_every,
            )
        else:
            null_permutations = permutation_summary(
                candidate_positions_by_acc=candidate_positions_by_acc,
                site_counts_by_acc=site_counts_by_acc,
                max_distance=args.max_distance,
                permutations=args.permutations,
                seed=args.seed,
                progress_every=args.progress_every,
            )
        null_curve = aggregate_null_curve(null_permutations)
        curve = observed_curve.merge(null_curve, how='left', on='distance_aa')
        curve = add_empirical_p_values(curve, permutation_df=null_permutations)
        curve['delta_prob_neighbor_within_x_vs_null'] = curve['prob_neighbor_within_x'] - curve['null_mean_prob_neighbor_within_x']
        curve['ratio_prob_neighbor_within_x_vs_null'] = curve['prob_neighbor_within_x'] / curve['null_mean_prob_neighbor_within_x'].replace(0, np.nan)
        curve['delta_other_methyl_args_within_x_vs_null'] = curve['mean_other_methyl_args_within_x'] - curve['null_mean_other_methyl_args_within_x']
        curve['ratio_other_methyl_args_within_x_vs_null'] = curve['mean_other_methyl_args_within_x'] / curve['null_mean_other_methyl_args_within_x'].replace(0, np.nan)
        curve['mean_other_modified_sites_within_x'] = curve['mean_other_methyl_args_within_x']
        curve['null_mean_other_modified_sites_within_x'] = curve['null_mean_other_methyl_args_within_x']
        curve['delta_other_modified_sites_within_x_vs_null'] = curve['delta_other_methyl_args_within_x_vs_null']
        curve['ratio_other_modified_sites_within_x_vs_null'] = curve['ratio_other_methyl_args_within_x_vs_null']
        curve['ptm_group'] = ptm_group
        curve_rows.append(curve)

        null_permutations['ptm_group'] = ptm_group
        permutation_rows.append(null_permutations)

        checkpoints = checkpoint_summary(curve, args.checkpoints)
        checkpoints['ptm_group'] = ptm_group
        checkpoint_rows.append(checkpoints)

        summary_rows.append({
            'ptm_group': ptm_group,
            'site_rows': len(group),
            'proteins_with_sites': len(observed_positions_by_acc),
            'candidate_proteins': len(candidate_positions_by_acc),
            'candidate_target_residues': int(sum(len(values) for values in candidate_positions_by_acc.values())),
            'restrict_context': args.restrict_context,
            'null_match_disorder': bool(args.null_match_disorder and intervals_by_acc is not None),
        })

    curves = pd.concat(curve_rows, ignore_index=True)
    checkpoints = pd.concat(checkpoint_rows, ignore_index=True)
    permutations = pd.concat(permutation_rows, ignore_index=True)
    summary = pd.DataFrame(summary_rows).sort_values('site_rows', ascending=False)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(summary, outdir / 'ptm_clustering_summary.tsv')
    save_table(curves, outdir / 'ptm_clustering_curves.tsv')
    save_table(checkpoints, outdir / 'ptm_clustering_checkpoints.tsv')
    save_table(permutations, outdir / 'ptm_clustering_null_permutations.tsv')

    plot_curve_family(
        curves=curves,
        value_col='prob_neighbor_within_x',
        ylabel='P(nearest modified residue <= X)',
        title='Cross-PTM nearest-neighbor probability curves',
        outpath=outdir / 'cross_ptm_nearest_neighbor_probability_curves',
    )
    plot_curve_family(
        curves=curves,
        value_col='ratio_prob_neighbor_within_x_vs_null',
        ylabel='Observed / null P(nearest modified residue <= X)',
        title='Cross-PTM nearest-neighbor enrichment vs null',
        outpath=outdir / 'cross_ptm_nearest_neighbor_enrichment_vs_null',
    )
    plot_curve_family(
        curves=curves,
        value_col='mean_other_modified_sites_within_x',
        ylabel='Mean other modified residues within X aa',
        title='Cross-PTM local modification density curves',
        outpath=outdir / 'cross_ptm_local_density_curves',
    )
    plot_curve_family(
        curves=curves,
        value_col='ratio_other_modified_sites_within_x_vs_null',
        ylabel='Observed / null mean modified residues within X aa',
        title='Cross-PTM local density enrichment vs null',
        outpath=outdir / 'cross_ptm_local_density_enrichment_vs_null',
    )

    print({
        'ptm_groups': sorted(set(curves['ptm_group'].astype(str))),
        'rows': len(curves),
        'output': str(outdir),
    })


if __name__ == '__main__':
    main()

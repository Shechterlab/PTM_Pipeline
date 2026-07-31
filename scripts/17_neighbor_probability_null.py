from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from clustering_utils import (
    merged_intervals_by_accession,
    nearest_neighbor_distances,
    normalize_disorder_intervals,
    position_is_in_idr,
    residue_positions,
    sample_positions,
    sample_positions_grouped,
)
from common import (
    BREWER_COLORS,
    add_q_values,
    apply_paper_style,
    canonical_site_table,
    parse_fasta,
    save_figure,
    save_table,
    style_axis,
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


def nearest_hist_only(positions_by_acc: dict[str, np.ndarray], max_distance: int) -> tuple[np.ndarray, int]:
    hist = np.zeros(max_distance + 1, dtype=np.int64)
    total_sites = 0
    for positions in positions_by_acc.values():
        positions = np.sort(np.asarray(positions, dtype=int))
        total_sites += len(positions)
        nearest = nearest_neighbor_distances(positions)
        if len(nearest) == 0:
            continue
        finite = nearest[np.isfinite(nearest)]
        if len(finite) == 0:
            continue
        finite = finite[(finite >= 1) & (finite <= max_distance)].astype(int)
        if len(finite):
            hist += np.bincount(finite, minlength=max_distance + 1)[: max_distance + 1]
    return hist, total_sites


def sample_candidate_positions(
    sites: pd.DataFrame,
    sequences: dict[str, str],
    intervals_by_acc: dict[str, pd.DataFrame] | None,
    restrict_context: str,
) -> tuple[dict[str, np.ndarray], dict[tuple[str, str], np.ndarray], dict[str, int], dict[tuple[str, str], int], dict[str, np.ndarray]]:
    methylated_proteins = sorted(set(sites['canonical_UniProtAC'].astype(str)))
    arg_positions_by_acc: dict[str, np.ndarray] = {}
    arg_positions_by_key: dict[tuple[str, str], np.ndarray] = {}
    for accession in methylated_proteins:
        if accession not in sequences:
            continue
        positions = residue_positions(sequences[accession], {'R'})
        if intervals_by_acc is None:
            if len(positions):
                arg_positions_by_acc[accession] = positions
            continue
        flags = position_is_in_idr(positions, intervals_by_acc.get(accession, pd.DataFrame(columns=['fragment_start', 'fragment_end'])))
        if restrict_context == 'IDR':
            subset = positions[flags]
            if len(subset):
                arg_positions_by_acc[accession] = subset
                arg_positions_by_key[(accession, 'IDR')] = subset
        elif restrict_context == 'non_IDR':
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

    site_counts_by_acc = sites.groupby('canonical_UniProtAC').size().to_dict()
    site_counts_by_key = {
        (str(accession), str(label)): int(count)
        for (accession, label), count in sites.groupby(['canonical_UniProtAC', 'idr_binary_class']).size().to_dict().items()
    } if intervals_by_acc is not None else {}
    observed_positions_by_acc = {
        accession: np.sort(group['corrected_position'].to_numpy(dtype=int))
        for accession, group in sites.groupby('canonical_UniProtAC')
    }
    return arg_positions_by_acc, arg_positions_by_key, site_counts_by_acc, site_counts_by_key, observed_positions_by_acc


def plot_probability_curve(stats: pd.DataFrame, outdir: Path) -> None:
    apply_paper_style()
    fig, ax = plt.subplots(figsize=(3.8, 4.0))
    ax.fill_between(
        stats['distance_aa'],
        100.0 * stats['null_q025_probability'],
        100.0 * stats['null_q975_probability'],
        alpha=0.14,
        color=BREWER_COLORS['light_gray'],
    )
    ax.plot(stats['distance_aa'], 100.0 * stats['null_mean_probability'], linewidth=2.0, linestyle='--', color=BREWER_COLORS['mid_gray'])
    ax.plot(stats['distance_aa'], 100.0 * stats['observed_probability'], linewidth=2.4, color=BREWER_COLORS['blue'])
    style_axis(ax, grid_axis='y')
    ax.set_xlim(1, int(stats['distance_aa'].max()))
    ax.set_xticks([10, 20, 30, 40])
    ax.set_xlabel('Within X aa')
    ax.set_ylabel('% of methyl-Arg sites with another nearby')
    ax.set_title('Nearest methyl-Arg neighbor probability')

    end_x = float(stats['distance_aa'].max())
    end_obs = 100.0 * float(stats.loc[stats['distance_aa'] == int(end_x), 'observed_probability'].iloc[0])
    end_null = 100.0 * float(stats.loc[stats['distance_aa'] == int(end_x), 'null_mean_probability'].iloc[0])
    ax.text(end_x + 0.5, end_obs, 'Observed', color=BREWER_COLORS['blue'], fontsize=9.0, va='center', ha='left')
    ax.text(end_x + 0.5, end_null, 'Within-protein null', color=BREWER_COLORS['dark_gray'], fontsize=9.0, va='center', ha='left')

    ax.set_ylim(0, max(62.0, 100.0 * float(stats['observed_probability'].max()) + 4.0))
    fig.tight_layout()
    save_figure(fig, outdir / 'nearest_neighbor_probability_curve')


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--integrated-sites', required=True)
    ap.add_argument('--canonical-fasta', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--disorder-intervals', default='')
    ap.add_argument('--null-match-disorder', action='store_true')
    ap.add_argument('--restrict-context', choices=['all', 'IDR', 'non_IDR'], default='all')
    ap.add_argument('--max-distance', type=int, default=40)
    ap.add_argument('--permutations', type=int, default=100000)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--progress-every', type=int, default=5000)
    args = ap.parse_args()

    sites = pd.read_csv(args.integrated_sites, sep='\t', low_memory=False)
    input_site_count = len(sites)
    sequences = parse_fasta(args.canonical_fasta)
    sites = validate_sites(sites, sequences)

    intervals_by_acc = None
    if args.disorder_intervals:
        disorder = normalize_disorder_intervals(pd.read_csv(args.disorder_intervals, sep='\t', low_memory=False))
        intervals_by_acc = merged_intervals_by_accession(disorder)
        sites = annotate_idr_binary(sites, intervals_by_acc)
        if args.restrict_context != 'all':
            sites = sites[sites['idr_binary_class'] == args.restrict_context].copy()
    elif args.null_match_disorder or args.restrict_context != 'all':
        raise ValueError('--disorder-intervals is required for --null-match-disorder or --restrict-context')

    arg_positions_by_acc, arg_positions_by_key, site_counts_by_acc, site_counts_by_key, observed_positions_by_acc = sample_candidate_positions(
        sites=sites,
        sequences=sequences,
        intervals_by_acc=intervals_by_acc,
        restrict_context=args.restrict_context,
    )

    observed_hist, observed_total = nearest_hist_only(observed_positions_by_acc, max_distance=args.max_distance)
    observed_cumulative = np.cumsum(observed_hist)[1:] / observed_total

    rng = np.random.default_rng(args.seed)
    null_probabilities = np.empty((args.permutations, args.max_distance), dtype=np.float32)
    for iteration in range(args.permutations):
        if args.null_match_disorder and intervals_by_acc is not None:
            sampled = sample_positions_grouped(arg_positions_by_key, site_counts_by_key, rng)
        else:
            sampled = sample_positions(arg_positions_by_acc, site_counts_by_acc, rng)
        hist, total_sites = nearest_hist_only(sampled, max_distance=args.max_distance)
        null_probabilities[iteration, :] = np.cumsum(hist)[1:] / total_sites
        if args.progress_every and (iteration + 1) % args.progress_every == 0:
            print(f'completed {iteration + 1} / {args.permutations} nearest-neighbor probability permutations')

    distances = np.arange(1, args.max_distance + 1, dtype=int)
    null_mean = null_probabilities.mean(axis=0)
    null_q025 = np.quantile(null_probabilities, 0.025, axis=0)
    null_q975 = np.quantile(null_probabilities, 0.975, axis=0)
    empirical_p = ((null_probabilities >= observed_cumulative).sum(axis=0) + 1) / (args.permutations + 1)

    stats = pd.DataFrame({
        'distance_aa': distances,
        'sites_total': observed_total,
        'observed_probability': observed_cumulative,
        'null_mean_probability': null_mean,
        'null_q025_probability': null_q025,
        'null_q975_probability': null_q975,
        'ratio_vs_null': observed_cumulative / null_mean,
        'delta_vs_null': observed_cumulative - null_mean,
        'empirical_p_greater': empirical_p,
    })
    stats = add_q_values(stats, p_col='empirical_p_greater', out_col='empirical_q_greater')

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(stats, outdir / 'nearest_neighbor_probability_curve_stats.tsv')
    save_table(
        stats[stats['distance_aa'].isin([5, 10, 20])].copy().assign(
            label=lambda df: df['distance_aa'].map(lambda x: f'<={int(x)} aa'),
            sites_with_neighbor_within_x=lambda df: np.rint(df['observed_probability'] * observed_total).astype(int),
        )[['distance_aa', 'label', 'sites_total', 'sites_with_neighbor_within_x', 'observed_probability', 'null_mean_probability', 'ratio_vs_null', 'empirical_p_greater', 'empirical_q_greater']],
        outdir / 'nearest_neighbor_checkpoint_summary.tsv',
    )
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
            'analysis_type': 'nearest_neighbor_probability_only',
        }]),
        outdir / 'nearest_neighbor_probability_qc.tsv',
    )
    plot_probability_curve(stats, outdir)
    print({
        'sites_analyzed': len(sites),
        'proteins_analyzed': len(observed_positions_by_acc),
        'permutations': args.permutations,
        'max_distance': args.max_distance,
        'output': str(outdir),
    })


if __name__ == '__main__':
    main()

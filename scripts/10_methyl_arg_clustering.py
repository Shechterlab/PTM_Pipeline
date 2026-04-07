from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import canonical_site_table, parse_fasta, save_figure, save_table


def residue_positions(sequence: str, residue: str = 'R') -> np.ndarray:
    seq = str(sequence).strip().upper()
    return np.array([idx for idx, aa in enumerate(seq, start=1) if aa == residue], dtype=int)


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


def nearest_neighbor_distances(positions: np.ndarray) -> np.ndarray:
    positions = np.sort(np.asarray(positions, dtype=int))
    n = len(positions)
    if n == 0:
        return np.array([], dtype=float)
    if n == 1:
        return np.array([np.inf], dtype=float)

    distances = np.full(n, np.inf, dtype=float)
    distances[0] = positions[1] - positions[0]
    distances[-1] = positions[-1] - positions[-2]
    if n > 2:
        distances[1:-1] = np.minimum(
            positions[1:-1] - positions[:-2],
            positions[2:] - positions[1:-1],
        )
    return distances


def ordered_pair_distance_histogram(positions: np.ndarray, max_distance: int) -> np.ndarray:
    positions = np.sort(np.asarray(positions, dtype=int))
    hist = np.zeros(max_distance + 1, dtype=int)
    if len(positions) <= 1 or max_distance < 1:
        return hist

    distances = np.abs(positions[:, None] - positions[None, :])
    mask = (distances > 0) & (distances <= max_distance)
    if mask.any():
        values = distances[mask].astype(int).ravel()
        hist += np.bincount(values, minlength=max_distance + 1)[: max_distance + 1]
    return hist


def histograms_from_positions(positions_by_acc: dict[str, np.ndarray], max_distance: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nearest_hist = np.zeros(max_distance + 1, dtype=int)
    ordered_pair_hist = np.zeros(max_distance + 1, dtype=int)
    per_site_rows = []

    for accession, positions in positions_by_acc.items():
        positions = np.sort(np.asarray(positions, dtype=int))
        nearest = nearest_neighbor_distances(positions)
        ordered_pair_hist += ordered_pair_distance_histogram(positions, max_distance=max_distance)
        for position, nearest_distance in zip(positions, nearest):
            finite_distance = nearest_distance if np.isfinite(nearest_distance) else pd.NA
            per_site_rows.append({
                'canonical_UniProtAC': accession,
                'corrected_position': int(position),
                'nearest_methyl_arg_distance': finite_distance,
                'has_neighbor_within_max_distance': bool(np.isfinite(nearest_distance) and nearest_distance <= max_distance),
            })
            if np.isfinite(nearest_distance) and nearest_distance <= max_distance:
                nearest_hist[int(nearest_distance)] += 1

    return nearest_hist, ordered_pair_hist, pd.DataFrame(per_site_rows)


def summarize_curve_from_histograms(
    nearest_hist: np.ndarray,
    ordered_pair_hist: np.ndarray,
    total_sites: int,
    max_distance: int,
) -> pd.DataFrame:
    distances = np.arange(1, max_distance + 1, dtype=int)
    cumulative_nearest = np.cumsum(nearest_hist)[1:]
    cumulative_pairs = np.cumsum(ordered_pair_hist)[1:]
    return pd.DataFrame({
        'distance_aa': distances,
        'sites_total': total_sites,
        'prob_neighbor_within_x': cumulative_nearest / total_sites if total_sites else np.nan,
        'mean_other_methyl_args_within_x': cumulative_pairs / total_sites if total_sites else np.nan,
        'sites_with_neighbor_within_x': cumulative_nearest,
        'ordered_neighbor_pairs_within_x': cumulative_pairs,
    })


def sample_positions(arg_positions_by_acc: dict[str, np.ndarray], methyl_counts_by_acc: dict[str, int], rng: np.random.Generator) -> dict[str, np.ndarray]:
    sampled = {}
    for accession, methyl_count in methyl_counts_by_acc.items():
        arg_positions = arg_positions_by_acc.get(accession, np.array([], dtype=int))
        if methyl_count <= 0 or len(arg_positions) == 0:
            continue
        if methyl_count >= len(arg_positions):
            sampled_positions = np.sort(arg_positions.copy())
        else:
            sampled_positions = np.sort(rng.choice(arg_positions, size=methyl_count, replace=False))
        sampled[accession] = sampled_positions.astype(int)
    return sampled


def permutation_summary(
    arg_positions_by_acc: dict[str, np.ndarray],
    methyl_counts_by_acc: dict[str, int],
    max_distance: int,
    permutations: int,
    seed: int,
    progress_every: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for iteration in range(1, permutations + 1):
        sampled = sample_positions(arg_positions_by_acc, methyl_counts_by_acc, rng)
        nearest_hist, ordered_pair_hist, _ = histograms_from_positions(sampled, max_distance=max_distance)
        total_sites = int(sum(len(values) for values in sampled.values()))
        summary = summarize_curve_from_histograms(
            nearest_hist=nearest_hist,
            ordered_pair_hist=ordered_pair_hist,
            total_sites=total_sites,
            max_distance=max_distance,
        )
        summary['permutation'] = iteration
        rows.append(summary)
        if progress_every and iteration % progress_every == 0:
            print(f'completed {iteration} / {permutations} clustering permutations')
    return pd.concat(rows, ignore_index=True)


def aggregate_null_curve(permutation_df: pd.DataFrame) -> pd.DataFrame:
    grouped = permutation_df.groupby('distance_aa', as_index=False)
    records = []
    for distance, group in grouped:
        probability = group['prob_neighbor_within_x']
        density = group['mean_other_methyl_args_within_x']
        records.append({
            'distance_aa': int(distance),
            'null_mean_prob_neighbor_within_x': float(probability.mean()),
            'null_sd_prob_neighbor_within_x': float(probability.std(ddof=1)) if len(group) > 1 else 0.0,
            'null_q025_prob_neighbor_within_x': float(probability.quantile(0.025)),
            'null_q975_prob_neighbor_within_x': float(probability.quantile(0.975)),
            'null_mean_other_methyl_args_within_x': float(density.mean()),
            'null_sd_other_methyl_args_within_x': float(density.std(ddof=1)) if len(group) > 1 else 0.0,
            'null_q025_other_methyl_args_within_x': float(density.quantile(0.025)),
            'null_q975_other_methyl_args_within_x': float(density.quantile(0.975)),
        })
    return pd.DataFrame(records)


def add_empirical_p_values(observed: pd.DataFrame, permutation_df: pd.DataFrame) -> pd.DataFrame:
    out = observed.copy()
    prob_p = []
    density_p = []
    for distance in out['distance_aa']:
        null_group = permutation_df[permutation_df['distance_aa'] == distance]
        null_prob = null_group['prob_neighbor_within_x'].to_numpy(dtype=float)
        null_density = null_group['mean_other_methyl_args_within_x'].to_numpy(dtype=float)
        obs_prob = float(out.loc[out['distance_aa'] == distance, 'prob_neighbor_within_x'].iloc[0])
        obs_density = float(out.loc[out['distance_aa'] == distance, 'mean_other_methyl_args_within_x'].iloc[0])
        prob_tail = (np.count_nonzero(null_prob >= obs_prob) + 1) / (len(null_prob) + 1)
        density_tail = (np.count_nonzero(null_density >= obs_density) + 1) / (len(null_density) + 1)
        prob_p.append(prob_tail)
        density_p.append(density_tail)
    out['empirical_p_greater_prob_neighbor_within_x'] = prob_p
    out['empirical_p_greater_other_methyl_args_within_x'] = density_p
    return out


def protein_density_summary(positions_by_acc: dict[str, np.ndarray], windows: list[int]) -> pd.DataFrame:
    rows = []
    for accession, positions in positions_by_acc.items():
        positions = np.sort(np.asarray(positions, dtype=int))
        methyl_site_count = len(positions)
        nearest = nearest_neighbor_distances(positions)
        row = {
            'canonical_UniProtAC': accession,
            'methyl_site_count': methyl_site_count,
            'nearest_distance_min': float(np.nanmin(nearest)) if np.isfinite(nearest).any() else pd.NA,
            'nearest_distance_median': float(np.nanmedian(nearest[np.isfinite(nearest)])) if np.isfinite(nearest).any() else pd.NA,
        }
        if methyl_site_count <= 1:
            for window in windows:
                row[f'max_methyl_sites_in_{window}aa_window'] = methyl_site_count
            rows.append(row)
            continue

        for window in windows:
            best = 0
            left = 0
            for right in range(methyl_site_count):
                while positions[right] - positions[left] > window:
                    left += 1
                best = max(best, right - left + 1)
            row[f'max_methyl_sites_in_{window}aa_window'] = int(best)
        rows.append(row)
    return pd.DataFrame(rows).sort_values(['methyl_site_count', 'nearest_distance_min'], ascending=[False, True])


def checkpoint_summary(curve_df: pd.DataFrame, checkpoints: list[int]) -> pd.DataFrame:
    selected = sorted(set(checkpoint for checkpoint in checkpoints if checkpoint > 0))
    return curve_df[curve_df['distance_aa'].isin(selected)].copy()


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
    args = ap.parse_args()

    sites = pd.read_csv(args.integrated_sites, sep='\t', low_memory=False)
    sequences = parse_fasta(args.canonical_fasta)
    sites = validate_sites(sites, sequences=sequences)

    methylated_proteins = sorted(set(sites['canonical_UniProtAC'].astype(str)))
    arg_positions_by_acc = {
        accession: residue_positions(sequences[accession], residue='R')
        for accession in methylated_proteins
        if accession in sequences
    }
    methyl_counts_by_acc = sites.groupby('canonical_UniProtAC').size().to_dict()
    observed_positions_by_acc = {
        accession: np.sort(group['corrected_position'].to_numpy(dtype=int))
        for accession, group in sites.groupby('canonical_UniProtAC')
    }

    observed_nearest_hist, observed_ordered_pair_hist, per_site = histograms_from_positions(
        observed_positions_by_acc,
        max_distance=args.max_distance,
    )
    observed_curve = summarize_curve_from_histograms(
        nearest_hist=observed_nearest_hist,
        ordered_pair_hist=observed_ordered_pair_hist,
        total_sites=len(sites),
        max_distance=args.max_distance,
    )

    permutations = permutation_summary(
        arg_positions_by_acc=arg_positions_by_acc,
        methyl_counts_by_acc=methyl_counts_by_acc,
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

    protein_summary = protein_density_summary(observed_positions_by_acc, windows=args.checkpoints)
    if 'substrate_genename' in sites.columns:
        gene_map = (
            sites[['canonical_UniProtAC', 'substrate_genename']]
            .drop_duplicates('canonical_UniProtAC')
            .rename(columns={'substrate_genename': 'gene'})
        )
        protein_summary = protein_summary.merge(gene_map, how='left', on='canonical_UniProtAC')

    per_site = per_site.merge(
        sites[['canonical_UniProtAC', 'corrected_position', 'substrate_genename', 'confidence_tier']].drop_duplicates(
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

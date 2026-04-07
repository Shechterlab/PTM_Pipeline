from __future__ import annotations

import numpy as np
import pandas as pd


def residue_positions(sequence: str, residues: set[str] | list[str] | tuple[str, ...]) -> np.ndarray:
    seq = str(sequence).strip().upper()
    residue_set = {str(residue).upper() for residue in residues}
    return np.array([idx for idx, aa in enumerate(seq, start=1) if aa in residue_set], dtype=int)


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
    n = len(positions)
    if n <= 1 or max_distance < 1:
        return hist

    for left in range(n):
        right = left + 1
        while right < n:
            distance = positions[right] - positions[left]
            if distance > max_distance:
                break
            hist[int(distance)] += 2
            right += 1
    return hist


def histograms_from_positions(
    positions_by_acc: dict[str, np.ndarray],
    max_distance: int,
    position_col: str = 'position',
    neighbor_distance_col: str = 'nearest_neighbor_distance',
    neighbor_flag_col: str = 'has_neighbor_within_max_distance',
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
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
                position_col: int(position),
                neighbor_distance_col: finite_distance,
                neighbor_flag_col: bool(np.isfinite(nearest_distance) and nearest_distance <= max_distance),
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


def sample_positions(
    candidate_positions_by_acc: dict[str, np.ndarray],
    site_counts_by_acc: dict[str, int],
    rng: np.random.Generator,
) -> dict[str, np.ndarray]:
    sampled = {}
    for accession, site_count in site_counts_by_acc.items():
        candidate_positions = candidate_positions_by_acc.get(accession, np.array([], dtype=int))
        if site_count <= 0 or len(candidate_positions) == 0:
            continue
        if site_count >= len(candidate_positions):
            sampled_positions = np.sort(candidate_positions.copy())
        else:
            sampled_positions = np.sort(rng.choice(candidate_positions, size=site_count, replace=False))
        sampled[accession] = sampled_positions.astype(int)
    return sampled


def permutation_summary(
    candidate_positions_by_acc: dict[str, np.ndarray],
    site_counts_by_acc: dict[str, int],
    max_distance: int,
    permutations: int,
    seed: int,
    progress_every: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for iteration in range(1, permutations + 1):
        sampled = sample_positions(candidate_positions_by_acc, site_counts_by_acc, rng)
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
    records = []
    for distance, group in permutation_df.groupby('distance_aa', as_index=False):
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
        site_count = len(positions)
        nearest = nearest_neighbor_distances(positions)
        row = {
            'canonical_UniProtAC': accession,
            'site_count': site_count,
            'nearest_distance_min': float(np.nanmin(nearest)) if np.isfinite(nearest).any() else pd.NA,
            'nearest_distance_median': float(np.nanmedian(nearest[np.isfinite(nearest)])) if np.isfinite(nearest).any() else pd.NA,
        }
        if site_count <= 1:
            for window in windows:
                row[f'max_sites_in_{window}aa_window'] = site_count
            rows.append(row)
            continue

        for window in windows:
            best = 0
            left = 0
            for right in range(site_count):
                while positions[right] - positions[left] > window:
                    left += 1
                best = max(best, right - left + 1)
            row[f'max_sites_in_{window}aa_window'] = int(best)
        rows.append(row)
    return pd.DataFrame(rows).sort_values(['site_count', 'nearest_distance_min'], ascending=[False, True])


def checkpoint_summary(curve_df: pd.DataFrame, checkpoints: list[int]) -> pd.DataFrame:
    selected = sorted(set(checkpoint for checkpoint in checkpoints if checkpoint > 0))
    return curve_df[curve_df['distance_aa'].isin(selected)].copy()

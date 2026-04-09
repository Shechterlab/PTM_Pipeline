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
    directed_histograms_from_position_maps,
    directed_permutation_summary,
    residue_positions,
    summarize_directed_curve_from_histograms,
)
from common import add_bh_q_values, canonicalize_uniprot_accession, parse_fasta, save_figure, save_table

DEFAULT_COMPARATOR_GROUPS = [
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



def build_position_maps(df: pd.DataFrame) -> dict[str, np.ndarray]:
    return {
        accession: np.sort(group['position'].to_numpy(dtype=int))
        for accession, group in df.groupby('canonical_UniProtAC')
    }



def low_complexity_mask(sequence: str, window: int = 12, max_unique: int = 3) -> np.ndarray:
    seq = str(sequence).strip().upper()
    if not seq:
        return np.zeros(0, dtype=bool)
    mask = np.zeros(len(seq), dtype=bool)
    if len(seq) <= window:
        mask[:] = len(set(seq)) <= max_unique
        return mask
    for start in range(0, len(seq) - window + 1):
        chunk = seq[start:start + window]
        if len(set(chunk)) <= max_unique:
            mask[start:start + window] = True
    return mask



def normalize_disorder_intervals(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    rename_map = {}
    for source, target in [('accession', 'canonical_UniProtAC'), ('Entry', 'canonical_UniProtAC'), ('start', 'fragment_start'), ('end', 'fragment_end')]:
        if source in out.columns and target not in out.columns:
            rename_map[source] = target
    out = out.rename(columns=rename_map)
    required = ['canonical_UniProtAC', 'fragment_start', 'fragment_end']
    missing = [col for col in required if col not in out.columns]
    if missing:
        raise ValueError(f'missing disorder interval columns: {missing}')
    out['canonical_UniProtAC'] = out['canonical_UniProtAC'].astype(str).map(canonicalize_uniprot_accession)
    out['fragment_start'] = pd.to_numeric(out['fragment_start'], errors='coerce')
    out['fragment_end'] = pd.to_numeric(out['fragment_end'], errors='coerce')
    out = out[out['fragment_start'].notna() & out['fragment_end'].notna()].copy()
    out['fragment_start'] = out['fragment_start'].astype(int)
    out['fragment_end'] = out['fragment_end'].astype(int)
    return out



def build_position_annotation_table(
    arg_sites: pd.DataFrame,
    sequences: dict[str, str],
    disorder_intervals: pd.DataFrame | None,
    low_complexity_window: int,
    low_complexity_max_unique: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    observed_by_acc = build_position_maps(arg_sites)
    disorder_masks: dict[str, np.ndarray] = {}
    if disorder_intervals is not None and not disorder_intervals.empty:
        for accession, group in disorder_intervals.groupby('canonical_UniProtAC'):
            seq_len = len(sequences.get(accession, ''))
            if seq_len == 0:
                continue
            mask = np.zeros(seq_len, dtype=bool)
            for row in group.itertuples(index=False):
                start = max(1, int(getattr(row, 'fragment_start'))) - 1
                end = min(seq_len, int(getattr(row, 'fragment_end')))
                mask[start:end] = True
            disorder_masks[accession] = mask

    low_complexity_masks = {accession: low_complexity_mask(seq, window=low_complexity_window, max_unique=low_complexity_max_unique) for accession, seq in sequences.items()}

    source_rows = []
    background_rows = []
    for accession, observed_positions in observed_by_acc.items():
        seq = sequences.get(accession, '')
        if not seq:
            continue
        arg_positions = residue_positions(seq, {'R'})
        observed_set = set(int(pos) for pos in observed_positions)
        disorder_mask = disorder_masks.get(accession, np.zeros(len(seq), dtype=bool))
        lc_mask = low_complexity_masks.get(accession, np.zeros(len(seq), dtype=bool))
        for position in arg_positions:
            is_observed = int(position) in observed_set
            row = {
                'canonical_UniProtAC': accession,
                'position': int(position),
                'is_disordered': bool(disorder_mask[position - 1]) if position - 1 < len(disorder_mask) else False,
                'is_local_low_complexity': bool(lc_mask[position - 1]) if position - 1 < len(lc_mask) else False,
            }
            if is_observed:
                source_rows.append(row)
            else:
                background_rows.append(row)
    return pd.DataFrame(source_rows), pd.DataFrame(background_rows)



def summarize_against_background(
    background_positions_by_acc: dict[str, np.ndarray],
    comparator_positions_by_acc: dict[str, np.ndarray],
    max_distance: int,
    exclude_self: bool,
) -> pd.DataFrame:
    nearest_hist, burden_hist, _ = directed_histograms_from_position_maps(
        source_positions_by_acc=background_positions_by_acc,
        target_positions_by_acc=comparator_positions_by_acc,
        max_distance=max_distance,
        source_position_col='background_arg_position',
        neighbor_distance_col='nearest_target_distance',
        neighbor_flag_col='has_target_within_max_distance',
        exclude_self=exclude_self,
    )
    total_background = int(sum(len(values) for values in background_positions_by_acc.values()))
    return summarize_directed_curve_from_histograms(
        nearest_hist=nearest_hist,
        burden_hist=burden_hist,
        total_sites=total_background,
        max_distance=max_distance,
    )



def position_maps_from_frame(df: pd.DataFrame) -> dict[str, np.ndarray]:
    if df.empty:
        return {}
    return {
        accession: np.sort(group['position'].to_numpy(dtype=int))
        for accession, group in df.groupby('canonical_UniProtAC')
    }



def strata_definitions(source_annotations: pd.DataFrame, background_annotations: pd.DataFrame) -> list[tuple[str, str, pd.DataFrame, pd.DataFrame]]:
    strata = [('all_sites', 'all_sites', source_annotations.copy(), background_annotations.copy())]
    has_disorder_signal = bool(source_annotations.get('is_disordered', pd.Series(dtype=bool)).any() or background_annotations.get('is_disordered', pd.Series(dtype=bool)).any())
    for value, label in [(True, 'disordered_source_sites'), (False, 'ordered_source_sites')]:
        if has_disorder_signal:
            strata.append((
                'disorder',
                label,
                source_annotations[source_annotations['is_disordered'] == value].copy(),
                background_annotations[background_annotations['is_disordered'] == value].copy(),
            ))
    for value, label in [(True, 'local_low_complexity_source_sites'), (False, 'non_low_complexity_source_sites')]:
        strata.append((
            'local_low_complexity',
            label,
            source_annotations[source_annotations['is_local_low_complexity'] == value].copy(),
            background_annotations[background_annotations['is_local_low_complexity'] == value].copy(),
        ))
    return strata



def plot_metric(curves: pd.DataFrame, value_col: str, ylabel: str, title: str, outpath: Path) -> None:
    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    for comparator, group in curves.groupby('comparator_ptm_group'):
        group = group.sort_values('distance_aa')
        ax.plot(group['distance_aa'], group[value_col], linewidth=2, label=comparator)
    ax.set_xlabel('Distance X from methyl-Arg (aa)')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    save_figure(fig, outpath)



def plot_checkpoint_heatmap(df: pd.DataFrame, value_col: str, significance_col: str, title: str, outpath: Path) -> None:
    if df.empty:
        return
    plot_df = df.copy()
    plot_df['row_label'] = plot_df['comparator_ptm_group']
    if 'stratum_label' in plot_df.columns and plot_df['stratum_label'].nunique() > 1:
        plot_df['row_label'] = plot_df['stratum_label'] + ' | ' + plot_df['comparator_ptm_group']
    heat = plot_df.pivot(index='row_label', columns='distance_aa', values=value_col)
    sig = plot_df.pivot(index='row_label', columns='distance_aa', values=significance_col)
    fig_height = max(4.8, 0.35 * len(heat.index) + 2.0)
    fig, ax = plt.subplots(figsize=(9.5, fig_height))
    im = ax.imshow(heat.to_numpy(dtype=float), aspect='auto')
    ax.set_xticks(np.arange(len(heat.columns)))
    ax.set_xticklabels(heat.columns)
    ax.set_yticks(np.arange(len(heat.index)))
    ax.set_yticklabels(heat.index)
    ax.set_xlabel('Distance checkpoint (aa)')
    ax.set_ylabel('Comparator / stratum')
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(value_col)
    for i in range(len(heat.index)):
        for j in range(len(heat.columns)):
            val = heat.iloc[i, j]
            star = ''
            q = sig.iloc[i, j]
            if pd.notna(q) and float(q) < 0.05:
                star = '*'
            ax.text(j, i, f'{val:.2f}{star}', ha='center', va='center', fontsize=8)
    fig.tight_layout()
    save_figure(fig, outpath)



def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--base-master', required=True)
    ap.add_argument('--integrated-arg-sites', required=True)
    ap.add_argument('--canonical-fasta', required=True)
    ap.add_argument('--disorder-intervals', default='')
    ap.add_argument('--outdir', default='results/cross_ptm_neighbors')
    ap.add_argument('--comparator-groups', nargs='*', default=list(DEFAULT_COMPARATOR_GROUPS))
    ap.add_argument('--max-distance', type=int, default=100)
    ap.add_argument('--permutations', type=int, default=250)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--progress-every', type=int, default=25)
    ap.add_argument('--checkpoints', nargs='*', type=int, default=[3, 5, 10, 20, 50])
    ap.add_argument('--low-complexity-window', type=int, default=12)
    ap.add_argument('--low-complexity-max-unique', type=int, default=3)
    args = ap.parse_args()

    sequences = parse_fasta(args.canonical_fasta)
    master = pd.read_csv(args.base_master, sep='\t', low_memory=False)
    integrated_arg = pd.read_csv(args.integrated_arg_sites, sep='\t', low_memory=False)

    arg_sites = normalize_arg_methyl_sites(integrated_arg, sequences=sequences)
    comparator_groups = [group for group in args.comparator_groups if group != 'Arg methylation']
    other_sites = normalize_master_sites(master, sequences=sequences, ptm_groups=comparator_groups)
    comparator_sites = pd.concat([arg_sites, other_sites], ignore_index=True, sort=False)

    disorder_intervals = None
    if args.disorder_intervals:
        disorder_intervals = normalize_disorder_intervals(pd.read_csv(args.disorder_intervals, sep='\t', low_memory=False))

    source_annotations, background_annotations = build_position_annotation_table(
        arg_sites=arg_sites,
        sequences=sequences,
        disorder_intervals=disorder_intervals,
        low_complexity_window=args.low_complexity_window,
        low_complexity_max_unique=args.low_complexity_max_unique,
    )

    curve_rows = []
    permutation_rows = []
    checkpoint_rows = []
    summary_rows = []
    per_site_rows = []

    for stratum_type, stratum_label, source_subset, background_subset in strata_definitions(source_annotations, background_annotations):
        if source_subset.empty:
            continue
        source_positions_by_acc = position_maps_from_frame(source_subset)
        background_positions_by_acc = position_maps_from_frame(background_subset)
        source_accessions = set(source_positions_by_acc)
        if not source_accessions:
            continue

        for comparator_ptm_group, group_all in comparator_sites.groupby('ptm_group'):
            group = group_all[group_all['canonical_UniProtAC'].isin(source_accessions)].copy()
            comparator_positions_by_acc = position_maps_from_frame(group)
            exclude_self = comparator_ptm_group == 'Arg methylation'
            nearest_hist, burden_hist, per_site = directed_histograms_from_position_maps(
                source_positions_by_acc=source_positions_by_acc,
                target_positions_by_acc=comparator_positions_by_acc,
                max_distance=args.max_distance,
                source_position_col='methyl_arg_position',
                neighbor_distance_col='nearest_target_ptm_distance',
                neighbor_flag_col='has_target_ptm_within_max_distance',
                exclude_self=exclude_self,
            )
            observed_curve = summarize_directed_curve_from_histograms(
                nearest_hist=nearest_hist,
                burden_hist=burden_hist,
                total_sites=int(len(source_subset)),
                max_distance=args.max_distance,
            )

            candidate_targets_by_acc = {
                accession: residue_positions(sequences[accession], PTM_GROUP_TARGET_RESIDUES[comparator_ptm_group])
                for accession in source_accessions
                if accession in sequences
            }
            target_site_counts_by_acc = group.groupby('canonical_UniProtAC').size().to_dict()
            null_permutations = directed_permutation_summary(
                source_positions_by_acc=source_positions_by_acc,
                candidate_target_positions_by_acc=candidate_targets_by_acc,
                target_site_counts_by_acc=target_site_counts_by_acc,
                max_distance=args.max_distance,
                permutations=args.permutations,
                seed=args.seed,
                progress_every=args.progress_every,
                exclude_self=exclude_self,
            )
            null_curve = aggregate_null_curve(null_permutations)
            matched_background_curve = summarize_against_background(
                background_positions_by_acc=background_positions_by_acc,
                comparator_positions_by_acc=comparator_positions_by_acc,
                max_distance=args.max_distance,
                exclude_self=False,
            )
            matched_background_curve = matched_background_curve.rename(columns={
                'prob_neighbor_within_x': 'matched_background_prob_neighbor_within_x',
                'mean_other_modified_sites_within_x': 'matched_background_mean_other_modified_sites_within_x',
                'sites_with_neighbor_within_x': 'matched_background_sites_with_neighbor_within_x',
                'ordered_neighbor_pairs_within_x': 'matched_background_ordered_neighbor_pairs_within_x',
            })

            curve = observed_curve.merge(null_curve, how='left', on='distance_aa').merge(
                matched_background_curve[['distance_aa', 'matched_background_prob_neighbor_within_x', 'matched_background_mean_other_modified_sites_within_x']],
                how='left',
                on='distance_aa',
            )
            curve = add_empirical_p_values(curve, permutation_df=null_permutations)
            curve['delta_prob_neighbor_within_x_vs_null'] = curve['prob_neighbor_within_x'] - curve['null_mean_prob_neighbor_within_x']
            curve['ratio_prob_neighbor_within_x_vs_null'] = curve['prob_neighbor_within_x'] / curve['null_mean_prob_neighbor_within_x'].replace(0, np.nan)
            curve['delta_prob_neighbor_within_x_vs_matched_background'] = curve['prob_neighbor_within_x'] - curve['matched_background_prob_neighbor_within_x']
            curve['ratio_prob_neighbor_within_x_vs_matched_background'] = curve['prob_neighbor_within_x'] / curve['matched_background_prob_neighbor_within_x'].replace(0, np.nan)
            curve['delta_other_modified_sites_within_x_vs_null'] = curve['mean_other_modified_sites_within_x'] - curve['null_mean_other_modified_sites_within_x']
            curve['ratio_other_modified_sites_within_x_vs_null'] = curve['mean_other_modified_sites_within_x'] / curve['null_mean_other_modified_sites_within_x'].replace(0, np.nan)
            curve['delta_other_modified_sites_within_x_vs_matched_background'] = curve['mean_other_modified_sites_within_x'] - curve['matched_background_mean_other_modified_sites_within_x']
            curve['ratio_other_modified_sites_within_x_vs_matched_background'] = curve['mean_other_modified_sites_within_x'] / curve['matched_background_mean_other_modified_sites_within_x'].replace(0, np.nan)
            curve['comparator_ptm_group'] = comparator_ptm_group
            curve['stratum_type'] = stratum_type
            curve['stratum_label'] = stratum_label
            curve_rows.append(curve)

            checkpoints = checkpoint_summary(curve, args.checkpoints)
            checkpoints['comparator_ptm_group'] = comparator_ptm_group
            checkpoints['stratum_type'] = stratum_type
            checkpoints['stratum_label'] = stratum_label
            checkpoint_rows.append(checkpoints)

            null_permutations['comparator_ptm_group'] = comparator_ptm_group
            null_permutations['stratum_type'] = stratum_type
            null_permutations['stratum_label'] = stratum_label
            permutation_rows.append(null_permutations)

            per_site['comparator_ptm_group'] = comparator_ptm_group
            per_site['stratum_type'] = stratum_type
            per_site['stratum_label'] = stratum_label
            per_site_rows.append(per_site)

            summary_rows.append({
                'stratum_type': stratum_type,
                'stratum_label': stratum_label,
                'comparator_ptm_group': comparator_ptm_group,
                'methyl_arg_sites': int(len(source_subset)),
                'methyl_arg_proteins': int(len(source_positions_by_acc)),
                'comparator_sites_on_methyl_arg_proteins': int(len(group)),
                'comparator_proteins_on_methyl_arg_proteins': int(group['canonical_UniProtAC'].nunique()),
                'matched_background_args': int(sum(len(values) for values in background_positions_by_acc.values())),
                'candidate_target_residues': int(sum(len(values) for values in candidate_targets_by_acc.values())),
                'exclude_self_in_nearest_neighbor': bool(exclude_self),
            })

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    curves = pd.concat(curve_rows, ignore_index=True)
    checkpoints = pd.concat(checkpoint_rows, ignore_index=True)
    permutations = pd.concat(permutation_rows, ignore_index=True)
    summary = pd.DataFrame(summary_rows).sort_values(['stratum_label', 'comparator_sites_on_methyl_arg_proteins'], ascending=[True, False])
    per_site = pd.concat(per_site_rows, ignore_index=True)

    curves = add_bh_q_values(curves, p_col='empirical_p_greater_prob_neighbor_within_x', q_col='global_empirical_q_greater_prob_neighbor_within_x')
    curves = add_bh_q_values(curves, p_col='empirical_p_greater_other_modified_sites_within_x', q_col='global_empirical_q_greater_other_modified_sites_within_x')
    checkpoints = curves[curves['distance_aa'].isin(sorted(set(checkpoint for checkpoint in args.checkpoints if checkpoint > 0)))].copy()

    save_table(summary, outdir / 'cross_ptm_neighbor_summary.tsv')
    save_table(curves, outdir / 'cross_ptm_neighbor_curves.tsv')
    save_table(checkpoints, outdir / 'cross_ptm_neighbor_checkpoints.tsv')
    save_table(permutations, outdir / 'cross_ptm_neighbor_null_permutations.tsv')
    save_table(per_site, outdir / 'cross_ptm_neighbor_per_site.tsv')
    save_table(source_annotations, outdir / 'methyl_arg_source_site_strata.tsv')
    save_table(background_annotations, outdir / 'non_methyl_arg_background_site_strata.tsv')

    all_curves = curves[curves['stratum_label'] == 'all_sites'].copy()
    if not all_curves.empty:
        plot_metric(
            curves=all_curves,
            value_col='prob_neighbor_within_x',
            ylabel='P(nearest comparator PTM <= X)',
            title='Methyl-Arg directed PTM-neighbor probability curves',
            outpath=outdir / 'cross_ptm_neighbor_probability_curves',
        )
        plot_metric(
            curves=all_curves,
            value_col='ratio_prob_neighbor_within_x_vs_null',
            ylabel='Observed / null nearest-neighbor probability',
            title='Methyl-Arg PTM-neighbor enrichment vs permutation null',
            outpath=outdir / 'cross_ptm_neighbor_enrichment_vs_null',
        )
        plot_metric(
            curves=all_curves,
            value_col='ratio_prob_neighbor_within_x_vs_matched_background',
            ylabel='Observed / matched-background nearest-neighbor probability',
            title='Methyl-Arg PTM-neighbor enrichment vs non-methyl-Arg background',
            outpath=outdir / 'cross_ptm_neighbor_enrichment_vs_background',
        )
        plot_metric(
            curves=all_curves,
            value_col='mean_other_modified_sites_within_x',
            ylabel='Mean comparator PTM sites within X aa',
            title='Methyl-Arg local heterotypic PTM burden',
            outpath=outdir / 'cross_ptm_neighbor_local_burden_curves',
        )
        plot_checkpoint_heatmap(
            df=checkpoints[checkpoints['stratum_label'] == 'all_sites'],
            value_col='ratio_prob_neighbor_within_x_vs_null',
            significance_col='global_empirical_q_greater_prob_neighbor_within_x',
            title='Methyl-Arg neighbor enrichment heatmap (all sites)',
            outpath=outdir / 'cross_ptm_neighbor_checkpoint_heatmap',
        )

    stratified_heatmap_df = checkpoints[checkpoints['distance_aa'] == min(args.checkpoints)].copy()
    if stratified_heatmap_df['stratum_label'].nunique() > 1:
        plot_checkpoint_heatmap(
            df=stratified_heatmap_df,
            value_col='ratio_prob_neighbor_within_x_vs_null',
            significance_col='global_empirical_q_greater_prob_neighbor_within_x',
            title=f'Methyl-Arg neighbor enrichment at {min(args.checkpoints)} aa across strata',
            outpath=outdir / f'cross_ptm_neighbor_stratified_{min(args.checkpoints)}aa_heatmap',
        )

    print({
        'comparators': sorted(set(curves['comparator_ptm_group'].astype(str))),
        'strata': sorted(set(curves['stratum_label'].astype(str))),
        'methyl_arg_sites': int(len(arg_sites)),
        'output': str(outdir),
    })


if __name__ == '__main__':
    main()

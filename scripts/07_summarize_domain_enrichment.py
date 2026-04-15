from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

from common import (
    BREWER_COLORS,
    add_q_values,
    annotate_barh,
    apply_paper_style,
    balanced_category_subset,
    classify_by_patterns,
    fisher_like_enrichment,
    format_p_value,
    load_json,
    log2_odds_ratio,
    parse_fasta,
    residue_background_from_sequences,
    save_figure,
    save_table,
    set_symmetric_xlim,
    signed_bar_colors,
    style_axis,
)


TYPE_PRIORITY = {
    'domain': 0,
    'repeat': 1,
    'homologous_superfamily': 2,
    'family': 3,
}

DOMAIN_EDGE_BIN_ORDER = [
    'Domain edge <=5 aa',
    'Domain edge 6-10 aa',
    'Domain edge 11-20 aa',
    'Domain edge 21-40 aa',
    'Domain-distal >40 aa',
]

DOMAIN_EDGE_BIN_COMPACT_ORDER = [
    'Domain edge <=20 aa',
    'Domain edge 21-40 aa',
    'Domain-distal >40 aa',
]


def interval_priority(interpro_type: str) -> int:
    return TYPE_PRIORITY.get(str(interpro_type), 9)


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
        }
    best, distance = choose_interval(intervals, position)
    if best is None:
        return {
            'domain_context_class': 'no_domain_annotation',
            'nearest_domain_class': 'No domain annotation',
            'nearest_domain_edge_distance': pd.NA,
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
    }


def annotate_background(background: pd.DataFrame, intervals_by_acc: dict[str, pd.DataFrame], boundary_window: int, ontology: dict) -> pd.DataFrame:
    out = background.copy()
    annotation_rows: list[pd.DataFrame] = []
    for accession, group in out.groupby('canonical_UniProtAC', sort=False):
        rows = []
        for row in group.itertuples():
            rows.append(
                classify_residue(
                    intervals_by_acc=intervals_by_acc,
                    accession=str(accession),
                    position=int(row.position),
                    boundary_window=boundary_window,
                    ontology=ontology,
                )
            )
        annotation_rows.append(pd.DataFrame(rows, index=group.index))
    annotations = pd.concat(annotation_rows, axis=0).sort_index() if annotation_rows else pd.DataFrame(index=out.index)
    return pd.concat([out, annotations], axis=1)


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


def domain_edge_bin(distance) -> str | pd.NA:
    if pd.isna(distance):
        return pd.NA
    distance = float(distance)
    if distance <= 5:
        return DOMAIN_EDGE_BIN_ORDER[0]
    if distance <= 10:
        return DOMAIN_EDGE_BIN_ORDER[1]
    if distance <= 20:
        return DOMAIN_EDGE_BIN_ORDER[2]
    if distance <= 40:
        return DOMAIN_EDGE_BIN_ORDER[3]
    return DOMAIN_EDGE_BIN_ORDER[4]


def domain_edge_bin_compact(distance) -> str | pd.NA:
    if pd.isna(distance):
        return pd.NA
    distance = float(distance)
    if distance <= 20:
        return DOMAIN_EDGE_BIN_COMPACT_ORDER[0]
    if distance <= 40:
        return DOMAIN_EDGE_BIN_COMPACT_ORDER[1]
    return DOMAIN_EDGE_BIN_COMPACT_ORDER[2]


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


def classify_disorder_positions(intervals: pd.DataFrame, positions: np.ndarray, boundary_window: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
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
    has_annotation[:] = True
    return classes, distances, has_annotation


def annotate_disorder_background(background: pd.DataFrame, intervals_by_acc: dict[str, pd.DataFrame], boundary_window: int) -> pd.DataFrame:
    out = background.copy()
    out['disorder_context_class'] = 'no_disorder_annotation'
    out['nearest_disorder_edge_distance'] = pd.NA
    out['has_disorder_annotation'] = False
    for accession, group in out.groupby('canonical_UniProtAC', sort=False):
        positions = pd.to_numeric(group['position'], errors='coerce').to_numpy(dtype=int)
        intervals = intervals_by_acc.get(str(accession), pd.DataFrame())
        classes, distances, has_annotation = classify_disorder_positions(intervals, positions, boundary_window)
        out.loc[group.index, 'disorder_context_class'] = classes
        out.loc[group.index, 'nearest_disorder_edge_distance'] = distances
        out.loc[group.index, 'has_disorder_annotation'] = has_annotation & (classes != 'no_disorder_annotation')
    return out


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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--annotated-sites', required=True)
    ap.add_argument('--interpro-intervals', required=True)
    ap.add_argument('--canonical-fasta', default='')
    ap.add_argument('--background-table', default='')
    ap.add_argument('--ontology', default='config/domain_class_ontology.json')
    ap.add_argument('--boundary-window', type=int, default=20)
    ap.add_argument('--include-types', nargs='*', default=['domain', 'repeat'])
    ap.add_argument('--disorder-sites', default='')
    ap.add_argument('--disorder-intervals', default='')
    ap.add_argument('--outdir', default='results/domain_context')
    args = ap.parse_args()

    annotated_path = Path(args.annotated_sites)
    interpro_path = Path(args.interpro_intervals)
    if not annotated_path.exists() or annotated_path.stat().st_size == 0:
        raise ValueError(f'annotated domain-context table is missing or empty: {annotated_path}')
    if not interpro_path.exists() or interpro_path.stat().st_size == 0:
        raise ValueError(f'InterPro interval table is missing or empty: {interpro_path}')

    annotated_sites = pd.read_csv(annotated_path, sep='\t', low_memory=False)
    intervals = pd.read_csv(interpro_path, sep='\t', low_memory=False)
    apply_paper_style()
    intervals['fragment_start'] = pd.to_numeric(intervals['fragment_start'], errors='coerce')
    intervals['fragment_end'] = pd.to_numeric(intervals['fragment_end'], errors='coerce')
    intervals = intervals[intervals['fragment_start'].notna() & intervals['fragment_end'].notna()].copy()
    if args.include_types:
        wanted = set(args.include_types)
        intervals = intervals[intervals['interpro_type'].astype(str).isin(wanted)].copy()
    intervals['fragment_start'] = intervals['fragment_start'].astype(int)
    intervals['fragment_end'] = intervals['fragment_end'].astype(int)
    intervals_by_acc = {
        accession: group.sort_values(['fragment_start', 'fragment_end', 'interpro_accession']).reset_index(drop=True)
        for accession, group in intervals.groupby('canonical_UniProtAC')
    }
    ontology = load_json(args.ontology)

    if 'nearest_domain_class' not in annotated_sites.columns:
        class_text = annotated_sites.fillna('').astype(str)[['nearest_interpro_name', 'nearest_interpro_type', 'nearest_member_database_accessions']].agg(' '.join, axis=1)
        annotated_sites['nearest_domain_class'] = class_text.map(lambda x: classify_by_patterns(x, ontology['domain_classes'], default='Other domain'))
    annotated_sites.loc[~annotated_sites['has_domain_annotation'].fillna(False), 'nearest_domain_class'] = 'No domain annotation'
    annotated_sites['domain_edge_bin'] = pd.NA
    annotated_sites['domain_edge_bin_compact'] = pd.NA
    annotated_domain_edge_mask = (
        annotated_sites['has_domain_annotation'].fillna(False) &
        annotated_sites['domain_context_class'].isin(['boundary', 'inter_domain_linker', 'distal']) &
        annotated_sites['nearest_domain_edge_distance'].notna()
    )
    annotated_sites.loc[annotated_domain_edge_mask, 'domain_edge_bin'] = (
        annotated_sites.loc[annotated_domain_edge_mask, 'nearest_domain_edge_distance'].map(domain_edge_bin)
    )
    annotated_sites.loc[annotated_domain_edge_mask, 'domain_edge_bin_compact'] = (
        annotated_sites.loc[annotated_domain_edge_mask, 'nearest_domain_edge_distance'].map(domain_edge_bin_compact)
    )

    disorder_intervals_by_acc: dict[str, pd.DataFrame] = {}
    if args.disorder_intervals:
        disorder_path = Path(args.disorder_intervals)
        if not disorder_path.exists() or disorder_path.stat().st_size == 0:
            raise ValueError(f'disorder interval table is missing or empty: {disorder_path}')
        disorder_intervals = normalize_disorder_intervals(pd.read_csv(disorder_path, sep='\t', low_memory=False))
        disorder_intervals_by_acc = {
            accession: group[['fragment_start', 'fragment_end']].reset_index(drop=True)
            for accession, group in disorder_intervals.groupby('canonical_UniProtAC')
        }

    if args.disorder_sites:
        disorder_sites_path = Path(args.disorder_sites)
        if not disorder_sites_path.exists() or disorder_sites_path.stat().st_size == 0:
            raise ValueError(f'disorder site table is missing or empty: {disorder_sites_path}')
        disorder_sites = pd.read_csv(disorder_sites_path, sep='\t', low_memory=False)
        required_cols = ['canonical_UniProtAC', 'corrected_position', 'disorder_context_class', 'nearest_disorder_edge_distance', 'has_disorder_annotation']
        missing = [col for col in required_cols if col not in disorder_sites.columns]
        if missing:
            raise ValueError(f'disorder site table missing columns: {missing}')
        disorder_sites = disorder_sites[required_cols].drop_duplicates(['canonical_UniProtAC', 'corrected_position']).copy()
        merge_cols = [col for col in ['disorder_context_class', 'nearest_disorder_edge_distance', 'has_disorder_annotation'] if col in annotated_sites.columns]
        if merge_cols:
            annotated_sites = annotated_sites.drop(columns=merge_cols)
        annotated_sites = annotated_sites.merge(
            disorder_sites,
            on=['canonical_UniProtAC', 'corrected_position'],
            how='left',
        )
    elif disorder_intervals_by_acc:
        position_col = 'corrected_position' if 'corrected_position' in annotated_sites.columns else 'position'
        disorder_frame = annotated_sites[['canonical_UniProtAC', position_col]].rename(columns={position_col: 'position'}).copy()
        disorder_frame = annotate_disorder_background(disorder_frame, disorder_intervals_by_acc, args.boundary_window)
        for col in ['disorder_context_class', 'nearest_disorder_edge_distance', 'has_disorder_annotation']:
            if col in annotated_sites.columns:
                annotated_sites = annotated_sites.drop(columns=[col])
        annotated_sites = pd.concat([annotated_sites, disorder_frame[['disorder_context_class', 'nearest_disorder_edge_distance', 'has_disorder_annotation']]], axis=1)

    if args.background_table:
        background_path = Path(args.background_table)
        if not background_path.exists() or background_path.stat().st_size == 0:
            raise ValueError(f'background domain-context table is missing or empty: {background_path}')
        methyl_protein_arg = pd.read_csv(background_path, sep='\t', low_memory=False)
    else:
        if not args.canonical_fasta:
            raise ValueError('--canonical-fasta is required when --background-table is not provided')
        sequences = parse_fasta(args.canonical_fasta)
        proteome_arg = residue_background_from_sequences(sequences, residue='R')
        methylated_proteins = sorted(set(annotated_sites['canonical_UniProtAC'].dropna().astype(str)))
        methyl_protein_arg = proteome_arg[proteome_arg['canonical_UniProtAC'].isin(methylated_proteins)].copy()
        methyl_site_keys = set(annotated_sites['site_key'].astype(str))
        methyl_protein_arg = methyl_protein_arg[~methyl_protein_arg['site_key'].isin(methyl_site_keys)].copy()
        methyl_protein_arg = annotate_background(methyl_protein_arg, intervals_by_acc, args.boundary_window, ontology)
    methyl_protein_arg['domain_edge_bin'] = pd.NA
    methyl_protein_arg['domain_edge_bin_compact'] = pd.NA
    background_domain_edge_mask = (
        methyl_protein_arg['domain_context_class'].isin(['boundary', 'inter_domain_linker', 'distal']) &
        methyl_protein_arg['nearest_domain_edge_distance'].notna()
    )
    methyl_protein_arg.loc[background_domain_edge_mask, 'domain_edge_bin'] = (
        methyl_protein_arg.loc[background_domain_edge_mask, 'nearest_domain_edge_distance'].map(domain_edge_bin)
    )
    methyl_protein_arg.loc[background_domain_edge_mask, 'domain_edge_bin_compact'] = (
        methyl_protein_arg.loc[background_domain_edge_mask, 'nearest_domain_edge_distance'].map(domain_edge_bin_compact)
    )
    if disorder_intervals_by_acc:
        methyl_protein_arg = annotate_disorder_background(methyl_protein_arg, disorder_intervals_by_acc, args.boundary_window)

    context_enrichment = fisher_like_enrichment(
        annotated_sites['domain_context_class'],
        methyl_protein_arg['domain_context_class'],
    )
    domain_class_enrichment = fisher_like_enrichment(
        annotated_sites['nearest_domain_class'],
        methyl_protein_arg['nearest_domain_class'],
    )
    edge_summary = annotated_sites[['nearest_domain_edge_distance', 'domain_context_class']].copy()
    observed_edge = annotated_sites['nearest_domain_edge_distance'].dropna().astype(float)
    background_edge = methyl_protein_arg['nearest_domain_edge_distance'].dropna().astype(float)
    edge_comparison = pd.DataFrame([{
        'observed_n': len(observed_edge),
        'background_n': len(background_edge),
        'observed_median': observed_edge.median() if not observed_edge.empty else pd.NA,
        'background_median': background_edge.median() if not background_edge.empty else pd.NA,
        'mannwhitney_u_p_value': float(mannwhitneyu(observed_edge, background_edge, alternative='two-sided').pvalue)
        if not observed_edge.empty and not background_edge.empty else pd.NA,
    }])
    edge_bin_enrichment = fisher_like_enrichment(
        annotated_sites[annotated_sites['domain_edge_bin'].notna()]['domain_edge_bin'],
        methyl_protein_arg[methyl_protein_arg['domain_edge_bin'].notna()]['domain_edge_bin'],
    )
    edge_bin_compact_enrichment = fisher_like_enrichment(
        annotated_sites[annotated_sites['domain_edge_bin_compact'].notna()]['domain_edge_bin_compact'],
        methyl_protein_arg[methyl_protein_arg['domain_edge_bin_compact'].notna()]['domain_edge_bin_compact'],
    )
    edge_curve = cumulative_edge_curve(
        annotated_sites.loc[annotated_domain_edge_mask, 'nearest_domain_edge_distance'],
        methyl_protein_arg.loc[background_domain_edge_mask, 'nearest_domain_edge_distance'],
        max_distance=40,
    )
    disordered_edge_bin_enrichment = pd.DataFrame()
    disordered_edge_bin_compact_enrichment = pd.DataFrame()
    disordered_edge_curve = pd.DataFrame()
    disordered_edge_qc = pd.DataFrame()
    if 'disorder_context_class' in annotated_sites.columns and 'disorder_context_class' in methyl_protein_arg.columns:
        annotated_disordered_edge_mask = annotated_domain_edge_mask & (annotated_sites['disorder_context_class'] == 'disordered')
        background_disordered_edge_mask = background_domain_edge_mask & (methyl_protein_arg['disorder_context_class'] == 'disordered')
        if int(annotated_disordered_edge_mask.sum()) > 0 and int(background_disordered_edge_mask.sum()) > 0:
            disordered_edge_bin_enrichment = fisher_like_enrichment(
                annotated_sites.loc[annotated_disordered_edge_mask, 'domain_edge_bin'],
                methyl_protein_arg.loc[background_disordered_edge_mask, 'domain_edge_bin'],
            )
            disordered_edge_bin_compact_enrichment = fisher_like_enrichment(
                annotated_sites.loc[annotated_disordered_edge_mask, 'domain_edge_bin_compact'],
                methyl_protein_arg.loc[background_disordered_edge_mask, 'domain_edge_bin_compact'],
            )
            disordered_edge_curve = cumulative_edge_curve(
                annotated_sites.loc[annotated_disordered_edge_mask, 'nearest_domain_edge_distance'],
                methyl_protein_arg.loc[background_disordered_edge_mask, 'nearest_domain_edge_distance'],
                max_distance=40,
            )
            disordered_edge_qc = pd.DataFrame([{
                'annotated_disordered_outside_domain_sites': int(annotated_disordered_edge_mask.sum()),
                'background_disordered_outside_domain_residues': int(background_disordered_edge_mask.sum()),
            }])

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(annotated_sites, outdir / 'sites_with_domain_context.tsv')
    save_table(context_enrichment, outdir / 'domain_context_enrichment.tsv')
    save_table(domain_class_enrichment, outdir / 'domain_class_enrichment.tsv')
    save_table(edge_summary, outdir / 'domain_edge_distance_summary.tsv')
    save_table(edge_comparison, outdir / 'domain_edge_distance_comparison.tsv')
    save_table(edge_bin_enrichment, outdir / 'domain_edge_bin_enrichment.tsv')
    save_table(edge_bin_compact_enrichment, outdir / 'domain_edge_bin_enrichment_compact.tsv')
    save_table(edge_curve, outdir / 'domain_edge_cumulative_enrichment.tsv')
    if not disordered_edge_bin_enrichment.empty:
        save_table(disordered_edge_bin_enrichment, outdir / 'domain_edge_bin_enrichment_disordered_only.tsv')
    if not disordered_edge_bin_compact_enrichment.empty:
        save_table(disordered_edge_bin_compact_enrichment, outdir / 'domain_edge_bin_enrichment_disordered_only_compact.tsv')
    if not disordered_edge_curve.empty:
        save_table(disordered_edge_curve, outdir / 'domain_edge_cumulative_enrichment_disordered_only.tsv')
    if not disordered_edge_qc.empty:
        save_table(disordered_edge_qc, outdir / 'domain_edge_disordered_only_qc.tsv')
    save_table(methyl_protein_arg, outdir / 'same_protein_nonmethyl_arginine_domain_background.tsv')
    save_table(
        pd.DataFrame([{
            'include_types': ';'.join(args.include_types) if args.include_types else 'all',
            'boundary_window': args.boundary_window,
            'annotated_sites': len(annotated_sites),
            'background_residues': len(methyl_protein_arg),
        }]),
        outdir / 'domain_analysis_qc.tsv',
    )

    plot_context = context_enrichment[
        context_enrichment['category'].isin(['in_domain', 'boundary', 'inter_domain_linker', 'distal', 'no_domain_annotation'])
    ].sort_values('log2_odds_ratio').copy()
    plot_context['display_category'] = plot_context['category'].map(domain_context_label)
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    context_palette = {
        'in_domain': BREWER_COLORS['blue'],
        'boundary': BREWER_COLORS['teal'],
        'inter_domain_linker': BREWER_COLORS['purple'],
        'distal': BREWER_COLORS['orange'],
        'no_domain_annotation': BREWER_COLORS['mid_gray'],
    }
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
    set_symmetric_xlim(ax, plot_context['log2_odds_ratio'], annotation_pad_ratio=0.75, center_on_zero=False)
    annotate_barh(
        ax,
        plot_context['display_category'],
        plot_context['log2_odds_ratio'],
        [f'n={n}, p={format_p_value(p)}, q={format_p_value(q)}' for n, p, q in zip(plot_context['target_count'], plot_context['p_value'], plot_context['q_value'])],
        fontsize=8.2,
    )
    fig.text(
        0.99,
        0.01,
        'Site-level counts using strict InterPro domain/repeat intervals.',
        ha='right',
        va='bottom',
        fontsize=8,
        color=BREWER_COLORS['dark_gray'],
    )
    fig.tight_layout()
    save_figure(fig, outdir / 'arg_methyl_domain_context_enrichment')

    plot_domain = balanced_category_subset(
        domain_class_enrichment,
        min_count=20,
        top_positive=7,
        top_negative=5,
        always_include=[],
    )
    plot_domain = plot_domain[~plot_domain['category'].isin({'Other domain', 'No domain annotation'})].copy()
    colors = signed_bar_colors(plot_domain['log2_odds_ratio'], positive=BREWER_COLORS['green'], negative=BREWER_COLORS['lavender'])
    colors = [BREWER_COLORS['mid_gray'] if cat in {'Other domain', 'No domain annotation'} else color for cat, color in zip(plot_domain['category'], colors)]
    fig2, ax2 = plt.subplots(figsize=(9.4, 6.2))
    plot_domain['display_category'] = plot_domain['category'].map(domain_class_label)
    ax2.barh(plot_domain['display_category'], plot_domain['log2_odds_ratio'], color=colors, edgecolor='white', linewidth=0.8)
    style_axis(ax2, zero='x')
    ax2.set_xlabel('log2(OR) vs arginines in methylated proteins')
    ax2.set_ylabel('Nearest domain class')
    ax2.set_title('Nearest domain-class enrichment around methylarginines')
    set_symmetric_xlim(ax2, plot_domain['log2_odds_ratio'], annotation_pad_ratio=0.82, center_on_zero=False)
    annotate_barh(
        ax2,
        plot_domain['display_category'],
        plot_domain['log2_odds_ratio'],
        [f'n={n}, p={format_p_value(p)}, q={format_p_value(q)}' for n, p, q in zip(plot_domain['target_count'], plot_domain['p_value'], plot_domain['q_value'])],
        fontsize=7.7,
    )
    fig2.text(
        0.99,
        0.01,
        'Site-level counts; selected enriched/depleted nearest-domain classes.',
        ha='right',
        va='bottom',
        fontsize=8,
        color=BREWER_COLORS['dark_gray'],
    )
    fig2.tight_layout()
    save_figure(fig2, outdir / 'arg_methyl_nearest_domain_class_enrichment')

    edge_plot = edge_bin_enrichment.copy()
    edge_plot['category'] = pd.Categorical(edge_plot['category'], categories=DOMAIN_EDGE_BIN_ORDER, ordered=True)
    edge_plot = edge_plot.sort_values('category').copy()
    edge_palette = {
        'Domain edge <=5 aa': BREWER_COLORS['teal'],
        'Domain edge 6-10 aa': BREWER_COLORS['light_green'],
        'Domain edge 11-20 aa': BREWER_COLORS['green'],
        'Domain edge 21-40 aa': BREWER_COLORS['gold'],
        'Domain-distal >40 aa': BREWER_COLORS['orange'],
    }
    fig3, ax3 = plt.subplots(figsize=(8.8, 5.4))
    ax3.barh(
        edge_plot['category'].astype(str),
        edge_plot['log2_odds_ratio'],
        color=[edge_palette.get(cat, BREWER_COLORS['mid_gray']) for cat in edge_plot['category'].astype(str)],
        edgecolor='white',
        linewidth=0.8,
    )
    style_axis(ax3, zero='x')
    ax3.set_xlabel('log2(OR) vs same-protein non-methyl arginines outside domains')
    ax3.set_ylabel('Outside-domain distance from nearest edge')
    ax3.set_title('Methylarginine proximity outside annotated domains')
    set_symmetric_xlim(ax3, edge_plot['log2_odds_ratio'], annotation_pad_ratio=0.62, center_on_zero=False)
    annotate_barh(
        ax3,
        edge_plot['category'].astype(str),
        edge_plot['log2_odds_ratio'],
        [f'n={n}, p={format_p_value(p)}, q={format_p_value(q)}' for n, p, q in zip(edge_plot['target_count'], edge_plot['p_value'], edge_plot['q_value'])],
        fontsize=8.0,
    )
    fig3.text(
        0.99,
        0.01,
        'Out-of-domain methylarginines only; compared with non-methyl arginines from the same methylated proteins.',
        ha='right',
        va='bottom',
        fontsize=8,
        color=BREWER_COLORS['dark_gray'],
    )
    fig3.tight_layout()
    save_figure(fig3, outdir / 'arg_methyl_domain_edge_bin_enrichment')

    if not disordered_edge_bin_enrichment.empty:
        disordered_edge_plot = disordered_edge_bin_enrichment.copy()
        disordered_edge_plot['category'] = pd.Categorical(disordered_edge_plot['category'], categories=DOMAIN_EDGE_BIN_ORDER, ordered=True)
        disordered_edge_plot = disordered_edge_plot.sort_values('category').copy()
        fig3b, ax3b = plt.subplots(figsize=(8.8, 5.4))
        ax3b.barh(
            disordered_edge_plot['category'].astype(str),
            disordered_edge_plot['log2_odds_ratio'],
            color=[edge_palette.get(cat, BREWER_COLORS['mid_gray']) for cat in disordered_edge_plot['category'].astype(str)],
            edgecolor='white',
            linewidth=0.8,
        )
        style_axis(ax3b, zero='x')
        ax3b.set_xlabel('log2(OR) vs disordered non-methyl arginines outside domains')
        ax3b.set_ylabel('Outside-domain distance from nearest edge')
        ax3b.set_title('Outside-domain proximity within disordered sequence')
        set_symmetric_xlim(ax3b, disordered_edge_plot['log2_odds_ratio'], annotation_pad_ratio=0.62, center_on_zero=False)
        annotate_barh(
            ax3b,
            disordered_edge_plot['category'].astype(str),
            disordered_edge_plot['log2_odds_ratio'],
            [f'n={n}, p={format_p_value(p)}, q={format_p_value(q)}' for n, p, q in zip(disordered_edge_plot['target_count'], disordered_edge_plot['p_value'], disordered_edge_plot['q_value'])],
            fontsize=8.0,
        )
        fig3b.text(
            0.99,
            0.01,
            'Restricted to methylarginine sites and background arginines that are already inside annotated IDRs.',
            ha='right',
            va='bottom',
            fontsize=8,
            color=BREWER_COLORS['dark_gray'],
        )
        fig3b.tight_layout()
        save_figure(fig3b, outdir / 'arg_methyl_domain_edge_bin_enrichment_disordered_only')

    if not disordered_edge_bin_compact_enrichment.empty:
        disordered_edge_plot_compact = disordered_edge_bin_compact_enrichment.copy()
        disordered_edge_plot_compact['category'] = pd.Categorical(
            disordered_edge_plot_compact['category'],
            categories=DOMAIN_EDGE_BIN_COMPACT_ORDER,
            ordered=True,
        )
        disordered_edge_plot_compact = disordered_edge_plot_compact.sort_values('category').copy()
        compact_palette = {
            'Domain edge <=20 aa': BREWER_COLORS['green'],
            'Domain edge 21-40 aa': BREWER_COLORS['gold'],
            'Domain-distal >40 aa': BREWER_COLORS['orange'],
        }
        fig3c, ax3c = plt.subplots(figsize=(7.8, 4.6))
        ax3c.barh(
            disordered_edge_plot_compact['category'].astype(str),
            disordered_edge_plot_compact['log2_odds_ratio'],
            color=[compact_palette.get(cat, BREWER_COLORS['mid_gray']) for cat in disordered_edge_plot_compact['category'].astype(str)],
            edgecolor='white',
            linewidth=0.8,
        )
        style_axis(ax3c, zero='x')
        ax3c.set_xlabel('log2(OR) vs disordered non-methyl arginines outside domains')
        ax3c.set_ylabel('Outside-domain distance from nearest edge')
        ax3c.set_title('Outside-domain proximity within disordered sequence')
        set_symmetric_xlim(ax3c, disordered_edge_plot_compact['log2_odds_ratio'], annotation_pad_ratio=0.52, center_on_zero=False)
        annotate_barh(
            ax3c,
            disordered_edge_plot_compact['category'].astype(str),
            disordered_edge_plot_compact['log2_odds_ratio'],
            [f'n={n}, p={format_p_value(p)}, q={format_p_value(q)}' for n, p, q in zip(disordered_edge_plot_compact['target_count'], disordered_edge_plot_compact['p_value'], disordered_edge_plot_compact['q_value'])],
            fontsize=8.0,
        )
        fig3c.text(
            0.99,
            0.01,
            'Disordered sites only; compact outside-domain bins for presentation.',
            ha='right',
            va='bottom',
            fontsize=8,
            color=BREWER_COLORS['dark_gray'],
        )
        fig3c.tight_layout()
        save_figure(fig3c, outdir / 'arg_methyl_domain_edge_bin_enrichment_disordered_only_compact')

        fraction_plot = disordered_edge_plot_compact.copy()
        fraction_plot['total_arginines_in_bin'] = fraction_plot['target_count'] + fraction_plot['background_count']
        fraction_plot['methylated_fraction_in_bin'] = fraction_plot['target_count'] / fraction_plot['total_arginines_in_bin']
        save_table(
            fraction_plot[['category', 'target_count', 'background_count', 'total_arginines_in_bin', 'methylated_fraction_in_bin', 'odds_ratio', 'log2_odds_ratio', 'p_value', 'q_value']],
            outdir / 'domain_edge_bin_enrichment_disordered_only_compact_fraction.tsv',
        )
        fig3d, ax3d = plt.subplots(figsize=(7.8, 4.6))
        ax3d.barh(
            fraction_plot['category'].astype(str),
            100.0 * fraction_plot['methylated_fraction_in_bin'],
            color=[compact_palette.get(cat, BREWER_COLORS['mid_gray']) for cat in fraction_plot['category'].astype(str)],
            edgecolor='white',
            linewidth=0.8,
        )
        style_axis(ax3d, grid_axis='x')
        ax3d.set_xlabel('% of arginines methylated in bin')
        ax3d.set_ylabel('Outside-domain distance from nearest edge')
        ax3d.set_title('Methylarginine frequency in disordered outside-domain bins')
        annotate_barh(
            ax3d,
            fraction_plot['category'].astype(str),
            100.0 * fraction_plot['methylated_fraction_in_bin'],
            [f"{n}/{t} ({100.0 * frac:.1f}%), q={format_p_value(q)}" for n, t, frac, q in zip(fraction_plot['target_count'], fraction_plot['total_arginines_in_bin'], fraction_plot['methylated_fraction_in_bin'], fraction_plot['q_value'])],
            fontsize=8.0,
        )
        fig3d.text(
            0.99,
            0.01,
            'Disordered outside-domain arginines only; denominator = methylated + non-methyl arginines in each bin.',
            ha='right',
            va='bottom',
            fontsize=8,
            color=BREWER_COLORS['dark_gray'],
        )
        fig3d.tight_layout()
        save_figure(fig3d, outdir / 'arg_methyl_domain_edge_bin_frequency_disordered_only_compact')

    fig4, ax4 = plt.subplots(figsize=(8.8, 4.8))
    ax4.plot(edge_curve['distance_aa'], edge_curve['log2_odds_ratio'], color=BREWER_COLORS['green'], linewidth=2.2)
    sig = edge_curve[edge_curve['q_value'] <= 0.05]
    if not sig.empty:
        ax4.scatter(sig['distance_aa'], sig['log2_odds_ratio'], color=BREWER_COLORS['orange'], s=18, zorder=3)
    style_axis(ax4, zero='y')
    ax4.set_xlabel('Outside-domain distance from nearest edge (<= X aa)')
    ax4.set_ylabel('log2(OR) vs same-protein non-methyl arginines')
    ax4.set_title('Cumulative enrichment outside domains by edge proximity')
    ax4.set_xlim(1, 40)
    for checkpoint in [5, 10, 20, 40]:
        row = edge_curve.loc[edge_curve['distance_aa'] == checkpoint]
        if row.empty:
            continue
        row = row.iloc[0]
        ax4.axvline(checkpoint, color=BREWER_COLORS['light_gray'], linestyle='--', linewidth=0.8)
        ax4.text(checkpoint, row['log2_odds_ratio'], f' {checkpoint}', fontsize=8, va='bottom', color=BREWER_COLORS['dark_gray'])
    fig4.text(
        0.99,
        0.01,
        'Uses out-of-domain sites with domain annotation only; points mark q<=0.05.',
        ha='right',
        va='bottom',
        fontsize=8,
        color=BREWER_COLORS['dark_gray'],
    )
    fig4.tight_layout()
    save_figure(fig4, outdir / 'arg_methyl_domain_edge_cumulative_enrichment')

    print({
        'annotated_sites': len(annotated_sites),
        'background_arginines_methylated_proteins': len(methyl_protein_arg),
        'outputs': str(outdir),
    })


if __name__ == '__main__':
    main()

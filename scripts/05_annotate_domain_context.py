from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


TYPE_PRIORITY = {
    'domain': 0,
    'repeat': 1,
    'homologous_superfamily': 2,
    'family': 3,
}


def interval_priority(interpro_type: str) -> int:
    return TYPE_PRIORITY.get(str(interpro_type), 9)


def edge_distance(position: int, start: int, end: int) -> tuple[int, str]:
    if start <= position <= end:
        return 0, 'inside'
    if position < start:
        return start - position, 'N_terminal_to_domain'
    return position - end, 'C_terminal_to_domain'


def choose_interval(intervals: pd.DataFrame, position: int) -> tuple[pd.Series | None, int, str, int]:
    if intervals.empty:
        return None, -1, '', 0

    overlaps = intervals[(intervals['fragment_start'] <= position) & (intervals['fragment_end'] >= position)].copy()
    if not overlaps.empty:
        overlaps['span'] = overlaps['fragment_end'] - overlaps['fragment_start']
        overlaps['priority'] = overlaps['interpro_type'].map(interval_priority)
        overlaps = overlaps.sort_values(['priority', 'span', 'fragment_start', 'interpro_accession'])
        best = overlaps.iloc[0]
        return best, 0, 'inside', len(overlaps)

    candidates = intervals.copy()
    distances = []
    sides = []
    for row in candidates.itertuples(index=False):
        distance, side = edge_distance(position, int(row.fragment_start), int(row.fragment_end))
        distances.append(distance)
        sides.append(side)
    candidates['distance'] = distances
    candidates['side'] = sides
    candidates['span'] = candidates['fragment_end'] - candidates['fragment_start']
    candidates['priority'] = candidates['interpro_type'].map(interval_priority)
    candidates = candidates.sort_values(['distance', 'priority', 'span', 'fragment_start', 'interpro_accession'])
    best = candidates.iloc[0]
    return best, int(best['distance']), str(best['side']), 0


def classify_row(row: pd.Series, intervals_by_accession: dict[str, pd.DataFrame], accession_col: str, position_col: str, boundary_window: int) -> dict:
    accession = str(row.get(accession_col, '') or '').strip()
    position_raw = row.get(position_col, '')
    position = pd.to_numeric(pd.Series([position_raw]), errors='coerce').iloc[0]
    if accession in {'', 'nan', 'None'} or pd.isna(position):
        return {
            'domain_context_class': 'no_position',
            'has_domain_annotation': False,
            'nearest_domain_edge_distance': pd.NA,
            'nearest_domain_side': '',
            'nearest_interpro_accession': '',
            'nearest_interpro_name': '',
            'nearest_interpro_type': '',
            'nearest_member_database_sources': '',
            'nearest_member_database_accessions': '',
            'nearest_domain_start': pd.NA,
            'nearest_domain_end': pd.NA,
            'overlapping_domain_count': 0,
        }

    position = int(position)
    intervals = intervals_by_accession.get(accession)
    if intervals is None or intervals.empty:
        return {
            'domain_context_class': 'no_domain_annotation',
            'has_domain_annotation': False,
            'nearest_domain_edge_distance': pd.NA,
            'nearest_domain_side': '',
            'nearest_interpro_accession': '',
            'nearest_interpro_name': '',
            'nearest_interpro_type': '',
            'nearest_member_database_sources': '',
            'nearest_member_database_accessions': '',
            'nearest_domain_start': pd.NA,
            'nearest_domain_end': pd.NA,
            'overlapping_domain_count': 0,
        }

    best, distance, side, overlap_count = choose_interval(intervals, position)
    if best is None:
        return {
            'domain_context_class': 'no_domain_annotation',
            'has_domain_annotation': False,
            'nearest_domain_edge_distance': pd.NA,
            'nearest_domain_side': '',
            'nearest_interpro_accession': '',
            'nearest_interpro_name': '',
            'nearest_interpro_type': '',
            'nearest_member_database_sources': '',
            'nearest_member_database_accessions': '',
            'nearest_domain_start': pd.NA,
            'nearest_domain_end': pd.NA,
            'overlapping_domain_count': 0,
        }

    if distance == 0:
        context_class = 'in_domain'
    elif distance <= boundary_window:
        context_class = 'boundary'
    else:
        context_class = 'distal'

    return {
        'domain_context_class': context_class,
        'has_domain_annotation': True,
        'nearest_domain_edge_distance': distance,
        'nearest_domain_side': side,
        'nearest_interpro_accession': best['interpro_accession'],
        'nearest_interpro_name': best['interpro_name'],
        'nearest_interpro_type': best['interpro_type'],
        'nearest_member_database_sources': best['member_database_sources'],
        'nearest_member_database_accessions': best['member_database_accessions'],
        'nearest_domain_start': int(best['fragment_start']),
        'nearest_domain_end': int(best['fragment_end']),
        'overlapping_domain_count': overlap_count,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--sites', required=True)
    ap.add_argument('--interpro-intervals', required=True)
    ap.add_argument('--outdir', default='results/domain_context')
    ap.add_argument('--accession-col', default='canonical_UniProtAC')
    ap.add_argument('--position-col', default='corrected_position')
    ap.add_argument('--boundary-window', type=int, default=20)
    ap.add_argument('--include-types', nargs='*', default=[])
    args = ap.parse_args()

    sites = pd.read_csv(args.sites, sep='\t', low_memory=False)
    intervals = pd.read_csv(args.interpro_intervals, sep='\t', low_memory=False)
    if args.include_types:
        wanted = set(args.include_types)
        intervals = intervals[intervals['interpro_type'].astype(str).isin(wanted)].copy()

    intervals['fragment_start'] = pd.to_numeric(intervals['fragment_start'], errors='coerce')
    intervals['fragment_end'] = pd.to_numeric(intervals['fragment_end'], errors='coerce')
    intervals = intervals[intervals['fragment_start'].notna() & intervals['fragment_end'].notna()].copy()
    intervals['fragment_start'] = intervals['fragment_start'].astype(int)
    intervals['fragment_end'] = intervals['fragment_end'].astype(int)

    intervals_by_accession = {
        accession: group.sort_values(['fragment_start', 'fragment_end', 'interpro_accession']).reset_index(drop=True)
        for accession, group in intervals.groupby('canonical_UniProtAC')
    }

    annotations = sites.apply(
        lambda row: pd.Series(
            classify_row(
                row,
                intervals_by_accession=intervals_by_accession,
                accession_col=args.accession_col,
                position_col=args.position_col,
                boundary_window=args.boundary_window,
            )
        ),
        axis=1,
    )
    annotated = pd.concat([sites, annotations], axis=1)
    annotated['domain_context_boundary_window'] = args.boundary_window
    annotated['domain_context_layer'] = 'InterPro_domain_context'

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    annotated_path = outdir / 'sites_with_domain_context.tsv'
    annotated.to_csv(annotated_path, sep='\t', index=False)

    summary = annotated['domain_context_class'].value_counts(dropna=False).rename_axis('domain_context_class').reset_index(name='site_count')
    summary['boundary_window'] = args.boundary_window
    summary.to_csv(outdir / 'domain_context_summary.tsv', sep='\t', index=False)
    print({
        'annotated_sites': len(annotated),
        'sites_with_domain_annotation': int(annotated['has_domain_annotation'].fillna(False).sum()),
        'output': str(annotated_path),
    })


if __name__ == '__main__':
    main()

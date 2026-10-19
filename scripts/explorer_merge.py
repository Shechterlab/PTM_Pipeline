"""Merge staged annotations into a Site Explorer export (the dict written to docs/data/methylarg_sites.json).

36_export_site_explorer_data.py calls these on the payload it builds. Each staging script (39 CPLM 4.0, 41 CARM1
knockout, 42 HPA subcellular locations) calls the same function through --explorer-json to update an existing export
in place, so both routes give the same file. The disorder intervals are the pipeline's own (PTM_MOBIDB_INTERVALS, read
by 06_annotate_disorder_context.py), so the website draws the regions its disorder classes were computed from. Fields added here are looked up by name (meta.protein_fields,
meta.site_fields), so their position does not matter.
"""
from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path

import pandas as pd
from scipy.stats import false_discovery_control, fisher_exact

from common import CPLM4_SOURCE_BITS, cplm4_evidence_mask

KME_GROUP = 'Lys methylation'
OTHER_PTM_SOURCES = ['Cross-PTM master table', 'CPLM 4.0 (Exp.)', 'CPLM 4.0 (Dat.)']
# HPA reliability scores, lowest first (index = code in the export).
HPA_RELIABILITY = ['Uncertain', 'Approved', 'Supported', 'Enhanced']


def field_index(data: dict, rows: str, name: str, empty) -> int:
    """Index of a named protein or site field, appended with `empty` on every row when the export lacks it."""
    names = data['meta']['protein_fields' if rows == 'proteins' else 'site_fields']
    if name not in names:
        names.append(name)
        for row in data[rows]:
            row.append(empty() if callable(empty) else empty)
    return names.index(name)


def update_export(path: str | Path, merge, *args) -> dict:
    """Load an export, apply one merge function and write it back; returns the merge function's summary."""
    path = Path(path)
    data = json.loads(path.read_text(encoding='utf-8'))
    summary = merge(data, *args)
    path.write_text(json.dumps(data, separators=(',', ':'), ensure_ascii=False), encoding='utf-8')
    return summary


def canonical_rows(data: dict):
    """(index, protein row) for proteins on their canonical accession (isoform entries carry their own numbering)."""
    for i, prot in enumerate(data['proteins']):
        if '-' not in prot[0]:
            yield i, prot


def merge_cplm4(data: dict, sites: pd.DataFrame, summary: dict) -> dict:
    """Lys methylation entries from 39_stage_cplm4_lys_methylation.py, replacing any CPLM entries already there."""
    meta = data['meta']
    field = meta['protein_fields'].index('other_ptms')
    group = meta['other_ptm_groups'].index(KME_GROUP)
    by_acc: dict[str, dict[int, int]] = defaultdict(dict)
    for acc, pos, evidence in sites[['canonical_UniProtAC', 'position', 'evidence']].itertuples(index=False):
        by_acc[acc][int(pos)] = by_acc[acc].get(int(pos), 0) | cplm4_evidence_mask(evidence)
    cplm_bits = sum(CPLM4_SOURCE_BITS.values())
    added = proteins_with = 0
    for prot in data['proteins']:
        entries: dict[tuple[int, int, str], int] = {}
        for item in prot[field]:
            mask = (item[3] if len(item) > 3 else 1) & ~cplm_bits
            if mask:
                entries[(item[0], item[1], item[2])] = mask
        hits = by_acc.get(prot[0], {}) if '-' not in prot[0] else {}
        for pos, mask in hits.items():
            key = (pos, group, 'K')
            entries[key] = entries.get(key, 0) | mask
        added += len(hits)
        proteins_with += bool(hits)
        prot[field] = [[pos, g, res, mask] for (pos, g, res), mask in sorted(entries.items())]
    meta['other_ptm_sources'] = OTHER_PTM_SOURCES
    meta['cplm4'] = {k: summary.get(k) for k in ('input_sha256', 'input_member')}
    meta['cplm4'].update({'validated_sites': len(sites), 'validated_proteins': int(sites['canonical_UniProtAC'].nunique()),
                          'explorer_sites': added, 'explorer_proteins': proteins_with})
    return {'explorer_sites': added, 'explorer_proteins': proteins_with}


def merge_carm1_ko(data: dict, sites: pd.DataFrame, summary: dict) -> dict:
    """CARM1 knockout response per site from 41_stage_carm1_ko_sites.py (site field carm1_ko).

    Each site gets [[cell line index, median log2 KO/WT, lowest p, peptides, peptides down, peptides up], ...] over the
    peptides that carry it; [] when no peptide was measured."""
    meta = data['meta']
    sites = sites.astype({'position': int, 'median_log2_ko_wt': float, 'min_p': float, 'peptides': int, 'peptides_down': int, 'peptides_up': int})
    cells = list(summary['cell_lines'])
    field = field_index(data, 'sites', 'carm1_ko', list)
    acc_of = {i: prot[0] for i, prot in canonical_rows(data)}
    site_at = {(acc_of[s[0]], s[1]): j for j, s in enumerate(data['sites']) if s[0] in acc_of}
    for s in data['sites']:
        s[field] = []
    matched = set()
    for row in sites.itertuples(index=False):
        j = site_at.get((row.canonical_UniProtAC, row.position))
        if j is None:
            continue
        matched.add(j)
        data['sites'][j][field].append([cells.index(row.cell_line), round(row.median_log2_ko_wt, 3), float(f'{row.min_p:.3g}'),
                                        row.peptides, row.peptides_down, row.peptides_up])
    for s in data['sites']:
        s[field].sort()
    placed = sites.drop_duplicates(['canonical_UniProtAC', 'position'])
    down = sites[sites['call'].eq('down')].drop_duplicates(['canonical_UniProtAC', 'position'])
    meta['carm1_ko'] = {
        'study': summary['study'], 'doi': summary['doi'], 'cell_lines': cells, 'rule': summary['rule'],
        'input_sha256': summary.get('input_sha256'),
        'placed_sites': len(placed),
        'placed_sites_down': len(down),
        'union_sites': len(matched),
        'union_sites_down': sum(1 for j in matched if any(c[4] for c in data['sites'][j][field])),
        'not_in_union': len(placed) - len(matched),
    }
    return {k: meta['carm1_ko'][k] for k in ('union_sites', 'union_sites_down', 'not_in_union')}


def merge_hpa(data: dict, locations: pd.DataFrame, summary: dict) -> dict:
    """HPA subcellular locations from 42_stage_hpa_subcellular.py (protein field hpa_locations) and, per location, a
    Fisher test of methylated against other HPA-profiled proteins (locations scored Approved or better)."""
    meta = data['meta']
    locations = locations.astype({'main': int})
    names = sorted(locations['location'].unique())
    field = field_index(data, 'proteins', 'hpa_locations', list)
    per_acc: dict[str, list] = defaultdict(list)
    for row in locations.itertuples(index=False):
        per_acc[row.canonical_UniProtAC].append([names.index(row.location), HPA_RELIABILITY.index(row.reliability), row.main])
    methylated = set()
    for _, prot in canonical_rows(data):
        prot[field] = sorted(per_acc.get(prot[0], []))
        methylated.add(prot[0])
    for prot in data['proteins']:
        if '-' in prot[0]:
            prot[field] = []
    confident = locations[locations['reliability'].ne('Uncertain')]
    universe = set(confident['canonical_UniProtAC'])
    meth_u = universe & methylated
    rows, pvals = [], []
    for k, name in enumerate(names):
        with_loc = set(confident.loc[confident['location'].eq(name), 'canonical_UniProtAC'])
        a = len(with_loc & meth_u)
        b = len(meth_u) - a
        c = len(with_loc - meth_u)
        d = len(universe) - len(meth_u) - c
        p = float(fisher_exact([[a, b], [c, d]]).pvalue)
        odds = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))   # Haldane-Anscombe corrected sample odds ratio
        rows.append([k, a, b, c, d, round(math.log2(odds), 3)])
        pvals.append(p)
    qvals = [float(q) for q in false_discovery_control(pvals, method='bh')] if pvals else []
    for row, p, q in zip(rows, pvals, qvals):
        row += [float(f'{p:.3g}'), float(f'{q:.3g}')]
    meta['hpa'] = {
        'version': summary.get('version', ''), 'input_sha256': summary.get('input_sha256', {}),
        'locations': names, 'reliability': HPA_RELIABILITY,
        'universe': len(universe), 'methylated': len(meth_u),
        'enrichment_fields': ['location', 'methylated with', 'methylated without', 'other with', 'other without', 'log2 OR', 'p', 'q'],
        'enrichment': rows,
    }
    return {'explorer_proteins_with_hpa': sum(1 for _, prot in canonical_rows(data) if prot[field]),
            'universe': len(universe), 'methylated_in_universe': len(meth_u)}


def merge_disorder_intervals(data: dict, intervals: pd.DataFrame, summary: dict) -> dict:
    """Disordered intervals per canonical protein (protein field disorder_intervals: [[start, end], ...])."""
    names = {'accession': 'canonical_UniProtAC', 'Entry': 'canonical_UniProtAC', 'start': 'fragment_start', 'end': 'fragment_end'}
    intervals = intervals.rename(columns={k: v for k, v in names.items() if k in intervals.columns and v not in intervals.columns})
    intervals = intervals.astype({'fragment_start': int, 'fragment_end': int})
    field = field_index(data, 'proteins', 'disorder_intervals', list)
    by_acc: dict[str, list] = defaultdict(list)
    for acc, a, b in intervals[['canonical_UniProtAC', 'fragment_start', 'fragment_end']].itertuples(index=False):
        by_acc[acc].append([a, b])
    for prot in data['proteins']:
        prot[field] = sorted(by_acc.get(prot[0], [])) if '-' not in prot[0] else []
    with_any = sum(1 for prot in data['proteins'] if prot[field])
    data['meta']['disorder_intervals'] = {'source': summary.get('source', 'MobiDB disorder intervals (PTM_MOBIDB_INTERVALS)'),
                                          'input_md5': summary.get('input_md5'), 'proteins': with_any}
    return {'explorer_proteins_with_intervals': with_any}

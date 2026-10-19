"""Merge staged annotations into a Site Explorer export (the dict written to docs/data/methylarg_sites.json).

36_export_site_explorer_data.py calls these on the payload it builds. Each staging script (39 CPLM 4.0, 41 CARM1
knockout, 42 HPA subcellular locations) calls the same function through --explorer-json to update an existing export
in place, so both routes give the same file. The disorder intervals are the pipeline's own (PTM_MOBIDB_INTERVALS, read
by 06_annotate_disorder_context.py), so the website draws the regions its disorder classes were computed from. The
disease-family membership comes from 34_regenerate_disease_ontology_figures.py, which calls the merge the same way. Fields added here are looked up by name (meta.protein_fields,
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


def label_symbols(labels: pd.DataFrame) -> dict[str, str]:
    """Accession -> the gene symbol the functional and disease analyses used for it (the first symbol of
    substrate_genename, else of genename, in protein_functional_labels.tsv; scripts/14 and 34 map these to Entrez)."""
    out = {}
    for acc, sub, gen in labels[['canonical_UniProtAC', 'substrate_genename', 'genename']].itertuples(index=False):
        sym = (str(sub or '').split(';', 1)[0].strip()) or (str(gen or '').split(';', 1)[0].strip())
        if sym and acc not in out:
            out[acc] = sym
    return out


def merge_disease_links(data: dict, evidence: dict, enriched: pd.DataFrame, labels: pd.DataFrame) -> dict:
    """Enriched DisGeNET terms per protein (evidence protein_disease_terms) and their families (protein field
    disease_family_mask), matched on the gene symbol the enrichment counted for that accession. Matching any symbol or
    alias linked proteins to terms of another gene whose symbol is among their aliases (GET3 to ARSA's terms).
    Isoform entries take their canonical accession's symbol. evidence['disease_terms'] must list `enriched` in order."""
    terms = evidence['disease_terms']
    if len(terms) != len(enriched) or any(t[0] != i for t, i in zip(terms, enriched['ID'])):
        raise ValueError('evidence disease_terms do not match the enrichment table')
    gene_terms: dict[str, list[int]] = defaultdict(list)
    for k, genes in enumerate(enriched['geneID']):
        for g in str(genes).split('/'):
            if g:
                gene_terms[g].append(k)
    symbols = label_symbols(labels)
    field = data['meta']['protein_fields'].index('disease_family_mask')
    links, changed = [], 0
    for prot in data['proteins']:
        ids = sorted(set(gene_terms.get(symbols.get(prot[0].split('-')[0], ''), [])))
        mask = 0
        for k in ids:
            if terms[k][2] >= 0:
                mask |= 1 << terms[k][2]
        changed += mask != prot[field]
        prot[field] = mask
        links.append(ids)
    before = sum(len(x) for x in evidence.get('protein_disease_terms', []))
    evidence['protein_disease_terms'] = links
    return {'term_links': sum(len(x) for x in links), 'term_links_before': before, 'proteins_with_changed_families': changed}


def merge_maron_peptide_calls(data: dict, peptides: pd.DataFrame, padj: float = 0.05) -> dict:
    """Append to each site's Maron Table S3 entries [IP, treatment, median change, lowest adj. p, peptide rows] the
    number of its peptides that fell and that rose at adj. p <= padj (maron_table_s3_sites.tsv). The site explorer calls
    a site decreased or increased only when its significant peptides agree; the median over all peptides can point the
    other way (one significant peptide up, several unchanged ones slightly down)."""
    meta = data['meta']
    sf = meta['site_fields']
    pk, pos_k, tk = sf.index('protein'), sf.index('position'), sf.index('maron_treatment')
    accs = [prot[0] for prot in data['proteins']]
    site_at = {(accs[row[pk]], row[pos_k]): j for j, row in enumerate(data['sites'])}
    ips, trs = meta['maron_ips'], meta['maron_treatments']
    counts: dict[tuple, list[int]] = defaultdict(lambda: [0, 0])
    for acc, pos, ip, tr, diff, q in peptides[['canonical_UniProtAC', 'position', 'ip', 'treatment', 'treatment_diff', 'treatment_padj']].itertuples(index=False):
        try:
            j, d, q = site_at.get((acc, int(pos))), float(diff), float(q)
        except (TypeError, ValueError):
            continue
        if j is None or ip not in ips or tr not in trs or math.isnan(d) or math.isnan(q) or q > padj or d == 0:
            continue
        counts[(j, ips.index(ip), trs.index(tr))][0 if d < 0 else 1] += 1
    disagree = 0
    for j, row in enumerate(data['sites']):
        for e in row[tk]:
            down, up = counts.get((j, e[0], e[1]), (0, 0))
            del e[5:]
            e += [down, up]
            if e[3] is not None and e[3] <= padj and e[2] is not None and e[2] != 0 and ((e[2] < 0 and not down) or (e[2] > 0 and not up)):
                disagree += 1
    meta['maron_treatment_fields'] = ['ip', 'treatment', 'median change', 'lowest adj. p', 'peptide rows',
                                      f'peptides down at adj. p <= {padj:g}', f'peptides up at adj. p <= {padj:g}']
    return {'entries_with_median_against_significant_peptides': disagree}


def merge_disease_families(data: dict, table: pd.DataFrame, summary: dict) -> dict:
    """Audited disease families whose gene unions hold each protein's gene, from 34_regenerate_disease_ontology_figures.py
    (the gene sets its family tests count): protein field disease_families, [index in meta.disease_families, number of
    the family's terms naming the gene] per family, and each family's term and gene counts in meta.disease_family_sets.
    Isoform entries take their canonical accession's families."""
    meta = data['meta']
    names = meta['disease_families']
    field = field_index(data, 'proteins', 'disease_families', list)
    per_acc: dict[str, dict[int, int]] = defaultdict(dict)
    unknown = set()
    for row in table.itertuples(index=False):
        if row.disease_family not in names:
            unknown.add(row.disease_family)
            continue
        k = names.index(row.disease_family)
        per_acc[row.canonical_UniProtAC][k] = max(per_acc[row.canonical_UniProtAC].get(k, 0), int(row.contributing_term_count))
    matched = 0
    for prot in data['proteins']:
        fams = per_acc.get(prot[0].split('-')[0], {})
        prot[field] = sorted([k, n] for k, n in fams.items())
        matched += bool(fams)
    info = {f['family']: f for f in summary.get('families', [])}
    meta['disease_family_sets'] = {
        'version': summary.get('version', ''), 'input_sha256': summary.get('input_sha256', {}),
        'fields': ['terms', 'genes', 'methylarginine_genes'],
        'families': [[int(info[n]['terms']), int(info[n]['genes']), int(info[n]['methylarginine_genes'])] if n in info else None
                     for n in names],
    }
    return {'proteins_with_a_family': matched, 'unknown_families': sorted(unknown)}


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


if __name__ == '__main__':
    # Re-apply the merges that read pipeline result tables to an existing export, without rebuilding it:
    #   python scripts/explorer_merge.py --results-root PTM_results --explorer-dir docs/data
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split('\n', 1)[0])
    ap.add_argument('--results-root', default='PTM_results')
    ap.add_argument('--explorer-dir', default='docs/data')
    args = ap.parse_args()
    root, out = Path(args.results_root), Path(args.explorer_dir)
    sites_path, ev_path = out / 'methylarg_sites.json', out / 'methylarg_evidence.json'
    sites, ev = json.loads(sites_path.read_text(encoding='utf-8')), json.loads(ev_path.read_text(encoding='utf-8'))
    read = lambda p: pd.read_csv(p, sep='\t', dtype=str, keep_default_na=False)
    disease = root / 'functional_class_union' / 'disease_enrichment'
    print('disease links:', merge_disease_links(sites, ev, read(disease / 'all_methyl_proteins_disease_enrichment.tsv'),
                                                read(root / 'functional_class_union' / 'protein_functional_labels.tsv')))
    print('Maron peptide calls:', merge_maron_peptide_calls(sites, read(root / 'maron_state_treatment' / 'maron_table_s3_sites.tsv')))
    for path, data in ((sites_path, sites), (ev_path, ev)):
        path.write_text(json.dumps(data, separators=(',', ':'), ensure_ascii=False), encoding='utf-8')

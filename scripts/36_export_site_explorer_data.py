#!/usr/bin/env python3
"""Encode a release (releases/VERSION/) as the Methylarginine Site Explorer's site and evidence data.

    python scripts/36_export_site_explorer_data.py --release releases/v1.1.0 --outdir docs/data

Writes docs/data/methylarg_sites.json (proteins, sites, vocabularies, overview tests) and methylarg_evidence.json (per-site
source records, disease terms), the format docs/assets/explorer.js reads: vocabularies, bit masks and arrays. Every value
comes from the release tables. Sites the release keeps pending (flag one_sided_window) are left out, as they are of
every analysis; the only other inputs are the MobiDB disorder intervals the protein maps draw
(--disorder-intervals, the ones the release's disorder classes were computed from) and the version notes of the
annotation layers in inputs/layers/. scripts/37 then adds the interaction data. --check decodes the files again and
compares them with the release tables.
"""
from __future__ import annotations

import argparse
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

SOURCES = ['Maron2021', 'ProMetheusDB', 'human_ptm_master']
SOURCE_OF_TABLE = {'maron2021_s5': 'Maron2021', 'prometheusdb_mmc2': 'ProMetheusDB', 'ptm_master': 'human_ptm_master'}
RAW_TABLE = {'maron2021_s5': 'maron_s5', 'prometheusdb_mmc2': 'prometheus_mmc2', 'ptm_master': 'human_ptm_master'}
STATES = ['', 'Rme1', 'Rme2', 'mixed']
STATE_SOURCES = [('methyl_states_reported', 'ProMetheusDB'), ('maron2021_s2_states', 'Maron 2021 Table S2'),
                 ('uniprot_states', 'UniProt'), ('uniprot_states_inferred', 'UniProt, by similarity (not counted)')]
DESIGN_LABEL = {'Rme1': 'Larsen 2016 (looked for Rme1 only)', 'Rme2': 'Wei 2020 (ADMA antibody, Rme2 only)'}
MOTIFS = ['RG', 'RGG', 'GRG', 'GAR', 'RXR', 'CARM1_proline_rich', 'CARM1_hydrophobic_proline', 'DR', 'RD', 'ER', 'RE']
MARON_IPS = ['MMA', 'ADMA', 'SDMA']
MARON_TREATMENTS = ['GSK591', 'MS023']
OTHER_PTM_GROUPS = ['Phosphorylation', 'Lys acylation', 'Ubiquitin/SUMO', 'Lys methylation']
PTM_GROUP = {'phospho': 0, 'lys_acyl': 1, 'ub_sumo': 2, 'lys_methyl': 3}
HPA_RELIABILITY = ['Uncertain', 'Approved', 'Supported', 'Enhanced']
CELL_LINES = ['MCF7', 'MDA-MB-231']
# The website's study names (docs/assets/explorer.js STUDY_DOI and the PhosphoSitePlus links key on them).
WEB_STUDY = {'PSP': 'PhosphoSitePlus'}
WINDOW_FLANK = 15
PADJ = 0.05


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, low_memory=False)


def L(value: str) -> list[str]:
    return value.split('; ') if value else []


def num(value, digits: int = 3):
    if value in ('', None):
        return None
    x = float(value)
    return None if math.isnan(x) else round(x, digits)


def sig(value, digits: int = 3):
    if value in ('', None):
        return None
    x = float(value)
    return None if math.isnan(x) else float(f'{x:.{digits}g}')


def integer(value):
    return None if value in ('', None) else int(float(value))


class Vocab:
    def __init__(self, values=None):
        self.values, self.lookup = [], {}
        for v in values or []:
            self.index(v)

    def index(self, value) -> int:
        if value not in self.lookup:
            self.lookup[value] = len(self.values)
            self.values.append(value)
        return self.lookup[value]

    def mask(self, labels) -> int:
        m = 0
        for x in labels:
            m |= 1 << self.index(x)
        return m


def overview_rows(stats: pd.DataFrame, analysis: str, extra=None) -> list[list]:
    rows = []
    for r in stats[stats['analysis'] == analysis].to_dict('records'):
        a, b, c, d = (int(r[k]) for k in 'abcd')
        row = [r['category'], a, round(a / (a + b), 4) if a + b else None, round(c / (c + d), 4) if c + d else None,
               num(r['log2_odds_ratio']), sig(r['q'])]
        rows.append(row + (extra(r) if extra else []))
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--release', required=True, help='the release folder, e.g. releases/v1.1.0')
    ap.add_argument('--inputs', default='inputs', help='frozen inputs (MobiDB intervals, layer notes)')
    ap.add_argument('--disorder-intervals', default='annotation/mobidb_disorder_intervals_2026-03-30.tsv', help='relative to --inputs')
    ap.add_argument('--panther-ontology', default='annotation/panther_19.0_protein_class_ontology.tsv', help='relative to --inputs')
    ap.add_argument('--panther-groups', default='config/panther_class_groups.json')
    ap.add_argument('--outdir', default='docs/data')
    ap.add_argument('--check', action='store_true', help='decode the written files and compare them with the release tables')
    args = ap.parse_args()
    rel, inputs, outdir = Path(args.release), Path(args.inputs), Path(args.outdir)
    info = json.loads((rel / 'release.json').read_text())
    v = info['version']
    t = lambda name: read(rel / name.format(v=v))  # noqa: E731

    sites = t('human_methylarginine_sites_{v}.tsv')
    proteins = t('human_methylarginine_proteins_{v}.tsv')
    studies = t('studies_{v}.tsv')
    stats = t('statistics_{v}.tsv')
    records = t('supporting/source_records_{v}.tsv')
    inhibitors = t('supporting/maron2021_inhibitors_{v}.tsv')
    carm1 = t('supporting/carm1_ko_{v}.tsv')
    disease_terms = t('supporting/disease_terms_{v}.tsv')
    disease_families = t('supporting/disease_families_{v}.tsv')
    density = t('supporting/density_windows_{v}.tsv')
    hpa = t('supporting/hpa_locations_{v}.tsv')
    other_ptms = t('supporting/other_ptm_sites_{v}.tsv')
    condensates = t('supporting/condensates_{v}.tsv')
    pending = sites['flags'].map(lambda f: 'one_sided_window' in L(f))
    n_pending = int(pending.sum())
    sites = sites[~pending].reset_index(drop=True)
    proteins = proteins[proteins['uniprot_accession'].isin(set(sites['uniprot_accession']))].reset_index(drop=True)

    # ---- Proteins ------------------------------------------------------------------------------------------------
    acc_list = list(proteins['uniprot_accession'])
    pidx = {a: i for i, a in enumerate(acc_list)}
    cond_vocab = Vocab(list(condensates['condensate_id']))
    families = list(disease_families['family'])
    family_idx = {f: k for k, f in enumerate(families)}
    term_family = [family_idx.get(f, -1) for f in disease_terms['family']]
    protein_terms = defaultdict(list)
    for k, accs in enumerate(disease_terms['uniprot_accessions']):
        for a in L(accs):
            protein_terms[a].append(k)
    dens = defaultdict(list)
    for r in density.to_dict('records'):
        dens[r['uniprot_accession']].append([int(r['window_aa']), int(r['start']), int(r['end']), int(r['sites']),
                                            num(r['log2_obs_over_null']), sig(r['q'])])
    ptms = defaultdict(list)
    for r in other_ptms.to_dict('records'):
        ptms[r['uniprot_accession']].append([int(r['position']), PTM_GROUP[r['ptm_type']], r['residue'], 1])
    hpa_names = sorted(set(hpa['location']) | set(stats.loc[stats['analysis'] == 'hpa_location', 'category']))
    hpa_of = defaultdict(list)
    for r in hpa.to_dict('records'):
        hpa_of[r['uniprot_accession']].append([hpa_names.index(r['location']), HPA_RELIABILITY.index(r['reliability']), int(r['main'] == 'true')])
    func_classes, func_subclasses = Vocab(), Vocab()
    protein_rows, protein_disease_terms = [], []
    for p in proteins.to_dict('records'):
        a = p['uniprot_accession']
        terms = sorted(protein_terms.get(a, []))
        fam_mask = 0
        for k in terms:
            if term_family[k] >= 0:
                fam_mask |= 1 << term_family[k]
        best = None
        if p['best_window_aa']:
            best = [int(p['best_window_aa']), int(p['best_window_start']), int(p['best_window_end']), int(p['best_window_sites']),
                    num(p['best_window_log2_obs_null']), sig(p['best_window_q'])]
        arch = None
        if p['architecture_class']:
            arch = [p['architecture_class'], [[min(c), max(c), len(c)] for c in
                                              ([int(x) for x in cl.split(',')] for cl in L(p['architecture_clusters']))]]
        local = None
        if p['nearest_spacing_min'] or p['max_sites_3aa']:
            local = [num(p['nearest_spacing_min'], 1), num(p['nearest_spacing_median'], 1)] + [integer(p[f'max_sites_{w}aa']) for w in (3, 5, 10, 20, 50)]
        gene = L(p['gene'])[0] if p['gene'] else ''
        aliases = sorted({g for g in L(p['gene_synonyms']) + L(p['reported_genes']) + L(p['gene'])[1:] if g != gene and ' ' not in g and ':' not in g})
        fams = dict(zip(L(p['disease_families']), map(int, L(p['disease_family_terms']))))
        protein_disease_terms.append(terms)
        protein_rows.append([
            a, gene, p['protein_name'], p['uniprot_entry_name'],
            func_classes.mask(L(p['functional_classes'])), func_subclasses.mask(L(p['functional_subclasses'])),
            int(p['protein_length']), int(p['arginine_count']), num(p['arg_enrichment_log2']), sig(p['arg_enrichment_q']),
            int('-' in a), sorted(cond_vocab.index(c) for c in L(p['condensates'])), fam_mask, int(p['als_ftd_family_gene'] == 'true'),
            local, best, dens.get(a, []), arch, sorted(ptms.get(a, [])), aliases, p['reference_sequence_md5'],
            sorted(hpa_of.get(a, [])), sorted([family_idx[f], n] for f, n in fams.items() if f in family_idx),
        ])

    # ---- Sites ----------------------------------------------------------------------------------------------------
    site_index = {s: j for j, s in enumerate(sites['site_id'])}
    site_of = records['placed_site'].where(records['placed_site'] != '', records['pending_site'])
    records = records.assign(_site=site_of)
    refs, studies_vocab = Vocab(), Vocab()
    web_name = {r['study_id']: WEB_STUDY.get(r['study_id'], r['label']) for r in studies.to_dict('records')}
    study_citation = {}
    for r in records.to_dict('records'):
        if r['source_table'] == 'maron2021_s5' and r['reference_text']:
            study_citation.setdefault(web_name[r['study']], r['reference_text'])
    by_site = defaultdict(list)
    for i, s in zip(records.index, records['_site']):
        if s:
            by_site[s].append(i)
    rec = records.to_dict('records')
    tr = defaultdict(list)
    for r in inhibitors.to_dict('records'):
        tr[r['site_id']].append([MARON_IPS.index(r['antibody']), MARON_TREATMENTS.index(r['inhibitor']), num(r['log2_diff_median']),
                                 sig(r['adj_p_min']), int(r['n_peptides']), int(r['n_decreased']), int(r['n_increased'])])
    ko = defaultdict(list)
    for r in carm1.to_dict('records'):
        ko[r['site_id']].append([CELL_LINES.index(r['cell_line']), num(r['log2_ko_over_wt']), sig(r['p']), int(r['n_peptides']),
                                 int(r['n_decreased']), int(r['n_increased'])])
    disorder_classes, domain_contexts, domain_classes = Vocab(), Vocab(), Vocab()
    interpro, interpro_idx = [], {}
    window_flags = {'missing': 0, 'center_not_r': 0}
    site_rows, evidence = [], []
    for s in sites.to_dict('records'):
        sid = s['site_id']
        state = {'': 0, 'Rme1': 1, 'Rme2': 2, 'Rme1_Rme2': 3}[s['methyl_state']]
        # Where the degree comes from: [source, degrees] for each source that gives one.
        state_detail = [[label, s[col]] for col, label in STATE_SOURCES if s[col]]
        state_detail += [[DESIGN_LABEL[x], x] for x in L(s['states_by_study_design'])]
        ipr = -1
        if s['nearest_domain_interpro']:
            if s['nearest_domain_interpro'] not in interpro_idx:
                interpro_idx[s['nearest_domain_interpro']] = len(interpro)
                interpro.append([s['nearest_domain_interpro'], s['nearest_domain_name']])
            ipr = interpro_idx[s['nearest_domain_interpro']]
        window = s['window_31']
        if window[WINDOW_FLANK] != 'R':
            window_flags['center_not_r'] += 1
        motif_mask = sum(1 << MOTIFS.index(m) for m in L(s['motifs']))
        ids = sorted(by_site[sid], key=lambda i: int(rec[i]['record_id']))
        site_refs = []
        for i in ids:
            ref = rec[i]['reference_text']
            if ref and refs.index(ref) not in site_refs:
                site_refs.append(refs.index(ref))
        evidence.append([[
            0, refs.index(rec[i]['reference_text']) if rec[i]['reference_text'] else -1, rec[i]['pmid'], rec[i]['evidence_sequence'],
            rec[i]['reported_enzyme'], rec[i]['reported_methyl_state'], rec[i]['reported_accession'], rec[i]['reported_position'],
            rec[i]['placement_code'], web_name[rec[i]['study']]] for i in ids])
        site_rows.append([
            pidx[s['uniprot_accession']], int(s['position']), state,
            sum(1 << SOURCES.index(SOURCE_OF_TABLE[x]) for x in L(s['source_tables'])),
            site_refs, disorder_classes.index(s['disorder']), integer(s['disorder_edge_distance']), domain_contexts.index(s['domain']),
            ipr, domain_classes.index(s['nearest_domain_class']), integer(s['domain_edge_distance']), s['domain_side'],
            window, motif_mask, int(window[WINDOW_FLANK] != 'R'), ';'.join(L(s['reported_positions'])),
            sum(1 << k for k, st in enumerate(('Rme1', 'Rme2')) if st in L(s['maron2021_s2_states'])),
            sum(1 << MARON_IPS.index(x) for x in L(s['maron2021_s3_antibodies'])),
            sorted(tr.get(sid, [])), studies_vocab.mask(web_name[x] for x in L(s['studies'])), sorted(ko.get(sid, [])),
            state_detail or None,
        ])
    # evidence rows: table and study codes as vocabularies
    source_tables, remap_statuses = Vocab(), Vocab()
    for rows, i_s in zip(evidence, [sorted(by_site[sid], key=lambda i: int(rec[i]['record_id'])) for sid in sites['site_id']]):
        for row, i in zip(rows, i_s):
            row[0] = source_tables.index(RAW_TABLE[rec[i]['source_table']])
            row[8] = remap_statuses.index(row[8])
            row[9] = studies_vocab.index(row[9])
    study_counts = Counter()
    study_unique = Counter()
    for row in site_rows:
        bits = [k for k in range(len(studies_vocab.values)) if row[19] >> k & 1]
        for k in bits:
            study_counts[k] += 1
        if len(bits) == 1:
            study_unique[bits[0]] += 1

    # ---- Overview, condensates, layers ----------------------------------------------------------------------------
    parent = defaultdict(Counter)
    for p in proteins.to_dict('records'):
        for sub in L(p['functional_subclasses']):
            for broad in L(p['functional_classes']):
                parent[sub][broad] += 1
    fam_stats = stats[stats['analysis'].isin(['disease_family', 'disease_family_additional'])].set_index('category')
    treatments = []
    for drug in MARON_TREATMENTS:
        for ip in ('SDMA', 'ADMA', 'MMA'):
            g = inhibitors[(inhibitors['inhibitor'] == drug) & (inhibitors['antibody'] == ip)]
            treatments.append([drug, ip, len(g), int((g['call'] == 'decreased').sum()), int((g['call'] == 'increased').sum()),
                               num(statistics.median(g['log2_diff_median'].astype(float))) if len(g) else None])
    cond_stats = stats[stats['analysis'] == 'condensate'].set_index('category')
    cond_rows = [[r['condensate_name'], r['condensate_id'], r['condensate_type'], int(r['n_members']), int(r['n_members_methylated']),
                  num(cond_stats.at[r['condensate_id'], 'log2_odds_ratio']) if r['condensate_id'] in cond_stats.index else None,
                  sig(cond_stats.at[r['condensate_id'], 'q']) if r['condensate_id'] in cond_stats.index else None]
                 for r in condensates.to_dict('records')]
    methylated = set(acc_list)
    in_condensate = {a for a, row in zip(acc_list, protein_rows) if row[11]}
    universe = int(cond_stats['a'].astype(int).iloc[0] + cond_stats['b'].astype(int).iloc[0] + cond_stats['c'].astype(int).iloc[0]
                   + cond_stats['d'].astype(int).iloc[0]) if len(cond_stats) else 0
    hpa_stats = stats[stats['analysis'] == 'hpa_location'].set_index('category')
    hpa_universe = int(hpa_stats.iloc[0][['a', 'b', 'c', 'd']].astype(int).sum()) if len(hpa_stats) else 0
    hpa_methyl = int(hpa_stats.iloc[0][['a', 'b']].astype(int).sum()) if len(hpa_stats) else 0
    layer = lambda name: json.loads((inputs / 'layers' / name).read_text()) if (inputs / 'layers' / name).exists() else {}  # noqa: E731
    hpa_note, ko_note = layer('hpa_summary.json'), layer('carm1_ko_summary.json')
    ko_sites = set(carm1['site_id'])
    dis = read(inputs / args.disorder_intervals)
    intervals = defaultdict(list)
    for a, s0, s1 in zip(dis['canonical_UniProtAC'], dis['fragment_start'], dis['fragment_end']):
        if a in pidx:
            intervals[a].append([int(float(s0)), int(float(s1))])
    for a, row in zip(acc_list, protein_rows):
        row.append(sorted(intervals.get(a, [])))

    # PANTHER protein classes: [group, top-level class (-1: none), the tested subclasses the protein is in]. The
    # release gives each protein's most specific classes; a protein is in every class above them too.
    ontology = read(inputs / args.panther_ontology).set_index('class_id')
    panther_config = json.loads(Path(args.panther_groups).read_text())
    by_name = {n: c for c, n in ontology['name'].items()}

    def above(name: str) -> list[str]:
        out, c = [], by_name[name]
        while c in ontology.index and c != 'PC00000':
            out.append(ontology.at[c, 'name'])
            c = ontology.at[c, 'parent_id']
        return out
    pg = [r['category'] for r in stats[stats['analysis'] == 'panther_group'].to_dict('records')]
    pc = [r['category'] for r in stats[stats['analysis'] == 'panther_class'].to_dict('records')]
    ps = [r['category'] for r in stats[stats['analysis'] == 'panther_subclass'].to_dict('records')]
    for p, row in zip(proteins.to_dict('records'), protein_rows):
        inside = {x for leaf in L(p['panther_classes']) for x in above(leaf)}
        row.append([pg.index(p['panther_group']), pc.index(p['panther_class']) if p['panther_class'] else -1,
                    sorted(k for k, x in enumerate(ps) if x in inside)])
    group_of = {ontology.at[c, 'name']: g['name'] for g in panther_config['groups'] for c in g['classes']}
    top_of = {n: ontology.at[ontology.at[by_name[n], 'top_level_id'], 'name'] for n in ps}
    totals = lambda r: [int(r['a']) + int(r['c'])]  # noqa: E731

    payload = {
        'meta': {
            'title': 'Human methylarginine site union',
            'generated': info['release_date'],
            'provenance': {
                'dataset_version': v,
                'pipeline_commit': info['code_commit'],
                'uncommitted_changes': False,
                'results_root': f'releases/{v}',
                'reference_sequences': {'file': info['reference']['fasta'], 'md5': info['reference']['fasta_md5'],
                                        'sequences': info['reference']['entries'], 'uniprot_release': info['reference']['release'],
                                        'uniprot_release_date': info['reference']['release_date']},
                'inputs': [{'file': x['file'].split('/')[-1], 'md5': x['md5']} for x in info['inputs'] if x['file'].startswith('inputs/sources/')],
                'validation': info.get('validation', {}),
                'changes': {k: info['counts'].get(k) for k in ('previous_sites', 'kept', 'kept_pending', 'removed', 'added')},
                'coordinates': 'Positions are 1-based on the reference sequences (canonical UniProt entries, initiator methionine = 1). '
                               'Protein field reference_md5 is the checksum of each reference sequence in UniProt\'s format; current UniProt '
                               'annotations line up with the positions only where the current sequence has the same checksum.',
            },
            'n_sites': len(site_rows), 'n_proteins': len(protein_rows), 'n_pending_not_shown': n_pending,
            'window_flank': WINDOW_FLANK, 'window_flags': window_flags,
            'sources': SOURCES,
            'studies': [[name, study_citation.get(name, ''), study_counts[k], study_unique[k]] for k, name in enumerate(studies_vocab.values)],
            'states': STATES, 'motifs': MOTIFS, 'maron_ips': MARON_IPS, 'maron_treatments': MARON_TREATMENTS,
            'other_ptm_groups': OTHER_PTM_GROUPS, 'other_ptm_sources': ['Cross-PTM master table'],
            'func_classes': func_classes.values, 'func_subclasses': func_subclasses.values,
            'disorder_classes': disorder_classes.values, 'domain_contexts': domain_contexts.values, 'domain_classes': domain_classes.values,
            'disease_families': families,
            # the eight audited families of the paper come first; the others are rule-based additions (scripts/disease_families.py)
            'disease_family_audited': [r['audited'] == 'true' for r in disease_families.to_dict('records')],
            'condensates': cond_rows,
            'condensate_fields': ['name', 'cdcode_id', 'type', 'members', 'methyl_members', 'log2_odds_ratio', 'q_value'],
            'condensate_summary': {'methyl_proteins_total': str(len(methylated)), 'condensate_universe_proteins': str(universe),
                                   'methyl_proteins_in_any_condensate': str(len(in_condensate)),
                                   'condensate_count_after_filters': str(len(cond_stats))},
            'refs': refs.values, 'interpro': interpro,
            'overview': {
                'functional_classes': overview_rows(stats, 'functional_class'),
                'motifs': overview_rows(stats, 'motif'),
                'disorder': overview_rows(stats, 'disorder'),
                'domain_context': overview_rows(stats, 'domain_context'),
                'domain_classes': overview_rows(stats, 'domain_class'),
                'disease_families': [[f, int(fam_stats.at[f, 'a']), int(fam_stats.at[f, 'a']) + int(fam_stats.at[f, 'c']),
                                      num(fam_stats.at[f, 'odds_ratio']), sig(fam_stats.at[f, 'q'])] for f in families],
                'treatments': treatments,
                # [category, methylated, fraction of methylated, fraction of others, log2 OR, q, (parent,) proteins]
                'panther_groups': overview_rows(stats, 'panther_group', totals),
                'panther_classes': overview_rows(stats, 'panther_class', lambda r: [group_of.get(r['category'], ''), *totals(r)]),
                'panther_subclasses': overview_rows(stats, 'panther_subclass', lambda r: [top_of[r['category']], *totals(r)]),
                'functional_subclasses': overview_rows(stats, 'functional_subclass',
                                                       lambda r: [parent[r['category']].most_common(1)[0][0] if parent[r['category']] else '']),
            },
            'protein_fields': ['accession', 'gene', 'protein_name', 'uniprot_id', 'func_class_mask', 'func_subclass_mask', 'length',
                               'arginine_count', 'shrinkage_log2_enrichment', 'arg_odds_q', 'is_isoform', 'condensates',
                               'disease_family_mask', 'als_ftd_audited_gene', 'local_density', 'best_density_window',
                               'significant_density_windows', 'cluster_architecture', 'other_ptms', 'gene_aliases', 'reference_md5',
                               'hpa_locations', 'disease_families', 'disorder_intervals', 'panther'],
            'panther': {'version': panther_config['panther_version'], 'unclassified': panther_config['unclassified'],
                        'groups': pg, 'classes': pc, 'subclasses': ps,
                        'universe': int(stats[stats['analysis'] == 'panther_group'][['a', 'c']].astype(int).sum().sum())},
            'site_fields': ['protein', 'position', 'state', 'source_mask', 'refs', 'disorder_class', 'disorder_edge_distance',
                            'domain_context', 'interpro', 'domain_class', 'domain_edge_distance', 'domain_side', 'window_31',
                            'motif_mask', 'window_center_mismatch', 'source_positions_original', 'maron_s2_state_mask',
                            'maron_ip_mask', 'maron_treatment', 'study_mask', 'carm1_ko', 'state_detail'],
            'maron_treatment_fields': ['ip', 'treatment', 'median change', 'lowest adj. p', 'peptide rows',
                                       f'peptides down at adj. p <= {PADJ:g}', f'peptides up at adj. p <= {PADJ:g}'],
            'carm1_ko': {'study': ko_note.get('study', 'Shishkova 2017'), 'doi': ko_note.get('doi', '10.1038/ncomms15571'),
                         'cell_lines': CELL_LINES, 'rule': ko_note.get('rule', ''), 'input_sha256': ko_note.get('input_sha256'),
                         'placed_sites': ko_note.get('placed_sites'), 'placed_sites_down': ko_note.get('placed_sites_down'),
                         'union_sites': len(ko_sites),
                         'union_sites_down': int(carm1.groupby('site_id')['n_decreased'].agg(lambda x: x.astype(int).sum() > 0).sum()),
                         'not_in_union': (ko_note.get('placed_sites') or 0) - len(ko_sites)},
            'hpa': {'version': hpa_note.get('version', ''), 'input_sha256': hpa_note.get('input_sha256', {}), 'locations': hpa_names,
                    'reliability': HPA_RELIABILITY, 'universe': hpa_universe, 'methylated': hpa_methyl,
                    'enrichment_fields': ['location', 'methylated with', 'methylated without', 'other with', 'other without', 'log2 OR', 'p', 'q'],
                    'enrichment': [[hpa_names.index(n), int(r['a']), int(r['b']), int(r['c']), int(r['d']), num(r['log2_odds_ratio']),
                                    sig(r['p']), sig(r['q'])] for n, r in hpa_stats.iterrows()]},
            'disease_family_sets': {'version': 'audited_v1', 'fields': ['terms', 'genes', 'methylarginine_genes'],
                                    'families': [[int(r['n_terms']), int(r['n_genes']), int(r['n_methylated_genes'])] for r in disease_families.to_dict('records')]},
            'disorder_intervals': {'source': f'MobiDB disorder intervals ({args.disorder_intervals})', 'proteins': sum(1 for r in protein_rows if r[-2])},
        },
        'proteins': protein_rows,
        'sites': site_rows,
    }
    evidence_payload = {
        'source_tables': source_tables.values, 'studies': studies_vocab.values, 'remap_statuses': remap_statuses.values,
        'evidence_fields': ['source_table', 'reference', 'pmid', 'peptide', 'enzyme', 'state', 'source_accession', 'source_position',
                            'remap_status', 'study'],
        'refs': refs.values, 'sites': evidence,
        'disease_terms': [[r['term_id'], r['term_name'], term_family[k], int(r['n_methylated']), sig(r['p_adjust'])]
                          for k, r in enumerate(disease_terms.to_dict('records'))],
        'protein_disease_terms': protein_disease_terms,
    }
    outdir.mkdir(parents=True, exist_ok=True)
    for name, data in (('methylarg_sites.json', payload), ('methylarg_evidence.json', evidence_payload)):
        path = outdir / name
        path.write_text(json.dumps(data, separators=(',', ':'), ensure_ascii=False), encoding='utf-8')
        print(f'wrote {path} ({path.stat().st_size / 1e6:.2f} MB)')
    print(f'{len(site_rows):,} sites, {len(protein_rows):,} proteins, {len(refs.values)} references, {len(evidence_payload["disease_terms"])} disease terms')
    if args.check:
        check(outdir, sites, proteins, records, disease_terms, web_name)


def check(outdir: Path, sites: pd.DataFrame, proteins: pd.DataFrame, records: pd.DataFrame, disease_terms: pd.DataFrame, web_name: dict) -> None:
    """Decode the written files and compare them with the release tables: every site, protein, record and term."""
    d = json.loads((outdir / 'methylarg_sites.json').read_text())
    e = json.loads((outdir / 'methylarg_evidence.json').read_text())
    m = d['meta']
    P = d['proteins']
    bad = []
    if [p[0] for p in P] != list(proteins['uniprot_accession']):
        bad.append('protein order')
    for s, r in zip(d['sites'], sites.to_dict('records')):
        if f'{P[s[0]][0]}:R{s[1]}' != r['site_id'] or s[12] != r['window_31'] \
                or [m['motifs'][k] for k in range(len(m['motifs'])) if s[13] >> k & 1] != L(r['motifs']) \
                or m['disorder_classes'][s[5]] != r['disorder'] or m['domain_contexts'][s[7]] != r['domain'] \
                or m['domain_classes'][s[9]] != r['nearest_domain_class'] \
                or {m['studies'][k][0] for k in range(len(m['studies'])) if s[19] >> k & 1} != {web_name[x] for x in L(r['studies'])}:
            bad.append(r['site_id'])
    n_records = sum(len(x) for x in e['sites'])
    shown = set(sites['site_id'])
    expected = int(records['placed_site'].isin(shown).sum() + records['pending_site'].isin(shown).sum())
    if n_records != expected:
        bad.append(f'evidence rows {n_records} vs {expected}')
    if len(d['sites']) != len(sites) or len(e['disease_terms']) != len(disease_terms):
        bad.append('counts')
    for p, r in zip(P, proteins.to_dict('records')):
        if p[6] != int(r['protein_length']) or p[20] != r['reference_sequence_md5']:
            bad.append(r['uniprot_accession'])
    print('check: ' + ('the site data match the release tables' if not bad else f'{len(bad)} differences, e.g. {bad[:5]}'))
    if bad:
        raise SystemExit(1)


if __name__ == '__main__':
    main()

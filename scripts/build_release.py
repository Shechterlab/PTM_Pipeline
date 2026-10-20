#!/usr/bin/env python3
"""Assemble a release of the human methylarginine sites from a pipeline run (hpc/build_release.sh runs it last).

Reads
  --results-root   the run: remapped/ (every source record and where it was placed), integrated/ (the union and the
                   pending sites), qc/remap_evidence_check/ (the independent check of every record)
  --previous       the previous release folder (releases/v1.0.0)
  --inputs         the frozen inputs (checked against inputs/MANIFEST.tsv by hpc/build_release.sh)
and writes --outdir (releases/VERSION), for version V and previous version P:

  human_methylarginine_sites_V.tsv      one row per site, including the sites kept pending (flag one_sided_window)
  human_methylarginine_sites_V.xlsx     the same table as a spreadsheet
  studies_V.tsv                         one row per study or database behind the sites
  supporting/source_records_V.tsv       one row per source record: what the study reported, where it is placed
  supporting/site_changes_P_to_V.tsv    one row per site of either version: kept, removed or added, and why
  supporting/maron2021_detection_V.tsv  Maron et al. 2021 Tables S2 and S3: which sites each table detects
  supporting/maron2021_inhibitors_V.tsv Maron et al. 2021 Table S3: response to the PRMT inhibitors
  supporting/carm1_ko_V.tsv             response to CARM1 knockout (Shishkova et al. 2017)
  human_methylarginine_proteins_V.tsv   one row per protein: identity, function, enrichment, clustering, context
  statistics_V.tsv                      every enrichment test the website or the paper reports, with its 2x2 counts
  supporting/disease_terms_V.tsv, density_windows_V.tsv, hpa_locations_V.tsv, other_ptm_sites_V.tsv
  data_dictionary.tsv                   every column of every table, with its codes
  release.json                          version, counts, inputs and files
  README.md, CHANGES.md, MD5SUMS

Codes and column definitions are in scripts/release_codes.py. Structural context and motifs use the pipeline's own
classifiers (scripts/05, 06, 07, 08), so they match the analyses.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))

import release_codes as rc  # noqa: E402
import release_docs  # noqa: E402
import release_tables  # noqa: E402
from common import classify_by_patterns  # noqa: E402


def load_script(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


remap = load_script('remap', '02_remap_arg_methyl_to_canonical.py')
domain_context = load_script('domain_context', '05_annotate_domain_context.py')
disorder_context = load_script('disorder_context', '06_annotate_disorder_context.py')
motif_odds = load_script('motif_odds', '08_motif_and_arg_odds.py')

# Frozen inputs read here, relative to --inputs.
INPUT = {
    'fasta': 'reference/uniprot_human_reviewed_2026_03.fasta',
    'entries': 'reference/uniprot_human_reviewed_2026_03_entries.tsv',
    'uniprot_release': 'reference/uniprot_release.json',
    'disorder': 'annotation/mobidb_disorder_intervals_2026-03-30.tsv',
    'domains': 'annotation/interpro_108.0_domain_intervals.tsv',
    'lift': 'annotation/interval_lift_to_uniprot_2026_03.tsv',
    'maron_s2': 'sources/maron2021_tableS2_methylarginine_ip_sites.xlsx',
    'maron_s3': 'sources/maron2021_tableS3_methylpeptides_prmt_inhibitors.xlsx',
    'carm1': 'layers/carm1_ko_sites.tsv',
    'uniprot_states': 'annotation/uniprot_2026_03_methylarginine_mod_res.tsv',
}
DOMAIN_ONTOLOGY = SCRIPTS.parent / 'config' / 'domain_class_ontology.json'


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, low_memory=False)


def clean(value) -> str:
    value = str(value if value is not None else '').strip()
    return '' if value.lower() in ('nan', 'none', '<na>') else value


def as_list(value: str) -> str:
    """A source field that holds several tokens separated by ';', as a release list."""
    return '; '.join(t.strip() for t in clean(value).split(';') if t.strip())


def key(acc: str, pos) -> str:
    return f'{acc}:R{int(float(pos))}'


def split_key(site: str) -> tuple[str, int]:
    acc, pos = site.split(':R')
    return acc, int(pos)


def join(values) -> str:
    return '; '.join(str(v) for v in values)


def ordered(values, order: list[str]) -> list[str]:
    values = set(values)
    unknown = values - set(order)
    if unknown:
        raise SystemExit(f'codes missing from scripts/release_codes.py: {sorted(unknown)}')
    return [v for v in order if v in values]


def study_of(row) -> str:
    """study_id of a source record."""
    if row['source_table'] == 'human_ptm_master':
        return {'unip': 'UniProt', 'iedb': 'IEDB'}[row['source']]
    if row['source_table'] == 'prometheus_mmc2':
        return 'ProMetheusDB'
    ref = row['reference'].strip()
    if ref.lower().startswith('this study'):
        return 'Maron2021'
    if 'PhosphoSitePlus' in ref:
        return 'PSP'
    author = ref.split(' et al.')[0].split()[-1]
    return author + re.findall(r'\((\d{4})\)', ref)[-1]


def window_of(row: pd.Series) -> str:
    """The sequence window a record carries besides its marked peptide (ProMetheusDB SeqWindow)."""
    return as_list(row.get('seq_window', ''))


def evidence_of(row: pd.Series) -> tuple[str, str]:
    """Kind and text of a record's sequence evidence, as the source gives it."""
    if clean(row.get('peptides', '')):
        return 'marked_peptide', as_list(row['peptides'])
    for field in ('seq_window', 'peptide_sequence'):
        text = as_list(row.get(field, ''))
        if text:
            length = int(float(row.get(f'{field}_window_length', 0) or 0))
            return (f'window_{length}' if length else 'peptide'), text
    return 'position_only', ''


# ---- The previous version, in this version's numbering ---------------------------------------------------------------

def reference_lift(inputs: Path) -> tuple[dict[str, int], set[str]]:
    """Residue shift per accession whose sequence gained N-terminal residues since the previous reference, and the
    accessions whose numbering cannot be carried over (scripts/00_lift_intervals_to_reference.py report)."""
    report = read_tsv(inputs / INPUT['lift'])
    shift = {r.accession: int(r.shift) for r in report.itertuples() if r.change == 'n_terminal_extension'}
    broken = set(report.loc[report['change'].isin(['other_change', 'not_in_new_reference']), 'accession'])
    return shift, broken


def lifted(site: str, shift: dict[str, int]) -> str:
    if not site:
        return ''
    acc, pos = split_key(site)
    return key(acc, pos + shift.get(acc, 0))


# ---- Records ---------------------------------------------------------------------------------------------------------

def record_reason(row: pd.Series, old: str, new: str, reference: dict, typed: pd.Series) -> str:
    """Reason code for a record whose placement differs between the previous version and this one (keys are lifted)."""
    status = row['remap_status']
    if status == 'sequence_peptide_start':
        return 'peptide_start'
    evidence = remap.row_evidence(typed)
    windows = [e for e in evidence if e['anchor_type'] in ('centered_window', 'truncated_window')]
    if new and not old:
        if any(re.search(r'^[_X-]|[_X-]$', t.strip()) for f in ('seq_window', 'peptide_sequence') for t in clean(row.get(f, '')).split(';') if t.strip()):
            return 'window_padding'
        # The previous code searched the row's accession and gene symbol; a record it could not place that now lands
        # on another entry was reached through the accession history, a gene synonym or the proteome-wide search.
        if status == 'sequence_unique_proteome':
            return 'entry_not_in_reference'
        reported = remap.canonicalize_uniprot_accession(row['substrate_UniProtAC_original'])
        if clean(row.get('source_accession_history', '')) or split_key(new)[0] != reported:
            return 'outdated_accession_or_symbol'
        return 'other'
    if old and not new:
        if any(e['offsets'] and len(e['sequence']) > e['offsets'][0] >= 0 and e['sequence'][e['offsets'][0]] == 'K' and e['anchor_type'] == 'centered_window' for e in windows):
            return 'lysine_window'
        # Evidence that is nowhere in the entry the previous version used (a source accession of another protein).
        old_seq = reference.get(split_key(old)[0], '')
        if evidence and old_seq and not any(remap.find_all_occurrences(old_seq, e['sequence']) for e in evidence):
            return 'evidence_mismatch'
        peptides = [e for e in evidence if e['anchor_type'] == 'peptide']
        if peptides and not windows and all(len(e['offsets']) > 1 for e in peptides) and status in ('unresolved_no_anchor', 'unresolved_evidence_elsewhere'):
            return 'multi_arginine_peptide'
        if status in ('unresolved_evidence_mismatch', 'unresolved_evidence_elsewhere'):
            return 'evidence_mismatch'
        if status == 'unresolved_ambiguous_hit':
            return 'ambiguous_evidence'
        return 'other'
    old_acc, new_acc = split_key(old)[0], split_key(new)[0]
    if old_acc != new_acc:
        return 'unreviewed_accession' if old_acc not in reference else 'wrong_protein'
    if status == 'sequence_position_verified':
        return 'peptide_recentred'
    if row['remap_anchor_type'] == 'truncated_window':
        return 'truncated_window'
    if not windows and any(e['anchor_type'] == 'peptide' for e in evidence):
        return 'peptide_recentred'
    return 'other'


def build_records(rows, typed, check, allowed, previous, pending_of, shift, reference) -> pd.DataFrame:
    out = []
    for i, r in rows.iterrows():
        study = study_of(r)
        kind, text = evidence_of(typed.loc[i])
        new = key(r['corrected_canonical_UniProtAC'], r['corrected_position']) if r['is_position_resolved'] == 'True' else ''
        old_published = previous.at[i, 'placed_site'] if i < len(previous) else ''
        old = lifted(old_published, shift)
        if i >= len(previous):
            change = 'new_record'
        elif new == old:
            change = 'same' if new else 'not_placed'
        else:
            change = 'moved' if new and old else 'newly_placed' if new else 'no_longer_placed'
        if change == 'same' and old != old_published:
            reason = 'reference_renumbered'
        elif change in ('same', 'not_placed', 'new_record'):
            reason = ''
        else:
            reason = record_reason(r, old, new, reference, typed.loc[i])
        c = check.loc[i]
        out.append({
            'record_id': i + 1,
            'source_table': rc.SOURCE_TABLE_CODE[r['source_table']],
            'study': study,
            'study_doi': dict((s[0], s[3]) for s in rc.STUDIES)[study],
            'reference_text': clean(r['reference']),
            'pmid': clean(r['pmid']),
            'reported_accession': as_list(r['substrate_UniProtAC_original']),
            'reported_gene': clean(r['substrate_genename']),
            'reported_position': clean(r['position']),
            'reported_methyl_state': clean(r['arg_methyl_state']),
            'reported_enzyme': clean(r.get('enzyme_genename', '')),
            'evidence_kind': kind,
            'evidence_sequence': text,
            'evidence_window': window_of(typed.loc[i]) if kind == 'marked_peptide' else '',
            'source_correction': clean(r.get('source_correction', '')),
            'evidence_as_reported': clean(r.get('peptide_sequence_reported', '')),
            'placed_site': new,
            'placement_code': r['remap_status'],
            'evidence_check': rc.VERDICT_CODE[c['verdict']],
            'accepted_exception': allowed.get((r['source_table'], c['reported_accession'], c['reported_position'], c['placed_site']), ''),
            'pending_site': pending_of.get(i + 1, ''),
            'previous_site': old_published,
            '_previous_lifted': old,
            'change_since_previous_version': change,
            'reason_code': reason,
        })
    return pd.DataFrame(out)


# ---- Maron et al. 2021 Tables S2 and S3, CARM1 knockout ------------------------------------------------------------

def place_maron_peptide(reference: dict, accession: str, peptide: str, begin, position: int) -> str:
    """Site key for a Maron S2/S3 site, or '' when its peptide does not support it on the reference entry: the peptide
    must sit at the reported start, or occur once in the entry (the offset then carries over), and cover an arginine."""
    acc = remap.canonicalize_uniprot_accession(remap.normalize_accession(str(accession).strip()))
    seq, peptide = reference.get(acc, ''), re.sub(r'[^A-Z]', '', str(peptide).upper())
    begin = pd.to_numeric(begin, errors='coerce')
    if not seq or not peptide or pd.isna(begin):
        return ''
    begin, offset = int(begin), position - int(begin)
    if not 0 <= offset < len(peptide) or peptide[offset] != 'R':
        return ''
    if seq[begin - 1:begin - 1 + len(peptide)] != peptide:
        hits = [m.start() for m in re.finditer(f'(?={peptide})', seq)]
        if len(hits) != 1:
            return ''
        begin = hits[0] + 1
    return key(acc, begin + offset)


def combine_calls(calls: list[str]) -> str:
    """Site call over rows (antibodies, or cell lines): decreased when a row is decreased and none increased or mixed;
    increased likewise; mixed when rows disagree or one is mixed; no_change when all are."""
    calls = set(calls)
    if not calls:
        return ''
    if calls <= {'no_change'}:
        return 'no_change'
    if 'mixed' in calls or {'decreased', 'increased'} <= calls:
        return 'mixed'
    return 'decreased' if 'decreased' in calls else 'increased'


def peptide_call(n_down: int, n_up: int) -> str:
    return 'mixed' if n_down and n_up else 'decreased' if n_down else 'increased' if n_up else 'no_change'


def maron_tables(inputs: Path, reference: dict, site_ids: set) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    s2 = pd.read_excel(inputs / INPUT['maron_s2'], sheet_name=0)
    s2_states, s2_peptides = defaultdict(set), defaultdict(set)
    counts = defaultdict(int)
    for _, r in s2.iterrows():
        mod = re.search(r'me([12])', str(r['Modification']))
        pos = pd.to_numeric(r['Position.In.Protein'], errors='coerce')
        site = place_maron_peptide(reference, r['Accession'], r['Sequence'], r['Begin.Pos'], int(pos)) if mod and pd.notna(pos) else ''
        counts['s2_rows'] += 1
        if site:
            counts['s2_rows_placed'] += 1
            counts['s2_rows_on_release_sites'] += site in site_ids
            s2_states[site].add(f'Rme{mod.group(1)}')
            s2_peptides[site].add((str(r['Sequence']), str(r['Modification'])))
    detected, measured = defaultdict(set), defaultdict(list)
    for sheet, drug in (('Control_vs_GSK591', 'GSK591'), ('Control_vs_MS023', 'MS023')):
        df = pd.read_excel(inputs / INPUT['maron_s3'], sheet_name=sheet)
        for _, r in df.iterrows():
            for pos, _degree in re.findall(r'R(\d+)me([12])', str(r['Modifications'])):
                counts[f's3_{drug}_site_rows'] += 1
                site = place_maron_peptide(reference, r['Accession'], r['Sequence'], r['Begin.Pos'], int(pos))
                if not site:
                    continue
                counts[f's3_{drug}_site_rows_placed'] += 1
                ip = clean(r['IP'])
                detected[(site, ip)].add((str(r['Sequence']), str(r['Modifications'])))
                diff, padj = pd.to_numeric(r[f'{drug}_vs_DMSO_diff'], errors='coerce'), pd.to_numeric(r[f'{drug}_vs_DMSO_p.adj'], errors='coerce')
                if pd.notna(diff) and pd.notna(padj):
                    measured[(site, ip, drug)].append((float(padj), float(diff)))
    detection = [{'site_id': s, 'maron_table': 'S2', 'antibody': '', 'methyl_states': join(ordered(states, list(rc.CODES['methyl_state']))),
                  'n_peptides': len(s2_peptides[s])} for s, states in s2_states.items() if s in site_ids]
    detection += [{'site_id': s, 'maron_table': 'S3', 'antibody': ip, 'methyl_states': '', 'n_peptides': len(peps)}
                  for (s, ip), peps in detected.items() if s in site_ids]
    inhibitors = []
    for (s, ip, drug), values in measured.items():
        if s not in site_ids:
            continue
        down = sum(p <= rc.PADJ_CUTOFF and d < 0 for p, d in values)
        up = sum(p <= rc.PADJ_CUTOFF and d > 0 for p, d in values)
        strongest = min(values)
        inhibitors.append({'site_id': s, 'antibody': ip, 'inhibitor': drug, 'n_peptides': len(values), 'n_decreased': down,
                           'n_increased': up, 'log2_diff_median': float(pd.Series([d for _, d in values]).median()),
                           'log2_diff_strongest': strongest[1], 'adj_p_min': strongest[0], 'call': peptide_call(down, up)})
    order = {c: k for k, c in enumerate(rc.ANTIBODY_ORDER)}
    detection = pd.DataFrame(detection, columns=[c for c, *_ in rc.MARON_DETECTION_COLUMNS])
    detection = detection.assign(_s=detection['site_id'].map(split_key), _a=detection['antibody'].map(lambda a: order.get(a, -1)))
    detection = detection.sort_values(['_s', 'maron_table', '_a']).drop(columns=['_s', '_a']).reset_index(drop=True)
    inhibitors = pd.DataFrame(inhibitors, columns=[c for c, *_ in rc.MARON_INHIBITOR_COLUMNS])
    inhibitors = inhibitors.assign(_s=inhibitors['site_id'].map(split_key), _a=inhibitors['antibody'].map(order))
    inhibitors = inhibitors.sort_values(['_s', 'inhibitor', '_a']).drop(columns=['_s', '_a']).reset_index(drop=True)
    return detection, inhibitors, dict(counts)


def carm1_table(inputs: Path, site_ids: set) -> pd.DataFrame:
    words = {'down': 'decreased', 'up': 'increased', 'unchanged': 'no_change', 'mixed': 'mixed'}
    layer = read_tsv(inputs / INPUT['carm1'])
    out = pd.DataFrame({
        'site_id': [key(a, p) for a, p in zip(layer['canonical_UniProtAC'], layer['position'])],
        'cell_line': layer['cell_line'],
        'n_peptides': layer['peptides'].astype(int),
        'log2_ko_over_wt': layer['median_log2_ko_wt'].astype(float),
        'p': layer['min_p'].astype(float),
        'n_decreased': layer['peptides_down'].astype(int),
        'n_increased': layer['peptides_up'].astype(int),
        'call': layer['call'].map(words),
    })
    if out['call'].isna().any():
        raise SystemExit(f"unknown CARM1 knockout calls: {sorted(set(layer['call']) - set(words))}")
    out = out[out['site_id'].isin(site_ids)]
    return out.assign(_s=out['site_id'].map(split_key)).sort_values(['_s', 'cell_line']).drop(columns='_s').reset_index(drop=True)


# ---- Structural context and motifs (the pipeline's own classifiers) ------------------------------------------------

def state_calls(sites: pd.DataFrame, records: pd.DataFrame, inputs: Path) -> pd.DataFrame:
    """The degree of methylation at each site, from every source that gives one: the records (ProMetheusDB), Maron
    Table S2, UniProt (experimental evidence counted, inferred annotations listed apart) and the studies that looked for
    one degree only (rc.STUDY_DESIGN_STATE)."""
    order = list(rc.CODES['methyl_state'])
    uniprot = read_tsv(inputs / INPUT['uniprot_states'])
    uniprot['site_id'] = [key(a, int(p)) for a, p in zip(uniprot['uniprot_accession'], uniprot['position'])]
    experimental = uniprot[uniprot['evidence_class'] == 'experimental'].groupby('site_id')['state'].agg(set)
    inferred = uniprot.groupby('site_id')['state'].agg(set)
    site_of = records['placed_site'].where(records['placed_site'] != '', records['pending_site'])
    design = records.assign(_site=site_of, _state=records['study'].map(rc.STUDY_DESIGN_STATE)).dropna(subset=['_state'])
    design = design[design['_site'] != ''].groupby('_site')['_state'].agg(set)
    out = pd.DataFrame(index=sites.index)
    out['uniprot_states'] = sites['site_id'].map(lambda s: join(ordered(experimental.get(s, set()), order)))
    out['uniprot_states_inferred'] = sites['site_id'].map(lambda s: join(ordered(inferred.get(s, set()) - experimental.get(s, set()), order)))
    out['states_by_study_design'] = sites['site_id'].map(lambda s: join(ordered(design.get(s, set()), order)))
    given = {'source_records': sites['methyl_states_reported'], 'maron2021_s2': sites['maron2021_s2_states'],
             'uniprot': out['uniprot_states'], 'study_design': out['states_by_study_design']}
    summary, evidence = [], []
    for i in sites.index:
        degrees, by = set(), []
        for name, column in given.items():
            found = {('Rme1' if x == 'Rme1' else 'Rme2') for x in column[i].split('; ') if x in order and x != 'Rme'}
            if found:
                degrees |= found
                by.append(name)
        summary.append('Rme1_Rme2' if len(degrees) == 2 else degrees.pop() if degrees else '')
        evidence.append(join(by))
    out['methyl_state'] = summary
    out['methyl_state_evidence'] = evidence
    return out


def context_annotations(sites: pd.DataFrame, inputs: Path) -> pd.DataFrame:
    frame = sites[['uniprot_accession', 'position']].rename(columns={'uniprot_accession': 'canonical_UniProtAC', 'position': 'corrected_position'})
    disorder = disorder_context.normalize_disorder_intervals(pd.read_csv(inputs / INPUT['disorder'], sep='\t', low_memory=False))
    by_acc = {a: g.reset_index(drop=True) for a, g in disorder.groupby('canonical_UniProtAC')}
    dis = disorder_context.annotate_frame(frame, by_acc, 'canonical_UniProtAC', 'corrected_position', rc.CONTEXT_WINDOW)

    intervals = pd.read_csv(inputs / INPUT['domains'], sep='\t', low_memory=False)
    intervals = intervals[intervals['interpro_type'].astype(str).isin({'domain', 'repeat'})].copy()  # scripts/05 default types
    for col in ('fragment_start', 'fragment_end'):
        intervals[col] = pd.to_numeric(intervals[col], errors='coerce')
    intervals = intervals.dropna(subset=['fragment_start', 'fragment_end'])
    intervals[['fragment_start', 'fragment_end']] = intervals[['fragment_start', 'fragment_end']].astype(int)
    domains_by_acc = {a: g.sort_values(['fragment_start', 'fragment_end', 'interpro_accession']).reset_index(drop=True)
                      for a, g in intervals.groupby('canonical_UniProtAC')}
    dom = pd.DataFrame([domain_context.classify_row(r, domains_by_acc, 'canonical_UniProtAC', 'corrected_position', rc.CONTEXT_WINDOW)
                        for _, r in frame.iterrows()], index=frame.index)
    ontology = json.loads(DOMAIN_ONTOLOGY.read_text())
    label = dom[['nearest_interpro_name', 'nearest_interpro_type', 'nearest_member_database_accessions']].fillna('').astype(str).agg(' '.join, axis=1)
    domain_class = [classify_by_patterns(t, ontology['domain_classes'], default='Other domain') if has else 'No domain annotation'
                    for t, has in zip(label, dom['has_domain_annotation'].fillna(False))]

    out = pd.DataFrame(index=sites.index)
    out['motifs'] = [join(m.removeprefix('motif_') for m, v in motif_odds.family_flags(w).items() if v) for w in sites['window_15']]
    out['disorder'] = dis['disorder_context_class']
    out['disorder_edge_distance'] = pd.to_numeric(dis['nearest_disorder_edge_distance'], errors='coerce').astype('Int64')
    out['domain'] = dom['domain_context_class']
    out['nearest_domain_interpro'] = dom['nearest_interpro_accession'].fillna('')
    out['nearest_domain_name'] = dom['nearest_interpro_name'].fillna('')
    out['nearest_domain_start'] = pd.to_numeric(dom['nearest_domain_start'], errors='coerce').astype('Int64')
    out['nearest_domain_end'] = pd.to_numeric(dom['nearest_domain_end'], errors='coerce').astype('Int64')
    out['nearest_domain_class'] = domain_class
    out['domain_side'] = dom['nearest_domain_side'].fillna('')
    out['domain_edge_distance'] = pd.to_numeric(dom['nearest_domain_edge_distance'], errors='coerce').astype('Int64')
    return out


# ---- Sites -----------------------------------------------------------------------------------------------------------

def build_sites(records, union, pending, reference, entry, previous_sites) -> pd.DataFrame:
    site_of = records['placed_site'].where(records['placed_site'] != '', records['pending_site'])
    observed = {key(a, p): o for a, p, o in zip(union['canonical_UniProtAC'], union['corrected_position'], union['arg_methyl_states_observed'])}
    proteome = remap.ProteomeIndex(reference)
    study_doi = {s[0]: s[3] for s in rc.STUDIES}
    rows = []
    for site, group in records.groupby(site_of, sort=False):
        if not site:
            continue
        acc, pos = split_key(site)
        seq = reference[acc]
        if seq[pos - 1] != 'R':
            raise SystemExit(f'{site} is not an arginine of the reference sequence')
        studies = ordered(group['study'], rc.STUDY_ORDER)
        states = ordered(clean(observed.get(site, '')).split(';') if site in observed else group['reported_methyl_state'], list(rc.CODES['methyl_state']) + [''])
        window_31 = ''.join(seq[k] if 0 <= k < len(seq) else '_' for k in range(pos - 16, pos + 15))
        shared = sorted({a for a, _ in proteome.find(seq[max(0, pos - 16):pos + 15]) if a != acc}) if len(seq[max(0, pos - 16):pos + 15]) >= 15 else []
        flags = []
        if site in pending:
            flags.append('one_sided_window')
        if (group['accepted_exception'] != '').any():
            flags.append('paralog_window')
        if site not in pending and group['placement_code'].str.contains('near_match').all():
            flags.append('near_match_only')
        if site not in pending and group['evidence_check'].eq('agrees_reported_arginine').all():
            flags.append('arginine_from_position')
        if shared:
            flags.append('shared_sequence')
        if len(studies) == 1:
            flags.append('single_study')
        reported = sorted({int(float(p)) for p in group['reported_position'] if clean(p)})
        in_previous = group['_previous_lifted'].eq(site).any() or site in previous_sites
        renumbered = group['reason_code'].eq('reference_renumbered').any()
        if in_previous:
            reasons = ['reference_renumbered'] if renumbered else []
        else:
            reasons = (['one_sided_window'] if site in pending else []) + sorted(set(group.loc[group['_previous_lifted'] != site, 'reason_code']) - {''})
        rows.append({
            'site_id': site,
            'uniprot_accession': acc,
            'uniprot_entry_name': entry(acc, 'Entry Name'),
            'gene': entry(acc, 'Gene Names (primary)'),
            'protein_name': re.split(r' [(\[]', entry(acc, 'Protein names'))[0],
            'position': pos,
            'residue': seq[pos - 1],
            'window_15': window_31[8:23],
            'window_31': window_31,
            'protein_length': len(seq),
            'reference_sequence_md5': hashlib.md5(seq.encode()).hexdigest().upper(),
            'methyl_states_reported': join(s for s in states if s),
            'n_studies': len(studies),
            'studies': join(studies),
            'study_dois': join(study_doi[s] for s in studies),
            'source_tables': join(ordered(group['source_table'], list(rc.SOURCE_TABLES))),
            'n_records': len(group),
            'record_ids': join(sorted(group['record_id'])),
            'reported_positions': join(reported) if any(p != pos for p in reported) else '',
            'flags': join(flags),
            'same_window_accessions': join(shared),
            'same_window_genes': join(entry(a, 'Gene Names (primary)').split(';')[0].strip() or a for a in shared),
            'in_previous_version': bool(in_previous),
            'change_since_previous_version': 'kept' if in_previous else 'added',
            'change_reasons': join(reasons),
        })
    sites = pd.DataFrame(rows)
    return sites.assign(_s=sites['site_id'].map(split_key)).sort_values('_s').drop(columns='_s').reset_index(drop=True)


# ---- Changes since the previous version ----------------------------------------------------------------------------

def build_changes(records, sites, previous_sites, shift, entry) -> pd.DataFrame:
    site_of = records['placed_site'].where(records['placed_site'] != '', records['pending_site'])
    now = set(sites['site_id'])
    before = {lifted(s, shift): s for s in previous_sites}
    pending = set(records['pending_site']) - {''}
    at_now = records.groupby(site_of).groups
    at_before = records.groupby('_previous_lifted').groups
    reported_gene = records.groupby('previous_site')['reported_gene'].agg(lambda g: next((x for x in g if x), ''))
    out = []
    for site in sorted(now | set(before), key=split_key):
        acc, pos = split_key(site)
        in_before, in_now = site in before, site in now
        if in_before and in_now:
            # The records that moved onto or off the site; with none of the previous ones left, its records were replaced.
            change = 'kept'
            gained = records.loc[at_now.get(site, [])]
            gained = gained[gained['_previous_lifted'] != site]
            lost = records.loc[at_before.get(site, [])]
            lost = lost[site_of[lost.index] != site]
            involved = pd.concat([gained, lost])
            kept_records = len(at_before.get(site, [])) - len(lost)
            codes = (['reference_renumbered'] if before[site] != site else []) + sorted(set(involved['reason_code']) - {'', 'reference_renumbered'})
        elif in_before:
            change = 'removed'
            rows = records.loc[at_before.get(site, [])]
            involved = rows[site_of[rows.index] != site]
            codes = sorted(set(involved['reason_code']) - {''})
        else:
            change = 'added'
            rows = records.loc[at_now.get(site, [])]
            involved = rows[rows['_previous_lifted'] != site]
            codes = (['one_sided_window'] if site in pending else []) + sorted(set(involved['reason_code']) - {''})
        if change == 'kept':
            elsewhere = sorted((set(site_of[lost.index]) | set(gained['previous_site'])) - {'', site, before[site]}, key=split_key)
            records_change = 'same' if involved.empty else 'replaced' if kept_records == 0 else 'changed'
        else:
            elsewhere = sorted(set(site_of[involved.index] if change == 'removed' else involved['previous_site']) - {'', site}, key=split_key)
            records_change = ''
        out.append({
            'site_id': site,
            'previous_site_id': before.get(site, ''),
            'gene': entry(acc, 'Gene Names (primary)') or reported_gene.get(before.get(site, ''), ''),
            'uniprot_accession': acc,
            'position': pos,
            'in_previous_version': in_before,
            'in_this_version': in_now,
            'change': change,
            'reason_codes': join(codes),
            'records_change': records_change,
            'records_now_at_or_previously_at': join(elsewhere),
            'record_ids': join(sorted(involved['record_id'])),
        })
    return pd.DataFrame(out)


# ---- Studies ---------------------------------------------------------------------------------------------------------

def build_studies(records, sites) -> pd.DataFrame:
    placed = sites[~sites['flags'].str.contains('one_sided_window')]
    per_site = placed['studies'].str.split('; ')
    out = []
    for sid, label, citation, doi in rc.STUDIES:
        group = records[records['study'] == sid]
        supported = per_site[per_site.map(lambda s: sid in s)]
        out.append({'study_id': sid, 'label': label, 'source_table': join(ordered(group['source_table'], list(rc.SOURCE_TABLES))),
                    'citation': citation, 'doi': doi, 'licence': rc.STUDY_LICENCE.get(sid, rc.MARON_S5),
                    'n_records': len(group), 'n_records_placed': int((group['placed_site'] != '').sum()),
                    'n_sites': len(supported), 'n_sites_only_here': int((supported.map(len) == 1).sum()),
                    'note': rc.STUDY_NOTE.get(sid, '')})
    unknown = set(records['study']) - set(rc.STUDY_ORDER)
    if unknown:
        raise SystemExit(f'studies missing from scripts/release_codes.py: {sorted(unknown)}')
    return pd.DataFrame(out)


def validation(root: Path, check: pd.DataFrame) -> dict:
    """The outcome of the build's two gates (hpc/build_release.sh stops on either)."""
    audit = root / 'qc' / 'remap_evidence_check' / 'union_site_residue_audit.tsv'
    failures = root / 'qc' / 'remap_evidence_check' / 'remap_evidence_check_failures.tsv'
    return {'evidence_check': {rc.VERDICT_CODE[k]: int(n) for k, n in check['verdict'].value_counts().items()},
            'rows_against_evidence_not_accepted': len(read_tsv(failures)) if failures.exists() else 0,
            'union_sites_not_arginine': len(read_tsv(audit)) if audit.exists() else 0}


def peptide_start_agreement(rows: pd.DataFrame, typed: pd.DataFrame, reference: dict, path: Path) -> dict:
    """How often counting from a known peptide start agrees with placements the pipeline made from plain peptides by
    other means (a single arginine, or the reported arginine inside the peptide)."""
    coordinates = remap.load_peptide_coordinates(path if path.exists() else None, reference)
    agree = disagree = 0
    for i in rows.index[(rows['is_position_resolved'] == 'True') & (rows['remap_status'] != 'sequence_peptide_start')]:
        peptides = [a for a in remap.build_anchors(typed.loc[i]) if a['anchor_type'] == 'peptide']
        acc = rows.at[i, 'corrected_canonical_UniProtAC']
        position, _ = remap.peptide_start_position(peptides, acc, int(float(rows.at[i, 'position'])), reference.get(acc, ''), coordinates) if peptides else (None, '')
        if position is not None:
            agree += position == int(float(rows.at[i, 'corrected_position']))
            disagree += position != int(float(rows.at[i, 'corrected_position']))
    return {'start_rule_agree': agree, 'start_rule_disagree': disagree}


# ---- Writing ---------------------------------------------------------------------------------------------------------

def interaction_tables(folder: Path, sites: pd.DataFrame) -> dict:
    """The tables scripts/37 wrote for these sites (--tables-out), with each site's interaction status: experimental
    when a crosslink, RBDmap fragment or PDB interface counts as evidence there, predicted_only when only a predicted
    interface or confident predicted contact does, none otherwise (IntAct records never set it, as on the website)."""
    t = {name: read_tsv(folder / f'{name}.tsv') for name in ('site_interactions', 'residue_evidence', 'protein_partners', 'complexes',
                                                            'interaction_statistics', 'protein_interaction_summary')}
    ev = t['site_interactions']
    analysed = sites[~sites['flags'].str.contains('one_sided_window', regex=False)]
    unknown = set(ev['site_id']) - set(analysed['site_id'])
    if unknown:
        raise SystemExit(f'{folder}: interaction rows for sites that are not analysed in the release: {sorted(unknown)[:5]}')
    if set(t['protein_interaction_summary']['uniprot_accession']) != set(analysed['uniprot_accession']):
        raise SystemExit(f'{folder}: the interaction tables are for another set of proteins')
    for frame, column in ((ev, 'counts_as_evidence'), (t['protein_partners'], 'partner_methylated')):
        frame[column] = frame[column] == 'true'
    counted = ev[ev['counts_as_evidence']]
    experimental = {'rna_crosslink', 'rbdmap_fragment', 'dna_crosslink', 'pdb_protein', 'pdb_rna', 'pdb_dna', 'pdb_self'}
    predicted = {'predicted_interface', 'predicted_contact'}
    kinds = counted.groupby('site_id')['evidence_type'].agg(set)
    status = {s: 'experimental' if k & experimental else 'predicted_only' if k & predicted else 'none' for s, k in kinds.items()}
    evidence = {s: join(ordered(k, rc.EVIDENCE_ORDER)) for s, k in kinds.items()}
    return {'sites': ev, 'residues': t['residue_evidence'], 'partners': t['protein_partners'], 'complexes': t['complexes'],
            'statistics': t['interaction_statistics'], 'proteins': t['protein_interaction_summary'], 'status': status, 'evidence': evidence}


def write_table(df: pd.DataFrame, columns: list[tuple], path: Path) -> None:
    """Write a release TSV in the declared column order and types: integers without decimals, numbers in up to 6
    significant digits, booleans true/false, empty for missing. Fails on a tab or newline in a value, or on ';' in a
    value of a column that is not a list."""
    out = pd.DataFrame(index=df.index)
    for name, kind, *_ in columns:
        col = df[name]
        if kind == 'integer':
            out[name] = [('' if pd.isna(v) or v == '' else str(int(v))) for v in col]
        elif kind == 'number':
            out[name] = [('' if pd.isna(v) or v == '' else f'{float(v):.6g}') for v in col]
        elif kind == 'boolean':
            out[name] = ['true' if bool(v) else 'false' for v in col]
        else:
            out[name] = ['' if v is None or (isinstance(v, float) and pd.isna(v)) or v is pd.NA else str(v) for v in col]
        bad = out[name].str.contains('[\t\n\r]', regex=True)
        if not kind.startswith('list'):
            bad |= out[name].str.contains(';', regex=False)
        if bad.any():
            raise SystemExit(f'{path.name}: column {name} has a separator inside a value: {out.loc[bad, name].iloc[0]!r}')
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, sep='\t', index=False, lineterminator='\n')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results-root', required=True, help='results root of the pipeline run (remapped/, integrated/, qc/)')
    ap.add_argument('--previous', required=True, help='folder of the previous release, e.g. releases/v1.0.0')
    ap.add_argument('--inputs', default='inputs', help='frozen inputs with MANIFEST.tsv')
    ap.add_argument('--allow', default='config/remap_evidence_allowlist.tsv')
    ap.add_argument('--version', required=True, help='version of this release, e.g. v1.1.0')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--date', default=dt.date.today().isoformat(), help='release date (default today)')
    ap.add_argument('--interactions', default='', help='the interaction tables scripts/37 wrote for these sites (--tables-out)')
    args = ap.parse_args()
    root, inputs, outdir, prev_dir = Path(args.results_root), Path(args.inputs), Path(args.outdir), Path(args.previous)
    version, previous = args.version, prev_dir.name

    reference = remap.parse_fasta(inputs / INPUT['fasta'])
    meta = read_tsv(inputs / INPUT['entries']).set_index('Entry')

    def entry(acc: str, field: str) -> str:
        return clean(meta[field].get(acc, ''))

    rows = read_tsv(root / 'remapped' / 'human_arg_methyl_source_rows_with_remap.tsv')
    typed = pd.read_csv(root / 'remapped' / 'human_arg_methyl_source_rows_with_remap.tsv', sep='\t', low_memory=False)
    typed = remap.annotate_window_lengths(typed) if 'seq_window_window_length' not in typed else typed
    check = read_tsv(root / 'qc' / 'remap_evidence_check' / 'remap_evidence_check_rows.tsv')
    union = read_tsv(root / 'integrated' / 'human_arg_methyl_union_dedup_by_site.tsv')
    pending_table = read_tsv(root / 'integrated' / 'pending_one_sided_window_sites.tsv')
    pending = {key(a, p) for a, p in zip(pending_table['canonical_UniProtAC'], pending_table['position'])}
    pending_of = {int(i): key(a, p) for a, p, ids in zip(pending_table['canonical_UniProtAC'], pending_table['position'], pending_table['record_ids'])
                  for i in ids.split('; ')}
    allow = pd.read_csv(args.allow, sep='\t', dtype=str, comment='#').fillna('')
    allowed = {(r.source_table, r.reported_accession, r.reported_position, r.placed_site): r.reason for r in allow.itertuples()}

    prev_sites = read_tsv(prev_dir / f'human_methylarginine_sites_{previous}.tsv')
    prev_records = read_tsv(prev_dir / 'supporting' / f'source_records_{previous}.tsv')
    shift, broken = reference_lift(inputs)
    stuck = sorted(s for s in prev_sites['site_id'] if split_key(s)[0] in broken)
    if stuck:
        raise SystemExit(f'previous sites on entries whose numbering cannot be carried over: {stuck}')

    # This version reads the previous version's source records in the same order, as reported (before the corrections
    # of source_correction), and after them the records read for the first time.
    n_prev = len(prev_records)
    head = rows.iloc[:n_prev]
    reported_peptide = head['peptide_sequence_reported'].map(clean).where(head['peptide_sequence_reported'].map(clean) != '', head['peptide_sequence'].map(clean)) \
        if 'peptide_sequence_reported' in head else head['peptide_sequence'].map(clean)
    same = (len(rows) >= n_prev and len(check) == len(rows)
            and (head['source_table'].map(rc.SOURCE_TABLE_CODE).values == prev_records['source_table'].values).all()
            and (head['substrate_UniProtAC_original'].map(clean).values == prev_records['reported_accession'].values).all()
            and (head['position'].map(clean).values == prev_records['reported_position'].values).all()
            and (reported_peptide.values == prev_records['peptide_sequence'].values).all()
            and (head['seq_window'].map(clean).values == prev_records['seq_window'].values).all())
    if not same:
        raise SystemExit('this run does not hold the previous release\'s source records in the same order')

    records = build_records(rows, typed, check, allowed, prev_records, pending_of, shift, reference)
    union_keys = {key(a, p) for a, p in zip(union['canonical_UniProtAC'], union['corrected_position'])}
    if union_keys != set(records['placed_site']) - {''}:
        raise SystemExit('the placed records and the union disagree')
    sites = build_sites(records, union, pending, reference, entry, {lifted(x, shift) for x in prev_sites['site_id']})
    site_ids = set(sites['site_id'])

    detection, inhibitors, maron_counts = maron_tables(inputs, reference, site_ids)
    s2 = detection[detection['maron_table'] == 'S2'].set_index('site_id')['methyl_states']
    s3 = detection[detection['maron_table'] == 'S3'].groupby('site_id')['antibody'].agg(lambda a: join(ordered(a, rc.ANTIBODY_ORDER)))
    sites['maron2021_s2_states'] = sites['site_id'].map(s2).fillna('')
    sites['maron2021_s3_antibodies'] = sites['site_id'].map(s3).fillna('')
    sites = pd.concat([sites, state_calls(sites, records, inputs)], axis=1)
    for drug in rc.CODES['inhibitor']:
        calls = inhibitors[inhibitors['inhibitor'] == drug].groupby('site_id')['call'].agg(list)
        sites[f'maron2021_{drug}_call'] = sites['site_id'].map(lambda s: combine_calls(calls.get(s, [])))
    carm1 = carm1_table(inputs, site_ids)
    carm1_calls = carm1.groupby('site_id')['call'].agg(list)
    sites['carm1_ko_call'] = sites['site_id'].map(lambda s: combine_calls(carm1_calls.get(s, [])))
    sites = pd.concat([sites, context_annotations(sites, inputs)], axis=1)
    ix = interaction_tables(Path(args.interactions), sites) if args.interactions else None
    analysed = ~sites['flags'].str.contains('one_sided_window', regex=False)   # sites kept pending are left out
    sites['interaction_status'] = sites['site_id'].map(ix['status']).fillna('none').where(analysed, '') if ix else ''
    sites['interaction_evidence'] = sites['site_id'].map(ix['evidence']).fillna('') if ix else ''
    sites = sites[[c for c, *_ in rc.SITE_COLUMNS]]

    changes = build_changes(records, sites, list(prev_sites['site_id']), shift, entry)
    studies = build_studies(records, sites)
    proteins, extra = release_tables.proteins_table(root, inputs, sites, records, entry, reference)
    statistics = release_tables.statistics(root, proteins)
    if ix:
        status = sites.groupby('uniprot_accession')['interaction_status'].agg(Counter)
        proteins['n_sites_experimental'] = proteins['uniprot_accession'].map(lambda a: status[a]['experimental'])
        proteins['n_sites_predicted_only'] = proteins['uniprot_accession'].map(lambda a: status[a]['predicted_only'])
        summary = ix['proteins'].set_index('uniprot_accession')
        for column in ('go_rna_binding', 'go_dna_binding', 'n_string_partners', 'n_partners_methylated', 'n_complexes'):
            proteins[column] = proteins['uniprot_accession'].map(summary[column]).fillna('')   # empty: only pending sites
        statistics = pd.concat([statistics, ix['statistics']], ignore_index=True)
    else:
        for column, *_ in rc.PROTEIN_COLUMNS[-7:]:
            proteins[column] = ''
    disease_terms = release_tables.disease_terms_table(root, proteins)

    files = {
        'sites': (f'human_methylarginine_sites_{version}.tsv', sites, rc.SITE_COLUMNS),
        'studies': (f'studies_{version}.tsv', studies, rc.STUDY_COLUMNS),
        'records': (f'supporting/source_records_{version}.tsv', records, rc.RECORD_COLUMNS),
        'changes': (f'supporting/site_changes_{previous}_to_{version}.tsv', changes, rc.CHANGE_COLUMNS),
        'maron_detection': (f'supporting/maron2021_detection_{version}.tsv', detection, rc.MARON_DETECTION_COLUMNS),
        'maron_inhibitors': (f'supporting/maron2021_inhibitors_{version}.tsv', inhibitors, rc.MARON_INHIBITOR_COLUMNS),
        'carm1_ko': (f'supporting/carm1_ko_{version}.tsv', carm1, rc.CARM1_COLUMNS),
        'proteins': (f'human_methylarginine_proteins_{version}.tsv', proteins, rc.PROTEIN_COLUMNS),
        'statistics': (f'statistics_{version}.tsv', statistics, rc.STATISTICS_COLUMNS),
        'disease_terms': (f'supporting/disease_terms_{version}.tsv', disease_terms, rc.DISEASE_TERM_COLUMNS),
        'density_windows': (f'supporting/density_windows_{version}.tsv', extra['density_windows'], rc.DENSITY_COLUMNS),
        'hpa_locations': (f'supporting/hpa_locations_{version}.tsv', extra['hpa'], rc.HPA_COLUMNS),
        'other_ptm_sites': (f'supporting/other_ptm_sites_{version}.tsv', extra['other_ptm_sites'], rc.OTHER_PTM_COLUMNS),
        'condensates': (f'supporting/condensates_{version}.tsv', release_tables.condensates_table(root, set(sites['uniprot_accession'])), rc.CONDENSATE_COLUMNS),
        'disease_families': (f'supporting/disease_families_{version}.tsv', release_tables.disease_families_table(root), rc.DISEASE_FAMILY_COLUMNS),
    }
    if ix:
        files.update({
            'site_interactions': (f'supporting/site_interactions_{version}.tsv', ix['sites'], rc.SITE_INTERACTION_COLUMNS),
            'residue_evidence': (f'supporting/residue_evidence_{version}.tsv', ix['residues'], rc.RESIDUE_EVIDENCE_COLUMNS),
            'protein_partners': (f'supporting/protein_partners_{version}.tsv', ix['partners'], rc.PARTNER_COLUMNS),
            'complexes': (f'supporting/complexes_{version}.tsv', ix['complexes'], rc.COMPLEX_COLUMNS),
        })
    for folder in (outdir / 'supporting', outdir):
        if folder.is_dir():
            for p in folder.iterdir():
                if p.is_file():
                    p.unlink()
    for name, frame, columns in files.values():
        write_table(frame, columns, outdir / name)
    import release_workbook
    release_workbook.write(outdir / f'human_methylarginine_sites_{version}.xlsx', read_tsv(outdir / files['sites'][0]),
                           rc.SITE_COLUMNS, title=f'Sites {version}')

    counts = release_docs.counts(sites, records, changes, studies, maron_counts, prev_sites)
    counts.update(peptide_start_agreement(rows, typed, reference, root / 'parsed_maron' / 'maron_peptide_coordinates.tsv'))
    context = {'version': version, 'previous': previous, 'files': files, 'counts': counts, 'inputs': inputs, 'root': root,
               'uniprot_release': json.loads((inputs / INPUT['uniprot_release']).read_text()), 'records': records,
               'sites': sites, 'changes': changes, 'studies': studies, 'release_date': args.date,
               'validation': validation(root, check)}
    release_docs.write_all(outdir, context)
    print(json.dumps({k: v for k, v in counts.items() if not isinstance(v, dict)}))


if __name__ == '__main__':
    main()

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
  supporting/disease_terms_V.tsv, disease_terms_tested_V.tsv, disease_ontology_V.tsv, density_windows_V.tsv, hpa_locations_V.tsv, other_ptm_sites_V.tsv
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
from common import classify_by_patterns, structure_context  # noqa: E402


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


def recommended_name(names: str) -> str:
    """UniProt's recommended name: its 'Protein names' text up to the first alternative name, a top-level (...) or [...]
    group that follows a space and is followed by a space or the end. Brackets that belong to the name stay: 'DNA
    (cytosine-5)-methyltransferase 1', 'Isocitrate dehydrogenase [NADP], mitochondrial', '[Pyruvate dehydrogenase
    (acetyl-transferring)] kinase isozyme 1, mitochondrial'."""
    depth, start = 0, -1
    for k, ch in enumerate(names):
        if ch in '([':
            if depth == 0:
                start = k
            depth += 1
        elif ch in ')]' and depth:
            depth -= 1
            if depth == 0 and start > 0 and names[start - 1] == ' ' and (k + 1 == len(names) or names[k + 1] == ' '):
                return names[:start - 1]
    return names


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
    if status in rc.PLACEMENT_NOT_LOCALISED:
        return 'not_localised'
    if status in rc.PLACEMENT_SEVERAL | rc.PLACEMENT_OTHER_PROTEIN | {'unresolved_study_not_found'} and row.get('site_localisation', ''):
        return 'not_unique_in_proteome'
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
        study = rc.study_of(r['source_table'], r['source'], r['reference'])
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
            'reported_sample': clean(r.get('study_sample', '')),
            'reported_enrichment': clean(r.get('study_enrichment', '')),
            'reported_detail': record_detail(r),
            'evidence_kind': kind,
            'evidence_sequence': text,
            'evidence_window': window_of(typed.loc[i]) if kind == 'marked_peptide' else '',
            'source_correction': clean(r.get('source_correction', '')),
            'evidence_as_reported': clean(r.get('peptide_sequence_reported', '')),
            'localised_by': clean(r.get('site_localisation', '')),
            'not_localised': '' if clean(r.get('site_localisation', '')) else clean(r.get('study_not_localised', '')),
            'candidate_sites': as_list(r.get('candidate_sites', '')),
            'reported_site': clean(r.get('reported_site', '')),
            **{t: clean(r.get(t, '')) for t in TEST_FIELDS},
            '_evidence_sites': clean(r.get('evidence_sites', '')),
            'placed_site': new,
            'placement_code': r['remap_status'],
            'evidence_check': rc.VERDICT_CODE[c['verdict']],
            'accepted_exception': allowed.get((r['source_table'], c['reported_accession'], c['reported_position'], c['placed_site']), ''),
            'pending_site': pending_of.get(i + 1, ''),
            'previous_site': old_published,
            '_previous_lifted': old,
            # where the previous version had the record, placed or (from v1.1.0) kept pending
            '_previous_at': old or lifted(previous.at[i, 'pending_site'] if 'pending_site' in previous and i < len(previous) else '', shift),
            'change_since_previous_version': change,
            'reason_code': reason,
            '_rme2_type': clean(r.get('study_rme2_type', '')),
        })
    return pd.DataFrame(out)


TEST_FIELDS = ['test_evidence_at_site', 'test_named_entry', 'test_unique_in_proteome', 'test_methyl_on_arginine', 'test_localised']
AMBIGUOUS = {**{c: 'several_sites' for c in rc.PLACEMENT_SEVERAL}, **{c: 'other_protein' for c in rc.PLACEMENT_OTHER_PROTEIN},
             'unresolved_study_not_found': 'not_in_reference'}


def lower_confidence_reason(localised: str, why: str, code: str) -> str:
    """Why a record that is not placed adds no site (CODES['lower_confidence']), empty when it says nothing."""
    if localised:
        return AMBIGUOUS.get(code, '')
    return why or ('several_residues' if code == 'unresolved_not_localised' else '')


def lower_confidence(records: pd.DataFrame, sites: set) -> pd.DataFrame:
    """The records that add no site but say where a methyl group could be (supporting/lower_confidence_evidence): not
    localised (a study's own record, or a compiled record no study table stands in for), or localised but not on one
    site of a protein the record names. Rows of the compiled tables a study's own table supersedes are left out."""
    keep = (records['placed_site'] == '') & (records['pending_site'] == '') & (records['placement_code'] != 'superseded_by_study_table')
    rows = records[keep].copy()
    rows['reason'] = [lower_confidence_reason(*v) for v in zip(rows['localised_by'], rows['not_localised'], rows['placement_code'])]
    rows = rows[rows['reason'] != '']
    rows['candidates_in_release'] = rows['candidate_sites'].map(lambda c: '; '.join(x for x in c.split('; ') if x in sites) if c else '')
    return rows[[c for c, *_ in rc.LOWER_CONFIDENCE_COLUMNS]]


def earlier_sites(prev_dir: Path, previous: str, prev_sites: pd.DataFrame, records: pd.DataFrame, shift: dict) -> tuple[dict, list]:
    """The sites of the earlier releases, in this release's numbering. Returns ({version: {site: record_ids}}, [(version,
    site as published, [sites now])]). The previous release's sites carry the records it placed there. The release
    before it is read from the previous release's site_changes table: its kept sites are sites of the previous release;
    a site it removed carries the records the table gives, at the sites they were moved to, or, when they were moved
    nowhere, at the site as published."""
    members, published = defaultdict(dict), []
    at = records[records['_previous_at'] != ''].groupby('_previous_at')['record_id'].agg(list)
    for site in prev_sites['site_id']:
        now = lifted(site, shift)
        members[previous][now] = at.get(now, [])
        published.append((previous, site, [now]))
    for path in sorted((prev_dir / 'supporting').glob(f'site_changes_*_to_{previous}.tsv')):
        version = path.name[len('site_changes_'):-len(f'_to_{previous}.tsv')]
        for r in read_tsv(path).to_dict('records'):
            if r['in_previous_version'] != 'true':
                continue
            if r['change'] == 'kept':
                published.append((version, r['previous_site_id'], [lifted(r['site_id'], shift)]))
                continue
            ids = [int(x) for x in r['record_ids'].split('; ') if x]
            now = [lifted(x, shift) for x in r['records_now_at_or_previously_at'].split('; ') if x] or [r['site_id']]
            for site in now:
                members[version][site] = sorted(set(members[version].get(site, [])) | set(ids))
            published.append((version, r['site_id'], now))
    return dict(members), published


def is_arginine(site: str, reference: dict) -> bool:
    acc, pos = split_key(site)
    return acc in reference and 0 < pos <= len(reference[acc]) and reference[acc][pos - 1] == 'R'


def lower_confidence_sites(records: pd.DataFrame, site_ids: set, reference: dict, members: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """The lower-confidence sites (CODES['tier']): every arginine that is not a site of this release and that a record
    not placed reports (reported_site, the compiled records a study's own table supersedes left out), or that an earlier
    release had (members). Returns the records exploded to these sites (one row per record and site, placed_site set
    to the site, for build_sites) and, per site, the tests of the record that passes most of them at it: a record's
    evidence_at_site is read again at each site (is the site among the arginines its evidence names), its other tests
    are its own. A site is counted in the analyses with lower-confidence sites when a record covers it with its own
    evidence and leaves the methyl group on an arginine (evidence_at_site pass or curated, methyl_on_arginine not
    fail)."""
    unplaced = (records['placed_site'] == '') & (records['pending_site'] == '') & (records['placement_code'] != 'superseded_by_study_table')
    at, versions = defaultdict(set), defaultdict(set)
    for i, site in zip(records.index[unplaced], records.loc[unplaced, 'reported_site']):
        if site and site not in site_ids:
            at[site].add(i)
    index_of = dict(zip(records['record_id'].astype(int), records.index))
    for version, sites in members.items():
        for site, ids in sites.items():
            if site not in site_ids and is_arginine(site, reference):
                at[site] |= {index_of[int(x)] for x in ids if int(x) in index_of}
                versions[site].add(version)
    rows, exploded = [], []
    for site, idx in at.items():
        tested = []
        for i in sorted(idx):
            r = records.loc[i]
            if not r['test_localised']:
                tests = dict.fromkeys(TEST_FIELDS, 'untested')
            else:
                evidence = set(r['_evidence_sites'].split('; ')) - {''}
                here = 'curated' if r['test_evidence_at_site'] == 'curated' and site == r['reported_site'] else 'pass' if site in evidence else 'fail'
                # named_entry at another site than the record reports: the same entry, or one of its reported accessions
                acc = split_key(site)[0]
                named = r['test_named_entry'] if site == r['reported_site'] else 'pass' if (
                    (r['test_named_entry'] == 'pass' and r['reported_site'].startswith(acc + ':R'))
                    or acc in {a.split('-')[0] for a in r['reported_accession'].split('; ')}) else 'fail'
                tests = {'test_evidence_at_site': here, 'test_named_entry': named, **{t: r[t] for t in TEST_FIELDS[2:]}}
            tested.append((i, tests))
        def rank(item):
            i, t = item
            return (t['test_evidence_at_site'] in ('pass', 'curated'), t['test_methyl_on_arginine'] != 'fail', t['test_named_entry'] == 'pass',
                    t['test_unique_in_proteome'] in ('pass', 'paralogs'), t['test_localised'] == 'pass',
                    records.at[i, 'placement_code'] != 'superseded_by_study_table', -int(records.at[i, 'record_id']))
        best = max(tested, key=rank)[1] if tested else dict.fromkeys(TEST_FIELDS, 'untested')
        reasons = {lower_confidence_reason(records.at[i, 'localised_by'], records.at[i, 'not_localised'], records.at[i, 'placement_code']) for i in idx}
        rows.append({'site_id': site,
                     'in_analyses_with_lower_confidence': any(t['test_evidence_at_site'] in ('pass', 'curated') and t['test_methyl_on_arginine'] != 'fail' for _, t in tested),
                     **best,
                     'tests_failed': join(k.removeprefix('test_') for k in TEST_FIELDS if best[k] == 'fail'),
                     'reasons': join(ordered(reasons - {''}, list(rc.CODES['lower_confidence']))),
                     'in_earlier_versions': join(sorted(versions.get(site, set())))})
        if idx:
            exploded.append(records.loc[sorted(idx)].assign(placed_site=site, pending_site=''))
    tests = pd.DataFrame(rows, columns=['site_id', 'in_analyses_with_lower_confidence', *TEST_FIELDS, 'tests_failed', 'reasons', 'in_earlier_versions'])
    exploded = pd.concat(exploded, ignore_index=True) if exploded else records.iloc[:0]
    return exploded, tests


def published_sites(published: list, site_ids: set, tests: pd.DataFrame) -> pd.DataFrame:
    """Every site of the earlier releases and what it is in this release (supporting/published_sites)."""
    failed = dict(zip(tests['site_id'], tests['tests_failed']))
    rows = []
    for version, site, now in published:
        status = 'site' if any(x in site_ids for x in now) else 'lower_confidence' if any(x in failed for x in now) else 'not_in_reference'
        rows.append({'version': version, 'published_site_id': site, 'site_ids': join(x for x in now if x in site_ids or x in failed),
                     'status': status, 'tests_failed': join(sorted({f for x in now if x not in site_ids for f in failed.get(x, '').split('; ') if f},
                                                                    key=list(rc.CODES['site_test']).index))})
    return pd.DataFrame(rows)


def record_detail(r: pd.Series) -> str:
    """What a staged study says about a record besides its sample and enrichment (scripts/00_stage_study_sites.py)."""
    parts = []
    if clean(r.get('study_rme2_type', '')):
        parts.append('neutral loss ' + clean(r['study_rme2_type']).replace('; ', ' and '))
    if clean(r.get('study_localisation', '')):
        parts.append('localisation probability ' + clean(r['study_localisation']))
    if clean(r.get('study_rows', '')):
        parts.append(clean(r['study_rows']) + (' source row' if clean(r['study_rows']) == '1' else ' source rows'))
    if clean(r.get('study_tables', '')):
        parts.append('tables ' + clean(r['study_tables']).replace('; ', ', '))
    if clean(r.get('study_note', '')):
        parts.append(clean(r['study_note']).replace('; ', ', '))
    return ', '.join(parts)


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
    out['structure_context'] = [structure_context(a, b, c) for a, b, c in zip(out['disorder'], out['domain'], out['domain_edge_distance'])]
    return out


# ---- Sites -----------------------------------------------------------------------------------------------------------

def build_sites(records, union, pending, reference, entry, previous_sites, placed_here: bool = True) -> pd.DataFrame:
    """placed_here False: the records are not placed on these sites (lower-confidence sites), so no placement flag."""
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
        # the paralog rule: placed only because the study named this one of several entries holding the sequence
        paralog = group['placement_code'].eq('sequence_paralog_named')
        paralogs = sorted({x.split(':R')[0] for e in group.loc[paralog, '_evidence_sites'] for x in e.split('; ') if x} - {acc}) \
            if placed_here and paralog.any() and not (~paralog).any() else []
        if paralogs:
            flags.append('paralog_assigned')
        if shared:
            flags.append('shared_sequence')
        if len(studies) == 1:
            flags.append('single_study')
        if studies == ['Ma2023']:
            flags.append('xenograft_only')
        reported = sorted({int(float(p)) for p in group['reported_position'] if clean(p)})
        in_previous = group['_previous_lifted'].eq(site).any() or site in previous_sites
        renumbered = group['reason_code'].eq('reference_renumbered').any()
        if in_previous:
            reasons = ['reference_renumbered'] if renumbered else []
        else:
            reasons = (['one_sided_window'] if site in pending else []) + sorted(set(group.loc[group['_previous_lifted'] != site, 'reason_code']) - {''})
            if group['change_since_previous_version'].eq('new_record').any():
                reasons = sorted(set(reasons) | {'new_record'})
        rows.append({
            'site_id': site,
            'uniprot_accession': acc,
            'uniprot_entry_name': entry(acc, 'Entry Name'),
            'gene': entry(acc, 'Gene Names (primary)'),
            'protein_name': recommended_name(entry(acc, 'Protein names')),
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
            'dimethyl_type_neutral_loss': join(ordered({t for x in group['_rme2_type'] for t in x.split('; ') if t}, list(rc.CODES['dimethyl_type']))),
            'flags': join(flags),
            'shared_with_paralogs': join(paralogs),
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
    at_before = records.groupby('_previous_at').groups
    reported_gene = records.groupby('previous_site')['reported_gene'].agg(lambda g: next((x for x in g if x), ''))
    out = []
    for site in sorted(now | set(before), key=split_key):
        acc, pos = split_key(site)
        in_before, in_now = site in before, site in now
        if in_before and in_now:
            # The records that moved onto or off the site; with none of the previous ones left, its records were replaced.
            change = 'kept'
            gained = records.loc[at_now.get(site, [])]
            gained = gained[gained['_previous_at'] != site]
            lost = records.loc[at_before.get(site, [])]
            lost = lost[site_of[lost.index] != site]
            involved = pd.concat([gained, lost])
            kept_records = len(at_before.get(site, [])) - len(lost)
            codes = (['reference_renumbered'] if before[site] != site else []) + sorted(set(involved['reason_code']) - {'', 'reference_renumbered'})
            if gained['change_since_previous_version'].eq('new_record').any():
                codes = sorted(set(codes) | {'new_record'})
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
            if involved['change_since_previous_version'].eq('new_record').any():
                codes = sorted(set(codes) | {'new_record'})
        if change == 'kept':
            elsewhere = sorted((set(site_of[lost.index]) | set(gained['previous_site'])) - {'', site, before[site]}, key=split_key)
            only_new = lost.empty and len(gained) and gained['change_since_previous_version'].eq('new_record').all()
            records_change = ('same' if involved.empty else 'replaced' if kept_records == 0 else 'records_added' if only_new
                              else 'changed')
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
        if group.empty:
            continue   # a study this build does not read (the study tables, when rebuilding v1.1.0)
        supported = per_site[per_site.map(lambda s: sid in s)]
        unplaced = group[group['placed_site'] == '']
        several = int(unplaced['placement_code'].isin(rc.PLACEMENT_SEVERAL).sum())
        other = int(unplaced['placement_code'].isin(rc.PLACEMENT_OTHER_PROTEIN).sum())
        not_localised = int(unplaced['placement_code'].isin(rc.PLACEMENT_NOT_LOCALISED).sum())
        out.append({'study_id': sid, 'label': label, 'source_table': join(ordered(group['source_table'], list(rc.SOURCE_TABLES))),
                    'citation': citation, 'doi': doi, 'licence': rc.STUDY_LICENCE.get(sid, rc.MARON_S5),
                    'n_records': len(group), 'n_records_placed': int((group['placed_site'] != '').sum()),
                    'n_sites': len(supported), 'n_sites_only_here': int((supported.map(len) == 1).sum()),
                    'n_sites_in_previous_version': int(placed.loc[supported.index, 'in_previous_version'].sum()),
                    'n_records_several_sites': several, 'n_records_other_protein': other, 'n_records_not_localised': not_localised,
                    'n_records_no_site': len(unplaced) - several - other - not_localised,
                    'note': rc.STUDY_NOTE.get(sid, ''), 'quality': rc.STUDY_QUALITY.get(sid, '')})
    unknown = set(records['study']) - set(rc.STUDY_ORDER)
    if unknown:
        raise SystemExit(f'studies missing from scripts/release_codes.py: {sorted(unknown)}')
    return pd.DataFrame(out)


def evidence_cores(record: pd.Series) -> list[tuple[str, int]]:
    """Each of the record's sequences (a merged record can carry several) as evidence_core reads it."""
    return [evidence_core(record, token) for token in clean(record['evidence_sequence']).split('; ')]


def evidence_core(record: pd.Series, token: str | None = None) -> tuple[str, int]:
    """The record's first sequence (or token) without marks or padding, and the 0-based offset of the site it names:
    the marked arginine of a marked peptide, the centre of a window, the arginine of a peptide with one; -1 when it
    names none."""
    if token is None:
        token = clean(record['evidence_sequence']).split('; ')[0]
    kind = record['evidence_kind']
    if kind == 'marked_peptide':
        plain, marked = remap.parse_modified_peptide(token)
        return plain, marked[0] - 1 if len(marked) == 1 else -1
    if kind.startswith('window_'):
        core, lead, _ = remap.split_window_padding(token)
        return core, int(kind.split('_')[1]) // 2 - lead
    if kind == 'peptide':
        arginines = [k for k, aa in enumerate(token) if aa == 'R']
        return token, arginines[0] if len(arginines) == 1 else -1
    return '', -1


def traced(sequence: str, offset: int, reference: str, start: int, flank: int = 6) -> list[str]:
    """Two aligned lines: the evidence with its site in brackets, and the reference around it (start is 1-based)."""
    site = start + offset
    lo, hi = max(1, start - flank), min(len(reference), start + len(sequence) - 1 + flank)
    ref = reference[lo - 1:site - 1] + '[' + reference[site - 1] + ']' + reference[site:hi]
    ev = ' ' * (start - lo) + sequence[:offset] + '[' + sequence[offset] + ']' + sequence[offset + 1:]
    return [f'evidence    {" " * len(str(lo))} {ev}', f'reference   {lo} {ref} {hi}']


def evidence_words(kind: str) -> str:
    return {'marked_peptide': 'marked peptide', 'peptide': 'peptide', 'position_only': 'position'}.get(
        kind, kind.replace('window_', 'window of ') + ' residues')


def placed_example(r: pd.Series, reference: dict, proteome) -> list[str]:
    """A placed record: what the source reported, where it is placed and why, and its evidence on the reference."""
    acc, pos = split_key(r['placed_site'])
    seq, off = evidence_core(r)
    lines = [f"record {r['record_id']}: reported {r['reported_accession']} at R{r['reported_position']}"
             + (f" ({r['reported_detail']})" if r['reported_detail'] else '')]
    if not seq:
        lo, hi = max(1, pos - 7), min(len(reference[acc]), pos + 7)
        return lines + [f"placed on {r['placed_site']} ({r['placement_code']}): no sequence, the reported residue is an arginine",
                        f'reference   {lo} {reference[acc][lo - 1:pos - 1]}[{reference[acc][pos - 1]}]{reference[acc][pos:hi]} {hi}']
    for s, o in evidence_cores(r):   # trace the sequence that names the placed site, not just the first
        if s and any(a == acc and (pos - start == o if o >= 0 else 0 <= pos - start < len(s)) for a, start in proteome.find(s)):
            seq, off = s, o
            break
    hits = proteome.find(seq)
    if off < 0:   # a peptide with several arginines: the occurrence that covers the site
        off = next((pos - start for a, start in hits if a == acc and 0 <= pos - start < len(seq)), -1)
    lines.append(f"placed on {r['placed_site']} ({r['placement_code']}; evidence check: {r['evidence_check']}); the "
                 f"{evidence_words(r['evidence_kind'])} " + ('occurs once' if len(hits) == 1 else f'occurs {len(hits)} times')
                 + ' in the reviewed proteome')
    return lines + (traced(seq, off, reference[acc], pos - off) if off >= 0 else [])


def code_reason(code: str) -> str:
    """A placement code's meaning without its leading 'not placed' words."""
    return re.sub(r'^not placed( \([^)]*\))?: ', '', rc.CODES['placement'].get(code, ''))


def unplaced_example(r: pd.Series, proteome) -> list[str]:
    """A record left out: its evidence, its placement code and meaning, and where the evidence occurs."""
    seq, off = evidence_core(r)
    hits = proteome.find(seq) if seq else []
    where = sorted({key(a, start + off) for a, start in hits}) if off >= 0 else sorted(f'{a}:{start}-{start + len(seq) - 1}' for a, start in hits)
    lines = [f"record {r['record_id']}: reported {r['reported_accession']} at R{r['reported_position']}"
             + (f", {evidence_words(r['evidence_kind'])} {r['evidence_sequence']}" if seq else ', no sequence'),
             f"not placed ({r['placement_code']}): {code_reason(r['placement_code'])}"]
    if seq and not where:
        lines.append('the evidence occurs nowhere in the reviewed proteome')
    elif seq:
        lines.append(f"the evidence {'names' if off >= 0 else 'occurs at'} {len(where)} "
                     f"{'site' if off >= 0 else 'place'}{'s' if len(where) != 1 else ''} in the reviewed proteome: "
                     + ', '.join(where[:6]) + (', ...' if len(where) > 6 else ''))
    return lines


def study_examples(records: pd.DataFrame, reference: dict) -> dict:
    """For each study, records traced from the source to the reviewed proteome, as text lines for STUDIES.md: a placed
    record (for a study read from its own table the one merged from most source rows, otherwise the first whose evidence check
    agrees), and records left out (for an added study the first whose evidence names several sites and the first whose
    one site is in a protein the record does not name, otherwise the first left out)."""
    proteome = remap.ProteomeIndex(reference)
    staged = set(rc.STUDY_TABLES.values())
    out = {}
    for sid in rc.STUDY_ORDER:
        g = records[records['study'] == sid]
        if g.empty:
            continue
        blocks = []
        placed = g[g['placed_site'] != '']
        if sid in staged:
            rows = pd.to_numeric(g['reported_detail'].str.extract(r'(\d+) source rows?\b')[0], errors='coerce').fillna(1)
            placed = placed.loc[rows[placed.index].sort_values(ascending=False, kind='stable').index]
        else:
            placed = pd.concat([placed[(placed['evidence_check'] == 'agrees') & (placed['evidence_kind'] != 'position_only')], placed])
        if len(placed):
            blocks.append(placed_example(placed.iloc[0], reference, proteome))
        left = g[(g['placed_site'] == '') & (g['pending_site'] == '')]
        codes = ['unresolved_study_not_unique', 'unresolved_study_other_entry'] if sid in staged else [None]
        for code in codes:
            pick = left if code is None else left[left['placement_code'] == code]
            if len(pick):
                blocks.append(unplaced_example(pick.iloc[0], proteome))
        out[sid] = blocks
    return out


def validation(root: Path, check: pd.DataFrame) -> dict:
    """The outcome of the build's two gates (hpc/build_release.sh stops on either)."""
    audit = root / 'qc' / 'remap_evidence_check' / 'union_site_residue_audit.tsv'
    failures = root / 'qc' / 'remap_evidence_check' / 'remap_evidence_check_failures.tsv'
    missing = [p.name for p in (audit, failures) if not p.exists()]
    if missing:   # a release records the gates' outcome; without their output it cannot say they passed
        raise SystemExit(f'gate output missing from {root / "qc" / "remap_evidence_check"}: {missing} (run the evidence check and audit)')
    return {'evidence_check': {rc.VERDICT_CODE[k]: int(n) for k, n in check['verdict'].value_counts().items()},
            'rows_against_evidence_not_accepted': len(read_tsv(failures)),
            'union_sites_not_arginine': len(read_tsv(audit))}


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
    # the RNA and DNA each site touches in structures (scripts/37 from v1.3.0)
    na = read_tsv(folder / 'nucleic_acid_contacts.tsv') if (folder / 'nucleic_acid_contacts.tsv').exists() else None
    if na is not None and set(na['site_id']) - set(analysed['site_id']):
        raise SystemExit(f'{folder}: nucleic-acid contacts for sites that are not analysed in the release')
    regions = read_tsv(folder / 'nucleic_acid_regions.tsv') if (folder / 'nucleic_acid_regions.tsv').exists() else None
    beyond = {name: read_tsv(folder / f'{name}.tsv') for name in ('arginine_nucleic_acid_distances', 'nucleic_acid_proximal_sites', 'nucleic_acid_near_sites', 'rna_bind_n_seq', 'dna_motifs', 'eclip')
              if (folder / f'{name}.tsv').exists()}
    for name in ('nucleic_acid_proximal_sites', 'nucleic_acid_near_sites'):
        if name in beyond and set(beyond[name]['site_id']) - set(analysed['site_id']):
            raise SystemExit(f'{folder}: {name} has sites that are not analysed in the release')
    return {'sites': ev, 'residues': t['residue_evidence'], 'partners': t['protein_partners'], 'complexes': t['complexes'],
            'statistics': t['interaction_statistics'], 'proteins': t['protein_interaction_summary'], 'status': status, 'evidence': evidence,
            'nucleic_acid': na, 'nucleic_acid_regions': regions, 'beyond': beyond}


def bound_columns(ix: dict) -> dict:
    """site_id -> the site table's bound_* values, from the RNA and DNA in direct contact with each site's side chain in
    structures (nucleic_acid_contacts), within 8 Å of it where there is none (nucleic_acid_proximal_sites), and in contact
    via a residue at i ± 1-5 where there is neither (nucleic_acid_near_sites). The bound sequence shown is the best
    supported: direct, then proximal, then via a flanking residue; then the most PDB entries, then the type's code order;
    its marked nucleotides are those of every entry showing it."""
    def hits(marked):
        out, k, i = set(), 0, 0
        while i < len(marked):
            if marked[i] == '[':
                out.add(k)
                i += 3
            else:
                i += 1
            k += 1
        return out

    groups = {}   # (site, window, type, tier) -> {tier, contact, offset, distance, entries, offsets, names}
    na = ix.get('nucleic_acid')
    for r in (na.itertuples() if na is not None else []):
        g = groups.setdefault((r.site_id, r.window, r.nucleic_acid_type, 0), {'tier': 0, 'contact': ('dna' if r.nucleic_acid == 'dna' else 'rna') + '_at', 'offset': 0,
                                                                          'distance': float(r.closest_distance_a), 'entries': set(), 'offsets': set(), 'names': Counter()})
        g['entries'].add(r.pdb_id)
        g['distance'] = min(g['distance'], float(r.closest_distance_a))
        g['offsets'] |= hits(r.window_marked)
        g['names'][r.molecule_name] += 1
    prox = ix['beyond'].get('nucleic_acid_proximal_sites')
    for r in (prox.itertuples() if prox is not None else []):
        # within 4 Å in the entry read: a direct contact PDBe-KB does not list at the interface
        tier = 1 if float(r.distance_a) > 4 else 0
        g = groups.setdefault((r.site_id, r.window, r.nucleic_acid_type, 'prox', tier), {'tier': tier, 'contact': r.nucleic_acid + ('_proximal' if tier else '_at'), 'offset': 0,
                                                                          'distance': float(r.distance_a), 'entries': 0, 'offsets': set(), 'names': Counter()})
        g['entries'] += int(r.n_entries)
        g['distance'] = min(g['distance'], float(r.distance_a))
        g['offsets'] |= hits(r.window_marked)
        g['names'][r.molecule_name] += int(r.n_entries)
    near = ix['beyond'].get('nucleic_acid_near_sites')
    for r in (near.itertuples() if near is not None else []):
        g = groups.setdefault((r.site_id, r.window, r.nucleic_acid_type, 2), {'tier': 2, 'contact': r.nucleic_acid + '_near', 'offset': int(r.residue_offset),
                                                                          'distance': None, 'entries': 0, 'offsets': set(), 'names': Counter()})
        g['entries'] += int(r.n_entries)
        g['offset'] = min(g['offset'], int(r.residue_offset), key=lambda v: (abs(v), v))
        g['offsets'] |= hits(r.window_marked)
        g['names'][r.molecule_name] += int(r.n_entries)
    by_site = defaultdict(list)
    for (site, window, ntype, *_), g in groups.items():
        g['n'] = len(g['entries']) if isinstance(g['entries'], set) else g['entries']
        by_site[site].append((window, ntype, g))
    types = list(rc.CODES['nucleic_acid_type'])
    out = {}
    for site, rows in by_site.items():
        window, ntype, g = min(rows, key=lambda x: (x[2]['tier'], -x[2]['n'], types.index(x[1]), x[0]))
        out[site] = {
            'bound_nucleic_acid': join(ordered({x[2]['contact'] for x in rows}, list(rc.CODES['bound_contact']))),
            'bound_nucleic_acid_types': join(ordered({x[1] for x in rows}, types)),
            'bound_sequence': ''.join(f'[{c}]' if k in g['offsets'] else c for k, c in enumerate(window)),
            'bound_sequence_offset': g['offset'],
            'bound_sequence_distance_a': g['distance'] if g['distance'] is not None else '',
            'bound_sequence_type': ntype,
            'bound_sequence_molecule': g['names'].most_common(1)[0][0],
            'bound_sequence_entries': g['n'],
        }
    return out


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
        if kind.startswith('list'):
            # items separated by '; ', none empty and none holding ';' itself
            bad |= out[name].map(lambda v: bool(v) and any(not x or ';' in x for x in v.split('; ')))
        else:
            bad |= out[name].str.contains(';', regex=False)
        if bad.any():
            raise SystemExit(f'{path.name}: column {name} has a separator inside a value: {out.loc[bad, name].iloc[0]!r}')
    release_docs.write_plain_tsv(out, path)


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
    ap.add_argument('--analyses-root', default='', help='where the analyses ran (default --results-root)')
    ap.add_argument('--with-lower-confidence', action='store_true',
                    help='write the view with lower-confidence sites: the site table holds both tiers (column tier), and the '
                         'proteins and statistics come from --analyses-root, the analyses of both (the website reads it)')
    ap.add_argument('--statistics-with-lower-confidence', default='',
                    help="the statistics table of that view, kept in the release as statistics_with_lower_confidence")
    args = ap.parse_args()
    root, inputs, outdir, prev_dir = Path(args.results_root), Path(args.inputs), Path(args.outdir), Path(args.previous)
    aroot = Path(args.analyses_root) if args.analyses_root else root
    version, previous = args.version, prev_dir.name

    reference = remap.parse_fasta(inputs / INPUT['fasta'])
    meta = read_tsv(inputs / INPUT['entries']).set_index('Entry')

    def entry(acc: str, field: str) -> str:
        return clean(meta[field].get(acc, ''))

    rows = read_tsv(root / 'remapped' / 'human_arg_methyl_source_rows_with_remap.tsv')
    typed = pd.read_csv(root / 'remapped' / 'human_arg_methyl_source_rows_with_remap.tsv', sep='\t', low_memory=False)
    typed = remap.annotate_window_lengths(typed) if 'seq_window_window_length' not in typed else typed
    check = read_tsv(root / 'qc' / 'remap_evidence_check' / 'remap_evidence_check_rows.tsv')
    # the evidence check writes one row per record, in the records' order: joined by position, so check it
    if len(check) != len(rows) or not (check['reported_accession'].values == rows['substrate_UniProtAC_original'].replace('nan', '').values).all() \
            or not (check['source_table'].values == rows['source_table'].values).all():
        raise SystemExit('the evidence check rows do not line up with the records (rerun scripts/qc_remap_evidence_check.py)')
    union = read_tsv(root / 'integrated' / 'human_arg_methyl_union_dedup_by_site.tsv')
    pending_table = read_tsv(root / 'integrated' / 'pending_one_sided_window_sites.tsv')
    pending = {key(a, p) for a, p in zip(pending_table['canonical_UniProtAC'], pending_table['position'])}
    pending_of = {int(i): key(a, p) for a, p, ids in zip(pending_table['canonical_UniProtAC'], pending_table['position'], pending_table['record_ids'])
                  for i in ids.split('; ')}
    allow = pd.read_csv(args.allow, sep='\t', dtype=str, comment='#').fillna('')
    allowed = {(r.source_table, r.reported_accession, r.reported_position, r.placed_site): r.reason for r in allow.itertuples()}

    prev_sites = read_tsv(prev_dir / f'human_methylarginine_sites_{previous}.tsv')
    prev_records = read_tsv(prev_dir / 'supporting' / f'source_records_{previous}.tsv')
    prev_info = json.loads((prev_dir / 'release.json').read_text()) if (prev_dir / 'release.json').exists() else {}
    if prev_info.get('reference', {}).get('fasta_md5') == hashlib.md5((inputs / INPUT['fasta']).read_bytes()).hexdigest():
        shift, broken = {}, set()   # the same reference sequences: the previous numbering is this one
    else:
        shift, broken = reference_lift(inputs)
    stuck = sorted(s for s in prev_sites['site_id'] if split_key(s)[0] in broken)
    if stuck:
        raise SystemExit(f'previous sites on entries whose numbering cannot be carried over: {stuck}')

    # This version reads the previous version's source records in the same order, as reported (before the corrections
    # of source_correction), and after them the records read for the first time.
    n_prev = len(prev_records)
    head = rows.iloc[:n_prev]
    same = (len(rows) >= n_prev and len(check) == len(rows)
            and (head['source_table'].map(rc.SOURCE_TABLE_CODE).values == prev_records['source_table'].values).all()
            and (head['substrate_UniProtAC_original'].map(clean).values == prev_records['reported_accession'].values).all()
            and (head['position'].map(clean).values == prev_records['reported_position'].values).all())
    if same and 'evidence_sequence' in prev_records:
        # a release of this format: the evidence as this run reads it
        evidence = [evidence_of(typed.loc[i]) for i in range(n_prev)]
        same = ([k for k, _ in evidence] == list(prev_records['evidence_kind'])
                and [t for _, t in evidence] == list(prev_records['evidence_sequence']))
    elif same:
        # v1.0.0: the peptide and window as the source tables give them
        reported_peptide = head['peptide_sequence_reported'].map(clean).where(head['peptide_sequence_reported'].map(clean) != '', head['peptide_sequence'].map(clean)) \
            if 'peptide_sequence_reported' in head else head['peptide_sequence'].map(clean)
        same = ((reported_peptide.values == prev_records['peptide_sequence'].values).all()
                and (head['seq_window'].map(clean).values == prev_records['seq_window'].values).all())
    if not same:
        raise SystemExit('this run does not hold the previous release\'s source records in the same order')

    records = build_records(rows, typed, check, allowed, prev_records, pending_of, shift, reference)
    union_keys = {key(a, p) for a, p in zip(union['canonical_UniProtAC'], union['corrected_position'])}
    if union_keys != set(records['placed_site']) - {''}:
        raise SystemExit('the placed records and the union disagree')
    previous_lifted = {lifted(x, shift) for x in prev_sites['site_id']}
    sites = build_sites(records, union, pending, reference, entry, previous_lifted)
    site_ids = set(sites['site_id'])
    # The lower-confidence sites go through the same annotation as the sites, then apart.
    members, published = earlier_sites(prev_dir, previous, prev_sites, records, shift)
    lc_records, lc_tests = lower_confidence_sites(records, site_ids, reference, members)
    lc_sites = build_sites(lc_records, union.iloc[:0], set(), reference, entry, previous_lifted, placed_here=False) if len(lc_records) else sites.iloc[:0]
    lc_tests = lc_tests[lc_tests['site_id'].isin(set(lc_sites['site_id']))]
    n_main = len(sites)
    sites = pd.concat([sites, lc_sites], ignore_index=True)
    all_ids = set(sites['site_id'])

    detection, inhibitors, _ = maron_tables(inputs, reference, all_ids)
    maron_counts = maron_tables(inputs, reference, site_ids)[2]
    s2 = detection[detection['maron_table'] == 'S2'].set_index('site_id')['methyl_states']
    s3 = detection[detection['maron_table'] == 'S3'].groupby('site_id')['antibody'].agg(lambda a: join(ordered(a, rc.ANTIBODY_ORDER)))
    sites['maron2021_s2_states'] = sites['site_id'].map(s2).fillna('')
    sites['maron2021_s3_antibodies'] = sites['site_id'].map(s3).fillna('')
    sites = pd.concat([sites, state_calls(sites, pd.concat([records, lc_records], ignore_index=True), inputs)], axis=1)
    for drug in rc.CODES['inhibitor']:
        calls = inhibitors[inhibitors['inhibitor'] == drug].groupby('site_id')['call'].agg(list)
        sites[f'maron2021_{drug}_call'] = sites['site_id'].map(lambda s: combine_calls(calls.get(s, [])))
    carm1 = carm1_table(inputs, all_ids)
    carm1_calls = carm1.groupby('site_id')['call'].agg(list)
    sites['carm1_ko_call'] = sites['site_id'].map(lambda s: combine_calls(carm1_calls.get(s, [])))
    sites = pd.concat([sites, context_annotations(sites, inputs)], axis=1)
    ix = interaction_tables(Path(args.interactions), sites if args.with_lower_confidence else sites.iloc[:n_main]) if args.interactions else None
    analysed = ~sites['flags'].str.contains('one_sided_window', regex=False)   # sites kept pending are left out
    sites['interaction_status'] = sites['site_id'].map(ix['status']).fillna('none').where(analysed, '') if ix else ''
    sites['interaction_evidence'] = sites['site_id'].map(ix['evidence']).fillna('') if ix else ''
    bound = bound_columns(ix) if ix else {}
    for c in rc.BOUND_SITE_COLUMNS:
        sites[c] = sites['site_id'].map(lambda s: bound[s][c] if s in bound else '')
    sites = sites[[c for c, *_ in rc.SITE_COLUMNS]]
    lc_sites = sites.iloc[n_main:].merge(lc_tests, on='site_id', how='left')
    lc_sites = lc_sites[[c for c, *_ in rc.LOWER_CONFIDENCE_SITE_COLUMNS]].reset_index(drop=True)
    published = published_sites(published, site_ids, lc_tests)
    if args.with_lower_confidence:
        # both tiers in one site table, the sites' tests all passed
        tier = pd.DataFrame({'tier': ['site'] * n_main + ['lower_confidence'] * len(lc_sites)})
        main_tests = pd.DataFrame({'site_id': sites['site_id'].iloc[:n_main], 'in_analyses_with_lower_confidence': True,
                                   **dict.fromkeys(TEST_FIELDS, 'pass'), 'tests_failed': '', 'reasons': '', 'in_earlier_versions': ''})
        tier_tests = pd.concat([main_tests, lc_tests.set_index('site_id').loc[sites['site_id'].iloc[n_main:]].reset_index()], ignore_index=True)
        sites = pd.concat([sites.reset_index(drop=True), tier, tier_tests.drop(columns='site_id')], axis=1)
        placed = pd.concat([records, lc_records], ignore_index=True)   # each record at every lower-confidence site it supports
    else:
        placed = records
        sites = sites.iloc[:n_main]
        detection, inhibitors, carm1 = (t[t['site_id'].isin(site_ids)] for t in (detection, inhibitors, carm1))

    main_sites = sites.iloc[:n_main]
    changes = build_changes(records, main_sites, list(prev_sites['site_id']), shift, entry)
    studies = build_studies(records, main_sites)
    proteins, extra = release_tables.proteins_table(aroot, inputs, sites, placed, entry, reference)
    statistics = release_tables.statistics(aroot, proteins)
    for column, *_ in rc.STATISTICS_COLUMNS:   # a column only some analyses fill (the 95% CI of an odds ratio)
        if column not in statistics.columns:
            statistics[column] = None
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
    disease_terms = release_tables.disease_terms_table(aroot, proteins)

    files = {
        'sites': (f'human_methylarginine_sites_{version}.tsv', sites, rc.VIEW_SITE_COLUMNS if args.with_lower_confidence else rc.SITE_COLUMNS),
        'lower_confidence_sites': (f'human_methylarginine_lower_confidence_sites_{version}.tsv', lc_sites, rc.LOWER_CONFIDENCE_SITE_COLUMNS),
        'studies': (f'studies_{version}.tsv', studies, rc.STUDY_COLUMNS),
        'records': (f'supporting/source_records_{version}.tsv', records, rc.RECORD_COLUMNS),
        'lower_confidence': (f'supporting/lower_confidence_evidence_{version}.tsv', lower_confidence(records, site_ids), rc.LOWER_CONFIDENCE_COLUMNS),
        'published_sites': (f'supporting/published_sites_{version}.tsv', published, rc.PUBLISHED_SITE_COLUMNS),
        'changes': (f'supporting/site_changes_{previous}_to_{version}.tsv', changes, rc.CHANGE_COLUMNS),
        'maron_detection': (f'supporting/maron2021_detection_{version}.tsv', detection, rc.MARON_DETECTION_COLUMNS),
        'maron_inhibitors': (f'supporting/maron2021_inhibitors_{version}.tsv', inhibitors, rc.MARON_INHIBITOR_COLUMNS),
        'carm1_ko': (f'supporting/carm1_ko_{version}.tsv', carm1, rc.CARM1_COLUMNS),
        'proteins': (f'human_methylarginine_proteins_{version}.tsv', proteins, rc.PROTEIN_COLUMNS),
        'statistics': (f'statistics_{version}.tsv', statistics, rc.STATISTICS_COLUMNS),
        'disease_terms': (f'supporting/disease_terms_{version}.tsv', disease_terms, rc.DISEASE_TERM_COLUMNS),
        'disease_terms_tested': (f'supporting/disease_terms_tested_{version}.tsv', release_tables.disease_terms_tested_table(aroot, proteins),
                                 rc.DISEASE_TERM_TESTED_COLUMNS),
        'disease_ontology': (f'supporting/disease_ontology_{version}.tsv', release_tables.disease_ontology_table(aroot, proteins), rc.DISEASE_ONTOLOGY_COLUMNS),
        'density_windows': (f'supporting/density_windows_{version}.tsv', extra['density_windows'], rc.DENSITY_COLUMNS),
        'hpa_locations': (f'supporting/hpa_locations_{version}.tsv', extra['hpa'], rc.HPA_COLUMNS),
        'other_ptm_sites': (f'supporting/other_ptm_sites_{version}.tsv', extra['other_ptm_sites'], rc.OTHER_PTM_COLUMNS),
        'condensates': (f'supporting/condensates_{version}.tsv', release_tables.condensates_table(aroot, set(sites['uniprot_accession'])), rc.CONDENSATE_COLUMNS),
        'disease_families': (f'supporting/disease_families_{version}.tsv', release_tables.disease_families_table(aroot), rc.DISEASE_FAMILY_COLUMNS),
    }
    if args.statistics_with_lower_confidence:
        files['statistics_with_lower_confidence'] = (f'statistics_with_lower_confidence_{version}.tsv',
                                                     read_tsv(Path(args.statistics_with_lower_confidence)), rc.STATISTICS_COLUMNS)
    if ix:
        files.update({
            'site_interactions': (f'supporting/site_interactions_{version}.tsv', ix['sites'], rc.SITE_INTERACTION_COLUMNS),
            'residue_evidence': (f'supporting/residue_evidence_{version}.tsv', ix['residues'], rc.RESIDUE_EVIDENCE_COLUMNS),
            'protein_partners': (f'supporting/protein_partners_{version}.tsv', ix['partners'], rc.PARTNER_COLUMNS),
            'complexes': (f'supporting/complexes_{version}.tsv', ix['complexes'], rc.COMPLEX_COLUMNS),
        })
        if ix['nucleic_acid'] is not None:
            files['nucleic_acid_contacts'] = (f'supporting/nucleic_acid_contacts_{version}.tsv', ix['nucleic_acid'], rc.NUCLEIC_ACID_CONTACT_COLUMNS)
        if ix['nucleic_acid_regions'] is not None:
            files['nucleic_acid_regions'] = (f'supporting/nucleic_acid_regions_{version}.tsv', ix['nucleic_acid_regions'], rc.NUCLEIC_ACID_REGION_COLUMNS)
        for name, columns in (('arginine_nucleic_acid_distances', rc.ARGININE_DISTANCE_COLUMNS), ('nucleic_acid_proximal_sites', rc.NUCLEIC_ACID_PROXIMAL_COLUMNS),
                              ('nucleic_acid_near_sites', rc.NUCLEIC_ACID_NEAR_COLUMNS), ('rna_bind_n_seq', rc.RNA_BIND_N_SEQ_COLUMNS),
                              ('dna_motifs', rc.DNA_MOTIF_COLUMNS), ('eclip', rc.ECLIP_COLUMNS)):
            if name in ix['beyond']:
                files[name] = (f'supporting/{name}_{version}.tsv', ix['beyond'][name], columns)
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

    counts = release_docs.counts(main_sites, records, changes, studies, maron_counts, prev_sites)
    counts['lower_confidence_sites'] = len(lc_sites)
    counts['lower_confidence_sites_analysed'] = int(lc_sites['in_analyses_with_lower_confidence'].astype(bool).sum())
    counts.update(peptide_start_agreement(rows, typed, reference, root / 'parsed_maron' / 'maron_peptide_coordinates.tsv'))
    context = {'version': version, 'previous': previous, 'files': files, 'counts': counts, 'inputs': inputs, 'root': root,
               'uniprot_release': json.loads((inputs / INPUT['uniprot_release']).read_text()), 'records': records,
               'sites': sites, 'changes': changes, 'studies': studies, 'release_date': args.date, 'previous_sites': prev_sites,
               'validation': validation(root, check), 'enrichment_tools': enrichment_tools(aroot),
               'study_examples': study_examples(records, reference)}
    release_docs.write_all(outdir, context)
    print(json.dumps({k: v for k, v in counts.items() if not isinstance(v, dict)}))


def enrichment_tools(root: Path) -> dict:
    """The R and package versions every enrichment run of the analyses recorded (scripts/enrich_sets.R); they must agree."""
    found = {}
    for path in sorted(root.glob('**/enricher/results/versions.tsv')):
        t = pd.read_csv(path, sep='\t', dtype=str)
        found[str(path.relative_to(root))] = dict(zip(t['package'], t['version']))
    versions = {tuple(sorted(v.items())) for v in found.values()}
    if len(versions) > 1:
        raise SystemExit(f'the enrichment runs used different tool versions: {found}')
    return dict(versions.pop()) if versions else {}


if __name__ == '__main__':
    main()

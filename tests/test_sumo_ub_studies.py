"""Tests of the Hendriks 2017 (SUMO-2) and Akimov 2018 (ubiquitination) staging (scripts/60, 61, window_sites.py) and of
their rows in releases/v1.5.0.

Run from the repository root: python -m pytest -q tests/test_sumo_ub_studies.py
"""
import csv
import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import window_sites as ws  # noqa: E402


def load(name, alias):
    spec = importlib.util.spec_from_file_location(alias, ROOT / 'scripts' / name)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


h = load('60_stage_hendriks2017_sumo_sites.py', 'stage60')
a = load('61_stage_akimov2018_ubiquitination_sites.py', 'stage61')
build = load('56_build_context_marks.py', 'build56')
MARKS = ROOT / 'inputs' / 'marks'
RELEASE = ROOT / 'releases' / 'v1.5.0'
REFERENCE = {'P00001': 'MAAAAAAAKCCCCCCCDDDDDDDKEEEEEEEFFF', 'P00002': 'MGGGGGGGKHHHHHHH', 'P00003': 'MKAAAAAAAKAAAAAAAK'}


def window(seq, pos, flank=7):
    padded = '_' * flank + seq + '_' * flank
    return padded[pos - 1:pos - 1 + 2 * flank + 1]


def test_placement_kept_secondary_relocated_rejected():
    motif = window(REFERENCE['P00001'], 9)
    assert ws.place([('P00001', 9)], [motif], REFERENCE, {})[:3] == ('kept', 'P00001', 9)
    assert ws.place([('Q99999', 9)], [motif], REFERENCE, {'Q99999': 'P00001'})[:3] == ('secondary_accession', 'P00001', 9)
    assert ws.place([('P00001-3', 2)], [motif], REFERENCE, {})[:3] == ('relocated', 'P00001', 9)          # isoform: by window
    assert ws.place([('P00001', 5)], [motif], REFERENCE, {})[:3] == ('relocated', 'P00001', 9)            # moved position
    assert ws.place([('A0A000', 9), ('P00001', 9)], [motif], REFERENCE, {})[:3] == ('kept', 'P00001', 9)  # later group member
    assert ws.place([('A0A000', 9)], [motif], REFERENCE, {})[0] == 'no_canonical_accession'
    assert ws.place([('P00002', 9)], [motif], REFERENCE, {})[0] == 'window_not_in_canonical'
    assert ws.place([('P00001', 9)], [None], REFERENCE, {})[0] == 'bad_window'
    # N-terminal padding must lie before the start
    nterm = window(REFERENCE['P00003'], 2)
    assert nterm.startswith('______') and ws.place([('P00003', 2)], [nterm], REFERENCE, {})[:3] == ('kept', 'P00003', 2)
    assert ws.hits(REFERENCE['P00003'], nterm) == [2]
    # a core found twice in the sequence and nowhere once: not unique
    assert ws.place([('P00004', 9)], ['AAAAAAAKAAAAAAA'], {'P00004': 'AAAAAAAKAAAAAAAKAAAAAAA'}, {})[0] == 'window_not_unique'


def test_core_and_group():
    w51 = 'X' * 25 + 'K' + 'Y' * 25
    assert ws.core(w51, 25) == 'XXXXXXXKYYYYYYY' and ws.core('X' * 51, 25) is None
    assert ws.group('A;B-2;C', '5;6;7', 'B-2', '6') == [('B-2', 6), ('A', 5), ('C', 7)]


def test_hendriks_toy_rows_merge_on_one_lysine():
    seq = REFERENCE['P00001']
    def row(protein, position, proteins, positions, prob, hela='+', u2os=None):
        r = {'Protein': protein, 'Position': position, 'Proteins': proteins, 'Positions within proteins': positions,
             h.WINDOW: '_' * 18 + window(seq, 9) + '_' * 18, 'Localization prob': prob, 'Amino acid': 'K'}
        for c, cols in h.DETECTION.items():
            for x in cols:
                r[x] = None
        r[h.DETECTION['HeLa'][0]], r[h.DETECTION['U2OS'][2]] = hela, u2os
        return r
    rows = [row('P00001', 9, 'P00001', '9', 0.9), row('P00001-2', 3, 'P00001-2', '3', 1.0, hela=None, u2os='+'),
            row('A0A000', 9, 'A0A000', '9', 1.0)]
    sites, bad, counts = h.stage(rows, REFERENCE, {})
    assert counts == {'kept': 1, 'relocated': 1, 'no_canonical_accession': 1}
    assert len(sites) == 1 and list(sites.columns) == h.COLUMNS and list(bad.columns) == h.REJECTED_COLUMNS
    r = sites.iloc[0]
    assert (r['position'], r['hendriks_localization_probability'], r['hendriks_cell_lines'], r['placement']) == (9, '1', 'HeLa;U2OS', 'kept;relocated')
    assert (r['state'], r['mark'], r['source'], r['evidence_class']) == ('sumo', 'sumoylation', 'Hendriks 2017', 'experimental')


@pytest.mark.parametrize('mod,mark,state,prob_col', [(h, 'sumoylation', 'sumo', 'hendriks_localization_probability'),
                                                     (a, 'ubiquitination', 'ub', 'akimov_localization_probability')])
def test_staged_tables_summaries_and_release(mod, mark, state, prob_col):
    from common import parse_fasta
    reference = parse_fasta(ROOT / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03.fasta')
    staged = pd.read_csv(MARKS / mod.OUT_SITES, sep='\t', dtype=str, keep_default_na=False)
    assert list(staged.columns) == mod.COLUMNS and not staged.duplicated(['canonical_UniProtAC', 'position']).any()
    assert set(staged['state']) == {state} and set(staged['residue']) == {'K'} and (staged[prob_col] != '').all()
    for r in staged.itertuples():
        assert reference[r.canonical_UniProtAC][int(r.position) - 1] == 'K'
    s = json.loads((MARKS / mod.OUT_SUMMARY).read_text(encoding='utf-8'))
    assert s['input_rows'] == sum(s['row_status'].values()) == s['rows_placed'] + sum(s['rejected_by_reason'].values())
    assert s['validated_sites'] == len(staged) == s['rows_placed'] - s['rows_merged_into_another_site']
    with open(MARKS / mod.OUT_REJECTED, encoding='utf-8') as fh:
        assert sum(1 for _ in fh) - 1 == sum(s['rejected_by_reason'].values())
    raw = {r['file']: r for r in csv.DictReader(open(MARKS / 'raw_files.tsv', encoding='utf-8'), delimiter='\t')}
    assert raw[mod.RAW]['md5'] == s['raw_md5'] and raw[mod.RAW]['url'] == mod.URL == s['url']
    # in the release: every staged site, with its columns; no other row carries them
    table = pd.read_csv(RELEASE / f'human_{mark}_sites_v1.5.0.tsv', sep='\t', dtype=str, keep_default_na=False)
    own = table[table['sources'].str.split('; ').map(lambda x: mod.SOURCE in x)]
    assert len(own) == len(staged) and not table.duplicated(['canonical_UniProtAC', 'position', 'state']).any()
    cols = build.SOURCE_COLUMNS[mod.SOURCE]
    merged = own.set_index(['canonical_UniProtAC', 'position'])[cols]
    assert merged.sort_index().equals(staged.set_index(['canonical_UniProtAC', 'position'])[cols].sort_index())
    assert set(table.loc[~table.index.isin(own.index), cols].sum(axis=1)) == {''}
    assert (own['evidence_class'] == 'experimental').all() and own['source_releases'].str.contains(mod.SOURCE_RELEASE, regex=False).all()


def test_akimov_has_no_n_terminal_sites():
    s = json.loads((MARKS / a.OUT_SUMMARY).read_text(encoding='utf-8'))
    assert 'not included' in s['n_terminal_sites']

"""Tests of the GO biological process staging (scripts/63_stage_go_bp.py) and of the staged tables in inputs/annotation/.

Run from the repository root: python -m pytest -q tests/test_go_bp.py
"""
import csv
import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
spec = importlib.util.spec_from_file_location('stage63', ROOT / 'scripts' / '63_stage_go_bp.py')
stage63 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(stage63)

ANNOTATION = ROOT / 'inputs' / 'annotation'
DATE = '2026-11-04'

# root <- A (is_a) <- B (is_a), C part_of A, D is_a B and part_of C, R regulates A (not followed), M a function,
# O obsolete, X has alt_id XA.
OBO = """format-version: 1.2
data-version: releases/2026-07-26

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process

[Term]
id: GO:0000001
name: A
namespace: biological_process
is_a: GO:0008150 ! biological_process

[Term]
id: GO:0000002
name: B
namespace: biological_process
is_a: GO:0000001 ! A

[Term]
id: GO:0000003
name: C
namespace: biological_process
relationship: part_of GO:0000001 ! A

[Term]
id: GO:0000004
name: D
namespace: biological_process
alt_id: GO:0000044
is_a: GO:0000002 ! B
relationship: part_of GO:0000003 ! C

[Term]
id: GO:0000005
name: R
namespace: biological_process
is_a: GO:0008150
relationship: regulates GO:0000001 ! A

[Term]
id: GO:0003674
name: molecular_function
namespace: molecular_function

[Term]
id: GO:0000009
name: O
namespace: biological_process
is_obsolete: true

[Typedef]
id: part_of
name: part of
"""


def gaf_row(acc, qualifier, go, code, aspect='P'):
    return '\t'.join(['UniProtKB', acc, 'G', qualifier, go, 'PMID:1', code, '', aspect, 'name', '', 'protein', 'taxon:9606',
                      '20260101', 'UniProt', '', '']) + '\n'


GAF = ['!gaf-version: 2.2\n', '!date-generated: 2026-11-04 15:05\n',
       gaf_row('P00001', 'involved_in', 'GO:0000044', 'IDA'),          # alt id of D
       gaf_row('P00001', 'involved_in', 'GO:0000004', 'IEA'),
       gaf_row('P00001', 'involved_in', 'GO:0000005', 'IEA'),          # R: regulates A, A not inherited
       gaf_row('P00001', 'enables', 'GO:0003674', 'IDA', aspect='F'),  # other aspect
       gaf_row('P00002', 'NOT|involved_in', 'GO:0000002', 'IMP'),      # NOT
       gaf_row('P00002', 'acts_upstream_of_or_within', 'GO:0000003', 'IEA'),
       gaf_row('P00002', 'involved_in', 'GO:0008150', 'ND'),           # ND dropped
       gaf_row('P00002', 'involved_in', 'GO:0000009', 'TAS'),          # obsolete
       gaf_row('Q99999', 'involved_in', 'GO:0000002', 'IDA'),          # secondary of P00003
       gaf_row('P99999', 'involved_in', 'GO:0000002', 'IDA')]          # not in reference


def staged():
    terms, header = stage63.read_obo(OBO.splitlines(True))
    parents, alt = stage63.bp_graph(terms)
    anc, depth = stage63.ancestors_of(parents), stage63.depths(parents)
    direct, counts, unknown, gaf_header, accs = stage63.read_gaf(GAF, {'P00001', 'P00002', 'P00003'}, {'Q99999': 'P00003'},
                                                                 terms, parents, alt)
    rows = stage63.build_rows(stage63.propagate(direct, anc))
    return terms, header, parents, depth, direct, counts, unknown, gaf_header, rows


def test_toy_propagation_and_filters():
    terms, header, parents, depth, direct, counts, unknown, gaf_header, rows = staged()
    assert header['data-version'] == 'releases/2026-07-26' and gaf_header['date-generated'].startswith('2026-11-04')
    assert 'GO:0003674' not in parents and 'GO:0000009' not in parents
    assert parents['GO:0000004'] == ['GO:0000002', 'GO:0000003'] and parents['GO:0000005'] == ['GO:0008150']
    assert depth['GO:0000004'] == (3, 3) and depth['GO:0000003'] == (2, 2) and depth['GO:0008150'] == (0, 0)
    assert direct == {('P00001', 'GO:0000004'): {'IDA', 'IEA'}, ('P00001', 'GO:0000005'): {'IEA'},
                      ('P00002', 'GO:0000003'): {'IEA'}, ('P00003', 'GO:0000002'): {'IDA'}}
    assert counts['bp_rows_not_qualifier'] == 1 and counts['bp_rows_nd_evidence'] == 1
    assert counts['bp_rows_not_in_reference'] == 1 and counts['bp_rows_alt_id_mapped'] == 1
    assert unknown == {('GO:0000009', 'obsolete'): 1}
    got = {(r[0], r[1]): r[2:] for r in rows}
    assert got[('P00001', 'GO:0000004')] == ['direct', 'IDA;IEA', 'IDA;IEA', 'true']
    assert got[('P00001', 'GO:0000001')] == ['propagated', '', 'IDA;IEA', 'true']     # via D, not via R (regulates)
    assert got[('P00001', 'GO:0008150')] == ['propagated', '', 'IDA;IEA', 'true']
    assert got[('P00001', 'GO:0000005')] == ['direct', 'IEA', 'IEA', 'false']
    assert got[('P00002', 'GO:0000001')] == ['propagated', '', 'IEA', 'false']        # part_of followed
    assert ('P00002', 'GO:0000002') not in got                                        # NOT not used
    assert sorted(t for a, t in got if a == 'P00001') == ['GO:0000001', 'GO:0000002', 'GO:0000003', 'GO:0000004',
                                                          'GO:0000005', 'GO:0008150']
    assert sorted(t for a, t in got if a == 'P00003') == ['GO:0000001', 'GO:0000002', 'GO:0008150']
    table = {r[0]: r for r in stage63.build_terms(rows, terms, parents, depth)}
    assert table['GO:0008150'][4:7] == [3, 2, 0] and table['GO:0000004'][4:8] == [1, 1, 1, 'GO:0000002;GO:0000003']


def read(name):
    with open(ANNOTATION / name, encoding='utf-8', newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


@pytest.mark.skipif(not (ANNOTATION / f'go_bp_human_{DATE}.tsv').exists(), reason='staged table not present')
def test_staged_tables_agree():
    rows = read(f'go_bp_human_{DATE}.tsv')
    terms = {r['go_id']: r for r in read(f'go_bp_terms_{DATE}.tsv')}
    summary = json.loads((ANNOTATION / f'go_bp_human_{DATE}_summary.json').read_text())
    assert list(rows[0]) == stage63.ROW_COLUMNS and list(next(iter(terms.values()))) == stage63.TERM_COLUMNS
    assert len(rows) == summary['rows'] and len(terms) == summary['terms_with_proteins']
    reference = {l.split('|')[1] for l in open(ROOT / 'inputs/reference/uniprot_human_reviewed_2026_03.fasta') if l.startswith('>')}
    proteins = {r['uniprot_accession'] for r in rows}
    assert proteins <= reference and len(proteins) == summary['proteins']
    assert len({(r['uniprot_accession'], r['go_id']) for r in rows}) == len(rows)
    for r in rows:
        assert r['go_id'] in terms
        assert (r['annotation'] == 'direct') == bool(r['direct_evidence_codes'])
        assert set(r['direct_evidence_codes'].split(';')) <= set(r['evidence_codes'].split(';')) or not r['direct_evidence_codes']
        assert r['non_iea'] == ('true' if set(r['evidence_codes'].split(';')) - {'IEA'} else 'false')
    # every protein has the root; counts in the terms table match the rows; parents of a protein's term are its terms too
    assert int(terms['GO:0008150']['n_proteins']) == len(proteins)
    by_term = {}
    for r in rows:
        by_term[r['go_id']] = by_term.get(r['go_id'], 0) + 1
    assert all(int(terms[t]['n_proteins']) == n for t, n in by_term.items())
    have = {(r['uniprot_accession'], r['go_id']) for r in rows}
    for a, t in sorted(have)[::40]:
        assert all((a, p) in have for p in terms[t]['parents'].split(';') if p)
    assert summary['files_md5']['rows'] == stage63.md5_file(ANNOTATION / f'go_bp_human_{DATE}.tsv')

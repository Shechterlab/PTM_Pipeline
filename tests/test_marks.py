"""Tests of the marks tables: scripts/54_stage_uniprot_ptm_sites.py and scripts/55_build_marks_release.py.

Run from the repository root: python -m pytest -q tests
"""
import csv
import importlib.util
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))


def load(name):
    spec = importlib.util.spec_from_file_location(name.replace('.py', '').replace('5', 'm5', 1), ROOT / 'scripts' / name)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


stage = load('54_stage_uniprot_ptm_sites.py')
build = load('55_build_marks_release.py')

REFERENCE = {'P00001': 'MKKRSTYKR', 'P00002': 'MKAK'}


def feature(term, pos, codes=('ECO:0000269',), pmid=None, isoform=None, end=None, modifier='EXACT'):
    loc = {'start': {'value': pos, 'modifier': modifier}, 'end': {'value': end or pos, 'modifier': 'EXACT'}}
    if isoform:
        loc['sequence'] = isoform
    evidences = [{'evidenceCode': c, 'source': 'PubMed', 'id': pmid} if pmid else {'evidenceCode': c} for c in codes]
    return {'type': 'Modified residue', 'description': term, 'location': loc, 'evidences': evidences}


def entry(acc, seq, *features):
    return {'primaryAccession': acc, 'sequence': {'value': seq}, 'features': list(features)}


def run(entries):
    return stage.stage(entries, REFERENCE, '2026_04')


def test_residue_check_drops_and_counts():
    tables, dropped = run([entry('P00001', REFERENCE['P00001'],
                                 feature('N6-methyllysine', 2),            # K: kept
                                 feature('N6-methyllysine', 5),            # S: dropped
                                 feature('Phosphoserine', 5),              # S: kept
                                 feature('Phosphothreonine', 6),           # T: kept
                                 feature('Phosphotyrosine', 6),            # T: dropped
                                 feature('Omega-N-methylarginine', 9),     # R: kept
                                 feature('Omega-N-methylarginine', 40))])  # past the end: dropped
    assert list(tables['methyllysine']['position']) == [2]
    assert list(tables['phosphorylation']['position']) == [5, 6]
    assert list(tables['methylarginine']['position']) == [9]
    assert dropped == {'residue_mismatch': 3}


def test_other_drops_are_counted_separately():
    tables, dropped = run([
        entry('P00001', REFERENCE['P00001'], feature('N6-methyllysine', 2, isoform='P00001-2'),
              feature('N6-methyllysine', 2, end=3), feature('N6-methyllysine', 2, modifier='UNKNOWN')),
        entry('P00002', 'MKAKX', feature('N6-methyllysine', 2)),     # sequence differs from the reference
        entry('Q99999', 'MKK', feature('N6-methyllysine', 2))])      # not in the reference
    assert all(len(t) == 0 for t in tables.values())
    assert dropped == {'on_an_isoform': 1, 'not_single_residue': 2, 'sequence_differs': 1, 'not_in_reference': 1}


def test_state_mapping():
    seq = REFERENCE['P00001']
    tables, _ = run([entry('P00001', seq, feature('N6-methyllysine', 2), feature('N6,N6-dimethyllysine', 2),
                           feature('N6,N6,N6-trimethyllysine', 3), feature('N6-methylated lysine', 8),
                           feature('Asymmetric dimethylarginine', 4), feature('Symmetric dimethylarginine', 9),
                           feature('Dimethylated arginine', 4), feature('N6-acetyllysine', 2),
                           feature('N6-acetyl-N6-methyllysine', 2), feature('(Microbial infection) Phosphoserine', 5),
                           feature('N6-methyllysine; alternate', 3))])
    ml = tables['methyllysine']
    assert list(zip(ml['position'], ml['state'], strict=True)) == [(2, 'me1'), (2, 'me2'), (3, 'me1'), (3, 'me3'), (8, 'me')]
    assert set(ml['mark']) == {'methyllysine'} and set(ml['residue']) == {'K'}
    ma = tables['methylarginine']
    assert sorted(zip(ma['position'], ma['state'], strict=True)) == [(4, 'me2'), (4, 'me2a'), (9, 'me2s')]
    assert len(tables['phosphorylation']) == 0     # a microbial-infection phosphoserine is not the entry's own mark


def test_evidence_class_and_pmids():
    assert stage.evidence_class(['ECO:0000250']) == 'by-similarity'
    assert stage.evidence_class(['ECO:0000250', 'ECO:0000269']) == 'experimental'
    assert stage.evidence_class(['ECO:0007744']) == 'experimental'
    assert stage.evidence_class(['ECO:0000305']) == 'by-similarity'
    assert stage.evidence_class([]) == 'by-similarity'
    tables, _ = run([entry('P00001', REFERENCE['P00001'], feature('N6-methyllysine', 2, ('ECO:0000269',), pmid='12345'),
                           feature('N6-methyllysine', 2, ('ECO:0000250',)), feature('N6-methyllysine', 3, ('ECO:0000250',)))])
    ml = tables['methyllysine']
    assert list(ml['evidence']) == ['ECO:0000250;ECO:0000269', 'ECO:0000250']
    assert list(ml['evidence_class']) == ['experimental', 'by-similarity']
    assert list(ml['pmids']) == ['12345', '']
    assert list(ml['source']) == ['UniProt'] * 2 and list(ml['source_release']) == ['2026_04'] * 2
    assert list(ml.columns) == stage.COLUMNS


def test_staging_is_deterministic():
    entries = [entry('P00001', REFERENCE['P00001'], feature('N6-methyllysine', 3), feature('N6-methyllysine', 2),
                     feature('Phosphoserine', 5)), entry('P00002', REFERENCE['P00002'], feature('N6-methyllysine', 2))]
    a, _ = run(entries)
    b, _ = run(list(reversed(entries)))
    for mark in a:
        assert a[mark].to_csv(sep='\t', index=False) == b[mark].to_csv(sep='\t', index=False)
    assert list(a['methyllysine']['canonical_UniProtAC']) == ['P00001', 'P00001', 'P00002']


def test_staged_tables_and_release_agree_with_manifest():
    manifest = list(csv.DictReader(open(ROOT / 'inputs' / 'MANIFEST.tsv', encoding='utf-8'), delimiter='\t'))
    rows = {r['file']: r for r in manifest if r['file'].startswith('marks/')}
    tables = [f for f in rows if f.endswith('_sites.tsv') or f.endswith('_crosscheck.tsv')]
    assert len(tables) >= 3
    for name in tables + ['marks/raw_files.tsv', 'marks/uniprot_ptm_sites_summary.json']:
        row = rows[name]
        assert row['redistribution'] == 'git'
        assert build.md5_file(ROOT / 'inputs' / name) == row['md5']
    for name in tables:
        row = rows[name]
        assert row['licence'] == 'CC BY 4.0' and 'UniProt' in row['citation'] and row['source'].startswith('https://rest.uniprot.org')
        with open(ROOT / 'inputs' / name, encoding='utf-8') as fh:
            assert sum(1 for _ in fh) - 1 == int(row['rows'])
            fh.seek(0)
            assert next(csv.reader(fh, delimiter='\t')) == stage.COLUMNS


def test_staged_residues_match_reference():
    reference = stage.parse_fasta(ROOT / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03.fasta')
    for path in (ROOT / 'inputs' / 'marks').glob('uniprot_*_methyllysine_sites.tsv'):
        for r in csv.DictReader(open(path, encoding='utf-8'), delimiter='\t'):
            assert r['residue'] == 'K' and reference[r['canonical_UniProtAC']][int(r['position']) - 1] == 'K'


def test_release_files_checksums_and_columns():
    for folder in (ROOT / 'releases' / 'marks').iterdir():
        sums = dict(line.split('  ', 1)[::-1] for line in (folder / 'MD5SUMS').read_text().splitlines())
        for name, md5 in ((n.strip(), m) for n, m in sums.items()):
            assert build.md5_file(folder / name) == md5
        for path in folder.glob('human_*_sites_*.tsv'):
            with open(path, encoding='utf-8') as fh:
                assert next(csv.reader(fh, delimiter='\t')) == build.COLUMNS
        assert 'CC BY 4.0' in (folder / 'README.md').read_text()


def test_release_licence_needs_every_source_open():
    assert {s['license'] for s in build.SOURCES.values()} == {'CC BY 4.0'}


def test_no_closed_source_is_referenced():
    banned = re.compile(r'cplm|phosphosite|iptmnet|dbptm|qptm', re.I)
    files = [ROOT / 'scripts' / '54_stage_uniprot_ptm_sites.py', ROOT / 'scripts' / '55_build_marks_release.py']
    files += [p for p in (ROOT / 'inputs' / 'marks').iterdir() if p.suffix in ('.tsv', '.json') and p.name != 'raw_files.tsv']
    for path in files:
        text = path.read_text(encoding='utf-8')
        if path.suffix == '.py':
            # the only mentions are the statement of what is excluded
            text = '\n'.join(ln for ln in text.splitlines() if not re.search(r'Excluded by design|PhosphoSitePlus, iPTMnet', ln))
        assert not banned.search(text), path
    rows = [r for r in csv.DictReader(open(ROOT / 'inputs' / 'MANIFEST.tsv', encoding='utf-8'), delimiter='\t')
            if r['file'].startswith('marks/') and r['file'].endswith('_sites.tsv')]
    assert all(not banned.search(r['source'] + r['derived_from'] + r['how_made_or_checked']) for r in rows)
    for path in (ROOT / 'releases' / 'marks').glob('*/human_*.tsv'):
        sources = {r['source'] for r in csv.DictReader(open(path, encoding='utf-8'), delimiter='\t')}
        assert sources <= set(build.SOURCES), path


def test_crosscheck_column(tmp_path):
    import pandas as pd
    (tmp_path / 'human_methylarginine_sites_v1.4.3.tsv').write_text('uniprot_accession\tposition\nA\t1\n')
    (tmp_path / 'human_methylarginine_lower_confidence_sites_v1.4.3.tsv').write_text('uniprot_accession\tposition\nA\t2\n')
    frame = pd.DataFrame({'canonical_UniProtAC': ['A', 'A', 'A'], 'position': ['1', '2', '3']})
    assert list(build.crosscheck(frame, tmp_path)['in_methylarginine_release']) == ['site', 'lower_confidence', 'no']


@pytest.mark.parametrize('name', ['54_stage_uniprot_ptm_sites.py', '55_build_marks_release.py'])
def test_scripts_use_american_spelling_and_no_banned_words(name):
    text = (ROOT / 'scripts' / name).read_text(encoding='utf-8')
    text = text.replace("'licence'", '')     # the MANIFEST column keeps its name
    assert not re.search(r'dtIDR|domain-tethered|domain flanks|licen[cs]e\b(?<!license)|colour|summaris|generalis', text, re.I)

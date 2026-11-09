"""Tests of the context-mark tables: scripts/54, 55 (staging) and 56 (releases/v1.5.0).

Run from the repository root: python -m pytest -q tests
"""
import csv
import importlib.util
import json
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
iptm = load('55_stage_iptmnet_ptm_sites.py')
build = load('56_build_context_marks.py')

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
    assert len(tables) == 6
    for name in tables + ['marks/raw_files.tsv', 'marks/uniprot_ptm_sites_summary.json', 'marks/iptmnet_ptm_sites_summary.json']:
        assert rows[name]['redistribution'] == 'git'
        assert build.md5_file(ROOT / 'inputs' / name) == rows[name]['md5']
    for name in tables:
        row = rows[name]
        with open(ROOT / 'inputs' / name, encoding='utf-8') as fh:
            assert sum(1 for _ in fh) - 1 == int(row['rows'])
            fh.seek(0)
            assert next(csv.reader(fh, delimiter='\t'))[:10] == stage.COLUMNS
        if '/uniprot_' in name:
            assert row['licence'] == 'CC BY 4.0' and row['source'].startswith('https://rest.uniprot.org')
        else:
            assert row['licence'].startswith('CC BY-NC-SA 4.0') and row['source'].startswith('https://research.bioinformatics.udel.edu/iptmnet_data/')
            assert 'iPTMnet' in row['citation']


def test_staged_residues_match_reference():
    reference = stage.parse_fasta(ROOT / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03.fasta')
    for path in (ROOT / 'inputs' / 'marks').glob('*_sites.tsv'):
        for r in csv.DictReader(open(path, encoding='utf-8'), delimiter='\t'):
            assert reference[r['canonical_UniProtAC']][int(r['position']) - 1] == r['residue']
            assert r['residue'] in {'methyllysine': 'K', 'phosphorylation': 'STY'}[r['mark']]


def test_iptmnet_staging():
    import pandas as pd
    fields = iptm.FIELDS
    def rec(ptm_type, source, acc, site, org='Homo sapiens (Human)', pmid=''):
        return dict(zip(fields, [ptm_type, source, acc, 'G', org, site, '', '', '', pmid], strict=True))
    records = pd.DataFrame([
        rec('METHYLATION', 'unip', 'P00001', 'K2'),                   # K: methyllysine, state me, class unspecified
        rec('METHYLATION', 'hprd', 'P00001-1', 'K2', pmid='11,5'),    # the -1 isoform is canonical; curated: experimental
        rec('METHYLATION', 'iedb', 'P00001', 'C1'),                   # not K or R
        rec('METHYLATION', 'unip', 'P00001', 'R9'),                   # methylarginine cross-check
        rec('PHOSPHORYLATION', 'rlim', 'P00001', 'S5', pmid='7'),     # text-mined
        rec('PHOSPHORYLATION', 'hprd', 'P00001', 'Y5'),               # residue mismatch (S)
        rec('PHOSPHORYLATION', 'hprd', 'P00001-2', 'S5'),             # isoform
        rec('PHOSPHORYLATION', 'hprd', 'Q99999', 'S5'),               # not in the reference
        rec('PHOSPHORYLATION', 'hprd', 'P00001', 'S5', org='Mus musculus (Mouse)'),
        rec('ACETYLATION', 'hprd', 'P00001', 'K2')])
    tables, dropped = iptm.stage(records, REFERENCE, '6.2')
    ml = tables['methyllysine']
    assert list(ml['state']) == ['me'] and list(ml['upstream_sources']) == ['hprd;unip']
    assert list(ml['evidence_class']) == ['experimental'] and list(ml['pmids']) == ['5;11']
    assert list(tables['methylarginine']['position']) == [9]
    assert list(tables['phosphorylation']['evidence_class']) == ['text-mined']
    assert dropped == {'other_residue': 1, 'residue_mismatch': 1, 'on_an_isoform': 1, 'not_in_reference': 1, 'not_human': 1, 'other_ptm_type': 1}
    assert iptm.evidence_class(['unip']) == 'unspecified' and iptm.evidence_class(['rlim', 'unip']) == 'text-mined'
    assert iptm.evidence_class(['unip', 'sign']) == 'experimental'


def test_merge_of_sources_is_deterministic_and_keeps_attribution():
    import pandas as pd
    u = pd.DataFrame([{**dict.fromkeys(stage.COLUMNS, ''), 'canonical_UniProtAC': 'P00001', 'position': '2', 'residue': 'K',
                       'mark': 'methyllysine', 'state': 'me1', 'evidence': 'ECO:0000250', 'evidence_class': 'by-similarity',
                       'source': 'UniProt', 'source_release': '2026_04', 'upstream_sources': ''}])
    i = pd.DataFrame([{**dict.fromkeys(stage.COLUMNS, ''), 'canonical_UniProtAC': 'P00001', 'position': '2', 'residue': 'K',
                       'mark': 'methyllysine', 'state': 'me1', 'evidence_class': 'experimental', 'pmids': '5',
                       'source': 'iPTMnet', 'source_release': '6.2', 'upstream_sources': 'hprd'},
                      {**dict.fromkeys(stage.COLUMNS, ''), 'canonical_UniProtAC': 'P00001', 'position': '3', 'residue': 'K',
                       'mark': 'methyllysine', 'state': 'me', 'evidence_class': 'unspecified', 'source': 'iPTMnet',
                       'source_release': '6.2', 'upstream_sources': 'unip'}])
    a, b = build.merge([u, i]), build.merge([i, u])
    assert a.to_csv(sep='\t', index=False) == b.to_csv(sep='\t', index=False)
    first = a.iloc[0]
    assert first['sources'] == 'UniProt; iPTMnet' and first['evidence_class'] == 'experimental'
    assert first['evidence'] == 'ECO:0000250' and first['iptmnet_sources'] == 'hprd' and first['source_releases'] == '2026_04; 6.2'
    assert a.iloc[1]['sources'] == 'iPTMnet'
    assert list(a.columns) == build.OUT_COLUMNS
    ov = build.overlap(a)
    assert ov['sites'] == {'both': 1, 'iPTMnet only': 1}


def test_release_v1_5_0_is_consistent():
    folder = ROOT / 'releases' / 'v1.5.0'
    for line in (folder / 'MD5SUMS').read_text().splitlines():
        md5, name = line.split('  ', 1)
        assert build.md5_file(folder / name) == md5, name
    release = json.loads((folder / 'release.json').read_text())
    previous = ROOT / 'releases' / 'v1.4.3'
    assert len(release['carried_over']) > 90
    for name, md5 in release['carried_over'].items():       # the methylarginine data: byte for byte as in v1.4.3
        assert build.md5_file(previous / name) == md5 == build.md5_file(folder / name), name
    for path in folder.glob('human_*_v1.5.0.tsv'):
        with open(path, encoding='utf-8') as fh:
            header = next(csv.reader(fh, delimiter='\t'))
        assert header[:11] == build.OUT_COLUMNS
    text = (folder / 'README.md').read_text() + (folder / 'CHANGES.md').read_text()
    assert 'unchanged from v1.4.3' in text and 'CC BY-NC-SA 4.0' in text and 'CC BY 4.0' in text
    assert release['license'].startswith('CC BY-NC-SA 4.0')
    assert not (ROOT / 'releases' / 'marks').exists()
    assert 'CC BY 4.0' not in (ROOT / 'LICENSE.md').read_text().split('## Code')[0].split('Data, figures')[-1]


def test_release_license_covers_every_source():
    # the release is CC BY-NC-SA 4.0: every source must be CC BY, CC BY-NC-SA or as open
    assert {s['license'] for s in build.SOURCES.values()} <= {'CC BY 4.0', 'CC BY-NC-SA 4.0'}
    assert set(build.SOURCES) == {'UniProt', 'iPTMnet'}


def test_no_closed_source_is_referenced():
    banned = re.compile(r'cplm|phosphosite|dbptm|qptm|\bpsp\b', re.I)
    files = [ROOT / 'scripts' / n for n in ('54_stage_uniprot_ptm_sites.py', '55_stage_iptmnet_ptm_sites.py', '56_build_context_marks.py')]
    files += [p for p in (ROOT / 'inputs' / 'marks').iterdir() if p.suffix in ('.tsv', '.json') and p.name != 'raw_files.tsv']
    files += list((ROOT / 'releases' / 'v1.5.0').glob('human_*_v1.5.0.tsv'))
    for path in files:
        text = path.read_text(encoding='utf-8')
        if path.suffix == '.py':    # only the statement of what is not used
            text = '\n'.join(ln for ln in text.splitlines() if not re.search(r'Not used|Berryhill2023 is added', ln))
        assert not banned.search(text), path
    rows = [r for r in csv.DictReader(open(ROOT / 'inputs' / 'MANIFEST.tsv', encoding='utf-8'), delimiter='\t')
            if r['file'].startswith('marks/')]
    assert all(not banned.search(r['source'] + r['derived_from'] + r['how_made_or_checked']) for r in rows)
    for path in (ROOT / 'releases' / 'v1.5.0').glob('human_*_v1.5.0.tsv'):
        for r in csv.DictReader(open(path, encoding='utf-8'), delimiter='\t'):
            assert set(r['sources'].split('; ')) <= set(build.SOURCES)
            assert not (set(r['iptmnet_sources'].split(';')) - set(iptm_upstream()))


def iptm_upstream():
    return ['hprd', 'unip', 'rlim', 'rlim+', 'sign', 'npro', 'iedb', 'pro', 'intact', '']


def test_crosscheck_column(tmp_path):
    import pandas as pd
    (tmp_path / 'human_methylarginine_sites_v1.4.3.tsv').write_text('uniprot_accession\tposition\nA\t1\n')
    (tmp_path / 'human_methylarginine_lower_confidence_sites_v1.4.3.tsv').write_text('uniprot_accession\tposition\nA\t2\n')
    frame = pd.DataFrame({'canonical_UniProtAC': ['A', 'A', 'A'], 'position': ['1', '2', '3']})
    assert list(build.crosscheck(frame, tmp_path)['in_methylarginine_release']) == ['site', 'lower_confidence', 'no']


@pytest.mark.parametrize('name', ['54_stage_uniprot_ptm_sites.py', '55_stage_iptmnet_ptm_sites.py', '56_build_context_marks.py'])
def test_scripts_use_american_spelling_and_no_banned_words(name):
    text = (ROOT / 'scripts' / name).read_text(encoding='utf-8')
    text = text.replace("'licence'", '')     # the MANIFEST column keeps its name
    assert not re.search(r'dtIDR|domain-tethered|domain flanks|licen[cs]e\b(?<!license)|colour|summaris|generalis', text, re.I)

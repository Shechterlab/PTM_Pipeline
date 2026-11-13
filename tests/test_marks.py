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
cplm = load('39_stage_cplm4_lys_methylation.py')
score = load('59_stage_iptmnet_scores.py')

def context_files(folder):
    return [folder / f'human_{m}_v1.5.0.tsv' for m in ('methyllysine_sites', 'phosphorylation_sites', 'lysine_acylation_sites', 'acetylation_sites', 'ubiquitination_sites',
                                                                       'sumoylation_sites', 'methylarginine_crosscheck')]


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
    assert len(tables) == 19 + len(cplm.ACYL_TYPES)
    for name in tables + ['marks/raw_files.tsv', 'marks/uniprot_ptm_sites_summary.json', 'marks/iptmnet_ptm_sites_summary.json',
                          'marks/cplm_ptm_sites_summary.json', 'marks/cplm_ubiquitination_summary.json', 'marks/cplm_sumoylation_summary.json',
                          *[f'marks/{c["summary"]}' for c in cplm.MODIFICATIONS.values() if c['mark'] == 'lysine_acylation']]:
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
        elif '/berryhill2023_' in name:
            assert 'CC BY 4.0' in row['licence'] and 'no separate notice' in row['licence'] and 'owner' not in row['licence']
            assert 'Berryhill' in row['citation'] and 'Sci Rep 13:377' in row['citation'] and '10.1038/s41598-022-27175-x' in row['citation']
        elif '/ochoa2020_' in name:
            assert 'not open access' in row['licence'] and 'exclusive licence to Springer Nature' in row['licence'] and 'owner' not in row['licence']
            assert 'Ochoa' in row['citation'] and 'Nat Biotechnol 38:365' in row['citation'] and '10.1038/s41587-019-0344-3' in row['citation']
            assert row['redistribution'] == 'git'
        elif '/hendriks2017_' in name or '/akimov2018_' in name:
            assert 'not open access' in row['licence'] and 'no separate notice' in row['licence'] and 'owner' not in row['licence']
            assert ('Nat Struct Mol Biol 24:325' in row['citation'] and '10.1038/nsmb.3366' in row['citation']) if 'hendriks' in name else (
                'Nat Struct Mol Biol 25:631' in row['citation'] and '10.1038/s41594-018-0084-y' in row['citation'])
            assert row['redistribution'] == 'git'
        elif '/cplm_' in name:
            assert row['licence'].startswith('none stated') and 'All Rights Reserved' in row['licence'] and 'pending' not in row['licence']
            assert row['source'].startswith('https://cplm.biocuckoo.cn/') and 'Zhang W' in row['citation'] and 'gkab849' in row['citation']
        else:
            assert row['licence'].startswith('CC BY-NC-SA 4.0') and row['source'].startswith('https://research.bioinformatics.udel.edu/iptmnet_data/')
            assert 'iPTMnet' in row['citation']


def test_iptmnet_score_is_staged_with_its_readme_wording_and_checksums():
    row = {r['file']: r for r in csv.DictReader(open(ROOT / 'inputs' / 'MANIFEST.tsv', encoding='utf-8'), delimiter='\t')}['marks/iptmnet_6_2_score.tsv']
    assert row['redistribution'] == 'git' and row['licence'].startswith('CC BY-NC-SA 4.0') and row['source'].endswith('/iptmnet_data/files/current/score.txt')
    assert build.md5_file(ROOT / 'inputs' / 'marks' / 'iptmnet_6_2_score.tsv') == row['md5']
    summary = json.loads((ROOT / 'inputs' / 'marks' / 'iptmnet_score_summary.json').read_text(encoding='utf-8'))
    assert 'The score is calculated based on all the evidence sources as shown on the website' in summary['readme_says']
    assert summary['raw_md5']['score.txt'] == '34d5a1be8f1150694f445d58245f2f30' and summary['score_range_in_download'] == [0, 4]
    raw = {r['file']: r for r in csv.DictReader(open(ROOT / 'inputs' / 'marks' / 'raw_files.tsv', encoding='utf-8'), delimiter='\t')}
    assert raw['iptmnet_score.txt']['md5'] == '34d5a1be8f1150694f445d58245f2f30' and len(raw['iptmnet_score.txt']['sha256']) == 64
    assert raw['iptmnet_ptm.txt']['md5'] == '058d6619f5c6cb650107695a257e6fbf'
    table = list(csv.DictReader(open(ROOT / 'inputs' / 'marks' / 'iptmnet_6_2_score.tsv', encoding='utf-8'), delimiter='\t'))
    assert len(table) == int(row['rows']) == summary['rows'] and set(r['iptmnet_score'] for r in table) <= {'0', '1', '2', '3', '4'}
    assert len({(r['canonical_UniProtAC'], r['position'], r['mark']) for r in table}) == len(table)
    reference = stage.parse_fasta(ROOT / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03.fasta')
    for r in table:
        assert reference[r['canonical_UniProtAC']][int(r['position']) - 1] == r['residue']


def test_iptmnet_score_staging_toy():
    import pandas as pd
    def rec(acc, site, enzyme, ptm, value):
        return dict(zip(score.FIELDS, [acc, site, enzyme, ptm, value], strict=True))
    records = pd.DataFrame([rec('P00001', 'K2', '', 'Methylation', '2'), rec('P00001-1', 'K2', '', 'Methylation', '0'),   # plain accession wins
                            rec('P00001-1', 'K3', '', 'Acetylation', '1'),                                              # only the -1 row: used
                            rec('P00001', 'S5', 'P00002', 'Phosphorylation', '3'),                                      # enzyme row: not used
                            rec('P00001', 'S5', '', 'Phosphorylation', '1'), rec('P00001', 'K2', '', 'S-Nitrosylation', '4'),
                            rec('P00001-2', 'K2', '', 'Methylation', '1'), rec('P99999', 'K2', '', 'Methylation', '1'),
                            rec('P00001', 'R2', '', 'Methylation', '1'), rec('P00001', 'S2', '', 'Phosphorylation', '1'),
                            rec('P00001', 'S', '', 'Phosphorylation', '1')])
    table, dropped = score.stage(records, {'P00001': 'MKKRSTYKR'}, '6.2')
    assert dropped['_conflicts'] == 1
    got = {(r.mark, r.position): r.iptmnet_score for r in table.itertuples()}
    assert got == {('methyllysine', 2): '2', ('lysine_acylation', 3): '1', ('phosphorylation', 5): '1'}
    assert list(table.columns) == score.COLUMNS and set(table['source']) == {'iPTMnet'}


def test_staged_residues_match_reference():
    reference = stage.parse_fasta(ROOT / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03.fasta')
    for path in (ROOT / 'inputs' / 'marks').glob('*_sites.tsv'):
        for r in csv.DictReader(open(path, encoding='utf-8'), delimiter='\t'):
            assert reference[r['canonical_UniProtAC']][int(r['position']) - 1] == r['residue']
            assert r['residue'] in {'methyllysine': 'K', 'phosphorylation': 'STY', 'lysine_acylation': 'K', 'ubiquitination': 'K', 'sumoylation': 'K'}[r['mark']]


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
        rec('N-GLYCOSYLATION', 'hprd', 'P00001', 'K2')])
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
    assert ov['sites'] == {'UniProt + iPTMnet': 1, 'iPTMnet only': 1}


def test_release_v1_5_0_is_consistent():
    folder = ROOT / 'releases' / 'v1.5.0'
    for line in (folder / 'MD5SUMS').read_text().splitlines():
        md5, name = line.split('  ', 1)
        assert build.md5_file(folder / name) == md5, name
    release = json.loads((folder / 'release.json').read_text())
    previous = ROOT / 'releases' / 'v1.4.3'
    assert len(release['carried_over']) > 90
    edited = set(release['carried_over_text_edited'])
    for name, info in release['carried_over'].items():       # the methylarginine data: as in v1.4.3, under this release's names
        assert name == info['from'].replace('_v1.4.3', '_v1.5.0')
        assert build.md5_file(previous / info['from']) == info['from_md5'] and build.md5_file(folder / name) == info['md5']
        assert (info['md5'] != info['from_md5']) == (name in edited), name
    assert len(edited) == 6 and 'human_methylarginine_sites_v1.5.0.xlsx' in edited      # strings that name the version
    for path in folder.rglob('*'):                            # every file carries this release's version, none another's
        if path.is_file() and re.search(r'_v\d+\.\d+\.\d+', path.name):
            assert path.name.endswith('_v1.5.0.tsv') or path.name.endswith('_v1.5.0.xlsx') or path.name == 'site_changes_v1.4.3_to_v1.5.0.tsv', path
    assert not list(folder.rglob('*_v1.4.3.*'))
    import pandas as pd
    changes = pd.read_csv(folder / 'supporting' / 'site_changes_v1.4.3_to_v1.5.0.tsv', sep='\t', dtype=str, keep_default_na=False)
    sites = pd.read_csv(folder / 'human_methylarginine_sites_v1.5.0.tsv', sep='\t', dtype=str, keep_default_na=False, usecols=['site_id'])
    assert set(changes['change']) == {'kept'} and set(changes['records_change']) == {'same'} and list(changes.columns) == build.CHANGE_COLUMNS
    assert set(changes['site_id']) == set(sites['site_id']) and len(changes) == len(sites)
    for path in context_files(folder):
        with open(path, encoding='utf-8') as fh:
            header = next(csv.reader(fh, delimiter='\t'))
        expected = [c if c != 'state' else 'acyl_type' for c in build.OUT_COLUMNS] if 'acylation' in path.name or 'acetylation' in path.name else build.OUT_COLUMNS
        assert header[:11] == expected
        if 'crosscheck' in path.name:
            assert header[11:] == ['in_methylarginine_release']
        else:       # the columns of the sources that contribute, then iPTMnet's score
            own = {'acylation': ['cplm_evidence'], 'acetylation': ['cplm_evidence'], 'ubiquitination': ['cplm_evidence', 'akimov_localization_probability', 'akimov_cell_lines'], 'sumoylation': ['cplm_evidence', 'hendriks_localization_probability', 'hendriks_cell_lines'], 'methyllysine': ['cplm_evidence', 'berryhill_site_probability', 'berryhill_cell_lines', 'berryhill_enrichment'],
                   'phosphorylation': ['functional_score', 'best_localization_probability', 'best_pep', 'biological_samples', 'spectral_counts']}
            assert header[11:] == next((v for k, v in own.items() if k in path.name), []) + ['iptmnet_score'], path.name
    text = (folder / 'README.md').read_text() + (folder / 'CHANGES.md').read_text()
    assert 'unchanged from v1.4.3' in text and 'CC BY-NC-SA 4.0' in text and 'CC BY 4.0' in text
    assert release['license'].startswith('CC BY-NC-SA 4.0')
    assert not (ROOT / 'releases' / 'marks').exists()
    assert 'CC BY 4.0' not in (ROOT / 'LICENSE.md').read_text().split('## Code')[0].split('Data, figures')[-1]


def test_release_license_covers_every_source():
    # the release is CC BY-NC-SA 4.0: every source must be CC BY, CC BY-NC-SA or as open
    # (CPLM, Berryhill 2023, Ochoa 2020, Hendriks 2017 and Akimov 2018 are the exceptions, each stating its own terms in its license entry, with no word on how the decision was taken)
    exceptions = {'CPLM', 'Berryhill 2023', 'Ochoa 2020', 'Hendriks 2017', 'Akimov 2018'}
    assert {s['license'] for k, s in build.SOURCES.items() if k not in exceptions} <= {'CC BY 4.0', 'CC BY-NC-SA 4.0'}
    assert all('owner' not in build.SOURCES[k]['license'] and 'pending' not in build.SOURCES[k]['license'] for k in exceptions)
    assert set(build.SOURCES) == {'UniProt', 'iPTMnet', 'CPLM', 'Berryhill 2023', 'Ochoa 2020', 'Hendriks 2017', 'Akimov 2018'} == set(build.SOURCE_ORDER)
    assert all('not open access' in build.SOURCES[s]['license'] and 'supplements carry no notice' in build.SOURCES[s]['license'] for s in ('Hendriks 2017', 'Akimov 2018'))
    assert 'CC BY 4.0' in build.SOURCES['Berryhill 2023']['license'] and 'supplements carry no separate notice' in build.SOURCES['Berryhill 2023']['license']
    assert 'not open access' in build.SOURCES['Ochoa 2020']['license'] and 'supplements carry no license' in build.SOURCES['Ochoa 2020']['license']
    # CPLM is the one source without a license: no license is stated
    assert build.SOURCES['CPLM']['license'] == 'none stated' and 'pending' not in build.SOURCES['CPLM']['citation']


def test_no_closed_source_is_referenced():
    banned = re.compile(r'phosphosite|dbptm|qptm|\bpsp\b', re.I)
    files = [ROOT / 'scripts' / n for n in ('39_stage_cplm4_lys_methylation.py', '54_stage_uniprot_ptm_sites.py', '55_stage_iptmnet_ptm_sites.py', '56_build_context_marks.py')]
    files += [p for p in (ROOT / 'inputs' / 'marks').iterdir() if p.suffix in ('.tsv', '.json') and p.name != 'raw_files.tsv'
              and not p.name.startswith('berryhill2023_')]   # peptide letters can spell a word
    files += list(context_files(ROOT / 'releases' / 'v1.5.0'))
    for path in files:
        text = path.read_text(encoding='utf-8')
        if path.suffix == '.py':    # only the statement of what is not used; the script's own name is not the closed database
            text = text.replace('58_stage_ochoa_phosphosites.py', '')
            text = '\n'.join(ln for ln in text.splitlines() if not re.search(r'Not used|Berryhill2023 is added', ln))
        assert not banned.search(text), path
    rows = [r for r in csv.DictReader(open(ROOT / 'inputs' / 'MANIFEST.tsv', encoding='utf-8'), delimiter='\t')
            if r['file'].startswith('marks/')]
    # the staging script's own file name (58_stage_ochoa_phosphosites.py) contains "phosphosite"; it is not the closed database
    assert all(not banned.search((r['source'] + r['derived_from'] + r['how_made_or_checked']).replace('58_stage_ochoa_phosphosites.py', ''))
               for r in rows)
    for path in context_files(ROOT / 'releases' / 'v1.5.0'):
        for r in csv.DictReader(open(path, encoding='utf-8'), delimiter='\t'):
            assert set(r['sources'].split('; ')) <= set(build.SOURCES)
            names = r['sources'].split('; ')
            assert ('CPLM' in names) == bool(r.get('cplm_evidence', ''))
            assert ('Ochoa 2020' in names) == bool(r.get('functional_score', ''))
            assert ('Berryhill 2023' in names) == bool(r.get('berryhill_cell_lines', ''))
            assert ('Hendriks 2017' in names) == bool(r.get('hendriks_localization_probability', ''))
            assert ('Akimov 2018' in names) == bool(r.get('akimov_localization_probability', ''))
            assert 'iPTMnet' in names or not r.get('iptmnet_score', '')
            assert not (set(r['iptmnet_sources'].split(';')) - set(iptm_upstream()))


def iptm_upstream():
    return ['hprd', 'unip', 'rlim', 'rlim+', 'sign', 'npro', 'iedb', 'pro', 'intact', '']


def test_crosscheck_column(tmp_path):
    import pandas as pd
    (tmp_path / 'human_methylarginine_sites_v1.4.3.tsv').write_text('uniprot_accession\tposition\nA\t1\n')
    (tmp_path / 'human_methylarginine_lower_confidence_sites_v1.4.3.tsv').write_text('uniprot_accession\tposition\nA\t2\n')
    frame = pd.DataFrame({'canonical_UniProtAC': ['A', 'A', 'A'], 'position': ['1', '2', '3']})
    assert list(build.crosscheck(frame, tmp_path)['in_methylarginine_release']) == ['site', 'lower_confidence', 'no']


@pytest.mark.parametrize('name', ['39_stage_cplm4_lys_methylation.py', '54_stage_uniprot_ptm_sites.py', '55_stage_iptmnet_ptm_sites.py', '56_build_context_marks.py',
                                  '57_stage_berryhill_lys_methylation.py', '58_stage_ochoa_phosphosites.py', '59_stage_iptmnet_scores.py',
                                  '60_stage_hendriks2017_sumo_sites.py', '61_stage_akimov2018_ubiquitination_sites.py', 'window_sites.py'])
def test_scripts_use_american_spelling_and_no_banned_words(name):
    text = (ROOT / 'scripts' / name).read_text(encoding='utf-8')
    text = text.replace("'licence'", '').replace('exclusive licence', '')     # the MANIFEST column keeps its name; the publisher's words
    assert not re.search(r'dtIDR|domain-tethered|domain flanks|licen[cs]e\b(?<!license)|colour|summaris|generalis', text, re.I)


def test_site_changes_generator_reproduces_the_pipelines_v1_4_3_table():
    import pandas as pd
    prev = ROOT / 'releases' / 'v1.4.3'
    old = (prev / 'supporting' / 'site_changes_v1.4.2_to_v1.4.3.tsv').read_text(encoding='utf-8')
    sites = pd.read_csv(prev / 'human_methylarginine_sites_v1.4.3.tsv', sep='\t', dtype=str, keep_default_na=False,
                        usecols=['site_id', 'gene', 'uniprot_accession', 'position'])
    prior = pd.read_csv(prev / 'supporting' / 'site_changes_v1.4.2_to_v1.4.3.tsv', sep='\t', dtype=str, keep_default_na=False, usecols=['site_id', 'gene'])
    assert build.site_changes(sites, prior).to_csv(sep='\t', index=False, lineterminator='\n') == old


def test_cplm_staging_pubmed_ids_and_evidence_classes():
    assert cplm.pubmed_ids('31072472;57987737;486334;2363903;31072472;abc;34300001;34300000') == ['2363903', '31072472', '34300000']
    import pandas as pd
    good = pd.DataFrame({'canonical_UniProtAC': ['P00002', 'P00001', 'P00001'], 'position': [3, 2, 4], 'evidence': ['Dat.', 'Exp.', 'Exp./Dat.'],
                         'references': ['57987737', '31072472;1', '']})
    table = cplm.stage_marks(good, 'CPLM 4.0 (2021-06-13)')
    assert list(table.columns) == cplm.MARKS_COLUMNS and list(table.columns[:10]) == stage.COLUMNS[:10]
    assert list(table['canonical_UniProtAC']) == ['P00001', 'P00001', 'P00002'] and set(table['state']) == {'me'}
    assert list(table['evidence_class']) == ['experimental', 'experimental', 'unspecified']
    assert list(table['pmids']) == ['31072472', '', ''] and list(table['cplm_evidence']) == ['Exp.', 'Exp./Dat.', 'Dat.']
    with pytest.raises(SystemExit):
        cplm.stage_marks(good.assign(evidence=['Exp.', 'Exp.', 'Other']), 'x')


def test_cplm_rows_in_the_release_and_the_merge():
    import pandas as pd
    folder = ROOT / 'releases' / 'v1.5.0'
    table = pd.read_csv(folder / 'human_methyllysine_sites_v1.5.0.tsv', sep='\t', dtype=str, keep_default_na=False)
    cp = table[table['sources'].str.contains('CPLM')]
    staged = pd.read_csv(ROOT / 'inputs' / 'marks' / 'cplm_4_0_methyllysine_sites.tsv', sep='\t', dtype=str, keep_default_na=False)
    assert len(staged) == len(cp) == 5356 and cp['canonical_UniProtAC'].nunique() == 2939
    assert not table.duplicated(['canonical_UniProtAC', 'position', 'state']).any()
    assert set(cp['cplm_evidence']) == {'Exp.', 'Exp./Dat.', 'Dat.'} and set(table.loc[~table['sources'].str.contains('CPLM'), 'cplm_evidence']) == {''}
    assert set(cp['state']) <= {'me'} and 'CPLM 4.0 (2021-06-13)' in set(cp['source_releases'].str.split('; ').str[-1])
    reference = stage.parse_fasta(ROOT / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03.fasta')
    for r in cp.itertuples():
        assert r.residue == 'K' and reference[r.canonical_UniProtAC][int(r.position) - 1] == 'K'
    # every CPLM Exp. or Exp./Dat. row is experimental, whatever else reports it
    assert set(cp.loc[cp['cplm_evidence'] != 'Dat.', 'evidence_class']) == {'experimental'}
    # the rows of the other sources are as they were: UniProt and iPTMnet rows without CPLM keep their class
    others = table[~table['sources'].str.contains('CPLM')]
    assert len(others) + len(cp[cp['sources'] != 'CPLM']) >= 756
    text = (folder / 'README.md').read_text() + (folder / 'CHANGES.md').read_text()
    assert 'CPLM states no license' in text and 'All Rights Reserved' in text and 'Zhang et al., Nucleic Acids Res 50:D451' in text


def test_merge_adds_cplm_to_sources_without_changing_other_rows():
    import pandas as pd
    base = dict.fromkeys(stage.COLUMNS, '')
    u = pd.DataFrame([{**base, 'canonical_UniProtAC': 'P00001', 'position': '2', 'residue': 'K', 'mark': 'methyllysine', 'state': 'me',
                       'evidence_class': 'by-similarity', 'evidence': 'ECO:0000250', 'source': 'UniProt', 'source_release': '2026_04'}])
    c = pd.DataFrame([{**base, 'canonical_UniProtAC': 'P00001', 'position': '2', 'residue': 'K', 'mark': 'methyllysine', 'state': 'me',
                       'evidence_class': 'experimental', 'pmids': '9;3', 'source': 'CPLM', 'source_release': 'CPLM 4.0 (2021-06-13)',
                       'cplm_evidence': 'Exp.'},
                      {**base, 'canonical_UniProtAC': 'P00002', 'position': '3', 'residue': 'K', 'mark': 'methyllysine', 'state': 'me',
                       'evidence_class': 'unspecified', 'source': 'CPLM', 'source_release': 'CPLM 4.0 (2021-06-13)', 'cplm_evidence': 'Dat.'}])
    out = build.merge([u.assign(upstream_sources='', cplm_evidence=''), c.assign(upstream_sources='')])
    assert list(out.columns) == build.OUT_COLUMNS + ['cplm_evidence']
    first = out.iloc[0]
    assert first['sources'] == 'UniProt; CPLM' and first['evidence_class'] == 'experimental' and first['evidence'] == 'ECO:0000250'
    assert first['pmids'] == '3;9' and first['source_releases'] == '2026_04; CPLM 4.0 (2021-06-13)' and first['cplm_evidence'] == 'Exp.'
    assert out.iloc[1]['sources'] == 'CPLM' and out.iloc[1]['evidence_class'] == 'unspecified'
    assert build.overlap(out)['sites'] == {'UniProt + CPLM': 1, 'CPLM only': 1}
    assert 'cplm_evidence' not in build.merge([u.assign(upstream_sources='', cplm_evidence='')]).columns


def test_acyl_terms_come_from_the_vocabulary():
    assert stage.acyl_type('N6-acetyllysine') == 'acetyl'
    assert stage.acyl_type('N6-(2-hydroxyisobutyryl)lysine') == '2-hydroxyisobutyryl'
    assert stage.acyl_type('N6-(beta-hydroxybutyryl)lysine') == 'beta-hydroxybutyryl'
    for term in ('N6-methyllysine', 'N6,N6-dimethyllysine', 'N6-acetyl-N6-methyllysine', 'N6-(ADP-ribosyl)lysine', 'N6-lipoyllysine',
                 'N6-biotinyllysine', 'N6-(pyridoxal phosphate)lysine', 'N6-(retinylidene)lysine', 'N6-carboxylysine',
                 '(Microbial infection) N6-acetyllysine', '5-hydroxylysine'):
        assert stage.acyl_type(term) is None, term
    tables, dropped = run([entry('P00001', REFERENCE['P00001'], feature('N6-acetyllysine', 2), feature('N6-succinyllysine; alternate', 2),
                                 feature('N6-acetyllysine', 5), feature('N6-lipoyllysine', 3))])
    acyl = tables['lysine_acylation']
    assert list(zip(acyl['position'], acyl['state'], strict=True)) == [(2, 'acetyl'), (2, 'succinyl')]
    assert dropped == {'residue_mismatch': 1} and set(acyl['mark']) == {'lysine_acylation'}


def test_iptmnet_acetylation_and_acylation_table_of_the_release():
    import pandas as pd
    records = pd.DataFrame([dict(zip(iptm.FIELDS, ['ACETYLATION', 'hprd', 'P00001', 'G', 'Homo sapiens (Human)', 'K2', '', '', '', '9'], strict=True)),
                            dict(zip(iptm.FIELDS, ['ACETYLATION', 'hprd', 'P00001', 'G', 'Homo sapiens (Human)', 'S5', '', '', '', ''], strict=True))])
    tables, dropped = iptm.stage(records, REFERENCE, '6.2')
    assert list(tables['lysine_acylation']['state']) == ['acetyl'] and dropped == {'other_residue': 1}
    folder = ROOT / 'releases' / 'v1.5.0'
    acyl = list(csv.DictReader(open(folder / 'human_lysine_acylation_sites_v1.5.0.tsv', encoding='utf-8'), delimiter='\t'))
    acetyl = list(csv.DictReader(open(folder / 'human_acetylation_sites_v1.5.0.tsv', encoding='utf-8'), delimiter='\t'))
    assert 'state' not in acyl[0] and 'acyl_type' in acyl[0]
    assert [r for r in acyl if r['acyl_type'] == 'acetyl'] == acetyl
    assert len({(r['canonical_UniProtAC'], r['position'], r['acyl_type']) for r in acyl}) == len(acyl)


def read_release(name):
    import pandas as pd
    return pd.read_csv(ROOT / 'releases' / 'v1.5.0' / name, sep='\t', dtype=str, keep_default_na=False)


def read_staged(name):
    import pandas as pd
    return pd.read_csv(ROOT / 'inputs' / 'marks' / name, sep='\t', dtype=str, keep_default_na=False)


def test_berryhill_rows_in_the_release():
    reference = stage.parse_fasta(ROOT / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03.fasta')
    table = read_release('human_methyllysine_sites_v1.5.0.tsv')
    staged = read_staged('berryhill2023_methyllysine_sites.tsv')
    bh = table[table['sources'].str.split('; ').map(lambda x: 'Berryhill 2023' in x)]
    assert len(staged) == len(bh) == 5072 and bh['canonical_UniProtAC'].nunique() == staged['canonical_UniProtAC'].nunique()
    assert not table.duplicated(['canonical_UniProtAC', 'position', 'state']).any()
    assert set(bh['state']) == {'me1', 'me2', 'me3'} and set(bh['residue']) == {'K'} and set(bh['mark']) == {'methyllysine'}
    assert (bh['evidence_class'] == 'experimental').all()
    for r in bh.itertuples():
        assert reference[r.canonical_UniProtAC][int(r.position) - 1] == 'K'
    merged = bh.set_index(['canonical_UniProtAC', 'position', 'state'])
    for r in staged.itertuples():
        row = merged.loc[(r.canonical_UniProtAC, r.position, r.state)]
        assert (row['berryhill_site_probability'], row['berryhill_cell_lines'], row['berryhill_enrichment']) == (
            r.berryhill_site_probability, r.berryhill_cell_lines, r.berryhill_enrichment)
        assert 'Sci Rep 13:377 (2023)' in row['source_releases']
    other = table[~table.index.isin(bh.index)]
    assert set(other['berryhill_site_probability'] + other['berryhill_cell_lines'] + other['berryhill_enrichment']) == {''}
    # a row no Berryhill row touches is what the other sources give: the merge without Berryhill has the same rows and values
    base = build.merge([read_staged(f.name) for f in sorted((ROOT / 'inputs' / 'marks').glob('*_methyllysine_sites.tsv'))
                        if not f.name.startswith('berryhill')]).set_index(['canonical_UniProtAC', 'position', 'state'])
    untouched = table[~table.set_index(['canonical_UniProtAC', 'position', 'state']).index.isin(staged.set_index(['canonical_UniProtAC', 'position', 'state']).index)]
    assert len(untouched) == len(base) - base.index.isin(staged.set_index(['canonical_UniProtAC', 'position', 'state']).index).sum()
    cols = build.OUT_COLUMNS + ['cplm_evidence']
    a = untouched[cols].set_index(['canonical_UniProtAC', 'position', 'state']).sort_index()
    b = base.reset_index().set_index(['canonical_UniProtAC', 'position', 'state']).loc[a.index, [c for c in cols if c not in ('canonical_UniProtAC', 'position', 'state')]].sort_index()
    assert a.equals(b)


def test_ochoa_rows_in_the_release():
    reference = stage.parse_fasta(ROOT / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03.fasta')
    table = read_release('human_phosphorylation_sites_v1.5.0.tsv')
    staged = read_staged('ochoa2020_phosphorylation_sites.tsv')
    oc = table[table['sources'].str.split('; ').map(lambda x: 'Ochoa 2020' in x)]
    assert len(staged) == len(oc) == 112405 and oc['canonical_UniProtAC'].nunique() == 10275
    assert not table.duplicated(['canonical_UniProtAC', 'position', 'state']).any() and set(table['state']) == {'phospho'}
    assert set(oc['residue']) <= {'S', 'T', 'Y'} and (oc['evidence_class'] == 'experimental').all()
    for r in oc.itertuples():
        assert reference[r.canonical_UniProtAC][int(r.position) - 1] == r.residue
    merged = oc.set_index(['canonical_UniProtAC', 'position'])
    for r in staged.sample(2000, random_state=1).itertuples():
        row = merged.loc[(r.canonical_UniProtAC, r.position)]
        assert [row[c] for c in build.SOURCE_COLUMNS['Ochoa 2020']] == [getattr(r, c) for c in build.SOURCE_COLUMNS['Ochoa 2020']]
    assert (oc['functional_score'] != '').all() and ((oc['functional_score'].astype(float) >= 0) & (oc['functional_score'].astype(float) <= 1)).all()
    other = table[~table.index.isin(oc.index)]
    assert set(other[build.SOURCE_COLUMNS['Ochoa 2020']].sum(axis=1)) == {''}
    assert 'Ochoa 2020' in set(table['sources'].str.split('; ').explode())


def test_iptmnet_score_in_the_release():
    staged = read_staged('iptmnet_6_2_score.tsv')
    for name, mark in (('human_methyllysine_sites_v1.5.0.tsv', 'methyllysine'), ('human_phosphorylation_sites_v1.5.0.tsv', 'phosphorylation'),
                       ('human_lysine_acylation_sites_v1.5.0.tsv', 'lysine_acylation'), ('human_acetylation_sites_v1.5.0.tsv', 'lysine_acylation'),
                       ('human_ubiquitination_sites_v1.5.0.tsv', 'ubiquitination'), ('human_sumoylation_sites_v1.5.0.tsv', 'sumoylation')):
        table = read_release(name)
        given = dict(zip(zip(staged['canonical_UniProtAC'][staged['mark'] == mark], staged['position'][staged['mark'] == mark], strict=True),
                         staged['iptmnet_score'][staged['mark'] == mark], strict=True))
        has_iptmnet = table['sources'].str.split('; ').map(lambda x: 'iPTMnet' in x)
        assert set(table.loc[~has_iptmnet, 'iptmnet_score']) <= {''}, name
        assert set(table['iptmnet_score']) <= {'', '0', '1', '2', '3', '4'}
        scored = table[table['iptmnet_score'] != '']
        assert len(scored) > 0.5 * has_iptmnet.sum(), name
        for r in scored.itertuples():
            assert given[(r.canonical_UniProtAC, r.position)] == r.iptmnet_score
        assert 'iptmnet_score' in table.columns
    ac, ph = read_release('human_acetylation_sites_v1.5.0.tsv'), read_release('human_lysine_acylation_sites_v1.5.0.tsv')
    assert ac.equals(ph[ph['acyl_type'] == 'acetyl'].reset_index(drop=True))
    # the score columns are in the data dictionary, and the methylarginine files have none
    dictionary = (ROOT / 'releases' / 'v1.5.0' / 'data_dictionary.tsv').read_text(encoding='utf-8')
    assert 'human_phosphorylation_sites_v1.5.0.tsv\tiptmnet_score' in dictionary and 'human_methylarginine_crosscheck_v1.5.0.tsv\tiptmnet_score' not in dictionary
    for col in ('functional_score', 'best_pep', 'berryhill_cell_lines', 'cplm_evidence'):
        assert col in dictionary


def test_license_notes_of_the_new_sources_in_the_release_text():
    folder = ROOT / 'releases' / 'v1.5.0'
    for name in ('README.md', 'CHANGES.md'):
        text = (folder / name).read_text(encoding='utf-8')
        for needle in ('Sci Rep 13:377', '10.1038/s41598-022-27175-x', 'CC BY 4.0', 'supplements carry no separate',
                       'Nat Biotechnol 38:365', '10.1038/s41587-019-0344-3', 'not open access', 'exclusive licence to Springer Nature America',
                       'The score is calculated based on all the evidence sources as shown on the website',
                       'Nat Struct Mol Biol 24:325-336', '10.1038/nsmb.3366', 'Nat Struct Mol Biol 25:631-640', '10.1038/s41594-018-0084-y'):
            assert needle in text, (name, needle)
        for phrase in ('pending permission', 'permission from', "owner's decision", "authors' permission", "authors' confirmation", 'when they answer'):
            assert phrase not in text, (name, phrase)
    release = json.loads((folder / 'release.json').read_text(encoding='utf-8'))
    assert {'Berryhill 2023', 'Ochoa 2020', 'Hendriks 2017', 'Akimov 2018'} <= set(release['sources']) and 'iptmnet_score' in release['inputs_md5']


def crosslink(term, pos, **kw):
    return {**feature(term, pos, **kw), 'type': 'Cross-link'}


def test_uniprot_crosslink_terms_ubiquitin_and_sumo_on_lysine():
    ub, sumo = 'Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in ubiquitin)', 'Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in SUMO2)'
    assert stage.crosslink_term(ub) == ('ubiquitination', 'ub', 'K') and stage.crosslink_term(ub + '; alternate') == ('ubiquitination', 'ub', 'K')
    for modifier in ('SUMO', 'SUMO1', 'SUMO3', 'SUMO1P1/SUMO5', 'SUMO2 and SUMO3'):
        assert stage.crosslink_term(sumo.replace('SUMO2', modifier)) == ('sumoylation', 'sumo', 'K'), modifier
    for term in ('Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in ISG15)', 'Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in NEDD8)',
                 '(Microbial infection) ' + ub, 'Glycyl cysteine thioester (Cys-Gly) (interchain with G-Cter in ubiquitin)',
                 'Peptide (Met-Gly) (interchain with G-Cter in ubiquitin)', 'Isoglutamyl lysine isopeptide (Lys-Gln) (interchain with Q-?)'):
        assert stage.crosslink_term(term) is None, term
    entries = [entry('P00001', REFERENCE['P00001'], crosslink(ub, 2, pmid='11'), crosslink(sumo, 2, codes=('ECO:0000250',)), crosslink(sumo, 2, pmid='12'),
                     crosslink(ub, 4), crosslink(sumo.replace('SUMO2', 'NEDD8'), 3), crosslink(ub, 5)),
               entry('P00002', 'MKAKX', crosslink(ub, 2))]
    tables, dropped = run(entries)
    u, m = tables['ubiquitination'], tables['sumoylation']
    assert list(zip(u['canonical_UniProtAC'], u['position'], u['state'], u['pmids'])) == [('P00001', 2, 'ub', '11')]
    assert dropped == {'residue_mismatch': 2, 'sequence_differs': 1}      # position 4 is R, position 5 is S, P00002's sequence is not the reference's
    assert list(m['state']) == ['sumo'] and list(m['evidence_class']) == ['experimental'] and m.loc[0, 'evidence'] == 'ECO:0000250;ECO:0000269'
    assert not tables['methyllysine'].shape[0] and tables['ubiquitination']['residue'].eq('K').all()


def test_iptmnet_ubiquitination_and_sumoylation_records():
    import pandas as pd
    human = 'Homo sapiens (Human)'
    rows = [('UBIQUITINATION', 'hprd', 'P00001', 'G', human, 'K2', '', '', '', '5'), ('UBIQUITINATION', 'unip', 'P00001', 'G', human, 'K2', '', '', '', ''),
            ('SUMOYLATION', 'sign', 'P00001-1', 'G', human, 'K3', '', '', '', '7'),
            ('UBIQUITINATION', 'hprd', 'P00001', 'G', human, 'S5', '', '', '', '5'), ('UBIQUITINATION', 'hprd', 'P00001', 'G', human, 'K4', '', '', '', '5'),
            ('UBIQUITINATION', 'hprd', 'P00001', 'G', 'Mus musculus (Mouse)', 'K2', '', '', '', '5')]
    tables, dropped = iptm.stage(pd.DataFrame(rows, columns=iptm.FIELDS), REFERENCE, '6.2')
    u, m = tables['ubiquitination'], tables['sumoylation']
    assert list(zip(u['position'], u['state'], u['upstream_sources'], u['evidence_class'], u['pmids'])) == [(2, 'ub', 'hprd;unip', 'experimental', '5')]
    assert list(zip(m['position'], m['state'], m['upstream_sources'])) == [(3, 'sumo', 'sign')]
    assert dropped == {'other_residue': 1, 'residue_mismatch': 1, 'not_human': 1}


def test_iptmnet_score_marks_ubiquitination_and_sumoylation():
    assert score.mark_of('Ubiquitination', 'K') == 'ubiquitination' and score.mark_of('Sumoylation', 'K') == 'sumoylation'
    assert score.mark_of('Ubiquitination', 'S') is None and score.mark_of('Sumoylation', 'R') is None


def test_cplm_modification_option_keeps_methylation_names_and_gives_others_their_own():
    assert cplm.DEFAULT_MODIFICATION == 'Methylation' and set(cplm.MODIFICATIONS) == {'Methylation', 'Ubiquitination', 'Sumoylation'} | set(cplm.ACYL_TYPES)
    me, ub = cplm.MODIFICATIONS['Methylation'], cplm.MODIFICATIONS['Ubiquitination']
    assert (me['table'], me['summary'], me['mark'], me['state'], me['tag']) == ('cplm_4_0_methyllysine_sites.tsv', 'cplm_ptm_sites_summary.json', 'methyllysine', 'me', 'methyl')
    assert (ub['table'], ub['mark'], ub['state']) == ('cplm_4_0_ubiquitination_sites.tsv', 'ubiquitination', 'ub')
    import pandas as pd
    good = pd.DataFrame({'canonical_UniProtAC': ['P00001'], 'position': [2], 'evidence': ['Exp.'], 'references': ['31072472']})
    table = cplm.stage_marks(good, 'CPLM 4.0 (2021-06-13)', ub)
    assert (table.loc[0, 'mark'], table.loc[0, 'state'], table.loc[0, 'residue']) == ('ubiquitination', 'ub', 'K')
    assert cplm.stage_marks(good, 'CPLM 4.0 (2021-06-13)').loc[0, 'mark'] == 'methyllysine'


def test_ubiquitination_and_sumoylation_tables_of_the_release():
    import pandas as pd
    folder = ROOT / 'releases' / 'v1.5.0'
    reference = stage.parse_fasta(ROOT / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03.fasta')
    for mark, state, sources in (('ubiquitination', 'ub', {'UniProt', 'iPTMnet', 'CPLM'}), ('sumoylation', 'sumo', {'UniProt', 'iPTMnet', 'CPLM'})):
        table = pd.read_csv(folder / f'human_{mark}_sites_v1.5.0.tsv', sep='\t', dtype=str, keep_default_na=False)
        assert set(table['state']) == {state} and set(table['mark']) == {mark} and set(table['residue']) == {'K'}
        assert not table.duplicated(['canonical_UniProtAC', 'position', 'state']).any()
        assert {x for v in table['sources'] for x in v.split('; ')} == sources
        for r in table.itertuples():
            assert reference[r.canonical_UniProtAC][int(r.position) - 1] == 'K'
        by_source = {s: table[table['sources'].str.contains(s)] for s in sources}
        staged = {'UniProt': f'uniprot_2026_04_{mark}_sites.tsv', 'iPTMnet': f'iptmnet_6_2_{mark}_sites.tsv', 'CPLM': f'cplm_4_0_{mark}_sites.tsv'}
        for s in sources:
            n = len(pd.read_csv(ROOT / 'inputs' / 'marks' / staged[s], sep='\t', dtype=str, keep_default_na=False))
            assert len(by_source[s]) == n, (mark, s)
    text = (folder / 'README.md').read_text(encoding='utf-8') + (folder / 'CHANGES.md').read_text(encoding='utf-8')
    for needle in ('All Rights Reserved', 'Zhang et al., Nucleic Acids Res 50:D451',
                   'CPLM states no license', 'The sumoylation table likewise includes'):
        assert needle in text, needle
    dictionary = (folder / 'data_dictionary.tsv').read_text(encoding='utf-8')
    for name in ('human_ubiquitination_sites_v1.5.0.tsv', 'human_sumoylation_sites_v1.5.0.tsv'):
        assert f'{name}\tiptmnet_score' in dictionary
    assert 'human_ubiquitination_sites_v1.5.0.tsv\tcplm_evidence' in dictionary and 'human_sumoylation_sites_v1.5.0.tsv\tcplm_evidence' in dictionary
    assert 'SUMOylation sites are not included' not in text


ACYL_COUNTS = {   # CPLM 4.0 acylation downloads: acyl type -> sites that passed validation (human lysines on the reference sequence)
    'acetyl': 52493, 'succinyl': 4130, 'crotonyl': 12021, 'malonyl': 4978, 'glutaryl': 49, '2-hydroxyisobutyryl': 8760,
    'beta-hydroxybutyryl': 2932, 'lactoyl': 71, 'propionyl': 44, 'butyryl': 192, 'benzoyl': 22, 'biotinyl': 9, 'lipoyl': 5, 'phosphoglyceryl': 38}


def test_cplm_acylation_staging_maps_every_download_to_one_acyl_type():
    acyl = {n: c for n, c in cplm.MODIFICATIONS.items() if c['mark'] == 'lysine_acylation'}
    assert set(acyl) == set(cplm.ACYL_TYPES) and {c['state'] for c in acyl.values()} == set(ACYL_COUNTS) and len(acyl) == 14
    assert len({c['table'] for c in acyl.values()}) == len({c['summary'] for c in acyl.values()}) == 14
    assert all(c['table'].endswith('_lysine_acylation_sites.tsv') for c in acyl.values())
    # the methylation, ubiquitination and sumoylation downloads keep their names and states
    assert [cplm.MODIFICATIONS[m]['table'] for m in ('Methylation', 'Ubiquitination', 'Sumoylation')] == [
        'cplm_4_0_methyllysine_sites.tsv', 'cplm_4_0_ubiquitination_sites.tsv', 'cplm_4_0_sumoylation_sites.tsv']
    # UniProt's names for the types that both have
    uniprot_types = {'acetyl', 'succinyl', '2-hydroxyisobutyryl', 'beta-hydroxybutyryl', 'lactoyl', 'glutaryl', 'crotonyl', 'benzoyl', 'butyryl', 'malonyl', 'propionyl'}
    assert uniprot_types <= set(ACYL_COUNTS)


def test_cplm_acylation_rows_in_the_release_by_type():
    import json

    import pandas as pd
    folder = ROOT / 'releases' / 'v1.5.0'
    table = read_release('human_lysine_acylation_sites_v1.5.0.tsv')
    acetyl = read_release('human_acetylation_sites_v1.5.0.tsv')
    reference = stage.parse_fasta(ROOT / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03.fasta')
    assert 'cplm_evidence' in table.columns and 'cplm_evidence' in acetyl.columns and list(table.columns)[-2:] == ['cplm_evidence', 'iptmnet_score']
    assert not table.duplicated(['canonical_UniProtAC', 'position', 'acyl_type']).any()
    assert set(table['acyl_type']) == set(ACYL_COUNTS) | {'methacryl', 'isonicotinyl'}
    cp = table[table['sources'].str.split('; ').map(lambda x: 'CPLM' in x)]
    assert set(cp['acyl_type']) == set(ACYL_COUNTS) and set(cp['cplm_evidence']) <= {'Exp.', 'Exp./Dat.', 'Dat.'} and '' not in set(cp['cplm_evidence'])
    assert set(table.loc[~table.index.isin(cp.index), 'cplm_evidence']) == {''}
    for r in table.itertuples():     # a lysine at the position, on the reference sequence
        assert r.residue == 'K' and reference[r.canonical_UniProtAC][int(r.position) - 1] == 'K'
    for name, c in cplm.MODIFICATIONS.items():
        if c['mark'] != 'lysine_acylation':
            continue
        staged = read_staged(c['table'])
        summary = json.loads((ROOT / 'inputs' / 'marks' / c['summary']).read_text(encoding='utf-8'))
        kind = cp[cp['acyl_type'] == c['state']]
        assert set(staged['state']) == {c['state']} and set(staged['source']) == {'CPLM'} and set(staged['mark']) == {'lysine_acylation'}
        assert len(staged) == len(kind) == summary['rows']['lysine_acylation'] and kind['canonical_UniProtAC'].nunique() == summary['proteins']
        assert summary['records_in_scope'] == len(staged) + sum(summary['dropped'].values()), name
        assert not kind.duplicated(['canonical_UniProtAC', 'position']).any() and len(staged) == ACYL_COUNTS[c['state']]
    assert acetyl.equals(table[table['acyl_type'] == 'acetyl'].reset_index(drop=True))
    assert 'CPLM states no license' in (folder / 'CHANGES.md').read_text() or 'CPLM states no license' in (folder / 'README.md').read_text()
    release = json.loads((folder / 'release.json').read_text(encoding='utf-8'))
    assert release['context_marks']['lysine_acylation']['rows'] == len(table)
    assert len(pd.read_csv(folder / 'human_acetylation_sites_v1.5.0.tsv', sep='\t', dtype=str, keep_default_na=False)) == len(acetyl)

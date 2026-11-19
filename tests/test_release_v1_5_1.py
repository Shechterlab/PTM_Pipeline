"""Release v1.5.1: complete, checksummed, the methylarginine data byte for byte those of v1.4.3, the context marks those of
v1.5.0, every site set's interaction tables, and the counts its README and release.json give."""
import hashlib
import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
REL = ROOT / 'releases' / 'v1.5.1'
pytestmark = pytest.mark.skipif(not (REL / 'release.json').exists(), reason='releases/v1.5.1 not built')


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def load(name):
    spec = importlib.util.spec_from_file_location(name.replace('.py', ''), ROOT / 'scripts' / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def release(version):
    return json.loads((ROOT / 'releases' / version / 'release.json').read_text(encoding='utf-8'))


@pytest.mark.parametrize('folder', ['.', 'analyses', 'analyses_with_lower_confidence'])
def test_checksums(folder):
    base = REL / folder
    for line in (base / 'MD5SUMS').read_text(encoding='utf-8').splitlines():
        digest, name = line.split('  ', 1)
        assert md5(base / name) == digest, f'{folder}/{name}'
    if folder == '.':
        listed = {ln.split('  ', 1)[1] for ln in (base / 'MD5SUMS').read_text(encoding='utf-8').splitlines()}
        files = {p.relative_to(base).as_posix() for p in base.rglob('*') if p.is_file() and p.name != 'MD5SUMS'}
        assert listed == files


def test_methylarginine_data_are_those_of_v1_4_3():
    """One-time check (DECISIONS_v1.5.1 D3): every methylarginine file v1.5.1 carries is, byte for byte, the v1.4.3 file it
    comes from, apart from the files that name their own version (STUDIES.md, the analyses READMEs and MD5SUMS, the workbook)."""
    r151, r150 = release('v1.5.1'), release('v1.5.0')
    edited = set(r151['carried_over_text_edited']) | {n.replace('_v1.5.0', '_v1.5.1') for n in r150['carried_over_text_edited']}
    checked = 0
    for name, info in r151['carried_over'].items():
        from150 = info['from']
        if from150 not in r150['carried_over'] or name in edited:
            continue
        v143 = ROOT / 'releases' / 'v1.4.3' / r150['carried_over'][from150]['from']
        assert md5(REL / name) == md5(v143), name
        checked += 1
    assert checked >= 90


def test_context_marks_are_those_of_v1_5_0():
    for name in ('methyllysine', 'phosphorylation', 'lysine_acylation', 'acetylation', 'ubiquitination', 'sumoylation'):
        assert md5(REL / f'human_{name}_sites_v1.5.1.tsv') == md5(ROOT / 'releases' / 'v1.5.0' / f'human_{name}_sites_v1.5.0.tsv'), name
    assert md5(REL / 'human_methylarginine_crosscheck_v1.5.1.tsv') == md5(ROOT / 'releases' / 'v1.5.0' / 'human_methylarginine_crosscheck_v1.5.0.tsv')


def test_every_site_set_has_its_interaction_tables_once():
    config = json.loads((ROOT / 'config' / 'site_sets.json').read_text(encoding='utf-8'))
    tables = load('65_carry_release.py').INTERACTION_TABLES
    for name in config['site_sets']:
        for t in tables:
            assert (REL / 'supporting' / f'{name}_{t}_v1.5.1.tsv').exists(), (name, t)
    for t in tables:   # the methylarginine tables under the site set's name only, the v1.4.3 bytes
        assert not (REL / 'supporting' / f'{t}_v1.5.1.tsv').exists()
        assert md5(REL / 'supporting' / f'methylarginine_{t}_v1.5.1.tsv') == md5(ROOT / 'releases' / 'v1.4.3' / 'supporting' / f'{t}_v1.4.3.tsv')


def test_no_cplm_in_the_release_tables():
    for path in list(REL.glob('human_*_sites_v1.5.1.tsv')):
        t = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, nrows=None)
        if 'sources' in t:
            assert not t['sources'].str.contains('CPLM').any(), path.name
        assert 'cplm_evidence' not in t.columns, path.name
    assert 'CPLM' not in release('v1.5.1')['sources'] and set(release('v1.5.1')['sources_excluded']) == {'CPLM'}
    readme = (REL / 'README.md').read_text(encoding='utf-8')
    assert '| CPLM' not in readme and 'CPLM 4.0 lysine methylation: Zhang' not in readme   # not listed as a layer or cited


def test_counts_in_release_json_and_readme():
    x65 = load('65_carry_release.py')
    config = json.loads((ROOT / 'config' / 'site_sets.json').read_text(encoding='utf-8'))
    rows = x65.count_marks(REL, 'v1.5.1', 'v1.5.0', config)
    assert rows == release('v1.5.1')['counts']
    by = {r['set']: r for r in rows}
    assert by['methylarginine']['sites'] == 11632 and by['phospho']['sites'] == 142415
    assert by['methyllysine']['sites'] == 5178 and by['lysine']['sites'] == 90876
    assert all(r['high_confidence'] <= r['sites'] and r['with_interaction_evidence'] <= r['sites'] for r in rows)
    readme = (REL / 'README.md').read_text(encoding='utf-8')
    for r in rows:
        assert f'| {r["sites"]:,} | {r["high_confidence"]:,} | {r["with_interaction_evidence"]:,}' in readme, r['set']


def test_high_confidence_rule():
    x65 = load('65_carry_release.py')
    t = pd.DataFrame({'canonical_UniProtAC': ['P1', 'P1', 'P2', 'P2', 'P3'], 'position': ['2', '2', '5', '5', '7'],
                      'evidence_class': ['experimental', 'experimental', 'experimental', 'by-similarity', 'experimental'],
                      'pmids': ['111', '', '111', '222', '111;222'], 'sources': ['UniProt', 'Hendriks 2017', 'UniProt', 'UniProt', 'iPTMnet']})
    units = x65.mark_units(t)
    # P1: PubMed 111 and Hendriks 2017 (its PubMed ID 28112733) -> 2; P2: the by-similarity row adds nothing -> 1; P3: 2 PubMed IDs
    assert units[('P1', '2')] == 2 and units[('P2', '5')] == 1 and units[('P3', '7')] == 2
    t2 = pd.DataFrame({'canonical_UniProtAC': ['P1'], 'position': ['2'], 'evidence_class': ['experimental'], 'pmids': ['28112733'], 'sources': ['UniProt; Hendriks 2017']})
    assert x65.mark_units(t2)[('P1', '2')] == 1     # a database citing Hendriks 2017 is not counted again
    # per modification (the lysine set): ubiquitin in one study and SUMO in another is one unit each; two studies of SUMO is two
    t3 = pd.DataFrame({'canonical_UniProtAC': ['P1', 'P1', 'P2', 'P2'], 'position': ['2', '2', '9', '9'], 'state': ['ub', 'sumo', 'sumo', 'sumo'],
                       'evidence_class': ['experimental'] * 4, 'pmids': ['', '', '111', ''], 'sources': ['Akimov 2018', 'Hendriks 2017', 'UniProt', 'Hendriks 2017']})
    pooled, per_mod = x65.mark_units(t3), x65.mark_units(t3, per_modification=True)
    assert pooled[('P1', '2')] == 2 and per_mod[('P1', '2')] == 1 and per_mod[('P2', '9')] == 2
    t4 = pd.DataFrame({'canonical_UniProtAC': ['P1', 'P1'], 'position': ['2', '2'], 'acyl_type': ['acetyl', 'succinyl'],
                       'evidence_class': ['experimental'] * 2, 'pmids': ['111', '222'], 'sources': ['UniProt', 'UniProt']})
    assert x65.mark_units(t4, per_modification=True)[('P1', '2')] == 1     # two acyl types are two modifications
    s = pd.DataFrame({'studies': ['Larsen2016; Sylvestersen2014', 'Larsen2016; Guo2014', 'PSP; UniProt; Guo2014', 'Guo2014; Fedoriw2019; PSP']})
    assert x65.arginine_high_confidence(s) == 2     # one laboratory; two studies; databases not counted; two studies


def test_american_spelling_keeps_identifiers():
    x65 = load('65_carry_release.py')
    text = 'the localisation of a tumour, the `localisation` column, "a localised title", [inputs/licences/](../../inputs/licences/a.pdf)'
    assert x65.americanize(text) == ('the localization of a tumor, the `localisation` column, "a localised title", '
                                     '[inputs/licences/](../../inputs/licences/a.pdf)')
    for name in ('README.md', 'CHANGES.md'):
        prose = x65._PROTECT.sub('', (REL / name).read_text(encoding='utf-8'))
        assert not x65._WORD.search(prose), name

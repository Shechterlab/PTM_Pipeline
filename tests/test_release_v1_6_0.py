"""Release v1.6.0: v1.5.1 carried forward unchanged, and the somatic changes at every candidate residue with the site sets
each residue belongs to (scripts/68_release_cancer_changes.py)."""
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
REL = ROOT / 'releases' / 'v1.6.0'
PREV = ROOT / 'releases' / 'v1.5.1'
LAYERS = ROOT / 'inputs' / 'layers'
RESIDUES = ['R', 'K', 'S', 'T', 'Y']
pytestmark = pytest.mark.skipif(not (REL / 'release.json').exists(), reason='releases/v1.6.0 not built')


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


def load(name):
    spec = importlib.util.spec_from_file_location(name.replace('.py', ''), ROOT / 'scripts' / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


RELEASE = json.loads((REL / 'release.json').read_text(encoding='utf-8')) if (REL / 'release.json').exists() else {}


@pytest.mark.parametrize('folder', ['.', 'analyses', 'analyses_with_lower_confidence'])
def test_checksums(folder):
    base = REL / folder
    lines = (base / 'MD5SUMS').read_text(encoding='utf-8').splitlines()
    for line in lines:
        digest, name = line.split('  ', 1)
        assert md5(base / name) == digest, f'{folder}/{name}'
    if folder == '.':
        listed = {ln.split('  ', 1)[1] for ln in lines}
        assert listed == {p.relative_to(base).as_posix() for p in base.rglob('*') if p.is_file() and p.name != 'MD5SUMS'}


def test_previous_release_is_carried_unchanged():
    """Every file of v1.5.1 but its own documents and site changes is here under the new name, the same bytes apart from the
    files that name their version."""
    edited = set(RELEASE['carried_over_text_edited'])
    assert edited <= {'STUDIES.md', 'analyses/README.md', 'analyses/MD5SUMS', 'analyses_with_lower_confidence/README.md',
                      'analyses_with_lower_confidence/MD5SUMS'} | {n for n in edited if n.endswith('.xlsx')}
    n = 0
    for path in PREV.rglob('*'):
        rel = path.relative_to(PREV).as_posix()
        if not path.is_file() or ('/' not in rel and rel in {'README.md', 'CHANGES.md', 'MD5SUMS', 'data_dictionary.tsv', 'release.json'}) \
                or path.name.startswith('site_changes_'):
            continue
        new = REL / rel.replace('_v1.5.1', '_v1.6.0')
        assert new.exists(), rel
        if rel.replace('_v1.5.1', '_v1.6.0') not in edited:
            assert md5(new) == md5(path), rel
            n += 1
    assert n > 100


@pytest.mark.parametrize('res', RESIDUES)
def test_change_tables_are_the_staged_layer_with_site_sets(res):
    for src in ('tcga_mc3', 'depmap_24q4'):
        staged = read(LAYERS / (f'{src}_arginine_changes.tsv' if res == 'R' else f'{src}_changes_{res}.tsv.gz'))
        t = read(REL / 'supporting' / f'{src}_changes_{res}_v1.6.0.tsv.gz')
        assert list(t.columns) == list(staged.columns) + ['site_sets']
        assert t.drop(columns='site_sets').equals(staged), (src, res)
        with open(REL / 'supporting' / f'{src}_changes_{res}_v1.6.0.tsv.gz', 'rb') as fh:
            assert fh.read(10)[4:8] == b'\x00\x00\x00\x00'     # gzip without a timestamp


def test_site_sets_name_the_sites_of_each_set():
    config = json.loads((ROOT / 'config' / 'site_sets.json').read_text(encoding='utf-8'))
    x65 = load('65_carry_release.py')
    # methylarginine: the sites with a change are those the site table itself counts (arginine layer columns)
    sites = read(REL / 'human_methylarginine_sites_v1.6.0.tsv')
    for src, col in (('tcga_mc3', 'tcga_changes'), ('depmap_24q4', 'depmap_changes')):
        t = read(REL / 'supporting' / f'{src}_changes_R_v1.6.0.tsv.gz')
        assert set(t['site_sets']) <= {'', 'methylarginine'}
        changed = set(zip(t.loc[t['site_sets'] == 'methylarginine', 'canonical_UniProtAC'], t.loc[t['site_sets'] == 'methylarginine', 'position']))
        expected = set(zip(sites.loc[sites[col].ne('') & sites[col].ne('0'), 'uniprot_accession'], sites.loc[sites[col].ne('') & sites[col].ne('0'), 'position']))
        assert changed == expected, src
    # the context marks: a row names a set exactly when that set's tables have the residue
    for name in ('methyllysine', 'lysine', 'phospho'):
        spec = config['site_sets'][name]
        known = set()
        for table in spec['tables']:
            s = read(REL / x65.release_name(table, 'v1.6.0'))
            s = s[s['residue'].isin(spec['residues'])]
            known |= set(zip(s['canonical_UniProtAC'], s['position'], s['residue']))
        for res in spec['residues']:
            t = read(REL / 'supporting' / f'tcga_mc3_changes_{res}_v1.6.0.tsv.gz')
            flagged = t['site_sets'].str.split(';').map(lambda v: name in v)
            assert all(((a, p, res) in known) == f for a, p, f in zip(t['canonical_UniProtAC'], t['position'], flagged)), (name, res)


def test_counts_in_release_json_readme_and_dictionary():
    readme = (REL / 'README.md').read_text(encoding='utf-8')
    per_set = {r['set']: r for r in RELEASE['cancer_changes']['per_site_set']}
    assert per_set['methylarginine']['tcga_mc3_sites_changed'] == 3100 and per_set['methylarginine']['depmap_24q4_sites_changed'] == 849
    for r in per_set.values():
        assert f'| {r["label"]} | {r["sites"]:,} | {r["tcga_mc3_sites_changed"]:,} | {r["depmap_24q4_sites_changed"]:,} |' in readme
    for key, s in RELEASE['cancer_changes']['per_residue'].items():
        src, res = key.split('/')
        assert s['rows'] == len(read(REL / 'supporting' / f'{src}_changes_{res}_v1.6.0.tsv.gz'))
    dd = read(REL / 'data_dictionary.tsv')
    for f in RELEASE['cancer_changes']['files']:
        path = REL / f
        opener = gzip.open if f.endswith('.gz') else open
        with opener(path, 'rt', encoding='utf-8') as fh:
            cols = fh.readline().rstrip('\n').split('\t')
        assert list(dd.loc[dd['file'] == f, 'column']) == cols, f
    assert not dd.duplicated(['file', 'column']).any()
    assert RELEASE['version'] == 'v1.6.0' and RELEASE['previous_version'] == 'v1.5.1'

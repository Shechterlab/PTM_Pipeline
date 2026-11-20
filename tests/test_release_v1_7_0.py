"""Release v1.7.0: v1.6.0 carried forward unchanged, and the per-site mutation propensity of its change tables
(scripts/71_release_propensity.py, scripts/propensity.py)."""
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
REL = ROOT / 'releases' / 'v1.7.0'
PREV = ROOT / 'releases' / 'v1.6.0'
LAYERS = ROOT / 'inputs' / 'layers'
RESIDUES = ['R', 'K', 'S', 'T', 'Y']
SOURCES = {'tcga_mc3': 'sample_id', 'depmap_24q4': 'model_id', 'clinvar': 'variation_id', 'gnomad_v4_1': 'variant_id'}
pytestmark = pytest.mark.skipif(not (REL / 'release.json').exists(), reason='releases/v1.7.0 not built')
RELEASE = json.loads((REL / 'release.json').read_text(encoding='utf-8')) if (REL / 'release.json').exists() else {}


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def read(path: Path, **kw) -> pd.DataFrame:
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, **kw)


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
    edited = set(RELEASE['carried_over_text_edited'])
    assert edited <= {'STUDIES.md', 'analyses/README.md', 'analyses/MD5SUMS', 'analyses_with_lower_confidence/README.md',
                      'analyses_with_lower_confidence/MD5SUMS'} | {n for n in edited if n.endswith('.xlsx')}
    n = 0
    for path in PREV.rglob('*'):
        rel = path.relative_to(PREV).as_posix()
        if not path.is_file() or ('/' not in rel and rel in {'README.md', 'CHANGES.md', 'MD5SUMS', 'data_dictionary.tsv', 'release.json'}) \
                or path.name.startswith('site_changes_'):
            continue
        new = REL / rel.replace('_v1.6.0', '_v1.7.0')
        assert new.exists(), rel
        if rel.replace('_v1.6.0', '_v1.7.0') not in edited:
            assert md5(new) == md5(path), rel
            n += 1
    assert n > 140


@pytest.fixture(scope='module')
def tables():
    return {src: read(REL / 'supporting' / f'site_propensity_{src}_v1.7.0.tsv.gz') for src in SOURCES}


@pytest.mark.parametrize('src', list(SOURCES))
def test_rows_are_known_sites_and_changed_residues_with_their_counts(src, tables):
    t = tables[src]
    for res in RESIDUES:
        ch = read(REL / 'supporting' / f'{src}_changes_{res}_v1.7.0.tsv.gz', usecols=['canonical_UniProtAC', 'position', SOURCES[src], 'site_sets'])
        k = ch.groupby(['canonical_UniProtAC', 'position'])[SOURCES[src]].nunique()
        sub = t[t['residue'] == res]
        got = dict(zip(zip(sub['canonical_UniProtAC'], sub['position']), sub['k'].astype(int)))
        # every changed residue is a row, with its count of distinct samples, lines, variations or variants
        assert all(got.get(key) == v for key, v in k.items()), (src, res)
        # the other rows are known sites without a change
        rest = sub[sub['k'] == '0']
        assert (rest['site_sets'] != '').all(), (src, res)
        # n is the sum of k over the protein's residues of the type (every residue with a change is written)
        tested = sub[sub['n'] != '']
        sums = sub.groupby('canonical_UniProtAC')['k'].apply(lambda s: s.astype(int).sum())
        assert (tested['n'].astype(int).to_numpy() == tested['canonical_UniProtAC'].map(sums).to_numpy()).all(), (src, res)
        assert len(sub) == RELEASE['propensity']['per_source_residue'][f'{src}/{res}']['rows']
    # known sites: every site of the site sets at this source's candidate residues is a row
    with gzip.open(REL / 'supporting' / f'site_propensity_{src}_v1.7.0.tsv.gz', 'rb') as fh:
        assert fh.read(1)
    assert (REL / 'supporting' / f'site_propensity_{src}_v1.7.0.tsv.gz').read_bytes()[4:8] == b'\x00\x00\x00\x00'


def test_expected_and_tests_are_consistent(tables):
    for src, t in tables.items():
        tested = t[t['p_upper'] != '']
        n, e, k = tested['n'].astype(float), tested['expected'].astype(float), tested['k'].astype(int)
        assert (e > 0).all() and (e <= n * (1 + 1e-6)).all()
        assert np.allclose(tested['oe'].astype(float), k / e, rtol=1e-4)
        assert np.allclose(tested['expected_unweighted'].astype(float), n / tested['residues_of_type'].astype(float), rtol=1e-4)
        q, p = tested['q_upper'].astype(float), tested['p_upper'].astype(float)
        assert ((q >= p * (1 - 1e-3)) & (q <= 1)).all()
        assert (t.loc[t['relative_mutability'] == '', 'p_upper'] == '').all()     # no weights, no test
        assert set(t['oe_status']) == {'descriptive' if src == 'clinvar' else 'test'}
        if src in ('tcga_mc3', 'depmap_24q4'):
            assert (t['k_no_hypermutated'].astype(int) <= t['k'].astype(int)).all()
            assert len(set(t['samples_profiled'])) == 1 and int(t['samples_profiled'].iloc[0]) > int(t['samples_profiled_no_hypermutated'].iloc[0])


def test_frequency_by_group_adds_up_to_k(tables):
    f = read(REL / 'supporting' / 'site_frequency_by_group_v1.7.0.tsv.gz')
    burden = read(REL / 'supporting' / 'sample_burden_v1.7.0.tsv')
    assert md5(REL / 'supporting' / 'sample_burden_v1.7.0.tsv') == md5(LAYERS / 'propensity_sample_burden.tsv')
    for src in ('tcga_mc3', 'depmap_24q4'):
        s = f[f['source'] == src]
        total = s.groupby(['canonical_UniProtAC', 'position'])['samples_mutated'].apply(lambda x: x.astype(int).sum())
        t = tables[src]
        t = t[t['k'] != '0']
        assert dict(zip(zip(t['canonical_UniProtAC'], t['position']), t['k'].astype(int))) == total.to_dict(), src
        den = burden[burden['source'] == src].groupby('group').size()
        assert (s['samples_profiled'].astype(int).to_numpy() == s['group'].map(den).to_numpy()).all()
        assert (s['samples_mutated'].astype(int) <= s['samples_profiled'].astype(int)).all()
    groups = read(REL / 'supporting' / 'group_denominators_v1.7.0.tsv')
    for src in ('tcga_mc3', 'depmap_24q4'):
        g = groups[groups['source'] == src]
        bb = burden[burden['source'] == src]
        assert g['samples_profiled'].astype(int).sum() == len(bb)
        assert g['samples_profiled_no_hypermutated'].astype(int).sum() == int((bb['hypermutated'] != 'True').sum())
        s = f[f['source'] == src].drop_duplicates('group')
        both = s.merge(g, on='group', suffixes=('', '_den'))
        assert len(both) == len(s) and (both['samples_profiled'] == both['samples_profiled_den']).all()
        assert (both['samples_profiled_no_hypermutated'] == both['samples_profiled_no_hypermutated_den']).all()
    projects = read(REL / 'supporting' / 'tcga_mc3_projects_v1.7.0.tsv').set_index('project_id')['tumour_samples']
    tg = groups[groups['source'] == 'tcga_mc3'].set_index('group')['samples_profiled']
    assert len(tg) == 33 and (tg == projects.reindex(tg.index)).all()


def test_site_set_rows_add_up(tables):
    sets = read(REL / 'supporting' / 'site_set_propensity_v1.7.0.tsv')
    assert len(sets) == 36 and set(sets['iterations']) == {'10000'}
    for _, r in sets.iterrows():
        t = tables[r['source']]
        sfx = '' if r['cohort'] == 'all' else '_no_hypermutated'
        sub = t[(t['residue'] == r['residue']) & (t['p_upper' + sfx] != '')]
        on = sub[sub['site_sets'].str.split(';').map(lambda s: r['site_set'] in s)]
        assert int(r['observed']) == on['k' + sfx].astype(int).sum(), (r['site_set'], r['source'], r['cohort'])
        assert int(r['sites']) == len(on) and int(r['proteins']) == on['canonical_UniProtAC'].nunique()
        assert np.isclose(float(r['expected']), on['expected' + sfx].astype(float).sum(), rtol=1e-4)
        assert 0 < float(r['p_upper']) <= 1 and 0 < float(r['p_lower']) <= 1


def test_controls():
    c = read(REL / 'supporting' / 'propensity_controls_v1.7.0.tsv')
    hot = c[(c['kind'] == 'recurrent hotspot') & (c['source'] == 'tcga_mc3')]
    assert len(hot) == 26 and (hot['result'] == 'pass').all()          # 13 hotspots, all samples and without the hypermutated
    onco = c[(c['family'] == 'every residue of the entry')]
    assert set(onco['canonical_UniProtAC']) == {'P68431', 'Q71DI3', 'P84243'}
    k36 = onco[(onco['control'].str.startswith('H3K36M')) & (onco['source'] == 'tcga_mc3') & (onco['cohort'] == 'all')]
    assert set(k36.loc[k36['result'] == 'pass', 'canonical_UniProtAC']) == {'P68431', 'P84243'}
    k27 = onco[(onco['control'].str.startswith('H3K27M')) & (onco['source'] == 'tcga_mc3') & (onco['cohort'] == 'all')]
    assert (k27['result'] == 'fail').all() and all(r for r in k27['reason'])
    assert set(c['result']) <= {'pass', 'fail', 'not_tested'} and (c['reason'] != '').all()
    g35 = onco[onco['control'].str.startswith('H3G34')]
    assert set(g35['residue']) == {'G'} and set(g35['position']) == {'35'}


def test_validation_agrees():
    v = read(REL / 'supporting' / 'propensity_validation_v1.7.0.tsv')
    assert list(v['source']) == list(SOURCES)
    for _, r in v.iterrows():
        assert int(r['units']) == 1000 and float(r['agreement_p05']) > 0.99 and float(r['max_abs_p_difference']) < 0.05


def test_rates_table_and_dictionary():
    rates = read(REL / 'supporting' / 'mutation_rates_v1.7.0.tsv')
    assert len(rates) == 96 * 4 + 104 and set(rates['source']) == {'tcga_mc3', 'depmap_24q4', 'gnomad_germline'}
    dd = read(REL / 'data_dictionary.tsv')
    assert not dd.duplicated(['file', 'column']).any()
    for f in RELEASE['propensity']['files']:
        opener = gzip.open if f.endswith('.gz') else open
        with opener(REL / f, 'rt', encoding='utf-8') as fh:
            cols = fh.readline().rstrip('\n').split('\t')
        assert list(dd.loc[dd['file'] == f, 'column']) == cols, f
    assert 'descriptive' in ' '.join(dd.loc[dd['file'] == 'supporting/site_propensity_clinvar_v1.7.0.tsv.gz', 'description'])
    assert RELEASE['version'] == 'v1.7.0' and RELEASE['previous_version'] == 'v1.6.0'
    assert 'mutation_rate_methylation_bins' in RELEASE['propensity']['germline_rates']
    readme = (REL / 'README.md').read_text(encoding='utf-8')
    assert 'mutation_rate_methylation_bins' in readme and '## New: mutation propensity' in readme and '\n## New: somatic' not in readme
    assert '\n## Somatic changes at every candidate residue' in readme

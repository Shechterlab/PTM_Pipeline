"""Release v1.8.0: v1.7.0 carried forward, per-site annotations (DepMap cell lines, AlphaMissense) for every site set, ClinVar
disease families and the scope of the site-set test (scripts/73_release_site_annotations.py, scripts/72_stage_alphamissense.py)."""
import gzip
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
REL = ROOT / 'releases' / 'v1.8.0'
PREV = ROOT / 'releases' / 'v1.7.0'
LAYERS = ROOT / 'inputs' / 'layers'
RESIDUES = ['R', 'K', 'S', 'T', 'Y']
pytestmark = pytest.mark.skipif(not (REL / 'release.json').exists(), reason='releases/v1.8.0 not built')
RELEASE = json.loads((REL / 'release.json').read_text(encoding='utf-8')) if (REL / 'release.json').exists() else {}
sys.path.insert(0, str(ROOT / 'scripts'))


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


def test_previous_release_is_carried_and_appended_tables_keep_their_columns():
    edited = set(RELEASE['carried_over_text_edited'])
    appended = set(RELEASE['appended_columns'])
    assert appended == {f'supporting/clinvar_changes_{r}_v1.8.0.tsv.gz' for r in RESIDUES} | {'supporting/site_set_propensity_v1.8.0.tsv'}
    n = 0
    for path in PREV.rglob('*'):
        rel = path.relative_to(PREV).as_posix()
        if not path.is_file() or ('/' not in rel and rel in {'README.md', 'CHANGES.md', 'MD5SUMS', 'data_dictionary.tsv', 'release.json'}) \
                or path.name.startswith('site_changes_'):
            continue
        new_rel = rel.replace('_v1.7.0', '_v1.8.0')
        assert (REL / new_rel).exists(), rel
        if new_rel in appended:
            a, b = read(path), read(REL / new_rel)
            extra = RELEASE['appended_columns'][new_rel]['appended_columns']
            assert list(b.columns) == list(a.columns) + extra and b.drop(columns=extra).equals(a), rel
        elif new_rel not in edited:
            assert md5(REL / new_rel) == md5(path), rel
            n += 1
    assert n > 150
    # the site tables keep their bytes: methylarginine since v1.4.3, the context marks since v1.5.0
    assert md5(REL / 'human_methylarginine_sites_v1.8.0.tsv') == md5(ROOT / 'releases' / 'v1.4.3' / 'human_methylarginine_sites_v1.4.3.tsv')
    for t in ('methyllysine', 'acetylation', 'lysine_acylation', 'ubiquitination', 'sumoylation', 'phosphorylation'):
        assert md5(REL / f'human_{t}_sites_v1.8.0.tsv') == md5(ROOT / 'releases' / 'v1.5.0' / f'human_{t}_sites_v1.5.0.tsv'), t


def test_site_annotations():
    config = json.loads((ROOT / 'config' / 'site_sets.json').read_text(encoding='utf-8'))
    am = read(LAYERS / 'alphamissense_residues.tsv.gz')
    am_at = {(a, p, r): (mn, mx, n) for a, p, r, mn, mx, n in zip(am['canonical_UniProtAC'], am['position'], am['residue'],
                                                                 am['am_mean_pathogenicity'], am['am_max_pathogenicity'], am['am_pathogenic'])}
    models = read(REL / 'supporting' / 'depmap_24q4_change_models_v1.8.0.tsv').set_index('model_id')['cell_line']
    dm = {r: read(REL / 'supporting' / f'depmap_24q4_changes_{r}_v1.8.0.tsv.gz') for r in RESIDUES}
    for name, spec in config['site_sets'].items():
        t = read(REL / 'supporting' / f'site_annotations_{name}_v1.8.0.tsv')
        assert not t.duplicated(['canonical_UniProtAC', 'position', 'residue']).any() and set(t['residue']) <= set(spec['residues'])
        assert len(t) == RELEASE['site_annotations']['per_set'][name]['sites']
        for r in spec['residues']:
            # every site of the set at this residue is a row (the per-site change tables name the sites they hold)
            ch = dm[r]
            flagged = ch[ch['site_sets'].str.split(';').map(lambda v: name in v)]
            want = flagged.groupby(['canonical_UniProtAC', 'position'])['model_id'].apply(lambda s: sorted(set(s)))
            got = t[t['residue'] == r].set_index(['canonical_UniProtAC', 'position'])
            for key, ms in want.items():
                row = got.loc[key]
                assert sorted(row['depmap_models'].split('|')) == ms and int(row['depmap_lines']) == len(ms), (name, key)
                assert row['depmap_cell_lines'].split('|') == [models[m] for m in row['depmap_models'].split('|')], (name, key)
            assert int((got['depmap_lines'].astype(int) > 0).sum()) == len(want), (name, r)
        sample = t.sample(n=min(500, len(t)), random_state=1)
        for a, p, r, mn, mx, n in zip(sample['canonical_UniProtAC'], sample['position'], sample['residue'], sample['am_mean_pathogenicity'],
                                      sample['am_max_pathogenicity'], sample['am_pathogenic']):
            assert (mn, mx, n) == am_at.get((a, p, r), ('', '', '')), (name, a, p)
    # the methylarginine counts agree with the site table's own DepMap columns (the arginine layer)
    s = read(REL / 'human_methylarginine_sites_v1.8.0.tsv', usecols=['uniprot_accession', 'position', 'depmap_change_lines'])
    a = read(REL / 'supporting' / 'site_annotations_methylarginine_v1.8.0.tsv')
    m = s.merge(a, left_on=['uniprot_accession', 'position'], right_on=['canonical_UniProtAC', 'position'])
    single = m[~m['depmap_change_lines'].str.contains(';')]
    assert len(m) == 11632 and (single['depmap_change_lines'].replace('', '0').astype(int) == single['depmap_lines'].astype(int)).all()


def test_clinvar_disease_families():
    from disease_families import classify_clinvar
    drop = {'not specified', 'not provided', 'see cases'}
    t = read(REL / 'supporting' / 'clinvar_changes_R_v1.8.0.tsv.gz', usecols=['conditions', 'germline_classification', 'disease_families'])
    for c, f in t.drop_duplicates('conditions').sample(n=2000, random_state=2)[['conditions', 'disease_families']].itertuples(index=False):
        names = {x.strip() for x in re.split(r'[|;]', c) if x.strip() and x.strip().lower() not in drop}
        want = ';'.join(sorted({classify_clinvar(n)[0] for n in names} - {'Other / unassigned'}))
        assert f == want, c
    assert 'Other / unassigned' not in set(';'.join(t['disease_families']).split(';'))
    lynch = t[t['conditions'].str.fullmatch('Lynch syndrome')]
    assert len(lynch) and set(lynch['disease_families']) == {'Cancer / neoplasia'}


def test_site_set_scope_columns():
    sets = read(REL / 'supporting' / 'site_set_propensity_v1.8.0.tsv')
    prev = read(PREV / 'supporting' / 'site_set_propensity_v1.7.0.tsv')
    assert sets.drop(columns=['sites_total', 'sites_weighted', 'k_total']).equals(prev)
    s = sets.set_index(['site_set', 'residue', 'source', 'cohort'])
    assert s.loc[('methylarginine', 'R', 'tcga_mc3', 'all'), ['sites_total', 'sites_weighted', 'k_total', 'observed']].tolist() == ['11632', '11235', '4516', '4462']
    assert s.loc[('methylarginine', 'R', 'gnomad_v4_1', 'all'), 'k_total'] == '20550'
    st, sw, sites = sets['sites_total'].astype(int), sets['sites_weighted'].astype(int), sets['sites'].astype(int)
    assert ((sites <= sw) & (sw <= st)).all() and (sets['observed'].astype(int) <= sets['k_total'].astype(int)).all()


def test_dictionary_title_and_release_json():
    dd = read(REL / 'data_dictionary.tsv')
    assert not dd.duplicated(['file', 'column']).any()
    files = list(RELEASE['site_annotations']['files']) + list(RELEASE['appended_columns'])
    for f in files:
        opener = gzip.open if f.endswith('.gz') else open
        with opener(REL / f, 'rt', encoding='utf-8') as fh:
            cols = fh.readline().rstrip('\n').split('\t')
        assert list(dd.loc[dd['file'] == f, 'column']) == cols, f
    assert RELEASE['version'] == 'v1.8.0' and RELEASE['previous_version'] == 'v1.7.0'
    assert all(w in RELEASE['title'] for w in ('arginine methylation', 'lysine methylation', 'acylations', 'ubiquitination', 'SUMOylation', 'phosphorylation'))
    readme = (REL / 'README.md').read_text(encoding='utf-8')
    assert readme.startswith('# Human post-translational modification sites') and 'intrinsically disordered' in readme
    assert RELEASE['propensity']['files'] and RELEASE['site_annotations']['alphamissense']['counts']['entries_matched'] == 19382

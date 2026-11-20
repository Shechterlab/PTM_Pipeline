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
    for src in ('tcga_mc3', 'depmap_24q4', 'clinvar', 'gnomad_v4_1'):
        staged = read(LAYERS / f'{src}_changes_{res}.tsv.gz')
        t = read(REL / 'supporting' / f'{src}_changes_{res}_v1.6.0.tsv.gz')
        assert list(t.columns) == list(staged.columns) + ['site_sets'] and list(t.columns)[-2:] == ['gene', 'site_sets']
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
    for r in RELEASE['germline_variants']['per_site_set']:
        assert r['clinvar_pathogenic_sites'] <= r['clinvar_sites_changed'] <= r['sites'] and r['gnomad_v4_1_sites_changed'] <= r['sites']
    for key, s in RELEASE['germline_variants']['per_residue'].items():
        src, res = key.split('/')
        assert s['rows'] == len(read(REL / 'supporting' / f'{src}_changes_{res}_v1.6.0.tsv.gz'))
    dd = read(REL / 'data_dictionary.tsv')
    for f in list(RELEASE['cancer_changes']['files']) + list(RELEASE['germline_variants']['files']):
        path = REL / f
        opener = gzip.open if f.endswith('.gz') else open
        with opener(path, 'rt', encoding='utf-8') as fh:
            cols = fh.readline().rstrip('\n').split('\t')
        assert list(dd.loc[dd['file'] == f, 'column']) == cols, f
    assert not dd.duplicated(['file', 'column']).any()
    assert RELEASE['version'] == 'v1.6.0' and RELEASE['previous_version'] == 'v1.5.1'


def test_rbdmap_pclap_line_corrected_and_previous_release_untouched():
    readme = (REL / 'README.md').read_text(encoding='utf-8')
    assert 'their extents are not republished' not in readme
    assert 'The methylarginine tables carry them, with their extents and sequences' in readme and 'PXD005739' in readme
    assert '## Corrected: RBDmap and pCLAP in the interaction notes' in (REL / 'CHANGES.md').read_text(encoding='utf-8')
    assert RELEASE['corrections']
    # the published v1.5.1 folder is left as archived
    assert 'their extents are not republished' in (PREV / 'README.md').read_text(encoding='utf-8')
    assert RELEASE['germline_variants']['gnomad']['license'].startswith('CC0')


def test_changes_not_shown():
    """supporting/variant_changes_excluded: the changes a source places on an entry whose residue there differs from the canonical
    sequence's (or lies beyond it), where either residue is a candidate residue; none of them is in a change table."""
    import sys
    sys.path.insert(0, str(ROOT / 'scripts'))
    from common import parse_fasta  # noqa: E402
    seqs = parse_fasta(ROOT / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03.fasta')
    ex = read(REL / 'supporting' / 'variant_changes_excluded_v1.6.0.tsv.gz')
    assert list(ex.columns) == ['canonical_UniProtAC', 'position', 'source', 'source_id', 'gene', 'source_ref', 'canonical_ref', 'alt', 'reason']
    assert set(ex['source']) == {'tcga_mc3', 'depmap_24q4', 'clinvar', 'gnomad_v4_1'} and set(ex['reason']) <= {'residue_differs', 'position_beyond_sequence'}
    assert (ex['source_ref'] != ex['canonical_ref']).all()
    assert (ex['source_ref'].isin(list('RKSTY')) | ex['canonical_ref'].isin(list('RKSTY'))).all()
    for a, p, c, why in zip(ex['canonical_UniProtAC'], ex['position'].astype(int), ex['canonical_ref'], ex['reason']):
        seq = seqs[a]
        assert (c == seq[p - 1]) if why == 'residue_differs' else (c == '' and p > len(seq))
    for src in ('tcga_mc3', 'depmap_24q4', 'clinvar', 'gnomad_v4_1'):
        who = {'tcga_mc3': 'sample_id', 'depmap_24q4': 'model_id', 'clinvar': 'variation_id', 'gnomad_v4_1': 'variant_id'}[src]
        kept = set()
        for res in RESIDUES:
            t = read(REL / 'supporting' / f'{src}_changes_{res}_v1.6.0.tsv.gz')
            kept |= set(zip(t['canonical_UniProtAC'], t['position'], t['alt'], t[who]))
        e = ex[ex['source'] == src]
        assert not kept & set(zip(e['canonical_UniProtAC'], e['position'], e['alt'], e['source_id'])), src
        if src in ('clinvar', 'gnomad_v4_1'):        # a variant the entry kept under another expression or transcript is not listed
            assert not {(a, w) for a, _, _, w in kept} & set(zip(e['canonical_UniProtAC'], e['source_id'])), src
    dd = read(REL / 'data_dictionary.tsv')
    assert list(dd.loc[dd['file'] == 'supporting/variant_changes_excluded_v1.6.0.tsv.gz', 'column']) == list(ex.columns)
    assert RELEASE['variant_changes_excluded']['rows'] == len(ex)

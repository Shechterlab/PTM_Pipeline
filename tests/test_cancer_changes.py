"""The changes at every candidate residue in TCGA tumors and DepMap cell lines (K, S, T, Y; scripts/66_stage_cancer_changes.py)
and in ClinVar (R, K, S, T, Y; scripts/67_stage_clinvar.py), and the per-protein backgrounds: on the residue they name, one row
per change and sample or variation, counted in the backgrounds as the tables count them, and as inputs/MANIFEST.tsv records them."""
import csv
import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
LAYERS = ROOT / 'inputs' / 'layers'
RESIDUES = ['K', 'S', 'T', 'Y']
pytestmark = pytest.mark.skipif(not (LAYERS / 'cancer_changes_summary.json').exists(), reason='cancer changes not staged')


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def fasta(path: Path) -> dict:
    seqs, acc = {}, ''
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            acc = line[1:].split('|')[1]
            seqs[acc] = []
        elif acc:
            seqs[acc].append(line.strip())
    return {a: ''.join(s) for a, s in seqs.items()}


CANONICAL = fasta(ROOT / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03.fasta')


def read(name: str) -> pd.DataFrame:
    return pd.read_csv(LAYERS / name, sep='\t', dtype=str, keep_default_na=False)


@pytest.mark.parametrize('res', RESIDUES)
def test_change_tables_are_on_their_residue_in_the_arginine_columns(res):
    for name, arginine, sample in ((f'tcga_mc3_changes_{res}.tsv.gz', 'tcga_mc3_arginine_changes.tsv', 'sample_id'),
                                   (f'depmap_24q4_changes_{res}.tsv.gz', 'depmap_24q4_arginine_changes.tsv', 'model_id')):
        t = read(name)
        assert list(t.columns) == list(pd.read_csv(LAYERS / arginine, sep='\t', nrows=0).columns), name
        assert len(t) > 10000, name
        assert all(CANONICAL[a][int(p) - 1] == res for a, p in zip(t['canonical_UniProtAC'], t['position'])), name
        assert not t.duplicated(['canonical_UniProtAC', 'position', 'alt', sample]).any(), name
        assert set(t['alt']) <= set('ACDEFGHIKLMNPQRSTVWY*') - {res}, name
    d = read(f'depmap_24q4_changes_{res}.tsv.gz')
    assert (d['protein_change'] == res + d['position'] + d['alt']).all()
    t = read(f'tcga_mc3_changes_{res}.tsv.gz')
    assert set(t['filter']) <= {'PASS', 'wga'} and set(t.loc[t['filter'] == 'wga', 'project_id']) <= {'TCGA-OV', 'TCGA-LAML'}


def test_background_counts_what_the_change_tables_hold():
    bg = read('cancer_change_protein_background.tsv.gz')
    assert set(bg['residue']) == {'R', 'K', 'S', 'T', 'Y', 'all'} and len(bg) == 6 * len(CANONICAL)
    bg = bg.set_index(['canonical_UniProtAC', 'residue']).astype(int)
    seq_counts = {(a, r): (len(s) if r == 'all' else s.count(r)) for a, s in CANONICAL.items() for r in ('R', 'K', 'S', 'T', 'Y', 'all')}
    assert all(bg.loc[k, 'residues'] == v for k, v in seq_counts.items())
    tables = {r: (f'tcga_mc3_changes_{r}.tsv.gz', f'depmap_24q4_changes_{r}.tsv.gz') for r in RESIDUES}
    tables['R'] = ('tcga_mc3_arginine_changes.tsv', 'depmap_24q4_arginine_changes.tsv')     # the arginine layer, by the same rules
    for res, (tname, dname) in tables.items():
        for name, who, cols in ((tname, 'sample_id', ['tcga_changes', 'tcga_samples', 'tcga_positions']),
                                (dname, 'model_id', ['depmap_changes', 'depmap_lines', 'depmap_positions'])):
            t = read(name)
            g = t.groupby('canonical_UniProtAC').agg(n=(who, 'size'), w=(who, 'nunique'), p=('position', 'nunique'))
            sub = bg.xs(res, level='residue')
            assert (sub.loc[g.index, cols].to_numpy() == g[['n', 'w', 'p']].to_numpy()).all(), (res, name)
            assert sub.drop(g.index)[cols[0]].eq(0).all(), (res, name)
    classes = bg.drop(index='all', level='residue').groupby(level=0).sum()
    total = bg.xs('all', level='residue')
    for c in ('tcga_changes', 'depmap_changes', 'residues'):
        assert (total[c] >= classes.loc[total.index, c]).all(), c


def test_summary_and_manifest():
    summary = json.loads((LAYERS / 'cancer_changes_summary.json').read_text(encoding='utf-8'))
    assert summary['inputs']['maf']['md5'] == '639ad8f8386e98dacc22e439188aa8fa'                 # the MAF of the arginine layer
    assert summary['inputs']['depmap_mutations']['md5'] == '7bdba347a1602fe96d5654a74d6e52f1'    # DepMap 24Q4, as scripts/44
    assert summary['denominators']['tcga_tumor_samples'] == 10294 and summary['denominators']['tcga_projects'] == 33
    manifest = {r['file']: r for r in csv.DictReader(open(ROOT / 'inputs' / 'MANIFEST.tsv', encoding='utf-8'), delimiter='\t')}
    names = [f'{p}_{r}.tsv.gz' for r in RESIDUES for p in ('tcga_mc3_changes', 'depmap_24q4_changes')]
    for name in names + ['cancer_change_protein_background.tsv.gz', 'cancer_changes_summary.json', 'raw_files.tsv',
                         'depmap_24q4_change_models.tsv']:
        row = manifest[f'layers/{name}']
        assert row['redistribution'] == 'git' and row['md5'] == md5(LAYERS / name), name
    for res in RESIDUES:
        n = summary['residues'][res]
        assert n['tcga_changes'] == len(read(f'tcga_mc3_changes_{res}.tsv.gz')) == int(manifest[f'layers/tcga_mc3_changes_{res}.tsv.gz']['rows'])
        assert n['md5'][f'tcga_mc3_changes_{res}.tsv.gz'] == md5(LAYERS / f'tcga_mc3_changes_{res}.tsv.gz')
    models = set(read('depmap_24q4_change_models.tsv')['model_id'])
    for res in RESIDUES:
        assert set(read(f'depmap_24q4_changes_{res}.tsv.gz')['model_id']) <= models, res
    assert set(read('depmap_24q4_arginine_changes.tsv')['model_id']) <= models
    raw = {r['file']: r for r in csv.DictReader(open(LAYERS / 'raw_files.tsv', encoding='utf-8'), delimiter='\t')}
    assert raw['uniprot_2026_04_human_reviewed_ensembl.tsv.gz']['md5'] == summary['inputs']['uniprot']['ensembl']['md5']
    assert raw['uniprot_human_reviewed_ensembl.tsv.gz']['version'].startswith('UniProt 2026_03')    # the arginine layer's row is kept


# ClinVar (scripts/67_stage_clinvar.py)
CLINVAR = ['R', 'K', 'S', 'T', 'Y']


@pytest.mark.skipif(not (LAYERS / 'clinvar_summary.json').exists(), reason='ClinVar not staged')
@pytest.mark.parametrize('res', CLINVAR)
def test_clinvar_tables_are_on_their_residue(res):
    t = read(f'clinvar_changes_{res}.tsv.gz')
    assert len(t) > 10000
    assert all(CANONICAL[a][int(p) - 1] == res for a, p in zip(t['canonical_UniProtAC'], t['position']))
    assert not t.duplicated(['canonical_UniProtAC', 'position', 'alt', 'variation_id']).any()
    assert set(t['alt']) <= set('ACDEFGHIKLMNPQRSTVWY*') - {res}
    assert t['refseq_protein'].str.match(r'^NP_\d+\.\d+$').all() and set(t['review_stars']) <= set('01234')
    assert not {'rcv', 'allele_id'} & set(t.columns)


@pytest.mark.skipif(not (LAYERS / 'clinvar_summary.json').exists(), reason='ClinVar not staged')
def test_clinvar_known_variants_and_background():
    t = read('clinvar_changes_R.tsv.gz').set_index(['canonical_UniProtAC', 'position', 'alt', 'variation_id'])
    r175h = t.loc[('P04637', '175', 'H', '12374')]                      # TP53 R175H
    assert r175h['germline_classification'] == 'Pathogenic' and r175h['oncogenicity'] == 'Oncogenic'
    bg = read('clinvar_protein_background.tsv.gz').set_index(['canonical_UniProtAC', 'residue']).astype(int)
    assert len(bg) == 6 * len(CANONICAL)
    summary = json.loads((LAYERS / 'clinvar_summary.json').read_text(encoding='utf-8'))
    manifest = {r['file']: r for r in csv.DictReader(open(ROOT / 'inputs' / 'MANIFEST.tsv', encoding='utf-8'), delimiter='\t')}
    for res in CLINVAR:
        x = read(f'clinvar_changes_{res}.tsv.gz')
        g = x.groupby('canonical_UniProtAC').agg(n=('variation_id', 'size'), p=('position', 'nunique'))
        sub = bg.xs(res, level='residue')
        assert (sub.loc[g.index, ['clinvar_variants', 'clinvar_positions']].to_numpy() == g[['n', 'p']].to_numpy()).all(), res
        assert summary['residues'][res]['variants'] == len(x) == int(manifest[f'layers/clinvar_changes_{res}.tsv.gz']['rows'])
        assert manifest[f'layers/clinvar_changes_{res}.tsv.gz']['md5'] == md5(LAYERS / f'clinvar_changes_{res}.tsv.gz')
    for name in ('clinvar_protein_background.tsv.gz', 'clinvar_summary.json'):
        assert manifest[f'layers/{name}']['md5'] == md5(LAYERS / name)
    assert (bg['clinvar_pathogenic'] <= bg['clinvar_variants']).all()

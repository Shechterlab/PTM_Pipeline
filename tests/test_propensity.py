"""The inputs of the per-site mutation propensity (scripts/70_stage_propensity.py) and its method (scripts/propensity.py)."""
import csv
import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
LAYERS = ROOT / 'inputs' / 'layers'
sys.path.insert(0, str(ROOT / 'scripts'))
import propensity as P  # noqa: E402
from common import parse_fasta  # noqa: E402

staged = pytest.mark.skipif(not (LAYERS / 'propensity_inputs_summary.json').exists(), reason='propensity inputs not staged')


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def read(name: str) -> pd.DataFrame:
    return pd.read_csv(LAYERS / name, sep='\t', dtype=str, keep_default_na=False)


@pytest.fixture(scope='module')
def cds():
    return read('propensity_cds.tsv.gz')


@pytest.fixture(scope='module')
def rates():
    return P.somatic_rates(read('propensity_somatic_spectra.tsv')), P.germline_rates(LAYERS / 'gnomad_v2_1_1_mutation_rates.txt.bgz')


def translate(s: str) -> str:
    return ''.join(P.CODE[s[i:i + 3]] for i in range(0, len(s), 3))


@staged
def test_coding_sequences_translate_to_the_canonical_sequences(cds):
    seqs = parse_fasta(ROOT / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03.fasta')
    for acc, s in zip(cds['canonical_UniProtAC'], cds['cds']):
        assert translate(s) == seqs[acc], acc
    status = read('propensity_cds_status.tsv')
    assert set(status['canonical_UniProtAC']) == set(seqs) and not status['canonical_UniProtAC'].duplicated().any()
    assert set(status.loc[status['status'] == 'matched', 'canonical_UniProtAC']) == set(cds['canonical_UniProtAC'])
    for matched, acc in zip(status['transcripts_matched'], status['canonical_UniProtAC']):
        if matched:
            rows = cds[cds['canonical_UniProtAC'] == acc]
            named = set(rows['transcript_id']) | {t for x in rows['same_cds_transcripts'] for t in x.split(';') if t}
            assert named == set(matched.split(';')), acc


def brute_force(row, key, somatic, germline):
    """w per residue, one substitution at a time, from the strings."""
    c = row['cds']
    left, right = [''] + list(c[:-1]), list(c[1:]) + ['']
    left[0], right[-1] = row['flank_5'], row['flank_3']
    for item in filter(None, row['junctions'].split(',')):
        j, xy = item.split(':')
        right[int(j) - 1], left[int(j)] = xy[0], xy[1]
    level = {int(a): int(b) for a, b in (x.split(':') for x in filter(None, row['cpg_methylation'].split(',')))}
    out = []
    for ci in range(len(c) // 3):
        cod, total = c[3 * ci:3 * ci + 3], 0.0
        for o in range(3):
            i = 3 * ci + o
            for a in 'ACGT':
                new = cod[:o] + a + cod[o + 1:]
                if a == cod[o] or P.CODE[new] == P.CODE[cod]:
                    continue
                ctx, ref, alt = P.collapse(left[i] + c[i] + right[i], a)
                k = P.CLASS_INDEX[(ctx, ref, alt)]
                if key == 'germline':
                    total += germline[k, level.get(i, 0) if ctx[1:] == 'CG' and alt == 'T' else 0]
                else:
                    total += somatic[key][k]
        out.append(total)
    return np.array(out)


@staged
def test_mutability_matches_a_substitution_by_substitution_count(cds, rates):
    somatic, germline = rates
    rng = np.random.default_rng(7)
    rows = list(rng.choice(len(cds), 25, replace=False)) + list(np.nonzero((cds['canonical_UniProtAC'] == 'P04637').to_numpy())[0])
    for i in rows:
        row = cds.iloc[i]
        got = P.codon_mutability(row, somatic, germline)
        for key in ('germline', 'tcga_mc3/all', 'tcga_mc3/no_hypermutated', 'depmap_24q4/all'):
            want = brute_force(row, key, somatic, germline)
            assert np.allclose(got[key], want, rtol=1e-12, atol=0), (row['transcript_id'], key)
    # TP53 R175 (CGC): its first C is a CpG transition in a methylated context, so it is more mutable than TP53's average arginine
    tp = cds[cds['canonical_UniProtAC'] == 'P04637'].iloc[0]
    w = P.codon_mutability(tp, somatic, germline)['germline']
    aa = translate(tp['cds'])
    assert aa[174] == 'R' and w[174] > np.mean([w[i] for i, a in enumerate(aa) if a == 'R'])


@staged
def test_germline_rate_table_is_gnomads():
    path = LAYERS / 'gnomad_v2_1_1_mutation_rates.txt.bgz'
    assert md5(path) == 'a089b6921f1061985af476aa3459a79e'           # GCS md5Hash oIm2kh8QYZha9HaqNFmnng== of the gnomAD object
    with gzip.open(path, 'rt') as fh:
        t = pd.read_csv(fh, sep='\t')
    assert len(t) == 104 and list(t.columns) == ['context', 'ref', 'alt', 'methylation_level', 'mu_snp']
    hi = t[t['methylation_level'] > 0]
    assert len(hi) == 8 and set(hi['context'].str[1:]) == {'CG'} and set(hi['alt']) == {'T'}
    r = P.germline_rates(path)
    assert r.shape == (96, 3) and (r[P.CPG_TRANSITION, 2] > r[P.CPG_TRANSITION, 0]).all()


@staged
def test_spectra_and_burden():
    summary = json.loads((LAYERS / 'propensity_inputs_summary.json').read_text(encoding='utf-8'))
    sp = read('propensity_somatic_spectra.tsv')
    b = read('propensity_sample_burden.tsv')
    assert len(sp) == 96 * 4 and sp.groupby(['source', 'cohort']).size().eq(96).all()
    assert sp.groupby(['context', 'ref', 'alt'])['opportunity'].nunique().eq(1).all()
    for src in ('tcga_mc3', 'depmap_24q4'):
        s, bb = sp[sp['source'] == src], b[b['source'] == src]
        hyper = bb['hypermutated'] == 'True'
        assert s.loc[s['cohort'] == 'all', 'observed'].astype(int).sum() == bb['snvs_in_coding_sequences'].astype(int).sum()
        assert s.loc[s['cohort'] == 'no_hypermutated', 'observed'].astype(int).sum() == bb.loc[~hyper, 'snvs_in_coding_sequences'].astype(int).sum()
        assert set(s.loc[s['cohort'] == 'no_hypermutated', 'samples'].astype(int)) == {int((~hyper).sum())}
        per_mb = bb['coding_calls'].astype(int) / (summary['coding_bp'] / 1e6)
        assert np.allclose(per_mb.round(4), bb['burden_per_mb'].astype(float))
        assert (hyper == (per_mb > 10)).all() and set(bb['threshold_per_mb']) == {'10.0'}
    t = b[b['source'] == 'tcga_mc3']
    assert len(t) == 10294 and t['group'].nunique() == 33 and not t['sample_id'].duplicated().any()
    assert summary['tcga_hypermutated'] == int((t['hypermutated'] == 'True').sum())
    d = b[b['source'] == 'depmap_24q4']
    # DepMap's mutation file keeps synonymous calls only where they touch a splice region
    assert len(d) == 1929 and d['silent'].astype(int).sum() < 0.001 * d['coding_calls'].astype(int).sum()


@staged
def test_h3_tables_hold_the_oncohistone_changes():
    seqs = parse_fasta(ROOT / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03.fasta')
    t = read('tcga_mc3_changes_h3_all_residues.tsv.gz')
    assert all(seqs[a][int(p) - 1] == r for a, p, r in zip(t['canonical_UniProtAC'], t['position'], t['residue']))
    k28m = t[(t['position'] == '28') & (t['alt'] == 'M')]
    assert set(zip(k28m['canonical_UniProtAC'], k28m['project_id'], k28m['gene'])) == {('P84243', 'TCGA-LGG', 'H3-3A'), ('P68431', 'TCGA-LAML', 'H3C4')}
    assert set(t['residue']) - set('RKSTY')                          # every residue, glycine 35 among them
    assert ((t['position'] == '35') & (t['residue'] == 'G') & t['alt'].isin(['R', 'V'])).sum() >= 5


@staged
def test_manifest_rows():
    manifest = {r['file']: r for r in csv.DictReader(open(ROOT / 'inputs' / 'MANIFEST.tsv', encoding='utf-8'), delimiter='\t')}
    for name in ('propensity_cds.tsv.gz', 'propensity_cds_status.tsv', 'gnomad_v2_1_1_mutation_rates.txt.bgz', 'propensity_somatic_spectra.tsv',
                 'propensity_sample_burden.tsv', 'tcga_mc3_changes_h3_all_residues.tsv.gz', 'depmap_24q4_changes_h3_all_residues.tsv.gz',
                 'propensity_inputs_summary.json'):
        assert manifest[f'layers/{name}']['md5'] == md5(LAYERS / name), name
    raw = {r['file']: r for r in csv.DictReader(open(LAYERS / 'raw_files.tsv', encoding='utf-8'), delimiter='\t')}
    assert raw['Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz']['md5'] == json.loads((LAYERS / 'propensity_inputs_summary.json').read_text())['inputs']['genome']['md5']
    assert raw['GRCh38_to_GRCh37.chain.gz']['md5'] == '9573fe0b0a4715f0dd4fa62d5bd7810d'        # Ensembl's CHECKSUMS
    assert (LAYERS / 'propensity_cds.tsv.gz').read_bytes()[4:8] == b'\x00\x00\x00\x00'         # gzip without a timestamp


def test_site_test_and_bh_on_a_small_case():
    res = pd.DataFrame({'canonical_UniProtAC': ['A'] * 4 + ['B'] * 2 + ['C'] * 2, 'position': [1, 2, 3, 4, 1, 2, 1, 2],
                        'w': [1.0, 1.0, 1.0, 3.0, 2.0, 2.0, np.nan, 1.0], 'k': [0, 1, 0, 5, 0, 0, 3, 0]})
    out = P.site_test(res, 'w', 'k', lower=True)
    assert list(out['n']) == [6, 6, 6, 6, 0, 0, 3, 3]
    assert np.allclose(out['expected'][:4], [1, 1, 1, 3]) and np.allclose(out['relative_mutability'][:4], [2 / 3, 2 / 3, 2 / 3, 2])
    assert out['p_upper'][4:].isna().all()                           # B: no change; C: no weight for every residue
    from scipy.stats import binom
    assert np.isclose(out['p_upper'][3], binom.sf(4, 6, 0.5)) and np.isclose(out['p_lower'][1], binom.cdf(1, 6, 1 / 6))
    p = np.array([0.01, 0.04, 0.03, np.nan, 0.2])
    assert np.allclose(P.bh(p)[[0, 1, 2, 4]], [0.04, 0.16 / 3, 0.16 / 3, 0.2]) and np.isnan(P.bh(p)[3])
    assert P.changes_needed(0, 10, 0.1, 1) == min(k for k in range(50) if binom.sf(k - 1, 10 + k, 0.1) <= 0.05)


def test_set_test_is_reproducible_and_calibrated():
    units = [(10, 0.2, 2), (5, 0.5, 3), (0, 0.3, 0)]
    a = P.set_test(units, np.random.default_rng(1), iterations=2000)
    b = P.set_test(units, np.random.default_rng(1), iterations=2000)
    assert a == b and a['observed'] == 5 and np.isclose(a['expected'], 4.5)
    assert 0.3 < a['p_upper'] < 0.7 and 0.3 < a['p_lower'] < 0.8

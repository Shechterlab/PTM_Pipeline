"""Tests of the Ochoa 2020 phosphosite staging (scripts/58). The staged table is in Git (published, cited)."""
import csv
import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
spec = importlib.util.spec_from_file_location('stage58', ROOT / 'scripts' / '58_stage_ochoa_phosphosites.py')
stage58 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(stage58)

REFERENCE = {'P00001': 'MKSTYKSS', 'P00002': 'MSAT'}
MARKS = ROOT / 'inputs' / 'marks'


def table(*rows):
    return pd.DataFrame([dict(zip(['uniprot', 'position', 'residue', 'functional_score', 'best_PEP', 'best_localization_prob',
                                   'biological_samples', 'spectralcounts'], r)) for r in rows])


def test_toy_classification_and_columns():
    t = table(('P00001', 3, 'S', 0.5, 0.001, 0.99, 2, 7),        # kept
              ('P00001', 4, 'T', 0.25, 0.0, 1.0, 1, 3),          # kept
              ('P00001', 2, 'S', 0.1, 0.0, 1.0, 1, 1),           # residue mismatch (K)
              ('P00002-2', 2, 'S', 0.1, 0.0, 1.0, 1, 1),         # isoform
              ('P99999', 2, 'S', 0.1, 0.0, 1.0, 1, 1),           # not in reference
              ('P00002', 9, 'S', 0.1, 0.0, 1.0, 1, 1),           # past the end
              ('P00002', 2, 'K', 0.1, 0.0, 1.0, 1, 1))           # not S/T/Y
    sites, bad, counts = stage58.stage(t, REFERENCE)
    assert list(sites.columns) == stage58.COLUMNS and list(bad.columns) == stage58.REJECTED_COLUMNS
    assert counts == {'kept': 2, 'residue_mismatch': 1, 'isoform_accession': 1, 'accession_not_in_reference': 1,
                      'bad_position': 1, 'bad_residue': 1}
    assert sites[['canonical_UniProtAC', 'position', 'residue']].values.tolist() == [['P00001', 3, 'S'], ['P00001', 4, 'T']]
    assert not sites.duplicated(['canonical_UniProtAC', 'position']).any()
    assert all(REFERENCE[a][p - 1] == r for a, p, r in zip(sites.canonical_UniProtAC, sites.position, sites.residue))
    assert sites.loc[0, 'functional_score'] == '0.5' and sites.loc[0, 'spectral_counts'] == '7'
    assert len(sites) + len(bad) == len(t)


def test_merge_stops_on_different_sites():
    a = pd.DataFrame({'uniprot': ['P00001'], 'position': [3], 'functional_score': [0.5]})
    b = pd.DataFrame({'uniprot': ['P00001'], 'position': [4], 'residue': ['T'], 'best_PEP': [0], 'best_localization_prob': [1],
                      'biological_samples': [1], 'spectralcounts': [1]})
    with pytest.raises(SystemExit):
        stage58.merge_tables(a, b)


def test_staged_table_matches_reference_and_summary():
    path = MARKS / stage58.OUT_SITES
    reference = stage58.parse_fasta(ROOT / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03.fasta')
    rows = list(csv.DictReader(open(path, encoding='utf-8'), delimiter='\t'))
    assert list(rows[0].keys()) == stage58.COLUMNS
    keys = [(r['canonical_UniProtAC'], r['position']) for r in rows]
    assert len(keys) == len(set(keys))
    for r in rows:
        assert r['residue'] in 'STY' and reference[r['canonical_UniProtAC']][int(r['position']) - 1] == r['residue']
        assert 0 <= float(r['functional_score']) <= 1 and r['source'] == 'Ochoa 2020' and r['evidence_class'] == 'experimental'
        assert r['mark'] == 'phosphorylation' and r['state'] == 'phospho' and 0 <= float(r['best_localization_probability']) <= 1
    summary = json.loads((MARKS / stage58.OUT_SUMMARY).read_text())
    assert summary['validated_sites'] == len(rows) == sum(summary['validated_sites_by_residue'].values())
    assert summary['validated_proteins'] == len({r['canonical_UniProtAC'] for r in rows})
    assert summary['input_sites'] == sum(summary['row_status'].values())
    assert summary['input_sites'] == summary['validated_sites'] + sum(summary['rejected_by_reason'].values())
    with open(MARKS / stage58.OUT_REJECTED, encoding='utf-8') as fh:
        assert sum(1 for _ in fh) - 1 == sum(summary['rejected_by_reason'].values())


def test_summary_and_manifest_in_git():
    summary = json.loads((MARKS / stage58.OUT_SUMMARY).read_text())
    assert summary['redistribution'] == 'git' and summary['relocated_sites'] == 0
    manifest = {r['file']: r for r in csv.DictReader(open(ROOT / 'inputs' / 'MANIFEST.tsv', encoding='utf-8'), delimiter='\t')}
    for name in (stage58.OUT_SITES, stage58.OUT_REJECTED):
        assert manifest['marks/' + name]['redistribution'] == 'git'
        assert 'exclusive licence to Springer Nature' in manifest['marks/' + name]['licence'] and "owner" not in manifest['marks/' + name]['licence']
        assert 'Nat Biotechnol 38:365' in manifest['marks/' + name]['citation']
        assert manifest['marks/' + name]['md5'] == summary['files_md5']['rejected' if 'rejected' in name else 'validated']
    assert manifest['marks/' + stage58.OUT_SUMMARY]['redistribution'] == 'git'
    raw = {r['file']: r for r in csv.DictReader(open(MARKS / 'raw_files.tsv', encoding='utf-8'), delimiter='\t')}
    assert raw[stage58.RAW_SCORES]['sha256'] == summary['raw']['scores']['sha256']
    assert raw[stage58.RAW_FEATURES]['sha256'] == summary['raw']['features']['sha256']

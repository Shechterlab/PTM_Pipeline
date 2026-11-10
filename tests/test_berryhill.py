"""Tests of the staged Berryhill et al. 2023 lysine methylation sites (scripts/57; Sci Rep 13:377, doi:10.1038/s41598-022-27175-x).

Run from the repository root: python -m pytest -q tests/test_berryhill.py
"""
import csv
import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
spec = importlib.util.spec_from_file_location('stage57', ROOT / 'scripts' / '57_stage_berryhill_lys_methylation.py')
stage = importlib.util.module_from_spec(spec)
spec.loader.exec_module(stage)

MARKS = ROOT / 'inputs' / 'marks'
REFERENCE = {'P00001': 'MAKKRSTYKRGGPEPTIDEKR', 'P00002': 'MAKAKAK', 'P00003': 'AKKAKBAKKAK'}


def row(acc, pos, peptide, pip, state='m1', prob=99.0, cell='HEK293T', enrichment='Kme1-A Rep 1'):
    return {'Kme_Site': f'{acc}-K{pos}-{state}', 'Site.Probability': prob, 'Peptide.Sequence': peptide, 'Position.in.Peptide': pip,
            'Modification': state, 'Protein.Accession': acc, 'Position.in.Protein': pos, 'Cell.Type': cell, 'Enrichment': enrichment}


def test_toy_staging_kept_relocated_rejected_and_merged_rows():
    rows = [row('P00001', 3, 'MAKK', 3),                                   # kept (K3)
            row('P00001', 3, 'MAKK', 3, cell='U2OS', prob='NA'),           # same site, other cell line, no probability
            row('P00001', 4, 'KRS', 1, state='m2'),                        # K4 given, peptide KRS starts at 4: kept
            row('P00001', 9, 'STYKR', 4),                                  # kept (S6 T7 Y8 K9)
            row('P00001', 10, 'PEPTIDEK', 8),                              # pos wrong, peptide unique: K20 (relocated)
            row('P00002', 4, 'AK', 2),                                     # position 4 is A and 'AK' occurs 3 times: not unique
            row('P00001', 5, 'ZZZ', 1),                                    # K not in peptide
            row('P00001', 5, 'RSTK', 9),                                   # bad peptide position
            row('P00004', 2, 'AKK', 2),                                    # accession not in reference
            row('P00001-2', 2, 'AKK', 2),                                  # isoform accession
            row('P00001', 2, 'AKK', 2, state='x'),                         # unknown state
            row('P00001', 2, 'WWWK', 4)]                                   # peptide absent
    sites, rejected, counts = stage.stage(rows, REFERENCE)
    assert counts == {'kept': 4, 'relocated': 1, 'peptide_not_unique': 1, 'peptide_position_not_K': 1, 'bad_position': 1,
                      'accession_not_in_reference': 1, 'isoform_accession': 1, 'bad_state': 1, 'peptide_not_in_canonical': 1}
    key = sites.set_index(['canonical_UniProtAC', 'position', 'state'])
    assert len(sites) == len(key) == 4
    merged = key.loc[('P00001', 3, 'me1')]
    assert merged['berryhill_cell_lines'] == 'HEK293T;U2OS' and merged['berryhill_site_probability'] == '99' and merged['relocated'] == 'no'
    moved = key.loc[('P00001', 20, 'me1')]
    assert moved['relocated'] == 'yes' and moved['source_position'] == '10'
    assert list(sites.columns) == stage.COLUMNS and set(sites['evidence_class']) == {'experimental'} and set(sites['source']) == {'Berryhill 2023'}
    assert len(rejected) == 7 and list(rejected.columns) == stage.REJECTED_COLUMNS


def test_probability_is_empty_when_every_row_is_na():
    sites, _, _ = stage.stage([row('P00001', 3, 'MAKK', 3, prob='NA')], REFERENCE)
    assert list(sites['berryhill_site_probability']) == ['']


def read_table(name):
    with open(MARKS / name, encoding='utf-8') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def test_staged_table_columns_residue_and_no_duplicates():
    reference = stage.parse_fasta(ROOT / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03.fasta')
    rows = read_table('berryhill2023_methyllysine_sites.tsv')
    with open(MARKS / 'berryhill2023_methyllysine_sites.tsv', encoding='utf-8') as fh:
        assert next(csv.reader(fh, delimiter='\t')) == stage.COLUMNS
    assert rows
    keys = [(r['canonical_UniProtAC'], r['position'], r['state']) for r in rows]
    assert len(keys) == len(set(keys))
    for r in rows:
        assert reference[r['canonical_UniProtAC']][int(r['position']) - 1] == 'K' and r['residue'] == 'K'
        assert r['mark'] == 'methyllysine' and r['state'] in {'me1', 'me2', 'me3'} and r['source'] == 'Berryhill 2023' and r['evidence_class'] == 'experimental'
        assert (r['relocated'] == 'yes') == bool(r['source_position'])
        assert r['berryhill_site_probability'] == '' or 90 <= float(r['berryhill_site_probability']) <= 100


def test_summary_counts_match_tables_and_manifest():
    summary = json.loads((MARKS / 'berryhill2023_methyllysine_summary.json').read_text(encoding='utf-8'))
    rows, rejected = read_table('berryhill2023_methyllysine_sites.tsv'), read_table('berryhill2023_methyllysine_rejected.tsv')
    assert summary['validated_sites'] == len(rows)
    assert summary['validated_proteins'] == len({r['canonical_UniProtAC'] for r in rows})
    assert summary['validated_residues'] == len({(r['canonical_UniProtAC'], r['position']) for r in rows})
    for state in ('me1', 'me2', 'me3'):
        assert summary['validated_sites_by_state'][state] == sum(r['state'] == state for r in rows)
        assert summary['validated_proteins_by_state'][state] == len({r['canonical_UniProtAC'] for r in rows if r['state'] == state})
    assert summary['relocated_sites'] == sum(r['relocated'] == 'yes' for r in rows)
    assert sum(summary['rejected_by_reason'].values()) == len(rejected)
    assert sum(summary['row_status'].values()) == summary['input_rows'] == 6576
    assert summary['row_status']['kept'] + summary['row_status']['relocated'] + len(rejected) == summary['input_rows']
    assert {r['reason'] for r in rejected} == set(summary['rejected_by_reason'])
    manifest = {r['file']: r for r in csv.DictReader(open(ROOT / 'inputs' / 'MANIFEST.tsv', encoding='utf-8'), delimiter='\t')}
    for name, kind in (('sites', 'validated'), ('rejected', 'rejected')):
        entry = manifest[f'marks/berryhill2023_methyllysine_{name}.tsv']
        assert entry['md5'] == stage.md5_file(MARKS / f'berryhill2023_methyllysine_{name}.tsv') == summary['files_md5'][kind]
        assert 'CC BY 4.0' in entry['licence'] and "owner's decision" in entry['licence'] and 'Sci Rep 13:377' in entry['citation']
    assert manifest['marks/berryhill2023_methyllysine_sites.tsv']['rows'] == str(len(rows))

"""Tests of the interaction records (scripts/interaction_sources.py) and of the PDB staging for site sets (scripts/53).

Run from the repository root: python -m pytest -q tests/test_site_interaction_records.py
"""
import importlib.util
import json
import sys
from pathlib import Path

import pytest

pd = pytest.importorskip('pandas')

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import interaction_sources as isrc  # noqa: E402


def load(name):
    spec = importlib.util.spec_from_file_location('stage53', ROOT / 'scripts' / name)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


stage = load('53_stage_interaction_sources.py')
CONFIG = json.loads((ROOT / 'config' / 'site_sets.json').read_text(encoding='utf-8'))
DATE = '2026-11-01'


# ---- the methylarginine release -------------------------------------------------------------------------------------

def release_records():
    rd = lambda p: pd.read_csv(ROOT / p, sep='\t', dtype=str, keep_default_na=False)  # noqa: E731
    ev = rd('releases/v1.4.3/supporting/site_interactions_v1.4.3.tsv')
    # the release table writes scores as numbers (0.3); the build passed scripts/37's text (0.30)
    ev['score'] = ev['score'].map(lambda v: f'{float(v):.2f}' if v else '')
    sites = rd('releases/v1.4.3/human_methylarginine_sites_v1.4.3.tsv')
    meta = rd('inputs/reference/uniprot_human_reviewed_2026_03_entries.tsv').set_index('Entry')
    genes = {a: str(g).split(';')[0].strip() for a, g in meta['Gene Names (primary)'].items() if str(g).strip()}
    return isrc.records(ev, sites, ROOT / 'inputs', genes)


def test_methylarginine_records_reproduce_the_release(tmp_path):
    ref = ROOT / 'releases' / 'v1.4.3' / 'supporting' / 'site_interaction_records_v1.4.3.tsv'
    if not ref.exists() or not (ROOT / 'inputs' / 'interactions' / f'pdb_entries_{isrc.DATE}.tsv').exists():
        pytest.skip('release v1.4.3 or its staged interaction inputs are not here')
    t = release_records()
    out = tmp_path / 'records.tsv'   # the string columns of build_release.write_table (release_docs.write_plain_tsv)
    with out.open('w', encoding='utf-8', newline='') as fh:
        fh.write('\t'.join(t.columns) + '\n')
        for row in t.astype(str).itertuples(index=False):
            fh.write('\t'.join(row) + '\n')
    assert out.read_bytes() == ref.read_bytes()


def test_parse_site_and_defaults():
    assert isrc.parse_site('P12345:K27') == ('P12345', 'K', 27)
    assert isrc.parse_site('P12345-2:R5') == ('P12345-2', 'R', 5)
    with pytest.raises(SystemExit):
        isrc.parse_site('P12345')
    assert isrc.ARGININE == (isrc.PDB_STATE, 'ARG', 'binds_methylated')
    assert 'binds_modified' not in isrc.CLASSES and 'Kac' not in isrc.STATES
    assert list(isrc.SITE_SET_CLASSES) == ['binds_modified', 'binds_unmodified', 'binds_other_residue', 'binds_state_not_read', 'nucleic_acid', 'self']
    vocab = isrc.site_set_state_vocabulary(CONFIG, ['S', 'T', 'Y'])
    assert {'serine', 'Sph', 'Tph', 'Yph', 'other_residue', 'not_read', 'not_recorded'} <= set(vocab) and all(vocab.values())
    assert isrc.site_states(CONFIG, ['K'])['K'].states == {'LYS': 'lysine', 'ALY': 'Kac', 'MLZ': 'Kme1', 'MLY': 'Kme2', 'M3L': 'Kme3'}


# ---- records of a site set from tiny staged tables -------------------------------------------------------------------

ENTRY = ['pdb_id', 'title', 'method', 'resolution_a', 'release_date', 'pubmed', 'doi', 'first_author', 'authors', 'year', 'journal', 'citation_title', 'ligands']
ENTITY = ['pdb_id', 'entity_id', 'type', 'description', 'uniprot_accessions', 'organisms', 'taxids']
RESIDUE = ['pdb_id', 'entity_id', 'uniprot_accession', 'position', 'residue', 'residue_name', 'residue_parent', 'entity_position', 'context_check', 'other_modified_residues']
EVIDENCE = ['site_id', 'evidence_type', 'localisation', 'partner_kind', 'partner_accession', 'partner_gene', 'source', 'structures', 'score', 'distance_a', 'detail',
            'counts_as_evidence']


def entry(pdb, pubmed='1234', doi='10.1/x'):
    return dict(zip(ENTRY, [pdb, f'title of {pdb}', 'X-RAY DIFFRACTION', '2.00', '2020-01-01', pubmed, doi, 'Smith', '3', '2020', 'Nature', 'a paper', '']))


def entity(pdb, ent, kind, desc, acc='', org='Homo sapiens'):
    return dict(zip(ENTITY, [pdb, ent, kind, desc, acc, org, '9606']))


def residue(pdb, ent, acc, pos, comp, check='matches'):
    return dict(zip(RESIDUE, [pdb, ent, acc, str(pos), comp, '', comp, str(pos), check, '']))


def ev_row(site, kind, partner_kind, acc='', gene='', source='PDBe-KB', structures='', localisation='at_residue', detail=''):
    return dict(zip(EVIDENCE, [site, kind, localisation, partner_kind, acc, gene, source, structures, '', '', detail, 'true']))


def write_tables(tmp_path, entries, entities, residues, names=()):
    folder = tmp_path / 'interactions'
    folder.mkdir()
    for name, rows, cols in (('pdb_entries', entries, ENTRY), ('pdb_entities', entities, ENTITY), ('pdb_site_residues', residues, RESIDUE),
                             ('uniprot_partner_names', names, ['accession', 'gene', 'organism', 'protein_name'])):
        pd.DataFrame(list(rows), columns=cols).to_csv(folder / f'site_set_{name}_{DATE}.tsv', sep='\t', index=False)
    return tmp_path


def by_source(t, pdb):
    return t[t['source_id'] == pdb].iloc[0]


def test_site_set_records_lysine(tmp_path):
    site = 'P11111:K27'
    entries = [entry(p) for p in ('1aaa', '1bbb', '1ccc', '1ddd', '1eee', '1rna')] + [entry('1fff', pubmed='', doi='')]
    entities = []
    for p in ('1aaa', '1bbb', '1ccc', '1ddd', '1fff'):
        entities += [entity(p, '1', 'polypeptide(L)', 'histone', 'P11111'), entity(p, '2', 'polypeptide(L)', 'reader', 'Q22222')]
    entities += [entity('1eee', '1', 'polypeptide(L)', 'histone', 'P11111'), entity('1eee', '2', 'polypeptide(L)', 'mouse reader', 'Q33333', 'Mus musculus'),
                 entity('1rna', '1', 'polypeptide(L)', 'histone', 'P11111'), entity('1rna', '2', 'polyribonucleotide', 'tRNA-Lys'),
                 entity('1rna', '3', 'polyribonucleotide', 'mRNA')]
    residues = [residue('1aaa', '1', 'P11111', 27, 'LYS'), residue('1bbb', '1', 'P11111', 27, 'ALY'), residue('1ccc', '1', 'P11111', 27, 'ARG'),
                residue('1eee', '1', 'P11111', 27, 'LYS'), residue('1rna', '1', 'P11111', 27, 'LYS'), residue('1fff', '1', 'P11111', 27, 'LYS')]
    names = [{'accession': 'Q33333', 'gene': 'Rdr', 'organism': 'Mus musculus', 'protein_name': 'Mouse reader'}]
    inputs = write_tables(tmp_path, entries, entities, residues, names)
    ev = pd.DataFrame([
        ev_row(site, 'pdb_protein', 'protein', 'Q22222', 'RDR2', structures='1aaa; 1bbb; 1ccc; 1ddd; 1fff'),
        ev_row(site, 'pdb_protein', 'protein', 'Q33333', '', structures='1eee'),
        ev_row(site, 'pdb_rna', 'rna', structures='1rna'),
        ev_row(site, 'pdb_self', 'self', structures='1aaa'),
        ev_row(site, 'rna_crosslink', 'rna', source='pepR-MS', localisation='in_peptide'),
        ev_row(site, 'capri_xl_peptide', 'rna', source='CAPRI', localisation='in_peptide', detail='20-30 AKQTARK'),
    ], columns=EVIDENCE)
    t = isrc.site_set_records(ev, inputs, {'Q22222': 'RDR2', 'P11111': 'H3'}, CONFIG, ['K'])
    assert list(t.columns) == [c for c, *_ in isrc.RECORD_COLUMNS]
    pr = t[t['evidence_types'] == 'pdb_protein']
    assert (by_source(pr, '1aaa')['record_class'], by_source(pr, '1aaa')['site_state']) == ('binds_unmodified', 'lysine')
    assert (by_source(pr, '1bbb')['record_class'], by_source(pr, '1bbb')['site_state']) == ('binds_modified', 'Kac')
    assert by_source(pr, '1bbb')['state_detail'].startswith('ALY modelled')
    assert (by_source(pr, '1ccc')['record_class'], by_source(pr, '1ccc')['site_state']) == ('binds_other_residue', 'other_residue')
    assert (by_source(pr, '1ddd')['record_class'], by_source(pr, '1ddd')['site_state']) == ('binds_state_not_read', 'not_read')
    assert by_source(pr, '1aaa')['partner_gene'] == 'RDR2' and by_source(pr, '1aaa')['partner_name'] == 'reader'
    assert by_source(pr, '1aaa')['flags'] == ''
    mouse = by_source(pr, '1eee')
    assert mouse['partner_gene'] == 'Rdr' and mouse['partner_species'] == 'Mus musculus' and mouse['flags'] == 'partner_other_species'
    assert by_source(pr, '1fff')['flags'] == 'no_publication'
    rna = t[t['evidence_types'] == 'pdb_rna'].iloc[0]
    assert (rna['record_class'], rna['partner_kind'], rna['partner_name']) == ('nucleic_acid', 'rna', 'mRNA, tRNA-Lys')
    assert rna['site_state'] == 'lysine'
    me = t[t['evidence_types'] == 'pdb_self'].iloc[0]
    assert (me['record_class'], me['partner_kind'], me['partner_accession']) == ('self', 'self', 'P11111')
    for kind, source in (('rna_crosslink', 'pepR-MS'), ('capri_xl_peptide', 'CAPRI')):
        x = t[t['evidence_types'] == kind].iloc[0]
        cite, doi, method = isrc.DATASETS[source]
        assert (x['record_class'], x['source'], x['method'], x['doi'], x['citation'], x['context']) == ('nucleic_acid', source, method, doi, cite, 'in peptide')
        assert (x['site_state'], x['flags'], x['partner_kind']) == ('not_recorded', '', 'rna')
    # the class order of the table: binds_modified, binds_unmodified, binds_other_residue, binds_state_not_read, nucleic_acid, self
    order = list(isrc.SITE_SET_CLASSES)
    assert [order.index(c) for c in t['record_class']] == sorted(order.index(c) for c in t['record_class'])


def test_site_set_records_phospho(tmp_path):
    entries = [entry('2aaa'), entry('2bbb')]
    entities = [entity(p, e, 'polypeptide(L)', d, a) for p in ('2aaa', '2bbb') for e, d, a in (('1', 'kinase substrate', 'P33333'), ('2', 'adaptor', 'Q44444'))]
    residues = [residue('2aaa', '1', 'P33333', 10, 'SEP'), residue('2bbb', '1', 'P33333', 10, 'THR')]   # a threonine where the site is a serine
    inputs = write_tables(tmp_path, entries, entities, residues)
    site = 'P33333:S10'
    ev = pd.DataFrame([ev_row(site, 'pdb_protein', 'protein', 'Q44444', 'ADP', structures='2aaa; 2bbb')], columns=EVIDENCE)
    t = isrc.site_set_records(ev, inputs, {'Q44444': 'ADP'}, CONFIG, ['S', 'T', 'Y'])
    a, b = by_source(t, '2aaa'), by_source(t, '2bbb')
    assert (a['record_class'], a['site_state']) == ('binds_modified', 'Sph')
    assert (b['record_class'], b['site_state']) == ('binds_other_residue', 'other_residue')
    # a site of a letter the set does not hold
    with pytest.raises(SystemExit):
        isrc.site_set_records(ev.assign(site_id='P33333:K10'), inputs, {}, CONFIG, ['S', 'T', 'Y'])


@pytest.mark.parametrize('kind', ['rbdmap_fragment', 'predicted_interface', 'pclap_peptide', 'intact_binding_region'])
def test_site_set_records_refuses_other_evidence(tmp_path, kind):
    inputs = write_tables(tmp_path, [entry('1aaa')], [], [])
    ev = pd.DataFrame([ev_row('P11111:K27', 'pdb_rna', 'rna', structures='1aaa'), ev_row('P11111:K27', kind, 'rna', source='RBDmap')], columns=EVIDENCE)
    with pytest.raises(SystemExit):
        isrc.site_set_records(ev, inputs, {}, CONFIG, ['K'])


def test_site_set_tables_use_the_latest_date(tmp_path):
    inputs = write_tables(tmp_path, [entry('1aaa')], [], [])
    newer = pd.DataFrame([entry('9zzz')], columns=ENTRY)
    newer.to_csv(tmp_path / 'interactions' / 'site_set_pdb_entries_2026-12-01.tsv', sep='\t', index=False)
    assert isrc.read_site_set(inputs, 'pdb_entries')['pdb_id'].tolist() == ['9zzz']
    assert isrc.read_site_set(inputs, 'pdb_entries', DATE)['pdb_id'].tolist() == ['1aaa']
    with pytest.raises(SystemExit):
        isrc.read_site_set(inputs, 'pdb_entities_missing')


# ---- scripts/53 in site-set mode --------------------------------------------------------------------------------------

ENTRY_ANSWER = {'data': {'entries': [{
    'rcsb_id': '1ABC', 'struct': {'title': 'a phosphopeptide complex'}, 'exptl': [{'method': 'X-RAY DIFFRACTION'}],
    'rcsb_entry_info': {'resolution_combined': [1.8]}, 'rcsb_accession_info': {'initial_release_date': '2020-02-03T00:00:00Z'},
    'rcsb_primary_citation': {'pdbx_database_id_PubMed': 123, 'pdbx_database_id_DOI': '10.1/abc', 'year': 2020, 'rcsb_authors': ['Doe J', 'Roe R'],
                              'journal_abbrev': 'Nature', 'title': 'a paper'},
    'nonpolymer_entities': [],
    'polymer_entities': [
        {'rcsb_id': '1ABC_1', 'rcsb_polymer_entity': {'pdbx_description': 'substrate'}, 'entity_poly': {'type': 'polypeptide(L)', 'pdbx_seq_one_letter_code': 'MAK(SEP)LKR'},
         'rcsb_entity_source_organism': [{'ncbi_scientific_name': 'Homo sapiens', 'ncbi_taxonomy_id': 9606}],
         'polymer_entity_instances': [{'rcsb_polymer_entity_instance_container_identifiers': {'asym_id': 'A', 'auth_asym_id': 'A'}}],
         'rcsb_polymer_entity_align': [{'reference_database_name': 'UniProt', 'reference_database_accession': 'P11111',
                                        'aligned_regions': [{'entity_beg_seq_id': 1, 'ref_beg_seq_id': 1, 'length': 7}]}]},
        {'rcsb_id': '1ABC_2', 'rcsb_polymer_entity': {'pdbx_description': '14-3-3'}, 'entity_poly': {'type': 'polypeptide(L)', 'pdbx_seq_one_letter_code': 'GGGGGG'},
         'rcsb_entity_source_organism': [{'ncbi_scientific_name': 'Homo sapiens', 'ncbi_taxonomy_id': 9606}],
         'polymer_entity_instances': [{'rcsb_polymer_entity_instance_container_identifiers': {'asym_id': 'B', 'auth_asym_id': 'B'}}],
         'rcsb_polymer_entity_align': [{'reference_database_name': 'UniProt', 'reference_database_accession': 'Q99999',
                                        'aligned_regions': [{'entity_beg_seq_id': 1, 'ref_beg_seq_id': 1, 'length': 6}]}]}]}]}}
CHEM_ANSWER = {'data': {'chem_comps': [{'chem_comp': {'id': 'SEP', 'name': 'PHOSPHOSERINE', 'type': 'L-peptide linking', 'one_letter_code': 'S',
                                                      'mon_nstd_parent_comp_id': ['SER']}}]}}


def test_stage_pdb_site_set_mode(tmp_path, monkeypatch):
    calls = []

    def http(url, data=None, headers=None, tries=4):
        calls.append(url)
        return CHEM_ANSWER if b'chem_comps' in data else ENTRY_ANSWER

    def forbidden(*a, **k):
        raise AssertionError('ModelServer must not be called in site-set mode')
    monkeypatch.setattr(stage, 'http', http)
    monkeypatch.setattr(stage, 'cached', lambda cache, name, fetch: fetch())
    monkeypatch.setattr(stage, 'residue_contacts', forbidden)
    monkeypatch.setattr(stage.urllib.request, 'urlopen', forbidden)
    tables, info = stage.stage_pdb(tmp_path, tmp_path, {'P11111': {4}}, {'P11111': 'MAKSLKR'}, ids=['1ABC'], methylarginine=False)
    assert set(tables) == {'pdb_entries', 'pdb_entities', 'pdb_site_residues'}
    assert calls and all(u == stage.RCSB for u in calls)
    e = tables['pdb_entries'].iloc[0]
    assert (e['pdb_id'], e['title'], e['method'], e['resolution_a'], e['release_date'], e['pubmed'], e['first_author'], e['year']) == \
        ('1abc', 'a phosphopeptide complex', 'X-RAY DIFFRACTION', '1.80', '2020-02-03', '123', 'Doe J', '2020')
    ents = tables['pdb_entities']
    assert ents[['entity_id', 'description', 'uniprot_accessions', 'organisms']].values.tolist() == [['1', 'substrate', 'P11111', 'Homo sapiens'],
                                                                                                      ['2', '14-3-3', 'Q99999', 'Homo sapiens']]
    st = tables['pdb_site_residues']
    assert len(st) == 1
    r = st.iloc[0]
    assert (r['pdb_id'], r['entity_id'], r['uniprot_accession'], r['position'], r['residue'], r['residue_parent'], r['context_check']) == \
        ('1abc', '1', 'P11111', '4', 'SEP', 'SER', 'matches')
    assert info['site_residues_by_parent'] == {'SER': 1} and 'site_residue_not_arginine' not in info


def test_pdbe_kb_entries_only(tmp_path):
    folder = tmp_path / 'interactions'
    folder.mkdir()
    pd.DataFrame([{'uniprot_accession': 'P11111', 'start': '1', 'end': '10', 'pdb_entries': '1abc, 2def'},
                  {'uniprot_accession': 'P11111', 'start': '50', 'end': '60', 'pdb_entries': '3ghi'}]).to_csv(folder / 'pdbe_kb_interface_residues_x.tsv', sep='\t', index=False)
    pd.DataFrame({'pdb_id': ['4jkl']}).to_csv(folder / 'pdb_arginine_contacts_x.tsv', sep='\t', index=False)
    assert stage.pdbe_kb_entries(tmp_path, {'P11111': {4}}) == {'1ABC', '2DEF'}
    assert stage.pdb_entries(tmp_path, {'P11111': {4}}) == ['1ABC', '2DEF', '4JKL']

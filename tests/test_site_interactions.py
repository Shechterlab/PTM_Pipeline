"""Tests of scripts/64_export_site_interactions.py (site-level interaction tables for any modified residue).

Run from the repository root: python -m pytest -q tests/test_site_interactions.py
"""
import csv
import importlib.util
import json
import sys
from pathlib import Path

import pytest

pytest.importorskip('scipy')

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))


def load(name, alias):
    spec = importlib.util.spec_from_file_location(alias, ROOT / 'scripts' / name)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


s64 = load('64_export_site_interactions.py', 'site_interactions_64')
x37 = s64.x37
s37a = load('37a_stage_interactions.py', 'stage_interactions_37a')

KINDS_NOT_WRITTEN = ('rbdmap_fragment', 'pclap_peptide', 'pclap_lysc_region', 'predicted_interface', 'predicted_contact',
                     'intact_mutation', 'intact_binding_region', 'intact_methylated')


def read_tsv(path):
    with open(path, encoding='utf-8', newline='') as fh:
        return list(csv.reader(fh, delimiter='\t'))


def write_tsv(path, header, rows):
    with open(path, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(header)
        w.writerows(rows)


# ----------------------------------------------------------------------------------------------------- 1. the release
def test_methylarginine_reproduces_v1_4_3(tmp_path):
    config = s64.load_config()
    spec = config['site_sets']['methylarginine']
    residues, tables = s64.residues_of(config, spec['residues'])
    topup = tmp_path / 'topup'
    src = s64.Source(ROOT / 'inputs' / 'interactions', topup, residues, tables, 2)
    proteome = s64.reference_proteome(src, ROOT / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03.fasta')
    sites = s64.read_sites(config, 'methylarginine')
    d = s64.compute(sites, src, proteome, spec['noun'], 2, 0, False)
    out = s64.write_tables(d, tmp_path / 'out')
    release = ROOT / 'releases' / 'v1.4.3' / 'supporting'

    for name in ('nucleic_acid_contacts', 'nucleic_acid_near_sites', 'nucleic_acid_proximal_sites'):
        assert out[name].read_bytes() == (release / f'{name}_v1.4.3.tsv').read_bytes(), name

    ref = read_tsv(release / 'site_interactions_v1.4.3.tsv')
    col = ref[0].index('evidence_type')
    kept = [ref[0]] + [r for r in ref[1:] if r[col] not in KINDS_NOT_WRITTEN]
    expected = ''.join('\t'.join(r) + '\n' for r in kept)
    assert out['site_interactions'].read_bytes() == expected.encode('utf-8')

    ref = read_tsv(release / 'residue_evidence_v1.4.3.tsv')
    got = read_tsv(out['residue_evidence'])
    assert got[0] == ref[0]
    col = ref[0].index('evidence_type')
    want = [tuple(r) for r in ref[1:] if r[col] not in KINDS_NOT_WRITTEN]
    got_rows = [tuple(r) for r in got[1:]]
    # every release row of a kind written is present, in the same relative order
    it = iter(got_rows)
    assert all(any(g == w for g in it) for w in want)
    wanted = set(want)
    extra = [r for r in got_rows if r not in wanted]
    # known difference: the release lists a protein's RBR-ID regions only when it has another kind of evidence
    assert {r[col] for r in extra} == {'rbrid_region'}
    assert len(extra) == 52

    assert not src.fetched
    assert not topup.exists() or not any(p.is_file() for p in topup.rglob('*'))


# ------------------------------------------------------------------------------------------------------- 2. the config
def test_config_matches_script_37():
    config = s64.load_config()
    r = config['residues']['R']
    assert frozenset(r['side_chain_atoms']) == frozenset(x37.SIDE_CHAIN)
    assert set(r['modified_components']) == set(x37.METHYL_R)
    residues, tables = s64.residues_of(config, ['R'])
    assert residues is x37.ARGININE
    assert tables == {'contacts': 'pdb_arginine_contacts_', 'chains': 'pdb_chain_arginines_',
                      'distances': 'pdb_arginine_nucleic_acid_distances_', 'chain_column': 'arginines'}

    expected_atoms = {'LYS', 'ALY', 'MLZ', 'MLY', 'M3L', 'SER', 'SEP', 'THR', 'TPO', 'TYR', 'PTR'}
    for letters in (['K'], ['S', 'T', 'Y']):
        residues, tables = s64.residues_of(config, letters)
        assert set(residues.atoms) == expected_atoms
        assert 'ARG' not in residues.atoms
        assert residues.atoms['LYS'] == frozenset({'CD', 'CE', 'NZ'})
        assert residues.atoms['PTR'] == frozenset({'CE1', 'CE2', 'CZ', 'OH'})
        assert (residues.contacts_folder, residues.proximity_folder) == ('residue_contacts', 'na_residues')
        assert tables == {'contacts': 'pdb_residue_contacts_', 'chains': 'pdb_chain_residues_',
                          'distances': 'pdb_residue_nucleic_acid_distances_', 'chain_column': 'residues'}
    with pytest.raises(SystemExit):
        s64.residues_of(config, ['R', 'K'])
    for spec in config['site_sets'].values():
        assert set(spec['residues']) <= set(config['residues'])


# ------------------------------------------------------------------------------------------- 3. synthetic fixtures
# One protein of 60 residues (odd positions E, even G) whose sites sit at these positions, each exercising one thing.
POS = {'protein': 5, 'self': 12, 'rna': 19, 'dna': 26, 'xl': 33, 'near': 40, 'prox': 47, 'dnaev': 54}
MISMATCH = 57
RNA_SEQ = 'ACGUACGUACGUACGUA'
DNA_SEQ = 'GGATCCGGATCCGGATC'
PROTEIN_PARTNER = ('P00009', 'PART', 'MEGEGEGEGEGEGEGEGEGE')
RESN = {'K': 'LYS', 'S': 'SER', 'T': 'THR', 'Y': 'TYR'}


def protein(acc, letter_at):
    seq = ['E' if i % 2 else 'G' for i in range(1, 61)]
    for p, letter in letter_at.items():
        seq[p - 1] = letter
    return ''.join(seq)


def write_fixture(folder, acc, letter_at):
    """The staged tables the steps read, for protein acc with residue letter_at {position: letter} (the sites)."""
    folder.mkdir(parents=True, exist_ok=True)
    seq = protein(acc, letter_at)
    L = lambda role: letter_at[POS[role]]  # noqa: E731
    resn = lambda role: RESN[L(role)]  # noqa: E731
    w = lambda name, header, rows=(): write_tsv(folder / name, header.split(), rows)  # noqa: E731
    p = POS
    w('pdbe_kb_interface_queries_2026-01-01.tsv', 'uniprot_accession answer saved', [[acc, 'interfaces', '2026-01-01']])
    w('pdbe_kb_interface_residues_2026-01-01.tsv', 'uniprot_accession start end partner partner_name pdb_entries', [
        [acc, p['protein'], p['protein'], 'P00009', 'Partner protein', '2prt'],
        [acc, MISMATCH, MISMATCH, 'P00009', 'Partner protein', '2prt'],   # a site whose residue is not the letter: dropped
        [acc, p['self'], p['self'], acc, 'Self', '3slf'],
        [acc, p['rna'], p['rna'], 'RNA', 'RNA', '1rna'],
        [acc, 42, 42, 'RNA', 'RNA', '1rna'],
        [acc, p['dna'], p['dna'], 'DNA', 'DNA', '1dna']])
    seg = lambda pdb, ch, sp, n: [pdb, ch, sp, 1, 1, n]  # noqa: E731
    w('sifts_pdb_chain_uniprot_human_2026-01-01.tsv', 'pdb_id chain uniprot_accession seqres_start uniprot_start uniprot_end', [
        seg('1dna', 'A', acc, 60), seg('1rna', 'A', acc, 60), seg('2prt', 'A', acc, 60), seg('2prt', 'B', 'P00009', 20),
        seg('3slf', 'A', acc, 60), seg('3slf', 'B', acc, 60)])
    w('pdb_nucleic_acid_chains_2026-01-01.tsv',
      'pdb_id chain entity_id nucleic_acid molecule_name source_organism sample_preparation entry_core_histones length sequence', [
          ['1dna', 'D', 2, 'dna', 'Fixture DNA', 'Homo sapiens', 'Synthetically obtained', 'false', len(DNA_SEQ), DNA_SEQ],
          ['1rna', 'C', 2, 'rna', 'Fixture RNA', 'Homo sapiens', 'Synthetically obtained', 'false', len(RNA_SEQ), RNA_SEQ]])
    # residues of the protein chains against RNA and DNA (a residue i+2 of the near site touches the RNA at 42)
    w('pdb_nucleic_acid_residue_contacts_2026-01-01.tsv',
      'pdb_id chain seqres_index author_number residue na_chain na_seqres_index nucleotide distance_a', [
          ['1dna', 'A', '', '', '', '', '', '', ''],
          ['1rna', 'A', 42, 42, 'GLY', 'C', 9, 'A', 3.5]])
    w('pdb_residue_contacts_2026-01-01.tsv',
      'pdb_id chain seqres_index residue contact_chain contact_seqres_index contact_author_number contact_residue distance_a', [
          ['1dna', 'A', p['dna'], resn('dna'), 'D', 4, 4, 'DT', 3.4],
          ['1rna', 'A', p['rna'], resn('rna'), 'C', 6, 6, 'C', 3.2],
          ['1rna', 'A', p['rna'], resn('rna'), 'C', 7, 7, 'G', 3.8],
          ['2prt', 'A', p['protein'], resn('protein'), 'B', 3, 3, 'ASP', 3.1],
          ['3slf', 'A', p['self'], resn('self'), 'B', p['self'], p['self'], resn('self'), 3.6],
          ['3slf', 'B', p['self'], resn('self'), 'A', p['self'], p['self'], resn('self'), 3.6]])
    w('pdb_chain_residues_2026-01-01.tsv', 'pdb_id chain n_residues residues', [
        ['1dna', 'A', 0, ''],
        ['1rna', 'A', 2, f"{p['near']}:{p['near']}:{resn('near')}; {p['prox']}:{p['prox']}:{resn('prox')}"]])
    w('pdb_residue_nucleic_acid_distances_2026-01-01.tsv',
      'pdb_id chain seqres_index na_chain closest_na_seqres_index closest_nucleotide side_chain_distance_a na_seqres_within_8a',
      [['1rna', 'A', p['prox'], 'C', 12, 'U', 6.1, '10; 11; 12; 13']])
    # evidence from cells
    x, d = p['xl'], p['dnaev']
    pep = seq[x - 2:x + 1]
    w('pepr_ms_v2_table_s2_rna_crosslinked_peptides.tsv', 'dataset uniprot_accession protein_ids position crosslinked_residue peptide modified_sequence rna_mass',
      [['Kramer 2014', acc, acc, x, L('xl'), pep, f'{pep[0]}{pep[1]}[324.0359]{pep[2]}', '324.0359']])
    w('pepr_ms_v2_table_s3_focused_sites.tsv', 'uniprot_accession protein_id focused_residues focused_positions')
    w('rbdmap_castello2016_rna_bound_fragments.tsv', 'protein_id gene fragment_start fragment_end fragment_sequence',
      [[acc, 'GENA', x - 5, x + 5, seq[x - 6:x + 5]]])
    w('capri_panhale2019_xl_peptides.tsv', 'gene uniprot_id peptide_start peptide_end xl_peptide', [['GENA', acc, x - 1, x + 1, pep]])
    w('capri_panhale2019_adj_peptides.tsv', 'gene uniprot_id peptide_start peptide_end adj_peptide')
    dpep = seq[d - 2:d + 1]
    w('trendel2025_dna_crosslink_localisation.tsv', 'uniprot_accession uniprot_id peptide_start sequence_with_localization',
      [[acc, acc, d - 1, dpep[0] + dpep[1].lower() + dpep[2]]])
    w('trendel2025_dna_crosslinked_peptides.tsv', 'uniprot_accession uniprot_id peptide_start peptide_end', [[acc, acc, d - 1, d + 1]])
    w('trendel2025_mcf7_dna_interactome.tsv', 'uniprot_accession uniprot_id')
    w('uniprot_human_dna_binding_regions_2026-01-01.tsv', 'uniprot_accession start end note', [[acc, d - 2, d + 2, '']])
    w('rbrid_oksuz2023_k562_positive_residues.tsv',
      'uniprot_accession gene position residue rbrid_positive_peptides best_p published_rule_peptides best_p_published summed_score residue_score',
      [[acc, 'GENA', x, L('xl'), 1, 0.001, 1, 0.001, 5.0, 4.2], [acc, 'GENA', x + 1, seq[x], 1, 0.001, 1, 0.001, 4.0, 3.1]])
    w('rbrid_oksuz2023_k562_coverage.tsv', 'uniprot_accession gene start end coverage length', [[acc, 'GENA', 1, 60, 'quantified', 60]])
    w('rbrid_oksuz2023_k562_segments.tsv',
      'uniprot_accession gene start end length max_residue_score mean_residue_score max_summed_score peptides rbrid_positive_peptides best_p '
      'published_rule_peptides protein_length protein_scored_residues residue_score_threshold',
      [[acc, 'GENA', x - 2, x + 3, 6, 4.2, 3.5, 6.0, 1, 1, 0.001, 1, 60, 6, 3]])
    w('rbrid_oksuz2023_k562_peptides.tsv', 'uniprot_accession gene start end peptide rbrid_positive log2_ratio p',
      [[acc, 'GENA', x - 1, x + 2, seq[x - 2:x + 2], 'true', -1.5, 0.001]])
    return seq


def run_fixture(tmp_path, monkeypatch, letters, letter_at, extra_sites=()):
    for name in ('fetch_interfaces', 'site_contacts', 'chain_na_contacts', 'chain_na_proximity', 'load_nucleic_chains'):
        def refuse(*a, **k):
            raise AssertionError('network')
        monkeypatch.setattr(s64.x37, name, refuse)
    acc = 'P00001'
    folder = tmp_path / 'staged'
    seq = write_fixture(folder, acc, letter_at)
    config = s64.load_config()
    residues, tables = s64.residues_of(config, letters)
    topup = tmp_path / 'topup'
    src = s64.Source(folder, topup, residues, tables, 1)
    proteome = {acc: ('GENA', '9606.ENSP1', seq), PROTEIN_PARTNER[0]: ('PART', '9606.ENSP9', PROTEIN_PARTNER[2])}
    sites = [(acc, pos, letter) for pos, letter in sorted(letter_at.items())] + [(acc, MISMATCH, letters[0])]
    d = s64.compute(sites, src, proteome)
    out = s64.write_tables(d, tmp_path / 'out')
    assert not src.fetched
    return {name: read_tsv(path) for name, path in out.items()}


def rows_of(table, **common):
    """The data rows of a table as dicts."""
    head, *body = table
    return [dict(zip(head, r)) for r in body]


def expected(letter_at, table, *rows):
    """Expected rows as dicts over the columns of the release table; unnamed columns are empty."""
    head = s64.release_columns()[table]
    return [{c: str(r.get(c, '')) for c, *_ in head} for r in rows]


def check_tables(t, letter_at):
    L = lambda role: letter_at[POS[role]]  # noqa: E731
    sid = lambda role: f'P00001:{L(role)}{POS[role]}'  # noqa: E731
    x, k = POS['xl'], L('xl')
    base = dict(source='PDBe-KB', counts_as_evidence='true', localisation='at_residue')
    assert rows_of(t['site_interactions']) == expected(letter_at, 'site_interactions',
        dict(base, site_id=sid('protein'), evidence_type='pdb_protein', partner_kind='protein', partner_accession='P00009', partner_gene='PART',
             structures='2prt', detail='1 PDB entries'),
        dict(base, site_id=sid('self'), evidence_type='pdb_self', partner_kind='self', structures='3slf', detail='1 PDB entries'),
        dict(base, site_id=sid('rna'), evidence_type='pdb_rna', partner_kind='rna', structures='1rna', detail='1 PDB entries'),
        dict(base, site_id=sid('dna'), evidence_type='pdb_dna', partner_kind='dna', structures='1dna', detail='1 PDB entries'),
        dict(base, site_id=sid('xl'), evidence_type='rna_crosslink', partner_kind='rna', source='Kramer2014', detail=f'G[{k}]G at {x} (324.0359 Da)'),
        dict(base, site_id=sid('xl'), evidence_type='capri_xl_peptide', localisation='in_peptide', partner_kind='rna', source='CAPRI',
             detail=f'{x - 1}-{x + 1} G{k}G'),
        dict(base, site_id=sid('xl'), evidence_type='rbrid_residue', partner_kind='rna', source='RBR-ID',
             detail=f'{x - 1}-{x + 2} G{k}GE (log2 +4SU/-4SU -1.50, p 0.001)'),
        dict(base, site_id=sid('dnaev'), evidence_type='dna_crosslink', partner_kind='dna', source='Trendel2025',
             detail='peptide 53-55 E%sE (lower case: the residues the crosslink may sit on)' % L('dnaev').lower()),
        dict(base, site_id=sid('dnaev'), evidence_type='dna_crosslink', localisation='in_peptide', partner_kind='dna', source='Trendel2025',
             counts_as_evidence='false',
             detail='inside a DNA-crosslinked peptide. peptide 53-55 E%sE (lower case: the residues the crosslink may sit on)' % L('dnaev').lower()),
        dict(base, site_id=sid('dnaev'), evidence_type='dna_binding_region', localisation='in_region', partner_kind='dna', source='UniProt',
             counts_as_evidence='false'))
    ev = lambda s, e, kind, loc, src, n='': dict(uniprot_accession='P00001', start=s, end=e, evidence_type=kind, localisation=loc, sources=src, n_partners=n)  # noqa: E731
    assert rows_of(t['residue_evidence']) == expected(letter_at, 'residue_evidence',
        ev(5, 5, 'pdb_protein', 'at_residue', 'PDBe-KB', 1), ev(12, 12, 'pdb_self', 'at_residue', 'PDBe-KB'),
        ev(19, 19, 'pdb_rna', 'at_residue', 'PDBe-KB'), ev(26, 26, 'pdb_dna', 'at_residue', 'PDBe-KB'),
        ev(31, 36, 'rbrid_region', 'in_region', 'RBR-ID'), ev(32, 34, 'capri_xl_peptide', 'in_peptide', 'CAPRI'),
        ev(33, 33, 'rna_crosslink', 'at_residue', 'Kramer2014'), ev(33, 34, 'rbrid_residue', 'in_region', 'RBR-ID'),
        ev(42, 42, 'pdb_rna', 'at_residue', 'PDBe-KB'), ev(52, 56, 'dna_binding_region', 'in_region', 'UniProt'),
        ev(53, 55, 'dna_crosslink', 'in_peptide', 'Trendel2025'), ev(54, 54, 'dna_crosslink', 'at_residue', 'Trendel2025'),
        ev(57, 57, 'pdb_protein', 'at_residue', 'PDBe-KB', 1))
    assert not any(r['site_id'].endswith(f':{L("protein")}{MISMATCH}') or r['site_id'].endswith(str(MISMATCH)) for r in rows_of(t['site_interactions']))
    # the RBDmap fragment (28-38) covers the RNA-crosslinked site and is not written, nor is any pCLAP row
    for name in ('site_interactions', 'residue_evidence'):
        assert not {r['evidence_type'] for r in rows_of(t[name])} & set(KINDS_NOT_WRITTEN)

    # the closest nucleotide is C6 (3.2 A) of the RNA, with G7 also contacted; T4 of the DNA
    na = dict(nucleic_acid_type='rna_oligo', molecule_name='Fixture RNA', source_organism='Homo sapiens', sample_preparation='Synthetically obtained')
    assert rows_of(t['nucleic_acid_contacts']) == expected(letter_at, 'nucleic_acid_contacts',
        dict(na, site_id=sid('rna'), pdb_id='1rna', chain='C', nucleic_acid='rna', chain_length=17, closest_nucleotide='C6', closest_distance_a='3.2',
             contacted_nucleotides='C6; G7', contacted_author_numbers='6; 7', contacted_distances_a='3.20; 3.80', window_start=-2,
             window='---ACGUACGUACGUAC', window_marked='---ACGUA[C][G]UACGUAC'),
        dict(na, site_id=sid('dna'), pdb_id='1dna', chain='D', nucleic_acid='dna', nucleic_acid_type='other_dna', molecule_name='Fixture DNA',
             chain_length=17, closest_nucleotide='T4', closest_distance_a='3.4', contacted_nucleotides='T4', contacted_author_numbers='4',
             contacted_distances_a='3.40', window_start=-4, window='-----GGATCCGGATCC', window_marked='-----GGA[T]CCGGATCC'))
    # the site at 40 touches nothing itself; residue 42 (i+2) touches A9. The site at 47 is also within 5 of 42 but is proximal, so it is
    # left out here
    assert rows_of(t['nucleic_acid_near_sites']) == expected(letter_at, 'nucleic_acid_near_sites',
        dict(na, site_id=sid('near'), nucleic_acid='rna', residue_offset=2, residue=42, pdb_id='1rna', n_entries=1, chain='C',
             centre_nucleotide='A9', window_start=1, window=RNA_SEQ, window_marked='ACGUACGU[A]CGUACGUA'))
    assert rows_of(t['nucleic_acid_proximal_sites']) == expected(letter_at, 'nucleic_acid_proximal_sites',
        dict(na, site_id=sid('prox'), nucleic_acid='rna', distance_a='6.1', pdb_id='1rna', n_entries=1, chain='C', closest_nucleotide='U12',
             window_start=4, window='UACGUACGUACGUA---', window_marked='UACGUA[C][G][U][A]CGUA---'))


def test_lysine_fixture(tmp_path, monkeypatch):
    letter_at = {p: 'K' for p in POS.values()}
    check_tables(run_fixture(tmp_path, monkeypatch, ['K'], letter_at), letter_at)


def test_phospho_fixture(tmp_path, monkeypatch):
    letter_at = dict(zip(POS.values(), 'SYSYSYSY'))
    check_tables(run_fixture(tmp_path, monkeypatch, ['S', 'T', 'Y'], letter_at), letter_at)


# ------------------------------------------------------------------------------------------------ 4. staging round trip
def test_stage_topup_residue_tables(tmp_path):
    cache = tmp_path / 'cache'
    (cache / 'residue_contacts').mkdir(parents=True)
    (cache / 'na_residues').mkdir()
    (cache / 'residue_contacts' / '1abc_A_3.json').write_text(json.dumps({'resn': 'LYS', 'hits': [['B', '3', '3', 'U', 3.1]]}))
    (cache / 'residue_contacts' / '1abc_A_9.json').write_text(json.dumps({'resn': 'SER', 'hits': []}))
    (cache / 'na_residues' / '1abc_A.json').write_text(json.dumps({'args': [[5, '5', 'LYS']], 'hits': [[5, 'B', 3, 'U', 6.1, [2, 3]]]}))
    out = tmp_path / 'staged'
    out.mkdir()
    manifest = tmp_path / 'MANIFEST.tsv'
    header = (ROOT / 'inputs' / 'MANIFEST.tsv').read_text(encoding='utf-8').splitlines()[0]
    manifest.write_text(header + '\n', encoding='utf-8')
    s37a.stage_topup(cache, out, 1, manifest, 'the test sites')
    names = sorted(p.name for p in out.glob('*.tsv'))
    stem = lambda prefix: [n for n in names if n.startswith(prefix)]  # noqa: E731
    assert len(stem('pdb_residue_contacts_')) == len(stem('pdb_chain_residues_')) == len(stem('pdb_residue_nucleic_acid_distances_')) == 1
    listed = [ln.split('\t')[0] for ln in manifest.read_text(encoding='utf-8').splitlines()[1:]]
    assert f'interactions/{stem("pdb_residue_contacts_")[0]}' in listed

    # Source also reads the interface tables; header only here
    write_tsv(out / 'pdbe_kb_interface_queries_2026-01-01.tsv', ['uniprot_accession', 'answer', 'saved'], [])
    config = s64.load_config()
    residues, tables = s64.residues_of(config, ['K'])
    src = s64.Source(out, tmp_path / 'topup', residues, tables, 1)
    assert src.contacts == {('1abc', 'A', 3): {'resn': 'LYS', 'hits': [('B', '3', '3', 'U', 3.1)]},
                            ('1abc', 'A', 9): {'resn': 'SER', 'hits': []}}
    assert src.proximity == {('1abc', 'A'): {'args': [[5, '5', 'LYS']], 'hits': [[5, 'B', 3, 'U', 6.1, [2, 3]]]}}
    # the same answers through the methods the steps call (no fetch)
    assert src.site_contacts('1abc', 'A', 3) == {'resn': 'LYS', 'hits': [['B', '3', '3', 'U', 3.1]]}
    assert src.chain_na_proximity('1abc', 'A', {'B'}) == src.proximity[('1abc', 'A')]
    assert not src.fetched

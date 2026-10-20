"""Edge cases of the remapping (scripts/02) and of its validation gates, on small hand-made sequences.

    python -m pytest tests/        or        python tests/test_remap_edge_cases.py
"""
from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))


def load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


remap = load('remap', '02_remap_arg_methyl_to_canonical.py')
evidence_check = load('evidence_check', 'qc_remap_evidence_check.py')

FILLER = 'LLKKAAPPQQEEDDSSTTVVWWYYLLKKAAPPQQEEDDSSTTVVWWYY'


def row(**fields) -> pd.Series:
    base = dict(substrate_UniProtAC='P1', substrate_UniProtAC_original='P1', canonical_UniProtAC='P1', substrate_genename='G',
                position=1, seq_window='', peptide_sequence='', peptides='', accession_is_isoform=False,
                seq_window_window_length=0, peptide_sequence_window_length=0, source_table='maron_s5', reference='X')
    base.update(fields)
    return pd.Series(base)


def place(r: pd.Series, sequences: dict, gene_map: dict | None = None) -> tuple[str, str]:
    out = remap.remap_row(r, sequences, gene_map or {})
    return (f"{out['corrected_canonical_UniProtAC']}:{out['corrected_site']}" if out['is_position_resolved'] else ''), out['remap_status']


# ---- windows cut short at a protein end ---------------------------------------------------------------------------

def test_n_terminal_window_reads_from_its_intact_end():
    seq = 'MAR' + 'GSTAR' + FILLER  # R3 and R8
    window = seq[:10]  # 15-residue window around R3, cut at the N-terminus
    for reported in (3, 4, 8):
        assert place(row(position=reported, peptide_sequence=window, peptide_sequence_window_length=15), {'P1': seq})[0] == 'P1:R3'


def test_c_terminal_window_reads_from_its_intact_end():
    seq = FILLER + 'GSRRAAAARGS'  # both readings of the last 10 residues are arginines: len-7 and len-2
    window = seq[-10:]  # 15-residue window around the R at len-2, cut at the C-terminus
    site = len(seq) - 2
    assert place(row(position=site - 6, peptide_sequence=window, peptide_sequence_window_length=15), {'P1': seq})[0] == f'P1:R{site}'


def test_window_inside_a_longer_entry_mature_chain():
    # Histone H4: sources number the mature chain (SGR3), UniProt counts the initiator methionine (R4).
    seq = 'MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG'
    assert place(row(position=3, peptide_sequence='SGRGKGGKGL', peptide_sequence_window_length=15), {'P1': seq})[0] == 'P1:R4'


def test_window_inside_a_longer_entry_with_two_readings_needs_the_reported_position():
    # 16 residues of a 31-residue window, inside the entry, with an arginine at both readings (window index 1 if cut
    # at its N-terminal side, 16 if cut at its C-terminal side): only the reported position decides.
    window = 'R' + 'A' * 14 + 'R'
    seq = 'MKKK' + window + FILLER
    kwargs = dict(peptide_sequence=window, peptide_sequence_window_length=31)
    assert place(row(position=5, **kwargs), {'P1': seq})[0] == 'P1:R5'  # this entry's numbering, n-cut reading
    assert place(row(position=20, **kwargs), {'P1': seq})[0] == 'P1:R20'  # this entry's numbering, c-cut reading
    assert place(row(position=1, **kwargs), {'P1': seq})[0] == 'P1:R5'  # source numbering from the window's first residue
    assert place(row(position=12, **kwargs), {'P1': seq})[0] == ''


def test_padding_gives_the_cut_side():
    seq = FILLER + 'VMSQGGRY'
    site = len(seq) - 1
    guo = row(position=site, peptide_sequence='VMSQGGRYXXXXX', peptide_sequence_window_length=13)
    assert place(guo, {'P1': seq}) == (f'P1:R{site}', 'sequence_unique_canonical')
    prometheus = row(source_table='prometheus_mmc2', position=5, seq_window='___MKEPRIFPRERP', seq_window_window_length=15)
    assert place(prometheus, {'P1': 'MKEPRIFPRERPAAQQ' + FILLER})[0] == 'P1:R5'


# ---- several candidate entries ------------------------------------------------------------------------------------

def test_equal_support_on_two_other_entries_stays_ambiguous():
    seq = 'MSTRAGGPLLR' + 'K' * 30
    r = row(substrate_UniProtAC='T0', substrate_UniProtAC_original='T0', canonical_UniProtAC='T0', position=13, peptide_sequence='GGPLLRK;TRAGG')
    for order in (['A1', 'B2'], ['B2', 'A1']):
        assert place(r, {'A1': seq, 'B2': seq}, {'G': order}) == ('', 'unresolved_ambiguous_hit')


def test_own_accession_breaks_a_tie():
    seq = 'MSTRAGGPLLR' + 'K' * 30
    r = row(substrate_UniProtAC='B2', substrate_UniProtAC_original='B2', canonical_UniProtAC='B2', position=13,
            peptide_sequence='GGPLLRK;TRAGG')
    assert place(r, {'A1': seq + 'GGPLLRK', 'B2': seq}, {'G': ['A1', 'B2']})[0] == 'B2:R11'


# ---- the reported position without placing evidence ----------------------------------------------------------------

def test_peptide_elsewhere_in_the_entry_does_not_back_the_reported_arginine():
    seq = 'RRR' + 'AAAAA' + 'R' + 'A' * 20
    assert place(row(position=9, peptide_sequence='RRR'), {'P1': seq}) == ('', 'unresolved_evidence_elsewhere')


def test_window_not_in_the_entry_backs_the_reported_arginine_only_when_it_nearly_matches():
    seq = FILLER[:20] + 'GSKPAGRGGAPLE' + FILLER
    site = 20 + 7
    near = 'GSKPAGRGGAPLD'  # one difference: a variant or sequence conflict
    far = 'WSKPWGRGGWPLD'  # four differences
    assert place(row(position=site, peptide_sequence=near, peptide_sequence_window_length=13), {'P1': seq}) == (f'P1:R{site}', 'same_accession_evidence_near_match')
    assert place(row(position=site, peptide_sequence=far, peptide_sequence_window_length=13), {'P1': seq}) == ('', 'unresolved_evidence_mismatch')


def test_short_peptide_needs_more_than_two_matching_residues():
    seq = 'AAAAARAAA' + 'A' * 20
    assert place(row(position=6, peptide_sequence='KRK'), {'P1': seq}) == ('', 'unresolved_evidence_mismatch')


def test_position_only_row_keeps_the_reported_arginine():
    seq = 'MAAAR' + FILLER
    assert place(row(source_table='human_ptm_master', reference='', position=5), {'P1': seq}) == ('P1:R5', 'same_accession_position_verified')


# ---- accessions and symbols that are no longer current (UniProt sec_ac.txt and gene synonyms) ------------------------

def place_history(r: pd.Series, sequences: dict, gene_map: dict | None = None, secondary: dict | None = None,
                  synonyms: dict | None = None) -> tuple[str, str]:
    out = remap.remap_row(r, sequences, gene_map or {}, secondary or {}, synonyms or {}, remap.ProteomeIndex(sequences))
    return (f"{out['corrected_canonical_UniProtAC']}:{out['corrected_site']}" if out['is_position_resolved'] else ''), out['remap_status']


WINDOW_SEQ = 'MSTRAGGPLLR' + 'K' * 30  # 13-residue window AGGPLLRKKKKKK centred on R11


def old_row(acc: str, gene: str = 'OLDSYM', **fields) -> pd.Series:
    return row(substrate_UniProtAC=acc, substrate_UniProtAC_original=acc, canonical_UniProtAC=acc, substrate_genename=gene, **fields)


def test_merged_accession_is_placed_on_its_current_entry_and_wins_a_tie():
    r = old_row('T9', position=11, peptide_sequence='AGGPLLRKKKKKK', peptide_sequence_window_length=13)
    seqs = {'A1': WINDOW_SEQ, 'P2': WINDOW_SEQ}
    # Without the history the window names two entries equally and nothing is placed; with it, P2 is the row's own.
    assert place_history(r, seqs, {'OLDSYM': ['A1', 'P2']}) == ('', 'unresolved_ambiguous_hit')
    assert place_history(r, seqs, {'OLDSYM': ['A1', 'P2']}, {'T9': ['P2']}) == ('P2:R11', 'sequence_unique_canonical')


def test_split_accession_goes_where_its_evidence_is():
    r = old_row('T8', position=11, peptide_sequence='AGGPLLRKKKKKK', peptide_sequence_window_length=13)
    seqs = {'P3': 'M' + 'A' * 40, 'P4': WINDOW_SEQ}
    assert place_history(r, seqs, secondary={'T8': ['P3', 'P4']}) == ('P4:R11', 'sequence_unique_candidate')


def test_split_accession_into_identical_entries_follows_the_gene_symbol():
    # HSPA1A P08107 was split into P0DMV8 (HSPA1A) and P0DMV9 (HSPA1B), whose sequences are identical.
    r = old_row('T8', gene='HSPA1A', position=11, peptide_sequence='AGGPLLRKKKKKK', peptide_sequence_window_length=13)
    seqs = {'P3': WINDOW_SEQ, 'P4': WINDOW_SEQ}
    assert place_history(r, seqs, {'HSPA1A': ['P3'], 'HSPA1B': ['P4']}, {'T8': ['P3', 'P4']}) == ('P3:R11', 'sequence_unique_canonical')
    assert place_history(old_row('T8', gene='OTHER', position=11, peptide_sequence='AGGPLLRKKKKKK', peptide_sequence_window_length=13),
                         seqs, {}, {'T8': ['P3', 'P4']}) == ('', 'unresolved_ambiguous_hit')


def test_outdated_symbol_is_found_through_its_synonym():
    r = old_row('T7', gene='MST4', position=11, peptide_sequence='AGGPLLRKKKKKK', peptide_sequence_window_length=13)
    seqs = {'P5': WINDOW_SEQ, 'Q1': 'M' + 'A' * 40}
    assert place_history(r, seqs, {'STK26': ['P5']}, synonyms={'MST4': ['P5']}) == ('P5:R11', 'sequence_unique_candidate')


def test_deleted_accession_is_placed_only_on_a_unique_proteome_match():
    r = old_row('T6', gene='GONE', position=11, peptide_sequence='AGGPLLRKKKKKK', peptide_sequence_window_length=13)
    assert place_history(r, {'P6': WINDOW_SEQ, 'Q2': 'M' + 'A' * 40}) == ('P6:R11', 'sequence_unique_proteome')
    assert place_history(r, {'P6': WINDOW_SEQ, 'Q2': WINDOW_SEQ})[0] == ''


def test_proteome_search_never_moves_a_row_to_a_paralog_of_its_named_gene():
    # ELMO1's window matching only ELMO2: the symbol names a current entry, so the row stays unplaced (issue 9).
    r = old_row('T5', gene='ELMO1', position=11, peptide_sequence='AGGPLLRKKKKKK', peptide_sequence_window_length=13)
    assert place_history(r, {'Q1': 'M' + 'A' * 40, 'Q2': WINDOW_SEQ}, {'ELMO1': ['Q1'], 'ELMO2': ['Q2']})[0] == ''


def test_merged_accession_without_evidence_keeps_no_position():
    # Numbering of the old entry can differ from the current one (A8MZ91 vs PRMT5 O14744), so a bare position is not moved.
    r = old_row('T9', source_table='human_ptm_master', reference='', position=5)
    assert place_history(r, {'P2': 'MAAAR' + FILLER}, secondary={'T9': ['P2']})[0] == ''


# ---- peptide starts from the source's coordinate table (Maron 2021 Tables S2/S3) ----------------------------------

COORD_PEPTIDE = 'GSRGSAAPRGS'  # arginines at offsets 2 and 8
COORD_SEQ = FILLER + COORD_PEPTIDE + FILLER
COORD_START = len(FILLER) + 1  # where the peptide starts in the entry; the source numbers it from 100


def write_coordinates(path: Path, rows: list[tuple[str, str, int]]) -> Path:
    pd.DataFrame([{'source': 'test', 'accession': a, 'gene': 'G', 'peptide': p, 'peptide_start': s, 'canonical_UniProtAC': a}
                  for a, p, s in rows]).to_csv(path, sep='\t', index=False)
    return path


def test_multi_arginine_peptide_is_counted_from_its_start(tmp_path):
    coords = remap.load_peptide_coordinates(write_coordinates(tmp_path / 'c.tsv', [('P1', COORD_PEPTIDE, 100)]), {'P1': COORD_SEQ})
    assert coords == {'P1': {(COORD_PEPTIDE, 100, COORD_START)}}  # found once in the entry, away from its stated start
    out = remap.remap_row(row(position=108, peptide_sequence=COORD_PEPTIDE), {'P1': COORD_SEQ}, {}, coordinates=coords)
    assert (out['remap_status'], out['corrected_position']) == ('sequence_peptide_start', COORD_START + 8)
    # a peptide inside the known one, or containing it, is counted from its own place
    out = remap.remap_row(row(position=108, peptide_sequence=COORD_PEPTIDE[1:]), {'P1': COORD_SEQ}, {}, coordinates=coords)
    assert out['corrected_position'] == COORD_START + 8
    longer = COORD_SEQ[COORD_START - 4:COORD_START + len(COORD_PEPTIDE)]
    out = remap.remap_row(row(position=102, peptide_sequence=longer), {'P1': COORD_SEQ}, {}, coordinates=coords)
    assert out['corrected_position'] == COORD_START + 2


def test_peptide_start_needs_an_arginine_at_the_reported_offset_and_a_known_peptide(tmp_path):
    coords = {'P1': {(COORD_PEPTIDE, 100, COORD_START)}}
    assert not remap.remap_row(row(position=105, peptide_sequence=COORD_PEPTIDE), {'P1': COORD_SEQ}, {}, coordinates=coords)['is_position_resolved']
    other = 'GSRGSAAPRGT'  # not the known peptide, nor inside or around it
    seq = COORD_SEQ + other
    assert not remap.remap_row(row(position=108, peptide_sequence=other), {'P1': seq}, {}, coordinates=coords)['is_position_resolved']
    # a peptide found twice in the entry and not at its stated start gives no coordinates
    assert remap.load_peptide_coordinates(write_coordinates(tmp_path / 'c.tsv', [('P1', COORD_PEPTIDE, 100)]),
                                          {'P1': COORD_SEQ + COORD_PEPTIDE}) == {}


def test_evidence_check_accepts_a_peptide_start_placement_only_with_the_coordinates(tmp_path):
    remapped = tmp_path / 'root' / 'remapped'
    remapped.mkdir(parents=True)
    pd.DataFrame([{
        'source_table': 'maron_s5', 'reference': 'This study.', 'substrate_UniProtAC_original': 'P12345', 'substrate_genename': 'GP12345',
        'position': 108, 'seq_window': '', 'peptide_sequence': COORD_PEPTIDE, 'peptides': '', 'corrected_canonical_UniProtAC': 'P12345',
        'corrected_position': COORD_START + 8, 'is_position_resolved': True, 'remap_status': 'sequence_peptide_start',
    }]).to_csv(remapped / 'human_arg_methyl_source_rows_with_remap.tsv', sep='\t', index=False)
    write_fasta(tmp_path / 'ref.fasta', {'P12345': COORD_SEQ})
    args = ['--results-root', str(tmp_path / 'root'), '--fasta', str(tmp_path / 'ref.fasta'), '--fail-on-contradictions']
    assert run('qc_remap_evidence_check.py', *args).returncode == 1  # a peptide with two arginines names neither
    coords = write_coordinates(tmp_path / 'c.tsv', [('P12345', COORD_PEPTIDE, 100)])
    checked = run('qc_remap_evidence_check.py', *args, '--peptide-coordinates', str(coords))
    assert checked.returncode == 0, checked.stdout + checked.stderr


# ---- the validation gates fail the run -----------------------------------------------------------------------------

def write_fasta(path: Path, sequences: dict) -> None:
    path.write_text(''.join(f'>sp|{acc}|{acc}_HUMAN Test OS=Homo sapiens GN=G{acc}\n{seq}\n' for acc, seq in sequences.items()))


def results_root(tmp_path: Path, placed_position: int) -> Path:
    seq = 'MAAAR' + 'GGRGG' + FILLER
    remapped = tmp_path / 'root' / 'remapped'
    remapped.mkdir(parents=True)
    pd.DataFrame([{
        'source_table': 'human_ptm_master', 'reference': '', 'substrate_UniProtAC_original': 'P12345', 'substrate_genename': 'GP12345',
        'position': 5, 'seq_window': '', 'peptide_sequence': '', 'peptides': '', 'corrected_canonical_UniProtAC': 'P12345',
        'corrected_position': placed_position, 'is_position_resolved': True, 'remap_status': 'same_accession_position_verified',
    }]).to_csv(remapped / 'human_arg_methyl_source_rows_with_remap.tsv', sep='\t', index=False)
    integrated = tmp_path / 'root' / 'integrated'
    integrated.mkdir()
    pd.DataFrame([{'canonical_UniProtAC': 'P12345', 'corrected_position': placed_position}]).to_csv(
        integrated / 'human_arg_methyl_union_dedup_by_site.tsv', sep='\t', index=False)
    write_fasta(tmp_path / 'ref.fasta', {'P12345': seq})
    return tmp_path / 'root'


def run(script: str, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(ROOT / 'scripts' / script), *args], capture_output=True, text=True)


def test_evidence_check_fails_on_a_contradiction(tmp_path):
    root = results_root(tmp_path, placed_position=8)  # placed on R8; the position-only evidence says R5
    args = ['--results-root', str(root), '--fasta', str(tmp_path / 'ref.fasta')]
    assert run('qc_remap_evidence_check.py', *args).returncode == 0  # report only
    gated = run('qc_remap_evidence_check.py', *args, '--fail-on-contradictions')
    assert gated.returncode == 1, gated.stdout + gated.stderr
    allow = tmp_path / 'allow.tsv'
    allow.write_text('source_table\treported_accession\treported_position\tplaced_site\treason\nhuman_ptm_master\tP12345\t5\tP12345:8\ttest\n')
    assert run('qc_remap_evidence_check.py', *args, '--fail-on-contradictions', '--allow', str(allow)).returncode == 0


def test_evidence_check_passes_an_agreeing_root(tmp_path):
    root = results_root(tmp_path, placed_position=5)
    assert run('qc_remap_evidence_check.py', '--results-root', str(root), '--fasta', str(tmp_path / 'ref.fasta'), '--fail-on-contradictions').returncode == 0


def test_residue_audit_fails_on_a_non_arginine(tmp_path):
    root = results_root(tmp_path, placed_position=6)  # residue 6 is G
    out = tmp_path / 'audit.tsv'
    audit = run('qc_remap_hand_check.py', '--audit-union', '--results-root', str(root), '--fasta', str(tmp_path / 'ref.fasta'), '--audit-out', str(out))
    assert audit.returncode == 1, audit.stdout + audit.stderr
    root = results_root(tmp_path / 'ok', placed_position=5)
    assert run('qc_remap_hand_check.py', '--audit-union', '--results-root', str(root), '--fasta', str(tmp_path / 'ok' / 'ref.fasta'),
               '--audit-out', str(out)).returncode == 0


if __name__ == '__main__':
    import inspect
    import tempfile
    failed = 0
    for name, test in sorted(globals().items()):
        if not name.startswith('test_'):
            continue
        try:
            if 'tmp_path' in inspect.signature(test).parameters:
                with tempfile.TemporaryDirectory() as tmp:
                    test(Path(tmp))
            else:
                test()
            print(f'ok    {name}')
        except AssertionError as error:
            failed += 1
            print(f'FAIL  {name}: {error}')
    sys.exit(1 if failed else 0)

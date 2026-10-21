#!/usr/bin/env python3
"""Stage the methylarginine records of each study and database from its own tables (inputs/sources/).

    python scripts/00_stage_study_sites.py --study STUDY --input TABLE [TABLE ...] [--reference FASTA] \
        [--entries ENTRIES --history HISTORY] --out inputs/sources/STUDY_methylarginine_records.tsv

STUDY is one of the keys of STUDIES below; the command that staged each table, with its inputs and their MD5s, is in
inputs/MANIFEST.tsv (how_made_or_checked, derived_from). --reference is the reviewed FASTA (lu2025, fedoriw2019, musiani2019),
--entries and --history the reference entries and accession history (guo2014, fedoriw2019, musiani2019: gene symbols to
entries).

Every record says how the methyl group is localised to its arginine (localised_by): localisation_score (the study's
score passes its threshold), one_methylatable_residue or every_arginine_dimethylated (the peptide's mass leaves one
arrangement, common.mass_localised, for studies that give no score), or curated_experimental (UniProt ECO:0000269).
Rows without one are dropped here and counted in the summary.

Writes the records each study reports, one per (first reported accession, reported position, unmarked peptide or
window, degree), with the position the source gives on that accession. Nothing is placed on a reference sequence here:
scripts/02_remap_arg_methyl_to_canonical.py does that from the accession, position, marked peptide and window. Source
rows that share a record are merged: cells, enrichments and tables are listed, the best localisation probability is
kept and the rows are counted. Lysine and other non-arginine marks are dropped, and so is an arginine a row lists as
both mono- and dimethylated. A site row with no usable position is dropped. The summary printed at the end counts the
rows read and what each filter dropped.

Hartel et al. 2019, Mol Cell Proteomics 18:2149, doi:10.1074/mcp.RA119.001625 (CC BY 4.0 from 2020-09-15 per Crossref,
exclusive licence to ASBMB per PMC; PRIDE PXD012357, CC0). Tables S1-S15: LN229 S1-S2 and HEK293T S3-S15, PSMs in the
identification tables S1-S7, peptide groups in the quantification tables S8-S15. Kept: ptmRS localisation >= 75 % and
Percolator q <= 0.01, the methods' thresholds (every row of the tables meets both). The position is the one the table
gives on the first master protein ('Gene Symbol Modification', paired in order with the methyl marks; in a
quantification table 'Protein Accession and Methyl Site' when the former numbers another isoform), checked against the
peptide start the quantification tables give ('Peptide Position in Master Proteins'). S1 and S2 give no positions:
their rows take the position S3-S15 give for the same first master protein, peptide and residue, and are dropped when
no other table gives one. rme2_type is the neutral-loss call for a dimethylarginine in the tables searched for neutral
losses (S3, S5, S6, S9, S11, S12, S15), per the paper's rules:
  ADMA / SDMA   the authors localised that loss to the residue ('Localized Modified R', 'Localized Neutral Loss'), or
                the peptide has one dimethylarginine and its class is ADMA (or SDMA with no MMA on the peptide, as an
                SDMA loss next to an MMA counts only when localised; Fig. S3 legend);
  ADMA; SDMA    both losses localised to the residue, a single-dimethylarginine peptide of class Mixed, a call on a
                quantification row whose ADMA and SDMA PSMs' mean Andromeda scores differ by < 30 (ambiguous; Fig. 4
                legend), or rows of one record that call the two types (a quantification row's call, which applies
                that rule to the PSMs of its enrichment, replaces their own calls);
  empty         no loss ion; a loss not localised on a peptide with several dimethylarginines; monomethylarginines;
                the tables without the search.
The note gives the PRMT1 knockdown result of each quantification row (log2 shPRMT1/shControl, q < 0.05 or ns), and
says when every HEK293T PSM of the record comes from shPRMT1 cells.

Geoghegan et al. 2015, Nat Commun 6:6758, doi:10.1038/ncomms7758 (CC BY 4.0). Supplementary Table 1: methyl-SILAC
paired PSMs from Jurkat cells and PBMCs, anti-MMA antibodies D5A12 (CST #8015) and MeR4100 (CST #8711), light and
heavy methyl masses, and positions in 'Site description'. The positions number the first accession, in alphabetical
order, of 'Protein' and 'Alternative.Proteins' (for every row the accession Supplementary Table 2 lists the site
under), not always 'Protein'; that accession leads reported_accession, followed by the others as written. Kept: the
methods' q.value <= 0.05 (5 % peptide FDR; with the methyl-SILAC pair the set's FDR is 0.93 %) and IP.Prob > 0.25
(every row meets them), and the localisation threshold ModLS > 8 or 'Unambiguous'. The methods also accept rows at
ModLS 8 or below whose every arginine is methylated, whose methylated arginines trypsin missed, or whose PTMscore is
>= 0.95; those rows do not fix which arginine carries which methyl group, and are dropped. The note gives the record's
best ModLS. localisation_probability is left empty: PTM.Localization.Probs gives 0 to many marked arginines that ModLS
localises (NRPAIAR(+28.03)GAAGGGGR(+14.02), ModLS 138, has 0 at all three arginines). Rows whose
protein list holds a decoy (DECOY1_, DECOY2_) are dropped; their site description is UNMAPPED or doubtful.

Sylvestersen et al. 2014, Mol Cell Proteomics 13:2072, doi:10.1074/mcp.O113.032748 (CC BY 4.0 from 2020-10-12 per
Crossref; PRIDE PXD000559, CC0). mmc1 file -4 (sheets 'All identified MMA sites' and 'Triplicate experiment - MMA') and
file -6 ('Dimethyl (KR)Sites'): MaxQuant site tables from HEK293T cells enriched with the RMe and RMeGG mono-methyl
arginine antibodies. Kept: localisation probability >= 0.75 (MaxQuant class I; the tables, like the abstract's 1,027
MMA sites, also list sites below it). The 13-residue window is kept as given; the peptide is the 'Modified sequence'
with the arginine whose flanks match the window marked (window only when no single arginine matches). The export
writes every L as I in both; they are kept as given.

Lu et al. 2025, EMBO Rep, doi:10.1038/s44319-025-00599-y (CC BY 4.0). Dataset EV1 sheet 'Arginine dimethylation
sites': the 743 dimethylarginine sites (accession and position, no peptide) from HeLa peptides enriched with the SMN
Tudor domain (Proteome Discoverer, 1 % FDR). The tables give no localisation score (SEQUEST HT and Percolator), and
most of the peptides hold several arginines, so a site is kept only where a PSM localises it by mass: a PSM of the
three replicates' PSM sheets (unambiguous, Percolator q <= 0.01) of the same protein, whose every arginine is
dimethylated, with a dimethylated R at the site's position in the reference sequence (--reference); the record takes
the peptide with most such PSMs. The other sites are dropped. EV1 has no list of monomethylarginine sites, only PSMs
without protein positions, so only the dimethyl list is staged. Sites of the 'Validated peptides' sheet (hmSILAC) are
noted.

Ma et al. 2023, J Am Soc Mass Spectrom, doi:10.1021/jasms.3c00154 (ACS; PRIDE PXD038203 tables CC0). PRIDE Table S1
sheet 'Identified ADMA sites': 403 MaxQuant sites from ADMA-antibody (PTMScan adme-R) enrichments of five breast cancer
patient-derived xenografts grown in mice, with the 31-residue window. Kept: localisation probability >= 0.75 (the
paper's threshold, with 1 % FDR; every row meets it). Five rows have their protein lists split over several cells;
they are read from the cells' content, not their header. The sheet names proteins by UniProt entry name (ROA1, ARI1B),
not gene symbol, so reported_gene is left empty.

Output columns: study_id, reported_accession, reported_gene, reported_position, peptide (this record's arginine marked
R(me1) or R(me2)), window, methyl_degree (Rme1, Rme2), rme2_type, cell_type, enrichment, localisation_probability (0-1),
n_source_rows, source_tables, note. Lists are joined with '; '.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from openpyxl import load_workbook

from common import mass_localised

COLUMNS = ['study_id', 'reported_accession', 'reported_gene', 'reported_position', 'peptide', 'window', 'methyl_degree',
           'rme2_type', 'cell_type', 'enrichment', 'localisation_probability', 'n_source_rows', 'source_tables', 'note',
           'localised_by']
# How a record localises its methyl group to the one arginine (every staged record does, or it is dropped):
#   localisation_score          the study's localisation score passes a standard threshold (MaxQuant localisation
#                               probability >= 0.75, ptmRS >= 75 %, Ascore/ModScore > 13, Geoghegan's ModLS > 8)
#   one_methylatable_residue    the peptide holds no other residue the search could methylate
#   every_arginine_dimethylated every arginine of the peptide carries a dimethyl group and the search methylated
#                               arginines only, so the mass allows no other arrangement
#   curated_experimental        UniProt annotates the residue with experimental evidence (ECO:0000269)
LOCALISED_BY = ('localisation_score', 'one_methylatable_residue', 'every_arginine_dimethylated', 'curated_experimental')
ACCESSION = re.compile(r'([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(-\d+)?')
DEGREE = {'Methyl': 'Rme1', 'Dimethyl': 'Rme2'}
MARK = {'Rme1': 'R(me1)', 'Rme2': 'R(me2)', '': 'R(methyl)'}  # '' : the degree is not given


def sheet_rows(path: Path, sheet: str, header_row: int) -> tuple[dict[str, int], list[tuple]]:
    """Header name -> column index, and the non-empty data rows below the header (1-based header_row)."""
    wb = load_workbook(path, read_only=True, data_only=True)
    rows = list(wb[sheet].iter_rows(values_only=True))
    header = {str(h).strip(): k for k, h in enumerate(rows[header_row - 1]) if h is not None}
    data = [r for r in rows[header_row:] if any(v not in (None, '') for v in r)]
    return header, data


def text(value) -> str:
    return '' if value is None else str(value).strip()


def accessions(value, sep: str = ';') -> list[str]:
    return [a.strip() for a in text(value).split(sep) if a.strip()]


def read_fasta(path: Path) -> dict[str, str]:
    """Sequences of a UniProt FASTA by accession."""
    seqs, acc = {}, ''
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            acc = line.split('|')[1] if '|' in line else line[1:].split()[0]
            seqs[acc] = []
        elif acc:
            seqs[acc].append(line.strip())
    return {a: ''.join(v) for a, v in seqs.items()}


def site_row(accs, gene, position, degree, cell, enrichment, table, peptide='', offset=None, window='', loc=None,
             rme2=None, note=(), localised='localisation_score') -> dict:
    return {'accessions': list(accs), 'gene': gene, 'position': position, 'degree': degree, 'cell': cell,
            'enrichment': enrichment, 'table': table, 'peptide': peptide, 'offset': offset, 'window': window,
            'loc': loc, 'rme2': rme2 or set(), 'note': list(note), 'localised': localised}


def open_workbook(path: Path):
    """A workbook read-only; some supplements carry custom document properties openpyxl cannot read (Lim 2020), which
    are skipped."""
    try:
        return load_workbook(path, read_only=True, data_only=True)
    except TypeError:
        import io
        import zipfile
        data = io.BytesIO()
        with zipfile.ZipFile(path) as zin, zipfile.ZipFile(data, 'w') as zout:
            for item in zin.infolist():
                if item.filename != 'docProps/custom.xml':
                    zout.writestr(item, zin.read(item.filename))
        data.seek(0)
        return load_workbook(data, read_only=True, data_only=True)


def gene_index(entries: Path, history: Path) -> dict[str, list[str]]:
    """Gene symbol (upper case) -> reviewed accessions: the primary gene names of the reference entries, and for a symbol
    no entry has as its primary name, the entries that list it as a synonym."""
    primary, synonym = defaultdict(set), defaultdict(set)
    with entries.open(encoding='utf-8') as fh:
        rows = csv.DictReader(fh, delimiter='\t')
        for r in rows:
            for g in re.split(r'[;\s]+', r.get('Gene Names (primary)') or ''):
                if g:
                    primary[g.upper()].add(r['Entry'])
    with history.open(encoding='utf-8') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['kind'] == 'synonym' and r['key']:
                synonym[r['key'].upper()].add(r['accession'])
    return {g: sorted(primary.get(g) or synonym.get(g, set())) for g in set(primary) | set(synonym)}


# --------------------------------------------------------------------------------------------------- Hartel 2019
HARTEL_SHEETS = [  # sheet, label, kind, enrichment, cell, searched for neutral losses
    ('LN229 SCX IDs', 'S1', 'id', 'SCX', 'LN229', False),
    ('LN229 MMA IAP IDs', 'S2', 'id', 'MMA IAP', 'LN229', False),
    ('SCX ID', 'S3', 'id', 'SCX', 'HEK293T', True),
    ('MMA IAP ID', 'S4', 'id', 'MMA IAP', 'HEK293T', False),
    ('ADMA IAP ID', 'S5', 'id', 'ADMA IAP', 'HEK293T', True),
    ('SDMA IAP ID', 'S6', 'id', 'SDMA IAP', 'HEK293T', True),
    ('PanK IAP ID', 'S7', 'id', 'PanK IAP', 'HEK293T', False),
    ('SCX Quant MMA', 'S8', 'quant', 'SCX', 'HEK293T', False),
    ('SCX Quant DMA', 'S9', 'quant', 'SCX', 'HEK293T', True),
    ('MMA IAP Quant', 'S10', 'quant', 'MMA IAP', 'HEK293T', False),
    ('ADMA IAP Quant', 'S11', 'quant', 'ADMA IAP', 'HEK293T', True),
    ('SDMA IAP Quant', 'S12', 'quant', 'SDMA IAP', 'HEK293T', True),
    ('SCX Quant K', 'S13', 'quant', 'SCX', 'HEK293T', False),
    ('PanK IAP Quant', 'S14', 'quant', 'PanK IAP', 'HEK293T', False),
    ('SCX Mixed Quant', 'S15', 'quant', 'SCX', 'HEK293T', True),
]
HARTEL_MIN_LOC, HARTEL_MAX_Q, ANDROMEDA_AMBIGUOUS = 75.0, 0.01, 30.0


def hartel_probabilities(value) -> list[float]:
    """'100,99.99' -> [100, 99.99]; an integer cell is Excel's reading of '100,100' as 100100 and is split in threes."""
    if value is None or value == '':
        return []
    if isinstance(value, int):
        s = str(value)
        head = len(s) % 3 or 3
        return [float(s[:head])] + [float(s[i:i + 3]) for i in range(head, len(s), 3)]
    if isinstance(value, float):
        return [value]
    return [float(p) for p in re.split('[,;]', value) if p.strip()]


def hartel_id_mods(value) -> list[tuple[str, int, str]]:
    """'R4(Methyl),R16(Dimethyl),M7(Oxidation)' -> [('R', 4, 'Methyl'), ('R', 16, 'Dimethyl')] (literal repeats once)."""
    out = []
    for m in re.finditer(r'([A-Z])(\d+)\((Methyl|Dimethyl|Trimethyl)\)', text(value)):
        mod = (m.group(1), int(m.group(2)), m.group(3))
        if mod not in out:
            out.append(mod)
    return out


def hartel_quant_mods(value) -> tuple[list[tuple[str, int, str, float | None]], bool]:
    """'2xDimethyl [R10(100);R12(99.2)];1xMethyl [R3(100)]' -> methyl marks with their probabilities, and whether
    each group holds the number of marks it states. A mark after a group's bracket belongs to that group
    ('3xDimethyl [R6(100);R10(100)];R8(100)')."""
    out, counts, kind = [], Counter(), None
    stated = {}
    for m in re.finditer(r'(\d+)x([A-Za-z]+) \[([^\]]*)\]|([A-Z])(\d+)(?:\(([\d.]+)\))?', text(value)):
        if m.group(1):
            kind = m.group(2)
            stated[kind] = stated.get(kind, 0) + int(m.group(1))
            items = re.findall(r'([A-Z])(\d+)(?:\(([\d.]+)\))?', m.group(3))
        else:
            items = [(m.group(4), m.group(5), m.group(6))]
        for res, off, prob in items:
            counts[kind] += 1
            if kind in ('Methyl', 'Dimethyl', 'Trimethyl'):
                out.append((res, int(off), kind, float(prob) if prob else None))
    return out, all(counts[k] == n for k, n in stated.items() if k in ('Methyl', 'Dimethyl', 'Trimethyl'))


def hartel_protein_sites(value) -> dict[str, list[int]]:
    """'P61978-2 1xMethyl [R256];P61978-3 1xMethyl [R232]' -> {'P61978-2': [256], 'P61978-3': [232]}."""
    out, current = defaultdict(list), None
    for m in re.finditer(r'(?:([A-Z0-9-]+) )?\d+x[A-Za-z]+ \[([^\]]*)\]', text(value)):
        current = m.group(1) or current
        out[current] += [int(p) for p in re.findall(r'[A-Z](\d+)', m.group(2))]
    return out


def hartel_localised(value) -> dict[int, set[str]]:
    """'R8 ADMA-y6*;R8 ADMA-y7*' or 'R12 ADMA, R14 ADMA' -> {peptide residue: {'ADMA'}}."""
    out = defaultdict(set)
    for m in re.finditer(r'R(\d+) (ADMA|SDMA)', text(value)):
        out[int(m.group(1))].add(m.group(2))
    return out


def hartel_nl_call(offset: int, localised: dict, nl_class: str, n_me2: int, has_mma: bool, andromeda) -> set[str]:
    calls = localised.get(offset, set())
    cl = nl_class.lower()
    if calls:
        call = set(calls)
    elif cl and n_me2 == 1 and cl in ('adma', '2xadma'):
        call = {'ADMA'}
    elif cl and n_me2 == 1 and cl == 'sdma' and not has_mma:
        call = {'SDMA'}
    elif cl and n_me2 == 1 and cl == 'mixed':
        call = {'ADMA', 'SDMA'}
    else:
        call = set()
    if call and isinstance(andromeda, (int, float)) and andromeda < ANDROMEDA_AMBIGUOUS:
        call = {'ADMA', 'SDMA'}
    return call


def condition(raw: str) -> str:
    if re.search(r'-(Cont|C1|C2)-', raw):
        return 'shControl'
    if re.search(r'-(sh|sh1|sh2)-', raw):
        return 'shPRMT1'
    return ''


def read_hartel(paths: list[Path], stats: Counter, ctx: dict) -> list[dict]:
    books = {}
    for path in paths:
        for name in load_workbook(path, read_only=True).sheetnames:
            books[name] = path
    out, borrowed = [], []
    lookup = defaultdict(set)  # (first master, peptide, residue) -> positions S3-S15 give
    for sheet, label, kind, enrichment, cell, nl_tested in HARTEL_SHEETS:
        if sheet not in books:
            raise SystemExit(f'Hartel sheet {sheet!r} ({label}) not in {[str(p) for p in paths]}')
        H, data = sheet_rows(books[sheet], sheet, 3)
        stats[f'rows_read:{label}'] += len(data)
        for r in data:
            if 'Methyl or Nonmethyl Peptide' in H and text(r[H['Methyl or Nonmethyl Peptide']]) != 'Methyl':
                stats['rows_skipped:nonmethyl_peptide'] += 1
                continue
            if float(r[H['Percolator.q.Value']]) > HARTEL_MAX_Q:
                stats['rows_dropped:q_above_0.01'] += 1
                continue
            peptide = text(r[H['Peptide Sequence']]).upper()
            masters = accessions(r[H['Master Protein Accessions']])
            if not masters or not all(ACCESSION.fullmatch(a) for a in masters):
                stats['rows_dropped:not_uniprot_accession'] += 1
                continue
            gene = text(r[H['Gene Symbol']]).split(';')[0].strip() if 'Gene Symbol' in H else ''
            note, starts = [], []
            if kind == 'id':
                mods = [(res, off, k, None) for res, off, k in hartel_id_mods(r[H['Modifications']])]
                residues = list(dict.fromkeys((res, off) for res, off, _, _ in mods))
                probs = hartel_probabilities(r[H['Localization Probability']])
                if len(probs) == len(mods):
                    mods = [(res, off, k, p) for (res, off, k, _), p in zip(mods, probs)]
                elif len(probs) == len(residues):
                    by_res = dict(zip(residues, probs))
                    mods = [(res, off, k, by_res[(res, off)]) for res, off, k, _ in mods]
                positions = []
                if 'Gene Symbol Modification' in H:
                    m = re.search(r'\s([\d,]+)\s*$', ' ' + text(r[H['Gene Symbol Modification']]))
                    positions = [int(x) for x in m.group(1).split(',') if x] if m else []
                if positions and len(positions) == len(mods):
                    pos_of = {(res, off): p for (res, off, _, _), p in zip(mods, positions)}
                elif positions and len(positions) == len(residues):
                    pos_of = dict(zip(residues, positions))
                else:
                    pos_of = {}
                    if 'Gene Symbol Modification' in H and mods:
                        stats['rows_position_list_mismatch'] += 1
                localised = hartel_localised(r[H['Localized Modified R']]) if nl_tested else {}
                nl_class = text(r[H['Neutral Loss Class']]) if nl_tested else ''
                andromeda = None
                if 'Raw File Origin' in H and cell == 'HEK293T':
                    note.append(condition(text(r[H['Raw File Origin']])))
            else:
                mods, complete = hartel_quant_mods(r[H['Modifications']])
                if not complete:
                    stats['rows_mod_count_mismatch'] += 1
                ppm = text(r[H['Peptide Position in Master Proteins']])
                current = None
                for token in ppm.split(';'):
                    m = re.fullmatch(r'\s*(\S+)?\s*\[(\d+)-(\d+)\]\s*', token)
                    if m:
                        current = m.group(1) or current
                        if current == masters[0]:
                            starts.append(int(m.group(2)))
                if 'Gene Symbol Modification' in H:
                    listed = [int(p) for p in re.findall(r'[A-Z](\d+)', text(r[H['Gene Symbol Modification']]).split(' ')[-1])]
                else:  # S13, S14 give 'O43929 1xMethyl [K7];P12345 ...' per master protein
                    first = re.split(r';(?=[A-Z0-9-]+ \d+x)', text(r[H['Protein Accession and Methyl Site']]))[0]
                    listed = sorted({int(p) for p in re.findall(r'[A-Z](\d+)', first.split(' ', 1)[-1])})
                ordered = sorted({(res, off) for res, off, _, _ in mods}, key=lambda m: m[1])
                pos_of = dict(zip(ordered, listed)) if len(listed) == len(ordered) else {}
                if not pos_of or not all(p - off + 1 in starts for (res, off), p in pos_of.items()):
                    # 'Gene Symbol Modification' can number another isoform; the per-protein column names the master
                    own = sorted(set(hartel_protein_sites(r[H['Protein Accession and Methyl Site']]).get(masters[0], [])))
                    if len(own) == len(ordered) and all(p - off + 1 in starts for (res, off), p in zip(ordered, own)):
                        pos_of = dict(zip(ordered, own))
                        stats['rows_positions_from_protein_accession_and_methyl_site'] += 1
                if mods and not pos_of:
                    stats['rows_position_list_mismatch'] += 1
                localised = hartel_localised(r[H['Localized Neutral Loss']]) if nl_tested else {}
                nl_class = text(r[H['Neutral Loss Class']]) if nl_tested else ''
                andromeda = r[H['Mean Andromeda Score difference between ADMA and SDMA PSMs']] if nl_tested else None
                fold = r[H['log2 Fold Change (shPRMT1 / shControl)']]
                if isinstance(fold, (int, float)):
                    significant = text(r[H['Significant Perseus t-test q < 0.05']]) == '+'
                    note.append(f'{label} PRMT1 KD log2 {fold:.2f} {"q<0.05" if significant else "ns"}')
            degrees = defaultdict(set)
            for res, off, k, _ in mods:
                degrees[(res, off)].add(k)
            both = {key for key, ks in degrees.items() if len(ks) > 1}
            n_me2 = sum(1 for res, off, k, _ in mods if res == 'R' and k == 'Dimethyl' and (res, off) not in both)
            has_mma = any(res == 'R' and k == 'Methyl' for res, off, k, _ in mods)
            dimethyl = {off for res, off, k, _ in mods if res == 'R' and k == 'Dimethyl'}
            if any(off not in dimethyl for off in localised):
                stats['rows_nl_localised_to_non_dimethyl_residue'] += 1
            # positions on the first master must agree with each other (ID tables) or with its peptide start
            consistent = True
            if pos_of:
                shifts = {p - off for (res, off), p in pos_of.items()}
                if kind == 'quant' and starts:
                    consistent = all(p - off + 1 in starts for (res, off), p in pos_of.items())
                elif len(shifts) > 1 or min(shifts) < 0:
                    consistent = False
            for res, off, k, prob in mods:
                if res != 'R' or k not in DEGREE:
                    stats['marks_dropped:lysine_or_trimethyl'] += 1
                    continue
                if (res, off) in both:
                    if k == 'Methyl':
                        stats['site_rows_dropped:same_arginine_me1_and_me2'] += 1
                    continue
                if not 1 <= off <= len(peptide) or peptide[off - 1] != 'R':
                    stats['site_rows_dropped:marked_residue_not_R'] += 1
                    continue
                if prob is None or prob < HARTEL_MIN_LOC:
                    stats['site_rows_dropped:localisation_below_75' if prob is not None else 'site_rows_dropped:no_localisation_probability'] += 1
                    continue
                if pos_of and not consistent:
                    stats['site_rows_dropped:position_inconsistent_with_peptide'] += 1
                    continue
                degree = DEGREE[k]
                rme2 = hartel_nl_call(off, localised, nl_class, n_me2, has_mma, andromeda) if nl_tested and degree == 'Rme2' else set()
                row = site_row(masters, gene, pos_of.get((res, off)), degree, cell, enrichment, label, peptide, off,
                               loc=prob / 100, rme2=rme2, note=[n for n in note if n])
                row['kd_condition'] = note[0] if kind == 'id' and note else ''
                row['quant'] = kind == 'quant'
                if row['position'] is None and label in ('S1', 'S2'):
                    borrowed.append(row)
                    continue
                if row['position'] is None:
                    stats['site_rows_dropped:no_position_given'] += 1
                    continue
                lookup[(masters[0], peptide, off)].add(row['position'])
                out.append(row)
    for row in borrowed:
        found = lookup.get((row['accessions'][0], row['peptide'], row['offset']), set())
        if len(found) == 1:
            row['position'] = next(iter(found))
            stats['site_rows_S1_S2_positioned_from_S3_S15'] += 1
            out.append(row)
        else:
            stats['site_rows_dropped:S1_S2_no_position_in_any_table' if not found else 'site_rows_dropped:S1_S2_several_positions'] += 1
    return out


def hartel_rme2_type(rows: list[dict]) -> set[str]:
    """The calls of a record's rows; a quantification row's call (made over the PSMs of its enrichment with the
    Andromeda-score rule) replaces the calls of that enrichment's PSMs."""
    calls = set()
    for enrichment in {r['enrichment'] for r in rows}:
        quant = set().union(*(r['rme2'] for r in rows if r['enrichment'] == enrichment and r['quant']))
        calls |= quant or set().union(*(r['rme2'] for r in rows if r['enrichment'] == enrichment and not r['quant']))
    return calls


def hartel_merge_note(rows: list[dict]) -> list[str]:
    """PRMT1 knockdown results, and whether every HEK293T PSM of the record came from shPRMT1 cells."""
    notes = [n for r in rows for n in r['note'] if n not in ('shControl', 'shPRMT1')]
    conditions = {r['kd_condition'] for r in rows if r['kd_condition']}
    if conditions == {'shPRMT1'}:
        notes.append('HEK293T PSMs from shPRMT1 cells only')
    return notes


# ------------------------------------------------------------------------------------------------ Geoghegan 2015
GEOGHEGAN_MAX_Q, GEOGHEGAN_MIN_IP, GEOGHEGAN_MIN_MODLS = 0.05, 0.25, 8
GEOGHEGAN_ANTIBODY = {'D5A12': ['anti-MMA antibody D5A12'], 'MeR4100': ['anti-MMA antibody MeR4100'],
                      'D5A12+MeR4100': ['anti-MMA antibody D5A12', 'anti-MMA antibody MeR4100']}


def geoghegan_degree(mass: float) -> str:
    """Light or heavy (13CD3) methyl: +14.02 / +18.04 mono, +28.03 / +36.08 di."""
    for degree, masses in (('Rme1', (14.02, 18.04)), ('Rme2', (28.03, 36.08))):
        if any(abs(mass - m) < 0.05 for m in masses):
            return degree
    return ''


def geoghegan_modls(modls) -> str:
    """The row's localisation as the note gives it ('ModLS unambiguous', 'ModLS 138'), empty at ModLS 8 or below."""
    if modls == 'Unambiguous':
        return 'ModLS unambiguous'
    return f'ModLS {float(modls):g}' if float(modls) > GEOGHEGAN_MIN_MODLS else ''


def geoghegan_merge_note(grp: list[dict]) -> list[str]:
    """The best localisation of the record's rows: unambiguous, else the highest ModLS."""
    def rank(note):
        score = note.split()[-1]
        return (score == 'unambiguous', 0.0 if score == 'unambiguous' else float(score))
    return [max((n for r in grp for n in r['note']), key=rank)]


def read_geoghegan(paths: list[Path], stats: Counter, ctx: dict) -> list[dict]:
    H, data = sheet_rows(paths[0], 'Supplementary Table 1', 2)
    stats['rows_read:ST1'] += len(data)
    out = []
    for r in data:
        if float(r[H['q.value']]) > GEOGHEGAN_MAX_Q:
            stats['rows_dropped:q_above_0.05'] += 1
            continue
        if float(r[H['IP.Prob']]) <= GEOGHEGAN_MIN_IP:
            stats['rows_dropped:ip_probability_0.25_or_less'] += 1
            continue
        protein = text(r[H['Protein']])
        accs = list(dict.fromkeys([protein] + accessions(r[H['Alternative.Proteins']])))
        if any(a.startswith('DECOY') for a in accs):
            stats['rows_dropped:decoy_in_protein_list'] += 1
            continue
        if not all(ACCESSION.fullmatch(a) for a in accs):
            stats['rows_dropped:not_uniprot_accession'] += 1
            continue
        numbered = min(accs)  # the protein 'Site description' numbers
        accs = [numbered] + [a for a in accs if a != numbered]
        if numbered != protein:
            stats['rows_numbered_on_alternative_protein'] += 1
        gene = re.search(r'\bGN=(\S+)', text(r[H['Protein.Description']]))
        gene = gene.group(1) if gene else ''
        residues, marks = '', {}
        for aa, mass in re.findall(r'([A-Z])(?:\(([+-][\d.]+)\))?', text(r[H['Modified.Peptide']])):
            residues += aa
            if mass:
                marks[len(residues)] = (aa, float(mass))
        peptide = text(r[H['Peptide']]).upper()
        if residues != peptide:
            stats['rows_dropped:modified_peptide_differs_from_peptide'] += 1
            continue
        sites = [(off, geoghegan_degree(mass)) for off, (aa, mass) in sorted(marks.items()) if aa == 'R' and geoghegan_degree(mass)]
        stats['marks_dropped:lysine_or_other'] += sum(1 for aa, mass in marks.values() if aa == 'K')
        description = [t for t in text(r[H['Site description']]).split(';') if t]
        reported = [re.fullmatch(r'(\d+)__(MMA|DMA)__([LH])', t) for t in description]
        if not sites:
            stats['rows_without_arginine_methyl'] += 1
            continue
        if not description or not all(reported) or len(reported) != len(sites):
            stats['site_rows_dropped:no_position_given' if not description or 'UNMAPPED' in description
                  else 'site_rows_dropped:site_description_mismatch'] += len(sites)
            continue
        if any((m.group(2) == 'MMA') != (degree == 'Rme1') for m, (off, degree) in zip(reported, sites)):
            stats['site_rows_dropped:site_description_degree_mismatch'] += len(sites)
            continue
        shifts = {int(m.group(1)) - off for m, (off, _) in zip(reported, sites)}
        if len(shifts) > 1 or min(shifts) < 0:  # e.g. rows with a methyllysine before the arginine
            stats['site_rows_dropped:position_inconsistent_with_peptide'] += len(sites)
            continue
        localised = geoghegan_modls(r[H['ModLS']])
        if not localised:
            stats['site_rows_dropped:modls_8_or_less'] += len(sites)
            continue
        cell = text(r[H['Cell']])
        enrichment = GEOGHEGAN_ANTIBODY[text(r[H['Antibody']])]
        for m, (off, degree) in zip(reported, sites):
            row = site_row(accs, gene, int(m.group(1)), degree, cell, '', 'ST1', peptide, off, note=[localised])
            row['enrichments'] = enrichment
            out.append(row)
    return out


# --------------------------------------------------------------------------------------------- Sylvestersen 2014
SYLVESTERSEN_MIN_LOC = 0.75
SYLVESTERSEN_SHEETS = [  # sheet, label, degree, position column, accession column, gene column
    ('All identified MMA sites', '032748-4 All identified MMA sites', 'Rme1', 'Positions within proteins', 'Uniprot ID', 'Gene name'),
    ('Triplicate experiment - MMA ', '032748-4 Triplicate experiment - MMA', 'Rme1', 'Positions within proteins', 'Protein', 'Gene name'),
    ('Dimethyl (KR)Sites', '032748-6 Dimethyl (KR)Sites', 'Rme2', 'Positions', 'Protein', 'Gene names'),
]
MAXQUANT_MARK = {'me': 'Rme1', 'di': 'Rme2'}


def maxquant_peptide(modified: str) -> tuple[str, dict[int, str]]:
    """'_(ac)ASEAR(me)GGIGAPPIQSAR(me)_' -> ('ASEARGGIGAPPIQSAR', {5: 'me', 17: 'me'}) (1-based residue: mark)."""
    residues, marks = '', {}
    for aa, mark in re.findall(r'([A-Z])(?:\((\w+)\))?', modified.strip('_').replace('(ac)', '')):
        residues += aa
        if mark:
            marks[len(residues)] = mark
    return residues, marks


def window_offset(peptide: str, candidates: list[int], window: str) -> int | None:
    """The one candidate residue (1-based) whose flanks in the peptide match the 13-residue window (I = L)."""
    half = len(window) // 2
    fits = []
    for off in candidates:
        ok = True
        for k in range(-half, half + 1):
            j = off - 1 + k
            if 0 <= j < len(peptide):
                w = window[half + k]
                if w == '_' or w.replace('I', 'L') != peptide[j].replace('I', 'L'):
                    ok = False
                    break
        if ok:
            fits.append(off)
    return fits[0] if len(fits) == 1 else None


def read_sylvestersen(paths: list[Path], stats: Counter, ctx: dict) -> list[dict]:
    books = {}
    for path in paths:
        for name in load_workbook(path, read_only=True).sheetnames:
            books[name] = path
    out = []
    for sheet, label, degree, pos_col, acc_col, gene_col in SYLVESTERSEN_SHEETS:
        if sheet not in books:
            raise SystemExit(f'Sylvestersen sheet {sheet!r} not in {[str(p) for p in paths]}')
        H, data = sheet_rows(books[sheet], sheet, 1)
        stats[f'rows_read:{label}'] += len(data)
        for r in data:
            m = re.fullmatch(r'([A-Z])(\d+)', text(r[H[pos_col]]))
            accs = accessions(r[H[acc_col]])
            if not m:
                stats['rows_dropped:no_position_given'] += 1
                continue
            if m.group(1) != 'R':
                stats['rows_dropped:lysine'] += 1
                continue
            if not accs or not all(ACCESSION.fullmatch(a) for a in accs):
                stats['rows_dropped:decoy_contaminant_or_not_uniprot'] += 1
                continue
            loc = float(r[H['Localization prob']])
            if loc < SYLVESTERSEN_MIN_LOC:
                stats['rows_dropped:localisation_below_0.75'] += 1
                continue
            window = text(r[H['Sequence window']])
            peptide, marks = maxquant_peptide(text(r[H['Modified sequence']]))
            candidates = [off for off, mk in marks.items() if peptide[off - 1] == 'R' and MAXQUANT_MARK.get(mk) == degree]
            off = window_offset(peptide, candidates, window) if window[len(window) // 2:len(window) // 2 + 1] == 'R' else None
            if off is None:
                stats['rows_window_only:site_not_matched_in_peptide'] += 1
            gene = text(r[H[gene_col]]).split(';')[0].strip()
            out.append(site_row(accs, gene, int(m.group(2)), degree, 'HEK293T', 'RMe and RMeGG MMA antibodies', label,
                                peptide if off else '', off, window=window, loc=loc))
    return out


# ------------------------------------------------------------------------------------------------------ Lu 2025
LU_MAX_Q = 0.01


def lu_psms(wb, stats: Counter) -> dict[str, Counter]:
    """Dimethylarginine PSMs of the three replicates: {accession: Counter((peptide, 1-based offset of the dimethyl R))},
    for unambiguous PSMs at Percolator q <= 0.01 (the methods' 1 % FDR) whose every arginine is dimethylated, under each
    accession of their protein group. The tables give no localisation score; in such a PSM the mass alone puts a
    dimethyl group on each arginine."""
    out = defaultdict(Counter)
    for name in ('PSMs of replicate1', 'PSMs of replicate2', 'PSMs of replicate3'):
        rows = list(wb[name].iter_rows(values_only=True))
        head = {h: k for k, h in enumerate(rows[1])}
        for r in rows[2:]:
            if text(r[head['PSM Ambiguity']]) != 'Unambiguous' or r[head['Percolator q-Value']] is None or float(r[head['Percolator q-Value']]) > LU_MAX_Q:
                stats['psm_rows_dropped:ambiguous_or_q'] += 1
                continue
            peptide = text(r[head['Annotated Sequence']]).upper()
            offsets = [int(x) for x in re.findall(r'R(\d+)\(Dimethyl\)', text(r[head['Modifications']]))]
            if not offsets:
                continue
            stats['psm_rows_read:dimethyl'] += 1
            if len(offsets) != peptide.count('R'):
                stats['psm_rows_dropped:not_every_arginine_dimethylated'] += 1
                continue
            for acc in accessions(r[head['Protein Accessions']]):
                for off in offsets:
                    out[acc.split('-')[0]][(peptide, off)] += 1
    return out


def read_lu(paths: list[Path], stats: Counter, ctx: dict) -> list[dict]:
    reference = ctx['reference']
    wb = load_workbook(paths[0], read_only=True, data_only=True)
    psms = lu_psms(wb, stats)
    rows = [text(r[0]) for r in wb['Arginine dimethylation sites'].iter_rows(values_only=True)]
    title, sites = rows[0], [v for v in rows[1:] if v]
    stats['rows_read:EV1 Arginine dimethylation sites'] += len(sites)
    stated = re.search(r'\((\d+)\)', title)
    if stated and int(stated.group(1)) != len(sites):
        raise SystemExit(f'Lu 2025: the sheet title states {stated.group(1)} sites, it lists {len(sites)}')
    validated = set()
    for r in list(wb['Validated peptides'].iter_rows(values_only=True))[2:]:
        m = re.match(r'(\S+) \d+x\w+ \[([^\]]*)\]', text(r[2]))
        if m:
            validated |= {(m.group(1), int(p)) for p in re.findall(r'R(\d+)', m.group(2))}
    out = []
    for v in sites:
        m = re.fullmatch(r'\(([^)]+)\)-R(\d+)', v)
        if not m:
            stats['rows_dropped:no_position_given'] += 1
            continue
        acc, pos = m.group(1), int(m.group(2))
        if not ACCESSION.fullmatch(acc):
            stats['rows_dropped:not_uniprot_accession'] += 1
            continue
        note = ['hmSILAC validated peptide (EV1)'] if (acc, pos) in validated else []
        # The site's peptide: a localised dimethylarginine PSM of the same protein whose marked R lies at this position
        # of the reference sequence (most PSMs, then the shorter peptide). Without one the site is dropped.
        seq, found = reference.get(acc.split('-')[0], ''), []
        for (peptide, off), n in psms.get(acc.split('-')[0], {}).items():
            start = pos - off
            if start >= 0 and seq[start:start + len(peptide)] == peptide:
                found.append((-n, len(peptide), peptide, off))
        if found:
            _, _, peptide, off = sorted(found)[0]
            stats['sites_with_psm_peptide'] += 1
            note.append(f'{-sorted(found)[0][0]} PSMs with every arginine dimethylated')
            localised = mass_localised(peptide, {k: 'Rme2' for k, aa in enumerate(peptide, 1) if aa == 'R'}, 'R')
            out.append(site_row([acc], '', pos, 'Rme2', 'HeLa', 'SMN Tudor', 'EV1 Arginine dimethylation sites; EV1 PSMs',
                                peptide=peptide, offset=off, note=note, localised=localised))
        else:
            stats['sites_dropped:no_localised_psm'] += 1
    return out


# ------------------------------------------------------------------------------------------------------ Ma 2023
MA_MIN_LOC = 0.75


def read_ma(paths: list[Path], stats: Counter, ctx: dict) -> list[dict]:
    H, data = sheet_rows(paths[0], 'Identified ADMA sites', 1)
    stats['rows_read:Table S1 Identified ADMA sites'] += len(data)
    first_t = H['T: T: Proteins']
    out = []
    for r in data:
        if text(r[H['C: C: Amino acid']]) != 'R':
            stats['rows_dropped:not_arginine'] += 1
            continue
        loc, position = float(r[H['N: N: Localization prob']]), int(r[H['N: N: Position']])
        # the T: columns hold proteins, positions, leading proteins, protein and window, whatever cells they spread over
        cells = [text(v) for v in r[first_t:] if text(v)]
        k = 0
        proteins = []
        while k < len(cells) and '|' in cells[k]:
            proteins.append(cells[k]); k += 1
        while k < len(cells) and re.fullmatch(r'[\d;]+', cells[k]):
            k += 1
        named = []
        while k < len(cells) and '|' in cells[k]:
            named.append(cells[k]); k += 1
        window = cells[k] if k < len(cells) else ''
        if len(named) < 2 or not re.fullmatch(r'[A-Z_]+', window):
            stats['rows_dropped:unreadable_protein_columns'] += 1
            continue
        protein, leading = named[-1], [a for cell in named[:-1] for a in cell.split(';')]
        if any(a.startswith(('REV__', 'CON__')) for a in proteins + leading + [protein]):
            stats['rows_dropped:decoy_or_contaminant'] += 1
            continue
        acc_of = lambda s: s.split('|')[1] if s.count('|') >= 2 else s
        accs = list(dict.fromkeys([acc_of(protein)] + [acc_of(a) for a in leading]))
        if not all(ACCESSION.fullmatch(a) for a in accs):
            stats['rows_dropped:not_uniprot_accession'] += 1
            continue
        if loc < MA_MIN_LOC:
            stats['rows_dropped:localisation_below_0.75'] += 1
            continue
        if not text(r[H['Sites']]).endswith(f'_{position}'):
            stats['rows_site_label_differs_from_position'] += 1
        note = ['patient-derived xenograft in mouse']
        out.append(site_row(accs, '', position, 'Rme2', 'breast cancer PDX',
                            'ADMA antibody (PTMScan adme-R)', 'Table S1 Identified ADMA sites', window=window, loc=loc,
                            note=note))
    return out


# ------------------------------------------------------------------------------------------------------ Guo 2014
GUO_MIN_MODSCORE, GUO_MAX_PPM = 13.0, 10.0
GUO_ANTIBODY = {'Me-R4-100': 'MMA antibody Me-R4-100', 'D5A12': 'MMA antibody R*GG (D5A12)', 'D4H5': 'ADMA antibody D4H5',
                'D6A8': 'ADMA antibody D6A8'}
MONO = {'G': 57.02146, 'A': 71.03711, 'S': 87.03203, 'P': 97.05276, 'V': 99.06841, 'T': 101.04768, 'C': 160.03065,
        'L': 113.08406, 'I': 113.08406, 'N': 114.04293, 'D': 115.02694, 'Q': 128.05858, 'K': 128.09496, 'E': 129.04259,
        'M': 131.04049, 'H': 137.05891, 'F': 147.06841, 'R': 156.10111, 'Y': 163.06333, 'W': 186.07931}  # C carbamidomethyl


def guo_degree(peptide: str, n_marks: int, n_ox: int, mz: float, charge: int) -> str:
    """'Rme1' or 'Rme2' when the precursor mass puts one or two methyl groups on every marked residue, else ''."""
    observed = mz * charge - charge * 1.007276
    units = (observed - sum(MONO[a] for a in peptide) - 18.010565 - n_ox * 15.99491) / 14.01565
    n = round(units)
    if abs(units - n) * 14.01565 > GUO_MAX_PPM * observed / 1e6:
        return ''
    return 'Rme1' if n == n_marks else 'Rme2' if n == 2 * n_marks else ''


def guo_peptide(text_value: str) -> tuple[str, list[int], int]:
    """'K.SYSGGTR*TSSGR.L' -> ('SYSGGTRTSSGR', [7], 0): residues, 1-based methyl marks (R* and K*), oxidised M (M#)."""
    core = text_value.split('.')[1] if text_value.count('.') >= 2 else text_value
    residues, marks, ox = '', [], 0
    for ch in core:
        if 'A' <= ch <= 'Z':
            residues += ch
        elif ch == '*':
            marks.append(len(residues))
        elif ch == '#':
            ox += 1
    return residues, marks, ox


def read_guo(paths: list[Path], stats: Counter, ctx: dict) -> list[dict]:
    """Table S3 lists the HCT116 sites (13-residue window, lower-case r at the site, RefSeq numbering); Table S2 the
    peptide-spectrum matches with the ModScore of each site. A site is kept where an HCT116 PSM of its gene, its
    peptide aligned with the site's window, localises the methyl group there at ModScore > 13 and its precursor mass
    gives the degree. The gene symbol names the reviewed entries; the RefSeq accession is kept in the note."""
    s2_path, s3_path = sorted(paths, key=lambda p: p.name)  # ...-3.xlsx is Table S2, ...-4.xlsx Table S3
    H2, rows2 = sheet_rows(s2_path, open_workbook(s2_path).sheetnames[0], 5)
    stats['rows_read:Table S2'] += len(rows2)
    psms = defaultdict(list)  # gene -> (peptide, offset, ModScore, degree, antibody)
    for r in rows2:
        if text(r[H2['sample']]) != 'HCT116':
            continue
        stats['psm_rows_read:HCT116'] += 1
        peptide, marks, ox = guo_peptide(text(r[H2['modScore_peptide']]) or text(r[H2['peptide']]))
        degree = guo_degree(peptide, len(marks), ox, float(r[H2['pre_m_z']]), int(float(r[H2['charge']])))
        if not degree:
            stats['psm_rows_dropped:degree_not_fixed_by_mass'] += 1
            continue
        scores = {}
        for k in range(1, 5):
            pos, score = text(r[H2[f'site_{k}_position']]), text(r[H2[f'site_{k}_score']])
            if pos and score:
                scores[int(float(pos))] = float(score)
        for off in marks:
            if peptide[off - 1] != 'R':
                continue
            if scores.get(off, 0) <= GUO_MIN_MODSCORE:
                stats['psm_sites_dropped:modscore_13_or_less'] += 1
                continue
            psms[text(r[H2['gene_symbol']]).upper()].append((peptide, off, scores[off], degree, text(r[H2['antibody ID']])))
    H3, rows3 = sheet_rows(s3_path, open_workbook(s3_path).sheetnames[0], 5)
    stats['rows_read:Table S3'] += len(rows3)
    out = []
    for r in rows3:
        gene, window = text(r[H3['Gene name']]), text(r[H3['-6/+6 peptide (lower case r is the methylated residue)']])
        position = text(r[H3['Site of Methylation']]).lstrip('§').strip()
        if not gene or not re.fullmatch(r'\d+', position) or len(window) != 13 or window[6] != 'r':
            stats['rows_dropped:no_gene_position_or_window'] += 1
            continue
        if text(r[H3['Site of Methylation']]).startswith('§'):
            stats['rows_previously_published_site'] += 1
        accs = ctx['genes'].get(gene.upper(), [])
        if not accs:
            stats['rows_dropped:gene_symbol_not_in_reference'] += 1
            continue
        fits = []
        for peptide, off, score, degree, antibody in psms.get(gene.upper(), []):
            if all(window[6 + k] in 'xX' or not 0 <= off - 1 + k < len(peptide) or peptide[off - 1 + k] == window[6 + k].upper()
                   for k in range(-6, 7)):
                fits.append((peptide, off, score, degree, antibody))
        if not fits:
            stats['rows_dropped:no_psm_localised_at_modscore_above_13'] += 1
            continue
        refseq = text(r[H3['Accession #']])
        for peptide, off, score, degree, antibody in fits:
            row = site_row(accs, gene, int(position), degree, 'HCT116', GUO_ANTIBODY.get(antibody, antibody), 'Table S2; Table S3',
                           peptide, off, window=window.upper().replace('X', '_'), note=[f'ModScore {score:g}', f'RefSeq {refseq}'])
            row['score'] = score
            out.append(row)
    return out


# ----------------------------------------------------------------------------------------------------- Fong 2019
FONG_MIN_LOC, FONG_MIN_SCORE = 0.75, 25.0
FONG_ENRICHMENT = {'mono': 'MMA antibody (PTMScan mme-RG)', 'sym': 'SDMA antibody (PTMScan sdme-R)',
                   'asym': 'ADMA antibody (PTMScan adme-R)'}


def read_fong(paths: list[Path], stats: Counter, ctx: dict) -> list[dict]:
    """Table S4 ('Redundant peptide list'): MaxQuant evidence rows of the NB4 SILAC experiments, one per methyl mark,
    with the leading protein, the position, LcPrb (localisation probability) and Andromeda score. Kept: arginine marks
    at LcPrb >= 0.75 (the paper keeps >= 0.5) and score >= 25 (the paper's)."""
    H, data = sheet_rows(paths[0], 'Redundant peptide list', 7)
    stats['rows_read:Table S4'] += len(data)
    out = []
    for r in data:
        if text(r[H['Residue']]) != 'R':
            stats['rows_dropped:not_arginine'] += 1
            continue
        accs = list(dict.fromkeys(a for a in accessions(r[H['Leading Proteins']]) + accessions(r[H['All Proteins']])
                                  if not a.startswith(('REV__', 'CON__'))))
        if not accs or not all(ACCESSION.fullmatch(a) for a in accs):
            stats['rows_dropped:decoy_contaminant_or_not_uniprot'] += 1
            continue
        loc, score = float(r[H['LcPrb']]), float(r[H['Score']])
        if loc < FONG_MIN_LOC or score < FONG_MIN_SCORE:
            stats['rows_dropped:localisation_below_0.75_or_score_below_25'] += 1
            continue
        degree = MAXQUANT_MARK.get(text(r[H['Modification']]))
        peptide, marks = maxquant_peptide(text(r[H['Modified sequence']]))
        window = text(r[H['Modification Window']])
        candidates = [off for off, mk in marks.items() if peptide[off - 1] == 'R' and MAXQUANT_MARK.get(mk) == degree]
        # the window picks the mark; 258 windows are a run of z and 15 empty, where the peptide's only such mark does
        off = (window_offset(peptide, candidates, window) if len(window) == 15 and window[7] == 'R'
               else candidates[0] if len(candidates) == 1 else None)
        if not degree or off is None:
            stats['rows_dropped:mark_not_matched_in_peptide'] += 1
            continue
        experiment = text(r[H['Experiment']])
        enrichment = next((v for k, v in FONG_ENRICHMENT.items() if f'_{k}_' in f'_{experiment}_'), '')
        out.append(site_row(accs, text(r[H['Gene Names']]).split(';')[0], int(float(r[H['Position']])), degree, 'NB4',
                            enrichment, 'Table S4', peptide, off, window=window, loc=loc))
    return out


# ------------------------------------------------------------------------------------------------------- Li 2021
LI_SHEETS = re.compile(r'(PRMT7|SGC3027|PRMT4|PRMT5)[- ]?regulated (me1|me2)', re.I)
LI_ENRICHMENT = {'PRMT7': 'MMA antibody', 'SGC3027': 'MMA antibody', 'PRMT4': 'MMA or ADMA antibody',
                 'PRMT5': 'MMA or SDMA antibody'}


def read_li(paths: list[Path], stats: Counter, ctx: dict) -> list[dict]:
    """Supplementary Data 1, 3, 5 and 6: the sites whose methylation fell at least two-fold (or vanished) on PRMT7,
    PRMT4 or PRMT5 knockdown or SGC3027 treatment of HEK293 cells, each with its accession, unmarked peptide and the
    position in the protein and the peptide; the 'Raw data' sheets are not read. The tables give no localisation score
    and the search allowed methylation of lysine as well as arginine, so a site is kept only where the peptide holds no
    other arginine or lysine (one_methylatable_residue)."""
    out = []
    for path in paths:
        wb = load_workbook(path, read_only=True, data_only=True)
        for ws in wb.worksheets:
            m = LI_SHEETS.search(ws.title.replace('SGC3027regulated', 'SGC3027 regulated'))
            if not m:
                continue
            perturbation, degree = m.group(1).upper(), {'me1': 'Rme1', 'me2': 'Rme2'}[m.group(2).lower()]
            H, data = sheet_rows(path, ws.title, 1)
            stats[f'rows_read:{ws.title}'] += len(data)
            for r in data:
                acc, peptide = text(r[H['Accession number']]), text(r[H['Peptide sequence']]).upper()
                try:
                    position, off = int(float(r[H['Rme1 position in protein']])), int(float(r[H['Rme1 position in peptide']]))
                except (TypeError, ValueError):
                    stats['rows_dropped:no_position_given'] += 1
                    continue
                if not ACCESSION.fullmatch(acc) or not 1 <= off <= len(peptide) or peptide[off - 1] != 'R':
                    stats['rows_dropped:not_uniprot_or_site_not_R'] += 1
                    continue
                localised = mass_localised(peptide, {off: degree}, 'RK')
                if not localised:
                    stats['site_rows_dropped:other_arginine_or_lysine_in_peptide'] += 1
                    continue
                out.append(site_row([acc], text(r[H['Gene name']]), position, degree, 'HEK293', LI_ENRICHMENT[perturbation],
                                    ws.title, peptide, off, note=[f'{perturbation} {"knockdown" if perturbation != "SGC3027" else "inhibitor"}'],
                                    localised=localised))
    return out


# ------------------------------------------------------------------------------------------------- Fedoriw 2019
FEDORIW_ENRICHMENT = {'MMA': 'MMA antibody (#12235)', 'ADMA': 'ADMA antibody (#13474)', 'SDMA': 'SDMA antibody (#13563)'}


def fedoriw_peptide(marked: str) -> tuple[str, dict[int, str]]:
    """'LR*GPGGPR@GGLGGGMR' -> ('LRGPGGPRGGLGGGMR', {2: 'Rme1', 8: 'Rme2'}): R* monomethyl, R@ dimethyl, M# oxidised
    (the marks every R* row of the MMA immunoprecipitations and every R@ row of the ADMA and SDMA ones carry)."""
    residues, marks = '', {}
    for ch in marked:
        if 'A' <= ch <= 'Z':
            residues += ch
        elif ch in '*@' and residues.endswith('R'):
            marks[len(residues)] = 'Rme1' if ch == '*' else 'Rme2'
    return residues, marks


def read_fedoriw(paths: list[Path], stats: Counter, ctx: dict) -> list[dict]:
    """Tables S5 (five DLBCL and pancreatic cancer cell lines, type I PRMT inhibitor) and S7 (type I PRMT and PRMT5
    inhibitors combined): peptides whose MMA, ADMA or SDMA changed, each with its gene symbol and the methyl marks; no
    accession, position or localisation score. The search methylated arginines only, so a mark is kept where the
    peptide holds no other arginine, or every arginine is dimethylated. Its position is the peptide's single occurrence
    in the reviewed entries of the gene."""
    out = []
    for path in paths:
        wb = load_workbook(path, read_only=True, data_only=True)
        for ws in wb.worksheets:
            rows = [r for r in ws.iter_rows(values_only=True)]
            header = next((k for k, r in enumerate(rows) if r and text(r[0]) == 'Modification'), None)
            if header is None:
                continue
            H = {text(h): k for k, h in enumerate(rows[header]) if h is not None}
            table = 'Table S7' if 'Combination (Log2 Fold Change)' in H else 'Table S5'
            cell = ws.title if table == 'Table S5' else 'not stated (Table S7)'
            for r in rows[header + 1:]:
                if not r or not text(r[H['Peptide']]):
                    continue
                stats[f'rows_read:{table} {ws.title}'] += 1
                ip, gene = text(r[H['Modification']]), text(r[H['Protein']])
                peptide, marks = fedoriw_peptide(text(r[H['Peptide']]))
                localised = mass_localised(peptide, marks, 'R')
                if not marks:
                    stats['rows_dropped:no_arginine_mark'] += 1
                    continue
                if not localised:
                    stats['rows_dropped:other_arginine_in_peptide'] += len(marks)
                    continue
                accs = ctx['genes'].get(gene.upper(), [])
                hits = [(a, m.start() + 1) for a in accs for m in re.finditer(f'(?={peptide})', ctx['reference'].get(a, ''))]
                if len(hits) != 1:
                    stats['rows_dropped:gene_symbol_not_in_reference' if not accs else 'rows_dropped:peptide_not_once_in_gene_entries'] += len(marks)
                    continue
                acc, start = hits[0]
                for off, degree in marks.items():
                    out.append(site_row([acc], gene, start + off - 1, degree, cell, FEDORIW_ENRICHMENT.get(ip, ip), table,
                                        peptide, off, note=[f'{ip} immunoprecipitation'], localised=localised))
    return out


# ------------------------------------------------------------------------------------------------------ Lim 2020
LIM_MIN_MODSCORE = 13.0
LIM_TABLES = {'Table S1': 'MMA antibody (#12235)', 'Table S2': 'ADMA antibody (#13474)'}


def read_lim(paths: list[Path], stats: Counter, ctx: dict) -> list[dict]:
    """Tables S1 (MMA) and S2 (ADMA): peptide-spectrum matches from pooled colorectal cancer tissue, each with the
    site(s), a 15-residue window per site (lower-case r at the site), the peptide with R* (monomethyl) and R@
    (asymmetric dimethyl) marks and a ModScore per replicate injection (one value per site, ';'-separated). A site is
    kept where a replicate's ModScore exceeds 13 (the paper's threshold). Some peptides write I for L; the placement
    restores them."""
    out = []
    for path in sorted(paths, key=lambda p: p.name):
        table = 'Table S1' if 'MOESM1' in path.name else 'Table S2'
        ws = open_workbook(path).worksheets[0]
        rows = [r for r in ws.iter_rows(values_only=True)]
        header = next(k for k, r in enumerate(rows) if r and text(r[0]) == 'Protein Type')
        names = [text(h) for h in rows[header]]
        col = {n: k for k, n in enumerate(names) if n}
        score_cols = [col['ModScore'], col['ModScore'] + 1]
        for r in rows[header + 2:]:
            if not r or not text(r[col['Site']]) or not text(r[col['Accession']]):
                continue
            stats[f'rows_read:{table}'] += 1
            acc = text(r[col['Accession']])
            if not ACCESSION.fullmatch(acc):
                stats['rows_dropped:not_uniprot_accession'] += 1
                continue
            sites = [t.strip() for t in text(r[col['Site']]).split(',')]
            windows = [t.strip() for t in re.split(r'[;,]', text(r[col['-7/+7 peptide']]))]
            marked = text(r[col['Peptide']])
            peptide, marks = fedoriw_peptide(marked.replace('#', ''))
            scores = [[float(x) for x in text(r[c]).split(';') if x.strip()] for c in score_cols if text(r[c])]
            if len(windows) != len(sites):
                stats['site_rows_dropped:windows_do_not_pair_with_sites'] += len(sites)
                continue
            for i, (site, window) in enumerate(zip(sites, windows)):
                per_site = [v[i] if len(v) == len(sites) else (v[0] if len(v) == 1 else None) for v in scores]
                best = max((x for x in per_site if x is not None), default=None)
                if not re.fullmatch(r'\d+', site) or len(window) != 15 or window[7] != 'r':
                    stats['site_rows_dropped:no_position_or_window'] += 1
                    continue
                if best is None or best <= LIM_MIN_MODSCORE:
                    stats['site_rows_dropped:modscore_13_or_less'] += 1
                    continue
                fits = [off for off in marks if all(window[7 + k] in 'xX' or not 0 <= off - 1 + k < len(peptide)
                                                    or peptide[off - 1 + k].replace('I', 'L') == window[7 + k].upper().replace('I', 'L')
                                                    for k in range(-7, 8))]
                if len(fits) != 1:
                    stats['site_rows_dropped:mark_not_matched_in_peptide'] += 1
                    continue
                off = fits[0]
                out.append(site_row([acc], text(r[col['Gene Name']]), int(site), marks[off], 'colorectal cancer tissue',
                                    LIM_TABLES[table], table, peptide, off, window=window.upper().replace('X', '_'),
                                    note=[f'ModScore {best:g}']))
    return out


# ------------------------------------------------------------------------------------------------------ Wei 2020
WEI_MIN_LOC = 0.75


def read_wei(paths: list[Path], stats: Counter, ctx: dict) -> list[dict]:
    """Table 1 ('ADMA sites'): MaxQuant sites from HPDE6c7 and PANC-1 cells enriched with the ADMA antibody, with the
    peptide and its per-residue probabilities, the position in peptide and protein, the 31-residue window and the
    localisation probability (the paper keeps > 0.75 and Andromeda score > 40; every row meets both). The 'Protein names'
    and 'Gene names' headers are swapped in the table; the gene is the shorter of the two cells."""
    H, data = sheet_rows(paths[0], 'ADMA sites', 3)
    stats['rows_read:Table 1'] += len(data)
    out = []
    for r in data:
        accs = accessions(r[H['Accessions']])
        if not accs or not all(ACCESSION.fullmatch(a) for a in accs):
            stats['rows_dropped:not_uniprot_accession'] += 1
            continue
        loc = float(r[H['Localization prob']])
        if loc <= WEI_MIN_LOC:
            stats['rows_dropped:localisation_0.75_or_less'] += 1
            continue
        peptide = re.sub(r'\([^)]*\)', '', text(r[H['Dimethyl (R) Probabilities']])).upper()
        off, position = int(float(r[H['Position in peptide']])), int(float(r[H['Position within protein']]))
        if not 1 <= off <= len(peptide) or peptide[off - 1] != 'R':
            stats['rows_dropped:site_not_R_in_peptide'] += 1
            continue
        window = text(r[H['Sequence window']]).split(';')[0]
        gene = min((text(r[H['Gene names']]), text(r[H['Protein names']])), key=len).split(';')[0]
        out.append(site_row(accs, gene, position, 'Rme2', 'HPDE6c7; PANC-1', 'ADMA antibody (#13474)', 'Table 1', peptide,
                            off, window=window, loc=loc))
    return out


# ---------------------------------------------------------------------------------------------- Maron 2021, own
MARON_DEGREE = {'me1': 'Rme1', 'me2': 'Rme2'}


PD_MODS = re.compile(r'\d+x(Methyl|Dimethyl)\s*\[([^\]]*)\]')


def pd_site_calls(path: Path, stats: Counter) -> dict[str, dict[int, dict[str, float]]]:
    """A Proteome Discoverer 2.4 peptide-group export of the lab's immunoprecipitations (Byonic search, ptmRS site
    probabilities): unmodified peptide -> {offset: {degree: best site probability}} for every arginine whose methyl
    group PD places. PD names a residue in the Modifications column only at site probability >= 75 ('R8(96)'); a methyl
    group it cannot place reads '[R]', and 'R/K' is a group it cannot tell between an arginine and a lysine."""
    wb = load_workbook(path, read_only=True, data_only=True)
    rows = wb.worksheets[0].iter_rows(values_only=True)
    H = {text(h): k for k, h in enumerate(next(rows))}
    calls = defaultdict(lambda: defaultdict(dict))
    for r in rows:
        if not any(v not in (None, '') for v in r):
            continue
        stats[f'pd_rows_read:{path.name}'] += 1
        peptide = re.sub(r'^\[.*?\]\.|\.\[.*?\]$', '', text(r[H['Annotated Sequence']])).upper()
        for kind, body in PD_MODS.findall(text(r[H['Modifications']])):
            for tok in (t.strip() for t in body.split(';') if t.strip()):
                m = re.fullmatch(r'R(\d+)\(([\d.]+)\)', tok)
                if not m:
                    stats['pd_methyl_groups_not_placed_on_one_arginine'] += 1
                    continue
                off, prob = int(m.group(1)), float(m.group(2))
                if prob < MARON_MIN_PTMRS or not 1 <= off <= len(peptide) or peptide[off - 1] != 'R':
                    stats['pd_methyl_groups_below_75_or_not_R'] += 1
                    continue
                degree = DEGREE[kind]
                calls[peptide][off][degree] = max(prob, calls[peptide][off].get(degree, 0.0))
    return calls


MARON_MIN_PTMRS = 75.0


def read_maron(paths: list[Path], stats: Counter, ctx: dict) -> list[dict]:
    """Tables S2 (each methylated residue with a peptide that carries it) and S3 (each peptide with all its methyl
    marks) of the A549 PTMScan data (sequential Rme1, Rme2a and Rme2s immunoprecipitation; Proteome Discoverer 2.4 with
    Byonic, 1 % FDR), with the lab's Proteome Discoverer peptide-group exports of the same immunoprecipitations ('IP
    MethylatedPeptide Groups', 'OG IP MethylatedPeptide Groups'; pd_site_calls). The tables give no localisation score
    and the search methylated lysines as well as arginines (mono- and dimethyl). A mark of the tables is kept where its
    peptide holds no other arginine or lysine (one_methylatable_residue), or where an export localises a methyl group to
    that arginine of the same peptide at ptmRS >= 75 % (localisation_score). When the export places the other degree
    there, the record takes the export's degree and the note says what the table gives. The exports' accessions (NCBI
    GI numbers in the earlier export) are not used: the tables name the protein. GluC peptides and the runs the exports
    do not hold are kept only by mass."""
    tables = [p for p in paths if p.name.startswith('maron2021_table')]
    calls = {}
    for path in sorted(set(paths) - set(tables), key=lambda p: p.name):
        for peptide, sites in pd_site_calls(path, stats).items():
            for off, degrees in sites.items():
                for degree, prob in degrees.items():
                    have = calls.setdefault(peptide, {}).setdefault(off, {})
                    have[degree] = max(prob, have.get(degree, 0.0))
    out = []
    for path in sorted(tables, key=lambda p: p.name):
        wb = load_workbook(path, read_only=True, data_only=True)
        for ws in wb.worksheets:
            H, data = sheet_rows(path, ws.title, 1)
            table = 'Table S2' if 'Position.In.Protein' in H else 'Table S3'
            stats[f'rows_read:{table} {ws.title}'] += len(data)
            for r in data:
                acc, peptide, begin = text(r[H['Accession']]), text(r[H['Sequence']]).upper(), int(float(r[H['Begin.Pos']]))
                mods = text(r[H['Modification' if table == 'Table S2' else 'Modifications']])
                marks = {int(m.group(1)) - begin + 1: MARON_DEGREE[m.group(2)] for m in re.finditer(r'R(\d+)(me1|me2)', mods)}
                marks = {off: d for off, d in marks.items() if 1 <= off <= len(peptide) and peptide[off - 1] == 'R'}
                if not marks or not ACCESSION.fullmatch(acc):
                    stats['rows_dropped:no_arginine_mark_or_accession'] += 1
                    continue
                by_mass = mass_localised(peptide, marks, 'RK')
                for off, degree in marks.items():
                    note, loc, localised = [f'{text(r[H["Protease"]])} peptide'], None, by_mass
                    if not localised:
                        placed = calls.get(peptide, {}).get(off, {})
                        if not placed:
                            stats['site_marks_dropped:not_localised_by_mass_or_ptmrs'] += 1
                            continue
                        if degree not in placed:
                            if len(placed) > 1:
                                stats['site_marks_dropped:ptmrs_places_both_degrees'] += 1
                                continue
                            note.append(f'{table} gives {MARK[degree][2:-1]}; the Proteome Discoverer export localises {MARK[next(iter(placed))][2:-1]}')
                            degree = next(iter(placed))
                            stats['site_marks_degree_from_ptmrs'] += 1
                        localised, loc = 'localisation_score', placed[degree] / 100
                        stats['site_marks_localised_by_ptmrs'] += 1
                    out.append(site_row([acc], text(r[H['Gene.Symbol']]), begin + off - 1, degree, 'A549',
                                        'Rme1, Rme2a and Rme2s antibodies (sequential)', table, peptide, off,
                                        note=note, loc=loc, localised=localised))
    return out


# ------------------------------------------------------------------------------------- Musiani 2019, Larsen 2016
MIN_LOC = 0.75   # MaxQuant localisation probability (class I sites)
MUSIANI_ENRICHMENT = {'Mono': 'monomethylarginine antibody', 'Sym': 'SDMA antibody', 'Input': 'no enrichment (input)'}


def musiani_site(peptide: str, marks: dict, window: str, degree: str) -> int | None:
    """The offset of the arginine a 15- or 13-residue window ('z' or '_' padding) centres on, among the peptide's marks
    of that degree."""
    candidates = [off for off, m in marks.items() if peptide[off - 1] == 'R' and MAXQUANT_MARK.get(m) == degree]
    return window_offset(peptide, candidates, window.replace('z', '_')) if candidates else None


def read_musiani(paths: list[Path], stats: Counter, ctx: dict) -> list[dict]:
    """The article's tables. Table S1 (Spreadsheet 1): every methyl site of the HeLa GSK591 SILAC experiments (MaxQuant;
    monomethylarginine and SDMA antibodies, and the input), one row per site and identification, with the leading protein,
    the position, the 15-residue window and the PTM localisation probability. Table S2: the heavy-methyl SILAC pairs of
    the orthogonal validation, Spreadsheet 1 with the leading protein list, LcPrb of the heavy and light peptides and a
    15-residue window (rows of class Mismatched, whose heavy and light peptides differ, are left out), and Spreadsheet 2
    with the gene only, a 13-residue window and LcPrb '> 0.75' (its entry is the one reviewed entry of the gene that holds
    the peptide once, --entries, --reference). Kept: arginine sites at localisation probability >= 0.75 (in Table S2,
    the best of the heavy and light peptides)."""
    out = []
    for path in paths:
        wb = load_workbook(path, read_only=True, data_only=True)
        if 'table_s1' in path.name.lower():
            rows = list(wb['Spreadsheet 1'].iter_rows(values_only=True))
            head = next(k for k, r in enumerate(rows) if r and text(r[0]) == 'Modified sequence')
            H = {text(h): k for k, h in enumerate(rows[head]) if h is not None}
            for r in rows[head + 1:]:
                if not r or not text(r[0]):
                    continue
                stats['rows_read:Table S1'] += 1
                if text(r[H['Residue']]) != 'R':
                    stats['rows_dropped:not_arginine'] += 1
                    continue
                loc = float(r[H['PTM Localization Probability']])
                if loc < MIN_LOC:
                    stats['rows_dropped:localisation_below_0.75'] += 1
                    continue
                acc, degree = text(r[H['Leading UniprotID']]), MAXQUANT_MARK[text(r[H['Modification']])]
                peptide, marks = maxquant_peptide(text(r[H['Modified sequence']]))
                window = text(r[H['Modification Window']])
                off = musiani_site(peptide, marks, window, degree)
                if not ACCESSION.fullmatch(acc) or off is None:
                    stats['rows_dropped:accession_unreadable_or_window_not_on_one_mark'] += 1
                    continue
                exp = text(r[H['Experiment']])
                out.append(site_row([acc], text(r[H['Gene Names']]).split(';')[0], int(float(r[H['Position in protein']])), degree,
                                    'HeLa', MUSIANI_ENRICHMENT.get(exp.split('_')[-1], exp), 'Table S1', peptide, off,
                                    window=window.replace('z', '_'), loc=loc, note=[f'SILAC {exp}']))
            continue
        for sheet in ('Spreadsheet 1', 'Spreadsheet 2'):
            rows = list(wb[sheet].iter_rows(values_only=True))
            head = next(k for k, r in enumerate(rows) if r and text(r[0]).startswith(('Sequence H', 'Modified Sequence H')))
            H = {text(h).lower(): k for k, h in enumerate(rows[head]) if h is not None}
            get = lambda r, name, side: r[H[f'{name} {side}'.lower()]] if f'{name} {side}'.lower() in H else None
            for r in rows[head + 1:]:
                if not r or not any(v not in (None, '') for v in r):
                    continue
                stats[f'rows_read:Table S2 {sheet}'] += 1
                if 'class' in H and text(r[H['class']]) == 'Mismatched':
                    stats['rows_dropped:heavy_and_light_mismatched'] += 1
                    continue
                side = next((x for x in ('H', 'L') if text(get(r, 'Residue', x)) == 'R'), None)
                if side is None:
                    stats['rows_dropped:not_arginine'] += 1
                    continue
                probs = [text(get(r, 'LcPrb', f'{x} (best)')) for x in ('H', 'L')]
                loc = max((float(v) for v in probs if re.fullmatch(r'[\d.]+', v)), default=None)
                if not (loc is not None and loc >= MIN_LOC or any(v == '> 0.75' for v in probs)):
                    stats['rows_dropped:localisation_below_0.75'] += 1
                    continue
                degree = MAXQUANT_MARK.get(text(get(r, 'Modification', side)))
                peptide, marks = maxquant_peptide(text(get(r, 'Modified Sequence', side)))
                window, pos = text(get(r, 'Window', side)), int(float(get(r, 'Position', side)))
                off = musiani_site(peptide, marks, window, degree) if degree else None
                if off is None:
                    stats['rows_dropped:window_not_on_one_mark'] += 1
                    continue
                gene = text(get(r, 'Gene', side)).split(';')[0]
                accs = [a for a in text(get(r, 'All Proteins', side) or get(r, 'All proteins', side)).split(';') if a]
                if not accs:   # Spreadsheet 2: the gene's reviewed entry that holds the peptide once
                    entries = ctx['genes'].get(gene.upper(), [])
                    hits = [(a, m.start() + 1) for a in entries for m in re.finditer(f'(?={peptide})', ctx['reference'].get(a, ''))]
                    if len(hits) != 1:
                        stats['rows_dropped:peptide_not_once_in_gene_entries'] += 1
                        continue
                    accs = [hits[0][0]]
                if not all(ACCESSION.fullmatch(a) for a in accs):
                    stats['rows_dropped:accession_unreadable'] += 1
                    continue
                out.append(site_row(accs, gene, pos, degree, 'HeLa', 'heavy-methyl SILAC, methylarginine antibodies',
                                    f'Table S2 {sheet}', peptide, off, window=window, loc=loc,
                                    note=['LcPrb > 0.75 (as given)'] if loc is None else []))
    return out


def read_larsen(paths: list[Path], stats: Counter, ctx: dict) -> list[dict]:
    """Table S1: the 8,900 monomethylarginine sites of the HEK293 experiments (MaxQuant), each with the leading protein,
    the position on it, the localisation probability and the 31-residue window. Kept: localisation probability >= 0.75
    (MaxQuant class I; the table also lists sites below it). The leading protein can be an isoform or an unreviewed
    entry; the gene names the reviewed entry, and the window places the site."""
    H, data = sheet_rows(paths[0], load_workbook(paths[0], read_only=True).sheetnames[0], 1)
    stats['rows_read:Table S1'] += len(data)
    out = []
    for r in data:
        acc, window = text(r[H['Protein']]), text(r[H['Sequence window']])
        if text(r[H['Amino acid']]) != 'R' or len(window) != 31 or window[15] != 'R':
            stats['rows_dropped:not_arginine_or_window_not_31_centred'] += 1
            continue
        if not ACCESSION.fullmatch(acc):
            stats['rows_dropped:not_uniprot_accession'] += 1
            continue
        loc = float(r[H['Localization prob']])
        if loc < MIN_LOC:
            stats['rows_dropped:localisation_below_0.75'] += 1
            continue
        gene = text(r[H['Gene names']]).split(';')[0]
        out.append(site_row([acc], gene, int(float(r[H['Position within proteins']])), 'Rme1', 'HEK293',
                            'methylarginine antibodies after high-pH fractionation', 'Table S1', window=window, loc=loc))
    return out


# ------------------------------------------------------------------------------------------------- UniProt, IEDB
def read_uniprot(paths: list[Path], stats: Counter, ctx: dict) -> list[dict]:
    """The methylated arginines UniProt 2026_03 annotates (MOD_RES, staged by scripts/00_stage_uniprot_methylarginine.py)
    with experimental evidence (ECO:0000269): curated from the literature, residue by residue. Annotations from large-scale
    proteomics (ECO:0007744) repeat the studies read here; by similarity (ECO:0000250) is not evidence of the site."""
    with paths[0].open(encoding='utf-8') as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    stats['rows_read:MOD_RES'] += len(rows)
    out = []
    for r in rows:
        if 'ECO:0000269' not in r['evidence']:
            stats[f'rows_dropped:evidence_{r["evidence_class"]}'] += 1
            continue
        state = r['state']
        degree = 'Rme1' if state == 'Rme1' else 'Rme2' if state.startswith('Rme2') else ''
        note = {'Rme2a': 'asymmetric (UniProt)', 'Rme2s': 'symmetric (UniProt)'}.get(state, '')
        out.append(site_row([r['uniprot_accession']], '', int(r['position']), degree, '', '', 'MOD_RES', note=[note] if note else [],
                            localised='curated_experimental'))
    return out


def read_iedb(paths: list[Path], stats: Counter, ctx: dict) -> list[dict]:
    """IEDB epitopes (HLA-presented peptides) with a methylated arginine: the 85 the iPTMnet-format PTM master table
    cites, as the IEDB query API returns them (iedb_epitopes_methylarginine.json): sequence, modified residue
    ('METH(R6)', degree not given) and the curated source antigen with the epitope's start. The search behind them is not
    known, so a site is kept only where the epitope holds no other arginine or lysine."""
    epitopes = json.loads(paths[0].read_text())
    stats['rows_read:epitopes'] += len(epitopes)
    out = []
    for e in epitopes:
        peptide = e['linear_sequence'].upper()
        mods = [int(m.group(1)) for d in e['structure_descriptions'] for m in re.finditer(r'METH\(R(\d+)\)', d)]
        antigens = [a for a in e['curated_source_antigens'] or [] if ACCESSION.fullmatch(a['accession'].split('.')[0])]
        if len(mods) != 1 or not antigens:
            stats['rows_dropped:not_one_arginine_mark_or_no_uniprot_antigen'] += 1
            continue
        off = mods[0]
        localised = mass_localised(peptide, {off: ''}, 'RK')
        if not localised:
            stats['rows_dropped:other_arginine_or_lysine_in_epitope'] += 1
            continue
        a = antigens[0]
        out.append(site_row([a['accession'].split('.')[0]], '', a['starting_position'] + off - 1, '', '', 'HLA-presented peptides',
                            'IEDB epitope', peptide, off, note=[f'IEDB epitope {e["structure_id"]}'], localised=localised))
    return out


# -------------------------------------------------------------------------------------------------------- merge
# study: (study_id, reader, the context it needs: reference FASTA, gene index, Maron 2021 Table S5)
STUDIES = {'hartel2019': ('Hartel2019', read_hartel, ()), 'geoghegan2015': ('Geoghegan2015', read_geoghegan, ()),
           'sylvestersen2014': ('Sylvestersen2014', read_sylvestersen, ()), 'lu2025': ('Lu2025', read_lu, ('reference',)),
           'ma2023': ('Ma2023', read_ma, ()), 'guo2014': ('Guo2014', read_guo, ('genes',)), 'fong2019': ('Fong2019', read_fong, ()),
           'li2021': ('Li2021', read_li, ()), 'fedoriw2019': ('Fedoriw2019', read_fedoriw, ('reference', 'genes')),
           'lim2020': ('Lim2020', read_lim, ()), 'wei2020': ('Wei2020', read_wei, ()), 'maron2021': ('Maron2021', read_maron, ()),
           'musiani2019': ('Musiani2019', read_musiani, ('reference', 'genes')), 'larsen2016': ('Larsen2016', read_larsen, ()),
           'uniprot': ('UniProt', read_uniprot, ()), 'iedb': ('IEDB', read_iedb, ())}


def table_key(label: str):
    return [int(t) if t.isdigit() else t for t in re.split(r'(\d+)', label)]


def fmt_prob(x: float | None) -> str:
    return '' if x is None else f'{x:.4f}'.rstrip('0').rstrip('.')


def merge(study_id: str, rows: list[dict]) -> list[dict]:
    groups = defaultdict(list)
    for row in rows:
        key = (row['accessions'][0], row['position'], row['peptide'] or row['window'], row['offset'] or 0, row['degree'])
        groups[key].append(row)
    records = []
    for (acc, pos, _, _, degree), grp in groups.items():
        accs = list(dict.fromkeys(a for r in grp for a in r['accessions']))
        peptide = grp[0]['peptide']
        marked = peptide[:grp[0]['offset'] - 1] + MARK[degree] + peptide[grp[0]['offset']:] if peptide else ''
        calls = hartel_rme2_type(grp) if study_id == 'Hartel2019' else set()
        locs = [r['loc'] for r in grp if r['loc'] is not None]
        enrich = sorted({e for r in grp for e in (r.get('enrichments') or [r['enrichment']]) if e})
        notes = (hartel_merge_note(grp) if study_id == 'Hartel2019' else geoghegan_merge_note(grp) if study_id == 'Geoghegan2015'
                 else [n for r in grp for n in r['note'] if not n.startswith('ModScore ')])
        scores = [float(n.split()[1]) for r in grp for n in r['note'] if n.startswith('ModScore ')]
        if scores:  # Guo 2014, Lim 2020: the best ModScore of the record's rows
            notes.insert(0, f'ModScore {max(scores):g}')
        records.append({
            'study_id': study_id, 'reported_accession': '; '.join(accs),
            'reported_gene': next((r['gene'] for r in grp if r['gene']), ''), 'reported_position': pos,
            'peptide': marked, 'window': next((r['window'] for r in grp if r['window']), ''), 'methyl_degree': degree,
            'rme2_type': '; '.join(sorted(calls)) if degree == 'Rme2' else '',
            'cell_type': '; '.join(sorted({r['cell'] for r in grp})), 'enrichment': '; '.join(enrich),
            'localisation_probability': fmt_prob(max(locs) if locs else None), 'n_source_rows': len(grp),
            'source_tables': '; '.join(sorted({r['table'] for r in grp}, key=table_key)),
            'note': '; '.join(dict.fromkeys(notes)),
            'localised_by': min((r['localised'] for r in grp), key=LOCALISED_BY.index),
        })
    records.sort(key=lambda r: (r['reported_accession'].split('; ')[0], r['reported_position'], r['methyl_degree'],
                                r['peptide'], r['window']))
    return records


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--study', required=True, choices=sorted(STUDIES))
    ap.add_argument('--input', required=True, nargs='+', help='the study supplementary table(s), as in the usage above')
    ap.add_argument('--out', required=True)
    ap.add_argument('--reference', default='', help='the reference FASTA (lu2025: pairs each listed site with its PSM peptide; fedoriw2019, musiani2019: places peptides named by gene only)')
    ap.add_argument('--entries', default='', help='the reference entries table (guo2014, fedoriw2019, musiani2019: gene symbols to entries)')
    ap.add_argument('--history', default='', help='the accession history table (gene synonyms, with --entries)')
    args = ap.parse_args()

    study_id, reader, needs = STUDIES[args.study]
    paths = [Path(p) for p in args.input]
    stats = Counter()
    ctx = {}
    for need, given in (('reference', args.reference), ('genes', args.entries and args.history)):
        if need in needs and not given:
            raise SystemExit(f'{args.study} needs --{ {"genes": "entries and --history"}.get(need, need) }')
    if 'reference' in needs:
        ctx['reference'] = read_fasta(Path(args.reference))
    if 'genes' in needs:
        ctx['genes'] = gene_index(Path(args.entries), Path(args.history))
    rows = reader(paths, stats, ctx)
    records = merge(study_id, rows)
    for rec in records:
        if any(ch in rec['peptide'] for ch in ',;'):
            raise SystemExit(f'peptide with a separator: {rec}')
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(records)

    summary = {
        'study': study_id,
        'inputs': {p.name: hashlib.md5(p.read_bytes()).hexdigest() for p in paths},
        'rows_read': {k.split(':', 1)[1]: v for k, v in stats.items() if k.startswith('rows_read:')},
        'filters_and_counts': {k: v for k, v in sorted(stats.items()) if not k.startswith('rows_read:')},
        'site_rows_kept': len(rows),
        'records_written': len(records),
        'by_degree': dict(sorted(Counter(r['methyl_degree'] for r in records).items())),
        'rme2_type': dict(sorted(Counter(r['rme2_type'] or '(empty)' for r in records if r['methyl_degree'] == 'Rme2').items())),
        'with_peptide': sum(1 for r in records if r['peptide']),
        'with_window': sum(1 for r in records if r['window']),
        'distinct_accession_position': len({(r['reported_accession'].split('; ')[0], r['reported_position']) for r in records}),
        'out': str(out),
    }
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()

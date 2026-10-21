#!/usr/bin/env python3
"""Stage the methylarginine records of five proteomics studies from their supplementary tables (inputs/sources/).

    python scripts/00_stage_study_sites.py --study hartel2019 \
        --input 153496_0_supp_351381_ptfjl4.xlsx 153496_0_supp_351382_ptfjl4.xlsx \
        --out inputs/sources/hartel2019_methylarginine_records.tsv
    python scripts/00_stage_study_sites.py --study geoghegan2015 --input ncomms7758-s2.xlsx \
        --out inputs/sources/geoghegan2015_methylarginine_records.tsv
    python scripts/00_stage_study_sites.py --study sylvestersen2014 --input mcp.O113.032748-4.xlsx mcp.O113.032748-6.xlsx \
        --out inputs/sources/sylvestersen2014_methylarginine_records.tsv
    python scripts/00_stage_study_sites.py --study lu2025 --input 44319_2025_599_MOESM3_ESM.xlsx \
        --out inputs/sources/lu2025_methylarginine_records.tsv
    python scripts/00_stage_study_sites.py --study ma2023 \
        --input Table_S1_Identified_and_quantified_ADMA_methylation_sites_06012022.xlsx \
        --out inputs/sources/ma2023_methylarginine_records.tsv

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
under), not always 'Protein'; that accession leads reported_accession, followed by the others as written. Kept, per
the methods: q.value <= 0.05 (5 % peptide FDR; with the methyl-SILAC pair the set's FDR is 0.93 %), IP.Prob > 0.25,
and ModLS > 8 or, below it, every arginine methylated, every methylated arginine missed by trypsin or every methylated
arginine with PTMscore >= 0.95 (every row meets them). PTM.Localization.Probs is reported as localisation_probability
but not used as a filter: it gives < 0.75 to many arginines of peptides whose every arginine is methylated. Rows whose
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
Tudor domain (Proteome Discoverer, 1 % FDR). EV1 has no list of monomethylarginine sites, only PSMs without protein
positions, so only the dimethyl list is staged. Sites of the 'Validated peptides' sheet (hmSILAC) are noted.

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

COLUMNS = ['study_id', 'reported_accession', 'reported_gene', 'reported_position', 'peptide', 'window', 'methyl_degree',
           'rme2_type', 'cell_type', 'enrichment', 'localisation_probability', 'n_source_rows', 'source_tables', 'note']
ACCESSION = re.compile(r'([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(-\d+)?')
DEGREE = {'Methyl': 'Rme1', 'Dimethyl': 'Rme2'}
MARK = {'Rme1': 'R(me1)', 'Rme2': 'R(me2)'}


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


def site_row(accs, gene, position, degree, cell, enrichment, table, peptide='', offset=None, window='', loc=None,
             rme2=None, note=()) -> dict:
    return {'accessions': list(accs), 'gene': gene, 'position': position, 'degree': degree, 'cell': cell,
            'enrichment': enrichment, 'table': table, 'peptide': peptide, 'offset': offset, 'window': window,
            'loc': loc, 'rme2': rme2 or set(), 'note': list(note)}


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


def read_hartel(paths: list[Path], stats: Counter) -> list[dict]:
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
                if prob is not None and prob < HARTEL_MIN_LOC:
                    stats['site_rows_dropped:localisation_below_75'] += 1
                    continue
                if pos_of and not consistent:
                    stats['site_rows_dropped:position_inconsistent_with_peptide'] += 1
                    continue
                degree = DEGREE[k]
                rme2 = hartel_nl_call(off, localised, nl_class, n_me2, has_mma, andromeda) if nl_tested and degree == 'Rme2' else set()
                row = site_row(masters, gene, pos_of.get((res, off)), degree, cell, enrichment, label, peptide, off,
                               loc=None if prob is None else prob / 100, rme2=rme2, note=[n for n in note if n])
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
GEOGHEGAN_MAX_Q, GEOGHEGAN_MIN_IP, GEOGHEGAN_MIN_MODLS, GEOGHEGAN_MIN_PTMSCORE = 0.05, 0.25, 8, 0.95
GEOGHEGAN_ANTIBODY = {'D5A12': ['anti-MMA antibody D5A12'], 'MeR4100': ['anti-MMA antibody MeR4100'],
                      'D5A12+MeR4100': ['anti-MMA antibody D5A12', 'anti-MMA antibody MeR4100']}


def geoghegan_degree(mass: float) -> str:
    """Light or heavy (13CD3) methyl: +14.02 / +18.04 mono, +28.03 / +36.08 di."""
    for degree, masses in (('Rme1', (14.02, 18.04)), ('Rme2', (28.03, 36.08))):
        if any(abs(mass - m) < 0.05 for m in masses):
            return degree
    return ''


def geoghegan_localised(modls, peptide: str, sites: list[int], probs: dict[int, float], protease: str) -> bool:
    """The paper's rule: ModLS > 8, or every arginine methylated, every methylated arginine missed by trypsin, or
    every methylated arginine with PTMscore >= 0.95."""
    if modls == 'Unambiguous' or float(modls) > GEOGHEGAN_MIN_MODLS:
        return True
    if all(off in sites for off, aa in enumerate(peptide, 1) if aa == 'R'):
        return True
    if protease == 'Trypsin' and all(off < len(peptide) for off in sites):
        return True
    return all(probs.get(off, 0) >= GEOGHEGAN_MIN_PTMSCORE for off in sites)


def read_geoghegan(paths: list[Path], stats: Counter) -> list[dict]:
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
        probs, i = {}, 0
        for aa, p in re.findall(r'([A-Z])(?:\(([\d.]+)\))?', text(r[H['PTM.Localization.Probs']])):
            i += 1
            if p:
                probs[i] = float(p)
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
        if not geoghegan_localised(r[H['ModLS']], peptide, [off for off, _ in sites], probs, text(r[H['Protease']])):
            stats['site_rows_dropped:modls_rule_not_met'] += len(sites)
            continue
        cell = text(r[H['Cell']])
        enrichment = GEOGHEGAN_ANTIBODY[text(r[H['Antibody']])]
        for m, (off, degree) in zip(reported, sites):
            row = site_row(accs, gene, int(m.group(1)), degree, cell, '', 'ST1', peptide, off, loc=probs.get(off))
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


def read_sylvestersen(paths: list[Path], stats: Counter) -> list[dict]:
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
def read_lu(paths: list[Path], stats: Counter) -> list[dict]:
    wb = load_workbook(paths[0], read_only=True, data_only=True)
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
        out.append(site_row([acc], '', pos, 'Rme2', 'HeLa', 'SMN Tudor', 'EV1 Arginine dimethylation sites', note=note))
    return out


# ------------------------------------------------------------------------------------------------------ Ma 2023
MA_MIN_LOC = 0.75


def read_ma(paths: list[Path], stats: Counter) -> list[dict]:
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


# -------------------------------------------------------------------------------------------------------- merge
STUDIES = {'hartel2019': ('Hartel2019', read_hartel), 'geoghegan2015': ('Geoghegan2015', read_geoghegan),
           'sylvestersen2014': ('Sylvestersen2014', read_sylvestersen), 'lu2025': ('Lu2025', read_lu),
           'ma2023': ('Ma2023', read_ma)}


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
        notes = hartel_merge_note(grp) if study_id == 'Hartel2019' else [n for r in grp for n in r['note']]
        records.append({
            'study_id': study_id, 'reported_accession': '; '.join(accs),
            'reported_gene': next((r['gene'] for r in grp if r['gene']), ''), 'reported_position': pos,
            'peptide': marked, 'window': next((r['window'] for r in grp if r['window']), ''), 'methyl_degree': degree,
            'rme2_type': '; '.join(sorted(calls)) if degree == 'Rme2' else '',
            'cell_type': '; '.join(sorted({r['cell'] for r in grp})), 'enrichment': '; '.join(enrich),
            'localisation_probability': fmt_prob(max(locs) if locs else None), 'n_source_rows': len(grp),
            'source_tables': '; '.join(sorted({r['table'] for r in grp}, key=table_key)),
            'note': '; '.join(dict.fromkeys(notes)),
        })
    records.sort(key=lambda r: (r['reported_accession'].split('; ')[0], r['reported_position'], r['methyl_degree'],
                                r['peptide'], r['window']))
    return records


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--study', required=True, choices=sorted(STUDIES))
    ap.add_argument('--input', required=True, nargs='+', help='the study supplementary table(s), as in the usage above')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    study_id, reader = STUDIES[args.study]
    paths = [Path(p) for p in args.input]
    stats = Counter()
    rows = reader(paths, stats)
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

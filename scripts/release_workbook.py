"""Annotated Excel workbook of a methylarginine release (called by scripts/build_release.py).

Sheets: About, Sites, Source records, Changes, Studies, Reason codes, Columns. Every table is an Excel table (filters,
banded rows) with a note on each column header; the Sites sheet groups its columns under coloured bands, links each
accession to UniProt, shows the sequence windows in a fixed-width font with the methylated arginine in bold red, and
highlights sites pending a decision or resting on a paralog window.
"""
from __future__ import annotations

import pandas as pd
from openpyxl import Workbook
from openpyxl.cell.rich_text import CellRichText, TextBlock
from openpyxl.cell.text import InlineFont
from openpyxl.comments import Comment
from openpyxl.formatting.rule import FormulaRule
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.table import Table, TableStyleInfo

MONO = 'Consolas'
LINK = Font(color='0563C1', underline='single')
BOLD = Font(bold=True)
WRAP = Alignment(wrap_text=True, vertical='top')
GROUP_FILL = {
    'Site': 'DDEBF7', 'Sequence': 'E2EFDA', 'Methylation evidence': 'FCE4D6', 'Context': 'EDE7F6',
    'Provenance': 'FFF2CC', 'Status': 'F8CBAD', 'Change': 'E7E6E6',
}
FILL = {name: PatternFill('solid', start_color=color, end_color=color) for name, color in {
    'pending': 'F8CBAD', 'paralog': 'FFF2CC', 'good': 'C6EFCE', 'bad': 'FFC7CE', 'warn': 'FFEB9C', 'down': 'DDEBF7',
}.items()}
WIDTH = {
    'site_id': 15, 'gene': 11, 'protein_name': 34, 'uniprot_accession': 12, 'uniprot_entry_name': 15, 'position': 9,
    'residue': 8, 'window_15': 21, 'window_31': 40, 'protein_length': 9, 'reference_sequence_md5': 34,
    'motifs': 26, 'nearest_domain': 40, 'studies': 30, 'study_dois': 40, 'source_tables': 24, 'record_ids': 16,
    'release_status': 14, 'flags': 28, 'same_window_in': 30, 'maron2021_inhibitor_detail': 60,
    'reference_text': 40, 'evidence_sequence': 34, 'placement': 44, 'evidence_check': 26, 'accepted_exception': 30,
    'citation': 48, 'doi': 30, 'terms': 40, 'meaning': 100, 'description': 110, 'records_now_at_or_previously_at': 30,
}


def cell_value(value):
    if value is None or value is pd.NA:
        return None
    if isinstance(value, float) and pd.isna(value):
        return None
    if hasattr(value, 'item'):  # numpy scalars
        value = value.item()
    return value


def centred(window: str) -> CellRichText:
    mid = len(window) // 2
    plain, site = InlineFont(rFont=MONO), InlineFont(rFont=MONO, b=True, color='C00000')
    return CellRichText(TextBlock(plain, window[:mid]), TextBlock(site, window[mid]), TextBlock(plain, window[mid + 1:]))


def add_table(ws, name: str, frame: pd.DataFrame, notes: dict, header_row: int, freeze: str, style: str = 'TableStyleMedium2'):
    """Write frame as an Excel table whose header sits on header_row, with a note on each header."""
    for k, column in enumerate(frame.columns, start=1):
        cell = ws.cell(row=header_row, column=k, value=column)
        if column in notes:
            cell.comment = Comment(f'{column}: {notes[column]}', 'release', width=320, height=140)
    for values in frame.itertuples(index=False):
        ws.append([cell_value(v) for v in values])
    last = get_column_letter(len(frame.columns))
    table = Table(displayName=name, ref=f'A{header_row}:{last}{header_row + max(len(frame), 1)}')
    table.tableStyleInfo = TableStyleInfo(name=style, showRowStripes=True)
    ws.add_table(table)
    ws.freeze_panes = freeze
    for k, column in enumerate(frame.columns, start=1):
        values = frame[column].astype(str)
        typical = int(values.str.len().quantile(0.9)) if len(values) else 0
        width = WIDTH.get(column, min(28, typical + 2))
        ws.column_dimensions[get_column_letter(k)].width = max(width, min(len(column) + 3, 32))
    return {column: get_column_letter(k) for k, column in enumerate(frame.columns, start=1)}


def highlight(ws, ref: str, formula: str, fill: str) -> None:
    ws.conditional_formatting.add(ref, FormulaRule(formula=[formula], fill=FILL[fill], stopIfTrue=False))


def link_column(ws, letter: str, first: int, last: int, url) -> None:
    for row in range(first, last + 1):
        cell = ws[f'{letter}{row}']
        if cell.value:
            cell.hyperlink = url(str(cell.value))
            cell.font = LINK


def write(path, version: str, previous: str, counts: dict, build: dict, sites: pd.DataFrame, records: pd.DataFrame,
          changes: pd.DataFrame, studies: pd.DataFrame, reasons: pd.DataFrame, columns: dict, flags: dict,
          source_about: dict) -> None:
    wb = Workbook()
    about = wb.active
    about.title = 'About'
    uniprot = lambda acc: f'https://www.uniprot.org/uniprotkb/{acc}/entry'  # noqa: E731

    # ---- Sites -------------------------------------------------------------------------------------------------
    ws = wb.create_sheet('Sites')
    groups = [g for _, g, _ in columns['Sites']]
    start = 1
    for k in range(1, len(groups) + 1):
        if k == len(groups) or groups[k] != groups[k - 1]:
            cell = ws.cell(row=1, column=start, value=groups[k - 1])
            cell.font, cell.alignment = BOLD, Alignment(horizontal='center')
            fill = PatternFill('solid', start_color=GROUP_FILL[groups[k - 1]], end_color=GROUP_FILL[groups[k - 1]])
            for j in range(start, k + 1):
                ws.cell(row=1, column=j).fill = fill
            if k > start:
                ws.merge_cells(start_row=1, start_column=start, end_row=1, end_column=k)
            start = k + 1
    shown = sites.copy()
    for column in ('window_15', 'window_31'):
        shown[column] = [centred(w) for w in shown[column]]
    col = add_table(ws, 'Sites', shown, {c: t for c, _, t in columns['Sites']}, header_row=2, freeze='C3')
    first, last = 3, len(sites) + 2
    link_column(ws, col['uniprot_accession'], first, last, uniprot)
    whole = f'A{first}:{get_column_letter(len(sites.columns))}{last}'
    highlight(ws, whole, f'ISNUMBER(SEARCH("pending",${col["release_status"]}{first}))', 'pending')
    highlight(ws, whole, f'ISNUMBER(SEARCH("paralog_window",${col["flags"]}{first}))', 'paralog')
    highlight(ws, f'{col["change_since_previous_version"]}{first}:{col["change_since_previous_version"]}{last}',
              f'{col["change_since_previous_version"]}{first}="added"', 'good')
    for c in ('maron2021_GSK591_PRMT5i', 'maron2021_MS023_typeIi'):
        rng = f'{col[c]}{first}:{col[c]}{last}'
        highlight(ws, rng, f'{col[c]}{first}="decreased"', 'down')
        highlight(ws, rng, f'{col[c]}{first}="increased"', 'bad')
    ws.sheet_properties.tabColor = '2F5597'

    # ---- Source records ----------------------------------------------------------------------------------------
    ws = wb.create_sheet('Source records')
    col = add_table(ws, 'SourceRecords', records, dict(columns['Source records']), header_row=1, freeze='B2')
    first, last = 2, len(records) + 1
    for row in range(first, last + 1):
        ws[f'{col["evidence_sequence"]}{row}'].font = Font(name=MONO)
    rng = f'{col["change_since_previous_version"]}{first}:{col["change_since_previous_version"]}{last}'
    highlight(ws, rng, f'{col["change_since_previous_version"]}{first}="newly placed"', 'good')
    highlight(ws, rng, f'{col["change_since_previous_version"]}{first}="no longer placed"', 'bad')
    highlight(ws, rng, f'{col["change_since_previous_version"]}{first}="moved"', 'warn')

    # ---- Changes -----------------------------------------------------------------------------------------------
    ws = wb.create_sheet(f'Changes since {previous}')
    col = add_table(ws, 'Changes', changes, dict(columns['Changes']), header_row=1, freeze='B2')
    first, last = 2, len(changes) + 1
    link_column(ws, col['uniprot_accession'], first, last, uniprot)
    rng = f'{col["change"]}{first}:{col["change"]}{last}'
    highlight(ws, rng, f'{col["change"]}{first}="added"', 'good')
    highlight(ws, rng, f'{col["change"]}{first}="removed"', 'bad')

    # ---- Studies, reason codes, columns ------------------------------------------------------------------------
    ws = wb.create_sheet('Studies')
    col = add_table(ws, 'Studies', studies, dict(columns['Studies']), header_row=1, freeze='B2')
    link_column(ws, col['doi'], 2, len(studies) + 1, lambda doi: f'https://doi.org/{doi}')
    ws = wb.create_sheet('Reason codes')
    col = add_table(ws, 'ReasonCodes', reasons, dict(columns['Reason codes']), header_row=1, freeze='B2')
    for row in range(2, len(reasons) + 2):
        ws[f'{col["meaning"]}{row}'].alignment = WRAP
    dictionary = pd.DataFrame(
        [{'sheet': 'Sites', 'column': c, 'group': g, 'description': t} for c, g, t in columns['Sites']]
        + [{'sheet': f'Changes since {previous}' if sheet == 'Changes' else sheet, 'column': c, 'group': '', 'description': t}
           for sheet in ('Source records', 'Changes', 'Studies', 'Reason codes') for c, t in columns[sheet]])
    ws = wb.create_sheet('Columns')
    col = add_table(ws, 'Columns', dictionary, {'description': 'what the column holds'}, header_row=1, freeze='C2')
    for row in range(2, len(dictionary) + 2):
        ws[f'{col["description"]}{row}'].alignment = WRAP

    # ---- About -------------------------------------------------------------------------------------------------
    ws = about
    ws.column_dimensions['A'].width = 30
    ws.column_dimensions['B'].width = 120
    ws['A1'] = f'Human methylarginine sites, {version}'
    ws['A1'].font = Font(bold=True, size=16)
    ws['A2'] = ('The corrected union of human protein arginine methylation sites, with where each site comes from, '
                'its sequence, its methylation evidence and its structural context.')
    ws['A2'].alignment = Alignment(wrap_text=False)
    rows: list[tuple] = [
        (),
        ('Paper', 'DeAngelo, Silverstein & Shechter, Cell Chem Biol 33:1404-1421 (2026)', 'https://doi.org/10.1016/j.chembiol.2026.09.011'),
        ('Replaces', f'{previous}: the published union (repository tag v1.0.0), 13,719 sites on 5,090 accessions'),
        (),
        ('In this version',),
        ('Sites', f"{counts['sites_placed']:,} on {counts['proteins_placed']:,} reviewed human proteins"),
        ('Pending a decision', f"{counts['sites_pending']} sites (flag pending_splice_junction; highlighted orange on the Sites sheet)"),
        ('Source records', f"{counts['records']:,}, of which {counts['records_placed']:,} placed"),
        (f'Since {previous}', f"{counts['kept']:,} sites kept, {counts['removed']:,} removed, {counts['added']:,} added"),
        ('Methylation evidence', f"degree known for {counts['with_state']:,} sites; Maron 2021 antibody and inhibitor data for {counts['with_maron_ip']:,}"),
        (),
        ('Sheets',),
        ('Sites', 'one row per site: start here. Hover over a column header for its definition.', "#'Sites'!A1"),
        ('Source records', 'one row per record: what each study reported, where it was placed, and why that changed', "#'Source records'!A1"),
        (f'Changes since {previous}', 'every site of either version: kept, removed or added, with the reason codes', f"#'Changes since {previous}'!A1"),
        ('Studies', 'citation, DOI and contribution of each study or database', "#'Studies'!A1"),
        ('Reason codes', 'why a placement changed, with the issue numbers of DATA_PROVENANCE_REVIEW.md', "#'Reason codes'!A1"),
        ('Columns', 'definition of every column of every sheet', "#'Columns'!A1"),
        (),
        ('Source tables',),
        *[(name, text) for name, text in source_about.items()],
        ('Maron 2021 Tables S2, S3', 'methylation-degree calls, antibodies and PRMT inhibitor responses; they annotate sites but add none'),
        (),
        ('Flags',),
        *[(name, f"{counts['flags'].get(name, 0):,} sites: {text}") for name, text in flags.items()],
        (),
        ('Built from',),
        ('Code', f"PTM_Pipeline commit {build['commit']} (hpc/build_release.sh, scripts/build_release.py)"),
        ('Inputs', f"{build['inputs']}/ in the repository, checked against {build['inputs']}/MANIFEST.tsv"),
        ('Reference sequences', 'UniProtKB reviewed human canonical (UP000005640), staged 2026-04-07; release not recorded, earlier than 2026_03'),
        ('Context', 'MobiDB disorder (2026-03-30), InterPro domains and repeats (2026-04-07), motifs (scripts/06, 05, 08)'),
        (),
        ('Licence', 'CC BY-NC-SA 4.0. PhosphoSitePlus data (www.phosphosite.org): CC BY-NC-SA 3.0, non-commercial. Cite the studies of the sites you use.'),
    ]
    for r in rows:
        ws.append(list(r[:2]))
        row = ws.max_row
        if len(r) == 1:
            ws.cell(row=row, column=1).font = Font(bold=True, size=12)
        elif r:
            ws.cell(row=row, column=1).font = BOLD
            ws.cell(row=row, column=2).alignment = WRAP
            if len(r) == 3:
                ws.cell(row=row, column=2).hyperlink = r[2]
                ws.cell(row=row, column=2).font = LINK
    ws.sheet_properties.tabColor = '548235'
    wb.save(path)

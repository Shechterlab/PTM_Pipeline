"""The release site table as a spreadsheet (called by scripts/build_release.py).

One sheet holding the same rows and columns as the TSV, as an Excel table (filters, banded rows) with the header and
the first two columns frozen. Accessions link to UniProt; sequence windows are in a fixed-width font with the
methylated arginine in bold red. Column definitions are in the release README.
"""
from __future__ import annotations

import pandas as pd
from openpyxl import Workbook
from openpyxl.cell.rich_text import CellRichText, TextBlock
from openpyxl.cell.text import InlineFont
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.table import Table, TableStyleInfo

MONO = 'Consolas'
WIDTH = {
    'site_id': 15, 'gene': 11, 'protein_name': 34, 'uniprot_accession': 13, 'uniprot_entry_name': 15, 'position': 10,
    'residue': 9, 'window_15': 21, 'window_31': 40, 'protein_length': 10, 'reference_sequence_md5': 35,
    'motifs': 26, 'nearest_domain': 40, 'studies': 32, 'study_dois': 40, 'record_ids': 16, 'flags': 28,
    'same_window_in': 30, 'maron2021_inhibitor_detail': 60,
}


def cell_value(value):
    if value is None or value is pd.NA or (isinstance(value, float) and pd.isna(value)):
        return None
    return value.item() if hasattr(value, 'item') else value


def centred(window: str) -> CellRichText:
    mid = len(window) // 2
    plain, site = InlineFont(rFont=MONO), InlineFont(rFont=MONO, b=True, color='C00000')
    return CellRichText(TextBlock(plain, window[:mid]), TextBlock(site, window[mid]), TextBlock(plain, window[mid + 1:]))


def write(path, sites: pd.DataFrame) -> None:
    wb = Workbook()
    ws = wb.active
    ws.title = 'Sites'
    shown = sites.copy()
    for column in ('window_15', 'window_31'):
        shown[column] = [centred(w) for w in shown[column]]
    ws.append(list(shown.columns))
    for values in shown.itertuples(index=False):
        ws.append([cell_value(v) for v in values])
    last_col, last_row = get_column_letter(len(shown.columns)), len(shown) + 1
    table = Table(displayName='Sites', ref=f'A1:{last_col}{last_row}')
    table.tableStyleInfo = TableStyleInfo(name='TableStyleMedium2', showRowStripes=True)
    ws.add_table(table)
    ws.freeze_panes = 'C2'
    letter = {column: get_column_letter(k) for k, column in enumerate(shown.columns, start=1)}
    for row in range(2, last_row + 1):
        cell = ws[f"{letter['uniprot_accession']}{row}"]
        cell.hyperlink = f'https://www.uniprot.org/uniprotkb/{cell.value}/entry'
        cell.font = Font(color='0563C1', underline='single')
    for column, col in letter.items():
        typical = int(sites[column].astype(str).str.len().quantile(0.9)) if len(sites) else 0
        ws.column_dimensions[col].width = max(WIDTH.get(column, min(28, typical + 2)), min(len(column) + 3, 32))
    wb.save(path)

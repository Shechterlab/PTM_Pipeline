"""The release site table as a spreadsheet (called by scripts/build_release.py).

One sheet with the same rows and columns as the TSV, as an Excel table (filters, banded rows):

- the header is coloured by topic (site, sequence, methylation evidence, context, interactions, sources, status,
  change), and each header cell carries the column's definition as a note;
- the header and the first two columns stay in view;
- integer columns are numbers, accessions link to UniProt, and sequence windows are in a fixed-width font with the
  site's arginine in bold red;
- every other row is shaded grey and pending sites orange (by formula, so the shading follows sorting and filtering).

The file's internal times are fixed, so the same table always gives the same bytes.
"""
from __future__ import annotations

import datetime
import re
import zipfile
from pathlib import Path

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
FIXED_TIME = datetime.datetime(2026, 1, 1)
WIDTH = {  # characters, for columns whose content is wider than their name
    'site_id': 18,
    'protein_name': 38, 'window_15': 21, 'window_31': 42, 'reference_sequence_md5': 36, 'motifs': 26,
    'nearest_domain_name': 34, 'nearest_domain_class': 24, 'studies': 32, 'study_dois': 32, 'record_ids': 18,
    'flags': 28, 'same_window_accessions': 26, 'same_window_genes': 26, 'change_reasons': 24,
    'interaction_evidence': 26, 'source_tables': 24, 'gene': 14, 'uniprot_entry_name': 20,
}
# Header colour per topic: (first column of the topic, fill).
TOPICS = [
    ('site_id', '155E75'),                  # site
    ('window_15', '0E7490'),                # sequence
    ('methyl_states_reported', '9A3412'),   # methylation evidence
    ('motifs', '3F6212'),                   # context
    ('interaction_status', '6D28D9'),       # interactions
    ('n_studies', '1E40AF'),                # sources
    ('flags', '854D0E'),                    # status
    ('in_previous_version', '374151'),      # change
]


def centred(window: str) -> CellRichText:
    mid = len(window) // 2
    plain, site = InlineFont(rFont=MONO), InlineFont(rFont=MONO, b=True, color='C00000')
    return CellRichText(TextBlock(plain, window[:mid]), TextBlock(site, window[mid]), TextBlock(plain, window[mid + 1:]))


def write(path, sites: pd.DataFrame, columns: list[tuple] | None = None, title: str = 'Sites') -> None:
    """sites: the TSV as read with dtype=str; columns: (name, type, codes, description), default the site columns of
    scripts/release_codes.py."""
    if columns is None:
        import release_codes
        columns = release_codes.SITE_COLUMNS
    kinds = {name: kind for name, kind, *_ in columns}
    descriptions = {name: description for name, _k, _c, description in columns}
    wb = Workbook()
    ws = wb.active
    ws.title = title
    ws.append(list(sites.columns))
    for row in sites.itertuples(index=False):
        values = []
        for name, value in zip(sites.columns, row):
            if value == '':
                values.append(None)
            elif kinds.get(name) == 'integer':
                values.append(int(value))
            elif kinds.get(name) == 'number':
                values.append(float(value))
            elif name in ('window_15', 'window_31'):
                values.append(centred(value))
            else:
                values.append(value)
        ws.append(values)
    last_col, last_row = get_column_letter(len(sites.columns)), len(sites) + 1
    table = Table(displayName=title.replace(' ', '_').replace('.', '_'), ref=f'A1:{last_col}{last_row}')
    table.tableStyleInfo = TableStyleInfo(name='TableStyleLight1', showRowStripes=False)
    ws.add_table(table)
    ws.freeze_panes = 'C2'
    letter = {name: get_column_letter(k) for k, name in enumerate(sites.columns, start=1)}

    starts = {first: fill for first, fill in TOPICS}
    fill = TOPICS[0][1]
    for name in sites.columns:
        fill = starts.get(name, fill)
        cell = ws[f'{letter[name]}1']
        cell.font = Font(bold=True, color='FFFFFF')
        cell.fill = PatternFill('solid', fgColor=fill)
        cell.alignment = Alignment(wrap_text=True, vertical='top')
        if descriptions.get(name):
            cell.comment = Comment(f'{name}: {descriptions[name]}', 'release', width=320, height=140)
    ws.row_dimensions[1].height = 32

    link = Font(color='0563C1', underline='single')
    for row in range(2, last_row + 1):
        cell = ws[f"{letter['uniprot_accession']}{row}"]
        cell.hyperlink = f'https://www.uniprot.org/uniprotkb/{cell.value}/entry'
        cell.font = link
    # Row shading by formula, so it follows the rows when the table is sorted or filtered: pending sites in orange,
    # every other row in light grey.
    rows = f'A2:{last_col}{last_row}'
    if 'status' in letter:
        ws.conditional_formatting.add(rows, FormulaRule(formula=[f'${letter["status"]}2="pending"'], stopIfTrue=True,
                                                        fill=PatternFill('solid', fgColor='FDEBD3', bgColor='FDEBD3')))
    ws.conditional_formatting.add(rows, FormulaRule(formula=['MOD(ROW(),2)=0'],
                                                    fill=PatternFill('solid', fgColor='EEF3F5', bgColor='EEF3F5')))
    for name, col in letter.items():
        typical = int(sites[name].str.len().quantile(0.9)) if len(sites) else 0
        ws.column_dimensions[col].width = max(WIDTH.get(name, min(typical + 2, 30)), len(name) + 3)
    wb.properties.creator = 'PTM_Pipeline'
    wb.properties.title = title
    wb.properties.created = wb.properties.modified = FIXED_TIME
    wb.save(path)
    fix_zip_times(Path(path))


def fix_zip_times(path: Path) -> None:
    with zipfile.ZipFile(path) as z:
        members = [(info, z.read(info.filename)) for info in z.infolist()]
    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as z:
        for info, data in members:
            if info.filename == 'docProps/core.xml':  # openpyxl stamps the save time here
                data = re.sub(rb'(<dcterms:modified[^>]*>)[^<]*', rb'\g<1>' + FIXED_TIME.strftime('%Y-%m-%dT%H:%M:%SZ').encode(), data)
            info.date_time = FIXED_TIME.timetuple()[:6]
            z.writestr(info, data)

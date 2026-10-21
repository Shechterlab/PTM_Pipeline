#!/usr/bin/env python3
"""Stage the CARM1 knockout arginine-methylation response (Shishkova et al. 2017) on the canonical human proteome.

Shishkova, Zeng, Liu, Kwiecien, Hebert, Coon & Xu, "Global mapping of CARM1 substrates defines enzyme specificity
and substrate recognition", Nat Commun 8:15571 (2017), doi:10.1038/ncomms15571 (CC BY 4.0). Supplementary Data 1 lists the
ADMA-enriched peptides quantified in parental (WT) and CARM1 knockout (KO) MCF7 and MDA-MB-231 cells, one sheet per
cell line: the peptide with each methylated arginine in lower case, its start and stop on the UniProt entry named,
three WT and three KO intensities, the log2 KO/WT fold change and a p value.

Each peptide is placed on the canonical sequence of its entry's base accession: it must occur there exactly once
(isoform-numbered peptides are placed this way too). A methylated arginine is called down when a peptide carrying it
has log2 KO/WT <= -1 with p < 0.05 (the twofold decrease the paper uses; of the 329 peptides it highlights as
putative CARM1 substrates and the 335 that meet this rule, 328 are the same) and up for log2 KO/WT >= 1 with
p < 0.05. Sites carried by peptides that disagree are called mixed.

Outputs (in --outdir, default $PTM_RESULTS_ROOT/carm1_ko)
  carm1_ko_peptides.tsv     every peptide with its placement status
  carm1_ko_sites.tsv        one row per canonical site and cell line
  carm1_ko_summary.json     input checksum, rule, counts

The release build reads the staged table (inputs/layers/): scripts/build_release.py calls each site and
scripts/36_export_site_explorer_data.py carries the calls to the website.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
from collections import Counter
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook

from common import parse_fasta, save_table

STUDY = 'Shishkova 2017'
DOI = '10.1038/ncomms15571'
LOG2_CUTOFF, P_CUTOFF = 1.0, 0.05
RULE = f'peptide log2 KO/WT <= -{LOG2_CUTOFF:g} (down) or >= {LOG2_CUTOFF:g} (up) with p < {P_CUTOFF:g}'


def read_sheets(path: Path) -> pd.DataFrame:
    """One row per peptide and cell line, from the columns named in each sheet's header row."""
    wb = load_workbook(path, read_only=True, data_only=True)
    frames = []
    for ws in wb.worksheets:
        rows = ws.iter_rows(values_only=True)
        header = [str(h).strip().lower() if h is not None else '' for h in next(rows)]
        def col(prefix):
            for k, h in enumerate(header):
                if h.startswith(prefix):
                    return k
            raise SystemExit(f'{path} sheet {ws.title}: no column starting with {prefix!r} in {header}')
        pep, start, acc, gene, fold, pval = col('peptide'), col('start'), col('uniprot'), col('gene'), col('fold change'), col('pvalue')
        for row in rows:
            if not row or not row[pep]:
                continue
            frames.append({'cell_line': ws.title.strip(), 'peptide': str(row[pep]).strip(), 'reported_start': row[start],
                           'source_UniProtAC': str(row[acc]).strip(), 'gene': str(row[gene] or '').strip(),
                           'log2_ko_wt': row[fold], 'p': row[pval]})
    return pd.DataFrame(frames)


def place(row, canonical: dict[str, str]) -> tuple[str, int | None]:
    """(status, 1-based start on the canonical sequence) for one peptide."""
    acc = row.source_UniProtAC.split('-')[0]
    seq = canonical.get(acc)
    if seq is None:
        return 'accession_not_in_canonical_fasta', None
    plain = row.peptide.upper()
    hits = [k for k in range(len(seq) - len(plain) + 1) if seq.startswith(plain, k)]
    if not hits:
        return 'peptide_not_in_canonical', None
    if len(hits) > 1:
        return 'peptide_ambiguous', None
    start = hits[0] + 1
    if '-' not in row.source_UniProtAC and str(row.reported_start) == str(start):
        return 'position_verified', start
    return 'placed_by_sequence', start


def call(log2: float, p: float) -> str:
    if p < P_CUTOFF and log2 <= -LOG2_CUTOFF:
        return 'down'
    if p < P_CUTOFF and log2 >= LOG2_CUTOFF:
        return 'up'
    return 'unchanged'


def main() -> None:
    results_root = Path(os.environ.get('PTM_RESULTS_ROOT', 'PTM_results'))
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--xlsx', default=os.environ.get('PTM_CARM1_KO_XLSX', ''),
                    help='Shishkova 2017 Supplementary Data 1 (41467_2017_BFncomms15571_MOESM422_ESM.xlsx; PTM_CARM1_KO_XLSX).')
    ap.add_argument('--fasta', default=os.environ.get('PTM_CANONICAL_FASTA', ''),
                    help='UniProt reviewed canonical FASTA (PTM_CANONICAL_FASTA).')
    ap.add_argument('--outdir', default=str(results_root / 'carm1_ko'))
    args = ap.parse_args()
    if not args.xlsx or not args.fasta:
        raise SystemExit('--xlsx and --fasta are required (or set PTM_CARM1_KO_XLSX and PTM_CANONICAL_FASTA)')

    xlsx = Path(args.xlsx)
    peptides = read_sheets(xlsx)
    canonical = parse_fasta(args.fasta)
    placed = [place(row, canonical) for row in peptides.itertuples(index=False)]
    peptides['placement'] = [p[0] for p in placed]
    peptides['canonical_start'] = [p[1] for p in placed]
    peptides['canonical_UniProtAC'] = peptides['source_UniProtAC'].str.split('-').str[0]
    peptides['log2_ko_wt'] = pd.to_numeric(peptides['log2_ko_wt'], errors='coerce')
    peptides['p'] = pd.to_numeric(peptides['p'], errors='coerce')
    peptides['call'] = [call(f, p) if pd.notna(f) and pd.notna(p) else 'not_quantified' for f, p in zip(peptides['log2_ko_wt'], peptides['p'])]
    peptides['methyl_positions'] = [
        ';'.join(str(int(s) + k) for k, ch in enumerate(pep) if ch == 'r') if pd.notna(s) else ''
        for pep, s in zip(peptides['peptide'], peptides['canonical_start'])
    ]

    site_rows = []
    kept = peptides[peptides['canonical_start'].notna() & peptides['call'].ne('not_quantified')]
    for row in kept.itertuples(index=False):
        for pos in row.methyl_positions.split(';'):
            site_rows.append({'canonical_UniProtAC': row.canonical_UniProtAC, 'position': int(pos), 'cell_line': row.cell_line,
                              'log2': row.log2_ko_wt, 'p': row.p, 'call': row.call})
    per_site = pd.DataFrame(site_rows)
    sites = []
    for (acc, pos, cell), grp in per_site.groupby(['canonical_UniProtAC', 'position', 'cell_line'], sort=True):
        calls = Counter(grp['call'])
        site_call = 'mixed' if calls['down'] and calls['up'] else 'down' if calls['down'] else 'up' if calls['up'] else 'unchanged'
        seq = canonical[acc]
        sites.append({'canonical_UniProtAC': acc, 'position': pos, 'residue': seq[pos - 1], 'cell_line': cell,
                      'median_log2_ko_wt': round(float(grp['log2'].median()), 3), 'min_p': float(grp['p'].min()),
                      'peptides': len(grp), 'peptides_down': calls['down'], 'peptides_up': calls['up'], 'call': site_call})
    sites = pd.DataFrame(sites)
    not_r = sites[sites['residue'].ne('R')]
    if len(not_r):
        raise SystemExit(f'{len(not_r)} placed methylation sites are not arginines on the canonical sequence:\n{not_r.head()}')

    outdir = Path(args.outdir)
    save_table(peptides, outdir / 'carm1_ko_peptides.tsv')
    save_table(sites, outdir / 'carm1_ko_sites.tsv')
    summary = {
        'generated': dt.date.today().isoformat(),
        'study': STUDY, 'doi': DOI, 'rule': RULE,
        'input': str(xlsx), 'input_sha256': hashlib.sha256(xlsx.read_bytes()).hexdigest(),
        'canonical_fasta': str(Path(args.fasta)),
        'cell_lines': sorted(peptides['cell_line'].unique(), key=list(dict.fromkeys(peptides['cell_line'])).index),
        'peptides': dict(Counter(peptides['cell_line'])),
        'placement': dict(Counter(peptides['placement'])),
        'peptide_calls': {cell: dict(Counter(grp['call'])) for cell, grp in peptides.groupby('cell_line')},
        'sites': {cell: dict(Counter(grp['call'])) for cell, grp in sites.groupby('cell_line')},
        'unique_sites': int(sites[['canonical_UniProtAC', 'position']].drop_duplicates().shape[0]),
    }
    (outdir / 'carm1_ko_summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()

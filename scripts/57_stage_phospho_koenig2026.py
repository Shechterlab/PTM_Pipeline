#!/usr/bin/env python3
"""Stage the phosphosites and site stoichiometries of the Koenig et al. 2026 cell-line phosphoproteome atlas.

    python scripts/57_stage_phospho_koenig2026.py --outdir STAGE --downloaded YYYY-MM-DD [--cache DIR]
    python scripts/apply_staged.py STAGE

Koenig et al. 2026 (Nat Struct Mol Biol 33:1561) measured the phosphoproteomes of 33 human cell lines, each in four
replicates under control conditions and after serum stimulation (Stim), and deposited the tables on Zenodo (record
21469258, CC BY 4.0). Two of its archives are read here; the kinase-substrate tables (Table 7, from PhosphoSitePlus) are
not.

Phosphosites_and_ProteinGroups_Tables.zip, Table2_Phosphosites_Processed_Data.xlsx (one sheet, 107,112 rows): one row
per site and multiplicity (PTM.Multiplicity: the site on singly, doubly, ... phosphorylated peptides), with Protein,
Gene, PTM.SiteAA, PTM.SiteLocation, PTM.SiteProbability (every row is class I, >= 0.75), PTM.FlankingRegion (15
residues centred on the site, '_' past the ends) and, per cell line, four control and four Stim replicate log2
intensities (Plate<k>_<line>_Ctrl_<r>, _Stim_<r>). A site is kept where the reference (--fasta, UniProt 2026_03) has its
residue at its position and its flanking region matches there; a site on an isoform, or on an accession that is not a
reference entry but whose gene names a single reviewed entry (GN=), is moved there when the flanking region matches at
exactly one position. Its rows over multiplicities become one. Each row of layers/phospho_koenig2026_sites.tsv:
uniprot_accession, position, residue, site_probability (the highest over its rows), cell_lines_detected (cell lines
with a quantified replicate, control or Stim, in any of its rows), cell_lines (';'-joined), multiplicities,
flanking_region, source_accession, placed_by (position, isoform_window or gene_window).

Site_Stoichiometry.zip, Table4_SiteStoichiometry_Table.xlsx (sheet FullTableWithStoichiometry_Ctrl, 94,555 rows): per
cell line and phosphopeptide, the site (CollapseKey <accession>_<residue><position>), the peptide (PEP.StrippedSequence)
and the site's place in it (Position_pept, 1-based), the percentage of the site's copies phosphorylated under control
(StoichiometryCtrl) and Stim (StoichiometryStim), and Status: 'OK', or 'Break rule' / 'No change' where the calculation
gave values outside 0-100 (negative, or above 100), which are left out. A row is kept where the peptide places its site
residue at the position on the reference. Each row of layers/phospho_koenig2026_stoichiometry.tsv: uniprot_accession,
position, residue, cell_line, stoichiometry_control_pct, stoichiometry_stim_pct (percent of copies phosphorylated,
0-100), peptide (one row per site, cell line and peptide: 1,047 sites of a cell line are measured on more than one
peptide).

Writes the two tables, phospho_koenig2026_summary.json, manifest_rows.tsv and raw_files_rows.tsv.
"""
from __future__ import annotations

import argparse
import io
import json
import re
import shutil
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from urllib.request import Request, urlopen

import openpyxl
import pandas as pd
from common import (
    canonicalize_uniprot_accession,
    fasta_genes,
    md5_file,
    parse_fasta,
    peptide_places,
    sha256_file,
    unique_gene_accessions,
    window_matches,
    window_positions,
    write_stage_rows,
)

RECORD = 'https://zenodo.org/api/records/21469258/files/{}/content'
SITES_ZIP = 'Phosphosites_and_ProteinGroups_Tables.zip'
SITES_MEMBER = 'Phosphosites_and_ProteinGroups_Tables/Table2_Phosphosites_Processed_Data.xlsx'
STOICH_ZIP = 'Site_Stoichiometry.zip'
STOICH_MEMBER = 'Site_Stoichiometry/Table4_SiteStoichiometry_Table.xlsx'
SITES = 'phospho_koenig2026_sites.tsv'
STOICH = 'phospho_koenig2026_stoichiometry.tsv'
SUMMARY = 'phospho_koenig2026_summary.json'
CITATION = ('Koenig C et al., A phosphoproteome atlas of human cell lines reveals the landscape of kinase activity, '
            'Nat Struct Mol Biol 33:1561-1574 (2026), doi:10.1038/s41594-026-01877-6; data doi:10.5281/zenodo.21469258')
LICENCE = 'CC BY 4.0 (Zenodo record 21469258 and the article)'
REPLICATE = re.compile(r'^Plate\d+_(.+)_(Ctrl|Stim)_\d+$')


def place_window(seqs, gene_acc, acc, pos, residue, flank, gene) -> tuple[str, int, str]:
    seq = seqs.get(acc)
    if seq is not None:
        if not (1 <= pos <= len(seq)) or seq[pos - 1] != residue:
            return '', 0, 'residue differs'
        return (acc, pos, 'position') if window_matches(seq, pos, flank) else ('', 0, 'flanking region differs')
    base = canonicalize_uniprot_accession(acc)
    tried = [(base, 'isoform_window')] if base != acc and base in seqs else []
    if gene in gene_acc and gene_acc[gene] != base:
        tried.append((gene_acc[gene], 'gene_window'))
    for ref, how in tried:
        hits = [p for p in window_positions(seqs[ref], flank) if seqs[ref][p - 1] == residue]
        if len(hits) == 1:
            return ref, hits[0], how
    return '', 0, 'isoform or gene: no single canonical position' if tried else 'not in reference'


def read_sites(zpath: Path, seqs, gene_acc) -> tuple[pd.DataFrame, Counter, int]:
    with zipfile.ZipFile(zpath) as z:
        book = openpyxl.load_workbook(io.BytesIO(z.read(SITES_MEMBER)), read_only=True)
    rows = book.worksheets[0].iter_rows(values_only=True)
    head = next(rows)
    ix = {h: i for i, h in enumerate(head)}
    lines = defaultdict(list)
    for i, h in enumerate(head):
        m = REPLICATE.match(str(h))
        if m:
            lines[m.group(1)].append(i)
    site = {}
    left, n = Counter(), 0
    for r in rows:
        if r[ix['Protein']] is None:
            continue
        n += 1
        if r[ix['PTM.SiteLocation']] is None:
            left['no site position'] += 1
            continue
        acc, residue = str(r[ix['Protein']]).strip(), str(r[ix['PTM.SiteAA']]).strip()
        pos, flank = int(r[ix['PTM.SiteLocation']]), str(r[ix['PTM.FlankingRegion']] or '')
        gene = str(r[ix['Gene']] or '').split(';')[0].strip()
        ref, rpos, how = place_window(seqs, gene_acc, acc, pos, residue, flank, gene)
        if not ref:
            left[how] += 1
            continue
        seen = {cl for cl, cols in lines.items() if any(r[c] is not None for c in cols)}
        s = site.setdefault((ref, rpos), {'residue': residue, 'prob': 0.0, 'lines': set(), 'mult': set(), 'flank': flank,
                                          'src': set(), 'how': how})
        s['prob'] = max(s['prob'], float(r[ix['PTM.SiteProbability']]))
        s['lines'] |= seen
        s['mult'].add(str(r[ix['PTM.Multiplicity']]))
        s['src'].add(acc)
    df = pd.DataFrame([{'uniprot_accession': a, 'position': p, 'residue': s['residue'], 'site_probability': s['prob'],
                        'cell_lines_detected': len(s['lines']), 'cell_lines': ';'.join(sorted(s['lines'])),
                        'multiplicities': ';'.join(sorted(s['mult'])), 'flanking_region': s['flank'],
                        'source_accession': ';'.join(sorted(s['src'])), 'placed_by': s['how']}
                       for (a, p), s in site.items()])
    return df.sort_values(['uniprot_accession', 'position']).reset_index(drop=True), left, n


def read_stoichiometry(zpath: Path, seqs) -> tuple[pd.DataFrame, Counter, int]:
    with zipfile.ZipFile(zpath) as z:
        tab = pd.read_excel(io.BytesIO(z.read(STOICH_MEMBER)), sheet_name='FullTableWithStoichiometry_Ctrl')
    left, rows = Counter(), []
    for r in tab.to_dict('records'):
        if r['Status'] != 'OK':
            left[f'status {r["Status"]}'] += 1
            continue
        acc, residue, pos = str(r['PG.ProteinGroups']).strip(), str(r['AA']).strip(), int(r['PhosPos'])
        seq = seqs.get(acc)
        if seq is None:
            left['not in reference'] += 1
            continue
        if not (1 <= pos <= len(seq)) or seq[pos - 1] != residue:
            left['residue differs'] += 1
            continue
        if not peptide_places(seq, pos, [(str(r['PEP.StrippedSequence']), [int(r['Position_pept'])])]):
            left['peptide does not place it'] += 1
            continue
        rows.append({'uniprot_accession': acc, 'position': pos, 'residue': residue, 'cell_line': r['CellLine'],
                     'stoichiometry_control_pct': round(float(r['StoichiometryCtrl']), 4),
                     'stoichiometry_stim_pct': round(float(r['StoichiometryStim']), 4),
                     'peptide': r['PEP.StrippedSequence']})
    df = pd.DataFrame(rows).drop_duplicates(['uniprot_accession', 'position', 'cell_line', 'peptide'])
    return df.sort_values(['uniprot_accession', 'position', 'cell_line', 'peptide']).reset_index(drop=True), left, len(tab)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--fasta', default='inputs/reference/uniprot_human_reviewed_2026_03.fasta')
    ap.add_argument('--cache', type=Path, help='read the two archives from here instead of downloading them')
    ap.add_argument('--manifest', default='inputs/MANIFEST.tsv', help='read only: its columns')
    ap.add_argument('--raw-files', default='inputs/layers/raw_files.tsv', help='read only: its columns')
    ap.add_argument('--downloaded', required=True)
    ap.add_argument('--outdir', required=True, type=Path)
    a = ap.parse_args()
    out = a.outdir
    out.mkdir(parents=True, exist_ok=True)
    paths = {}
    for name in (SITES_ZIP, STOICH_ZIP):
        paths[name] = (a.cache / name) if a.cache else out / name
        if not a.cache:
            req = Request(RECORD.format(name), headers={'User-Agent': 'PTM-pipeline-stager/1.0'})
            with urlopen(req, timeout=3600) as r, open(paths[name], 'wb') as fh:
                shutil.copyfileobj(r, fh, 1 << 20)
    seqs = parse_fasta(a.fasta)
    gene_acc = unique_gene_accessions(fasta_genes(a.fasta))
    sites, left_sites, n_sites = read_sites(paths[SITES_ZIP], seqs, gene_acc)
    stoich, left_stoich, n_stoich = read_stoichiometry(paths[STOICH_ZIP], seqs)
    sites.to_csv(out / SITES, sep='\t', index=False)
    stoich.to_csv(out / STOICH, sep='\t', index=False)
    summary = {
        'script': 'scripts/57_stage_phospho_koenig2026.py', 'downloaded': a.downloaded,
        'files': {n: {'md5': md5_file(p), 'url': RECORD.format(n)} for n, p in paths.items()},
        'reference': {'fasta': a.fasta, 'md5': md5_file(a.fasta)},
        'sites': {'source_rows': n_sites, 'rows': len(sites), 'proteins': int(sites.uniprot_accession.nunique()),
                  'by_residue': {k: int(v) for k, v in sites.residue.value_counts().items()},
                  'cell_lines_detected_quantiles': {q: float(sites.cell_lines_detected.quantile(q)) for q in (0.25, 0.5, 0.75)},
                  'in_all_33_cell_lines': int((sites.cell_lines_detected == 33).sum()),
                  'placed_by': {k: int(v) for k, v in sites.placed_by.value_counts().items()}, 'left_out': dict(left_sites)},
        'stoichiometry': {'source_rows': n_stoich, 'rows': len(stoich),
                          'sites': int(len(stoich.drop_duplicates(['uniprot_accession', 'position']))),
                          'cell_lines': int(stoich.cell_line.nunique()),
                          'median_control_pct': float(stoich.stoichiometry_control_pct.median()),
                          'median_stim_pct': float(stoich.stoichiometry_stim_pct.median()), 'left_out': dict(left_stoich)},
    }
    (out / SUMMARY).write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    s, t = summary['sites'], summary['stoichiometry']
    derived = '; '.join(f'{n} (MD5 {summary["files"][n]["md5"]})' for n in paths)
    mrows = [
        {'file': f'layers/{SITES}', 'rows': len(sites),
         'content': (f'Class I phosphosites (localisation >= 0.75) of 33 human cell lines placed on the 2026_03 reference: '
                     f'{len(sites):,} sites on {s["proteins"]:,} proteins, each with the cell lines where it was quantified'),
         'derived_from': f'{SITES_MEMBER} in {derived.split("; ")[0]}',
         'how_made_or_checked': 'scripts/57_stage_phospho_koenig2026.py: residue at the position of the reference and the '
                                '15-residue flanking region matching there; isoform and gene sites moved only to a single '
                                'canonical position'},
        {'file': f'layers/{STOICH}', 'rows': len(stoich),
         'content': (f'Phosphorylation stoichiometry (percent of copies, control and serum-stimulated) per cell line: '
                     f'{len(stoich):,} rows, {t["sites"]:,} sites in {t["cell_lines"]} cell lines, rows of status OK only'),
         'derived_from': f'{STOICH_MEMBER} in {derived.split("; ")[1]}',
         'how_made_or_checked': 'scripts/57_stage_phospho_koenig2026.py: status OK; residue at the position of the '
                                'reference, where the row\'s peptide places it'},
        {'file': f'layers/{SUMMARY}', 'rows': '', 'derived_from': '',
         'content': 'Archives, MD5s, counts and what was left out of the Koenig 2026 phosphosite and stoichiometry staging',
         'how_made_or_checked': 'scripts/57_stage_phospho_koenig2026.py'},
    ]
    for r in mrows:
        p = out / r['file'].split('/')[-1]
        r.update({'bytes': p.stat().st_size, 'md5': md5_file(p), 'redistribution': 'git', 'citation': CITATION,
                  'source': 'Derived by scripts/57_stage_phospho_koenig2026.py from ' + '; '.join(RECORD.format(n) for n in paths),
                  'licence': LICENCE, 'version_or_date': f'Zenodo record 21469258 (2026-07-21), downloaded {a.downloaded}'})
    rrows = [{'file': n, 'url': RECORD.format(n), 'bytes': p.stat().st_size, 'md5': md5_file(p), 'sha256': sha256_file(p),
              'downloaded': a.downloaded, 'version': 'Zenodo record 21469258', 'kept': 'zenodo'} for n, p in paths.items()]
    write_stage_rows(out, mrows, rrows, a.manifest, a.raw_files)
    print(json.dumps({k: summary[k] for k in ('sites', 'stoichiometry')}, indent=1))


if __name__ == '__main__':
    main()

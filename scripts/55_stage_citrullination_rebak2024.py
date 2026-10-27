#!/usr/bin/env python3
"""Stage the human citrullination sites of Rebak et al. 2024 (Supplementary Table 3) on the reference.

    python scripts/55_stage_citrullination_rebak2024.py --outdir STAGE --downloaded YYYY-MM-DD [--cache DIR]
    python scripts/apply_staged.py STAGE

Rebak et al. 2024 (Nat Struct Mol Biol 31:977) mapped citrullination in HL-60 cells and HL-60-derived neutrophil-like
cells (NLC), with and without calcium ionophore, and in a second study with the PADI4 inhibitor GSK484 (1, 5 and 20 uM)
against its own control (NLC_Cal30_GSKstudy). Their workbook (Supplementary Tables 1-15, 173 MB) lists the 14,056 sites
in sheet Supp_Table3_Cit_Sites (two header rows; the second names the columns): UniProt accession, the position on the
canonical sequence (Position_canonical), the modified peptides (citrullines marked '(ci)'), localisation probability,
the mode of identification per condition (By MS/MS, By matching or empty) and the intensity of each replicate. The sheet
has no column of regulation or response to the inhibitor; the replicate intensities of the GSK484 study are kept so one
can be computed. The article is CC BY 4.0.

A site is kept where the reference (--fasta, UniProt 2026_03) has an arginine at its position and one of its peptides
places a '(ci)' arginine there; a site on an isoform accession, or whose accession is not a reference entry but whose
gene names a single reviewed entry (GN=), is moved there when its peptides place a citrulline at exactly one position.
The others are counted and left out.

Each row: uniprot_accession, position, residue (R), localization_probability, total_intensity, identified_in (the
conditions with an identification, ';'-joined, as the source names them: HL60, HL60_Cal15, NLC, NLC_Cal15, NLC_Cal30,
NLC_Cal30_GSKstudy, GSK_484_1uM, GSK_484_5uM, GSK_484_20uM), then the GSK484 study's replicate intensities
(gsk_control_r1-r4 = NLC_Cal30_GSKstudy, gsk484_1uM_r1-r4, gsk484_5uM_r1-r4, gsk484_20uM_r1-r4; empty = not
quantified), modified_peptides, source_accession, placed_by (position, isoform_peptide or gene_peptide).

Writes layers/citrullination_rebak2024_sites.tsv, citrullination_rebak2024_summary.json, manifest_rows.tsv and
raw_files_rows.tsv.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from urllib.request import Request, urlopen

import pandas as pd
from common import (
    canonicalize_uniprot_accession,
    fasta_genes,
    marked_peptide,
    md5_file,
    parse_fasta,
    peptide_places,
    peptide_positions,
    sha256_file,
    unique_gene_accessions,
    write_stage_rows,
)

URL = ('https://static-content.springer.com/esm/art%3A10.1038%2Fs41594-024-01214-9/MediaObjects/'
       '41594_2024_1214_MOESM3_ESM.xlsx')
RAW = '41594_2024_1214_MOESM3_ESM.xlsx'
SHEET = 'Supp_Table3_Cit_Sites'
OUT = 'citrullination_rebak2024_sites.tsv'
SUMMARY = 'citrullination_rebak2024_summary.json'
CITATION = ('Rebak AS et al., A quantitative and site-specific atlas of the citrullinome reveals widespread existence of '
            'citrullination and insights into PADI4 substrates, Nat Struct Mol Biol 31:977-995 (2024), '
            'doi:10.1038/s41594-024-01214-9, Supplementary Table 3')
LICENCE = 'CC BY 4.0 (article and supplementary tables)'
CONDITIONS = ('HL60', 'HL60_Cal15', 'NLC', 'NLC_Cal15', 'NLC_Cal30', 'NLC_Cal30_GSKstudy', 'GSK_484_1uM', 'GSK_484_5uM',
              'GSK_484_20uM')
GSK = {'NLC_Cal30_GSKstudy': 'gsk_control', 'GSK_484_1uM': 'gsk484_1uM', 'GSK_484_5uM': 'gsk484_5uM',
       'GSK_484_20uM': 'gsk484_20uM'}


def place(seqs, gene_acc, acc, pos, gene, peptides) -> tuple[str, int, str]:
    seq = seqs.get(acc)
    if seq is not None:
        if not (1 <= pos <= len(seq)) or seq[pos - 1] != 'R':
            return '', 0, 'residue differs'
        return (acc, pos, 'position') if peptide_places(seq, pos, peptides) else ('', 0, 'no peptide places it')
    base = canonicalize_uniprot_accession(acc)
    tried = [(base, 'isoform_peptide')] if base != acc and base in seqs else []
    if gene in gene_acc and gene_acc[gene] != base:
        tried.append((gene_acc[gene], 'gene_peptide'))
    for ref, how in tried:
        hits = sorted(p for p in peptide_positions(seqs[ref], peptides) if seqs[ref][p - 1] == 'R')
        if len(hits) == 1:
            return ref, hits[0], how
    return '', 0, 'isoform or gene: no single canonical position' if tried else 'not in reference'


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--fasta', default='inputs/reference/uniprot_human_reviewed_2026_03.fasta')
    ap.add_argument('--cache', type=Path, help='read the workbook from here instead of downloading it')
    ap.add_argument('--manifest', default='inputs/MANIFEST.tsv', help='read only: its columns')
    ap.add_argument('--raw-files', default='inputs/layers/raw_files.tsv', help='read only: its columns')
    ap.add_argument('--downloaded', required=True)
    ap.add_argument('--outdir', required=True, type=Path)
    a = ap.parse_args()
    out = a.outdir
    out.mkdir(parents=True, exist_ok=True)
    path = (a.cache / RAW) if a.cache else out / RAW
    if not a.cache:
        with urlopen(Request(URL, headers={'User-Agent': 'PTM-pipeline-stager/1.0'}), timeout=1800) as r:
            path.write_bytes(r.read())
    seqs = parse_fasta(a.fasta)
    gene_acc = unique_gene_accessions(fasta_genes(a.fasta))
    tab = pd.read_excel(path, sheet_name=SHEET, header=1, dtype=str, keep_default_na=False)
    missing = [c for c in ('Protein (Uniprot)', 'Gene names', 'Position_canonical', 'Modified sequences identified',
                           'Localization prob', *CONDITIONS) if c not in tab.columns]
    if missing:
        raise SystemExit(f'{SHEET}: columns not found: {missing}')
    rows, left = [], Counter()
    for r in tab.to_dict('records'):
        acc = r['Protein (Uniprot)'].split(';')[0].strip()
        try:
            pos = int(float(r['Position_canonical']))
        except ValueError:
            left['no position'] += 1
            continue
        peptides = [marked_peptide(p, 'ci') for p in r['Modified sequences identified'].split(';') if p.strip()]
        peptides = [p for p in peptides if p[1]]
        if not peptides:
            left['no (ci) peptide'] += 1
            continue
        ref, rpos, how = place(seqs, gene_acc, acc, pos, r['Gene names'].split(';')[0].strip(), peptides)
        if not ref:
            left[how] += 1
            continue
        row = {'uniprot_accession': ref, 'position': rpos, 'residue': 'R',
               'localization_probability': r['Localization prob'],
               'total_intensity': r['Total intensity across all conditions'],
               'identified_in': ';'.join(c for c in CONDITIONS if r[c].strip())}
        for src, name in GSK.items():
            for k in range(1, 5):
                v = r.get(f'{src}_R{k}', '').strip()
                row[f'{name}_r{k}'] = '' if v in ('', '0') else v
        row.update({'modified_peptides': r['Modified sequences identified'], 'source_accession': r['Protein (Uniprot)'],
                    'placed_by': how})
        rows.append(row)
    df = pd.DataFrame(rows)
    dup = int(df.duplicated(['uniprot_accession', 'position']).sum())
    df = df.drop_duplicates(['uniprot_accession', 'position']).sort_values(['uniprot_accession', 'position'])
    df.to_csv(out / OUT, sep='\t', index=False)
    summary = {
        'script': 'scripts/55_stage_citrullination_rebak2024.py', 'downloaded': a.downloaded,
        'source': {'file': RAW, 'url': URL, 'md5': md5_file(path), 'sheet': SHEET, 'rows': len(tab)},
        'reference': {'fasta': a.fasta, 'md5': md5_file(a.fasta)},
        'rows': len(df), 'proteins': int(df.uniprot_accession.nunique()),
        'identified_in': {c: int(sum(c in x.split(';') for x in df.identified_in)) for c in CONDITIONS},
        'placed_by': {k: int(v) for k, v in df.placed_by.value_counts().items()},
        'duplicates_after_placing': dup, 'left_out': dict(left),
    }
    (out / SUMMARY).write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    mrows = [
        {'file': f'layers/{OUT}', 'rows': len(df),
         'content': (f'Citrullinated arginines of HL-60 and neutrophil-like cells, with the replicate intensities of the '
                     f'PADI4 inhibitor (GSK484) study, placed on the 2026_03 reference: {len(df):,} sites on '
                     f'{summary["proteins"]:,} proteins'),
         'derived_from': f'{RAW} sheet {SHEET} (MD5 {summary["source"]["md5"]})',
         'how_made_or_checked': 'scripts/55_stage_citrullination_rebak2024.py: an arginine at the position of the '
                                'reference on which a (ci) peptide places; isoform and gene sites moved only to a single '
                                'canonical position'},
        {'file': f'layers/{SUMMARY}', 'rows': '', 'derived_from': '',
         'content': 'Source file, MD5, counts by condition and what was left out of the Rebak 2024 citrullination staging',
         'how_made_or_checked': 'scripts/55_stage_citrullination_rebak2024.py'},
    ]
    for r in mrows:
        p = out / r['file'].split('/')[-1]
        r.update({'bytes': p.stat().st_size, 'md5': md5_file(p), 'redistribution': 'git', 'citation': CITATION,
                  'source': f'Derived by scripts/55_stage_citrullination_rebak2024.py from {URL}', 'licence': LICENCE,
                  'version_or_date': f'Supplementary Tables of the published article, downloaded {a.downloaded}'})
    rrows = [{'file': RAW, 'url': URL, 'bytes': path.stat().st_size, 'md5': md5_file(path), 'sha256': sha256_file(path),
              'downloaded': a.downloaded, 'version': 'Rebak et al. 2024 Supplementary Tables 1-15', 'kept': 'zenodo'}]
    write_stage_rows(out, mrows, rrows, a.manifest, a.raw_files)
    print(json.dumps({k: summary[k] for k in ('rows', 'proteins', 'identified_in', 'placed_by', 'left_out')}, indent=1))


if __name__ == '__main__':
    main()

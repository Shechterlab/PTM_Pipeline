#!/usr/bin/env python3
"""Stage the human ADP-ribosylation sites of Hendriks et al. 2021 (Supplementary Data 2) on the reference.

    python scripts/56_stage_adpr_hendriks2021.py --outdir STAGE --downloaded YYYY-MM-DD [--cache DIR]
    python scripts/apply_staged.py STAGE

Hendriks et al. 2021 (Nat Commun 12:5893) enriched ADP-ribosylated peptides from human cells with the Af1521
macrodomain and with two pan-ADPr antibodies, and list the 2,758 unique sites they found in Supplementary Data 2 (sheet
Hendriks_DataS2; n = 4 cell culture replicates): UniProt accession, residue, position, a 17-residue sequence window, the
modified peptides, localisation probability, and replicate intensities per condition and enrichment, in columns named
<condition>_<enrichment>_<replicate> with conditions CTRL, ARH3 and PARG and enrichments AB1, AB2 (the antibodies) and
Af1521. The article is CC BY 4.0.

A site is kept where the reference (--fasta, UniProt 2026_03) has its residue at its position and one of its peptides
(the residue marked '(ADPr)') places it there, or, where the row lists no peptide, its window matches there. A site on
an isoform accession (P12345-2) is moved to the canonical sequence when its peptides or window place it at exactly one
position there; the others are counted and left out.

Each row: uniprot_accession, position, residue, localization_probability, ks_motif (the source's KS column: a lysine
before the serine), intensity, msms_count, enrichments (the enrichments with intensity in any replicate, ';'-joined),
conditions (likewise, by condition), modified_peptides, source_accession (as listed), placed_by (position, or
isoform_window / isoform_peptide when moved).

Writes layers/adpr_hendriks2021_sites.tsv, adpr_hendriks2021_summary.json, manifest_rows.tsv and raw_files_rows.tsv.
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
    marked_peptide,
    md5_file,
    parse_fasta,
    peptide_places,
    peptide_positions,
    sha256_file,
    window_matches,
    window_positions,
    write_stage_rows,
)

URL = ('https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-021-26172-4/MediaObjects/'
       '41467_2021_26172_MOESM5_ESM.xlsx')
RAW = '41467_2021_26172_MOESM5_ESM.xlsx'
SHEET = 'Hendriks_DataS2'
OUT = 'adpr_hendriks2021_sites.tsv'
SUMMARY = 'adpr_hendriks2021_summary.json'
CITATION = ('Hendriks IA et al., The regulatory landscape of the human HPF1- and ARH3-dependent ADP-ribosylome, '
            'Nat Commun 12:5893 (2021), doi:10.1038/s41467-021-26172-4, Supplementary Data 2')
LICENCE = 'CC BY 4.0 (article and supplementary data)'
CONDITIONS = ('CTRL', 'ARH3', 'PARG')
ENRICHMENTS = ('AB1', 'AB2', 'Af1521')


def place(seqs: dict[str, str], acc: str, pos: int, residue: str, window: str, peptides: list) -> tuple[str, int, str]:
    """(accession, position, how) of a site on the reference, or ('', 0, reason) when it does not place."""
    seq = seqs.get(acc)
    if seq is not None:
        if not (1 <= pos <= len(seq)) or seq[pos - 1] != residue:
            return '', 0, 'residue differs'
        if peptides:
            return (acc, pos, 'position') if peptide_places(seq, pos, peptides) else ('', 0, 'no peptide places it')
        return (acc, pos, 'position') if window_matches(seq, pos, window) else ('', 0, 'window differs')
    base = canonicalize_uniprot_accession(acc)
    seq = seqs.get(base)
    if base == acc or seq is None:
        return '', 0, 'not in reference'
    hits = sorted(p for p in (peptide_positions(seq, peptides) if peptides else set(window_positions(seq, window)))
                  if seq[p - 1] == residue)
    if len(hits) == 1:
        return base, hits[0], 'isoform_peptide' if peptides else 'isoform_window'
    return '', 0, 'isoform: no single canonical position'


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
        with urlopen(Request(URL, headers={'User-Agent': 'PTM-pipeline-stager/1.0'}), timeout=600) as r:
            path.write_bytes(r.read())
    seqs = parse_fasta(a.fasta)
    tab = pd.read_excel(path, sheet_name=SHEET, dtype=str, keep_default_na=False)
    groups = {(c, e): [col for col in tab.columns if col.startswith(f'{c}_{e}_')] for c in CONDITIONS for e in ENRICHMENTS}
    rows, left = [], Counter()
    for r in tab.to_dict('records'):
        acc, residue = r['Protein (Uniprot ID)'].strip(), r['Amino acid'].strip()
        try:
            pos = int(float(r['Position']))
        except ValueError:
            left['no position'] += 1
            continue
        peptides = [marked_peptide(p, 'ADPr') for p in r['Modified sequences'].split(';') if p.strip()]
        peptides = [p for p in peptides if p[1]]
        ref, rpos, how = place(seqs, acc, pos, residue, r['Sequence window'], peptides)
        if not ref:
            left[how] += 1
            continue
        seen = {k for k, cols in groups.items() if any(r[c] not in ('', '0') for c in cols)}
        rows.append({'uniprot_accession': ref, 'position': rpos, 'residue': residue,
                     'localization_probability': r['Localization probability'], 'ks_motif': r['KS'],
                     'intensity': r['Intensity'], 'msms_count': r['MS/MS count'],
                     'enrichments': ';'.join(e for e in ENRICHMENTS if any(k[1] == e for k in seen)),
                     'conditions': ';'.join(c for c in CONDITIONS if any(k[0] == c for k in seen)),
                     'modified_peptides': r['Modified sequences'], 'source_accession': acc, 'placed_by': how})
    df = pd.DataFrame(rows)
    dup = int(df.duplicated(['uniprot_accession', 'position']).sum())
    df = df.drop_duplicates(['uniprot_accession', 'position']).sort_values(['uniprot_accession', 'position'])
    df.to_csv(out / OUT, sep='\t', index=False)
    summary = {
        'script': 'scripts/56_stage_adpr_hendriks2021.py', 'downloaded': a.downloaded,
        'source': {'file': RAW, 'url': URL, 'md5': md5_file(path), 'sheet': SHEET, 'rows': len(tab)},
        'reference': {'fasta': a.fasta, 'md5': md5_file(a.fasta)},
        'rows': len(df), 'proteins': int(df.uniprot_accession.nunique()),
        'by_residue': {k: int(v) for k, v in df.residue.value_counts().items()},
        'placed_by': {k: int(v) for k, v in df.placed_by.value_counts().items()},
        'duplicates_after_placing': dup, 'left_out': dict(left),
    }
    (out / SUMMARY).write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    by = summary['by_residue']
    mrows = [
        {'file': f'layers/{OUT}', 'rows': len(df),
         'content': (f'ADP-ribosylation sites of human cells (Af1521 and two pan-ADPr antibodies) placed on the 2026_03 '
                     f'reference: {len(df):,} sites on {summary["proteins"]:,} proteins ('
                     + ', '.join(f'{v:,} {k}' for k, v in by.items()) + ')'),
         'derived_from': f'{RAW} (MD5 {summary["source"]["md5"]})',
         'how_made_or_checked': 'scripts/56_stage_adpr_hendriks2021.py: residue at the position of the reference; an '
                                '(ADPr) peptide places it there (window where no peptide); isoform sites moved only to '
                                'a single canonical position'},
        {'file': f'layers/{SUMMARY}', 'rows': '', 'derived_from': '',
         'content': 'Source file, MD5, counts by residue and what was left out of the Hendriks 2021 ADP-ribosylation staging',
         'how_made_or_checked': 'scripts/56_stage_adpr_hendriks2021.py'},
    ]
    for r in mrows:
        p = out / r['file'].split('/')[-1]
        r.update({'bytes': p.stat().st_size, 'md5': md5_file(p), 'redistribution': 'git', 'citation': CITATION,
                  'source': f'Derived by scripts/56_stage_adpr_hendriks2021.py from {URL}', 'licence': LICENCE,
                  'version_or_date': f'Supplementary Data 2 of the published article, downloaded {a.downloaded}'})
    rrows = [{'file': RAW, 'url': URL, 'bytes': path.stat().st_size, 'md5': md5_file(path), 'sha256': sha256_file(path),
              'downloaded': a.downloaded, 'version': 'Hendriks et al. 2021 Supplementary Data 2', 'kept': 'zenodo'}]
    write_stage_rows(out, mrows, rrows, a.manifest, a.raw_files)
    print(json.dumps({k: summary[k] for k in ('rows', 'proteins', 'by_residue', 'placed_by', 'left_out')}, indent=1))


if __name__ == '__main__':
    main()

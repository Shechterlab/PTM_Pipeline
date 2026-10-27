#!/usr/bin/env python3
"""Stage the lysine acetylation stoichiometries of Hansen et al. 2019 (Supplementary Data 1) on the reference.

    python scripts/58_stage_acetyl_stoichiometry_hansen2019.py --outdir STAGE --downloaded YYYY-MM-DD [--cache DIR]
    python scripts/apply_staged.py STAGE

Hansen et al. 2019 (Nat Commun 10:1055) measured the fraction of each acetylated lysine's copies that carry the acetyl
group in HeLa cells, by partial chemical acetylation with an isotope-labelled reagent. Supplementary Data 1 (one
workbook): sheet 1c lists 6,829 sites with their 'Stoichiometry' (a fraction, 0-1: the median is 0.000235, the paper's
0.02%), the leading protein and position, every protein and position the peptides match, and a 31-residue window; sheet
1a lists the acetylated peptides ('Modified sequence', lysines marked '(ac)') of each site by the same protein and
position, and writes '>10%' where a peptide's stoichiometry was above the range it could measure. The article is
CC BY 4.0.

A site is placed on the reference (--fasta, UniProt 2026_03) by the first protein of its list that is a reference
entry with a lysine at the listed position on which one of the site's peptides places its '(ac)' lysine; failing that,
on the canonical sequence of an isoform, or the single reviewed entry of its gene (GN=), where its peptides place an
acetyl-lysine at exactly one position. The others are counted and left out.

Each row: uniprot_accession, position, residue (K), stoichiometry (fraction), above_measurable_range (True where sheet
1a gives '>10%' for a peptide of the site: the stoichiometry is then a lower bound), peptides (';'-joined),
source_protein and source_position (the leading protein and position, as listed), placed_by (position,
isoform_peptide or gene_peptide), sequence_window.

Writes layers/acetyl_stoichiometry_hansen2019.tsv, acetyl_stoichiometry_hansen2019_summary.json, manifest_rows.tsv and
raw_files_rows.tsv.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
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

URL = ('https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-09024-0/MediaObjects/'
       '41467_2019_9024_MOESM4_ESM.xlsx')
RAW = '41467_2019_9024_MOESM4_ESM.xlsx'
OUT = 'acetyl_stoichiometry_hansen2019.tsv'
SUMMARY = 'acetyl_stoichiometry_hansen2019_summary.json'
CITATION = ('Hansen BK et al., Analysis of human acetylation stoichiometry defines mechanistic constraints on protein '
            'regulation, Nat Commun 10:1055 (2019), doi:10.1038/s41467-019-09024-0, Supplementary Data 1')
LICENCE = 'CC BY 4.0 (article and supplementary data)'


def place(seqs, gene_acc, proteins, positions, gene, peptides) -> tuple[str, int, str]:
    for acc, p in zip(proteins, positions):
        seq = seqs.get(acc)
        if seq is not None and p and 1 <= p <= len(seq) and seq[p - 1] == 'K' and peptide_places(seq, p, peptides):
            return acc, p, 'position'
    tried = []
    for acc in proteins:
        base = canonicalize_uniprot_accession(acc)
        if base != acc and base in seqs:
            tried.append((base, 'isoform_peptide'))
    if gene in gene_acc:
        tried.append((gene_acc[gene], 'gene_peptide'))
    for acc, how in tried:
        hits = sorted(p for p in peptide_positions(seqs[acc], peptides) if seqs[acc][p - 1] == 'K')
        if len(hits) == 1:
            return acc, hits[0], how
    return '', 0, ('no listed protein, isoform or gene places it' if tried or any(a in seqs for a in proteins)
                   else 'not in reference')


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
    gene_acc = unique_gene_accessions(fasta_genes(a.fasta))
    peps = pd.read_excel(path, sheet_name='Supplementary Data 1a', dtype=str, keep_default_na=False)
    sites = pd.read_excel(path, sheet_name='Supplementary Data 1c', dtype=str, keep_default_na=False)
    pep_of, above = defaultdict(set), set()
    for r in peps.to_dict('records'):
        key = (r['Protein'], r['Position within protein'])
        pep_of[key].add(r['Modified sequence'])
        if '%' in r['Median Stoichiometry']:
            above.add(key)
    rows, left = [], Counter()
    for r in sites.to_dict('records'):
        key = (r['Protein'], r['Position within Protein'])
        peptides = [marked_peptide(p, 'ac') for p in sorted(pep_of.get(key, ()))]
        peptides = [p for p in peptides if p[1]]
        if not peptides:
            left['no acetyl peptide in sheet 1a'] += 1
            continue
        proteins = [x.strip() for x in r['Proteins'].split(';')]
        positions = [int(x) if x.strip().isdigit() else 0 for x in r['Positions within Proteins'].split(';')]
        gene = r['Gene names'].split(';')[0].strip()
        acc, pos, how = place(seqs, gene_acc, proteins, positions, gene, peptides)
        if not acc:
            left[how] += 1
            continue
        rows.append({'uniprot_accession': acc, 'position': pos, 'residue': 'K', 'stoichiometry': r['Stoichiometry'],
                     'above_measurable_range': key in above, 'peptides': ';'.join(sorted(pep_of[key])),
                     'source_protein': r['Protein'], 'source_position': r['Position within Protein'], 'placed_by': how,
                     'sequence_window': r['Sequence Window']})
    df = pd.DataFrame(rows)
    dup = int(df.duplicated(['uniprot_accession', 'position']).sum())
    df = df.sort_values(['uniprot_accession', 'position', 'stoichiometry']).drop_duplicates(
        ['uniprot_accession', 'position'], keep='last')
    df.to_csv(out / OUT, sep='\t', index=False)
    st = df.stoichiometry.astype(float)
    summary = {
        'script': 'scripts/58_stage_acetyl_stoichiometry_hansen2019.py', 'downloaded': a.downloaded,
        'source': {'file': RAW, 'url': URL, 'md5': md5_file(path), 'sites_sheet_rows': len(sites), 'peptide_sheet_rows': len(peps)},
        'reference': {'fasta': a.fasta, 'md5': md5_file(a.fasta)},
        'rows': len(df), 'proteins': int(df.uniprot_accession.nunique()),
        'stoichiometry_fraction': {'median': float(st.median()), 'q90': float(st.quantile(0.9)), 'max': float(st.max()),
                                   'above_measurable_range': int(df.above_measurable_range.sum())},
        'placed_by': {k: int(v) for k, v in df.placed_by.value_counts().items()},
        'duplicates_after_placing': dup, 'left_out': dict(left),
    }
    (out / SUMMARY).write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    mrows = [
        {'file': f'layers/{OUT}', 'rows': len(df),
         'content': (f'Lysine acetylation stoichiometry in HeLa cells (fraction of copies acetylated) placed on the '
                     f'2026_03 reference: {len(df):,} sites on {summary["proteins"]:,} proteins, median '
                     f'{summary["stoichiometry_fraction"]["median"]:.6f}; {summary["stoichiometry_fraction"]["above_measurable_range"]} '
                     f'above the measurable range (lower bounds)'),
         'derived_from': f'{RAW} (MD5 {summary["source"]["md5"]})',
         'how_made_or_checked': 'scripts/58_stage_acetyl_stoichiometry_hansen2019.py: a lysine at the position of the '
                                'reference on which an (ac) peptide of the site places; isoform and gene sites moved only to '
                                'a single canonical position'},
        {'file': f'layers/{SUMMARY}', 'rows': '', 'derived_from': '',
         'content': 'Source file, MD5, stoichiometry quantiles and what was left out of the Hansen 2019 staging',
         'how_made_or_checked': 'scripts/58_stage_acetyl_stoichiometry_hansen2019.py'},
    ]
    for r in mrows:
        p = out / r['file'].split('/')[-1]
        r.update({'bytes': p.stat().st_size, 'md5': md5_file(p), 'redistribution': 'git', 'citation': CITATION,
                  'source': f'Derived by scripts/58_stage_acetyl_stoichiometry_hansen2019.py from {URL}', 'licence': LICENCE,
                  'version_or_date': f'Supplementary Data 1 of the published article, downloaded {a.downloaded}'})
    rrows = [{'file': RAW, 'url': URL, 'bytes': path.stat().st_size, 'md5': md5_file(path), 'sha256': sha256_file(path),
              'downloaded': a.downloaded, 'version': 'Hansen et al. 2019 Supplementary Data 1', 'kept': 'zenodo'}]
    write_stage_rows(out, mrows, rrows, a.manifest, a.raw_files)
    print(json.dumps({k: summary[k] for k in ('rows', 'proteins', 'stoichiometry_fraction', 'placed_by', 'left_out')}, indent=1))


if __name__ == '__main__':
    main()

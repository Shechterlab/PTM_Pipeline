"""Peptide start positions from Maron et al. 2021 Tables S2 and S3, in the numbering Maron's Table S5 uses.

Table S5 gives each site's position and peptide but not where the peptide starts, so a peptide with several arginines
does not say which one it reports when the numbering differs from the reference entry. Tables S2 and S3 give the start
(Begin.Pos) of the same peptides; scripts/02 and scripts/qc_remap_evidence_check.py use it to count the reported
position from the peptide's start (--peptide-coordinates).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from common import canonicalize_uniprot_accession, normalize_accession, save_table

S3_SHEETS = ('Control_vs_GSK591', 'Control_vs_MS023')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--table-s2', required=True)
    ap.add_argument('--table-s3', required=True)
    ap.add_argument('--outdir', default='results/parsed_maron')
    args = ap.parse_args()

    frames = [pd.read_excel(args.table_s2, sheet_name=0).assign(source='Maron2021_TableS2')]
    frames += [pd.read_excel(args.table_s3, sheet_name=sheet).assign(source='Maron2021_TableS3') for sheet in S3_SHEETS]
    out = pd.concat([pd.DataFrame({
        'source': df['source'],
        'accession': df['Accession'].astype(str).str.strip(),
        'gene': df['Gene.Symbol'].astype(str).str.strip(),
        'peptide': df['Sequence'].astype(str).str.strip().str.upper(),
        'peptide_start': pd.to_numeric(df['Begin.Pos'], errors='coerce'),
    }) for df in frames], ignore_index=True)
    out['canonical_UniProtAC'] = out['accession'].map(normalize_accession).map(canonicalize_uniprot_accession)
    out = out[out['peptide_start'].notna() & out['peptide'].str.fullmatch(r'[A-Z]+')].copy()
    out['peptide_start'] = out['peptide_start'].astype(int)
    out = out.drop_duplicates(['canonical_UniProtAC', 'peptide', 'peptide_start', 'source']).sort_values(
        ['canonical_UniProtAC', 'peptide_start', 'peptide', 'source'])

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(out, outdir / 'maron_peptide_coordinates.tsv')
    print({'peptides': len(out), 'proteins': out['canonical_UniProtAC'].nunique()})


if __name__ == '__main__':
    main()

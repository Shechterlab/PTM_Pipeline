from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from common import (
    accession_is_isoform,
    canonical_gene_name,
    canonicalize_uniprot_accession,
    normalize_accession,
    save_table,
)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True)
    ap.add_argument('--outdir', default='results/parsed_maron')
    ap.add_argument('--exclude-pattern', default='Larsen')
    args = ap.parse_args()

    df = pd.read_excel(args.input, sheet_name=0)
    out = pd.DataFrame()
    out['substrate_UniProtAC_original'] = df['Accession'].astype(str).str.strip()
    out['substrate_UniProtAC'] = out['substrate_UniProtAC_original'].map(normalize_accession)
    out['canonical_UniProtAC'] = out['substrate_UniProtAC'].map(canonicalize_uniprot_accession)
    out['accession_is_isoform'] = out['substrate_UniProtAC'].map(accession_is_isoform)
    out['substrate_genename'] = df['Gene Symbol'].map(canonical_gene_name)
    out['position'] = pd.to_numeric(df['Position In Protein'], errors='coerce')
    out['residue'] = 'R'
    out['site'] = 'R' + out['position'].astype('Int64').astype(str)
    out['peptide_sequence'] = df['Peptide Sequence'].astype(str)
    out['reference'] = df['Reference'].astype(str)
    out['source'] = 'Maron2021_TableS5'
    out['source_family'] = 'Maron2021'
    out['ptm_type'] = 'METHYLATION'
    out['ptm_group'] = 'Arg methylation'
    out['organism'] = 'Homo sapiens (Human)'
    out['site_key'] = out['substrate_UniProtAC'] + ':' + out['site']
    out['site_key_canonical_naive'] = out['canonical_UniProtAC'] + ':' + out['site']
    out['is_arg_methyl'] = True
    out['data_layer'] = 'external_methyl_table'
    out['arg_methyl_state'] = ''
    out['arg_methyl_state_known'] = False
    out = out[out['position'].notna() & out['substrate_UniProtAC'].ne('')].copy()
    if args.exclude_pattern:
        out = out[~out['reference'].str.contains(args.exclude_pattern, case=False, na=False)].copy()
    out['position'] = out['position'].astype(int)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(out, outdir / 'maron_s5_arg_methyl_sites.tsv')
    summary = pd.DataFrame([{
        'rows': len(out),
        'proteins': out['substrate_UniProtAC'].nunique(),
        'unique_sites': out['site_key'].nunique(),
    }])
    save_table(summary, outdir / 'summary.tsv')
    print(summary.to_dict(orient='records')[0])


if __name__ == '__main__':
    main()

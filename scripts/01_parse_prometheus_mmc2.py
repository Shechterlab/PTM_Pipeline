from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from common import (
    accession_is_isoform,
    canonicalize_uniprot_accession,
    clean_peptide_field,
    normalize_accession,
    normalize_arg_methyl_state,
    save_table,
)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True)
    ap.add_argument('--sheet', default='RK-methylsites')
    ap.add_argument('--outdir', default='results/parsed_prometheus')
    args = ap.parse_args()

    df = pd.read_excel(args.input, sheet_name=args.sheet, header=6)
    df = df.rename(columns={
        'LeadProt': 'substrate_UniProtAC',
        'GeneName': 'substrate_genename',
        'RES': 'residue',
        'POS': 'position',
        'MOD': 'mod_state',
        'SeqWindow': 'seq_window',
        'Peptides': 'peptides',
    })
    out = df[df['residue'].astype(str).str.upper().eq('R')].copy()
    out['substrate_UniProtAC_original'] = out['substrate_UniProtAC'].astype(str).str.strip()
    out['substrate_UniProtAC'] = out['substrate_UniProtAC_original'].map(normalize_accession)
    out['canonical_UniProtAC'] = out['substrate_UniProtAC'].map(canonicalize_uniprot_accession)
    out['accession_is_isoform'] = out['substrate_UniProtAC'].map(accession_is_isoform)
    out['substrate_genename'] = out['substrate_genename'].astype(str).str.strip()
    out['position'] = pd.to_numeric(out['position'], errors='coerce')
    out = out[out['position'].notna() & out['substrate_UniProtAC'].ne('')].copy()
    out['position'] = out['position'].astype(int)
    out['residue'] = 'R'
    out['site'] = 'R' + out['position'].astype(str)
    out['seq_window'] = out['seq_window'].astype(str).str.strip()
    out['peptides'] = out['peptides'].map(clean_peptide_field)
    out['arg_methyl_state'] = out['mod_state'].map(normalize_arg_methyl_state)
    out['arg_methyl_state_known'] = out['arg_methyl_state'].astype(str).ne('')
    out['source'] = 'ProMetheusDB_mmc2'
    out['source_family'] = 'ProMetheusDB'
    out['ptm_type'] = 'METHYLATION'
    out['ptm_group'] = 'Arg methylation'
    out['organism'] = 'Homo sapiens (Human)'
    out['site_key'] = out['substrate_UniProtAC'] + ':' + out['site']
    out['site_key_canonical_naive'] = out['canonical_UniProtAC'] + ':' + out['site']
    out['is_arg_methyl'] = True
    out['data_layer'] = 'external_methyl_table'

    keep = [
        'substrate_UniProtAC_original', 'substrate_UniProtAC', 'canonical_UniProtAC', 'accession_is_isoform',
        'substrate_genename', 'residue', 'position', 'site', 'mod_state', 'arg_methyl_state',
        'arg_methyl_state_known', 'seq_window', 'peptides', 'source', 'source_family', 'ptm_type',
        'ptm_group', 'organism', 'site_key', 'site_key_canonical_naive', 'is_arg_methyl', 'data_layer'
    ]
    out = out[keep].copy()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(out, outdir / 'prometheus_mmc2_arg_methyl_sites.tsv')
    summary = pd.DataFrame([{
        'rows': len(out),
        'proteins': out['substrate_UniProtAC'].nunique(),
        'unique_sites': out['site_key'].nunique(),
    }])
    save_table(summary, outdir / 'summary.tsv')
    print(summary.to_dict(orient='records')[0])


if __name__ == '__main__':
    main()

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd

from common import clean_aa_sequence, normalize_arg_methyl_state, save_table



def inferred_residue_from_window(window: str) -> str:
    seq = clean_aa_sequence(window)
    if not seq:
        return ''
    midpoint = len(seq) // 2
    if 0 <= midpoint < len(seq):
        return seq[midpoint]
    return ''



def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True)
    ap.add_argument('--outdir', default='results/parsed_dbptm')
    args = ap.parse_args()

    path = Path(args.input)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8', errors='replace') as fh:
        df = pd.read_csv(fh, sep='\t', header=None)
    cols = [f'col{i}' for i in range(df.shape[1])]
    df.columns = cols

    out = pd.DataFrame()
    out['substrate_UniProtAC'] = df['col0'].astype(str).str.strip()
    out['position'] = pd.to_numeric(df['col1'], errors='coerce')
    out['ptm_annotation'] = df['col2'].astype(str)
    out['seq_window'] = df['col3'].astype(str) if df.shape[1] >= 4 else ''
    out['source'] = 'dbPTM_experimental'
    out['source_family'] = 'dbPTM'
    out['data_source_label'] = 'dbPTM_experimental'
    out['position'] = out['position'].astype('Int64')
    out = out[out['position'].notna()].copy()
    out['position'] = out['position'].astype(int)
    out['residue'] = out['seq_window'].map(inferred_residue_from_window)
    annotation_text = out['ptm_annotation'].astype(str).str.lower()
    keep = (out['residue'] == 'R') | annotation_text.str.contains('argin', na=False) | annotation_text.str.contains('rme', na=False)
    out = out[keep].copy()
    out['residue'] = 'R'
    out['site'] = out['residue'] + out['position'].astype(str)
    out['site_key'] = out['substrate_UniProtAC'] + ':' + out['site']
    out['ptm_group'] = 'Arg methylation'
    out['ptm_type'] = 'METHYLATION'
    out['mod_state'] = out['ptm_annotation'].map(normalize_arg_methyl_state)
    out['arg_methyl_state'] = out['mod_state']
    out['organism'] = 'unknown_from_file'
    out['raw_columns'] = df.astype(str).agg('|'.join, axis=1)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(out, outdir / 'dbptm_arg_methyl_sites.tsv')
    save_table(out, outdir / 'dbptm_parsed_preview.tsv')
    print({'rows_arg_methyl': len(out), 'columns_in_input': df.shape[1]})


if __name__ == '__main__':
    main()

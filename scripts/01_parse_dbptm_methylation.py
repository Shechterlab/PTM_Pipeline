from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd

from common import save_table


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True)
    ap.add_argument('--outdir', default='results/parsed_dbptm')
    args = ap.parse_args()

    path = Path(args.input)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8', errors='replace') as fh:
        df = pd.read_csv(fh, sep='\t', header=None)
    # Expected dbPTM experimental export per site: UniProt AC/ID, position, PTM type, window sequence, optional PMID columns.
    # This parser keeps the first four columns robustly and preserves the remainder for inspection.
    cols = [f'col{i}' for i in range(df.shape[1])]
    df.columns = cols
    out = pd.DataFrame()
    out['substrate_UniProtAC'] = df['col0'].astype(str).str.strip()
    out['position'] = pd.to_numeric(df['col1'], errors='coerce')
    out['ptm_annotation'] = df['col2'].astype(str)
    out['seq_window'] = df['col3'].astype(str) if df.shape[1] >= 4 else ''
    out['source'] = 'dbPTM_experimental'
    out['source_family'] = 'dbPTM'
    out['position'] = out['position'].astype('Int64')
    out = out[out['position'].notna()].copy()
    out['position'] = out['position'].astype(int)
    out['residue'] = out['seq_window'].str[10:11].replace('', 'R')
    out['site'] = out['residue'] + out['position'].astype(str)
    out['site_key'] = out['substrate_UniProtAC'] + ':' + out['site']
    out['ptm_type'] = 'METHYLATION'
    out['organism'] = 'unknown_from_file'
    out['raw_columns'] = df.astype(str).agg('|'.join, axis=1)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(out, outdir / 'dbptm_parsed_preview.tsv')
    print({'rows': len(out), 'columns_in_input': df.shape[1]})


if __name__ == '__main__':
    main()

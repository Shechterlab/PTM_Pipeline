#!/usr/bin/env python3
"""Move a reviewed staging folder into inputs/.

    python scripts/apply_staged.py STAGE_DIR [STAGE_DIR ...] [--inputs inputs]

A staging script run for review (into a scratch folder, without --manifest) writes its tables there with
manifest_rows.tsv (their rows for inputs/MANIFEST.tsv) and raw_files_rows.tsv (their downloads, for
inputs/layers/raw_files.tsv). Once the tables are reviewed, this copies each table the manifest rows mark 'git' to its
path under inputs/ (its MD5 must equal the row's), checks the MD5 of any other listed file the folder holds, and adds
the rows to inputs/MANIFEST.tsv and inputs/layers/raw_files.tsv (common.update_manifest, common.upsert_tsv_rows), so the
tables in inputs/ are the bytes that were reviewed. Files kept out of Git ('zenodo', 'checksum_only') are not copied.
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pandas as pd

from common import md5_file, refresh_manifest_row, update_manifest, upsert_tsv_rows


def apply(stage: Path, inputs: Path) -> list[str]:
    rows = pd.read_csv(stage / 'manifest_rows.tsv', sep='\t', dtype=str, keep_default_na=False).to_dict('records')
    done = []
    for r in rows:
        src = stage / Path(r['file']).name
        if r['redistribution'] == 'git':
            if not src.exists():
                raise SystemExit(f'{stage}: {src.name} is listed for Git but not in the folder')
            if md5_file(src) != r['md5']:
                raise SystemExit(f'{src}: MD5 differs from its manifest row (stage it again)')
            dest = inputs / r['file']
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, dest)
            done.append(r['file'])
        elif src.exists() and md5_file(src) != r['md5']:
            raise SystemExit(f'{src}: MD5 differs from its manifest row')
    update_manifest(inputs / 'MANIFEST.tsv', rows)
    raw = stage / 'raw_files_rows.tsv'
    if raw.exists():
        target = inputs / 'layers' / 'raw_files.tsv'
        upsert_tsv_rows(target, pd.read_csv(raw, sep='\t', dtype=str, keep_default_na=False).to_dict('records'))
        refresh_manifest_row(inputs / 'MANIFEST.tsv', target, 'layers/raw_files.tsv')
    return done


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('stage', nargs='+', help='reviewed staging folders')
    ap.add_argument('--inputs', default='inputs')
    args = ap.parse_args()
    for stage in args.stage:
        for f in apply(Path(stage), Path(args.inputs)):
            print(f'{stage}: {f}')


if __name__ == '__main__':
    main()

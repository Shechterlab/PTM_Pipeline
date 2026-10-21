"""Download the raw files behind inputs/layers/ and check each against the MD5 the release used.

    python scripts/00_fetch_static_sources.py --outdir PTM_data/external [--only hpa_subcellular_zip ...]

Most of these URLs are not versioned: they serve the source's current file (the Human Protein Atlas its latest
release), so a later download can differ from the one inputs/layers/ was staged from. Every download is compared with
its row in --checksums (inputs/layers/raw_files.tsv); a file that differs is kept, reported, and the script exits 1.
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import pandas as pd
import requests

URLS = {
    # CPLM 4.0 lysine methylation (Zhang et al., Nucleic Acids Res 2022); stage with 39_stage_cplm4_lys_methylation.py.
    'cplm4_methylation_zip': 'https://cplm.biocuckoo.cn/Download/Methylation.zip',
    # CARM1 knockout ADMA peptides (Shishkova et al., Nat Commun 2017, Supplementary Data 1, CC BY 4.0); 41_stage_carm1_ko_sites.py.
    'carm1_ko_xlsx': 'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fncomms15571/MediaObjects/41467_2017_BFncomms15571_MOESM422_ESM.xlsx',
    # Human Protein Atlas subcellular locations and its UniProt mapping (CC BY 4.0); 42_stage_hpa_subcellular.py.
    'hpa_subcellular_zip': 'https://www.proteinatlas.org/download/tsv/subcellular_location.tsv.zip',
    'hpa_atlas_zip': 'https://www.proteinatlas.org/download/proteinatlas.tsv.zip',
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir', default='PTM_data/external')
    ap.add_argument('--only', nargs='*', default=[])
    ap.add_argument('--checksums', default=str(Path(__file__).resolve().parents[1] / 'inputs' / 'layers' / 'raw_files.tsv'))
    args = ap.parse_args()

    known = pd.read_csv(args.checksums, sep='\t', dtype=str, keep_default_na=False).set_index('url')
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    names = args.only if args.only else list(URLS)
    differs = []
    for name in names:
        url = URLS[name]
        dest = outdir / Path(url).name
        print(f'downloading {name} -> {dest}')
        md5 = hashlib.md5()
        with requests.get(url, stream=True, timeout=120) as r:
            r.raise_for_status()
            with open(dest, 'wb') as fh:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        fh.write(chunk)
                        md5.update(chunk)
        if url not in known.index:
            print(f'  MD5 {md5.hexdigest()} (not in {Path(args.checksums).name})')
        elif md5.hexdigest() != known.at[url, 'md5']:
            differs.append(name)
            print(f'  MD5 {md5.hexdigest()} differs from the file of {known.at[url, "downloaded"]} '
                  f'({known.at[url, "version"]}, MD5 {known.at[url, "md5"]}): the source has changed it')
        else:
            print(f'  MD5 matches the file of {known.at[url, "downloaded"]} ({known.at[url, "version"]})')
    if differs:
        sys.exit(f'{len(differs)} downloads differ from the files inputs/layers/ was staged from: {differs}')


if __name__ == '__main__':
    main()

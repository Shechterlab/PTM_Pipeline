from __future__ import annotations

import argparse
from pathlib import Path

import requests

URLS = {
    'dbptm_methylation_gz': 'https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Methylation.gz',
    'dbptm_phosphorylation_gz': 'https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Phosphorylation.gz',
    'dbptm_acetylation_gz': 'https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Acetylation.gz',
    'dbptm_ubiquitination_gz': 'https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Ubiquitination.gz',
    'dbptm_sumoylation_gz': 'https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Sumoylation.gz',
    'dbptm_citrullination_gz': 'https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Citrullination.gz',
    'dbptm_n_glycosylation_gz': 'https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/N-linked_Glycosylation.gz',
    'dbptm_o_glycosylation_gz': 'https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/O-linked_Glycosylation.gz',
    'dbptm_myristoylation_gz': 'https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Myristoylation.gz',
    # CPLM 4.0 lysine methylation (Zhang et al., Nucleic Acids Res 2022); stage with 39_stage_cplm4_lys_methylation.py.
    'cplm4_methylation_zip': 'https://cplm.biocuckoo.cn/Download/Methylation.zip',
    # CARM1 knockout ADMA peptides (Shishkova et al., Nat Commun 2017, Supplementary Data 1, CC BY 4.0); 41_stage_carm1_ko_sites.py.
    'carm1_ko_xlsx': 'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fncomms15571/MediaObjects/41467_2017_BFncomms15571_MOESM422_ESM.xlsx',
    # Human Protein Atlas subcellular locations and its UniProt mapping (CC BY 4.0); 42_stage_hpa_subcellular.py.
    'hpa_subcellular_zip': 'https://www.proteinatlas.org/download/tsv/subcellular_location.tsv.zip',
    'hpa_atlas_zip': 'https://www.proteinatlas.org/download/proteinatlas.tsv.zip',
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--outdir', default='PTM_data/external')
    ap.add_argument('--only', nargs='*', default=[])
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    names = args.only if args.only else list(URLS)
    for name in names:
        url = URLS[name]
        dest = outdir / Path(url).name
        print(f'downloading {name} -> {dest}')
        with requests.get(url, stream=True, timeout=120) as r:
            r.raise_for_status()
            with open(dest, 'wb') as fh:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        fh.write(chunk)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Stage Human Protein Atlas subcellular locations for the canonical human proteome.

The HPA subcellular section (Thul et al., Science 356:eaal3321, 2017, doi:10.1126/science.aal3321; CC BY 4.0,
proteinatlas.org/about/licence) maps proteins to 49 locations by immunofluorescence, with a reliability score per
location (Enhanced, Supported, Approved, Uncertain). Two downloads are read:

  subcellular_location.tsv(.zip)   Ensembl gene, main and additional locations, the locations at each reliability
  proteinatlas.tsv(.zip)           Ensembl gene to UniProt accession(s)

A location is kept for every reviewed canonical accession (--fasta) that HPA maps to the gene. Genes without a UniProt
accession, or whose accessions are not in the FASTA, are counted in the summary and left out.

Outputs (in --outdir, default $PTM_RESULTS_ROOT/hpa)
  hpa_subcellular_locations.tsv   one row per accession and location: main (1) or additional (0), reliability
  hpa_summary.json                input checksums and file dates, counts

With --explorer-json, the table is merged into an existing Site Explorer export (protein field hpa_locations, and a
per-location test of methylated against other profiled proteins in meta.hpa); 36_export_site_explorer_data.py applies
the same merge (explorer_merge.py) when it builds the export.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import os
import zipfile
from collections import Counter
from pathlib import Path

import pandas as pd

from common import parse_fasta, save_table
from explorer_merge import HPA_RELIABILITY, merge_hpa, update_export


def read_table(path: Path) -> tuple[list[dict], str, str]:
    """Rows of a TSV given plain or as the HPA zip; with its sha256 and, for a zip, the member's date."""
    raw = path.read_bytes()
    digest, member = hashlib.sha256(raw).hexdigest(), ''
    if zipfile.is_zipfile(io.BytesIO(raw)):
        with zipfile.ZipFile(io.BytesIO(raw)) as zf:
            names = [n for n in zf.namelist() if n.endswith('.tsv')]
            if len(names) != 1:
                raise SystemExit(f'expected one .tsv in {path}, found {names}')
            info = zf.getinfo(names[0])
            member = f'{info.filename} ({dt.datetime(*info.date_time).date().isoformat()})'
            raw = zf.read(names[0])
    return list(csv.DictReader(io.StringIO(raw.decode('utf-8')), delimiter='\t')), digest, member


def main() -> None:
    results_root = Path(os.environ.get('PTM_RESULTS_ROOT', 'PTM_results'))
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--subcellular', default=os.environ.get('PTM_HPA_SUBCELLULAR', ''),
                    help='HPA subcellular_location.tsv or its zip (PTM_HPA_SUBCELLULAR).')
    ap.add_argument('--atlas', default=os.environ.get('PTM_HPA_ATLAS', ''),
                    help='HPA proteinatlas.tsv or its zip, for the Ensembl-to-UniProt mapping (PTM_HPA_ATLAS).')
    ap.add_argument('--version', default=os.environ.get('PTM_HPA_VERSION', ''),
                    help='HPA release the files come from, as shown on proteinatlas.org/about/download (e.g. 25.1).')
    ap.add_argument('--fasta', default=os.environ.get('PTM_CANONICAL_FASTA', ''),
                    help='UniProt reviewed canonical FASTA (PTM_CANONICAL_FASTA).')
    ap.add_argument('--outdir', default=str(results_root / 'hpa'))
    ap.add_argument('--explorer-json', default='',
                    help='Optionally merge the table into this Site Explorer export in place.')
    args = ap.parse_args()
    if not args.subcellular or not args.atlas or not args.fasta:
        raise SystemExit('--subcellular, --atlas and --fasta are required (or set PTM_HPA_SUBCELLULAR, PTM_HPA_ATLAS and PTM_CANONICAL_FASTA)')

    sub_rows, sub_sha, sub_member = read_table(Path(args.subcellular))
    atlas_rows, atlas_sha, atlas_member = read_table(Path(args.atlas))
    canonical = parse_fasta(args.fasta)
    accessions = {row['Ensembl']: [a.strip() for a in row['Uniprot'].split(',') if a.strip()] for row in atlas_rows}

    out, genes = [], Counter()
    for row in sub_rows:
        accs = accessions.get(row['Gene'], [])
        kept = [a for a in accs if a in canonical]
        genes['no_uniprot' if not accs else 'not_in_canonical_fasta' if not kept else 'kept'] += 1
        main = {loc for loc in row['Main location'].split(';') if loc}
        reliability = {loc: level for level in HPA_RELIABILITY for loc in row[level].split(';') if loc}
        for acc in kept:
            for loc in sorted(main | set(row['Additional location'].split(';')) - {''}):
                out.append({'canonical_UniProtAC': acc, 'gene': row['Gene name'], 'ensembl_gene': row['Gene'], 'location': loc,
                            'main': int(loc in main), 'reliability': reliability.get(loc, row['Reliability'])})
    table = pd.DataFrame(out).drop_duplicates(['canonical_UniProtAC', 'location']).sort_values(['canonical_UniProtAC', 'location'])

    outdir = Path(args.outdir)
    save_table(table, outdir / 'hpa_subcellular_locations.tsv')
    summary = {
        'generated': dt.date.today().isoformat(),
        'version': args.version,
        'input': {'subcellular': str(Path(args.subcellular)), 'atlas': str(Path(args.atlas))},
        'input_sha256': {'subcellular': sub_sha, 'atlas': atlas_sha},
        'input_member': {'subcellular': sub_member, 'atlas': atlas_member},
        'canonical_fasta': str(Path(args.fasta)),
        'genes': dict(genes),
        'proteins': int(table['canonical_UniProtAC'].nunique()),
        'locations': int(table['location'].nunique()),
        'reliability': dict(Counter(table['reliability'])),
    }
    if args.explorer_json:
        summary['explorer'] = update_export(args.explorer_json, merge_hpa, table, summary)
    (outdir / 'hpa_summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()

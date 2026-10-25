#!/usr/bin/env python3
"""Fetch the TCGA PanCancer Atlas protein changes at arginines from the cBioPortal web API.

cBioPortal (Cerami et al. 2012, Cancer Discov 2:401, doi:10.1158/2159-8290.CD-12-0095; Gao et al. 2013, Sci Signal
6:pl1, doi:10.1126/scisignal.2004088; de Bruijn et al. 2023, Cancer Res 83:3861, doi:10.1158/0008-5472.CAN-23-0816)
serves the 32 TCGA PanCancer Atlas studies (Hoadley et al. 2018, Cell 173:291, doi:10.1016/j.cell.2018.03.022), whose
mutations are the MC3 calls (Ellrott et al. 2018, Cell Syst 6:271, doi:10.1016/j.cels.2018.03.002). Its data are
under the ODC Open Database License, with attribution to the original studies
(https://docs.cbioportal.org/user-guide/faq/).

For every study named *_tcga_pan_can_atlas_2018, the mutations of its sequenced samples are fetched for every Entrez
gene that org.Hs.eg.db names (--symbol-entrez, batches of 2,000), and the missense and nonsense changes at an arginine
(protein change R<position><residue or *>) are kept. Positions are on the transcript cBioPortal annotates
(refseqMrnaId); scripts/44_stage_depmap.py keeps a change only where the UniProt canonical sequence has an arginine at
that position.

Outputs (in --outdir)
  tcga_pancan_arginine_changes_raw.tsv   one row per sample and change, as cBioPortal gives it
  tcga_pancan_studies.tsv                each study: name, sequenced samples, citation, PubMed IDs
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import re
import time
import urllib.request
from pathlib import Path

from common import own_entrez

API = 'https://www.cbioportal.org/api'
CHANGE = re.compile(r'^R(\d+)([A-Z*])$')
KEEP_TYPES = {'Missense_Mutation', 'Nonsense_Mutation'}
FIELDS = ['study_id', 'sample_id', 'patient_id', 'entrez_id', 'hugo_symbol', 'protein_change', 'mutation_type',
          'refseq_mrna_id', 'ncbi_build', 'chr', 'start_position', 'reference_allele', 'variant_allele',
          'tumor_alt_count', 'tumor_ref_count']


def call(path: str, body: dict | None = None, tries: int = 5):
    url = API + path
    for k in range(tries):
        try:
            req = urllib.request.Request(url, data=None if body is None else json.dumps(body).encode(),
                                         headers={'Content-Type': 'application/json', 'Accept': 'application/json'})
            with urllib.request.urlopen(req, timeout=600) as r:
                return json.load(r)
        except Exception as e:  # network or server error: wait and try again
            if k == tries - 1:
                raise SystemExit(f'{url}: {e}')
            time.sleep(2 ** (k + 1))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--symbol-entrez', default='inputs/disease/orghsegdb_3.22.0_symbol_entrez.tsv')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--batch', type=int, default=2000)
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    genes = sorted({int(e) for _, e in own_entrez(args.symbol_entrez)})
    studies = sorted((s for s in call('/studies?projection=DETAILED') if s['studyId'].endswith('_tcga_pan_can_atlas_2018')), key=lambda s: s['studyId'])
    rows, meta = [], []
    for s in studies:
        sid = s['studyId']
        lists = {x['sampleListId']: x for x in call(f'/studies/{sid}/sample-lists?projection=DETAILED')}
        sl = f'{sid}_sequenced' if f'{sid}_sequenced' in lists else f'{sid}_all'
        n0 = len(rows)
        for i in range(0, len(genes), args.batch):
            for m in call(f'/molecular-profiles/{sid}_mutations/mutations/fetch?projection=DETAILED',
                          {'entrezGeneIds': genes[i:i + args.batch], 'sampleListId': sl}):
                pc = m.get('proteinChange') or ''
                if m.get('mutationType') not in KEEP_TYPES or not CHANGE.match(pc):
                    continue
                rows.append({'study_id': sid, 'sample_id': m['sampleId'], 'patient_id': m['patientId'], 'entrez_id': m['entrezGeneId'],
                             'hugo_symbol': m['gene']['hugoGeneSymbol'], 'protein_change': pc, 'mutation_type': m['mutationType'],
                             'refseq_mrna_id': m.get('refseqMrnaId') or '', 'ncbi_build': m.get('ncbiBuild') or '', 'chr': m.get('chr') or '',
                             'start_position': m.get('startPosition'), 'reference_allele': m.get('referenceAllele') or '',
                             'variant_allele': m.get('variantAllele') or '', 'tumor_alt_count': m.get('tumorAltCount'),
                             'tumor_ref_count': m.get('tumorRefCount')})
        meta.append({'study_id': sid, 'name': s['name'], 'cancer_type': (s.get('cancerType') or {}).get('name', ''),
                     'sample_list': sl, 'sequenced_samples': len(lists[sl]['sampleIds']), 'citation': s.get('citation') or '',
                     'pmid': s.get('pmid') or '', 'arginine_changes': len(rows) - n0})
        print(sid, meta[-1]['sequenced_samples'], 'samples', meta[-1]['arginine_changes'], 'arginine changes', flush=True)
    rows.sort(key=lambda r: (r['study_id'], r['hugo_symbol'], r['sample_id'], r['protein_change']))
    with open(out / 'tcga_pancan_arginine_changes_raw.tsv', 'w', newline='') as f:
        w = csv.DictWriter(f, FIELDS, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    with open(out / 'tcga_pancan_studies.tsv', 'w', newline='') as f:
        w = csv.DictWriter(f, list(meta[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(meta)
    (out / 'tcga_pancan_fetch.json').write_text(json.dumps({'fetched': dt.date.today().isoformat(), 'api': API, 'studies': len(meta),
                                                            'genes_queried': len(genes), 'arginine_changes': len(rows)}, indent=2) + '\n')


if __name__ == '__main__':
    main()

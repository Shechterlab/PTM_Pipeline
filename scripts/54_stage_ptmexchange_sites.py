#!/usr/bin/env python3
"""Stage the PTMeXchange human site tables (acetylation, methylation, SUMOylation, ubiquitination) on the reference.

    python scripts/54_stage_ptmexchange_sites.py --outdir STAGE --downloaded YYYY-MM-DD [--cache DIR]
    python scripts/apply_staged.py STAGE

PTMeXchange (ProteomeXchange; Kalyuzhnyy et al.) reprocesses public mass-spectrometry datasets with one pipeline and
grades every site Gold, Silver or Bronze by its false localisation rate. Its human "GSB" tables (FILES below) give the
protein, the position on the UniProt sequences it searched (2024-04), the residue, a decoy flag and, per dataset, the
modified peptide with the modified position. The PRIDE deposits of these builds are CC0 (PRIDE API licence field).

A site is kept where it is not a decoy, the reference sequence (--fasta, UniProt 2026_03) has the residue at the
position, and at least one of its peptides, stripped of modifications, sits in the reference sequence so that its
modified residue falls on that position ('peptide'); a site whose residue matches but none of whose peptides place it
is counted and left out. Each row: uniprot_accession, position, residue, ptm (acetyl_K, methyl_R, methyl_K, sumo_K,
ubiquitin_K), confidence (Gold, Silver, Bronze), psm_count (5% FLR, where the table gives it), datasets (PXD
accessions with a peptide), file.

Writes layers/ptmexchange_human_sites.tsv, ptmexchange_summary.json, manifest_rows.tsv and raw_files_rows.tsv.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import re
from collections import Counter
from pathlib import Path
from urllib.request import Request, urlopen

import pandas as pd

from common import md5_file, parse_fasta

BASE = 'https://www.proteomexchange.org/ptmexchange/Files/'
FILES = {  # file: (ptm by residue, PRIDE deposit)
    'Human_GSB_Acetyl_all_prots_1106.csv.gz': ({'K': 'acetyl_K'}, 'PXD079218'),
    'Human_GSB_Methyl_all_prots_0406.csv.gz': ({'R': 'methyl_R', 'K': 'methyl_K'}, 'PXD076372'),
    'Human_GSB_SUMO_all_prots_0924.csv.gz': ({'K': 'sumo_K'}, 'PXD068489'),
    'Human_GSB_Ubiq_all_prots_1106.csv.gz': ({'K': 'ubiquitin_K'}, 'PXD068989'),
}
CITATION = ('PTMeXchange human builds, PRIDE PXD079218 (acetylation), PXD076372 (methylation), PXD068489 (SUMOylation), '
            'PXD068989 (ubiquitination); ProteomeXchange, https://www.proteomexchange.org/ptmexchange/')
LICENCE = 'CC0 (PRIDE deposits PXD079218, PXD076372, PXD068489, PXD068989: PRIDE API licence field)'
PEPTIDE = re.compile(r'^(.*)-(\d+)$')
MOD = re.compile(r'\[[^\]]*\]')


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def placed(seq: str, pos: int, peptides: list[str]) -> bool:
    """Does any peptide (modifications stripped, '-k' its modified position) put its residue k on pos?"""
    for p in peptides:
        m = PEPTIDE.match(p)
        if not m:
            continue
        plain, k = MOD.sub('', m.group(1)), int(m.group(2))
        start = seq.find(plain)
        while start >= 0:
            if start + k == pos:
                return True
            start = seq.find(plain, start + 1)
    return False


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--fasta', default='inputs/reference/uniprot_human_reviewed_2026_03.fasta')
    ap.add_argument('--cache', type=Path, help='read the files from here instead of downloading them')
    ap.add_argument('--manifest', default='inputs/MANIFEST.tsv', help='read only: its columns')
    ap.add_argument('--raw-files', default='inputs/layers/raw_files.tsv', help='read only: its columns')
    ap.add_argument('--downloaded', required=True)
    ap.add_argument('--outdir', required=True, type=Path)
    a = ap.parse_args()
    out = a.outdir
    out.mkdir(parents=True, exist_ok=True)
    seqs = parse_fasta(a.fasta)
    rows, left, raw = [], Counter(), []
    for name, (ptm_of, pxd) in FILES.items():
        path = (a.cache / name) if a.cache else out / name
        if not a.cache:
            with urlopen(Request(BASE + name, headers={'User-Agent': 'PTM-pipeline-stager/1.0'}), timeout=600) as r:
                path.write_bytes(r.read())
        raw.append((name, path))
        tab = pd.read_csv(io.BytesIO(gzip.decompress(path.read_bytes())), dtype=str, keep_default_na=False)
        pos_col = 'Protein_pos' if 'Protein_pos' in tab.columns else 'Protein position'
        pep_cols = [c for c in tab.columns if c.lower().endswith('_peptide_mod_pos')]
        psm_col = next((c for c in tab.columns if c.startswith('Sum_of_PSM_counts')), None)
        for r in tab.to_dict('records'):
            if r.get('Decoy_mod', '0') not in ('0', ''):
                left[f'{name}: decoy'] += 1
                continue
            parts = r['Protein'].split('|')
            acc = parts[1] if len(parts) > 1 else parts[0]
            res = r['PTM_residue']
            if res not in ptm_of:
                left[f'{name}: residue {res}'] += 1
                continue
            try:
                pos = int(float(r[pos_col]))
            except ValueError:
                left[f'{name}: no position'] += 1
                continue
            seq = seqs.get(acc)
            if seq is None:
                left[f'{name}: not in reference'] += 1
                continue
            if not (1 <= pos <= len(seq)) or seq[pos - 1] != res:
                left[f'{name}: residue differs'] += 1
                continue
            peps = {c.split('_')[0]: r[c] for c in pep_cols if r[c]}
            if not placed(seq, pos, list(peps.values())):
                left[f'{name}: no peptide places it'] += 1
                continue
            rows.append({'uniprot_accession': acc, 'position': pos, 'residue': res, 'ptm': ptm_of[res],
                         'confidence': r['PTM_FLR_category'], 'psm_count': r.get(psm_col, '') if psm_col else '',
                         'datasets': ';'.join(sorted(peps)), 'file': name})
    df = pd.DataFrame(rows).drop_duplicates(['uniprot_accession', 'position', 'ptm'])
    df = df.sort_values(['uniprot_accession', 'position', 'ptm']).reset_index(drop=True)
    df.to_csv(out / 'ptmexchange_human_sites.tsv', sep='\t', index=False)
    summary = {
        'script': 'scripts/54_stage_ptmexchange_sites.py', 'downloaded': a.downloaded,
        'files': {n: {'md5': md5_file(p), 'pride': FILES[n][1]} for n, p in raw},
        'reference': {'fasta': a.fasta, 'md5': md5_file(a.fasta)},
        'rows': len(df), 'by_ptm': {k: int(v) for k, v in df.ptm.value_counts().items()},
        'by_ptm_confidence': {f'{k[0]} {k[1]}': int(v) for k, v in df.groupby(['ptm', 'confidence']).size().items()},
        'left_out': dict(left),
    }
    (out / 'ptmexchange_summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    mcols = open(a.manifest, encoding='utf-8').readline().rstrip('\n').split('\t')
    rcols = open(a.raw_files, encoding='utf-8').readline().rstrip('\n').split('\t')
    by = summary['by_ptm']
    source = 'Derived by scripts/54_stage_ptmexchange_sites.py from ' + '; '.join(BASE + n for n in FILES)
    version = f'PTMeXchange human builds (file dates in the names), downloaded {a.downloaded}'
    mrows = [
        {'file': 'layers/ptmexchange_human_sites.tsv', 'rows': len(df),
         'content': (f'PTMeXchange reprocessed human sites placed on the 2026_03 reference by residue and peptide: '
                     f'{by.get("acetyl_K", 0):,} acetyl-K, {by.get("ubiquitin_K", 0):,} ubiquitin-K, {by.get("sumo_K", 0):,} SUMO-K, '
                     f'{by.get("methyl_R", 0):,} methyl-R and {by.get("methyl_K", 0):,} methyl-K, each Gold, Silver or Bronze'),
         'derived_from': '; '.join(f'{n} (MD5 {summary["files"][n]["md5"]})' for n in FILES),
         'how_made_or_checked': 'scripts/54_stage_ptmexchange_sites.py: not a decoy; residue at the position of the reference; '
                                'a peptide of the site places it there'},
        {'file': 'layers/ptmexchange_summary.json', 'rows': '', 'derived_from': '',
         'content': 'Files, MD5s, counts by PTM and confidence, and what was left out of the PTMeXchange staging',
         'how_made_or_checked': 'scripts/54_stage_ptmexchange_sites.py'},
    ]
    for r in mrows:
        p = out / r['file'].split('/')[-1]
        r.update({'bytes': p.stat().st_size, 'md5': md5_file(p), 'redistribution': 'git', 'source': source,
                  'version_or_date': version, 'citation': CITATION, 'licence': LICENCE})
    pd.DataFrame([{c: r.get(c, '') for c in mcols} for r in mrows]).to_csv(out / 'manifest_rows.tsv', sep='\t', index=False)
    rrows = [{'file': n, 'url': BASE + n, 'bytes': p.stat().st_size, 'md5': md5_file(p), 'sha256': sha256_file(p),
              'downloaded': a.downloaded, 'version': f'PTMeXchange human build, PRIDE {FILES[n][1]}', 'kept': 'zenodo'}
             for n, p in raw]
    pd.DataFrame([{c: r.get(c, '') for c in rcols} for r in rrows]).to_csv(out / 'raw_files_rows.tsv', sep='\t', index=False)
    print(json.dumps({k: summary[k] for k in ('rows', 'by_ptm', 'left_out')}, indent=1))


if __name__ == '__main__':
    main()

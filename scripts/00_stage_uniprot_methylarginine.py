#!/usr/bin/env python3
"""Stage UniProt's methylarginine annotations (MOD_RES) of the reviewed human entries (inputs/annotation/).

    python scripts/00_stage_uniprot_methylarginine.py --fasta inputs/reference/uniprot_human_reviewed_2026_03.fasta \
        --release 2026_03 --raw-out uniprot_human_reviewed_mod_res_2026_03.tsv \
        --out inputs/annotation/uniprot_2026_03_methylarginine_mod_res.tsv

Downloads (or, with --raw, reads) every reviewed human entry with a modified residue, with its sequence, from the
UniProt REST API (CC BY 4.0), and refuses to go on unless the API serves --release. Keeps the arginine methylation
terms, each on its own row, at positions that are an arginine of an entry whose sequence equals the --fasta one:

  Omega-N-methylarginine         Rme1
  Asymmetric dimethylarginine    Rme2a
  Symmetric dimethylarginine     Rme2s
  Dimethylated arginine          Rme2
  Omega-N-methylated arginine    Rme (degree not stated)

with the evidence codes and their class: experimental (ECO:0000269 experimental evidence, ECO:0007744 combinatorial
computational and experimental, i.e. large-scale proteomics), by_similarity (ECO:0000250) or other. Features
annotated on an isoform's sequence (40 in 2026_03) are counted in the report, not kept.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from urllib.request import Request, urlopen

import pandas as pd

URL = ('https://rest.uniprot.org/uniprotkb/stream?query=(organism_id:9606)%20AND%20(reviewed:true)%20AND%20(ft_mod_res:*)'
       '&fields=accession,ft_mod_res,sequence&format=tsv')
TERMS = {'Omega-N-methylarginine': 'Rme1', 'Asymmetric dimethylarginine': 'Rme2a', 'Symmetric dimethylarginine': 'Rme2s',
         'Dimethylated arginine': 'Rme2', 'Omega-N-methylated arginine': 'Rme'}
EXPERIMENTAL = {'ECO:0000269', 'ECO:0007744'}
FEATURE = re.compile(r'MOD_RES (\d+); /note="([^"]*)"(?:; /evidence="([^"]*)")?')
# A feature on an isoform's sequence (MOD_RES P12345-2:249) is not on the canonical one: counted, not kept.
ISOFORM_FEATURE = re.compile(r'MOD_RES [A-Z0-9]+-\d+:\d+; /note="([^";]*)')


def read_fasta(path: str) -> dict[str, str]:
    seqs, acc = {}, None
    for line in open(path):
        if line.startswith('>'):
            acc = line.split('|')[1]
            seqs[acc] = []
        elif acc:
            seqs[acc].append(line.strip())
    return {a: ''.join(s) for a, s in seqs.items()}


def evidence_class(codes: list[str]) -> str:
    if any(c in EXPERIMENTAL for c in codes):
        return 'experimental'
    if 'ECO:0000250' in codes:
        return 'by_similarity'
    return 'other'


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--release', required=True, help='the UniProt release the reference is, e.g. 2026_03')
    ap.add_argument('--raw', help='read this download instead of fetching')
    ap.add_argument('--raw-out', help='keep the download here')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    if args.raw:
        raw_path = args.raw
    else:
        with urlopen(Request(URL, headers={'User-Agent': 'PTM-pipeline-stager/1.0'}), timeout=600) as response:
            release = response.headers.get('X-UniProt-Release', '')
            if release != args.release:
                sys.exit(f'the UniProt API serves release {release!r}, not {args.release}')
            data = response.read()
        raw_path = args.raw_out or 'uniprot_mod_res.tsv'
        Path(raw_path).write_bytes(data)
        print(f'UniProt release {release}: {raw_path}')
    raw = pd.read_csv(raw_path, sep='\t', dtype=str, keep_default_na=False)
    reference = read_fasta(args.fasta)

    rows, skipped = [], {'not_in_reference': 0, 'sequence_differs': 0, 'not_arginine': 0, 'on_an_isoform': 0}
    for r in raw.itertuples(index=False):
        acc, features, sequence = r[0], r[1], r[2]
        skipped['on_an_isoform'] += sum(note.strip() in TERMS for note in ISOFORM_FEATURE.findall(features))
        terms = [(int(pos), note.split(';')[0].strip(), 'alternate' in note, [e.split('|')[0] for e in ev.split(', ') if e])
                 for pos, note, ev in FEATURE.findall(features)]
        terms = [t for t in terms if t[1] in TERMS]
        if not terms:
            continue
        if acc not in reference:
            skipped['not_in_reference'] += len(terms)
            continue
        if sequence != reference[acc]:
            skipped['sequence_differs'] += len(terms)
            continue
        for pos, term, alternate, codes in terms:
            if sequence[pos - 1] != 'R':
                skipped['not_arginine'] += 1
                continue
            rows.append({'uniprot_accession': acc, 'position': pos, 'term': term, 'state': TERMS[term], 'alternate': alternate,
                         'evidence': ';'.join(sorted(set(codes))), 'evidence_class': evidence_class(codes)})
    out = pd.DataFrame(rows).sort_values(['uniprot_accession', 'position', 'state']).reset_index(drop=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep='\t', index=False)
    sites = out.drop_duplicates(['uniprot_accession', 'position'])
    print(f'{len(out):,} methylarginine annotations at {len(sites):,} arginines of {out["uniprot_accession"].nunique():,} entries; '
          f'by evidence: {out["evidence_class"].value_counts().to_dict()}; left out: {skipped}')


if __name__ == '__main__':
    main()

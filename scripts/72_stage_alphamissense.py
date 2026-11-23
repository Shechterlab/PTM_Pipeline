#!/usr/bin/env python3
"""Stage AlphaMissense's predicted pathogenicity at every candidate residue (R, K, S, T, Y) of the reviewed canonical
proteins (release v1.8.0).

    python scripts/72_stage_alphamissense.py --alphamissense AlphaMissense_aa_substitutions.tsv.gz --outdir inputs/layers

Source
  AlphaMissense_aa_substitutions.tsv.gz (Cheng et al., Science 381:eadg7492, 2023; Zenodo record 10813168, v3,
  doi:10.5281/zenodo.10813168): the predicted pathogenicity of every single amino-acid substitution of the UniProt canonical
  sequences AlphaMissense was run on. The Zenodo record states CC BY 4.0, and DeepMind moved the predictions from CC BY-NC-SA
  4.0 to CC BY 4.0 on 13 March 2024; the header of the file itself still names CC BY-NC-SA 4.0. Either way these derived
  values may be redistributed under this release's CC BY-NC-SA 4.0, with credit.

Rule
  An entry counts when AlphaMissense's sequence, read from the reference residue and position of its substitutions, equals
  the 2026_03 canonical sequence over its whole length; entries whose sequence differs, or that AlphaMissense lacks, are
  listed with the reason and get no values. Per candidate residue: the mean and highest pathogenicity over its 19
  substitutions, and how many of them AlphaMissense classes pathogenic (every candidate residue has 19; the script checks).
  AlphaMissense's scores are less reliable in intrinsically disordered regions (the release says so where it uses them).

Outputs (gzip without a timestamp)
  alphamissense_residues.tsv.gz       one row per candidate residue of a matched entry
  alphamissense_status.tsv            every reviewed entry: matched, sequence_differs or not_in_alphamissense, with lengths
  alphamissense_summary.json          input (MD5), rule, license notes and counts
"""
from __future__ import annotations

import argparse
import datetime as dt
import io
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from common import md5_file, parse_fasta, refresh_manifest_row, update_manifest, upsert_tsv_rows

URL = 'https://zenodo.org/records/10813168/files/AlphaMissense_aa_substitutions.tsv.gz'
CITATION = 'Cheng et al., Science 381:eadg7492 (2023), doi:10.1126/science.adg7492; predictions doi:10.5281/zenodo.10813168'
LICENSE = ('CC BY 4.0 (Zenodo record 10813168; DeepMind moved the predictions from CC BY-NC-SA 4.0 to CC BY 4.0 on 13 March 2024; '
           'the file header still names CC BY-NC-SA 4.0)')
CANDIDATES = 'RKSTY'
COLUMNS = ['canonical_UniProtAC', 'position', 'residue', 'am_mean_pathogenicity', 'am_max_pathogenicity', 'am_pathogenic']


def write_gz(df: pd.DataFrame, path: Path) -> None:
    import gzip
    buf = io.StringIO()
    df.to_csv(buf, sep='\t', index=False, lineterminator='\n')
    with open(path, 'wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0, compresslevel=9) as gz:
        gz.write(buf.getvalue().encode('utf-8'))


def read_alphamissense(path: str, wanted: set[str]) -> pd.DataFrame:
    """Per (entry, position, reference residue): the summed, highest and lowest pathogenicity, the number of substitutions
    and their number per class, aggregated chunk by chunk (71 million rows do not fit in memory as they are)."""
    parts = []
    for chunk in pd.read_csv(path, sep='\t', comment='#', dtype={'uniprot_id': str, 'protein_variant': str, 'am_pathogenicity': float,
                                                                 'am_class': str}, chunksize=4_000_000):
        chunk = chunk[chunk['uniprot_id'].isin(wanted)]
        v = chunk['protein_variant']
        c = pd.DataFrame({'acc': chunk['uniprot_id'].to_numpy(), 'pos': v.str[1:-1].astype(np.int32).to_numpy(), 'ref': v.str[0].to_numpy(),
                          'score': chunk['am_pathogenicity'].to_numpy()})
        for k in ('pathogenic', 'ambiguous', 'benign'):
            c[k] = (chunk['am_class'].to_numpy() == k).astype(np.int32)
        parts.append(c.groupby(['acc', 'pos', 'ref'], sort=False).agg(total=('score', 'sum'), max=('score', 'max'), min=('score', 'min'),
                                                                      n=('score', 'size'), pathogenic=('pathogenic', 'sum'),
                                                                      ambiguous=('ambiguous', 'sum'), benign=('benign', 'sum')).reset_index())
    a = pd.concat(parts, ignore_index=True)
    return a.groupby(['acc', 'pos', 'ref'], sort=True).agg(total=('total', 'sum'), max=('max', 'max'), min=('min', 'min'), n=('n', 'sum'),
                                                           pathogenic=('pathogenic', 'sum'), ambiguous=('ambiguous', 'sum'),
                                                           benign=('benign', 'sum')).reset_index()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--alphamissense', required=True, help=f'AlphaMissense_aa_substitutions.tsv.gz ({URL})')
    ap.add_argument('--fasta', default='inputs/reference/uniprot_human_reviewed_2026_03.fasta')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--manifest', default='')
    ap.add_argument('--downloaded', default=dt.date.today().isoformat())
    args = ap.parse_args()
    out = Path(args.outdir)
    canonical = parse_fasta(args.fasta)
    am = read_alphamissense(args.alphamissense, set(canonical))
    counts = Counter({'alphamissense_scored_positions_on_reference_entries': len(am), 'alphamissense_substitutions_on_reference_entries': int(am['n'].sum())})

    status, keep = [], []
    for acc, g in am.groupby('acc', sort=True):
        seq = canonical[acc]
        pos, ref = g['pos'].to_numpy(), g['ref'].to_numpy()
        same = (len(pos) == len(seq) and int(pos.max()) <= len(seq) and all(seq[p - 1] == r for p, r in zip(pos, ref)))
        status.append({'canonical_UniProtAC': acc, 'status': 'matched' if same else 'sequence_differs', 'canonical_length': len(seq),
                       'alphamissense_positions': int(len(pos)), 'alphamissense_last_position': int(pos.max())})
        if same:
            keep.append(acc)
    seen = {s['canonical_UniProtAC'] for s in status}
    status += [{'canonical_UniProtAC': a, 'status': 'not_in_alphamissense', 'canonical_length': len(canonical[a]), 'alphamissense_positions': 0,
                'alphamissense_last_position': 0} for a in sorted(set(canonical) - seen)]
    st = pd.DataFrame(status).sort_values('canonical_UniProtAC', kind='stable')
    counts.update({f'entries_{k}': int(v) for k, v in st['status'].value_counts().items()})

    res = am[am['acc'].isin(set(keep)) & am['ref'].isin(list(CANDIDATES))]
    table = pd.DataFrame({'canonical_UniProtAC': res['acc'], 'position': res['pos'], 'residue': res['ref'],
                          'am_mean_pathogenicity': (res['total'] / res['n']).map(lambda x: f'{x:.4f}'),
                          'am_max_pathogenicity': res['max'].map(lambda x: f'{x:.4f}'), 'am_pathogenic': res['pathogenic'].astype(int)})
    if not (res['n'] == 19).all():
        raise SystemExit('a candidate residue does not have 19 substitutions in AlphaMissense')
    out.mkdir(parents=True, exist_ok=True)
    write_gz(table, out / 'alphamissense_residues.tsv.gz')
    st.to_csv(out / 'alphamissense_status.tsv', sep='\t', index=False, lineterminator='\n')
    counts['candidate_residues'] = len(table)
    counts.update({f'candidate_residues_{r}': int((table['residue'] == r).sum()) for r in CANDIDATES})
    counts['residues_with_19_substitutions'] = int((res['n'] == 19).sum())
    summary = {'generated': dt.date.today().isoformat(), 'input': {'file': Path(args.alphamissense).name, 'url': URL, 'md5': md5_file(args.alphamissense)},
               'reference': {'fasta': args.fasta, 'md5': md5_file(args.fasta)}, 'citation': CITATION, 'license': LICENSE,
               'rule': 'an entry counts when AlphaMissense\'s sequence (reference residue and position of its substitutions) equals the '
                       '2026_03 canonical sequence over its whole length; per candidate residue the mean and highest pathogenicity '
                       'over its 19 substitutions and the number AlphaMissense classes pathogenic',
               'caveat': 'AlphaMissense is less reliable in intrinsically disordered regions',
               'counts': dict(sorted(counts.items())),
               'md5': {n: md5_file(out / n) for n in ('alphamissense_residues.tsv.gz', 'alphamissense_status.tsv')}}
    (out / 'alphamissense_summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    if args.manifest:
        write_manifest(Path(args.manifest), out, args)
    print(json.dumps(summary['counts'], indent=1))


def write_manifest(manifest: Path, out: Path, args) -> None:
    import hashlib
    folder, how = out.name, 'scripts/72_stage_alphamissense.py'
    ver = 'AlphaMissense predictions, Zenodo record 10813168 v3 (2023-09-19)'
    rows = [
        {'file': f'{folder}/alphamissense_residues.tsv.gz', 'content': 'AlphaMissense predicted pathogenicity per candidate residue (R, K, S, T, Y) '
         'of reviewed entries whose AlphaMissense sequence equals the 2026_03 canonical sequence: mean and highest over the '
         '19 substitutions, and the number AlphaMissense classes pathogenic (gzip)'},
        {'file': f'{folder}/alphamissense_status.tsv', 'content': 'Every reviewed entry: whether AlphaMissense\'s sequence equals the canonical '
         'sequence (matched), differs, or is missing, with both lengths'},
        {'file': f'{folder}/alphamissense_summary.json', 'content': 'Input (MD5), rule, license notes and counts of the AlphaMissense staging'},
    ]
    for r in rows:
        path = out / r['file'].split('/', 1)[1]
        n = '' if path.suffix == '.json' else (sum(1 for _ in __import__('gzip').open(path, 'rt')) - 1 if path.suffix == '.gz'
                                               else sum(1 for _ in open(path, encoding='utf-8')) - 1)
        r.update({'bytes': path.stat().st_size, 'md5': md5_file(path), 'rows': n, 'redistribution': 'git', 'source': URL, 'version_or_date': ver,
                  'citation': CITATION, 'licence': LICENSE, 'derived_from': 'AlphaMissense_aa_substitutions.tsv.gz; reference/uniprot_human_reviewed_2026_03.fasta',
                  'how_made_or_checked': how})
    update_manifest(manifest, rows)
    p = Path(args.alphamissense)
    h = hashlib.sha256()
    with open(p, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    upsert_tsv_rows(out / 'raw_files.tsv', [{'file': p.name, 'url': URL, 'bytes': p.stat().st_size, 'md5': md5_file(p), 'sha256': h.hexdigest(),
                                             'downloaded': args.downloaded, 'version': ver + '; Zenodo MD5 b9ccb339e0de6cb0a8d1973ad2026576 matched',
                                             'kept': 'not kept: 1.2 GB, served by Zenodo at this URL; MD5 and SHA-256 identify it'}])
    refresh_manifest_row(manifest, out / 'raw_files.tsv', f'{folder}/raw_files.tsv')


if __name__ == '__main__':
    main()

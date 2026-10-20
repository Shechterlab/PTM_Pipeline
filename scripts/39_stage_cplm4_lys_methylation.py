#!/usr/bin/env python3
"""Stage and validate CPLM 4.0 lysine methylation sites against the canonical human proteome.

CPLM 4.0 (Zhang et al., Nucleic Acids Res 50:D451-D459, 2022, doi:10.1093/nar/gkab849; cplm.biocuckoo.cn)
compiles lysine modification sites. Its Methylation download (Download/Methylation.zip, one headerless
tab-separated Methylation.txt) has nine columns:

  CPLM id, UniProt accession, position, modification, gene, organism, protein sequence, evidence, references

Rows are kept when they are for the chosen organism (Homo sapiens), describe methylation, sit on a lysine of
CPLM's own sequence and on a lysine of the same accession and position in the canonical FASTA, and either the
two sequences are identical or the +/-7 residue windows around the lysine agree. Nothing is remapped: a site
whose sequence moved is reported in the rejected table with its reason. The evidence label ("Exp.", "Dat.",
"Exp./Dat.") and the references field are carried through verbatim; some "Dat." rows hold identifiers that
are not PubMed IDs.

Outputs (in --outdir, default $PTM_RESULTS_ROOT/cplm4)
  cplm4_human_lys_methyl_validated.tsv   one row per canonical accession and position
  cplm4_human_lys_methyl_rejected.tsv    input rows that failed a check, with the reason
  cplm4_summary.json                     input checksum, counts per check, sites and proteins kept

With --explorer-json, the validated sites are also merged into an existing Site Explorer export
(docs/data/methylarg_sites.json) as Lys methylation entries with their CPLM source bits, replacing any CPLM
entries already there; 36_export_site_explorer_data.py applies the same merge (explorer_merge.py) when it builds the export.
"""
from __future__ import annotations

import argparse
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
from explorer_merge import merge_cplm4, update_export

CPLM4_COLUMNS = ['cplm_id', 'source_UniProtAC', 'position', 'modification', 'gene', 'organism', 'cplm_sequence',
                 'evidence', 'references']
WINDOW_FLANK = 7


def read_cplm(path: Path) -> tuple[pd.DataFrame, str, str]:
    """Read a CPLM download given as the zip, the extracted text file, or a directory holding either."""
    if path.is_dir():
        for name in ('Methylation.txt', 'Methylation.zip'):
            if (path / name).exists():
                path = path / name
                break
        else:
            raise SystemExit(f'no Methylation.txt or Methylation.zip in {path}')
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    member = ''
    if zipfile.is_zipfile(io.BytesIO(raw)):
        with zipfile.ZipFile(io.BytesIO(raw)) as zf:
            names = [n for n in zf.namelist() if n.lower().endswith('.txt')]
            if len(names) != 1:
                raise SystemExit(f'expected one .txt file in {path}, found {names}')
            info = zf.getinfo(names[0])
            member = f'{info.filename} ({dt.datetime(*info.date_time).isoformat()})'
            raw = zf.read(names[0])
    df = pd.read_csv(io.BytesIO(raw), sep='\t', header=None, dtype=str, keep_default_na=False, quoting=3)
    if df.shape[1] != len(CPLM4_COLUMNS):
        raise SystemExit(f'{path}: expected {len(CPLM4_COLUMNS)} columns, found {df.shape[1]}')
    df.columns = CPLM4_COLUMNS
    return df, digest, member


def window(seq: str, pos: int) -> str:
    return seq[max(0, pos - 1 - WINDOW_FLANK):pos + WINDOW_FLANK]


def check_row(row, canonical: dict[str, str]) -> tuple[str, str]:
    """Return (status, detail): status is 'identical_sequence' or 'window_match' when kept, else a reason."""
    pos = row.position_int
    own = row.cplm_sequence.strip().upper()
    if pos is None or pos < 1:
        return 'bad_position', ''
    if not (pos <= len(own) and own[pos - 1] == 'K'):
        return 'cplm_sequence_not_K', own[pos - 1] if pos <= len(own) else 'past end'
    seq = canonical.get(row.source_UniProtAC)
    if seq is None:
        return 'accession_not_in_canonical_fasta', ''
    if not (pos <= len(seq) and seq[pos - 1] == 'K'):
        return 'canonical_not_K', seq[pos - 1] if pos <= len(seq) else 'past end'
    if own == seq:
        return 'identical_sequence', ''
    if window(own, pos) == window(seq, pos):
        return 'window_match', ''
    return 'window_mismatch', f'{window(own, pos)} vs {window(seq, pos)}'


def main() -> None:
    results_root = Path(os.environ.get('PTM_RESULTS_ROOT', 'PTM_results'))
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--cplm', default=os.environ.get('PTM_CPLM4_METHYL', ''),
                    help='CPLM 4.0 Methylation.zip, Methylation.txt, or a directory holding one (PTM_CPLM4_METHYL).')
    ap.add_argument('--fasta', default=os.environ.get('PTM_CANONICAL_FASTA', ''),
                    help='UniProt reviewed canonical FASTA (PTM_CANONICAL_FASTA).')
    ap.add_argument('--organism', default='Homo sapiens')
    ap.add_argument('--outdir', default=str(results_root / 'cplm4'))
    ap.add_argument('--explorer-json', default='',
                    help='Optionally merge the validated sites into this Site Explorer export in place.')
    args = ap.parse_args()
    if not args.cplm or not args.fasta:
        raise SystemExit('--cplm and --fasta are required (or set PTM_CPLM4_METHYL and PTM_CANONICAL_FASTA)')

    df, digest, member = read_cplm(Path(args.cplm))
    canonical = parse_fasta(args.fasta)
    n_input = len(df)
    df = df[df['organism'].eq(args.organism) & df['modification'].str.contains('Methylation', case=False)].copy()
    n_scope = len(df)
    df['position_int'] = [int(p) if p.strip().isdigit() else None for p in df['position']]
    checks = [check_row(row, canonical) for row in df.itertuples(index=False)]
    df['validation'] = [c[0] for c in checks]
    df['validation_detail'] = [c[1] for c in checks]
    kept = df['validation'].isin(['identical_sequence', 'window_match'])

    rejected = df[~kept].drop(columns=['cplm_sequence', 'position_int'])
    good = df[kept].copy()
    good['canonical_UniProtAC'] = good['source_UniProtAC']
    good['position'] = good['position_int'].astype(int)
    good['residue'] = 'K'
    good['site'] = 'K' + good['position'].astype(str)
    good = (
        good.groupby(['canonical_UniProtAC', 'position'], as_index=False, sort=True)
        .agg(residue=('residue', 'first'), site=('site', 'first'), gene=('gene', 'first'),
             cplm_id=('cplm_id', lambda v: ';'.join(sorted(set(v)))),
             modification=('modification', 'first'),
             evidence=('evidence', lambda v: '/'.join(e for e in ('Exp.', 'Dat.') if any(e in x.split('/') for x in v))),
             references=('references', lambda v: ';'.join(dict.fromkeys(t for x in v for t in x.split(';') if t))),
             validation=('validation', 'first'))
    )

    outdir = Path(args.outdir)
    save_table(good, outdir / 'cplm4_human_lys_methyl_validated.tsv')
    save_table(rejected, outdir / 'cplm4_human_lys_methyl_rejected.tsv')
    summary = {
        'generated': dt.date.today().isoformat(),
        'input': str(Path(args.cplm)),
        'input_sha256': digest,
        'input_member': member,
        'canonical_fasta': str(Path(args.fasta)),
        'canonical_sequences': len(canonical),
        'organism': args.organism,
        'input_rows': n_input,
        'rows_in_scope': n_scope,
        'checks': dict(Counter(df['validation'])),
        'validated_sites': len(good),
        'validated_proteins': int(good['canonical_UniProtAC'].nunique()),
        'validated_evidence': dict(Counter(good['evidence'])),
    }
    if args.explorer_json:
        summary.update(update_export(args.explorer_json, merge_cplm4, good, summary))
    (outdir / 'cplm4_summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()

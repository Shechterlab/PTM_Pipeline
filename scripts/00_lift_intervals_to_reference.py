#!/usr/bin/env python3
"""Move per-protein interval tables (InterPro, MobiDB) from the sequences they were computed on to a newer reference.

InterPro 108.0 and the MobiDB export of 2026-03-30 give coordinates on the UniProt sequences of early 2026 (the
2026-04-07 staging). The reference of release v1.1.0 is UniProt 2026_03, where 13 reviewed entries have a different
sequence. For each interval table this script:

  - keeps every entry whose sequence is the same in both;
  - shifts an entry whose new sequence is the old one with residues added at the N-terminus;
  - keeps the coordinates of an entry whose new sequence has the same length (single-residue changes) or extends the
    old one at the C-terminus;
  - drops the intervals of every other changed entry, and of entries the new reference no longer has;
  - updates `protein_length` where the table has it.

    python scripts/00_lift_intervals_to_reference.py --old-fasta OLD.fasta --new-fasta inputs/reference/....fasta \
        --tables inputs/annotation/interpro_108.0_domain_intervals.tsv ... --report lift_report.tsv

Tables are rewritten in place; the report lists what happened to each changed entry.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def read_fasta(path: Path) -> dict[str, str]:
    seqs, acc, parts = {}, None, []
    for line in open(path, encoding='utf-8'):
        if line.startswith('>'):
            if acc:
                seqs[acc] = ''.join(parts)
            acc, parts = line.split('|')[1], []
        else:
            parts.append(line.strip())
    if acc:
        seqs[acc] = ''.join(parts)
    return seqs


def mapping(old: str, new: str) -> tuple[str, int]:
    """How coordinates on old map onto new: (kind, shift)."""
    if old == new:
        return 'same', 0
    if new.endswith(old):
        return 'n_terminal_extension', len(new) - len(old)
    if len(new) == len(old):
        return 'substitutions', 0
    if new.startswith(old):
        return 'c_terminal_extension', 0
    return 'other_change', 0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--old-fasta', required=True, help='the sequences the tables were computed on')
    ap.add_argument('--new-fasta', required=True, help='the reference to move them to')
    ap.add_argument('--tables', nargs='+', required=True)
    ap.add_argument('--report', required=True)
    args = ap.parse_args()
    old, new = read_fasta(Path(args.old_fasta)), read_fasta(Path(args.new_fasta))
    changed = {a: mapping(old[a], new[a]) for a in old.keys() & new.keys() if old[a] != new[a]}
    gone = sorted(old.keys() - new.keys())

    report = []
    for table in args.tables:
        path = Path(table)
        df = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, low_memory=False)
        acc = df['canonical_UniProtAC']
        keep = ~acc.isin(gone)
        for a, (kind, shift) in sorted(changed.items()):
            rows = acc == a
            if not rows.any():
                continue
            ends = df.loc[rows, 'fragment_end'].astype(int)
            if ends.max() > len(old[a]):
                raise SystemExit(f'{path.name}: {a} has intervals beyond its old length; it was not computed on --old-fasta')
            if kind == 'other_change':
                keep &= ~rows
            elif shift:
                for col in ('fragment_start', 'fragment_end'):
                    df.loc[rows, col] = (df.loc[rows, col].astype(int) + shift).astype(str)
            report.append({'table': path.name, 'accession': a, 'old_length': len(old[a]), 'new_length': len(new[a]),
                           'change': kind, 'shift': shift, 'intervals': int(rows.sum()),
                           'action': 'dropped' if kind == 'other_change' else f'shifted by {shift}' if shift else 'kept'})
        for a in gone:
            n = int((acc == a).sum())
            if n:
                report.append({'table': path.name, 'accession': a, 'old_length': len(old[a]), 'new_length': 0,
                               'change': 'not_in_new_reference', 'shift': 0, 'intervals': n, 'action': 'dropped'})
        df = df[keep]
        if 'protein_length' in df:
            lengths = df['canonical_UniProtAC'].map(lambda a: str(len(new[a])) if a in new else '')
            df = df.assign(protein_length=lengths)
        df.to_csv(path, sep='\t', index=False, lineterminator='\n')
    pd.DataFrame(report).to_csv(args.report, sep='\t', index=False, lineterminator='\n')
    print(pd.DataFrame(report).groupby(['table', 'action']).size().to_string())


if __name__ == '__main__':
    main()

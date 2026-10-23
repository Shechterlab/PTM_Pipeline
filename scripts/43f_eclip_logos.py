#!/usr/bin/env python3
"""Sequence logos of the RNA crosslinked in cells, from the eCLIP 5-mer enrichment tables scripts/43e staged.

For each eCLIP experiment (ENCODE reproducible peaks, ENCORE reproducible enriched windows), the 5-mers enriched in the
top peaks or windows are aligned into one logo:

  1. 5-mers enriched against the background (--background: 'other', the top peaks or windows of the other RBPs of the
     same source reweighted to the experiment's mix of transcript regions, or 'flank', the same length on each side of
     each peak) with R >= --min-r, z >= --min-z and at least --min-count occurrences in the foreground, ranked by R;
     at most --top of them.
  2. The first is placed at offsets 0-4. Each next one is placed at the offset (-2 to +2) where it agrees best with the
     logo so far (the most frequent base of each column): at least 3 overlapping positions with at most one mismatch,
     ties to the smaller shift. One that fits nowhere is left out.
  3. Each placed 5-mer adds its weight R - 1 to its base in each column; columns holding under 30% of the largest
     column's weight are trimmed from the ends; each column is then made frequencies.

An experiment with fewer than --min-kmers placed 5-mers gets no logo. The logo summarises the enriched 5-mers; it is
not a binding model fitted to the peaks.

Writes docs/data/methylarg_eclip_logos.json: {meta, logos: [[experiment, source, uniprot_accession, cell_line,
[[A, C, G, U] per column], [[kmer, R, offset] per placed 5-mer]]]}.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
from collections import defaultdict
from pathlib import Path

BASES = 'ACGU'


def read(path: Path) -> list[dict]:
    with open(path, encoding='utf-8') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def place(cols: dict, kmer: str) -> int | None:
    """The offset of kmer against the logo so far (cols: position -> [A, C, G, U] weights), or None."""
    cons = {p: BASES[max(range(4), key=lambda b: w[b])] for p, w in cols.items()}
    best = None
    for off in (0, -1, 1, -2, 2):
        over = [(cons[off + i], c) for i, c in enumerate(kmer) if off + i in cons]
        if len(over) < 3:
            continue
        match = sum(a == b for a, b in over)
        if match < len(over) - 1:
            continue
        key = (match, -abs(off))
        if best is None or key > best[0]:
            best = (key, off)
    return None if best is None else best[1]


def logo(kmers: list[tuple[str, float]]) -> tuple[list[list[float]], list[list]]:
    cols: dict[int, list[float]] = {}
    used = []
    for kmer, r in kmers:
        off = 0 if not cols else place(cols, kmer)
        if off is None:
            continue
        for i, c in enumerate(kmer):
            cols.setdefault(off + i, [0.0, 0.0, 0.0, 0.0])[BASES.index(c)] += r - 1
        used.append([kmer, round(r, 2), off])
    if not cols:
        return [], used
    lo, hi = min(cols), max(cols)
    top = max(sum(w) for w in cols.values())
    while lo < hi and sum(cols.get(lo, [0])) < 0.3 * top:
        lo += 1
    while hi > lo and sum(cols.get(hi, [0])) < 0.3 * top:
        hi -= 1
    pfm = []
    for p in range(lo, hi + 1):
        w = cols.get(p, [0.0] * 4)
        t = sum(w) or 1.0
        pfm.append([round(x / t, 3) for x in w])
    shift = lo
    return pfm, [[k, r, o - shift] for k, r, o in used]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--inputs', default='inputs/nucleic_acid', type=Path)
    ap.add_argument('--out', default='docs/data/methylarg_eclip_logos.json', type=Path)
    ap.add_argument('--background', choices=['other', 'flank'], default='other')
    ap.add_argument('--min-r', type=float, default=1.5)
    ap.add_argument('--min-z', type=float, default=3.0)
    ap.add_argument('--min-count', type=int, default=10)
    ap.add_argument('--top', type=int, default=10)
    ap.add_argument('--min-kmers', type=int, default=3)
    a = ap.parse_args()

    exps = {}
    for r in read(a.inputs / 'encode_eclip_experiments.tsv'):
        exps[r['encode_experiment']] = ('ENCODE', r['uniprot_accession'], r['cell_line'])
    for r in read(a.inputs / 'encore_eclip_experiments.tsv'):
        exps[r['encore_experiment']] = ('ENCORE', r['uniprot_accession'], r['cell_line'])
    R, Z = 'R_' + a.background, 'z_' + a.background
    kmers = defaultdict(list)
    for name, col in (('encode_eclip_5mer_enrichment.tsv', 'encode_experiment'), ('encore_eclip_5mer_enrichment.tsv', 'encore_experiment')):
        for r in read(a.inputs / name):
            if float(r[R]) >= a.min_r and float(r[Z]) >= a.min_z and int(r['n_fg']) >= a.min_count:
                kmers[r[col]].append((r['kmer'], float(r[R])))
    out, none = [], 0
    for exp, (src, acc, cell) in sorted(exps.items()):
        ks = sorted(kmers.get(exp, []), key=lambda x: (-x[1], x[0]))[:a.top]
        pfm, used = logo(ks)
        if len(used) < a.min_kmers or not pfm:
            none += 1
            continue
        out.append([exp, src, acc, cell, pfm, used])
    meta = {
        'generated': dt.date.today().isoformat(),
        'script': 'scripts/43f_eclip_logos.py',
        'inputs': ['encode_eclip_5mer_enrichment.tsv', 'encore_eclip_5mer_enrichment.tsv', 'encode_eclip_experiments.tsv', 'encore_eclip_experiments.tsv'],
        'background': a.background,
        'rule': f'5-mers with R >= {a.min_r}, z >= {a.min_z} and at least {a.min_count} foreground occurrences against the {a.background} background, ranked by R, at most {a.top}, '
                f'aligned to the first (offsets -2 to +2, at least 3 overlapping positions, at most 1 mismatch) and weighted by R - 1; a logo needs {a.min_kmers} placed 5-mers',
        'fields': ['experiment', 'source (ENCODE or ENCORE)', 'uniprot_accession', 'cell line', '[[A, C, G, U] frequencies per column]', '[[kmer, R, column of its first base]]'],
        'experiments': len(exps), 'with_logo': len(out), 'without_logo': none,
    }
    a.out.write_text(json.dumps({'meta': meta, 'logos': out}, separators=(',', ':')), encoding='utf-8')
    print(f'{len(out)} of {len(exps)} experiments with a logo -> {a.out}')


if __name__ == '__main__':
    main()

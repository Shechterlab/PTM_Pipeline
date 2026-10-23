#!/usr/bin/env python3
"""InterPro domains of the explorer's proteins, for the protein map's "Domains (InterPro)" row.

From the InterPro intervals staged for the release (inputs/annotation/interpro_108.0_domain_intervals.tsv, the ones the
domain context of the analyses uses), for every canonical protein of the site data (docs/data/methylarg_sites_all.json,
which holds the lower-confidence tier's proteins too):

  - entries of type domain and repeat; a protein with neither keeps its homologous superfamilies;
  - each entry's member-database hits merged into one interval wherever they overlap (Pfam, SMART, PROSITE and the
    others place an entry's boundaries a few residues apart);
  - where two entries cover the same stretch (each 80% or more of the other), the one on more of these proteins labels
    it and the other is kept as a note (InterPro's specific subfamily entries beside the general domain).

Writes docs/data/methylarg_interpro.json: {meta, names: [[InterPro accession, name, type]], prot: {accession:
[[start, end, index into names, [indices of entries covering the same stretch]?], ...]}}.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
from collections import defaultdict
from pathlib import Path

KEEP = ('domain', 'repeat')


def merge(spans: list[tuple[int, int]]) -> list[list[int]]:
    out: list[list[int]] = []
    for a, b in sorted(spans):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--intervals', default='inputs/annotation/interpro_108.0_domain_intervals.tsv', type=Path)
    ap.add_argument('--sites-json', default='docs/data/methylarg_sites_all.json', type=Path)
    ap.add_argument('--out', default='docs/data/methylarg_interpro.json', type=Path)
    ap.add_argument('--release', default='InterPro 108.0')
    a = ap.parse_args()

    d = json.loads(a.sites_json.read_text(encoding='utf-8'))
    accs = {p[0] for p in d['proteins'] if '-' not in p[0]}
    hits: dict[str, dict[str, list]] = defaultdict(lambda: defaultdict(list))
    names: dict[str, tuple[str, str]] = {}
    with open(a.intervals, encoding='utf-8') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            acc = r['canonical_UniProtAC']
            if acc not in accs or not r['fragment_start'] or not r['fragment_end']:
                continue
            ipr, kind = r['interpro_accession'], r['interpro_type']
            names[ipr] = (r['interpro_name'], kind)
            hits[acc][ipr].append((int(r['fragment_start']), int(r['fragment_end'])))
    # how many proteins each entry is on: of two entries covering the same stretch, the more widely used one labels it
    use = defaultdict(int)
    for acc in hits:
        for ipr in hits[acc]:
            use[ipr] += 1
    index: dict[str, int] = {}
    prot: dict[str, list] = {}
    for acc in sorted(hits):
        entries = hits[acc]
        keep = [i for i in entries if names[i][1] in KEEP] or [i for i in entries if names[i][1] == 'homologous_superfamily']
        spans = sorted(((s, e, ipr) for ipr in keep for s, e in merge(entries[ipr])), key=lambda x: (-use[x[2]], x[2], x[0]))
        kept: list[list] = []
        for s, e, ipr in spans:
            # the same stretch under another entry (each covers 80% or more of the other) becomes a note on it
            same = next((r for r in kept if min(e, r[1]) - max(s, r[0]) + 1 >= 0.8 * max(e - s + 1, r[1] - r[0] + 1)), None)
            k = index.setdefault(ipr, len(index))
            if same is not None:
                if k != same[2] and k not in same[3]:
                    same[3].append(k)
                continue
            kept.append([s, e, k, []])
        if kept:
            prot[acc] = sorted([s, e, k] + ([also] if also else []) for s, e, k, also in kept)
    order = sorted(index, key=index.get)
    meta = {
        'generated': dt.date.today().isoformat(),
        'script': 'scripts/36c_export_interpro_domains.py',
        'release': a.release,
        'source': str(a.intervals),
        'rule': 'InterPro entries of type domain or repeat (homologous superfamilies for a protein with neither), each entry\'s member-database hits merged where they overlap; an entry covering the same stretch as a more widely used one (80% reciprocal overlap) is a note on it',
        'proteins': len(prot), 'of_proteins': len(accs), 'entries': len(order),
    }
    a.out.write_text(json.dumps({'meta': meta, 'names': [[i, names[i][0], names[i][1]] for i in order], 'prot': prot}, separators=(',', ':')), encoding='utf-8')
    print(f'{len(prot):,} of {len(accs):,} proteins with InterPro domains ({len(order):,} entries) -> {a.out}')


if __name__ == '__main__':
    main()

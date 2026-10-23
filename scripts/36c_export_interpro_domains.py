#!/usr/bin/env python3
"""InterPro domains and families for the protein map, and a light record of every reviewed protein without a site.

From the InterPro intervals staged for the release (inputs/annotation/interpro_108.0_domain_intervals.tsv, the ones the
domain context of the analyses uses, and interpro_108.0_all_entry_intervals.tsv for the family entries), per protein:

  - domains: entries of type domain and repeat; a protein with neither keeps its homologous superfamilies here;
  - families: entries of type family, and the homologous superfamilies not already shown as domains;
  - each entry's member-database hits merged into one interval wherever they overlap (Pfam, SMART, PROSITE and the
    others place an entry's boundaries a few residues apart);
  - where two entries of one kind cover the same stretch (each 80% or more of the other), the one on more reviewed
    human proteins labels it and the other is kept as a note (InterPro's specific subfamily entries beside the general
    one).

Writes two files:

  docs/data/methylarg_interpro.json, for the proteins of the site data (docs/data/methylarg_sites_all.json, which
  holds the lower-confidence tier's proteins too): {meta, names: [[InterPro accession, name, type]], prot: {accession:
  [[start, end, index into names, [indices of entries covering the same stretch]?], ...]}, fam: {accession: [...]}}.

  docs/data/methylarg_proteome.json and docs/data/proteome/pNN.json, for the reviewed human proteins of UniProt
  2026_03 with no site in the default view (docs/data/methylarg_sites.json), which the explorer shows on a light page.
  The index is small, for search: {meta, shards, prot: [[accession, gene, lower-confidence sites]]}. Each protein's
  record is in shard shard(accession) (shard() below; the explorer computes the same): {names, prot: {accession:
  [protein name, length, disorder [[start, end]], domains, families, CD-CODE condensate ids]}}, with domains and
  families indexing that shard's own names.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
from collections import defaultdict
from pathlib import Path

DOMAIN = ('domain', 'repeat')
SUPER = 'homologous_superfamily'


def merge(spans: list[tuple[int, int]]) -> list[list[int]]:
    out: list[list[int]] = []
    for a, b in sorted(spans):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def read_hits(path: Path, types: set[str], hits: dict, names: dict) -> None:
    with open(path, encoding='utf-8') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['interpro_type'] not in types or not r['fragment_start'] or not r['fragment_end']:
                continue
            ipr = r['interpro_accession']
            names[ipr] = (r['interpro_name'], r['interpro_type'])
            hits[r['canonical_UniProtAC']][ipr].append((int(r['fragment_start']), int(r['fragment_end'])))


def collapse(entries: dict, keep: list[str], use: dict, index: dict) -> list:
    """The kept entries' merged intervals, one per stretch: [[start, end, k, [also k]?], ...] by start."""
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
    return sorted([s, e, k] + ([also] if also else []) for s, e, k, also in kept)


def shard(acc: str, n: int) -> int:
    """The detail file of an accession: a hash of its characters, as the explorer computes it."""
    h = 0
    for c in acc:
        h = (h * 31 + ord(c)) % 1000003
    return h % n


def split(entries: dict, names: dict) -> tuple[list[str], list[str]]:
    """The entries of the domain row and of the family row."""
    dom = [i for i in entries if names[i][1] in DOMAIN]
    sup = [i for i in entries if names[i][1] == SUPER]
    fam = [i for i in entries if names[i][1] == 'family']
    return (dom, fam + sup) if dom else (sup, fam)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--intervals', default='inputs/annotation/interpro_108.0_domain_intervals.tsv', type=Path)
    ap.add_argument('--all-intervals', default='inputs/annotation/interpro_108.0_all_entry_intervals.tsv', type=Path)
    ap.add_argument('--entries', default='inputs/reference/uniprot_human_reviewed_2026_03_entries.tsv', type=Path)
    ap.add_argument('--disorder', default='inputs/annotation/mobidb_disorder_intervals_2026-03-30.tsv', type=Path)
    ap.add_argument('--cdcode', default='inputs/annotation/cdcode_protein_membership_v2.2.tsv', type=Path)
    ap.add_argument('--sites-json', default='docs/data/methylarg_sites.json', type=Path)
    ap.add_argument('--sites-all-json', default='docs/data/methylarg_sites_all.json', type=Path)
    ap.add_argument('--out', default='docs/data/methylarg_interpro.json', type=Path)
    ap.add_argument('--proteome-out', default='docs/data/methylarg_proteome.json', type=Path)
    ap.add_argument('--shards', type=int, default=64, help='detail files under <proteome-out dir>/proteome/')
    ap.add_argument('--release', default='InterPro 108.0')
    a = ap.parse_args()

    dall = json.loads(a.sites_all_json.read_text(encoding='utf-8'))
    dhi = json.loads(a.sites_json.read_text(encoding='utf-8'))
    site_accs = {p[0] for p in dall['proteins'] if '-' not in p[0]}
    listed = {p[0] for p in dhi['proteins']}
    lc_sites: dict[str, int] = defaultdict(int)
    for s in dall['sites']:
        acc = dall['proteins'][s[0]][0]
        if acc not in listed:
            lc_sites[acc] += 1

    hits: dict[str, dict[str, list]] = defaultdict(lambda: defaultdict(list))
    names: dict[str, tuple[str, str]] = {}
    read_hits(a.intervals, set(DOMAIN) | {SUPER}, hits, names)
    read_hits(a.all_intervals, {'family'}, hits, names)
    # how many reviewed proteins each entry is on: of two entries covering the same stretch, the more widely used one
    # labels it
    use: dict[str, int] = defaultdict(int)
    for acc in hits:
        for ipr in hits[acc]:
            use[ipr] += 1

    # the proteins of the site data
    index: dict[str, int] = {}
    prot: dict[str, list] = {}
    fam: dict[str, list] = {}
    for acc in sorted(site_accs & set(hits)):
        dom_ids, fam_ids = split(hits[acc], names)
        if dom_ids:
            prot[acc] = collapse(hits[acc], dom_ids, use, index)
        if fam_ids:
            fam[acc] = collapse(hits[acc], fam_ids, use, index)
    order = sorted(index, key=index.get)
    meta = {
        'generated': dt.date.today().isoformat(),
        'script': 'scripts/36c_export_interpro_domains.py',
        'release': a.release,
        'source': [str(a.intervals), str(a.all_intervals)],
        'rule': 'Domains: InterPro entries of type domain or repeat (homologous superfamilies for a protein with neither). Families: entries of type family and the other homologous superfamilies. '
                'Each entry\'s member-database hits merged where they overlap; an entry covering the same stretch as a more widely used one of its kind (80% reciprocal overlap) is a note on it',
        'proteins': len(prot), 'with_families': len(fam), 'of_proteins': len(site_accs), 'entries': len(order),
    }
    a.out.write_text(json.dumps({'meta': meta, 'names': [[i, names[i][0], names[i][1]] for i in order], 'prot': prot, 'fam': fam}, separators=(',', ':')), encoding='utf-8')
    print(f'{len(prot):,} of {len(site_accs):,} proteins with InterPro domains, {len(fam):,} with families ({len(order):,} entries) -> {a.out}')

    # every reviewed protein without a site in the default view
    idr: dict[str, list] = defaultdict(list)
    with open(a.disorder, encoding='utf-8') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            idr[r['canonical_UniProtAC']].append((int(r['fragment_start']), int(r['fragment_end'])))
    cd: dict[str, list] = defaultdict(list)
    with open(a.cdcode, encoding='utf-8') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['condensate_id'] not in cd[r['uniprotkb_ac']]:
                cd[r['uniprotkb_ac']].append(r['condensate_id'])
    shards = [{'index': {}, 'prot': {}} for _ in range(a.shards)]
    index_rows = []
    with open(a.entries, encoding='utf-8') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            acc = r['Entry']
            if acc in listed:
                continue
            sh = shards[shard(acc, a.shards)]
            dom_ids, fam_ids = split(hits[acc], names) if acc in hits else ([], [])
            index_rows.append([acc, r['Gene Names (primary)'] or '', lc_sites.get(acc, 0)])
            sh['prot'][acc] = [r['Protein names'].split(' (')[0] if r['Protein names'] else '', int(r['Length']), [list(x) for x in merge(idr.get(acc, []))],
                               collapse(hits[acc], dom_ids, use, sh['index']) if dom_ids else [],
                               collapse(hits[acc], fam_ids, use, sh['index']) if fam_ids else [], cd.get(acc, [])]
    pmeta = {
        'generated': dt.date.today().isoformat(),
        'script': 'scripts/36c_export_interpro_domains.py',
        'reference': str(a.entries), 'disorder': str(a.disorder), 'condensates': str(a.cdcode), 'interpro': meta['rule'],
        'rule': 'Reviewed human proteins of UniProt 2026_03 with no methylarginine site in the default view; lower-confidence sites counted from ' + str(a.sites_all_json),
        'proteins': len(index_rows), 'with_lower_confidence_sites': sum(1 for x in index_rows if x[2]),
        'fields': ['accession', 'gene', 'lower_confidence_sites'],
        'shard_fields': ['protein_name', 'length', 'disorder_intervals', 'interpro_domains', 'interpro_families', 'cdcode_condensates'],
        'shard_rule': 'h = 0; for each character c: h = (h * 31 + code(c)) mod 1000003; shard = h mod shards',
    }
    a.proteome_out.write_text(json.dumps({'meta': pmeta, 'shards': a.shards, 'prot': index_rows}, separators=(',', ':')), encoding='utf-8')
    sdir = a.proteome_out.parent / 'proteome'
    sdir.mkdir(exist_ok=True)
    for k, sh in enumerate(shards):
        order_k = sorted(sh['index'], key=sh['index'].get)
        (sdir / f'p{k:02d}.json').write_text(json.dumps({'names': [[i, names[i][0], names[i][1]] for i in order_k], 'prot': sh['prot']}, separators=(',', ':')), encoding='utf-8')
    print(f'{len(index_rows):,} reviewed proteins without a site ({pmeta["with_lower_confidence_sites"]:,} with lower-confidence sites only) -> {a.proteome_out} and {a.shards} files in {sdir}')


if __name__ == '__main__':
    main()

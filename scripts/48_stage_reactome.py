#!/usr/bin/env python3
"""Stage Reactome release 97 pathways of the reviewed canonical human proteome.

Reactome (reactome.org) maps UniProt accessions to the pathways their proteins take part in. Licence of the data: CC0
(https://reactome.org/license, section 1 c: "All data in the Reactome database and files derived from that data are
licensed under the Creative Commons Public Domain Dedication (CC0)"). Cite Ragueneau et al., The Reactome Knowledgebase
2026, Nucleic Acids Res (2025), doi:10.1093/nar/gkaf1223 (reactome.org/cite).

Downloads read (--reactome folder; https://download.reactome.org/97/, all species, filtered here to Homo sapiens)
  UniProt2Reactome_All_Levels.txt   accession, pathway, URL, name, evidence code, species: every level of the hierarchy
  UniProt2Reactome.txt              the same for the lowest-level pathway (diagram or subset) only
  ReactomePathways.txt              pathway, name, species
  ReactomePathwaysRelation.txt      parent pathway, child pathway
  api_events_toplevel_9606.json     optional: ContentService /data/pathways/top/9606, to check the top-level pathways

An isoform accession (P12345-2) is mapped to its entry (P12345); rows whose entry is not a reviewed entry of --fasta are
left out and counted. One row per entry and pathway: its evidence codes (TAS, IEA, as Reactome gives them; both when
Reactome lists the pair with each), whether the pathway is one of the entry's lowest-level pathways (UniProt2Reactome.txt),
and the Reactome accessions behind the row. The top-level pathways are those without a parent in the relation file.

Outputs (in --outdir)
  reactome_human_protein_pathways.tsv   canonical_UniProtAC, pathway_id, evidence, lowest_level, reactome_accessions
  reactome_human_pathways.tsv           every human pathway: name, parents, top-level pathway(s) (ids and names), depth
                                        (fewest steps below a top-level pathway), children, n_reviewed_proteins (all
                                        levels) and n_reviewed_proteins_lowest
  reactome_summary.json                 release, input MD5s, counts, what was left out and why
  manifest_rows.tsv                     the rows these files would take in inputs/MANIFEST.tsv (its columns)
  raw_files_rows.tsv                    the rows of the downloads for inputs/layers/raw_files.tsv (its columns)
Nothing outside --outdir is written: --manifest and --raw-files are read only, for their columns.

    python scripts/48_stage_reactome.py --reactome <download folder> --outdir <folder>
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict, deque
from pathlib import Path

import pandas as pd

from common import md5_file, parse_fasta, save_table

RELEASE = '97'
SPECIES = 'Homo sapiens'
PREFIX = 'R-HSA-'
BASE = f'https://download.reactome.org/{RELEASE}'
# The downloads: URL, the Last-Modified and ETag the server gave (2026-10-25). An ETag "md5-N" is Amazon S3's multipart
# tag (8 MiB parts), checked by recomputing it.
DOWNLOADS = {
    'UniProt2Reactome_All_Levels.txt': (f'{BASE}/UniProt2Reactome_All_Levels.txt', '2026-06-24 19:18:42 GMT', '6ff6a10c934261389750f2bb8dffb427-15'),
    'UniProt2Reactome.txt': (f'{BASE}/UniProt2Reactome.txt', '2026-06-24 19:18:40 GMT', '7d51554a754f6ce155c705d5e78d914c-6'),
    'ReactomePathways.txt': (f'{BASE}/ReactomePathways.txt', '2026-06-24 19:18:39 GMT', '63953262f948a43232e01ac18c9606d2'),
    'ReactomePathwaysRelation.txt': (f'{BASE}/ReactomePathwaysRelation.txt', '2026-06-24 19:18:40 GMT', '4cd76055c4ea7d18308d5993f160533d'),
    'api_events_toplevel_9606.json': ('https://reactome.org/ContentService/data/pathways/top/9606', '', ''),
}
MAP_COLUMNS = ['accession', 'pathway_id', 'url', 'name', 'evidence', 'species']
LICENCE = ('CC0 (https://reactome.org/license, 1 c: "All data in the Reactome database and files derived from that data are licensed '
           'under the Creative Commons Public Domain Dedication (CC0)")')
MANIFEST_LICENCE = 'CC0 1.0 (https://reactome.org/license, section 1 c, data)'
CITATION = 'Ragueneau et al., The Reactome Knowledgebase 2026, Nucleic Acids Res (2025), doi:10.1093/nar/gkaf1223; Reactome release 97'
OUT_PAIRS = 'reactome_human_protein_pathways.tsv'
OUT_PATHWAYS = 'reactome_human_pathways.tsv'
OUT_SUMMARY = 'reactome_summary.json'


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def s3_etag(path: Path, multipart: bool, part: int = 8 * 1024 * 1024) -> str:
    if not multipart:
        return md5_file(path)
    digests = []
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(part), b''):
            digests.append(hashlib.md5(block).digest())
    return hashlib.md5(b''.join(digests)).hexdigest() + f'-{len(digests)}'


def write_rows(path: Path, columns: list[str], rows: list[dict]) -> None:
    """A TSV with the given columns, written as common.upsert_tsv_rows writes rows (no quoting)."""
    lines = ['\t'.join(columns)]
    for r in rows:
        cells = [str(r.get(c, '')) for c in columns]
        if any(ch in c for c in cells for ch in '\t\n\r'):
            raise SystemExit(f'{path.name}: a value holds a tab or line break')
        lines.append('\t'.join(cells))
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def read_mapping(path: Path, canonical: dict, counts: dict) -> pd.DataFrame:
    """Human rows of a UniProt2Reactome file, with entry (isoforms mapped to their entry) and only reference entries."""
    t = pd.read_csv(path, sep='\t', header=None, names=MAP_COLUMNS, dtype=str, keep_default_na=False, quoting=3)
    h = t[t['species'] == SPECIES].copy()
    if not h['pathway_id'].str.startswith(PREFIX).all() or t.loc[t['species'] != SPECIES, 'pathway_id'].str.startswith(PREFIX).any():
        raise SystemExit(f'{path.name}: species and {PREFIX} pathway identifiers disagree')
    h['entry'] = h['accession'].str.split('-').str[0]
    iso = h['accession'].str.contains('-')
    keep = h['entry'].isin(canonical)
    counts.update({'rows_all_species': int(len(t)), 'rows_human': int(len(h)), 'human_accessions': int(h['accession'].nunique()),
                   'isoform_accessions': int(h.loc[iso, 'accession'].nunique()),
                   'isoform_accessions_of_reference_entries': int(h.loc[iso & keep, 'accession'].nunique()),
                   'reference_entries_reached_by_isoforms': int(h.loc[iso & keep, 'entry'].nunique()),
                   'reference_entries_reached_only_by_isoforms': len(set(h.loc[iso & keep, 'entry']) - set(h.loc[~iso & keep, 'entry'])),
                   'isoform_rows': int(iso.sum()),
                   'accessions_not_in_reference': int(h.loc[~keep, 'accession'].nunique()), 'rows_not_in_reference': int((~keep).sum()),
                   'rows_kept': int(keep.sum()), 'evidence_rows': dict(sorted(Counter(h.loc[keep, 'evidence']).items()))})
    return h[keep]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--reactome', required=True, help='folder with the four Reactome files (and optionally api_events_toplevel_9606.json)')
    ap.add_argument('--fasta', default='inputs/reference/uniprot_human_reviewed_2026_03.fasta')
    ap.add_argument('--manifest', default='inputs/MANIFEST.tsv', help='read only: its columns')
    ap.add_argument('--raw-files', default='inputs/layers/raw_files.tsv', help='read only: its columns')
    ap.add_argument('--downloaded', default='2026-10-25', help='the day the downloads were made')
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    src, out = Path(args.reactome), Path(args.outdir)
    canonical = parse_fasta(args.fasta)

    counts_all, counts_low = {}, {}
    allv = read_mapping(src / 'UniProt2Reactome_All_Levels.txt', canonical, counts_all)
    low = read_mapping(src / 'UniProt2Reactome.txt', canonical, counts_low)
    paths = pd.read_csv(src / 'ReactomePathways.txt', sep='\t', header=None, names=['pathway_id', 'name', 'species'], dtype=str,
                        keep_default_na=False, quoting=3)
    paths = paths[paths['species'] == SPECIES]
    if paths['pathway_id'].duplicated().any() or not paths['pathway_id'].str.startswith(PREFIX).all():
        raise SystemExit('ReactomePathways.txt: duplicated or non-human human pathway identifiers')
    name = dict(zip(paths['pathway_id'], paths['name']))
    rel = pd.read_csv(src / 'ReactomePathwaysRelation.txt', sep='\t', header=None, names=['parent', 'child'], dtype=str, keep_default_na=False)
    rel = rel[rel['parent'].str.startswith(PREFIX) | rel['child'].str.startswith(PREFIX)]
    if not (rel['parent'].str.startswith(PREFIX) & rel['child'].str.startswith(PREFIX)).all():
        raise SystemExit('ReactomePathwaysRelation.txt: a relation crosses species')
    unknown = (set(rel['parent']) | set(rel['child']) | set(allv['pathway_id']) | set(low['pathway_id'])) - set(name)
    if unknown:
        raise SystemExit(f'{len(unknown)} pathway identifiers are not in ReactomePathways.txt, e.g. {sorted(unknown)[:3]}')
    names_differ = int((allv['name'] != allv['pathway_id'].map(name)).sum() + (low['name'] != low['pathway_id'].map(name)).sum())

    # hierarchy
    parents, children = defaultdict(set), defaultdict(set)
    for p, c in zip(rel['parent'], rel['child']):
        parents[c].add(p)
        children[p].add(c)
    top = sorted(p for p in name if not parents[p])
    depth = {p: 0 for p in top}
    queue = deque(top)
    while queue:
        p = queue.popleft()
        for c in sorted(children[p]):
            if c not in depth:
                depth[c] = depth[p] + 1
                queue.append(c)
    if len(depth) != len(name):
        raise SystemExit('the pathway hierarchy has a cycle or a pathway no top-level pathway reaches')

    def ancestors(p: str, memo: dict = {}) -> frozenset:  # noqa: B006 - a cache
        if p not in memo:
            memo[p] = frozenset(set().union(*(ancestors(q) | {q} for q in parents[p])) if parents[p] else set())
        return memo[p]
    top_of = {p: sorted(a for a in ancestors(p) | {p} if a in set(top)) for p in name}

    # one row per entry and pathway
    ev = allv.groupby(['entry', 'pathway_id'])['evidence'].agg(lambda s: ';'.join(sorted(set(s))))
    accs = allv.groupby(['entry', 'pathway_id'])['accession'].agg(lambda s: ';'.join(sorted(set(s))))
    low_pairs = set(zip(low['entry'], low['pathway_id']))
    pairs = pd.DataFrame({'evidence': ev, 'reactome_accessions': accs}).reset_index().rename(columns={'entry': 'canonical_UniProtAC'})
    pairs['lowest_level'] = [(a, p) in low_pairs for a, p in zip(pairs['canonical_UniProtAC'], pairs['pathway_id'])]
    pairs = pairs[['canonical_UniProtAC', 'pathway_id', 'evidence', 'lowest_level', 'reactome_accessions']]
    pairs = pairs.sort_values(['canonical_UniProtAC', 'pathway_id']).reset_index(drop=True)
    pairs['lowest_level'] = pairs['lowest_level'].map({True: 'true', False: 'false'})

    # checks: the lowest-level pairs are among all levels, and all levels are the lowest-level pairs and their ancestors
    all_pairs = set(zip(pairs['canonical_UniProtAC'], pairs['pathway_id']))
    closure = {(a, q) for a, p in low_pairs for q in ancestors(p) | {p}}
    # pathway table
    n_all = pairs.groupby('pathway_id')['canonical_UniProtAC'].nunique()
    n_low = pairs[pairs['lowest_level'] == 'true'].groupby('pathway_id')['canonical_UniProtAC'].nunique()
    prow = []
    for p in sorted(name):
        prow.append({'pathway_id': p, 'name': name[p], 'parents': ';'.join(sorted(parents[p])), 'top_level': ';'.join(top_of[p]),
                     'top_level_names': ';'.join(name[t] for t in top_of[p]), 'depth': depth[p], 'n_children': len(children[p]),
                     'n_reviewed_proteins': int(n_all.get(p, 0)), 'n_reviewed_proteins_lowest': int(n_low.get(p, 0))})
    ptab = pd.DataFrame(prow)

    api_check = {}
    jp = src / 'api_events_toplevel_9606.json'
    if jp.exists():
        api = {x['stId'] for x in json.loads(jp.read_text(encoding='utf-8'))}
        api_check = {'content_service_top_level': len(api), 'same_as_relation_roots': api == set(top), 'md5': md5_file(jp)}

    save_table(pairs, out / OUT_PAIRS)
    save_table(ptab, out / OUT_PATHWAYS)
    sizes = ptab['n_reviewed_proteins']
    lsz = ptab['n_reviewed_proteins_lowest']
    summary = {
        'reactome': {'release': RELEASE, 'release_source': 'https://reactome.org/ContentService/data/database/version answered 97 on '
                     '2026-10-25; https://reactome.org/download-data rewrites its links to download.reactome.org/<that version>/. A folder '
                     'download.reactome.org/98/ (files dated 2026-10-22) exists but the site did not serve it as current',
                     'licence': LICENCE, 'citation': CITATION, 'downloaded': args.downloaded},
        'inputs': {n: {'url': u, 'last_modified': lm, 'md5': md5_file(src / n),
                       **({'etag': et, 'etag_matches': s3_etag(src / n, '-' in et) == et} if et else {})}
                   for n, (u, lm, et) in DOWNLOADS.items() if (src / n).exists()},
        'reference': {'fasta_md5': md5_file(args.fasta), 'reference_entries': len(canonical)},
        'all_levels_file': counts_all,
        'lowest_level_file': counts_low,
        'names_differing_from_ReactomePathways': names_differ,
        'pairs': int(len(pairs)), 'pairs_lowest_level': int((pairs['lowest_level'] == 'true').sum()),
        'pairs_by_evidence': dict(sorted(Counter(pairs['evidence']).items())),
        'pairs_lowest_level_not_in_all_levels': len(low_pairs - all_pairs),
        'pairs_all_levels_not_lowest_or_ancestor': len(all_pairs - closure),
        'pairs_ancestor_of_lowest_not_in_all_levels': len(closure - all_pairs),
        'entries_with_a_pathway': int(pairs['canonical_UniProtAC'].nunique()),
        'entries_without_a_pathway': len(canonical) - int(pairs['canonical_UniProtAC'].nunique()),
        'pathways_human': int(len(ptab)), 'top_level_pathways': len(top), 'top_level_check': api_check,
        'pathways_with_several_parents': int(sum(1 for p in name if len(parents[p]) > 1)),
        'pathways_under_several_top_level': int(sum(1 for p in name if len(top_of[p]) > 1)),
        'max_depth': int(max(depth.values())),
        'pathways_with_reviewed_proteins': int((sizes > 0).sum()),
        'pathways_10_to_500_reviewed_proteins_all_levels': int(((sizes >= 10) & (sizes <= 500)).sum()),
        'pathways_10_to_500_reviewed_proteins_lowest_level': int(((lsz >= 10) & (lsz <= 500)).sum()),
        'pathways_over_500_reviewed_proteins_all_levels': int((sizes > 500).sum()),
        'top_level_sizes': {name[t]: int(n_all.get(t, 0)) for t in top},
    }
    (out / OUT_SUMMARY).write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')

    mcols = open(args.manifest, encoding='utf-8').readline().rstrip('\n').split('\t')
    rcols = open(args.raw_files, encoding='utf-8').readline().rstrip('\n').split('\t')
    source = '; '.join(DOWNLOADS[n][0] for n in ('UniProt2Reactome_All_Levels.txt', 'UniProt2Reactome.txt', 'ReactomePathways.txt',
                                                  'ReactomePathwaysRelation.txt'))
    version = f'Reactome release {RELEASE} (files dated 2026-06-24), downloaded {args.downloaded}'
    how = ('scripts/48_stage_reactome.py: Homo sapiens rows; isoform accessions mapped to their entry; only entries of '
           'reference/uniprot_human_reviewed_2026_03.fasta')
    mrows = [
        {'file': f'layers/{OUT_PAIRS}', 'content': 'Reactome release 97 pathways of reviewed canonical human proteins, every level of '
         'the hierarchy: one row per entry and pathway, with Reactome\'s evidence codes (TAS, IEA), whether it is one of the entry\'s '
         'lowest-level pathways, and the Reactome accessions (entry or isoforms) behind it', 'rows': len(pairs),
         'derived_from': 'UniProt2Reactome_All_Levels.txt; UniProt2Reactome.txt', 'how_made_or_checked': how},
        {'file': f'layers/{OUT_PATHWAYS}', 'content': 'The Reactome release 97 human pathways: name, parents, top-level pathway(s), depth, '
         'children, and the reviewed canonical proteins in each (all levels, and lowest level)', 'rows': len(ptab),
         'derived_from': 'ReactomePathways.txt; ReactomePathwaysRelation.txt; UniProt2Reactome_All_Levels.txt; UniProt2Reactome.txt',
         'how_made_or_checked': how},
        {'file': f'layers/{OUT_SUMMARY}', 'content': 'Release, input MD5s, counts and what was left out of the Reactome staging',
         'rows': '', 'derived_from': '', 'how_made_or_checked': 'scripts/48_stage_reactome.py'},
    ]
    for r in mrows:
        p = out / r['file'].split('/')[-1]
        r.update({'bytes': p.stat().st_size, 'md5': md5_file(p), 'redistribution': 'git', 'source': source, 'version_or_date': version,
                  'citation': CITATION, 'licence': MANIFEST_LICENCE})
    write_rows(out / 'manifest_rows.tsv', mcols, mrows)
    rrows = []
    for n, (u, lm, et) in DOWNLOADS.items():
        p = src / n
        if p.exists():
            what = (f'Reactome release {RELEASE} (server date {lm}; S3 ETag {et} matches)' if et else
                    f'Reactome ContentService answer for release {RELEASE}, top-level human pathways (a check)')
            rrows.append({'file': n, 'url': u, 'bytes': p.stat().st_size, 'md5': md5_file(p), 'sha256': sha256_file(p),
                          'downloaded': args.downloaded, 'version': what, 'kept': 'release_asset'})
    bad = [n for n, v in summary['inputs'].items() if v.get('etag_matches') is False]
    if bad:
        raise SystemExit(f'the ETag of {bad} does not match: download again')
    write_rows(out / 'raw_files_rows.tsv', rcols, rrows)
    print(json.dumps({k: summary[k] for k in ('pairs', 'pairs_lowest_level', 'pairs_by_evidence', 'entries_with_a_pathway', 'pathways_human',
                                              'top_level_pathways', 'pathways_10_to_500_reviewed_proteins_all_levels')}, indent=2))


if __name__ == '__main__':
    main()

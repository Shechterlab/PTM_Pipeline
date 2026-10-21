#!/usr/bin/env python3
"""PANTHER protein classes of the methylated proteins against the other reviewed human proteins.

    python scripts/03b_panther_class_enrichment.py --integrated-sites WORK/integrated/human_arg_methyl_union_dedup_by_site.tsv \
        --uniprot-metadata inputs/reference/uniprot_human_reviewed_2026_03_entries.tsv \
        --protein-classes inputs/annotation/panther_19.0_protein_classes.tsv \
        --ontology inputs/annotation/panther_19.0_protein_class_ontology.tsv \
        --groups config/panther_class_groups.json --outdir WORK/panther_class

The universe is the reviewed entries of --uniprot-metadata; the targets are the entries with a site of the union. Each
protein PANTHER classifies has one top-level class, and config/panther_class_groups.json puts each top-level class in
one group, so the groups (with the proteins PANTHER leaves unclassified) and the top-level classes (likewise) each
partition the universe. A protein belongs to a class below the top level when PANTHER gives it that class or one under
it. Three analyses, each an over-representation test of every category with clusterProfiler::enricher (one-sided
hypergeometric test, fold enrichment, Benjamini-Hochberg within the analysis; scripts/enrichment.py), the background
being every reviewed protein (the unclassified ones form their own category):

  panther_group     the groups of the overview, and the unclassified proteins
  panther_class     the 24 top-level classes, and the unclassified proteins
  panther_subclass  every class below the top level with at least --min-proteins proteins in the universe

Writes panther_class_enrichment.tsv (one row per test, with the 2x2 counts and each category's parent: the group of a top-level class, the top-level class of a subclass) and
protein_panther_classes.tsv (one row per protein of the universe). Needs Rscript with clusterProfiler on PATH.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd

from enrichment import Job, run


def read(path) -> pd.DataFrame:
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


def tests(analysis: str, members: dict[str, set], parents: dict[str, str], ids: dict[str, str], universe: set, targets: set,
          workdir: Path, min_size: int = 1) -> list[dict]:
    """clusterProfiler::enricher over the categories, the background every protein of the universe."""
    sets = pd.DataFrame([(name, acc) for name, group in members.items() for acc in group], columns=['term', 'gene'])
    res = run([Job(analysis, targets, universe, sets=sets, whole_universe=True, min_size=min_size)], workdir)[analysis]
    rows = []
    for r in res.to_dict('records'):
        rows.append({'analysis': analysis, 'category': r['ID'], 'category_id': ids.get(r['ID'], ''), 'parent': parents.get(r['ID'], ''),
                     'a': r['a'], 'b': r['b'], 'c': r['c'], 'd': r['d'], 'fold_enrichment': r['FoldEnrichment'],
                     'log2_fold_enrichment': math.log2(r['FoldEnrichment']) if r['FoldEnrichment'] > 0 else None,
                     'p_value': r['pvalue'], 'q_value': r['p.adjust'], 'method': 'enricher'})
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--integrated-sites', required=True)
    ap.add_argument('--uniprot-metadata', required=True)
    ap.add_argument('--protein-classes', required=True)
    ap.add_argument('--ontology', required=True)
    ap.add_argument('--groups', required=True)
    ap.add_argument('--min-proteins', type=int, default=10)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    universe = set(read(args.uniprot_metadata)['Entry'])
    targets = set(read(args.integrated_sites)['canonical_UniProtAC']) & universe
    ontology = read(args.ontology).set_index('class_id')
    name = ontology['name'].to_dict()
    parent = ontology['parent_id'].to_dict()
    top_of = ontology['top_level_id'].to_dict()
    config = json.loads(Path(args.groups).read_text())
    unclassified = config['unclassified']
    group_of = {c: g['name'] for g in config['groups'] for c in g['classes']}
    tops = set(ontology.loc[ontology['parent_id'] == 'PC00000'].index)
    if set(group_of) != tops or sum(len(g['classes']) for g in config['groups']) != len(tops):
        raise SystemExit(f'{args.groups} must put each top-level class in exactly one group: '
                         f'missing {sorted(tops - set(group_of))}, unknown {sorted(set(group_of) - tops)}')

    def ancestors(c: str) -> list[str]:
        out = []
        while c in parent and c != 'PC00000':
            out.append(c)
            c = parent[c]
        return out

    classes = read(args.protein_classes).set_index('uniprot_accession')
    rows = []
    for acc in sorted(universe):
        assigned = classes.at[acc, 'protein_class_ids'].split(';') if acc in classes.index and classes.at[acc, 'protein_class_ids'] else []
        top = sorted({top_of[c] for c in assigned})
        if len(top) > 1:
            raise SystemExit(f'{acc} has several top-level classes: {top}')
        closure = sorted({a for c in assigned for a in ancestors(c)})
        # The most specific classes: those no other assigned class lies under.
        leaves = [c for c in assigned if not any(c != o and c in ancestors(o) for o in assigned)]
        rows.append({'uniprot_accession': acc, 'methylated': acc in targets,
                     'panther_group': group_of[top[0]] if top else unclassified,
                     'panther_class_id': top[0] if top else '', 'panther_class': name[top[0]] if top else '',
                     'panther_subclass_ids': ';'.join(c for c in closure if c not in tops),
                     'panther_subclasses': ';'.join(name[c] for c in closure if c not in tops),
                     'panther_leaf_classes': ';'.join(name[c] for c in leaves),
                     'panther_family': classes.at[acc, 'panther_family'] if acc in classes.index else '',
                     'panther_subfamily_name': classes.at[acc, 'subfamily_name'] if acc in classes.index else ''})
    proteins = pd.DataFrame(rows)
    proteins.to_csv(out / 'protein_panther_classes.tsv', sep='\t', index=False)

    by = lambda col: {k: set(g['uniprot_accession']) for k, g in proteins.groupby(col)}  # noqa: E731
    groups = by('panther_group')
    group_members = {g['name']: groups.get(g['name'], set()) for g in config['groups']} | {unclassified: groups.get(unclassified, set())}
    top_members = {name[c]: set(proteins.loc[proteins['panther_class_id'] == c, 'uniprot_accession']) for c in sorted(tops, key=lambda c: name[c])}
    top_members[unclassified] = groups.get(unclassified, set())
    sub_members = {}
    for c in ontology.index:
        if c in tops or c == 'PC00000':
            continue
        members = set(proteins.loc[proteins['panther_subclass_ids'].str.split(';').map(lambda ids, c=c: c in ids), 'uniprot_accession'])
        if len(members) >= args.min_proteins:
            sub_members[name[c]] = members
    if len(set(name[c] for c in ontology.index)) != len(ontology):
        raise SystemExit('two PANTHER classes share a name')
    ids = {v: k for k, v in name.items()}
    work = out / 'enricher'
    result = (tests('panther_group', group_members, {}, {}, universe, targets, work)
              + tests('panther_class', top_members, {name[c]: group_of[c] for c in tops}, ids, universe, targets, work)
              + tests('panther_subclass', sub_members, {name[c]: name[top_of[c]] for c in ontology.index if c not in tops and c != 'PC00000'},
                      ids, universe, targets, work, min_size=args.min_proteins))
    pd.DataFrame(result).to_csv(out / 'panther_class_enrichment.tsv', sep='\t', index=False)
    # The groups and the top-level classes must partition the universe. Checked on the categories themselves:
    # clusterProfiler leaves out a category without a methylated protein, so the tested rows need not add up.
    for analysis, members in (('panther_group', group_members), ('panther_class', top_members)):
        listed = [acc for group in members.values() for acc in group]
        if len(listed) != len(set(listed)) or set(listed) != universe:
            raise SystemExit(f'{analysis} does not partition the universe: {len(listed)} listed, {len(set(listed))} distinct, '
                             f'{len(universe)} in the universe, {len(set(listed) - universe)} outside it')
    classified = int((proteins['panther_class_id'] != '').sum())
    print(f'{len(universe):,} reviewed proteins, {len(targets):,} methylated; {classified:,} with a PANTHER class; '
          f'{len(group_members)} groups, {len(top_members)} top-level rows, {len(sub_members)} subclasses of at least {args.min_proteins} proteins, '
          f'{sum(r["analysis"] == "panther_subclass" for r in result)} of them with a methylated protein (clusterProfiler tests only those)')


if __name__ == '__main__':
    main()

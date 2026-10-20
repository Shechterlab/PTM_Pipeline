#!/usr/bin/env python3
"""Stage the PANTHER 19.0 protein classes of the reviewed human proteins (inputs/annotation/).

    python scripts/00_stage_panther_classes.py --classification PTHR19.0_human --protein-class Protein_Class_19.0 \
        --relationship Protein_class_relationship \
        --entries inputs/reference/uniprot_human_reviewed_2026_03_entries.tsv --outdir inputs/annotation

The three files are PANTHER's downloads (data CC BY 4.0, https://www.pantherdb.org/tou.jsp):
  https://data.pantherdb.org/ftp/sequence_classifications/19.0/PANTHER_Sequence_Classification_files/PTHR19.0_human
  https://data.pantherdb.org/PANTHER19.0/ontology/Protein_Class_19.0
  https://data.pantherdb.org/PANTHER19.0/ontology/Protein_class_relationship
(the protein class ontology itself was last revised for PANTHER 17.0, as its header says).

Writes
  panther_19.0_protein_classes.tsv         one row per reviewed entry of --entries that PANTHER classifies into a
                                           family: its family, subfamily and protein class ids (empty when PANTHER
                                           gives no protein class)
  panther_19.0_protein_class_ontology.tsv  the classes in use (lines PANTHER comments out are left out), each with
                                           its parent and its top-level class (the ancestor just below the root,
                                           PC00000 "protein class")
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = 'PC00000'


def read_ontology(protein_class: str, relationship: str) -> pd.DataFrame:
    classes = pd.read_csv(protein_class, sep='\t', comment='!', header=None, dtype=str, keep_default_na=False,
                          names=['class_id', 'code', 'name', 'definition'])
    classes = classes[(classes['class_id'] != '') & ~classes['class_id'].str.startswith('#')].copy()
    rel = pd.read_csv(relationship, sep='\t', comment='!', header=None, dtype=str, keep_default_na=False,
                      names=['child', 'child_name', 'parent', 'parent_name', 'order'])
    rel = rel[(rel['child'] != '') & ~rel['child'].str.startswith('#')]
    parent = dict(zip(rel['child'], rel['parent']))
    known = set(classes['class_id'])
    missing = sorted((set(parent) | set(parent.values())) - known)
    if missing:
        raise SystemExit(f'relationship names classes the class list does not define: {missing}')

    def top(c: str) -> str:
        seen = []
        while parent.get(c, ROOT) != ROOT:
            if c in seen:
                raise SystemExit(f'cycle in the class relationship at {c}')
            seen.append(c)
            c = parent[c]
        return c

    classes['parent_id'] = classes['class_id'].map(lambda c: parent.get(c, ''))
    classes['top_level_id'] = classes['class_id'].map(lambda c: '' if c == ROOT else top(c))
    return classes[['class_id', 'name', 'parent_id', 'top_level_id', 'code', 'definition']]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--classification', required=True, help='PTHR19.0_human')
    ap.add_argument('--protein-class', required=True, help='Protein_Class_19.0')
    ap.add_argument('--relationship', required=True, help='Protein_class_relationship')
    ap.add_argument('--entries', required=True, help='reviewed UniProt entries (Entry column)')
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    ontology = read_ontology(args.protein_class, args.relationship)
    ontology.to_csv(out / 'panther_19.0_protein_class_ontology.tsv', sep='\t', index=False)

    # Columns: gene identifier, UniProt accession, gene symbol, subfamily id, family name, subfamily name, GO molecular
    # function, biological process, cellular component, protein class, pathway.
    raw = pd.read_csv(args.classification, sep='\t', header=None, dtype=str, keep_default_na=False)
    if raw.shape[1] != 11:
        raise SystemExit(f'{args.classification}: expected 11 columns, found {raw.shape[1]}')
    reviewed = set(pd.read_csv(args.entries, sep='\t', dtype=str, keep_default_na=False)['Entry'])
    raw = raw[raw[1].isin(reviewed)]
    if raw[1].duplicated().any():
        raise SystemExit('an accession is classified twice')
    known = set(ontology['class_id'])
    ids = raw[9].map(lambda s: [x.split('#')[1] for x in s.split(';') if '#' in x])
    unknown = sorted({c for row in ids for c in row} - known)
    if unknown:
        raise SystemExit(f'protein classes not in the ontology: {unknown}')
    table = pd.DataFrame({
        'uniprot_accession': raw[1].values,
        'gene_symbol': raw[2].values,
        'panther_family': raw[3].str.split(':').str[0].values,
        'panther_subfamily': raw[3].values,
        'family_name': raw[4].values,
        'subfamily_name': raw[5].values,
        'protein_class_ids': ids.map(';'.join).values,
    }).sort_values('uniprot_accession')
    table.to_csv(out / 'panther_19.0_protein_classes.tsv', sep='\t', index=False)
    with_class = int((table['protein_class_ids'] != '').sum())
    print(f'{len(table):,} of {len(reviewed):,} reviewed entries in PANTHER families, {with_class:,} with a protein class; '
          f'{len(ontology):,} classes')


if __name__ == '__main__':
    main()

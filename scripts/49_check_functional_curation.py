#!/usr/bin/env python3
"""Check the curated functional experiments (inputs/curation/functional_experiments.tsv) and summarise them.

    python scripts/49_check_functional_curation.py [--table ...] [--fasta ...] [--sites RELEASE_SITE_TABLE]

Each row must have the columns of inputs/curation/README.md in that order, a new FE id, a reviewed accession of the
reference proteome, an arginine at its position on that sequence (or no position and a region), codes from the README,
a DOI and a PubMed id. A row read from the abstract alone may leave experiment_kind empty. With --sites, the summary
says which curated sites are sites of that release. Exits 1 on any failure, so a staging or build step can call it.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd

from common import parse_fasta

COLUMNS = ['experiment_id', 'uniprot_accession', 'gene', 'position', 'region', 'species', 'system', 'writer', 'experiment_kind', 'perturbation',
           'perturbation_note', 'readout', 'direction', 'finding', 'figure', 'reference', 'doi', 'pmid', 'read_from', 'curated_on', 'note']
KINDS = {'methylated_vs_unmethylated', 'site_mutation', 'prmt_inhibitor', 'prmt_knockdown', 'prmt_knockout', 'correlation', 'prediction'}
DIRECTIONS = {'increased', 'decreased', 'no_change'}
READ_FROM = {'full_text', 'abstract'}
SPECIES = {'human', 'mouse', 'human protein'}


def problems(t: pd.DataFrame, seqs: dict) -> list[str]:
    bad = []
    if list(t.columns) != COLUMNS:
        return [f'columns differ from the README: {list(t.columns)}']
    for dup in t.loc[t['experiment_id'].duplicated(), 'experiment_id']:
        bad.append(f'{dup}: id used twice')
    for r in t.to_dict('records'):
        e, acc, pos = r['experiment_id'], r['uniprot_accession'], r['position']
        why = []
        if not re.fullmatch(r'FE\d{4}', e):
            why.append('id is not FE and four digits')
        if acc not in seqs:
            why.append(f'{acc} is not a canonical entry of the reference')
        elif pos:
            if not pos.isdigit() or not 0 < int(pos) <= len(seqs[acc]) or seqs[acc][int(pos) - 1] != 'R':
                why.append(f'position {pos} is not an arginine of {acc}')
        elif not r['region']:
            why.append('neither a position nor a region')
        if r['experiment_kind'] not in KINDS and not (r['experiment_kind'] == '' and r['read_from'] == 'abstract'):
            why.append(f'experiment_kind {r["experiment_kind"]!r}')
        for col, allowed in (('direction', DIRECTIONS), ('read_from', READ_FROM), ('species', SPECIES)):
            if r[col] not in allowed:
                why.append(f'{col} {r[col]!r}')
        if not re.fullmatch(r'10\.\d{4,9}/\S+', r['doi']):
            why.append(f'doi {r["doi"]!r}')
        if not r['pmid'].isdigit():
            why.append(f'pmid {r["pmid"]!r}')
        if not re.fullmatch(r'\d{4}-\d{2}-\d{2}', r['curated_on']):
            why.append(f'curated_on {r["curated_on"]!r}')
        for col in ('perturbation', 'readout', 'finding', 'reference'):
            if not r[col]:
                why.append(f'{col} empty')
        bad += [f'{e}: {w}' for w in why]
    return bad


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--table', default='inputs/curation/functional_experiments.tsv')
    ap.add_argument('--fasta', default='inputs/reference/uniprot_human_reviewed_2026_03.fasta')
    ap.add_argument('--sites', default='', help="a release's site table, to say which curated sites are release sites")
    args = ap.parse_args()
    t = pd.read_csv(args.table, sep='\t', dtype=str, keep_default_na=False)
    bad = problems(t, parse_fasta(args.fasta))
    if bad:
        print('\n'.join(bad), file=sys.stderr)
        raise SystemExit(f'{len(bad)} problems in {args.table}')
    sited = t[t['position'] != '']
    summary = {'rows': len(t), 'proteins': t['uniprot_accession'].nunique(), 'papers': t['doi'].nunique(),
               'rows_with_a_site': len(sited), 'rows_read_from_abstract': int((t['read_from'] == 'abstract').sum()),
               'by_kind': t['experiment_kind'].replace('', 'not_read').value_counts().sort_index().to_dict()}
    if args.sites:
        release = set(pd.read_csv(args.sites, sep='\t', dtype=str, usecols=['site_id'])['site_id'])
        curated = sorted(set(sited['uniprot_accession'] + ':R' + sited['position']))
        summary['curated_sites'] = curated
        summary['curated_sites_not_in_release'] = [s for s in curated if s not in release]
        summary['release'] = Path(args.sites).name
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()

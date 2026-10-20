#!/usr/bin/env python3
"""Write releases/v1.0.0/: the published union (DeAngelo et al. 2026) as readable tables, for change tracking.

The published results are the folder PTM_results/ at git tag v1.0.0 (commit 0574c61). This script reads four of its
tables and writes them, unchanged in content, as one site table and one source-record table:

    git worktree add /tmp/v100 v1.0.0
    python scripts/build_release_v1_0_0.py --published-root /tmp/v100/PTM_results --outdir releases/v1.0.0

Nothing is recomputed: every value is the published one, renamed to the column names of later releases.
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd

SOURCE_TABLE = {'maron_s5': 'maron2021_s5', 'prometheus_mmc2': 'prometheusdb_mmc2', 'human_ptm_master': 'ptm_master'}


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, low_memory=False)


def clean(value: str) -> str:
    return '' if str(value).strip().lower() in ('nan', 'none', '<na>') else str(value).strip()


def site_key(acc: str, pos: str) -> str:
    return f'{acc}:R{int(float(pos))}'


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--published-root', required=True, help='PTM_results/ at tag v1.0.0')
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    root, outdir = Path(args.published_root), Path(args.outdir)
    (outdir / 'supporting').mkdir(parents=True, exist_ok=True)

    union = read(root / 'integrated' / 'human_arg_methyl_union_dedup_by_site.tsv')
    rows = read(root / 'remapped' / 'human_arg_methyl_source_rows_with_remap.tsv')
    disorder = read(root / 'disorder_context' / 'sites_with_disorder_context.tsv')
    domain = read(root / 'domain_context' / 'sites_with_domain_context.tsv')
    dis = {site_key(a, p): c for a, p, c in zip(disorder['canonical_UniProtAC'], disorder['corrected_position'], disorder['disorder_context_class'])}
    dom = {site_key(a, p): c for a, p, c in zip(domain['canonical_UniProtAC'], domain['corrected_position'], domain['domain_context_class'])}

    sites = pd.DataFrame({
        'site_id': [site_key(a, p) for a, p in zip(union['canonical_UniProtAC'], union['corrected_position'])],
        'uniprot_accession': union['canonical_UniProtAC'],
        'gene': union['substrate_genename'].map(clean),
        'protein_name': union['protein_name'].map(clean),
        'position': union['corrected_position'].astype(float).astype(int),
        'methyl_states_reported': union['arg_methyl_states_observed'].map(lambda s: '; '.join(x for x in clean(s).split(';') if x)),
        'source_families': union['source_families_exact'].map(lambda s: '; '.join(x for x in clean(s).split(';') if x)),
        'n_records': union['n_rows_union'].astype(int),
    })
    sites['disorder'] = sites['site_id'].map(dis).fillna('')
    sites['domain'] = sites['site_id'].map(dom).fillna('')
    sites = sites.sort_values(['uniprot_accession', 'position']).reset_index(drop=True)

    resolved = rows['is_position_resolved'] == 'True'
    records = pd.DataFrame({
        'record_id': range(1, len(rows) + 1),
        'source_table': rows['source_table'].map(SOURCE_TABLE),
        'source': rows['source'].map(clean),
        'reference': rows['reference'].map(clean),
        'pmid': rows['pmid'].map(clean),
        'reported_accession': rows['substrate_UniProtAC_original'].map(clean),
        'reported_gene': rows['substrate_genename'].map(clean),
        'reported_position': rows['position'].map(clean),
        'reported_methyl_state': rows['arg_methyl_state'].map(clean),
        'peptide_sequence': rows['peptide_sequence'].map(clean),
        'seq_window': rows['seq_window'].map(clean),
        'marked_peptides': rows['peptides'].map(clean),
        'placement_code': rows['remap_status'],
        'placed_site': [site_key(a, p) if r else '' for a, p, r in zip(rows['corrected_canonical_UniProtAC'], rows['corrected_position'], resolved)],
    })

    sites.to_csv(outdir / 'human_methylarginine_sites_v1.0.0.tsv', sep='\t', index=False, lineterminator='\n')
    records.to_csv(outdir / 'supporting' / 'source_records_v1.0.0.tsv', sep='\t', index=False, lineterminator='\n')
    placed = set(records['placed_site']) - {''}
    if placed != set(sites['site_id']):
        raise SystemExit('the placed records and the union disagree')
    print({'sites': len(sites), 'proteins': sites['uniprot_accession'].nunique(), 'records': len(records),
           'records_placed': int(resolved.sum())})


if __name__ == '__main__':
    main()

#!/usr/bin/env python
"""Build the per-gene table used by the disease-ontology, ClinGen and TDP-43 tests.

One row per reviewed human UniProt entry with an Entrez Gene ID. Columns:
methylarginine membership (any evidence in the repo's site union), site count,
protein length, PaxDb whole-organism abundance, and an MS-detectable flag
(detected in any PaxDb 5.0 cell-line proteome, or carrying a methylarginine site).
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get('PTM_DATA_ROOT', ROOT / 'PTM_data'))
RESULTS = Path(os.environ.get('PTM_RESULTS_ROOT', ROOT / 'PTM_results'))

# PaxDb 5.0 cell-line proteomes used for the MS-detectable background.
CELL_LINES = ['A549', 'GAMG', 'HEK293', 'HeLa', 'HepG2', 'Jurkat', 'K562', 'LnCap', 'MCF7', 'RKO', 'U2OS']
CELLLINE_FILES = (
    [f'9606-{c}_Geiger_2012.txt' for c in CELL_LINES]
    + [f'9606-iBAQ_{c}_Geiger_2012_uniprot.txt' for c in CELL_LINES]
    + [
        '9606-hela_Nagaraj_2011.txt',
        '9606-hela_Nagaraj_2011_cpc.txt',
        '9606-hela_Nagaraj_2011_iBAQ.txt',
        '9606-PXD000271_Kulak_Nat_Methods_2014_HeLa.txt',
        '9606-PXD145877_Mueller_Nature_2020_Homo_sapiens_HeLa.txt',
    ]
)


def read_paxdb(path: Path) -> pd.DataFrame:
    # Some files carry a third (raw count) column; only the ID and ppm abundance are used.
    df = pd.read_csv(path, sep='\t', comment='#', header=None, usecols=[0, 1], names=['string_id', 'abundance'])
    return df[df['abundance'] > 0]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--uniprot', type=Path, default=DATA / 'context' / 'uniprot_human_reviewed.tsv',
                    help='UniProt reviewed human TSV: Entry, Entry Name, Gene Names (primary), Length, GeneID')
    ap.add_argument('--union', type=Path, default=RESULTS / 'integrated' / 'human_arg_methyl_union_dedup_by_site.tsv')
    ap.add_argument('--paxdb-dir', type=Path, default=DATA / 'paxdb')
    ap.add_argument('--out', type=Path, default=RESULTS / 'disease_gene_tests' / 'gene_table.tsv')
    args = ap.parse_args()

    up = pd.read_csv(args.uniprot, sep='\t')
    up.columns = ['uniprot', 'entry_name', 'symbol', 'length', 'geneid']
    up = up.dropna(subset=['geneid'])
    up['entrez'] = up['geneid'].str.split(';').str[0].astype(int)
    acc2entrez = dict(zip(up['uniprot'], up['entrez']))
    sym2entrez = up.dropna(subset=['symbol']).drop_duplicates('symbol').set_index('symbol')['entrez'].to_dict()
    # One row per Entrez ID (keep the longest reviewed protein).
    up = up.sort_values('length', ascending=False).drop_duplicates('entrez').drop(columns='geneid')

    # Sites are assigned to genes by reviewed accession; unreviewed (TrEMBL) accessions fall back
    # to their gene symbol, matching the symbol-based mapping of the existing disease analysis.
    union = pd.read_csv(args.union, sep='\t', usecols=['canonical_UniProtAC', 'substrate_genename'],
                        low_memory=False)
    by_acc = union['canonical_UniProtAC'].map(acc2entrez)
    by_sym = union['substrate_genename'].map(sym2entrez)
    union['entrez'] = by_acc.fillna(by_sym)
    union['mapped_by'] = np.where(by_acc.notna(), 'reviewed_accession',
                                  np.where(by_sym.notna(), 'gene_symbol', 'unmapped'))
    sites = union.dropna(subset=['entrez']).groupby('entrez').size().rename('arg_methyl_sites')
    sites.index = sites.index.astype(int)
    up = up.merge(sites, left_on='entrez', right_index=True, how='left')
    up['arg_methyl_sites'] = up['arg_methyl_sites'].fillna(0).astype(int)
    up['methyl_any'] = up['arg_methyl_sites'] > 0
    prot_map = union.drop_duplicates('canonical_UniProtAC')['mapped_by'].value_counts()

    links = pd.read_csv(args.paxdb_dir / 'human_paxdb_uniprot_links.tsv', sep='\t', header=None,
                        names=['string_id', 'entry_name'])
    whole = read_paxdb(args.paxdb_dir / '9606-WHOLE_ORGANISM-integrated.txt').merge(links, on='string_id')
    whole = whole.groupby('entry_name')['abundance'].max()
    up['paxdb_whole_organism_ppm'] = up['entry_name'].map(whole).fillna(0.0)
    up['log10_abundance'] = np.log10(up['paxdb_whole_organism_ppm'] + 0.01)
    up['log10_length'] = np.log10(up['length'])

    detected = set()
    per_file = []
    for fname in CELLLINE_FILES:
        df = read_paxdb(args.paxdb_dir / fname).merge(links, on='string_id')
        detected |= set(df['entry_name'])
        per_file.append((fname, df['entry_name'].nunique()))
    up['detected_cellline'] = up['entry_name'].isin(detected)
    up['ms_detectable'] = up['detected_cellline'] | up['methyl_any']

    up = up.sort_values('symbol').reset_index(drop=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    up.to_csv(args.out, sep='\t', index=False)

    summary = pd.DataFrame([
        ('reviewed_uniprot_with_entrez', len(up)),
        ('union_canonical_proteins', int(union['canonical_UniProtAC'].nunique())),
        ('union_proteins_mapped_by_reviewed_accession', int(prot_map.get('reviewed_accession', 0))),
        ('union_proteins_mapped_by_gene_symbol', int(prot_map.get('gene_symbol', 0))),
        ('union_proteins_unmapped', int(prot_map.get('unmapped', 0))),
        ('methyl_any_genes', int(up['methyl_any'].sum())),
        ('detected_in_paxdb_celllines', int(up['detected_cellline'].sum())),
        ('ms_detectable_background', int(up['ms_detectable'].sum())),
        ('methyl_any_not_detected_in_celllines', int((up['methyl_any'] & ~up['detected_cellline']).sum())),
    ] + [(f'paxdb_{f}', n) for f, n in per_file], columns=['metric', 'value'])
    summary.to_csv(args.out.with_name('gene_table_summary.tsv'), sep='\t', index=False)
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()

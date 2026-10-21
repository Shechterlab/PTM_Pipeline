#!/usr/bin/env python3
"""Stage Human Protein Atlas RNA expression for the reviewed human entries (inputs/annotation/).

    python scripts/00_stage_hpa_rna_expression.py --celline rna_celline.tsv.zip --tissue rna_tissue_consensus.tsv.zip \
        --atlas proteinatlas.tsv.zip --entries inputs/reference/uniprot_human_reviewed_2026_03_entries.tsv \
        --version 25.1 --out inputs/annotation/hpa_25.1_rna_expression.tsv

The three files are HPA downloads (CC BY 4.0, https://www.proteinatlas.org/about/licence):
  https://www.proteinatlas.org/download/tsv/rna_celline.tsv.zip            nTPM of each gene in 1,206 cell lines
  https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip   consensus nTPM of each gene in 51 tissues
  https://www.proteinatlas.org/download/proteinatlas.tsv.zip               Ensembl gene to UniProt accession(s)

Writes one row per reviewed accession that HPA maps to a gene with RNA data: its gene, the median and mean nTPM over
the cell lines and the mean and maximum consensus nTPM over the tissues. Where HPA maps one accession to several genes,
the row is the gene with the highest cell-line median (genes lists them all). The disease family test
(scripts/14) stratifies by the cell-line median, the measure that best separates the proteins found methylated from the
others: most methylarginine sites were found by mass spectrometry of cell lines.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import zipfile
from pathlib import Path

import pandas as pd


def read_member(path: str, usecols: list[str]) -> tuple[pd.DataFrame, str]:
    raw = Path(path).read_bytes()
    with zipfile.ZipFile(io.BytesIO(raw)) as zf:
        names = [n for n in zf.namelist() if n.endswith('.tsv')]
        if len(names) != 1:
            raise SystemExit(f'expected one .tsv in {path}, found {names}')
        table = pd.read_csv(zf.open(names[0]), sep='\t', dtype=str, keep_default_na=False, usecols=usecols)
        info = zf.getinfo(names[0])
    return table, f'{names[0]} ({info.date_time[0]:04d}-{info.date_time[1]:02d}-{info.date_time[2]:02d}), sha256 {hashlib.sha256(raw).hexdigest()}'


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--celline', required=True, help='rna_celline.tsv.zip')
    ap.add_argument('--tissue', required=True, help='rna_tissue_consensus.tsv.zip')
    ap.add_argument('--atlas', required=True, help='proteinatlas.tsv.zip')
    ap.add_argument('--entries', required=True, help='reviewed UniProt entries (Entry column)')
    ap.add_argument('--version', required=True, help='the HPA version of the downloads, e.g. 25.1')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    cell, cell_src = read_member(args.celline, ['Gene', 'Cell line', 'nTPM'])
    tissue, tissue_src = read_member(args.tissue, ['Gene', 'Tissue', 'nTPM'])
    atlas, atlas_src = read_member(args.atlas, ['Gene', 'Ensembl', 'Uniprot'])
    for t in (cell, tissue):
        t['nTPM'] = t['nTPM'].astype(float)
    c = cell.groupby('Gene')['nTPM'].agg(celline_median_ntpm='median', celline_mean_ntpm='mean', cell_lines='size')
    t = tissue.groupby('Gene')['nTPM'].agg(tissue_mean_ntpm='mean', tissue_max_ntpm='max')
    genes = c.join(t, how='outer').reset_index().rename(columns={'Gene': 'ensembl_gene'})
    genes = genes.merge(atlas.rename(columns={'Ensembl': 'ensembl_gene', 'Gene': 'gene'}), on='ensembl_gene', how='left')
    genes['uniprot_accession'] = genes['Uniprot'].fillna('').str.split(', ')
    rows = genes.explode('uniprot_accession')
    reviewed = set(pd.read_csv(args.entries, sep='\t', dtype=str, keep_default_na=False)['Entry'])
    rows = rows[rows['uniprot_accession'].isin(reviewed)]
    all_genes = rows.groupby('uniprot_accession')['ensembl_gene'].agg(lambda s: ';'.join(sorted(set(s))))
    rows = rows.sort_values(['uniprot_accession', 'celline_median_ntpm'], ascending=[True, False], na_position='last')
    out = rows.drop_duplicates('uniprot_accession').copy()
    out['genes'] = out['uniprot_accession'].map(all_genes)
    out['cell_lines'] = out['cell_lines'].fillna(0).astype(int)
    for col in ('celline_median_ntpm', 'celline_mean_ntpm', 'tissue_mean_ntpm', 'tissue_max_ntpm'):
        out[col] = out[col].round(3)
    out = out[['uniprot_accession', 'gene', 'ensembl_gene', 'genes', 'celline_median_ntpm', 'celline_mean_ntpm', 'cell_lines',
               'tissue_mean_ntpm', 'tissue_max_ntpm']].sort_values('uniprot_accession')
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep='\t', index=False, lineterminator='\n')
    print(f'HPA {args.version}: {len(out):,} of {len(reviewed):,} reviewed entries with RNA data '
          f'({int(out["celline_median_ntpm"].notna().sum()):,} with cell-line nTPM, {int((out["genes"].str.count(";") > 0).sum())} mapped to several genes)')
    print(f'  {cell_src}\n  {tissue_src}\n  {atlas_src}')


if __name__ == '__main__':
    main()

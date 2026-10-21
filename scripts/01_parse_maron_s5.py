"""Parse Maron et al. 2021 Table S5 (compiled methylarginine sites) into one row per source record.

Two faults of the table are corrected, and each corrected row says so (source_correction; the reported values stay in
*_reported columns):

- nan_restored: 14 Larsen 2016 windows are 28 residues long where the study's windows are 31, each missing the
  tripeptide NAN (no Larsen window holds it; the letters were read as a missing value when the table was made). With
  --canonical-fasta, a window is restored when inserting NAN at exactly one place gives the 31-residue window centred
  on the reported position of the reported entry.
- accession_na: 4 rows give no accession ("NA", the table's mark for a missing value), only a gene symbol and a
  window. Written to a separate file
  (maron_s5_rows_accession_na.tsv), they are placed from their gene symbol and window, after every other record, so
  the record numbers of the published version stay the same.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from common import (
    accession_is_isoform,
    canonical_gene_name,
    canonicalize_uniprot_accession,
    normalize_accession,
    save_table,
)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--exclude-pattern', default='')
    ap.add_argument('--canonical-fasta', default='', help='reference sequences, to restore the windows that lost NAN')
    args = ap.parse_args()

    df = pd.read_excel(args.input, sheet_name=0)
    out = pd.DataFrame()
    out['substrate_UniProtAC_original'] = df['Accession'].astype(str).str.strip()
    out['substrate_UniProtAC'] = out['substrate_UniProtAC_original'].map(normalize_accession)
    out['canonical_UniProtAC'] = out['substrate_UniProtAC'].map(canonicalize_uniprot_accession)
    out['accession_is_isoform'] = out['substrate_UniProtAC'].map(accession_is_isoform)
    out['substrate_genename'] = df['Gene Symbol'].map(canonical_gene_name)
    out['position'] = pd.to_numeric(df['Position In Protein'], errors='coerce')
    out['residue'] = 'R'
    out['site'] = 'R' + out['position'].astype('Int64').astype(str)
    out['peptide_sequence'] = df['Peptide Sequence'].astype(str)
    out['reference'] = df['Reference'].astype(str)
    out['source'] = 'Maron2021_TableS5'
    out['source_family'] = 'Maron2021'
    out['ptm_type'] = 'METHYLATION'
    out['ptm_group'] = 'Arg methylation'
    out['organism'] = 'Homo sapiens (Human)'
    out['site_key'] = out['substrate_UniProtAC'] + ':' + out['site']
    out['site_key_canonical_naive'] = out['canonical_UniProtAC'] + ':' + out['site']
    out['is_arg_methyl'] = True
    out['data_layer'] = 'external_methyl_table'
    out['arg_methyl_state'] = ''
    out['arg_methyl_state_known'] = False
    out['peptide_sequence_reported'] = ''
    out['source_correction'] = ''
    if args.canonical_fasta:
        restore_nan(out, read_fasta(Path(args.canonical_fasta)))
    na = df['Accession'].isna() & out['position'].notna() & out['substrate_genename'].ne('')
    na_rows = out[na].copy()
    na_rows['source_correction'] = 'accession_na'
    na_rows['substrate_UniProtAC_original'] = ''
    na_rows['substrate_UniProtAC'] = ''
    na_rows['canonical_UniProtAC'] = ''
    na_rows['site_key'] = ':' + na_rows['site']
    na_rows['site_key_canonical_naive'] = na_rows['site_key']
    out = out[out['position'].notna() & out['substrate_UniProtAC'].ne('') & ~na].copy()
    if args.exclude_pattern:
        out = out[~out['reference'].str.contains(args.exclude_pattern, case=False, na=False)].copy()
    out['position'] = out['position'].astype(int)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(out, outdir / 'maron_s5_arg_methyl_sites.tsv')
    if len(na_rows):
        na_rows['position'] = na_rows['position'].astype(int)
        save_table(na_rows, outdir / 'maron_s5_rows_accession_na.tsv')
    print(f"{int(out['source_correction'].eq('nan_restored').sum())} windows with NAN restored; {len(na_rows)} rows with accession NA")
    summary = pd.DataFrame([{
        'rows': len(out),
        'proteins': out['substrate_UniProtAC'].nunique(),
        'unique_sites': out['site_key'].nunique(),
    }])
    save_table(summary, outdir / 'summary.tsv')
    print(summary.to_dict(orient='records')[0])


def read_fasta(path: Path) -> dict[str, str]:
    seqs, acc = {}, None
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            acc = line[1:].split()[0].split('|')[1] if '|' in line else line[1:].split()[0]
            seqs[acc] = ''
        elif acc:
            seqs[acc] += line.strip()
    return seqs


def restore_nan(out: pd.DataFrame, seqs: dict[str, str]) -> None:
    """Larsen 2016 windows of Table S5 that lost the tripeptide NAN: restored in place (see the module docstring)."""
    larsen = out['reference'].str.contains('Larsen', na=False)
    full = out.loc[larsen, 'peptide_sequence'].str.len().mode().iloc[0]
    for i in out.index[larsen & out['peptide_sequence'].str.len().eq(full - 3)]:
        pep, pos, acc = out.at[i, 'peptide_sequence'], out.at[i, 'position'], out.at[i, 'canonical_UniProtAC']
        seq, flank = seqs.get(acc, ''), (full - 1) // 2
        if not seq or pd.isna(pos) or pos - 1 - flank < 0:
            continue
        want = seq[int(pos) - 1 - flank:int(pos) + flank]
        fits = [k for k in range(len(pep) + 1) if pep[:k] + 'NAN' + pep[k:] == want]
        if len(fits) == 1:
            out.at[i, 'peptide_sequence_reported'] = pep
            out.at[i, 'peptide_sequence'] = want
            out.at[i, 'source_correction'] = 'nan_restored'


if __name__ == '__main__':
    main()

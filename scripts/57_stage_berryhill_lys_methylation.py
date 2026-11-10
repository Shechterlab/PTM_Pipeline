#!/usr/bin/env python3
"""Stage the human lysine methylation sites of Berryhill et al. 2023 as a methyllysine table (inputs/marks/).

    python scripts/57_stage_berryhill_lys_methylation.py --xlsx PTM_data/marks/berryhill2023_41598_2022_27175_MOESM2_ESM.xlsx \
        --fasta inputs/reference/uniprot_human_reviewed_2026_03.fasta --out-dir inputs/marks --downloaded 2026-11-10

Source: Berryhill CA, Hanquier JN, Doud EH, Cordeiro-Spinetti E, Dickson BM, Rothbart SB, Mosley AL. Global lysine
methylome profiling using systematically characterized affinity reagents. Sci Rep 13:377 (2023),
doi:10.1038/s41598-022-27175-x. Licence: article CC BY 4.0; the supplementary files carry no separate notice; included
by the owner's decision (2026-11-10), without waiting for the authors' confirmation. Cite the paper wherever the data is
described.

Reads only Supplementary Table 5 (sheet "Sup Table 5" of Supplementary Information 2, MOESM2_ESM.xlsx): one row per
site, cell line and enrichment, with Kme_Site, Site.Probability, Peptide.Sequence, Position.in.Peptide, Protein.Accession,
Position.in.Protein, Modification (m1, m2, m3), Cell.Type and Enrichment. The other columns (PSM counts, motif,
description) and the other sheets are not used; Sup Tables 4 and 7 carry PhosphoSitePlus information and are left out.
The paper filters peptides "for high confidence identification and modification site probability" but gives no number
in the text; the table's own Site.Probability is 90.1-100, with "NA" for the PTMScan Rep 2 rows, which carry no
localisation probability. The column is kept (berryhill_site_probability, the highest value
over the rows of a site; empty when every row is NA), not filtered on.

Each row is checked on the reference proteome (--fasta, canonical reviewed accessions):
  kept       the canonical sequence has K at Position.in.Protein and the peptide lies on the canonical sequence at the
             offset that puts its K (Position.in.Peptide) there
  relocated  the position or peptide does not fit, but the peptide occurs exactly once in the canonical sequence and has
             K at Position.in.Peptide: the site moves to that K (the source position is kept in source_position;
             e.g. histone numbering with or without the initiator Met)
  rejected   everything else, with a reason: accession_not_in_reference, isoform_accession, bad_position, bad_state,
             peptide_position_not_K, peptide_not_unique, peptide_not_in_canonical
States: m1 me1, m2 me2, m3 me3. Rows of one site (accession, position, state) from several cell lines and enrichments
are one row, with berryhill_cell_lines, berryhill_enrichment and peptides ';'-joined. Columns: the ten of the other mark tables
(evidence empty, evidence_class experimental (mass spectrometry), pmids empty, source 'Berryhill 2023') and
berryhill_site_probability (the highest Site.Probability over the rows of the site), berryhill_cell_lines,
berryhill_enrichment, peptides, relocated (yes/no), source_position. The release (scripts/56) takes the first three; the
peptides and the relocation flags stay in inputs/marks/.
The staged table is named `..._sites.tsv`, so that scripts/56's glob `inputs/marks/*_methyllysine_sites.tsv` picks it up:
the release merges it with the other sources.
Outputs are sorted, so running again gives the same bytes. The download and the staged files have rows in
inputs/marks/raw_files.tsv and inputs/MANIFEST.tsv (the download itself stays out of Git).
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from datetime import date
from pathlib import Path

import openpyxl
import pandas as pd

from common import md5_file, parse_fasta, refresh_manifest_row, update_manifest, upsert_tsv_rows

URL = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-022-27175-x/MediaObjects/41598_2022_27175_MOESM2_ESM.xlsx'
ARTICLE_URL = 'https://www.nature.com/articles/s41598-022-27175-x'
SHEET = 'Sup Table 5'
CITATION = 'Berryhill et al., Sci Rep 13:377 (2023), doi:10.1038/s41598-022-27175-x'
LICENSE = "article CC BY 4.0; supplementary files carry no separate notice; included by the owner's decision (2026-11-10)"
SOURCE = 'Berryhill 2023'
SOURCE_RELEASE = 'Sci Rep 13:377 (2023)'
STATES = {'m1': 'me1', 'm2': 'me2', 'm3': 'me3'}
NEEDED = ['Kme_Site', 'Site.Probability', 'Peptide.Sequence', 'Position.in.Peptide', 'Modification', 'Protein.Accession',
          'Position.in.Protein', 'Cell.Type', 'Enrichment']
COLUMNS = ['canonical_UniProtAC', 'position', 'residue', 'mark', 'state', 'evidence', 'evidence_class', 'pmids', 'source',
           'source_release', 'berryhill_site_probability', 'berryhill_cell_lines', 'berryhill_enrichment', 'peptides', 'relocated', 'source_position']
REJECTED_COLUMNS = ['Kme_Site', 'Protein.Accession', 'Position.in.Protein', 'Modification', 'Peptide.Sequence',
                    'Position.in.Peptide', 'Site.Probability', 'Cell.Type', 'Enrichment', 'reason', 'detail']
OUT_SITES = 'berryhill2023_methyllysine_sites.tsv'
OUT_REJECTED = 'berryhill2023_methyllysine_rejected.tsv'
OUT_SUMMARY = 'berryhill2023_methyllysine_summary.json'
CONTENT = ('Berryhill et al. 2023 (Sci Rep 13:377, doi:10.1038/s41598-022-27175-x) human lysine methylation sites '
           '(Supplementary Table 5; HEK293T and U2OS, antibody-enriched and global proteomics), checked on the reference '
           'proteome, one row per residue and state, with the site probability, cell lines, enrichments and peptides')


def to_int(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return int(number) if number == int(number) else None


def to_prob(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def classify(row: dict, reference: dict[str, str]) -> tuple[str, int | None, str]:
    """(status, canonical position, detail): status kept, relocated, or a rejection reason."""
    acc = str(row['Protein.Accession']).strip()
    if '-' in acc:
        return 'isoform_accession', None, acc
    if row['Modification'] not in STATES:
        return 'bad_state', None, str(row['Modification'])
    seq = reference.get(acc)
    if seq is None:
        return 'accession_not_in_reference', None, acc
    pos, pip = to_int(row['Position.in.Protein']), to_int(row['Position.in.Peptide'])
    peptide = str(row['Peptide.Sequence']).strip().upper()
    if pos is None or pos < 1 or pip is None or not 1 <= pip <= len(peptide):
        return 'bad_position', None, f"{row['Position.in.Protein']}/{row['Position.in.Peptide']}"
    if peptide[pip - 1] != 'K':
        return 'peptide_position_not_K', None, f'{peptide[pip - 1]} at {pip} of {peptide}'
    start = pos - pip
    if pos <= len(seq) and seq[pos - 1] == 'K' and start >= 0 and seq[start:start + len(peptide)] == peptide:
        return 'kept', pos, ''
    starts, i = [], seq.find(peptide)
    while i != -1:
        starts.append(i)
        i = seq.find(peptide, i + 1)
    if len(starts) == 1:
        return 'relocated', starts[0] + pip, f'{pos} -> {starts[0] + pip}'
    if len(starts) > 1:
        return 'peptide_not_unique', None, f'{len(starts)} matches of {peptide}'
    there = seq[pos - 1] if pos <= len(seq) else 'past end'
    return 'peptide_not_in_canonical', None, f'{peptide}; canonical has {there} at {pos}'


def stage(rows: list[dict], reference: dict[str, str]) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """The validated sites, the rejected rows and the counts, from the rows of Sup Table 5 (dicts keyed by header)."""
    kept, rejected, status_counts = [], [], Counter()
    for row in rows:
        status, pos, detail = classify(row, reference)
        status_counts[status] += 1
        if pos is None:
            rejected.append({**{c: row[c] for c in REJECTED_COLUMNS[:-2]}, 'reason': status, 'detail': detail})
            continue
        prob = to_prob(row['Site.Probability'])
        kept.append({'canonical_UniProtAC': str(row['Protein.Accession']).strip(), 'position': pos,
                     'state': STATES[row['Modification']], 'prob': prob, 'cell': str(row['Cell.Type']),
                     'enrichment': str(row['Enrichment']), 'peptide': str(row['Peptide.Sequence']).strip().upper(),
                     'relocated': status == 'relocated', 'source_position': to_int(row['Position.in.Protein'])})
    table = pd.DataFrame(kept, columns=['canonical_UniProtAC', 'position', 'state', 'prob', 'cell', 'enrichment', 'peptide',
                                        'relocated', 'source_position'])

    def joined(values):
        return ';'.join(sorted({str(v) for v in values if str(v) != ''}))

    sites = (table.groupby(['canonical_UniProtAC', 'position', 'state'], as_index=False)
             .agg(site_probability_max=('prob', 'max'), cell_lines=('cell', joined), enrichments=('enrichment', joined),
                  peptides=('peptide', joined), relocated=('relocated', 'any'),
                  source_position=('source_position', lambda v: joined(v))))
    sites['residue'], sites['mark'] = 'K', 'methyllysine'
    sites['evidence'], sites['evidence_class'], sites['pmids'] = '', 'experimental', ''
    sites['source'], sites['source_release'] = SOURCE, SOURCE_RELEASE
    sites['site_probability_max'] = [f'{p:g}' if pd.notna(p) else '' for p in sites['site_probability_max']]
    sites['source_position'] = [sp if rel else '' for sp, rel in zip(sites['source_position'], sites['relocated'])]
    sites['relocated'] = ['yes' if r else 'no' for r in sites['relocated']]
    sites = sites.rename(columns={'site_probability_max': 'berryhill_site_probability', 'cell_lines': 'berryhill_cell_lines',
                                  'enrichments': 'berryhill_enrichment'})
    sites = sites[COLUMNS].sort_values(['canonical_UniProtAC', 'position', 'state']).reset_index(drop=True)
    bad = pd.DataFrame(rejected, columns=REJECTED_COLUMNS)
    if len(bad):
        bad = bad.sort_values(['reason', 'Protein.Accession', 'Position.in.Protein', 'Modification', 'Cell.Type', 'Enrichment'],
                              key=lambda c: c.astype(str)).reset_index(drop=True)
    return sites, bad, dict(status_counts)


def read_sheet(path: Path) -> list[dict]:
    sheet = openpyxl.load_workbook(path, read_only=True)[SHEET]
    it = sheet.iter_rows(values_only=True)
    head = list(next(it))
    missing = [c for c in NEEDED if c not in head]
    if missing:
        raise SystemExit(f'{path} {SHEET}: columns missing: {missing}')
    return [dict(zip(head, r)) for r in it if r[0] is not None]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--xlsx', required=True, help='Supplementary Information 2 (41598_2022_27175_MOESM2_ESM.xlsx)')
    ap.add_argument('--fasta', required=True, help='the reference proteome (canonical reviewed accessions)')
    ap.add_argument('--out-dir', default='inputs/marks')
    ap.add_argument('--manifest', default='inputs/MANIFEST.tsv')
    ap.add_argument('--raw-name', default='berryhill2023_41598_2022_27175_MOESM2_ESM.xlsx', help='name of the download in raw_files.tsv')
    ap.add_argument('--downloaded', default='', help='date of the download for raw_files.tsv (default: today)')
    args = ap.parse_args()

    xlsx, out_dir = Path(args.xlsx), Path(args.out_dir)
    raw_bytes = xlsx.read_bytes()
    raw_md5, raw_sha = hashlib.md5(raw_bytes).hexdigest(), hashlib.sha256(raw_bytes).hexdigest()
    rows = read_sheet(xlsx)
    reference = parse_fasta(args.fasta)
    sites, rejected, status_counts = stage(rows, reference)

    out_dir.mkdir(parents=True, exist_ok=True)
    sites.to_csv(out_dir / OUT_SITES, sep='\t', index=False, lineterminator='\n')
    rejected.to_csv(out_dir / OUT_REJECTED, sep='\t', index=False, lineterminator='\n')
    probs = pd.to_numeric(pd.Series([r['Site.Probability'] for r in rows]), errors='coerce')
    summary = {
        'source': f'{CITATION}; {SHEET} of Supplementary Information 2',
        'article_url': ARTICLE_URL,
        'url': URL,
        'licence': LICENSE,
        'reference_fasta': args.fasta,
        'reference_fasta_md5': md5_file(args.fasta),
        'raw_file': args.raw_name,
        'raw_md5': raw_md5,
        'raw_sha256': raw_sha,
        'input_rows': len(rows),
        'input_rows_by_state': dict(sorted(Counter(r['Modification'] for r in rows).items())),
        'row_status': dict(sorted(status_counts.items())),
        'rejected_by_reason': {k: v for k, v in sorted(status_counts.items()) if k not in ('kept', 'relocated')},
        'site_probability_note': ('the paper filters on "high confidence identification and modification site probability" '
                                  'without a number in the text; the column is kept, not filtered. Site.Probability is NA '
                                  f'in {int(probs.isna().sum())} input rows (all PTMScan Rep 2 rows); the others range '
                                  f'{probs.min():g}-{probs.max():g}'),
        'validated_sites': len(sites),
        'validated_proteins': int(sites['canonical_UniProtAC'].nunique()),
        'validated_residues': int(sites[['canonical_UniProtAC', 'position']].drop_duplicates().shape[0]),
        'validated_sites_by_state': {k: int(v) for k, v in sites['state'].value_counts().sort_index().items()},
        'validated_proteins_by_state': {k: int(v) for k, v in sites.groupby('state')['canonical_UniProtAC'].nunique().sort_index().items()},
        'relocated_sites': int((sites['relocated'] == 'yes').sum()),
        'sites_with_probability': int((sites['berryhill_site_probability'] != '').sum()),
        'sites_by_cell_lines': {k: int(v) for k, v in sites['berryhill_cell_lines'].value_counts().sort_index().items()},
        'files': {'validated': OUT_SITES, 'rejected': OUT_REJECTED},
        'files_md5': {'validated': md5_file(out_dir / OUT_SITES), 'rejected': md5_file(out_dir / OUT_REJECTED)},
    }
    (out_dir / OUT_SUMMARY).write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n', encoding='utf-8')

    raw_list = out_dir / 'raw_files.tsv'
    when = args.downloaded or date.today().isoformat()
    upsert_tsv_rows(raw_list, [{
        'file': args.raw_name, 'url': URL, 'bytes': xlsx.stat().st_size, 'md5': raw_md5, 'sha256': raw_sha, 'downloaded': when,
        'version': f'{CITATION}; Supplementary Information 2 (Sup Tables 3-7: antibody tables, Kme sites); only Sup Table 5 is used',
        'kept': f"not in Git; Sup Table 5 is staged as marks/{OUT_SITES} ({LICENSE})"}])

    command = (f'python scripts/57_stage_berryhill_lys_methylation.py --xlsx PTM_data/marks/{args.raw_name} '
               f'--fasta {args.fasta} --out-dir inputs/marks --downloaded {when}')
    check = ('; K and peptide checked on the reference sequence, unique peptide match relocates, others rejected with a reason; '
             f'counts in marks/{OUT_SUMMARY}')
    rows_out = [
        {'file': f'marks/{OUT_SITES}', 'bytes': (out_dir / OUT_SITES).stat().st_size, 'md5': md5_file(out_dir / OUT_SITES),
         'rows': len(sites), 'redistribution': 'git', 'content': CONTENT + f': {len(sites):,} rows', 'source': URL,
         'version_or_date': f'{SOURCE_RELEASE}; downloaded {when}', 'citation': CITATION, 'licence': LICENSE,
         'derived_from': f'{args.raw_name} (MD5 {raw_md5}), sheet {SHEET}', 'how_made_or_checked': command + check},
        {'file': f'marks/{OUT_REJECTED}', 'bytes': (out_dir / OUT_REJECTED).stat().st_size, 'md5': md5_file(out_dir / OUT_REJECTED),
         'rows': len(rejected), 'redistribution': 'git',
         'content': f'Rows of Berryhill et al. 2023 Supplementary Table 5 that failed the reference check, with the reason: {len(rejected):,} rows',
         'source': URL, 'version_or_date': f'{SOURCE_RELEASE}; downloaded {when}', 'citation': CITATION, 'licence': LICENSE,
         'derived_from': f'{args.raw_name} (MD5 {raw_md5}), sheet {SHEET}', 'how_made_or_checked': command + check},
        {'file': f'marks/{OUT_SUMMARY}', 'bytes': (out_dir / OUT_SUMMARY).stat().st_size, 'md5': md5_file(out_dir / OUT_SUMMARY),
         'redistribution': 'git', 'content': 'Counts, drops by reason, site-probability note and the MD5/SHA-256 of the input and MD5 of both tables of the Berryhill 2023 staging',
         'source': 'the tables above', 'version_or_date': SOURCE_RELEASE, 'citation': CITATION, 'licence': LICENSE,
         'derived_from': 'marks/raw_files.tsv', 'how_made_or_checked': command},
    ]
    update_manifest(args.manifest, rows_out)
    refresh_manifest_row(args.manifest, raw_list, 'marks/raw_files.tsv', redistribution='git',
                         content='URL, size, MD5, SHA-256 and date of every raw download behind inputs/marks/',
                         source='the URLs it lists', derived_from='',
                         how_made_or_checked='written by the scripts/54, 55, 57, 58 and 59 staging runs')
    print(json.dumps({k: summary[k] for k in ('input_rows', 'row_status', 'validated_sites', 'validated_proteins',
                                               'validated_sites_by_state', 'validated_proteins_by_state', 'relocated_sites')}, indent=1))


if __name__ == '__main__':
    main()

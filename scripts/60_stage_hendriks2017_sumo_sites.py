#!/usr/bin/env python3
"""Stage the human SUMO-2 sites of Hendriks et al. 2017 as a sumoylation table (inputs/marks/).

    python scripts/60_stage_hendriks2017_sumo_sites.py --xlsx PTM_data/marks/41594_2017_BFnsmb3366_MOESM30_ESM.xlsx \
        --fasta inputs/reference/uniprot_human_reviewed_2026_03.fasta \
        --history inputs/reference/uniprot_accession_history_2026_03.tsv --out-dir inputs/marks --downloaded 2026-11-13

Source: Hendriks IA, Lyon D, Young C, Jensen LJ, Vertegaal ACO, Nielsen ML. Site-specific mapping of the human SUMO
proteome reveals co-modification with phosphorylation. Nat Struct Mol Biol 24:325-336 (2017), doi:10.1038/nsmb.3366
(PubMed 28112733; not in PMC; PRIDE PXD005296). The article is not open access: its page gives "2017 Springer Nature
America, Inc." as the copyright, and the supplementary file carries no license or notice of its own. The derived table
(accession, position, localization probability, cell lines) is published in the release, cited. Cite the paper wherever
the data is described. The download stays out of Git (raw_files.tsv has its MD5 and SHA-256).

Reads only Supplementary Table 1 (MOESM30_ESM.xlsx, sheet Hendriks_TableS1; a two-row header, the second row naming
the columns): "A list of all 40,765 identified SUMOylation sites", one row per SUMO-2 acceptor lysine, as published
(the paper's own localization and FDR filters already applied; nothing is filtered here). Columns used: Protein and
Position (the leading protein), Proteins and Positions within proteins (the protein group), Sequence window (51 AA,
Unique), Localization prob and the eight "(MS/MS detection)" columns (HeLa and U2OS, each control, BTZ, Heat and MG132).
The other columns (quantification, structural context, peptides) and Supplementary Tables 2-11 are not used.

Each row is placed on the reference proteome (--fasta, canonical reviewed accessions; --history, its secondary
accessions) by scripts/window_sites.py: kept where a protein of the group is a canonical entry with K at the given
position and the window agreeing (+/-7 residues), secondary_accession where that needs the current primary accession,
relocated where the window's core occurs once in the canonical sequence of a protein of the group (mostly sites whose
leading protein is an isoform); otherwise rejected with a reason. Rows that land on the same canonical lysine are one
row. Columns: the ten of the other mark tables (state sumo, evidence empty, evidence_class experimental (mass
spectrometry), pmids empty, source 'Hendriks 2017', source_release 'Nat Struct Mol Biol 24:325 (2017)') and
  hendriks_localization_probability  the highest Localization prob of the rows of the site (MaxQuant, 0-1)
  hendriks_cell_lines                HeLa and/or U2OS: the cell lines with an MS/MS detection of the site
  source_accession, source_position  the leading protein and position of the row(s), ';'-joined
  placement                          kept, secondary_accession or relocated (';'-joined if the rows differ)
The release (scripts/56) takes the first two. The staged table is named `hendriks2017_sumoylation_sites.tsv`, which
scripts/56's glob `inputs/marks/*_sumoylation_sites.tsv` picks up. Outputs are sorted, so running again gives the same bytes.
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
from window_sites import REJECTIONS, core, group, place, read_secondary

ARTICLE_URL = 'https://www.nature.com/articles/nsmb.3366'
URL = 'https://static-content.springer.com/esm/art%3A10.1038%2Fnsmb.3366/MediaObjects/41594_2017_BFnsmb3366_MOESM30_ESM.xlsx'
RAW = '41594_2017_BFnsmb3366_MOESM30_ESM.xlsx'
SHEET = 'Hendriks_TableS1'
CITATION = 'Hendriks et al., Nat Struct Mol Biol 24:325-336 (2017), doi:10.1038/nsmb.3366'
LICENSE = ('article page: copyright "2017 Springer Nature America, Inc."; not open access; '
           'the supplementary file carries no separate notice')
SOURCE = 'Hendriks 2017'
SOURCE_RELEASE = 'Nat Struct Mol Biol 24:325 (2017)'
WINDOW, CENTER = 'Sequence window (51 AA, Unique)', 25
DETECTION = {'HeLa': ['HeLa (MS/MS detection)', 'HeLa BTZ (MS/MS detection)', 'HeLa Heat (MS/MS detection)', 'HeLa MG132 (MS/MS detection)'],
             'U2OS': ['U2OS (MS/MS detection)', 'U2OS BTZ (MS/MS detection)', 'U2OS Heat (MS/MS detection)', 'U2OS MG132 (MS/MS detection)']}
NEEDED = ['Protein', 'Position', 'Proteins', 'Positions within proteins', WINDOW, 'Localization prob', 'Amino acid',
          *[c for cols in DETECTION.values() for c in cols]]
COLUMNS = ['canonical_UniProtAC', 'position', 'residue', 'mark', 'state', 'evidence', 'evidence_class', 'pmids', 'source',
           'source_release', 'hendriks_localization_probability', 'hendriks_cell_lines', 'source_accession', 'source_position', 'placement']
REJECTED_COLUMNS = ['Protein', 'Position', 'Proteins', 'Localization prob', 'window', 'reason', 'detail']
OUT_SITES = 'hendriks2017_sumoylation_sites.tsv'
OUT_REJECTED = 'hendriks2017_sumoylation_rejected.tsv'
OUT_SUMMARY = 'hendriks2017_sumoylation_summary.json'
CONTENT = ('Hendriks et al. 2017 (Nat Struct Mol Biol 24:325, doi:10.1038/nsmb.3366) human SUMO-2 acceptor lysines '
           '(Supplementary Table 1; HeLa and U2OS), placed on the reference proteome, one row per lysine, with the '
           'localization probability and cell lines')


def joined(values) -> str:
    return ';'.join(sorted({str(v) for v in values if str(v) != ''}))


def stage(rows: list[dict], reference: dict[str, str], secondary: dict[str, str]) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """The placed sites, the rejected rows and the counts, from the rows of Supplementary Table 1 (dicts keyed by header)."""
    kept, rejected, counts = [], [], Counter()
    for row in rows:
        motif = core(str(row[WINDOW] or ''), CENTER)
        pairs = group(row['Proteins'], row['Positions within proteins'], row['Protein'], row['Position'])
        status, acc, pos, detail = place(pairs, [motif], reference, secondary)
        counts[status] += 1
        if acc is None:
            rejected.append({'Protein': row['Protein'], 'Position': row['Position'], 'Proteins': row['Proteins'],
                             'Localization prob': row['Localization prob'], 'window': row[WINDOW], 'reason': status, 'detail': detail})
            continue
        cells = [c for c, cols in DETECTION.items() if any(str(row[x] or '').strip() == '+' for x in cols)]
        kept.append({'canonical_UniProtAC': acc, 'position': pos, 'prob': float(row['Localization prob']),
                     'cells': ';'.join(cells), 'source_accession': str(row['Protein']).strip(),
                     'source_position': str(row['Position']).strip(), 'placement': status})
    table = pd.DataFrame(kept, columns=['canonical_UniProtAC', 'position', 'prob', 'cells', 'source_accession', 'source_position', 'placement'])
    sites = (table.groupby(['canonical_UniProtAC', 'position'], as_index=False)
             .agg(prob=('prob', 'max'), cells=('cells', lambda v: joined(c for x in v for c in x.split(';'))),
                  source_accession=('source_accession', joined), source_position=('source_position', joined),
                  placement=('placement', joined)))
    sites['residue'], sites['mark'], sites['state'] = 'K', 'sumoylation', 'sumo'
    sites['evidence'], sites['evidence_class'], sites['pmids'] = '', 'experimental', ''
    sites['source'], sites['source_release'] = SOURCE, SOURCE_RELEASE
    sites['hendriks_localization_probability'] = [f'{p:g}' for p in sites['prob']]
    sites['hendriks_cell_lines'] = sites['cells']
    sites = sites[COLUMNS].sort_values(['canonical_UniProtAC', 'position']).reset_index(drop=True)
    bad = pd.DataFrame(rejected, columns=REJECTED_COLUMNS)
    if len(bad):
        bad = bad.sort_values(['reason', 'Protein', 'Position'], key=lambda c: c.astype(str)).reset_index(drop=True)
    return sites, bad, dict(counts)


def read_sheet(path: Path) -> list[dict]:
    it = openpyxl.load_workbook(path, read_only=True, data_only=True)[SHEET].iter_rows(values_only=True)
    next(it)                                  # the first header row names column groups
    head = list(next(it))
    missing = [c for c in NEEDED if c not in head]
    if missing:
        raise SystemExit(f'{path} {SHEET}: columns missing: {missing}')
    keep = [head.index(c) for c in NEEDED]
    return [dict(zip(NEEDED, (r[i] for i in keep))) for r in it if r[0] is not None]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--xlsx', required=True, help=f'Supplementary Table 1 ({RAW})')
    ap.add_argument('--fasta', required=True, help='the reference proteome (canonical reviewed accessions)')
    ap.add_argument('--history', required=True, help='secondary accessions of the reference (uniprot_accession_history_<release>.tsv)')
    ap.add_argument('--out-dir', default='inputs/marks')
    ap.add_argument('--manifest', default='inputs/MANIFEST.tsv')
    ap.add_argument('--downloaded', default='', help='date of the download for raw_files.tsv (default: today)')
    args = ap.parse_args()

    xlsx, out_dir = Path(args.xlsx), Path(args.out_dir)
    data = xlsx.read_bytes()
    raw_md5, raw_sha = hashlib.md5(data).hexdigest(), hashlib.sha256(data).hexdigest()
    rows = read_sheet(xlsx)
    if {r['Amino acid'] for r in rows} != {'K'}:
        raise SystemExit(f'{SHEET}: rows other than lysine')
    reference = parse_fasta(args.fasta)
    secondary = read_secondary(args.history, reference)
    sites, rejected, counts = stage(rows, reference, secondary)

    out_dir.mkdir(parents=True, exist_ok=True)
    sites.to_csv(out_dir / OUT_SITES, sep='\t', index=False, lineterminator='\n')
    rejected.to_csv(out_dir / OUT_REJECTED, sep='\t', index=False, lineterminator='\n')
    probs = pd.Series([float(r['Localization prob']) for r in rows])
    placed = sum(v for k, v in counts.items() if k not in REJECTIONS)
    summary = {
        'source': f'{CITATION}; Supplementary Table 1 (sheet {SHEET})',
        'article_url': ARTICLE_URL,
        'url': URL,
        'license': LICENSE,
        'redistribution': 'git',
        'reference_fasta': args.fasta,
        'reference_fasta_md5': md5_file(args.fasta),
        'accession_history': args.history,
        'accession_history_md5': md5_file(args.history),
        'raw_file': RAW,
        'raw_md5': raw_md5,
        'raw_sha256': raw_sha,
        'input_rows': len(rows),
        'input_leading_proteins': len({r['Protein'] for r in rows}),
        'input_rows_leading_isoform': sum('-' in str(r['Protein']) for r in rows),
        'row_status': dict(sorted(counts.items())),
        'rows_placed': placed,
        'rejected_by_reason': {k: v for k, v in sorted(counts.items()) if k in REJECTIONS},
        'rows_merged_into_another_site': placed - len(sites),
        'localization_probability_note': ('the table is used as published (the paper filters its sites itself); Localization prob '
                                          f'ranges {probs.min():g}-{probs.max():g} over the input rows, {int((probs < 0.75).sum())} below 0.75'),
        'validated_sites': len(sites),
        'validated_proteins': int(sites['canonical_UniProtAC'].nunique()),
        'sites_by_placement': {k: int(v) for k, v in sites['placement'].value_counts().sort_index().items()},
        'sites_by_cell_lines': {k: int(v) for k, v in sites['hendriks_cell_lines'].value_counts().sort_index().items()},
        'files': {'validated': OUT_SITES, 'rejected': OUT_REJECTED},
        'files_md5': {'validated': md5_file(out_dir / OUT_SITES), 'rejected': md5_file(out_dir / OUT_REJECTED)},
    }
    (out_dir / OUT_SUMMARY).write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n', encoding='utf-8')

    when = args.downloaded or date.today().isoformat()
    raw_list = out_dir / 'raw_files.tsv'
    upsert_tsv_rows(raw_list, [{
        'file': RAW, 'url': URL, 'bytes': len(data), 'md5': raw_md5, 'sha256': raw_sha, 'downloaded': when,
        'version': f'{CITATION}; Supplementary Table 1 (all 40,765 SUMOylation sites)',
        'kept': f'not in Git; staged as marks/{OUT_SITES} ({LICENSE})'}])

    command = (f'python scripts/60_stage_hendriks2017_sumo_sites.py --xlsx PTM_data/marks/{RAW} --fasta {args.fasta} '
               f'--history {args.history} --out-dir inputs/marks --downloaded {when}')
    check = ('; K and +/-7 window checked on the reference sequence, isoform and moved sites placed by a unique window match, '
             f'others rejected with a reason; counts in marks/{OUT_SUMMARY}')
    common_fields = {'source': URL, 'version_or_date': f'{SOURCE_RELEASE}; downloaded {when}', 'citation': CITATION,
                     'licence': LICENSE, 'derived_from': f'{RAW} (MD5 {raw_md5}), sheet {SHEET}', 'how_made_or_checked': command + check}
    update_manifest(args.manifest, [
        {'file': f'marks/{OUT_SITES}', 'bytes': (out_dir / OUT_SITES).stat().st_size, 'md5': md5_file(out_dir / OUT_SITES),
         'rows': len(sites), 'redistribution': 'git', 'content': CONTENT + f': {len(sites):,} rows', **common_fields},
        {'file': f'marks/{OUT_REJECTED}', 'bytes': (out_dir / OUT_REJECTED).stat().st_size, 'md5': md5_file(out_dir / OUT_REJECTED),
         'rows': len(rejected), 'redistribution': 'git',
         'content': f'Rows of Hendriks et al. 2017 Supplementary Table 1 that could not be placed on the reference, with the reason: {len(rejected):,} rows',
         **common_fields},
        {'file': f'marks/{OUT_SUMMARY}', 'bytes': (out_dir / OUT_SUMMARY).stat().st_size, 'md5': md5_file(out_dir / OUT_SUMMARY),
         'redistribution': 'git', 'content': 'Counts, drops by reason and the MD5/SHA-256 of the input and MD5 of both tables of the Hendriks 2017 staging',
         **{**common_fields, 'source': 'the tables above', 'derived_from': 'marks/raw_files.tsv', 'how_made_or_checked': command}},
    ])
    refresh_manifest_row(args.manifest, raw_list, 'marks/raw_files.tsv', redistribution='git',
                         content='URL, size, MD5, SHA-256 and date of every raw download behind inputs/marks/',
                         source='the URLs it lists', derived_from='',
                         how_made_or_checked='written by the scripts/54, 55, 57, 58, 59, 60 and 61 staging runs')
    print(json.dumps({k: summary[k] for k in ('input_rows', 'row_status', 'rows_placed', 'rows_merged_into_another_site',
                                               'validated_sites', 'validated_proteins', 'sites_by_placement')}, indent=1))


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Stage the human lysine ubiquitination sites of Akimov et al. 2018 (UbiSite) as a ubiquitination table (inputs/marks/).

    python scripts/61_stage_akimov2018_ubiquitination_sites.py --xlsx PTM_data/marks/41594_2018_84_MOESM3_ESM.xlsx \
        --fasta inputs/reference/uniprot_human_reviewed_2026_03.fasta \
        --history inputs/reference/uniprot_accession_history_2026_03.tsv --out-dir inputs/marks --downloaded 2026-11-13

Source: Akimov V, Barrio-Hernandez I, Hansen SVF, Hallenborg P, Pedersen AK, Bekker-Jensen DB, Puglia M,
Christensen SDK, Vanselow JT, Nielsen MM, Kratchmarova I, Kelstrup CD, Olsen JV, Blagoev B. UbiSite approach for
comprehensive mapping of lysine and N-terminal ubiquitination sites. Nat Struct Mol Biol 25:631-640 (2018),
doi:10.1038/s41594-018-0084-y (PubMed 29967540; not in PMC). The article is not open access: its page gives
"2018 The Author(s)" as the copyright, and the supplementary file carries no license or notice of its own. The derived
table (accession, position, localization probability, cell lines) is published in the release, cited. Cite the paper
wherever the data is described. The download stays out of Git (raw_files.tsv has its MD5 and SHA-256).

Reads Supplementary Table 1 (MOESM3_ESM.xlsx, "Lysine ubiquitination sites of the UbiSite-enriched ubiquitinomes of
Hep2 and Jurkat cells"). Its sheet "Combined Ubi sites" is the site list (one row per lysine, the sites of the
LysC+trypsin and the LysC-only digests combined, as published; nothing is filtered here): id.Tryp/id.LysC,
Protein.*, Position.*, Localization.prob.*, Sequence.window.* of each digest, and HEP2* and JURK* (1 where the site has
an intensity value in that cell line). The protein group of each site (Proteins, Positions.within.proteins) is read
from the digest's own sheet ("LysC+Trypsin Ubi sites", "LysC only Ubi sites") by its id. The sheet "Comparaison with
Udeshi 2013", the other columns, and the ISG15 and NEDD8 sites listed below the LysC-only sites (after an empty row,
under headers of their own) are not used. N-terminal ubiquitination sites (Supplementary Table 6, peptides of
N-terminally ubiquitinated proteins) are left out: the table of this release is of lysine sites only, and that file is
not downloaded.

Each site is placed on the reference proteome (--fasta, canonical reviewed accessions; --history, its secondary
accessions) by scripts/window_sites.py (the window is 31 residues, centre 16; kept, secondary_accession or relocated,
otherwise rejected with a reason), from its trypsin entry first and from its LysC entry when that fails. Sites that
land on the same canonical lysine are one row. Columns: the ten of the other mark tables (state ub, evidence empty,
evidence_class experimental (mass spectrometry), pmids empty, source 'Akimov 2018', source_release
'Nat Struct Mol Biol 25:631 (2018)') and
  akimov_localization_probability  the highest localization probability of the site over both digests (MaxQuant, 0-1)
  akimov_cell_lines                Hep2 and/or Jurkat: the cell lines in which the site has an intensity value
                                   (HEP2*, JURK*); empty where it has none
  source_accession, source_position  the leading protein and position, ';'-joined
  placement                        kept, secondary_accession or relocated
The release (scripts/56) takes the first two. The staged table is named `akimov2018_ubiquitination_sites.tsv`, which
scripts/56's glob `inputs/marks/*_ubiquitination_sites.tsv` picks up. Outputs are sorted, so running again gives the same bytes.
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
from window_sites import REJECTIONS, core, group, place, read_secondary, to_int

ARTICLE_URL = 'https://www.nature.com/articles/s41594-018-0084-y'
URL = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41594-018-0084-y/MediaObjects/41594_2018_84_MOESM3_ESM.xlsx'
RAW = '41594_2018_84_MOESM3_ESM.xlsx'
SHEET = 'Combined Ubi sites'
DIGESTS = {'Tryp': ('LysC+Trypsin Ubi sites', 2), 'LysC': ('LysC only Ubi sites', 0)}     # sheet, rows above the header
CITATION = 'Akimov et al., Nat Struct Mol Biol 25:631-640 (2018), doi:10.1038/s41594-018-0084-y'
LICENSE = ('article page: copyright "2018 The Author(s)"; not open access; '
           'the supplementary file carries no separate notice')
SOURCE = 'Akimov 2018'
SOURCE_RELEASE = 'Nat Struct Mol Biol 25:631 (2018)'
CENTER = 15
CELLS = {'Hep2': 'HEP2*', 'Jurkat': 'JURK*'}
NEEDED = [f'{c}.{d}' for d in DIGESTS for c in ('id', 'Protein', 'Position', 'Localization.prob', 'Sequence.window')] + list(CELLS.values())
NEEDED_DIGEST = ['id', 'Proteins', 'Positions.within.proteins', 'Amino.acid']
COLUMNS = ['canonical_UniProtAC', 'position', 'residue', 'mark', 'state', 'evidence', 'evidence_class', 'pmids', 'source',
           'source_release', 'akimov_localization_probability', 'akimov_cell_lines', 'source_accession', 'source_position', 'placement']
REJECTED_COLUMNS = ['id.Tryp', 'id.LysC', 'protein', 'position', 'localization_probability', 'window', 'reason', 'detail']
OUT_SITES = 'akimov2018_ubiquitination_sites.tsv'
OUT_REJECTED = 'akimov2018_ubiquitination_rejected.tsv'
OUT_SUMMARY = 'akimov2018_ubiquitination_summary.json'
CONTENT = ('Akimov et al. 2018 (Nat Struct Mol Biol 25:631, doi:10.1038/s41594-018-0084-y) human lysine ubiquitination '
           'sites (UbiSite; Supplementary Table 1, combined sites; Hep2 and Jurkat), placed on the reference proteome, one row '
           'per lysine, with the localization probability and cell lines; N-terminal sites not included')


def joined(values) -> str:
    return ';'.join(sorted({str(v) for v in values if str(v) != ''}))


def prob(value) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def stage(rows: list[dict], groups: dict[str, dict[int, dict]], reference: dict[str, str],
          secondary: dict[str, str]) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """The placed sites, the rejected rows and the counts, from the rows of the combined sheet and the protein groups of
    each digest's sheet (digest -> id -> row)."""
    kept, rejected, counts = [], [], Counter()
    for row in rows:
        first = None
        for d in DIGESTS:
            sid = to_int(row[f'id.{d}'])
            if sid is None:
                continue
            g = groups[d].get(sid, {})
            pairs = group(g.get('Proteins'), g.get('Positions.within.proteins'), row[f'Protein.{d}'], row[f'Position.{d}'])
            motifs = [core(w, CENTER) for w in str(row[f'Sequence.window.{d}'] or '').split(';')]
            result = place(pairs, motifs, reference, secondary)
            first = first or (d, result)
            if result[1] is not None:
                first = (d, result)
                break
        d, (status, acc, pos, detail) = first
        counts[status] += 1
        probs = [p for p in (prob(row['Localization.prob.Tryp']), prob(row['Localization.prob.LysC'])) if p is not None]
        if acc is None:
            ids = {k: '' if to_int(row[k]) is None else str(to_int(row[k])) for k in ('id.Tryp', 'id.LysC')}
            rejected.append({**ids, 'protein': row[f'Protein.{d}'],
                             'position': row[f'Position.{d}'], 'localization_probability': f'{max(probs):g}' if probs else '',
                             'window': row[f'Sequence.window.{d}'], 'reason': status, 'detail': detail})
            continue
        cells = [c for c, col in CELLS.items() if to_int(row[col]) == 1]
        kept.append({'canonical_UniProtAC': acc, 'position': pos, 'prob': max(probs) if probs else None, 'cells': ';'.join(cells),
                     'source_accession': str(row[f'Protein.{d}']).strip(), 'source_position': str(to_int(row[f'Position.{d}'])),
                     'placement': status})
    table = pd.DataFrame(kept, columns=['canonical_UniProtAC', 'position', 'prob', 'cells', 'source_accession', 'source_position', 'placement'])
    sites = (table.groupby(['canonical_UniProtAC', 'position'], as_index=False)
             .agg(prob=('prob', 'max'), cells=('cells', lambda v: joined(c for x in v for c in x.split(';'))),
                  source_accession=('source_accession', joined), source_position=('source_position', joined),
                  placement=('placement', joined)))
    sites['residue'], sites['mark'], sites['state'] = 'K', 'ubiquitination', 'ub'
    sites['evidence'], sites['evidence_class'], sites['pmids'] = '', 'experimental', ''
    sites['source'], sites['source_release'] = SOURCE, SOURCE_RELEASE
    sites['akimov_localization_probability'] = [f'{p:g}' if pd.notna(p) else '' for p in sites['prob']]
    sites['akimov_cell_lines'] = sites['cells']
    sites = sites[COLUMNS].sort_values(['canonical_UniProtAC', 'position']).reset_index(drop=True)
    bad = pd.DataFrame(rejected, columns=REJECTED_COLUMNS)
    if len(bad):
        bad = bad.sort_values(['reason', 'protein', 'position'], key=lambda c: c.astype(str)).reset_index(drop=True)
    return sites, bad, dict(counts)


def read_rows(book, sheet: str, skip: int, needed: list[str], first_block: bool = False) -> list[dict]:
    """The rows below the header, as dicts of the needed columns; with first_block, only the rows up to the first empty row."""
    it = book[sheet].iter_rows(values_only=True)
    for _ in range(skip):
        next(it)
    head = list(next(it))
    missing = [c for c in needed if c not in head]
    if missing:
        raise SystemExit(f'{sheet}: columns missing: {missing}')
    keep = [head.index(c) for c in needed]
    out = []
    for r in it:
        if first_block and all(c is None for c in r):
            break
        out.append(dict(zip(needed, (r[i] for i in keep))))
    return out


def read_book(path: Path) -> tuple[list[dict], dict[str, dict[int, dict]], dict]:
    """The site rows of the combined sheet (rows with an id of either digest; the footnote rows have none), the protein
    groups of each digest by id, and counts. A digest sheet is read down to its first empty row: below it, the LysC-only
    sheet lists ISG15 and NEDD8 sites under headers of their own, which are not ubiquitination sites and are not used."""
    book = openpyxl.load_workbook(path, read_only=True, data_only=True)
    rows = [r for r in read_rows(book, SHEET, 0, NEEDED) if to_int(r['id.Tryp']) is not None or to_int(r['id.LysC']) is not None]
    groups, info = {}, {}
    for d, (sheet, skip) in DIGESTS.items():
        digest = read_rows(book, sheet, skip, NEEDED_DIGEST, first_block=True)
        if any(to_int(r['id']) is None for r in digest):
            raise SystemExit(f'{sheet}: a row without an id')
        if {r['Amino.acid'] for r in digest} != {'K'}:
            raise SystemExit(f'{sheet}: rows other than lysine')
        groups[d] = {to_int(r['id']): r for r in digest}
        if len(groups[d]) != len(digest):
            raise SystemExit(f'{sheet}: an id is listed twice')
        info[sheet] = len(digest)
    return rows, groups, info


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
    rows, groups, digest_rows = read_book(xlsx)
    reference = parse_fasta(args.fasta)
    secondary = read_secondary(args.history, reference)
    sites, rejected, counts = stage(rows, groups, reference, secondary)

    out_dir.mkdir(parents=True, exist_ok=True)
    sites.to_csv(out_dir / OUT_SITES, sep='\t', index=False, lineterminator='\n')
    rejected.to_csv(out_dir / OUT_REJECTED, sep='\t', index=False, lineterminator='\n')
    best = pd.Series([max(p for p in (prob(r['Localization.prob.Tryp']), prob(r['Localization.prob.LysC'])) if p is not None) for r in rows])
    placed = sum(v for k, v in counts.items() if k not in REJECTIONS)
    by_digest = Counter('both' if to_int(r['id.Tryp']) is not None and to_int(r['id.LysC']) is not None
                        else 'LysC+trypsin only' if to_int(r['id.Tryp']) is not None else 'LysC only' for r in rows)
    summary = {
        'source': f'{CITATION}; Supplementary Table 1 (sheet {SHEET}; protein groups from the sheets {", ".join(s for s, _ in DIGESTS.values())})',
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
        'input_rows_by_digest': dict(sorted(by_digest.items())),
        'digest_sheet_rows': digest_rows,
        'n_terminal_sites': 'not included: Supplementary Table 6 (peptides of N-terminally ubiquitinated proteins) is not used; lysine sites only',
        'row_status': dict(sorted(counts.items())),
        'rows_placed': placed,
        'rejected_by_reason': {k: v for k, v in sorted(counts.items()) if k in REJECTIONS},
        'rows_merged_into_another_site': placed - len(sites),
        'localization_probability_note': ('the table is used as published (the paper filters its sites itself); the best localization '
                                          f'probability of a site ranges {best.min():g}-{best.max():g}, {int((best < 0.75).sum())} sites below 0.75'),
        'validated_sites': len(sites),
        'validated_proteins': int(sites['canonical_UniProtAC'].nunique()),
        'sites_by_placement': {k: int(v) for k, v in sites['placement'].value_counts().sort_index().items()},
        'sites_by_cell_lines': {(k or 'none'): int(v) for k, v in sites['akimov_cell_lines'].value_counts().sort_index().items()},
        'files': {'validated': OUT_SITES, 'rejected': OUT_REJECTED},
        'files_md5': {'validated': md5_file(out_dir / OUT_SITES), 'rejected': md5_file(out_dir / OUT_REJECTED)},
    }
    (out_dir / OUT_SUMMARY).write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n', encoding='utf-8')

    when = args.downloaded or date.today().isoformat()
    raw_list = out_dir / 'raw_files.tsv'
    upsert_tsv_rows(raw_list, [{
        'file': RAW, 'url': URL, 'bytes': len(data), 'md5': raw_md5, 'sha256': raw_sha, 'downloaded': when,
        'version': f'{CITATION}; Supplementary Table 1 (lysine ubiquitination sites); the Udeshi 2013 comparison sheet is not used',
        'kept': f'not in Git; staged as marks/{OUT_SITES} ({LICENSE})'}])

    command = (f'python scripts/61_stage_akimov2018_ubiquitination_sites.py --xlsx PTM_data/marks/{RAW} --fasta {args.fasta} '
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
         'content': f'Sites of Akimov et al. 2018 Supplementary Table 1 that could not be placed on the reference, with the reason: {len(rejected):,} rows',
         **common_fields},
        {'file': f'marks/{OUT_SUMMARY}', 'bytes': (out_dir / OUT_SUMMARY).stat().st_size, 'md5': md5_file(out_dir / OUT_SUMMARY),
         'redistribution': 'git', 'content': 'Counts, drops by reason and the MD5/SHA-256 of the input and MD5 of both tables of the Akimov 2018 staging',
         **{**common_fields, 'source': 'the tables above', 'derived_from': 'marks/raw_files.tsv', 'how_made_or_checked': command}},
    ])
    refresh_manifest_row(args.manifest, raw_list, 'marks/raw_files.tsv', redistribution='git',
                         content='URL, size, MD5, SHA-256 and date of every raw download behind inputs/marks/',
                         source='the URLs it lists', derived_from='',
                         how_made_or_checked='written by the scripts/54, 55, 57, 58, 59, 60 and 61 staging runs')
    print(json.dumps({k: summary[k] for k in ('input_rows', 'input_rows_by_digest', 'row_status', 'rows_placed', 'rows_merged_into_another_site',
                                               'validated_sites', 'validated_proteins', 'sites_by_placement', 'sites_by_cell_lines')}, indent=1))


if __name__ == '__main__':
    main()

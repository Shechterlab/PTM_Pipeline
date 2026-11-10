#!/usr/bin/env python3
"""Stage the human phosphosites of Ochoa et al. 2020 as a phosphorylation table with functional scores (inputs/marks/).

    python scripts/58_stage_ochoa_phosphosites.py --scores PTM_data/marks/41587_2019_344_MOESM5_ESM.xlsx \
        --features PTM_data/marks/41587_2019_344_MOESM4_ESM.xlsx \
        --fasta inputs/reference/uniprot_human_reviewed_2026_03.fasta --out-dir inputs/marks --downloaded 2026-11-10

Source: Ochoa D, Jarnuczak AF, Viéitez C, Gehre M, Soucheray M, Mateus A, Kleefeldt AA, Hill A, Garcia-Alonso L,
Stein F, Krogan NJ, Savitski MM, Swaney DL, Vizcaíno JA, Noh KM, Beltrao P. The functional landscape of the human
phosphoproteome. Nat Biotechnol 38:365-373 (2020), doi:10.1038/s41587-019-0344-3. Included by the owner's decision
(2026-11-10). The article is NOT open access: its page says "© 2019 The Author(s), under exclusive licence to Springer
Nature America, Inc." and gives the supplements no licence of their own. The staged tables are therefore not
redistributed: MANIFEST redistribution = checksum_only, and inputs/marks/.gitignore lists them. This script rebuilds them
from the two downloads (same bytes). The counts file (summary JSON) is in Git. Cite the paper wherever the data is described.

Reads two sheets, each with one row per site (116,258 sites), keyed by `uniprot` and `position`:
  Supplementary Table 3 (MOESM5_ESM.xlsx, sheet functional_score): functional_score (0-1). It has no residue and no
      evidence columns, so Supplementary Table 2 is read for those.
  Supplementary Table 2 (MOESM4_ESM.xlsx, sheet annotated_phosphoproteome): only residue, best_PEP,
      best_localization_prob, biological_samples and spectralcounts. Its other columns, and its other sheets (one of
      which lists known regulatory sites taken from a third-party database), are not read.
The two sheets must hold the same sites; otherwise the script stops.

Each site is checked on the reference proteome (--fasta, canonical reviewed accessions):
  kept      the accession is a canonical entry of the reference, the table's residue is S, T or Y, and the canonical
            sequence has that residue at the position
  rejected  everything else, with a reason: isoform_accession, accession_not_in_reference, bad_residue (not S/T/Y),
            bad_position, residue_mismatch
Relocation by sequence window is not possible: neither table gives a peptide or window around the site, so a site
whose residue does not match is rejected, never moved.
Columns of the validated table: the ten of the other mark tables (state phospho, evidence empty, evidence_class
high-throughput, pmids empty, source Ochoa2020, source_release 'Nat Biotechnol 38:365 (2020)') and
  functional_score          Supplementary Table 3, 0-1, higher means more likely functional
  best_localization_prob    the best localisation probability over the observations of the site (MaxQuant)
  best_PEP                  the best posterior error probability of a peptide carrying the site
  biological_samples        the number of biological samples (experimental conditions) that observed the site
  spectral_counts           the total spectral counts of the site
The validated table is named `..._validated.tsv`, not `..._sites.tsv`, so that scripts/56's glob
`inputs/marks/*_phosphorylation_sites.tsv` does not pick it up before the release merge step adds it deliberately.
Outputs are sorted, so running again gives the same bytes.
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

ARTICLE_URL = 'https://www.nature.com/articles/s41587-019-0344-3'
BASE = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-019-0344-3/MediaObjects/'
URL_SCORES = BASE + '41587_2019_344_MOESM5_ESM.xlsx'
URL_FEATURES = BASE + '41587_2019_344_MOESM4_ESM.xlsx'
RAW_SCORES, RAW_FEATURES = '41587_2019_344_MOESM5_ESM.xlsx', '41587_2019_344_MOESM4_ESM.xlsx'
SHEET_SCORES, SHEET_FEATURES = 'functional_score', 'annotated_phosphoproteome'
CITATION = 'Ochoa et al., Nat Biotechnol 38:365-373 (2020), doi:10.1038/s41587-019-0344-3'
LICENSE = ('article page: "© 2019 The Author(s), under exclusive licence to Springer Nature America, Inc."; not open access; '
           "supplements carry no separate notice; included by the owner's decision (2026-11-10)")
SOURCE = 'Ochoa2020'
SOURCE_RELEASE = 'Nat Biotechnol 38:365 (2020)'
RESIDUES = 'STY'
NEEDED_SCORES = ['uniprot', 'position', 'functional_score']
NEEDED_FEATURES = ['uniprot', 'position', 'residue', 'best_PEP', 'best_localization_prob', 'biological_samples', 'spectralcounts']
COLUMNS = ['canonical_UniProtAC', 'position', 'residue', 'mark', 'state', 'evidence', 'evidence_class', 'pmids', 'source',
           'source_release', 'functional_score', 'best_localization_prob', 'best_PEP', 'biological_samples', 'spectral_counts']
REJECTED_COLUMNS = ['uniprot', 'position', 'residue', 'functional_score', 'reason', 'detail']
OUT_SITES = 'ochoa2020_phosphorylation_validated.tsv'
OUT_REJECTED = 'ochoa2020_phosphorylation_rejected.tsv'
OUT_SUMMARY = 'ochoa2020_phosphorylation_summary.json'
CONTENT = ('Ochoa et al. 2020 (Nat Biotechnol 38:365, doi:10.1038/s41587-019-0344-3) human phosphoserine, phosphothreonine and '
           'phosphotyrosine sites (Supplementary Tables 3 and 2), checked on the reference proteome, one row per residue, with the '
           'functional score, best localisation probability, best PEP, number of biological samples and spectral counts')


def num(value, fmt='.9g') -> str:
    return '' if value is None or value == 'NA' or pd.isna(value) else format(float(value), fmt)


def classify(acc: str, position, residue, reference: dict[str, str]) -> tuple[str, str]:
    """(status, detail): kept, or the reason for the rejection."""
    acc = str(acc).strip()
    if '-' in acc:
        return 'isoform_accession', acc
    seq = reference.get(acc)
    if seq is None:
        return 'accession_not_in_reference', acc
    if residue not in tuple(RESIDUES):
        return 'bad_residue', str(residue)
    try:
        pos = int(position)
    except (TypeError, ValueError):
        return 'bad_position', str(position)
    if pos < 1 or pos > len(seq):
        return 'bad_position', f'{pos} of {len(seq)}'
    if seq[pos - 1] != residue:
        return 'residue_mismatch', f'table {residue}, canonical {seq[pos - 1]} at {pos}'
    return 'kept', ''


def stage(table: pd.DataFrame, reference: dict[str, str]) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """The validated sites, the rejected sites and the counts, from the merged tables 3 and 2 (one row per site, the
    columns uniprot, position, residue, functional_score, best_PEP, best_localization_prob, biological_samples, spectralcounts)."""
    kept, rejected, counts = [], [], Counter()
    for r in table.itertuples(index=False):
        status, detail = classify(r.uniprot, r.position, r.residue, reference)
        counts[status] += 1
        if status != 'kept':
            rejected.append({'uniprot': r.uniprot, 'position': r.position, 'residue': r.residue,
                             'functional_score': num(r.functional_score), 'reason': status, 'detail': detail})
            continue
        kept.append({'canonical_UniProtAC': str(r.uniprot).strip(), 'position': int(r.position), 'residue': r.residue,
                     'mark': 'phosphorylation', 'state': 'phospho', 'evidence': '', 'evidence_class': 'high-throughput',
                     'pmids': '', 'source': SOURCE, 'source_release': SOURCE_RELEASE,
                     'functional_score': num(r.functional_score), 'best_localization_prob': num(r.best_localization_prob),
                     'best_PEP': num(r.best_PEP), 'biological_samples': num(r.biological_samples, '.0f'),
                     'spectral_counts': num(r.spectralcounts, '.0f')})
    sites = pd.DataFrame(kept, columns=COLUMNS)
    sites = sites.sort_values(['canonical_UniProtAC', 'position']).reset_index(drop=True)
    bad = pd.DataFrame(rejected, columns=REJECTED_COLUMNS)
    if len(bad):
        bad = bad.sort_values(['reason', 'uniprot', 'position'], key=lambda c: c.astype(str)).reset_index(drop=True)
    return sites, bad, dict(counts)


def read_sheet(path: Path, sheet: str, needed: list[str]) -> pd.DataFrame:
    it = openpyxl.load_workbook(path, read_only=True)[sheet].iter_rows(values_only=True)
    head = list(next(it))
    missing = [c for c in needed if c not in head]
    if missing:
        raise SystemExit(f'{path} {sheet}: columns missing: {missing}')
    keep = [head.index(c) for c in needed]
    return pd.DataFrame([[r[i] for i in keep] for r in it if r[0] is not None], columns=needed)


def merge_tables(scores: pd.DataFrame, features: pd.DataFrame) -> pd.DataFrame:
    """Tables 3 and 2 joined on uniprot and position; stops unless both hold the same sites, each once."""
    for name, t in (('Supplementary Table 3', scores), ('Supplementary Table 2', features)):
        if t.duplicated(['uniprot', 'position']).any():
            raise SystemExit(f'{name}: a site is listed twice')
    merged = scores.merge(features, on=['uniprot', 'position'], how='outer', indicator=True, validate='one_to_one')
    if (merged['_merge'] != 'both').any():
        raise SystemExit(f'Supplementary Tables 3 and 2 hold different sites: {merged["_merge"].value_counts().to_dict()}')
    return merged.drop(columns='_merge')


def quantiles(values: pd.Series) -> dict:
    q = values.quantile([0, .01, .05, .25, .5, .75, .95, .99, 1])
    return {f'{k:g}': round(float(v), 6) for k, v in q.items()}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--scores', required=True, help=f'Supplementary Table 3 ({RAW_SCORES})')
    ap.add_argument('--features', required=True, help=f'Supplementary Table 2 ({RAW_FEATURES})')
    ap.add_argument('--fasta', required=True, help='the reference proteome (canonical reviewed accessions)')
    ap.add_argument('--out-dir', default='inputs/marks')
    ap.add_argument('--manifest', default='inputs/MANIFEST.tsv')
    ap.add_argument('--downloaded', default='', help='date of the download for raw_files.tsv (default: today)')
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    raws = {}
    for key, path, name in (('scores', Path(args.scores), RAW_SCORES), ('features', Path(args.features), RAW_FEATURES)):
        data = path.read_bytes()
        raws[key] = {'name': name, 'path': path, 'bytes': len(data), 'md5': hashlib.md5(data).hexdigest(),
                     'sha256': hashlib.sha256(data).hexdigest()}
    table = merge_tables(read_sheet(raws['scores']['path'], SHEET_SCORES, NEEDED_SCORES),
                         read_sheet(raws['features']['path'], SHEET_FEATURES, NEEDED_FEATURES))
    reference = parse_fasta(args.fasta)
    sites, rejected, counts = stage(table, reference)

    out_dir.mkdir(parents=True, exist_ok=True)
    sites.to_csv(out_dir / OUT_SITES, sep='\t', index=False, lineterminator='\n')
    rejected.to_csv(out_dir / OUT_REJECTED, sep='\t', index=False, lineterminator='\n')
    score = pd.to_numeric(sites['functional_score'])
    summary = {
        'source': f'{CITATION}; Supplementary Table 3 (functional scores) and Supplementary Table 2 (residue and evidence columns)',
        'article_url': ARTICLE_URL,
        'urls': {'scores': URL_SCORES, 'features': URL_FEATURES},
        'licence': LICENSE,
        'redistribution': 'checksum_only',
        'reference_fasta': args.fasta,
        'reference_fasta_md5': md5_file(args.fasta),
        'raw': {k: {f: v[f] for f in ('name', 'md5', 'sha256')} for k, v in raws.items()},
        'input_sites': len(table),
        'input_proteins': int(table['uniprot'].nunique()),
        'input_sites_by_residue': {k: int(v) for k, v in sorted(table['residue'].value_counts().items())},
        'row_status': dict(sorted(counts.items())),
        'rejected_by_reason': {k: v for k, v in sorted(counts.items()) if k != 'kept'},
        'relocated_sites': 0,
        'relocation_note': 'no peptide or window is given, so a site whose residue does not match the canonical sequence is rejected, not moved',
        'validated_sites': len(sites),
        'validated_proteins': int(sites['canonical_UniProtAC'].nunique()),
        'validated_sites_by_residue': {k: int(v) for k, v in sites['residue'].value_counts().sort_index().items()},
        'validated_proteins_by_residue': {k: int(v) for k, v in sites.groupby('residue')['canonical_UniProtAC'].nunique().sort_index().items()},
        'functional_score_quantiles': quantiles(score),
        'functional_score_ge_0.5': int((score >= 0.5).sum()),
        'functional_score_ge_0.5_fraction': round(float((score >= 0.5).mean()), 6),
        'files': {'validated': OUT_SITES, 'rejected': OUT_REJECTED},
        'files_md5': {'validated': md5_file(out_dir / OUT_SITES), 'rejected': md5_file(out_dir / OUT_REJECTED)},
    }
    (out_dir / OUT_SUMMARY).write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n', encoding='utf-8')

    when = args.downloaded or date.today().isoformat()
    raw_list = out_dir / 'raw_files.tsv'
    upsert_tsv_rows(raw_list, [
        {'file': raws[k]['name'], 'url': u, 'bytes': raws[k]['bytes'], 'md5': raws[k]['md5'], 'sha256': raws[k]['sha256'],
         'downloaded': when, 'version': f'{CITATION}; {what}',
         'kept': f"not in Git; staged as marks/{OUT_SITES} (checksum_only: {LICENSE})"}
        for k, u, what in (('scores', URL_SCORES, 'Supplementary Table 3 (functional scores)'),
                           ('features', URL_FEATURES, 'Supplementary Table 2 (annotated features); only residue and evidence columns are used'))])

    command = (f'python scripts/58_stage_ochoa_phosphosites.py --scores PTM_data/marks/{RAW_SCORES} --features PTM_data/marks/{RAW_FEATURES} '
               f'--fasta {args.fasta} --out-dir inputs/marks --downloaded {when}')
    check = ('; residue S/T/Y checked on the reference sequence, mismatches rejected with a reason (no window to relocate by); '
             f'counts in marks/{OUT_SUMMARY}')
    derived = f'{RAW_SCORES} (MD5 {raws["scores"]["md5"]}), {RAW_FEATURES} (MD5 {raws["features"]["md5"]})'
    common_fields = {'source': URL_SCORES, 'version_or_date': f'{SOURCE_RELEASE}; downloaded {when}', 'citation': CITATION,
                     'licence': LICENSE, 'derived_from': derived, 'how_made_or_checked': command + check}
    update_manifest(args.manifest, [
        {'file': f'marks/{OUT_SITES}', 'bytes': (out_dir / OUT_SITES).stat().st_size, 'md5': md5_file(out_dir / OUT_SITES),
         'rows': len(sites), 'redistribution': 'checksum_only', 'content': CONTENT + f': {len(sites):,} rows', **common_fields},
        {'file': f'marks/{OUT_REJECTED}', 'bytes': (out_dir / OUT_REJECTED).stat().st_size, 'md5': md5_file(out_dir / OUT_REJECTED),
         'rows': len(rejected), 'redistribution': 'checksum_only',
         'content': f'Sites of Ochoa et al. 2020 Supplementary Table 3 that failed the reference check, with the reason: {len(rejected):,} rows',
         **common_fields},
        {'file': f'marks/{OUT_SUMMARY}', 'bytes': (out_dir / OUT_SUMMARY).stat().st_size, 'md5': md5_file(out_dir / OUT_SUMMARY),
         'redistribution': 'git', 'content': 'Counts, drops by reason, functional-score quantiles and the MD5/SHA-256 of the inputs and MD5 of both tables of the Ochoa 2020 staging',
         **{**common_fields, 'source': 'the tables above', 'derived_from': 'marks/raw_files.tsv', 'how_made_or_checked': command}},
    ])
    refresh_manifest_row(args.manifest, raw_list, 'marks/raw_files.tsv', redistribution='git',
                         content='URL, size, MD5, SHA-256 and date of every raw download behind inputs/marks/',
                         source='the URLs it lists', derived_from='',
                         how_made_or_checked='written by the scripts/54, 55 and 58 staging runs')
    print(json.dumps({k: summary[k] for k in ('input_sites', 'input_proteins', 'row_status', 'validated_sites', 'validated_proteins',
                                               'validated_sites_by_residue', 'validated_proteins_by_residue',
                                               'functional_score_quantiles', 'functional_score_ge_0.5_fraction')}, indent=1))


if __name__ == '__main__':
    main()

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from common import canonicalize_uniprot_accession, normalize_accession, save_table


def iter_json_array(path: str | Path, chunk_size: int = 1024 * 1024):
    decoder = json.JSONDecoder()
    buffer = ''
    started = False
    path = Path(path)
    with path.open('rt', encoding='utf-8') as fh:
        while True:
            chunk = fh.read(chunk_size)
            eof = chunk == ''
            buffer += chunk
            while True:
                buffer = buffer.lstrip()
                if not started:
                    if not buffer:
                        break
                    if not buffer.startswith('['):
                        raise ValueError(f'{path} does not contain a top-level JSON array')
                    buffer = buffer[1:]
                    started = True
                    continue
                if not buffer:
                    break
                if buffer.startswith(']'):
                    return
                try:
                    obj, idx = decoder.raw_decode(buffer)
                except json.JSONDecodeError:
                    if eof:
                        raise
                    break
                yield obj
                buffer = buffer[idx:].lstrip()
                if buffer.startswith(','):
                    buffer = buffer[1:]
                    continue
                if buffer.startswith(']'):
                    return
            if eof:
                break


def load_accessions(path: str | Path) -> set[str]:
    path = Path(path)
    if path.suffix.lower() == '.tsv':
        df = pd.read_csv(path, sep='\t', low_memory=False)
        for column in ['Entry', 'canonical_UniProtAC', 'substrate_UniProtAC', 'accession']:
            if column in df.columns:
                return set(df[column].dropna().astype(str).map(canonicalize_uniprot_accession))
        raise ValueError(f'no accession column found in {path}')
    if path.suffix.lower() in {'.txt', '.list'}:
        return {canonicalize_uniprot_accession(line.strip()) for line in path.read_text().splitlines() if line.strip()}
    raise ValueError(f'unsupported accession source: {path}')


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True)
    ap.add_argument('--outdir', default='data/context')
    ap.add_argument('--annotation-key', default='prediction-disorder-th_50')
    ap.add_argument('--proteome', default='UP000005640')
    ap.add_argument('--reviewed-only', action='store_true')
    ap.add_argument('--accessions', default='')
    ap.add_argument('--chunk-size', type=int, default=1024 * 1024)
    args = ap.parse_args()

    wanted_accessions = load_accessions(args.accessions) if args.accessions else None
    rows = []
    protein_rows = []
    total_entries = 0
    kept_entries = 0
    annotated_entries = 0

    for entry in iter_json_array(args.input, chunk_size=args.chunk_size):
        total_entries += 1
        accession = canonicalize_uniprot_accession(entry.get('acc', ''))
        if not accession:
            continue
        if args.proteome and str(entry.get('proteome', '')).strip() != args.proteome:
            continue
        if args.reviewed_only and not bool(entry.get('reviewed', False)):
            continue
        if wanted_accessions is not None and accession not in wanted_accessions:
            continue
        kept_entries += 1
        annotation = entry.get(args.annotation_key)
        if not isinstance(annotation, dict):
            continue
        regions = annotation.get('regions') or []
        if not regions:
            continue
        annotated_entries += 1
        protein_rows.append({
            'canonical_UniProtAC': accession,
            'mobidb_annotation_key': args.annotation_key,
            'region_count': len(regions),
            'content_fraction': annotation.get('content_fraction', pd.NA),
            'content_count': annotation.get('content_count', pd.NA),
            'source_id': annotation.get('source_id', ''),
        })
        for region_index, region in enumerate(regions, start=1):
            if not isinstance(region, list) or len(region) < 2:
                continue
            start = pd.to_numeric(pd.Series([region[0]]), errors='coerce').iloc[0]
            end = pd.to_numeric(pd.Series([region[1]]), errors='coerce').iloc[0]
            if pd.isna(start) or pd.isna(end):
                continue
            rows.append({
                'canonical_UniProtAC': accession,
                'fragment_start': int(start),
                'fragment_end': int(end),
                'mobidb_annotation_key': args.annotation_key,
                'mobidb_region_index': region_index,
                'source_id': annotation.get('source_id', ''),
                'content_fraction': annotation.get('content_fraction', pd.NA),
                'content_count': annotation.get('content_count', pd.NA),
            })

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    intervals = pd.DataFrame(rows).sort_values(['canonical_UniProtAC', 'fragment_start', 'fragment_end']) if rows else pd.DataFrame(columns=[
        'canonical_UniProtAC', 'fragment_start', 'fragment_end', 'mobidb_annotation_key',
        'mobidb_region_index', 'source_id', 'content_fraction', 'content_count',
    ])
    proteins = pd.DataFrame(protein_rows).sort_values(['canonical_UniProtAC']) if protein_rows else pd.DataFrame(columns=[
        'canonical_UniProtAC', 'mobidb_annotation_key', 'region_count', 'content_fraction', 'content_count', 'source_id',
    ])
    save_table(intervals, outdir / 'mobidb_human_reviewed_disorder_intervals.tsv')
    save_table(proteins, outdir / 'mobidb_human_reviewed_disorder_protein_summary.tsv')

    summary = pd.DataFrame([{
        'input': str(Path(args.input)),
        'annotation_key': args.annotation_key,
        'proteome_filter': args.proteome,
        'reviewed_only': bool(args.reviewed_only),
        'accession_filter_count': len(wanted_accessions) if wanted_accessions is not None else pd.NA,
        'entries_total': total_entries,
        'entries_after_filters': kept_entries,
        'entries_with_annotation': annotated_entries,
        'interval_rows': len(intervals),
    }])
    save_table(summary, outdir / 'mobidb_human_reviewed_disorder_summary.tsv')
    print(summary.to_dict(orient='records')[0])


if __name__ == '__main__':
    main()

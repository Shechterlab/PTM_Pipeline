from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from common import canonicalize_uniprot_accession, save_table


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
    ap.add_argument('--outdir', default='PTM_data/context')
    ap.add_argument('--annotation-key', default='prediction-disorder-th_50')
    ap.add_argument('--proteome', default='UP000005640')
    ap.add_argument('--reviewed-only', action='store_true')
    ap.add_argument('--accessions', default='')
    ap.add_argument('--chunk-size', type=int, default=1024 * 1024)
    ap.add_argument('--fasta', default='', help='reference FASTA: entries whose MobiDB sequence differs are left out (their '
                    'regions number another sequence) and listed in mobidb_entries_other_sequence.tsv')
    ap.add_argument('--patch', default='', help='instead of reading an export: check the staged intervals table PATCH '
                    'against --fasta and replace the entries numbered on another sequence (see patch_staged)')
    args = ap.parse_args()
    reference = read_fasta(Path(args.fasta)) if args.fasta else {}
    if args.patch:
        if not reference:
            raise SystemExit('--patch needs --fasta')
        patch_staged(Path(args.patch), reference, args.annotation_key)
        return

    wanted_accessions = load_accessions(args.accessions) if args.accessions else None
    other_sequence = []
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
        if reference and entry.get('sequence', '') != reference.get(accession, entry.get('sequence', '')):
            other_sequence.append({'canonical_UniProtAC': accession, 'mobidb_length': len(entry.get('sequence', '')),
                                   'reference_length': len(reference[accession])})
            continue
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
    if reference:
        save_table(pd.DataFrame(other_sequence, columns=['canonical_UniProtAC', 'mobidb_length', 'reference_length']),
                   outdir / 'mobidb_entries_other_sequence.tsv')
    print(summary.to_dict(orient='records')[0], f'entries on another sequence: {len(other_sequence)}')


def read_fasta(path: Path) -> dict[str, str]:
    seqs, acc = {}, None
    for line in path.open(encoding='utf-8'):
        line = line.strip()
        if line.startswith('>'):
            parts = line[1:].split('|')
            acc = parts[1] if len(parts) > 2 else line[1:].split()[0]
            seqs[acc] = ''
        elif acc:
            seqs[acc] += line
    return seqs


def implied_length_range(count, fraction) -> tuple[float, float]:
    """The sequence lengths a MobiDB entry's content_count and content_fraction (rounded to 3 decimals) allow."""
    count, fraction = float(count), float(fraction)
    return count / (fraction + 0.0005), (count / (fraction - 0.0005) if fraction > 0.0005 else float('inf'))


def patch_staged(path: Path, reference: dict, key: str) -> None:
    """An already staged table whose export is not kept: an entry whose regions run past the reference sequence, or
    whose content_count and content_fraction imply another length, numbers another sequence. Each such entry is
    fetched from the MobiDB API (https://mobidb.org/api/download?acc=...&format=json): its regions replace the staged
    ones when the API's sequence is the reference sequence, otherwise the entry is dropped. The table is rewritten in
    place and the entries are listed in mobidb_entries_patched.tsv next to it."""
    import datetime as dt
    import time
    import urllib.request
    table = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    wrong = []
    for acc, g in table.groupby('canonical_UniProtAC', sort=False):
        ref = reference.get(acc, '')
        lo, hi = implied_length_range(g['content_count'].iloc[0], g['content_fraction'].iloc[0])
        if not ref or g['fragment_end'].astype(int).max() > len(ref) or not lo - 0.5 <= len(ref) <= hi + 0.5:
            wrong.append(acc)
    rows, report = [], []
    for acc in wrong:
        url = f'https://mobidb.org/api/download?acc={acc}&format=json'
        for attempt in range(4):
            try:
                with urllib.request.urlopen(url, timeout=60) as response:
                    body = response.read().decode().strip()
                entry = json.loads(body) if body else {}   # MobiDB answers nothing for an entry it no longer has
                break
            except OSError:
                if attempt == 3:
                    raise
                time.sleep(2 ** (attempt + 1))
        entry = (entry[0] if entry else {}) if isinstance(entry, list) else entry
        same = entry.get('sequence', '') == reference.get(acc)
        annotation = entry.get(key) if same else None
        regions = (annotation or {}).get('regions') or []
        report.append({'canonical_UniProtAC': acc, 'reference_length': len(reference.get(acc, '')),
                       'api_length': len(entry.get('sequence', '')), 'outcome': 'replaced' if same else 'dropped',
                       'regions': len(regions)})
        for k, (a, b) in enumerate(regions, start=1):
            rows.append({'canonical_UniProtAC': acc, 'fragment_start': str(a), 'fragment_end': str(b), 'mobidb_annotation_key': key,
                         'mobidb_region_index': str(k), 'source_id': annotation.get('source_id', '') or '',
                         'content_fraction': str(annotation.get('content_fraction', '')), 'content_count': str(annotation.get('content_count', ''))})
    out = pd.concat([table[~table['canonical_UniProtAC'].isin(wrong)], pd.DataFrame(rows, columns=table.columns)], ignore_index=True)
    out = out.sort_values(['canonical_UniProtAC', 'fragment_start', 'fragment_end'], key=lambda c: c.astype(int) if c.name != 'canonical_UniProtAC' else c)
    out.to_csv(path, sep='\t', index=False, lineterminator='\n')
    rep = pd.DataFrame(report)
    rep['fetched'] = dt.date.today().isoformat()
    rep.to_csv(path.parent / 'mobidb_entries_patched.tsv', sep='\t', index=False, lineterminator='\n')
    print(f'{len(wrong)} entries numbered on another sequence: {(rep["outcome"] == "replaced").sum() if len(rep) else 0} replaced from the API, '
          f'{(rep["outcome"] == "dropped").sum() if len(rep) else 0} dropped -> {path}')


if __name__ == '__main__':
    main()

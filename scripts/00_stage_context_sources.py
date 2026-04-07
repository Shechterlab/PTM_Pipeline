from __future__ import annotations

import argparse
import gzip
import io
import json
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO

import pandas as pd
import requests


UNIPROT_QUERY = '(proteome:UP000005640) AND (reviewed:true)'
UNIPROT_FASTA_URL = 'https://rest.uniprot.org/uniprotkb/stream'
UNIPROT_TSV_URL = 'https://rest.uniprot.org/uniprotkb/stream'
INTERPRO_ENTRY_URL = 'https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/reviewed/{accession}'
INTERPRO_ENTRY_LIST_URL = 'https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list'
INTERPRO_PROTEIN2IPR_URL = 'https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/protein2ipr.dat.gz'
DEFAULT_DOMAIN_TYPES = ('domain', 'homologous_superfamily', 'repeat')
HEADERS = {'User-Agent': 'PTM-pipeline-context-stager/1.0'}


def download_file(url: str, params: dict[str, str], dest: Path, timeout: int = 300) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, params=params, stream=True, timeout=timeout, headers=HEADERS) as response:
        response.raise_for_status()
        with open(dest, 'wb') as fh:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    fh.write(chunk)


def parse_accessions_from_tsv(path: Path) -> list[str]:
    df = pd.read_csv(path, sep='\t')
    if 'Entry' not in df.columns:
        raise ValueError(f'expected UniProt metadata TSV with Entry column: {path}')
    accessions = df['Entry'].dropna().astype(str).str.strip()
    return sorted(set(acc for acc in accessions if acc))


def load_accessions(path: Path) -> list[str]:
    if path.suffix.lower() in {'.tsv', '.txt', '.csv'}:
        if path.suffix.lower() == '.tsv':
            df = pd.read_csv(path, sep='\t', low_memory=False)
        elif path.suffix.lower() == '.csv':
            df = pd.read_csv(path, low_memory=False)
        else:
            lines = [line.strip() for line in path.read_text().splitlines()]
            return sorted(set(line for line in lines if line and not line.startswith('#')))
        for column in ['canonical_UniProtAC', 'substrate_UniProtAC', 'Entry', 'accession']:
            if column in df.columns:
                accessions = df[column].dropna().astype(str).str.strip()
                return sorted(set(acc for acc in accessions if acc))
        raise ValueError(f'no accession column found in {path}')
    if path.suffix.lower() in {'.fa', '.fasta'}:
        accessions = []
        for line in path.read_text().splitlines():
            if not line.startswith('>'):
                continue
            header = line[1:].strip()
            parts = header.split('|')
            if len(parts) >= 2:
                accessions.append(parts[1].strip())
        return sorted(set(acc for acc in accessions if acc))
    raise ValueError(f'unsupported accession source: {path}')


def load_protein_lengths(path: Path) -> dict[str, int]:
    df = pd.read_csv(path, sep='\t', usecols=['Entry', 'Length'])
    df['Entry'] = df['Entry'].astype(str).str.strip()
    df['Length'] = pd.to_numeric(df['Length'], errors='coerce')
    df = df[df['Entry'] != '']
    df = df[df['Length'].notna()].copy()
    df['Length'] = df['Length'].astype(int)
    return dict(zip(df['Entry'], df['Length']))


def normalize_interpro_type(value: str) -> str:
    return str(value or '').strip().lower().replace(' ', '_').replace('-', '_')


def infer_member_database_source(member_accession: str) -> str:
    text = str(member_accession or '').strip()
    if not text:
        return ''
    if ':' in text:
        return text.split(':', 1)[0]
    match = re.match(r'[A-Za-z]+', text)
    return match.group(0) if match else ''


@contextmanager
def open_text_stream(source: str | Path, timeout: int = 300) -> Iterator[TextIO]:
    source_text = str(source)
    if source_text.startswith(('http://', 'https://')):
        with requests.get(source_text, stream=True, timeout=timeout, headers=HEADERS) as response:
            response.raise_for_status()
            if source_text.endswith('.gz'):
                with gzip.GzipFile(fileobj=response.raw) as gz_fh:
                    with io.TextIOWrapper(gz_fh, encoding='utf-8') as text_fh:
                        yield text_fh
            else:
                with io.TextIOWrapper(response.raw, encoding='utf-8') as text_fh:
                    yield text_fh
        return

    path = Path(source)
    if path.suffix.lower() == '.gz':
        with gzip.open(path, 'rt', encoding='utf-8') as text_fh:
            yield text_fh
    else:
        with path.open('rt', encoding='utf-8') as text_fh:
            yield text_fh


def fetch_interpro_payload(accession: str, cache_dir: Path, page_size: int, force: bool, timeout: int, retries: int) -> dict:
    cache_path = cache_dir / f'{accession}.json'
    if cache_path.exists() and not force:
        return json.loads(cache_path.read_text())

    url = INTERPRO_ENTRY_URL.format(accession=accession)
    results: list[dict] = []
    count = 0
    last_error: Exception | None = None

    for attempt in range(retries):
        try:
            next_url = url
            params: dict[str, str] | None = {'page_size': str(page_size)}
            first_page = True
            results = []
            while next_url:
                with requests.get(next_url, params=params if first_page else None, timeout=timeout, headers=HEADERS) as response:
                    response.raise_for_status()
                    payload = response.json()
                count = int(payload.get('count', 0) or 0)
                results.extend(payload.get('results', []))
                next_url = payload.get('next')
                first_page = False
                params = None
            break
        except requests.HTTPError as exc:
            status_code = exc.response.status_code if exc.response is not None else None
            if status_code == 404:
                count = 0
                results = []
                last_error = None
                break
            last_error = exc
        except Exception as exc:  # pragma: no cover - defensive network retry path
            last_error = exc
        time.sleep(min(30, 2 ** attempt))
    else:
        raise RuntimeError(f'failed to fetch InterPro payload for {accession}: {last_error}')

    payload = {
        'accession': accession,
        'count': count,
        'results': results,
    }
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(payload))
    return payload


def member_database_fields(member_databases: dict | None) -> tuple[str, str]:
    if not member_databases:
        return '', ''
    sources = sorted(str(source) for source in member_databases)
    ids = []
    for source, members in member_databases.items():
        for accession in sorted(members):
            ids.append(f'{source}:{accession}')
    return ';'.join(sources), ';'.join(ids)


def payload_to_rows(payload: dict, default_domain_types: set[str]) -> list[dict]:
    accession = str(payload.get('accession', '')).upper()
    rows: list[dict] = []
    for result in payload.get('results', []):
        metadata = result.get('metadata', {})
        member_sources, member_accessions = member_database_fields(metadata.get('member_databases'))
        interpro_type = normalize_interpro_type(metadata.get('type', ''))
        for protein in result.get('proteins', []):
            protein_length = protein.get('protein_length')
            for location_index, location in enumerate(protein.get('entry_protein_locations', []), start=1):
                fragments = location.get('fragments') or []
                for fragment_index, fragment in enumerate(fragments, start=1):
                    rows.append({
                        'canonical_UniProtAC': accession,
                        'protein_length': protein_length,
                        'interpro_accession': metadata.get('accession', ''),
                        'interpro_name': metadata.get('name', ''),
                        'interpro_type': interpro_type,
                        'interpro_source_database': metadata.get('source_database', ''),
                        'member_database_sources': member_sources,
                        'member_database_accessions': member_accessions,
                        'fragment_start': int(fragment.get('start')),
                        'fragment_end': int(fragment.get('end')),
                        'fragment_status': fragment.get('dc-status', ''),
                        'location_index': location_index,
                        'fragment_index': fragment_index,
                        'location_representative': bool(location.get('representative', False)),
                        'location_model': location.get('model', ''),
                        'location_score': location.get('score', ''),
                        'domain_track_default': interpro_type in default_domain_types,
                        'interpro_interval_source': 'api',
                    })
    return rows


def load_entry_type_map(source: str | Path, timeout: int) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with open_text_stream(source, timeout=timeout) as fh:
        header_seen = False
        for raw_line in fh:
            line = raw_line.rstrip('\n')
            if not line:
                continue
            if not header_seen:
                header_seen = True
                if line.startswith('ENTRY_AC\tENTRY_TYPE'):
                    continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            mapping[parts[0].strip()] = normalize_interpro_type(parts[1])
    return mapping


def build_rows_from_bulk(
    accessions: list[str],
    protein_lengths: dict[str, int],
    entry_type_map: dict[str, str],
    source: str | Path,
    default_domain_types: set[str],
    timeout: int,
    progress_every: int,
) -> tuple[list[dict], int, int, int]:
    accession_set = set(accessions)
    if not accession_set:
        return [], 0, 0, 0
    min_target = min(accession_set)
    max_target = max(accession_set)
    proteins_with_entries: set[str] = set()
    rows: list[dict] = []
    lines_scanned = 0

    with open_text_stream(source, timeout=timeout) as fh:
        for raw_line in fh:
            lines_scanned += 1
            if progress_every and lines_scanned % progress_every == 0:
                print(f'scanned {lines_scanned:,} protein2ipr rows; kept {len(rows):,}')

            parts = raw_line.rstrip('\n').split('\t')
            if len(parts) < 6:
                continue
            accession = parts[0].strip().upper()
            if accession < min_target:
                continue
            if accession > max_target:
                break
            if accession not in accession_set:
                continue

            interpro_accession = parts[1].strip()
            interpro_type = entry_type_map.get(interpro_accession, '')
            member_accession = parts[3].strip()
            try:
                fragment_start = int(parts[4])
                fragment_end = int(parts[5])
            except ValueError:
                continue

            proteins_with_entries.add(accession)
            rows.append({
                'canonical_UniProtAC': accession,
                'protein_length': protein_lengths.get(accession, pd.NA),
                'interpro_accession': interpro_accession,
                'interpro_name': parts[2].strip(),
                'interpro_type': interpro_type,
                'interpro_source_database': 'InterPro',
                'member_database_sources': infer_member_database_source(member_accession),
                'member_database_accessions': member_accession,
                'fragment_start': fragment_start,
                'fragment_end': fragment_end,
                'fragment_status': '',
                'location_index': pd.NA,
                'fragment_index': pd.NA,
                'location_representative': pd.NA,
                'location_model': '',
                'location_score': '',
                'domain_track_default': interpro_type in default_domain_types,
                'interpro_interval_source': 'bulk_protein2ipr',
            })

    proteins_without_entries = len(accession_set - proteins_with_entries)
    return rows, len(proteins_with_entries), proteins_without_entries, lines_scanned


def build_rows_from_bulk_pandas(
    accessions: list[str],
    protein_lengths: dict[str, int],
    entry_type_map: dict[str, str],
    source: str | Path,
    default_domain_types: set[str],
    chunk_size: int,
    progress_every: int,
) -> tuple[list[dict], int, int, int]:
    source_path = Path(source)
    accession_set = {str(accession).strip().upper() for accession in accessions if str(accession).strip()}
    if not accession_set:
        return [], 0, 0, 0

    column_names = [
        'canonical_UniProtAC',
        'interpro_accession',
        'interpro_name',
        'member_database_accessions',
        'fragment_start',
        'fragment_end',
    ]
    rows: list[dict] = []
    proteins_with_entries: set[str] = set()
    lines_scanned = 0

    for chunk in pd.read_csv(
        source_path,
        sep='\t',
        header=None,
        names=column_names,
        usecols=range(6),
        compression='infer',
        dtype=str,
        chunksize=max(1, chunk_size),
        low_memory=False,
    ):
        lines_scanned += len(chunk)
        if progress_every and lines_scanned % progress_every < len(chunk):
            print(f'scanned {lines_scanned:,} protein2ipr rows; kept {len(rows):,}')

        chunk['canonical_UniProtAC'] = chunk['canonical_UniProtAC'].astype(str).str.strip().str.upper()
        filtered = chunk[chunk['canonical_UniProtAC'].isin(accession_set)].copy()
        if filtered.empty:
            continue

        filtered['fragment_start'] = pd.to_numeric(filtered['fragment_start'], errors='coerce')
        filtered['fragment_end'] = pd.to_numeric(filtered['fragment_end'], errors='coerce')
        filtered = filtered[filtered['fragment_start'].notna() & filtered['fragment_end'].notna()].copy()
        if filtered.empty:
            continue

        filtered['fragment_start'] = filtered['fragment_start'].astype(int)
        filtered['fragment_end'] = filtered['fragment_end'].astype(int)
        filtered['interpro_accession'] = filtered['interpro_accession'].astype(str).str.strip()
        filtered['interpro_name'] = filtered['interpro_name'].astype(str).str.strip()
        filtered['member_database_accessions'] = filtered['member_database_accessions'].astype(str).str.strip()
        filtered['protein_length'] = filtered['canonical_UniProtAC'].map(protein_lengths)
        filtered['interpro_type'] = filtered['interpro_accession'].map(entry_type_map).fillna('')
        filtered['interpro_source_database'] = 'InterPro'
        filtered['member_database_sources'] = filtered['member_database_accessions'].map(infer_member_database_source)
        filtered['fragment_status'] = ''
        filtered['location_index'] = pd.NA
        filtered['fragment_index'] = pd.NA
        filtered['location_representative'] = pd.NA
        filtered['location_model'] = ''
        filtered['location_score'] = ''
        filtered['domain_track_default'] = filtered['interpro_type'].isin(default_domain_types)
        filtered['interpro_interval_source'] = 'bulk_protein2ipr'

        proteins_with_entries.update(filtered['canonical_UniProtAC'].tolist())
        rows.extend(
            filtered[[
                'canonical_UniProtAC',
                'protein_length',
                'interpro_accession',
                'interpro_name',
                'interpro_type',
                'interpro_source_database',
                'member_database_sources',
                'member_database_accessions',
                'fragment_start',
                'fragment_end',
                'fragment_status',
                'location_index',
                'fragment_index',
                'location_representative',
                'location_model',
                'location_score',
                'domain_track_default',
                'interpro_interval_source',
            ]].to_dict(orient='records')
        )

    proteins_without_entries = len(accession_set - proteins_with_entries)
    return rows, len(proteins_with_entries), proteins_without_entries, lines_scanned


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--outdir', default='data/context')
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--page-size', type=int, default=200)
    ap.add_argument('--timeout', type=int, default=120)
    ap.add_argument('--retries', type=int, default=5)
    ap.add_argument('--max-accessions', type=int, default=0)
    ap.add_argument('--force-interpro', action='store_true')
    ap.add_argument('--force-uniprot', action='store_true')
    ap.add_argument('--skip-interpro', action='store_true')
    ap.add_argument('--accessions', default='')
    ap.add_argument('--domain-types', nargs='*', default=list(DEFAULT_DOMAIN_TYPES))
    ap.add_argument('--interpro-source', choices=['auto', 'api', 'bulk'], default='auto')
    ap.add_argument('--interpro-protein2ipr', default='')
    ap.add_argument('--interpro-entry-list', default='')
    ap.add_argument('--bulk-parser', choices=['auto', 'pandas', 'stream'], default='auto')
    ap.add_argument('--bulk-chunk-size', type=int, default=1_000_000)
    ap.add_argument('--bulk-progress-every', type=int, default=5_000_000)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cache_dir = outdir / 'interpro_cache'

    fasta_path = outdir / 'uniprot_human_reviewed_canonical.fasta'
    metadata_path = outdir / 'uniprot_human_reviewed_canonical_metadata.tsv'

    if args.force_uniprot or not fasta_path.exists():
        print(f'downloading UniProt FASTA -> {fasta_path}')
        download_file(
            UNIPROT_FASTA_URL,
            {'format': 'fasta', 'query': UNIPROT_QUERY},
            fasta_path,
            timeout=max(args.timeout, 300),
        )
    else:
        print(f'reusing UniProt FASTA -> {fasta_path}')
    if args.force_uniprot or not metadata_path.exists():
        print(f'downloading UniProt metadata -> {metadata_path}')
        download_file(
            UNIPROT_TSV_URL,
            {'format': 'tsv', 'query': UNIPROT_QUERY, 'fields': 'accession,id,gene_primary,protein_name,length'},
            metadata_path,
            timeout=max(args.timeout, 300),
        )
    else:
        print(f'reusing UniProt metadata -> {metadata_path}')

    if args.accessions:
        accessions = load_accessions(Path(args.accessions))
    else:
        accessions = parse_accessions_from_tsv(metadata_path)
    if args.max_accessions:
        accessions = accessions[:args.max_accessions]

    protein_lengths = load_protein_lengths(metadata_path)
    default_domain_types = {normalize_interpro_type(value) for value in args.domain_types}
    interpro_source_mode = 'skipped'
    interpro_entry_list_source = ''
    interpro_protein2ipr_source = ''
    bulk_lines_scanned = 0

    rows: list[dict] = []
    proteins_with_entries = 0
    proteins_without_entries = 0
    if not args.skip_interpro:
        interpro_source_mode = args.interpro_source
        if interpro_source_mode == 'auto':
            local_bulk = outdir / 'protein2ipr.dat.gz'
            if args.interpro_protein2ipr or local_bulk.exists() or len(accessions) > 250:
                interpro_source_mode = 'bulk'
            else:
                interpro_source_mode = 'api'

        if interpro_source_mode == 'bulk':
            entry_list_path = outdir / 'entry.list'
            if args.interpro_entry_list:
                interpro_entry_list_source = args.interpro_entry_list
            else:
                if args.force_interpro or not entry_list_path.exists():
                    print(f'downloading InterPro entry metadata -> {entry_list_path}')
                    download_file(INTERPRO_ENTRY_LIST_URL, {}, entry_list_path, timeout=max(args.timeout, 300))
                else:
                    print(f'reusing InterPro entry metadata -> {entry_list_path}')
                interpro_entry_list_source = str(entry_list_path)

            local_bulk = outdir / 'protein2ipr.dat.gz'
            if args.interpro_protein2ipr:
                interpro_protein2ipr_source = args.interpro_protein2ipr
            elif local_bulk.exists():
                interpro_protein2ipr_source = str(local_bulk)
            else:
                interpro_protein2ipr_source = INTERPRO_PROTEIN2IPR_URL

            print(f'loading InterPro entry types from {interpro_entry_list_source}')
            entry_type_map = load_entry_type_map(interpro_entry_list_source, timeout=max(args.timeout, 300))
            bulk_parser = args.bulk_parser
            if bulk_parser == 'auto':
                bulk_parser = 'pandas' if not str(interpro_protein2ipr_source).startswith(('http://', 'https://')) else 'stream'
            print(f'loading InterPro protein intervals from {interpro_protein2ipr_source} using {bulk_parser} parser')
            if bulk_parser == 'pandas':
                rows, proteins_with_entries, proteins_without_entries, bulk_lines_scanned = build_rows_from_bulk_pandas(
                    accessions=accessions,
                    protein_lengths=protein_lengths,
                    entry_type_map=entry_type_map,
                    source=interpro_protein2ipr_source,
                    default_domain_types=default_domain_types,
                    chunk_size=args.bulk_chunk_size,
                    progress_every=args.bulk_progress_every,
                )
            else:
                rows, proteins_with_entries, proteins_without_entries, bulk_lines_scanned = build_rows_from_bulk(
                    accessions=accessions,
                    protein_lengths=protein_lengths,
                    entry_type_map=entry_type_map,
                    source=interpro_protein2ipr_source,
                    default_domain_types=default_domain_types,
                    timeout=max(args.timeout, 300),
                    progress_every=args.bulk_progress_every,
                )
        else:
            print(f'fetching InterPro entries for {len(accessions)} accessions via API')
            with ThreadPoolExecutor(max_workers=max(1, args.workers)) as executor:
                futures = {
                    executor.submit(
                        fetch_interpro_payload,
                        accession,
                        cache_dir,
                        args.page_size,
                        args.force_interpro,
                        args.timeout,
                        args.retries,
                    ): accession
                    for accession in accessions
                }
                for future in as_completed(futures):
                    payload = future.result()
                    protein_rows = payload_to_rows(payload, default_domain_types)
                    if protein_rows:
                        proteins_with_entries += 1
                        rows.extend(protein_rows)
                    else:
                        proteins_without_entries += 1
                    if (proteins_with_entries + proteins_without_entries) % 250 == 0:
                        print(f'processed {proteins_with_entries + proteins_without_entries} / {len(accessions)} accessions')

    interpro_all = pd.DataFrame(rows)
    if interpro_all.empty:
        interpro_all = pd.DataFrame(columns=[
            'canonical_UniProtAC', 'protein_length', 'interpro_accession', 'interpro_name', 'interpro_type',
            'interpro_source_database', 'member_database_sources', 'member_database_accessions',
            'fragment_start', 'fragment_end', 'fragment_status', 'location_index', 'fragment_index',
            'location_representative', 'location_model', 'location_score', 'domain_track_default',
            'interpro_interval_source',
        ])
    interpro_all = interpro_all.sort_values([
        'canonical_UniProtAC', 'fragment_start', 'fragment_end', 'interpro_type', 'interpro_accession'
    ])
    interpro_all_path = outdir / 'interpro_human_reviewed_all_entry_intervals.tsv'
    interpro_all.to_csv(interpro_all_path, sep='\t', index=False)

    domain_like = interpro_all[interpro_all['domain_track_default']].copy()
    domain_like_path = outdir / 'interpro_human_reviewed_domain_like_intervals.tsv'
    domain_like.to_csv(domain_like_path, sep='\t', index=False)

    summary = pd.DataFrame([{
        'uniprot_query': UNIPROT_QUERY,
        'fasta_path': str(fasta_path),
        'metadata_path': str(metadata_path),
        'accessions_requested': len(accessions),
        'proteins_with_interpro_entries': proteins_with_entries,
        'proteins_without_interpro_entries': proteins_without_entries,
        'interpro_rows_all': len(interpro_all),
        'interpro_rows_domain_like': len(domain_like),
        'domain_types_default': ';'.join(sorted(default_domain_types)),
        'interpro_cache_dir': str(cache_dir),
        'uniprot_fasta_url': UNIPROT_FASTA_URL,
        'interpro_entry_url_template': INTERPRO_ENTRY_URL,
        'interpro_source_mode': interpro_source_mode,
        'interpro_entry_list_source': interpro_entry_list_source,
        'interpro_protein2ipr_source': interpro_protein2ipr_source,
        'bulk_lines_scanned': bulk_lines_scanned,
    }])
    summary.to_csv(outdir / 'context_source_summary.tsv', sep='\t', index=False)
    print(summary.to_dict(orient='records')[0])


if __name__ == '__main__':
    main()

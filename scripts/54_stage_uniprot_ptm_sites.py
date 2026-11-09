#!/usr/bin/env python3
"""Stage UniProt's methyllysine, phosphorylation and methylarginine annotations as one table per mark (inputs/marks/).

    python scripts/54_stage_uniprot_ptm_sites.py --fasta inputs/reference/uniprot_human_reviewed_2026_03.fasta \
        --raw-dir PTM_data/marks --out-dir inputs/marks

Reads the "Modified residue" features of every reviewed human UniProtKB/Swiss-Prot entry (UniProt REST API, CC BY 4.0;
the JSON stream of `reviewed:true AND organism_id:9606`, fields accession, ft_mod_res and sequence). With --raw it
reads an earlier download instead of fetching. The API serves the current release only, which may be newer than the
reference proteome: the release served is recorded (`source_release`) and an entry whose sequence differs from the
reference is left out whole, so a position is never moved between sequences.

Marks, terms and states (a term is the feature description before its first ';'):

  methyllysine     N6-methyllysine me1; N6,N6-dimethyllysine me2; N6,N6,N6-trimethyllysine me3;
                   N6-methylated lysine me (degree not stated)                                  residue K
  phosphorylation  Phosphoserine S, Phosphothreonine T, Phosphotyrosine Y: state phospho
  methylarginine   Omega-N-methylarginine me1; Asymmetric dimethylarginine me2a; Symmetric dimethylarginine me2s;
                   Dimethylated arginine me2; Omega-N-methylated arginine me                    residue R
                   (a cross-check only: the methylarginine release has its own table)

A residue carrying several states has one row per state. Keys: `canonical_UniProtAC` (an entry of --fasta) and
`position` on its sequence. Columns of every table:

  canonical_UniProtAC, position, residue, mark, state,
  evidence          the ECO codes of the feature, sorted, ';'-joined
  evidence_class    experimental: ECO:0000269 (published experiment) or ECO:0007744 (combinatorial computational and
                    experimental evidence, i.e. large-scale proteomics) among the codes. by-similarity: every other
                    feature (ECO:0000250 by similarity, and ECO:0000255, 0000303, 0000305 or no code)
  pmids             the PubMed IDs of the feature's evidence, sorted, ';'-joined
  source            UniProt
  source_release    the UniProt release served

Dropped and counted in the summary JSON: not_in_reference (accession absent from --fasta), sequence_differs (entry
sequence is not the reference's), on_an_isoform (feature on an isoform's sequence), not_single_residue (a range or
uncertain position), residue_mismatch (the reference has another residue at the position: `common.residue_is`).
Outputs are sorted, so running again on the same download gives the same bytes. Every download is listed in
inputs/marks/raw_files.tsv and each staged file has its row in inputs/MANIFEST.tsv (redistribution git: CC BY 4.0).
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from datetime import date
from pathlib import Path
from urllib.request import Request, urlopen

import pandas as pd

from common import md5_file, parse_fasta, refresh_manifest_row, residue_is, update_manifest, upsert_tsv_rows

URL = ('https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:9606'
       '&fields=accession,ft_mod_res,sequence&format=json')
CITATION = 'UniProt Consortium, Nucleic Acids Res 53:D609-D617 (2025), doi:10.1093/nar/gkae1010'
LICENSE = 'CC BY 4.0'
# feature description (before ';') -> (mark, state, residue)
TERMS = {
    'N6-methyllysine': ('methyllysine', 'me1', 'K'),
    'N6,N6-dimethyllysine': ('methyllysine', 'me2', 'K'),
    'N6,N6,N6-trimethyllysine': ('methyllysine', 'me3', 'K'),
    'N6-methylated lysine': ('methyllysine', 'me', 'K'),
    'Phosphoserine': ('phosphorylation', 'phospho', 'S'),
    'Phosphothreonine': ('phosphorylation', 'phospho', 'T'),
    'Phosphotyrosine': ('phosphorylation', 'phospho', 'Y'),
    'Omega-N-methylarginine': ('methylarginine', 'me1', 'R'),
    'Asymmetric dimethylarginine': ('methylarginine', 'me2a', 'R'),
    'Symmetric dimethylarginine': ('methylarginine', 'me2s', 'R'),
    'Dimethylated arginine': ('methylarginine', 'me2', 'R'),
    'Omega-N-methylated arginine': ('methylarginine', 'me', 'R'),
}
EXPERIMENTAL = {'ECO:0000269', 'ECO:0007744'}
COLUMNS = ['canonical_UniProtAC', 'position', 'residue', 'mark', 'state', 'evidence', 'evidence_class', 'pmids',
           'source', 'source_release']
OUT_NAME = {'methyllysine': 'uniprot_{rel}_methyllysine_sites.tsv', 'phosphorylation': 'uniprot_{rel}_phosphorylation_sites.tsv',
            'methylarginine': 'uniprot_{rel}_methylarginine_crosscheck.tsv'}
CONTENT = {
    'methyllysine': 'UniProt N6-methyllysine, N6,N6-dimethyllysine and N6,N6,N6-trimethyllysine (and N6-methylated lysine, degree not stated) annotations of the reviewed human entries, one row per residue and state, with ECO evidence codes, their class (experimental or by-similarity) and PubMed IDs',
    'phosphorylation': 'UniProt phosphoserine, phosphothreonine and phosphotyrosine annotations of the reviewed human entries, one row per residue, with ECO evidence codes, their class (experimental or by-similarity) and PubMed IDs',
    'methylarginine': 'UniProt arginine methylation annotations (omega-N-methylarginine, asymmetric and symmetric dimethylarginine, dimethylated and omega-N-methylated arginine) of the reviewed human entries: a cross-check only, the methylarginine release has its own table',
}


def evidence_class(codes) -> str:
    """experimental if any code is ECO:0000269 or ECO:0007744, else by-similarity (ECO:0000250 and every other code)."""
    return 'experimental' if any(c in EXPERIMENTAL for c in codes) else 'by-similarity'


def feature_term(description: str) -> str:
    return description.split(';')[0].strip()


def stage(entries: list[dict], reference: dict[str, str], release: str) -> tuple[dict[str, pd.DataFrame], dict]:
    """The per-mark tables and the counts of what was dropped, from the `results` of the UniProt JSON."""
    rows, dropped = [], Counter()
    for entry in entries:
        acc = entry['primaryAccession']
        sequence = entry.get('sequence', {}).get('value', '')
        for feature in entry.get('features', []):
            if feature.get('type') != 'Modified residue':
                continue
            term = TERMS.get(feature_term(feature.get('description', '')))
            if term is None:
                continue
            mark, state, residue = term
            loc = feature['location']
            if 'sequence' in loc:
                dropped['on_an_isoform'] += 1
                continue
            start, end = loc['start'], loc['end']
            if start['value'] != end['value'] or start['modifier'] != 'EXACT' or end['modifier'] != 'EXACT':
                dropped['not_single_residue'] += 1
                continue
            if acc not in reference:
                dropped['not_in_reference'] += 1
                continue
            if sequence != reference[acc]:
                dropped['sequence_differs'] += 1
                continue
            pos = int(start['value'])
            if not residue_is(reference, acc, pos, residue):
                dropped['residue_mismatch'] += 1
                continue
            evidences = feature.get('evidences', [])
            codes = sorted({e['evidenceCode'] for e in evidences})
            pmids = sorted({e['id'] for e in evidences if e.get('source') == 'PubMed' and e.get('id')}, key=int)
            rows.append({'canonical_UniProtAC': acc, 'position': pos, 'residue': residue, 'mark': mark, 'state': state,
                         'evidence': ';'.join(codes), 'evidence_class': evidence_class(codes), 'pmids': ';'.join(pmids),
                         'source': 'UniProt', 'source_release': release})
    table = pd.DataFrame(rows, columns=COLUMNS)
    # a feature listed twice for the residue and state (two notes) is one row; evidence and PMIDs are merged
    table = (table.groupby(['canonical_UniProtAC', 'position', 'residue', 'mark', 'state', 'source', 'source_release'], as_index=False)
             .agg(evidence=('evidence', lambda v: ';'.join(sorted({c for x in v for c in x.split(';') if c}))),
                  pmids=('pmids', lambda v: ';'.join(sorted({c for x in v for c in x.split(';') if c}, key=int)))))
    table['evidence_class'] = [evidence_class(e.split(';')) for e in table['evidence']]
    table = table[COLUMNS].sort_values(['canonical_UniProtAC', 'position', 'state']).reset_index(drop=True)
    return {m: table[table['mark'] == m].reset_index(drop=True) for m in OUT_NAME}, dict(dropped)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fasta', required=True, help='the reference proteome (canonical reviewed accessions)')
    ap.add_argument('--raw', help='read this download (and its .release file) instead of fetching')
    ap.add_argument('--raw-dir', default='PTM_data/marks', help='where the download is kept (not in Git)')
    ap.add_argument('--out-dir', default='inputs/marks')
    ap.add_argument('--manifest', default='inputs/MANIFEST.tsv')
    ap.add_argument('--downloaded', default='', help='date of the download for raw_files.tsv (default: today when fetching)')
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    if args.raw:
        raw_path = Path(args.raw)
        release = raw_path.with_suffix('.release').read_text().strip()
    else:
        with urlopen(Request(URL, headers={'User-Agent': 'PTM-pipeline-stager/1.0'}), timeout=1800) as response:
            release = response.headers.get('X-UniProt-Release', '')
            if not release:
                raise SystemExit('the UniProt API gave no X-UniProt-Release header')
            data = response.read()
        raw_path = Path(args.raw_dir) / f'uniprot_human_reviewed_mod_res_{release}.json'
        raw_path.parent.mkdir(parents=True, exist_ok=True)
        raw_path.write_bytes(data)
        raw_path.with_suffix('.release').write_text(release + '\n')
    raw_md5 = md5_file(raw_path)
    raw_sha = hashlib.sha256(raw_path.read_bytes()).hexdigest()
    entries = json.loads(raw_path.read_text(encoding='utf-8'))['results']
    reference = parse_fasta(args.fasta)

    tables, dropped = stage(entries, reference, release)
    out_dir.mkdir(parents=True, exist_ok=True)
    files = {}
    for mark, table in tables.items():
        name = OUT_NAME[mark].format(rel=release)
        table.to_csv(out_dir / name, sep='\t', index=False, lineterminator='\n')
        files[mark] = name
    summary = {
        'source': 'UniProtKB/Swiss-Prot reviewed human, "Modified residue" features (UniProt REST API)',
        'url': URL,
        'uniprot_release_served': release,
        'reference_fasta': args.fasta,
        'reference_fasta_md5': md5_file(args.fasta),
        'raw_file': raw_path.name,
        'raw_md5': raw_md5,
        'entries_in_download': len(entries),
        'dropped': dropped,
        'rows': {m: len(t) for m, t in tables.items()},
        'rows_by_state': {m: {k: int(v) for k, v in t['state'].value_counts().sort_index().items()} for m, t in tables.items()},
        'rows_by_evidence_class': {m: {k: int(v) for k, v in t['evidence_class'].value_counts().sort_index().items()}
                                   for m, t in tables.items()},
        'files': files,
        'files_md5': {m: md5_file(out_dir / f) for m, f in files.items()},
    }
    summary_name = 'uniprot_ptm_sites_summary.json'
    (out_dir / summary_name).write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n', encoding='utf-8')

    # raw_files.tsv: a row that is already there for the same bytes keeps its download date
    raw_list = out_dir / 'raw_files.tsv'
    head = 'file\turl\tbytes\tmd5\tsha256\tdownloaded\tversion\tkept'
    if not raw_list.exists():
        raw_list.write_text(head + '\n', encoding='utf-8')
    old = pd.read_csv(raw_list, sep='\t', dtype=str, keep_default_na=False)
    old = old[(old['file'] == raw_path.name) & (old['md5'] == raw_md5)]
    when = old['downloaded'].iloc[0] if len(old) else (args.downloaded or date.today().isoformat())
    upsert_tsv_rows(raw_list, [{'file': raw_path.name, 'url': URL, 'bytes': raw_path.stat().st_size, 'md5': raw_md5,
                                'sha256': raw_sha, 'downloaded': when, 'version': f'UniProt release {release} (X-UniProt-Release header)',
                                'kept': 'checksum_only (30 MB; re-download is the source of the tables, which are in Git)'}])

    command = (f'python scripts/54_stage_uniprot_ptm_sites.py --fasta {args.fasta} '
               f'--raw-dir PTM_data/marks --out-dir inputs/marks')
    rows = []
    for mark, name in files.items():
        path = out_dir / name
        rows.append({'file': f'marks/{name}', 'bytes': path.stat().st_size, 'md5': md5_file(path), 'rows': len(tables[mark]),
                     'redistribution': 'git', 'content': CONTENT[mark] + f': {len(tables[mark]):,} rows',
                     'source': URL, 'version_or_date': f'UniProt release {release} (X-UniProt-Release header); downloaded {when}',
                     'citation': CITATION, 'licence': LICENSE, 'derived_from': f'{raw_path.name} (MD5 {raw_md5})',
                     'how_made_or_checked': command + '; residue checked on the reference sequence, entries whose sequence differs dropped; counts in marks/uniprot_ptm_sites_summary.json'})
    path = out_dir / summary_name
    rows.append({'file': f'marks/{summary_name}', 'bytes': path.stat().st_size, 'md5': md5_file(path),
                 'redistribution': 'git', 'content': 'Counts, release served, drops and the MD5 of the input and of every table of the UniProt PTM-site staging',
                 'source': 'the tables above', 'version_or_date': f'UniProt release {release}', 'derived_from': 'marks/raw_files.tsv',
                 'how_made_or_checked': command})
    update_manifest(args.manifest, rows)
    # raw_files.tsv has its own row, refreshed last (it holds the MD5 of every raw download)
    refresh_manifest_row(args.manifest, raw_list, 'marks/raw_files.tsv', redistribution='git',
                         content='URL, size, MD5, SHA-256 and date of every raw download behind inputs/marks/',
                         source='the URLs it lists', derived_from='', how_made_or_checked='written by the scripts/54 and scripts/55 staging runs')
    print(json.dumps({k: summary[k] for k in ('uniprot_release_served', 'dropped', 'rows', 'rows_by_state', 'rows_by_evidence_class')}, indent=1))


if __name__ == '__main__':
    main()

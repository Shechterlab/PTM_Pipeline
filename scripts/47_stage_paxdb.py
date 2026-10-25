#!/usr/bin/env python3
"""Stage PaxDb 6.1 protein abundance (Homo sapiens, whole organism, integrated) for the reviewed canonical proteome.

PaxDb (pax-db.org; von Mering lab, SIB / University of Zurich) gives protein abundance in parts per million (ppm) of the
proteome, one value per STRING protein (STRING v12.0 identifiers for PaxDb 6.x). Its integrated datasets are weighted
averages of the datasets of one organism or tissue; this script stages the human whole-organism one
(9606-WHOLE_ORGANISM-integrated.txt), to be used as a detection covariate (abundant proteins are the ones mass
spectrometry finds most). Licence: CC BY 4.0 (https://pax-db.org/license). Cite Huang et al., Mol Cell Proteomics
22:100640 (2023), doi:10.1016/j.mcpro.2023.100640 (PaxDb 5.0, the publication PaxDb's licence page links for credit).

Downloads read (--paxdb folder)
  9606.zip                           every human dataset of PaxDb 6.1 (https://pax-db.org/downloads/6.1/datasets/9606.zip)
  paxdb-uniprot-crossreferences.txt  PaxDb's UniProt accession -> STRING protein table, all species
  paxdb-uniprot-links-v6.1.zip       PaxDb's STRING protein -> UniProt entry name table (optional; only compared, not used)

Mapping to canonical_UniProtAC (one value per reviewed entry of --fasta)
  1. PaxDb's own cross-reference: the entry's accession -> STRING protein (paxdb-uniprot-crossreferences.txt).
  2. Entries step 1 does not reach: UniProt's own STRING cross-reference of the entry (release 2026_03,
     --uniprot-string, the table scripts/37a staged), when no other entry took that STRING protein in step 1.
  An entry given several STRING proteins with a value takes the largest (none in PaxDb 6.1); a STRING protein given to
  several entries gives its value to each (none in 6.1). Entry names (paxdb-uniprot-links) are not used: UniProt
  reassigns them between releases; the summary counts how often they agree. Isoforms are not involved: PaxDb
  cross-references entry accessions only. An entry encoded by several genes (identical histone or haemoglobin
  copies, say) takes the value of the one STRING protein the cross-reference names, not the sum of its copies (PaxDb
  lists each copy as a protein of its own; the summary counts the copies and whether they carry the same value).

The percentile of an entry is 100 x its average rank (ascending abundance, ties averaged) / the number of staged
entries: the lowest has 100/n, the highest 100.

Outputs (in --outdir)
  paxdb_human_integrated_abundance.tsv   one row per reviewed entry with a value: canonical_UniProtAC, gene (UniProt
                                         primary gene names), paxdb_gene_name, string_id, abundance_ppm (as PaxDb
                                         writes it), percentile, mapped_by (paxdb_xref or uniprot_xref)
  paxdb_summary.json                     versions, input MD5s, dataset header, the per-tissue integrated datasets of
                                         the same download (listed, not staged), counts, what was left out and why
  manifest_rows.tsv                      the rows these two files would take in inputs/MANIFEST.tsv (its columns)
  raw_files_rows.tsv                     the rows of the downloads for inputs/layers/raw_files.tsv (its columns)
Nothing outside --outdir is written: --manifest and --raw-files are read only, for their columns.

    python scripts/47_stage_paxdb.py --paxdb <download folder> --outdir <folder>
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import zipfile
from collections import Counter
from pathlib import Path

import pandas as pd

from common import md5_file, parse_fasta, save_table

VERSION = '6.1'
STRING_VERSION = '12.0'
TAXON = '9606'
DATASET = f'{TAXON}/{TAXON}-WHOLE_ORGANISM-integrated.txt'
COLUMNS = ['gene_name', 'string_external_id', 'abundance']      # the column header PaxDb 6.1 writes (after '#')
BASE = f'https://pax-db.org/downloads/{VERSION}'
# The downloads: URL and the Last-Modified the server gave (2026-10-25)
DOWNLOADS = {
    '9606.zip': (f'{BASE}/datasets/9606.zip', '2026-07-09 11:19:44 GMT'),
    'paxdb-uniprot-crossreferences.txt': (f'{BASE}/paxdb-uniprot-crossreferences.txt', '2026-07-09 11:18:14 GMT'),
    'paxdb-uniprot-links-v6.1.zip': (f'{BASE}/paxdb-uniprot-links-v6.1.zip', '2026-07-09 11:20:28 GMT'),
}
LICENCE = ('CC BY 4.0 (https://pax-db.org/license: "Creative Commons Attribution 4.0 International (CC BY 4.0) applies to all '
           'copyrightable parts of our database")')
MANIFEST_LICENCE = 'CC BY 4.0 (https://pax-db.org/license)'
CITATION = ('Huang et al., Mol Cell Proteomics 22:100640 (2023), doi:10.1016/j.mcpro.2023.100640 (PaxDb 5.0); PaxDb 6.1, '
            'https://pax-db.org')
OUT_TABLE = 'paxdb_human_integrated_abundance.tsv'
OUT_SUMMARY = 'paxdb_summary.json'


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read_dataset(text: str, name: str) -> tuple[dict, pd.DataFrame]:
    """A PaxDb dataset file: its '#key: value' header and its rows (gene_name, string_id, abundance as written)."""
    meta, rows, head = {}, [], None
    for line in text.splitlines():
        if not line.strip():
            continue
        if line.startswith('#'):
            body = line[1:]
            if '\t' in body:
                head = body.split('\t')
            elif ':' in body:
                k, v = body.split(':', 1)
                meta[k.strip()] = v.strip()
            continue
        rows.append(line.split('\t'))
    if head != COLUMNS:
        raise SystemExit(f'{name}: column header {head}, expected {COLUMNS}')
    bad = [r for r in rows if len(r) != 3]
    if bad:
        raise SystemExit(f'{name}: {len(bad)} rows without 3 fields, e.g. {bad[0]}')
    table = pd.DataFrame(rows, columns=['paxdb_gene_name', 'string_id', 'abundance_ppm'])
    table['value'] = table['abundance_ppm'].astype(float)
    if table['string_id'].duplicated().any():
        raise SystemExit(f'{name}: a STRING protein is listed twice')
    return meta, table


def write_rows(path: Path, columns: list[str], rows: list[dict]) -> None:
    """A TSV with the given columns, written as common.upsert_tsv_rows writes rows (no quoting)."""
    lines = ['\t'.join(columns)]
    for r in rows:
        cells = [str(r.get(c, '')) for c in columns]
        if any(ch in c for c in cells for ch in '\t\n\r'):
            raise SystemExit(f'{path.name}: a value holds a tab or line break')
        lines.append('\t'.join(cells))
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def weights(text: str) -> dict[str, float]:
    return {k: float(v) for k, v in (w.rsplit(':', 1) for w in text.split(';') if w.strip())}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--paxdb', required=True, help='folder with 9606.zip and paxdb-uniprot-crossreferences.txt (and optionally '
                    'paxdb-uniprot-links-v6.1.zip)')
    ap.add_argument('--fasta', default='inputs/reference/uniprot_human_reviewed_2026_03.fasta')
    ap.add_argument('--entries', default='inputs/reference/uniprot_human_reviewed_2026_03_entries.tsv')
    ap.add_argument('--uniprot-string', default='inputs/interactions/uniprot_human_reviewed_2026-10-16.tsv',
                    help="UniProt 2026_03 reviewed entries with UniProt's STRING cross-reference (column string_id)")
    ap.add_argument('--manifest', default='inputs/MANIFEST.tsv', help='read only: its columns')
    ap.add_argument('--raw-files', default='inputs/layers/raw_files.tsv', help='read only: its columns')
    ap.add_argument('--downloaded', default='2026-10-25', help='the day the downloads were made')
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    src, out = Path(args.paxdb), Path(args.outdir)
    canonical = parse_fasta(args.fasta)
    entries = pd.read_csv(args.entries, sep='\t', dtype=str, keep_default_na=False).set_index('Entry')
    entries = entries[entries.index.isin(canonical)]
    genes_of = {a: [g.strip() for g in n.split(';') if g.strip()] for a, n in entries['Gene Names (primary)'].items()}

    # the whole-organism dataset, and the other integrated datasets of the same download (listed only)
    with zipfile.ZipFile(src / '9606.zip') as zf:
        members = sorted(i.filename for i in zf.infolist() if not i.is_dir())
        info = zf.getinfo(DATASET)
        meta, data = read_dataset(zf.read(DATASET).decode('utf-8'), DATASET)
        tissues = []
        for name in members:
            if name.endswith('-integrated.txt') and name != DATASET:
                m, t = read_dataset(zf.read(name).decode('utf-8'), name)
                tissues.append({'file': name.split('/')[-1], 'organ': m.get('organ', ''), 'name': m.get('name', ''),
                                'score': m.get('score', ''), 'coverage': m.get('coverage', ''), 'proteins': int(len(t)),
                                'weighted_datasets': len(weights(m.get('weights', '')))})
    value = dict(zip(data['string_id'], data['value']))

    # 1. PaxDb's cross-reference
    xref = pd.read_csv(src / 'paxdb-uniprot-crossreferences.txt', sep='\t', header=None, names=['accession', 'paxdb_id'], dtype=str,
                       keep_default_na=False)
    xref = xref[xref['paxdb_id'].str.startswith(TAXON + '-')]
    xref = xref.assign(string_id=TAXON + '.' + xref['paxdb_id'].str[len(TAXON) + 1:])
    if xref['accession'].str.contains('-').any():
        raise SystemExit('the PaxDb cross-reference holds isoform accessions: decide how to map them')
    ref_xref = xref[xref['accession'].isin(canonical)]
    step1 = ref_xref.groupby('accession')['string_id'].agg(lambda s: sorted(set(s))).to_dict()
    taken = {s for ids in step1.values() for s in ids}
    # 2. UniProt's own STRING cross-reference, for the entries step 1 does not reach
    up = pd.read_csv(args.uniprot_string, sep='\t', dtype=str, keep_default_na=False)
    up = up[up['uniprot_accession'].isin(canonical) & (up['string_id'] != '')]
    up_ids = {a: sorted({s.strip() for s in v.split(';') if s.strip()}) for a, v in zip(up['uniprot_accession'], up['string_id'])}
    agree = sum(1 for a, ids in step1.items() if a in up_ids and set(ids) & set(up_ids[a]))
    disagree = sum(1 for a, ids in step1.items() if a in up_ids and not set(ids) & set(up_ids[a]))
    step2, step2_taken = {}, 0
    for a, ids in sorted(up_ids.items()):
        if a in step1:
            continue
        free = [s for s in ids if s not in taken]
        step2_taken += len(ids) - len(free)
        if free:
            step2[a] = free

    rows, left = [], Counter()
    several, shared = 0, Counter()
    for acc in sorted(canonical):
        how, ids = ('paxdb_xref', step1[acc]) if acc in step1 else ('uniprot_xref', step2[acc]) if acc in step2 else ('', [])
        with_value = [s for s in ids if s in value]
        if not ids:
            left['entries_not_cross_referenced_to_string'] += 1
            continue
        if not with_value:
            left[f'entries_{how}_string_protein_not_in_dataset'] += 1
            continue
        if len(with_value) > 1:
            several += 1
        best = sorted(with_value, key=lambda s: (-value[s], s))[0]
        shared[best] += 1
        rows.append({'canonical_UniProtAC': acc, 'gene': '; '.join(genes_of.get(acc, [])), 'string_id': best, 'mapped_by': how})
    table = pd.DataFrame(rows).merge(data[['string_id', 'paxdb_gene_name', 'abundance_ppm', 'value']], on='string_id', how='left')
    pct = table['value'].rank(method='average', pct=True) * 100
    table['percentile'] = [f'{p:.2f}' for p in pct]
    table = table[['canonical_UniProtAC', 'gene', 'paxdb_gene_name', 'string_id', 'abundance_ppm', 'percentile', 'mapped_by']]
    table = table.sort_values('canonical_UniProtAC').reset_index(drop=True)

    # what of the dataset is not used, and why
    used = set(table['string_id'])
    by_id = xref.groupby('string_id')['accession'].agg(list).to_dict()
    unused = Counter()
    for s in data['string_id']:
        if s in used:
            continue
        accs = by_id.get(s, [])
        unused['not_in_paxdb_cross_reference' if not accs else 'cross_referenced_to_no_reference_entry'
               if not any(a in canonical for a in accs) else 'other'] += 1
    # entries encoded by several genes: dataset proteins of their other gene copies (by PaxDb gene name), not used
    gene_entry = {}
    for acc in table['canonical_UniProtAC']:
        for g in genes_of.get(acc, []):
            gene_entry.setdefault(g, set()).add(acc)
    staged_value = dict(zip(table['canonical_UniProtAC'], table['abundance_ppm'].astype(float)))
    copies = data[~data['string_id'].isin(used) & data['paxdb_gene_name'].isin(gene_entry)]
    same = sum(1 for g, v in zip(copies['paxdb_gene_name'], copies['value'])
               if any(abs(staged_value[a] - v) <= 1e-9 * max(v, 1e-12) for a in gene_entry[g]))
    copy_entries = sorted({a for g in copies['paxdb_gene_name'] for a in gene_entry[g]})
    gene_ok = sum(1 for a, g in zip(table['canonical_UniProtAC'], table['paxdb_gene_name']) if g in genes_of.get(a, []))

    # PaxDb's entry-name table, compared (not used)
    links = {}
    lp = src / 'paxdb-uniprot-links-v6.1.zip'
    if lp.exists():
        with zipfile.ZipFile(lp) as zf:
            name = [n for n in zf.namelist() if n.endswith('.tsv')][0]
            lk = pd.read_csv(io.BytesIO(zf.read(name)), sep='\t', header=None, names=['string_id', 'entry_name'], dtype=str,
                             keep_default_na=False)
        lk = lk[lk['string_id'].str.startswith(TAXON + '.')]
        acc_of_name = {n: a for a, n in entries['Entry Name'].items()}
        lk = lk.assign(acc=lk['entry_name'].map(acc_of_name))
        mine = dict(zip(table['canonical_UniProtAC'], table['string_id']))
        hit = lk.dropna(subset=['acc'])
        hit_ds = hit[hit['string_id'].isin(value)]
        links = {'human_rows': int(len(lk)), 'entry_name_is_a_reference_entry_name': int(len(hit)),
                 'dataset_proteins_named': int(len(hit_ds)),
                 'agree_with_staged': int(sum(mine.get(a) == s for a, s in zip(hit_ds['acc'], hit_ds['string_id']))),
                 'differ_from_staged': int(sum(a in mine and mine[a] != s for a, s in zip(hit_ds['acc'], hit_ds['string_id']))),
                 'name_of_an_entry_not_staged': int(sum(a not in mine for a in hit_ds['acc'])),
                 'input_md5': md5_file(lp)}

    save_table(table, out / OUT_TABLE)
    v = table['abundance_ppm'].astype(float)
    summary = {
        'paxdb': {'version': VERSION, 'string_version': STRING_VERSION, 'version_source': 'https://pax-db.org (site config '
                  'PAXDB_VERSION 6.1; "What\'s new": v6.0 -> v6.1 (May 2026)); https://pax-db.org/downloads/latest/ lists the 6.1 files',
                  'licence': LICENCE, 'citation': CITATION, 'downloaded': args.downloaded},
        'inputs': {n: {'url': u, 'last_modified': lm, 'md5': md5_file(src / n)} for n, (u, lm) in DOWNLOADS.items() if (src / n).exists()},
        'reference': {'fasta_md5': md5_file(args.fasta), 'entries_md5': md5_file(args.entries),
                      'uniprot_string_md5': md5_file(args.uniprot_string), 'reference_entries': len(canonical)},
        'dataset': {'member': DATASET, 'member_date': '%04d-%02d-%02d' % info.date_time[:3], 'member_bytes': info.file_size,
                    'member_md5': hashlib.md5(zipfile.ZipFile(src / '9606.zip').read(DATASET)).hexdigest(),
                    **{k: meta.get(k, '') for k in ('id', 'name', 'score', 'coverage', 'organ', 'integrated', 'publication_year')},
                    'weights': weights(meta.get('weights', '')), 'proteins': int(len(data)),
                    'ppm_total': round(float(data['value'].sum()), 1)},
        'download_members': len(members),
        'other_integrated_datasets': {'count': len(tissues), 'staged': False,
                                      'why': 'each covers a subset of the proteome and the set would be ~50 columns or ~0.5 M rows; '
                                             'only the whole-organism dataset is staged', 'datasets': tissues},
        'mapping': {
            'paxdb_cross_reference_human_rows': int(len(xref)), 'paxdb_cross_reference_reference_entries': len(step1),
            'paxdb_vs_uniprot_cross_reference': {'agree': agree, 'disagree': disagree,
                                                 'paxdb_only': sum(1 for a in step1 if a not in up_ids),
                                                 'uniprot_only': sum(1 for a in up_ids if a not in step1)},
            'uniprot_cross_reference_entries_added': len(step2),
            'uniprot_cross_reference_string_ids_already_taken': step2_taken,
            'entries_with_several_string_proteins_with_value': several,
            'string_proteins_given_to_several_entries': sum(1 for n in shared.values() if n > 1),
            'staged_entries_whose_paxdb_gene_name_is_a_uniprot_primary_gene': gene_ok,
            'paxdb_uniprot_links_entry_names': links,
        },
        'staged_entries': int(len(table)),
        'staged_by': dict(Counter(table['mapped_by'])),
        'entries_without_value': len(canonical) - int(len(table)),
        'entries_without_value_by_reason': dict(sorted(left.items())),
        'dataset_proteins_not_staged': dict(sorted(unused.items())),
        'unstaged_proteins_named_like_a_staged_entry_gene': {
            'dataset_proteins': int(len(copies)), 'entries': len(copy_entries), 'same_value_as_the_entry': same,
            'note': 'dataset proteins left out whose PaxDb gene name is a primary gene of a staged entry, mostly further gene '
                    'copies of an entry (e.g. the H4C genes of P62805, HBA1 of P69905): the entry keeps the value of its '
                    'cross-referenced protein, not the sum of its copies'},
        'abundance_ppm_quantiles_staged': {q: float(f'{v.quantile(x):.6g}') for q, x in
                                           (('min', 0), ('q25', .25), ('median', .5), ('q75', .75), ('max', 1))},
    }
    (out / OUT_SUMMARY).write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')

    # rows for inputs/MANIFEST.tsv and inputs/layers/raw_files.tsv (written here, not there)
    mcols = open(args.manifest, encoding='utf-8').readline().rstrip('\n').split('\t')
    rcols = open(args.raw_files, encoding='utf-8').readline().rstrip('\n').split('\t')
    source = f'{DOWNLOADS["9606.zip"][0]} (member {DATASET}); {DOWNLOADS["paxdb-uniprot-crossreferences.txt"][0]}'
    version = (f'PaxDb {VERSION} (May 2026; STRING v{STRING_VERSION} protein identifiers; dataset id {meta.get("id", "")}, files dated '
               f'2026-07-09), downloaded {args.downloaded}')
    how = ('scripts/47_stage_paxdb.py: PaxDb\'s UniProt cross-reference, else UniProt 2026_03\'s STRING cross-reference '
           '(interactions/uniprot_human_reviewed_2026-10-16.tsv); one STRING protein per entry; checked against '
           'reference/uniprot_human_reviewed_2026_03.fasta')
    mrows = [
        {'file': f'layers/{OUT_TABLE}', 'content': 'PaxDb 6.1 protein abundance (ppm) of reviewed canonical human proteins in the '
         'whole-organism integrated dataset (H. sapiens): one row per entry with a value, with its STRING protein, the PaxDb and '
         'UniProt gene names, the percentile among the staged entries and which cross-reference mapped it',
         'derived_from': '9606.zip (9606-WHOLE_ORGANISM-integrated.txt); paxdb-uniprot-crossreferences.txt; '
                         'interactions/uniprot_human_reviewed_2026-10-16.tsv', 'how_made_or_checked': how},
        {'file': f'layers/{OUT_SUMMARY}', 'content': 'Versions, input MD5s, dataset header and weights, the other integrated '
         'datasets of the download, mapping counts and what was left out of the PaxDb staging', 'derived_from': '',
         'how_made_or_checked': 'scripts/47_stage_paxdb.py'},
    ]
    for r in mrows:
        p = out / r['file'].split('/')[-1]
        r.update({'bytes': p.stat().st_size, 'md5': md5_file(p), 'rows': '' if p.suffix == '.json' else len(table),
                  'redistribution': 'git', 'source': source, 'version_or_date': version, 'citation': CITATION, 'licence': MANIFEST_LICENCE})
    write_rows(out / 'manifest_rows.tsv', mcols, mrows)
    rrows = []
    for n, (u, lm) in DOWNLOADS.items():
        p = src / n
        if not p.exists():
            continue
        what = {'9606.zip': f'PaxDb {VERSION}, all human datasets (server date {lm})',
                'paxdb-uniprot-crossreferences.txt': f'PaxDb {VERSION} UniProt cross-reference, all species (server date {lm})',
                'paxdb-uniprot-links-v6.1.zip': f'PaxDb {VERSION} UniProt entry-name links, all species (server date {lm}); compared only'}[n]
        rrows.append({'file': n, 'url': u, 'bytes': p.stat().st_size, 'md5': md5_file(p), 'sha256': sha256_file(p),
                      'downloaded': args.downloaded, 'version': what, 'kept': 'release_asset'})
    write_rows(out / 'raw_files_rows.tsv', rcols, rrows)
    print(json.dumps({k: summary[k] for k in ('staged_entries', 'staged_by', 'entries_without_value', 'entries_without_value_by_reason',
                                              'dataset_proteins_not_staged', 'unstaged_proteins_named_like_a_staged_entry_gene')}, indent=2))


if __name__ == '__main__':
    main()

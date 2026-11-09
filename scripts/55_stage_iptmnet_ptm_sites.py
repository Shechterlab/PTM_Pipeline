#!/usr/bin/env python3
"""Stage iPTMnet's human phosphorylation and lysine and arginine methylation records as one table per mark (inputs/marks/).

    python scripts/55_stage_iptmnet_ptm_sites.py --fasta inputs/reference/uniprot_human_reviewed_2026_03.fasta \
        --raw-dir PTM_data/marks --out-dir inputs/marks

Reads ptm.txt (and readme.txt, for the release) of https://research.bioinformatics.udel.edu/iptmnet_data/files/current/
(iPTMnet, Protein Information Resource, CC BY-NC-SA 4.0 per its license page and its download readme; the readme's
last line also mentions CC BY 4.0, and this repository follows the stricter terms). With --raw it reads a kept download.

ptm.txt has one record per upstream database and site, without a header: ptm_type, source, substrate_UniProtAC,
substrate_genename, organism, site (residue letter and position), enzyme_UniProtAC, enzyme_genename, note, pmid. Kept:
organism Homo sapiens; PHOSPHORYLATION at S, T or Y (mark phosphorylation, state phospho); METHYLATION at K (methyllysine)
or R (methylarginine, a cross-check only). iPTMnet gives no methylation degree: state is `me`.

Keys: `canonical_UniProtAC` of --fasta and `position` on its sequence. An accession with the isoform suffix -1 is the
canonical entry. Dropped and counted in the summary JSON: not_human; other_ptm_type; other_residue (a methylation at C, Q,
E, ... or a phosphorylation not at S, T or Y); on_an_isoform (suffix -2 and up); not_in_reference; residue_mismatch (the
reference has another residue at the position: `common.residue_is`); unparsed_site.

Records of the same residue are merged. Columns: canonical_UniProtAC, position, residue, mark, state, evidence (empty:
iPTMnet has no ECO codes), evidence_class, pmids, source=iPTMnet, source_release, upstream_sources.

  upstream_sources  the upstream databases of iPTMnet's `source` column that report the residue, sorted, ';'-joined
                    (hprd, unip = UniProt, rlim and rlim+ = RLIMS-P text mining, sign = SIGNOR, npro = neXtProt,
                    iedb, pro = Protein Ontology, intact); iPTMnet asks that they are cited when it is redistributed
  pmids             the PubMed IDs of those records
  evidence_class    the kind of record, not a new judgment of its evidence. experimental: a curated record of
                    hprd, sign, npro, pro, intact or iedb (curated from an experiment, with a PubMed ID). text-mined:
                    only RLIMS-P (rlim, rlim+) reports it. unspecified: only unip reports it (iPTMnet's copy of
                    UniProt carries no ECO code; the UniProt table has the codes). A residue several kinds report
                    takes experimental, then text-mined, then unspecified.

Outputs are sorted, so a second run on the same download gives the same bytes. Downloads are listed in
inputs/marks/raw_files.tsv; each table has its MANIFEST row (CC BY-NC-SA 4.0, redistribution git).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter
from datetime import date
from pathlib import Path
from urllib.request import Request, urlopen

import pandas as pd

from common import md5_file, parse_fasta, refresh_manifest_row, residue_is, update_manifest, upsert_tsv_rows

BASE = 'https://research.bioinformatics.udel.edu/iptmnet_data/files/current/'
CITATION = ('Huang H, Arighi CN, Ross KE, et al. iPTMnet: an integrated resource for protein post-translational modification '
            'network discovery. Nucleic Acids Res 46:D542-D550 (2018), doi:10.1093/nar/gkx1104')
LICENSE = 'CC BY-NC-SA 4.0 (iPTMnet license page and download readme; the readme also mentions CC BY 4.0, the stricter terms are followed)'
COLUMNS = ['canonical_UniProtAC', 'position', 'residue', 'mark', 'state', 'evidence', 'evidence_class', 'pmids', 'source',
           'source_release', 'upstream_sources']
FIELDS = ['ptm_type', 'source', 'acc', 'gene', 'organism', 'site', 'enzyme_acc', 'enzyme_gene', 'note', 'pmid']
HUMAN = 'Homo sapiens (Human)'
CURATED = {'hprd', 'sign', 'npro', 'pro', 'intact', 'iedb'}
TEXT_MINED = {'rlim', 'rlim+'}
OUT_NAME = {'methyllysine': 'iptmnet_{rel}_methyllysine_sites.tsv', 'phosphorylation': 'iptmnet_{rel}_phosphorylation_sites.tsv',
            'methylarginine': 'iptmnet_{rel}_methylarginine_crosscheck.tsv'}
CONTENT = {
    'methyllysine': 'iPTMnet human lysine methylation records (degree not stated), one row per residue, with the upstream database and PubMed IDs of each',
    'phosphorylation': 'iPTMnet human phosphoserine, phosphothreonine and phosphotyrosine records, one row per residue, with the upstream database and PubMed IDs of each',
    'methylarginine': 'iPTMnet human arginine methylation records (degree not stated): a cross-check only, the methylarginine release has its own table',
}


def evidence_class(upstream) -> str:
    upstream = set(upstream)
    if upstream & CURATED:
        return 'experimental'
    if upstream & TEXT_MINED:
        return 'text-mined'
    return 'unspecified'


def mark_of(ptm_type: str, letter: str):
    """(mark, state) of a record, or None when it is not a mark staged here."""
    if ptm_type == 'PHOSPHORYLATION' and letter in 'STY':
        return 'phosphorylation', 'phospho'
    if ptm_type == 'METHYLATION' and letter == 'K':
        return 'methyllysine', 'me'
    if ptm_type == 'METHYLATION' and letter == 'R':
        return 'methylarginine', 'me'
    return None


def stage(records: pd.DataFrame, reference: dict[str, str], release: str) -> tuple[dict[str, pd.DataFrame], dict]:
    """The per-mark tables and the counts of what was dropped, from the rows of ptm.txt."""
    rows, dropped = [], Counter()
    for r in records.itertuples(index=False):
        if r.organism != HUMAN:
            dropped['not_human'] += 1
            continue
        if r.ptm_type not in ('PHOSPHORYLATION', 'METHYLATION'):
            dropped['other_ptm_type'] += 1
            continue
        m = re.fullmatch(r'([A-Z])(\d+)', r.site)
        if not m:
            dropped['unparsed_site'] += 1
            continue
        letter, pos = m.group(1), int(m.group(2))
        mark = mark_of(r.ptm_type, letter)
        if mark is None:
            dropped['other_residue'] += 1
            continue
        acc = r.acc[:-2] if r.acc.endswith('-1') else r.acc
        if '-' in acc:
            dropped['on_an_isoform'] += 1
            continue
        if acc not in reference:
            dropped['not_in_reference'] += 1
            continue
        if not residue_is(reference, acc, pos, letter):
            dropped['residue_mismatch'] += 1
            continue
        pmids = [p for p in re.split(r'[;,]', r.pmid) if p.strip().isdigit()]
        rows.append((acc, pos, letter, mark[0], mark[1], r.source, ';'.join(p.strip() for p in pmids)))
    frame = pd.DataFrame(rows, columns=['canonical_UniProtAC', 'position', 'residue', 'mark', 'state', 'upstream', 'pmids'])
    merged = []
    for key, grp in frame.groupby(['canonical_UniProtAC', 'position', 'residue', 'mark', 'state'], sort=True):
        upstream = sorted(set(grp['upstream']))
        pmids = sorted({p for v in grp['pmids'] for p in v.split(';') if p}, key=int)
        merged.append({'canonical_UniProtAC': key[0], 'position': key[1], 'residue': key[2], 'mark': key[3], 'state': key[4],
                       'evidence': '', 'evidence_class': evidence_class(upstream), 'pmids': ';'.join(pmids),
                       'source': 'iPTMnet', 'source_release': release, 'upstream_sources': ';'.join(upstream)})
    table = pd.DataFrame(merged, columns=COLUMNS)
    return {m: table[table['mark'] == m].reset_index(drop=True) for m in OUT_NAME}, dict(dropped)


def fetch(url: str) -> tuple[bytes, str]:
    with urlopen(Request(url, headers={'User-Agent': 'PTM-pipeline-stager/1.0'}), timeout=1800) as response:
        return response.read(), response.headers.get('Last-Modified', '')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--raw-dir', default='PTM_data/marks', help='where the downloads are kept (not in Git)')
    ap.add_argument('--raw', help='read this ptm.txt (with iptmnet_readme.txt beside it) instead of fetching')
    ap.add_argument('--out-dir', default='inputs/marks')
    ap.add_argument('--manifest', default='inputs/MANIFEST.tsv')
    ap.add_argument('--downloaded', default='', help='date of the download for raw_files.tsv (default: today when fetching)')
    args = ap.parse_args()

    out_dir, raw_dir = Path(args.out_dir), Path(args.raw_dir)
    if args.raw:
        ptm_path = Path(args.raw)
        readme_path = ptm_path.with_name('iptmnet_readme.txt')
        modified = {}
    else:
        raw_dir.mkdir(parents=True, exist_ok=True)
        ptm_path, readme_path = raw_dir / 'iptmnet_ptm.txt', raw_dir / 'iptmnet_readme.txt'
        data, modified['ptm.txt'] = fetch(BASE + 'ptm.txt')
        ptm_path.write_bytes(data)
        data, modified['readme.txt'] = fetch(BASE + 'readme.txt')
        readme_path.write_bytes(data)
    readme = readme_path.read_text(encoding='utf-8')
    found = re.search(r'for release ([0-9][0-9.]*)', readme)
    if not found:
        raise SystemExit('the iPTMnet readme names no release')
    release = found.group(1)

    records = pd.read_csv(ptm_path, sep='\t', header=None, names=FIELDS, dtype=str, keep_default_na=False)
    reference = parse_fasta(args.fasta)
    tables, dropped = stage(records, reference, release)
    out_dir.mkdir(parents=True, exist_ok=True)
    files = {}
    for mark, table in tables.items():
        name = OUT_NAME[mark].format(rel=release.replace('.', '_'))
        table.to_csv(out_dir / name, sep='\t', index=False, lineterminator='\n')
        files[mark] = name
    ptm_md5, readme_md5 = md5_file(ptm_path), md5_file(readme_path)
    humans = records[records['organism'] == HUMAN]
    summary = {
        'source': 'iPTMnet ptm.txt (Protein Information Resource)', 'url': BASE + 'ptm.txt', 'iptmnet_release': release,
        'license': LICENSE, 'reference_fasta': args.fasta, 'reference_fasta_md5': md5_file(args.fasta),
        'raw_files': {'ptm.txt': ptm_path.name, 'readme.txt': readme_path.name}, 'raw_md5': {'ptm.txt': ptm_md5, 'readme.txt': readme_md5},
        'records_in_download': len(records), 'human_records': len(humans), 'dropped': dropped,
        'rows': {m: len(t) for m, t in tables.items()},
        'rows_by_state': {m: {k: int(v) for k, v in t['state'].value_counts().sort_index().items()} for m, t in tables.items()},
        'rows_by_evidence_class': {m: {k: int(v) for k, v in t['evidence_class'].value_counts().sort_index().items()} for m, t in tables.items()},
        'rows_by_upstream_source': {m: {k: int(v) for k, v in Counter(u for x in t['upstream_sources'] for u in x.split(';')).most_common()}
                                    for m, t in tables.items()},
        'files': files, 'files_md5': {m: md5_file(out_dir / f) for m, f in files.items()},
    }
    summary_name = 'iptmnet_ptm_sites_summary.json'
    (out_dir / summary_name).write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n', encoding='utf-8')

    raw_list = out_dir / 'raw_files.tsv'
    old = pd.read_csv(raw_list, sep='\t', dtype=str, keep_default_na=False)
    raw_rows = []
    for name, path, md5 in (('ptm.txt', ptm_path, ptm_md5), ('readme.txt', readme_path, readme_md5)):
        same = old[(old['file'] == path.name) & (old['md5'] == md5)]
        when = same['downloaded'].iloc[0] if len(same) else (args.downloaded or date.today().isoformat())
        stamp = f'; Last-Modified {modified[name]}' if modified.get(name) else ''
        version = f'iPTMnet release {release}' + stamp
        if len(same) and not stamp:
            version = same['version'].iloc[0]     # a re-run on the kept download keeps the row as it was
        raw_rows.append({'file': path.name, 'url': BASE + name, 'bytes': path.stat().st_size, 'md5': md5,
                         'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'downloaded': when,
                         'version': version,
                         'kept': 'release_asset (Zenodo archive of the release; 44 MB, not in Git)'})
    upsert_tsv_rows(raw_list, raw_rows)

    command = (f'python scripts/55_stage_iptmnet_ptm_sites.py --fasta {args.fasta} --raw-dir PTM_data/marks --out-dir inputs/marks')
    rows = []
    for mark, name in files.items():
        path = out_dir / name
        rows.append({'file': f'marks/{name}', 'bytes': path.stat().st_size, 'md5': md5_file(path), 'rows': len(tables[mark]),
                     'redistribution': 'git', 'content': CONTENT[mark] + f': {len(tables[mark]):,} rows',
                     'source': BASE + 'ptm.txt', 'version_or_date': f'iPTMnet release {release}; downloaded {raw_rows[0]["downloaded"]}',
                     'citation': CITATION, 'licence': LICENSE, 'derived_from': f'{ptm_path.name} (MD5 {ptm_md5})',
                     'how_made_or_checked': command + '; residue checked on the reference sequence; counts in marks/' + summary_name})
    path = out_dir / summary_name
    rows.append({'file': f'marks/{summary_name}', 'bytes': path.stat().st_size, 'md5': md5_file(path), 'redistribution': 'git',
                 'content': 'Counts, release, drops and the MD5 of the input and of every table of the iPTMnet staging',
                 'source': 'the tables above', 'version_or_date': f'iPTMnet release {release}', 'derived_from': 'marks/raw_files.tsv',
                 'how_made_or_checked': command})
    update_manifest(args.manifest, rows)
    refresh_manifest_row(args.manifest, raw_list, 'marks/raw_files.tsv')
    print(json.dumps({k: summary[k] for k in ('iptmnet_release', 'human_records', 'dropped', 'rows', 'rows_by_state',
                                              'rows_by_evidence_class')}, indent=1))


if __name__ == '__main__':
    main()

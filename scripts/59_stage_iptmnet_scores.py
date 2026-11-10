#!/usr/bin/env python3
"""Stage iPTMnet's site scores (score.txt) for the sites of the context-mark tables (inputs/marks/).

    python scripts/59_stage_iptmnet_scores.py --fasta inputs/reference/uniprot_human_reviewed_2026_03.fasta \
        --raw-dir PTM_data/marks --out-dir inputs/marks

Reads score.txt of https://research.bioinformatics.udel.edu/iptmnet_data/files/current/ (iPTMnet, Protein Information
Resource, CC BY-NC-SA 4.0 per its license page and its download readme; this repository follows the stricter terms) and
the readme.txt beside it. With --raw it reads a kept download (readme: iptmnet_readme.txt beside it).

What the score is. The only description iPTMnet gives in the download is the readme's: "score.txt -- score information:
substrate_UniProtAC, site, enzyme_UniProtAC, ptm_type, score. Note: The score is calculated based on all the evidence
sources as shown on the website." The readme does not say how. In the file the score is an integer from 0 to 4. It is
iPTMnet's own score and is carried as given: it is not a probability, and it is not comparable between PTM types.

score.txt has no header and five tab-separated fields; the enzyme field is empty on the rows that score the site itself
(827,718 of 851,361 rows) and holds the kinase or enzyme on the rows that score an enzyme-site relation (those are not
used). Site rows are kept for the marks of the release: Phosphorylation at S, T or Y (phosphorylation), Methylation at K
(methyllysine) and Acetylation at K (lysine_acylation, the acetyl rows); other PTM types and residues are dropped, and
counted in the summary JSON with: not_in_reference (the accession is not a canonical entry of --fasta), on_an_isoform
(suffix -2 and up), residue_mismatch (the reference has another residue at the position: `common.residue_is`),
unparsed_site. Accessions with the suffix -1 are the canonical entry, as in script 55. score.txt scores a site under
both the plain accession (794,961 site rows) and the -1 accession (32,757): the plain accession's scores are almost always the
higher (the -1 rows mostly carry 0 or 1), so the plain accession's score is kept and the -1 one is used only where there is no
plain row (`plain_and_suffix_1_accession_conflicts` in the summary counts the sites where the two differ).

Columns of the staged table: canonical_UniProtAC, position, residue, mark, iptmnet_score, source (iPTMnet),
source_release. Script 56 adds iptmnet_score to the rows of the release tables whose sources include iPTMnet; a row of
another source only has it empty. The name of the staged table does not end in `_sites.tsv`, so script 56 does not read it as
a source of sites. Outputs are sorted, so a second run on the same download gives the same bytes.
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
FIELDS = ['acc', 'site', 'enzyme_acc', 'ptm_type', 'score']
COLUMNS = ['canonical_UniProtAC', 'position', 'residue', 'mark', 'iptmnet_score', 'source', 'source_release']
OUT_NAME = 'iptmnet_{rel}_score.tsv'
OUT_SUMMARY = 'iptmnet_score_summary.json'
README_SAYS = ('score.txt -- score information (substrate_UniProtAC, site, enzyme_UniProtAC, ptm_type, score). '
               'Note: The score is calculated based on all the evidence sources as shown on the website.')
CONTENT = ("iPTMnet's own score of each human site of the phosphorylation, lysine methylation and lysine acetylation tables "
           "(score.txt rows without an enzyme; integer 0-4; the readme says only that it is calculated from all the evidence sources shown on the website)")


def mark_of(ptm_type: str, letter: str):
    if ptm_type == 'Phosphorylation' and letter in 'STY':
        return 'phosphorylation'
    if ptm_type == 'Methylation' and letter == 'K':
        return 'methyllysine'
    if ptm_type == 'Acetylation' and letter == 'K':
        return 'lysine_acylation'
    return None


def stage(records: pd.DataFrame, reference: dict[str, str], release: str) -> tuple[pd.DataFrame, dict]:
    """The staged scores and the counts, from the rows of score.txt."""
    dropped, conflicts = Counter(), 0
    best: dict[tuple, tuple[bool, str]] = {}
    for r in records.itertuples(index=False):
        if r.enzyme_acc != '':
            dropped['enzyme_row'] += 1
            continue
        m = re.fullmatch(r'([A-Z])(\d+)', r.site)
        if not m:
            dropped['unparsed_site'] += 1
            continue
        letter, pos = m.group(1), int(m.group(2))
        mark = mark_of(r.ptm_type, letter)
        if mark is None:
            dropped['other_ptm_type_or_residue'] += 1
            continue
        plain = not r.acc.endswith('-1')
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
        key = (acc, pos, letter, mark)
        if key in best:
            was_plain, was = best[key]
            if was != r.score:
                conflicts += 1
            if was_plain or not plain:
                continue
        best[key] = (plain, r.score)
    rows = [{'canonical_UniProtAC': k[0], 'position': k[1], 'residue': k[2], 'mark': k[3], 'iptmnet_score': v[1],
             'source': 'iPTMnet', 'source_release': release} for k, v in best.items()]
    table = pd.DataFrame(rows, columns=COLUMNS)
    table['_p'] = table['position'].astype(int)
    table = table.sort_values(['mark', 'canonical_UniProtAC', '_p']).drop(columns='_p').reset_index(drop=True)
    return table, {**dict(dropped), '_conflicts': conflicts}


def fetch(url: str) -> tuple[bytes, str]:
    with urlopen(Request(url, headers={'User-Agent': 'PTM-pipeline-stager/1.0'}), timeout=1800) as response:
        return response.read(), response.headers.get('Last-Modified', '')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--raw-dir', default='PTM_data/marks', help='where the downloads are kept (not in Git)')
    ap.add_argument('--raw', help='read this score.txt (with iptmnet_readme.txt beside it) instead of fetching')
    ap.add_argument('--out-dir', default='inputs/marks')
    ap.add_argument('--manifest', default='inputs/MANIFEST.tsv')
    ap.add_argument('--downloaded', default='', help='date of the download for raw_files.tsv (default: today when fetching)')
    args = ap.parse_args()

    out_dir, raw_dir = Path(args.out_dir), Path(args.raw_dir)
    modified = {}
    if args.raw:
        score_path = Path(args.raw)
        readme_path = score_path.with_name('iptmnet_readme.txt')
    else:
        raw_dir.mkdir(parents=True, exist_ok=True)
        score_path, readme_path = raw_dir / 'iptmnet_score.txt', raw_dir / 'iptmnet_readme.txt'
        data, modified['score.txt'] = fetch(BASE + 'score.txt')
        score_path.write_bytes(data)
        data, modified['readme.txt'] = fetch(BASE + 'readme.txt')
        readme_path.write_bytes(data)
    readme = readme_path.read_text(encoding='utf-8')
    found = re.search(r'for release ([0-9][0-9.]*)', readme)
    if not found:
        raise SystemExit('the iPTMnet readme names no release')
    release = found.group(1)
    if 'score is calculated based on all the evidence sources' not in readme:
        raise SystemExit('the iPTMnet readme no longer describes the score as before: read it and update README_SAYS')

    records = pd.read_csv(score_path, sep='\t', header=None, names=FIELDS, dtype=str, keep_default_na=False)
    if not records['score'].str.fullmatch(r'\d+').all():
        raise SystemExit('score.txt has a score that is not an integer')
    reference = parse_fasta(args.fasta)
    table, dropped = stage(records, reference, release)
    conflicts = dropped.pop('_conflicts')
    out_dir.mkdir(parents=True, exist_ok=True)
    name = OUT_NAME.format(rel=release.replace('.', '_'))
    table.to_csv(out_dir / name, sep='\t', index=False, lineterminator='\n')
    score_md5, readme_md5 = md5_file(score_path), md5_file(readme_path)
    scores = records['score'].astype(int)
    summary = {
        'source': 'iPTMnet score.txt (Protein Information Resource)', 'url': BASE + 'score.txt', 'iptmnet_release': release,
        'license': LICENSE, 'readme_says': README_SAYS, 'reference_fasta': args.fasta, 'reference_fasta_md5': md5_file(args.fasta),
        'raw_files': {'score.txt': score_path.name, 'readme.txt': readme_path.name},
        'raw_md5': {'score.txt': score_md5, 'readme.txt': readme_md5},
        'records_in_download': len(records), 'site_records_without_enzyme': int((records['enzyme_acc'] == '').sum()),
        'score_range_in_download': [int(scores.min()), int(scores.max())], 'dropped': dropped,
        'plain_and_suffix_1_accession_conflicts': conflicts,
        'rows': len(table), 'rows_by_mark': {k: int(v) for k, v in table['mark'].value_counts().sort_index().items()},
        'rows_by_score': {m: {k: int(v) for k, v in g['iptmnet_score'].value_counts().sort_index().items()}
                          for m, g in table.groupby('mark')},
        'file': name, 'file_md5': md5_file(out_dir / name),
    }
    (out_dir / OUT_SUMMARY).write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n', encoding='utf-8')

    raw_list = out_dir / 'raw_files.tsv'
    old = pd.read_csv(raw_list, sep='\t', dtype=str, keep_default_na=False)
    same = old[(old['file'] == score_path.name) & (old['md5'] == score_md5)]
    when = same['downloaded'].iloc[0] if len(same) else (args.downloaded or date.today().isoformat())
    stamp = f'; Last-Modified {modified["score.txt"]}' if modified.get('score.txt') else ''
    version = f'iPTMnet release {release}' + stamp
    if len(same) and not stamp:
        version = same['version'].iloc[0]
    upsert_tsv_rows(raw_list, [{
        'file': score_path.name, 'url': BASE + 'score.txt', 'bytes': score_path.stat().st_size, 'md5': score_md5,
        'sha256': hashlib.sha256(score_path.read_bytes()).hexdigest(), 'downloaded': when, 'version': version,
        'kept': 'release_asset (Zenodo archive of the release; 26 MB, not in Git)'}])

    command = f'python scripts/59_stage_iptmnet_scores.py --fasta {args.fasta} --raw-dir PTM_data/marks --out-dir inputs/marks'
    path, spath = out_dir / name, out_dir / OUT_SUMMARY
    update_manifest(args.manifest, [
        {'file': f'marks/{name}', 'bytes': path.stat().st_size, 'md5': md5_file(path), 'rows': len(table), 'redistribution': 'git',
         'content': CONTENT + f': {len(table):,} rows', 'source': BASE + 'score.txt',
         'version_or_date': f'iPTMnet release {release}; downloaded {when}', 'citation': CITATION, 'licence': LICENSE,
         'derived_from': f'{score_path.name} (MD5 {score_md5})',
         'how_made_or_checked': command + f'; residue checked on the reference sequence; counts in marks/{OUT_SUMMARY}'},
        {'file': f'marks/{OUT_SUMMARY}', 'bytes': spath.stat().st_size, 'md5': md5_file(spath), 'redistribution': 'git',
         'content': 'Counts, release, drops, what the readme says the score is, and the MD5 of the input and of the table of the iPTMnet score staging',
         'source': 'the table above', 'version_or_date': f'iPTMnet release {release}', 'derived_from': 'marks/raw_files.tsv',
         'how_made_or_checked': command}])
    refresh_manifest_row(args.manifest, raw_list, 'marks/raw_files.tsv')
    print(json.dumps({k: summary[k] for k in ('iptmnet_release', 'records_in_download', 'dropped', 'rows', 'rows_by_mark',
                                              'plain_and_suffix_1_accession_conflicts')}, indent=1))


if __name__ == '__main__':
    main()

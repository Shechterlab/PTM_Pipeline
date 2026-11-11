#!/usr/bin/env python3
"""Stage and validate CPLM 4.0 lysine methylation (or ubiquitination, SUMOylation) sites against the canonical human proteome.

CPLM 4.0 (Zhang et al., Nucleic Acids Res 50:D451-D459, 2022, doi:10.1093/nar/gkab849; cplm.biocuckoo.cn)
compiles lysine modification sites. Each modification has its own download (--modification Methylation, the default,
Ubiquitination or Sumoylation: Download/<Modification>.zip, one headerless tab-separated <Modification>.txt of the same
nine columns); what follows says "Methylation" for all three. Its Methylation download (Download/Methylation.zip, one headerless
tab-separated Methylation.txt) has nine columns:

  CPLM id, UniProt accession, position, modification, gene, organism, protein sequence, evidence, references

Rows are kept when they are for the chosen organism (Homo sapiens), describe methylation, sit on a lysine of
CPLM's own sequence and on a lysine of the same accession and position in the canonical FASTA, and either the
two sequences are identical or the +/-7 residue windows around the lysine agree. Nothing is remapped: a site
whose sequence moved is reported in the rejected table with its reason. The evidence label ("Exp.", "Dat.",
"Exp./Dat.") and the references field are carried through verbatim; some "Dat." rows hold identifiers that
are not PubMed IDs.

Outputs (in --outdir, default $PTM_RESULTS_ROOT/cplm4)
  cplm4_human_lys_methyl_validated.tsv   one row per canonical accession and position
  cplm4_human_lys_methyl_rejected.tsv    input rows that failed a check, with the reason
  cplm4_summary.json                     input checksum, counts per check, sites and proteins kept

With --modification Ubiquitination or Sumoylation the same checks and outputs apply under other names: cplm4_human_lys_ubiquitination_validated.tsv
(and _rejected.tsv, cplm4_summary.json), and with --marks-dir cplm_4_0_ubiquitination_sites.tsv (mark ubiquitination, state ub) and
cplm_ubiquitination_summary.json, or cplm_4_0_sumoylation_sites.tsv (mark sumoylation, state sumo) and cplm_sumoylation_summary.json.
The Methylation outputs are unchanged by the option. Sumoylation is included in release v1.5.0; script 56 reads every
inputs/marks/*_sumoylation_sites.tsv.

With --marks-dir the validated sites are also staged for the context-mark release (inputs/marks/, read by
56_build_context_marks.py) in the format of the other marks:

    python scripts/39_stage_cplm4_lys_methylation.py --cplm PTM_data/marks/Methylation.zip \\
        --fasta inputs/reference/uniprot_human_reviewed_2026_03.fasta --marks-dir inputs/marks --downloaded 2026-10-18

  cplm_4_0_methyllysine_sites.tsv   the columns of the other staged tables, state `me` (CPLM gives no degree), source=CPLM,
                                    source_release=CPLM 4.0 (<date of Methylation.txt>), and two more columns:
                                    cplm_evidence (the label as CPLM gives it: Exp., Exp./Dat., Dat.) and cplm_references
                                    (CPLM's reference field, verbatim)
  cplm_ptm_sites_summary.json       counts, input checksums, the MD5 of the table
  raw_files.tsv, ../MANIFEST.tsv    rows for the download and the table

  evidence_class  Exp. and Exp./Dat.: experimental (CPLM curates these from experimental papers); Dat.: unspecified (the
                  record comes from another database and CPLM gives no experiment)
  pmids           the references that are PubMed IDs. CPLM's field also holds identifiers of other databases and its own
                  (Dat. rows); a number is taken as a PubMed ID when it has seven or eight digits and is at most
                  34,300,000, the PubMed IDs assigned when Methylation.txt was made (2021-06-13). Shorter and larger numbers
                  stay in cplm_references only.

CPLM 4.0 states no license (its page has a copyright notice, "All Rights Reserved"). Its human lysine methylation, ubiquitination and SUMOylation sites are
included in the release, cited (Zhang et al., Nucleic Acids Res 50:D451-D459, 2022, doi:10.1093/nar/gkab849).
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import io
import json
import os
import re
import zipfile
from collections import Counter
from pathlib import Path

import pandas as pd

from common import md5_file, parse_fasta, refresh_manifest_row, save_table, update_manifest, upsert_tsv_rows

CPLM4_COLUMNS = ['cplm_id', 'source_UniProtAC', 'position', 'modification', 'gene', 'organism', 'cplm_sequence',
                 'evidence', 'references']
WINDOW_FLANK = 7
MARKS_COLUMNS = ['canonical_UniProtAC', 'position', 'residue', 'mark', 'state', 'evidence', 'evidence_class', 'pmids', 'source',
                 'source_release', 'cplm_evidence', 'cplm_references']
URL_BASE = 'https://cplm.biocuckoo.cn/Download/'
DEFAULT_MODIFICATION = 'Methylation'
# per CPLM download: the tag of the validated/rejected tables, the staged mark and state, the staged table and summary names,
# the environment variable of --cplm, and what the manifest calls the sites
MODIFICATIONS = {
    'Methylation': {'tag': 'methyl', 'mark': 'methyllysine', 'state': 'me', 'table': 'cplm_4_0_methyllysine_sites.tsv',
                    'summary': 'cplm_ptm_sites_summary.json', 'env': 'PTM_CPLM4_METHYL',
                    'what': 'lysine methylation sites (degree not stated)'},
    'Ubiquitination': {'tag': 'ubiquitination', 'mark': 'ubiquitination', 'state': 'ub', 'table': 'cplm_4_0_ubiquitination_sites.tsv',
                       'summary': 'cplm_ubiquitination_summary.json', 'env': 'PTM_CPLM4_UBIQUITINATION', 'what': 'lysine ubiquitination sites'},
    'Sumoylation': {'tag': 'sumoylation', 'mark': 'sumoylation', 'state': 'sumo', 'table': 'cplm_4_0_sumoylation_sites.tsv',
                    'summary': 'cplm_sumoylation_summary.json', 'env': 'PTM_CPLM4_SUMOYLATION', 'what': 'lysine SUMOylation sites'},
}
CITATION = ('Zhang W, Tan X, Lin S, et al. CPLM 4.0: an updated database with rich annotations for protein lysine modifications. '
            'Nucleic Acids Res 50:D451-D459 (2022), doi:10.1093/nar/gkab849')
LICENSE = 'none stated (cplm.biocuckoo.cn carries a copyright notice, "All Rights Reserved")'
LICENSE_UBIQUITINATION = LICENSE
LICENSE_SUMOYLATION = LICENSE
MAX_PMID = 34_300_000    # the highest PubMed ID assigned by the date of Methylation.txt (2021-06-13)
MIN_PMID = 1_000_000
EVIDENCE_CLASS = {'Exp.': 'experimental', 'Exp./Dat.': 'experimental', 'Dat.': 'unspecified'}


def read_cplm(path: Path, modification: str = DEFAULT_MODIFICATION) -> tuple[pd.DataFrame, str, str]:
    """Read a CPLM download given as the zip, the extracted text file, or a directory holding either."""
    if path.is_dir():
        for name in (f'{modification}.txt', f'{modification}.zip'):
            if (path / name).exists():
                path = path / name
                break
        else:
            raise SystemExit(f'no {modification}.txt or {modification}.zip in {path}')
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    member = ''
    if zipfile.is_zipfile(io.BytesIO(raw)):
        with zipfile.ZipFile(io.BytesIO(raw)) as zf:
            names = [n for n in zf.namelist() if n.lower().endswith('.txt')]
            if len(names) != 1:
                raise SystemExit(f'expected one .txt file in {path}, found {names}')
            info = zf.getinfo(names[0])
            member = f'{info.filename} ({dt.datetime(*info.date_time).isoformat()})'
            raw = zf.read(names[0])
    df = pd.read_csv(io.BytesIO(raw), sep='\t', header=None, dtype=str, keep_default_na=False, quoting=3)
    if df.shape[1] != len(CPLM4_COLUMNS):
        raise SystemExit(f'{path}: expected {len(CPLM4_COLUMNS)} columns, found {df.shape[1]}')
    df.columns = CPLM4_COLUMNS
    return df, digest, member


def pubmed_ids(references: str) -> list[str]:
    """The references of a CPLM row that are PubMed IDs (see the module docstring), sorted, without repeats."""
    ids = {int(t) for t in references.split(';') if re.fullmatch(r'\d{7,8}', t) and MIN_PMID <= int(t) <= MAX_PMID}
    return [str(i) for i in sorted(ids)]


def stage_marks(good: pd.DataFrame, release: str, cfg: dict = MODIFICATIONS[DEFAULT_MODIFICATION]) -> pd.DataFrame:
    """The validated sites as a staged table of inputs/marks/ (one row per residue, one state)."""
    out = pd.DataFrame({
        'canonical_UniProtAC': good['canonical_UniProtAC'], 'position': good['position'].astype(int), 'residue': 'K',
        'mark': cfg['mark'], 'state': cfg['state'], 'evidence': '',
        'evidence_class': good['evidence'].map(EVIDENCE_CLASS),
        'pmids': [';'.join(pubmed_ids(r)) for r in good['references']], 'source': 'CPLM', 'source_release': release,
        'cplm_evidence': good['evidence'], 'cplm_references': good['references']})
    if out['evidence_class'].isna().any():
        raise SystemExit(f'unknown CPLM evidence label: {sorted(set(good["evidence"][out["evidence_class"].isna()]))}')
    return out.sort_values(['canonical_UniProtAC', 'position'], kind='stable')[MARKS_COLUMNS].reset_index(drop=True)


def window(seq: str, pos: int) -> str:
    return seq[max(0, pos - 1 - WINDOW_FLANK):pos + WINDOW_FLANK]


def check_row(row, canonical: dict[str, str]) -> tuple[str, str]:
    """Return (status, detail): status is 'identical_sequence' or 'window_match' when kept, else a reason."""
    pos = row.position_int
    own = row.cplm_sequence.strip().upper()
    if pos is None or pos < 1:
        return 'bad_position', ''
    if not (pos <= len(own) and own[pos - 1] == 'K'):
        return 'cplm_sequence_not_K', own[pos - 1] if pos <= len(own) else 'past end'
    seq = canonical.get(row.source_UniProtAC)
    if seq is None:
        return 'accession_not_in_canonical_fasta', ''
    if not (pos <= len(seq) and seq[pos - 1] == 'K'):
        return 'canonical_not_K', seq[pos - 1] if pos <= len(seq) else 'past end'
    if own == seq:
        return 'identical_sequence', ''
    if window(own, pos) == window(seq, pos):
        return 'window_match', ''
    return 'window_mismatch', f'{window(own, pos)} vs {window(seq, pos)}'


def write_marks(args, good: pd.DataFrame, summary: dict, rejected: pd.DataFrame, raw: Path, member: str) -> None:
    """The staged table, its summary, the download's row in raw_files.tsv and the MANIFEST rows."""
    mod = args.modification
    cfg = MODIFICATIONS[mod]
    marks_name, marks_summary, mark, state = cfg['table'], cfg['summary'], cfg['mark'], cfg['state']
    url = URL_BASE + f'{mod}.zip'
    terms = {'Ubiquitination': LICENSE_UBIQUITINATION, 'Sumoylation': LICENSE_SUMOYLATION}.get(mod, LICENSE)
    marks = Path(args.marks_dir)
    stamp = re.search(r'\((\d{4}-\d{2}-\d{2})T', member)
    if not stamp:
        raise SystemExit(f'no date in the member of the download: {member!r}')
    release = f'CPLM 4.0 ({stamp.group(1)})'
    table = stage_marks(good, release, cfg)
    marks.mkdir(parents=True, exist_ok=True)
    table.to_csv(marks / marks_name, sep='\t', index=False, lineterminator='\n')
    if raw.is_dir():
        raw = raw / (f'{mod}.zip' if (raw / f'{mod}.zip').exists() else f'{mod}.txt')
    raw_md5 = md5_file(raw)
    pm = table['pmids'].map(bool)
    out = {
        'source': f'CPLM 4.0 {mod} download (cplm.biocuckoo.cn)', 'url': url, 'cplm_release': release, 'license': terms,
        'reference_fasta': args.fasta, 'reference_fasta_md5': md5_file(args.fasta),
        'raw_files': {f'{mod}.zip': raw.name}, 'raw_md5': {f'{mod}.zip': raw_md5},
        'raw_sha256': {f'{mod}.zip': summary['input_sha256']}, 'raw_member': member,
        'records_in_download': summary['input_rows'], 'records_in_scope': summary['rows_in_scope'],
        'dropped': {k: v for k, v in sorted(summary['checks'].items()) if k not in ('identical_sequence', 'window_match')},
        'checks': dict(sorted(summary['checks'].items())),
        'rows': {mark: len(table)}, 'proteins': int(table['canonical_UniProtAC'].nunique()),
        'rows_by_state': {mark: {state: len(table)}},
        'rows_by_evidence_class': {mark: {k: int(v) for k, v in table['evidence_class'].value_counts().sort_index().items()}},
        'rows_by_cplm_evidence': {k: int(v) for k, v in table['cplm_evidence'].value_counts().sort_index().items()},
        'rows_with_pubmed_ids': int(pm.sum()), 'pubmed_id_rule': f'{len(str(MIN_PMID))}-8 digits, {MIN_PMID:,} to {MAX_PMID:,}',
        'files': {mark: marks_name}, 'files_md5': {mark: md5_file(marks / marks_name)},
    }
    (marks / marks_summary).write_text(json.dumps(out, indent=2, sort_keys=True) + '\n', encoding='utf-8')

    raw_list = marks / 'raw_files.tsv'
    old = pd.read_csv(raw_list, sep='\t', dtype=str, keep_default_na=False)
    same = old[(old['file'] == raw.name) & (old['md5'] == raw_md5)]
    when = same['downloaded'].iloc[0] if len(same) else (args.downloaded or dt.date.today().isoformat())
    upsert_tsv_rows(raw_list, [{
        'file': raw.name, 'url': url, 'bytes': raw.stat().st_size, 'md5': raw_md5, 'sha256': summary['input_sha256'],
        'downloaded': when, 'version': f'CPLM 4.0; member {member}',
        'kept': f'not in Git, not redistributed ({raw.stat().st_size / 1e6:.1f} MB; download it again from the URL and check the SHA-256); the staged table is in Git'}])

    command = (f'python scripts/39_stage_cplm4_lys_methylation.py --cplm PTM_data/marks/{raw.name} --fasta {args.fasta} '
               f'--marks-dir inputs/marks' + ('' if mod == DEFAULT_MODIFICATION else f' --modification {mod}'))
    manifest_rows = [
        {'file': f'marks/{marks_name}', 'bytes': (marks / marks_name).stat().st_size, 'md5': md5_file(marks / marks_name), 'rows': len(table),
         'redistribution': 'git',
         'content': f'CPLM 4.0 human {cfg["what"]}, one row per residue, with CPLM\'s evidence label (Exp., Exp./Dat., Dat.) and references: {len(table):,} rows',
         'source': url, 'version_or_date': f'{release}; downloaded {when}', 'citation': CITATION, 'licence': terms,
         'derived_from': f'{raw.name} (MD5 {raw_md5})',
         'how_made_or_checked': command + '; sequence and residue checked on the reference sequence; counts in marks/' + marks_summary},
        {'file': f'marks/{marks_summary}', 'bytes': (marks / marks_summary).stat().st_size, 'md5': md5_file(marks / marks_summary),
         'redistribution': 'git', 'content': 'Counts, release, drops and the MD5 of the input and of the table of the CPLM staging',
         'source': 'the table above', 'version_or_date': release, 'derived_from': 'marks/raw_files.tsv', 'how_made_or_checked': command}]
    update_manifest(args.manifest, manifest_rows)
    refresh_manifest_row(args.manifest, raw_list, 'marks/raw_files.tsv')


def main() -> None:
    results_root = Path(os.environ.get('PTM_RESULTS_ROOT', 'PTM_results'))
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--modification', default=DEFAULT_MODIFICATION, choices=sorted(MODIFICATIONS),
                    help='which CPLM download: Methylation (default), Ubiquitination or Sumoylation')
    ap.add_argument('--cplm', default='',
                    help='CPLM 4.0 <Modification>.zip, <Modification>.txt, or a directory holding one '
                         '(default: the environment variable PTM_CPLM4_METHYL, PTM_CPLM4_UBIQUITINATION or PTM_CPLM4_SUMOYLATION).')
    ap.add_argument('--fasta', default=os.environ.get('PTM_CANONICAL_FASTA', ''),
                    help='UniProt reviewed canonical FASTA (PTM_CANONICAL_FASTA).')
    ap.add_argument('--organism', default='Homo sapiens')
    ap.add_argument('--outdir', default=str(results_root / 'cplm4'))
    ap.add_argument('--marks-dir', default='', help='also stage the validated sites here (inputs/marks) for the context-mark release')
    ap.add_argument('--manifest', default='inputs/MANIFEST.tsv')
    ap.add_argument('--downloaded', default='', help='date of the download for raw_files.tsv (default: today)')
    args = ap.parse_args()
    cfg = MODIFICATIONS[args.modification]
    args.cplm = args.cplm or os.environ.get(cfg['env'], '')
    if not args.cplm or not args.fasta:
        raise SystemExit(f'--cplm and --fasta are required (or set {cfg["env"]} and PTM_CANONICAL_FASTA)')

    df, digest, member = read_cplm(Path(args.cplm), args.modification)
    canonical = parse_fasta(args.fasta)
    n_input = len(df)
    df = df[df['organism'].eq(args.organism) & df['modification'].str.contains(args.modification, case=False)].copy()
    n_scope = len(df)
    df['position_int'] = [int(p) if p.strip().isdigit() else None for p in df['position']]
    checks = [check_row(row, canonical) for row in df.itertuples(index=False)]
    df['validation'] = [c[0] for c in checks]
    df['validation_detail'] = [c[1] for c in checks]
    kept = df['validation'].isin(['identical_sequence', 'window_match'])

    rejected = df[~kept].drop(columns=['cplm_sequence', 'position_int'])
    good = df[kept].copy()
    good['canonical_UniProtAC'] = good['source_UniProtAC']
    good['position'] = good['position_int'].astype(int)
    good['residue'] = 'K'
    good['site'] = 'K' + good['position'].astype(str)
    good = (
        good.groupby(['canonical_UniProtAC', 'position'], as_index=False, sort=True)
        .agg(residue=('residue', 'first'), site=('site', 'first'), gene=('gene', 'first'),
             cplm_id=('cplm_id', lambda v: ';'.join(sorted(set(v)))),
             modification=('modification', 'first'),
             evidence=('evidence', lambda v: '/'.join(e for e in ('Exp.', 'Dat.') if any(e in x.split('/') for x in v))),
             references=('references', lambda v: ';'.join(dict.fromkeys(t for x in v for t in x.split(';') if t))),
             validation=('validation', 'first'))
    )

    outdir = Path(args.outdir)
    save_table(good, outdir / f'cplm4_human_lys_{cfg["tag"]}_validated.tsv')
    save_table(rejected, outdir / f'cplm4_human_lys_{cfg["tag"]}_rejected.tsv')
    summary = {
        'generated': dt.date.today().isoformat(),
        'input': str(Path(args.cplm)),
        'input_sha256': digest,
        'input_member': member,
        'canonical_fasta': str(Path(args.fasta)),
        'canonical_sequences': len(canonical),
        'organism': args.organism,
        'input_rows': n_input,
        'rows_in_scope': n_scope,
        'checks': dict(Counter(df['validation'])),
        'validated_sites': len(good),
        'validated_proteins': int(good['canonical_UniProtAC'].nunique()),
        'validated_evidence': dict(Counter(good['evidence'])),
    }
    (outdir / 'cplm4_summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))
    if args.marks_dir:
        write_marks(args, good, summary, rejected, Path(args.cplm), member)


if __name__ == '__main__':
    main()

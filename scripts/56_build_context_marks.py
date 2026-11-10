#!/usr/bin/env python3
"""Write releases/v1.5.0: the methylarginine data of v1.4.3, carried over unchanged, and the context-mark tables.

    python scripts/56_build_context_marks.py --version v1.5.0 --release-date 2026-11-09

The methylarginine tables of v1.5.0 are the v1.4.3 files, copied under this release's own names (_v1.4.3 becomes
_v1.5.0, as every release names all its files) with the data unchanged (the build, hpc/build_release.sh, needs R and the
network; this release changes no methylarginine data). Only the strings that name the version inside a file change, as when
v1.4.3 followed v1.4.2: the title and the statistics file name in the two analyses READMEs, the file names in their MD5SUMS,
the title of STUDIES.md and the title and sheet name of the workbook. supporting/site_changes_v1.4.3_to_v1.5.0.tsv is written
in the format of build_release.build_changes: every site kept, none added or removed. Beside them, for each
context mark (methyllysine, phosphorylation), `human_<mark>_sites_<version>.tsv` merges the staged tables
`inputs/marks/*_<mark>_sites.tsv` of every source (UniProt, iPTMnet, and CPLM 4.0 for methyllysine), one row per residue and state, and
`human_methylarginine_crosscheck_<version>.tsv` merges the sources' methylarginine records with a column that says whether
the arginine is in the methylarginine data. Also written: data_dictionary.tsv (v1.4.3's, with the new tables), README.md,
CHANGES.md, release.json and MD5SUMS. Nothing in the output depends on the date but --release-date, and the commit named
is the last one that changed scripts 39, 54, 55 or 56, so a second run gives the same bytes.

Merged columns:
  canonical_UniProtAC, position, residue, mark, state
  sources          the sources that report the residue and state, in order: UniProt; iPTMnet; CPLM
  evidence_class   the best class of any source: experimental, then text-mined, then by-similarity, then unspecified
                   (CPLM: its labels Exp. and Exp./Dat. are experimental, Dat. is unspecified)
  evidence         UniProt's ECO codes ('' where only iPTMnet reports it)
  iptmnet_sources  iPTMnet's upstream databases for the residue (hprd, unip, rlim, rlim+, sign, npro, iedb, pro, intact)
  pmids            the PubMed IDs of all sources
  source_releases  the release of each source in `sources`, same order
  cplm_evidence    (methyllysine only) CPLM's own evidence label for the residue: Exp., Exp./Dat. or Dat. ('' where CPLM does not report it)

CPLM 4.0 states no license. The owner decided on 2026-11-09 to include its human lysine methylation sites, cited, pending the
authors' permission; their terms, when given, are followed. CPLM gives no methylation degree: its rows are state `me`.
"""
from __future__ import annotations

import argparse
import io
import json
import re
import subprocess
import zipfile
from pathlib import Path

import pandas as pd

from common import md5_file

MARKS = ['methyllysine', 'phosphorylation']
SOURCE_ORDER = ['UniProt', 'iPTMnet', 'CPLM']
CARRIED_FROM = 'v1.4.3'
NOT_CARRIED = {'README.md', 'CHANGES.md', 'MD5SUMS', 'data_dictionary.tsv', 'release.json'}
# files that name their release inside, and whose naming strings change with the version
VERSIONED_TEXT = {'analyses/README.md', 'analyses/MD5SUMS', 'analyses_with_lower_confidence/README.md',
                  'analyses_with_lower_confidence/MD5SUMS', 'STUDIES.md'}
CHANGE_COLUMNS = ['site_id', 'previous_site_id', 'gene', 'uniprot_accession', 'position', 'in_previous_version',
                  'in_this_version', 'change', 'reason_codes', 'records_change', 'records_now_at_or_previously_at', 'record_ids']
CLASS_RANK = ['experimental', 'text-mined', 'by-similarity', 'unspecified']
OUT_COLUMNS = ['canonical_UniProtAC', 'position', 'residue', 'mark', 'state', 'sources', 'evidence_class', 'evidence',
               'iptmnet_sources', 'pmids', 'source_releases']
CPLM_COLUMN = 'cplm_evidence'     # after OUT_COLUMNS, in the tables that have CPLM rows
SOURCES = {
    'UniProt': {'license': 'CC BY 4.0', 'url': 'https://www.uniprot.org/help/license',
                'citation': 'UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. Nucleic Acids Res '
                            '53:D609-D617 (2025). doi:10.1093/nar/gkae1010'},
    'iPTMnet': {'license': 'CC BY-NC-SA 4.0', 'url': 'https://research.bioinformatics.udel.edu/iptmnet/license',
                'citation': 'Huang H, Arighi CN, Ross KE, et al. iPTMnet: an integrated resource for protein post-translational '
                            'modification network discovery. Nucleic Acids Res 46:D542-D550 (2018). doi:10.1093/nar/gkx1104'},
    'CPLM': {'license': 'none stated; included pending the authors\' permission', 'url': 'https://cplm.biocuckoo.cn/',
             'citation': 'Zhang W, Tan X, Lin S, et al. CPLM 4.0: an updated database with rich annotations for protein lysine '
                         'modifications. Nucleic Acids Res 50:D451-D459 (2022). doi:10.1093/nar/gkab849 (included pending the authors\' permission)'},
    # Berryhill2023 is added only once the license of its supplements is confirmed (inputs/marks/README.md).
}
UPSTREAM = {'hprd': 'HPRD', 'unip': 'UniProt', 'rlim': 'RLIMS-P text mining', 'rlim+': 'RLIMS-P text mining', 'sign': 'SIGNOR',
            'npro': 'neXtProt', 'iedb': 'IEDB', 'pro': 'Protein Ontology', 'intact': 'IntAct'}
CODES = {
    'state': 'methyllysine: me1, me2, me3, me (degree not stated); phosphorylation: phospho; methylarginine: me1, me2a, me2s, me2 (type not stated), me (degree not stated)',
    'evidence_class': 'experimental, text-mined, by-similarity, unspecified',
    'mark': 'methyllysine, phosphorylation, methylarginine',
    'sources': 'UniProt, iPTMnet, CPLM',
    'cplm_evidence': 'Exp., Exp./Dat., Dat.',
}
DICTIONARY = [
    ('canonical_UniProtAC', 'string', 'reviewed UniProtKB accession of the canonical sequence in the reference proteome (UniProt 2026_03, UP000005640)'),
    ('position', 'integer', 'residue number on that sequence, 1-based'),
    ('residue', 'code', 'one-letter code of the residue, checked against the reference sequence (K for methyllysine, S, T or Y for phosphorylation, R for methylarginine)'),
    ('mark', 'code', 'the modification'),
    ('state', 'code', 'the state reported for the residue; a residue with several states has one row per state. iPTMnet and CPLM do not give a methylation degree (state me)'),
    ('sources', 'list<code>', 'the sources that report the residue and state, separated by "; "'),
    ('evidence_class', 'code', 'the best class of any source. experimental: UniProt ECO:0000269 or ECO:0007744, an iPTMnet record curated from an experiment (upstream hprd, sign, npro, pro, intact, iedb), or a CPLM record labelled Exp. or Exp./Dat. (CPLM curates these from experimental papers). text-mined: only RLIMS-P reports it. by-similarity: UniProt records with no experimental code (mainly ECO:0000250). unspecified: only iPTMnet\'s copy of UniProt (upstream unip), which has no ECO code, or only a CPLM record labelled Dat. (taken from another database; CPLM names no experiment). A site can hide by-similarity and unspecified rows by filtering on this column'),
    ('evidence', 'list<string>', 'UniProt ECO evidence codes, separated by ";" (empty where only iPTMnet reports the residue)'),
    ('iptmnet_sources', 'list<code>', 'the upstream databases behind the iPTMnet record, separated by ";": hprd (HPRD), unip (UniProt), rlim and rlim+ (RLIMS-P text mining), sign (SIGNOR), npro (neXtProt), iedb (IEDB), pro (Protein Ontology), intact (IntAct); iPTMnet asks that they are cited'),
    ('pmids', 'list<string>', 'PubMed IDs of the records, separated by ";", when given. CPLM\'s reference field also holds identifiers that are not PubMed IDs; only seven- and eight-digit numbers up to 34,300,000 (the PubMed IDs assigned by the date of CPLM\'s file) are taken as PubMed IDs'),
    ('source_releases', 'list<string>', 'the release of each source in sources, same order (UniProt release, iPTMnet release, CPLM release and the date of its file)'),
    ('cplm_evidence', 'code', 'methyllysine table only: CPLM\'s own evidence label for the residue, as CPLM gives it. Exp.: from an experimental paper; Exp./Dat.: both an experimental paper and another database; Dat.: another database only. Empty where CPLM does not report the residue. CPLM 4.0 states no license; its sites are included by the owner\'s decision pending the authors\' permission, cited to Zhang et al., Nucleic Acids Res 50:D451 (2022)'),
    ('in_methylarginine_release', 'code', f'crosscheck table only: whether the arginine is a site, or a lower-confidence site, in the methylarginine data of this release (unchanged from {CARRIED_FROM}), or not at all. Codes: site, lower_confidence, no'),
]


def read_stage(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    if 'upstream_sources' not in frame:
        frame['upstream_sources'] = ''
    if CPLM_COLUMN not in frame:
        frame[CPLM_COLUMN] = ''
    return frame


def best_class(classes) -> str:
    return min(set(classes), key=CLASS_RANK.index)


def union(grp: pd.DataFrame, col: str) -> list[str]:
    return sorted({v for x in grp[col] for v in x.split(';') if v})


def merge(frames: list[pd.DataFrame]) -> pd.DataFrame:
    """One row per residue and state from the staged rows of every source."""
    table = pd.concat(frames, ignore_index=True)
    for col in ('upstream_sources', CPLM_COLUMN):
        if col not in table:
            table[col] = ''
    table = table.fillna('')
    rows = []
    for key, grp in table.groupby(['canonical_UniProtAC', 'position', 'residue', 'mark', 'state'], sort=False):
        present = [s for s in SOURCE_ORDER if s in set(grp['source'])]
        release = {s: grp.loc[grp['source'] == s, 'source_release'].iloc[0] for s in present}
        rows.append({'canonical_UniProtAC': key[0], 'position': key[1], 'residue': key[2], 'mark': key[3], 'state': key[4],
                     'sources': '; '.join(present), 'evidence_class': best_class(grp['evidence_class']),
                     'evidence': ';'.join(union(grp, 'evidence')), 'iptmnet_sources': ';'.join(union(grp, 'upstream_sources')),
                     'pmids': ';'.join(sorted(union(grp, 'pmids'), key=int)), 'source_releases': '; '.join(release[s] for s in present),
                     CPLM_COLUMN: ';'.join(union(grp, CPLM_COLUMN))})
    out = pd.DataFrame(rows, columns=OUT_COLUMNS + [CPLM_COLUMN])
    if 'CPLM' not in set(table['source']):
        out = out.drop(columns=CPLM_COLUMN)
    out['_pos'] = out['position'].astype(int)
    return out.sort_values(['canonical_UniProtAC', '_pos', 'state']).drop(columns='_pos').reset_index(drop=True)


def crosscheck(merged: pd.DataFrame, release_dir: Path, version: str = CARRIED_FROM) -> pd.DataFrame:
    def keys(name: str) -> set[tuple[str, str]]:
        t = pd.read_csv(release_dir / name, sep='\t', dtype=str, keep_default_na=False, usecols=['uniprot_accession', 'position'])
        return set(zip(t['uniprot_accession'], t['position'], strict=True))
    sites = keys(f'human_methylarginine_sites_{version}.tsv')
    lower = keys(f'human_methylarginine_lower_confidence_sites_{version}.tsv')
    out = merged.copy()
    out['in_methylarginine_release'] = ['site' if k in sites else 'lower_confidence' if k in lower else 'no'
                                        for k in zip(out['canonical_UniProtAC'], out['position'], strict=True)]
    return out


def combo(sources: str) -> str:
    """'UniProt only', or the sources joined by ' + ' in the order of SOURCE_ORDER."""
    s = [x for x in SOURCE_ORDER if x in set(sources.split('; '))]
    return f'{s[0]} only' if len(s) == 1 else ' + '.join(s)


def combo_order(label: str) -> tuple:
    names = label.replace(' only', '').split(' + ')
    return (len(names), [SOURCE_ORDER.index(n) for n in names])


def overlap(table: pd.DataFrame) -> dict:
    """Rows and residues (any state) by which sources report them: 'UniProt only', 'UniProt + iPTMnet', ..."""
    rows = table['sources'].map(combo).value_counts()
    sites = table.groupby(['canonical_UniProtAC', 'position'])['sources'].agg(lambda v: '; '.join(sorted({x for s in v for x in s.split('; ')})))
    sites = sites.map(combo).value_counts()
    return {'rows': {k: int(v) for k, v in rows.items()}, 'sites': {k: int(v) for k, v in sites.items()}}


def counts_table(table: pd.DataFrame) -> list[str]:
    classes = [c for c in CLASS_RANK if c in set(table['evidence_class'])]
    combos = sorted(set(table['sources'].map(combo)), key=combo_order)
    lines = ['| State | ' + ' | '.join(classes) + ' | ' + ' | '.join(combos) + ' | Rows |', '|---|' + '---:|' * (len(classes) + len(combos) + 1)]
    def line(name, grp):
        by = grp['sources'].map(combo).value_counts()
        return (f'| {name} | ' + ' | '.join(f'{(grp["evidence_class"] == c).sum():,}' for c in classes)
                + ' | ' + ' | '.join(f'{by.get(c, 0):,}' for c in combos) + f' | {len(grp):,} |')
    for state, grp in table.groupby('state'):
        lines.append(line(state, grp))
    lines.append(line('all', table))
    return lines


def renamed(name: str, version: str) -> str:
    """A carried file's name with this release's version (the version in a name is the release's own)."""
    return name.replace(f'_{CARRIED_FROM}', f'_{version}')


def retitle_workbook(data: bytes, version: str) -> bytes:
    """The workbook with its title and sheet name ('Sites v1.4.3') at this version; every other part as it was."""
    out = io.BytesIO()
    with zipfile.ZipFile(io.BytesIO(data)) as zin, zipfile.ZipFile(out, 'w') as zout:
        for info in zin.infolist():
            part = zin.read(info.filename)
            if info.filename in ('docProps/core.xml', 'xl/workbook.xml'):
                part = part.replace(f'Sites {CARRIED_FROM}'.encode(), f'Sites {version}'.encode())
            zout.writestr(info, part, compress_type=info.compress_type)
    return out.getvalue()


def site_changes(sites: pd.DataFrame, previous_changes: pd.DataFrame) -> pd.DataFrame:
    """supporting/site_changes_<previous>_to_<version>.tsv for a release whose sites are the previous release's: every
    site kept, with its records the same (build_release.build_changes gives these rows when nothing changes). A gene the
    entry lacks is the one the previous table took from the records."""
    reported = dict(zip(previous_changes['site_id'], previous_changes['gene'], strict=True))
    out = sites[['site_id', 'gene', 'uniprot_accession', 'position']].copy()
    out['gene'] = [g or reported.get(sid, '') for sid, g in zip(out['site_id'], out['gene'], strict=True)]
    out['previous_site_id'] = out['site_id']
    out['in_previous_version'], out['in_this_version'], out['change'] = 'true', 'true', 'kept'
    out['reason_codes'], out['records_change'], out['records_now_at_or_previously_at'], out['record_ids'] = '', 'same', '', ''
    out['_p'] = out['position'].astype(int)
    out = out.sort_values(['uniprot_accession', '_p'], kind='stable')
    return out[CHANGE_COLUMNS].reset_index(drop=True)


def code_commit() -> str:
    paths = ['scripts/39_stage_cplm4_lys_methylation.py', 'scripts/54_stage_uniprot_ptm_sites.py', 'scripts/55_stage_iptmnet_ptm_sites.py', 'scripts/56_build_context_marks.py']
    try:
        dirty = subprocess.run(['git', 'status', '--porcelain', '--', *paths], capture_output=True, text=True, check=True).stdout
        sha = subprocess.run(['git', 'log', '-1', '--format=%H', '--', *paths], capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return 'unknown'
    return sha + (' (scripts modified)' if dirty.strip() else '')


def demote(text: str) -> str:
    out, fence = [], False
    for line in text.splitlines():
        if line.startswith('```'):
            fence = not fence
        out.append('#' + line if line.startswith('#') and not fence else line)
    return '\n'.join(out) + '\n'


def write_docs(out: Path, version, release_date, tables, cross, summaries, commit) -> dict:
    files = {m: f'human_{m}_sites_{version}.tsv' for m in MARKS}
    files['methylarginine'] = f'human_methylarginine_crosscheck_{version}.tsv'
    allt = {**tables, 'methylarginine': cross}
    releases = {s: set() for s in SOURCE_ORDER}
    for t in allt.values():
        for names, rels in zip(t['sources'], t['source_releases'], strict=True):
            for name, rel in zip(names.split('; '), rels.split('; '), strict=True):
                releases[name].add(rel)
    releases = {s: sorted(v) for s, v in releases.items()}
    old = (out / f'README_{CARRIED_FROM}.tmp').read_text(encoding='utf-8')
    L = [f'# Human methylarginine sites and context marks, {version}', '',
         f'{version} keeps the methylarginine data of {CARRIED_FROM} unchanged (`human_methylarginine_*`, '
         '`statistics_*`, `studies_*`, `supporting/` and `analyses*` are the v1.4.3 files, renamed `_v1.5.0` and otherwise unchanged apart from the '
         'strings that name the version inside a few of them; see CHANGES.md) and adds '
         'tables of context marks: the lysine methylation sites (UniProt, iPTMnet and CPLM 4.0) and phosphorylation sites (UniProt and iPTMnet) of the human reviewed proteome, and the '
         'arginine methylation records of UniProt and iPTMnet as a cross-check. Sites for each mark read their own table and show the '
         'others as context. Positions of the new tables are on UniProt release 2026_03, the reference of the methylarginine data.', '',
         '## Context marks', '', '| File | One row per |', '|---|---|']
    for m in MARKS:
        t = tables[m]
        L.append(f'| `{files[m]}` | {m} residue and state ({len(t):,} rows on {t[["canonical_UniProtAC", "position"]].drop_duplicates().shape[0]:,} '
                 f'residues of {t["canonical_UniProtAC"].nunique():,} proteins) |')
    L += [f'| `{files["methylarginine"]}` | arginine methylation record of UniProt or iPTMnet ({len(cross):,} rows), with whether the methylarginine data has the arginine. A cross-check, not a second methylarginine table |', '',
          'Columns are in `data_dictionary.tsv`. Each row says which sources report it (`sources`), the best `evidence_class`, UniProt\'s ECO codes, '
          'iPTMnet\'s upstream databases (`iptmnet_sources`), the PubMed IDs and, in the methyllysine table, CPLM\'s own evidence label (`cplm_evidence`: Exp., Exp./Dat., Dat.). A site can hide `by-similarity` (UniProt ECO:0000250 and other '
          'non-experimental records) and `unspecified` rows by filtering on `evidence_class`.', '',
          '### Sources and licenses', '', '| Source | Release | License | Terms |', '|---|---|---|---|']
    for s in SOURCE_ORDER:
        L.append(f'| {s} | {", ".join(releases[s])} | {SOURCES[s]["license"]} | {SOURCES[s]["url"]} |')
    L += ['', 'The release is CC BY-NC-SA 4.0 ([LICENSE.md](../../LICENSE.md)). UniProt (CC BY 4.0) is compatible with it. iPTMnet\'s license page says '
          'CC BY-NC-SA 4.0; its download readme also mentions CC BY 4.0 in its last line. The release follows the stricter terms: share-alike and '
          'non-commercial. iPTMnet asks that the upstream databases of its `source` column are cited when its data are redistributed: each row keeps them '
          '(`iptmnet_sources`), and they are HPRD, UniProt, RLIMS-P text mining, SIGNOR, neXtProt, IEDB, Protein Ontology and IntAct (terms are theirs). '
          'Not used: PhosphoSitePlus, dbPTM and qPTM. The supplements of Berryhill et al. 2023 are not used: they state no license of their own.', '',
          '**CPLM 4.0 states no license** (cplm.biocuckoo.cn carries a copyright notice, "All Rights Reserved"). Its human lysine methylation sites are included by the '
          'owner\'s decision of 2026-11-09, cited (Zhang et al., Nucleic Acids Res 50:D451-D459, 2022, doi:10.1093/nar/gkab849), **pending the authors\' permission**; '
          'the terms the authors give, when they answer, are followed. The CC BY-NC-SA 4.0 of this release does not grant rights in CPLM\'s records beyond what the authors allow. '
          'CPLM contributes only to the methyllysine table, as the rows whose `sources` name CPLM (state `me`, with `cplm_evidence`); a reader who needs the table without CPLM removes that source from `sources`, '
          'and the rows that only CPLM reports.', '',
          '### Counts', '',
          f'Rows by state, evidence class and the sources that report them (iPTMnet {summaries["iptmnet"]["iptmnet_release"]} and CPLM give no methylation degree, so their rows are state `me` and do not meet UniProt\'s `me1`, `me2` and `me3` rows).', '']
    for m in MARKS:
        ov = overlap(tables[m])
        L += [f'**{m}** ({len(tables[m]):,} rows). Residues by source, any state: ' + ', '.join(f'{k} {v:,}' for k, v in sorted(ov['sites'].items())) + '.', ''] + counts_table(tables[m]) + ['']
    ov = overlap(cross)
    L += ['**methylarginine cross-check** (' + f'{len(cross):,} rows). Residues by source, any state: ' + ', '.join(f'{k} {v:,}' for k, v in sorted(ov['sites'].items())) + '.', ''] + counts_table(cross) + ['']
    arg = cross.drop_duplicates(['canonical_UniProtAC', 'position'])['in_methylarginine_release'].value_counts().sort_index()
    L += ['Arginines of the cross-check in the methylarginine data: ' + ', '.join(f'{k} {v:,}' for k, v in arg.items()) + '.', '']
    unip_only = {m: int(tables[m]['iptmnet_sources'].eq('unip').sum()) for m in MARKS}
    cp = tables['methyllysine']
    cp = cp[cp['sources'].str.contains('CPLM')]
    L += [f'iPTMnet carries a copy of UniProt: {unip_only["methyllysine"]:,} methyllysine and {unip_only["phosphorylation"]:,} phosphorylation rows have `unip` as their only upstream database, '
          'so "UniProt + iPTMnet" counts them as agreement of two sources when they are one.', '',
          f'CPLM 4.0 evidence of the {len(cp):,} methyllysine rows it reports (`cplm_evidence`): ' + ', '.join(f'{k} {v:,}' for k, v in cp['cplm_evidence'].value_counts().sort_index().items())
          + '. CPLM compiles records from other databases too (`Dat.`), so a residue that CPLM and UniProt or iPTMnet both report may rest on the same original record; the tables only record that the residue is shared.', '',
          '### Dropped in staging', '',
          'Each count is of source records. UniProt: ' + ', '.join(f'{v:,} {k.replace("_", " ")}' for k, v in sorted(summaries['uniprot']['dropped'].items())) +
          f'; the entries\' sequences were compared with the reference (UniProt {summaries["uniprot"]["uniprot_release_served"]} against 2026_03). iPTMnet: ' +
          ', '.join(f'{v:,} {k.replace("_", " ")}' for k, v in sorted(summaries['iptmnet']['dropped'].items())) +
          f'. CPLM: {summaries["cplm"]["records_in_download"]:,} rows in the download, {summaries["cplm"]["records_in_scope"]:,} human methylation rows, of which ' +
          ', '.join(f'{v:,} {k.replace("_", " ")}' for k, v in sorted(summaries['cplm']['dropped'].items())) + ' (nothing remapped).', '',
          '### How the context marks were made', '',
          f'`python scripts/54_stage_uniprot_ptm_sites.py`, `python scripts/55_stage_iptmnet_ptm_sites.py` and `python scripts/39_stage_cplm4_lys_methylation.py --marks-dir inputs/marks` (inputs in `inputs/marks/`, downloads listed in `inputs/marks/raw_files.tsv`), then '
          f'`python scripts/56_build_context_marks.py --version {version} --release-date {release_date}` at commit `{commit}` of https://github.com/Shechterlab/PTM_Pipeline.', '',
          '### Cite', '', 'The paper below, and the sources you use:', '']
    L += [f'- {s}: {SOURCES[s]["citation"]}' for s in SOURCE_ORDER]
    L += ['', f'## Methylarginine data (unchanged from {CARRIED_FROM})', '',
          f'What follows is the README of {CARRIED_FROM}, with its headings one level down. The files it names are in this folder with `_v1.5.0` in place of `_v1.4.3`.', '']
    old = re.sub(rf'_{CARRIED_FROM}_to_', f'_{version}_to_', old)
    old = old.replace(f'_v1.4.2_to_{CARRIED_FROM}', f'_{CARRIED_FROM}_to_{version}')
    old = re.sub(rf'_{re.escape(CARRIED_FROM)}(?=\.(?:tsv|xlsx))', f'_{version}', old)
    (out / 'README.md').write_text('\n'.join(L) + '\n' + demote(old), encoding='utf-8')
    return files


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--version', required=True)
    ap.add_argument('--release-date', required=True)
    ap.add_argument('--out', help='default releases/<version>')
    ap.add_argument('--inputs', default='inputs/marks')
    ap.add_argument('--carry-from', default=f'releases/{CARRIED_FROM}')
    args = ap.parse_args()
    inputs, out, src = Path(args.inputs), Path(args.out or f'releases/{args.version}'), Path(args.carry_from)
    out.mkdir(parents=True, exist_ok=True)
    version = args.version

    carried, edited = {}, []
    for path in sorted(p for p in src.rglob('*') if p.is_file()):
        rel = path.relative_to(src)
        if len(rel.parts) == 1 and rel.name in NOT_CARRIED:
            continue
        if rel.name.startswith('site_changes_'):
            continue     # the changes since the release before are rewritten below
        data = path.read_bytes()
        if rel.as_posix() in VERSIONED_TEXT:
            data = data.replace(CARRIED_FROM.encode(), version.encode())
        elif rel.suffix == '.xlsx':
            data = retitle_workbook(data, version)
        new_rel = renamed(rel.as_posix(), version)
        dest = out / new_rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not dest.exists() or dest.read_bytes() != data:
            dest.write_bytes(data)
        carried[new_rel] = {'from': rel.as_posix(), 'from_md5': md5_file(path), 'md5': md5_file(dest)}
        if carried[new_rel]['from_md5'] != carried[new_rel]['md5']:
            edited.append(new_rel)
    # every site is kept: the changes table in the format of the pipeline's
    sites = pd.read_csv(src / f'human_methylarginine_sites_{CARRIED_FROM}.tsv', sep='\t', dtype=str, keep_default_na=False,
                        usecols=['site_id', 'gene', 'uniprot_accession', 'position'])
    prev_changes = pd.read_csv(src / 'supporting' / f'site_changes_v1.4.2_to_{CARRIED_FROM}.tsv', sep='\t', dtype=str,
                               keep_default_na=False, usecols=['site_id', 'gene'])
    changes_name = f'supporting/site_changes_{CARRIED_FROM}_to_{version}.tsv'
    site_changes(sites, prev_changes).to_csv(out / changes_name, sep='\t', index=False, lineterminator='\n')

    tables, used = {}, {}
    for mark in MARKS + ['methylarginine']:
        suffix = 'crosscheck' if mark == 'methylarginine' else 'sites'
        paths = sorted(inputs.glob(f'*_{mark}_{suffix}.tsv'))
        if not paths:
            raise SystemExit(f'no staged table for {mark} in {inputs}')
        frames = [read_stage(p) for p in paths]
        unknown = {s for f in frames for s in f['source']} - set(SOURCES)
        if unknown:
            raise SystemExit(f'source {sorted(unknown)} has no license entry in SOURCES: check its terms first')
        tables[mark] = merge(frames)
        used[mark] = {p.name: md5_file(p) for p in paths}
    cross = crosscheck(tables.pop('methylarginine'), out, version)

    summaries = {'uniprot': json.loads((inputs / 'uniprot_ptm_sites_summary.json').read_text(encoding='utf-8')),
                 'iptmnet': json.loads((inputs / 'iptmnet_ptm_sites_summary.json').read_text(encoding='utf-8')),
                 'cplm': json.loads((inputs / 'cplm_ptm_sites_summary.json').read_text(encoding='utf-8'))}
    for m in MARKS:
        tables[m].to_csv(out / f'human_{m}_sites_{version}.tsv', sep='\t', index=False, lineterminator='\n')
    cross.to_csv(out / f'human_methylarginine_crosscheck_{version}.tsv', sep='\t', index=False, lineterminator='\n')

    # data dictionary: v1.4.3's, then the new tables
    new = {**{m: f'human_{m}_sites_{version}.tsv' for m in MARKS}, 'methylarginine': f'human_methylarginine_crosscheck_{version}.tsv'}
    lines = []
    for line in (src / 'data_dictionary.tsv').read_text(encoding='utf-8').rstrip('\n').split('\n'):
        name, rest = line.split('\t', 1)
        name = name.replace(f'_v1.4.2_to_{CARRIED_FROM}', f'_{CARRIED_FROM}_to_{version}')
        lines.append(renamed(name, version) + '\t' + rest)
    for name in new.values():
        for column, kind, description in DICTIONARY:
            if column == 'in_methylarginine_release' and 'crosscheck' not in name:
                continue
            if column == 'cplm_evidence' and 'methyllysine' not in name:
                continue
            if column in ('evidence_class', 'state', 'mark', 'sources'):
                values = CODES[column]
            else:
                values = ''
            lines.append(f'{name}\t{column}\t{kind}\t{values}\t\t{description}')
    (out / 'data_dictionary.tsv').write_text('\n'.join(lines) + '\n', encoding='utf-8')

    old_readme = src / 'README.md'
    (out / f'README_{CARRIED_FROM}.tmp').write_text(old_readme.read_text(encoding='utf-8'), encoding='utf-8')
    commit = code_commit()
    write_docs(out, version, args.release_date, tables, cross, summaries, commit)
    (out / f'README_{CARRIED_FROM}.tmp').unlink()

    counts = {m: {'rows': len(t), 'residues': int(t[['canonical_UniProtAC', 'position']].drop_duplicates().shape[0]),
                  'proteins': int(t['canonical_UniProtAC'].nunique()), 'by_state': {k: int(v) for k, v in t['state'].value_counts().sort_index().items()},
                  'by_evidence_class': {k: int(v) for k, v in t['evidence_class'].value_counts().sort_index().items()}, 'by_source': overlap(t)}
              for m, t in {**tables, 'methylarginine_crosscheck': cross}.items()}
    changes = [f'# Changes from {CARRIED_FROM} to {version}', '',
               f'{version} adds tables of context marks next to the methylarginine data. **The methylarginine data is unchanged from {CARRIED_FROM}**: '
               f'its {len(carried)} files are the v1.4.3 files under this release\'s names, with the data not rebuilt (`release.json`, `carried_over`, gives the v1.4.3 name and MD5 of each), because the build '
               '(`hpc/build_release.sh`) needs R with clusterProfiler and DOSE and the network, and no methylarginine input changed. Sites, proteins, '
               'records, statistics and layers are those of v1.4.3.', '',
               '### Files renamed', '',
               f'Every file name carries its release\'s version, so each of the {len(carried)} carried files is renamed `_{CARRIED_FROM}` to `_{version}` (for example '
               f'`human_methylarginine_sites_{CARRIED_FROM}.tsv` is now `human_methylarginine_sites_{version}.tsv`), its content unchanged except where it names its version: '
               + ', '.join(f'`{n}`' for n in edited) + f'. `supporting/site_changes_{CARRIED_FROM}_to_{version}.tsv` replaces v1.4.3\'s (v1.4.2 to v1.4.3); it lists every one of the '
               f'{len(sites):,} sites as kept, none added or removed, with its records the same.', '',
               '| | v1.4.3 | ' + version + ' |', '|---|---:|---:|']
    rel_prev = json.loads((src / 'release.json').read_text(encoding='utf-8'))['counts']
    for label, key in (('Sites', 'sites'), ('Proteins', 'proteins'), ('Sites kept pending', 'pending'), ('Source records', 'records'), ('Studies and databases', 'studies')):
        changes.append(f'| {label} | {rel_prev[key]:,} | {rel_prev[key]:,} |')
    changes += ['', '## New: context marks', '',
                'Versioned tables of other marks, for sites of those marks to read (`README.md` has the sources, licenses, counts and how to cite). '
                'The sources are UniProt and iPTMnet, which are openly licensed, and, for methyllysine only, CPLM 4.0, which is not (below).', '']
    for m in MARKS:
        c = counts[m]
        changes.append(f'- `human_{m}_sites_{version}.tsv`: {c["rows"]:,} rows on {c["residues"]:,} residues of {c["proteins"]:,} proteins; by state '
                       + ', '.join(f'{k} {v:,}' for k, v in c['by_state'].items()) + '; by evidence class ' + ', '.join(f'{k} {v:,}' for k, v in c['by_evidence_class'].items())
                       + '; residues by source ' + ', '.join(f'{k} {v:,}' for k, v in sorted(c['by_source']['sites'].items())) + '.')
    c = counts['methylarginine_crosscheck']
    changes.append(f'- `human_methylarginine_crosscheck_{version}.tsv`: {c["rows"]:,} arginine methylation rows of UniProt and iPTMnet with whether the methylarginine data has the arginine; residues by source '
                   + ', '.join(f'{k} {v:,}' for k, v in sorted(c['by_source']['sites'].items())) + '.')
    cc = summaries['cplm']
    changes += ['', 'Sources: UniProt (release ' + ', '.join(sorted({summaries["uniprot"]["uniprot_release_served"]})) + ', CC BY 4.0), iPTMnet (release ' + summaries['iptmnet']['iptmnet_release'] +
                ', CC BY-NC-SA 4.0) and CPLM 4.0 (below). The release stays CC BY-NC-SA 4.0. iPTMnet\'s download readme also mentions CC BY 4.0; the stricter terms are followed.', '',
                '### CPLM 4.0: license note', '',
                f'The methyllysine table includes all {cc["rows"]["methyllysine"]:,} validated human lysine methylation sites of CPLM 4.0 ({cc["proteins"]:,} proteins; Zhang et al., Nucleic Acids Res 50:D451-D459, 2022, '
                f'doi:10.1093/nar/gkab849; `Methylation.txt` of {cc["cplm_release"].split("(")[1].rstrip(")")}). **CPLM states no license**: cplm.biocuckoo.cn carries only a copyright notice ("All Rights Reserved"). '
                'The sites are included by the owner\'s decision of 2026-11-09 (confirmed 2026-11-10), cited, **pending the authors\' permission**; the terms the authors give, when they answer, are followed. The CC BY-NC-SA 4.0 of the release does not grant rights in CPLM\'s records beyond what the authors allow. '
                'Each CPLM row has `CPLM` in `sources` and CPLM\'s own evidence label in `cplm_evidence` (' + ', '.join(f'{k} {v:,}' for k, v in cc['rows_by_cplm_evidence'].items()) +
                '), so a reader can filter on it or leave CPLM out. CPLM rows are state `me` (CPLM gives no degree); `evidence_class` is experimental for Exp. and Exp./Dat., and unspecified for Dat. '
                f'(a record CPLM took from another database). `pmids` holds the references that are PubMed IDs ({cc["rows_with_pubmed_ids"]:,} rows; the rest of CPLM\'s reference field is identifiers of other databases). '
                f'Dropped on validation, nothing remapped: ' + ', '.join(f'{v:,} {k.replace("_", " ")}' for k, v in sorted(cc['dropped'].items())) + '. '
                'No lysine acylation is part of this release.', '',
                'Code: `scripts/39_stage_cplm4_lys_methylation.py`, `scripts/54_stage_uniprot_ptm_sites.py`, `scripts/55_stage_iptmnet_ptm_sites.py`, `scripts/56_build_context_marks.py`; staged tables in `inputs/marks/`.']
    (out / 'CHANGES.md').write_text('\n'.join(changes) + '\n', encoding='utf-8')

    release = {'title': 'Human methylarginine sites and context marks', 'version': version, 'previous_version': CARRIED_FROM,
               'release_date': args.release_date, 'repository': 'https://github.com/Shechterlab/PTM_Pipeline',
               'license': 'CC BY-NC-SA 4.0 (LICENSE.md); sources keep their own terms (inputs/MANIFEST.tsv); CPLM 4.0 states no license and is included pending the authors\' permission',
               'methylarginine_data': f'unchanged from {CARRIED_FROM}: files carried over byte for byte, not rebuilt',
               'methylarginine_code_commit': json.loads((src / 'release.json').read_text(encoding='utf-8')).get('code_commit'),
               'context_marks_code_commit': commit, 'carried_over': carried, 'carried_over_text_edited': edited,
               'site_changes': changes_name,
               'sources': {s: {'license': SOURCES[s]['license']} for s in SOURCE_ORDER},
               'inputs_md5': used, 'context_marks': counts}
    (out / 'release.json').write_text(json.dumps(release, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    names = sorted(p.relative_to(out).as_posix() for p in out.rglob('*') if p.is_file() and p.name != 'MD5SUMS')
    (out / 'MD5SUMS').write_text(''.join(f'{md5_file(out / n)}  {n}\n' for n in names), encoding='utf-8')
    print(json.dumps({m: {k: v for k, v in c.items() if k != 'residues' and k != 'proteins'} for m, c in counts.items()}, indent=1))


if __name__ == '__main__':
    main()

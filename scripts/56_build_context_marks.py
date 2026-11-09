#!/usr/bin/env python3
"""Write releases/v1.5.0: the methylarginine data of v1.4.3, carried over unchanged, and the context-mark tables.

    python scripts/56_build_context_marks.py --version v1.5.0 --release-date 2026-11-09

The methylarginine tables of v1.5.0 are the v1.4.3 files, copied byte for byte under their own names (the build,
hpc/build_release.sh, needs R and the network; this release changes no methylarginine data). Beside them, for each
context mark (methyllysine, phosphorylation), `human_<mark>_sites_<version>.tsv` merges the staged tables
`inputs/marks/*_<mark>_sites.tsv` of every source (UniProt, iPTMnet), one row per residue and state, and
`human_methylarginine_crosscheck_<version>.tsv` merges the sources' methylarginine records with a column that says whether
the arginine is in the methylarginine data. Also written: data_dictionary.tsv (v1.4.3's, with the new tables), README.md,
CHANGES.md, release.json and MD5SUMS. Nothing in the output depends on the date but --release-date, and the commit named
is the last one that changed scripts 54 to 56, so a second run gives the same bytes.

Merged columns:
  canonical_UniProtAC, position, residue, mark, state
  sources          the sources that report the residue and state, in order: UniProt; iPTMnet
  evidence_class   the best class of any source: experimental, then text-mined, then by-similarity, then unspecified
  evidence         UniProt's ECO codes ('' where only iPTMnet reports it)
  iptmnet_sources  iPTMnet's upstream databases for the residue (hprd, unip, rlim, rlim+, sign, npro, iedb, pro, intact)
  pmids            the PubMed IDs of all sources
  source_releases  the release of each source in `sources`, same order
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path

import pandas as pd

from common import md5_file

MARKS = ['methyllysine', 'phosphorylation']
SOURCE_ORDER = ['UniProt', 'iPTMnet']
CARRIED_FROM = 'v1.4.3'
NOT_CARRIED = {'README.md', 'CHANGES.md', 'MD5SUMS', 'data_dictionary.tsv', 'release.json'}
CLASS_RANK = ['experimental', 'text-mined', 'by-similarity', 'unspecified']
OUT_COLUMNS = ['canonical_UniProtAC', 'position', 'residue', 'mark', 'state', 'sources', 'evidence_class', 'evidence',
               'iptmnet_sources', 'pmids', 'source_releases']
SOURCES = {
    'UniProt': {'license': 'CC BY 4.0', 'url': 'https://www.uniprot.org/help/license',
                'citation': 'UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. Nucleic Acids Res '
                            '53:D609-D617 (2025). doi:10.1093/nar/gkae1010'},
    'iPTMnet': {'license': 'CC BY-NC-SA 4.0', 'url': 'https://research.bioinformatics.udel.edu/iptmnet/license',
                'citation': 'Huang H, Arighi CN, Ross KE, et al. iPTMnet: an integrated resource for protein post-translational '
                            'modification network discovery. Nucleic Acids Res 46:D542-D550 (2018). doi:10.1093/nar/gkx1104'},
    # Berryhill2023 is added only once the license of its supplements is confirmed (inputs/marks/README.md).
}
UPSTREAM = {'hprd': 'HPRD', 'unip': 'UniProt', 'rlim': 'RLIMS-P text mining', 'rlim+': 'RLIMS-P text mining', 'sign': 'SIGNOR',
            'npro': 'neXtProt', 'iedb': 'IEDB', 'pro': 'Protein Ontology', 'intact': 'IntAct'}
CODES = {
    'state': 'methyllysine: me1, me2, me3, me (degree not stated); phosphorylation: phospho; methylarginine: me1, me2a, me2s, me2 (type not stated), me (degree not stated)',
    'evidence_class': 'experimental, text-mined, by-similarity, unspecified',
    'mark': 'methyllysine, phosphorylation, methylarginine',
    'sources': 'UniProt, iPTMnet',
}
DICTIONARY = [
    ('canonical_UniProtAC', 'string', 'reviewed UniProtKB accession of the canonical sequence in the reference proteome (UniProt 2026_03, UP000005640)'),
    ('position', 'integer', 'residue number on that sequence, 1-based'),
    ('residue', 'code', 'one-letter code of the residue, checked against the reference sequence (K for methyllysine, S, T or Y for phosphorylation, R for methylarginine)'),
    ('mark', 'code', 'the modification'),
    ('state', 'code', 'the state reported for the residue; a residue with several states has one row per state. iPTMnet does not give a methylation degree (state me)'),
    ('sources', 'list<code>', 'the sources that report the residue and state, separated by "; "'),
    ('evidence_class', 'code', 'the best class of any source. experimental: UniProt ECO:0000269 or ECO:0007744, or an iPTMnet record curated from an experiment (upstream hprd, sign, npro, pro, intact, iedb). text-mined: only RLIMS-P reports it. by-similarity: UniProt records with no experimental code (mainly ECO:0000250). unspecified: only iPTMnet\'s copy of UniProt (upstream unip), which has no ECO code. A site can hide by-similarity and unspecified rows by filtering on this column'),
    ('evidence', 'list<string>', 'UniProt ECO evidence codes, separated by ";" (empty where only iPTMnet reports the residue)'),
    ('iptmnet_sources', 'list<code>', 'the upstream databases behind the iPTMnet record, separated by ";": hprd (HPRD), unip (UniProt), rlim and rlim+ (RLIMS-P text mining), sign (SIGNOR), npro (neXtProt), iedb (IEDB), pro (Protein Ontology), intact (IntAct); iPTMnet asks that they are cited'),
    ('pmids', 'list<string>', 'PubMed IDs of the records, separated by ";", when given'),
    ('source_releases', 'list<string>', 'the release of each source in sources, same order (UniProt release, iPTMnet release)'),
    ('in_methylarginine_release', 'code', f'crosscheck table only: whether the arginine is a site, or a lower-confidence site, in the methylarginine data of this release (unchanged from {CARRIED_FROM}), or not at all. Codes: site, lower_confidence, no'),
]


def read_stage(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    if 'upstream_sources' not in frame:
        frame['upstream_sources'] = ''
    return frame


def best_class(classes) -> str:
    return min(set(classes), key=CLASS_RANK.index)


def union(grp: pd.DataFrame, col: str) -> list[str]:
    return sorted({v for x in grp[col] for v in x.split(';') if v})


def merge(frames: list[pd.DataFrame]) -> pd.DataFrame:
    """One row per residue and state from the staged rows of every source."""
    table = pd.concat(frames, ignore_index=True)
    rows = []
    for key, grp in table.groupby(['canonical_UniProtAC', 'position', 'residue', 'mark', 'state'], sort=False):
        present = [s for s in SOURCE_ORDER if s in set(grp['source'])]
        release = {s: grp.loc[grp['source'] == s, 'source_release'].iloc[0] for s in present}
        rows.append({'canonical_UniProtAC': key[0], 'position': key[1], 'residue': key[2], 'mark': key[3], 'state': key[4],
                     'sources': '; '.join(present), 'evidence_class': best_class(grp['evidence_class']),
                     'evidence': ';'.join(union(grp, 'evidence')), 'iptmnet_sources': ';'.join(union(grp, 'upstream_sources')),
                     'pmids': ';'.join(sorted(union(grp, 'pmids'), key=int)), 'source_releases': '; '.join(release[s] for s in present)})
    out = pd.DataFrame(rows, columns=OUT_COLUMNS)
    out['_pos'] = out['position'].astype(int)
    return out.sort_values(['canonical_UniProtAC', '_pos', 'state']).drop(columns='_pos').reset_index(drop=True)


def crosscheck(merged: pd.DataFrame, release_dir: Path) -> pd.DataFrame:
    def keys(name: str) -> set[tuple[str, str]]:
        t = pd.read_csv(release_dir / name, sep='\t', dtype=str, keep_default_na=False, usecols=['uniprot_accession', 'position'])
        return set(zip(t['uniprot_accession'], t['position'], strict=True))
    sites = keys(f'human_methylarginine_sites_{CARRIED_FROM}.tsv')
    lower = keys(f'human_methylarginine_lower_confidence_sites_{CARRIED_FROM}.tsv')
    out = merged.copy()
    out['in_methylarginine_release'] = ['site' if k in sites else 'lower_confidence' if k in lower else 'no'
                                        for k in zip(out['canonical_UniProtAC'], out['position'], strict=True)]
    return out


def overlap(table: pd.DataFrame) -> dict:
    """Rows and residues (any state) by which sources report them."""
    def label(sources):
        s = set(sources.split('; '))
        return 'both' if len(s) == 2 else f'{next(iter(s))} only'
    rows = table['sources'].map(label).value_counts()
    sites = table.groupby(['canonical_UniProtAC', 'position'])['sources'].agg(lambda v: '; '.join(sorted({x for s in v for x in s.split('; ')})))
    sites = sites.map(lambda s: 'both' if '; ' in s else f'{s} only').value_counts()
    return {'rows': {k: int(v) for k, v in rows.items()}, 'sites': {k: int(v) for k, v in sites.items()}}


def counts_table(table: pd.DataFrame) -> list[str]:
    classes = [c for c in CLASS_RANK if c in set(table['evidence_class'])]
    lines = ['| State | ' + ' | '.join(classes) + ' | UniProt only | iPTMnet only | both | Rows |', '|---|' + '---:|' * (len(classes) + 4)]
    def line(name, grp):
        by = grp['sources'].map(lambda s: 'both' if '; ' in s else f'{s} only').value_counts()
        return (f'| {name} | ' + ' | '.join(f'{(grp["evidence_class"] == c).sum():,}' for c in classes)
                + f' | {by.get("UniProt only", 0):,} | {by.get("iPTMnet only", 0):,} | {by.get("both", 0):,} | {len(grp):,} |')
    for state, grp in table.groupby('state'):
        lines.append(line(state, grp))
    lines.append(line('all', table))
    return lines


def code_commit() -> str:
    paths = ['scripts/54_stage_uniprot_ptm_sites.py', 'scripts/55_stage_iptmnet_ptm_sites.py', 'scripts/56_build_context_marks.py']
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
         f'{version} keeps the methylarginine data of {CARRIED_FROM} unchanged (every file below under `human_methylarginine_*`, '
         '`statistics_*`, `studies_*`, `supporting/` and `analyses*` is the v1.4.3 file, byte for byte, under its own name) and adds '
         'tables of context marks: the lysine methylation and phosphorylation sites of the human reviewed proteome, and the '
         'arginine methylation records of UniProt and iPTMnet as a cross-check. Sites for each mark read their own table and show the '
         'others as context. Positions of the new tables are on UniProt release 2026_03, the reference of the methylarginine data.', '',
         '## Context marks', '', '| File | One row per |', '|---|---|']
    for m in MARKS:
        t = tables[m]
        L.append(f'| `{files[m]}` | {m} residue and state ({len(t):,} rows on {t[["canonical_UniProtAC", "position"]].drop_duplicates().shape[0]:,} '
                 f'residues of {t["canonical_UniProtAC"].nunique():,} proteins) |')
    L += [f'| `{files["methylarginine"]}` | arginine methylation record of UniProt or iPTMnet ({len(cross):,} rows), with whether the methylarginine data has the arginine. A cross-check, not a second methylarginine table |', '',
          'Columns are in `data_dictionary.tsv`. Each row says which sources report it (`sources`), the best `evidence_class`, UniProt\'s ECO codes, '
          'iPTMnet\'s upstream databases (`iptmnet_sources`) and the PubMed IDs. A site can hide `by-similarity` (UniProt ECO:0000250 and other '
          'non-experimental records) and `unspecified` rows by filtering on `evidence_class`.', '',
          '### Sources and licenses', '', '| Source | Release | License | Terms |', '|---|---|---|---|']
    for s in SOURCE_ORDER:
        L.append(f'| {s} | {", ".join(releases[s])} | {SOURCES[s]["license"]} | {SOURCES[s]["url"]} |')
    L += ['', 'The release is CC BY-NC-SA 4.0 ([LICENSE.md](../../LICENSE.md)). UniProt (CC BY 4.0) is compatible with it. iPTMnet\'s license page says '
          'CC BY-NC-SA 4.0; its download readme also mentions CC BY 4.0 in its last line. The release follows the stricter terms: share-alike and '
          'non-commercial. iPTMnet asks that the upstream databases of its `source` column are cited when its data are redistributed: each row keeps them '
          '(`iptmnet_sources`), and they are HPRD, UniProt, RLIMS-P text mining, SIGNOR, neXtProt, IEDB, Protein Ontology and IntAct (terms are theirs). '
          'Not used: CPLM (all rights reserved), PhosphoSitePlus, dbPTM and qPTM. The supplements of Berryhill et al. 2023 are not used: they state no license of their own.', '',
          '### Counts', '',
          f'Rows by state, evidence class and source (rows reported by UniProt only, iPTMnet only, or both; iPTMnet {summaries["iptmnet"]["iptmnet_release"]} gives no methylation degree, so its rows are state `me` and do not meet UniProt\'s `me1`, `me2` and `me3` rows).', '']
    for m in MARKS:
        ov = overlap(tables[m])
        L += [f'**{m}** ({len(tables[m]):,} rows). Residues by source, any state: ' + ', '.join(f'{k} {v:,}' for k, v in sorted(ov['sites'].items())) + '.', ''] + counts_table(tables[m]) + ['']
    ov = overlap(cross)
    L += ['**methylarginine cross-check** (' + f'{len(cross):,} rows). Residues by source, any state: ' + ', '.join(f'{k} {v:,}' for k, v in sorted(ov['sites'].items())) + '.', ''] + counts_table(cross) + ['']
    arg = cross.drop_duplicates(['canonical_UniProtAC', 'position'])['in_methylarginine_release'].value_counts().sort_index()
    L += ['Arginines of the cross-check in the methylarginine data: ' + ', '.join(f'{k} {v:,}' for k, v in arg.items()) + '.', '']
    unip_only = {m: int(tables[m]['iptmnet_sources'].eq('unip').sum()) for m in MARKS}
    L += [f'iPTMnet carries a copy of UniProt: {unip_only["methyllysine"]:,} methyllysine and {unip_only["phosphorylation"]:,} phosphorylation rows have `unip` as their only upstream database, '
          'so "both" counts them as agreement of two sources when they are one.', '',
          '### Dropped in staging', '',
          'Each count is of source records. UniProt: ' + ', '.join(f'{v:,} {k.replace("_", " ")}' for k, v in sorted(summaries['uniprot']['dropped'].items())) +
          f'; the entries\' sequences were compared with the reference (UniProt {summaries["uniprot"]["uniprot_release_served"]} against 2026_03). iPTMnet: ' +
          ', '.join(f'{v:,} {k.replace("_", " ")}' for k, v in sorted(summaries['iptmnet']['dropped'].items())) + '.', '',
          '### How the context marks were made', '',
          f'`python scripts/54_stage_uniprot_ptm_sites.py` and `python scripts/55_stage_iptmnet_ptm_sites.py` (inputs in `inputs/marks/`, downloads listed in `inputs/marks/raw_files.tsv`), then '
          f'`python scripts/56_build_context_marks.py --version {version} --release-date {release_date}` at commit `{commit}` of https://github.com/Shechterlab/PTM_Pipeline.', '',
          '### Cite', '', 'The paper below, and the sources you use:', '']
    L += [f'- {s}: {SOURCES[s]["citation"]}' for s in SOURCE_ORDER]
    L += ['', f'## Methylarginine data (unchanged from {CARRIED_FROM})', '',
          f'What follows is the README of {CARRIED_FROM}, with its headings one level down. The files it names are in this folder under those names.', '']
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

    carried = {}
    for path in sorted(p for p in src.rglob('*') if p.is_file()):
        rel = path.relative_to(src)
        if len(rel.parts) == 1 and rel.name in NOT_CARRIED:
            continue
        dest = out / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not dest.exists() or md5_file(dest) != md5_file(path):
            shutil.copyfile(path, dest)
        carried[rel.as_posix()] = md5_file(dest)

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
    cross = crosscheck(tables.pop('methylarginine'), out)

    summaries = {'uniprot': json.loads((inputs / 'uniprot_ptm_sites_summary.json').read_text(encoding='utf-8')),
                 'iptmnet': json.loads((inputs / 'iptmnet_ptm_sites_summary.json').read_text(encoding='utf-8'))}
    for m in MARKS:
        tables[m].to_csv(out / f'human_{m}_sites_{version}.tsv', sep='\t', index=False, lineterminator='\n')
    cross.to_csv(out / f'human_methylarginine_crosscheck_{version}.tsv', sep='\t', index=False, lineterminator='\n')

    # data dictionary: v1.4.3's, then the new tables
    new = {**{m: f'human_{m}_sites_{version}.tsv' for m in MARKS}, 'methylarginine': f'human_methylarginine_crosscheck_{version}.tsv'}
    lines = (src / 'data_dictionary.tsv').read_text(encoding='utf-8').rstrip('\n').split('\n')
    for name in new.values():
        for column, kind, description in DICTIONARY:
            if column == 'in_methylarginine_release' and 'crosscheck' not in name:
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
               f'its {len(carried)} files are copied byte for byte (their MD5 are in `release.json`, `carried_over`), not rebuilt, because the build '
               '(`hpc/build_release.sh`) needs R with clusterProfiler and DOSE and the network, and no methylarginine input changed. Sites, proteins, '
               'records, statistics and layers are those of v1.4.3.', '',
               '| | v1.4.3 | ' + version + ' |', '|---|---:|---:|']
    rel_prev = json.loads((src / 'release.json').read_text(encoding='utf-8'))['counts']
    for label, key in (('Sites', 'sites'), ('Proteins', 'proteins'), ('Sites kept pending', 'pending'), ('Source records', 'records'), ('Studies and databases', 'studies')):
        changes.append(f'| {label} | {rel_prev[key]:,} | {rel_prev[key]:,} |')
    changes += ['', '## New: context marks', '',
                'Versioned tables of other marks, from openly licensed sources, for sites of those marks to read (`README.md` has the sources, licenses, counts and how to cite).', '']
    for m in MARKS:
        c = counts[m]
        changes.append(f'- `human_{m}_sites_{version}.tsv`: {c["rows"]:,} rows on {c["residues"]:,} residues of {c["proteins"]:,} proteins; by state '
                       + ', '.join(f'{k} {v:,}' for k, v in c['by_state'].items()) + '; by evidence class ' + ', '.join(f'{k} {v:,}' for k, v in c['by_evidence_class'].items())
                       + '; residues by source ' + ', '.join(f'{k} {v:,}' for k, v in sorted(c['by_source']['sites'].items())) + '.')
    c = counts['methylarginine_crosscheck']
    changes.append(f'- `human_methylarginine_crosscheck_{version}.tsv`: {c["rows"]:,} arginine methylation rows of UniProt and iPTMnet with whether the methylarginine data has the arginine; residues by source '
                   + ', '.join(f'{k} {v:,}' for k, v in sorted(c['by_source']['sites'].items())) + '.')
    changes += ['', 'Sources: UniProt (release ' + ', '.join(sorted({summaries["uniprot"]["uniprot_release_served"]})) + ', CC BY 4.0) and iPTMnet (release ' + summaries['iptmnet']['iptmnet_release'] +
                ', CC BY-NC-SA 4.0). The release stays CC BY-NC-SA 4.0. iPTMnet\'s download readme also mentions CC BY 4.0; the stricter terms are followed.', '',
                'Code: `scripts/54_stage_uniprot_ptm_sites.py`, `scripts/55_stage_iptmnet_ptm_sites.py`, `scripts/56_build_context_marks.py`; staged tables in `inputs/marks/`.']
    (out / 'CHANGES.md').write_text('\n'.join(changes) + '\n', encoding='utf-8')

    release = {'title': 'Human methylarginine sites and context marks', 'version': version, 'previous_version': CARRIED_FROM,
               'release_date': args.release_date, 'repository': 'https://github.com/Shechterlab/PTM_Pipeline',
               'license': 'CC BY-NC-SA 4.0 (LICENSE.md); sources keep their own terms (inputs/MANIFEST.tsv)',
               'methylarginine_data': f'unchanged from {CARRIED_FROM}: files carried over byte for byte, not rebuilt',
               'methylarginine_code_commit': json.loads((src / 'release.json').read_text(encoding='utf-8')).get('code_commit'),
               'context_marks_code_commit': commit, 'carried_over': carried,
               'sources': {s: {'license': SOURCES[s]['license']} for s in SOURCE_ORDER},
               'inputs_md5': used, 'context_marks': counts}
    (out / 'release.json').write_text(json.dumps(release, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    names = sorted(p.relative_to(out).as_posix() for p in out.rglob('*') if p.is_file() and p.name != 'MD5SUMS')
    (out / 'MD5SUMS').write_text(''.join(f'{md5_file(out / n)}  {n}\n' for n in names), encoding='utf-8')
    print(json.dumps({m: {k: v for k, v in c.items() if k != 'residues' and k != 'proteins'} for m, c in counts.items()}, indent=1))


if __name__ == '__main__':
    main()

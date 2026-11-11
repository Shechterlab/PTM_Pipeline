#!/usr/bin/env python3
"""Write releases/v1.5.0: the methylarginine data of v1.4.3, carried over unchanged, and the context-mark tables.

    python scripts/56_build_context_marks.py --version v1.5.0 --release-date 2026-11-09

The methylarginine tables of v1.5.0 are the v1.4.3 files, copied under this release's own names (_v1.4.3 becomes
_v1.5.0, as every release names all its files) with the data unchanged (the build, hpc/build_release.sh, needs R and the
network; this release changes no methylarginine data). Only the strings that name the version inside a file change, as when
v1.4.3 followed v1.4.2: the title and the statistics file name in the two analyses READMEs, the file names in their MD5SUMS,
the title of STUDIES.md and the title and sheet name of the workbook. supporting/site_changes_v1.4.3_to_v1.5.0.tsv is written
in the format of build_release.build_changes: every site kept, none added or removed. Beside them, for each
context mark (methyllysine, phosphorylation, lysine_acylation, and acetylation, the acetyl rows of the last; ubiquitination; sumoylation), `human_<mark>_sites_<version>.tsv` merges the staged tables
`inputs/marks/*_<mark>_sites.tsv` of every source (UniProt, iPTMnet, CPLM 4.0 and Berryhill 2023 for methyllysine, Ochoa 2020 for phosphorylation, UniProt, iPTMnet and CPLM 4.0 for ubiquitination and sumoylation), one row per residue and state
(acylation: acyl type, in `acyl_type` in place of `state`), and
`human_methylarginine_crosscheck_<version>.tsv` merges the sources' methylarginine records with a column that says whether
the arginine is in the methylarginine data. Also written: data_dictionary.tsv (v1.4.3's, with the new tables), README.md,
CHANGES.md, release.json and MD5SUMS. Nothing in the output depends on the date but --release-date, and the commit named
is the last one that changed scripts 39, 54, 55, 56, 57, 58 or 59, so a second run gives the same bytes.

Merged columns:
  canonical_UniProtAC, position, residue, mark, state
  sources          the sources that report the residue and state, in order: UniProt; iPTMnet; CPLM; Berryhill 2023; Ochoa 2020
  evidence_class   the best class of any source: experimental, then text-mined, then by-similarity, then unspecified
                   (CPLM: its labels Exp. and Exp./Dat. are experimental, Dat. is unspecified; Berryhill 2023 and Ochoa 2020: experimental, mass spectrometry)
  evidence         UniProt's ECO codes ('' where only iPTMnet reports it)
  iptmnet_sources  iPTMnet's upstream databases for the residue (hprd, unip, rlim, rlim+, sign, npro, iedb, pro, intact)
  pmids            the PubMed IDs of all sources
  source_releases  the release of each source in `sources`, same order
  cplm_evidence    (methyllysine, ubiquitination and sumoylation) CPLM's own evidence label for the residue: Exp., Exp./Dat. or Dat. ('' where CPLM does not report it)
  berryhill_site_probability, berryhill_cell_lines, berryhill_enrichment   (methyllysine only) from Berryhill 2023; '' on other rows
  functional_score, best_localization_probability, best_pep, biological_samples, spectral_counts   (phosphorylation only) from
                   Ochoa 2020; '' on other rows
  iptmnet_score    (methyllysine, phosphorylation, acylation, ubiquitination, sumoylation) iPTMnet's own score of the residue, 0 to 4 (scripts/59), on the rows whose
                   sources include iPTMnet and that score.txt has; '' otherwise

CPLM 4.0 states no license. The owner decided on 2026-11-09 to include its human lysine methylation sites, and on 2026-11-10 its
ubiquitination and SUMOylation sites on the same terms, cited, pending the authors' permission; their terms, when given, are followed. CPLM gives no
methylation degree (its methylation rows are state `me`); its ubiquitination rows are state `ub` and its SUMOylation rows state `sumo`.
Ubiquitination and sumoylation tables: states `ub` and `sumo`, from UniProt's cross-links to ubiquitin and to any SUMO paralog (the paralog is not kept), iPTMnet's
UBIQUITINATION and SUMOYLATION records and CPLM's ubiquitination and SUMOylation, all on lysine.
Berryhill 2023: article CC BY 4.0, supplements without a separate notice, included by the owner's decision of 2026-11-10. Ochoa 2020: article
not open access, supplements without a license, the derived table published by the owner's decision of 2026-11-10.
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

MARKS = ['methyllysine', 'phosphorylation', 'lysine_acylation', 'ubiquitination', 'sumoylation']
ACYL = {'lysine_acylation', 'acetylation'}          # their `state` column is written as `acyl_type`
# acetylation: the acetyl rows of lysine_acylation, as a table of its own
FILE_MARKS = ['methyllysine', 'phosphorylation', 'lysine_acylation', 'acetylation', 'ubiquitination', 'sumoylation']
SOURCE_ORDER = ['UniProt', 'iPTMnet', 'CPLM', 'Berryhill 2023', 'Ochoa 2020']
# columns a source adds to the tables it contributes to (empty on the rows of other sources), in this order after OUT_COLUMNS
SOURCE_COLUMNS = {'CPLM': ['cplm_evidence'],
                  'Berryhill 2023': ['berryhill_site_probability', 'berryhill_cell_lines', 'berryhill_enrichment'],
                  'Ochoa 2020': ['functional_score', 'best_localization_probability', 'best_pep', 'biological_samples', 'spectral_counts']}
EXTRA_COLUMNS = [c for cols in SOURCE_COLUMNS.values() for c in cols]
SCORE_COLUMN = 'iptmnet_score'    # last column of the context-mark tables, filled for the rows iPTMnet reports (scripts/59)
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
    'Berryhill 2023': {'license': 'article CC BY 4.0; supplements carry no separate notice; included by the owner\'s decision',
                       'url': 'https://www.nature.com/articles/s41598-022-27175-x',
                       'citation': 'Berryhill CA, Hanquier JN, Doud EH, Cordeiro-Spinetti E, Dickson BM, Rothbart SB, Mosley AL. Global lysine '
                                   'methylome profiling using systematically characterized affinity reagents. Sci Rep 13:377 (2023). '
                                   'doi:10.1038/s41598-022-27175-x'},
    'Ochoa 2020': {'license': 'article not open access (exclusive license to Springer Nature America); supplements carry no license; included by the owner\'s decision',
                   'url': 'https://www.nature.com/articles/s41587-019-0344-3',
                   'citation': 'Ochoa D, Jarnuczak AF, Viéitez C, et al. The functional landscape of the human phosphoproteome. '
                               'Nat Biotechnol 38:365-373 (2020). doi:10.1038/s41587-019-0344-3'},
}
UPSTREAM = {'hprd': 'HPRD', 'unip': 'UniProt', 'rlim': 'RLIMS-P text mining', 'rlim+': 'RLIMS-P text mining', 'sign': 'SIGNOR',
            'npro': 'neXtProt', 'iedb': 'IEDB', 'pro': 'Protein Ontology', 'intact': 'IntAct'}
CODES = {
    'state': 'methyllysine: me1, me2, me3, me (degree not stated); phosphorylation: phospho; ubiquitination: ub; sumoylation: sumo; methylarginine: me1, me2a, me2s, me2 (type not stated), me (degree not stated)',
    'evidence_class': 'experimental, text-mined, by-similarity, unspecified',
    'mark': 'methyllysine, phosphorylation, lysine_acylation, ubiquitination, sumoylation, methylarginine',
    'sources': ', '.join(SOURCE_ORDER),
    'cplm_evidence': 'Exp., Exp./Dat., Dat.',
}
DICTIONARY = [
    ('canonical_UniProtAC', 'string', 'reviewed UniProtKB accession of the canonical sequence in the reference proteome (UniProt 2026_03, UP000005640)'),
    ('position', 'integer', 'residue number on that sequence, 1-based'),
    ('residue', 'code', 'one-letter code of the residue, checked against the reference sequence (K for methyllysine, lysine acylation, ubiquitination and sumoylation, S, T or Y for phosphorylation, R for methylarginine)'),
    ('mark', 'code', 'the modification (the acylation tables have lysine_acylation)'),
    ('acyl_type', 'code', 'acylation tables only, in place of state: the acyl group on the lysine, as UniProt names it in its N6-<acyl>lysine terms (acetyl, succinyl, 2-hydroxyisobutyryl, beta-hydroxybutyryl, lactoyl, glutaryl, crotonyl, methacryl, isonicotinyl, benzoyl, butyryl, malonyl, propionyl); iPTMnet gives acetyl only. A residue with several acyl types has one row per type'),
    ('state', 'code', 'the state reported for the residue; a residue with several states has one row per state. iPTMnet and CPLM do not give a methylation degree (state me); UniProt and Berryhill 2023 give me1, me2 and me3. Ubiquitination rows are state ub and sumoylation rows state sumo (UniProt\'s SUMO paralog is not kept)'),
    ('sources', 'list<code>', 'the sources that report the residue and state, separated by "; "'),
    ('evidence_class', 'code', 'the best class of any source. experimental: UniProt ECO:0000269 or ECO:0007744, an iPTMnet record curated from an experiment (upstream hprd, sign, npro, pro, intact, iedb), a CPLM record labelled Exp. or Exp./Dat. (CPLM curates these from experimental papers), a Berryhill 2023 site (identified by mass spectrometry in HEK293T or U2OS cells) or an Ochoa 2020 site (identified by mass spectrometry). text-mined: only RLIMS-P reports it. by-similarity: UniProt records with no experimental code (mainly ECO:0000250). unspecified: only iPTMnet\'s copy of UniProt (upstream unip), which has no ECO code, or only a CPLM record labelled Dat. (taken from another database; CPLM names no experiment). A site can hide by-similarity and unspecified rows by filtering on this column'),
    ('evidence', 'list<string>', 'UniProt ECO evidence codes, separated by ";" (empty where only iPTMnet reports the residue)'),
    ('iptmnet_sources', 'list<code>', 'the upstream databases behind the iPTMnet record, separated by ";": hprd (HPRD), unip (UniProt), rlim and rlim+ (RLIMS-P text mining), sign (SIGNOR), npro (neXtProt), iedb (IEDB), pro (Protein Ontology), intact (IntAct); iPTMnet asks that they are cited'),
    ('pmids', 'list<string>', 'PubMed IDs of the records, separated by ";", when given. CPLM\'s reference field also holds identifiers that are not PubMed IDs; only seven- and eight-digit numbers up to 34,300,000 (the PubMed IDs assigned by the date of CPLM\'s file) are taken as PubMed IDs'),
    ('source_releases', 'list<string>', 'the release of each source in sources, same order (UniProt release, iPTMnet release, CPLM release and the date of its file)'),
    ('cplm_evidence', 'code', 'methyllysine, ubiquitination and sumoylation tables only: CPLM\'s own evidence label for the residue, as CPLM gives it. Exp.: from an experimental paper; Exp./Dat.: both an experimental paper and another database; Dat.: another database only. Empty where CPLM does not report the residue. CPLM 4.0 states no license; its methylation sites (owner\'s decision of 2026-11-09) and ubiquitination and SUMOylation sites (2026-11-10) are included pending the authors\' permission, cited to Zhang et al., Nucleic Acids Res 50:D451 (2022)'),
    ('berryhill_site_probability', 'number', 'methyllysine table only: the highest site-localization probability (percent; 90.1 to 100 in Supplementary Table 5 of Berryhill et al., Sci Rep 13:377, 2023) over the cell lines and enrichments that identified the residue and state; empty where the residue is not a Berryhill 2023 row, and where the paper gives none (its PTMScan Rep 2 rows)'),
    ('berryhill_cell_lines', 'list<string>', 'methyllysine table only: the cell lines (HEK293T, U2OS) in which Berryhill et al. 2023 identified the residue and state, separated by ";"; empty for other sources'),
    ('berryhill_enrichment', 'list<string>', 'methyllysine table only: the enrichments (the antibody or PTMScan reagents of Berryhill et al. 2023, or none for global proteomics) that identified the residue and state, separated by ";"; empty for other sources'),
    ('functional_score', 'number', 'phosphorylation table only: the functional score of Ochoa et al., Nat Biotechnol 38:365 (2020), Supplementary Table 3, 0 to 1 (higher means the phosphorylation site is more likely functional, as predicted from the paper\'s features); empty where the residue is not an Ochoa 2020 row'),
    ('best_localization_probability', 'number', 'phosphorylation table only: the best localization probability of the phosphorylation site over its observations in the data compiled by Ochoa et al. 2020 (MaxQuant), 0 to 1; empty for other sources'),
    ('best_pep', 'number', 'phosphorylation table only: the best posterior error probability of a peptide that carries the phosphorylation site (Ochoa et al. 2020; lower is better); empty for other sources'),
    ('biological_samples', 'integer', 'phosphorylation table only: the number of biological samples (experimental conditions) in which Ochoa et al. 2020 observed the phosphorylation site; empty for other sources'),
    ('spectral_counts', 'integer', 'phosphorylation table only: the total spectral counts of the phosphorylation site in the data of Ochoa et al. 2020; empty for other sources'),
    ('iptmnet_score', 'integer', 'methyllysine, phosphorylation, acylation, ubiquitination and sumoylation tables: iPTMnet 6.2\'s own score of the residue (score.txt, rows without an enzyme; 0 to 4). The iPTMnet readme says only that the score is calculated from all the evidence sources shown on its website. Given only on rows whose sources include iPTMnet, and empty where its score file has no row for the residue; it is iPTMnet\'s and does not take the other sources of this table into account'),
    ('in_methylarginine_release', 'code', f'crosscheck table only: whether the arginine is a site, or a lower-confidence site, in the methylarginine data of this release (unchanged from {CARRIED_FROM}), or not at all. Codes: site, lower_confidence, no'),
]


def read_stage(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    if 'upstream_sources' not in frame:
        frame['upstream_sources'] = ''
    for col in EXTRA_COLUMNS:
        if col not in frame:
            frame[col] = ''
    return frame


def best_class(classes) -> str:
    return min(set(classes), key=CLASS_RANK.index)


def union(grp: pd.DataFrame, col: str) -> list[str]:
    return sorted({v for x in grp[col] for v in x.split(';') if v})


def merge(frames: list[pd.DataFrame]) -> pd.DataFrame:
    """One row per residue and state from the staged rows of every source."""
    table = pd.concat(frames, ignore_index=True)
    for col in ['upstream_sources'] + EXTRA_COLUMNS:
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
                     **{c: ';'.join(union(grp, c)) for c in EXTRA_COLUMNS}})
    out = pd.DataFrame(rows, columns=OUT_COLUMNS + EXTRA_COLUMNS)
    out = out[OUT_COLUMNS + [c for src, cols in SOURCE_COLUMNS.items() if src in set(table['source']) for c in cols]]
    out['_pos'] = out['position'].astype(int)
    return out.sort_values(['canonical_UniProtAC', '_pos', 'state']).drop(columns='_pos').reset_index(drop=True)


def add_iptmnet_scores(table: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    """The table with `iptmnet_score` after its other columns: iPTMnet's score of the residue (scripts/59) on the rows whose
    sources include iPTMnet, empty on the others and where iPTMnet's score file has no row for the residue."""
    given = dict(zip(zip(scores['canonical_UniProtAC'], scores['position'], scores['mark'], strict=True), scores[SCORE_COLUMN], strict=True))
    out = table.copy()
    keys = zip(out['canonical_UniProtAC'], out['position'], out['mark'], strict=True)
    out[SCORE_COLUMN] = [given.get(k, '') if 'iPTMnet' in s.split('; ') else '' for k, s in zip(keys, out['sources'], strict=True)]
    return out


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
    paths = ['scripts/39_stage_cplm4_lys_methylation.py', 'scripts/54_stage_uniprot_ptm_sites.py', 'scripts/55_stage_iptmnet_ptm_sites.py',
             'scripts/57_stage_berryhill_lys_methylation.py', 'scripts/58_stage_ochoa_phosphosites.py', 'scripts/59_stage_iptmnet_scores.py',
             'scripts/56_build_context_marks.py']
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


def pct(n, d) -> str:
    return f'{n:,} of {d:,}'


def source_notes(tables: dict, summaries: dict) -> list[str]:
    """The paragraphs on Berryhill 2023, Ochoa 2020 and the iPTMnet score: license, citation, what the columns hold, counts."""
    b, o, sc = summaries['berryhill'], summaries['ochoa'], summaries['score']
    ml, ph = tables['methyllysine'], tables['phosphorylation']
    bh_rows, oc_rows = ml['sources'].str.contains('Berryhill 2023'), ph['sources'].str.contains('Ochoa 2020')
    lines = []
    for name in ('methyllysine', 'phosphorylation', 'lysine_acylation', 'ubiquitination', 'sumoylation'):
        t = tables[name]
        ip = t['sources'].str.split('; ').map(lambda x: 'iPTMnet' in x)
        lines.append(f'{name.replace("_", " ")} {pct(int((t.loc[ip, SCORE_COLUMN] != "").sum()), int(ip.sum()))}')
    return [
        '**Berryhill 2023** (Berryhill CA, Hanquier JN, Doud EH, Cordeiro-Spinetti E, Dickson BM, Rothbart SB, Mosley AL. "Global lysine methylome profiling '
        'using systematically characterized affinity reagents", Sci Rep 13:377, 2023, doi:10.1038/s41598-022-27175-x). The article is CC BY 4.0; its supplements carry no separate '
        'license notice. Its Supplementary Table 5 is included by the owner\'s decision of 2026-11-10, without waiting for the authors\' confirmation, cited. '
        f'It contributes {int(bh_rows.sum()):,} methyllysine rows (states me1, me2, me3 kept; `evidence_class` experimental, as the sites are identified by mass spectrometry in HEK293T or U2OS cells; '
        f'{int(ml.loc[bh_rows, "berryhill_site_probability"].ne("").sum()):,} with a site probability) on {ml.loc[bh_rows, "canonical_UniProtAC"].nunique():,} proteins. Of the {b["input_rows"]:,} rows of the table, '
        f'{b["relocated_sites"]:,} sites were moved to the lysine of their peptide on the canonical sequence and {sum(b["rejected_by_reason"].values()):,} rows were rejected '
        f'(' + ', '.join(f'{v:,} {k.replace("_", " ")}' for k, v in sorted(b['rejected_by_reason'].items())) + '). The other sheets of the supplements, which carry third-party information, are not used.', '',
        '**Ochoa 2020** (Ochoa D, Jarnuczak AF, Viéitez C, et al. "The functional landscape of the human phosphoproteome", Nat Biotechnol 38:365-373, 2020, doi:10.1038/s41587-019-0344-3). '
        'The article is **not open access** ("© 2019 The Author(s), under exclusive licence to Springer Nature America, Inc."), and its supplements carry no license of their own. '
        'The derived table (accession, position, residue, functional score, best localization probability, best PEP, biological samples, spectral counts) is published here by the owner\'s decision of 2026-11-10, cited. '
        'The CC BY-NC-SA 4.0 of this release does not grant rights in Ochoa et al.\'s data beyond what Springer Nature and the authors allow; a reader who needs the phosphorylation table without it removes that source from `sources`, '
        f'and the rows that only it reports. It contributes {int(oc_rows.sum()):,} phosphorylation rows (state `phospho`; `evidence_class` experimental, mass spectrometry) on {ph.loc[oc_rows, "canonical_UniProtAC"].nunique():,} proteins; '
        f'of its {o["input_sites"]:,} sites {sum(o["rejected_by_reason"].values()):,} were rejected on the reference check (' + ', '.join(f'{v:,} {k.replace("_", " ")}' for k, v in sorted(o['rejected_by_reason'].items())) +
        '; no window is given, so none was relocated).', '',
        f'**iPTMnet score** (`iptmnet_score`). iPTMnet\'s download has a `score.txt` (iPTMnet {sc["iptmnet_release"]}, {sc["url"]}); its readme says only "The score is calculated based on all the evidence sources as shown on the website". '
        f'It is an integer from 0 to 4 and is iPTMnet\'s own: this release carries it unchanged and does not explain it. `score.txt` scores the site (the rows without an enzyme) under the plain accession and, less often, the accession with the suffix -1 (its isoform 1); '
        'the plain accession\'s score is used. A row has the score when its `sources` include iPTMnet and `score.txt` has a row for the residue: ' + '; '.join(lines) + '.',
    ]


def write_docs(out: Path, version, release_date, tables, cross, summaries, commit) -> dict:
    files = {m: f'human_{m}_sites_{version}.tsv' for m in FILE_MARKS}
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
         'tables of context marks: the lysine methylation sites (UniProt, iPTMnet, CPLM 4.0 and Berryhill 2023), phosphorylation sites (UniProt, iPTMnet and Ochoa 2020) lysine acylation sites (UniProt and iPTMnet; acetylation also as a table of its own), lysine ubiquitination sites and lysine SUMOylation sites (each UniProt, iPTMnet and CPLM 4.0) of the human reviewed proteome, and the '
         'arginine methylation records of UniProt and iPTMnet as a cross-check. Sites for each mark read their own table and show the '
         'others as context. Positions of the new tables are on UniProt release 2026_03, the reference of the methylarginine data.', '',
         '## Context marks', '', '| File | One row per |', '|---|---|']
    for m in FILE_MARKS:
        t = tables[m]
        L.append(f'| `{files[m]}` | {m.replace("_", " ")} residue and {"acyl type" if m in ACYL else "state"} ({len(t):,} rows on {t[["canonical_UniProtAC", "position"]].drop_duplicates().shape[0]:,} '
                 f'residues of {t["canonical_UniProtAC"].nunique():,} proteins) |')
    L += [f'| `{files["methylarginine"]}` | arginine methylation record of UniProt or iPTMnet ({len(cross):,} rows), with whether the methylarginine data has the arginine. A cross-check, not a second methylarginine table |', '',
          'Columns are in `data_dictionary.tsv`. Each row says which sources report it (`sources`), the best `evidence_class`, UniProt\'s ECO codes, '
          'iPTMnet\'s upstream databases (`iptmnet_sources`), the PubMed IDs and, in the methyllysine, ubiquitination and sumoylation tables, CPLM\'s own evidence label (`cplm_evidence`: Exp., Exp./Dat., Dat.), in its own columns the figures of Berryhill 2023 (`berryhill_*`) and of Ochoa 2020 (`functional_score`, `best_localization_probability`, `best_pep`, `biological_samples`, `spectral_counts`) and, in the methyllysine, phosphorylation, acylation, ubiquitination and sumoylation tables, iPTMnet\'s own score (`iptmnet_score`); each of these is empty on the rows of other sources. A site can hide `by-similarity` (UniProt ECO:0000250 and other '
          'non-experimental records) and `unspecified` rows by filtering on `evidence_class`.', '',
          '### Sources and licenses', '', '| Source | Release | License | Terms |', '|---|---|---|---|']
    for s in SOURCE_ORDER:
        L.append(f'| {s} | {", ".join(releases[s])} | {SOURCES[s]["license"]} | {SOURCES[s]["url"]} |')
    L += ['', 'The release is CC BY-NC-SA 4.0 ([LICENSE.md](../../LICENSE.md)). UniProt (CC BY 4.0) is compatible with it. iPTMnet\'s license page says '
          'CC BY-NC-SA 4.0; its download readme also mentions CC BY 4.0 in its last line. The release follows the stricter terms: share-alike and '
          'non-commercial. iPTMnet asks that the upstream databases of its `source` column are cited when its data are redistributed: each row keeps them '
          '(`iptmnet_sources`), and they are HPRD, UniProt, RLIMS-P text mining, SIGNOR, neXtProt, IEDB, Protein Ontology and IntAct (terms are theirs). '
          'Not used: PhosphoSitePlus, dbPTM and qPTM. ', '',
          '**CPLM 4.0 states no license** (cplm.biocuckoo.cn carries a copyright notice, "All Rights Reserved"). Its human lysine methylation sites are included by the '
          'owner\'s decision of 2026-11-09, and its human lysine ubiquitination and SUMOylation sites by the owner\'s decision of 2026-11-10 on the same terms, cited (Zhang et al., Nucleic Acids Res 50:D451-D459, 2022, doi:10.1093/nar/gkab849), **pending the authors\' permission**; '
          'the terms the authors give, when they answer, are followed. The CC BY-NC-SA 4.0 of this release does not grant rights in CPLM\'s records beyond what the authors allow. '
          'CPLM contributes only to the methyllysine, ubiquitination and sumoylation tables, as the rows whose `sources` name CPLM (state `me`, `ub` or `sumo`, with `cplm_evidence`); a reader who needs a table without CPLM removes that source from `sources`, '
          'and the rows that only CPLM reports.', '',
          *source_notes(tables, summaries), '',
          '### Counts', '',
          f'Rows by state, evidence class and the sources that report them (iPTMnet {summaries["iptmnet"]["iptmnet_release"]} and CPLM give no methylation degree, so their rows are state `me` and do not meet UniProt\'s `me1`, `me2` and `me3` rows; ubiquitination rows are state `ub` and sumoylation rows state `sumo`, UniProt\'s SUMO paralog not kept).', '']
    for m in MARKS:
        ov = overlap(tables[m])
        L += [f'**{m}** ({len(tables[m]):,} rows). Residues by source, any state: ' + ', '.join(f'{k} {v:,}' for k, v in sorted(ov['sites'].items())) + '.', ''] + counts_table(tables[m]) + ['']
    ov = overlap(cross)
    L += ['**methylarginine cross-check** (' + f'{len(cross):,} rows). Residues by source, any state: ' + ', '.join(f'{k} {v:,}' for k, v in sorted(ov['sites'].items())) + '.', ''] + counts_table(cross) + ['']
    arg = cross.drop_duplicates(['canonical_UniProtAC', 'position'])['in_methylarginine_release'].value_counts().sort_index()
    L += ['Arginines of the cross-check in the methylarginine data: ' + ', '.join(f'{k} {v:,}' for k, v in arg.items()) + '.', '']
    unip_only = {m: int(tables[m]['iptmnet_sources'].eq('unip').sum()) for m in MARKS}
    unip_only_text = ', '.join(f'{unip_only[m]:,} {m.replace("_", " ")}' for m in MARKS)
    cp = tables['methyllysine']
    cp = cp[cp['sources'].str.contains('CPLM')]
    cu = tables['ubiquitination']
    cu = cu[cu['sources'].str.contains('CPLM')]
    cs = tables['sumoylation']
    cs = cs[cs['sources'].str.contains('CPLM')]
    L += [f'iPTMnet carries a copy of UniProt: {unip_only_text} rows have `unip` as their only upstream database, '
          'so "UniProt + iPTMnet" counts them as agreement of two sources when they are one.', '',
          f'CPLM 4.0 evidence of the {len(cp):,} methyllysine rows it reports (`cplm_evidence`): ' + ', '.join(f'{k} {v:,}' for k, v in cp['cplm_evidence'].value_counts().sort_index().items())
          + f'; of the {len(cu):,} ubiquitination rows it reports: ' + ', '.join(f'{k} {v:,}' for k, v in cu['cplm_evidence'].value_counts().sort_index().items())
          + f'; of the {len(cs):,} sumoylation rows it reports: ' + ', '.join(f'{k} {v:,}' for k, v in cs['cplm_evidence'].value_counts().sort_index().items())
          + '. CPLM compiles records from other databases too (`Dat.`), so a residue that CPLM and UniProt or iPTMnet both report may rest on the same original record; the tables only record that the residue is shared.', '',
          '### Dropped in staging', '',
          'Each count is of source records. UniProt: ' + ', '.join(f'{v:,} {k.replace("_", " ")}' for k, v in sorted(summaries['uniprot']['dropped'].items())) +
          f'; the entries\' sequences were compared with the reference (UniProt {summaries["uniprot"]["uniprot_release_served"]} against 2026_03). iPTMnet: ' +
          ', '.join(f'{v:,} {k.replace("_", " ")}' for k, v in sorted(summaries['iptmnet']['dropped'].items())) +
          f'. CPLM: {summaries["cplm"]["records_in_download"]:,} rows in the methylation download, {summaries["cplm"]["records_in_scope"]:,} human methylation rows, of which ' +
          ', '.join(f'{v:,} {k.replace("_", " ")}' for k, v in sorted(summaries['cplm']['dropped'].items())) + f' (nothing remapped); {summaries["cplm_ub"]["records_in_download"]:,} rows in the ubiquitination download, {summaries["cplm_ub"]["records_in_scope"]:,} human ubiquitination rows, of which ' +
          ', '.join(f'{v:,} {k.replace("_", " ")}' for k, v in sorted(summaries['cplm_ub']['dropped'].items())) + f' (nothing remapped); {summaries["cplm_su"]["records_in_download"]:,} rows in the sumoylation download, {summaries["cplm_su"]["records_in_scope"]:,} human sumoylation rows, of which ' +
          ', '.join(f'{v:,} {k.replace("_", " ")}' for k, v in sorted(summaries['cplm_su']['dropped'].items())) + ' (nothing remapped). '
          f'Berryhill 2023: {summaries["berryhill"]["input_rows"]:,} rows of Supplementary Table 5, ' + ', '.join(f'{v:,} {k.replace("_", " ")}' for k, v in sorted(summaries['berryhill']['rejected_by_reason'].items())) +
          f' (and {summaries["berryhill"]["relocated_sites"]:,} sites relocated by their peptide). Ochoa 2020: {summaries["ochoa"]["input_sites"]:,} sites, ' + ', '.join(f'{v:,} {k.replace("_", " ")}' for k, v in sorted(summaries['ochoa']['rejected_by_reason'].items())) +
          f' (nothing relocated). iPTMnet scores: {summaries["score"]["records_in_download"]:,} rows of score.txt, ' + ', '.join(f'{v:,} {k.replace("_", " ")}' for k, v in sorted(summaries['score']['dropped'].items())) + '.', '',
          '### How the context marks were made', '',
          f'`python scripts/54_stage_uniprot_ptm_sites.py`, `python scripts/55_stage_iptmnet_ptm_sites.py` and `python scripts/39_stage_cplm4_lys_methylation.py --marks-dir inputs/marks` (also with `--modification Ubiquitination` and `--modification Sumoylation`), `python scripts/57_stage_berryhill_lys_methylation.py`, `python scripts/58_stage_ochoa_phosphosites.py` and `python scripts/59_stage_iptmnet_scores.py` (inputs in `inputs/marks/`, downloads listed in `inputs/marks/raw_files.tsv`), then '
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
    score_paths = sorted(inputs.glob('iptmnet_*_score.tsv'))
    if len(score_paths) != 1:
        raise SystemExit(f'expected one staged iPTMnet score table (scripts/59) in {inputs}, found {len(score_paths)}')
    scores = pd.read_csv(score_paths[0], sep='\t', dtype=str, keep_default_na=False)
    for m in MARKS:
        tables[m] = add_iptmnet_scores(tables[m], scores)
    used['iptmnet_score'] = {score_paths[0].name: md5_file(score_paths[0])}
    cross = crosscheck(tables.pop('methylarginine'), out, version)

    summaries = {'uniprot': json.loads((inputs / 'uniprot_ptm_sites_summary.json').read_text(encoding='utf-8')),
                 'iptmnet': json.loads((inputs / 'iptmnet_ptm_sites_summary.json').read_text(encoding='utf-8')),
                 'cplm': json.loads((inputs / 'cplm_ptm_sites_summary.json').read_text(encoding='utf-8')),
                 'cplm_ub': json.loads((inputs / 'cplm_ubiquitination_summary.json').read_text(encoding='utf-8')),
                 'cplm_su': json.loads((inputs / 'cplm_sumoylation_summary.json').read_text(encoding='utf-8')),
                 'berryhill': json.loads((inputs / 'berryhill2023_methyllysine_summary.json').read_text(encoding='utf-8')),
                 'ochoa': json.loads((inputs / 'ochoa2020_phosphorylation_summary.json').read_text(encoding='utf-8')),
                 'score': json.loads((inputs / 'iptmnet_score_summary.json').read_text(encoding='utf-8'))}
    tables['acetylation'] = tables['lysine_acylation'][tables['lysine_acylation']['state'] == 'acetyl'].reset_index(drop=True)
    for m in FILE_MARKS:
        t = tables[m].rename(columns={'state': 'acyl_type'}) if m in ACYL else tables[m]
        t.to_csv(out / f'human_{m}_sites_{version}.tsv', sep='\t', index=False, lineterminator='\n')
    cross.to_csv(out / f'human_methylarginine_crosscheck_{version}.tsv', sep='\t', index=False, lineterminator='\n')

    # data dictionary: v1.4.3's, then the new tables
    new = {**{m: f'human_{m}_sites_{version}.tsv' for m in FILE_MARKS}, 'methylarginine': f'human_methylarginine_crosscheck_{version}.tsv'}
    lines = []
    for line in (src / 'data_dictionary.tsv').read_text(encoding='utf-8').rstrip('\n').split('\n'):
        name, rest = line.split('\t', 1)
        name = name.replace(f'_v1.4.2_to_{CARRIED_FROM}', f'_{CARRIED_FROM}_to_{version}')
        lines.append(renamed(name, version) + '\t' + rest)
    headers = {name: list((cross if mark == 'methylarginine' else tables[mark].rename(columns={'state': 'acyl_type'} if mark in ACYL else {})).columns)
               for mark, name in new.items()}
    for mark, name in new.items():
        for column, kind, description in DICTIONARY:
            if column not in headers[name]:
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
                'The sources are UniProt and iPTMnet, which are openly licensed; CPLM 4.0 (methyllysine, ubiquitination and sumoylation), which states no license; Berryhill 2023 (methyllysine), whose article is CC BY 4.0 and whose supplements carry no notice; '
                'and Ochoa 2020 (phosphorylation), whose article is not open access (all below, with the iptmnet_score column).', '']
    for m in FILE_MARKS:
        c = counts[m]
        changes.append(f'- `human_{m}_sites_{version}.tsv`: {c["rows"]:,} rows on {c["residues"]:,} residues of {c["proteins"]:,} proteins; by {"acyl type" if m in ACYL else "state"} '
                       + ', '.join(f'{k} {v:,}' for k, v in c['by_state'].items()) + '; by evidence class ' + ', '.join(f'{k} {v:,}' for k, v in c['by_evidence_class'].items())
                       + '; residues by source ' + ', '.join(f'{k} {v:,}' for k, v in sorted(c['by_source']['sites'].items())) + '.')
    c = counts['methylarginine_crosscheck']
    changes.append(f'- `human_methylarginine_crosscheck_{version}.tsv`: {c["rows"]:,} arginine methylation rows of UniProt and iPTMnet with whether the methylarginine data has the arginine; residues by source '
                   + ', '.join(f'{k} {v:,}' for k, v in sorted(c['by_source']['sites'].items())) + '.')
    cc, cu, cs = summaries['cplm'], summaries['cplm_ub'], summaries['cplm_su']
    changes += ['', 'Sources: UniProt (release ' + ', '.join(sorted({summaries["uniprot"]["uniprot_release_served"]})) + ', CC BY 4.0), iPTMnet (release ' + summaries['iptmnet']['iptmnet_release'] +
                ', CC BY-NC-SA 4.0) and CPLM 4.0, Berryhill 2023 and Ochoa 2020 (below). The release stays CC BY-NC-SA 4.0. iPTMnet\'s download readme also mentions CC BY 4.0; the stricter terms are followed.', '',
                '### CPLM 4.0: license note', '',
                f'The methyllysine table includes all {cc["rows"]["methyllysine"]:,} validated human lysine methylation sites of CPLM 4.0 ({cc["proteins"]:,} proteins; Zhang et al., Nucleic Acids Res 50:D451-D459, 2022, '
                f'doi:10.1093/nar/gkab849; `Methylation.txt` of {cc["cplm_release"].split("(")[1].rstrip(")")}). **CPLM states no license**: cplm.biocuckoo.cn carries only a copyright notice ("All Rights Reserved"). '
                'The sites are included by the owner\'s decision of 2026-11-09 (confirmed 2026-11-10), cited, **pending the authors\' permission**; the terms the authors give, when they answer, are followed. The CC BY-NC-SA 4.0 of the release does not grant rights in CPLM\'s records beyond what the authors allow. '
                'Each CPLM row has `CPLM` in `sources` and CPLM\'s own evidence label in `cplm_evidence` (' + ', '.join(f'{k} {v:,}' for k, v in cc['rows_by_cplm_evidence'].items()) +
                '), so a reader can filter on it or leave CPLM out. CPLM rows are state `me` (CPLM gives no degree); `evidence_class` is experimental for Exp. and Exp./Dat., and unspecified for Dat. '
                f'(a record CPLM took from another database). `pmids` holds the references that are PubMed IDs ({cc["rows_with_pubmed_ids"]:,} rows; the rest of CPLM\'s reference field is identifiers of other databases). '
                f'Dropped on validation, nothing remapped: ' + ', '.join(f'{v:,} {k.replace("_", " ")}' for k, v in sorted(cc['dropped'].items())) + '. '
                f'The ubiquitination table likewise includes all {cu["rows"]["ubiquitination"]:,} validated human lysine ubiquitination sites of CPLM 4.0 ({cu["proteins"]:,} proteins; `Ubiquitination.txt` of {cu["cplm_release"].split("(")[1].rstrip(")")}), '
                'by the owner\'s decision of 2026-11-10 on the same terms as the methylation sites: cited (Zhang et al., Nucleic Acids Res 50:D451-D459, 2022, doi:10.1093/nar/gkab849), **pending the authors\' permission**, the authors\' terms followed when given. '
                'Its rows have `CPLM` in `sources`, state `ub`, and CPLM\'s evidence label in `cplm_evidence` (' + ', '.join(f'{k} {v:,}' for k, v in cu['rows_by_cplm_evidence'].items()) + '); '
                f'`pmids` holds the PubMed IDs among CPLM\'s references ({cu["rows_with_pubmed_ids"]:,} rows). Dropped on validation, nothing remapped: ' + ', '.join(f'{v:,} {k.replace("_", " ")}' for k, v in sorted(cu['dropped'].items())) + '. '
                f'The sumoylation table likewise includes all {cs["rows"]["sumoylation"]:,} validated human lysine SUMOylation sites of CPLM 4.0 ({cs["proteins"]:,} proteins; `Sumoylation.txt` of {cs["cplm_release"].split("(")[1].rstrip(")")}), '
                'by the owner\'s decision of 2026-11-10 on the same terms: cited (Zhang et al., Nucleic Acids Res 50:D451-D459, 2022, doi:10.1093/nar/gkab849), **pending the authors\' permission**, the authors\' terms followed when given. '
                'Its rows have `CPLM` in `sources`, state `sumo`, and CPLM\'s evidence label in `cplm_evidence` (' + ', '.join(f'{k} {v:,}' for k, v in cs['rows_by_cplm_evidence'].items()) + '); '
                f'`pmids` holds the PubMed IDs among CPLM\'s references ({cs["rows_with_pubmed_ids"]:,} rows). Dropped on validation, nothing remapped: ' + ', '.join(f'{v:,} {k.replace("_", " ")}' for k, v in sorted(cs['dropped'].items())) + '. '
                'CPLM contributes no lysine acylation.', '',
                '### Berryhill 2023, Ochoa 2020 and the iPTMnet score: license notes', '', *source_notes(tables, summaries), '',
                '### Ubiquitination and SUMOylation', '',
                'The two tables hold lysine sites (state `ub` and `sumo`) of UniProt (its cross-link features: the Lys-Gly isopeptide to ubiquitin, or to any SUMO paralog, the paralog not kept), iPTMnet (UBIQUITINATION and SUMOYLATION at K) and CPLM 4.0 (above). '
                'UniProt and iPTMnet are openly licensed (CC BY 4.0, CC BY-NC-SA 4.0). iPTMnet carries a copy of UniProt, so a residue reported by both may be one record.', '',
                '### Lysine acylation', '',
                'The two acylation tables hold the lysine acylation (13 acyl types of UniProt\'s N6-acyl lysine terms, and iPTMnet\'s acetylation) of UniProt and iPTMnet, with `acyl_type` in place of `state`; '
                '`human_acetylation_sites_' + version + '.tsv` repeats the acetyl rows because they are over half of them. Both sources are openly licensed (UniProt CC BY 4.0, iPTMnet CC BY-NC-SA 4.0).', '',
                'Code: `scripts/39_stage_cplm4_lys_methylation.py`, `scripts/54_stage_uniprot_ptm_sites.py`, `scripts/55_stage_iptmnet_ptm_sites.py`, `scripts/57_stage_berryhill_lys_methylation.py`, `scripts/58_stage_ochoa_phosphosites.py`, `scripts/59_stage_iptmnet_scores.py`, `scripts/56_build_context_marks.py`; staged tables in `inputs/marks/`.']
    (out / 'CHANGES.md').write_text('\n'.join(changes) + '\n', encoding='utf-8')

    release = {'title': 'Human methylarginine sites and context marks', 'version': version, 'previous_version': CARRIED_FROM,
               'release_date': args.release_date, 'repository': 'https://github.com/Shechterlab/PTM_Pipeline',
               'license': 'CC BY-NC-SA 4.0 (LICENSE.md); sources keep their own terms (inputs/MANIFEST.tsv); CPLM 4.0 states no license and its methylation, ubiquitination and SUMOylation sites are included pending the authors\' permission; Berryhill 2023 supplements carry no notice (article CC BY 4.0); Ochoa 2020 is not open access; both included by the owner\'s decision',
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

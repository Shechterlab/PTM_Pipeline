#!/usr/bin/env python3
"""Stage the PRMT7 knockout arginine-monomethylation response (Szewczyk et al. 2020) on the canonical human proteome.

Szewczyk, Ishikawa, Organ, ..., Nara & Barsyte-Lovejoy, "Pharmacological inhibition of PRMT7 links arginine
monomethylation to the cellular stress response", Nat Commun 11:2396 (2020), doi:10.1038/s41467-020-16271-z. The article
is CC BY 4.0 (its Rights and permissions section); the publisher's pages state no licence for the Supplementary Data
files themselves, so the staged tables are checksum_only until that is settled. Raw MS files: MassIVE MSV000084773
(doi:10.25345/C5GQ37), not used.

Supplementary Data 2 (41467_2020_16271_MOESM5_ESM.xlsx, sheet IP_All_Hits_Sites) lists the 346 arginine-monomethylated
peptides of a monomethylarginine (Rme1) immunoprecipitation of SILAC-labelled HCT116 cells, PRMT7 knockout (heavy) over
parental (light): each named by its leading UniProt accession(s), gene and Rme1 site, GENE (Rn), with the limma log2 H/L
over four replicates, its P value and Benjamini-Hochberg adjusted P value, the replicate H/L ratios, the paper's call
(Significantly Depleted) and whether PhosphoSitePlus had the site. It gives no peptide sequence. The paper tested the
inhibitor SGC3027 (and SGC3027N) by western blots of HSP70 only, so the knockout is the one comparison. Sheet1 of the
same file (Supplementary Table 4) and Source Data Fig. 2c repeat rows of this sheet. Supplementary Data 3 (MOESM6, the
input proteome) gives the FASTA header, and so the sequence version, of each entry in the database the paper searched
(UniProt 2018_08).

Each row is placed at its reported position on the canonical sequence of its leading (first) accession, and kept when
that residue is an arginine (common.residue_is): position_verified when the entry has the sequence version that was
searched. An entry the input proteome does not hold is checked against the sequence UniProt had in release 2018_08
(UniSave, --history, fetched once with --fetch-history): position_verified_unisave when that sequence is the reference
sequence, lifted when it differs and the 15 residues around the arginine occur once in the reference sequence (the
position then moves with them), and dropped otherwise. The other accessions of a shared peptide are listed, not placed
(the position is numbered on the leading entry). Rows of one entry whose replicate H/L ratios are
all identical are one quantified peptide carrying those methylarginines together.

A peptide is decreased at log2 KO/WT <= -1 with adjusted P < 0.01: the paper's rule (Methods, "Significant hits were
called as H/L ratio of <-1 ... and a Benjamini-Hochberg adjusted p-value of <0.01 (n = 4)"; it gives its Significantly
Depleted column exactly, and no row sits at log2 = -1). It is increased at log2 >= 1 with adjusted P < 0.01: the paper
calls depletion only, and the increase is the same rule turned round, as scripts/41 does for CARM1. A site takes the
call of its peptides by the shared rule (scripts/peptide_calls.py), and the call of its single-site peptides alone.

Outputs (in --outdir)
  prmt7_ko_peptides.tsv     every row of Supplementary Data 2 with its placement, values and call
  prmt7_ko_sites.tsv        one row per canonical site and comparison
  prmt7_ko_summary.json     input MD5s, rule, counts, the rows left out and why
  manifest_rows.tsv         the staged tables' rows for inputs/MANIFEST.tsv
  raw_files_rows.tsv        the downloads' rows for inputs/layers/raw_files.tsv
With --manifest, the rows also go into that MANIFEST.tsv and the raw_files.tsv in --outdir (as scripts/44 does).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook

import peptide_calls as pc
from common import md5_file, parse_fasta, refresh_manifest_row, residue_is, save_table, update_manifest, upsert_tsv_rows

STUDY = 'Szewczyk 2020'
DOI = '10.1038/s41467-020-16271-z'
CITATION = 'Szewczyk et al., Nat Commun 11:2396 (2020), doi:10.1038/s41467-020-16271-z'
COMPARISON, CELL_LINE = 'PRMT7_KO_vs_WT', 'HCT116'
LOG2_CUTOFF, Q_CUTOFF = 1.0, 0.01
RULE = (f'peptide log2 KO/WT <= -{LOG2_CUTOFF:g} (decreased) or >= {LOG2_CUTOFF:g} (increased) with BH-adjusted p < {Q_CUTOFF:g}; '
        'sites by scripts/peptide_calls.py')
SOURCE_RULE = ('Significant hits were called as H/L ratio of <-1 (knockout cells (H) relative to control (L)) and a Benjamini-Hochberg '
               'adjusted p-value of <0.01 (n = 4) (Methods; Fig. 2c legend; Supplementary Table 4)')
SITE = re.compile(r'^(.*) \(R(\d+)\)$')
SD2_SHEET = 'IP_All_Hits_Sites'
SD2_COLUMNS = {'protein': 'protein', 'gene': 'gene names', 'log2': 'logfc', 'p': 'p.value', 'q': 'adj.p.val', 'site': 'site',
               'depleted': 'significantly depleted', 'psp': 'previously_ observed_ phosphositedb'}
RATIOS = ['ratio_hl_r1', 'ratio_hl_r2', 'ratio_hl_r3', 'ratio_hl_r4']
KEPT = ('position_verified', 'position_verified_unisave', 'lifted')
SEARCHED_RELEASE = '2018_08'   # the UniProt release the paper searched (Supplementary Data 3)
UNISAVE = 'https://rest.uniprot.org/unisave/'

# What inputs/MANIFEST.tsv and inputs/layers/raw_files.tsv say of the staged tables and their downloads.
ESM = 'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-020-16271-z/MediaObjects/'
SD2_NAME, SD3_NAME = '41467_2020_16271_MOESM5_ESM.xlsx', '41467_2020_16271_MOESM6_ESM.xlsx'
LICENCE = ('CC BY 4.0: the article\'s licence (nature.com Rights and permissions; Crossref), taken to cover its Supplementary Data, '
           'which state none of their own (the lab\'s decision of 2026-10-25, as for the CARM1 knockout of Shishkova et al. 2017)')
REDISTRIBUTION = 'git'
HOW = ('scripts/45_stage_prmt7_szewczyk2020.py; checked against reference/uniprot_human_reviewed_2026_03.fasta: the reported position on '
       'the sequence version searched (the 2018_08 FASTA headers of Supplementary Data 3, or UniSave for entries they lack, '
       'layers/prmt7_ko_searched_sequences.tsv), kept where it is an arginine there and maps to the reference sequence')
MANIFEST_ROWS = {
    'prmt7_ko_peptides.tsv': ('Rme1 peptides quantified by monomethylarginine immunoprecipitation of SILAC-labelled PRMT7 knockout and parental '
        'HCT116 cells (Supplementary Data 2, one row each, named by Rme1 site), placed at the reported position on the canonical sequence, with '
        f'log2 KO/WT, p, BH-adjusted p and call (decreased: log2 <= -1, increased: >= 1, both at adjusted p < {Q_CUTOFF:g})'),
    'prmt7_ko_sites.tsv': ('One row per canonical methylarginine and comparison (PRMT7 KO vs WT, HCT116): its peptides, median log2 KO/WT, lowest '
        'p and adjusted p, and the call by scripts/peptide_calls.py over all its peptides and over its single-site peptides'),
    'prmt7_ko_searched_sequences.tsv': ('The sequence UniProt release 2018_08 had for each entry of Supplementary Data 2 that the input '
        'proteome of Supplementary Data 3 does not list (UniSave: entry version, sequence version, sequence), to check and move '
        'their positions to the reference sequence'),
    'prmt7_ko_summary.json': ('Counts, rule and input MD5s of the PRMT7 knockout staging, and the rows left out and why (file names only; '
        'counts and checksums, no source values)'),
}
RAW_FILES = {SD2_NAME: 'Supplementary Data 2 (Identification of monomethyl arginine peptides in PRMT7 WT and KO cells); publisher file '
                       'Last-Modified 2023-12-04',
             SD3_NAME: 'Supplementary Data 3 (Input peptide level analysis in PRMT7 WT and KO cells); publisher file Last-Modified 2023-12-04'}
RAW_KEPT = 'release_asset (CC BY 4.0, the article\'s licence)'


def number(v) -> float | None:
    """A cell as a number; 'NaN', 'NA' and empty cells are None."""
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(x) else x


def read_sd2(path: Path) -> pd.DataFrame:
    """One row per peptide of Supplementary Data 2, from the columns named in the sheet's header row."""
    wb = load_workbook(path, read_only=True, data_only=True)
    if SD2_SHEET not in wb.sheetnames:
        raise SystemExit(f'{path}: no sheet {SD2_SHEET!r} in {wb.sheetnames}')
    rows = wb[SD2_SHEET].iter_rows(values_only=True)
    header = [str(h).strip().lower() if h is not None else '' for h in next(rows)]
    missing = [c for c in [*SD2_COLUMNS.values(), *RATIOS] if c not in header]
    if missing:
        raise SystemExit(f'{path} sheet {SD2_SHEET}: no columns {missing} in {header}')
    at = {k: header.index(c) for k, c in SD2_COLUMNS.items()}
    out = []
    for k, row in enumerate(rows, 2):
        if not row or row[at['site']] is None:
            continue
        m = SITE.match(str(row[at['site']]).strip())
        if not m:
            raise SystemExit(f'{path} row {k}: site {row[at["site"]]!r} is not GENE (Rn)')
        accs = [a.strip() for a in str(row[at['protein']]).split('/') if a.strip()]
        ratios = tuple(number(row[header.index(c)]) for c in RATIOS)
        out.append({'source_row': k, 'reported_site': str(row[at['site']]).strip(), 'source_UniProtAC': str(row[at['protein']]).strip(),
                    'leading': accs[0], 'other_accessions': ';'.join(accs[1:]), 'gene': str(row[at['gene']] or '').strip(),
                    'reported_position': int(m.group(2)), 'log2_ko_wt': number(row[at['log2']]), 'p': number(row[at['p']]),
                    'q': number(row[at['q']]), 'n_replicates': sum(r is not None for r in ratios), 'ratios': ratios,
                    'significantly_depleted': str(row[at['depleted']] or '').strip(),
                    'phosphositeplus_reported': str(row[at['psp']] or '').strip()})
    return pd.DataFrame(out)


def searched_versions(path: Path) -> dict[str, int]:
    """Sequence version (SV=) of each entry in the FASTA headers of Supplementary Data 3, the database the paper searched."""
    wb = load_workbook(path, read_only=True, data_only=True)
    ws = wb.worksheets[1]
    rows = ws.iter_rows(values_only=True)
    header = [str(h).strip() if h is not None else '' for h in next(rows)]
    if 'Fasta headers' not in header:
        raise SystemExit(f'{path} sheet {ws.title}: no column Fasta headers')
    col = header.index('Fasta headers')
    versions, conflicts = {}, set()
    for row in rows:
        for part in str(row[col] or '').split(';'):
            m = re.match(r'\s*>?(?:sp|tr)\|([^|]+)\|.*\bSV=(\d+)', part)
            if m:
                acc, sv = m.group(1), int(m.group(2))
                if versions.setdefault(acc, sv) != sv:
                    conflicts.add(acc)
    if conflicts:
        raise SystemExit(f'{path}: entries with two sequence versions: {sorted(conflicts)}')
    return versions


def reference_versions(path: Path) -> dict[str, int]:
    """Sequence version (SV=) of each entry of the reference FASTA."""
    out = {}
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('>'):
                m = re.search(r'\bSV=(\d+)', line)
                out[line[1:].split('|')[1]] = int(m.group(1)) if m else None
    return out


def fetch_history(accessions: list[str], release: str = SEARCHED_RELEASE) -> pd.DataFrame:
    """For each entry, the UniSave entry version current in UniProt release `release`, and its sequence. Releases before
    2010 are numbered otherwise ('15.0/55.0') and come before any 'YYYY_MM' release."""
    from urllib.request import urlopen
    dated = lambda r: re.fullmatch(r'\d{4}_\d{2}', r) is not None  # noqa: E731
    rows = []
    for acc in sorted(set(accessions)):
        with urlopen(f'{UNISAVE}{acc}?format=json', timeout=120) as fh:
            versions = json.load(fh)['results']
        before = [v for v in versions if not dated(v['firstRelease'].split('/')[0]) or v['firstRelease'].split('/')[0] <= release]
        if not before:
            rows.append({'accession': acc, 'release': release, 'entry_version': '', 'sequence_version': '', 'sequence': ''})
            continue
        v = max(before, key=lambda x: x['entryVersion'])
        last = v['lastRelease'].split('/')[0]
        if dated(last) and last < release:   # the entry was gone by then
            rows.append({'accession': acc, 'release': release, 'entry_version': '', 'sequence_version': '', 'sequence': ''})
            continue
        with urlopen(f'{UNISAVE}{acc}?format=fasta&versions={v["entryVersion"]}', timeout=120) as fh:
            seq = ''.join(line.strip() for line in fh.read().decode().splitlines() if line and not line.startswith('>'))
        rows.append({'accession': acc, 'release': release, 'entry_version': v['entryVersion'], 'sequence_version': v['sequenceVersion'], 'sequence': seq})
    return pd.DataFrame(rows)


def place(row, canonical: dict[str, str], searched: dict[str, int], reference: dict[str, int], history: dict[str, str]) -> tuple[str, int | None]:
    """Placement of one row at its reported position: (status, position on the leading entry's canonical sequence)."""
    acc, pos = row.leading, row.reported_position
    if acc not in canonical:
        return 'accession_not_in_canonical_fasta', None
    if acc in searched and searched[acc] == reference.get(acc):
        return ('position_verified', pos) if residue_is(canonical, acc, pos) else ('residue_not_arginine', None)
    old = history.get(acc, '')
    if not old:
        return 'no_searched_sequence', None
    if not (0 < pos <= len(old) and old[pos - 1] == 'R'):
        return 'not_arginine_in_searched_sequence', None
    if old == canonical[acc]:
        return 'position_verified_unisave', pos
    window = old[max(0, pos - 8):pos + 7]
    hits = [m.start() for m in re.finditer(f'(?={re.escape(window)})', canonical[acc])]
    if len(hits) != 1:
        return 'not_lifted', None
    return 'lifted', hits[0] + (pos - 1 - max(0, pos - 8)) + 1


def call(log2: float | None, q: float | None) -> str:
    """One peptide by the shared rule: changed at adjusted p < Q_CUTOFF with |log2| >= LOG2_CUTOFF."""
    if log2 is None or q is None:
        return 'not_quantified'
    return pc.site_response([(q, log2, 1)], Q_CUTOFF, LOG2_CUTOFF, strict=True)['call']


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--sd2', required=True, help=f'Szewczyk 2020 Supplementary Data 2 ({SD2_NAME})')
    ap.add_argument('--sd3', required=True, help=f'Szewczyk 2020 Supplementary Data 3 ({SD3_NAME}), for the sequence versions searched')
    ap.add_argument('--fasta', default='inputs/reference/uniprot_human_reviewed_2026_03.fasta')
    ap.add_argument('--history', default='', help='the 2018_08 sequences of the entries Supplementary Data 3 lacks (prmt7_ko_searched_sequences.tsv)')
    ap.add_argument('--fetch-history', action='store_true', help='fetch them from UniSave into --history first (network)')
    ap.add_argument('--downloaded', default='2026-10-25', help='the date the two files were downloaded (raw_files.tsv)')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--manifest', default='', help='inputs/MANIFEST.tsv: also write the rows there (and into raw_files.tsv in --outdir); '
                    'empty to leave both alone')
    args = ap.parse_args()
    sd2, sd3, fasta, out = Path(args.sd2), Path(args.sd3), Path(args.fasta), Path(args.outdir)

    peptides = read_sd2(sd2)
    canonical = parse_fasta(fasta)
    searched, reference = searched_versions(sd3), reference_versions(fasta)
    need = sorted({a for a in peptides['leading'] if a in canonical and not (a in searched and searched[a] == reference.get(a))})
    if args.fetch_history:
        if not args.history:
            raise SystemExit('--fetch-history needs --history, the file to write')
        save_table(fetch_history(need), Path(args.history))
    hist = pd.read_csv(args.history, sep='\t', dtype=str, keep_default_na=False) if args.history else pd.DataFrame(columns=['accession', 'sequence'])
    history = dict(zip(hist['accession'], hist['sequence']))
    missing = [a for a in need if a not in history]
    if args.history and missing:
        raise SystemExit(f'{args.history} lacks {missing}: fetch it again (--fetch-history)')
    where = [place(r, canonical, searched, reference, history) for r in peptides.itertuples(index=False)]
    peptides['placement'] = [w[0] for w in where]
    peptides['canonical_position'] = [w[1] for w in where]
    placed = peptides['placement'].isin(KEPT)
    peptides['canonical_UniProtAC'] = peptides['leading'].where(placed, '')
    peptides['searched_sequence_version'] = [str(searched[a]) if a in searched else '' for a in peptides['leading']]
    peptides['reference_sequence_version'] = [str(reference[a]) if reference.get(a) is not None else '' for a in peptides['leading']]
    peptides['methyl_positions'] = [str(int(p)) if k else '' for p, k in zip(peptides['canonical_position'], placed)]
    peptides['call'] = [call(f, q) for f, q in zip(peptides['log2_ko_wt'], peptides['q'])]

    # One quantified peptide can carry several of the sites the table names: those rows of the entry have the same
    # replicate ratios (and so the same log2, p and q).
    same = defaultdict(list)
    for r in peptides.itertuples(index=False):
        if r.n_replicates:
            same[(r.source_UniProtAC, r.ratios)].append(r.reported_position)
    group = [same.get((a, k), [p]) if n else [p] for a, k, p, n in
             zip(peptides['source_UniProtAC'], peptides['ratios'], peptides['reported_position'], peptides['n_replicates'])]
    peptides['n_methyl_sites'] = [len(g) for g in group]
    peptides['shares_quantification_with'] = [';'.join(f'R{x}' for x in sorted(g) if x != p) for g, p in zip(group, peptides['reported_position'])]
    multi = peptides[peptides['n_methyl_sites'] > 1]
    for _, g in multi.groupby([multi['source_UniProtAC'], multi['ratios'].map(str)]):
        if g[['log2_ko_wt', 'p', 'q']].nunique().max() > 1:
            raise SystemExit(f'rows with the same replicate ratios but different log2, p or q:\n{g}')

    # The rule against the paper's own column.
    depleted = peptides['significantly_depleted'].eq('Yes')
    disagree = peptides[peptides['call'].eq('decreased') != depleted]
    at_boundary = int((peptides['log2_ko_wt'].abs() == LOG2_CUTOFF).sum())

    # Sites: each placed, quantified row gives its site one peptide.
    groups = defaultdict(list)
    for r in peptides[placed & peptides['call'].ne('not_quantified')].itertuples(index=False):
        groups[(r.canonical_UniProtAC, int(r.canonical_position))].append(r)
    sites = []
    for (acc, pos), rows in sorted(groups.items()):
        res = pc.site_response([(r.q, r.log2_ko_wt, r.n_methyl_sites) for r in rows], Q_CUTOFF, LOG2_CUTOFF, strict=True)
        sites.append({'canonical_UniProtAC': acc, 'position': pos, 'residue': canonical[acc][pos - 1], 'comparison': COMPARISON,
                      'cell_line': CELL_LINE, 'n_peptides': res['n_peptides'], 'log2_ko_wt_median': res['log2_median'],
                      'log2_ko_wt_strongest': res['log2_strongest'], 'p_min': min(r.p for r in rows), 'q_min': res['adj_p_min'],
                      'n_decreased': res['n_decreased'], 'n_increased': res['n_increased'],
                      'n_multi_site_peptides': res['n_multi_site_peptides'], 'call': res['call'],
                      'call_single_site_peptides': res['call_single_site_peptides'],
                      'significantly_depleted_in_source': 'Yes' if any(r.significantly_depleted == 'Yes' for r in rows) else 'No',
                      'placement': ';'.join(sorted({r.placement for r in rows}))})
    sites = pd.DataFrame(sites)
    if sites['residue'].ne('R').any():
        raise SystemExit(f'placed sites that are not arginines:\n{sites[sites["residue"].ne("R")]}')

    columns = ['comparison', 'cell_line', 'source_row', 'reported_site', 'source_UniProtAC', 'gene', 'reported_position', 'log2_ko_wt', 'p',
               'q', 'n_replicates', 'significantly_depleted', 'phosphositeplus_reported', 'placement', 'canonical_UniProtAC',
               'other_accessions', 'searched_sequence_version', 'reference_sequence_version', 'methyl_positions', 'n_methyl_sites',
               'shares_quantification_with', 'call']
    peptides = peptides.assign(comparison=COMPARISON, cell_line=CELL_LINE)[columns]
    save_table(peptides, out / 'prmt7_ko_peptides.tsv')
    save_table(sites, out / 'prmt7_ko_sites.tsv')
    if args.history:
        # the searched sequences are staged beside the tables, so the staging can be run again offline
        hp = Path(args.history)
        if hp.resolve() != (out / 'prmt7_ko_searched_sequences.tsv').resolve():
            save_table(pd.read_csv(hp, sep='\t', dtype=str, keep_default_na=False), out / 'prmt7_ko_searched_sequences.tsv')
    q_bh = pc.bh(peptides['p'].tolist())
    summary = {
        'study': STUDY, 'doi': DOI, 'article_licence': 'CC BY 4.0', 'raw_data': 'MassIVE MSV000084773 (doi:10.25345/C5GQ37), not used',
        'comparison': COMPARISON, 'cell_line': CELL_LINE, 'rule': RULE, 'source_rule': SOURCE_RULE,
        'inputs': {sd2.name: {'md5': md5_file(sd2), 'sheet': SD2_SHEET}, sd3.name: {'md5': md5_file(sd3), 'use': 'FASTA header sequence versions'},
                   fasta.name: {'md5': md5_file(fasta)}},
        'rows': len(peptides), 'entries': int(peptides['source_UniProtAC'].nunique()),
        'placement': dict(sorted(Counter(peptides['placement']).items())),
        'dropped': {why: g['reported_site'].tolist() for why, g in peptides[~peptides['placement'].isin(KEPT)].groupby('placement')},
        'shared_peptides': {r.reported_site: {'placed_on': r.canonical_UniProtAC, 'other_accessions': r.other_accessions}
                            for r in peptides[peptides['other_accessions'] != ''].itertuples(index=False)},
        'multi_site_peptides': {'peptides': int(multi.groupby([multi['source_UniProtAC'], multi['ratios'].map(str)]).ngroups),
                                'rows': len(multi)},
        'peptide_calls': dict(sorted(Counter(peptides['call']).items())),
        'source_significantly_depleted': int(depleted.sum()),
        'rule_disagrees_with_source': disagree['reported_site'].tolist(),
        'rows_at_log2_cutoff': at_boundary,
        'adj_p_minus_bh_of_p_max_abs': round(max(abs(a - b) for a, b in zip(peptides['q'], q_bh)), 6),
        'sites': dict(sorted(Counter(sites['call']).items())),
        'sites_single_site_call': dict(sorted(Counter(sites['call_single_site_peptides'].replace('', 'none')).items())),
        'unique_sites': len(sites),
        'downloaded': args.downloaded,
    }
    (out / 'prmt7_ko_summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')

    # Their rows for inputs/MANIFEST.tsv and inputs/layers/raw_files.tsv.
    mrows = []
    for name, content in MANIFEST_ROWS.items():
        path = out / name
        n = '' if name.endswith('.json') else len(path.read_text(encoding='utf-8').splitlines()) - 1
        mrows.append({'file': f'layers/{name}', 'bytes': path.stat().st_size, 'md5': md5_file(path), 'rows': n,
                      'redistribution': 'git' if name.endswith('.json') else REDISTRIBUTION, 'content': content,
                      'source': ESM + SD2_NAME + '; ' + ESM + SD3_NAME,
                      'version_or_date': f'Supplementary Data 2 and 3, downloaded {args.downloaded}', 'citation': CITATION,
                      'licence': LICENCE, 'derived_from': f'{SD2_NAME}; {SD3_NAME}',
                      'how_made_or_checked': HOW if not name.endswith('.json') else 'scripts/45_stage_prmt7_szewczyk2020.py', 'original_name': ''})
    raw = [{'file': p.name, 'url': ESM + p.name, 'bytes': p.stat().st_size, 'md5': md5_file(p), 'sha256': sha256_file(p),
            'downloaded': args.downloaded, 'version': RAW_FILES[p.name], 'kept': RAW_KEPT} for p in (sd2, sd3)]
    if not args.manifest:   # staged elsewhere for review: the rows to add once the tables move into inputs/layers/
        for rows, name, cols in ((mrows, 'manifest_rows.tsv', list(mrows[0])), (raw, 'raw_files_rows.tsv', list(raw[0]))):
            save_table(pd.DataFrame(rows, columns=cols), out / name)
    else:
        update_manifest(args.manifest, mrows)
        upsert_tsv_rows(out / 'raw_files.tsv', raw)
        refresh_manifest_row(args.manifest, out / 'raw_files.tsv', 'layers/raw_files.tsv')
    print(json.dumps({k: v for k, v in summary.items() if k not in ('dropped', 'shared_peptides', 'inputs')}, indent=2))


if __name__ == '__main__':
    main()

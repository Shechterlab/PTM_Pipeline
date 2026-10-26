#!/usr/bin/env python3
"""Stage the PRMT1 knockdown arginine-methylation response (Hartel et al. 2019) on the canonical human proteome.

Hartel, Chew, Qin, Xu & Graham, "Deep protein methylation profiling by combined chemical and immunoaffinity approaches
reveals novel PRMT1 targets", Mol Cell Proteomics 18:2149-2164 (2019), doi:10.1074/mcp.RA119.001625 (CC BY 4.0 from
2020-09-15 per Crossref; PRIDE PXD012357, CC0). The quantification tables S8-S15, all in 153496_0_supp_351381_ptfjl4.xlsx
of the supplement mmc1.zip (the other file, 153496_0_supp_351382_ptfjl4.xlsx, holds Table S1 only), list the methyl
peptide groups quantified label-free (Minora, Proteome Discoverer 2.2) in HEK293T cells stably expressing an shRNA
against PRMT1 or a non-targeting control, two biological replicates each, per enrichment: high-pH SCX (S8 MMA, S9 DMA,
S15 peptides with both, S13 lysine) and immunoaffinity (S10 MMA, S11 ADMA, S12 SDMA, S14 pan-methyl-K). Each row gives
the peptide, its methyl marks with their ptmRS probabilities (every arginine mark is >= 91.5 %), the master proteins
and the peptide's span on them, the four normalised log2 abundances, log2 shPRMT1/shControl, and the Perseus test:
-log10 p, q, t and the flag 'Significant Perseus t-test q < 0.05'.

The test, as the methods give it ('Identification and Quantitation of Peptides'; 'Experimental Design and Statistical
Rationale'): log2 abundances, normalised to each sample's median, two-sample t test with permutation-based FDR in
Perseus 1.6.2.3, "a q-value FDR of 0.05 and S0 value of 0.5"; "Statistical significance was given to differential
peptide pairs with q-value < 0.05" and ptmRS localisation >= 75 %; only peptides quantified in all four samples were
tested. p and q are kept as the tables give them (p = 10^-(-log10 p)); nothing is recomputed. The lysine tables S13 and
S14 are read and counted, not staged (one lysine peptide changes, HMGN2 K40).

Each peptide is placed on the canonical sequence of the base accession of its first master protein (the accession the
Hartel records of scripts/00_stage_study_sites.py take): it must occur there exactly once. Its marks are read as
scripts/00 reads them (hartel_quant_mods: the residue on the peptide, the degree, the ptmRS probability). A placement is
position_verified when the first master is the canonical entry and the table's own start for it ('Peptide Position in
Master Proteins') is the start found, else placed_by_sequence. The table's own sites on the first master ('Protein
Accession and Methyl Site'), as offsets from its start there, must be the marks' offsets ('Modifications'); a peptide
whose two columns disagree (3 rows: S9 SARGTRGGGFGAHGRGR, S9 GEGGEFSVDRPIIDRPIRGR, S10 LSPLLPGGGGRGSR) is kept as
sites_differ_from_table and gives no site, as scripts/00 drops the same rows (position_inconsistent_with_peptide,
position_list_mismatch). A peptide changes at q < 0.05 (decreased for log2 < 0, increased for > 0; no
fold-change minimum, as S0 already enters the test); calls over the peptides carrying a site follow
scripts/peptide_calls.py (site_response: decreased, increased, mixed, no_significant_change, and the call from the
site's single-site peptides). call_rme1 and call_rme2 apply the same rule to the peptides that carry the site as
monomethyl or dimethyl arginine. call_source_flag applies it to the paper's own flag instead of q: the flag marks one
S10 peptide at q = 0.0504 (GGGFGGNDNFGRGGNFSGR, R12 me1, log2 1.87) that the q rule leaves unchanged (its site, HNRNPA1
R218, changes through another peptide, so the site calls agree). The flag reproduces the paper's peptide counts (119 MMA:
5 SCX + 114 IAP, 61 up and 58 down; DMA 2 SCX, 4 ADMA, 3 SDMA; 4 mixed). With 2 against 2 samples and S0 = 0.5 the
Perseus q weighs the fold change: 44 of the 132 flagged peptides have nominal p >= 0.05, and BH over the tables' p
values would keep 2 peptides (the summary counts both); a reader of these calls should know they are the paper's.

The peptide table keeps table_row (1-based among the table's data rows), the marks' offsets on the peptide,
methyl_positions and methyl_degrees (aligned, ';'), ptmRS, the four log2 abundances as given, log2_kd_control, p, q,
t, paper_significant and the peptide's call; other_master_entries lists the other entries a shared peptide's masters
name (the peptide is placed on the first only, as the records are).

Outputs (in --outdir, default $PTM_RESULTS_ROOT/hartel2019_prmt1_kd)
  hartel2019_prmt1_kd_peptides.tsv   every peptide of S8-S12 and S15 with its placement, marks, log2, p, q and call
  hartel2019_prmt1_kd_sites.tsv      one row per canonical methylarginine: peptide counts and the calls
  hartel2019_prmt1_kd_summary.json   input checksums, test, rule, counts, rows left out and why, the paper's counts
  manifest_rows.tsv                  the three tables' rows for inputs/MANIFEST.tsv (as layers/...)
  raw_files_rows.tsv                 the input's row for inputs/layers/raw_files.tsv (pointing at its MANIFEST row)
Nothing outside --outdir is written: --manifest is read only, to check the inputs' MD5s and take their URL and date.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

import peptide_calls as pc
from common import md5_file, parse_fasta, residue_is, save_table

SCRIPTS = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location('stage_study_sites', SCRIPTS / '00_stage_study_sites.py')
s00 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s00)

STUDY = 'Hartel 2019'
DOI = '10.1074/mcp.RA119.001625'
Q_CUTOFF, MIN_ABS_LOG2, STRICT = 0.05, 0.0, True
RULE = (f'peptide Perseus permutation-based FDR q < {Q_CUTOFF:g} as given (S0 = 0.5 inside the test; no fold-change '
        'minimum): decreased for log2 shPRMT1/shControl < 0, increased for > 0; sites by scripts/peptide_calls.site_response')
TEST = ('Perseus 1.6.2.3 two-sample t test with permutation-based FDR, q < 0.05, S0 = 0.5, on log2 abundances normalised to '
        'sample medians; peptides quantified in all four samples (2 shControl, 2 shPRMT1); ptmRS localisation >= 75 % '
        '(methods: Identification and Quantitation of Peptides; Experimental Design and Statistical Rationale)')
ARG_TABLES = ('S8', 'S9', 'S10', 'S11', 'S12', 'S15')
PLACED = ('position_verified', 'placed_by_sequence')   # the placements whose methyl positions make sites
LYS_TABLES = ('S13', 'S14')
COL = {'fold': 'log2 Fold Change (shPRMT1 / shControl)', 'q': 'Perseus t-test q value', 'mlogp': '-log10(Perseus p-value)',
       't': 'Perseus t-test t-statistic', 'flag': 'Significant Perseus t-test q < 0.05'}
LFQ = ('LFQ shControl_rep1', 'LFQ shControl_rep2', 'LFQ shPRMT1_rep1', 'LFQ shPRMT1_rep2')
DEGREE = {'Methyl': 'Rme1', 'Dimethyl': 'Rme2'}
# What the paper reports (abstract, Results, Discussion), to compare the staged counts with.
PAPER = {'arginine_sites_changed': 127, 'proteins': 78, 'mma_changed': 119, 'mma_changed_scx': 5, 'mma_changed_iap': 114,
         'mma_peptides_increased': 61, 'mma_peptides_decreased': 58, 'dma_peptides_increased': {'S9': 2, 'S11': 4, 'S12': 3},
         'mixed_peptides_changed': 4, 'lysine_sites_changed': 1,
         'note': ('The Results count the 119 as peptides ("Of the total 119 significantly changing methyl peptides, 5 came from SCX '
                  'enrichment, and 114 from IAP enrichment"); the Discussion calls them sites ("119 of 127, or 93.7%")')}

LICENCE = 'CC BY 4.0 (Crossref, from 2020-09-15); PRIDE PXD012357 CC0'
CITATION = 'Hartel NG et al., Mol Cell Proteomics 18:2149-2164 (2019), doi:10.1074/mcp.RA119.001625'
HOW = ('scripts/46_stage_hartel2019_prmt1_kd.py; peptides placed on reference/uniprot_human_reviewed_2026_03.fasta '
       '(a site is kept only where the canonical sequence has an arginine)')
MANIFEST_ROWS = {
    'hartel2019_prmt1_kd_peptides.tsv': (
        'Methyl peptides quantified in HEK293T shPRMT1 and shControl cells (Hartel et al. 2019 Tables S8-S12, S15), each placed on the '
        'canonical sequence, with methyl positions and degree, log2 shPRMT1/shControl, Perseus p and q as given, the paper\'s flag and the '
        'call (q < 0.05)'),
    'hartel2019_prmt1_kd_sites.tsv': (
        'One row per canonical methylarginine of those peptides: peptide counts, median log2, lowest q and the calls (decreased, increased, '
        'mixed, no_significant_change; from single-site peptides; by degree; by the paper\'s flag)'),
    'hartel2019_prmt1_kd_summary.json': (
        'Input checksums, test, rule, counts, the rows left out and why, and the paper\'s counts, of the PRMT1 knockdown staging'),
}


def read_tables(paths: list[Path]) -> tuple[dict[str, tuple[dict, list]], dict[str, str]]:
    """{table label: (header, rows)} for S8-S15, found by sheet name in whichever input holds it, and {label: file name}."""
    books = {}
    for path in paths:
        for name in s00.load_workbook(path, read_only=True).sheetnames:
            books[name] = path
    tables, source = {}, {}
    for sheet, label, kind, _enrichment, _cell, _nl in s00.HARTEL_SHEETS:
        if kind != 'quant':
            continue
        if sheet not in books:
            raise SystemExit(f'Hartel sheet {sheet!r} ({label}) not in {[str(p) for p in paths]}')
        tables[label] = s00.sheet_rows(books[sheet], sheet, 3)
        source[label] = books[sheet].name
        missing = [c for c in (*COL.values(), *LFQ, 'Peptide Sequence', 'Modifications', 'Master Protein Accessions',
                               'Peptide Position in Master Proteins', 'Protein Accession and Methyl Site') if c not in tables[label][0]]
        if missing:
            raise SystemExit(f'{label} ({sheet}): columns missing {missing}')
    return tables, source


def reported_starts(value, accession: str) -> list[int]:
    """Starts that 'Peptide Position in Master Proteins' gives on one master ('P09651 [197-218]; P0DMV8 [5-26]')."""
    starts, current = [], None
    for token in s00.text(value).split(';'):
        m = re.fullmatch(r'\s*(\S+)?\s*\[(\d+)-(\d+)\]\s*', token)
        if m:
            current = m.group(1) or current
            if current == accession:
                starts.append(int(m.group(2)))
    return starts


def place(peptide: str, masters: list[str], ppm, canonical: dict[str, str]) -> tuple[str, str, int | None]:
    """(status, canonical accession, 1-based start) of one peptide on its first master's canonical sequence."""
    acc = masters[0].split('-')[0]
    seq = canonical.get(acc)
    if seq is None:
        return 'accession_not_in_canonical_fasta', acc, None
    hits = [m.start() for m in re.finditer(f'(?={peptide})', seq)]
    if not hits:
        return 'peptide_not_in_canonical', acc, None
    if len(hits) > 1:
        return 'peptide_ambiguous', acc, None
    start = hits[0] + 1
    if '-' not in masters[0] and start in reported_starts(ppm, masters[0]):
        return 'position_verified', acc, start
    return 'placed_by_sequence', acc, start


def peptide_rows(tables: dict, canonical: dict[str, str], dropped: Counter, checks: Counter) -> list[dict]:
    rows = []
    for label in ARG_TABLES:
        H, data = tables[label]
        sheet = next(s for s, lab, *_ in s00.HARTEL_SHEETS if lab == label)
        enrichment = next(e for _s, lab, _k, e, *_ in s00.HARTEL_SHEETS if lab == label)
        for k, r in enumerate(data, start=1):
            peptide = s00.text(r[H['Peptide Sequence']]).upper()
            masters = s00.accessions(r[H['Master Protein Accessions']])
            mods, complete = s00.hartel_quant_mods(r[H['Modifications']])
            if not complete:
                checks['modification_count_mismatch'] += 1
            marks = sorted({(off, DEGREE[kind], prob) for res, off, kind, prob in mods if res == 'R' and kind in DEGREE})
            other = [m for m in mods if not (m[0] == 'R' and m[2] in DEGREE)]
            if other:
                checks['non_arginine_methyl_marks_ignored'] += len(other)
            offsets = [off for off, _, _ in marks]
            if len(set(offsets)) != len(offsets):
                raise SystemExit(f'{label} row {k}: an arginine marked with two degrees: {r[H["Modifications"]]}')
            if any(not 1 <= off <= len(peptide) or peptide[off - 1] != 'R' for off in offsets):
                raise SystemExit(f'{label} row {k}: a mark is not on an R of {peptide}: {r[H["Modifications"]]}')
            status, acc, start = place(peptide, masters, r[H['Peptide Position in Master Proteins']], canonical) if masters else \
                ('no_master_protein', '', None)
            positions = [start + off - 1 for off in offsets] if start else []
            if start and any(not residue_is(canonical, acc, p, 'R') for p in positions):
                raise SystemExit(f'{label} row {k}: a placed methyl position is not an R of {acc}')
            if start:
                if '-' not in masters[0] and status == 'placed_by_sequence':
                    checks['canonical_master_start_differs_from_table'] += 1
                # the table's sites on the first master, as offsets from the table's start there: must be the marks' offsets
                given = sorted(set(s00.hartel_protein_sites(r[H['Protein Accession and Methyl Site']]).get(masters[0], [])))
                starts = reported_starts(r[H['Peptide Position in Master Proteins']], masters[0])
                if not given or not starts:
                    checks['sites_not_checked:no_table_sites_for_first_master'] += 1
                elif any([g - s + 1 for g in given] == offsets for s in starts):
                    checks['sites_match_table'] += 1
                else:
                    checks['sites_differ_from_table'] += 1
                    status = 'sites_differ_from_table'
            fold, q, mlogp = (r[H[COL[c]]] for c in ('fold', 'q', 'mlogp'))
            if not all(isinstance(v, (int, float)) for v in (fold, q, mlogp)):
                raise SystemExit(f'{label} row {k}: fold change, q or -log10 p not a number')
            p = 10 ** -float(mlogp)
            call = pc.site_response([(float(q), float(fold), len(offsets))], Q_CUTOFF, MIN_ABS_LOG2, strict=STRICT)['call']
            others = sorted({m.split('-')[0] for m in masters} - {acc})
            rows.append({
                'table': label, 'sheet': sheet, 'enrichment': enrichment, 'table_row': k, 'peptide': peptide,
                'modifications': s00.text(r[H['Modifications']]), 'master_accessions': ';'.join(masters),
                'gene': s00.text(r[H['Gene Symbol']]).split(';')[0].strip() if 'Gene Symbol' in H else '',
                'reported_start': ';'.join(str(s) for s in reported_starts(r[H['Peptide Position in Master Proteins']], masters[0])) if masters else '',
                'reported_sites': s00.text(r[H['Protein Accession and Methyl Site']]),
                'canonical_UniProtAC': acc, 'placement': status, 'canonical_start': start if start else '',
                'other_master_entries': ';'.join(others),
                'methyl_offsets': ';'.join(str(o) for o in offsets),
                'methyl_positions': ';'.join(str(p) for p in positions),
                'methyl_degrees': ';'.join(d for _, d, _ in marks),
                'ptmrs': ';'.join('' if pr is None else f'{pr:g}' for _, _, pr in marks),
                'n_methylarginines': len(offsets),
                **{f'log2_{c.split(" ", 1)[1]}': r[H[c]] for c in LFQ},
                'log2_kd_control': float(fold), 'p': p, 'q': float(q), 't': r[H[COL['t']]],
                'paper_significant': s00.text(r[H[COL['flag']]]) == '+', 'call': call,
            })
            if status not in PLACED:
                dropped[f'peptides_without_sites:{status}'] += 1
                if rows[-1]['paper_significant']:
                    dropped[f'flagged_peptides_without_sites:{status}'] += 1
    return rows


def site_rows(peptides: pd.DataFrame, canonical: dict[str, str]) -> pd.DataFrame:
    groups = defaultdict(list)
    for r in peptides[peptides['placement'].isin(PLACED)].itertuples(index=False):
        for pos, degree in zip(r.methyl_positions.split(';'), r.methyl_degrees.split(';')):
            groups[(r.canonical_UniProtAC, int(pos))].append((r.q, r.log2_kd_control, r.n_methylarginines, degree,
                                                              r.paper_significant, r.table, r.gene))
    out = []
    for (acc, pos), vals in sorted(groups.items()):
        res = pc.site_response([(q, d, n) for q, d, n, *_ in vals], Q_CUTOFF, MIN_ABS_LOG2, strict=STRICT)
        by_degree = {deg: [(q, d, n) for q, d, n, g, *_ in vals if g == deg] for deg in ('Rme1', 'Rme2')}
        flag_down = sum(f and d < 0 for _, d, _, _, f, _, _ in vals)
        flag_up = sum(f and d > 0 for _, d, _, _, f, _, _ in vals)
        out.append({
            'canonical_UniProtAC': acc, 'position': pos, 'residue': canonical[acc][pos - 1],
            'gene': next((g for *_, g in vals if g), ''),
            'degrees': ';'.join(d for d in ('Rme1', 'Rme2') if by_degree[d]),
            'tables': ';'.join(t for t in ARG_TABLES if any(v[5] == t for v in vals)),
            'n_peptides': res['n_peptides'], 'n_decreased': res['n_decreased'], 'n_increased': res['n_increased'],
            'n_multi_site_peptides': res['n_multi_site_peptides'], 'log2_median': round(res['log2_median'], 3),
            'log2_strongest': res['log2_strongest'], 'q_min': res['adj_p_min'],
            'call': res['call'], 'call_single_site_peptides': res['call_single_site_peptides'],
            'call_rme1': pc.site_response(by_degree['Rme1'], Q_CUTOFF, MIN_ABS_LOG2, strict=STRICT)['call'],
            'call_rme2': pc.site_response(by_degree['Rme2'], Q_CUTOFF, MIN_ABS_LOG2, strict=STRICT)['call'],
            'call_source_flag': pc.peptide_call(flag_down, flag_up),
        })
    return pd.DataFrame(out)


def lysine_counts(tables: dict) -> dict:
    out = {}
    for label in LYS_TABLES:
        H, data = tables[label]
        flagged = [r for r in data if s00.text(r[H[COL['flag']]]) == '+']
        out[label] = {'rows': len(data), 'flagged_significant': len(flagged),
                      'flagged': [f'{s00.text(r[H["Gene Symbol"]])} {s00.text(r[H["Protein Accession and Methyl Site"]])}' for r in flagged]}
    return out


def manifest_inputs(manifest: Path | None, paths: list[Path]) -> dict[str, dict]:
    """Each input's row in inputs/MANIFEST.tsv (sources/raw/hartel2019/<name>), checked by MD5."""
    if manifest is None:
        return {}
    table = pd.read_csv(manifest, sep='\t', dtype=str, keep_default_na=False)
    rows = {}
    for p in paths:
        rel = f'sources/raw/hartel2019/{p.name}'
        hit = table[table['file'] == rel]
        if len(hit) != 1:
            raise SystemExit(f'{manifest}: no single row {rel}')
        row = hit.iloc[0].to_dict()
        if row['md5'] != md5_file(p):
            raise SystemExit(f'{p}: MD5 {md5_file(p)} differs from its MANIFEST row ({row["md5"]})')
        rows[p.name] = row
    return rows


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def write_rows(outdir: Path, paths: list[Path], tables_from: dict[str, str], inputs_rows: dict[str, dict]) -> None:
    """manifest_rows.tsv (columns of inputs/MANIFEST.tsv) and raw_files_rows.tsv (columns of inputs/layers/raw_files.tsv)."""
    used = [p for p in paths if p.name in tables_from.values()]
    derived = '; '.join(f'{p.name} (MD5 {md5_file(p)})' for p in used)
    source = '; '.join(sorted({inputs_rows[p.name]['source'] for p in used})) if inputs_rows else ''
    cols = ['file', 'bytes', 'md5', 'rows', 'redistribution', 'content', 'source', 'version_or_date', 'citation', 'licence',
            'derived_from', 'how_made_or_checked', 'original_name']
    rows = []
    for name, content in MANIFEST_ROWS.items():
        path = outdir / name
        n = '' if name.endswith('.json') else sum(1 for _ in open(path, encoding='utf-8')) - 1
        rows.append({'file': f'layers/{name}', 'bytes': path.stat().st_size, 'md5': md5_file(path), 'rows': n, 'redistribution': 'git',
                     'content': content, 'source': source, 'version_or_date': 'Supplementary Tables S8-S15 as published',
                     'citation': CITATION, 'licence': LICENCE, 'derived_from': derived,
                     'how_made_or_checked': HOW if not name.endswith('.json') else 'scripts/46_stage_hartel2019_prmt1_kd.py',
                     'original_name': ''})
    pd.DataFrame(rows, columns=cols).to_csv(outdir / 'manifest_rows.tsv', sep='\t', index=False, lineterminator='\n')
    raw = []
    for p in used:
        row = inputs_rows.get(p.name, {})
        date = re.search(r'Downloaded (\d{4}-\d\d-\d\d)', row.get('how_made_or_checked', ''))
        raw.append({'file': p.name, 'url': f'{row.get("source", "")} (member {p.name})', 'bytes': p.stat().st_size, 'md5': md5_file(p),
                    'sha256': sha256_file(p), 'downloaded': date.group(1) if date else '',
                    'version': 'Supplementary Tables S2-S15 (mmc1.zip), as published',
                    'kept': f'zenodo, as sources/raw/hartel2019/{p.name} (its inputs/MANIFEST.tsv row)'})
    pd.DataFrame(raw, columns=['file', 'url', 'bytes', 'md5', 'sha256', 'downloaded', 'version', 'kept']).to_csv(
        outdir / 'raw_files_rows.tsv', sep='\t', index=False, lineterminator='\n')


def main() -> None:
    results_root = Path(os.environ.get('PTM_RESULTS_ROOT', 'PTM_results'))
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--input', nargs='+', required=True,
                    help='Hartel 2019 supplement(s); the one holding Tables S8-S15 is 153496_0_supp_351381_ptfjl4.xlsx.')
    ap.add_argument('--fasta', default=os.environ.get('PTM_CANONICAL_FASTA', ''),
                    help='UniProt reviewed canonical FASTA (inputs/reference/uniprot_human_reviewed_2026_03.fasta; PTM_CANONICAL_FASTA).')
    ap.add_argument('--manifest', default='', help='inputs/MANIFEST.tsv, read only: checks each input against its sources/raw/hartel2019 row.')
    ap.add_argument('--outdir', default=str(results_root / 'hartel2019_prmt1_kd'))
    args = ap.parse_args()
    if not args.fasta:
        raise SystemExit('--fasta is required (or set PTM_CANONICAL_FASTA)')

    paths = [Path(p) for p in args.input]
    inputs_rows = manifest_inputs(Path(args.manifest) if args.manifest else None, paths)
    tables, tables_from = read_tables(paths)
    canonical = parse_fasta(args.fasta)
    dropped, checks = Counter(), Counter()
    peptides = pd.DataFrame(peptide_rows(tables, canonical, dropped, checks))
    sites = site_rows(peptides, canonical)
    if sites['residue'].ne('R').any():
        raise SystemExit('placed methylation sites that are not arginines on the canonical sequence')

    outdir = Path(args.outdir)
    save_table(peptides, outdir / 'hartel2019_prmt1_kd_peptides.tsv')
    save_table(sites, outdir / 'hartel2019_prmt1_kd_sites.tsv')

    placed = peptides[peptides['placement'].isin(PLACED)]
    changed_sites = sites[sites['call'] != pc.NO_CHANGE]
    flagged_sites = sites[sites['call_source_flag'] != pc.NO_CHANGE]
    disagree = peptides[(peptides['call'] != pc.NO_CHANGE) != peptides['paper_significant']]
    shared = placed[placed['other_master_entries'] != '']
    fractions = placed.groupby(['table', 'peptide', 'modifications']).size()
    summary = {
        'study': STUDY, 'doi': DOI, 'pride': 'PXD012357', 'licence': LICENCE,
        'comparison': 'HEK293T shPRMT1 / shControl (stable shRNA), 2 biological replicates each, per enrichment',
        'test': TEST, 'rule': RULE,
        'p_and_q': 'as given in Tables S8-S15: p = 10^-(-log10(Perseus p-value)), q = Perseus t-test q value; nothing recomputed',
        'inputs': {p.name: {'md5': md5_file(p), 'tables': [t for t, f in tables_from.items() if f == p.name],
                            'manifest_row': f'sources/raw/hartel2019/{p.name}' if p.name in inputs_rows else ''} for p in paths},
        'canonical_fasta': {'file': Path(args.fasta).name, 'md5': md5_file(args.fasta), 'sequences': len(canonical)},
        'rows_read': {t: len(tables[t][1]) for t in (*ARG_TABLES, *LYS_TABLES)},
        'not_staged': {'lysine_tables': lysine_counts(tables),
                       'why': 'lysine methylation (S13 SCX mono/di/tri-methyl K, S14 pan-methyl-K IAP): no methylarginine'},
        'placement': dict(sorted(Counter(peptides['placement']).items())),
        'dropped': dict(sorted(dropped.items())),
        'checks': dict(sorted(checks.items())),
        'placed_peptides_shared_with_other_entries': int(len(shared)),
        'same_peptide_and_marks_in_two_scx_fractions': int((fractions > 1).sum()),
        'peptides': {t: int((peptides['table'] == t).sum()) for t in ARG_TABLES},
        'peptide_calls': {t: dict(sorted(Counter(g['call']).items())) for t, g in peptides.groupby('table', sort=False)},
        'peptides_paper_flag': {t: int(g['paper_significant'].sum()) for t, g in peptides.groupby('table', sort=False)},
        'q_and_p': {
            'flagged_peptides_with_p_at_least_0.05': int((peptides['paper_significant'] & (peptides['p'] >= 0.05)).sum()),
            'flagged_peptides_p_max': float(peptides.loc[peptides['paper_significant'], 'p'].max()),
            'unflagged_peptides_with_p_below_0.05': int((~peptides['paper_significant'] & (peptides['p'] < 0.05)).sum()),
            'peptides_with_q_0': int((peptides['q'] == 0).sum()),
            'bh_over_p_within_table_below_0.05': {t: sum(q < Q_CUTOFF for q in pc.bh(g['p'].tolist()))
                                                  for t, g in peptides.groupby('table', sort=False)},
            'note': ('Perseus q (permutation-based FDR with S0 = 0.5, 2 against 2 samples) is not BH of p: it weighs the fold '
                     'change, so peptides with large log2 pass at nominal p above 0.05, and BH over the given p would keep few'),
        },
        'flag_and_q_disagree': [f'{r.table} {r.peptide} {r.modifications} q={r.q:g} flag={"+" if r.paper_significant else ""}'
                                for r in disagree.itertuples(index=False)],
        'sites': int(len(sites)), 'site_proteins': int(sites['canonical_UniProtAC'].nunique()),
        'site_calls': dict(sorted(Counter(sites['call']).items())),
        'site_calls_single_site_peptides': dict(sorted(Counter(sites['call_single_site_peptides'].replace('', '(none)')).items())),
        'site_calls_rme1': dict(sorted(Counter(sites['call_rme1'].replace('', '(not Rme1)')).items())),
        'site_calls_rme2': dict(sorted(Counter(sites['call_rme2'].replace('', '(not Rme2)')).items())),
        'site_calls_source_flag': dict(sorted(Counter(sites['call_source_flag']).items())),
        'changed_sites': {'q_rule': int(len(changed_sites)), 'q_rule_proteins': int(changed_sites['canonical_UniProtAC'].nunique()),
                          'paper_flag': int(len(flagged_sites)), 'paper_flag_proteins': int(flagged_sites['canonical_UniProtAC'].nunique())},
        'paper': PAPER,
    }
    (outdir / 'hartel2019_prmt1_kd_summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    write_rows(outdir, paths, tables_from, inputs_rows)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()

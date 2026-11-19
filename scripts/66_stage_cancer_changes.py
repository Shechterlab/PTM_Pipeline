#!/usr/bin/env python3
"""Stage the somatic protein changes at every candidate residue of the context-mark site sets (lysine, serine, threonine,
tyrosine), in TCGA tumors and DepMap cell lines, and each protein's background counts.

The arginine layer (scripts/44b_extract_tcga_mc3.py, scripts/44_stage_depmap.py) keeps the changes at every arginine of a
canonical sequence, a methylarginine site or not; the website then shows them at the sites. This script applies the same
rules to the residues the other site sets modify (K for the lysine and methyllysine sets, S, T and Y for phosphorylation),
so a site table can be read against every candidate residue or against its known sites only, and adds the counts a site's
changes can be compared with: for each protein and residue, how many residues of that kind it has and how many changes
fall on them, and the same over all its residues.

Sources (the files of the arginine layer, same MD5)
  TCGA MC3 public MAF (mc3.v0.2.8.PUBLIC.maf.gz, GDC open access; Ellrott et al. 2018, Cell Syst 6:271) and the GDC cases
    of TCGA (each patient's project).
  DepMap Public 24Q4 OmicsSomaticMutations.csv (CC BY 4.0; doi:10.25452/figshare.plus.27993248), and the genes its CRISPR
    screens hold (the header of CRISPRGeneEffect.csv), which decide an entry's Entrez gene as in scripts/44.
  UniProt reviewed human entries: their Ensembl transcripts and displayed isoforms (--uniprot-ensembl, --uniprot-isoforms;
    the release is read from the files' download record, --uniprot-release), and the reference proteome of the release
    tables (UniProt 2026_03, inputs/reference/).

Rules (as for arginines)
  TCGA: FILTER PASS, and in TCGA-OV and TCGA-LAML also FILTER wga alone (Ellrott et al. 2018). A missense or stop change is
    read on a transcript UniProt assigns to a reviewed entry's displayed isoform (all_effects); when the variant falls in no
    such transcript, the Missense_Mutation or Nonsense_Mutation of the MAF transcript (HGVSp_Short) is placed on the entry
    its SWISSPROT field names, else on the entries of its gene symbol (the Entrez gene of the entry's primary gene name in
    org.Hs.eg.db). A change counts only where the 2026_03 canonical sequence has the changed residue at its position.
    One row per tumor sample and change (accession, position, new residue).
  DepMap: a call on the canonical isoform (UniprotID unlabelled or -1; else the entry VepSwissprot names) whose residue the
    canonical sequence has at its position; one row per cell line and change.
  Background: every change kept by these rules, at any residue.

Outputs (in --outdir; gzip without a timestamp, so the same input gives the same bytes)
  tcga_mc3_changes_<X>.tsv.gz              per residue X: one row per tumor sample and change at a canonical X, the
                                           columns of tcga_mc3_arginine_changes.tsv (alt * is a stop)
  depmap_24q4_changes_<X>.tsv.gz           per residue X: one row per cell line and change, the columns of
                                           depmap_24q4_arginine_changes.tsv
  cancer_change_protein_background.tsv.gz  per reference entry and residue class (R, K, S, T, Y, and all residues): the
                                           residues of that class in the canonical sequence (its length for all), and the
                                           TCGA and DepMap changes on them, the samples or lines carrying one and the
                                           positions changed
  depmap_24q4_change_models.tsv            every cell line with a placed change (any residue): name, Oncotree lineage and
                                           primary disease (DepMap Model.csv), with --depmap-models
  cancer_changes_summary.json              inputs (MD5), the rules, the denominators (tumor samples, cell lines with
                                           mutation calls) and the counts of what was kept and dropped
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import importlib.util
import io
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from common import (entrez_by_symbol, entrez_of_entries, md5_file, own_entrez, parse_fasta, refresh_manifest_row, update_manifest,
                    upsert_tsv_rows)

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location('extract_tcga_mc3', HERE / '44b_extract_tcga_mc3.py')
s44b = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s44b)
_spec = importlib.util.spec_from_file_location('stage_depmap', HERE / '44_stage_depmap.py')
s44 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s44)

AA = 'ACDEFGHIKLMNPQRSTVWY'
THREE = s44b.THREE | {'Arg': 'R'}     # scripts/44b reads only changes from arginine, so its table lacks Arg
CLASSES = ['R', 'K', 'S', 'T', 'Y']
EFFECT = re.compile(r'^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter)$')
SHORT = re.compile(rf'^p\.([{AA}])(\d+)([{AA}*])$')
DEPMAP = re.compile(rf'^p\.([{AA}])(\d+)([{AA}]|Ter)$')
TCGA_COLUMNS = ['canonical_UniProtAC', 'position', 'alt', 'project_id', 'sample_id', 'transcript_id', 'mapped_by', 'filter',
                'tumor_alt_count', 'tumor_ref_count']
DEPMAP_COLUMNS = ['canonical_UniProtAC', 'position', 'alt', 'protein_change', 'model_id', 'allele_fraction', 'hotspot', 'likely_lof',
                  'am_pathogenicity', 'am_class']
BACKGROUND_COLUMNS = ['canonical_UniProtAC', 'residue', 'residues', 'tcga_changes', 'tcga_samples', 'tcga_positions',
                      'depmap_changes', 'depmap_lines', 'depmap_positions']


def write_gz(df: pd.DataFrame, path: Path) -> None:
    """A TSV, gzip-compressed with mtime 0 and no file name in the header, so the bytes depend only on the table."""
    buf = io.StringIO()
    df.to_csv(buf, sep='\t', index=False, lineterminator='\n')
    with open(path, 'wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0, compresslevel=9) as gz:
        gz.write(buf.getvalue().encode('utf-8'))


class Placer:
    """Places a change read on the MAF's own transcript: the entry its SWISSPROT field names, else the entries of its gene
    symbol, as scripts/44_stage_depmap.py does for arginines."""

    def __init__(self, fasta: str, entries: str, symbol_entrez: str, canonical: dict, screened: set[str]):
        self.entry_acc = {}
        for line in open(fasta, encoding='utf-8'):
            if line.startswith('>'):
                parts = line[1:].split('|')
                if len(parts) >= 3:
                    self.entry_acc[parts[2].split()[0]] = parts[1]
        self.entrez_of, self.own = entrez_by_symbol(symbol_entrez), own_entrez(symbol_entrez)
        # an entry's Entrez gene among those DepMap screened; an entry whose gene was not screened takes the Entrez gene
        # org.Hs.eg.db gives its name as its own (scripts/44_stage_depmap.py)
        acc_entrez, acc_gene, _ = entrez_of_entries(entries, symbol_entrez, allowed=screened, reference=canonical)
        self.accs_of = defaultdict(list)
        for acc, e in acc_entrez.items():
            self.accs_of[e].append(acc)
        for acc, g in acc_gene.items():
            if acc not in acc_entrez:
                for e in self.entrez_of.get(g, []):
                    if (g, e) in self.own:
                        self.accs_of[e].append(acc)
        self.canonical = canonical

    def candidates(self, swissprot: str, symbol: str) -> list[str]:
        if swissprot in self.entry_acc:
            return [self.entry_acc[swissprot]]
        es = [e for e in self.entrez_of.get(symbol, []) if (symbol, e) in self.own] or self.entrez_of.get(symbol, [])
        return sorted({a for e in es for a in self.accs_of.get(e, [])})


def on_residue(canonical: dict, acc: str, pos: int, ref: str) -> bool:
    seq = canonical.get(acc)
    return bool(seq) and 1 <= pos <= len(seq) and seq[pos - 1] == ref


def tcga(args, canonical: dict, placer: Placer, residues: list[str], counts: Counter):
    """Kept TCGA changes: (rows at the residues asked for, every placed change as (acc, ref, pos, alt, sample), tumor samples
    per project)."""
    project_of = {}
    with open(args.cases, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            project_of[r['submitter_id']] = r['project.project_id']
    canon, xref = s44b.canonical_transcripts(args.uniprot_ensembl, args.uniprot_isoforms)
    counts.update({f'uniprot_{k}': v for k, v in xref.items()})
    samples = defaultdict(set)
    raw = []          # every kept change on a canonical transcript or placed from the MAF transcript, before the residue check
    with gzip.open(args.maf, 'rt', encoding='utf-8', newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        ix = {c: k for k, c in enumerate(next(reader))}
        for r in reader:
            counts['tcga_maf_rows'] += 1
            barcode = r[ix['Tumor_Sample_Barcode']]
            patient, sample = barcode[:12], barcode[:15]
            project = project_of.get(patient, '')
            samples[project].add(sample)
            filt = r[ix['FILTER']]
            if not s44b.kept(filt, project):
                continue
            if not project:
                counts['tcga_kept_calls_patient_not_in_gdc'] += 1
                continue
            common = (project, sample, filt, r[ix['Hugo_Symbol']], r[ix['Chromosome']], r[ix['Start_Position']], r[ix['Tumor_Seq_Allele2']],
                      r[ix['t_alt_count']], r[ix['t_ref_count']])
            effects = [e.split(',') for e in r[ix['all_effects']].split(';') if e]
            on_canonical = [e for e in effects if len(e) >= 4 and e[3] in canon]
            if on_canonical:
                for e in on_canonical:
                    m = EFFECT.match(e[2])
                    if (not m or m.group(1) not in THREE or m.group(1) == 'Ter' or m.group(3) not in THREE
                            or not any(c in e[1] for c in s44b.KEEP_CONSEQUENCES)):
                        continue
                    ref, pos, alt = THREE[m.group(1)], int(m.group(2)), THREE[m.group(3)]
                    if alt == ref:
                        continue
                    for acc in sorted(canon[e[3]]):
                        raw.append((acc, ref, pos, alt, e[3], 'canonical_transcript') + common)
                continue
            if r[ix['Variant_Classification']] not in s44b.KEEP_TYPES:
                continue
            m = SHORT.match(r[ix['HGVSp_Short']])
            if not m:
                continue
            ref, pos, alt = m.group(1), int(m.group(2)), m.group(3)
            cands = placer.candidates(r[ix['SWISSPROT']].replace('.', ''), r[ix['Hugo_Symbol']])
            if not cands:
                counts['tcga_maf_transcript_changes_gene_without_entry'] += 1
            for acc in cands:
                raw.append((acc, ref, pos, alt, r[ix['Transcript_ID']], 'maf_transcript') + common)
    # as scripts/44b then 44: one transcript per (sample, entry, change, variant) -- the last by transcript ID -- then the
    # residue check, then one row per (entry, position, new residue, sample)
    best = {}
    for x in sorted(raw, key=lambda x: x[4]):
        best[(x[7], x[0], x[9], x[1], x[2], x[3], x[10], x[11], x[12])] = x
    rows, placed, seen = [], [], set()
    for x in sorted(best.values(), key=lambda x: (x[9], x[0], x[2], x[1], x[3], x[7], x[10], int(x[11]), x[12])):
        acc, ref, pos, alt = x[0], x[1], x[2], x[3]
        if not on_residue(canonical, acc, pos, ref):
            counts['tcga_changes_residue_not_on_canonical'] += 1
            continue
        key = (acc, pos, alt, x[7])
        if key in seen:
            continue
        seen.add(key)
        placed.append((acc, ref, pos, alt, x[7]))
        if ref in residues:
            rows.append({'canonical_UniProtAC': acc, 'position': pos, 'alt': alt, 'project_id': x[6], 'sample_id': x[7],
                         'transcript_id': x[4], 'mapped_by': x[5], 'filter': x[8], 'tumor_alt_count': x[13], 'tumor_ref_count': x[14]})
    counts['tcga_changes_placed_all_residues'] = len(placed)
    return rows, placed, samples


def depmap(args, canonical: dict, residues: list[str], counts: Counter):
    mut = pd.read_csv(args.depmap_mutations, dtype=str, keep_default_na=False,
                      usecols=['ProteinChange', 'UniprotID', 'VepSwissprot', 'ModelID', 'AF', 'Hotspot', 'LikelyLoF', 'AMPathogenicity', 'AMClass'])
    lines = mut['ModelID'].nunique()
    parts = mut['ProteinChange'].str.extract(DEPMAP)
    m = mut[parts[0].notna()].copy()
    m[['ref', 'position', 'alt']] = parts[parts[0].notna()].to_numpy()
    m['alt'] = m['alt'].replace('Ter', '*')
    m = m[m['alt'] != m['ref']]
    counts['depmap_missense_and_stop_calls'] = len(m)
    iso = m['UniprotID'].str.extract(r'^([A-Z0-9]+?)(?:-(\d+))?$')
    base = iso[0].where(iso[0].notna() & (iso[0] != ''), m['VepSwissprot'].str.split('.').str[0])
    keep = iso[1].isna() | (iso[1] == '1')
    counts['depmap_calls_on_other_isoforms'] = int((~keep).sum())
    m, base = m[keep].copy(), base[keep]
    m['canonical_UniProtAC'] = base
    counts['depmap_calls_not_in_reference'] = int((~m['canonical_UniProtAC'].isin(canonical)).sum())
    m = m[m['canonical_UniProtAC'].isin(canonical)]
    m['position'] = m['position'].astype(int)
    ok = [on_residue(canonical, a, p, r) for a, p, r in zip(m['canonical_UniProtAC'], m['position'], m['ref'])]
    counts['depmap_calls_residue_not_on_canonical'] = int(len(m) - sum(ok))
    m = m[ok].drop_duplicates(['canonical_UniProtAC', 'position', 'alt', 'ModelID'])
    counts['depmap_changes_placed_all_residues'] = len(m)
    placed = list(zip(m['canonical_UniProtAC'], m['ref'], m['position'], m['alt'], m['ModelID']))
    t = m[m['ref'].isin(residues)]
    rows = pd.DataFrame({'canonical_UniProtAC': t['canonical_UniProtAC'], 'position': t['position'], 'alt': t['alt'],
                         'protein_change': t['ref'] + t['position'].astype(str) + t['alt'], 'model_id': t['ModelID'],
                         'allele_fraction': t['AF'], 'hotspot': t['Hotspot'] == 'True', 'likely_lof': t['LikelyLoF'] == 'True',
                         'am_pathogenicity': t['AMPathogenicity'], 'am_class': t['AMClass'], 'ref': t['ref']})
    return rows, placed, lines


def background(canonical: dict, tplaced, dplaced) -> pd.DataFrame:
    agg = defaultdict(lambda: [0, set(), set(), 0, set(), set()])
    for src, placed in ((0, tplaced), (3, dplaced)):
        for acc, ref, pos, alt, who in placed:
            for cls in ([ref] if ref in CLASSES else []) + ['all']:
                a = agg[(acc, cls)]
                a[src] += 1
                a[src + 1].add(who)
                a[src + 2].add(pos)
    rows = []
    for acc in sorted(canonical):
        seq = canonical[acc]
        for cls in CLASSES + ['all']:
            a = agg.get((acc, cls), [0, set(), set(), 0, set(), set()])
            rows.append([acc, cls, len(seq) if cls == 'all' else seq.count(cls), a[0], len(a[1]), len(a[2]), a[3], len(a[4]), len(a[5])])
    return pd.DataFrame(rows, columns=BACKGROUND_COLUMNS)


TC_LICENSE = s44.TC_LIC.replace('licence', 'license')     # the MANIFEST column keeps its name
NAMES = {'R': 'arginine', 'K': 'lysine', 'S': 'serine', 'T': 'threonine', 'Y': 'tyrosine'}
HOW = ('scripts/66_stage_cancer_changes.py, the rules of the arginine layer (scripts/44b, scripts/44; rerun on arginines it gives '
       'tcga_mc3_arginine_changes.tsv and depmap_24q4_arginine_changes.tsv byte for byte); checked against '
       'reference/uniprot_human_reviewed_2026_03.fasta (a change counts only where the canonical sequence has the changed residue)')


def gz_rows(path: Path) -> int:
    with gzip.open(path, 'rt', encoding='utf-8') as fh:
        return sum(1 for _ in fh) - 1


def sha256(p: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(p, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def write_manifest(manifest: Path, out: Path, args, residues: list[str], downloaded: str) -> None:
    """The staged tables' rows in inputs/MANIFEST.tsv, and the downloads this staging adds in the folder's raw_files.tsv
    (the MC3 MAF and OmicsSomaticMutations.csv are those of scripts/44, same MD5, and keep their rows)."""
    folder = out.name
    tc_ver = f'MC3 public MAF v0.2.8 (GDC publication page, downloaded 2026-10-25, same MD5); GDC cases fetched {downloaded}'
    tc_src = s44.TC_SRC.replace('scripts/44b_extract_tcga_mc3.py', 'scripts/66_stage_cancer_changes.py')
    rows = []
    for res in residues:
        name = NAMES[res]
        rows.append({'file': f'{folder}/tcga_mc3_changes_{res}.tsv.gz', 'content': f'TCGA MC3 somatic missense and stop changes at every {name} '
                     f'({res}) of reviewed canonical proteins, a known site or not: one row per tumor sample and change (residue * a stop), calls '
                     'with FILTER PASS (OV and LAML also wga alone, as Ellrott et al. 2018), with the Ensembl transcript, its FILTER and the '
                     'tumor read counts; the columns of tcga_mc3_arginine_changes.tsv (gzip)', 'source': tc_src, 'version_or_date': tc_ver,
                     'citation': s44.TC_CIT, 'licence': TC_LICENSE, 'derived_from': 'mc3.v0.2.8.PUBLIC.maf.gz; gdc_tcga_cases (fetched ' + downloaded + ')'})
        rows.append({'file': f'{folder}/depmap_24q4_changes_{res}.tsv.gz', 'content': f'DepMap somatic missense and stop changes at every {name} '
                     f'({res}) of reviewed canonical proteins, a known site or not: one row per cell line and change, with allele fraction, '
                     'hotspot, likely loss of function and AlphaMissense as DepMap gives them; the columns of depmap_24q4_arginine_changes.tsv (gzip)',
                     'source': s44.DM_SRC, 'version_or_date': s44.DM_VER, 'citation': s44.DM_CIT + '; AlphaMissense: Cheng et al., Science 381:eadg7492 (2023)',
                     'licence': 'CC BY 4.0', 'derived_from': 'OmicsSomaticMutations.csv'})
    rows.append({'file': f'{folder}/cancer_change_protein_background.tsv.gz', 'content': 'Per reviewed canonical entry and residue class (R, K, '
                 'S, T, Y, all): the residues of that class in the sequence (its length for all), and the TCGA MC3 and DepMap missense and stop '
                 'changes on them, the tumor samples or cell lines carrying one and the positions changed, every change counted by the rules '
                 'of the change tables (gzip)', 'source': tc_src + '; ' + s44.DM_SRC, 'version_or_date': tc_ver + '; ' + s44.DM_VER,
                 'citation': s44.TC_CIT + '; ' + s44.DM_CIT, 'licence': TC_LICENSE + '; CC BY 4.0 (DepMap)',
                 'derived_from': 'mc3.v0.2.8.PUBLIC.maf.gz; OmicsSomaticMutations.csv; reference/uniprot_human_reviewed_2026_03.fasta'})
    for r in rows:
        path = out / r['file'].split('/', 1)[1]
        r.update({'bytes': path.stat().st_size, 'md5': md5_file(path), 'rows': gz_rows(path), 'redistribution': 'git', 'how_made_or_checked': HOW})
    models = out / 'depmap_24q4_change_models.tsv'
    if models.exists():
        rows.append({'file': f'{folder}/depmap_24q4_change_models.tsv', 'bytes': models.stat().st_size, 'md5': md5_file(models),
                     'rows': sum(1 for _ in open(models, encoding='utf-8')) - 1, 'redistribution': 'git',
                     'content': 'Every DepMap cell line with a missense or stop change kept by the rules of the change tables (any residue): '
                     'its name, Oncotree lineage and primary disease', 'source': s44.DM_SRC, 'version_or_date': s44.DM_VER,
                     'citation': s44.DM_CIT, 'licence': 'CC BY 4.0', 'derived_from': 'Model.csv; OmicsSomaticMutations.csv', 'how_made_or_checked': HOW})
    summary = out / 'cancer_changes_summary.json'
    rows.append({'file': f'{folder}/cancer_changes_summary.json', 'bytes': summary.stat().st_size, 'md5': md5_file(summary), 'rows': '',
                 'redistribution': 'git', 'content': 'Inputs (MD5), rules, denominators (TCGA tumor samples, DepMap lines with mutation calls) '
                 'and counts of the staging of the changes at K, S, T and Y and of the per-protein background',
                 'source': 'the downloads above', 'version_or_date': tc_ver + '; ' + s44.DM_VER, 'citation': '', 'licence': '',
                 'derived_from': '', 'how_made_or_checked': 'scripts/66_stage_cancer_changes.py'})
    update_manifest(manifest, rows)
    raw = []
    for path, url, ver in ((Path(args.cases), s44b.CASES_URL, f'GDC cases of TCGA, fetched {downloaded}'),
                           (Path(args.uniprot_ensembl), s44b.UNIPROT_ENSEMBL_URL, f'UniProt {args.uniprot_release}, Ensembl cross-references'),
                           (Path(args.uniprot_isoforms), s44b.UNIPROT_ISOFORMS_URL, f'UniProt {args.uniprot_release}, alternative products'),
                           (Path(args.crispr_gene_effect), 'https://ndownloader.figshare.com/files/51064667 (bytes 0 to the end of the first line)',
                            'DepMap Public 24Q4, the header line of CRISPRGeneEffect.csv (the genes screened)')):
        raw.append({'file': path.name, 'url': url, 'bytes': path.stat().st_size, 'md5': md5_file(path), 'sha256': sha256(path),
                    'downloaded': downloaded, 'version': ver, 'kept': 'release_asset'})
    upsert_tsv_rows(out / 'raw_files.tsv', raw)
    refresh_manifest_row(manifest, out / 'raw_files.tsv', f'{folder}/raw_files.tsv')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--maf', required=True, help='mc3.v0.2.8.PUBLIC.maf.gz')
    ap.add_argument('--cases', required=True, help='GDC cases of TCGA (TSV: submitter_id, project.project_id, project.name)')
    ap.add_argument('--depmap-mutations', required=True, help='DepMap 24Q4 OmicsSomaticMutations.csv')
    ap.add_argument('--crispr-gene-effect', required=True, help='DepMap 24Q4 CRISPRGeneEffect.csv, or its header line: the genes screened, '
                    'which decide an entry\'s Entrez gene as in scripts/44')
    ap.add_argument('--depmap-models', default='', help='DepMap 24Q4 Model.csv: write depmap_24q4_change_models.tsv')
    ap.add_argument('--uniprot-ensembl', required=True, help='UniProt reviewed human: Entry, Entry Name, Ensembl (TSV, gzip)')
    ap.add_argument('--uniprot-isoforms', required=True, help='UniProt reviewed human: Entry, Alternative products (TSV, gzip)')
    ap.add_argument('--uniprot-release', required=True, help='the UniProt release the two UniProt files were downloaded from')
    ap.add_argument('--fasta', default='inputs/reference/uniprot_human_reviewed_2026_03.fasta')
    ap.add_argument('--entries', default='inputs/reference/uniprot_human_reviewed_2026_03_entries.tsv')
    ap.add_argument('--symbol-entrez', default='inputs/disease/orghsegdb_3.22.0_symbol_entrez.tsv')
    ap.add_argument('--residues', default='KSTY', help='residues to write change tables for (default KSTY)')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--manifest', default='', help='inputs/MANIFEST.tsv: write the rows of the tables staged (and of the new downloads in '
                    'raw_files.tsv beside them), as scripts/44 does; empty to leave both alone')
    ap.add_argument('--downloaded', default=dt.date.today().isoformat(), help='download date of the GDC cases and UniProt files')
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    residues = [r for r in args.residues.upper() if r in AA]
    canonical = parse_fasta(args.fasta)
    counts = Counter()

    screened = {m.group(1) for c in pd.read_csv(args.crispr_gene_effect, nrows=0).columns if (m := re.match(r'^.+ \((\d+)\)$', str(c)))}
    placer = Placer(args.fasta, args.entries, args.symbol_entrez, canonical, screened)
    trows, tplaced, samples = tcga(args, canonical, placer, residues, counts)
    drows, dplaced, lines = depmap(args, canonical, residues, counts)

    tt = pd.DataFrame(trows, columns=TCGA_COLUMNS)
    ref_of = {(a, p): canonical[a][p - 1] for a, p in zip(tt['canonical_UniProtAC'], tt['position'])}
    written = {}
    for res in residues:
        t = tt[[ref_of[(a, p)] == res for a, p in zip(tt['canonical_UniProtAC'], tt['position'])]]
        t = t.sort_values(['canonical_UniProtAC', 'position', 'alt', 'project_id', 'sample_id'], kind='stable')
        d = drows[drows['ref'] == res].drop(columns='ref').sort_values(['canonical_UniProtAC', 'position', 'alt', 'model_id'], kind='stable')
        write_gz(t, out / f'tcga_mc3_changes_{res}.tsv.gz')
        write_gz(d, out / f'depmap_24q4_changes_{res}.tsv.gz')
        written[res] = {'tcga_changes': len(t), 'tcga_positions': int(t[['canonical_UniProtAC', 'position']].drop_duplicates().shape[0]),
                        'tcga_samples': int(t['sample_id'].nunique()), 'depmap_changes': len(d),
                        'depmap_positions': int(d[['canonical_UniProtAC', 'position']].drop_duplicates().shape[0]),
                        'depmap_lines': int(d['model_id'].nunique()),
                        'md5': {f'tcga_mc3_changes_{res}.tsv.gz': md5_file(out / f'tcga_mc3_changes_{res}.tsv.gz'),
                                f'depmap_24q4_changes_{res}.tsv.gz': md5_file(out / f'depmap_24q4_changes_{res}.tsv.gz')}}
    bg = background(canonical, tplaced, dplaced)
    write_gz(bg, out / 'cancer_change_protein_background.tsv.gz')
    if args.depmap_models:
        models = pd.read_csv(args.depmap_models, dtype=str, keep_default_na=False).set_index('ModelID')
        used = sorted({x[4] for x in dplaced})
        mt = pd.DataFrame({'model_id': used, 'cell_line': [models['CellLineName'].get(m, '') for m in used],
                           'lineage': [models['OncotreeLineage'].get(m, '') for m in used],
                           'primary_disease': [models['OncotreePrimaryDisease'].get(m, '') for m in used]})
        mt.to_csv(out / 'depmap_24q4_change_models.tsv', sep='\t', index=False, lineterminator='\n')
    summary = {
        'generated': dt.date.today().isoformat(),
        'inputs': {'maf': {'file': Path(args.maf).name, 'url': s44b.MAF_URL, 'md5': md5_file(args.maf)},
                   'cases': {'file': Path(args.cases).name, 'url': s44b.CASES_URL, 'md5': md5_file(args.cases)},
                   'depmap_mutations': {'file': Path(args.depmap_mutations).name, 'release': 'DepMap Public 24Q4',
                                        'doi': '10.25452/figshare.plus.27993248', 'license': 'CC BY 4.0', 'md5': md5_file(args.depmap_mutations)},
                   'uniprot': {'release': args.uniprot_release,
                               'ensembl': {'file': Path(args.uniprot_ensembl).name, 'url': s44b.UNIPROT_ENSEMBL_URL, 'md5': md5_file(args.uniprot_ensembl)},
                               'isoforms': {'file': Path(args.uniprot_isoforms).name, 'url': s44b.UNIPROT_ISOFORMS_URL, 'md5': md5_file(args.uniprot_isoforms)}},
                   'crispr_gene_effect_header': {'file': Path(args.crispr_gene_effect).name, 'genes': len(screened),
                                                 'md5_of_file_read': md5_file(args.crispr_gene_effect),
                                                 'note': 'the header line of CRISPRGeneEffect.csv (DepMap 24Q4, MD5 of the whole file 6edf7ade09b9b34199210b559d4745d3)'},
                   'depmap_models': ({'file': Path(args.depmap_models).name, 'md5': md5_file(args.depmap_models)} if args.depmap_models else None),
                   'reference': {'fasta': args.fasta, 'md5': md5_file(args.fasta)}},
        'rule': {'tcga': 'FILTER PASS; in TCGA-OV and TCGA-LAML also FILTER wga alone (Ellrott et al. 2018); a missense or stop change '
                         'on a transcript UniProt assigns to a reviewed entry\'s displayed isoform (all_effects); where the variant falls '
                         'in no such transcript, a Missense_Mutation or Nonsense_Mutation of the MAF transcript (HGVSp_Short), placed on '
                         'the entry SWISSPROT names, else on the entries of the gene symbol; kept where the 2026_03 canonical sequence '
                         'has the changed residue at the position; one row per (entry, position, new residue, tumor sample)',
                 'depmap': 'a missense or stop call on the canonical isoform (UniprotID unlabelled or -1, else VepSwissprot) whose '
                           'residue the 2026_03 canonical sequence has at the position; one row per (entry, position, new residue, line)',
                 'background': 'every change kept by these rules, at any residue, counted per entry and residue class'},
        'denominators': {'tcga_tumor_samples': sum(len(v) for k, v in samples.items() if k),
                         'tcga_projects': len([k for k in samples if k]), 'depmap_lines_with_mutation_calls': int(lines)},
        'residues': written,
        'background_md5': md5_file(out / 'cancer_change_protein_background.tsv.gz'),
        'counts': dict(sorted(counts.items())),
    }
    (out / 'cancer_changes_summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    if args.manifest:
        write_manifest(Path(args.manifest), out, args, residues, args.downloaded)
    print(json.dumps({'denominators': summary['denominators'], 'residues': {k: {x: v[x] for x in v if x != 'md5'} for k, v in written.items()}}, indent=1))


if __name__ == '__main__':
    main()

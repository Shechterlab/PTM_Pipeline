#!/usr/bin/env python3
"""The annotation layers (inputs/layers/): data from other resources that the release joins to its sites and proteins.

Each layer is registered once, in LAYERS, with its staged files and what the build does with them:

  files            its tables in inputs/layers/, each listed in inputs/MANIFEST.tsv
  site_columns     columns it adds to the site table (release_codes.SITE_COLUMNS takes them in), and site_values,
                   the function that fills them
  protein_columns  likewise for the protein table (protein_values)
  analyses         the analyses it adds to the statistics table, and stats, the function that makes their rows
  website          the function that writes its website files from a release

    python scripts/release_layers.py check --inputs inputs
        every registered file is present and matches MANIFEST.tsv (hpc/build_release.sh runs it first)
    python scripts/release_layers.py website --release releases/v1.4.0 --inputs inputs --site-out docs/data
        each layer's website files, from the release (hpc/build_release.sh runs it last)

docs/ADDING_A_SOURCE.md says how to add a layer. The layers registered before this file existed (CPLM 4.0, the CARM1
knockout, the Human Protein Atlas) are listed with their files only; their columns and tests are made where they were
before (build_release.carm1_table, release_tables.hpa_rows, scripts/36).
"""
from __future__ import annotations

import argparse
import functools
import hashlib
import json
import math
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import pandas as pd

HERE = Path(__file__).resolve().parent


@dataclass
class Layer:
    key: str
    title: str
    files: list[str]                       # paths under inputs/layers/
    staging: list[str]                     # the scripts that stage them, in order
    licence: str
    citation: str
    site_columns: list[tuple] = field(default_factory=list)     # (name, type, codes, description), as release_codes
    protein_columns: list[tuple] = field(default_factory=list)
    codes: dict[str, dict] = field(default_factory=dict)        # the code lists its columns use
    analyses: dict[str, str] = field(default_factory=dict)      # analysis code -> description (statistics table)
    site_values: Callable | None = None     # (sites, inputs) -> DataFrame of site_columns on the index of sites
    protein_values: Callable | None = None  # (proteins, inputs) -> DataFrame of protein_columns on the index of proteins
    stats: Callable | None = None           # (sites with a boolean column counted, inputs, sequences, workdir) -> rows
    website: Callable | None = None         # (release folder, version, inputs, site out) -> None
    topic: str = ''                         # header fill of its site columns in the workbook (release_workbook.TOPICS)


def read_layer(inputs, name: str) -> pd.DataFrame:
    return pd.read_csv(Path(inputs) / 'layers' / name, sep='\t', dtype=str, keep_default_na=False, low_memory=False)


def join(values) -> str:
    return '; '.join(str(v) for v in values)


# ---- Cancer cell lines and tumours (DepMap 24Q4, PRISM 24Q2, TCGA PanCancer Atlas) -----------------------------------

DEPENDENT = 0.5        # DepMap's dependency probability at which a line depends on a gene
MOST_LINES = 0.5       # dependency classes: the share of screened lines that depend on the gene
SOME_LINES = 0.05
CHANGES = (  # source, staged table, what counts once, site columns (changes, count per change)
    ('tcga', 'tcga_pancan_arginine_changes.tsv', 'sample_id', 'tcga_changes', 'tcga_change_tumours'),
    ('depmap', 'depmap_24q4_arginine_changes.tsv', 'model_id', 'depmap_changes', 'depmap_change_lines'),
)
CHANGE_CATEGORY = {'tcga': 'tcga_tumours', 'depmap': 'depmap_cell_lines'}


@functools.lru_cache(maxsize=4)
def arginine_changes(inputs: str) -> dict[str, dict[tuple[str, int], dict[str, int]]]:
    """Per source, {(accession, position): {residue: tumour samples or cell lines with that change}}."""
    out = {}
    for src, name, who, *_ in CHANGES:
        t = read_layer(inputs, name)
        n = t.groupby(['canonical_UniProtAC', 'position', 'alt'])[who].nunique()
        per = defaultdict(dict)
        for (acc, pos, alt), k in n.items():
            per[(acc, int(pos))][alt] = int(k)
        out[src] = dict(per)
    return out


@functools.lru_cache(maxsize=4)
def gene_effect(inputs: str) -> pd.DataFrame:
    g = read_layer(inputs, 'depmap_24q4_gene_effect_summary.tsv').set_index('canonical_UniProtAC')
    share = g['dependent_lines'].astype(int) / g['lines_with_probability'].astype(int).where(lambda x: x > 0)
    ce = g['common_essential'] == 'True'
    g['dependency_class'] = ['common_essential' if c else 'most_lines' if s >= MOST_LINES else 'some_lines' if s >= SOME_LINES else 'few_or_none'
                             for c, s in zip(ce, share.fillna(0))]
    return g


def change_text(pos: int, alt: str) -> str:
    return f'R{pos}{alt}'


def cancer_site_values(sites: pd.DataFrame, inputs) -> pd.DataFrame:
    changes = arginine_changes(str(inputs))
    out = pd.DataFrame(index=sites.index)
    keys = [(a, int(p)) for a, p in zip(sites['uniprot_accession'], sites['position'])]
    for src, _name, _who, col, count_col in CHANGES:
        per = [sorted(changes[src].get(k, {}).items(), key=lambda x: (-x[1], x[0])) for k in keys]
        out[col] = [join(change_text(k[1], alt) for alt, _ in x) for k, x in zip(keys, per)]
        out[count_col] = [join(n for _, n in x) for x in per]
    return out


def cancer_protein_values(proteins: pd.DataFrame, inputs) -> pd.DataFrame:
    g = gene_effect(str(inputs))
    entry = proteins['uniprot_accession'].str.split('-').str[0]   # an isoform takes its entry's gene
    pick = lambda col, f=str: [f(g.at[a, col]) if a in g.index else None for a in entry]  # noqa: E731
    return pd.DataFrame({'depmap_gene_effect_median': pick('median', float), 'depmap_lines_screened': pick('lines', int),
                         'depmap_dependent_lines': pick('dependent_lines', int), 'depmap_dependency_class': pick('dependency_class')},
                        index=proteins.index)


def change_tables(sites: pd.DataFrame, inputs, sequences: dict) -> dict[str, list[tuple[str, tuple[int, int, int, int]]]]:
    """Per source, one 2x2 table per methylated protein whose gene DepMap screened: methylarginines changed (in at least one
    tumour or line), methylarginines not changed, other arginines changed, other arginines not changed. The counted sites
    are the methylarginines; the sites of the table not counted (kept pending, or lower-confidence sites the analyses
    leave out) are in neither group."""
    changes = arginine_changes(str(inputs))
    screened = set(gene_effect(str(inputs)).index)
    pos = sites.assign(position=sites['position'].astype(int))
    counted = pos[pos['counted']].groupby('uniprot_accession')['position'].agg(set).to_dict()
    left_out = pos[~pos['counted']].groupby('uniprot_accession')['position'].agg(set).to_dict()
    out = {}
    for src, *_ in CHANGES:
        hit = changes[src]
        tables = []
        for acc in sorted(a for a in counted if a in screened and a in sequences):
            meth, skip = counted[acc], left_out.get(acc, set())
            rs = [k + 1 for k, ch in enumerate(sequences[acc]) if ch == 'R']
            mm = [r for r in rs if r in meth]
            oo = [r for r in rs if r not in meth and r not in skip]
            a, c = sum((acc, r) in hit for r in mm), sum((acc, r) in hit for r in oo)
            tables.append((acc, (a, len(mm) - a, c, len(oo) - c)))
        out[src] = tables
    return out


def bh(ps: list[float]) -> list[float]:
    """Benjamini-Hochberg, as stats::p.adjust."""
    order = sorted(range(len(ps)), key=lambda i: ps[i])
    q, prev = [0.0] * len(ps), 1.0
    for rank in range(len(ps), 0, -1):
        i = order[rank - 1]
        prev = min(prev, ps[i] * len(ps) / rank)
        q[i] = prev
    return q


def cancer_stats(sites: pd.DataFrame, inputs, sequences: dict, workdir: Path) -> list[dict]:
    from scipy.stats import fisher_exact
    from common import mantel_haenszel_ci
    from enrichment import Job, run
    rows = []
    # Dependency classes: methylated proteins (by their entry's gene) against the other screened entries.
    g = gene_effect(str(inputs))
    universe = set(g.index)
    target = {a.split('-')[0] for a in sites.loc[sites['counted'], 'uniprot_accession']} & universe
    sets = pd.DataFrame({'term': g['dependency_class'], 'gene': g.index})
    res = run([Job('depmap_dependency', target, universe, sets=sets, whole_universe=True)], Path(workdir) / 'enricher')['depmap_dependency']
    order = list(CANCER.codes['dependency_class'])
    for r in sorted(res.to_dict('records'), key=lambda r: order.index(r['ID'])):
        rows.append({'analysis': 'depmap_dependency', 'category': r['ID'], 'stratum': 'all', 'a': r['a'], 'b': r['b'], 'c': r['c'], 'd': r['d'],
                     'fold_enrichment': r['FoldEnrichment'], 'log2_fold_enrichment': math.log2(r['FoldEnrichment']),
                     'p': r['pvalue'], 'q': r['p.adjust'], 'method': 'enricher',
                     'background': 'other reviewed human proteins whose gene DepMap 24Q4 screened (CRISPR, Chronos)', 'n_strata': ''})
    # Changes at methylarginines against the other arginines of the same proteins: pooled (Fisher) and within protein
    # (Mantel-Haenszel, a protein a stratum).
    pooled, within = [], []
    for src, tables in change_tables(sites, inputs, sequences).items():
        a, b, c, d = (sum(t[k] for _, t in tables) for k in range(4))
        _, p = fisher_exact([[a, b], [c, d]])
        odds = (a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5))
        pooled.append({'analysis': 'cancer_change', 'category': CHANGE_CATEGORY[src], 'stratum': 'all', 'a': a, 'b': b, 'c': c, 'd': d,
                       'odds_ratio': odds, 'log2_odds_ratio': math.log2(odds), 'p': float(p), 'method': 'fisher',
                       'background': 'other arginines of the methylated proteins whose gene DepMap screened', 'n_strata': ''})
        mh, lo, hi, mp, used = mantel_haenszel_ci([t for _, t in tables])
        within.append({'analysis': 'cancer_change', 'category': CHANGE_CATEGORY[src], 'stratum': 'all', 'a': a, 'b': b, 'c': c, 'd': d,
                       'odds_ratio': mh, 'log2_odds_ratio': math.log2(mh) if mh and mh > 0 and math.isfinite(mh) else None,
                       'odds_ratio_ci_low': lo, 'odds_ratio_ci_high': hi, 'p': mp, 'method': 'mantel_haenszel',
                       'background': 'other arginines of the same protein, each protein a stratum', 'n_strata': used})
    for group in (pooled, within):
        for r, q in zip(group, bh([r['p'] for r in group])):
            r['q'] = q
    return rows + pooled + within


def cancer_website(release: Path, version: str, inputs: Path, site_out: Path) -> None:
    subprocess.run([sys.executable, str(HERE / '44c_export_depmap.py'), '--layers', str(inputs / 'layers'), '--out', str(site_out / 'depmap'),
                    '--sites', str(release / f'human_methylarginine_sites_{version}.tsv'),
                    '--statistics', str(release / f'statistics_{version}.tsv'),
                    '--fasta', str(inputs / 'reference' / 'uniprot_human_reviewed_2026_03.fasta')], check=True)


CANCER = Layer(
    key='cancer', title='Cancer cell lines and tumours',
    files=['depmap_24q4_gene_effect_summary.tsv', 'depmap_24q4_arginine_changes.tsv', 'depmap_24q4_models.tsv',
           'depmap_24q4_prmt_gene_effect.tsv', 'prism_24q2_prmt_inhibitors.tsv', 'tcga_pancan_arginine_changes.tsv',
           'tcga_pancan_studies.tsv', 'depmap_summary.json'],
    staging=['scripts/44b_fetch_tcga_arginine_changes.py', 'scripts/44_stage_depmap.py'],
    licence='DepMap Public 24Q4 and PRISM Repurposing 24Q2: CC BY 4.0; cBioPortal (TCGA PanCancer Atlas): ODC Open Database License',
    citation='DepMap Public 24Q4 (doi:10.25452/figshare.plus.27993248; Tsherniak et al. 2017, doi:10.1016/j.cell.2017.06.010; Chronos, '
             'Dempster et al. 2021, doi:10.1186/s13059-021-02540-7); PRISM Repurposing 24Q2 (doi:10.6084/m9.figshare.25917643; Corsello et al. '
             '2020, doi:10.1038/s43018-019-0018-6); TCGA PanCancer Atlas (Hoadley et al. 2018, doi:10.1016/j.cell.2018.03.022; MC3 calls, '
             'Ellrott et al. 2018, doi:10.1016/j.cels.2018.03.002) through cBioPortal (Cerami et al. 2012, doi:10.1158/2159-8290.CD-12-0095; '
             'Gao et al. 2013, doi:10.1126/scisignal.2004088; de Bruijn et al. 2023, doi:10.1158/0008-5472.CAN-23-0816)',
    site_columns=[
        ('tcga_changes', 'list<string>', '', 'missense and stop changes of the arginine in TCGA PanCancer Atlas tumours (MC3 calls through cBioPortal, '
         '32 studies, 10,443 sequenced samples), most tumours first, * a stop. A change counts where the UniProt canonical sequence has this '
         'arginine at the position cBioPortal gives (its RefSeq transcript), so a change numbered on a transcript that differs is not here '
         '(inputs/layers/README.md)'),
        ('tcga_change_tumours', 'list<integer>', '', 'for each of these changes, the tumour samples with it, same order'),
        ('depmap_changes', 'list<string>', '', 'missense and stop changes of the arginine in DepMap 24Q4 cell lines (OmicsSomaticMutations, on the '
         'canonical UniProt isoform), most lines first'),
        ('depmap_change_lines', 'list<integer>', '', 'for each of these changes, the cell lines with it, same order'),
    ],
    protein_columns=[
        ('depmap_gene_effect_median', 'number', '', 'median CRISPR gene effect (Chronos) of the protein\'s gene across the DepMap 24Q4 lines that '
         'screened it: 0 no effect, -1 the median of commonly essential genes; empty when DepMap did not screen the gene (an isoform '
         'takes its entry\'s gene)'),
        ('depmap_lines_screened', 'integer', '', 'lines with a gene effect for the gene'),
        ('depmap_dependent_lines', 'integer', '', f'lines whose dependency probability for the gene is at least {DEPENDENT}'),
        ('depmap_dependency_class', 'code', 'dependency_class', 'the gene\'s dependency class in DepMap 24Q4, empty when not screened'),
    ],
    codes={'dependency_class': {
        'common_essential': 'DepMap 24Q4 calls the gene commonly essential (CRISPRInferredCommonEssentials: it ranks among the strongest '
                            'dependencies of at least 90% of lines)',
        'most_lines': f'not commonly essential, and at least {MOST_LINES:.0%} of the screened lines depend on it',
        'some_lines': f'not commonly essential, and {SOME_LINES:.0%} to under {MOST_LINES:.0%} of the lines depend on it (selective dependencies)',
        'few_or_none': f'fewer than {SOME_LINES:.0%} of the lines depend on it',
    }},
    analyses={
        'depmap_dependency': 'methylated proteins (by the entry\'s gene) in a DepMap 24Q4 dependency class (code dependency_class), against the '
                             'other reviewed human proteins whose gene DepMap screened. The classes partition the screened entries. Essential '
                             'proteins are abundant, and abundant proteins are the ones mass spectrometry finds most, so methylated proteins '
                             'are expected to be enriched among commonly essential genes whatever methylation does',
        'cancer_change': 'methylarginines with a missense or stop change in at least one TCGA PanCancer Atlas tumour (category tcga_tumours) or '
                         'DepMap 24Q4 cell line (depmap_cell_lines), against the other arginines of the methylated proteins whose gene DepMap '
                         'screened: pooled (method fisher) and within each protein (method mantel_haenszel, each protein a stratum). Sites of '
                         'the table that the analyses leave out are in neither group. Arginine codons CGN hold a CpG, the dinucleotide that '
                         'mutates most often in tumours, so codon use as well as selection shapes the comparison',
    },
    site_values=cancer_site_values, protein_values=cancer_protein_values, stats=cancer_stats, website=cancer_website,
    topic='9F1239',
)

LAYERS = [
    Layer('cplm4', 'CPLM 4.0 lysine methylation', ['cplm4_summary.json'], ['scripts/39_stage_cplm4_lys_methylation.py'],
          'All rights reserved (© The CUCKOO Workgroup; no licence stated): the tables stay out of Git, checksum_only in MANIFEST.tsv',
          'Zhang et al., Nucleic Acids Res 50:D451 (2022)'),
    Layer('carm1_ko', 'CARM1 knockout response (Shishkova et al. 2017)', ['carm1_ko_peptides.tsv', 'carm1_ko_sites.tsv', 'carm1_ko_summary.json'],
          ['scripts/41_stage_carm1_ko_sites.py'], 'CC BY 4.0', 'Shishkova et al., Nat Commun 8:15571 (2017)'),
    Layer('hpa', 'Human Protein Atlas 25.1 subcellular locations', ['hpa_subcellular_locations.tsv', 'hpa_summary.json'],
          ['scripts/42_stage_hpa_subcellular.py'], 'CC BY 4.0', 'Thul et al., Science 356:eaal3321 (2017)'),
    CANCER,
]


def site_columns() -> list[tuple]:
    return [c for layer in LAYERS for c in layer.site_columns]


def protein_columns() -> list[tuple]:
    return [c for layer in LAYERS for c in layer.protein_columns]


def add_site_values(sites: pd.DataFrame, inputs) -> pd.DataFrame:
    parts = [layer.site_values(sites, inputs) for layer in LAYERS if layer.site_values]
    return pd.concat([sites, *parts], axis=1) if parts else sites


def add_protein_values(proteins: pd.DataFrame, inputs) -> pd.DataFrame:
    parts = [layer.protein_values(proteins, inputs) for layer in LAYERS if layer.protein_values]
    out = pd.concat([proteins, *parts], axis=1) if parts else proteins
    out.attrs = proteins.attrs
    return out


def statistics(sites: pd.DataFrame, inputs, sequences: dict, workdir: Path) -> list[dict]:
    """sites: uniprot_accession, position and counted (the analyses count the site) of every site of the table."""
    return [r for layer in LAYERS if layer.stats for r in layer.stats(sites, inputs, sequences, Path(workdir) / 'layers' / layer.key)]


# ---- Command line ------------------------------------------------------------------------------------------------------

def check(inputs: Path) -> list[str]:
    manifest = pd.read_csv(inputs / 'MANIFEST.tsv', sep='\t', dtype=str, keep_default_na=False).set_index('file')
    bad = []
    for layer in LAYERS:
        for name in layer.files:
            rel, path = f'layers/{name}', inputs / 'layers' / name
            if rel not in manifest.index:
                bad.append(f'{layer.key}: {rel} is not in MANIFEST.tsv')
            elif not path.exists():
                bad.append(f'{layer.key}: {rel} is missing')
            elif hashlib.md5(path.read_bytes()).hexdigest() != manifest.at[rel, 'md5']:
                bad.append(f'{layer.key}: {rel} differs from MANIFEST.tsv (stage it again with {", ".join(layer.staging)})')
    return bad


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    c = sub.add_parser('check')
    c.add_argument('--inputs', default='inputs')
    w = sub.add_parser('website')
    w.add_argument('--release', required=True)
    w.add_argument('--inputs', default='inputs')
    w.add_argument('--site-out', default='docs/data')
    args = ap.parse_args()
    inputs = Path(args.inputs)
    if args.cmd == 'check':
        bad = check(inputs)
        if bad:
            raise SystemExit('the layer inputs do not match MANIFEST.tsv:\n  ' + '\n  '.join(bad))
        print(f'{sum(len(layer.files) for layer in LAYERS)} layer files of {len(LAYERS)} layers match MANIFEST.tsv')
    else:
        release = Path(args.release)
        version = json.loads((release / 'release.json').read_text())['version']
        for layer in LAYERS:
            if layer.website:
                layer.website(release, version, inputs, Path(args.site_out))
                print(f'{layer.key}: website files written')


if __name__ == '__main__':
    main()

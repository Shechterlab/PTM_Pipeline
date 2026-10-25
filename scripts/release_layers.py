#!/usr/bin/env python3
"""The annotation layers (inputs/layers/): data from other resources that the release joins to its sites and proteins.

Each layer is registered once, in LAYERS, with its staged files and what the build does with them:

  files            its tables in inputs/layers/, each listed in inputs/MANIFEST.tsv
  site_columns     columns it adds to the site table (release_codes.SITE_COLUMNS takes them in), and site_values,
                   the function that fills them
  protein_columns  likewise for the protein table (protein_values)
  analyses         the analyses it adds to the statistics table, and stats, the function that makes their rows
  supporting       tables it adds to the release's supporting/ folder, with their columns (the data dictionary and the
                   README take them in)
  website          the function that writes its website files from a release

    python scripts/release_layers.py check --inputs inputs
        every registered file is present and matches MANIFEST.tsv (hpc/build_release.sh runs it first)
    python scripts/release_layers.py website --release releases/v1.4.0 --inputs inputs --site-out docs/data
        each layer's website files, from the release (hpc/build_release.sh runs it last)

ADDING_A_SOURCE.md says how to add a layer. The layers registered before this file existed (CPLM 4.0, the CARM1
knockout, the Human Protein Atlas) are listed with their files only; their columns and tests are made where they were
before (build_release.carm1_table, release_tables.hpa_rows, scripts/36). The writer perturbations (knockouts, knockdowns
and inhibitors) are one layer, each dataset one entry in PERTURBATIONS.
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
    files: list[str]                       # paths under inputs/<folder>/
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
    elsewhere: str = ''                     # what it adds that is made outside this file (layers registered before it)
    supporting: Callable | None = None      # (sites, inputs, version) -> {key: (path in the release, frame, columns)}
    supporting_rows: dict[str, str] = field(default_factory=dict)  # key -> what one row of that table is (README)
    folder: str = 'layers'                  # where its files are, under inputs/
    since: str = ''                         # the release that added it (CHANGES.md describes it there)
    summary: str = ''                       # what it adds, for CHANGES.md
    changes: Callable | None = None         # (files of the release) -> sentences with its counts, for CHANGES.md


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


def cancer_changes(files: dict) -> list[str]:
    sites = files['sites'][1]
    t, d = (sites['tcga_changes'] != ''), (sites['depmap_changes'] != '')
    return [f'{int(t.sum()):,} sites have a change in at least one tumour and {int(d.sum()):,} in at least one cell line.']


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
    topic='9F1239', since='v1.4.0',
    summary='Cancer genetics beside the sites, not evidence about methylation: how much DepMap\'s CRISPR-screened cancer cell lines '
            'depend on each methylated protein (protein columns `depmap_*`), and the missense and stop changes of each methylarginine '
            'in TCGA PanCancer Atlas tumours and DepMap cell lines (site columns `tcga_changes`, `depmap_changes`). Two tests: '
            'dependency classes of the methylated proteins (analysis `depmap_dependency`) and changes at methylarginines against the '
            'other arginines of their proteins (analysis `cancer_change`).',
    changes=lambda files: cancer_changes(files),
)

def release_table(release: Path, name: str) -> pd.DataFrame:
    return pd.read_csv(Path(release) / name, sep='\t', dtype=str, keep_default_na=False)


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, separators=(',', ':'), ensure_ascii=False), encoding='utf-8')


# ---- Writer perturbations: knockouts, knockdowns and inhibitors -------------------------------------------------------

RELEASE_Q, RELEASE_LOG2 = 0.05, 1.0   # the release's rule for every perturbation: |log2| >= 1 at q < 0.05 (as the CARM1 knockout)
SOURCE_CALL = {'down': 'decreased', 'up': 'increased', 'unchanged': 'no_significant_change', 'mixed': 'mixed',
               'decreased': 'decreased', 'increased': 'increased', 'no_significant_change': 'no_significant_change'}


@dataclass
class Perturbation:
    key: str             # its code (CODES['perturbation'])
    label: str           # short name, as the website shows it
    description: str
    enzyme: str
    kind: str            # knockout, knockdown or inhibitor
    cell_line: str
    study: str
    doi: str
    peptides: str        # its staged peptide table in inputs/layers/ (canonical_UniProtAC, methyl_positions)
    log2: str            # the column of log2 treated / control
    p: str               # the column of the peptide p value, '' when the study gives none
    q: str               # the column of an adjusted p value to use as given, or 'bh': p adjusted here
    source_call: str     # the column of each peptide's call by the study's own rule
    source_rule: str
    select: dict = field(default_factory=dict)       # column -> value: the rows of this arm in a shared table
    bh_groups: list = field(default_factory=list)    # columns within which p is adjusted ('bh')


PERTURBATIONS = [
    Perturbation('carm1_ko_mcf7', 'CARM1 KO, MCF7', 'CARM1 knockout against parental MCF7 cells, ADMA immunoprecipitation, label-free '
                 '(Shishkova et al. 2017)', 'CARM1', 'knockout', 'MCF7', 'Shishkova et al., Nat Commun 8:15571 (2017)', '10.1038/ncomms15571',
                 'carm1_ko_peptides.tsv', 'log2_ko_wt', 'p', 'bh', 'call', 'twofold at nominal p < 0.05 (the study\'s rule)',
                 select={'cell_line': 'MCF7'}),
    Perturbation('carm1_ko_mdamb231', 'CARM1 KO, MDA-MB-231', 'CARM1 knockout against parental MDA-MB-231 cells, ADMA '
                 'immunoprecipitation, label-free (Shishkova et al. 2017)', 'CARM1', 'knockout', 'MDA-MB-231',
                 'Shishkova et al., Nat Commun 8:15571 (2017)', '10.1038/ncomms15571', 'carm1_ko_peptides.tsv', 'log2_ko_wt', 'p', 'bh', 'call',
                 'twofold at nominal p < 0.05 (the study\'s rule)', select={'cell_line': 'MDA-MB-231'}),
    Perturbation('prmt7_ko_hct116', 'PRMT7 KO, HCT116', 'PRMT7 knockout against parental HCT116 cells, monomethylarginine '
                 'immunoprecipitation, SILAC, 4 replicates (Szewczyk et al. 2020)', 'PRMT7', 'knockout', 'HCT116',
                 'Szewczyk et al., Nat Commun 11:2396 (2020)', '10.1038/s41467-020-16271-z', 'prmt7_ko_peptides.tsv', 'log2_ko_wt', 'p', 'q', 'call',
                 'twofold at BH-adjusted p < 0.01 (the study\'s rule; its adjusted p values are used as given)'),
    Perturbation('prmt1_kd_hek293t', 'PRMT1 KD, HEK293T', 'PRMT1 shRNA knockdown against a non-targeting shRNA in HEK293T cells, '
                 'label-free, 2 against 2 biological replicates, SCX and antibody enrichments (Hartel et al. 2019)', 'PRMT1', 'knockdown',
                 'HEK293T', 'Hartel et al., Mol Cell Proteomics 18:2149-2164 (2019)', '10.1074/mcp.RA119.001625',
                 'hartel2019_prmt1_kd_peptides.tsv', 'log2_kd_control', 'p', 'bh', 'call',
                 'Perseus t-test with permutation-based FDR q < 0.05 and S0 = 0.5 (the study\'s rule; 2 against 2 samples allow few '
                 'distinct permutations)', bh_groups=['table']),
    Perturbation('prmt5i_epz015666_hela', 'PRMT5i EPZ015666, HeLa', 'PRMT5 inhibitor EPZ015666 against untreated HeLa cells, SDMA '
                 'antibody after ulilysin digestion, SILAC (Lu et al. 2024); the treated channel is read from the data, as the '
                 'readable text does not state it', 'PRMT5', 'inhibitor', 'HeLa', 'Lu et al., J Proteome Res 23:1014-1027 (2024)',
                 '10.1021/acs.jproteome.3c00724', 'lu2024_prmt5_peptides.tsv', 'log2_inhibitor_control', 'p', 'bh', 'call',
                 'twofold at nominal p < 0.05, inferred: it reproduces 254 of the 255 sites of the study\'s Table 1.4 that map to the '
                 'reference and adds none (the study states no threshold in the text that could be read)'),
]


@functools.lru_cache(maxsize=4)
def perturbation_rows(inputs: str) -> pd.DataFrame:
    """One row per site and perturbation: the release's call (|log2| >= 1 at q < 0.05, by scripts/peptide_calls.py over the
    peptides carrying the site, q adjusted here where the study gives p only) and the call by the study's own rule."""
    import peptide_calls as pc
    rows = []
    for d in PERTURBATIONS:
        pep = read_layer(inputs, d.peptides)
        for col, value in d.select.items():
            pep = pep[pep[col] == value]
        pep = pep.copy()
        pep['rel_log2'] = pd.to_numeric(pep[d.log2], errors='coerce')
        pep['rel_p'] = pd.to_numeric(pep[d.p], errors='coerce') if d.p else float('nan')
        pep = pep[pep['rel_log2'].notna()]
        if d.q == 'bh':   # over every quantified peptide of the comparison, placed on the reference or not (as v1.3.1's CARM1 table)
            pep['rel_q'] = float('nan')
            tested = pep[pep['rel_p'].notna()]
            for _, g in (tested.groupby(d.bh_groups) if d.bh_groups else [((), tested)]):
                pep.loc[g.index, 'rel_q'] = pc.bh(g['rel_p'].tolist())
        else:
            pep['rel_q'] = pd.to_numeric(pep[d.q], errors='coerce')
        pep = pep[(pep['canonical_UniProtAC'] != '') & (pep['methyl_positions'] != '')]
        groups = defaultdict(list)
        for r in pep.itertuples(index=False):
            positions = [int(float(x)) for x in str(r.methyl_positions).split(';') if x]
            for pos in positions:
                groups[(r.canonical_UniProtAC, pos)].append((r.rel_q, r.rel_log2, len(positions), r.rel_p, SOURCE_CALL.get(getattr(r, d.source_call), '')))
        for (acc, pos), vals in sorted(groups.items()):
            res = pc.site_response([(q, l2, n) for q, l2, n, _, _ in vals], RELEASE_Q, RELEASE_LOG2, strict=True)
            tested = any(not pd.isna(q) for q, *_ in vals)
            ps = [p for *_, p, _ in vals if not pd.isna(p)]
            rows.append({'site_id': f'{acc}:R{pos}', 'perturbation': d.key, 'n_peptides': res['n_peptides'],
                         'log2_treated_over_control': round(res['log2_median'], 3), 'p': min(ps) if ps else None,
                         'q': res['adj_p_min'] if tested else None, 'n_decreased': res['n_decreased'], 'n_increased': res['n_increased'],
                         'n_multi_site_peptides': res['n_multi_site_peptides'], 'call': res['call'] if tested else 'not_tested',
                         'call_single_site_peptides': (res['call_single_site_peptides'] or '') if tested else '',
                         'call_source_rule': pc.combine_calls(c for *_, c in vals)})
    return pd.DataFrame(rows)


def perturbation_site_values(sites: pd.DataFrame, inputs) -> pd.DataFrame:
    t = perturbation_rows(str(inputs))
    order = {d.key: k for k, d in enumerate(PERTURBATIONS)}
    by = defaultdict(list)
    for r in t.itertuples(index=False):
        by[r.site_id].append((order[r.perturbation], r.perturbation, r.call, r.call_source_rule))
    ids = sites['uniprot_accession'] + ':R' + sites['position'].astype(int).astype(str)
    vals = [sorted(by.get(s, [])) for s in ids]
    return pd.DataFrame({'perturbations': [join(x[1] for x in v) for v in vals],
                         'perturbation_calls': [join(x[2] for x in v) for v in vals],
                         'perturbation_calls_source_rule': [join(x[3] or 'not_tested' for x in v) for v in vals]}, index=sites.index)


PERTURBATION_COLUMNS = [
    ('site_id', 'string', '', 'the site'),
    ('perturbation', 'code', 'perturbation', 'the knockout, knockdown or inhibitor, and the cells'),
    ('n_peptides', 'integer', '', 'quantified peptides carrying the site'),
    ('log2_treated_over_control', 'number', '', 'median log2 perturbed / control over those peptides: descriptive, not the effect of a test'),
    ('p', 'number', '', 'the lowest p value over the peptides, as the study gives it'),
    ('q', 'number', '', f'the lowest q over the peptides: the study\'s adjusted p where it gives one, else its p adjusted here (Benjamini-Hochberg, within the comparison); empty when the study gives no p'),
    ('n_decreased', 'integer', '', f'peptides at log2 <= -{RELEASE_LOG2:g} and q < {RELEASE_Q:g}'),
    ('n_increased', 'integer', '', f'peptides at log2 >= {RELEASE_LOG2:g} and q < {RELEASE_Q:g}'),
    ('n_multi_site_peptides', 'integer', '', 'of the peptides, those carrying another methylarginine too'),
    ('call', 'code', 'perturbation_call', f'the release\'s rule for every perturbation: a peptide changes at |log2| >= {RELEASE_LOG2:g} and q < {RELEASE_Q:g}, combined over the peptides carrying the site (scripts/peptide_calls.py): peptide-level evidence'),
    ('call_single_site_peptides', 'code', 'perturbation_call', 'the same over the peptides carrying no other methylarginine; empty when every peptide carries another'),
    ('call_source_rule', 'code', 'perturbation_call', 'the call by the study\'s own rule (code perturbation), for comparison'),
]


def perturbation_supporting(sites: pd.DataFrame, inputs, version: str) -> dict:
    t = perturbation_rows(str(inputs))
    ids = set(sites['uniprot_accession'] + ':R' + sites['position'].astype(int).astype(str))
    return {'writer_perturbations': (f'supporting/writer_perturbations_{version}.tsv', t[t['site_id'].isin(ids)].reset_index(drop=True),
                                     PERTURBATION_COLUMNS)}


def perturbation_changes(files: dict) -> list[str]:
    t = files['writer_perturbations'][1]
    out = []
    for d in PERTURBATIONS:
        x = t[t['perturbation'] == d.key]
        if len(x):
            dec, inc, differ = int((x['call'] == 'decreased').sum()), int((x['call'] == 'increased').sum()), int((x['call'] != x['call_source_rule']).sum())
            out.append(f'{d.label} ({d.study}): {len(x):,} sites quantified, {dec:,} decreased and {inc:,} increased by the release\'s rule; '
                       f'the study\'s rule calls {differ:,} of them differently.')
    return out


def perturbation_website(release: Path, version: str, inputs: Path, site_out: Path) -> None:
    """SITE_OUT/perturbations.json: the datasets, and for each site of the release its rows of supporting/writer_perturbations."""
    t = release_table(release, f'supporting/writer_perturbations_{version}.tsv')
    keys = [d.key for d in PERTURBATIONS]
    calls = list(PERTURBATION.codes['perturbation_call'])
    num = lambda v, nd=3: None if v == '' else float(f'{float(v):.{nd}g}')  # noqa: E731
    sites = defaultdict(list)
    for r in t.itertuples(index=False):
        sites[r.site_id].append([keys.index(r.perturbation), calls.index(r.call), calls.index(r.call_source_rule) if r.call_source_rule else -1,
                                 num(r.log2_treated_over_control), num(r.q), int(r.n_peptides), int(r.n_decreased), int(r.n_increased),
                                 int(r.n_multi_site_peptides), calls.index(r.call_single_site_peptides) if r.call_single_site_peptides else -1])
    payload = {'fields': ['perturbation', 'call', 'call_source_rule', 'log2_treated_over_control', 'q', 'n_peptides', 'n_decreased',
                          'n_increased', 'n_multi_site_peptides', 'call_single_site_peptides'],
               'calls': calls, 'rule': f'|log2| >= {RELEASE_LOG2:g} at q < {RELEASE_Q:g}, combined over the peptides carrying the site',
               'perturbations': [[d.key, d.label, d.description, d.enzyme, d.kind, d.cell_line, d.study, d.doi, d.source_rule] for d in PERTURBATIONS],
               'version': version, 'sites': dict(sorted(sites.items()))}
    write_json(site_out / 'perturbations.json', payload)


PERTURBATION = Layer(
    key='writer_perturbations', title='Writer perturbations (PRMT knockouts, knockdowns and inhibitors)',
    files=['prmt7_ko_peptides.tsv', 'prmt7_ko_sites.tsv', 'prmt7_ko_searched_sequences.tsv', 'prmt7_ko_summary.json',
           'hartel2019_prmt1_kd_peptides.tsv', 'hartel2019_prmt1_kd_sites.tsv', 'hartel2019_prmt1_kd_summary.json',
           'lu2024_prmt5_peptides.tsv', 'lu2024_prmt5_sites.tsv', 'lu2024_prmt5_summary.json'],
    staging=['scripts/41_stage_carm1_ko_sites.py', 'scripts/45_stage_prmt7_szewczyk2020.py', 'scripts/46_stage_hartel2019_prmt1_kd.py',
             'scripts/51_stage_lu2024_prmt5.py'],
    licence='CC BY 4.0 (Shishkova et al. 2017; Szewczyk et al. 2020; Hartel et al. 2019, from 2020-09-15); CC BY-NC 4.0 (Lu et al. 2024 supporting information)',
    citation='; '.join(sorted({f'{d.study}, doi:{d.doi}' for d in PERTURBATIONS})),
    site_columns=[
        ('perturbations', 'list<code>', 'perturbation', 'the knockouts, knockdowns and inhibitors that quantified a peptide carrying the '
         'site, in code order (supporting/writer_perturbations)'),
        ('perturbation_calls', 'list<code>', 'perturbation_call', f'for each of them, the release\'s call (|log2| >= {RELEASE_LOG2:g} at q < '
         f'{RELEASE_Q:g}, the same rule for every perturbation): peptide-level evidence, same order'),
        ('perturbation_calls_source_rule', 'list<code>', 'perturbation_call', 'for each of them, the call by the study\'s own rule, same order'),
    ],
    codes={'perturbation': {d.key: f'{d.description}; the study\'s rule: {d.source_rule}' for d in PERTURBATIONS},
           'perturbation_call': {'decreased': 'a peptide carrying the site was lower with the perturbation: at least one decreased and none increased',
                                 'increased': 'a peptide carrying the site was higher with the perturbation: at least one increased and none decreased',
                                 'mixed': 'peptides changed in opposite directions',
                                 'no_significant_change': 'quantified, and no peptide changed by the rule: not evidence that the site is unchanged',
                                 'not_tested': 'quantified, but the study gives no p value, so the rule cannot be applied'}},
    site_values=perturbation_site_values, supporting=perturbation_supporting, website=perturbation_website, since='v1.4.0',
    summary='Knockouts, knockdowns and inhibitors of the PRMTs that quantified methylated peptides, called by one rule for all '
            f'(a peptide changes at |log2| >= {RELEASE_LOG2:g} and q < {RELEASE_Q:g}; a site is decreased when a peptide carrying it '
            'decreased and none increased: scripts/peptide_calls.py, as the CARM1 knockout of v1.3.1), with each study\'s own rule '
            'beside it (`call_source_rule`). Peptide-level evidence with no total-protein correction. Site columns `perturbations`, '
            '`perturbation_calls`, `perturbation_calls_source_rule`; one row per site and perturbation in `supporting/writer_perturbations`. '
            'The CARM1 knockout rows equal `supporting/carm1_ko`, which stays.',
    changes=lambda files: perturbation_changes(files),
    supporting_rows={'writer_perturbations': 'site and knockout, knockdown or inhibitor (PRMT1, PRMT5, PRMT7, CARM1): its peptides, '
                                             'change, and the call by the release\'s rule and by the study\'s'},
    topic='5B21B6',
)


# ---- Protein abundance (PaxDb 6.1) ------------------------------------------------------------------------------------

@functools.lru_cache(maxsize=4)
def paxdb_abundance(inputs: str) -> pd.DataFrame:
    return read_layer(inputs, 'paxdb_human_integrated_abundance.tsv').set_index('canonical_UniProtAC')


def paxdb_protein_values(proteins: pd.DataFrame, inputs) -> pd.DataFrame:
    p = paxdb_abundance(str(inputs))
    entry = proteins['uniprot_accession'].str.split('-').str[0]   # an isoform takes its entry's value
    pick = lambda col: [float(p.at[a, col]) if a in p.index else None for a in entry]  # noqa: E731
    return pd.DataFrame({'paxdb_abundance_ppm': pick('abundance_ppm'), 'paxdb_percentile': pick('percentile')}, index=proteins.index)


def paxdb_changes(files: dict) -> list[str]:
    p = files['proteins'][1]
    v = pd.to_numeric(p['paxdb_percentile'], errors='coerce')
    return [f'{int(v.notna().sum()):,} of {len(p):,} proteins have a value; their median percentile is {v.median():.0f}.']


def paxdb_website(release: Path, version: str, site_out: Path) -> None:
    """SITE_OUT/paxdb.json: abundance and percentile of each entry of the release's protein table."""
    p = release_table(release, f'human_methylarginine_proteins_{version}.tsv')
    p = p[p['paxdb_abundance_ppm'] != '']
    prot = {a.split('-')[0]: [float(f'{float(v):.4g}'), float(c)] for a, v, c in zip(p['uniprot_accession'], p['paxdb_abundance_ppm'], p['paxdb_percentile'])}
    write_json(site_out / 'paxdb.json', {'version': version, 'fields': ['abundance_ppm', 'percentile'], 'prot': dict(sorted(prot.items()))})


PAXDB = Layer(
    key='paxdb', title='PaxDb 6.1 protein abundance',
    files=['paxdb_human_integrated_abundance.tsv', 'paxdb_summary.json'],
    staging=['scripts/47_stage_paxdb.py'],
    licence='CC BY 4.0',
    citation='PaxDb 6.1, H. sapiens whole-organism integrated dataset (pax-db.org); Huang et al., Mol Cell Proteomics 22:100640 (2023), '
             'doi:10.1016/j.mcpro.2023.100640',
    protein_columns=[
        ('paxdb_abundance_ppm', 'number', '', 'protein abundance in parts per million of the proteome, PaxDb 6.1 H. sapiens whole-organism '
         'integrated dataset (a weighted average of 7 whole-organism datasets, PeptideAtlas 2011 weighing 0.75), through PaxDb\'s '
         'UniProt cross-reference to one STRING v12.0 protein (else UniProt\'s STRING cross-reference); empty when PaxDb has no value '
         'for the entry (an isoform takes its entry\'s value). For a detection covariate, not a measurement in the cells the sites '
         'come from'),
        ('paxdb_percentile', 'number', '', 'percentile of that abundance among the 18,805 reviewed entries with a value (ties averaged; '
         '100 the most abundant)'),
    ],
    protein_values=paxdb_protein_values, website=lambda release, version, inputs, site_out: paxdb_website(release, version, site_out),
    since='v1.4.0',
    summary='Protein abundance from PaxDb 6.1 (whole-organism integrated dataset) for each protein (`paxdb_abundance_ppm`, `paxdb_percentile`): '
            'a covariate for detection, since mass spectrometry finds abundant proteins most.',
    changes=lambda files: paxdb_changes(files),
)


# ---- Reactome pathways (release 97) -----------------------------------------------------------------------------------

REACTOME_MIN, REACTOME_MAX = 10, 500   # pathway sizes tested, in reviewed proteins with at least one Reactome pathway


@functools.lru_cache(maxsize=4)
def reactome(inputs: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    pairs = read_layer(inputs, 'reactome_human_protein_pathways.tsv')
    pathways = read_layer(inputs, 'reactome_human_pathways.tsv').set_index('pathway_id')
    return pairs, pathways


def reactome_protein_values(proteins: pd.DataFrame, inputs) -> pd.DataFrame:
    pairs, pathways = reactome(str(inputs))
    top = pathways.index[pathways['depth'] == '0']
    in_top = pairs[pairs['pathway_id'].isin(top)].groupby('canonical_UniProtAC')['pathway_id'].agg(lambda s: sorted(pathways.loc[s, 'name']))
    n_low = pairs[pairs['lowest_level'] == 'true'].groupby('canonical_UniProtAC')['pathway_id'].nunique()
    entry = proteins['uniprot_accession'].str.split('-').str[0]
    return pd.DataFrame({'reactome_top_level': [join(in_top[a]) if a in in_top.index else '' for a in entry],
                         'reactome_n_pathways': [int(n_low[a]) if a in n_low.index else None for a in entry]}, index=proteins.index)


def reactome_stats(sites: pd.DataFrame, inputs, sequences: dict, workdir: Path) -> list[dict]:
    from enrichment import Job, run
    pairs, pathways = reactome(str(inputs))
    universe = set(pairs['canonical_UniProtAC'])            # reviewed entries with at least one pathway
    target = {a.split('-')[0] for a in sites.loc[sites['counted'], 'uniprot_accession']} & universe
    top = set(pathways.index[pathways['depth'] == '0'])
    sets = pairs[['pathway_id', 'canonical_UniProtAC']]     # all levels: a pathway holds the proteins of its sub-pathways
    names = pathways['name'].rename_axis('pathway_id').reset_index()
    jobs = [Job('reactome_pathway', target, universe, sets=sets, names=names, min_size=REACTOME_MIN, max_size=REACTOME_MAX),
            Job('reactome_top_level', target, universe, sets=sets[sets['pathway_id'].isin(top)], names=names)]
    res = run(jobs, Path(workdir) / 'enricher')
    background = 'other reviewed human proteins with at least one Reactome 97 pathway'
    rows = []
    for job in ('reactome_top_level', 'reactome_pathway'):
        for r in sorted(res[job].to_dict('records'), key=lambda r: (r['pvalue'], r['ID'])):
            rows.append({'analysis': job, 'category': r['ID'], 'stratum': 'all', 'a': r['a'], 'b': r['b'], 'c': r['c'], 'd': r['d'],
                         'fold_enrichment': r['FoldEnrichment'], 'log2_fold_enrichment': math.log2(r['FoldEnrichment']),
                         'p': r['pvalue'], 'q': r['p.adjust'], 'method': 'enricher', 'background': background, 'n_strata': ''})
    return rows


REACTOME_LISTED_Q = 0.05   # pathways below the top level the website lists: q at most this


def reactome_changes(files: dict) -> list[str]:
    st = files['statistics'][1]
    q = pd.to_numeric(st['q'], errors='coerce')
    fold = pd.to_numeric(st['fold_enrichment'], errors='coerce')
    top, low = st['analysis'] == 'reactome_top_level', st['analysis'] == 'reactome_pathway'
    return [f'{int((top & (q <= 0.05) & (fold > 1)).sum())} of the {int(top.sum())} top-level pathways, and {int((low & (q <= 0.05) & (fold > 1)).sum()):,} '
            f'of the {int(low.sum()):,} pathways of {REACTOME_MIN} to {REACTOME_MAX} proteins at any level, are enriched at q <= 0.05.']


def reactome_website(release: Path, version: str, inputs: Path, site_out: Path) -> None:
    """SITE_OUT/reactome.json: the top-level pathways, each protein's, and the release statistics of the top level and of
    the pathways below it at q <= 0.05 (for the Overview: the whole dataset's tests and a selection's retest)."""
    _, pathways = reactome(str(inputs))
    top = pathways[pathways['depth'] == '0'].sort_values('name')
    index = {name: k for k, name in enumerate(top['name'])}
    p = release_table(release, f'human_methylarginine_proteins_{version}.tsv')
    prot = {}
    for a, names in zip(p['uniprot_accession'], p['reactome_top_level']):
        if names:
            prot[a.split('-')[0]] = sorted(index[x] for x in names.split('; '))
    fields = ['analysis', 'category', 'a', 'b', 'c', 'd', 'fold_enrichment', 'p', 'q']
    st = release_table(release, f'statistics_{version}.tsv')
    st = st[(st['analysis'] == 'reactome_top_level') | ((st['analysis'] == 'reactome_pathway') & (pd.to_numeric(st['q']) <= REACTOME_LISTED_Q))]
    num = lambda f, v: int(v) if f in 'abcd' else float(f'{float(v):.4g}') if f in ('fold_enrichment', 'p', 'q') else v  # noqa: E731
    listed = sorted(set(st.loc[st['analysis'] == 'reactome_pathway', 'category']))
    write_json(site_out / 'reactome.json', {
        'version': version, 'top': [[i, n] for i, n in zip(top.index, top['name'])], 'prot': dict(sorted(prot.items())),
        'stats_fields': fields, 'stats': [[num(f, r[f]) for f in fields] for r in st.to_dict('records')],
        'pathways': {i: [pathways.at[i, 'name'], pathways.at[i, 'top_level_names'], int(pathways.at[i, 'n_reviewed_proteins'])] for i in listed},
        'listed_q': REACTOME_LISTED_Q, 'sizes': [REACTOME_MIN, REACTOME_MAX]})


REACTOME = Layer(
    key='reactome', title='Reactome 97 pathways',
    files=['reactome_human_protein_pathways.tsv', 'reactome_human_pathways.tsv', 'reactome_summary.json'],
    staging=['scripts/48_stage_reactome.py'],
    licence='CC0 (Reactome data; reactome.org/license)',
    citation='Reactome release 97; Ragueneau et al., Nucleic Acids Res 54:D673-D681 (2026), doi:10.1093/nar/gkaf1223',
    protein_columns=[
        ('reactome_top_level', 'list<string>', '', 'the Reactome 97 top-level pathways the protein is in, at any level below them, sorted '
         '(an isoform takes its entry\'s); empty when Reactome gives the entry no pathway'),
        ('reactome_n_pathways', 'integer', '', 'the entry\'s lowest-level Reactome 97 pathways (UniProt2Reactome); empty when Reactome gives the entry none'),
    ],
    analyses={
        'reactome_pathway': f'methylated proteins (by entry) in a Reactome 97 pathway, against the other reviewed human proteins with at '
                            f'least one Reactome pathway; a pathway holds the proteins of all its sub-pathways (UniProt2Reactome_All_Levels), '
                            f'and pathways of {REACTOME_MIN} to {REACTOME_MAX} such proteins are tested. Pathways nest, so neighbouring '
                            f'tests share proteins (inputs/layers/reactome_human_pathways.tsv gives each pathway\'s parents). Methylated '
                            f'proteins are much more abundant than the others (PaxDb), and abundant proteins are the ones mass spectrometry '
                            f'finds most, so pathways of abundant proteins are expected to be enriched whatever methylation does',
        'reactome_top_level': 'the same for the 29 Reactome top-level pathways, whatever their size (a protein can be in several)',
    },
    protein_values=reactome_protein_values, stats=reactome_stats, website=reactome_website, since='v1.4.0',
    summary='Reactome release 97 pathways of each protein (`reactome_top_level`, `reactome_n_pathways`), and two over-representation '
            'tests of the methylated proteins against the reviewed proteins with a pathway (analyses `reactome_top_level`, `reactome_pathway`; '
            'clusterProfiler enricher, as the other protein sets).',
    changes=lambda files: reactome_changes(files),
)


# ---- Curated functional experiments (inputs/curation/) ----------------------------------------------------------------

CURATION_COLUMNS = [
    ('experiment_id', 'string', '', 'FE and four digits, never reused'),
    ('site_id', 'string', '', 'the site, accession:R position; empty when the experiment is not resolved to one site'),
    ('in_release', 'boolean', '', 'the site is a site of this release'),
    ('uniprot_accession', 'string', '', 'reviewed human accession'), ('gene', 'string', '', 'its primary gene name'),
    ('position', 'integer', '', 'the arginine; empty when the experiment is not resolved to one site'),
    ('region', 'string', '', 'for a row without a site, which arginines the experiment concerns'),
    ('species', 'string', '', 'of the cells or protein in the experiment'), ('system', 'string', '', 'cells or organism, or the in vitro system'),
    ('writer', 'string', '', 'the PRMT the paper names as methylating the site, empty when none'),
    ('experiment_kind', 'code', 'experiment_kind', 'what kind of experiment; empty when only the abstract was read'),
    ('perturbation', 'string', '', 'what was done, as the paper names it'), ('perturbation_note', 'string', '', 'what else the perturbation changes'),
    ('readout', 'string', '', 'what was measured'), ('direction', 'code', 'direction', 'of the readout under the perturbation'),
    ('finding', 'string', '', 'one sentence, close to the paper\'s words'), ('figure', 'string', '', 'where in the paper'),
    ('reference', 'string', '', 'the paper'), ('doi', 'string', '', 'its DOI'), ('pmid', 'string', '', 'its PubMed id'),
    ('read_from', 'code', 'read_from', 'what the curator read'), ('curated_on', 'string', '', 'date of curation'), ('note', 'string', '', 'anything else'),
]


def curation_supporting(sites: pd.DataFrame, inputs, version: str) -> dict:
    t = pd.read_csv(Path(inputs) / 'curation' / 'functional_experiments.tsv', sep='\t', dtype=str, keep_default_na=False)
    ids = set(sites['uniprot_accession'] + ':R' + sites['position'].astype(int).astype(str))
    t['site_id'] = [f'{a}:R{p}' if p else '' for a, p in zip(t['uniprot_accession'], t['position'])]
    t['in_release'] = t['site_id'].isin(ids)
    return {'functional_experiments': (f'supporting/functional_experiments_{version}.tsv', t, CURATION_COLUMNS)}


def curation_website(release: Path, version: str, site_out: Path) -> None:
    """SITE_OUT/curation.json: the curated experiments as the release table has them."""
    t = release_table(release, f'supporting/functional_experiments_{version}.tsv')
    write_json(site_out / 'curation.json', {'version': version, 'fields': list(t.columns), 'rows': t.values.tolist()})


CURATION = Layer(
    key='curation', title='Curated functional experiments', folder='curation',
    files=['functional_experiments.tsv'], staging=['curated by hand; scripts/49_check_functional_curation.py'],
    licence='CC BY-NC-SA 4.0 (this repository\'s data licence); the facts are from the papers cited', citation='the paper of each row',
    codes={'experiment_kind': {'methylated_vs_unmethylated': 'the same protein or peptide with and without (or with less) methylation',
                               'site_mutation': 'the arginine mutated (R to K keeps the charge, R to A removes it; both remove the methyl group)',
                               'prmt_inhibitor': 'a PRMT inhibitor in cells', 'prmt_knockdown': 'a PRMT knocked down',
                               'prmt_knockout': 'a PRMT knocked out', 'correlation': 'methylation and a readout measured together, no perturbation',
                               'prediction': 'a computational prediction'},
           'direction': {'increased': 'the readout rose', 'decreased': 'the readout fell', 'no_change': 'the readout did not change'},
           'read_from': {'full_text': 'the full text was read', 'abstract': 'only the abstract was read'}},
    supporting=curation_supporting, website=lambda release, version, inputs, site_out: curation_website(release, version, site_out),
    since='v1.4.0',
    summary='Experiments from the literature on the methylation of a site or protein (`supporting/functional_experiments`), curated by hand '
            'from the papers (inputs/curation/README.md): what was done, what was measured and which way it changed. Kept apart from the '
            'site table and every statistic.',
    changes=lambda files: [f"{len(files['functional_experiments'][1]):,} experiments so far, from "
                           f"{files['functional_experiments'][1]['doi'].nunique()} papers."],
    supporting_rows={'functional_experiments': 'curated experiment on a site or a protein\'s methylation: perturbation, readout, '
                                               'direction and reference (kept apart from every statistic)'},
)

LAYERS = [
    Layer('cplm4', 'CPLM 4.0 lysine methylation', ['cplm4_summary.json'], ['scripts/39_stage_cplm4_lys_methylation.py'],
          'All rights reserved (© The CUCKOO Workgroup; no licence stated): the tables stay out of Git, checksum_only in MANIFEST.tsv',
          'Zhang et al., Nucleic Acids Res 50:D451 (2022)', elsewhere='the lysine methylation sites the website draws on the protein maps (scripts/36)'),
    Layer('carm1_ko', 'CARM1 knockout response (Shishkova et al. 2017)', ['carm1_ko_peptides.tsv', 'carm1_ko_sites.tsv', 'carm1_ko_summary.json'],
          ['scripts/41_stage_carm1_ko_sites.py'], 'CC BY 4.0', 'Shishkova et al., Nat Commun 8:15571 (2017)',
          elsewhere='site column `carm1_ko_call` and `supporting/carm1_ko` (scripts/build_release.py)'),
    Layer('hpa', 'Human Protein Atlas 25.1 subcellular locations', ['hpa_subcellular_locations.tsv', 'hpa_summary.json'],
          ['scripts/42_stage_hpa_subcellular.py'], 'CC BY 4.0', 'Thul et al., Science 356:eaal3321 (2017)',
          elsewhere='protein column `hpa_locations`, `supporting/hpa_locations` and analysis `hpa_location` (scripts/release_tables.py)'),
    CANCER,
    PERTURBATION,
    PAXDB,
    REACTOME,
    CURATION,
]


def readme_section() -> list[str]:
    """The release README's section on the layers (scripts/release_docs.py)."""
    code = lambda names: ', '.join(f'`{c}`' for c in names)  # noqa: E731
    lines = ['## Annotation layers', '',
             'Data from other resources, joined to the sites and proteins by accession and position; none adds a site. Each is '
             'staged into `inputs/layers/` (its files listed, with checksums, in `inputs/MANIFEST.tsv`) and registered in '
             '`scripts/release_layers.py`, which the build reads for the columns and tests a layer adds '
             '(`ADDING_A_SOURCE.md`).', '',
             '| Layer | Staged by | Adds | Licence |', '|---|---|---|---|']
    for layer in LAYERS:
        adds = ([f'site columns {code(c[0] for c in layer.site_columns)}'] if layer.site_columns else []) \
            + ([f'protein columns {code(c[0] for c in layer.protein_columns)}'] if layer.protein_columns else []) \
            + ([f'analyses {code(layer.analyses)} in the statistics table'] if layer.analyses else []) \
            + ([f'supporting tables {code(f"supporting/{k}" for k in layer.supporting_rows)}'] if layer.supporting_rows else []) \
            + ([layer.elsewhere] if layer.elsewhere else [])
        lines.append(f'| {layer.title} | {code(Path(s).name for s in layer.staging)} | {"; ".join(adds)} | {layer.licence} |')
    lines += ['', 'Cite them as: ' + '; '.join(f'{layer.title}: {layer.citation}' for layer in LAYERS) + '.', '']
    return lines


def changes_section(files: dict, version: str, previous: str) -> list[str]:
    """CHANGES.md: the layers added after the previous release, with their counts in this one."""
    v = lambda x: tuple(int(k) for k in x.lstrip('v').split('.'))  # noqa: E731
    new = [layer for layer in LAYERS if layer.since and v(previous) < v(layer.since) <= v(version)]
    if not new:
        return []
    out = ['## New annotation layers', '',
           f'{version} adds {len(new)} layers of annotation beside the sites. None adds, removes or moves a site; their columns, tables '
           'and tests are in `README.md`, `data_dictionary.tsv` and the statistics table, and scripts/release_layers.py builds them.', '']
    for layer in new:
        out += [f'### {layer.title}', '', layer.summary, '']
        if layer.changes:
            out += [f'- {x}' for x in layer.changes(files)] + ['']
        out += [f'Source: {layer.citation}. Licence: {layer.licence}.', '']
    return out


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


def supporting_tables(sites: pd.DataFrame, inputs, version: str) -> dict:
    """{key: (path in the release, frame, columns)} of every layer's supporting tables (build_release.py writes them)."""
    return {k: v for layer in LAYERS if layer.supporting for k, v in layer.supporting(sites, inputs, version).items()}


def supporting_file_rows() -> dict[str, str]:
    """What one row of each layer's supporting table is (release_docs.FILE_ROWS)."""
    return {k: v for layer in LAYERS for k, v in layer.supporting_rows.items()}


def statistics(sites: pd.DataFrame, inputs, sequences: dict, workdir: Path) -> list[dict]:
    """sites: uniprot_accession, position and counted (the analyses count the site) of every site of the table."""
    return [r for layer in LAYERS if layer.stats for r in layer.stats(sites, inputs, sequences, Path(workdir) / 'layers' / layer.key)]


# ---- Command line ------------------------------------------------------------------------------------------------------

def check(inputs: Path) -> list[str]:
    manifest = pd.read_csv(inputs / 'MANIFEST.tsv', sep='\t', dtype=str, keep_default_na=False).set_index('file')
    bad = []
    for layer in LAYERS:
        for name in layer.files:
            rel, path = f'{layer.folder}/{name}', inputs / layer.folder / name
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

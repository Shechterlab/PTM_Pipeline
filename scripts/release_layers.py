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


# ---- Cancer cell lines and tumours (DepMap 24Q4, PRISM 24Q2, TCGA MC3 from the GDC) -----------------------------------

DEPENDENT = 0.5        # DepMap's dependency probability at which a line depends on a gene
MOST_LINES = 0.5       # dependency classes: the share of screened lines that depend on the gene
SOME_LINES = 0.05
CHANGES = (  # source, staged table, what counts once, site columns (changes, count per change)
    ('tcga', 'tcga_mc3_arginine_changes.tsv', 'sample_id', 'tcga_changes', 'tcga_change_tumours'),
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
           'depmap_24q4_prmt_gene_effect.tsv', 'prism_24q2_prmt_inhibitors.tsv', 'tcga_mc3_arginine_changes.tsv',
           'tcga_mc3_projects.tsv', 'depmap_summary.json'],
    staging=['scripts/44b_extract_tcga_mc3.py', 'scripts/44_stage_depmap.py'],
    licence='DepMap Public 24Q4 and PRISM Repurposing 24Q2: CC BY 4.0; TCGA MC3 calls: GDC open-access data (NIH Genomic Data Sharing '
            'policy, no licence stated)',
    citation='DepMap Public 24Q4 (doi:10.25452/figshare.plus.27993248; Tsherniak et al. 2017, doi:10.1016/j.cell.2017.06.010; Chronos, '
             'Dempster et al. 2021, doi:10.1186/s13059-021-02540-7); PRISM Repurposing 24Q2 (doi:10.6084/m9.figshare.25917643; Corsello et al. '
             '2020, doi:10.1038/s43018-019-0018-6); the results shown here are in whole or part based upon data generated by the TCGA Research '
             'Network (https://www.cancer.gov/tcga): MC3 somatic calls (Ellrott et al. 2018, doi:10.1016/j.cels.2018.03.002), the calls of '
             'the PanCancer Atlas (Hoadley et al. 2018, doi:10.1016/j.cell.2018.03.022), from the GDC (mc3.v0.2.8.PUBLIC.maf.gz)',
    site_columns=[
        ('tcga_changes', 'list<string>', '', 'missense and stop changes of the arginine in TCGA tumours (MC3 public calls from the GDC, 33 '
         'projects, 10,294 tumour samples; FILTER PASS, and in OV and LAML also wga alone, as MC3 used them), most tumours first, * a stop. '
         'A change is read on a transcript UniProt assigns to the canonical isoform where the variant falls in one, else on the MAF\'s '
         'transcript, and counts only where the canonical sequence has this arginine at that position (inputs/layers/README.md)'),
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
        'cancer_change': 'methylarginines with a missense or stop change in at least one TCGA tumour (MC3 calls, category tcga_tumours) or '
                         'DepMap 24Q4 cell line (depmap_cell_lines), against the other arginines of the methylated proteins whose gene DepMap '
                         'screened: pooled (method fisher) and within each protein (method mantel_haenszel, each protein a stratum). Sites of '
                         'the table that the analyses leave out are in neither group. Arginine codons CGN hold a CpG, the dinucleotide that '
                         'mutates most often in tumours, so codon use as well as selection shapes the comparison',
    },
    site_values=cancer_site_values, protein_values=cancer_protein_values, stats=cancer_stats, website=cancer_website,
    topic='9F1239', since='v1.4.0',
    summary='Cancer genetics beside the sites, not evidence about methylation: how much DepMap\'s CRISPR-screened cancer cell lines '
            'depend on each methylated protein (protein columns `depmap_*`), and the missense and stop changes of each methylarginine '
            'in TCGA tumours (the MC3 calls of the PanCancer Atlas, from the GDC) and DepMap cell lines (site columns `tcga_changes`, '
            '`depmap_changes`). Two tests: '
            'dependency classes of the methylated proteins (analysis `depmap_dependency`) and changes at methylarginines against the '
            'other arginines of their proteins (analysis `cancer_change`).',
    changes=lambda files: cancer_changes(files),
)

def release_table(release: Path, name: str) -> pd.DataFrame:
    return pd.read_csv(Path(release) / name, sep='\t', dtype=str, keep_default_na=False)


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, separators=(',', ':'), ensure_ascii=False), encoding='utf-8')


# ---- Writer perturbations: knockouts, knockdowns and inhibitors -------------------------------------------------------

RELEASE_Q, RELEASE_LOG2 = 0.05, 1.0   # the uniform call: twofold (|log2| >= 1) at q < 0.05, as for the CARM1 knockout
SOURCE_CALL = {'down': 'decreased', 'up': 'increased', 'unchanged': 'no_significant_change', 'mixed': 'mixed',
               'decreased': 'decreased', 'increased': 'increased', 'no_significant_change': 'no_significant_change'}
LI2021 = ('Li WJ et al., Nat Commun 12:1946 (2021)', '10.1038/s41467-021-21963-1')
LI2021_RULE = 'more than twofold lower methylation ratio (H/L < 0.5) in one SILAC experiment, listed per site'


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
    q: str               # the column of an adjusted p value to use as given, 'bh': p adjusted here, '': none (no twofold,
                         # q < 0.05 call)
    source_call: str     # the column of each peptide's published call
    source_rule: str     # the published criteria
    select: dict = field(default_factory=dict)       # column -> value: the rows of this arm in a shared table
    bh_groups: list = field(default_factory=list)    # columns within which p is adjusted ('bh')
    primary: str = 'twofold_q05'                     # which call is `call`: 'twofold_q05' or 'published'
    primary_why: str = ''                            # why the published call, when it is `call`
    source_sites: tuple = ()                         # (staged site table, column): the study's own call per site, where it calls
                                                     # sites rather than peptides (keyed by canonical_UniProtAC, position and select)
    localised: tuple = ()                            # (column, value): a peptide counts as single-site only when its one position
                                                     # has this value (positions aligned with methyl_positions)
    changed_only: bool = False                       # the study lists its changed sites only: no row is not "unchanged"


PERTURBATIONS = [
    Perturbation('carm1_ko_mcf7', 'CARM1 KO, MCF7', 'CARM1 knockout against parental MCF7 cells, ADMA immunoprecipitation, label-free '
                 '(Shishkova et al. 2017)', 'CARM1', 'knockout', 'MCF7', 'Shishkova et al., Nat Commun 8:15571 (2017)', '10.1038/ncomms15571',
                 'carm1_ko_peptides.tsv', 'log2_ko_wt', 'p', 'bh', 'call', 'twofold at nominal p < 0.05',
                 select={'cell_line': 'MCF7'}),
    Perturbation('carm1_ko_mdamb231', 'CARM1 KO, MDA-MB-231', 'CARM1 knockout against parental MDA-MB-231 cells, ADMA '
                 'immunoprecipitation, label-free (Shishkova et al. 2017)', 'CARM1', 'knockout', 'MDA-MB-231',
                 'Shishkova et al., Nat Commun 8:15571 (2017)', '10.1038/ncomms15571', 'carm1_ko_peptides.tsv', 'log2_ko_wt', 'p', 'bh', 'call',
                 'twofold at nominal p < 0.05', select={'cell_line': 'MDA-MB-231'}),
    Perturbation('prmt7_ko_hct116', 'PRMT7 KO, HCT116', 'PRMT7 knockout against parental HCT116 cells, monomethylarginine '
                 'immunoprecipitation, SILAC, 4 replicates (Szewczyk et al. 2020)', 'PRMT7', 'knockout', 'HCT116',
                 'Szewczyk et al., Nat Commun 11:2396 (2020)', '10.1038/s41467-020-16271-z', 'prmt7_ko_peptides.tsv', 'log2_ko_wt', 'p', 'q', 'call',
                 'twofold at BH-adjusted p < 0.01 (its adjusted p values are used as given)'),
    Perturbation('prmt1_kd_hek293t', 'PRMT1 KD, HEK293T', 'PRMT1 shRNA knockdown against a non-targeting shRNA in HEK293T cells, '
                 'label-free, 2 against 2 biological replicates, SCX and antibody enrichments (Hartel et al. 2019)', 'PRMT1', 'knockdown',
                 'HEK293T', 'Hartel et al., Mol Cell Proteomics 18:2149-2164 (2019)', '10.1074/mcp.RA119.001625',
                 'hartel2019_prmt1_kd_peptides.tsv', 'log2_kd_control', 'p', 'bh', 'call',
                 'Perseus t-test with permutation-based FDR q < 0.05 and S0 = 0.5', bh_groups=['table'], primary='published',
                 primary_why='2 against 2 samples, too few for a Benjamini-Hochberg-adjusted t-test, so the authors used a '
                             'permutation-based FDR'),
    Perturbation('prmt5i_epz015666_hela', 'PRMT5i EPZ015666, HeLa', 'PRMT5 inhibitor EPZ015666 against untreated HeLa cells, SDMA '
                 'antibody after ulilysin digestion, SILAC (Lu et al. 2024), the treated channel inferred from the data', 'PRMT5', 'inhibitor', 'HeLa', 'Lu et al., J Proteome Res 23:1014-1027 (2024)',
                 '10.1021/acs.jproteome.3c00724', 'lu2024_prmt5_peptides.tsv', 'log2_inhibitor_control', 'p', 'bh', 'call',
                 'twofold at nominal p < 0.05 (not stated in the paper, this threshold reproduces 254 of the 255 sites of its Table 1.4 '
                 'that map to the reference, and adds none)'),
    Perturbation('prmt7_kd_hek293', 'PRMT7 KD, HEK293', 'PRMT7 siRNA knockdown against control siRNA in HEK293 cells, '
                 'monomethylarginine antibody, SILAC, one experiment (Supplementary Data 1, repeat 1, Li et al. 2021)', 'PRMT7', 'knockdown',
                 'HEK293', *LI2021, 'li2021_peptides.tsv', 'log2_treated_control', '', '', 'call_source', LI2021_RULE,
                 select={'comparison': 'PRMT7_si'}, primary='published', source_sites=('li2021_sites.tsv', 'call_source'),
                 localised=('localised_by', 'search'),
                 primary_why='one SILAC experiment, with no test across replicates (the Proteome Discoverer adjusted p is not one)'),
    Perturbation('prmt7i_sgc3027_hek293', 'PRMT7i SGC3027, HEK293', 'PRMT7 inhibitor SGC3027 against untreated HEK293 cells, '
                 'monomethylarginine antibody, SILAC, one experiment (Supplementary Data 3, Li et al. 2021). Heavy is read as treated, '
                 'as the paper reads its ratios, though one sentence of its Methods says the reverse', 'PRMT7', 'inhibitor', 'HEK293', *LI2021,
                 'li2021_peptides.tsv', 'log2_treated_control', '', '', 'call_source', LI2021_RULE, select={'comparison': 'PRMT7_SGC3027'},
                 primary='published', source_sites=('li2021_sites.tsv', 'call_source'), localised=('localised_by', 'search'),
                 primary_why='one SILAC experiment, with no test across replicates'),
    Perturbation('carm1_kd_hek293', 'CARM1 KD, HEK293', 'CARM1 (PRMT4) siRNA knockdown against control siRNA in HEK293 cells, '
                 'monomethylarginine then ADMA antibodies, SILAC, one experiment (Supplementary Data 5, Li et al. 2021)', 'CARM1', 'knockdown',
                 'HEK293', *LI2021, 'li2021_peptides.tsv', 'log2_treated_control', '', '', 'call_source', LI2021_RULE,
                 select={'comparison': 'PRMT4_si'}, primary='published', source_sites=('li2021_sites.tsv', 'call_source'),
                 localised=('localised_by', 'search'),
                 primary_why='one SILAC experiment, with no p value'),
    Perturbation('prmt5_kd_hek293', 'PRMT5 KD, HEK293', 'PRMT5 siRNA knockdown against control siRNA in HEK293 cells, '
                 'monomethylarginine then SDMA antibodies, SILAC, one experiment (Supplementary Data 6, Li et al. 2021): the decreased '
                 'sites only, as the study lists them, so a site without a row may have been measured unchanged', 'PRMT5', 'knockdown',
                 'HEK293', *LI2021, 'li2021_peptides.tsv', 'log2_treated_control', '', '', 'call_source', LI2021_RULE,
                 select={'comparison': 'PRMT5_si'}, primary='published', source_sites=('li2021_sites.tsv', 'call_source'),
                 localised=('localised_by', 'search'), changed_only=True,
                 primary_why='one SILAC experiment, with no p value, and only the decreased sites listed'),
]


@functools.lru_cache(maxsize=4)
def perturbation_rows(inputs: str) -> pd.DataFrame:
    """One row per site and perturbation: the twofold, q < 0.05 call (|log2| >= 1 at q < 0.05, scripts/peptide_calls.py over
    the peptides carrying the site, q adjusted here where the study gives p only) and the published call; `call` is the one
    the perturbation names."""
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
        elif d.q:
            pep['rel_q'] = pd.to_numeric(pep[d.q], errors='coerce')
        else:
            pep['rel_q'] = float('nan')
        pep = pep[(pep['canonical_UniProtAC'] != '') & (pep['methyl_positions'] != '')]
        site_call = {}
        if d.source_sites:
            st = read_layer(inputs, d.source_sites[0])
            for col, value in d.select.items():
                st = st[st[col] == value]
            site_call = {(a, int(p)): SOURCE_CALL.get(c, '') for a, p, c in zip(st['canonical_UniProtAC'], st['position'], st[d.source_sites[1]])}
        groups = defaultdict(list)
        for r in pep.to_dict('records'):
            positions = [int(float(x)) for x in str(r['methyl_positions']).split(';') if x]
            how = str(r[d.localised[0]]).split(';') if d.localised else []
            single = len(positions) == 1 and (not d.localised or how[:1] == [d.localised[1]])
            for pos in positions:
                groups[(r['canonical_UniProtAC'], pos)].append((r['rel_q'], r['rel_log2'], len(positions), r['rel_p'],
                                                                SOURCE_CALL.get(r[d.source_call], ''), single))
        for (acc, pos), vals in sorted(groups.items()):
            res = pc.site_response([(q, l2, 1 if single else 2) for q, l2, _, _, _, single in vals], RELEASE_Q, RELEASE_LOG2, strict=True)
            tested = any(not pd.isna(q) for q, *_ in vals)
            ps = [p for _, _, _, p, _, _ in vals if not pd.isna(p)]
            release = res['call'] if tested else 'not_tested'
            release_single = (res['call_single_site_peptides'] or '') if tested else ''
            source = site_call.get((acc, pos), '') if d.source_sites else pc.combine_calls(c for *_, c, _ in vals)
            source_single = pc.combine_calls(c for *_, c, single in vals if single)
            if d.primary == 'published':
                call, single_call = source or 'not_tested', source_single
                n_down = sum(c == 'decreased' for *_, c, _ in vals)
                n_up = sum(c == 'increased' for *_, c, _ in vals)
            else:
                call, single_call, n_down, n_up = release, release_single, res['n_decreased'], res['n_increased']
            rows.append({'site_id': f'{acc}:R{pos}', 'perturbation': d.key, 'n_peptides': len(vals),
                         'log2_treated_over_control': round(res['log2_median'], 3), 'p': min(ps) if ps else None,
                         'q': res['adj_p_min'] if tested else None, 'n_decreased': n_down, 'n_increased': n_up,
                         'n_multi_site_peptides': sum(n > 1 for _, _, n, _, _, _ in vals), 'call': call, 'call_criteria': d.primary,
                         'call_single_site_peptides': single_call, 'call_twofold_q05': release, 'call_published': source})
    return pd.DataFrame(rows)


def perturbation_site_values(sites: pd.DataFrame, inputs) -> pd.DataFrame:
    t = perturbation_rows(str(inputs))
    order = {d.key: k for k, d in enumerate(PERTURBATIONS)}
    by = defaultdict(list)
    for r in t.itertuples(index=False):
        by[r.site_id].append((order[r.perturbation], r.perturbation, r.call, r.call_twofold_q05, r.call_published or 'not_tested'))
    ids = sites['uniprot_accession'] + ':R' + sites['position'].astype(int).astype(str)
    vals = [sorted(by.get(s, [])) for s in ids]
    return pd.DataFrame({'perturbations': [join(x[1] for x in v) for v in vals],
                         'perturbation_calls': [join(x[2] for x in v) for v in vals],
                         'perturbation_calls_twofold_q05': [join(x[3] for x in v) for v in vals],
                         'perturbation_calls_published': [join(x[4] for x in v) for v in vals]}, index=sites.index)


PERTURBATION_COLUMNS = [
    ('site_id', 'string', '', 'the site'),
    ('perturbation', 'code', 'perturbation', 'the knockout, knockdown or inhibitor, and the cells'),
    ('n_peptides', 'integer', '', 'quantified peptides carrying the site'),
    ('log2_treated_over_control', 'number', '', 'median log2 perturbed / control over those peptides: descriptive, not the effect of a test'),
    ('p', 'number', '', 'the lowest p value over the peptides, as the study gives it'),
    ('q', 'number', '', 'the lowest q over the peptides: the study\'s adjusted p where it gives one, else its p adjusted here (Benjamini-Hochberg, within the comparison); empty when the study gives no p'),
    ('n_decreased', 'integer', '', 'peptides called decreased by the criteria of `call`'),
    ('n_increased', 'integer', '', 'peptides called increased by the criteria of `call`'),
    ('n_multi_site_peptides', 'integer', '', 'of the peptides, those carrying another methylarginine too'),
    ('call', 'code', 'perturbation_call', 'the call: `call_twofold_q05` or `call_published`, as `call_criteria` says (the perturbation\'s code says why); peptide-level evidence'),
    ('call_criteria', 'code', 'perturbation_criteria', 'which of the two calls `call` is'),
    ('call_single_site_peptides', 'code', 'perturbation_call', 'the same over the peptides that carry no other methylarginine and whose '
     'one mark the search localised; empty when there are none'),
    ('call_twofold_q05', 'code', 'perturbation_call', f'the same criteria for every perturbation: a peptide changes at |log2| >= {RELEASE_LOG2:g} and q < {RELEASE_Q:g}, '
     'combined over the peptides carrying the site (scripts/peptide_calls.py); not_tested where the study gives no p value'),
    ('call_published', 'code', 'perturbation_call', 'by the published criteria (code perturbation); empty where the study makes no call '
     '(as for a dimethylated site in a monomethylation list)'),
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
        if not len(x):
            continue
        dec, inc = int((x['call'] == 'decreased').sum()), int((x['call'] == 'increased').sum())
        rule = 'the published criteria' if d.primary == 'published' else 'twofold at q < 0.05'
        tested = x[(x['call_twofold_q05'] != 'not_tested') & (x['call_published'] != '')]
        other = (f'; the two calls differ at {int((tested["call_twofold_q05"] != tested["call_published"]).sum()):,} of the '
                 f'{len(tested):,} sites both call') if len(tested) else '; no p values, so no twofold, q < 0.05 call'
        out.append(f'{d.label} ({d.study}): {len(x):,} sites{" listed" if d.changed_only else " quantified"}, {dec:,} decreased and '
                   f'{inc:,} increased by {rule}{other}.')
    return out


def perturbation_website(release: Path, version: str, inputs: Path, site_out: Path) -> None:
    """SITE_OUT/perturbations.json: the datasets, and for each site of the release its rows of supporting/writer_perturbations."""
    t = release_table(release, f'supporting/writer_perturbations_{version}.tsv')
    keys = [d.key for d in PERTURBATIONS]
    calls = list(PERTURBATION.codes['perturbation_call'])
    ix = lambda c: calls.index(c) if c else -1  # noqa: E731
    num = lambda v, nd=3: None if v == '' else float(f'{float(v):.{nd}g}')  # noqa: E731
    sites = defaultdict(list)
    for r in t.itertuples(index=False):
        sites[r.site_id].append([keys.index(r.perturbation), ix(r.call), ix(r.call_twofold_q05), ix(r.call_published),
                                 num(r.log2_treated_over_control), num(r.q), int(r.n_peptides), int(r.n_decreased), int(r.n_increased),
                                 int(r.n_multi_site_peptides), ix(r.call_single_site_peptides)])
    payload = {'fields': ['perturbation', 'call', 'call_twofold_q05', 'call_published', 'log2_treated_over_control', 'q', 'n_peptides',
                          'n_decreased', 'n_increased', 'n_multi_site_peptides', 'call_single_site_peptides'],
               'calls': calls, 'rule': f'|log2| >= {RELEASE_LOG2:g} at q < {RELEASE_Q:g}, combined over the peptides carrying the site',
               'perturbations': [[d.key, d.label, d.description, d.enzyme, d.kind, d.cell_line, d.study, d.doi, d.source_rule, d.primary,
                                  d.primary_why, int(d.changed_only)] for d in PERTURBATIONS],
               'version': version, 'sites': dict(sorted(sites.items()))}
    write_json(site_out / 'perturbations.json', payload)


def perturbation_code(d: Perturbation) -> str:
    primary = (f'Called by the published criteria ({d.primary_why}), with the twofold, q < 0.05 call beside it' if d.primary == 'published'
               else 'Called at twofold and q < 0.05, with the published call beside it')
    return f'{d.description}. Published criteria: {d.source_rule}. {primary}'


PERTURBATION = Layer(
    key='writer_perturbations', title='Writer perturbations (PRMT knockouts, knockdowns and inhibitors)',
    files=['prmt7_ko_peptides.tsv', 'prmt7_ko_sites.tsv', 'prmt7_ko_searched_sequences.tsv', 'prmt7_ko_summary.json',
           'hartel2019_prmt1_kd_peptides.tsv', 'hartel2019_prmt1_kd_sites.tsv', 'hartel2019_prmt1_kd_summary.json',
           'lu2024_prmt5_peptides.tsv', 'lu2024_prmt5_sites.tsv', 'lu2024_prmt5_summary.json',
           'li2021_peptides.tsv', 'li2021_sites.tsv', 'li2021_summary.json'],
    staging=['scripts/41_stage_carm1_ko_sites.py', 'scripts/45_stage_prmt7_szewczyk2020.py', 'scripts/46_stage_hartel2019_prmt1_kd.py',
             'scripts/51_stage_lu2024_prmt5.py', 'scripts/50_stage_li2021_perturbations.py'],
    licence='CC BY 4.0 (Shishkova et al. 2017; Szewczyk et al. 2020; Hartel et al. 2019, from 2020-09-15; Li et al. 2021, its Supplementary '
            'Data taken as covered by the article\'s licence); CC BY-NC 4.0 (Lu et al. 2024 supporting information)',
    citation='; '.join(sorted({f'{d.study}, doi:{d.doi}' for d in PERTURBATIONS})),
    site_columns=[
        ('perturbations', 'list<code>', 'perturbation', 'the knockouts, knockdowns and inhibitors that quantified a peptide carrying the '
         'site (or, for one that lists only its changed sites, listed it), in code order (supporting/writer_perturbations)'),
        ('perturbation_calls', 'list<code>', 'perturbation_call', 'for each of them, the call (twofold at q < 0.05, or the published call '
         'where its code says so): peptide-level evidence, same order'),
        ('perturbation_calls_twofold_q05', 'list<code>', 'perturbation_call', f'for each of them, the call at |log2| >= {RELEASE_LOG2:g} and q < {RELEASE_Q:g}, same order'),
        ('perturbation_calls_published', 'list<code>', 'perturbation_call', 'for each of them, the published call, same order'),
    ],
    codes={'perturbation': {d.key: perturbation_code(d) for d in PERTURBATIONS},
           'perturbation_call': {'decreased': 'a peptide carrying the site was lower with the perturbation: at least one decreased and none increased',
                                 'increased': 'a peptide carrying the site was higher with the perturbation: at least one increased and none decreased',
                                 'mixed': 'peptides changed in opposite directions',
                                 'no_significant_change': 'quantified, and no peptide changed by the rule: not evidence that the site is unchanged',
                                 'not_tested': 'quantified, but not called: the study gives no p value (twofold, q < 0.05 call) or '
                                               'makes no call on the site (published call)'},
           'perturbation_criteria': {'twofold_q05': f'the same criteria for every perturbation: |log2| >= {RELEASE_LOG2:g} at q < {RELEASE_Q:g}',
                                     'published': 'the published criteria, where the perturbation\'s code says why'}},
    site_values=perturbation_site_values, supporting=perturbation_supporting, website=perturbation_website, since='v1.4.0',
    summary='Knockouts, knockdowns and inhibitors of the PRMTs that quantified methylated peptides. Each site is called twice: at '
            f'twofold and q < {RELEASE_Q:g} (a peptide changes at |log2| >= {RELEASE_LOG2:g} and q < {RELEASE_Q:g}; a site is decreased when '
            'a peptide carrying it decreased and none increased: scripts/peptide_calls.py, as for the CARM1 knockout), and by the '
            'published criteria. `call` is the twofold, q < 0.05 call, except for the PRMT1 knockdown (2 against 2 samples) and the four '
            'Li et al. 2021 arms (one SILAC experiment each, no test), where it is the published call. Peptide-level evidence with no '
            'total-protein correction. Site columns `perturbations`, `perturbation_calls`, `perturbation_calls_twofold_q05`, '
            '`perturbation_calls_published`; one row per site and perturbation in `supporting/writer_perturbations`. The CARM1 knockout '
            'rows equal `supporting/carm1_ko`, whose column `call_source_rule` is renamed `call_published`.',
    changes=lambda files: perturbation_changes(files),
    supporting_rows={'writer_perturbations': 'site and knockout, knockdown or inhibitor (PRMT1, PRMT5, PRMT7, CARM1): its peptides, '
                                             'change, and the twofold, q < 0.05 and published calls'},
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
    st = st[st['method'] == 'enricher']   # the tests themselves, not the same sets within abundance deciles
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
    st = st[st['method'] == 'enricher']   # the within-abundance-decile rows reach the site through scripts/36 (abundance_strata)
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
                            f'least one Reactome pathway. A pathway holds the proteins of all its sub-pathways (UniProt2Reactome_All_Levels), '
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
                               'site_mutation': 'the arginine mutated (R to K keeps the charge and R to A removes it, and both remove the methyl group)',
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

# ---- Maps of the cell (scripts/cell_maps.py) ---------------------------------------------------------------------------

import cell_maps  # noqa: E402

CELL = Layer(
    key='cell', title='Maps of the cell',
    files=cell_maps.FILES, staging=['scripts/52_stage_cell_maps.py'],
    licence='Human Cell Map CC BY 4.0; organelle IP BSD-3; Kustatscher et al. 2014 CC BY; Barutcu et al. 2022 Table S1 (Elsevier): '
            'the calls derived from it only; GEO and ENCODE data without restrictions; GENCODE open',
    citation='Go et al., Nature 595:120 (2021), doi:10.1038/s41586-021-03592-2; Hein et al., Cell 188:1137 (2025), '
             'doi:10.1016/j.cell.2024.11.028; Kustatscher et al., EMBO J 33:648 (2014), doi:10.1002/embj.201387614; Barutcu et al., Mol Cell '
             '82:1035 (2022), doi:10.1016/j.molcel.2021.12.010; Fazal et al., Cell 178:473 (2019), doi:10.1016/j.cell.2019.05.027; Chen et al., '
             'J Cell Biol 217:4025 (2018), doi:10.1083/jcb.201807108; Rao et al., Cell 159:1665 (2014), doi:10.1016/j.cell.2014.11.021; '
             'ENCODE Project Consortium, Nature 583:699 (2020), doi:10.1038/s41586-020-2493-4',
    protein_columns=[
        ('hcm_compartment', 'string', '', 'the Human Cell Map compartment of the entry (Go et al. 2021, HEK293 BioID, NMF localisation); empty when the map does not localise it'),
        ('hcm_domain', 'string', '', 'its Human Cell Map SAFE domain, the finer neighbourhood of proteins that label one another; empty when SAFE gives it none'),
        ('organelle_ip_compartment', 'string', '', 'its organelle IP cluster (Hein et al. 2025, HEK293T), or unclassified; empty when organelle IP did not detect it'),
        ('chromatin_probability', 'number', '', 'its interphase chromatin probability (Kustatscher et al. 2014): how reliably it comes down with chromatin; empty when not scored'),
        ('nuclear_body_proteomes', 'list<string>', '', 'the nuclear speckles, lamina and PML bodies whose APEX2 proximity proteome holds it (Barutcu et al. 2022, '
         'SAINT AvgP >= 0.9; speckles with two of the three markers); none when detected in none; empty when not detected'),
        ('eclip_rna_compartments', 'list<string>', '', 'where the RNAs its eCLIP experiments bind gather, as map: compartment (q <= 0.05 and odds ratio above 1 '
         'within expression and gene-length deciles; supporting/eclip_rna_compartments); none when a map tests it and they gather nowhere; '
         'empty when no map tests it (no eCLIP data, or fewer than 20 bound RNAs placed)'),
    ],
    codes={'eclip_map': cell_maps.MAP_CODES},
    analyses={
        'cell_hcm_compartment': 'methylated proteins (by entry) in a Human Cell Map compartment (Go et al. 2021, HEK293 BioID, NMF), against the other '
                                'proteins the map localises (compartments of at least 5 proteins)',
        'cell_hcm_domain': 'the same for the Human Cell Map\'s SAFE domains, against the other proteins placed in one',
        'cell_organelle_ip': 'the same for the organelle IP clusters (Hein et al. 2025, HEK293T), against the other proteins it classifies',
        'cell_nuclear_body_proteome': 'the same for the nuclear speckle, lamina and PML-body APEX2 proteomes (Barutcu et al. 2022, SAINT AvgP >= 0.9), '
                                      'against the other proteins those experiments detected',
        'cell_chromatin': 'the same for four bins of the interphase chromatin probability (Kustatscher et al. 2014: >= 0.8 chromatin, 0.5-0.8, '
                          '0.2-0.5, < 0.2 not chromatin), against the other proteins scored',
        'cell_region_hcm': 'methylated proteins in a region of the cell (nucleus, cytoplasm outside the organelles, organelles and membranes), each '
                           'the Human Cell Map compartments it holds (scripts/cell_maps.py, REGION_OF, the nuclear envelope with the membranes), '
                           'against the other proteins the map localises',
        'cell_region_organelle_ip': 'the same for the regions of the organelle IP clusters, against the other proteins it classifies',
        'cell_region_hpa': 'the same for the regions of the Human Protein Atlas locations (cilia and sperm structures in none), against the other '
                           'proteins HPA locates (Approved or better)',
        'eclip_rna_apexseq': 'methylated proteins among the proteins with eCLIP data (ENCODE, ENCORE) whose bound RNAs gather in an APEX-seq compartment '
                             '(Fazal et al. 2019, HEK293T, supporting/eclip_rna_compartments), against the other eCLIP\'d proteins it tests. Method '
                             'mann_whitney compares the size of every tested protein\'s odds ratio, called or not',
        'eclip_rna_nuclear_body': 'the same for the nuclear bodies of nuclear-body APEX2-seq (Barutcu et al. 2022, Flp-In HEK293)',
        'eclip_gene_tsaseq': 'the same for the genes of the bound RNAs near speckles or the lamina by TSA-seq (Chen et al. 2018, K562)',
        'eclip_gene_hic': 'the same for the Hi-C subcompartment of the genes of the bound RNAs (Rao et al. 2014, GM12878)',
    },
    protein_values=cell_maps.protein_values, stats=cell_maps.stats, supporting=cell_maps.supporting, website=cell_maps.website,
    supporting_rows={'eclip_rna_compartments': 'protein with eCLIP data, map and compartment: whether the RNAs it binds gather there, within '
                                               'expression and gene-length deciles (every protein and compartment a map tests)'},
    since='v1.4.2',
    summary='Where the methylated proteins sit in the cell, by four maps beside the Human Protein Atlas: the Human Cell Map and its finer '
            'domains, organelle IP, the nuclear speckle, lamina and PML-body proteomes, and chromatin probability (protein columns '
            '`hcm_compartment`, `hcm_domain`, `organelle_ip_compartment`, `nuclear_body_proteomes`, `chromatin_probability`). Each map tests '
            'the methylated proteins against the proteins it places, with enricher and within PaxDb abundance deciles. And where the RNAs '
            'that proteins with eCLIP data bind gather, by two RNA maps (APEX-seq, nuclear-body APEX2-seq) and two DNA maps (TSA-seq, Hi-C '
            'subcompartments), each protein\'s bound RNAs against the other RNAs a map places within expression and gene-length deciles '
            '(`supporting/eclip_rna_compartments`, column `eclip_rna_compartments`); the methylated proteins are compared with the other '
            'eCLIP\'d proteins in each compartment. The website showed these maps as a preview before this release.',
    changes=cell_maps.changes,
)

# ---- Provenance of the interaction records (scripts/interaction_sources.py) --------------------------------------------

import interaction_sources  # noqa: E402

SOURCES = Layer(
    key='interaction_sources', title='Provenance of the interaction records', folder='interactions',
    files=interaction_sources.FILES + [f'intact_{n}_human_{interaction_sources.INTACT_FEATURES}.tsv' for n in ('methylarginines', 'mutations')],
    staging=['scripts/53_stage_interaction_sources.py', 'scripts/37a_stage_interactions.py'],
    licence='PDB data CC0 1.0 (RCSB Data API and ModelServer); IntAct CC BY 4.0; PSI-MI and PSI-MOD (OLS) CC BY 4.0; PubMed citation data '
            'under the NLM Terms and Conditions (courtesy of the U.S. National Library of Medicine)',
    citation='Burley et al., Nucleic Acids Res 51:D488 (2023); Sehnal et al., Nucleic Acids Res 49:W431 (2021); del Toro et al., Nucleic '
             'Acids Res 50:D648 (2022); Montecchi-Palazzi et al., Nat Biotechnol 26:864 (2008); McLaughlin et al., Nucleic Acids Res 53:D1318 '
             '(2025); Sayers et al., Nucleic Acids Res 53:D20 (2025), doi:10.1093/nar/gkae979',
    site_columns=interaction_sources.SITE_COLUMNS, topic='5B21B6',
    codes={'interaction_class': interaction_sources.CLASSES, 'record_partner': interaction_sources.PARTNER_KINDS,
           'site_state': interaction_sources.STATES, 'record_flag': interaction_sources.FLAGS},
    elsewhere='its table and site columns are filled by scripts/build_release.py from supporting/site_interactions',
    supporting_rows={'site_interaction_records': 'interaction record at a site: its PDB entry, IntAct interaction, study or model, method, '
                                                 'species, paper, the state of the site in the experiment and what the record shows'},
    website=interaction_sources.website, since='v1.4.3',
    summary='Where every interaction record at a site comes from, and what it shows there. Each PDB row of '
            '`supporting/site_interactions` is split into its entries and each IntAct row into its interactions; the crosslinks and '
            'predicted models keep one row each. Every record gets its method, the species of both proteins, its paper, the state of the '
            'site in that experiment (the residue a structure models; the methylarginine IntAct curates there; a mutation) and the other '
            'marks of the protein, a class (writes, removes, binds methylated, binds unmodified, a mutation changing the binding, a binding '
            'region, another reaction, RNA or DNA, not a protein, another copy, predicted) and flags for what to check '
            '(`supporting/site_interaction_records`). The PDB entries that model a site as a methylarginine are read in full, '
            'including peptides that PDBe-KB does not map to the human protein.',
    changes=interaction_sources.changes,
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
    CELL,
    SOURCES,
]


def readme_section() -> list[str]:
    """The release README's section on the layers (scripts/release_docs.py)."""
    code = lambda names: ', '.join(f'`{c}`' for c in names)  # noqa: E731
    lines = ['## Annotation layers', '',
             'Data from other resources, joined to the sites and proteins by accession and position; none adds a site. Each is '
             'staged into `inputs/layers/` (the provenance of the interaction records into `inputs/interactions/`; its files listed, with '
             'checksums, in `inputs/MANIFEST.tsv`) and registered in '
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
           (f'{version} adds {len(new)} layers of annotation beside the sites. None adds, removes or moves a site; their columns, tables '
            'and tests are in `README.md`, `data_dictionary.tsv` and the statistics table, and scripts/release_layers.py builds them.' if len(new) > 1 else
            f'{version} adds a layer of annotation beside the sites. It adds, removes or moves no site; its columns, tables and tests are in '
            '`README.md`, `data_dictionary.tsv` and the statistics table, and scripts/release_layers.py builds them.'), '']
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
    """The layer files against MANIFEST.tsv, and every code meaning and study text free of the release's separators (';'
    between list items and in data_dictionary.tsv value_meanings; tab and newline anywhere), which the build would
    otherwise meet only when it writes the tables."""
    import release_codes as rc
    manifest = pd.read_csv(inputs / 'MANIFEST.tsv', sep='\t', dtype=str, keep_default_na=False).set_index('file')
    bad = [f'code {k}.{c} holds ";", a tab or a newline' for k, codes in rc.CODES.items() for c, m in codes.items()
           if any(x in str(m) for x in (';', '\t', '\n'))]
    bad += [f'{name} {k} holds ";", a tab or a newline' for name in ('STUDY_LICENCE', 'STUDY_NOTE', 'STUDY_QUALITY')
            for k, v in getattr(rc, name, {}).items() if any(x in str(v) for x in (';', '\t', '\n'))]
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
            raise SystemExit('the layer inputs or texts are not ready for a release:\n  ' + '\n  '.join(bad))
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

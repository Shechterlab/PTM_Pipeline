#!/usr/bin/env python
"""Gene-set tests: are ALS genes / TDP-43 targets more likely to encode methylarginine proteins?

1. ClinGen panels. Genes classified Definitive/Strong/Moderate, grouped by curation panel (GCEP);
   panels with >=10 genes. For each panel: fraction of genes with methylarginine (any evidence),
   rank among panels, and an adjusted odds ratio from logistic regression
       methyl ~ in_panel + rna_binding + log10_length + log10_abundance
   fitted on all curated disease genes (controls for the bias toward well-studied genes).
2. TDP-43 targets from the external cryptic-junction compendium (Shechterlab/cryptic_exon,
   external_tdp43_compendium_gene_summary.hg38.tsv). Same model, fitted on the MS-detectable
   background, with the ClinGen ALS panel as a positive comparator under the same background.
"""
from __future__ import annotations

import argparse
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get('PTM_DATA_ROOT', ROOT / 'PTM_data'))
RESULTS = Path(os.environ.get('PTM_RESULTS_ROOT', ROOT / 'PTM_results'))
OUT = RESULTS / 'disease_gene_tests'

ALS_PANEL = 'Amyotrophic Lateral Sclerosis Spectrum Disorders Gene Curation Expert Panel'
VALID = {'Definitive', 'Strong', 'Moderate'}
COVARIATES = ['rna_binding', 'log10_length', 'log10_abundance']


def load_genes(out_dir: Path) -> pd.DataFrame:
    g = pd.read_csv(out_dir / 'gene_table.tsv', sep='\t')
    rb = pd.read_csv(out_dir / 'gene_rna_binding_go0003723.tsv', sep='\t')
    g = g.merge(rb, on='entrez')
    for c in ['methyl_any', 'ms_detectable', 'rna_binding']:
        g[c] = g[c].astype(str).eq('True')
    g['heavy_methyl'] = g['arg_methyl_sites'] >= 3
    return g


def adjusted_or(df: pd.DataFrame, indicator: str, outcome: str = 'methyl_any') -> dict:
    X = sm.add_constant(df[[indicator] + COVARIATES].astype(float))
    y = df[outcome].astype(float)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        fit = sm.Logit(y, X).fit(disp=0)
    lo, hi = np.clip(fit.conf_int().loc[indicator], -50, 50)  # perfect separation in tiny panels
    return {'adj_OR': float(np.exp(fit.params[indicator])), 'adj_OR_lo95': float(np.exp(lo)),
            'adj_OR_hi95': float(np.exp(hi)), 'adj_p': float(fit.pvalues[indicator])}


def crude(df: pd.DataFrame, indicator: str, outcome: str = 'methyl_any') -> dict:
    a = int((df[indicator] & df[outcome]).sum())
    b = int((df[indicator] & ~df[outcome]).sum())
    c = int((~df[indicator] & df[outcome]).sum())
    d = int((~df[indicator] & ~df[outcome]).sum())
    odds, p = fisher_exact([[a, b], [c, d]])
    return {'n_set': a + b, 'n_methyl': a, 'frac_methyl': a / max(a + b, 1),
            'background_frac_methyl': c / max(c + d, 1), 'fisher_OR': odds, 'fisher_p': p}


def clingen_panels(genes: pd.DataFrame, clingen_csv: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    cg = pd.read_csv(clingen_csv, skiprows=[0, 1, 2, 3, 5])
    cg = cg[cg['CLASSIFICATION'].isin(VALID)]
    cg = cg.merge(genes[['symbol', 'entrez']], left_on='GENE SYMBOL', right_on='symbol')
    panel_genes = cg.groupby('GCEP')['entrez'].apply(set)
    panel_genes = panel_genes[panel_genes.map(len) >= 10]
    curated = genes[genes['entrez'].isin(set(cg['entrez']))].copy()
    curated[['entrez', 'symbol', 'methyl_any', 'arg_methyl_sites']].to_csv(
        OUT / 'clingen_curated_genes.tsv', sep='\t', index=False)

    rows = []
    for panel, members in panel_genes.items():
        curated['in_panel'] = curated['entrez'].isin(members)
        row = {'panel': panel.replace(' Gene Curation Expert Panel', '')}
        for outcome in ['methyl_any', 'heavy_methyl']:
            cr = crude(curated, 'in_panel', outcome)
            ad = adjusted_or(curated, 'in_panel', outcome)
            prefix = '' if outcome == 'methyl_any' else 'heavy_'
            row.update({f'{prefix}{k}': v for k, v in {**cr, **ad}.items()})
        row['n_rna_binding'] = int(curated.loc[curated['in_panel'], 'rna_binding'].sum())
        row['n_methyl_not_rna_binding'] = int((curated['in_panel'] & curated['methyl_any'] & ~curated['rna_binding']).sum())
        row['methyl_genes'] = ','.join(sorted(curated.loc[curated['in_panel'] & curated['methyl_any'], 'symbol']))
        rows.append(row)
    res = pd.DataFrame(rows)
    res['rank_frac_methyl'] = res['frac_methyl'].rank(ascending=False, method='min').astype(int)
    res['rank_adj_OR'] = res['adj_OR'].rank(ascending=False, method='min').astype(int)
    res['heavy_rank_frac'] = res['heavy_frac_methyl'].rank(ascending=False, method='min').astype(int)
    res = res.sort_values('frac_methyl', ascending=False)
    als = curated['entrez'].isin(panel_genes[ALS_PANEL])
    return res, genes.assign(als_clingen=genes['entrez'].isin(panel_genes[ALS_PANEL])), len(curated), als


def tdp43_tests(genes: pd.DataFrame, compendium: Path) -> pd.DataFrame:
    comp = pd.read_csv(compendium, sep='\t')
    comp = comp.merge(genes[['symbol', 'entrez']], left_on='gene_name', right_on='symbol')
    # Gene sets are defined from the compendium's own flags, fixed before looking at methylation.
    sets = {
        'All compendium host genes (any TDP-43-associated junction)': comp,
        'TDP-43-sensitive cryptic events': comp[comp['has_tdp43_sensitive'].astype(str).eq('True')],
        'Pathology-discriminating events (AUC >= 0.8)': comp[comp['max_auc'] >= 0.8],
        'NMD-rescued or dark cryptic events': comp[comp['has_nmd_or_dark'].astype(str).eq('True')],
        'Down after TDP-43 knockdown (>=2 datasets)': comp[comp['significant_down_dataset_count'] >= 2],
    }
    rows = []
    for bg_name, bg in [('MS-detectable', genes[genes['ms_detectable']]), ('Reviewed proteome', genes)]:
        bg = bg.copy()
        named = {**{k: set(v['entrez']) for k, v in sets.items()},
                 'Inherited ALS genes (ClinGen ALS panel)': set(genes.loc[genes['als_clingen'], 'entrez'])}
        for name, members in named.items():
            bg['in_set'] = bg['entrez'].isin(members)
            for outcome in ['methyl_any', 'heavy_methyl']:
                rows.append({'gene_set': name, 'background': bg_name, 'outcome': outcome,
                             **crude(bg, 'in_set', outcome), **adjusted_or(bg, 'in_set', outcome)})
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out-dir', type=Path, default=OUT)
    ap.add_argument('--clingen', type=Path, default=DATA / 'clingen' / 'gene_validity.csv')
    ap.add_argument('--tdp43-compendium', type=Path,
                    default=ROOT.parent / 'cryptic_exon' / 'cryptic_sf3b2_screen' / 'results'
                    / 'external_tdp43_compendium_gene_summary.hg38.tsv')
    args = ap.parse_args()

    genes = load_genes(args.out_dir)
    panels, genes, n_curated, _ = clingen_panels(genes, args.clingen)
    panels.to_csv(args.out_dir / 'clingen_panel_methylarginine.tsv', sep='\t', index=False)
    tdp = tdp43_tests(genes, args.tdp43_compendium)
    tdp.to_csv(args.out_dir / 'tdp43_compendium_methylarginine.tsv', sep='\t', index=False)

    als = panels[panels['panel'].str.startswith('Amyotrophic')].iloc[0]
    print(f'Curated disease genes (Definitive/Strong/Moderate, in universe): {n_curated}; panels >=10 genes: {len(panels)}')
    print(f"ALS panel: {als['n_methyl']}/{als['n_set']} methyl; rank {als['rank_frac_methyl']} of {len(panels)}; "
          f"adj OR {als['adj_OR']:.2f} ({als['adj_OR_lo95']:.2f}-{als['adj_OR_hi95']:.2f}), p={als['adj_p']:.3g}; "
          f"heavy rank {als['heavy_rank_frac']}")
    print(panels[['panel', 'n_set', 'n_methyl', 'frac_methyl', 'adj_OR', 'adj_p', 'rank_frac_methyl']].head(8).to_string(index=False))
    cols = ['gene_set', 'background', 'n_set', 'n_methyl', 'frac_methyl', 'background_frac_methyl', 'adj_OR', 'adj_OR_lo95', 'adj_OR_hi95', 'adj_p']
    print(tdp[tdp['outcome'].eq('methyl_any')][cols].to_string(index=False))


if __name__ == '__main__':
    main()

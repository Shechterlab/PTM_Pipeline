#!/usr/bin/env python
"""Export the data behind the ALS talk slides (script 41) as one JSON file for the editable
PowerPoint builder (scripts/43_build_als_talk_pptx.js). Same inputs and numbers as script 41."""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get('PTM_DATA_ROOT', ROOT / 'PTM_data'))
RESULTS = Path(os.environ.get('PTM_RESULTS_ROOT', ROOT / 'PTM_results'))
D = RESULTS / 'disease_gene_tests'
OUT = D / 'als_talk_figures' / 'als_talk_slide_data.json'

ALS_KEGG = {'Amyotrophic lateral sclerosis', 'Pathways of neurodegeneration - multiple diseases'}
PROTEINS = [('FUS', 'P35637', 'FUS'), ('EWSR1', 'Q01844', 'EWSR1'), ('TAF15', 'Q92804', 'TAF15'),
            ('TARDBP', 'Q13148', 'TDP-43'), ('HNRNPA1', 'P09651', 'hnRNPA1'),
            ('HNRNPA2B1', 'P22626', 'hnRNPA2B1'), ('MATR3', 'P43243', 'MATR3')]
KEGG_SHOWCASE = ['TARDBP', 'FUS', 'TAF15', 'EWSR1', 'HNRNPA1', 'HNRNPA2B1', 'MATR3', 'ATXN2', 'VCP', 'SOD1',
                 'UBQLN2', 'SQSTM1', 'G3BP1', 'NEFH', 'PRPH', 'DCTN1', 'SETX', 'TUBA4A']
DO_SHOWCASE = ['FUS', 'TAF15', 'EWSR1', 'TARDBP', 'HNRNPA1', 'MATR3', 'VCP', 'SOD1', 'ANXA11', 'NEFM',
               'DYNC1H1', 'BAG3', 'HSPB1', 'UBA1']


def kegg() -> dict:
    kg = pd.read_csv(D / 'kegg_2022_rerun' / 'enrichKEGG_ptm_pipeline_union.tsv', sep='\t')
    kg['count'] = kg['GeneRatio'].str.split('/').str[0].astype(int)
    kg['gene_ratio'] = kg['count'] / kg['GeneRatio'].str.split('/').str[1].astype(int)
    by_p = kg.sort_values('pvalue').reset_index(drop=True)
    by_p['rank'] = np.arange(1, len(by_p) + 1)
    als = by_p[by_p['Description'] == 'Amyotrophic lateral sclerosis'].iloc[0]
    genes = set(als['geneID'].split('/'))
    rows = lambda df: [{'name': r['Description'].replace(' - multiple diseases', ''), 'count': int(r['count']),
                        'neglog10fdr': round(float(-np.log10(r['p.adjust'])), 2),
                        'gene_ratio': round(float(r['gene_ratio']), 4), 'fdr': float(r['p.adjust']),
                        'highlight': r['Description'] in ALS_KEGG}
                       for r in df.to_dict('records')]
    return {
        'n_query': int(als['GeneRatio'].split('/')[1]), 'n_significant': len(kg),
        'top_by_p': rows(by_p.head(12)),
        'top_by_gene_ratio': rows(kg.sort_values('gene_ratio', ascending=False).head(20)),
        'als': {'rank': int(als['rank']), 'count': int(als['Count']), 'size': int(als['BgRatio'].split('/')[0]),
                'fold': round(float(als['FoldEnrichment']), 2), 'fdr': float(als['p.adjust'])},
        'showcase': [g for g in KEGG_SHOWCASE if g in genes],
    }


def clingen() -> dict:
    panels = pd.read_csv(D / 'clingen_panel_methylarginine.tsv', sep='\t')
    cur = pd.read_csv(D / 'clingen_curated_genes.tsv', sep='\t')
    genes = pd.read_csv(D / 'gene_table.tsv', sep='\t')
    cur['methyl_any'] = cur['methyl_any'].astype(str).eq('True')
    genes['methyl_any'] = genes['methyl_any'].astype(str).eq('True')
    als = panels[panels['panel'].str.startswith('Amyotrophic')].iloc[0]
    p = panels.sort_values('frac_methyl', ascending=False)
    return {
        'bars': [
            {'label': 'ALS-causing proteins', 'pct': round(100 * float(als['frac_methyl']), 1),
             'n': f"{int(als['n_methyl'])} of {int(als['n_set'])}"},
            {'label': 'All disease proteins', 'pct': round(100 * float(cur['methyl_any'].mean()), 1),
             'n': f"{int(cur['methyl_any'].sum()):,} of {len(cur):,}"},
            {'label': 'All human proteins', 'pct': round(100 * float(genes['methyl_any'].mean()), 1),
             'n': f"{int(genes['methyl_any'].sum()):,} of {len(genes):,}"},
        ],
        'panels': [{'panel': r.panel, 'pct': round(100 * float(r.frac_methyl), 1), 'n': int(r.n_set),
                    'als': r.panel.startswith('Amyotrophic')} for r in p.itertuples()],
        'als_adj_or': [round(float(als[c]), 2) for c in ('adj_OR', 'adj_OR_lo95', 'adj_OR_hi95')],
        'als_adj_p': float(als['adj_p']),
        'als_methyl_genes': als['methyl_genes'].split(','),
    }


def do() -> dict:
    do_ = pd.read_csv(D / 'dose_DO_reviewed_bg_full.tsv', sep='\t')
    show = pd.concat([do_.head(10), do_[do_['Description'].isin(['motor neuron disease', 'amyotrophic lateral sclerosis'])
                                         & (do_['rank'] > 10)]])
    mnd = do_[do_['Description'] == 'motor neuron disease'].iloc[0]
    det = pd.read_csv(D / 'dose_DO_ms_detectable_bg_full.tsv', sep='\t')
    mnd_det = det[det['Description'] == 'motor neuron disease'].iloc[0]
    return {
        'n_terms': len(do_), 'n_query': int(mnd['GeneRatio'].split('/')[1]),
        'rows': [{'name': r['Description'][0].upper() + r['Description'][1:], 'rank': int(r['rank']),
                  'count': int(r['GeneRatio'].split('/')[0]), 'size': int(r['BgRatio'].split('/')[0]),
                  'neglog10fdr': round(float(-np.log10(r['p.adjust'])), 2), 'fdr': float(r['p.adjust']),
                  'highlight': r['Description'] in ('motor neuron disease', 'amyotrophic lateral sclerosis')}
                 for r in show.to_dict('records')],
        'mnd': {'rank': int(mnd['rank']), 'count': int(mnd['Count']), 'fdr': float(mnd['p.adjust'])},
        'mnd_ms_detectable': {'rank': int(mnd_det['rank']), 'fdr': float(mnd_det['p.adjust'])},
        'showcase': [g for g in DO_SHOWCASE if g in set(mnd['geneSymbol'].split('/'))],
    }


def proteins() -> list[dict]:
    union = pd.read_csv(RESULTS / 'integrated' / 'human_arg_methyl_union_dedup_by_site.tsv', sep='\t', low_memory=False,
                        usecols=['canonical_UniProtAC', 'corrected_position', 'reference_any'])
    out = []
    for gene, acc, label in PROTEINS:
        d = json.load(open(DATA / 'context' / 'uniprot_features' / f'{acc}.json'))
        feats = []
        for f in d.get('features', []):
            s, e, desc = f['location']['start']['value'], f['location']['end']['value'], f.get('description', '')
            if f['type'] in ('Domain', 'Zinc finger') and desc != 'IQ':
                feats.append({'kind': 'domain', 'start': s, 'end': e,
                              'label': 'RRM' if desc.startswith('RRM') else ('ZnF' if f['type'] == 'Zinc finger' else desc)})
            elif (f['type'] == 'Compositional bias' and ('Gly' in desc or 'Low complexity' in desc)) or \
                    (f['type'] == 'Region' and ('RGG' in desc or 'D-R-' in desc)):
                feats.append({'kind': 'lcd', 'start': s, 'end': e, 'label': desc})
        sites = union[union['canonical_UniProtAC'] == acc]
        studies = {}
        for r in sites.itertuples():
            refs = [x.strip() for x in str(r.reference_any).split(';') if x.strip() and x.strip() != 'nan']
            studies[int(r.corrected_position)] = len(set(refs))
        out.append({'gene': gene, 'label': label, 'uniprot': acc, 'length': d['sequence']['length'],
                    'features': feats, 'sites': sorted(studies), 'studies_per_site': studies})
    return out


PRMT_PROTEINS = [('FUS', 'FUS'), ('EWSR1', 'EWSR1'), ('TAF15', 'TAF15'), ('HNRNPA1', 'hnRNPA1'),
                 ('HNRNPA2B1', 'hnRNPA2B1'), ('HNRNPA3', 'hnRNPA3'), ('G3BP1', 'G3BP1')]


def condensates() -> dict:
    """Maron 2021 Fig. 7 analogue: condensate / stress-granule membership (CD-CODE) of methylated ALS proteins."""
    g = pd.read_csv(D / 'gene_table.tsv', sep='\t')
    panels = pd.read_csv(D / 'clingen_panel_methylarginine.tsv', sep='\t')
    als_genes = panels[panels['panel'].str.startswith('Amyotrophic')].iloc[0]['methyl_genes'].split(',')
    mem = pd.read_csv(RESULTS / 'condensates' / 'condensate_membership.tsv', sep='\t')
    any_c = set(mem['canonical_UniProtAC'])
    sg = set(mem.loc[mem['condensate_name'].str.contains('stress granule', case=False), 'canonical_UniProtAC'])
    als = set(g.loc[g['symbol'].isin(als_genes), 'uniprot'])
    meth = set(g.loc[g['methyl_any'].astype(str).eq('True'), 'uniprot'])
    prot = set(g['uniprot'])
    groups = [('Methylated ALS-causing proteins', als), ('All methylated proteins', meth), ('All human proteins', prot)]
    return {
        'groups': [{'label': lab, 'n': len(S), 'any_pct': round(100 * len(S & any_c) / len(S), 1),
                    'sg_pct': round(100 * len(S & sg) / len(S), 1)} for lab, S in groups],
        'als_sg_proteins': sorted(g.loc[g['uniprot'].isin(als & sg), 'symbol']),
    }


def prmt() -> dict:
    """Maron 2021 Table S3: site-level log2(treated/DMSO) in MMA/ADMA/SDMA PTMScan after GSK591 or MS023."""
    t = pd.read_csv(RESULTS / 'maron_state_treatment' / 'maron_table_s3_site_treatment_collapsed.tsv', sep='\t')
    cols = [(tr, ip) for tr in ('GSK591', 'MS023') for ip in ('MMA', 'ADMA', 'SDMA')]
    rows = []
    for gene, label in PRMT_PROTEINS + [(None, 'All sites (proteome)')]:
        sub = t if gene is None else t[t['gene'] == gene]
        cells = []
        for tr, ip in cols:
            c = sub[(sub['treatment'] == tr) & (sub['ip'] == ip)]
            cells.append({'n': int(len(c)), 'median': None if c.empty else round(float(c['treatment_diff_median'].median()), 2),
                          'sig': int((c['treatment_padj_min'] < 0.05).sum())})
        rows.append({'label': label, 'cells': cells})
    return {'columns': [f'{tr} {ip}' for tr, ip in cols], 'rows': rows,
            'tdp43_detected': bool((t['gene'] == 'TARDBP').any())}


def main() -> None:
    data = {'kegg': kegg(), 'clingen': clingen(), 'do': do(), 'proteins': proteins(),
            'condensates': condensates(), 'prmt': prmt()}
    OUT.write_text(json.dumps(data, indent=1))
    print('wrote', OUT)
    print({p['label']: len(p['sites']) for p in data['proteins']}, 'TDP-43 studies', data['proteins'][3]['studies_per_site'])


if __name__ == '__main__':
    main()

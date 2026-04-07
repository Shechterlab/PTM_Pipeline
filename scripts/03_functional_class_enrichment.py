from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from common import canonicalize_uniprot_accession, classify_by_patterns, fisher_like_enrichment, load_json, save_table


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--base-master', required=True)
    ap.add_argument('--integrated-sites', required=True)
    ap.add_argument('--ontology', default='config/functional_ontology.json')
    ap.add_argument('--outdir', default='results/functional_class_union')
    args = ap.parse_args()

    master = pd.read_csv(args.base_master, sep='\t', low_memory=False)
    sites = pd.read_csv(args.integrated_sites, sep='\t', low_memory=False)
    ontology = load_json(args.ontology)

    master['canonical_UniProtAC'] = master['substrate_UniProtAC'].astype(str).map(canonicalize_uniprot_accession)
    ann_cols = [c for c in ['canonical_UniProtAC', 'substrate_UniProtAC', 'substrate_genename', 'UniProtID', 'protein_name', 'genename'] if c in master.columns]
    protein_meta = master[ann_cols].drop_duplicates('canonical_UniProtAC').copy()
    protein_meta['annotation_text'] = protein_meta.fillna('').astype(str).agg(' '.join, axis=1)
    protein_meta['functional_class_broad'] = protein_meta['annotation_text'].map(lambda x: classify_by_patterns(x, ontology['broad_classes']))
    protein_meta['functional_class_subclass'] = protein_meta['annotation_text'].map(lambda x: classify_by_patterns(x, ontology['subclasses']))

    target_column = 'canonical_UniProtAC' if 'canonical_UniProtAC' in sites.columns else 'substrate_UniProtAC'
    target_acc = set(sites[target_column].dropna().astype(str))
    protein_meta['is_arg_methyl_protein'] = protein_meta['canonical_UniProtAC'].isin(target_acc)
    target = protein_meta[protein_meta['is_arg_methyl_protein']].copy()

    broad = fisher_like_enrichment(target['functional_class_broad'], protein_meta['functional_class_broad'])
    sub = fisher_like_enrichment(target['functional_class_subclass'], protein_meta['functional_class_subclass'])

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(protein_meta, outdir / 'protein_functional_classes.tsv')
    save_table(broad, outdir / 'arg_methyl_functional_class_broad_enrichment.tsv')
    save_table(sub, outdir / 'arg_methyl_functional_class_subclass_enrichment.tsv')

    pb = broad[broad['target_count'] >= 10].head(10).sort_values('odds_ratio')
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    ax.barh(pb['category'], pb['odds_ratio'])
    ax.set_xlabel('Odds ratio vs all PTM proteins')
    ax.set_ylabel('Broad class')
    ax.set_title('Integrated arg-methylome: broad functional-class enrichment')
    for y, v, n in zip(pb['category'], pb['odds_ratio'], pb['target_count']):
        ax.text(v, y, f'  n={n}, OR={v:.1f}', va='center', fontsize=9)
    fig.tight_layout()
    fig.savefig(outdir / 'arg_methyl_functional_class_broad_enrichment.png', dpi=300)
    fig.savefig(outdir / 'arg_methyl_functional_class_broad_enrichment.svg')

    ps = sub[sub['target_count'] >= 5].head(12).sort_values('odds_ratio')
    fig2, ax2 = plt.subplots(figsize=(8.9, 6.3))
    ax2.barh(ps['category'], ps['odds_ratio'])
    ax2.set_xlabel('Odds ratio vs all PTM proteins')
    ax2.set_ylabel('Subclass')
    ax2.set_title('Integrated arg-methylome: review-grade subclass enrichment')
    for y, v, n in zip(ps['category'], ps['odds_ratio'], ps['target_count']):
        ax2.text(v, y, f'  n={n}, OR={v:.1f}', va='center', fontsize=9)
    fig2.tight_layout()
    fig2.savefig(outdir / 'arg_methyl_functional_class_subclass_enrichment.png', dpi=300)
    fig2.savefig(outdir / 'arg_methyl_functional_class_subclass_enrichment.svg')


if __name__ == '__main__':
    main()

from __future__ import annotations

import argparse
from pathlib import Path
import re

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import fisher_exact

from common import (
    add_q_values,
    canonicalize_uniprot_accession,
    format_p_value,
    load_json,
    log2_odds_ratio,
    save_figure,
    save_table,
)


def build_annotation_text(master: pd.DataFrame) -> pd.DataFrame:
    ann_cols = [
        c for c in [
            'canonical_UniProtAC',
            'substrate_UniProtAC',
            'substrate_genename',
            'UniProtID',
            'protein_name',
            'genename',
        ] if c in master.columns
    ]
    protein_meta = master[ann_cols].drop_duplicates('canonical_UniProtAC').copy()
    protein_meta['annotation_text'] = protein_meta.fillna('').astype(str).agg(' '.join, axis=1).str.lower()
    return protein_meta


def compile_label_defs(label_defs: list[dict]) -> list[dict]:
    compiled = []
    for label_def in label_defs:
        compiled.append({
            'label': label_def['label'],
            'patterns': [re.compile(pattern, flags=re.IGNORECASE) for pattern in label_def.get('patterns', [])],
        })
    return compiled


def labels_for_text(text: str, label_defs: list[dict]) -> list[str]:
    out = []
    for label_def in label_defs:
        patterns = label_def.get('patterns', [])
        if any(pattern.search(text) for pattern in patterns):
            out.append(label_def['label'])
    return out


def expand_membership(protein_meta: pd.DataFrame, label_defs: list[dict], family_name: str) -> pd.DataFrame:
    rows = []
    for row in protein_meta.itertuples(index=False):
        labels = labels_for_text(str(row.annotation_text), label_defs)
        if not labels:
            labels = ['Other / unclassified']
        for label in labels:
            rows.append({
                'canonical_UniProtAC': row.canonical_UniProtAC,
                'family': family_name,
                'label': label,
            })
    return pd.DataFrame(rows).drop_duplicates()


def odds_ratio(a: int, b: int, c: int, d: int) -> float:
    return ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))


def label_enrichment(memberships: pd.DataFrame, target_acc: set[str], family_name: str) -> pd.DataFrame:
    family = memberships[memberships['family'] == family_name].copy()
    universe = set(family['canonical_UniProtAC'].astype(str))
    background_acc = universe - target_acc
    rows = []
    for label, group in family.groupby('label'):
        label_acc = set(group['canonical_UniProtAC'].astype(str))
        a = len(label_acc & target_acc)
        b = len(target_acc - label_acc)
        c = len(label_acc & background_acc)
        d = len(background_acc - label_acc)
        rows.append({
            'category': label,
            'target_count': a,
            'target_fraction': a / len(target_acc) if target_acc else pd.NA,
            'background_count': c,
            'background_fraction': c / len(background_acc) if background_acc else pd.NA,
            'odds_ratio': odds_ratio(a, b, c, d),
            'log2_odds_ratio': log2_odds_ratio(odds_ratio(a, b, c, d)),
            'p_value': float(fisher_exact([[a, b], [c, d]], alternative='two-sided').pvalue),
            'family': family_name,
        })
    out = add_q_values(pd.DataFrame(rows))
    return out.sort_values(['odds_ratio', 'target_count'], ascending=[False, False]).reset_index(drop=True)


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
    protein_meta = build_annotation_text(master)
    target_column = 'canonical_UniProtAC' if 'canonical_UniProtAC' in sites.columns else 'substrate_UniProtAC'
    target_acc = set(sites[target_column].dropna().astype(str).map(canonicalize_uniprot_accession))

    broad_membership = expand_membership(protein_meta, compile_label_defs(ontology['broad_labels']), 'broad')
    subclass_membership = expand_membership(protein_meta, compile_label_defs(ontology['subclass_labels']), 'subclass')
    memberships = pd.concat([broad_membership, subclass_membership], ignore_index=True)

    protein_summary = protein_meta[['canonical_UniProtAC']].copy()
    broad_map = broad_membership.groupby('canonical_UniProtAC')['label'].apply(lambda s: ';'.join(sorted(set(s)))).rename('functional_class_broad')
    sub_map = subclass_membership.groupby('canonical_UniProtAC')['label'].apply(lambda s: ';'.join(sorted(set(s)))).rename('functional_class_subclass')
    protein_summary = protein_summary.merge(broad_map, how='left', on='canonical_UniProtAC')
    protein_summary = protein_summary.merge(sub_map, how='left', on='canonical_UniProtAC')
    protein_summary = protein_summary.merge(protein_meta, how='left', on='canonical_UniProtAC')
    protein_summary['is_arg_methyl_protein'] = protein_summary['canonical_UniProtAC'].isin(target_acc)

    broad = label_enrichment(memberships, target_acc, 'broad')
    sub = label_enrichment(memberships, target_acc, 'subclass')

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(protein_summary, outdir / 'protein_functional_labels.tsv')
    save_table(broad_membership, outdir / 'protein_functional_broad_membership.tsv')
    save_table(subclass_membership, outdir / 'protein_functional_subclass_membership.tsv')
    save_table(broad, outdir / 'arg_methyl_functional_class_broad_enrichment.tsv')
    save_table(sub, outdir / 'arg_methyl_functional_class_subclass_enrichment.tsv')

    pb = broad[broad['target_count'] >= 15].head(12).sort_values('odds_ratio')
    fig, ax = plt.subplots(figsize=(9.0, 6.2))
    ax.barh(pb['category'], pb['log2_odds_ratio'])
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlabel('log2(OR) vs all proteins')
    ax.set_ylabel('Broad class')
    ax.set_title('Integrated arg-methylome: multi-label broad enrichment')
    for y, v, n, p, q in zip(pb['category'], pb['log2_odds_ratio'], pb['target_count'], pb['p_value'], pb['q_value']):
        ax.text(v, y, f'  n={n}, p={format_p_value(p)}, q={format_p_value(q)}', va='center', ha='left' if v >= 0 else 'right', fontsize=8.5)
    fig.tight_layout()
    save_figure(fig, outdir / 'arg_methyl_functional_class_broad_enrichment')

    ps = sub[sub['target_count'] >= 8].head(14).sort_values('odds_ratio')
    fig2, ax2 = plt.subplots(figsize=(9.8, 7.2))
    ax2.barh(ps['category'], ps['log2_odds_ratio'])
    ax2.axvline(0, color='black', linewidth=0.8)
    ax2.set_xlabel('log2(OR) vs all proteins')
    ax2.set_ylabel('Subclass')
    ax2.set_title('Integrated arg-methylome: multi-label RNA-centric subclasses')
    for y, v, n, p, q in zip(ps['category'], ps['log2_odds_ratio'], ps['target_count'], ps['p_value'], ps['q_value']):
        ax2.text(v, y, f'  n={n}, p={format_p_value(p)}, q={format_p_value(q)}', va='center', ha='left' if v >= 0 else 'right', fontsize=8.2)
    fig2.tight_layout()
    save_figure(fig2, outdir / 'arg_methyl_functional_class_subclass_enrichment')


if __name__ == '__main__':
    main()

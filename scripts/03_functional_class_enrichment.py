from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError
from scipy.stats import fisher_exact

from common import (
    add_bh_q_values,
    canonicalize_uniprot_accession,
    format_p_value,
    load_json,
    log2_odds_ratio,
    save_figure,
    save_table,
)





def load_interpro_terms(path: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep='\t', low_memory=False)
    except (EmptyDataError, FileNotFoundError):
        return pd.DataFrame(columns=['canonical_UniProtAC', 'interpro_annotation_text'])
    if df.empty or 'canonical_UniProtAC' not in df.columns:
        return pd.DataFrame(columns=['canonical_UniProtAC', 'interpro_annotation_text'])
    text_cols = [col for col in ['interpro_name', 'interpro_type', 'member_database_accessions', 'interpro_accession'] if col in df.columns]
    if not text_cols:
        return pd.DataFrame(columns=['canonical_UniProtAC', 'interpro_annotation_text'])
    merged = df[['canonical_UniProtAC'] + text_cols].fillna('').astype(str)
    merged['interpro_annotation_text'] = merged[text_cols].agg(' '.join, axis=1)
    merged = merged.groupby('canonical_UniProtAC', as_index=False)['interpro_annotation_text'].agg(' '.join)
    return merged

DEFAULT_MULTI_LABEL_CATEGORIES = [
    {'label': 'General RNA-binding / RNP', 'panel': 'rna_centric', 'patterns': [r'rna binding', r'rna-binding', r'ribonucle', r'\brbm\d*\b', r'\bhnrnp', r'\brnp\b', r'\brbmxl', r'\bpcbp\d*\b']},
    {'label': 'RRM / KH RNA-binding proteins', 'panel': 'rna_centric', 'patterns': [r'rna recognition motif', r'\brrm\b', r'kh domain', r'k homology', r'\bkh\b', r'quaking', r'igf2bp', r'khsrp', r'\bqki\b']},
    {'label': 'hnRNP / RG-rich RBPs', 'panel': 'rna_centric', 'patterns': [r'hnrnp', r'heterogeneous nuclear ribonucleoprotein', r'rg-rich', r'rgg', r'rg motif', r'fused in sarcoma', r'\bfus\b', r'ewsr1', r'taf15', r'\brbm\d+\b', r'\bprrc2']},
    {'label': 'Spliceosome / splicing factors', 'panel': 'rna_centric', 'patterns': [r'splice', r'snrnp', r'prpf', r'sf3', r'u2af', r'serine/arginine-rich', r'\blsm\b', r'sm protein', r'\bsnrp', r'\bsrsf\d+\b', r'\brbmx\b']},
    {'label': 'RNA helicases / remodeling ATPases', 'panel': 'rna_centric', 'patterns': [r'rna helicase', r'helicase', r'dead-box', r'\bddx\d+\b', r'\bdhx\d+\b', r'deah', r'dexd', r'ski2']},
    {'label': 'mRNA processing / export / 3-prime end', 'panel': 'rna_centric', 'patterns': [r'mrna processing', r'mrna export', r'cleavage and polyadenylation', r'polyadenyl', r'pcf11', r'thoc', r'alyref', r'exportin', r'nxf1', r'cpsf', r'cstf', r'cfim', r'cleavage factor', r'3-prime end']},
    {'label': 'RNA surveillance / decay / exosome', 'panel': 'rna_centric', 'patterns': [r'exosome', r'nonsense-mediated decay', r'upf\d', r'ski complex', r'rna decay', r'deadenyl', r'decapping', r'dcp\d', r'xrn\d', r'lsm1-7']},
    {'label': 'Ribosome biogenesis / nucleolar RNP', 'panel': 'rna_centric', 'patterns': [r'ribosome biogenesis', r'nucleolar', r'\bsno', r'rrna processing', r'ribosomal rna', r'fibrillarin', r'\bnop\d*\b', r'dyskerin', r'nucleophosmin', r'gar1', r'ncl\b']},
    {'label': 'Translation / ribosome / initiation factors', 'panel': 'rna_centric', 'patterns': [r'ribosomal protein', r'translation', r'elongation factor', r'initiation factor', r'\beif\d', r'\beef\d', r'trna synthetase', r'rps\d', r'rpl\d']},
    {'label': 'Stress granule / paraspeckle / P-body factors', 'panel': 'rna_centric', 'patterns': [r'stress granule', r'paraspeck', r'p-body', r'processing body', r'\bnono\b', r'\bsfpq\b', r'\brbm14\b', r'tardbp', r'fmr1', r'g3bp\d', r'tia1', r'tiar', r'phase separation', r'condensate']},
    {'label': 'Chromatin / transcription regulators', 'panel': 'non_rna', 'patterns': [r'chromatin', r'histone', r'nucleosome', r'transcription', r'polymerase', r'mediator', r'epigen', r'bromodomain', r'pwwp', r'tudor', r'reader', r'writer', r'remodel', r'mll', r'kmt2', r'cbx\d']},
    {'label': 'DNA repair / genome maintenance', 'panel': 'non_rna', 'patterns': [r'dna repair', r'damage response', r'double-strand', r'homologous recombination', r'replication', r'fork', r'brca', r'53bp', r'rad51', r'mre11', r'atm', r'atr', r'fanconi', r'genome maintenance']},
    {'label': 'Ubiquitin / proteostasis', 'panel': 'non_rna', 'patterns': [r'ubiquitin', r'proteasome', r'sumo', r'autophagy', r'chaperone', r'heat shock', r'protein quality control', r'ubiquitin ligase', r'deubiquitin']},
    {'label': 'Cytoskeleton / trafficking / membrane', 'panel': 'non_rna', 'patterns': [r'actin', r'microtub', r'myosin', r'dynein', r'kinesin', r'vesicle', r'membrane', r'golgi', r'endosome', r'mitochond', r'coatomer', r'clathrin']},
    {'label': 'Signaling / catalytic enzymes', 'panel': 'non_rna', 'patterns': [r'kinase', r'phosphatase', r'gtpase', r'atpase', r'metabolic', r'enzyme', r'cataly', r'transferase', r'hydrolase', r'ligase']},
    {'label': 'Organelle / metabolism / other enzymes', 'panel': 'non_rna', 'patterns': [r'mitochond', r'peroxis', r'lysosom', r'golgi', r'endoplasmic reticulum', r'metabolism', r'oxidoreductase', r'synthetase', r'dehydrogenase']},
]


def normalize_categories(ontology: dict) -> list[dict]:
    blocks = ontology.get('multi_label_categories')
    if isinstance(blocks, list) and blocks:
        return blocks
    return DEFAULT_MULTI_LABEL_CATEGORIES


def match_any(text: str, patterns: list[str]) -> bool:
    lowered = str(text).lower()
    return any(re.search(pattern, lowered) for pattern in patterns)


def compute_binary_enrichment(target_flags: pd.Series, universe_flags: pd.Series, label: str, target_total: int, universe_total: int) -> dict:
    a = int(target_flags.sum())
    b = int(target_total - a)
    inclusive = int(universe_flags.sum())
    c = int(inclusive - a)
    d = int((universe_total - target_total) - c)
    orr = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    return {
        'category': label,
        'target_count': a,
        'target_fraction': a / target_total if target_total else np.nan,
        'background_count': c,
        'background_fraction': c / (universe_total - target_total) if universe_total > target_total else np.nan,
        'inclusive_universe_count': inclusive,
        'inclusive_universe_fraction': inclusive / universe_total if universe_total else np.nan,
        'odds_ratio': orr,
        'log2_odds_ratio': log2_odds_ratio(orr),
        'p_value': float(fisher_exact([[a, b], [c, d]], alternative='two-sided').pvalue),
    }


def classify_proteins(protein_meta: pd.DataFrame, categories: list[dict]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    columns = {}
    for block in categories:
        label = block['label']
        mask = protein_meta['annotation_text'].map(lambda text: match_any(text, block['patterns']))
        columns[label] = mask
        matched = protein_meta.loc[mask, ['canonical_UniProtAC']].copy()
        matched['functional_category'] = label
        matched['panel'] = block.get('panel', 'all')
        rows.append(matched)
    membership = protein_meta[['canonical_UniProtAC']].copy()
    for label, mask in columns.items():
        membership[label] = mask.to_numpy()
    if columns:
        assigned_any = np.column_stack([mask.to_numpy(dtype=bool) for mask in columns.values()]).any(axis=1)
        rna_masks = [columns[block['label']].to_numpy(dtype=bool) for block in categories if block.get('panel') == 'rna_centric' and block['label'] in columns]
        assigned_rna = np.column_stack(rna_masks).any(axis=1) if rna_masks else np.zeros(len(protein_meta), dtype=bool)
    else:
        assigned_any = np.zeros(len(protein_meta), dtype=bool)
        assigned_rna = np.zeros(len(protein_meta), dtype=bool)
    membership['Other / unclassified'] = ~assigned_any
    membership['Any RNA-centric label'] = assigned_rna
    rows.append(protein_meta.loc[~assigned_any, ['canonical_UniProtAC']].assign(functional_category='Other / unclassified', panel='other'))
    rows.append(protein_meta.loc[assigned_rna, ['canonical_UniProtAC']].assign(functional_category='Any RNA-centric label', panel='summary'))
    long_df = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=['canonical_UniProtAC', 'functional_category', 'panel'])
    return membership, long_df


def plot_panel(df: pd.DataFrame, title: str, xlabel: str, outpath: Path, min_count: int = 15, top_n: int = 14) -> None:
    if df.empty:
        return
    plot_df = df[df['target_count'] >= min_count].sort_values(['q_value_bh', 'odds_ratio', 'target_count'], ascending=[True, False, False]).head(top_n)
    if plot_df.empty:
        return
    plot_df = plot_df.sort_values('log2_odds_ratio')
    fig, ax = plt.subplots(figsize=(10.2, max(5.2, 0.42 * len(plot_df) + 1.8)))
    ax.barh(plot_df['category'], plot_df['log2_odds_ratio'])
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Functional category')
    ax.set_title(title)
    for y, v, n, q in zip(plot_df['category'], plot_df['log2_odds_ratio'], plot_df['target_count'], plot_df['q_value_bh']):
        ax.text(v, y, f'  n={n}, q={format_p_value(q)}', va='center', ha='left' if v >= 0 else 'right', fontsize=9)
    fig.tight_layout()
    save_figure(fig, outpath)


def plot_membership_distribution(protein_meta: pd.DataFrame, outpath: Path) -> None:
    if protein_meta.empty:
        return
    counts = protein_meta['functional_category_count'].value_counts().sort_index().rename_axis('label_count').reset_index(name='protein_count')
    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    ax.bar(counts['label_count'].astype(str), counts['protein_count'])
    ax.set_xlabel('Number of assigned multi-label categories per protein')
    ax.set_ylabel('Proteins')
    ax.set_title('Multi-label annotation coverage across PTM proteins')
    fig.tight_layout()
    save_figure(fig, outpath)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--base-master', required=True)
    ap.add_argument('--integrated-sites', required=True)
    ap.add_argument('--ontology', default='config/functional_ontology.json')
    ap.add_argument('--interpro-intervals', default='')
    ap.add_argument('--outdir', default='results/functional_class_union')
    args = ap.parse_args()

    master = pd.read_csv(args.base_master, sep='\t', low_memory=False)
    sites = pd.read_csv(args.integrated_sites, sep='\t', low_memory=False)
    ontology = load_json(args.ontology)
    categories = normalize_categories(ontology)

    master['canonical_UniProtAC'] = master['substrate_UniProtAC'].astype(str).map(canonicalize_uniprot_accession)
    ann_cols = [c for c in ['canonical_UniProtAC', 'substrate_UniProtAC', 'substrate_genename', 'UniProtID', 'protein_name', 'genename'] if c in master.columns]
    protein_meta = master[ann_cols].drop_duplicates('canonical_UniProtAC').copy()
    protein_meta['annotation_text'] = protein_meta.fillna('').astype(str).agg(' '.join, axis=1)
    if args.interpro_intervals:
        interpro_terms = load_interpro_terms(args.interpro_intervals)
        protein_meta = protein_meta.merge(interpro_terms, how='left', on='canonical_UniProtAC')
        protein_meta['annotation_text'] = (protein_meta['annotation_text'].fillna('') + ' ' + protein_meta['interpro_annotation_text'].fillna('')).str.strip()

    target_column = 'canonical_UniProtAC' if 'canonical_UniProtAC' in sites.columns else 'substrate_UniProtAC'
    target_acc = set(sites[target_column].dropna().astype(str))
    protein_meta['is_arg_methyl_protein'] = protein_meta['canonical_UniProtAC'].isin(target_acc)

    membership, long_membership = classify_proteins(protein_meta, categories)
    protein_meta = protein_meta.merge(membership, how='left', on='canonical_UniProtAC')
    protein_meta['Other / unclassified'] = protein_meta['Other / unclassified'].fillna(True)
    category_cols = [col for col in membership.columns if col != 'canonical_UniProtAC']
    protein_meta['functional_category_count'] = protein_meta[category_cols].fillna(False).sum(axis=1).to_numpy(dtype=int)

    target = protein_meta[protein_meta['is_arg_methyl_protein']].copy()
    universe_total = len(protein_meta)
    target_total = len(target)
<<<<<<< HEAD

    enrichment_rows = []
    panel_map = {block['label']: block.get('panel', 'all') for block in categories}
    panel_map['Other / unclassified'] = 'other'
    panel_map['Any RNA-centric label'] = 'summary'
    for category in category_cols:
        enrichment = compute_binary_enrichment(target[category].fillna(False).astype(bool), protein_meta[category].fillna(False).astype(bool), category, target_total, universe_total)
        enrichment['panel'] = panel_map.get(category, 'all')
        enrichment_rows.append(enrichment)
    enrichment_df = add_bh_q_values(pd.DataFrame(enrichment_rows), p_col='p_value', q_col='q_value_bh')
    enrichment_df = enrichment_df.sort_values(['panel', 'q_value_bh', 'odds_ratio', 'target_count'], ascending=[True, True, False, False])

=======

    enrichment_rows = []
    panel_map = {block['label']: block.get('panel', 'all') for block in categories}
    panel_map['Other / unclassified'] = 'other'
    panel_map['Any RNA-centric label'] = 'summary'
    for category in category_cols:
        enrichment = compute_binary_enrichment(target[category].fillna(False).astype(bool), protein_meta[category].fillna(False).astype(bool), category, target_total, universe_total)
        enrichment['panel'] = panel_map.get(category, 'all')
        enrichment_rows.append(enrichment)
    enrichment_df = add_bh_q_values(pd.DataFrame(enrichment_rows), p_col='p_value', q_col='q_value_bh')
    enrichment_df = enrichment_df.sort_values(['panel', 'q_value_bh', 'odds_ratio', 'target_count'], ascending=[True, True, False, False])

>>>>>>> 4fbe8b2dc495fcf875e4e868ffb2aca4033251f1
    counts = long_membership.merge(protein_meta[['canonical_UniProtAC', 'is_arg_methyl_protein']], how='left', on='canonical_UniProtAC')
    counts_summary = counts.groupby(['functional_category', 'panel', 'is_arg_methyl_protein'], as_index=False).size().rename(columns={'size': 'protein_count'})
    coverage = pd.DataFrame([
        {'metric': 'universe_proteins', 'value': int(universe_total)},
        {'metric': 'arg_methyl_union_proteins', 'value': int(target_total)},
        {'metric': 'arg_methyl_proteins_with_any_rna_label', 'value': int(target['Any RNA-centric label'].fillna(False).sum())},
        {'metric': 'arg_methyl_proteins_unclassified', 'value': int(target['Other / unclassified'].fillna(False).sum())},
        {'metric': 'all_ptm_proteins_with_any_rna_label', 'value': int(protein_meta['Any RNA-centric label'].fillna(False).sum())},
        {'metric': 'all_ptm_proteins_unclassified', 'value': int(protein_meta['Other / unclassified'].fillna(False).sum())},
        {'metric': 'median_labels_per_arg_methyl_protein', 'value': float(target['functional_category_count'].median()) if not target.empty else np.nan},
        {'metric': 'median_labels_per_all_ptm_protein', 'value': float(protein_meta['functional_category_count'].median()) if not protein_meta.empty else np.nan},
    ])

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(protein_meta, outdir / 'protein_functional_multilabel_annotations.tsv')
    save_table(long_membership, outdir / 'protein_functional_multilabel_long.tsv')
    save_table(enrichment_df, outdir / 'arg_methyl_functional_multilabel_enrichment.tsv')
    save_table(counts_summary, outdir / 'functional_multilabel_category_counts.tsv')
    save_table(coverage, outdir / 'functional_multilabel_annotation_coverage.tsv')

<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 4fbe8b2dc495fcf875e4e868ffb2aca4033251f1
    plot_panel(
        enrichment_df[enrichment_df['panel'].isin(['rna_centric', 'summary'])].copy(),
        title='Integrated arg-methylome: RNA-centric multi-label enrichment',
        xlabel='log2(OR) vs all PTM proteins',
        outpath=outdir / 'arg_methyl_functional_multilabel_rna_centric_enrichment',
    )
    plot_panel(
        enrichment_df[enrichment_df['panel'].isin(['non_rna', 'other'])].copy(),
        title='Integrated arg-methylome: non-RNA multi-label enrichment',
        xlabel='log2(OR) vs all PTM proteins',
        outpath=outdir / 'arg_methyl_functional_multilabel_non_rna_enrichment',
    )
    plot_membership_distribution(protein_meta, outdir / 'functional_multilabel_membership_distribution')
<<<<<<< HEAD
=======
    pb = broad[broad['target_count'] >= 10].head(10).sort_values('odds_ratio')
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    ax.barh(pb['category'], pb['log2_odds_ratio'])
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlabel('log2(OR) vs all PTM proteins')
    ax.set_ylabel('Broad class')
    ax.set_title('Integrated arg-methylome: broad functional-class enrichment')
    for y, v, n, q in zip(pb['category'], pb['log2_odds_ratio'], pb['target_count'], pb['q_value_bh']):
        ax.text(v, y, f'  n={n}, q={format_p_value(q)}', va='center', ha='left' if v >= 0 else 'right', fontsize=9)
    fig.tight_layout()
    save_figure(fig, outdir / 'arg_methyl_functional_class_broad_enrichment')

    ps = sub[sub['target_count'] >= 5].head(12).sort_values('odds_ratio')
    fig2, ax2 = plt.subplots(figsize=(8.9, 6.3))
    ax2.barh(ps['category'], ps['log2_odds_ratio'])
    ax2.axvline(0, color='black', linewidth=0.8)
    ax2.set_xlabel('log2(OR) vs all PTM proteins')
    ax2.set_ylabel('Subclass')
    ax2.set_title('Integrated arg-methylome: review-grade subclass enrichment')
    for y, v, n, q in zip(ps['category'], ps['log2_odds_ratio'], ps['target_count'], ps['q_value_bh']):
        ax2.text(v, y, f'  n={n}, q={format_p_value(q)}', va='center', ha='left' if v >= 0 else 'right', fontsize=9)
    fig2.tight_layout()
    save_figure(fig2, outdir / 'arg_methyl_functional_class_subclass_enrichment')
>>>>>>> 2c1a8ff7d8603628918f8ffc773cdcc98c903bff
=======
>>>>>>> 4fbe8b2dc495fcf875e4e868ffb2aca4033251f1


if __name__ == '__main__':
    main()

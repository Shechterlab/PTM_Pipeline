from __future__ import annotations

import argparse
from pathlib import Path
import re

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import fisher_exact

from common import (
    BREWER_COLORS,
    add_q_values,
    annotate_barh,
    apply_paper_style,
    balanced_category_subset,
    canonicalize_uniprot_accession,
    format_p_value,
    load_json,
    log2_odds_ratio,
    save_figure,
    save_table,
    set_symmetric_xlim,
    signed_bar_colors,
    style_axis,
)

SUBCLASS_TO_BROAD = {
    'RRM / KH RBPs': 'RNA-binding / RNP',
    'RNA-associated zinc finger / post-transcriptional': 'RNA-binding / RNP',
    'hnRNP / RG-rich RBPs': 'RNA-binding / RNP',
    'Spliceosome / splicing': 'RNA-binding / RNP',
    'RNA helicases': 'RNA-binding / RNP',
    "mRNA processing / export / 3'-end": 'RNA-binding / RNP',
    'Stress granule / paraspeckle / P-body': 'RNA-binding / RNP',
    'Poly(A) / stability RBPs': 'RNA-binding / RNP',
    'Nucleolar / ribosome biogenesis': 'Translation / ribosome / biogenesis',
    'Translation / initiation / elongation': 'Translation / ribosome / biogenesis',
    'Chromatin readers': 'Chromatin / transcription',
    'Chromatin writers / remodelers': 'Chromatin / transcription',
    'Transcription / mediator / polymerase': 'Chromatin / transcription',
    'DNA-binding transcription factors / developmental regulators': 'Chromatin / transcription',
    'DNA repair / DDR': 'DNA repair / genome maintenance',
    'Nuclear scaffold / pore / lamina': 'Nuclear scaffold / pore / lamina',
    'Cytoskeleton / trafficking': 'Cytoplasmic / trafficking / cytoskeleton',
    'Endomembrane / coat trafficking': 'Cytoplasmic / trafficking / cytoskeleton',
    'Junction / polarity / scaffold': 'Cytoplasmic / trafficking / cytoskeleton',
    'Adhesion / intermediate filament scaffolds': 'Cytoplasmic / trafficking / cytoskeleton',
    'Membrane / receptor / extracellular': 'Membrane / receptor / extracellular',
    'Mitochondrial / metabolic enzymes': 'Mitochondrial / metabolism',
    'Cell cycle / centrosome / spindle': 'Cell cycle / centrosome / spindle',
    'Proteostasis / chaperone / proteasome': 'Proteostasis / ubiquitin / chaperone',
    'Poorly characterized / specialized': 'Poorly characterized / specialized',
    'Signaling / enzymes': 'Signaling / enzymes',
}


def extract_gene_tokens(values: list[str]) -> list[str]:
    tokens: list[str] = []
    seen: set[str] = set()
    reserved = {'NAN', 'NONE', 'NAME', 'ORFNAMES', 'SYNONYMS', 'GENE', 'PROTEIN'}
    for value in values:
        text = str(value or '').strip()
        if text.lower() in {'', 'nan', 'none'}:
            continue
        for part in re.split(r'[;\s]+', text):
            part = part.strip().upper().strip(',')
            if not part or part in reserved or part in seen:
                continue
            if not re.fullmatch(r'[A-Z0-9][A-Z0-9-]{1,14}', part):
                continue
            tokens.append(part)
            seen.add(part)
    return tokens


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


def aggregate_unique_text(series: pd.Series) -> str:
    values = [str(v).strip() for v in series if str(v).strip() not in {'', 'nan', 'None'}]
    if not values:
        return ''
    return ' ; '.join(dict.fromkeys(values))


def split_semicolon_values(series: pd.Series) -> str:
    values: list[str] = []
    seen: set[str] = set()
    for raw in series:
        text = str(raw).strip()
        if text in {'', 'nan', 'None'}:
            continue
        for piece in text.split(';'):
            piece = piece.strip()
            if piece and piece not in seen:
                values.append(piece)
                seen.add(piece)
    return '; '.join(values)


def build_interpro_protein_metadata(interpro: pd.DataFrame, domain_ontology: dict | None = None) -> pd.DataFrame:
    if interpro is None or interpro.empty:
        return pd.DataFrame(columns=['canonical_UniProtAC', 'interpro_annotation_text', 'interpro_domain_classes'])
    out = interpro.copy()
    out['canonical_UniProtAC'] = out['canonical_UniProtAC'].astype(str).map(canonicalize_uniprot_accession)
    out = out[out['canonical_UniProtAC'].notna() & out['canonical_UniProtAC'].ne('')].copy()

    class_defs = compile_label_defs(domain_ontology['domain_classes']) if domain_ontology else []
    class_text = out[['interpro_name', 'interpro_type']].fillna('').astype(str).agg(' '.join, axis=1)
    out['interpro_domain_class'] = class_text.map(lambda x: '; '.join(labels_for_text(x, class_defs)) if class_defs else '')

    protein = out.groupby('canonical_UniProtAC', as_index=False).agg({
        'interpro_name': aggregate_unique_text,
        'interpro_accession': aggregate_unique_text,
        'interpro_type': aggregate_unique_text,
        'interpro_domain_class': split_semicolon_values,
    })
    informative_classes = protein['interpro_domain_class'].fillna('').str.split(';').map(
        lambda xs: '; '.join(
            x.strip()
            for x in xs
            if x.strip() and x.strip() not in {'Other domain', 'Repeat / superfamily'}
        )
    )
    protein['interpro_domain_classes'] = informative_classes
    protein['interpro_annotation_text'] = protein[
        ['interpro_name', 'interpro_accession', 'interpro_type', 'interpro_domain_classes']
    ].fillna('').astype(str).agg(' '.join, axis=1).str.strip()
    return protein[['canonical_UniProtAC', 'interpro_annotation_text', 'interpro_domain_classes']]


def build_combined_protein_metadata(
    master: pd.DataFrame,
    sites: pd.DataFrame,
    metadata: pd.DataFrame | None = None,
    interpro: pd.DataFrame | None = None,
    domain_ontology: dict | None = None,
) -> pd.DataFrame:
    frames = []

    master_cols = [c for c in ['canonical_UniProtAC', 'substrate_genename', 'UniProtID', 'protein_name', 'genename'] if c in master.columns]
    if master_cols:
        frames.append(master[master_cols].copy())

    site_cols = [c for c in ['canonical_UniProtAC', 'substrate_genename', 'UniProtID', 'protein_name'] if c in sites.columns]
    if site_cols:
        frames.append(sites[site_cols].copy())

    if metadata is not None and not metadata.empty:
        meta = metadata.copy()
        rename_map = {
            'Entry': 'canonical_UniProtAC',
            'Entry Name': 'UniProtID',
            'Gene Names (primary)': 'substrate_genename',
            'Protein names': 'protein_name',
        }
        meta = meta.rename(columns=rename_map)
        meta_cols = [c for c in ['canonical_UniProtAC', 'substrate_genename', 'UniProtID', 'protein_name'] if c in meta.columns]
        if meta_cols:
            frames.append(meta[meta_cols].copy())

    protein_meta = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=['canonical_UniProtAC'])
    protein_meta['canonical_UniProtAC'] = protein_meta['canonical_UniProtAC'].astype(str).map(canonicalize_uniprot_accession)
    protein_meta = protein_meta[protein_meta['canonical_UniProtAC'].notna() & (protein_meta['canonical_UniProtAC'] != '')].copy()
    agg_map = {
        col: aggregate_unique_text
        for col in ['substrate_genename', 'UniProtID', 'protein_name', 'genename']
        if col in protein_meta.columns
    }
    protein_meta = protein_meta.groupby('canonical_UniProtAC', as_index=False).agg(agg_map)
    for col in ['substrate_genename', 'UniProtID', 'protein_name', 'genename']:
        if col not in protein_meta.columns:
            protein_meta[col] = ''
    if interpro is not None and not interpro.empty:
        protein_meta = protein_meta.merge(
            build_interpro_protein_metadata(interpro, domain_ontology=domain_ontology),
            how='left',
            on='canonical_UniProtAC',
        )
    else:
        protein_meta['interpro_annotation_text'] = ''
        protein_meta['interpro_domain_classes'] = ''
    protein_meta['annotation_text'] = protein_meta[
        ['canonical_UniProtAC', 'substrate_genename', 'UniProtID', 'protein_name', 'genename', 'interpro_annotation_text', 'interpro_domain_classes']
    ].fillna('').astype(str).agg(' '.join, axis=1).str.lower()
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


def add_broad_fallback_from_subclass(broad_membership: pd.DataFrame, subclass_membership: pd.DataFrame) -> pd.DataFrame:
    fallback = subclass_membership.copy()
    fallback['label'] = fallback['label'].map(SUBCLASS_TO_BROAD)
    fallback = fallback[fallback['label'].notna()].copy()
    fallback['family'] = 'broad'
    combined = pd.concat([broad_membership, fallback[['canonical_UniProtAC', 'family', 'label']]], ignore_index=True).drop_duplicates()
    return combined


def drop_other_when_informative(membership: pd.DataFrame) -> pd.DataFrame:
    informative = membership[membership['label'] != 'Other / unclassified'].groupby('canonical_UniProtAC').size()
    if informative.empty:
        return membership
    informative_acc = set(informative.index.astype(str))
    keep = ~(
        membership['canonical_UniProtAC'].astype(str).isin(informative_acc) &
        membership['label'].eq('Other / unclassified')
    )
    return membership[keep].copy()


def propagate_labels_by_gene_consensus(protein_meta: pd.DataFrame, membership: pd.DataFrame) -> pd.DataFrame:
    meta_cols = [c for c in ['canonical_UniProtAC', 'substrate_genename', 'genename'] if c in protein_meta.columns]
    gene_frame = protein_meta[meta_cols].drop_duplicates('canonical_UniProtAC').copy()
    merged = gene_frame.merge(membership, how='inner', on='canonical_UniProtAC')
    informative = merged[merged['label'] != 'Other / unclassified'].copy()
    if informative.empty:
        return membership

    gene_to_labels: dict[str, set[str]] = {}
    for row in informative.itertuples(index=False):
        labels = [x.strip() for x in str(row.label).split(';') if x.strip() and x.strip() != 'Other / unclassified']
        if not labels:
            continue
        genes = extract_gene_tokens([
            getattr(row, 'substrate_genename', ''),
            getattr(row, 'genename', ''),
        ])
        for gene in genes:
            gene_to_labels.setdefault(gene, set()).update(labels)

    informative_acc = set(informative['canonical_UniProtAC'].astype(str))
    rescue_rows = []
    for row in gene_frame.itertuples(index=False):
        acc = str(row.canonical_UniProtAC)
        if acc in informative_acc:
            continue
        genes = extract_gene_tokens([
            getattr(row, 'substrate_genename', ''),
            getattr(row, 'genename', ''),
        ])
        rescue_labels: set[str] = set()
        for gene in genes:
            rescue_labels.update(gene_to_labels.get(gene, set()))
        for label in sorted(rescue_labels):
            rescue_rows.append({
                'canonical_UniProtAC': acc,
                'family': membership['family'].iloc[0],
                'label': label,
            })
    if not rescue_rows:
        return membership
    return pd.concat([membership, pd.DataFrame(rescue_rows)], ignore_index=True).drop_duplicates()


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
    ap.add_argument('--uniprot-metadata', default='')
    ap.add_argument('--interpro-intervals', default='')
    ap.add_argument('--ontology', default='config/functional_ontology.json')
    ap.add_argument('--domain-ontology', default='config/domain_class_ontology.json')
    ap.add_argument('--outdir', default='results/functional_class_union')
    args = ap.parse_args()

    master = pd.read_csv(args.base_master, sep='\t', low_memory=False)
    sites = pd.read_csv(args.integrated_sites, sep='\t', low_memory=False)
    ontology = load_json(args.ontology)
    domain_ontology = load_json(args.domain_ontology) if args.domain_ontology else None
    apply_paper_style()
    metadata = None
    if args.uniprot_metadata:
        metadata = pd.read_csv(args.uniprot_metadata, sep='\t', low_memory=False)
    interpro = None
    if args.interpro_intervals:
        interpro = pd.read_csv(args.interpro_intervals, sep='\t', low_memory=False)

    master['canonical_UniProtAC'] = master['substrate_UniProtAC'].astype(str).map(canonicalize_uniprot_accession)
    protein_meta = build_combined_protein_metadata(master, sites, metadata=metadata, interpro=interpro, domain_ontology=domain_ontology)
    target_column = 'canonical_UniProtAC' if 'canonical_UniProtAC' in sites.columns else 'substrate_UniProtAC'
    target_acc = set(sites[target_column].dropna().astype(str).map(canonicalize_uniprot_accession))

    broad_membership = expand_membership(protein_meta, compile_label_defs(ontology['broad_labels']), 'broad')
    subclass_membership = expand_membership(protein_meta, compile_label_defs(ontology['subclass_labels']), 'subclass')
    subclass_membership = propagate_labels_by_gene_consensus(protein_meta, subclass_membership)
    subclass_membership = drop_other_when_informative(subclass_membership)
    broad_membership = add_broad_fallback_from_subclass(broad_membership, subclass_membership)
    broad_membership = propagate_labels_by_gene_consensus(protein_meta, broad_membership)
    broad_membership = drop_other_when_informative(broad_membership)
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
    save_table(protein_summary[protein_summary['functional_class_broad'].fillna('').str.contains('Other / unclassified')], outdir / 'other_unclassified_broad_proteins.tsv')
    save_table(protein_summary[protein_summary['functional_class_subclass'].fillna('').str.contains('Other / unclassified')], outdir / 'other_unclassified_subclass_proteins.tsv')
    save_table(broad, outdir / 'arg_methyl_functional_class_broad_enrichment.tsv')
    save_table(sub, outdir / 'arg_methyl_functional_class_subclass_enrichment.tsv')

    pb = broad[broad['target_count'] >= 10].sort_values('log2_odds_ratio').copy()
    broad_colors = signed_bar_colors(pb['log2_odds_ratio'])
    broad_colors = [BREWER_COLORS['mid_gray'] if cat == 'Other / unclassified' else color for cat, color in zip(pb['category'], broad_colors)]
    fig, ax = plt.subplots(figsize=(9.2, 6.4))
    ax.barh(pb['category'], pb['log2_odds_ratio'], color=broad_colors, edgecolor='white', linewidth=0.8)
    style_axis(ax, zero='x')
    ax.set_xlabel('log2(OR) vs all proteins')
    ax.set_ylabel('Broad class')
    ax.set_title('Integrated arg-methylome: multi-label broad enrichment')
    set_symmetric_xlim(ax, pb['log2_odds_ratio'], annotation_pad_ratio=0.65, center_on_zero=False)
    annotate_barh(
        ax,
        pb['category'],
        pb['log2_odds_ratio'],
        [f'n={n}, p={format_p_value(p)}, q={format_p_value(q)}' for n, p, q in zip(pb['target_count'], pb['p_value'], pb['q_value'])],
        fontsize=8.2,
    )
    fig.text(
        0.99,
        0.01,
        'Protein-level multi-label counts; totals are not additive.',
        ha='right',
        va='bottom',
        fontsize=8,
        color=BREWER_COLORS['dark_gray'],
    )
    fig.tight_layout()
    save_figure(fig, outdir / 'arg_methyl_functional_class_broad_enrichment')

    ps = balanced_category_subset(sub, min_count=8, top_positive=10, top_negative=6)
    subclass_colors = signed_bar_colors(ps['log2_odds_ratio'], positive=BREWER_COLORS['green'], negative=BREWER_COLORS['purple'])
    subclass_colors = [BREWER_COLORS['mid_gray'] if cat == 'Other / unclassified' else color for cat, color in zip(ps['category'], subclass_colors)]
    fig2, ax2 = plt.subplots(figsize=(10.2, 7.4))
    ax2.barh(ps['category'], ps['log2_odds_ratio'], color=subclass_colors, edgecolor='white', linewidth=0.8)
    style_axis(ax2, zero='x')
    ax2.set_xlabel('log2(OR) vs all proteins')
    ax2.set_ylabel('Subclass')
    ax2.set_title('Integrated arg-methylome: multi-label RNA-centric subclasses')
    set_symmetric_xlim(ax2, ps['log2_odds_ratio'], annotation_pad_ratio=0.58, center_on_zero=False)
    annotate_barh(
        ax2,
        ps['category'],
        ps['log2_odds_ratio'],
        [f'n={n}, p={format_p_value(p)}, q={format_p_value(q)}' for n, p, q in zip(ps['target_count'], ps['p_value'], ps['q_value'])],
        fontsize=7.8,
    )
    fig2.text(
        0.99,
        0.01,
        'Protein-level multi-label counts; totals are not additive.',
        ha='right',
        va='bottom',
        fontsize=8,
        color=BREWER_COLORS['dark_gray'],
    )
    fig2.tight_layout()
    save_figure(fig2, outdir / 'arg_methyl_functional_class_subclass_enrichment')


if __name__ == '__main__':
    main()

from __future__ import annotations

import argparse
from collections import defaultdict, deque
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist

from common import BREWER_COLORS, apply_paper_style, save_figure, save_table


def choose_gene_name(row: pd.Series) -> str:
    for col in ['substrate_genename', 'gene_name', 'genename']:
        value = str(row.get(col, '') or '').strip()
        if value and value not in {'nan', 'None'}:
            return value.split(';')[0].strip()
    return str(row.get('canonical_UniProtAC', ''))


def cluster_order(matrix: pd.DataFrame, axis: int) -> list[int]:
    if axis == 0:
        data = matrix.to_numpy(dtype=float)
    else:
        data = matrix.to_numpy(dtype=float).T
    if data.shape[0] <= 1:
        return list(range(data.shape[0]))
    condensed = pdist(data, metric='hamming')
    if condensed.size == 0:
        return list(range(data.shape[0]))
    return leaves_list(linkage(condensed, method='average')).tolist()


def connected_components(edges: pd.DataFrame) -> pd.DataFrame:
    graph: dict[str, set[str]] = defaultdict(set)
    for row in edges.itertuples(index=False):
        graph[str(row.protein_a)].add(str(row.protein_b))
        graph[str(row.protein_b)].add(str(row.protein_a))

    seen: set[str] = set()
    rows: list[dict] = []
    component_id = 0
    for node in sorted(graph):
        if node in seen:
            continue
        component_id += 1
        queue = deque([node])
        members: list[str] = []
        seen.add(node)
        while queue:
            current = queue.popleft()
            members.append(current)
            for neighbor in graph[current]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
        for member in sorted(members):
            rows.append({
                'community_id': component_id,
                'protein_label': member,
                'community_size': len(members),
            })
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--integrated-sites', required=True)
    ap.add_argument('--condensate-enrichment', required=True)
    ap.add_argument('--condensate-membership', required=True)
    ap.add_argument('--functional-labels', required=True)
    ap.add_argument('--outdir', default='results/assembly_clustering')
    ap.add_argument('--q-cutoff', type=float, default=0.1)
    ap.add_argument('--min-member-count', type=int, default=6)
    ap.add_argument('--top-assemblies', type=int, default=14)
    ap.add_argument('--top-proteins', type=int, default=40)
    args = ap.parse_args()

    apply_paper_style()
    sites = pd.read_csv(args.integrated_sites, sep='\t', low_memory=False)
    enrich = pd.read_csv(args.condensate_enrichment, sep='\t', low_memory=False)
    membership = pd.read_csv(args.condensate_membership, sep='\t', low_memory=False)
    labels = pd.read_csv(args.functional_labels, sep='\t', low_memory=False)

    site_counts = (
        sites.groupby('canonical_UniProtAC', dropna=False)
        .agg(methyl_site_count=('site_key', 'nunique'), substrate_genename=('substrate_genename', 'first'))
        .reset_index()
    )
    labels = labels[['canonical_UniProtAC', 'functional_class_broad', 'functional_class_subclass']].drop_duplicates()
    protein_meta = site_counts.merge(labels, how='left', on='canonical_UniProtAC')
    protein_meta['gene_name'] = protein_meta.apply(choose_gene_name, axis=1)
    protein_meta['protein_label'] = protein_meta['gene_name'] + ' (' + protein_meta['canonical_UniProtAC'].astype(str) + ')'

    selected_assemblies = enrich[
        (enrich['member_count'] >= args.min_member_count) &
        (enrich['q_value'] <= args.q_cutoff)
    ].copy()
    if selected_assemblies.empty:
        selected_assemblies = enrich[enrich['member_count'] >= args.min_member_count].copy()
    selected_assemblies = selected_assemblies.sort_values(['log2_odds_ratio', 'methyl_protein_count'], ascending=[False, False]).head(args.top_assemblies).copy()

    membership = membership[membership['condensate_id'].isin(set(selected_assemblies['condensate_id'].astype(str)))].copy()
    methyl_proteins = set(site_counts['canonical_UniProtAC'].astype(str))
    membership = membership[membership['canonical_UniProtAC'].astype(str).isin(methyl_proteins)].copy()

    protein_membership_counts = (
        membership.groupby('canonical_UniProtAC', dropna=False)
        .agg(selected_assembly_count=('condensate_id', 'nunique'))
        .reset_index()
    )
    protein_meta = protein_meta.merge(protein_membership_counts, how='left', on='canonical_UniProtAC')
    protein_meta['selected_assembly_count'] = protein_meta['selected_assembly_count'].fillna(0).astype(int)
    selected_proteins = (
        protein_meta[protein_meta['selected_assembly_count'] > 0]
        .sort_values(['selected_assembly_count', 'methyl_site_count', 'gene_name'], ascending=[False, False, True])
        .head(args.top_proteins)
        .copy()
    )

    membership = membership[membership['canonical_UniProtAC'].isin(set(selected_proteins['canonical_UniProtAC'].astype(str)))].copy()
    membership = membership.merge(
        selected_assemblies[['condensate_id', 'condensate_name', 'log2_odds_ratio', 'q_value', 'member_count', 'methyl_protein_count']],
        how='left',
        on=['condensate_id', 'condensate_name'],
    )
    membership = membership.merge(
        selected_proteins[['canonical_UniProtAC', 'protein_label', 'methyl_site_count', 'functional_class_broad']],
        how='left',
        on='canonical_UniProtAC',
    )

    matrix = (
        membership.assign(value=1)
        .pivot_table(index='protein_label', columns='condensate_name', values='value', aggfunc='max', fill_value=0)
        .astype(int)
    )
    row_order = cluster_order(matrix, axis=0)
    col_order = cluster_order(matrix, axis=1)
    matrix = matrix.iloc[row_order, col_order]

    selected_proteins = selected_proteins.set_index('protein_label').loc[matrix.index].reset_index()
    assembly_meta = selected_assemblies.set_index('condensate_name').loc[matrix.columns].reset_index()

    edges = []
    assembly_sets = {
        protein: set(matrix.columns[matrix.loc[protein] > 0].tolist())
        for protein in matrix.index
    }
    protein_labels = list(matrix.index)
    for i, protein_a in enumerate(protein_labels):
        for protein_b in protein_labels[i + 1:]:
            shared = assembly_sets[protein_a] & assembly_sets[protein_b]
            if not shared:
                continue
            union = assembly_sets[protein_a] | assembly_sets[protein_b]
            edges.append({
                'protein_a': protein_a,
                'protein_b': protein_b,
                'shared_assembly_count': len(shared),
                'shared_assemblies': ';'.join(sorted(shared)),
                'jaccard_similarity': len(shared) / len(union) if union else np.nan,
            })
    edge_df = pd.DataFrame(edges).sort_values(['shared_assembly_count', 'jaccard_similarity'], ascending=[False, False]) if edges else pd.DataFrame(columns=['protein_a', 'protein_b', 'shared_assembly_count', 'shared_assemblies', 'jaccard_similarity'])
    community_df = connected_components(edge_df) if not edge_df.empty else pd.DataFrame(columns=['community_id', 'protein_label', 'community_size'])

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(selected_assemblies, outdir / 'selected_assemblies.tsv')
    save_table(selected_proteins, outdir / 'selected_proteins.tsv')
    save_table(matrix.reset_index(), outdir / 'protein_assembly_membership_matrix.tsv')
    save_table(edge_df, outdir / 'protein_assembly_shared_membership_edges.tsv')
    save_table(community_df, outdir / 'protein_assembly_communities.tsv')

    fig, ax = plt.subplots(figsize=(12.5, max(6.0, 0.28 * len(matrix.index) + 2.2)))
    ax.imshow(matrix.to_numpy(dtype=float), aspect='auto', cmap=plt.cm.Greys, interpolation='nearest', vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(matrix.columns)))
    ax.set_xticklabels(matrix.columns, rotation=45, ha='right', fontsize=9)
    ax.set_yticks(np.arange(len(matrix.index)))
    ax.set_yticklabels(matrix.index, fontsize=8.4)
    ax.set_xlabel('Enriched assemblies / condensates')
    ax.set_ylabel('Top methylarginine proteins')
    ax.set_title('Clustering of methylarginine proteins by enriched assembly membership')
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(matrix.columns), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(matrix.index), 1), minor=True)
    ax.grid(which='minor', color='white', linestyle='-', linewidth=0.6)
    ax.tick_params(which='minor', bottom=False, left=False)
    fig.text(
        0.99,
        0.01,
        'Rows: top methyl proteins by shared enriched-assembly membership. Columns: enriched CD-CODE assemblies.',
        ha='right',
        va='bottom',
        fontsize=8,
        color=BREWER_COLORS['dark_gray'],
    )
    fig.tight_layout()
    save_figure(fig, outdir / 'protein_assembly_membership_clustering')

    print({
        'selected_assemblies': len(selected_assemblies),
        'selected_proteins': len(selected_proteins),
        'shared_edges': len(edge_df),
        'outdir': str(outdir),
    })


if __name__ == '__main__':
    main()

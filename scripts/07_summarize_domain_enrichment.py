from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from common import (
    classify_by_patterns,
    fisher_like_enrichment,
    format_p_value,
    load_json,
    parse_fasta,
    residue_background_from_sequences,
    save_figure,
    save_table,
)


TYPE_PRIORITY = {
    'domain': 0,
    'repeat': 1,
    'homologous_superfamily': 2,
    'family': 3,
}


def interval_priority(interpro_type: str) -> int:
    return TYPE_PRIORITY.get(str(interpro_type), 9)


def edge_distance(position: int, start: int, end: int) -> tuple[int, str]:
    if start <= position <= end:
        return 0, 'inside'
    if position < start:
        return start - position, 'N_terminal_to_domain'
    return position - end, 'C_terminal_to_domain'


def choose_interval(intervals: pd.DataFrame, position: int) -> tuple[pd.Series | None, int]:
    if intervals.empty:
        return None, -1
    overlaps = intervals[(intervals['fragment_start'] <= position) & (intervals['fragment_end'] >= position)].copy()
    if not overlaps.empty:
        overlaps['span'] = overlaps['fragment_end'] - overlaps['fragment_start']
        overlaps['priority'] = overlaps['interpro_type'].map(interval_priority)
        overlaps = overlaps.sort_values(['priority', 'span', 'fragment_start', 'interpro_accession'])
        return overlaps.iloc[0], 0
    candidates = intervals.copy()
    candidates['distance'] = candidates.apply(lambda r: edge_distance(position, int(r['fragment_start']), int(r['fragment_end']))[0], axis=1)
    candidates['span'] = candidates['fragment_end'] - candidates['fragment_start']
    candidates['priority'] = candidates['interpro_type'].map(interval_priority)
    candidates = candidates.sort_values(['distance', 'priority', 'span', 'fragment_start', 'interpro_accession'])
    best = candidates.iloc[0]
    return best, int(best['distance'])


def is_inter_domain_linker(intervals: pd.DataFrame, position: int) -> bool:
    if len(intervals) < 2:
        return False
    ordered = intervals[['fragment_start', 'fragment_end']].sort_values(['fragment_start', 'fragment_end']).drop_duplicates()
    previous_end = None
    for row in ordered.itertuples(index=False):
        start = int(row.fragment_start)
        end = int(row.fragment_end)
        if previous_end is not None and previous_end < position < start:
            return True
        previous_end = max(previous_end or end, end)
    return False


def classify_residue(intervals_by_acc: dict[str, pd.DataFrame], accession: str, position: int, boundary_window: int, ontology: dict) -> dict:
    intervals = intervals_by_acc.get(accession, pd.DataFrame())
    if intervals.empty:
        return {
            'domain_context_class': 'no_domain_annotation',
            'nearest_domain_class': 'No domain annotation',
            'nearest_domain_edge_distance': pd.NA,
        }
    best, distance = choose_interval(intervals, position)
    if best is None:
        return {
            'domain_context_class': 'no_domain_annotation',
            'nearest_domain_class': 'No domain annotation',
            'nearest_domain_edge_distance': pd.NA,
        }
    if distance == 0:
        context_class = 'in_domain'
    elif distance <= boundary_window:
        context_class = 'boundary'
    elif is_inter_domain_linker(intervals, position):
        context_class = 'inter_domain_linker'
    else:
        context_class = 'distal'
    label_text = ' '.join([
        str(best.get('interpro_name', '')),
        str(best.get('interpro_type', '')),
        str(best.get('member_database_accessions', '')),
    ])
    return {
        'domain_context_class': context_class,
        'nearest_domain_class': classify_by_patterns(label_text, ontology['domain_classes'], default='Other domain'),
        'nearest_domain_edge_distance': distance,
    }


def annotate_background(background: pd.DataFrame, intervals_by_acc: dict[str, pd.DataFrame], boundary_window: int, ontology: dict) -> pd.DataFrame:
    annotations = background.apply(
        lambda row: pd.Series(
            classify_residue(
                intervals_by_acc=intervals_by_acc,
                accession=str(row['canonical_UniProtAC']),
                position=int(row['position']),
                boundary_window=boundary_window,
                ontology=ontology,
            )
        ),
        axis=1,
    )
    return pd.concat([background, annotations], axis=1)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--annotated-sites', required=True)
    ap.add_argument('--interpro-intervals', required=True)
    ap.add_argument('--canonical-fasta', required=True)
    ap.add_argument('--ontology', default='config/domain_class_ontology.json')
    ap.add_argument('--boundary-window', type=int, default=20)
    ap.add_argument('--outdir', default='results/domain_context')
    args = ap.parse_args()

    annotated_sites = pd.read_csv(args.annotated_sites, sep='\t', low_memory=False)
    intervals = pd.read_csv(args.interpro_intervals, sep='\t', low_memory=False)
    intervals['fragment_start'] = pd.to_numeric(intervals['fragment_start'], errors='coerce')
    intervals['fragment_end'] = pd.to_numeric(intervals['fragment_end'], errors='coerce')
    intervals = intervals[intervals['fragment_start'].notna() & intervals['fragment_end'].notna()].copy()
    intervals['fragment_start'] = intervals['fragment_start'].astype(int)
    intervals['fragment_end'] = intervals['fragment_end'].astype(int)
    intervals_by_acc = {
        accession: group.sort_values(['fragment_start', 'fragment_end', 'interpro_accession']).reset_index(drop=True)
        for accession, group in intervals.groupby('canonical_UniProtAC')
    }
    ontology = load_json(args.ontology)

    if 'nearest_domain_class' not in annotated_sites.columns:
        class_text = annotated_sites.fillna('').astype(str)[['nearest_interpro_name', 'nearest_interpro_type', 'nearest_member_database_accessions']].agg(' '.join, axis=1)
        annotated_sites['nearest_domain_class'] = class_text.map(lambda x: classify_by_patterns(x, ontology['domain_classes'], default='Other domain'))

    sequences = parse_fasta(args.canonical_fasta)
    proteome_arg = residue_background_from_sequences(sequences, residue='R')
    methylated_proteins = sorted(set(annotated_sites['canonical_UniProtAC'].dropna().astype(str)))
    methyl_protein_arg = proteome_arg[proteome_arg['canonical_UniProtAC'].isin(methylated_proteins)].copy()
    methyl_protein_arg = annotate_background(methyl_protein_arg, intervals_by_acc, args.boundary_window, ontology)

    context_enrichment = fisher_like_enrichment(
        annotated_sites['domain_context_class'],
        methyl_protein_arg['domain_context_class'],
    )
    domain_class_enrichment = fisher_like_enrichment(
        annotated_sites['nearest_domain_class'],
        methyl_protein_arg['nearest_domain_class'],
    )
    edge_summary = annotated_sites[['nearest_domain_edge_distance', 'domain_context_class']].copy()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(context_enrichment, outdir / 'domain_context_enrichment.tsv')
    save_table(domain_class_enrichment, outdir / 'domain_class_enrichment.tsv')
    save_table(edge_summary, outdir / 'domain_edge_distance_summary.tsv')

    plot_context = context_enrichment[context_enrichment['category'].isin(['in_domain', 'boundary', 'inter_domain_linker', 'distal'])].sort_values('odds_ratio')
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    ax.barh(plot_context['category'], plot_context['log2_odds_ratio'])
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlabel('log2(OR) vs arginines in methylated proteins')
    ax.set_ylabel('Domain context')
    ax.set_title('Methylarginine enrichment by domain context')
    for y, v, n, q in zip(plot_context['category'], plot_context['log2_odds_ratio'], plot_context['target_count'], plot_context['q_value_bh']):
        ax.text(v, y, f'  n={n}, q={format_p_value(q)}', va='center', ha='left' if v >= 0 else 'right', fontsize=9)
    fig.tight_layout()
    save_figure(fig, outdir / 'arg_methyl_domain_context_enrichment')

    plot_domain = domain_class_enrichment[domain_class_enrichment['target_count'] >= 20].head(12).sort_values('odds_ratio')
    fig2, ax2 = plt.subplots(figsize=(9.0, 6.0))
    ax2.barh(plot_domain['category'], plot_domain['log2_odds_ratio'])
    ax2.axvline(0, color='black', linewidth=0.8)
    ax2.set_xlabel('log2(OR) vs arginines in methylated proteins')
    ax2.set_ylabel('Nearest domain class')
    ax2.set_title('Nearest domain-class enrichment around methylarginines')
    for y, v, n, q in zip(plot_domain['category'], plot_domain['log2_odds_ratio'], plot_domain['target_count'], plot_domain['q_value_bh']):
        ax2.text(v, y, f'  n={n}, q={format_p_value(q)}', va='center', ha='left' if v >= 0 else 'right', fontsize=9)
    fig2.tight_layout()
    save_figure(fig2, outdir / 'arg_methyl_nearest_domain_class_enrichment')

    print({
        'annotated_sites': len(annotated_sites),
        'background_arginines_methylated_proteins': len(methyl_protein_arg),
        'outputs': str(outdir),
    })


if __name__ == '__main__':
    main()

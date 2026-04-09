from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import (
    canonical_site_table,
    fisher_like_enrichment,
    format_p_value,
    parse_fasta,
    residue_background_from_sequences,
    save_figure,
    save_table,
)


def normalize_disorder_intervals(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    rename_map = {}
    for source, target in [
        ('accession', 'canonical_UniProtAC'),
        ('Entry', 'canonical_UniProtAC'),
        ('start', 'fragment_start'),
        ('end', 'fragment_end'),
    ]:
        if source in out.columns and target not in out.columns:
            rename_map[source] = target
    out = out.rename(columns=rename_map)
    required = ['canonical_UniProtAC', 'fragment_start', 'fragment_end']
    missing = [col for col in required if col not in out.columns]
    if missing:
        raise ValueError(f'missing disorder interval columns: {missing}')
    out['canonical_UniProtAC'] = out['canonical_UniProtAC'].astype(str).str.strip()
    out['fragment_start'] = pd.to_numeric(out['fragment_start'], errors='coerce')
    out['fragment_end'] = pd.to_numeric(out['fragment_end'], errors='coerce')
    out = out[out['fragment_start'].notna() & out['fragment_end'].notna()].copy()
    out['fragment_start'] = out['fragment_start'].astype(int)
    out['fragment_end'] = out['fragment_end'].astype(int)
    return out.sort_values(['canonical_UniProtAC', 'fragment_start', 'fragment_end']).reset_index(drop=True)


def classify_positions(intervals: pd.DataFrame, positions: np.ndarray, boundary_window: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = len(positions)
    classes = np.full(n, 'no_disorder_annotation', dtype=object)
    distances = np.full(n, np.nan, dtype=float)
    has_annotation = np.zeros(n, dtype=bool)
    if intervals.empty or n == 0:
        return classes, distances, has_annotation

    starts = intervals['fragment_start'].to_numpy(dtype=int)
    ends = intervals['fragment_end'].to_numpy(dtype=int)
    min_distance = np.full(n, np.iinfo(np.int32).max, dtype=int)
    in_interval = np.zeros(n, dtype=bool)

    for start, end in zip(starts, ends):
        overlap = (positions >= start) & (positions <= end)
        in_interval |= overlap
        dist = np.minimum(np.abs(positions - start), np.abs(positions - end))
        dist[overlap] = 0
        min_distance = np.minimum(min_distance, dist)

    classes[in_interval] = 'disordered'
    distances[in_interval] = 0
    boundary = (~in_interval) & (min_distance <= boundary_window)
    classes[boundary] = 'disorder_boundary'
    distances[boundary] = min_distance[boundary]
    ordered = (~in_interval) & (min_distance > boundary_window)
    classes[ordered] = 'ordered'
    distances[ordered] = min_distance[ordered]
    has_annotation[:] = True
    return classes, distances, has_annotation


def annotate_frame(df: pd.DataFrame, intervals_by_acc: dict[str, pd.DataFrame], accession_col: str, position_col: str, boundary_window: int) -> pd.DataFrame:
    out = df.copy()
    out[accession_col] = out[accession_col].astype(str).str.strip()
    out[position_col] = pd.to_numeric(out[position_col], errors='coerce')
    out['disorder_context_class'] = 'no_position'
    out['nearest_disorder_edge_distance'] = pd.NA
    out['has_disorder_annotation'] = False

    valid_mask = (~out[accession_col].isin(['', 'nan', 'None'])) & out[position_col].notna()
    if valid_mask.any():
        valid = out.loc[valid_mask, [accession_col, position_col]].copy()
        valid[position_col] = valid[position_col].astype(int)
        for accession, idx in valid.groupby(accession_col).groups.items():
            positions = valid.loc[idx, position_col].to_numpy(dtype=int)
            intervals = intervals_by_acc.get(accession, pd.DataFrame())
            classes, distances, has_annotation = classify_positions(intervals, positions, boundary_window)
            out.loc[idx, 'disorder_context_class'] = classes
            out.loc[idx, 'nearest_disorder_edge_distance'] = distances
            out.loc[idx, 'has_disorder_annotation'] = has_annotation & (classes != 'no_disorder_annotation')

    out['disorder_context_boundary_window'] = boundary_window
    out['disorder_context_layer'] = 'disorder_context'
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--sites', required=True)
    ap.add_argument('--disorder-intervals', required=True)
    ap.add_argument('--canonical-fasta', default='')
    ap.add_argument('--outdir', default='results/disorder_context')
    ap.add_argument('--accession-col', default='canonical_UniProtAC')
    ap.add_argument('--position-col', default='corrected_position')
    ap.add_argument('--boundary-window', type=int, default=20)
    args = ap.parse_args()

    sites = pd.read_csv(args.sites, sep='\t', low_memory=False)
    sites = canonical_site_table(sites, accession_col=args.accession_col, position_col=args.position_col)
    intervals = normalize_disorder_intervals(pd.read_csv(args.disorder_intervals, sep='\t', low_memory=False))
    intervals_by_acc = {
        accession: group.reset_index(drop=True)
        for accession, group in intervals.groupby('canonical_UniProtAC')
    }

    annotated_sites = annotate_frame(
        sites,
        intervals_by_acc=intervals_by_acc,
        accession_col='canonical_UniProtAC',
        position_col='corrected_position',
        boundary_window=args.boundary_window,
    )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(annotated_sites, outdir / 'sites_with_disorder_context.tsv')

    summary = annotated_sites['disorder_context_class'].value_counts(dropna=False).rename_axis('disorder_context_class').reset_index(name='site_count')
    summary['boundary_window'] = args.boundary_window
    save_table(summary, outdir / 'disorder_context_summary.tsv')

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    plot_df = summary[summary['disorder_context_class'] != 'no_position'].copy()
    ax.bar(plot_df['disorder_context_class'], plot_df['site_count'])
    ax.set_ylabel('Methylarginine sites')
    ax.set_title('Methylarginine disorder-context distribution')
    ax.tick_params(axis='x', rotation=20)
    fig.tight_layout()
    save_figure(fig, outdir / 'arg_methyl_disorder_context_distribution')

    if args.canonical_fasta:
        sequences = parse_fasta(args.canonical_fasta)
        proteome_arg = residue_background_from_sequences(sequences, residue='R')
        proteome_arg = proteome_arg.rename(columns={'position': 'corrected_position'})
        proteome_arg = annotate_frame(
            proteome_arg,
            intervals_by_acc=intervals_by_acc,
            accession_col='canonical_UniProtAC',
            position_col='corrected_position',
            boundary_window=args.boundary_window,
        )
        methylated_proteins = sorted(set(annotated_sites['canonical_UniProtAC'].dropna().astype(str)))
        methyl_protein_arg = proteome_arg[proteome_arg['canonical_UniProtAC'].isin(methylated_proteins)].copy()
        save_table(
            fisher_like_enrichment(annotated_sites['disorder_context_class'], proteome_arg['disorder_context_class']),
            outdir / 'idr_enrichment_vs_all_arginines.tsv',
        )
        save_table(
            fisher_like_enrichment(annotated_sites['disorder_context_class'], methyl_protein_arg['disorder_context_class']),
            outdir / 'idr_enrichment_vs_arginines_in_methylated_proteins.tsv',
        )



        compare = pd.read_csv(outdir / 'idr_enrichment_vs_arginines_in_methylated_proteins.tsv', sep='\t')
        compare = compare[compare['category'].isin(['disordered', 'disorder_boundary', 'ordered'])].sort_values('odds_ratio')
        fig2, ax2 = plt.subplots(figsize=(7.4, 4.8))
        ax2.barh(compare['category'], compare['log2_odds_ratio'])
        ax2.axvline(0, color='black', linewidth=0.8)
        ax2.set_xlabel('log2(OR) vs arginines in methylated proteins')
        ax2.set_ylabel('Disorder class')
        ax2.set_title('Methylarginine enrichment by disorder context')
        for y, v, n, q in zip(compare['category'], compare['log2_odds_ratio'], compare['target_count'], compare['q_value_bh']):
            ax2.text(v, y, f'  n={n}, q={format_p_value(q)}', va='center', ha='left' if v >= 0 else 'right', fontsize=9)
        fig2.tight_layout()
        save_figure(fig2, outdir / 'arg_methyl_idr_enrichment')

    print({
        'annotated_sites': len(annotated_sites),
        'sites_with_disorder_annotation': int(annotated_sites['has_disorder_annotation'].fillna(False).sum()),
        'output': str(outdir / 'sites_with_disorder_context.tsv'),
    })


if __name__ == '__main__':
    main()

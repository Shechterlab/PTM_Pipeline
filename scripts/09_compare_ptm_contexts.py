from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import (
    BREWER_COLORS,
    annotate_barh,
    apply_paper_style,
    canonicalize_uniprot_accession,
    fisher_like_enrichment,
    parse_fasta,
    q_threshold_note,
    residue_background_from_sequences,
    save_figure,
    save_table,
    set_symmetric_xlim,
    signed_bar_colors,
    style_axis,
)


DEFAULT_PTM_GROUPS = [
    'Arg methylation',
    'Phosphorylation',
    'Lys acylation',
    'Lys methylation',
    'Ubiquitin/SUMO',
]

PTM_GROUP_TARGET_RESIDUES = {
    'Arg methylation': {'R'},
    'Phosphorylation': {'S', 'T', 'Y'},
    'Lys acylation': {'K'},
    'Lys methylation': {'K'},
    'Ubiquitin/SUMO': {'K'},
}


def normalize_master_sites(master: pd.DataFrame, sequences: dict[str, str], ptm_groups: list[str]) -> pd.DataFrame:
    out = master[master['ptm_group'].isin(ptm_groups)].copy()
    out['canonical_UniProtAC'] = out['substrate_UniProtAC'].astype(str).map(canonicalize_uniprot_accession)
    out['position'] = pd.to_numeric(out['position'], errors='coerce')
    out = out[out['position'].notna()].copy()
    out['position'] = out['position'].astype(int)
    out['residue'] = out['residue'].astype(str).str.upper().str[:1]
    out = out[
        out.apply(
            lambda row: row['residue'] in PTM_GROUP_TARGET_RESIDUES.get(str(row['ptm_group']), {str(row['residue'])}),
            axis=1,
        )
    ].copy()

    def valid_row(row: pd.Series) -> bool:
        seq = sequences.get(str(row['canonical_UniProtAC']), '')
        pos = int(row['position'])
        return bool(seq) and 1 <= pos <= len(seq) and seq[pos - 1] == row['residue']

    out = out[out.apply(valid_row, axis=1)].copy()
    out['site'] = out['residue'] + out['position'].astype(str)
    out['site_key'] = out['canonical_UniProtAC'] + ':' + out['site']
    out = out.drop_duplicates(['ptm_group', 'canonical_UniProtAC', 'residue', 'position']).copy()
    return out


def normalize_arg_methyl_sites(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out['canonical_UniProtAC'] = out['canonical_UniProtAC'].astype(str).map(canonicalize_uniprot_accession)
    out['position'] = pd.to_numeric(out['corrected_position'], errors='coerce')
    out = out[out['position'].notna()].copy()
    out['position'] = out['position'].astype(int)
    out['residue'] = 'R'
    out['ptm_group'] = 'Arg methylation'
    out['site'] = out['residue'] + out['position'].astype(str)
    out['site_key'] = out['canonical_UniProtAC'] + ':' + out['site']
    return out.drop_duplicates(['ptm_group', 'canonical_UniProtAC', 'residue', 'position']).copy()


def normalize_intervals(df: pd.DataFrame, start_col: str = 'fragment_start', end_col: str = 'fragment_end') -> pd.DataFrame:
    out = df.copy()
    out['canonical_UniProtAC'] = out['canonical_UniProtAC'].astype(str).str.strip()
    out[start_col] = pd.to_numeric(out[start_col], errors='coerce')
    out[end_col] = pd.to_numeric(out[end_col], errors='coerce')
    out = out[out[start_col].notna() & out[end_col].notna()].copy()
    out[start_col] = out[start_col].astype(int)
    out[end_col] = out[end_col].astype(int)
    return out


def classify_positions(intervals_df: pd.DataFrame, positions: np.ndarray, boundary_window: int) -> np.ndarray:
    n = len(positions)
    classes = np.full(n, 'no_disorder_annotation', dtype=object)
    if intervals_df.empty or n == 0:
        return classes
    starts = intervals_df['fragment_start'].to_numpy(dtype=int)
    ends = intervals_df['fragment_end'].to_numpy(dtype=int)
    min_distance = np.full(n, np.iinfo(np.int32).max, dtype=int)
    in_interval = np.zeros(n, dtype=bool)
    for start, end in zip(starts, ends):
        overlap = (positions >= start) & (positions <= end)
        in_interval |= overlap
        dist = np.minimum(np.abs(positions - start), np.abs(positions - end))
        dist[overlap] = 0
        min_distance = np.minimum(min_distance, dist)
    classes[in_interval] = 'disordered'
    boundary = (~in_interval) & (min_distance <= boundary_window)
    classes[boundary] = 'disorder_boundary'
    ordered = (~in_interval) & (min_distance > boundary_window)
    classes[ordered] = 'ordered'
    return classes


def annotate_disorder(df: pd.DataFrame, intervals: pd.DataFrame, boundary_window: int) -> pd.DataFrame:
    intervals_by_acc = {acc: group[['fragment_start', 'fragment_end']].reset_index(drop=True) for acc, group in intervals.groupby('canonical_UniProtAC')}
    out = df.copy()
    out['disorder_context_class'] = 'no_disorder_annotation'
    for accession, idx in out.groupby('canonical_UniProtAC').groups.items():
        positions = out.loc[idx, 'position'].to_numpy(dtype=int)
        out.loc[idx, 'disorder_context_class'] = classify_positions(
            intervals_by_acc.get(accession, pd.DataFrame(columns=['fragment_start', 'fragment_end'])),
            positions,
            boundary_window=boundary_window,
        )
    return out


def add_binary_idr_label(df: pd.DataFrame, source_col: str = 'disorder_context_class') -> pd.DataFrame:
    out = df.copy()
    out['idr_binary_class'] = pd.NA
    out.loc[out[source_col] == 'disordered', 'idr_binary_class'] = 'IDR'
    out.loc[out[source_col].isin(['ordered', 'disorder_boundary']), 'idr_binary_class'] = 'non_IDR'
    return out


def interval_priority(interpro_type: str) -> int:
    return {'domain': 0, 'repeat': 1, 'homologous_superfamily': 2, 'family': 3}.get(str(interpro_type), 9)


def annotate_domain(df: pd.DataFrame, intervals: pd.DataFrame, boundary_window: int) -> pd.DataFrame:
    intervals_by_acc = {
        acc: group.sort_values(['fragment_start', 'fragment_end', 'interpro_accession']).reset_index(drop=True)
        for acc, group in intervals.groupby('canonical_UniProtAC')
    }

    out = df.copy()
    out['domain_context_class'] = 'no_domain_annotation'
    for accession, idx in out.groupby('canonical_UniProtAC').groups.items():
        group = intervals_by_acc.get(str(accession), pd.DataFrame())
        if group.empty:
            continue
        positions = out.loc[idx, 'position'].to_numpy(dtype=int)
        classes = np.full(len(positions), 'distal', dtype=object)
        min_distance = np.full(len(positions), np.iinfo(np.int32).max, dtype=int)
        in_domain = np.zeros(len(positions), dtype=bool)
        starts = group['fragment_start'].to_numpy(dtype=int)
        ends = group['fragment_end'].to_numpy(dtype=int)
        for start, end in zip(starts, ends):
            overlap = (positions >= start) & (positions <= end)
            in_domain |= overlap
            dist = np.minimum(np.abs(positions - start), np.abs(positions - end))
            dist[overlap] = 0
            min_distance = np.minimum(min_distance, dist)
        classes[min_distance <= boundary_window] = 'boundary'
        classes[in_domain] = 'in_domain'
        out.loc[idx, 'domain_context_class'] = classes
    return out


def make_context_background(sequences: dict[str, str], residues: list[str]) -> pd.DataFrame:
    frames = [residue_background_from_sequences(sequences, residue=residue) for residue in sorted(set(residues))]
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--base-master', required=True)
    ap.add_argument('--integrated-arg-sites', required=True)
    ap.add_argument('--canonical-fasta', required=True)
    ap.add_argument('--disorder-intervals', default='')
    ap.add_argument('--interpro-intervals', default='')
    ap.add_argument('--boundary-window', type=int, default=20)
    ap.add_argument('--ptm-groups', nargs='*', default=list(DEFAULT_PTM_GROUPS))
    ap.add_argument('--make-detailed-context-plots', action='store_true')
    ap.add_argument('--outdir', default='results/ptm_compare')
    args = ap.parse_args()

    sequences = parse_fasta(args.canonical_fasta)
    master = pd.read_csv(args.base_master, sep='\t', low_memory=False)
    arg_sites = pd.read_csv(args.integrated_arg_sites, sep='\t', low_memory=False)
    apply_paper_style()

    ptm_sites = normalize_master_sites(master, sequences=sequences, ptm_groups=[g for g in args.ptm_groups if g != 'Arg methylation'])
    arg_sites = normalize_arg_methyl_sites(arg_sites)
    ptm_sites = pd.concat([ptm_sites, arg_sites], ignore_index=True, sort=False)
    input_qc = (
        ptm_sites.groupby('ptm_group', dropna=False)
        .agg(
            site_count=('site_key', 'nunique'),
            protein_count=('canonical_UniProtAC', 'nunique'),
            residues=('residue', lambda s: ';'.join(sorted(set(s.astype(str))))),
        )
        .reset_index()
    )

    rows = []
    backgrounds = []
    disorder_bg_cache: dict[tuple[str, ...], pd.DataFrame] = {}
    domain_bg_cache: dict[tuple[str, ...], pd.DataFrame] = {}

    if args.disorder_intervals:
        disorder = normalize_intervals(pd.read_csv(args.disorder_intervals, sep='\t', low_memory=False))
        ptm_sites = annotate_disorder(ptm_sites, disorder, boundary_window=args.boundary_window)
        ptm_sites = add_binary_idr_label(ptm_sites)
        for ptm_group, group in ptm_sites.groupby('ptm_group'):
            residues = sorted(set(group['residue'].astype(str)))
            residue_key = tuple(residues)
            if residue_key not in disorder_bg_cache:
                bg = make_context_background(sequences, residues)
                bg = annotate_disorder(bg, disorder, boundary_window=args.boundary_window)
                bg = add_binary_idr_label(bg)
                disorder_bg_cache[residue_key] = bg
            bg = disorder_bg_cache[residue_key].copy()
            bg = bg[~bg['site_key'].isin(set(group['site_key'].astype(str)))].copy()
            enrich = fisher_like_enrichment(group['disorder_context_class'], bg['disorder_context_class'])
            enrich['ptm_group'] = ptm_group
            enrich['context_layer'] = 'disorder'
            rows.append(enrich)
            backgrounds.append({'ptm_group': ptm_group, 'context_layer': 'disorder', 'background_residues': len(bg)})

        binary_rows = []
        binary_counts = []
        for ptm_group, group in ptm_sites.groupby('ptm_group'):
            residues = sorted(set(group['residue'].astype(str)))
            residue_key = tuple(residues)
            if residue_key not in disorder_bg_cache:
                bg = make_context_background(sequences, residues)
                bg = annotate_disorder(bg, disorder, boundary_window=args.boundary_window)
                bg = add_binary_idr_label(bg)
                disorder_bg_cache[residue_key] = bg
            bg = disorder_bg_cache[residue_key].copy()
            bg = bg[~bg['site_key'].isin(set(group['site_key'].astype(str)))].copy()
            group_binary = group[group['idr_binary_class'].notna()].copy()
            bg_binary = bg[bg['idr_binary_class'].notna()].copy()
            enrich = fisher_like_enrichment(group_binary['idr_binary_class'], bg_binary['idr_binary_class'])
            enrich['ptm_group'] = ptm_group
            enrich['comparison'] = 'binary_idr_vs_non_idr'
            binary_rows.append(enrich)
            binary_counts.append({
                'ptm_group': ptm_group,
                'classified_site_count': len(group_binary),
                'idr_site_count': int((group_binary['idr_binary_class'] == 'IDR').sum()),
                'background_classified_site_count': len(bg_binary),
                'background_idr_site_count': int((bg_binary['idr_binary_class'] == 'IDR').sum()),
            })
        if binary_rows:
            binary = pd.concat(binary_rows, ignore_index=True)
            save_table(binary, Path(args.outdir) / 'ptm_idr_binary_comparison.tsv')
            save_table(pd.DataFrame(binary_counts), Path(args.outdir) / 'ptm_idr_binary_counts.tsv')

            plot_df = binary[binary['category'] == 'IDR'].copy()
            if not plot_df.empty:
                plot_df = plot_df.merge(pd.DataFrame(binary_counts), on='ptm_group', how='left')
                plot_df['observed_idr_fraction_pct'] = 100.0 * plot_df['idr_site_count'] / plot_df['classified_site_count']
                plot_df['background_idr_fraction_pct'] = 100.0 * plot_df['background_idr_site_count'] / plot_df['background_classified_site_count']
                plot_df = plot_df.sort_values('log2_odds_ratio')

                fig_idr, ax_idr = plt.subplots(figsize=(9.2, 5.6))
                colors = signed_bar_colors(
                    plot_df['log2_odds_ratio'],
                    positive=BREWER_COLORS['dark_gray'],
                    negative=BREWER_COLORS['light_gray'],
                )
                ax_idr.barh(plot_df['ptm_group'], plot_df['log2_odds_ratio'], color=colors, edgecolor='white', linewidth=0.8)
                style_axis(ax_idr, zero='x')
                ax_idr.set_xlabel('log2(OR) for IDR vs target-residue proteome background')
                ax_idr.set_ylabel('PTM class')
                ax_idr.set_title('Cross-PTM IDR enrichment')
                set_symmetric_xlim(ax_idr, plot_df['log2_odds_ratio'], annotation_pad_ratio=0.34, center_on_zero=False)
                annotate_barh(
                    ax_idr,
                    plot_df['ptm_group'],
                    plot_df['log2_odds_ratio'],
                    [f'{n}/{t}' for n, t in zip(plot_df['idr_site_count'], plot_df['classified_site_count'])],
                    fontsize=8.7,
                )
                fig_idr.text(0.99, 0.01, q_threshold_note(plot_df['q_value']), ha='right', va='bottom', fontsize=8, color=BREWER_COLORS['dark_gray'])
                fig_idr.tight_layout()
                save_figure(fig_idr, Path(args.outdir) / 'cross_ptm_idr_binary_comparison')

                fig_frac, ax_frac = plt.subplots(figsize=(9.2, 5.6))
                ax_frac.barh(
                    plot_df['ptm_group'],
                    plot_df['observed_idr_fraction_pct'],
                    color=BREWER_COLORS['mid_gray'],
                    edgecolor='white',
                    linewidth=0.8,
                )
                ax_frac.scatter(
                    plot_df['background_idr_fraction_pct'],
                    plot_df['ptm_group'],
                    color=BREWER_COLORS['dark_gray'],
                    marker='D',
                    s=28,
                    zorder=3,
                    label='Target-residue proteome background',
                )
                style_axis(ax_frac)
                ax_frac.set_xlabel('% of classified sites in IDRs')
                ax_frac.set_ylabel('PTM class')
                ax_frac.set_title('Cross-PTM IDR site burden')
                ax_frac.legend(loc='lower right', frameon=False)
                annotate_barh(
                    ax_frac,
                    plot_df['ptm_group'],
                    plot_df['observed_idr_fraction_pct'],
                    [f'{n}/{t}' for n, t in zip(plot_df['idr_site_count'], plot_df['classified_site_count'])],
                    fontsize=8.0,
                )
                fig_frac.tight_layout()
                save_figure(fig_frac, Path(args.outdir) / 'cross_ptm_idr_site_fraction')

    if args.interpro_intervals:
        interpro = normalize_intervals(pd.read_csv(args.interpro_intervals, sep='\t', low_memory=False))
        ptm_sites = annotate_domain(ptm_sites, interpro, boundary_window=args.boundary_window)
        for ptm_group, group in ptm_sites.groupby('ptm_group'):
            residues = sorted(set(group['residue'].astype(str)))
            residue_key = tuple(residues)
            if residue_key not in domain_bg_cache:
                bg = make_context_background(sequences, residues)
                bg = annotate_domain(bg, interpro, boundary_window=args.boundary_window)
                domain_bg_cache[residue_key] = bg
            bg = domain_bg_cache[residue_key].copy()
            bg = bg[~bg['site_key'].isin(set(group['site_key'].astype(str)))].copy()
            enrich = fisher_like_enrichment(group['domain_context_class'], bg['domain_context_class'])
            enrich['ptm_group'] = ptm_group
            enrich['context_layer'] = 'domain'
            rows.append(enrich)
            backgrounds.append({'ptm_group': ptm_group, 'context_layer': 'domain', 'background_residues': len(bg)})

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(ptm_sites, outdir / 'ptm_sites_with_contexts.tsv')
    save_table(input_qc, outdir / 'ptm_input_qc.tsv')
    if rows:
        comparison = pd.concat(rows, ignore_index=True)
        save_table(comparison, outdir / 'ptm_context_comparison.tsv')
        save_table(pd.DataFrame(backgrounds), outdir / 'ptm_context_backgrounds.tsv')

        if args.make_detailed_context_plots and (comparison['context_layer'] == 'disorder').any():
            plot_df = comparison[
                (comparison['context_layer'] == 'disorder') &
                (comparison['category'].isin(['disordered', 'disorder_boundary', 'ordered']))
            ].copy()
            if not plot_df.empty:
                pivot = plot_df.pivot(index='ptm_group', columns='category', values='log2_odds_ratio').fillna(0.0)
                fig, ax = plt.subplots(figsize=(9.5, 5.5))
                for category in pivot.columns:
                    ax.plot(pivot.index, pivot[category], marker='o', label=category)
                    category_df = plot_df[plot_df['category'] == category]
                    for _, row in category_df.iterrows():
                        ax.text(row['ptm_group'], row['log2_odds_ratio'], f" n={int(row['target_count'])}", fontsize=7.0, va='bottom' if row['log2_odds_ratio'] >= 0 else 'top')
                ax.axhline(0, color='black', linewidth=0.8)
                ax.set_ylabel('log2(OR) vs target-residue proteome background')
                ax.set_xlabel('PTM class')
                ax.set_title('Cross-PTM disorder-context comparison')
                ax.tick_params(axis='x', rotation=25)
                ax.legend()
                fig.tight_layout()
                save_figure(fig, outdir / 'cross_ptm_disorder_context_comparison')

        if args.make_detailed_context_plots and (comparison['context_layer'] == 'domain').any():
            plot_df = comparison[
                (comparison['context_layer'] == 'domain') &
                (comparison['category'].isin(['in_domain', 'boundary', 'distal']))
            ].copy()
            if not plot_df.empty:
                pivot = plot_df.pivot(index='ptm_group', columns='category', values='log2_odds_ratio').fillna(0.0)
                fig2, ax2 = plt.subplots(figsize=(9.5, 5.5))
                for category in pivot.columns:
                    ax2.plot(pivot.index, pivot[category], marker='o', label=category)
                    category_df = plot_df[plot_df['category'] == category]
                    for _, row in category_df.iterrows():
                        ax2.text(row['ptm_group'], row['log2_odds_ratio'], f" n={int(row['target_count'])}", fontsize=7.0, va='bottom' if row['log2_odds_ratio'] >= 0 else 'top')
                ax2.axhline(0, color='black', linewidth=0.8)
                ax2.set_ylabel('log2(OR) vs target-residue proteome background')
                ax2.set_xlabel('PTM class')
                ax2.set_title('Cross-PTM domain-context comparison')
                ax2.tick_params(axis='x', rotation=25)
                ax2.legend()
                fig2.tight_layout()
                save_figure(fig2, outdir / 'cross_ptm_domain_context_comparison')

    print({
        'ptm_site_rows': len(ptm_sites),
        'ptm_groups': sorted(set(ptm_sites['ptm_group'].astype(str))),
        'output': str(outdir),
    })


if __name__ == '__main__':
    main()

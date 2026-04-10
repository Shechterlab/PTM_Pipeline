from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import (
    canonicalize_uniprot_accession,
    fisher_like_enrichment,
    format_p_value,
    parse_fasta,
    residue_background_from_sequences,
    save_figure,
    save_table,
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

    def classify(row: pd.Series) -> str:
        group = intervals_by_acc.get(str(row['canonical_UniProtAC']), pd.DataFrame())
        if group.empty:
            return 'no_domain_annotation'
        position = int(row['position'])
        overlaps = group[(group['fragment_start'] <= position) & (group['fragment_end'] >= position)].copy()
        if not overlaps.empty:
            return 'in_domain'
        group = group.copy()
        group['distance'] = group.apply(
            lambda r: min(abs(int(r['fragment_start']) - position), abs(int(r['fragment_end']) - position)),
            axis=1,
        )
        group['priority'] = group['interpro_type'].map(interval_priority)
        group = group.sort_values(['distance', 'priority', 'fragment_start', 'interpro_accession'])
        nearest = int(group.iloc[0]['distance'])
        if nearest <= boundary_window:
            return 'boundary'
        return 'distal'

    out = df.copy()
    out['domain_context_class'] = out.apply(classify, axis=1)
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
    ap.add_argument('--outdir', default='results/ptm_compare')
    args = ap.parse_args()

    sequences = parse_fasta(args.canonical_fasta)
    master = pd.read_csv(args.base_master, sep='\t', low_memory=False)
    arg_sites = pd.read_csv(args.integrated_arg_sites, sep='\t', low_memory=False)

    ptm_sites = normalize_master_sites(master, sequences=sequences, ptm_groups=[g for g in args.ptm_groups if g != 'Arg methylation'])
    arg_sites = normalize_arg_methyl_sites(arg_sites)
    ptm_sites = pd.concat([ptm_sites, arg_sites], ignore_index=True, sort=False)

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
        if binary_rows:
            binary = pd.concat(binary_rows, ignore_index=True)
            save_table(binary, Path(args.outdir) / 'ptm_idr_binary_comparison.tsv')

            plot_df = binary[binary['category'] == 'IDR'].copy()
            if not plot_df.empty:
                plot_df = plot_df.sort_values('odds_ratio', ascending=False)
                fig_idr, ax_idr = plt.subplots(figsize=(8.8, 5.2))
                ax_idr.bar(plot_df['ptm_group'], plot_df['log2_odds_ratio'])
                ax_idr.axhline(0, color='black', linewidth=0.8)
                ax_idr.set_ylabel('log2(OR) for IDR vs target-residue proteome background')
                ax_idr.set_xlabel('PTM class')
                ax_idr.set_title('Cross-PTM IDR enrichment')
                ax_idr.tick_params(axis='x', rotation=25)
                for x, v, p, q in zip(plot_df['ptm_group'], plot_df['log2_odds_ratio'], plot_df['p_value'], plot_df['q_value']):
                    ax_idr.text(x, v, f'p={format_p_value(p)}\nq={format_p_value(q)}', ha='center', va='bottom' if v >= 0 else 'top', fontsize=7.5, rotation=90)
                fig_idr.tight_layout()
                save_figure(fig_idr, Path(args.outdir) / 'cross_ptm_idr_binary_comparison')

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
    if rows:
        comparison = pd.concat(rows, ignore_index=True)
        save_table(comparison, outdir / 'ptm_context_comparison.tsv')
        save_table(pd.DataFrame(backgrounds), outdir / 'ptm_context_backgrounds.tsv')

        if (comparison['context_layer'] == 'disorder').any():
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
                        ax.text(row['ptm_group'], row['log2_odds_ratio'], f" p={format_p_value(row['p_value'])}, q={format_p_value(row['q_value'])}", fontsize=6.5, va='bottom' if row['log2_odds_ratio'] >= 0 else 'top')
                ax.axhline(0, color='black', linewidth=0.8)
                ax.set_ylabel('log2(OR) vs target-residue proteome background')
                ax.set_xlabel('PTM class')
                ax.set_title('Cross-PTM disorder-context comparison')
                ax.tick_params(axis='x', rotation=25)
                ax.legend()
                fig.tight_layout()
                save_figure(fig, outdir / 'cross_ptm_disorder_context_comparison')

        if (comparison['context_layer'] == 'domain').any():
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
                        ax2.text(row['ptm_group'], row['log2_odds_ratio'], f" p={format_p_value(row['p_value'])}, q={format_p_value(row['q_value'])}", fontsize=6.5, va='bottom' if row['log2_odds_ratio'] >= 0 else 'top')
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

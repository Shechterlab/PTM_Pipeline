from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import fisher_exact

from common import (
    BREWER_COLORS,
    add_q_values,
    annotate_barh,
    apply_paper_style,
    canonicalize_uniprot_accession,
    log2_odds_ratio,
    q_threshold_note,
    save_figure,
    save_table,
    set_symmetric_xlim,
    signed_bar_colors,
    style_axis,
)


def normalize_membership(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    accession_cols = ['canonical_UniProtAC', 'uniprot_acc', 'uniprotkb_ac', 'UniProt', 'Entry', 'protein_accession']
    condensate_id_cols = ['condensate_id', 'CDCODE_ID', 'cdcode_id', 'condensate']
    condensate_name_cols = ['condensate_name', 'name', 'condensate_label', 'CDCODE_name']

    def first_existing(candidates: list[str]) -> str:
        for col in candidates:
            if col in out.columns:
                return col
        return ''

    acc_col = first_existing(accession_cols)
    id_col = first_existing(condensate_id_cols)
    name_col = first_existing(condensate_name_cols)
    if not acc_col or not id_col:
        raise ValueError('staged condensate membership must include an accession column and condensate id column')

    normalized = pd.DataFrame()
    normalized['canonical_UniProtAC'] = out[acc_col].astype(str).map(canonicalize_uniprot_accession)
    normalized['condensate_id'] = out[id_col].astype(str).str.strip()
    normalized['condensate_name'] = out[name_col].astype(str).str.strip() if name_col else normalized['condensate_id']
    normalized = normalized[
        normalized['canonical_UniProtAC'].ne('') &
        normalized['condensate_id'].ne('')
    ].drop_duplicates()
    return normalized


def normalize_protein_species(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    accession_cols = ['canonical_UniProtAC', 'uniprot_acc', 'uniprotkb_ac', 'UniProt', 'Entry', 'protein_accession', 'uniprot_id']
    species_tax_cols = ['species_taxon_id', 'taxon_id']
    species_name_cols = ['species_name', 'organism_name']

    def first_existing(candidates: list[str]) -> str:
        for col in candidates:
            if col in out.columns:
                return col
        return ''

    acc_col = first_existing(accession_cols)
    tax_col = first_existing(species_tax_cols)
    name_col = first_existing(species_name_cols)
    if not acc_col:
        raise ValueError('protein table must include a UniProt accession column')
    normalized = pd.DataFrame()
    normalized['canonical_UniProtAC'] = out[acc_col].astype(str).map(canonicalize_uniprot_accession)
    normalized['species_taxon_id'] = out[tax_col] if tax_col else pd.NA
    normalized['species_name'] = out[name_col].astype(str) if name_col else ''
    return normalized.drop_duplicates()


def normalize_condensate_species(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    id_cols = ['condensate_id', 'CDCODE_ID', 'cdcode_id', 'condensate', 'uid']
    name_cols = ['condensate_name', 'name', 'condensate_label', 'CDCODE_name']
    species_tax_cols = ['species_taxon_id', 'taxon_id']
    species_name_cols = ['species_name', 'organism_name']

    def first_existing(candidates: list[str]) -> str:
        for col in candidates:
            if col in out.columns:
                return col
        return ''

    id_col = first_existing(id_cols)
    if not id_col:
        raise ValueError('condensate table must include a condensate id column')
    normalized = pd.DataFrame()
    normalized['condensate_id'] = out[id_col].astype(str).str.strip()
    name_col = first_existing(name_cols)
    normalized['condensate_name'] = out[name_col].astype(str).str.strip() if name_col else normalized['condensate_id']
    tax_col = first_existing(species_tax_cols)
    normalized['species_taxon_id'] = out[tax_col] if tax_col else pd.NA
    species_name_col = first_existing(species_name_cols)
    normalized['species_name'] = out[species_name_col].astype(str) if species_name_col else ''
    normalized['condensate_type'] = out['condensate_type'].astype(str) if 'condensate_type' in out.columns else ''
    return normalized.drop_duplicates()


def odds_ratio(a: int, b: int, c: int, d: int) -> float:
    return ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--integrated-sites', required=True)
    ap.add_argument('--staged-membership', default='')
    ap.add_argument('--proteins-table', default='')
    ap.add_argument('--condensates-table', default='')
    ap.add_argument('--outdir', default='results/condensates')
    ap.add_argument('--min-condensate-size', type=int, default=5)
    ap.add_argument('--include-synthetic', action='store_true')
    args = ap.parse_args()

    if not args.staged_membership:
        raise ValueError(
            '--staged-membership is required and should point to a protein-to-condensate membership table '
            '(for example protein2cdcode_v2.2.tsv or cdcode_staged_membership.tsv)'
        )

    sites = pd.read_csv(args.integrated_sites, sep='\t', low_memory=False)
    membership_sep = '\t' if str(args.staged_membership).endswith('.tsv') else ','
    memberships = normalize_membership(pd.read_csv(args.staged_membership, sep=membership_sep, low_memory=False))

    if args.proteins_table:
        proteins_sep = '\t' if str(args.proteins_table).endswith('.tsv') else ','
        proteins = normalize_protein_species(pd.read_csv(args.proteins_table, sep=proteins_sep, low_memory=False))
        human_acc = set(
            proteins.loc[
                proteins['species_taxon_id'].astype(str).eq('9606') |
                proteins['species_name'].astype(str).str.contains('Homo sapiens', case=False, na=False),
                'canonical_UniProtAC',
            ].dropna().astype(str)
        )
        memberships = memberships[memberships['canonical_UniProtAC'].isin(human_acc)].copy()

    if args.condensates_table:
        condensates_sep = '\t' if str(args.condensates_table).endswith('.tsv') else ','
        condensates = normalize_condensate_species(pd.read_csv(args.condensates_table, sep=condensates_sep, low_memory=False))
        human_cond = condensates.loc[
            condensates['species_taxon_id'].astype(str).eq('9606') |
            condensates['species_name'].astype(str).str.contains('Homo sapiens', case=False, na=False),
            ['condensate_id', 'condensate_name', 'condensate_type'],
        ].drop_duplicates()
        if not args.include_synthetic:
            human_cond = human_cond[~human_cond['condensate_type'].astype(str).str.contains('synthetic', case=False, na=False)].copy()
        human_cond = human_cond.sort_values(['condensate_id', 'condensate_name']).drop_duplicates(subset=['condensate_id'], keep='first')
        memberships = memberships.merge(
            human_cond,
            how='inner',
            on='condensate_id',
            suffixes=('', '_catalog'),
        )
        if 'condensate_name_catalog' in memberships.columns:
            memberships['condensate_name'] = (
                memberships['condensate_name_catalog']
                .fillna(memberships['condensate_name'])
                .astype(str)
                .str.strip()
                .replace({'': pd.NA})
                .fillna(memberships['condensate_name'])
            )
            memberships = memberships.drop(columns=['condensate_name_catalog'])

    methyl_proteins = set(sites['canonical_UniProtAC'].dropna().astype(str).map(canonicalize_uniprot_accession))
    universe = set(memberships['canonical_UniProtAC'].astype(str))
    methyl_in_universe = methyl_proteins & universe
    global_methyl_fraction = (len(methyl_in_universe) / len(universe)) if universe else 0.0

    rows = []
    for (condensate_id, condensate_name), group in memberships.groupby(['condensate_id', 'condensate_name']):
        members = set(group['canonical_UniProtAC'].astype(str))
        if len(members) < args.min_condensate_size:
            continue
        a = len(members & methyl_proteins)
        b = len(methyl_in_universe) - a
        c = len(members) - a
        d = len(universe - members) - b
        current_or = odds_ratio(a, b, c, d)
        rows.append({
            'condensate_id': condensate_id,
            'condensate_name': condensate_name,
            'member_count': len(members),
            'methyl_protein_count': a,
            'methyl_fraction_within_condensate': a / len(members) if members else 0.0,
            'expected_methyl_count_at_global_rate': len(members) * global_methyl_fraction,
            'odds_ratio': current_or,
            'log2_odds_ratio': log2_odds_ratio(current_or),
            'p_value': float(fisher_exact([[a, b], [c, d]], alternative='two-sided').pvalue),
        })

    enrich = add_q_values(pd.DataFrame(rows)).sort_values(['odds_ratio', 'methyl_protein_count'], ascending=[False, False])
    summary = pd.DataFrame([{
        'methyl_proteins_total': len(methyl_proteins),
        'condensate_universe_proteins': len(universe),
        'methyl_proteins_in_any_condensate': len(methyl_in_universe),
        'fraction_methyl_proteins_in_any_condensate': (len(methyl_in_universe) / len(methyl_proteins)) if methyl_proteins else 0.0,
        'fraction_condensate_universe_proteins_methylated': global_methyl_fraction,
        'condensate_count_after_filters': int(enrich['condensate_id'].nunique()) if not enrich.empty else 0,
        'include_synthetic': bool(args.include_synthetic),
    }])

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(memberships, outdir / 'condensate_membership.tsv')
    save_table(enrich, outdir / 'condensate_enrichment.tsv')
    save_table(summary, outdir / 'condensate_summary.tsv')

    apply_paper_style()
    positive_plot_df = enrich[(enrich['member_count'] >= 10) & (enrich['q_value'] <= 0.05)].sort_values(['log2_odds_ratio', 'methyl_protein_count'], ascending=[False, False]).head(15).sort_values('log2_odds_ratio')
    if not positive_plot_df.empty:
        fig0, ax0 = plt.subplots(figsize=(10.0, 7.2))
        colors0 = signed_bar_colors(positive_plot_df['log2_odds_ratio'])
        ax0.barh(positive_plot_df['condensate_name'], positive_plot_df['log2_odds_ratio'], color=colors0, edgecolor='white', linewidth=0.8)
        style_axis(ax0, zero='x')
        ax0.set_xlabel('log2(OR) vs human condensate protein universe')
        ax0.set_ylabel('Condensate')
        ax0.set_title('Condensate enrichment of methylarginine proteins')
        set_symmetric_xlim(ax0, positive_plot_df['log2_odds_ratio'], annotation_pad_ratio=0.36, center_on_zero=False)
        annotate_barh(
            ax0,
            positive_plot_df['condensate_name'],
            positive_plot_df['log2_odds_ratio'],
            [f'{nm}/{nt}' for nm, nt in zip(positive_plot_df['methyl_protein_count'], positive_plot_df['member_count'])],
            fontsize=8.8,
        )
        fig0.text(0.99, 0.01, q_threshold_note(positive_plot_df['q_value']), ha='right', va='bottom', fontsize=8, color=BREWER_COLORS['dark_gray'])
        fig0.tight_layout()
        save_figure(fig0, outdir / 'condensate_enrichment')

    plot_df = enrich[(enrich['member_count'] >= 10) & (enrich['q_value'] <= 0.05)].sort_values(['log2_odds_ratio', 'methyl_protein_count'], ascending=[False, False]).head(12).sort_values('log2_odds_ratio')
    if not plot_df.empty:
        fig, ax = plt.subplots(figsize=(10.2, 7.6))
        colors = signed_bar_colors(plot_df['log2_odds_ratio'])
        ax.barh(plot_df['condensate_name'], plot_df['log2_odds_ratio'], color=colors, edgecolor='white', linewidth=0.8)
        style_axis(ax, zero='x')
        ax.set_xlabel('log2(OR) vs human condensate protein universe')
        ax.set_ylabel('Condensate')
        ax.set_title('Condensate enrichment of methylarginine proteins')
        set_symmetric_xlim(ax, plot_df['log2_odds_ratio'], annotation_pad_ratio=0.38, center_on_zero=False)
        annotate_barh(
            ax,
            plot_df['condensate_name'],
            plot_df['log2_odds_ratio'],
            [f'{nm}/{nt}' for nm, nt in zip(plot_df['methyl_protein_count'], plot_df['member_count'])],
            fontsize=8.8,
        )
        fig.text(0.99, 0.01, q_threshold_note(plot_df['q_value']), ha='right', va='bottom', fontsize=8, color=BREWER_COLORS['dark_gray'])
        fig.tight_layout()
        save_figure(fig, outdir / 'condensate_enrichment_top_bottom')

        fig2, ax2 = plt.subplots(figsize=(10.2, 7.6))
        frac_pct = plot_df['methyl_fraction_within_condensate'] * 100.0
        baseline_pct = global_methyl_fraction * 100.0
        colors2 = [
            BREWER_COLORS['dark_gray'] if frac >= baseline_pct else BREWER_COLORS['mid_gray']
            for frac in frac_pct
        ]
        ax2.barh(plot_df['condensate_name'], frac_pct, color=colors2, edgecolor='white', linewidth=0.8)
        ax2.axvline(baseline_pct, color=BREWER_COLORS['dark_gray'], linewidth=1.0, linestyle='--')
        style_axis(ax2)
        ax2.set_xlabel('% condensate proteins in arg-methylome')
        ax2.set_ylabel('Condensate')
        ax2.set_title('Condensate methyl-protein burden')
        for y, frac, nm, nt in zip(plot_df['condensate_name'], frac_pct, plot_df['methyl_protein_count'], plot_df['member_count']):
            ax2.text(frac, y, f'  {nm}/{nt}', va='center', ha='left', fontsize=8.0)
        fig2.text(
            0.99,
            0.01,
            f'Dashed line: global condensate-universe baseline = {baseline_pct:.1f}%',
            ha='right',
            va='bottom',
            fontsize=8,
            color=BREWER_COLORS['dark_gray'],
        )
        fig2.tight_layout()
        save_figure(fig2, outdir / 'condensate_fraction_methylated_top_bottom')


if __name__ == '__main__':
    main()

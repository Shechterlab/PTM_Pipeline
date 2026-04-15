from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from common import (
    BREWER_COLORS,
    add_q_values,
    apply_paper_style,
    canonicalize_uniprot_accession,
    normalize_accession,
    parse_fasta,
    save_figure,
    save_table,
    sequence_window,
    style_axis,
)


MOD_RE = re.compile(r'R(\d+)me([12])')
HYDROPHOBIC = set('AILMFWVY')


def residue_at(window: str, relative_position: int) -> str:
    center = len(window) // 2
    index = center + relative_position
    if index < 0 or index >= len(window):
        return '_'
    return window[index]


def family_flags(window: str) -> dict[str, bool]:
    neighborhood = [residue_at(window, i) for i in range(-5, 6) if i != 0]
    near3 = [residue_at(window, i) for i in range(-3, 4) if i != 0]
    return {
        'RG': residue_at(window, 1) == 'G',
        'RGG': residue_at(window, 1) == 'G' and residue_at(window, 2) == 'G',
        'GRG': residue_at(window, -1) == 'G' and residue_at(window, 1) == 'G',
        'GAR': residue_at(window, -2) == 'G' and residue_at(window, -1) == 'A',
        'RXR': residue_at(window, -2) == 'R' or residue_at(window, 2) == 'R',
        'DR': residue_at(window, -1) == 'D',
        'ER': residue_at(window, -1) == 'E',
        'RD': residue_at(window, 1) == 'D',
        'RE': residue_at(window, 1) == 'E',
        'CARM1_proline_rich': 'P' in neighborhood,
        'CARM1_hydrophobic_proline': ('P' in neighborhood) and any(aa in HYDROPHOBIC for aa in near3),
    }


def first_nonempty(series: pd.Series):
    for value in series:
        if pd.notna(value) and str(value).strip() and str(value).strip().lower() != 'nan':
            return value
    return pd.NA


def site_key_frame(df: pd.DataFrame) -> pd.Series:
    return df['canonical_UniProtAC'].astype(str) + ':' + df['position'].astype(int).astype(str)


def odds_ratio(a: int, b: int, c: int, d: int) -> float:
    return ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))


def enrich_binary(target: pd.Series, target_label: str, background: pd.Series, background_label: str) -> pd.DataFrame:
    rows = []
    values = sorted(set(target.astype(str)).union(set(background.astype(str))))
    for value in values:
        a = int((target.astype(str) == value).sum())
        b = int(len(target) - a)
        c = int((background.astype(str) == value).sum())
        d = int(len(background) - c)
        rows.append({
            'category': value,
            'target_label': target_label,
            'background_label': background_label,
            'target_count': a,
            'target_fraction': a / len(target) if len(target) else np.nan,
            'background_count': c,
            'background_fraction': c / len(background) if len(background) else np.nan,
            'odds_ratio': odds_ratio(a, b, c, d),
            'log2_odds_ratio': np.log2(odds_ratio(a, b, c, d)),
            'p_value': float(fisher_exact([[a, b], [c, d]], alternative='two-sided').pvalue),
        })
    return add_q_values(pd.DataFrame(rows))


def parse_modifications(mod_text: str) -> list[tuple[int, str]]:
    out = []
    for position, degree in MOD_RE.findall(str(mod_text)):
        out.append((int(position), f'me{degree}'))
    return out


def normalize_acc(accession: str) -> str:
    return canonicalize_uniprot_accession(normalize_accession(str(accession).strip()))


def attach_window(df: pd.DataFrame, sequences: dict[str, str], flank: int = 7) -> pd.DataFrame:
    out = df.copy()
    out['canonical_UniProtAC'] = out['accession'].map(normalize_acc)
    out['sequence_window'] = ''
    out['sequence_validated'] = False
    for idx, row in out.iterrows():
        acc = str(row['canonical_UniProtAC'])
        pos = int(row['position'])
        seq = sequences.get(acc, '')
        if not seq or pos < 1 or pos > len(seq) or seq[pos - 1] != 'R':
            continue
        out.at[idx, 'sequence_window'] = sequence_window(seq, pos, flank)
        out.at[idx, 'sequence_validated'] = True
    return out


def collapse_s2_site_state(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    keys = ['canonical_UniProtAC', 'position', 'methyl_state']
    base = (
        df.sort_values(['sequence_validated'], ascending=[False])
        .drop_duplicates(keys)
        .copy()
    )
    obs = df.groupby(keys, dropna=False).size().rename('observation_rows').reset_index()
    base = base.merge(obs, on=keys, how='left')
    base['site_key'] = site_key_frame(base)
    return base


def collapse_s3_site_treatment(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    grouped = (
        df.groupby(['canonical_UniProtAC', 'gene', 'position', 'ip', 'treatment'], dropna=False)
        .agg(
            accession=('accession', first_nonempty),
            peptide_rows=('position', 'size'),
            sequence_window=('sequence_window', first_nonempty),
            sequence_validated=('sequence_validated', 'max'),
            me1_row_fraction=('methyl_state', lambda x: (x == 'me1').mean()),
            me2_row_fraction=('methyl_state', lambda x: (x == 'me2').mean()),
            methyl_states_observed=('methyl_state', lambda x: ';'.join(sorted(set(str(v) for v in x if pd.notna(v) and str(v).strip())))),
            treatment_diff_median=('treatment_diff', 'median'),
            treatment_diff_mean=('treatment_diff', 'mean'),
            treatment_padj_min=('treatment_padj', 'min'),
            treatment_pval_min=('treatment_pval', 'min'),
            peptide_sequence_example=('peptide_sequence', first_nonempty),
        )
        .reset_index()
    )
    grouped['site_key'] = site_key_frame(grouped)
    return grouped


def collapse_s3_site_ip(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    grouped = (
        df.groupby(['canonical_UniProtAC', 'gene', 'position', 'ip'], dropna=False)
        .agg(
            accession=('accession', first_nonempty),
            treatment_count=('treatment', 'nunique'),
            treatments=('treatment', lambda x: ';'.join(sorted(set(str(v) for v in x if pd.notna(v) and str(v).strip())))),
            peptide_rows_total=('peptide_rows', 'sum'),
            sequence_window=('sequence_window', first_nonempty),
            sequence_validated=('sequence_validated', 'max'),
            me1_row_fraction=('me1_row_fraction', 'mean'),
            me2_row_fraction=('me2_row_fraction', 'mean'),
        )
        .reset_index()
    )
    grouped['site_key'] = site_key_frame(grouped)
    return grouped


def s2_sites(path: str | Path) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name=0)
    out = pd.DataFrame({
        'accession': df['Accession'].astype(str).str.strip(),
        'gene': df['Gene.Symbol'].astype(str).str.strip(),
        'position': pd.to_numeric(df['Position.In.Protein'], errors='coerce'),
        'modification': df['Modification'].astype(str).str.strip(),
        'disordered': df['Disordered'].astype(str).map({'True': True, 'False': False}),
        'peptide_sequence': df['Sequence'].astype(str).str.strip(),
        'table': 'Table_S2',
    })
    out = out[out['position'].notna()].copy()
    out['position'] = out['position'].astype(int)
    out['methyl_state'] = out['modification'].str.extract(r'(me[12])', expand=False).fillna('')
    return out


def s3_sites(path: str | Path) -> pd.DataFrame:
    rows = []
    for sheet, treatment in [('Control_vs_GSK591', 'GSK591'), ('Control_vs_MS023', 'MS023')]:
        df = pd.read_excel(path, sheet_name=sheet)
        diff_col = f'{treatment}_vs_DMSO_diff'
        padj_col = f'{treatment}_vs_DMSO_p.adj'
        pval_col = f'{treatment}_vs_DMSO_p.val'
        for _, row in df.iterrows():
            mods = parse_modifications(row.get('Modifications', ''))
            if not mods:
                continue
            for pos, state in mods:
                rows.append({
                    'accession': str(row.get('Accession', '')).strip(),
                    'gene': str(row.get('Gene.Symbol', '')).strip(),
                    'peptide_sequence': str(row.get('Sequence', '')).strip(),
                    'position': pos,
                    'methyl_state': state,
                    'ip': str(row.get('IP', '')).strip(),
                    'treatment': treatment,
                    'treatment_diff': pd.to_numeric(row.get(diff_col), errors='coerce'),
                    'treatment_padj': pd.to_numeric(row.get(padj_col), errors='coerce'),
                    'treatment_pval': pd.to_numeric(row.get(pval_col), errors='coerce'),
                    'begin_pos': pd.to_numeric(row.get('Begin.Pos'), errors='coerce'),
                    'end_pos': pd.to_numeric(row.get('End.Pos'), errors='coerce'),
                    'modifications': str(row.get('Modifications', '')).strip(),
                    'sheet': sheet,
                    'table': 'Table_S3',
                })
    out = pd.DataFrame(rows)
    out = out[out['position'].notna()].copy()
    out['position'] = out['position'].astype(int)
    return out


def motif_by_group(df: pd.DataFrame, group_col: str) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()
    motif_cols = sorted(family_flags(df['sequence_window'].iloc[0]).keys())
    flags = pd.DataFrame([family_flags(w) for w in df['sequence_window']])
    data = pd.concat([df.reset_index(drop=True), flags], axis=1)
    rows = []
    groups = sorted(set(data[group_col].dropna().astype(str)))
    for group in groups:
        target = data[data[group_col].astype(str) == group]
        background = data[data[group_col].astype(str) != group]
        if target.empty or background.empty:
            continue
        for motif in motif_cols:
            a = int(target[motif].sum())
            b = int(len(target) - a)
            c = int(background[motif].sum())
            d = int(len(background) - c)
            rows.append({
                'group_col': group_col,
                'group_value': group,
                'motif_family': motif,
                'target_count': a,
                'target_fraction': a / len(target),
                'background_count': c,
                'background_fraction': c / len(background),
                'odds_ratio': odds_ratio(a, b, c, d),
                'log2_odds_ratio': np.log2(odds_ratio(a, b, c, d)),
                'p_value': float(fisher_exact([[a, b], [c, d]], alternative='two-sided').pvalue),
            })
    return add_q_values(pd.DataFrame(rows)) if rows else pd.DataFrame()


def treatment_summary(df: pd.DataFrame, padj_cutoff: float) -> pd.DataFrame:
    out = df.copy()
    out['sig_down'] = out['treatment_padj_min'].notna() & (out['treatment_padj_min'] <= padj_cutoff) & (out['treatment_diff_median'] < 0)
    out['sig_up'] = out['treatment_padj_min'].notna() & (out['treatment_padj_min'] <= padj_cutoff) & (out['treatment_diff_median'] > 0)
    rows = []
    for (treatment, ip), group in out.groupby(['treatment', 'ip'], dropna=False):
        rows.append({
            'treatment': treatment,
            'ip': ip,
            'unique_site_ip_treatment': len(group),
            'unique_canonical_sites': group['site_key'].nunique(),
            'peptide_rows_total': int(group['peptide_rows'].sum()),
            'median_diff': group['treatment_diff_median'].median(),
            'mean_diff': group['treatment_diff_mean'].mean(),
            'sig_down_sites': int(group['sig_down'].sum()),
            'sig_up_sites': int(group['sig_up'].sum()),
            'fraction_sig_down_sites': float(group['sig_down'].mean()),
            'fraction_sig_up_sites': float(group['sig_up'].mean()),
        })
    return pd.DataFrame(rows).sort_values(['treatment', 'fraction_sig_down_sites'], ascending=[True, False])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--table-s2', required=True)
    ap.add_argument('--table-s3', required=True)
    ap.add_argument('--canonical-fasta', required=True)
    ap.add_argument('--outdir', default='PTM_results/maron_state_treatment')
    ap.add_argument('--padj-cutoff', type=float, default=0.05)
    args = ap.parse_args()

    apply_paper_style()
    sequences = parse_fasta(args.canonical_fasta)

    s2 = attach_window(s2_sites(args.table_s2), sequences)
    s3 = attach_window(s3_sites(args.table_s3), sequences)
    s2_site_state = collapse_s2_site_state(s2)
    s2_valid = s2_site_state[s2_site_state['sequence_validated']].copy()
    s3_site_treatment = collapse_s3_site_treatment(s3)
    s3_site_ip = collapse_s3_site_ip(s3_site_treatment)
    s3_valid = s3_site_ip[s3_site_ip['sequence_validated']].copy()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(s2, outdir / 'maron_table_s2_sites.tsv')
    save_table(s3, outdir / 'maron_table_s3_sites.tsv')
    save_table(s2_site_state, outdir / 'maron_table_s2_site_state_collapsed.tsv')
    save_table(s3_site_treatment, outdir / 'maron_table_s3_site_treatment_collapsed.tsv')
    save_table(s3_site_ip, outdir / 'maron_table_s3_site_ip_collapsed.tsv')

    qc = pd.DataFrame([{
        's2_state_observation_rows': len(s2),
        's2_unique_site_state_observations': len(s2_site_state),
        's2_validated_unique_site_state_observations': len(s2_valid),
        's2_unique_canonical_sites': s2_site_state['site_key'].nunique(),
        's2_sites_observed_in_both_me1_me2': int((s2_site_state.groupby('site_key')['methyl_state'].nunique() > 1).sum()),
        's3_peptide_rows': len(s3),
        's3_unique_site_ip_treatment': len(s3_site_treatment),
        's3_unique_site_ip': len(s3_site_ip),
        's3_validated_unique_site_ip': len(s3_valid),
        's3_unique_canonical_sites': s3_site_ip['site_key'].nunique(),
        's3_sites_observed_in_multiple_ip_classes': int((s3_site_ip.groupby('site_key')['ip'].nunique() > 1).sum()),
        's3_sites_with_multiple_states_within_ip': int((s3.groupby(['canonical_UniProtAC', 'position', 'ip'])['methyl_state'].nunique() > 1).sum()),
        'padj_cutoff': args.padj_cutoff,
    }])
    save_table(qc, outdir / 'maron_state_treatment_qc.tsv')

    s2_state_summary = (
        s2_site_state.groupby('methyl_state', dropna=False)
        .agg(
            unique_site_state_observations=('methyl_state', 'size'),
            proteins=('canonical_UniProtAC', 'nunique'),
            unique_canonical_sites=('site_key', 'nunique'),
            observation_rows=('observation_rows', 'sum'),
            fraction_disordered=('disordered', 'mean'),
        )
        .reset_index()
    )
    save_table(s2_state_summary, outdir / 'maron_s2_state_summary.tsv')
    save_table(
        s2_site_state.groupby('site_key')['methyl_state'].nunique().value_counts().sort_index().rename_axis('states_per_site').reset_index(name='site_count'),
        outdir / 'maron_s2_state_overlap_summary.tsv',
    )

    if not s2_valid.empty:
        s2_state_motifs = motif_by_group(s2_valid, 'methyl_state')
        save_table(s2_state_motifs, outdir / 'maron_s2_state_motif_enrichment.tsv')

    s3_ip_summary = (
        s3_site_ip.groupby('ip', dropna=False)
        .agg(
            unique_site_ip=('ip', 'size'),
            proteins=('canonical_UniProtAC', 'nunique'),
            unique_canonical_sites=('site_key', 'nunique'),
            peptide_rows_total=('peptide_rows_total', 'sum'),
            mean_me1_row_fraction=('me1_row_fraction', 'mean'),
            mean_me2_row_fraction=('me2_row_fraction', 'mean'),
        )
        .reset_index()
        .sort_values('unique_site_ip', ascending=False)
    )
    save_table(s3_ip_summary, outdir / 'maron_s3_ip_summary.tsv')
    save_table(
        s3_site_ip.groupby('site_key')['ip'].nunique().value_counts().sort_index().rename_axis('ip_classes_per_site').reset_index(name='site_count'),
        outdir / 'maron_s3_ip_overlap_summary.tsv',
    )

    s3_treatment = treatment_summary(s3_site_treatment, args.padj_cutoff)
    save_table(s3_treatment, outdir / 'maron_s3_treatment_ip_summary.tsv')

    if not s3_valid.empty:
        s3_ip_motifs = motif_by_group(s3_valid, 'ip')
        save_table(s3_ip_motifs, outdir / 'maron_s3_ip_motif_enrichment.tsv')
        adma_sdma = s3_valid[s3_valid['ip'].isin(['ADMA', 'SDMA'])].copy()
        if not adma_sdma.empty:
            save_table(motif_by_group(adma_sdma, 'ip'), outdir / 'maron_s3_adma_sdma_motif_enrichment.tsv')

    # Site-level inhibitor response figure by IP.
    if not s3_treatment.empty:
        plot_df = s3_treatment.copy()
        plot_df['label'] = plot_df['ip'] + ' | ' + plot_df['treatment']
        fig, ax = plt.subplots(figsize=(7.8, 4.8))
        ax.bar(
            plot_df['label'],
            plot_df['fraction_sig_down_sites'],
            color=[BREWER_COLORS['orange'] if t == 'GSK591' else BREWER_COLORS['blue'] for t in plot_df['treatment']],
            edgecolor='white',
            linewidth=0.8,
        )
        style_axis(ax, grid_axis='y')
        ax.set_ylabel(f'Fraction of unique sites significantly decreased (adj p<={args.padj_cutoff:g})')
        ax.set_xlabel('IP and treatment contrast')
        ax.set_title('Maron site-level depletion by inhibitor and methyl-Arg IP')
        ax.tick_params(axis='x', rotation=25)
        fig.tight_layout()
        save_figure(fig, outdir / 'maron_s3_treatment_ip_depletion')

    # Direct ADMA vs SDMA motif panel if present.
    if outdir.joinpath('maron_s3_ip_motif_enrichment.tsv').exists():
        motif_df = pd.read_csv(outdir / 'maron_s3_ip_motif_enrichment.tsv', sep='\t', low_memory=False)
        motif_df = motif_df[motif_df['group_value'].isin(['ADMA', 'SDMA']) & motif_df['motif_family'].isin(['RG', 'RGG', 'GRG', 'DR', 'ER', 'RD', 'RE'])].copy()
        if not motif_df.empty:
            plot_df = motif_df.sort_values(['group_value', 'log2_odds_ratio'])
            fig2, ax2 = plt.subplots(figsize=(8.4, 5.4))
            xpos = np.arange(len(plot_df['motif_family'].unique()))
            motif_order = ['RG', 'RGG', 'GRG', 'DR', 'ER', 'RD', 'RE']
            width = 0.38
            adma = plot_df[plot_df['group_value'] == 'ADMA'].set_index('motif_family').reindex(motif_order)
            sdma = plot_df[plot_df['group_value'] == 'SDMA'].set_index('motif_family').reindex(motif_order)
            ax2.bar(xpos - width / 2, adma['log2_odds_ratio'], width=width, color=BREWER_COLORS['orange'], edgecolor='white', linewidth=0.8, label='ADMA IP')
            ax2.bar(xpos + width / 2, sdma['log2_odds_ratio'], width=width, color=BREWER_COLORS['blue'], edgecolor='white', linewidth=0.8, label='SDMA IP')
            style_axis(ax2, zero='y', grid_axis='y')
            ax2.set_xticks(xpos)
            ax2.set_xticklabels(motif_order)
            ax2.set_ylabel('log2(OR) vs other IP classes')
            ax2.set_xlabel('Motif family')
            ax2.set_title('Maron IP-specific motif enrichment')
            ax2.legend(frameon=False)
            fig2.tight_layout()
            save_figure(fig2, outdir / 'maron_s3_adma_sdma_motif_comparison')

    print({
        's2_state_observation_rows': len(s2),
        's2_unique_site_state_observations': len(s2_site_state),
        's2_validated_unique_site_state_observations': len(s2_valid),
        's3_peptide_rows': len(s3),
        's3_unique_site_ip_treatment': len(s3_site_treatment),
        's3_validated_unique_site_ip': len(s3_valid),
        'output': str(outdir),
    })


if __name__ == '__main__':
    main()

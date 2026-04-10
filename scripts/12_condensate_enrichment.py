from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import fisher_exact

from common import add_q_values, canonicalize_uniprot_accession, format_p_value, log2_odds_ratio, save_figure, save_table


def normalize_membership(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    accession_cols = ['canonical_UniProtAC', 'uniprot_acc', 'UniProt', 'Entry', 'protein_accession']
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


def odds_ratio(a: int, b: int, c: int, d: int) -> float:
    return ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--integrated-sites', required=True)
    ap.add_argument('--staged-membership', default='')
    ap.add_argument('--outdir', default='results/condensates')
    ap.add_argument('--min-condensate-size', type=int, default=5)
    args = ap.parse_args()

    if not args.staged_membership:
        raise ValueError('only --staged-membership is supported in this version')

    sites = pd.read_csv(args.integrated_sites, sep='\t', low_memory=False)
    memberships = normalize_membership(pd.read_csv(args.staged_membership, sep='\t', low_memory=False))
    methyl_proteins = set(sites['canonical_UniProtAC'].dropna().astype(str).map(canonicalize_uniprot_accession))
    universe = set(memberships['canonical_UniProtAC'].astype(str))

    rows = []
    for (condensate_id, condensate_name), group in memberships.groupby(['condensate_id', 'condensate_name']):
        members = set(group['canonical_UniProtAC'].astype(str))
        if len(members) < args.min_condensate_size:
            continue
        a = len(members & methyl_proteins)
        b = len(methyl_proteins & universe) - a
        c = len(members)
        d = len(universe - members)
        rows.append({
            'condensate_id': condensate_id,
            'condensate_name': condensate_name,
            'member_count': len(members),
            'methyl_protein_count': a,
            'odds_ratio': odds_ratio(a, b, c, d),
            'log2_odds_ratio': log2_odds_ratio(odds_ratio(a, b, c, d)),
            'p_value': float(fisher_exact([[a, b], [c, d]], alternative='two-sided').pvalue),
        })

    enrich = add_q_values(pd.DataFrame(rows)).sort_values(['odds_ratio', 'methyl_protein_count'], ascending=[False, False])

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(memberships, outdir / 'condensate_membership.tsv')
    save_table(enrich, outdir / 'condensate_enrichment.tsv')

    plot_df = enrich.head(15).sort_values('odds_ratio')
    if not plot_df.empty:
        fig, ax = plt.subplots(figsize=(9.5, 6.4))
        ax.barh(plot_df['condensate_name'], plot_df['log2_odds_ratio'])
        ax.axvline(0, color='black', linewidth=0.8)
        ax.set_xlabel('log2(OR) vs condensate protein universe')
        ax.set_ylabel('Condensate')
        ax.set_title('Condensate enrichment of methylarginine proteins')
        for y, v, n, p, q in zip(plot_df['condensate_name'], plot_df['log2_odds_ratio'], plot_df['methyl_protein_count'], plot_df['p_value'], plot_df['q_value']):
            ax.text(v, y, f'  n={n}, p={format_p_value(p)}, q={format_p_value(q)}', va='center', ha='left' if v >= 0 else 'right', fontsize=8)
        fig.tight_layout()
        save_figure(fig, outdir / 'condensate_enrichment')


if __name__ == '__main__':
    main()

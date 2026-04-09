from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import save_figure, save_table



def split_sources(text: str) -> list[str]:
    return [token for token in str(text).split(';') if token and token != 'nan']



def source_combo_table(dedup: pd.DataFrame) -> pd.DataFrame:
    combo = dedup.copy()
    combo['source_combo'] = combo['source_families_exact'].map(lambda x: ';'.join(sorted(split_sources(x))))
    out = combo['source_combo'].value_counts().rename_axis('source_combo').reset_index(name='site_count')
    out['source_family_count'] = out['source_combo'].map(lambda x: len(split_sources(x)))
    return out.sort_values(['source_family_count', 'site_count'], ascending=[False, False])



def overlap_matrix(dedup: pd.DataFrame) -> pd.DataFrame:
    families = sorted({fam for text in dedup['source_families_exact'].fillna('') for fam in split_sources(text)})
    rows = []
    for fam_a in families:
        for fam_b in families:
            mask = dedup['source_families_exact'].fillna('').map(lambda text: fam_a in split_sources(text) and fam_b in split_sources(text))
            rows.append({'source_a': fam_a, 'source_b': fam_b, 'site_count': int(mask.sum())})
    return pd.DataFrame(rows)



def plot_overlap_heatmap(overlap: pd.DataFrame, outpath: Path) -> None:
    pivot = overlap.pivot(index='source_a', columns='source_b', values='site_count').fillna(0).astype(int)
    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    im = ax.imshow(pivot.to_numpy(), aspect='auto')
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=30, ha='right')
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    ax.set_title('Exact site overlap across methyl-Arg source families')
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label('Shared exact sites')
    for i in range(len(pivot.index)):
        for j in range(len(pivot.columns)):
            ax.text(j, i, str(int(pivot.iloc[i, j])), ha='center', va='center', fontsize=8)
    fig.tight_layout()
    save_figure(fig, outpath)



def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--all-rows', required=True)
    ap.add_argument('--dedup-sites', required=True)
    ap.add_argument('--outdir', default='results/review_audit')
    args = ap.parse_args()

    all_rows = pd.read_csv(args.all_rows, sep='\t', low_memory=False)
    dedup = pd.read_csv(args.dedup_sites, sep='\t', low_memory=False)

    combo = source_combo_table(dedup)
    overlap = overlap_matrix(dedup)
    source_counts = all_rows['source_family'].value_counts().rename_axis('source_family').reset_index(name='row_count')
    protein_counts = all_rows.groupby('source_family')['substrate_UniProtAC'].nunique().rename('protein_count').reset_index()
    source_counts = source_counts.merge(protein_counts, how='left', on='source_family')
    support = dedup['source_family_count_exact'].value_counts().rename_axis('source_family_count_exact').reset_index(name='site_count')
<<<<<<< HEAD
<<<<<<< HEAD
=======
    confidence = dedup['confidence_tier'].value_counts().rename_axis('confidence_tier').reset_index(name='site_count') if 'confidence_tier' in dedup.columns else pd.DataFrame()
>>>>>>> 2c1a8ff7d8603628918f8ffc773cdcc98c903bff
=======
>>>>>>> 4fbe8b2dc495fcf875e4e868ffb2aca4033251f1
    top_proteins = dedup.groupby('substrate_UniProtAC').size().rename('site_count').reset_index().sort_values('site_count', ascending=False).head(25)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(combo, outdir / 'source_combo_counts.tsv')
    save_table(overlap, outdir / 'source_overlap_matrix.tsv')
    save_table(source_counts, outdir / 'source_row_and_protein_counts.tsv')
    save_table(support, outdir / 'site_support_distribution.tsv')
    save_table(top_proteins, outdir / 'top_proteins_by_union_site_count.tsv')
<<<<<<< HEAD
<<<<<<< HEAD
=======
    if not confidence.empty:
        save_table(confidence, outdir / 'confidence_tier_distribution.tsv')
>>>>>>> 2c1a8ff7d8603628918f8ffc773cdcc98c903bff
=======
>>>>>>> 4fbe8b2dc495fcf875e4e868ffb2aca4033251f1

    fig1, ax1 = plt.subplots(figsize=(8.8, 5.2))
    plot_combo = combo.head(12).sort_values(['source_family_count', 'site_count'])
    ax1.barh(plot_combo['source_combo'], plot_combo['site_count'])
    ax1.set_xlabel('Deduplicated methyl-Arg sites')
    ax1.set_ylabel('Exact source combination')
    ax1.set_title('Largest exact-support source combinations in the union methylarginiome')
    fig1.tight_layout()
    save_figure(fig1, outdir / 'source_combo_barplot')

    plot_overlap_heatmap(overlap, outdir / 'source_overlap_heatmap')

    fig2, ax2 = plt.subplots(figsize=(6.4, 4.6))
    ax2.bar(support['source_family_count_exact'].astype(str), support['site_count'])
    ax2.set_xlabel('Exact supporting source families')
    ax2.set_ylabel('Deduplicated methyl-Arg sites')
    ax2.set_title('Exact support distribution across source families')
    fig2.tight_layout()
    save_figure(fig2, outdir / 'site_support_distribution')

<<<<<<< HEAD
<<<<<<< HEAD

=======
    if not confidence.empty:
        fig3, ax3 = plt.subplots(figsize=(8.2, 4.8))
        ax3.barh(confidence['confidence_tier'], confidence['site_count'])
        ax3.set_xlabel('Deduplicated methyl-Arg sites')
        ax3.set_ylabel('Confidence tier')
        ax3.set_title('Cross-resource recurrence tier distribution')
        fig3.tight_layout()
        save_figure(fig3, outdir / 'confidence_tier_distribution')
>>>>>>> 2c1a8ff7d8603628918f8ffc773cdcc98c903bff
=======

>>>>>>> 4fbe8b2dc495fcf875e4e868ffb2aca4033251f1

    print({
        'all_rows': int(len(all_rows)),
        'dedup_sites': int(len(dedup)),
        'source_families': sorted(set(all_rows['source_family'].astype(str))),
        'output': str(outdir),
    })


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Do methylated lysines sit near methylated arginines more often than the other lysines of the same proteins?

Every lysine of every canonical protein in the methylarginine union is one observation: the outcome is whether
it is methylated (Lys methylation entries of the Site Explorer export: the cross-PTM master table and the
sequence-validated CPLM 4.0 sites from 39_stage_cplm4_lys_methylation.py), the exposure whether a
methylarginine of the union lies within +/-W residues. Odds ratios are Mantel-Haenszel estimates over strata of

  protein                                   (each lysine compared with lysines of its own protein)
  protein x disorder                        (MobiDB consensus disorder, prediction-disorder-th_50)
  protein x disorder x arginines within W   (lysines with the same number of nearby arginines, so that
                                             "near a methylarginine" is not just "in an arginine-rich stretch")

with the CMH test of no association. The distance profile compares lysines whose nearest methylarginine is
1-5, 6-10, 11-20 or 21-50 residues away with lysines more than 50 residues from any (or in proteins whose
methylarginines are all further), within protein x disorder: enrichment confined to the first bins is the
scale of one proteolytic peptide, where co-detection of both marks on the same peptide can inflate it.

Both marks are mostly from mass spectrometry, so a study that reports both on one peptide inflates co-occurrence.
With --cplm4-sites, one more comparison keeps only "independent" methyl-K: CPLM 4.0 experimental sites (evidence
Exp.) none of whose references is a study or database that also contributed methylarginines to the union
(UNION_STUDY_PMIDS, including PhosphoSitePlus, which feeds both).

Inputs: the explorer export (docs/data/methylarg_sites.json), the canonical FASTA its positions refer to, and
disorder intervals (canonical_UniProtAC, fragment_start, fragment_end; PTM_MOBIDB_INTERVALS). Proteins whose
FASTA length differs from the export's are skipped, and so are methylarginines without an R and Lys-methyl
entries without a K in the FASTA.
"""
from __future__ import annotations

import argparse
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2
from statsmodels.stats.contingency_tables import StratifiedTable

from common import BREWER_COLORS, apply_paper_style, cmh_statistic, parse_fasta, save_figure, save_table, style_axis

KME_GROUP = 'Lys methylation'
# Same histone families as the explorer's histone numbering (core and linker histones and their variants).
HISTONE_GENE = re.compile(r'^(H1-\d+|H2A[A-Z0-9-]*|H2B[A-Z0-9-]*|H3C\d+|H3-\d+[A-Z]?|H3Y\d|H4C\d+|H4-\d+|MACROH2A\d|CENPA)$')
# PubMed IDs of the studies and databases behind the methylarginine union (from their DOIs, PubMed esearch).
UNION_STUDY_PMIDS = {
    'Maron 2021': '34505004', 'ProMetheusDB': '35577067', 'PhosphoSitePlus': '25514926', 'Larsen 2016': '27577262',
    'Musiani 2019': '30940768', 'Fong 2019': '31408619', 'Guo 2014': '24129315', 'Fedoriw 2019': '31257072',
    'Wei 2020': '33344439', 'Lim 2020': '32467672', 'Li 2021': '33782401',
}
DIST_BINS = [(1, 5), (6, 10), (11, 20), (21, 50)]
FAR = 50


def r_bin(n: int) -> str:
    return '0' if n == 0 else '1' if n == 1 else '2' if n == 2 else '3-4' if n <= 4 else '5+'


def load_disorder(path: Path) -> dict[str, list[tuple[int, int]]]:
    df = pd.read_csv(path, sep='\t')
    df = df.rename(columns={'accession': 'canonical_UniProtAC', 'Entry': 'canonical_UniProtAC', 'start': 'fragment_start', 'end': 'fragment_end'})
    out: dict[str, list[tuple[int, int]]] = {}
    for acc, a, b in df[['canonical_UniProtAC', 'fragment_start', 'fragment_end']].itertuples(index=False):
        out.setdefault(str(acc), []).append((int(a), int(b)))
    return out


def lysine_table(export: dict, seqs: dict[str, str], disorder: dict, window: int) -> tuple[pd.DataFrame, dict]:
    meta = export['meta']
    pf = meta['protein_fields']
    i_acc, i_gene, i_len, i_ptm = (pf.index(k) for k in ('accession', 'gene', 'length', 'other_ptms'))
    sf = meta['site_fields']
    s_prot, s_pos, s_mis = (sf.index(k) for k in ('protein', 'position', 'window_center_mismatch'))
    group = meta['other_ptm_groups'].index(KME_GROUP)
    methyl_r: dict[int, list[int]] = {}
    for s in export['sites']:
        if not s[s_mis]:
            methyl_r.setdefault(s[s_prot], []).append(s[s_pos])
    rows = []
    qc = {'proteins_in_export': len(export['proteins']), 'skipped_isoform': 0, 'skipped_no_sequence': 0, 'skipped_length_mismatch': 0,
          'methylR_not_R': 0, 'methylK_not_K': 0, 'proteins_without_disorder_record': 0}
    for i, prot in enumerate(export['proteins']):
        acc = prot[i_acc]
        if '-' in acc:
            qc['skipped_isoform'] += 1
            continue
        seq = seqs.get(acc)
        if not seq:
            qc['skipped_no_sequence'] += 1
            continue
        if prot[i_len] and prot[i_len] != len(seq):
            qc['skipped_length_mismatch'] += 1
            continue
        mr = sorted({p for p in methyl_r.get(i, []) if p <= len(seq) and seq[p - 1] == 'R'})
        qc['methylR_not_R'] += len(set(methyl_r.get(i, []))) - len(mr)
        if not mr:
            continue
        kme: dict[int, int] = {}
        for item in prot[i_ptm]:
            if item[1] != group:
                continue
            pos, mask = item[0], (item[3] if len(item) > 3 else 1)
            if pos <= len(seq) and seq[pos - 1] == 'K':
                kme[pos] = kme.get(pos, 0) | mask
            else:
                qc['methylK_not_K'] += 1
        dis = disorder.get(acc)
        if dis is None:
            qc['proteins_without_disorder_record'] += 1
        disordered = np.zeros(len(seq) + 2, dtype=bool)
        for a, b in dis or []:
            disordered[max(1, a):min(len(seq), b) + 1] = True
        mr_arr = np.array(mr)
        gene = prot[i_gene] or ''
        for pos in (k + 1 for k, c in enumerate(seq) if c == 'K'):
            d = np.abs(mr_arr - pos)
            rows.append((acc, gene, pos, kme.get(pos, 0), int(d.min()), int((d <= window).sum()), bool(disordered[pos]), dis is not None))
    df = pd.DataFrame(rows, columns=['accession', 'gene', 'position', 'kme_sources', 'nearest_methylR', 'methylR_within_window',
                                     'disordered', 'has_disorder_record'])
    df['methylK'] = df['kme_sources'] > 0
    df['histone'] = df['gene'].map(lambda g: bool(HISTONE_GENE.match(g)))
    return df, qc


def mh_test(df: pd.DataFrame, exposure: pd.Series, outcome: pd.Series, strata: list[str]) -> dict:
    """Mantel-Haenszel odds ratio of outcome for exposure, within strata; only informative strata are kept."""
    key = df[strata].astype(str).agg('|'.join, axis=1)
    t = pd.DataFrame({'k': key, 'e': exposure.astype(bool).values, 'o': outcome.astype(bool).values})
    cells = t.groupby(['k', 'e', 'o']).size().unstack(['e', 'o'], fill_value=0)
    for c in [(True, True), (True, False), (False, True), (False, False)]:
        if c not in cells.columns:
            cells[c] = 0
    a, b, c, d = (cells[(True, True)], cells[(True, False)], cells[(False, True)], cells[(False, False)])
    keep = ((a + b) > 0) & ((c + d) > 0) & ((a + c) > 0) & ((b + d) > 0)
    tables = [np.array([[w, x], [y, z]]) for w, x, y, z in zip(a[keep], b[keep], c[keep], d[keep])]
    res = {'strata': ' x '.join(strata), 'informative_strata': int(keep.sum()),
           'lysines': int(len(t)), 'methylK': int(t['o'].sum()), 'exposed': int(t['e'].sum()),
           'methylK_exposed': int((t['o'] & t['e']).sum()),
           'pct_methylK_exposed': 100 * (t['o'] & t['e']).sum() / max(1, t['e'].sum()),
           'pct_methylK_unexposed': 100 * (t['o'] & ~t['e']).sum() / max(1, (~t['e']).sum())}
    if not tables:
        return res | {'mh_or': np.nan, 'ci_low': np.nan, 'ci_high': np.nan, 'cmh_chi2': np.nan, 'log10_p': np.nan}
    st = StratifiedTable(tables)
    lo, hi = st.oddsratio_pooled_confint()
    stat = cmh_statistic([(t[0, 0], t[0, 1], t[1, 0], t[1, 1]) for t in tables])
    # log10 p from the chi-square (1 df) survival function, which stays finite where p underflows to 0.
    return res | {'mh_or': st.oddsratio_pooled, 'ci_low': lo, 'ci_high': hi, 'cmh_chi2': stat, 'log10_p': chi2.logsf(stat, 1) / np.log(10)}


def plot(models: pd.DataFrame, prof: pd.DataFrame, window: int, path_prefix: Path) -> None:
    """Forest plot: (a) odds ratio by distance to the nearest methylarginine, (b) the +/-window models."""
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter

    apply_paper_style()
    a = prof.dropna(subset=['mh_or'])
    a_labels = [f'{r} aa'.replace('-', '\u2013') for r in a['nearest_methylR']]
    b = models[(models['window'] == window) & (models['analysis'] != 'window sensitivity')].drop_duplicates('analysis')
    names = {
        'protein': 'Within protein', 'protein x disorder': '+ disorder', 'protein x disorder x arginines nearby': '+ arginines nearby (main)',
        'without histones': 'Main, without histones', 'evidence: master table': 'Methyl-K: PTM master table',
        'evidence: CPLM 4.0 Exp.': 'Methyl-K: CPLM Exp.', 'evidence: CPLM 4.0 Dat. only': 'Methyl-K: CPLM Dat. only',
        'evidence: independent CPLM Exp. (no union study cited)': 'Methyl-K: CPLM Exp., no shared study',
        'independent CPLM Exp., without histones': 'Same, without histones',
    }
    b_labels = [names.get(n, n) for n in b['analysis']]
    fig, axes = plt.subplots(2, 1, figsize=(6.6, 6.2), sharex=True, gridspec_kw={'height_ratios': [len(a), len(b)], 'hspace': 0.32})
    ink, muted, mark = '#111111', BREWER_COLORS['dark_gray'], BREWER_COLORS['blue']
    for ax, df, labels, title in [(axes[0], a, a_labels, 'a  Lysines by distance to the nearest methylarginine (vs > 50 aa)'),
                                  (axes[1], b, b_labels, f'b  Methylarginine within \u00b1{window} aa: strata and evidence')]:
        y = np.arange(len(df))[::-1]
        ax.axvline(1, color=muted, linewidth=0.8, linestyle=(0, (3, 2)))
        ax.hlines(y, df['ci_low'], df['ci_high'], color=mark, linewidth=2)
        ax.plot(df['mh_or'], y, 'o', color=mark, markersize=6.5, markeredgecolor='white', markeredgewidth=1)
        for yy, (_, r) in zip(y, df.iterrows()):
            ax.text(9.6, yy, f"{r['mh_or']:.2f} ({r['ci_low']:.2f}\u2013{r['ci_high']:.2f})  n={int(r['methylK_exposed'])}",
                    va='center', ha='left', fontsize=8.5, color=muted, clip_on=False)
        ax.set_yticks(y, labels)
        ax.set_title(title, loc='left', fontsize=10.5)
        ax.set_xscale('log')
        ax.set_xlim(0.8, 9)
        ax.xaxis.set_minor_formatter(NullFormatter())
        style_axis(ax, grid_axis='x', force_grid=True)
        ax.spines['left'].set_visible(False)
        ax.tick_params(axis='y', length=0)
    axes[1].set_xticks([1, 2, 3, 5, 8], ['1', '2', '3', '5', '8'])
    axes[1].set_xlabel('Mantel\u2013Haenszel odds ratio of lysine methylation (95% CI)')
    fig.text(0.0, -0.02, 'Strata: protein x MobiDB disorder (a); protein x disorder x arginines within the window (b, from "main"). '
             'n: methylated lysines in the exposed group.', fontsize=7.5, color=muted, ha='left')
    save_figure(fig, path_prefix)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--sites-json', default='docs/data/methylarg_sites.json')
    ap.add_argument('--fasta', default=os.environ.get('PTM_CANONICAL_FASTA', ''))
    ap.add_argument('--disorder-intervals', default=os.environ.get('PTM_MOBIDB_INTERVALS', ''))
    ap.add_argument('--cplm4-sites', default='',
                    help='Validated CPLM 4.0 table from 39_stage_cplm4_lys_methylation.py, for the independent methyl-K test.')
    ap.add_argument('--window', type=int, default=10)
    ap.add_argument('--outdir', default=str(Path(os.environ.get('PTM_RESULTS_ROOT', 'PTM_results')) / 'methylk_methylr'))
    args = ap.parse_args()
    if not args.fasta or not args.disorder_intervals:
        raise SystemExit('--fasta and --disorder-intervals are required (or PTM_CANONICAL_FASTA and PTM_MOBIDB_INTERVALS)')

    export = json.loads(Path(args.sites_json).read_text(encoding='utf-8'))
    sources = export['meta'].get('other_ptm_sources', ['Cross-PTM master table'])
    seqs = parse_fasta(args.fasta)
    disorder = load_disorder(Path(args.disorder_intervals))
    w = args.window
    lys, qc = lysine_table(export, seqs, disorder, w)
    near = lambda df, win: (df['nearest_methylR'] <= win)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Arginines within +/-win of each lysine, binned, for the arginine-density strata.
    def arg_counts(win: int) -> pd.Series:
        out = []
        for acc, grp in lys.groupby('accession', sort=False):
            seq = seqs[acc]
            for pos in grp['position']:
                lo, hi = max(1, pos - win), min(len(seq), pos + win)
                out.append(seq[lo - 1:hi].count('R'))
        return pd.Series(out, index=lys.index)

    rbins = {win: arg_counts(win).map(r_bin) for win in sorted({5, w, 20})}
    dis_ok = lys['has_disorder_record']
    models = []

    def run(name, df, exposure, outcome, strata, extra=None):
        r = mh_test(df, exposure, outcome, strata)
        models.append({'analysis': name, **(extra or {}), **r})

    base = lys.assign(arg_bin=rbins[w])
    full = ['accession', 'disordered', 'arg_bin']
    run('protein', base, near(base, w), base['methylK'], ['accession'], {'window': w})
    run('protein x disorder', base[dis_ok], near(base[dis_ok], w), base[dis_ok]['methylK'], ['accession', 'disordered'], {'window': w})
    run('protein x disorder x arginines nearby', base[dis_ok], near(base[dis_ok], w), base[dis_ok]['methylK'], full, {'window': w})
    for win in sorted({5, w, 20}):
        b = lys.assign(arg_bin=rbins[win])[dis_ok]
        run('window sensitivity', b, near(b, win), b['methylK'], full, {'window': win})
    nh = base[dis_ok & ~base['histone']]
    run('without histones', nh, near(nh, w), nh['methylK'], full, {'window': w})
    # Evidence subsets: each compares one source's methyl-K with unmethylated lysines (other sources' sites left out).
    subsets = [('master table', 1), ('CPLM 4.0 Exp.', 2), ('CPLM 4.0 Dat. only', None)]
    for label, bit in subsets:
        src = base['kme_sources']
        hit = (src & 2 == 0) & (src & 4 > 0) if bit is None else (src & bit > 0)
        sub = base[dis_ok & (hit | ~base['methylK'])]
        run('evidence: ' + label, sub, near(sub, w), hit[sub.index], full, {'window': w})
    if args.cplm4_sites:
        cplm = pd.read_csv(args.cplm4_sites, sep='\t', dtype=str)
        shared = set(UNION_STUDY_PMIDS.values())
        indep = cplm[cplm['evidence'].str.split('/').map(lambda e: 'Exp.' in e)
                     & ~cplm['references'].fillna('').str.split(';').map(lambda refs: bool(shared & set(refs)))]
        keys = set(zip(indep['canonical_UniProtAC'], indep['position'].astype(int)))
        hit = pd.Series([(a, p) in keys for a, p in zip(base['accession'], base['position'])], index=base.index) & base['methylK']
        sub = base[dis_ok & (hit | ~base['methylK'])]
        run('evidence: independent CPLM Exp. (no union study cited)', sub, near(sub, w), hit[sub.index], full, {'window': w})
        sub_nh = sub[~sub['histone']]
        run('independent CPLM Exp., without histones', sub_nh, near(sub_nh, w), hit[sub_nh.index], full, {'window': w})
    # Distance profile against lysines > FAR residues from any methylarginine, within protein x disorder.
    prof = []
    far = base['nearest_methylR'] > FAR
    for a, b in DIST_BINS:
        inb = base['nearest_methylR'].between(a, b)
        sub = base[dis_ok & (inb | far)]
        r = mh_test(sub, inb[sub.index], sub['methylK'], ['accession', 'disordered'])
        prof.append({'nearest_methylR': f'{a}-{b}', **r})
    prof.append({'nearest_methylR': f'>{FAR} (reference)', 'lysines': int((dis_ok & far).sum()), 'methylK': int((dis_ok & far & base['methylK']).sum())})

    models_df = pd.DataFrame(models)
    prof_df = pd.DataFrame(prof)
    save_table(models_df, outdir / 'cmh_models.tsv')
    save_table(prof_df, outdir / 'distance_profile.tsv')
    save_table(lys, outdir / 'lysines.tsv')
    summary = {'window': w, 'other_ptm_sources': sources, 'qc': qc,
               'lysines': int(len(lys)), 'proteins': int(lys['accession'].nunique()),
               'methylK': int(lys['methylK'].sum()), 'proteins_with_methylK': int(lys.loc[lys['methylK'], 'accession'].nunique())}
    (outdir / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    plot(models_df, prof_df, w, outdir / 'methylk_methylr_cooccurrence')
    pd.set_option('display.width', 220)
    cols = ['analysis', 'window', 'informative_strata', 'lysines', 'methylK', 'methylK_exposed', 'pct_methylK_exposed', 'pct_methylK_unexposed', 'mh_or', 'ci_low', 'ci_high', 'cmh_chi2', 'log10_p']
    print(json.dumps(summary, indent=2))
    print(models_df[cols].to_string(index=False, float_format=lambda v: f'{v:.3g}'))
    print(prof_df[['nearest_methylR', 'informative_strata', 'lysines', 'methylK', 'methylK_exposed', 'mh_or', 'ci_low', 'ci_high', 'cmh_chi2', 'log10_p']].to_string(index=False, float_format=lambda v: f'{v:.3g}'))


if __name__ == '__main__':
    main()

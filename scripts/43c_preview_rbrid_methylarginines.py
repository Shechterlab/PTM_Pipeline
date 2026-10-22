#!/usr/bin/env python3
"""Preview: do methylated arginines sit in or near RBR-ID RNA-binding regions more often than the other arginines of
the same proteins? (K562 RBR-ID, Oksuz et al. 2023, staged by scripts/43c_stage_rbrid.py.)

    python -I scripts/43c_preview_rbrid_methylarginines.py [--sites releases/v1.2.0/human_methylarginine_sites_v1.2.0.tsv]
        [--threshold 3.0] [--out preview.tsv]

Two definitions of an RNA-binding region (RBR):
  segment    runs of residues whose residue score is >= --threshold (the staged segments). Residue scores and segments
             are recomputed from the staged peptide table with the staging script's own functions, so other
             thresholds can be tried; at the staged threshold they must equal the staged segment table (checked).
  published  residues under an RBR-ID+ peptide (p < 0.05, depleted; He et al. 2016, Oksuz et al. 2023), read from the
             staged positive-residue table.
Backgrounds (the detectability control), read from the staged coverage table:
  all        every arginine of the analysed proteins;
  measured   arginines under any quantified peptide (coverage 'quantified');
  scored     arginines under a peptide with an RBR-ID score (coverage 'scored'), for segments;
  tested     arginines under a peptide with a t or Fisher test (coverage 'tested'), for the published definition.

Analysis set: proteins with at least one scored residue ("quantified in RBR-ID"); every arginine of their canonical
sequence; exposure = the arginine is a methylarginine of the site table (uniprot_accession, position).
  (a) arginine inside an RBR; (b) within +/-5 residues of one; (c) within +/-10.
      Mantel-Haenszel odds ratio with each protein a stratum, Robins-Breslow-Greenland 95% CI, Cochran-Mantel-Haenszel
      chi-square with continuity correction (R's mantelhaen.test default), and p_rotation: the methylation pattern of
      each protein rotated by a random offset along its arginines (2,000 times), two-sided on log OR.
  (a') (a) without arginines in RG or GR dipeptides (methylated or not).
  (d) residue scores of scored methylated vs other arginines within proteins: mean of the per-protein differences of
      means, and van Elteren's stratified rank-sum statistic (normal p), each with permutation p values (labels
      shuffled within protein; and the methylation pattern rotated along each protein's scored arginines, which keeps
      the spacing of clustered sites). (d') the same without RG/GR arginines.
Residues near each other share peptides and smoothing, and methylarginines cluster, so residues are not independent:
the chi-square, normal and label-shuffle p values are optimistic; the rotation p values are the more cautious ones.
"""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]


def load_stager():
    spec = importlib.util.spec_from_file_location('stage_rbrid', Path(__file__).with_name('43c_stage_rbrid.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class Strata:
    """Rows ordered by protein, then position. Rotating a protein's labels by a random offset along its own rows keeps
    how many of its arginines are methylated and how they are spaced (clusters stay clusters)."""

    def __init__(self, stratum: np.ndarray):
        self.codes, uniq = pd.factorize(stratum)
        self.k = len(uniq)
        self.size = np.bincount(self.codes, minlength=self.k)
        self.first = np.r_[0, np.cumsum(self.size)[:-1]]
        self.rank = np.arange(len(stratum)) - self.first[self.codes]

    def rotate(self, lab: np.ndarray, rng) -> np.ndarray:
        shift = rng.integers(0, self.size)[self.codes]
        return lab[self.first[self.codes] + (self.rank + shift) % self.size[self.codes]]

    def counts(self, e: np.ndarray, o: np.ndarray):
        bc = lambda w: np.bincount(self.codes, weights=w.astype(float), minlength=self.k)  # noqa: E731
        return bc(e & o), bc(e & ~o), bc(~e & o), bc(~e & ~o)


def mantel_haenszel(stratum: np.ndarray, exposed: np.ndarray, outcome: np.ndarray, n_perm: int = 0, seed: int = 1) -> dict:
    """Mantel-Haenszel OR over proteins with both methylated and other arginines; RGB 95% CI; CMH chi-square p with
    continuity correction; and a rotation-permutation p (two-sided on log OR)."""
    st = Strata(stratum)
    a, b, c, d = st.counts(exposed, outcome)
    keep = ((a + b) > 0) & ((c + d) > 0)
    a, b, c, d = a[keep], b[keep], c[keep], d[keep]
    n = a + b + c + d
    R, S = a * d / n, b * c / n
    P, Q = (a + d) / n, (b + c) / n
    res = {'proteins': int(keep.sum()), 'informative_proteins': int((((a + c) > 0) & ((b + d) > 0)).sum()),
           'a': int(a.sum()), 'b': int(b.sum()), 'c': int(c.sum()), 'd': int(d.sum())}
    if R.sum() == 0 or S.sum() == 0:
        return {**res, 'or_mh': np.nan, 'ci_low': np.nan, 'ci_high': np.nan, 'p': np.nan, 'p_rotation': np.nan}
    orr = R.sum() / S.sum()
    var = ((P * R).sum() / (2 * R.sum() ** 2) + (P * S + Q * R).sum() / (2 * R.sum() * S.sum())
           + (Q * S).sum() / (2 * S.sum() ** 2))
    se = np.sqrt(var)
    m = n > 1
    ea = ((a + b) * (a + c) / n)[m]
    va = ((a + b) * (c + d) * (a + c) * (b + d) / (n ** 2 * (n - 1)))[m]
    chi = (abs(a[m].sum() - ea.sum()) - 0.5) ** 2 / va.sum()
    p_rot = np.nan
    if n_perm:
        rng = np.random.default_rng(seed)
        obs, ge = abs(np.log(orr)), 0
        for _ in range(n_perm):
            pa, pb, pc, pd_ = (x[keep] for x in st.counts(st.rotate(exposed, rng), outcome))
            r_, s_ = (pa * pd_ / n).sum(), (pb * pc / n).sum()
            ge += (not (r_ > 0 and s_ > 0)) or abs(np.log(r_ / s_)) >= obs - 1e-12
        p_rot = (ge + 1) / (n_perm + 1)
    return {**res, 'or_mh': orr, 'ci_low': np.exp(np.log(orr) - 1.96 * se), 'ci_high': np.exp(np.log(orr) + 1.96 * se),
            'p': float(stats.chi2.sf(chi, 1)), 'p_rotation': p_rot}


def within_protein_scores(stratum, methyl, score, n_perm: int, seed: int) -> dict:
    """Mean over proteins of (mean score of methylated - mean of other covered arginines), with permutation p values,
    and van Elteren's stratified Wilcoxon test."""
    df = pd.DataFrame({'s': stratum, 'm': methyl, 'x': score})
    g = df.groupby('s')['m'].agg(['sum', 'size'])
    keep = g[(g['sum'] > 0) & (g['sum'] < g['size'])].index
    df = df[df.s.isin(keep)].reset_index(drop=True)  # already ordered by protein and position
    st = Strata(df.s.to_numpy())
    codes, k, size = st.codes, st.k, st.size
    x, m = df.x.to_numpy(), df.m.to_numpy()
    r = df.groupby(codes)['x'].rank().to_numpy()  # within-protein ranks, ties averaged
    nm = np.bincount(codes, weights=m.astype(float), minlength=k)
    tot = np.bincount(codes, weights=x, minlength=k)
    w = 1.0 / (size + 1)
    ve_mean = nm * (size + 1) / 2.0

    def per_protein(lab):
        sm = np.bincount(codes, weights=x * lab, minlength=k)
        return sm / nm - (tot - sm) / (size - nm)

    def ve_num(lab):  # van Elteren: sum over proteins of (rank sum of methylated - its expectation) / (n + 1)
        return float((w * (np.bincount(codes, weights=r * lab, minlength=k) - ve_mean)).sum())

    lab0 = m.astype(float)
    per = per_protein(lab0)
    obs, obs_ve = float(per.mean()), ve_num(lab0)
    rng = np.random.default_rng(seed)
    ge = {'shuffle': 0, 'rotation': 0, 've_shuffle': 0, 've_rotation': 0}
    for _ in range(n_perm):
        order = np.lexsort((rng.random(len(x)), codes))
        lab = np.empty(len(x))
        lab[order] = lab0  # each protein's own labels, randomly rearranged within it
        lab_r = st.rotate(lab0, rng)
        for name, L in [('shuffle', lab), ('rotation', lab_r)]:
            ge[name] += abs(float(per_protein(L).mean())) >= abs(obs) - 1e-12
            ge['ve_' + name] += abs(ve_num(L)) >= abs(obs_ve) - 1e-12
    den = 0.0
    for _, grp in df.groupby(codes):
        n_, mm = len(grp), int(grp.m.sum())
        _, t = np.unique(grp.x.to_numpy(), return_counts=True)
        den += (1.0 / (n_ + 1)) ** 2 * mm * (n_ - mm) / 12.0 * ((n_ + 1) - (t ** 3 - t).sum() / (n_ * (n_ - 1)))
    z = obs_ve / np.sqrt(den) if den > 0 else np.nan
    P = lambda c: (c + 1) / (n_perm + 1)  # noqa: E731
    return {'proteins': k, 'methylated_arginines': int(m.sum()), 'other_arginines': int(len(m) - m.sum()),
            'mean_within_protein_difference': obs, 'median_within_protein_difference': float(np.median(per)),
            'proteins_methylated_higher': int((per > 0).sum()), 'proteins_methylated_lower': int((per < 0).sum()),
            'p_mean_label_shuffle': P(ge['shuffle']), 'p_mean_rotation': P(ge['rotation']),
            'van_elteren_z': z, 'van_elteren_p': float(2 * stats.norm.sf(abs(z))),
            'p_van_elteren_label_shuffle': P(ge['ve_shuffle']), 'p_van_elteren_rotation': P(ge['ve_rotation'])}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--peptides', type=Path, default=ROOT / 'inputs' / 'interactions' / 'rbrid_oksuz2023_k562_peptides.tsv')
    ap.add_argument('--segments', type=Path, default=ROOT / 'inputs' / 'interactions' / 'rbrid_oksuz2023_k562_segments.tsv')
    ap.add_argument('--positive-residues', type=Path,
                    default=ROOT / 'inputs' / 'interactions' / 'rbrid_oksuz2023_k562_positive_residues.tsv')
    ap.add_argument('--coverage', type=Path, default=ROOT / 'inputs' / 'interactions' / 'rbrid_oksuz2023_k562_coverage.tsv')
    ap.add_argument('--fasta', type=Path, default=ROOT / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03.fasta')
    ap.add_argument('--sites', type=Path, default=ROOT / 'releases' / 'v1.2.0' / 'human_methylarginine_sites_v1.2.0.tsv')
    ap.add_argument('--threshold', type=float, default=3.0, help='the staged threshold; others for sensitivity')
    ap.add_argument('--also-thresholds', default='1,2,5', help='(a) repeated at these thresholds')
    ap.add_argument('--score', default='rbrid_score', help='peptide score column (rbrid_score_5v5 for the 5 x 5 test)')
    ap.add_argument('--permutations', type=int, default=10000, help='for (d)')
    ap.add_argument('--mh-permutations', type=int, default=2000, help='rotations for the p_rotation of (a)-(c)')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out', type=Path)
    args = ap.parse_args()
    st = load_stager()

    seqs = {}
    acc = None
    for line in open(args.fasta, encoding='utf-8'):
        if line.startswith('>'):
            acc = line.split('|')[1]
            seqs[acc] = []
        else:
            seqs[acc].append(line.strip())
    seqs = {a: ''.join(v) for a, v in seqs.items()}
    lengths = {a: len(s) for a, s in seqs.items()}
    pep = pd.read_csv(args.peptides, sep='\t', keep_default_na=False, na_values=[''])
    prof = st.residue_profiles(pep['uniprot_accession'].to_numpy(), pep['start'].to_numpy(), pep['end'].to_numpy(),
                               pep[args.score].to_numpy(dtype=float), lengths)

    sites = pd.read_csv(args.sites, sep='\t', usecols=['uniprot_accession', 'position', 'residue'], dtype={'position': int})
    ok = sites.apply(lambda r: r.uniprot_accession in seqs and r.position <= len(seqs[r.uniprot_accession])
                     and seqs[r.uniprot_accession][r.position - 1] == 'R', axis=1)
    print(f'sites: {len(sites):,} in {sites.uniprot_accession.nunique():,} proteins; {int((~ok).sum())} not an R of the '
          f'reference sequence (left out)')
    methyl = set(zip(sites.uniprot_accession[ok], sites.position[ok]))

    # measured background (staged coverage intervals) and the published-cut-off residues
    covsets = {}
    for kind, g in pd.read_csv(args.coverage, sep='\t').groupby('coverage'):
        d_ = {}
        for a, s0, e0 in zip(g.uniprot_accession, g.start, g.end):
            d_.setdefault(a, []).append((s0, e0))
        covsets[kind] = {a: np.array(v) for a, v in d_.items()}
    posres = pd.read_csv(args.positive_residues, sep='\t', usecols=['uniprot_accession', 'position'])
    pos_by = {a: np.sort(g.position.to_numpy()) for a, g in posres.groupby('uniprot_accession')}

    def inside(iv, p):
        if iv is None:
            return np.zeros(len(p), bool)
        return ((p[:, None] >= iv[None, :, 0]) & (p[:, None] <= iv[None, :, 1])).any(1)

    rows = []
    for a in sorted(prof):
        s = seqs[a]
        sc = dict(zip(prof[a][0].tolist(), prof[a][2].tolist()))
        rpos = np.array([i + 1 for i, ch in enumerate(s) if ch == 'R'], dtype=int)
        if not len(rpos):
            continue
        fl = {k: inside(covsets.get(k, {}).get(a), rpos) for k in ('quantified', 'tested', 'scored')}
        for j, p in enumerate(rpos.tolist()):
            rows.append((a, p, (a, p) in methyl, fl['quantified'][j], fl['tested'][j], fl['scored'][j], sc.get(p, np.nan),
                         (p < len(s) and s[p] == 'G') or (p > 1 and s[p - 2] == 'G')))
    R = pd.DataFrame(rows, columns=['acc', 'pos', 'methyl', 'measured', 'tested', 'scored', 'score', 'rg_gr'])
    if not (R.scored == R.score.notna()).all():
        raise SystemExit('the scored coverage intervals and the residue scores disagree')
    n_meth_prot = R[R.methyl].acc.nunique()
    print(f'proteins quantified (scored residues): {len(prof):,}; with a methylarginine: {n_meth_prot:,} of '
          f'{len({a for a, _ in methyl}):,} methylated proteins; methylarginines in them: {int(R.methyl.sum()):,} of '
          f'{len(methyl):,}; other arginines: {int((~R.methyl).sum()):,}')
    for k in ('measured', 'tested', 'scored'):
        print(f'  {k}: methylarginines {int((R.methyl & R[k]).sum()):,} ({(R[k][R.methyl]).mean():.1%}), other arginines '
              f'{int((~R.methyl & R[k]).sum()):,} ({(R[k][~R.methyl]).mean():.1%})')

    def segment_distances(threshold):
        segs = st.find_segments(prof, threshold)
        by = {}
        for a, s0, e0, _ in segs:
            by.setdefault(a, []).append((s0, e0))
        dist = np.full(len(R), np.inf)
        for a, grp in R.groupby('acc'):
            if a in by:
                iv = np.array(by[a])
                p = grp.pos.to_numpy()[:, None]
                dist[grp.index.to_numpy()] = np.maximum(0, np.maximum(iv[:, 0][None, :] - p, p - iv[:, 1][None, :])).min(1)
        return segs, dist

    def point_distances(points: dict):
        dist = np.full(len(R), np.inf)
        for a, grp in R.groupby('acc'):
            q = points.get(a)
            if q is not None and len(q):
                p = grp.pos.to_numpy()
                i = np.clip(np.searchsorted(q, p), 1, len(q)) - 1
                j = np.clip(i + 1, 0, len(q) - 1)
                dist[grp.index.to_numpy()] = np.minimum(np.abs(q[i] - p), np.abs(q[j] - p))
        return dist

    segs, dist_seg = segment_distances(args.threshold)
    if args.segments.exists():
        staged = pd.read_csv(args.segments, sep='\t')
        same = sorted((a, s0, e0) for a, s0, e0, _ in segs) == sorted(zip(staged.uniprot_accession, staged.start, staged.end))
        print(f'segments at {args.threshold}: {len(segs):,}; equal to the staged table: {same}')
    dist_pub = point_distances(pos_by)
    print(f'published-cut-off residues in these proteins: {len(posres[posres.uniprot_accession.isin(prof)]):,}')

    stratum, meth = R.acc.to_numpy(), R.methyl.to_numpy()
    nrg = ~R.rg_gr.to_numpy()

    def mask(bg):
        return np.ones(len(R), bool) if bg == 'all' else R[bg].to_numpy()

    def mh(m, outcome):
        return mantel_haenszel(stratum[m], meth[m], outcome[m], args.mh_permutations, args.seed)

    out = []
    for rbr, dist, own, thr in [('segment', dist_seg, 'scored', f'score >= {args.threshold:g}'),
                                ('published', dist_pub, 'tested', 'RBR-ID+ peptide (p < 0.05, depleted)')]:
        for label, k in [('(a) inside', 0), ('(b) within +/-5', 5), ('(c) within +/-10', 10)]:
            for bg in ('all', 'measured', own):
                out.append({'rbr': rbr, 'rule': thr, 'test': label, 'arginines': bg, **mh(mask(bg), dist <= k)})
        for bg in ('all', 'measured', own):
            out.append({'rbr': rbr, 'rule': thr, 'test': "(a') inside, RG/GR excluded", 'arginines': bg,
                        **mh(mask(bg) & nrg, dist <= 0)})
    for T in [float(x) for x in args.also_thresholds.split(',') if x]:
        _, dT = segment_distances(T)
        for bg in ('all', 'measured', 'scored'):
            out.append({'rbr': 'segment', 'rule': f'score >= {T:g}', 'test': '(a) inside', 'arginines': bg,
                        **mh(mask(bg), dT <= 0)})
    res = pd.DataFrame(out)
    pd.set_option('display.width', 250)
    pd.set_option('display.max_colwidth', 40)
    print(res.to_string(index=False, float_format=lambda v: f'{v:.4g}'))
    sc_ = R[R.scored]
    d = within_protein_scores(sc_.acc.to_numpy(), sc_.methyl.to_numpy(), sc_.score.to_numpy(), args.permutations, args.seed)
    print('(d) residue scores, scored arginines, within protein:')
    for k_, v in d.items():
        print(f'    {k_}: {v:.4g}' if isinstance(v, float) else f'    {k_}: {v}')
    sc2 = sc_[~sc_.rg_gr]
    d2 = within_protein_scores(sc2.acc.to_numpy(), sc2.methyl.to_numpy(), sc2.score.to_numpy(), args.permutations, args.seed)
    print("(d') the same without RG/GR arginines:")
    for k_, v in d2.items():
        print(f'    {k_}: {v:.4g}' if isinstance(v, float) else f'    {k_}: {v}')
    if args.out:
        res.to_csv(args.out, sep='\t', index=False)
        pd.DataFrame([{'test': '(d) residue score', **d}, {'test': "(d') residue score, no RG/GR", **d2}]).to_csv(
            args.out.with_name(args.out.stem + '_scores.tsv'), sep='\t', index=False)


if __name__ == '__main__':
    main()

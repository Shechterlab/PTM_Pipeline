"""The per-site mutation propensity of the change tables (release v1.7.0): codon mutability, the within-protein binomial test,
the site-set test, the positive controls and the permutation check. Read by scripts/71_release_propensity.py and the tests.

Mutability. For each coding base of a protein's coding sequence (inputs/layers/propensity_cds.tsv.gz) and each substitution
that changes its residue (missense or stop), the rate of that substitution in its trinucleotide (transcript strand; at an
exon edge the intron neighbor): for TCGA and DepMap the source's own spectrum (propensity_somatic_spectra.tsv), for gnomAD
and ClinVar gnomAD's germline rates (mutation_rate_methylation_bins), with the CpG methylation level for a CpG transition.
A residue's mutability w is the sum over its codon's three bases. An entry with several coding sequences (several genes,
as the histones) takes their mean.

Test. For a source, a residue type X and a protein with n changes at its X residues (k_i at residue i; a change is a tumor
sample, cell line, ClinVar variation or gnomAD variant), expected_i = n w_i / sum w over the protein's X residues, and
k_i ~ Binomial(n, w_i / sum w). p_upper is the chance of k_i or more, p_lower of k_i or fewer. Benjamini-Hochberg q within
the source, cohort and residue type, over every residue tested (every X residue of a protein with n >= 1 and weights).
"""
from __future__ import annotations

import gzip
import io
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import binom

B = 'ACGT'
IDX = {b: i for i, b in enumerate(B)}
COMP = str.maketrans('ACGTN', 'TGCAN')
CODON_ORDER = 'TCAG'
CODE = dict(zip((a + b + c for a in CODON_ORDER for b in CODON_ORDER for c in CODON_ORDER),
                'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
SEED = 20261120
ITERATIONS = 10_000
ALPHA = 0.05


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def collapse(context: str, alt: str) -> tuple[str, str, str]:
    if context[1] in 'AG':
        return rc(context), context[1].translate(COMP), alt.translate(COMP)
    return context, context[1], alt


CLASSES = [(l + m + r, m, a) for m in 'CT' for l in B for r in B for a in B if a != m]
CLASS_INDEX = {c: k for k, c in enumerate(CLASSES)}
CLASS_OF = np.full(65 * 4, -1, dtype=np.int64)          # (trinucleotide code * 4 + alt) -> class; code 64 = unknown
for _t in range(64):
    _ctx = B[_t // 16] + B[(_t // 4) % 4] + B[_t % 4]
    for _a in range(4):
        if B[_a] != _ctx[1]:
            CLASS_OF[_t * 4 + _a] = CLASS_INDEX[collapse(_ctx, B[_a])]
CPG_TRANSITION = np.array([c[0][1:] == 'CG' and c[2] == 'T' for c in CLASSES])
NS = np.zeros((65, 3, 4), dtype=bool)                   # codon code 64 = a codon with an unknown base
for _c in range(64):
    _codon = B[_c // 16] + B[(_c // 4) % 4] + B[_c % 4]
    for _o in range(3):
        for _a in range(4):
            if B[_a] != _codon[_o]:
                NS[_c, _o, _a] = CODE[_codon] != CODE[_codon[:_o] + B[_a] + _codon[_o + 1:]]
LUT = np.full(256, 4, dtype=np.int64)
for _b, _i in IDX.items():
    LUT[ord(_b)] = _i


def somatic_rates(spectra: pd.DataFrame) -> dict[str, np.ndarray]:
    """'<source>/<cohort>' -> rate per class (96)."""
    out = {}
    for (src, cohort), g in spectra.groupby(['source', 'cohort'], sort=False):
        r = np.zeros(len(CLASSES))
        for ctx, ref, alt, rate in zip(g['context'], g['ref'], g['alt'], g['rate']):
            r[CLASS_INDEX[(ctx, ref, alt)]] = float(rate) if rate != '' else 0.0
        out[f'{src}/{cohort}'] = r
    return out


def germline_rates(path) -> np.ndarray:
    """gnomAD mutation_rate_methylation_bins -> rate per class (96) and methylation level (0, 1, 2); a class without levels
    1 and 2 (any but a CpG transition) takes its level-0 rate there."""
    with gzip.open(path, 'rt') as fh:
        t = pd.read_csv(fh, sep='\t', dtype={'context': str, 'ref': str, 'alt': str, 'methylation_level': int, 'mu_snp': float})
    r = np.full((len(CLASSES), 3), np.nan)
    for ctx, ref, alt, lv, mu in zip(t['context'], t['ref'], t['alt'], t['methylation_level'], t['mu_snp']):
        if ctx[1] != ref:
            raise ValueError(f'context {ctx} does not center on {ref}')
        r[CLASS_INDEX[collapse(ctx, alt)], lv] = mu
    if np.isnan(r[:, 0]).any():
        raise ValueError('a class has no level-0 rate')
    for lv in (1, 2):
        missing = np.isnan(r[:, lv])
        if (missing & CPG_TRANSITION).any():
            raise ValueError(f'a CpG transition has no level-{lv} rate')
        r[missing, lv] = r[missing, 0]
    return r


def contexts(row) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Base codes, left and right neighbor codes (transcript strand, genomic) and CpG methylation level of each coding base."""
    codes = LUT[np.frombuffer(row['cds'].encode(), dtype=np.uint8)]
    n = len(codes)
    left, right = np.empty(n, dtype=np.int64), np.empty(n, dtype=np.int64)
    left[1:], right[:-1] = codes[:-1], codes[1:]
    left[0], right[-1] = IDX.get(row['flank_5'], 4), IDX.get(row['flank_3'], 4)
    for item in filter(None, row['junctions'].split(',')):
        j, pair = item.split(':')
        j = int(j)
        right[j - 1], left[j] = IDX.get(pair[0], 4), IDX.get(pair[1], 4)
    level = np.zeros(n, dtype=np.int64)
    for item in filter(None, row['cpg_methylation'].split(',')):
        j, lv = item.split(':')
        level[int(j)] = int(lv)
    return codes, left, right, level


def codon_mutability(row, somatic: dict[str, np.ndarray], germline: np.ndarray) -> dict[str, np.ndarray]:
    """Per residue of one coding sequence: w under each somatic spectrum and under the germline rates."""
    codes, left, right, level = contexts(row)
    n = len(codes)
    known = (codes < 4) & (left < 4) & (right < 4)
    tri = np.where(known, left * 16 + codes * 4 + right, 64)
    cls = CLASS_OF[tri[:, None] * 4 + np.arange(4)[None, :]]                  # (n, 4), -1 where alt = base or unknown
    cod = codes.reshape(-1, 3)
    cc = np.where((cod < 4).all(axis=1), cod[:, 0] * 16 + cod[:, 1] * 4 + cod[:, 2], 64)
    ns = NS[np.repeat(cc, 3), np.tile(np.arange(3), len(cc))] & (cls >= 0)
    safe = np.where(cls >= 0, cls, 0)
    out = {}
    for key, rate in somatic.items():
        out[key] = (np.where(ns, rate[safe], 0.0)).reshape(-1, 3, 4).sum(axis=(1, 2))
    lv = np.where(CPG_TRANSITION[safe], level[:, None], 0)
    out['germline'] = (np.where(ns, germline[safe, lv], 0.0)).reshape(-1, 3, 4).sum(axis=(1, 2))
    return out


def mutability(cds: pd.DataFrame, somatic: dict[str, np.ndarray], germline: np.ndarray) -> dict[str, dict[str, np.ndarray]]:
    """accession -> {spectrum key -> w per residue (the mean over the entry's coding sequences)}."""
    parts = defaultdict(list)
    for _, row in cds.iterrows():
        parts[row['canonical_UniProtAC']].append(codon_mutability(row, somatic, germline))
    out = {}
    for acc, ws in parts.items():
        out[acc] = {k: np.mean([w[k] for w in ws], axis=0) for k in ws[0]}
    return out


def bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values (NaN in, NaN out; NaNs not counted)."""
    q = np.full(len(p), np.nan)
    ok = ~np.isnan(p)
    m = int(ok.sum())
    if not m:
        return q
    pv = p[ok]
    order = np.argsort(pv, kind='stable')
    ranked = pv[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def site_test(res: pd.DataFrame, w_col: str, k_col: str, lower: bool = False) -> pd.DataFrame:
    """res: one row per candidate residue (canonical_UniProtAC, position, w_col, k_col) of one residue type. Adds n, residues,
    relative mutability, expected, O/E (weighted and unweighted), p_upper (and p_lower) and q within the rows tested."""
    g = res.groupby('canonical_UniProtAC', sort=False)
    n = g[k_col].transform('sum').to_numpy(dtype=float)
    count = g[k_col].transform('size').to_numpy(dtype=float)
    w = res[w_col].to_numpy(dtype=float)
    wsum = g[w_col].transform('sum').to_numpy(dtype=float)
    has_w = g[w_col].transform(lambda s: s.notna().all() and (s > 0).any()).to_numpy(dtype=bool)
    k = res[k_col].to_numpy(dtype=float)
    pi = np.where(has_w, w / np.where(wsum > 0, wsum, 1), np.nan)
    tested = has_w & (n >= 1)
    expected = np.where(tested, n * pi, np.nan)
    p_up = np.where(tested, binom.sf(k - 1, n, np.nan_to_num(pi)), np.nan)
    out = pd.DataFrame({'n': n.astype(int), 'residues_of_type': count.astype(int),
                        'relative_mutability': np.where(has_w, w / np.where(wsum > 0, wsum, 1) * count, np.nan),
                        'expected': expected, 'oe': np.where(tested & (expected > 0), k / np.where(expected > 0, expected, 1), np.nan),
                        'expected_unweighted': np.where(n >= 1, n / count, np.nan),
                        'oe_unweighted': np.where(n >= 1, k / np.where(n >= 1, n / count, 1), np.nan),
                        'p_upper': p_up, 'q_upper': bh(p_up)}, index=res.index)
    if lower:
        p_lo = np.where(tested, binom.cdf(k, n, np.nan_to_num(pi)), np.nan)
        out['p_lower'], out['q_lower'] = p_lo, bh(p_lo)
    return out


def changes_needed(k: int, n: int, pi: float, m: int, alpha: float = ALPHA) -> int | None:
    """The fewest changes at the residue (the protein's other changes unchanged) whose p_upper is at most alpha / m: enough
    for q < alpha whatever the other residues' p (Bonferroni bound)."""
    if not (pi > 0) or m < 1:
        return None
    other = n - k
    for kk in range(0, 100_000):
        # the residue's share stays pi of the protein's weight; n grows with kk
        if binom.sf(kk - 1, other + kk, pi) <= alpha / m:
            return kk
    return None


def set_test(units: list[tuple[int, float, int]], rng: np.random.Generator, iterations: int = ITERATIONS) -> dict:
    """units: per protein (n, share of the protein's weight on the set's residues, changes observed at them). The total at
    the set's residues under the null is a sum of Binomial(n, share); p from `iterations` draws."""
    n = np.array([u[0] for u in units], dtype=np.int64)
    share = np.array([u[1] for u in units], dtype=float)
    obs = int(sum(u[2] for u in units))
    expected = float((n * share).sum())
    sims = np.zeros(iterations, dtype=np.int64)
    for start in range(0, len(n), 500):
        sl = slice(start, start + 500)
        sims += rng.binomial(n[sl][:, None], share[sl][:, None], size=(len(n[sl]), iterations)).sum(axis=0)
    return {'observed': obs, 'expected': expected, 'oe': obs / expected if expected > 0 else np.nan,
            'p_upper': (1 + int((sims >= obs).sum())) / (iterations + 1), 'p_lower': (1 + int((sims <= obs).sum())) / (iterations + 1)}


def permutation_check(units: list[tuple[np.ndarray, np.ndarray]], rng: np.random.Generator, iterations: int = ITERATIONS) -> dict:
    """units: per protein and residue type, (k per residue, w per residue). Each protein's changes redrawn over its residues
    in proportion to w (multinomial), `iterations` times; the per-residue permutation p against the binomial p."""
    pb, pp = [], []
    for k, w in units:
        n = int(k.sum())
        pi = w / w.sum()
        sims = rng.multinomial(n, pi, size=iterations)
        pp.append((1 + (sims >= k[None, :]).sum(axis=0)) / (iterations + 1))
        pb.append(binom.sf(k - 1, n, pi))
    pb, pp = np.concatenate(pb), np.concatenate(pp)
    qb, qp = bh(pb), bh(pp)
    cb, cp = qb < ALPHA, qp < ALPHA
    nb, np_ = pb < ALPHA, pp < ALPHA
    diff = np.abs(pb - pp)
    return {'units': len(units), 'residues': int(len(pb)), 'binomial_q05': int(cb.sum()), 'permutation_q05': int(cp.sum()),
            'both_q05': int((cb & cp).sum()), 'agreement': float((cb == cp).mean()) if len(pb) else np.nan,
            'binomial_p05': int(nb.sum()), 'permutation_p05': int(np_.sum()), 'both_p05': int((nb & np_).sum()),
            'agreement_p05': float((nb == np_).mean()) if len(pb) else np.nan,
            'median_abs_p_difference': float(np.median(diff)) if len(pb) else np.nan,
            'max_abs_p_difference': float(diff.max()) if len(pb) else np.nan,
            'max_abs_p_difference_p_above_0_001': float(diff[pb > 1e-3].max()) if (pb > 1e-3).any() else np.nan}


def fmt(x, digits: int = 6) -> str:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return ''
    return f'{x:.{digits}g}'


def write_gz(df: pd.DataFrame, path) -> None:
    buf = io.StringIO()
    df.to_csv(buf, sep='\t', index=False, lineterminator='\n')
    with open(path, 'wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0, compresslevel=9) as gz:
        gz.write(buf.getvalue().encode('utf-8'))

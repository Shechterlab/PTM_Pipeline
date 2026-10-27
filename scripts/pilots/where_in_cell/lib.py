"""Shared pieces of the where-in-the-cell pilots: reviewed entries, methylated proteins, PaxDb, and the per-compartment
tests (one-sided hypergeometric with fold enrichment among the dataset's localised proteins, and the Mantel-Haenszel odds
ratio within PaxDb deciles of the same universe, two-sided; BH within each test family)."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import hypergeom

def setup(repo):
    sys.path.insert(0, str(Path(repo) / 'scripts'))
    import abundance_strata as A
    from common import mantel_haenszel_ci, benjamini_hochberg
    entries = pd.read_csv(Path(repo) / 'inputs/reference/uniprot_human_reviewed_2026_03_entries.tsv', sep='\t', dtype=str, keep_default_na=False)
    by_symbol = {}
    for e, g in zip(entries['Entry'], entries['Gene Names (primary)']):
        for s in g.split(';'):
            if s.strip(): by_symbol.setdefault(s.strip(), set()).add(e)
    prot = pd.read_csv(Path(repo) / 'releases/v1.4.1/human_methylarginine_proteins_v1.4.1.tsv', sep='\t', dtype=str, keep_default_na=False)
    return dict(A=A, mh=mantel_haenszel_ci, bh=benjamini_hochberg, reviewed=set(entries['Entry']), by_symbol=by_symbol,
                meth={a.split('-')[0] for a in prot['uniprot_accession']}, pax=A.abundance(Path(repo) / 'inputs')['protein'],
                gene={a.split('-')[0]: g for a, g in zip(prot['uniprot_accession'], prot['gene'])})

def resolve(ctx, ids, symbol=''):
    """The reviewed entry an identifier list names (first reviewed accession), else the one entry its gene symbol names."""
    for a in ids:
        a = a.strip().split('-')[0]
        if a in ctx['reviewed']: return a
    hits = ctx['by_symbol'].get(symbol, set())
    return next(iter(hits)) if len(hits) == 1 else ''

def tests(ctx, universe, sets, min_size=5):
    meth, pax = ctx['meth'], ctx['pax']
    universe = sorted(universe)
    is_t = np.array([u in meth for u in universe]); n, N = int(is_t.sum()), len(universe)
    codes, layers = np.unique(ctx['A'].strata(universe, pax).to_numpy(dtype=str), return_inverse=True)
    rows = []
    for name, members in sets.items():
        ins = np.array([u in members for u in universe]); M = int(ins.sum()); k = int((ins & is_t).sum())
        if M < min_size: continue
        tables = [(int((s & ins & is_t).sum()), int((s & is_t & ~ins).sum()), int((s & ins & ~is_t).sum()), int((s & ~ins & ~is_t).sum()))
                  for s in (layers == i for i in range(len(codes)))]
        odds, lo, hi, pmh, used = ctx['mh'](tables)
        raw = (k * (N - M - n + k)) / ((n - k) * (M - k)) if (n - k) and (M - k) else float('inf')
        rows.append(dict(compartment=name, members=M, methylated=k, share=k / M, fold=(k / n) / (M / N),
                         p=float(hypergeom.sf(k - 1, N, M, n)) if k else 1.0, odds_raw=raw, odds_mh=odds, lo_mh=lo, hi_mh=hi, p_mh=pmh))
    t = pd.DataFrame(rows)
    t['q'] = ctx['bh'](t['p']); t['q_mh'] = ctx['bh'](t['p_mh'].fillna(1))
    return t.sort_values('fold', ascending=False)

def show(t, title):
    print('\n==', title)
    print(t[['compartment', 'members', 'methylated', 'share', 'fold', 'q', 'odds_raw', 'odds_mh', 'q_mh']].to_string(index=False, float_format=lambda v: f'{v:.3g}'))

"""The one rule that turns peptide-level changes into a call for each methylarginine the peptides carry.

Used by scripts/build_release.py (Maron et al. 2021 Table S3 inhibitor responses, the CARM1 knockout of Shishkova et
al. 2017) and scripts/16_maron_state_treatment_context.py, so the release, the website and the analyses agree.

A peptide changes when its adjusted p value is at most the cut-off (and, where given, its |log2| reaches a minimum).
A site's call over its peptides: decreased when at least one peptide decreased and none increased, increased
likewise, mixed when peptides changed in opposite directions, no_significant_change when none changed. The median
log2 over all peptides is descriptive only: it is not the effect of any test, and it can have the opposite sign to the
call (one peptide significantly up, three not significantly down).

A peptide reports the methyl groups it carries together. When it carries several methylarginines, its change cannot be
given to one of them, so each site also gets the call from its single-site peptides alone (peptides carrying no other
methylarginine), empty when it has none: that call, unlike the other, can be read as a change at this site.
"""
from __future__ import annotations

import math

NO_CHANGE = 'no_significant_change'


def bh(p: list) -> list:
    """Benjamini-Hochberg q values (None stays None)."""
    idx = [k for k, v in enumerate(p) if v is not None and not (isinstance(v, float) and math.isnan(v))]
    q = [None] * len(p)
    m = len(idx)
    prev = 1.0
    for rank, k in reversed(list(enumerate(sorted(idx, key=lambda k: p[k]), 1))):
        prev = min(prev, p[k] * m / rank)
        q[k] = prev
    return q


def peptide_call(n_down: int, n_up: int) -> str:
    return 'mixed' if n_down and n_up else 'decreased' if n_down else 'increased' if n_up else NO_CHANGE


def combine_calls(calls) -> str:
    """Call over rows (antibodies, cell lines): decreased when a row is decreased and none increased or mixed; increased
    likewise; mixed when rows disagree or one is mixed; no_significant_change when all are; empty with no row."""
    calls = {c for c in calls if c}
    if not calls:
        return ''
    if calls <= {NO_CHANGE}:
        return NO_CHANGE
    if 'mixed' in calls or {'decreased', 'increased'} <= calls:
        return 'mixed'
    return 'decreased' if 'decreased' in calls else 'increased'


def site_response(peptides: list, cutoff: float, min_abs_log2: float = 0.0, strict: bool = False) -> dict:
    """peptides: [(adjusted p, log2 change, methylarginines on the peptide), ...] for one site and condition. A peptide
    changes at p <= cutoff (p < cutoff when strict) with |log2| >= min_abs_log2."""
    missing = lambda v: v is None or (isinstance(v, float) and math.isnan(v))  # noqa: E731
    peptides = [(None if missing(q) else q, d, n) for q, d, n in peptides]
    changed = lambda q, d: q is not None and (q < cutoff if strict else q <= cutoff) and abs(d) >= min_abs_log2  # noqa: E731
    down = sum(changed(q, d) and d < 0 for q, d, _ in peptides)
    up = sum(changed(q, d) and d > 0 for q, d, _ in peptides)
    single = [(q, d) for q, d, n in peptides if n == 1]
    s_down = sum(changed(q, d) and d < 0 for q, d in single)
    s_up = sum(changed(q, d) and d > 0 for q, d in single)
    ds = sorted(d for _, d, _ in peptides)
    mid = len(ds) // 2
    median = (ds[mid] if len(ds) % 2 else (ds[mid - 1] + ds[mid]) / 2) if ds else None
    best = min(peptides, key=lambda x: (x[0] is None, x[0] if x[0] is not None else 0)) if peptides else (None, None, 0)
    return {'n_peptides': len(peptides), 'n_decreased': down, 'n_increased': up,
            'n_multi_site_peptides': sum(n > 1 for _, _, n in peptides), 'log2_median': median,
            'log2_strongest': best[1], 'adj_p_min': best[0], 'call': peptide_call(down, up) if peptides else '',
            'call_single_site_peptides': peptide_call(s_down, s_up) if single else ''}

"""Place lysine sites of a mass-spectrometry site table (MaxQuant style) on the canonical reviewed sequences of the reference.

Used by scripts/60 (Hendriks 2017, SUMO-2) and scripts/61 (Akimov 2018, ubiquitin). A source row names a leading protein
and its position, the other proteins of its protein group with theirs, and a sequence window centred on the lysine. The
row is placed by the +/-7 residue core of its window ('_' pads past either end of the sequence), as the CPLM staging
(scripts/39) checks its sites:
  kept                 the leading protein, or a later protein of the group, is a canonical accession of the reference,
                       which has K at the given position with the window core agreeing
  secondary_accession  as kept, after a secondary (merged) accession is replaced by its current primary accession
  relocated            no given position fits, but the core occurs exactly once in the canonical sequence of a protein of
                       the group (an isoform's canonical entry, or a sequence that changed since the study): the site
                       moves to that K
The group is tried leading protein first (then the other proteins of the leading protein's entry, then the rest, in the
order given): first every given position, then every core search, and the first fit wins. Otherwise the row is rejected:
  no_canonical_accession   no protein of the group is (or maps by secondary accession to) a canonical reference entry
  window_not_unique        the core occurs more than once in a candidate sequence and nowhere exactly once
  window_not_in_canonical  the core occurs in no candidate sequence
  bad_window               the window has no K at its centre
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

FLANK = 7
REJECTIONS = ('no_canonical_accession', 'window_not_unique', 'window_not_in_canonical', 'bad_window')


def read_secondary(history: str | Path, reference: dict[str, str]) -> dict[str, str]:
    """Secondary accession -> current primary accession, from the reference's accession history (kind 'secondary'),
    for accessions that are not themselves canonical entries and map to exactly one entry of the reference."""
    table = pd.read_csv(history, sep='\t', dtype=str, keep_default_na=False)
    table = table[(table['kind'] == 'secondary') & table['accession'].isin(set(reference)) & ~table['key'].isin(set(reference))]
    targets = table.groupby('key')['accession'].agg(lambda v: sorted(set(v)))
    return {k: v[0] for k, v in targets.items() if len(v) == 1}


def core(window: str, center: int) -> str | None:
    """The +/-7 residues around the centre of a window, or None when the centre is not K or the window is too short."""
    window = (window or '').strip().upper()
    if not 0 <= center < len(window) or window[center] != 'K' or center < FLANK or center + FLANK >= len(window):
        return None
    return window[center - FLANK:center + FLANK + 1]


def agrees(seq: str, pos: int, motif: str) -> bool:
    """Whether the core lies on the sequence with its centre at pos (1-based): letters equal, '_' past an end."""
    if not 1 <= pos <= len(seq) or seq[pos - 1] != 'K':
        return False
    for k, ch in enumerate(motif):
        j = pos - 1 - FLANK + k
        inside = 0 <= j < len(seq)
        if ch == '_':
            if inside:
                return False
        elif not inside or seq[j] != ch:
            return False
    return True


def hits(seq: str, motif: str) -> list[int]:
    """Every position (1-based) of the sequence at which the core agrees."""
    left = len(motif) - len(motif.lstrip('_'))
    stripped = motif.strip('_')
    found, i = [], seq.find(stripped)
    while i != -1:
        pos = i + (FLANK - left) + 1
        if agrees(seq, pos, motif):
            found.append(pos)
        i = seq.find(stripped, i + 1)
    return found


def to_int(value) -> int | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return int(number) if number == int(number) else None


def group(proteins, positions, leading=None, leading_position=None) -> list[tuple[str, int | None]]:
    """The (accession, position) pairs of a row, leading protein first, from ';'-joined lists."""
    accs = [a.strip() for a in str(proteins or '').split(';')]
    poss = [to_int(p) for p in str(positions or '').split(';')]
    pairs = [(a, poss[i] if i < len(poss) else None) for i, a in enumerate(accs) if a]
    if leading:
        lead = (str(leading).strip(), to_int(leading_position))
        pairs = [lead] + [p for p in pairs if p != lead]
    seen, out = set(), []
    for p in pairs:
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def place(pairs: list[tuple[str, int | None]], motifs: list[str], reference: dict[str, str],
          secondary: dict[str, str]) -> tuple[str, str | None, int | None, str]:
    """(status, canonical accession, position, detail) of one source row: see the module docstring."""
    motifs = [m for m in motifs if m]
    if not motifs:
        return 'bad_window', None, None, 'no window with K at its centre'
    lead = pairs[0][0].split('-')[0] if pairs else ''
    ordered = [p for p in pairs if p[0].split('-')[0] == lead] + [p for p in pairs if p[0].split('-')[0] != lead]
    candidates = []
    for acc, pos in ordered:
        base = acc.split('-')[0]
        primary = base if base in reference else secondary.get(base)
        if primary is not None:
            candidates.append((acc, pos, primary))
    if not candidates:
        return 'no_canonical_accession', None, None, ';'.join(a for a, _ in pairs)[:200]
    for acc, pos, primary in candidates:      # a given position that fits
        if '-' not in acc and pos is not None and any(agrees(reference[primary], pos, m) for m in motifs):
            return ('kept' if primary == acc else 'secondary_accession'), primary, pos, '' if primary == acc else f'{acc} -> {primary}'
    not_unique = []
    for acc, pos, primary in candidates:      # else the core, once in a canonical sequence of the group
        found = sorted({p for m in motifs for p in hits(reference[primary], m)})
        if len(found) == 1:
            return 'relocated', primary, found[0], f'{acc} {pos} -> {primary} {found[0]}'
        if len(found) > 1:
            not_unique.append(f'{primary}: {len(found)} matches')
    if not_unique:
        return 'window_not_unique', None, None, '; '.join(dict.fromkeys(not_unique))[:200]
    return 'window_not_in_canonical', None, None, f'{motifs[0]} in {";".join(a for a, _ in pairs)[:160]}'

#!/usr/bin/env python3
"""Check every remapped source row against its own evidence, with code independent of scripts/02.

For each source row it derives the site its evidence supports on a reviewed UniProt entry and compares that with
where scripts/02 placed the row. Evidence is read per source: a source whose tokens mostly share one window length
(13, 15 or 31) gives centred windows, other sources give peptides, ProMetheusDB also marks the modified residue, and
the master table gives positions only.

Each row also gets a drift class, which says whether the source position still fits the current sequence:
  A  fits: the window cut around the reported position equals the source window; or the peptide covers the reported
     position and that residue is R; or (positions only) the residue is R
  B  the evidence is in the same entry, elsewhere            C  the reported accession has no sequence in --fasta
  D  the evidence is not in that entry                       E  the reported residue is not R

A row placed on the one site its evidence implies agrees. It reads 'agrees: reported arginine' when the evidence covers
other arginines of the entry too (a peptide with several, or a window found more than once), so that the reported
position picks the site, and 'agrees: position only' when the row gives no sequence.

    python scripts/qc_remap_evidence_check.py --results-root PTM_results --fasta canonical_human.fasta
    python scripts/qc_remap_evidence_check.py --results-root PTM_results_remapfix --fasta ... --extra-fasta unreviewed.fasta

A window cut short at a protein end is read from its intact end (padding with _ or X names the end); one inside a
longer entry gives the reading that lands on an arginine, the reported position deciding between two. Evidence that
is not in the reported entry still supports the reported arginine when, laid with its site there, it differs at no
more than 2 residues and no more than a fifth of them.

With --peptide-coordinates (scripts/01_parse_maron_peptide_coordinates.py) a plain peptide with several arginines also
names the arginine at the reported position counted from the peptide's start, when the table gives the start of that
peptide, or of one containing it or inside it, on the same entry.

With --fail-on-contradictions the exit status is 1 when a row is placed against its evidence, without evidence, or on an
accession outside --fasta, unless --allow lists it (TSV: source_table, reported_accession, reported_position,
placed_site, reason). Such rows are written to remap_evidence_check_failures.tsv.

Writes <outdir>/remap_evidence_check_rows.tsv and remap_evidence_check_summary.tsv and prints the summary.
"""
from __future__ import annotations

import argparse
import bisect
import gzip
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

WINDOW_LENGTHS = (13, 15, 31)
PADDING = '_-X'
NEAR_MAX_MISMATCHES, NEAR_MAX_FRACTION = 2, 0.2
FAIL_VERDICTS = ('contradicts evidence', 'placed without evidence', 'filed under an accession outside --fasta')
AGREES = {'evidence': 'agrees', 'position': 'agrees: reported arginine', 'position_only': 'agrees: position only'}
ACC = re.compile(r'([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(-\d+)?')


def read_fasta(path: str) -> tuple[dict, dict]:
    seqs, genes, acc, buf = {}, {}, None, []
    with (gzip.open if str(path).endswith('.gz') else open)(path, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if acc:
                    seqs[acc] = ''.join(buf)
                acc, buf = line[1:].split('|')[1] if '|' in line else line[1:].split()[0], []
                m = re.search(r' GN=(\S+)', line)
                genes[acc] = m.group(1).upper() if m else ''
            else:
                buf.append(line.strip())
    if acc:
        seqs[acc] = ''.join(buf)
    return seqs, genes


def letters(x) -> str:
    return re.sub(r'[^A-Z]', '', str(x).upper()) if str(x) not in ('', 'nan', 'None') else ''


def occurrences(seq: str, tok: str) -> list[int]:
    return [m.start() for m in re.finditer(f'(?={re.escape(tok)})', seq)] if seq and tok else []


class Proteome:
    """Substring search over all reviewed sequences at once."""

    def __init__(self, seqs: dict):
        self.order = sorted(seqs)
        self.offsets, parts, off = [], [], 0
        for a in self.order:
            self.offsets.append(off)
            parts.append(seqs[a] + '#')
            off += len(seqs[a]) + 1
        self.text = ''.join(parts)

    def find(self, tok: str) -> list[tuple[str, int]]:
        out, i = [], self.text.find(tok)
        while i != -1:
            k = bisect.bisect_right(self.offsets, i) - 1
            out.append((self.order[k], i - self.offsets[k]))
            i = self.text.find(tok, i + 1)
        return out


def marked_peptides(field) -> list[tuple[str, list[int]]]:
    """Marked peptides such as GVEYGSGR(di)GHLGTFEGR (ProMetheusDB, the study tables, R(methyl) where the degree is not
    given): the plain sequence and 0-based methyl-R offsets."""
    out = []
    for pep in re.split(r'[;,]', str(field)):
        plain, idx = '', []
        for m in re.finditer(r'([A-Z])(\([a-z]+\))?', pep):
            if m.group(2) in ('(me)', '(di)', '(methyl)') and m.group(1) == 'R':
                idx.append(len(plain))
            plain += m.group(1)
        if plain:
            out.append((plain, idx))
    return out


def window_lengths(rows: pd.DataFrame, field: str) -> pd.Series:
    """Full window length of each row's source for this field (0 for a peptide source)."""
    key = rows['source_table'].astype(str) + '|' + rows.get('reference', pd.Series('', index=rows.index)).fillna('').astype(str)
    toks = rows[field].map(lambda v: [letters(t) for t in str(v).split(';') if letters(t)]) if field in rows else pd.Series([[]] * len(rows), index=rows.index)
    out = {}
    for k, g in toks.groupby(key):
        flat = [len(t) for ts in g for t in ts]
        if flat:
            top, n = Counter(flat).most_common(1)[0]
            out[k] = top if top in WINDOW_LENGTHS and n / len(flat) >= 0.5 else 0
    return key.map(out).fillna(0).astype(int)


def window_token(tok, w: int) -> tuple[str, int | None]:
    """Residues of a window and the 0-based offset of its site; None for a window cut short without padding."""
    raw = re.sub(r'[^A-Z_\-]', '', str(tok).upper())
    core = raw.strip(PADDING)
    lead, trail = len(raw) - len(raw.lstrip(PADDING)), len(raw) - len(raw.rstrip(PADDING))
    if len(core) == w:
        return core, w // 2
    if (lead or trail) and lead + len(core) + trail == w:
        return core, w // 2 - lead
    return core, None


def window_sites(s: str, i: int, t: str, w: int, off: int | None, pos: int) -> list[int]:
    """0-based site supported by an occurrence of window t at offset i of s (empty if none or undecided)."""
    if off is not None:
        offs = [off]
    else:
        n_cut, c_cut = len(t) - 1 - w // 2, w // 2
        at_n, at_c = i == 0, i + len(t) == len(s)
        offs = range(max(0, n_cut), min(len(t) - 1, c_cut) + 1) if at_n and at_c else [n_cut] if at_n else [c_cut] if at_c else [n_cut, c_cut]
    sites = [i + o for o in offs if 0 <= o < len(t) and s[i + o] == 'R']
    if len(sites) > 1:
        sites = [c for c in sites if c == pos - 1 or (off is None and c - i == len(t) - 1 - w // 2 and c - i + 1 == pos)]
    return sites if len(sites) == 1 else []


def near_match(s: str, pos: int, tokens: list[tuple[str, list[int]]]) -> bool:
    """Does one token, laid with one of its arginine site offsets on pos, match s within the near-match limits?"""
    for t, offs in tokens:
        for o in offs:
            if 0 <= o < len(t) and t[o] == 'R':
                st = pos - 1 - o
                mm = sum(1 for k, aa in enumerate(t) if not 0 <= st + k < len(s) or s[st + k] != aa)
                if mm <= NEAR_MAX_MISMATCHES and mm <= NEAR_MAX_FRACTION * len(t):
                    return True
    return False


def read_coordinates(path: str | None, seqs: dict) -> dict:
    """Per entry: (peptide, 1-based start in the source numbering, 1-based start in the entry). A peptide counts when it
    lies at its stated start in the entry or occurs there once."""
    coords = defaultdict(set)
    if not path:
        return coords
    for r in pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False).itertuples():
        s, t, start = seqs.get(r.canonical_UniProtAC, ''), letters(r.peptide), int(r.peptide_start)
        at = [start - 1] if s[start - 1:start - 1 + len(t)] == t else occurrences(s, t)
        if s and t and len(at) == 1:
            coords[r.canonical_UniProtAC].add((t, start, at[0] + 1))
    return coords


def coordinate_sites(s: str, acc: str, pos: int, peptides: list[str], coords: dict) -> set[int]:
    """0-based arginines that plain peptides name when counted from a start the coordinate table gives."""
    out = set()
    for t in peptides:
        for known, src, ref in coords.get(acc, ()):
            if t in known:
                k = known.index(t)
                src_t, ref_t = src + k, ref + k
            elif known in t:
                k = t.index(known)
                src_t, ref_t = src - k, ref - k
            else:
                continue
            off = pos - src_t
            if 0 <= off < len(t) and t[off] == 'R' and s[ref_t - 1:ref_t - 1 + len(t)] == t:
                out.add(ref_t - 1 + off)
    return out


def check_row(r, seqs: dict, genes: dict, proteome: Proteome, coords: dict | None = None) -> tuple[str, set, str]:
    """Drift class, the evidence-implied sites (accession, position) on reviewed entries, and what singles the site
    out: 'evidence', 'position' (the evidence covers other arginines of the entry too, and the reported position picks
    one) or 'position_only' (no sequence)."""
    m = ACC.search(str(r.substrate_UniProtAC_original))
    acc = (m.group(1) + (m.group(2) or '')) if m else ''
    s, pos = seqs.get(acc, ''), int(r.position)
    at_r = bool(s) and 0 < pos <= len(s) and s[pos - 1] == 'R'
    windows = [(t, r.wlen_seq_window) for t in str(r.seq_window).split(';') if r.wlen_seq_window and letters(t)]
    windows += [(t, r.wlen_peptide_sequence) for t in str(r.peptide_sequence).split(';') if r.wlen_peptide_sequence and letters(t)]
    windows = [(core, off, w) for t, w in windows for core, off in [window_token(t, w)] if core and len(core) <= w]
    peptides = [letters(t) for t in str(r.peptide_sequence).split(';') if letters(t) and not r.wlen_peptide_sequence]
    marks = marked_peptides(r.peptides) if str(r.peptides) not in ('', 'nan') else []

    if not windows and not peptides and not marks:  # positions only
        return ('C' if not s else 'A' if at_r else 'E'), ({(acc, pos)} if at_r else set()), 'position_only'

    cls, implied, basis = 'C', set(), 'evidence'
    if s:
        sites = {c for t, off, w in windows for i in occurrences(s, t) for c in window_sites(s, i, t, w, off, pos)}
        sites |= {i + j for p, idx in marks for i in occurrences(s, p) for j in idx}
        covers = [(i, t) for t in peptides for i in occurrences(s, t) if i <= pos - 1 < i + len(t)]
        if pos - 1 in sites or (covers and at_r):
            cls, implied = 'A', {(acc, pos)}
            # A window or marked peptide found more than once, or a peptide with several arginines, names more than
            # the reported arginine: the reported position picks it.
            named = sites or {i + k for t in peptides for i in occurrences(s, t) for k, aa in enumerate(t) if aa == 'R'}
            basis = 'evidence' if named == {pos - 1} else 'position'
        else:
            # A peptide with a single arginine names it wherever the peptide lies.
            sites |= {i + t.index('R') for t in peptides if t.count('R') == 1 for i in occurrences(s, t)}
            sites |= coordinate_sites(s, acc, pos, peptides, coords or {})
            implied = {(acc, c + 1) for c in sites}
            found = any(occurrences(s, t) for t, _, _ in windows) or any(occurrences(s, p) for p, _ in marks) or any(occurrences(s, t) for t in peptides)
            cls = 'E' if covers else 'B' if found else 'D'
            if not found:
                tokens = [(t, [off] if off is not None else [len(t) - 1 - w // 2, w // 2]) for t, off, w in windows]
                tokens += [(p, idx) for p, idx in marks] + [(t, [k for k, aa in enumerate(t) if aa == 'R']) for t in peptides]
                if near_match(s, pos, tokens):
                    implied = {(acc, pos)}
    if acc not in genes or cls not in ('A', 'B') or not implied:
        # Reported entry unusable or not matching: search the evidence across the reviewed proteome.
        cand = set()
        for t, off, w in windows:
            for a, i in proteome.find(t):
                cand |= {(a, c + 1) for c in window_sites(seqs[a], i, t, w, off, pos)}
        for p, idx in marks:
            cand |= {(a, i + j + 1) for a, i in proteome.find(p) for j in idx}
        for t in peptides:
            src = [i for i in occurrences(s, t) if i <= pos - 1 < i + len(t)]
            for a, i in proteome.find(t):
                if src and seqs[a][i + pos - 1 - src[0]] == 'R':
                    cand.add((a, i + pos - src[0]))
                elif t.count('R') == 1:
                    cand.add((a, i + t.index('R') + 1))
        same_gene = {c for c in cand if genes.get(c[0]) == str(r.substrate_genename).upper()}
        implied = same_gene or cand or implied
    return cls, implied, basis


def read_allow(path: str | None) -> set:
    if not path:
        return set()
    allow = pd.read_csv(path, sep='\t', dtype=str, comment='#').fillna('')
    return set(zip(allow.source_table, allow.reported_accession, allow.reported_position, allow.placed_site))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results-root', default='PTM_results')
    ap.add_argument('--fasta', required=True, help='reviewed human UniProt FASTA used by the run')
    ap.add_argument('--extra-fasta', action='append', default=[], help='more sequences (unreviewed or isoform entries) for the drift class')
    ap.add_argument('--outdir', help='default: <results-root>/qc/remap_evidence_check')
    ap.add_argument('--fail-on-contradictions', action='store_true', help='exit 1 when a row is placed against or without its evidence')
    ap.add_argument('--peptide-coordinates', help='peptide starts in the source numbering (scripts/01_parse_maron_peptide_coordinates.py)')
    ap.add_argument('--allow', help='TSV of accepted exceptions (source_table, reported_accession, reported_position, placed_site, reason)')
    args = ap.parse_args()
    root = Path(args.results_root)
    outdir = Path(args.outdir) if args.outdir else root / 'qc' / 'remap_evidence_check'
    outdir.mkdir(parents=True, exist_ok=True)

    reviewed, genes = read_fasta(args.fasta)
    seqs = dict(reviewed)
    for p in args.extra_fasta:
        extra, _ = read_fasta(p)
        seqs.update({a: v for a, v in extra.items() if a not in reviewed})
    proteome = Proteome(reviewed)
    coords = read_coordinates(args.peptide_coordinates, reviewed)

    rows = pd.read_csv(root / 'remapped' / 'human_arg_methyl_source_rows_with_remap.tsv', sep='\t', low_memory=False)
    for field in ('seq_window', 'peptide_sequence', 'peptides', 'reference'):
        if field not in rows:
            rows[field] = ''
    rows['wlen_seq_window'] = window_lengths(rows, 'seq_window')
    rows['wlen_peptide_sequence'] = window_lengths(rows, 'peptide_sequence')
    kinds = rows.groupby(['source_table', rows['reference'].fillna('').astype(str).str[:40]])[['wlen_seq_window', 'wlen_peptide_sequence']].first()
    print('window length per source (0 = peptides or positions only):')
    print(kinds.to_string(), '\n')

    out = []
    for r in rows.itertuples():
        cls, implied, basis = check_row(r, seqs, genes, proteome, coords)
        placed = f'{r.corrected_canonical_UniProtAC}:{int(r.corrected_position)}' if bool(r.is_position_resolved) and pd.notna(r.corrected_position) else ''
        evidence = sorted(f'{a}:{p}' for a, p in implied if a in reviewed)
        if placed and placed.split(':')[0] not in reviewed:
            verdict = 'filed under an accession outside --fasta'
        elif not evidence:
            verdict = 'agrees: no site' if not placed else 'placed without evidence'
        elif placed in evidence:
            verdict = 'one of several evidence sites' if len(evidence) > 1 else AGREES[basis]
        elif placed:
            verdict = 'contradicts evidence'
        else:
            verdict = 'dropped, evidence has a site' if len(evidence) == 1 else 'evidence ambiguous'
        study = str(r.reference).split(' et al')[0] if str(r.reference) not in ('', 'nan') else str(r.source_table)
        out.append({'source_table': r.source_table, 'study': study, 'reported_accession': r.substrate_UniProtAC_original, 'reported_position': int(r.position),
                    'drift_class': cls, 'evidence_sites': ';'.join(evidence), 'placed_site': placed, 'remap_status': r.remap_status, 'verdict': verdict})
    res = pd.DataFrame(out)
    res.to_csv(outdir / 'remap_evidence_check_rows.tsv', sep='\t', index=False)
    summary = res.groupby(['verdict']).size().rename('rows').reset_index()
    summary.to_csv(outdir / 'remap_evidence_check_summary.tsv', sep='\t', index=False)
    print(summary.to_string(index=False), '\n')
    print('drift class by study (A = reported position fits the current sequence):')
    print(pd.crosstab(res['study'].str[:24], res['drift_class'], margins=True).to_string(), '\n')
    print('rows that contradict their evidence, by study:')
    print(res[res.verdict.isin(FAIL_VERDICTS)].groupby(['study', 'verdict']).size().to_string())
    print(f'\nrows: {outdir / "remap_evidence_check_rows.tsv"}')

    allow = read_allow(args.allow)
    keys = zip(res.source_table.astype(str), res.reported_accession.astype(str), res.reported_position.astype(str), res.placed_site.astype(str))
    failures = res[res.verdict.isin(FAIL_VERDICTS) & ~pd.Series([k in allow for k in keys], index=res.index)]
    failures.to_csv(outdir / 'remap_evidence_check_failures.tsv', sep='\t', index=False)
    allowed = int(res.verdict.isin(FAIL_VERDICTS).sum()) - len(failures)
    print(f'\nrows placed against, or without, their evidence: {len(failures)} (plus {allowed} accepted in --allow)')
    if len(failures):
        print(failures.head(40).to_string(index=False))
        print(f'all: {outdir / "remap_evidence_check_failures.tsv"}')
        if args.fail_on_contradictions:
            sys.exit(1)


if __name__ == '__main__':
    main()

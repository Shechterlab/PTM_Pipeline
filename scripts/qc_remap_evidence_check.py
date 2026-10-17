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

    python scripts/qc_remap_evidence_check.py --results-root PTM_results --fasta canonical_human.fasta
    python scripts/qc_remap_evidence_check.py --results-root PTM_results_remapfix --fasta ... --extra-fasta unreviewed.fasta

Writes <outdir>/remap_evidence_check_rows.tsv and remap_evidence_check_summary.tsv and prints the summary.
"""
from __future__ import annotations

import argparse
import bisect
import gzip
import re
from collections import Counter
from pathlib import Path

import pandas as pd

WINDOW_LENGTHS = (13, 15, 31)
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
    """ProMetheusDB peptides such as GVEYGSGR(di)GHLGTFEGR: the plain sequence and 0-based methyl-R offsets."""
    out = []
    for pep in re.split(r'[;,]', str(field)):
        plain, idx = '', []
        for m in re.finditer(r'([A-Z])(\([a-z]+\))?', pep):
            if m.group(2) in ('(me)', '(di)') and m.group(1) == 'R':
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


def check_row(r, seqs: dict, genes: dict, proteome: Proteome) -> tuple[str, set]:
    """Drift class and the evidence-implied sites (accession, position) on reviewed entries."""
    m = ACC.search(str(r.substrate_UniProtAC_original))
    acc = (m.group(1) + (m.group(2) or '')) if m else ''
    s, pos = seqs.get(acc, ''), int(r.position)
    at_r = bool(s) and 0 < pos <= len(s) and s[pos - 1] == 'R'
    windows = [(t, r.wlen_seq_window) for t in str(r.seq_window).split(';') if r.wlen_seq_window and letters(t)]
    windows += [(t, r.wlen_peptide_sequence) for t in str(r.peptide_sequence).split(';') if r.wlen_peptide_sequence and letters(t)]
    windows = [(letters(t), w) for t, w in windows if len(letters(t)) <= w]
    peptides = [letters(t) for t in str(r.peptide_sequence).split(';') if letters(t) and not r.wlen_peptide_sequence]
    marks = marked_peptides(r.peptides) if str(r.peptides) not in ('', 'nan') else []

    if not windows and not peptides and not marks:  # positions only
        return ('C' if not s else 'A' if at_r else 'E'), ({(acc, pos)} if at_r else set())

    cls, implied = 'C', set()
    if s:
        fits = any(t == s[max(0, pos - 1 - w // 2):min(len(s), pos + w // 2)] for t, w in windows) and at_r
        fits = fits or any(i + j == pos - 1 for p, idx in marks for i in occurrences(s, p) for j in idx)
        covers = [(i, t) for t in peptides for i in occurrences(s, t) if i <= pos - 1 < i + len(t)]
        if fits or (covers and at_r):
            cls, implied = 'A', {(acc, pos)}
        elif covers:
            cls = 'E'
        else:
            for t, w in windows:
                for i in occurrences(s, t):
                    c = i + w // 2 if len(t) == w or i + len(t) == len(s) else i + len(t) - 1 - w // 2 if i == 0 else None
                    if c is not None and 0 <= c < len(s) and s[c] == 'R':
                        implied.add((acc, c + 1))
            implied |= {(acc, i + j + 1) for p, idx in marks for i in occurrences(s, p) for j in idx}
            found = any(occurrences(s, t) for t, _ in windows) or any(occurrences(s, p) for p, _ in marks) or any(occurrences(s, t) for t in peptides)
            cls = 'B' if found else 'D'
    if acc not in genes or cls not in ('A', 'B') or not implied:
        # Reported entry unusable: search the evidence across the reviewed proteome.
        cand = set()
        for t, w in windows:
            if len(t) == w:
                cand |= {(a, i + w // 2 + 1) for a, i in proteome.find(t) if seqs[a][i + w // 2] == 'R'}
        for p, idx in marks:
            cand |= {(a, i + j + 1) for a, i in proteome.find(p) for j in idx}
        for t in peptides:
            src = [i for i in occurrences(s, t) if i <= pos - 1 < i + len(t)]
            for a, i in proteome.find(t):
                if src and seqs[a][i + pos - 1 - src[0]] == 'R':
                    cand.add((a, i + pos - src[0]))
                elif not src and t.count('R') == 1:
                    cand.add((a, i + t.index('R') + 1))
        same_gene = {c for c in cand if genes.get(c[0]) == str(r.substrate_genename).upper()}
        implied = same_gene or cand
    return cls, implied


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results-root', default='PTM_results')
    ap.add_argument('--fasta', required=True, help='reviewed human UniProt FASTA used by the run')
    ap.add_argument('--extra-fasta', action='append', default=[], help='more sequences (unreviewed or isoform entries) for the drift class')
    ap.add_argument('--outdir', help='default: <results-root>/qc/remap_evidence_check')
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
        cls, implied = check_row(r, seqs, genes, proteome)
        placed = f'{r.corrected_canonical_UniProtAC}:{int(r.corrected_position)}' if bool(r.is_position_resolved) and pd.notna(r.corrected_position) else ''
        evidence = sorted(f'{a}:{p}' for a, p in implied if a in reviewed)
        if len(evidence) > 1:
            verdict = 'one of several evidence sites' if placed in evidence else 'evidence ambiguous'
        elif not evidence:
            verdict = 'agrees: no site' if not placed else 'placed without evidence'
        else:
            verdict = 'agrees' if placed == evidence[0] else 'dropped, evidence has a site' if not placed else 'contradicts evidence'
        if placed and placed.split(':')[0] not in reviewed:
            verdict = 'filed under an accession outside --fasta'
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
    print(res[res.verdict.isin(['contradicts evidence', 'filed under an accession outside --fasta'])].groupby(['study', 'verdict']).size().to_string())
    print(f'\nrows: {outdir / "remap_evidence_check_rows.tsv"}')


if __name__ == '__main__':
    main()

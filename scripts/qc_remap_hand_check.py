#!/usr/bin/env python3
"""Print the evidence behind a protein's remapped methylarginine sites, for checking by eye.

For every source row of the protein it shows the study, the accession and position the source reported, the
residue UniProt has at that position, where the row's peptide or window occurs in the canonical sequence, and where
the pipeline placed the site. The sequence around the sites is printed with the reported (^) and placed (*)
positions marked. Reads the committed remap and union tables; sequences come from UniProt (or --fasta).

    python scripts/qc_remap_hand_check.py SNRPD1            # gene symbol or UniProt accession
    python scripts/qc_remap_hand_check.py examples          # the cases behind the 2026-10 remap fixes
    python scripts/qc_remap_hand_check.py P62314 --results-root /path/to/PTM_results_remapfix
    python scripts/qc_remap_hand_check.py --audit-union     # is every union site an arginine in UniProt?
    python scripts/qc_remap_hand_check.py --audit-union --fasta <the run's FASTA>   # exit 1 if not (pipeline gate)
"""
from __future__ import annotations

import argparse
import re
import sys
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

import pandas as pd

EXAMPLES = [
    ('SNRPD1', 'Peptides placed on the arginine nearest their middle: Maron reports R104-R119, the union keeps R98-R106 only.'),
    ('A2A274', 'Position computed on the reviewed ACO2 entry (Q99798) but filed under the unreviewed accession A2A274.'),
    ('Q93074', 'A PTCD1 row filed under MED12 (Q93074): its window belongs to PTCD1 (O75127), so MED12 "R235" is not a reported site.'),
    ('TTC28', 'Musiani rows whose 15-residue window is centred on a lysine; the pipeline made an arginine site from the nearest R.'),
    ('MRE11', 'Script 08 tested MRE11 as two proteins (MRE11, MRE11A); see per_protein_arg_odds.tsv.'),
]


def fetch_sequence(acc: str, cache: dict, fasta: dict | None) -> str:
    if acc in cache:
        return cache[acc]
    seq = ''
    if fasta is not None:
        seq = fasta.get(acc, '')
    else:
        try:
            with urllib.request.urlopen(f'https://rest.uniprot.org/uniprotkb/{acc}.fasta', timeout=30) as r:
                seq = ''.join(line.strip() for line in r.read().decode().splitlines()[1:])
        except Exception:  # noqa: BLE001 - obsolete or merged accessions return nothing
            seq = ''
    cache[acc] = seq
    return seq


def read_fasta(path: str) -> dict:
    out, acc, buf = {}, None, []
    for line in open(path):
        if line.startswith('>'):
            if acc:
                out[acc] = ''.join(buf)
            acc, buf = line[1:].split('|')[1] if '|' in line else line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if acc:
        out[acc] = ''.join(buf)
    return out


def tokens(row) -> list[str]:
    out = []
    for field in ('peptides', 'seq_window', 'peptide_sequence'):
        for t in str(row.get(field, '') or '').split(';'):
            t = re.sub(r'[^A-Za-z]', '', re.sub(r'\[.*?\]|\(.*?\)', '', t)).upper()
            if t and t != 'NAN':
                out.append(t)
    return out


def show(query: str, rows: pd.DataFrame, union: pd.DataFrame, root: Path, cache: dict, fasta: dict | None) -> None:
    q = query.upper()
    hit = rows[(rows['substrate_UniProtAC'].astype(str).str.upper() == q) | (rows['canonical_UniProtAC'].astype(str).str.upper() == q)
               | (rows['corrected_canonical_UniProtAC'].astype(str).str.upper() == q) | (rows['substrate_genename'].astype(str).str.upper() == q)]
    if hit.empty:
        print(f'No source rows for {query}.')
        return
    accs = sorted(set(hit['canonical_UniProtAC'].dropna().astype(str)) | set(hit['corrected_canonical_UniProtAC'].dropna().astype(str))
                  | {a for c in hit['candidate_canonical_accessions'].dropna().astype(str) for a in c.split(';') if a})
    print(f'Accessions involved: {", ".join(accs)}')
    for acc in accs:
        u = sorted(union[union['canonical_UniProtAC'] == acc]['corrected_position'].astype(int))
        seq = fetch_sequence(acc, cache, fasta)
        print(f'  {acc}: {len(seq) or "not in UniProt (obsolete or merged)"}{" aa" if seq else ""}; union sites: {", ".join(f"R{p}" for p in u) or "none"}')
    print()
    print(f'{"study":<18} {"src acc":<11} {"gene":<7} {"rep":>5} {"res@rep":>7}  {"peptide/window found at":<34} {"placed":<17} {"res@placed":>10}  status')
    marks = {}
    for _, r in hit.sort_values(['position']).iterrows():
        acc = str(r['canonical_UniProtAC'])
        seq = fetch_sequence(acc, cache, fasta)
        rep = int(r['position'])
        res = seq[rep - 1] if seq and 0 < rep <= len(seq) else '-'
        found = []
        cands = [a for a in str(r.get('candidate_canonical_accessions', '') or '').split(';') if a and a != 'nan']
        for t in tokens(r):
            for target in dict.fromkeys([acc, str(r.get('corrected_canonical_UniProtAC', '') or '')] + cands):
                s = fetch_sequence(target, cache, fasta) if target and target != 'nan' else ''
                for m in re.finditer(f'(?={re.escape(t)})', s):
                    found.append(f'{target}:{m.start() + 1}-{m.start() + len(t)}')
        placed_acc = str(r.get('corrected_canonical_UniProtAC', '') or '')
        placed = f'{placed_acc}:R{int(r["corrected_position"])}' if pd.notna(r.get('corrected_position')) else 'dropped'
        pseq = fetch_sequence(placed_acc, cache, fasta) if placed != 'dropped' else ''
        pres = pseq[int(r['corrected_position']) - 1] if pseq and int(r['corrected_position']) <= len(pseq) else '-'
        study = r.get('reference')
        study = str(study) if isinstance(study, str) and study else str(r.get('source', ''))
        study = 'Maron 2021' if study.startswith('This study') else study.split(' et al.')[0].split(',')[0][:17]
        print(f'{study:<18} {str(r["substrate_UniProtAC"]):<11} {str(r["substrate_genename"])[:7]:<7} {rep:>5} {res:>7}  {";".join(sorted(set(found)))[:34]:<34} {placed:<17} {pres:>10}  {r["remap_status"]}')
        marks.setdefault(acc, set()).add(('^', rep))
        if pd.notna(r.get('corrected_position')):
            marks.setdefault(placed_acc, set()).add(('*', int(r['corrected_position'])))
    # Sequence ruler with reported (^) and placed (*) positions.
    for acc, ms in marks.items():
        seq = fetch_sequence(acc, cache, fasta)
        if not seq:
            continue
        lo = max(1, min(p for _, p in ms) - 5)
        hi = min(len(seq), max(p for _, p in ms) + 5)
        if hi - lo > 120:
            continue
        print(f'\n{acc} {lo}-{hi}   ^ reported by a source   * where the pipeline placed a site')
        print(f'{"":>7}{seq[lo - 1:hi]}')
        line = [' '] * (hi - lo + 1)
        for kind, p in ms:
            if lo <= p <= hi:
                i = p - lo
                line[i] = '+' if line[i] not in (' ', kind) else kind
        print(f'{lo:>6} {"".join(line)}')
    odds = root / 'motif_arg_odds' / 'per_protein_arg_odds.tsv'
    if odds.exists():
        t = pd.read_csv(odds, sep='\t')
        sub = t[t['canonical_UniProtAC'].isin(accs)]
        if len(sub) > 1 or (len(sub) and sub['canonical_UniProtAC'].duplicated().any()):
            print('\nper_protein_arg_odds.tsv rows (one protein should have one row):')
            print(sub[['canonical_UniProtAC', 'substrate_genename', 'methyl_site_count', 'arginine_count', 'q_value']].to_string(index=False))
    print()


def fetch_many(accs: list[str], cache: dict) -> None:
    """Fill the cache with current UniProt sequences, 80 accessions per request; missing accessions stay ''."""
    todo = [a for a in accs if a not in cache]
    for i in range(0, len(todo), 80):
        q = ' OR '.join(f'accession:{a}' for a in todo[i:i + 80])
        url = 'https://rest.uniprot.org/uniprotkb/stream?' + urllib.parse.urlencode({'query': q, 'format': 'fasta'})
        with urllib.request.urlopen(url, timeout=120) as r:
            acc, buf = None, []
            for line in r.read().decode().splitlines():
                if line.startswith('>'):
                    if acc:
                        cache[acc] = ''.join(buf)
                    acc, buf = line[1:].split('|')[1], []
                else:
                    buf.append(line.strip())
            if acc:
                cache[acc] = ''.join(buf)
    for a in todo:
        cache.setdefault(a, '')


def audit_union(union: pd.DataFrame, cache: dict, fasta: dict | None, out: Path) -> int:
    """Check every union site against the sequences; returns the number of sites that are not an arginine."""
    accs = sorted(set(union['canonical_UniProtAC'].astype(str)))
    if fasta is None:
        fetch_many(accs, cache)
    rows, counts = [], Counter()
    for acc, pos in zip(union['canonical_UniProtAC'].astype(str), union['corrected_position'].astype(int)):
        seq = fasta.get(acc, '') if fasta is not None else cache.get(acc, '')
        if not seq:
            kind = 'accession not in --fasta' if fasta is not None else 'accession not in current UniProt'
        elif pos > len(seq):
            kind = 'position beyond the sequence end'
        elif seq[pos - 1] != 'R':
            kind = f'residue is {seq[pos - 1]}, not R'
        else:
            counts['arginine'] += 1
            continue
        counts[kind.split(',')[0] if kind.startswith('residue') else kind] += 1
        rows.append({'accession': acc, 'position': pos, 'problem': kind,
                     'context': seq[max(0, pos - 8):pos + 7] if seq and pos <= len(seq) else ''})
    print(f'{len(union):,} union sites checked against {"the --fasta" if fasta is not None else "current UniProt"} sequences:')
    for k, v in counts.most_common():
        print(f'  {v:7,}  {k}')
    pd.DataFrame(rows, columns=['accession', 'position', 'problem', 'context']).to_csv(out, sep='\t', index=False)
    print(f'sites that are not an arginine, or cannot be checked: {out}')
    return len(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('query', nargs='*', help='gene symbols or UniProt accessions, or "examples"')
    ap.add_argument('--audit-union', action='store_true', help='check that every site of the union is an arginine in UniProt')
    ap.add_argument('--audit-out', default='union_site_residue_audit.tsv')
    ap.add_argument('--results-root', default='PTM_results')
    ap.add_argument('--fasta', help='use this UniProt FASTA instead of fetching sequences')
    args = ap.parse_args()
    root = Path(args.results_root)
    rows = pd.read_csv(root / 'remapped' / 'human_arg_methyl_source_rows_with_remap.tsv', sep='\t', low_memory=False)
    union = pd.read_csv(root / 'integrated' / 'human_arg_methyl_union_dedup_by_site.tsv', sep='\t', low_memory=False)
    fasta = read_fasta(args.fasta) if args.fasta else None
    cache: dict = {}
    if args.audit_union:
        problems = audit_union(union, cache, fasta, Path(args.audit_out))
        if not args.query:
            sys.exit(1 if problems else 0)
    queries = [(g, note) for g, note in EXAMPLES] if args.query == ['examples'] else [(q, '') for q in args.query]
    for q, note in queries:
        print('=' * 100)
        print(q + (f': {note}' if note else ''))
        print('=' * 100)
        show(q, rows, union, root, cache, fasta)
    sys.stdout.flush()


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""List the sites a release keeps pending a policy decision: windows that match the entry on one side of the site only.

scripts/02_remap_arg_methyl_to_canonical.py leaves a window record unplaced when its window does not match the entry
at the reported arginine. Some of these windows match one full flank (at least 6 residues) and the site exactly, but
not the other flank in full, and no other reviewed entry holds the window. Some are windows of another isoform of the
protein (a different exon after the site). A window whose numbering is that of another isoform as well is kept the
same way at the arginine where one complete flank (15 residues and the arginine) matches the entry, when that happens
at exactly one arginine (FMR1 R606: a window reported at R564, from an isoform that skips an exon). The release keeps
these sites, flagged one_sided_window, until the policy is decided.

    python scripts/02_find_one_sided_window_sites.py --results-root WORK --fasta inputs/reference/....fasta

Writes WORK/integrated/pending_one_sided_window_sites.tsv: one row per site (accession, position) with the source records
(record_id = row number of remapped/human_arg_methyl_source_rows_with_remap.tsv, from 1).
"""
from __future__ import annotations

import argparse
import importlib.util
from collections import defaultdict
from pathlib import Path

import pandas as pd

SCRIPTS = Path(__file__).resolve().parent


def load_script(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


remap = load_script('remap', '02_remap_arg_methyl_to_canonical.py')
MIN_FLANK = 6


def one_sided_window_sites(typed: pd.DataFrame, reference: dict, union_keys: set) -> dict[tuple[str, int], list[int]]:
    proteome = remap.ProteomeIndex(reference)
    out: dict[tuple[str, int], list[int]] = defaultdict(list)
    for i in typed.index[typed['remap_status'] == 'unresolved_evidence_mismatch']:
        r = typed.loc[i]
        acc = remap.canonicalize_uniprot_accession(r['canonical_UniProtAC'])
        seq, pos = reference.get(acc, ''), int(r['position'])
        if not (seq and 0 < pos <= len(seq) and seq[pos - 1] == 'R') or (acc, pos) in union_keys:
            continue
        for item in remap.row_evidence(r):
            if item['anchor_type'] not in ('centered_window', 'truncated_window'):
                continue
            t = item['sequence']
            if proteome.find(t):
                continue
            offsets = item['offsets'] if item['offsets'] is not None else [len(t) - item['full_length'] // 2 - 1, item['full_length'] // 2]
            for o in offsets:
                if not (0 <= o < len(t) and t[o] == 'R'):
                    continue
                start = pos - 1 - o
                ok = [0 <= start + k < len(seq) and seq[start + k] == aa for k, aa in enumerate(t)]
                if min(o, len(t) - 1 - o) >= MIN_FLANK and (all(ok[:o + 1]) or all(ok[o:])):
                    out[(acc, pos)].append(i)
                    break
            else:
                if item['anchor_type'] == 'centered_window':
                    elsewhere = flank_elsewhere(t, seq)
                    if elsewhere and (acc, elsewhere) not in union_keys:
                        out[(acc, elsewhere)].append(i)
                        break
                continue
            break
    return dict(out)


def flank_elsewhere(t: str, seq: str) -> int | None:
    """The one arginine of seq where a complete flank of the window t (15 residues and the central arginine) matches."""
    o = len(t) // 2
    if len(t) != 31 or t[o] != 'R':
        return None
    hits = set()
    for frag, shift in ((t[:o + 1], o), (t[o:], 0)):
        start = seq.find(frag)
        while start != -1:
            hits.add(start + 1 + shift)
            start = seq.find(frag, start + 1)
    return hits.pop() if len(hits) == 1 else None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results-root', required=True)
    ap.add_argument('--fasta', required=True)
    args = ap.parse_args()
    root = Path(args.results_root)
    reference = remap.parse_fasta(args.fasta)
    typed = pd.read_csv(root / 'remapped' / 'human_arg_methyl_source_rows_with_remap.tsv', sep='\t', low_memory=False)
    typed = remap.annotate_window_lengths(typed) if 'seq_window_window_length' not in typed else typed
    union = pd.read_csv(root / 'integrated' / 'human_arg_methyl_union_dedup_by_site.tsv', sep='\t', low_memory=False)
    union_keys = {(a, int(float(p))) for a, p in zip(union['canonical_UniProtAC'], union['corrected_position'])}
    found = one_sided_window_sites(typed, reference, union_keys)
    table = pd.DataFrame([{'canonical_UniProtAC': a, 'position': p, 'record_ids': '; '.join(str(i + 1) for i in ids)}
                          for (a, p), ids in sorted(found.items())], columns=['canonical_UniProtAC', 'position', 'record_ids'])
    table.to_csv(root / 'integrated' / 'pending_one_sided_window_sites.tsv', sep='\t', index=False, lineterminator='\n')
    print(f'{len(table)} pending sites from {sum(len(v) for v in found.values())} records')


if __name__ == '__main__':
    main()

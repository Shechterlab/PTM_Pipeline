#!/usr/bin/env python3
"""Preview: do methylated proteins carry disordered acidic tracts, and do their methylarginines sit near them?

    python scripts/52_preview_acidic_tracts.py [--sites releases/v1.4.1/human_methylarginine_sites_v1.4.1.tsv]
        [--min-runs 5,6,8] [--nears 10,30,60,100] [--out preview.tsv]

Reads only staged inputs: the reference proteome, MobiDB disorder intervals, PaxDb abundance and a release's sites.

An acidic tract is an uninterrupted run of at least --min-run D or E (the rule of the MicroDomains prospective list),
called disordered when at least half of its residues lie in a MobiDB disorder interval. Only disordered tracts count.

  (A) proteins: methylated (a site in --sites) against the other reviewed proteins, outcome = at least one disordered
      acidic tract. Fisher's odds ratio (0.5 added to each cell), and the Mantel-Haenszel odds ratio with each PaxDb
      abundance decile a stratum (plus one for proteins without a PaxDb value), as scripts/abundance_strata.py does.
  (B) residues, in cis: the disordered arginines of methylated proteins that carry a disordered tract, exposure = the
      arginine is a methylarginine, outcome = within --near residues of a disordered tract. Mantel-Haenszel with each
      protein a stratum, its Robins-Breslow-Greenland CI, and p_rotation (the methylation pattern rotated along each
      protein's arginines), from scripts/43c_preview_rbrid_methylarginines.py.
  (B') (B) without arginines in RG or GR dipeptides (methylated or not), since PRMT1 prefers RG context.
Every combination of --min-runs and --nears is run; (A) depends on the tract length only.
"""
from __future__ import annotations

import argparse
import importlib.util
import re
from pathlib import Path

import numpy as np
import pandas as pd

from common import mantel_haenszel_ci, parse_fasta

ROOT = Path(__file__).resolve().parents[1]


def load_rbrid_preview():
    spec = importlib.util.spec_from_file_location('preview_rbrid', Path(__file__).with_name('43c_preview_rbrid_methylarginines.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def disorder_masks(path: Path, seqs: dict[str, str]) -> dict[str, np.ndarray]:
    dis = pd.read_csv(path, sep='\t')
    out = {}
    for acc, g in dis.groupby('canonical_UniProtAC'):
        if acc not in seqs:
            continue
        m = np.zeros(len(seqs[acc]) + 1, bool)  # 1-based
        for a, b in zip(g.fragment_start, g.fragment_end):
            m[max(1, a):min(len(seqs[acc]), b) + 1] = True
        out[acc] = m
    return out


def disordered_tracts(seqs: dict[str, str], masks: dict[str, np.ndarray], min_run: int) -> dict[str, list[tuple[int, int]]]:
    out = {}
    for acc, s in seqs.items():
        m = masks.get(acc)
        out[acc] = [(t.start() + 1, t.end()) for t in re.finditer('[DE]{%d,}' % min_run, s)
                    if m is not None and m[t.start() + 1:t.end() + 1].mean() >= 0.5]
    return out


def protein_test(seqs, tracts, methylated: set[str], percentile: pd.Series) -> dict:
    df = pd.DataFrame({'acc': list(seqs)})
    df['tract'] = df.acc.map(lambda a: bool(tracts[a]))
    df['methylated'] = df.acc.isin(methylated)
    pct = df.acc.map(percentile)
    df['stratum'] = pd.qcut(pct, 10, labels=False, duplicates='drop').fillna(-1)

    def cells(g):
        return (int((g.methylated & g.tract).sum()), int((g.methylated & ~g.tract).sum()),
                int((~g.methylated & g.tract).sum()), int((~g.methylated & ~g.tract).sum()))
    a, b, c, d = cells(df)
    odds, lo, hi, p, used = mantel_haenszel_ci([cells(g) for _, g in df.groupby('stratum')])
    return {'test': 'A_proteins', 'a': a, 'b': b, 'c': c, 'd': d,
            'target_fraction': a / (a + b), 'background_fraction': c / (c + d),
            'or_raw': (a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5)),
            'or_mh': odds, 'ci_low': lo, 'ci_high': hi, 'p': p, 'p_rotation': np.nan, 'strata': used}


def residue_rows(seqs, masks, tracts, sites: dict[str, set[int]], near: int) -> pd.DataFrame:
    rows = []
    for acc in sorted(sites):
        s, m, dt = seqs.get(acc), masks.get(acc), tracts.get(acc)
        if not s or m is None or not dt:
            continue
        for i, ch in enumerate(s, 1):
            if ch != 'R' or not m[i]:
                continue
            rg = (i > 1 and s[i - 2] == 'G') or (i < len(s) and s[i] == 'G')
            rows.append((acc, i, i in sites[acc], any(x - near <= i <= y + near for x, y in dt), rg))
    return pd.DataFrame(rows, columns=['acc', 'position', 'methylated', 'near_tract', 'rg_gr'])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--fasta', type=Path, default=ROOT / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03.fasta')
    ap.add_argument('--disorder', type=Path, default=ROOT / 'inputs' / 'annotation' / 'mobidb_disorder_intervals_2026-03-30.tsv')
    ap.add_argument('--paxdb', type=Path, default=ROOT / 'inputs' / 'layers' / 'paxdb_human_integrated_abundance.tsv')
    ap.add_argument('--sites', type=Path, default=ROOT / 'releases' / 'v1.4.1' / 'human_methylarginine_sites_v1.4.1.tsv')
    ap.add_argument('--min-runs', default='5,6,8')
    ap.add_argument('--nears', default='10,30,60,100')
    ap.add_argument('--rotations', type=int, default=2000)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out', type=Path)
    args = ap.parse_args()

    mh_rot = load_rbrid_preview().mantel_haenszel
    seqs = parse_fasta(args.fasta)
    masks = disorder_masks(args.disorder, seqs)
    percentile = pd.read_csv(args.paxdb, sep='\t').set_index('canonical_UniProtAC')['percentile']
    st = pd.read_csv(args.sites, sep='\t', usecols=['uniprot_accession', 'position'])
    sites = st.groupby('uniprot_accession')['position'].apply(lambda x: set(map(int, x))).to_dict()

    out = []
    for min_run in map(int, args.min_runs.split(',')):
        tracts = disordered_tracts(seqs, masks, min_run)
        n_tracts = sum(len(v) for v in tracts.values())
        res = protein_test(seqs, tracts, set(sites), percentile)
        out.append({'min_run': min_run, 'near': np.nan, 'disordered_tracts': n_tracts, **res})
        for near in map(int, args.nears.split(',')):
            r = residue_rows(seqs, masks, tracts, sites, near)
            for test, sub in [('B_cis', r), ('B_cis_without_RG_GR', r[~r.rg_gr])]:
                sub = sub.reset_index(drop=True)
                res = mh_rot(sub.acc.to_numpy(), sub.methylated.to_numpy(), sub.near_tract.to_numpy(),
                             n_perm=args.rotations, seed=args.seed)
                a, b, c, d = res['a'], res['b'], res['c'], res['d']
                out.append({'min_run': min_run, 'near': near, 'disordered_tracts': n_tracts, 'test': test,
                            'a': a, 'b': b, 'c': c, 'd': d, 'target_fraction': a / (a + b), 'background_fraction': c / (c + d),
                            'or_mh': res['or_mh'], 'ci_low': res['ci_low'], 'ci_high': res['ci_high'],
                            'p': res['p'], 'p_rotation': res['p_rotation'], 'strata': res['informative_proteins']})
    df = pd.DataFrame(out)
    with pd.option_context('display.width', 200, 'display.max_columns', 30):
        print(df.to_string(index=False, float_format=lambda x: f'{x:.3g}'))
    if args.out:
        df.to_csv(args.out, sep='\t', index=False)


if __name__ == '__main__':
    main()

"""Pilot: methylated vs other proteins the Human Cell Map (Go et al. 2021) localises, per NMF compartment and per SAFE
domain; one-sided hypergeometric with fold enrichment (as the release's enricher tests), and the Mantel-Haenszel odds
ratio within PaxDb abundance deciles of the same universe; beside it the same sets against the reviewed proteome."""
import sys, math, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import hypergeom
repo, dl, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
sys.path.insert(0, str(repo / 'scripts'))
import abundance_strata as A
from common import mantel_haenszel_ci, benjamini_hochberg

entries = pd.read_csv(repo / 'inputs/reference/uniprot_human_reviewed_2026_03_entries.tsv', sep='\t', dtype=str, keep_default_na=False)
reviewed = set(entries['Entry'])
by_symbol = {}
for e, g in zip(entries['Entry'], entries['Gene Names (primary)']):
    for s in g.split(';'):
        s = s.strip()
        if s: by_symbol.setdefault(s, set()).add(e)
prot = pd.read_csv(repo / 'releases/v1.4.1/human_methylarginine_proteins_v1.4.1.tsv', sep='\t', dtype=str, keep_default_na=False)
meth = {a.split('-')[0] for a in prot['uniprot_accession']}
pax = A.abundance(repo / 'inputs')['protein']

p = pd.read_csv(dl / 'preys-v1.txt', sep='\t', dtype=str, keep_default_na=False)
def acc_of(row):
    a = row['UniProt'].split('-')[0]
    if a in reviewed: return a
    hits = by_symbol.get(row['symbol'], set())
    return next(iter(hits)) if len(hits) == 1 else ''
p['acc'] = p.apply(acc_of, axis=1)
lost = p[p['acc'] == '']
p = p[p['acc'] != ''].drop_duplicates('acc')
print(f'preys {len(p) + len(lost)}; on a reviewed entry {len(p)}; not mapped {len(lost)} (e.g. {list(lost.symbol[:6])})')
ns = pd.read_csv(dl / 'not-significant-v1.txt', sep='\t', dtype=str, keep_default_na=False)
ns_acc = {next(iter(by_symbol[s])) for s in ns['prey'] if len(by_symbol.get(s, ())) == 1}
detected = set(p['acc']) | ns_acc
print(f'detected (significant + not significant preys, on reviewed entries): {len(detected)}; localised: {len(p)}')
print(f'methylated: of localised {len(set(p.acc) & meth)} ({len(set(p.acc) & meth)/len(p):.1%}); of detected {len(detected & meth)} ({len(detected & meth)/len(detected):.1%}); of reviewed proteome {len(meth & reviewed)}/{len(reviewed)} ({len(meth & reviewed)/len(reviewed):.1%})')

def tests(universe, sets, label):
    universe = sorted(universe)
    U = np.array(universe); is_t = np.array([u in meth for u in universe]); n, N = int(is_t.sum()), len(universe)
    lay = A.strata(universe, pax).to_numpy(dtype=str); codes, layers = np.unique(lay, return_inverse=True)
    rows = []
    for name, members in sets.items():
        ins = np.array([u in members for u in universe]); M = int(ins.sum()); k = int((ins & is_t).sum())
        if M < 5: continue
        pval = float(hypergeom.sf(k - 1, N, M, n)) if k else 1.0
        fe = (k / n) / (M / N) if n and M else float('nan')
        tables = []
        for i in range(len(codes)):
            s = layers == i
            a = int((s & ins & is_t).sum()); b = int((s & is_t & ~ins).sum()); c = int((s & ins & ~is_t).sum()); d = int((s & ~ins & ~is_t).sum())
            tables.append((a, b, c, d))
        odds, lo, hi, pmh, used = mantel_haenszel_ci(tables)
        raw = (k * (N - M - n + k)) / ((n - k) * (M - k)) if (n - k) and (M - k) else float('inf')
        rows.append(dict(background=label, compartment=name, members=M, methylated=k, share=k / M, fold=fe, p=pval, odds_raw=raw, odds_mh=odds, lo_mh=lo, hi_mh=hi, p_mh=pmh))
    t = pd.DataFrame(rows)
    t['q'] = benjamini_hochberg(t['p']); t['q_mh'] = benjamini_hochberg(t['p_mh'].fillna(1))
    return t

nmf = {lab: set(g['acc']) for lab, g in p.groupby('MMF localization')}
safe = {lab: set(g['acc']) for lab, g in p.groupby('SAFE localization') if lab not in ('no prediction', '-')}
res = pd.concat([tests(set(p.acc), nmf, 'localised (NMF)'), tests(reviewed, nmf, 'reviewed proteome (NMF)'),
                 tests(set(p.loc[~p['SAFE localization'].isin(['no prediction', '-']), 'acc']), safe, 'localised (SAFE)')], ignore_index=True)
res.to_csv(out / 'hcm_pilot.tsv', sep='\t', index=False)
pd.set_option('display.width', 220)
for bg in res.background.unique():
    t = res[res.background == bg].sort_values('fold', ascending=False)
    print('\n==', bg)
    print(t[['compartment', 'members', 'methylated', 'share', 'fold', 'q', 'odds_raw', 'odds_mh', 'q_mh']].to_string(index=False, float_format=lambda v: f'{v:.3g}'))
# per-protein placement for the protein pages
json.dump({r.acc: [r['MMF localization'], r['SAFE localization']] for _, r in p.iterrows()}, open(out / 'hcm_protein.json', 'w'))

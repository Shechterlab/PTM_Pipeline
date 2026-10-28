"""Pilot: the speckle, lamina and PML-body proximity proteomes of Barutcu et al. 2022 (Mol Cell 82:1035; APEX2-MS in
Flp-In HEK293, Table S1: SAINT scores of every prey with each APEX2 marker; raw data MassIVE MSV000087649).

A protein belongs to nuclear speckles when SAINT's average probability (AvgP) is >= 0.9 with at least two of the three
speckle markers (SRSF1, SRSF7, RNPS1), to the nuclear lamina when it is >= 0.9 with LMNA or LMNB1, and to PML bodies
when it is >= 0.9 with SP100. The proteins the map places are every prey the experiment detected. Methylated proteins
are compared with the other detected proteins (lib.tests: hypergeometric, and within PaxDb abundance deciles).

    python pilot_apexms.py REPO mmc2.xlsx OUT     -> OUT/apexms_pilot.tsv, OUT/apexms_protein.json
"""
import sys, json
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
import lib

repo, table, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
ctx = lib.setup(repo)
t = pd.read_excel(table, header=1)
t['acc'] = [lib.resolve(ctx, [], g) for g in t['Prey Gene'].astype(str)]
lost = t.loc[t['acc'] == '', 'Prey Gene'].nunique()
t = t[t['acc'] != '']
sure = t[t['AvgP'] >= 0.9]
by = lambda baits: sure[sure['Bait'].isin(baits)].groupby('acc')['Bait'].nunique()  # noqa: E731
speck = by(['SRSF1', 'SRSF7', 'RNPS1'])
sets = {'Nuclear speckles': set(speck[speck >= 2].index), 'Nuclear lamina': set(by(['LMNA', 'LMNB1']).index), 'PML bodies': set(by(['SP100']).index)}
universe = set(t['acc'])
print(f'preys on a reviewed entry {len(universe)} (unmapped symbols {lost}); ' + ', '.join(f'{k} {len(v)}' for k, v in sets.items()))
res = lib.tests(ctx, universe, sets)
lib.show(res, 'APEX2-MS (Barutcu 2022)')
res.to_csv(out / 'apexms_pilot.tsv', sep='\t', index=False)
json.dump({a: sorted(k for k, v in sets.items() if a in v) for a in sorted(universe)}, open(out / 'apexms_protein.json', 'w'))

import sys, json
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(sys.argv[3])))
import lib
ctx = lib.setup(sys.argv[1])
f = sys.argv[2]
for sheet, ratio in (('Input_DMSO_vs_GSK591', 'GSK/DMSO'), ('Input_DMSO_vs_MS023', 'MS023/DMSO'), ('polyA_DMSO_vs_GSK591', 'GSK/DMSO'), ('polyA_DMSO_vs_MS023', 'MS023/DMSO')):
    d = pd.read_excel(f, sheet_name=sheet)
    d['acc'] = [lib.resolve(ctx, str(a).split(';'), str(g)) for a, g in zip(d['Accession'], d['Gene Symbol'])]
    d = d[d.acc != ''].copy()
    d['lp'] = pd.to_numeric(d['p-value'], errors='coerce'); d['fc'] = pd.to_numeric(d[ratio], errors='coerce')
    d = d.dropna(subset=['lp', 'fc']).drop_duplicates('acc')
    sets = {'up (log2 >= 1, p <= 0.05)': set(d.acc[(d.fc >= 1) & (d.lp >= 1.301)]), 'down (log2 <= -1, p <= 0.05)': set(d.acc[(d.fc <= -1) & (d.lp >= 1.301)])}
    m = set(d.acc) & ctx['meth']
    print(f'\n== {sheet}: {len(d)} proteins quantified, {len(m)} methylated ({len(m)/len(d):.0%}); up {len(sets[list(sets)[0]])}, down {len(sets[list(sets)[1]])}')
    t = lib.tests(ctx, set(d.acc), sets, min_size=3)
    print(t[['compartment', 'members', 'methylated', 'share', 'fold', 'q', 'odds_raw', 'odds_mh', 'q_mh']].to_string(index=False, float_format=lambda v: f'{v:.3g}'))
    top = d[d.acc.isin(sets[list(sets)[0]] & m)].sort_values('fc', ascending=False)['Gene Symbol'].head(15).tolist()
    print('   methylated, up:', top)

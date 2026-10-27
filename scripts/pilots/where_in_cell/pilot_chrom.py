import sys, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
sys.path.insert(0, str(Path(__file__).parent))
import lib
repo, f, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
ctx = lib.setup(repo)
t = pd.read_excel(f, sheet_name=0, dtype={'Majority protein IDs': str})
print(t[['Training function', 'Lab internal category']].describe())
t['acc'] = [lib.resolve(ctx, str(ids).split(';')) for ids in t['Majority protein IDs']]
lost = (t['acc'] == '').sum()
t = t[t['acc'] != ''].sort_values('Interphase chromatin probability', ascending=False).drop_duplicates('acc')
icp = t.set_index('acc')['Interphase chromatin probability'].astype(float)
m = ctx['meth']; ism = icp.index.isin(list(m))
print(f'proteins {len(t) + lost}; on a reviewed entry {len(t)}; unmapped {lost}; methylated {ism.sum()} ({ism.mean():.1%})')
print(f'median ICP methylated {icp[ism].median():.3f} vs other {icp[~ism].median():.3f}; Mann-Whitney p {mannwhitneyu(icp[ism], icp[~ism]).pvalue:.2g}')
bins = {'ICP >= 0.8 (chromatin)': set(icp[icp >= 0.8].index), 'ICP 0.5-0.8': set(icp[(icp >= 0.5) & (icp < 0.8)].index),
        'ICP 0.2-0.5': set(icp[(icp >= 0.2) & (icp < 0.5)].index), 'ICP < 0.2 (not chromatin)': set(icp[icp < 0.2].index),
        'ICP >= 0.5': set(icp[icp >= 0.5].index)}
res = lib.tests(ctx, set(icp.index), bins)
lib.show(res, 'chromatin probability (Kustatscher 2014), among scored proteins')
res.to_csv(out / 'chrom_pilot.tsv', sep='\t', index=False)
json.dump({a: round(float(v), 3) for a, v in icp.items()}, open(out / 'chrom_protein.json', 'w'))

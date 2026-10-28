import sys, numpy as np, pandas as pd
from scipy.stats import fisher_exact
sys.path.insert(0, '/home/user/PTM_Pipeline/scripts')
from common import mantel_haenszel_ci
counts, rmats = sys.argv[1:3]
c = pd.read_csv(counts, sep='\t', index_col=0, low_memory=False).apply(pd.to_numeric, errors='coerce').fillna(0)
c.columns = [x.replace('_Rep', '-rep').replace('Rep', 'rep') for x in c.columns]
cpm = c / c.sum() * 1e6
inputs = [x for x in cpm.columns if 'input' in x.lower()]
expr = cpm[inputs].mean(axis=1); expressed = expr >= 1
enr = lambda b, r: np.log2((cpm[f'{b}-rep{r}'] + 0.5) / (cpm[f'{b}input-rep{r}'] + 0.5))
sp = []
for b in ('SRSF1', 'SRSF7', 'RNPS1'):
    d = [enr(b, r) - enr('NLS', r) for r in (1, 2)]
    sp.append(set(cpm.index[expressed & (d[0] >= 1) & (d[1] >= 1)]))
speckle = {g for g in set().union(*sp) if sum(g in s for s in sp) >= 2}
r = pd.read_csv(rmats, sep='\t'); r = r.loc[:, ~r.columns.duplicated()]
r = r[r.geneSymbol.isin(cpm.index[expressed])].copy()
# per event: reads supporting it (power), significance and direction
tot = lambda col: r[col].astype(str).str.split(',').apply(lambda v: sum(int(x) for x in v if x.isdigit()))
r['reads'] = tot('IJC_SAMPLE_1') + tot('SJC_SAMPLE_1') + tot('IJC_SAMPLE_2') + tot('SJC_SAMPLE_2')
r['sig'] = (r.FDR < 0.05) & (r.IncLevelDifference.abs() >= 0.1)
r['up'] = r.sig & (r.IncLevelDifference < 0)
r['speckle'] = r.geneSymbol.isin(speckle)
r['expr_decile'] = pd.qcut(np.log10(r.geneSymbol.map(expr) + 1), 10, labels=False, duplicates='drop')
r['read_decile'] = pd.qcut(np.log10(r.reads + 1), 10, labels=False, duplicates='drop')
def mh(df, outcome, strata):
    tables = []
    for _, g in df.groupby(strata):
        a = int((g.speckle & g[outcome]).sum()); b = int((~g.speckle & g[outcome]).sum()); cc = int((g.speckle & ~g[outcome]).sum()); d = int((~g.speckle & ~g[outcome]).sum())
        tables.append((a, b, cc, d))
    return mantel_haenszel_ci(tables)
print(f'events {len(r)}; significant {r.sig.sum()} (more retained {r.up.sum()}); speckle-gene events {r.speckle.sum()}')
a = int((r.speckle & r.sig).sum()); b = int((~r.speckle & r.sig).sum()); cc = int((r.speckle & ~r.sig).sum()); d = int((~r.speckle & ~r.sig).sum())
print(f'any significant change: speckle {a}/{a+cc} ({a/(a+cc):.1%}) vs other {b}/{b+d} ({b/(b+d):.1%}); OR {fisher_exact([[a,b],[cc,d]])[0]:.2f}')
for s in ('expr_decile', 'read_decile'):
    o, lo, hi, p, n = mh(r, 'sig', s); print(f'  within {s}: OR {o:.2f} ({lo:.2f}-{hi:.2f}) p {p:.2g}')
sg = r[r.sig]
a = int((sg.speckle & sg.up).sum()); b = int((~sg.speckle & sg.up).sum()); cc = int((sg.speckle & ~sg.up).sum()); d = int((~sg.speckle & ~sg.up).sum())
print(f'direction among significant events: more retained with knockdown in speckle genes {a}/{a+cc} ({a/(a+cc):.1%}) vs other genes {b}/{b+d} ({b/(b+d):.1%}); OR {fisher_exact([[a,b],[cc,d]])[0]:.2f} p {fisher_exact([[a,b],[cc,d]])[1]:.2g}')
o, lo, hi, p, n = mh(r, 'up', 'read_decile'); print(f'more retained (vs all other events) within read-depth deciles: OR {o:.2f} ({lo:.2f}-{hi:.2f}) p {p:.2g}')
print('median input CPM speckle genes', round(expr[list(speckle)].median(), 1), 'vs expressed', round(expr[expressed].median(), 1))

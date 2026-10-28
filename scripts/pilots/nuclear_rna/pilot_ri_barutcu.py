"""Pilot: do transcripts with introns PRMT5 knockdown retains sit at nuclear speckles (Barutcu et al. 2022 APEX2-seq)?
Gene level. Domain association per bait: log2 enrichment of the pulldown over its own input (CPM, both replicates), minus
the same for the nucleoplasmic NLS-APEX2 control; a gene is associated when that difference is >= 1 in both replicates.
Background: genes with a retained-intron event tested by rMATS and expressed in the APEX inputs (CPM >= 1)."""
import sys
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
counts, rmats, out = sys.argv[1:4]
c = pd.read_csv(counts, sep='\t', index_col=0, low_memory=False).apply(pd.to_numeric, errors='coerce').fillna(0)
c.columns = [x.replace('_Rep', '-rep').replace('_rep', '-rep').replace('Rep', 'rep') for x in c.columns]
cpm = c / c.sum() * 1e6
def col(name, rep):
    for cand in (f'{name}-rep{rep}',):
        if cand in cpm: return cpm[cand]
    raise KeyError(name, rep)
inputs = [x for x in cpm.columns if 'input' in x.lower()]
expressed = cpm[inputs].mean(axis=1) >= 1
def enr(bait, rep):
    return np.log2((col(bait, rep) + 0.5) / (col(bait + 'input', rep) + 0.5))
baits = {'Speckle (SRSF1)': 'SRSF1', 'Speckle (SRSF7)': 'SRSF7', 'Speckle (RNPS1)': 'RNPS1', 'Lamina (LMNA)': 'LMNA', 'Nucleolus (FBL)': 'FBL',
         'PML body (PML)': 'PML', 'PML body (SP100)': 'SP100', 'Cajal body (WRAP53)': 'WRAP53', 'Gem/Cajal (SMN2)': 'SMN2', 'Histone locus (NPAT)': 'NPAT', 'Sam68 body (SAM68)': 'SAM68'}
assoc = {}
for label, b in baits.items():
    try:
        d = [enr(b, r) - enr('NLS', r) for r in (1, 2)]
    except KeyError as e:
        print('skip', label, e); continue
    assoc[label] = set(cpm.index[expressed & (d[0] >= 1) & (d[1] >= 1)])
sp = [assoc[k] for k in assoc if k.startswith('Speckle')]
assoc['Speckle (2 of 3 markers)'] = {g for g in set().union(*sp) if sum(g in s for s in sp) >= 2}
r = pd.read_csv(rmats, sep='\t'); r = r.loc[:, ~r.columns.duplicated()]
tested = set(r.geneSymbol) & set(cpm.index[expressed])
up = set(r[(r.FDR < 0.05) & (r.IncLevelDifference <= -0.1)].geneSymbol) & tested     # more retained in SAMPLE_2 (knockdown)
down = set(r[(r.FDR < 0.05) & (r.IncLevelDifference >= 0.1)].geneSymbol) & tested
print(f'genes: expressed in APEX inputs {expressed.sum()}; with an RI event tested {len(tested)}; more retention with knockdown {len(up)}; less {len(down)}')
rows = []
for label, s in assoc.items():
    s = s & tested
    for name, hits in (('more retained', up), ('less retained', down)):
        a = len(hits & s); b = len(hits - s); cc = len(s - hits); d = len(tested - hits - s)
        orr, p = fisher_exact([[a, b], [cc, d]])
        rows.append(dict(domain=label, genes_in_domain=len(s), change=name, in_domain=a, of=len(hits), share=a / max(1, len(hits)), background_share=len(s) / len(tested), odds_ratio=orr, p=p))
t = pd.DataFrame(rows)
from statsmodels.stats.multitest import multipletests
t['q'] = multipletests(t['p'], method='fdr_bh')[1]
pd.set_option('display.width', 200)
print(t.to_string(index=False, float_format=lambda v: f'{v:.3g}'))
t.to_csv(out, sep='\t', index=False)

"""Preview data for the website's "Where in the cell" card (not a release file): per-dataset tests from the pilots, the HPA
rows of the v1.4.1 statistics, the compartment grid, the nucleus parts and each protein's placement."""
import sys, json, math
from pathlib import Path
import pandas as pd
repo, sc, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
l2 = lambda v: None if v is None or (isinstance(v, float) and math.isnan(v)) or v <= 0 else (99 if math.isinf(v) else round(math.log2(v), 3))
sig = lambda v: None if v is None or (isinstance(v, float) and math.isnan(v)) else float(f'{v:.3g}')

def rows_of(t):
    return [[r.compartment, int(r.members), int(r.methylated), l2(r.fold), sig(r.q), l2(r.odds_mh), l2(r.lo_mh), l2(r.hi_mh), sig(r.q_mh)]
            for r in t.itertuples()]
hcm = pd.read_csv(sc / 'hcm_pilot.tsv', sep='\t')
tests = {'hcm': rows_of(hcm[hcm.background == 'localised (NMF)']), 'hcm_safe': rows_of(hcm[hcm.background == 'localised (SAFE)']),
         'oip': rows_of(pd.read_csv(sc / 'oip_pilot.tsv', sep='\t')),
         'chep': rows_of(pd.read_csv(sc / 'chrom_pilot.tsv', sep='\t').query("compartment != 'ICP >= 0.5'"))}
st = pd.read_csv(repo / 'releases/v1.4.1/statistics_v1.4.1.tsv', sep='\t', dtype=str, keep_default_na=False)
num = lambda x: float(x) if x not in ('', None) else None
en = st[(st.analysis == 'hpa_location') & (st.method == 'enricher')]
mh = st[(st.analysis == 'hpa_location') & (st.method == 'mantel_haenszel')].set_index('category')
tests['hpa'] = sorted([[r.category, int(r.a) + int(r.c), int(r.a), l2(num(r.fold_enrichment)), sig(num(r.q)),
                        l2(num(mh.at[r.category, 'odds_ratio'])), l2(num(mh.at[r.category, 'odds_ratio_ci_low'])), l2(num(mh.at[r.category, 'odds_ratio_ci_high'])), sig(num(mh.at[r.category, 'q']))]
                       for r in en.itertuples() if int(r.a) + int(r.c) >= 20], key=lambda x: -(x[3] or -9))
hcm_n = hcm[hcm.background == 'localised (NMF)']
datasets = [
    dict(key='hcm', name='Human Cell Map', method='Proximity labelling', how='BioID with 192 markers; each protein in its NMF compartment', cells='HEK293',
         ref='Go et al. 2021', doi='10.1038/s41586-021-03592-2', licence='CC BY 4.0', url='https://humancellmap.org'),
    dict(key='hcm_safe', name='Human Cell Map, finer domains', method='Proximity labelling', how='the same map, SAFE domains (proteins without a domain left out)', cells='HEK293',
         ref='Go et al. 2021', doi='10.1038/s41586-021-03592-2', licence='CC BY 4.0', url='https://humancellmap.org'),
    dict(key='oip', name='Organelle IP', method='Organelle pull-down', how='organelle immunoprecipitation; each protein in its cluster (unclassified left out)', cells='HEK293T',
         ref='Hein et al. 2025', doi='10.1016/j.cell.2024.11.028', licence='BSD-3 (analysis repository)', url='https://github.com/czbiohub-sf/Organelle_IP_analyses_and_figures'),
    dict(key='chep', name='Chromatin probability', method='Fractionation', how='chromatin enrichment (ChEP) in many conditions, as a probability per protein', cells='human cell lines',
         ref='Kustatscher et al. 2014', doi='10.1002/embj.201387614', licence='CC BY', url='https://doi.org/10.1002/embj.201387614'),
    dict(key='hpa', name='Human Protein Atlas', method='Imaging', how='immunofluorescence; release statistics (hpa_location), locations of 20 or more proteins', cells='many cell lines',
         ref='HPA 25.1', doi='10.1126/science.aal3321', licence='CC BY 4.0', url='https://www.proteinatlas.org/humanproteome/subcellular'),
]
for d in datasets:
    t = tests[d['key']]
    if d['key'] == 'hpa':
        r0 = en.iloc[0]; d['universe'] = int(r0.a) + int(r0.b) + int(r0.c) + int(r0.d); d['methylated'] = int(r0.a) + int(r0.b)
hp, op, cp = (json.load(open(sc / f)) for f in ('hcm_protein.json', 'oip_protein.json', 'chrom_protein.json'))
prot = pd.read_csv(repo / 'releases/v1.4.1/human_methylarginine_proteins_v1.4.1.tsv', sep='\t', dtype=str, keep_default_na=False)
meth = {a.split('-')[0] for a in prot.uniprot_accession}
for d in datasets:
    keys = {'hcm': hp, 'hcm_safe': {a: v for a, v in hp.items() if v[1] not in ('no prediction', '-')}, 'oip': {a: v for a, v in op.items() if v != 'unclassified'}, 'chep': cp}.get(d['key'])
    if keys is not None:
        d['universe'] = len(keys); d['methylated'] = len(set(keys) & meth)
cats = {'hcm': sorted({v[0] for v in hp.values()}), 'hcm_safe': sorted({v[1] for v in hp.values()}), 'oip': sorted(set(op.values()))}
proteins = {}
for a in sorted(set(hp) | set(op) | set(cp)):
    h = hp.get(a)
    proteins[a] = [cats['hcm'].index(h[0]) if h else -1, cats['hcm_safe'].index(h[1]) if h else -1,
                   cats['oip'].index(op[a]) if a in op else -1, round(cp[a] * 100) if a in cp else -1]
# the compartment grid: [group, [[row, {dataset: its category}, nucleus part or '']]]
G = lambda **k: k
groups = [
    ['Nucleus', [
        ['Chromatin', G(hcm='chromatin', chep='ICP >= 0.8 (chromatin)'), 'chromatin'],
        ['Nucleoplasm', G(hcm='nucleoplasm', oip='nucleus', hpa='Nucleoplasm'), 'nucleoplasm'],
        ['Nucleolus', G(hcm='nucleolus', oip='nucleolus', hpa='Nucleoli'), 'nucleolus'],
        ['Nuclear speckles', G(hcm_safe='splicing speckles, nucleolus', hpa='Nuclear speckles'), 'speckles'],
        ['Paraspeckles', G(hcm_safe='paraspeckles'), 'paraspeckles'],
        ['Nuclear bodies (PML, Cajal …)', G(hcm='nuclear body', hpa='Nuclear bodies'), 'bodies'],
        ['Nuclear envelope and pores', G(hcm='nuclear outer membrane-ER membrane network', hcm_safe='nuclear pore', hpa='Nuclear membrane'), 'envelope'],
        ['Mitotic chromosome', G(hpa='Mitotic chromosome'), '']]],
    ['Cytoplasm', [
        ['RNA granules (stress granules, P-bodies)', G(hcm='cytoplasmic ribonucleoprotein granule', hcm_safe='cytoplasmic ribonucleoprotein granule', oip='stress_granule'), 'granules'],
        ['P-bodies', G(oip='p-body', hcm_safe='centrosome enriched, P-body enriched'), 'pbodies'],
        ['Translation', G(oip='translation'), 'ribosomes'],
        ['Cytosol', G(hcm='actin cytoskeleton, cytosol', oip='cytosol', hpa='Cytosol'), 'cytosol'],
        ['Actin cytoskeleton', G(oip='actin_cytoskeleton', hcm_safe='actin cytoskeleton', hpa='Actin filaments'), 'actin'],
        ['Microtubules', G(hcm='microtubule cytoskeleton', hpa='Microtubules'), 'microtubules'],
        ['Centrosome', G(hcm='centrosome', oip='centrosome', hpa='Centrosome'), 'centrosome'],
        ['Proteasome', G(oip='proteasome'), '']]],
    ['Organelles and membranes', [
        ['Mitochondria', G(hcm='mitochondrial matrix', oip='mitochondrion', hpa='Mitochondria'), 'mito'],
        ['Endoplasmic reticulum', G(hcm='ER membrane', oip='ER', hpa='Endoplasmic reticulum'), 'er'],
        ['Golgi', G(hcm='Golgi apparatus', oip='Golgi', hpa='Golgi apparatus'), 'golgi'],
        ['Endosomes and lysosomes', G(hcm='endosome, lysosome', oip='lysosome', hpa='Vesicles'), 'endolyso'],
        ['Peroxisomes', G(hcm='mitochondrial outer membrane, peroxisome', oip='peroxisome'), 'perox'],
        ['Plasma membrane', G(hcm='plasma membrane', oip='plasma_membrane', hpa='Plasma membrane'), 'pm']]],
]
have = {k: {r[0] for r in v} for k, v in tests.items()}
for g in groups:
    for row in g[1]:
        for k, c in list(row[1].items()):
            assert c in have[k], (k, c)
json.dump(dict(preview=True, generated='pilot', datasets=datasets, tests=tests, test_fields=['category', 'members', 'methylated', 'log2 fold enrichment', 'q', 'log2 OR within abundance deciles', 'CI low', 'CI high', 'q'],
               groups=groups, categories=cats, proteins=proteins, protein_fields=['hcm NMF index', 'hcm SAFE index', 'oip index', 'chromatin probability x100']),
          open(out, 'w'), separators=(',', ':'))
print({d['key']: (d['universe'], d['methylated']) for d in datasets}, len(proteins), out.stat().st_size)

"""Preview data for the website's "Where in the cell" card (not a release file): per-dataset tests from the pilots, the HPA
rows of the v1.4.1 statistics, the compartment grid, the parts of the drawings and each protein's placement, with its
PaxDb decile in each map so the website can repeat the tests within abundance deciles for a selection."""
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
         'chep': rows_of(pd.read_csv(sc / 'chrom_pilot.tsv', sep='\t').query("compartment != 'ICP >= 0.5'")),
         'apexms': rows_of(pd.read_csv(sc / 'apexms_pilot.tsv', sep='\t'))}
st = pd.read_csv(repo / 'releases/v1.4.1/statistics_v1.4.1.tsv', sep='\t', dtype=str, keep_default_na=False)
num = lambda x: float(x) if x not in ('', None) else None
en = st[(st.analysis == 'hpa_location') & (st.method == 'enricher')]
mh = st[(st.analysis == 'hpa_location') & (st.method == 'mantel_haenszel')].set_index('category')
tests['hpa'] = sorted([[r.category, int(r.a) + int(r.c), int(r.a), l2(num(r.fold_enrichment)), sig(num(r.q)),
                        l2(num(mh.at[r.category, 'odds_ratio'])), l2(num(mh.at[r.category, 'odds_ratio_ci_low'])), l2(num(mh.at[r.category, 'odds_ratio_ci_high'])), sig(num(mh.at[r.category, 'q']))]
                       for r in en.itertuples() if int(r.a) + int(r.c) >= 20], key=lambda x: -(x[3] or -9))
hcm_n = hcm[hcm.background == 'localised (NMF)']
datasets = [
    dict(key='hcm', name='Human Cell Map', method='Proximity labelling', how='proximity labelling from 192 baits; each protein sits in the compartment whose baits label it most', cells='HEK293',
         ref='Go et al. 2021', doi='10.1038/s41586-021-03592-2', licence='CC BY 4.0', url='https://humancellmap.org'),
    dict(key='hcm_safe', name='Human Cell Map, finer domains', method='Proximity labelling', how='the same map at a finer scale: neighbourhoods of proteins that label one another (SAFE domains)', cells='HEK293',
         ref='Go et al. 2021', doi='10.1038/s41586-021-03592-2', licence='CC BY 4.0', url='https://humancellmap.org'),
    dict(key='oip', name='Organelle IP', method='Organelle pull-down', how='organelles pulled out of the cell one at a time; each protein belongs to the organelle it comes down with', cells='HEK293T',
         ref='Hein et al. 2025', doi='10.1016/j.cell.2024.11.028', licence='BSD-3 (analysis repository)', url='https://github.com/czbiohub-sf/Organelle_IP_analyses_and_figures'),
    dict(key='chep', name='Chromatin probability', method='Fractionation', how='how reliably a protein comes down with chromatin across many conditions, as a probability', cells='human cell lines',
         ref='Kustatscher et al. 2014', doi='10.1002/embj.201387614', licence='CC BY', url='https://doi.org/10.1002/embj.201387614'),
    dict(key='apexms', name='Speckle and lamina proteomes', method='Proximity labelling', how='APEX2 on speckle markers (SRSF1, SRSF7, RNPS1), lamins (LMNA, LMNB1) and a PML-body marker (SP100) labels the proteins around them; a protein belongs where SAINT gives it a probability of at least 0.9 (speckles: with two of the three markers)', cells='HEK293',
         ref='Barutcu et al. 2022', doi='10.1016/j.molcel.2021.12.010', licence='Mol Cell Table S1 (calls derived from it); raw data MassIVE MSV000087649', url='https://doi.org/10.1016/j.molcel.2021.12.010'),
    dict(key='hpa', name='Human Protein Atlas', method='Imaging', how='antibody staining in many cell lines; the locations seen with approved or better reliability (those holding 20 or more proteins)', cells='many cell lines',
         ref='HPA 25.1', doi='10.1126/science.aal3321', licence='CC BY 4.0', url='https://www.proteinatlas.org/humanproteome/subcellular'),
]
for d in datasets:
    t = tests[d['key']]
    if d['key'] == 'hpa':
        r0 = en.iloc[0]; d['universe'] = int(r0.a) + int(r0.b) + int(r0.c) + int(r0.d); d['methylated'] = int(r0.a) + int(r0.b)
hp, op, cp, ap = (json.load(open(sc / f)) for f in ('hcm_protein.json', 'oip_protein.json', 'chrom_protein.json', 'apexms_protein.json'))
prot = pd.read_csv(repo / 'releases/v1.4.1/human_methylarginine_proteins_v1.4.1.tsv', sep='\t', dtype=str, keep_default_na=False)
meth = {a.split('-')[0] for a in prot.uniprot_accession}
universes = {'hcm': set(hp), 'hcm_safe': {a for a, v in hp.items() if v[1] not in ('no prediction', '-')}, 'oip': {a for a, v in op.items() if v != 'unclassified'}, 'chep': set(cp), 'apexms': set(ap)}
for d in datasets:
    if d['key'] in universes:
        d['universe'] = len(universes[d['key']]); d['methylated'] = len(universes[d['key']] & meth)
# The release's HPA test: every protein with a location scored Approved or better (release_tables.hpa_rows)
hpa_all = pd.read_csv(repo / 'inputs/layers/hpa_subcellular_locations.tsv', sep='\t', dtype=str, keep_default_na=False)
confident = hpa_all[hpa_all['reliability'] != 'Uncertain']
hpa_sets = confident.groupby('canonical_UniProtAC')['location'].agg(lambda s: sorted(set(s))).to_dict()
hpa_ds = next(d for d in datasets if d['key'] == 'hpa')
assert len(hpa_sets) == hpa_ds['universe'], (len(hpa_sets), hpa_ds['universe'])
universes['hpa'] = set(hpa_sets)
cats = {'hcm': sorted({v[0] for v in hp.values()}), 'hcm_safe': sorted({v[1] for v in hp.values()}), 'oip': sorted(set(op.values())),
        'hpa': sorted(set(confident['location'])), 'apexms': sorted({c for v in ap.values() for c in v})}
# Each protein's PaxDb decile among the proteins each map places, as the tests stratify (abundance_strata.strata):
# one character per map of MAPS, 0-9 (the least to the most abundant), x without a PaxDb value, - not placed.
sys.path.insert(0, str(Path(__file__).parent))
import lib
ctx = lib.setup(repo)
MAPS = ['hcm', 'hcm_safe', 'oip', 'chep', 'hpa', 'apexms']
dec = {k: ctx['A'].strata(sorted(universes[k]), ctx['pax']).to_dict() for k in MAPS}
code = lambda k, a: '-' if a not in universes[k] else 'x' if dec[k][a] == 'none' else str(int(dec[k][a]) - 1)
# chromatin probability x1000, kept on its side of each bin edge the pilot uses (0.2, 0.5, 0.8)
def icp_code(v):
    c = round(v * 1000)
    for edge in (200, 500, 800):
        if v < edge / 1000 and c >= edge: c = edge - 1
    return c
proteins = {}
for a in sorted(set(hp) | set(op) | set(cp) | set(hpa_sets) | set(ap)):
    h = hp.get(a)
    proteins[a] = [cats['hcm'].index(h[0]) if h else -1, cats['hcm_safe'].index(h[1]) if h else -1,
                   cats['oip'].index(op[a]) if a in op else -1, icp_code(cp[a]) if a in cp else -1,
                   [cats['hpa'].index(x) for x in hpa_sets.get(a, [])], ''.join(code(k, a) for k in MAPS), [cats['apexms'].index(x) for x in ap.get(a, [])]]
# the compartment grid: [group, [[row, {dataset: its category}, nucleus part or '']]]
G = lambda **k: k
groups = [
    ['Nucleus', [
        ['Chromatin', G(hcm='chromatin', chep='ICP >= 0.8 (chromatin)'), 'chromatin'],
        ['Nucleoplasm', G(hcm='nucleoplasm', oip='nucleus', hpa='Nucleoplasm'), 'nucleoplasm'],
        ['Nucleolus', G(hcm='nucleolus', oip='nucleolus', hpa='Nucleoli'), 'nucleolus'],
        ['Nuclear speckles', G(hcm_safe='splicing speckles, nucleolus', hpa='Nuclear speckles', apexms='Nuclear speckles'), 'speckles'],
        ['Paraspeckles', G(hcm_safe='paraspeckles'), 'paraspeckles'],
        ['Nuclear bodies (PML, Cajal …)', G(hcm='nuclear body', hpa='Nuclear bodies'), 'bodies'],
        ['PML bodies', G(apexms='PML bodies'), 'pml'],
        ['Nuclear envelope and lamina', G(hcm='nuclear outer membrane-ER membrane network', hpa='Nuclear membrane', apexms='Nuclear lamina'), 'envelope'],
        ['Nuclear pores', G(hcm_safe='nuclear pore'), 'pores'],
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
# ---- RNA and DNA maps (pilot_layers.py): where the RNAs each protein with eCLIP data binds gather (APEX-seq,
# nuclear-body APEX2-seq) and where their genes lie (TSA-seq, Hi-C subcompartments), and the methylated proteins
# against the others in each compartment. Written beside the file as subcellular_na.json: the card states the
# comparison in a sentence, and each protein's page lists its own calls.
lt = pd.read_csv(sc / 'layer_tests.tsv', sep='\t')
lp = json.load(open(sc / 'layer_proteins.json'))
geo = lambda g: f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={g}'  # noqa: E731
na_ds = {
    'apex': dict(key='apex', name='APEX-seq', method='RNA proximity labelling', cells='HEK293T', ref='Fazal et al. 2019', doi='10.1016/j.cell.2019.05.027',
                 how='APEX2 placed in eight compartments labels the RNAs around it; an RNA belongs to a compartment when it is labelled more than in the control (FDR < 0.05)',
                 licence='GEO GSE116008, no restrictions', url=geo('GSE116008'), expression='the mean of the control libraries'),
    'nb': dict(key='nb', name='Nuclear-body APEX2-seq', method='RNA proximity labelling', cells='HEK293', ref='Barutcu et al. 2022', doi='10.1016/j.molcel.2021.12.010',
               how='APEX2 on markers of speckles, the lamina and nuclear bodies; an RNA belongs to a body when it is labelled at least twice as much as by nucleoplasmic APEX2 in both replicates',
               licence='GEO GSE176439, no restrictions', url=geo('GSE176439'), expression='the mean of the input libraries'),
    'tsa': dict(key='tsa', name='TSA-seq', method='DNA distance mapping', cells='K562', ref='Chen et al. 2018', doi='10.1083/jcb.201807108',
                how='tyramide signal spreading from speckles (SON) and from the lamina (lamin B) marks DNA by its distance; a gene lies near when its signal is in the top tenth of the genome',
                licence='GEO GSE81553, no restrictions', url=geo('GSE81553'), expression='K562 total RNA-seq (ENCODE ENCSR109IQO)'),
    'sub': dict(key='sub', name='Hi-C subcompartments', method='Chromosome contacts (Hi-C)', cells='GM12878', ref='Rao et al. 2014', doi='10.1016/j.cell.2014.11.021',
                how='the contacts each stretch of the genome makes divide it into six subcompartments: A1 and A2 active (A1 beside speckles), B1 Polycomb-repressed, '
                    'B2 at the nucleolus and the lamina with pericentromeric heterochromatin, B3 at the lamina, B4 a few clusters on chromosome 19; a gene takes the subcompartment at its middle',
                licence='GEO GSE63525, no restrictions', url=geo('GSE63525'), expression='GM12878 total RNA-seq (ENCODE ENCSR000AEE)'),
}
NA_MAPS = ['apex', 'nb', 'tsa', 'sub']
na_cats = {k: lp['categories'][k] for k in NA_MAPS}
na_tested = {k: {a for a, x in lp['proteins'].items() if k in x['tested']} for k in NA_MAPS}
def na_rows(k):
    t = lt[lt['map'] == k]
    m2 = lambda v: None if pd.isna(v) else 99 if v == math.inf else -99 if v == -math.inf else round(float(v), 3)  # noqa: E731
    return [r + [m2(x.log2_or_median_methylated), m2(x.log2_or_median_other), sig(x.q_size)]
            for r, x in zip(rows_of(t), t.itertuples())]
na_tests = {k: na_rows(k) for k in NA_MAPS}
for k in NA_MAPS:
    na_ds[k]['universe'] = len(na_tested[k]); na_ds[k]['methylated'] = len(na_tested[k] & meth); na_ds[k]['placed'] = lp['placed'][k]
# each protein's calls in each map, null where the map tests it not (fewer than 20 of its bound RNAs placed)
na_prot = {a: [x['name'], [[[na_cats[k].index(c), v, q] for c, v, q in x['calls'].get(k, []) if c in na_cats[k]] if k in x['tested'] else None for k in NA_MAPS]]
           for a, x in sorted(lp['proteins'].items()) if x['tested']}
json.dump(dict(preview=True, maps=NA_MAPS, datasets=na_ds, tests=na_tests, categories=na_cats, proteins=na_prot,
               test_fields=['compartment', "eCLIP'd proteins whose bound RNAs gather there", 'methylated', 'log2 fold enrichment', 'q', 'log2 OR within PaxDb deciles', 'CI low', 'CI high', 'q',
                            "median log2 OR of the methylated proteins' bound RNAs", 'of the other proteins', 'q (Mann-Whitney, BH per map)'],
               protein_fields=['name', 'per map of maps: [category index, log2 OR of its bound RNAs within expression and gene-length deciles, q] where q <= 0.05; null where the map does not test it']),
          open(out.parent / 'subcellular_na.json', 'w'), separators=(',', ':'))
print('RNA/DNA maps', {k: (na_ds[k]['universe'], na_ds[k]['methylated']) for k in NA_MAPS}, len(na_prot), (out.parent / 'subcellular_na.json').stat().st_size)

json.dump(dict(preview=True, generated='pilot', datasets=datasets, tests=tests, test_fields=['category', 'members', 'methylated', 'log2 fold enrichment', 'q', 'log2 OR within abundance deciles', 'CI low', 'CI high', 'q'],
               groups=groups, categories=cats, proteins=proteins, decile_maps=MAPS,
               protein_fields=['hcm NMF index', 'hcm SAFE index', 'oip index', 'chromatin probability x1000', 'hpa location indices (Approved or better)',
                               'PaxDb decile among the proteins each map of decile_maps places (0-9; x no value; - not placed)', 'apexms category indices']),
          open(out, 'w'), separators=(',', ':'))
print({d['key']: (d['universe'], d['methylated']) for d in datasets}, len(proteins), out.stat().st_size)

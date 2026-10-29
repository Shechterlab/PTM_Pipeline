"""The cell-map layer (scripts/release_layers.py CELL): where methylated proteins sit in the cell, and where the RNAs that
RNA-binding proteins bind gather.

Protein maps (staged by scripts/52_stage_cell_maps.py into inputs/layers/): the Human Cell Map's compartments and its
finer SAFE domains, organelle IP, the speckle, lamina and PML-body proteomes of Barutcu et al. 2022, and the interphase
chromatin probability of Kustatscher et al. 2014 (in four bins). The Human Protein Atlas, the fifth, is tested where it
always was (analysis hpa_location, scripts/release_tables.py). Each map is tested against the proteins it places: the
methylated proteins among the members of each compartment, with clusterProfiler enricher (one-sided hypergeometric, BH
within the map), and again within PaxDb abundance deciles (scripts/abundance_strata.py, TESTS).

RNA and DNA maps: APEX-seq, nuclear-body APEX2-seq, TSA-seq and Hi-C subcompartments place RNAs (genes). Each protein with
eCLIP data (ENCODE and ENCORE, inputs/nucleic_acid/: the genes with most reproducible peaks or windows in each
experiment, pooled over the protein's experiments) is a unit. Its bound RNAs gather in a compartment when a map places
them there more often than the other RNAs it places, within deciles of expression in the map's cells crossed with
deciles of gene length (GENCODE v19, whose missing genes are one more length stratum; Mantel-Haenszel, BH over every
protein and compartment of the map; q <= 0.05 and odds ratio above 1). Gene length matters: eCLIP target lists favour long genes, whose introns collect more peaks, and
compartments differ in it (the RNAs nuclear-body APEX2-seq places at the lamina come from genes of median 319 kb, its
speckle RNAs from genes of 12 kb). A map tests the proteins with at least 20 bound RNAs it places, in the compartments
where it places at least 20 RNAs. The methylated proteins are then compared with the other eCLIP'd proteins a map tests,
as the protein maps are (enricher, and within PaxDb deciles), and, because a call depends on how many bound RNAs a map
places, the size of every protein's odds ratio, called or not, is compared too (two-sided Mann-Whitney U, BH per map).
"""
from __future__ import annotations

import functools
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

MIN_SET = 5      # proteins a compartment needs to be tested (enricher's smallest set)
MIN_BOUND = 20   # bound RNAs a map must place for a protein to be tested in it, and RNAs a compartment must hold
Q = 0.05
CHEP_BINS = [('ICP >= 0.8 (chromatin)', 0.8, math.inf), ('ICP 0.5-0.8', 0.5, 0.8), ('ICP 0.2-0.5', 0.2, 0.5), ('ICP < 0.2 (not chromatin)', -math.inf, 0.2)]

# analysis -> the website's key, its name, the staged table, and the background its rows name
PROTEIN_MAPS = {
    'cell_hcm_compartment': ('hcm', 'Human Cell Map compartments', 'cell_map_hcm_go2021.tsv', 'the other proteins the Human Cell Map localises (HEK293 BioID, NMF compartments)'),
    'cell_hcm_domain': ('hcm_safe', 'Human Cell Map domains', 'cell_map_hcm_go2021.tsv', 'the other proteins the Human Cell Map places in a SAFE domain'),
    'cell_organelle_ip': ('oip', 'Organelle IP', 'cell_map_organelle_ip_hein2025.tsv', 'the other proteins organelle IP classifies (HEK293T)'),
    'cell_nuclear_body_proteome': ('apexms', 'Speckle, lamina and PML-body proteomes', 'cell_map_nuclear_bodies_barutcu2022.tsv',
                                   'the other proteins the speckle, lamina and PML-body APEX2-MS experiments detected (Flp-In HEK293)'),
    'cell_chromatin': ('chep', 'Chromatin probability', 'cell_map_chromatin_kustatscher2014.tsv', 'the other proteins with an interphase chromatin probability'),
}
GENE_MAPS = {
    'eclip_rna_apexseq': ('apex', 'APEX-seq', 'rna_map_apexseq_fazal2019.tsv', 'expression_cpm'),
    'eclip_rna_nuclear_body': ('nb', 'Nuclear-body APEX2-seq', 'rna_map_nuclear_bodies_barutcu2022.tsv', 'expression_cpm'),
    'eclip_gene_tsaseq': ('tsa', 'TSA-seq', 'dna_map_tsaseq_chen2018.tsv', 'expression_tpm'),
    'eclip_gene_hic': ('sub', 'Hi-C subcompartments', 'dna_map_hic_rao2014.tsv', 'expression_tpm'),
}
FILES = sorted({m[2] for m in PROTEIN_MAPS.values()} | {m[2] for m in GENE_MAPS.values()}) + ['cell_maps_summary.json']


def read(inputs, name: str) -> pd.DataFrame:
    return pd.read_csv(Path(inputs) / 'layers' / name, sep='\t', dtype=str, keep_default_na=False)


def entry(acc: str) -> str:
    return str(acc).split('-')[0]


def groups(frame: pd.DataFrame, column: str) -> dict[str, set]:
    return {c: set(g['canonical_UniProtAC']) for c, g in frame.groupby(column)}


@functools.lru_cache(maxsize=4)
def protein_maps(inputs: str) -> dict[str, tuple[set, dict[str, set]]]:
    """analysis -> (the proteins the map places, {compartment: its members})."""
    out = {}
    h = read(inputs, 'cell_map_hcm_go2021.tsv')
    out['cell_hcm_compartment'] = (set(h['canonical_UniProtAC']), groups(h, 'nmf_compartment'))
    hs = h[h['safe_domain'] != '']
    out['cell_hcm_domain'] = (set(hs['canonical_UniProtAC']), groups(hs, 'safe_domain'))
    o = read(inputs, 'cell_map_organelle_ip_hein2025.tsv')
    o = o[o['compartment'] != 'unclassified']
    out['cell_organelle_ip'] = (set(o['canonical_UniProtAC']), groups(o, 'compartment'))
    b = read(inputs, 'cell_map_nuclear_bodies_barutcu2022.tsv')
    sets = defaultdict(set)
    for a, cs in zip(b['canonical_UniProtAC'], b['compartments']):
        for c in cs.split('; '):
            if c and c != 'none':
                sets[c].add(a)
    out['cell_nuclear_body_proteome'] = (set(b['canonical_UniProtAC']), dict(sets))
    c = read(inputs, 'cell_map_chromatin_kustatscher2014.tsv')
    icp = c['interphase_chromatin_probability'].astype(float)
    out['cell_chromatin'] = (set(c['canonical_UniProtAC']), {name: set(c.loc[(icp >= lo) & (icp < hi), 'canonical_UniProtAC']) for name, lo, hi in CHEP_BINS})
    return out


@functools.lru_cache(maxsize=4)
def eclip_targets(inputs: str) -> tuple[dict[str, set], dict[str, str]]:
    """The genes each protein's eCLIP experiments bind (Ensembl IDs without version), and the protein's name there."""
    nuc = Path(inputs) / 'nucleic_acid'
    targets, name = defaultdict(set), {}
    for src, key, sym in (('encode', 'encode_experiment', 'target'), ('encore', 'encore_experiment', 'rbp')):
        exp = pd.read_csv(nuc / f'{src}_eclip_experiments.tsv', sep='\t', dtype=str, keep_default_na=False)
        for a, n in zip(exp['uniprot_accession'], exp[sym]):
            if a:
                name.setdefault(entry(a), n)
        acc_of = exp.set_index(key)['uniprot_accession']
        tg = pd.read_csv(nuc / f'{src}_eclip_target_genes.tsv', sep='\t', dtype=str, keep_default_na=False)
        for e, g in zip(tg[key], tg['gene_id']):
            a = entry(acc_of.get(e, ''))
            if a:
                targets[a].add(str(g).split('.')[0])
    return dict(targets), name


@functools.lru_cache(maxsize=4)
def bound_rna(inputs: str) -> dict[str, dict]:
    """analysis -> {tested: proteins with >= MIN_BOUND bound RNAs placed, rows: one per protein and tested compartment
    (accession, compartment, placed, in_compartment, odds_ratio, ci_low, ci_high, p, q), calls: {compartment: proteins
    whose bound RNAs gather there}, categories: the compartments tested, sizes: RNAs placed in each, placed: RNAs placed}."""
    import abundance_strata as A
    from common import benjamini_hochberg, mantel_haenszel_ci
    targets, _ = eclip_targets(inputs)
    out = {}
    for analysis, (_, _, name, ecol) in GENE_MAPS.items():
        t = read(inputs, name)
        universe = list(t['gene_id'])
        if universe != sorted(universe):
            raise SystemExit(f'cell_maps: {name} is not sorted by gene_id')
        pos = {g: i for i, g in enumerate(universe)}
        expr = pd.Series(t[ecol].astype(float).to_numpy(), index=universe)
        glen = pd.Series(pd.to_numeric(t['gene_length'], errors='coerce').to_numpy(), index=universe)
        by_expr = A.strata(universe, expr).to_numpy(dtype=str)
        by_len = A.strata(universe, glen).to_numpy(dtype=str)
        _, layer = np.unique(np.char.add(np.char.add(by_expr, '|'), by_len), return_inverse=True)
        k = int(layer.max()) + 1
        n_all = np.bincount(layer, minlength=k)
        cats_of = [[c for c in cs.split('; ') if c] for cs in t['compartments']]
        in_cat = {c: np.array([c in cs for cs in cats_of]) for c in sorted({c for cs in cats_of for c in cs})}
        sizes = {c: int(v.sum()) for c, v in in_cat.items()}
        cats = [c for c, v in in_cat.items() if v.sum() >= MIN_BOUND]
        rows, tested = [], set()
        for acc in sorted(targets):
            idx = np.array(sorted(pos[g] for g in targets[acc] if g in pos), dtype=int)
            if len(idx) < MIN_BOUND:
                continue
            tested.add(acc)
            sel = np.zeros(len(universe), dtype=bool)
            sel[idx] = True
            n_sel = np.bincount(layer[sel], minlength=k)
            for c in cats:
                n_cat = np.bincount(layer[in_cat[c]], minlength=k)
                both = np.bincount(layer[sel & in_cat[c]], minlength=k)
                tables = [(int(both[i]), int(n_sel[i] - both[i]), int(n_cat[i] - both[i]), int(n_all[i] - n_sel[i] - n_cat[i] + both[i])) for i in range(k)]
                odds, lo, hi, p, _ = mantel_haenszel_ci(tables)
                rows.append({'accession': acc, 'compartment': c, 'placed': len(idx), 'in_compartment': int(both.sum()),
                             'odds_ratio': odds, 'ci_low': lo, 'ci_high': hi, 'p': p})
        qs = benjamini_hochberg([r['p'] if r['p'] is not None else 1.0 for r in rows])
        calls = defaultdict(set)
        for r, q in zip(rows, qs):
            r['q'] = float(q)
            if q <= Q and r['odds_ratio'] is not None and r['odds_ratio'] > 1:
                calls[r['compartment']].add(r['accession'])
        out[analysis] = {'tested': tested, 'rows': rows, 'calls': dict(calls), 'categories': cats, 'sizes': sizes, 'placed': len(universe)}
    return out


def log2_or(odds) -> float:
    return math.log2(odds) if odds > 0 else -math.inf   # inf stays inf


def size_tests(b: dict, target: set) -> list[dict]:
    """Per compartment with at least MIN_SET calls: every tested protein's log2 odds ratio, methylated against the
    others (two-sided Mann-Whitney U, BH over the map's compartments)."""
    from scipy.stats import mannwhitneyu
    from common import benjamini_hochberg
    size = defaultdict(dict)
    for r in b['rows']:
        if r['odds_ratio'] is not None:
            size[r['compartment']][r['accession']] = log2_or(r['odds_ratio'])
    out = []
    for c in sorted(x for x, v in b['calls'].items() if len(v) >= MIN_SET):
        yes = [v for a, v in size[c].items() if a in target]
        no = [v for a, v in size[c].items() if a not in target]
        p = float(mannwhitneyu(yes, no).pvalue) if yes and no else math.nan
        med = lambda v: float(np.median(v)) if v else math.nan  # noqa: E731
        out.append({'compartment': c, 'n_methylated': len(yes), 'n_other': len(no), 'median_methylated': med(yes), 'median_other': med(no), 'p': p})
    qs = benjamini_hochberg([1.0 if math.isnan(r['p']) else r['p'] for r in out]) if out else []
    for r, q in zip(out, qs):
        r['q'] = float(q)
    return out


def target_of(sites: pd.DataFrame) -> set:
    return {entry(a) for a in sites.loc[sites['counted'], 'uniprot_accession']}


def stats(sites: pd.DataFrame, inputs, sequences: dict, workdir: Path) -> list[dict]:
    from enrichment import Job, run
    target = target_of(sites)
    pm, br = protein_maps(str(inputs)), bound_rna(str(inputs))
    jobs, background = [], {}
    term2gene = lambda sets: pd.DataFrame([(c, a) for c in sorted(sets) for a in sorted(sets[c])], columns=['term', 'gene'])  # noqa: E731
    for analysis, (_, _, _, bg) in PROTEIN_MAPS.items():
        universe, sets = pm[analysis]
        jobs.append(Job(analysis, target & universe, universe, sets=term2gene(sets), whole_universe=True, min_size=MIN_SET))
        background[analysis] = bg
    for analysis, (_, label, _, _) in GENE_MAPS.items():
        b = br[analysis]
        jobs.append(Job(analysis, target & b['tested'], b['tested'], sets=term2gene(b['calls']), whole_universe=True, min_size=MIN_SET))
        background[analysis] = (f'the other proteins with eCLIP data that {label} tests (at least {MIN_BOUND} bound RNAs placed). A set is the '
                                f'proteins whose bound RNAs gather in the compartment')
    res = run(jobs, Path(workdir) / 'enricher')
    rows = []
    for job in jobs:
        for r in sorted(res[job.job].to_dict('records'), key=lambda r: (r['pvalue'], r['ID'])):
            rows.append({'analysis': job.job, 'category': r['ID'], 'stratum': 'all', 'a': r['a'], 'b': r['b'], 'c': r['c'], 'd': r['d'],
                         'fold_enrichment': r['FoldEnrichment'], 'log2_fold_enrichment': math.log2(r['FoldEnrichment']) if r['FoldEnrichment'] > 0 else None,
                         'p': r['pvalue'], 'q': r['p.adjust'], 'method': 'enricher', 'background': background[job.job], 'n_strata': ''})
    for analysis, (_, label, _, _) in GENE_MAPS.items():
        for r in size_tests(br[analysis], target):
            rows.append({'analysis': analysis, 'category': r['compartment'], 'stratum': 'all', 'a': r['n_methylated'], 'b': None, 'c': r['n_other'], 'd': None,
                         'p': r['p'], 'q': r['q'], 'method': 'mann_whitney', 'n_strata': '',
                         'background': f'the log2 odds ratio of the bound RNAs of each other protein {label} tests gathering in the compartment, '
                                       f'called or not (supporting/eclip_rna_compartments). a and c count the proteins compared'})
    return rows


def protein_values(proteins: pd.DataFrame, inputs) -> pd.DataFrame:
    inputs = str(inputs)
    h = read(inputs, 'cell_map_hcm_go2021.tsv').set_index('canonical_UniProtAC')
    o = read(inputs, 'cell_map_organelle_ip_hein2025.tsv').set_index('canonical_UniProtAC')['compartment']
    c = read(inputs, 'cell_map_chromatin_kustatscher2014.tsv').set_index('canonical_UniProtAC')['interphase_chromatin_probability']
    b = read(inputs, 'cell_map_nuclear_bodies_barutcu2022.tsv').set_index('canonical_UniProtAC')['compartments']
    gathers, tested = defaultdict(list), set()
    for analysis, (_, label, _, _) in GENE_MAPS.items():
        x = bound_rna(inputs)[analysis]
        tested |= x['tested']
        for comp, accs in sorted(x['calls'].items()):
            for a in accs:
                gathers[a].append(f'{label}: {comp}')
    e = [entry(a) for a in proteins['uniprot_accession']]
    get = lambda s, a: s.get(a, '')  # noqa: E731
    return pd.DataFrame({
        'hcm_compartment': [h.at[a, 'nmf_compartment'] if a in h.index else '' for a in e],
        'hcm_domain': [h.at[a, 'safe_domain'] if a in h.index else '' for a in e],
        'organelle_ip_compartment': [get(o, a) for a in e],
        'chromatin_probability': [float(c[a]) if a in c.index else None for a in e],
        'nuclear_body_proteomes': [get(b, a) for a in e],
        'eclip_rna_compartments': ['; '.join(gathers[a]) if gathers.get(a) else 'none' if a in tested else '' for a in e],
    }, index=proteins.index)


SUPPORTING_COLUMNS = [
    ('uniprot_accession', 'string', '', 'reviewed accession of the protein with eCLIP data'),
    ('protein', 'string', '', 'its name in the eCLIP experiments'),
    ('methylated', 'boolean', '', 'the protein has a site in this release'),
    ('map', 'code', 'eclip_map', 'the map of RNA or DNA'),
    ('compartment', 'string', '', 'the compartment (a subcompartment for Hi-C)'),
    ('bound_rnas_placed', 'integer', '', 'the protein\'s bound RNAs (genes) the map places'),
    ('bound_rnas_in_compartment', 'integer', '', 'of these, those the map places in the compartment'),
    ('odds_ratio', 'number', '', 'Mantel-Haenszel odds ratio of its bound RNAs being in the compartment against the other RNAs the map places, '
     'within deciles of expression in the map\'s cells crossed with deciles of gene length (GENCODE v19; the genes v19 lacks '
     'are one more length stratum within each expression decile)'),
    ('odds_ratio_ci_low', 'number', '', 'lower bound of its 95% confidence interval (Robins-Breslow-Greenland)'),
    ('odds_ratio_ci_high', 'number', '', 'upper bound'),
    ('p', 'number', '', 'continuity-corrected Cochran-Mantel-Haenszel test, two-sided'),
    ('q', 'number', '', 'Benjamini-Hochberg over every protein and compartment of the map'),
    ('gathers', 'boolean', '', 'its bound RNAs gather in the compartment: q <= 0.05 and odds ratio above 1'),
]
MAP_CODES = {'apexseq': 'APEX-seq (Fazal et al. 2019, HEK293T, poly(A) RNA)', 'nuclear_body': 'nuclear-body APEX2-seq (Barutcu et al. 2022, Flp-In HEK293)',
             'tsaseq': 'TSA-seq, the bound RNAs\' genes (Chen et al. 2018, K562)', 'hic': 'Hi-C subcompartments, the bound RNAs\' genes (Rao et al. 2014, GM12878)'}
MAP_CODE = {'eclip_rna_apexseq': 'apexseq', 'eclip_rna_nuclear_body': 'nuclear_body', 'eclip_gene_tsaseq': 'tsaseq', 'eclip_gene_hic': 'hic'}


def supporting(sites: pd.DataFrame, inputs, version: str) -> dict:
    # the proteins the tests count as methylated: those with a site the analyses count (not only pending ones)
    methylated = {entry(a) for a in sites.loc[~sites['flags'].str.contains('one_sided_window', regex=False), 'uniprot_accession']}
    _, names = eclip_targets(str(inputs))
    rows = []
    for analysis in GENE_MAPS:
        for r in bound_rna(str(inputs))[analysis]['rows']:
            rows.append({'uniprot_accession': r['accession'], 'protein': names.get(r['accession'], ''), 'methylated': r['accession'] in methylated,
                         'map': MAP_CODE[analysis], 'compartment': r['compartment'], 'bound_rnas_placed': r['placed'],
                         'bound_rnas_in_compartment': r['in_compartment'], 'odds_ratio': r['odds_ratio'], 'odds_ratio_ci_low': r['ci_low'],
                         'odds_ratio_ci_high': r['ci_high'], 'p': r['p'], 'q': r['q'],
                         'gathers': r['q'] <= Q and r['odds_ratio'] is not None and r['odds_ratio'] > 1})
    t = pd.DataFrame(rows, columns=[c[0] for c in SUPPORTING_COLUMNS])
    return {'eclip_rna_compartments': (f'supporting/eclip_rna_compartments_{version}.tsv', t, SUPPORTING_COLUMNS)}


def changes(files: dict) -> list[str]:
    st = files['statistics'][1]
    num = lambda x: pd.to_numeric(x, errors='coerce')  # noqa: E731
    out = []
    pm = st[st['analysis'].isin(list(PROTEIN_MAPS))]
    en = pm[pm['method'] == 'enricher']
    mh = pm[pm['method'] == 'mantel_haenszel']
    both = en.merge(mh[['analysis', 'category', 'odds_ratio', 'q']].rename(columns={'odds_ratio': 'oa', 'q': 'qa'}), on=['analysis', 'category'], how='left')
    up = (num(both['q']) <= Q) & (num(both['fold_enrichment']) > 1)
    held = up & (num(both['qa']) <= Q) & (num(both['oa']) > 1)
    parts = []
    for analysis, (_, label, _, _) in PROTEIN_MAPS.items():
        g = both['analysis'] == analysis
        parts.append(f'{label}, {int((up & g).sum())} of {int(g.sum())} ({int((held & g).sum())} of them also within abundance deciles)')
    out.append('Compartments where the methylated proteins are enriched at q <= 0.05 among the proteins each map places: ' + '; '.join(parts) + '.')
    # the strongest differences among proteins of similar abundance, either way
    mh = mh.assign(o=num(mh['odds_ratio']), qq=num(mh['q']))
    name = lambda r: f"{PROTEIN_MAPS[r['analysis']][1]}: {r['category']} ({r['o']:.2g})"  # noqa: E731
    more = mh[(mh['qq'] <= Q) & (mh['o'] > 1)].sort_values('o', ascending=False).head(6)
    fewer = mh[(mh['qq'] <= Q) & (mh['o'] < 1)].sort_values('o').head(6)
    if len(more) or len(fewer):
        out.append('Among proteins of similar abundance (odds ratio within PaxDb deciles, q <= 0.05), the methylated proteins sit more often in '
                   + ('; '.join(name(r) for _, r in more.iterrows()) or 'no compartment') + ', and less often in '
                   + ('; '.join(name(r) for _, r in fewer.iterrows()) or 'none') + '.')
    gm = st[st['analysis'].isin(list(GENE_MAPS))]
    en, adj, mw = (gm[gm['method'] == m] for m in ('enricher', 'mantel_haenszel', 'mann_whitney'))
    n_all = int((num(en['a']) + num(en['b']) + num(en['c']) + num(en['d'])).max()) if len(en) else 0
    n_meth = int((num(en['a']) + num(en['b'])).max()) if len(en) else 0
    sep = lambda t: int((num(t['q']) <= Q).sum())  # noqa: E731
    low = lambda t: f'lowest q {num(t["q"]).min():.2g}' if len(t) else 'none tested'  # noqa: E731
    out.append(f'RNAs bound by proteins with eCLIP data: of the up to {n_all} such proteins a map tests, {n_meth} carry a methylarginine. '
               f'{sep(adj) or "No"} compartment{"s" if sep(adj) != 1 else ""} separate{"" if sep(adj) != 1 else "s"} them from the others within '
               f'abundance deciles ({low(adj)}), and {sep(mw) or "none"} by the size of their odds ratios ({low(mw)}).')
    return out


# ---- website: SITE_OUT/subcellular.json and subcellular_na.json ------------------------------------------------------

GEO = lambda g: f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={g}'  # noqa: E731
DATASETS = [
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
NA_DATASETS = {
    'apex': dict(key='apex', name='APEX-seq', method='RNA proximity labelling', cells='HEK293T', ref='Fazal et al. 2019', doi='10.1016/j.cell.2019.05.027',
                 how='APEX2 placed in eight compartments labels the RNAs around it; an RNA belongs to a compartment when it is labelled more than in the control (FDR < 0.05)',
                 licence='GEO GSE116008, no restrictions', url=GEO('GSE116008'), expression='the mean of the control libraries'),
    'nb': dict(key='nb', name='Nuclear-body APEX2-seq', method='RNA proximity labelling', cells='HEK293', ref='Barutcu et al. 2022', doi='10.1016/j.molcel.2021.12.010',
               how='APEX2 on markers of speckles, the lamina and nuclear bodies; an RNA belongs to a body when it is labelled at least twice as much as by nucleoplasmic APEX2 in both replicates',
               licence='GEO GSE176439, no restrictions', url=GEO('GSE176439'), expression='the mean of the input libraries'),
    'tsa': dict(key='tsa', name='TSA-seq', method='DNA distance mapping', cells='K562', ref='Chen et al. 2018', doi='10.1083/jcb.201807108',
                how='tyramide signal spreading from speckles (SON) and from the lamina (lamin B) marks DNA by its distance; a gene lies near when its signal is in the top tenth of the genome',
                licence='GEO GSE81553, no restrictions', url=GEO('GSE81553'), expression='K562 total RNA-seq (ENCODE ENCSR109IQO)'),
    'sub': dict(key='sub', name='Hi-C subcompartments', method='Chromosome contacts (Hi-C)', cells='GM12878', ref='Rao et al. 2014', doi='10.1016/j.cell.2014.11.021',
                how='the contacts each stretch of the genome makes divide it into six subcompartments: A1 and A2 active (A1 beside speckles), B1 Polycomb-repressed, '
                    'B2 at the nucleolus and the lamina with pericentromeric heterochromatin, B3 at the lamina, B4 a few clusters on chromosome 19; a gene takes the subcompartment at its middle',
                licence='GEO GSE63525, no restrictions', url=GEO('GSE63525'), expression='GM12878 total RNA-seq (ENCODE ENCSR000AEE)'),
}
G = dict
# the compartment grid: [group, [[row, {map: its category}, part of the drawings or '']]]
GRID = [
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


def l2(v):
    """log2 for the website, as abundance_strata.web_log2: +-99 where infinite or 0, None where there is none."""
    import abundance_strata as A
    return A.web_log2(v, 3)


def sig(v):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(v) else float(f'{v:.3g}')


def test_rows(st: pd.DataFrame, analysis: str, members: dict[str, int], methylated: dict[str, int]) -> list[list]:
    """[category, members, methylated, log2 fold, q, log2 OR within deciles, CI low, CI high, q] for every compartment with
    at least MIN_SET members, from the release's statistics (nulls where enricher did not test it: no methylated member)."""
    en = st[(st['analysis'] == analysis) & (st['method'] == 'enricher')].set_index('category')
    mh = st[(st['analysis'] == analysis) & (st['method'] == 'mantel_haenszel')].set_index('category')
    rows = []
    for c, m in members.items():
        if m < MIN_SET:
            continue
        if c in en.index:
            r = en.loc[c]
            if int(r['a']) != methylated[c] or int(r['a']) + int(r['c']) != m:
                raise SystemExit(f'cell_maps: {analysis} {c}: the release counts {r["a"]} of {int(r["a"]) + int(r["c"])}, the staged map {methylated[c]} of {m}')
            x = mh.loc[c] if c in mh.index else None
            rows.append([c, m, methylated[c], l2(r['fold_enrichment']), sig(r['q'])] +
                        ([l2(x['odds_ratio']), l2(x['odds_ratio_ci_low']), l2(x['odds_ratio_ci_high']), sig(x['q'])] if x is not None else [None] * 4))
        else:
            rows.append([c, m, methylated[c], None, None, None, None, None, None])
    return sorted(rows, key=lambda r: -(r[3] if r[3] is not None else -999))


def website(release: Path, version: str, inputs: Path, site_out: Path) -> None:
    import abundance_strata as A
    from release_layers import release_table, write_json
    inputs = Path(inputs)
    st = release_table(release, f'statistics_{version}.tsv')
    prot = release_table(release, f'human_methylarginine_proteins_{version}.tsv')
    meth = {entry(a) for a, n in zip(prot['uniprot_accession'], prot['n_sites']) if int(n or 0) > 0}   # as the tests count them
    pm = protein_maps(str(inputs))
    universes, tests, cats = {}, {}, {}
    for analysis, (key, _, _, _) in PROTEIN_MAPS.items():
        universe, sets = pm[analysis]
        universes[key] = universe
        tests[key] = test_rows(st, analysis, {c: len(v) for c, v in sorted(sets.items())}, {c: len(v & meth) for c, v in sorted(sets.items())})
    # the Human Protein Atlas, as the release tests it (analysis hpa_location: every protein with a location scored Approved or better)
    hpa_all = pd.read_csv(inputs / 'layers' / 'hpa_subcellular_locations.tsv', sep='\t', dtype=str, keep_default_na=False)
    confident = hpa_all[hpa_all['reliability'] != 'Uncertain']
    hpa_sets = confident.groupby('canonical_UniProtAC')['location'].agg(lambda s: sorted(set(s))).to_dict()
    universes['hpa'] = set(hpa_sets)
    by_loc = defaultdict(set)
    for a, locs in hpa_sets.items():
        for x in locs:
            by_loc[x].add(a)
    tests['hpa'] = [r for r in test_rows(st, 'hpa_location', {c: len(v) for c, v in sorted(by_loc.items())}, {c: len(v & meth) for c, v in sorted(by_loc.items())})
                    if r[1] >= 20 and r[3] is not None]
    datasets = [dict(d) for d in DATASETS]
    for d in datasets:
        d['universe'], d['methylated'] = len(universes[d['key']]), len(universes[d['key']] & meth)
    h = read(inputs, 'cell_map_hcm_go2021.tsv').set_index('canonical_UniProtAC')
    o = read(inputs, 'cell_map_organelle_ip_hein2025.tsv').set_index('canonical_UniProtAC')['compartment']
    c = read(inputs, 'cell_map_chromatin_kustatscher2014.tsv').set_index('canonical_UniProtAC')['interphase_chromatin_probability'].astype(float)
    b = read(inputs, 'cell_map_nuclear_bodies_barutcu2022.tsv').set_index('canonical_UniProtAC')['compartments']
    ap = {a: [x for x in v.split('; ') if x and x != 'none'] for a, v in b.items()}
    cats = {'hcm': sorted(set(h['nmf_compartment'])), 'hcm_safe': sorted({x if x else 'no prediction' for x in h['safe_domain']}), 'oip': sorted(set(o)),
            'hpa': sorted(set(confident['location'])), 'apexms': sorted({x for v in ap.values() for x in v})}
    # each protein's PaxDb decile among the proteins each map places, as the tests stratify (abundance_strata.strata):
    # one character per map of MAPS, 0-9 (the least to the most abundant), x without a PaxDb value, - not placed
    maps = ['hcm', 'hcm_safe', 'oip', 'chep', 'hpa', 'apexms']
    pax = A.abundance(inputs)['protein']
    dec = {k: A.strata(sorted(universes[k]), pax).to_dict() for k in maps}
    code = lambda k, a: '-' if a not in universes[k] else 'x' if dec[k][a] == 'none' else str(int(dec[k][a]) - 1)  # noqa: E731

    def icp_code(v):   # chromatin probability x1000, kept on its side of each bin edge (0.2, 0.5, 0.8)
        x = round(v * 1000)
        for edge in (200, 500, 800):
            if v < edge / 1000 and x >= edge:
                x = edge - 1
        return x
    proteins = {}
    for a in sorted(set(h.index) | set(o.index) | set(c.index) | set(hpa_sets) | set(ap)):
        hr = h.loc[a] if a in h.index else None
        proteins[a] = [cats['hcm'].index(hr['nmf_compartment']) if hr is not None else -1,
                       cats['hcm_safe'].index(hr['safe_domain'] or 'no prediction') if hr is not None else -1,
                       cats['oip'].index(o[a]) if a in o.index else -1, icp_code(float(c[a])) if a in c.index else -1,
                       [cats['hpa'].index(x) for x in hpa_sets.get(a, [])], ''.join(code(k, a) for k in maps),
                       [cats['apexms'].index(x) for x in ap.get(a, [])]]
    have = {k: {r[0] for r in v} for k, v in tests.items()}
    for _, rows in GRID:
        for row in rows:
            for k, cat in row[1].items():
                if cat not in have[k]:
                    raise SystemExit(f'cell_maps: the grid names {k} {cat!r}, which the release does not test')
    write_json(site_out / 'subcellular.json', {
        'version': version, 'datasets': datasets, 'tests': tests,
        'test_fields': ['category', 'members', 'methylated', 'log2 fold enrichment', 'q', 'log2 OR within abundance deciles', 'CI low', 'CI high', 'q'],
        'groups': GRID, 'categories': cats, 'proteins': proteins, 'decile_maps': maps,
        'protein_fields': ['hcm NMF index', 'hcm SAFE index', 'oip index', 'chromatin probability x1000', 'hpa location indices (Approved or better)',
                           'PaxDb decile among the proteins each map of decile_maps places (0-9; x no value; - not placed)', 'apexms category indices']})
    # RNA and DNA maps: the release's tests of the methylated eCLIP'd proteins against the others, the medians of their
    # odds ratios (from supporting/eclip_rna_compartments), and each protein's calls
    sup = release_table(release, f'supporting/eclip_rna_compartments_{version}.tsv')
    sup['l2'] = [log2_or(float(x)) if x != '' else None for x in sup['odds_ratio']]
    na_maps = [GENE_MAPS[a][0] for a in GENE_MAPS]
    na_ds, na_tests, na_cats = {}, {}, {}
    br = bound_rna(str(inputs))
    for analysis, (key, _, _, _) in GENE_MAPS.items():
        x = br[analysis]
        s = sup[sup['map'] == MAP_CODE[analysis]]
        tested = set(s['uniprot_accession'])
        na_cats[key] = x['categories']
        d = dict(NA_DATASETS[key])
        d.update(universe=len(tested), methylated=len(tested & meth), placed=x['placed'])
        na_ds[key] = d
        calls = s[s['gathers'] == 'true'].groupby('compartment')['uniprot_accession'].agg(set).to_dict()
        rows = test_rows(st, analysis, {cmp: len(v) for cmp, v in sorted(calls.items())}, {cmp: len(v & meth) for cmp, v in sorted(calls.items())})
        mw = st[(st['analysis'] == analysis) & (st['method'] == 'mann_whitney')].set_index('category')
        med = lambda v: None if not len(v) else (lambda m: 99 if m == math.inf else -99 if m == -math.inf else round(float(m), 3))(float(np.median(v)))  # noqa: E731
        out = []
        for r in rows:
            g = s[(s['compartment'] == r[0]) & s['l2'].notna()]
            out.append(r + [med(g.loc[g['uniprot_accession'].isin(meth), 'l2'].tolist()), med(g.loc[~g['uniprot_accession'].isin(meth), 'l2'].tolist()),
                            sig(mw.at[r[0], 'q']) if r[0] in mw.index else None])
        na_tests[key] = out
    names = dict(zip(sup['uniprot_accession'], sup['protein']))
    na_prot = {}
    for a in sorted(set(sup['uniprot_accession'])):
        per = []
        for analysis, (key, _, _, _) in GENE_MAPS.items():
            s = sup[(sup['map'] == MAP_CODE[analysis]) & (sup['uniprot_accession'] == a)]
            if not len(s):
                per.append(None)
                continue
            hit = s[(pd.to_numeric(s['q']) <= Q) & s['l2'].notna()]
            per.append([[na_cats[key].index(r.compartment), 99 if r.l2 == math.inf else -99 if r.l2 == -math.inf else round(float(r.l2), 3), sig(r.q)]
                        for r in hit.itertuples()])
        na_prot[a] = [names[a] or a, per]
    write_json(site_out / 'subcellular_na.json', {
        'version': version, 'maps': na_maps, 'datasets': na_ds, 'tests': na_tests, 'categories': na_cats, 'proteins': na_prot,
        'test_fields': ['compartment', "eCLIP'd proteins whose bound RNAs gather there", 'methylated', 'log2 fold enrichment', 'q', 'log2 OR within PaxDb deciles', 'CI low', 'CI high', 'q',
                        "median log2 OR of the methylated proteins' bound RNAs", 'of the other proteins', 'q (Mann-Whitney, BH per map)'],
        'protein_fields': ['name', 'per map of maps: [category index, log2 OR of its bound RNAs within expression and gene-length deciles, q] where q <= 0.05; null where the map does not test it']})
    print(f'cell maps: subcellular.json ({len(proteins):,} proteins), subcellular_na.json ({len(na_prot)} proteins with eCLIP data)')

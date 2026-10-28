"""Do methylated proteins bind RNAs that gather in particular compartments, more than other RNA-binding proteins do?
(pilot for the "Where in the cell" card's RNA and DNA layers; not part of a release)

Each protein with eCLIP data (ENCODE and ENCORE: the genes with most reproducible peaks or windows in each experiment,
pooled over its experiments; inputs/nucleic_acid/) is a unit. Its bound RNAs gather in a compartment when they are
placed there more often than the other RNAs (genes) a map places, within deciles of expression in the map's cells
crossed with deciles of gene length (GENCODE v19; 100 strata; Mantel-Haenszel; BH over every protein and compartment of
the map; q <= 0.05 and odds ratio > 1). Gene length matters: eCLIP target lists favour long genes, whose introns collect
more peaks, and compartments differ in it (the nuclear-body lamina RNAs come from genes of median 319 kb, the speckle
RNAs from genes of 12 kb). Compartments where a map places fewer than 20 RNAs (genes) are not tested. Then, as for every
protein-level test of the release, the methylated proteins are compared with the other eCLIP'd proteins a map can test
(those with at least 20 bound RNAs it places), within PaxDb abundance deciles (lib.tests). A call depends on how many
bound RNAs a map places, so the size of each protein's own odds ratio is compared too, whether or not it is called: the
methylated proteins' log2 odds ratios against the others' (two-sided Mann-Whitney U; BH per map).

The maps:
RNA
- APEX-seq (Fazal et al. 2019, Cell; HEK293T, poly(A) RNA; GEO GSE116008): an RNA belongs to a compartment when the
  compartment's APEX2 labelling exceeds its control at FDR < 0.05 (the study's DESeq2 table). Expression: mean CPM of
  the control libraries; RNAs with CPM >= 1 are the ones the map places.
- Nuclear-body APEX2-seq (Barutcu et al. 2022, Mol Cell; Flp-In HEK293, ribosome-depleted RNA; GEO GSE176439): an RNA
  belongs to a body when the marker's labelling over its input exceeds the nucleoplasmic NLS-APEX2 control's by at least
  twofold in both replicates; speckles need two of the three speckle markers (SRSF1, SRSF7, RNPS1) and PML bodies both
  PML and SP100. Expression: mean CPM of the inputs (>= 1).
DNA
- TSA-seq (Chen et al. 2018, J Cell Biol; K562; GEO GSE81553): a gene lies near nuclear speckles (or the nuclear lamina)
  when its SON (lamin B) TSA-seq signal, the mean of the two replicates over the gene and at least 20 kb, is in the top
  tenth of the genome's 20 kb bins. Genes expressed in K562 (total RNA-seq TPM >= 1, ENCODE ENCSR109IQO, the mean of the
  two replicates) are the ones the map places; that TPM is also their expression.
- Hi-C subcompartments (Rao et al. 2014, Cell; GM12878; GEO GSE63525, GSE63525_GM12878_subcompartments.bed.gz): the
  subcompartment (A1, A2, B1, B2, B3, B4) at the middle of a gene expressed in GM12878 (total RNA-seq TPM >= 1, ENCODE
  ENCSR000AEE, the mean of the two replicates; also its expression).

Genes are joined by Ensembl gene ID (version dropped); nuclear-body counts carry symbols only and are joined by a symbol
that names a single gene.

    python pilot_layers.py REPO DL OUT
DL holds apex/, barutcu/ (GSE176439_NuclearBody_APEX_seq_gene_counts.txt.gz), tsa/ (the GSE81553 bigWigs and
gencode.v19.annotation.gtf.gz), k562/ (ENCFF611MXW.tsv, ENCFF384BFE.tsv), hic/ (GSE63525_GM12878_subcompartments.bed.gz)
and gm12878/ (ENCFF723ICA.tsv, ENCFF362RMV.tsv). Writes OUT/layer_tests.tsv and OUT/layer_proteins.json.
"""
import sys, json, gzip, re
from collections import defaultdict
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
import lib

repo, dl, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
ctx = lib.setup(repo)
nuc = Path(repo) / 'inputs' / 'nucleic_acid'
strip = lambda g: str(g).split('.')[0]  # noqa: E731

# ---- eCLIP targets: protein (reviewed accession) -> Ensembl gene IDs; and the gene symbols they name
names = {}
targets = defaultdict(set)
rbp_name = {}
for src, key, sym in (('encode', 'encode_experiment', 'target'), ('encore', 'encore_experiment', 'rbp')):
    exp = pd.read_csv(nuc / f'{src}_eclip_experiments.tsv', sep='\t', dtype=str, keep_default_na=False)
    for a, n in zip(exp['uniprot_accession'], exp[sym]):
        if a:
            rbp_name.setdefault(a.split('-')[0], n)
    ex = exp.set_index(key)['uniprot_accession']
    tg = pd.read_csv(nuc / f'{src}_eclip_target_genes.tsv', sep='\t', dtype=str, keep_default_na=False)
    for e, g, n in zip(tg[key], tg['gene_id'], tg['gene_name']):
        acc = ex.get(e, '').split('-')[0]
        if acc:
            targets[acc].add(strip(g)); names.setdefault(strip(g), n)
meth = ctx['meth']
print(f'eCLIP proteins {len(targets)}: methylated {len(set(targets) & meth)}, without a methylarginine {len(set(targets) - meth)}')

genes = {}   # Ensembl gene ID -> {'cats': {map: [categories]}, 'expr': {map: value}}
def put(g, key, cats, expr):
    x = genes.setdefault(g, {'cats': {}, 'expr': {}})
    x['cats'][key] = sorted(set(x['cats'].get(key, [])) | set(cats))
    x['expr'][key] = max(expr, x['expr'].get(key, -1.0))

# ---- APEX-seq (Fazal 2019)
APEX = {'NLS': 'Nucleus', 'NIK': 'Nucleolus', 'LMA': 'Nuclear lamina', 'NUCP': 'Nuclear pore', 'NES': 'Cytosol', 'ERM': 'ER membrane',
        'OMM': 'Outer mitochondrial membrane', 'MITO': 'Mitochondrial matrix'}
de = pd.read_csv(dl / 'apex' / 'GSE116008_APEX_seq_deseq2_all_data_polyA_hek293t.txt.gz', sep='\t').set_index('Ensembl_Gene')
cnt = pd.read_csv(dl / 'apex' / 'GSE116008_APEX_seq_counts_all_data_polyA_hek293t.txt.gz', sep='\t', index_col=0).apply(pd.to_numeric, errors='coerce').fillna(0)
cpm = cnt / cnt.sum() * 1e6
expr = cpm[[c for c in cpm.columns if re.search(r'-C\d$', c)]].mean(axis=1)
for ens, row in de.iterrows():
    e = float(expr.get(ens, 0.0))
    if e < 1:
        continue
    if isinstance(row['Common_Gene'], str):
        names.setdefault(strip(ens), row['Common_Gene'])
    put(strip(ens), 'apex', [n for code, n in APEX.items() if row[code + '.18C.P'] < 0.05 and row[code + '.18C'] > 0], e)

# ---- symbols -> a single Ensembl gene (GENCODE v19, APEX-seq's Ensembl 84, eCLIP's GENCODE v50)
gtf = []
with gzip.open(dl / 'tsa' / 'gencode.v19.annotation.gtf.gz', 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.split('\t')
        if f[2] == 'gene':
            gtf.append((f[0], int(f[3]) - 1, int(f[4]), re.search(r'gene_id "([^".]+)', f[8]).group(1), re.search(r'gene_name "([^"]+)"', f[8]).group(1)))
by_symbol = defaultdict(set)
for _, _, _, g, n in gtf:
    by_symbol[n].add(g)
for ens, n in de['Common_Gene'].items():
    by_symbol[str(n)].add(strip(ens))
for g, n in list(names.items()):
    by_symbol[n].add(g)
for _, _, _, g, n in gtf:
    names.setdefault(g, n)
glen = pd.Series({g: b - a for _, a, b, g, _ in gtf})  # gene length (bp); genes GENCODE v19 lacks have none

# ---- Nuclear-body APEX2-seq (Barutcu 2022)
c = pd.read_csv(dl / 'barutcu' / 'GSE176439_NuclearBody_APEX_seq_gene_counts.txt.gz', sep='\t', index_col=0, low_memory=False).apply(pd.to_numeric, errors='coerce').fillna(0)
c.columns = [x.replace('_Rep', '-rep').replace('_rep', '-rep').replace('Rep', 'rep') for x in c.columns]
bcpm = c / c.sum() * 1e6
bexpr = bcpm[[x for x in bcpm.columns if 'input' in x.lower()]].mean(axis=1)
col = lambda name, rep: bcpm[f'{name}-rep{rep}']  # noqa: E731
enr = lambda bait, rep: np.log2((col(bait, rep) + 0.5) / (col(bait + 'input', rep) + 0.5))  # noqa: E731
over = {b: (enr(b, 1) - enr('NLS', 1) >= 1) & (enr(b, 2) - enr('NLS', 2) >= 1) for b in ('SRSF1', 'SRSF7', 'RNPS1', 'LMNA', 'FBL', 'PML', 'SP100', 'WRAP53', 'SMN2', 'NPAT', 'SAM68')}
nb = {'Nuclear speckles': (over['SRSF1'].astype(int) + over['SRSF7'] + over['RNPS1']) >= 2, 'Nuclear lamina': over['LMNA'], 'Nucleolus': over['FBL'],
      'PML bodies': over['PML'] & over['SP100'], 'Cajal bodies': over['WRAP53'], 'Gems': over['SMN2'], 'Histone locus bodies': over['NPAT'], 'Sam68 bodies': over['SAM68']}
n_nb = 0
for sym, e in bexpr.items():
    g = by_symbol.get(str(sym), set())
    if len(g) != 1 or e < 1:
        continue
    put(next(iter(g)), 'nb', [k for k, v in nb.items() if v.get(sym, False)], float(e)); n_nb += 1

# ---- TSA-seq (Chen 2018), genes expressed in K562
import pyBigWig
tracks = {'son': ['GSE81553_SON_TSA-Seq_Condition2.bw', 'GSE81553_SON_TSA-Seq_Condition2.Replicate.bw'],
          'lmnb': ['GSE81553_LaminB_TSA-Seq_Condition2.bw', 'GSE81553_LaminB_TSA-Seq_Conditoin2_Replicate.bw']}
bw = {k: [pyBigWig.open(str(dl / 'tsa' / f)) for f in fs] for k, fs in tracks.items()}
CHROMS = [f'chr{i}' for i in range(1, 23)] + ['chrX']
BIN = 20000
def mean(k, chrom, a, b):
    v = [x for x in (w.stats(chrom, a, b, type='mean')[0] for w in bw[k]) if x is not None]
    return float(np.mean(v)) if v else np.nan
edges = {}
for k in tracks:
    vals = []
    for chrom in CHROMS:
        nbins = bw[k][0].chroms()[chrom] // BIN
        reps = [np.array([np.nan if x is None else x for x in w.stats(chrom, 0, nbins * BIN, type='mean', nBins=nbins)], dtype=float) for w in bw[k]]
        vals.append(np.nanmean(np.vstack(reps), axis=0))
    v = np.concatenate(vals); edges[k] = float(np.quantile(v[np.isfinite(v)], 0.9))
def encode_tpm(folder, files):
    t = pd.concat([pd.read_csv(dl / folder / f, sep='\t', usecols=['gene_id', 'TPM']).set_index('gene_id')['TPM'] for f in files], axis=1).mean(axis=1)
    t.index = [strip(i) for i in t.index]
    return t.groupby(level=0).max()
tpm = encode_tpm('k562', ('ENCFF611MXW.tsv', 'ENCFF384BFE.tsv'))
n_tsa = 0
for chrom, a, b, g, n in gtf:
    e = float(tpm.get(g, 0.0))
    if chrom not in CHROMS or e < 1:
        continue
    if b - a < BIN:
        m = (a + b) // 2; a, b = max(0, m - BIN // 2), m + BIN // 2
    s, l = mean('son', chrom, a, b), mean('lmnb', chrom, a, b)
    if not (np.isfinite(s) and np.isfinite(l)):
        continue
    put(g, 'tsa', (['Near nuclear speckles'] if s >= edges['son'] else []) + (['Near the nuclear lamina'] if l >= edges['lmnb'] else []), e); n_tsa += 1
# ---- Hi-C subcompartments (Rao 2014), genes expressed in GM12878
sub = pd.read_csv(dl / 'hic' / 'GSE63525_GM12878_subcompartments.bed.gz', sep='\t', header=None, usecols=[0, 1, 2, 3], names=['chrom', 'a', 'b', 'sub'])
sub = sub[sub['sub'].isin(['A1', 'A2', 'B1', 'B2', 'B3', 'B4'])]
by_chrom = {c: (g['a'].to_numpy(), g['b'].to_numpy(), g['sub'].to_numpy()) for c, g in sub.sort_values(['chrom', 'a']).groupby('chrom')}
gm = encode_tpm('gm12878', ('ENCFF723ICA.tsv', 'ENCFF362RMV.tsv'))
SUBNAME = {'A1': 'A1 (speckle-associated)', 'A2': 'A2 (active)', 'B1': 'B1 (Polycomb)', 'B2': 'B2 (perinucleolar, pericentromeric)', 'B3': 'B3 (lamina-associated)', 'B4': 'B4 (chromosome 19 KZFP clusters)'}
n_sub = 0
for chrom, a, b, g, n in gtf:
    e = float(gm.get(g, 0.0))
    if chrom not in by_chrom or e < 1:
        continue
    starts, ends, lab = by_chrom[chrom]; m = (a + b) // 2; i = np.searchsorted(starts, m, side='right') - 1
    if i >= 0 and m < ends[i]:
        put(g, 'sub', [SUBNAME[lab[i]]], e); n_sub += 1
print(f'genes placed: APEX-seq {sum("apex" in g["cats"] for g in genes.values())}, nuclear bodies {n_nb}, TSA-seq {n_tsa}, Hi-C subcompartments {n_sub}; genes {len(genes)}')

# ---- each eCLIP'd protein: do the RNAs it binds gather in each compartment? (within expression deciles; BH per map)
MAPS = ('apex', 'nb', 'tsa', 'sub')
per = defaultdict(dict)   # accession -> map -> [[compartment, log2 OR, q], ...] where q <= 0.05
tested = defaultdict(set)  # map -> accessions with at least 20 bound RNAs it places
calls, tested_cats, sizes = {}, {}, {}
size = defaultdict(lambda: defaultdict(dict))  # map -> compartment -> accession -> log2 OR (-inf for 0, inf for infinite)
for key in MAPS:
    universe = sorted(g for g, x in genes.items() if key in x['cats'])
    pos = {g: i for i, g in enumerate(universe)}
    by_expr = ctx['A'].strata(universe, pd.Series({g: genes[g]['expr'][key] for g in universe})).to_numpy(dtype=str)
    by_len = ctx['A'].strata(universe, glen).to_numpy(dtype=str)
    codes, layer = np.unique(np.char.add(np.char.add(by_expr, '|'), by_len), return_inverse=True)
    K = len(codes); n_all = np.bincount(layer, minlength=K)
    in_cat = {c: np.array([c in genes[g]['cats'][key] for g in universe]) for c in sorted({c for g in universe for c in genes[g]['cats'][key]})}
    sizes[key] = {c: int(v.sum()) for c, v in in_cat.items()}
    print(f'{key}: RNAs (genes) placed in each compartment ' + ', '.join(f'{c} {n}' for c, n in sizes[key].items()))
    cats = tested_cats[key] = [c for c, v in in_cat.items() if v.sum() >= 20]
    rows = []
    for acc, t in targets.items():
        idx = np.array(sorted(pos[g] for g in t if g in pos), dtype=int)
        if len(idx) < 20:
            continue
        tested[key].add(acc)
        sel = np.zeros(len(universe), dtype=bool); sel[idx] = True
        n_sel = np.bincount(layer[sel], minlength=K)
        for c in cats:
            n_cat = np.bincount(layer[in_cat[c]], minlength=K); both = np.bincount(layer[sel & in_cat[c]], minlength=K)
            tables = [(int(both[i]), int(n_sel[i] - both[i]), int(n_cat[i] - both[i]), int(n_all[i] - n_sel[i] - n_cat[i] + both[i])) for i in range(K)]
            odds, lo, hi, p, used = ctx['mh'](tables)
            rows.append([acc, c, odds, p])
            if odds is not None:
                size[key][c][acc] = float(np.log2(odds)) if odds > 0 else -np.inf
    q = ctx['bh']([r[3] if r[3] is not None else 1.0 for r in rows])
    calls[key] = defaultdict(set)
    for (acc, c, odds, p), qq in zip(rows, q):
        if qq <= 0.05 and odds is not None and odds > 0:
            per[acc].setdefault(key, []).append([c, round(float(np.log2(odds)), 3) if np.isfinite(odds) else 99, float(f'{qq:.3g}')])
            if odds > 1:
                calls[key][c].add(acc)
    print(f'{key}: {len(tested[key])} proteins tested; calls per compartment ' + ', '.join(f'{c} {len(v)}' for c, v in sorted(calls[key].items())))

# ---- the methylated proteins against the other eCLIP'd proteins each map tests, within PaxDb abundance deciles
from scipy.stats import mannwhitneyu
out_rows = []
for key in MAPS:
    t = lib.tests(ctx, tested[key], {c: v for c, v in calls[key].items()}, min_size=5)
    t.insert(0, 'map', key)
    med = lambda v: float(np.median(v)) if len(v) else np.nan  # noqa: E731
    yes = {c: [v for a, v in size[key][c].items() if a in meth] for c in t['compartment']}
    no = {c: [v for a, v in size[key][c].items() if a not in meth] for c in t['compartment']}
    t['log2_or_median_methylated'] = [med(yes[c]) for c in t['compartment']]
    t['log2_or_median_other'] = [med(no[c]) for c in t['compartment']]
    t['p_size'] = [float(mannwhitneyu(yes[c], no[c]).pvalue) if yes[c] and no[c] else np.nan for c in t['compartment']]
    t['q_size'] = ctx['bh'](t['p_size'].fillna(1))
    print(t[['compartment', 'log2_or_median_methylated', 'log2_or_median_other', 'p_size', 'q_size']].to_string(index=False, float_format=lambda v: f'{v:.3g}'))
    out_rows.append(t)
    lib.show(t, f'{key}: methylated among the proteins whose bound RNAs gather in each compartment ({len(tested[key] & meth)} of {len(tested[key])} methylated)')
pd.concat(out_rows).to_csv(out / 'layer_tests.tsv', sep='\t', index=False)
json.dump({'proteins': {a: {'name': rbp_name.get(a, a), 'tested': sorted(k for k in MAPS if a in tested[k]), 'calls': per.get(a, {})} for a in sorted(targets)},
           'categories': tested_cats, 'sizes': sizes,
           'placed': {k: sum(k in x['cats'] for x in genes.values()) for k in MAPS}}, open(out / 'layer_proteins.json', 'w'))

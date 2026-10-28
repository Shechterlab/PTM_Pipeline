"""Which compartments hold the RNAs that methylated proteins bind (pilot for the "Where in the cell" card's RNA and DNA
layers; not part of a release).

The RNAs bound by a protein are its eCLIP targets: the genes with most reproducible peaks (ENCODE, K562 and HepG2) or
windows (ENCORE) in each experiment (inputs/nucleic_acid/encode_eclip_target_genes.tsv, encore_eclip_target_genes.tsv),
pooled over its experiments. The compartments are those of three maps:

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

Genes are joined by Ensembl gene ID (version dropped); nuclear-body counts carry symbols only and are joined by a symbol
that names a single gene. Each map compares the RNAs bound by methylated proteins with the other RNAs (genes) it places,
within deciles of their expression in the same cells (lib.tests); the RNAs bound by eCLIP'd proteins without a
methylarginine are compared the same way, for reference.

    python pilot_layers.py REPO DL OUT
DL holds apex/, barutcu/ (GSE176439_NuclearBody_APEX_seq_gene_counts.txt.gz), tsa/ (the GSE81553 bigWigs and
gencode.v19.annotation.gtf.gz) and k562/ (ENCFF611MXW.tsv, ENCFF384BFE.tsv). Writes OUT/layer_tests.tsv and
OUT/layer_genes.json.
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
tpm = pd.concat([pd.read_csv(dl / 'k562' / f, sep='\t', usecols=['gene_id', 'TPM']).set_index('gene_id')['TPM'] for f in ('ENCFF611MXW.tsv', 'ENCFF384BFE.tsv')], axis=1).mean(axis=1)
tpm.index = [strip(i) for i in tpm.index]
tpm = tpm.groupby(level=0).max()
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
print(f'genes placed: APEX-seq {sum("apex" in g["cats"] for g in genes.values())}, nuclear bodies {n_nb}, TSA-seq {n_tsa}; genes {len(genes)}; TSA top-tenth edges {edges}')

# ---- tests: RNAs bound by methylated proteins (and, for reference, by eCLIP'd proteins without a methylarginine)
bound = {'methylated': set().union(*[targets[a] for a in targets if a in meth]), 'other': set().union(*[targets[a] for a in targets if a not in meth])}
rows = []
for key in ('apex', 'nb', 'tsa'):
    universe = {g for g, x in genes.items() if key in x['cats']}
    values = pd.Series({g: genes[g]['expr'][key] for g in universe})
    sets = defaultdict(set)
    for g in universe:
        for cat in genes[g]['cats'][key]:
            sets[cat].add(g)
    for who, s in bound.items():
        t = lib.tests(ctx, universe, sets, min_size=10, values=values, target=s)
        t.insert(0, 'bound_by', who); t.insert(0, 'map', key)
        rows.append(t)
        lib.show(t, f'{key}: RNAs bound by {who} proteins ({len(s & universe)} of {len(universe)})')
pd.concat(rows).to_csv(out / 'layer_tests.tsv', sep='\t', index=False)
json.dump({'genes': genes, 'names': {g: names.get(g, g) for g in genes}, 'targets': {a: sorted(s & set(genes)) for a, s in targets.items()}, 'rbps': rbp_name}, open(out / 'layer_genes.json', 'w'))

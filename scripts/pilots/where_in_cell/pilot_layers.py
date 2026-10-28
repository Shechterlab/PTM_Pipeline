"""Where the transcripts and the genes of methylated proteins sit (pilot for the "Where in the cell" card's RNA and DNA
layers; not part of a release).

RNA
- APEX-seq (Fazal et al. 2019, Cell; HEK293T, poly(A) RNA; GEO GSE116008): a transcript is enriched in a compartment
  when the compartment's APEX2 labelling exceeds its control at FDR < 0.05 (the study's DESeq2 table). Expression: mean
  CPM of the control libraries; transcripts with CPM >= 1 are the ones the map places.
- Nuclear-body APEX2-seq (Barutcu et al. 2022, Mol Cell; Flp-In HEK293, ribosome-depleted RNA; GEO GSE176439): a
  transcript is associated with a body when the bait's labelling over its input exceeds the nucleoplasmic NLS-APEX2
  control's by at least twofold in both replicates; speckles need two of the three speckle baits (SRSF1, SRSF7, RNPS1)
  and PML bodies both PML and SP100. Expression: mean CPM of the inputs (>= 1).
DNA
- TSA-seq (Chen et al. 2018, J Cell Biol; K562; GEO GSE81553): a protein-coding gene is near nuclear speckles (or the
  nuclear lamina) when its SON (lamin B) TSA-seq signal, the mean of the two replicates over the gene and at least 20 kb,
  is in the top tenth of the genome's 20 kb bins. Expression: K562 total RNA-seq TPM (ENCODE ENCSR109IQO, mean of the two
  replicates).

Each gene is placed on its reviewed UniProt entry by its gene symbol (lib.resolve). The tests (lib.tests) compare the
genes of methylated proteins with the other genes each dataset places, within deciles of their expression in the same
cells.

    python pilot_layers.py REPO DL OUT
DL holds apex/, barutcu/ (GSE176439_NuclearBody_APEX_seq_gene_counts.txt.gz), tsa/ (the GSE81553 bigWigs and
gencode.v19.annotation.gtf.gz) and k562/ (ENCFF611MXW.tsv, ENCFF384BFE.tsv). Writes OUT/layer_tests.tsv and
OUT/layer_genes.json.
"""
import sys, json, gzip, re
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
import lib

repo, dl, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
ctx = lib.setup(repo)
acc_of = lambda symbol: lib.resolve(ctx, [], symbol) if isinstance(symbol, str) and symbol else ''  # noqa: E731
genes = {}   # accession -> {'cats': {map: [categories]}, 'expr': {map: value}}
def put(acc, key, cats, expr):
    g = genes.setdefault(acc, {'cats': {}, 'expr': {}})
    g['cats'][key] = sorted(set(g['cats'].get(key, [])) | set(cats))
    g['expr'][key] = max(expr, g['expr'].get(key, -1.0))

# ---- APEX-seq (Fazal 2019)
APEX = {'NLS': 'Nucleus', 'NIK': 'Nucleolus', 'LMA': 'Nuclear lamina', 'NUCP': 'Nuclear pore', 'NES': 'Cytosol', 'ERM': 'ER membrane',
        'OMM': 'Outer mitochondrial membrane', 'MITO': 'Mitochondrial matrix'}
de = pd.read_csv(dl / 'apex' / 'GSE116008_APEX_seq_deseq2_all_data_polyA_hek293t.txt.gz', sep='\t')
cnt = pd.read_csv(dl / 'apex' / 'GSE116008_APEX_seq_counts_all_data_polyA_hek293t.txt.gz', sep='\t', index_col=0)
cnt = cnt.apply(pd.to_numeric, errors='coerce').fillna(0)
cpm = cnt / cnt.sum() * 1e6
expr = cpm[[c for c in cpm.columns if re.search(r'-C\d$', c)]].mean(axis=1)
de = de.set_index('Ensembl_Gene')
n_apex = 0
for ens, row in de.iterrows():
    e = float(expr.get(ens, 0.0))
    acc = acc_of(row['Common_Gene'])
    if not acc or e < 1:
        continue
    cats = [name for code, name in APEX.items() if row[code + '.18C.P'] < 0.05 and row[code + '.18C'] > 0]
    put(acc, 'apex', cats, e); n_apex += 1

# ---- Nuclear-body APEX2-seq (Barutcu 2022)
c = pd.read_csv(dl / 'barutcu' / 'GSE176439_NuclearBody_APEX_seq_gene_counts.txt.gz', sep='\t', index_col=0, low_memory=False).apply(pd.to_numeric, errors='coerce').fillna(0)
c.columns = [x.replace('_Rep', '-rep').replace('_rep', '-rep').replace('Rep', 'rep') for x in c.columns]
bcpm = c / c.sum() * 1e6
inputs = [x for x in bcpm.columns if 'input' in x.lower()]
bexpr = bcpm[inputs].mean(axis=1)
col = lambda name, rep: bcpm[f'{name}-rep{rep}']  # noqa: E731
enr = lambda bait, rep: np.log2((col(bait, rep) + 0.5) / (col(bait + 'input', rep) + 0.5))  # noqa: E731
over = {b: (enr(b, 1) - enr('NLS', 1) >= 1) & (enr(b, 2) - enr('NLS', 2) >= 1) for b in ('SRSF1', 'SRSF7', 'RNPS1', 'LMNA', 'FBL', 'PML', 'SP100', 'WRAP53', 'SMN2', 'NPAT', 'SAM68')}
nb = {'Nuclear speckles': (over['SRSF1'].astype(int) + over['SRSF7'] + over['RNPS1']) >= 2, 'Nuclear lamina': over['LMNA'], 'Nucleolus': over['FBL'],
      'PML bodies': over['PML'] & over['SP100'], 'Cajal bodies': over['WRAP53'], 'Gems': over['SMN2'], 'Histone locus bodies': over['NPAT'], 'Sam68 bodies': over['SAM68']}
n_nb = 0
for sym, e in bexpr.items():
    acc = acc_of(sym)
    if not acc or e < 1:
        continue
    put(acc, 'nb', [k for k, v in nb.items() if v.get(sym, False)], float(e)); n_nb += 1

# ---- TSA-seq (Chen 2018)
import pyBigWig
tsa = dl / 'tsa'
tracks = {'son': ['GSE81553_SON_TSA-Seq_Condition2.bw', 'GSE81553_SON_TSA-Seq_Condition2.Replicate.bw'],
          'lmnb': ['GSE81553_LaminB_TSA-Seq_Condition2.bw', 'GSE81553_LaminB_TSA-Seq_Conditoin2_Replicate.bw']}
bw = {k: [pyBigWig.open(str(tsa / f)) for f in fs] for k, fs in tracks.items()}
CHROMS = [f'chr{i}' for i in range(1, 23)] + ['chrX']
BIN = 20000
def mean(k, chrom, a, b):
    v = [w.stats(chrom, a, b, type='mean')[0] for w in bw[k]]
    v = [x for x in v if x is not None]
    return float(np.mean(v)) if v else np.nan
edges = {}
for k in tracks:
    vals = []
    for chrom in CHROMS:
        L = bw[k][0].chroms()[chrom]; nb_ = L // BIN
        reps = [np.array([np.nan if x is None else x for x in w.stats(chrom, 0, nb_ * BIN, type='mean', nBins=nb_)], dtype=float) for w in bw[k]]
        vals.append(np.nanmean(np.vstack(reps), axis=0))
    v = np.concatenate(vals); v = v[np.isfinite(v)]
    edges[k] = float(np.quantile(v, 0.9))
tpm = pd.concat([pd.read_csv(dl / 'k562' / f, sep='\t', usecols=['gene_id', 'TPM']).set_index('gene_id')['TPM'] for f in ('ENCFF611MXW.tsv', 'ENCFF384BFE.tsv')], axis=1).mean(axis=1)
tpm.index = [str(i).split('.')[0] for i in tpm.index]
tpm = tpm.groupby(level=0).max()
n_tsa = 0
with gzip.open(tsa / 'gencode.v19.annotation.gtf.gz', 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.split('\t')
        if f[2] != 'gene' or f[0] not in CHROMS or 'gene_type "protein_coding"' not in f[8]:
            continue
        name = re.search(r'gene_name "([^"]+)"', f[8]).group(1); ens = re.search(r'gene_id "([^".]+)', f[8]).group(1)
        acc = acc_of(name)
        if not acc:
            continue
        a, b = int(f[3]) - 1, int(f[4])
        if b - a < BIN:
            m = (a + b) // 2; a, b = max(0, m - BIN // 2), m + BIN // 2
        s, l = mean('son', f[0], a, b), mean('lmnb', f[0], a, b)
        if not (np.isfinite(s) and np.isfinite(l)):
            continue
        cats = (['Near nuclear speckles'] if s >= edges['son'] else []) + (['Near the nuclear lamina'] if l >= edges['lmnb'] else [])
        put(acc, 'tsa', cats, float(tpm.get(ens, np.nan)) if ens in tpm.index else -1.0); n_tsa += 1
print(f'genes placed: APEX-seq {n_apex}, nuclear bodies {n_nb}, TSA-seq {n_tsa}; accessions {len(genes)}; TSA top-decile edges {edges}')

# ---- tests: the genes of methylated proteins against the other genes each dataset places, within expression deciles
rows = []
for key in ('apex', 'nb', 'tsa'):
    universe = {a for a, g in genes.items() if key in g['cats']}
    values = pd.Series({a: genes[a]['expr'][key] for a in universe}).where(lambda s: s >= 0)
    sets = {}
    for a in universe:
        for cat in genes[a]['cats'][key]:
            sets.setdefault(cat, set()).add(a)
    t = lib.tests(ctx, universe, sets, values=values)
    t.insert(0, 'map', key)
    rows.append(t)
    lib.show(t, key)
pd.concat(rows).to_csv(out / 'layer_tests.tsv', sep='\t', index=False)
json.dump(genes, open(out / 'layer_genes.json', 'w'))

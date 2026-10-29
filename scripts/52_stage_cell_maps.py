#!/usr/bin/env python3
"""Stage the maps of the cell for the cell-map layer (scripts/cell_maps.py): where proteins sit, by four protein maps,
and where RNAs gather and genes lie, by two RNA maps and two DNA maps.

Protein maps (one row per reviewed entry the map places; an identifier list names its first reviewed accession, else
the one entry its gene symbol names; an entry named by several rows keeps the first that names its accession, else the
first that names its symbol, and for chromatin its highest probability. Several rows name one accession only for
identical copies, such as the histone H3 and H4 genes and SMN1 and SMN2 in the Human Cell Map)
  cell_map_hcm_go2021.tsv          Human Cell Map (Go et al. 2021, Nature 595:120; HEK293, BioID from 192 baits): each
                                   localised prey's NMF compartment and its SAFE domain (empty where SAFE gives none).
                                   humancellmap.org preys-v1.txt; CC BY 4.0
  cell_map_organelle_ip_hein2025.tsv  organelle IP (Hein et al. 2025, Cell 188:1137; HEK293T): each protein group's
                                   cluster annotation ('unclassified' kept). Organelle_IP_analyses_and_figures
                                   data/labels/cluster_annotation_Dec6.csv; BSD-3
  cell_map_chromatin_kustatscher2014.tsv  interphase chromatin probability (Kustatscher et al. 2014, EMBO J 33:648,
                                   Table EV1 = embj0033-0648-sd8.xlsx). CC BY
  cell_map_nuclear_bodies_barutcu2022.tsv  the speckle, lamina and PML-body APEX2-MS proteomes (Barutcu et al. 2022, Mol
                                   Cell 82:1035, Table S1; Flp-In HEK293, SAINT): a prey belongs to nuclear speckles when
                                   its average probability (AvgP) is >= 0.9 with at least two of SRSF1, SRSF7 and RNPS1, to
                                   the nuclear lamina with LMNA or LMNB1, to PML bodies with SP100; every prey detected is
                                   placed ('none' where it belongs to none). Preys are named by gene symbol only. The
                                   table holds these calls, not the supplementary table (Elsevier; not redistributed)

RNA and DNA maps (one row per gene the map places, by Ensembl gene ID without version; gene_length from GENCODE v19,
empty for genes v19 lacks; compartments '; '-separated, empty for placed genes in none)
  rna_map_apexseq_fazal2019.tsv    APEX-seq (Fazal et al. 2019, Cell 178:473; HEK293T, poly(A); GEO GSE116008): an RNA
                                   belongs to a compartment when its labelling exceeds the control at FDR < 0.05 with a
                                   positive log2 change (the study's DESeq2 table, 18C columns). expression_cpm: mean CPM
                                   of the control libraries; RNAs with CPM >= 1 are placed
  rna_map_nuclear_bodies_barutcu2022.tsv  nuclear-body APEX2-seq (Barutcu et al. 2022; Flp-In HEK293, ribosome-depleted;
                                   GEO GSE176439): a body's RNA is labelled at least twice as much over its input as by
                                   nucleoplasmic NLS-APEX2 in both replicates (speckles: two of SRSF1, SRSF7, RNPS1; PML
                                   bodies: both PML and SP100). Counts are by symbol: a symbol naming a single gene
                                   (GENCODE v19, the APEX-seq table, or the staged eCLIP target genes) is placed.
                                   expression_cpm: mean CPM of the inputs (>= 1)
  dna_map_tsaseq_chen2018.tsv      TSA-seq (Chen et al. 2018, J Cell Biol 217:4025; K562; GEO GSE81553, condition 2,
                                   two replicates): son_tsa and laminb_tsa, the mean signal of the replicates over the gene
                                   (at least 20 kb about its middle); a gene lies near speckles (the lamina) when that is
                                   in the top tenth of the genome's 20 kb bins (chromosomes 1-22 and X). Genes expressed in
                                   K562 (ENCODE ENCSR109IQO total RNA-seq, mean TPM of ENCFF611MXW and ENCFF384BFE >= 1)
  dna_map_hic_rao2014.tsv          Hi-C subcompartments (Rao et al. 2014, Cell 159:1665; GM12878; GEO GSE63525): the
                                   subcompartment at the middle of a gene expressed in GM12878 (ENCODE ENCSR000AEE, mean
                                   TPM of ENCFF723ICA and ENCFF362RMV >= 1); hg19 coordinates, as GENCODE v19
  cell_maps_summary.json           counts, the TSA-seq edges, versions and the MD5 of every input

Outputs go to --outdir, with manifest_rows.tsv and raw_files_rows.tsv (the rows these files take in inputs/MANIFEST.tsv
and inputs/layers/raw_files.tsv). --install then copies the tables into inputs/layers/ and puts the rows in place.

    python scripts/52_stage_cell_maps.py --downloads PTM_data/external/cell_maps --outdir /tmp/cell_maps [--install]

Needs pyBigWig. Running it again on the same downloads gives the same bytes.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import shutil
import zipfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from common import md5_file, refresh_manifest_row, update_manifest, upsert_tsv_rows

GEO = 'https://ftp.ncbi.nlm.nih.gov/geo/series'
ENCODE = 'https://www.encodeproject.org/files/{0}/@@download/{0}.tsv'
OIP_COMMIT = '65f22ac57deeab5e0d5b8959049e88e36c3ad876'
# file: (URL, what it is, how it is kept)
DOWNLOADS = {
    'preys-v1.txt': ('https://humancellmap.org/resources/downloads/preys-v1.txt', 'Human Cell Map v1 preys (NMF and SAFE localisations)',
                     'release_asset (CC BY 4.0)'),
    'cluster_annotation_Dec6.csv': (f'https://raw.githubusercontent.com/czbiohub-sf/Organelle_IP_analyses_and_figures/{OIP_COMMIT}/data/labels/cluster_annotation_Dec6.csv',
                                    f'organelle IP cluster annotation, repository commit {OIP_COMMIT[:7]} (2026-03-10)', 'release_asset (BSD-3)'),
    'supp.zip': ('https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3983682/supplementaryFiles',
                 'Kustatscher et al. 2014 supplementary files (Europe PMC); member embj0033-0648-sd8.xlsx read', 'release_asset (CC BY, the article\'s licence)'),
    '1-s2.0-S1097276521010728-mmc2.xlsx': ('https://ars.els-cdn.com/content/image/1-s2.0-S1097276521010728-mmc2.xlsx',
                                           'Barutcu et al. 2022 Table S1 (APEX2-MS SAINT scores)', 'checksum_only (Elsevier; the article is not open access)'),
    'GSE116008_APEX_seq_deseq2_all_data_polyA_hek293t.txt.gz': (f'{GEO}/GSE116nnn/GSE116008/suppl/GSE116008_APEX_seq_deseq2_all_data_polyA_hek293t.txt.gz',
                                                                'APEX-seq DESeq2 table, poly(A), HEK293T', 'release_asset (GEO)'),
    'GSE116008_APEX_seq_counts_all_data_polyA_hek293t.txt.gz': (f'{GEO}/GSE116nnn/GSE116008/suppl/GSE116008_APEX_seq_counts_all_data_polyA_hek293t.txt.gz',
                                                                'APEX-seq counts, poly(A), HEK293T', 'release_asset (GEO)'),
    'GSE176439_NuclearBody_APEX_seq_gene_counts.txt.gz': (f'{GEO}/GSE176nnn/GSE176439/suppl/GSE176439_NuclearBody_APEX_seq_gene_counts.txt.gz',
                                                          'nuclear-body APEX2-seq gene counts', 'release_asset (GEO)'),
    'GSE81553_SON_TSA-Seq_Condition2.bw': (f'{GEO}/GSE81nnn/GSE81553/suppl/GSE81553_SON_TSA-Seq_Condition2.bw', 'SON TSA-seq, K562, condition 2', 'not kept: served by GEO at this URL (open access); MD5 and SHA-256 identify it'),
    'GSE81553_SON_TSA-Seq_Condition2.Replicate.bw': (f'{GEO}/GSE81nnn/GSE81553/suppl/GSE81553_SON_TSA-Seq_Condition2.Replicate.bw', 'SON TSA-seq, K562, condition 2, replicate', 'not kept: served by GEO at this URL (open access); MD5 and SHA-256 identify it'),
    'GSE81553_LaminB_TSA-Seq_Condition2.bw': (f'{GEO}/GSE81nnn/GSE81553/suppl/GSE81553_LaminB_TSA-Seq_Condition2.bw', 'lamin B TSA-seq, K562, condition 2', 'not kept: served by GEO at this URL (open access); MD5 and SHA-256 identify it'),
    'GSE81553_LaminB_TSA-Seq_Conditoin2_Replicate.bw': (f'{GEO}/GSE81nnn/GSE81553/suppl/GSE81553_LaminB_TSA-Seq_Conditoin2_Replicate.bw', 'lamin B TSA-seq, K562, condition 2, replicate (GEO\'s spelling)', 'not kept: served by GEO at this URL (open access); MD5 and SHA-256 identify it'),
    'GSE63525_GM12878_subcompartments.bed.gz': (f'{GEO}/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_subcompartments.bed.gz', 'GM12878 Hi-C subcompartments (hg19)', 'release_asset (GEO)'),
    'gencode.v19.annotation.gtf.gz': ('https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz', 'GENCODE release 19 (GRCh37) gene annotation', 'release_asset (GENCODE, open)'),
    'ENCFF611MXW.tsv': (ENCODE.format('ENCFF611MXW'), 'K562 total RNA-seq gene quantifications (ENCSR109IQO), replicate 1', 'release_asset (ENCODE, open)'),
    'ENCFF384BFE.tsv': (ENCODE.format('ENCFF384BFE'), 'K562 total RNA-seq gene quantifications (ENCSR109IQO), replicate 2', 'release_asset (ENCODE, open)'),
    'ENCFF723ICA.tsv': (ENCODE.format('ENCFF723ICA'), 'GM12878 total RNA-seq gene quantifications (ENCSR000AEE), replicate 1', 'release_asset (ENCODE, open)'),
    'ENCFF362RMV.tsv': (ENCODE.format('ENCFF362RMV'), 'GM12878 total RNA-seq gene quantifications (ENCSR000AEE), replicate 2', 'release_asset (ENCODE, open)'),
}
HCM_CITE = 'Go et al., A proximity-dependent biotinylation map of a human cell, Nature 595:120-124 (2021), doi:10.1038/s41586-021-03592-2'
OIP_CITE = 'Hein et al., Global organelle profiling reveals subcellular localization and remodeling at proteome scale, Cell 188:1137-1155 (2025), doi:10.1016/j.cell.2024.11.028'
CHEP_CITE = 'Kustatscher et al., Proteomics of a fuzzy organelle: interphase chromatin, EMBO J 33:648-664 (2014), doi:10.1002/embj.201387614'
BARUTCU_CITE = 'Barutcu et al., Systematic mapping of nuclear domain-associated transcripts reveals speckles and lamina as hubs of functionally distinct retained introns, Mol Cell 82:1035-1052 (2022), doi:10.1016/j.molcel.2021.12.010'
FAZAL_CITE = 'Fazal et al., Atlas of subcellular RNA localization revealed by APEX-seq, Cell 178:473-490 (2019), doi:10.1016/j.cell.2019.05.027'
CHEN_CITE = 'Chen et al., Mapping 3D genome organization relative to nuclear compartments using TSA-Seq as a cytological ruler, J Cell Biol 217:4025-4048 (2018), doi:10.1083/jcb.201807108'
RAO_CITE = 'Rao et al., A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping, Cell 159:1665-1680 (2014), doi:10.1016/j.cell.2014.11.021'
ENCODE_CITE = 'ENCODE Project Consortium, Nature 583:699-710 (2020), doi:10.1038/s41586-020-2493-4'

APEX = {'NLS': 'Nucleus', 'NIK': 'Nucleolus', 'LMA': 'Nuclear lamina', 'NUCP': 'Nuclear pore', 'NES': 'Cytosol', 'ERM': 'ER membrane',
        'OMM': 'Outer mitochondrial membrane', 'MITO': 'Mitochondrial matrix'}
SUBNAME = {'A1': 'A1 (speckle-associated)', 'A2': 'A2 (active)', 'B1': 'B1 (Polycomb)', 'B2': 'B2 (perinucleolar, pericentromeric)',
           'B3': 'B3 (lamina-associated)', 'B4': 'B4 (chromosome 19 KZFP clusters)'}
CHROMS = [f'chr{i}' for i in range(1, 23)] + ['chrX']
BIN = 20000
TSA = {'son': ['GSE81553_SON_TSA-Seq_Condition2.bw', 'GSE81553_SON_TSA-Seq_Condition2.Replicate.bw'],
       'lmnb': ['GSE81553_LaminB_TSA-Seq_Condition2.bw', 'GSE81553_LaminB_TSA-Seq_Conditoin2_Replicate.bw']}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def strip(g) -> str:
    return str(g).split('.')[0]


class Entries:
    """Reviewed entries and the symbols that name them (UniProt primary gene names)."""
    def __init__(self, path: Path):
        e = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
        self.reviewed = set(e['Entry'])
        self.gene = dict(zip(e['Entry'], e['Gene Names (primary)']))
        self.by_symbol = defaultdict(set)
        for a, g in zip(e['Entry'], e['Gene Names (primary)']):
            for s in g.split(';'):
                if s.strip():
                    self.by_symbol[s.strip()].add(a)

    def resolve(self, ids, symbol: str = '') -> tuple[str, str]:
        """(entry, how): the first reviewed accession of ids, else the one entry symbol names; ('', '') otherwise."""
        for a in ids:
            a = a.strip().split('-')[0]
            if a in self.reviewed:
                return a, 'accession'
        hits = self.by_symbol.get(symbol, set())
        return (next(iter(hits)), 'symbol') if len(hits) == 1 else ('', '')


def first_per_entry(t: pd.DataFrame) -> pd.DataFrame:
    """One row per entry: the first (in file order) matched by accession, else the first matched by gene symbol. A protein
    group without a reviewed accession (a fragment, say) can come before the group that names the entry."""
    t = t.assign(_row=np.arange(len(t)), _symbol=(t['mapped_by'] != 'accession'))
    t = t[t['canonical_UniProtAC'] != ''].sort_values(['_symbol', '_row']).drop_duplicates('canonical_UniProtAC').sort_values('_row')
    return t.drop(columns=['_row', '_symbol'])


def protein_maps(dl: Path, ent: Entries) -> tuple[dict[str, pd.DataFrame], dict]:
    out, counts = {}, {}
    # Human Cell Map: the localised preys
    p = pd.read_csv(dl / 'preys-v1.txt', sep='\t', dtype=str, keep_default_na=False)
    got = [ent.resolve([u], s) for u, s in zip(p['UniProt'], p['symbol'])]
    p['canonical_UniProtAC'], p['mapped_by'] = [g[0] for g in got], [g[1] for g in got]
    counts['hcm'] = {'preys': len(p), 'not_mapped': int((p['canonical_UniProtAC'] == '').sum())}
    p = first_per_entry(p)
    safe = p['SAFE localization'].where(~p['SAFE localization'].isin(['no prediction', '-']), '')
    out['cell_map_hcm_go2021.tsv'] = pd.DataFrame({'canonical_UniProtAC': p['canonical_UniProtAC'], 'gene': p['canonical_UniProtAC'].map(ent.gene),
                                                   'hcm_symbol': p['symbol'], 'nmf_compartment': p['MMF localization'], 'safe_domain': safe,
                                                   'mapped_by': p['mapped_by']})
    counts['hcm'].update(staged=len(p), with_safe_domain=int((safe != '').sum()))
    # organelle IP
    t = pd.read_csv(dl / 'cluster_annotation_Dec6.csv', dtype=str, keep_default_na=False)
    got = [ent.resolve(ids.split(';'), g) for ids, g in zip(t['Majority protein IDs'], t['Gene_name_canonical'])]
    t['canonical_UniProtAC'], t['mapped_by'] = [g[0] for g in got], [g[1] for g in got]
    counts['organelle_ip'] = {'protein_groups': len(t), 'not_mapped': int((t['canonical_UniProtAC'] == '').sum())}
    t = first_per_entry(t)
    out['cell_map_organelle_ip_hein2025.tsv'] = pd.DataFrame({'canonical_UniProtAC': t['canonical_UniProtAC'], 'gene': t['canonical_UniProtAC'].map(ent.gene),
                                                              'oip_gene': t['Gene_name_canonical'], 'compartment': t['cluster_annotation'], 'mapped_by': t['mapped_by']})
    counts['organelle_ip'].update(staged=len(t), classified=int((t['cluster_annotation'] != 'unclassified').sum()))
    # chromatin probability
    with zipfile.ZipFile(dl / 'supp.zip') as zf:
        c = pd.read_excel(zf.open('embj0033-0648-sd8.xlsx'), sheet_name=0, dtype={'Majority protein IDs': str})
    c['canonical_UniProtAC'] = [ent.resolve(str(ids).split(';'))[0] for ids in c['Majority protein IDs']]
    counts['chromatin'] = {'protein_groups': len(c), 'not_mapped': int((c['canonical_UniProtAC'] == '').sum())}
    c = c[c['canonical_UniProtAC'] != ''].sort_values('Interphase chromatin probability', ascending=False, kind='mergesort').drop_duplicates('canonical_UniProtAC')
    c = c.sort_values('canonical_UniProtAC')
    out['cell_map_chromatin_kustatscher2014.tsv'] = pd.DataFrame({'canonical_UniProtAC': c['canonical_UniProtAC'], 'gene': c['canonical_UniProtAC'].map(ent.gene),
                                                                  'interphase_chromatin_probability': c['Interphase chromatin probability'].astype(float)})
    counts['chromatin']['staged'] = len(c)
    # Barutcu APEX2-MS: speckles, lamina, PML bodies
    b = pd.read_excel(dl / '1-s2.0-S1097276521010728-mmc2.xlsx', header=1)
    b['canonical_UniProtAC'] = [ent.resolve([], g)[0] for g in b['Prey Gene'].astype(str)]
    counts['nuclear_bodies'] = {'preys': int(b['Prey Gene'].nunique()), 'symbols_not_mapped': int(b.loc[b['canonical_UniProtAC'] == '', 'Prey Gene'].nunique())}
    b = b[b['canonical_UniProtAC'] != '']
    sure = b[b['AvgP'] >= 0.9]
    by = lambda baits: sure[sure['Bait'].isin(baits)].groupby('canonical_UniProtAC')['Bait'].nunique()  # noqa: E731
    speck = by(['SRSF1', 'SRSF7', 'RNPS1'])
    sets = {'Nuclear speckles': set(speck[speck >= 2].index), 'Nuclear lamina': set(by(['LMNA', 'LMNB1']).index), 'PML bodies': set(by(['SP100']).index)}
    acc = sorted(set(b['canonical_UniProtAC']))
    out['cell_map_nuclear_bodies_barutcu2022.tsv'] = pd.DataFrame({'canonical_UniProtAC': acc, 'gene': [ent.gene[a] for a in acc],
                                                                   'compartments': ['; '.join(sorted(k for k, v in sets.items() if a in v)) or 'none' for a in acc]})
    counts['nuclear_bodies'].update(staged=len(acc), **{k: len(v) for k, v in sets.items()})
    return out, counts


def encode_tpm(dl: Path, files) -> pd.Series:
    t = pd.concat([pd.read_csv(dl / f, sep='\t', usecols=['gene_id', 'TPM']).set_index('gene_id')['TPM'] for f in files], axis=1).mean(axis=1)
    t.index = [strip(i) for i in t.index]
    return t.groupby(level=0).max()


def gene_maps(dl: Path, nucleic_acid: Path) -> tuple[dict[str, pd.DataFrame], dict]:
    import pyBigWig
    names = {}
    # eCLIP target gene names first (they name genes for the nuclear-body counts, which carry symbols only)
    for src, key in (('encode', 'encode_experiment'), ('encore', 'encore_experiment')):
        tg = pd.read_csv(nucleic_acid / f'{src}_eclip_target_genes.tsv', sep='\t', dtype=str, keep_default_na=False)
        for g, n in zip(tg['gene_id'], tg['gene_name']):
            names.setdefault(strip(g), n)
    genes = {}

    def put(g, key, cats, expr, extra=None):
        x = genes.setdefault(g, {'cats': {}, 'expr': {}, 'extra': {}})
        x['cats'][key] = sorted(set(x['cats'].get(key, [])) | set(cats))
        x['expr'][key] = max(expr, x['expr'].get(key, -1.0))
        if extra:
            x['extra'][key] = extra
    # APEX-seq
    de = pd.read_csv(dl / 'GSE116008_APEX_seq_deseq2_all_data_polyA_hek293t.txt.gz', sep='\t').set_index('Ensembl_Gene')
    cnt = pd.read_csv(dl / 'GSE116008_APEX_seq_counts_all_data_polyA_hek293t.txt.gz', sep='\t', index_col=0).apply(pd.to_numeric, errors='coerce').fillna(0)
    cpm = cnt / cnt.sum() * 1e6
    expr = cpm[[c for c in cpm.columns if re.search(r'-C\d$', c)]].mean(axis=1)
    for ens, row in de.iterrows():
        e = float(expr.get(ens, 0.0))
        if e < 1:
            continue
        if isinstance(row['Common_Gene'], str):
            names.setdefault(strip(ens), row['Common_Gene'])
        put(strip(ens), 'apex', [n for code, n in APEX.items() if row[code + '.18C.P'] < 0.05 and row[code + '.18C'] > 0], e)
    # symbols -> one Ensembl gene (GENCODE v19, APEX-seq's table, the eCLIP target genes)
    gtf = []
    with gzip.open(dl / 'gencode.v19.annotation.gtf.gz', 'rt') as fh:
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
    glen = {g: b - a for _, a, b, g, _ in gtf}
    # nuclear-body APEX2-seq
    c = pd.read_csv(dl / 'GSE176439_NuclearBody_APEX_seq_gene_counts.txt.gz', sep='\t', index_col=0, low_memory=False).apply(pd.to_numeric, errors='coerce').fillna(0)
    c.columns = [x.replace('_Rep', '-rep').replace('_rep', '-rep').replace('Rep', 'rep') for x in c.columns]
    bcpm = c / c.sum() * 1e6
    bexpr = bcpm[[x for x in bcpm.columns if 'input' in x.lower()]].mean(axis=1)
    col = lambda name, rep: bcpm[f'{name}-rep{rep}']  # noqa: E731
    enr = lambda bait, rep: np.log2((col(bait, rep) + 0.5) / (col(bait + 'input', rep) + 0.5))  # noqa: E731
    over = {b: (enr(b, 1) - enr('NLS', 1) >= 1) & (enr(b, 2) - enr('NLS', 2) >= 1)
            for b in ('SRSF1', 'SRSF7', 'RNPS1', 'LMNA', 'FBL', 'PML', 'SP100', 'WRAP53', 'SMN2', 'NPAT', 'SAM68')}
    nb = {'Nuclear speckles': (over['SRSF1'].astype(int) + over['SRSF7'] + over['RNPS1']) >= 2, 'Nuclear lamina': over['LMNA'], 'Nucleolus': over['FBL'],
          'PML bodies': over['PML'] & over['SP100'], 'Cajal bodies': over['WRAP53'], 'Gems': over['SMN2'], 'Histone locus bodies': over['NPAT'], 'Sam68 bodies': over['SAM68']}
    n_nb_symbols = 0
    for sym, e in bexpr.items():
        g = by_symbol.get(str(sym), set())
        if len(g) != 1 or e < 1:
            continue
        put(next(iter(g)), 'nb', [k for k, v in nb.items() if v.get(sym, False)], float(e))
        n_nb_symbols += 1
    # TSA-seq, genes expressed in K562
    bw = {k: [pyBigWig.open(str(dl / f)) for f in fs] for k, fs in TSA.items()}

    def mean(k, chrom, a, b):
        v = [x for x in (w.stats(chrom, a, b, type='mean')[0] for w in bw[k]) if x is not None]
        return float(np.mean(v)) if v else np.nan
    edges = {}
    for k in TSA:
        vals = []
        for chrom in CHROMS:
            nbins = bw[k][0].chroms()[chrom] // BIN
            reps = [np.array([np.nan if x is None else x for x in w.stats(chrom, 0, nbins * BIN, type='mean', nBins=nbins)], dtype=float) for w in bw[k]]
            vals.append(np.nanmean(np.vstack(reps), axis=0))
        v = np.concatenate(vals)
        edges[k] = float(np.quantile(v[np.isfinite(v)], 0.9))
    tpm = encode_tpm(dl, ('ENCFF611MXW.tsv', 'ENCFF384BFE.tsv'))
    for chrom, a, b, g, n in gtf:
        e = float(tpm.get(g, 0.0))
        if chrom not in CHROMS or e < 1:
            continue
        if b - a < BIN:
            m = (a + b) // 2
            a, b = max(0, m - BIN // 2), m + BIN // 2
        s, lm = mean('son', chrom, a, b), mean('lmnb', chrom, a, b)
        if not (np.isfinite(s) and np.isfinite(lm)):
            continue
        put(g, 'tsa', (['Near nuclear speckles'] if s >= edges['son'] else []) + (['Near the nuclear lamina'] if lm >= edges['lmnb'] else []), e,
            {'son_tsa': s, 'laminb_tsa': lm})
    # Hi-C subcompartments, genes expressed in GM12878
    sub = pd.read_csv(dl / 'GSE63525_GM12878_subcompartments.bed.gz', sep='\t', header=None, usecols=[0, 1, 2, 3], names=['chrom', 'a', 'b', 'sub'])
    sub = sub[sub['sub'].isin(list(SUBNAME))]
    by_chrom = {ch: (g['a'].to_numpy(), g['b'].to_numpy(), g['sub'].to_numpy()) for ch, g in sub.sort_values(['chrom', 'a']).groupby('chrom')}
    gm = encode_tpm(dl, ('ENCFF723ICA.tsv', 'ENCFF362RMV.tsv'))
    for chrom, a, b, g, n in gtf:
        e = float(gm.get(g, 0.0))
        if chrom not in by_chrom or e < 1:
            continue
        starts, ends, lab = by_chrom[chrom]
        m = (a + b) // 2
        i = np.searchsorted(starts, m, side='right') - 1
        if i >= 0 and m < ends[i]:
            put(g, 'sub', [SUBNAME[lab[i]]], e)
    out, counts = {}, {'nuclear_body_symbols_placed': n_nb_symbols, 'tsa_edges_top_tenth': edges}
    spec = {'apex': ('rna_map_apexseq_fazal2019.tsv', 'expression_cpm'), 'nb': ('rna_map_nuclear_bodies_barutcu2022.tsv', 'expression_cpm'),
            'tsa': ('dna_map_tsaseq_chen2018.tsv', 'expression_tpm'), 'sub': ('dna_map_hic_rao2014.tsv', 'expression_tpm')}
    for key, (name, ecol) in spec.items():
        gs = sorted(g for g, x in genes.items() if key in x['cats'])
        t = pd.DataFrame({'gene_id': gs, 'gene_name': [names.get(g, '') for g in gs],
                          'gene_length': [glen.get(g) for g in gs], ecol: [genes[g]['expr'][key] for g in gs]})
        if key == 'tsa':
            t['son_tsa'] = [genes[g]['extra'][key]['son_tsa'] for g in gs]
            t['laminb_tsa'] = [genes[g]['extra'][key]['laminb_tsa'] for g in gs]
        t['compartments'] = ['; '.join(genes[g]['cats'][key]) for g in gs]
        t['gene_length'] = t['gene_length'].astype('Int64')
        out[name] = t
        cats = defaultdict(int)
        for g in gs:
            for x in genes[g]['cats'][key]:
                cats[x] += 1
        counts[name] = {'genes_placed': len(gs), 'by_compartment': dict(sorted(cats.items()))}
    return out, counts


def write_table(df: pd.DataFrame, path: Path) -> None:
    for c in df.columns:
        if df[c].dtype == object and df[c].astype(str).str.contains('[\t\n\r]', regex=True).any():
            raise SystemExit(f'{path.name}: column {c} holds a tab or line break')
    df.to_csv(path, sep='\t', index=False, lineterminator='\n')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--downloads', required=True, help='folder with the files of DOWNLOADS')
    ap.add_argument('--entries', default='inputs/reference/uniprot_human_reviewed_2026_03_entries.tsv')
    ap.add_argument('--nucleic-acid', default='inputs/nucleic_acid', help='the staged eCLIP target genes (their gene names)')
    ap.add_argument('--manifest', default='inputs/MANIFEST.tsv')
    ap.add_argument('--raw-files', default='inputs/layers/raw_files.tsv')
    ap.add_argument('--layers', default='inputs/layers', help='where --install copies the tables')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--install', action='store_true', help='copy the tables into --layers and put their rows in --manifest and --raw-files')
    args = ap.parse_args()
    dl, out = Path(args.downloads), Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    missing = [f for f in DOWNLOADS if not (dl / f).exists()]
    if missing:
        raise SystemExit(f'missing downloads in {dl}: {missing}')
    ent = Entries(Path(args.entries))
    tables, pcounts = protein_maps(dl, ent)
    gtables, gcounts = gene_maps(dl, Path(args.nucleic_acid))
    tables.update(gtables)
    for name, df in tables.items():
        write_table(df, out / name)
    inputs_md5 = {f: md5_file(dl / f) for f in DOWNLOADS}
    # the repository's own tables it read, under fixed names whatever path was given
    repo_md5 = {'reference/' + Path(args.entries).name: md5_file(args.entries)}
    repo_md5.update({f'nucleic_acid/{src}_eclip_target_genes.tsv': md5_file(Path(args.nucleic_acid) / f'{src}_eclip_target_genes.tsv') for src in ('encode', 'encore')})
    summary = {'staged_by': 'scripts/52_stage_cell_maps.py', 'protein_maps': pcounts, 'gene_maps': gcounts,
               'repository_inputs_md5': repo_md5, 'inputs_md5': inputs_md5, 'tables': {n: len(df) for n, df in tables.items()}}
    (out / 'cell_maps_summary.json').write_text(json.dumps(summary, indent=1, sort_keys=True) + '\n')
    downloaded = {f: datetime.fromtimestamp((dl / f).resolve().stat().st_mtime, timezone.utc).strftime('%Y-%m-%d') for f in DOWNLOADS}
    raw = [{'file': f, 'url': u, 'bytes': (dl / f).stat().st_size, 'md5': inputs_md5[f], 'sha256': sha256_file(dl / f),
            'downloaded': downloaded[f], 'version': what, 'kept': kept} for f, (u, what, kept) in DOWNLOADS.items()]
    rcols = Path(args.raw_files).read_text().split('\n')[0].split('\t')
    pd.DataFrame(raw)[rcols].to_csv(out / 'raw_files_rows.tsv', sep='\t', index=False, lineterminator='\n')
    day = max(downloaded.values())
    geo = lambda *s: '; '.join(DOWNLOADS[x][0] for x in s)  # noqa: E731
    how = 'scripts/52_stage_cell_maps.py (its docstring gives every rule); inputs and their MD5 in cell_maps_summary.json'
    meta = {
        'cell_map_hcm_go2021.tsv': ('Human Cell Map (HEK293 BioID): the NMF compartment and SAFE domain of each localised prey, on reviewed canonical entries',
                                    geo('preys-v1.txt'), f'Human Cell Map v1 preys, downloaded {downloaded["preys-v1.txt"]}', HCM_CITE, 'CC BY 4.0 (humancellmap.org)'),
        'cell_map_organelle_ip_hein2025.tsv': ('Organelle IP (HEK293T): the cluster annotation of each protein group, on reviewed canonical entries',
                                               geo('cluster_annotation_Dec6.csv'), f'repository commit {OIP_COMMIT[:7]}, downloaded {downloaded["cluster_annotation_Dec6.csv"]}', OIP_CITE,
                                               'BSD-3-Clause (the analysis repository)'),
        'cell_map_chromatin_kustatscher2014.tsv': ('Interphase chromatin probability of each protein (chromatin enrichment for proteomics), on reviewed canonical entries',
                                                   geo('supp.zip'), f'Table EV1 (embj0033-0648-sd8.xlsx), downloaded {downloaded["supp.zip"]}', CHEP_CITE, 'CC BY (the article\'s licence)'),
        'cell_map_nuclear_bodies_barutcu2022.tsv': ('Nuclear speckle, lamina and PML-body proteomes (APEX2-MS, Flp-In HEK293): each prey detected and the bodies it belongs to at SAINT AvgP >= 0.9',
                                                    geo('1-s2.0-S1097276521010728-mmc2.xlsx'), f'Mol Cell Table S1, downloaded {downloaded["1-s2.0-S1097276521010728-mmc2.xlsx"]}', BARUTCU_CITE,
                                                    'calls derived from Table S1 (Elsevier); the table itself is not redistributed'),
        'rna_map_apexseq_fazal2019.tsv': ('APEX-seq (HEK293T, poly(A)): each RNA with CPM >= 1 in the controls, its CPM, gene length and the compartments whose labelling it exceeds the control in (FDR < 0.05)',
                                          geo('GSE116008_APEX_seq_deseq2_all_data_polyA_hek293t.txt.gz', 'GSE116008_APEX_seq_counts_all_data_polyA_hek293t.txt.gz', 'gencode.v19.annotation.gtf.gz'),
                                          f'GEO GSE116008, GENCODE 19; downloaded {day}', f'{FAZAL_CITE}; GENCODE release 19', 'GEO (no restrictions); GENCODE (open)'),
        'rna_map_nuclear_bodies_barutcu2022.tsv': ('Nuclear-body APEX2-seq (Flp-In HEK293): each RNA with CPM >= 1 in the inputs, its CPM, gene length and the bodies whose marker labels it at least twice as much as NLS-APEX2',
                                                   geo('GSE176439_NuclearBody_APEX_seq_gene_counts.txt.gz', 'gencode.v19.annotation.gtf.gz'),
                                                   f'GEO GSE176439, GENCODE 19; downloaded {day}', f'{BARUTCU_CITE}; GENCODE release 19', 'GEO (no restrictions); GENCODE (open)'),
        'dna_map_tsaseq_chen2018.tsv': ('TSA-seq (K562): each gene expressed in K562 (TPM >= 1), its SON and lamin B TSA-seq signal and whether it lies near speckles or the lamina (top tenth of the genome\'s 20 kb bins)',
                                        geo(*TSA['son'], *TSA['lmnb'], 'ENCFF611MXW.tsv', 'ENCFF384BFE.tsv', 'gencode.v19.annotation.gtf.gz'),
                                        f'GEO GSE81553, ENCODE ENCSR109IQO, GENCODE 19; downloaded {day}', f'{CHEN_CITE}; {ENCODE_CITE}; GENCODE release 19',
                                        'GEO (no restrictions); ENCODE (open); GENCODE (open)'),
        'dna_map_hic_rao2014.tsv': ('Hi-C subcompartments (GM12878): the subcompartment at the middle of each gene expressed in GM12878 (TPM >= 1)',
                                    geo('GSE63525_GM12878_subcompartments.bed.gz', 'ENCFF723ICA.tsv', 'ENCFF362RMV.tsv', 'gencode.v19.annotation.gtf.gz'),
                                    f'GEO GSE63525, ENCODE ENCSR000AEE, GENCODE 19; downloaded {day}', f'{RAO_CITE}; {ENCODE_CITE}; GENCODE release 19',
                                    'GEO (no restrictions); ENCODE (open); GENCODE (open)'),
        'cell_maps_summary.json': ('Counts, the TSA-seq edges and the MD5 of every input of the cell-map staging', 'the downloads above', f'downloaded {day}', '', ''),
    }
    mcols = Path(args.manifest).read_text().split('\n')[0].split('\t')
    mrows = []
    for name, (content, source, version, citation, licence) in meta.items():
        p = out / name
        rows = sum(1 for _ in open(p)) - 1 if name.endswith('.tsv') else ''
        mrows.append({'file': f'layers/{name}', 'bytes': p.stat().st_size, 'md5': md5_file(p), 'rows': rows, 'redistribution': 'git',
                      'content': content, 'source': source, 'version_or_date': version, 'citation': citation, 'licence': licence,
                      'derived_from': 'layers/raw_files.tsv', 'how_made_or_checked': how, 'original_name': ''})
    pd.DataFrame(mrows)[mcols].to_csv(out / 'manifest_rows.tsv', sep='\t', index=False, lineterminator='\n')
    print(json.dumps({'protein_maps': pcounts, 'tables': summary['tables'], 'tsa_edges': gcounts['tsa_edges_top_tenth']}, indent=1))
    if args.install:
        layers = Path(args.layers)
        for name in meta:
            shutil.copyfile(out / name, layers / name)
        update_manifest(args.manifest, mrows)
        upsert_tsv_rows(args.raw_files, raw)
        refresh_manifest_row(args.manifest, args.raw_files, 'layers/raw_files.tsv')
        print(f'installed {len(meta)} files into {layers} and their rows into {args.manifest} and {args.raw_files}')


if __name__ == '__main__':
    main()

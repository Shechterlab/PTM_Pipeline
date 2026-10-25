#!/usr/bin/env python3
"""Stage DepMap Public 24Q4, PRISM Repurposing 24Q2 and the TCGA PanCancer Atlas arginine changes for the website.

Cell-line genetics, not methylation: DepMap's CRISPR screens knock out whole genes, so the gene effect and the common-
essential call are about a protein; the residue-level data are the somatic protein changes at arginines, in cell lines
(DepMap) and in tumours (TCGA, scripts/44b_fetch_tcga_arginine_changes.py).

Sources
  DepMap Public 24Q4 (Broad Institute; CC BY 4.0; doi:10.25452/figshare.plus.27993248): CRISPRGeneEffect.csv (Chronos
    gene effect; 0 = no effect, -1 = the median of common essential genes), CRISPRGeneDependency.csv (probability that a
    line depends on the gene), CRISPRInferredCommonEssentials.csv, AchillesCommonEssentialControls.csv, Model.csv,
    OmicsSomaticMutations.csv, OmicsAbsoluteCNGene.csv (PureCN absolute copy number; MTAP at 0 copies is deleted).
  PRISM Repurposing Public 24Q2 (CC BY 4.0; doi:10.6084/m9.figshare.25917643): log2 viability of each cell line after 5
    days at 2.5 uM (Extended_Primary_Data_Matrix), for the PRMT inhibitors it holds.
  TCGA PanCancer Atlas via cBioPortal (ODC Open Database License): the raw changes scripts/44b fetched.

A UniProt entry takes the Entrez gene of its primary gene name in org.Hs.eg.db (--symbol-entrez, as the disease tests
do); when the name has several Entrez genes, the one DepMap screened. A protein change counts at an arginine of the
canonical sequence only when it is a missense or stop change at an R, its transcript is the canonical isoform (DepMap:
UniprotID unlabelled or -1; TCGA: the transcript cBioPortal annotates) and the canonical sequence has an R at its
position.

Outputs (in --outdir)
  depmap_24q4_gene_effect_summary.tsv   per UniProt entry: lines screened, gene effect quantiles and histogram, lines
                                        with dependency probability >= 0.5, common essential, the 3 most dependent
                                        lineages
  depmap_24q4_arginine_changes.tsv      per cell line and change at a canonical arginine
  depmap_24q4_models.tsv                the cell lines: name, lineage, disease, which data they have, MTAP copy number
  depmap_24q4_prmt_gene_effect.tsv      gene effect and dependency probability of each PRMT in each line
  prism_24q2_prmt_inhibitors.tsv        log2 viability per line and PRMT inhibitor
  tcga_pancan_arginine_changes.tsv      per tumour sample and change at a canonical arginine (the sample ID names the
                                        patient; residue * is a stop)
  tcga_pancan_studies.tsv               the 32 studies: name, sequenced samples, citation, PubMed IDs
  depmap_summary.json                   inputs (MD5), versions, counts, what was dropped and why
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from common import parse_fasta, save_table

PRMTS = ['PRMT1', 'PRMT2', 'PRMT3', 'CARM1', 'PRMT5', 'PRMT6', 'PRMT7', 'PRMT8', 'PRMT9']
# PRISM's PRMT inhibitors (its compound list names them; the target is the enzyme each was developed against)
PRISM_PRMT = {'BRD-K70446295': ('EPZ015666', 'PRMT5'), 'BRD-K00003421': ('GSK3326595', 'PRMT5'),
              'BRD-K30421675': ('JNJ-64619178', 'PRMT5'), 'BRD-K33586289': ('EPZ020411', 'PRMT6'),
              'BRD-K86718786': ('SGC-707', 'PRMT3')}
HIST_EDGES = np.round(np.arange(-3.0, 1.0001, 0.2), 1)     # 20 bins; values outside go to the end bins
MIN_LINEAGE = 10                                            # lines a lineage needs to be ranked
DEPMAP_CHANGE = re.compile(r'^p\.R(\d+)([ACDEFGHIKLMNPQSTVWY]|Ter)$')
TCGA_CHANGE = re.compile(r'^R(\d+)([ACDEFGHIKLMNPQSTVWY*])$')


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def gene_columns(cols) -> dict[str, str]:
    """'SYMBOL (ENTREZ)' column names by Entrez ID."""
    out = {}
    for c in cols:
        m = re.match(r'^(.+) \((\d+)\)$', str(c))
        if m:
            out[m.group(2)] = c
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--depmap', required=True, help='folder with the DepMap 24Q4 files')
    ap.add_argument('--prism', required=True, help='folder with the PRISM Repurposing 24Q2 files')
    ap.add_argument('--tcga', default='', help='folder with scripts/44b output (tcga_pancan_arginine_changes_raw.tsv)')
    ap.add_argument('--fasta', default='inputs/reference/uniprot_human_reviewed_2026_03.fasta')
    ap.add_argument('--entries', default='inputs/reference/uniprot_human_reviewed_2026_03_entries.tsv')
    ap.add_argument('--symbol-entrez', default='inputs/disease/orghsegdb_3.22.0_symbol_entrez.tsv')
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    dep, prism, out = Path(args.depmap), Path(args.prism), Path(args.outdir)
    canonical = parse_fasta(args.fasta)
    dropped = Counter()

    # UniProt entry -> Entrez gene
    entries = pd.read_csv(args.entries, sep='\t', dtype=str, keep_default_na=False)
    sym = pd.read_csv(args.symbol_entrez, sep='\t', dtype=str, keep_default_na=False)
    entrez_of = sym.groupby('symbol')['entrez_id'].agg(list).to_dict()
    own = set(zip(sym.loc[sym['entrez_symbol'] == 'true', 'symbol'], sym.loc[sym['entrez_symbol'] == 'true', 'entrez_id']))

    effect = pd.read_csv(dep / 'CRISPRGeneEffect.csv', index_col=0)
    dependency = pd.read_csv(dep / 'CRISPRGeneDependency.csv', index_col=0).reindex(index=effect.index, columns=effect.columns)
    screened = gene_columns(effect.columns)
    acc_entrez, acc_gene = {}, {}
    for acc, names in zip(entries['Entry'], entries['Gene Names (primary)']):
        if acc not in canonical or not names:
            continue
        g = names.split(';')[0].strip()
        cands = entrez_of.get(g, [])
        hit = [e for e in cands if e in screened]
        if len(hit) > 1:
            hit = [e for e in hit if (g, e) in own] or hit
        acc_gene[acc] = g
        if len(hit) == 1:
            acc_entrez[acc] = hit[0]
        elif cands:
            dropped['entries_with_entrez_not_screened' if not hit else 'entries_with_several_screened_entrez'] += 1
        else:
            dropped['entries_without_entrez'] += 1

    # gene effect summaries
    models = pd.read_csv(dep / 'Model.csv', dtype=str, keep_default_na=False).set_index('ModelID')
    lineage = models['OncotreeLineage'].reindex(effect.index).fillna('')
    lin_n = lineage[lineage != ''].value_counts()
    ranked_lineages = [l for l, n in lin_n.items() if n >= MIN_LINEAGE]
    common = {l.strip() for l in (dep / 'CRISPRInferredCommonEssentials.csv').read_text().splitlines()[1:] if l.strip()}
    controls = {l.strip() for l in (dep / 'AchillesCommonEssentialControls.csv').read_text().splitlines()[1:] if l.strip()}
    rows = []
    for acc, e in sorted(acc_entrez.items()):
        col = screened[e]
        v = effect[col]
        ok = v.notna()
        x = v[ok].to_numpy()
        if not len(x):
            dropped['entries_screened_without_values'] += 1
            continue
        p = dependency[col][ok]
        hist = np.histogram(np.clip(x, HIST_EDGES[0], HIST_EDGES[-1] - 1e-9), bins=HIST_EDGES)[0]
        lin = v[ok].groupby(lineage[ok]).agg(['mean', 'size'])
        lin = lin[lin.index.isin(ranked_lineages)].sort_values('mean').head(3)
        q = np.quantile(x, [0.05, 0.25, 0.5, 0.75, 0.95])
        rows.append({'canonical_UniProtAC': acc, 'gene': acc_gene[acc], 'entrez_id': e, 'depmap_gene': col,
                     'lines': len(x), 'mean': round(float(x.mean()), 3), 'sd': round(float(x.std(ddof=1)), 3) if len(x) > 1 else '',
                     'min': round(float(x.min()), 3), 'q05': round(q[0], 3), 'q25': round(q[1], 3), 'median': round(q[2], 3),
                     'q75': round(q[3], 3), 'q95': round(q[4], 3), 'max': round(float(x.max()), 3),
                     'lines_with_probability': int(p.notna().sum()), 'dependent_lines': int((p >= 0.5).sum()),
                     'common_essential': col in common, 'essential_control': col in controls,
                     'histogram': ';'.join(map(str, hist)),
                     'most_dependent_lineages': ';'.join(f'{l}|{int(r["size"])}|{r["mean"]:.3f}' for l, r in lin.iterrows())})
    summary_table = pd.DataFrame(rows)

    # PRMT gene effect per line
    prmt = []
    for g in PRMTS:
        cols = [c for c in effect.columns if str(c).startswith(g + ' (')]
        if not cols:
            dropped[f'prmt_not_screened_{g}'] += 1
            continue
        c = cols[0]
        for mid, ge, pr in zip(effect.index, effect[c], dependency[c]):
            if pd.notna(ge):
                prmt.append({'model_id': mid, 'gene': g, 'gene_effect': round(float(ge), 4), 'dependency_probability': '' if pd.isna(pr) else round(float(pr), 4)})
    prmt = pd.DataFrame(prmt)

    # MTAP absolute copy number
    cn_head = pd.read_csv(dep / 'OmicsAbsoluteCNGene.csv', nrows=0).columns
    mtap_col = [c for c in cn_head if str(c).startswith('MTAP (')]
    mtap = pd.read_csv(dep / 'OmicsAbsoluteCNGene.csv', usecols=[cn_head[0]] + mtap_col, index_col=0)[mtap_col[0]] if mtap_col else pd.Series(dtype=float)

    # DepMap protein changes at canonical arginines
    mut = pd.read_csv(dep / 'OmicsSomaticMutations.csv', dtype=str, keep_default_na=False,
                      usecols=['ProteinChange', 'VariantInfo', 'UniprotID', 'VepSwissprot', 'ModelID', 'AF', 'Hotspot', 'LikelyLoF', 'AMPathogenicity', 'AMClass'])
    mutated_models = mut['ModelID'].nunique()
    m = mut[mut['ProteinChange'].str.match(DEPMAP_CHANGE)].copy()
    dropped['depmap_calls_at_R'] = len(m)
    iso = m['UniprotID'].str.extract(r'^([A-Z0-9]+?)(?:-(\d+))?$')
    base = iso[0].where(iso[0].notna() & (iso[0] != ''), m['VepSwissprot'].str.split('.').str[0])
    keep_iso = iso[1].isna() | (iso[1] == '1')
    dropped['depmap_calls_on_other_isoforms'] = int((~keep_iso).sum())
    m, base = m[keep_iso], base[keep_iso]
    m['canonical_UniProtAC'] = base
    dropped['depmap_calls_not_in_reference'] = int((~m['canonical_UniProtAC'].isin(canonical)).sum())
    m = m[m['canonical_UniProtAC'].isin(canonical)]
    m['position'] = m['ProteinChange'].str.extract(r'p\.R(\d+)')[0].astype(int)
    on_r = [p <= len(canonical[a]) and canonical[a][p - 1] == 'R' for a, p in zip(m['canonical_UniProtAC'], m['position'])]
    dropped['depmap_calls_not_R_on_canonical'] = int(len(m) - sum(on_r))
    m = m[on_r]
    alt = m['ProteinChange'].str.extract(r'p\.R\d+(.+)$')[0].replace('Ter', '*')
    dchanges = pd.DataFrame({'canonical_UniProtAC': m['canonical_UniProtAC'], 'position': m['position'], 'alt': alt,
                             'protein_change': 'R' + m['position'].astype(str) + alt, 'model_id': m['ModelID'],
                             'allele_fraction': m['AF'], 'hotspot': m['Hotspot'] == 'True', 'likely_lof': m['LikelyLoF'] == 'True',
                             'am_pathogenicity': m['AMPathogenicity'], 'am_class': m['AMClass']})
    dchanges = dchanges.drop_duplicates(['canonical_UniProtAC', 'position', 'alt', 'model_id']).sort_values(['canonical_UniProtAC', 'position', 'alt', 'model_id'])

    # PRISM
    pmat = pd.read_csv(prism / 'Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv', index_col=0)
    plist = pd.read_csv(prism / 'Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv', dtype=str, keep_default_na=False)
    pr = []
    for rid in pmat.index:
        key = '-'.join(str(rid).replace('BRD:', '').split('-')[:2])
        if key not in PRISM_PRMT:
            continue
        name, target = PRISM_PRMT[key]
        screen = ';'.join(sorted(set(plist.loc[plist['IDs'] == rid, 'screen'])))
        for mid, v in pmat.loc[rid].items():
            if pd.notna(v):
                pr.append({'model_id': mid, 'compound': name, 'broad_id': rid, 'screen': screen, 'target': target, 'log2_viability': round(float(v), 4)})
    pr = pd.DataFrame(pr)
    missing = sorted({v[0] for v in PRISM_PRMT.values()} - set(pr['compound']))
    if missing:
        raise SystemExit(f'PRISM matrix lacks {missing}')

    # TCGA changes at canonical arginines
    tchanges, tcga_meta = pd.DataFrame(), {}
    if args.tcga:
        t = pd.read_csv(Path(args.tcga) / 'tcga_pancan_arginine_changes_raw.tsv', sep='\t', dtype=str, keep_default_na=False)
        studies = pd.read_csv(Path(args.tcga) / 'tcga_pancan_studies.tsv', sep='\t', dtype=str, keep_default_na=False)
        dropped['tcga_changes_at_R'] = len(t)
        accs_of = {}
        for acc, e in acc_entrez.items():
            accs_of.setdefault(e, []).append(acc)
        # entries whose gene DepMap did not screen still take the Entrez gene of their name
        for acc, g in acc_gene.items():
            if acc not in acc_entrez:
                for e in entrez_of.get(g, []):
                    if (g, e) in own:
                        accs_of.setdefault(e, []).append(acc)
        t = t[t['protein_change'].str.match(TCGA_CHANGE)]
        t['position'] = t['protein_change'].str.extract(r'^R(\d+)')[0].astype(int)
        rows = []
        for r in t.itertuples(index=False):
            hits = [a for a in accs_of.get(r.entrez_id, []) if r.position <= len(canonical[a]) and canonical[a][r.position - 1] == 'R']
            if not accs_of.get(r.entrez_id):
                dropped['tcga_changes_gene_without_entry'] += 1
            elif not hits:
                dropped['tcga_changes_not_R_on_canonical'] += 1
            for a in hits:
                rows.append({'canonical_UniProtAC': a, 'position': r.position, 'alt': r.protein_change[-1], 'study_id': r.study_id,
                             'sample_id': r.sample_id, 'refseq_mrna_id': r.refseq_mrna_id, 'tumor_alt_count': r.tumor_alt_count,
                             'tumor_ref_count': r.tumor_ref_count})
        tchanges = pd.DataFrame(rows).drop_duplicates(['canonical_UniProtAC', 'position', 'alt', 'sample_id']).sort_values(['canonical_UniProtAC', 'position', 'alt', 'study_id', 'sample_id'])
        fetch = json.loads((Path(args.tcga) / 'tcga_pancan_fetch.json').read_text())
        tcga_meta = {'fetched': fetch['fetched'], 'api': fetch['api'], 'studies': len(studies),
                     'sequenced_samples': int(studies['sequenced_samples'].astype(int).sum()),
                     'raw_md5': md5(Path(args.tcga) / 'tcga_pancan_arginine_changes_raw.tsv'),
                     'studies_md5': md5(Path(args.tcga) / 'tcga_pancan_studies.tsv')}
        save_table(studies, out / 'tcga_pancan_studies.tsv')

    # models
    used = set(effect.index) | set(dchanges['model_id']) | set(pr['model_id'])
    mtab = pd.DataFrame({'model_id': sorted(used)})
    mtab['cell_line'] = mtab['model_id'].map(models['CellLineName']).fillna('')
    mtab['lineage'] = mtab['model_id'].map(models['OncotreeLineage']).fillna('')
    mtab['primary_disease'] = mtab['model_id'].map(models['OncotreePrimaryDisease']).fillna('')
    mtab['crispr'] = mtab['model_id'].isin(effect.index)
    mtab['mutation_calls'] = mtab['model_id'].isin(set(mut['ModelID']))
    mtab['prism'] = mtab['model_id'].isin(set(pr['model_id']))
    mtab['mtap_absolute_cn'] = mtab['model_id'].map(mtap).map(lambda v: '' if pd.isna(v) else round(float(v), 2))

    save_table(summary_table, out / 'depmap_24q4_gene_effect_summary.tsv')
    save_table(dchanges, out / 'depmap_24q4_arginine_changes.tsv')
    save_table(mtab, out / 'depmap_24q4_models.tsv')
    save_table(prmt, out / 'depmap_24q4_prmt_gene_effect.tsv')
    save_table(pr, out / 'prism_24q2_prmt_inhibitors.tsv')
    if len(tchanges):
        save_table(tchanges, out / 'tcga_pancan_arginine_changes.tsv')
    names = ['CRISPRGeneEffect.csv', 'CRISPRGeneDependency.csv', 'CRISPRInferredCommonEssentials.csv', 'AchillesCommonEssentialControls.csv',
             'Model.csv', 'OmicsSomaticMutations.csv', 'OmicsAbsoluteCNGene.csv']
    pnames = ['Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv', 'Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv']
    summary = {
        'generated': dt.date.today().isoformat(),
        'depmap': {'release': 'DepMap Public 24Q4', 'doi': '10.25452/figshare.plus.27993248', 'licence': 'CC BY 4.0',
                   'input_md5': {n: md5(dep / n) for n in names}},
        'prism': {'release': 'PRISM Repurposing Public 24Q2', 'doi': '10.6084/m9.figshare.25917643', 'licence': 'CC BY 4.0',
                  'dose_uM': 2.5, 'days': 5, 'input_md5': {n: md5(prism / n) for n in pnames}},
        'tcga': tcga_meta,
        'crispr_lines': int(len(effect.index)), 'crispr_genes': int(len(screened)), 'common_essential_genes': len(common),
        'entries_with_gene_effect': int(len(summary_table)), 'entries_common_essential': int(summary_table['common_essential'].sum()),
        'lineages_ranked': {l: int(lin_n[l]) for l in ranked_lineages}, 'histogram_edges': HIST_EDGES.tolist(),
        'mutation_call_lines': int(mutated_models), 'depmap_arginine_changes': int(len(dchanges)),
        'tcga_arginine_changes': int(len(tchanges)), 'prism_lines': int(pr['model_id'].nunique()),
        'prism_compounds': {n: int((pr['compound'] == n).sum()) for n in sorted(set(pr['compound']))},
        'mtap_deleted_lines': int((mtap == 0).sum()), 'mtap_lines_with_cn': int(mtap.notna().sum()),
        'dropped': dict(sorted(dropped.items())),
    }
    (out / 'depmap_summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in summary.items() if k not in ('lineages_ranked', 'histogram_edges')}, indent=2))


if __name__ == '__main__':
    main()

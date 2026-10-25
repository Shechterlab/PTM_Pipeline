#!/usr/bin/env python3
"""Export the staged DepMap, PRISM and TCGA layers (scripts/44_stage_depmap.py, inputs/layers/) for the website.

Writes docs/data/depmap/:
  pNN.json    per UniProt accession, in the 64 files of the other per-protein data (shard as scripts/36c):
                g: [lines, mean, median, q05, q25, q75, q95, min, max, dependent lines, lines with a probability,
                    common essential (0/1), histogram counts, [[lineage index, lines, mean], ...]]
                m: DepMap changes at arginines, [[position, residue, [cell line indices], hotspot, AlphaMissense x100 or
                    null, likely loss of function], ...]
                t: TCGA changes at arginines, [[position, residue, [[study index, samples], ...]], ...]
  meta.json   releases, licences, histogram edges, lineages, cell lines [id, name, lineage index, MTAP absolute copy
              number], TCGA studies, counts, and the tests of whether a methylated arginine is changed more often than
              the other arginines of its protein (--sites): pooled Fisher and within-protein Mantel-Haenszel
  prmt.json   the PRMT panel: each PRMT's gene effect across the CRISPR lines and by MTAP status (0 copies = deleted),
              and each PRISM PRMT inhibitor's log2 viability, by MTAP status and against its target's gene effect

Quantiles are [5th, 25th, 50th, 75th, 95th percentile, lines]; MTAP groups are compared with a two-sided Mann-Whitney
test, an inhibitor's viability with its target's gene effect by Spearman correlation over the lines with both.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr

from common import mantel_haenszel_ci, parse_fasta

SHARDS = 64


def shard_of(acc: str) -> int:
    h = 0
    for c in acc:
        h = (h * 31 + ord(c)) % 1000003
    return h % SHARDS


def num(v, nd=3):
    return None if pd.isna(v) or v == '' else round(float(v), nd)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--layers', default='inputs/layers')
    ap.add_argument('--out', default='docs/data/depmap')
    ap.add_argument('--sites', default='releases/v1.3.1/human_methylarginine_sites_v1.3.1.tsv', help='the release site table')
    ap.add_argument('--fasta', default='inputs/reference/uniprot_human_reviewed_2026_03.fasta')
    args = ap.parse_args()
    L, out = Path(args.layers), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    summary = json.loads((L / 'depmap_summary.json').read_text())
    read = lambda name: pd.read_csv(L / name, sep='\t', dtype=str, keep_default_na=False)  # noqa: E731
    models = read('depmap_24q4_models.tsv')
    lineages = sorted({x for x in models['lineage'] if x})
    lin_ix = {x: k for k, x in enumerate(lineages)}
    mix = {m: k for k, m in enumerate(models['model_id'])}

    prot: dict[str, dict] = {}
    g = read('depmap_24q4_gene_effect_summary.tsv')
    for r in g.itertuples(index=False):
        lins = []
        for part in filter(None, r.most_dependent_lineages.split(';')):
            name, n, mean = part.split('|')
            lins.append([lin_ix[name], int(n), round(float(mean), 3)])
        prot.setdefault(r.canonical_UniProtAC, {})['g'] = [
            int(r.lines), num(r.mean), num(r.median), num(r.q05), num(r.q25), num(r.q75), num(r.q95), num(r.min), num(r.max),
            int(r.dependent_lines), int(r.lines_with_probability), int(r.common_essential == 'True'),
            [int(x) for x in r.histogram.split(';')], lins]

    d = read('depmap_24q4_arginine_changes.tsv')
    for (acc, pos, alt), grp in d.groupby(['canonical_UniProtAC', 'position', 'alt'], sort=True):
        am = pd.to_numeric(grp['am_pathogenicity'], errors='coerce').max()
        prot.setdefault(acc, {}).setdefault('m', []).append([
            int(pos), alt, sorted(mix[m] for m in grp['model_id']), int((grp['hotspot'] == 'True').any()),
            None if pd.isna(am) else int(round(am * 100)), int((grp['likely_lof'] == 'True').any())])

    studies = pd.DataFrame()
    if (L / 'tcga_pancan_arginine_changes.tsv').exists():
        studies = read('tcga_pancan_studies.tsv')
        six = {s: k for k, s in enumerate(studies['study_id'])}
        t = read('tcga_pancan_arginine_changes.tsv')
        for (acc, pos, alt), grp in t.groupby(['canonical_UniProtAC', 'position', 'alt'], sort=True):
            per = grp.groupby('study_id')['sample_id'].nunique()
            prot.setdefault(acc, {}).setdefault('t', []).append([int(pos), alt, sorted([six[s], int(n)] for s, n in per.items())])
    for v in prot.values():
        for k in ('m', 't'):
            if k in v:
                v[k].sort(key=lambda x: (x[0], x[1]))

    shards = [dict() for _ in range(SHARDS)]
    for acc in sorted(prot):
        shards[shard_of(acc)][acc] = prot[acc]
    for k, s in enumerate(shards):
        (out / f'p{k:02d}.json').write_text(json.dumps({'prot': s}, separators=(',', ':')), encoding='utf-8')

    cn = lambda v: None if v == '' else float(v)  # noqa: E731
    # Is a methylated arginine changed (in at least one line or tumour) more often than the other arginines of the same
    # proteins? Pooled over the methylated proteins (Fisher), and within protein (Mantel-Haenszel, a protein a stratum).
    canonical = parse_fasta(args.fasta)
    sites = pd.read_csv(args.sites, sep='\t', dtype=str, keep_default_na=False, usecols=['uniprot_accession', 'position'])
    sites = sites[sites['uniprot_accession'].isin(canonical)]
    meth = set(zip(sites['uniprot_accession'], sites['position'].astype(int)))
    prots = sorted(set(sites['uniprot_accession']))
    tests = {}
    for key, name, who in (('depmap', 'depmap_24q4_arginine_changes.tsv', 'model_id'), ('tcga', 'tcga_pancan_arginine_changes.tsv', 'sample_id')):
        if not (L / name).exists():
            continue
        hit = set(map(tuple, read(name)[['canonical_UniProtAC', 'position']].drop_duplicates().assign(position=lambda x: x['position'].astype(int)).to_numpy().tolist()))
        a = b = c = d = 0
        tables = []
        for acc in prots:
            rs = [k + 1 for k, ch in enumerate(canonical[acc]) if ch == 'R']
            mm = [r for r in rs if (acc, r) in meth]
            oo = [r for r in rs if (acc, r) not in meth]
            a1, c1 = sum((acc, r) in hit for r in mm), sum((acc, r) in hit for r in oo)
            a, b, c, d = a + a1, b + len(mm) - a1, c + c1, d + len(oo) - c1
            tables.append((a1, len(mm) - a1, c1, len(oo) - c1))
        orr, p = fisher_exact([[a, b], [c, d]])
        mh, lo, hi, mp, used = mantel_haenszel_ci(tables)
        tests[key] = {'meth_changed': a, 'meth': a + b, 'other_changed': c, 'other': c + d, 'pooled_or': round(orr, 3), 'pooled_p': float(f'{p:.3g}'),
                      'mh_or': round(mh, 3), 'mh_lo': round(lo, 3), 'mh_hi': round(hi, 3), 'mh_p': float(f'{mp:.3g}'), 'proteins': used}
    meta = {
        'depmap': {'release': summary['depmap']['release'], 'doi': summary['depmap']['doi'], 'licence': summary['depmap']['licence'],
                   'crispr_lines': summary['crispr_lines'], 'crispr_genes': summary['crispr_genes'],
                   'common_essential_genes': summary['common_essential_genes'], 'mutation_call_lines': summary['mutation_call_lines'],
                   'mtap_deleted_lines': summary['mtap_deleted_lines'], 'mtap_lines_with_cn': summary['mtap_lines_with_cn']},
        'prism': {k: summary['prism'][k] for k in ('release', 'doi', 'licence', 'dose_uM', 'days')} | {'lines': summary['prism_lines']},
        'tcga': {k: summary['tcga'][k] for k in ('fetched', 'studies', 'sequenced_samples')} if summary.get('tcga') else {},
        'hist_edges': summary['histogram_edges'], 'lineages': lineages,
        'lineage_lines': summary['lineages_ranked'],
        'models': [[r.model_id, r.cell_line, lin_ix.get(r.lineage, -1), cn(r.mtap_absolute_cn)] for r in models.itertuples(index=False)],
        'studies': [[r.study_id.replace('_tcga_pan_can_atlas_2018', '').upper(), r.name.replace(' (TCGA, PanCancer Atlas)', ''), int(r.sequenced_samples)]
                    for r in studies.itertuples(index=False)] if len(studies) else [],
        'tests': tests, 'sites_table': Path(args.sites).name,
    }
    (out / 'meta.json').write_text(json.dumps(meta, separators=(',', ':')), encoding='utf-8')

    # the PRMT panel
    mtap = {r.model_id: cn(r.mtap_absolute_cn) for r in models.itertuples(index=False)}
    deleted = {m for m, v in mtap.items() if v == 0}
    intact = {m for m, v in mtap.items() if v is not None and v > 0}

    def quant(v):
        v = np.asarray([x for x in v if x is not None and not pd.isna(x)], dtype=float)
        return [round(float(x), 3) for x in np.quantile(v, [0.05, 0.25, 0.5, 0.75, 0.95])] + [int(len(v))] if len(v) else None

    def split(series: pd.Series) -> dict:
        dl, it = series[series.index.isin(deleted)], series[series.index.isin(intact)]
        return {'all': quant(series), 'mtap_deleted': quant(dl), 'mtap_intact': quant(it),
                'mtap_p': float(f'{mannwhitneyu(dl, it).pvalue:.3g}') if len(dl) and len(it) else None}

    pg = read('depmap_24q4_prmt_gene_effect.tsv')
    effect_of = {g: grp.set_index('model_id')['gene_effect'].astype(float) for g, grp in pg.groupby('gene')}
    genes = {}
    for gene, ge in effect_of.items():
        dp = pd.to_numeric(pg.loc[pg['gene'] == gene, 'dependency_probability'], errors='coerce')
        genes[gene] = split(ge) | {'dependent_lines': int((dp >= 0.5).sum())}
    pr = read('prism_24q2_prmt_inhibitors.tsv')
    prism = {}
    for comp, grp in pr.groupby('compound'):
        v = grp.set_index('model_id')['log2_viability'].astype(float)
        tgt = grp['target'].iloc[0]
        rec = {'target': tgt, 'screen': grp['screen'].iloc[0], 'broad_id': grp['broad_id'].iloc[0]} | split(v)
        if tgt in effect_of:
            j = v.index.intersection(effect_of[tgt].index)
            r = spearmanr(v[j], effect_of[tgt][j])
            rec |= {'rho': round(float(r.statistic), 3), 'rho_p': float(f'{r.pvalue:.3g}'), 'rho_lines': int(len(j))}
        rec['most_sensitive'] = [[mix[m], round(float(x), 3)] for m, x in v.sort_values().head(8).items()]
        prism[comp] = rec
    (out / 'prmt.json').write_text(json.dumps({'genes': genes, 'prism': prism}, separators=(',', ':')), encoding='utf-8')
    sizes = {p.name: p.stat().st_size for p in out.glob('*.json')}
    print(len(prot), 'proteins;', sum(sizes.values()) // 1024, 'KB in', len(sizes), 'files; largest', max(sizes.values()) // 1024, 'KB')


if __name__ == '__main__':
    main()

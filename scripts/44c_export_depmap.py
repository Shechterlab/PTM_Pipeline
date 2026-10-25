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
              number], TCGA studies, counts, and stats: the release's statistics rows of the layer (analyses
              depmap_dependency and cancer_change, scripts/release_layers.py), read from --statistics, or made with the
              same code from --sites when that release has none
  genes.json  every screened entry, for the Overview's comparisons: accessions, median gene effect x1000, common
              essential (a string of 0 and 1), dependency class (a string of indices into meta.dependency_classes),
              dependent lines, lines screened
  prmt.json   the PRMT panel: each PRMT's gene effect across the CRISPR lines and by MTAP status (0 copies = deleted),
              and each PRISM PRMT inhibitor's log2 viability, by MTAP status and against its target's gene effect

Quantiles are [5th, 25th, 50th, 75th, 95th percentile, lines]; MTAP groups are compared with a two-sided Mann-Whitney
test, an inhibitor's viability with its target's gene effect by Spearman correlation over the lines with both.
"""
from __future__ import annotations

import argparse
import json
import math
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

import release_layers
from common import parse_fasta

SHARDS = 64


def shard_of(acc: str) -> int:
    h = 0
    for c in acc:
        h = (h * 31 + ord(c)) % 1000003
    return h % SHARDS


def num(v, nd=3):
    return None if pd.isna(v) or v == '' else round(float(v), nd)


STAT_FIELDS = ['analysis', 'category', 'method', 'a', 'b', 'c', 'd', 'odds_ratio', 'odds_ratio_ci_low', 'odds_ratio_ci_high',
               'fold_enrichment', 'p', 'q', 'n_strata', 'background']


def layer_stats(args, layers: Path) -> dict:
    """The layer's statistics rows as the release has them (--statistics, or the statistics table beside --sites), or,
    for a release made before the layer was registered, made with the same code (scripts/release_layers.py) from its site
    table, the sites kept pending left out (the enrichment test then needs R, as in the build)."""
    analyses = set(release_layers.CANCER.analyses)
    sites_path = Path(args.sites)
    path = Path(args.statistics) if args.statistics else sites_path.with_name(sites_path.name.replace('human_methylarginine_sites_', 'statistics_'))
    rows, source = [], ''
    if path.exists():
        st = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
        rows = st[st['analysis'].isin(analyses)].to_dict('records')
        source = path.name
    if not rows:
        sites = pd.read_csv(sites_path, sep='\t', dtype=str, keep_default_na=False)
        frame = pd.DataFrame({'uniprot_accession': sites['uniprot_accession'], 'position': sites['position'].astype(int),
                              'counted': ~sites['flags'].str.contains('one_sided_window', regex=False)})
        with tempfile.TemporaryDirectory() as work:
            rows = [{k: ('' if v is None or (isinstance(v, float) and math.isnan(v)) else v) for k, v in r.items()}
                    for r in release_layers.cancer_stats(frame, layers.parent, parse_fasta(args.fasta), Path(work))]
        source = f'scripts/release_layers.py on {sites_path.name}'

    def val(k, v):
        if k in ('analysis', 'category', 'method', 'background'):
            return v
        if v in ('', None):
            return None
        x = float(v)
        return int(x) if k in 'abcd' or k == 'n_strata' else float(f'{x:.4g}')
    return {'rows': [[val(k, r.get(k, '')) for k in STAT_FIELDS] for r in rows], 'from': source, 'fields': STAT_FIELDS}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--layers', default='inputs/layers')
    ap.add_argument('--out', default='docs/data/depmap')
    ap.add_argument('--sites', default='releases/v1.3.1/human_methylarginine_sites_v1.3.1.tsv', help='the release site table')
    ap.add_argument('--statistics', default='', help="the release's statistics table (default: next to --sites, when it exists)")
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

    gs = g.sort_values('canonical_UniProtAC')
    classes = list(release_layers.CANCER.codes['dependency_class'])
    cls = release_layers.gene_effect(str(L.parent))['dependency_class']
    (out / 'genes.json').write_text(json.dumps({
        'acc': gs['canonical_UniProtAC'].tolist(), 'median': [int(round(float(v) * 1000)) for v in gs['median']],
        'ce': ''.join('1' if v == 'True' else '0' for v in gs['common_essential']),
        'cls': ''.join(str(classes.index(cls[a])) for a in gs['canonical_UniProtAC']),
        'dependent': gs['dependent_lines'].astype(int).tolist(), 'lines': gs['lines'].astype(int).tolist()}, separators=(',', ':')), encoding='utf-8')

    shards = [dict() for _ in range(SHARDS)]
    for acc in sorted(prot):
        shards[shard_of(acc)][acc] = prot[acc]
    for k, s in enumerate(shards):
        (out / f'p{k:02d}.json').write_text(json.dumps({'prot': s}, separators=(',', ':')), encoding='utf-8')

    cn = lambda v: None if v == '' else float(v)  # noqa: E731
    stats = layer_stats(args, L)
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
        'stats': stats['rows'], 'stats_fields': stats['fields'], 'stats_from': stats['from'], 'sites_table': Path(args.sites).name,
        'dependency_classes': [[k, v] for k, v in release_layers.CANCER.codes['dependency_class'].items()],
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

#!/usr/bin/env python3
"""Compare the analyses of several releases, each rerun with the same code (releases/VERSION/analyses/).

    python scripts/compare_release_analyses.py --analyses v1.0.0=releases/v1.0.0/analyses \
        v1.1.0=releases/v1.1.0/analyses v1.2.0=releases/v1.2.0/analyses --out releases/v1.2.0/analyses

Writes comparison_<versions>.tsv (every test of the statistics tables side by side: counts, effect, q and the call,
enriched or depleted at q <= 0.05 or not significant) and COMPARISON.md (the headline numbers of each analysis and
the tests whose call differs between releases).
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd

Q = 0.05
# The analyses compared, in the order of the report, with what a row is.
ANALYSES = [('functional_class', 'functional class'), ('panther_group', 'PANTHER group'), ('panther_class', 'PANTHER class'),
            ('motif', 'motif family at the site'), ('disorder', 'disorder context'), ('domain_context', 'domain context'),
            ('domain_class', 'nearest domain class'), ('condensate', 'CD-CODE condensate'), ('hpa_location', 'HPA location'),
            ('disease_family', 'disease family (audited)'), ('functional_subclass', 'functional subclass'),
            ('panther_subclass', 'PANTHER subclass'), ('disease_family_additional', 'disease family (additional)')]


def num(v):
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


def effect(r) -> float | None:
    """log2 of the fold enrichment (enricher) or of the odds ratio (Fisher)."""
    return num(r['log2_fold_enrichment']) if num(r['log2_fold_enrichment']) is not None else num(r['log2_odds_ratio'])


def call(e, q) -> str:
    if q is None or e is None or q > Q:
        return 'ns'
    return 'enriched' if e > 0 else 'depleted'


def headline(root: Path) -> dict:
    """The descriptive numbers of one release's analyses."""
    h = {}
    d = read(root / 'disorder' / 'disorder_context_summary.tsv').set_index('disorder_context_class')['site_count'].astype(int)
    h['sites analysed'] = int(d.sum())
    h['in a disordered region'] = d.get('disordered', 0) / d.sum()
    h['at a disorder boundary'] = d.get('disorder_boundary', 0) / d.sum()
    dc = read(root / 'domains' / 'domain_context_summary.tsv').set_index('domain_context_class')['site_count'].astype(int)
    h['inside a domain'] = dc.get('in_domain', 0) / dc.sum()
    m = read(root / 'motifs' / 'motif_family_enrichment_vs_all_proteome_arginines.tsv').set_index('motif_family')
    for fam in ('RG', 'RGG', 'GAR', 'RXR', 'CARM1_proline_rich'):
        if fam in m.index:
            h[f'{fam} at the site'] = float(m.at[fam, 'target_fraction'])
            h[f'{fam}, odds ratio vs all arginines'] = float(m.at[fam, 'odds_ratio'])
    nn = read(root / 'clustering' / 'nearest_neighbor_summary_5_10_20.tsv').set_index('distance_aa')
    for dist in ('5', '10', '20'):
        if dist in nn.index:
            obs, null = float(nn.at[dist, 'prob_neighbor_within_x']), float(nn.at[dist, 'null_mean_prob_neighbor_within_x'])
            h[f'another site within {dist} aa'] = obs
            h[f'another site within {dist} aa, obs/null'] = obs / null
    dw = read(root / 'clustering' / 'arg_methyl_multiscale_density_window_summary.tsv').set_index('window_aa')
    h['proteins with a dense window (any width)'] = int(read(root / 'clustering' / 'arg_methyl_multiscale_density_significant_regions.tsv')['canonical_UniProtAC'].nunique())
    for w in ('10', '20', '50'):
        if w in dw.index:
            h[f'proteins with a dense {w}-aa window'] = int(dw.at[w, 'significant_proteins'])
    cs = read(root / 'condensates' / 'condensate_summary.tsv').iloc[0]
    h['methylated proteins in a condensate'] = float(cs['fraction_methyl_proteins_in_any_condensate'])
    return h


def fmt(v) -> str:
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return ''
    if isinstance(v, float):
        return f'{v:.1%}' if 0 <= v <= 1 else f'{v:.2f}'
    return f'{v:,}' if isinstance(v, int) else str(v)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--analyses', nargs='+', required=True, help='VERSION=FOLDER, oldest first')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    roots = dict(a.split('=', 1) for a in args.analyses)
    versions = list(roots)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    tables = {}
    for v, root in roots.items():
        s = read(Path(root) / f'statistics_{v}.tsv')
        s = s[s['method'] != 'mantel_haenszel']
        s['key'] = s['analysis'] + '\t' + s['category'] + '\t' + s['stratum']
        tables[v] = s.drop_duplicates('key').set_index('key')
    keys = list(dict.fromkeys(k for v in versions for k in tables[v].index))
    rows = []
    for k in keys:
        analysis, category, stratum = k.split('\t')
        row = {'analysis': analysis, 'category': category, 'stratum': stratum}
        calls = []
        for v in versions:
            t = tables[v]
            if k in t.index:
                r = t.loc[k]
                e, q = effect(r), num(r['q'])
                row.update({f'a_{v}': r['a'], f'log2_effect_{v}': '' if e is None else f'{e:.3g}', f'q_{v}': r['q'], f'call_{v}': call(e, q)})
                calls.append(call(e, q))
            else:
                row.update({f'a_{v}': '', f'log2_effect_{v}': '', f'q_{v}': '', f'call_{v}': 'not tested'})
                calls.append('not tested')
        row['calls_agree'] = 'true' if len(set(c for c in calls if c != 'not tested')) <= 1 else 'false'
        rows.append(row)
    comp = pd.DataFrame(rows)
    name = '_'.join(versions)
    comp.to_csv(out / f'comparison_{name}.tsv', sep='\t', index=False)

    heads = {v: headline(Path(root)) for v, root in roots.items()}
    lines = [f'# Analyses of {", ".join(versions)} compared', '',
             'Each release\'s sites were analysed with the same code (`hpc/release_analyses.sh`; the commit is in each '
             '`analyses/README.md`), so the differences below come from the sites, not the methods. A test is called '
             f'enriched or depleted at BH q <= {Q} (the effect is the log2 fold enrichment, or the log2 odds ratio). '
             f'`comparison_{name}.tsv` has every test side by side.', '',
             '## Headline numbers', '', '| | ' + ' | '.join(versions) + ' |', '|---|' + '---:|' * len(versions)]
    for k in heads[versions[0]]:
        lines.append(f'| {k} | ' + ' | '.join(fmt(heads[v].get(k)) for v in versions) + ' |')
    lines += ['', '## Tests by analysis', '',
              '| Analysis | Tested | ' + ' | '.join(f'Enriched {v}' for v in versions) + ' | Same call in all | Spearman r of effects (first vs last) |',
              '|---|---:|' + '---:|' * len(versions) + '---:|---:|']
    for analysis, label in ANALYSES:
        g = comp[comp['analysis'] == analysis]
        if g.empty:
            continue
        enriched = [int((g[f'call_{v}'] == 'enriched').sum()) for v in versions]
        e0 = pd.to_numeric(g[f'log2_effect_{versions[0]}'], errors='coerce')
        e1 = pd.to_numeric(g[f'log2_effect_{versions[-1]}'], errors='coerce')
        both = e0.notna() & e1.notna()
        rho = e0[both].corr(e1[both], method='spearman') if both.sum() >= 3 else float('nan')
        lines.append(f'| {label} | {len(g):,} | ' + ' | '.join(f'{x:,}' for x in enriched) +
                     f" | {int((g['calls_agree'] == 'true').sum()):,} | {'' if math.isnan(rho) else f'{rho:.2f}'} |")
    for analysis, label in ANALYSES[:10]:
        g = comp[comp['analysis'] == analysis].copy()
        if g.empty:
            continue
        lines += ['', f'### {label[0].upper() + label[1:]}', '',
                  '| Category | ' + ' | '.join(f'{v}: log2 effect (q)' for v in versions) + ' |', '|---|' + '---|' * len(versions)]
        g['_sort'] = pd.to_numeric(g[f'log2_effect_{versions[-1]}'], errors='coerce').fillna(-99)
        changed = g['calls_agree'] == 'false'
        show = g[changed | g[[f'call_{v}' for v in versions]].isin(['enriched', 'depleted']).any(axis=1)] if analysis in ('condensate', 'hpa_location') else g
        for r in show.sort_values('_sort', ascending=False).to_dict('records'):
            cells = []
            for v in versions:
                if r[f'call_{v}'] == 'not tested':
                    cells.append('not tested')
                    continue
                q = num(r[f'q_{v}'])
                mark = '**' if r[f'call_{v}'] in ('enriched', 'depleted') else ''
                cells.append(f"{mark}{r[f'log2_effect_{v}']}{mark} ({'' if q is None else f'{q:.2g}'})")
            lines.append(f"| {r['category']}{' †' if r['calls_agree'] == 'false' else ''} | " + ' | '.join(cells) + ' |')
        if analysis in ('condensate', 'hpa_location'):
            lines.append(f'\n(Only the categories called in some release are listed; {len(g) - len(show)} others are not significant in any.)')
    lines += ['', '† the call (enriched, depleted, not significant) differs between releases. Bold: q <= 0.05.', '']
    (out / 'COMPARISON.md').write_text('\n'.join(lines))
    print(f'{len(comp):,} tests -> {out}')


if __name__ == '__main__':
    main()

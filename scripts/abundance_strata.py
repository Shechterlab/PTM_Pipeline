"""The protein-set tests of a release again, with the proteins split by abundance: a sensitivity test for detection.

Mass spectrometry finds abundant proteins most often, so the methylated proteins are mostly abundant ones (the median
PaxDb percentile of the methylated proteins is 74), and a protein set rich in abundant proteins can look enriched among
them whatever their methylation does. This test asks whether each enrichment holds among proteins of similar abundance.

For every protein-level test the release runs with clusterProfiler enricher (scripts/enrichment.py), it reads the test's own inputs (target, universe, sets in <root>/<folder>/enricher/<job>/; a test whose sets each have their own
background, binding_protein, has one job per set), checks that
they give the release's 2x2 counts for every set it tested (it stops otherwise), and adds one row per set:

  - strata: the deciles of PaxDb 6.1 abundance (whole-organism integrated dataset, inputs/layers/
    paxdb_human_integrated_abundance.tsv) among the test's universe, plus one stratum for the identifiers without a
    PaxDb value;
  - the Mantel-Haenszel odds ratio across them with its Robins-Breslow-Greenland 95% CI and the continuity-corrected
    two-sided CMH test (common.mantel_haenszel_ci, as the within-protein tests), BH-adjusted within the analysis;
  - the same a, b, c, d as the release's row (the counts do not change), method mantel_haenszel, and the test's
    background with ", each PaxDb abundance decile a stratum" added.

The disease families count genes (Entrez): a gene takes the highest PaxDb percentile of the reviewed entries whose
primary gene symbol org.Hs.eg.db maps to it, as scripts/14 maps them. The rows sit beside the release's own; nothing
the release reports is replaced. The residue tests, already within proteins, are not repeated. The proteins found binding
RNA or DNA in cells (binding_protein) are: four of its experiments are tested among the proteins they measured, which
controls for detection, but the other three against the whole reviewed proteome, and abundance can matter beyond
detection too.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

from common import benjamini_hochberg, entrez_by_symbol, mantel_haenszel_ci

DECILES = 10
SUFFIX = ', each PaxDb abundance decile a stratum'
# The website's log2 odds ratio for an infinite odds ratio (every set member in the target, in every informative
# stratum), and minus it for an odds ratio of 0: JSON has no infinity.
WEB_INF = 99
# analysis: (folder of its enricher job under the analyses root, or 'interactions' for the interactions root; job, or one
# job per category where it holds {category}; identifiers)
TESTS = {
    'panther_group': ('panther_class', 'panther_group', 'protein'),
    'panther_class': ('panther_class', 'panther_class', 'protein'),
    'panther_subclass': ('panther_class', 'panther_subclass', 'protein'),
    'functional_class': ('functional_class_union', 'functional_broad', 'protein'),
    'functional_subclass': ('functional_class_union', 'functional_subclass', 'protein'),
    'reactome_top_level': ('layers/reactome', 'reactome_top_level', 'protein'),
    'reactome_pathway': ('layers/reactome', 'reactome_pathway', 'protein'),
    'hpa_location': ('hpa', 'hpa_location', 'protein'),
    'condensate': ('condensates', 'condensate', 'protein'),
    'condensate_proteome': ('condensates', 'condensate_proteome', 'protein'),
    'depmap_dependency': ('layers/cancer', 'depmap_dependency', 'protein'),
    'go_binding': ('interactions', 'go_binding', 'protein'),
    'cell_hcm_compartment': ('layers/cell', 'cell_hcm_compartment', 'protein'),
    'cell_hcm_domain': ('layers/cell', 'cell_hcm_domain', 'protein'),
    'cell_organelle_ip': ('layers/cell', 'cell_organelle_ip', 'protein'),
    'cell_nuclear_body_proteome': ('layers/cell', 'cell_nuclear_body_proteome', 'protein'),
    'cell_chromatin': ('layers/cell', 'cell_chromatin', 'protein'),
    'cell_region_hcm': ('layers/cell', 'cell_region_hcm', 'protein'),
    'cell_region_organelle_ip': ('layers/cell', 'cell_region_organelle_ip', 'protein'),
    'cell_region_hpa': ('layers/cell', 'cell_region_hpa', 'protein'),
    'eclip_rna_apexseq': ('layers/cell', 'eclip_rna_apexseq', 'protein'),
    'eclip_rna_nuclear_body': ('layers/cell', 'eclip_rna_nuclear_body', 'protein'),
    'eclip_gene_tsaseq': ('layers/cell', 'eclip_gene_tsaseq', 'protein'),
    'eclip_gene_hic': ('layers/cell', 'eclip_gene_hic', 'protein'),
    'disease_family': ('disease', 'disease_family', 'gene'),
    'disease_family_additional': ('disease', 'disease_family_additional', 'gene'),
    # one clusterProfiler job per experiment, each against the proteins it measured or the reviewed proteome (scripts/37)
    'binding_protein': ('interactions', 'binding_protein_{category}', 'protein'),
}
# the job's term for a release category, where scripts/37 renames it (rna_exp -> rna_binding_experimental)
TERM_OF = {'go_binding': lambda c: c.replace('_binding_experimental', '_exp').replace('_binding', '')}


def abundance(inputs: Path) -> dict[str, pd.Series]:
    """PaxDb percentile by reviewed entry, and by Entrez gene (the highest of the entries its symbol maps to)."""
    pax = pd.read_csv(Path(inputs) / 'layers' / 'paxdb_human_integrated_abundance.tsv', sep='\t', dtype=str, keep_default_na=False)
    by_entry = pd.to_numeric(pax.set_index('canonical_UniProtAC')['percentile'], errors='coerce').dropna()
    entries = pd.read_csv(Path(inputs) / 'reference' / 'uniprot_human_reviewed_2026_03_entries.tsv', sep='\t', dtype=str, keep_default_na=False)
    entrez_of = entrez_by_symbol(Path(inputs) / 'disease' / 'orghsegdb_3.22.0_symbol_entrez.tsv')
    by_gene = {}
    for entry, symbol in zip(entries['Entry'], entries['Gene Names (primary)']):
        v = by_entry.get(entry)
        if v is None or not symbol:
            continue
        for e in entrez_of.get(symbol.split(';')[0].strip(), []):
            by_gene[e] = max(by_gene.get(e, -1.0), float(v))
    return {'protein': by_entry, 'gene': pd.Series(by_gene, dtype=float)}


def strata(universe: list[str], values: pd.Series) -> pd.Series:
    """Each identifier's abundance decile among the universe (1 the least abundant), or 'none' without a value."""
    v = values.reindex(universe)
    have = v.notna()
    out = pd.Series('none', index=universe, dtype=object)
    if have.any():
        edges = np.quantile(v[have].to_numpy(), np.linspace(0, 1, DECILES + 1))
        out[have] = np.clip(np.searchsorted(edges, v[have].to_numpy(), side='right'), 1, DECILES).astype(str)
    return out


def job_inputs(job_dir: Path, whole_universe: bool) -> tuple[list[str], set[str], dict[str, set[str]]]:
    """A job's universe, target and sets as clusterProfiler enricher tests them (scripts/enrich_sets.R): the universe is
    the job's universe file, or its identifiers in at least one set when the job does not test the whole universe."""
    lines = lambda name: [x for x in (job_dir / name).read_text().split('\n') if x]  # noqa: E731
    sets = pd.read_csv(job_dir / 'sets.tsv', sep='\t', dtype=str, keep_default_na=False)
    universe = lines('universe.txt')
    if not whole_universe:
        annotated = set(sets['gene'])
        universe = [u for u in universe if u in annotated]
    members = set(universe)
    return universe, set(lines('target.txt')) & members, sets[sets['gene'].isin(members)].groupby('term')['gene'].agg(set).to_dict()


def whole_universe_of(job_dir: Path) -> bool:
    """The job's whole_universe setting, from the jobs table scripts/enrichment.py wrote beside it."""
    for table in job_dir.parent.glob('jobs_*.tsv'):
        jobs = pd.read_csv(table, sep='\t', dtype=str, keep_default_na=False)
        hit = jobs[jobs['job'] == job_dir.name]
        if len(hit):
            return hit['whole_universe'].iloc[-1] == 'true'
    raise SystemExit(f'abundance_strata: no jobs table names {job_dir}')


def job_tests(job_dir: Path, values: pd.Series, tests: list[tuple[str, tuple[int, int, int, int]]], label: str) -> list[tuple]:
    """For each (term, the a, b, c, d its test counted) of one enricher job, the Mantel-Haenszel test over the abundance
    deciles of the job's universe: (odds ratio, CI low, CI high, p, informative strata). Stops unless the job's inputs
    give the same a, b, c, d."""
    universe, target, sets = job_inputs(job_dir, whole_universe_of(job_dir))
    code, layers = np.unique(strata(universe, values).to_numpy(dtype=str), return_inverse=True)
    k = len(code)
    pos = {u: i for i, u in enumerate(universe)}
    is_target = np.array([u in target for u in universe])
    n_target, n_all = np.bincount(layers[is_target], minlength=k), np.bincount(layers, minlength=k)
    out = []
    for term, counts in tests:
        in_set = np.zeros(len(universe), dtype=bool)
        in_set[[pos[u] for u in sets.get(term, ())]] = True
        n_both, n_set = np.bincount(layers[is_target & in_set], minlength=k), np.bincount(layers[in_set], minlength=k)
        # per stratum: a in the set and the target, b the target only, c the set only, d neither
        tables = [(int(n_both[i]), int(n_target[i] - n_both[i]), int(n_set[i] - n_both[i]),
                   int(n_all[i] - n_target[i] - n_set[i] + n_both[i])) for i in range(k)]
        got = tuple(sum(t[j] for t in tables) for j in range(4))
        if got != tuple(counts):
            raise SystemExit(f'abundance_strata: {label} {term}: the inputs in {job_dir} give a, b, c, d = {got}, the test {tuple(counts)}')
        out.append(mantel_haenszel_ci(tables))
    return out


def web_log2(odds, digits: int = 4):
    """log2 of an odds ratio for the website's JSON: +-WEB_INF where it is infinite or 0, None where there is none."""
    try:
        v = float(odds)
    except (TypeError, ValueError):
        return None
    if math.isnan(v) or v < 0:
        return None
    return -WEB_INF if v == 0 else WEB_INF if math.isinf(v) else round(math.log2(v), digits)


def bh(ps: list) -> list:
    """Benjamini-Hochberg q of p values, None where p is None."""
    return [None if math.isnan(q) else float(q) for q in benjamini_hochberg([float('nan') if p is None else p for p in ps])]


def rows(statistics: pd.DataFrame, analyses_root: Path, interactions_root: Path | None, inputs: Path) -> list[dict]:
    """The abundance-stratified rows for the protein-level tests in statistics (a release's or an analyses folder's)."""
    values = abundance(inputs)
    out = []
    for analysis, (folder, job, kind) in TESTS.items():
        tested = statistics[(statistics['analysis'] == analysis) & (statistics['method'] == 'enricher')].to_dict('records')
        if not tested:
            continue
        root = interactions_root if folder == 'interactions' else (Path(analyses_root) / folder if analyses_root else None)
        counts = [tuple(int(float(r[x])) for x in 'abcd') for r in tested]
        # the rows of each job (one job for most tests; one per category where the job names {category})
        by_job: dict[str, list[int]] = {}
        for k, r in enumerate(tested):
            by_job.setdefault(job.format(category=r['category']), []).append(k)
        res = [None] * len(tested)
        for name, idx in by_job.items():
            job_dir = Path(root) / 'enricher' / name if root else None
            if job_dir is None or not (job_dir / 'target.txt').exists():
                raise SystemExit(f'abundance_strata: the release has {analysis} rows but their inputs are not in {job_dir}')
            for k, x in zip(idx, job_tests(job_dir, values[kind], [(TERM_OF.get(analysis, str)(tested[k]['category']), counts[k]) for k in idx], analysis)):
                res[k] = x
        # BH within the analysis, over all its sets whatever job tested them
        for r, (a, b, c, d), (odds, lo, hi, p, used), q in zip(tested, counts, res, bh([x[3] for x in res])):
            out.append({'analysis': analysis, 'category': r['category'], 'stratum': 'all', 'a': a, 'b': b, 'c': c, 'd': d,
                        'odds_ratio': odds, 'log2_odds_ratio': math.log2(odds) if odds and odds > 0 and math.isfinite(odds) else None,
                        'odds_ratio_ci_low': lo, 'odds_ratio_ci_high': hi, 'p': p, 'method': 'mantel_haenszel',
                        'background': r['background'] + SUFFIX, 'n_strata': used, 'q': q})
    return out


SINCE = 'v1.4.1'
LABELS = {'panther_group': 'PANTHER class groups', 'panther_class': 'PANTHER protein classes', 'panther_subclass': 'PANTHER subclasses',
          'functional_class': 'Functional classes', 'functional_subclass': 'Functional subclasses',
          'reactome_top_level': 'Reactome top-level pathways', 'reactome_pathway': 'Reactome pathways (10 to 500 proteins)',
          'hpa_location': 'HPA locations', 'condensate': 'CD-CODE condensates, against other condensate proteins',
          'condensate_proteome': 'CD-CODE condensates, against the proteome', 'depmap_dependency': 'DepMap dependency classes',
          'go_binding': 'GO RNA and DNA binding', 'disease_family': 'Disease families (audited)', 'disease_family_additional': 'Disease families (additional)',
          'binding_protein': 'Proteins found binding RNA or DNA in cells'}


def changes_section(files: dict, version: str, previous: str) -> list[str]:
    """CHANGES.md: the abundance-stratified tests, in the release that adds them, with how many enrichments hold."""
    v = lambda x: tuple(int(k) for k in x.lstrip('v').split('.'))  # noqa: E731
    stats = files['statistics'][1]
    adjusted = stats[(stats['method'] == 'mantel_haenszel') & stats['background'].astype(str).str.endswith(SUFFIX)]
    if not (v(previous) < v(SINCE) <= v(version)) or not len(adjusted):
        return []
    num = lambda x: pd.to_numeric(x, errors='coerce')  # noqa: E731
    base = stats[stats['method'] == 'enricher'].assign(fe=lambda t: num(t['fold_enrichment']), qq=lambda t: num(t['q']))
    both = base.merge(adjusted.assign(oa=lambda t: num(t['odds_ratio']), qa=lambda t: num(t['q']))[['analysis', 'category', 'oa', 'qa']],
                      on=['analysis', 'category'])
    pax = num(files['proteins'][1]['paxdb_percentile']).dropna()
    out = ['## Protein-level tests within abundance deciles', '',
           f'Mass spectrometry finds abundant proteins most, and the methylated proteins are mostly abundant ones (median PaxDb '
           f'percentile {pax.median():.0f} for the {len(pax):,} with a value), so a protein set rich in abundant proteins can look '
           'enriched among them whatever their methylation does. Each protein-level over-representation test is run again with '
           'its background split into PaxDb 6.1 abundance deciles (method `mantel_haenszel`, the background ending '
           '"each PaxDb abundance decile a stratum"; README.md, Statistics). The rows sit beside the tests they repeat, with the '
           'same sets and counts; no other value of the release changes. Both directions are tested, so a set can be depleted '
           'among proteins of similar abundance though enriched without the strata.', '',
           '| Analysis | Sets tested | Enriched at q <= 0.05 | Of these, enriched within abundance deciles | Depleted within abundance deciles |',
           '|---|---:|---:|---:|---:|']
    for analysis, label in LABELS.items():
        g = both[both['analysis'] == analysis]
        if not len(g):
            continue
        enriched = (g['qq'] <= 0.05) & (g['fe'] > 1)
        held = enriched & (g['qa'] <= 0.05) & (g['oa'] > 1)
        depleted = (g['qa'] <= 0.05) & (g['oa'] < 1)
        out.append(f'| {label} | {len(g):,} | {int(enriched.sum()):,} | {int(held.sum()):,} | {int(depleted.sum()):,} |')
    return out + ['']

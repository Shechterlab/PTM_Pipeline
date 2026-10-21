"""Per-protein, statistics and supporting tables of a release, from the analyses of a pipeline run (called by
scripts/build_release.py after the site table is built).

Reads the analysis outputs that hpc/build_release.sh writes under the run's results root (the same folders as the
published PTM_results/) and the frozen inputs, and returns tables in the release's column format (scripts/release_codes.py):

  proteins            one row per protein of the site table
  statistics          every enrichment test the website or the paper reports, with its 2x2 counts
  disease_terms       enriched DisGeNET terms (DOSE enrichDGN: p, BH-adjusted p and q all <= 0.05)
  disease_terms_tested     every DisGeNET term DOSE tested, enriched or not
  disease_ontology    every Disease Ontology term tested
  density_windows     significant methylarginine density windows (scripts/28)
  hpa_locations       Human Protein Atlas locations of the proteins
  other_ptm_sites     other modifications on the proteins, from the PTM master table
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

import release_codes as rc
from disease_families import ALL_FAMILIES

OTHER_PTM = {'Phosphorylation': 'phospho', 'Lys acylation': 'lys_acyl', 'Ubiquitin/SUMO': 'ub_sumo', 'Lys methylation': 'lys_methyl'}
TARGET_RESIDUES = {'phospho': set('STY'), 'lys_acyl': {'K'}, 'ub_sumo': {'K'}, 'lys_methyl': {'K'}}


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, low_memory=False)


def num(value):
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(x) else x


def join(values) -> str:
    return '; '.join(str(v) for v in values)


def split(value: str, sep: str = ';') -> list[str]:
    return [p.strip() for p in str(value or '').split(sep) if p.strip()]


# ---- Statistics ------------------------------------------------------------------------------------------------------

def effect(r: dict) -> dict:
    """The effect size columns of a statistics row: clusterProfiler's fold enrichment for its tests (method enricher),
    the odds ratio for the 2x2 tests."""
    if 'fold_enrichment' in r:
        return {'fold_enrichment': num(r['fold_enrichment']), 'log2_fold_enrichment': num(r['log2_fold_enrichment']),
                'odds_ratio': None, 'log2_odds_ratio': None, 'method': r.get('method') or 'enricher'}
    return {'odds_ratio': num(r['odds_ratio']), 'log2_odds_ratio': num(r['log2_odds_ratio']), 'method': 'fisher'}


def from_fractions(df: pd.DataFrame, label: str, analysis: str, background: str, stratum: str = 'all') -> list[dict]:
    """Rows of an enrichment table that gives counts and fractions (target_count, target_fraction, background_count,
    background_fraction): the 2x2 table is recovered from the totals the fractions imply."""
    def total(count_col, frac_col):
        values = {round(float(c) / float(f)) for c, f in zip(df[count_col], df[frac_col]) if num(f) and float(f) > 0}
        if len(values) != 1:
            raise SystemExit(f'{analysis}: the fractions imply more than one total: {sorted(values)[:5]}')
        return values.pop()
    n_target, n_background = total('target_count', 'target_fraction'), total('background_count', 'background_fraction')
    rows = []
    for r in df.to_dict('records'):
        a, c = int(float(r['target_count'])), int(float(r['background_count']))
        rows.append({'analysis': analysis, 'category': r[label], 'stratum': stratum, 'a': a, 'b': n_target - a, 'c': c,
                     'd': n_background - c, **effect(r), 'p': num(r['p_value']), 'q': num(r['q_value']), 'background': background, 'n_strata': ''})
    return rows


def panther_rows(root: Path) -> list[dict]:
    rows = []
    for r in read(root / 'panther_class' / 'panther_class_enrichment.tsv').to_dict('records'):
        rows.append({'analysis': r['analysis'], 'category': r['category'], 'stratum': 'all',
                     **{k: int(r[k]) for k in 'abcd'}, **effect(r), 'p': num(r['p_value']), 'q': num(r['q_value']),
                     'background': 'other reviewed human proteins', 'n_strata': ''})
    return rows


def statistics(root: Path, proteins: pd.DataFrame) -> pd.DataFrame:
    rows = panther_rows(root)
    rows += from_fractions(read(root / 'functional_class_union' / 'arg_methyl_functional_class_broad_enrichment.tsv'), 'category',
                           'functional_class', 'other reviewed human proteins')
    rows += from_fractions(read(root / 'functional_class_union' / 'arg_methyl_functional_class_subclass_enrichment.tsv'), 'category',
                           'functional_subclass', 'other reviewed human proteins')
    rows += from_fractions(read(root / 'motif_arg_odds' / 'motif_family_enrichment_vs_methylated_protein_arginines.tsv'), 'motif_family',
                           'motif', 'other arginines of the methylated proteins')
    rows += from_fractions(read(root / 'disorder_context' / 'idr_enrichment_vs_arginines_in_methylated_proteins.tsv'), 'category',
                           'disorder', 'other arginines of the methylated proteins')
    rows += from_fractions(read(root / 'domain_context' / 'domain_context_enrichment.tsv'), 'category',
                           'domain_context', 'other arginines of the methylated proteins')
    rows += from_fractions(read(root / 'domain_context' / 'domain_class_enrichment.tsv'), 'category',
                           'domain_class', 'other arginines of the methylated proteins')
    # The disease families (scripts/14): clusterProfiler::enricher against the genes with a DisGeNET term.
    for r in read(root / 'disease' / 'disease_families.tsv').to_dict('records'):
        a, c = int(r['methylarginine_protein_count']), int(r['background_count'])
        fold = float(r['fold_enrichment'])
        rows.append({'analysis': 'disease_family' if r['audited'] == 'True' else 'disease_family_additional',
                     'category': r['disease_family'], 'stratum': 'all', 'a': a,
                     'b': int(r['methylarginine_denominator']) - a, 'c': c, 'd': int(r['background_denominator']) - c,
                     'fold_enrichment': fold, 'log2_fold_enrichment': math.log2(fold), 'p': num(r['p_value']), 'q': num(r['BH_adjusted_p']),
                     'method': 'enricher', 'background': 'other reviewed human genes (Entrez) with a DisGeNET term', 'n_strata': ''})
    rows += condensate_rows(root)
    rows += hpa_rows(proteins, root / 'hpa' / 'enricher')
    return pd.DataFrame(rows)


def condensate_rows(root: Path) -> list[dict]:
    membership = read(root / 'condensates' / 'condensate_membership.tsv')
    enrichment = read(root / 'condensates' / 'condensate_enrichment.tsv')
    summary = read(root / 'condensates' / 'condensate_summary.tsv').iloc[0]
    universe = int(summary['condensate_universe_proteins'])
    methyl_in_universe = round(float(summary['fraction_condensate_universe_proteins_methylated']) * universe)
    rows = []
    for r in enrichment.to_dict('records'):
        a, members = int(r['methyl_protein_count']), int(r['member_count'])
        rows.append({'analysis': 'condensate', 'category': r['condensate_id'], 'stratum': 'all', 'a': a,
                     'b': methyl_in_universe - a, 'c': members - a, 'd': universe - members - (methyl_in_universe - a),
                     **effect(r), 'p': num(r['p_value']), 'q': num(r['q_value']),
                     'background': 'other human proteins in CD-CODE condensates (export of 2026-03-18)', 'n_strata': ''})
    del membership
    return rows


def hpa_rows(proteins: pd.DataFrame, workdir: Path) -> list[dict]:
    """Per HPA location, methylated against other HPA-profiled proteins, over locations scored Approved or better:
    clusterProfiler::enricher (one-sided hypergeometric, fold enrichment, BH; scripts/enrichment.py)."""
    from enrichment import Job, run
    hpa = proteins.attrs['hpa_all']
    confident = hpa[hpa['reliability'] != 'Uncertain']
    universe = set(confident['canonical_UniProtAC'])
    methylated = universe & set(proteins['uniprot_accession'])
    res = run([Job('hpa_location', methylated, universe, sets=confident[['location', 'canonical_UniProtAC']], whole_universe=True)], workdir)['hpa_location']
    return [{'analysis': 'hpa_location', 'category': r['ID'], 'stratum': 'all', 'a': r['a'], 'b': r['b'], 'c': r['c'], 'd': r['d'],
             'fold_enrichment': r['FoldEnrichment'], 'log2_fold_enrichment': math.log2(r['FoldEnrichment']),
             'p': r['pvalue'], 'q': r['p.adjust'], 'method': 'enricher',
             'background': 'other proteins with an HPA 25.1 location (Approved or better)', 'n_strata': ''} for r in res.to_dict('records')]


# ---- Proteins ----------------------------------------------------------------------------------------------------------

def proteins_table(root: Path, inputs: Path, sites: pd.DataFrame, records: pd.DataFrame, entry, sequences: dict) -> tuple[pd.DataFrame, dict]:
    import hashlib
    accs = sorted(set(sites['uniprot_accession']))
    pending = sites['flags'].str.contains('one_sided_window')
    history = read(inputs / 'reference' / 'uniprot_accession_history_2026_03.tsv')
    synonyms = history[history['kind'] == 'synonym'].groupby('accession')['key'].agg(lambda k: sorted(set(k))).to_dict()
    site_of = records['placed_site'].where(records['placed_site'] != '', records['pending_site'])
    reported = records.assign(acc=site_of.str.split(':').str[0]).groupby('acc')['reported_gene'].agg(lambda g: Counter(x for x in g if x))

    functional = read(root / 'functional_class_union' / 'protein_functional_labels.tsv').drop_duplicates('canonical_UniProtAC').set_index('canonical_UniProtAC')
    panther = read(root / 'panther_class' / 'protein_panther_classes.tsv').set_index('uniprot_accession')
    odds = read(root / 'motif_arg_odds' / 'per_protein_arg_odds.tsv').drop_duplicates('canonical_UniProtAC').set_index('canonical_UniProtAC')
    density = read(root / 'methyl_arg_clustering' / 'protein_local_density_summary.tsv').drop_duplicates('canonical_UniProtAC').set_index('canonical_UniProtAC')
    best = read(root / 'methyl_arg_clustering' / 'arg_methyl_multiscale_density_best_region_per_protein.tsv').drop_duplicates('canonical_UniProtAC').set_index('canonical_UniProtAC')
    significant = read(root / 'methyl_arg_clustering' / 'arg_methyl_multiscale_density_significant_regions.tsv')
    architecture = read(root / 'domain_context' / 'domain_cluster_architecture_classes_ge3.tsv').drop_duplicates('canonical_UniProtAC').set_index('canonical_UniProtAC')
    membership = read(root / 'condensates' / 'condensate_membership.tsv')
    condensates = membership.groupby('canonical_UniProtAC')['condensate_id'].agg(lambda c: sorted(set(c))).to_dict()
    hpa_all = read(inputs / 'layers' / 'hpa_subcellular_locations.tsv')
    hpa = hpa_all[hpa_all['canonical_UniProtAC'].isin(accs)]
    family_rows = read(root / 'disease' / 'disease_family_proteins.tsv')
    families = defaultdict(dict)
    for r in family_rows.to_dict('records'):
        families[r['canonical_UniProtAC']][r['disease_family']] = int(r['contributing_term_count'])
    family_order = ALL_FAMILIES
    other = other_ptm_sites(inputs, set(accs), sequences)
    other_counts = other.groupby(['uniprot_accession', 'ptm_type']).size()

    rows = []
    for acc in accs:
        group = sites[sites['uniprot_accession'] == acc]
        seq = sequences[acc]
        gene = entry(acc, 'Gene Names (primary)')
        rep = [g for g, _ in reported.get(acc, Counter()).most_common() if g not in split(gene)]
        f = functional.loc[acc] if acc in functional.index else None
        pc = panther.loc[acc] if acc in panther.index else None   # None: a protein the analyses did not see (lower-confidence sites only)
        o = odds.loc[acc] if acc in odds.index else None
        d = density.loc[acc] if acc in density.index else None
        b = best.loc[acc] if acc in best.index else None
        a = architecture.loc[acc] if acc in architecture.index else None
        locs = hpa[hpa['canonical_UniProtAC'] == acc].sort_values(['main', 'location'], ascending=[False, True])
        fams = families.get(acc, {})
        studies = sorted({s for x in group['studies'] for s in split(x)}, key=rc.STUDY_ORDER.index)
        rows.append({
            'uniprot_accession': acc,
            'uniprot_entry_name': entry(acc, 'Entry Name'),
            'gene': gene,
            'gene_synonyms': join(s for s in synonyms.get(acc, []) if s not in split(gene)),
            'reported_genes': join(rep),
            'protein_name': group['protein_name'].iloc[0],
            'protein_length': len(seq),
            'reference_sequence_md5': hashlib.md5(seq.encode()).hexdigest().upper(),
            'n_sites': int((~pending[group.index]).sum()),
            'n_pending': int(pending[group.index].sum()),
            'studies': join(studies),
            'panther_group': pc['panther_group'] if pc is not None else '',
            'panther_class': pc['panther_class'] if pc is not None else '',
            'panther_classes': join(split(pc['panther_leaf_classes'])) if pc is not None else '',
            'panther_family': pc['panther_family'] if pc is not None else '',
            'functional_classes': join(split(f['functional_class_broad'])) if f is not None else '',
            'functional_subclasses': join(split(f['functional_class_subclass'])) if f is not None else '',
            'arginine_count': seq.count('R'),
            'arg_enrichment_log2': num(o['shrinkage_log2_enrichment']) if o is not None else None,
            'arg_enrichment_q': num(o['q_value']) if o is not None else None,
            'nearest_spacing_min': num(d['nearest_distance_min']) if d is not None else None,
            'nearest_spacing_median': num(d['nearest_distance_median']) if d is not None else None,
            **{f'max_sites_{w}aa': (num(d[f'max_methyl_sites_in_{w}aa_window']) if d is not None else None) for w in (3, 5, 10, 20, 50)},
            'best_window_aa': num(b['window_aa']) if b is not None else None,
            'best_window_start': num(b['best_window_start']) if b is not None else None,
            'best_window_end': num(b['best_window_end']) if b is not None else None,
            'best_window_sites': num(b['observed_max_methyl_sites']) if b is not None else None,
            'best_window_log2_obs_null': num(b['log2_observed_over_null_mean']) if b is not None else None,
            'best_window_q': num(b['q_value']) if b is not None else None,
            'n_dense_windows': int((significant['canonical_UniProtAC'] == acc).sum()),
            'architecture_class': a['architecture_class'] if a is not None else '',
            'architecture_clusters': join(','.join(c.split(',')) for c in split(a['best_cluster_positions'])) if a is not None else '',
            'condensates': join(condensates.get(acc, [])),
            'hpa_locations': join(locs['location']),
            'disease_families': join(x for x in family_order if x in fams),
            'disease_family_terms': join(fams[x] for x in family_order if x in fams),
            'als_ftd_family_gene': 'ALS / motor neuron / FTD' in fams,
            **{f'n_{t}': int(other_counts.get((acc, t), 0)) for t in ('phospho', 'lys_acyl', 'ub_sumo', 'lys_methyl')},
            # the symbol the disease enrichment counted for this protein (scripts/14: the UniProt primary gene, else the
            # labels' first symbol)
            '_label_symbol': gene.split(';')[0].strip() or ((split(f['substrate_genename']) or split(f['genename']) or [''])[0] if f is not None else ''),
        })
    table = pd.DataFrame(rows)
    table.attrs['hpa_all'] = hpa_all
    hpa_out = hpa.rename(columns={'canonical_UniProtAC': 'uniprot_accession'}).assign(main=hpa['main'] == '1')
    hpa_out = hpa_out.sort_values(['uniprot_accession', 'main', 'location'], ascending=[True, False, True])
    extra = {'hpa': hpa_out, 'other_ptm_sites': other,
             'density_windows': density_windows(significant, set(accs))}
    return table, extra


def density_windows(significant: pd.DataFrame, accs: set) -> pd.DataFrame:
    df = significant[significant['canonical_UniProtAC'].isin(accs)]
    out = pd.DataFrame({
        'uniprot_accession': df['canonical_UniProtAC'],
        'window_aa': df['window_aa'].astype(int),
        'start': df['best_window_start'].astype(int),
        'end': df['best_window_end'].astype(int),
        'sites': df['observed_max_methyl_sites'].astype(int),
        'arginines': df['best_window_arg_count'].astype(int),
        'null_mean_sites': df['null_mean_max_sites'].astype(float),
        'log2_obs_over_null': df['log2_observed_over_null_mean'].astype(float),
        'p': df['empirical_p_greater'].astype(float),
        'q': df['q_value'].astype(float),
        'site_positions': df['best_window_methyl_positions'].map(lambda s: join(split(s, ','))),
    })
    return out.sort_values(['uniprot_accession', 'window_aa', 'start']).reset_index(drop=True)


def other_ptm_sites(inputs: Path, accs: set, sequences: dict) -> pd.DataFrame:
    """Other modifications on the union's proteins from the PTM master table: rows on the canonical accession itself
    (isoform-numbered rows are left out), at a residue of the right kind that the reference sequence has there."""
    master = read(inputs / 'sources' / 'ptm_master_iptmnet_human.tsv')
    master = master[master['ptm_group'].isin(OTHER_PTM) & master['substrate_UniProtAC'].isin(accs)].copy()
    master['ptm_type'] = master['ptm_group'].map(OTHER_PTM)
    master['position'] = pd.to_numeric(master['position'], errors='coerce')
    master = master[master['position'].notna()].copy()
    master['position'] = master['position'].astype(int)
    master['residue'] = master['residue'].str.upper().str[:1]
    ok = [r in TARGET_RESIDUES[t] and 0 < p <= len(sequences.get(a, '')) and sequences[a][p - 1] == r
          for a, p, r, t in zip(master['substrate_UniProtAC'], master['position'], master['residue'], master['ptm_type'])]
    master = master[ok]
    out = (master.groupby(['substrate_UniProtAC', 'position', 'residue', 'ptm_type'])
           .agg(evidence=('source', lambda s: join(sorted(set(s)))), n_records=('source', 'size')).reset_index()
           .rename(columns={'substrate_UniProtAC': 'uniprot_accession'}))
    out['source'] = 'ptm_master'
    order = {t: k for k, t in enumerate(OTHER_PTM.values())}
    out = out.assign(_t=out['ptm_type'].map(order)).sort_values(['uniprot_accession', 'position', '_t']).drop(columns='_t')
    return out[['uniprot_accession', 'position', 'residue', 'ptm_type', 'source', 'evidence', 'n_records']].reset_index(drop=True)


# ---- Disease terms -------------------------------------------------------------------------------------------------------

def disease_terms_table(root: Path, proteins: pd.DataFrame) -> pd.DataFrame:
    terms = read(root / 'disease' / 'disease_terms.tsv')
    family = read(root / 'disease' / 'disease_term_families.tsv').set_index('term_id')['primary_family']
    by_symbol = defaultdict(list)
    for acc, sym in zip(proteins['uniprot_accession'], proteins['_label_symbol']):
        by_symbol[sym].append(acc)
    return pd.DataFrame({
        'term_id': terms['ID'],
        'term_name': terms['Description'],
        'family': terms['ID'].map(family).replace({'Other / unassigned': ''}),
        'n_genes_term': terms['BgRatio'].str.split('/').str[0].astype(int),
        'n_methylated': terms['Count'].astype(int),
        'fold_enrichment': terms['fold_enrichment'].astype(float),
        'p': terms['pvalue'].astype(float),
        'p_adjust': terms['p.adjust'].astype(float),
        'q': terms['qvalue'].astype(float),
        'genes': terms['geneID'].map(lambda g: join(split(g, '/'))),
        'uniprot_accessions': terms['geneID'].map(lambda g: join(sorted({a for s in split(g, '/') for a in by_symbol.get(s, [])}))),
    })


def condensates_table(root: Path, accs: set) -> pd.DataFrame:
    """Human CD-CODE condensates (export of 2026-03-18, synthetic ones left out; scripts/12) with their members."""
    m = read(root / 'condensates' / 'condensate_membership.tsv')
    rows = []
    for cid, g in m.groupby('condensate_id', sort=False):
        members = set(g['canonical_UniProtAC'])
        rows.append({'condensate_id': cid, 'condensate_name': g['condensate_name'].iloc[0], 'condensate_type': g['condensate_type'].iloc[0],
                     'n_members': len(members), 'n_members_methylated': len(members & accs)})
    return pd.DataFrame(rows)


def disease_families_table(root: Path) -> pd.DataFrame:
    f = read(root / 'disease' / 'disease_families.tsv')
    return pd.DataFrame({'family': f['disease_family'], 'audited': f['audited'] == 'True', 'n_terms': f['contributing_ontology_terms'].astype(int),
                         'n_genes': f['family_gene_count'].astype(int), 'n_methylated_genes': f['methylarginine_protein_count'].astype(int),
                         'fold_enrichment': f['fold_enrichment'].astype(float), 'q': f['BH_adjusted_p'].astype(float)})


def _accessions_of(proteins: pd.DataFrame):
    by_symbol = defaultdict(list)
    for acc, sym in zip(proteins['uniprot_accession'], proteins['_label_symbol']):
        by_symbol[sym].append(acc)
    return lambda g: join(sorted({a for s in split(g, '/') for a in by_symbol.get(s, [])}))


def _clean_name(names: pd.Series) -> pd.Series:
    # ";" is the release's list separator: a name that holds one is written with ","
    return names.str.replace('; ', ', ', regex=False).str.replace(';', ',', regex=False)


def disease_terms_tested_table(root: Path, proteins: pd.DataFrame) -> pd.DataFrame:
    """Every DisGeNET term DOSE enrichDGN tested, enriched or not (scripts/14 disease_terms_tested.tsv)."""
    t = read(root / 'disease' / 'disease_terms_tested.tsv')
    accs = _accessions_of(proteins)
    return pd.DataFrame({
        'term_id': t['ID'], 'term_name': _clean_name(t['Description']), 'family': t['family'],
        'n_genes_term': t['M'].astype(int), 'n_methylated': t['k'].astype(int), 'fold_enrichment': t['FoldEnrichment'].astype(float),
        'p': t['pvalue'].astype(float), 'p_adjust': t['p.adjust'].astype(float), 'q': t['qvalue'].astype(float), 'enriched': t['enriched'] == 'True',
        'genes': t['geneID'].map(lambda g: join(split(g, '/'))), 'uniprot_accessions': t['geneID'].map(accs),
    })


def disease_ontology_table(root: Path, proteins: pd.DataFrame) -> pd.DataFrame:
    """Every Disease Ontology term tested (scripts/14 disease_ontology.tsv)."""
    t = read(root / 'disease' / 'disease_ontology.tsv')
    accs = _accessions_of(proteins)
    return pd.DataFrame({
        'term_id': t['ID'], 'term_name': _clean_name(t['Description']), 'parents': t['parents'].map(lambda v: join(split(v, ';'))),
        'n_genes_term': t['M'].astype(int), 'n_methylated': t['k'].astype(int), 'fold_enrichment': t['FoldEnrichment'].astype(float),
        'p': t['pvalue'].astype(float), 'q': t['p.adjust'].astype(float),
        'genes': t['geneID'].map(lambda g: join(split(g, '/'))), 'uniprot_accessions': t['geneID'].map(accs),
    })

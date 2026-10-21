#!/usr/bin/env python3
"""Disease enrichment of the methylarginine proteins: DisGeNET terms, and eight audited disease families.

The term test is the one of DOSE 4.4.0 enrichDGN (Yu et al., Bioinformatics 31:608, 2015), as the published analysis
ran it in R (scripts/14_disease_enrichment.R), computed here from tables exported once from R
(scripts/00_export_disease_tables.R): DOSE's DisGeNET term-gene table and org.Hs.eg.db's symbol-to-Entrez map.

- Universe: the Entrez genes of the primary gene symbols of every reviewed protein in --labels
  (scripts/03 protein_functional_labels.tsv), restricted to genes with a DisGeNET term.
- Target: the Entrez genes of the methylarginine proteins (is_arg_methyl_protein).
- A term is tested when it has 10 to 300 universe genes and at least one target gene. p = P(X >= k) for the
  hypergeometric distribution (k target genes in the term, M universe genes in the term, N universe genes, n target
  genes with a term); Benjamini-Hochberg over the tested terms; q values as DOSE computes them (qvalue, lambda 0.05).
- Kept: p, adjusted p and q at most 0.05.

Families: each DisGeNET term is assigned to at most one family (scripts/disease_families.py); a family's genes are the
union of its terms' universe genes. Each family is tested against every reviewed gene (the universe without the
restriction to genes with a term), with Benjamini-Hochberg over the eight audited families and apart over the
additional, rule-based families (column audited False), in four ways:

- uncorrected: a two-sided Fisher exact test of target against the other genes (as scripts/34, and v1.1.0).
- corrected for expression (--expression): a Mantel-Haenszel odds ratio and test (chi-square with continuity
  correction) over strata of the gene's expression, the deciles of its HPA cell-line median nTPM (the highest of its
  reviewed proteins), and one stratum for genes without HPA RNA data. Proteins expressed highly in cell lines are the
  ones mass spectrometry finds methylated, and they are also the better studied ones.
- corrected for annotation: the same over strata of how many DisGeNET terms outside the family name the gene (none,
  then the quintiles of the others). Well-studied genes carry more disease terms of every kind, so every family gene
  set favours them; counting only the terms outside the family keeps the family's own membership out of it.
- corrected: both, over the joint strata (the default on the website from v1.1.1).

Terms: every term with 10 to 300 universe genes and at least one target gene is also tested with the corrected test
(strata of expression and of the gene's other terms), Benjamini-Hochberg over them; disease_terms_corrected.tsv lists
them all, enriched or not.

    python scripts/14_disease_enrichment.py --labels RESULTS/functional_class_union/protein_functional_labels.tsv \\
        --term-genes inputs/disease/dose_4.4.0_disgenet_term_genes.tsv \\
        --symbol-entrez inputs/disease/orghsegdb_3.22.0_symbol_entrez.tsv \\
        --expression inputs/annotation/hpa_25.1_rna_expression.tsv --outdir RESULTS/disease

Writes to --outdir:
  disease_terms.tsv             enriched terms (DOSE's columns, gene symbols in geneID)
  disease_terms_corrected.tsv   every tested term with the corrected test, and whether DOSE finds it enriched
  disease_families.tsv          the family tests: the eight audited families, then the additional ones
  disease_family_proteins.tsv   for each methylarginine protein, the families its gene is in
  disease_term_families.tsv     every DisGeNET term with its family assignment
  disease_input_summary.tsv     proteins, symbols and Entrez genes of the universe and the target
"""
from __future__ import annotations

import argparse
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2, fisher_exact, hypergeom

from disease_families import ADDITIONAL_FAMILIES, ALL_FAMILIES, AUDITED_FAMILIES, classify

MIN_SIZE, MAX_SIZE, CUTOFF = 10, 300, 0.05


def primary_gene(value: str) -> str:
    return str(value or '').split(';')[0].strip()


def is_true(value) -> bool:
    return str(value).strip().lower() in {'1', 'true', 't', 'yes', 'y'}


def bh(p: list[float]) -> list[float]:
    """R's p.adjust(p, 'BH')."""
    m = len(p)
    order = sorted(range(m), key=lambda i: p[i], reverse=True)
    out, running = [0.0] * m, 1.0
    for rank, i in enumerate(order):
        running = min(running, p[i] * m / (m - rank))
        out[i] = min(1.0, running)
    return out


def qvalues(p: list[float], lam: float = 0.05) -> list[float]:
    """qvalue(p, lambda = 0.05) as DOSE calls it: with one lambda, pi0 = #(p >= lambda) / (m (1 - lambda)), at most 1."""
    pi0 = min(1.0, sum(x >= lam for x in p) / (len(p) * (1 - lam)))
    return [pi0 * q for q in bh(p)]


def mantel_haenszel(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> tuple[float, float]:
    """Mantel-Haenszel odds ratio over strata of 2x2 tables (a: target in the set, b: target outside it, c: other genes
    in the set, d: other genes outside it) and the p of the Mantel-Haenszel chi-square with continuity correction."""
    n = a + b + c + d
    ok = n > 1
    a, b, c, d, n = (x[ok].astype(float) for x in (a, b, c, d, n))
    num, den = float(np.sum(a * d / n)), float(np.sum(b * c / n))
    expected = np.sum((a + b) * (a + c) / n)
    variance = np.sum((a + b) * (c + d) * (a + c) * (b + d) / (n * n * (n - 1)))
    p = float(chi2.sf((abs(a.sum() - expected) - 0.5) ** 2 / variance, 1)) if variance > 0 else 1.0
    return (num / den if den else math.inf), p


def quantile_bins(values: np.ndarray, k: int, missing: np.ndarray | None = None) -> np.ndarray:
    """0..k-1 by the k-quantiles of values (over the non-missing ones); k for the missing."""
    have = ~missing if missing is not None else np.ones(len(values), bool)
    edges = np.quantile(values[have], np.linspace(0, 1, k + 1)[1:-1])
    out = np.searchsorted(edges, values, side='right')
    return np.where(have, out, k)


def depth_bins(depth: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """0 for genes no other term names, then 1..5 by the quintiles of the reference counts above 0."""
    edges = np.quantile(reference[reference > 0], [0.2, 0.4, 0.6, 0.8])
    return np.where(depth > 0, 1 + np.searchsorted(edges, depth, side='right'), 0)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--labels', required=True)
    ap.add_argument('--term-genes', required=True)
    ap.add_argument('--symbol-entrez', required=True)
    ap.add_argument('--expression', required=True, help='HPA RNA expression per reviewed accession (scripts/00_stage_hpa_rna_expression.py)')
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    labels = pd.read_csv(args.labels, sep='\t', dtype=str, keep_default_na=False)
    labels['gene_symbol'] = [primary_gene(a) or primary_gene(b) for a, b in zip(labels['substrate_genename'], labels.get('genename', [''] * len(labels)))]
    labels = labels[labels['gene_symbol'] != '']
    target_rows = labels[labels['is_arg_methyl_protein'].map(is_true)]

    sym = pd.read_csv(args.symbol_entrez, sep='\t', dtype=str, keep_default_na=False)
    entrez_of = sym.groupby('symbol')['entrez_id'].agg(list).to_dict()
    name_of = dict(zip(sym.loc[sym['entrez_symbol'] == 'true', 'entrez_id'], sym.loc[sym['entrez_symbol'] == 'true', 'symbol']))

    def genes(symbols) -> list[str]:
        out, seen = [], set()
        for s in dict.fromkeys(symbols):
            for e in entrez_of.get(s, []):
                if e not in seen:
                    seen.add(e)
                    out.append(e)
        return out

    universe_list = genes(labels['gene_symbol'])
    targeted = set(genes(target_rows['gene_symbol']))
    target_list = [g for g in universe_list if g in targeted]  # in universe order, as DOSE receives them from R
    summary = pd.DataFrame([
        {'target': 'Universe', 'protein_count': labels['canonical_UniProtAC'].nunique(),
         'symbol_count': labels['gene_symbol'].nunique(), 'entrez_count': len(universe_list)},
        {'target': 'All methyl proteins', 'protein_count': target_rows['canonical_UniProtAC'].nunique(),
         'symbol_count': target_rows['gene_symbol'].nunique(), 'entrez_count': len(target_list)},
    ])
    summary.to_csv(outdir / 'disease_input_summary.tsv', sep='\t', index=False, lineterminator='\n')

    tg = pd.read_csv(args.term_genes, sep='\t', dtype=str, keep_default_na=False, quoting=3)
    term_name = dict(zip(tg['term_id'], tg['term_name']))
    term_genes_all = tg.groupby('term_id')['entrez_id'].agg(set).to_dict()
    gene_terms = tg.groupby('entrez_id')['term_id'].agg(set).to_dict()

    # ---- Terms (DOSE enricher_internal) ---------------------------------------------------------------------------
    ext = set(tg['entrez_id']) & set(universe_list)
    query = [g for g in target_list if g in gene_terms]       # target genes with a term (DOSE's qExtID2TermID)
    n, N = len(query), len(ext)
    q_terms = defaultdict(list)
    for g in query:
        for t in gene_terms[g]:
            if g in ext:
                q_terms[t].append(g)
    rows = []
    for t in sorted(set().union(*(gene_terms[g] for g in query))):
        M = len(term_genes_all[t] & ext)
        if not MIN_SIZE <= M <= MAX_SIZE:
            continue
        k = len(q_terms.get(t, []))
        mu = M * n / N
        sigma = mu * (N - n) * (N - M) / N / (N - 1)
        rows.append({'ID': t, 'Description': term_name[t], 'GeneRatio': f'{k}/{n}', 'BgRatio': f'{M}/{N}',
                     'RichFactor': k / M, 'FoldEnrichment': k / M * N / n, 'zScore': (k - mu) / math.sqrt(sigma),
                     'pvalue': float(hypergeom.sf(k - 1, N, M, n)),
                     'geneID': '/'.join(name_of.get(g, g) for g in q_terms.get(t, [])), 'Count': k})
    terms = pd.DataFrame(rows)
    terms['p.adjust'] = bh(terms['pvalue'].tolist())
    terms['qvalue'] = qvalues(terms['pvalue'].tolist())
    terms = terms[(terms['pvalue'] <= CUTOFF) & (terms['p.adjust'] <= CUTOFF) & (terms['qvalue'] <= CUTOFF)]
    terms = terms.assign(_c=-terms['Count']).sort_values(['p.adjust', '_c'], kind='stable').drop(columns='_c')
    terms['gene_ratio_num'] = terms['Count'] / n
    terms['bg_ratio_num'] = terms['BgRatio'].map(lambda r: int(r.split('/')[0]) / int(r.split('/')[1]))
    terms['fold_enrichment'] = terms['gene_ratio_num'] / terms['bg_ratio_num']
    terms = terms[['ID', 'Description', 'GeneRatio', 'BgRatio', 'RichFactor', 'FoldEnrichment', 'zScore', 'pvalue',
                   'p.adjust', 'qvalue', 'geneID', 'Count', 'gene_ratio_num', 'bg_ratio_num', 'fold_enrichment']]
    terms.to_csv(outdir / 'disease_terms.tsv', sep='\t', index=False, lineterminator='\n')

    # ---- Families (scripts/34) --------------------------------------------------------------------------------------
    universe, target = set(universe_list), set(target_list)
    background = universe - target
    audit, family_terms = [], defaultdict(list)
    for t in sorted(term_name):
        fam, status, method, reason, secondary = classify(term_name[t])
        audit.append({'term_id': t, 'term_name': term_name[t], 'primary_family': fam, 'audit_status': status,
                      'assignment_method': method, 'secondary_rule_matches': secondary, 'override_reason': reason,
                      'genes_in_reviewed_universe': len(term_genes_all[t] & universe)})
        if status == 'included':
            family_terms[fam].append(t)
    pd.DataFrame(audit).to_csv(outdir / 'disease_term_families.tsv', sep='\t', index=False, lineterminator='\n')
    # Strata for the corrected tests: expression (deciles of the HPA cell-line median nTPM, and no data) and annotation
    # depth (how many DisGeNET terms name the gene, counted outside the family or term under test).
    genes_arr = np.array(universe_list)
    is_target = np.isin(genes_arr, list(target))
    expression = pd.read_csv(args.expression, sep='\t', dtype={'uniprot_accession': str}, keep_default_na=False, na_values=[''])
    expr_of = dict(zip(expression['uniprot_accession'], expression['celline_median_ntpm']))
    gene_expr: dict[str, float] = {}
    for acc, symbol in zip(labels['canonical_UniProtAC'], labels['gene_symbol']):
        v = expr_of.get(acc)
        if v is not None and not pd.isna(v):
            for e in entrez_of.get(symbol, []):
                gene_expr[e] = max(gene_expr.get(e, -1.0), float(v))
    expr = np.array([gene_expr.get(g, np.nan) for g in universe_list])
    e_bin = quantile_bins(np.log1p(np.nan_to_num(expr)), 10, np.isnan(expr))
    fam_of_term = {t: f for f, ts in family_terms.items() for t in ts}
    gene_term_count = np.array([len(gene_terms.get(g, ())) for g in universe_list])
    gene_fam_count = [defaultdict(int) for _ in universe_list]
    for i, g in enumerate(universe_list):
        for t in gene_terms.get(g, ()):
            gene_fam_count[i][fam_of_term.get(t, '')] += 1

    def stratified(x: np.ndarray, strata: np.ndarray) -> tuple[float, float, int]:
        k = int(strata.max()) + 1
        count = lambda m: np.bincount(strata[m], minlength=k)  # noqa: E731
        odds, p = mantel_haenszel(count(x & is_target), count(~x & is_target), count(x & ~is_target), count(~x & ~is_target))
        return odds, p, int((np.bincount(strata, minlength=k) > 1).sum())

    fam_rows, fam_genes = [], {}
    for fam in ALL_FAMILIES:
        g = set().union(*(term_genes_all[t] & universe for t in family_terms[fam])) if family_terms[fam] else set()
        fam_genes[fam] = g
        a, c = len(target & g), len(background & g)
        b, d = len(target) - a, len(background) - c
        odds = (a * d) / (b * c) if b * c else ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        row = {'disease_family': fam, 'contributing_ontology_terms': len(family_terms[fam]),
               'methylarginine_protein_count': a, 'background_count': c, 'family_gene_count': a + c,
               'methylarginine_denominator': len(target), 'background_denominator': len(background),
               'odds_ratio': odds, 'p_value': float(fisher_exact([[a, b], [c, d]]).pvalue),
               'audited': fam in AUDITED_FAMILIES}
        x = np.isin(genes_arr, list(g))
        d_bin = depth_bins(gene_term_count - np.array([fc.get(fam, 0) for fc in gene_fam_count]),
                           gene_term_count - np.array([fc.get(fam, 0) for fc in gene_fam_count]))
        for name, strata in (('expression', e_bin), ('annotation', d_bin), ('corrected', e_bin * 10 + d_bin)):
            row[f'odds_ratio_{name}'], row[f'p_value_{name}'], row[f'n_strata_{name}'] = stratified(x, strata)
        fam_rows.append(row)
    for names in (AUDITED_FAMILIES, ADDITIONAL_FAMILIES):
        part = [r for r in fam_rows if r['disease_family'] in names]
        for suffix in ('', '_expression', '_annotation', '_corrected'):
            for r, q in zip(part, bh([r['p_value' + suffix] for r in part])):
                r['BH_adjusted_p' + suffix] = q
    pd.DataFrame(fam_rows).sort_values(['audited', 'p_value'], ascending=[False, True], kind='stable').to_csv(outdir / 'disease_families.tsv', sep='\t', index=False, lineterminator='\n')

    # ---- Terms, corrected ---------------------------------------------------------------------------------------------
    # Members of the term count its other terms (one fewer); everyone else keeps their count.
    K = 11 * 10
    base = e_bin * 10 + depth_bins(gene_term_count, gene_term_count)
    N0 = np.bincount(base, minlength=K).astype(float)
    Y0 = np.bincount(base, weights=is_target, minlength=K)
    index_of = {g: i for i, g in enumerate(universe_list)}
    dose_padj = dict(zip(terms['ID'], terms['p.adjust']))
    corrected = []
    for t in sorted(term_genes_all):
        idx = np.array(sorted(index_of[g] for g in term_genes_all[t] if g in index_of), dtype=int)
        if not MIN_SIZE <= len(idx) <= MAX_SIZE or not is_target[idx].any():
            continue
        ym = is_target[idx].astype(float)
        old, new = base[idx], e_bin[idx] * 10 + depth_bins(gene_term_count[idx] - 1, gene_term_count)
        N = N0 - np.bincount(old, minlength=K) + np.bincount(new, minlength=K)
        Y = Y0 - np.bincount(old, weights=ym, minlength=K) + np.bincount(new, weights=ym, minlength=K)
        In = np.bincount(new, minlength=K).astype(float)
        A = np.bincount(new, weights=ym, minlength=K)
        odds, p = mantel_haenszel(A, Y - A, In - A, N - Y - In + A)
        corrected.append({'term_id': t, 'term_name': term_name[t], 'family': fam_of_term.get(t, ''), 'genes': len(idx),
                          'methylated': int(ym.sum()), 'odds_ratio': odds, 'p_value': p,
                          'dose_enriched': t in dose_padj, 'dose_p_adjust': dose_padj.get(t, ''),
                          'methylated_genes': '/'.join(name_of.get(universe_list[i], universe_list[i]) for i in idx if is_target[i])})
    corrected = pd.DataFrame(corrected)
    corrected['BH_adjusted_p'] = bh(corrected['p_value'].tolist())
    corrected = corrected.sort_values(['p_value', 'term_id'], kind='stable')[['term_id', 'term_name', 'family', 'genes', 'methylated', 'odds_ratio', 'p_value',
                                                                           'BH_adjusted_p', 'dose_enriched', 'dose_p_adjust', 'methylated_genes']]
    corrected.to_csv(outdir / 'disease_terms_corrected.tsv', sep='\t', index=False, lineterminator='\n')
    prot = []
    for acc, s in zip(target_rows['canonical_UniProtAC'], target_rows['gene_symbol']):
        ents = set(entrez_of.get(s, []))
        for fam in ALL_FAMILIES:
            hit = sorted(ents & fam_genes[fam])
            if hit:
                prot.append({'canonical_UniProtAC': acc, 'gene_symbol': s, 'entrez_ids': ';'.join(hit), 'disease_family': fam,
                             'contributing_term_count': len({t for t in family_terms[fam] for e in hit if e in term_genes_all[t]})})
    pd.DataFrame(prot, columns=['canonical_UniProtAC', 'gene_symbol', 'entrez_ids', 'disease_family', 'contributing_term_count']).to_csv(
        outdir / 'disease_family_proteins.tsv', sep='\t', index=False, lineterminator='\n')
    sig = corrected[corrected['BH_adjusted_p'] <= CUTOFF]
    print(f'corrected tests: {len(corrected)} terms, {int((sig["odds_ratio"] > 1).sum())} enriched and {int((sig["odds_ratio"] < 1).sum())} depleted '
          f'(BH <= {CUTOFF}); {int(np.isnan(expr).sum())} universe genes without HPA RNA data')
    print(f'{len(terms)} enriched terms; target {len(target)} of {len(universe)} genes; '
          f'{len({r["canonical_UniProtAC"] for r in prot})} proteins in a family')


if __name__ == '__main__':
    main()

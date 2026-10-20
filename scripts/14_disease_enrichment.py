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
union of its terms' universe genes, and each family is tested with a two-sided Fisher exact test of target against the
other universe genes, with Benjamini-Hochberg over the eight families (as scripts/34).

    python scripts/14_disease_enrichment.py --labels RESULTS/functional_class_union/protein_functional_labels.tsv \\
        --term-genes inputs/disease/dose_4.4.0_disgenet_term_genes.tsv \\
        --symbol-entrez inputs/disease/orghsegdb_3.22.0_symbol_entrez.tsv --outdir RESULTS/disease

Writes to --outdir:
  disease_terms.tsv             enriched terms (DOSE's columns, gene symbols in geneID)
  disease_families.tsv          the eight family tests
  disease_family_proteins.tsv   for each methylarginine protein, the families its gene is in
  disease_term_families.tsv     every DisGeNET term with its family assignment
  disease_input_summary.tsv     proteins, symbols and Entrez genes of the universe and the target
"""
from __future__ import annotations

import argparse
import math
from collections import defaultdict
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact, hypergeom

from disease_families import FAMILY_RULES, classify

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


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--labels', required=True)
    ap.add_argument('--term-genes', required=True)
    ap.add_argument('--symbol-entrez', required=True)
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
    fam_rows, fam_genes = [], {}
    for fam, _ in FAMILY_RULES:
        g = set().union(*(term_genes_all[t] & universe for t in family_terms[fam])) if family_terms[fam] else set()
        fam_genes[fam] = g
        a, c = len(target & g), len(background & g)
        b, d = len(target) - a, len(background) - c
        odds = (a * d) / (b * c) if b * c else ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        fam_rows.append({'disease_family': fam, 'contributing_ontology_terms': len(family_terms[fam]),
                         'methylarginine_protein_count': a, 'background_count': c, 'family_gene_count': a + c,
                         'methylarginine_denominator': len(target), 'background_denominator': len(background),
                         'odds_ratio': odds, 'p_value': float(fisher_exact([[a, b], [c, d]]).pvalue)})
    for r, q in zip(fam_rows, bh([r['p_value'] for r in fam_rows])):
        r['BH_adjusted_p'] = q
    pd.DataFrame(fam_rows).sort_values('p_value', kind='stable').to_csv(outdir / 'disease_families.tsv', sep='\t', index=False, lineterminator='\n')
    prot = []
    for acc, s in zip(target_rows['canonical_UniProtAC'], target_rows['gene_symbol']):
        ents = set(entrez_of.get(s, []))
        for fam, _ in FAMILY_RULES:
            hit = sorted(ents & fam_genes[fam])
            if hit:
                prot.append({'canonical_UniProtAC': acc, 'gene_symbol': s, 'entrez_ids': ';'.join(hit), 'disease_family': fam,
                             'contributing_term_count': len({t for t in family_terms[fam] for e in hit if e in term_genes_all[t]})})
    pd.DataFrame(prot, columns=['canonical_UniProtAC', 'gene_symbol', 'entrez_ids', 'disease_family', 'contributing_term_count']).to_csv(
        outdir / 'disease_family_proteins.tsv', sep='\t', index=False, lineterminator='\n')
    print(f'{len(terms)} enriched terms; target {len(target)} of {len(universe)} genes; '
          f'{len({r["canonical_UniProtAC"] for r in prot})} proteins in a family')


if __name__ == '__main__':
    main()

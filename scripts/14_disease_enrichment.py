#!/usr/bin/env python3
"""Disease enrichment of the methylarginine proteins, with published tools (scripts/enrichment.py runs them in R):

- DisGeNET terms: clusterProfiler::enricher on DOSE 4.4.0's DisGeNET term-gene table (--term-genes, exported by
  scripts/00_export_disease_tables.R), the test DOSE::enrichDGN runs on the same data (Yu et al., Bioinformatics 31:608,
  2015), which the published analysis used (scripts/14_disease_enrichment.R) and later DOSE releases dropped: the target
  is the Entrez genes of the methylarginine proteins, the universe the Entrez genes of every reviewed protein (the
  background is those with a DisGeNET term, as in DOSE), terms of 10 to 300 genes. Enriched: pvalue, p.adjust and qvalue
  at most 0.05, as the published run cut them.
- Disease Ontology: clusterProfiler::enricher on the gene sets DOSE's enrichDO tests (each Disease Ontology term with
  the genes of the terms under it; inputs/disease/hdo_<release>_term_genes.tsv, exported by
  scripts/00_export_disease_ontology.R; the same results as enrichDO), the same target and universe (the background is
  the universe genes with a Disease Ontology term), terms of at least 10 genes and no upper limit, so the broad classes
  are tested too.
- Disease families: each DisGeNET term is assigned to at most one family (scripts/disease_families.py); a family's
  genes are the union of its terms' genes. clusterProfiler::enricher tests the families against the genes with a
  DisGeNET term (DOSE's background for the terms), the eight audited families apart from the additional, rule-based
  ones (column audited False).

All three are one-sided over-representation tests (hypergeometric), with Benjamini-Hochberg within each.

    python scripts/14_disease_enrichment.py --labels RESULTS/functional_class_union/protein_functional_labels.tsv \\
        --term-genes inputs/disease/dose_4.4.0_disgenet_term_genes.tsv \\
        --symbol-entrez inputs/disease/orghsegdb_3.22.0_symbol_entrez.tsv \\
        --ontology-genes inputs/disease/hdo_2026-09-30_term_genes.tsv --ontology-terms inputs/disease/hdo_2026-09-30_terms.tsv \\
        --outdir RESULTS/disease

Writes to --outdir:
  disease_terms.tsv             enriched DisGeNET terms (DOSE's columns, gene symbols in geneID)
  disease_terms_tested.tsv      every DisGeNET term DOSE tested, enriched or not, with its family
  disease_ontology.tsv          every Disease Ontology term tested, with its parents and the methylated genes in it
  disease_families.tsv          the family tests: the eight audited families, then the additional ones
  disease_family_proteins.tsv   for each methylarginine protein, the families its gene is in
  disease_background_proteins.tsv  for each methylarginine protein, whether its gene has a DisGeNET term (the background
                                of the family tests)
  disease_term_families.tsv     every DisGeNET term with its family assignment
  disease_input_summary.tsv     proteins, symbols and Entrez genes of the universe and the target
Needs Rscript with clusterProfiler and DOSE on PATH.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd

from common import entrez_by_symbol
from disease_families import ADDITIONAL_FAMILIES, ALL_FAMILIES, AUDITED_FAMILIES, classify
from enrichment import Job, run

MIN_SIZE, MAX_SIZE, CUTOFF = 10, 300, 0.05


def primary_gene(value: str) -> str:
    return str(value or '').split(';')[0].strip()


def is_true(value) -> bool:
    return str(value).strip().lower() in {'1', 'true', 't', 'yes', 'y'}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--labels', required=True)
    ap.add_argument('--term-genes', required=True, help="DOSE's DisGeNET term-gene table (for the family gene sets)")
    ap.add_argument('--symbol-entrez', required=True)
    ap.add_argument('--ontology-genes', required=True, help='Disease Ontology term-gene table (scripts/00_export_disease_ontology.R)')
    ap.add_argument('--ontology-terms', required=True, help='Disease Ontology terms with names and parents')
    ap.add_argument('--uniprot-metadata', default='',
                    help="UniProt entries table: each protein counts under its UniProt primary gene (else the labels' first symbol)")
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    labels = pd.read_csv(args.labels, sep='\t', dtype=str, keep_default_na=False)
    # The protein's UniProt primary gene: the labels' first symbol comes from the source tables first (IGF2 for IGF2-AS).
    uniprot_gene = {}
    if args.uniprot_metadata:
        meta = pd.read_csv(args.uniprot_metadata, sep='\t', dtype=str, keep_default_na=False)
        uniprot_gene = {a: primary_gene(g) for a, g in zip(meta['Entry'], meta['Gene Names (primary)'])}
    labels['gene_symbol'] = [uniprot_gene.get(acc, '') or primary_gene(a) or primary_gene(b) for acc, a, b in
                             zip(labels['canonical_UniProtAC'], labels['substrate_genename'], labels.get('genename', [''] * len(labels)))]
    labels = labels[labels['gene_symbol'] != '']
    target_rows = labels[labels['is_arg_methyl_protein'].map(is_true)]

    entrez_of = entrez_by_symbol(args.symbol_entrez)
    # A gene is named by the symbol of the protein that brought it in, the symbol the release maps back to proteins.
    name_of = {}
    for s in labels['gene_symbol']:
        for e in entrez_of.get(s, []):
            name_of.setdefault(e, s)

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
    target_list = [g for g in universe_list if g in targeted]
    universe, target = set(universe_list), set(target_list)
    unmapped = lambda rows: int(rows.loc[~rows['gene_symbol'].isin(entrez_of), 'gene_symbol'].nunique())  # noqa: E731
    summary = pd.DataFrame([
        {'target': 'Universe', 'protein_count': labels['canonical_UniProtAC'].nunique(),
         'symbol_count': labels['gene_symbol'].nunique(), 'entrez_count': len(universe_list), 'symbols_without_entrez': unmapped(labels)},
        {'target': 'All methyl proteins', 'protein_count': target_rows['canonical_UniProtAC'].nunique(),
         'symbol_count': target_rows['gene_symbol'].nunique(), 'entrez_count': len(target_list), 'symbols_without_entrez': unmapped(target_rows)},
    ])
    summary.to_csv(outdir / 'disease_input_summary.tsv', sep='\t', index=False, lineterminator='\n')
    symbols = lambda ids: '/'.join(name_of.get(g, g) for g in str(ids or '').split('/') if g)  # noqa: E731

    # ---- Family assignment of every DisGeNET term --------------------------------------------------------------------
    tg = pd.read_csv(args.term_genes, sep='\t', dtype=str, keep_default_na=False, quoting=3)
    term_name = dict(zip(tg['term_id'], tg['term_name']))
    term_genes_all = tg.groupby('term_id')['entrez_id'].agg(set).to_dict()
    audit, family_terms = [], defaultdict(list)
    for t in sorted(term_name):
        fam, status, method, reason, secondary = classify(term_name[t])
        audit.append({'term_id': t, 'term_name': term_name[t], 'primary_family': fam, 'audit_status': status,
                      'assignment_method': method, 'secondary_rule_matches': secondary, 'override_reason': reason,
                      'genes_in_reviewed_universe': len(term_genes_all[t] & universe)})
        if status == 'included':
            family_terms[fam].append(t)
    pd.DataFrame(audit).to_csv(outdir / 'disease_term_families.tsv', sep='\t', index=False, lineterminator='\n')
    family_of = {t: f for f, ts in family_terms.items() for t in ts}
    fam_genes = {f: set().union(*(term_genes_all[t] & universe for t in family_terms[f])) if family_terms[f] else set() for f in ALL_FAMILIES}
    annotated = set(tg['entrez_id']) & universe       # DOSE's background for the DisGeNET terms

    # ---- The tests ----------------------------------------------------------------------------------------------------
    do_genes = pd.read_csv(args.ontology_genes, sep='\t', dtype=str, keep_default_na=False)
    do_terms = pd.read_csv(args.ontology_terms, sep='\t', dtype=str, keep_default_na=False)
    fam_sets = lambda names: pd.DataFrame([(f, g) for f in names for g in fam_genes[f]], columns=['term', 'gene'])  # noqa: E731
    jobs = [
        Job('disgenet_terms', target, universe, sets=tg[['term_id', 'entrez_id']], names=tg[['term_id', 'term_name']], min_size=MIN_SIZE, max_size=MAX_SIZE),
        Job('disease_ontology', target, universe, sets=do_genes[['term_id', 'entrez_id']], names=do_terms[['term_id', 'name']], min_size=MIN_SIZE),
        Job('disease_family', target & annotated, annotated, sets=fam_sets(AUDITED_FAMILIES), whole_universe=True),
        Job('disease_family_additional', target & annotated, annotated, sets=fam_sets(ADDITIONAL_FAMILIES), whole_universe=True),
    ]
    res = run(jobs, outdir / 'enricher')

    dgn = res['disgenet_terms'].copy()
    dgn['geneID'] = dgn['geneID'].map(symbols)
    dgn['family'] = dgn['ID'].map(family_of).fillna('')
    dgn['enriched'] = (dgn['pvalue'] <= CUTOFF) & (dgn['p.adjust'] <= CUTOFF) & (dgn['qvalue'] <= CUTOFF)
    dgn.sort_values(['p.adjust', 'Count'], ascending=[True, False], kind='stable').to_csv(outdir / 'disease_terms_tested.tsv', sep='\t', index=False, lineterminator='\n')
    terms = dgn[dgn['enriched']].copy()
    terms['gene_ratio_num'] = terms['k'] / terms['n']
    terms['bg_ratio_num'] = terms['M'] / terms['N']
    terms['fold_enrichment'] = terms['FoldEnrichment']
    terms = terms.sort_values(['p.adjust', 'Count'], ascending=[True, False], kind='stable')
    terms[['ID', 'Description', 'GeneRatio', 'BgRatio', 'RichFactor', 'FoldEnrichment', 'pvalue', 'p.adjust', 'qvalue', 'geneID', 'Count',
           'gene_ratio_num', 'bg_ratio_num', 'fold_enrichment']].to_csv(outdir / 'disease_terms.tsv', sep='\t', index=False, lineterminator='\n')

    do = res['disease_ontology'].copy()
    parents = dict(zip(do_terms['term_id'], do_terms['parents']))
    do['parents'] = do['ID'].map(parents).fillna('')
    do['geneID'] = do['geneID'].map(symbols)
    do.sort_values(['p.adjust', 'Count'], ascending=[True, False], kind='stable')[
        ['ID', 'Description', 'parents', 'GeneRatio', 'BgRatio', 'k', 'n', 'M', 'N', 'FoldEnrichment', 'pvalue', 'p.adjust', 'qvalue', 'geneID']
    ].to_csv(outdir / 'disease_ontology.tsv', sep='\t', index=False, lineterminator='\n')

    fam_rows = []
    for job, audited in (('disease_family', True), ('disease_family_additional', False)):
        for r in res[job].to_dict('records'):
            fam_rows.append({'disease_family': r['ID'], 'contributing_ontology_terms': len(family_terms[r['ID']]),
                             'methylarginine_protein_count': r['k'], 'background_count': r['M'] - r['k'], 'family_gene_count': r['M'],
                             'methylarginine_denominator': r['n'], 'background_denominator': r['N'] - r['n'],
                             'fold_enrichment': r['FoldEnrichment'], 'p_value': r['pvalue'], 'BH_adjusted_p': r['p.adjust'],
                             'audited': audited, 'method': 'enricher'})
    fams = pd.DataFrame(fam_rows)
    fams.sort_values(['audited', 'p_value'], ascending=[False, True], kind='stable').to_csv(outdir / 'disease_families.tsv', sep='\t', index=False, lineterminator='\n')

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
    pd.DataFrame([{'canonical_UniProtAC': acc, 'gene_symbol': s, 'disgenet_background': bool(set(entrez_of.get(s, [])) & annotated)}
                  for acc, s in zip(target_rows['canonical_UniProtAC'], target_rows['gene_symbol'])]).to_csv(
        outdir / 'disease_background_proteins.tsv', sep='\t', index=False, lineterminator='\n')
    print(f'DisGeNET (enricher): {len(dgn)} terms tested, {len(terms)} enriched; Disease Ontology (enricher): {len(do)} terms tested, '
          f'{int((do["p.adjust"] <= CUTOFF).sum())} with p.adjust <= {CUTOFF}; families: {int((fams["BH_adjusted_p"] <= CUTOFF).sum())} of {len(fams)}; '
          f'target {len(target)} of {len(universe)} genes; {len({r["canonical_UniProtAC"] for r in prot})} proteins in a family')


if __name__ == '__main__':
    main()

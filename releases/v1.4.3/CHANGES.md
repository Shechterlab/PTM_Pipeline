# Changes from v1.4.2 to v1.4.3

| | v1.4.2 | v1.4.3 |
|---|---:|---:|
| Sites | 11,632 | 11,632 |
| Proteins | 4,439 | 4,439 |
| Sites kept pending | 0 | 0 |
| Source records | 67,516 | 67,516 |
| Source records placed | 38,325 | 38,325 |
| Studies and databases | 19 | 19 |

Of the 11,632 sites of v1.4.2 (pending ones included), 11,632 are kept and 0 removed; 0 sites were added. The records of v1.4.2 are read in the same order and all are placed as before. 

### Xenograft-only sites

Ma et al. 2023 measured breast cancer patient-derived xenografts grown in mice. The 46 sites only they report carry the flag `xenograft_only`: where the human and mouse sequences agree, the peptide can come from the mouse stroma.

## New annotation layers

v1.4.3 adds a layer of annotation beside the sites. It adds, removes or moves no site; its columns, tables and tests are in `README.md`, `data_dictionary.tsv` and the statistics table, and scripts/release_layers.py builds them.

### Provenance of the interaction records

Where every interaction record at a site comes from, and what it shows there. Each PDB row of `supporting/site_interactions` is split into its entries and each IntAct row into its interactions; the crosslinks and predicted models keep one row each. Every record gets its method, the species of both proteins, its paper, the state of the site in that experiment (the residue a structure models; the methylarginine IntAct curates there; a mutation) and the other marks of the protein, a class (writes, removes, binds methylated, binds unmodified, a mutation changing the binding, a binding region, another reaction, RNA or DNA, not a protein, another copy, predicted) and flags for what to check (`supporting/site_interaction_records`). The PDB entries that model a site as a methylarginine are read in full, including peptides that PDBe-KB does not map to the human protein.

- Every interaction record at a site now says where it comes from (supporting/site_interaction_records, 52,060 rows for 4,701 sites): the PDB entry, IntAct interaction, crosslinking study or predicted model, its method, the species of both proteins, the paper (2,677 papers behind the experimental records), the state of the site in the experiment and the other marks of the protein.
- What the records show: methylates the site 10; demethylates it 4; binds it methylated 230; binds it unmodified 16,797; binds a mutant or variant there 2; at an interface whose residue could not be read 32; mutating it changes the binding 51; mutating it leaves the binding unchanged 111; mutated, effect not stated 3; binds a region holding it 1,114; another reaction on the region 39; RNA or DNA 17,710; not a protein partner 12; another copy of the protein 3,125; predicted only 12,820.
- 18 sites have a partner recorded binding them methylated (230 records from 199 PDB entries and 11 IntAct interactions). 219 PDB entries model a site methylated and are not among those PDBe-KB lists for it: peptides without a human UniProt mapping, orthologues identical around the site, and the other histone variants sharing a tail. Their 234 records are new: binds it methylated 204, RNA or DNA 30.
- Flags for the reader to check (column flags): partner_other_species 908; site_protein_other_species 337; matched_by_sequence 219; state_not_reported 25; paper_names_other_arginine 4; title_names_paralogue 6; no_publication 2,315.
- Site columns interaction_classes, partners_binding_methylated, methylating_enzymes (6 sites), demethylating_enzymes and interaction_papers summarise them.

Source: Burley et al., Nucleic Acids Res 51:D488 (2023); Sehnal et al., Nucleic Acids Res 49:W431 (2021); del Toro et al., Nucleic Acids Res 50:D648 (2022); Montecchi-Palazzi et al., Nat Biotechnol 26:864 (2008); McLaughlin et al., Nucleic Acids Res 53:D1318 (2025); Sayers et al., Nucleic Acids Res 53:D20 (2025), doi:10.1093/nar/gkae979. Licence: PDB data CC0 1.0 (RCSB Data API and ModelServer); IntAct CC BY 4.0; PSI-MI and PSI-MOD (OLS) CC BY 4.0; PubMed citation data under the NLM Terms and Conditions (courtesy of the U.S. National Library of Medicine).

## Corrected: statistics

A consistency audit of every test led to these corrections. The sites are unchanged; only rows of the statistics tables change.

- **Proteins each in-cell experiment found binding RNA or DNA** (`binding_protein`, 7 experiments). Up to v1.4.2 they were Fisher's exact tests; they are now clusterProfiler enricher tests (one-sided hypergeometric), like every other protein-set test, with the same counts and fold enrichment. The p of 4 of the 7 changes, and none crosses q <= 0.05 either way. Each is also repeated within PaxDb abundance deciles (7 new rows, method mantel_haenszel), as the other protein-set tests have been since v1.4.1.
- **95% confidence intervals of the within-protein interaction tests.** 82 Mantel-Haenszel rows of `pdb_interface`, `predicted_interface`, `crosslink`, `sequence_proximity`, `nucleic_acid_distance` and `intact` had no interval; 82 now carry the Robins-Breslow-Greenland interval of the odds ratio (`odds_ratio_ci_low`, `odds_ratio_ci_high`), as every other Mantel-Haenszel row does. Their odds ratios and p are unchanged.
- **`sequence_proximity`**: the pCLAP category is now `pclap_peptide`, as in `crosslink` (10 rows; was `pclap_lysc_region`). The values are unchanged.
- **The pooled residue tests count the sites the within-protein tests count.** The motif and disorder tests (method fisher) counted the 5 sites of MIF, a protein whose every arginine is a site and so has no other arginine to compare them with; the domain tests and the within-protein tests already left it out. The p of 12 rows changes slightly, and none crosses q <= 0.05 either way.
- **The view with lower-confidence sites** (`statistics_with_lower_confidence`). The lower-confidence sites this view does not count (`in_analyses_with_lower_confidence` false) were left out of the methylarginines but still counted among the other arginines of the residue tests, and their proteins counted as methylated in `go_binding` and `binding_protein`. They are now left out of both sides (the lower-confidence view's binding_protein rows also become enricher tests, as above). Of 4,662 rows the two versions share, 248 have a different p (sequence_proximity 60, crosslink 50, domain_class 26, pdb_interface 24, motif 16, structure_context 13, nucleic_acid_distance 12, predicted_interface 12, domain_context 9, go_binding 8, binding_protein 7, disorder 6, intact 5), and 2 became significant at q <= 0.05 (structure_context domain_disordered (mantel_haenszel); binding_protein rna_pclap (enricher)) and none stopped being significant.

The website reads the same tables, so its charts and tables change with them.

## Still open

- **Paralog windows** (11 sites, flag `paralog_window`), accepted in `config/remap_evidence_allowlist.tsv`.
- **Sites without a readable localisation** (Site specificity in the CHANGES.md of v1.2.0): mostly sites whose study gives them below its localisation threshold, or whose search placed the methyl group differently from one confident spectrum to the next.

## Tracing a change

- `supporting/site_changes_v1.4.2_to_v1.4.3.tsv`: one row per site of either version, with its change and the records involved.
- `supporting/source_records_v1.4.3.tsv`: one row per source record, with what the study reported, where each version placed it and the evidence check.

Their columns are listed in `README.md` and `data_dictionary.tsv`.

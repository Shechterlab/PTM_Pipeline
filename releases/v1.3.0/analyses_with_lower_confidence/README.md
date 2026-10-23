# Analyses of v1.3.0 with its lower-confidence sites, rerun with the current pipeline

The analyses of the release, run on its 13,445 sites on 4,749 proteins with the code of commit `4d9dbbb+uncommitted` (2026-10-23), so that every release is analysed the same way (hpc/release_analyses.sh). The sites kept pending (flag one_sided_window) are left out, as in the release build. The lower-confidence sites whose own evidence covers the arginine and leaves the methyl group on an arginine (in_analyses_with_lower_confidence in the release's lower-confidence site table) are counted as sites; the others are left out.

`statistics_v1.3.0.tsv` holds every test in one table (columns as in the release's statistics table: analysis, category, the 2x2 counts a-d, effect, p, BH q, background). The folders hold the result tables:

| File | Holds |
|---|---|
| `functional_classes/arg_methyl_functional_class_broad_enrichment.tsv` | functional classes, methylated proteins vs the reviewed proteome |
| `functional_classes/arg_methyl_functional_class_subclass_enrichment.tsv` | functional subclasses |
| `panther/panther_class_enrichment.tsv` | PANTHER protein classes |
| `motifs/motif_family_enrichment_vs_all_proteome_arginines.tsv` | motif families at the site, vs every arginine of the proteome |
| `motifs/motif_family_enrichment_vs_methylated_protein_arginines.tsv` | motif families, vs the other arginines of the methylated proteins |
| `motifs/positional_amino_acid_enrichment.tsv` | residues at -7..+7 of the site |
| `disorder/disorder_context_summary.tsv` | sites in MobiDB disordered regions |
| `disorder/idr_enrichment_vs_all_arginines.tsv` | disorder, vs every arginine |
| `disorder/idr_enrichment_vs_arginines_in_methylated_proteins.tsv` | disorder, vs the other arginines of the methylated proteins |
| `disorder/idr_edge_bin_enrichment_vs_all_disordered_arginines.tsv` | distance to the edge of the disordered region |
| `domains/domain_context_summary.tsv` | sites inside, near and away from InterPro domains |
| `domains/domain_context_enrichment.tsv` | domain context, vs the other arginines of the methylated proteins |
| `domains/domain_class_enrichment.tsv` | nearest domain class |
| `domains/domain_edge_distance_summary.tsv` | distance to the nearest domain edge |
| `domains/domain_cluster_architecture_summary_ge3.tsv` | domain-cluster architectures (3 or more sites) |
| `other_ptms/ptm_context_comparison.tsv` | disorder and domain context of methylarginines vs other PTM sites |
| `other_ptms/ptm_idr_binary_comparison.tsv` | disordered fraction, methylarginines vs other PTM sites |
| `clustering/nearest_neighbor_summary_5_10_20.tsv` | another methylarginine within 5, 10, 20 residues, vs a same-protein arginine null |
| `clustering/nearest_neighbor_curve_vs_null.tsv` | nearest-neighbour distance curve vs the null |
| `clustering/arg_methyl_multiscale_density_window_summary.tsv` | dense windows by window size |
| `clustering/arg_methyl_multiscale_density_significant_regions.tsv` | significant dense regions |
| `condensates/condensate_enrichment.tsv` | CD-CODE condensates |
| `condensates/condensate_summary.tsv` | condensate universe and methylated fraction |
| `disease/disease_families.tsv` | disease families (DisGeNET terms) |
| `disease/disease_terms.tsv` | DisGeNET terms |
| `disease/disease_ontology.tsv` | Disease Ontology terms |
| `disease/disease_input_summary.tsv` | genes and terms tested |

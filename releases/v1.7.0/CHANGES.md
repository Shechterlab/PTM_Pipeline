# Changes from v1.6.0 to v1.7.0

Everything of v1.6.0 is in v1.7.0 under this release's names, its data unchanged. v1.7.0 adds the mutation propensity of every candidate residue (the decisions of the owner of 20 November 2026: hypermutated samples above 10 coding calls per Mb, with and without them; gnomAD's germline rates for gnomAD and ClinVar; the oncohistone controls over every residue of H3; every known site and every changed residue; ClinVar descriptive; a failed control reported, not blocking).

## New: mutation propensity at every candidate residue

v1.7.0 asks, for each source, whether a residue is changed more (or less) often than its codon's mutability predicts, compared with the other residues of its type in the same protein. It reads the change tables of the release (TCGA, DepMap, ClinVar, gnomAD) and adds:

- `supporting/site_propensity_<source>_v1.7.0.tsv.gz` (tcga_mc3, depmap_24q4, clinvar, gnomad_v4_1): every candidate residue (R, K, S, T, Y) that is a known site of a site set, changed or not, and every other candidate residue with a change: observed changes k, expected changes, O/E, binomial p and BH q, with an unweighted expectation beside it; TCGA and DepMap also without the hypermutated samples; gnomAD also the lower tail (depletion);
- `supporting/site_frequency_by_group_v1.7.0.tsv.gz`: samples changed and profiled per residue and TCGA project or DepMap lineage (for waterfall and oncoprint views; the change tables keep the sample of every change);
- `supporting/sample_burden_v1.7.0.tsv`: every TCGA sample and DepMap line, its coding calls, burden per Mb and whether it is hypermutated;
- `supporting/mutation_rates_v1.7.0.tsv`: the spectra and rates the mutability comes from;
- `supporting/site_set_propensity_v1.7.0.tsv`: per site set, residue, source and cohort, the changes at the set's sites within its proteins against expectation;
- `supporting/propensity_controls_v1.7.0.tsv` and `supporting/propensity_validation_v1.7.0.tsv`: positive controls and a permutation check.

### Method

- **Observed (k).** Per residue: TCGA tumor samples, DepMap cell lines, ClinVar variations or gnomAD variants with a missense or stop change there, from the change tables. n is the sum of k over the protein's residues of the same type.
- **Mutability (w).** For each base of the residue's codon and each substitution that changes the residue (missense or stop), the rate of that substitution in its trinucleotide; w is their sum. The codons are those of the Ensembl 105 transcripts UniProt assigns to the entry's displayed isoform whose coding sequence translates to the 2026_03 canonical sequence (18,717 of 20,416 entries; the others have no weights, and their rows leave the weighted columns empty). A base at an exon edge takes its intron neighbor. An entry encoded by several genes (the histones) takes the mean over its coding sequences.
- **Rates.** TCGA and DepMap: each source's own spectrum, the residue-changing SNVs per pyrimidine-centered trinucleotide class in those coding sequences over the positions where the class would change the residue. gnomAD and ClinVar: gnomAD's germline rates, gnomAD mutation_rate_methylation_bins (gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/model/mutation_rate_methylation_bins.txt.bgz; the gnomAD v2.1.1 flagship paper, Karczewski et al., Nature 581:434, 2020; CC0), with the CpG methylation level of each CpG base from gnomAD's methylation table of the same paper (GRCh37, `methylation.txt.bgz`; positions lifted from GRCh38 with Ensembl's chain) binned as gnomAD's constraint code bins it (level 2 above a mean of 0.6, level 1 above 0.2). gnomAD v4 publishes its methylation sites but no per-context rate table; this table is the one gnomAD publishes.
- **Test.** expected = n × w / (the sum of w over the protein's residues of the type); k ~ Binomial(n, w / sum w). p_upper is the chance of k or more; for gnomAD also p_lower, k or fewer. Benjamini–Hochberg q within each source, cohort and residue type, over every residue tested (every residue of a protein with n ≥ 1 and weights), not only the rows written.
- **Hypermutated samples.** More than 10 coding calls per Mb (Campbell et al., Cell 171:1042 (2017), doi:10.1016/j.cell.2017.09.048), over the 32.3 Mb of the coding sequences above: 1,298 of 10,294 TCGA samples and 449 of 1,929 DepMap lines. Columns ending `_no_hypermutated` repeat the test without them, with spectra from the other samples. DepMap's mutation file keeps synonymous calls only where they touch a splice region, so its burden counts protein-changing calls almost only.
- **ClinVar.** ClinVar's variants reflect clinical testing and submission more than mutation. Its O/E uses the germline rates as gnomAD's does, but is descriptive: `oe_status` says `descriptive` on its rows, and the data dictionary says so.
- **Site sets.** For each set, residue and source, over the proteins with a site of the set: the changes at the set's sites against the sum over the proteins of n times the sites' share of the protein's mutability. p from 10,000 draws, each protein's changes redrawn over its residues by mutability (a fixed seed per row). The same code serves any set of residues.

| Source | Residue | Rows | Residues tested | q < 0.05, more changes than expected | of which known sites | q < 0.05, fewer (gnomAD) |
|---|---|---:|---:|---:|---:|---:|
| TCGA MC3 | R | 189,173 | 575,198 | 222 | 2 |  |
| TCGA MC3 | K | 143,077 | 546,770 | 15 | 3 |  |
| TCGA MC3 | S | 209,844 | 840,495 | 103 | 14 |  |
| TCGA MC3 | T | 81,274 | 497,494 | 14 | 5 |  |
| TCGA MC3 | Y | 33,253 | 213,446 | 5 | 1 |  |
| DepMap 24Q4 | R | 55,440 | 451,714 | 108 | 0 |  |
| DepMap 24Q4 | K | 111,095 | 386,888 | 55 | 11 |  |
| DepMap 24Q4 | S | 136,448 | 655,992 | 181 | 19 |  |
| DepMap 24Q4 | T | 50,571 | 360,913 | 143 | 8 |  |
| DepMap 24Q4 | Y | 21,563 | 150,347 | 30 | 1 |  |
| ClinVar | R | 233,669 | 575,136 | 8 | 0 |  |
| ClinVar | K | 163,913 | 547,931 | 1 | 0 |  |
| ClinVar | S | 234,566 | 831,775 | 6 | 1 |  |
| ClinVar | T | 141,454 | 524,884 | 0 | 0 |  |
| ClinVar | Y | 59,094 | 235,449 | 0 | 0 |  |
| gnomAD v4.1 | R | 501,621 | 601,257 | 76 | 5 | 1 |
| gnomAD v4.1 | K | 397,321 | 614,937 | 0 | 0 | 0 |
| gnomAD v4.1 | S | 636,874 | 888,113 | 1 | 0 | 0 |
| gnomAD v4.1 | T | 425,659 | 567,902 | 0 | 0 | 0 |
| gnomAD v4.1 | Y | 198,507 | 284,658 | 1 | 0 | 0 |

### Site sets

A p of 9.999e-05 is the smallest 10,000 draws give ((1 + 0) / 10,001).

| Site set | Residue | Source | Cohort | Proteins | Sites | Observed | Expected | O/E | O/E (unweighted) | p (more) | p (fewer) |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Methylarginine (R) | R | TCGA MC3 | all | 4,148 | 10,944 | 4,462 | 5091.52 | 0.876359 | 0.888175 | 1 | 9.999e-05 |
| Methylarginine (R) | R | TCGA MC3 | no_hypermutated | 3,615 | 9,923 | 1,192 | 1417.96 | 0.840646 | 0.862487 | 1 | 9.999e-05 |
| Methylarginine (R) | R | DepMap 24Q4 | all | 3,171 | 8,571 | 1,034 | 1165.18 | 0.887417 | 0.886618 | 1 | 9.999e-05 |
| Methylarginine (R) | R | DepMap 24Q4 | no_hypermutated | 2,067 | 5,917 | 317 | 372.739 | 0.850462 | 0.864253 | 0.9994 | 0.0006999 |
| Methylarginine (R) | R | ClinVar | all | 4,102 | 10,770 | 5,979 | 5667.09 | 1.05504 | 1.02285 | 9.999e-05 | 1 |
| Methylarginine (R) | R | gnomAD v4.1 | all | 4,293 | 11,210 | 20,460 | 18953.6 | 1.07948 | 1.04109 | 9.999e-05 | 1 |
| Methyllysine (K) | K | TCGA MC3 | all | 2,426 | 4,475 | 452 | 507.256 | 0.891069 | 0.871399 | 0.9965 | 0.0046 |
| Methyllysine (K) | K | TCGA MC3 | no_hypermutated | 1,716 | 3,406 | 128 | 158.393 | 0.808115 | 0.790068 | 0.9974 | 0.0033 |
| Methyllysine (K) | K | DepMap 24Q4 | all | 1,727 | 3,213 | 171 | 185.103 | 0.92381 | 0.914526 | 0.8773 | 0.1371 |
| Methyllysine (K) | K | DepMap 24Q4 | no_hypermutated | 1,027 | 2,089 | 74 | 70.7614 | 1.04577 | 1.02688 | 0.3529 | 0.6968 |
| Methyllysine (K) | K | ClinVar | all | 2,415 | 4,464 | 808 | 884.604 | 0.913403 | 0.909997 | 0.9964 | 0.004 |
| Methyllysine (K) | K | gnomAD v4.1 | all | 2,716 | 4,956 | 4,179 | 4321.48 | 0.967031 | 0.965907 | 0.9885 | 0.0119 |
| Lysine acylation, ubiquitination and SUMOylation (K) | K | TCGA MC3 | all | 8,333 | 79,669 | 7,494 | 7557.22 | 0.991635 | 0.967144 | 0.8218 | 0.1823 |
| Lysine acylation, ubiquitination and SUMOylation (K) | K | TCGA MC3 | no_hypermutated | 5,071 | 55,196 | 2,117 | 2103.37 | 1.00648 | 0.988779 | 0.3506 | 0.6594 |
| Lysine acylation, ubiquitination and SUMOylation (K) | K | DepMap 24Q4 | all | 5,429 | 56,386 | 2,718 | 2769.09 | 0.981548 | 0.972087 | 0.8931 | 0.1119 |
| Lysine acylation, ubiquitination and SUMOylation (K) | K | DepMap 24Q4 | no_hypermutated | 2,921 | 33,578 | 968 | 990.807 | 0.976982 | 0.968315 | 0.8331 | 0.1771 |
| Lysine acylation, ubiquitination and SUMOylation (K) | K | ClinVar | all | 8,572 | 80,663 | 15,284 | 15530.6 | 0.984121 | 0.981697 | 0.9935 | 0.006999 |
| Lysine acylation, ubiquitination and SUMOylation (K) | K | gnomAD v4.1 | all | 10,170 | 88,582 | 72,657 | 72561.4 | 1.00132 | 0.997978 | 0.3204 | 0.6812 |
| Phosphorylation (S, T, Y) | S | TCGA MC3 | all | 10,256 | 94,548 | 14,033 | 14180.3 | 0.98961 | 0.997537 | 0.926 | 0.07519 |
| Phosphorylation (S, T, Y) | S | TCGA MC3 | no_hypermutated | 7,920 | 82,362 | 3,948 | 3986.99 | 0.99022 | 0.996447 | 0.7718 | 0.234 |
| Phosphorylation (S, T, Y) | S | DepMap 24Q4 | all | 7,601 | 72,277 | 4,180 | 4326.8 | 0.966073 | 0.967767 | 0.9963 | 0.0039 |
| Phosphorylation (S, T, Y) | S | DepMap 24Q4 | no_hypermutated | 4,764 | 52,712 | 1,397 | 1483.44 | 0.941733 | 0.94271 | 0.9971 | 0.0033 |
| Phosphorylation (S, T, Y) | S | ClinVar | all | 10,203 | 94,081 | 22,921 | 22485.3 | 1.01938 | 1.02013 | 0.0003 | 0.9999 |
| Phosphorylation (S, T, Y) | S | gnomAD v4.1 | all | 11,156 | 98,141 | 108,538 | 104819 | 1.03548 | 1.03544 | 9.999e-05 | 1 |
| Phosphorylation (S, T, Y) | T | TCGA MC3 | all | 6,374 | 26,606 | 2,771 | 2804.82 | 0.987941 | 0.985882 | 0.7648 | 0.2397 |
| Phosphorylation (S, T, Y) | T | TCGA MC3 | no_hypermutated | 4,168 | 19,317 | 904 | 888.797 | 1.01711 | 1.01246 | 0.2904 | 0.7237 |
| Phosphorylation (S, T, Y) | T | DepMap 24Q4 | all | 4,461 | 18,878 | 1,257 | 1174.32 | 1.0704 | 1.06918 | 0.0048 | 0.9953 |
| Phosphorylation (S, T, Y) | T | DepMap 24Q4 | no_hypermutated | 2,513 | 11,860 | 484 | 423.358 | 1.14324 | 1.14154 | 0.0006999 | 0.9994 |
| Phosphorylation (S, T, Y) | T | ClinVar | all | 7,060 | 28,461 | 8,972 | 8190.65 | 1.0954 | 1.08353 | 9.999e-05 | 1 |
| Phosphorylation (S, T, Y) | T | gnomAD v4.1 | all | 7,800 | 30,260 | 37,670 | 35232.4 | 1.06919 | 1.06422 | 9.999e-05 | 1 |
| Phosphorylation (S, T, Y) | Y | TCGA MC3 | all | 2,913 | 7,182 | 904 | 890.92 | 1.01468 | 1.01211 | 0.3145 | 0.6993 |
| Phosphorylation (S, T, Y) | Y | TCGA MC3 | no_hypermutated | 1,746 | 4,668 | 363 | 376.408 | 0.96438 | 0.978508 | 0.8021 | 0.2171 |
| Phosphorylation (S, T, Y) | Y | DepMap 24Q4 | all | 2,010 | 5,035 | 490 | 485.727 | 1.0088 | 1.00031 | 0.4141 | 0.6078 |
| Phosphorylation (S, T, Y) | Y | DepMap 24Q4 | no_hypermutated | 881 | 2,301 | 152 | 150.339 | 1.01105 | 1.01441 | 0.454 | 0.5799 |
| Phosphorylation (S, T, Y) | Y | ClinVar | all | 3,356 | 8,045 | 2,628 | 2547.82 | 1.03147 | 1.03679 | 0.0347 | 0.9673 |
| Phosphorylation (S, T, Y) | Y | gnomAD v4.1 | all | 4,326 | 9,773 | 10,070 | 9785.1 | 1.02912 | 1.03642 | 0.0004 | 0.9997 |

### Controls

The oncohistone controls are computed over every residue of the H3 entries (glycine 34 is not a candidate residue), from the TCGA and DepMap changes at every residue of H3 (UniProt counts the initiator methionine: H3K27 is K28, H3K36 is K37, H3G34 is G35). The recurrent hotspots (Chang et al., Nat Biotechnol 34:155 (2016), doi:10.1038/nbt.3391) are read from the per-site tables. A control passes at q < 0.05; a failed control does not hold the release, and its row says why.

TCGA, all samples: H3K27M (UniProt K28) on P68431 fail, H3K36M (UniProt K37) on P68431 pass, H3G34R/V (UniProt G35) on P68431 fail, H3K27M (UniProt K28) on Q71DI3 fail, H3K36M (UniProt K37) on Q71DI3 fail, H3G34R/V (UniProt G35) on Q71DI3 fail, H3K27M (UniProt K28) on P84243 fail, H3K36M (UniProt K37) on P84243 pass, H3G34R/V (UniProt G35) on P84243 fail; hotspots 13 of 13 pass.

| Control | Entry | Source | Cohort | Family | k | n | Expected | q | Result |
|---|---|---|---|---|---:|---:|---:|---:|---|
| H3K27M (UniProt K28) | P68431 (H3C1;H3C10;H3C11;H3C12;H3C2;H3C3;H3C4;H3C6;H3C7;H3C8) | TCGA MC3 | all | every residue of the entry | 2 | 269 | 1.2505 | 1 | fail |
| H3K36M (UniProt K37) | P68431 (H3C1;H3C10;H3C11;H3C12;H3C2;H3C3;H3C4;H3C6;H3C7;H3C8) | TCGA MC3 | all | every residue of the entry | 7 | 269 | 1.3899 | 0.01929 | pass |
| H3G34R/V (UniProt G35) | P68431 (H3C1;H3C10;H3C11;H3C12;H3C2;H3C3;H3C4;H3C6;H3C7;H3C8) | TCGA MC3 | all | every residue of the entry | 4 | 269 | 2.03149 | 1 | fail |
| H3K27M (UniProt K28) | P68431 (H3C1;H3C10;H3C11;H3C12;H3C2;H3C3;H3C4;H3C6;H3C7;H3C8) | TCGA MC3 | no_hypermutated | every residue of the entry | 1 | 110 | 0.460419 | 1 | fail |
| H3K36M (UniProt K37) | P68431 (H3C1;H3C10;H3C11;H3C12;H3C2;H3C3;H3C4;H3C6;H3C7;H3C8) | TCGA MC3 | no_hypermutated | every residue of the entry | 5 | 110 | 0.52026 | 0.02609 | pass |
| H3G34R/V (UniProt G35) | P68431 (H3C1;H3C10;H3C11;H3C12;H3C2;H3C3;H3C4;H3C6;H3C7;H3C8) | TCGA MC3 | no_hypermutated | every residue of the entry | 0 | 110 | 0.921413 | 1 | fail |
| H3K27M (UniProt K28) | P68431 (H3C1;H3C10;H3C11;H3C12;H3C2;H3C3;H3C4;H3C6;H3C7;H3C8) | DepMap 24Q4 | all | every residue of the entry | 4 | 162 | 0.778023 | 0.3659 | fail |
| H3K36M (UniProt K37) | P68431 (H3C1;H3C10;H3C11;H3C12;H3C2;H3C3;H3C4;H3C6;H3C7;H3C8) | DepMap 24Q4 | all | every residue of the entry | 0 | 162 | 0.828625 | 1 | fail |
| H3G34R/V (UniProt G35) | P68431 (H3C1;H3C10;H3C11;H3C12;H3C2;H3C3;H3C4;H3C6;H3C7;H3C8) | DepMap 24Q4 | all | every residue of the entry | 0 | 162 | 1.36638 | 1 | fail |
| H3K27M (UniProt K28) | P68431 (H3C1;H3C10;H3C11;H3C12;H3C2;H3C3;H3C4;H3C6;H3C7;H3C8) | DepMap 24Q4 | no_hypermutated | every residue of the entry | 1 | 73 | 0.377836 | 1 | fail |
| H3K36M (UniProt K37) | P68431 (H3C1;H3C10;H3C11;H3C12;H3C2;H3C3;H3C4;H3C6;H3C7;H3C8) | DepMap 24Q4 | no_hypermutated | every residue of the entry | 0 | 73 | 0.402767 | 1 | fail |
| H3G34R/V (UniProt G35) | P68431 (H3C1;H3C10;H3C11;H3C12;H3C2;H3C3;H3C4;H3C6;H3C7;H3C8) | DepMap 24Q4 | no_hypermutated | every residue of the entry | 0 | 73 | 0.6457 | 1 | fail |
| H3K27M (UniProt K28) | Q71DI3 (H3C13;H3C14;H3C15) | TCGA MC3 | all | every residue of the entry | 0 | 29 | 0.171723 | 1 | fail |
| H3K36M (UniProt K37) | Q71DI3 (H3C13;H3C14;H3C15) | TCGA MC3 | all | every residue of the entry | 1 | 29 | 0.170634 | 1 | fail |
| H3G34R/V (UniProt G35) | Q71DI3 (H3C13;H3C14;H3C15) | TCGA MC3 | all | every residue of the entry | 1 | 29 | 0.28142 | 1 | fail |
| H3K27M (UniProt K28) | Q71DI3 (H3C13;H3C14;H3C15) | TCGA MC3 | no_hypermutated | every residue of the entry | 0 | 16 | 0.0854914 | 1 | fail |
| H3K36M (UniProt K37) | Q71DI3 (H3C13;H3C14;H3C15) | TCGA MC3 | no_hypermutated | every residue of the entry | 1 | 16 | 0.0843074 | 1 | fail |
| H3G34R/V (UniProt G35) | Q71DI3 (H3C13;H3C14;H3C15) | TCGA MC3 | no_hypermutated | every residue of the entry | 0 | 16 | 0.164374 | 1 | fail |
| H3K27M (UniProt K28) | Q71DI3 (H3C13;H3C14;H3C15) | DepMap 24Q4 | all | every residue of the entry | 0 | 0 |  |  | not_tested |
| H3K36M (UniProt K37) | Q71DI3 (H3C13;H3C14;H3C15) | DepMap 24Q4 | all | every residue of the entry | 0 | 0 |  |  | not_tested |
| H3G34R/V (UniProt G35) | Q71DI3 (H3C13;H3C14;H3C15) | DepMap 24Q4 | all | every residue of the entry | 0 | 0 |  |  | not_tested |
| H3K27M (UniProt K28) | Q71DI3 (H3C13;H3C14;H3C15) | DepMap 24Q4 | no_hypermutated | every residue of the entry | 0 | 0 |  |  | not_tested |
| H3K36M (UniProt K37) | Q71DI3 (H3C13;H3C14;H3C15) | DepMap 24Q4 | no_hypermutated | every residue of the entry | 0 | 0 |  |  | not_tested |
| H3G34R/V (UniProt G35) | Q71DI3 (H3C13;H3C14;H3C15) | DepMap 24Q4 | no_hypermutated | every residue of the entry | 0 | 0 |  |  | not_tested |
| H3K27M (UniProt K28) | P84243 (H3-3A;H3-3B) | TCGA MC3 | all | every residue of the entry | 1 | 62 | 0.292345 | 1 | fail |
| H3K36M (UniProt K37) | P84243 (H3-3A;H3-3B) | TCGA MC3 | all | every residue of the entry | 6 | 62 | 0.398159 | 0.0004309 | pass |
| H3G34R/V (UniProt G35) | P84243 (H3-3A;H3-3B) | TCGA MC3 | all | every residue of the entry | 4 | 62 | 0.511944 | 0.1204 | fail |
| H3K27M (UniProt K28) | P84243 (H3-3A;H3-3B) | TCGA MC3 | no_hypermutated | every residue of the entry | 1 | 28 | 0.116933 | 1 | fail |
| H3K36M (UniProt K37) | P84243 (H3-3A;H3-3B) | TCGA MC3 | no_hypermutated | every residue of the entry | 3 | 28 | 0.163893 | 0.04004 | pass |
| H3G34R/V (UniProt G35) | P84243 (H3-3A;H3-3B) | TCGA MC3 | no_hypermutated | every residue of the entry | 4 | 28 | 0.228903 | 0.01063 | pass |
| H3K27M (UniProt K28) | P84243 (H3-3A;H3-3B) | DepMap 24Q4 | all | every residue of the entry | 3 | 20 | 0.0965065 | 0.01638 | pass |
| H3K36M (UniProt K37) | P84243 (H3-3A;H3-3B) | DepMap 24Q4 | all | every residue of the entry | 1 | 20 | 0.120106 | 1 | fail |
| H3G34R/V (UniProt G35) | P84243 (H3-3A;H3-3B) | DepMap 24Q4 | all | every residue of the entry | 1 | 20 | 0.185748 | 1 | fail |
| H3K27M (UniProt K28) | P84243 (H3-3A;H3-3B) | DepMap 24Q4 | no_hypermutated | every residue of the entry | 2 | 11 | 0.0559618 | 0.1878 | fail |
| H3K36M (UniProt K37) | P84243 (H3-3A;H3-3B) | DepMap 24Q4 | no_hypermutated | every residue of the entry | 1 | 11 | 0.0695954 | 1 | fail |
| H3G34R/V (UniProt G35) | P84243 (H3-3A;H3-3B) | DepMap 24Q4 | no_hypermutated | every residue of the entry | 1 | 11 | 0.101963 | 1 | fail |
| IDH1 R132 | O75874 (IDH1) | TCGA MC3 | all | arginines of every protein (R) | 459 | 466 | 53.4221 | 0 | pass |
| IDH2 R172 | P48735 (IDH2) | TCGA MC3 | all | arginines of every protein (R) | 26 | 35 | 0.392092 | 1.411e-38 | pass |
| TP53 R175 | P04637 (TP53) | TCGA MC3 | all | arginines of every protein (R) | 154 | 1140 | 55.668 | 1.47e-24 | pass |
| TP53 R248 | P04637 (TP53) | TCGA MC3 | all | arginines of every protein (R) | 207 | 1140 | 34.6524 | 1.034e-88 | pass |
| TP53 R273 | P04637 (TP53) | TCGA MC3 | all | arginines of every protein (R) | 252 | 1140 | 55.736 | 5.079e-85 | pass |
| PTEN R130 | P60484 (PTEN) | TCGA MC3 | all | arginines of every protein (R) | 98 | 191 | 21.1827 | 8.357e-38 | pass |
| FBXW7 R465 | Q969H0 (FBXW7) | TCGA MC3 | all | arginines of every protein (R) | 59 | 220 | 7.70967 | 8.578e-30 | pass |
| SF3B1 K700 | O75533 (SF3B1) | TCGA MC3 | all | lysines of every protein (K) | 14 | 33 | 0.290185 | 1.835e-15 | pass |
| CTNNB1 S37 | P35222 (CTNNB1) | TCGA MC3 | all | serines of every protein (S) | 42 | 130 | 3.79707 | 9.704e-27 | pass |
| CTNNB1 S45 | P35222 (CTNNB1) | TCGA MC3 | all | serines of every protein (S) | 40 | 130 | 3.79707 | 1.985e-24 | pass |
| U2AF1 S34 | Q01081 (U2AF1) | TCGA MC3 | all | serines of every protein (S) | 30 | 31 | 1.94322 | 5.029e-30 | pass |
| FGFR3 S249 | P22607 (FGFR3) | TCGA MC3 | all | serines of every protein (S) | 41 | 55 | 0.750391 | 5.161e-59 | pass |
| ERBB2 S310 | P04626 (ERBB2) | TCGA MC3 | all | serines of every protein (S) | 44 | 58 | 0.814226 | 2.131e-63 | pass |
| H3K27M (UniProt K28) | P68431 (H3C1;H3C10;H3C11;H3C12;H3C2;H3C3;H3C4;H3C6;H3C7;H3C8) | TCGA MC3 | all | lysines of every protein (K) | 2 | 19 | 2.10362 | 1 | fail |
| H3K36M (UniProt K37) | P68431 (H3C1;H3C10;H3C11;H3C12;H3C2;H3C3;H3C4;H3C6;H3C7;H3C8) | TCGA MC3 | all | lysines of every protein (K) | 7 | 19 | 2.33812 | 1 | fail |
| H3K27M (UniProt K28) | Q71DI3 (H3C13;H3C14;H3C15) | TCGA MC3 | all | lysines of every protein (K) | 0 | 3 | 0.385007 | 1 | fail |
| H3K36M (UniProt K37) | Q71DI3 (H3C13;H3C14;H3C15) | TCGA MC3 | all | lysines of every protein (K) | 1 | 3 | 0.382564 | 1 | fail |
| H3K27M (UniProt K28) | P84243 (H3-3A;H3-3B) | TCGA MC3 | all | lysines of every protein (K) | 1 | 10 | 1.0448 | 1 | fail |
| H3K36M (UniProt K37) | P84243 (H3-3A;H3-3B) | TCGA MC3 | all | lysines of every protein (K) | 6 | 10 | 1.42297 | 1 | fail |
| IDH1 R132 | O75874 (IDH1) | TCGA MC3 | no_hypermutated | arginines of every protein (R) | 442 | 442 | 46.6491 | 0 | pass |
| IDH2 R172 | P48735 (IDH2) | TCGA MC3 | no_hypermutated | arginines of every protein (R) | 26 | 33 | 0.336177 | 8.512e-41 | pass |
| TP53 R175 | P04637 (TP53) | TCGA MC3 | no_hypermutated | arginines of every protein (R) | 133 | 932 | 42.7489 | 1.484e-25 | pass |
| TP53 R248 | P04637 (TP53) | TCGA MC3 | no_hypermutated | arginines of every protein (R) | 185 | 932 | 33.4394 | 3.139e-74 | pass |
| TP53 R273 | P04637 (TP53) | TCGA MC3 | no_hypermutated | arginines of every protein (R) | 218 | 932 | 51.7067 | 3.308e-68 | pass |
| PTEN R130 | P60484 (PTEN) | TCGA MC3 | no_hypermutated | arginines of every protein (R) | 49 | 81 | 8.75569 | 3.016e-22 | pass |
| FBXW7 R465 | Q969H0 (FBXW7) | TCGA MC3 | no_hypermutated | arginines of every protein (R) | 34 | 123 | 5.34761 | 9.509e-14 | pass |
| SF3B1 K700 | O75533 (SF3B1) | TCGA MC3 | no_hypermutated | lysines of every protein (K) | 14 | 27 | 0.264196 | 3.269e-16 | pass |
| CTNNB1 S37 | P35222 (CTNNB1) | TCGA MC3 | no_hypermutated | serines of every protein (S) | 35 | 111 | 2.84429 | 3.227e-23 | pass |
| CTNNB1 S45 | P35222 (CTNNB1) | TCGA MC3 | no_hypermutated | serines of every protein (S) | 33 | 111 | 2.84429 | 7.965e-21 | pass |
| U2AF1 S34 | Q01081 (U2AF1) | TCGA MC3 | no_hypermutated | serines of every protein (S) | 23 | 23 | 1.34317 | 6.138e-24 | pass |
| FGFR3 S249 | P22607 (FGFR3) | TCGA MC3 | no_hypermutated | serines of every protein (S) | 29 | 37 | 0.533898 | 1.038e-40 | pass |
| ERBB2 S310 | P04626 (ERBB2) | TCGA MC3 | no_hypermutated | serines of every protein (S) | 21 | 28 | 0.384052 | 1.985e-28 | pass |
| H3K27M (UniProt K28) | P68431 (H3C1;H3C10;H3C11;H3C12;H3C2;H3C3;H3C4;H3C6;H3C7;H3C8) | TCGA MC3 | no_hypermutated | lysines of every protein (K) | 1 | 11 | 1.00436 | 1 | fail |
| H3K36M (UniProt K37) | P68431 (H3C1;H3C10;H3C11;H3C12;H3C2;H3C3;H3C4;H3C6;H3C7;H3C8) | TCGA MC3 | no_hypermutated | lysines of every protein (K) | 5 | 11 | 1.1349 | 1 | fail |
| H3K27M (UniProt K28) | Q71DI3 (H3C13;H3C14;H3C15) | TCGA MC3 | no_hypermutated | lysines of every protein (K) | 0 | 2 | 0.213665 | 1 | fail |
| H3K36M (UniProt K37) | Q71DI3 (H3C13;H3C14;H3C15) | TCGA MC3 | no_hypermutated | lysines of every protein (K) | 1 | 2 | 0.210706 | 1 | fail |
| H3K27M (UniProt K28) | P84243 (H3-3A;H3-3B) | TCGA MC3 | no_hypermutated | lysines of every protein (K) | 1 | 6 | 0.541454 | 1 | fail |
| H3K36M (UniProt K37) | P84243 (H3-3A;H3-3B) | TCGA MC3 | no_hypermutated | lysines of every protein (K) | 3 | 6 | 0.758905 | 1 | fail |
| IDH1 R132 | O75874 (IDH1) | DepMap 24Q4 | all | arginines of every protein (R) | 11 | 12 | 1.06783 | 5.746e-07 | pass |
| IDH2 R172 | P48735 (IDH2) | DepMap 24Q4 | all | arginines of every protein (R) | 3 | 3 | 0.0661306 | 0.04652 | pass |
| TP53 R175 | P04637 (TP53) | DepMap 24Q4 | all | arginines of every protein (R) | 44 | 391 | 20.9221 | 0.02114 | pass |
| TP53 R248 | P04637 (TP53) | DepMap 24Q4 | all | arginines of every protein (R) | 105 | 391 | 12.7826 | 2.728e-58 | pass |
| TP53 R273 | P04637 (TP53) | DepMap 24Q4 | all | arginines of every protein (R) | 96 | 391 | 17.757 | 1.661e-37 | pass |
| PTEN R130 | P60484 (PTEN) | DepMap 24Q4 | all | arginines of every protein (R) | 16 | 43 | 3.1175 | 0.0002606 | pass |
| FBXW7 R465 | Q969H0 (FBXW7) | DepMap 24Q4 | all | arginines of every protein (R) | 14 | 56 | 1.94069 | 7.017e-05 | pass |
| SF3B1 K700 | O75533 (SF3B1) | DepMap 24Q4 | all | lysines of every protein (K) | 3 | 11 | 0.103541 | 0.3961 | fail |
| CTNNB1 S37 | P35222 (CTNNB1) | DepMap 24Q4 | all | serines of every protein (S) | 9 | 23 | 0.548852 | 2.027e-05 | pass |
| CTNNB1 S45 | P35222 (CTNNB1) | DepMap 24Q4 | all | serines of every protein (S) | 7 | 23 | 0.548852 | 0.004913 | pass |
| U2AF1 S34 | Q01081 (U2AF1) | DepMap 24Q4 | all | serines of every protein (S) |  |  |  |  | fail |
| FGFR3 S249 | P22607 (FGFR3) | DepMap 24Q4 | all | serines of every protein (S) | 4 | 7 | 0.105377 | 0.009475 | pass |
| ERBB2 S310 | P04626 (ERBB2) | DepMap 24Q4 | all | serines of every protein (S) | 5 | 10 | 0.137495 | 0.000913 | pass |
| H3K27M (UniProt K28) | P68431 (H3C1;H3C10;H3C11;H3C12;H3C2;H3C3;H3C4;H3C6;H3C7;H3C8) | DepMap 24Q4 | all | lysines of every protein (K) | 4 | 15 | 1.30569 | 1 | fail |
| H3K36M (UniProt K37) | P68431 (H3C1;H3C10;H3C11;H3C12;H3C2;H3C3;H3C4;H3C6;H3C7;H3C8) | DepMap 24Q4 | all | lysines of every protein (K) | 0 | 15 | 1.39061 | 1 | fail |
| H3K27M (UniProt K28) | Q71DI3 (H3C13;H3C14;H3C15) | DepMap 24Q4 | all | lysines of every protein (K) |  |  |  |  | fail |
| H3K36M (UniProt K37) | Q71DI3 (H3C13;H3C14;H3C15) | DepMap 24Q4 | all | lysines of every protein (K) |  |  |  |  | fail |
| H3K27M (UniProt K28) | P84243 (H3-3A;H3-3B) | DepMap 24Q4 | all | lysines of every protein (K) | 3 | 6 | 0.521576 | 1 | fail |
| H3K36M (UniProt K37) | P84243 (H3-3A;H3-3B) | DepMap 24Q4 | all | lysines of every protein (K) | 1 | 6 | 0.649121 | 1 | fail |
| IDH1 R132 | O75874 (IDH1) | DepMap 24Q4 | no_hypermutated | arginines of every protein (R) | 10 | 10 | 0.898938 | 5.081e-07 | pass |
| IDH2 R172 | P48735 (IDH2) | DepMap 24Q4 | no_hypermutated | arginines of every protein (R) | 3 | 3 | 0.0688509 | 0.05683 | fail |
| TP53 R175 | P04637 (TP53) | DepMap 24Q4 | no_hypermutated | arginines of every protein (R) | 30 | 266 | 12.1427 | 0.0348 | pass |
| TP53 R248 | P04637 (TP53) | DepMap 24Q4 | no_hypermutated | arginines of every protein (R) | 77 | 266 | 10.2434 | 9.124e-40 | pass |
| TP53 R273 | P04637 (TP53) | DepMap 24Q4 | no_hypermutated | arginines of every protein (R) | 60 | 266 | 11.4293 | 1.205e-21 | pass |
| PTEN R130 | P60484 (PTEN) | DepMap 24Q4 | no_hypermutated | arginines of every protein (R) | 4 | 16 | 1.14854 | 1 | fail |
| FBXW7 R465 | Q969H0 (FBXW7) | DepMap 24Q4 | no_hypermutated | arginines of every protein (R) | 5 | 21 | 0.657448 | 0.6958 | fail |
| SF3B1 K700 | O75533 (SF3B1) | DepMap 24Q4 | no_hypermutated | lysines of every protein (K) | 3 | 4 | 0.0398995 | 0.02664 | pass |
| CTNNB1 S37 | P35222 (CTNNB1) | DepMap 24Q4 | no_hypermutated | serines of every protein (S) | 4 | 11 | 0.253227 | 0.2772 | fail |
| CTNNB1 S45 | P35222 (CTNNB1) | DepMap 24Q4 | no_hypermutated | serines of every protein (S) | 4 | 11 | 0.253227 | 0.2772 | fail |
| U2AF1 S34 | Q01081 (U2AF1) | DepMap 24Q4 | no_hypermutated | serines of every protein (S) |  |  |  |  | fail |
| FGFR3 S249 | P22607 (FGFR3) | DepMap 24Q4 | no_hypermutated | serines of every protein (S) | 3 | 5 | 0.0738873 | 0.1375 | fail |
| ERBB2 S310 | P04626 (ERBB2) | DepMap 24Q4 | no_hypermutated | serines of every protein (S) | 5 | 5 | 0.0664205 | 4.811e-06 | pass |
| H3K27M (UniProt K28) | P68431 (H3C1;H3C10;H3C11;H3C12;H3C2;H3C3;H3C4;H3C6;H3C7;H3C8) | DepMap 24Q4 | no_hypermutated | lysines of every protein (K) | 1 | 8 | 0.657523 | 1 | fail |
| H3K36M (UniProt K37) | P68431 (H3C1;H3C10;H3C11;H3C12;H3C2;H3C3;H3C4;H3C6;H3C7;H3C8) | DepMap 24Q4 | no_hypermutated | lysines of every protein (K) | 0 | 8 | 0.700908 | 1 | fail |
| H3K27M (UniProt K28) | Q71DI3 (H3C13;H3C14;H3C15) | DepMap 24Q4 | no_hypermutated | lysines of every protein (K) |  |  |  |  | fail |
| H3K36M (UniProt K37) | Q71DI3 (H3C13;H3C14;H3C15) | DepMap 24Q4 | no_hypermutated | lysines of every protein (K) |  |  |  |  | fail |
| H3K27M (UniProt K28) | P84243 (H3-3A;H3-3B) | DepMap 24Q4 | no_hypermutated | lysines of every protein (K) | 2 | 3 | 0.247673 | 1 | fail |
| H3K36M (UniProt K37) | P84243 (H3-3A;H3-3B) | DepMap 24Q4 | no_hypermutated | lysines of every protein (K) | 1 | 3 | 0.308012 | 1 | fail |

### Permutation check

On up to 1,000 proteins per source (protein and residue-type units with at least two changes, drawn with a fixed seed), each unit's changes were redrawn over its residues by mutability 10,000 times, and the permutation p compared with the binomial p. A residue's count under that redrawing is binomial, so the two should agree up to sampling error; the permutation p cannot go below 1 / 10,001, so it cannot reach q < 0.05 over tens of thousands of residues, and the calls are compared at p < 0.05.

| Source | Units | Residues | Binomial p < 0.05 | Permutation p < 0.05 | Both | Agreement (p < 0.05) | Binomial q < 0.05 | Largest p difference |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| TCGA MC3 | 1,000 | 42,373 | 747 | 744 | 734 | 0.9995 | 5 | 0.01732 |
| DepMap 24Q4 | 1,000 | 52,590 | 1,006 | 1,010 | 986 | 0.9992 | 14 | 0.01347 |
| ClinVar | 1,000 | 37,321 | 557 | 550 | 538 | 0.9992 | 0 | 0.01721 |
| gnomAD v4.1 | 1,000 | 34,760 | 1,401 | 1,399 | 1,378 | 0.9987 | 0 | 0.02253 |

Sizes: `site_propensity_tcga_mc3_v1.7.0.tsv.gz` 18.1 MB, `site_propensity_depmap_24q4_v1.7.0.tsv.gz` 8.9 MB, `site_propensity_clinvar_v1.7.0.tsv.gz` 14.5 MB, `site_propensity_gnomad_v4_1_v1.7.0.tsv.gz` 41.4 MB, `site_frequency_by_group_v1.7.0.tsv.gz` 3.5 MB.

## New files

- `supporting/mutation_rates_v1.7.0.tsv`
- `supporting/propensity_controls_v1.7.0.tsv`
- `supporting/propensity_validation_v1.7.0.tsv`
- `supporting/sample_burden_v1.7.0.tsv`
- `supporting/site_frequency_by_group_v1.7.0.tsv.gz`
- `supporting/site_propensity_clinvar_v1.7.0.tsv.gz`
- `supporting/site_propensity_depmap_24q4_v1.7.0.tsv.gz`
- `supporting/site_propensity_gnomad_v4_1_v1.7.0.tsv.gz`
- `supporting/site_propensity_tcga_mc3_v1.7.0.tsv.gz`
- `supporting/site_set_propensity_v1.7.0.tsv`

## Files renamed

Every carried file is renamed `_v1.6.0` to `_v1.7.0`, its content unchanged except where it names its version: `STUDIES.md`, `analyses/MD5SUMS`, `analyses/README.md`, `analyses_with_lower_confidence/MD5SUMS`, `analyses_with_lower_confidence/README.md`, `human_methylarginine_sites_v1.7.0.xlsx`. `supporting/site_changes_v1.6.0_to_v1.7.0.tsv` replaces v1.6.0's site changes table; it lists every one of the 11,632 methylarginine sites as kept.

Code: `scripts/70_stage_propensity.py` (coding sequences, rates, burden, H3 changes), `scripts/propensity.py` (the method), `scripts/71_release_propensity.py` (this release), `scripts/65_carry_release.py` (carrying the previous release), at commit `a4f91cff557c3f732a11198ee91938db0b5f62c0 (scripts modified)`.

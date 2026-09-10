# Audited disease-ontology consolidation

## Methods and denominators

The analysis used the reviewed-human universe and identifier mapping of the existing pipeline. After SYMBOL-to-Entrez mapping, the methylarginine set contained **4,683 genes** and the non-methylarginine background contained **15,391 genes** (20,074 total). Each ontology term was assigned to at most one primary family using first-match precedence. Exact curated overrides were applied afterward and are visible in the audit table. Family gene sets are unions of unique reviewed-universe genes across assigned terms. Each gene is counted once per family. Two-sided Fisher exact tests were recalculated from each family union, followed by BH correction across the eight final family tests.

## Consolidation decisions

Precedence is ALS/motor-neuron/FTD, Alzheimer/tau/dementia, Parkinsonism, ataxia/cerebellar degeneration, myopathy/muscle/cardiomyopathy, cancer/neoplasia, immune/inflammatory, then developmental/congenital. This prevents FTD labels from entering general dementia. The ALS abbreviation is matched as a standalone token rather than as the substring `als`; this removes false matches inside unrelated words. Ataxia telangiectasia is explicitly excluded from the plotted ataxia family because its defining placement here is genome instability; the audit retains that decision. Four FTD-plus-inclusion-body-myopathy labels are explicitly assigned to ALS/motor-neuron/FTD. Unmatched terms remain `Other / unassigned` and are not plotted. Term assignment never used enrichment magnitude or significance.

## ALS / motor neuron / FTD comparison

The audited family contains **492 methylarginine genes among 1,474 family genes**, OR **1.723**, two-sided p **8.65e-20**, BH-adjusted p **1.15e-19**. The prior regex result was 617/1880 genes, OR 1.698, BH q 9.31e-23. The reduction is driven mainly by correcting the unbounded `als` substring rule, which had admitted unrelated labels, plus mutually exclusive precedence and the documented exact-name overrides. The tested universe is unchanged.

The contributing-term table is ranked by methylarginine-gene count. The contributing-gene table gives every family gene, methylarginine membership, and source terms.

## Family overlap sensitivity

Disease families are not biologically independent. The five largest overlap coefficients are:

- ALS / motor neuron / FTD vs Cancer / neoplasia: 1410 shared genes; Jaccard 0.103; overlap coefficient 0.957 (flagged).
- Alzheimer / tau / dementia vs Cancer / neoplasia: 3245 shared genes; Jaccard 0.236; overlap coefficient 0.945 (flagged).
- Cancer / neoplasia vs Parkinsonism: 1983 shared genes; Jaccard 0.145; overlap coefficient 0.942 (flagged).
- Cancer / neoplasia vs Myopathy / muscle / cardiomyopathy: 2615 shared genes; Jaccard 0.190; overlap coefficient 0.935 (flagged).
- Cancer / neoplasia vs Immune / inflammatory: 6024 shared genes; Jaccard 0.428; overlap coefficient 0.924 (flagged).

## Plot-table agreement and interpretation

Both figures are generated directly from their exact exported TSV tables in the same script. Overview annotations use family-union counts and family-level BH-adjusted p-values. The focused panel retains source term identifiers and names and displays the existing clusterProfiler `p.adjust` values; Storey `qvalue` is retained only for provenance. These results show enrichment of disease-associated gene sets among methylarginine proteins. They do not establish that methylarginine causes disease or that the ontology families are independent mechanisms.

## Recommended legend language

**Disease ontology enrichment among methylarginine proteins.** Disease terms were assigned to mutually exclusive primary families by prespecified label precedence plus documented exact-name overrides. Genes were unioned within each family and counted once. Odds ratios compare mapped methylarginine proteins with the remaining reviewed-human background. P values are from two-sided Fisher exact tests, with Benjamini-Hochberg adjustment across the eight plotted families. The focused panel shows selected direct neurodegenerative terms and uses the original clusterProfiler BH-adjusted `p.adjust` values.

## Top ALS / motor neuron / FTD contributors

### Ontology terms

- Amyotrophic Lateral Sclerosis (C0002736): 345 methylarginine genes among 1051 reviewed-universe genes.
- Frontotemporal dementia (C0338451): 117 methylarginine genes among 300 reviewed-universe genes.
- Motor Neuron Disease (C0085084): 81 methylarginine genes among 183 reviewed-universe genes.
- Frontotemporal Lobar Degeneration (C0751072): 68 methylarginine genes among 187 reviewed-universe genes.
- Motor neuron atrophy (C4024896): 63 methylarginine genes among 131 reviewed-universe genes.
- Amyotrophic Lateral Sclerosis, Sporadic (C1862941): 55 methylarginine genes among 164 reviewed-universe genes.
- Amyotrophic Lateral Sclerosis, Familial (C4551993): 50 methylarginine genes among 127 reviewed-universe genes.
- AMYOTROPHIC LATERAL SCLEROSIS 1 (C1862939): 49 methylarginine genes among 156 reviewed-universe genes.
- GRN-related frontotemporal dementia (C3811918): 42 methylarginine genes among 106 reviewed-universe genes.
- Bulbar palsy (C4082299): 23 methylarginine genes among 48 reviewed-universe genes.

### Methylarginine genes

- SOD1 (6647): present in 27 contributing terms.
- TARDBP (23435): present in 23 contributing terms.
- FUS (2521): present in 22 contributing terms.
- VCP (7415): present in 21 contributing terms.
- SQSTM1 (8878): present in 19 contributing terms.
- CHMP2B (25978): present in 15 contributing terms.
- DCTN1 (1639): present in 14 contributing terms.
- SETX (23064): present in 12 contributing terms.
- ALS2 (57679): present in 11 contributing terms.
- ATXN2 (6311): present in 11 contributing terms.
- MATR3 (9782): present in 11 contributing terms.
- HNRNPA1 (3178): present in 10 contributing terms.
- PFN1 (5216): present in 10 contributing terms.
- UBQLN2 (29978): present in 10 contributing terms.
- GFAP (2670): present in 9 contributing terms.

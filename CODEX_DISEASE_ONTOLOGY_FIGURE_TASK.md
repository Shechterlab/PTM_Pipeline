# Task: regenerate methylarginine disease-ontology figures with audited term consolidation

Work in the `PTM_pipeline` repository. Inspect the existing disease-enrichment analysis and regenerate a concise, manuscript-ready disease-ontology overview, plus a focused neurodegeneration panel emphasizing ALS, FTD, and related biology.

## Existing inputs and code

- Main analysis: `scripts/14_disease_enrichment.R`
- Complete term-level results: `PTM_results/functional_class_union/disease_enrichment/all_methyl_proteins_disease_enrichment.tsv`
- Existing full dot plot: `PTM_results/functional_class_union/disease_enrichment/all_methyl_proteins_disease_dotplot.*`
- Existing manually consolidated results: `PTM_results/functional_class_union/disease_enrichment/all_methyl_proteins_disease_theme_summary.tsv`
- Existing consolidated plot: `PTM_results/functional_class_union/disease_enrichment/all_methyl_proteins_disease_theme_summary.*`
- Existing neurodegeneration terms: `PTM_results/functional_class_union/disease_enrichment/all_methyl_proteins_neurodegenerative_focus.tsv`
- Compact focused panel and plotting code already drafted:
  - `scripts/33_plot_neurodegenerative_disease_focus.py`
  - `PTM_results/functional_class_union/disease_enrichment/all_methyl_proteins_als_ftd_neurodegeneration_compact.*`

Preserve existing files unless deliberately superseding them with clearly named versioned outputs.

## Scientific objective

The full ontology plot is dominated by redundant congenital and dysmorphic terms and does not communicate the neurodegenerative signal. Produce:

1. A consolidated overview containing approximately 6–8 interpretable disease families.
2. A focused panel containing direct neurodegenerative terms, including ALS, FTD, motor-neuron disease, and tauopathy terms when supported by the current results.

The analysis must show enrichment of disease-associated gene sets among methylarginine proteins. Do not imply that methylarginine causes disease.

## Consolidation requirements

Audit the keyword rules currently defined in `scripts/14_disease_enrichment.R`. The current rules can assign the same term to multiple themes. Replace or augment them with an explicit, reproducible term-to-family mapping that has these properties:

- Each ontology term is assigned to at most one primary disease family in the overview.
- Use explicit precedence rules for potentially overlapping matches.
- Keep ALS, motor-neuron disease, and FTD together in one primary family.
- Prevent FTD terms from also entering a general dementia family.
- Prevent terms such as ataxia telangiectasia from entering both ataxia and genome-instability families without an explicit documented decision.
- Separate label matching from the final curated overrides. Curated inclusions and exclusions must be visible in a small configuration table or clearly structured code object.
- Retain an `Other/unassigned` audit status, but do not necessarily plot it.
- Do not select terms merely because they yield a favorable enrichment result.

Suggested overview families:

- ALS / motor neuron / FTD
- Alzheimer / tau / dementia
- Parkinsonism
- Ataxia / cerebellar degeneration
- Myopathy / muscle / cardiomyopathy
- Cancer / neoplasia
- Developmental / congenital syndromes
- Immune / inflammatory

Modify this list if the term audit supports a clearer non-overlapping classification. Document every material change.

## Statistical requirements

- Use the same reviewed-human universe and mapped identifiers as the existing analysis unless inspection reveals an error.
- For each consolidated family, union the unique genes from all assigned ontology terms.
- Count each gene only once within a family.
- Recalculate a two-sided Fisher exact test from the family-level gene union.
- Apply Benjamini–Hochberg correction across the final family-level tests.
- Label the adjusted value unambiguously as `BH-adjusted p` or `q_BH`.
- For the focused term-level panel, use the existing clusterProfiler `p.adjust` column, not the separate Storey `qvalue` column, unless there is a documented reason to change this choice.
- Report at least: family or term, number of contributing ontology terms, methylarginine-protein count, background count, odds ratio or fold enrichment, raw p-value, and BH-adjusted p-value.
- Confirm and report the analyzed methylarginine-set denominator and background denominator after identifier mapping.

## Sensitivity and audit checks

Perform and summarize these checks:

1. Verify that no overview ontology term is assigned to more than one primary family.
2. Export the full term-to-family mapping, including unmatched and excluded terms and the reason for curated overrides.
3. Compare the new ALS/motor-neuron/FTD result with the existing result (currently OR approximately 1.70 and BH q approximately 9.3e-23 from 617 methylarginine proteins among 1,880 theme genes). Explain any difference.
4. List the top contributing terms and genes for the ALS/motor-neuron/FTD family.
5. Quantify gene overlap among the plotted disease families. Flag pairs with substantial overlap rather than treating the families as biologically independent.
6. Confirm that all plotted values exactly match the exported plotting table.

## Figure design

Generate a clean overview and focused panel that can be used separately or combined.

### Consolidated overview

- Horizontal dot, lollipop, or compact forest-style plot.
- Plot odds ratio or fold enrichment on the x-axis with a null reference at 1.
- Show the methylarginine-protein count and BH-adjusted p-value as concise annotations.
- Use a restrained palette compatible with the existing review figures.
- Avoid encoding nearly identical information simultaneously through bar length, dot size, and color.
- Keep labels short and legible at final manuscript size.
- Do not title the figure with a mechanistic or causal claim.

### Focused neurodegeneration panel

- Include direct terms supported by the results, prioritizing motor neuron disease, FTD, familial ALS, and tauopathies.
- Shorten display labels without altering the source ontology names in the exported table.
- Retain source term names and identifiers for provenance.
- Use `p.adjust` for displayed BH-adjusted significance.

## Required outputs

Write outputs under `PTM_results/functional_class_union/disease_enrichment/` with clear names. Produce:

- Consolidated overview in SVG, PDF, and 300-dpi PNG.
- Focused neurodegeneration panel in SVG, PDF, and 300-dpi PNG.
- Optional two-panel composite in SVG, PDF, and 300-dpi PNG if it remains legible.
- Exact plotting tables as TSV files.
- Full term-to-family audit mapping as TSV.
- Family gene-overlap matrix or long-format TSV.
- ALS/motor-neuron/FTD contributing-term and contributing-gene tables.
- A concise Markdown analysis note describing methods, denominators, consolidation decisions, sensitivity results, and recommended figure-legend language.

Update or add a reproducible script under `scripts/`. Integrate it with the existing pipeline only if doing so does not disrupt other outputs.

## Validation before completion

- Render and visually inspect every final figure for clipping, overlap, tiny text, and excessive whitespace.
- Confirm SVG/PDF vector text and shapes render correctly.
- Confirm PNG dimensions and resolution.
- Check every TSV for nonzero dimensions, expected headers, unique keys where appropriate, and agreement with plotted annotations.
- Report all created or modified files.
- Do not overwrite unrelated user changes.

## Final response

Lead with the recommended figure and show its PNG. Briefly state what changed from the previous consolidation, give the key ALS/FTD result, link all final artifacts, and identify any remaining interpretive limitations specific to the data.

# Curated functional experiments

`functional_experiments.tsv` holds one row per experiment and site: what a paper did to a methylarginine (or to a
protein's methylation) and what changed. It is curated by hand from the papers. It is kept apart from the site table
and from every statistic, because a curated list is not a sample of sites. A site can be curated without being a site
of the release: METTL14 R255 is an example, found in mouse cells. `scripts/49_check_functional_curation.py` checks the
table; the roadmap (v1.4.0, item 4) describes the plan for it.

## Columns

| Column | Holds |
|---|---|
| `experiment_id` | `FE` and four digits, never reused |
| `uniprot_accession` | the reviewed human accession, a canonical entry of `../reference/uniprot_human_reviewed_2026_03.fasta` |
| `gene` | its primary gene name |
| `position` | the arginine on that sequence; empty when the experiment is not resolved to one site |
| `region` | for a row without a position, which arginines the experiment concerns |
| `species` | of the cells or protein in the experiment: `human`, `mouse`, `human protein` (recombinant, in vitro) |
| `system` | cells or organism and how they were made, or the in vitro system |
| `writer` | the PRMT the paper names as methylating the site, empty when none |
| `experiment_kind` | a code, below; empty only when the row was read from the abstract alone |
| `perturbation` | what was done, as the paper names it |
| `perturbation_note` | what else the perturbation changes |
| `readout` | what was measured, with its assay where the paper gives it |
| `direction` | of the readout under the perturbation: `increased`, `decreased` or `no_change` |
| `finding` | one sentence, as close to the paper's own words as clarity allows |
| `figure` | where in the paper |
| `reference`, `doi`, `pmid` | the paper |
| `read_from` | `full_text` or `abstract` |
| `curated_on` | date of curation |
| `note` | anything else, such as how a mouse position maps to the human sequence |

## Experiment kinds

| Code | Meaning |
|---|---|
| `methylated_vs_unmethylated` | the same protein or peptide with and without (or with less) methylation, for example from bacteria, from cells treated with a methyltransferase inhibitor, or methylated in vitro |
| `site_mutation` | the arginine mutated; R to K keeps the charge, R to A removes the charge and the side chain, and both remove the methyl group |
| `prmt_inhibitor` | a PRMT inhibitor in cells |
| `prmt_knockdown` | a PRMT knocked down |
| `prmt_knockout` | a PRMT knocked out |
| `correlation` | methylation and a readout measured together, without a perturbation |
| `prediction` | a computational prediction |

## Rules

- **What a row may hold.** Only what the paper shows. Never infer the assay, the direction or the writer.
- **Abstract-only rows.** A row read from the abstract alone says so (`read_from`) and leaves `experiment_kind`
  empty until the full text is read.
- **Mouse positions.** A position from a mouse experiment is given on the human sequence only when the residues
  around it are identical, with the check in `note`.
- **Readouts that are not site-specific.** A PRMT knockdown, knockout or inhibitor changes all of that enzyme's
  substrates, so its readout is not specific to the site. `perturbation_note` says so.

# Human PTM sites by mark, v0.1.0

One table per modification of the human reviewed proteome, from openly licensed sources only, for the sites of the marks family (each with its own questions) to read. The methylarginine sites are not here: they are in [`releases/v1.4.3/`](../../v1.4.3/), whose license is unchanged. This is a first version (v0): UniProt annotations only.

Positions are on UniProt release 2026_03 (reviewed human reference proteome UP000005640, canonical accessions), the reference of every release of this repository. The annotations are from UniProt release 2026_04; an entry whose sequence differs from the reference is left out, and so is any annotation whose residue is not the expected one (dropped in staging: 641 on an isoform, 2 residue mismatch).

## Files

| File | One row per |
|---|---|
| `human_methyllysine_sites_v0.1.0.tsv` | methyllysine site, state and source (533 rows on 397 residues of 176 proteins) |
| `human_phosphorylation_sites_v0.1.0.tsv` | phosphorylation site, state and source (40,928 rows on 40,928 residues of 8,052 proteins) |
| `human_methylarginine_uniprot_crosscheck_v0.1.0.tsv` | UniProt arginine methylation annotation (1,669 rows), with whether the methylarginine release has the arginine. A cross-check, not the methylarginine table |
| `data_dictionary.tsv` | column of each table |
| `release.json` | sources, inputs and their MD5 |
| `MD5SUMS` | checksums of the files |

## Sources and licenses

| Source | License | Used for |
|---|---|---|
| UniProt ([terms](https://www.uniprot.org/help/license)) | CC BY 4.0 | methyllysine, phosphorylation |

**License of this release: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).** Every source above is CC BY 4.0, which allows redistribution and adaptation under CC BY 4.0 with attribution, so the release does not need the share-alike and non-commercial terms of the methylarginine releases. That is a choice for this folder only: the methylarginine release stays CC BY-NC-SA 4.0 (see `LICENSE.md`). A source whose terms are not CC BY, or as open, is not added to this release. Excluded by design: CPLM, PhosphoSitePlus, iPTMnet, dbPTM and qPTM (all rights reserved, non-commercial or unclear terms).

## Evidence classes

- `experimental`: a UniProt record with ECO:0000269 (published experiment) or ECO:0007744 (combinatorial computational and experimental evidence: large-scale proteomics, which is most of UniProt's phosphorylation).
- `by-similarity`: the other UniProt records, mainly ECO:0000250 (by similarity to another species' entry), also ECO:0000255, 0000303 and 0000305, and records with no code. Filter these out to keep measured sites.
- `high-throughput`: a proteomics study that is not a UniProt record. None in this version.

## Counts: methyllysine

| State | by-similarity | experimental | Rows |
|---|---|---|---|
| me | 18 | 10 | 28 |
| me1 | 75 | 181 | 256 |
| me2 | 39 | 88 | 127 |
| me3 | 58 | 64 | 122 |
| all | 190 | 343 | 533 |

## Counts: phosphorylation

| State | by-similarity | experimental | Rows |
|---|---|---|---|
| phospho | 10,371 | 30,557 | 40,928 |
| all | 10,371 | 30,557 | 40,928 |

## Counts: methylarginine cross-check

| State | by-similarity | experimental | Rows |
|---|---|---|---|
| me | 2 | 3 | 5 |
| me1 | 483 | 502 | 985 |
| me2 | 33 | 38 | 71 |
| me2a | 245 | 308 | 553 |
| me2s | 34 | 21 | 55 |
| all | 797 | 872 | 1,669 |

Arginines of the cross-check table in the methylarginine release v1.4.3: lower_confidence 47, no 308, site 1,127.

## How to cite

DeAngelo JD, Silverstein AM, Shechter D. PTM_Pipeline: human PTM sites by mark, v0.1.0. https://github.com/Shechterlab/PTM_Pipeline (Zenodo: [10.5281/zenodo.22696798](https://doi.org/10.5281/zenodo.22696798)), CC BY 4.0, and the sources you use:

- UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. Nucleic Acids Res 53:D609-D617 (2025). doi:10.1093/nar/gkae1010

## Rebuilding

```bash
python scripts/54_stage_uniprot_ptm_sites.py --fasta inputs/reference/uniprot_human_reviewed_2026_03.fasta
python scripts/55_build_marks_release.py --version v0.1.0 --out releases/marks/v0.1.0
```

The UniProt API serves its current release only; `--raw` of script 54 re-stages from a kept download (`inputs/marks/raw_files.tsv` lists its MD5). A published version's folder is never rewritten: a change is a new version.

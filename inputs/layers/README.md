# Layer inputs

These tables hold three annotation layers. The website shows them on protein and site pages; none of them adds
methylarginine sites to the union.

| Layer | Made by | Files |
|---|---|---|
| CPLM 4.0 lysine methylation | `scripts/39_stage_cplm4_lys_methylation.py` | `cplm4_human_lys_methyl_validated.tsv`, `cplm4_human_lys_methyl_rejected.tsv`, `cplm4_summary.json` |
| CARM1 knockout response (Shishkova et al. 2017) | `scripts/41_stage_carm1_ko_sites.py` | `carm1_ko_peptides.tsv`, `carm1_ko_sites.tsv`, `carm1_ko_summary.json` |
| Human Protein Atlas 25.1 subcellular locations | `scripts/42_stage_hpa_subcellular.py` | `hpa_subcellular_locations.tsv`, `hpa_summary.json` |

- `MANIFEST.tsv` describes each file: rows, MD5, what it holds, source, version, citation, licence and how it
  was made.
- `raw_files.tsv` lists the downloads the tables were made from, with URL, size, MD5, SHA-256 and download date.
- The summary JSON files are the scripts' own summaries, with paths shortened to file names. Their input
  checksums equal those stored in the website data (`docs/data/methylarg_sites.json`, `meta.cplm4`, `meta.hpa`,
  `meta.carm1_ko`).

## Checks

On 20 October 2026 the three scripts were run again from the downloads in `raw_files.tsv`. All five tables came
out byte-identical, both times:

- with the FASTA the scripts originally used (the UniProt 2026_03 reviewed human reference proteome,
  20,416 sequences);
- with `../interactions/uniprot_human_reviewed_2026-10-16.fasta`, which is in Git.

## What each table holds

- **`cplm4_human_lys_methyl_validated.tsv`** (5,356 sites on 2,939 proteins): CPLM 4.0 human lysine methylation
  sites whose CPLM sequence matches the canonical UniProt sequence. The match is either an identical sequence, or
  the same ±7-residue window around the lysine. Nothing is remapped. The CPLM id, evidence label (`Exp.`, `Dat.`,
  `Exp./Dat.`) and references are given as CPLM gives them.
- **`cplm4_human_lys_methyl_rejected.tsv`** (559 rows): the rows that failed, with the reason:
  - 494: the accession is not in the reference proteome;
  - 60: no lysine at that position in the canonical sequence;
  - 5: the window differs.
- **`carm1_ko_peptides.tsv`** (515 peptides): ADMA peptides quantified in parental and CARM1 knockout MCF7 and
  MDA-MB-231 cells. Each peptide is placed on the canonical sequence and carries its log2 KO/WT, p value and call:
  - down: log2 ≤ −1 and p < 0.05;
  - up: log2 ≥ 1 and p < 0.05.
- **`carm1_ko_sites.tsv`** (603 rows, 431 sites): one row per methylarginine and cell line, with:
  - median log2 KO/WT;
  - lowest p;
  - peptide counts;
  - the call (down, up, unchanged, mixed).
- **`hpa_subcellular_locations.tsv`** (27,150 rows on 13,301 proteins): one row per protein and HPA location,
  with:
  - whether it is a main (1) or additional (0) location;
  - the HPA reliability (Enhanced, Supported, Approved, Uncertain).

## Held tables

The two CPLM tables are kept out of Git (`.gitignore`) until the PI decides whether they may be redistributed
(see section 9 of `RELEASE_SPEC_v1.1.0.md`). CPLM does not state a licence. The CPLM summary is in Git, and
`MANIFEST.tsv` gives the MD5 of both tables.

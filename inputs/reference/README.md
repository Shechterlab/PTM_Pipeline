# UniProt reference: release 2026_03

The pipeline uses one UniProt reference: the reviewed human entries of the reference proteome UP000005640, from
UniProt release 2026_03 (released 2 September 2026). It covers site positions, sequence checks, layers and
interaction data.

| File | What it is |
|---|---|
| `uniprot_human_reviewed_2026_03.fasta` | The 20,416 canonical sequences. |
| `uniprot_human_reviewed_2026_03_entries.tsv` | Entry, entry name, primary gene, protein names and length of the same entries. |
| `uniprot_accession_history_2026_03.tsv` | Secondary accessions (70,627) and gene-name synonyms (22,537) of these entries. scripts/02 uses them to reach the current entry from an old accession or gene symbol. |
| `uniprot_release.json` | Release, query, download time and MD5 of each file. |
| `v1.1.0_sites_changed_in_2026_03.tsv` | The v1.1.0 sites whose protein changed in this release (see below). |

`MANIFEST.tsv` lists each file with its source URL, MD5, licence (CC BY 4.0) and how it was made.

## How the release is known

The release is read from the `X-UniProt-Release` header of UniProt's response (`2026_03`) on the day of each
download. The same FASTA was downloaded on 18 and 20 October 2026, and the two files are byte-identical (MD5
`c75aad9d197cbba00f51063f65a96ac4`). UniProt has served release 2026_03 since 2 September 2026, so the other
downloads made in between are also 2026_03:

- **The interaction export's UniProt download** (16 October 2026, query `organism_id:9606 AND reviewed:true`):
  it has 15 more entries than the reference proteome. Every entry the two share has the identical sequence, and
  `../interactions/uniprot_human_reviewed_2026-10-16.tsv` marks each entry `in_reference_proteome`.
- **The layers** (`../layers/`): they were made with this FASTA itself.
- **The accession history:** rebuilt on 20 October 2026 against this FASTA. It is identical (MD5) to the copy on
  the release branch.

## What changes for release v1.1.0

Release v1.1.0 (branch `claude/epic-franklin-oni9ft`) placed its 13,368 sites on an earlier FASTA,
`inputs/reference/uniprot_human_reviewed_2026-04-07.fasta`. Its release was not recorded. Compared with 2026_03:

- **13 entries have a different sequence.** v1.1.0 has sites on two of them:
  - POLR2A (P24928): one residue changes (A1835T), outside every site's window, so R426, R1603 and R1810 keep
    their positions;
  - RAB12 (Q6IQ22): 2026_03 adds 96 residues at the N-terminus, so site R10 becomes **R106**.
- **One entry is gone (O43930) and one is new (C0HMG9).** Neither carries a v1.1.0 site.
- **Every other v1.1.0 site** is on a protein whose sequence is identical in both releases.

To move v1.1.0 to 2026_03, the release session should:

1. use these files in place of the 2026-04-07 FASTA and entries;
2. renumber RAB12 R10 to R106 (site ID, windows, the source records that place it, and the distances to its
   domain and disorder boundaries);
3. update `protein_length` and `reference_sequence_md5` for POLR2A and RAB12;
4. re-run the release checks.

`v1.1.0_sites_changed_in_2026_03.tsv` gives the new values.

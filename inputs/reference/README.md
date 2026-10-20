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

`../MANIFEST.tsv` (one manifest for all of `inputs/`) lists each file with its source URL, MD5, licence (CC BY 4.0) and how it was made.

## How the release is known

The release is read from the `X-UniProt-Release` header of UniProt's response (`2026_03`) on the day of each
download. The same FASTA was downloaded on 18 and 20 October 2026, and the two files are byte-identical (MD5
`c75aad9d197cbba00f51063f65a96ac4`). UniProt has served release 2026_03 since 2 September 2026, so the other
downloads made in between are also 2026_03:

- **The interaction export's UniProt download** (16 October 2026, query `organism_id:9606 AND reviewed:true`):
  it has 15 more entries than the reference proteome. Every entry the two share has the identical sequence, and
  `../interactions/uniprot_human_reviewed_2026-10-16.tsv` marks each entry `in_reference_proteome`.
- **The layers** (`../layers/`): they were made with this FASTA itself.
- **The accession history:** rebuilt on 20 October 2026 against this FASTA. It is identical (MD5) to the copy the
  first v1.1.0 build used.

## Earlier reference sequences

Release v1.0.0 (the published union) used UniProt sequences staged on 2026-04-07, whose release was not recorded;
that FASTA is in the Zenodo archive. Compared with 2026_03, 13 entries have a different sequence, one entry is gone
(O43930) and one is new (C0HMG9). For the union this matters in two places:

- POLR2A (P24928): one residue changes (A1835T), outside every site's window, so its sites keep their positions;
- RAB12 (Q6IQ22): 2026_03 adds 96 residues at the N-terminus, so the site published as R10 is R106 here.

The InterPro and MobiDB intervals in `../annotation/` were computed on the earlier sequences and moved to these
(`../annotation/interval_lift_to_uniprot_2026_03.tsv`).

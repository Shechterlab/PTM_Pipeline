# Interaction inputs

These tables hold the data behind the website's interaction views (`website/data/methylarg_interactions.json` and
`website/data/methylarg_pred_iface.json`).

`scripts/37_export_interaction_data.py` downloaded them into a cache between 16 and 19 October 2026
(`data/interactions_cache`: 36,550 files, 658 MB, not in Git). `scripts/37a_stage_interactions.py` then turned
the cache into one plain table per source.

- `../MANIFEST.tsv` (one manifest for all of `inputs/`) describes each table: rows, MD5, what it holds, source URL, version or download date, citation,
  licence, how it was made, and which cached files it came from.
- `raw_cache_files.tsv` lists every cached file with its size, MD5 and the day it was saved.

## Checks

`scripts/37b_check_staged_interactions.py` was run on 20 October 2026 and passed 30 of 30 checks.

- **Sources:** for every source, the data scripts/37 reads from the cache is the same as the data rebuilt from
  these tables alone. Network access is switched off during the check.
- **IntAct and crosslinks:** the IntAct counts and the crosslink sequence check (26,012 matched, 68 dropped)
  equal the values stored in the website data.
- **Measurements:**
  - the predicted-model contacts equal the website's contact gallery (6,267 contacts);
  - every interface-geometry site on the website (1,655) has its row here.

Separately, scripts/37 was re-run from the cache with network access off. It rebuilt both website files
unchanged except for the generation date.

## Conventions

- **Format:** UTF-8, tab-separated, one header row.
- **Values:**
  - lists are joined with `; `;
  - booleans are `true` or `false`;
  - an empty cell means no value.
- **Positions:** canonical UniProt numbering. The exception is `seqres` columns, which hold the residue index in
  the PDB chain (mmCIF `label_seq_id`).
- **File names:** each name carries the source version, or the date of the download.
- **Filtering:** a table is filtered only where its source covers far more than the pipeline can use (other
  species, unreviewed entries, scores below the lowest cut-off the pipeline applies). The analysis choices stay in
  scripts/37, such as the score cut-offs the website shows and the sequence checks.

## Tables

### Proteins

| File | Rows | What it holds |
|---|---:|---|
| `uniprot_human_reviewed_2026-10-16.tsv` | 20,431 | Reviewed human UniProt entries, each with: <ul><li>primary gene name(s);</li><li>STRING identifier;</li><li>sequence length and MD5;</li><li>whether it belongs to the InterPro Tudor (IPR002999), WD40 repeat (IPR001680) or PHD finger (IPR001965, IPR019787) entries.</li></ul> |
| `uniprot_human_reviewed_2026-10-16.fasta` | 20,431 | The canonical sequences scripts/37 checks crosslinks, IntAct features and RBDmap fragments against. |
| `uniprot_human_regions_2026-10-16.tsv` | 72,666 | UniProt disordered regions (Region features noted "Disordered"), domains, repeats and zinc fingers. |
| `uniprot_human_dna_binding_regions_2026-10-16.tsv` | 725 | UniProt DNA-binding regions. |
| `quickgo_human_rna_dna_binding_2026-10-16.tsv` | 27,656 | GO annotations to RNA binding (GO:0003723) or DNA binding (GO:0003677), or to their descendants. One row per protein, term and evidence code. |

**UniProt release.** These UniProt tables are release 2026_03, the release the whole pipeline uses (see
[`../reference/`](../reference/README.md)).

- **This download:** made on 16 October 2026 with `organism_id:9606 AND reviewed:true`. It has 15 more entries
  than the reference proteome.
- **Shared entries:** each of the 20,416 entries in both sets has the identical sequence, and
  `in_reference_proteome` marks them.
- **Rebuilding on the reference:** `scripts/37_export_interaction_data.py --reference-fasta
  inputs/reference/uniprot_human_reviewed_2026_03.fasta` restricts the export to the reference entries. It stops
  if any downloaded sequence differs from the reference.

### Protein-protein networks and complexes

| File | Rows | What it holds |
|---|---:|---|
| `string_v12.5_human_physical_links_score400.tsv` | 414,284 | STRING v12.5 physical-subnetwork links between two reviewed proteins with combined score ≥ 400, each pair once. The website shows the links ≥ 700 that touch a methylated protein. |
| `complex_portal_human_complexes_2026-10-16.tsv` | 88,953 | Complex Portal human complexes: 2,498 curated and 15,284 predicted (hu.MAP 3.0). One row per member slot and accession. A slot that can be filled by any of several interchangeable proteins has one row per protein. |

STRING identifiers are matched to UniProt through UniProt's own STRING cross-reference. When more than one entry
lists the same identifier, the first entry in UniProt's order is used. Links to identifiers that match no
reviewed entry are left out.

### Experimental structures (PDB)

| File | Rows | What it holds |
|---|---:|---|
| `pdbe_kb_interface_residues_2026-10-16.tsv` | 214,124 | PDBe-KB interface residues: stretches of a protein that contact a partner (a protein, RNA or DNA) in at least one PDB entry, with the entries. Consecutive residues with the same partner and the same PDB entries are joined into one start-end row. |
| `pdbe_kb_interface_queries_2026-10-16.tsv` | 5,090 | Which proteins were looked up in PDBe-KB, and whether the answer listed interfaces. A protein that is not listed was not looked up. |
| `sifts_observed_ranges_human_2026-10-16.tsv` | 17,609 | Residues of reviewed human proteins that are modelled in at least one PDB chain, as merged ranges. |
| `sifts_pdb_chain_uniprot_human_2026-10-16.tsv` | 225,706 | SIFTS mapping of PDB chains to reviewed human proteins. Each segment has its first SEQRES index and its UniProt range. |
| `pdb_methylarginine_residues_2026-10-16.tsv` | 302 | Methylated arginines modelled in PDB entries: 2MR (symmetric dimethyl), DA2 (asymmetric dimethyl) and NMM (monomethyl). |
| `pdb_arginine_contacts_2026-10-16.tsv` | 14,989 | Residues of other chains within 4 Å of the side chain of an interface arginine (first model). A row with empty contact columns is an arginine that was examined and touches nothing. |
| `pdb_arginine_interface_environments_2026-10-16.tsv` | 8,176 | Geometry of each interface arginine from the atoms within 10 Å:<ul><li>relative solvent accessibility alone and in the complex;</li><li>Levy core, rim or support class;</li><li>surface buried by nucleic acid and by protein;</li><li>base or backbone hydrogen bonds;</li><li>protein hydrogen bonds.</li></ul> |

The contacts and environments come from the RCSB ModelServer. They were measured by scripts/37 (`site_contacts`,
`analyze_environment`) on the neighbourhoods it saved, and only the measurements are kept here.

### Predicted complexes

| File | Rows | What it holds |
|---|---:|---|
| `protvar_predicted_interfaces_2024.05.28_pdockq0.23.tsv` | 108,851 | ProtVar interface residues (within 5 Å) of AlphaFold-Multimer models of human protein pairs with pDockQ ≥ 0.23, the lowest cut-off the website uses. Residues are listed as ProtVar gives them (`R169,A203,…`). |
| `protvar_served_models_af2complexes_v1.tsv` | 68,756 | Pairs whose model the ProtVar API serves, with chains, pDockQ and interface residues. The residue numbers are written as runs (`3-37,39,43`). |
| `predicted_model_arginine_contacts_2026-10-16.tsv` | 6,267 | Methylarginines within 5 Å of the partner in ProtVar models (pDockQ ≥ 0.5) and AlphaFold DB human heterodimers. Each row gives:<ul><li>closest distance;</li><li>salt bridge, hydrogen bond, cation-π, clash;</li><li>pLDDT;</li><li>the partner residues.</li></ul> |
| `alphafold_db_human_heterodimers_2026-10-16.tsv` | 3,095 | AlphaFold DB human heterodimer models that include a protein of the published union, with chains, ipTM and pDockQ. |
| `alphafold_db_complex_queries_2026-10-16.tsv` | 5,088 | Which proteins were looked up in the AlphaFold DB complex API, and how many human heterodimers each answer listed. |

### Curated interaction features (IntAct)

| File | Rows | What it holds |
|---|---:|---|
| `intact_mutations_human_2026-10-19.tsv` | 61,726 | Mutations with their effect on an interaction (no effect, disrupting, decreasing, increasing, causing). |
| `intact_binding_regions_human_2026-10-19.tsv` | 192,541 | Binding regions: sufficient, necessary, binding-associated, direct. |
| `intact_methylarginines_human_2026-10-19.tsv` | 30 | Methylated arginines curated on interaction participants. |

These tables come from the IntAct release finished on 12 January 2026. They keep reviewed canonical human
accessions only.

- **`sequence_check`** compares IntAct's curated sequence with the canonical UniProt sequence. scripts/37 uses only
  the rows marked `matches` or `no curated sequence`.
- **`curated_sequence`** is filled in only where the check fails. Where it matches, the curated sequence is the
  canonical sequence at `ranges`.

### RNA and DNA crosslinks

| File | Rows | What it holds |
|---|---:|---|
| `pepr_ms_v2_table_s2_rna_crosslinked_peptides.tsv` | 16,799 | pepR-MS compendium, Table S2: RNA-crosslinked peptide-spectrum matches from five datasets, re-analysed with one pipeline. Each row includes the residue carrying the nucleotide mass. |
| `pepr_ms_v2_table_s3_focused_sites.tsv` | 2,211 | pepR-MS compendium, Table S3: the focused crosslinking site(s) of each crosslink. |
| `rbdmap_castello2016_rna_bound_fragments.tsv` | 1,610 | RBDmap (HeLa) high-confidence RNA-bound fragments, as Table S1 gives them. scripts/37 keeps 1,376: 83 of these it moves to the place the fragment matches the sequence, and it drops 234. |
| `trendel2025_dna_crosslink_localisation.tsv` | 1,688 | DNA-crosslinked peptides (MCF7). Lower-case letters mark the residues the crosslink may sit on. |
| `trendel2025_dna_crosslinked_peptides.tsv` | 688 | DNA-crosslinked peptide ranges. |
| `trendel2025_mcf7_dna_interactome.tsv` | 1,191 | Proteins in the MCF7 DNA interactome. |

### Condensates

| File | Rows | What it holds |
|---|---:|---|
| `cdcode_2.3_condensate_members_2026-10-16.tsv` | 6,992 | CD-CODE 2.3 members of the 147 condensates the pipeline uses. Each row gives:<ul><li>confidence score;</li><li>roles;</li><li>kinds of experimental evidence (`unspecified` where CD-CODE gives an empty value);</li><li>PubMed IDs.</li></ul> |

## Licences

Each table carries its source's licence in `../MANIFEST.tsv`. Licences verified on 20 October 2026:

- **CC BY 4.0, CC0 or EMBL-EBI terms:** UniProt, GO, STRING, Complex Portal, PDBe-KB and SIFTS, the PDB, ProtVar,
  AlphaFold DB and IntAct.
- **RBDmap** (Castello et al. 2016): CC BY 4.0, from the open access article in
  [PMC5003815](https://pmc.ncbi.nlm.nih.gov/articles/PMC5003815/).
- **Trendel et al. 2025:** CC BY-NC 4.0, from the open access article.
- **CD-CODE:** CC BY-NC 4.0, as stated on [cd-code.org](https://cd-code.org/).
- **pepR-MS:** CC BY-NC.

The non-commercial terms are compatible with the repository's CC BY-NC-SA 4.0 for data.

## What is not here

- **The raw downloads.** These are:
  - bulk files;
  - one API answer per protein or structure;
  - the atoms around each examined residue;
  - the journal spreadsheets.

  `raw_cache_files.tsv` records each one's MD5.
- **Data for the sites and proteins that release v1.1.0 adds.** This cache was built for the published union of
  13,719 sites. Of the 13,368 v1.1.0 sites, 873 are not in that union, and 94 v1.1.0 proteins were never looked up
  in PDBe-KB or AlphaFold DB.
  - The proteome-wide tables already cover them: UniProt, GO, STRING, Complex Portal, SIFTS, the ProtVar interface
    table, IntAct and the crosslink tables.
  - The per-site structure measurements and those per-protein lookups need a top-up download, dated separately.

## Rebuilding

```
python scripts/37a_stage_interactions.py --cache data/interactions_cache --out inputs/interactions
python scripts/37b_check_staged_interactions.py --cache data/interactions_cache --tables inputs/interactions
```

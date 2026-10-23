# Changes from v1.2.0 to v1.3.0

| | v1.2.0 | v1.3.0 |
|---|---:|---:|
| Sites | 11,601 | 11,601 |
| Proteins | 4,432 | 4,432 |
| Sites kept pending | 0 | 0 |
| Source records | 66,660 | 66,660 |
| Source records placed | 37,518 | 37,518 |
| Studies and databases | 18 | 18 |

Of the 11,601 sites of v1.2.0 (pending ones included), 11,601 are kept and 0 removed; 0 sites were added. The records of v1.2.0 are read in the same order and all are placed as before. 

### Xenograft-only sites

Ma et al. 2023 measured breast cancer patient-derived xenografts grown in mice. The 46 sites only they report carry the flag `xenograft_only`: where the human and mouse sequences agree, the peptide can come from the mouse stroma.

## RNA and DNA at methylarginines

The sites are those of v1.2.0. This version adds the RNA and DNA in contact with them in PDB structures, read in every entry PDBe-KB lists for a site with RNA or DNA, and what their proteins bind in vitro and in cells.

- `supporting/nucleic_acid_contacts_v1.3.0.tsv`: 256 methylarginines in direct contact with RNA and 90 with DNA (a nucleotide within 4 Å of the side chain), 14,212 site–chain pairs in 1,272 entries: the chain (molecule name, source organism and preparation method: Synthetically obtained for a designed oligonucleotide), the nucleotides contacted and the sequence window of 8 nucleotides on each side of the nearest one, so the windows of different sites align. The windows are the sequences deposited in each structure, often chosen by its authors, so they record contacts, not sequence preference; no consensus is derived from them.
- `supporting/nucleic_acid_proximal_sites_v1.3.0.tsv`: 157 sites with RNA and 63 with DNA within 8 Å of the arginine side chain in a complex where it has no direct contact (every entry of the complex analysed): the minimum heavy-atom distance and the sequence window around the nearest nucleotide.
- `supporting/arginine_nucleic_acid_distances_v1.3.0.tsv`: the minimum heavy-atom distance from the side chain to RNA and to DNA of every arginine of the methylated proteins modelled in an entry with nucleic acid (931 methylarginines, 12,416 other arginines), and in `interaction_statistics` (analysis nucleic_acid_distance) whether methylarginines lie within 4 or 8 Å more often than the other arginines of the same proteins.
- `supporting/nucleic_acid_near_sites_v1.3.0.tsv`: 135 sites with RNA and 68 with DNA in contact via a flanking residue (one at i ± 1–5 within 4 Å, in a complex where the arginine side chain is not within 8 Å): the bound sequence there, centred on the nucleotide the closest such residue contacts.
- `supporting/nucleic_acid_regions_v1.3.0.tsv`: 4,217 runs of residues within 4 Å of RNA or DNA on 606 proteins (289 with methylarginines inside), read once per distinct complex (the same bound sequences), with the segments each contacts.
- `supporting/rna_bind_n_seq_v1.3.0.tsv`: ENCODE RNA Bind-n-Seq (in vitro): the most enriched 5-mers and logos, for 76 methylated proteins.
- `supporting/htr_selex_v1.3.0.tsv`: HTR-SELEX (in vitro, Jolma et al., Genome Res 30:962, 2020, CC BY-NC 4.0): linear and stem-loop RNA motifs, for 53 methylated proteins.
- `supporting/dna_motifs_v1.3.0.tsv`: JASPAR 2026 CORE vertebrates DNA-binding profiles (CC BY 4.0), for 183 methylated proteins.
- `supporting/eclip_v1.3.0.tsv`: eCLIP in K562 and HepG2 cells: peaks by transcript region (GENCODE v50 basic) and the genes with the most peaks, for 186 methylated proteins (137 from ENCODE, 86 from ENCORE (Her et al., Cell Systems 2026, GEO GSE315407)).

The website shows them on protein pages (RNA and DNA interactions: where the protein binds, the bound sequences in contact with each methylarginine with a logo and a type filter, its preferences and eCLIP targets), on a site's card, in the site tables (Bound RNA/DNA, with the search term bound:), and for a selection in the Interactions view. Each chain is typed (rRNA, tRNA, snRNA, mRNA, other ncRNA, viral, designed oligo; nucleosomal, telomeric or other DNA) and the site table gains bound_* columns. Inputs: `inputs/interactions/pdb_nucleic_acid_*` (PDB data, CC0) and `inputs/nucleic_acid/` (scripts/43_stage_nucleic_acid_binding.py). RNAcompete (CISBP-RNA) states no licence and is not included.

## RNA crosslinking in living cells

Besides the crosslinked residues (pepR-MS compendium) and RBDmap fragments of v1.2.0, three experiments in living cells are added to `supporting/residue_evidence`, each with its own evidence type:

- CAPRI (Panhale et al., Nat Commun 2019, doi:10.1038/s41467-019-10585-3, HEK293): peptides crosslinked to RNA and the peptides released next to them: 1,678 stretches on 381 proteins (`evidence_type` capri_xl_peptide, capri_adj_peptide).
- pCLAP (Mullari et al., J Proteome Res 16:2762, 2017, doi:10.1021/acs.jproteome.7b00042, HEK293): RNA-binding peptides and the Lys-C peptides around them: 7,070 stretches on 1,172 proteins (`evidence_type` pclap_peptide, pclap_lysc_region).
- RBR-ID (Oksuz et al., Mol Cell 83:2449, 2023, doi:10.1016/j.molcel.2023.06.012, K562 nuclei): residues under a peptide depleted when RNA is crosslinked, and regions of residue score >= 3: 6,700 stretches on 1,062 proteins (`evidence_type` rbrid_residue, rbrid_region).

A peptide experiment can only report arginines in peptides it measured, so `interaction_statistics` (analysis crosslink) also tests each one among those arginines alone (categories ending in _among_measured, from the measured peptides staged in `inputs/interactions/`): methylarginines against the other arginines of the same proteins, each protein a stratum: pCLAP OR 1.26 (p 0.032, 307 proteins); CAPRI OR 1.37 (p 0.011, 165 proteins); RBDmap OR 1.67 (p 0.00077, 150 proteins); RBR-ID OR 1.43 (p 1.5e-08, 661 proteins). `sequence_proximity` tests whether methylarginines lie within a given distance along the sequence of each kind of RNA-binding residue more often than the other arginines of the same proteins.

## Structural context

Each site gains `structure_context`: disorder (MobiDB) and domains (InterPro) in one class, from inside a folded domain to a disordered region more than 30 residues from any domain, with disordered regions binned by the distance to the nearest domain edge (within 10, 11-30, more than 30). `statistics_v1.3.0.tsv` (analysis structure_context) compares methylarginines with the other arginines of the same proteins in each class. The disorder and domain columns and their tests are unchanged.

## Condensates against the proteome

Each CD-CODE condensate is also tested against the reviewed human proteome (analysis condensate_proteome): 30 of 40 condensates are enriched at q <= 0.05. The test against the other condensate proteins (analysis condensate) is unchanged; as about half of all condensate proteins are methylated, it tests which condensates are enriched relative to other condensates; few are.

## Still open

- **Paralog windows** (11 sites, flag `paralog_window`), accepted in `config/remap_evidence_allowlist.tsv`.
- **Sites without a readable localisation** (Site specificity in the CHANGES.md of v1.2.0): mostly sites whose study gives them below its localisation threshold, or whose search placed the methyl group differently from one confident spectrum to the next.

## Tracing a change

- `supporting/site_changes_v1.2.0_to_v1.3.0.tsv`: one row per site of either version, with its change and the records involved.
- `supporting/source_records_v1.3.0.tsv`: one row per source record, with what the study reported, where each version placed it and the evidence check.

Their columns are listed in `README.md` and `data_dictionary.tsv`.

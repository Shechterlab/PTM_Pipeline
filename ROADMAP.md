# Roadmap

What is planned after v1.3.1, in order of priority. The aim of the next release is to say, for a site, which enzyme
methylates it, under what conditions, and with what consequence, with each kind of evidence kept apart.

## v1.4.0

### Work packages

1. **PRMT7 and monomethylarginine.** Szewczyk et al. 2020, *Pharmacological inhibition of PRMT7 links arginine
   monomethylation to the cellular stress response*, Nat Commun 11:2396,
   [doi:10.1038/s41467-020-16271-z](https://doi.org/10.1038/s41467-020-16271-z): PRMT7 inhibition (SGC3027) and
   knockout, with Rme1 immunoprecipitation; raw data MassIVE
   [MSV000084773](https://www.omicsdi.org/dataset/massive/MSV000084773). A new writer-perturbation layer beside the
   Maron et al. 2021 inhibitors and the CARM1 knockout, under the same peptide rule (`scripts/peptide_calls.py`).
2. **PRMT5 substrates from ETD.** Lu et al. 2024, *ETD-Based Proteomic Profiling Improves Arginine Methylation
   Identification and Reveals Novel PRMT5 Substrates*, J Proteome Res 23:1014–1027,
   [doi:10.1021/acs.jproteome.3c00724](https://doi.org/10.1021/acs.jproteome.3c00724): sDMA sites localised with
   electron-transfer dissociation in HeLa cells, 138 of them new and proposed as PRMT5 substrates. Staged as a study
   (its sites, under the five site tests) and, where its tables give them, its PRMT5-dependent changes as a
   perturbation layer.
3. **Perturbations already measured by included studies.** Hartel et al. 2019, Mol Cell Proteomics 18:2149–2164,
   [doi:10.1074/mcp.RA119.001625](https://doi.org/10.1074/mcp.RA119.001625), PRIDE
   [PXD012357](https://www.ebi.ac.uk/pride/archive/projects/PXD012357), is staged for its sites only; its PRMT1
   knockdown in HEK293T (127 sites changed at q < 0.05, 119 of them Rme1) is in its supplementary tables, which
   `inputs/` already holds as raw downloads. Read the other included studies the same way for knockdowns, knockouts and
   inhibitor arms.
4. **A curated layer of functional experiments.** One row per experiment and site: the protein and site, the
   experiment (methylated against unmethylated protein or peptide; a site mutation, with R→K and R→A noted as removing
   the charge or the residue as well as the methyl group; a PRMT inhibitor, knockdown or knockout; a correlation or a
   prediction), the readout and its direction, the cell or system, and the reference. Kept apart from the site table
   and from every statistic, since curation is not a sample of sites. First entries, as examples of the kinds:
   - METTL14 R255, methylated by PRMT1; the R255K mutant lowers global mRNA m6A in mouse ES cells (Liu et al. 2021,
     Nat Commun 12:3780, [doi:10.1038/s41467-021-24035-6](https://doi.org/10.1038/s41467-021-24035-6));
   - FUS, whose arginine methylation reduces its phase separation, with transportin as a chaperone (Hofweber et al.
     2018, Cell 173:706–719, [doi:10.1016/j.cell.2018.03.004](https://doi.org/10.1016/j.cell.2018.03.004); Qamar et
     al. 2018, Cell 173:720–734, [doi:10.1016/j.cell.2018.03.056](https://doi.org/10.1016/j.cell.2018.03.056)).

### Annotation sources

5. **PhosphoSitePlus, directly.** The sites now come from the 2021 download compiled in Maron et al. 2021 Table S5,
   without PhosphoSitePlus's evidence categories. A direct, dated download
   ([phosphosite.org](https://www.phosphosite.org)) would give them, and the curated-only sites (site field
   `curated_only`) could then be split by evidence.
6. **Protein abundance.** Matched total proteomes come first, where a study has them: they correct a peptide's change
   for its protein's change in the same samples. Otherwise [PaxDb](https://pax-db.org) integrated abundances, as a
   detection covariate for the enrichment tests (abundant proteins are more often reported as methylated).
7. **Cell-line context.** [DepMap](https://depmap.org/portal/) dependency and the CCLE proteomes, for the PRMTs and
   for methylated proteins, as protein-page context.
8. **Pathways.** [Reactome](https://reactome.org) as a further protein-set annotation, under the same
   over-representation test and with the same reported-site caveat.

Each source is staged into `inputs/` with its terms in `inputs/MANIFEST.tsv` before it is used; a source whose terms do
not allow redistribution is kept `checksum_only`, as now.

### Analysis

- Leave-one-study-out versions of the main enrichment and within-protein tests, so a result that rests on one study
  shows it.
- Measured backgrounds wherever a study gives the peptides it measured (as RBR-ID, pCLAP and RBDmap already do for the
  RNA-binding tests).
- A detection control for the CAPRI and RBDmap tests, and SDMA against ADMA where the neutral-loss calls allow it.

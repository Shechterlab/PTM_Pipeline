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
   `curated_only`) could then be split by evidence. Its current
   [terms](https://www.phosphosite.org/psp/legal/terms-and-conditions) give no right to download or redistribute data
   without a written agreement with Cell Signaling Technology, so this waits for one; the same terms bear on the
   PhosphoSitePlus-derived sites the releases carry now, which is for the lab to settle with them.
6. **Protein abundance.** Matched total proteomes come first, where a study has them: they correct a peptide's change
   for its protein's change in the same samples. Otherwise [PaxDb](https://pax-db.org) integrated abundances, as a
   detection covariate for the enrichment tests (abundant proteins are more often reported as methylated).
7. **Cell-line dependency ([DepMap](https://depmap.org/portal/)).** DepMap Public 24Q4 from Figshare (CC BY 4.0,
   [doi:10.25452/figshare.plus.27993248](https://doi.org/10.25452/figshare.plus.27993248)), the latest release
   published there; later releases are downloaded behind a browser check from the portal, under its own terms. PRISM
   Repurposing Public 24Q2 (CC BY 4.0, [doi:10.6084/m9.figshare.25917643](https://doi.org/10.6084/m9.figshare.25917643))
   gives single-dose (2.5 µM, 5-day) viability across its cell lines for five PRMT inhibitors: three of PRMT5
   (EPZ015666, GSK3326595, JNJ-64619178), EPZ020411 (PRMT6) and SGC-707 (PRMT3); it has no GSK591, MS023 or
   MRTX1719. DepMap's screens are of whole genes, so the residue-level use is the somatic mutations it calls
   (`OmicsSomaticMutations`), with TCGA from [cBioPortal](https://www.cbioportal.org) beside them: a first, pooled count on
   v1.3.1 finds 771 of 11,601 methylated arginines changed in at least one cell line (6.6%) against 6.7% of the other
   arginines of the same proteins (odds ratio 1.00), so the card reports mutations at a site as observations, not as a
   test result. Uses, each labelled as cell-line genetics, not methylation:
   - a **protein-page card**, on light pages too: the CRISPR gene effect (Chronos) across cell lines, whether the gene
     is commonly essential, and the CCLE protein level where measured;
   - a **writer layer**: the dependencies on each PRMT, including the known PRMT5 dependency of MTAP-deleted lines
     (Kryukov et al. 2016, Science 351:1214, [doi:10.1126/science.aad5214](https://doi.org/10.1126/science.aad5214);
     Mavrakis et al. 2016, Science 351:1208, [doi:10.1126/science.aad5944](https://doi.org/10.1126/science.aad5944)),
     and which methylated proteins' dependencies correlate with a PRMT's across lines, as a suggestion of a pathway
     beside the perturbation layers, never as evidence that the PRMT methylates the site; and, on each PRMT's page,
     the PRISM viability of cell lines treated with its inhibitors;
   - a **site-level mutation row**: the cell lines (DepMap) and tumours (TCGA) whose protein carries a change at the
     arginine, with AlphaMissense where DepMap gives it;
   - a **protein-set test**: whether methylated proteins are more often commonly essential than other proteins, with
     abundance as a covariate, since essential and abundant proteins are also the most often detected.
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

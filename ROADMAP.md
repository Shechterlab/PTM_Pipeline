# Roadmap

What is planned after v1.3.1, in order of priority. The aim of the next release is to say, for a site, which enzyme
methylates it, under what conditions, and with what consequence, with each kind of evidence kept apart.

## v1.4.0

### Work packages

Status on 25 October 2026 (branch `feature/depmap`): each source below is downloaded, its licence read and its staging
script written; the staged tables wait in a scratch folder for the decisions noted against them, and none is in
`inputs/` yet except the curated layer and DepMap.


1. **PRMT7 and monomethylarginine.** Szewczyk et al. 2020, *Pharmacological inhibition of PRMT7 links arginine
   monomethylation to the cellular stress response*, Nat Commun 11:2396,
   [doi:10.1038/s41467-020-16271-z](https://doi.org/10.1038/s41467-020-16271-z): PRMT7 inhibition (SGC3027) and
   knockout, with Rme1 immunoprecipitation; raw data MassIVE
   [MSV000084773](https://www.omicsdi.org/dataset/massive/MSV000084773). A new writer-perturbation layer beside the
   Maron et al. 2021 inhibitors and the CARM1 knockout, under the same peptide rule (`scripts/peptide_calls.py`).
   *Staged* (`scripts/45`): the paper's proteomics is one comparison, PRMT7 knockout against wild-type HCT116 (Rme1
   immunoprecipitation, SILAC); SGC3027 was tested by western blot only. 343 of 346 peptides placed (99 by residue
   alone, to be checked against the sequence versions searched); the paper's rule (log2 <= -1, BH-adjusted p < 0.01)
   gives its 29 depleted peptides. Open: the licence of the supplementary files (none stated beside the article's CC BY).
   Li et al. 2021, already included, has the SGC3027 and PRMT7 knockdown arms with statistics: stage them next.
2. **PRMT5 substrates from ETD.** Lu et al. 2024, *ETD-Based Proteomic Profiling Improves Arginine Methylation
   Identification and Reveals Novel PRMT5 Substrates*, J Proteome Res 23:1014–1027,
   [doi:10.1021/acs.jproteome.3c00724](https://doi.org/10.1021/acs.jproteome.3c00724): sDMA sites localised with
   electron-transfer dissociation in HeLa cells, 138 of them new and proposed as PRMT5 substrates. Staged as a study
   (its sites, under the five site tests) and, where its tables give them, its PRMT5-dependent changes as a
   perturbation layer.
   *Staged* (reader for `scripts/00`, study `lu2024`): 856 records, 401 localised sites (ptmRS >= 98 %), 20 new; the
   same lab as `lu2025` but a different enrichment and deposit. The PRMT5 inhibitor arm (EPZ015666, SILAC) is staged
   with an inferred rule (twofold, p < 0.05; heavy read as treated), which reproduces 254 of the 255 sites of its Table
   1.4. Open: the supporting information came from ACS's figshare deposit (CC BY-NC 4.0) after pubs.acs.org refused.
3. **Perturbations already measured by included studies.** Hartel et al. 2019, Mol Cell Proteomics 18:2149–2164,
   [doi:10.1074/mcp.RA119.001625](https://doi.org/10.1074/mcp.RA119.001625), PRIDE
   [PXD012357](https://www.ebi.ac.uk/pride/archive/projects/PXD012357), is staged for its sites only; its PRMT1
   knockdown in HEK293T (127 sites changed at q < 0.05, 119 of them Rme1) is in its supplementary tables, which
   `inputs/` already holds as raw downloads. Read the other included studies the same way for knockdowns, knockouts and
   inhibitor arms.
   *Staged* (`scripts/46`): Hartel's 119 are changed peptides, not Rme1 sites (the tables' flag gives them exactly);
   751 sites, 129 changed. Open: the paper's Perseus permutation FDR on 2 against 2 samples keeps 132 peptides, the
   release's BH rule 2. A survey of the other included studies found 40 perturbation arms
   (Li 2021, Fong 2019, Musiani 2019, Larsen 2016 and Fedoriw 2019 the most useful).
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
   *Started*: `inputs/curation/functional_experiments.tsv` with its README and checker (`scripts/49`), six rows from
   these three papers.

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
   *Staged* (`scripts/47`): PaxDb 6.1 (CC BY 4.0), 18,805 entries; methylated proteins have a median of 4.1 ppm
   against 0.16 ppm for the other reviewed proteins. The whole-organism set leans on a 2011 PeptideAtlas build: a
   cell-line set may suit the covariate better.
7. **Cell-line dependency ([DepMap](https://depmap.org/portal/)).** DepMap Public 24Q4 from Figshare (CC BY 4.0,
   [doi:10.25452/figshare.plus.27993248](https://doi.org/10.25452/figshare.plus.27993248)), the latest release
   published there; later releases are downloaded behind a browser check from the portal, under its own terms. PRISM
   Repurposing Public 24Q2 (CC BY 4.0, [doi:10.6084/m9.figshare.25917643](https://doi.org/10.6084/m9.figshare.25917643))
   gives single-dose (2.5 µM, 5-day) viability across its cell lines for five PRMT inhibitors: three of PRMT5
   (EPZ015666, GSK3326595, JNJ-64619178), EPZ020411 (PRMT6) and SGC-707 (PRMT3); it has no GSK591, MS023 or
   MRTX1719. DepMap's screens are of whole genes, so the residue-level use is the somatic mutations it calls
   (`OmicsSomaticMutations`), with TCGA PanCancer Atlas from [cBioPortal](https://www.cbioportal.org) beside them. A
   preview is on the `feature/depmap` branch (`scripts/44*`, `inputs/layers/`); the release build now carries it
   (`scripts/release_layers.py`: four site columns, four protein columns, analyses `depmap_dependency` and
   `cancer_change`). Open: the TCGA columns derive from cBioPortal data under the ODbL, whose share-alike clause
   reaches a derived database; keep them on the website only, publish them under the ODbL, or take MC3 from the GDC. Counting missense and stop changes on the
   canonical sequence, 7.3% of the 11,601 methylated arginines of v1.3.1 are changed in at least one DepMap cell line and
   28.3% in at least one TCGA tumour, against 6.9% and 30.6% of the other arginines of the same proteins: within protein
   (Mantel-Haenszel), odds ratio 0.95 (0.88-1.03) in cell lines and 0.89 (0.85-0.93) in tumours, not adjusted for codon
   context (CpG-containing CGN codons mutate most). So the page reports mutations at a site as observations. Uses, each
   labelled as cancer genetics, not methylation:
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
   *Staged* (`scripts/48`): Reactome 97 (CC0), 11,433 entries in 2,883 pathways; release 98 appeared on 22 October.

Each source is staged into `inputs/` with its terms in `inputs/MANIFEST.tsv` before it is used; a source whose terms do
not allow redistribution is kept `checksum_only`, as now.

### Analysis

- Leave-one-study-out versions of the main enrichment and within-protein tests, so a result that rests on one study
  shows it.
- Measured backgrounds wherever a study gives the peptides it measured (as RBR-ID, pCLAP and RBDmap already do for the
  RNA-binding tests).
- A detection control for the CAPRI and RBDmap tests, and SDMA against ADMA where the neutral-loss calls allow it.

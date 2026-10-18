# Methylarginine Interaction-Regulation Pipeline

This repository contains the analysis code associated with DeAngelo JD, Silverstein AM, Shechter D. "Protein arginine methylation rewires molecular interactions." *Cell Chemical Biology* 33(10):1404–1421 (2026). [doi:10.1016/j.chembiol.2026.09.011](https://doi.org/10.1016/j.chembiol.2026.09.011)

An interactive, searchable version of the methylarginine site union lives in [`docs/`](docs/) (see [Site Explorer](#site-explorer)).

The pipeline assembles, remaps, annotates, and analyzes a human methylarginine site union with an emphasis on questions that should be broadly useful to investigators studying methylarginine biology: source provenance, canonical residue mapping, sequence motifs, disorder/domain context, local site density, condensate/RNP-associated classes, and representative substrate architectures.

The repository is intentionally code-only. Large source tables, downloaded annotations, and generated outputs are staged outside Git.

## Repository Layout

| Path | Purpose |
| --- | --- |
| `scripts/00_*` | Static source download/staging helpers for UniProt, InterPro, and MobiDB context. |
| `scripts/01_*` | Source parsers for methylarginine datasets. |
| `scripts/02_*` | Canonical remapping and source integration. |
| `scripts/03_*` to `scripts/16_*` | Main analysis modules: source provenance, functional classes, motifs, disorder/domain context, PTM comparisons, clustering, condensates, GO/disease supplements, and Maron state/treatment checks. |
| `scripts/17_*` to `scripts/32_*` | Extended analyses and summary-figure utilities, including local-density nulls, domain-edge length nulls, architecture schematics, and 2D nearest-neighbor density panels. |
| `scripts/common.py`, `scripts/clustering_utils.py` | Shared utilities. |
| `scripts/36_export_site_explorer_data.py` | Exports the site union and its annotations as JSON for the static site explorer. |
| `docs/` | Static Methylarginine Site Explorer (HTML/CSS/JS plus exported JSON data), served by GitHub Pages or Netlify. |
| `config/` | JSON ontologies used for functional and domain-class annotation. |
| `envs/` | Conda environment definition for the full Python/R workflow. |
| `hpc/` | Environment, hygiene, and Slurm runner scripts. |

## Staged Inputs

By default, [`hpc/ptm_env.sh`](hpc/ptm_env.sh) resolves local or HPC inputs from `PTM_data/` and outputs from `PTM_results/`. You can override any path with the exported `PTM_*` variables defined there.

The expected local staging layout is:

```text
PTM_data/
  base/
    human_ptm_master.tsv
  external/
    Table_S5_Compiled_Methylarginine_Data.xlsx
    mmc2.xlsx
    mobidb_human_2026-03-30.json
  context/
    uniprot_human_reviewed_canonical.fasta
    uniprot_human_reviewed_canonical_metadata.tsv
    protein2ipr.dat.gz
    interpro_human_reviewed_all_entry_intervals.tsv
    interpro_human_reviewed_domain_like_intervals.tsv
    mobidb_human_reviewed_disorder_intervals.tsv
  condensates/
    cdcode_staged_membership.tsv
    proteins_202603181653.csv
    condensates_202603181647.csv
```

Only the first five files listed below are required for the core methylarginine union and sequence-confirmed motif/clustering analyses. InterPro, MobiDB, CD-CODE, GO, disease, and Maron state/treatment inputs enable additional annotation or supplement modules.

Required inputs for the main methylarginine union:

| Variable | Default path |
| --- | --- |
| `PTM_MASTER_TSV` | `PTM_data/base/human_ptm_master.tsv` |
| `PTM_MARON_XLSX` | `PTM_data/external/Table_S5_Compiled_Methylarginine_Data.xlsx` |
| `PTM_PROMETHEUS_MMC2` | `PTM_data/external/mmc2.xlsx` |
| `PTM_CANONICAL_FASTA` | `PTM_data/context/uniprot_human_reviewed_canonical.fasta` |
| `PTM_UNIPROT_METADATA` | `PTM_data/context/uniprot_human_reviewed_canonical_metadata.tsv` |

Recommended context inputs:

| Variable | Default path |
| --- | --- |
| `PTM_INTERPRO_BULK` | `PTM_data/context/protein2ipr.dat.gz` |
| `PTM_INTERPRO_DOMAIN_INTERVALS` | `PTM_data/context/interpro_human_reviewed_domain_like_intervals.tsv` |
| `PTM_MOBIDB_JSON` | `PTM_data/external/mobidb_human_2026-03-30.json` |
| `PTM_MOBIDB_INTERVALS` | `PTM_data/context/mobidb_human_reviewed_disorder_intervals.tsv` |

Optional modules use staged CD-CODE condensate files and Maron Table S2/S3 workbooks if available. These files are not redistributed in this repository.

### Data Staging

Create the staging directories and source the path configuration:

```bash
mkdir -p PTM_data/base PTM_data/external PTM_data/context PTM_data/condensates PTM_results
. hpc/ptm_env.sh
```

Place user-supplied source tables in the staged paths expected by [`hpc/ptm_env.sh`](hpc/ptm_env.sh):

```text
$PTM_MASTER_TSV          Cross-PTM human PTM master table used for background comparisons.
$PTM_MARON_XLSX          Maron Table S5 compiled methylarginine workbook.
$PTM_PROMETHEUS_MMC2     ProMetheusDB mmc2 methylarginine workbook.
$PTM_MOBIDB_JSON         Optional MobiDB human JSON export for disorder intervals.
```

Then stage public sequence and metadata context:

```bash
python scripts/00_stage_context_sources.py \
  --outdir "$PTM_CONTEXT_ROOT" \
  --skip-interpro
```

If a MobiDB JSON export is available, convert it to intervals:

```bash
python scripts/00_stage_mobidb_disorder.py \
  --input "$PTM_MOBIDB_JSON" \
  --accessions "$PTM_MASTER_TSV" \
  --outdir "$PTM_CONTEXT_ROOT"
```

For InterPro domain intervals, the preferred public workflow is to download `protein2ipr.dat.gz` into `$PTM_INTERPRO_BULK` and parse only proteins observed in the remapped methylarginine union. If you have not remapped yet, run the main pipeline once with `PTM_INTERPRO_MODE=api` or omit domain modules until `PTM_results/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv` exists.

```bash
curl -L -o "$PTM_INTERPRO_BULK" \
  https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/protein2ipr.dat.gz

python scripts/00_stage_context_sources.py \
  --outdir "$PTM_CONTEXT_ROOT" \
  --accessions "$PTM_RESULTS_ROOT/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv" \
  --interpro-source bulk \
  --interpro-protein2ipr "$PTM_INTERPRO_BULK" \
  --bulk-parser pandas
```

Validate the staged paths before a full run:

```bash
. hpc/ptm_env.sh
for f in \
  "$PTM_MASTER_TSV" \
  "$PTM_MARON_XLSX" \
  "$PTM_PROMETHEUS_MMC2" \
  "$PTM_CANONICAL_FASTA" \
  "$PTM_UNIPROT_METADATA"
do
  test -s "$f" || { echo "missing or empty: $f"; exit 1; }
done
```

To stage data outside the repository, set `PTM_DATA_ROOT` and `PTM_RESULTS_ROOT` before sourcing the environment script:

```bash
export PTM_DATA_ROOT=/path/to/PTM_data
export PTM_RESULTS_ROOT=/path/to/PTM_results
. hpc/ptm_env.sh
```

## Environment

Create the Conda environment:

```bash
conda env create -n ptm_pipeline -f envs/ptm_pipeline_hpc.yml
conda activate ptm_pipeline
```

For local plotting on shared or restricted filesystems, set writable cache directories:

```bash
export HOME=/tmp/ptm_home
export XDG_CACHE_HOME=/tmp/ptm_cache
export MPLCONFIGDIR=/tmp/ptm_mpl
mkdir -p "$HOME" "$XDG_CACHE_HOME" "$MPLCONFIGDIR"
```

If Arial is not installed but you want Matplotlib outputs to use it, point `PTM_FONT_DIRS` to one or more directories containing `.ttf` files, separated by `:` on Linux/macOS.

## Main Pipeline

Source the path configuration:

```bash
. hpc/ptm_env.sh
```

Run the full Slurm workflow on HPC:

```bash
bash hpc/create_ptm_env.sh envs/ptm_pipeline_hpc.yml ptm_pipeline
sbatch hpc/slurm_run_integration.sh "$PWD"
```

The Slurm runner stages context, parses source tables, remaps methylarginine sites, integrates the union, annotates motif/disorder/domain context, runs clustering and cross-PTM comparisons, and runs optional condensate/GO/disease modules when their inputs are available.

Useful run-time overrides:

```bash
PTM_INTERPRO_MODE=bulk sbatch hpc/slurm_run_integration.sh "$PWD"
PTM_INTERPRO_MODE=api sbatch hpc/slurm_run_integration.sh "$PWD"
PTM_CLUSTER_RESTRICT_CONTEXT=IDR PTM_CROSS_PTM_CLUSTER_RESTRICT_CONTEXT=IDR sbatch hpc/slurm_run_integration.sh "$PWD"
PTM_MARON_EXCLUDE_PATTERN=Larsen sbatch hpc/slurm_run_integration.sh "$PWD"
```

For a local, stepwise rerun, use the commands in [`hpc/slurm_run_integration.sh`](hpc/slurm_run_integration.sh) as the authoritative ordered recipe.

### Rerun after a pipeline fix and compare with the published results

```bash
cd /path/to/PTM_Pipeline && git pull
sbatch hpc/slurm_rerun_remap_fix.sh "$PWD"                       # new results in ../PTM_results_remapfix
# or choose the location:  sbatch hpc/slurm_rerun_remap_fix.sh "$PWD" /path/to/PTM_results_remapfix
```

The job reuses the staged inputs under `PTM_DATA_ROOT` (so no newer UniProt or InterPro release enters the comparison), runs `hpc/slurm_run_integration.sh` and `scripts/run_extended_analyses.sh` into the new root, rebuilds the explorer's site data there (`site_explorer_data/`), and finishes with:

```bash
python scripts/38_compare_results_roots.py --old PTM_results --new ../PTM_results_remapfix \
  --outdir ../PTM_results_remapfix/rerun_comparison
```

`rerun_comparison/report.md` lists the site-union difference (sites and proteins added or removed, per-source counts), every changed headline `summary.tsv`, p/q values that cross 0.05 and effect sizes that change sign; `file_summary.tsv` covers every table and figure, `changes_detail.tsv` the cell-level changes and `site_union_diff.tsv` each added or removed site. The published root is never written to. Options: `PTM_RERUN_INTERACTIONS=1` also rebuilds the interaction data (needs internet), `PTM_RERUN_EXTENDED=0` skips the extended analyses. `PTM_SKIP_HYGIENE=1` (the default in this job) skips the repo hygiene guard, which fails on a checkout that tracks `PTM_results/`. To publish the new results in the explorer, copy `site_explorer_data/methylarg_sites.json` and `methylarg_evidence.json` to `docs/data/` and rerun `scripts/37_export_interaction_data.py`. The scripts themselves expose explicit `--help` arguments and write tables/figures under `PTM_RESULTS_ROOT`.

## Extended Analyses and Summary Figures

After the main outputs exist, run the extended analysis and summary-figure workflow:

```bash
bash scripts/run_extended_analyses.sh "$PWD"
```

This runner produces:

| Output folder | Main contents |
| --- | --- |
| `PTM_results/summary_figures/` | Representative substrate architecture SVGs, interaction-redistribution schematic, and density/domain-adjacency schematic. |
| `PTM_results/domain_context/` | IDR-length-conditioned domain-edge null, signed domain-edge cumulative density, clustered architecture classifications. |
| `PTM_results/methyl_arg_clustering/` | Nearest-neighbor plots, multiscale rolling-density nulls, pairwise occupancy heatmaps, and 2D k-nearest methylarginine density panels. |

Runtime controls:

```bash
PTM_EXTENDED_PERMUTATIONS=1000 PTM_EXTENDED_2D_PERMUTATIONS=250 bash scripts/run_extended_analyses.sh "$PWD"
```

## Key Output Families

| Folder | Interpretation |
| --- | --- |
| `integrated/` | Canonically remapped methylarginine union, provenance, support summaries, accession audit. |
| `qc/source_provenance/` | Source-support and accession-level provenance summaries. These are useful audit outputs, not primary biological figures. |
| `functional_class_union/` | Multi-label RNA/RNP-oriented functional class enrichment and membership tables. |
| `motif_arg_odds/` | Positional amino-acid enrichment, RG/RGG/GAR motif summaries, acidic asymmetry, and per-protein arginine odds. |
| `disorder_context/` | Disorder/IDR occupancy and IDR-edge localization summaries. |
| `domain_context/` | Domain placement, nearest-domain class enrichment, domain-edge proximity summaries, and length-null controls. |
| `methyl_arg_clustering/` | Local methylarginine spacing, rolling-density, and same-protein arginine-null clustering analyses. |
| `ptm_compare/` | Cross-PTM disorder/domain context comparisons. |
| `condensates/`, `assembly_clustering/` | Optional CD-CODE condensate enrichment and shared-assembly clustering. |
| `maron_state_treatment/` | Optional Maron-only state/treatment supplement analyses. |
| `summary_figures/` | Lightweight SVG schematics for communicating representative architectures and interaction-redistribution concepts. |

## Site Explorer

`docs/` is a static, client-side web app for querying the methylarginine union by gene symbol or synonym, UniProt accession, protein name, site (`FUS R216`) or pasted gene list. It filters by study-level provenance, methylation state, Maron methyl-IP (Rme1, Rme2a, Rme2s) and PRMT inhibitor response, motif, disorder and domain context, clustering, nearby PTMs, CD-CODE condensates, functional class/subclass and disease family; a bar under the search box always shows the active constraints with one-click removal. Each protein page has a protein map (methylarginine lollipops coloured by state, methyl-IP, PRMT inhibitor response or number of reports, over UniProt domains, disorder, composition, RG/RGG-rich segments, UniProt variants and other PTMs) that exports as SVG or PNG; a 3D view of the AlphaFold DB model or any mapped PDB entry with the methylarginines highlighted (entries that model a methylated arginine are flagged), where any residue can be clicked to zoom to it and show its 5 Å contacts, with options for style, colouring, surface, ambient-occlusion lighting, other chains, ligands and site labels; a Save PDF button that prints the whole protein page as a report (light palette, vector map, 3D snapshot, all table columns); per-site pLDDT, AlphaMissense and UniProt variant/modification annotations; an Interactions card (STRING physical partner network coloured by each partner's own methylation, Complex Portal complexes, Gene Ontology RNA/DNA binding with QuickGO evidence, ENCODE eCLIP/ChIP-seq and ChIP-Atlas, PDB entries with RNA or DNA, and methylarginines at PDB interfaces); domain, condensate, ontology and disease panels; and per-site source evidence. The Interactions tab ranks complexes by methylated members, counts methylated proteins bound by each PRMT, Tudor-domain, WD40-repeat and PHD-finger protein, and tests whether methylated arginines sit at protein, RNA or DNA interfaces more often than other arginines of the same protein (Cochran–Mantel–Haenszel, also within ordered regions and within domains). Protein pages keep protein-specific evidence first: CD-CODE confidence, role and evidence for each condensate membership, UniProt disease involvement and Open Targets association scores, residue-level contacts of each methylarginine in PDB structures, and an optional second ring of PRMTs and Tudor proteins bound by the protein's partners; dataset-level statistics are shown separately. The other-PTM sites (from the cross-PTM master table) are labelled as non-exhaustive. A switch at the top of each protein page flips all context cards between protein-specific evidence and dataset context. The overview condensate chart can count only CD-CODE memberships at a chosen confidence (≥ 3 of 5 by default), and the interface card adds a residue-environment analysis (Levy core/rim/support, base vs backbone hydrogen bonds with RNA/DNA) and an RG/RGG-context test against same-protein controls. Sites modelled as methylated arginine in PDB entries are flagged and filterable, and contact entries open the structure at the site. Putative contacts from predicted complexes come from ProtVar (methylarginines within 5 Å of the partner in AlphaFold-Multimer models of human pairs, pDockQ ≥ 0.23): sites are tagged and filterable, the overview repeats the within-protein interface test on these models (including disordered regions), and protein pages list the partners with the model and AlphaFold DB complexes (live API); either model opens in the 3D view, which measures each methylarginine's distance to the partner. The ProtVar API serves a later prediction run than its 2024 5 Å interface table, so a site the table places at the interface can sit further away in the viewable model; the viewer says so when that happens, and pairs without a served model show no View button. UniProt, AlphaFold DB, RCSB PDB and PDBe data are fetched live in the browser. Overview charts export as SVG, PNG or TSV. Styling follows shechterlab.org.

The overview leads with the functional ontology and disease families, then the interactions findings; the Documentation tab holds the methods behind every chart (each chart links to its section). The Figures tab builds figures from the current search and filters: the selected sites aligned at the methylarginine with a frequency or information-content sequence logo (SVG, PNG or FASTA), and the selected proteins on one residue scale with UniProt domains and disordered regions. Protein pages add a conservation card that loads Ensembl Compara orthologue alignments (mouse to yeast, live) and shows the residue aligned to each methylarginine, with an exportable alignment around any site. The Interactions tab opens with an evidence map (an UpSet table of how RNA crosslinks, RBDmap fragments, DNA crosslinks, PDB protein and RNA/DNA interfaces, predicted interfaces and confident predicted contacts combine on the selected sites; totals and combinations are clickable filters, shareable as `evAll`/`evExact` in the URL) with the proteins carrying most sites of a chosen type, and complex rows expand to their members. Whenever a search or filter is set, the Proteins and Sites views show what sets the selection apart: site features over- or under-represented against the other sites (Fisher, BH), each clickable to narrow the selection. Every exported figure, including 3D snapshots and printed pages, ends with a credit line naming the site, its address and DeAngelo, Silverstein & Shechter (2026). The site's name is read from `<meta name="application-name">` in `docs/index.html` (change it there, in `<title>` and in the page heading to rename the site) and its address from the browser, so credits follow a rename or a move.

The **Interactions** tab gathers the site-wide interaction analyses. It opens with key findings written from the exported statistics (so the wording follows the numbers when the data change), then one collapsible section each for RNA/DNA crosslinks in cells, PDB interfaces, predicted complexes, partners, complexes and Gene Ontology; the overview keeps a short summary that links there. Partner filters cover PRMTs, Tudor-domain (InterPro IPR002999), WD40-repeat (IPR001680) and PHD-finger (IPR001965, IPR019787) proteins, other methylated proteins and "any other protein" (a partner in none of these classes), and the specific-partner box takes any gene symbol. Residues that touch another copy of the same protein in a PDB structure are labelled "self-contact" (homodimer or homo-oligomer). The sequence card exports an SVG or PNG figure titled with the gene, protein name and UniProt accession, with a legend.

Regenerate the data after rerunning the pipeline:

```bash
. hpc/ptm_env.sh
python scripts/36_export_site_explorer_data.py \
  --results-root "$PTM_RESULTS_ROOT" \
  --uniprot-metadata "$PTM_UNIPROT_METADATA" \
  --outdir docs/data
```

This writes `docs/data/methylarg_sites.json` (search and filter index) and `docs/data/methylarg_evidence.json` (per-row source evidence and disease terms, loaded when a protein is opened). Optional inputs that are missing are skipped with a note.

Then export the interaction context (downloads are cached in `data/interactions_cache/`, which Git ignores; the first run makes about 5,000 PDBe requests and takes roughly 25 minutes, later runs reuse the cache, and `--refresh` downloads everything again):

```bash
python scripts/37_export_interaction_data.py
```

This writes `docs/data/methylarg_interactions.json` from STRING v12.5 physical links (CC BY 4.0), Complex Portal curated and predicted complexes (CC0), the UniProt reviewed human proteome and QuickGO GO annotations (CC BY 4.0), PDBe-KB interface residues and SIFTS observed segments and chain mapping (EMBL-EBI), UniProt disordered-region and domain features, CD-CODE condensate membership with per-protein confidence (cd-code.org API; licence not stated on the site) and residue-level contacts from the RCSB ModelServer (PDB data, CC0), RCSB-flagged PDB entries that model a methylated arginine (2MR, DA2, NMM, AGM; PDBe modified-residue API), and ProtVar predicted interfaces (`2024.05.28_interface_summary_5A.tsv.gz`, plus `af2complexes_interfaces_v1.tsv.gz` to record which pairs the ProtVar API serves a model for; EMBL-EBI, CC BY 4.0). The contact stage makes about 7,000 ModelServer requests on its first run (run it alone with `--fetch-contacts-only`), and the interface-geometry stage about 8,000 more for 10 Å residue environments (`--fetch-environments-only`).

The predicted-contact gallery is part of the same export: it measures every methylarginine within 5 Å of the partner in served ProtVar models (pDockQ ≥ 0.5) and in AlphaFold DB human heterodimers (salt bridges ≤ 4 Å, hydrogen bonds ≤ 3.5 Å, cation–π ≤ 5 Å, clashes < 2.2 Å, pLDDT of both residues). It downloads a few thousand model files on its first run and keeps only the per-site results (`data/interactions_cache/protvar_measure`, `afdb_measure`); run it alone with `--fetch-models-only`.

The same export adds residue-level RNA and DNA contact in living human cells:
- **RNA-crosslinked residues:** from the pepR-MS supplement (Sha, Kuster & Trendel, bioRxiv doi:10.1101/2025.09.05.674461, CC BY-NC). Its Tables S2–S3 re-analyse Kramer 2014, XRNAX 2019, RBS-ID 2020 and pRBS-ID 2021 with pepR-MS.
- **RNA-bound fragments:** RBDmap (Castello et al. 2016, Oxford Research Archive).
- **DNA-crosslinked peptides and residues, and the DNA interactome:** Trendel et al. 2025, *Cell* (doi:10.1016/j.cell.2025.04.037).
- **DNA-binding regions:** UniProt.

bioRxiv often refuses scripted downloads (HTTP 429). If it does, save the pepR-MS supplementary zip as `data/interactions_cache/crosslinks/pepr_media-1.zip`.

**Licence.** Code is MIT, and data, figures and text are CC BY-NC-SA 4.0; see [LICENSE.md](LICENSE.md). The share-alike, non-commercial terms follow from PhosphoSitePlus, DisGeNET and pepR-MS.

**Data licences.** The Documentation page of the explorer lists every source with its version, licence and citation. Two sources carry non-commercial terms: PhosphoSitePlus® sites (8,175 sites in the union, via Maron et al. 2021 Table S5) are licensed under CC BY-NC-SA 3.0 (cite Hornbeck et al. 2015, doi:10.1093/nar/gku1267, and acknowledge www.phosphosite.org), and DisGeNET-derived disease terms are CC BY-NC-SA 4.0. pepR-MS crosslinks are CC BY-NC. Downloads that include these data carry the same terms.

Gene symbols follow current UniProt; symbols reported by the source tables (some outdated or belonging to a paralog) are kept as search aliases. If the metadata TSV also has UniProt's `Gene Names` column, its synonyms become aliases too. The committed data used the same query as `00_stage_context_sources.py` with `gene_names` added to the fields:

```bash
curl -G https://rest.uniprot.org/uniprotkb/stream \
  --data-urlencode "format=tsv" \
  --data-urlencode "query=(proteome:UP000005640) AND (reviewed:true)" \
  --data-urlencode "fields=accession,id,gene_primary,gene_names,protein_name,length" \
  -o uniprot_human_reviewed_metadata_with_synonyms.tsv
```

Preview locally (the page fetches its JSON, so open it through a server rather than from disk):

```bash
cd docs && python3 -m http.server 8000
# then open http://localhost:8000
```

Deploy either way; neither needs a build step:

- **GitHub Pages:** Settings → Pages → Build and deployment → Source: *Deploy from a branch*, Branch: `main`, Folder: `/docs`. The site appears at `https://shechterlab.github.io/PTM_Pipeline/`.
- **Netlify:** Add new site → Import an existing project → this repository. `netlify.toml` already sets the publish directory to `docs` with no build command.

## Current Interpretive Guardrails

- The integrated union is the primary analysis universe; source support is retained as provenance/audit metadata rather than used as a primary confidence-tier scheme.
- Methylarginine state labels are useful where source-supported, but the full union should not be overinterpreted as a clean `Rme1` versus `Rme2a` versus `Rme2s` state-resolved dataset.
- Domain adjacency is descriptive in this dataset. The IDR-length-conditioned null shows that simple domain-edge proximity is largely explained by IDR geometry/length.
- The stronger architectural signal is short-range methylarginine clustering/local density relative to same-protein arginine positions.
- Remapping (`scripts/02`) was corrected in four ways; rerun the pipeline from `scripts/02` (`hpc/slurm_rerun_remap_fix.sh`) to regenerate the union and every downstream result.
  1. Sources give either centred site windows (Guo 13 residues; PhosphoSitePlus, Musiani, Fong and ProMetheusDB 15; Larsen and Wei 31) or peptides that contain the site somewhere (Maron's own data, Fedoriw, Li, Lim). Earlier runs read every token as a centred window and, when its middle residue was not an arginine, took the arginine nearest the middle. For peptides that ignored the reported position: 1,773 rows whose reported arginine lies inside the peptide in the current sequence (1,767 on proteins whose sequence has not changed since the source was published) were moved to another arginine, and multi-site peptides collapsed onto one site (SNRPD1 R104–R119 onto R98–R106; TAF15 R467, R483 and R490 from one Fedoriw peptide onto R475). Which kind a source gives is now read from the data (a source whose tokens mostly share one window length gives windows), the reported position of a peptide is kept when it is an arginine inside the peptide (`sequence_position_verified`), a peptide places a site by itself only when it has a single arginine, and a window cut short is read from the protein end it touches (Guo pads such windows with X and ProMetheusDB with _, which names the end). A short window that lies inside a longer entry (the source protein began or ended there: another start codon, the mature histone chain, the ubiquitin C-terminus of RPS27A/UBA52, an isoform) gives the reading that lands on an arginine; when both readings do, only one the reported position names is taken.
  2. When the sequence matched a different accession than the row's own (usually the reviewed entry behind an unreviewed or mislabelled source accession), the position was computed on the matched protein but filed under the row's accession. It is now filed under the matched accession (`is_accession_remapped`). Earlier unions carried 857 rows on unreviewed accessions in reviewed-entry numbering (578 of those sites were already in the union under the reviewed entry, so they were counted twice) and 49 rows on the wrong reviewed protein (PTCD1 R235 filed as MED12 R235).
  3. A full-length window whose middle residue is not an arginine gives no site. Maron Table S5 carries 268 Musiani windows centred on a lysine at the reported position; earlier runs turned 131 of them into arginine sites a few residues away.
  4. The reported arginine stands on its own only for rows without sequence evidence (the master table's position-only records, `same_accession_position_verified`). A row whose peptide or window lies elsewhere in the entry is not placed (`unresolved_evidence_elsewhere`); one whose evidence is not in the entry must match it with its site on the reported arginine at no more than 2 residues and no more than a fifth of them (variants, sequence conflicts; `same_accession_evidence_near_match`), or it is not placed (`unresolved_evidence_mismatch`). Earlier code kept the reported arginine whenever it was one, which placed 90 rows against their own evidence, among them isoform-only sites filed on unrelated canonical arginines (FMR1 R562/R564, PML R632/R717). When two entries other than the row's own support a site equally, the row stays ambiguous instead of going to whichever was listed first.
  Replaying the remapping on the committed source rows with current UniProt sequences, the old code reproduces the published remapping on 28,958 of 28,960 rows (the two others are RAB12, whose sequence UniProt extended in September 2026) and the corrected code gives 13,324 sites on 4,706 reviewed proteins: 12,481 sites kept, 843 added and 1,238 removed, 81 reviewed proteins entering and 35 leaving (UniProt release 2026_03; the HPC rerun uses the staged FASTA of the published run). Of Maron's own 1,829 sites (630 of them new in that paper) the published union held 1,455 (288 of the new ones) at the reported position and the corrected union holds 1,812 (620); Maron Table S3 SDMA-IP sites rise from 372 to 454 of 458.
- To check these by hand: `python scripts/qc_remap_hand_check.py examples` prints, for SNRPD1, ACO2 (A2A274), MED12/PTCD1, TTC28 and MRE11, each source row's reported position, the residue UniProt has there, where its peptide or window occurs, and where the pipeline placed the site; any gene or accession works in place of `examples`. `python scripts/qc_remap_hand_check.py --audit-union` checks that every union site is an arginine in current UniProt: in the published union 13,024 of 13,719 are, 233 are not or fall beyond the sequence end, and 462 sit on accessions UniProt no longer carries; in the corrected union all 13,324 are. With `--fasta` (the run's own FASTA) the audit exits 1 when a site is not an arginine; it cannot see a site placed on the wrong arginine, which is what the evidence check is for.
- To check every row: `python scripts/qc_remap_evidence_check.py --results-root PTM_results --fasta <reviewed FASTA>` derives, with code independent of `scripts/02`, the site each source row's own window, peptide or position supports, compares it with the placed site, and says whether the reported position still fits the current sequence (drift class A–E; `--extra-fasta` adds unreviewed and isoform sequences). On the published run 24,838 of 28,960 rows agree with their evidence, 2,142 contradict it, 238 are placed without evidence and 901 are filed under accessions outside the reviewed FASTA; on the corrected code 27,828 agree and 22 contradict it, all Larsen or Fedoriw sequences that match a paralog exactly and the reported protein within 2 residues (issue 9 of `DATA_PROVENANCE_REVIEW.md`, listed with reasons in `config/remap_evidence_allowlist.tsv`). With `--fail-on-contradictions` the check exits 1 on any row placed against or without its evidence that the `--allow` list does not name. 90% of rows (26,204) have a reported position that still fits the current sequence; most of the rest are Larsen 2016 windows on sequences or accessions UniProt has since changed, which both versions remap from the window.
- Per-protein arginine odds (`scripts/08`) now group sites by accession only. They were grouped by accession and gene symbol, which split 138 proteins whose sites carry an alias (MRE11/MRE11A, 13 sites tested as 7 and 6) into separate, weaker tests; on the published union this changes the proteins at q ≤ 0.05 from 807 to 832.
- The stars on the positional heatmaps (`scripts/08`, `arg_methyl_positional_enrichment_*`) now use BH q-values (*q < 0.05, **q < 1e-5, ***q < 1e-20), as the Figure 3A legend of DeAngelo et al. (2026) states and as the printed panel shows. They used raw p-values, which differ from the printed panel at one cell (P +1: p = 0.048, q = 0.071).
- Both checks are pipeline steps: `hpc/slurm_run_integration.sh` runs them right after the union is built and stops the run if either fails. `python -m pytest tests/` (or `python tests/test_remap_edge_cases.py`) runs the remapping edge cases and the gates on small hand-made sequences.
- Provenance: `scripts/00_stage_context_sources.py` writes `uniprot_release.json` next to the staged FASTA (release and date from UniProt's response headers, checksums; a FASTA staged earlier gets its checksum and an unknown release). `scripts/36_export_site_explorer_data.py --canonical-fasta ... --dataset-version ...` puts the dataset version, pipeline commit, UniProt release, input and FASTA checksums and the gate results in `meta.provenance`, and each protein's reference-sequence MD5 (UniProt's `sequence.md5` format) in the protein rows. The explorer compares that checksum with the live UniProt entry and, when they differ, does not draw current UniProt positional annotations against the sites; `--strict` stops the export when an input table or a reference sequence is missing.
- Motif and clustering modules can analyze fewer sites than the integrated union because they require sequence-confirmed canonical residues in the staged FASTA. Use each module's audit TSV when reporting counts.

## Hygiene

Before publishing or archiving, run:

```bash
bash hpc/check_repo_hygiene.sh "$PWD"
git status --short --ignored
```

The hygiene check fails on merge-conflict markers, tracked local app artifacts, tracked bytecode, tracked archives, tracked staged data/results, or tracked files larger than 10 MB.

Ignored local paths include `PTM_data/`, `PTM_results/`, bytecode caches, local archives, editor artifacts, and cluster logs. The tracked repository should contain only source code, lightweight configuration, environment files, and documentation.

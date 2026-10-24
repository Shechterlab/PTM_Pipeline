# PTM_Pipeline: human methylarginine sites

The human protein arginine methylation sites of DeAngelo, Silverstein & Shechter, *Protein arginine methylation
rewires molecular interactions*, Cell Chem Biol 33:1404–1421 (2026),
[doi:10.1016/j.chembiol.2026.09.011](https://doi.org/10.1016/j.chembiol.2026.09.011): the data, the code that builds
them, and the [Methylarginine Site Explorer](https://shechterlab.github.io/PTM_Pipeline/) website.

## Versions

| Version | Sites | Data | Website |
|---|---|---|---|
| **v1.3.0** (current) | 11,601 sites on 4,432 proteins, and 2,406 lower-confidence sites | [`releases/v1.3.0/`](releases/v1.3.0/) ([what changed](releases/v1.3.0/CHANGES.md)) | [`docs/`](docs/) |
| v1.2.0 | 11,601 sites on 4,432 proteins, and 2,406 lower-confidence sites | [`releases/v1.2.0/`](releases/v1.2.0/) ([what changed](releases/v1.2.0/CHANGES.md)) | [`docs/v1.2.0/`](docs/v1.2.0/) |
| v1.1.0 | 13,365 sites on 4,723 proteins, and 6 pending | [`releases/v1.1.0/`](releases/v1.1.0/) ([what changed](releases/v1.1.0/CHANGES.md)) | [`docs/v1.1.0/`](docs/v1.1.0/) |
| v1.0.0 (published with the paper) | 13,719 sites on 5,090 accessions | [`releases/v1.0.0/`](releases/v1.0.0/) | [`docs/v1.0.0/`](docs/v1.0.0/) |

v1.3.0 keeps the sites of v1.2.0 and adds what methylated arginines contact: RNA and DNA in PDB structures (with the bound sequences), RNA-binding evidence from living cells (crosslinked residues, RBDmap, CAPRI, pCLAP and RBR-ID, each also tested among the arginines its experiment measured), and what their proteins bind in vitro and in cells (RNA Bind-n-Seq, HTR-SELEX, JASPAR, ENCODE and ENCORE eCLIP); the website opens each page with a picture and its answer. v1.2.0 reads every study from its own tables and keeps as a site only an arginine that one of its records passes five tests for: the record's peptide or window covers it, in a protein the record names, at one place in the reviewed proteome (or once in each of several paralogs, when the site goes on the copy most records name and is flagged `paralog_assigned`), the methyl group can only be on an arginine, and the study localises it there. The 2,406 arginines that a study reports, or an earlier release had, and that no record passes every test for are kept as lower-confidence sites with their test results; the website shows them on request, and a second run of the analyses counts the 1,844 that are surely methylarginine (`statistics_with_lower_confidence`). Every release's analyses, rerun with the same code, are in its `analyses/` folder, and [`releases/comparison/`](releases/comparison/) compares them. v1.1.0 uses the sources and methods of the paper, mapped to UniProt release 2026_03, with corrections to how some source records were placed on UniProt sequences. The repository as it was before this reorganisation, with the full
analysis output of the published run, is the branch
[`v1.0.0_archive`](https://github.com/Shechterlab/PTM_Pipeline/tree/v1.0.0_archive); the tag `v1.0.0` and the branch
`release/v1.0.0` mark the published run. [CONTRIBUTING.md](CONTRIBUTING.md) gives the branches (`main` is what the
website serves, `develop` the working branch, `release/vX.Y.Z` each published version, frozen) and how a change becomes
a release.

Each release folder has a README that lists its tables and every column, and checksums; from v1.1.0 it also has a data
dictionary and a record of what changed since the release before.

## Layout

| Folder | Holds |
|---|---|
| [`releases/`](releases/) | the site union, one folder per version |
| [`docs/`](docs/) | the Methylarginine Site Explorer (GitHub Pages); `docs/v1.0.0/` is the published version and `docs/v1.1.0/` and `docs/v1.2.0/` the releases since, unchanged |
| [`inputs/`](inputs/) | every file a release and the website are built from, with `MANIFEST.tsv` |
| [`provenance/`](provenance/) | where the sites come from, in one page |
| [`scripts/`](scripts/) | the pipeline, numbered by step |
| [`config/`](config/) | the functional and domain-class ontologies, and the accepted paralog placements |
| [`hpc/`](hpc/) | the release build, and the cluster scripts of the published run |
| [`envs/`](envs/) | the conda environment, and the resolved Python packages of the latest build |

## Data staging

Staging and building are separate steps. Staging downloads the sources and writes the tables a build reads into
[`inputs/`](inputs/) (`scripts/00_*`, `37a`, `43*`); [`inputs/README.md`](inputs/README.md) says how each folder was
staged and with which script, and [`inputs/MANIFEST.tsv`](inputs/MANIFEST.tsv) gives every file's source, version or
download date, licence, MD5, how it was made, and whether it is redistributed: `git` (in this repository, about 410 MB),
`zenodo` (the raw downloads behind the tables, journal supplements, bulk files and API answers, in the Zenodo archive of
the release) or `checksum_only` (kept out because of its terms; the manifest has its MD5 and source). A build reads only
the `git` files and checks each against its MD5, with two exceptions it derives when they are absent:

- `inputs/interactions/pclap_mullari2017_peptides.tsv`, from the MaxQuant output of PRIDE
  [PXD005739](https://www.ebi.ac.uk/pride/archive/projects/PXD005739) (`scripts/43d_stage_pclap.py`);
- `inputs/interactions/rbdmap_castello2016_measured_peptides.tsv`, from the peptide and quantification tables of PRIDE
  [PXD000883](https://www.ebi.ac.uk/pride/archive/projects/PXD000883) (`scripts/43f_*`).

Both are PRIDE datasets under the EMBL-EBI terms of use without a licence from their submitters, so they are not
redistributed. A derived table must match its MD5 in the manifest, or the build stops.

## Building a release

```bash
conda env create -n ptm_pipeline -f envs/ptm_pipeline_hpc.yml
conda activate ptm_pipeline
PREVIOUS=releases/v1.3.0 bash hpc/build_release.sh v1.3.1 WORK_ROOT
```

The environment file gives version ranges. The exact Python packages of the latest build are in
[`envs/cloud_build_requirements.txt`](envs/cloud_build_requirements.txt), and each release's `release.json` records the
Python and package versions that wrote it (`python_environment`) and the R and package versions of its enrichment tests
(`enrichment_tools`).

The build checks every input against `inputs/MANIFEST.tsv`, reads each study from its own table in `inputs/sources/`,
places each record on the reference sequences, keeps as sites only the records that pass the five site tests (the
paralog rule included), keeps the arginines that fail as lower-confidence sites, stops if a record is placed against
its own evidence or a site is not an arginine, runs the analyses twice (on the sites, and with the lower-confidence
sites that are surely methylarginine), and writes `releases/VERSION/` (with the changes traced from `PREVIOUS`
and both analyses' result tables) and the website's data in `docs/data/` (the view with lower-confidence sites as
`methylarg_*_all.json`). It writes the checksums last and checks them (`md5sum -c MD5SUMS` in the release folder and
in each analyses folder). v1.1.0, which placed records without localisation, is
rebuilt with `PREVIOUS=releases/v1.0.0 STUDIES='' LOCALISED_ONLY=0 bash hpc/build_release.sh v1.1.0`; the header of
`hpc/build_release.sh` lists the options. It takes about three hours on one core.

What a build needs besides the repository:

- **R** (clusterProfiler and DOSE, in the environment file) for the tests of protein and gene sets. A work root that
  already holds their results for the same inputs, settings and R script (`scripts/enrichment.py` checks a SHA-256 of
  all three) reuses them, so a rebuild that changes nothing those tests read runs without R.
- **The network** only to derive the two PRIDE tables above, and when `inputs/interactions/` lacks interaction data for
  a site (a protein or structure the staged tables do not cover): `scripts/37` then fetches it into
  `WORK_ROOT/interactions_topup` and the build stops before writing the release, with the command that stages it
  (`scripts/37a_stage_interactions.py --topup`). With both PRIDE tables in place and nothing to top up, a build runs
  offline.
- **Storage**: about 0.7 GB for the work root, and about 0.2 GB of downloads in `PTM_data/external/` when the PRIDE
  tables are derived (the pCLAP MaxQuant output is 193 MB; the four RBDmap tables 19 MB).

## Scripts

| Scripts | Step |
|---|---|
| `00_*` | staging inputs: UniProt, InterPro, MobiDB, PANTHER, the disease tables |
| `01_*` | parsing the source tables |
| `02_*` | placing records on reference sequences, the union, the sites kept pending |
| `qc_*`, `38_*` | the evidence check, the residue audit, and comparing two results roots |
| `03_*` to `16_*` | analyses: functional classes (the paper's and PANTHER's), domains, disorder, motifs, other PTMs, clustering, condensates, disease, Maron et al. 2021 summaries |
| `17_*` to `35_*` | extended analyses and the paper's figures |
| `36_*`, `37_*` | the website's data: sites, evidence and interactions (`37a` stages the interaction inputs) |
| `39_*` to `42_*` | the website's annotation layers (CARM1 knockout, Human Protein Atlas; CPLM 4.0 for v1.0.0) and the lysine and arginine methylation co-occurrence |
| `enrichment.py`, `enrich_sets.R` | the over-representation tests of protein and gene sets (clusterProfiler `enricher`, DOSE `enrichDGN`), which the analyses call |
| `build_release.py`, `release_*.py` | writing a release folder |

`hpc/slurm_run_integration.sh` and `scripts/run_extended_analyses.sh` are the workflow of the published run with its
figures; they read staged files from `PTM_data/` and write `PTM_results/` (paths in `hpc/ptm_env.sh`).

## Cite and licence

Cite the paper and the studies that reported the sites you use (each release's studies table lists them). Code is MIT;
data and content are CC BY-NC-SA 4.0; third-party data keep their own terms, and PhosphoSitePlus®, DisGeNET, CD-CODE,
the pepR-MS compendium and the Trendel et al. 2025 DNA crosslinks are non-commercial. See [LICENSE.md](LICENSE.md).

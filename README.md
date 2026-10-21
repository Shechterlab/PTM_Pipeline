# PTM_Pipeline: human methylarginine sites

The human protein arginine methylation sites of DeAngelo, Silverstein & Shechter, *Protein arginine methylation
rewires molecular interactions*, Cell Chem Biol 33:1404–1421 (2026),
[doi:10.1016/j.chembiol.2026.09.011](https://doi.org/10.1016/j.chembiol.2026.09.011): the data, the code that builds
them, and the [Methylarginine Site Explorer](https://shechterlab.github.io/PTM_Pipeline/) website.

## Versions

| Version | Sites | Data | Website |
|---|---|---|---|
| **v1.1.0** (current) | 13,365 sites on 4,723 proteins, and 6 pending | [`releases/v1.1.0/`](releases/v1.1.0/) ([what changed](releases/v1.1.0/CHANGES.md)) | [`docs/`](docs/) |
| v1.0.0 (published with the paper) | 13,719 sites on 5,090 accessions | [`releases/v1.0.0/`](releases/v1.0.0/) | [`docs/v1.0.0/`](docs/v1.0.0/) |

v1.1.0 uses the sources and methods of the paper, mapped to UniProt release 2026_03, with corrections to how some
source records were placed on UniProt sequences. The repository as it was before this reorganisation, with the full
analysis output of the published run, is the branch
[`v1.0.0_archive`](https://github.com/Shechterlab/PTM_Pipeline/tree/v1.0.0_archive) (and tag `v1.0.0`).

Each release folder has a README that lists its tables and every column, and checksums; from v1.1.0 it also has a data
dictionary and a record of what changed since the release before.

## Layout

| Folder | Holds |
|---|---|
| [`releases/`](releases/) | the site union, one folder per version |
| [`docs/`](docs/) | the Methylarginine Site Explorer (GitHub Pages); `docs/v1.0.0/` is the published version, unchanged |
| [`inputs/`](inputs/) | every file a release and the website are built from, with `MANIFEST.tsv` |
| [`provenance/`](provenance/) | where the sites come from, in one page |
| [`scripts/`](scripts/) | the pipeline, numbered by step |
| [`config/`](config/) | the functional and domain-class ontologies, and the accepted paralog placements |
| [`hpc/`](hpc/) | the release build, and the cluster scripts of the published run |
| [`envs/`](envs/) | the conda environment |

## Data staging

Everything a build reads is staged once into [`inputs/`](inputs/); a build downloads nothing.
[`inputs/README.md`](inputs/README.md) says how each folder was staged and with which script, and
[`inputs/MANIFEST.tsv`](inputs/MANIFEST.tsv) gives every file's source, version or download date, licence, MD5 and how
it was made. The raw downloads behind the tables (journal supplements, bulk files, API answers) go to the Zenodo
archive of the release.

## Building a release

```bash
conda env create -n ptm_pipeline -f envs/ptm_pipeline_hpc.yml
conda activate ptm_pipeline
bash hpc/build_release.sh v1.2.0
```

The build checks every input against `inputs/MANIFEST.tsv`, reads each study from its own table in `inputs/sources/`,
places each record on the reference sequences, keeps as sites only the records whose source localises the methyl group
to one arginine that the peptide names once in the reviewed proteome, stops if a record is placed against its own
evidence or a site is not an arginine, runs the analyses, and writes `releases/v1.2.0/` (with the changes traced from
`releases/v1.1.0/`) and the website's data in `docs/data/`. v1.1.0, which placed records without localisation, is
rebuilt with `PREVIOUS=releases/v1.0.0 STUDIES='' LOCALISED_ONLY=0 bash hpc/build_release.sh v1.1.0`; the header of
`hpc/build_release.sh` lists the options. It takes about half an hour on one core. The tests of
protein and gene sets run in R (clusterProfiler and DOSE, in the environment file), so `Rscript` must be on the path.

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

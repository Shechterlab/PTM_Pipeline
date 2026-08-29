# Methylarginine Interaction-Regulation Pipeline

This repository contains the analysis code associated with DeAngelo, Silverstein, and Shechter, "Molecular Interaction Regulation by Arginine Methylation", invited for *Cell Chemical Biology* (2026).

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
| `scripts/33_*` | Export helper: the non-methyl PTM slice that falls inside an external flank table. |
| `scripts/common.py`, `scripts/clustering_utils.py` | Shared utilities. |
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

For a local, stepwise rerun, use the commands in [`hpc/slurm_run_integration.sh`](hpc/slurm_run_integration.sh) as the authoritative ordered recipe. The scripts themselves expose explicit `--help` arguments and write tables/figures under `PTM_RESULTS_ROOT`.

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

## Current Interpretive Guardrails

- The integrated union is the primary analysis universe; source support is retained as provenance/audit metadata rather than used as a primary confidence-tier scheme.
- Methylarginine state labels are useful where source-supported, but the full union should not be overinterpreted as a clean `Rme1` versus `Rme2a` versus `Rme2s` state-resolved dataset.
- Domain adjacency is descriptive in this dataset. The IDR-length-conditioned null shows that simple domain-edge proximity is largely explained by IDR geometry/length.
- The stronger architectural signal is short-range methylarginine clustering/local density relative to same-protein arginine positions.
- Motif and clustering modules can analyze fewer sites than the integrated union because they require sequence-confirmed canonical residues in the staged FASTA. Use each module's audit TSV when reporting counts.

## Hygiene

Before publishing or archiving, run:

```bash
bash hpc/check_repo_hygiene.sh "$PWD"
git status --short --ignored
```

The hygiene check fails on merge-conflict markers, tracked local app artifacts, tracked bytecode, tracked archives, tracked staged data/results, or tracked files larger than 10 MB.

Ignored local paths include `PTM_data/`, `PTM_results/`, bytecode caches, local archives, editor artifacts, and cluster logs. The tracked repository should contain only source code, lightweight configuration, environment files, and documentation.

## Exports for downstream projects

The cross-PTM union assembled here carries phosphorylation, acylation,
ubiquitin/SUMO and lysine methylation with enzyme attribution, which is useful
to projects working on a different coordinate universe. `scripts/33_*` writes
the slice of that union falling inside an externally supplied set of protein
windows, so the consumer needs neither this repository's staged inputs nor a
second copy of the source databases.

```bash
python scripts/33_export_flank_ptm_annotation.py \
    --flanks /path/to/flank_table.tsv \
    --output PTM_results/exports/ptm_sites_in_domain_flanks.tsv
```

The flank table must carry `canonical_UniProtAC`, `flank_start` and
`flank_end`. Arginine methylation is deliberately not exported: a consumer
studying methylarginine already has it, and re-exporting would put the same
evidence in two places under different provenance. CK2 and the proline-directed
kinases are labelled from the attributed enzyme, never from sequence.

The current consumer is the MicroDomains repository, whose acidic
domain-tethered IDR analysis uses the phosphorylation slice to ask which kinase
family tracks acidic domain-adjacent segments.

#!/bin/bash
#SBATCH -J ptm_integrate
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 16:00:00

CONDA_BASE=/gs/gsfs0/hpc01/rhel8/apps/conda3
CONDA_EXE="$CONDA_BASE/bin/conda"

if [ -x "$CONDA_EXE" ]; then
  eval "$("$CONDA_EXE" shell.bash hook)"
  CONDA_HOOK_STATUS=$?
else
  . "$CONDA_BASE/bin/activate"
  CONDA_HOOK_STATUS=$?
fi

ENV_NAME=${PTM_ENV_NAME:-ptm_pipeline}
ENV_PREFIX=${PTM_ENV_PREFIX:-}
INTERPRO_MODE=${PTM_INTERPRO_MODE:-auto}
INTERPRO_API_WORKERS=${PTM_INTERPRO_API_WORKERS:-8}
INTERPRO_API_PAGE_SIZE=${PTM_INTERPRO_API_PAGE_SIZE:-200}
INTERPRO_API_TIMEOUT=${PTM_INTERPRO_API_TIMEOUT:-120}
INTERPRO_API_RETRIES=${PTM_INTERPRO_API_RETRIES:-5}

run_step() {
  echo "[ptm_pipeline] $1"
  shift
  "$@"
  status=$?
  if [ "$status" -ne 0 ]; then
    echo "[ptm_pipeline] step failed with exit code $status"
    echo "[ptm_pipeline] command: $*"
    exit "$status"
  fi
}

if [ "$CONDA_HOOK_STATUS" -ne 0 ]; then
  echo "[ptm_pipeline] failed to initialize Conda shell hook"
  exit "$CONDA_HOOK_STATUS"
fi

if [ -n "$ENV_PREFIX" ]; then
  conda activate "$ENV_PREFIX"
  ACTIVATE_STATUS=$?
else
  conda activate "$ENV_NAME"
  ACTIVATE_STATUS=$?
fi

if [ "$ACTIVATE_STATUS" -ne 0 ]; then
  echo "[ptm_pipeline] failed to activate Conda environment"
  echo "[ptm_pipeline] create one with: bash hpc/create_ptm_env.sh envs/ptm_pipeline_hpc.yml $ENV_NAME"
  exit "$ACTIVATE_STATUS"
fi

BASE_DIR=${1:-$PWD}
export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib-${SLURM_JOB_ID:-manual}"
mkdir -p "$MPLCONFIGDIR"
mkdir -p "$BASE_DIR/data/context"

run_step "stage UniProt context" python "$BASE_DIR/scripts/00_stage_context_sources.py" \
  --outdir "$BASE_DIR/data/context" \
  --skip-interpro
run_step "parse Maron S5" python "$BASE_DIR/scripts/01_parse_maron_s5.py" \
  --input "$BASE_DIR/data/external/Table_S5_Compiled_Methylarginine_Data.xlsx" \
  --outdir "$BASE_DIR/results/parsed_maron"
run_step "parse ProMetheus mmc2" python "$BASE_DIR/scripts/01_parse_prometheus_mmc2.py" \
  --input "$BASE_DIR/data/external/mmc2.xlsx" \
  --outdir "$BASE_DIR/results/parsed_prometheus"
run_step "remap to canonical coordinates" python "$BASE_DIR/scripts/02_remap_arg_methyl_to_canonical.py" \
  --base-master "$BASE_DIR/data/base/human_ptm_master.tsv" \
  --maron "$BASE_DIR/results/parsed_maron/maron_s5_arg_methyl_sites.tsv" \
  --prometheus "$BASE_DIR/results/parsed_prometheus/prometheus_mmc2_arg_methyl_sites.tsv" \
  --canonical-fasta "$BASE_DIR/data/context/uniprot_human_reviewed_canonical.fasta" \
  --uniprot-metadata "$BASE_DIR/data/context/uniprot_human_reviewed_canonical_metadata.tsv" \
  --outdir "$BASE_DIR/results/remapped"
run_step "integrate remapped sites" python "$BASE_DIR/scripts/02_integrate_arg_methyl_sources.py" \
  --premerged-remapped "$BASE_DIR/results/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv" \
  --outdir "$BASE_DIR/results/integrated"
run_step "functional enrichment" python "$BASE_DIR/scripts/03_functional_class_enrichment.py" \
  --base-master "$BASE_DIR/data/base/human_ptm_master.tsv" \
  --integrated-sites "$BASE_DIR/results/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --ontology "$BASE_DIR/config/functional_ontology.json" \
  --outdir "$BASE_DIR/results/functional_class_union"
run_step "example figures" python "$BASE_DIR/scripts/04_example_figures.py" \
  --integrated-sites "$BASE_DIR/results/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --outdir "$BASE_DIR/results/example_figures"
run_step "motif enrichment and Arg odds" python "$BASE_DIR/scripts/08_motif_and_arg_odds.py" \
  --integrated-sites "$BASE_DIR/results/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --canonical-fasta "$BASE_DIR/data/context/uniprot_human_reviewed_canonical.fasta" \
  --outdir "$BASE_DIR/results/motif_arg_odds"
run_step "methyl-Arg clustering" python "$BASE_DIR/scripts/10_methyl_arg_clustering.py" \
  --integrated-sites "$BASE_DIR/results/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --canonical-fasta "$BASE_DIR/data/context/uniprot_human_reviewed_canonical.fasta" \
  --outdir "$BASE_DIR/results/methyl_arg_clustering"
run_step "cross-PTM clustering" python "$BASE_DIR/scripts/11_compare_ptm_clustering.py" \
  --base-master "$BASE_DIR/data/base/human_ptm_master.tsv" \
  --integrated-arg-sites "$BASE_DIR/results/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --canonical-fasta "$BASE_DIR/data/context/uniprot_human_reviewed_canonical.fasta" \
  --outdir "$BASE_DIR/results/ptm_clustering"

if [[ -f "$BASE_DIR/data/context/mobidb_human_reviewed_disorder_intervals.tsv" ]]; then
  run_step "annotate disorder context" python "$BASE_DIR/scripts/06_annotate_disorder_context.py" \
    --sites "$BASE_DIR/results/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
    --disorder-intervals "$BASE_DIR/data/context/mobidb_human_reviewed_disorder_intervals.tsv" \
    --canonical-fasta "$BASE_DIR/data/context/uniprot_human_reviewed_canonical.fasta" \
    --outdir "$BASE_DIR/results/disorder_context"
else
  echo "Skipping disorder-context staging: $BASE_DIR/data/context/mobidb_human_reviewed_disorder_intervals.tsv not found"
fi

if [[ "$INTERPRO_MODE" = "api" ]]; then
  run_step "stage InterPro domain intervals via API" python "$BASE_DIR/scripts/00_stage_context_sources.py" \
    --outdir "$BASE_DIR/data/context" \
    --accessions "$BASE_DIR/results/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv" \
    --interpro-source api \
    --workers "$INTERPRO_API_WORKERS" \
    --page-size "$INTERPRO_API_PAGE_SIZE" \
    --timeout "$INTERPRO_API_TIMEOUT" \
    --retries "$INTERPRO_API_RETRIES"
elif [[ "$INTERPRO_MODE" = "bulk" || ( "$INTERPRO_MODE" = "auto" && -f "$BASE_DIR/data/context/protein2ipr.dat.gz" ) ]]; then
  if [[ -f "$BASE_DIR/data/context/protein2ipr.dat.gz" ]]; then
    run_step "stage InterPro domain intervals" python "$BASE_DIR/scripts/00_stage_context_sources.py" \
      --outdir "$BASE_DIR/data/context" \
      --accessions "$BASE_DIR/results/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv" \
      --interpro-source bulk \
      --interpro-protein2ipr "$BASE_DIR/data/context/protein2ipr.dat.gz" \
      --bulk-parser pandas
  else
    echo "Skipping InterPro bulk staging: $BASE_DIR/data/context/protein2ipr.dat.gz not found"
  fi
fi

if [[ -f "$BASE_DIR/data/context/interpro_human_reviewed_domain_like_intervals.tsv" ]]; then
  run_step "annotate domain context" python "$BASE_DIR/scripts/05_annotate_domain_context.py" \
    --sites "$BASE_DIR/results/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
    --interpro-intervals "$BASE_DIR/data/context/interpro_human_reviewed_domain_like_intervals.tsv" \
    --position-col corrected_position \
    --outdir "$BASE_DIR/results/domain_context"
  run_step "summarize domain enrichment" python "$BASE_DIR/scripts/07_summarize_domain_enrichment.py" \
    --annotated-sites "$BASE_DIR/results/domain_context/sites_with_domain_context.tsv" \
    --interpro-intervals "$BASE_DIR/data/context/interpro_human_reviewed_domain_like_intervals.tsv" \
    --canonical-fasta "$BASE_DIR/data/context/uniprot_human_reviewed_canonical.fasta" \
    --outdir "$BASE_DIR/results/domain_context"
else
  echo "Skipping domain-context annotation: $BASE_DIR/data/context/interpro_human_reviewed_domain_like_intervals.tsv not found"
fi

if [[ -f "$BASE_DIR/data/context/mobidb_human_reviewed_disorder_intervals.tsv" || -f "$BASE_DIR/data/context/interpro_human_reviewed_domain_like_intervals.tsv" ]]; then
  COMPARE_ARGS=(
    --base-master "$BASE_DIR/data/base/human_ptm_master.tsv"
    --integrated-arg-sites "$BASE_DIR/results/integrated/human_arg_methyl_union_dedup_by_site.tsv"
    --canonical-fasta "$BASE_DIR/data/context/uniprot_human_reviewed_canonical.fasta"
    --outdir "$BASE_DIR/results/ptm_compare"
  )
  if [[ -f "$BASE_DIR/data/context/mobidb_human_reviewed_disorder_intervals.tsv" ]]; then
    COMPARE_ARGS+=(--disorder-intervals "$BASE_DIR/data/context/mobidb_human_reviewed_disorder_intervals.tsv")
  fi
  if [[ -f "$BASE_DIR/data/context/interpro_human_reviewed_domain_like_intervals.tsv" ]]; then
    COMPARE_ARGS+=(--interpro-intervals "$BASE_DIR/data/context/interpro_human_reviewed_domain_like_intervals.tsv")
  fi
  run_step "cross-PTM context comparison" python "$BASE_DIR/scripts/09_compare_ptm_contexts.py" "${COMPARE_ARGS[@]}"
fi

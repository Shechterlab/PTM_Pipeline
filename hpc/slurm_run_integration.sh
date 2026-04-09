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
RUN_CROSS_PTM_NEIGHBORS=${PTM_RUN_CROSS_PTM_NEIGHBORS:-0}
RUN_REVIEW_AUDIT=${PTM_RUN_REVIEW_AUDIT:-1}
RUN_CONDENSATE=${PTM_RUN_CONDENSATE:-auto}
CDCODE_DIR=${PTM_CDCODE_DIR:-}
CDCODE_PROTEINS=${PTM_CDCODE_PROTEINS:-}
CDCODE_PROTEIN_MAP=${PTM_CDCODE_PROTEIN_MAP:-}
CDCODE_CONDENSATES=${PTM_CDCODE_CONDENSATES:-}
CONDENSATE_STAGED=${PTM_CONDENSATE_STAGED:-}
CONDENSATE_OUTDIR=${PTM_CONDENSATE_OUTDIR:-}
CONDENSATE_EXPERIMENTAL_ONLY=${PTM_CONDENSATE_EXPERIMENTAL_ONLY:-0}
CONDENSATE_MIN_CONFIDENCE=${PTM_CONDENSATE_MIN_CONFIDENCE:-}
MARON_EXCLUDE_PATTERN=${PTM_MARON_EXCLUDE_PATTERN:-}
METHYL_CLUSTER_PERMUTATIONS=${PTM_METHYL_CLUSTER_PERMUTATIONS:-100000}
METHYL_CLUSTER_PROGRESS_EVERY=${PTM_METHYL_CLUSTER_PROGRESS_EVERY:-1000}

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
mkdir -p "$BASE_DIR/data/context" "$BASE_DIR/data/external/cd_code"

if [ -z "$CDCODE_DIR" ]; then
  CDCODE_DIR="$BASE_DIR/data/external/cd_code"
fi
if [ -z "$CDCODE_PROTEINS" ]; then
  CDCODE_PROTEINS="$CDCODE_DIR/proteins_202603181653.csv"
fi
if [ -z "$CDCODE_PROTEIN_MAP" ]; then
  CDCODE_PROTEIN_MAP="$CDCODE_DIR/protein2cdcode_v2.2.tsv"
fi
if [ -z "$CDCODE_CONDENSATES" ]; then
  CDCODE_CONDENSATES="$CDCODE_DIR/condensates_202603181647.csv"
fi
if [ -z "$CONDENSATE_STAGED" ]; then
  CONDENSATE_STAGED="$BASE_DIR/data/context/cd_code_condensate_proteins.tsv"
fi
if [ -z "$CONDENSATE_OUTDIR" ]; then
  CONDENSATE_OUTDIR="$BASE_DIR/results/condensate_enrichment"
fi

run_step "stage UniProt context" python "$BASE_DIR/scripts/00_stage_context_sources.py" \
  --outdir "$BASE_DIR/data/context" \
  --skip-interpro
run_step "parse Maron S5" python "$BASE_DIR/scripts/01_parse_maron_s5.py" \
  --input "$BASE_DIR/data/external/Table_S5_Compiled_Methylarginine_Data.xlsx" \
  --exclude-pattern "$MARON_EXCLUDE_PATTERN" \
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
if [[ ! -s "$BASE_DIR/data/context/interpro_human_reviewed_domain_like_intervals.tsv" ]] && [[ "$INTERPRO_MODE" = "api" ]]; then
  run_step "stage InterPro domain intervals via API" python "$BASE_DIR/scripts/00_stage_context_sources.py" \
    --outdir "$BASE_DIR/data/context" \
    --accessions "$BASE_DIR/results/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv" \
    --interpro-source api \
    --workers "$INTERPRO_API_WORKERS" \
    --page-size "$INTERPRO_API_PAGE_SIZE" \
    --timeout "$INTERPRO_API_TIMEOUT" \
    --retries "$INTERPRO_API_RETRIES"
elif [[ ! -s "$BASE_DIR/data/context/interpro_human_reviewed_domain_like_intervals.tsv" ]] && [[ "$INTERPRO_MODE" = "bulk" || ( "$INTERPRO_MODE" = "auto" && -f "$BASE_DIR/data/context/protein2ipr.dat.gz" ) ]]; then
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
FUNCTIONAL_ARGS=(
  --base-master "$BASE_DIR/data/base/human_ptm_master.tsv"
  --integrated-sites "$BASE_DIR/results/integrated/human_arg_methyl_union_dedup_by_site.tsv"
  --ontology "$BASE_DIR/config/functional_ontology.json"
  --outdir "$BASE_DIR/results/functional_class_union"
)
if [[ -s "$BASE_DIR/data/context/interpro_human_reviewed_domain_like_intervals.tsv" ]]; then
  FUNCTIONAL_ARGS+=(--interpro-intervals "$BASE_DIR/data/context/interpro_human_reviewed_domain_like_intervals.tsv")
fi
run_step "functional enrichment" python "$BASE_DIR/scripts/03_functional_class_enrichment.py" "${FUNCTIONAL_ARGS[@]}"
run_step "example figures" python "$BASE_DIR/scripts/04_example_figures.py" \
  --integrated-sites "$BASE_DIR/results/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --outdir "$BASE_DIR/results/example_figures"
MOTIF_ARGS=(
  --integrated-sites "$BASE_DIR/results/integrated/human_arg_methyl_union_dedup_by_site.tsv"
  --canonical-fasta "$BASE_DIR/data/context/uniprot_human_reviewed_canonical.fasta"
  --outdir "$BASE_DIR/results/motif_arg_odds"
)
if [[ -f "$BASE_DIR/data/context/mobidb_human_reviewed_disorder_intervals.tsv" ]]; then
  MOTIF_ARGS+=(--disorder-intervals "$BASE_DIR/data/context/mobidb_human_reviewed_disorder_intervals.tsv")
fi
run_step "motif enrichment and protein shrinkage model" python "$BASE_DIR/scripts/08_motif_and_arg_odds.py" "${MOTIF_ARGS[@]}"
run_step "methyl-Arg clustering" python "$BASE_DIR/scripts/10_methyl_arg_clustering.py" \
  --integrated-sites "$BASE_DIR/results/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --canonical-fasta "$BASE_DIR/data/context/uniprot_human_reviewed_canonical.fasta" \
  --permutations "$METHYL_CLUSTER_PERMUTATIONS" \
  --progress-every "$METHYL_CLUSTER_PROGRESS_EVERY" \
  --outdir "$BASE_DIR/results/methyl_arg_clustering"
run_step "cross-PTM clustering" python "$BASE_DIR/scripts/11_compare_ptm_clustering.py" \
  --base-master "$BASE_DIR/data/base/human_ptm_master.tsv" \
  --integrated-arg-sites "$BASE_DIR/results/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --canonical-fasta "$BASE_DIR/data/context/uniprot_human_reviewed_canonical.fasta" \
  --outdir "$BASE_DIR/results/ptm_clustering"

if [ "$RUN_CROSS_PTM_NEIGHBORS" = "1" ]; then
  CROSS_NEIGHBOR_ARGS=(
    --base-master "$BASE_DIR/data/base/human_ptm_master.tsv"
    --integrated-arg-sites "$BASE_DIR/results/integrated/human_arg_methyl_union_dedup_by_site.tsv"
    --canonical-fasta "$BASE_DIR/data/context/uniprot_human_reviewed_canonical.fasta"
    --outdir "$BASE_DIR/results/cross_ptm_neighbors"
  )
  if [[ -f "$BASE_DIR/data/context/mobidb_human_reviewed_disorder_intervals.tsv" ]]; then
    CROSS_NEIGHBOR_ARGS+=(--disorder-intervals "$BASE_DIR/data/context/mobidb_human_reviewed_disorder_intervals.tsv")
  fi
  run_step "cross-PTM neighbor enrichment" python "$BASE_DIR/scripts/12_cross_ptm_neighbor_enrichment.py" "${CROSS_NEIGHBOR_ARGS[@]}"
fi

if [[ -f "$BASE_DIR/data/context/mobidb_human_reviewed_disorder_intervals.tsv" ]]; then
  run_step "annotate disorder context" python "$BASE_DIR/scripts/06_annotate_disorder_context.py" \
    --sites "$BASE_DIR/results/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
    --disorder-intervals "$BASE_DIR/data/context/mobidb_human_reviewed_disorder_intervals.tsv" \
    --canonical-fasta "$BASE_DIR/data/context/uniprot_human_reviewed_canonical.fasta" \
    --outdir "$BASE_DIR/results/disorder_context"
else
  echo "Skipping disorder-context staging: $BASE_DIR/data/context/mobidb_human_reviewed_disorder_intervals.tsv not found"
fi

if [[ ! -s "$BASE_DIR/data/context/interpro_human_reviewed_domain_like_intervals.tsv" ]] && [[ "$INTERPRO_MODE" = "api" ]]; then
  run_step "stage InterPro domain intervals via API" python "$BASE_DIR/scripts/00_stage_context_sources.py" \
    --outdir "$BASE_DIR/data/context" \
    --accessions "$BASE_DIR/results/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv" \
    --interpro-source api \
    --workers "$INTERPRO_API_WORKERS" \
    --page-size "$INTERPRO_API_PAGE_SIZE" \
    --timeout "$INTERPRO_API_TIMEOUT" \
    --retries "$INTERPRO_API_RETRIES"
elif [[ ! -s "$BASE_DIR/data/context/interpro_human_reviewed_domain_like_intervals.tsv" ]] && [[ "$INTERPRO_MODE" = "bulk" || ( "$INTERPRO_MODE" = "auto" && -f "$BASE_DIR/data/context/protein2ipr.dat.gz" ) ]]; then
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

if [[ -s "$BASE_DIR/data/context/interpro_human_reviewed_domain_like_intervals.tsv" ]]; then
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

if [[ -f "$BASE_DIR/data/context/mobidb_human_reviewed_disorder_intervals.tsv" || -s "$BASE_DIR/data/context/interpro_human_reviewed_domain_like_intervals.tsv" ]]; then
  COMPARE_ARGS=(
    --base-master "$BASE_DIR/data/base/human_ptm_master.tsv"
    --integrated-arg-sites "$BASE_DIR/results/integrated/human_arg_methyl_union_dedup_by_site.tsv"
    --canonical-fasta "$BASE_DIR/data/context/uniprot_human_reviewed_canonical.fasta"
    --outdir "$BASE_DIR/results/ptm_compare"
  )
  if [[ -f "$BASE_DIR/data/context/mobidb_human_reviewed_disorder_intervals.tsv" ]]; then
    COMPARE_ARGS+=(--disorder-intervals "$BASE_DIR/data/context/mobidb_human_reviewed_disorder_intervals.tsv")
  fi
  if [[ -s "$BASE_DIR/data/context/interpro_human_reviewed_domain_like_intervals.tsv" ]]; then
    COMPARE_ARGS+=(--interpro-intervals "$BASE_DIR/data/context/interpro_human_reviewed_domain_like_intervals.tsv")
  fi
  run_step "cross-PTM context comparison" python "$BASE_DIR/scripts/09_compare_ptm_contexts.py" "${COMPARE_ARGS[@]}"
fi

if [ "$RUN_REVIEW_AUDIT" = "1" ] && [[ -f "$BASE_DIR/results/integrated/human_arg_methyl_union_all_rows.tsv" ]] && [[ -f "$BASE_DIR/results/integrated/human_arg_methyl_union_dedup_by_site.tsv" ]]; then
  run_step "review audit figures" python "$BASE_DIR/scripts/14_union_source_audit_and_review_figures.py" \
    --all-rows "$BASE_DIR/results/integrated/human_arg_methyl_union_all_rows.tsv" \
    --dedup-sites "$BASE_DIR/results/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
    --outdir "$BASE_DIR/results/review_audit"
fi

if [ "$RUN_CONDENSATE" = "1" ] || [ "$RUN_CONDENSATE" = "true" ] || [ "$RUN_CONDENSATE" = "yes" ] || [ "$RUN_CONDENSATE" = "auto" ]; then
  HAVE_STAGE=0
  HAVE_RAW=0
  if [[ -f "$CONDENSATE_STAGED" ]]; then
    HAVE_STAGE=1
  fi
  if [[ -f "$CDCODE_PROTEINS" ]] && [[ -f "$CDCODE_PROTEIN_MAP" ]] && [[ -f "$CDCODE_CONDENSATES" ]]; then
    HAVE_RAW=1
  fi
  if [ "$RUN_CONDENSATE" = "1" ] || [ "$RUN_CONDENSATE" = "true" ] || [ "$RUN_CONDENSATE" = "yes" ] || [ "$HAVE_STAGE" -eq 1 ] || [ "$HAVE_RAW" -eq 1 ]; then
    if [ "$HAVE_STAGE" -ne 1 ]; then
      run_step "stage CD-CODE condensates" python "$BASE_DIR/scripts/00_stage_cd_code_condensates.py" \
        --cdcode-proteins "$CDCODE_PROTEINS" \
        --cdcode-protein-map "$CDCODE_PROTEIN_MAP" \
        --cdcode-condensates "$CDCODE_CONDENSATES" \
        --out "$CONDENSATE_STAGED"
    fi
    COND_ARGS=(
      --base-master "$BASE_DIR/data/base/human_ptm_master.tsv"
      --integrated-arg-sites "$BASE_DIR/results/integrated/human_arg_methyl_union_dedup_by_site.tsv"
      --canonical-fasta "$BASE_DIR/data/context/uniprot_human_reviewed_canonical.fasta"
      --condensate-proteins "$CONDENSATE_STAGED"
      --outdir "$CONDENSATE_OUTDIR"
    )
    if [[ -f "$BASE_DIR/data/context/mobidb_human_reviewed_disorder_intervals.tsv" ]]; then
      COND_ARGS+=(--disorder-intervals "$BASE_DIR/data/context/mobidb_human_reviewed_disorder_intervals.tsv")
    fi
    if [ "$CONDENSATE_EXPERIMENTAL_ONLY" = "1" ]; then
      COND_ARGS+=(--experimental-only)
    fi
    if [ -n "$CONDENSATE_MIN_CONFIDENCE" ]; then
      COND_ARGS+=(--min-confidence-score "$CONDENSATE_MIN_CONFIDENCE")
    fi
    run_step "condensate PTM enrichment" python "$BASE_DIR/scripts/13_condensate_ptm_enrichment.py" "${COND_ARGS[@]}"
  else
    echo "Skipping condensate enrichment: no staged table or raw CD-CODE files found under $CDCODE_DIR"
  fi
fi

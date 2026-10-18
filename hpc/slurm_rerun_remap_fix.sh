#!/bin/bash
#SBATCH -J ptm_remapfix
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 24:00:00

# Rerun the pipeline after the scripts/02 remapping fixes into a fresh results root, leave the published results
# untouched, then compare the two roots table by table. The main pipeline stops if a source row is placed against its
# own evidence or a union site is not an arginine (see hpc/slurm_run_integration.sh); the site explorer data are only
# exported once every table they read has been rebuilt from the new union.
#
#   sbatch hpc/slurm_rerun_remap_fix.sh /path/to/PTM_Pipeline [NEW_RESULTS_ROOT]
#
# Old root: PTM_OLD_RESULTS_ROOT, else PTM_RESULTS_ROOT as hpc/ptm_env.sh resolves it (the published run).
# New root: the second argument, else PTM_NEW_RESULTS_ROOT, else <project>/PTM_results_remapfix.
# The staged inputs under PTM_DATA_ROOT are reused as they are (no fresh UniProt or InterPro download), so the
# remapping fix is the only change between the two runs. Clustering permutations and seeds keep their defaults.
#
# Options (environment):
#   PTM_RERUN_OVERWRITE=1      allow writing into a non-empty new root
#   PTM_RERUN_EXTENDED=0       skip scripts/run_extended_analyses.sh
#   PTM_RERUN_INTERACTIONS=1   also rebuild the explorer's interaction data (scripts/37; needs internet and its cache)
#   PTM_SKIP_HYGIENE=0         run the repo hygiene guard (it fails on a checkout that tracks PTM_results/)
#   PTM_ORG_SQLITE=path        org.Hs.eg.db SQLite for scripts/34 (default: the one installed with the R package)
#   PTM_MARON_TABLE_S2/_S3     Maron Table S2/S3 to rebuild maron_state_treatment/ (default: copy it from the old root)

set -eo pipefail

BASE_DIR=$(cd "${1:-$PWD}" && pwd)
NEW_ARG=${2:-}

CONDA_BASE=${PTM_CONDA_BASE:-/gs/gsfs0/hpc01/rhel8/apps/conda3}
if [ -x "$CONDA_BASE/bin/conda" ]; then
  eval "$("$CONDA_BASE/bin/conda" shell.bash hook)"
elif [ -f "$CONDA_BASE/bin/activate" ]; then
  . "$CONDA_BASE/bin/activate"
fi
if command -v conda >/dev/null 2>&1; then
  if [ -n "${PTM_ENV_PREFIX:-}" ]; then conda activate "$PTM_ENV_PREFIX"; else conda activate "${PTM_ENV_NAME:-ptm_pipeline}"; fi
fi

. "$BASE_DIR/hpc/ptm_env.sh"
OLD_ROOT=$(cd "${PTM_OLD_RESULTS_ROOT:-$PTM_RESULTS_ROOT}" && pwd)
NEW_ROOT=${NEW_ARG:-${PTM_NEW_RESULTS_ROOT:-$PTM_PROJECT_ROOT/PTM_results_remapfix}}
mkdir -p "$NEW_ROOT"
NEW_ROOT=$(cd "$NEW_ROOT" && pwd)

fail() { echo "[ptm_remapfix] $*"; exit 1; }

[ "$OLD_ROOT" != "$NEW_ROOT" ] || fail "the new results root must differ from the published one ($OLD_ROOT)"
[ -s "$OLD_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" ] || fail "no published union under $OLD_ROOT"
if [ -n "$(ls -A "$NEW_ROOT" 2>/dev/null)" ] && [ "${PTM_RERUN_OVERWRITE:-0}" != "1" ]; then
  fail "$NEW_ROOT is not empty; set PTM_RERUN_OVERWRITE=1 to reuse it"
fi
# Without the staged FASTA and metadata the first step would download a newer UniProt release and mix it into the comparison.
for f in "$PTM_CANONICAL_FASTA" "$PTM_UNIPROT_METADATA" "$PTM_MASTER_TSV" "$PTM_MARON_XLSX" "$PTM_PROMETHEUS_MMC2"; do
  [ -s "$f" ] || fail "staged input missing: $f"
done
grep -q "unresolved_evidence_mismatch" "$BASE_DIR/scripts/02_remap_arg_methyl_to_canonical.py" || fail "scripts/02 does not contain the remapping fixes"

{
  echo "started: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "code: $BASE_DIR"
  echo "commit: $(git -C "$BASE_DIR" rev-parse HEAD 2>/dev/null || echo unknown)"
  echo "uncommitted changes:"; git -C "$BASE_DIR" status --short 2>/dev/null | grep -v '^?? ' || true
  echo "old results root: $OLD_ROOT"
  echo "new results root: $NEW_ROOT"
  echo "data root: $PTM_DATA_ROOT"
  for f in "$PTM_CANONICAL_FASTA" "$PTM_UNIPROT_METADATA" "$PTM_MASTER_TSV" "$PTM_MARON_XLSX" "$PTM_PROMETHEUS_MMC2"; do
    echo "input: $(md5sum "$f")"
  done
  if [ -f "$PTM_CPLM4_METHYL" ]; then echo "input: $(md5sum "$PTM_CPLM4_METHYL")"; fi
  release_file="$(dirname "$PTM_CANONICAL_FASTA")/uniprot_release.json"
  echo "UniProt release record: $(cat "$release_file" 2>/dev/null | tr -d '\n' || echo 'none (FASTA staged before releases were recorded)')"
} > "$NEW_ROOT/RUN_INFO.txt"
cat "$NEW_ROOT/RUN_INFO.txt"

export PTM_RESULTS_ROOT="$NEW_ROOT"
export PTM_SKIP_HYGIENE=${PTM_SKIP_HYGIENE:-1}

echo "[ptm_remapfix] main pipeline (same steps and options as hpc/slurm_run_integration.sh)"
bash "$BASE_DIR/hpc/slurm_run_integration.sh" "$BASE_DIR"

if [ "${PTM_RERUN_EXTENDED:-1}" = "1" ]; then
  echo "[ptm_remapfix] extended analyses"
  PYTHON_BIN=${PYTHON_BIN:-python} bash "$BASE_DIR/scripts/run_extended_analyses.sh" "$BASE_DIR"
fi

# The disease-family analysis (scripts/34, then the figures of 35 and 33) counts methylated proteins, which the remap
# changes, and the site explorer reads its tables, so it is rebuilt here from the new root's labels and scripts/14 output.
DISEASE_DIR="$NEW_ROOT/functional_class_union/disease_enrichment"
[ -s "$DISEASE_DIR/disease_input_summary.tsv" ] || fail "no disease enrichment under $DISEASE_DIR (scripts/14 needs Rscript with clusterProfiler and DOSE)"
echo "[ptm_remapfix] disease-family analysis and figures"
Rscript "$BASE_DIR/scripts/34_export_dgn_mapping.R" "$DISEASE_DIR/dgn_term_gene.tsv"
ORG_SQLITE=${PTM_ORG_SQLITE:-$(Rscript -e 'cat(system.file("extdata", "org.Hs.eg.sqlite", package = "org.Hs.eg.db"))')}
[ -s "$ORG_SQLITE" ] || fail "org.Hs.eg.db SQLite not found; set PTM_ORG_SQLITE"
python "$BASE_DIR/scripts/34_regenerate_disease_ontology_figures.py" \
  --labels "$NEW_ROOT/functional_class_union/protein_functional_labels.tsv" \
  --term-gene "$DISEASE_DIR/dgn_term_gene.tsv" \
  --org-sqlite "$ORG_SQLITE" \
  --term-results "$DISEASE_DIR/all_methyl_proteins_disease_enrichment.tsv" \
  --input-summary "$DISEASE_DIR/disease_input_summary.tsv" \
  --outdir "$DISEASE_DIR"
python "$BASE_DIR/scripts/35_plot_disease_overview_clusterprofiler.py" --disease-dir "$DISEASE_DIR"
python "$BASE_DIR/scripts/33_plot_neurodegenerative_disease_focus.py" --disease-dir "$DISEASE_DIR"

# scripts/16 reads only Maron Table S2/S3 and the FASTA (same staged FASTA as the old root), not the union; the site
# explorer joins its tables to the new union by accession and position.
if [ -n "${PTM_MARON_TABLE_S2:-}" ] && [ -n "${PTM_MARON_TABLE_S3:-}" ]; then
  echo "[ptm_remapfix] Maron Table S2/S3 context"
  python "$BASE_DIR/scripts/16_maron_state_treatment_context.py" \
    --table-s2 "$PTM_MARON_TABLE_S2" --table-s3 "$PTM_MARON_TABLE_S3" \
    --canonical-fasta "$PTM_CANONICAL_FASTA" --outdir "$NEW_ROOT/maron_state_treatment"
else
  [ -d "$OLD_ROOT/maron_state_treatment" ] || fail "no maron_state_treatment/ in $OLD_ROOT; set PTM_MARON_TABLE_S2 and PTM_MARON_TABLE_S3"
  rm -rf "$NEW_ROOT/maron_state_treatment"
  cp -r "$OLD_ROOT/maron_state_treatment" "$NEW_ROOT/maron_state_treatment"
  echo "copied: maron_state_treatment/ from $OLD_ROOT (scripts/16 output; independent of the remap)" >> "$NEW_ROOT/RUN_INFO.txt"
fi

echo "[ptm_remapfix] site explorer data"
if [ -e "$PTM_CPLM4_METHYL" ]; then
  python "$BASE_DIR/scripts/39_stage_cplm4_lys_methylation.py" \
    --cplm "$PTM_CPLM4_METHYL" \
    --fasta "$PTM_CANONICAL_FASTA" \
    --outdir "$NEW_ROOT/cplm4"
else
  echo "[ptm_remapfix] no CPLM 4.0 download at $PTM_CPLM4_METHYL; Lys methylation from the master table only"
fi
python "$BASE_DIR/scripts/36_export_site_explorer_data.py" \
  --results-root "$NEW_ROOT" \
  --uniprot-metadata "$PTM_UNIPROT_METADATA" \
  --canonical-fasta "$PTM_CANONICAL_FASTA" \
  --dataset-version "${PTM_DATASET_VERSION:-remapfix-$(date -u +%Y%m%d)}" \
  --strict \
  --outdir "$NEW_ROOT/site_explorer_data"
if [ -s "$PTM_MOBIDB_INTERVALS" ]; then
  CPLM_ARGS=()
  [ -s "$NEW_ROOT/cplm4/cplm4_human_lys_methyl_validated.tsv" ] && CPLM_ARGS=(--cplm4-sites "$NEW_ROOT/cplm4/cplm4_human_lys_methyl_validated.tsv")
  python "$BASE_DIR/scripts/40_methylk_methylr_cooccurrence.py" \
    --sites-json "$NEW_ROOT/site_explorer_data/methylarg_sites.json" \
    --fasta "$PTM_CANONICAL_FASTA" \
    --disorder-intervals "$PTM_MOBIDB_INTERVALS" \
    "${CPLM_ARGS[@]}" \
    --outdir "$NEW_ROOT/methylk_methylr"
fi
if [ "${PTM_RERUN_INTERACTIONS:-0}" = "1" ]; then
  python "$BASE_DIR/scripts/37_export_interaction_data.py" \
    --sites-json "$NEW_ROOT/site_explorer_data/methylarg_sites.json" \
    --out "$NEW_ROOT/site_explorer_data/methylarg_interactions.json" \
    --cache "${PTM_INTERACTIONS_CACHE:-$BASE_DIR/data/interactions_cache}"
fi

echo "[ptm_remapfix] comparing $OLD_ROOT with $NEW_ROOT"
python "$BASE_DIR/scripts/38_compare_results_roots.py" \
  --old "$OLD_ROOT" \
  --new "$NEW_ROOT" \
  --outdir "$NEW_ROOT/rerun_comparison"

# Report only: the new root already passed the gate in the main pipeline; the old root shows what the fixes changed.
echo "[ptm_remapfix] checking every source row against its own evidence (old and new root)"
for root in "$OLD_ROOT" "$NEW_ROOT"; do
  python "$BASE_DIR/scripts/qc_remap_evidence_check.py" \
    --results-root "$root" \
    --fasta "$PTM_CANONICAL_FASTA" \
    --allow "${PTM_REMAP_ALLOWLIST:-$BASE_DIR/config/remap_evidence_allowlist.tsv}" \
    --outdir "$NEW_ROOT/rerun_comparison/evidence_check_$(basename "$root")" | tee "$NEW_ROOT/rerun_comparison/evidence_check_$(basename "$root").txt"
done

echo "finished: $(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$NEW_ROOT/RUN_INFO.txt"
echo "[ptm_remapfix] done; read $NEW_ROOT/rerun_comparison/report.md"

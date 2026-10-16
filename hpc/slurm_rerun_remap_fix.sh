#!/bin/bash
#SBATCH -J ptm_remapfix
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 24:00:00

# Rerun the pipeline after the scripts/02 remapping fix (remap_status sequence_position_verified) into a fresh
# results root, leave the published results untouched, then compare the two roots table by table.
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
grep -q "sequence_position_verified" "$BASE_DIR/scripts/02_remap_arg_methyl_to_canonical.py" || fail "scripts/02 does not contain the remapping fix"

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

# Scripts 16 (Maron Table S2/S3 only) and 33-35 (disease figures from the protein list, hard-coded to
# <repo>/PTM_results) do not read site positions, so they are not rerun.

echo "[ptm_remapfix] site explorer data"
python "$BASE_DIR/scripts/36_export_site_explorer_data.py" \
  --results-root "$NEW_ROOT" \
  --uniprot-metadata "$PTM_UNIPROT_METADATA" \
  --outdir "$NEW_ROOT/site_explorer_data"
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

echo "finished: $(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$NEW_ROOT/RUN_INFO.txt"
echo "[ptm_remapfix] done; read $NEW_ROOT/rerun_comparison/report.md"

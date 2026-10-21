#!/bin/bash
# Run the release analyses on a release's sites with the current code, and keep them in the release folder.
#
#   bash hpc/release_analyses.sh VERSION [WORK_ROOT]      e.g. bash hpc/release_analyses.sh v1.1.0
#
# Reads releases/VERSION/human_methylarginine_sites_VERSION.tsv (for v1.0.0, through
# releases/v1.1.0/supporting/site_changes_v1.0.0_to_v1.1.0.tsv, which places its sites on the current reference), runs
# hpc/analyses.sh on them, and writes releases/VERSION/analyses/ (scripts/release_analyses.py collect). The inputs must
# match inputs/MANIFEST.tsv; python and Rscript must be on PATH.
set -eo pipefail
VERSION=${1:?usage: hpc/release_analyses.sh VERSION [WORK_ROOT]}
CODE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
IN="$CODE/inputs"
FASTA="$IN/reference/uniprot_human_reviewed_2026_03.fasta"
ENTRIES="$IN/reference/uniprot_human_reviewed_2026_03_entries.tsv"
WORK=${2:-$(mktemp -d)}
SITES="${SITES:-$CODE/releases/$VERSION/human_methylarginine_sites_$VERSION.tsv}"
OUT="${ANALYSES_OUT:-$CODE/releases/$VERSION/analyses}"
MAP=''
[ "$VERSION" = v1.0.0 ] && MAP="$CODE/releases/v1.1.0/supporting/site_changes_v1.0.0_to_v1.1.0.tsv"
[ -f "$SITES" ] || { echo "[analyses] no site table $SITES"; exit 1; }
. "$CODE/hpc/analyses.sh"
echo "[analyses] $VERSION: checking inputs against $IN/MANIFEST.tsv"
check_inputs
python "$CODE/scripts/release_analyses.py" prepare --sites "$SITES" ${MAP:+--map "$MAP"} --entries "$ENTRIES" --work "$WORK"
echo "[analyses] $VERSION: analyses"
run_analyses
python "$CODE/scripts/release_analyses.py" collect --sites "$SITES" ${MAP:+--map "$MAP"} --work "$WORK" --inputs "$IN" \
  --version "$VERSION" --out "$OUT"
echo "[analyses] $VERSION: done -> $OUT"

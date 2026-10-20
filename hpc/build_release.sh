#!/bin/bash
# Build the corrected union from the frozen inputs and write the release tables.
#
#   bash hpc/build_release.sh RELEASE_NAME [WORK_ROOT]
#
# Runs, with the python on PATH: the two parsers, the remapping (with the accession history), the union, the evidence
# check and the residue audit (either stops the build), then scripts/build_release.py into release/RELEASE_NAME/.
# Every input comes from inputs/ and must match inputs/MANIFEST.tsv. WORK_ROOT (default: a temporary directory) holds
# the intermediate results root.

set -eo pipefail

NAME=${1:?usage: hpc/build_release.sh RELEASE_NAME [WORK_ROOT]}
CODE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
IN="$CODE/inputs"
WORK=${2:-$(mktemp -d)}
OUT="$CODE/release/$NAME"
mkdir -p "$WORK"

echo "[release] checking inputs against $IN/MANIFEST.tsv"
python - "$IN" <<'EOF'
import hashlib, sys
from pathlib import Path
import pandas as pd
inputs = Path(sys.argv[1])
bad = [r.file for r in pd.read_csv(inputs / 'MANIFEST.tsv', sep='\t').itertuples()
       if not (inputs / r.file).exists() or hashlib.md5((inputs / r.file).read_bytes()).hexdigest() != r.md5]
if bad:
    sys.exit(f'inputs differ from MANIFEST.tsv: {bad}')
EOF

{
  echo "built: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "commit: $(git -C "$CODE" rev-parse HEAD 2>/dev/null || echo unknown)"
  echo "uncommitted changes:"; git -C "$CODE" status --short 2>/dev/null | grep -v '^?? ' || true
  echo "python: $(python -c 'import sys, pandas, numpy; print(sys.version.split()[0], "pandas", pandas.__version__, "numpy", numpy.__version__)')"
  echo "inputs: $IN (all match MANIFEST.tsv)"
} > "$WORK/RUN_INFO.txt"

echo "[release] parsing source tables"
python "$CODE/scripts/01_parse_maron_s5.py" --input "$IN/maron2021_TableS5_Compiled_Methylarginine_Data.xlsx" --outdir "$WORK/parsed_maron"
python "$CODE/scripts/01_parse_prometheus_mmc2.py" --input "$IN/promteheusdb_supplement_mmc2.xlsx" --outdir "$WORK/parsed_prometheus"

echo "[release] placing records on the reference sequences"
python "$CODE/scripts/02_remap_arg_methyl_to_canonical.py" \
  --base-master "$IN/human_ptm_master.tsv" \
  --maron "$WORK/parsed_maron/maron_s5_arg_methyl_sites.tsv" \
  --prometheus "$WORK/parsed_prometheus/prometheus_mmc2_arg_methyl_sites.tsv" \
  --canonical-fasta "$IN/uniprot_human_reviewed_canonical.fasta" \
  --uniprot-metadata "$IN/uniprot_human_reviewed_canonical_metadata.tsv" \
  --accession-history "$IN/uniprot_accession_history.tsv" \
  --outdir "$WORK/remapped" > "$WORK/remap.log"
tail -1 "$WORK/remap.log"
python "$CODE/scripts/02_integrate_arg_methyl_sources.py" \
  --premerged-remapped "$WORK/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv" \
  --outdir "$WORK/integrated" > /dev/null

echo "[release] checking every record against its evidence, and every site's residue"
python "$CODE/scripts/qc_remap_evidence_check.py" --results-root "$WORK" \
  --fasta "$IN/uniprot_human_reviewed_canonical.fasta" \
  --fail-on-contradictions --allow "$CODE/config/remap_evidence_allowlist.tsv" > "$WORK/evidence_check.log"
grep "rows placed against" "$WORK/evidence_check.log"
python "$CODE/scripts/qc_remap_hand_check.py" --audit-union --results-root "$WORK" \
  --fasta "$IN/uniprot_human_reviewed_canonical.fasta" \
  --audit-out "$WORK/qc/remap_evidence_check/union_site_residue_audit.tsv"

echo "[release] writing $OUT"
python "$CODE/scripts/build_release.py" --release-root "$WORK" --published-root "$CODE/PTM_results" \
  --inputs "$IN" --allow "$CODE/config/remap_evidence_allowlist.tsv" --release-name "$NAME" --outdir "$OUT"
echo "[release] done: $OUT"

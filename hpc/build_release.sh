#!/bin/bash
# Build the corrected union from the frozen inputs and write a release.
#
#   bash hpc/build_release.sh VERSION [WORK_ROOT]        e.g. bash hpc/build_release.sh v1.1.0
#
# Runs, with the python on PATH:
#   1. the parsers, the remapping (with the accession history and Maron's peptide starts), the union, and the sites
#      whose window crosses a splice junction (kept as pending);
#   2. the evidence check and the residue audit: either stops the build;
#   3. scripts/build_release.py, which writes releases/VERSION/.
# Every input comes from inputs/ and must match inputs/MANIFEST.tsv (the files it lists as in Git). WORK_ROOT (default:
# a temporary directory) holds the intermediate results. The previous version is releases/v1.0.0/.

set -eo pipefail

VERSION=${1:?usage: hpc/build_release.sh VERSION [WORK_ROOT]}
CODE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
IN="$CODE/inputs"
WORK=${2:-$(mktemp -d)}
OUT="$CODE/releases/$VERSION"
FASTA="$IN/reference/uniprot_human_reviewed_2026_03.fasta"
ENTRIES="$IN/reference/uniprot_human_reviewed_2026_03_entries.tsv"
mkdir -p "$WORK"

echo "[release] checking inputs against $IN/MANIFEST.tsv"
python - "$IN" <<'EOF'
import hashlib, sys
from pathlib import Path
import pandas as pd
inputs = Path(sys.argv[1])
manifest = pd.read_csv(inputs / 'MANIFEST.tsv', sep='\t', dtype=str, keep_default_na=False)
bad = [r.file for r in manifest[manifest['redistribution'] == 'git'].itertuples()
       if not (inputs / r.file).exists() or hashlib.md5((inputs / r.file).read_bytes()).hexdigest() != r.md5]
if bad:
    sys.exit(f'inputs differ from MANIFEST.tsv: {bad}')
EOF

{
  echo "built: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "commit: $(git -C "$CODE" rev-parse HEAD 2>/dev/null || echo unknown)"
  echo "uncommitted changes:"; git -C "$CODE" status --short 2>/dev/null | grep -v '^?? ' || true
  echo "python: $(python -c 'import sys, pandas, numpy, openpyxl; print(sys.version.split()[0], "pandas", pandas.__version__, "numpy", numpy.__version__, "openpyxl", openpyxl.__version__)')"
  echo "inputs: $IN (all match MANIFEST.tsv)"
} > "$WORK/RUN_INFO.txt"

echo "[release] parsing source tables"
python "$CODE/scripts/01_parse_maron_s5.py" --input "$IN/sources/maron2021_tableS5_compiled_methylarginine_sites.xlsx" --outdir "$WORK/parsed_maron"
python "$CODE/scripts/01_parse_prometheus_mmc2.py" --input "$IN/sources/massignani2022_prometheusdb_mmc2_methylsites.xlsx" --outdir "$WORK/parsed_prometheus"
python "$CODE/scripts/01_parse_maron_peptide_coordinates.py" --table-s2 "$IN/sources/maron2021_tableS2_methylarginine_ip_sites.xlsx" \
  --table-s3 "$IN/sources/maron2021_tableS3_methylpeptides_prmt_inhibitors.xlsx" --outdir "$WORK/parsed_maron"
COORDS="$WORK/parsed_maron/maron_peptide_coordinates.tsv"

echo "[release] placing records on the reference sequences"
python "$CODE/scripts/02_remap_arg_methyl_to_canonical.py" \
  --base-master "$IN/sources/ptm_master_iptmnet_human.tsv" \
  --maron "$WORK/parsed_maron/maron_s5_arg_methyl_sites.tsv" \
  --prometheus "$WORK/parsed_prometheus/prometheus_mmc2_arg_methyl_sites.tsv" \
  --canonical-fasta "$FASTA" \
  --uniprot-metadata "$ENTRIES" \
  --accession-history "$IN/reference/uniprot_accession_history_2026_03.tsv" \
  --peptide-coordinates "$COORDS" \
  --outdir "$WORK/remapped" > "$WORK/remap.log"
tail -1 "$WORK/remap.log"
python "$CODE/scripts/02_integrate_arg_methyl_sources.py" \
  --premerged-remapped "$WORK/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv" \
  --outdir "$WORK/integrated" > /dev/null
python "$CODE/scripts/02_find_junction_window_sites.py" --results-root "$WORK" --fasta "$FASTA"

echo "[release] checking every record against its evidence, and every site's residue"
python "$CODE/scripts/qc_remap_evidence_check.py" --results-root "$WORK" --fasta "$FASTA" --peptide-coordinates "$COORDS" \
  --fail-on-contradictions --allow "$CODE/config/remap_evidence_allowlist.tsv" > "$WORK/evidence_check.log"
grep "rows placed against" "$WORK/evidence_check.log"
python "$CODE/scripts/qc_remap_hand_check.py" --audit-union --results-root "$WORK" --fasta "$FASTA" \
  --audit-out "$WORK/qc/remap_evidence_check/union_site_residue_audit.tsv"

echo "[release] writing $OUT"
python "$CODE/scripts/build_release.py" --results-root "$WORK" --previous "$CODE/releases/v1.0.0" --inputs "$IN" \
  --allow "$CODE/config/remap_evidence_allowlist.tsv" --version "$VERSION" --outdir "$OUT"
echo "[release] done: $OUT"

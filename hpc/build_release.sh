#!/bin/bash
# Build the corrected union from the frozen inputs and write a release.
#
#   bash hpc/build_release.sh VERSION [WORK_ROOT]        e.g. bash hpc/build_release.sh v1.2.0
#
# Runs, with the python on PATH:
#   1. the parsers, the placement of every record on the reference sequences (with the accession history, Maron's
#      peptide starts and the study tables of inputs/sources/; a record adds a site only where its source localises the
#      methyl group to one arginine, scripts/02 --localised-only), the union, and the sites whose windows match the entry
#      on one side of the arginine only (kept pending a decision);
#   2. the evidence check and the residue audit: either stops the build;
#   3. the analyses the release and the website report, on the union (the sites kept pending are left out): functional
#      classes, motifs and arginine enrichment, disorder and domain context, other PTMs, clustering and density,
#      domain-cluster architectures, condensates, the Maron 2021 summaries and disease enrichment; the tests of protein
#      and gene sets run clusterProfiler and DOSE in R (scripts/enrich_sets.R), so Rscript must be on PATH;
#   4. scripts/build_release.py, a first pass: the sites, and the lower-confidence sites (arginines a study reports, or an
#      earlier release had, whose records each fail a site test);
#   5. the analyses again, on the sites and the lower-confidence sites whose own evidence covers the arginine and leaves
#      the methyl group on an arginine (WORK_ROOT/with_lower_confidence), and the view with both tiers
#      (build_release.py --with-lower-confidence), which the website reads when the reader switches them on;
#   6. the website's data from the release and from that view (scripts/36) and the interaction evidence of their sites
#      (scripts/37, from inputs/interactions/), then scripts/build_release.py again with the interaction tables, which
#      writes releases/VERSION/ with the statistics of the second run and both analyses' result tables (analyses/,
#      analyses_with_lower_confidence/); the website data go to docs/data/ (SITE_OUT), the view's as methylarg_*_all.json;
#   7. the annotation layers (scripts/release_layers.py): checked against MANIFEST.tsv first; each build_release.py pass
#      adds their site and protein columns and statistics, and their website files are written from the release last.
# Every input comes from inputs/ and must match inputs/MANIFEST.tsv (the files it lists as in Git, and the pCLAP and
# RBDmap peptides, which scripts/43d and scripts/43f derive from PRIDE when they are absent). WORK_ROOT (default: a temporary directory) holds the
# intermediate results. The build needs the network for those derivations and when inputs/interactions/ lacks data for a site: scripts/37 then fetches it into WORK_ROOT/interactions_topup, and the build stops, before
# writing the release, and says how to stage it.
#
# Environment (optional): PREVIOUS, the release the changes are traced from (default releases/v1.2.0); RELEASE_OUT, where
# the release goes (default releases/VERSION); SITE_OUT, where the website data go (default docs/data); STUDIES, the
# study tables read (default all sixteen, inputs/sources/<study>_methylarginine_records.tsv); LOCALISED_ONLY (default 1),
# 0 to place records without localisation as v1.1.0 did. v1.1.0 itself: PREVIOUS=releases/v1.0.0 STUDIES=''
# LOCALISED_ONLY=0.

set -eo pipefail

VERSION=${1:?usage: hpc/build_release.sh VERSION [WORK_ROOT]}
CODE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
IN="$CODE/inputs"
WORK=${2:-$(mktemp -d)}
OUT=${RELEASE_OUT:-$CODE/releases/$VERSION}
PREVIOUS=${PREVIOUS:-$CODE/releases/v1.2.0}
STUDIES=${STUDIES-larsen2016 maron2021 guo2014 musiani2019 fong2019 li2021 fedoriw2019 uniprot lim2020 wei2020 iedb hartel2019 geoghegan2015 sylvestersen2014 lu2025 ma2023}
LOCALISED_ONLY=${LOCALISED_ONLY:-1}
SITE_OUT=${SITE_OUT:-$CODE/docs/data}
FASTA="$IN/reference/uniprot_human_reviewed_2026_03.fasta"
ENTRIES="$IN/reference/uniprot_human_reviewed_2026_03_entries.tsv"
mkdir -p "$WORK"
if [ "$(realpath -m "$PREVIOUS")" = "$(realpath -m "$OUT")" ]; then
  echo "[release] PREVIOUS and the output are the same folder ($OUT): set PREVIOUS to the release before $VERSION"; exit 1
fi

echo "[release] checking inputs against $IN/MANIFEST.tsv"
. "$CODE/hpc/analyses.sh"
check_inputs
# the annotation layers the release joins to its sites and proteins (scripts/release_layers.py LAYERS)
python "$CODE/scripts/release_layers.py" check --inputs "$IN"

# Tables derived from PRIDE and kept out of Git (EMBL-EBI terms of use): the pCLAP peptides (scripts/43d, from the
# MaxQuant output) and the RBDmap measured peptides (scripts/43f, from the peptide and quantification tables). When one
# is absent, its script derives it (download into PTM_data/external/) in WORK, leaving the files in Git alone. Either way
# each must match MANIFEST.tsv, or scripts/37 would test it silently changed or leave it out.
derive_input() {  # derive_input TABLE SCRIPT [ARGS...]
  local table=$1 script=$2 name
  shift 2
  name=$(basename "$table" .tsv)
  if [ ! -f "$IN/$table" ]; then
    echo "[release] deriving $IN/$table ($script)"
    mkdir -p "$WORK/logs" "$WORK/derived/$name/interactions"
    python "$CODE/$script" "$@" --out "$WORK/derived/$name/interactions" --manifest '' > "$WORK/logs/stage_$name.log" 2>&1 \
      || { tail -20 "$WORK/logs/stage_$name.log"; echo "[release] could not derive $table (log: $WORK/logs/stage_$name.log)"; exit 1; }
    cp "$WORK/derived/$name/$table" "$IN/$table"
  fi
  python - "$IN" "$table" "$script" <<'PYEOF'
import hashlib, sys
from pathlib import Path
import pandas as pd
inputs, name, script = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
manifest = pd.read_csv(inputs / 'MANIFEST.tsv', sep='\t', dtype=str, keep_default_na=False).set_index('file')
if hashlib.md5((inputs / name).read_bytes()).hexdigest() != manifest.at[name, 'md5']:
    sys.exit(f'{name} differs from MANIFEST.tsv: move it aside and build again to derive it ({script})')
PYEOF
}
derive_input interactions/pclap_mullari2017_peptides.tsv scripts/43d_stage_pclap.py
derive_input interactions/rbdmap_castello2016_measured_peptides.tsv scripts/43f_stage_peptide_backgrounds.py --only rbdmap

{
  echo "built: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "commit: $(git -C "$CODE" rev-parse HEAD 2>/dev/null || echo unknown)"
  echo "uncommitted changes:"; git -C "$CODE" status --short 2>/dev/null | grep -v '^?? ' || true
  echo "python: $(python -c 'import sys, pandas, numpy, openpyxl; print(sys.version.split()[0], "pandas", pandas.__version__, "numpy", numpy.__version__, "openpyxl", openpyxl.__version__)')"
  echo "inputs: $IN (all match MANIFEST.tsv)"
} > "$WORK/RUN_INFO.txt"

echo "[release] parsing source tables"
python "$CODE/scripts/01_parse_maron_s5.py" --input "$IN/sources/maron2021_tableS5_compiled_methylarginine_sites.xlsx" \
  --canonical-fasta "$FASTA" --outdir "$WORK/parsed_maron"
python "$CODE/scripts/01_parse_prometheus_mmc2.py" --input "$IN/sources/massignani2022_prometheusdb_mmc2_methylsites.xlsx" --outdir "$WORK/parsed_prometheus"
python "$CODE/scripts/01_parse_maron_peptide_coordinates.py" --table-s2 "$IN/sources/maron2021_tableS2_methylarginine_ip_sites.xlsx" \
  --table-s3 "$IN/sources/maron2021_tableS3_methylpeptides_prmt_inhibitors.xlsx" --outdir "$WORK/parsed_maron"
COORDS="$WORK/parsed_maron/maron_peptide_coordinates.tsv"
# The study tables staged by scripts/00_stage_study_sites.py; their records follow those of the three compiled tables,
# in this order, so that record_ids stay those of earlier releases.
STUDY_RECORDS=""
for study in $STUDIES; do
  STUDY_RECORDS="$STUDY_RECORDS $IN/sources/${study}_methylarginine_records.tsv"
done
LOCALISED=""
[ "$LOCALISED_ONLY" = 1 ] && LOCALISED="--localised-only"

echo "[release] placing records on the reference sequences"
python "$CODE/scripts/02_remap_arg_methyl_to_canonical.py" \
  --base-master "$IN/sources/ptm_master_iptmnet_human.tsv" \
  --maron "$WORK/parsed_maron/maron_s5_arg_methyl_sites.tsv" \
  --prometheus "$WORK/parsed_prometheus/prometheus_mmc2_arg_methyl_sites.tsv" \
  --maron-extra "$WORK/parsed_maron/maron_s5_rows_accession_na.tsv" \
  --canonical-fasta "$FASTA" \
  --uniprot-metadata "$ENTRIES" \
  --accession-history "$IN/reference/uniprot_accession_history_2026_03.tsv" \
  --peptide-coordinates "$COORDS" \
  --study-records $STUDY_RECORDS $LOCALISED \
  --outdir "$WORK/remapped" > "$WORK/remap.log"
tail -1 "$WORK/remap.log"
python "$CODE/scripts/02_integrate_arg_methyl_sources.py" \
  --premerged-remapped "$WORK/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv" \
  --outdir "$WORK/integrated" > /dev/null
python "$CODE/scripts/02_find_one_sided_window_sites.py" --results-root "$WORK" --fasta "$FASTA"

echo "[release] checking every record against its evidence, and every site's residue"
python "$CODE/scripts/qc_remap_evidence_check.py" --results-root "$WORK" --fasta "$FASTA" --peptide-coordinates "$COORDS" \
  --fail-on-contradictions --allow "$CODE/config/remap_evidence_allowlist.tsv" > "$WORK/evidence_check.log" \
  || { tail -45 "$WORK/evidence_check.log"; echo "[release] records placed against their evidence (log: $WORK/evidence_check.log)"; exit 1; }
grep "rows placed against" "$WORK/evidence_check.log"
python "$CODE/scripts/qc_remap_hand_check.py" --audit-union --results-root "$WORK" --fasta "$FASTA" \
  --audit-out "$WORK/qc/remap_evidence_check/union_site_residue_audit.tsv"

echo "[release] analyses on the union"
run_analyses

release() {
  python "$CODE/scripts/build_release.py" --results-root "$WORK" --previous "$PREVIOUS" --inputs "$IN" \
    --allow "$CODE/config/remap_evidence_allowlist.tsv" --version "$VERSION" "$@"
}
# The first pass goes to WORK, not to OUT: the website data and the interaction evidence are made from it, and OUT is
# written once, with the interaction evidence.
PASS1="$WORK/release_pass1"
echo "[release] first pass: $PASS1"
release --outdir "$PASS1" > "$WORK/logs/release_first_pass.log"
SITES_TABLE="$PASS1/human_methylarginine_sites_$VERSION.tsv"
LOWER_TABLE="$PASS1/human_methylarginine_lower_confidence_sites_$VERSION.tsv"

# The view with lower-confidence sites: the analyses again, on the sites and the lower-confidence sites whose own
# evidence covers the arginine and leaves the methyl group on an arginine; its site table holds both tiers. The website
# reads it when the reader switches lower-confidence sites on, and the release keeps its statistics.
WORK_LC="$WORK/with_lower_confidence"
echo "[release] analyses with lower-confidence sites: $WORK_LC"
python "$CODE/scripts/release_analyses.py" prepare --sites "$SITES_TABLE" --lower-confidence "$LOWER_TABLE" --entries "$ENTRIES" \
  --work "$WORK_LC"
( WORK="$WORK_LC"; run_analyses )
VIEW1="$WORK/view_pass1"
release --outdir "$VIEW1" --with-lower-confidence --analyses-root "$WORK_LC" > "$WORK/logs/view_first_pass.log"

echo "[release] website data and interaction evidence"
SITE="$WORK/site_data"
SITE_ALL="$WORK/site_data_with_lower_confidence"
step "site data (scripts/36)" python "$CODE/scripts/36_export_site_explorer_data.py" --release "$PASS1" --inputs "$IN" \
  --outdir "$SITE" --check
step "site data with lower-confidence sites (scripts/36)" python "$CODE/scripts/36_export_site_explorer_data.py" \
  --release "$VIEW1" --inputs "$IN" --outdir "$SITE_ALL" --check
rm -rf "$WORK/interactions_topup"
interactions() {   # $1 site data folder, $2 interaction tables out
  python "$CODE/scripts/37_export_interaction_data.py" \
    --sites-json "$1/methylarg_sites.json" --out "$1/methylarg_interactions.json" \
    --pred-iface-out "$1/methylarg_pred_iface.json" --tables "$IN/interactions" --topup-cache "$WORK/interactions_topup" \
    --tables-out "$2" --reference-fasta "$FASTA" --nucleic-acid "$IN/nucleic_acid" \
    --measured-sites "$CODE/releases/v1.0.0/human_methylarginine_sites_v1.0.0.tsv" --unlisted-out "$1/methylarg_binding_unlisted.json"
}
# both tiers first: when they lack data, so may the sites, and one top-up covers both
step "interaction evidence with lower-confidence sites (scripts/37)" interactions "$SITE_ALL" "$WORK/interactions_with_lower_confidence"
if [ -z "$(find "$WORK/interactions_topup" -type f 2>/dev/null | head -1)" ]; then
  step "interaction evidence (scripts/37)" interactions "$SITE" "$WORK/interactions"
fi
if [ -n "$(find "$WORK/interactions_topup" -type f 2>/dev/null | head -1)" ]; then
  echo "[release] inputs/interactions/ lacks data for some sites; scripts/37 fetched it into $WORK/interactions_topup."
  echo "          Stage it with: python scripts/37a_stage_interactions.py --topup $WORK/interactions_topup, then build again."
  exit 1
fi

echo "[release] writing $OUT with the interaction evidence"
VIEW="$WORK/view"
release --outdir "$VIEW" --with-lower-confidence --analyses-root "$WORK_LC" --interactions "$WORK/interactions_with_lower_confidence" \
  > "$WORK/logs/view.log"
release --outdir "$OUT" --interactions "$WORK/interactions" --statistics-with-lower-confidence "$VIEW/statistics_$VERSION.tsv"
# The result tables of both analyses, in the release folder (scripts/release_analyses.py collect)
python "$CODE/scripts/release_analyses.py" collect --sites "$OUT/human_methylarginine_sites_$VERSION.tsv" --work "$WORK" \
  --inputs "$IN" --version "$VERSION" --out "$OUT/analyses" > /dev/null
python "$CODE/scripts/release_analyses.py" collect --sites "$OUT/human_methylarginine_sites_$VERSION.tsv" \
  --lower-confidence "$OUT/human_methylarginine_lower_confidence_sites_$VERSION.tsv" --work "$WORK_LC" \
  --inputs "$IN" --version "$VERSION" --out "$OUT/analyses_with_lower_confidence" > /dev/null
# The release's checksums last: collect has just rewritten analyses*/ (their README names the commit), so the top-level
# MD5SUMS written with the release is out of date until it is written again here. Every manifest is then verified.
python -c "import sys; sys.path.insert(0, '$CODE/scripts'); from pathlib import Path; import release_docs; release_docs.write_md5sums(Path('$OUT'))"
for sums in "$OUT/MD5SUMS" "$OUT/analyses/MD5SUMS" "$OUT/analyses_with_lower_confidence/MD5SUMS"; do
  (cd "$(dirname "$sums")" && md5sum -c --quiet MD5SUMS) || { echo "[release] checksums do not match in $(dirname "$sums")"; exit 1; }
done
echo "[release] checksums verified: $OUT"
for f in sites evidence disease interactions; do
  cp "$SITE/methylarg_$f.json" "$SITE_OUT/methylarg_$f.json"
  cp "$SITE_ALL/methylarg_$f.json" "$SITE_OUT/methylarg_${f}_all.json"
done
# InterPro domains and families of the explorer's proteins (scripts/36c), one file for both tiers; and every other reviewed
# protein (an index and 64 record files under proteome/) for the light pages of proteins without a site
step "InterPro domains and the proteins without a site (scripts/36c)" python "$CODE/scripts/36c_export_interpro_domains.py" \
  --intervals "$IN/annotation/interpro_108.0_domain_intervals.tsv" --all-intervals "$IN/annotation/interpro_108.0_all_entry_intervals.tsv" \
  --entries "$IN/reference/uniprot_human_reviewed_2026_03_entries.tsv" --disorder "$IN/annotation/mobidb_disorder_intervals_2026-03-30.tsv" \
  --cdcode "$IN/annotation/cdcode_protein_membership_v2.2.tsv" --sites-json "$SITE_OUT/methylarg_sites.json" --sites-all-json "$SITE_OUT/methylarg_sites_all.json" \
  --out-dir "$SITE_OUT" --pred-iface "$SITE_ALL/methylarg_pred_iface.json" --proteome-out "$SITE_OUT/methylarg_proteome.json" \
  --binding "$SITE/methylarg_binding_unlisted.json"
# eCLIP logos (scripts/43f), one file for both tiers: the 5-mers enriched in each experiment's peaks or windows, aligned
step "eCLIP logos (scripts/43f)" python "$CODE/scripts/43f_eclip_logos.py" --inputs "$IN/nucleic_acid" --out "$SITE_OUT/methylarg_eclip_logos.json"
# Each annotation layer's website files, from the release (scripts/release_layers.py: DepMap, PRISM and TCGA to
# SITE_OUT/depmap/ by scripts/44c, with the release's statistics of the layer)
step "annotation layers' website files (scripts/release_layers.py)" python "$CODE/scripts/release_layers.py" website --release "$OUT" \
  --inputs "$IN" --site-out "$SITE_OUT"
echo "[release] done: $OUT and $SITE_OUT/"

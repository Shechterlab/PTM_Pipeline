#!/bin/bash
# Build the corrected union from the frozen inputs and write a release.
#
#   bash hpc/build_release.sh VERSION [WORK_ROOT]        e.g. bash hpc/build_release.sh v1.1.0
#
# Runs, with the python on PATH:
#   1. the parsers, the remapping (with the accession history, Maron's peptide starts and the studies staged in
#      inputs/sources/), the union, and the sites
#      whose windows match the entry on one side of the arginine only (kept pending a decision);
#   2. the evidence check and the residue audit: either stops the build;
#   3. the analyses the release and the website report, on the union (the sites kept pending are left out): functional
#      classes, motifs and arginine enrichment, disorder and domain context, other PTMs, clustering and density,
#      domain-cluster architectures, condensates, the Maron 2021 summaries and disease enrichment; the tests of protein
#      and gene sets run clusterProfiler and DOSE in R (scripts/enrich_sets.R), so Rscript must be on PATH;
#   4. scripts/build_release.py, which writes releases/VERSION/;
#   5. the website's data from the release (scripts/36) and the interaction evidence of its sites (scripts/37, from
#      inputs/interactions/), then scripts/build_release.py again with the interaction tables; the website data go to
#      docs/data/ (SITE_OUT).
# Every input comes from inputs/ and must match inputs/MANIFEST.tsv (the files it lists as in Git). WORK_ROOT (default:
# a temporary directory) holds the intermediate results. The build needs the network only when inputs/interactions/
# lacks data for a site: scripts/37 then fetches it into WORK_ROOT/interactions_topup, and the build stops, before
# writing the release, and says how to stage it.
#
# Environment (optional): PREVIOUS, the release the changes are traced from (default releases/v1.1.0); RELEASE_OUT, where
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
PREVIOUS=${PREVIOUS:-$CODE/releases/v1.1.0}
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
U="$WORK/integrated/human_arg_methyl_union_dedup_by_site.tsv"
DOMAINS="$IN/annotation/interpro_108.0_domain_intervals.tsv"
DISORDER="$IN/annotation/mobidb_disorder_intervals_2026-03-30.tsv"
MASTER="$IN/sources/ptm_master_iptmnet_human.tsv"
export PTM_SKIP_PDF=1 MPLCONFIGDIR="$WORK/.matplotlib"
step() {
  local name="$1" log
  shift
  log="$WORK/logs/$(printf '%s' "$name" | tr -cs 'A-Za-z0-9' _).log"
  echo "  $name"
  "$@" > "$log" 2>&1 || { tail -20 "$log"; echo "[release] failed: $name (log: $log)"; exit 1; }
}
mkdir -p "$WORK/logs"
step "functional classes (scripts/03)" python "$CODE/scripts/03_functional_class_enrichment.py" --base-master "$MASTER" \
  --integrated-sites "$U" --uniprot-metadata "$ENTRIES" --interpro-intervals "$DOMAINS" \
  --ontology "$CODE/config/functional_ontology.json" --domain-ontology "$CODE/config/domain_class_ontology.json" \
  --outdir "$WORK/functional_class_union"
step "PANTHER protein classes (scripts/03b)" python "$CODE/scripts/03b_panther_class_enrichment.py" --integrated-sites "$U" \
  --uniprot-metadata "$ENTRIES" --protein-classes "$IN/annotation/panther_19.0_protein_classes.tsv" \
  --ontology "$IN/annotation/panther_19.0_protein_class_ontology.tsv" --groups "$CODE/config/panther_class_groups.json" \
  --outdir "$WORK/panther_class"
step "motifs and arginine enrichment (scripts/08)" python "$CODE/scripts/08_motif_and_arg_odds.py" --integrated-sites "$U" \
  --canonical-fasta "$FASTA" --outdir "$WORK/motif_arg_odds"
step "disorder context (scripts/06)" python "$CODE/scripts/06_annotate_disorder_context.py" --sites "$U" \
  --disorder-intervals "$DISORDER" --canonical-fasta "$FASTA" --outdir "$WORK/disorder_context"
step "domain context (scripts/05)" python "$CODE/scripts/05_annotate_domain_context.py" --sites "$U" \
  --interpro-intervals "$DOMAINS" --position-col corrected_position --outdir "$WORK/domain_context"
step "domain enrichment (scripts/07)" python "$CODE/scripts/07_summarize_domain_enrichment.py" \
  --annotated-sites "$WORK/domain_context/sites_with_domain_context.tsv" --interpro-intervals "$DOMAINS" \
  --canonical-fasta "$FASTA" --ontology "$CODE/config/domain_class_ontology.json" \
  --disorder-sites "$WORK/disorder_context/sites_with_disorder_context.tsv" --disorder-intervals "$DISORDER" \
  --outdir "$WORK/domain_context"
step "other PTMs (scripts/09)" python "$CODE/scripts/09_compare_ptm_contexts.py" --base-master "$MASTER" \
  --integrated-arg-sites "$U" --canonical-fasta "$FASTA" --disorder-intervals "$DISORDER" --interpro-intervals "$DOMAINS" \
  --outdir "$WORK/ptm_compare"
# Clustering: 100 permutations, seed 42, as in the published run (its clustering_qc.tsv records them). Multiscale
# density: 10,000 permutations per protein (from v1.2.0, 1,000 before), each protein seeded from 42 and its accession.
step "clustering (scripts/10)" python "$CODE/scripts/10_methyl_arg_clustering.py" --integrated-sites "$U" \
  --canonical-fasta "$FASTA" --permutations 100 --disorder-intervals "$DISORDER" --null-match-disorder \
  --outdir "$WORK/methyl_arg_clustering"
step "multiscale density (scripts/28)" python "$CODE/scripts/28_multiscale_methylarg_density_null.py" \
  --sites "$WORK/domain_context/sites_with_domain_context.tsv" --fasta "$FASTA" --permutations 10000 \
  --outdir "$WORK/methyl_arg_clustering"
step "domain-cluster architectures (scripts/25)" python "$CODE/scripts/25_classify_domain_cluster_architectures.py" \
  --clusters "$WORK/domain_context/domain_edge_class_20aa_cluster_examples.tsv" \
  --site-context "$WORK/domain_context/sites_with_domain_context.tsv" --outdir "$WORK/domain_context" --cluster-threshold 3
step "condensates (scripts/12)" python "$CODE/scripts/12_condensate_enrichment.py" --integrated-sites "$U" \
  --staged-membership "$IN/annotation/cdcode_protein_membership_v2.2.tsv" \
  --proteins-table "$IN/annotation/cdcode_proteins_2026-03-18.csv" \
  --condensates-table "$IN/annotation/cdcode_condensates_2026-03-18.csv" --reference-fasta "$FASTA" \
  --accession-history "$IN/reference/uniprot_accession_history_2026_03.tsv" --outdir "$WORK/condensates"
step "Maron 2021 summaries (scripts/16)" python "$CODE/scripts/16_maron_state_treatment_context.py" \
  --table-s2 "$IN/sources/maron2021_tableS2_methylarginine_ip_sites.xlsx" \
  --table-s3 "$IN/sources/maron2021_tableS3_methylpeptides_prmt_inhibitors.xlsx" --canonical-fasta "$FASTA" \
  --outdir "$WORK/maron_state_treatment"
step "disease enrichment (scripts/14)" python "$CODE/scripts/14_disease_enrichment.py" \
  --labels "$WORK/functional_class_union/protein_functional_labels.tsv" \
  --term-genes "$IN/disease/dose_4.4.0_disgenet_term_genes.tsv" \
  --symbol-entrez "$IN/disease/orghsegdb_3.22.0_symbol_entrez.tsv" \
  --ontology-genes "$IN/disease/hdo_2026-09-30_term_genes.tsv" --ontology-terms "$IN/disease/hdo_2026-09-30_terms.tsv" \
  --uniprot-metadata "$ENTRIES" \
  --outdir "$WORK/disease"

release() {
  python "$CODE/scripts/build_release.py" --results-root "$WORK" --previous "$PREVIOUS" --inputs "$IN" \
    --allow "$CODE/config/remap_evidence_allowlist.tsv" --version "$VERSION" "$@"
}
# The first pass goes to WORK, not to OUT: the website data and the interaction evidence are made from it, and OUT is
# written once, with the interaction evidence.
PASS1="$WORK/release_pass1"
echo "[release] first pass: $PASS1"
release --outdir "$PASS1" > "$WORK/logs/release_first_pass.log"

echo "[release] website data and interaction evidence"
SITE="$WORK/site_data"
step "site data (scripts/36)" python "$CODE/scripts/36_export_site_explorer_data.py" --release "$PASS1" --inputs "$IN" \
  --outdir "$SITE" --check
rm -rf "$WORK/interactions_topup"
step "interaction evidence (scripts/37)" python "$CODE/scripts/37_export_interaction_data.py" \
  --sites-json "$SITE/methylarg_sites.json" --out "$SITE/methylarg_interactions.json" \
  --pred-iface-out "$SITE/methylarg_pred_iface.json" --tables "$IN/interactions" --topup-cache "$WORK/interactions_topup" \
  --tables-out "$WORK/interactions" --reference-fasta "$FASTA" \
  --measured-sites "$CODE/releases/v1.0.0/human_methylarginine_sites_v1.0.0.tsv"
if [ -n "$(find "$WORK/interactions_topup" -type f 2>/dev/null | head -1)" ]; then
  echo "[release] inputs/interactions/ lacks data for some sites; scripts/37 fetched it into $WORK/interactions_topup."
  echo "          Stage it with: python scripts/37a_stage_interactions.py --topup $WORK/interactions_topup, then build again."
  exit 1
fi

echo "[release] writing $OUT with the interaction evidence"
release --outdir "$OUT" --interactions "$WORK/interactions"
cp "$SITE"/methylarg_sites.json "$SITE"/methylarg_evidence.json "$SITE"/methylarg_interactions.json \
  "$SITE"/methylarg_pred_iface.json "$SITE_OUT/"
echo "[release] done: $OUT and $SITE_OUT/"

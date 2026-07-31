#!/bin/bash

set -eo pipefail

ROOT_DIR=${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
PYTHON_BIN=${PYTHON_BIN:-python3}
EXTENDED_PERMUTATIONS=${PTM_EXTENDED_PERMUTATIONS:-1000}
EXTENDED_2D_PERMUTATIONS=${PTM_EXTENDED_2D_PERMUTATIONS:-250}

. "$ROOT_DIR/hpc/ptm_env.sh"

export MPLCONFIGDIR=${MPLCONFIGDIR:-${TMPDIR:-/tmp}/ptm_mpl_extended}
export XDG_CACHE_HOME=${XDG_CACHE_HOME:-${TMPDIR:-/tmp}/ptm_cache_extended}
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME" "$PTM_RESULTS_ROOT"

run_step() {
  echo "[ptm_extended] $1"
  shift
  "$@"
}

require_file() {
  if [ ! -s "$1" ]; then
    echo "[ptm_extended] required file missing or empty: $1"
    exit 1
  fi
}

run_step "repo hygiene guard" bash "$PTM_CODE_ROOT/hpc/check_repo_hygiene.sh" "$PTM_CODE_ROOT"

require_file "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv"
require_file "$PTM_RESULTS_ROOT/domain_context/sites_with_domain_context.tsv"
require_file "$PTM_RESULTS_ROOT/domain_context/same_protein_nonmethyl_arginine_domain_background.tsv"
require_file "$PTM_RESULTS_ROOT/methyl_arg_clustering/per_site_nearest_methyl_arg_distance.tsv"
require_file "$PTM_RESULTS_ROOT/methyl_arg_clustering/nearest_neighbor_curve_vs_null.tsv"
require_file "$PTM_CANONICAL_FASTA"

NEIGHBOR_ARGS=(
  --integrated-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv"
  --canonical-fasta "$PTM_CANONICAL_FASTA"
  --outdir "$PTM_RESULTS_ROOT/methyl_arg_clustering"
  --permutations "$EXTENDED_PERMUTATIONS"
)
if [ -s "$PTM_MOBIDB_INTERVALS" ]; then
  NEIGHBOR_ARGS+=(--disorder-intervals "$PTM_MOBIDB_INTERVALS" --null-match-disorder)
fi
run_step "nearest-neighbor probability null" "$PYTHON_BIN" "$PTM_CODE_ROOT/scripts/17_neighbor_probability_null.py" "${NEIGHBOR_ARGS[@]}"

if [ -d "$PTM_RESULTS_ROOT/functional_class_union/clusterprofiler_go" ]; then
  run_step "format clusterProfiler GO panels" "$PYTHON_BIN" "$PTM_CODE_ROOT/scripts/18_format_clusterprofiler_outputs.py" \
    --indir "$PTM_RESULTS_ROOT/functional_class_union/clusterprofiler_go"
fi

run_step "format selected summary figures" "$PYTHON_BIN" "$PTM_CODE_ROOT/scripts/19_format_summary_figures.py" \
  --results-root "$PTM_RESULTS_ROOT"

run_step "representative methylarginine architecture figures" "$PYTHON_BIN" "$PTM_CODE_ROOT/scripts/20_draw_representative_methylarginine_architectures.py" \
  --site-table "$PTM_RESULTS_ROOT/domain_context/sites_with_domain_context.tsv" \
  --remapped-rows "$PTM_RESULTS_ROOT/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv" \
  --interpro "$PTM_INTERPRO_DOMAIN_INTERVALS" \
  --disorder "$PTM_MOBIDB_INTERVALS" \
  --metadata "$PTM_UNIPROT_METADATA" \
  --fasta "$PTM_CANONICAL_FASTA" \
  --outdir "$PTM_RESULTS_ROOT/summary_figures"

run_step "interaction redistribution concept figure" "$PYTHON_BIN" "$PTM_CODE_ROOT/scripts/21_draw_interaction_redistribution_model.py" \
  --outdir "$PTM_RESULTS_ROOT/summary_figures"

run_step "density and domain-adjacency schematic" "$PYTHON_BIN" "$PTM_CODE_ROOT/scripts/22_draw_density_and_adjacency_schematic.py" \
  --outdir "$PTM_RESULTS_ROOT/summary_figures"

run_step "domain-edge length-null analysis" "$PYTHON_BIN" "$PTM_CODE_ROOT/scripts/23_domain_edge_length_null.py" \
  --sites "$PTM_RESULTS_ROOT/domain_context/sites_with_domain_context.tsv" \
  --background "$PTM_RESULTS_ROOT/domain_context/same_protein_nonmethyl_arginine_domain_background.tsv" \
  --disorder "$PTM_MOBIDB_INTERVALS" \
  --interpro "$PTM_INTERPRO_DOMAIN_INTERVALS" \
  --outdir "$PTM_RESULTS_ROOT/domain_context"

run_step "plot domain-edge length-null analysis" "$PYTHON_BIN" "$PTM_CODE_ROOT/scripts/24_plot_domain_edge_length_null.py" \
  --outdir "$PTM_RESULTS_ROOT/domain_context"

for threshold in 3 4; do
  run_step "classify clustered domain architectures, threshold ${threshold}" "$PYTHON_BIN" "$PTM_CODE_ROOT/scripts/25_classify_domain_cluster_architectures.py" \
    --clusters "$PTM_RESULTS_ROOT/domain_context/domain_edge_class_20aa_cluster_examples.tsv" \
    --site-context "$PTM_RESULTS_ROOT/domain_context/sites_with_domain_context.tsv" \
    --outdir "$PTM_RESULTS_ROOT/domain_context" \
    --cluster-threshold "$threshold"

  run_step "plot signed domain-edge cluster distances, threshold ${threshold}" "$PYTHON_BIN" "$PTM_CODE_ROOT/scripts/26_plot_domain_cluster_signed_distance.py" \
    --clusters "$PTM_RESULTS_ROOT/domain_context/domain_edge_class_20aa_cluster_examples.tsv" \
    --site-context "$PTM_RESULTS_ROOT/domain_context/sites_with_domain_context.tsv" \
    --outdir "$PTM_RESULTS_ROOT/domain_context" \
    --cluster-threshold "$threshold"
done

run_step "signed domain-edge cumulative density" "$PYTHON_BIN" "$PTM_CODE_ROOT/scripts/27_plot_signed_domain_edge_cumulative_density.py" \
  --site-context "$PTM_RESULTS_ROOT/domain_context/sites_with_domain_context.tsv" \
  --length-null-sites "$PTM_RESULTS_ROOT/domain_context/domain_edge_idr_length_null_sites.tsv" \
  --background "$PTM_RESULTS_ROOT/domain_context/same_protein_nonmethyl_arginine_domain_background.tsv" \
  --interpro "$PTM_INTERPRO_DOMAIN_INTERVALS" \
  --outdir "$PTM_RESULTS_ROOT/domain_context"

run_step "multiscale methylarginine density null" "$PYTHON_BIN" "$PTM_CODE_ROOT/scripts/28_multiscale_methylarg_density_null.py" \
  --sites "$PTM_RESULTS_ROOT/domain_context/sites_with_domain_context.tsv" \
  --fasta "$PTM_CANONICAL_FASTA" \
  --outdir "$PTM_RESULTS_ROOT/methyl_arg_clustering" \
  --permutations "$EXTENDED_PERMUTATIONS"

run_step "nearest-neighbor density plot" "$PYTHON_BIN" "$PTM_CODE_ROOT/scripts/29_plot_nearest_neighbor_density.py" \
  --per-site "$PTM_RESULTS_ROOT/methyl_arg_clustering/per_site_nearest_methyl_arg_distance.tsv" \
  --curve "$PTM_RESULTS_ROOT/methyl_arg_clustering/nearest_neighbor_curve_vs_null.tsv" \
  --outdir "$PTM_RESULTS_ROOT/methyl_arg_clustering"

run_step "pairwise occupancy heatmap" "$PYTHON_BIN" "$PTM_CODE_ROOT/scripts/30_pairwise_occupancy_heatmap.py" \
  --sites "$PTM_RESULTS_ROOT/domain_context/sites_with_domain_context.tsv" \
  --fasta "$PTM_CANONICAL_FASTA" \
  --outdir "$PTM_RESULTS_ROOT/methyl_arg_clustering" \
  --permutations "$EXTENDED_2D_PERMUTATIONS"

run_step "k-nearest 2D density" "$PYTHON_BIN" "$PTM_CODE_ROOT/scripts/31_plot_knearest_2d_density.py" \
  --sites "$PTM_RESULTS_ROOT/domain_context/sites_with_domain_context.tsv" \
  --fasta "$PTM_CANONICAL_FASTA" \
  --outdir "$PTM_RESULTS_ROOT/methyl_arg_clustering" \
  --permutations "$EXTENDED_2D_PERMUTATIONS"

run_step "amino-acid-resolution k-nearest 2D density" "$PYTHON_BIN" "$PTM_CODE_ROOT/scripts/32_plot_knearest_aa_resolution_density.py" \
  --sites "$PTM_RESULTS_ROOT/domain_context/sites_with_domain_context.tsv" \
  --fasta "$PTM_CANONICAL_FASTA" \
  --outdir "$PTM_RESULTS_ROOT/methyl_arg_clustering" \
  --permutations "$EXTENDED_2D_PERMUTATIONS"

echo "[ptm_extended] extended-analysis outputs written under $PTM_RESULTS_ROOT"

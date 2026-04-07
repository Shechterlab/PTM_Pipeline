#!/bin/bash
#SBATCH -J ptm_integrate
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 16:00:00

set -euo pipefail

. /gs/gsfs0/hpc01/rhel8/apps/conda3/bin/activate
conda activate base

BASE_DIR=${1:-$PWD}
export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib-${SLURM_JOB_ID:-manual}"
mkdir -p "$MPLCONFIGDIR"
mkdir -p "$BASE_DIR/data/context"

python "$BASE_DIR/scripts/00_stage_context_sources.py" \
  --outdir "$BASE_DIR/data/context" \
  --skip-interpro
python "$BASE_DIR/scripts/01_parse_maron_s5.py" \
  --input "$BASE_DIR/data/external/Table_S5_Compiled_Methylarginine_Data.xlsx" \
  --outdir "$BASE_DIR/results/parsed_maron"
python "$BASE_DIR/scripts/01_parse_prometheus_mmc2.py" \
  --input "$BASE_DIR/data/external/mmc2.xlsx" \
  --outdir "$BASE_DIR/results/parsed_prometheus"
python "$BASE_DIR/scripts/02_remap_arg_methyl_to_canonical.py" \
  --base-master "$BASE_DIR/data/base/human_ptm_master.tsv" \
  --maron "$BASE_DIR/results/parsed_maron/maron_s5_arg_methyl_sites.tsv" \
  --prometheus "$BASE_DIR/results/parsed_prometheus/prometheus_mmc2_arg_methyl_sites.tsv" \
  --canonical-fasta "$BASE_DIR/data/context/uniprot_human_reviewed_canonical.fasta" \
  --uniprot-metadata "$BASE_DIR/data/context/uniprot_human_reviewed_canonical_metadata.tsv" \
  --outdir "$BASE_DIR/results/remapped"
python "$BASE_DIR/scripts/02_integrate_arg_methyl_sources.py" \
  --premerged-remapped "$BASE_DIR/results/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv" \
  --outdir "$BASE_DIR/results/integrated"
python "$BASE_DIR/scripts/03_functional_class_enrichment.py" \
  --base-master "$BASE_DIR/data/base/human_ptm_master.tsv" \
  --integrated-sites "$BASE_DIR/results/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --ontology "$BASE_DIR/config/functional_ontology.json" \
  --outdir "$BASE_DIR/results/functional_class_union"
python "$BASE_DIR/scripts/04_example_figures.py" \
  --integrated-sites "$BASE_DIR/results/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --outdir "$BASE_DIR/results/example_figures"

if [[ -f "$BASE_DIR/data/context/protein2ipr.dat.gz" ]]; then
  python "$BASE_DIR/scripts/00_stage_context_sources.py" \
    --outdir "$BASE_DIR/data/context" \
    --accessions "$BASE_DIR/results/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv" \
    --interpro-source bulk \
    --interpro-protein2ipr "$BASE_DIR/data/context/protein2ipr.dat.gz" \
    --bulk-parser pandas
  python "$BASE_DIR/scripts/05_annotate_domain_context.py" \
    --sites "$BASE_DIR/results/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
    --interpro-intervals "$BASE_DIR/data/context/interpro_human_reviewed_domain_like_intervals.tsv" \
    --position-col corrected_position \
    --outdir "$BASE_DIR/results/domain_context"
else
  echo "Skipping InterPro domain-context staging: $BASE_DIR/data/context/protein2ipr.dat.gz not found"
fi

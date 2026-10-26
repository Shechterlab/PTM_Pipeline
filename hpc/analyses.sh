# The analyses of a release, on the site union WORK/integrated/human_arg_methyl_union_dedup_by_site.tsv (sourced by
# hpc/build_release.sh and hpc/release_analyses.sh, so every release is analysed with the same steps). Needs CODE, IN,
# FASTA, ENTRIES and WORK, and python and Rscript on PATH; each step's output goes to WORK/<folder>, its log to
# WORK/logs/.

# Stops unless every input MANIFEST.tsv lists as in Git matches its MD5.
check_inputs() {
  python - "$IN" <<'PYEOF'
import hashlib, sys
from pathlib import Path
import pandas as pd
inputs = Path(sys.argv[1])
manifest = pd.read_csv(inputs / 'MANIFEST.tsv', sep='\t', dtype=str, keep_default_na=False)
bad = [r.file for r in manifest[manifest['redistribution'] == 'git'].itertuples()
       if not (inputs / r.file).exists() or hashlib.md5((inputs / r.file).read_bytes()).hexdigest() != r.md5]
if bad:
    sys.exit(f'inputs differ from MANIFEST.tsv: {bad}')
PYEOF
}

run_analyses() {
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
  # scripts/25 derives its cluster examples from the site context only when the file is absent: remove it, or a second
  # run in the same WORK_ROOT reads the sites of the first
  rm -f "$WORK/domain_context/domain_edge_class_20aa_cluster_examples.tsv"
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

}

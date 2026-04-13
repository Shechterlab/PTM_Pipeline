#!/bin/bash

ptm_choose_default_root() {
  primary="$1"
  secondary="$2"
  legacy="$3"
  if [ -n "$primary" ] && [ -e "$primary" ]; then
    printf '%s\n' "$primary"
    return
  fi
  if [ -n "$secondary" ] && [ -e "$secondary" ]; then
    printf '%s\n' "$secondary"
    return
  fi
  printf '%s\n' "$legacy"
}

PTM_ENV_SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PTM_CODE_ROOT=${PTM_CODE_ROOT:-$(cd "$PTM_ENV_SCRIPT_DIR/.." && pwd)}
PTM_PROJECT_ROOT=${PTM_PROJECT_ROOT:-$(cd "$PTM_CODE_ROOT/.." && pwd)}

DEFAULT_DATA_ROOT=$(ptm_choose_default_root "$PTM_PROJECT_ROOT/PTM_data" "$PTM_CODE_ROOT/PTM_data" "$PTM_CODE_ROOT/data")
DEFAULT_RESULTS_ROOT=$(ptm_choose_default_root "$PTM_PROJECT_ROOT/PTM_results" "$PTM_CODE_ROOT/PTM_results" "$PTM_CODE_ROOT/results")

PTM_DATA_ROOT=${PTM_DATA_ROOT:-$DEFAULT_DATA_ROOT}
PTM_RESULTS_ROOT=${PTM_RESULTS_ROOT:-$DEFAULT_RESULTS_ROOT}

PTM_BASE_ROOT=${PTM_BASE_ROOT:-$PTM_DATA_ROOT/base}
PTM_EXTERNAL_ROOT=${PTM_EXTERNAL_ROOT:-$PTM_DATA_ROOT/external}
PTM_CONTEXT_ROOT=${PTM_CONTEXT_ROOT:-$PTM_DATA_ROOT/context}

PTM_MASTER_TSV=${PTM_MASTER_TSV:-$PTM_BASE_ROOT/human_ptm_master.tsv}
PTM_MARON_XLSX=${PTM_MARON_XLSX:-$PTM_EXTERNAL_ROOT/Table_S5_Compiled_Methylarginine_Data.xlsx}
PTM_PROMETHEUS_MMC2=${PTM_PROMETHEUS_MMC2:-$PTM_EXTERNAL_ROOT/mmc2.xlsx}
PTM_MOBIDB_JSON=${PTM_MOBIDB_JSON:-$PTM_EXTERNAL_ROOT/mobidb_human_2026-03-30.json}
PTM_CANONICAL_FASTA=${PTM_CANONICAL_FASTA:-$PTM_CONTEXT_ROOT/uniprot_human_reviewed_canonical.fasta}
PTM_UNIPROT_METADATA=${PTM_UNIPROT_METADATA:-$PTM_CONTEXT_ROOT/uniprot_human_reviewed_canonical_metadata.tsv}
PTM_MOBIDB_INTERVALS=${PTM_MOBIDB_INTERVALS:-$PTM_CONTEXT_ROOT/mobidb_human_reviewed_disorder_intervals.tsv}
PTM_INTERPRO_BULK=${PTM_INTERPRO_BULK:-$PTM_CONTEXT_ROOT/protein2ipr.dat.gz}
PTM_INTERPRO_ALL_INTERVALS=${PTM_INTERPRO_ALL_INTERVALS:-$PTM_CONTEXT_ROOT/interpro_human_reviewed_all_entry_intervals.tsv}
PTM_INTERPRO_DOMAIN_INTERVALS=${PTM_INTERPRO_DOMAIN_INTERVALS:-$PTM_CONTEXT_ROOT/interpro_human_reviewed_domain_like_intervals.tsv}
PTM_CDCODE_STAGED=${PTM_CDCODE_STAGED:-$(ptm_choose_default_root "$PTM_DATA_ROOT/condensates/cdcode_staged_membership.tsv" "$PTM_DATA_ROOT/condensates/protein2cdcode_v2.2.tsv" "$PTM_DATA_ROOT/condensates/cdcode_staged_membership.tsv")}
PTM_CDCODE_PROTEINS=${PTM_CDCODE_PROTEINS:-$(ptm_choose_default_root "$PTM_DATA_ROOT/condensates/proteins_202603181653.csv" "" "$PTM_DATA_ROOT/condensates/proteins_202603181653.csv")}
PTM_CDCODE_CONDENSATES=${PTM_CDCODE_CONDENSATES:-$(ptm_choose_default_root "$PTM_DATA_ROOT/condensates/condensates_202603181647.csv" "" "$PTM_DATA_ROOT/condensates/condensates_202603181647.csv")}

if [ -z "$PTM_CODE_ROOT" ] || [ -z "$PTM_DATA_ROOT" ] || [ -z "$PTM_RESULTS_ROOT" ]; then
  echo "[ptm_env] PTM_CODE_ROOT, PTM_DATA_ROOT, and PTM_RESULTS_ROOT must not be empty"
  return 1 2>/dev/null || exit 1
fi

mkdir -p "$PTM_DATA_ROOT" "$PTM_RESULTS_ROOT" "$PTM_CONTEXT_ROOT"

export PTM_CODE_ROOT
export PTM_PROJECT_ROOT
export PTM_DATA_ROOT
export PTM_RESULTS_ROOT
export PTM_BASE_ROOT
export PTM_EXTERNAL_ROOT
export PTM_CONTEXT_ROOT
export PTM_MASTER_TSV
export PTM_MARON_XLSX
export PTM_PROMETHEUS_MMC2
export PTM_MOBIDB_JSON
export PTM_CANONICAL_FASTA
export PTM_UNIPROT_METADATA
export PTM_MOBIDB_INTERVALS
export PTM_INTERPRO_BULK
export PTM_INTERPRO_ALL_INTERVALS
export PTM_INTERPRO_DOMAIN_INTERVALS
export PTM_CDCODE_STAGED
export PTM_CDCODE_PROTEINS
export PTM_CDCODE_CONDENSATES

if [ "${PTM_DEBUG_PATHS:-0}" = "1" ]; then
  echo "[ptm_env] PTM_CODE_ROOT=$PTM_CODE_ROOT"
  echo "[ptm_env] PTM_DATA_ROOT=$PTM_DATA_ROOT"
  echo "[ptm_env] PTM_RESULTS_ROOT=$PTM_RESULTS_ROOT"
  echo "[ptm_env] PTM_CONTEXT_ROOT=$PTM_CONTEXT_ROOT"
fi

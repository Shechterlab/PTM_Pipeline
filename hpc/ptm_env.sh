#!/bin/bash

export PTM_PROJECT_ROOT="/gs/gsfs0/home/dshecht1/projects/2026-04_ArgReview"
export PTM_DATA_ROOT="${PTM_DATA_ROOT:-$PTM_PROJECT_ROOT/PTM_data}"
export PTM_RESULTS_ROOT="${PTM_RESULTS_ROOT:-$PTM_PROJECT_ROOT/PTM_results}"

export PTM_BASE_MASTER="${PTM_BASE_MASTER:-$PTM_DATA_ROOT/base/human_ptm_master.tsv}"
export PTM_MARON_XLSX="${PTM_MARON_XLSX:-$PTM_DATA_ROOT/external/Table_S5_Compiled_Methylarginine_Data.xlsx}"
export PTM_PROMETHEUS_XLSX="${PTM_PROMETHEUS_XLSX:-$PTM_DATA_ROOT/external/ProMetheus_Arg_Methyl.xlsx}"

export PTM_CANONICAL_FASTA="${PTM_CANONICAL_FASTA:-$PTM_DATA_ROOT/context/uniprot_human_reviewed_canonical.fasta}"
export PTM_CANONICAL_METADATA="${PTM_CANONICAL_METADATA:-$PTM_DATA_ROOT/context/uniprot_human_reviewed_canonical_metadata.tsv}"
export PTM_DISORDER_INTERVALS="${PTM_DISORDER_INTERVALS:-$PTM_DATA_ROOT/context/mobidb_human_reviewed_disorder_intervals.tsv}"
export PTM_INTERPRO_INTERVALS="${PTM_INTERPRO_INTERVALS:-$PTM_DATA_ROOT/context/interpro_human_reviewed_domain_like_intervals.tsv}"
export PTM_INTERPRO_PROTEIN2IPR="${PTM_INTERPRO_PROTEIN2IPR:-$PTM_DATA_ROOT/context/protein2ipr.dat.gz}"
export PTM_INTERPRO_ENTRY_LIST="${PTM_INTERPRO_ENTRY_LIST:-$PTM_DATA_ROOT/context/entry.list}"

export PTM_CDCODE_PROTEINS="${PTM_CDCODE_PROTEINS:-$PTM_DATA_ROOT/external/cd_code/proteins_202603181653.csv}"
export PTM_CDCODE_PROTEIN_MAP="${PTM_CDCODE_PROTEIN_MAP:-$PTM_DATA_ROOT/external/cd_code/protein2cdcode_v2.2.tsv}"
export PTM_CDCODE_CONDENSATES="${PTM_CDCODE_CONDENSATES:-$PTM_DATA_ROOT/external/cd_code/condensates_202603181647.csv}"
export PTM_CDCODE_STAGED="${PTM_CDCODE_STAGED:-$PTM_DATA_ROOT/context/cd_code_condensate_proteins.tsv}"

export PTM_MARON_EXCLUDE_PATTERN="${PTM_MARON_EXCLUDE_PATTERN:-}"
export PTM_METHYL_CLUSTER_PERMUTATIONS="${PTM_METHYL_CLUSTER_PERMUTATIONS:-100000}"
export PTM_RUN_CONDENSATE="${PTM_RUN_CONDENSATE:-1}"
export PTM_RUN_CROSS_PTM="${PTM_RUN_CROSS_PTM:-0}"

export PTM_CONDENSATE_EXPERIMENTAL_ONLY="${PTM_CONDENSATE_EXPERIMENTAL_ONLY:-0}"
export PTM_CONDENSATE_MIN_CONFIDENCE="${PTM_CONDENSATE_MIN_CONFIDENCE:-}"
export PTM_CONDA_ENV="${PTM_CONDA_ENV:-}"
export PTM_PYTHON_BIN="${PTM_PYTHON_BIN:-python}"

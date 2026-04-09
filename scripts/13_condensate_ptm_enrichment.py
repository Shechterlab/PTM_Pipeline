from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import ConvergenceWarning, PerfectSeparationWarning

from common import add_bh_q_values, canonicalize_uniprot_accession, fisher_like_enrichment, parse_fasta, save_figure, save_table

DEFAULT_PTM_GROUPS = [
    'Arg methylation',
    'Phosphorylation',
    'Lys acylation',
    'Lys methylation',
    'Ubiquitin/SUMO',
]

PTM_GROUP_TARGET_RESIDUES = {
    'Arg methylation': {'R'},
    'Phosphorylation': {'S', 'T', 'Y'},
    'Lys acylation': {'K'},
    'Lys methylation': {'K'},
    'Ubiquitin/SUMO': {'K'},
}

ROLE_PRIORITY = {
    'driver': 0,
    'scaffold': 1,
    'member': 2,
    'client': 3,
    'unspecified': 4,
    'non_condensate': 5,
}



def normalize_master_sites(master: pd.DataFrame, sequences: dict[str, str], ptm_groups: list[str]) -> pd.DataFrame:
    out = master[master['ptm_group'].isin(ptm_groups)].copy()
    out['canonical_UniProtAC'] = out['substrate_UniProtAC'].astype(str).map(canonicalize_uniprot_accession)
    out['position'] = pd.to_numeric(out['position'], errors='coerce')
    out = out[out['position'].notna()].copy()
    out['position'] = out['position'].astype(int)
    out['residue'] = out['residue'].astype(str).str.upper().str[:1]
    out = out[
        out.apply(
            lambda row: row['residue'] in PTM_GROUP_TARGET_RESIDUES.get(str(row['ptm_group']), {str(row['residue'])}),
            axis=1,
        )
    ].copy()

    def valid_row(row: pd.Series) -> bool:
        seq = sequences.get(str(row['canonical_UniProtAC']), '')
        pos = int(row['position'])
        return bool(seq) and 1 <= pos <= len(seq) and seq[pos - 1] == row['residue']

    out = out[out.apply(valid_row, axis=1)].copy()
    out['site'] = out['residue'] + out['position'].astype(str)
    out['site_key'] = out['canonical_UniProtAC'] + ':' + out['site']
    return out.drop_duplicates(['ptm_group', 'canonical_UniProtAC', 'residue', 'position']).copy()



def normalize_arg_methyl_sites(df: pd.DataFrame, sequences: dict[str, str]) -> pd.DataFrame:
    out = df.copy()
    out['canonical_UniProtAC'] = out['canonical_UniProtAC'].astype(str).map(canonicalize_uniprot_accession)
    out['position'] = pd.to_numeric(out['corrected_position'], errors='coerce')
    out = out[out['position'].notna()].copy()
    out['position'] = out['position'].astype(int)
    out['residue'] = 'R'
    out['ptm_group'] = 'Arg methylation'

    def valid_row(row: pd.Series) -> bool:
        seq = sequences.get(str(row['canonical_UniProtAC']), '')
        pos = int(row['position'])
        return bool(seq) and 1 <= pos <= len(seq) and seq[pos - 1] == 'R'

    out = out[out.apply(valid_row, axis=1)].copy()
    out['site'] = out['residue'] + out['position'].astype(str)
    out['site_key'] = out['canonical_UniProtAC'] + ':' + out['site']
    return out.drop_duplicates(['ptm_group', 'canonical_UniProtAC', 'residue', 'position']).copy()



def low_complexity_fraction(sequence: str, window: int = 12, max_unique: int = 3) -> float:
    seq = str(sequence).strip().upper()
    if not seq:
        return np.nan
    if len(seq) <= window:
        return float(len(set(seq)) <= max_unique)
    mask = np.zeros(len(seq), dtype=bool)
    for start in range(0, len(seq) - window + 1):
        chunk = seq[start:start + window]
        if len(set(chunk)) <= max_unique:
            mask[start:start + window] = True
    return float(mask.mean())



def disorder_fraction_table(disorder: pd.DataFrame, sequences: dict[str, str]) -> pd.DataFrame:
    if disorder.empty:
        return pd.DataFrame(columns=['canonical_UniProtAC', 'disorder_fraction'])
    rows = []
    for accession, group in disorder.groupby('canonical_UniProtAC'):
        length = len(sequences.get(accession, ''))
        if length == 0:
            continue
        covered = np.zeros(length, dtype=bool)
        for row in group.itertuples(index=False):
            start = max(1, int(getattr(row, 'fragment_start'))) - 1
            end = min(length, int(getattr(row, 'fragment_end')))
            covered[start:end] = True
        rows.append({'canonical_UniProtAC': accession, 'disorder_fraction': float(covered.mean())})
    return pd.DataFrame(rows)



def parse_bool(value) -> bool:
    text = str(value).strip().lower()
    return text in {'true', '1', 'yes', 'y', 't'}



def normalize_role(value: str) -> str:
    text = str(value).strip().lower()
    if text in {'', 'nan', 'none'}:
        return 'unspecified'
    if 'driver' in text:
        return 'driver'
    if 'scaffold' in text:
        return 'scaffold'
    if 'member' in text:
        return 'member'
    if 'client' in text:
        return 'client'
    return text.replace(' ', '_')



def choose_role(values: pd.Series, is_condensate: bool) -> str:
    roles = sorted({normalize_role(value) for value in values if str(value).strip() not in {'', 'nan', 'None'}})
    if not roles:
        return 'unspecified' if is_condensate else 'non_condensate'
    return sorted(roles, key=lambda role: ROLE_PRIORITY.get(role, 99))[0]



def parse_condensate_table(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    rename_map = {}
    for column in out.columns:
        key = str(column).strip().lower()
        if key in {'uniprot', 'uniprotac', 'uniprot_acc', 'accession', 'protein_accession', 'canonical_uniprotac'}:
            rename_map[column] = 'canonical_UniProtAC'
        elif key in {'is_condensate', 'condensate', 'condensate_associated'}:
            rename_map[column] = 'is_condensate'
        elif key in {'driver_or_member', 'functional_type', 'role'}:
            rename_map[column] = 'condensate_role'
        elif key in {'confidence_score', 'confidence'}:
            rename_map[column] = 'confidence_score'
        elif key in {'is_experimental', 'experimental'}:
            rename_map[column] = 'is_experimental'
    out = out.rename(columns=rename_map)
    if 'canonical_UniProtAC' not in out.columns:
        raise ValueError('condensate table must include a UniProt accession column')
    out['canonical_UniProtAC'] = out['canonical_UniProtAC'].astype(str).map(canonicalize_uniprot_accession)
    out = out[out['canonical_UniProtAC'].ne('')].copy()
    if 'is_condensate' not in out.columns:
        out['is_condensate'] = True
    out['is_condensate'] = out['is_condensate'].map(parse_bool)
    if 'condensate_role' not in out.columns:
        out['condensate_role'] = 'unspecified'
    out['condensate_role'] = out['condensate_role'].map(normalize_role)
    if 'is_experimental' in out.columns:
        out['is_experimental'] = out['is_experimental'].map(parse_bool)
    else:
        out['is_experimental'] = False
    if 'confidence_score' in out.columns:
        out['confidence_score'] = pd.to_numeric(out['confidence_score'], errors='coerce')
    else:
        out['confidence_score'] = np.nan

    agg = out.groupby('canonical_UniProtAC', as_index=False).agg(
        is_condensate=('is_condensate', 'max'),
        is_experimental=('is_experimental', 'max'),
        confidence_score=('confidence_score', 'max'),
        condensate_role=('condensate_role', lambda s: choose_role(s, bool(np.any(pd.Series(s).map(lambda x: normalize_role(x) != 'non_condensate'))))),
    )
    agg.loc[~agg['is_condensate'], 'condensate_role'] = 'non_condensate'
    return agg



def build_protein_table(
    sequences: dict[str, str],
    ptm_sites: pd.DataFrame,
    condensate_table: pd.DataFrame,
    disorder_fraction: pd.DataFrame,
) -> pd.DataFrame:
    proteins = pd.DataFrame({'canonical_UniProtAC': sorted(sequences)})
    proteins['protein_length'] = proteins['canonical_UniProtAC'].map(lambda acc: len(sequences.get(acc, '')))
    proteins['low_complexity_fraction'] = proteins['canonical_UniProtAC'].map(lambda acc: low_complexity_fraction(sequences.get(acc, '')))
    proteins = proteins.merge(condensate_table, how='left', on='canonical_UniProtAC')
    proteins['is_condensate'] = proteins['is_condensate'].where(proteins['is_condensate'].notna(), False).astype(bool)
    proteins['is_experimental'] = proteins['is_experimental'].where(proteins['is_experimental'].notna(), False).astype(bool)
    proteins['condensate_role'] = proteins['condensate_role'].fillna('non_condensate')
    proteins['confidence_score'] = pd.to_numeric(proteins['confidence_score'], errors='coerce')
    proteins = proteins.merge(disorder_fraction, how='left', on='canonical_UniProtAC')
    proteins['disorder_fraction'] = pd.to_numeric(proteins['disorder_fraction'], errors='coerce').fillna(0.0)

    for ptm_group, residues in PTM_GROUP_TARGET_RESIDUES.items():
        proteins[f'candidate_residue_count__{ptm_group}'] = proteins['canonical_UniProtAC'].map(
            lambda acc: sum(sequences.get(acc, '').count(residue) for residue in residues)
        )
    site_counts = ptm_sites.groupby(['ptm_group', 'canonical_UniProtAC']).size().rename('site_count').reset_index()
    for ptm_group in sorted(PTM_GROUP_TARGET_RESIDUES):
        group_counts = site_counts[site_counts['ptm_group'] == ptm_group][['canonical_UniProtAC', 'site_count']].copy()
        group_counts = group_counts.rename(columns={'site_count': f'site_count__{ptm_group}'})
        proteins = proteins.merge(group_counts, how='left', on='canonical_UniProtAC')
        proteins[f'site_count__{ptm_group}'] = proteins[f'site_count__{ptm_group}'].fillna(0).astype(int)
        proteins[f'has_ptm__{ptm_group}'] = proteins[f'site_count__{ptm_group}'] > 0
        candidate = proteins[f'candidate_residue_count__{ptm_group}'].replace(0, np.nan)
        proteins[f'site_density_per_candidate__{ptm_group}'] = proteins[f'site_count__{ptm_group}'] / candidate
    return proteins



def fit_adjusted_model(proteins: pd.DataFrame, ptm_group: str) -> dict:
    feature_cols = [
        f'has_ptm__{ptm_group}',
        'protein_length_log10',
        f'candidate_residue_count_log10__{ptm_group}',
        'disorder_fraction',
        'low_complexity_fraction',
    ]
    model_df = proteins[['is_condensate'] + feature_cols].dropna().copy()
    if model_df.empty:
        return {'ptm_group': ptm_group, 'model_status': 'no_rows'}
    X = sm.add_constant(model_df[feature_cols].astype(float), has_constant='add')
    y = model_df['is_condensate'].astype(float)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=PerfectSeparationWarning)
            warnings.simplefilter('ignore', category=ConvergenceWarning)
            warnings.simplefilter('ignore', category=RuntimeWarning)
            fit = sm.Logit(y, X).fit(disp=False, maxiter=300)
        coef = fit.params[f'has_ptm__{ptm_group}']
        conf = fit.conf_int().loc[f'has_ptm__{ptm_group}']
        return {
            'ptm_group': ptm_group,
            'model_status': 'ok',
            'model_type': 'logit',
            'n_rows_model': int(len(model_df)),
            'adjusted_log_odds': float(coef),
            'adjusted_odds_ratio': float(np.exp(coef)),
            'adjusted_ci_low': float(np.exp(conf.iloc[0])),
            'adjusted_ci_high': float(np.exp(conf.iloc[1])),
            'adjusted_p_value': float(fit.pvalues[f'has_ptm__{ptm_group}']),
        }
    except Exception as exc:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', category=PerfectSeparationWarning)
                warnings.simplefilter('ignore', category=ConvergenceWarning)
                warnings.simplefilter('ignore', category=RuntimeWarning)
                fit = sm.Logit(y, X).fit_regularized(alpha=1e-6, L1_wt=0.0, disp=False, maxiter=500)
            coef = float(fit.params[f'has_ptm__{ptm_group}'])
            return {
                'ptm_group': ptm_group,
                'model_status': f'penalized_fallback:{type(exc).__name__}',
                'model_type': 'logit_regularized',
                'n_rows_model': int(len(model_df)),
                'adjusted_log_odds': coef,
                'adjusted_odds_ratio': float(np.exp(coef)),
                'adjusted_ci_low': np.nan,
                'adjusted_ci_high': np.nan,
                'adjusted_p_value': np.nan,
            }
        except Exception as exc2:
            return {
                'ptm_group': ptm_group,
                'model_status': f'failed:{type(exc2).__name__}',
                'model_type': 'none',
                'n_rows_model': int(len(model_df)),
                'adjusted_log_odds': np.nan,
                'adjusted_odds_ratio': np.nan,
                'adjusted_ci_low': np.nan,
                'adjusted_ci_high': np.nan,
                'adjusted_p_value': np.nan,
            }



def plot_adjusted_models(adjusted: pd.DataFrame, outpath: Path) -> None:
    if 'adjusted_odds_ratio' not in adjusted.columns:
        return
    plot_df = adjusted[adjusted['adjusted_odds_ratio'].notna()].copy()
    if plot_df.empty:
        return
    plot_df = plot_df.sort_values('adjusted_odds_ratio')
    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    y = np.arange(len(plot_df))
    ax.scatter(plot_df['adjusted_odds_ratio'], y)
    ci_mask = plot_df['adjusted_ci_low'].notna() & plot_df['adjusted_ci_high'].notna()
    if ci_mask.any():
        ax.hlines(
            y[ci_mask.to_numpy()],
            plot_df.loc[ci_mask, 'adjusted_ci_low'],
            plot_df.loc[ci_mask, 'adjusted_ci_high'],
        )
    ax.axvline(1.0, color='black', linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df['ptm_group'])
    ax.set_xlabel('Adjusted odds ratio for condensate association')
    ax.set_ylabel('PTM class')
    ax.set_title('Condensate association after covariate adjustment')
    for yi, row in zip(y, plot_df.itertuples(index=False)):
        q = getattr(row, 'adjusted_q_value_bh', np.nan)
        label = f'q={q:.3g}' if pd.notna(q) else str(row.model_status)
        ax.text(float(row.adjusted_odds_ratio), yi, f'  {label}', va='center', fontsize=8)
    fig.tight_layout()
    save_figure(fig, outpath)




def read_auto_table(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine='python')


def detect_column(df: pd.DataFrame, aliases: set[str]) -> str | None:
    lowered = {str(column).strip().lower(): column for column in df.columns}
    for alias in aliases:
        if alias in lowered:
            return lowered[alias]
    for key, original in lowered.items():
        if any(alias in key for alias in aliases):
            return original
    return None



def split_multi_value_cell(value) -> list[str]:
    text = str(value).strip()
    if text in {'', 'nan', 'None', 'null'}:
        return []
    parts = [part.strip() for part in __import__('re').split(r'[;,|]', text)]
    return [part for part in parts if part]



def build_condensate_table_from_cd_code_raw(
    proteins_path: str,
    protein_map_path: str = '',
    condensates_path: str = '',
) -> pd.DataFrame:
    proteins_raw = read_auto_table(proteins_path)
    proteins = proteins_raw.copy()

    acc_col = detect_column(proteins, {
        'uniprot', 'uniprotac', 'uniprot_acc', 'accession', 'protein_accession', 'canonical_uniprotac'
    })
    if acc_col is None:
        raise ValueError('CD-CODE proteins file must include a UniProt accession-like column')
    proteins['canonical_UniProtAC'] = proteins[acc_col].astype(str).map(canonicalize_uniprot_accession)
    proteins = proteins[proteins['canonical_UniProtAC'].ne('')].copy()

    role_col = detect_column(proteins, {'functional_type', 'driver_or_member', 'role'})
    exp_col = detect_column(proteins, {'is_experimental', 'experimental'})
    conf_col = detect_column(proteins, {'confidence_score', 'confidence'})
    cdcode_id_col = detect_column(proteins, {'cdcode_id', 'condensate_id', 'condensate_ids', 'condensates'})

    direct = proteins[['canonical_UniProtAC']].copy()
    direct['is_condensate'] = True
    direct['condensate_role'] = proteins[role_col].map(normalize_role) if role_col else 'unspecified'
    direct['is_experimental'] = proteins[exp_col].map(parse_bool) if exp_col else False
    direct['confidence_score'] = pd.to_numeric(proteins[conf_col], errors='coerce') if conf_col else np.nan
    direct_rows = [direct]

    if protein_map_path:
        mapping_raw = read_auto_table(protein_map_path)
        mapping = mapping_raw.copy()
        map_acc_col = detect_column(mapping, {
            'uniprot', 'uniprotac', 'uniprot_acc', 'accession', 'protein_accession', 'canonical_uniprotac'
        })
        map_cdcode_col = detect_column(mapping, {'cdcode_id', 'condensate_id', 'condensate_ids', 'cdcode', 'condensate'})
        if map_acc_col is None or map_cdcode_col is None:
            raise ValueError('protein2cdcode mapping must include both UniProt accession and condensate/CD-CODE ID columns')
        mapping['canonical_UniProtAC'] = mapping[map_acc_col].astype(str).map(canonicalize_uniprot_accession)
        mapping = mapping[mapping['canonical_UniProtAC'].ne('')].copy()
        mapping['condensate_id'] = mapping[map_cdcode_col].astype(str)
        mapping['condensate_id'] = mapping['condensate_id'].map(split_multi_value_cell)
        mapping = mapping.explode('condensate_id')
        mapping['condensate_id'] = mapping['condensate_id'].astype(str).str.strip()
        mapping = mapping[mapping['condensate_id'].ne('')].copy()

        edge = mapping[['canonical_UniProtAC', 'condensate_id']].drop_duplicates().copy()
        edge['is_condensate'] = True
        edge['condensate_role'] = 'unspecified'
        edge['is_experimental'] = False
        edge['confidence_score'] = np.nan

        if condensates_path:
            condensates_raw = read_auto_table(condensates_path)
            condensates = condensates_raw.copy()
            condensate_id_col = detect_column(condensates, {'cdcode_id', 'condensate_id', 'id'})
            if condensate_id_col is None:
                raise ValueError('condensates file must include a condensate/CD-CODE ID column')
            rename_map = {condensate_id_col: 'condensate_id'}
            role2 = detect_column(condensates, {'functional_type', 'driver_or_member', 'role'})
            exp2 = detect_column(condensates, {'is_experimental', 'experimental'})
            conf2 = detect_column(condensates, {'confidence_score', 'confidence'})
            if role2 is not None:
                rename_map[role2] = 'condensate_role'
            if exp2 is not None:
                rename_map[exp2] = 'is_experimental'
            if conf2 is not None:
                rename_map[conf2] = 'confidence_score'
            condensates = condensates.rename(columns=rename_map)
            keep_cols = ['condensate_id']
            for col in ['condensate_role', 'is_experimental', 'confidence_score']:
                if col in condensates.columns:
                    keep_cols.append(col)
            condensates = condensates[keep_cols].copy()
            if 'condensate_role' in condensates.columns:
                condensates['condensate_role'] = condensates['condensate_role'].map(normalize_role)
            if 'is_experimental' in condensates.columns:
                condensates['is_experimental'] = condensates['is_experimental'].map(parse_bool)
            if 'confidence_score' in condensates.columns:
                condensates['confidence_score'] = pd.to_numeric(condensates['confidence_score'], errors='coerce')
            edge = edge.merge(condensates, how='left', on='condensate_id', suffixes=('', '_meta'))
            for col in ['condensate_role', 'is_experimental', 'confidence_score']:
                meta_col = f'{col}_meta'
                if meta_col in edge.columns:
                    edge[col] = edge[meta_col].where(edge[meta_col].notna(), edge[col])
                    edge = edge.drop(columns=[meta_col])

        direct_rows.append(edge[['canonical_UniProtAC', 'is_condensate', 'condensate_role', 'is_experimental', 'confidence_score']])

    combined = pd.concat(direct_rows, ignore_index=True, sort=False)
    return parse_condensate_table(combined)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--base-master', required=True)
    ap.add_argument('--integrated-arg-sites', required=True)
    ap.add_argument('--canonical-fasta', required=True)
    ap.add_argument('--condensate-proteins', default='')
    ap.add_argument('--cdcode-proteins', default='')
    ap.add_argument('--cdcode-protein-map', default='')
    ap.add_argument('--cdcode-condensates', default='')
    ap.add_argument('--disorder-intervals', default='')
    ap.add_argument('--outdir', default='results/condensate_enrichment')
    ap.add_argument('--ptm-groups', nargs='*', default=list(DEFAULT_PTM_GROUPS))
    ap.add_argument('--experimental-only', action='store_true')
    ap.add_argument('--min-confidence-score', type=float, default=np.nan)
    args = ap.parse_args()

    sequences = parse_fasta(args.canonical_fasta)
    master = pd.read_csv(args.base_master, sep='\t', low_memory=False)
    arg_sites_in = pd.read_csv(args.integrated_arg_sites, sep='\t', low_memory=False)
    if args.condensate_proteins:
        condensate_table = parse_condensate_table(read_auto_table(args.condensate_proteins))
    elif args.cdcode_proteins:
        condensate_table = build_condensate_table_from_cd_code_raw(
            proteins_path=args.cdcode_proteins,
            protein_map_path=args.cdcode_protein_map,
            condensates_path=args.cdcode_condensates,
        )
    else:
        raise ValueError('provide either --condensate-proteins or --cdcode-proteins')
    if args.experimental_only and 'is_experimental' in condensate_table.columns:
        condensate_table = condensate_table[condensate_table['is_experimental'].fillna(False)].copy()
    if pd.notna(args.min_confidence_score) and 'confidence_score' in condensate_table.columns:
        keep = condensate_table['confidence_score'].fillna(-np.inf) >= float(args.min_confidence_score)
        condensate_table = condensate_table[keep].copy()
    disorder_fraction = pd.DataFrame(columns=['canonical_UniProtAC', 'disorder_fraction'])
    if args.disorder_intervals:
        disorder = pd.read_csv(args.disorder_intervals, sep='\t', low_memory=False)
        disorder['canonical_UniProtAC'] = disorder['canonical_UniProtAC'].astype(str).map(canonicalize_uniprot_accession)
        disorder_fraction = disorder_fraction_table(disorder, sequences=sequences)

    other_groups = [group for group in args.ptm_groups if group != 'Arg methylation']
    master_sites = normalize_master_sites(master, sequences=sequences, ptm_groups=other_groups)
    arg_sites = normalize_arg_methyl_sites(arg_sites_in, sequences=sequences)
    ptm_sites = pd.concat([master_sites, arg_sites], ignore_index=True, sort=False)
    ptm_sites = ptm_sites[ptm_sites['ptm_group'].isin(args.ptm_groups)].copy()

    proteins = build_protein_table(
        sequences=sequences,
        ptm_sites=ptm_sites,
        condensate_table=condensate_table,
        disorder_fraction=disorder_fraction,
    )
    proteins['protein_length_log10'] = np.log10(proteins['protein_length'].clip(lower=1))
    for ptm_group in args.ptm_groups:
        proteins[f'candidate_residue_count_log10__{ptm_group}'] = np.log10(proteins[f'candidate_residue_count__{ptm_group}'].clip(lower=1))

    enrichment_rows = []
    adjusted_rows = []
    role_rows = []
    density_rows = []
    for ptm_group in args.ptm_groups:
        target = proteins.loc[proteins[f'has_ptm__{ptm_group}'], 'is_condensate'].map({True: 'condensate', False: 'non_condensate'})
        universe = proteins['is_condensate'].map({True: 'condensate', False: 'non_condensate'})
        enrich = fisher_like_enrichment(target, universe)
        enrich['ptm_group'] = ptm_group
        enrichment_rows.append(enrich)
        adjusted_rows.append(fit_adjusted_model(proteins, ptm_group))

        role_target = proteins.loc[proteins[f'has_ptm__{ptm_group}'], 'condensate_role'].fillna('non_condensate')
        role_universe = proteins['condensate_role'].fillna('non_condensate')
        role_enrich = fisher_like_enrichment(role_target, role_universe)
        role_enrich['ptm_group'] = ptm_group
        role_rows.append(role_enrich)

        density_rows.append({
            'ptm_group': ptm_group,
            'mean_density_condensate': float(proteins.loc[proteins['is_condensate'], f'site_density_per_candidate__{ptm_group}'].fillna(0.0).mean()),
            'mean_density_non_condensate': float(proteins.loc[~proteins['is_condensate'], f'site_density_per_candidate__{ptm_group}'].fillna(0.0).mean()),
            'median_density_condensate': float(proteins.loc[proteins['is_condensate'], f'site_density_per_candidate__{ptm_group}'].fillna(0.0).median()),
            'median_density_non_condensate': float(proteins.loc[~proteins['is_condensate'], f'site_density_per_candidate__{ptm_group}'].fillna(0.0).median()),
        })

    enrichment = pd.concat(enrichment_rows, ignore_index=True)
    role_enrichment = pd.concat(role_rows, ignore_index=True)
    adjusted = pd.DataFrame(adjusted_rows)
    adjusted = add_bh_q_values(adjusted, p_col='adjusted_p_value', q_col='adjusted_q_value_bh')

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(proteins, outdir / 'protein_level_ptm_condensate_table.tsv')
    save_table(enrichment, outdir / 'condensate_ptm_enrichment.tsv')
    save_table(adjusted, outdir / 'condensate_ptm_adjusted_logistic_models.tsv')
    save_table(role_enrichment, outdir / 'condensate_role_ptm_enrichment.tsv')
    save_table(pd.DataFrame(density_rows), outdir / 'condensate_ptm_density_summary.tsv')

    plot_adjusted_models(adjusted, outdir / 'condensate_ptm_adjusted_odds_ratios')

    summary = pd.DataFrame([{
        'proteins_total': len(proteins),
        'condensate_proteins': int(proteins['is_condensate'].sum()),
        'ptm_groups': ';'.join(args.ptm_groups),
        'experimental_only': bool(args.experimental_only),
        'min_confidence_score': '' if pd.isna(args.min_confidence_score) else float(args.min_confidence_score),
        'output': str(outdir),
    }])
    save_table(summary, outdir / 'summary.tsv')
    print(summary.to_dict(orient='records')[0])


if __name__ == '__main__':
    main()

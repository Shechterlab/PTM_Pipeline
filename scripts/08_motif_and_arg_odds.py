from __future__ import annotations

import argparse
import math
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import beta as beta_dist
from scipy.stats import fisher_exact

from common import (
    add_bh_q_values,
    canonical_site_table,
    format_p_value,
    log2_odds_ratio,
    parse_fasta,
    residue_background_from_sequences,
    save_figure,
    save_table,
    sequence_window,
)

try:
    import statsmodels.api as sm
except ImportError:  # pragma: no cover
    sm = None


HYDROPHOBIC = set('AILMFWVY')
CENTERED_KMER_OFFSETS = [(-1, 0), (0, 1), (-1, 0, 1), (0, 1, 2), (-1, 0, 1, 2), (-2, -1, 0), (-2, -1, 0, 1), (-2, -1, 0, 1, 2)]
FEATURE_COLUMNS = [
    'feature_rg_plus1',
    'feature_rgg_plus12',
    'feature_grg_center',
    'feature_rxr_pm2',
    'feature_proline_flank3',
    'feature_hydrophobic_flank3',
    'feature_acidic_within2',
    'feature_acidic_within4',
    'feature_low_complexity_site',
    'feature_disordered_site',
    'feature_gly_fraction_flank3',
    'feature_gly_fraction_flank7',
    'feature_arg_fraction_flank7',
    'feature_pro_fraction_flank7',
    'feature_unique_aa_flank7',
    'feature_unique_aa_flank15',
]


def residue_at(window: str, relative_position: int) -> str:
    center = len(window) // 2
    index = center + relative_position
    if index < 0 or index >= len(window):
        return '_'
    return window[index]


def flank_residues(window: str, offsets: list[int]) -> list[str]:
    return [residue_at(window, offset) for offset in offsets if residue_at(window, offset) != '_']


def family_flags(window: str) -> dict[str, bool]:
    neighborhood = flank_residues(window, [i for i in range(-5, 6) if i != 0])
    near3 = flank_residues(window, [i for i in range(-3, 4) if i != 0])
    return {
        'motif_RG': residue_at(window, 1) == 'G',
        'motif_RGG': residue_at(window, 1) == 'G' and residue_at(window, 2) == 'G',
        'motif_GRG': residue_at(window, -1) == 'G' and residue_at(window, 1) == 'G',
        'motif_GAR': residue_at(window, -2) == 'G' and residue_at(window, -1) == 'A',
        'motif_RXR': residue_at(window, -2) == 'R' or residue_at(window, 2) == 'R',
        'motif_CARM1_proline_rich': 'P' in neighborhood,
        'motif_CARM1_hydrophobic_proline': ('P' in neighborhood) and any(aa in HYDROPHOBIC for aa in near3),
        'motif_PRMT5_acidic_DR_like': residue_at(window, -1) in {'D', 'E'} or residue_at(window, 1) in {'D', 'E'},
        'motif_PRMT5_acidic_within2': any(residue_at(window, i) in {'D', 'E'} for i in [-2, -1, 1, 2]),
    }


def exclusive_motif_family(window: str) -> str:
    if residue_at(window, 1) == 'G' and (residue_at(window, 2) == 'G' or residue_at(window, -1) == 'G'):
        return 'RGG / GRG box-like'
    if residue_at(window, 1) == 'G':
        return 'RG-only'
    if residue_at(window, -2) == 'R' or residue_at(window, 2) == 'R':
        return 'RXR-like'
    if 'P' in flank_residues(window, [i for i in range(-4, 5) if i != 0]) and any(aa in HYDROPHOBIC for aa in flank_residues(window, [-3, -2, -1, 1, 2, 3])):
        return 'CARM1 P/hydrophobic'
    if 'P' in flank_residues(window, [i for i in range(-4, 5) if i != 0]):
        return 'Proline-rich'
    if any(residue_at(window, i) in {'D', 'E'} for i in [-1, 1]):
        return 'Acidic DR-like'
    if any(residue_at(window, i) in {'D', 'E'} for i in [-2, -1, 1, 2]):
        return 'Acidic flank'
    if sum(1 for aa in flank_residues(window, [-2, -1, 1, 2]) if aa == 'R') >= 2:
        return 'Arg-rich'
    return 'Other / unmatched'


def low_complexity_mask(sequence: str, window: int = 12, max_unique: int = 3) -> np.ndarray:
    seq = str(sequence).strip().upper()
    if not seq:
        return np.zeros(0, dtype=bool)
    mask = np.zeros(len(seq), dtype=bool)
    if len(seq) <= window:
        mask[:] = len(set(seq)) <= max_unique
        return mask
    for start in range(0, len(seq) - window + 1):
        chunk = seq[start:start + window]
        if len(set(chunk)) <= max_unique:
            mask[start:start + window] = True
    return mask


def normalize_disorder_intervals(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    rename_map = {}
    for source, target in [('accession', 'canonical_UniProtAC'), ('Entry', 'canonical_UniProtAC'), ('start', 'fragment_start'), ('end', 'fragment_end')]:
        if source in out.columns and target not in out.columns:
            rename_map[source] = target
    out = out.rename(columns=rename_map)
    required = ['canonical_UniProtAC', 'fragment_start', 'fragment_end']
    missing = [col for col in required if col not in out.columns]
    if missing:
        raise ValueError(f'missing disorder interval columns: {missing}')
    out['fragment_start'] = pd.to_numeric(out['fragment_start'], errors='coerce')
    out['fragment_end'] = pd.to_numeric(out['fragment_end'], errors='coerce')
    out = out[out['fragment_start'].notna() & out['fragment_end'].notna()].copy()
    out['fragment_start'] = out['fragment_start'].astype(int)
    out['fragment_end'] = out['fragment_end'].astype(int)
    return out


def odds_ratio(a: int, b: int, c: int, d: int) -> float:
    return ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))


def motif_enrichment(target: pd.DataFrame, background: pd.DataFrame, motif_columns: list[str], background_label: str) -> pd.DataFrame:
    rows = []
    for column in motif_columns:
        a = int(target[column].sum())
        b = int(len(target) - a)
        c = int(background[column].sum())
        d = int(len(background) - c)
        rows.append({
            'motif_family': column.replace('motif_', ''),
            'target_count': a,
            'target_fraction': a / len(target) if len(target) else np.nan,
            'background_count': c,
            'background_fraction': c / len(background) if len(background) else np.nan,
            'odds_ratio': odds_ratio(a, b, c, d),
            'log2_odds_ratio': log2_odds_ratio(odds_ratio(a, b, c, d)),
            'p_value': float(fisher_exact([[a, b], [c, d]], alternative='two-sided').pvalue),
            'background_model': background_label,
        })
    out = add_bh_q_values(pd.DataFrame(rows), p_col='p_value', q_col='q_value_bh')
    return out.sort_values(['q_value_bh', 'odds_ratio', 'target_count'], ascending=[True, False, False])


def categorical_enrichment(target: pd.Series, background: pd.Series, label_col: str = 'category', background_label: str = '') -> pd.DataFrame:
    target = target.fillna('NA').astype(str)
    background = background.fillna('NA').astype(str)
    rows = []
    for label in sorted(set(target) | set(background)):
        a = int((target == label).sum())
        b = int(len(target) - a)
        c = int((background == label).sum())
        d = int(len(background) - c)
        rows.append({
            label_col: label,
            'target_count': a,
            'target_fraction': a / len(target) if len(target) else np.nan,
            'background_count': c,
            'background_fraction': c / len(background) if len(background) else np.nan,
            'odds_ratio': odds_ratio(a, b, c, d),
            'log2_odds_ratio': log2_odds_ratio(odds_ratio(a, b, c, d)),
            'p_value': float(fisher_exact([[a, b], [c, d]], alternative='two-sided').pvalue),
            'background_model': background_label,
        })
    out = add_bh_q_values(pd.DataFrame(rows), p_col='p_value', q_col='q_value_bh')
    return out.sort_values(['q_value_bh', 'odds_ratio', 'target_count'], ascending=[True, False, False])


def positional_enrichment(target_windows: pd.Series, background_windows: pd.Series) -> pd.DataFrame:
    rows = []
    if target_windows.empty or background_windows.empty:
        return pd.DataFrame()
    flank = len(str(target_windows.iloc[0])) // 2
    amino_acids = sorted(set(''.join(target_windows.astype(str).tolist() + background_windows.astype(str).tolist())) - {'_'})
    for rel in range(-flank, flank + 1):
        if rel == 0:
            continue
        target_chars = target_windows.astype(str).map(lambda x: residue_at(x, rel))
        background_chars = background_windows.astype(str).map(lambda x: residue_at(x, rel))
        for aa in amino_acids:
            a = int((target_chars == aa).sum())
            b = int(len(target_chars) - a)
            c = int((background_chars == aa).sum())
            d = int(len(background_chars) - c)
            rows.append({
                'relative_position': rel,
                'amino_acid': aa,
                'target_count': a,
                'background_count': c,
                'odds_ratio': odds_ratio(a, b, c, d),
                'log2_odds_ratio': math.log2(odds_ratio(a, b, c, d)),
                'p_value': float(fisher_exact([[a, b], [c, d]], alternative='two-sided').pvalue),
            })
    out = add_bh_q_values(pd.DataFrame(rows), p_col='p_value', q_col='q_value_bh')
    return out


def top_protein_label(row: pd.Series) -> str:
    gene = str(row.get('substrate_genename', '') or '').strip()
    accession = str(row.get('canonical_UniProtAC', '') or '').strip()
    if gene:
        return f'{gene} ({accession})'
    return accession


def estimate_beta_prior(successes: np.ndarray, totals: np.ndarray) -> tuple[float, float, float]:
    successes = np.asarray(successes, dtype=float)
    totals = np.asarray(totals, dtype=float)
    mask = totals > 0
    successes = successes[mask]
    totals = totals[mask]
    if len(totals) == 0:
        return 1.0, 1.0, 0.5
    mu = float(successes.sum() / totals.sum())
    mu = min(max(mu, 1e-8), 1 - 1e-8)
    prop = successes / totals
    observed_var = float(np.var(prop, ddof=1)) if len(prop) > 1 else 0.0
    sampling_var = float(mu * (1 - mu) * np.mean(1 / totals))
    latent_var = max(observed_var - sampling_var, 1e-8)
    prior_strength = min(max((mu * (1 - mu) / latent_var) - 1.0, 2.0), 500.0)
    alpha = mu * prior_strength
    beta = (1 - mu) * prior_strength
    return float(alpha), float(beta), mu


def build_disorder_masks(sequences: dict[str, str], disorder_intervals: pd.DataFrame | None) -> dict[str, np.ndarray]:
    masks: dict[str, np.ndarray] = {}
    if disorder_intervals is None or disorder_intervals.empty:
        return masks
    for accession, group in disorder_intervals.groupby('canonical_UniProtAC'):
        seq_len = len(sequences.get(accession, ''))
        if seq_len == 0:
            continue
        mask = np.zeros(seq_len, dtype=bool)
        for row in group.itertuples(index=False):
            start = max(1, int(getattr(row, 'fragment_start'))) - 1
            end = min(seq_len, int(getattr(row, 'fragment_end')))
            mask[start:end] = True
        masks[accession] = mask
    return masks


def protein_covariates(sequences: dict[str, str], disorder_intervals: pd.DataFrame | None) -> pd.DataFrame:
    rows = []
    disorder_masks = build_disorder_masks(sequences, disorder_intervals)
    for accession, sequence in sequences.items():
        seq = str(sequence)
        if not seq:
            continue
        lc_mask = low_complexity_mask(seq)
        dis_mask = disorder_masks.get(accession, np.zeros(len(seq), dtype=bool))
        rows.append({
            'canonical_UniProtAC': accession,
            'protein_length': len(seq),
            'low_complexity_fraction': float(lc_mask.mean()) if len(lc_mask) else 0.0,
            'disorder_fraction': float(dis_mask.mean()) if len(dis_mask) else np.nan,
        })
    return pd.DataFrame(rows)


def fit_expected_rate_model(protein_counts: pd.DataFrame) -> tuple[pd.Series, str]:
    if sm is None:
        return pd.Series(protein_counts['global_methyl_rate'].to_numpy(), index=protein_counts.index), 'global_rate_fallback'
    predictors = ['low_complexity_fraction']
    if 'disorder_fraction' in protein_counts.columns and protein_counts['disorder_fraction'].notna().any():
        predictors.append('disorder_fraction')
    design = protein_counts[predictors].copy()
    if design.empty:
        return pd.Series(protein_counts['global_methyl_rate'].to_numpy(), index=protein_counts.index), 'global_rate_fallback'
    for col in design.columns:
        if design[col].notna().any():
            design[col] = design[col].fillna(float(design[col].median()))
        else:
            design[col] = 0.0
    design = sm.add_constant(design, has_constant='add')
    try:
        model = sm.GLM(
            protein_counts['methyl_site_count'] / protein_counts['arginine_count'],
            design,
            family=sm.families.Binomial(),
            var_weights=protein_counts['arginine_count'],
        )
        fit = model.fit(maxiter=200, disp=False)
        pred = pd.Series(np.clip(fit.predict(design), 1e-8, 1 - 1e-8), index=protein_counts.index)
        return pred, 'glm_covariate_adjusted'
    except Exception:
        return pd.Series(protein_counts['global_methyl_rate'].to_numpy(), index=protein_counts.index), 'global_rate_fallback'


def compute_site_features(df: pd.DataFrame, sequences: dict[str, str], disorder_masks: dict[str, np.ndarray]) -> pd.DataFrame:
    out = df.copy()
    unique_aa_7 = []
    unique_aa_15 = []
    low_complexity = []
    disordered = []
    pos_col = 'position' if 'position' in out.columns else 'corrected_position'
    for row in out.itertuples(index=False):
        accession = str(getattr(row, 'canonical_UniProtAC'))
        position = int(getattr(row, pos_col))
        window7 = str(getattr(row, 'window_7'))
        window15 = str(getattr(row, 'window_15'))
        sequence = sequences.get(accession, '')
        lc_mask = low_complexity_mask(sequence) if sequence else np.zeros(0, dtype=bool)
        dis_mask = disorder_masks.get(accession, np.zeros(len(sequence), dtype=bool)) if sequence else np.zeros(0, dtype=bool)
        index0 = position - 1
        low_complexity.append(bool(index0 >= 0 and index0 < len(lc_mask) and lc_mask[index0]))
        disordered.append(bool(index0 >= 0 and index0 < len(dis_mask) and dis_mask[index0]))
        flank7 = [aa for aa in window7 if aa != '_']
        flank15 = [aa for aa in window15 if aa != '_']
        unique_aa_7.append(len(set(flank7)))
        unique_aa_15.append(len(set(flank15)))
    out['site_low_complexity'] = low_complexity
    out['site_disordered'] = disordered
    out['feature_low_complexity_site'] = out['site_low_complexity'].astype(int)
    out['feature_disordered_site'] = out['site_disordered'].astype(int)
    out['feature_rg_plus1'] = out['window_7'].astype(str).map(lambda w: int(residue_at(w, 1) == 'G'))
    out['feature_rgg_plus12'] = out['window_7'].astype(str).map(lambda w: int(residue_at(w, 1) == 'G' and residue_at(w, 2) == 'G'))
    out['feature_grg_center'] = out['window_7'].astype(str).map(lambda w: int(residue_at(w, -1) == 'G' and residue_at(w, 1) == 'G'))
    out['feature_rxr_pm2'] = out['window_7'].astype(str).map(lambda w: int(residue_at(w, -2) == 'R' or residue_at(w, 2) == 'R'))
    out['feature_proline_flank3'] = out['window_7'].astype(str).map(lambda w: int('P' in flank_residues(w, [-3, -2, -1, 1, 2, 3])))
    out['feature_hydrophobic_flank3'] = out['window_7'].astype(str).map(lambda w: int(any(aa in HYDROPHOBIC for aa in flank_residues(w, [-3, -2, -1, 1, 2, 3]))))
    out['feature_acidic_within2'] = out['window_7'].astype(str).map(lambda w: int(any(residue_at(w, i) in {'D', 'E'} for i in [-2, -1, 1, 2])))
    out['feature_acidic_within4'] = out['window_15'].astype(str).map(lambda w: int(any(residue_at(w, i) in {'D', 'E'} for i in [-4, -3, -2, -1, 1, 2, 3, 4])))
    out['feature_gly_fraction_flank3'] = out['window_7'].astype(str).map(lambda w: float(np.mean([aa == 'G' for aa in flank_residues(w, [-3, -2, -1, 1, 2, 3])])) if flank_residues(w, [-3, -2, -1, 1, 2, 3]) else 0.0)
    out['feature_gly_fraction_flank7'] = out['window_15'].astype(str).map(lambda w: float(np.mean([aa == 'G' for aa in flank_residues(w, [i for i in range(-7, 8) if i != 0])])) if flank_residues(w, [i for i in range(-7, 8) if i != 0]) else 0.0)
    out['feature_arg_fraction_flank7'] = out['window_15'].astype(str).map(lambda w: float(np.mean([aa == 'R' for aa in flank_residues(w, [i for i in range(-7, 8) if i != 0])])) if flank_residues(w, [i for i in range(-7, 8) if i != 0]) else 0.0)
    out['feature_pro_fraction_flank7'] = out['window_15'].astype(str).map(lambda w: float(np.mean([aa == 'P' for aa in flank_residues(w, [i for i in range(-7, 8) if i != 0])])) if flank_residues(w, [i for i in range(-7, 8) if i != 0]) else 0.0)
    out['feature_unique_aa_flank7'] = unique_aa_7
    out['feature_unique_aa_flank15'] = unique_aa_15
    out['exclusive_motif_family'] = out['window_7'].astype(str).map(exclusive_motif_family)
    return out


def choose_matched_controls(sites: pd.DataFrame, background: pd.DataFrame) -> pd.DataFrame:
    if sites.empty or background.empty:
        return background.head(0).copy()
    chosen_parts = []
    site_pos_col = 'position' if 'position' in sites.columns else 'corrected_position'
    bg_pos_col = 'position' if 'position' in background.columns else 'corrected_position'
    for accession, site_group in sites.groupby('canonical_UniProtAC'):
        bg_group = background[background['canonical_UniProtAC'] == accession].copy()
        if bg_group.empty:
            continue
        used = set()
        selected = []
        strata_cols = ['site_low_complexity', 'site_disordered']
        for site_row in site_group.sort_values(site_pos_col).itertuples(index=False):
            candidates = bg_group.copy()
            for col in strata_cols:
                value = getattr(site_row, col, None)
                if col in candidates.columns and candidates[col].nunique() > 1:
                    subset = candidates[candidates[col] == value].copy()
                    if not subset.empty:
                        candidates = subset
            candidates = candidates[~candidates['site_key'].isin(used)].copy()
            if candidates.empty:
                candidates = bg_group[~bg_group['site_key'].isin(used)].copy()
            if candidates.empty:
                candidates = bg_group.copy()
            candidates['distance'] = (candidates[bg_pos_col] - int(getattr(site_row, site_pos_col))).abs()
            pick = candidates.sort_values(['distance', bg_pos_col, 'site_key']).iloc[0]
            selected.append(pick.drop(labels=['distance']).to_dict())
            used.add(str(pick['site_key']))
        if selected:
            chosen_parts.append(pd.DataFrame(selected))
    if not chosen_parts:
        return background.head(0).copy()
    matched = pd.concat(chosen_parts, ignore_index=True)
    matched['matched_background'] = True
    return matched


def centered_kmer_value(window: str, offsets: tuple[int, ...]) -> str:
    letters = [residue_at(window, offset) for offset in offsets]
    if '_' in letters:
        return ''
    return ''.join(letters)


def centered_kmer_enrichment(target: pd.DataFrame, background: pd.DataFrame, min_target_count: int = 25) -> pd.DataFrame:
    rows = []
    if target.empty or background.empty:
        return pd.DataFrame()
    background_windows = background['window_15'].astype(str)
    for offsets in CENTERED_KMER_OFFSETS:
        label = 'motif_' + '_'.join(f'{offset:+d}' for offset in offsets)
        target_vals = target['window_15'].astype(str).map(lambda w: centered_kmer_value(w, offsets))
        background_vals = background_windows.map(lambda w: centered_kmer_value(w, offsets))
        all_kmers = sorted(set(target_vals) | set(background_vals))
        for kmer in all_kmers:
            if not kmer:
                continue
            a = int((target_vals == kmer).sum())
            if a < min_target_count:
                continue
            b = int(len(target_vals) - a)
            c = int((background_vals == kmer).sum())
            d = int(len(background_vals) - c)
            rows.append({
                'kmer_definition': label,
                'kmer': kmer,
                'target_count': a,
                'background_count': c,
                'target_fraction': a / len(target_vals),
                'background_fraction': c / len(background_vals),
                'odds_ratio': odds_ratio(a, b, c, d),
                'log2_odds_ratio': log2_odds_ratio(odds_ratio(a, b, c, d)),
                'p_value': float(fisher_exact([[a, b], [c, d]], alternative='two-sided').pvalue),
            })
    if not rows:
        return pd.DataFrame()
    out = add_bh_q_values(pd.DataFrame(rows), p_col='p_value', q_col='q_value_bh')
    return out.sort_values(['q_value_bh', 'odds_ratio', 'target_count'], ascending=[True, False, False])


def fit_site_feature_model(target: pd.DataFrame, background: pd.DataFrame) -> pd.DataFrame:
    if target.empty or background.empty:
        return pd.DataFrame(columns=['feature', 'coefficient', 'p_value', 'q_value_bh', 'model'])
    combined = pd.concat([
        target.assign(is_methyl=1),
        background.assign(is_methyl=0),
    ], ignore_index=True)
    rows = []
    if sm is not None:
        design = combined[FEATURE_COLUMNS].copy()
        for col in design.columns:
            design[col] = pd.to_numeric(design[col], errors='coerce').fillna(0.0)
        design = sm.add_constant(design, has_constant='add')
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                fit = sm.GLM(combined['is_methyl'], design, family=sm.families.Binomial()).fit(maxiter=200, disp=False)
            params = fit.params.drop('const', errors='ignore')
            pvals = fit.pvalues.drop('const', errors='ignore')
            for feature, coef in params.items():
                rows.append({'feature': feature, 'coefficient': float(coef), 'p_value': float(pvals.get(feature, np.nan)), 'model': 'multivariable_glm'})
            out = add_bh_q_values(pd.DataFrame(rows), p_col='p_value', q_col='q_value_bh')
            return out.sort_values(['q_value_bh', 'coefficient'], ascending=[True, False])
        except Exception:
            rows = []
    for feature in FEATURE_COLUMNS:
        t = pd.to_numeric(target[feature], errors='coerce').fillna(0.0)
        b = pd.to_numeric(background[feature], errors='coerce').fillna(0.0)
        if set(np.unique(pd.concat([t, b]))) <= {0.0, 1.0}:
            a = int(t.sum())
            b_not = int(len(t) - a)
            c = int(b.sum())
            d = int(len(b) - c)
            coef = math.log(odds_ratio(a, b_not, c, d))
            pval = float(fisher_exact([[a, b_not], [c, d]], alternative='two-sided').pvalue)
        else:
            coef = float(t.mean() - b.mean())
            pval = np.nan
        rows.append({'feature': feature, 'coefficient': coef, 'p_value': pval, 'model': 'univariate_fallback'})
    out = add_bh_q_values(pd.DataFrame(rows), p_col='p_value', q_col='q_value_bh')
    return out.sort_values(['q_value_bh', 'coefficient'], ascending=[True, False])


def plot_barh_with_stats(plot_df: pd.DataFrame, y_col: str, x_col: str, outpath: Path, xlabel: str, title: str, stat_col: str, count_col: str, max_items: int = 12, min_count: int = 20) -> None:
    if plot_df.empty:
        return
    use = plot_df[plot_df[count_col] >= min_count].head(max_items).sort_values(x_col)
    if use.empty:
        return
    fig, ax = plt.subplots(figsize=(10.0, max(5.2, 0.42 * len(use) + 1.8)))
    ax.barh(use[y_col], use[x_col])
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Category')
    ax.set_title(title)
    for y, x, n, q in zip(use[y_col], use[x_col], use[count_col], use[stat_col]):
        ax.text(x, y, f'  n={n}, q={format_p_value(q)}', va='center', ha='left' if x >= 0 else 'right', fontsize=8.5)
    fig.tight_layout()
    save_figure(fig, outpath)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--integrated-sites', required=True)
    ap.add_argument('--canonical-fasta', required=True)
    ap.add_argument('--disorder-intervals', default='')
    ap.add_argument('--outdir', default='results/motif_arg_odds')
    ap.add_argument('--window-flank', type=int, default=7)
    ap.add_argument('--window-flank-wide', type=int, default=15)
    ap.add_argument('--min-arg-count', type=int, default=10)
    ap.add_argument('--min-site-count', type=int, default=2)
    ap.add_argument('--min-kmer-count', type=int, default=25)
    args = ap.parse_args()

    sites = pd.read_csv(args.integrated_sites, sep='\t', low_memory=False)
    sites = canonical_site_table(sites, accession_col='canonical_UniProtAC', position_col='corrected_position', residue_col='residue')
    sequences = parse_fasta(args.canonical_fasta)
    sites = sites[sites['canonical_UniProtAC'].isin(sequences)].copy()
    sites['window_15'] = sites.apply(lambda row: sequence_window(sequences[row['canonical_UniProtAC']], int(row['corrected_position']), args.window_flank), axis=1)
    sites['window_31'] = sites.apply(lambda row: sequence_window(sequences[row['canonical_UniProtAC']], int(row['corrected_position']), args.window_flank_wide), axis=1)
    sites = sites.rename(columns={'window_15': 'window_15', 'window_31': 'window_31'})
    # Rename for fixed feature helper expectations.
    sites['window_7'] = sites['window_15']
    sites['window_15'] = sites['window_31']

    proteome_arg = residue_background_from_sequences(sequences, residue='R', flank=args.window_flank_wide)
    proteome_arg = proteome_arg.rename(columns={f'window_{2 * args.window_flank_wide + 1}': 'window_15'})
    proteome_arg['window_7'] = proteome_arg['window_15'].astype(str).map(lambda w: w[args.window_flank_wide - args.window_flank: args.window_flank_wide + args.window_flank + 1])

    methylated_proteins = sorted(set(sites['canonical_UniProtAC'].dropna().astype(str)))
    methyl_protein_arg = proteome_arg[proteome_arg['canonical_UniProtAC'].isin(methylated_proteins)].copy()
    methyl_site_keys = set(sites['site_key'].astype(str))
    background_nonmethyl = methyl_protein_arg[~methyl_protein_arg['site_key'].isin(methyl_site_keys)].copy()

    disorder_intervals = normalize_disorder_intervals(pd.read_csv(args.disorder_intervals, sep='\t', low_memory=False)) if args.disorder_intervals else None
    disorder_masks = build_disorder_masks(sequences, disorder_intervals)

    sites = compute_site_features(sites, sequences, disorder_masks)
    background_nonmethyl = compute_site_features(background_nonmethyl, sequences, disorder_masks)
    matched_background = choose_matched_controls(sites, background_nonmethyl)

    motif_columns = sorted([c for c in sites.columns if c.startswith('motif_')])
    if not motif_columns:
        motif_frame = pd.DataFrame([family_flags(window) for window in sites['window_7']])
        sites = pd.concat([sites.reset_index(drop=True), motif_frame], axis=1)
        bg_frame = pd.DataFrame([family_flags(window) for window in background_nonmethyl['window_7']])
        background_nonmethyl = pd.concat([background_nonmethyl.reset_index(drop=True), bg_frame], axis=1)
        matched_frame = pd.DataFrame([family_flags(window) for window in matched_background['window_7']]) if not matched_background.empty else pd.DataFrame(columns=sites.columns)
        matched_background = pd.concat([matched_background.reset_index(drop=True), matched_frame], axis=1) if not matched_background.empty else matched_background
        motif_columns = sorted([c for c in sites.columns if c.startswith('motif_')])

    motif_vs_methyl_proteins = motif_enrichment(sites, background_nonmethyl, motif_columns, 'arginines_in_methylated_proteins')
    motif_vs_matched = motif_enrichment(sites, matched_background, motif_columns, 'protein_matched_nonmethyl_arginines') if not matched_background.empty else pd.DataFrame()
    exclusive_vs_matched = categorical_enrichment(sites['exclusive_motif_family'], matched_background['exclusive_motif_family'], label_col='exclusive_motif_family', background_label='protein_matched_nonmethyl_arginines') if not matched_background.empty else pd.DataFrame()
    positional = positional_enrichment(sites['window_15'], matched_background['window_15']) if not matched_background.empty else pd.DataFrame()
    centered_kmers = centered_kmer_enrichment(sites, matched_background, min_target_count=args.min_kmer_count) if not matched_background.empty else pd.DataFrame()
    feature_model = fit_site_feature_model(sites, matched_background) if not matched_background.empty else pd.DataFrame()

    arg_counts = residue_background_from_sequences(sequences, residue='R')[['canonical_UniProtAC', 'site_key']].groupby('canonical_UniProtAC').size().rename('arginine_count').reset_index()
    protein_counts = sites.groupby(['canonical_UniProtAC', 'substrate_genename'], as_index=False).size().rename(columns={'size': 'methyl_site_count'})
    protein_counts = protein_counts.merge(arg_counts, how='left', on='canonical_UniProtAC')
    covariates = protein_covariates(sequences, disorder_intervals)
    protein_counts = protein_counts.merge(covariates, how='left', on='canonical_UniProtAC')
    protein_counts['global_methyl_rate'] = float(protein_counts['methyl_site_count'].sum() / protein_counts['arginine_count'].sum())
    alpha0, beta0, global_rate = estimate_beta_prior(protein_counts['methyl_site_count'].to_numpy(), protein_counts['arginine_count'].to_numpy())
    protein_counts['posterior_alpha'] = protein_counts['methyl_site_count'] + alpha0
    protein_counts['posterior_beta'] = protein_counts['arginine_count'] - protein_counts['methyl_site_count'] + beta0
    protein_counts['posterior_mean_rate'] = protein_counts['posterior_alpha'] / (protein_counts['posterior_alpha'] + protein_counts['posterior_beta'])
    protein_counts['posterior_ci_low'] = beta_dist.ppf(0.025, protein_counts['posterior_alpha'], protein_counts['posterior_beta'])
    protein_counts['posterior_ci_high'] = beta_dist.ppf(0.975, protein_counts['posterior_alpha'], protein_counts['posterior_beta'])
    protein_counts['posterior_prob_above_global'] = 1.0 - beta_dist.cdf(global_rate, protein_counts['posterior_alpha'], protein_counts['posterior_beta'])
    protein_counts['expected_rate'], expected_model = fit_expected_rate_model(protein_counts)
    protein_counts['expected_rate_model'] = expected_model
    protein_counts['shrinkage_adjusted_log2_enrichment'] = np.log2(protein_counts['posterior_mean_rate'] / protein_counts['expected_rate'])
    protein_counts['posterior_log2_enrichment_vs_global'] = np.log2(protein_counts['posterior_mean_rate'] / global_rate)
    protein_counts['shrinkage_rank_score'] = protein_counts['shrinkage_adjusted_log2_enrichment'] * np.sqrt(protein_counts['arginine_count']) * (protein_counts['posterior_prob_above_global'] - 0.5)
    protein_counts['label'] = protein_counts.apply(top_protein_label, axis=1)
    protein_counts = protein_counts[(protein_counts['arginine_count'] >= args.min_arg_count) & (protein_counts['methyl_site_count'] >= args.min_site_count)].copy()
    protein_counts = protein_counts.sort_values(['shrinkage_rank_score', 'posterior_prob_above_global', 'methyl_site_count'], ascending=[False, False, False])

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(sites, outdir / 'sequence_windows.tsv')
    save_table(background_nonmethyl, outdir / 'background_nonmethyl_arginines_same_proteins.tsv')
    save_table(matched_background, outdir / 'matched_nonmethyl_arginines_same_proteins.tsv')
    save_table(motif_vs_methyl_proteins, outdir / 'motif_family_feature_enrichment_vs_methylated_protein_arginines.tsv')
    save_table(motif_vs_matched, outdir / 'motif_family_feature_enrichment_vs_matched_nonmethyl_arginines.tsv')
    save_table(exclusive_vs_matched, outdir / 'exclusive_motif_family_enrichment_vs_matched_nonmethyl_arginines.tsv')
    save_table(positional, outdir / 'matched_control_positional_amino_acid_enrichment.tsv')
    save_table(centered_kmers, outdir / 'de_novo_centered_kmer_enrichment.tsv')
    save_table(feature_model, outdir / 'motif_feature_model.tsv')
    save_table(protein_counts, outdir / 'per_protein_shrinkage_model.tsv')

    plot_motif = exclusive_vs_matched[exclusive_vs_matched['exclusive_motif_family'] != 'Other / unmatched'].copy() if not exclusive_vs_matched.empty else pd.DataFrame()
    plot_motif = plot_motif.sort_values(['q_value_bh', 'odds_ratio', 'target_count'], ascending=[True, False, False]).head(10)
    plot_barh_with_stats(
        plot_motif,
        y_col='exclusive_motif_family',
        x_col='log2_odds_ratio',
        outpath=outdir / 'arg_methyl_exclusive_motif_family_enrichment',
        xlabel='log2(OR) vs protein-matched non-methyl Arg controls',
        title='Exclusive methylarginine motif-family enrichment',
        stat_col='q_value_bh',
        count_col='target_count',
        max_items=10,
        min_count=15,
    )

    if not centered_kmers.empty:
        plot_kmers = centered_kmers.head(12).sort_values('log2_odds_ratio')
        figk, axk = plt.subplots(figsize=(10.2, max(5.5, 0.42 * len(plot_kmers) + 1.7)))
        labels = plot_kmers['kmer'] + ' [' + plot_kmers['kmer_definition'].str.replace('motif_', '', regex=False) + ']'
        axk.barh(labels, plot_kmers['log2_odds_ratio'])
        axk.axvline(0, color='black', linewidth=0.8)
        axk.set_xlabel('log2(OR) vs protein-matched non-methyl Arg controls')
        axk.set_ylabel('Centered k-mer')
        axk.set_title('De novo enriched centered motifs around methylarginines')
        for y, x, n, q in zip(labels, plot_kmers['log2_odds_ratio'], plot_kmers['target_count'], plot_kmers['q_value_bh']):
            axk.text(x, y, f'  n={n}, q={format_p_value(q)}', va='center', ha='left' if x >= 0 else 'right', fontsize=8.3)
        figk.tight_layout()
        save_figure(figk, outdir / 'arg_methyl_de_novo_centered_kmers')

    heatmap_df = positional.copy()
    selected_aa = (
        heatmap_df.groupby('amino_acid')['target_count'].sum().sort_values(ascending=False).head(12).index.tolist()
    ) if not heatmap_df.empty else []
    if selected_aa:
        matrix = heatmap_df[heatmap_df['amino_acid'].isin(selected_aa)].pivot(index='amino_acid', columns='relative_position', values='log2_odds_ratio').fillna(0.0)
        q_matrix = heatmap_df[heatmap_df['amino_acid'].isin(selected_aa)].pivot(index='amino_acid', columns='relative_position', values='q_value_bh')
        q_matrix = q_matrix.reindex(index=matrix.index, columns=matrix.columns)
        p_floor = np.nextafter(0.0, 1.0)
        neglog10_q = -np.log10(q_matrix.fillna(1.0).clip(lower=p_floor))
        fig2, (ax2, ax2b) = plt.subplots(ncols=2, figsize=(15.0, 5.8), gridspec_kw={'width_ratios': [1.0, 1.0]})
        im = ax2.imshow(matrix.to_numpy(), aspect='auto', cmap='coolwarm', vmin=-3, vmax=3)
        ax2.set_yticks(range(len(matrix.index)))
        ax2.set_yticklabels(matrix.index)
        ax2.set_xticks(range(len(matrix.columns)))
        ax2.set_xticklabels(matrix.columns)
        ax2.set_xlabel('Position relative to methyl-Arg')
        ax2.set_title('Matched-control positional amino-acid enrichment')
        for row_idx, aa in enumerate(matrix.index):
            for col_idx, rel in enumerate(matrix.columns):
                qval = q_matrix.loc[aa, rel]
                if pd.isna(qval):
                    continue
                mark = '***' if qval < 1e-10 else '**' if qval < 1e-3 else '*' if qval < 0.05 else ''
                if mark:
                    ax2.text(col_idx, row_idx, mark, ha='center', va='center', fontsize=7, color='black')
        fig2.colorbar(im, ax=ax2, label='log2(OR)')
        im_p = ax2b.imshow(neglog10_q.to_numpy(), aspect='auto', cmap='viridis')
        ax2b.set_yticks(range(len(neglog10_q.index)))
        ax2b.set_yticklabels(neglog10_q.index)
        ax2b.set_xticks(range(len(neglog10_q.columns)))
        ax2b.set_xticklabels(neglog10_q.columns)
        ax2b.set_xlabel('Position relative to methyl-Arg')
        ax2b.set_title('Matched-control positional BH-corrected significance')
        fig2.colorbar(im_p, ax=ax2b, label='-log10(q)')
        fig2.tight_layout()
        save_figure(fig2, outdir / 'arg_methyl_positional_enrichment_heatmap')

    if not feature_model.empty:
        feat_plot = feature_model.head(12).sort_values('coefficient')
        figf, axf = plt.subplots(figsize=(9.8, max(5.2, 0.42 * len(feat_plot) + 1.8)))
        axf.barh(feat_plot['feature'].str.replace('feature_', '', regex=False), feat_plot['coefficient'])
        axf.axvline(0, color='black', linewidth=0.8)
        axf.set_xlabel('Model coefficient (log-odds scale when multivariable)')
        axf.set_ylabel('Feature')
        axf.set_title('Discriminative motif/context features for methylarginines')
        for y, x, q in zip(feat_plot['feature'].str.replace('feature_', '', regex=False), feat_plot['coefficient'], feat_plot['q_value_bh']):
            axf.text(x, y, f'  q={format_p_value(q)}', va='center', ha='left' if x >= 0 else 'right', fontsize=8.3)
        figf.tight_layout()
        save_figure(figf, outdir / 'arg_methyl_feature_model_coefficients')

    top_score = protein_counts.head(20).sort_values('shrinkage_adjusted_log2_enrichment')
    fig3, ax3 = plt.subplots(figsize=(10.6, 7.3))
    ax3.barh(top_score['label'], top_score['shrinkage_adjusted_log2_enrichment'])
    ax3.axvline(0, color='black', linewidth=0.8)
    ax3.set_xlabel('Shrinkage-adjusted log2 enrichment vs expected methylation rate')
    ax3.set_ylabel('Protein')
    ax3.set_title('Top proteins by shrinkage-adjusted methylarginine enrichment')
    for y, v, n, prob in zip(top_score['label'], top_score['shrinkage_adjusted_log2_enrichment'], top_score['methyl_site_count'], top_score['posterior_prob_above_global']):
        ax3.text(v, y, f'  n={n}, postP={prob:.3f}', va='center', ha='left' if v >= 0 else 'right', fontsize=8)
    fig3.tight_layout()
    save_figure(fig3, outdir / 'top_proteins_by_shrinkage_adjusted_enrichment')

    fig4, ax4 = plt.subplots(figsize=(7.4, 5.4))
    ax4.scatter(protein_counts['methyl_site_count'], protein_counts['shrinkage_adjusted_log2_enrichment'], alpha=0.55)
    ax4.set_xlabel('Methylarginine site count')
    ax4.set_ylabel('Shrinkage-adjusted log2 enrichment')
    ax4.set_title('Protein-level methylarginine burden and shrinkage-adjusted enrichment')
    fig4.tight_layout()
    save_figure(fig4, outdir / 'protein_site_count_vs_shrinkage_adjusted_enrichment')

    print({
        'methyl_sites': len(sites),
        'background_arginines_same_proteins': len(background_nonmethyl),
        'matched_controls': len(matched_background),
        'proteins_with_shrinkage_model': len(protein_counts),
        'beta_prior_alpha': round(alpha0, 4),
        'beta_prior_beta': round(beta0, 4),
        'expected_rate_model': expected_model,
        'output': str(outdir),
    })


if __name__ == '__main__':
    main()

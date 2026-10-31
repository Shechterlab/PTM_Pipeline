from __future__ import annotations

import json
import math
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import chi2, fisher_exact


BREWER_COLORS = {
    'blue': '#377eb8',
    'orange': '#d95f02',
    'green': '#1b9e77',
    'purple': '#7570b3',
    'magenta': '#e7298a',
    'brown': '#a6761d',
    'red': '#e41a1c',
    'teal': '#66c2a5',
    'salmon': '#fc8d62',
    'lavender': '#8da0cb',
    'light_green': '#a6d854',
    'gold': '#ffd92f',
    'gray': '#6b6b6b',
    'light_gray': '#d9d9d9',
    'mid_gray': '#bdbdbd',
    'dark_gray': '#4d4d4d',
}

_ARIAL_REGISTERED = False


def save_table(df: pd.DataFrame, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() == '.tsv':
        df.to_csv(path, sep='\t', index=False)
    else:
        df.to_csv(path, index=False)


def save_figure(fig, path_prefix: str | Path, dpi: int = 300) -> None:
    path_prefix = Path(path_prefix)
    path_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path_prefix.with_suffix('.png'), dpi=dpi, bbox_inches='tight', facecolor='white')
    if str(os.environ.get('PTM_SKIP_PDF', '')).lower() in {'1', 'true', 'yes'}:
        return
    try:
        fig.savefig(path_prefix.with_suffix('.pdf'), bbox_inches='tight', facecolor='white')
    except PermissionError:
        # Some mounted-drive files can remain externally locked; keep the PNG/TSV outputs fresh.
        pass


def apply_paper_style() -> None:
    import matplotlib as mpl
    from matplotlib import font_manager

    global _ARIAL_REGISTERED
    if not _ARIAL_REGISTERED:
        font_dirs = [path for path in os.environ.get('PTM_FONT_DIRS', '').split(os.pathsep) if path]
        for font_dir in font_dirs:
            for font_path in Path(font_dir).glob('*.ttf'):
                try:
                    font_manager.fontManager.addfont(str(font_path))
                except RuntimeError:
                    pass
        _ARIAL_REGISTERED = True

    preferred_sans = ['Liberation Sans', 'DejaVu Sans', 'Arial', 'Helvetica', 'sans-serif']
    if os.environ.get('PTM_FONT_DIRS'):
        preferred_sans = ['Arial', *[font for font in preferred_sans if font != 'Arial']]

    mpl.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': preferred_sans,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'axes.titlesize': 12.5,
        'axes.titleweight': 'regular',
        'axes.titlepad': 7,
        'axes.labelsize': 11,
        'axes.labelcolor': '#111111',
        'axes.edgecolor': '#111111',
        'axes.linewidth': 0.85,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'xtick.color': '#111111',
        'ytick.color': '#111111',
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'xtick.major.size': 4,
        'ytick.major.size': 4,
        'legend.fontsize': 9.5,
        'legend.frameon': False,
        'grid.color': '#d9d9d9',
        'grid.linewidth': 0.6,
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
        'savefig.facecolor': 'white',
    })


def style_axis(ax, zero: str | None = None, grid_axis: str | None = None, force_grid: bool = False) -> None:
    wants_grid = grid_axis is not None and (force_grid or len(ax.lines) > 0 or len(ax.collections) > 0)
    if wants_grid:
        ax.grid(axis=grid_axis, color=BREWER_COLORS['light_gray'], linestyle='-', linewidth=0.6)
    ax.set_axisbelow(True)
    if zero == 'x':
        ax.axvline(0, color='#111111', linewidth=0.85)
    elif zero == 'y':
        ax.axhline(0, color='#111111', linewidth=0.85)
    ax.spines['left'].set_color('#111111')
    ax.spines['bottom'].set_color('#111111')
    ax.tick_params(colors='#111111')


def signed_bar_colors(values: Iterable[float], positive: str | None = None, negative: str | None = None, neutral: str | None = None) -> list[str]:
    positive = positive or BREWER_COLORS['dark_gray']
    negative = negative or BREWER_COLORS['mid_gray']
    neutral = neutral or BREWER_COLORS['light_gray']
    out = []
    for value in values:
        if pd.isna(value):
            out.append(neutral)
        elif float(value) > 0:
            out.append(positive)
        elif float(value) < 0:
            out.append(negative)
        else:
            out.append(neutral)
    return out


def set_symmetric_xlim(ax, values: Iterable[float], annotation_pad_ratio: float = 0.32, center_on_zero: bool = False) -> None:
    numeric = pd.Series(list(values), dtype='float').replace([np.inf, -np.inf], np.nan).dropna()
    if numeric.empty:
        return
    vmin = float(numeric.min())
    vmax = float(numeric.max())
    if center_on_zero:
        extent = max(abs(vmin), abs(vmax))
        extent = extent if extent > 0 else 1.0
        pad = extent * annotation_pad_ratio
        ax.set_xlim(-(extent + pad), extent + pad)
    else:
        span = max(vmax - vmin, 1.0)
        pad = span * annotation_pad_ratio
        ax.set_xlim(vmin - pad, vmax + pad)


def annotate_barh(ax, y_values: Iterable, x_values: Iterable[float], labels: Iterable[str], fontsize: float = 8.0, color: str | None = None) -> None:
    import matplotlib.patheffects as pe

    color = color or BREWER_COLORS['dark_gray']
    xmin, xmax = ax.get_xlim()
    span = max(xmax - xmin, 1.0)
    offset = span * 0.02
    for y, x, label in zip(y_values, x_values, labels):
        if pd.isna(x):
            continue
        x = float(x)
        if x >= 0:
            xpos = x + offset
            ha = 'left'
        else:
            xpos = x - offset
            ha = 'right'
        ax.text(
            xpos,
            y,
            str(label),
            va='center',
            ha=ha,
            fontsize=fontsize,
            color=color,
            clip_on=False,
            path_effects=[pe.withStroke(linewidth=2.6, foreground='white')],
        )


def count_over_total_labels(
    df: pd.DataFrame,
    target_col: str = 'target_count',
    background_col: str = 'background_count',
    include_percent: bool = False,
) -> list[str]:
    labels: list[str] = []
    for _, row in df.iterrows():
        target = row.get(target_col, np.nan)
        background = row.get(background_col, np.nan)
        if pd.notna(target) and pd.notna(background):
            target_i = int(target)
            total_i = int(target + background)
            if include_percent and total_i > 0:
                labels.append(f'{target_i}/{total_i} ({100.0 * target_i / total_i:.1f}%)')
            else:
                labels.append(f'{target_i}/{total_i}')
        elif pd.notna(target):
            labels.append(f'n={int(target)}')
        else:
            labels.append('')
    return labels


def balanced_category_subset(
    df: pd.DataFrame,
    value_col: str = 'log2_odds_ratio',
    count_col: str = 'target_count',
    min_count: int = 0,
    top_positive: int = 6,
    top_negative: int = 4,
    always_include: Iterable[str] | None = None,
) -> pd.DataFrame:
    plot_df = df.copy()
    if count_col in plot_df.columns:
        plot_df = plot_df[plot_df[count_col] >= min_count].copy()
    pos = plot_df[plot_df[value_col] >= 0].sort_values(value_col, ascending=False).head(top_positive)
    neg = plot_df[plot_df[value_col] < 0].sort_values(value_col, ascending=True).head(top_negative)
    frames = [pos, neg]
    if always_include:
        frames.append(plot_df[plot_df['category'].isin(set(always_include))])
    out = pd.concat(frames, ignore_index=True).drop_duplicates(subset=['category'])
    return out.sort_values(value_col).reset_index(drop=True)


def load_json(path: str | Path) -> dict:
    return json.loads(Path(path).read_text())


def parse_fasta(path: str | Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    current_acc = ''
    for raw_line in Path(path).read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith('>'):
            header = line[1:]
            parts = header.split('|')
            if len(parts) >= 2:
                current_acc = parts[1].strip()
            else:
                current_acc = header.split()[0]
            sequences[current_acc] = []
            continue
        if current_acc:
            sequences[current_acc].append(line.strip())
    return {acc: ''.join(chunks) for acc, chunks in sequences.items()}


def classify_by_patterns(text: str, blocks: list[dict], default: str = 'Other / unclassified') -> str:
    t = str(text).lower()
    for block in blocks:
        if any(re.search(p, t) for p in block['patterns']):
            return block['label']
    return default


def log2_odds_ratio(value: float) -> float:
    value = float(value)
    if value <= 0 or math.isnan(value):
        return math.nan
    return math.log2(value)


def format_p_value(value) -> str:
    if pd.isna(value):
        return 'NA'
    value = float(value)
    if value == 0.0:
        return '0'
    if value < 1e-3:
        return f'{value:.1e}'
    return f'{value:.3f}'


def format_compact_number(value, scientific_threshold: float = 1e-2, significant_digits: int = 1) -> str:
    if pd.isna(value):
        return 'NA'
    value = float(value)
    if value == 0.0:
        return '0'
    abs_value = abs(value)
    if abs_value < scientific_threshold:
        exponent = int(math.floor(math.log10(abs_value)))
        mantissa = abs_value / (10 ** exponent)
        decimals = max(significant_digits - 1, 0)
        mantissa = round(mantissa, decimals)
        if mantissa >= 10:
            mantissa /= 10
            exponent += 1
        mantissa_str = f'{mantissa:.{decimals}f}'.rstrip('0').rstrip('.')
        return rf'${mantissa_str}\times10^{{{exponent}}}$'
    if abs_value < 0.1:
        return f'{value:.2f}'
    if abs_value < 1:
        return f'{value:.1f}'
    return f'{value:.0f}'


def compact_q_label(value, prefix: str = 'q', significant_digits: int = 1) -> str:
    if pd.isna(value):
        return f'{prefix}=NA'
    value = float(value)
    if value == 0.0:
        return f'{prefix}\u22480'
    body = format_compact_number(value, significant_digits=significant_digits)
    return f'{prefix}={body}'


def q_threshold_note(values: Iterable[float], prefix: str = 'All shown BH q') -> str:
    series = pd.Series(list(values), dtype='float').replace([np.inf, -np.inf], np.nan).dropna()
    if series.empty:
        return ''
    max_q = float(series.max())
    if max_q == 0.0:
        return f'{prefix}\u22480'
    if max_q <= 0.05:
        return f'{prefix}<0.05'
    return f'{prefix}\u2264{format_compact_number(max_q, significant_digits=1)}'


def benjamini_hochberg(p_values: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(p_values), dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return out
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    n = len(ranked)
    adjusted = np.empty(n, dtype=float)
    cumulative = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        bh = ranked[i] * n / rank
        cumulative = min(cumulative, bh)
        adjusted[i] = min(cumulative, 1.0)
    restored = np.empty(n, dtype=float)
    restored[order] = adjusted
    out[valid] = restored
    return out


def add_q_values(df: pd.DataFrame, p_col: str = 'p_value', out_col: str = 'q_value') -> pd.DataFrame:
    out = df.copy()
    if p_col not in out.columns:
        out[out_col] = np.nan
        return out
    out[out_col] = benjamini_hochberg(out[p_col].astype(float).to_numpy())
    return out


# Structural context: MobiDB disorder and InterPro domains in one class per arginine, from inside a folded domain to deep
# in a disordered region (scripts/07 tests it; build_release writes it per site). An arginine in a disordered region
# outside domains is binned by its distance to the nearest domain edge; on a protein without domains it is far.
STRUCTURE_NEAR, STRUCTURE_MID = 10, 30
STRUCTURE_CONTEXTS = ['domain_folded', 'domain_disordered', 'ordered_outside_domain', 'idr_near_domain', 'idr_mid_domain',
                      'idr_far_from_domain', 'no_disorder_annotation']


def structure_context(disorder_class, domain_class, domain_edge_distance) -> str:
    """The class of STRUCTURE_CONTEXTS from a site's disorder class (scripts/06), domain context (scripts/05) and
    distance to the nearest domain edge."""
    if str(disorder_class) == 'no_disorder_annotation':
        return 'no_disorder_annotation'
    disordered = str(disorder_class) == 'disordered'
    if str(domain_class) == 'in_domain':
        return 'domain_disordered' if disordered else 'domain_folded'
    if not disordered:
        return 'ordered_outside_domain'
    d = pd.to_numeric(domain_edge_distance, errors='coerce')
    if pd.isna(d) or d > STRUCTURE_MID:
        return 'idr_far_from_domain'
    return 'idr_near_domain' if d <= STRUCTURE_NEAR else 'idr_mid_domain'


def per_protein_background(analysis: str, tested, background: pd.DataFrame, col: str | None = None,
                           flags: list[str] | None = None) -> pd.DataFrame:
    """The background of a residue test, per protein: for every protein whose sites the test counts (tested), the number
    of its other arginines (category '*'), and of those in each class of col (missing: 'Other / unclassified', as
    fisher_like_enrichment counts them) or with each flag column true (named without its 'motif_' prefix). The website
    repeats the test for a selection of sites from these counts: the selected methylarginines against the other
    arginines of the same proteins."""
    tested = sorted({str(a) for a in tested})
    acc = background['canonical_UniProtAC'].astype(str)
    keep = acc.isin(tested)
    parts = [pd.DataFrame({'canonical_UniProtAC': tested, 'category': '*',
                           'other_arginines': acc[keep].value_counts().reindex(tested, fill_value=0).to_numpy()})]
    if col is not None:
        cats = background.loc[keep, col].fillna('Other / unclassified').astype(str)
        parts.append(pd.DataFrame({'canonical_UniProtAC': acc[keep], 'category': cats})
                     .groupby(['canonical_UniProtAC', 'category']).size().rename('other_arginines').reset_index())
    for flag in flags or []:
        hit = keep & background[flag].fillna(False).astype(bool)
        parts.append(acc[hit].value_counts().rename('other_arginines').rename_axis('canonical_UniProtAC').reset_index()
                     .assign(category=flag.removeprefix('motif_')))
    out = pd.concat(parts, ignore_index=True)
    out.insert(0, 'analysis', analysis)
    return out[['analysis', 'canonical_UniProtAC', 'category', 'other_arginines']].sort_values(
        ['canonical_UniProtAC', 'category'], kind='stable').reset_index(drop=True)


# Sites of a release the analyses leave out altogether: neither methylarginines nor background arginines. scripts/
# release_analyses.py prepare lists them beside the union it writes (WORK/integrated/): the lower-confidence sites the
# view with lower-confidence sites does not count (in_analyses_with_lower_confidence false) and the sites kept pending,
# as build_release.py and release_layers.change_tables leave them out of both groups.
LEFT_OUT_SITES = 'human_arg_methyl_union_left_out.tsv'


def left_out_site_keys(path: str | Path) -> set[str]:
    """The site_keys (accession:Rposition) of the sites the analyses leave out, from the LEFT_OUT_SITES table beside a
    union (WORK/integrated/<union>) or, for a table a later step wrote from the union (WORK/<step>/<table>), in
    WORK/integrated/. Empty when there is no such table (a pipeline run without one leaves no site out)."""
    p = Path(path)
    for candidate in (p.parent / LEFT_OUT_SITES, p.parent.parent / 'integrated' / LEFT_OUT_SITES):
        if candidate.exists():
            return set(pd.read_csv(candidate, sep='\t', dtype=str, keep_default_na=False)['site_key'])
    return set()


def with_background(sites: pd.DataFrame, background: pd.DataFrame, col: str = 'canonical_UniProtAC') -> pd.DataFrame:
    """The sites a same-protein residue test counts: those on a protein with at least one other arginine in background
    (not a protein whose every arginine is a site, MIF in v1.1.0 with 5 sites, nor one outside the FASTA), so that the
    pooled test (Fisher) and the within-protein test (Mantel-Haenszel) of a class count the same sites."""
    return sites[sites[col].isin(set(background[col]))]


# Within-protein tests: Mantel-Haenszel with each protein a stratum (scripts/07, scripts/37).
def informative_strata(tables):
    """Strata that carry information for a Mantel-Haenszel analysis (no empty margin, at least two observations)."""
    return [(a, b, c, d) for a, b, c, d in tables if a + b + c + d >= 2 and (a + b) and (c + d) and (a + c) and (b + d)]


def cmh_statistic(tables) -> float:
    """The continuity-corrected Cochran-Mantel-Haenszel statistic of 2x2 tables (a, b, c, d): (|O - E| - 0.5)^2 / V,
    with the bracket floored at 0, so that a count within 0.5 of its expectation gives 0 (p = 1). statsmodels'
    test_null_odds(correction=True) leaves the bracket negative there, and the square then grows as V shrinks: a sparse
    test with nothing observed and almost nothing expected came out significant. NaN when V is 0."""
    o = e = v = 0.0
    for a, b, c, d in tables:
        n = a + b + c + d
        o += a
        e += (a + b) * (a + c) / n
        v += (a + b) * (c + d) * (a + c) * (b + d) / (n * n * (n - 1))
    return max(0.0, abs(o - e) - 0.5) ** 2 / v if v > 0 else float('nan')


def mantel_haenszel(tables):
    """Mantel-Haenszel common odds ratio (statsmodels StratifiedTable) and the continuity-corrected CMH test (cmh_statistic).

    Returns (odds, p, informative strata); odds may be inf when no stratum has b*c > 0."""
    strata = informative_strata(tables)
    if not strata:
        return None, None, 0
    from statsmodels.stats.contingency_tables import StratifiedTable
    st = StratifiedTable([np.array([[a, b], [c, d]], dtype=float) for a, b, c, d in strata])
    den = sum(b * c / (a + b + c + d) for a, b, c, d in strata)
    odds = float(st.oddsratio_pooled) if den > 0 else float('inf')
    # scipy's survival function, which stays accurate where 1 - cdf underflows (below ~1e-16)
    stat = cmh_statistic(strata)
    return odds, (float(chi2.sf(stat, 1)) if math.isfinite(stat) else None), len(strata)


def mantel_haenszel_ci(tables):
    """Mantel-Haenszel OR with its Robins-Breslow-Greenland 95% CI (statsmodels): (odds, lo, hi, p, informative strata)."""
    odds, p, used = mantel_haenszel(tables)
    if odds is None or not math.isfinite(odds) or odds <= 0:
        return odds, None, None, p, used
    strata = informative_strata(tables)
    from statsmodels.stats.contingency_tables import StratifiedTable
    st = StratifiedTable([np.array([[a, b], [c, d]], dtype=float) for a, b, c, d in strata])
    lo, hi = st.oddsratio_pooled_confint(alpha=0.05)
    return odds, (float(lo) if math.isfinite(lo) else None), (float(hi) if math.isfinite(hi) else None), p, used


def fisher_like_enrichment(target: pd.Series, universe: pd.Series) -> pd.DataFrame:
    target = target.fillna('Other / unclassified')
    universe = universe.fillna('Other / unclassified')
    rows = []
    target_n = len(target)
    universe_n = len(universe)
    for cat in sorted(set(target.unique()) | set(universe.unique())):
        a = int((target == cat).sum())
        b = int(target_n - a)
        c = int((universe == cat).sum())
        d = int(universe_n - c)
        # Haldane-Anscombe correction
        orr = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        rows.append({
            'category': cat,
            'target_count': a,
            'target_fraction': a / target_n if target_n else np.nan,
            'background_count': c,
            'background_fraction': c / universe_n if universe_n else np.nan,
            'odds_ratio': orr,
            'log2_odds_ratio': log2_odds_ratio(orr),
            'p_value': float(fisher_exact([[a, b], [c, d]], alternative='two-sided').pvalue),
        })
    out = add_q_values(pd.DataFrame(rows))
    return out.sort_values(['odds_ratio', 'target_count'], ascending=[False, False])


def normalize_accession(acc: str) -> str:
    acc = str(acc).strip()
    if acc in {'', 'nan', 'None'}:
        return ''
    return acc.split(';')[0].strip()


def canonicalize_uniprot_accession(acc: str) -> str:
    acc = normalize_accession(acc)
    if not acc:
        return ''
    return re.sub(r'-\d+$', '', acc)


def accession_is_isoform(acc: str) -> bool:
    acc = normalize_accession(acc)
    if not acc:
        return False
    return canonicalize_uniprot_accession(acc) != acc


def canonical_gene_name(g: str) -> str:
    g = str(g).strip()
    if g in {'', 'nan', 'None'}:
        return ''
    return g.split(';')[0].strip()


def source_family_from_base_source(source: str) -> str:
    source = str(source)
    if source == 'unip':
        return 'UniProt'
    return 'iPTMnet'


def extract_site_fields(residue: str, position) -> tuple[str, int | None, str]:
    res = str(residue).strip().upper()[:1]
    pos = pd.to_numeric(pd.Series([position]), errors='coerce').iloc[0]
    if pd.isna(pos):
        return res, None, ''
    pos = int(pos)
    return res, pos, f'{res}{pos}'


def clean_peptide_field(x: str) -> str:
    s = str(x)
    s = s.replace('[', '').replace(']', '').replace("'", '')
    s = s.replace('(ox)', '')
    return s.strip()


def clean_aa_sequence(x: str) -> str:
    return ''.join(ch for ch in str(x).upper() if 'A' <= ch <= 'Z')


def split_multi_value_field(x: str, separators: str = ';,') -> list[str]:
    if str(x) in {'', 'nan', 'None'}:
        return []
    pattern = '[' + re.escape(separators) + ']'
    return [token.strip() for token in re.split(pattern, str(x)) if token.strip()]


def centered_anchor_indices(sequence: str, residue: str = 'R') -> list[int]:
    seq = clean_aa_sequence(sequence)
    if not seq:
        return []
    midpoint = (len(seq) + 1) // 2
    if 1 <= midpoint <= len(seq) and seq[midpoint - 1] == residue:
        return [midpoint]
    residue_positions = [idx + 1 for idx, aa in enumerate(seq) if aa == residue]
    if not residue_positions:
        return []
    distances = [abs(idx - midpoint) for idx in residue_positions]
    min_distance = min(distances)
    return [idx for idx in residue_positions if abs(idx - midpoint) == min_distance]


def parse_modified_peptide(x: str, residue: str = 'R') -> tuple[str, list[int]]:
    text = str(x).strip()
    clean = []
    modified_positions: list[int] = []
    i = 0
    while i < len(text):
        char = text[i]
        upper = char.upper()
        if 'A' <= upper <= 'Z':
            clean.append(upper)
            aa_index = len(clean)
            i += 1
            if i < len(text) and text[i] == '(':
                end = text.find(')', i)
                if end == -1:
                    end = i
                annotation = text[i + 1:end].strip().lower()
                if upper == residue and annotation not in {'', 'ox'}:
                    modified_positions.append(aa_index)
                i = end + 1
            continue
        i += 1
    return ''.join(clean), modified_positions


def mass_localised(peptide: str, marks: dict[int, str], methylatable: str) -> str:
    """How a peptide's mass alone localises its methyl groups, or '' when it does not. peptide: the plain sequence;
    methylatable: the residues the study's search could methylate (R, or RK); marks: 1-based residue -> 'Rme1' or
    'Rme2' for the methylated arginines. 'one_methylatable_residue': the peptide holds one residue the search could
    methylate; 'every_arginine_dimethylated': every such residue is an arginine and every one is dimethylated, so the
    mass leaves no other arrangement."""
    candidates = [k for k, aa in enumerate(peptide, 1) if aa in methylatable]
    if len(candidates) == 1:
        return 'one_methylatable_residue'
    if candidates and all(peptide[k - 1] == 'R' and marks.get(k) == 'Rme2' for k in candidates):
        return 'every_arginine_dimethylated'
    return ''


def normalize_arg_methyl_state(state: str) -> str:
    s = str(state).strip().lower()
    if s in {'', 'nan', 'none'}:
        return ''
    if s in {'me', 'mono', 'mono-methyl', 'monomethyl', 'monomethylation', 'mma', 'rme1'}:
        return 'Rme1'
    if s in {'di', 'di-methyl', 'dimethyl', 'dimethylation', 'rme2'}:
        return 'Rme2'
    if any(token in s for token in ['asymmetric', 'adma', 'rme2a']):
        return 'Rme2a'
    if any(token in s for token in ['symmetric', 'sdma', 'rme2s']):
        return 'Rme2s'
    if 'mono' in s:
        return 'Rme1'
    if 'di' in s or 'dimethyl' in s:
        return 'Rme2'
    return str(state).strip()


def build_fuzzy_clusters(df: pd.DataFrame, tolerance: int = 2) -> pd.DataFrame:
    rows = []
    required = ['substrate_UniProtAC', 'residue', 'position', 'source_family']
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f'missing columns for clustering: {missing}')
    for (acc, residue), g in df.dropna(subset=['position']).groupby(['substrate_UniProtAC', 'residue']):
        g = g.sort_values(['position', 'source_family']).reset_index(drop=True)
        parent = list(range(len(g)))

        def find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(a: int, b: int) -> None:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra

        pos = g['position'].astype(int).tolist()
        fam = g['source_family'].astype(str).tolist()
        for i in range(len(g)):
            for j in range(i + 1, min(i + 50, len(g))):
                if pos[j] - pos[i] > tolerance:
                    break
                if fam[i] != fam[j]:
                    union(i, j)

        comps: dict[int, list[int]] = defaultdict(list)
        for i in range(len(g)):
            comps[find(i)].append(i)

        for idxs in comps.values():
            sub = g.iloc[idxs].copy()
            srcs = sorted(set(sub['source_family'].astype(str)))
            rows.append({
                'substrate_UniProtAC': acc,
                'residue': residue,
                'cluster_pos_min': int(sub['position'].min()),
                'cluster_pos_max': int(sub['position'].max()),
                'cluster_pos_rep': int(round(float(sub['position'].median()))),
                'cluster_site_rep': f'{residue}{int(round(float(sub["position"].median())))}',
                'source_family_count_fuzzy': len(srcs),
                'source_families_fuzzy': ';'.join(srcs),
                'rows_in_cluster': len(sub),
            })
    return pd.DataFrame(rows)


def summarize_sources(values: Iterable[str]) -> str:
    return ';'.join(sorted(set(str(v) for v in values if str(v) not in {'', 'nan', 'None'})))


def consensus_or_mixed(values: Iterable[str], default: str = '') -> str:
    uniq = sorted(set(str(v) for v in values if str(v) not in {'', 'nan', 'None'}))
    if not uniq:
        return default
    if len(uniq) == 1:
        return uniq[0]
    return 'mixed'


def sequence_window(sequence: str, position: int, flank: int, pad: str = '_') -> str:
    seq = clean_aa_sequence(sequence)
    if not seq or position < 1:
        return ''
    left = max(0, position - flank - 1)
    right = min(len(seq), position + flank)
    window = seq[left:right]
    missing_left = flank - (position - 1 - left)
    missing_right = flank - (right - position)
    return (pad * missing_left) + window + (pad * missing_right)


def residue_background_from_sequences(
    sequences: dict[str, str],
    residue: str,
    accessions: Iterable[str] | None = None,
    flank: int = 0,
) -> pd.DataFrame:
    residue = str(residue).upper()
    selected = set(accessions) if accessions is not None else None
    rows = []
    for acc, sequence in sequences.items():
        if selected is not None and acc not in selected:
            continue
        seq = clean_aa_sequence(sequence)
        for idx, aa in enumerate(seq, start=1):
            if aa != residue:
                continue
            row = {
                'canonical_UniProtAC': acc,
                'position': idx,
                'residue': residue,
                'site': f'{residue}{idx}',
                'site_key': f'{acc}:{residue}{idx}',
            }
            if flank > 0:
                row[f'window_{2 * flank + 1}'] = sequence_window(seq, idx, flank=flank)
            rows.append(row)
    return pd.DataFrame(rows)


def canonical_site_table(
    df: pd.DataFrame,
    accession_col: str = 'canonical_UniProtAC',
    position_col: str = 'corrected_position',
    residue_col: str = 'residue',
) -> pd.DataFrame:
    out = df.copy()
    out[accession_col] = out[accession_col].astype(str).map(canonicalize_uniprot_accession)
    out[position_col] = pd.to_numeric(out[position_col], errors='coerce')
    out = out[out[position_col].notna()].copy()
    out[position_col] = out[position_col].astype(int)
    out[residue_col] = out[residue_col].astype(str).str.upper().str[:1]
    out['site'] = out[residue_col] + out[position_col].astype(str)
    out['site_key'] = out[accession_col] + ':' + out['site']
    return out


# ---- Staging inputs ----------------------------------------------------------------------------------------------------
# What every staging script needs: the MD5 of what it wrote, its rows in inputs/MANIFEST.tsv (and a layer's raw downloads
# in raw_files.tsv), the Entrez gene of a UniProt entry, and whether a position of a reference sequence holds a residue.

def md5_file(path: str | Path) -> str:
    import hashlib
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def upsert_tsv_rows(path: str | Path, rows: list[dict], fmt=str) -> None:
    """Put rows (dicts keyed by the file's columns; a missing column is left empty) into a TSV keyed by its first column:
    a row whose key is already there replaces it in place, a new row goes after the last row of its folder (the key up to
    its first '/'), or at the end. The other rows and their order stay as they are, so re-staging the same files gives the
    same file."""
    path = Path(path)
    lines = path.read_text(encoding='utf-8').splitlines()
    head = lines[0].split('\t')
    key = head[0]
    text = {str(r[key]): '\t'.join(fmt(r.get(c, '')) for c in head) for r in rows}
    body, done = [], set()
    for ln in lines[1:]:
        k = ln.split('\t')[0]
        if k in text:
            body.append(text[k])
            done.add(k)
        else:
            body.append(ln)
    for r in rows:
        k = str(r[key])
        if k in done:
            continue
        folder = k.split('/')[0] + '/' if '/' in k else None
        last = max((i for i, ln in enumerate(body) if folder and ln.startswith(folder)), default=len(body) - 1)
        body.insert(last + 1, text[k])
        done.add(k)
    path.write_text('\n'.join([lines[0]] + body) + '\n', encoding='utf-8')


def update_manifest(manifest: str | Path, rows: list[dict], fmt=str) -> None:
    """The rows of the files a staging run wrote, into inputs/MANIFEST.tsv (upsert_tsv_rows)."""
    upsert_tsv_rows(manifest, rows, fmt)


def refresh_manifest_row(manifest: str | Path, path: str | Path, rel: str, **fields) -> None:
    """The row of a file a staging run rewrote in part (a folder's raw_files.tsv) in inputs/MANIFEST.tsv: its bytes, MD5
    and rows as they are now, its other fields kept unless given."""
    path = Path(path)
    table = pd.read_csv(manifest, sep='\t', dtype=str, keep_default_na=False)
    old = table[table['file'] == rel]
    row = old.iloc[0].to_dict() if len(old) else {'file': rel, 'redistribution': 'git'}
    with open(path, encoding='utf-8') as fh:
        n = sum(1 for _ in fh) - 1
    row.update({'bytes': path.stat().st_size, 'md5': md5_file(path), 'rows': n}, **fields)
    update_manifest(manifest, [row])


def entrez_by_symbol(path: str | Path) -> dict[str, list[str]]:
    """Gene symbol -> its Entrez genes in org.Hs.eg.db (inputs/disease/orghsegdb_*_symbol_entrez.tsv), in file order."""
    sym = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    return sym.groupby('symbol', sort=False)['entrez_id'].agg(list).to_dict()


def own_entrez(path: str | Path) -> set[tuple[str, str]]:
    """The (symbol, Entrez) pairs where org.Hs.eg.db gives the gene that symbol as its own (entrez_symbol true)."""
    sym = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    own = sym[sym['entrez_symbol'] == 'true']
    return set(zip(own['symbol'], own['entrez_id']))


def entrez_of_entries(entries: str | Path, symbol_entrez: str | Path, allowed: set[str] | None = None,
                      reference: dict | None = None) -> tuple[dict[str, str], dict[str, str], dict[str, int]]:
    """Each reviewed UniProt entry's Entrez gene, as the disease tests map genes: the first of its primary gene names
    (UniProt entries table) to its Entrez genes in org.Hs.eg.db. With allowed (e.g. the genes a screen holds), only those
    count; of several, the one org.Hs.eg.db gives the name as its own. Returns (accession -> Entrez, accession -> gene
    name, counts of the entries left without one, by reason); reference limits the entries to its accessions."""
    table = pd.read_csv(entries, sep='\t', dtype=str, keep_default_na=False)
    of, own = entrez_by_symbol(symbol_entrez), own_entrez(symbol_entrez)
    acc_entrez, acc_gene, dropped = {}, {}, defaultdict(int)
    for acc, names in zip(table['Entry'], table['Gene Names (primary)']):
        if (reference is not None and acc not in reference) or not names:
            continue
        g = names.split(';')[0].strip()
        cands = of.get(g, [])
        hit = [e for e in cands if allowed is None or e in allowed]
        if len(hit) > 1:
            hit = [e for e in hit if (g, e) in own] or hit
        acc_gene[acc] = g
        if len(hit) == 1:
            acc_entrez[acc] = hit[0]
        elif not cands:
            dropped['without_entrez'] += 1
        else:
            dropped['entrez_not_allowed' if not hit else 'several_entrez'] += 1
    return acc_entrez, acc_gene, dict(dropped)


def residue_is(sequences: dict[str, str], accession: str, position: int, residue: str = 'R') -> bool:
    """Whether the reference sequence of accession has residue at position (1-based)."""
    seq = sequences.get(accession)
    return bool(seq) and 1 <= int(position) <= len(seq) and seq[int(position) - 1] == residue

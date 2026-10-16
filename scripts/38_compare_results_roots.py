#!/usr/bin/env python3
"""Compare two PTM_results roots file by file, for example the published run and a rerun after a fix.

    python scripts/38_compare_results_roots.py --old PTM_results --new ../PTM_results_remapfix \
        --outdir ../PTM_results_remapfix/rerun_comparison

Tables (.tsv, .csv) are matched row by row on a key found from their columns (site keys, accession + position,
or the leading text columns), else compared in order when their shapes agree. For each table it reports rows
added and removed, cells changed, the largest numeric change, p/q values that cross the significance threshold
and effect sizes whose sign flips. Figures (.png, .pdf, .svg) and other files are compared by content hash.

Writes:
  file_summary.tsv     one row per file in either root
  changes_detail.tsv   significance flips, sign flips, largest numeric changes and changed keys
  site_union_diff.tsv  sites added to or removed from integrated/human_arg_methyl_union_dedup_by_site.tsv
  report.md            the union difference, headline summary tables, flips, and changed files by folder
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

TABLE_EXT = {'.tsv': '\t', '.csv': ','}
UNION = Path('integrated/human_arg_methyl_union_dedup_by_site.tsv')
KEY_CANDIDATES = [
    ['corrected_site_key'], ['site_key_canonical'], ['canonical_UniProtAC', 'corrected_position'],
    ['canonical_UniProtAC', 'position'], ['UniProtAC', 'position'], ['accession', 'position'], ['site_key'],
]
SIG_COL = re.compile(r'(^|[_\s])(p|q|pval|pvalue|p_value|padj|p_adj|p\.adjust|qvalue|q_value|fdr|bh)([_\s]|$)', re.I)
EFFECT_COL = re.compile(r'(log2|log_?odds|odds_?ratio|(^|_)or($|_)|fold|enrichment|effect|delta|z_?score|(^|_)z($|_))', re.I)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def read_table(path: Path) -> pd.DataFrame | None:
    try:
        return pd.read_csv(path, sep=TABLE_EXT[path.suffix], low_memory=False, dtype=str, keep_default_na=False)
    except Exception:  # noqa: BLE001 - malformed or non-tabular text is compared by hash instead
        return None


def numeric(series: pd.Series) -> pd.Series | None:
    """The column as floats when at least 90% of its non-empty values parse, else None."""
    values = series.replace('', np.nan)
    parsed = pd.to_numeric(values, errors='coerce')
    filled = values.notna().sum()
    if filled == 0 or parsed.notna().sum() < 0.9 * filled:
        return None
    return parsed


def find_key(old: pd.DataFrame, new: pd.DataFrame) -> list[str] | None:
    common = [c for c in old.columns if c in new.columns]
    for cols in KEY_CANDIDATES:
        if all(c in common for c in cols) and not old.duplicated(cols).any() and not new.duplicated(cols).any():
            return cols
    text = [c for c in common if numeric(old[c]) is None and numeric(new[c]) is None]
    for k in range(1, min(3, len(text)) + 1):
        cols = text[:k]
        if not old.duplicated(cols).any() and not new.duplicated(cols).any():
            return cols
    return None


def key_strings(df: pd.DataFrame, cols: list[str]) -> pd.Series:
    return df[cols].astype(str).agg(' | '.join, axis=1)


def compare_tables(rel: str, old: pd.DataFrame, new: pd.DataFrame, alpha: float, details: list[dict]) -> dict:
    out = {'old_rows': len(old), 'new_rows': len(new), 'columns_added': ';'.join(c for c in new.columns if c not in old.columns),
           'columns_removed': ';'.join(c for c in old.columns if c not in new.columns)}
    key = find_key(old, new)
    if key:
        ko, kn = key_strings(old, key), key_strings(new, key)
        o, n = old.set_axis(ko, axis=0), new.set_axis(kn, axis=0)
        added, removed = sorted(set(kn) - set(ko)), sorted(set(ko) - set(kn))
        shared = [k for k in ko if k in set(kn)]
        o, n = o.loc[shared], n.loc[shared]
        out.update(match='key: ' + '+'.join(key), rows_added=len(added), rows_removed=len(removed))
        for label, keys in (('row added', added), ('row removed', removed)):
            for k in keys[:25]:
                details.append({'file': rel, 'kind': label, 'row': k, 'column': '', 'old': '', 'new': ''})
    elif len(old) == len(new):
        o, n = old.reset_index(drop=True), new.reset_index(drop=True)
        out.update(match='row order', rows_added=0, rows_removed=0)
    else:
        out.update(match='none (no unique key, different row counts)', rows_added='', rows_removed='')
        return out

    cells = 0
    max_abs, max_col = 0.0, ''
    sig_flips = sign_flips = 0
    for col in [c for c in o.columns if c in n.columns and c not in (key or [])]:
        a, b = o[col], n[col]
        changed = a.ne(b)
        if not changed.any():
            continue
        na, nb = numeric(a), numeric(b)
        if na is not None and nb is not None:
            diff = (nb - na).abs()
            both = na.notna() & nb.notna()
            real = changed & ~(both & np.isclose(na, nb, rtol=1e-9, atol=1e-12))
            cells += int(real.sum())
            if both.any() and diff[both].max() > max_abs:
                max_abs, max_col = float(diff[both].max()), col
            if SIG_COL.search(col):
                flip = both & ((na <= alpha) != (nb <= alpha))
                sig_flips += int(flip.sum())
                for idx in flip[flip].index[:50]:
                    details.append({'file': rel, 'kind': f'significance flip (alpha {alpha})', 'row': str(idx), 'column': col, 'old': a[idx], 'new': b[idx]})
            if EFFECT_COL.search(col):
                # Sign flips of effect sizes; log-scale columns flip around 0, ratios around 1.
                centre = 1.0 if re.search(r'(ratio|fold|(^|_)or($|_))', col, re.I) and not re.search(r'log', col, re.I) else 0.0
                flip = both & (np.sign(na - centre) * np.sign(nb - centre) < 0)
                sign_flips += int(flip.sum())
                for idx in flip[flip].index[:50]:
                    details.append({'file': rel, 'kind': 'effect sign flip', 'row': str(idx), 'column': col, 'old': a[idx], 'new': b[idx]})
            top = diff[real & both].sort_values(ascending=False).head(3)
            for idx, v in top.items():
                details.append({'file': rel, 'kind': 'numeric change', 'row': str(idx), 'column': col, 'old': a[idx], 'new': b[idx], 'abs_diff': v})
        else:
            cells += int(changed.sum())
    out.update(cells_changed=cells, max_abs_diff=round(max_abs, 6) if max_abs else '', max_abs_diff_column=max_col,
               significance_flips=sig_flips, effect_sign_flips=sign_flips)
    return out


def union_diff(old_root: Path, new_root: Path) -> tuple[pd.DataFrame, list[str]]:
    lines: list[str] = []
    po, pn = old_root / UNION, new_root / UNION
    if not (po.exists() and pn.exists()):
        return pd.DataFrame(), ['Site union not found in both roots.']
    o, n = read_table(po), read_table(pn)
    key = ['canonical_UniProtAC', 'corrected_position']
    so = set(map(tuple, o[key].values)); sn = set(map(tuple, n[key].values))
    gene = {}
    for df in (o, n):
        gene.update(dict(zip(df['canonical_UniProtAC'], df.get('substrate_genename', df['canonical_UniProtAC']))))
    rows = [{'change': 'added', 'accession': a, 'gene': gene.get(a, ''), 'position': int(p)} for a, p in sn - so] + \
           [{'change': 'removed', 'accession': a, 'gene': gene.get(a, ''), 'position': int(p)} for a, p in so - sn]
    diff = pd.DataFrame(rows).sort_values(['gene', 'accession', 'position', 'change']) if rows else pd.DataFrame(columns=['change', 'accession', 'gene', 'position'])
    po_, pn_ = Counter(a for a, _ in so), Counter(a for a, _ in sn)
    lines += ['## Site union', '',
              f'| | Old | New |', '|---|---|---|',
              f'| Sites | {len(so):,} | {len(sn):,} |',
              f'| Proteins | {len(po_):,} | {len(pn_):,} |',
              f'| Sites added | | {len(sn - so):,} |', f'| Sites removed | | {len(so - sn):,} |', '']
    changed = [(gene.get(a, ''), a, po_[a], pn_[a]) for a in set(po_) | set(pn_) if po_[a] != pn_[a]]
    changed.sort(key=lambda r: -abs(r[3] - r[2]))
    lines += [f'{len({a for a, _ in so ^ sn}):,} proteins have a different site set; {len(changed):,} have a different site count. Largest changes:', '',
              '| Gene | Accession | Old sites | New sites |', '|---|---|---|---|']
    lines += [f'| {g} | {a} | {x} | {y} |' for g, a, x, y in changed[:20]] + ['']
    if 'source_labels_exact' in o.columns:
        def per_source(df):
            c = Counter()
            for v in df['source_labels_exact']:
                for s in set(filter(None, str(v).split(';'))):
                    c[s.strip()] += 1
            return c
        co, cn = per_source(o), per_source(n)
        lines += ['Sites per source label:', '', '| Source | Old | New | Change |', '|---|---|---|---|']
        lines += [f'| {s} | {co[s]:,} | {cn[s]:,} | {cn[s] - co[s]:+,} |' for s in sorted(set(co) | set(cn), key=lambda s: -cn[s])] + ['']
    return diff, lines


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--old', required=True, help='results root of the reference run')
    ap.add_argument('--new', required=True, help='results root of the new run')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--alpha', type=float, default=0.05, help='threshold for p/q significance flips')
    ap.add_argument('--skip', nargs='*', default=['rerun_comparison', 'site_explorer_data'], help='top-level folders to ignore')
    args = ap.parse_args()
    old_root, new_root, outdir = Path(args.old).resolve(), Path(args.new).resolve(), Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    def files(root: Path) -> set[str]:
        out = set()
        for p in root.rglob('*'):
            rel = p.relative_to(root)
            if p.is_file() and rel.parts[0] not in args.skip and p.name != 'RUN_INFO.txt' and not (outdir == p.parent or outdir in p.parents):
                out.add(str(rel))
        return out

    fo, fn = files(old_root), files(new_root)
    summary, details = [], []
    for rel in sorted(fo | fn):
        row = {'file': rel, 'folder': str(Path(rel).parent), 'type': Path(rel).suffix.lstrip('.')}
        if rel not in fn:
            row['status'] = 'only in old'
        elif rel not in fo:
            row['status'] = 'only in new'
        else:
            po, pn = old_root / rel, new_root / rel
            if sha256(po) == sha256(pn):
                row['status'] = 'identical'
            elif po.suffix in TABLE_EXT:
                o, n = read_table(po), read_table(pn)
                if o is None or n is None:
                    row['status'] = 'changed (unreadable as a table)'
                else:
                    row.update(compare_tables(rel, o, n, args.alpha, details))
                    row['status'] = 'changed' if (row.get('cells_changed') or row.get('rows_added') or row.get('rows_removed')
                                                   or row['old_rows'] != row['new_rows'] or row['columns_added'] or row['columns_removed']) else 'identical values'
            elif po.suffix == '.json':
                try:
                    row['status'] = 'identical values' if json.loads(po.read_text()) == json.loads(pn.read_text()) else 'changed'
                except ValueError:
                    row['status'] = 'changed'
            else:
                row['status'] = 'changed (content hash)'
        summary.append(row)

    summ = pd.DataFrame(summary)
    summ.to_csv(outdir / 'file_summary.tsv', sep='\t', index=False)
    det = pd.DataFrame(details, columns=['file', 'kind', 'row', 'column', 'old', 'new', 'abs_diff'])
    det.to_csv(outdir / 'changes_detail.tsv', sep='\t', index=False)
    diff, union_lines = union_diff(old_root, new_root)
    diff.to_csv(outdir / 'site_union_diff.tsv', sep='\t', index=False)

    lines = ['# Results comparison', '', f'Old: `{old_root}`  ', f'New: `{new_root}`', '']
    for root, label in ((new_root, 'New run'),):
        info = root / 'RUN_INFO.txt'
        if info.exists():
            lines += [f'{label} info:', '', '```', info.read_text().strip(), '```', '']
    counts = summ['status'].value_counts()
    lines += ['## Files', '', '| Status | Files |', '|---|---|'] + [f'| {s} | {c} |' for s, c in counts.items()] + ['']
    lines += union_lines

    # Small headline tables (summary.tsv and the like) shown side by side where they changed.
    lines += ['## Headline summary tables', '']
    for rel in sorted(r for r in fo & fn if Path(r).name.endswith('summary.tsv')):
        o, n = read_table(old_root / rel), read_table(new_root / rel)
        if o is None or n is None or o.equals(n) or len(o) > 40 or len(n) > 40:
            continue
        lines += [f'### `{rel}`', '']
        if len(o) == 1 and len(n) == 1:
            lines += ['| Field | Old | New |', '|---|---|---|']
            lines += [f'| {c} | {o.at[0, c] if c in o else ""} | {n.at[0, c] if c in n else ""} |' for c in dict.fromkeys(list(o.columns) + list(n.columns))
                      if (o.at[0, c] if c in o else None) != (n.at[0, c] if c in n else None)]
        else:
            key = find_key(o, n)
            if key:
                m = o.merge(n, on=key, how='outer', suffixes=(' (old)', ' (new)'))
                show = [c for c in m.columns if c.endswith(' (old)') and m[c].ne(m[c.replace(' (old)', ' (new)')]).any()]
                cols = key + [x for c in show for x in (c, c.replace(' (old)', ' (new)'))]
                changed = m[m[[c for c in cols if c not in key]].astype(str).apply(lambda r: any(r.iloc[i] != r.iloc[i + 1] for i in range(0, len(r), 2)), axis=1)] if show else m.iloc[0:0]
                lines += ['| ' + ' | '.join(cols) + ' |', '|' + '---|' * len(cols)]
                lines += ['| ' + ' | '.join(str(v) for v in r) + ' |' for r in changed[cols].head(40).values.tolist()]
            else:
                lines += ['(row order changed; see file_summary.tsv)']
        lines += ['']

    for kind, title in (('significance flip', f'Significance flips (alpha {args.alpha})'), ('effect sign flip', 'Effect sign flips')):
        sub = det[det['kind'].str.startswith(kind)]
        lines += [f'## {title}', '']
        if sub.empty:
            lines += ['None.', '']
            continue
        # Files with few flips (study-wide tables) first; long per-protein tables last.
        per_file = sub.groupby('file').size().sort_values()
        sub = sub.assign(_n=sub['file'].map(per_file)).sort_values(['_n', 'file'], kind='stable')
        lines += [f'{len(sub):,} cells in {sub["file"].nunique():,} files; tables with the fewest flips first (first 60 rows; all in changes_detail.tsv).', '',
                  '| File | Flips |', '|---|---|'] + [f'| {f} | {n} |' for f, n in per_file.items()] + ['',
                  '| File | Row | Column | Old | New |', '|---|---|---|---|---|']
        lines += [f'| {r.file} | {r.row} | {r.column} | {r.old} | {r.new} |' for r in sub.head(60).itertuples()] + ['']

    changed = summ[summ['status'].str.startswith('changed') | summ['status'].str.startswith('only')]
    lines += ['## Changed files by folder', '', '| Folder | Changed | Of | Significance flips | Effect sign flips |', '|---|---|---|---|---|']
    for folder, g in summ.groupby('folder'):
        ch = g[g['status'].str.startswith('changed') | g['status'].str.startswith('only')]
        if len(ch):
            sf = pd.to_numeric(g.get('significance_flips', pd.Series(dtype=float)), errors='coerce').fillna(0).sum()
            ef = pd.to_numeric(g.get('effect_sign_flips', pd.Series(dtype=float)), errors='coerce').fillna(0).sum()
            lines += [f'| {folder} | {len(ch)} | {len(g)} | {int(sf)} | {int(ef)} |']
    lines += ['', 'Figures are compared by content hash only; regenerate and view them for visual differences.', '']
    (outdir / 'report.md').write_text('\n'.join(lines))
    print(f'{len(summ):,} files compared: ' + ', '.join(f'{c} {s}' for s, c in counts.items()))
    print(f'report: {outdir / "report.md"}')


if __name__ == '__main__':
    main()

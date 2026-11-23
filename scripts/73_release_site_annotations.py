#!/usr/bin/env python3
"""Write a release that carries the previous release forward and adds, for every site set, a per-site table of the DepMap
cell lines changed at each site and AlphaMissense's predicted pathogenicity there; the disease families of each ClinVar
variation; and the scope of the site-set propensity test.

    python scripts/73_release_site_annotations.py --version v1.8.0 --previous v1.7.0 --release-date 2026-11-23

It carries every file of releases/<previous> under this release's names, as scripts/65_carry_release.py does, and then:

  supporting/site_annotations_<set>_<version>.tsv   new, per site set (config/site_sets.json): one row per site (accession,
                                                    position, residue), the DepMap 24Q4 cell lines with a missense or stop
                                                    change there (names and model IDs, from the change tables and
                                                    depmap_24q4_change_models), and AlphaMissense's mean and highest
                                                    pathogenicity over the residue's 19 substitutions and how many it
                                                    classes pathogenic (inputs/layers/alphamissense_residues.tsv.gz,
                                                    scripts/72). The site tables themselves are carried byte for byte (the
                                                    methylarginine table unchanged since v1.4.3): join on accession,
                                                    position and residue.
  supporting/clinvar_changes_<X>_<version>.tsv.gz   the carried table with one column appended, disease_families: the
                                                    variation's conditions split on '|' and ';' (without 'not specified',
                                                    'not provided' and 'see cases'), each classed by
                                                    scripts/disease_families.py classify_clinvar(), 'Other / unassigned'
                                                    dropped, the families sorted and joined by ';'
  supporting/site_set_propensity_<version>.tsv      the carried table with three columns appended, sites_total,
                                                    sites_weighted and k_total (scripts/propensity.py set_scope), computed
                                                    from the carried per-site tables, which hold every known site

Every other carried file keeps its bytes, and the appended tables keep every earlier column as it was. It then writes
data_dictionary.tsv (the previous one under the new names, with the new columns and tables), README.md, CHANGES.md,
release.json and MD5SUMS. Reader-facing prose is in American spelling.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import propensity as P  # noqa: E402
from common import md5_file  # noqa: E402
from disease_families import CLINVAR_FAMILIES, classify_clinvar  # noqa: E402


def _load(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, HERE / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


x65 = _load('carry65', '65_carry_release.py')
x68 = _load('release68', '68_release_cancer_changes.py')
LAYERS = ROOT / 'inputs' / 'layers'
RESIDUES = ['R', 'K', 'S', 'T', 'Y']
SOURCES = {'tcga_mc3': 'TCGA MC3', 'depmap_24q4': 'DepMap 24Q4', 'clinvar': 'ClinVar', 'gnomad_v4_1': 'gnomAD v4.1'}
DROP_CONDITIONS = {'not specified', 'not provided', 'see cases'}
UNASSIGNED = 'Other / unassigned'
TITLE = ('Human post-translational modification sites (arginine methylation, lysine methylation, acetylation and other acylations, '
         'ubiquitination, SUMOylation, phosphorylation) with interaction context, cancer and germline variants, and mutation propensity')
ANNOTATION_COLUMNS = ['canonical_UniProtAC', 'position', 'residue', 'depmap_lines', 'depmap_cell_lines', 'depmap_models',
                      'am_mean_pathogenicity', 'am_max_pathogenicity', 'am_pathogenic']
AM_NOTE = ('AlphaMissense (Cheng et al., Science 381:eadg7492, 2023; predictions doi:10.5281/zenodo.10813168, CC BY 4.0 by the Zenodo '
           'record, the file header still naming CC BY-NC-SA 4.0) predicts the pathogenicity of each amino-acid substitution; its '
           'scores are less reliable in intrinsically disordered regions, where many modification sites lie')


def disease_families(conditions: str, cache: dict) -> str:
    if conditions not in cache:
        names = {x.strip() for x in re.split(r'[|;]', conditions) if x.strip() and x.strip().lower() not in DROP_CONDITIONS}
        fams = {classify_clinvar(n)[0] for n in names} - {UNASSIGNED}
        cache[conditions] = ';'.join(sorted(fams))
    return cache[conditions]


def code_commit() -> str:
    paths = ['scripts/73_release_site_annotations.py', 'scripts/72_stage_alphamissense.py', 'scripts/disease_families.py', 'scripts/propensity.py',
             'scripts/65_carry_release.py', 'scripts/68_release_cancer_changes.py', 'config/site_sets.json']
    try:
        dirty = subprocess.run(['git', 'status', '--porcelain', '--', *paths], capture_output=True, text=True, check=True, cwd=ROOT).stdout
        sha = subprocess.run(['git', 'log', '-1', '--format=%H', '--', *paths], capture_output=True, text=True, check=True, cwd=ROOT).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return 'unknown'
    return sha + (' (scripts modified)' if dirty.strip() else '')


def dictionary_rows(version: str, config: dict) -> tuple[list[str], dict]:
    """New tables' rows, and rows to insert after the last row of an appended table (file -> rows)."""
    ref = 'reviewed UniProtKB accession (canonical sequence, UniProt 2026_03)'
    rows = []
    for name in config['site_sets']:
        f = f'supporting/site_annotations_{name}_{version}.tsv'
        for c in (('canonical_UniProtAC', 'string', '', '', ref),
                  ('position', 'integer', '', '', 'position of the site in the canonical sequence'),
                  ('residue', 'string', '', '', 'the residue of the site'),
                  ('depmap_lines', 'integer', '', '', 'DepMap 24Q4 cell lines with a missense or stop change at the site (the change tables)'),
                  ('depmap_cell_lines', 'string', '', '', 'their names (DepMap Model.csv, via depmap_24q4_change_models), sorted, separated by | (one '
                   'name, RS4;11, holds a semicolon)'),
                  ('depmap_models', 'string', '', '', 'their DepMap model IDs, in the same order, separated by |'),
                  ('am_mean_pathogenicity', 'number', '0-1', '', 'AlphaMissense predicted pathogenicity, mean over the residue\'s 19 substitutions; '
                   'empty where AlphaMissense\'s sequence differs from the canonical one or the entry is missing. ' + AM_NOTE),
                  ('am_max_pathogenicity', 'number', '0-1', '', 'the highest AlphaMissense pathogenicity of the 19 substitutions'),
                  ('am_pathogenic', 'integer', '0-19', '', 'how many of the 19 substitutions AlphaMissense classes pathogenic')):
            rows.append('\t'.join([f, *c]))
    after = {}
    fams = '; '.join(sorted(set(x for x in __import__('disease_families').ALL_FAMILIES + CLINVAR_FAMILIES) - {UNASSIGNED}))
    for res in RESIDUES:
        after[f'supporting/clinvar_changes_{res}_{version}.tsv.gz'] = [
            '\t'.join([f'supporting/clinvar_changes_{res}_{version}.tsv.gz', 'disease_families', 'string', fams, '',
                       'the disease families of the variation\'s conditions: conditions split on | and ; without "not specified", "not provided" '
                       'and "see cases", each classed by scripts/disease_families.py classify_clinvar() (ClinVar rules, then the DisGeNET '
                       'rules of classify()), "Other / unassigned" left out, sorted, separated by semicolons; empty when no condition gets a family'])]
    f = f'supporting/site_set_propensity_{version}.tsv'
    after[f] = ['\t'.join([f, *c]) for c in (
        ('sites_total', 'integer', '', '', 'the set\'s sites of this residue type, on every protein'),
        ('sites_weighted', 'integer', '', '', 'of those, the sites on proteins with weights (a coding sequence that translates to the canonical '
         'sequence): the sites observed and expected cover (sites counts those of them on proteins with a change)'),
        ('k_total', 'integer', '', '', 'the sum of k over all sites_total, proteins without weights included (observed is the sum over sites_weighted)'))]
    return rows, after


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--version', required=True)
    ap.add_argument('--previous', required=True)
    ap.add_argument('--release-date', required=True)
    ap.add_argument('--config', default=str(ROOT / 'config' / 'site_sets.json'))
    args = ap.parse_args()
    version, previous = args.version, args.previous
    prev, out = ROOT / 'releases' / previous, ROOT / 'releases' / version
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    config = json.loads(Path(args.config).read_text(encoding='utf-8'))
    arginine_sets = [n for n, s in config['site_sets'].items() if set(s['residues']) <= {'R'}]
    x56 = _load('build56', '56_build_context_marks.py')

    carried, edited = x65.carry(prev, out, previous, version, arginine_sets)
    sites = x65.read(out / f'human_methylarginine_sites_{version}.tsv')
    prev_changes = next((prev / 'supporting').glob(f'site_changes_*_to_{previous}.tsv'))
    changes_name = f'supporting/site_changes_{previous}_to_{version}.tsv'
    x56.site_changes(sites, x65.read(prev_changes)).to_csv(out / changes_name, sep='\t', index=False, lineterminator='\n')
    index = x68.site_index(out, version, config)

    # DepMap cell lines per changed residue
    models = x65.read(out / f'supporting/depmap_24q4_change_models_{version}.tsv').set_index('model_id')['cell_line']
    lines_at = {}
    for res in RESIDUES:
        t = x65.read(out / f'supporting/depmap_24q4_changes_{res}_{version}.tsv.gz')
        for (acc, pos), g in t.groupby(['canonical_UniProtAC', 'position'], sort=False):
            pairs = sorted({(models.get(m, ''), m) for m in g['model_id']})
            lines_at[(acc, int(pos), res)] = pairs
    am = x65.read(LAYERS / 'alphamissense_residues.tsv.gz')
    am_at = {(a, int(p), r): (mn, mx, n) for a, p, r, mn, mx, n in zip(am['canonical_UniProtAC'], am['position'], am['residue'],
                                                                     am['am_mean_pathogenicity'], am['am_max_pathogenicity'], am['am_pathogenic'])}
    files, ann_stats = {}, {}
    for name, spec in config['site_sets'].items():
        keys = sorted(k for k, v in index.items() if name in v.split(';'))
        rows = []
        for acc, pos, res in keys:
            pairs = lines_at.get((acc, pos, res), [])
            mn, mx, n = am_at.get((acc, pos, res), ('', '', ''))
            rows.append({'canonical_UniProtAC': acc, 'position': pos, 'residue': res, 'depmap_lines': len(pairs),
                         'depmap_cell_lines': '|'.join(c for c, _ in pairs), 'depmap_models': '|'.join(m for _, m in pairs),
                         'am_mean_pathogenicity': mn, 'am_max_pathogenicity': mx, 'am_pathogenic': n})
        t = pd.DataFrame(rows, columns=ANNOTATION_COLUMNS)
        fname = f'supporting/site_annotations_{name}_{version}.tsv'
        t.to_csv(out / fname, sep='\t', index=False, lineterminator='\n')
        files[fname] = {'rows': len(t), 'md5': md5_file(out / fname)}
        ann_stats[name] = {'label': x65.SET_LABELS.get(name, name), 'sites': len(t), 'sites_with_depmap_change': int((t['depmap_lines'] > 0).sum()),
                           'sites_with_alphamissense': int((t['am_mean_pathogenicity'] != '').sum()),
                           'sites_am_mean_above_0_564': int(pd.to_numeric(t['am_mean_pathogenicity'], errors='coerce').gt(0.564).sum())}

    # ClinVar disease families, appended
    cache, fam_stats, appended = {}, {}, {}
    for res in RESIDUES:
        fname = f'supporting/clinvar_changes_{res}_{version}.tsv.gz'
        t = x65.read(out / fname)
        t['disease_families'] = [disease_families(c, cache) for c in t['conditions']]
        P.write_gz(t, out / fname)
        appended[fname] = {'appended_columns': ['disease_families'], 'from_md5': carried[fname]['from_md5'], 'md5': md5_file(out / fname)}
        patho = t['germline_classification'].str.match(x68.PATHOGENIC)
        named = t['conditions'].map(lambda c: any(x.strip() and x.strip().lower() not in DROP_CONDITIONS for x in re.split(r'[|;]', c)))
        fam_stats[res] = {'rows': len(t), 'pathogenic_rows': int(patho.sum()), 'pathogenic_rows_with_named_condition': int((patho & named).sum()),
                          'pathogenic_rows_with_family': int((patho & (t['disease_families'] != '')).sum())}
    # the scope of the site-set test, appended
    fname = f'supporting/site_set_propensity_{version}.tsv'
    sets = x65.read(out / fname)
    tables = {src: x65.read(out / f'supporting/site_propensity_{src}_{version}.tsv.gz') for src in SOURCES}
    scope = []
    for _, r in sets.iterrows():
        t = tables[r['source']]
        sfx = '' if r['cohort'] == 'all' else '_no_hypermutated'
        sub = t[t['residue'] == r['residue']]
        in_set = sub['site_sets'].str.split(';').map(lambda v: r['site_set'] in v).to_numpy()
        s = P.set_scope(in_set, (sub['relative_mutability' + sfx] != '').to_numpy(), sub['k' + sfx].astype(int).to_numpy())
        if int(sub['k' + sfx].astype(int).to_numpy()[in_set & (sub['relative_mutability' + sfx] != '').to_numpy()].sum()) != int(r['observed']):
            raise SystemExit(f'site-set scope: observed differs for {r["site_set"]} {r["source"]} {r["cohort"]}')
        scope.append(s)
    for c in ('sites_total', 'sites_weighted', 'k_total'):
        sets[c] = [s[c] for s in scope]
    sets.to_csv(out / fname, sep='\t', index=False, lineterminator='\n')
    appended[fname] = {'appended_columns': ['sites_total', 'sites_weighted', 'k_total'], 'from_md5': carried[fname]['from_md5'], 'md5': md5_file(out / fname)}

    # data dictionary
    new_rows, after = dictionary_rows(version, config)
    lines = (prev / 'data_dictionary.tsv').read_text(encoding='utf-8').rstrip('\n').split('\n')
    dict_lines = [lines[0]]
    for line in lines[1:]:
        name, rest = line.split('\t', 1)
        name = changes_name if name.startswith('supporting/site_changes_') else x65.new_name(name, previous, version, arginine_sets)
        dict_lines.append(name + '\t' + rest)
    for f, extra in after.items():
        last = max(i for i, ln in enumerate(dict_lines) if ln.split('\t', 1)[0] == f)
        dict_lines[last + 1:last + 1] = extra
    dict_lines += new_rows
    (out / 'data_dictionary.tsv').write_text('\n'.join(dict_lines) + '\n', encoding='utf-8')

    write_docs(args, prev, out, version, previous, arginine_sets, carried, edited, changes_name, sites, files, appended, ann_stats, fam_stats, config)
    print(json.dumps({'files': files, 'annotations': ann_stats, 'clinvar_families': fam_stats}, indent=1))


def write_docs(args, prev, out, version, previous, arginine_sets, carried, edited, changes_name, sites, files, appended, ann_stats, fam_stats, config):
    rows = x65.count_marks(out, version, previous, config)
    prev_release = json.loads((prev / 'release.json').read_text(encoding='utf-8'))
    amsum = json.loads((LAYERS / 'alphamissense_summary.json').read_text(encoding='utf-8'))
    commit = code_commit()
    ann_table = ['| Site set | Sites | with a DepMap cell line changed | with AlphaMissense | mean of the 19 substitutions above 0.564 |', '|---|---:|---:|---:|---:|']
    for name, s in ann_stats.items():
        ann_table.append(f'| {s["label"]} | {s["sites"]:,} | {s["sites_with_depmap_change"]:,} | {s["sites_with_alphamissense"]:,} | {s["sites_am_mean_above_0_564"]:,} |')
    fam_named = sum(v['pathogenic_rows_with_named_condition'] for v in fam_stats.values())
    fam_mapped = sum(v['pathogenic_rows_with_family'] for v in fam_stats.values())
    c = amsum['counts']
    new_section = [
        '## New: per-site annotations, ClinVar disease families, and the scope of the site-set test', '',
        f'- **`supporting/site_annotations_<set>_{version}.tsv`** (methylarginine, methyllysine, lysine, phospho): one row per site of the '
        'set (accession, position, residue). Each row gives the DepMap 24Q4 cell lines with a missense or stop change at the site (count, '
        'names and model IDs), and AlphaMissense\'s mean and highest predicted pathogenicity over the residue\'s 19 substitutions, with how '
        'many of them it classes pathogenic. The site tables themselves keep their bytes (the methylarginine table is unchanged since v1.4.3, '
        'the context-mark tables since v1.5.0); join the annotations to them on accession, position and residue.',
        f'- **AlphaMissense** counts for an entry when its sequence equals the 2026_03 canonical sequence over its whole length: '
        f'{c.get("entries_matched", 0):,} entries do, {c.get("entries_sequence_differs", 0):,} differ and {c.get("entries_not_in_alphamissense", 0):,} are '
        f'missing (`inputs/layers/alphamissense_status.tsv`). {AM_NOTE}. 0.564 is AlphaMissense\'s threshold for calling a single '
        'substitution likely pathogenic (Cheng et al. 2023); a residue whose mean over all 19 lies above it is one where most substitutions are '
        'predicted damaging. Cell-line names and model IDs are separated by | (one name, RS4;11, holds a semicolon).',
        f'- **`disease_families`** is appended to `supporting/clinvar_changes_<residue>_{version}.tsv.gz`. Each variation\'s conditions are '
        'classed by `classify_clinvar()` in `scripts/disease_families.py`: ClinVar-specific rules first, then the DisGeNET rules of `classify()`, '
        'which is unchanged (so the disease analyses of earlier releases stay reproducible). On ClinVar names `classify()` alone misfiles '
        'common conditions (Lynch and Li-Fraumeni syndromes as developmental, retinal dystrophies as neurodegeneration, the developmental and '
        'epileptic encephalopathies as craniofacial, malignant hyperthermia as cancer); the ClinVar rules correct these and add a family, '
        f'Epilepsy and seizures. Of the {fam_named:,} pathogenic or likely pathogenic rows with a named condition, {fam_mapped:,} '
        f'({100 * fam_mapped / max(fam_named, 1):.1f}%) get a family.',
        f'- **`sites_total`, `sites_weighted`, `k_total`** are appended to `supporting/site_set_propensity_{version}.tsv`: the test covers '
        'sites_weighted of the sites_total sites and observed of the k_total changes (the rest lie on proteins without weights). Every '
        'earlier column keeps its values.', '',
        *ann_table, '',
        'The title of the release now names every mark it holds.', '']
    carried_text = re.sub(r'(^|\n)## New: (.)', lambda m: f'{m.group(1)}## {m.group(2).upper()}', x65.carried_readme(prev, previous, version, arginine_sets).rstrip('\n'))
    R = [f'# {TITLE}, {version}', '',
         f'{version} holds everything of {previous} under this release\'s names, its data unchanged except for columns appended to the '
         'ClinVar change tables and the site-set propensity table, and adds per-site annotations for every site set. See CHANGES.md.', '',
         *new_section, carried_text]
    (out / 'README.md').write_text(x65.americanize('\n'.join(R)) + '\n', encoding='utf-8')
    C = [f'# Changes from {previous} to {version}', '',
         f'Everything of {previous} is in {version} under this release\'s names. The owner\'s decisions of 23 November 2026: the per-site '
         'DepMap cell lines and AlphaMissense scores go in separate per-site tables (the site tables keep their bytes); ClinVar gets disease '
         'families; the site-set table gets the scope of its test; the title names every mark. AlphaMissense is a new source, so this is a '
         'minor release.', '',
         *new_section,
         '## New files', '', *[f'- `{f}`' for f in sorted(files)], '',
         '## Columns appended', '', *[f'- `{f}`: ' + ', '.join(f'`{c}`' for c in v['appended_columns']) for f, v in sorted(appended.items())], '',
         '## Files renamed', '',
         f'Every carried file is renamed `_{previous}` to `_{version}`, its content unchanged except where it names its version: '
         + ', '.join(f'`{e}`' for e in edited) + f', and the tables above with columns appended. `{changes_name}` replaces {previous}\'s site '
         f'changes table; it lists every one of the {len(sites):,} methylarginine sites as kept.', '',
         'Code: `scripts/72_stage_alphamissense.py` (AlphaMissense per residue), `scripts/disease_families.py` (`classify_clinvar`), '
         '`scripts/propensity.py` (`set_scope`), `scripts/73_release_site_annotations.py` (this release), `scripts/65_carry_release.py`, at '
         f'commit `{commit}`.']
    (out / 'CHANGES.md').write_text(x65.americanize('\n'.join(C)) + '\n', encoding='utf-8')

    def rename(obj):
        return json.loads(json.dumps(obj).replace(f'_{previous}.', f'_{version}.'))

    release = {k: rename(v) for k, v in prev_release.items() if k not in ('carried_over', 'carried_over_text_edited', 'site_changes', 'counts', 'corrections')}
    release.update({'title': TITLE, 'version': version, 'previous_version': previous, 'release_date': args.release_date, 'release_code_commit': commit,
                    'carried_over': carried, 'carried_over_text_edited': edited, 'appended_columns': appended, 'site_changes': changes_name, 'counts': rows,
                    'site_annotations': {'files': files, 'per_set': ann_stats, 'alphamissense': {'input': amsum['input'], 'license': amsum['license'],
                                                                                                 'citation': amsum['citation'], 'rule': amsum['rule'],
                                                                                                 'caveat': amsum['caveat'], 'counts': amsum['counts']}},
                    'clinvar_disease_families': {'rule': 'conditions split on | and ;, without not specified, not provided and see cases; classify_clinvar(); '
                                                         'Other / unassigned left out; sorted; ";"-joined', 'per_residue': fam_stats},
                    'corrections': []})
    (out / 'release.json').write_text(json.dumps(release, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    names = sorted(p.relative_to(out).as_posix() for p in out.rglob('*') if p.is_file() and p.name != 'MD5SUMS')
    (out / 'MD5SUMS').write_text(''.join(f'{md5_file(out / n)}  {n}\n' for n in names), encoding='utf-8')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Write a release that carries the previous release forward and adds the per-site mutation propensity of its change tables.

    python scripts/71_release_propensity.py --version v1.7.0 --previous v1.6.0 --release-date 2026-11-21

It carries every file of releases/<previous> under this release's names, as scripts/65_carry_release.py does (the data
unchanged), and adds to supporting/, from the carried change tables and the inputs staged by scripts/70_stage_propensity.py
(inputs/layers/propensity_*, gnomad_v2_1_1_mutation_rates.txt.bgz, *_changes_h3_all_residues.tsv.gz):

  site_propensity_<source>_<version>.tsv.gz   per source (tcga_mc3, depmap_24q4, clinvar, gnomad_v4_1): every candidate
                                              residue (R, K, S, T, Y) that is a known site of a site set, changed or not,
                                              and every other candidate residue with at least one change: observed and
                                              expected changes under the residue's codon mutability, O/E, binomial p and
                                              BH q; TCGA and DepMap also without the hypermutated samples
  site_frequency_by_group_<version>.tsv.gz    per residue changed and TCGA project or DepMap lineage: samples changed and
                                              samples profiled
  sample_burden_<version>.tsv                 per TCGA sample and DepMap line: coding calls, burden per Mb, hypermutated
  mutation_rates_<version>.tsv                the somatic spectra and the germline rates the mutability is computed from
  site_set_propensity_<version>.tsv           per site set, residue, source and cohort: changes at the set's sites within
                                              its proteins against expectation, with a sampled p
  propensity_controls_<version>.tsv           positive controls: the oncohistone residues of H3 (over every residue of H3)
                                              and recurrent cancer hotspots, pass or fail with the reason
  propensity_validation_<version>.tsv         the binomial p against a within-protein permutation, on 1,000 proteins

The method is in scripts/propensity.py. It then writes data_dictionary.tsv (the previous one under the new names, and the
new tables), README.md (what is new, then the previous release's README from its first section on), CHANGES.md,
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
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import propensity as P  # noqa: E402
from common import md5_file, parse_fasta  # noqa: E402


def _load(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, HERE / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


x65 = _load('carry65', '65_carry_release.py')
x68 = _load('release68', '68_release_cancer_changes.py')
LAYERS = ROOT / 'inputs' / 'layers'
RESIDUES = ['R', 'K', 'S', 'T', 'Y']
NAMES = x68.NAMES
SOURCES = {'tcga_mc3': 'TCGA MC3', 'depmap_24q4': 'DepMap 24Q4', 'clinvar': 'ClinVar', 'gnomad_v4_1': 'gnomAD v4.1'}
SOMATIC = ('tcga_mc3', 'depmap_24q4')
WHO = x68.WHO
GROUP = {'tcga_mc3': 'project_id', 'depmap_24q4': None}
COHORTS = {'tcga_mc3': ('all', 'no_hypermutated'), 'depmap_24q4': ('all', 'no_hypermutated'), 'clinvar': ('all',), 'gnomad_v4_1': ('all',)}
SUFFIX = {'all': '', 'no_hypermutated': '_no_hypermutated'}
OE_STATUS = {'tcga_mc3': 'test', 'depmap_24q4': 'test', 'gnomad_v4_1': 'test', 'clinvar': 'descriptive'}
GNOMAD_RATES = ('gnomAD mutation_rate_methylation_bins (gs://gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/model/'
                'mutation_rate_methylation_bins.txt.bgz; the gnomAD v2.1.1 flagship paper, Karczewski et al., Nature 581:434, 2020; CC0)')
CHANG = 'Chang et al., Nat Biotechnol 34:155 (2016), doi:10.1038/nbt.3391'
CAMPBELL = 'Campbell et al., Cell 171:1042 (2017), doi:10.1016/j.cell.2017.09.048'
# recurrent somatic hotspots at candidate residues (Chang et al. 2016, cancerhotspots.org), UniProt canonical numbering
HOTSPOTS = [('IDH1 R132', 'O75874', 132, 'R'), ('IDH2 R172', 'P48735', 172, 'R'), ('TP53 R175', 'P04637', 175, 'R'),
            ('TP53 R248', 'P04637', 248, 'R'), ('TP53 R273', 'P04637', 273, 'R'), ('PTEN R130', 'P60484', 130, 'R'),
            ('FBXW7 R465', 'Q969H0', 465, 'R'), ('SF3B1 K700', 'O75533', 700, 'K'), ('CTNNB1 S37', 'P35222', 37, 'S'),
            ('CTNNB1 S45', 'P35222', 45, 'S'), ('U2AF1 S34', 'Q01081', 34, 'S'), ('FGFR3 S249', 'P22607', 249, 'S'),
            ('ERBB2 S310', 'P04626', 310, 'S')]
# the oncohistone residues of H3 (histone numbering; UniProt counts the initiator methionine): label, UniProt position,
# residue, the named new residues
ONCOHISTONES = [('H3K27M', 28, 'K', {'M'}), ('H3K36M', 37, 'K', {'M'}), ('H3G34R/V', 35, 'G', {'R', 'V'})]
# the H3 entries the oncohistone mutations are found in: H3.1 (H3C1-H3C12), H3.2 (H3C13-H3C15), H3.3 (H3-3A, H3-3B)
H3_ENTRIES = ('P68431', 'Q71DI3', 'P84243')
VALIDATION_UNITS = 1000


def read_changes(out: Path, version: str, src: str, res: str) -> pd.DataFrame:
    cols = ['canonical_UniProtAC', 'position', WHO[src], 'gene'] + ([GROUP[src]] if GROUP.get(src) else [])
    t = pd.read_csv(out / f'supporting/{src}_changes_{res}_{version}.tsv.gz', sep='\t', dtype=str, keep_default_na=False, usecols=cols)
    t['position'] = t['position'].astype(int)
    return t


def k_counts(t: pd.DataFrame, who: str, drop: set | None = None) -> pd.Series:
    if drop:
        t = t[~t[who].isin(drop)]
    return t.groupby(['canonical_UniProtAC', 'position'])[who].nunique()


def candidates(canonical: dict, res: str) -> pd.DataFrame:
    accs, pos = [], []
    for acc in sorted(canonical):
        seq = canonical[acc]
        ps = [i + 1 for i, c in enumerate(seq) if c == res]
        accs += [acc] * len(ps)
        pos += ps
    return pd.DataFrame({'canonical_UniProtAC': accs, 'position': np.array(pos, dtype=np.int64)})


def w_column(frame: pd.DataFrame, mut: dict, key: str) -> np.ndarray:
    w = np.full(len(frame), np.nan)
    for acc, idx in frame.groupby('canonical_UniProtAC', sort=False).indices.items():
        m = mut.get(acc)
        if m is not None:
            w[idx] = m[key][frame['position'].to_numpy()[idx] - 1]
    return w


def spectrum_key(src: str, cohort: str) -> str:
    return 'germline' if src not in SOMATIC else f'{src}/{cohort}'


def code_commit() -> str:
    paths = ['scripts/71_release_propensity.py', 'scripts/propensity.py', 'scripts/70_stage_propensity.py', 'scripts/65_carry_release.py',
             'scripts/68_release_cancer_changes.py', 'config/site_sets.json']
    try:
        dirty = subprocess.run(['git', 'status', '--porcelain', '--', *paths], capture_output=True, text=True, check=True, cwd=ROOT).stdout
        sha = subprocess.run(['git', 'log', '-1', '--format=%H', '--', *paths], capture_output=True, text=True, check=True, cwd=ROOT).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return 'unknown'
    return sha + (' (scripts modified)' if dirty.strip() else '')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--version', required=True)
    ap.add_argument('--previous', required=True)
    ap.add_argument('--release-date', required=True)
    ap.add_argument('--config', default=str(ROOT / 'config' / 'site_sets.json'))
    ap.add_argument('--fasta', default=str(ROOT / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03.fasta'))
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

    # inputs
    canonical = parse_fasta(args.fasta)
    index = x68.site_index(out, version, config)
    cds = pd.read_csv(LAYERS / 'propensity_cds.tsv.gz', sep='\t', dtype=str, keep_default_na=False)
    spectra = pd.read_csv(LAYERS / 'propensity_somatic_spectra.tsv', sep='\t', dtype=str, keep_default_na=False)
    somatic = P.somatic_rates(spectra)
    germline = P.germline_rates(LAYERS / 'gnomad_v2_1_1_mutation_rates.txt.bgz')
    mut = P.mutability(cds, somatic, germline)
    burden = pd.read_csv(LAYERS / 'propensity_sample_burden.tsv', sep='\t', dtype=str, keep_default_na=False)
    hyper = {src: set(burden.loc[(burden['source'] == src) & (burden['hypermutated'] == 'True'), 'sample_id']) for src in SOMATIC}
    profiled = {src: {'all': int((burden['source'] == src).sum())} for src in SOMATIC}
    for src in SOMATIC:
        profiled[src]['no_hypermutated'] = profiled[src]['all'] - len(hyper[src])
    staging = json.loads((LAYERS / 'propensity_inputs_summary.json').read_text(encoding='utf-8'))
    gene_of = defaultdict(set)
    for acc, g in zip(cds['canonical_UniProtAC'], cds['gene']):
        if g:
            gene_of[acc].add(g)

    # per-site tables
    files, stats, keep = {}, {}, {}
    found = {}           # (source, cohort, accession, position) -> the test of that residue, for the controls
    family = {}          # (source, cohort, residue) -> residues tested
    for src in SOURCES:
        parts = []
        for res in RESIDUES:
            t = read_changes(out, version, src, res)
            for acc, g in zip(t['canonical_UniProtAC'], t['gene']):
                if g and acc not in gene_of:
                    gene_of.setdefault(('fallback', acc), set()).add(g)
            frame = candidates(canonical, res)
            frame['residue'] = res
            frame['site_sets'] = [index.get((a, int(p), res), '') for a, p in zip(frame['canonical_UniProtAC'], frame['position'])]
            for cohort in COHORTS[src]:
                sfx = SUFFIX[cohort]
                k = k_counts(t, WHO[src], hyper[src] if cohort == 'no_hypermutated' else None)
                frame['k' + sfx] = pd.Series(list(zip(frame['canonical_UniProtAC'], frame['position']))).map(k).fillna(0).astype(int).to_numpy()
                frame['_w' + sfx] = w_column(frame, mut, spectrum_key(src, cohort))
                tested = P.site_test(frame, '_w' + sfx, 'k' + sfx, lower=(src == 'gnomad_v4_1'))
                for c in tested.columns:
                    frame[c + sfx] = tested[c]
                if src in SOMATIC:
                    frame['samples_profiled' + sfx] = profiled[src][cohort]
                ok = frame['p_upper' + sfx].notna()
                family[(src, cohort, res)] = int(ok.sum())
                st = stats.setdefault(f'{src}/{res}', {})
                st[cohort] = {'candidate_residues': len(frame), 'tested': int(ok.sum()), 'with_change': int((frame['k' + sfx] > 0).sum()),
                              'q_upper_below_0_05': int((frame['q_upper' + sfx] < P.ALPHA).sum()),
                              'known_sites_q_upper_below_0_05': int(((frame['q_upper' + sfx] < P.ALPHA) & (frame['site_sets'] != '')).sum())}
                if src == 'gnomad_v4_1':
                    st[cohort]['q_lower_below_0_05'] = int((frame['q_lower'] < P.ALPHA).sum())
            keep[(src, res)] = frame
            parts.append(frame[(frame['site_sets'] != '') | (frame['k'] > 0)])
            stats[f'{src}/{res}']['rows'] = int(parts[-1].shape[0])
            for cohort in COHORTS[src]:
                sfx = SUFFIX[cohort]
                sub = frame[frame['p_upper' + sfx].notna()]
                for c in HOTSPOT_KEYS(res):
                    hit = sub[(sub['canonical_UniProtAC'] == c[0]) & (sub['position'] == c[1])]
                    if len(hit):
                        r = hit.iloc[0]
                        found[(src, cohort, c[0], c[1])] = {'k': int(r['k' + sfx]), 'n': int(r['n' + sfx]), 'expected': float(r['expected' + sfx]),
                                                            'oe': float(r['oe' + sfx]), 'p': float(r['p_upper' + sfx]), 'q': float(r['q_upper' + sfx]),
                                                            'pi': float(r['expected' + sfx]) / float(r['n' + sfx]), 'm': family[(src, cohort, res)]}
        table = pd.concat(parts, ignore_index=True)
        table['gene'] = [';'.join(sorted(gene_of.get(a) or gene_of.get(('fallback', a)) or [])) for a in table['canonical_UniProtAC']]
        table['oe_status'] = OE_STATUS[src]
        cols = ['canonical_UniProtAC', 'position', 'residue', 'gene', 'site_sets', 'residues_of_type']
        for cohort in COHORTS[src]:
            sfx = SUFFIX[cohort]
            cols += [c + sfx for c in ('relative_mutability', 'k', 'n', 'expected', 'oe', 'expected_unweighted', 'oe_unweighted', 'p_upper', 'q_upper')]
            if src == 'gnomad_v4_1':
                cols += ['p_lower', 'q_lower']
            if src in SOMATIC:
                cols.append('samples_profiled' + sfx)
        cols.append('oe_status')
        table = table[cols].copy()
        for c in cols:
            if c.startswith(('relative_mutability', 'expected', 'oe', 'p_', 'q_')) and c != 'oe_status':
                digits = 4 if c.startswith(('p_', 'q_')) else 6
                table[c] = [P.fmt(x, digits) for x in table[c].to_numpy(dtype=float)]
        table = table.sort_values(['canonical_UniProtAC', 'position'], kind='stable')
        name = f'supporting/site_propensity_{src}_{version}.tsv.gz'
        P.write_gz(table, out / name)
        files[name] = {'rows': len(table), 'md5': md5_file(out / name), 'bytes': (out / name).stat().st_size}

    # frequency per TCGA project and DepMap lineage
    lineage = dict(zip(burden.loc[burden['source'] == 'depmap_24q4', 'sample_id'], burden.loc[burden['source'] == 'depmap_24q4', 'group']))
    denom = {src: burden[burden['source'] == src].groupby('group')['sample_id'].size().to_dict() for src in SOMATIC}
    denom_nh = {src: burden[(burden['source'] == src) & (burden['hypermutated'] != 'True')].groupby('group')['sample_id'].size().to_dict() for src in SOMATIC}
    freq = []
    for src in SOMATIC:
        who = WHO[src]
        for res in RESIDUES:
            t = read_changes(out, version, src, res)
            t['group'] = t['project_id'] if src == 'tcga_mc3' else t[who].map(lineage).fillna('')
            t['_hyper'] = t[who].isin(hyper[src])
            g = t.groupby(['canonical_UniProtAC', 'position', 'group'])
            f = g[who].nunique().rename('samples_mutated').to_frame()
            f['samples_mutated_no_hypermutated'] = t[~t['_hyper']].groupby(['canonical_UniProtAC', 'position', 'group'])[who].nunique()
            f = f.reset_index()
            f['samples_mutated_no_hypermutated'] = f['samples_mutated_no_hypermutated'].fillna(0).astype(int)
            f['residue'], f['source'] = res, src
            f['samples_profiled'] = f['group'].map(denom[src]).fillna(0).astype(int)
            f['samples_profiled_no_hypermutated'] = f['group'].map(denom_nh[src]).fillna(0).astype(int)
            freq.append(f)
    freq = pd.concat(freq, ignore_index=True)[['canonical_UniProtAC', 'position', 'residue', 'source', 'group', 'samples_mutated', 'samples_profiled',
                                               'samples_mutated_no_hypermutated', 'samples_profiled_no_hypermutated']]
    freq = freq.sort_values(['canonical_UniProtAC', 'position', 'source', 'group'], kind='stable')
    name = f'supporting/site_frequency_by_group_{version}.tsv.gz'
    P.write_gz(freq, out / name)
    files[name] = {'rows': len(freq), 'md5': md5_file(out / name), 'bytes': (out / name).stat().st_size}
    # every group's denominators, those without a changed sample included (for k = 0 per group)
    den = pd.DataFrame([{'source': src, 'group': g, 'samples_profiled': n, 'samples_profiled_no_hypermutated': int(denom_nh[src].get(g, 0))}
                        for src in SOMATIC for g, n in sorted(denom[src].items())])
    name = f'supporting/group_denominators_{version}.tsv'
    den.to_csv(out / name, sep='\t', index=False, lineterminator='\n')
    files[name] = {'rows': len(den), 'md5': md5_file(out / name)}

    # burden and rates
    name = f'supporting/sample_burden_{version}.tsv'
    shutil.copyfile(LAYERS / 'propensity_sample_burden.tsv', out / name)
    files[name] = {'rows': len(burden), 'md5': md5_file(out / name), 'from': 'inputs/layers/propensity_sample_burden.tsv'}
    rates = spectra.assign(methylation_level='')[['source', 'cohort', 'context', 'ref', 'alt', 'methylation_level', 'observed', 'opportunity', 'rate']]
    g_rows = []
    for k, (ctx, ref, alt) in enumerate(P.CLASSES):
        for lv in ((0, 1, 2) if P.CPG_TRANSITION[k] else (0,)):
            g_rows.append({'source': 'gnomad_germline', 'cohort': '', 'context': ctx, 'ref': ref, 'alt': alt, 'methylation_level': str(lv),
                           'observed': '', 'opportunity': '', 'rate': f'{germline[k, lv]:.4e}'})
    rates = pd.concat([rates, pd.DataFrame(g_rows)], ignore_index=True)
    name = f'supporting/mutation_rates_{version}.tsv'
    rates.to_csv(out / name, sep='\t', index=False, lineterminator='\n')
    files[name] = {'rows': len(rates), 'md5': md5_file(out / name)}

    # site sets
    set_rows = []
    i = 0
    for set_name, spec in config['site_sets'].items():
        for res in sorted(spec['residues']):
            for src in SOURCES:
                frame = keep[(src, res)]
                in_set = frame['site_sets'].str.split(';').map(lambda s: set_name in s).to_numpy()
                for cohort in COHORTS[src]:
                    sfx = SUFFIX[cohort]
                    tested = frame['p_upper' + sfx].notna().to_numpy()
                    sub = frame[tested]
                    flag = in_set[tested]
                    units, unw = [], []
                    n_sites = 0
                    for acc, idx in sub.groupby('canonical_UniProtAC', sort=False).indices.items():
                        f = flag[idx]
                        if not f.any():
                            continue
                        w = sub['_w' + sfx].to_numpy()[idx]
                        k = sub['k' + sfx].to_numpy()[idx]
                        n = int(k.sum())
                        units.append((n, float(w[f].sum() / w.sum()), int(k[f].sum())))
                        unw.append(n * f.sum() / len(f))
                        n_sites += int(f.sum())
                    rng = np.random.default_rng(P.SEED + i)
                    r = P.set_test(units, rng) if units else {'observed': 0, 'expected': np.nan, 'oe': np.nan, 'p_upper': np.nan, 'p_lower': np.nan}
                    eu = float(sum(unw)) if units else np.nan
                    set_rows.append({'site_set': set_name, 'label': x65.SET_LABELS.get(set_name, set_name), 'residue': res, 'source': src, 'cohort': cohort,
                                     'proteins': len(units), 'sites': n_sites,
                                     'residues_in_proteins': int(sum(len(idx) for acc, idx in sub.groupby('canonical_UniProtAC', sort=False).indices.items()
                                                                     if flag[idx].any())),
                                     'observed': r['observed'], 'expected': P.fmt(r['expected']), 'oe': P.fmt(r['oe']), 'expected_unweighted': P.fmt(eu),
                                     'oe_unweighted': P.fmt(r['observed'] / eu if eu and eu > 0 else np.nan), 'p_upper': P.fmt(r['p_upper'], 4),
                                     'p_lower': P.fmt(r['p_lower'], 4), 'iterations': P.ITERATIONS, 'seed': P.SEED + i, 'oe_status': OE_STATUS[src],
                                     **P.set_scope(in_set, frame['relative_mutability' + sfx].notna().to_numpy(), frame['k' + sfx].to_numpy())})
                    i += 1
    set_table = pd.DataFrame(set_rows)
    name = f'supporting/site_set_propensity_{version}.tsv'
    set_table.to_csv(out / name, sep='\t', index=False, lineterminator='\n')
    files[name] = {'rows': len(set_table), 'md5': md5_file(out / name)}

    # controls
    controls = control_rows(canonical, mut, hyper, staging, found, gene_of)
    name = f'supporting/propensity_controls_{version}.tsv'
    controls.to_csv(out / name, sep='\t', index=False, lineterminator='\n')
    files[name] = {'rows': len(controls), 'md5': md5_file(out / name)}

    # permutation check
    val_rows = []
    for j, src in enumerate(SOURCES):
        units = []
        for res in RESIDUES:
            frame = keep[(src, res)]
            sub = frame[frame['p_upper'].notna()]
            for acc, idx in sub.groupby('canonical_UniProtAC', sort=False).indices.items():
                k = sub['k'].to_numpy()[idx]
                if k.sum() >= 2:
                    units.append((res, acc, k.astype(np.int64), sub['_w'].to_numpy()[idx]))
        rng = np.random.default_rng(P.SEED + 1000 + j)
        pick = rng.choice(len(units), size=min(VALIDATION_UNITS, len(units)), replace=False)
        pick.sort()
        r = P.permutation_check([(units[x][2], units[x][3]) for x in pick], rng)
        val_rows.append({'source': src, 'cohort': 'all', 'units_available': len(units), **{k: (P.fmt(v, 4) if isinstance(v, float) else v) for k, v in r.items()},
                         'iterations': P.ITERATIONS, 'seed': P.SEED + 1000 + j})
    validation = pd.DataFrame(val_rows)
    name = f'supporting/propensity_validation_{version}.tsv'
    validation.to_csv(out / name, sep='\t', index=False, lineterminator='\n')
    files[name] = {'rows': len(validation), 'md5': md5_file(out / name)}

    write_docs(args, prev, out, version, previous, arginine_sets, carried, edited, changes_name, sites, files, stats, set_table, controls,
               validation, staging, burden, hyper, profiled, cds, mut, canonical)
    print(json.dumps({'files': files, 'stats': stats}, indent=1))


def HOTSPOT_KEYS(res: str) -> list[tuple[str, int]]:
    out = [(acc, pos) for _, acc, pos, r in HOTSPOTS if r == res]
    if res == 'K':
        out += [(acc, pos) for acc in H3_ENTRIES for _, pos, r, _ in ONCOHISTONES if r == 'K']
    return out


def control_rows(canonical: dict, mut: dict, hyper: dict, staging: dict, found: dict, gene_of: dict) -> pd.DataFrame:
    rows = []

    def result(q: float, k: int, need) -> tuple[str, str]:
        if q == q and q < P.ALPHA:
            return 'pass', f'q < {P.ALPHA}'
        if need is None:
            return 'fail', 'not testable'
        return 'fail', f'{k} change{"s" if k != 1 else ""} at the residue; at least {need} would pass at this n (p <= {P.ALPHA} / residues tested)'

    # oncohistones, over every residue of each H3 entry
    tables = {'tcga_mc3': pd.read_csv(LAYERS / 'tcga_mc3_changes_h3_all_residues.tsv.gz', sep='\t', dtype=str, keep_default_na=False),
              'depmap_24q4': pd.read_csv(LAYERS / 'depmap_24q4_changes_h3_all_residues.tsv.gz', sep='\t', dtype=str, keep_default_na=False)}
    for acc in H3_ENTRIES:
        if acc not in staging['h3_entries']:
            raise SystemExit(f'{acc} is not among the staged H3 entries')
        seq = canonical[acc]
        hits = [(lab, pos, r, alts) for lab, pos, r, alts in ONCOHISTONES if len(seq) >= pos and seq[pos - 1] == r]
        if not hits:
            continue
        for src, t in tables.items():
            t = t[t['canonical_UniProtAC'] == acc].copy()
            t['position'] = t['position'].astype(int)
            for cohort in ('all', 'no_hypermutated'):
                tt = t[~t[WHO[src]].isin(hyper[src])] if cohort == 'no_hypermutated' else t
                k_all = tt.groupby('position')[WHO[src]].nunique()
                frame = pd.DataFrame({'canonical_UniProtAC': acc, 'position': np.arange(1, len(seq) + 1)})
                frame['k'] = frame['position'].map(k_all).fillna(0).astype(int)
                m = mut.get(acc)
                frame['w'] = m[spectrum_key(src, cohort)] if m is not None else np.nan
                res = P.site_test(frame, 'w', 'k')
                m_tested = int(res['p_upper'].notna().sum())
                for lab, pos, r, alts in hits:
                    i = pos - 1
                    k, n = int(frame['k'].iloc[i]), int(res['n'].iloc[i])
                    named = int(tt[(tt['position'] == pos) & tt['alt'].isin(alts)][WHO[src]].nunique())
                    pi = float(res['expected'].iloc[i]) / n if n else float('nan')
                    q = float(res['q_upper'].iloc[i])
                    need = P.changes_needed(k, n, pi, m_tested) if n else None
                    outcome, why = result(q, k, need) if n else ('not_tested', 'no change at any residue of the entry in this source and cohort')
                    rows.append({'control': f'{lab} (UniProt {r}{pos})', 'kind': 'oncohistone', 'canonical_UniProtAC': acc,
                                 'gene': ';'.join(sorted(gene_of.get(acc, []))), 'position': pos, 'residue': r,
                                 'named_change': '/'.join(f'{r}{pos}{a}' for a in sorted(alts)), 'source': src, 'cohort': cohort,
                                 'family': 'every residue of the entry', 'residues_tested': m_tested, 'k': k, 'k_named_change': named, 'n': n,
                                 'expected': P.fmt(float(res['expected'].iloc[i])), 'oe': P.fmt(float(res['oe'].iloc[i])),
                                 'p_upper': P.fmt(float(res['p_upper'].iloc[i]), 4), 'q_upper': P.fmt(q, 4), 'changes_needed': '' if need is None else need,
                                 'result': outcome, 'reason': why})
    # the per-site tables' own result at the hotspots and the H3 lysines
    for src in SOMATIC:
        for cohort in ('all', 'no_hypermutated'):
            for lab, acc, pos, r in HOTSPOTS + [(f'{lab} (UniProt {r}{pos})', acc, pos, r) for acc in H3_ENTRIES
                                                 for lab, pos, r, _ in ONCOHISTONES if r == 'K']:
                seq = canonical.get(acc, '')
                if len(seq) < pos or seq[pos - 1] != r:
                    rows.append({'control': lab, 'kind': 'recurrent hotspot' if 'UniProt' not in lab else 'oncohistone', 'canonical_UniProtAC': acc,
                                 'gene': ';'.join(sorted(gene_of.get(acc, []))), 'position': pos, 'residue': r, 'named_change': '', 'source': src,
                                 'cohort': cohort, 'family': '', 'result': 'not_tested', 'reason': 'the canonical sequence has another residue there'})
                    continue
                f = found.get((src, cohort, acc, pos))
                if f is None:
                    rows.append({'control': lab, 'kind': 'recurrent hotspot' if 'UniProt' not in lab else 'oncohistone', 'canonical_UniProtAC': acc,
                                 'gene': ';'.join(sorted(gene_of.get(acc, []))), 'position': pos, 'residue': r, 'named_change': '', 'source': src,
                                 'cohort': cohort, 'family': f'{NAMES[r]}s of every protein ({r})', 'result': 'fail',
                                 'reason': 'not tested: no change on the protein\'s residues of this type, or no coding sequence'})
                    continue
                need = P.changes_needed(f['k'], f['n'], f['pi'], f['m'])
                outcome, why = result(f['q'], f['k'], need)
                rows.append({'control': lab, 'kind': 'recurrent hotspot' if 'UniProt' not in lab else 'oncohistone', 'canonical_UniProtAC': acc,
                             'gene': ';'.join(sorted(gene_of.get(acc, []))), 'position': pos, 'residue': r, 'named_change': '', 'source': src,
                             'cohort': cohort, 'family': f'{NAMES[r]}s of every protein ({r})', 'residues_tested': f['m'], 'k': f['k'],
                             'k_named_change': '', 'n': f['n'], 'expected': P.fmt(f['expected']), 'oe': P.fmt(f['oe']), 'p_upper': P.fmt(f['p'], 4),
                             'q_upper': P.fmt(f['q'], 4), 'changes_needed': '' if need is None else need, 'result': outcome, 'reason': why})
    cols = ['control', 'kind', 'canonical_UniProtAC', 'gene', 'position', 'residue', 'named_change', 'source', 'cohort', 'family', 'residues_tested',
            'k', 'k_named_change', 'n', 'expected', 'oe', 'p_upper', 'q_upper', 'changes_needed', 'result', 'reason']
    out = pd.DataFrame(rows).reindex(columns=cols)
    for c in ('residues_tested', 'k', 'k_named_change', 'n', 'changes_needed'):
        out[c] = [str(int(x)) if isinstance(x, (int, float, np.integer, np.floating)) and x == x else (x if isinstance(x, str) else '') for x in out[c]]
    return out.fillna('')


def dictionary_rows(version: str) -> list[str]:
    ref = 'reviewed UniProtKB accession (canonical sequence, UniProt 2026_03)'
    unit = {'tcga_mc3': 'tumor samples', 'depmap_24q4': 'cell lines', 'clinvar': 'ClinVar variations', 'gnomad_v4_1': 'gnomAD variants'}
    rows = []
    for src in SOURCES:
        f = f'supporting/site_propensity_{src}_{version}.tsv.gz'
        spec = 'its own spectrum (mutation_rates, source ' + src + ')' if src in SOMATIC else 'gnomAD\'s germline rates (mutation_rates, source gnomad_germline)'
        cols = [('canonical_UniProtAC', 'string', '', '', ref), ('position', 'integer', '', '', 'position in the canonical sequence'),
                ('residue', 'string', 'R; K; S; T; Y', '', 'the residue'),
                ('gene', 'string', '', '', 'HGNC approved symbol of the gene(s) whose coding sequences give the weights, separated by semicolons '
                 '(an entry without one: the genes of its changes)'),
                ('site_sets', 'string', '', '', 'the site sets of this release with a site at this accession, position and residue '
                 '(config/site_sets.json), separated by semicolons; empty: a candidate residue, not a known site'),
                ('residues_of_type', 'integer', '', '', 'residues of this type in the canonical sequence')]
        for cohort in COHORTS[src]:
            sfx = SUFFIX[cohort]
            who = unit[src] + (' (hypermutated ones left out: more than 10 coding calls per Mb, sample_burden)' if sfx else '')
            cols += [('relative_mutability' + sfx, 'number', '', '', f'the residue\'s codon mutability under {spec} over the mean of the '
                      'protein\'s residues of this type (1 = average); empty: the entry has no coding sequence that translates to its canonical sequence'),
                     ('k' + sfx, 'integer', '', '', f'{who} with a missense or stop change at the residue (the change tables)'),
                     ('n' + sfx, 'integer', '', '', 'the sum of k over the protein\'s residues of this type'),
                     ('expected' + sfx, 'number', '', '', 'n times the residue\'s share of the mutability of the protein\'s residues of this type'),
                     ('oe' + sfx, 'number', '', '', 'k / expected' + (' (descriptive: ClinVar records clinical testing and submission, not '
                      'mutation; oe_status)' if src == 'clinvar' else '')),
                     ('expected_unweighted' + sfx, 'number', '', '', 'n / residues_of_type (every residue of the type equally mutable)'),
                     ('oe_unweighted' + sfx, 'number', '', '', 'k / expected_unweighted'),
                     ('p_upper' + sfx, 'number', '', '', 'P(X >= k), X ~ Binomial(n, the residue\'s share of the mutability); empty where not tested '
                      '(n = 0, or no coding sequence)'),
                     ('q_upper' + sfx, 'number', '', '', 'Benjamini-Hochberg q of p_upper over every residue tested of this source, cohort and residue type')]
            if src == 'gnomad_v4_1':
                cols += [('p_lower', 'number', '', '', 'P(X <= k): fewer variants than expected (depletion)'),
                         ('q_lower', 'number', '', '', 'Benjamini-Hochberg q of p_lower, as q_upper')]
            if src in SOMATIC:
                cols.append(('samples_profiled' + sfx, 'integer', '', '', f'{unit[src]} profiled' + (' that are not hypermutated' if sfx else '')))
        cols.append(('oe_status', 'string', 'test; descriptive', 'test = a test of mutation propensity; descriptive = ClinVar, whose '
                     'variants reflect clinical testing and submission more than mutation, so its O/E describes and does not test',
                     'how to read oe and p'))
        rows += ['\t'.join([f, *c]) for c in cols]
    f = f'supporting/site_frequency_by_group_{version}.tsv.gz'
    for c in (('canonical_UniProtAC', 'string', '', '', ref), ('position', 'integer', '', '', 'position in the canonical sequence'),
              ('residue', 'string', 'R; K; S; T; Y', '', 'the residue'), ('source', 'string', 'tcga_mc3; depmap_24q4', '', 'the change table'),
              ('group', 'string', '', '', 'TCGA project, or the cell line\'s Oncotree lineage (DepMap Model.csv)'),
              ('samples_mutated', 'integer', '', '', 'tumor samples or cell lines of the group with a change at the residue'),
              ('samples_profiled', 'integer', '', '', 'tumor samples or cell lines of the group profiled'),
              ('samples_mutated_no_hypermutated', 'integer', '', '', 'samples_mutated, the hypermutated ones left out'),
              ('samples_profiled_no_hypermutated', 'integer', '', '', 'samples_profiled, the hypermutated ones left out')):
        rows.append('\t'.join([f, *c]))
    f = f'supporting/group_denominators_{version}.tsv'
    for c in (('source', 'string', 'tcga_mc3; depmap_24q4', '', 'TCGA MC3 or DepMap 24Q4'),
              ('group', 'string', '', '', 'TCGA project, or Oncotree lineage (empty: a DepMap line without one in Model.csv)'),
              ('samples_profiled', 'integer', '', '', 'tumor samples or cell lines of the group profiled (sample_burden), changed or not'),
              ('samples_profiled_no_hypermutated', 'integer', '', '', 'of those, not hypermutated')):
        rows.append('\t'.join([f, *c]))
    f = f'supporting/sample_burden_{version}.tsv'
    for c in (('source', 'string', 'tcga_mc3; depmap_24q4', '', 'TCGA MC3 or DepMap 24Q4'),
              ('sample_id', 'string', '', '', 'TCGA tumor sample (barcode to the sample, 15 characters) or DepMap model ID'),
              ('group', 'string', '', '', 'TCGA project or Oncotree lineage'),
              ('coding_calls', 'integer', '', '', 'TCGA: kept MC3 calls of a coding Variant_Classification (Missense_Mutation, Nonsense_Mutation, '
               'Nonstop_Mutation, Translation_Start_Site, Silent, Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins); DepMap: calls with '
               'a coding consequence (OmicsSomaticMutations.csv keeps synonymous calls only where they touch a splice region)'),
              ('non_silent_calls', 'integer', '', '', 'coding_calls less silent'), ('missense', 'integer', '', '', 'missense calls'),
              ('stop', 'integer', '', '', 'stop-gained calls'), ('silent', 'integer', '', '', 'silent calls (TCGA only)'),
              ('indel', 'integer', '', '', 'frameshift and in-frame insertions and deletions'),
              ('other', 'integer', '', '', 'stop lost, start lost, other protein-altering calls'),
              ('snvs_in_coding_sequences', 'integer', '', '', 'residue-changing SNVs in the coding sequences of the weights (the spectrum\'s count)'),
              ('burden_per_mb', 'number', '', '', 'coding_calls per Mb of coding sequence (32.3 Mb: the coding sequences of the weights, each position once)'),
              ('hypermutated', 'boolean', 'True; False', '', f'burden_per_mb above threshold_per_mb ({CAMPBELL})'),
              ('threshold_per_mb', 'number', '', '', 'the threshold')):
        rows.append('\t'.join([f, *c]))
    f = f'supporting/mutation_rates_{version}.tsv'
    for c in (('source', 'string', 'tcga_mc3; depmap_24q4; gnomad_germline', '', 'the spectrum: TCGA and DepMap from their own calls; gnomad_germline, '
               'gnomAD\'s mutation_rate_methylation_bins, the rates of gnomAD and ClinVar'),
              ('cohort', 'string', 'all; no_hypermutated', '', 'the samples of a somatic spectrum (empty for the germline rates)'),
              ('context', 'string', '', '', 'trinucleotide, pyrimidine-centered (the strand with C or T in the middle)'),
              ('ref', 'string', 'C; T', '', 'the middle base'), ('alt', 'string', 'A; C; G; T', '', 'the new base'),
              ('methylation_level', 'integer', '0; 1; 2', '0 = CpG methylation mean at most 0.2 or not CpG; 1 = above 0.2; 2 = above 0.6',
               'germline rates only; levels 1 and 2 only for CpG transitions (gnomAD constraint code)'),
              ('observed', 'integer', '', '', 'somatic: residue-changing SNVs of the class in the coding sequences'),
              ('opportunity', 'integer', '', '', 'somatic: positions where the class would change the residue (each genomic position once)'),
              ('rate', 'number', '', '', 'somatic: observed / opportunity (relative); germline: gnomAD\'s mu_snp per generation')):
        rows.append('\t'.join([f, *c]))
    f = f'supporting/site_set_propensity_{version}.tsv'
    for c in (('site_set', 'string', '', '', 'site set (config/site_sets.json)'), ('label', 'string', '', '', 'its name'),
              ('residue', 'string', 'R; K; S; T; Y', '', 'residue type'), ('source', 'string', '', '', 'the change table'),
              ('cohort', 'string', 'all; no_hypermutated', '', 'the samples (somatic sources)'),
              ('proteins', 'integer', '', '', 'proteins tested with at least one site of the set'),
              ('sites', 'integer', '', '', 'the set\'s sites on them'), ('residues_in_proteins', 'integer', '', '', 'residues of the type on them'),
              ('observed', 'integer', '', '', 'the sum of k over the set\'s sites'),
              ('expected', 'number', '', '', 'the sum over the proteins of n times the sites\' share of the protein\'s mutability'),
              ('oe', 'number', '', '', 'observed / expected'), ('expected_unweighted', 'number', '', '', 'as expected, every residue equally mutable'),
              ('oe_unweighted', 'number', '', '', 'observed / expected_unweighted'),
              ('p_upper', 'number', '', '', 'share of draws, each protein\'s n changes redrawn over its residues by mutability, with at least the '
               'observed at the sites ((1 + count) / (iterations + 1))'),
              ('p_lower', 'number', '', '', 'the same, at most the observed'), ('iterations', 'integer', '', '', 'draws'),
              ('seed', 'integer', '', '', 'NumPy default_rng seed'), ('oe_status', 'string', 'test; descriptive', '', 'as in the per-site tables'),
              ('sites_total', 'integer', '', '', 'the set\'s sites of this residue type, on every protein'),
              ('sites_weighted', 'integer', '', '', 'of those, the sites on proteins with weights (a coding sequence that translates to the canonical '
               'sequence): the sites observed and expected cover (sites counts those of them on proteins with a change)'),
              ('k_total', 'integer', '', '', 'the sum of k over all sites_total, proteins without weights included (observed is the sum over '
               'sites_weighted)')):
        rows.append('\t'.join([f, *c]))
    f = f'supporting/propensity_controls_{version}.tsv'
    for c in (('control', 'string', '', '', 'the control (histone numbering for the oncohistones, with the UniProt position)'),
              ('kind', 'string', 'oncohistone; recurrent hotspot', '', f'oncohistone residue of H3, or a recurrent cancer hotspot ({CHANG})'),
              ('canonical_UniProtAC', 'string', '', '', ref), ('gene', 'string', '', '', 'HGNC approved symbol(s)'),
              ('position', 'integer', '', '', 'UniProt position'), ('residue', 'string', '', '', 'the residue'),
              ('named_change', 'string', '', '', 'the oncohistone change(s)'), ('source', 'string', 'tcga_mc3; depmap_24q4', '', 'the change table'),
              ('cohort', 'string', 'all; no_hypermutated', '', 'the samples'),
              ('family', 'string', '', '', 'the residues tested together: every residue of the H3 entry (the oncohistone rows from '
               'the H3 tables at every residue), or the residues of the type in every protein (the per-site tables)'),
              ('residues_tested', 'integer', '', '', 'residues in the family (the BH family)'), ('k', 'integer', '', '', 'samples or lines changed at the residue'),
              ('k_named_change', 'integer', '', '', 'of those, with the named change'), ('n', 'integer', '', '', 'changes at the protein\'s residues of the family'),
              ('expected', 'number', '', '', 'expected changes'), ('oe', 'number', '', '', 'k / expected'), ('p_upper', 'number', '', '', 'binomial p'),
              ('q_upper', 'number', '', '', 'BH q within the family'),
              ('changes_needed', 'integer', '', '', 'the fewest changes at the residue, the rest unchanged, with p_upper at most 0.05 / residues_tested '
               '(enough for q < 0.05)'),
              ('result', 'string', 'pass; fail; not_tested', 'pass = q_upper < 0.05', 'the outcome'), ('reason', 'string', '', '', 'why')):
        rows.append('\t'.join([f, *c]))
    f = f'supporting/propensity_validation_{version}.tsv'
    for c in (('source', 'string', '', '', 'the change table'), ('cohort', 'string', '', '', 'all samples'),
              ('units_available', 'integer', '', '', 'protein and residue-type units tested with at least 2 changes'),
              ('units', 'integer', '', '', 'units drawn (up to 1,000)'), ('residues', 'integer', '', '', 'their residues'),
              ('binomial_q05', 'integer', '', '', 'residues with binomial BH q < 0.05 (BH within the units drawn)'),
              ('permutation_q05', 'integer', '', '', 'residues with permutation BH q < 0.05 (each unit\'s changes redrawn over its residues by '
               'mutability, multinomial; p = (1 + draws at or above k) / (iterations + 1))'),
              ('both_q05', 'integer', '', '', 'residues with both (the permutation p is at least 1 / (iterations + 1), so over tens of thousands of '
               'residues it cannot reach q < 0.05: compare the nominal calls)'),
              ('agreement', 'number', '', '', 'share of residues on which the two q < 0.05 calls agree'),
              ('binomial_p05', 'integer', '', '', 'residues with binomial p < 0.05'), ('permutation_p05', 'integer', '', '', 'residues with permutation p < 0.05'),
              ('both_p05', 'integer', '', '', 'residues with both'), ('agreement_p05', 'number', '', '', 'share of residues on which the two p < 0.05 calls agree'),
              ('median_abs_p_difference', 'number', '', '', 'median |binomial p - permutation p|'),
              ('max_abs_p_difference', 'number', '', '', 'largest |binomial p - permutation p|'),
              ('max_abs_p_difference_p_above_0_001', 'number', '', '', 'the same, where binomial p > 0.001'),
              ('iterations', 'integer', '', '', 'draws per unit'), ('seed', 'integer', '', '', 'NumPy default_rng seed')):
        rows.append('\t'.join([f, *c]))
    return rows


def write_docs(args, prev, out, version, previous, arginine_sets, carried, edited, changes_name, sites, files, stats, set_table, controls,
               validation, staging, burden, hyper, profiled, cds, mut, canonical) -> None:
    config = json.loads(Path(args.config).read_text(encoding='utf-8'))
    # data dictionary: the previous one under the new names, then the new tables
    lines = (prev / 'data_dictionary.tsv').read_text(encoding='utf-8').rstrip('\n').split('\n')
    dict_lines = [lines[0]]
    for line in lines[1:]:
        name, rest = line.split('\t', 1)
        name = changes_name if name.startswith('supporting/site_changes_') else x65.new_name(name, previous, version, arginine_sets)
        dict_lines.append(name + '\t' + rest)
    dict_lines += dictionary_rows(version)
    (out / 'data_dictionary.tsv').write_text('\n'.join(dict_lines) + '\n', encoding='utf-8')

    rows = x65.count_marks(out, version, previous, config)
    prev_release = json.loads((prev / 'release.json').read_text(encoding='utf-8'))
    commit = code_commit()
    mb = staging['coding_bp'] / 1e6
    matched = staging['entries_by_status'].get('matched', 0)
    n_entries = staging['entries']
    tcga_h, dm_h = len(hyper['tcga_mc3']), len(hyper['depmap_24q4'])

    def per_source_table() -> list[str]:
        t = ['| Source | Residue | Rows | Residues tested | q < 0.05, more changes than expected | of which known sites | q < 0.05, fewer (gnomAD) |',
             '|---|---|---:|---:|---:|---:|---:|']
        for src in SOURCES:
            for res in RESIDUES:
                s = stats[f'{src}/{res}']
                a = s['all']
                lower = f'{a["q_lower_below_0_05"]:,}' if 'q_lower_below_0_05' in a else ''
                t.append(f'| {SOURCES[src]} | {res} | {s["rows"]:,} | {a["tested"]:,} | {a["q_upper_below_0_05"]:,} | {a["known_sites_q_upper_below_0_05"]:,} | {lower} |')
        return t

    def set_md() -> list[str]:
        t = ['| Site set | Residue | Source | Cohort | Proteins | Sites | Observed | Expected | O/E | O/E (unweighted) | p (more) | p (fewer) |',
             '|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
        for _, r in set_table.iterrows():
            t.append(f'| {r["label"]} | {r["residue"]} | {SOURCES[r["source"]]} | {r["cohort"]} | {r["proteins"]:,} | {r["sites"]:,} | {r["observed"]:,} | '
                     f'{r["expected"]} | {r["oe"]} | {r["oe_unweighted"]} | {r["p_upper"]} | {r["p_lower"]} |')
        return t

    def controls_md() -> list[str]:
        t = ['| Control | Entry | Source | Cohort | Family | k | n | Expected | q | Result |', '|---|---|---|---|---|---:|---:|---:|---:|---|']
        for _, r in controls.iterrows():
            t.append(f'| {r["control"]} | {r["canonical_UniProtAC"]} ({r["gene"]}) | {SOURCES[r["source"]]} | {r["cohort"]} | {r["family"]} | {r["k"]} | '
                     f'{r["n"]} | {r["expected"]} | {r["q_upper"]} | {r["result"]} |')
        return t

    def validation_md() -> list[str]:
        t = ['| Source | Units | Residues | Binomial p < 0.05 | Permutation p < 0.05 | Both | Agreement (p < 0.05) | Binomial q < 0.05 | Largest p difference |',
             '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
        for _, r in validation.iterrows():
            t.append(f'| {SOURCES[r["source"]]} | {r["units"]:,} | {r["residues"]:,} | {r["binomial_p05"]:,} | {r["permutation_p05"]:,} | {r["both_p05"]:,} | '
                     f'{r["agreement_p05"]} | {r["binomial_q05"]:,} | {r["max_abs_p_difference"]} |')
        return t

    onco = controls[(controls['kind'] == 'oncohistone') & (controls['family'] == 'every residue of the entry') & (controls['source'] == 'tcga_mc3')
                    & (controls['cohort'] == 'all')]
    hot = controls[(controls['kind'] == 'recurrent hotspot') & (controls['source'] == 'tcga_mc3') & (controls['cohort'] == 'all')]
    sizes = ', '.join(f'`{Path(k).name}` {v["bytes"] / 1e6:.1f} MB' for k, v in files.items() if 'bytes' in v)
    new_section = [
        '## New: mutation propensity at every candidate residue', '',
        f'{version} asks, for each source, whether a residue is changed more (or less) often than its codon\'s mutability predicts, '
        'compared with the other residues of its type in the same protein. It reads the change tables of the release (TCGA, DepMap, '
        'ClinVar, gnomAD) and adds:', '',
        f'- `supporting/site_propensity_<source>_{version}.tsv.gz` (tcga_mc3, depmap_24q4, clinvar, gnomad_v4_1): every candidate residue '
        '(R, K, S, T, Y) that is a known site of a site set, changed or not, and every other candidate residue with a change: observed '
        'changes k, expected changes, O/E, binomial p and BH q, with an unweighted expectation beside it; TCGA and DepMap also without '
        'the hypermutated samples; gnomAD also the lower tail (depletion);',
        f'- `supporting/site_frequency_by_group_{version}.tsv.gz`: samples changed and profiled per residue and TCGA project or DepMap '
        'lineage (for waterfall and oncoprint views; the change tables keep the sample of every change);',
        f'- `supporting/group_denominators_{version}.tsv`: the samples profiled per TCGA project and DepMap lineage, all and not '
        'hypermutated, every group included (so a group without a changed sample shows k = 0);',
        f'- `supporting/sample_burden_{version}.tsv`: every TCGA sample and DepMap line, its coding calls, burden per Mb and whether it is hypermutated;',
        f'- `supporting/mutation_rates_{version}.tsv`: the spectra and rates the mutability comes from;',
        f'- `supporting/site_set_propensity_{version}.tsv`: per site set, residue, source and cohort, the changes at the set\'s sites '
        'within its proteins against expectation;',
        f'- `supporting/propensity_controls_{version}.tsv` and `supporting/propensity_validation_{version}.tsv`: positive controls and a '
        'permutation check.', '',
        '### Method', '',
        '- **Observed (k).** Per residue: TCGA tumor samples, DepMap cell lines, ClinVar variations or gnomAD variants with a missense or '
        'stop change there, from the change tables. n is the sum of k over the protein\'s residues of the same type.',
        '- **Mutability (w).** For each base of the residue\'s codon and each substitution that changes the residue (missense or stop), '
        'the rate of that substitution in its trinucleotide; w is their sum. The codons are those of the Ensembl 105 transcripts UniProt '
        f'assigns to the entry\'s displayed isoform whose coding sequence translates to the 2026_03 canonical sequence ({matched:,} of '
        f'{n_entries:,} entries; the others have no weights, and their rows leave the weighted columns empty). A base at an exon edge '
        'takes its intron neighbor. An entry encoded by several genes (the histones) takes the mean over its coding sequences.',
        '- **Rates.** TCGA and DepMap: each source\'s own spectrum, the residue-changing SNVs per pyrimidine-centered trinucleotide class '
        'in those coding sequences over the positions where the class would change the residue. gnomAD and ClinVar: gnomAD\'s germline '
        f'rates, {GNOMAD_RATES}, with the CpG methylation level of each CpG base from gnomAD\'s methylation table of the same paper '
        '(GRCh37, `methylation.txt.bgz`; positions lifted from GRCh38 with Ensembl\'s chain) binned as gnomAD\'s constraint code bins it '
        '(level 2 above a mean of 0.6, level 1 above 0.2). gnomAD v4 publishes its methylation sites but no per-context rate table; '
        'this table is the one gnomAD publishes.',
        '- **Test.** expected = n × w / (the sum of w over the protein\'s residues of the type); k ~ Binomial(n, w / sum w). p_upper is '
        'the chance of k or more; for gnomAD also p_lower, k or fewer. Benjamini–Hochberg q within each source, cohort and residue type, '
        'over every residue tested (every residue of a protein with n ≥ 1 and weights), not only the rows written.',
        f'- **Hypermutated samples.** More than 10 coding calls per Mb ({CAMPBELL}), over the {mb:.1f} Mb of the coding sequences above: '
        f'{tcga_h:,} of {profiled["tcga_mc3"]["all"]:,} TCGA samples and {dm_h:,} of {profiled["depmap_24q4"]["all"]:,} DepMap lines. '
        'Columns ending `_no_hypermutated` repeat the test without them, with spectra from the other samples. DepMap\'s mutation file keeps '
        'synonymous calls only where they touch a splice region, so its burden counts protein-changing calls almost only.',
        '- **ClinVar.** ClinVar\'s variants reflect clinical testing and submission more than mutation. Its O/E uses the germline '
        'rates as gnomAD\'s does, but is descriptive: `oe_status` says `descriptive` on its rows, and the data dictionary says so.',
        '- **Site sets.** For each set, residue and source, over the proteins with a site of the set: the changes at the set\'s sites '
        'against the sum over the proteins of n times the sites\' share of the protein\'s mutability. p from 10,000 draws, each protein\'s '
        'changes redrawn over its residues by mutability (a fixed seed per row). The same code serves any set of residues. '
        '`sites_total`, `sites_weighted` and `k_total` give the scope: the test covers sites_weighted of sites_total sites and '
        'observed of k_total changes (the rest lie on proteins without weights).', '',
        *per_source_table(), '',
        '### Site sets', '', 'A p of 9.999e-05 is the smallest 10,000 draws give ((1 + 0) / 10,001).', '', *set_md(), '',
        '### Controls', '',
        'The oncohistone controls are computed over every residue of the H3 entries (glycine 34 is not a candidate residue), from the '
        'TCGA and DepMap changes at every residue of H3 (UniProt counts the initiator methionine: H3K27 is K28, H3K36 is K37, H3G34 is '
        f'G35). The recurrent hotspots ({CHANG}) are read from the per-site tables. A control passes at q < 0.05; a failed control '
        'does not hold the release, and its row says why.', '',
        f'TCGA, all samples: {", ".join(f"{r.control} on {r.canonical_UniProtAC} {r.result}" for r in onco.itertuples())}; hotspots '
        f'{int((hot["result"] == "pass").sum())} of {len(hot)} pass.', '',
        *controls_md(), '',
        '### Permutation check', '',
        'On up to 1,000 proteins per source (protein and residue-type units with at least two changes, drawn with a fixed seed), each '
        'unit\'s changes were redrawn over its residues by mutability 10,000 times, and the permutation p compared with the binomial p. '
        'A residue\'s count under that redrawing is binomial, so the two should agree up to sampling error; the permutation p cannot go '
        'below 1 / 10,001, so it cannot reach q < 0.05 over tens of thousands of residues, and the calls are compared at p < 0.05.', '',
        *validation_md(), '',
        f'Sizes: {sizes}.', '']
    carried_text = re.sub(r'(^|\n)## New: (.)', lambda m: f'{m.group(1)}## {m.group(2).upper()}',
                          x65.carried_readme(prev, previous, version, arginine_sets).rstrip('\n'))
    R = [f'# Human methylarginine sites, context marks, site-level interaction evidence, cancer and germline variants and their propensity, {version}', '',
         f'{version} holds everything of {previous} under this release\'s names, its data unchanged, and adds the mutation propensity '
         'of every candidate residue in each source: observed against expected changes under codon mutability, within each protein, '
         'with tests at the residue and site-set level, frequencies per tumor type and lineage, and each sample\'s mutation burden. '
         'See CHANGES.md for what changed.', '',
         *new_section,
         carried_text]
    (out / 'README.md').write_text(x65.americanize('\n'.join(R)) + '\n', encoding='utf-8')
    C = [f'# Changes from {previous} to {version}', '',
         f'Everything of {previous} is in {version} under this release\'s names, its data unchanged. {version} adds the mutation '
         'propensity of every candidate residue (the decisions of the owner of 20 November 2026: hypermutated samples above 10 coding '
         'calls per Mb, with and without them; gnomAD\'s germline rates for gnomAD and ClinVar; the oncohistone controls over every '
         'residue of H3; every known site and every changed residue; ClinVar descriptive; a failed control reported, not blocking).', '',
         *new_section,
         '## New files', '', *[f'- `{f}`' for f in sorted(files)], '',
         '## Files renamed', '',
         f'Every carried file is renamed `_{previous}` to `_{version}`, its content unchanged except where it names its version: '
         + ', '.join(f'`{e}`' for e in edited) + '. '
         f'`{changes_name}` replaces {previous}\'s site changes table; it lists every one of the {len(sites):,} methylarginine sites as kept.', '',
         'Code: `scripts/70_stage_propensity.py` (coding sequences, rates, burden, H3 changes), `scripts/propensity.py` (the method), '
         '`scripts/71_release_propensity.py` (this release), `scripts/65_carry_release.py` (carrying the previous release), at '
         f'commit `{commit}`.']
    (out / 'CHANGES.md').write_text(x65.americanize('\n'.join(C)) + '\n', encoding='utf-8')

    def rename(obj):
        return json.loads(json.dumps(obj).replace(f'_{previous}.', f'_{version}.'))

    release = {'title': 'Human methylarginine sites, context marks, site-level interaction evidence, cancer and germline variants and their propensity',
               'version': version, 'previous_version': previous, 'release_date': args.release_date,
               'repository': 'https://github.com/Shechterlab/PTM_Pipeline', 'license': prev_release.get('license'),
               'methylarginine_data': prev_release.get('methylarginine_data'), 'methylarginine_code_commit': prev_release.get('methylarginine_code_commit'),
               'context_marks_code_commit': prev_release.get('context_marks_code_commit'), 'release_code_commit': commit,
               'carried_over': carried, 'carried_over_text_edited': edited, 'site_changes': changes_name,
               'sources': prev_release.get('sources'), 'sources_excluded': prev_release.get('sources_excluded'),
               'context_marks': prev_release.get('context_marks'), 'interaction_tables': rename(prev_release.get('interaction_tables')),
               'counts': rows, 'cancer_changes': rename(prev_release.get('cancer_changes')), 'germline_variants': rename(prev_release.get('germline_variants')),
               'variant_changes_excluded': rename(prev_release.get('variant_changes_excluded')),
               'propensity': {'files': files, 'per_source_residue': stats, 'site_sets': set_table.to_dict(orient='records'),
                              'controls': controls.to_dict(orient='records'), 'validation': validation.to_dict(orient='records'),
                              'germline_rates': GNOMAD_RATES, 'hypermutated': {'threshold_per_mb': 10, 'coding_mb': round(mb, 4), 'citation': CAMPBELL,
                                                                               'tcga_mc3': tcga_h, 'depmap_24q4': dm_h},
                              'entries_with_weights': matched, 'entries': n_entries, 'staging_inputs': staging['inputs'], 'rules': staging['rules'],
                              'seed': P.SEED, 'iterations': P.ITERATIONS},
               'corrections': []}
    (out / 'release.json').write_text(json.dumps(release, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    names = sorted(p.relative_to(out).as_posix() for p in out.rglob('*') if p.is_file() and p.name != 'MD5SUMS')
    (out / 'MD5SUMS').write_text(''.join(f'{md5_file(out / n)}  {n}\n' for n in names), encoding='utf-8')


if __name__ == '__main__':
    main()

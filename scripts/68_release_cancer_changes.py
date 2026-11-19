#!/usr/bin/env python3
"""Write a release that carries the previous release forward and adds the somatic changes, in TCGA tumors and DepMap cell
lines, at every candidate residue of the site sets.

    python scripts/68_release_cancer_changes.py --version v1.6.0 --previous v1.5.1 --release-date 2026-11-19

It carries every file of releases/<previous> under this release's names, as scripts/65_carry_release.py does (the data
unchanged; the strings that name the version inside STUDIES.md, the analyses READMEs and MD5SUMS and the workbook title
change), and adds to supporting/, from the staged layer (inputs/layers/, scripts/44b, 44 and 66):

  tcga_mc3_changes_<X>_<version>.tsv.gz      per residue X (R, K, S, T, Y): one row per tumor sample and change at an X of a
                                             reviewed canonical protein, a known site or not, and site_sets
  depmap_24q4_changes_<X>_<version>.tsv.gz   the same for DepMap cell lines, and site_sets
  cancer_change_protein_background_<version>.tsv.gz   per entry and residue class: residues, changes, samples or lines,
                                             positions (the staged bytes)
  tcga_mc3_projects_<version>.tsv            the TCGA projects and their tumor samples (the staged bytes)
  depmap_24q4_change_models_<version>.tsv    the DepMap cell lines of the change tables (the staged bytes)

site_sets names the site sets of this release (config/site_sets.json) that have a site at the row's accession and
position, residue included, separated by semicolons; empty, the residue is a candidate but no known site. A website can
then show a site set's changes at every candidate residue, or at its known sites only. The arginine rows are the
arginine layer's tables of earlier releases' inputs (scripts/44), gzip-compressed; the methylarginine release files
themselves are carried unchanged.

It then writes supporting/site_changes_<previous>_to_<version>.tsv, data_dictionary.tsv (the previous one under the new
names, and the new tables), README.md (what is new, then the previous release's README from its first section on),
CHANGES.md, release.json and MD5SUMS. Reader-facing prose is in American spelling.
"""
from __future__ import annotations

import argparse
import gzip
import importlib.util
import io
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from common import md5_file  # noqa: E402


def _load(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, HERE / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


x65 = _load('carry65', '65_carry_release.py')
LAYERS = ROOT / 'inputs' / 'layers'
RESIDUES = ['R', 'K', 'S', 'T', 'Y']
SOURCES = {'tcga_mc3': 'TCGA MC3', 'depmap_24q4': 'DepMap 24Q4'}
NAMES = {'R': 'arginine', 'K': 'lysine', 'S': 'serine', 'T': 'threonine', 'Y': 'tyrosine'}
COPIED = {'cancer_change_protein_background.tsv.gz': 'cancer_change_protein_background_{v}.tsv.gz',
          'tcga_mc3_projects.tsv': 'tcga_mc3_projects_{v}.tsv', 'depmap_24q4_change_models.tsv': 'depmap_24q4_change_models_{v}.tsv'}


def staged(src: str, res: str) -> Path:
    return LAYERS / (f'{src}_arginine_changes.tsv' if res == 'R' else f'{src}_changes_{res}.tsv.gz')


def released(src: str, res: str, version: str) -> str:
    return f'supporting/{src}_changes_{res}_{version}.tsv.gz'


def write_gz(df: pd.DataFrame, path: Path) -> None:
    buf = io.StringIO()
    df.to_csv(buf, sep='\t', index=False, lineterminator='\n')
    with open(path, 'wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0, compresslevel=9) as gz:
        gz.write(buf.getvalue().encode('utf-8'))


def site_index(out: Path, version: str, config: dict) -> dict:
    """(accession, position, residue) -> the site sets with a site there, from this release's site tables."""
    index: dict = {}
    for name, spec in config['site_sets'].items():
        letters = set(spec['residues'])
        cols = spec['columns']
        if letters == {'R'}:
            t = x65.read(out / f'human_methylarginine_sites_{version}.tsv')
            keys = zip(t[cols['accession']], t[cols['position']].astype(int), ['R'] * len(t))
        else:
            parts = []
            for table in spec['tables']:
                t = x65.read(out / x65.release_name(table, version))
                parts.append(t[[cols['accession'], cols['position'], cols['residue']]])
            t = pd.concat(parts)
            t = t[t[cols['residue']].isin(letters)]
            keys = zip(t[cols['accession']], t[cols['position']].astype(int), t[cols['residue']])
        for key in set(keys):
            index.setdefault(key, set()).add(name)
    return {k: ';'.join(sorted(v)) for k, v in index.items()}


def dictionary_rows(version: str) -> list[str]:
    ref = 'reviewed UniProtKB accession of the canonical sequence in the reference proteome (UniProt 2026_03, UP000005640)'
    site_sets = ('the site sets of this release with a site at this accession and position (config/site_sets.json), separated by '
                 'semicolons; empty: a candidate residue, not a known site')
    alt = ('string', 'ACDEFGHIKLMNPQRSTVWY*', '', 'the residue after the change; * a stop')
    rows = []
    for res in RESIDUES:
        f = released('tcga_mc3', res, version)
        for col, typ, vals, means, desc in (
                ('canonical_UniProtAC', 'string', '', '', ref),
                ('position', 'integer', '', '', f'position of the {NAMES[res]} on the canonical sequence'),
                ('alt', *alt),
                ('project_id', 'string', '', '', f'TCGA project of the tumor (supporting/tcga_mc3_projects_{version}.tsv)'),
                ('sample_id', 'string', '', '', 'TCGA tumor sample: patient barcode and sample type, e.g. TCGA-02-0003-01'),
                ('transcript_id', 'string', '', '', 'Ensembl transcript the change was read on'),
                ('mapped_by', 'string', 'canonical_transcript; maf_transcript',
                 'canonical_transcript = read on a transcript UniProt assigns to the canonical isoform (all_effects of the MAF); '
                 'maf_transcript = read on the MAF\'s own transcript and placed through the entry it names or the gene symbol',
                 'how the change was placed on the canonical sequence'),
                ('filter', 'string', 'PASS; wga', 'PASS = passed every MC3 filter; wga = flagged only as from whole-genome-amplified '
                 'DNA, kept in TCGA-OV and TCGA-LAML as Ellrott et al. 2018 kept them', 'MC3 FILTER of the call'),
                ('tumor_alt_count', 'integer', '', '', 'tumor reads with the variant'),
                ('tumor_ref_count', 'integer', '', '', 'tumor reads with the reference allele'),
                ('site_sets', 'string', '', '', site_sets)):
            rows.append('\t'.join([f, col, typ, vals, means, desc]))
        f = released('depmap_24q4', res, version)
        for col, typ, vals, means, desc in (
                ('canonical_UniProtAC', 'string', '', '', ref),
                ('position', 'integer', '', '', f'position of the {NAMES[res]} on the canonical sequence'),
                ('alt', *alt),
                ('protein_change', 'string', '', '', 'residue, position and new residue, e.g. K120E or K79*'),
                ('model_id', 'string', '', '', f'DepMap model ID (supporting/depmap_24q4_change_models_{version}.tsv)'),
                ('allele_fraction', 'number', '', '', 'variant allele fraction in the cell line, as DepMap gives it'),
                ('hotspot', 'boolean', 'True; False', '', 'DepMap marks the change a hotspot'),
                ('likely_lof', 'boolean', 'True; False', '', 'DepMap marks the change likely loss of function'),
                ('am_pathogenicity', 'number', '', '', 'AlphaMissense pathogenicity (0 to 1), as DepMap gives it; empty for stops'),
                ('am_class', 'string', 'likely_benign; ambiguous; likely_pathogenic', '', 'AlphaMissense class, as DepMap gives it'),
                ('site_sets', 'string', '', '', site_sets)):
            rows.append('\t'.join([f, col, typ, vals, means, desc]))
    f = f'supporting/cancer_change_protein_background_{version}.tsv.gz'
    for col, typ, vals, desc in (
            ('canonical_UniProtAC', 'string', '', ref),
            ('residue', 'string', 'R; K; S; T; Y; all', 'residue class; all = every residue of the protein'),
            ('residues', 'integer', '', 'residues of that class in the canonical sequence (its length for all)'),
            ('tcga_changes', 'integer', '', 'TCGA changes on them (rows of the change tables; for all, every residue)'),
            ('tcga_samples', 'integer', '', 'tumor samples with one of those changes'),
            ('tcga_positions', 'integer', '', 'positions changed in TCGA'),
            ('depmap_changes', 'integer', '', 'DepMap changes on them'),
            ('depmap_lines', 'integer', '', 'cell lines with one of those changes'),
            ('depmap_positions', 'integer', '', 'positions changed in DepMap')):
        rows.append('\t'.join([f, col, typ, vals, '', desc]))
    f = f'supporting/tcga_mc3_projects_{version}.tsv'
    for col, typ, desc in (('project_id', 'string', 'TCGA project'), ('name', 'string', 'project name'),
                           ('tumour_samples', 'integer', 'tumor samples of the project in the MC3 MAF (every one counts as sequenced)'),
                           ('tumour_samples_with_kept_call', 'integer', 'of those, the samples with a call kept by the MC3 filters'),
                           ('arginine_changes', 'integer', 'arginine changes extracted for the project (the arginine layer)')):
        rows.append('\t'.join([f, col, typ, '', '', desc]))
    f = f'supporting/depmap_24q4_change_models_{version}.tsv'
    for col, desc in (('model_id', 'DepMap model ID'), ('cell_line', 'cell line name'), ('lineage', 'Oncotree lineage'),
                      ('primary_disease', 'Oncotree primary disease')):
        rows.append('\t'.join([f, col, 'string', '', '', desc]))
    return rows


def code_commit() -> str:
    paths = ['scripts/68_release_cancer_changes.py', 'scripts/65_carry_release.py', 'scripts/66_stage_cancer_changes.py', 'config/site_sets.json']
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

    # the change tables, with site_sets
    index = site_index(out, version, config)
    stats, files = {}, {}
    for src in SOURCES:
        for res in RESIDUES:
            t = x65.read(staged(src, res))
            t['site_sets'] = [index.get((a, int(p), res), '') for a, p in zip(t['canonical_UniProtAC'], t['position'])]
            name = released(src, res, version)
            write_gz(t, out / name)
            files[name] = {'from': staged(src, res).relative_to(ROOT).as_posix(), 'from_md5': md5_file(staged(src, res)), 'md5': md5_file(out / name)}
            known = t[t['site_sets'] != '']
            stats[f'{src}/{res}'] = {'rows': len(t), 'positions': int(t[['canonical_UniProtAC', 'position']].drop_duplicates().shape[0]),
                                     'rows_at_known_sites': len(known),
                                     'known_sites_changed': int(known[['canonical_UniProtAC', 'position']].drop_duplicates().shape[0])}
    for src_name, dest in COPIED.items():
        name = 'supporting/' + dest.format(v=version)
        shutil.copyfile(LAYERS / src_name, out / name)
        files[name] = {'from': f'inputs/layers/{src_name}', 'from_md5': md5_file(LAYERS / src_name), 'md5': md5_file(out / name)}

    # per site set: known sites, and those with a change
    per_set = []
    for name, spec in config['site_sets'].items():
        letters = sorted(spec['residues'])
        n_sites = sum(1 for k, v in index.items() if name in v.split(';'))
        row = {'set': name, 'label': x65.SET_LABELS.get(name, name), 'residues': ''.join(letters), 'sites': n_sites}
        for src in SOURCES:
            changed = set()
            for res in letters:
                t = x65.read(out / released(src, res, version))
                t = t[t['site_sets'].str.split(';').map(lambda s: name in s)]
                changed |= set(zip(t['canonical_UniProtAC'], t['position']))
            row[f'{src}_sites_changed'] = len(changed)
        per_set.append(row)

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
    staging = json.loads((LAYERS / 'cancer_changes_summary.json').read_text(encoding='utf-8'))
    den = staging['denominators']
    commit = code_commit()

    # the new tables, for README and CHANGES
    def residue_table() -> list[str]:
        out_lines = ['| Residue | TCGA changes | at known sites | DepMap changes | at known sites |', '|---|---:|---:|---:|---:|']
        for res in RESIDUES:
            tc, dm = stats[f'tcga_mc3/{res}'], stats[f'depmap_24q4/{res}']
            out_lines.append(f'| {NAMES[res].capitalize()} ({res}) | {tc["rows"]:,} | {tc["rows_at_known_sites"]:,} | {dm["rows"]:,} | {dm["rows_at_known_sites"]:,} |')
        return out_lines

    def set_table() -> list[str]:
        out_lines = ['| Site set | Sites | Sites changed in TCGA tumors | Sites changed in DepMap cell lines |', '|---|---:|---:|---:|']
        for r in per_set:
            out_lines.append(f'| {r["label"]} | {r["sites"]:,} | {r["tcga_mc3_sites_changed"]:,} | {r["depmap_24q4_sites_changed"]:,} |')
        return out_lines

    new_section = [
        '## New: somatic changes at every candidate residue', '',
        f'`supporting/tcga_mc3_changes_<residue>_{version}.tsv.gz` and `supporting/depmap_24q4_changes_<residue>_{version}.tsv.gz` '
        'list the somatic missense and stop changes at every arginine, lysine, serine, threonine and tyrosine of the reviewed '
        'canonical proteins, a known site or not: in TCGA tumors (the MC3 calls of the PanCancer Atlas, '
        f'{den["tcga_tumor_samples"]:,} tumor samples in {den["tcga_projects"]} projects) and in DepMap 24Q4 cell lines '
        f'({den["depmap_lines_with_mutation_calls"]:,} lines with mutation calls). One row per tumor sample or cell line and change. '
        'The column `site_sets` names the site sets of this release with a site at that residue (empty: a candidate residue, not '
        'a known site), so a set\'s changes can be read at every candidate residue or at its known sites only.', '',
        *residue_table(), '',
        *set_table(), '',
        f'`supporting/cancer_change_protein_background_{version}.tsv.gz` gives, for each protein and residue class (and over all '
        'its residues), how many residues it has and the changes, samples or lines and positions on them, so a site\'s changes '
        'can be compared with those of the other residues of its kind in the same protein. Arginines are changed far more often '
        'than the other residues (their codons hold CpG dinucleotides, which mutate often), so such comparisons are best made '
        'within a residue class.', '',
        'The rules are those of the arginine layer of earlier releases (`scripts/44b_extract_tcga_mc3.py`, '
        '`scripts/44_stage_depmap.py`): TCGA calls with FILTER PASS (in OV and LAML also `wga` alone), read on a transcript UniProt '
        'assigns to the canonical isoform, else on the MAF\'s own transcript placed through its entry or gene symbol; DepMap calls '
        'on the canonical isoform; a change counts only where the 2026_03 canonical sequence has that residue. '
        '`scripts/66_stage_cancer_changes.py` stages the K, S, T and Y tables and, rerun on arginines, writes the arginine tables '
        'byte for byte. TCGA data are GDC open access (NIH Genomic Data Sharing policy; acknowledge TCGA and the MC3 publication, '
        'Ellrott et al., Cell Syst 6:271, 2018, and make no attempt to identify participants); DepMap 24Q4 is CC BY 4.0 '
        '(doi:10.25452/figshare.plus.27993248).', '',
        f'`supporting/tcga_mc3_projects_{version}.tsv` names the TCGA projects and their tumor samples, and '
        f'`supporting/depmap_24q4_change_models_{version}.tsv` the cell lines (name, lineage, primary disease).', '']

    R = [f'# Human methylarginine sites, context marks, site-level interaction evidence and cancer changes, {version}', '',
         f'{version} holds everything of {previous} under this release\'s names, its data unchanged: the methylarginine data of '
         'v1.4.3, the tables of context marks and the site-level interaction tables of every site set. It adds the somatic changes '
         'in TCGA tumors and DepMap cell lines at every candidate residue of the site sets, with the site sets each residue belongs '
         'to, and per-protein background counts. See CHANGES.md for what changed.', '',
         *new_section,
         x65.carried_readme(prev, previous, version, arginine_sets).rstrip('\n')]
    (out / 'README.md').write_text(x65.americanize('\n'.join(R)) + '\n', encoding='utf-8')

    new_files = sorted(files)
    C = [f'# Changes from {previous} to {version}', '',
         f'Everything of {previous} is in {version} under this release\'s names, its data unchanged. {version} adds the somatic changes '
         'at every candidate residue of the site sets, in TCGA tumors and DepMap cell lines.', '',
         *new_section,
         '## New files', '', *[f'- `{f}`' for f in new_files], '',
         '## Files renamed', '',
         f'Every carried file is renamed `_{previous}` to `_{version}`, its content unchanged except where it names its version: '
         + ', '.join(f'`{e}`' for e in edited) + '. '
         f'`{changes_name}` replaces {previous}\'s site changes table; it lists every one of the {len(sites):,} methylarginine sites as kept.', '',
         'Code: `scripts/66_stage_cancer_changes.py` (the K, S, T and Y change tables and the background), '
         '`scripts/68_release_cancer_changes.py` (this release), `scripts/65_carry_release.py` (carrying the previous release), at '
         f'commit `{commit}`.']
    (out / 'CHANGES.md').write_text(x65.americanize('\n'.join(C)) + '\n', encoding='utf-8')

    release = {'title': 'Human methylarginine sites, context marks, site-level interaction evidence and cancer changes', 'version': version,
               'previous_version': previous, 'release_date': args.release_date, 'repository': 'https://github.com/Shechterlab/PTM_Pipeline',
               'license': prev_release.get('license'), 'methylarginine_data': 'unchanged from v1.4.3: files carried over byte for byte, not rebuilt',
               'methylarginine_code_commit': prev_release.get('methylarginine_code_commit'),
               'context_marks_code_commit': prev_release.get('context_marks_code_commit'), 'release_code_commit': commit,
               'carried_over': carried, 'carried_over_text_edited': edited, 'site_changes': changes_name,
               'sources': prev_release.get('sources'), 'sources_excluded': prev_release.get('sources_excluded'),
               'context_marks': prev_release.get('context_marks'), 'interaction_tables': prev_release.get('interaction_tables') and {
                   x65.new_name(k, previous, version, arginine_sets): md5_file(out / x65.new_name(k, previous, version, arginine_sets))
                   for k in prev_release['interaction_tables']},
               'counts': rows,
               'cancer_changes': {'files': files, 'per_residue': stats, 'per_site_set': per_set, 'denominators': den,
                                  'staging_inputs': staging['inputs'], 'rules': staging['rule']}}
    (out / 'release.json').write_text(json.dumps(release, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    names = sorted(p.relative_to(out).as_posix() for p in out.rglob('*') if p.is_file() and p.name != 'MD5SUMS')
    (out / 'MD5SUMS').write_text(''.join(f'{md5_file(out / n)}  {n}\n' for n in names), encoding='utf-8')
    print(json.dumps({'per_residue': stats, 'per_site_set': per_set}, indent=1))


if __name__ == '__main__':
    main()

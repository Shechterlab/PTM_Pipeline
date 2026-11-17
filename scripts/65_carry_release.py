#!/usr/bin/env python3
"""Write a release that carries the previous release forward and holds the site-level interaction tables of every site set.

    python scripts/65_carry_release.py --version v1.5.1 --previous v1.5.0 --release-date 2026-11-17

Run it after scripts/64 has written the interaction tables into releases/<version>/supporting/ (one file per site set
and table, <set>_<table>_<version>.tsv; the arginine sets' tables copied there with --carry-over). It then:

  carries every file of releases/<previous> under this release's names (_<previous> becomes _<version>), the data
  unchanged. The strings that name the version inside STUDIES.md and the two analyses READMEs and MD5SUMS change, as
  does the workbook's title; the analyses MD5SUMS are checksummed again. The previous release's own documents
  (README.md, CHANGES.md, MD5SUMS, release.json, data_dictionary.tsv) and its site_changes table are not carried. The
  methylarginine interaction tables of the previous release (supporting/<table>_<previous>.tsv) are the arginine site
  set's tables of this release (supporting/methylarginine_<table>_<version>.tsv): the build checks they are the same
  bytes and does not write them twice;
  writes supporting/site_changes_<previous>_to_<version>.tsv (every site kept), data_dictionary.tsv (the previous
  one under the new names, and the interaction tables of every site set), README.md (the per-mark counts, the
  interaction tables, then the previous release's README from its first section on), CHANGES.md, release.json and
  MD5SUMS. Reader-facing prose is in American spelling; file names, column names, codes and quoted titles keep theirs.

Counts per site set (README, release.json):
  sites                      distinct residues in the site set's tables (config/site_sets.json; the main site table
                             for methylarginine)
  high-confidence sites      the website's control, not a release tier (release tables keep every row). For the context
                             marks (DECISIONS_v1.5.1 D5, as ptmcontext's site files): only experimental rows give units; a
                             unit is a PubMed ID in pmids, or a direct dataset (Berryhill 2023, Ochoa 2020, Hendriks 2017 as
                             PubMed ID 28112733, Akimov 2018 as 29967540, so a database citing the paper is not counted
                             again); databases count only through the papers they cite; a residue's units are pooled over
                             all its rows; high confidence is two units or more. For methylarginine, its own control,
                             unchanged: two or more studies in `studies`, a laboratory's studies counted once, PhosphoSitePlus,
                             UniProt, IEDB and ProMetheusDB not counted.
  with interaction evidence  sites with a row of <set>_site_interactions that counts as evidence (counts_as_evidence true);
                             for methylarginine also the count from the kinds of evidence the generic step computes
"""
from __future__ import annotations

import argparse
import importlib.util
import io
import json
import re
import subprocess
import sys
import tempfile
import zipfile
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


NOT_CARRIED = {'README.md', 'CHANGES.md', 'MD5SUMS', 'data_dictionary.tsv', 'release.json'}
VERSIONED_TEXT = {'analyses/README.md', 'analyses/MD5SUMS', 'analyses_with_lower_confidence/README.md',
                  'analyses_with_lower_confidence/MD5SUMS', 'STUDIES.md'}
INTERACTION_TABLES = ('site_interactions', 'residue_evidence', 'nucleic_acid_contacts', 'nucleic_acid_near_sites',
                      'nucleic_acid_proximal_sites', 'site_interaction_records')
# the kinds of evidence scripts/64 computes for any site set (methylarginine's tables hold more: predicted interfaces and
# contacts, IntAct features, RBDmap and pCLAP)
GENERIC_EVIDENCE = {'pdb_protein', 'pdb_rna', 'pdb_dna', 'pdb_self', 'rna_crosslink', 'capri_xl_peptide', 'capri_adj_peptide',
                    'rbrid_residue', 'dna_crosslink', 'dna_binding_region'}
# DECISIONS_v1.5.1 D5 (ptmcontext sites/methyllysine.yaml, lysine.yaml and phospho.yaml, control hc)
DIRECT_UNITS = {'Berryhill 2023': 'Berryhill 2023', 'Ochoa 2020': 'Ochoa 2020', 'Hendriks 2017': '28112733', 'Akimov 2018': '29967540'}
# methylarginine's own control (ptmcontext sites/methylarginine.yaml, control hc), unchanged
ARG_NOT_COUNTED = {'PSP', 'UniProt', 'IEDB', 'ProMetheusDB'}
ARG_GROUPS = {'Larsen2016': 'Nielsen laboratory', 'Sylvestersen2014': 'Nielsen laboratory', 'Lu2024': 'Lu laboratory',
              'Lu2025': 'Lu laboratory', 'Musiani2019': 'Bonaldi laboratory', 'Fong2019': 'Bonaldi laboratory'}
LYSINE_MARKS = (('acetylation', 'acetylation'), ('lysine_acylation', 'lysine acylation (all acyl types)'),
                ('ubiquitination', 'ubiquitination'), ('sumoylation', 'SUMOylation'))
SET_LABELS = {'methylarginine': 'Methylarginine (R)', 'methyllysine': 'Methyllysine (K)',
              'lysine': 'Lysine acylation, ubiquitination and SUMOylation (K)', 'phospho': 'Phosphorylation (S, T, Y)'}

# British spellings in prose and their American forms (DECISIONS_v1.5.1 D6). Code spans, quoted text and links are left as
# they are, so column names (`localisation`), codes and the titles of papers keep their spelling.
SPELLING = {
    'licence': 'license', 'licences': 'licenses', 'colour': 'color', 'colours': 'colors',
    'coloured': 'colored', 'localised': 'localized', 'localise': 'localize', 'localises': 'localizes',
    'localisation': 'localization', 'localisations': 'localizations', 'analysed': 'analyzed', 'analyse': 'analyze',
    'analysing': 'analyzing', 'behaviour': 'behavior', 'centre': 'center', 'centred': 'centered', 'tumour': 'tumor',
    'tumours': 'tumors', 'organisation': 'organization', 'organisations': 'organizations', 'recognised': 'recognized',
    'recognise': 'recognize', 'summarised': 'summarized', 'summarise': 'summarize', 'catalogue': 'catalog',
    'catalogued': 'cataloged', 'modelled': 'modeled', 'modelling': 'modeling', 'labelled': 'labeled', 'labelling': 'labeling',
    'favour': 'favor', 'favoured': 'favored', 'neighbour': 'neighbor', 'neighbours': 'neighbors', 'neighbouring': 'neighboring',
    'grey': 'gray', 'artefact': 'artifact', 'artefacts': 'artifacts', 'characterised': 'characterized',
    'characterise': 'characterize', 'normalised': 'normalized', 'normalise': 'normalize', 'standardised': 'standardized',
    'prioritised': 'prioritized', 'minimise': 'minimize', 'minimised': 'minimized', 'optimised': 'optimized',
    'fibre': 'fiber', 'fibres': 'fibers', 'judgement': 'judgment', 'travelled': 'traveled', 'cancelled': 'canceled',
    'signalling': 'signaling', 'synthesised': 'synthesized',
    'emphasise': 'emphasize', 'utilise': 'utilize', 'utilised': 'utilized', 'acknowledgement': 'acknowledgment',
    'realised': 'realized', 'harmonised': 'harmonized', 'categorised': 'categorized', 'visualised': 'visualized',
    'visualise': 'visualize', 'randomised': 'randomized', 'stabilised': 'stabilized', 'mobilised': 'mobilized',
    'generalise': 'generalize', 'generalised': 'generalized', 'specialised': 'specialized', 'parametrised': 'parametrized',
}
_PROTECT = re.compile(r'`[^`\n]*`|"[^"\n]*"|“[^”\n]*”|https?://\S+|\]\([^)]*\)')
_WORD = re.compile(r'\b(' + '|'.join(sorted(SPELLING, key=len, reverse=True)) + r')\b', re.I)


def americanize(text: str) -> str:
    """Prose in American spelling: SPELLING applied outside code spans, quoted text and link targets."""
    def fix(m):
        word = m.group(0)
        new = SPELLING[word.lower()]
        return new.upper() if word.isupper() else new[0].upper() + new[1:] if word[0].isupper() else new
    out, pos = [], 0
    for m in _PROTECT.finditer(text):
        out.append(_WORD.sub(fix, text[pos:m.start()]))
        out.append(m.group(0))
        pos = m.end()
    out.append(_WORD.sub(fix, text[pos:]))
    return ''.join(out)


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


def new_name(rel: str, previous: str, version: str, arginine_sets: list[str]) -> str:
    """A carried file's name in this release: its version, and the arginine set's interaction tables under the set's name."""
    for t in INTERACTION_TABLES:
        if rel == f'supporting/{t}_{previous}.tsv':
            return f'supporting/{arginine_sets[0]}_{t}_{version}.tsv'
    return rel.replace(f'_{previous}', f'_{version}')


def retitle_workbook(data: bytes, previous: str, version: str) -> bytes:
    out = io.BytesIO()
    with zipfile.ZipFile(io.BytesIO(data)) as zin, zipfile.ZipFile(out, 'w') as zout:
        for info in zin.infolist():
            part = zin.read(info.filename)
            if info.filename in ('docProps/core.xml', 'xl/workbook.xml'):
                part = part.replace(f'Sites {previous}'.encode(), f'Sites {version}'.encode())
            zout.writestr(info, part, compress_type=info.compress_type)
    return out.getvalue()


def refresh_md5sums(path: Path) -> None:
    lines = [ln for ln in path.read_text(encoding='utf-8').splitlines() if ln.strip()]
    names = [ln.split('  ', 1)[1] for ln in lines]
    path.write_text(''.join(f'{md5_file(path.parent / n)}  {n}\n' for n in names), encoding='utf-8')


def carry(prev: Path, out: Path, previous: str, version: str, arginine_sets: list[str]) -> tuple[dict, list]:
    carried = {}
    for path in sorted(p for p in prev.rglob('*') if p.is_file()):
        rel = path.relative_to(prev).as_posix()
        if ('/' not in rel and rel in NOT_CARRIED) or Path(rel).name.startswith('site_changes_'):
            continue
        data = path.read_bytes()
        if rel in VERSIONED_TEXT:
            data = data.replace(previous.encode(), version.encode())
        elif rel.endswith('.xlsx'):
            data = retitle_workbook(data, previous, version)
        dest_rel = new_name(rel, previous, version, arginine_sets)
        dest = out / dest_rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest_rel.startswith(f'supporting/{arginine_sets[0]}_') and dest.exists():
            if dest.read_bytes() != data:
                raise SystemExit(f'{dest_rel} (scripts/64 --carry-over) differs from {prev.name}/{rel}')
        elif not dest.exists() or dest.read_bytes() != data:
            dest.write_bytes(data)
        carried[dest_rel] = {'from': rel, 'from_md5': md5_file(path), 'md5': md5_file(dest)}
    for sums in sorted(out.glob('*/MD5SUMS')):
        refresh_md5sums(sums)
        rel = sums.relative_to(out).as_posix()
        if rel in carried:
            carried[rel]['md5'] = md5_file(sums)
    edited = sorted(k for k, v in carried.items() if v['md5'] != v['from_md5'])
    return carried, edited


def mark_units(table: pd.DataFrame) -> pd.Series:
    """Per residue (accession, position): its distinct units, pooled over its rows (D5)."""
    exp = table[table['evidence_class'] == 'experimental']
    units: dict = {}
    for acc, pos, pmids, sources in zip(exp['canonical_UniProtAC'], exp['position'], exp['pmids'], exp['sources'], strict=True):
        got = units.setdefault((acc, pos), set())
        got.update(x.strip() for x in pmids.split(';') if x.strip())
        got.update(DIRECT_UNITS[s] for s in sources.split('; ') if s in DIRECT_UNITS)
    return pd.Series({k: len(v) for k, v in units.items()}, dtype='int64')


def arginine_high_confidence(sites: pd.DataFrame) -> int:
    def units(studies: str) -> int:
        names = {x.strip() for x in studies.split(';') if x.strip()} - ARG_NOT_COUNTED
        return len({ARG_GROUPS.get(n, n) for n in names})
    return int(sites['studies'].map(units).ge(2).sum())


def evidence_sites(path: Path, generic_only: bool = False) -> set:
    t = read(path)
    t = t[t['counts_as_evidence'] == 'true']
    if generic_only:
        t = t[t['evidence_type'].isin(GENERIC_EVIDENCE)]
    return set(t['site_id'])


def site_ids(table: pd.DataFrame, residue_col: str = 'residue') -> set:
    return set(table['canonical_UniProtAC'] + ':' + table[residue_col] + table['position'])


def count_marks(out: Path, version: str, previous: str, config: dict) -> list[dict]:
    """One row per site set (and, for lysine, per table): sites, high-confidence sites, sites with interaction evidence."""
    rows = []
    for name, spec in config['site_sets'].items():
        ev_path = out / 'supporting' / f'{name}_site_interactions_{version}.tsv'
        ev = evidence_sites(ev_path)
        if set(spec['residues']) == {'R'}:
            sites = read(out / f'human_methylarginine_sites_{version}.tsv')
            ids = set(sites['site_id'])
            rows.append({'set': name, 'label': SET_LABELS.get(name, name), 'sites': len(ids), 'high_confidence': arginine_high_confidence(sites),
                         'with_interaction_evidence': len(ev & ids), 'with_generic_interaction_evidence': len(evidence_sites(ev_path, True) & ids),
                         'high_confidence_rule': 'methylarginine control: two or more studies, a laboratory counted once, PhosphoSitePlus, UniProt, IEDB and ProMetheusDB not counted'})
            continue
        tables = {Path(t).name.split('_sites_')[0].removeprefix('human_'): read(out / Path(t).name.replace(f'_{previous}.', f'_{version}.'))
                  for t in spec['tables']}
        allt = pd.concat(tables.values(), ignore_index=True)
        units = mark_units(allt)
        ids = site_ids(allt)
        rows.append({'set': name, 'label': SET_LABELS.get(name, name), 'sites': len(ids), 'high_confidence': int((units >= 2).sum()),
                     'with_interaction_evidence': len(ev & ids), 'high_confidence_rule': 'D5: two or more units (PubMed IDs or direct datasets) over experimental rows'})
        if name == 'lysine':
            for mark, label in LYSINE_MARKS:
                t = tables[mark]
                u = mark_units(t)
                mids = site_ids(t)
                rows.append({'set': f'{name}/{mark}', 'label': f'  of which {label}', 'sites': len(mids), 'high_confidence': int((u >= 2).sum()),
                             'with_interaction_evidence': len(ev & mids), 'high_confidence_rule': 'D5, units pooled within this table'})
    return rows


def interaction_section(out: Path, version: str, previous: str, config: dict) -> str:
    """The interaction part of the README, from scripts/64's release notes (headings one level down, its title and
    first paragraph replaced by this README's own)."""
    x64 = _load('export64', '64_export_site_interactions.py')
    with tempfile.TemporaryDirectory() as tmp:
        t = Path(tmp)
        (t / 'supporting').mkdir()
        for f in sorted((out / 'supporting').glob(f'*_{version}.tsv')):
            if any(f.name.endswith(f'_{tb}_{version}.tsv') for tb in INTERACTION_TABLES) and not f.name.startswith('site_changes_'):
                (t / 'supporting' / f.name).symlink_to(f)
        arg = [n for n, s in config['site_sets'].items() if set(s['residues']) <= {'R'}]
        x64.write_release_notes(t, version, config, previous if arg else '')
        text = (t / 'README.md').read_text(encoding='utf-8')
    body = text[text.index('The tables are made by'):]
    body = re.sub(r'^(#+) ', lambda m: '#' + m.group(1) + ' ', body, flags=re.M)
    return body.replace(f'`releases/{previous}/data_dictionary.tsv` and its README describe them',
                        '`data_dictionary.tsv` and the methylarginine section below describe them')


def carried_readme(prev: Path, previous: str, version: str, arginine_sets: list[str]) -> str:
    text = (prev / 'README.md').read_text(encoding='utf-8')
    text = text[text.index('\n## ') + 1:]
    text = re.sub(rf'_{re.escape(previous)}(?=\.(?:tsv|xlsx))', f'_{version}', text)
    text = text.replace(f'site_changes_v1.4.3_to_{version}.tsv', f'site_changes_{previous}_to_{version}.tsv')
    text = re.sub(r'(?<![\w])(' + '|'.join(INTERACTION_TABLES) + rf')_{re.escape(version)}\.tsv', rf'{arginine_sets[0]}_\1_{version}.tsv', text)
    text = text.replace(f'`_{previous}` in place of', f'`_{version}` in place of').replace(f'renamed `_{previous}`', f'renamed `_{version}`')
    return text


def dictionary(prev: Path, out: Path, previous: str, version: str, config: dict, arginine_sets: list[str]) -> str:
    x_src = _load('sources_vocab', 'interaction_sources.py')
    lines = (prev / 'data_dictionary.tsv').read_text(encoding='utf-8').rstrip('\n').split('\n')
    header, body = lines[0], lines[1:]
    out_lines, arg_rows = [header], {}
    for line in body:
        name, rest = line.split('\t', 1)
        if name.startswith('supporting/site_changes_'):
            name = f'supporting/site_changes_{previous}_to_{version}.tsv'
        else:
            name = new_name(name, previous, version, arginine_sets)
        out_lines.append(name + '\t' + rest)
        for t in INTERACTION_TABLES:
            if name == f'supporting/{arginine_sets[0]}_{t}_{version}.tsv':
                arg_rows.setdefault(t, []).append(rest.split('\t'))
    for name, spec in config['site_sets'].items():
        if name in arginine_sets:
            continue
        letters = sorted(spec['residues'])
        for t in INTERACTION_TABLES:
            path = out / 'supporting' / f'{name}_{t}_{version}.tsv'
            if not path.exists():
                continue
            cols = path.read_text(encoding='utf-8').split('\n', 1)[0].split('\t')
            by_col = {r[0]: r for r in arg_rows.get(t, [])}
            for col in cols:
                r = list(by_col.get(col, [col, 'string', '', '', '']))
                if col == 'record_class':
                    r[2], r[3] = '; '.join(x_src.SITE_SET_CLASSES), '; '.join(f'{k} = {v}' for k, v in x_src.SITE_SET_CLASSES.items())
                elif col == 'site_state':
                    voc = x_src.site_set_state_vocabulary(config, letters)
                    r[2], r[3] = '; '.join(voc), '; '.join(f'{k} = {v}' for k, v in voc.items())
                elif col == 'evidence_type':
                    r[2] = '; '.join(sorted(set(read(path)['evidence_type'])))
                    r[3] = ''
                out_lines.append(f'supporting/{name}_{t}_{version}.tsv\t' + '\t'.join(r))
    return '\n'.join(out_lines) + '\n'


def code_commit() -> str:
    paths = ['scripts/65_carry_release.py', 'scripts/64_export_site_interactions.py', 'scripts/56_build_context_marks.py', 'config/site_sets.json']
    try:
        dirty = subprocess.run(['git', 'status', '--porcelain', '--', *paths], capture_output=True, text=True, check=True, cwd=ROOT).stdout
        sha = subprocess.run(['git', 'log', '-1', '--format=%H', '--', *paths], capture_output=True, text=True, check=True, cwd=ROOT).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return 'unknown'
    return sha + (' (scripts modified)' if dirty.strip() else '')


def counts_table(rows: list[dict]) -> list[str]:
    lines = ['| Site set | Sites | High-confidence sites | Sites with interaction evidence |', '|---|---:|---:|---:|']
    for r in rows:
        ev = f'{r["with_interaction_evidence"]:,}'
        if 'with_generic_interaction_evidence' in r:
            ev += f' ({r["with_generic_interaction_evidence"]:,} from the kinds of evidence of the other sets)'
        lines.append(f'| {r["label"].strip() if not r["label"].startswith("  ") else "&nbsp;&nbsp;" + r["label"].strip()} | {r["sites"]:,} | {r["high_confidence"]:,} | {ev} |')
    return lines


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--version', required=True)
    ap.add_argument('--previous', required=True)
    ap.add_argument('--release-date', required=True)
    ap.add_argument('--config', default=str(ROOT / 'config' / 'site_sets.json'))
    args = ap.parse_args()
    version, previous = args.version, args.previous
    prev, out = ROOT / 'releases' / previous, ROOT / 'releases' / version
    config = json.loads(Path(args.config).read_text(encoding='utf-8'))
    arginine_sets = [n for n, s in config['site_sets'].items() if set(s['residues']) <= {'R'}]
    missing = [f'{n}_site_interactions_{version}.tsv' for n in config['site_sets'] if not (out / 'supporting' / f'{n}_site_interactions_{version}.tsv').exists()]
    if missing:
        raise SystemExit(f'interaction tables missing in {out}/supporting (run scripts/64 first): {missing}')
    x56 = _load('build56', '56_build_context_marks.py')

    carried, edited = carry(prev, out, previous, version, arginine_sets)
    sites = read(out / f'human_methylarginine_sites_{version}.tsv')
    prev_changes = next((prev / 'supporting').glob(f'site_changes_*_to_{previous}.tsv'))
    changes_name = f'supporting/site_changes_{previous}_to_{version}.tsv'
    x56.site_changes(sites, read(prev_changes)).to_csv(out / changes_name, sep='\t', index=False, lineterminator='\n')
    (out / 'data_dictionary.tsv').write_text(dictionary(prev, out, previous, version, config, arginine_sets), encoding='utf-8')

    rows = count_marks(out, version, previous, config)
    prev_release = json.loads((prev / 'release.json').read_text(encoding='utf-8'))
    interaction_files = {f.relative_to(out).as_posix(): md5_file(f) for f in sorted((out / 'supporting').glob(f'*_{version}.tsv'))
                         if any(f.name.endswith(f'_{t}_{version}.tsv') for t in INTERACTION_TABLES)}
    new_sets = [n for n in config['site_sets'] if n not in arginine_sets]
    commit = code_commit()

    # README
    R = [f'# Human methylarginine sites, context marks and site-level interaction evidence, {version}', '',
         f'{version} holds everything of {previous}, the methylarginine data of v1.4.3 and the tables of context marks, under this '
         f'release\'s names and otherwise unchanged, and adds site-level interaction tables for the context marks '
         f'({", ".join(new_sets)}), made the way the methylarginine release makes its own. The methylarginine interaction tables are '
         f'here too, carried over byte for byte, so one release holds the tables of every site set. {previous} was not released '
         'on its own. See CHANGES.md for what changed and why.', '',
         '## Counts', '',
         'Sites are the distinct residues of each site set\'s tables. High-confidence sites are the website\'s control, not a tier '
         'of the release: the release tables keep every row. For the context marks a site is high confidence when at least two '
         'independent studies report it from experimental records: each PubMed ID counts once, a direct dataset (Berryhill 2023, '
         'Ochoa 2020, Hendriks 2017, Akimov 2018) counts as one study, and so as the paper a database cites, and UniProt and iPTMnet '
         'count only through the papers they cite; a residue\'s studies are pooled over all its rows. For methylarginine it is the '
         'existing control: two or more studies, a laboratory\'s studies counted once, PhosphoSitePlus, UniProt, IEDB and '
         'ProMetheusDB not counted. Sites with interaction evidence have a row of their set\'s `site_interactions` table that counts '
         'as evidence; the methylarginine tables hold more kinds of evidence (predicted interfaces and contacts, IntAct, RBDmap, '
         'pCLAP), so their count from only the kinds the other sets have is given too.', '',
         *counts_table(rows), '',
         '## Site-level interaction evidence', '', interaction_section(out, version, previous, config).rstrip('\n'), '',
         carried_readme(prev, previous, version, arginine_sets).rstrip('\n'), '',
         '## Checksums', '', 'Run `md5sum -c MD5SUMS` in this folder, and in `analyses/` and `analyses_with_lower_confidence/`.']
    (out / 'README.md').write_text(americanize('\n'.join(R)) + '\n', encoding='utf-8')

    # CHANGES
    x56_counts = prev_release.get('context_marks', {})
    C = [f'# Changes from {previous} to {version}', '',
         f'{version} is the release; {previous} was prepared but not released on its own. Everything of {previous} is in {version} '
         f'under this release\'s names, its data unchanged, and {version} adds site-level interaction tables for the context marks.', '',
         '## New: site-level interaction evidence for the context marks', '',
         f'For {", ".join(new_sets)}, `supporting/<site set>_<table>_{version}.tsv` gives what the methylarginine release gives for '
         'its sites: PDB interfaces with each entry\'s method, paper and the state of the residue it models, RNA and DNA contacts '
         'within 4 Å in structures with the bound sequences, proximal and nearby nucleic-acid contacts, and RNA-binding evidence '
         'from cells (`scripts/64_export_site_interactions.py`, the steps of `scripts/37` with the modified residue as a parameter).', '',
         *counts_table(rows), '',
         '## The methylarginine interaction tables, carried over', '',
         f'The six methylarginine interaction tables of {previous} (and of v1.4.3) are carried over byte for byte under the names '
         'the other site sets use, so that every set\'s tables sit side by side:', '',
         '| In ' + previous + ' | In ' + version + ' |', '|---|---|']
    C += [f'| `supporting/{t}_{previous}.tsv` | `supporting/{arginine_sets[0]}_{t}_{version}.tsv` |' for t in INTERACTION_TABLES]
    C += ['', '## Files renamed', '',
          f'Every carried file is renamed `_{previous}` to `_{version}`, its content unchanged except where it names its version: '
          + ', '.join(f'`{e}`' for e in edited) + '. '
          f'`{changes_name}` replaces {previous}\'s site changes table; it lists every one of the {len(sites):,} methylarginine sites as kept.', '',
          f'## Before the release: CPLM 4.0 taken out of {previous}', '',
          f'CPLM 4.0 states no license ("All Rights Reserved") and compiles sites from other databases as well as from papers. The owner '
          f'decided on 17 November 2026 to leave it out, and {previous}, not yet released, was rebuilt without it (pull request #27). '
          'Ochoa 2020, Berryhill 2023, Hendriks 2017 and Akimov 2018 stay. The context-mark tables of this release:', '',
          '| Table | Rows with CPLM | Rows in this release |', '|---|---:|---:|']
    before = {'methyllysine': 10889, 'lysine_acylation': 88816, 'acetylation': 53177, 'ubiquitination': 125771, 'sumoylation': 48593, 'phosphorylation': 142415}
    for m, n in before.items():
        now = x56_counts.get(m, {}).get('rows')
        C.append(f'| `human_{m}_sites_{version}.tsv` | {n:,} | {now:,} |')
    C += ['', 'The site-level interaction tables of methyllysine and lysine were rebuilt from these tables; each holds exactly the rows '
          'of the sites still in the set.', '',
          f'## Fixed: the analyses checksums of {previous}', '',
          f'`analyses/MD5SUMS` and `analyses_with_lower_confidence/MD5SUMS` listed their README.md with the checksum of the v1.4.3 '
          f'text, while the README names the release\'s own version; they are checksummed again here (and in {previous}).', '',
          'Code: `scripts/64_export_site_interactions.py` (interaction tables), `scripts/56_build_context_marks.py` (context marks), '
          f'`scripts/65_carry_release.py` (this release), at commit `{commit}`.']
    (out / 'CHANGES.md').write_text(americanize('\n'.join(C)) + '\n', encoding='utf-8')

    release = {'title': 'Human methylarginine sites, context marks and site-level interaction evidence', 'version': version,
               'previous_version': previous, 'release_date': args.release_date, 'repository': 'https://github.com/Shechterlab/PTM_Pipeline',
               'license': prev_release.get('license'), 'methylarginine_data': 'unchanged from v1.4.3: files carried over byte for byte, not rebuilt',
               'methylarginine_code_commit': prev_release.get('methylarginine_code_commit'),
               'context_marks_code_commit': prev_release.get('context_marks_code_commit'), 'release_code_commit': commit,
               'carried_over': carried, 'carried_over_text_edited': edited, 'site_changes': changes_name,
               'sources': prev_release.get('sources'), 'sources_excluded': prev_release.get('sources_excluded'),
               'context_marks': prev_release.get('context_marks'), 'interaction_tables': interaction_files, 'counts': rows}
    (out / 'release.json').write_text(json.dumps(release, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    names = sorted(p.relative_to(out).as_posix() for p in out.rglob('*') if p.is_file() and p.name != 'MD5SUMS')
    (out / 'MD5SUMS').write_text(''.join(f'{md5_file(out / n)}  {n}\n' for n in names), encoding='utf-8')
    print(json.dumps(rows, indent=1))


if __name__ == '__main__':
    main()

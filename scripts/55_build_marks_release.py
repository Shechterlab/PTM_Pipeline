#!/usr/bin/env python3
"""Build a marks release: one open table per mark from the staged tables of inputs/marks/.

    python scripts/55_build_marks_release.py --version v0.1.0 --out releases/marks/v0.1.0

For every mark (methyllysine, phosphorylation) it concatenates the staged tables named `*_<mark>_sites.tsv` in
inputs/marks/ (one per source: UniProt now, others as they are staged), sorted, into
`human_<mark>_sites_<version>.tsv`, and writes the UniProt methylarginine annotations as
`human_methylarginine_uniprot_crosscheck_<version>.tsv` with a column that says whether the methylarginine release
(releases/v1.4.3) has the arginine as a site, a lower-confidence site, or not. The methylarginine sites themselves stay in
that release. Also written: `data_dictionary.tsv`, `README.md` (sources, licenses, counts by state and evidence
class, how to cite), `release.json` and `MD5SUMS`. Nothing in the output depends on the date, so a second run on the
same inputs gives the same bytes.

The release is CC BY 4.0, which every source allows. A source whose terms are not CC BY (or as open) is not added here.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from common import md5_file

COLUMNS = ['canonical_UniProtAC', 'position', 'residue', 'mark', 'state', 'evidence', 'evidence_class', 'pmids',
           'source', 'source_release']
MARKS = ['methyllysine', 'phosphorylation']
CHECKED_RELEASE = 'v1.4.3'
# Terms of each source name allowed in the table. A source not listed here stops the build.
SOURCES = {
    'UniProt': {'license': 'CC BY 4.0', 'citation': 'UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. '
                'Nucleic Acids Res 53:D609-D617 (2025). doi:10.1093/nar/gkae1010',
                'url': 'https://www.uniprot.org/help/license'},
    # Berryhill2023 is added here only once the license of its supplements is confirmed (inputs/marks/README.md).
}
DICTIONARY = [
    ('canonical_UniProtAC', 'string', '', 'reviewed UniProtKB accession of the canonical sequence in the reference proteome (UniProt 2026_03, UP000005640)'),
    ('position', 'integer', '', 'residue number on that sequence, 1-based'),
    ('residue', 'string', 'K, S, T, Y, R', 'one-letter code of the residue; checked against the reference sequence'),
    ('mark', 'string', 'methyllysine, phosphorylation, methylarginine', 'the modification'),
    ('state', 'string', 'methyllysine: me1, me2, me3, me (degree not stated); phosphorylation: phospho; methylarginine: me1, me2a, me2s, me2 (type not stated), me (degree not stated)',
     'the state reported for the residue; a residue with several states has one row per state'),
    ('evidence', 'string', 'ECO codes, ";"-separated', 'ECO evidence codes of the source record (UniProt: ECO:0000269 published experiment, ECO:0007744 combinatorial computational and experimental evidence, ECO:0000250 by similarity, ECO:0000255, ECO:0000303, ECO:0000305)'),
    ('evidence_class', 'string', 'experimental, by-similarity, high-throughput',
     'experimental: ECO:0000269 or ECO:0007744 among the codes. by-similarity: UniProt records with no such code (ECO:0000250 and others). high-throughput: a proteomics study that is not a UniProt record (none in v0.1.0). A site can hide by-similarity rows by filtering on this column'),
    ('pmids', 'string', 'PubMed IDs, ";"-separated', 'PubMed IDs of the source record, when it gives them'),
    ('source', 'string', 'UniProt', 'database or study the row comes from'),
    ('source_release', 'string', '', 'release of the source (UniProt release)'),
    ('in_methylarginine_release', 'string', 'site, lower_confidence, no',
     f'crosscheck table only: whether the methylarginine release {CHECKED_RELEASE} has the arginine as a site, as a lower-confidence site, or not at all'),
]


def read_stage(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


def combine(frames: list[pd.DataFrame]) -> pd.DataFrame:
    table = pd.concat(frames, ignore_index=True)[COLUMNS]
    table['_pos'] = table['position'].astype(int)
    table = table.sort_values(['canonical_UniProtAC', '_pos', 'state', 'source']).drop(columns='_pos')
    return table.reset_index(drop=True)


def crosscheck(stage: pd.DataFrame, release_dir: Path) -> pd.DataFrame:
    """The UniProt methylarginine rows with whether the methylarginine release has the arginine."""
    def keys(name: str) -> set[tuple[str, str]]:
        t = pd.read_csv(release_dir / name, sep='\t', dtype=str, keep_default_na=False, usecols=['uniprot_accession', 'position'])
        return set(zip(t['uniprot_accession'], t['position'], strict=True))
    sites = keys(f'human_methylarginine_sites_{CHECKED_RELEASE}.tsv')
    lower = keys(f'human_methylarginine_lower_confidence_sites_{CHECKED_RELEASE}.tsv')
    out = stage.copy()
    out['in_methylarginine_release'] = [
        'site' if k in sites else 'lower_confidence' if k in lower else 'no'
        for k in zip(out['canonical_UniProtAC'], out['position'], strict=True)]
    return out


def counts_table(table: pd.DataFrame) -> list[str]:
    lines = ['| State | ' + ' | '.join(sorted(table['evidence_class'].unique())) + ' | Rows |',
             '|---|' + '---|' * (table['evidence_class'].nunique() + 1)]
    classes = sorted(table['evidence_class'].unique())
    for state, grp in table.groupby('state'):
        lines.append(f'| {state} | ' + ' | '.join(f'{(grp["evidence_class"] == c).sum():,}' for c in classes) + f' | {len(grp):,} |')
    lines.append('| all | ' + ' | '.join(f'{(table["evidence_class"] == c).sum():,}' for c in classes) + f' | {len(table):,} |')
    return lines


def write_readme(path: Path, version: str, tables: dict[str, pd.DataFrame], cross: pd.DataFrame, files: dict[str, str],
                 summary: dict) -> None:
    sources = sorted({s for t in tables.values() for s in t['source'].unique()})
    rel = sorted({r for t in tables.values() for r in t['source_release'].unique()})
    lines = [f'# Human PTM sites by mark, {version}', '',
             'One table per modification of the human reviewed proteome, from openly licensed sources only, for the sites '
             'of the marks family (each with its own questions) to read. The methylarginine sites are not here: they are in '
             f'[`releases/{CHECKED_RELEASE}/`](../../{CHECKED_RELEASE}/), whose license is unchanged. This is a first version (v0): '
             'UniProt annotations only.', '',
             'Positions are on UniProt release 2026_03 (reviewed human reference proteome UP000005640, canonical accessions), '
             f'the reference of every release of this repository. The annotations are from UniProt release {", ".join(rel)}; an entry '
             'whose sequence differs from the reference is left out, and so is any annotation whose residue is not the expected '
             f'one ({summary["dropped_total"]}).', '',
             '## Files', '', '| File | One row per |', '|---|---|']
    for mark in MARKS:
        lines.append(f'| `{files[mark]}` | {mark} site, state and source ({len(tables[mark]):,} rows on '
                     f'{tables[mark][["canonical_UniProtAC", "position"]].drop_duplicates().shape[0]:,} residues of '
                     f'{tables[mark]["canonical_UniProtAC"].nunique():,} proteins) |')
    lines += [f'| `{files["methylarginine"]}` | UniProt arginine methylation annotation ({len(cross):,} rows), with whether the methylarginine release has the arginine. A cross-check, not the methylarginine table |',
              '| `data_dictionary.tsv` | column of each table |', '| `release.json` | sources, inputs and their MD5 |',
              '| `MD5SUMS` | checksums of the files |', '',
              '## Sources and licenses', '', '| Source | License | Used for |', '|---|---|---|']
    for s in sources:
        lines.append(f'| {s} ([terms]({SOURCES[s]["url"]})) | {SOURCES[s]["license"]} | '
                     + ', '.join(m for m in MARKS if s in set(tables[m]['source'])) + ' |')
    lines += ['', '**License of this release: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).** Every source above is '
              'CC BY 4.0, which allows redistribution and adaptation under CC BY 4.0 with attribution, so the release does not '
              'need the share-alike and non-commercial terms of the methylarginine releases. That is a choice for this '
              'folder only: the methylarginine release stays CC BY-NC-SA 4.0 (see `LICENSE.md`). A source whose terms are not '
              'CC BY, or as open, is not added to this release. Excluded by design: CPLM, PhosphoSitePlus, iPTMnet, dbPTM and qPTM '
              '(all rights reserved, non-commercial or unclear terms).', '',
              '## Evidence classes', '',
              '- `experimental`: a UniProt record with ECO:0000269 (published experiment) or ECO:0007744 (combinatorial computational '
              'and experimental evidence: large-scale proteomics, which is most of UniProt\'s phosphorylation).',
              '- `by-similarity`: the other UniProt records, mainly ECO:0000250 (by similarity to another species\' entry), also '
              'ECO:0000255, 0000303 and 0000305, and records with no code. Filter these out to keep measured sites.',
              '- `high-throughput`: a proteomics study that is not a UniProt record. None in this version.', '']
    for mark in MARKS:
        lines += [f'## Counts: {mark}', ''] + counts_table(tables[mark]) + ['']
    lines += ['## Counts: methylarginine cross-check', ''] + counts_table(cross) + ['',
              'Arginines of the cross-check table in the methylarginine release ' + CHECKED_RELEASE + ': '
              + ', '.join(f'{k} {v:,}' for k, v in cross.drop_duplicates(['canonical_UniProtAC', 'position'])['in_methylarginine_release'].value_counts().sort_index().items()) + '.', '',
              '## How to cite', '',
              f'DeAngelo JD, Silverstein AM, Shechter D. PTM_Pipeline: human PTM sites by mark, {version}. '
              'https://github.com/Shechterlab/PTM_Pipeline (Zenodo: [10.5281/zenodo.22696798](https://doi.org/10.5281/zenodo.22696798)), '
              'CC BY 4.0, and the sources you use:', '']
    lines += [f'- {SOURCES[s]["citation"]}' for s in sources]
    lines += ['', '## Rebuilding', '',
              '```bash',
              'python scripts/54_stage_uniprot_ptm_sites.py --fasta inputs/reference/uniprot_human_reviewed_2026_03.fasta',
              f'python scripts/55_build_marks_release.py --version {version} --out releases/marks/{version}',
              '```', '',
              'The UniProt API serves its current release only; `--raw` of script 54 re-stages from a kept download (`inputs/marks/raw_files.tsv` '
              'lists its MD5). A published version\'s folder is never rewritten: a change is a new version.', '']
    path.write_text('\n'.join(lines), encoding='utf-8')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--version', required=True, help='for example v0.1.0')
    ap.add_argument('--out', required=True)
    ap.add_argument('--inputs', default='inputs/marks')
    ap.add_argument('--methylarginine-release', default=f'releases/{CHECKED_RELEASE}')
    args = ap.parse_args()
    inputs, out = Path(args.inputs), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    tables, used = {}, {}
    for mark in MARKS:
        paths = sorted(inputs.glob(f'*_{mark}_sites.tsv'))
        if not paths:
            raise SystemExit(f'no staged table for {mark} in {inputs}')
        frames = [read_stage(p) for p in paths]
        for frame in frames:
            unknown = set(frame['source']) - set(SOURCES)
            if unknown:
                raise SystemExit(f'source {sorted(unknown)} has no license entry in SOURCES: check its terms first')
        tables[mark] = combine(frames)
        used[mark] = {p.name: md5_file(p) for p in paths}
    cross_stage = sorted(inputs.glob('*_methylarginine_crosscheck.tsv'))
    if len(cross_stage) != 1:
        raise SystemExit(f'expected one methylarginine crosscheck table in {inputs}, found {len(cross_stage)}')
    cross = crosscheck(combine([read_stage(cross_stage[0])]), Path(args.methylarginine_release))
    used['methylarginine'] = {cross_stage[0].name: md5_file(cross_stage[0])}

    files = {m: f'human_{m}_sites_{args.version}.tsv' for m in MARKS}
    files['methylarginine'] = f'human_methylarginine_uniprot_crosscheck_{args.version}.tsv'
    for mark in MARKS:
        tables[mark].to_csv(out / files[mark], sep='\t', index=False, lineterminator='\n')
    cross.to_csv(out / files['methylarginine'], sep='\t', index=False, lineterminator='\n')

    stage_summary = json.loads((inputs / 'uniprot_ptm_sites_summary.json').read_text(encoding='utf-8'))
    dropped = stage_summary['dropped']
    dropped_text = ', '.join(f'{v:,} {k.replace("_", " ")}' for k, v in sorted(dropped.items())) or 'none dropped'
    summary = {'dropped_total': f'dropped in staging: {dropped_text}'}

    dictionary = ['file\tcolumn\ttype\tvalues\tdescription']
    for name in files.values():
        for column, kind, values, description in DICTIONARY:
            if column == 'in_methylarginine_release' and 'crosscheck' not in name:
                continue
            dictionary.append(f'{name}\t{column}\t{kind}\t{values}\t{description}')
    (out / 'data_dictionary.tsv').write_text('\n'.join(dictionary) + '\n', encoding='utf-8')
    write_readme(out / 'README.md', args.version, tables, cross, files, summary)

    release = {
        'version': args.version, 'license': 'CC BY 4.0',
        'reference': 'UniProt 2026_03, reviewed human reference proteome UP000005640 (inputs/reference/)',
        'sources': {s: SOURCES[s]['license'] for s in sorted({s for t in tables.values() for s in t['source'].unique()})},
        'inputs_md5': used,
        'rows': {m: len(t) for m, t in tables.items()} | {'methylarginine_crosscheck': len(cross)},
        'methylarginine_release_checked': CHECKED_RELEASE,
    }
    (out / 'release.json').write_text(json.dumps(release, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    names = sorted(p.name for p in out.iterdir() if p.is_file() and p.name != 'MD5SUMS')
    (out / 'MD5SUMS').write_text(''.join(f'{md5_file(out / n)}  {n}\n' for n in names), encoding='utf-8')
    print(json.dumps(release['rows']))


if __name__ == '__main__':
    main()

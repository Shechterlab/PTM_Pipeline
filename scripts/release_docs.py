"""Documents of a release (called by scripts/build_release.py): data_dictionary.tsv, release.json, README.md,
CHANGES.md and MD5SUMS. Counts are computed from the tables just written, and the column sections come from the
column definitions in scripts/release_codes.py, so the documents cannot drift from the data.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import pandas as pd

import release_codes as rc
from disease_families import ADDITIONAL_FAMILIES, AUDITED_FAMILIES

REPO = Path(__file__).resolve().parent.parent
GITHUB = 'https://github.com/Shechterlab/PTM_Pipeline'
PAPER = ('DeAngelo JD, Silverstein AM, Shechter D. Protein arginine methylation rewires molecular interactions. '
         'Cell Chem Biol 33:1404-1421 (2026), doi:10.1016/j.chembiol.2026.09.011')
FEATURED = {'FUS': 'P35637', 'EWSR1': 'Q01844', 'FMR1': 'Q06787', 'Histone H3.1': 'P68431', 'SmB/SNRPB': 'P14678',
            'SmD1/SNRPD1': 'P62314', 'SmD3/SNRPD3': 'P62318'}


def split_key(site: str) -> tuple[str, int]:
    acc, pos = site.split(':R')
    return acc, int(pos)


def clean_int(x) -> int | None:
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def counts(sites, records, changes, studies, maron_counts, prev_sites) -> dict:
    pending_rows = sites['flags'].str.contains('one_sided_window')
    placed = sites[~pending_rows]
    flags = sites['flags'].str.split('; ').explode()
    out = {
        'rows': len(sites), 'sites': len(placed), 'proteins': placed['uniprot_accession'].nunique(),
        'pending': int(pending_rows.sum()),
        'proteins_all_rows': sites['uniprot_accession'].nunique(),
        'records': len(records), 'records_placed': int((records['placed_site'] != '').sum()),
        'previous_sites': len(prev_sites), 'previous_proteins': prev_sites['uniprot_accession'].nunique(),
        'previous_records_placed': int((records['previous_site'] != '').sum()),
        'kept': int((changes['change'] == 'kept').sum()), 'removed': int((changes['change'] == 'removed').sum()),
        'added': int((changes['change'] == 'added').sum()),
        'renumbered': int(changes['reason_codes'].str.contains('reference_renumbered').sum()),
        'kept_pending': int(changes.loc[changes['change'] == 'kept', 'site_id'].isin(sites.loc[pending_rows, 'site_id']).sum()),
        'kept_records_changed': int((changes['records_change'] == 'changed').sum()),
        'kept_records_replaced': int((changes['records_change'] == 'replaced').sum()),
        'flags': {k: int(v) for k, v in flags[flags != ''].value_counts().items()},
        'studies': len(studies),
        'with_maron_s3': int((sites['maron2021_s3_antibodies'] != '').sum()),
        'with_state': int((placed['methyl_state'] != '').sum()),
        'state': {k: int(v) for k, v in placed['methyl_state'].replace('', 'none').value_counts().items()},
        'state_evidence': {k: int(v) for k, v in placed['methyl_state_evidence'].str.split('; ').explode().replace('', pd.NA).dropna().value_counts().items()},
        'no_domain': int((placed['domain'] == 'no_domain_annotation').sum()),
        'peptide_start_records': int((records['placement_code'] == 'sequence_peptide_start').sum()),
    }
    out.update({f'maron_{k}': v for k, v in maron_counts.items()})
    # Domain context of kept sites, and the proteins the published InterPro staging left without intervals.
    before = prev_sites.set_index('site_id')['domain']
    previous_of = changes.loc[changes['change'] == 'kept'].set_index('site_id')['previous_site_id']
    now = placed.set_index('site_id')['domain']
    kept = [s for s in previous_of.index if s in now.index]
    out['domain_changed_kept'] = int(sum(before[previous_of[s]] != now[s] for s in kept))
    by_protein = prev_sites.groupby('uniprot_accession')['domain'].agg(set)
    unannotated = {a for a, d in by_protein.items() if d == {'no_domain_annotation'}}
    annotated_now = set(now.index.map(lambda s: split_key(s)[0])[now != 'no_domain_annotation'])
    out['previously_unannotated_now_annotated'] = len(unannotated & annotated_now)
    return out


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def commit() -> str:
    try:
        return subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True, text=True, cwd=REPO).stdout.strip() or 'unknown'
    except OSError:
        return 'unknown'


def dictionary(files: dict) -> pd.DataFrame:
    rows = []
    for name, _frame, columns in files.values():
        for column, kind, codes, description in columns:
            meanings = rc.CODES.get(codes, {})
            rows.append({'file': name, 'column': column, 'type': kind, 'values': '; '.join(meanings),
                         'value_meanings': '; '.join(f'{c} = {m}' for c, m in meanings.items()), 'description': description})
    return pd.DataFrame(rows)


def n(value) -> str:
    return f'{value:,}' if isinstance(value, int) else str(value)


def write_all(outdir: Path, ctx: dict) -> None:
    files, c, version, previous = ctx['files'], ctx['counts'], ctx['version'], ctx['previous']
    dictionary_columns = [('file', 'string', '', 'the release file, relative to the release folder'),
                          ('column', 'string', '', 'the column'),
                          ('type', 'code', '', 'integer, number, string, boolean, code, or a list of these'),
                          ('values', 'list<string>', '', 'the codes a coded column may hold'),
                          ('value_meanings', 'list<string>', '', 'code = meaning, for each code'),
                          ('description', 'string', '', 'what the column holds')]
    files = dict(files)
    files['dictionary'] = ('data_dictionary.tsv', None, dictionary_columns)
    table = dictionary(files)
    table.to_csv(outdir / 'data_dictionary.tsv', sep='\t', index=False, lineterminator='\n')
    check_codes(outdir, files)
    write_release_json(outdir, ctx, files)
    write_readme(outdir, ctx, files)
    write_changes(outdir, ctx, files)
    write_md5sums(outdir)


def check_codes(outdir: Path, files: dict) -> None:
    """Every value of a coded column is a code of the dictionary."""
    for name, _frame, columns in files.values():
        data = pd.read_csv(outdir / name, sep='\t', dtype=str, keep_default_na=False)
        for column, kind, codes, _ in columns:
            if not codes:
                continue
            values = data[column].str.split('; ').explode() if kind.startswith('list') else data[column]
            unknown = set(values) - set(rc.CODES[codes]) - {''}
            if unknown:
                raise SystemExit(f'{name}: column {column} holds values that are not codes: {sorted(unknown)[:5]}')


def write_release_json(outdir: Path, ctx: dict, files: dict) -> None:
    manifest = pd.read_csv(ctx['inputs'] / 'MANIFEST.tsv', sep='\t', dtype=str, keep_default_na=False)
    record = {
        'title': 'Human methylarginine sites',
        'version': ctx['version'],
        'previous_version': ctx['previous'],
        'code_commit': commit(),
        'repository': GITHUB,
        'release_date': ctx['release_date'],
        'paper': PAPER,
        'doi': None,
        'licence': 'CC BY-NC-SA 4.0 (LICENSE.md); sources keep their own terms (studies file, inputs/MANIFEST.tsv)',
        'counts': {k: v for k, v in ctx['counts'].items()},
        'validation': ctx['validation'],
        'reference': {'proteome': 'UP000005640', 'reviewed': True, **{k: ctx['uniprot_release'][k] for k in ('release', 'release_date', 'entries', 'fasta', 'fasta_md5')}},
        'annotation': {'InterPro': '108.0', 'MobiDB': 'export of 2026-03-30', 'domain_ontology': 'config/domain_class_ontology.json'},
        'phosphositeplus': {'downloaded': '2021-05-22', 'via': 'Maron et al. 2021 Table S5',
                            'terms': 'CC BY-NC-SA 3.0 Unported, inputs/licences/phosphositeplus_terms_2021-07-11_web-archive.pdf'},
        'files': [],
        'inputs': [{'file': f'inputs/{r.file}', 'md5': r.md5, 'version_or_date': r.version_or_date, 'licence': r.licence,
                    'redistribution': r.redistribution} for r in manifest.itertuples()],
    }
    for name, *_ in files.values():
        path = outdir / name
        with open(path, encoding='utf-8') as fh:
            rows = sum(1 for _ in fh) - 1
        record['files'].append({'file': name, 'rows': rows, 'md5': md5(path)})
    xlsx = outdir / f"human_methylarginine_sites_{ctx['version']}.xlsx"
    record['files'].append({'file': xlsx.name, 'rows': ctx['counts']['rows'], 'md5': md5(xlsx)})
    (outdir / 'release.json').write_text(json.dumps(record, indent=1, ensure_ascii=False) + '\n')


# What one row of each release file is, in the order the README lists them.
FILE_ROWS = {
    'proteins': 'protein: names, site counts, functional classes, arginine enrichment, site density and clusters, condensates, locations, disease families, other PTMs, interactions',
    'studies': 'study or database: citation, DOI, terms, counts',
    'statistics': 'enrichment test behind the summaries: classes, motifs, disorder, domains, disease families, condensates, locations, interactions',
    'records': 'source record: what the study reported, its evidence, where each version placed it',
    'changes': 'site of either version: kept, removed or added, and why',
    'maron_detection': 'site detected in Maron et al. 2021 Table S2 or S3 (per antibody)',
    'maron_inhibitors': 'site, antibody and PRMT inhibitor in Maron et al. 2021 Table S3',
    'carm1_ko': 'site and cell line in the CARM1 knockout data of Shishkova et al. 2017',
    'site_interactions': 'piece of interaction evidence at a site: crosslinks, PDB interfaces, predicted interfaces and contacts, IntAct records',
    'residue_evidence': 'run of residues of a protein with the same interaction evidence',
    'protein_partners': 'protein and STRING physical partner (score >= 0.7)',
    'complexes': 'Complex Portal or hu.MAP complex with a methylated member',
    'disease_terms': 'enriched DisGeNET term (DOSE enrichDGN)',
    'disease_families': 'disease family: terms and genes',
    'density_windows': 'protein and window width with significantly dense sites (q <= 0.05)',
    'condensates': 'CD-CODE condensate: name, type, members and methylated members',
    'hpa_locations': 'protein and Human Protein Atlas location',
    'other_ptm_sites': 'other modification on a protein of the release (PTM master table)',
}


def column_section(name: str, columns: list[tuple]) -> list[str]:
    lines = ['', f'### `{name}`', '', '| Column | Type | Description |', '|---|---|---|']
    for column, kind, codes, description in columns:
        codes_text = f" Codes: {', '.join(f'`{k}`' for k in rc.CODES[codes])}." if codes and len(rc.CODES[codes]) <= 14 else (' Codes: see `data_dictionary.tsv`.' if codes else '')
        lines.append(f'| `{column}` | {kind.replace("<", "&lt;").replace(">", "&gt;")} | {description}{"." if description and not description.endswith(".") else ""}{codes_text} |')
    return lines


def write_readme(outdir: Path, ctx: dict, files: dict) -> None:
    c, version, previous = ctx['counts'], ctx['version'], ctx['previous']
    studies = ctx['studies']
    ur = ctx['uniprot_release']
    lines = [
        f'# Human methylarginine sites, {version}',
        '',
        f"{n(c['sites'])} sites on {n(c['proteins'])} proteins, and {n(c['pending'])} pending sites, "
        f"{n(c['rows'])} rows: the human protein arginine methylation sites of {PAPER}, with corrections to how some "
        f'source records were placed on UniProt sequences. [CHANGES.md](CHANGES.md) says what changed since {previous} '
        f'(the published version) and why.',
        '',
        f'Every site has its sequence, methylation evidence, structural context, the studies that report it, and the '
        f'source records behind it. Positions are on UniProt release {ur["release"]} (reviewed human reference proteome '
        f'UP000005640).',
        '',
        '## Files',
        '',
        '| File | One row per |',
        '|---|---|',
        f"| `{files['sites'][0]}` | site, including the {n(c['pending'])} sites kept pending (flag `one_sided_window`) |",
        f"| `human_methylarginine_sites_{version}.xlsx` | the same table as a spreadsheet |",
    ]
    rows = dict(FILE_ROWS, sites=None)
    lines += [f"| `{files[k][0]}` | {text} |" for k, text in rows.items() if text and k in files]
    lines += [
        '| `data_dictionary.tsv` | column of every table, with its codes and their meanings |',
        '| `release.json` | the release: version, counts, reference, inputs, files with MD5 |',
        '| `CHANGES.md` | what changed since the previous version |',
        '| `MD5SUMS` | checksums (`md5sum -c MD5SUMS`) |',
        '',
        'Tables are UTF-8, tab-separated, one header row. Lists are separated by `; `, booleans are `true`/`false`, '
        'and a missing value is empty. Coded columns hold codes only; `data_dictionary.tsv` gives their meanings.',
        '',
        '## Cite',
        '',
        f'Cite {PAPER}, and the studies of the sites you use (`studies` column, `{files["studies"][0]}`):',
        '',
        '| Study or database | Citation | DOI | Records | Sites | Only here |',
        '|---|---|---|---:|---:|---:|',
    ]
    lines += [f'| {s.label}{"" if s.study_id in s.label else f" (`{s.study_id}`)"} | {s.citation} | [{s.doi}](https://doi.org/{s.doi}) | {s.n_records:,} | {s.n_sites:,} | {s.n_sites_only_here:,} |'
              for s in studies.itertuples()]
    lines += [
        '',
        'The site records come from three source tables: Maron et al. 2021 Table S5 (sites compiled from ten studies), '
        'ProMetheusDB (Massignani et al. 2022) and the arginine methylation rows of an iPTMnet-format PTM master table '
        '(UniProt and IEDB records). Maron et al. 2021 Tables S2 and S3 and the CARM1 knockout data (Shishkova et al. '
        '2017, doi:10.1038/ncomms15571) annotate sites but add none.',
        '',
        '## Licence',
        '',
        'The release is CC BY-NC-SA 4.0 ([LICENSE.md](../../LICENSE.md)); each study keeps its own terms (`licence` in '
        'the studies file). The release includes PhosphoSitePlus® (PSP, www.phosphosite.org) records compiled in Maron et '
        'al. 2021 Table S5. PhosphoSite, created by Cell Signaling Technology, was licensed under CC BY-NC-SA 3.0 '
        'Unported in 2021, when those records were downloaded '
        '([inputs/licences/](../../inputs/licences/phosphositeplus_terms_2021-07-11_web-archive.pdf)). Those terms ask '
        'that PSP be cited (Hornbeck et al., Nucleic Acids Res 43:D512-D520, 2015), that bulk downloads list PSP as the '
        'source (the `studies` column does), and that web pages showing PSP-derived sites link them to the PSP protein '
        'page (www.phosphosite.org/uniprotAccAction?id=<UniProt accession>).',
        '',
        'The annotation and interaction sources keep their own terms too (`licence` in `inputs/MANIFEST.tsv`). Three are '
        'non-commercial: the DisGeNET terms behind the disease tables (CC BY-NC-SA 4.0), CD-CODE condensates (CC BY-NC 4.0) '
        'and the pepR-MS RNA crosslink compendium (CC BY-NC).',
        '',
        '## How it was built',
        '',
        f"`bash hpc/build_release.sh {version}` at commit `{commit()}` of {GITHUB}:",
        '',
        '1. checks every input against `inputs/MANIFEST.tsv`;',
        '2. parses the source tables and places each record on the reference sequence '
        '(`scripts/02_remap_arg_methyl_to_canonical.py`);',
        '3. stops if any record is placed against its own evidence (`scripts/qc_remap_evidence_check.py`; the 22 '
        'accepted paralog rows of `config/remap_evidence_allowlist.tsv` aside) or any site is not an arginine;',
        '4. runs the analyses on the union: functional classes (`scripts/03`), motifs and arginine enrichment (`08`), '
        'disorder and domain context (`06`, `05`, `07`), other PTMs (`09`), site density and clusters (`10`, `28`, `25`), '
        'condensates (`12`), Maron et al. 2021 summaries (`16`) and disease enrichment (`14`);',
        '5. encodes the release for the website (`scripts/36`) and computes the interaction evidence of its sites from '
        '`inputs/interactions/` (`scripts/37`);',
        '6. writes this folder (`scripts/build_release.py`).',
        '',
        f'- **Reference sequences:** UniProt release {ur["release"]} ({ur["release_date"]}), the {ur["entries"]:,} reviewed '
        f'entries of the human reference proteome UP000005640 (`inputs/reference/`). `reference_sequence_md5` tells '
        'whether a later UniProt entry still has the same numbering.',
        '- **Context:** MobiDB predicted disorder (export of 2026-03-30), InterPro 108.0 domains and repeats, sequence '
        'motifs: the classifiers of `scripts/06`, `scripts/05`, `scripts/07` and `scripts/08`. The intervals were '
        'computed on earlier UniProt sequences and moved to these (`inputs/annotation/interval_lift_to_uniprot_2026_03.tsv`).',
        f"- **Previous version:** `releases/{previous}/` (the published union, git tag {previous}).",
        '',
        '## Columns',
        '',
        'The same definitions, with every code and its meaning, are in `data_dictionary.tsv`.',
    ]
    for k in ['sites'] + list(FILE_ROWS):
        if k in files:
            lines += column_section(files[k][0], files[k][2])
    (outdir / 'README.md').write_text('\n'.join(lines) + '\n')


def write_changes(outdir: Path, ctx: dict, files: dict) -> None:
    c, version, previous = ctx['counts'], ctx['version'], ctx['previous']
    records, changes, sites = ctx['records'], ctx['changes'], ctx['sites']
    reason_records = records['reason_code'].replace('', pd.NA).dropna().value_counts()
    code_sites = changes.assign(code=changes['reason_codes'].str.split('; ')).explode('code')
    gene_of = sites.set_index('site_id')['gene']
    started = sorted(set(records.loc[records['placement_code'] == 'sequence_peptide_start', 'placed_site']), key=lambda s: (gene_of.get(s, ''), split_key(s)[1]))
    start_sites = ', '.join(f"{gene_of.get(s, '') or split_key(s)[0]} R{split_key(s)[1]}" for s in started)
    agreement = (f"all {c['start_rule_agree']:,} placements the pipeline made from plain peptides by other means"
                 if not c['start_rule_disagree'] else
                 f"{c['start_rule_agree']:,} placements made from plain peptides by other means, and disagrees with {c['start_rule_disagree']:,}")
    pending = sites[sites['flags'].str.contains('one_sided_window')]
    gene_or_acc = {r.site_id: (r.gene.split(';')[0] or r.uniprot_accession) for r in pending.itertuples()}
    shifted = [f"{gene_or_acc[r.pending_site]} R{split_key(r.pending_site)[1]} (reported at R{int(float(r.reported_position))})"
               for r in records[records['pending_site'] != ''].itertuples()
               if clean_int(r.reported_position) not in (None, split_key(r.pending_site)[1])]
    corrected = records[records['source_correction'] != '']
    nan_rows = corrected[corrected['source_correction'] == 'nan_restored']
    na_rows = corrected[corrected['source_correction'] == 'accession_na']
    lines = [
        f'# Changes from {previous} to {version}',
        '',
        f'{version} uses the sources and methods of DeAngelo et al. (2026), mapped to UniProt release '
        f"{ctx['uniprot_release']['release']}, with corrections to how some source records were placed on UniProt "
        f'sequences. {previous} is the published union (git tag {previous}, `releases/{previous}/`).',
        '',
        f'| | {previous} | {version} |',
        '|---|---:|---:|',
        f"| Sites | {n(c['previous_sites'])} | {n(c['sites'])}, and {n(c['pending'])} pending |",
        f"| Proteins | {n(c['previous_proteins'])} | {n(c['proteins'])} |",
        f"| Source records placed | {n(c['previous_records_placed'])} of {n(c['records'])} | {n(c['records_placed'])} of {n(c['records'])} |",
        '',
        f"Of the {n(c['previous_sites'])} published sites, {n(c['kept'])} are kept (sites of {version} too, "
        f"{n(c['kept_pending'])} of them pending) and {n(c['removed'])} were removed; {n(c['added'])} sites were added. "
        'A site moves when its records now point at another arginine: the old one is removed and the new one added. '
        f"A kept site can still have gained or lost records: {n(c['kept_records_changed'])} did, and at "
        f"{n(c['kept_records_replaced'])} none of the previous records remain, other records now naming the same "
        'arginine (supporting/site_changes, column records_change).'
        + (f" One kept site has a new number only: RAB12 R10 is R106 in UniProt {ctx['uniprot_release']['release']}, which adds 96 residues at its N-terminus." if c['renumbered'] == 1 else ''),
        '',
        '## What was corrected',
        '',
        'Each record whose placement changed carries one reason code.',
        '',
        '| Reason code | What changed | Records | Sites removed | Sites added |',
        '|---|---|---:|---:|---:|',
    ]
    for code, text in rc.CODES['reason'].items():
        r = int(reason_records.get(code, 0))
        removed = int(((code_sites['code'] == code) & (code_sites['change'] == 'removed')).sum())
        added = int(((code_sites['code'] == code) & (code_sites['change'] == 'added')).sum())
        if (r or removed or added) and code not in ('reference_renumbered', 'one_sided_window'):
            lines.append(f"| `{code}` | {text} | {r:,} | {removed:,} | {added:,} |")
    lines += [
        '',
        'A site can carry several codes when its records changed for different reasons.',
        '',
        '## Corrections to the source tables',
        '',
        f"- **Windows missing NAN** ({len(nan_rows)} records, `source_correction` `nan_restored`). Maron et al. 2021 Table "
        'S5 gives these Larsen 2016 windows 28 residues long, where every other Larsen window has 31, and no Larsen window '
        'there holds the tripeptide NAN: the letters were read as a missing value when the table was made. The parser '
        'restores NAN where inserting it at exactly one place gives the window centred on the reported position of the '
        f"reported entry; the window as the table gives it is kept in `evidence_as_reported`. {previous} placed "
        f"{int((nan_rows['previous_site'] != '').sum())} of these records; this version places "
        f"{int((nan_rows['placed_site'] != '').sum())}.",
        f"- **Records without an accession** ({len(na_rows)} records, `source_correction` `accession_na`). Table S5 gives "
        f"them a gene symbol and a window, and the accession NA; {previous} left them out. They are placed from their gene "
        f"symbol and window ({int((na_rows['placed_site'] != '').sum())} placed) and numbered after the records of "
        f'{previous}, so that every other record keeps its record_id.',
        '',
        '## Proteins featured in the paper',
        '',
        f'| Protein | Sites in {previous} | Sites in {version} | Removed | Added |',
        '|---|---:|---:|---|---|',
    ]
    pending_ids = set(pending['site_id'])
    for name, acc in FEATURED.items():
        g = changes[changes['uniprot_accession'] == acc]
        mark = lambda s: s.split(':')[1] + (' (pending)' if s in pending_ids else '')  # noqa: E731
        removed = ', '.join(map(mark, g.loc[g['change'] == 'removed', 'site_id'])) or '-'
        added = ', '.join(map(mark, g.loc[g['change'] == 'added', 'site_id'])) or '-'
        now = g[g['in_this_version'].astype(bool)]
        n_pending = int(now['site_id'].isin(pending_ids).sum())
        lines.append(f"| {name} ({acc}) | {int(g['in_previous_version'].sum())} | {len(now) - n_pending}{f', and {n_pending} pending' if n_pending else ''} | {removed} | {added} |")
    lines += [
        '',
        '## Decided in this version',
        '',
        '- **Peptides with several arginines.** A peptide does not say which of its arginines is methylated, and the '
        'position the source reports may follow an older numbering. Where Maron et al. 2021 Tables S2 and S3 give the '
        "peptide's start in the same numbering, the reported position is counted from that start (reason "
        f"`peptide_start`): {c['peptide_start_records']} records ({int(reason_records.get('peptide_start', 0))} of them placed differently in {previous}), on {start_sites}. Without such a start the record "
        f'stays unplaced (reason `multi_arginine_peptide`). Where a start is known, counting from it agrees with {agreement}.',
        f"- **One UniProt reference, release {ctx['uniprot_release']['release']}.** {previous} used sequences staged on "
        '2026-04-07 whose UniProt release was not recorded. Between them, 13 reviewed entries changed sequence; two carry '
        'sites: POLR2A (one substitution, outside every window; its sites keep their numbers) and RAB12, whose site R10 '
        'is R106 now.',
        '- **Domain context for every protein.** The published analyses used InterPro intervals staged for only part of '
        'the union, so some proteins had no domain annotation. This version uses InterPro 108.0 for all reviewed '
        f"entries: {c['previously_unannotated_now_annotated']:,} proteins whose published sites were all \"no domain "
        f"annotation\" now have a domain or repeat, and the domain context of {c['domain_changed_kept']:,} kept sites changes.",
        '- **Protein classes from PANTHER.** The paper\'s functional classes are keyword classes that overlap and leave '
        'much of the proteome out. This version adds the PANTHER 19.0 protein classes (scripts/03b): each protein\'s '
        'top-level class, grouped for the overview (config/panther_class_groups.json) so that the groups and the proteins '
        'PANTHER leaves unclassified partition the reviewed proteome, with tests of the groups, the 24 top-level classes '
        'and their subclasses (statistics: panther_group, panther_class, panther_subclass). The paper\'s classes are kept.',
        f'- **Disease families.** The ALS / motor neuron / FTD rule now also takes TDP-43 proteinopathy, Pick disease '
        '(not Niemann-Pick) and anterior horn cell disease, which the published rule missed (5 DisGeNET terms, among '
        'them the enriched terms TDP-43 Proteinopathies, Pick Disease of the Brain and Degeneration of anterior horn '
        f'cells); the other {len(AUDITED_FAMILIES) - 1} audited families are unchanged. The '
        f'DisGeNET terms none of them takes now go to {len(ADDITIONAL_FAMILIES)} additional, rule-based families '
        '(scripts/disease_families.py ADDITIONAL_FAMILY_RULES: mostly clinical signs, by organ system), tested and '
        'BH-adjusted apart (statistics: disease_family_additional). They are not audited.',
        f"- **Methylation degree at each site** (`methyl_state`, {n(c['with_state'])} of {n(c['sites'])} sites). "
        f"{previous} gave a degree only from ProMetheusDB. This version also counts Maron et al. 2021 Table S2, UniProt's "
        'modified-residue annotations with experimental evidence (inferred ones are listed apart) and the two studies that '
        'looked for one degree only: Larsen 2016 (monomethylation) and Wei 2020 (ADMA antibody, dimethylation). Mass '
        'spectrometry tells mono- from dimethylation, but asymmetric from symmetric mostly by the antibody used, so '
        '`methyl_state` gives Rme1, Rme2 or both, and the columns it is drawn from keep the detail: '
        + ', '.join(f"{n(v)} {k.replace('_', ' + ')}" for k, v in c['state'].items() if k != 'none') + '.',
        f'- **PRMT inhibitor calls** (Maron et al. 2021 Table S3). A peptide changes when its adjusted p value is at most '
        f'{rc.PADJ_CUTOFF}. Each site, antibody and inhibitor gets a call over its peptides, and each site and inhibitor a '
        f'call over its antibodies (`{files["maron_inhibitors"][0]}`). The same combination gives the CARM1 knockout call '
        'over the two cell lines.',
        '',
        '## Still open',
        '',
        f"- **Windows that match on one side only** ({n(c['pending'])} sites, flag `one_sided_window`; "
        f"{int(pending['in_previous_version'].sum())} of them in {previous}): the source window matches the entry exactly on "
        'one side of the arginine and not on the other, and no other reviewed entry holds it. Three of these windows occur '
        'in full in another UniProt isoform of the protein (GANAB R185, TCF7L2 R126, SNX12 R153; checked 2026-10-20), so '
        'they come from a splice variant. '
        + (('Some are reported at another position, in the numbering of an isoform, and match one complete flank (15 '
            'residues and the arginine) here: ' + ', '.join(shifted) + '. ') if shifted else '')
        + 'The pipeline leaves these records unplaced; the release keeps the sites, flagged and left out of every '
        'analysis, until the policy is decided: '
        + ', '.join(f"{g} R{split_key(s)[1]}" for g, s in zip(pending['gene'], pending['site_id'])) + '.',
        f"- **Paralog windows** ({c['flags'].get('paralog_window', 0)} sites, flag `paralog_window`): a window or "
        'peptide matches another reviewed entry exactly and the named entry within 2 residues. The other entry is a '
        'paralog, except at RBM39 R109, where the 6-residue peptide SALRGK reads SAIRGK in RBM39 (I and L have the same '
        'mass) and also occurs in CCDC13. The records stay on the entry the source named, accepted in '
        '`config/remap_evidence_allowlist.tsv`.',
        f"- **Arginine picked by the reported position** ({c['flags'].get('arginine_from_position', 0)} sites, flag "
        '`arginine_from_position`): every record of the site has a peptide with several arginines, or a window or peptide '
        'found more than once in the entry, so the evidence alone does not single out this arginine; the positions the '
        'sources report do.',
        '',
        '## Tracing a change',
        '',
        f"- `{files['changes'][0]}`: one row per site of either version, with its change, reason codes, the records "
        'involved and, for a removed site, where those records are now (for an added site, where they were; for a kept '
        'site, the records that moved onto or off it, and both).',
        f"- `{files['records'][0]}`: one row per source record, with what the study reported, the evidence sequence, where "
        'each version placed it, the independent evidence check and the reason code.',
        '',
        'Their columns are listed in `README.md` and `data_dictionary.tsv`.',
    ]
    (outdir / 'CHANGES.md').write_text('\n'.join(lines) + '\n')


def write_md5sums(outdir: Path) -> None:
    with open(outdir / 'MD5SUMS', 'w') as fh:
        for p in sorted(outdir.rglob('*')):
            if p.is_file() and p.name != 'MD5SUMS':
                fh.write(f'{md5(p)}  {p.relative_to(outdir).as_posix()}\n')

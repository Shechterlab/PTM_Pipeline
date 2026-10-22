"""Documents of a release (called by scripts/build_release.py): data_dictionary.tsv, release.json, README.md,
CHANGES.md and MD5SUMS. Counts are computed from the tables just written, and the column sections come from the
column definitions in scripts/release_codes.py, so the documents cannot drift from the data.
"""
from __future__ import annotations

import hashlib
import json
import re
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
    """HEAD, with '+uncommitted' when tracked files differ from it (the release was built from a working tree)."""
    try:
        head = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True, text=True, cwd=REPO).stdout.strip() or 'unknown'
        dirty = subprocess.run(['git', 'status', '--porcelain', '--untracked-files=no'], capture_output=True, text=True, cwd=REPO).stdout.strip()
        return head + ('+uncommitted' if dirty else '')
    except OSError:
        return 'unknown'


def write_plain_tsv(table: pd.DataFrame, path: Path) -> None:
    """A table of strings as plain TSV, no quoting (the values hold no tab or newline)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as fh:
        fh.write('\t'.join(table.columns) + '\n')
        for row in table.itertuples(index=False):
            fh.write('\t'.join(row) + '\n')


def dictionary(files: dict) -> pd.DataFrame:
    rows = []
    for name, _frame, columns in files.values():
        for column, kind, codes, description in columns:
            meanings = rc.CODES.get(codes, {})
            if any(';' in m for m in meanings.values()):
                raise SystemExit(f'a meaning of the {codes} codes holds ";", the separator of value_meanings')
            rows.append({'file': name, 'column': column, 'type': kind, 'values': '; '.join(meanings),
                         'value_meanings': '; '.join(f'{c} = {m}' for c, m in meanings.items()), 'description': description})
    return pd.DataFrame(rows)


def n(value) -> str:
    return f'{value:,}' if isinstance(value, int) else str(value)


def write_all(outdir: Path, ctx: dict) -> None:
    files, previous = ctx['files'], ctx['previous']
    dictionary_columns = [('file', 'string', '', 'the release file, relative to the release folder'),
                          ('column', 'string', '', 'the column'),
                          ('type', 'code', '', 'integer, number, string, boolean, code, or a list of these'),
                          ('values', 'list<string>', '', 'the codes a coded column may hold'),
                          ('value_meanings', 'list<string>', '', 'code = meaning, for each code'),
                          ('description', 'string', '', 'what the column holds')]
    files = dict(files)
    files['dictionary'] = ('data_dictionary.tsv', None, dictionary_columns)
    table = dictionary(files)
    if table.apply(lambda col: col.str.contains('[\t\n\r]', regex=True)).any().any():
        raise SystemExit('data_dictionary.tsv: a tab or newline inside a value')
    write_plain_tsv(table, outdir / 'data_dictionary.tsv')
    check_codes(outdir, files)
    write_release_json(outdir, ctx, files)
    write_readme(outdir, ctx, files)
    if previous == 'v1.0.0':
        write_changes(outdir, ctx, files)
    else:
        write_changes_since_release(outdir, ctx, files)
    write_studies(outdir, ctx)
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
        # The tools of the enrichment tests (scripts/enrich_sets.R), as the analyses recorded them.
        'enrichment_tools': ctx.get('enrichment_tools', {}),
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
    'lower_confidence_sites': 'lower-confidence site: an arginine a study reports as methylated, or an earlier release had as a site, whose records each fail a site test, with the results of the tests',
    'studies': 'study or database: citation, DOI, terms, counts',
    'statistics': 'enrichment test behind the summaries: classes, motifs, disorder, domains, disease families, condensates, locations, interactions',
    'statistics_with_lower_confidence': 'the same tests, run again with the lower-confidence sites whose own evidence covers the arginine and leaves the methyl group on an arginine',
    'records': 'source record: what the study reported, its evidence, where each version placed it',
    'lower_confidence': 'source record that adds no site (not localised, or not on one site of a protein it names), with where its methyl group could be',
    'changes': 'site of either version: kept, removed or added, and why',
    'published_sites': 'site of v1.0.0 or v1.1.0: what it is in this release (a site, a lower-confidence site, or not on the reference)',
    'maron_detection': 'site detected in Maron et al. 2021 Table S2 or S3 (per antibody)',
    'maron_inhibitors': 'site, antibody and PRMT inhibitor in Maron et al. 2021 Table S3',
    'carm1_ko': 'site and cell line in the CARM1 knockout data of Shishkova et al. 2017',
    'site_interactions': 'piece of interaction evidence at a site: crosslinks, PDB interfaces, predicted interfaces and contacts, IntAct records',
    'residue_evidence': 'run of residues of a protein with the same interaction evidence',
    'protein_partners': 'protein and STRING physical partner (score >= 0.7)',
    'complexes': 'Complex Portal or hu.MAP complex with a methylated member',
    'nucleic_acid_contacts': 'site and RNA or DNA chain within 4 Å of its side chain in a PDB structure: the contacted nucleotides and the bound sequence around them',
    'nucleic_acid_proximal_sites': 'site, complex and kind of nucleic acid in which the arginine side chain is within 8 Å, not 4 Å, of RNA or DNA: the bound sequence around the closest nucleotide',
    'nucleic_acid_near_sites': 'site and complex in which a residue at i ± 1-5, not the arginine side chain, is within 4 Å of RNA or DNA (the side chain not within 8 Å): the bound sequence there',
    'nucleic_acid_regions': 'run of residues of a protein within 4 Å of RNA or DNA in PDB structures: the bound segments and the methylarginines inside it',
    'rna_bind_n_seq': 'methylated protein and ENCODE RNA Bind-n-Seq experiment: its most enriched 5-mers and logos (in vitro)',
    'dna_motifs': 'methylated protein and JASPAR 2026 DNA-binding profile',
    'eclip': 'methylated protein and ENCODE eCLIP experiment: peaks by transcript region and the genes with most peaks (in cells)',
    'disease_terms': 'enriched DisGeNET term (DOSE enrichDGN)',
    'disease_terms_tested': 'DisGeNET term DOSE enrichDGN tested, enriched or not',
    'disease_ontology': 'Disease Ontology term tested (clusterProfiler enricher on the gene sets of DOSE enrichDO)',
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
    added = [label for sid, label, *_ in rc.STUDIES if sid in set(rc.STUDIES_ADDED) and sid in set(ctx['records']['study'])]
    localised = 'localised_by' in ctx['records'].columns and ctx['records']['placement_code'].isin(rc.PLACEMENT_NOT_LOCALISED).any()
    lines = [
        f'# Human methylarginine sites, {version}',
        '',
        f"{n(c['sites'])} sites on {n(c['proteins'])} proteins"
        + (f", and {n(c['pending'])} pending sites" if c['pending'] else '') + f", {n(c['rows'])} rows: the human protein arginine "
        f'methylation sites of {PAPER}, with corrections to how some source records were placed on UniProt sequences (v1.1.0)'
        + ('; from v1.2.0 every study is read from its own tables and a site is kept only where its methyl group is '
           'localised to one arginine' if localised else '')
        + (f" and the sites of {len(added)} more studies ({', '.join(added)})" if added else '')
        + f'. [CHANGES.md](CHANGES.md) says what changed since {previous}'
        + (' (the published version)' if previous == 'v1.0.0' else '') + ' and why.',
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
        '| `STUDIES.md` | each study and database behind the sites: where its records come from, how they were filtered and placed, and records traced from the source to the site |',
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
        'The site records come from Maron et al. 2021 Table S5 (sites compiled from ten studies), '
        'ProMetheusDB (Massignani et al. 2022) and the arginine methylation rows of an iPTMnet-format PTM master table '
        '(UniProt and IEDB records)'
        + (f", and the supplementary tables of {', '.join(added)}, each staged into one table of records "
           '(`scripts/00_stage_study_sites.py`, `inputs/sources/`)' if added else '')
        + '. Maron et al. 2021 Tables S2 and S3 and the CARM1 knockout data (Shishkova et al. '
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
        'The annotation and interaction sources keep their own terms too (`licence` in `inputs/MANIFEST.tsv`). Four are '
        'non-commercial: the DisGeNET terms behind the disease tables (CC BY-NC-SA 4.0), CD-CODE condensates (CC BY-NC 4.0), '
        'the pepR-MS RNA crosslink compendium (CC BY-NC) and the DNA crosslinks of Trendel et al. 2025 (CC BY-NC 4.0). Maron '
        'et al. 2021 and ProMetheusDB supplements are CC BY-NC-ND 4.0: the release redistributes records derived from them, '
        'not the tables themselves.',
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
        f"- **Previous version:** `releases/{previous}/` (git tag {previous}).",
        '',
        '## Statistics',
        '',
        f"`{files['statistics'][0]}` holds every test behind the summaries: its analysis and category, the 2x2 counts a, b, "
        'c and d, the effect, p, q and the method:',
        '',
        '- **Protein and gene sets** (method `enricher`): PANTHER and functional classes, CD-CODE condensates, Human Protein '
        'Atlas locations, disease families and Gene Ontology RNA and DNA binding; the complexes, Disease Ontology terms and '
        'DisGeNET terms, tested the same way, are in `supporting/`. Over-representation analysis in R with clusterProfiler::enricher, and '
        'DOSE::enrichDGN for the DisGeNET terms (`scripts/enrich_sets.R`; '
        + ', '.join(f'{k} {v}' for k, v in ctx.get('enrichment_tools', {}).items()) + '). '
        'For a set holding M of the N members of the background, where n of the N are methylarginine-containing proteins '
        '(or their genes) and k of those are in the set, p is the one-sided hypergeometric probability of k or more, '
        '`fold_enrichment` is (k/n)/(M/N), and q is the Benjamini-Hochberg adjusted p over the sets of the analysis; '
        'a = k, b = n - k, c = M - k and d = N - M - n + k. A set is tested when it holds a methylarginine-containing '
        'protein and its size is within the analysis\'s range. Backgrounds: the reviewed human proteome (classes, Gene '
        'Ontology, complexes), the proteins of any human CD-CODE condensate (condensates), the proteins the Human Protein '
        'Atlas profiled (locations), the reviewed human genes with a DisGeNET term (disease families and DisGeNET terms) '
        'and those with a Disease Ontology term (Disease Ontology).',
        '- **Arginines** (methods `fisher` and `mantel_haenszel`): motifs, disorder, domain context and the interaction '
        'tests compare methylated arginines with the other arginines of the same proteins: two-sided Fisher exact test '
        '(scipy) with the sample odds ratio after adding 0.5 to each cell, and for the interaction tests also the '
        'Mantel-Haenszel common odds ratio with each protein a stratum (statsmodels); q is Benjamini-Hochberg within the '
        'analysis and method.',
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
        '- **Protein-set tests.** Every test of protein or gene sets (PANTHER and functional classes, CD-CODE condensates, '
        'Human Protein Atlas locations, disease families, Gene Ontology binding, complexes) is an over-representation '
        'analysis with clusterProfiler::enricher, and the DisGeNET terms with DOSE::enrichDGN (scripts/enrich_sets.R; '
        'statistics: method enricher, with fold enrichment): p is the one-sided hypergeometric probability, fold enrichment '
        '(k/n)/(M/N), q Benjamini-Hochberg over the sets of the analysis (README, Statistics). The disease families are tested '
        'against the reviewed human genes with a DisGeNET term. The Human Disease Ontology is added, with the gene sets of '
        'DOSE enrichDO (supporting/disease_ontology), and every DisGeNET term tested is listed (supporting/disease_terms_tested).',
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


def write_changes_since_release(outdir: Path, ctx: dict, files: dict) -> None:
    """CHANGES.md from a previous release of this format (v1.1.0 on): the studies added, what happened to the previous
    sites and records, and what changed in the analyses."""
    c, version, previous = ctx['counts'], ctx['version'], ctx['previous']
    records, changes, sites, studies = ctx['records'], ctx['changes'], ctx['sites'], ctx['studies']
    placed = sites[~sites['flags'].str.contains('one_sided_window')]
    prev = ctx['previous_sites']
    prev_pending = prev['flags'].str.contains('one_sided_window')
    old = records[records['change_since_previous_version'] != 'new_record']
    new = records[records['change_since_previous_version'] == 'new_record']
    added_ids = set(changes.loc[changes['change'] == 'added', 'site_id'])
    prev_proteins = set(changes.loc[changes['change'] != 'added', 'site_id'].map(lambda x: split_key(x)[0]))
    first_new = int(new['record_id'].min()) if len(new) else None
    # what v1.2.0 changed (studies read from their own tables, site tests, corrections) is told against v1.1.0 only
    since_v11 = previous == 'v1.1.0'
    lines = [
        f'# Changes from {previous} to {version}',
        '',
        f'| | {previous} | {version} |',
        '|---|---:|---:|',
        f"| Sites | {n(int((~prev_pending).sum()))} | {n(c['sites'])} |",
        f"| Proteins | {n(prev.loc[~prev_pending, 'uniprot_accession'].nunique())} | {n(c['proteins'])} |",
        f"| Sites kept pending | {n(int(prev_pending.sum()))} | {n(c['pending'])} |",
        f"| Source records | {n(len(old))} | {n(len(records))} |",
        f"| Source records placed | {n(int((old['previous_site'] != '').sum()))} | {n(c['records_placed'])} |",
        f"| Studies and databases | {n(len(old['study'].unique()))} | {n(len(studies))} |",
        '',
        f"Of the {n(c['previous_sites'])} sites of {previous} (pending ones included), {n(c['kept'])} are kept and {n(c['removed'])} removed; "
        f"{n(c['added'])} sites were added. The records of {previous} are read in the same order and "
        + ('all are placed as before. ' if (old['change_since_previous_version'].isin(['same', 'not_placed'])).all() else
           f"{n(int((~old['change_since_previous_version'].isin(['same', 'not_placed'])).sum()))} of them are placed differently "
           '(supporting/source_records, columns change_since_previous_version and reason_code). ')
        + (f'The records read for the first time (each study\'s own table) follow them, from record_id {first_new:,}.' if first_new else ''),
        '',
    ]
    # ---- Site specificity: the records without localisation, and the sites of the previous version they alone held
    if since_v11 and 'localised_by' in records.columns and records['placement_code'].isin(rc.PLACEMENT_NOT_LOCALISED).any():
        code_n = records['placement_code'].value_counts()
        unloc = records[records['placement_code'] == 'unresolved_not_localised']
        removed = changes[changes['change'] == 'removed']
        prev_studies = prev.set_index('site_id')['studies'].str.split('; ')
        rem_studies = removed['previous_site_id'].map(prev_studies)

        def support(st) -> str:
            others = set(st) - {'PSP'}
            if not others:
                return 'PhosphoSitePlus only'
            if len(others) == 1:
                return dict((x[0], x[1]) for x in rc.STUDIES)[next(iter(others))] + ' (with or without PhosphoSitePlus)'
            if others == {'Larsen2016', 'Maron2021'}:
                return 'Larsen 2016 and Maron 2021 (with or without PhosphoSitePlus)'
            return 'two or more studies'
        by_support = rem_studies.dropna().map(support).value_counts()
        localised_elsewhere = records[(records['localised_by'] != '') & records['placement_code'].isin(rc.PLACEMENT_SEVERAL | rc.PLACEMENT_OTHER_PROTEIN | {'unresolved_study_not_found'})]
        not_unique_sites = set(localised_elsewhere['_previous_at']) & set(removed['previous_site_id'])
        lines += [
            '## Site specificity',
            '',
            'A record adds a site only where its source localises the methyl group to one arginine: a localisation score over '
            "the study's threshold (MaxQuant localisation probability, ptmRS, Ascore or ModScore, ModLS), a peptide whose mass "
            'leaves one arrangement (no other residue the search could methylate, or every such residue a dimethylated '
            'arginine), a UniProt annotation with experimental evidence from the literature, or a PhosphoSitePlus site, '
            'trusted as curated (`localised_by` in `supporting/source_records`). Its peptide, or without one its window, must '
            'also name exactly one arginine in the reviewed proteome, in a protein the record names.',
            '',
            'Every study is read from its own tables (`inputs/sources/<study>_methylarginine_records.tsv`, staged by '
            '`scripts/00_stage_study_sites.py`), so the rows of the compiled tables give way to them: '
            f"{n(int(code_n.get('superseded_by_study_table', 0)))} records of Maron et al. 2021 Table S5 and the PTM master table "
            'are not placed (`superseded_by_study_table`), the study\'s own records standing in for them. '
            f"{n(len(unloc))} records carry no localisation at all (`unresolved_not_localised`): "
            + ', '.join(f"{n(v)} {k}" for k, v in unloc['study'].map(dict((x[0], x[1]) for x in rc.STUDIES)).value_counts().items())
            + '. All records stay in `supporting/source_records`.',
            '',
            f"The {n(len(removed))} sites of {previous} that no localised record supports are removed. By the sources that reported "
            f'them in {previous}:',
            '',
            '| Reported in ' + previous + ' by | Sites removed |',
            '|---|---:|',
        ] + [f'| {k} | {v:,} |' for k, v in by_support.items()] + [
            '',
            f"{n(len(not_unique_sites))} of them have a localised record whose peptide or window names several arginines in "
            'the reviewed proteome (paralogs, repeats) or one in a protein the record does not name.',
            '',
            "The others were reported through the compiled tables by studies whose own tables do not localise them: the "
            "study's table gives the site below its localisation threshold, gives no score for a peptide that leaves "
            "several residues the search could methylate, or does not list it. Their records stay in "
            "`supporting/source_records` (`superseded_by_study_table`, or `unresolved_not_localised` for ProMetheusDB).",
            '',
        ]
        # ---- Lower-confidence sites: the arginines reported or published that no record passes every test for
        lc_sites = files['lower_confidence_sites'][1] if 'lower_confidence_sites' in files else None
        published = files['published_sites'][1] if 'published_sites' in files else None
        if lc_sites is not None and len(lc_sites):
            analysed = int(lc_sites['in_analyses_with_lower_confidence'].astype(str).str.lower().eq('true').sum())
            failed = lc_sites['tests_failed'].str.split('; ').explode()
            lines += [
                '## Lower-confidence sites',
                '',
                f"{n(len(lc_sites))} arginines that a study reports as methylated, or that v1.0.0 or v1.1.0 had as sites, are not sites "
                'of this release: each of their records fails at least one site test. They are kept as lower-confidence sites '
                f"(`human_methylarginine_lower_confidence_sites_{ctx['version']}.tsv`), with the results of the record that passes most "
                "tests, and are never in this release's counts and tests. The analyses were run a second time with the "
                f"{n(analysed)} whose own evidence covers the arginine and leaves the methyl group on an arginine "
                f"(`statistics_with_lower_confidence_{ctx['version']}.tsv`, `analyses_with_lower_confidence/`). The website shows "
                'them, hatched, when the reader switches Lower-confidence sites on, and its counts and charts then come from that '
                'second run.',
                '',
                "Every record has passed its study's identification filter (false discovery rate or score) before it is read. It "
                'is then put to five tests at the arginine its study reports; a site needs a record that passes all five.',
                '',
                '| Test | Passes when | Lower-confidence sites whose best record fails it |',
                '|---|---|---:|',
            ] + [f"| `{t}` | {rc.CODES['site_test'][t]} | {n(int((failed == t).sum()))} |" for t in rc.CODES['site_test']] + ['']
            if published is not None and len(published):
                status = list(rc.CODES['published_status'])
                tab = pd.crosstab(published['version'], published['status']).reindex(columns=status, fill_value=0)
                lines += [
                    'What became of every site of the earlier releases (`supporting/published_sites`):',
                    '',
                    '| Release | Sites | ' + ' | '.join(f'`{x}`' for x in status) + ' |',
                    '|---|---:|' + '---:|' * len(status),
                ] + [f'| {v} | {n(int(row.sum()))} | ' + ' | '.join(n(int(row[x])) for x in status) + ' |' for v, row in tab.iterrows()] + [
                    '',
                    '`' + '`, `'.join(status) + '`: ' + '; '.join(rc.CODES['published_status'][x] for x in status) + '.',
                    '',
                ]
        # ---- Lower-confidence evidence: the records that add no site but say where a methyl group could be
        ambiguous = rc.PLACEMENT_SEVERAL | rc.PLACEMENT_OTHER_PROTEIN | {'unresolved_study_not_found'}
        lower = records[(records['placed_site'] == '') & (records['pending_site'] == '')
                        & (records['placement_code'] != 'superseded_by_study_table')
                        & ((records['localised_by'] == '') | records['placement_code'].isin(ambiguous))]
        if 'candidate_sites' in lower.columns and len(lower):
            label = {**{c: 'several_sites' for c in rc.PLACEMENT_SEVERAL}, **{c: 'other_protein' for c in rc.PLACEMENT_OTHER_PROTEIN},
                     'unresolved_study_not_found': 'not_in_reference'}
            why = [label.get(c, '') if loc else (nl or 'several_residues')
                   for loc, nl, c in zip(lower['localised_by'], lower['not_localised'], lower['placement_code'])]
            counts = pd.Series(why).value_counts()
            lines += [
                '## Lower-confidence evidence',
                '',
                f"{n(len(lower))} records add no site but say where a methyl group could be: they are not localised, or localised "
                'but not to one site of a protein they name. They are listed with their candidate sites in '
                '`supporting/lower_confidence_evidence`, never counted as sites or used in a test, and the website shows them '
                'only on request (protein map, row Lower-confidence evidence).',
                '',
                '| Why | Records |',
                '|---|---:|',
            ] + [f"| `{k}`: {rc.CODES['lower_confidence'][k]} | {v:,} |" for k, v in counts.items() if k] + ['']
        staged = [x for x in rc.STUDIES if x[0] in set(rc.STUDY_TABLES.values()) and x[0] not in set(rc.STUDIES_ADDED)]
        site_studies = placed.set_index('site_id')['studies'].str.split('; ')
        prev_by = prev_studies
        rows_ = []
        for sid, label, citation, doi in staged:
            g = records[(records['study'] == sid) & (records['source_table'].isin([k for k, (raw, _) in rc.SOURCE_TABLES.items() if raw in rc.STUDY_TABLES]))]
            if g.empty:
                continue
            here = {x for x in g['placed_site'] if x} & set(site_studies.index)
            before = {x for x, st in prev_by.items() if sid in st}
            rows_.append(f'| [{label}](https://doi.org/{doi}) | {len(g):,} | {int((g["placed_site"] != "").sum()):,} | {len(here):,} | {len(before):,} |')
        if rows_:
            lines += [
                '## Studies read from their own tables',
                '',
                'The studies compiled in Maron et al. 2021 Table S5, and UniProt and IEDB, read from their own tables. Records: '
                'those the staging keeps (localised, and within the study\'s false discovery rate); sites: those they support '
                f'now, beside the sites the study supported in {previous} through the compiled tables.',
                '',
                f'| Study | Records | Placed | Sites | Sites in {previous} |',
                '|---|---:|---:|---:|---:|',
            ] + rows_ + ['', '`STUDIES.md` gives, for each, the tables read, the thresholds, and records traced from the source table to the '
                         'reviewed proteome.', '']
    added = [s for s in rc.STUDIES if s[0] in set(new['study']) and s[0] in set(rc.STUDIES_ADDED)]
    if added:
        lines += [
            f'## Studies added ({len(added)})',
            '',
            'Each study\'s supplementary table is staged like the others: one record per reported peptide and arginine, kept when '
            'it passes the study\'s false discovery rate and the methyl group is localised to that arginine, and placed by the '
            'rule above.',
            '',
            '| Study | Sample | Records | Placed | Sites | Sites new to the union | Proteins new to the union |',
            '|---|---|---:|---:|---:|---:|---:|',
        ]
        site_studies = placed.set_index('site_id')['studies'].str.split('; ')
        for sid, label, citation, doi in added:
            g = new[new['study'] == sid]
            here = {s for s in g['placed_site'] if s}
            here_placed = here & set(site_studies.index)
            new_sites = here_placed & added_ids
            new_prot = {split_key(x)[0] for x in new_sites} - prev_proteins
            sample = '; '.join(sorted({x for v in g['reported_sample'] for x in v.split('; ') if x}))
            lines.append(f'| [{label}](https://doi.org/{doi}) | {sample} | {len(g):,} | {int((g["placed_site"] != "").sum()):,} | '
                         f'{len(here_placed):,} | {len(new_sites):,} | {len(new_prot):,} |')
        only_new = placed[placed['site_id'].isin(added_ids)]
        lines += [
            '',
            f"{n(len(only_new))} sites are new to the union, on {n(only_new['uniprot_accession'].nunique())} proteins "
            f"({n(len({a for a in only_new['uniprot_accession']} - prev_proteins))} of them new).",
            '',
        ]
    nl = placed[placed['dimethyl_type_neutral_loss'] != '']
    if since_v11 and len(nl):
        ab = nl['maron2021_s3_antibodies'].str.split('; ').map(lambda v: {x for x in v if x in ('ADMA', 'SDMA')})
        call = nl['dimethyl_type_neutral_loss'].str.split('; ').map(set)
        both = [(a, b) for a, b in zip(call, ab) if b]
        same = sum(a == b for a, b in both)
        part = sum(bool(a & b) and a != b for a, b in both)
        lines += [
            '### ADMA and SDMA from neutral losses',
            '',
            f"Hartel et al. 2019 tell asymmetric from symmetric dimethylarginine by the neutral loss in the fragment "
            f"spectra: dimethylamine (45.058 Da) for ADMA, methylamine (31.042 Da) for SDMA. `dimethyl_type_neutral_loss` "
            f"gives the call at {n(len(nl))} sites: "
            + ', '.join(f"{v} {k}" for k, v in nl['dimethyl_type_neutral_loss'].value_counts().items()) + '. '
            f"Maron et al. 2021 Table S3 says which antibodies recovered a site (`maron2021_s3_antibodies`). At the "
            f"{n(len(both))} sites with both a neutral-loss call and an ADMA or SDMA immunoprecipitation, the two agree at "
            f"{n(same)}, overlap at {n(part)} and differ at {n(len(both) - same - part)}. Both columns are kept as reported.",
            '',
        ]
    xeno = int(placed['flags'].str.contains('xenograft_only').sum())
    if xeno:
        lines += [
            '### Xenograft-only sites',
            '',
            f"Ma et al. 2023 measured breast cancer patient-derived xenografts grown in mice. The {n(xeno)} sites only "
            'they report carry the flag `xenograft_only`: where the human and mouse sequences agree, the peptide can come '
            'from the mouse stroma.',
            '',
        ]
    lines += nucleic_acid_section(ctx, files, previous)
    if since_v11:
        lines += [
        '## Methylation degree',
        '',
        f"`methyl_state` gives a degree at {n(c['with_state'])} of {n(c['sites'])} sites: "
        + ', '.join(f"{n(v)} {k.replace('_', ' + ')}" for k, v in c['state'].items() if k != 'none') + '. '
        'The studies read from their own tables count through their records (`methyl_states_reported`).',
        '',
        '## Protein-set tests',
        '',
        '- Every test of protein or gene sets (PANTHER and functional classes, CD-CODE condensates, Human Protein Atlas '
        'locations, disease families, Gene Ontology binding, complexes) is an over-representation analysis with '
        'clusterProfiler::enricher, and the DisGeNET terms with DOSE::enrichDGN (scripts/enrich_sets.R; statistics: method '
        'enricher, with fold enrichment): p is the one-sided hypergeometric probability, fold enrichment (k/n)/(M/N), q '
        'Benjamini-Hochberg over the sets of the analysis (README, Statistics). The disease families are tested against the '
        'reviewed human genes with a DisGeNET term.',
        '- The Human Disease Ontology is added, with the gene sets of DOSE enrichDO (`supporting/disease_ontology`), and every '
        'DisGeNET term tested is listed (`supporting/disease_terms_tested`).',
        '- The analyses run on the sites of this version.',
        '',
        '## Corrections to the analyses',
        '',
        '- Functional classes: short gene tokens of `config/functional_ontology.json` matched inside other words (`atr` in '
        'matrix and tetratricopeptide, `actin` in interacting, `lamin` in filamin, `nop` in aminopeptidase); they now match '
        'whole words or word starts. Every protein is labelled from the same text (UniProt and InterPro); the names in '
        'the source tables, which list mostly methylated proteins, are no longer added to theirs. The classes called '
        'enriched do not change.',
        '- Condensates: the background holds reviewed entries only, secondary accessions counted as their current entry.',
        '- Disease enrichment: a protein counts under its UniProt primary gene (the source tables\' first symbol could be a '
        'neighbouring gene, e.g. IGF2 for IGF2-AS), and the Disease Ontology table names genes by symbol.',
        '- Density windows: a window of W residues spans W residues (it spanned W+1), each protein\'s null is seeded from '
        'its accession, and the null takes 10,000 permutations (1,000 before).',
        '- Motifs: the proteome background leaves out the methylated arginines. Per-protein arginine enrichment: q over the '
        'proteins tested.',
        '- Disorder: a protein end is not an edge of a disordered region; MobiDB entries numbered on another sequence than '
        'the reference were replaced from the MobiDB API or dropped (`inputs/annotation/mobidb_entries_patched_2026-10-21.tsv`).',
        '- Human Protein Atlas: an accession several HPA genes map to takes the best reliability and the main location of any.',
        '- Protein names keep brackets that belong to the name (DNA (cytosine-5)-methyltransferase 1); gene synonyms are '
        'split at the semicolons between genes.',
        '- Interaction evidence: predicted contacts and interface residues whose residues disagree with the reference '
        'sequence (a chain of another protein, another numbering) and ProtVar pairs with a non-human partner are left out; '
        'PDB contacts in fusion chains name the protein of the segment that holds the residue; CD-CODE members listed '
        'under isoforms count once, with their best confidence.',
        '',
        ]
    lines += [
        '## Still open',
        '',
    ] + ([f"- **Windows that match on one side only** ({n(c['pending'])} sites, flag `one_sided_window`): kept, flagged and left "
          'out of every analysis.'] if c['pending'] else []) + ([
        f"- **Paralog windows** ({c['flags'].get('paralog_window', 0)} sites, flag `paralog_window`), accepted in "
        '`config/remap_evidence_allowlist.tsv`.'] if c['flags'].get('paralog_window', 0) else []) + [
        '- **Sites without a readable localisation** (' + ('above' if since_v11 else 'Site specificity in the CHANGES.md of v1.2.0') + '): mostly sites '
        'whose study gives them below its localisation threshold, or whose search placed the methyl group differently from one '
        'confident spectrum to the next.',
        '',
        '## Tracing a change',
        '',
        f"- `{files['changes'][0]}`: one row per site of either version, with its change and the records involved.",
        f"- `{files['records'][0]}`: one row per source record, with what the study reported, where each version placed it "
        'and the evidence check.',
        '',
        'Their columns are listed in `README.md` and `data_dictionary.tsv`.',
    ]
    (outdir / 'CHANGES.md').write_text('\n'.join(lines) + '\n')


def nucleic_acid_section(ctx: dict, files: dict, previous: str) -> list:
    """CHANGES.md: the RNA and DNA in contact with methylarginines, directly or via a flanking residue, and what their proteins bind (v1.3.0 on)."""
    if 'nucleic_acid_contacts' not in files or previous not in ('v1.0.0', 'v1.1.0', 'v1.2.0'):
        return []
    at, near = files['nucleic_acid_contacts'][1], files.get('nucleic_acid_near_sites', (None, None))[1]
    prox = files.get('nucleic_acid_proximal_sites', (None, None))[1]
    regions = files.get('nucleic_acid_regions', (None, None))[1]
    kind = lambda df, k: int(df.loc[df['nucleic_acid'].isin(['rna', 'hybrid'] if k == 'rna' else ['dna']), 'site_id'].nunique())  # noqa: E731
    out = [
        '## RNA and DNA at methylarginines',
        '',
        f'The sites are those of {previous}. This version adds the RNA and DNA in contact with them in PDB structures, read in every entry '
        'PDBe-KB lists for a site with RNA or DNA, and what their proteins bind in vitro and in cells.',
        '',
        f"- `{files['nucleic_acid_contacts'][0]}`: {n(kind(at, 'rna'))} methylarginines in direct contact with RNA and {n(kind(at, 'dna'))} with DNA "
        f"(a nucleotide within 4 Å of the side chain), {n(len(at))} site-chain pairs in {n(at['pdb_id'].nunique())} entries: the "
        'chain (molecule, organism, how it was made: Synthetically obtained for a designed oligonucleotide), the nucleotides '
        'contacted and the bound sequence 8 nucleotides each side of the closest one, so the windows of any sites align.',
    ]
    if prox is not None:
        out.append(f"- `{files['nucleic_acid_proximal_sites'][0]}`: {n(kind(prox, 'rna'))} sites with RNA and {n(kind(prox, 'dna'))} with DNA "
                   'within 8 Å of the arginine side chain in a complex where it has no direct contact (one entry read per distinct '
                   'complex): the distance and the bound sequence around the closest nucleotide.')
    if near is not None:
        out.append(f"- `{files['nucleic_acid_near_sites'][0]}`: {n(kind(near, 'rna'))} sites with RNA and {n(kind(near, 'dna'))} with DNA "
                   'in contact via a flanking residue (one at i ± 1-5 within 4 Å, in a complex where the arginine side chain is not within 8 Å): '
                   'the bound sequence there, centred on the nucleotide the closest such residue contacts.')
    if regions is not None:
        has = regions['methylarginines'].astype(str).str.len() > 0
        out.append(f"- `{files['nucleic_acid_regions'][0]}`: {n(len(regions))} runs of residues within 4 Å of RNA or DNA on "
                   f"{n(regions['uniprot_accession'].nunique())} proteins ({n(int(has.sum()))} with methylarginines inside), read once per "
                   'distinct complex (the same bound sequences), with the segments each contacts.')
    for key, text in (('rna_bind_n_seq', 'ENCODE RNA Bind-n-Seq (in vitro): the most enriched 5-mers and logos'),
                      ('dna_motifs', 'JASPAR 2026 CORE vertebrates DNA-binding profiles (CC BY 4.0)'),
                      ('eclip', 'ENCODE eCLIP in K562 and HepG2 cells: peaks by transcript region (GENCODE v50 basic) and the genes with most')):
        if key in files:
            df = files[key][1]
            out.append(f"- `{files[key][0]}`: {text}, for {n(df['uniprot_accession'].nunique())} methylated proteins.")
    out += ['', 'The website shows them on protein pages (RNA and DNA interactions: where the protein binds, the bound sequences '
            'in contact with each methylarginine with a logo and a type filter, its preferences and eCLIP targets), on a site\'s card, in '
            'the site tables (Bound RNA/DNA, with the search term bound:), and for a selection in the Interactions view.'
            ' Each chain is typed (rRNA, tRNA, snRNA, mRNA, other ncRNA, viral, designed oligo; nucleosomal, telomeric or other '
            'DNA) and the site table gains bound_* columns. Inputs: `inputs/interactions/pdb_nucleic_acid_*` (PDB data, CC0) and '
            '`inputs/nucleic_acid/` (scripts/43_stage_nucleic_acid_binding.py). RNAcompete (CISBP-RNA) states no licence and '
            'is not included.', '']
    return out


def anchor(label: str) -> str:
    """The anchor GitHub gives a Markdown heading."""
    return re.sub(r'[^a-z0-9 -]', '', label.lower()).replace(' ', '-')


def write_studies(outdir: Path, ctx: dict) -> None:
    """STUDIES.md: for every study and database behind the sites, where its records come from, what they hold and how
    they were filtered (release_codes STUDY_NOTE, STUDY_QUALITY), what became of them (placed, and why the others were
    not), and records traced from the source table to the reviewed proteome (build_release.study_examples)."""
    version, studies, records = ctx['version'], ctx['studies'], ctx['records']
    examples = ctx.get('study_examples', {})
    st = studies.set_index('study_id')
    order = [x for x in rc.STUDIES if x[0] in st.index]
    words = {code: text for code, (_, text) in rc.SOURCE_TABLES.items()}
    lines = [
        f'# Studies and databases behind the sites, {version}',
        '',
        'Every site rests on records: rows of the source tables, each crediting the study or database that reported the '
        'site. A record is placed on a reviewed UniProt entry by its own sequence evidence (a window centred on the site, '
        'or a peptide with or without the site marked), and only by its reported position when it gives no sequence; every '
        'placed record is then checked against that evidence, and every site is an arginine of the reference sequence. The '
        'records of a study read from its own table add a site only when the study localised the methyl group to one '
        'arginine and the peptide (or, without one, the window) names exactly one arginine in the whole reviewed proteome, '
        'in a protein the record names; the rows of the compiled tables (Maron et al. 2021 Table S5, the PTM master table) '
        'carry no localisation and add none.',
        '',
        'Each section gives where the study\'s records come from, what they hold, what became of them, and records traced '
        'from the source table to the reviewed proteome: a placed record, with its evidence aligned to the reference '
        'sequence (brackets mark the site), and a record left out, with the reason. The placement codes are explained in '
        '`data_dictionary.tsv`; every record is in `supporting/source_records`.',
        '',
        '| Study or database | Records from | Records | Placed | Sites | Reported by no other source |',
        '|---|---|---:|---:|---:|---:|',
    ]
    for sid, label, citation, doi in order:
        r = st.loc[sid]
        lines.append(f"| [{label}](#{anchor(label)}) | {'; '.join(words.get(t, t) for t in str(r['source_table']).split('; '))} | "
                     f"{n(int(r['n_records']))} | {n(int(r['n_records_placed']))} | {n(int(r['n_sites']))} | {n(int(r['n_sites_only_here']))} |")
    lines.append('')
    for sid, label, citation, doi in order:
        r = st.loc[sid]
        g = records[records['study'] == sid]
        left = g[(g['placed_site'] == '') & (g['pending_site'] == '')]['placement_code'].value_counts()
        note = r['note'][:1].upper() + r['note'][1:] if r['note'] else ''
        lines += [f'## {label}', '', f'{citation}, [doi:{doi}](https://doi.org/{doi}).' + (f' {note}.' if note else ''), '',
                  '**Records from** ' + '; '.join(words.get(t, t) for t in str(r['source_table']).split('; '))
                  + f". **Terms:** {r['licence']}.", '']
        if r['quality']:
            lines += [f"> {r['quality']}", '']
        pending = int(((g['placed_site'] == '') & (g['pending_site'] != '')).sum())
        lines.append(f"{n(int(r['n_records']))} record{'s' if int(r['n_records']) != 1 else ''}: {n(int(r['n_records_placed']))} placed on "
                     f"{n(int(r['n_sites']))} site{'s' if int(r['n_sites']) != 1 else ''} ({n(int(r['n_sites_in_previous_version']))} of "
                     f"them in the previous version, {n(int(r['n_sites_only_here']))} reported by no other source)"
                     + (f', {n(pending)} pending' if pending else '') + '.'
                     + (' Not placed: ' + '; '.join(f'{n(int(k))} `{code}` ({re.sub(r"^not placed( [(][^)]*[)])?: ", "", rc.CODES["placement"].get(code, ""))})' for code, k in left.items()) + '.'
                        if len(left) else ' Every record is placed.'))
        lines.append('')
        for block in examples.get(sid, []):
            lines += ['```text', *block, '```', '']
    (outdir / 'STUDIES.md').write_text('\n'.join(lines))


def write_md5sums(outdir: Path) -> None:
    with open(outdir / 'MD5SUMS', 'w') as fh:
        for p in sorted(outdir.rglob('*')):
            if p.is_file() and p.name != 'MD5SUMS':
                fh.write(f'{md5(p)}  {p.relative_to(outdir).as_posix()}\n')

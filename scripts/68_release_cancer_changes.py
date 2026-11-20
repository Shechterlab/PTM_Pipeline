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
SOMATIC = {'tcga_mc3': 'TCGA MC3', 'depmap_24q4': 'DepMap 24Q4'}
GERMLINE = {'clinvar': 'ClinVar', 'gnomad_v4_1': 'gnomAD v4.1'}
SOURCES = SOMATIC | GERMLINE
WHO = {'tcga_mc3': 'sample_id', 'depmap_24q4': 'model_id', 'clinvar': 'variation_id', 'gnomad_v4_1': 'variant_id'}
PATHOGENIC = r'^(Pathogenic|Likely pathogenic|Pathogenic/Likely pathogenic)(;|$|,)'
# the line of the carried interaction notes (scripts/64) that said the RBDmap and pCLAP extents are not republished: the
# methylarginine tables carry them, with their sequences, credited (owner's decision, 19 November 2026)
WRONG_LINE = '- RBDmap fragments and pCLAP peptides: their extents are not republished.'
RIGHT_LINE = ('- RBDmap fragments and pCLAP peptides: `scripts/64` does not read them for these sets. The methylarginine tables carry '
              'them, with their extents and sequences: RBDmap from Castello et al., Mol Cell 63:696 (2016), doi:10.1016/j.molcel.2016.06.029 '
              '(CC BY 4.0); pCLAP, the RNA-binding peptides of Mullari et al., J Proteome Res 16:2762 (2017), doi:10.1021/acs.jproteome.7b00042, '
              'called by the pipeline from PRIDE PXD005739 with the rule of Mullari et al., Nat Commun 14:4348 (2023).')
NAMES = {'R': 'arginine', 'K': 'lysine', 'S': 'serine', 'T': 'threonine', 'Y': 'tyrosine'}
COPIED = {'cancer_change_protein_background.tsv.gz': 'cancer_change_protein_background_{v}.tsv.gz',
          'tcga_mc3_projects.tsv': 'tcga_mc3_projects_{v}.tsv', 'depmap_24q4_change_models.tsv': 'depmap_24q4_change_models_{v}.tsv',
          'clinvar_protein_background.tsv.gz': 'clinvar_protein_background_{v}.tsv.gz',
          'gnomad_protein_background.tsv.gz': 'gnomad_protein_background_{v}.tsv.gz'}


def staged(src: str, res: str) -> Path:
    # every family and residue from its staged table; the TCGA and DepMap arginine rows from scripts/66 (the rows of the arginine
    # layer, byte for byte, with the gene), the arginine layer tables themselves left as they are
    return LAYERS / f'{src}_changes_{res}.tsv.gz'


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
                ('gene', 'string', '', '', 'HGNC approved symbol of the gene of the transcript the change was read on (its Ensembl 105 gene; else the symbol the MC3 MAF gives, resolved through HGNC previous symbols; else the Ensembl gene ID)'),
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
                ('gene', 'string', '', '', 'HGNC approved symbol of the gene (DepMap\'s Ensembl gene, else its HugoSymbol)'),
                ('site_sets', 'string', '', '', site_sets)):
            rows.append('\t'.join([f, col, typ, vals, means, desc]))
        f = released('clinvar', res, version)
        for col, typ, vals, means, desc in (
                ('canonical_UniProtAC', 'string', '', '', ref),
                ('position', 'integer', '', '', f'position of the {NAMES[res]} on the canonical sequence'),
                ('alt', *alt),
                ('variation_id', 'string', '', '', 'ClinVar variation ID (https://www.ncbi.nlm.nih.gov/clinvar/variation/<ID>/ lists its records)'),
                ('refseq_protein', 'string', '', '', 'RefSeq protein the change is written on, a displayed isoform of the entry'),
                ('variant_type', 'string', '', '', 'ClinVar variant type (single nucleotide variant, ...)'),
                ('germline_classification', 'string', '', '', 'ClinVar aggregate germline classification (Pathogenic, Likely pathogenic, Uncertain significance, ...)'),
                ('clinsig_simple', 'integer', '1; 0; -1', '1 = at least one current record pathogenic or likely pathogenic; 0 = none; -1 = no classification', 'ClinVar ClinSigSimple'),
                ('review_status', 'string', '', '', 'ClinVar review status of the classification'),
                ('review_stars', 'integer', '0; 1; 2; 3; 4', '4 = practice guideline; 3 = reviewed by expert panel; 2 = criteria provided, multiple submitters, no conflicts; 1 = criteria provided, single submitter or conflicting classifications; 0 = otherwise', 'ClinVar review stars'),
                ('number_submitters', 'integer', '', '', 'submitters of the variation'),
                ('origin', 'string', '', '', 'ClinVar OriginSimple (germline, somatic, ...)'),
                ('conditions', 'string', '', '', 'ClinVar PhenotypeList, separated by | as ClinVar gives them'),
                ('last_evaluated', 'string', '', '', 'date the classification was last evaluated'),
                ('somatic_clinical_impact', 'string', '', '', 'ClinVar somatic clinical impact classification'),
                ('oncogenicity', 'string', '', '', 'ClinVar oncogenicity classification'),
                ('dbsnp', 'string', '', '', 'dbSNP ID'),
                ('gene', 'string', '', '', 'HGNC approved symbol of the gene ClinVar gives the RefSeq protein\'s expression (its NCBI Gene ID; else its symbol, else variant_summary\'s)'),
                ('site_sets', 'string', '', '', site_sets)):
            rows.append('\t'.join([f, col, typ, vals, means, desc]))
        f = released('gnomad_v4_1', res, version)
        for col, typ, vals, means, desc in (
                ('canonical_UniProtAC', 'string', '', '', ref),
                ('position', 'integer', '', '', f'position of the {NAMES[res]} on the canonical sequence'),
                ('alt', *alt),
                ('variant_id', 'string', '', '', 'gnomAD variant, GRCh38 chromosome-position-reference-alternate (gnomad.broadinstitute.org/variant/<ID>?dataset=gnomad_r4)'),
                ('dbsnp', 'string', '', '', 'dbSNP ID'),
                ('ac', 'integer', '', '', 'allele count in gnomAD v4.1 exomes'),
                ('an', 'integer', '', '', 'allele number (alleles genotyped)'),
                ('af', 'number', '', '', 'allele frequency'),
                ('nhomalt', 'integer', '', '', 'homozygous individuals'),
                ('grpmax', 'string', '', '', 'genetic ancestry group with the highest allele frequency'),
                ('af_grpmax', 'number', '', '', 'allele frequency in that group'),
                ('faf95_max', 'number', '', '', 'highest filtering allele frequency (Poisson 95% CI) across groups'),
                ('consequence', 'string', 'missense_variant; stop_gained', '', 'VEP consequence of the change'),
                ('gene', 'string', '', '', 'HGNC approved symbol of the Ensembl 105 gene of the transcript VEP read the change on; several separated by semicolons when genes overlap; the Ensembl gene ID for a gene without an HGNC symbol'),
                ('site_sets', 'string', '', '', site_sets)):
            rows.append('\t'.join([f, col, typ, vals, means, desc]))
    for name, label in (('clinvar', 'ClinVar'), ('gnomad', 'gnomAD v4.1')):
        f = f'supporting/{name}_protein_background_{version}.tsv.gz'
        cols = [('canonical_UniProtAC', 'string', '', ref), ('residue', 'string', 'R; K; S; T; Y; all', 'residue class; all = every residue of the protein'),
                ('residues', 'integer', '', 'residues of that class in the canonical sequence (its length for all)')]
        cols += ([('clinvar_variants', 'integer', '', 'ClinVar variants on them'),
                  ('clinvar_pathogenic', 'integer', '', 'of those, classified pathogenic or likely pathogenic (germline)'),
                  ('clinvar_positions', 'integer', '', 'positions changed')] if name == 'clinvar' else
                 [('gnomad_variants', 'integer', '', 'gnomAD v4.1 exome variants on them'), ('gnomad_positions', 'integer', '', 'positions changed'),
                  ('gnomad_allele_count', 'integer', '', 'summed allele count of those variants')])
        for col, typ, vals, desc in cols:
            rows.append('\t'.join([f, col, typ, vals, '', desc]))
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
    paths = ['scripts/68_release_cancer_changes.py', 'scripts/65_carry_release.py', 'scripts/66_stage_cancer_changes.py',
             'scripts/67_stage_clinvar.py', 'scripts/69_stage_gnomad.py', 'config/site_sets.json']
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
            if src == 'clinvar':
                patho = known[known['germline_classification'].str.match(PATHOGENIC)]
                stats[f'{src}/{res}'] |= {'pathogenic_rows': int(t['germline_classification'].str.match(PATHOGENIC).sum()),
                                          'pathogenic_rows_at_known_sites': len(patho),
                                          'known_sites_with_pathogenic': int(patho[['canonical_UniProtAC', 'position']].drop_duplicates().shape[0])}
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
            if src == 'clinvar':
                patho = set()
                for res in letters:
                    t = x65.read(out / released(src, res, version))
                    t = t[t['site_sets'].str.split(';').map(lambda s: name in s) & t['germline_classification'].str.match(PATHOGENIC)]
                    patho |= set(zip(t['canonical_UniProtAC'], t['position']))
                row['clinvar_pathogenic_sites'] = len(patho)
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
    clinvar = json.loads((LAYERS / 'clinvar_summary.json').read_text(encoding='utf-8'))
    gnomad = json.loads((LAYERS / 'gnomad_summary.json').read_text(encoding='utf-8'))
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

    def germline_residue_table() -> list[str]:
        out_lines = ['| Residue | ClinVar variants | at known sites | of which pathogenic or likely pathogenic | gnomAD variants | at known sites |',
                     '|---|---:|---:|---:|---:|---:|']
        for res in RESIDUES:
            cv, gn = stats[f'clinvar/{res}'], stats[f'gnomad_v4_1/{res}']
            out_lines.append(f'| {NAMES[res].capitalize()} ({res}) | {cv["rows"]:,} | {cv["rows_at_known_sites"]:,} | '
                             f'{cv["pathogenic_rows_at_known_sites"]:,} | {gn["rows"]:,} | {gn["rows_at_known_sites"]:,} |')
        return out_lines

    def germline_set_table() -> list[str]:
        out_lines = ['| Site set | Sites | Sites with a ClinVar variant | with a pathogenic or likely pathogenic one | Sites with a gnomAD variant |',
                     '|---|---:|---:|---:|---:|']
        for r in per_set:
            out_lines.append(f'| {r["label"]} | {r["sites"]:,} | {r["clinvar_sites_changed"]:,} | {r["clinvar_pathogenic_sites"]:,} | '
                             f'{r["gnomad_v4_1_sites_changed"]:,} |')
        return out_lines

    germline_section = [
        '## New: germline variants at every candidate residue (ClinVar, gnomAD)', '',
        f'`supporting/clinvar_changes_<residue>_{version}.tsv.gz` and `supporting/gnomad_v4_1_changes_<residue>_{version}.tsv.gz` list '
        'the human germline variants that change an arginine, lysine, serine, threonine or tyrosine of a reviewed canonical '
        'protein (a substitution or a stop), a known site or not, with `site_sets` as in the cancer tables. ClinVar (files of '
        f'{clinvar["inputs"]["clinvar_files_dated"]}): one row per variation, with its germline classification, review status and '
        'stars, submitters, origin, conditions, somatic clinical impact and oncogenicity; the variant counts where one of its '
        'protein changes is written on a RefSeq protein UniProt assigns to the entry\'s displayed isoform. gnomAD v4.1 (730,947 '
        'exomes): one row per variant that passes every gnomAD filter, with its allele count and number, frequency, homozygotes, '
        'the ancestry group with the highest frequency and the highest filtering allele frequency; the change is VEP\'s, on a '
        'transcript UniProt assigns to the displayed isoform. Either way a variant counts only where the 2026_03 canonical sequence '
        'has that residue. As in the cancer tables, `gene` names the gene of each variant.', '',
        *germline_residue_table(), '',
        *germline_set_table(), '',
        f'`supporting/clinvar_protein_background_{version}.tsv.gz` and `supporting/gnomad_protein_background_{version}.tsv.gz` give, for '
        'each protein and residue class, the residues and the variants on them (and, for ClinVar, those pathogenic or likely '
        'pathogenic; for gnomAD, the summed allele count).', '',
        'NCBI places no restrictions on the use or distribution of ClinVar data (submitters may claim rights in what they '
        'submitted); cite ClinVar (Landrum et al., Nucleic Acids Res 46:D1062, 2018). The gnomAD primary data are CC0; cite gnomAD v4 '
        '(Chen et al., Nature 625:92, 2024). gnomAD\'s SpliceAI and CADD annotations, which carry non-commercial terms, are not used. '
        'gnomAD was read from Google\'s public copy (`gs://gcp-public-data--gnomad/release/4.1/vcf/exomes/`), streamed '
        '(`scripts/69_stage_gnomad.py`); ClinVar from its FTP tables (`scripts/67_stage_clinvar.py`).', '']

    new_section = [
        '## New: somatic changes at every candidate residue', '',
        f'`supporting/tcga_mc3_changes_<residue>_{version}.tsv.gz` and `supporting/depmap_24q4_changes_<residue>_{version}.tsv.gz` '
        'list the somatic missense and stop changes at every arginine, lysine, serine, threonine and tyrosine of the reviewed '
        'canonical proteins, a known site or not: in TCGA tumors (the MC3 calls of the PanCancer Atlas, '
        f'{den["tcga_tumor_samples"]:,} tumor samples in {den["tcga_projects"]} projects) and in DepMap 24Q4 cell lines '
        f'({den["depmap_lines_with_mutation_calls"]:,} lines with mutation calls). One row per tumor sample or cell line and change. '
        'The column `site_sets` names the site sets of this release with a site at that residue (empty: a candidate residue, not '
        'a known site), so a set\'s changes can be read at every candidate residue or at its known sites only. The column `gene` '
        'names the gene each change lies in (its current HGNC approved symbol, for all four sources), which matters where several genes encode one protein: the 14 genes '
        'of histone H4 (P62805), for example, are pooled under one accession.', '',
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

    carried_text = x65.carried_readme(prev, previous, version, arginine_sets).rstrip('\n')
    fixed_line = WRONG_LINE in carried_text
    carried_text = carried_text.replace(WRONG_LINE, RIGHT_LINE)
    R = [f'# Human methylarginine sites, context marks, site-level interaction evidence, and cancer and germline variants, {version}', '',
         f'{version} holds everything of {previous} under this release\'s names, its data unchanged: the methylarginine data of '
         'v1.4.3, the tables of context marks and the site-level interaction tables of every site set. It adds, at every candidate '
         'residue of the site sets, the somatic changes in TCGA tumors and DepMap cell lines and the germline variants of ClinVar '
         'and gnomAD, each with the site sets the residue belongs to, and per-protein background counts. See CHANGES.md for what '
         'changed.', '',
         *new_section,
         *germline_section,
         carried_text]
    (out / 'README.md').write_text(x65.americanize('\n'.join(R)) + '\n', encoding='utf-8')

    new_files = sorted(files)
    C = [f'# Changes from {previous} to {version}', '',
         f'Everything of {previous} is in {version} under this release\'s names, its data unchanged. {version} adds the somatic changes '
         'at every candidate residue of the site sets, in TCGA tumors and DepMap cell lines, and the germline variants of ClinVar '
         'and gnomAD there.', '',
         *new_section,
         *germline_section,
         '## Corrected: RBDmap and pCLAP in the interaction notes', '',
         f'The interaction notes of {previous}\'s README said, under what the tables of methyllysine, lysine and phospho lack, that the '
         'extents of RBDmap fragments and pCLAP peptides are not republished. They are: the methylarginine interaction tables carry '
         'each RBDmap fragment and pCLAP peptide with its extent and sequence, credited (RBDmap: Castello et al. 2016, CC BY 4.0; pCLAP: '
         'Mullari et al. 2017, called from PRIDE PXD005739 with the rule of Mullari et al. 2023). This README says so'
         + ('' if fixed_line else ' (the line was not found to replace)') + f'; {previous}\'s published folder is left as archived.', '',
         '## New files', '', *[f'- `{f}`' for f in new_files], '',
         '## Files renamed', '',
         f'Every carried file is renamed `_{previous}` to `_{version}`, its content unchanged except where it names its version: '
         + ', '.join(f'`{e}`' for e in edited) + '. '
         f'`{changes_name}` replaces {previous}\'s site changes table; it lists every one of the {len(sites):,} methylarginine sites as kept.', '',
         'Code: `scripts/66_stage_cancer_changes.py` (the K, S, T and Y change tables and the background), '
         '`scripts/67_stage_clinvar.py` (ClinVar), `scripts/69_stage_gnomad.py` (gnomAD), '
         '`scripts/68_release_cancer_changes.py` (this release), `scripts/65_carry_release.py` (carrying the previous release), at '
         f'commit `{commit}`.']
    (out / 'CHANGES.md').write_text(x65.americanize('\n'.join(C)) + '\n', encoding='utf-8')

    som = {k: v for k, v in files.items() if not any(s in k for s in ('clinvar', 'gnomad'))}
    germ = {k: v for k, v in files.items() if any(s in k for s in ('clinvar', 'gnomad'))}
    release = {'title': 'Human methylarginine sites, context marks, site-level interaction evidence, and cancer and germline variants', 'version': version,
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
               'cancer_changes': {'files': som, 'per_residue': {k: v for k, v in stats.items() if k.split('/')[0] in SOMATIC},
                                  'per_site_set': per_set, 'denominators': den, 'staging_inputs': staging['inputs'], 'rules': staging['rule']},
               'germline_variants': {'files': germ, 'per_residue': {k: v for k, v in stats.items() if k.split('/')[0] in GERMLINE},
                                     'per_site_set': per_set,
                                     'clinvar': {'inputs': clinvar['inputs'], 'rule': clinvar['rule'], 'license': clinvar['license']},
                                     'gnomad': {'release': gnomad['release'], 'inputs': gnomad['inputs'], 'rule': gnomad['rule'],
                                                'license': gnomad['license'], 'citation': gnomad['citation']}},
               'corrections': ['README: RBDmap and pCLAP extents are carried in the methylarginine tables (the previous README said they '
                               'were not republished)'] if fixed_line else []}
    (out / 'release.json').write_text(json.dumps(release, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    names = sorted(p.relative_to(out).as_posix() for p in out.rglob('*') if p.is_file() and p.name != 'MD5SUMS')
    (out / 'MD5SUMS').write_text(''.join(f'{md5_file(out / n)}  {n}\n' for n in names), encoding='utf-8')
    print(json.dumps({'per_residue': stats, 'per_site_set': per_set}, indent=1))


if __name__ == '__main__':
    main()

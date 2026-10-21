#!/usr/bin/env python3
"""Run the release analyses on any release's sites, and collect them into the release folder.

    python scripts/release_analyses.py prepare --sites releases/v1.1.0/human_methylarginine_sites_v1.1.0.tsv \
        --entries inputs/reference/uniprot_human_reviewed_2026_03_entries.tsv --work WORK
    (hpc/analyses.sh then runs the analyses on WORK/integrated/; hpc/release_analyses.sh does both)
    python scripts/release_analyses.py collect --work WORK --inputs inputs --version v1.1.0 --sites ... --out releases/v1.1.0/analyses

prepare writes the site union the analysis scripts read (WORK/integrated/human_arg_methyl_union_dedup_by_site.tsv) from a
release's site table: one row per site, on the reviewed reference sequences, without the sites kept pending (flag
one_sided_window). v1.0.0 numbered its sites on an earlier UniProt release; --map (site_changes_v1.0.0_to_v1.1.0.tsv)
takes each v1.0.0 site to where v1.1.0 placed its records on the reference, and leaves out the sites v1.1.0 could not
place. collect writes the analyses' result tables, the statistics table of every test (release_tables.statistics), a
README and MD5SUMS into --out, so every release can be compared with the same code.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_release as br  # noqa: E402
import release_codes as rc  # noqa: E402
import release_tables as rt  # noqa: E402

CODE = Path(__file__).resolve().parents[1]
# The result tables kept, by folder: (analysis output folder, file, what it holds).
KEEP = [
    ('functional_classes', 'functional_class_union', 'arg_methyl_functional_class_broad_enrichment.tsv', 'functional classes, methylated proteins vs the reviewed proteome'),
    ('functional_classes', 'functional_class_union', 'arg_methyl_functional_class_subclass_enrichment.tsv', 'functional subclasses'),
    ('panther', 'panther_class', 'panther_class_enrichment.tsv', 'PANTHER protein classes'),
    ('motifs', 'motif_arg_odds', 'motif_family_enrichment_vs_all_proteome_arginines.tsv', 'motif families at the site, vs every arginine of the proteome'),
    ('motifs', 'motif_arg_odds', 'motif_family_enrichment_vs_methylated_protein_arginines.tsv', 'motif families, vs the other arginines of the methylated proteins'),
    ('motifs', 'motif_arg_odds', 'positional_amino_acid_enrichment.tsv', 'residues at -7..+7 of the site'),
    ('disorder', 'disorder_context', 'disorder_context_summary.tsv', 'sites in MobiDB disordered regions'),
    ('disorder', 'disorder_context', 'idr_enrichment_vs_all_arginines.tsv', 'disorder, vs every arginine'),
    ('disorder', 'disorder_context', 'idr_enrichment_vs_arginines_in_methylated_proteins.tsv', 'disorder, vs the other arginines of the methylated proteins'),
    ('disorder', 'disorder_context', 'idr_edge_bin_enrichment_vs_all_disordered_arginines.tsv', 'distance to the edge of the disordered region'),
    ('domains', 'domain_context', 'domain_context_summary.tsv', 'sites inside, near and away from InterPro domains'),
    ('domains', 'domain_context', 'domain_context_enrichment.tsv', 'domain context, vs the other arginines of the methylated proteins'),
    ('domains', 'domain_context', 'domain_class_enrichment.tsv', 'nearest domain class'),
    ('domains', 'domain_context', 'domain_edge_distance_summary.tsv', 'distance to the nearest domain edge'),
    ('domains', 'domain_context', 'domain_cluster_architecture_summary_ge3.tsv', 'domain-cluster architectures (3 or more sites)'),
    ('other_ptms', 'ptm_compare', 'ptm_context_comparison.tsv', 'disorder and domain context of methylarginines vs other PTM sites'),
    ('other_ptms', 'ptm_compare', 'ptm_idr_binary_comparison.tsv', 'disordered fraction, methylarginines vs other PTM sites'),
    ('clustering', 'methyl_arg_clustering', 'nearest_neighbor_summary_5_10_20.tsv', 'another methylarginine within 5, 10, 20 residues, vs a same-protein arginine null'),
    ('clustering', 'methyl_arg_clustering', 'nearest_neighbor_curve_vs_null.tsv', 'nearest-neighbour distance curve vs the null'),
    ('clustering', 'methyl_arg_clustering', 'arg_methyl_multiscale_density_window_summary.tsv', 'dense windows by window size'),
    ('clustering', 'methyl_arg_clustering', 'arg_methyl_multiscale_density_significant_regions.tsv', 'significant dense regions'),
    ('condensates', 'condensates', 'condensate_enrichment.tsv', 'CD-CODE condensates'),
    ('condensates', 'condensates', 'condensate_summary.tsv', 'condensate universe and methylated fraction'),
    ('disease', 'disease', 'disease_families.tsv', 'disease families (DisGeNET terms)'),
    ('disease', 'disease', 'disease_terms.tsv', 'DisGeNET terms'),
    ('disease', 'disease', 'disease_ontology.tsv', 'Disease Ontology terms'),
    ('disease', 'disease', 'disease_input_summary.tsv', 'genes and terms tested'),
]


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def release_sites(sites_path: Path, map_path: str) -> pd.DataFrame:
    """(accession, position, gene, studies) of the release's sites on the reference, pending sites left out."""
    sites = pd.read_csv(sites_path, sep='\t', dtype=str, keep_default_na=False)
    if 'flags' in sites.columns:
        sites = sites[~sites['flags'].str.contains('one_sided_window')]
    studies = sites['studies'] if 'studies' in sites.columns else sites['source_families'].str.replace(';', '; ')
    out = pd.DataFrame({'site_id': sites['site_id'], 'gene': sites['gene'], 'studies': studies})
    if map_path:
        ch = pd.read_csv(map_path, sep='\t', dtype=str, keep_default_na=False)
        to = dict(zip(ch['previous_site_id'], ch['site_id'].where(ch['in_this_version'] == 'true', '')))
        out['site_id'] = out['site_id'].map(lambda s: to.get(s, ''))
        out = out[out['site_id'] != ''].drop_duplicates('site_id')
    acc_pos = out['site_id'].str.split(':R', expand=True)
    return out.assign(accession=acc_pos[0], position=acc_pos[1])


def prepare(args) -> None:
    sites = release_sites(Path(args.sites), args.map)
    entries = pd.read_csv(args.entries, sep='\t', dtype=str, keep_default_na=False).set_index('Entry')
    gene = sites['accession'].map(lambda a: entries['Gene Names (primary)'].get(a, '')).where(lambda g: g != '', sites['gene'])
    fam = sites['studies'].str.replace('; ', ';')
    n = fam.str.count(';') + 1
    union = pd.DataFrame({
        'substrate_UniProtAC': sites['accession'], 'canonical_UniProtAC': sites['accession'], 'substrate_genename': gene,
        'protein_name': sites['accession'].map(lambda a: entries['Protein names'].get(a, '')),
        'UniProtID': sites['accession'].map(lambda a: entries['Entry Name'].get(a, '')),
        'position': sites['position'], 'corrected_position': sites['position'], 'site_key': sites['site_id'],
        'site': 'R' + sites['position'], 'residue': 'R', 'source_families_exact': fam, 'source_families_fuzzy': fam,
        'source_family_count_exact': n.astype(str), 'source_family_count_fuzzy': n.astype(str)})
    out = Path(args.work) / 'integrated' / 'human_arg_methyl_union_dedup_by_site.tsv'
    out.parent.mkdir(parents=True, exist_ok=True)
    union.to_csv(out, sep='\t', index=False)
    print(f'{len(union):,} sites on {union["canonical_UniProtAC"].nunique():,} proteins -> {out}')


def collect(args) -> None:
    work, out = Path(args.work), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    sites = release_sites(Path(args.sites), args.map)
    proteins = pd.DataFrame({'uniprot_accession': sorted(set(sites['accession']))})
    proteins.attrs['hpa_all'] = rt.read(Path(args.inputs) / 'layers' / 'hpa_subcellular_locations.tsv')
    stats = rt.statistics(work, proteins)
    for name, *_ in rc.STATISTICS_COLUMNS:
        if name not in stats.columns:
            stats[name] = ''
    br.write_table(stats, rc.STATISTICS_COLUMNS, out / f'statistics_{args.version}.tsv')
    kept = []
    for folder, src, name, what in KEEP:
        (out / folder).mkdir(exist_ok=True)
        shutil.copyfile(work / src / name, out / folder / name)
        kept.append((f'{folder}/{name}', what))
    commit = subprocess.run(['git', '-C', str(CODE), 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(['git', '-C', str(CODE), 'status', '--porcelain', '--untracked-files=no'], capture_output=True, text=True).stdout.strip()
    n_sites, n_prot = len(sites), sites['accession'].nunique()
    lines = [f'# Analyses of {args.version}, rerun with the current pipeline', '',
             f'The analyses of the release, run on its {n_sites:,} sites on {n_prot:,} proteins with the code of commit '
             f'`{commit}{"+uncommitted" if dirty else ""}` ({dt.date.today().isoformat()}), so that every release is analysed the same '
             'way (hpc/release_analyses.sh). ' + (
                 'v1.0.0 numbered its sites on an earlier UniProt release: each site is taken to where v1.1.0 placed its records on '
                 'UniProt 2026_03, and the sites v1.1.0 could not place are left out. The analyses of the published run are in the '
                 'v1.0.0_archive branch. ' if args.map else '') +
             'The sites kept pending (flag one_sided_window) are left out, as in the release build.', '',
             f'`statistics_{args.version}.tsv` holds every test in one table (columns as in the release\'s statistics table: analysis, '
             'category, the 2x2 counts a-d, effect, p, BH q, background). The folders hold the result tables:', '',
             '| File | Holds |', '|---|---|'] + [f'| `{f}` | {w} |' for f, w in kept] + ['']
    (out / 'README.md').write_text('\n'.join(lines))
    files = sorted(p for p in out.rglob('*') if p.is_file() and p.name != 'MD5SUMS')
    (out / 'MD5SUMS').write_text(''.join(f'{md5(p)}  {p.relative_to(out)}\n' for p in files))
    print(f'{len(files)} files -> {out}')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    for name in ('prepare', 'collect'):
        p = sub.add_parser(name)
        p.add_argument('--sites', required=True, help="the release's site table")
        p.add_argument('--map', default='', help='v1.0.0 only: site_changes_v1.0.0_to_v1.1.0.tsv')
        p.add_argument('--work', required=True)
        if name == 'prepare':
            p.add_argument('--entries', required=True)
        else:
            p.add_argument('--inputs', required=True)
            p.add_argument('--version', required=True)
            p.add_argument('--out', required=True)
    args = ap.parse_args()
    {'prepare': prepare, 'collect': collect}[args.cmd](args)


if __name__ == '__main__':
    main()

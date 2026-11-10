#!/usr/bin/env python3
"""Stage UniProt's methyllysine, phosphorylation, acylation, methylarginine, ubiquitination and SUMOylation annotations as one table per mark (inputs/marks/).

    python scripts/54_stage_uniprot_ptm_sites.py --fasta inputs/reference/uniprot_human_reviewed_2026_03.fasta \
        --raw-dir PTM_data/marks --out-dir inputs/marks

Reads the "Modified residue" and "Cross-link" features of every reviewed human UniProtKB/Swiss-Prot entry (UniProt REST
API, CC BY 4.0; two JSON streams of `reviewed:true AND organism_id:9606`: fields accession, ft_mod_res and sequence, and
fields accession, ft_crosslnk and sequence; both must be of the same release). With --raw and --raw-crosslnk it reads
earlier downloads instead of fetching. The API serves the current release only, which may be newer than the
reference proteome: the release served is recorded (`source_release`) and an entry whose sequence differs from the
reference is left out whole, so a position is never moved between sequences.

Marks, terms and states (a term is the feature description before its first ';'):

  methyllysine     N6-methyllysine me1; N6,N6-dimethyllysine me2; N6,N6,N6-trimethyllysine me3;
                   N6-methylated lysine me (degree not stated)                                  residue K
  phosphorylation  Phosphoserine S, Phosphothreonine T, Phosphotyrosine Y: state phospho
  lysine_acylation every N6-acyl lysine term the download holds, read from its feature descriptions: 'N6-<acyl>lysine' with
                   the acyl group in plain or bracketed form ending in -yl or -oyl (N6-acetyllysine, N6-succinyllysine,
                   N6-(2-hydroxyisobutyryl)lysine, ...); state = the acyl type (acetyl, succinyl, 2-hydroxyisobutyryl, ...),
                   residue K. Not acylation: N6-methyl terms, N6-(ADP-ribosyl), N6-lipoyl, N6-biotinyl (cofactors and
                   ADP-ribose), N6-(pyridoxal phosphate), N6-(retinylidene), N6-carboxy and microbial-infection terms
  ubiquitination   "Cross-link" features "Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in ubiquitin)": state ub, residue K.
                   Not taken: ubiquitin on other residues (Cys-Gly thioester, Ser-Gly ester, Met-Gly peptide), the other
                   ubiquitin-like modifiers (ISG15, NEDD8, UFM1, ATG12) and the features marked "(Microbial infection)"
  sumoylation      the same isopeptide "... interchain with G-Cter in SUMO", SUMO1, SUMO2, SUMO3, SUMO1P1/SUMO5 and the
                   combinations of them: state sumo (the paralog is not kept: iPTMnet and CPLM do not give it), residue K
  methylarginine   Omega-N-methylarginine me1; Asymmetric dimethylarginine me2a; Symmetric dimethylarginine me2s;
                   Dimethylated arginine me2; Omega-N-methylated arginine me                    residue R
                   (a cross-check only: the methylarginine release has its own table)

A residue carrying several states has one row per state. Keys: `canonical_UniProtAC` (an entry of --fasta) and
`position` on its sequence. Columns of every table:

  canonical_UniProtAC, position, residue, mark, state,
  evidence          the ECO codes of the feature, sorted, ';'-joined
  evidence_class    experimental: ECO:0000269 (published experiment) or ECO:0007744 (combinatorial computational and
                    experimental evidence, i.e. large-scale proteomics) among the codes. by-similarity: every other
                    feature (ECO:0000250 by similarity, and ECO:0000255, 0000303, 0000305 or no code)
  pmids             the PubMed IDs of the feature's evidence, sorted, ';'-joined
  source            UniProt
  source_release    the UniProt release served

Dropped and counted in the summary JSON: not_in_reference (accession absent from --fasta), sequence_differs (entry
sequence is not the reference's), on_an_isoform (feature on an isoform's sequence), not_single_residue (a range or
uncertain position), residue_mismatch (the reference has another residue at the position: `common.residue_is`).
Outputs are sorted, so running again on the same download gives the same bytes. Every download is listed in
inputs/marks/raw_files.tsv and each staged file has its row in inputs/MANIFEST.tsv (redistribution git: CC BY 4.0).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter
from datetime import date
from pathlib import Path
from urllib.request import Request, urlopen

import pandas as pd

from common import md5_file, parse_fasta, refresh_manifest_row, residue_is, update_manifest, upsert_tsv_rows

URL = ('https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:9606'
       '&fields=accession,ft_mod_res,sequence&format=json')
CITATION = 'UniProt Consortium, Nucleic Acids Res 53:D609-D617 (2025), doi:10.1093/nar/gkae1010'
LICENSE = 'CC BY 4.0'
# feature description (before ';') -> (mark, state, residue)
ACYL_TERM = re.compile(r'N6-\(?([^()]*?(?:yl|oyl))\)?lysine')
NOT_ACYL = {'methyl', 'ADP-ribosyl', 'lipoyl', 'biotinyl'}
TERMS = {
    'N6-methyllysine': ('methyllysine', 'me1', 'K'),
    'N6,N6-dimethyllysine': ('methyllysine', 'me2', 'K'),
    'N6,N6,N6-trimethyllysine': ('methyllysine', 'me3', 'K'),
    'N6-methylated lysine': ('methyllysine', 'me', 'K'),
    'Phosphoserine': ('phosphorylation', 'phospho', 'S'),
    'Phosphothreonine': ('phosphorylation', 'phospho', 'T'),
    'Phosphotyrosine': ('phosphorylation', 'phospho', 'Y'),
    'Omega-N-methylarginine': ('methylarginine', 'me1', 'R'),
    'Asymmetric dimethylarginine': ('methylarginine', 'me2a', 'R'),
    'Symmetric dimethylarginine': ('methylarginine', 'me2s', 'R'),
    'Dimethylated arginine': ('methylarginine', 'me2', 'R'),
    'Omega-N-methylated arginine': ('methylarginine', 'me', 'R'),
}
CROSSLINK_MARKS = {'ubiquitination', 'sumoylation'}
EXPERIMENTAL = {'ECO:0000269', 'ECO:0007744'}
CROSSLINK = re.compile(r'Glycyl lysine isopeptide \(Lys-Gly\) \(interchain with G-Cter in ([^;)]+)')
URL_CROSSLNK = ('https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:9606'
                '&fields=accession,ft_crosslnk,sequence&format=json')
COLUMNS = ['canonical_UniProtAC', 'position', 'residue', 'mark', 'state', 'evidence', 'evidence_class', 'pmids',
           'source', 'source_release']
def acyl_type(term: str):
    """The acyl type of an 'N6-<acyl>lysine' term (acetyl, succinyl, 2-hydroxyisobutyryl, ...), else None."""
    m = ACYL_TERM.fullmatch(term)
    return m.group(1) if m and m.group(1) not in NOT_ACYL and 'methyl' not in m.group(1) else None


def crosslink_term(description: str):
    """(mark, state, residue) of a "Cross-link" feature that is ubiquitin or SUMO on a lysine (see the module docstring), else None."""
    m = CROSSLINK.match(description)
    if not m:
        return None
    modifier = m.group(1).strip()
    if modifier == 'ubiquitin':
        return 'ubiquitination', 'ub', 'K'
    if modifier.startswith('SUMO'):
        return 'sumoylation', 'sumo', 'K'
    return None


OUT_NAME = {'methyllysine': 'uniprot_{rel}_methyllysine_sites.tsv', 'phosphorylation': 'uniprot_{rel}_phosphorylation_sites.tsv',
            'lysine_acylation': 'uniprot_{rel}_lysine_acylation_sites.tsv',
            'methylarginine': 'uniprot_{rel}_methylarginine_crosscheck.tsv',
            'ubiquitination': 'uniprot_{rel}_ubiquitination_sites.tsv', 'sumoylation': 'uniprot_{rel}_sumoylation_sites.tsv'}
CONTENT = {
    'ubiquitination': 'UniProt lysine ubiquitination annotations (cross-link "Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in ubiquitin)") of the reviewed human entries, one row per residue, with ECO evidence codes, their class and PubMed IDs',
    'sumoylation': 'UniProt lysine SUMOylation annotations (cross-link "Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in SUMO", any paralog, paralog not kept) of the reviewed human entries, one row per residue, with ECO evidence codes, their class and PubMed IDs',
    'lysine_acylation': 'UniProt N6-acyl lysine annotations (acetyl, succinyl, 2-hydroxyisobutyryl, beta-hydroxybutyryl, lactoyl, glutaryl, crotonyl, methacryl, isonicotinyl, benzoyl, butyryl, malonyl, propionyl: every N6-acyl term of the download) of the reviewed human entries, one row per residue and acyl type, with ECO evidence codes, their class and PubMed IDs',
    'methyllysine': 'UniProt N6-methyllysine, N6,N6-dimethyllysine and N6,N6,N6-trimethyllysine (and N6-methylated lysine, degree not stated) annotations of the reviewed human entries, one row per residue and state, with ECO evidence codes, their class (experimental or by-similarity) and PubMed IDs',
    'phosphorylation': 'UniProt phosphoserine, phosphothreonine and phosphotyrosine annotations of the reviewed human entries, one row per residue, with ECO evidence codes, their class (experimental or by-similarity) and PubMed IDs',
    'methylarginine': 'UniProt arginine methylation annotations (omega-N-methylarginine, asymmetric and symmetric dimethylarginine, dimethylated and omega-N-methylated arginine) of the reviewed human entries: a cross-check only, the methylarginine release has its own table',
}


def evidence_class(codes) -> str:
    """experimental if any code is ECO:0000269 or ECO:0007744, else by-similarity (ECO:0000250 and every other code)."""
    return 'experimental' if any(c in EXPERIMENTAL for c in codes) else 'by-similarity'


def feature_term(description: str) -> str:
    return description.split(';')[0].strip()


def stage(entries: list[dict], reference: dict[str, str], release: str) -> tuple[dict[str, pd.DataFrame], dict]:
    """The per-mark tables and the counts of what was dropped, from the `results` of the UniProt JSON."""
    rows, dropped = [], Counter()
    for entry in entries:
        acc = entry['primaryAccession']
        sequence = entry.get('sequence', {}).get('value', '')
        for feature in entry.get('features', []):
            kind = feature.get('type')
            if kind == 'Cross-link':
                term = crosslink_term(feature.get('description', ''))
            elif kind == 'Modified residue':
                name = feature_term(feature.get('description', ''))
                term = TERMS.get(name)
                if term is None and acyl_type(name):
                    term = ('lysine_acylation', acyl_type(name), 'K')
            else:
                continue
            if term is None:
                continue
            mark, state, residue = term
            loc = feature['location']
            if 'sequence' in loc:
                dropped['on_an_isoform'] += 1
                continue
            start, end = loc['start'], loc['end']
            if start['value'] != end['value'] or start['modifier'] != 'EXACT' or end['modifier'] != 'EXACT':
                dropped['not_single_residue'] += 1
                continue
            if acc not in reference:
                dropped['not_in_reference'] += 1
                continue
            if sequence != reference[acc]:
                dropped['sequence_differs'] += 1
                continue
            pos = int(start['value'])
            if not residue_is(reference, acc, pos, residue):
                dropped['residue_mismatch'] += 1
                continue
            evidences = feature.get('evidences', [])
            codes = sorted({e['evidenceCode'] for e in evidences})
            pmids = sorted({e['id'] for e in evidences if e.get('source') == 'PubMed' and e.get('id')}, key=int)
            rows.append({'canonical_UniProtAC': acc, 'position': pos, 'residue': residue, 'mark': mark, 'state': state,
                         'evidence': ';'.join(codes), 'evidence_class': evidence_class(codes), 'pmids': ';'.join(pmids),
                         'source': 'UniProt', 'source_release': release})
    table = pd.DataFrame(rows, columns=COLUMNS)
    # a feature listed twice for the residue and state (two notes) is one row; evidence and PMIDs are merged
    table = (table.groupby(['canonical_UniProtAC', 'position', 'residue', 'mark', 'state', 'source', 'source_release'], as_index=False)
             .agg(evidence=('evidence', lambda v: ';'.join(sorted({c for x in v for c in x.split(';') if c}))),
                  pmids=('pmids', lambda v: ';'.join(sorted({c for x in v for c in x.split(';') if c}, key=int)))))
    table['evidence_class'] = [evidence_class(e.split(';')) for e in table['evidence']]
    table = table[COLUMNS].sort_values(['canonical_UniProtAC', 'position', 'state']).reset_index(drop=True)
    return {m: table[table['mark'] == m].reset_index(drop=True) for m in OUT_NAME}, dict(dropped)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fasta', required=True, help='the reference proteome (canonical reviewed accessions)')
    ap.add_argument('--raw', help='read this "Modified residue" download (and its .release file) instead of fetching')
    ap.add_argument('--raw-crosslnk', help='read this "Cross-link" download (and its .release file) instead of fetching '
                                           '(default with --raw: the --raw name with mod_res replaced by crosslnk)')
    ap.add_argument('--raw-dir', default='PTM_data/marks', help='where the download is kept (not in Git)')
    ap.add_argument('--out-dir', default='inputs/marks')
    ap.add_argument('--manifest', default='inputs/MANIFEST.tsv')
    ap.add_argument('--downloaded', default='', help='date of the download for raw_files.tsv (default: today when fetching)')
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    downloads = {}    # stem -> (path, release, url)
    if args.raw:
        raw_cross = Path(args.raw_crosslnk) if args.raw_crosslnk else Path(args.raw.replace('mod_res', 'crosslnk'))
        for stem, path, url in (('mod_res', Path(args.raw), URL), ('crosslnk', raw_cross, URL_CROSSLNK)):
            if not path.exists():
                raise SystemExit(f'{path} not found (--raw gives the "Modified residue" download, --raw-crosslnk the "Cross-link" one)')
            downloads[stem] = (path, path.with_suffix('.release').read_text().strip(), url)
    else:
        for stem, url in (('mod_res', URL), ('crosslnk', URL_CROSSLNK)):
            with urlopen(Request(url, headers={'User-Agent': 'PTM-pipeline-stager/1.0'}), timeout=1800) as response:
                release = response.headers.get('X-UniProt-Release', '')
                if not release:
                    raise SystemExit('the UniProt API gave no X-UniProt-Release header')
                data = response.read()
            path = Path(args.raw_dir) / f'uniprot_human_reviewed_{stem}_{release}.json'
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
            path.with_suffix('.release').write_text(release + '\n')
            downloads[stem] = (path, release, url)
    release = downloads['mod_res'][1]
    if downloads['crosslnk'][1] != release:
        raise SystemExit(f'the two downloads are of different UniProt releases: {release} and {downloads["crosslnk"][1]}')
    entries = []
    for path, _, _ in downloads.values():
        entries += json.loads(path.read_text(encoding='utf-8'))['results']
    reference = parse_fasta(args.fasta)

    tables, dropped = stage(entries, reference, release)
    out_dir.mkdir(parents=True, exist_ok=True)
    files = {}
    for mark, table in tables.items():
        name = OUT_NAME[mark].format(rel=release)
        table.to_csv(out_dir / name, sep='\t', index=False, lineterminator='\n')
        files[mark] = name
    raws = {stem: {'path': path, 'md5': md5_file(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'url': url}
            for stem, (path, _, url) in downloads.items()}
    of_mark = {m: 'crosslnk' if m in CROSSLINK_MARKS else 'mod_res' for m in OUT_NAME}    # the download a mark's rows come from
    summary = {
        'source': 'UniProtKB/Swiss-Prot reviewed human, "Modified residue" and "Cross-link" features (UniProt REST API)',
        'url': URL, 'url_crosslnk': URL_CROSSLNK,
        'uniprot_release_served': release,
        'reference_fasta': args.fasta,
        'reference_fasta_md5': md5_file(args.fasta),
        'raw_file': raws['mod_res']['path'].name,
        'raw_md5': raws['mod_res']['md5'],
        'raw_crosslnk_file': raws['crosslnk']['path'].name,
        'raw_crosslnk_md5': raws['crosslnk']['md5'],
        'entries_in_download': len(entries),
        'dropped': dropped,
        'rows': {m: len(t) for m, t in tables.items()},
        'rows_by_state': {m: {k: int(v) for k, v in t['state'].value_counts().sort_index().items()} for m, t in tables.items()},
        'rows_by_evidence_class': {m: {k: int(v) for k, v in t['evidence_class'].value_counts().sort_index().items()}
                                   for m, t in tables.items()},
        'files': files,
        'files_md5': {m: md5_file(out_dir / f) for m, f in files.items()},
    }
    summary_name = 'uniprot_ptm_sites_summary.json'
    (out_dir / summary_name).write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n', encoding='utf-8')

    # raw_files.tsv: a row that is already there for the same bytes keeps its download date
    raw_list = out_dir / 'raw_files.tsv'
    head = 'file\turl\tbytes\tmd5\tsha256\tdownloaded\tversion\tkept'
    if not raw_list.exists():
        raw_list.write_text(head + '\n', encoding='utf-8')
    old = pd.read_csv(raw_list, sep='\t', dtype=str, keep_default_na=False)
    raw_rows, when = [], {}
    for stem, r in raws.items():
        same = old[(old['file'] == r['path'].name) & (old['md5'] == r['md5'])]
        when[stem] = same['downloaded'].iloc[0] if len(same) else (args.downloaded or date.today().isoformat())
        raw_rows.append({'file': r['path'].name, 'url': r['url'], 'bytes': r['path'].stat().st_size, 'md5': r['md5'],
                         'sha256': r['sha256'], 'downloaded': when[stem], 'version': f'UniProt release {release} (X-UniProt-Release header)',
                         'kept': f'release_asset (Zenodo archive of the release; {round(r["path"].stat().st_size / 1e6)} MB, not in Git)'})
    upsert_tsv_rows(raw_list, raw_rows)

    command = (f'python scripts/54_stage_uniprot_ptm_sites.py --fasta {args.fasta} '
               f'--raw-dir PTM_data/marks --out-dir inputs/marks')
    rows = []
    for mark, name in files.items():
        path = out_dir / name
        rows.append({'file': f'marks/{name}', 'bytes': path.stat().st_size, 'md5': md5_file(path), 'rows': len(tables[mark]),
                     'redistribution': 'git', 'content': CONTENT[mark] + f': {len(tables[mark]):,} rows',
                     'source': raws[of_mark[mark]]['url'],
                     'version_or_date': f'UniProt release {release} (X-UniProt-Release header); downloaded {when[of_mark[mark]]}',
                     'citation': CITATION, 'licence': LICENSE,
                     'derived_from': f'{raws[of_mark[mark]]["path"].name} (MD5 {raws[of_mark[mark]]["md5"]})',
                     'how_made_or_checked': command + '; residue checked on the reference sequence, entries whose sequence differs dropped; counts in marks/uniprot_ptm_sites_summary.json'})
    path = out_dir / summary_name
    rows.append({'file': f'marks/{summary_name}', 'bytes': path.stat().st_size, 'md5': md5_file(path),
                 'redistribution': 'git', 'content': 'Counts, release served, drops and the MD5 of the input and of every table of the UniProt PTM-site staging',
                 'source': 'the tables above', 'version_or_date': f'UniProt release {release}', 'derived_from': 'marks/raw_files.tsv',
                 'how_made_or_checked': command})
    update_manifest(args.manifest, rows)
    # raw_files.tsv has its own row, refreshed last (it holds the MD5 of every raw download)
    refresh_manifest_row(args.manifest, raw_list, 'marks/raw_files.tsv', redistribution='git',
                         content='URL, size, MD5, SHA-256 and date of every raw download behind inputs/marks/',
                         source='the URLs it lists', derived_from='', how_made_or_checked='written by the scripts/54 and scripts/55 staging runs')
    print(json.dumps({k: summary[k] for k in ('uniprot_release_served', 'dropped', 'rows', 'rows_by_state', 'rows_by_evidence_class')}, indent=1))


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Stage every UniProt modified residue (MOD_RES) of the reviewed human entries, with the enzyme UniProt names.

    python scripts/53_stage_uniprot_mod_res.py --release 2026_03 --outdir STAGE [--raw uniprot_mod_res.tsv]
    python scripts/apply_staged.py STAGE

The download, the release check and the feature parsing are scripts/00_stage_uniprot_methylarginine.py's (one REST
stream of accession, MOD_RES and sequence; it stops unless the API serves --release); this keeps every term rather than
the arginine methylation ones. A feature is kept where its entry's sequence equals the reference one (--fasta) and the
residue the term names is there; features on isoform sequences and on other sequences are counted and left out.

Each row: uniprot_accession, position, residue, term (UniProt's controlled vocabulary, before the first ';'),
ptm_class (ptm_class() below; 'other' for the rest), enzyme (the 'by X' parts of the note, e.g. CK2, PKA, autocatalysis;
';'-joined), alternate (UniProt's 'alternate': another modification is annotated at the same residue), evidence (ECO
codes), evidence_class (experimental = ECO:0000269 or ECO:0007744, by_similarity = ECO:0000250, other).

Writes layers/uniprot_2026_03_mod_res.tsv, uniprot_mod_res_summary.json, manifest_rows.tsv and raw_files_rows.tsv.
UniProt: CC BY 4.0.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import re
import sys
from collections import Counter
from pathlib import Path
from urllib.request import Request, urlopen

import pandas as pd

from common import md5_file, parse_fasta

CITATION = 'UniProt Consortium, Nucleic Acids Res 53:D609-D617 (2025), doi:10.1093/nar/gkae1010'
LICENCE = 'CC BY 4.0'
OUT_SUMMARY = 'uniprot_mod_res_summary.json'

# term (UniProt ptmlist) -> class, tested in this order; anything else is 'other'
ACYL_STEMS = ('succinyl', 'crotonyl', 'lactoyl', 'malonyl', 'glutaryl', 'hydroxyisobutyryl', 'hydroxybutyryl', 'butyryl',
              'propionyl', 'benzoyl', 'formyl', 'methacryl', 'isonicotinyl')
METHYL_R = ('Omega-N-methylarginine', 'Asymmetric dimethylarginine', 'Symmetric dimethylarginine', 'Dimethylated arginine',
            'Omega-N-methylated arginine')
METHYL_K = ('N6-methyllysine', 'N6,N6-dimethyllysine', 'N6,N6,N6-trimethyllysine')
CK2 = re.compile(r'\bCK2\b')


# the residue each class must sit on (others are checked against the residue the term names, when it names one)
CLASS_RESIDUES = {'phospho': 'STYHDC', 'acetyl_K': 'K', 'acyl_K': 'K', 'methyl_R': 'R', 'methyl_K': 'K', 'glutamyl': 'E',
                  'glycyl': 'E', 'citrulline': 'R'}


def load_stager():
    spec = importlib.util.spec_from_file_location('stage_uniprot_methylarginine',
                                                  Path(__file__).with_name('00_stage_uniprot_methylarginine.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def ptm_class(term: str) -> str:
    if 'ADP-ribosyl' in term:
        return 'adp_ribosyl'
    if term.startswith('Phospho'):
        return 'phospho'
    if term == 'N6-acetyllysine':
        return 'acetyl_K'
    if term.startswith('N6-') and term.endswith('lysine') and any(x in term for x in ACYL_STEMS):
        return 'acyl_K'
    if term.startswith(('N-acetyl', 'N2-acetyl')):
        return 'acetyl_Nterm'
    if term in METHYL_R:
        return 'methyl_R'
    if term in METHYL_K:
        return 'methyl_K'
    if term == '5-glutamyl polyglutamate':
        return 'glutamyl'
    if term in ('5-glutamyl glycine', 'Polyglycine'):
        return 'glycyl'
    if term == 'Citrulline':
        return 'citrulline'
    return 'other'


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--fasta', default='inputs/reference/uniprot_human_reviewed_2026_03.fasta')
    ap.add_argument('--release', required=True, help='the UniProt release the reference is, e.g. 2026_03')
    ap.add_argument('--raw', help='read this download instead of fetching')
    ap.add_argument('--manifest', default='inputs/MANIFEST.tsv', help='read only: its columns')
    ap.add_argument('--raw-files', default='inputs/layers/raw_files.tsv', help='read only: its columns')
    ap.add_argument('--downloaded', required=True, help='the day of the download, YYYY-MM-DD')
    ap.add_argument('--outdir', required=True, type=Path)
    a = ap.parse_args()
    st = load_stager()
    out = a.outdir
    out.mkdir(parents=True, exist_ok=True)
    raw_name = f'uniprot_human_reviewed_mod_res_{a.release}.tsv'
    raw_path = Path(a.raw) if a.raw else out / raw_name
    if not a.raw:
        with urlopen(Request(st.URL, headers={'User-Agent': 'PTM-pipeline-stager/1.0'}), timeout=900) as response:
            release = response.headers.get('X-UniProt-Release', '')
            if release != a.release:
                sys.exit(f'the UniProt API serves release {release!r}, not {a.release}')
            raw_path.write_bytes(response.read())
    raw = pd.read_csv(raw_path, sep='\t', dtype=str, keep_default_na=False)
    reference = parse_fasta(a.fasta)

    rows, left = [], Counter()
    for acc, features, sequence in raw.itertuples(index=False):
        left['on_an_isoform'] += len(st.ISOFORM_FEATURE.findall(features))
        feats = st.FEATURE.findall(features)
        if not feats:
            continue
        if acc not in reference:
            left['not_in_reference'] += len(feats)
            continue
        if sequence != reference[acc]:
            left['sequence_differs'] += len(feats)
            continue
        for pos, note, ev in feats:
            pos = int(pos)
            parts = [p.strip() for p in note.split(';')]
            term = parts[0]
            cls = ptm_class(term)
            residue = sequence[pos - 1] if 1 <= pos <= len(sequence) else ''
            if not residue or (cls in CLASS_RESIDUES and residue not in CLASS_RESIDUES[cls]):
                left['residue_not_as_named'] += 1
                continue
            codes = [e.split('|')[0] for e in ev.split(', ') if e]
            rows.append({'uniprot_accession': acc, 'position': pos, 'residue': residue, 'term': term, 'ptm_class': cls,
                         'enzyme': ';'.join(p[3:] for p in parts[1:] if p.startswith('by ')),
                         'alternate': 'alternate' in parts[1:], 'evidence': ';'.join(sorted(set(codes))),
                         'evidence_class': st.evidence_class(codes)})
    tab = pd.DataFrame(rows).sort_values(['uniprot_accession', 'position', 'term']).reset_index(drop=True)
    out_name = f'uniprot_{a.release}_mod_res.tsv'
    tab.to_csv(out / out_name, sep='\t', index=False)
    ck2 = tab[(tab.ptm_class == 'phospho') & tab.enzyme.str.contains(CK2)]
    summary = {
        'release': a.release, 'script': 'scripts/53_stage_uniprot_mod_res.py', 'raw': {'file': raw_name, 'md5': md5_file(raw_path)},
        'reference': {'fasta': a.fasta, 'md5': md5_file(a.fasta)},
        'rows': len(tab), 'entries': int(tab.uniprot_accession.nunique()),
        'by_class': {k: int(v) for k, v in tab.ptm_class.value_counts().items()},
        'by_evidence_class': {k: int(v) for k, v in tab.evidence_class.value_counts().items()},
        'phospho_with_enzyme': int(((tab.ptm_class == 'phospho') & (tab.enzyme != '')).sum()),
        'ck2_sites': int(len(ck2)), 'ck2_entries': int(ck2.uniprot_accession.nunique()),
        'top_enzymes': dict(Counter(e for v in tab[tab.ptm_class == 'phospho'].enzyme for e in v.split(';') if e).most_common(25)),
        'left_out': dict(left),
    }
    (out / OUT_SUMMARY).write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')

    mcols = open(a.manifest, encoding='utf-8').readline().rstrip('\n').split('\t')
    rcols = open(a.raw_files, encoding='utf-8').readline().rstrip('\n').split('\t')
    source = ('Derived by scripts/53_stage_uniprot_mod_res.py from the UniProt REST API: '
              'https://rest.uniprot.org/uniprotkb/stream?query=(organism_id:9606) AND (reviewed:true) AND (ft_mod_res:*)'
              '&fields=accession,ft_mod_res,sequence&format=tsv')
    version = f'UniProt release {a.release} (X-UniProt-Release header); downloaded {a.downloaded}'
    how = (f'scripts/53_stage_uniprot_mod_res.py --release {a.release}: features of entries whose sequence equals '
           'reference/uniprot_human_reviewed_2026_03.fasta, on the residue the term names')
    by = summary['by_class']
    mrows = [
        {'file': f'layers/{out_name}', 'rows': len(tab), 'derived_from': f'{raw_name} (MD5 {summary["raw"]["md5"]})',
         'content': (f'Every UniProt {a.release} modified residue (MOD_RES) of the reviewed human entries, with its class, the '
                     f'enzyme UniProt names and the evidence codes: {len(tab):,} annotations on {summary["entries"]:,} entries '
                     f'({by.get("phospho", 0):,} phospho, {summary["ck2_sites"]} of them by CK2; {by.get("acetyl_K", 0):,} N6-acetyllysine; '
                     f'{by.get("glutamyl", 0)} polyglutamylation; {by.get("adp_ribosyl", 0)} ADP-ribosylation; {by.get("citrulline", 0)} citrulline)'),
         'how_made_or_checked': how},
        {'file': f'layers/{OUT_SUMMARY}', 'rows': '', 'derived_from': '', 'how_made_or_checked': 'scripts/53_stage_uniprot_mod_res.py',
         'content': 'Release, input MD5s, counts by class and evidence, CK2 counts and what was left out of the UniProt MOD_RES staging'},
    ]
    for r in mrows:
        p = out / r['file'].split('/')[-1]
        r.update({'bytes': p.stat().st_size, 'md5': md5_file(p), 'redistribution': 'git', 'source': source,
                  'version_or_date': version, 'citation': CITATION, 'licence': LICENCE, 'original_name': raw_name})
    pd.DataFrame([{c: r.get(c, '') for c in mcols} for r in mrows]).to_csv(out / 'manifest_rows.tsv', sep='\t', index=False)
    rrows = [{'file': raw_name, 'url': st.URL, 'bytes': raw_path.stat().st_size, 'md5': md5_file(raw_path),
              'sha256': sha256_file(raw_path), 'downloaded': a.downloaded, 'version': f'UniProt release {a.release}',
              'kept': 'zenodo'}]
    pd.DataFrame([{c: r.get(c, '') for c in rcols} for r in rrows]).to_csv(out / 'raw_files_rows.tsv', sep='\t', index=False)
    print(json.dumps({k: summary[k] for k in ('rows', 'entries', 'by_class', 'ck2_sites', 'ck2_entries', 'left_out')}, indent=1))


if __name__ == '__main__':
    main()

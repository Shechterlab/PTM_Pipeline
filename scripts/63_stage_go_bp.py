#!/usr/bin/env python3
"""Stage the Gene Ontology Biological Process annotations of the reviewed human proteins (inputs/annotation/).

    python -I scripts/63_stage_go_bp.py --gaf PTM_data/annotation/goa_human.gaf.gz --obo PTM_data/annotation/go-basic.obo \
        --fasta inputs/reference/uniprot_human_reviewed_2026_03.fasta \
        --history inputs/reference/uniprot_accession_history_2026_03.tsv --out-dir inputs/annotation --downloaded 2026-11-13 --go-build 2026-08-05

Sources (Gene Ontology Consortium; GO and its annotations are CC BY 4.0, https://geneontology.org/docs/go-citation-policy/):
  https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz   the EBI GOA human annotation file (GAF 2.2; its
      header gives date-generated and the GO version it was checked against); the current file has no versioned URL,
      so the download goes into the release's Zenodo archive
  https://current.geneontology.org/ontology/go-basic.obo              go-basic (no cross-ontology links, acyclic over
      is_a and part_of), for term names and the graph; the same file stays at the permanent release URL (GO_RELEASE_URL)

What is kept
  * GAF rows of aspect P (biological process), DB UniProtKB, whose accession is an entry of the reference proteome
    (--fasta), or a secondary accession of exactly one such entry (--history), which then takes its place;
  * every qualifier but those with NOT (involved_in, acts_upstream_of and its variants); NOT rows are dropped and
    counted, and a positive row is kept even where the same protein has a NOT row to the same term (counted);
  * evidence ND (no biological data, an annotation to the root that says nothing is known) is dropped and counted;
  * GO ids are read through go-basic: an alt_id stands for its term; ids that go-basic does not have, or has as
    obsolete (the GAF may be checked against a newer GO than the current go-basic), are dropped and listed in the summary.

Propagation (true-path rule): every direct term gives the protein all its ancestors over is_a and part_of within
biological_process (regulates, positively_regulates and negatively_regulates are not followed), the root GO:0008150
included.

Writes, named by the GAF's date-generated
  go_bp_human_<date>.tsv        one row per protein and BP term (direct or inherited):
      uniprot_accession       reviewed accession (UniProt 2026_03 reference)
      go_id                   the term (its name is in go_bp_terms_<date>.tsv; left out here to keep the table at
                              about half the size)
      annotation              direct (the protein is annotated to this term itself) or propagated (only to descendants)
      direct_evidence_codes   evidence codes of the direct annotations to this term (';'-joined, sorted; empty when propagated)
      evidence_codes          evidence codes of all annotations at or below this term that give it to the protein
      non_iea                 true when evidence_codes has a code other than IEA: the term survives dropping IEA
  go_bp_terms_<date>.tsv        one row per BP term with at least one protein:
      go_id, go_name
      depth_min, depth_max    shortest and longest is_a/part_of path to the root (root 0)
      n_proteins              proteins with the term, propagated, any evidence
      n_proteins_non_iea      the same without IEA (non_iea true)
      n_proteins_direct       proteins annotated to the term itself
      parents                 is_a and part_of parents (';'-joined)
  go_bp_human_<date>_summary.json  counts, versions, drops, the IEA share and the MD5 of every input and output
Outputs are sorted, so running again on the same downloads gives the same bytes.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import md5_file, parse_fasta, refresh_manifest_row, update_manifest  # noqa: E402

GAF_URL = 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz'
OBO_URL = 'https://current.geneontology.org/ontology/go-basic.obo'
GO_RELEASE_URL = 'https://release.geneontology.org/{build}/ontology/go-basic.obo'
ROOT = 'GO:0008150'
NAMESPACE = 'biological_process'
EDGES = ('is_a', 'part_of')
CITATION = ('Gene Ontology Consortium, Genetics 224:iyad031 (2023), doi:10.1093/genetics/iyad031; '
            'Ashburner et al., Nat Genet 25:25 (2000), doi:10.1038/75556')
LICENCE = 'CC BY 4.0 (Gene Ontology and GO annotations, https://geneontology.org/docs/go-citation-policy/)'
ROW_COLUMNS = ['uniprot_accession', 'go_id', 'annotation', 'direct_evidence_codes', 'evidence_codes', 'non_iea']
TERM_COLUMNS = ['go_id', 'go_name', 'depth_min', 'depth_max', 'n_proteins', 'n_proteins_non_iea', 'n_proteins_direct', 'parents']


# ----------------------------------------------------------------------------- ontology
def read_obo(lines) -> tuple[dict, dict]:
    """Terms of an OBO file: {id: {name, namespace, obsolete, parents (is_a and part_of), alt_ids}}, and the header tags."""
    header, terms, cur, in_term = {}, {}, None, False
    for raw in lines:
        line = raw.rstrip('\n')
        if line.startswith('['):
            in_term = line == '[Term]'
            cur = None
            continue
        if not line or line.startswith('!'):
            continue
        tag, _, value = line.partition(': ')
        if cur is None and not in_term:
            header.setdefault(tag, value)
            continue
        value = value.split(' ! ')[0].strip()
        if tag == 'id':
            cur = terms.setdefault(value, {'name': '', 'namespace': '', 'obsolete': False, 'parents': set(), 'alt_ids': []})
        elif cur is None:
            continue
        elif tag == 'name':
            cur['name'] = value
        elif tag == 'namespace':
            cur['namespace'] = value
        elif tag == 'is_obsolete':
            cur['obsolete'] = value == 'true'
        elif tag == 'alt_id':
            cur['alt_ids'].append(value)
        elif tag == 'is_a':
            cur['parents'].add(value.split()[0])
        elif tag == 'relationship':
            rel, target = value.split()[:2]
            if rel in EDGES:
                cur['parents'].add(target)
    return terms, header


def bp_graph(terms: dict) -> tuple[dict, dict]:
    """parents of each live BP term (BP parents only), and the alt_id -> id map."""
    bp = {t for t, v in terms.items() if v['namespace'] == NAMESPACE and not v['obsolete']}
    parents = {t: sorted(p for p in terms[t]['parents'] if p in bp) for t in bp}
    alt = {a: t for t, v in terms.items() for a in v['alt_ids']}
    return parents, alt


def ancestors_of(parents: dict) -> dict:
    """{term: frozenset of its ancestors} (not the term itself)."""
    memo: dict[str, frozenset] = {}

    def anc(t, stack=()):
        if t in memo:
            return memo[t]
        if t in stack:
            raise SystemExit(f'cycle over is_a/part_of at {t}')
        out = set()
        for p in parents[t]:
            out.add(p)
            out |= anc(p, stack + (t,))
        memo[t] = frozenset(out)
        return memo[t]

    sys.setrecursionlimit(10000)
    for t in parents:
        anc(t)
    return memo


def depths(parents: dict) -> dict:
    """{term: (shortest, longest) path length to the root}."""
    memo: dict[str, tuple[int, int]] = {}

    def d(t):
        if t not in memo:
            if t == ROOT:
                memo[t] = (0, 0)
            elif not parents[t]:
                raise SystemExit(f'{t} has no path to {ROOT}')
            else:
                ds = [d(p) for p in parents[t]]
                memo[t] = (1 + min(x for x, _ in ds), 1 + max(y for _, y in ds))
        return memo[t]

    for t in parents:
        d(t)
    return memo


# ----------------------------------------------------------------------------- annotations
def read_secondary(path, reference) -> dict:
    """secondary accession -> the one reference entry that lists it."""
    to = defaultdict(set)
    with open(path, encoding='utf-8') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['kind'] == 'secondary' and r['accession'] in reference:
                to[r['key']].add(r['accession'])
    return {k: next(iter(v)) for k, v in to.items() if len(v) == 1 and k not in reference}


def read_gaf(lines, reference: set, secondary: dict, terms: dict, parents: dict, alt: dict):
    """Positive BP annotations: {(accession, term): set of evidence codes}, and counts."""
    direct = defaultdict(set)
    negative = set()
    counts = Counter()
    accession_status = {}
    unknown = Counter()
    header = {}
    for line in lines:
        if line.startswith('!'):
            m = re.match(r'!\s*([\w-]+):\s*(.*)', line.strip())
            if m:
                header.setdefault(m.group(1), m.group(2))
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 15:
            counts['gaf_rows_malformed'] += 1
            continue
        if f[8] != 'P':
            counts['gaf_rows_other_aspect'] += 1
            continue
        counts['gaf_rows_bp'] += 1
        if f[0] != 'UniProtKB':
            counts['bp_rows_not_uniprotkb'] += 1
            continue
        acc = f[1]
        if acc in reference:
            accession_status[acc] = 'reference'
        elif acc in secondary:
            accession_status[acc] = 'secondary'
            acc = secondary[acc]
        else:
            accession_status[acc] = 'not_in_reference'
            counts['bp_rows_not_in_reference'] += 1
            continue
        go = alt.get(f[4], f[4])
        if go != f[4]:
            counts['bp_rows_alt_id_mapped'] += 1
        if go not in parents:
            why = 'obsolete' if go in terms and terms[go]['obsolete'] else 'not_bp' if go in terms else 'not_in_go_basic'
            unknown[(f[4], why)] += 1
            counts[f'bp_rows_term_{why}'] += 1
            continue
        if f[6] == 'ND':
            counts['bp_rows_nd_evidence'] += 1
            continue
        if 'NOT' in f[3].split('|'):
            counts['bp_rows_not_qualifier'] += 1
            negative.add((acc, go))
            continue
        counts['bp_rows_kept'] += 1
        counts['bp_rows_kept_iea'] += f[6] == 'IEA'
        counts[f'qualifier:{f[3]}'] += 1
        direct[(acc, go)].add(f[6])
    counts['positive_with_not_row_same_term'] = sum(1 for k in negative if k in direct)
    gaf_accessions = Counter(accession_status.values())
    return direct, counts, unknown, header, gaf_accessions


def propagate(direct: dict, ancestors: dict) -> dict:
    """{accession: {term: [direct codes, all codes]}} over the true-path rule."""
    out: dict[str, dict] = defaultdict(dict)
    for (acc, go), codes in direct.items():
        prot = out[acc]
        for t in (go, *ancestors[go]):
            entry = prot.setdefault(t, [set(), set()])
            entry[1] |= codes
        prot[go][0] |= codes
    return out


def build_rows(prop: dict) -> list[list]:
    rows = []
    for acc in sorted(prop):
        for go in sorted(prop[acc]):
            d, a = prop[acc][go]
            rows.append([acc, go, 'direct' if d else 'propagated', ';'.join(sorted(d)), ';'.join(sorted(a)),
                         'true' if a - {'IEA'} else 'false'])
    return rows


def build_terms(rows: list[list], terms: dict, parents: dict, depth: dict) -> list[list]:
    n, n_ok, n_direct = Counter(), Counter(), Counter()
    for acc, go, ann, _, _, ok in rows:
        n[go] += 1
        n_ok[go] += ok == 'true'
        n_direct[go] += ann == 'direct'
    return [[go, terms[go]['name'], depth[go][0], depth[go][1], n[go], n_ok[go], n_direct[go], ';'.join(parents[go])]
            for go in sorted(n)]


def write_tsv(path: Path, columns, rows) -> None:
    with open(path, 'w', encoding='utf-8', newline='') as fh:
        fh.write('\t'.join(columns) + '\n')
        for r in rows:
            fh.write('\t'.join(str(x) for x in r) + '\n')


def open_text(path: Path):
    return gzip.open(path, 'rt', encoding='utf-8') if path.suffix == '.gz' else open(path, encoding='utf-8')


def digest(path: Path) -> tuple[int, str, str]:
    data = path.read_bytes()
    return len(data), hashlib.md5(data).hexdigest(), hashlib.sha256(data).hexdigest()


# ----------------------------------------------------------------------------- main
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--gaf', required=True, help='goa_human.gaf.gz')
    ap.add_argument('--obo', required=True, help='go-basic.obo')
    ap.add_argument('--fasta', required=True, help='the reference proteome (canonical reviewed accessions)')
    ap.add_argument('--history', required=True, help='uniprot_accession_history_<release>.tsv')
    ap.add_argument('--out-dir', default='inputs/annotation')
    ap.add_argument('--manifest', default='inputs/MANIFEST.tsv')
    ap.add_argument('--downloaded', default='', help='date of the downloads (default: today)')
    ap.add_argument('--go-build', default='', help="GO release build of go-basic (the x-amz-meta-build-number of the download), for its permanent URL")
    args = ap.parse_args()
    gaf, obo, out = Path(args.gaf), Path(args.obo), Path(args.out_dir)
    when = args.downloaded or date.today().isoformat()

    with open_text(obo) as fh:
        terms, obo_header = read_obo(fh)
    parents, alt = bp_graph(terms)
    ancestors, depth = ancestors_of(parents), depths(parents)
    reference = set(parse_fasta(args.fasta))
    secondary = read_secondary(args.history, reference)
    with open_text(gaf) as fh:
        direct, counts, unknown, gaf_header, gaf_accessions = read_gaf(fh, reference, secondary, terms, parents, alt)
    gaf_date = gaf_header.get('date-generated', '').split()[0]
    if not gaf_date:
        raise SystemExit('the GAF header has no date-generated')
    go_release = obo_header.get('data-version', '').replace('releases/', '')

    prop = propagate(direct, ancestors)
    rows = build_rows(prop)
    term_rows = build_terms(rows, terms, parents, depth)

    out.mkdir(parents=True, exist_ok=True)
    f_rows, f_terms = f'go_bp_human_{gaf_date}.tsv', f'go_bp_terms_{gaf_date}.tsv'
    f_summary = f'go_bp_human_{gaf_date}_summary.json'
    write_tsv(out / f_rows, ROW_COLUMNS, rows)
    write_tsv(out / f_terms, TERM_COLUMNS, term_rows)

    # IEA share
    ann_codes = Counter(c for codes in direct.values() for c in codes)
    direct_pairs = len(direct)
    iea_only_direct = sum(1 for codes in direct.values() if codes == {'IEA'})
    iea_only_rows = sum(1 for r in rows if r[5] == 'false')
    proteins = len(prop)
    proteins_non_iea = len({r[0] for r in rows if r[5] == 'true'})
    gaf_n, gaf_md5, gaf_sha = digest(gaf)
    obo_n, obo_md5, obo_sha = digest(obo)
    summary = {
        'source': {'gaf': GAF_URL, 'obo': OBO_URL},
        'citation': CITATION,
        'license': LICENCE,
        'gaf': {'file': gaf.name, 'bytes': gaf_n, 'md5': gaf_md5, 'sha256': gaf_sha, 'downloaded': when,
                'gaf_version': gaf_header.get('gaf-version', ''), 'date_generated': gaf_header.get('date-generated', ''),
                'go_version_checked_against': gaf_header.get('go-version', '')},
        'obo': {'file': obo.name, 'bytes': obo_n, 'md5': obo_md5, 'sha256': obo_sha, 'downloaded': when,
                'data_version': obo_header.get('data-version', ''), 'format_version': obo_header.get('format-version', '')},
        'reference_fasta': args.fasta, 'reference_fasta_md5': md5_file(args.fasta), 'reference_entries': len(reference),
        'accession_history': args.history, 'accession_history_md5': md5_file(args.history),
        'propagation': 'is_a and part_of within biological_process (true-path rule); regulates relations not followed; root included',
        'go_basic_bp_terms_live': len(parents),
        'gaf_rows': dict(sorted(counts.items())),
        'gaf_bp_accessions': dict(sorted(gaf_accessions.items())),
        'dropped_term_ids': {f'{t} ({why})': n for (t, why), n in sorted(unknown.items())},
        'direct_annotations_by_evidence_code': dict(sorted(ann_codes.items())),
        'direct_protein_term_pairs': direct_pairs,
        'direct_pairs_iea_only': iea_only_direct,
        'direct_pairs_iea_only_share': round(iea_only_direct / direct_pairs, 4),
        'bp_rows_kept_iea_share': round(counts['bp_rows_kept_iea'] / counts['bp_rows_kept'], 4),
        'direct_pair_codes_iea_share': round(ann_codes['IEA'] / sum(ann_codes.values()), 4),
        'rows': len(rows),
        'rows_direct': sum(1 for r in rows if r[2] == 'direct'),
        'rows_propagated': sum(1 for r in rows if r[2] == 'propagated'),
        'rows_iea_only': iea_only_rows,
        'rows_iea_only_share': round(iea_only_rows / len(rows), 4),
        'proteins': proteins,
        'proteins_share_of_reference': round(proteins / len(reference), 4),
        'proteins_with_non_iea_term': proteins_non_iea,
        'terms_with_proteins': len(term_rows),
        'files': {'rows': f_rows, 'terms': f_terms},
        'files_md5': {'rows': md5_file(out / f_rows), 'terms': md5_file(out / f_terms)},
    }
    (out / f_summary).write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n', encoding='utf-8')

    # MANIFEST
    folder = out.name
    command = (f'python -I scripts/63_stage_go_bp.py --gaf goa_human.gaf.gz --obo go-basic.obo --fasta {args.fasta} '
               f'--history {args.history} --out-dir inputs/{folder} --downloaded {when}'
               + (f' --go-build {args.go_build}' if args.go_build else ''))
    version = (f'GOA human GAF {gaf_header.get("gaf-version", "")} generated {gaf_header.get("date-generated", "")} '
               f'(checked against GO {gaf_header.get("go-version", "").split("/releases/")[-1].split("/")[0]}); '
               f'go-basic {go_release}; downloaded {when}')
    derived = f'goa_human.gaf.gz (MD5 {gaf_md5}); go-basic.obo (MD5 {obo_md5}); reference/uniprot_human_reviewed_2026_03.fasta'
    common_fields = {'source': f'{GAF_URL}; {OBO_URL}', 'version_or_date': version, 'citation': CITATION, 'licence': LICENCE,
                     'derived_from': derived, 'redistribution': 'git'}
    check = (f'; aspect P, NOT qualifiers ({counts["bp_rows_not_qualifier"]:,} rows) and ND evidence ({counts["bp_rows_nd_evidence"]:,}) dropped, accessions of the 2026_03 reference '
             f'(secondary accessions mapped through the accession history), terms through go-basic (alt_id mapped, '
             f'{sum(unknown.values()):,} rows to {len(unknown):,} obsolete or missing ids dropped), propagated over is_a and part_of; '
             f'counts in {folder}/{f_summary}')
    iea_note = (f'{iea_only_rows:,} of {len(rows):,} rows ({iea_only_rows / len(rows):.1%}) rest on IEA alone (non_iea false); '
                f'{counts["bp_rows_kept_iea"]:,} of {counts["bp_rows_kept"]:,} kept GAF rows ({counts["bp_rows_kept_iea"] / counts["bp_rows_kept"]:.1%}) are IEA; '
                f'{iea_only_direct:,} of {direct_pairs:,} direct protein-term pairs ({iea_only_direct / direct_pairs:.1%}) have IEA only; '
                'term names are in the terms table')
    rows_manifest = [
        {'file': f'{folder}/{f_rows}', 'bytes': (out / f_rows).stat().st_size, 'md5': md5_file(out / f_rows), 'rows': len(rows),
         'content': (f'GO biological process terms of {proteins:,} reviewed human proteins (UniProt 2026_03 reference), one row per '
                     f'protein and term, direct or propagated (true-path rule over is_a and part_of), with evidence codes; {iea_note}'),
         'how_made_or_checked': command + check, **common_fields},
        {'file': f'{folder}/{f_terms}', 'bytes': (out / f_terms).stat().st_size, 'md5': md5_file(out / f_terms), 'rows': len(term_rows),
         'content': (f'The {len(term_rows):,} GO biological process terms with at least one reviewed human protein: name, shortest and '
                     'longest depth below the root, proteins with the term (propagated; without IEA; direct), is_a/part_of parents'),
         'how_made_or_checked': command + check, **common_fields},
        {'file': f'{folder}/{f_summary}', 'bytes': (out / f_summary).stat().st_size, 'md5': md5_file(out / f_summary), 'rows': '',
         'content': 'Versions, MD5/SHA-256 of both downloads, rows kept and dropped by reason, the IEA share and the MD5 of both GO BP tables',
         'how_made_or_checked': command, **{**common_fields, 'source': 'the tables above'}},
        {'file': 'goa_human.gaf.gz', 'bytes': gaf_n, 'md5': gaf_md5, 'rows': '', 'redistribution': 'zenodo',
         'content': f'EBI GOA human annotation file (GAF {gaf_header.get("gaf-version", "")}), all aspects; SHA-256 {gaf_sha}',
         'source': GAF_URL, 'version_or_date': f'generated {gaf_header.get("date-generated", "")}; downloaded {when}',
         'citation': CITATION, 'licence': LICENCE, 'derived_from': '',
         'how_made_or_checked': ('Downloaded as is (gzip -t passes); the current file has no versioned URL, so it is kept in '
                                 'the Zenodo archive of the release'), 'original_name': 'goa_human.gaf.gz'},
        {'file': 'go-basic.obo', 'bytes': obo_n, 'md5': obo_md5, 'rows': '', 'redistribution': 'checksum_only' if args.go_build else 'zenodo',
         'content': f'Gene Ontology go-basic, release {go_release}; SHA-256 {obo_sha}',
         'source': OBO_URL, 'version_or_date': f'{obo_header.get("data-version", "")}; downloaded {when}',
         'citation': CITATION, 'licence': LICENCE, 'derived_from': '',
         'how_made_or_checked': (f'Kept at its permanent GO release URL {GO_RELEASE_URL.format(build=args.go_build)} (same size and ETag '
                                 'as the current download)' if args.go_build else 'Downloaded from the current GO URL'),
         'original_name': 'go-basic.obo'},
    ]
    update_manifest(args.manifest, rows_manifest)
    print(json.dumps({k: summary[k] for k in ('gaf_rows', 'proteins', 'proteins_with_non_iea_term', 'rows', 'rows_direct', 'rows_propagated',
                                               'rows_iea_only_share', 'bp_rows_kept_iea_share', 'direct_pairs_iea_only_share', 'terms_with_proteins',
                                               'dropped_term_ids')}, indent=1))


if __name__ == '__main__':
    main()

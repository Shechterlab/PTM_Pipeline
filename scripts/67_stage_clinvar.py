#!/usr/bin/env python3
"""Stage the ClinVar variants that change a candidate residue of the site sets (arginine, lysine, serine, threonine,
tyrosine) on the canonical sequence, with their classification, and each protein's background counts.

ClinVar (NCBI) holds what submitters report about human variants: germline classifications (pathogenic to benign), and
since 2024 somatic clinical impact and oncogenicity. NCBI "places no restrictions on the use or distribution of the data"
(https://www.ncbi.nlm.nih.gov/home/about/policies/), though submitters may claim rights in what they submitted.

Sources
  hgvs4variation.txt.gz (ClinVar FTP, tab_delimited): every HGVS expression of each variant, among them its protein changes
    on RefSeq proteins (NP_).
  variant_summary.txt.gz (ClinVar FTP): each variant's aggregate classifications, review status, submitters, conditions
    and origin (the GRCh38 rows).
  UniProt reviewed human entries: their RefSeq proteins and displayed isoforms (--uniprot-refseq, --uniprot-isoforms), and
    the reference proteome of the release tables (UniProt 2026_03, inputs/reference/).

Rules
  A variant counts at a residue when one of its protein changes is a single-residue substitution or a stop (p.Arg123His,
  p.Arg123Ter, p.Arg123*) on a RefSeq protein UniProt assigns to a reviewed entry's displayed isoform (a RefSeq protein
  without an isoform label counts for an entry without isoforms), and the 2026_03 canonical sequence has that residue at
  that position. One row per (entry, position, new residue, ClinVar variation).
  Review stars follow ClinVar (https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/): practice guideline 4, reviewed
  by expert panel 3, criteria provided by multiple submitters without conflicts 2, criteria provided by a single
  submitter or with conflicting classifications 1, otherwise 0.

Outputs (in --outdir; gzip without a timestamp)
  clinvar_changes_<X>.tsv.gz             per residue X (R, K, S, T, Y): one row per entry, position, new residue and
                                         variation, with the RefSeq protein, the variant type, the germline classification,
                                         its simple form (1 = pathogenic or likely pathogenic in some submission), review
                                         status and stars, submitters, origin, conditions, last evaluated, the somatic
                                         clinical impact and oncogenicity classifications, and the dbSNP ID (the
                                         variation's ClinVar page, by its ID, lists its RCV records), and gene: the HGNC
                                         approved symbol of the gene ClinVar gives the expression (NCBI Gene ID)
  clinvar_protein_background.tsv.gz      per reference entry and residue class (R, K, S, T, Y, all): the residues of that
                                         class, and the ClinVar variants on them, those classified pathogenic or likely
                                         pathogenic (germline), and the positions they change
  clinvar_summary.json                   inputs (MD5, ClinVar file dates), rules and counts
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import io
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from common import md5_file, parse_fasta, refresh_manifest_row, update_manifest, upsert_tsv_rows
from ensembl_genes import HGNC, HGNC_URL

AA = 'ACDEFGHIKLMNPQRSTVWY'
CLASSES = ['R', 'K', 'S', 'T', 'Y']
THREE = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
         'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'}
PROTEIN = re.compile(r'^(NP_\d+)(?:\.\d+)?:p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|\*)\)?$')
STARS = {'practice guideline': 4, 'reviewed by expert panel': 3, 'criteria provided, multiple submitters, no conflicts': 2,
         'criteria provided, single submitter': 1, 'criteria provided, conflicting classifications': 1}
COLUMNS = ['canonical_UniProtAC', 'position', 'alt', 'variation_id', 'refseq_protein', 'variant_type', 'germline_classification',
           'clinsig_simple', 'review_status', 'review_stars', 'number_submitters', 'origin', 'conditions', 'last_evaluated',
           'somatic_clinical_impact', 'oncogenicity', 'dbsnp', 'gene']     # the variation's page lists its RCV records
BACKGROUND_COLUMNS = ['canonical_UniProtAC', 'residue', 'residues', 'clinvar_variants', 'clinvar_pathogenic', 'clinvar_positions']
PATHOGENIC = re.compile(r'^(Pathogenic|Likely pathogenic|Pathogenic/Likely pathogenic)(;|$|,)', re.I)
URLS = {'hgvs4variation': 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/hgvs4variation.txt.gz',
        'variant_summary': 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz',
        'uniprot_refseq': ('https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:9606'
                           '&fields=accession,id,xref_refseq&format=tsv&compressed=true')}


def write_gz(df: pd.DataFrame, path: Path) -> None:
    buf = io.StringIO()
    df.to_csv(buf, sep='\t', index=False, lineterminator='\n')
    with open(path, 'wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0, compresslevel=9) as gz:
        gz.write(buf.getvalue().encode('utf-8'))


def canonical_refseq(refseq_tsv: str, isoforms_tsv: str) -> tuple[dict[str, set[str]], Counter]:
    """RefSeq protein (no version) -> the reviewed entries whose displayed isoform it is (as scripts/44b reads Ensembl)."""
    displayed, counts = {}, Counter()
    with gzip.open(isoforms_tsv, 'rt', encoding='utf-8') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            txt = r['Alternative products (isoforms)']
            if 'IsoId=' in txt:
                displayed[r['Entry']] = {x.strip() for d in re.findall(r'IsoId=([^;]+); Sequence=Displayed', txt) for x in d.split(',')}
    of = defaultdict(set)
    with gzip.open(refseq_tsv, 'rt', encoding='utf-8') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            acc = r['Entry']
            for item in filter(None, (x.strip() for x in r['RefSeq'].split(';'))):
                m = re.match(r'^(NP_\d+)(?:\.\d+)?(?: \[([A-Z0-9-]+)\])?$', item)
                if not m:
                    counts['refseq_xref_not_np'] += 1
                    continue
                np_, iso = m.groups()
                if (iso and iso in displayed.get(acc, set())) or (not iso and acc not in displayed):
                    of[np_].add(acc)
                    counts['canonical_refseq_proteins'] += 1
                else:
                    counts['other_isoform_refseq_proteins'] += 1
    return of, counts


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--hgvs4variation', required=True)
    ap.add_argument('--variant-summary', required=True)
    ap.add_argument('--uniprot-refseq', required=True, help='UniProt reviewed human: Entry, Entry Name, RefSeq (TSV, gzip)')
    ap.add_argument('--uniprot-isoforms', required=True, help='UniProt reviewed human: Entry, Alternative products (TSV, gzip)')
    ap.add_argument('--uniprot-release', required=True)
    ap.add_argument('--clinvar-date', required=True, help='the date of the ClinVar files (their Last-Modified on the FTP)')
    ap.add_argument('--hgnc', required=True, help=f'HGNC complete set ({HGNC_URL}): current approved symbols')
    ap.add_argument('--fasta', default='inputs/reference/uniprot_human_reviewed_2026_03.fasta')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--manifest', default='', help='inputs/MANIFEST.tsv: write the rows of the staged files and of the downloads')
    ap.add_argument('--downloaded', default=dt.date.today().isoformat())
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    canonical = parse_fasta(args.fasta)
    canon, counts = canonical_refseq(args.uniprot_refseq, args.uniprot_isoforms)

    # protein changes on canonical RefSeq proteins, checked on the 2026_03 sequence
    changes = {}
    with gzip.open(args.hgvs4variation, 'rt', encoding='utf-8') as fh:
        head = None
        for line in fh:
            if line.startswith('#Symbol'):
                head = line[1:].rstrip('\n').split('\t')
                ix = {c: k for k, c in enumerate(head)}
                continue
            if head is None or line.startswith('#'):
                continue
            r = line.rstrip('\n').split('\t')
            m = PROTEIN.match(r[ix['ProteinExpression']])
            if not m:
                continue
            counts['protein_substitutions_or_stops'] += 1
            np_, ref3, pos, alt3 = m.group(1), m.group(2), int(m.group(3)), m.group(4)
            if ref3 not in THREE or (alt3 not in THREE and alt3 not in ('Ter', '*')):
                continue
            ref, alt = THREE[ref3], ('*' if alt3 in ('Ter', '*') else THREE[alt3])
            if alt == ref or np_ not in canon:
                counts['on_other_refseq_proteins' if np_ not in canon else 'synonymous'] += 1
                continue
            for acc in sorted(canon[np_]):
                seq = canonical.get(acc, '')
                if not (1 <= pos <= len(seq) and seq[pos - 1] == ref):
                    counts['residue_not_on_canonical'] += 1
                    continue
                key = (acc, pos, alt, r[ix['VariationID']])
                if key not in changes:
                    changes[key] = (ref, r[ix['AlleleID']], r[ix['ProteinExpression']].split(':')[0], r[ix['Symbol']], r[ix['GeneID']])
    counts['changes_placed_all_residues'] = len(changes)

    # each variation's classifications (GRCh38 rows of variant_summary)
    wanted = {k[3] for k in changes}
    info = {}
    with gzip.open(args.variant_summary, 'rt', encoding='utf-8') as fh:
        reader = csv.reader(fh, delimiter='\t')
        h = next(reader)
        ix = {c.lstrip('#'): k for k, c in enumerate(h)}
        for r in reader:
            vid = r[ix['VariationID']]
            if vid not in wanted or (vid in info and r[ix['Assembly']] != 'GRCh38'):
                continue
            rs = r[ix['RS# (dbSNP)']]
            info[vid] = {'variant_type': r[ix['Type']], 'germline_classification': r[ix['ClinicalSignificance']],
                         'clinsig_simple': r[ix['ClinSigSimple']], 'review_status': r[ix['ReviewStatus']],
                         'review_stars': STARS.get(r[ix['ReviewStatus']], 0), 'number_submitters': r[ix['NumberSubmitters']],
                         'origin': r[ix['OriginSimple']], 'conditions': r[ix['PhenotypeList']], 'last_evaluated': r[ix['LastEvaluated']],
                         'somatic_clinical_impact': r[ix['SomaticClinicalImpact']], 'oncogenicity': r[ix['Oncogenicity']],
                         'dbsnp': '' if rs in ('-1', '') else 'rs' + rs, 'rcv': r[ix['RCVaccession']], 'summary_gene': r[ix['GeneSymbol']]}
    counts['variations_without_summary'] = len(wanted - set(info))

    rows = []
    hgnc = HGNC(args.hgnc)
    for (acc, pos, alt, vid), (ref, allele, np_, symbol, gene_id) in changes.items():
        if vid not in info:
            continue
        # the HGNC approved symbol of the NCBI gene ClinVar gives the RefSeq protein's expression, else of its symbol, else of
        # the variation's gene in variant_summary
        gene = hgnc.symbol(entrez=gene_id, name=symbol if symbol not in ('', '-') else info[vid]['summary_gene'])
        rows.append({'canonical_UniProtAC': acc, 'position': pos, 'alt': alt, 'variation_id': vid, 'allele_id': allele,
                     'refseq_protein': np_, 'ref': ref, 'gene': gene} | info[vid])
    t = pd.DataFrame(rows)
    t['variation_id_num'] = t['variation_id'].astype(int)
    t = t.sort_values(['canonical_UniProtAC', 'position', 'alt', 'variation_id_num'], kind='stable')
    written = {}
    for res in CLASSES:
        sub = t[t['ref'] == res][COLUMNS]
        write_gz(sub, out / f'clinvar_changes_{res}.tsv.gz')
        written[res] = {'variants': len(sub), 'positions': int(sub[['canonical_UniProtAC', 'position']].drop_duplicates().shape[0]),
                        'pathogenic_or_likely': int(sub['germline_classification'].str.match(PATHOGENIC).sum()),
                        'md5': md5_file(out / f'clinvar_changes_{res}.tsv.gz')}

    # background, over every placed variant with a summary
    agg = defaultdict(lambda: [0, 0, set()])
    for r in t.itertuples(index=False):
        patho = bool(PATHOGENIC.match(r.germline_classification))
        for cls in ([r.ref] if r.ref in CLASSES else []) + ['all']:
            a = agg[(r.canonical_UniProtAC, cls)]
            a[0] += 1
            a[1] += patho
            a[2].add(r.position)
    bg = []
    for acc in sorted(canonical):
        seq = canonical[acc]
        for cls in CLASSES + ['all']:
            a = agg.get((acc, cls), [0, 0, set()])
            bg.append([acc, cls, len(seq) if cls == 'all' else seq.count(cls), a[0], a[1], len(a[2])])
    write_gz(pd.DataFrame(bg, columns=BACKGROUND_COLUMNS), out / 'clinvar_protein_background.tsv.gz')

    summary = {
        'generated': dt.date.today().isoformat(),
        'inputs': {'hgvs4variation': {'file': Path(args.hgvs4variation).name, 'url': URLS['hgvs4variation'], 'md5': md5_file(args.hgvs4variation)},
                   'variant_summary': {'file': Path(args.variant_summary).name, 'url': URLS['variant_summary'], 'md5': md5_file(args.variant_summary)},
                   'clinvar_files_dated': args.clinvar_date,
                   'uniprot': {'release': args.uniprot_release,
                               'refseq': {'file': Path(args.uniprot_refseq).name, 'url': URLS['uniprot_refseq'], 'md5': md5_file(args.uniprot_refseq)},
                               'isoforms': {'file': Path(args.uniprot_isoforms).name, 'md5': md5_file(args.uniprot_isoforms)}},
                   'hgnc': {'file': Path(args.hgnc).name, 'url': HGNC_URL, 'md5': md5_file(args.hgnc), 'license': 'CC0'},
                   'reference': {'fasta': args.fasta, 'md5': md5_file(args.fasta)}},
        'license': 'NCBI places no restrictions on the use or distribution of ClinVar data; submitters may claim rights in what they submitted',
        'rule': 'a single-residue substitution or stop on a RefSeq protein UniProt assigns to a reviewed entry\'s displayed isoform, '
                'kept where the 2026_03 canonical sequence has the residue at the position; classifications from the GRCh38 rows of '
                'variant_summary; one row per (entry, position, new residue, variation)',
        'residues': written,
        'background_md5': md5_file(out / 'clinvar_protein_background.tsv.gz'),
        'counts': dict(sorted(counts.items())),
    }
    (out / 'clinvar_summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    if args.manifest:
        write_manifest(Path(args.manifest), out, args)
    print(json.dumps({k: {x: v[x] for x in v if x != 'md5'} for k, v in written.items()}, indent=1))


def gz_rows(path: Path) -> int:
    with gzip.open(path, 'rt', encoding='utf-8') as fh:
        return sum(1 for _ in fh) - 1


def sha256(p: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(p, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def write_manifest(manifest: Path, out: Path, args) -> None:
    folder = out.name
    src = f'{URLS["hgvs4variation"]}; {URLS["variant_summary"]}'
    ver = f'ClinVar files dated {args.clinvar_date}, downloaded {args.downloaded}; UniProt {args.uniprot_release} RefSeq cross-references'
    cit = 'Landrum et al., Nucleic Acids Res 46:D1062 (2018), doi:10.1093/nar/gkx1153 (ClinVar)'
    lic = 'NCBI places no restrictions on the use or distribution of ClinVar data (https://www.ncbi.nlm.nih.gov/home/about/policies/); submitters may claim rights in their data'
    how = ('scripts/67_stage_clinvar.py; checked against reference/uniprot_human_reviewed_2026_03.fasta (a variant counts only where '
           'the canonical sequence has the changed residue at its position)')
    names = {'R': 'arginine', 'K': 'lysine', 'S': 'serine', 'T': 'threonine', 'Y': 'tyrosine'}
    rows = []
    for res in CLASSES:
        rows.append({'file': f'{folder}/clinvar_changes_{res}.tsv.gz', 'content': f'ClinVar variants that change a {names[res]} ({res}) of a '
                     'reviewed canonical protein (a substitution or a stop on its RefSeq protein): one row per entry, position, new residue '
                     'and variation, with the germline classification, review status and stars, submitters, origin, conditions, the '
                     'somatic clinical impact and oncogenicity, and the dbSNP ID (gzip)', 'source': src, 'version_or_date': ver,
                     'citation': cit, 'licence': lic, 'derived_from': 'hgvs4variation.txt.gz; variant_summary.txt.gz'})
    rows.append({'file': f'{folder}/clinvar_protein_background.tsv.gz', 'content': 'Per reviewed canonical entry and residue class (R, K, S, '
                 'T, Y, all): the residues of that class, the ClinVar variants on them, those classified pathogenic or likely pathogenic, '
                 'and the positions changed (gzip)', 'source': src, 'version_or_date': ver, 'citation': cit, 'licence': lic,
                 'derived_from': 'hgvs4variation.txt.gz; variant_summary.txt.gz; reference/uniprot_human_reviewed_2026_03.fasta'})
    for r in rows:
        path = out / r['file'].split('/', 1)[1]
        r.update({'bytes': path.stat().st_size, 'md5': md5_file(path), 'rows': gz_rows(path), 'redistribution': 'git', 'how_made_or_checked': how})
    s = out / 'clinvar_summary.json'
    rows.append({'file': f'{folder}/clinvar_summary.json', 'bytes': s.stat().st_size, 'md5': md5_file(s), 'rows': '', 'redistribution': 'git',
                 'content': 'Inputs (MD5, ClinVar file date), rule and counts of the ClinVar staging', 'source': 'the downloads above',
                 'version_or_date': ver, 'citation': cit, 'licence': lic, 'derived_from': '', 'how_made_or_checked': 'scripts/67_stage_clinvar.py'})
    update_manifest(manifest, rows)
    raw = []
    for path, url, v in ((Path(args.hgvs4variation), URLS['hgvs4variation'], f'ClinVar, dated {args.clinvar_date}'),
                         (Path(args.variant_summary), URLS['variant_summary'], f'ClinVar, dated {args.clinvar_date}'),
                         (Path(args.uniprot_refseq), URLS['uniprot_refseq'], f'UniProt {args.uniprot_release}, RefSeq cross-references'),
                         (Path(args.hgnc), HGNC_URL, 'HGNC complete set, file of 2026-11-13 (CC0)')):
        raw.append({'file': path.name, 'url': url, 'bytes': path.stat().st_size, 'md5': md5_file(path), 'sha256': sha256(path),
                    'downloaded': args.downloaded, 'version': v, 'kept': 'release_asset'})
    upsert_tsv_rows(out / 'raw_files.tsv', raw)
    refresh_manifest_row(manifest, out / 'raw_files.tsv', f'{folder}/raw_files.tsv')


if __name__ == '__main__':
    main()

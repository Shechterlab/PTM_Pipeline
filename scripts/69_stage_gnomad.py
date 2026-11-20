#!/usr/bin/env python3
"""Stage the gnomAD v4.1 exome variants that change a candidate residue of the site sets (arginine, lysine, serine,
threonine, tyrosine) on the canonical sequence, with their allele counts, and each protein's background counts.

gnomAD v4.1 (Chen et al., Nature 625:92, 2024; the v4.1 release of April 2024) aggregates 730,947 exomes. "The primary data
from the gnomAD exomes and genomes are available free of restrictions under the Creative Commons Zero Public Domain
Dedication" (https://github.com/broadinstitute/gnomad-browser/blob/main/browser/about/policies/terms.md). Its SpliceAI
annotations (CC BY-NC 4.0) and CADD scores (non-commercial) are not read. The files are Google's public copy
(gs://gcp-public-data--gnomad/release/4.1/vcf/exomes/, read anonymously over HTTPS); no project or credential is needed.

    python scripts/69_stage_gnomad.py fetch --chrom chr21 --outdir WORK          # one chromosome, streamed
    python scripts/69_stage_gnomad.py stage --partials WORK --outdir inputs/layers --manifest inputs/MANIFEST.tsv ...

fetch streams a chromosome's sites VCF (curl | gzip -dc | grep), never storing it, and keeps the variants that pass every
gnomAD filter (FILTER PASS) and that VEP calls missense or stop gained on some transcript: for each, its allele counts and
the VEP entries of those consequences (transcript, protein position, residues). It writes WORK/<chrom>.tsv.gz and, when
the stream ended cleanly, WORK/<chrom>.done with the bytes read; a broken stream is retried.

stage reads the partials and keeps, per variant, the change on a transcript UniProt assigns to a reviewed entry's displayed
isoform (as scripts/44b reads Ensembl), where the 2026_03 canonical sequence has the residue VEP names at that position.

Outputs (in --outdir; gzip without a timestamp)
  gnomad_v4_1_changes_<X>.tsv.gz         per residue X (R, K, S, T, Y): one row per entry, position, new residue and
                                         variant (GRCh38 chrom-pos-ref-alt), with dbSNP ID, allele count and number,
                                         allele frequency, homozygotes, the group with the highest frequency and its
                                         frequency, and the highest filtering allele frequency (95%), and gene: the HGNC
                                         approved symbol of the Ensembl 105 gene of the transcript VEP read the change on
                                         (';' between overlapping genes)
  gnomad_protein_background.tsv.gz       per reference entry and residue class (R, K, S, T, Y, all): residues, variants,
                                         positions changed, and the summed allele count
  gnomad_summary.json                    inputs (bytes and MD5 read per chromosome, file dates), rules and counts
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import importlib.util
import io
import json
import re
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common import md5_file, parse_fasta, refresh_manifest_row, update_manifest, upsert_tsv_rows  # noqa: E402
from ensembl_genes import GTF_105, HGNC, HGNC_URL, transcript_genes  # noqa: E402

BASE = 'https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.{chrom}.vcf.bgz'
CHROMS = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
CLASSES = ['R', 'K', 'S', 'T', 'Y']
AA = set('ACDEFGHIKLMNPQRSTVWY')
KEYS = ('AC', 'AN', 'AF', 'nhomalt', 'grpmax', 'AF_grpmax', 'fafmax_faf95_max')
PARTIAL = ['chrom', 'pos', 'ref', 'alt', 'rsid'] + list(KEYS) + ['transcript', 'protein_position', 'ref_aa', 'alt_aa', 'consequence']
COLUMNS = ['canonical_UniProtAC', 'position', 'alt', 'variant_id', 'dbsnp', 'ac', 'an', 'af', 'nhomalt', 'grpmax', 'af_grpmax',
           'faf95_max', 'consequence', 'gene']
BACKGROUND_COLUMNS = ['canonical_UniProtAC', 'residue', 'residues', 'gnomad_variants', 'gnomad_positions', 'gnomad_allele_count']
LICENSE = ('CC0 1.0 (gnomAD primary data, https://github.com/broadinstitute/gnomad-browser/blob/main/browser/about/policies/terms.md); '
           'SpliceAI (CC BY-NC 4.0) and CADD (non-commercial) annotations not used')
CITATION = 'Chen et al., Nature 625:92 (2024), doi:10.1038/s41586-023-06045-0 (gnomAD v4)'


def info_value(info: str, key: str) -> str:
    m = re.search(rf'(?:^|;){key}=([^;]*)', info)
    return m.group(1) if m else ''


def parse_line(line: str, csq_ix: dict, keep_tx: set) -> list[list[str]]:
    f = line.rstrip('\n').split('\t', 8)
    if len(f) < 8 or f[6] != 'PASS':
        return []
    info = f[7]
    vep = info_value(info, 'vep')
    vals = [info_value(info, k) for k in KEYS]
    rows = []
    for entry in vep.split(','):
        e = entry.split('|')
        if len(e) <= max(csq_ix.values()):
            continue
        if e[csq_ix['Feature']].split('.')[0] not in keep_tx:   # transcripts of reviewed entries' displayed isoforms only
            continue
        cons = e[csq_ix['Consequence']]
        if 'missense_variant' not in cons and 'stop_gained' not in cons:
            continue
        aa = e[csq_ix['Amino_acids']].split('/')
        pos = e[csq_ix['Protein_position']]
        if len(aa) != 2 or len(aa[0]) != 1 or len(aa[1]) != 1 or not pos.isdigit():
            continue
        rows.append([f[0], f[1], f[3], f[4], f[2]] + vals + [e[csq_ix['Feature']].split('.')[0], pos, aa[0], aa[1], cons])
    return rows


def write_gz(df: pd.DataFrame, path: Path) -> None:
    buf = io.StringIO()
    df.to_csv(buf, sep='\t', index=False, lineterminator='\n')
    with open(path, 'wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0, compresslevel=9) as gz:
        gz.write(buf.getvalue().encode('utf-8'))


def fetch(chrom: str, outdir: Path, keep_tx: set, retries: int = 3) -> None:
    url = BASE.format(chrom=chrom)
    outdir.mkdir(parents=True, exist_ok=True)
    for attempt in range(1, retries + 1):
        md5 = hashlib.md5()
        counter = {'bytes': 0}
        curl = subprocess.Popen(['curl', '-sS', '--fail', url], stdout=subprocess.PIPE)
        gz = subprocess.Popen(['gzip', '-dc'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        grep = subprocess.Popen(['grep', '-F', '-e', 'missense_variant', '-e', 'stop_gained', '-e', '##INFO=<ID=vep,'],
                                stdin=gz.stdout, stdout=subprocess.PIPE, text=True)
        gz.stdout.close()

        def pump():
            for block in iter(lambda: curl.stdout.read(1 << 20), b''):
                md5.update(block)
                counter['bytes'] += len(block)
                gz.stdin.write(block)
            gz.stdin.close()

        import threading
        t = threading.Thread(target=pump, daemon=True)
        t.start()
        csq_ix, n_rows, lines = None, 0, 0
        tmp = outdir / f'{chrom}.tsv.gz.part'
        with open(tmp, 'wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0, compresslevel=6) as gzout:
            gzout.write(('\t'.join(PARTIAL) + '\n').encode())
            for line in grep.stdout:
                if line.startswith('##INFO=<ID=vep,'):
                    fmt = line.split('Format: ', 1)[1].split('"', 1)[0].split('|')
                    csq_ix = {k: i for i, k in enumerate(fmt)}
                    continue
                if csq_ix is None or line.startswith('#'):
                    continue
                lines += 1
                for r in parse_line(line, csq_ix, keep_tx):
                    gzout.write(('\t'.join(r) + '\n').encode())
                    n_rows += 1
        t.join()
        codes = (curl.wait(), gz.wait(), grep.wait())
        if codes[0] == 0 and codes[1] == 0 and codes[2] in (0, 1) and csq_ix is not None:
            tmp.replace(outdir / f'{chrom}.tsv.gz')
            (outdir / f'{chrom}.done').write_text(json.dumps({'url': url, 'bytes': counter['bytes'], 'md5': md5.hexdigest(),
                                                              'lines_with_missense_or_stop': lines, 'entries': n_rows,
                                                              'fetched': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')}) + '\n')
            print(chrom, counter['bytes'], 'bytes;', lines, 'variant lines;', n_rows, 'entries', flush=True)
            return
        print(f'{chrom}: attempt {attempt} failed (curl, gzip, grep exit {codes}); retrying', file=sys.stderr, flush=True)
        time.sleep(30 * attempt)
    raise SystemExit(f'{chrom}: failed after {retries} attempts')


def canonical_transcripts(ensembl_tsv: str, isoforms_tsv: str) -> dict:
    spec = importlib.util.spec_from_file_location('extract_tcga_mc3', HERE / '44b_extract_tcga_mc3.py')
    s44b = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(s44b)
    return s44b.canonical_transcripts(ensembl_tsv, isoforms_tsv)[0]


def stage(args) -> None:
    """Per chromosome (a variant lies on one): the changes on canonical transcripts, checked on the 2026_03 sequence; the
    residue tables are written per chromosome to a work folder and joined, sorted, at the end."""
    partials, out = Path(args.partials), Path(args.outdir)
    missing = [c for c in CHROMS if not (partials / f'{c}.done').exists()]
    if missing:
        raise SystemExit(f'chromosomes not fetched: {missing}')
    canonical = parse_fasta(args.fasta)
    canon = canonical_transcripts(args.uniprot_ensembl, args.uniprot_isoforms)
    tx_acc = pd.DataFrame([(tx, a) for tx, accs in canon.items() for a in sorted(accs)], columns=['transcript', 'canonical_UniProtAC'])
    tx_gene, hgnc = transcript_genes(args.ensembl_gtf), HGNC(args.hgnc)
    gene_of = {tx: hgnc.symbol(*eg) or eg[0] for tx, eg in tx_gene.items()}     # the Ensembl gene for a gene without an HGNC symbol
    work = partials / 'staged'
    work.mkdir(exist_ok=True)
    counts = Counter()
    bg_counts, bg_positions = [], []
    for c in CHROMS:
        p = pd.read_csv(partials / f'{c}.tsv.gz', sep='\t', dtype=str, keep_default_na=False)
        counts['vep_entries_on_canonical_transcripts'] += len(p)
        p = p.merge(tx_acc, on='transcript')
        p = p[p['ref_aa'].isin(AA) & (p['alt_aa'].isin(AA) | (p['alt_aa'] == '*')) & (p['alt_aa'] != p['ref_aa'])]
        p['position'] = p['protein_position'].astype(int)
        seq_res = [canonical.get(a, '')[q - 1] if 1 <= q <= len(canonical.get(a, '')) else '' for a, q in zip(p['canonical_UniProtAC'], p['position'])]
        ok = p['ref_aa'].to_numpy() == pd.Series(seq_res, index=p.index).to_numpy()
        counts['residue_not_on_canonical'] += int((~ok).sum())
        p = p[ok]
        p['variant_id'] = p['chrom'].str.removeprefix('chr') + '-' + p['pos'] + '-' + p['ref'] + '-' + p['alt']
        p = p.drop(columns=['ref', 'alt'])          # the alleles live in variant_id; alt is the new residue below
        p = p.rename(columns={'alt_aa': 'alt_res', 'rsid': 'dbsnp', 'AC': 'ac', 'AN': 'an', 'AF': 'af', 'nhomalt': 'nhomalt',
                              'grpmax': 'grpmax', 'AF_grpmax': 'af_grpmax', 'fafmax_faf95_max': 'faf95_max'})
        # the gene: Ensembl 105 name of the transcript VEP read the change on; the names of every such transcript of the same
        # change, separated by ';' (overlapping genes)
        p['gene'] = p['transcript'].map(gene_of).fillna('')
        key = ['canonical_UniProtAC', 'position', 'alt_res', 'variant_id']
        genes = p.groupby(key)['gene'].agg(lambda s: ';'.join(sorted({g for g in s if g}))).rename('genes').reset_index()
        p = p.drop_duplicates(key).drop(columns='gene').merge(genes, on=key, how='left').rename(columns={'genes': 'gene'})
        counts['changes_without_ensembl_105_gene'] += int((p['gene'] == '').sum())
        p['dbsnp'] = p['dbsnp'].replace('.', '')
        p['consequence'] = ['stop_gained' if a == '*' else 'missense_variant' for a in p['alt_res']]
        counts['changes_placed_all_residues'] += len(p)
        cls = p['ref_aa'].where(p['ref_aa'].isin(CLASSES), '')
        for frame in (p.assign(cls=cls)[cls != ''], p.assign(cls='all')):
            g = frame.groupby(['canonical_UniProtAC', 'cls']).agg(n=('variant_id', 'size'), ac=('ac', lambda s: pd.to_numeric(s, errors='coerce').fillna(0).astype(int).sum()))
            bg_counts.append(g)
            bg_positions.append(frame[['canonical_UniProtAC', 'cls', 'position']].drop_duplicates())
        for res in CLASSES:
            sub = p[p['ref_aa'] == res].rename(columns={'alt_res': 'alt'})[COLUMNS]
            assert sub['alt'].isin(AA | {'*'}).all() and not sub.columns.duplicated().any()
            sub.to_csv(work / f'{res}_{c}.tsv', sep='\t', index=False, lineterminator='\n')
        print(c, len(p), 'changes', flush=True)
    written = {}
    for res in CLASSES:
        t = pd.concat([pd.read_csv(work / f'{res}_{c}.tsv', sep='\t', dtype=str, keep_default_na=False) for c in CHROMS], ignore_index=True)
        t['_pos'] = t['position'].astype(int)
        t = t.sort_values(['canonical_UniProtAC', '_pos', 'alt', 'variant_id'], kind='stable').drop(columns='_pos')
        write_gz(t, out / f'gnomad_v4_1_changes_{res}.tsv.gz')
        written[res] = {'variants': len(t), 'positions': int(t[['canonical_UniProtAC', 'position']].drop_duplicates().shape[0]),
                        'bytes': (out / f'gnomad_v4_1_changes_{res}.tsv.gz').stat().st_size,
                        'md5': md5_file(out / f'gnomad_v4_1_changes_{res}.tsv.gz')}
        del t
    g = pd.concat(bg_counts).groupby(level=[0, 1]).sum()
    pos = pd.concat(bg_positions).drop_duplicates().groupby(['canonical_UniProtAC', 'cls']).size()
    bg = []
    for acc in sorted(canonical):
        seq = canonical[acc]
        for cls in CLASSES + ['all']:
            n, ac = (int(g.at[(acc, cls), 'n']), int(g.at[(acc, cls), 'ac'])) if (acc, cls) in g.index else (0, 0)
            bg.append([acc, cls, len(seq) if cls == 'all' else seq.count(cls), n, int(pos.get((acc, cls), 0)), ac])
    write_gz(pd.DataFrame(bg, columns=BACKGROUND_COLUMNS), out / 'gnomad_protein_background.tsv.gz')
    fetched = {c: json.loads((partials / f'{c}.done').read_text()) for c in CHROMS}
    summary = {'generated': dt.date.today().isoformat(), 'release': 'gnomAD v4.1 exomes (730,947 exomes; files dated 2024-04-17)',
               'license': LICENSE, 'citation': CITATION,
               'inputs': {'vcfs': fetched, 'uniprot': {'release': args.uniprot_release, 'ensembl_md5': md5_file(args.uniprot_ensembl),
                                                       'isoforms_md5': md5_file(args.uniprot_isoforms)},
                          'ensembl_gtf': {'file': Path(args.ensembl_gtf).name, 'url': GTF_105, 'md5': md5_file(args.ensembl_gtf)},
                          'hgnc': {'file': Path(args.hgnc).name, 'url': HGNC_URL, 'md5': md5_file(args.hgnc), 'license': 'CC0'},
                          'reference': {'fasta': args.fasta, 'md5': md5_file(args.fasta)}},
               'rule': 'FILTER PASS; a missense or stop-gained VEP consequence on a transcript UniProt assigns to a reviewed entry\'s '
                       'displayed isoform, kept where the 2026_03 canonical sequence has the residue VEP names at the protein '
                       'position; one row per (entry, position, new residue, variant)',
               'residues': written, 'background_md5': md5_file(out / 'gnomad_protein_background.tsv.gz'), 'counts': dict(sorted(counts.items()))}
    (out / 'gnomad_summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    if args.manifest:
        write_manifest(Path(args.manifest), out, args, fetched)
    print(json.dumps({k: {x: v[x] for x in v if x != 'md5'} for k, v in written.items()}, indent=1))


def write_manifest(manifest: Path, out: Path, args, fetched: dict) -> None:
    folder = out.name
    src = 'https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/ (gnomad.exomes.v4.1.sites.<chrom>.vcf.bgz, streamed)'
    ver = f'gnomAD v4.1 exomes (files dated 2024-04-17), read {args.downloaded}; UniProt {args.uniprot_release} Ensembl cross-references'
    how = ('scripts/69_stage_gnomad.py; checked against reference/uniprot_human_reviewed_2026_03.fasta (a variant counts only where '
           'the canonical sequence has the residue VEP names at its position)')
    names = {'R': 'arginine', 'K': 'lysine', 'S': 'serine', 'T': 'threonine', 'Y': 'tyrosine'}
    rows = []
    for res in CLASSES:
        rows.append({'file': f'{folder}/gnomad_v4_1_changes_{res}.tsv.gz', 'content': f'gnomAD v4.1 exome variants (FILTER PASS) that change a '
                     f'{names[res]} ({res}) of a reviewed canonical protein (missense or stop gained): one row per entry, position, new '
                     'residue and variant, with dbSNP ID, allele count and number, frequency, homozygotes, the group with the highest '
                     'frequency and the highest filtering allele frequency (gzip)', 'source': src, 'version_or_date': ver,
                     'citation': CITATION, 'licence': LICENSE, 'derived_from': 'gnomad.exomes.v4.1.sites.<chrom>.vcf.bgz (24 files)'})
    rows.append({'file': f'{folder}/gnomad_protein_background.tsv.gz', 'content': 'Per reviewed canonical entry and residue class (R, K, S, T, '
                 'Y, all): residues, gnomAD v4.1 exome variants on them, positions changed and summed allele count (gzip)', 'source': src,
                 'version_or_date': ver, 'citation': CITATION, 'licence': LICENSE,
                 'derived_from': 'gnomad.exomes.v4.1.sites.<chrom>.vcf.bgz; reference/uniprot_human_reviewed_2026_03.fasta'})
    for r in rows:
        path = out / r['file'].split('/', 1)[1]
        with gzip.open(path, 'rt', encoding='utf-8') as fh:
            n = sum(1 for _ in fh) - 1
        r.update({'bytes': path.stat().st_size, 'md5': md5_file(path), 'rows': n, 'redistribution': 'git', 'how_made_or_checked': how})
    s = out / 'gnomad_summary.json'
    rows.append({'file': f'{folder}/gnomad_summary.json', 'bytes': s.stat().st_size, 'md5': md5_file(s), 'rows': '', 'redistribution': 'git',
                 'content': 'Inputs (bytes and MD5 read per chromosome), rule and counts of the gnomAD staging', 'source': src,
                 'version_or_date': ver, 'citation': CITATION, 'licence': LICENSE, 'derived_from': '', 'how_made_or_checked': 'scripts/69_stage_gnomad.py'})
    update_manifest(manifest, rows)
    raw = [{'file': Path(v['url']).name, 'url': v['url'], 'bytes': v['bytes'], 'md5': v['md5'], 'sha256': '',
            'downloaded': v['fetched'][:10], 'version': 'gnomAD v4.1 exomes, sites VCF (streamed, not stored)',
            'kept': 'not kept: served by Google Cloud public datasets at this URL (CC0); MD5 of the bytes read identifies it'} for v in fetched.values()]
    gtf = Path(args.ensembl_gtf)
    raw.append({'file': gtf.name, 'url': GTF_105, 'bytes': gtf.stat().st_size, 'md5': md5_file(gtf), 'sha256': '', 'downloaded': args.downloaded,
                'version': 'Ensembl release 105 GTF (Ensembl CHECKSUMS sum 40862 49638 matched); transcript gene names',
                'kept': 'not kept: served by Ensembl at this URL; MD5 identifies it'})
    hg = Path(args.hgnc)
    raw.append({'file': hg.name, 'url': HGNC_URL, 'bytes': hg.stat().st_size, 'md5': md5_file(hg), 'sha256': '', 'downloaded': args.downloaded,
                'version': 'HGNC complete set, file of 2026-11-13 (CC0)', 'kept': 'release_asset'})
    upsert_tsv_rows(out / 'raw_files.tsv', raw)
    refresh_manifest_row(manifest, out / 'raw_files.tsv', f'{folder}/raw_files.tsv')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    f = sub.add_parser('fetch')
    f.add_argument('--chrom', required=True, choices=CHROMS)
    f.add_argument('--outdir', required=True)
    f.add_argument('--uniprot-ensembl', required=True)
    f.add_argument('--uniprot-isoforms', required=True)
    s = sub.add_parser('stage')
    s.add_argument('--partials', required=True)
    s.add_argument('--uniprot-ensembl', required=True)
    s.add_argument('--uniprot-isoforms', required=True)
    s.add_argument('--uniprot-release', required=True)
    s.add_argument('--hgnc', required=True, help=f'HGNC complete set ({HGNC_URL}): current approved symbols')
    s.add_argument('--ensembl-gtf', required=True, help=f'Ensembl 105 GTF ({GTF_105}): each transcript\'s gene name')
    s.add_argument('--fasta', default='inputs/reference/uniprot_human_reviewed_2026_03.fasta')
    s.add_argument('--outdir', required=True)
    s.add_argument('--manifest', default='')
    s.add_argument('--downloaded', default=dt.date.today().isoformat())
    args = ap.parse_args()
    if args.cmd == 'fetch':
        fetch(args.chrom, Path(args.outdir), set(canonical_transcripts(args.uniprot_ensembl, args.uniprot_isoforms)))
    else:
        stage(args)


if __name__ == '__main__':
    main()

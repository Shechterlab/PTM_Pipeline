#!/usr/bin/env python3
"""Stage ENCORE eCLIP (Her et al., Cell Systems 2026; GEO SuperSeries GSE315407) beside ENCODE eCLIP (inputs/nucleic_acid/).

    python -I scripts/43e_stage_encore.py --cache <download folder> --eclip-cache <the cache scripts/43 eclip used> \
        [--out inputs/nucleic_acid] [--skip-encode-kmers]

What ENCORE is, as GEO (SubSeries GSE315347 and GSE315406 of GSE315407) and the Skipper paper (Boyle et al., Cell Genomics
3:100317, 2023) describe it: eCLIP of RNA-binding proteins in K562 and HepG2 cells, two biological replicates each with a
size-matched input, reads mapped to hg38 with STAR and analysed with Skipper v1.99.0 (GEO). Skipper, as its paper describes
it, tiles annotated transcripts (GENCODE; transcripts below 1 TPM in the cell line left out) into windows of at most 100 nt,
counts read starts in IP and input, and tests each window for enrichment of IP over input with a beta-binomial model that
learns overdispersion and GC bias (windows in GC bins); a window is reproducible when it passes FDR 0.2 in both replicates
(GEO: only uniquely mapped reads for these files). Each experiment has one *.reproducible_enriched_windows.tsv.gz: coordinates (0-based start, end
exclusive; checked below against GENCODE feature starts), strand, Skipper's feature type and gene, read sums of IP and
input, the log2 odds ratio of enrichment in each replicate (min, mean, max) and p and q values. Every file is sorted by
the mean log2 odds ratio, which is the window's score here.

Tables (inputs/nucleic_acid/, each with its row in inputs/MANIFEST.tsv):
  encore_eclip_experiments.tsv      one row per experiment and reviewed UniProt entry: GEO series, file and samples, cell line,
                                    number of windows, and the ENCODE experiment with the same Yeo-lab number (ENCODE alias
                                    gene-yeo:<number>, same cell line) when ENCODE has one
  encore_eclip_window_regions.tsv   windows per GENCODE v50 basic region of the window's middle, placed as scripts/43 places
                                    ENCODE eCLIP peaks (CDS, UTR3, UTR5, noncoding_exon, intron, intergenic)
  encore_eclip_window_features.tsv  windows per Skipper feature type, as the files give it
  encore_eclip_target_genes.tsv     the 200 genes (GENCODE v50) with most windows per experiment, as for ENCODE eCLIP
  encore_eclip_top_windows.tsv      the 200 windows with the highest mean log2 odds ratio per experiment, with the RNA sequence
                                    and the number of ENCORE experiments in which the same window is reproducible: the files' top
                                    windows include windows with no input reads that many unrelated RBPs share, which the
                                    ranking keeps (it is the files' own) and this count flags
  encore_eclip_5mer_enrichment.tsv  5-mer enrichment of each experiment's top windows (below)
  encode_eclip_5mer_enrichment.tsv  the same, by the same method, for the ENCODE eCLIP reproducible peaks scripts/43 staged
  eclip_5mer_backgrounds.tsv        per experiment of either source: sequences and 5-mer positions behind those numbers, and
                                    the top 5-mers
  encore_raw_files.tsv              every download, with URL, size, MD5 and date

RBP to reviewed UniProt entry (release 2026_03): the target of the ENCODE experiment with the same Yeo-lab number, by its
UniProt cross-references (as scripts/43 eclip); else the primary gene symbol (when a symbol has several entries, the one
named <SYMBOL>_HUMAN: RBM10, not the MINAS-60 microprotein of the same gene); else the one reviewed human entry UniProt's
gene_exact search finds for the name.

RNA sequence: from the GRCh38 primary assembly of GENCODE release 50 (EBI; MD5 checked against GENCODE's MD5SUMS), read
0-based end-exclusive, reverse-complemented for the minus strand, written 5' to 3' in U alphabet (N where the assembly has
N). Checked: the GT..AG ends of a sample of GENCODE v50 introns read through the same code (share in the MANIFEST row).

5-mer enrichment, per experiment (the same code for both sources):
  foreground   the experiment's top K sequences (K = 1,000 by default; all when fewer): ENCORE windows ranked by the mean
               log2 odds ratio, ENCODE peaks by signalValue (column 7, log2 fold enrichment), then -log10 p (column 8).
               Every overlapping 5-mer inside a sequence is counted (none containing N); n_k of 5-mer k, N in all.
  flank        for each foreground sequence of length L, the L nt immediately 5' and the L nt immediately 3' of it on the
               same strand, less every nucleotide inside any reproducible window (or peak) of the same experiment; 5-mers
               counted inside what remains.
  other        the foreground sequences of every other experiment of the same source whose RBP differs (no UniProt
               entry in common), by GENCODE v50 region of the sequence's middle, reweighted to the experiment's own mix:
               q_k = sum over regions r of w_r (o_rk + 0.5) / (O_r + 512), w_r the share of the foreground's 5-mer
               positions in region r, o_rk and O_r the background counts in r.
  R            (n_k + 0.5) / (N + 512) over the background frequency ((m_k + 0.5) / (M + 512) for flank, q_k for other),
               as RNA Bind-n-Seq's R is a frequency ratio.
  z            two-proportion z: (n_k/N - p_bg) / sqrt(p (1 - p) (1/N + 1/M)), p the pooled frequency, M the background's
               5-mer positions (for other, those of the regions the foreground has). Overlapping 5-mers and neighbouring
               windows are not independent, so z ranks k-mers; it is not a calibrated p-value.
The flank background keeps the gene and its local composition but not UV crosslinking's preference for U; the other
background removes what all eCLIP shares (crosslinking, region, windows many RBPs share), but depends on which RBPs the
panel holds. eclip_5mer_backgrounds.tsv lists the five 5-mers with the highest z per background (ranking by R alone puts
rare 5-mers with a few counts first).
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import importlib.util
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

_spec = importlib.util.spec_from_file_location('stage43', Path(__file__).resolve().parent / '43_stage_nucleic_acid_binding.py')
s43 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s43)
log, md5_of = s43.log, s43.md5_of

SUPERSERIES = 'GSE315407'
GEO_TEXT = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=brief'
GEO_SOFT = 'https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{acc}/soft/{acc}_family.soft.gz'
GENCODE_DIR = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/'
GENOME = GENCODE_DIR + 'GRCh38.primary_assembly.genome.fa.gz'
UNIPROT_GENE = ('https://rest.uniprot.org/uniprotkb/search?query=%28gene_exact%3A{g}%29+AND+%28organism_id%3A9606%29+AND+'
                '%28reviewed%3Atrue%29&fields=accession,gene_primary,gene_synonym&format=tsv')
WINDOWS_RE = re.compile(r'(GSE\d+)_(.+)_(K562|HepG2)_(\d+)\.reproducible_enriched_windows\.tsv\.gz$')
LIBRARY_RE = re.compile(r'(.+)_(K562|HepG2)_(\d+)_(IP|INPUT)_rep(\d)$')
REGIONS = s43.REGIONS + ['intergenic']
TOP_WINDOWS, TOP_GENES, KMER_WINDOWS = 200, 200, 1000
K, NK = 5, 4 ** 5
ALPHABET = np.frombuffer(b'ACGU' + b'N' * 252, dtype=np.uint8)
TO_CODES = bytes.maketrans(b'ACGTacgt', b'\x00\x01\x02\x03\x00\x01\x02\x03')
KMERS = [''.join('ACGU'[(i >> (2 * (K - 1 - j))) & 3] for j in range(K)) for i in range(NK)]
CITATION = ('Her et al., Cell Syst 17:101588 (2026), doi:10.1016/j.cels.2026.101588; Skipper: Boyle et al., Cell Genom 3:100317 (2023), '
            'doi:10.1016/j.xgen.2023.100317')
LICENCE = ('GEO: public, no licence stated by the submitters (NCBI FTP README: all data there public and unrestricted); '
           'article CC BY 4.0 (Crossref, version of record, from 2026-03-23)')


# ---- GEO --------------------------------------------------------------------------------------------------------------

def geo_fields(dl, acc: str) -> dict:
    """The !Series_ fields of a GEO series (brief text form), each a list of values."""
    fields = defaultdict(list)
    for line in dl.get(GEO_TEXT.format(acc=acc), f'geo/{acc}_brief.txt', 'cache').read_text(encoding='utf-8').splitlines():
        if line.startswith('!Series_'):
            k, _, v = line[len('!Series_'):].partition(' = ')
            fields[k.strip()].append(v.strip())
    return fields


def geo_samples(path: Path) -> tuple[dict, Counter]:
    """(Yeo-lab number -> {('IP' | 'INPUT', replicate): GSM}, the data-processing lines) from a family SOFT file."""
    by_number, processing, cur = defaultdict(dict), Counter(), None
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line.startswith('^SAMPLE'):
                cur = line.split('=', 1)[1].strip()
            elif cur and line.startswith('!Sample_description') and 'Library name:' in line:
                m = LIBRARY_RE.match(line.split('Library name:', 1)[1].strip())
                if not m:
                    raise SystemExit(f'{cur}: library name not understood: {line.strip()}')
                by_number[m.group(3)][(m.group(4), int(m.group(5)))] = cur
            elif cur and line.startswith('!Sample_data_processing'):
                processing[line.split('=', 1)[1].strip()] += 1
    return by_number, processing


# ---- genome -----------------------------------------------------------------------------------------------------------

class Genome:
    """GRCh38 chromosomes (chr*) of a FASTA, kept as one file of base codes per chromosome (A C G T = 0 1 2 3, anything
    else >= 4) beside the cache, made once and read through memory maps."""

    def __init__(self, fasta_gz: Path, folder: Path):
        folder.mkdir(parents=True, exist_ok=True)
        marker, source = folder / 'source.md5', md5_of(fasta_gz)
        if not (marker.exists() and marker.read_text().split()[0] == source):
            self._convert(fasta_gz, folder)
            marker.write_text(f'{source}  {fasta_gz.name}\n')
        self.seq = {p.stem: np.memmap(p, dtype=np.uint8, mode='r') for p in sorted(folder.glob('chr*.codes'))}

    @staticmethod
    def _convert(fasta_gz: Path, folder: Path):
        name, lines = None, []

        def flush():
            if name and name.startswith('chr'):
                (folder / f'{name}.codes').write_bytes(b''.join(lines).translate(TO_CODES))
                log(f'  genome: {name} {sum(map(len, lines)):,} nt')

        with gzip.open(fasta_gz, 'rb') as fh:
            for line in fh:
                if line.startswith(b'>'):
                    flush()
                    name, lines = line[1:].split()[0].decode(), []
                else:
                    lines.append(line.rstrip(b'\r\n'))
        flush()

    def length(self, chrom: str) -> int:
        return len(self.seq[chrom]) if chrom in self.seq else 0

    def codes(self, chrom: str, start: int, end: int, strand: str) -> np.ndarray:
        """Codes of [start, end) (0-based) read 5' to 3' on the strand."""
        a = np.array(self.seq[chrom][max(start, 0):min(end, self.length(chrom))])
        if strand == '-':
            a = np.where(a < 4, 3 - a, 4).astype(np.uint8)[::-1]
        return a


def rna(codes: np.ndarray) -> str:
    return ALPHABET[codes].tobytes().decode()


def kmer_counts(arrays: list) -> np.ndarray:
    """Counts of every 5-mer (index base 4, A C G U = 0 1 2 3, first base most significant) inside the sequences, none with N."""
    if not arrays:
        return np.zeros(NK, dtype=np.int64)
    sep = np.full(1, 4, dtype=np.uint8)
    s = np.concatenate([x for a in arrays for x in (a, sep)]).astype(np.int64)
    if len(s) < K:
        return np.zeros(NK, dtype=np.int64)
    bad = np.concatenate([[0], np.cumsum(s >= 4)])
    s = np.minimum(s, 3)
    idx = sum(s[j:len(s) - K + 1 + j] << (2 * (K - 1 - j)) for j in range(K))
    ok = (bad[K:] - bad[:-K]) == 0
    return np.bincount(idx[ok], minlength=NK)


class Mask:
    """Merged intervals per (chromosome, strand); free() gives the parts of an interval outside them."""

    def __init__(self, intervals):
        by = defaultdict(list)
        for chrom, strand, a, b in intervals:
            by[(chrom, strand)].append((a, b))
        self.m = {}
        for key, iv in by.items():
            iv.sort()
            merged = []
            for a, b in iv:
                if merged and a <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], b)
                else:
                    merged.append([a, b])
            arr = np.array(merged, dtype=np.int64)
            self.m[key] = (arr[:, 0], arr[:, 1])

    def free(self, chrom, strand, a, b):
        if a >= b:
            return []
        st, en = self.m.get((chrom, strand), (np.empty(0, np.int64), np.empty(0, np.int64)))
        pieces, cur = [], a
        for k in range(int(np.searchsorted(en, a, 'right')), int(np.searchsorted(st, b, 'left'))):
            if st[k] > cur:
                pieces.append((cur, min(int(st[k]), b)))
            cur = max(cur, int(en[k]))
        if cur < b:
            pieces.append((cur, b))
        return pieces


def kmer_profile(genome: Genome, fg: list, mask: Mask) -> tuple[np.ndarray, np.ndarray, int]:
    """For foreground sequences [(chromosome, strand, start, end, region index)]: 5-mer counts per region (6 x 1024), counts
    in the flanks (1024) and the number of sequences read."""
    by_region, flanks, n = defaultdict(list), [], 0
    for chrom, strand, a, b, region in fg:
        if chrom not in genome.seq:
            continue
        n += 1
        by_region[region].append(genome.codes(chrom, a, b, strand))
        span = b - a
        for fa, fb in ((a - span, a), (b, b + span)):
            for pa, pb in mask.free(chrom, strand, max(fa, 0), min(fb, genome.length(chrom))):
                flanks.append(genome.codes(chrom, pa, pb, strand))
    fg_r = np.zeros((len(REGIONS), NK), dtype=np.int64)
    for r, arrays in by_region.items():
        fg_r[r] = kmer_counts(arrays)
    return fg_r, kmer_counts(flanks), n


def two_prop_z(n, N, bg_freq, M):
    p = (n + bg_freq * M) / (N + M)
    se = np.sqrt(np.maximum(p * (1 - p) * (1 / N + 1 / M), 1e-300))
    return (n / N - bg_freq) / se


def enrichment(fg_r: np.ndarray, flank: np.ndarray, pool_r: np.ndarray) -> dict:
    """R and z of every 5-mer against the flank and the region-weighted pool of other RBPs (see the module docstring)."""
    n = fg_r.sum(0)
    N = int(n.sum())
    out = {'n': n, 'N': N, 'M_flank': int(flank.sum()), 'M_other': 0}
    if N == 0:
        return out
    f = (n + 0.5) / (N + NK / 2)
    M = int(flank.sum())
    if M:
        g = (flank + 0.5) / (M + NK / 2)
        out['R_flank'], out['z_flank'] = f / g, two_prop_z(n, N, flank / M, M)
    w = fg_r.sum(1) / N
    Or = pool_r.sum(1)
    q = sum(w[r] * (pool_r[r] + 0.5) / (Or[r] + NK / 2) for r in range(len(REGIONS)) if w[r] > 0)
    Mo = int(Or[w > 0].sum())
    if Mo:
        out['R_other'], out['z_other'], out['M_other'] = f / q, two_prop_z(n, N, q, Mo), Mo
    return out


def top_kmers(e: dict, which: str, n: int = 5) -> str:
    """The n 5-mers with the highest z (ranking by R alone favours rare 5-mers with a few counts), as kmer:R."""
    if f'R_{which}' not in e:
        return ''
    R, z = e[f'R_{which}'], e[f'z_{which}']
    order = [k for k in np.argsort(-z, kind='stable') if z[k] > 0][:n]
    return ';'.join(f'{KMERS[k]}:{R[k]:.2f}' for k in order)


def kmer_rows(key: str, e: dict) -> list:
    if not e['N']:
        return []
    get = lambda name, k, nd: round(float(e[name][k]), nd) if name in e else ''  # noqa: E731
    return [[key, KMERS[k], int(e['n'][k]), get('R_flank', k, 3), get('z_flank', k, 2), get('R_other', k, 3), get('z_other', k, 2)]
            for k in range(NK)]


def pools(profiles: dict, rbp_of: dict) -> dict:
    """For each experiment: the region counts of all other experiments whose RBP differs."""
    total = sum(p[0] for p in profiles.values())
    out = {}
    for exp in profiles:
        same = [x for x in profiles if rbp_of[x] & rbp_of[exp]]
        out[exp] = total - sum(profiles[x][0] for x in same)
    return out


# ---- checks -----------------------------------------------------------------------------------------------------------

def splice_check(genome: Genome, segs: dict, sample: int = 20000) -> tuple[int, int]:
    """(introns read, of which GU..AG on the transcript) for an even sample of GENCODE v50 basic introns."""
    introns = sorted({(c, s, a, b) for (c, s), rows in segs.items() for a, b, region, _ in rows if region == 4 and c in genome.seq})
    step = max(1, len(introns) // sample)
    good = n = 0
    for c, s, a, b in introns[::step]:
        # the intron is [a, b], 1-based: its first two and last two nucleotides on the transcript
        five, three = ((a - 1, a + 1), (b - 2, b)) if s == '+' else ((b - 2, b), (a - 1, a + 1))
        n += 1
        good += rna(genome.codes(c, *five, s)) == 'GU' and rna(genome.codes(c, *three, s)) == 'AG'
    return n, good


def start_check(starts_by_key: dict, windows: pd.DataFrame) -> tuple[int, int]:
    """Windows whose start + 1, and whose start, is the first base of a GENCODE v50 segment (exon part or intron)."""
    plus1 = same = 0
    for (chrom, strand), d in windows.groupby(['chr', 'strand']):
        arr = starts_by_key.get((chrom, strand))
        if arr is None:
            continue
        plus1 += int(np.isin(d['start'].to_numpy() + 1, arr).sum())
        same += int(np.isin(d['start'].to_numpy(), arr).sum())
    return plus1, same


# ---- ENCORE -----------------------------------------------------------------------------------------------------------

def entry_names(entries: Path) -> dict:
    with open(entries, encoding='utf-8') as fh:
        return {r['Entry']: r['Entry Name'] for r in csv.DictReader(fh, delimiter='\t')}


def by_symbol(ref, names: dict, label: str) -> set:
    """The reviewed entries of a primary gene symbol; when there are several (RBM10: the protein and the MINAS-60
    microprotein of the same gene), the one named <SYMBOL>_HUMAN if there is exactly one."""
    accs = ref.symbol(label)
    named = {a for a in accs if names.get(a, '').upper() == f'{label.upper()}_HUMAN'}
    return named if len(accs) > 1 and len(named) == 1 else accs


def encore_experiments(dl, ref, names: dict) -> tuple[list, dict]:
    """The window files of the SubSeries with their samples, RBP mapping and ENCODE experiment; and the GEO fields."""
    sup = geo_fields(dl, SUPERSERIES)
    subs = [v.split(':', 1)[1].strip() for v in sup['relation'] if v.startswith('SuperSeries of:')]
    q = '&'.join(['type=Experiment', 'assay_title=eCLIP', 'limit=all', 'format=json'] +
                 [f'field={f}' for f in ('accession', 'aliases', 'status', 'target.label', 'target.genes.dbxrefs', 'target.genes.symbol',
                                         'biosample_ontology.term_name')])
    by_yeo = defaultdict(list)
    for e in dl.json(f'{s43.ENCODE}/search/?{q}', 'encode/eclip_experiments_all_statuses.json')['@graph']:
        for a in e.get('aliases') or []:
            m = re.match(r'gene-yeo:(\d+)$', a)
            if m:
                by_yeo[(m.group(1), (e.get('biosample_ontology') or {}).get('term_name', ''))].append(e)
    exps, meta = [], {'series': {}, 'processing': Counter()}
    for acc in subs:
        f = geo_fields(dl, acc)
        meta['series'][acc] = f
        stub = acc[:-3] + 'nnn'
        samples, processing = geo_samples(dl.get(GEO_SOFT.format(stub=stub, acc=acc), f'geo/{acc}_family.soft.gz', 'cache'))
        meta['processing'] += processing
        for url in f['supplementary_file']:
            m = WINDOWS_RE.search(url)
            if not m:
                continue
            gse, label, cell, number = m.groups()
            s = samples.get(number, {})
            hits = sorted(by_yeo.get((number, cell), []), key=lambda e: (e.get('status') != 'released', e['accession']))
            enc = hits[0] if hits else None
            accs, how = set(), ''
            if enc:
                accs, how = ref.encode_target(enc.get('target'))
                how = 'encode_' + how
            if not accs:
                accs, how = by_symbol(ref, names, label), 'symbol'
                if accs != ref.symbol(label):
                    how = 'symbol (entry named for it)'
            if not accs:
                p = dl.get(UNIPROT_GENE.format(g=label), f'uniprot/gene_exact_{label}.tsv', 'cache')
                with open(p, encoding='utf-8') as fh:
                    found = {r['Entry'] for r in csv.DictReader(fh, delimiter='\t')} & ref.accs
                accs, how = (found, 'uniprot_gene_name') if len(found) == 1 else (set(), 'unmapped')
            if enc and how.startswith('encode_') and ref.symbol(label) and ref.symbol(label) != accs:
                log(f'  note: {label} ({number}) maps to {sorted(accs)} through ENCODE and {sorted(ref.symbol(label))} by symbol')
            exps.append({'key': f'{label}_{cell}_{number}', 'gse': gse, 'file': m.group(0), 'url': url.replace('ftp://', 'https://', 1),
                         'label': label, 'cell': cell, 'number': number, 'accs': sorted(accs), 'how': how,
                         'ip': ';'.join(s.get(('IP', r), '') for r in (1, 2)), 'input': ';'.join(s.get(('INPUT', r), '') for r in (1, 2)),
                         'encode': enc['accession'] if enc else '', 'encode_status': enc.get('status', '') if enc else ''})
    return sorted(exps, key=lambda x: x['key']), meta


WINDOW_COLUMNS = ['chr', 'start', 'end', 'strand', 'feature_type_top', 'gene_name', 'gene_id', 'gene_type_top', 'input_sum', 'clip_sum',
                  'enrichment_l2or_min', 'enrichment_l2or_mean', 'enrichment_l2or_max', 'q_max']
CHROM_IX: dict = {}


def window_keys(d: pd.DataFrame) -> np.ndarray:
    """One integer per window: chromosome, start, width and strand."""
    c = d['chr'].map(lambda s: CHROM_IX.setdefault(s, len(CHROM_IX))).to_numpy(np.int64)
    width = (d['end'] - d['start']).to_numpy(np.int64)
    if (width < 0).any() or (width > 255).any():
        raise SystemExit('a window wider than 255 nt')
    return ((c * 2 ** 28 + d['start'].to_numpy(np.int64)) * 256 + width) * 2 + (d['strand'] == '-').to_numpy(np.int64)


def stage_encore(dl, ref, names, genome, segs, genes, n_kmer: int, checks: dict):
    exps, meta = encore_experiments(dl, ref, names)
    starts = defaultdict(list)
    for key, rows in segs.items():
        starts[key] = np.unique(np.array([a for a, _, _, _ in rows], dtype=np.int64))
    erows, rrows, frows, grows, wrows, profiles, all_keys = [], [], [], [], [], {}, []
    plus1 = same = 0
    for x in exps:
        path = dl.get(x['url'], f"windows/{x['file']}", 'cache')
        text = ['chr', 'strand', 'feature_type_top', 'gene_name', 'gene_id', 'gene_type_top']
        d = pd.read_csv(path, sep='\t', usecols=WINDOW_COLUMNS, dtype={c: str for c in text})
        d[text] = d[text].fillna('')
        d['row'] = np.arange(len(d))
        d = d.sort_values(['enrichment_l2or_mean', 'enrichment_l2or_min', 'row'], ascending=[False, False, True], kind='stable').reset_index(drop=True)
        a, b = check = start_check(starts, d)
        plus1, same = plus1 + a, same + b
        mids = (d['start'] + (d['end'] - d['start']) // 2 + 1).tolist()
        placed = s43.place_all(segs, genes, list(zip(d['chr'], d['strand'], mids)))
        region = [5 if h is None else h[0] for h in placed]
        d['region'] = region
        x['n'] = len(d)
        counts = Counter(region)
        rrows += [[x['key'], REGIONS[r], counts.get(r, 0)] for r in range(len(REGIONS))]
        frows += [[x['key'], ft, n] for ft, n in sorted(Counter(d['feature_type_top']).items(), key=lambda kv: (-kv[1], kv[0]))]
        per_gene = defaultdict(lambda: [0, float('-inf'), set()])
        for h, score in zip(placed, d['enrichment_l2or_mean']):
            if h is None:
                continue
            pg = per_gene[h[1]]
            pg[0] += 1
            pg[1] = max(pg[1], score)
            pg[2].add(h[0])
        top = sorted(per_gene.items(), key=lambda kv: (-kv[1][0], -kv[1][1], genes[kv[0]][1]))[:TOP_GENES]
        grows += [[x['key'], genes[g][0], genes[g][1], genes[g][2], v[0], round(v[1], 3), ';'.join(s43.REGIONS[r] for r in sorted(v[2]))] for g, v in top]
        k = window_keys(d)
        all_keys.append(np.unique(k))
        for rank, (r, key) in enumerate(zip(d.head(TOP_WINDOWS).itertuples(index=False), k[:TOP_WINDOWS]), 1):
            seq = rna(genome.codes(r.chr, r.start, r.end, r.strand)) if r.chr in genome.seq else ''
            wrows.append([x['key'], rank, r.chr, r.start, r.end, r.strand, r.gene_name, r.gene_id, r.gene_type_top, r.feature_type_top,
                          REGIONS[r.region], round(r.enrichment_l2or_mean, 3), round(r.enrichment_l2or_min, 3), round(r.enrichment_l2or_max, 3),
                          f'{r.q_max:.3g}', r.clip_sum, r.input_sum, key, seq])
        mask = Mask(zip(d['chr'], d['strand'], d['start'], d['end']))
        head = d.head(n_kmer)
        fg = list(zip(head['chr'], head['strand'], head['start'], head['end'], head['region']))
        profiles[x['key']] = kmer_profile(genome, fg, mask)
        log(f"  {x['key']}: {len(d):,} windows; starts at a GENCODE segment start as 0-based {check[0]:,}, as 1-based {check[1]:,}")
    checks['encore_starts'] = (plus1, same)
    # in how many experiments each top window (same coordinates and strand) is a reproducible enriched window
    uniq, n_exp = np.unique(np.concatenate(all_keys), return_counts=True)
    for w in wrows:
        w[-2] = int(n_exp[np.searchsorted(uniq, w[-2])])
    checks['shared'] = (sum(1 for w in wrows if w[-2] >= 10), sum(1 for w in wrows if w[-3] == 0), len(wrows))
    rbp_of = {x['key']: set(x['accs']) or {x['label']} for x in exps}
    pool = pools(profiles, rbp_of)
    krows, brows = [], []
    for x in exps:
        fg_r, flank, n = profiles[x['key']]
        e = enrichment(fg_r, flank, pool[x['key']])
        krows += kmer_rows(x['key'], e)
        brows.append(['ENCORE', x['key'], ';'.join(x['accs']), x['label'], x['cell'], x['n'], n, e['N'], e['M_flank'], e['M_other'],
                      sum(1 for y in exps if not rbp_of[y['key']] & rbp_of[x['key']]), top_kmers(e, 'other'), top_kmers(e, 'flank')])
    for x in exps:
        for acc in x['accs'] or ['']:
            erows.append([x['key'], x['gse'], x['file'], x['label'], acc, x['how'], x['cell'], x['ip'], x['input'], x['n'], x['encode'], x['encode_status']])
    return exps, meta, erows, rrows, frows, grows, wrows, krows, brows


# ---- ENCODE peaks -----------------------------------------------------------------------------------------------------

def encode_kmers(eclip_cache: Path, folder: Path, genome: Genome, segs, genes, n_kmer: int):
    """The same 5-mer enrichment for the ENCODE eCLIP reproducible peaks of encode_eclip_experiments.tsv, from the files
    scripts/43 downloaded (each checked against the MD5 raw_files.tsv records)."""
    with open(folder / 'encode_eclip_experiments.tsv', encoding='utf-8') as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    with open(folder / 'raw_files.tsv', encoding='utf-8') as fh:
        md5 = {r['file']: r['md5'] for r in csv.DictReader(fh, delimiter='\t')}
    exps = {}
    for r in rows:
        x = exps.setdefault(r['encode_experiment'], {'file': r['encode_file'], 'accs': set(), 'target': r['target'], 'cell': r['cell_line']})
        x['accs'].add(r['uniprot_accession'])
    tops, masks, info = {}, {}, {}
    for exp, x in sorted(exps.items()):
        name = f"eclip/{x['file']}.bed.gz"
        path = eclip_cache / name
        if md5_of(path) != md5.get(name):
            raise SystemExit(f'{path}: MD5 differs from raw_files.tsv')
        peaks = []
        with gzip.open(path, 'rt') as fh:
            for line in fh:
                c = line.rstrip('\n').split('\t')
                if len(c) >= 8:
                    peaks.append((c[0], c[5], int(c[1]), int(c[2]), float(c[6]), float(c[7])))
        order = sorted(range(len(peaks)), key=lambda i: (-peaks[i][4], -peaks[i][5], peaks[i][0], peaks[i][2], peaks[i][1]))[:n_kmer]
        tops[exp] = [peaks[i][:4] for i in order]
        masks[exp] = Mask(p[:4] for p in peaks)
        info[exp] = len(peaks)
    # one sweep over GENCODE for every experiment's top peaks
    flat = [(exp, p) for exp in sorted(tops) for p in tops[exp]]
    placed = s43.place_all(segs, genes, [(c, s, a + (b - a) // 2 + 1) for _, (c, s, a, b) in flat])
    fgs = defaultdict(list)
    for (exp, (c, s, a, b)), h in zip(flat, placed):
        fgs[exp].append((c, s, a, b, 5 if h is None else h[0]))
    profiles = {exp: kmer_profile(genome, fgs[exp], masks[exp]) for exp in sorted(tops)}
    rbp_of = {exp: x['accs'] for exp, x in exps.items()}
    pool = pools(profiles, rbp_of)
    krows, brows = [], []
    for exp, x in sorted(exps.items()):
        fg_r, flank, n = profiles[exp]
        e = enrichment(fg_r, flank, pool[exp])
        krows += kmer_rows(exp, e)
        brows.append(['ENCODE', exp, ';'.join(sorted(x['accs'])), x['target'], x['cell'], info[exp], n, e['N'], e['M_flank'], e['M_other'],
                      sum(1 for y in exps if not rbp_of[y] & rbp_of[exp]), top_kmers(e, 'other'), top_kmers(e, 'flank')])
    return krows, brows, len(exps)


# ---- main -------------------------------------------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    root = Path(__file__).resolve().parents[1]
    ap.add_argument('--cache', required=True, help='download folder for ENCORE, the genome and GENCODE')
    ap.add_argument('--eclip-cache', default='PTM_data/external/nucleic_acid', help='the cache scripts/43 eclip read (eclip/*.bed.gz)')
    ap.add_argument('--out', default=str(root / 'inputs' / 'nucleic_acid'))
    ap.add_argument('--entries', default=str(root / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03_entries.tsv'))
    ap.add_argument('--manifest', default=str(root / 'inputs' / 'MANIFEST.tsv'))
    ap.add_argument('--kmer-windows', type=int, default=KMER_WINDOWS)
    ap.add_argument('--skip-encode-kmers', action='store_true')
    args = ap.parse_args()
    cache = Path(args.cache)
    dl, ref, out = s43.Downloads(cache), s43.Reference(Path(args.entries)), s43.Out(Path(args.out))
    today = dt.date.today().isoformat()

    sums = {ln.split()[1]: ln.split()[0] for ln in dl.get(GENCODE_DIR + 'MD5SUMS', 'gencode/MD5SUMS', 'cache').read_text().splitlines() if ln.strip()}
    gtf = dl.get(s43.GENCODE, 'gencode/' + Path(s43.GENCODE).name, 'cache')
    fasta = dl.get(GENOME, 'gencode/' + Path(GENOME).name, 'cache')
    for p in (gtf, fasta):
        if md5_of(p) != sums.get(p.name):
            raise SystemExit(f'{p}: MD5 differs from GENCODE MD5SUMS')
    genome = Genome(fasta, cache / 'genome_codes')
    segs, genes = s43.load_gencode(gtf)
    checks = {'splice': splice_check(genome, segs)}
    log(f"GENCODE v50 introns read: {checks['splice'][0]:,}, GU..AG {checks['splice'][1]:,}")

    names = entry_names(Path(args.entries))
    exps, meta, erows, rrows, frows, grows, wrows, krows, brows = stage_encore(dl, ref, names, genome, segs, genes, args.kmer_windows, checks)
    iso = lambda s: dt.datetime.strptime(s.replace('Public on ', ''), '%b %d %Y').date().isoformat()  # noqa: E731
    dates = {acc: (iso(f['status'][0]), iso(f['last_update_date'][0])) for acc, f in meta['series'].items()}
    series = ' and '.join(sorted(meta['series']))
    pipeline = [p for p in meta['processing'] if 'Skipper' in p or p.startswith('Assembly')]
    sp_n, sp_good = checks['splice']
    st1, st0 = checks['encore_starts']
    sh_10, sh_zero, sh_rows = checks['shared']
    version = (f'GEO {SUPERSERIES} ({series}; public {", ".join(sorted({d[0] for d in dates.values()}))}, last updated '
               f'{", ".join(sorted({d[1] for d in dates.values()}))}; {"; ".join(pipeline)}), read {today}')
    base = dict(source=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={SUPERSERIES} (SubSeries {series}, *.reproducible_enriched_windows.tsv.gz)',
                version_or_date=version, citation=CITATION, licence=LICENCE,
                derived_from=f'{len(exps)} GEO window files (encore_raw_files.tsv)')
    gencode_meta = dict(base, source=base['source'] + f'; {s43.GENCODE}', version_or_date=version + '; GENCODE release 50 basic',
                        licence=LICENCE + '; GENCODE open access')
    mapped = Counter(x['how'] for x in exps)
    renamed = sorted({f"{x['label']} as {';'.join(x['accs']) or 'none'} ({x['how']})" for x in exps if not ref.symbol(x['label'])})
    out.table('encore_eclip_experiments.tsv', ['encore_experiment', 'geo_series', 'geo_file', 'rbp', 'uniprot_accession', 'mapped_by', 'cell_line',
                                               'geo_samples_ip', 'geo_samples_input', 'n_windows', 'encode_experiment', 'encode_status'], erows,
              content='ENCORE eCLIP experiments (Skipper reproducible enriched windows, K562 or HepG2) with the RBP as the file names it, its reviewed '
                      'UniProt entry, the GEO samples (IP replicates 1;2, size-matched input replicates 1;2), the number of windows, and the ENCODE '
                      'eCLIP experiment with the same Yeo-lab number (alias gene-yeo:<number>, same cell line) when ENCODE has one, with its status',
              how_made_or_checked=('scripts/43e_stage_encore.py; RBP mapped through that ENCODE experiment\'s target (its UniProt cross-references, '
                                   'as scripts/43 eclip), else by primary gene symbol, else by the one reviewed human entry UniProt finds for the '
                                   'name (gene_exact): ' + ', '.join(f'{k} {v}' for k, v in sorted(mapped.items())) + '; names that are not a '
                                   'primary gene symbol: ' + ', '.join(renamed)),
              **base)
    out.table('encore_eclip_window_regions.tsv', ['encore_experiment', 'region', 'n_windows'], rrows,
              content='Where the middle of each ENCORE reproducible enriched window falls on GENCODE v50 basic transcripts of its strand: CDS, UTR3, UTR5, '
                      'noncoding_exon, intron (that order when it falls in more than one), or intergenic; placed as scripts/43 places ENCODE eCLIP peaks',
              how_made_or_checked='scripts/43e_stage_encore.py (scripts/43 load_gencode, place_all)', **gencode_meta)
    out.table('encore_eclip_window_features.tsv', ['encore_experiment', 'feature_type', 'n_windows'], frows,
              content='ENCORE reproducible enriched windows per Skipper feature type, as the files give it (feature_type_top: CDS, CDS_START, CDS_STOP, '
                      'UTR5, UTR3, INTRON, splice-site windows SS5/SS3/SSB _ADJ and _PROX, EXON_MRNA, EXON_LNCRNA, EXON_SMALL, EXON_PSEUDO, PRIMIRNA ...)',
              how_made_or_checked='scripts/43e_stage_encore.py', **base)
    out.table('encore_eclip_target_genes.tsv', ['encore_experiment', 'gene_id', 'gene_name', 'gene_type', 'n_windows', 'max_l2or_mean', 'regions'], grows,
              content=f'The {TOP_GENES} genes (GENCODE v50 basic, gene of the window middle) with most ENCORE reproducible enriched windows per experiment '
                      '(ties by the highest mean log2 odds ratio), with the regions the windows fall in; as encode_eclip_target_genes.tsv',
              how_made_or_checked='scripts/43e_stage_encore.py', **gencode_meta)
    out.table('encore_eclip_top_windows.tsv', ['encore_experiment', 'rank', 'chrom', 'start', 'end', 'strand', 'gene_name', 'gene_id', 'gene_type',
                                               'feature_type', 'region', 'l2or_mean', 'l2or_min', 'l2or_max', 'q_max', 'clip_sum', 'input_sum',
                                               'n_experiments', 'rna_sequence'], wrows,
              content=f'The {TOP_WINDOWS} ENCORE reproducible enriched windows per experiment with the highest mean log2 odds ratio of IP over input '
                      '(the order of the files; ties by the lower replicate\'s), with coordinates (GRCh38, 0-based start, end exclusive), strand, the '
                      'file\'s gene, gene type and Skipper feature type, the GENCODE v50 region of the window middle, the replicates\' log2 odds ratios, '
                      'the larger replicate q value, the IP and input read sums the file gives, the number of ENCORE experiments (of '
                      f'{len(exps)}) in which the same window is a reproducible enriched window (a window many unrelated RBPs share is not specific to '
                      'this one), and the RNA sequence of the window (5\' to 3\' on its strand, U alphabet) from the GRCh38 primary assembly of '
                      'GENCODE release 50',
              source=gencode_meta['source'] + f'; {GENOME}', version_or_date=gencode_meta['version_or_date'], citation=CITATION,
              licence=gencode_meta['licence'], derived_from=base['derived_from'] + '; GRCh38.primary_assembly.genome.fa.gz',
              how_made_or_checked=(f'scripts/43e_stage_encore.py (Genome); genome MD5 checked against GENCODE MD5SUMS; {sp_good:,} of {sp_n:,} '
                                   f'GENCODE v50 introns read GU..AG through the same code; {st1:,} window starts + 1 fall on a GENCODE v50 '
                                   f'segment start against {st0:,} starts as given (so the files are 0-based); of the {sh_rows:,} rows, '
                                   f'{sh_zero:,} have no input reads (input_sum 0) and {sh_10:,} are windows reproducible in 10 or more '
                                   'experiments'))
    kmeta = dict(content='', source=gencode_meta['source'] + f'; {GENOME}', version_or_date=gencode_meta['version_or_date'], citation=CITATION,
                 licence=gencode_meta['licence'])
    method = (f'5-mer enrichment of the top {args.kmer_windows:,} windows (by mean log2 odds ratio; all when fewer) per experiment, every overlapping '
              '5-mer on the RNA strand counted (n_fg), against two backgrounds: flank (the same length 5\' and 3\' of each window on its strand, less '
              'any reproducible window of the experiment) and other (the top windows of the other RBPs of the source, reweighted to the experiment\'s '
              'mix of GENCODE v50 regions); R = frequency ratio with 0.5 pseudocounts, z = two-proportion z (a ranking, not a calibrated p); U alphabet')
    out.table('encore_eclip_5mer_enrichment.tsv', ['encore_experiment', 'kmer', 'n_fg', 'R_flank', 'z_flank', 'R_other', 'z_other'], krows,
              **dict(kmeta, content='ENCORE ' + method), derived_from=base['derived_from'] + '; GRCh38.primary_assembly.genome.fa.gz',
              how_made_or_checked='scripts/43e_stage_encore.py (kmer_profile, enrichment); method in the script\'s docstring')
    if not args.skip_encode_kmers:
        ek, eb, n_enc = encode_kmers(Path(args.eclip_cache), Path(args.out), genome, segs, genes, args.kmer_windows)
        out.table('encode_eclip_5mer_enrichment.tsv', ['encode_experiment', 'kmer', 'n_fg', 'R_flank', 'z_flank', 'R_other', 'z_other'], ek,
                  content=('ENCODE eCLIP reproducible peaks: ' + method.replace('windows (by mean log2 odds ratio', 'peaks (by signalValue, then -log10 p')
                           .replace('of each window', 'of each peak').replace('any reproducible window', 'any reproducible peak')
                           .replace('top windows of the other RBPs', 'top peaks of the other RBPs')),
                  source=f'{s43.ENCODE} (eCLIP reproducible peaks, as encode_eclip_experiments.tsv); {s43.GENCODE}; {GENOME}',
                  version_or_date=f'ENCODE files of raw_files.tsv; GENCODE release 50; computed {today}',
                  citation='Van Nostrand et al., Nature 583:711 (2020), doi:10.1038/s41586-020-2077-3; Frankish et al., Nucleic Acids Res (GENCODE)',
                  licence=s43.ENCODE_TERMS + '; GENCODE open access',
                  derived_from=f'{n_enc} ENCODE peak files (raw_files.tsv, MD5 checked); GRCh38.primary_assembly.genome.fa.gz',
                  how_made_or_checked='scripts/43e_stage_encore.py encode_kmers (the same code as the ENCORE 5-mers)')
        brows += eb
    out.table('eclip_5mer_backgrounds.tsv', ['source', 'experiment', 'uniprot_accession', 'rbp', 'cell_line', 'n_sequences_all', 'n_sequences_fg',
                                             'fg_positions', 'flank_positions', 'other_positions', 'other_experiments', 'top_other', 'top_flank'], brows,
              content=('Per eCLIP experiment (ENCORE windows, ENCODE peaks) behind the 5-mer tables: reproducible windows or peaks, those used as '
                       'foreground, and the 5-mer positions of the foreground, the flanks and the other-RBP background (regions of the foreground only), '
                       'the number of other experiments pooled, and the five 5-mers with the highest z against each background, with their R '
                       '(kmer:R)'),
              source=kmeta['source'], version_or_date=kmeta['version_or_date'], citation=CITATION + '; Van Nostrand et al., Nature 583:711 (2020)',
              licence=LICENCE + '; ' + s43.ENCODE_TERMS + '; GENCODE open access', derived_from='the 5-mer tables',
              how_made_or_checked='scripts/43e_stage_encore.py')

    raw = Path(args.out) / 'encore_raw_files.tsv'
    rows = sorted(dl.rows())
    with open(raw, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['file', 'url', 'bytes', 'md5', 'downloaded', 'kept'])
        w.writerows(rows)
    out.manifest.append({'file': f'{Path(args.out).name}/{raw.name}', 'bytes': raw.stat().st_size, 'md5': md5_of(raw), 'rows': len(rows),
                         'redistribution': 'git', 'content': 'Every file scripts/43e downloaded (GEO series records, SOFT files and window files; the '
                                                             'ENCODE alias search; UniProt name lookups; GENCODE release 50 GTF, genome and MD5SUMS), with '
                                                             'its URL, size, MD5 and date; the files themselves are not kept here',
                         'source': 'the sources above', 'version_or_date': today, 'how_made_or_checked': 'scripts/43e_stage_encore.py'})
    out.write_manifest(Path(args.manifest))
    n_acc = len({a for x in exps for a in x['accs']})
    log(f'ENCORE: {len(exps)} experiments of {n_acc} proteins; mapped by {dict(mapped)}')


if __name__ == '__main__':
    main()

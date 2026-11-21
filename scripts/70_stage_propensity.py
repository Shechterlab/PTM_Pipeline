#!/usr/bin/env python3
"""Stage the inputs of the per-site mutation propensity (release v1.7.0): each reviewed protein's coding sequence with the
genomic neighbors and CpG methylation of its bases, the mutation rates (TCGA and DepMap spectra, gnomAD's germline rates),
each tumor sample's and cell line's coding mutation burden, and the TCGA and DepMap changes at every residue of the histone
H3 proteins (the oncohistone controls).

    python scripts/70_stage_propensity.py --genome Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
        --chain GRCh38_to_GRCh37.chain.gz --methylation methylation.txt.bgz --rates mutation_rate_methylation_bins.txt.bgz \
        --maf mc3.v0.2.8.PUBLIC.maf.gz --cases gdc_tcga_cases.tsv --depmap-mutations OmicsSomaticMutations.csv \
        --depmap-models Model.csv --crispr-gene-effect CRISPRGeneEffect_header.csv --ensembl-gtf Homo_sapiens.GRCh38.105.gtf.gz \
        --hgnc hgnc_complete_set.txt --uniprot-ensembl ... --uniprot-isoforms ... --uniprot-release 2026_04 --outdir inputs/layers

Coding sequences
  For each reviewed entry, the Ensembl 105 transcripts UniProt assigns to its displayed isoform (scripts/44b, as scripts/66
  reads them). A transcript counts when its coding sequence (the GTF's CDS records on the GRCh38 primary assembly; the stop
  codon is not in them) translates with the standard code to the entry's 2026_03 canonical sequence. Entries with no such
  transcript are listed with the reason and get no weights. Each base's neighbors are read on the genome, so a base at the
  edge of an exon has its intron neighbor.
  CpG methylation (for the germline rates): at each C followed by G and each G preceded by C (transcript strand, genomic
  neighbors), the mean methylation in gnomAD's methylation table (the v2.1.1 flagship paper's, GRCh37:
  papers/2019-flagship-lof/v1.0/model/methylation.txt.bgz) at the position lifted to GRCh37 with Ensembl's GRCh38_to_GRCh37
  chain, binned as gnomAD's constraint code bins it (gnomad_lof, constraint_utils/constraint_basics.py, prepare_ht: level 2
  above 0.6, level 1 above 0.2, else 0). A position without a lift or a value is level 0, as there.

Rates
  Germline: gnomAD's mutation_rate_methylation_bins (Karczewski et al., Nature 581:434, 2020; gnomAD data are CC0), the
    mutation rate per trinucleotide, substitution and methylation level. Copied with its bytes unchanged.
  Somatic: per trinucleotide and substitution (pyrimidine-centered, 96 classes), the SNVs that change the residue (missense
    or stop) in the coding sequences above, over the substitutions there that would (each genomic position once; a position
    counts for a substitution if it changes the residue in any of its transcripts). TCGA: calls kept by the rules of
    scripts/66 (FILTER PASS; OV and LAML also wga alone), placed on GRCh38 through the lift of the coding positions.
    DepMap: the SNVs of OmicsSomaticMutations.csv (GRCh38). Each source with all its samples, and without the hypermutated
    ones.

Burden
  Per tumor sample (every TCGA sample in the MAF with a GDC project, as the denominators of scripts/66) and cell line (every
  DepMap model with mutation calls): coding calls per Mb of coding sequence. Coding calls: TCGA, the kept calls of a coding
  Variant_Classification (Missense_Mutation, Nonsense_Mutation, Nonstop_Mutation, Translation_Start_Site, Silent,
  Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins); DepMap, the calls with a coding consequence (missense,
  stop gained or lost, start lost, frameshift, in-frame insertion or deletion, protein altering). OmicsSomaticMutations.csv
  keeps synonymous calls only where they touch a splice region, so DepMap burdens count protein-changing calls almost only;
  the TCGA table also gives its non-silent count. The Mb are the summed length of the coding sequences above (each genomic
  position once). Hypermutated: more than 10 per Mb (Campbell et al., Cell 171:1042, 2017).

H3
  The TCGA and DepMap changes at every residue (not only R, K, S, T and Y) of the entries named Histone H3*, by the rules of
  scripts/66, for the controls at H3K27, H3K36 and H3G34 (UniProt numbering counts the initiator methionine: K28, K37, G35).

Outputs (in --outdir; gzip without a timestamp, so the same input gives the same bytes)
  propensity_cds.tsv.gz                 one row per matched coding sequence (transcripts with the same coding exons share a
                                        row): entry, transcripts, gene, chromosome, strand, coding sequence (transcript strand,
                                        no stop codon), the bases before and after it, the intron bases at each exon junction,
                                        and the CpG methylation levels above 0
  propensity_cds_status.tsv             every reviewed entry: the transcripts UniProt assigns, those that match, and why
                                        none does
  gnomad_v2_1_1_mutation_rates.txt.bgz  gnomAD's rate table, the bytes downloaded
  propensity_somatic_spectra.tsv        per source, cohort (all; no_hypermutated) and the 96 classes: SNVs observed,
                                        opportunities, rate
  propensity_sample_burden.tsv          per TCGA sample and DepMap line: group, coding calls by kind, burden per Mb,
                                        hypermutated or not
  tcga_mc3_changes_h3_all_residues.tsv.gz, depmap_24q4_changes_h3_all_residues.tsv.gz   the H3 changes at every residue
  propensity_inputs_summary.json        inputs (MD5), rules and counts
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import importlib.util
import io
import json
import re
import shutil
import subprocess
from argparse import Namespace
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from common import md5_file, parse_fasta, refresh_manifest_row, update_manifest, upsert_tsv_rows
from ensembl_genes import GTF_105, HGNC, HGNC_URL, transcript_genes

HERE = Path(__file__).resolve().parent


def _load(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, HERE / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


s66 = _load('stage_cancer_changes', '66_stage_cancer_changes.py')
s44b, s44 = s66.s44b, s66.s44

GENOME_URL = 'https://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
CHAIN_URL = 'https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh38_to_GRCh37.chain.gz'
GNOMAD_MODEL = 'https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/model/'
METHYLATION_URL = GNOMAD_MODEL + 'methylation.txt.bgz'
RATES_URL = GNOMAD_MODEL + 'mutation_rate_methylation_bins.txt.bgz'
BINNING = ('gnomad_lof constraint_utils/constraint_basics.py, prepare_ht (https://github.com/macarthur-lab/gnomad_lof): '
           'methylation_level 2 where a CpG transition\'s methylation MEAN > 0.6, 1 where > 0.2, else 0')
GNOMAD_CIT = 'Karczewski et al., Nature 581:434 (2020), doi:10.1038/s41586-020-2308-7'
CAMPBELL = 'Campbell et al., Cell 171:1042 (2017), doi:10.1016/j.cell.2017.09.048'
HYPERMUTATED_PER_MB = 10.0
CHROMS = [str(c) for c in range(1, 23)] + ['X', 'Y']
B = 'ACGT'
IDX = {b: i for i, b in enumerate(B)}
COMP = str.maketrans('ACGTN', 'TGCAN')
CODON_ORDER = 'TCAG'
CODE = dict(zip((a + b + c for a in CODON_ORDER for b in CODON_ORDER for c in CODON_ORDER),
                'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
TCGA_CODING = {'Missense_Mutation': 'missense', 'Nonsense_Mutation': 'stop', 'Nonstop_Mutation': 'other', 'Translation_Start_Site': 'other',
               'Silent': 'silent', 'Frame_Shift_Del': 'indel', 'Frame_Shift_Ins': 'indel', 'In_Frame_Del': 'indel', 'In_Frame_Ins': 'indel'}
DEPMAP_CODING = [('missense_variant', 'missense'), ('stop_gained', 'stop'), ('frameshift_variant', 'indel'), ('inframe_deletion', 'indel'),
                 ('inframe_insertion', 'indel'), ('stop_lost', 'other'), ('start_lost', 'other'), ('protein_altering_variant', 'other'),
                 ('synonymous_variant', 'silent')]
BURDEN_COLUMNS = ['source', 'sample_id', 'group', 'coding_calls', 'non_silent_calls', 'missense', 'stop', 'silent', 'indel', 'other',
                  'snvs_in_coding_sequences', 'burden_per_mb', 'hypermutated', 'threshold_per_mb']
SPECTRA_COLUMNS = ['source', 'cohort', 'samples', 'context', 'ref', 'alt', 'observed', 'opportunity', 'rate']


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def write_gz(df: pd.DataFrame, path: Path) -> None:
    buf = io.StringIO()
    df.to_csv(buf, sep='\t', index=False, lineterminator='\n')
    with open(path, 'wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0, compresslevel=9) as gz:
        gz.write(buf.getvalue().encode('utf-8'))


def translate(cds: str) -> str:
    return ''.join(CODE.get(cds[i:i + 3], 'X') for i in range(0, len(cds) - 2, 3))


def collapse(context: str, alt: str) -> tuple[str, str, str]:
    """Pyrimidine-centered (context, ref, alt)."""
    if context[1] in 'AG':
        return rc(context), context[1].translate(COMP), alt.translate(COMP)
    return context, context[1], alt


# the 96 classes, in a fixed order, and the class of each (genome-oriented trinucleotide code, alt base)
CLASSES = [(l + m + r, m, a) for m in 'CT' for l in B for r in B for a in B if a != m]
CLASS_INDEX = {c: k for k, c in enumerate(CLASSES)}
CLASS_OF = np.full(64 * 4, -1, dtype=np.int16)
for _t in range(64):
    _ctx = B[_t // 16] + B[(_t // 4) % 4] + B[_t % 4]
    for _a in range(4):
        if B[_a] != _ctx[1]:
            CLASS_OF[_t * 4 + _a] = CLASS_INDEX[collapse(_ctx, B[_a])]


def read_genome(path: str) -> dict[str, np.ndarray]:
    """GRCh38 chromosomes 1-22, X, Y as uint8 codes (A C G T = 0 1 2 3, anything else 4)."""
    lut = np.full(256, 4, dtype=np.uint8)
    for b, i in IDX.items():
        lut[ord(b)] = lut[ord(b.lower())] = i
    out, name, parts = {}, None, []

    def flush():
        if name in CHROMS:
            out[name] = lut[np.frombuffer(b''.join(parts), dtype=np.uint8)]

    with gzip.open(path, 'rb') as fh:
        for line in fh:
            if line.startswith(b'>'):
                flush()
                name, parts = line[1:].split()[0].decode(), []
            elif name in CHROMS:
                parts.append(line.rstrip(b'\n'))
    flush()
    return out


class Chain:
    """Point lift-over with a UCSC-format chain file (Ensembl's GRCh38_to_GRCh37); 1-based positions in and out."""

    def __init__(self, path: str):
        blocks = defaultdict(list)
        with gzip.open(path, 'rt') as fh:
            head = None
            for line in fh:
                f = line.split()
                if not f:
                    continue
                if f[0] == 'chain':
                    head = {'t': f[2], 'ts': int(f[5]), 'q': f[7], 'qsize': int(f[8]), 'qstrand': f[9], 'qs': int(f[10])}
                    t, q = head['ts'], head['qs']
                    continue
                size = int(f[0])
                blocks[head['t']].append((t, t + size, head['q'], q, head['qstrand'] == '-', head['qsize']))
                if len(f) == 3:
                    t, q = t + size + int(f[1]), q + size + int(f[2])
        self.by = {}
        self.overlaps = 0
        for chrom, bl in blocks.items():
            bl.sort()
            a = np.array([(x[0], x[1], x[3], x[4], x[5]) for x in bl], dtype=np.int64)
            self.overlaps += int((a[1:, 0] < a[:-1, 1]).sum())
            self.by[chrom] = (a, [x[2] for x in bl])

    def lift(self, chrom: str, pos: np.ndarray):
        """(GRCh37 chromosome per position or '', GRCh37 position, strand flipped) for 1-based GRCh38 positions."""
        n = len(pos)
        qchrom, qpos, flip = np.array([''] * n, dtype=object), np.zeros(n, dtype=np.int64), np.zeros(n, dtype=bool)
        if chrom not in self.by or not n:
            return qchrom, qpos, flip
        a, names = self.by[chrom]
        p0 = pos - 1
        k = np.searchsorted(a[:, 0], p0, side='right') - 1
        ok = (k >= 0) & (p0 < a[np.clip(k, 0, None), 1])
        for back in (1, 2, 3):                  # blocks of different chains can overlap: try the few starting just before
            k2 = k - back
            alt = ~ok & (k2 >= 0) & (p0 < a[np.clip(k2, 0, None), 1])
            k, ok = np.where(alt, k2, k), ok | alt
        kk = k[ok]
        off = p0[ok] - a[kk, 0]
        q = a[kk, 2] + off
        minus = a[kk, 3].astype(bool)
        q = np.where(minus, a[kk, 4] - 1 - q, q)
        qpos[ok], flip[ok] = q + 1, minus
        qchrom[ok] = np.array(names, dtype=object)[kk]
        return qchrom, qpos, flip


def gtf_cds(gtf: str, wanted: set[str]) -> dict[str, dict]:
    tx = {}
    rx = re.compile(r'(gene_id|transcript_id|gene_name) "([^"]+)"')
    with gzip.open(gtf, 'rt', encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if f[2] != 'CDS':
                continue
            at = dict(rx.findall(f[8]))
            t = at.get('transcript_id', '')
            if t not in wanted:
                continue
            d = tx.setdefault(t, {'chrom': f[0], 'strand': f[6], 'gene_id': at.get('gene_id', ''), 'gene_name': at.get('gene_name', ''), 'exons': []})
            d['exons'].append((int(f[3]), int(f[4])))
    return tx


def build_models(args, canonical: dict, genome: dict, counts: Counter):
    """Matched transcripts (with positions) and every entry's status."""
    canon, xref = s44b.canonical_transcripts(args.uniprot_ensembl, args.uniprot_isoforms)
    counts.update({f'uniprot_{k}': v for k, v in xref.items()})
    of_acc = defaultdict(set)
    for t, accs in canon.items():
        for a in accs:
            of_acc[a].add(t)
    tx = gtf_cds(args.ensembl_gtf, set(canon))
    models, status = [], []
    for acc in sorted(canonical):
        assigned = sorted(of_acc.get(acc, ()))
        matched, why, same = [], Counter(), {}
        for t in assigned:
            d = tx.get(t)
            if d is None:
                why['not_in_ensembl_105'] += 1
                continue
            if d['chrom'] not in genome:
                why['not_on_chromosomes_1_22_X_Y'] += 1
                continue
            exons = sorted(d['exons'])
            pos = np.concatenate([np.arange(s, e + 1, dtype=np.int64) for s, e in exons])
            if d['strand'] == '-':
                pos = pos[::-1]
            codes = genome[d['chrom']][pos - 1]
            if d['strand'] == '-':
                codes = np.where(codes < 4, 3 - codes, 4)
            cds = ''.join('ACGTN'[c] for c in codes)
            protein = translate(cds) if len(cds) % 3 == 0 else ''
            if protein.endswith('*') and protein[:-1] == canonical[acc]:
                cds, pos, protein = cds[:-3], pos[:-3], protein[:-1]
            if protein != canonical[acc]:
                why['translation_differs'] += 1
                continue
            matched.append(t)
            key = (d['chrom'], d['strand'], tuple(exons))
            if key in same:                     # another transcript with the same coding exons: one coding sequence, both named
                same[key]['transcripts'].append(t)
                continue
            same[key] = {'acc': acc, 'transcript': t, 'transcripts': [t], 'gene_id': d['gene_id'], 'gene_name': d['gene_name'],
                         'chrom': d['chrom'], 'strand': d['strand'], 'cds': cds, 'pos': pos}
            models.append(same[key])
        if matched:
            st = 'matched'
        elif not assigned:
            st = 'no_transcript_assigned'
        else:
            st = ';'.join(sorted(why))
        counts[f'entries_{st}' if st in ('matched', 'no_transcript_assigned') else 'entries_no_transcript_matches'] += 1
        status.append({'canonical_UniProtAC': acc, 'status': st, 'transcripts_assigned': ';'.join(assigned), 'transcripts_matched': ';'.join(matched)})
    return models, status


def neighbors(m: dict, genome: dict) -> tuple[np.ndarray, np.ndarray]:
    """Transcript-strand left and right neighbor codes of every coding base, read on the genome."""
    g = genome[m['chrom']]
    p = m['pos']
    if m['strand'] == '+':
        left, right = g[p - 2], g[p]
    else:
        left, right = g[p], g[p - 2]
        left, right = np.where(left < 4, 3 - left, 4), np.where(right < 4, 3 - right, 4)
    return left, right


def methylation_levels(path: str, wanted: dict[str, np.ndarray]) -> dict[str, dict[int, float]]:
    """MEAN methylation at the GRCh37 positions wanted (per chromosome), streamed from gnomAD's table."""
    keys = {c: np.unique(v) for c, v in wanted.items() if len(v)}
    out = {c: {} for c in keys}
    proc = subprocess.Popen(['bash', '-c', f'gzip -dc "{path}" | tr ":" "\\t"'], stdout=subprocess.PIPE)
    for chunk in pd.read_csv(proc.stdout, sep='\t', names=['chrom', 'pos', 'mean'], skiprows=1, dtype={'chrom': str, 'pos': np.int64, 'mean': float},
                             chunksize=5_000_000):
        for chrom, sub in chunk.groupby('chrom', sort=False):
            if chrom not in keys:
                continue
            hit = np.isin(sub['pos'].to_numpy(), keys[chrom])
            if hit.any():
                out[chrom].update(zip(sub['pos'].to_numpy()[hit].tolist(), sub['mean'].to_numpy()[hit].tolist()))
    proc.wait()
    return out


def level(mean: float | None) -> int:
    if mean is None:
        return 0
    return 2 if mean > 0.6 else 1 if mean > 0.2 else 0


# NS[codon_code, offset, alt] = 1 where replacing the codon's base at offset by alt changes the residue (missense or stop)
NS = np.zeros((64, 3, 4), dtype=bool)
for _c in range(64):
    _codon = B[_c // 16] + B[(_c // 4) % 4] + B[_c % 4]
    for _o in range(3):
        for _a in range(4):
            if B[_a] != _codon[_o]:
                _new = _codon[:_o] + B[_a] + _codon[_o + 1:]
                NS[_c, _o, _a] = CODE[_codon] != CODE[_new]


def coding_positions(models: list[dict], genome: dict):
    """Per chromosome, the unique coding positions (sorted), each one's genome-oriented trinucleotide code, and a 4-bit mask of
    the genome-oriented alt bases that change the residue in any transcript."""
    per = defaultdict(list)
    for m in models:
        codes = np.array([IDX.get(b, 4) for b in m['cds']], dtype=np.int64)
        cod = codes.reshape(-1, 3)
        cc = cod[:, 0] * 16 + cod[:, 1] * 4 + cod[:, 2]
        ns = NS[np.repeat(cc, 3), np.tile(np.arange(3), len(cc))]        # (L, 4), transcript-strand alts
        if m['strand'] == '-':
            ns = ns[:, ::-1]                                               # alt t on the transcript is 3 - t on the genome
        mask = (ns * (1 << np.arange(4))).sum(axis=1)
        per[m['chrom']].append((m['pos'], mask))
    out = {}
    for chrom, parts in per.items():
        pos = np.concatenate([p for p, _ in parts])
        mask = np.concatenate([k for _, k in parts])
        order = np.argsort(pos, kind='stable')
        pos, mask = pos[order], mask[order]
        u, start = np.unique(pos, return_index=True)
        merged = np.bitwise_or.reduceat(mask, start)
        g = genome[chrom]
        tri = g[u - 2].astype(np.int64) * 16 + g[u - 1].astype(np.int64) * 4 + g[u].astype(np.int64)
        valid = (g[u - 2] < 4) & (g[u - 1] < 4) & (g[u] < 4)
        out[chrom] = {'pos': u, 'mask': merged, 'tri': np.where(valid, tri, -1)}
    return out


def opportunity(cp: dict) -> np.ndarray:
    opp = np.zeros(len(CLASSES), dtype=np.int64)
    for d in cp.values():
        ok = d['tri'] >= 0
        for a in range(4):
            sel = ok & ((d['mask'] >> a) & 1).astype(bool)
            cls = CLASS_OF[d['tri'][sel] * 4 + a]
            opp += np.bincount(cls[cls >= 0], minlength=len(CLASSES))
    return opp


def snv_class(cp: dict, chrom: str, pos: int, ref: str, alt: str) -> int:
    """The class of a genome-oriented GRCh38 SNV if it lies in a coding position and changes the residue there; -2 if the
    reference base differs from the genome's; -1 otherwise."""
    d = cp.get(chrom)
    if d is None or ref not in IDX or alt not in IDX:
        return -1
    k = np.searchsorted(d['pos'], pos)
    if k >= len(d['pos']) or d['pos'][k] != pos or d['tri'][k] < 0:
        return -1
    tri = int(d['tri'][k])
    if (tri // 4) % 4 != IDX[ref]:
        return -2
    if not (int(d['mask'][k]) >> IDX[alt]) & 1:
        return -1
    return int(CLASS_OF[tri * 4 + IDX[alt]])


def tcga_pass(args, cp: dict, lifted: dict, counts: Counter):
    """Per TCGA sample: project, coding calls by kind, and the residue-changing SNVs per class."""
    project_of = {}
    with open(args.cases, newline='', encoding='utf-8') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            project_of[r['submitter_id']] = r['project.project_id']
    group, calls, classes = {}, defaultdict(dict), defaultdict(Counter)
    with gzip.open(args.maf, 'rt', encoding='utf-8', newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        ix = {c: k for k, c in enumerate(next(reader))}
        for r in reader:
            barcode = r[ix['Tumor_Sample_Barcode']]
            sample, project = barcode[:15], project_of.get(barcode[:12], '')
            if not project:
                continue
            group[sample] = project
            if not s44b.kept(r[ix['FILTER']], project):
                continue
            chrom, start = r[ix['Chromosome']], int(r[ix['Start_Position']])
            ref, alt = r[ix['Reference_Allele']], r[ix['Tumor_Seq_Allele2']]
            kind = TCGA_CODING.get(r[ix['Variant_Classification']])
            if kind:
                calls[sample][(chrom, start, ref, alt)] = kind
            if r[ix['Variant_Type']] != 'SNP' or len(ref) != 1 or len(alt) != 1:
                continue
            hit = lifted.get(chrom)
            if hit is None:
                continue
            k = np.searchsorted(hit['key'], start)
            if k >= len(hit['key']) or hit['key'][k] != start:
                continue
            c38, p38, flip = hit['chrom'][k], int(hit['pos'][k]), bool(hit['flip'][k])
            if flip:
                ref, alt = ref.translate(COMP), alt.translate(COMP)
            cls = snv_class(cp, c38, p38, ref, alt)
            if cls == -2:
                counts['tcga_snv_reference_differs_from_grch38'] += 1
            elif cls >= 0:
                classes[sample][(c38, p38, alt)] = cls
    return group, calls, classes


def lift_coding(cp: dict, chain: Chain, counts: Counter) -> dict:
    """GRCh37 chromosome -> sorted GRCh37 positions of the coding positions, with their GRCh38 chromosome, position and strand."""
    parts = defaultdict(list)
    for chrom, d in cp.items():
        qc, qp, fl = chain.lift(chrom, d['pos'])
        ok = qc != ''
        counts['coding_positions'] += len(d['pos'])
        counts['coding_positions_lifted_to_grch37'] += int(ok.sum())
        for c in sorted(set(qc[ok])):
            s = ok & (qc == c)
            parts[c].append((qp[s], np.array([chrom] * int(s.sum()), dtype=object), d['pos'][s], fl[s]))
    out = {}
    for c, ps in parts.items():
        key = np.concatenate([p[0] for p in ps])
        order = np.argsort(key, kind='stable')
        out[c] = {'key': key[order], 'chrom': np.concatenate([p[1] for p in ps])[order], 'pos': np.concatenate([p[2] for p in ps])[order],
                  'flip': np.concatenate([p[3] for p in ps])[order]}
    return out


def depmap_pass(args, cp: dict, counts: Counter):
    mut = pd.read_csv(args.depmap_mutations, dtype=str, keep_default_na=False, usecols=['Chrom', 'Pos', 'Ref', 'Alt', 'VariantType', 'VariantInfo', 'ModelID'])
    mut = mut.drop_duplicates(['ModelID', 'Chrom', 'Pos', 'Ref', 'Alt'])
    calls, classes = defaultdict(dict), defaultdict(Counter)
    for model, chrom, pos, ref, alt, vtype, info in zip(mut['ModelID'], mut['Chrom'], mut['Pos'], mut['Ref'], mut['Alt'], mut['VariantType'], mut['VariantInfo']):
        terms = set(info.split('&'))
        kind = next((k for t, k in DEPMAP_CODING if t in terms), None)
        if kind:
            calls[model][(chrom, pos, ref, alt)] = kind
        if vtype != 'SNV' or len(ref) != 1 or len(alt) != 1:
            continue
        c = chrom[3:] if chrom.startswith('chr') else chrom
        cls = snv_class(cp, c, int(pos), ref, alt)
        if cls == -2:
            counts['depmap_snv_reference_differs_from_grch38'] += 1
        elif cls >= 0:
            classes[model][(c, int(pos), alt)] = cls
    models = sorted(set(mut['ModelID']))
    return models, calls, classes


def burden_rows(source: str, groups: dict, calls: dict, classes: dict, mb: float) -> list[dict]:
    rows = []
    for sid in sorted(groups):
        kinds = Counter(calls.get(sid, {}).values())
        coding = sum(kinds.values())
        per_mb = coding / mb
        rows.append({'source': source, 'sample_id': sid, 'group': groups[sid], 'coding_calls': coding, 'non_silent_calls': coding - kinds['silent'],
                     'missense': kinds['missense'], 'stop': kinds['stop'], 'silent': kinds['silent'], 'indel': kinds['indel'], 'other': kinds['other'],
                     'snvs_in_coding_sequences': len(classes.get(sid, {})), 'burden_per_mb': round(per_mb, 4),
                     'hypermutated': per_mb > HYPERMUTATED_PER_MB, 'threshold_per_mb': HYPERMUTATED_PER_MB})
    return rows


def spectra_rows(source: str, classes: dict, hyper: set, samples: list, opp: np.ndarray) -> list[dict]:
    rows = []
    for cohort, keep in (('all', samples), ('no_hypermutated', [s for s in samples if s not in hyper])):
        obs = np.zeros(len(CLASSES), dtype=np.int64)
        for s in keep:
            for cls in classes.get(s, {}).values():
                obs[cls] += 1
        for k, (ctx, ref, alt) in enumerate(CLASSES):
            rows.append({'source': source, 'cohort': cohort, 'samples': len(keep), 'context': ctx, 'ref': ref, 'alt': alt, 'observed': int(obs[k]),
                         'opportunity': int(opp[k]), 'rate': f'{obs[k] / opp[k]:.6e}' if opp[k] else ''})
    return rows


def model_rows(models: list[dict], genome: dict, chain: Chain, args, tx_gene: dict, hgnc: HGNC, counts: Counter) -> pd.DataFrame:
    """The coding-sequence table: flanks, junction bases and CpG methylation levels."""
    cpg = []          # (model index, cds index, chrom, GRCh38 position)
    flank = []
    for i, m in enumerate(models):
        left, right = neighbors(m, genome)
        cds = m['cds']
        cds_codes = np.array([IDX.get(b, 4) for b in cds])
        jumps = np.nonzero(np.abs(np.diff(m['pos'])) != 1)[0] + 1      # first coding index of each later exon
        junctions = ','.join(f'{j}:{"ACGTN"[right[j - 1]]}{"ACGTN"[left[j]]}' for j in jumps)
        flank.append(("ACGTN"[left[0]], "ACGTN"[right[-1]], junctions))
        is_c = (cds_codes == 1) & (right == 2)
        is_g = (cds_codes == 2) & (left == 1)
        for j in np.nonzero(is_c | is_g)[0]:
            cpg.append((i, int(j), m['chrom'], int(m['pos'][j])))
    counts['cpg_coding_bases'] = len(cpg)
    wanted, lift_of = defaultdict(list), {}
    by_chrom = defaultdict(list)
    for x in cpg:
        by_chrom[x[2]].append(x[3])
    for chrom, ps in by_chrom.items():
        ps = np.unique(np.array(ps, dtype=np.int64))
        qc, qp, _ = chain.lift(chrom, ps)
        for p, c, q in zip(ps.tolist(), qc.tolist(), qp.tolist()):
            if c:
                lift_of[(chrom, p)] = (c, q)
                wanted[c].append(q)
    counts['cpg_coding_positions'] = sum(len(np.unique(v)) for v in by_chrom.values())
    counts['cpg_coding_positions_lifted'] = len(lift_of)
    means = methylation_levels(args.methylation, {c: np.array(v, dtype=np.int64) for c, v in wanted.items()})
    levels = defaultdict(list)
    found = Counter()
    for i, j, chrom, p in cpg:
        q = lift_of.get((chrom, p))
        mean = means.get(q[0], {}).get(q[1]) if q else None
        lv = level(mean)
        found['with_value' if mean is not None else 'no_value'] += 1
        found[f'level_{lv}'] += 1
        if lv:
            levels[i].append(f'{j}:{lv}')
    counts.update({f'cpg_coding_bases_{k}': v for k, v in found.items()})
    rows = []
    for i, m in enumerate(models):
        ensg, name = tx_gene.get(m['transcript'], (m['gene_id'], m['gene_name']))
        rows.append({'canonical_UniProtAC': m['acc'], 'transcript_id': m['transcript'], 'same_cds_transcripts': ';'.join(m['transcripts'][1:]),
                     'gene_id': m['gene_id'],
                     'gene': hgnc.symbol(ensg, name=name) or name or ensg, 'chromosome': m['chrom'], 'strand': m['strand'],
                     'cds_start': int(m['pos'].min()), 'cds_end': int(m['pos'].max()), 'cds': m['cds'], 'flank_5': flank[i][0], 'flank_3': flank[i][1],
                     'junctions': flank[i][2], 'cpg_methylation': ','.join(levels.get(i, []))})
    return pd.DataFrame(rows)


def h3_tables(args, canonical: dict, names: dict, counts: Counter, tx_gene: dict, hgnc: HGNC):
    h3 = sorted(a for a, n in names.items() if n.startswith('Histone H3'))
    screened = {m.group(1) for c in pd.read_csv(args.crispr_gene_effect, nrows=0).columns if (m := re.match(r'^.+ \((\d+)\)$', str(c)))}
    placer = s66.Placer(args.fasta, args.entries, args.symbol_entrez, canonical, screened)
    every = list(s66.AA)
    c66 = Counter()
    trows, _, _, _ = s66.tcga(args, canonical, placer, every, c66, tx_gene, hgnc)
    t = pd.DataFrame(trows, columns=s66.TCGA_COLUMNS)
    t = t[t['canonical_UniProtAC'].isin(h3)].copy()
    t.insert(2, 'residue', [canonical[a][p - 1] for a, p in zip(t['canonical_UniProtAC'], t['position'])])
    t = t.sort_values(['canonical_UniProtAC', 'position', 'alt', 'project_id', 'sample_id'], kind='stable')
    drows, _, _, _ = s66.depmap(args, canonical, every, c66, hgnc)
    d = drows[drows['canonical_UniProtAC'].isin(h3)].rename(columns={'ref': 'residue'})
    d = d[['canonical_UniProtAC', 'position', 'residue'] + [c for c in d.columns if c not in ('canonical_UniProtAC', 'position', 'residue')]]
    d = d.sort_values(['canonical_UniProtAC', 'position', 'alt', 'model_id'], kind='stable')
    counts['h3_entries'] = len(h3)
    counts['h3_tcga_changes'] = len(t)
    counts['h3_depmap_changes'] = len(d)
    return h3, t, d


def sha256(p: Path) -> str:
    return s66.sha256(p)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--genome', required=True, help=f'Ensembl 105 GRCh38 primary assembly ({GENOME_URL})')
    ap.add_argument('--chain', required=True, help=f'Ensembl GRCh38 to GRCh37 chain ({CHAIN_URL})')
    ap.add_argument('--methylation', required=True, help=f'gnomAD methylation table, GRCh37 ({METHYLATION_URL})')
    ap.add_argument('--rates', required=True, help=f'gnomAD mutation_rate_methylation_bins ({RATES_URL})')
    ap.add_argument('--maf', required=True)
    ap.add_argument('--cases', required=True)
    ap.add_argument('--depmap-mutations', required=True)
    ap.add_argument('--depmap-models', required=True)
    ap.add_argument('--crispr-gene-effect', required=True)
    ap.add_argument('--ensembl-gtf', required=True)
    ap.add_argument('--hgnc', required=True)
    ap.add_argument('--uniprot-ensembl', required=True)
    ap.add_argument('--uniprot-isoforms', required=True)
    ap.add_argument('--uniprot-release', required=True)
    ap.add_argument('--fasta', default='inputs/reference/uniprot_human_reviewed_2026_03.fasta')
    ap.add_argument('--entries', default='inputs/reference/uniprot_human_reviewed_2026_03_entries.tsv')
    ap.add_argument('--symbol-entrez', default='inputs/disease/orghsegdb_3.22.0_symbol_entrez.tsv')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--manifest', default='')
    ap.add_argument('--downloaded', default=dt.date.today().isoformat())
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    counts = Counter()
    canonical = parse_fasta(args.fasta)
    names = {}
    for line in open(args.fasta, encoding='utf-8'):
        if line.startswith('>'):
            m = re.match(r'^>\w+\|([^|]+)\|\S+ (.*?) OS=', line)
            if m:
                names[m.group(1)] = m.group(2)
    tx_gene, hgnc = transcript_genes(args.ensembl_gtf), HGNC(args.hgnc)

    genome = read_genome(args.genome)
    models, status = build_models(args, canonical, genome, counts)
    chain = Chain(args.chain)
    counts['chain_block_overlaps'] = chain.overlaps
    cds = model_rows(models, genome, chain, args, tx_gene, hgnc, counts)
    write_gz(cds, out / 'propensity_cds.tsv.gz')
    pd.DataFrame(status).to_csv(out / 'propensity_cds_status.tsv', sep='\t', index=False, lineterminator='\n')
    shutil.copyfile(args.rates, out / 'gnomad_v2_1_1_mutation_rates.txt.bgz')

    cp = coding_positions(models, genome)
    coding_bp = int(sum(len(d['pos']) for d in cp.values()))
    mb = coding_bp / 1e6
    opp = opportunity(cp)
    lifted = lift_coding(cp, chain, counts)
    groups, tcalls, tclasses = tcga_pass(args, cp, lifted, counts)
    dmodels, dcalls, dclasses = depmap_pass(args, cp, counts)
    lineage = pd.read_csv(args.depmap_models, dtype=str, keep_default_na=False).set_index('ModelID')['OncotreeLineage']
    dgroups = {m: lineage.get(m, '') for m in dmodels}
    burden = burden_rows('tcga_mc3', groups, tcalls, tclasses, mb) + burden_rows('depmap_24q4', dgroups, dcalls, dclasses, mb)
    bt = pd.DataFrame(burden, columns=BURDEN_COLUMNS)
    bt.to_csv(out / 'propensity_sample_burden.tsv', sep='\t', index=False, lineterminator='\n')
    hyper = set(bt.loc[bt['hypermutated'], 'sample_id'])
    spectra = spectra_rows('tcga_mc3', tclasses, hyper, sorted(groups), opp) + spectra_rows('depmap_24q4', dclasses, hyper, dmodels, opp)
    pd.DataFrame(spectra, columns=SPECTRA_COLUMNS).to_csv(out / 'propensity_somatic_spectra.tsv', sep='\t', index=False, lineterminator='\n')

    h3, ht, hd = h3_tables(args, canonical, names, counts, tx_gene, hgnc)
    write_gz(ht, out / 'tcga_mc3_changes_h3_all_residues.tsv.gz')
    write_gz(hd, out / 'depmap_24q4_changes_h3_all_residues.tsv.gz')

    st = pd.DataFrame(status)
    summary = {
        'generated': dt.date.today().isoformat(),
        'inputs': {k: {'file': Path(getattr(args, k)).name, 'md5': md5_file(getattr(args, k))}
                   for k in ('genome', 'chain', 'methylation', 'rates', 'maf', 'cases', 'depmap_mutations', 'depmap_models', 'crispr_gene_effect',
                             'ensembl_gtf', 'hgnc', 'uniprot_ensembl', 'uniprot_isoforms', 'fasta')},
        'urls': {'genome': GENOME_URL, 'chain': CHAIN_URL, 'methylation': METHYLATION_URL, 'rates': RATES_URL, 'ensembl_gtf': GTF_105, 'hgnc': HGNC_URL},
        'uniprot_release': args.uniprot_release,
        'rules': {
            'coding_sequence': 'Ensembl 105 transcripts UniProt assigns to the displayed isoform (scripts/44b canonical_transcripts) whose GTF CDS '
                               'records on the GRCh38 primary assembly translate with the standard code to the 2026_03 canonical sequence',
            'methylation': 'gnomAD v2.1.1 flagship methylation table (GRCh37) at the GRCh38 position lifted with Ensembl GRCh38_to_GRCh37; ' + BINNING,
            'germline_rates': 'gnomAD mutation_rate_methylation_bins (' + GNOMAD_CIT + '), CC0, bytes unchanged',
            'somatic_spectra': 'residue-changing SNVs in the coding positions over the residue-changing substitutions there, 96 pyrimidine-centered '
                               'classes; TCGA kept calls (scripts/66 rules) placed on GRCh38 through the lifted coding positions; DepMap SNVs',
            'burden': f'coding calls per Mb of coding sequence ({coding_bp} bp, each GRCh38 position once); hypermutated above '
                      f'{HYPERMUTATED_PER_MB:g} per Mb ({CAMPBELL}); DepMap OmicsSomaticMutations.csv keeps synonymous calls only where they touch a splice region',
            'h3': 'TCGA and DepMap changes at every residue of the entries named Histone H3*, by the rules of scripts/66'},
        'coding_bp': coding_bp,
        'entries': len(canonical),
        'entries_by_status': dict(Counter(st['status']).most_common()),
        'coding_sequences': len(cds),
        'matched_transcripts': int(sum(len(m['transcripts']) for m in models)),
        'tcga_samples': len(groups), 'tcga_hypermutated': int(bt[(bt['source'] == 'tcga_mc3') & bt['hypermutated']].shape[0]),
        'depmap_lines': len(dmodels), 'depmap_hypermutated': int(bt[(bt['source'] == 'depmap_24q4') & bt['hypermutated']].shape[0]),
        'h3_entries': h3,
        'counts': dict(sorted(counts.items())),
        'md5': {p.name: md5_file(p) for p in sorted(out.glob('propensity_*')) if p.suffix != '.json'} |
               {n: md5_file(out / n) for n in ('gnomad_v2_1_1_mutation_rates.txt.bgz', 'tcga_mc3_changes_h3_all_residues.tsv.gz',
                                               'depmap_24q4_changes_h3_all_residues.tsv.gz')},
    }
    (out / 'propensity_inputs_summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    if args.manifest:
        write_manifest(Path(args.manifest), out, args)
    print(json.dumps({k: summary[k] for k in ('coding_bp', 'entries_by_status', 'coding_sequences', 'matched_transcripts', 'tcga_samples', 'tcga_hypermutated',
                                              'depmap_lines', 'depmap_hypermutated', 'counts')}, indent=1))


def write_manifest(manifest: Path, out: Path, args) -> None:
    folder = out.name
    how = 'scripts/70_stage_propensity.py'
    ens_ver = 'Ensembl release 105 (GRCh38 primary assembly, GTF); UniProt 2026_03 canonical sequences'
    tc_src = s44.TC_SRC.replace('scripts/44b_extract_tcga_mc3.py', how)
    rows = [
        {'file': f'{folder}/propensity_cds.tsv.gz', 'content': 'Per reviewed entry, each Ensembl 105 transcript UniProt assigns to its displayed isoform '
         'whose coding sequence translates to the 2026_03 canonical sequence: gene, chromosome, strand, coding sequence (no stop codon), the genomic '
         'bases before and after it and at each exon junction, and the CpG methylation level (gnomAD v2.1.1 methylation, binned 1 or 2) of its '
         'CpG bases (gzip)', 'source': f'{GENOME_URL}; {GTF_105}; {METHYLATION_URL}; {CHAIN_URL}', 'version_or_date': ens_ver,
         'citation': 'Ensembl (Cunningham et al., Nucleic Acids Res 50:D988, 2022); ' + GNOMAD_CIT, 'licence': 'Ensembl: no restrictions; gnomAD: CC0',
         'derived_from': 'Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz; Homo_sapiens.GRCh38.105.gtf.gz; methylation.txt.bgz; GRCh38_to_GRCh37.chain.gz'},
        {'file': f'{folder}/propensity_cds_status.tsv', 'content': 'Every reviewed entry: the Ensembl transcripts UniProt assigns to its displayed '
         'isoform, those whose Ensembl 105 coding sequence translates to the canonical sequence, and why none does',
         'source': GTF_105, 'version_or_date': ens_ver, 'citation': 'Ensembl', 'licence': 'no restrictions', 'derived_from': 'Homo_sapiens.GRCh38.105.gtf.gz'},
        {'file': f'{folder}/gnomad_v2_1_1_mutation_rates.txt.bgz', 'content': 'gnomAD mutation_rate_methylation_bins: germline mutation rate per '
         'trinucleotide, substitution and CpG methylation level (the bytes downloaded)', 'source': RATES_URL, 'version_or_date': 'gnomAD v2.1.1 '
         'flagship paper resources v1.0 (object of 2020-08-20)', 'citation': GNOMAD_CIT, 'licence': 'CC0 (gnomAD)', 'derived_from': ''},
        {'file': f'{folder}/propensity_somatic_spectra.tsv', 'content': 'TCGA MC3 and DepMap 24Q4 somatic spectra, all samples and without the '
         'hypermutated: residue-changing SNVs in the coding sequences of propensity_cds per pyrimidine-centered trinucleotide class, the '
         'residue-changing substitutions there, and their ratio', 'source': tc_src + '; ' + s44.DM_SRC, 'version_or_date': s44.TC_VER + '; ' + s44.DM_VER,
         'citation': s44.TC_CIT + '; ' + s44.DM_CIT, 'licence': s66.TC_LICENSE + '; CC BY 4.0 (DepMap)', 'derived_from': 'mc3.v0.2.8.PUBLIC.maf.gz; OmicsSomaticMutations.csv'},
        {'file': f'{folder}/propensity_sample_burden.tsv', 'content': 'Per TCGA tumor sample and DepMap cell line: project or Oncotree lineage, '
         'coding calls by kind, burden per Mb of coding sequence and the hypermutated flag (more than 10 per Mb, Campbell et al. 2017)',
         'source': tc_src + '; ' + s44.DM_SRC, 'version_or_date': s44.TC_VER + '; ' + s44.DM_VER, 'citation': s44.TC_CIT + '; ' + s44.DM_CIT + '; ' + CAMPBELL,
         'licence': s66.TC_LICENSE + '; CC BY 4.0 (DepMap)', 'derived_from': 'mc3.v0.2.8.PUBLIC.maf.gz; OmicsSomaticMutations.csv; Model.csv'},
        {'file': f'{folder}/tcga_mc3_changes_h3_all_residues.tsv.gz', 'content': 'TCGA MC3 missense and stop changes at every residue of the '
         'histone H3 entries, the rules and columns of tcga_mc3_changes_<X>.tsv.gz, with the residue (gzip)', 'source': tc_src,
         'version_or_date': s44.TC_VER, 'citation': s44.TC_CIT, 'licence': s66.TC_LICENSE, 'derived_from': 'mc3.v0.2.8.PUBLIC.maf.gz'},
        {'file': f'{folder}/depmap_24q4_changes_h3_all_residues.tsv.gz', 'content': 'DepMap 24Q4 missense and stop changes at every residue of the '
         'histone H3 entries, the rules and columns of depmap_24q4_changes_<X>.tsv.gz, with the residue (gzip)', 'source': s44.DM_SRC,
         'version_or_date': s44.DM_VER, 'citation': s44.DM_CIT, 'licence': 'CC BY 4.0', 'derived_from': 'OmicsSomaticMutations.csv'},
        {'file': f'{folder}/propensity_inputs_summary.json', 'content': 'Inputs (MD5), rules and counts of the staging of the propensity inputs',
         'source': 'the downloads above', 'version_or_date': ens_ver, 'citation': '', 'licence': '', 'derived_from': ''},
    ]
    for r in rows:
        path = out / r['file'].split('/', 1)[1]
        n = '' if path.suffix in ('.json', '.bgz') else (s66.gz_rows(path) if path.suffix == '.gz' else sum(1 for _ in open(path, encoding='utf-8')) - 1)
        r.update({'bytes': path.stat().st_size, 'md5': md5_file(path), 'rows': n, 'redistribution': 'git', 'how_made_or_checked': how})
    update_manifest(manifest, rows)
    raw = []
    for path, url, ver, kept in (
            (Path(args.genome), GENOME_URL, 'Ensembl release 105, GRCh38 primary assembly', 'not kept: 881 MB, served by Ensembl at this URL; MD5 and SHA-256 identify it'),
            (Path(args.chain), CHAIN_URL, 'Ensembl assembly mapping, GRCh38 to GRCh37', 'release_asset'),
            (Path(args.methylation), METHYLATION_URL, 'gnomAD v2.1.1 flagship paper resources v1.0, methylation (GRCh37)',
             'not kept: 247 MB, served by gnomAD at this URL (CC0); MD5 and SHA-256 identify it'),
            (Path(args.rates), RATES_URL, 'gnomAD v2.1.1 flagship paper resources v1.0, mutation_rate_methylation_bins', 'in Git as layers/gnomad_v2_1_1_mutation_rates.txt.bgz (same bytes)')):
        raw.append({'file': path.name, 'url': url, 'bytes': path.stat().st_size, 'md5': md5_file(path), 'sha256': sha256(path),
                    'downloaded': args.downloaded, 'version': ver, 'kept': kept})
    upsert_tsv_rows(out / 'raw_files.tsv', raw)
    refresh_manifest_row(manifest, out / 'raw_files.tsv', f'{folder}/raw_files.tsv')


if __name__ == '__main__':
    main()

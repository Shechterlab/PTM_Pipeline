#!/usr/bin/env python3
"""Stage RBR-ID of K562 nuclei (Oksuz et al. 2023; PRIDE PXD035484) as peptide and RNA-binding-region tables.

    python -I scripts/43c_stage_rbrid.py --cache <download dir> --out inputs/interactions \
        [--fasta inputs/reference/uniprot_human_reviewed_2026_03.fasta] [--threshold 3.0] [--min-runs 3]

RBR-ID (He et al., Mol Cell 64:416, 2016, doi:10.1016/j.molcel.2016.09.034) crosslinks RNA to protein in cells
pulsed or not with 4-thiouridine (4SU) and irradiated with 312 nm UV. A peptide crosslinked to RNA changes mass, so
the signal of the unmodified peptide drops in the +4SU samples. Oksuz et al. (Mol Cell 83:2449, 2023,
doi:10.1016/j.molcel.2023.06.012; preprint bioRxiv doi:10.1101/2022.09.27.509776) ran it on K562 nuclei: 5
biological replicates per condition, each acquired in 3 technical replicates, DIA on a Q Exactive HF, searched with
DIA-NN 1.8.1. PRIDE holds the 30 raw files and one DIA-NN report.tsv (2.02 GB) for all 30 runs; the licence is CC0.

What the authors describe, and how this script follows it
---------------------------------------------------------
Oksuz et al. 2023, STAR Methods, "Bioinformatic analysis of the RBR-ID data" (author manuscript PMC10529847):
  - charge states of the same peptide summed; abundances normalised by the median peptide intensity of each run;
  - normalised abundances in cells treated or not with 4SU compared by unpaired, two-sided Student's t tests;
  - peptides missing in all 5 x 3 runs of one treatment tested with Fisher's exact test on detection frequency;
  - p values of both tests adjusted together (Benjamini-Hochberg);
  - "RBR-ID+" peptides: p < 0.05 and log2 fold change < 0 (their Figure S5 analysis).
He et al. 2016 (Experimental Procedures and Supplemental Experimental Procedures, mmc10):
  - depletion = log2(mean normalised +4SU intensity / mean normalised -4SU intensity);
  - score = -log2(+4SU / -4SU) * (log10 p)^2;
  - "For residue level RBR-ID scores, we summed the RBR-ID score of each peptide overlapping any given amino acid
    and smoothed the resulting curve using Friedman's 'super smoother' (Friedman, 1984)."
    The super smoother here is a port of R's stats::supsmu (span = "cv", bass = 0), checked against R.

Choices this script makes where the papers are silent (each is a command-line option)
--------------------------------------------------------------------------------------
  - DIA-NN filter: precursor Q.Value <= 0.01, Lib.Q.Value <= 0.01 and Lib.PG.Q.Value <= 0.01, the filter DIA-NN
    itself applies to its matrices after match-between-runs (--reanalyse). The report holds rows above these.
  - Peptide = stripped sequence: all precursors (charge states, oxidised Met, N-terminal acetyl) summed per run,
    as MaxQuant's peptides table (He et al. 2016) sums them.
  - The t test is on the normalised abundances on their linear scale, as both papers word it (--scale log2 runs it
    on log2 values instead), over the 15 runs of each condition (technical replicates as the authors counted them,
    "5 x 3"). The 5 x 5 test on the means of technical replicates is given beside it (p_5v5).
  - Minimum runs: a t test needs >= --min-runs (3) quantified runs in each condition. A peptide absent from all 15
    runs of one condition gets Fisher's test (as the authors); one seen in 1-2 runs of a condition gets no test.
  - Fisher-tested peptides have no finite ratio, so they get no score and add nothing to residue scores; they keep
    their p, q and direction, and count as RBR-ID+ when seen only in controls with p < 0.05.
  - Score uses the unadjusted p (He et al. 2016); q (BH) is given beside it.
  - Residue scores: only residues covered by at least one scored peptide get a value (He et al. plot the others as
    gaps), and the smoother runs over those positions in one pass per protein.
  - RBR segment: maximal run of consecutive covered residues whose smoothed residue score is >= --threshold (3.0).
    The papers give no residue-level cut-off. A score of 3 is, for one peptide alone, a 3.4-fold depletion at
    p = 0.05, 1.7-fold at p = 0.01 or 1.26-fold at p = 0.001 (overlapping peptides add up). The default was set by a
    rule fixed before any comparison with methylarginines: the lowest of 1, 2, 3 and 5 at which the label-permutation
    estimate of the residue false-discovery rate (below) is under 10%. With the defaults the estimates were 19% at 1,
    10.8% at 2, 7.3% at 3 and 4.2% at 5; the script prints them and the manifest keeps the one at the threshold.
  - Permutation null: the 10 samples split 5 vs 5 with 2 or 3 samples swapped between the groups (the 200 balanced
    splits; technical replicates stay together), --null-permutations (20) of them drawn with --seed (1), each scored as
    the real split. Estimated FDR at a threshold = mean null fraction / observed fraction of scored residues at or
    above it. Swapped splits keep a fifth of any real effect, so the estimate errs high.

Run-to-condition design (worked out here; PRIDE has no SDRF or experimental-design file)
-------------------------------------------------------------------------------------------
  - Raw files: K562_nuclei_DIA_no_4SU_{1..5}.t{1,2,3}.raw and K562_nuclei_DIA_4SU_{6..10}.t{1,2,3}.raw.
  - The DIA-NN log (DIA.DIA-NN.search.log, on PRIDE) loads 30 mzML runs: 220131_RL_RWT190_Nuclei_DIA_{1..10},
    220206_RL_RWT190_Nuclei_DIA_{1..10} and 220226_190_{1..10}_2, i.e. samples 1-10 acquired on three dates.
  - Methods: 5 biological replicates per condition, each acquired in three technical replicates (5 x 3).
  So the sample number N in a run name is the raw files' sample number: N <= 5 is -4SU, N >= 6 is +4SU, and the
  three dates are the three technical replicates (which date is t1, t2 or t3 is not recorded and not needed). The
  script checks that the 30 runs are exactly these, and prints checks from the data. On this report:
  - each run's most correlated run has the same condition in 30 of 30 runs, the same date in 21 and is of the same
    sample in 4 (the date weighs more than the sample);
  - significant changes (p < 0.05) are 9,218 depletions and 10,493 enrichments in +4SU. Median normalisation lifts
    the rest when abundant RNA-binding proteins are depleted;
  - RBR-ID+ peptides: 26.6% over RRM domains, 29.7% elsewhere in RRM proteins, 13.2% in other proteins.
  Not printed but checked: HNRNPC MIAGQVLDINLAAEPK (74-89), the peptide He et al. 2016 show depleted in mouse cells
  (same sequence and numbering), is depleted in samples 6-10 (log2 -0.40, p = 0.004), as are most peptides of its RRM
  (VFIGNLNTLVVK, IVGCSVHK, NARAAVAGEDGR: p < 0.02), while its peptides from residue 188 on are not (p > 0.1).

Search caveat: the DIA-NN search allowed carbamidomethyl-Cys (fixed), Met oxidation and protein N-terminal acetyl
only (log: --unimod4, --var-mod UniMod:35 and UniMod:1). No arginine methylation was searched, so a peptide carrying
a methylarginine is quantified only through its unmethylated form, and its depletion says nothing about the methyl
form. The script counts the UniMod ids in Modified.Sequence to confirm it.

Residue scores are staged only at or above the published cut-off
------------------------------------------------------------------
Neither paper sets a cut-off on residue scores. The published cut-off is the peptide one: He et al. 2016,
Experimental Procedures, "Peptides with log2(fold-change) < 0 and p value < 0.05 (Student's t test) were considered
primary hits"; Oksuz et al. 2023, STAR Methods, "RBR-ID+ peptides (p-val<0.05, log2FC<0)". The residue table
therefore holds a residue, with its residue scores, only when at least one RBR-ID+ peptide covers it (unadjusted
p < 0.05 and depleted in +4SU; a peptide absent from every +4SU run counts with Fisher p < 0.05, the test Oksuz et
al. use for such peptides). He et al. also removed peptides that passed the same cut-off without UV; this dataset has
no no-UV samples, so that filter cannot be applied. The full profile (about 0.7 million scored residues) is not
staged: residue_profiles() recomputes it from the peptide table in seconds. The measured background (which residues,
and so which arginines, any peptide covered) is staged as intervals in the coverage table.

Outputs (in --out):
  rbrid_oksuz2023_k562_peptides.tsv            one row per peptide placed on a reviewed canonical sequence
  rbrid_oksuz2023_k562_positive_residues.tsv   residues under an RBR-ID+ peptide, with their residue scores
  rbrid_oksuz2023_k562_coverage.tsv            measured background: quantified, tested and scored intervals
  rbrid_oksuz2023_k562_segments.tsv            RBR segments (residue score >= --threshold) with summarised scores
  rbrid_oksuz2023_k562_runs.tsv                the 30 runs: date, sample, condition, raw files, normalisation median
  raw_rbrid_files.tsv                          the PRIDE files read, with size, MD5 and date
Their rows in --manifest (inputs/MANIFEST.tsv) are replaced, or added after the last interactions/ row.
The report itself is streamed and not kept: the per-run x precursor aggregate goes to the cache
(report_precursor_runs.tsv.gz, report_precursors.tsv.gz, report_stream_facts.json), and peptides dropped by the
sequence placement to rbrid_dropped_peptides.tsv there.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import math
import re
import time
import urllib.request
import warnings
from collections import Counter, defaultdict, deque
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

UA = 'PTM_Pipeline-stager/1.0 (https://github.com/Shechterlab/PTM_Pipeline)'
PRIDE_API = 'https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD035484'
PRIDE_FTP = 'https://ftp.pride.ebi.ac.uk/pride/data/archive/2023/07/PXD035484'
REPORT = 'report.tsv'
LOG = 'DIA.DIA-NN.search.log'
REPORT_COLUMNS = ['Run', 'Protein.Ids', 'Genes', 'Modified.Sequence', 'Stripped.Sequence', 'Precursor.Id',
                  'Precursor.Charge', 'Q.Value', 'Global.Q.Value', 'Lib.Q.Value', 'PG.Q.Value', 'Lib.PG.Q.Value',
                  'Global.PG.Q.Value', 'Proteotypic', 'Precursor.Quantity', 'Precursor.Normalised']
RUN_PATTERNS = [re.compile(r'^220131_RL_RWT190_Nuclei_DIA_(\d+)$'), re.compile(r'^220206_RL_RWT190_Nuclei_DIA_(\d+)$'),
                re.compile(r'^220226_190_(\d+)_2$')]
BATCHES = ['220131', '220206', '220226']
CITATION = ('Oksuz et al., Mol Cell 83:2449 (2023), doi:10.1016/j.molcel.2023.06.012; RBR-ID: He et al., Mol Cell '
            '64:416 (2016), doi:10.1016/j.molcel.2016.09.034')
LICENCE = 'CC0 (PRIDE project licence, "Creative Commons Public Domain (CC0)", read from the PRIDE API)'
MANIFEST_COLUMNS = ['file', 'bytes', 'md5', 'rows', 'redistribution', 'content', 'source', 'version_or_date', 'citation',
                    'licence', 'derived_from', 'how_made_or_checked', 'original_name']
PEPTIDES = 'rbrid_oksuz2023_k562_peptides.tsv'
SEGMENTS = 'rbrid_oksuz2023_k562_segments.tsv'
RAW_FILES = 'raw_rbrid_files.tsv'
RUNS = 'rbrid_oksuz2023_k562_runs.tsv'
RESIDUES = 'rbrid_oksuz2023_k562_positive_residues.tsv'
COVERAGE = 'rbrid_oksuz2023_k562_coverage.tsv'
POSITIVE_RULE = ('a residue under at least one RBR-ID+ peptide: log2(+4SU/-4SU) < 0 and unadjusted p < 0.05 (Student t '
                 'test), the primary-hit cut-off of He et al. 2016 (Experimental Procedures) and the RBR-ID+ definition '
                 'of Oksuz et al. 2023 (STAR Methods, Figure S5 analysis), or, for a peptide in no +4SU run, Fisher p '
                 '< 0.05 (Oksuz et al. 2023 test such peptides with Fisher\'s exact test). Neither paper sets a '
                 'cut-off on residue scores')


def log(*a):
    print(*a, flush=True)


def md5_of(path: Path) -> str:
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path) -> Path:
    """Download once into the cache (cached thereafter)."""
    if dest.exists() and dest.stat().st_size:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + '.part')
    for k in range(5):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': UA, 'Accept': '*/*'})
            with urllib.request.urlopen(req, timeout=300) as r, open(part, 'wb') as fh:
                while chunk := r.read(1 << 20):
                    fh.write(chunk)
            part.replace(dest)
            return dest
        except Exception as e:  # noqa: BLE001 - retried, then raised
            if k == 4:
                raise
            log(f'  retry {url}: {e}')
            time.sleep(2 ** (k + 1))
    return dest


# ----------------------------------------------------------------------------------------------------------------
# Streaming the 2 GB report
# ----------------------------------------------------------------------------------------------------------------
class RangedStream(io.RawIOBase):
    """A URL read front to back in fixed blocks, each an HTTP Range request; several blocks are fetched at once
    (one connection to PRIDE gives about 1.6 MB/s) and handed on in file order. A failed block is fetched again.
    Every byte handed on goes through MD5 and SHA-1 (PRIDE lists SHA-1 checksums)."""

    def __init__(self, url: str, size: int, block: int = 16 << 20, workers: int = 6):
        super().__init__()
        self.url, self.size, self.pos, self.workers = url, size, 0, workers
        self.md5, self.sha1 = hashlib.md5(), hashlib.sha1()
        self.retries = 0
        self.blocks = [(s, min(size, s + block) - 1) for s in range(0, size, block)]
        self.pool = ThreadPoolExecutor(workers)
        self.pending, self.submitted = deque(), 0
        self.cur, self.off = b'', 0
        self.t0, self.reported = time.time(), 0

    def readable(self):
        return True

    def _fetch(self, start: int, end: int) -> bytes:
        for k in range(10):
            try:
                req = urllib.request.Request(self.url, headers={'User-Agent': UA, 'Range': f'bytes={start}-{end}'})
                with urllib.request.urlopen(req, timeout=180) as r:
                    if r.status != 206:
                        raise IOError(f'HTTP status {r.status} for a range request')
                    data = r.read()
                if len(data) != end - start + 1:
                    raise IOError(f'block {start}-{end}: {len(data)} bytes')
                return data
            except Exception as e:  # noqa: BLE001 - fetched again, then raised
                self.retries += 1
                if k == 9:
                    raise
                log(f'  block at byte {start:,} again: {e}')
                time.sleep(min(60, 2 ** (k + 1)))
        raise IOError('unreachable')

    def readinto(self, b):
        if self.off >= len(self.cur):
            while self.submitted < len(self.blocks) and len(self.pending) < 2 * self.workers:
                self.pending.append(self.pool.submit(self._fetch, *self.blocks[self.submitted]))
                self.submitted += 1
            if not self.pending:
                return 0
            self.cur, self.off = self.pending.popleft().result(), 0
            self.md5.update(self.cur)
            self.sha1.update(self.cur)
            self.pos += len(self.cur)
            if self.pos - self.reported >= 200_000_000:
                self.reported = self.pos
                log(f'  {self.pos / 1e9:.2f} of {self.size / 1e9:.2f} GB read ({self.pos / (time.time() - self.t0) / 1e6:.1f} MB/s)')
        n = min(len(b), len(self.cur) - self.off)
        b[:n] = self.cur[self.off:self.off + n]
        self.off += n
        return n

    def close(self):
        for f in self.pending:
            f.cancel()
        self.pool.shutdown(wait=False)
        super().close()


def aggregate_report(cache: Path, url: str, size: int, sha1_expected: str | None) -> dict:
    """Stream report.tsv once and keep per run x precursor rows (cache/report_precursor_runs.tsv.gz) and one row of
    annotations per precursor (cache/report_precursors.tsv.gz). Returns the stream's facts (MD5, SHA-1, row count,
    UniMod ids seen). The report is not written to disk."""
    facts_path = cache / 'report_stream_facts.json'
    runs_path, prec_path = cache / 'report_precursor_runs.tsv.gz', cache / 'report_precursors.tsv.gz'
    if facts_path.exists() and runs_path.exists() and prec_path.exists():
        return json.loads(facts_path.read_text())
    log(f'Streaming {url} ({size / 1e9:.2f} GB); only the aggregate is kept')
    raw = RangedStream(url, size)
    buf = io.BufferedReader(raw, buffer_size=8 << 20)
    unimods, rows, runs, conflicts = Counter(), 0, Counter(), 0
    prec_info = {}
    tmp_runs = runs_path.with_name(runs_path.name + '.part')
    num_cols = ['Precursor.Quantity', 'Precursor.Normalised', 'Q.Value', 'Global.Q.Value', 'Lib.Q.Value', 'PG.Q.Value',
                'Lib.PG.Q.Value', 'Global.PG.Q.Value']
    with gzip.open(tmp_runs, 'wt', encoding='utf-8', compresslevel=5) as out:
        out.write('\t'.join(['Run', 'Precursor.Id'] + num_cols) + '\n')
        reader = pd.read_csv(buf, sep='\t', usecols=REPORT_COLUMNS, dtype={c: str for c in REPORT_COLUMNS[:7] + ['Proteotypic']},
                             chunksize=400_000, low_memory=False)
        for chunk in reader:
            rows += len(chunk)
            runs.update(chunk['Run'].value_counts().to_dict())
            for ms, n in chunk['Modified.Sequence'].value_counts().items():
                for u in re.findall(r'UniMod:(\d+)', ms):
                    unimods[u] += n
            ann = chunk.drop_duplicates('Precursor.Id')[['Precursor.Id', 'Stripped.Sequence', 'Modified.Sequence',
                                                         'Precursor.Charge', 'Protein.Ids', 'Genes', 'Proteotypic']]
            for r in ann.itertuples(index=False):
                key = r[0]
                val = tuple('' if (isinstance(v, float) and math.isnan(v)) else str(v) for v in r[1:])
                old = prec_info.get(key)
                if old is None:
                    prec_info[key] = val
                elif old != val:  # DIA-NN annotates a precursor the same way in every run; counted, first kept
                    conflicts += 1
            chunk[['Run', 'Precursor.Id'] + num_cols].to_csv(out, sep='\t', header=False, index=False, float_format='%.6g')
            log(f'  {rows:,} rows, {len(prec_info):,} precursors')
    buf.close()
    if raw.pos != size:
        raise IOError(f'read {raw.pos} of {size} bytes')
    sha1, md5 = raw.sha1.hexdigest(), raw.md5.hexdigest()
    if sha1_expected and sha1 != sha1_expected:
        raise IOError(f'SHA-1 {sha1} differs from PRIDE {sha1_expected}')
    tmp_runs.replace(runs_path)
    with gzip.open(prec_path, 'wt', encoding='utf-8') as fh:
        fh.write('\t'.join(['Precursor.Id', 'Stripped.Sequence', 'Modified.Sequence', 'Precursor.Charge', 'Protein.Ids',
                            'Genes', 'Proteotypic']) + '\n')
        for k in sorted(prec_info):
            fh.write('\t'.join((k,) + prec_info[k]) + '\n')
    facts = {'url': url, 'bytes': size, 'md5': md5, 'sha1': sha1, 'sha1_pride': sha1_expected, 'rows': rows,
             'runs': dict(sorted(runs.items())), 'precursors': len(prec_info), 'unimod_counts': dict(unimods),
             'precursor_annotation_conflicts': conflicts,
             'resumes': raw.retries, 'streamed_on': dt.date.today().isoformat()}
    facts_path.write_text(json.dumps(facts, indent=1))
    log(f'  done: {rows:,} rows, MD5 {md5}, SHA-1 {sha1} (PRIDE {sha1_expected}), {raw.retries} resumes')
    return facts


def pride_metadata(cache: Path) -> tuple[dict, list]:
    proj = json.loads(fetch(PRIDE_API, cache / 'pride_project.json').read_text())
    files = json.loads(fetch(PRIDE_API + '/files?pageSize=500&page=0', cache / 'pride_files.json').read_text())
    return proj, files


# ----------------------------------------------------------------------------------------------------------------
# Design: which run is which sample and condition
# ----------------------------------------------------------------------------------------------------------------
RAW_NAME = re.compile(r'^K562_nuclei_DIA_(no_4SU|4SU)_(\d+)\.t([123])\.raw$')


def raw_files_by_sample(files: list) -> dict:
    """The 30 raw files on PRIDE: samples 1-5 without 4SU, 6-10 with, each in technical replicates t1-t3."""
    by = defaultdict(list)
    for f in files:
        if f['fileName'].endswith('.raw'):
            m = RAW_NAME.match(f['fileName'])
            if not m or (m.group(1) == 'no_4SU') != (int(m.group(2)) <= 5):
                raise SystemExit(f'unexpected raw file {f["fileName"]}')
            by[int(m.group(2))].append(f['fileName'])
    if sorted(by) != list(range(1, 11)) or any(len(v) != 3 for v in by.values()):
        raise SystemExit('expected raw files for samples 1-10, three each')
    return {k: sorted(v) for k, v in by.items()}


def search_log(path: Path) -> tuple[list, str]:
    """Runs the DIA-NN log loads, and its command line."""
    text = path.read_text(encoding='utf-8', errors='replace')
    runs = sorted(set(re.findall(r'Loading run \S*/([^/\s]+)\.mzML', text)))
    cmd = next((ln.strip() for ln in text.splitlines() if ln.startswith('diann-')), '')
    return runs, cmd


def run_design(report_runs: list, log_runs: list, raw: dict) -> pd.DataFrame:
    if sorted(report_runs) != sorted(log_runs):
        raise SystemExit('the runs in report.tsv differ from the runs the DIA-NN log loads')
    rows = []
    for run in report_runs:
        for date, pat in zip(BATCHES, RUN_PATTERNS):
            m = pat.match(run)
            if m:
                s = int(m.group(1))
                rows.append({'run': run, 'acquisition_date': '20' + date, 'sample': s,
                             'condition': 'minus_4SU' if s <= 5 else 'plus_4SU', 'raw_files_of_sample': '; '.join(raw[s])})
                break
        else:
            raise SystemExit(f'run {run!r} matches none of the three acquisition series')
    d = pd.DataFrame(rows)
    if len(d) != 30 or d.groupby(['acquisition_date', 'sample']).size().max() != 1 or set(d['sample']) != set(range(1, 11)):
        raise SystemExit('expected 30 runs: samples 1-10 acquired on each of three dates')
    return d.sort_values(['sample', 'acquisition_date']).reset_index(drop=True)


# ----------------------------------------------------------------------------------------------------------------
# Peptide intensities and statistics
# ----------------------------------------------------------------------------------------------------------------
def peptide_matrix(cache: Path, runs: list, qcut: float):
    """Per run x stripped sequence: sum of Precursor.Quantity over the precursors that pass the filter."""
    prec = pd.read_csv(cache / 'report_precursors.tsv.gz', sep='\t', dtype=str, keep_default_na=False)
    rr = pd.read_csv(cache / 'report_precursor_runs.tsv.gz', sep='\t',
                     usecols=['Run', 'Precursor.Id', 'Precursor.Quantity', 'Q.Value', 'Lib.Q.Value', 'Lib.PG.Q.Value'],
                     dtype={'Run': str, 'Precursor.Id': str})
    n_rows = len(rr)
    keep = ((rr['Q.Value'] <= qcut) & (rr['Lib.Q.Value'] <= qcut) & (rr['Lib.PG.Q.Value'] <= qcut)
            & (rr['Precursor.Quantity'] > 0))
    rr = rr.loc[keep]
    rr = rr.assign(peptide=rr['Precursor.Id'].map(prec.set_index('Precursor.Id')['Stripped.Sequence']))
    mat = rr.groupby(['peptide', 'Run'])['Precursor.Quantity'].sum().unstack('Run').reindex(columns=runs)
    used = prec[prec['Precursor.Id'].isin(set(rr['Precursor.Id']))]

    def join(v):
        return ';'.join(sorted({x for s in v for x in s.split(';') if x}))

    ann = used.groupby('Stripped.Sequence').agg(protein_ids=('Protein.Ids', join), genes=('Genes', join),
                                                precursors=('Precursor.Id', 'size')).reindex(mat.index)
    info = {'report_rows': n_rows, 'rows_passing_filter': int(keep.sum()), 'precursors_passing': int(len(used)),
            'peptides': int(len(mat))}
    return mat, ann, info


def _ttest(A: np.ndarray, B: np.ndarray, min_n: int) -> np.ndarray:
    """Student's two-sided t test (pooled variance) row by row over the values present (NaN = missing); NaN where a
    group has fewer than min_n values or both groups have no spread."""
    nA, nB = (~np.isnan(A)).sum(1), (~np.isnan(B)).sum(1)
    with warnings.catch_warnings(), np.errstate(all='ignore'):
        warnings.simplefilter('ignore', RuntimeWarning)
        mA, mB = np.nanmean(A, 1), np.nanmean(B, 1)
        vA, vB = np.nanvar(A, 1, ddof=1), np.nanvar(B, 1, ddof=1)
        df = nA + nB - 2
        sp2 = ((nA - 1) * vA + (nB - 1) * vB) / df
        t = (mA - mB) / np.sqrt(sp2 * (1.0 / nA + 1.0 / nB))
    ok = (nA >= min_n) & (nB >= min_n) & (sp2 > 0)
    p = np.full(len(A), np.nan)
    p[ok] = 2 * stats.t.sf(np.abs(t[ok]), df[ok])
    return p


def bh(p: np.ndarray) -> np.ndarray:
    q = np.full(len(p), np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    if not len(pv):
        return q
    order = np.argsort(pv)
    ranked = pv[order] * len(pv) / np.arange(1, len(pv) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(len(pv))
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def _nanmean_rows(X):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return np.nanmean(X, 1)


def peptide_stats(norm: np.ndarray, plus: np.ndarray, sample: np.ndarray, min_runs: int, min_samples: int, scale: str,
                  full: bool = True) -> dict:
    """norm: peptides x runs of median-normalised intensities (NaN = not quantified); plus: +4SU runs; sample: the
    sample of each run. Returns the per-peptide columns (He et al. 2016 ratio and score, Oksuz et al. 2023 tests)."""
    P, M = norm[:, plus], norm[:, ~plus]
    n_p, n_m = (~np.isnan(P)).sum(1), (~np.isnan(M)).sum(1)
    mean_p, mean_m = _nanmean_rows(P), _nanmean_rows(M)
    with np.errstate(all='ignore'):
        log2r = np.log2(mean_p / mean_m)
        tv = np.log2(norm) if scale == 'log2' else norm
    p = _ttest(tv[:, plus], tv[:, ~plus], min_runs)
    test = np.where(np.isnan(p), '', 't').astype(object)
    with np.errstate(all='ignore'):
        score = np.where(test == 't', -log2r * np.log10(p) ** 2, np.nan)
    out = {'runs_minus_4su': n_m, 'runs_plus_4su': n_p, 'mean_minus_4su': mean_m, 'mean_plus_4su': mean_p,
           'log2_ratio': np.where(test == 't', log2r, np.nan), 'rbrid_score': score}
    if not full:
        return out
    NP, NM = int(plus.sum()), int((~plus).sum())
    cache = {}
    for k in np.flatnonzero((n_p == 0) ^ (n_m == 0)):
        key = (int(n_p[k]), int(n_m[k]))
        if key not in cache:
            cache[key] = stats.fisher_exact([[key[0], NP - key[0]], [key[1], NM - key[1]]])[1]
        p[k], test[k] = cache[key], 'fisher'
    q = bh(p)
    positive = (p < 0.05) & (((test == 't') & (log2r < 0)) | ((test == 'fisher') & (n_p == 0)))
    # 5 x 5: technical replicates averaged per sample first
    samples = sorted(set(sample.tolist()))
    S = np.column_stack([_nanmean_rows(norm[:, sample == s]) for s in samples])
    splus = np.array([bool(plus[sample == s][0]) for s in samples])
    with np.errstate(all='ignore'):
        Sv = np.log2(S) if scale == 'log2' else S
        p5 = _ttest(Sv[:, splus], Sv[:, ~splus], min_samples)
        log2r5 = np.log2(_nanmean_rows(S[:, splus]) / _nanmean_rows(S[:, ~splus]))
        score5 = np.where(np.isnan(p5), np.nan, -log2r5 * np.log10(p5) ** 2)
    out.update({'test': test, 'p': p, 'q': q, 'rbrid_positive': positive,
                'samples_minus_4su': (~np.isnan(S[:, ~splus])).sum(1), 'samples_plus_4su': (~np.isnan(S[:, splus])).sum(1),
                'log2_ratio_5v5': np.where(np.isnan(p5), np.nan, log2r5), 'p_5v5': p5, 'rbrid_score_5v5': score5})
    return out


# ----------------------------------------------------------------------------------------------------------------
# Placing peptides on the reviewed canonical sequences
# ----------------------------------------------------------------------------------------------------------------
class Reference:
    def __init__(self, fasta: Path, entries: Path, history: Path):
        seq, acc = defaultdict(list), None
        with open(fasta, encoding='utf-8') as fh:
            for line in fh:
                if line.startswith('>'):
                    acc = line.split('|')[1]
                else:
                    seq[acc].append(line.strip())
        self.seq = {a: ''.join(v) for a, v in seq.items()}
        self.gene_of, self.gene = {}, defaultdict(set)
        with open(entries, encoding='utf-8') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                names = (r.get('Gene Names (primary)') or '').replace(';', ' ').split()
                self.gene_of[r['Entry']] = '; '.join(names)
                for g in names:
                    self.gene[g.upper()].add(r['Entry'])
        self.secondary, self.synonym = defaultdict(set), defaultdict(set)
        with open(history, encoding='utf-8') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                if r['kind'] == 'secondary':
                    self.secondary[r['key']].add(r['accession'])
                elif r['kind'] == 'synonym':
                    self.synonym[r['key'].upper()].add(r['accession'])


ROUTES = ['accession', 'secondary_accession', 'gene', 'gene_synonym']
DROP_REASONS = {'repeated': 'occurs more than once in the canonical sequence',
                'absent': 'not in the canonical sequence of the given accessions or genes',
                'no_entry': 'no reviewed reference entry for the given accessions or genes'}


def place_peptides(ann: pd.DataFrame, ref: Reference) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Each peptide goes to one reviewed canonical sequence that contains it exactly once. Candidates, in order:
    the accessions DIA-NN gives (isoform suffix removed); their current entry when they are secondary accessions;
    the reference entries of the genes DIA-NN gives; of those genes as synonyms. The first route with a match wins.
    Several matches: as Oksuz et al. 2023, the entry with most matching peptides in the dataset (then alphabetical);
    the others go to also_in."""
    placed, dropped = [], []
    for pep, pids, genes in zip(ann.index, ann['protein_ids'], ann['genes']):
        given = list(dict.fromkeys(x.split('-')[0] for x in str(pids).split(';') if x))
        gnames = [g.upper() for g in str(genes).split(';') if g]
        cand = {'accession': [a for a in given if a in ref.seq],
                'secondary_accession': sorted({c for a in given if a not in ref.seq for c in ref.secondary.get(a, ())}),
                'gene': sorted({c for g in gnames for c in ref.gene.get(g, ())}),
                'gene_synonym': sorted({c for g in gnames for c in ref.synonym.get(g, ())})}
        any_cand, repeated, hit = False, False, None
        for route in ROUTES:
            hits = []
            for a in cand[route]:
                any_cand = True
                s = ref.seq[a]
                i = s.find(pep)
                if i < 0:
                    continue
                if s.find(pep, i + 1) >= 0:
                    repeated = True
                    continue
                hits.append((a, i + 1))
            if hits:
                hit = (route, hits)
                break
        if hit is None:
            dropped.append({'peptide': pep, 'protein_ids': pids, 'genes': genes,
                            'reason': DROP_REASONS['repeated' if repeated else 'absent' if any_cand else 'no_entry']})
        else:
            placed.append((pep, hit[0], hit[1]))
    count = Counter(a for _, _, hits in placed for a, _ in hits)
    rows = []
    for pep, route, hits in placed:
        best = min(hits, key=lambda h: (-count[h[0]], h[0]))
        rows.append({'peptide': pep, 'uniprot_accession': best[0], 'start': best[1], 'end': best[1] + len(pep) - 1,
                     'mapped_by': route, 'also_in': '; '.join(sorted(a for a, _ in hits if a != best[0]))})
    return pd.DataFrame(rows), pd.DataFrame(dropped, columns=['peptide', 'protein_ids', 'genes', 'reason'])


# ----------------------------------------------------------------------------------------------------------------
# Residue scores (He et al. 2016): summed peptide scores, smoothed with Friedman's super smoother
# ----------------------------------------------------------------------------------------------------------------
SPANS = (0.05, 0.2, 0.5)  # tweeter, midrange, woofer: R's supsmu


def _running_lines(x, y, span, vsmlsq, cv):
    """Friedman's running-lines smoother as R's stats ppr.f 'smooth' (unit weights, sorted distinct x, not periodic):
    a least-squares line through the 2*ibw+1 nearest points (window clamped at the ends); with cv, also each point's
    cross-validated absolute residual."""
    n = len(x)
    ibw = max(2, int(0.5 * span * n + 0.5))
    it = min(2 * ibw + 1, n)
    lo = np.clip(np.arange(n) - ibw, 0, n - it)
    hi = lo + it
    z = np.zeros(1)
    cx, cy = np.concatenate([z, np.cumsum(x)]), np.concatenate([z, np.cumsum(y)])
    cxx, cxy = np.concatenate([z, np.cumsum(x * x)]), np.concatenate([z, np.cumsum(x * y)])
    xm, ym = (cx[hi] - cx[lo]) / it, (cy[hi] - cy[lo]) / it
    var = (cxx[hi] - cxx[lo]) - it * xm * xm
    cvar = (cxy[hi] - cxy[lo]) - it * xm * ym
    ok = var > vsmlsq
    safe = np.where(ok, var, 1.0)
    smo = np.where(ok, cvar / safe, 0.0) * (x - xm) + ym
    if not cv:
        return smo, None
    a = 1.0 - (1.0 / it + np.where(ok, (x - xm) ** 2 / safe, 0.0))
    acvr = np.zeros(n)
    good = a > 0
    acvr[good] = np.abs(y[good] - smo[good]) / a[good]
    for j in np.flatnonzero(~good):  # as the Fortran: the previous point's value
        acvr[j] = acvr[j - 1] if j > 0 else 0.0
    return smo, acvr


def supsmu(x, y) -> np.ndarray:
    """R's stats::supsmu(x, y) (span = "cv", bass = 0, unit weights) for sorted distinct x. Checked against R 4.5.3 on
    18 profiles: equal within 1e-8, except where the cross-validated residuals of two spans tie exactly (flat stretches
    of a summed profile), where the span R picks follows rounding noise (largest difference 0.002 on a value of 3.8)."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    n = len(x)
    if n == 0 or x[-1] <= x[0]:
        return np.full(n, y.mean() if n else np.nan)
    i, j = max(n // 4, 1), max(3 * (n // 4), 1)
    scale = x[j - 1] - x[i - 1]
    while scale <= 0:
        j, i = (j + 1 if j < n else j), (i - 1 if i > 1 else i)
        scale = x[j - 1] - x[i - 1]
    vsmlsq = (1e-3 * scale) ** 2
    fits, resid = [], []
    for s in SPANS:
        sm, cvr = _running_lines(x, y, s, vsmlsq, True)
        fits.append(sm)
        resid.append(_running_lines(x, cvr, SPANS[1], vsmlsq, False)[0])
    best = np.asarray(SPANS)[np.argmin(np.vstack(resid), axis=0)]
    sp = np.clip(_running_lines(x, best, SPANS[1], vsmlsq, False)[0], SPANS[0], SPANS[2])
    f = sp - SPANS[1]
    up, dn = f / (SPANS[2] - SPANS[1]), -f / (SPANS[1] - SPANS[0])
    mix = np.where(f >= 0, (1 - up) * fits[1] + up * fits[2], (1 - dn) * fits[1] + dn * fits[0])
    return _running_lines(x, mix, SPANS[0], vsmlsq, False)[0]


def residue_profiles(acc: np.ndarray, start: np.ndarray, end: np.ndarray, score: np.ndarray, lengths: dict) -> dict:
    """Per protein: covered positions (1-based), summed peptide scores and the smoothed residue score. Only residues
    under at least one scored peptide get a value; the smoother runs over those positions."""
    ok = np.isfinite(score)
    acc, start, end, score = acc[ok], start[ok], end[ok], score[ok]
    order = np.argsort(acc, kind='stable')
    acc, start, end, score = acc[order], start[order], end[order], score[order]
    bounds = np.flatnonzero(np.r_[True, acc[1:] != acc[:-1], True])
    out = {}
    for a, b in zip(bounds[:-1], bounds[1:]):
        L = lengths[acc[a]]
        add, cov = np.zeros(L + 2), np.zeros(L + 2)
        np.add.at(add, start[a:b], score[a:b])
        np.add.at(add, end[a:b] + 1, -score[a:b])
        np.add.at(cov, start[a:b], 1)
        np.add.at(cov, end[a:b] + 1, -1)
        covered = np.cumsum(cov)[1:L + 1] > 0
        pos = np.flatnonzero(covered) + 1
        summed = np.cumsum(add)[1:L + 1][covered]
        out[acc[a]] = (pos, summed, supsmu(pos, summed))
    return out


def find_segments(profiles: dict, threshold: float) -> list:
    """Maximal runs of consecutive covered residues with smoothed score >= threshold: (acc, start, end, indices)."""
    segs = []
    for acc, (pos, summed, sm) in profiles.items():
        idx = np.flatnonzero(sm >= threshold)
        if not len(idx):
            continue
        cut = np.flatnonzero((np.diff(idx) != 1) | (np.diff(pos[idx]) != 1)) + 1
        for g in np.split(idx, cut):
            segs.append((acc, int(pos[g[0]]), int(pos[g[-1]]), g))
    return segs


def segment_table(segs: list, profiles: dict, table: pd.DataFrame, ref: Reference) -> pd.DataFrame:
    scored = table[np.isfinite(table['rbrid_score'])]
    by_acc = {a: g for a, g in scored.groupby('uniprot_accession')}
    rows = []
    for acc, s, e, idx in segs:
        pos, summed, sm = profiles[acc]
        g = by_acc[acc]
        over = g[(g['start'] <= e) & (g['end'] >= s)]
        rows.append({'uniprot_accession': acc, 'gene': ref.gene_of.get(acc, ''), 'start': s, 'end': e, 'length': e - s + 1,
                     'max_residue_score': sm[idx].max(), 'mean_residue_score': sm[idx].mean(),
                     'max_summed_score': summed[idx].max(), 'peptides': len(over),
                     'rbrid_positive_peptides': int(over['rbrid_positive'].sum()), 'best_p': over['p'].min(),
                     'protein_length': len(ref.seq[acc]), 'protein_scored_residues': len(pos)})
    return pd.DataFrame(rows).sort_values(['uniprot_accession', 'start']).reset_index(drop=True)


def positive_residues(table: pd.DataFrame, profiles: dict, ref: Reference) -> pd.DataFrame:
    """Residues that pass the published cut-off: under at least one RBR-ID+ peptide (POSITIVE_RULE). Each keeps its
    residue scores (summed and smoothed, from all scored peptides over it; empty when only a Fisher-tested peptide
    covers it). No residue outside this rule is staged."""
    rows = []
    for acc, g in table[table['rbrid_positive'].astype(bool)].groupby('uniprot_accession'):
        seq = ref.seq[acc]
        L = len(seq)
        cnt, best = np.zeros(L + 2, dtype=int), np.full(L + 2, np.inf)
        for s, e, p in zip(g['start'], g['end'], g['p']):
            cnt[s:e + 1] += 1
            best[s:e + 1] = np.minimum(best[s:e + 1], p)
        prof = profiles.get(acc)
        sc = dict(zip(prof[0].tolist(), zip(prof[1].tolist(), prof[2].tolist()))) if prof else {}
        for i in np.flatnonzero(cnt[1:L + 1] > 0) + 1:
            summed, smooth = sc.get(int(i), (np.nan, np.nan))
            rows.append({'uniprot_accession': acc, 'gene': ref.gene_of.get(acc, ''), 'position': int(i),
                         'residue': seq[i - 1], 'rbrid_positive_peptides': int(cnt[i]), 'best_p': best[i],
                         'summed_score': summed, 'residue_score': smooth})
    return pd.DataFrame(rows)


def coverage_intervals(table: pd.DataFrame, ref: Reference) -> pd.DataFrame:
    """The measured background, per protein, as merged intervals (overlapping or adjacent peptides joined):
    quantified = under any staged peptide (quantified in at least one run after the q-value filter);
    tested     = under a peptide with a t or Fisher test, where an RBR-ID+ call was possible;
    scored     = under a peptide with an RBR-ID score (t test), where a residue score exists."""
    rows = []
    for kind, mask in [('quantified', np.ones(len(table), bool)), ('tested', table['test'].isin(['t', 'fisher']).to_numpy()),
                       ('scored', np.isfinite(table['rbrid_score'].to_numpy(dtype=float)))]:
        for acc, g in table[mask].groupby('uniprot_accession'):
            iv = sorted(zip(g['start'], g['end']))
            cs, ce = iv[0]
            for s, e in iv[1:]:
                if s <= ce + 1:
                    ce = max(ce, e)
                else:
                    rows.append((acc, cs, ce, kind))
                    cs, ce = s, e
            rows.append((acc, cs, ce, kind))
    df = pd.DataFrame(rows, columns=['uniprot_accession', 'start', 'end', 'coverage'])
    df.insert(1, 'gene', df['uniprot_accession'].map(ref.gene_of))
    df['length'] = df['end'] - df['start'] + 1
    return df.sort_values(['uniprot_accession', 'coverage', 'start']).reset_index(drop=True)


# ----------------------------------------------------------------------------------------------------------------
# Checks of the design from the data, and a label-permutation null for the residue threshold
# ----------------------------------------------------------------------------------------------------------------
def design_checks(norm: np.ndarray, design: pd.DataFrame, st: dict, table: pd.DataFrame, interpro: Path | None) -> list:
    out = []
    full = ~np.isnan(norm).any(1)
    C = np.corrcoef(np.log2(norm[full]).T)
    np.fill_diagonal(C, -np.inf)
    smp, date = design['sample'].to_numpy(), design['acquisition_date'].to_numpy()
    same_sample = sum(bool(np.all(smp[np.argsort(C[i])[::-1][:2]] == smp[i])) for i in range(len(smp)))
    same_date = sum(bool(date[np.argmax(C[i])] == date[i]) for i in range(len(smp)))
    cond = design['condition'].to_numpy()
    same_cond = sum(bool(cond[np.argmax(C[i])] == cond[i]) for i in range(len(smp)))
    out.append(f'runs whose two most correlated runs (log2 intensities of the {int(full.sum()):,} peptides quantified in '
               f'all 30) are the other two runs of the same sample: {same_sample} of 30; whose most correlated run has the '
               f'same acquisition date: {same_date} of 30; the same condition: {same_cond} of 30')
    sig = (st['test'] == 't') & (st['p'] < 0.05)
    dep = int((st['log2_ratio'][sig] < 0).sum())
    out.append(f't-tested peptides with p < 0.05: {int(sig.sum()):,}; depleted in +4SU {dep:,}, enriched {int(sig.sum()) - dep:,}')
    if interpro and interpro.exists():
        iv = pd.read_csv(interpro, sep='\t', usecols=['canonical_UniProtAC', 'interpro_accession', 'fragment_start', 'fragment_end'])
        rrm = iv[iv['interpro_accession'] == 'IPR000504']
        spans = defaultdict(list)
        for a, s, e in zip(rrm['canonical_UniProtAC'], rrm['fragment_start'], rrm['fragment_end']):
            spans[a].append((int(s), int(e)))
        t = table[np.isfinite(table['rbrid_score'])]
        in_rrm = np.array([any(s <= pe and e >= ps for s, e in spans.get(a, ())) for a, ps, pe in
                           zip(t['uniprot_accession'], t['start'], t['end'])])
        prot = t['uniprot_accession'].isin(spans).to_numpy()
        for name, m in [('overlapping an RRM (IPR000504)', in_rrm), ('elsewhere in RRM proteins', prot & ~in_rrm),
                        ('in other proteins', ~prot)]:
            sub = t[m]
            out.append(f'scored peptides {name}: {len(sub):,}, RBR-ID+ {sub["rbrid_positive"].mean():.1%}, '
                       f'median score {sub["rbrid_score"].median():.3f}')
    return out


def permutation_null(norm, design, row_of, table, lengths, args, observed: dict, grid) -> pd.DataFrame:
    """Balanced label permutations at the level of biological samples (technical replicates kept together): 2 or 3
    of the 5 samples swapped between the groups (200 such splits), --null-permutations of them at random. Each is
    scored exactly as the real split, and the fraction of scored residues at or above each threshold is compared."""
    import itertools
    sample = design['sample'].to_numpy()
    minus, plus = set(range(1, 6)), set(range(6, 11))
    splits = [set(c) for c in itertools.combinations(range(1, 11), 5) if len(set(c) & minus) in (2, 3)]
    rng = np.random.default_rng(args.seed)
    picks = rng.choice(len(splits), size=min(args.null_permutations, len(splits)), replace=False)
    obs = np.concatenate([v[2] for v in observed.values()])
    res = {T: [] for T in grid}
    segs_null = []
    acc, st_, en = table['uniprot_accession'].to_numpy(), table['start'].to_numpy(), table['end'].to_numpy()
    for k, i in enumerate(picks):
        plus_p = np.isin(sample, sorted(splits[i]))
        sc = peptide_stats(norm, plus_p, sample, args.min_runs, args.min_samples, args.scale, full=False)['rbrid_score'][row_of]
        prof = residue_profiles(acc, st_, en, sc, lengths)
        vals = np.concatenate([v[2] for v in prof.values()])
        for T in grid:
            res[T].append(float((vals >= T).mean()))
        segs_null.append(len(find_segments(prof, args.threshold)))
    rows = []
    for T in grid:
        o, nl = float((obs >= T).mean()), float(np.mean(res[T]))
        rows.append({'threshold': T, 'observed_fraction': o, 'null_mean_fraction': nl, 'null_max_fraction': max(res[T]),
                     'estimated_fdr': min(1.0, nl / o) if o else float('nan')})
    out = pd.DataFrame(rows)
    out.attrs['segments_null'] = segs_null
    out.attrs['splits'] = [sorted(splits[i]) for i in picks]
    return out


# ----------------------------------------------------------------------------------------------------------------
# Writing
# ----------------------------------------------------------------------------------------------------------------
def fmt(v):
    if isinstance(v, (bool, np.bool_)):
        return 'true' if v else 'false'
    if isinstance(v, (float, np.floating)):
        return '' if not np.isfinite(v) else f'{v:.4g}'
    return '' if v is None else str(v)


def write_tsv(df: pd.DataFrame, path: Path) -> int:
    with open(path, 'w', encoding='utf-8', newline='') as fh:
        fh.write('\t'.join(df.columns) + '\n')
        for row in df.itertuples(index=False):
            fh.write('\t'.join(fmt(v) for v in row) + '\n')
    return len(df)


def update_manifest(manifest: Path, rows: list):
    """Replace the rows of these files, or add them after the last interactions/ row."""
    lines = manifest.read_text(encoding='utf-8').splitlines()
    head = lines[0].split('\t')
    new = {r['file']: '\t'.join(str(r.get(c, '')) for c in head) for r in rows}
    body, done = [], set()
    for ln in lines[1:]:
        f = ln.split('\t')[0]
        if f in new:
            body.append(new[f])
            done.add(f)
        else:
            body.append(ln)
    rest = [new[r['file']] for r in rows if r['file'] not in done]
    if rest:
        last = max((i for i, ln in enumerate(body) if ln.startswith('interactions/')), default=len(body) - 1)
        body[last + 1:last + 1] = rest
    manifest.write_text('\n'.join([lines[0]] + body) + '\n', encoding='utf-8')


def main():
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cache', type=Path, required=True, help='download folder (PRIDE metadata, log, report aggregate)')
    ap.add_argument('--out', type=Path, default=root / 'inputs' / 'interactions')
    ap.add_argument('--fasta', type=Path, default=root / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03.fasta')
    ap.add_argument('--entries', type=Path, default=root / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03_entries.tsv')
    ap.add_argument('--history', type=Path, default=root / 'inputs' / 'reference' / 'uniprot_accession_history_2026_03.tsv')
    ap.add_argument('--interpro', type=Path, default=root / 'inputs' / 'annotation' / 'interpro_108.0_domain_intervals.tsv',
                    help='only for the RRM check printed at the end')
    ap.add_argument('--manifest', default=str(root / 'inputs' / 'MANIFEST.tsv'), help="'' to leave the manifest alone")
    ap.add_argument('--qvalue', type=float, default=0.01, help='Q.Value, Lib.Q.Value and Lib.PG.Q.Value cut-off')
    ap.add_argument('--min-runs', type=int, default=3, help='quantified runs per condition for the t test (of 15)')
    ap.add_argument('--min-samples', type=int, default=3, help='quantified samples per condition for the 5 x 5 test')
    ap.add_argument('--scale', choices=['linear', 'log2'], default='linear', help='scale of the values the t test compares')
    ap.add_argument('--threshold', type=float, default=3.0, help='smoothed residue score that makes an RBR residue')
    ap.add_argument('--null-permutations', type=int, default=20)
    ap.add_argument('--seed', type=int, default=1)
    args = ap.parse_args()
    cache, out = args.cache, args.out
    cache.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)

    proj, files = pride_metadata(cache)
    if 'CC0' not in proj.get('license', ''):
        raise SystemExit(f'PRIDE licence is {proj.get("license")!r}, not CC0: check before staging')
    rep = next(f for f in files if f['fileName'] == REPORT)
    log_path = fetch(f'{PRIDE_FTP}/{LOG}', cache / LOG)
    facts = aggregate_report(cache, f'{PRIDE_FTP}/{REPORT}', int(rep['fileSizeBytes']), rep.get('checksum'))
    log_runs, cmd = search_log(log_path)
    design = run_design(list(facts['runs']), log_runs, raw_files_by_sample(files))
    log('DIA-NN command line:', cmd)
    log('UniMod ids in Modified.Sequence (rows):', facts['unimod_counts'])
    methyl = {u: n for u, n in facts['unimod_counts'].items() if u in ('34', '36', '37')}
    log('  methyl (UniMod:34, 36, 37):', methyl or 'none')

    mat, ann, info = peptide_matrix(cache, list(design['run']), args.qvalue)
    X = mat.to_numpy(dtype=float)
    med = np.nanmedian(X, 0)
    norm = X / med
    design['peptides_quantified'] = (~np.isnan(X)).sum(0)
    design['median_peptide_intensity'] = med
    plus = (design['condition'] == 'plus_4SU').to_numpy()
    sample = design['sample'].to_numpy()
    st = peptide_stats(norm, plus, sample, args.min_runs, args.min_samples, args.scale)
    pep_all = pd.DataFrame(st, index=mat.index)
    log(f'{info["report_rows"]:,} report rows; {info["rows_passing_filter"]:,} pass the q-value filter; '
        f'{info["precursors_passing"]:,} precursors, {info["peptides"]:,} peptides')
    log('tests:', dict(Counter(pep_all['test'].replace('', 'none'))))

    ref = Reference(args.fasta, args.entries, args.history)
    placed, dropped = place_peptides(ann, ref)
    table = placed.merge(ann.rename(columns={'protein_ids': 'diann_protein_ids', 'genes': 'diann_genes'}),
                         left_on='peptide', right_index=True).merge(pep_all, left_on='peptide', right_index=True)
    row_of = mat.index.get_indexer(table['peptide'])
    table.insert(1, 'gene', table['uniprot_accession'].map(ref.gene_of))
    lengths = {a: len(s) for a, s in ref.seq.items()}
    log(f'placed {len(table):,} peptides on {table["uniprot_accession"].nunique():,} entries:',
        dict(Counter(table['mapped_by'])), f'| dropped {len(dropped):,}:', dict(Counter(dropped['reason'])))

    # residue scores come from the peptide scores as written (4 significant digits), so the staged table reproduces them
    table['rbrid_score'] = [float(fmt(v)) if np.isfinite(v) else np.nan for v in table['rbrid_score']]
    profiles = residue_profiles(table['uniprot_accession'].to_numpy(), table['start'].to_numpy(), table['end'].to_numpy(),
                                table['rbrid_score'].to_numpy(), lengths)
    segs = find_segments(profiles, args.threshold)
    seg = segment_table(segs, profiles, table, ref)
    scored_res = sum(len(v[0]) for v in profiles.values())
    log(f'{len(profiles):,} proteins with scored residues ({scored_res:,} residues); {len(seg):,} segments >= '
        f'{args.threshold} in {seg["uniprot_accession"].nunique():,} proteins ({int(seg["length"].sum()):,} residues)')
    res = positive_residues(table, profiles, ref)
    res_regions = int(((res['uniprot_accession'] != res['uniprot_accession'].shift())
                       | (res['position'] != res['position'].shift() + 1)).sum())
    log(f'published cut-off ({POSITIVE_RULE.split(":")[0]}): {len(res):,} residues in {res_regions:,} runs of consecutive '
        f'residues, {res["uniprot_accession"].nunique():,} proteins; {int(res["residue_score"].notna().sum()):,} with a '
        f'residue score; {int((res["residue_score"] >= args.threshold).sum()):,} also >= {args.threshold}')
    cov = coverage_intervals(table, ref)
    for kind, g in cov.groupby('coverage'):
        log(f'coverage {kind}: {len(g):,} intervals, {int(g["length"].sum()):,} residues, {g["uniprot_accession"].nunique():,} proteins')
    checks = design_checks(norm, design, st, table, args.interpro)
    for c in checks:
        log('check:', c)
    null = None
    if args.null_permutations > 0:
        grid = sorted({0.5, 1.0, 2.0, 3.0, 5.0, args.threshold})
        null = permutation_null(norm, design, row_of, table, lengths, args, profiles, grid)
        log(null.to_string(index=False))
        log(f'segments at {args.threshold} under the null: {null.attrs["segments_null"]}')

    # ---- tables
    cols = ['uniprot_accession', 'gene', 'start', 'end', 'peptide', 'mapped_by', 'also_in', 'diann_protein_ids',
            'diann_genes', 'precursors', 'runs_minus_4su', 'runs_plus_4su', 'mean_minus_4su', 'mean_plus_4su', 'log2_ratio',
            'test', 'p', 'q', 'rbrid_score', 'rbrid_positive', 'samples_minus_4su', 'samples_plus_4su', 'log2_ratio_5v5',
            'p_5v5', 'rbrid_score_5v5']
    table = table[cols].sort_values(['uniprot_accession', 'start', 'end']).reset_index(drop=True)
    n_pep = write_tsv(table, out / PEPTIDES)
    n_seg = write_tsv(seg, out / SEGMENTS)
    n_res = write_tsv(res, out / RESIDUES)
    n_cov = write_tsv(cov, out / COVERAGE)
    n_run = write_tsv(design, out / RUNS)
    if len(dropped):
        write_tsv(dropped, cache / 'rbrid_dropped_peptides.tsv')
    raw_rows = [{'path': f'PXD035484/{REPORT}', 'url': f'{PRIDE_FTP}/{REPORT}', 'bytes': facts['bytes'], 'md5': facts['md5'],
                 'saved': facts['streamed_on'],
                 'note': f'streamed, not kept; SHA-1 {facts["sha1"]} equals the checksum PRIDE lists'}]
    for name, url in [(LOG, f'{PRIDE_FTP}/{LOG}'), ('pride_project.json', PRIDE_API),
                      ('pride_files.json', PRIDE_API + '/files?pageSize=500&page=0')]:
        p = cache / name
        raw_rows.append({'path': f'PXD035484/{name}', 'url': url, 'bytes': p.stat().st_size, 'md5': md5_of(p),
                         'saved': dt.date.fromtimestamp(p.stat().st_mtime).isoformat(), 'note': 'kept in the download folder'})
    write_tsv(pd.DataFrame(raw_rows), out / RAW_FILES)

    # ---- manifest rows
    if args.manifest:
        n_prot = table['uniprot_accession'].nunique()
        n_scored = int(np.isfinite(table['rbrid_score']).sum())
        how = (f'scripts/43c_stage_rbrid.py: report filtered at Q.Value, Lib.Q.Value and Lib.PG.Q.Value <= {args.qvalue}; '
               f'precursors summed per stripped sequence and run; median-normalised per run; Student t test on the '
               f'{args.scale} normalised values, 15 vs 15 runs, at least {args.min_runs} runs per condition; Fisher test '
               f'when absent from one condition; BH q; score = -log2(+4SU/-4SU) x (log10 p)^2 (He et al. 2016). Placed by '
               f'exact unique match to inputs/reference/uniprot_human_reviewed_2026_03.fasta via the accession given, '
               f'else the gene ({len(dropped):,} peptides dropped)')
        src = f'{PRIDE_FTP}/{REPORT}'
        ver = f'PRIDE PXD035484 (public 2023-07-02), streamed {facts["streamed_on"]}'
        rows = [
            {'file': f'{out.name}/{PEPTIDES}', 'rows': n_pep,
             'content': (f'RBR-ID (K562 nuclei, 4-thiouridine and 312 nm UV vs UV alone; DIA, DIA-NN 1.8.1): {n_pep:,} '
                         f'peptides on {n_prot:,} reviewed canonical entries with start and end, runs quantified per '
                         f'condition (of 15), mean normalised intensities, log2 +4SU/-4SU, p (t or Fisher), BH q, RBR-ID '
                         f'score ({n_scored:,} scored), RBR-ID+ call (p < 0.05, depleted), and the 5 x 5 test on sample '
                         f'means. The search did not allow methylarginine'),
             'derived_from': f'PXD035484/{REPORT}', 'how_made_or_checked': how, 'original_name': REPORT},
            {'file': f'{out.name}/{RESIDUES}', 'rows': n_res,
             'content': (f'RBR-ID residues at or above the published cut-off, {POSITIVE_RULE}: {n_res:,} residues '
                         f'({res_regions:,} runs of consecutive residues) in {res["uniprot_accession"].nunique():,} '
                         f'proteins, each with the number of RBR-ID+ peptides over it, their best p, and its residue scores '
                         f'(summed peptide scores and the super-smoothed score of He et al. 2016; empty where only a '
                         f'Fisher-tested peptide covers it). No other residue is staged'),
             'derived_from': f'{out.name}/{PEPTIDES}',
             'how_made_or_checked': 'scripts/43c_stage_rbrid.py positive_residues', 'original_name': ''},
            {'file': f'{out.name}/{COVERAGE}', 'rows': n_cov,
             'content': ('The measured background of RBR-ID per protein, as merged residue intervals: quantified (under any '
                         'staged peptide), tested (under a peptide with a t or Fisher test, where an RBR-ID+ call was '
                         'possible) and scored (under a peptide with an RBR-ID score). An arginine is measured when its '
                         'position lies in a quantified interval'),
             'derived_from': f'{out.name}/{PEPTIDES}',
             'how_made_or_checked': 'scripts/43c_stage_rbrid.py coverage_intervals (overlapping or adjacent peptides joined)',
             'original_name': ''},
            {'file': f'{out.name}/{SEGMENTS}', 'rows': n_seg,
             'content': (f'RBR-ID RNA-binding regions: runs of consecutive residues whose residue score (peptide scores '
                         f'summed per residue, smoothed with Friedman\'s super smoother; He et al. 2016) is >= '
                         f'{args.threshold}, with the maximum and mean residue score, the peptides over them and the best p. '
                         f'The threshold is this script\'s choice (no published residue cut-off): the lowest of 1, 2, 3, 5 '
                         f'with a label-permutation FDR estimate under 10%'),
             'derived_from': f'{out.name}/{PEPTIDES}',
             'how_made_or_checked': ('scripts/43c_stage_rbrid.py (residue_profiles, find_segments; the super smoother is a '
                                     'port of R stats::supsmu checked against R 4.5.3)'
                                     + (f'; label-permutation estimate of the residue FDR at {args.threshold}: '
                                        f'{null.loc[null.threshold == args.threshold, "estimated_fdr"].iloc[0]:.2f}'
                                        if null is not None else '')),
             'original_name': ''},
            {'file': f'{out.name}/{RUNS}', 'rows': n_run,
             'content': ('The 30 runs of report.tsv with acquisition date, sample (1-5 without 4SU, 6-10 with), condition, '
                         'the raw files of that sample, peptides quantified and the median peptide intensity used to '
                         'normalise'),
             'derived_from': f'PXD035484/{LOG}; PXD035484/pride_files.json',
             'how_made_or_checked': ('scripts/43c_stage_rbrid.py run_design: sample number in the run name = raw-file '
                                     'sample number; the three dates are the technical replicates'),
             'original_name': ''},
            {'file': f'{out.name}/{RAW_FILES}', 'rows': len(raw_rows),
             'content': 'The PRIDE files behind the rbrid_oksuz2023 tables, with their size, MD5 and date',
             'derived_from': '', 'how_made_or_checked': 'scripts/43c_stage_rbrid.py', 'original_name': ''},
        ]
        for r in rows:
            p = out / Path(r['file']).name
            r.update({'bytes': p.stat().st_size, 'md5': md5_of(p), 'redistribution': 'git', 'source': src,
                      'version_or_date': ver, 'citation': CITATION, 'licence': LICENCE})
        rows[-1].update({'source': f'{PRIDE_FTP}/', 'version_or_date': facts['streamed_on'], 'citation': '', 'licence': ''})
        update_manifest(Path(args.manifest), rows)
    log(f'wrote {out / PEPTIDES} ({n_pep:,}), {out / RESIDUES} ({n_res:,}), {out / COVERAGE} ({n_cov:,}), '
        f'{out / SEGMENTS} ({n_seg:,}), {out / RUNS} ({n_run}), {out / RAW_FILES}')


if __name__ == '__main__':
    main()



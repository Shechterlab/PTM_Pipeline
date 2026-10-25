#!/usr/bin/env python3
"""Stage the pCLAP RNA-binding peptides of Mullari et al. 2017 (HEK293) from the PRIDE MaxQuant search (inputs/interactions/).

    python -I scripts/43d_stage_pclap.py --cache PTM_data/external/pclap --out inputs/interactions [--preview [--preview-out FILE]]
    python -I scripts/43d_stage_pclap.py --no-stage --preview [--sites releases/v1.2.0/human_methylarginine_sites_v1.2.0.tsv]

pCLAP (peptide cross-linking and affinity purification; Mullari, Lyon, Jensen & Nielsen, J Proteome Res 16:2762, 2017):
cells are UV-crosslinked, lysed and digested with Lys-C, so only the Lys-C peptides crosslinked to RNA stay attached to it;
poly(A) RNA is pulled down with oligo(dT), eluted, digested with RNase and trypsin, and the tryptic peptides are measured
(label-free, single shots). A peptide counts as RNA-binding when it is enriched in UV-crosslinked (XL) over
non-crosslinked (No XL) pull-downs. The RNA-binding region it stands for is its Lys-C peptide.

Neither the article (ACS, not open access) nor PRIDE gives the list of RNA-binding peptides, so it is derived here from
the deposited MaxQuant output (PXD005739, MaxQuantOutput.tar.gz: 4 XL and 4 No XL runs of 17 August 2016):

  * reverse hits and contaminants removed; intensities log2-transformed (0 = missing);
  * kept for testing: at least 3 valid values out of 4 in at least one group;
  * missing values imputed per run from a normal distribution 1.8 standard deviations below the run's mean, 0.3
    standard deviations wide (Perseus defaults);
  * XL against No XL: two-sided Student's t-test, with s0 = 0.1 in the statistic and a permutation-based FDR of 5 %
    (all 70 splits of the eight runs into two groups of four); RNA-binding = FDR-significant with XL > No XL.

That is the rule Mullari et al., Nat Commun 14:4348 (2023) give for pCLAP's "Filter 1" ("as previously described", citing
the 2017 article); the 2017 article itself could not be read here. Imputation is random, so the test is run on
--imputations draws (seeded) and a peptide is called RNA-binding when it passes in at least half of them. As a check, the
call is compared at gene level with the 2023 article's Supplementary Data 1 (CC BY 4.0), which marks which of its 526 mouse
brain RBPs the 2017 HEK293 pCLAP found; that file is downloaded into the cache and not staged.

The peptide table carries the authors' MaxQuant intensities, and PRIDE gives the dataset only the EMBL-EBI terms of use
(no licence from the submitters), so it is not redistributed: MANIFEST redistribution = checksum_only, and
inputs/interactions/.gitignore keeps it out of Git. Run this script to rebuild it (same MD5).

Each peptide is placed on the canonical sequence of inputs/reference/uniprot_human_reviewed_2026_03.fasta by exact,
unique match: the leading razor protein, then the other accessions MaxQuant gives (a secondary accession followed to its
entry), then the gene names (primary names, then synonyms; one entry only); everything else is dropped and counted. The
Lys-C peptide is found as the authors' tool does it (github.com/dblyon/pCLAP_Tryp2LysC: the in-silico Lys-C peptide of
at least 7 residues, with at most 6 missed cleavages, that contains the tryptic peptide and has the fewest missed
cleavages).

The MaxQuant search allowed only Acetyl (Protein N-term) and Oxidation (M) as variable modifications: a peptide that
carries a methylarginine could be identified only in its unmethylated form.

--preview compares the methylarginines of a release with the other arginines of the same proteins (Mantel-Haenszel odds
ratios, one stratum per protein, among proteins with at least one RNA-binding peptide). It writes nothing in inputs/.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import io
import itertools
import json
import math
import re
import tarfile
import time
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

UA = 'PTM_Pipeline-stager/1.0 (https://github.com/Shechterlab/PTM_Pipeline)'
PRIDE_ACC = 'PXD005739'
PRIDE_API = f'https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{PRIDE_ACC}'
TARBALL = 'MaxQuantOutput.tar.gz'
TARBALL_URL = f'https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/06/{PRIDE_ACC}/{TARBALL}'
MEMBERS = ['MaxQuantOutput/peptides.txt', 'MaxQuantOutput/parameters.txt', 'MaxQuantOutput/summary.txt']
# Brain-pCLAP (2023, CC BY 4.0) Supplementary Data 1: which of its 526 mouse brain RBPs the 2017 HEK293 pCLAP found (check only)
BRAIN_SD1_URL = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-39936-x/MediaObjects/41467_2023_39936_MOESM3_ESM.xlsx'
BRAIN_SD1 = 'mullari2023_moesm3_supplementary_data_1.xlsx'
TABLE = 'pclap_mullari2017_peptides.tsv'
RAW_TABLE = 'raw_pclap_files.tsv'
CITATION = 'Mullari et al., J Proteome Res 16:2762 (2017), doi:10.1021/acs.jproteome.7b00042'
RULE_CITATION = 'Mullari et al., Nat Commun 14:4348 (2023), doi:10.1038/s41467-023-39936-x'
# Perseus defaults and the 2023 article's "Filter 1"
MIN_VALID = 3
IMPUTE_SHIFT, IMPUTE_WIDTH = 1.8, 0.3
S0, FDR = 0.1, 0.05
# The authors' tryptic-to-Lys-C mapping (pCLAP_Tryp2LysC, tryptic2LysC.py)
LYSC_MAX_MISSED, LYSC_MIN_LEN = 6, 7


def log(*a):
    print(*a, flush=True)


def md5_bytes(b: bytes) -> str:
    return hashlib.md5(b).hexdigest()


def file_hash(path: Path, algo: str = 'md5') -> str:
    h = hashlib.new(algo)
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path) -> Path:
    """Download once into the cache (kept as is on later runs)."""
    if dest.exists() and dest.stat().st_size:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + '.part')
    for k in range(5):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': UA, 'Accept': '*/*'})
            with urllib.request.urlopen(req, timeout=600) as r, open(part, 'wb') as fh:
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


# ---- inputs ---------------------------------------------------------------------------------------------------------------

def read_members(tarball: Path) -> dict[str, bytes]:
    """The MaxQuant tables the staging reads, straight from the archive (nothing is extracted to disk)."""
    out = {}
    with tarfile.open(tarball, 'r:gz') as t:
        for name in MEMBERS:
            m = t.getmember(name)
            if not m.isreg():
                raise SystemExit(f'{name}: not a regular file in {tarball.name}')
            out[name] = t.extractfile(m).read()
    return out


def parameters(raw: bytes) -> dict:
    p = {}
    for line in raw.decode('utf-8').splitlines()[1:]:
        k, _, v = line.partition('\t')
        p.setdefault(k, v)
    return p


def runs(raw: bytes) -> list[dict]:
    rows = list(csv.DictReader(io.StringIO(raw.decode('utf-8')), delimiter='\t'))
    return [r for r in rows if r.get('Enzyme')]  # the per-experiment and total rows have no search settings


class Reference:
    def __init__(self, fasta: Path, entries: Path, history: Path):
        self.seq, acc = {}, None
        for line in fasta.read_text(encoding='utf-8').splitlines():
            if line.startswith('>'):
                acc = line.split('|')[1]
                self.seq[acc] = []
            elif acc:
                self.seq[acc].append(line.strip())
        self.seq = {a: ''.join(s) for a, s in self.seq.items()}
        self.gene, self.by_gene = {}, defaultdict(set)
        with open(entries, encoding='utf-8') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                genes = (r.get('Gene Names (primary)') or '').replace(';', ' ').split()
                self.gene[r['Entry']] = ';'.join(genes)
                for g in genes:
                    self.by_gene[g.upper()].add(r['Entry'])
        self.secondary, self.synonym = defaultdict(set), defaultdict(set)
        with open(history, encoding='utf-8') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                if r['accession'] not in self.seq:
                    continue
                if r['kind'] == 'secondary':
                    self.secondary[r['key']].add(r['accession'])
                elif r['kind'] == 'synonym':
                    self.synonym[r['key'].upper()].add(r['accession'])

    def unique_start(self, acc: str, pep: str):
        """0-based start of the peptide in the canonical sequence if it occurs there exactly once, else None."""
        s = self.seq.get(acc, '')
        i = s.find(pep)
        if i < 0 or s.find(pep, i + 1) >= 0:
            return None
        return i

    def place(self, pep: str, razor: str, proteins: list[str], genes: list[str]):
        """(accession, 0-based start, how) or (None, None, reason)."""
        tried = set()
        for k, given in enumerate([razor] + [p for p in proteins if p != razor]):
            base = given.split('-')[0]
            accs, how = ({base}, 'leading_razor_accession' if k == 0 else 'other_given_accession') if base in self.seq \
                else (self.secondary.get(base, set()), 'secondary_accession')
            if len(accs) != 1:
                continue
            acc = next(iter(accs))
            tried.add(acc)
            i = self.unique_start(acc, pep)
            if i is not None:
                return acc, i, how
        hits, how = set(), ''
        for g in genes:
            accs = self.by_gene.get(g.upper()) or self.synonym.get(g.upper()) or set()
            tried |= accs
            for acc in accs:
                if self.unique_start(acc, pep) is not None:
                    hits.add(acc)
                    how = 'gene_primary' if self.by_gene.get(g.upper()) else 'gene_synonym'
        if len(hits) == 1:
            acc = next(iter(hits))
            return acc, self.unique_start(acc, pep), how
        if len(hits) > 1:
            return None, None, 'gene_matches_several_entries'
        if not tried:
            return None, None, 'no_reviewed_entry'
        if any(self.seq[a].count(pep) > 1 for a in tried):
            return None, None, 'more_than_one_position'
        return None, None, 'not_in_canonical_sequence'


def lysc_region(seq: str, start: int, end: int):
    """The Lys-C peptide holding seq[start:end] as pCLAP_Tryp2LysC chooses it: in-silico Lys-C peptides (cleavage after
    every K) of >= 7 residues with <= 6 missed cleavages, the one with the fewest missed cleavages (counted as the tool
    does: Ks in the peptide minus one), ties to the alphabetically first sequence. Returns 0-based (start, end, missed)."""
    cuts = sorted({0, len(seq), *[m.end() for m in re.finditer('K', seq)]})
    best = None
    for i, a in enumerate(cuts):
        if a > start:
            break
        for b in cuts[i + 1:i + LYSC_MAX_MISSED + 2]:
            if b < end or b - a < LYSC_MIN_LEN:
                continue
            frag = seq[a:b]
            key = (max(frag.count('K') - 1, 0), frag)
            if best is None or key < best[0]:
                best = (key, a, b)
    return (best[1], best[2], best[0][0]) if best else (None, None, None)


# ---- the XL against No XL test --------------------------------------------------------------------------------------------

def statistic(y: np.ndarray, g: np.ndarray):
    """Two-sample Student's t on rows of y (group g against the rest): difference, its standard error."""
    a, b = y[:, g], y[:, ~g]
    na, nb = a.shape[1], b.shape[1]
    sp2 = ((na - 1) * a.var(1, ddof=1) + (nb - 1) * b.var(1, ddof=1)) / (na + nb - 2)
    return a.mean(1) - b.mean(1), np.sqrt(sp2 * (1 / na + 1 / nb))


def permutation_q(d_obs: np.ndarray, d_null: np.ndarray, n_perm: int) -> np.ndarray:
    """Two-sided permutation FDR (expected null count at or above |d| over the observed count), as q-values."""
    obs = np.abs(d_obs)
    null = np.sort(np.abs(d_null).ravel())
    order = np.argsort(-obs, kind='stable')
    so = obs[order]
    exp_false = (null.size - np.searchsorted(null, so, side='left')) / n_perm
    called = np.searchsorted(-so, -so, side='right')  # rows with |d| >= this one
    fdr = np.minimum(exp_false / called, 1.0)
    q_sorted = np.minimum.accumulate(fdr[::-1])[::-1]
    q = np.empty_like(q_sorted)
    q[order] = q_sorted
    return q


def xl_test(log2: np.ndarray, is_xl: np.ndarray, n_imputations: int, seed: int):
    """Imputation, t-test and permutation FDR, repeated over seeded imputations. Returns per-row arrays."""
    n = log2.shape[0]
    splits = [np.isin(np.arange(log2.shape[1]), c) for c in itertools.combinations(range(log2.shape[1]), int(is_xl.sum()))]
    valid = ~np.isnan(log2)
    mu = np.nanmean(log2, 0)
    sd = np.nanstd(log2, 0, ddof=1)
    passes = np.zeros((n, n_imputations), bool)
    diffs, ps, qs = (np.zeros((n, n_imputations)) for _ in range(3))
    for k in range(n_imputations):
        rng = np.random.default_rng(seed + k)
        y = log2.copy()
        for j in range(y.shape[1]):
            miss = ~valid[:, j]
            y[miss, j] = rng.normal(mu[j] - IMPUTE_SHIFT * sd[j], IMPUTE_WIDTH * sd[j], miss.sum())
        diff, se = statistic(y, is_xl)
        d_obs = diff / (se + S0)
        d_null = []
        for g in splits:
            dg, sg = statistic(y, g)
            d_null.append(dg / (sg + S0))
        d_null = np.stack(d_null)
        q = permutation_q(d_obs, d_null, len(splits))
        t = diff / se
        p = 2 * stats.t.sf(np.abs(t), is_xl.sum() + (~is_xl).sum() - 2)
        diffs[:, k], ps[:, k], qs[:, k] = diff, p, q
        passes[:, k] = (q <= FDR) & (diff > 0)
    return {'diff': np.median(diffs, 1), 'p': np.median(ps, 1), 'q': np.median(qs, 1), 'pass_fraction': passes.mean(1),
            'passes_per_imputation': passes.sum(0)}


def check_against_2023(xlsx: Path, rows: list, human_genes: set) -> dict:
    """Protein-level check of the reproduced call: brain-pCLAP's Supplementary Data 1 marks which of its 526 brain RBPs
    (mouse gene symbols) the 2017 HEK293 pCLAP found ('HEK293-pCLAP', '+'). Compared through identical symbols."""
    raw = pd.read_excel(xlsx, sheet_name='Comparison with other studies', header=None)
    row = next(i for i in range(len(raw)) if 'HEK293-pCLAP' in [str(v).strip() for v in raw.iloc[i].tolist()])
    col = [str(v).strip() for v in raw.iloc[row].tolist()].index('HEK293-pCLAP')
    body = raw.iloc[row + 1:]
    genes, flag = body.iloc[:, 0].astype(str).str.strip(), body.iloc[:, col].astype(str).str.strip()
    ok = genes.ne('nan') & genes.ne('')
    genes, flag = genes[ok].str.upper(), flag[ok]
    bind = {g for r in rows if r[-1] == 'true' for g in r[1].upper().split(';') if g}
    meas = {g for r in rows for g in r[1].upper().split(';') if g}
    out = {'brain_rbps': int(len(genes))}
    for label, sel in (('flagged_2017', flag.eq('+')), ('not_flagged', ~flag.eq('+'))):
        gs = set(genes[sel]) & human_genes
        out[label] = {'same_human_symbol': len(gs), 'rna_binding_here': len(gs & bind),
                      'measured_not_rna_binding_here': len((gs & meas) - bind), 'not_measured_here': len(gs - meas)}
    return out


# ---- staging --------------------------------------------------------------------------------------------------------------

def stage(args) -> dict:
    cache = Path(args.cache)
    project_json = fetch(PRIDE_API, cache / 'pride_project.json')
    files_json = fetch(f'{PRIDE_API}/files?pageSize=100', cache / 'pride_files.json')
    tarball = fetch(TARBALL_URL, cache / TARBALL)
    project = json.loads(project_json.read_text(encoding='utf-8'))
    listed = {f['fileName']: f for f in json.loads(files_json.read_text(encoding='utf-8'))}
    sha1 = file_hash(tarball, 'sha1')
    if listed.get(TARBALL, {}).get('checksum') != sha1:
        raise SystemExit(f'{TARBALL}: SHA-1 {sha1} differs from PRIDE\'s {listed.get(TARBALL, {}).get("checksum")}')
    log(f'{TARBALL}: {tarball.stat().st_size:,} bytes, SHA-1 matches PRIDE; PRIDE licence: {project.get("license")!r}')

    members = read_members(tarball)
    par = parameters(members['MaxQuantOutput/parameters.txt'])
    rr = runs(members['MaxQuantOutput/summary.txt'])
    var_mods = sorted({m for r in rr for m in r['Variable modifications'].split(';') if m})
    methyl_r = [m for m in var_mods if re.search(r'methyl', m, re.I) and re.search(r'\([^)]*R[^)]*\)', m)]
    exps = [r['Experiment'] for r in rr]
    xl = sorted(e for e in exps if re.fullmatch(r'XL \d+', e))
    noxl = sorted(e for e in exps if re.fullmatch(r'No XL \d+', e))
    if len(xl) != 4 or len(noxl) != 4 or len(exps) != 8:
        raise SystemExit(f'unexpected experiments: {exps}')
    log(f'MaxQuant {par.get("Version")}; enzyme {sorted({r["Enzyme"] for r in rr})}; variable modifications {var_mods}; '
        f'fixed {par.get("Fixed modifications")}; match between runs {par.get("Match between runs")}; '
        f'methylarginine searched: {bool(methyl_r)}')

    pep = pd.read_csv(io.BytesIO(members['MaxQuantOutput/peptides.txt']), sep='\t', dtype=str, keep_default_na=False)
    n_all = len(pep)
    pep = pep[(pep['Reverse'] != '+') & (pep['Potential contaminant'] != '+')].reset_index(drop=True)
    cols = [f'Intensity {e}' for e in xl + noxl]
    inten = pep[cols].apply(pd.to_numeric).to_numpy(float)
    measured = (inten > 0).any(1)
    pep, inten = pep[measured].reset_index(drop=True), inten[measured]
    log2 = np.where(inten > 0, np.log2(np.where(inten > 0, inten, 1)), np.nan)
    is_xl = np.array([c.startswith('Intensity XL') for c in cols])
    n_valid_xl, n_valid_noxl = (~np.isnan(log2[:, is_xl])).sum(1), (~np.isnan(log2[:, ~is_xl])).sum(1)
    tested = (n_valid_xl >= MIN_VALID) | (n_valid_noxl >= MIN_VALID)
    res = xl_test(log2[tested], is_xl, args.imputations, args.seed)
    log(f'peptides: {n_all:,} in peptides.txt, {len(pep):,} measured after removing reverse hits and contaminants, '
        f'{tested.sum():,} tested; RNA-binding per imputation: min {res["passes_per_imputation"].min()}, '
        f'max {res["passes_per_imputation"].max()}')

    ref = Reference(Path(args.fasta), Path(args.entries), Path(args.history))
    rows, dropped, how_count = [], Counter(), Counter()
    test_index = np.full(len(pep), -1)
    test_index[np.flatnonzero(tested)] = np.arange(tested.sum())
    for i, r in pep.iterrows():
        seq = r['Sequence']
        acc, s, how = ref.place(seq, r['Leading razor protein'], [p for p in r['Proteins'].split(';') if p],
                                [g for g in r['Gene names'].split(';') if g])
        j = test_index[i]
        binding = bool(j >= 0 and res['pass_fraction'][j] >= 0.5)
        if acc is None:
            dropped[(how, binding)] += 1
            continue
        how_count[how] += 1
        ls, le, lmc = lysc_region(ref.seq[acc], s, s + len(seq))
        idt = [r[f'Identification type {e}'] for e in xl + noxl]
        rows.append([
            acc, ref.gene.get(acc, ''), s + 1, s + len(seq), seq,
            '' if ls is None else ls + 1, '' if le is None else le, '' if lmc is None else lmc,
            how, r['Leading razor protein'], r['Proteins'], r['Gene names'], r['Unique (Groups)'],
            r['Missed cleavages'], r['PEP'], r['Score'], r['MS/MS Count'],
            *[int(v) for v in inten[i]],
            sum(t == 'By MS/MS' for t in idt[:4]), sum(t == 'By matching' for t in idt[:4]),
            sum(t == 'By MS/MS' for t in idt[4:]), sum(t == 'By matching' for t in idt[4:]),
            int(n_valid_xl[i]), int(n_valid_noxl[i]), 'true' if j >= 0 else 'false',
            *(['', '', '', ''] if j < 0 else [f'{res["diff"][j]:.4f}', f'{res["p"][j]:.4g}', f'{res["q"][j]:.4g}',
                                              f'{res["pass_fraction"][j]:.2f}']),
            'true' if binding else 'false'])
    rows.sort(key=lambda x: (x[0], x[2], x[3], x[4]))
    columns = ['uniprot_accession', 'gene', 'start', 'end', 'sequence', 'lysc_start', 'lysc_end', 'lysc_missed_cleavages',
               'mapped_by', 'given_leading_razor_protein', 'given_proteins', 'given_gene_names', 'unique_groups',
               'missed_cleavages', 'pep', 'score', 'msms_count',
               *[f'intensity_{e.lower().replace(" ", "_")}' for e in xl + noxl],
               'xl_runs_by_msms', 'xl_runs_by_matching', 'noxl_runs_by_msms', 'noxl_runs_by_matching',
               'xl_valid_values', 'noxl_valid_values', 'tested', 'log2_xl_minus_noxl', 't_test_p', 'permutation_q',
               'pass_fraction', 'rna_binding']
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / TABLE
    with open(path, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(columns)
        w.writerows(rows)
    n_bind = sum(r[-1] == 'true' for r in rows)
    n_bind_all = int((res['pass_fraction'] >= 0.5).sum())
    log(f'{TABLE}: {len(rows):,} peptides on {len({r[0] for r in rows}):,} proteins, {n_bind:,} RNA-binding on '
        f'{len({r[0] for r in rows if r[-1] == "true"}):,} proteins; {path.stat().st_size / 1e6:.2f} MB')
    log(f'  placed by: {dict(how_count)}')
    by_reason = Counter()
    for (h, b), v in dropped.items():
        by_reason[h] += v
    log(f'  dropped: {sum(dropped.values()):,} ({sum(v for (h, b), v in dropped.items() if b):,} of the {n_bind_all:,} '
        f'RNA-binding): {dict(by_reason)}')
    brain_xlsx = fetch(BRAIN_SD1_URL, cache / BRAIN_SD1)
    human_genes = {g.upper() for gg in ref.gene.values() for g in gg.split(';') if g}
    check = check_against_2023(brain_xlsx, rows, human_genes)
    f17, n17 = check['flagged_2017'], check['not_flagged']
    log(f'check against brain-pCLAP Supplementary Data 1: {json.dumps(check)}')

    saved = lambda p: dt.date.fromtimestamp(p.stat().st_mtime).isoformat()  # noqa: E731
    raw = [[TARBALL, tarball.stat().st_size, file_hash(tarball), saved(tarball)]]
    raw += [[f'{TARBALL}:{m}', len(b), md5_bytes(b), saved(tarball)] for m, b in members.items()]
    raw += [[p.name, p.stat().st_size, file_hash(p), saved(p)] for p in (project_json, files_json, brain_xlsx)]
    raw_path = out / RAW_TABLE
    with open(raw_path, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['path', 'bytes', 'md5', 'saved'])
        w.writerows(raw)

    day = saved(tarball)
    search_day = dt.datetime.strptime(par.get('Date of writing', '').split()[0], '%m/%d/%Y').date().isoformat()
    summary = {
        'peptides_in_file': n_all, 'measured': len(pep), 'tested': int(tested.sum()), 'rna_binding_before_mapping': n_bind_all,
        'kept': len(rows), 'kept_rna_binding': n_bind, 'dropped': sum(dropped.values()), 'dropped_by_reason': dict(by_reason),
        'dropped_rna_binding': sum(v for (h, b), v in dropped.items() if b), 'placed_by': dict(how_count),
        'proteins': len({r[0] for r in rows}), 'proteins_rna_binding': len({r[0] for r in rows if r[-1] == 'true'}),
        'variable_modifications': var_mods, 'methylarginine_searched': bool(methyl_r), 'pride_licence': project.get('license'),
        'maxquant_version': par.get('Version'), 'passes_per_imputation': res['passes_per_imputation'].tolist(),
        'check_against_brain_pclap_2023': check}
    manifest_rows = [
        {'file': f'{out.name}/{TABLE}', 'bytes': path.stat().st_size, 'md5': file_hash(path), 'rows': len(rows),
         'redistribution': 'checksum_only',
         'content': f'pCLAP (HEK293; UV crosslinking, Lys-C digest, oligo(dT) capture of poly(A) RNA, RNase and trypsin) human '
                    f'peptides measured in 4 UV-crosslinked (XL) and 4 non-crosslinked (No XL) pull-downs, from the PRIDE MaxQuant '
                    f'search: each peptide on the canonical sequence (exact, unique match) with its Lys-C peptide, the intensity '
                    f'and identification type (MS/MS or matching) in each run, and the XL against No XL test reproduced from the '
                    f'authors\' rule (>= {MIN_VALID} of 4 valid values in a group; Perseus imputation; two-sided t-test with '
                    f's0 {S0}; permutation FDR {FDR:.0%}; XL > No XL; RNA-binding when it passes in at least half of '
                    f'{args.imputations} seeded imputations). {n_bind:,} RNA-binding peptides on '
                    f'{summary["proteins_rna_binding"]:,} proteins; the others are the measured background. The search did not '
                    f'allow methylarginine',
         'source': f'{TARBALL_URL} (MaxQuantOutput/peptides.txt)',
         'version_or_date': f'PRIDE {PRIDE_ACC} (public {project.get("publicationDate")}; MaxQuant {par.get("Version")}, '
                            f'search written {search_day}), downloaded {day}',
         'citation': f'{CITATION}; rule from {RULE_CITATION}',
         'licence': f'{project.get("license")} (PRIDE {PRIDE_ACC}: no licence from the submitters); the article is not open '
                    f'access (ACS; not AuthorChoice, not in PMC). Not redistributed',
         'derived_from': f'{TARBALL}: MaxQuantOutput/peptides.txt, parameters.txt, summary.txt ({RAW_TABLE})',
         'how_made_or_checked': f'scripts/43d_stage_pclap.py (tarball SHA-1 checked against PRIDE; seed {args.seed}; '
                                f'{len(rows):,} peptides kept, {sum(dropped.values()):,} dropped: '
                                + ', '.join(f'{v:,} {k.replace("_", " ")}' for k, v in sorted(by_reason.items()))
                                + f'). Checked against the 2017 HEK293-pCLAP column of {RULE_CITATION.split(",")[0]} 2023 '
                                  f'Supplementary Data 1 (gene level, its 526 brain RBPs): of {f17["same_human_symbol"]} genes marked '
                                  f'found, {f17["rna_binding_here"]} have an RNA-binding peptide here; of {n17["same_human_symbol"]} '
                                  f'not marked, {n17["rna_binding_here"]} do',
         'original_name': 'peptides.txt'},
        {'file': f'{out.name}/{RAW_TABLE}', 'bytes': raw_path.stat().st_size, 'md5': file_hash(raw_path), 'rows': len(raw),
         'redistribution': 'git',
         'content': f'The files behind {TABLE}: the PRIDE {PRIDE_ACC} MaxQuant output archive, the three tables read from '
                    f'it, PRIDE\'s project and file listings (licence; the archive\'s SHA-1) and the brain-pCLAP supplement used '
                    f'as a check, with their size, MD5 and date',
         'source': f'{TARBALL_URL}; {PRIDE_API}; {BRAIN_SD1_URL}', 'version_or_date': day, 'citation': '', 'licence': '',
         'derived_from': '',
         'how_made_or_checked': 'scripts/43d_stage_pclap.py', 'original_name': ''}]
    if args.manifest:
        write_manifest(Path(args.manifest), manifest_rows)
    return summary


def write_manifest(manifest: Path, new: list[dict]):
    from common import update_manifest
    update_manifest(manifest, new)
    log(f'{manifest.name}: {len(new)} rows written')


# ---- preview: methylarginines against the other arginines of the same proteins ---------------------------------------------

def mantel_haenszel(tables: list[tuple[int, int, int, int]]) -> dict:
    """Mantel-Haenszel odds ratio over 2x2 strata (a: methylated & in, b: methylated & out, c: other & in, d: other & out),
    Robins-Breslow-Greenland 95 % CI and the Cochran-Mantel-Haenszel test (no continuity correction)."""
    A = B = C = D = 0
    sr = ss = spr = sps_qr = sqs = 0.0
    sa = se = sv = 0.0
    informative = 0
    for a, b, c, d in tables:
        A, B, C, D = A + a, B + b, C + c, D + d
        n = a + b + c + d
        if n < 2:
            continue
        r, s = a * d / n, b * c / n
        p, q = (a + d) / n, (b + c) / n
        sr, ss = sr + r, ss + s
        spr, sps_qr, sqs = spr + p * r, sps_qr + p * s + q * r, sqs + q * s
        v = (a + b) * (c + d) * (a + c) * (b + d) / (n * n * (n - 1))
        if v > 0:
            informative += 1
        sa, se, sv = sa + a, se + (a + b) * (a + c) / n, sv + v
    out = {'a': A, 'b': B, 'c': C, 'd': D, 'informative_proteins': informative, 'proteins': len(tables)}
    if sr > 0 and ss > 0:
        orr = sr / ss
        var = spr / (2 * sr ** 2) + sps_qr / (2 * sr * ss) + sqs / (2 * ss ** 2)
        out.update(OR=orr, ci_low=math.exp(math.log(orr) - 1.96 * math.sqrt(var)), ci_high=math.exp(math.log(orr) + 1.96 * math.sqrt(var)))
    else:
        out.update(OR=float('nan'), ci_low=float('nan'), ci_high=float('nan'))
    out['chi2'] = (sa - se) ** 2 / sv if sv > 0 else float('nan')
    out['p'] = float(stats.chi2.sf(out['chi2'], 1)) if sv > 0 else float('nan')
    return out


def preview(args) -> list[dict]:
    ref_seq = Reference(Path(args.fasta), Path(args.entries), Path(args.history)).seq
    pep = pd.read_csv(Path(args.out) / TABLE, sep='\t', dtype=str, keep_default_na=False)
    pep[['start', 'end']] = pep[['start', 'end']].astype(int)
    sites = pd.read_csv(args.sites, sep='\t', dtype=str, keep_default_na=False, usecols=['uniprot_accession', 'position', 'residue'])
    methyl_d = defaultdict(set)
    off = 0
    for acc, pos in zip(sites['uniprot_accession'], sites['position'].astype(int)):
        s = ref_seq.get(acc, '')
        if pos <= len(s) and s[pos - 1] == 'R':
            methyl_d[acc].add(pos)
        else:
            off += 1
    methyl = dict(methyl_d)

    def intervals(df, a, b):
        d = defaultdict(list)
        for acc, s, e in zip(df['uniprot_accession'], df[a], df[b]):
            if s != '' and e != '':
                d[acc].append((int(s), int(e)))
        return dict(d)

    bind = pep[pep['rna_binding'] == 'true']
    spans = intervals(bind, 'start', 'end')
    lysc = intervals(bind, 'lysc_start', 'lysc_end')
    covered = intervals(pep, 'start', 'end')

    def dist(p, iv):
        return min((0 if s <= p <= e else (s - p if p < s else p - e) for s, e in iv), default=10 ** 9)

    def no_rg_gr(seq, p):
        return not (seq[p - 2:p - 1] == 'G' or seq[p:p + 1] == 'G')

    sp, cv, lc = (lambda a: spans.get(a, [])), (lambda a: covered.get(a, [])), (lambda a: lysc.get(a, []))
    # each test: which arginines enter (True = all) and which count as "in"
    tests = {
        'a_inside_peptide': lambda acc, p, seq: (True, dist(p, sp(acc)) == 0),
        'b_within_5': lambda acc, p, seq: (True, dist(p, sp(acc)) <= 5),
        'c_within_10': lambda acc, p, seq: (True, dist(p, sp(acc)) <= 10),
        'a_inside_peptide_no_RG_GR': lambda acc, p, seq: (no_rg_gr(seq, p), dist(p, sp(acc)) == 0),
        'a_inside_peptide_covered_by_measured_peptide': lambda acc, p, seq: (dist(p, cv(acc)) == 0, dist(p, sp(acc)) == 0),
        'inside_lysc_region_of_rna_binding_peptide': lambda acc, p, seq: (True, dist(p, lc(acc)) == 0),
    }
    results = []
    for name, f in tests.items():
        tables = []
        for acc in sorted(spans):
            seq, m = ref_seq[acc], methyl.get(acc, set())
            t = [0, 0, 0, 0]
            for i, aa in enumerate(seq):
                if aa != 'R':
                    continue
                keep, inside = f(acc, i + 1, seq)
                if keep:
                    t[(0 if (i + 1) in m else 2) + (0 if inside else 1)] += 1
            tables.append(tuple(t))
        r = mantel_haenszel(tables)
        r['test'] = name
        results.append(r)
    prot_bind = set(spans)
    m_in = [(acc, p) for acc in sorted(prot_bind & set(methyl)) for p in methyl[acc]]
    m_cov = [(a, p) for a in sorted(set(covered) & set(methyl)) for p in methyl[a]]
    cov = {
        'proteins_with_rna_binding_peptide': len(prot_bind),
        'proteins_with_measured_peptide': len(covered),
        'methylated_proteins_in_release': len(methyl),
        'methylarginines_in_release': sum(len(v) for v in methyl.values()),
        'release_sites_not_R_on_reference': off,
        'methylated_proteins_with_rna_binding_peptide': len(prot_bind & set(methyl)),
        'methylarginines_in_proteins_with_rna_binding_peptide': len(m_in),
        'methylarginines_inside_rna_binding_peptide': sum(dist(p, sp(a)) == 0 for a, p in m_in),
        'methylarginines_within_5': sum(dist(p, sp(a)) <= 5 for a, p in m_in),
        'methylarginines_within_10': sum(dist(p, sp(a)) <= 10 for a, p in m_in),
        'methylarginines_inside_lysc_region': sum(dist(p, lc(a)) == 0 for a, p in m_in),
        'methylated_proteins_with_any_measured_peptide': len(set(covered) & set(methyl)),
        'methylarginines_in_proteins_with_any_measured_peptide': len(m_cov),
        'methylarginines_covered_by_measured_peptide': sum(dist(p, cv(a)) == 0 for a, p in m_cov),
    }
    log(json.dumps(cov, indent=1))
    log('test\ta\tb\tc\td\tOR\t95% CI\tp\tinformative proteins / proteins')
    for r in results:
        log(f'{r["test"]}\t{r["a"]}\t{r["b"]}\t{r["c"]}\t{r["d"]}\t{r["OR"]:.3f}\t{r["ci_low"]:.3f}-{r["ci_high"]:.3f}\t'
            f'{r["p"]:.3g}\t{r["informative_proteins"]} / {r["proteins"]}')
    if args.preview_out:
        Path(args.preview_out).write_text(json.dumps({'coverage': cov, 'tests': results, 'sites': str(args.sites)}, indent=1) + '\n')
    return results


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cache', default=str(root / 'PTM_data' / 'external' / 'pclap'), help='download folder (kept between runs)')
    ap.add_argument('--out', default=str(root / 'inputs' / 'interactions'))
    ap.add_argument('--manifest', default=str(root / 'inputs' / 'MANIFEST.tsv'), help="'' to leave the manifest alone")
    ap.add_argument('--fasta', default=str(root / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03.fasta'))
    ap.add_argument('--entries', default=str(root / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03_entries.tsv'))
    ap.add_argument('--history', default=str(root / 'inputs' / 'reference' / 'uniprot_accession_history_2026_03.tsv'))
    ap.add_argument('--imputations', type=int, default=25, help='seeded imputations; RNA-binding = passes in at least half')
    ap.add_argument('--seed', type=int, default=2017)
    ap.add_argument('--no-stage', action='store_true', help='only run the preview on an existing table')
    ap.add_argument('--preview', action='store_true', help='compare methylarginines with the other arginines (prints)')
    ap.add_argument('--sites', default=str(root / 'releases' / 'v1.2.0' / 'human_methylarginine_sites_v1.2.0.tsv'))
    ap.add_argument('--preview-out', default='', help='also write the preview as JSON here (outside inputs/)')
    args = ap.parse_args()
    if not args.no_stage:
        log(json.dumps(stage(args), indent=1))
    if args.preview:
        preview(args)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Stage what is known of the RNA and DNA that methylated proteins bind, beyond structures (inputs/nucleic_acid/).

    python scripts/43_stage_nucleic_acid_binding.py --cache PTM_data/external/nucleic_acid --out inputs/nucleic_acid \
        [--only rbns jaspar eclip htrselex]

Four sources, each mapped to the reviewed human reference (--entries, release 2026_03):

  RNA Bind-n-Seq (ENCODE; Dominguez et al., Mol Cell 70:854, 2018, and ENCODE4): what a purified RBP picks from random
      RNA in vitro. Kept: the 5-mer enrichment R at the protein concentration where the top 5-mer is most enriched,
      and the logos ENCODE gives as position weight matrices.
  JASPAR 2026 CORE vertebrates, non-redundant (CC BY 4.0): DNA-binding profiles, with each matrix's UniProt
      accessions, species, TF class and family from the release's database dump. A matrix of another species is
      given to the human entry of the same gene symbol; a heterodimer's matrix (A::B) to both proteins.
  ENCODE eCLIP (Van Nostrand et al., Nature 583:711, 2020; ENCODE4): reproducible peaks (GRCh38, replicates 1 and 2)
      of RBPs in K562 and HepG2 cells, placed by their middle on GENCODE v50 basic transcripts: per experiment the
      peaks in CDS, 3' UTR, 5' UTR, non-coding exon, intron or outside genes (that order when a peak falls in more
      than one), and the 200 genes with most peaks.
  HTR-SELEX (Jolma et al., Genome Res 30:962, 2020; CC BY-NC 4.0): what a purified RBP picks in vitro from random RNA
      of 40 nucleotides, over several selection cycles, as position frequency matrices of linear motifs (Supplemental
      Table S2) and of motifs within a stem-loop (Table S3), primary and secondary, with the experiment and its
      replicate. Taken from the Europe PMC copy of the supplements; each RBP mapped by its gene symbol (two renamed since).

ENCODE data have no use restrictions (encodeproject.org/help/citing-encode). RNAcompete (CISBP-RNA) states no licence,
so it is not staged; the website links to it. Every download is listed with its MD5 in raw_files.tsv, and each table
gets its row in inputs/MANIFEST.tsv.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import heapq
import json
import re
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

UA = 'PTM_Pipeline-stager/1.0 (https://github.com/Shechterlab/PTM_Pipeline)'
ENCODE = 'https://www.encodeproject.org'
JASPAR_PFM = 'https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt'
JASPAR_SQL = 'https://jaspar.elixir.no/download/database/JASPAR2026.sql.gz'
GENCODE = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.basic.annotation.gtf.gz'
ENCODE_TERMS = 'ENCODE data: no restrictions on use (encodeproject.org/help/citing-encode)'
# Jolma et al. 2020 supplements through Europe PMC (PMC's own links answer scripts with a browser check)
JOLMA_SUPPL = 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7397871/supplementaryFiles'
JOLMA_TABLES = 'supp_gr.258848.119_Supplemental_Tables_S1-S8.xlsx'
# Symbols of the paper that HGNC has renamed since (UniProt primary gene names, release 2026_03)
JOLMA_RENAMED = {'CSDA': 'YBX3', 'HNRPLL': 'HNRNPLL'}
# Where an eCLIP peak's middle falls, best first when it falls in more than one.
REGIONS = ['CDS', 'UTR3', 'UTR5', 'noncoding_exon', 'intron']
TOP_GENES = 200
MANIFEST_COLUMNS = ['file', 'bytes', 'md5', 'rows', 'redistribution', 'content', 'source', 'version_or_date', 'citation', 'licence',
                    'derived_from', 'how_made_or_checked', 'original_name']


def log(*a):
    print(*a, flush=True)


def md5_of(path: Path) -> str:
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


class Downloads:
    """Downloads into the cache, each once, listed for raw_files.tsv."""

    def __init__(self, cache: Path):
        self.cache = cache
        self.seen = {}

    def get(self, url: str, name: str, kept: str) -> Path:
        dest = self.cache / name
        if not (dest.exists() and dest.stat().st_size):
            dest.parent.mkdir(parents=True, exist_ok=True)
            for k in range(5):
                try:
                    req = urllib.request.Request(url, headers={'User-Agent': UA, 'Accept': '*/*'})
                    with urllib.request.urlopen(req, timeout=300) as r, open(dest.with_suffix(dest.suffix + '.part'), 'wb') as fh:
                        while chunk := r.read(1 << 20):
                            fh.write(chunk)
                    dest.with_suffix(dest.suffix + '.part').replace(dest)
                    break
                except Exception as e:  # noqa: BLE001 - retried, then raised
                    if k == 4:
                        raise
                    log(f'  retry {url[:100]}: {e}')
                    time.sleep(2 ** (k + 1))
        self.seen[name] = (url, dest, kept)
        return dest

    def json(self, url: str, name: str) -> dict:
        return json.loads(self.get(url, name, 'cache').read_text())

    def rows(self):
        for name, (url, dest, kept) in sorted(self.seen.items()):
            yield [name, url, dest.stat().st_size, md5_of(dest), dt.date.fromtimestamp(dest.stat().st_mtime).isoformat(), kept]


class Reference:
    def __init__(self, entries: Path):
        self.accs, self.by_gene = set(), defaultdict(set)
        with open(entries, encoding='utf-8') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                self.accs.add(r['Entry'])
                for g in (r.get('Gene Names (primary)') or '').replace(';', ' ').split():
                    self.by_gene[g.upper()].add(r['Entry'])

    def symbol(self, s: str) -> set:
        return set(self.by_gene.get((s or '').upper(), ()))

    def encode_target(self, target: dict) -> tuple[set, str]:
        """An ENCODE target's reviewed accessions: its UniProtKB cross-references in the reference, else its symbols."""
        accs, how = set(), 'uniprot'
        for g in (target or {}).get('genes') or []:
            up = {x.split(':', 1)[1] for x in g.get('dbxrefs') or [] if x.startswith('UniProtKB:')} & self.accs
            if up:
                accs |= up
            else:
                accs |= self.symbol(g.get('symbol', ''))
                how = 'symbol'
        return accs, how


class Out:
    """Tables written to the inputs folder, with their MANIFEST rows."""

    def __init__(self, folder: Path):
        self.folder = folder
        folder.mkdir(parents=True, exist_ok=True)
        self.manifest = []

    def table(self, name: str, columns: list, rows: list, **meta) -> int:
        path = self.folder / name
        with open(path, 'w', encoding='utf-8', newline='') as fh:
            w = csv.writer(fh, delimiter='\t', lineterminator='\n')
            w.writerow(columns)
            for r in rows:
                if len(r) != len(columns):
                    raise ValueError(f'{name}: {len(r)} values for {len(columns)} columns')
                w.writerow(['' if v is None else v for v in r])
        self.manifest.append({'file': f'{self.folder.name}/{name}', 'bytes': path.stat().st_size, 'md5': md5_of(path), 'rows': len(rows),
                              'redistribution': 'git', **meta})
        log(f'{name}: {len(rows):,} rows, {path.stat().st_size / 1e6:.1f} MB')
        return len(rows)

    def write_manifest(self, manifest: Path):
        lines = manifest.read_text(encoding='utf-8').splitlines()
        head = lines[0].split('\t')
        names = {m['file'] for m in self.manifest}
        keep = [ln for ln in lines[1:] if ln.split('\t')[0] not in names]
        keep += ['\t'.join(str(m.get(c, '')) for c in head) for m in self.manifest]
        manifest.write_text('\n'.join([lines[0]] + keep) + '\n', encoding='utf-8')


# ---- RNA Bind-n-Seq ------------------------------------------------------------------------------------------------------

def stage_rbns(dl: Downloads, ref: Reference, out: Out, workers: int):
    q = '&'.join(['type=Experiment', 'assay_title=RNA+Bind-n-Seq', 'status=released', 'limit=all', 'format=json'] +
                 [f'field={f}' for f in ('accession', 'target.label', 'target.genes.dbxrefs', 'target.genes.symbol', 'award.rfa', 'date_released')])
    exps = {e['accession']: e for e in dl.json(f'{ENCODE}/search/?{q}', 'rbns/experiments.json')['@graph'] if e.get('target')}
    q = '&'.join(['type=File', 'assay_title=RNA+Bind-n-Seq', 'status=released', 'output_type=enrichment', 'limit=all', 'format=json'] +
                 [f'field={f}' for f in ('accession', 'dataset', 'submitted_file_name', 'href', 'md5sum', 'file_size')])
    files = dl.json(f'{ENCODE}/search/?{q}', 'rbns/files.json')['@graph']
    kmer_files, pwm_files = {}, []
    for f in files:
        exp = f['dataset'].split('/')[2]
        if exp not in exps:
            continue
        name = f.get('submitted_file_name', '').rsplit('/', 1)[-1]
        if re.search(r'enrichment_R\.5mers', name):
            kmer_files[exp] = f
        m = re.search(r'_(\d+)mer_logo(\d+)\.PWM$', name)
        if m:
            pwm_files.append((exp, f, int(m.group(1)), int(m.group(2))))

    def fetch(f):
        p = dl.get(ENCODE + f['href'], f"rbns/{f['accession']}.tsv", 'cache')
        if md5_of(p) != f['md5sum']:
            raise SystemExit(f"{f['accession']}: MD5 differs from ENCODE's")
        return p

    with ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(fetch, list(kmer_files.values()) + [f for _, f, _, _ in pwm_files]))
    who = {}
    for exp, e in exps.items():
        accs, how = ref.encode_target(e['target'])
        who[exp] = (sorted(accs), e['target'].get('label', ''), how, (e.get('award') or {}).get('rfa', ''))
    rows, chosen = [], []
    for exp, f in sorted(kmer_files.items()):
        accs, label, how, rfa = who[exp]
        lines = (dl.cache / f"rbns/{f['accession']}.tsv").read_text().splitlines()
        concs = lines[0].split('\t')[1:]
        table = [(ln.split('\t')[0], [float(x) for x in ln.split('\t')[1:]]) for ln in lines[1:] if ln.strip()]
        col = max(range(len(concs)), key=lambda c: max(v[c] for _, v in table))
        chosen.append((exp, label, concs[col].strip()))
        for acc in accs:
            rows += [[acc, label, exp, f['accession'], concs[col].strip(), kmer, round(v[col], 4)] for kmer, v in table]
    out.table('rbns_5mer_enrichment.tsv', ['uniprot_accession', 'target', 'encode_experiment', 'encode_file', 'concentration', 'kmer', 'R'], rows,
              content='RNA Bind-n-Seq 5-mer enrichment R (frequency in the RNA bound to the RBP over the input) of every 5-mer, at the protein '
                      'concentration where the top 5-mer is most enriched; ENCODE target mapped to reviewed UniProt by its cross-references',
              source=f'{ENCODE} (assay RNA Bind-n-Seq, output enrichment, *_enrichment_R.5mers)', version=f'ENCODE, read {dt.date.today().isoformat()}',
              citation='Dominguez et al., Mol Cell 70:854 (2018), doi:10.1016/j.molcel.2018.05.001; ENCODE Project Consortium, Nature 583:699 (2020)',
              licence=ENCODE_TERMS, derived_from=f'{len(kmer_files)} ENCODE enrichment files (raw_files.tsv)',
              how_made_or_checked='scripts/43_stage_nucleic_acid_binding.py rbns; each file checked against ENCODE\'s MD5')
    rows = []
    for exp, f, k, logo in sorted(pwm_files, key=lambda x: (x[0], x[2], x[3])):
        accs, label, how, rfa = who[exp]
        lines = [ln.split() for ln in (dl.cache / f"rbns/{f['accession']}.tsv").read_text().splitlines()]
        start = next(i for i, ln in enumerate(lines) if ln and ln[0] == 'PO') + 1
        for acc in accs:
            rows += [[acc, label, exp, f['accession'], k, logo, int(ln[0]) + 1, *[round(float(x), 4) for x in ln[1:5]]] for ln in lines[start:] if len(ln) >= 5]
    out.table('rbns_logos.tsv', ['uniprot_accession', 'target', 'encode_experiment', 'encode_file', 'k', 'logo', 'position', 'A', 'C', 'G', 'U'], rows,
              content='RNA Bind-n-Seq logos as ENCODE gives them (position weight matrices of the enriched k-mers, logo 0 the strongest), '
                      'one row per position',
              source=f'{ENCODE} (assay RNA Bind-n-Seq, output enrichment, *_logo*.PWM)', version=f'ENCODE, read {dt.date.today().isoformat()}',
              citation='Dominguez et al., Mol Cell 70:854 (2018), doi:10.1016/j.molcel.2018.05.001', licence=ENCODE_TERMS,
              derived_from=f'{len(pwm_files)} ENCODE PWM files (raw_files.tsv)', how_made_or_checked='scripts/43_stage_nucleic_acid_binding.py rbns')
    n = len({a for exp in kmer_files for a in who[exp][0]})
    log(f'RBNS: {len(kmer_files)} experiments with 5-mer enrichment ({n} proteins), {len(pwm_files)} logos')


# ---- JASPAR ------------------------------------------------------------------------------------------------------------

def sql_values(text: str):
    """The values of one INSERT statement's VALUES(...): numbers and '...' strings ('' an escaped quote)."""
    out, i = [], 0
    while i < len(text):
        c = text[i]
        if c == "'":
            j, buf = i + 1, []
            while True:
                if text[j] == "'" and j + 1 < len(text) and text[j + 1] == "'":
                    buf.append("'")
                    j += 2
                elif text[j] == "'":
                    break
                else:
                    buf.append(text[j])
                    j += 1
            out.append(''.join(buf))
            i = j + 1
        elif c in ', ':
            i += 1
        else:
            j = i
            while j < len(text) and text[j] != ',':
                j += 1
            v = text[i:j].strip()
            out.append(None if v == 'NULL' else v)
            i = j
    return out


def stage_jaspar(dl: Downloads, ref: Reference, out: Out):
    pfm = dl.get(JASPAR_PFM, 'jaspar/' + Path(JASPAR_PFM).name, 'cache').read_text()
    sql = dl.get(JASPAR_SQL, 'jaspar/' + Path(JASPAR_SQL).name, 'cache')
    tables = defaultdict(list)
    with gzip.open(sql, 'rt', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            m = re.match(r"INSERT INTO [`\"]?(MATRIX|MATRIX_PROTEIN|MATRIX_SPECIES|MATRIX_ANNOTATION)[`\"]? VALUES\((.*)\);\s*$", line)
            if m:
                tables[m.group(1)].append(sql_values(m.group(2)))
    internal = {(r[2], str(r[3])): int(r[0]) for r in tables['MATRIX']}
    prot, species, ann = defaultdict(list), defaultdict(list), defaultdict(dict)
    for mid, acc in tables['MATRIX_PROTEIN']:
        prot[int(mid)].append(acc)
    for mid, tax in tables['MATRIX_SPECIES']:
        species[int(mid)].append(str(tax))
    for _, mid, tag, val in tables['MATRIX_ANNOTATION']:
        ann[int(mid)][tag] = val
    mrows, prows, unmapped = [], [], []
    for block in re.split(r'^>', pfm, flags=re.M)[1:]:
        lines = block.strip().splitlines()
        mid_s, name = lines[0].split('\t')[0].strip(), lines[0].split('\t')[1].strip()
        base, version = mid_s.split('.')
        mid = internal.get((base, version))
        if mid is None:
            raise SystemExit(f'{mid_s}: not in the JASPAR 2026 database dump')
        counts = {ln[0]: [float(x) for x in re.findall(r'[\d.]+', ln[1:])] for ln in lines[1:]}
        accs, how = set(), 'accession'
        if '9606' in species[mid]:
            accs = set(prot[mid]) & ref.accs
        if not accs:
            how = 'symbol'
            for part in name.split('::'):
                accs |= ref.symbol(part)
        if not accs:
            unmapped.append(mid_s)
        a = ann[mid]
        for acc in sorted(accs):
            mrows.append([mid_s, name, acc, how, ','.join(species[mid]), a.get('class', '').replace(';', ','), a.get('family', '').replace(';', ','), len(counts['A']),
                          int(round(sum(c[0] for c in counts.values())))])
        prows += [[mid_s, k + 1, *[int(round(counts[b][k])) for b in 'ACGT']] for k in range(len(counts['A']))]
    meta = dict(source=JASPAR_PFM + '; ' + JASPAR_SQL, version='JASPAR 2026 (release 11)',
                citation='Rauluseviciute et al., Nucleic Acids Res 52:D174 (2024) (JASPAR 2024); JASPAR 2026 release',
                licence='CC BY 4.0 (jaspar.elixir.no/downloads)')
    out.table('jaspar2026_matrices.tsv', ['matrix_id', 'name', 'uniprot_accession', 'mapped_by', 'species_taxids', 'tf_class', 'tf_family', 'length', 'sites'], mrows,
              content='JASPAR 2026 CORE vertebrates non-redundant DNA-binding profiles with the human reviewed entry of each: the matrix\'s own '
                      'UniProt accession for a human matrix, else the entry of the same gene symbol (each partner of a heterodimer A::B)',
              derived_from='JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt; JASPAR2026.sql.gz (MATRIX, MATRIX_PROTEIN, MATRIX_SPECIES, MATRIX_ANNOTATION)',
              how_made_or_checked=f'scripts/43_stage_nucleic_acid_binding.py jaspar; {len(unmapped)} matrices with no human entry: {", ".join(unmapped[:12])}', **meta)
    out.table('jaspar2026_pfms.tsv', ['matrix_id', 'position', 'A', 'C', 'G', 'T'], prows,
              content='The position frequency matrices (site counts per position) of the same profiles',
              derived_from='JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt', how_made_or_checked='scripts/43_stage_nucleic_acid_binding.py jaspar', **meta)
    log(f'JASPAR: {len({r[0] for r in mrows})} matrices on {len({r[2] for r in mrows})} proteins; {len(unmapped)} with no human entry')


# ---- HTR-SELEX ---------------------------------------------------------------------------------------------------------

def cell_number(v):
    """A spreadsheet cell as a number (the sheets are read as text), None when empty or not a number."""
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if x == x else None


def htrselex_blocks(sheet):
    """The motifs of a supplemental table of Jolma et al. 2020: (field: value, [[A, C, G, U] per position]) for each block
    under the column header row: a row naming the motif, then its A, C, G and U rows."""
    rows = sheet.values.tolist()
    head = next(i for i, r in enumerate(rows) if r[0] == 'Gene Symbol' and r[1] == 'Barcode')
    fields = [str(v).strip() for v in rows[head] if isinstance(v, str)]
    i = head + 1
    while i < len(rows):
        if [str(rows[i + k][0]).strip() if i + k < len(rows) else '' for k in range(1, 5)] == ['A', 'C', 'G', 'U']:
            meta = {f: rows[i][k] for k, f in enumerate(fields)}
            cols = [k for k in range(1, len(rows[i + 1])) if cell_number(rows[i + 1][k]) is not None]
            yield meta, [[cell_number(rows[i + b][k]) or 0.0 for b in range(1, 5)] for k in cols]
            i += 5
        else:
            i += 1


def stage_htrselex(dl: Downloads, ref: Reference, out: Out):
    import io
    import zipfile

    import pandas as pd
    z = dl.get(JOLMA_SUPPL, 'htrselex/jolma2020_PMC7397871_supplementary_files.zip', 'cache')
    with zipfile.ZipFile(z) as zf:
        book = pd.ExcelFile(io.BytesIO(zf.read(JOLMA_TABLES)))
    mrows, prows, unmapped = [], [], set()
    text = lambda v: '' if v is None or (isinstance(v, float) and v != v) else str(v).strip()  # noqa: E731
    for sheet, kind in (('S2-PWMs-linear', 'linear'), ('S3- PWMs-structured', 'structured')):
        for k, (meta, cols) in enumerate(htrselex_blocks(book.parse(sheet, header=None, dtype=object)), 1):
            gene = text(meta['Gene Symbol'])
            symbol = JOLMA_RENAMED.get(gene, gene)
            accs = ref.symbol(symbol)
            if not accs:
                unmapped.add(gene)
            model = f'jolma2020_{kind}_{k:03d}'
            for acc in sorted(accs):
                mrows.append([model, gene, acc, 'renamed symbol' if symbol != gene else 'symbol', kind, text(meta.get('Motif type')),
                              text(meta.get('seed')), text(meta.get('Mono or multimer')), text(meta.get('Experiment')), text(meta.get('cycle')),
                              text(meta.get('Barcode')), text(meta.get('Batch')), text(next((v for f, v in meta.items() if f.startswith('Number of RBDs')), '')),
                              text(meta.get('Note')), len(cols)])
            # Table S2 gives counts, S3 frequencies: both as frequencies per position
            prows += [[model, j + 1, *[round(v / (sum(c) or 1), 4) for v in c]] for j, c in enumerate(cols)]
    meta = dict(source=JOLMA_SUPPL + ' (' + JOLMA_TABLES + ', sheets S2 and S3)', version='Genome Res 30:962 (2020), Supplemental Tables S2 and S3',
                citation='Jolma et al., Genome Res 30:962 (2020), doi:10.1101/gr.258848.119', licence='CC BY-NC 4.0 (the article and its supplements, PMC7397871)')
    out.table('htrselex_jolma2020_models.tsv', ['model_id', 'gene_symbol', 'uniprot_accession', 'mapped_by', 'motif_kind', 'motif_type', 'seed', 'oligomer',
                                                 'experiment', 'cycle', 'barcode', 'batch', 'rbd_count', 'note', 'length'], mrows,
              content='HTR-SELEX RNA-binding models of human RBPs (in vitro selection from random 40-nucleotide RNA): linear motifs and motifs within a '
                      'stem-loop, primary or secondary, with the experiment (or its replicate), the selection cycle used and the RBP\'s reviewed entry',
              derived_from=JOLMA_TABLES + ' sheets S2-PWMs-linear and S3- PWMs-structured',
              how_made_or_checked=f'scripts/43_stage_nucleic_acid_binding.py htrselex; mapped by gene symbol ({", ".join(f"{a} as {b}" for a, b in JOLMA_RENAMED.items())}); '
                                  f'unmapped: {", ".join(sorted(unmapped)) or "none"}', **meta)
    out.table('htrselex_jolma2020_pfms.tsv', ['model_id', 'position', 'A', 'C', 'G', 'U'], prows,
              content='The position frequency matrices of the same models (Table S2 counts and Table S3 frequencies, both as frequencies per position)',
              derived_from=JOLMA_TABLES, how_made_or_checked='scripts/43_stage_nucleic_acid_binding.py htrselex', **meta)
    log(f'HTR-SELEX: {len({r[0] for r in mrows})} models on {len({r[2] for r in mrows})} proteins; unmapped: {sorted(unmapped)}')


# ---- eCLIP -------------------------------------------------------------------------------------------------------------

def gtf_attrs(text: str) -> dict:
    return {k: v for k, v in re.findall(r'(\w+) "([^"]*)"', text)}


def load_gencode(path: Path):
    """Segments of GENCODE basic transcripts by (chromosome, strand): (start, end, region index, gene) sorted by start,
    and the genes (id, name, type)."""
    tx, exons, cds, genes, gene_ix = {}, defaultdict(list), defaultdict(list), [], {}
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if f[2] not in ('transcript', 'exon', 'CDS'):
                continue
            a = gtf_attrs(f[8])
            t = a['transcript_id']
            if f[2] == 'transcript':
                g = a['gene_id']
                if g not in gene_ix:
                    gene_ix[g] = len(genes)
                    genes.append((g, a.get('gene_name', ''), a.get('gene_type', '')))
                tx[t] = (f[0], f[6], gene_ix[g])
            elif f[2] == 'exon':
                exons[t].append((int(f[3]), int(f[4])))
            else:
                cds[t].append((int(f[3]), int(f[4])))
    segs = defaultdict(list)
    for t, (chrom, strand, g) in tx.items():
        ex = sorted(exons[t])
        cs = sorted(cds[t])
        for (s, e) in ex:
            if not cs:
                segs[(chrom, strand)].append((s, e, 3, g))
                continue
            lo, hi = cs[0][0], cs[-1][1]
            parts = [(s, min(e, lo - 1), 'before'), (max(s, lo), min(e, hi), 'CDS'), (max(s, hi + 1), e, 'after')]
            for a_, b_, where in parts:
                if a_ > b_:
                    continue
                if where == 'CDS':
                    region = 0
                elif (where == 'before') == (strand == '+'):
                    region = 2   # 5' UTR
                else:
                    region = 1   # 3' UTR
                segs[(chrom, strand)].append((a_, b_, region, g))
        for (s1, e1), (s2, e2) in zip(ex, ex[1:]):
            if s2 > e1 + 1:
                segs[(chrom, strand)].append((e1 + 1, s2 - 1, 4, g))
    for rows in segs.values():
        rows.sort()
    return segs, genes


def place_all(segs: dict, genes: list, peaks: list) -> list:
    """For each (chromosome, strand, position): the best region (REGIONS order) covering it and its gene (a protein-coding
    gene first), or None outside genes. One sweep per chromosome and strand over the sorted segments and positions."""
    res = [None] * len(peaks)
    by_key = defaultdict(list)
    for i, (chrom, strand, pos) in enumerate(peaks):
        by_key[(chrom, strand)].append((pos, i))
    for key, plist in by_key.items():
        rows, heap, k = segs.get(key, []), [], 0
        for pos, i in sorted(plist):
            while k < len(rows) and rows[k][0] <= pos:
                s, e, region, g = rows[k]
                heapq.heappush(heap, (e, region, g))
                k += 1
            while heap and heap[0][0] < pos:
                heapq.heappop(heap)
            if heap:
                _, region, g = min(heap, key=lambda x: (x[1], genes[x[2]][2] != 'protein_coding', x[2]))
                res[i] = (region, g)
    return res


def stage_eclip(dl: Downloads, ref: Reference, out: Out, workers: int):
    q = '&'.join(['type=Experiment', 'assay_title=eCLIP', 'status=released', 'limit=all', 'format=json'] +
                 [f'field={f}' for f in ('accession', 'target.label', 'target.genes.dbxrefs', 'target.genes.symbol', 'biosample_ontology.term_name',
                                         'award.rfa', 'date_released')])
    exps = {e['accession']: e for e in dl.json(f'{ENCODE}/search/?{q}', 'eclip/experiments.json')['@graph'] if e.get('target')}
    q = '&'.join(['type=File', 'assay_title=eCLIP', 'status=released', 'output_type=peaks', 'file_type=bed+narrowPeak', 'assembly=GRCh38', 'limit=all',
                  'format=json'] + [f'field={f}' for f in ('accession', 'dataset', 'biological_replicates', 'href', 'md5sum', 'file_size', 'genome_annotation')])
    files = {}
    for f in dl.json(f'{ENCODE}/search/?{q}', 'eclip/files.json')['@graph']:
        exp = f['dataset'].split('/')[2]
        if exp in exps and sorted(f.get('biological_replicates') or []) == [1, 2]:
            if exp in files:
                raise SystemExit(f'{exp}: more than one reproducible peak file')
            files[exp] = f

    def fetch(f):
        p = dl.get(ENCODE + f['href'], f"eclip/{f['accession']}.bed.gz", 'cache')
        if md5_of(p) != f['md5sum']:
            raise SystemExit(f"{f['accession']}: MD5 differs from ENCODE's")
        return p

    with ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(fetch, files.values()))
    gtf = dl.get(GENCODE, 'gencode/' + Path(GENCODE).name, 'cache')
    segs, genes = load_gencode(gtf)
    peaks, owner = [], []   # every peak's middle (1-based), and its experiment and signalValue
    for exp, f in sorted(files.items()):
        with gzip.open(dl.cache / f"eclip/{f['accession']}.bed.gz", 'rt') as fh:
            for line in fh:
                c = line.rstrip('\n').split('\t')
                if len(c) >= 7:
                    start, end = int(c[1]), int(c[2])
                    peaks.append((c[0], c[5], start + (end - start) // 2 + 1))
                    owner.append((exp, float(c[6])))
    placed = place_all(segs, genes, peaks)
    by_exp = defaultdict(list)
    for (exp, signal), hit in zip(owner, placed):
        by_exp[exp].append((signal, hit))
    erows, rrows, grows = [], [], []
    for exp, f in sorted(files.items()):
        e = exps[exp]
        accs, how = ref.encode_target(e['target'])
        n, regions, per_gene = 0, Counter(), defaultdict(lambda: [0, 0.0, set()])
        for signal, hit in by_exp[exp]:
            n += 1
            if hit is None:
                regions['intergenic'] += 1
                continue
            region, g = hit
            regions[REGIONS[region]] += 1
            pg = per_gene[g]
            pg[0] += 1
            pg[1] = max(pg[1], signal)
            pg[2].add(REGIONS[region])
        cell = (e.get('biosample_ontology') or {}).get('term_name', '')
        for acc in sorted(accs):
            erows.append([exp, f['accession'], acc, e['target'].get('label', ''), how, cell, (e.get('award') or {}).get('rfa', ''), e.get('date_released', ''),
                          n, f.get('genome_annotation', '')])
        rrows += [[exp, r, regions.get(r, 0)] for r in REGIONS + ['intergenic']]
        top = sorted(per_gene.items(), key=lambda kv: (-kv[1][0], -kv[1][1], genes[kv[0]][1]))[:TOP_GENES]
        grows += [[exp, genes[g][0], genes[g][1], genes[g][2], v[0], round(v[1], 3), ';'.join(r for r in REGIONS if r in v[2])] for g, v in top]
    meta = dict(source=f'{ENCODE} (assay eCLIP, output peaks, bed narrowPeak, GRCh38, biological replicates 1 and 2); {GENCODE}',
                version=f'ENCODE, read {dt.date.today().isoformat()}; GENCODE release 50 basic',
                citation='Van Nostrand et al., Nature 583:711 (2020), doi:10.1038/s41586-020-2077-3; Frankish et al., Nucleic Acids Res (GENCODE)',
                licence=ENCODE_TERMS + '; GENCODE open access')
    out.table('encode_eclip_experiments.tsv', ['encode_experiment', 'encode_file', 'uniprot_accession', 'target', 'mapped_by', 'cell_line', 'award', 'date_released',
                                               'n_peaks', 'genome_annotation'], erows,
              content='ENCODE eCLIP experiments and the reproducible peak file of each, with the RBP\'s reviewed UniProt entry and its number of peaks',
              derived_from=f'{len(files)} ENCODE peak files (raw_files.tsv)', how_made_or_checked='scripts/43_stage_nucleic_acid_binding.py eclip; each file checked against ENCODE\'s MD5', **meta)
    out.table('encode_eclip_peak_regions.tsv', ['encode_experiment', 'region', 'n_peaks'], rrows,
              content='Where the middle of each reproducible eCLIP peak falls on GENCODE v50 basic transcripts of its strand: CDS, UTR3, UTR5, '
                      'noncoding_exon, intron (that order when it falls in more than one), or intergenic',
              derived_from='the peak files; gencode.v50.basic.annotation.gtf.gz', how_made_or_checked='scripts/43_stage_nucleic_acid_binding.py eclip (load_gencode, place)', **meta)
    out.table('encode_eclip_target_genes.tsv', ['encode_experiment', 'gene_id', 'gene_name', 'gene_type', 'n_peaks', 'max_signal', 'regions'], grows,
              content=f'The {TOP_GENES} genes with most reproducible eCLIP peaks per experiment (ties by the highest peak signalValue, column 7 of '
                      'narrowPeak), with the regions the peaks fall in',
              derived_from='the peak files; gencode.v50.basic.annotation.gtf.gz', how_made_or_checked='scripts/43_stage_nucleic_acid_binding.py eclip', **meta)
    log(f'eCLIP: {len(files)} experiments of {len({r[2] for r in erows})} proteins')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    root = Path(__file__).resolve().parents[1]
    ap.add_argument('--cache', default='PTM_data/external/nucleic_acid')
    ap.add_argument('--out', default=str(root / 'inputs' / 'nucleic_acid'))
    ap.add_argument('--entries', default=str(root / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03_entries.tsv'))
    ap.add_argument('--manifest', default=str(root / 'inputs' / 'MANIFEST.tsv'))
    ap.add_argument('--only', nargs='*', default=['rbns', 'jaspar', 'eclip', 'htrselex'])
    ap.add_argument('--workers', type=int, default=6)
    args = ap.parse_args()
    dl, ref, out = Downloads(Path(args.cache)), Reference(Path(args.entries)), Out(Path(args.out))
    if 'rbns' in args.only:
        stage_rbns(dl, ref, out, args.workers)
    if 'jaspar' in args.only:
        stage_jaspar(dl, ref, out)
    if 'eclip' in args.only:
        stage_eclip(dl, ref, out, args.workers)
    if 'htrselex' in args.only:
        stage_htrselex(dl, ref, out)
    raw = Path(args.out) / 'raw_files.tsv'
    old = {}
    if raw.exists():
        with open(raw, encoding='utf-8') as fh:
            old = {r[0]: r for r in csv.reader(fh, delimiter='\t')}
    rows = {r[0]: r for r in dl.rows()}
    keep = {k: v for k, v in old.items() if k != 'file' and k not in rows}
    with open(raw, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['file', 'url', 'bytes', 'md5', 'downloaded', 'kept'])
        for k in sorted({**keep, **rows}):
            w.writerow(({**keep, **rows})[k])
    out.manifest.append({'file': f'{Path(args.out).name}/raw_files.tsv', 'bytes': raw.stat().st_size, 'md5': md5_of(raw),
                         'rows': len(keep) + len(rows), 'redistribution': 'git',
                         'content': 'Every file scripts/43 downloaded (ENCODE files and searches, JASPAR, GENCODE, the Jolma 2020 supplements), with its URL, size, MD5 and date; '
                                    'the files themselves are not kept here', 'source': 'the sources above',
                         'how_made_or_checked': 'scripts/43_stage_nucleic_acid_binding.py'})
    out.write_manifest(Path(args.manifest))


if __name__ == '__main__':
    main()

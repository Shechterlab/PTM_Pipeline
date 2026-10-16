#!/usr/bin/env python3
"""Export interaction context for the Site Explorer to docs/data/methylarg_interactions.json.

Run after 36_export_site_explorer_data.py (it reads docs/data/methylarg_sites.json for the protein and
site lists). Public sources are downloaded once into --cache (default data/interactions_cache, which Git
ignores) and reused; pass --refresh to download them again.

Sources
  STRING v12.5 human physical links (CC BY 4.0) .......... physical partners and scores
  Complex Portal human complexes, curated + predicted (CC0) complex membership
  UniProt reviewed human proteome (CC BY 4.0) .............. gene symbols, STRING ids, sequences, Tudor domains
  QuickGO / GOA (CC BY 4.0) ................................ RNA binding (GO:0003723) and DNA binding (GO:0003677)
  PDBe-KB interface residues (EMBL-EBI) .................... residues at protein, RNA and DNA interfaces in the PDB
  SIFTS observed UniProt segments (EMBL-EBI) ............... structural coverage, for the interface enrichment test
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import io
import json
import math
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

STRING_VERSION = '12.5'
STRING_LINKS = f'https://stringdb-downloads.org/download/protein.physical.links.v{STRING_VERSION}/9606.protein.physical.links.v{STRING_VERSION}.txt.gz'
COMPLEX_URLS = {
    'curated': 'https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606.tsv',
    'predicted': 'https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606_predicted.tsv',
}
UNIPROT_STREAM = 'https://rest.uniprot.org/uniprotkb/stream'
QUICKGO_DOWNLOAD = 'https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch'
PDBE_IFACE = 'https://www.ebi.ac.uk/pdbe/graph-api/uniprot/interface_residues/{acc}'
SIFTS_OBSERVED = 'https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/uniprot_segments_observed.tsv.gz'

# The nine human protein arginine methyltransferases.
PRMT_GENES = ['PRMT1', 'PRMT2', 'PRMT3', 'CARM1', 'PRMT5', 'PRMT6', 'PRMT7', 'PRMT8', 'PRMT9']
TUDOR_INTERPRO = 'IPR002999'
GO_TERMS = {'rna': 'GO:0003723', 'dna': 'GO:0003677'}
# GO evidence codes counted as experimental support.
GO_EXPERIMENTAL = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP', 'HGI', 'HEP'}

# Partner flag bits (shared with docs/assets/explorer.js).
F_PRMT, F_TUDOR, F_RNA, F_DNA, F_RNA_EXP, F_DNA_EXP = 1, 2, 4, 8, 16, 32
# Site interface bits.
I_PROTEIN, I_SELF, I_RNA, I_DNA = 1, 2, 4, 8

ACC_RE = re.compile(r'^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$')
UA = 'PTM_Pipeline-site-explorer/1.0 (https://github.com/Shechterlab/PTM_Pipeline)'


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def fetch(url: str, dest: Path, refresh: bool, params: dict | None = None, headers: dict | None = None, tries: int = 4) -> Path:
    if dest.exists() and dest.stat().st_size > 0 and not refresh:
        return dest
    if params:
        url = url + '?' + urllib.parse.urlencode(params)
    dest.parent.mkdir(parents=True, exist_ok=True)
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': UA, **(headers or {})})
            with urllib.request.urlopen(req, timeout=300) as r, open(dest.with_suffix(dest.suffix + '.part'), 'wb') as fh:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            dest.with_suffix(dest.suffix + '.part').replace(dest)
            log(f'  downloaded {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)')
            return dest
        except Exception as e:  # noqa: BLE001 - retried, then raised
            if k == tries - 1:
                raise
            log(f'  retry {url[:90]}: {e}')
            time.sleep(2 ** (k + 1))
    return dest


def open_text(path: Path):
    return io.TextIOWrapper(gzip.open(path, 'rb'), encoding='utf-8') if path.suffix == '.gz' else open(path, encoding='utf-8')


def base_acc(acc: str) -> str:
    return acc.split('-')[0]


# ----------------------------------------------------------------------------- sources
def load_uniprot(cache: Path, refresh: bool):
    """Reviewed human proteome: accession -> (gene, STRING id, sequence); plus Tudor-domain accessions."""
    path = fetch(UNIPROT_STREAM, cache / 'uniprot_human_reviewed.tsv.gz', refresh,
                 {'query': '(organism_id:9606) AND (reviewed:true)', 'fields': 'accession,gene_primary,xref_string,sequence',
                  'format': 'tsv', 'compressed': 'true'})
    proteome = {}
    with open_text(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            sid = (row.get('STRING') or '').split(';')[0].strip()
            proteome[row['Entry']] = (row.get('Gene Names (primary)') or '', sid, row.get('Sequence') or '')
    tud = fetch(UNIPROT_STREAM, cache / 'uniprot_human_tudor.tsv', refresh,
                {'query': f'(organism_id:9606) AND (reviewed:true) AND (xref:interpro-{TUDOR_INTERPRO})', 'fields': 'accession', 'format': 'tsv'})
    tudor = {line.split('\t')[0].strip() for line in open(tud, encoding='utf-8').read().splitlines()[1:] if line.strip()}
    log(f'UniProt: {len(proteome):,} reviewed human entries, {len(tudor)} with a Tudor domain')
    return proteome, tudor


def load_go(cache: Path, refresh: bool):
    """QuickGO annotations to GO term or descendants: accession -> {'rna': set(codes), 'dna': set(codes)}."""
    out = defaultdict(lambda: {'rna': set(), 'dna': set()})
    for key, term in GO_TERMS.items():
        path = fetch(QUICKGO_DOWNLOAD, cache / f'quickgo_{key}.tsv', refresh, {
            'goId': term, 'goUsage': 'descendants', 'goUsageRelationships': 'is_a,part_of,occurs_in',
            'taxonId': '9606', 'taxonUsage': 'exact', 'geneProductType': 'protein',
            'selectedFields': 'geneProductId,goId,goEvidence', 'downloadLimit': '200000'}, headers={'Accept': 'text/tsv'})
        n = 0
        with open(path, encoding='utf-8') as fh:
            for row in csv.reader(fh, delimiter='\t'):
                if not row or row[0].startswith('GENE PRODUCT'):
                    continue
                acc = base_acc(row[1])
                code = row[-1].strip()
                out[acc][key].add(code)
                n += 1
        log(f'QuickGO {term}: {n:,} annotations')
    return out


def load_string(cache: Path, refresh: bool, sid_to_acc: dict, keep: set, min_score: int):
    """Physical links with a methylated protein on at least one side: {(a, b): score} with a < b."""
    path = fetch(STRING_LINKS, cache / Path(STRING_LINKS).name, refresh)
    edges = {}
    with open_text(path) as fh:
        next(fh)
        for line in fh:
            s1, s2, sc = line.split()
            sc = int(sc)
            if sc < min_score:
                continue
            a, b = sid_to_acc.get(s1), sid_to_acc.get(s2)
            if not a or not b or a == b or (a not in keep and b not in keep):
                continue
            key = (a, b) if a < b else (b, a)
            if sc > edges.get(key, 0):
                edges[key] = sc
    log(f'STRING v{STRING_VERSION}: {len(edges):,} physical links (score >= {min_score}) touching a methylated protein')
    return edges


def parse_members(field: str):
    """Complex Portal participant list -> list of member slots, each a list of UniProt accessions."""
    slots = []
    for part in field.split('|'):
        part = re.sub(r'\(\d+\)$', '', part.strip())
        if not part:
            continue
        ids = part.strip('[]').split(',') if part.startswith('[') else [part]
        accs = []
        for i in ids:
            a = base_acc(i.strip())
            a = a.split('_')[0]
            if ACC_RE.match(a) and a not in accs:
                accs.append(a)
        if accs:
            slots.append(accs)
    return slots


def load_complexes(cache: Path, refresh: bool):
    out = []
    for kind, url in COMPLEX_URLS.items():
        path = fetch(url, cache / Path(url).name, refresh)
        with open(path, encoding='utf-8') as fh:
            rows = csv.reader(fh, delimiter='\t')
            header = next(rows)
            col = header.index('Expanded participant list') if 'Expanded participant list' in header else 4
            for row in rows:
                if len(row) <= col:
                    continue
                slots = parse_members(row[col] or row[4])
                if slots:
                    out.append({'id': row[0], 'name': row[1], 'predicted': kind == 'predicted', 'slots': slots})
    log(f'Complex Portal: {sum(not c["predicted"] for c in out):,} curated and {sum(c["predicted"] for c in out):,} predicted human complexes')
    return out


def fetch_interfaces(cache: Path, accs: list, refresh: bool, workers: int):
    """PDBe-KB interface residues for each accession, cached one JSON per accession ({} when none)."""
    folder = cache / 'pdbe_interfaces'
    folder.mkdir(parents=True, exist_ok=True)

    def one(acc):
        dest = folder / f'{acc}.json'
        if dest.exists() and not refresh:
            return
        url = PDBE_IFACE.format(acc=acc)
        for k in range(5):
            try:
                req = urllib.request.Request(url, headers={'User-Agent': UA})
                with urllib.request.urlopen(req, timeout=60) as r:
                    dest.write_bytes(r.read())
                return
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    dest.write_text('{}')
                    return
                time.sleep(2 ** k)
            except Exception:  # noqa: BLE001 - network hiccup, retried
                time.sleep(2 ** k)
        log(f'  PDBe interface lookup failed for {acc}')

    todo = [a for a in accs if refresh or not (folder / f'{a}.json').exists()]
    log(f'PDBe interfaces: {len(accs) - len(todo):,} cached, fetching {len(todo):,}')
    done = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for _ in pool.map(one, todo):
            done += 1
            if done % 250 == 0:
                log(f'  {done:,}/{len(todo):,}')
    out = {}
    for acc in accs:
        p = folder / f'{acc}.json'
        if not p.exists():
            continue
        try:
            d = json.loads(p.read_text() or '{}')
        except json.JSONDecodeError:
            continue
        rec = d.get(acc) if isinstance(d, dict) else None
        if rec and rec.get('data'):
            out[acc] = rec['data']
    log(f'PDBe interfaces: {len(out):,} proteins have interface residues')
    return out


def load_observed(cache: Path, refresh: bool, keep: set):
    """UniProt positions observed in any PDB chain (SIFTS), for the proteins in keep."""
    path = fetch(SIFTS_OBSERVED, cache / Path(SIFTS_OBSERVED).name, refresh)
    obs = defaultdict(list)
    with open_text(path) as fh:
        rows = csv.reader(fh, delimiter='\t')
        header = None
        for row in rows:
            if not row or row[0].startswith('#'):
                continue
            if header is None:
                header = [h.strip() for h in row]
                ix = {h: i for i, h in enumerate(header)}
                continue
            acc = row[ix['SP_PRIMARY']]
            if acc in keep:
                try:
                    obs[acc].append((int(row[ix['SP_BEG']]), int(row[ix['SP_END']])))
                except ValueError:
                    pass
    return obs


# ----------------------------------------------------------------------------- statistics
def log_choose(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(k, N, K, n):
    """P(X >= k) for X ~ Hypergeometric(N population, K successes, n draws)."""
    hi = min(K, n)
    if k > hi:
        return 0.0
    lo = max(0, n - (N - K))
    k = max(k, lo)
    denom = log_choose(N, n)
    terms = [log_choose(K, i) + log_choose(N - K, n - i) - denom for i in range(k, hi + 1)]
    m = max(terms)
    return min(1.0, math.exp(m) * sum(math.exp(t - m) for t in terms))


def fisher_2x2(a, b, c, d):
    """Two-sided Fisher exact p for [[a, b], [c, d]] and the sample odds ratio (Haldane-corrected)."""
    n1, n2, m1 = a + b, c + d, a + c
    N = n1 + n2
    lo, hi = max(0, m1 - n2), min(m1, n1)
    denom = log_choose(N, m1)
    probs = {x: math.exp(log_choose(n1, x) + log_choose(n2, m1 - x) - denom) for x in range(lo, hi + 1)}
    p_obs = probs[a]
    p = min(1.0, sum(v for v in probs.values() if v <= p_obs * (1 + 1e-7)))
    odds = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    return p, odds


def bh(pvals):
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [1.0] * n
    prev = 1.0
    for rank, i in reversed(list(enumerate(order, start=1))):
        prev = min(prev, pvals[i] * n / rank)
        q[i] = prev
    return q


def sig(x, d=3):
    return None if x is None or (isinstance(x, float) and math.isnan(x)) else float(f'{x:.{d}g}')


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--sites-json', default='docs/data/methylarg_sites.json')
    ap.add_argument('--out', default='docs/data/methylarg_interactions.json')
    ap.add_argument('--cache', default='data/interactions_cache')
    ap.add_argument('--string-min-score', type=int, default=700, help='STRING physical score floor (0-1000); 700 = high confidence. Protein pages fetch lower-confidence links live.')
    ap.add_argument('--workers', type=int, default=6, help='parallel PDBe requests')
    ap.add_argument('--refresh', action='store_true', help='download every source again')
    ap.add_argument('--fetch-interfaces-only', action='store_true', help='only fill the PDBe interface cache')
    args = ap.parse_args()
    cache = Path(args.cache)

    explorer = json.loads(Path(args.sites_json).read_text())
    prots, sites = explorer['proteins'], explorer['sites']
    nP = len(prots)
    meth_accs = [p[0] for p in prots]
    canonical = [p[0] for p in prots if not p[10] and '-' not in p[0]]

    if args.fetch_interfaces_only:
        fetch_interfaces(cache, canonical, args.refresh, args.workers)
        return

    proteome, tudor = load_uniprot(cache, args.refresh)
    go = load_go(cache, args.refresh)
    sid_to_acc = {}
    for acc, (_, sid, _) in proteome.items():
        if sid and sid not in sid_to_acc:
            sid_to_acc[sid] = acc
    meth_base = {}
    for i, acc in enumerate(meth_accs):
        meth_base.setdefault(base_acc(acc), i)
    edges = load_string(cache, args.refresh, sid_to_acc, set(meth_base), args.string_min_score)
    complexes = load_complexes(cache, args.refresh)
    ifaces = fetch_interfaces(cache, canonical, args.refresh, args.workers)

    # Protein vocabulary: the explorer's proteins first (same order), then every other protein referenced.
    vocab = list(meth_accs)
    index = {base_acc(a): i for i, a in reversed(list(enumerate(vocab)))}

    def vid(acc):
        acc = base_acc(acc)
        if acc not in index:
            index[acc] = len(vocab)
            vocab.append(acc)
        return index[acc]

    gene_of = {acc: g for acc, (g, _, _) in proteome.items()}
    prmt_accs = {acc for acc, g in gene_of.items() if g in PRMT_GENES}

    # STRING partners per explorer protein.
    partners = defaultdict(list)
    for (a, b), sc in edges.items():
        for x, y in ((a, b), (b, a)):
            if x in meth_base:
                partners[meth_base[x]].append((vid(y), sc))
    # Writers and readers are always in the vocabulary, even without links.
    for acc, g in gene_of.items():
        if g in PRMT_GENES or acc in tudor:
            vid(acc)
    # Partner indices are delta-encoded in ascending order: [d1, score1, d2, score2, ...].
    string_rows = []
    for i in range(nP):
        lst = sorted(partners.get(i, []))
        row, last = [], 0
        for k, sc in lst:
            row += [k - last, sc]
            last = k
        string_rows.append(row)

    # Complexes: membership, methylated share and enrichment against the reviewed proteome.
    N = len(proteome)
    K = sum(1 for a in proteome if a in meth_base)
    cx_rows, cx_p = [], []
    for c in complexes:
        slots = [[a for a in s if a in proteome] or s for s in c['slots']]
        n = len(slots)
        k = sum(1 for s in slots if any(a in meth_base for a in s))
        if k == 0:
            # Tested for the BH correction, but only complexes with a methylated member are exported.
            cx_p.append(hypergeom_sf(0, N, K, n) if n >= 2 else 1.0)
            cx_rows.append(None)
            continue
        members = [vid(s[0]) if len(s) == 1 else [vid(a) for a in s] for s in slots]
        p = hypergeom_sf(k, N, K, n) if n >= 2 else 1.0
        fold = (k / n) / (K / N) if n else 0
        cx_rows.append([c['id'], c['name'], int(c['predicted']), members, n, k, sig(fold), sig(p)])
        cx_p.append(p)
    kinds = [int(c['predicted']) for c in complexes]
    for kind in (0, 1):
        idx = [i for i, kd in enumerate(kinds) if kd == kind]
        for i, q in zip(idx, bh([cx_p[i] for i in idx])):
            if cx_rows[i] is not None:
                cx_rows[i].append(sig(q))
    n_curated, n_predicted = kinds.count(0), kinds.count(1)
    cx_rows = [r for r in cx_rows if r is not None]

    # Interfaces at methylarginine sites, plus all-arginine counts for the enrichment test.
    site_rows = []
    observed = load_observed(cache, args.refresh, set(ifaces))
    enrich = {'rna': [0, 0, 0, 0], 'dna': [0, 0, 0, 0], 'protein': [0, 0, 0, 0], 'any': [0, 0, 0, 0]}
    by_prot_sites = defaultdict(list)
    for j, s in enumerate(sites):
        by_prot_sites[s[0]].append(j)
    for i, acc in enumerate(meth_accs):
        data = ifaces.get(acc)
        if not data:
            continue
        res_partners = defaultdict(dict)   # position -> {partner key: set(pdb ids)}
        for part in data:
            pacc = part.get('accession') or ''
            if pacc in ('RNA', 'DNA'):
                key = pacc
            elif ACC_RE.match(base_acc(pacc)):
                key = base_acc(pacc)
            else:
                continue
            for r in part.get('residues', []):
                if r.get('indexType', 'UNIPROT') != 'UNIPROT':
                    continue
                pdbs = {e.get('pdbId') for e in r.get('interactingPDBEntries') or [] if e.get('pdbId')}
                for pos in range(int(r['startIndex']), int(r['endIndex']) + 1):
                    res_partners[pos].setdefault(key, set()).update(pdbs)

        def mask_of(keys):
            m = 0
            for key in keys:
                m |= I_RNA if key == 'RNA' else I_DNA if key == 'DNA' else I_SELF if key == base_acc(acc) else I_PROTEIN
            return m

        site_pos = {sites[j][1]: j for j in by_prot_sites[i]}
        for pos, j in site_pos.items():
            hit = res_partners.get(pos)
            if not hit:
                continue
            plist, pdbs = [], set()
            for key, ids in hit.items():
                plist.append(-1 if key == 'RNA' else -2 if key == 'DNA' else vid(key))
                pdbs.update(ids)
            site_rows.append([j, mask_of(hit), plist, sorted(pdbs)[:8], len(pdbs)])

        # Enrichment: methylated vs other observed arginines of this protein at each interface type.
        seq = proteome.get(acc, ('', '', ''))[2]
        segs = observed.get(acc, [])
        if not seq or not segs:
            continue
        cov = set()
        for a, b in segs:
            cov.update(range(a, b + 1))
        for pos in cov:
            if pos > len(seq) or seq[pos - 1] != 'R':
                continue
            m = mask_of(res_partners.get(pos, {}))
            meth = pos in site_pos
            for key, bit in (('rna', I_RNA), ('dna', I_DNA), ('protein', I_PROTEIN | I_SELF), ('any', I_RNA | I_DNA | I_PROTEIN | I_SELF)):
                at = bool(m & bit)
                enrich[key][(0 if meth else 2) + (0 if at else 1)] += 1
    enrich_out = {}
    for key, (a, b, c, d) in enrich.items():
        p, odds = fisher_2x2(a, b, c, d)
        enrich_out[key] = [a, b, c, d, sig(math.log2(odds)), sig(p)]

    # Flags for every protein in the vocabulary.
    flags = []
    for acc in vocab:
        b = base_acc(acc)
        f = 0
        if b in prmt_accs:
            f |= F_PRMT
        if b in tudor:
            f |= F_TUDOR
        g = go.get(b)
        if g:
            if g['rna']:
                f |= F_RNA
            if g['rna'] & GO_EXPERIMENTAL:
                f |= F_RNA_EXP
            if g['dna']:
                f |= F_DNA
            if g['dna'] & GO_EXPERIMENTAL:
                f |= F_DNA_EXP
        flags.append(f)

    # GO RNA/DNA binding among methylated vs other reviewed human proteins.
    go_enrich = {}
    for key, bit in (('rna', F_RNA), ('rna_exp', F_RNA_EXP), ('dna', F_DNA), ('dna_exp', F_DNA_EXP)):
        a = b = c = d = 0
        for acc in proteome:
            g = go.get(acc)
            codes = (g['rna'] if key.startswith('rna') else g['dna']) if g else set()
            has = bool(codes & GO_EXPERIMENTAL) if key.endswith('_exp') else bool(codes)
            if acc in meth_base:
                a, b = a + has, b + (not has)
            else:
                c, d = c + has, d + (not has)
        p, odds = fisher_2x2(a, b, c, d)
        go_enrich[key] = [a, a + b, c, c + d, sig(math.log2(odds)), sig(p)]

    # STRING ids (numeric ENSP part) let the browser ask STRING for edges among a protein's partners.
    def ensp(acc):
        sid = proteome.get(base_acc(acc), ('', '', ''))[1]
        m = re.match(r'9606\.ENSP(\d+)$', sid)
        return int(m.group(1)) if m else 0

    out = {
        'meta': {
            'generated': dt.date.today().isoformat(),
            'string_version': STRING_VERSION,
            'string_min_score': args.string_min_score,
            'proteome_size': N,
            'complexes_tested': {'curated': n_curated, 'predicted': n_predicted},
            'methylated_in_proteome': K,
            'prmt_genes': PRMT_GENES,
            'flag_bits': {'prmt': F_PRMT, 'tudor': F_TUDOR, 'rna': F_RNA, 'dna': F_DNA, 'rna_exp': F_RNA_EXP, 'dna_exp': F_DNA_EXP},
            'iface_bits': {'protein': I_PROTEIN, 'self': I_SELF, 'rna': I_RNA, 'dna': I_DNA},
            'go_experimental_codes': sorted(GO_EXPERIMENTAL),
            'sources': {
                'string': f'STRING v{STRING_VERSION} physical links (CC BY 4.0)',
                'complex_portal': 'Complex Portal human curated and predicted complexes (CC0)',
                'uniprot': 'UniProt reviewed human proteome (CC BY 4.0)',
                'go': 'QuickGO/GOA annotations to GO:0003723 and GO:0003677 incl. descendants (CC BY 4.0)',
                'pdbe': 'PDBe-KB interface residues and SIFTS observed segments (EMBL-EBI)',
            },
            'fields': {
                'proteins': ['accession', 'gene', 'flags', 'ensp'],
                'string': 'per explorer protein: [delta partner index, score, ...] with partner indices ascending',
                'complexes': ['id', 'name', 'predicted', 'members (index or list of alternatives)', 'n', 'methylated', 'fold', 'p', 'q'],
                'iface': ['site', 'mask', 'partners (-1 RNA, -2 DNA, else protein index)', 'pdb ids (first 8)', 'pdb count'],
                'enrich': '[methylated at, methylated not at, other R at, other R not at, log2 OR, p] over structurally observed arginines',
                'go_enrich': '[methylated with term, methylated reviewed, other with term, other reviewed, log2 OR, p]',
            },
        },
        'proteins': [[a, gene_of.get(base_acc(a), ''), flags[k], ensp(a)] for k, a in enumerate(vocab)],
        'string': string_rows,
        'complexes': cx_rows,
        'iface': site_rows,
        'enrich': enrich_out,
        'go_enrich': go_enrich,
    }
    # Explorer proteins keep their own gene labels.
    for k, p in enumerate(prots):
        out['proteins'][k][1] = p[1] or out['proteins'][k][1]
    Path(args.out).write_text(json.dumps(out, separators=(',', ':'), ensure_ascii=False))
    log(f'wrote {args.out}: {len(vocab):,} proteins, {sum(len(r) // 2 for r in string_rows):,} partner links, '
        f'{len(cx_rows):,} complexes with methylated members, {len(site_rows):,} sites at PDB interfaces '
        f'({Path(args.out).stat().st_size / 1e6:.2f} MB)')


if __name__ == '__main__':
    main()

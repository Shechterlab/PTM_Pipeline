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
  UniProt disordered regions and domains (CC BY 4.0) ....... strata for the interface test
  CD-CODE condensate membership (cd-code.org API) .......... per-protein confidence score, role and evidence
  RCSB ModelServer + SIFTS chain mapping ................... residues each interface methylarginine touches

Stages that make many requests cache every answer: the PDBe interface lookup (about 5,000 calls) and the
structure contacts (about 5,500 calls). Run them alone first with --fetch-interfaces-only / --fetch-contacts-only.
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
SIFTS_CHAINS = 'https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz'
CDCODE_API = 'https://cd-code.org/api'
MODELSERVER = 'https://models.rcsb.org/v1/{pdb}/residueSurroundings'
RCSB_SEARCH = 'https://search.rcsb.org/rcsbsearch/v2/query'
PDBE_MODRES = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/modified_AA_or_NA/{pdb}'

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
# CD-CODE role and evidence bits (shared with docs/assets/explorer.js).
ROLE_BITS = {'member': 1, 'driver': 2, 'client': 4, 'regulator': 8}
ROLE_OTHER = 16
EVIDENCE_BITS = {'in_vivo': 1, 'in_cellulo': 2, 'in_vitro': 4, 'frap': 8, 'colocalization': 16, 'mass_spectrometry': 32, 'emergent properties': 64}
EVIDENCE_OTHER = 128
# Residue classes for structure contacts.
NUC_RNA = {'A', 'C', 'G', 'U', 'I', 'PSU', '5MC', '7MG', 'OMG', 'OMC', 'OMU', '2MG', 'M2G', 'H2U', '1MA', '5MU', '4SU'}
NUC_DNA = {'DA', 'DC', 'DG', 'DT', 'DI', 'DU', '5CM'}
METHYL_R = {'2MR', 'DA2', 'NMM', 'AGM'}
AA3 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L',
       'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V', 'SEP': 'pS', 'TPO': 'pT',
       'PTR': 'pY', 'MLY': 'Kme2', 'M3L': 'Kme3', 'MLZ': 'Kme1', 'ALY': 'Kac', '2MR': 'Rme2', 'DA2': 'Rme2a', 'NMM': 'Rme1', 'AGM': 'Rme1'}
WATERS = {'HOH', 'WAT', 'DOD'}
SIDE_CHAIN = {'NE', 'CZ', 'NH1', 'NH2', 'CD'}
CONTACT_CUTOFF = 4.0
# Interface geometry: residue environments within ENV_RADIUS, solvent accessibility by Shrake-Rupley.
ENV_RADIUS = 10.0
VDW = {'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8, 'P': 1.8, 'SE': 1.9}
PROBE = 1.4
ARG_MAX_ASA = 274.0        # theoretical maximum for Arg (Tien et al. 2013, PLoS One 8:e80635)
NA_BACKBONE = {'P', 'OP1', 'OP2', 'OP3', "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"}
GUANIDINO_N = {'NE', 'NH1', 'NH2'}
HBOND = 3.5
LEVY = ['none', 'support', 'rim', 'core']

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


def load_uniprot_regions(cache: Path, refresh: bool):
    """Reviewed human proteome: disordered regions (MobiDB-lite) and domain/repeat/zinc-finger intervals."""
    path = fetch(UNIPROT_STREAM, cache / 'uniprot_human_regions.tsv.gz', refresh,
                 {'query': '(organism_id:9606) AND (reviewed:true)', 'fields': 'accession,ft_region,ft_domain,ft_repeat,ft_zn_fing',
                  'format': 'tsv', 'compressed': 'true'})
    feat = re.compile(r'(REGION|DOMAIN|REPEAT|ZN_FING) [<?]?(\d+)\.\.[>?]?(\d+);(?: /note="([^"]*)")?')
    disordered, domains = defaultdict(list), defaultdict(list)
    with open_text(path) as fh:
        for row in csv.reader(fh, delimiter='\t'):
            if not row or row[0] == 'Entry':
                continue
            for m in feat.finditer('\t'.join(row[1:])):
                kind, a, b, note = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4) or ''
                if kind == 'REGION':
                    if note == 'Disordered':
                        disordered[row[0]].append((a, b))
                else:
                    domains[row[0]].append((a, b))
    return disordered, domains


def load_sifts_chains(cache: Path, refresh: bool, keep_pdbs: set):
    """SIFTS chain-level mapping for the PDB entries in keep_pdbs: pdb -> list of (chain, acc, res_beg, sp_beg, sp_end).

    RES_BEG is the SEQRES index (mmCIF label_seq_id), which exists even when the segment's first residue
    is unobserved and has no author number."""
    path = fetch(SIFTS_CHAINS, cache / Path(SIFTS_CHAINS).name, refresh)
    out = defaultdict(list)
    with open_text(path) as fh:
        header = None
        for row in csv.reader(fh, delimiter='\t'):
            if not row or row[0].startswith('#'):
                continue
            if header is None:
                header = {h.strip(): i for i, h in enumerate(row)}
                continue
            pdb = row[header['PDB']].lower()
            if pdb not in keep_pdbs:
                continue
            beg = row[header['RES_BEG']]
            if not re.fullmatch(r'\d+', beg):
                continue
            try:
                out[pdb].append((row[header['CHAIN']], row[header['SP_PRIMARY']], int(beg), int(row[header['SP_BEG']]), int(row[header['SP_END']])))
            except ValueError:
                pass
    return out


def load_cdcode(cache: Path, refresh: bool, uids: list, workers: int):
    """CD-CODE member lists with per-protein confidence, role and evidence, one cached JSON per condensate."""
    folder = cache / 'cdcode'
    folder.mkdir(parents=True, exist_ok=True)
    release = fetch(f'{CDCODE_API}/releases/latest', folder / 'release.json', refresh)

    def one(uid):
        fetch(f'{CDCODE_API}/condensates/{uid}/proteins', folder / f'{uid}.json', refresh, {'per_page': '5000'})

    with ThreadPoolExecutor(max_workers=min(4, workers)) as pool:
        list(pool.map(one, uids))
    out = {}
    for uid in uids:
        p = folder / f'{uid}.json'
        try:
            rows = json.loads(p.read_text()).get('data', [])
        except (OSError, json.JSONDecodeError):
            continue
        members = {}
        for r in rows:
            acc = base_acc((r.get('protein') or {}).get('uniprot_id') or '')
            if not acc:
                continue
            roles = evid = 0
            pmids = []
            for a in (r.get('new_annotations') or []) + (r.get('old_annotations') or []):
                an = a.get('annotations') or {}
                t = an.get('functional_type')
                for x in (t if isinstance(t, list) else [t] if t else []):
                    roles |= ROLE_BITS.get(x, ROLE_OTHER)
                for e in an.get('exp_evidence') or []:
                    evid |= EVIDENCE_BITS.get(e, EVIDENCE_OTHER)
                for pm in ([a.get('publication')] if a.get('publication') else []) + list(a.get('publications') or []):
                    pm = str(pm)
                    if pm.isdigit() and pm not in pmids:
                        pmids.append(pm)
            members[acc] = (r.get('confidence_score') or 0, roles, evid, pmids)
        out[uid] = members
    try:
        version = json.loads(release.read_text()).get('latest_release_version', '')
    except (OSError, json.JSONDecodeError):
        version = ''
    log(f'CD-CODE {version}: member lists for {len(out)} condensates')
    return out, version


def parse_cif_atoms(text: str):
    """Minimal mmCIF atom_site reader: list of dicts with the columns needed for contacts."""
    lines = text.splitlines()
    cols, atoms, in_loop = [], [], False
    for line in lines:
        if line.startswith('_atom_site.'):
            cols.append(line.split('.', 1)[1].strip())
            in_loop = True
            continue
        if in_loop:
            if line.startswith(('ATOM', 'HETATM')):
                parts = line.split()
                if len(parts) >= len(cols):
                    atoms.append(dict(zip(cols, parts)))
            elif line.startswith(('#', 'loop_', '_')) and atoms:
                break
    return atoms


def site_contacts(cache: Path, pdb: str, chain: str, label: int, refresh: bool):
    """Residues of other chains within CONTACT_CUTOFF of an arginine side chain, cached as JSON.

    The residue is addressed by author chain and SEQRES index (label_seq_id). Hits are
    [author chain, label_seq_id, author residue number, residue name, distance]."""
    folder = cache / 'contacts_seqres'
    folder.mkdir(parents=True, exist_ok=True)
    dest = folder / f'{pdb}_{chain}_{label}.json'
    if dest.exists() and not refresh:
        try:
            return json.loads(dest.read_text())
        except json.JSONDecodeError:
            pass
    url = MODELSERVER.format(pdb=pdb) + '?' + urllib.parse.urlencode(
        {'auth_asym_id': chain, 'label_seq_id': label, 'radius': CONTACT_CUTOFF + 1.5, 'model_nums': 1, 'encoding': 'cif'})
    text = None
    for k in range(4):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': UA})
            with urllib.request.urlopen(req, timeout=90) as r:
                text = r.read().decode('utf-8', 'replace')
            break
        except urllib.error.HTTPError as e:
            if e.code in (400, 404):
                text = ''
                break
            time.sleep(2 ** k)
        except Exception:  # noqa: BLE001 - network hiccup, retried
            time.sleep(2 ** k)
    if text is None:
        return None
    atoms = parse_cif_atoms(text)
    me = [a for a in atoms if a.get('auth_asym_id') == chain and a.get('label_seq_id') == str(label)]
    resn = me[0].get('label_comp_id') if me else ''
    side = [(float(a['Cartn_x']), float(a['Cartn_y']), float(a['Cartn_z'])) for a in me if a.get('label_atom_id') in SIDE_CHAIN]
    out = {'resn': resn, 'hits': []}
    if side and resn in ({'ARG'} | METHYL_R):
        seen = {}
        for a in atoms:
            ch, comp = a.get('auth_asym_id'), a.get('label_comp_id')
            if ch == chain or comp in WATERS:
                continue
            x, y, z = float(a['Cartn_x']), float(a['Cartn_y']), float(a['Cartn_z'])
            d = min(math.dist((x, y, z), q) for q in side)
            if d > CONTACT_CUTOFF:
                continue
            key = (ch, a.get('label_seq_id'), a.get('auth_seq_id'), comp)
            if key not in seen or d < seen[key]:
                seen[key] = d
        out['hits'] = [[ch, lab, seq, comp, round(d, 2)] for (ch, lab, seq, comp), d in sorted(seen.items(), key=lambda kv: kv[1])]
    dest.write_text(json.dumps(out))
    return out


def rg_boxes(seq: str):
    """RG/RGG-rich segments: runs of >=3 RG dipeptides with <=12 aa between consecutive RG arginines (as in the explorer)."""
    out, cur = [], None
    for k in range(len(seq) - 1):
        if seq[k] != 'R' or seq[k + 1] != 'G':
            continue
        pos = k + 1
        if cur and pos - cur[1] <= 12:
            cur = [cur[0], pos, cur[2] + 1]
        else:
            if cur and cur[2] >= 3:
                out.append((cur[0], cur[1] + 1))
            cur = [pos, pos, 1]
    if cur and cur[2] >= 3:
        out.append((cur[0], cur[1] + 1))
    return out


def residue_environment(cache: Path, pdb: str, chain: str, label: int, refresh: bool):
    """All atoms of residues within ENV_RADIUS of a residue (RCSB ModelServer), cached as gzipped JSON rows
    [author chain, label_seq_id, author residue number, residue name, atom name, element, x, y, z]."""
    folder = cache / 'environments'
    folder.mkdir(parents=True, exist_ok=True)
    dest = folder / f'{pdb}_{chain}_{label}.json.gz'
    if dest.exists() and not refresh:
        try:
            return json.loads(gzip.decompress(dest.read_bytes()))
        except (OSError, json.JSONDecodeError):
            pass
    url = MODELSERVER.format(pdb=pdb) + '?' + urllib.parse.urlencode(
        {'auth_asym_id': chain, 'label_seq_id': label, 'radius': ENV_RADIUS, 'model_nums': 1, 'encoding': 'cif'})
    text = None
    for k in range(4):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': UA})
            with urllib.request.urlopen(req, timeout=120) as r:
                text = r.read().decode('utf-8', 'replace')
            break
        except urllib.error.HTTPError as e:
            if e.code in (400, 404):
                text = ''
                break
            time.sleep(2 ** k)
        except Exception:  # noqa: BLE001 - network hiccup, retried
            time.sleep(2 ** k)
    if text is None:
        return None
    rows = []
    for a in parse_cif_atoms(text):
        el = (a.get('type_symbol') or '').upper()
        if el in ('H', 'D') or a.get('label_comp_id') in WATERS or a.get('label_alt_id') not in (None, '.', 'A'):
            continue
        try:
            rows.append([a.get('auth_asym_id'), a.get('label_seq_id'), a.get('auth_seq_id'), a.get('label_comp_id'), a.get('label_atom_id'), el,
                         round(float(a['Cartn_x']), 2), round(float(a['Cartn_y']), 2), round(float(a['Cartn_z']), 2)])
        except (KeyError, ValueError):
            continue
    dest.write_bytes(gzip.compress(json.dumps(rows, separators=(',', ':')).encode()))
    return rows


_SPHERE = None


def sphere_points(n=96):
    global _SPHERE
    if _SPHERE is None:
        import numpy as np
        k = np.arange(n) + 0.5
        phi = np.arccos(1 - 2 * k / n)
        theta = np.pi * (1 + 5 ** 0.5) * k
        _SPHERE = np.stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], axis=1)
    return _SPHERE


def residue_sasa(target, occluders):
    """Shrake-Rupley SASA (Å²) of the target atoms given occluding atoms (both lists of (element, xyz))."""
    import numpy as np
    if not target:
        return 0.0
    pts = sphere_points()
    allc = np.array([xyz for _, xyz in occluders]) if occluders else np.zeros((0, 3))
    allr = np.array([VDW.get(el, 1.8) + PROBE for el, _ in occluders]) if occluders else np.zeros(0)
    total = 0.0
    for el, xyz in target:
        r = VDW.get(el, 1.8) + PROBE
        sp = np.asarray(xyz) + pts * r
        if len(allc):
            d0 = np.linalg.norm(allc - np.asarray(xyz), axis=1)
            near = (d0 < r + allr) & (d0 > 1e-3)
            if near.any():
                dd = np.linalg.norm(sp[:, None, :] - allc[near][None, :, :], axis=2)
                free = ~(dd < allr[near][None, :]).any(axis=1)
            else:
                free = np.ones(len(sp), bool)
        else:
            free = np.ones(len(sp), bool)
        total += 4 * math.pi * r * r * free.mean()
    return total


def analyze_environment(rows, chain: str, label: int):
    """Levy interface class, relative accessibility and nucleic-acid / protein hydrogen-bond partners of one arginine."""
    me = [r for r in rows if r[0] == chain and r[1] == str(label)]
    if not me or me[0][3] not in ({'ARG'} | METHYL_R):
        return None
    me_ids = {id(r) for r in me}
    own = [r for r in rows if r[0] == chain and id(r) not in me_ids]
    na = [r for r in rows if r[0] != chain and (r[3] in NUC_RNA or r[3] in NUC_DNA)]
    prot = [r for r in rows if r[0] != chain and r[3] in AA3]
    other = [r for r in rows if r[0] != chain]
    tgt = [(r[5], r[6:9]) for r in me]
    occ = lambda rs: [(r[5], r[6:9]) for r in rs]
    # The residue's own atoms occlude each other (an atom never occludes itself).
    s_mono = float(residue_sasa(tgt, occ(me + own)))
    s_cplx = float(residue_sasa(tgt, occ(me + own + other)))
    s_na = float(residue_sasa(tgt, occ(me + own + na))) if na else s_mono
    s_prot = float(residue_sasa(tgt, occ(me + own + prot))) if prot else s_mono
    r_mono, r_cplx = s_mono / ARG_MAX_ASA, s_cplx / ARG_MAX_ASA
    if s_mono - s_cplx < 1.0:
        levy = 'none'
    elif r_mono < 0.25:
        levy = 'support'
    elif r_cplx < 0.25:
        levy = 'core'
    else:
        levy = 'rim'
    gn = [r for r in me if r[4] in GUANIDINO_N]
    na_type = 'none'
    base_hb, bb_hb, prot_hb = set(), set(), set()
    for r in na + prot:
        if r[5] not in ('N', 'O'):
            continue
        d = min(math.dist(r[6:9], g[6:9]) for g in gn) if gn else 99
        if d > HBOND:
            continue
        if r in prot:
            prot_hb.add((r[0], r[1], r[3]))
        elif r[4] in NA_BACKBONE:
            bb_hb.add((r[0], r[2], r[3]))
        else:
            base_hb.add((r[0], r[2], r[3]))
    if na:
        sc = [r for r in me if r[4] in SIDE_CHAIN]
        close = any(min(math.dist(a[6:9], b[6:9]) for b in sc) <= CONTACT_CUTOFF for a in na) if sc else False
        na_type = 'base' if base_hb else 'backbone' if bb_hb else 'contact' if close else 'none'
    return {
        'resn': me[0][3], 'levy': levy, 'rasa_mono': round(r_mono, 3), 'rasa_cplx': round(r_cplx, 3),
        'd_na': round(s_mono - s_na, 1), 'd_prot': round(s_mono - s_prot, 1), 'na_type': na_type,
        'base_hb': sorted({c + n for _, n, c in [(x[0], x[1], x[2]) for x in base_hb]}),
        'n_prot_hb': len(prot_hb),
    }


def load_modelled_methyl_r(cache: Path, refresh: bool, workers: int):
    """PDB entries that model a methylated arginine (2MR, DA2, NMM, AGM), with chain, SEQRES index and residue name."""
    folder = cache / 'modelled_methyl_r'
    folder.mkdir(parents=True, exist_ok=True)
    q = {'query': {'type': 'terminal', 'service': 'text', 'parameters': {'attribute': 'rcsb_polymer_entity_container_identifiers.chem_comp_monomers', 'operator': 'in', 'value': sorted(METHYL_R)}},
         'return_type': 'entry', 'request_options': {'return_all_hits': True}}
    hits = fetch(RCSB_SEARCH, folder / 'entries.json', refresh, {'json': json.dumps(q)})
    entries = sorted(x['identifier'].lower() for x in json.loads(hits.read_text()).get('result_set', []))

    def one(pdb):
        dest = folder / f'{pdb}.json'
        if dest.exists() and not refresh:
            return
        try:
            fetch(PDBE_MODRES.format(pdb=pdb), dest, refresh)
        except Exception:  # noqa: BLE001 - a missing entry is skipped
            dest.write_text('{}')

    with ThreadPoolExecutor(max_workers=min(4, workers)) as pool:
        list(pool.map(one, entries))
    out = []
    for pdb in entries:
        try:
            rows = json.loads((folder / f'{pdb}.json').read_text()).get(pdb, [])
        except (OSError, json.JSONDecodeError):
            continue
        for r in rows:
            if r.get('chem_comp_id') in METHYL_R and r.get('residue_number') is not None:
                out.append((pdb, r.get('chain_id'), int(r['residue_number']), r['chem_comp_id']))
    log(f'PDB-modelled methylarginines: {len(out)} residues in {len(entries)} entries')
    return out


def cdcode_threshold_enrichment(cd, proteome, meth_base, uids, thresholds=(1, 2, 3, 4)):
    """Methylated share of each condensate counting only members at or above a CD-CODE confidence score."""
    N = len(proteome)
    K = sum(1 for a in proteome if a in meth_base)
    out = {}
    for t in thresholds:
        rows, pv = [], []
        for ci, uid in enumerate(uids):
            members = {a for a, m in cd.get(uid, {}).items() if m[0] >= t and a in proteome}
            n = len(members)
            if n < 5:
                continue
            m = sum(1 for a in members if a in meth_base)
            p, odds = fisher_2x2(m, n - m, K - m, (N - K) - (n - m))
            rows.append([ci, n, m, sig(math.log2(odds)), sig(p)])
            pv.append(p)
        for r, q in zip(rows, bh(pv)):
            r.append(sig(q))
        out[str(t)] = rows
    return out


def mantel_haenszel(tables):
    """Mantel-Haenszel common odds ratio and CMH test (continuity-corrected) over 2x2 strata [[a, b], [c, d]]."""
    num = den = sa = sea = sva = 0.0
    used = 0
    for a, b, c, d in tables:
        n = a + b + c + d
        if n < 2 or (a + b) == 0 or (c + d) == 0 or (a + c) == 0 or (b + d) == 0:
            continue
        used += 1
        num += a * d / n
        den += b * c / n
        sa += a
        sea += (a + b) * (a + c) / n
        sva += (a + b) * (c + d) * (a + c) * (b + d) / (n * n * (n - 1))
    if not used or sva <= 0:
        return None, None, 0
    odds = (num / den) if den > 0 else float('inf')
    stat = (abs(sa - sea) - 0.5) ** 2 / sva
    p = math.erfc(math.sqrt(max(stat, 0) / 2))
    return odds, p, used


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
    ap.add_argument('--fetch-contacts-only', action='store_true', help='only fill the structure-contact cache (needs the interface cache)')
    ap.add_argument('--contact-entries', type=int, default=2, help='PDB entries examined per site and partner type')
    ap.add_argument('--fetch-environments-only', action='store_true', help='only fill the 10 Å residue-environment cache (interface geometry)')
    ap.add_argument('--controls-per-protein', type=int, default=12, help='non-methylated interface arginines sampled per protein for the geometry comparison')
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
    prmt_tudor = {acc for acc, (g, _, _) in proteome.items() if g in PRMT_GENES} | set(tudor)
    edges = load_string(cache, args.refresh, sid_to_acc, set(meth_base) | prmt_tudor, args.string_min_score)
    complexes = load_complexes(cache, args.refresh)
    ifaces = fetch_interfaces(cache, canonical, args.refresh, args.workers)
    disordered, domains = load_uniprot_regions(cache, args.refresh)

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
    observed = load_observed(cache, args.refresh, set(meth_base))
    enrich = {'rna': [0, 0, 0, 0], 'dna': [0, 0, 0, 0], 'protein': [0, 0, 0, 0], 'any': [0, 0, 0, 0]}
    IFACE_TYPES = (('rna', I_RNA), ('dna', I_DNA), ('protein', I_PROTEIN | I_SELF), ('any', I_RNA | I_DNA | I_PROTEIN | I_SELF))
    STRATA = ('all', 'ordered', 'domain')
    strata = {(t, f): defaultdict(lambda: [0, 0, 0, 0]) for t, _ in IFACE_TYPES for f in STRATA}
    contact_jobs = []   # (site j, partner key, pdb, own accession, position)
    env_jobs = []       # (methylated 1/0, protein index, position, site j or -1, pdb, accession)
    RG_CONTEXTS = ('segment', 'dipeptide', 'other')
    rg_strata = {(c, t): defaultdict(lambda: [0, 0, 0, 0]) for c in RG_CONTEXTS for t in ('rna', 'any')}
    rg_cover = {c: [0, 0] for c in RG_CONTEXTS}   # methylated sites in context: [all, structurally observed]
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
                for pdb in sorted(ids, reverse=True)[:args.contact_entries]:
                    contact_jobs.append((j, key, pdb, acc, pos))
            site_rows.append([j, mask_of(hit), plist, sorted(pdbs)[:8], len(pdbs)])
            if pdbs:
                env_jobs.append((1, i, pos, j, max(pdbs), acc))
        # Controls: non-methylated arginines of the same protein at an interface, evenly spaced, capped.
        seq0 = proteome.get(acc, ('', '', ''))[2]
        if any(pos in res_partners for pos in site_pos) and seq0:
            ctrl = sorted(pos for pos in res_partners if pos not in site_pos and pos <= len(seq0) and seq0[pos - 1] == 'R')
            if len(ctrl) > args.controls_per_protein:
                step = len(ctrl) / args.controls_per_protein
                ctrl = [ctrl[int(k * step)] for k in range(args.controls_per_protein)]
            for pos in ctrl:
                ids = set().union(*res_partners[pos].values())
                if ids:
                    env_jobs.append((0, i, pos, -1, max(ids), acc))

        # Enrichment: methylated vs other observed arginines of this protein at each interface type.
        seq = proteome.get(acc, ('', '', ''))[2]
        segs = observed.get(acc, [])
        if not seq or not segs:
            continue
        cov = set()
        for a, b in segs:
            cov.update(range(a, b + 1))
        dis = disordered.get(acc, [])
        dom = domains.get(acc, [])
        boxes = rg_boxes(seq)
        rg_ctx = lambda pos: 'segment' if any(a <= pos <= b for a, b in boxes) else 'dipeptide' if pos < len(seq) and seq[pos] == 'G' else 'other'
        for pos in cov:
            if pos > len(seq) or seq[pos - 1] != 'R':
                continue
            m = mask_of(res_partners.get(pos, {}))
            meth = pos in site_pos
            ordered = not any(a <= pos <= b for a, b in dis)
            in_dom = any(a <= pos <= b for a, b in dom)
            ctx = rg_ctx(pos)
            for key, bit in (('rna', I_RNA), ('any', I_RNA | I_DNA | I_PROTEIN | I_SELF)):
                rg_strata[(ctx, key)][acc][(0 if meth else 2) + (0 if m & bit else 1)] += 1
            for key, bit in IFACE_TYPES:
                cell = (0 if meth else 2) + (0 if m & bit else 1)
                enrich[key][cell] += 1
                strata[(key, 'all')][acc][cell] += 1
                if ordered:
                    strata[(key, 'ordered')][acc][cell] += 1
                if in_dom:
                    strata[(key, 'domain')][acc][cell] += 1
    enrich_out = {}
    for key, (a, b, c, d) in enrich.items():
        p, odds = fisher_2x2(a, b, c, d)
        enrich_out[key] = [a, b, c, d, sig(math.log2(odds)), sig(p)]
    # Stratified versions: pooled within a residue class, and Mantel-Haenszel with each protein as a stratum,
    # so the comparison is always against arginines of the same protein.
    enrich_strata = {}
    for key, _ in IFACE_TYPES:
        enrich_strata[key] = {}
        for f in STRATA:
            tabs = list(strata[(key, f)].values())
            a, b, c, d = (sum(t[k] for t in tabs) for k in range(4))
            p, odds = fisher_2x2(a, b, c, d) if a + b + c + d else (None, None)
            mh_or, mh_p, used = mantel_haenszel(tabs)
            enrich_strata[key][f] = {
                'pooled': [a, b, c, d, sig(math.log2(odds)) if odds else None, sig(p)],
                'mh': [sig(math.log2(mh_or)) if mh_or and math.isfinite(mh_or) else None, sig(mh_p), used],
            }

    # How many methylated arginines in each RG context are structurally observed at all (all canonical proteins).
    for i, acc in enumerate(meth_accs):
        if prots[i][10] or '-' in acc:
            continue
        seq = proteome.get(acc, ('', '', ''))[2]
        if not seq:
            continue
        boxes = rg_boxes(seq)
        segs = observed.get(acc, [])
        for j in by_prot_sites[i]:
            pos = sites[j][1]
            if pos > len(seq) or seq[pos - 1] != 'R':
                continue
            c = 'segment' if any(a <= pos <= b for a, b in boxes) else 'dipeptide' if pos < len(seq) and seq[pos] == 'G' else 'other'
            rg_cover[c][0] += 1
            rg_cover[c][1] += any(a <= pos <= b for a, b in segs)

    # RG/RGG context: methylated vs other arginines of the same protein in each sequence context.
    rg_out = {}
    for c in RG_CONTEXTS:
        rg_out[c] = {'cover': rg_cover[c]}
        for key in ('rna', 'any'):
            tabs = list(rg_strata[(c, key)].values())
            a, b, cc, d = (sum(t[k] for t in tabs) for k in range(4))
            p, odds = fisher_2x2(a, b, cc, d) if a + b + cc + d else (None, None)
            mh_or, mh_p, used = mantel_haenszel(tabs)
            rg_out[c][key] = {'pooled': [a, b, cc, d, sig(math.log2(odds)) if odds else None, sig(p)],
                              'mh': [sig(math.log2(mh_or)) if mh_or and math.isfinite(mh_or) else None, sig(mh_p), used]}

    # Interface geometry from 10 Å residue environments: Levy core/rim/support and nucleic-acid hydrogen bonds,
    # methylated interface arginines vs sampled non-methylated interface arginines of the same proteins.
    geo_sites, geo_summary = [], {}
    if env_jobs:
        echains = load_sifts_chains(cache, args.refresh, {pdb for *_, pdb, _ in env_jobs})
        units = []
        for meth, i, pos, j, pdb, acc in env_jobs:
            cands = sorted({(ch, res_beg + (pos - sp_beg)) for ch, sp, res_beg, sp_beg, sp_end in echains.get(pdb, [])
                            if base_acc(sp) == acc and sp_beg <= pos <= sp_end})
            if cands:
                units.append((meth, i, pos, j, pdb, cands[0][0], cands[0][1]))
        log(f'Interface geometry: {len(units):,} residue environments ({sum(u[0] for u in units):,} methylated)')
        done = 0

        def run_env(u):
            rows = residue_environment(cache, u[4], u[5], u[6], args.refresh)
            return u, (analyze_environment(rows, u[5], u[6]) if rows else None)

        results = []
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            for u, res in pool.map(run_env, units):
                results.append((u, res))
                done += 1
                if done % 500 == 0:
                    log(f'  {done:,}/{len(units):,}')
        if args.fetch_environments_only:
            return
        levy_tabs = defaultdict(lambda: [0, 0, 0, 0])
        base_tabs = defaultdict(lambda: [0, 0, 0, 0])
        levy_counts = {1: {k: 0 for k in LEVY}, 0: {k: 0 for k in LEVY}}
        na_counts = {1: {'base': 0, 'backbone': 0, 'contact': 0}, 0: {'base': 0, 'backbone': 0, 'contact': 0}}
        for (meth, i, pos, j, pdb, ch, lab), res in results:
            if not res:
                continue
            levy_counts[meth][res['levy']] += 1
            if res['levy'] != 'none':
                levy_tabs[i][(0 if meth else 2) + (0 if res['levy'] == 'core' else 1)] += 1
            if res['na_type'] != 'none':
                na_counts[meth][res['na_type']] += 1
                base_tabs[i][(0 if meth else 2) + (0 if res['na_type'] == 'base' else 1)] += 1
            if meth:
                geo_sites.append([j, pdb, ch, LEVY.index(res['levy']), res['rasa_mono'], res['rasa_cplx'], res['na_type'], res['base_hb'], res['d_na'], res['d_prot'], res['n_prot_hb']])

        def summary(tabs):
            tl = list(tabs.values())
            a, b, c, d = (sum(t[k] for t in tl) for k in range(4))
            p, odds = fisher_2x2(a, b, c, d) if a + b + c + d else (None, None)
            mh_or, mh_p, used = mantel_haenszel(tl)
            return {'pooled': [a, b, c, d, sig(math.log2(odds)) if odds else None, sig(p)],
                    'mh': [sig(math.log2(mh_or)) if mh_or and math.isfinite(mh_or) else None, sig(mh_p), used]}

        geo_summary = {'levy_counts': levy_counts, 'na_counts': na_counts, 'core': summary(levy_tabs), 'base': summary(base_tabs),
                       'controls_per_protein': args.controls_per_protein}
        log(f'Interface geometry: {len(geo_sites):,} methylated sites classified')
    elif args.fetch_environments_only:
        return

    # Methylated arginines modelled as such in PDB entries, mapped to UniProt through SIFTS SEQRES ranges.
    modelled = load_modelled_methyl_r(cache, args.refresh, args.workers)
    mchains = load_sifts_chains(cache, args.refresh, {m[0] for m in modelled})
    site_index = {(base_acc(meth_accs[s[0]]), s[1]): j for j, s in enumerate(sites)}
    modelled_sites = defaultdict(set)
    modelled_other = defaultdict(set)
    for pdb, ch, lab, comp in modelled:
        for c2, sp, res_beg, sp_beg, sp_end in mchains.get(pdb, []):
            if c2 != ch:
                continue
            u = sp_beg + (lab - res_beg)
            if not (sp_beg <= u <= sp_end):
                continue
            j = site_index.get((base_acc(sp), u))
            if j is not None:
                modelled_sites[j].add((pdb, comp))
            else:
                modelled_other[base_acc(sp)].add((u, pdb, comp))
    modelled_rows = [[j, sorted({p2 for p2, _ in v}), sorted({c for _, c in v})] for j, v in sorted(modelled_sites.items())]
    modelled_gaps = [[vid(a), sorted({u for u, _, _ in v}), sorted({p2 for _, p2, _ in v})[:6]] for a, v in sorted(modelled_other.items()) if a in proteome]
    log(f'PDB-modelled methylarginines: {len(modelled_rows)} union sites, {sum(len(r[1]) for r in modelled_gaps)} positions not in the union')

    # Residues each interface methylarginine touches, from the structures themselves (RCSB ModelServer).
    contact_rows = []
    if contact_jobs:
        chains = load_sifts_chains(cache, args.refresh, {pdb for _, _, pdb, _, _ in contact_jobs})
        requests_ = defaultdict(list)   # (pdb, chain, author residue) -> [(site, partner key)]
        for j, key, pdb, acc, pos in contact_jobs:
            # The site's residue in at most two copies of the protein in the entry.
            cands = sorted({(ch, res_beg + (pos - sp_beg)) for ch, sp, res_beg, sp_beg, sp_end in chains.get(pdb, [])
                            if base_acc(sp) == acc and sp_beg <= pos <= sp_end})
            for ch, auth in cands[:2]:
                requests_[(pdb, ch, auth)].append((j, key))
        todo = sorted(requests_)
        log(f'Structure contacts: {len(todo):,} residue neighbourhoods to read')
        done = 0
        results = {}

        def run(req):
            return req, site_contacts(cache, req[0], req[1], req[2], args.refresh)

        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            for req, res in pool.map(run, todo):
                results[req] = res
                done += 1
                if done % 500 == 0:
                    log(f'  {done:,}/{len(todo):,}')
        chain_acc = {(pdb, ch): base_acc(sp) for pdb, rows in chains.items() for ch, sp, *_ in rows}
        chain_map = defaultdict(list)
        for pdb, rows in chains.items():
            for ch, sp, beg, sp_beg, sp_end in rows:
                chain_map[(pdb, ch)].append((beg, sp_beg, sp_end))
        by_site = defaultdict(dict)   # j -> {(key, pdb): set(labels)}
        modelled = {}
        for (pdb, ch, auth), owners in requests_.items():
            res = results.get((pdb, ch, auth))
            if not res:
                continue
            for j, _ in owners:
                own = meth_accs[sites[j][0]]
                if res.get('resn') in METHYL_R:
                    modelled[j] = res['resn']
                for hch, hlab, hseq, comp, dist in res['hits']:
                    if comp in NUC_RNA:
                        k, label = -1, comp + hseq
                    elif comp in NUC_DNA:
                        k, label = -2, comp[1:] + hseq if comp.startswith('D') and len(comp) == 2 else comp + hseq
                    elif comp in AA3:
                        pacc = chain_acc.get((pdb, hch))
                        unp = None
                        if re.fullmatch(r'\d+', hlab or ''):
                            for res_beg, sp_beg, sp_end in chain_map.get((pdb, hch), []):
                                u = sp_beg + (int(hlab) - res_beg)
                                if sp_beg <= u <= sp_end:
                                    unp = u
                                    break
                        k = vid(pacc) if pacc else -4
                        label = AA3[comp] + (str(unp) if unp is not None else hch + ':' + hseq)
                    else:
                        k, label = -3, comp
                    by_site[j].setdefault((k, pdb), set()).add(label)
        for j, groups in by_site.items():
            rows = []
            for (k, pdb), labels in sorted(groups.items(), key=lambda kv: (0 if kv[0][0] in (-1, -2) else 1, kv[0])):
                rows.append([k, pdb, sorted(labels, key=lambda x: (re.sub(r'\d', '', x), int(re.sub(r'\D', '', x) or 0)))[:8]])
            contact_rows.append([j, rows, modelled.get(j, '')])
        log(f'Structure contacts: {len(contact_rows):,} sites with residue-level contacts')
    if args.fetch_contacts_only:
        return

    # CD-CODE: this protein's own membership confidence, role and evidence in each condensate.
    cond_meta = explorer['meta'].get('condensates', [])
    cd_uids = [c[1] for c in cond_meta]
    cd, cd_version = load_cdcode(cache, args.refresh, cd_uids, args.workers)
    cond_conf = cdcode_threshold_enrichment(cd, proteome, meth_base, cd_uids)
    cd_rows = []
    for i, acc in enumerate(meth_accs):
        b = base_acc(acc)
        for ci, uid in enumerate(cd_uids):
            m = cd.get(uid, {}).get(b)
            if m:
                score, roles, evid, pmids = m
                cd_rows.append([i, ci, score, roles, evid, len(pmids), pmids])

    # STRING partners of PRMTs and Tudor-domain proteins, for writer-substrate-reader paths.
    wr_lists = defaultdict(set)
    for (a, b), sc in edges.items():
        if a in prmt_tudor:
            wr_lists[a].add(b)
        if b in prmt_tudor:
            wr_lists[b].add(a)
    wr_rows = [[vid(w), sorted(vid(x) for x in xs)] for w, xs in sorted(wr_lists.items())]

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
            'cdcode_version': cd_version,
            'role_bits': dict(ROLE_BITS, other=ROLE_OTHER),
            'evidence_bits': dict(EVIDENCE_BITS, other=EVIDENCE_OTHER),
            'contact_cutoff': CONTACT_CUTOFF,
            'contact_entries': args.contact_entries,
            'go_experimental_codes': sorted(GO_EXPERIMENTAL),
            'sources': {
                'string': f'STRING v{STRING_VERSION} physical links (CC BY 4.0)',
                'complex_portal': 'Complex Portal human curated and predicted complexes (CC0)',
                'uniprot': 'UniProt reviewed human proteome (CC BY 4.0)',
                'go': 'QuickGO/GOA annotations to GO:0003723 and GO:0003677 incl. descendants (CC BY 4.0)',
                'pdbe': 'PDBe-KB interface residues and SIFTS observed segments and chain mapping (EMBL-EBI)',
                'cdcode': f'CD-CODE {cd_version} condensate membership via cd-code.org API',
                'contacts': 'RCSB ModelServer residue surroundings (wwPDB data, CC0)',
                'regions': 'UniProt disordered regions (MobiDB-lite) and domain, repeat and zinc-finger features (CC BY 4.0)',
            },
            'fields': {
                'proteins': ['accession', 'gene', 'flags', 'ensp'],
                'string': 'per explorer protein: [delta partner index, score, ...] with partner indices ascending',
                'complexes': ['id', 'name', 'predicted', 'members (index or list of alternatives)', 'n', 'methylated', 'fold', 'p', 'q'],
                'iface': ['site', 'mask', 'partners (-1 RNA, -2 DNA, else protein index)', 'pdb ids (first 8)', 'pdb count'],
                'enrich': '[methylated at, methylated not at, other R at, other R not at, log2 OR, p] over structurally observed arginines',
                'go_enrich': '[methylated with term, methylated reviewed, other with term, other reviewed, log2 OR, p]',
                'enrich_strata': '{type: {all|ordered|domain: {pooled: [a, b, c, d, log2 OR, p], mh: [log2 MH OR, CMH p, informative proteins]}}}',
                'contacts': ['site', '[partner (-1 RNA, -2 DNA, -3 ligand, -4 unmapped chain, else protein index), pdb, residue labels]', 'methylated residue modelled (2MR/DA2/NMM/AGM)'],
                'cd_members': ['protein', 'condensate index (meta.condensates of the sites file)', 'confidence 1-5', 'role bits', 'evidence bits', 'publications', 'all PMIDs'],
                'wr': ['PRMT or Tudor protein index', 'its STRING physical partners (indices)'],
                'rg_iface': '{segment|dipeptide|other: {cover: [methylated sites, structurally observed], rna|any: {pooled, mh}}}',
                'geometry': ['site', 'pdb', 'chain', 'Levy class (0 none, 1 support, 2 rim, 3 core)', 'rASA alone', 'rASA in complex', 'nucleic-acid contact (base/backbone/contact/none)', 'bases hydrogen-bonded', 'ASA buried by nucleic acid', 'ASA buried by protein', 'hydrogen bonds to protein'],
                'geometry_summary': '{levy_counts, na_counts: {1 methylated, 0 control}, core|base: {pooled, mh}}',
                'modelled': ['site', 'PDB entries modelling a methylated arginine there', 'residue names'],
                'modelled_gaps': ['protein index', 'UniProt positions modelled as methylarginine but not in the union', 'PDB entries'],
                'cond_conf': '{min stars: [[condensate index, members, methylated, log2 OR, p, q], ...]}',
            },
        },
        'proteins': [[a, gene_of.get(base_acc(a), ''), flags[k], ensp(a)] for k, a in enumerate(vocab)],
        'string': string_rows,
        'complexes': cx_rows,
        'iface': site_rows,
        'enrich': enrich_out,
        'go_enrich': go_enrich,
        'enrich_strata': enrich_strata,
        'contacts': contact_rows,
        'cd_members': cd_rows,
        'wr': wr_rows,
        'rg_iface': rg_out,
        'geometry': geo_sites,
        'geometry_summary': geo_summary,
        'modelled': modelled_rows,
        'modelled_gaps': modelled_gaps,
        'cond_conf': cond_conf,
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

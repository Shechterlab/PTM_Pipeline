#!/usr/bin/env python3
"""Export interaction context for the Site Explorer to docs/data/methylarg_interactions.json.

Run after 36_export_site_explorer_data.py (it reads docs/data/methylarg_sites.json for the protein and
site lists). Public sources are downloaded once into --cache (default PTM_data/interactions_cache, which Git
ignores) and reused; pass --refresh to download them again.

Sources
  STRING v12.5 human physical links (CC BY 4.0) .......... physical partners and scores
  Complex Portal human complexes, curated + predicted (CC0) complex membership
  UniProt reviewed human proteome (CC BY 4.0) .............. gene symbols, STRING ids, sequences, Tudor, WD40 and PHD domains
  QuickGO / GOA (CC BY 4.0) ................................ RNA binding (GO:0003723) and DNA binding (GO:0003677)
  PDBe-KB interface residues (EMBL-EBI) .................... residues at protein, RNA and DNA interfaces in the PDB
  SIFTS observed UniProt segments (EMBL-EBI) ............... structural coverage, for the interface enrichment test
  UniProt disordered regions and domains (CC BY 4.0) ....... strata for the interface test
  CD-CODE condensate membership (cd-code.org API) .......... per-protein confidence score, role and evidence
  RCSB ModelServer + SIFTS chain mapping ................... residues each interface methylarginine touches
  ProtVar predicted interfaces (EMBL-EBI, CC BY 4.0) ....... AlphaFold-Multimer interface residues of human pairs
  IntAct features (EMBL-EBI, CC BY 4.0) .................... mutations that change an interaction, binding regions, methylated residues

The protein-set enrichments (Gene Ontology RNA and DNA binding, complexes, condensates by CD-CODE confidence) are
over-representation tests run with clusterProfiler's enricher (scripts/enrich_sets.R, so Rscript must be on PATH).

It also writes docs/data/methylarg_pred_iface.json (--pred-iface-out): every residue of an explorer protein at a ProtVar
predicted interface, as runs with the number of partners and the best pDockQ, which the explorer loads only when a
figure draws the predicted lane.

Stages that make many requests cache every answer: the PDBe interface lookup (about 5,000 calls) and the
structure contacts (about 5,500 calls). Run them alone first with --fetch-interfaces-only / --fetch-contacts-only.
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
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps
from statsmodels.stats.contingency_tables import StratifiedTable

from enrichment import Job, run as run_enrichment
from interaction_inputs import merge_member

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
# ProtVar predicted interfaces: AlphaFold-Multimer models of human binary interactions, residues within 5 Å
# of the partner chain, in UniProt numbering (CC BY 4.0).
PROTVAR_IFACE = 'https://ftp.ebi.ac.uk/pub/databases/ProtVar/predictions/interfaces/2024.05.28_interface_summary_5A.tsv.gz'
PRED_TIERS = (0.23, 0.5)
# Models served by the ProtVar API (/prediction/interaction/{A}/{B}/model) are a later prediction run, listed with
# 10 Å interfaces in this table; the 2024 5 Å table above was computed on earlier models that the API often lacks.
PROTVAR_SERVED = 'https://ftp.ebi.ac.uk/pub/databases/ProtVar/predictions/interfaces/af2complexes_interfaces_v1.tsv.gz'
PROTVAR_MODEL = 'https://www.ebi.ac.uk/ProtVar/api/prediction/interaction/{a}/{b}/model'
AFDB_COMPLEX = 'https://alphafold.ebi.ac.uk/api/complex/{acc}'
AFDB_PREDICTION = 'https://alphafold.ebi.ac.uk/api/prediction/{id}'
# Predicted-contact gallery: every methylarginine within 5 Å of the partner in a served model is measured.
GALLERY_PDOCKQ = 0.5          # ProtVar models measured (served-model pDockQ floor)
SALT_BRIDGE, HBOND_PRED, CATION_PI, CLASH, CONTACT_KEEP = 4.0, 3.5, 5.0, 2.2, 5.0
G_SB, G_HB, G_CATPI, G_CLASH = 1, 2, 4, 8
ARG_N = ('NE', 'NH1', 'NH2')
ACIDIC_O = {('ASP', 'OD1'), ('ASP', 'OD2'), ('GLU', 'OE1'), ('GLU', 'OE2')}
# Residue-level nucleic-acid crosslinks (mass spectrometry, human cells).
# pepR-MS supplement (Sha, Kuster & Trendel, bioRxiv 2025.09.05.674461, CC BY-NC): Table S2 re-analyses Kramer 2014,
# XRNAX (Trendel 2019), RBS-ID (Bae 2020) and pRBS-ID (Bae 2021) with the pepR-MS pipeline; Table S3 lists focused sites.
# bioRxiv rate-limits scripts: if the download fails, save the supplementary zip as crosslinks/pepr_media-1.zip.
PEPR_ZIP = 'https://www.biorxiv.org/content/biorxiv/early/2026/09/21/2025.09.05.674461/DC1/embed/media-1.zip?download=true'
TRENDEL_MMC = 'https://ars.els-cdn.com/content/image/1-s2.0-S0092867425005070-mmc{k}.xlsx'   # Trendel et al. 2025 Cell
# RBDmap (Castello et al. 2016 Mol Cell, HeLa): RNA-bound proteolytic fragments, Oxford Research Archive dataset.
RBDMAP_S1 = 'https://ora.ox.ac.uk/objects/uuid:092f4366-cac9-4ed7-8283-89dc6032834a/files/ma75a36ddd2e91c194d65041242e5a683'
UNIPROT_DNA_BIND = 'https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:9606+AND+ft_dna_bind:*&fields=accession,ft_dna_bind&format=tsv'
RNA_XL_SETS = [('Kramer_2014', 'Kramer 2014'), ('Trendel_2019', 'XRNAX 2019'), ('Bae_2020', 'RBS-ID 2020'), ('Bae_2021', 'pRBS-ID 2021'), ('this study', 'pepR-MS 2025')]
RNA_FOCUSED = 32                       # pepR-MS Table S3 focused crosslinking site
D_XL, D_PEP, D_BIND, D_CAND = 1, 2, 4, 8  # DNA: crosslink localised to this residue alone, inside a crosslinked peptide,
                                          # UniProt DNA-binding region, one of several candidate residues of a crosslink
XL_NEAR = 3                            # "near" = within this many residues of a crosslinked residue
RINGS = {'PHE': ('CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'), 'TYR': ('CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'), 'TRP': ('CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2')}
PDBE_MODRES = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/modified_AA_or_NA/{pdb}'
# IntAct curated features of interaction participants (EMBL-EBI, CC BY 4.0): mutations and their effect on the
# interaction, binding regions mapped for a partner, and modified residues of the participants.
INTACT_FEATURES = 'https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/features/{name}.tsv'
INTACT_RELEASE = 'https://ftp.ebi.ac.uk/pub/databases/intact/current/various/release_time.txt'
INTACT_EFFECTS = ['no effect', 'disrupts', 'decreases', 'increases', 'causes', 'unspecified']   # mutation effect codes
INTACT_EFFECT_MI = {'MI:2226': 0, 'MI:0573': 1, 'MI:1128': 1, 'MI:1129': 1, 'MI:0119': 2, 'MI:1133': 2, 'MI:1130': 2,
                    'MI:0382': 3, 'MI:1132': 3, 'MI:1131': 3, 'MI:2227': 4}
INTACT_REGIONS = ['sufficient', 'necessary', 'binding-associated', 'direct']   # binding region type codes
INTACT_REGION_MI = {'MI:0442': 0, 'MI:0429': 1, 'MI:0117': 2, 'MI:1125': 3}
INTACT_PEPTIDE_MAX = 20     # a binding region at most this long (peptide-level mapping) counts as evidence at its residues
INTACT_MUT_MAX = 3          # mutations spanning at most this many residues are kept per residue
INTACT_PARTNER = re.compile(r'\(([A-Za-z_]+):([^()|]+?)\(psi-mi:"(MI:\d+)"\(([^)]*)\)\)')

# The nine human protein arginine methyltransferases.
PRMT_GENES = ['PRMT1', 'PRMT2', 'PRMT3', 'CARM1', 'PRMT5', 'PRMT6', 'PRMT7', 'PRMT8', 'PRMT9']
# Reader-domain classes marked on partners (InterPro): Tudor, WD40 repeat, PHD finger (SMART and Pfam/PROSITE entries).
DOMAIN_CLASSES = {'tudor': ('IPR002999',), 'wd40': ('IPR001680',), 'phd': ('IPR001965', 'IPR019787')}
GO_TERMS = {'rna': 'GO:0003723', 'dna': 'GO:0003677'}
# GO evidence codes counted as experimental support.
GO_EXPERIMENTAL = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP', 'HGI', 'HEP'}

# Partner flag bits (shared with docs/assets/explorer.js).
F_PRMT, F_TUDOR, F_RNA, F_DNA, F_RNA_EXP, F_DNA_EXP, F_WD40, F_PHD = 1, 2, 4, 8, 16, 32, 64, 128
F_DOMAIN = {'tudor': F_TUDOR, 'wd40': F_WD40, 'phd': F_PHD}
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
# AGM (5-methylarginine) is C-methylated, not a guanidino N-methylarginine, so it is not counted.
METHYL_R = {'2MR', 'DA2', 'NMM'}
AA3 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L',
       'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V', 'SEP': 'pS', 'TPO': 'pT',
       'PTR': 'pY', 'MLY': 'Kme2', 'M3L': 'Kme3', 'MLZ': 'Kme1', 'ALY': 'Kac', '2MR': 'Rme2', 'DA2': 'Rme2a', 'NMM': 'Rme1'}
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


def read_fasta(path: Path) -> dict:
    """UniProt FASTA -> {accession: sequence}; the accession is the second |-field of the header."""
    out, acc = {}, None
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('>'):
                head = line[1:].split()[0]
                acc = head.split('|')[1] if head.count('|') >= 2 else head
                out[acc] = []
            elif acc:
                out[acc].append(line.strip())
    return {a: ''.join(v) for a, v in out.items()}


def open_text(path: Path):
    return io.TextIOWrapper(gzip.open(path, 'rb'), encoding='utf-8') if path.suffix == '.gz' else open(path, encoding='utf-8')


def base_acc(acc: str) -> str:
    return acc.split('-')[0]


# ----------------------------------------------------------------------------- sources
def load_uniprot(cache: Path, refresh: bool):
    """Reviewed human proteome: accession -> (gene, STRING id, sequence); plus {class: accessions} for DOMAIN_CLASSES."""
    path = fetch(UNIPROT_STREAM, cache / 'uniprot_human_reviewed.tsv.gz', refresh,
                 {'query': '(organism_id:9606) AND (reviewed:true)', 'fields': 'accession,gene_primary,xref_string,sequence',
                  'format': 'tsv', 'compressed': 'true'})
    proteome = {}
    with open_text(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            sid = (row.get('STRING') or '').split(';')[0].strip()
            proteome[row['Entry']] = (row.get('Gene Names (primary)') or '', sid, row.get('Sequence') or '')
    classes = {}
    for key, entries in DOMAIN_CLASSES.items():
        xref = ' OR '.join(f'xref:interpro-{e}' for e in entries)
        path = fetch(UNIPROT_STREAM, cache / f'uniprot_human_{key}.tsv', refresh,
                     {'query': f'(organism_id:9606) AND (reviewed:true) AND ({xref})', 'fields': 'accession', 'format': 'tsv'})
        classes[key] = {line.split('\t')[0].strip() for line in open(path, encoding='utf-8').read().splitlines()[1:] if line.strip()}
    log(f'UniProt: {len(proteome):,} reviewed human entries; ' + ', '.join(f'{len(v)} {k}' for k, v in classes.items()))
    return proteome, classes


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
        if accs and accs not in slots:
            # A protein listed twice (two PRO chains of one precursor, or a repeated accession) is one member.
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
        failed.append(acc)

    failed = []
    todo = [a for a in accs if refresh or not (folder / f'{a}.json').exists()]
    log(f'PDBe interfaces: {len(accs) - len(todo):,} cached, fetching {len(todo):,}')
    done = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for _ in pool.map(one, todo):
            done += 1
            if done % 250 == 0:
                log(f'  {done:,}/{len(todo):,}')
    if failed:   # stop rather than leave these proteins' interfaces out
        raise SystemExit(f'PDBe-KB did not answer for {len(failed)} proteins after 5 tries (e.g. {failed[0]}); run the build again')
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
            merge_member(members, acc, (r.get('confidence_score') or 0, roles, evid, pmids))
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
    if text is None:   # the server failed every try: stop rather than leave the evidence out
        raise SystemExit(f'ModelServer did not answer for {pdb} chain {chain} residue {label} after 4 tries; run the build again')
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
    if text is None:   # the server failed every try: stop rather than leave the evidence out
        raise SystemExit(f'ModelServer did not answer for {pdb} chain {chain} residue {label} after 4 tries; run the build again')
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
    """PDB entries that model a methylated arginine (2MR, DA2, NMM), with chain, SEQRES index and residue name."""
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


def cdcode_threshold_enrichment(cd, proteome, meth_base, uids, workdir, thresholds=(1, 2, 3, 4)):
    """Methylated share of each condensate counting only members at or above a CD-CODE confidence score, against the
    proteins in any human CD-CODE condensate (at any score), the background scripts/12 uses for all memberships:
    clusterProfiler enricher, one job per score, condensates with at least five counted members.
    {score: [[condensate index, members, methylated members, log2 fold enrichment, p, q], ...]}"""
    universe = {a for uid in uids for a in cd.get(uid, {}) if a in proteome}
    target = {a for a in universe if a in meth_base}
    jobs = []
    for t in thresholds:
        sets = pd.DataFrame([(ci, a) for ci, uid in enumerate(uids) for a, m in cd.get(uid, {}).items() if m[0] >= t and a in proteome],
                            columns=['term', 'gene'])
        jobs.append(Job(f'condensate_confidence_{t}', target, universe, sets=sets, whole_universe=True, min_size=5))
    res = run_enrichment(jobs, workdir)
    out = {}
    for t in thresholds:
        r = res[f'condensate_confidence_{t}']
        rows = [[int(x.ID), int(x.M), int(x.k), sig(math.log2(x.FoldEnrichment)), sig(x.pvalue), sig(x['p.adjust'])] for _, x in r.iterrows()]
        out[str(t)] = sorted(rows)
    return out


def load_protvar(cache: Path, refresh: bool, keep: set):
    """Predicted heterodimer interfaces touching a protein in keep: acc -> list of (partner, pDockQ, R positions, all positions, is_first)."""
    path = fetch(PROTVAR_IFACE, cache / Path(PROTVAR_IFACE).name, refresh)
    out = defaultdict(list)
    n = 0
    with open_text(path) as fh:
        rows = csv.reader(fh, delimiter='\t')
        header = next(rows)
        ix = {h: k for k, h in enumerate(header)}
        for row in rows:
            try:
                score = float(row[ix['pdockq']])
            except (ValueError, IndexError):
                continue
            a, b = base_acc(row[ix['uniprot_id1']]), base_acc(row[ix['uniprot_id2']])
            if a == b or score < PRED_TIERS[0] or (a not in keep and b not in keep):
                continue
            for acc, other, field, first in ((a, b, 'ifresid1', True), (b, a, 'ifresid2', False)):
                if acc not in keep:
                    continue
                res = [x for x in (row[ix[field]] or '').split(',') if x]
                pos_all = {int(x[1:]) for x in res if x[1:].isdigit()}
                pos_r = {int(x[1:]) for x in res if x[:1] == 'R' and x[1:].isdigit()}
                if pos_all:
                    out[acc].append((other, score, pos_r, pos_all, first))
                    n += 1
    log(f'ProtVar predicted interfaces: {n:,} heterodimer interface sides (pDockQ >= {PRED_TIERS[0]}) on {len(out):,} methylated proteins')
    return out


def predicted_interface_runs(pred: dict) -> dict:
    """acc -> [[start, end, partners, best pDockQ x100], ...]: runs of consecutive residues at the predicted interface
    with any partner (pDockQ at or above the lowest tier), with the most partners and the best pDockQ along the run."""
    out = {}
    for acc, sides in pred.items():
        partners, best = defaultdict(set), {}
        for other, score, _, pos_all, _ in sides:
            for pos in pos_all:
                partners[pos].add(other)
                best[pos] = max(best.get(pos, 0.0), score)
        runs = []
        for pos in sorted(partners):
            n, q = len(partners[pos]), round(best[pos] * 100)
            if runs and pos == runs[-1][1] + 1:
                runs[-1][1:] = [pos, max(runs[-1][2], n), max(runs[-1][3], q)]
            else:
                runs.append([pos, pos, n, q])
        out[acc] = runs
    return out


def load_protvar_served(cache: Path, refresh: bool):
    """Ordered (uniprot_id1, uniprot_id2) pairs that the ProtVar API serves a model for."""
    path = fetch(PROTVAR_SERVED, cache / Path(PROTVAR_SERVED).name, refresh)
    out = set()
    with open_text(path) as fh:
        rows = csv.reader(fh, delimiter='\t')
        ix = {h: k for k, h in enumerate(next(rows))}
        for row in rows:
            out.add((row[ix['uniprot_id1']], row[ix['uniprot_id2']]))
    log(f'ProtVar served models: {len(out):,} pairs')
    return out


def get_bytes(url: str, tries: int = 4):
    """Body of a GET (gunzipped when needed), or None on 404."""
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': UA, 'Accept-Encoding': 'gzip'})
            with urllib.request.urlopen(req, timeout=180) as r:
                data = r.read()
            return gzip.decompress(data) if data[:2] == b'\x1f\x8b' else data
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if k == tries - 1:
                raise
        except Exception:  # noqa: BLE001 - retried, then raised
            if k == tries - 1:
                raise
        time.sleep(2 ** (k + 1))
    return None


def parse_model_chains(text: str):
    """Heavy atoms of the first model: chain -> residue number -> [(resn, atom, x, y, z, b)]."""
    chains = defaultdict(lambda: defaultdict(list))
    for line in text.splitlines():
        if line.startswith('ENDMDL'):
            break
        if not line.startswith('ATOM') or line[16] not in ' A':
            continue
        name = line[12:16].strip()
        elem = line[76:78].strip() or name[:1]
        if elem == 'H':
            continue
        chains[line[21]][int(line[22:26])].append((line[17:20], name, float(line[30:38]), float(line[38:46]), float(line[46:54]), float(line[60:66])))
    return chains


def measure_methyl_r(chains, chain: str, positions):
    """Arginines of `chain` at `positions` against every other chain of a predicted complex.

    Returns pos -> [min distance, flags, pLDDT of the arginine, lowest pLDDT of the partner residues named,
    [partner chain, residue label with bond codes]] for arginines within CONTACT_KEEP of a partner atom."""
    import numpy as np
    me = chains.get(chain) or {}
    meta, xyz = [], []
    rings = []
    ca_b = {}
    for c, residues in chains.items():
        if c == chain:
            continue
        for resi, atoms in residues.items():
            for resn, name, x, y, z, b in atoms:
                meta.append((c, resi, resn, name))
                xyz.append((x, y, z))
                if name == 'CA':
                    ca_b[(c, resi)] = b
            ring = RINGS.get(atoms[0][0])
            if ring:
                pts = [(x, y, z) for resn, name, x, y, z, b in atoms if name in ring]
                if len(pts) == 6:
                    rings.append((c, resi, atoms[0][0], np.mean(np.array(pts), axis=0)))
    if not xyz:
        return {}
    P = np.array(xyz)
    out = {}
    for pos in positions:
        atoms = me.get(pos)
        if not atoms or atoms[0][0] != 'ARG':
            continue
        X = np.array([(x, y, z) for _, _, x, y, z, _ in atoms])
        D = np.sqrt(((X[:, None, :] - P[None, :, :]) ** 2).sum(-1))
        dmin = float(D.min())
        if dmin > CONTACT_KEEP:
            continue
        names = [a[1] for a in atoms]
        flags = G_CLASH if dmin < CLASH else 0
        bonds = defaultdict(set)
        for ai, an in enumerate(names):
            if an not in ARG_N:
                continue
            for pj in np.nonzero(D[ai] <= HBOND_PRED)[0]:
                c, resi, resn, pn = meta[pj]
                if pn.startswith('O') or (resn == 'HIS' and pn in ('ND1', 'NE2')):
                    bonds[(c, resi, resn)].add('hb')
            for pj in np.nonzero(D[ai] <= SALT_BRIDGE)[0]:
                c, resi, resn, pn = meta[pj]
                if (resn, pn) in ACIDIC_O:
                    bonds[(c, resi, resn)].add('sb')
        if 'CZ' in names:
            cz = X[names.index('CZ')]
            for c, resi, resn, centre in rings:
                if float(np.sqrt(((cz - centre) ** 2).sum())) <= CATION_PI:
                    bonds[(c, resi, resn)].add('pi')
        kinds = set().union(*bonds.values()) if bonds else set()
        flags |= (G_SB if 'sb' in kinds else 0) | (G_HB if 'hb' in kinds else 0) | (G_CATPI if 'pi' in kinds else 0)
        if not bonds:
            c, resi, resn, _ = meta[int(D.min(axis=0).argmin())]
            bonds[(c, resi, resn)] = set()
        ca = next((b for _, n, *_rest, b in atoms if n == 'CA'), atoms[0][5])
        partner_b = [ca_b.get((c, resi)) for c, resi, _ in bonds]
        labels = [[c, AA3.get(resn, 'X') + str(resi) + (':' + '+'.join(sorted(k)) if k else '')] for (c, resi, resn), k in sorted(bonds.items(), key=lambda kv: kv[0][1])]
        out[pos] = [round(dmin, 1), flags, round(ca), round(min(b for b in partner_b if b is not None)) if any(b is not None for b in partner_b) else None, labels[:6]]
    return out


def measure_protvar_models(cache: Path, refresh: bool, workers: int, sites_by_acc: dict):
    """Methylarginines at the interface of served ProtVar models (pDockQ >= GALLERY_PDOCKQ): one cached JSON per pair."""
    path = fetch(PROTVAR_SERVED, cache / Path(PROTVAR_SERVED).name, refresh)
    jobs = []
    with open_text(path) as fh:
        rows = csv.reader(fh, delimiter='\t')
        ix = {h: k for k, h in enumerate(next(rows))}
        for row in rows:
            a, b = row[ix['uniprot_id1']], row[ix['uniprot_id2']]
            if '-' in a or '-' in b or a == b:
                continue
            try:
                sc = float(row[ix['pdockq']])
            except ValueError:
                continue
            if sc < GALLERY_PDOCKQ:
                continue
            want = {}
            for acc, ch, col in ((a, row[ix['chain1']] or 'A', 'ifresid1'), (b, row[ix['chain2']] or 'B', 'ifresid2')):
                near = sites_by_acc.get(acc, set()) & {int(v) for v in row[ix[col]].split(',') if v.isdigit()}
                if near:
                    want[ch] = (acc, sorted(near))
            if want:
                jobs.append((a, b, sc, want))
    folder = cache / 'protvar_measure'
    folder.mkdir(parents=True, exist_ok=True)
    done = [0]

    def one(job):
        a, b, sc, want = job
        dest = folder / f'{a}_{b}.json'
        if dest.exists() and not refresh:
            return
        data = get_bytes(PROTVAR_MODEL.format(a=a, b=b))
        res = {}
        if data:
            chains = parse_model_chains(data.decode('utf-8', 'replace'))
            for ch, (acc, pos) in want.items():
                m = measure_methyl_r(chains, ch, pos)
                if m:
                    res[ch] = {str(k): v for k, v in m.items()}
        dest.write_text(json.dumps({'served': bool(data), 'chains': res}))
        done[0] += 1
        if done[0] % 250 == 0:
            log(f'  ProtVar models measured: {done[0]:,}')

    todo = [j for j in jobs if refresh or not (folder / f'{j[0]}_{j[1]}.json').exists()]
    log(f'ProtVar served models (pDockQ >= {GALLERY_PDOCKQ}) with a methylarginine within 10 Å: {len(jobs):,} pairs, {len(todo):,} to measure')
    with ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(one, todo))
    out = []
    for a, b, sc, want in jobs:
        f = folder / f'{a}_{b}.json'
        if not f.exists():
            continue
        rec = json.loads(f.read_text())
        for ch, sites in rec.get('chains', {}).items():
            acc = want[ch][0]
            other = b if acc == a else a
            for pos, m in sites.items():
                out.append({'acc': acc, 'pos': int(pos), 'partner': other, 'src': 0, 'score': sc, 'pdockq': sc, 'm': m, 'ref': 1 if acc == a else 2, 'chain': ch})
    return out


def measure_afdb_complexes(cache: Path, refresh: bool, workers: int, sites_by_acc: dict):
    """Human heterodimers in AlphaFold DB that include a methylated protein; every methylarginine is measured."""
    lists = cache / 'afdb_complex'
    preds = cache / 'afdb_prediction'
    folder = cache / 'afdb_measure'
    for d in (lists, preds, folder):
        d.mkdir(parents=True, exist_ok=True)
    accs = sorted(sites_by_acc)

    def get_list(acc):
        dest = lists / f'{acc}.json'
        if dest.exists() and not refresh:
            return
        data = get_bytes(AFDB_COMPLEX.format(acc=acc))
        dest.write_text(data.decode() if data else '[]')

    with ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(get_list, accs))
    complexes = {}
    for acc in accs:
        try:
            rows = json.loads((lists / f'{acc}.json').read_text())
        except (OSError, json.JSONDecodeError):
            continue
        for x in rows if isinstance(rows, list) else []:
            if x.get('assemblyType') == 'Hetero' and set(x.get('taxId') or []) == {9606} and x.get('modelEntityId'):
                complexes[x['modelEntityId']] = x
    log(f'AlphaFold DB: {len(complexes):,} human heterodimer models with a methylated protein')
    done = [0]

    def one(mid):
        dest = folder / f'{mid}.json'
        if dest.exists() and not refresh:
            return
        pf = preds / f'{mid}.json'
        if not pf.exists() or refresh:
            data = get_bytes(AFDB_PREDICTION.format(id=mid))
            pf.write_text(data.decode() if data else '[]')
        chains_meta = json.loads(pf.read_text())
        res = {}
        mine = [c for c in chains_meta if c.get('uniprotAccession') in sites_by_acc]
        if mine and chains_meta[0].get('pdbUrl'):
            data = get_bytes(chains_meta[0]['pdbUrl'])
            if data:
                chains = parse_model_chains(data.decode('utf-8', 'replace'))
                for c in mine:
                    m = measure_methyl_r(chains, c['chainId'], sorted(sites_by_acc[c['uniprotAccession']]))
                    if m:
                        res[c['chainId']] = {str(k): v for k, v in m.items()}
        dest.write_text(json.dumps({'chains': res, 'meta': [[c.get('chainId'), c.get('uniprotAccession')] for c in chains_meta]}))
        done[0] += 1
        if done[0] % 250 == 0:
            log(f'  AlphaFold DB models measured: {done[0]:,}')

    with ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(one, sorted(complexes)))
    out = []
    for mid, x in complexes.items():
        f = folder / f'{mid}.json'
        if not f.exists():
            continue
        rec = json.loads(f.read_text())
        acc_of = {c: a for c, a in rec.get('meta', [])}
        for ch, sites in rec.get('chains', {}).items():
            acc = acc_of.get(ch)
            others = [a for c, a in rec.get('meta', []) if c != ch]
            if not acc or not others:
                continue
            for pos, m in sites.items():
                out.append({'acc': acc, 'pos': int(pos), 'partner': others[0], 'src': 1, 'score': x.get('complexPredictionAccuracy_ipTM'),
                            'pdockq': x.get('complexPredictionAccuracy_pDockQ'), 'm': m, 'ref': mid, 'chain': ch})
    return out


def load_crosslink_sequences(seq: dict, proteome: dict) -> dict:
    """The crosslinked peptides and RBDmap fragments, for display; read from the frozen tables (use_frozen_tables)."""
    log('crosslink peptides are shown only when the frozen tables are read (--tables)')
    return {'rna': {}, 'dna': {}, 'rbd': {}}


def load_crosslinks(cache: Path, refresh: bool, seq: dict, seq_meta: dict | None = None):
    """Residue-level RNA and DNA crosslinks in canonical UniProt numbering, each checked against the sequence.

    Returns rna[acc][pos] -> bitmask (one bit per dataset, RNA_FOCUSED for pepR-MS focused sites), rbd[acc] ->
    [(start, end)] RBDmap RNA-bound fragments, dna[acc][pos] -> D_XL, dna_peps[acc] -> [(start, end)],
    dna_interactome (set of accessions), dna_bind[acc] -> [(start, end)], and the sequence check {matched, dropped}."""
    import openpyxl
    import zipfile
    folder = cache / 'crosslinks'
    folder.mkdir(parents=True, exist_ok=True)
    rna, dna, rbd = defaultdict(dict), defaultdict(dict), defaultdict(list)
    rna_n = defaultdict(Counter)   # crosslinked peptide-spectrum rows per residue (Table S2, all datasets)
    # A crosslink (dataset, peptide, RNA mass) is often listed at several candidate positions; residues named alone
    # by at least one crosslink are "single-localisation" residues.
    rna_cands, rna_single = defaultdict(set), defaultdict(set)
    dna_peps, dna_bind, dna_interactome = defaultdict(list), defaultdict(list), set()
    checked = Counter()

    def ok(acc, pos, aa):
        s = seq.get(acc)
        good = bool(s) and 0 < pos <= len(s) and s[pos - 1] == aa
        checked['ok' if good else 'mismatch'] += 1
        return good

    zpath = folder / 'pepr_media-1.zip'
    try:
        fetch(PEPR_ZIP, zpath, refresh)
    except Exception as e:  # noqa: BLE001 - bioRxiv often answers scripts with 429
        if not zpath.exists():
            log(f'pepR-MS supplement not available ({e}); save it as {zpath} to include RNA crosslinks')
    if zpath.exists():
        with zipfile.ZipFile(zpath) as zf:
            names = {Path(n).name: n for n in zf.namelist()}
            wb = openpyxl.load_workbook(io.BytesIO(zf.read(names['Table_S2.xlsx'])), read_only=True)
            for bit, (sheet, _) in enumerate(RNA_XL_SETS):
                rows = wb[sheet].iter_rows(values_only=True)
                ix = {h: k for k, h in enumerate(next(rows))}
                for r in rows:
                    acc = base_acc(str(r[ix['Protein ID']]).split(';')[0])
                    try:
                        pos = int(r[ix['Position']])
                    except (TypeError, ValueError):
                        continue
                    # The crosslinked residue is the one carrying the nucleotide mass in the modified sequence.
                    m = re.search(r'([A-Z])\[(\d+\.\d+)\]', re.sub(r'C\[57\.02\d*\]|M\[15\.99\d*\]', lambda x: x.group(0)[0], str(r[ix['Modified Sequence']])))
                    if m and ok(acc, pos, m.group(1)):
                        rna[acc][pos] = rna[acc].get(pos, 0) | (1 << bit)
                        rna_n[acc][pos] += 1
                        rna_cands[(bit, acc, str(r[ix['Peptide Sequence']]), str(r[ix['RNA Mass']]))].add(pos)
            wb = openpyxl.load_workbook(io.BytesIO(zf.read(names['Table_S3.xlsx'])), read_only=True)
            rows = wb['Crosslinking site'].iter_rows(values_only=True)
            ix = {h: k for k, h in enumerate(next(rows))}
            for r in rows:
                acc = base_acc(str(r[ix['Protein ID']]))
                focused = list(zip(json.loads(r[ix['Focused Site AA']] or '[]'), json.loads(r[ix['Focused Site Position']] or '[]')))
                for aa, pos in focused:
                    if ok(acc, int(pos), aa):
                        rna[acc][int(pos)] = rna[acc].get(int(pos), 0) | RNA_FOCUSED
                        if len(focused) == 1:
                            rna_single[acc].add(int(pos))
            for (_, acc, _, _), positions in rna_cands.items():
                if len(positions) == 1:
                    rna_single[acc].update(positions)
            log(f'RNA crosslinks: {sum(len(v) for v in rna.values()):,} residues, {sum(len(v) for v in rna_single.values()):,} named alone by at least one crosslink')

    # RBDmap high-confidence RNA-bound fragments (RBDpep), kept where the fragment matches the canonical sequence.
    rpath = folder / 'rbdmap_castello2016_table_s1.xlsx'
    try:
        fetch(RBDMAP_S1, rpath, refresh)
    except Exception as e:  # noqa: BLE001
        log(f'RBDmap Table S1 not available ({e})')
    if rpath.exists():
        gene_acc = {v[0]: a for a, v in sorted(seq_meta.items())} if seq_meta else {}
        wb = openpyxl.load_workbook(rpath, read_only=True)
        rows = wb['RBDpep'].iter_rows(values_only=True)
        head = next(rows)
        ix = {str(h).strip(): k for k, h in enumerate(head)}
        frag_col = next(k for h, k in ix.items() if h.startswith('LysC/ArgC'))
        kept = moved = dropped = 0
        for r in rows:
            frag = str(r[frag_col] or '')
            acc = base_acc(str(r[ix['ProtID']]))
            if acc not in seq:
                acc = gene_acc.get(str(r[ix['Symbol']]), acc)
            s = seq.get(acc, '')
            try:
                a, b = int(r[ix['fragmentStart']]), int(r[ix['fragmentStop']])
            except (TypeError, ValueError):
                dropped += 1
                continue
            if not s or not frag:
                dropped += 1
                continue
            if s[a - 1:b] != frag:
                hit = s.find(frag)
                if hit < 0 or s.find(frag, hit + 1) >= 0:
                    dropped += 1
                    continue
                a, b = hit + 1, hit + len(frag)
                moved += 1
            rbd[acc].append((a, b))
            kept += 1
        log(f'RBDmap: {kept:,} RNA-bound fragments kept ({moved:,} re-located by sequence), {dropped:,} dropped')

    for k in (1, 2):
        dest = folder / f'trendel2025_mmc{k}.xlsx'
        try:
            fetch(TRENDEL_MMC.format(k=k), dest, refresh)
        except Exception as e:  # noqa: BLE001
            log(f'Trendel 2025 table S{k} not available ({e})')
    if (folder / 'trendel2025_mmc2.xlsx').exists():
        wb = openpyxl.load_workbook(folder / 'trendel2025_mmc2.xlsx', read_only=True)
        rows = wb['4ST-xl-peptides localization'].iter_rows(values_only=True)
        ix = {h: k for k, h in enumerate(next(rows))}
        for r in rows:
            pep, acc = str(r[ix['Sequence with localization']]), base_acc(str(r[ix['Uniprot ID']]))
            try:
                start = int(r[ix['Protein position start']])
            except (TypeError, ValueError):
                continue
            # Lower case marks the residues the crosslink may sit on; most peptides leave several candidates.
            bit = D_XL if sum(ch.islower() for ch in pep) == 1 else D_CAND
            for k, ch in enumerate(pep):
                if ch.islower() and ok(acc, start + k, ch.upper()):
                    dna[acc][start + k] = dna[acc].get(start + k, 0) | bit
        rows = wb['4ST-xl-peptides combined'].iter_rows(values_only=True)
        ix = {h: k for k, h in enumerate(next(rows))}
        for r in rows:
            try:
                dna_peps[base_acc(str(r[ix['Uniprot ID']]))].append((int(r[ix['Protein position start']]), int(r[ix['Protein position end']])))
            except (TypeError, ValueError):
                continue
    if (folder / 'trendel2025_mmc1.xlsx').exists():
        wb = openpyxl.load_workbook(folder / 'trendel2025_mmc1.xlsx', read_only=True)
        rows = wb['DNA interactome & subproteomes'].iter_rows(values_only=True)
        ix = {h: k for k, h in enumerate(next(rows))}
        for r in rows:
            if r[ix['MCF7 DNA interactome']] == 'DNA interactome':
                dna_interactome.add(base_acc(str(r[ix['Uniprot ID']])))
    path = fetch(UNIPROT_DNA_BIND, folder / 'uniprot_dna_bind.tsv', refresh)
    with open(path, encoding='utf-8') as fh:
        next(fh)
        for line in fh:
            acc, _, feats = line.rstrip('\n').partition('\t')
            for a, b in re.findall(r'DNA_BIND (\d+)\.\.(\d+)', feats):
                dna_bind[acc].append((int(a), int(b)))
    log(f'Crosslinks: RNA {sum(len(v) for v in rna.values()):,} residues in {len(rna):,} proteins; DNA {sum(len(v) for v in dna.values()):,} residues in {len(dna):,} proteins, '
        f'{len(dna_interactome):,} DNA-interactome proteins; UniProt DNA-binding regions in {len(dna_bind):,} proteins; sequence check {dict(checked)}')
    return rna, rna_n, rbd, dna, dna_peps, dna_interactome, dna_bind, rna_single, {'matched': checked['ok'], 'dropped': checked['mismatch']}


def load_intact(cache: Path, refresh: bool, seq: dict):
    """IntAct features on reviewed human proteins whose curated sequence matches the canonical sequence at their range.

    Returns ({accession: [(kind, code, start, end, resulting, partners, pmid, interaction AC)]}, counts, release) with kind
    0 mutation (code: INTACT_EFFECTS), 1 binding region (INTACT_REGIONS), 2 methylated arginine (code 0); partners are
    (database, identifier, participant type) of the other participants, or the protein itself for a homomer."""
    feats, counts = defaultdict(list), Counter()
    for kind, name in ((0, 'mutations'), (1, 'bindings_regions'), (2, 'ptms')):
        path = fetch(INTACT_FEATURES.format(name=name), cache / 'intact' / f'{name}.tsv', refresh)
        with open(path, encoding='utf-8') as fh:
            next(fh)
            for line in fh:
                x = line.rstrip('\n').split('\t')
                if len(x) < 16 or 'taxid:9606' not in x[10] or not x[7].startswith('uniprotkb:'):
                    continue
                acc = x[7][len('uniprotkb:'):]
                if kind == 2:
                    label = x[5].lower()
                    if 'arginine' not in label or 'methyl' not in label:
                        continue
                    code = 0
                else:
                    mi = re.search(r'MI:\d+', x[5])
                    code = (INTACT_EFFECT_MI.get(mi.group(0), 5) if kind == 0 else INTACT_REGION_MI.get(mi.group(0))) if mi else None
                    if code is None:
                        counts[f'{name}:other_type'] += 1
                        continue
                s = seq.get(acc)
                if s is None:
                    counts[f'{name}:' + ('isoform_or_chain' if '-' in acc else 'not_reviewed_human')] += 1
                    continue
                ranges = x[2].split(',')
                spans = [tuple(map(int, r.split('-'))) for r in ranges if re.fullmatch(r'\d+-\d+', r)]
                if len(spans) != len(ranges):
                    counts[f'{name}:undetermined_range'] += 1
                    continue
                original = x[3].split(',') if x[3] not in ('', '-') else []
                if any(b > len(s) or a < 1 or (k < len(original) and s[a - 1:b] != original[k]) for k, (a, b) in enumerate(spans)):
                    counts[f'{name}:sequence_mismatch'] += 1
                    continue
                if kind == 2 and any(s[a - 1:b] != 'R' for a, b in spans):
                    counts[f'{name}:not_arginine'] += 1
                    continue
                parts = INTACT_PARTNER.findall(x[11])
                own = [p for p in parts if p[0] == 'uniprotkb' and p[1].split('-')[0] == acc]
                others = [p for p in parts if not (p[0] == 'uniprotkb' and p[1].split('-')[0] == acc)]
                partners = tuple((db, ident.strip('"'), ptype) for db, ident, _, ptype in others) if others else (('uniprotkb', acc, 'protein'),) if len(own) > 1 or len(parts) == 1 else ()
                pmid = next((int(m.group(1)) for m in re.finditer(r'pubmed:(\d+)', x[12])), 0)
                resulting = x[4] if kind == 0 and x[4] not in ('', '-') else ''
                for k, (a, b) in enumerate(spans):
                    feats[acc].append((kind, code, a, b, resulting.split(',')[k] if resulting and len(resulting.split(',')) == len(spans) else resulting, partners, pmid, x[14]))
                counts[f'{name}:kept'] += 1
    release = ''
    try:
        release = fetch(INTACT_RELEASE, cache / 'intact' / 'release_time.txt', refresh).read_text(encoding='utf-8').strip()
    except Exception as e:  # noqa: BLE001 - the release note is informative only
        log(f'IntAct release note not available ({e})')
    log('IntAct features: ' + ', '.join(f'{k} {v:,}' for k, v in sorted(counts.items())))
    return feats, dict(counts), release


def informative_strata(tables):
    """Strata that carry information for a Mantel-Haenszel analysis (no empty margin, at least two observations)."""
    return [(a, b, c, d) for a, b, c, d in tables if a + b + c + d >= 2 and (a + b) and (c + d) and (a + c) and (b + d)]


def mantel_haenszel(tables):
    """Mantel-Haenszel common odds ratio and the continuity-corrected CMH test (statsmodels StratifiedTable).

    Returns (odds, p, informative strata); odds may be inf when no stratum has b*c > 0."""
    strata = informative_strata(tables)
    if not strata:
        return None, None, 0
    st = StratifiedTable([np.array([[a, b], [c, d]], dtype=float) for a, b, c, d in strata])
    den = sum(b * c / (a + b + c + d) for a, b, c, d in strata)
    odds = float(st.oddsratio_pooled) if den > 0 else float('inf')
    # statsmodels reports 1 - cdf, which underflows below ~1e-16; take its statistic to scipy's survival function.
    return odds, float(sps.chi2.sf(st.test_null_odds(correction=True).statistic, 1)), len(strata)


def mantel_haenszel_ci(tables):
    """Mantel-Haenszel OR with its Robins-Breslow-Greenland 95% CI (statsmodels): (odds, lo, hi, p, informative strata)."""
    odds, p, used = mantel_haenszel(tables)
    if odds is None or not math.isfinite(odds) or odds <= 0:
        return odds, None, None, p, used
    strata = informative_strata(tables)
    st = StratifiedTable([np.array([[a, b], [c, d]], dtype=float) for a, b, c, d in strata])
    lo, hi = st.oddsratio_pooled_confint(alpha=0.05)
    return odds, (float(lo) if math.isfinite(lo) else None), (float(hi) if math.isfinite(hi) else None), p, used


def woolf_ci(a, b, c, d):
    """95% CI of the sample odds ratio (scipy, normal approximation on the log scale)."""
    with np.errstate(divide='ignore', invalid='ignore'):   # a zero cell gives an open interval, reported as missing
        ci = sps.contingency.odds_ratio([[a, b], [c, d]], kind='sample').confidence_interval(0.95)
    return float(ci.low), float(ci.high)


def or_summary(tables):
    """Pooled Fisher OR and within-stratum MH OR, each with 95% CI and p: {'pooled': [a, b, c, d, OR, lo, hi, p], 'mh': [OR, lo, hi, p, used]}."""
    tables = list(tables)
    a, b, c, d = (sum(x[k] for x in tables) for k in range(4))
    if not (a + b) or not (c + d):
        return {'pooled': [a, b, c, d, None, None, None, None], 'mh': [None, None, None, None, 0], 'proteins': len(tables)}
    pv, _ = fisher_2x2(a, b, c, d)
    odds = a * d / (b * c) if b * c else float('inf')   # sample odds ratio, the centre of the Woolf interval
    lo, hi = woolf_ci(a, b, c, d)
    mo, mlo, mhi, mp, used = mantel_haenszel_ci(tables)
    return {'pooled': [a, b, c, d, sig(odds), sig(lo), sig(hi), sig(pv)], 'mh': [sig(mo), sig(mlo), sig(mhi), sig(mp), used], 'proteins': len(tables)}


def mann_whitney(x, y):
    """Two-sided Mann-Whitney U test (scipy): (U for x, p)."""
    if not len(x) or not len(y):
        return None, None
    res = sps.mannwhitneyu(x, y, alternative='two-sided')
    return float(res.statistic), float(res.pvalue)


# ----------------------------------------------------------------------------- statistics
def fisher_2x2(a, b, c, d):
    """Two-sided Fisher exact p (scipy) for [[a, b], [c, d]] and the Haldane-Anscombe corrected sample odds ratio."""
    p = float(sps.fisher_exact([[a, b], [c, d]]).pvalue)
    return p, ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))


def bh(pvals):
    """Benjamini-Hochberg adjusted p values (scipy)."""
    return [float(q) for q in sps.false_discovery_control(pvals, method='bh')] if len(pvals) else []


def sig(x, d=3):
    """Round to d significant digits; None for missing, NaN or infinite values (JSON has no Infinity)."""
    return None if x is None or (isinstance(x, float) and not math.isfinite(x)) else float(f'{x:.{d}g}')


# ----------------------------------------------------------------------------- main
# ----------------------------------------------------------------------------- frozen tables (--tables)
class Measured:
    """A residue environment already measured (inputs/interactions/pdb_arginine_interface_environments_*): its result."""

    def __init__(self, result):
        self.result = result


def use_frozen_tables(folder: Path, topup: Path, measured_sites: dict):
    """Read the frozen tables in folder (scripts/interaction_inputs.py) instead of the download cache. Anything the
    tables do not cover (proteins never looked up, structure neighbourhoods and model positions never measured) is
    fetched as before, into topup. measured_sites: {accession: positions} measured when the tables were made."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import interaction_inputs as ii
    t = ii.Tables(folder)
    g = globals()
    orig = {k: g[k] for k in ('fetch_interfaces', 'site_contacts', 'residue_environment', 'analyze_environment')}
    manifest = folder.parent / 'MANIFEST.tsv'
    intact_release = ''
    if manifest.exists():
        with open(manifest, encoding='utf-8') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                if r['file'].startswith('interactions/intact_mutations_human_'):
                    intact_release = r['version_or_date']
    g['load_uniprot'] = lambda cache, refresh: ii.uniprot(t)
    g['load_go'] = lambda cache, refresh: ii.go(t)
    g['load_string'] = lambda cache, refresh, sid_to_acc, keep, min_score: ii.string(t, keep, min_score, sid_to_acc)
    g['load_complexes'] = lambda cache, refresh: ii.complexes(t)
    g['load_observed'] = lambda cache, refresh, keep: ii.observed(t, keep)
    g['load_uniprot_regions'] = lambda cache, refresh: ii.regions(t)
    g['load_sifts_chains'] = lambda cache, refresh, keep_pdbs: ii.sifts_chains(t, keep_pdbs)
    g['load_modelled_methyl_r'] = lambda cache, refresh, workers: ii.modelled_methyl_r(t)
    reference = {a: v[2] for a, v in ii.uniprot(t)[0].items()}
    g['load_protvar'] = lambda cache, refresh, keep: ii.protvar(t, keep, PRED_TIERS[0], reference)
    g['load_protvar_served'] = lambda cache, refresh: {(r['uniprot_id1'], r['uniprot_id2']) for r in ii.protvar_served_rows(t)}
    g['load_cdcode'] = lambda cache, refresh, uids, workers: ii.cdcode(t, ROLE_BITS, ROLE_OTHER, EVIDENCE_BITS, EVIDENCE_OTHER)

    def load_intact_(cache, refresh, seq):
        feats, counts = ii.intact(t, INTACT_EFFECT_MI, INTACT_REGION_MI)
        return feats, counts, intact_release
    g['load_intact'] = load_intact_
    g['load_crosslinks'] = lambda cache, refresh, seq, seq_meta=None: ii.crosslinks(t, seq, seq_meta, RNA_XL_SETS, RNA_FOCUSED, D_XL, D_CAND)
    g['load_crosslink_sequences'] = lambda seq, proteome: ii.crosslink_sequences(t, seq, proteome, RNA_XL_SETS)

    queried = ii.queried_interfaces(t)

    def fetch_interfaces_(cache, accs, refresh, workers):
        out = ii.interfaces(t, [a for a in accs if a in queried])
        missing = [a for a in accs if a not in queried]
        if missing:
            log(f'PDBe interfaces: {len(missing):,} proteins not in the frozen tables, fetched into {topup}')
            out.update(orig['fetch_interfaces'](topup, missing, False, workers))
        return out
    g['fetch_interfaces'] = fetch_interfaces_

    contacts = ii.contacts(t)

    def site_contacts_(cache, pdb, chain, label, refresh):
        hit = contacts.get((pdb, chain, int(label)))
        return {'resn': hit['resn'], 'hits': [list(h) for h in hit['hits']]} if hit else orig['site_contacts'](topup, pdb, chain, label, False)
    g['site_contacts'] = site_contacts_

    envs = ii.environments(t)

    def residue_environment_(cache, pdb, chain, label, refresh):
        key = (pdb, chain, int(label))
        return Measured(envs[key]) if key in envs else orig['residue_environment'](topup, pdb, chain, label, False)
    g['residue_environment'] = residue_environment_
    g['analyze_environment'] = lambda rows, chain, label: rows.result if isinstance(rows, Measured) else orig['analyze_environment'](rows, chain, label)

    measured = ii.model_contacts(t, G_CLASH, G_SB, G_HB, G_CATPI, reference)
    measured_at = ii.model_positions(t)
    served_rows = ii.protvar_served_rows(t)
    heterodimers = ii.afdb_heterodimers(t)
    afdb_queried = ii.afdb_queried(t)

    def partner_agrees(m, partner):
        """A measure's partner residues agree with the partner's reference sequence (else the chains are mixed up)."""
        labels = [lab.split(':')[0] for _, lab in m[4]]   # 'E98:hb+sb' -> 'E98'
        return partner in reference and (not labels or ii.labels_agree(labels, reference[partner]))

    def measure_new(folder_name, key, data, want):
        """Measure positions the tables do not have, in a model downloaded now; cached in topup/folder_name/key.json."""
        dest = topup / folder_name / f'{key}.json'
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists():
            return json.loads(dest.read_text()).get('chains', {})
        raw = data() if callable(data) else data
        res = {}
        if raw:
            chains = parse_model_chains(raw.decode('utf-8', 'replace'))
            for ch, (acc, pos) in want.items():
                m = measure_methyl_r(chains, ch, pos)
                if m:
                    res[ch] = {str(k): v for k, v in m.items()}
        dest.write_text(json.dumps({'served': bool(raw), 'chains': res, 'positions': {ch: [acc, pos] for ch, (acc, pos) in want.items()}}))
        return res

    def measure_protvar_models_(cache, refresh, workers, sites_by_acc):
        jobs = []
        for r in served_rows:
            a, b = r['uniprot_id1'], r['uniprot_id2']
            if '-' in a or '-' in b or a == b or a not in reference or b not in reference:   # human reviewed pairs only
                continue
            sc = float(r['pdockq'])
            if sc < GALLERY_PDOCKQ:
                continue
            want, new = {}, {}
            for acc, ch, col in ((a, r['chain1'] or 'A', 'ifresid1_ranges'), (b, r['chain2'] or 'B', 'ifresid2_ranges')):
                near = sites_by_acc.get(acc, set()) & ii.expand(r[col])
                if near:
                    want[ch] = (acc, sorted(near))
                    fresh = sorted(near - measured_sites.get(acc, set()) - measured_at.get(('ProtVar', f'{a}_{b}', acc), set()))
                    if fresh:
                        new[ch] = (acc, fresh)
            if want:
                jobs.append((a, b, sc, want, new))
        todo = [j for j in jobs if j[4]]
        log(f'ProtVar served models (pDockQ >= {GALLERY_PDOCKQ}) with a methylarginine at the interface: {len(jobs):,} pairs; '
            f'{len(todo):,} have positions the tables do not, measured into {topup}')
        fresh = {}
        with ThreadPoolExecutor(max_workers=workers) as pool:
            for (a, b, *_), res in zip(todo, pool.map(lambda j: measure_new('protvar_measure', f'{j[0]}_{j[1]}', lambda: get_bytes(PROTVAR_MODEL.format(a=j[0], b=j[1])), j[4]), todo)):
                fresh[(a, b)] = res
        out = []
        for a, b, sc, want, new in jobs:
            chains = {ch: dict(v) for ch, v in measured.get(('ProtVar', f'{a}_{b}'), {}).items()}
            for ch, sites in fresh.get((a, b), {}).items():
                chains.setdefault(ch, {}).update(sites)
            for ch, sites in chains.items():
                if ch not in want:
                    continue
                acc = want[ch][0]
                other = b if acc == a else a
                for pos, m in sites.items():
                    if int(pos) in want[ch][1] and partner_agrees(m, other):
                        out.append({'acc': acc, 'pos': int(pos), 'partner': other, 'src': 0, 'score': sc, 'pdockq': sc, 'm': m, 'ref': 1 if acc == a else 2, 'chain': ch})
        return out
    g['measure_protvar_models'] = measure_protvar_models_

    def measure_afdb_complexes_(cache, refresh, workers, sites_by_acc):
        accs = sorted(sites_by_acc)
        lists = topup / 'afdb_complex'
        lists.mkdir(parents=True, exist_ok=True)
        # an accession whose saved query listed more heterodimers than the tables hold models for (models staged
        # without their accessions before the uniprotAccession key was read) is queried again
        held = Counter(a for x in heterodimers.values() for a in x['accessions'])
        new_accs = [a for a in accs if a not in afdb_queried or held[a] < afdb_queried[a]]

        def get_list(acc):
            dest = lists / f'{acc}.json'
            if not dest.exists():
                data = get_bytes(AFDB_COMPLEX.format(acc=acc))
                dest.write_text(data.decode() if data else '[]')
        with ThreadPoolExecutor(max_workers=workers) as pool:
            list(pool.map(get_list, new_accs))
        models = {mid: dict(x, scores=(x['iptm'], x['pdockq'])) for mid, x in heterodimers.items() if x['accessions']}
        for acc in new_accs:
            for x in ii.read_json(lists / f'{acc}.json', []):
                if isinstance(x, dict) and x.get('assemblyType') == 'Hetero' and set(x.get('taxId') or []) == {9606} and x.get('modelEntityId'):
                    listed = x.get('uniprotAccession') or x.get('uniprotAccessions') or []
                    models.setdefault(x['modelEntityId'], {'accessions': [listed] if isinstance(listed, str) else list(listed), 'chains': [],
                                                          'scores': (x.get('complexPredictionAccuracy_ipTM'), x.get('complexPredictionAccuracy_pDockQ'))})
        models = {mid: x for mid, x in models.items() if set(x['accessions']) & set(accs)}
        preds = topup / 'afdb_prediction'
        preds.mkdir(parents=True, exist_ok=True)

        def meta_of(mid, x):
            if x['chains']:
                return list(zip(x['chains'], x['accessions']))
            pf = preds / f'{mid}.json'
            if not pf.exists():
                data = get_bytes(AFDB_PREDICTION.format(id=mid))
                pf.write_text(data.decode() if data else '[]')
            return [(c.get('chainId'), c.get('uniprotAccession')) for c in ii.read_json(pf, [])]

        def job(item):
            mid, x = item
            meta = meta_of(mid, x)
            new = {ch: (acc, sorted(sites_by_acc[acc] - (measured_sites.get(acc, set()) if mid in heterodimers else set())
                                    - measured_at.get(('AlphaFold DB', mid, acc), set())))
                   for ch, acc in meta if acc in sites_by_acc}
            new = {ch: v for ch, v in new.items() if v[1]}
            if not new:
                return mid, meta, {}

            def model_bytes():
                pf = preds / f'{mid}.json'
                if not pf.exists():
                    data = get_bytes(AFDB_PREDICTION.format(id=mid))
                    pf.write_text(data.decode() if data else '[]')
                info = ii.read_json(pf, [])
                return get_bytes(info[0]['pdbUrl']) if info and info[0].get('pdbUrl') else None
            return mid, meta, measure_new('afdb_measure', mid, model_bytes, new)
        log(f'AlphaFold DB: {len(models):,} human heterodimer models with a methylated protein '
            f'({len(new_accs):,} proteins not in the tables looked up into {topup})')
        with ThreadPoolExecutor(max_workers=workers) as pool:
            results = list(pool.map(job, sorted(models.items())))
        out = []
        for mid, meta, fresh in results:
            chains = {ch: dict(v) for ch, v in measured.get(('AlphaFold DB', mid), {}).items()}
            for ch, sites in fresh.items():
                chains.setdefault(ch, {}).update(sites)
            acc_of = dict(meta)
            iptm, pdockq = models[mid]['scores']
            for ch, sites in chains.items():
                acc = acc_of.get(ch)
                others = [a for c, a in meta if c != ch]
                if not acc or not others or acc not in sites_by_acc:
                    continue
                for pos, m in sites.items():
                    if int(pos) in sites_by_acc[acc] and partner_agrees(m, others[0]):
                        out.append({'acc': acc, 'pos': int(pos), 'partner': others[0], 'src': 1, 'score': iptm, 'pdockq': pdockq,
                                    'm': m, 'ref': mid, 'chain': ch})
        return out
    g['measure_afdb_complexes'] = measure_afdb_complexes_


# ----------------------------------------------------------------------------- release tables (--tables-out)
XL_SOURCE = ['Kramer2014', 'XRNAX', 'RBS-ID', 'pRBS-ID', 'pepR-MS', 'pepR-MS_focused']
# Categories of the release statistics for the crosslink and IntAct tests.
XL_TEST = {'rna_xl': 'rna_crosslink', 'rna_xl_single': 'rna_crosslink_single_residue', 'rna_near': 'rna_crosslink_within_3',
           'rna_rbdmap': 'rbdmap_fragment', 'dna_xl': 'dna_crosslink', 'dna_cand': 'dna_crosslink_or_candidate',
           'dna_pep': 'dna_crosslinked_peptide', 'dna_bind': 'dna_binding_region', 'rna_xl_nonrg': 'rna_crosslink_outside_rg',
           'rna_rbdmap_nonrg': 'rbdmap_fragment_outside_rg'}
INTACT_TEST = {'mutation': 'mutation_changes_interaction', 'peptide': 'binding_region', 'peptide_disordered': 'binding_region_disordered'}


def runs(positions):
    """Consecutive positions as [start, end] runs."""
    out = []
    for q in sorted(positions):
        if out and q == out[-1][1] + 1:
            out[-1][1] = q
        else:
            out.append([q, q])
    return out


def write_release_tables(outdir: Path, d: dict):
    """The interaction tables of a release, from what main() computed: site_interactions, residue_evidence,
    protein_partners, complexes, interaction_statistics and protein_interaction_summary (one row per protein: GO
    binding, partners, complexes). Written in the release's conventions; scripts/build_release.py adds them."""
    outdir.mkdir(parents=True, exist_ok=True)
    prots, sites, vocab, gene_of = d['prots'], d['sites'], d['vocab'], d['gene_of']
    sid = lambda j: f'{prots[sites[j][0]][0]}:R{sites[j][1]}'  # noqa: E731
    gene = lambda a: (gene_of.get(base_acc(a), '') or '').split(';')[0].strip()  # noqa: E731

    def write(name, header, rows):
        with open(outdir / name, 'w', encoding='utf-8', newline='') as fh:
            w = csv.writer(fh, delimiter='\t', lineterminator='\n')
            w.writerow(header)
            for r in rows:
                w.writerow(['true' if v is True else 'false' if v is False else '' if v is None else v for v in r])

    def partner(k, own):
        if k == -1:
            return 'rna', '', ''
        if k == -2:
            return 'dna', '', ''
        if k < 0:
            return 'other', '', ''
        a = base_acc(vocab[k])
        return ('self' if a == base_acc(own) else 'protein'), a, gene(a)

    # One row per piece of evidence at a site. counts_as_evidence: the row is one of the website's evidence classes
    # (the evidence map); interaction_status in the site table then counts only the crosslink, PDB and predicted rows.
    rows = []
    for j, key, pdbs in d['iface_detail']:
        own = prots[sites[j][0]][0]
        if key in ('RNA', 'DNA'):
            ev, kind, pa, pg = f'pdb_{key.lower()}', key.lower(), '', ''
        elif base_acc(key) == base_acc(own):
            ev, kind, pa, pg = 'pdb_self', 'self', '', ''
        else:
            ev, kind, pa, pg = 'pdb_protein', 'protein', base_acc(key), gene(key)
        rows.append([sid(j), ev, 'at_residue', kind, pa, pg, 'PDBe-KB', '; '.join(pdbs), '', '', f'{len(pdbs)} PDB entries', True])
    xs = d.get('xl_seq') or {'rna': {}, 'dna': {}, 'rbd': {}}

    def marked(pep, k):
        return pep[:k] + '[' + pep[k] + ']' + pep[k + 1:]

    def rna_peptides(acc, positions, bit):
        # Crosslinked peptides of one dataset at these residues: the crosslinked residue in brackets, with its RNA masses.
        return ', '.join(f'{marked(pep, k)} at {q} ({masses} Da)' for q in positions for b, pep, k, masses in xs['rna'].get(acc, {}).get(q, []) if b == bit)

    for j, rb, nb, db, fb, single in d['xl_rows']:
        acc, pos = prots[sites[j][0]][0], sites[j][1]
        for k, src in enumerate(XL_SOURCE):
            if rb >> k & 1:
                rows.append([sid(j), 'rna_crosslink', 'at_residue' if single else 'candidate_residue', 'rna', '', '', src, '', '', '',
                             rna_peptides(acc, [pos], k), True])
            elif not rb and nb >> k & 1:
                close = [q for q in range(pos - XL_NEAR, pos + XL_NEAR + 1) if q != pos]
                rows.append([sid(j), 'rna_crosslink', 'within_3', 'rna', '', '', src, '', '', '',
                             '. '.join(x for x in ('a crosslinked residue within 3 residues', rna_peptides(acc, close, k)) if x), False])
        if fb:
            frags = [f'{a}-{b} {frag}' for a, b, frag in xs['rbd'].get(acc, []) if a <= pos <= b]
            rows.append([sid(j), 'rbdmap_fragment', 'in_fragment', 'rna', '', '', 'RBDmap', '', '', '', ', '.join(frags), True])
        dna_peps = ', '.join(f'{a}-{a + len(pep) - 1} {pep}' for a, pep in xs['dna'].get(acc, []) if a <= pos <= a + len(pep) - 1)
        dna_note = (f'peptide {dna_peps} (lower case: the residues the crosslink may sit on)') if dna_peps else ''
        if db & D_XL:
            rows.append([sid(j), 'dna_crosslink', 'at_residue', 'dna', '', '', 'Trendel2025', '', '', '', dna_note, True])
        if db & D_CAND:
            rows.append([sid(j), 'dna_crosslink', 'candidate_residue', 'dna', '', '', 'Trendel2025', '', '', '', dna_note, True])
        if db & D_PEP:
            rows.append([sid(j), 'dna_crosslink', 'in_peptide', 'dna', '', '', 'Trendel2025', '', '', '',
                         '. '.join(x for x in ('inside a DNA-crosslinked peptide', dna_note) if x), False])
        if db & D_BIND:
            rows.append([sid(j), 'dna_binding_region', 'in_region', 'dna', '', '', 'UniProt', '', '', '', '', False])
    for j, hits in d['pred_all']:
        own = prots[sites[j][0]][0]
        for k, score, first, served in hits:
            kind, pa, pg = partner(k, own)
            rows.append([sid(j), 'predicted_interface', 'at_residue', kind, pa, pg, 'ProtVar', '', f'{score / 100:.2f}', '',
                         'model served by ProtVar' if served else '', True])
    for g in d['gallery']:
        j, k, src, score, pdockq, pl_r, pl_p, dmin, flags, labels, ref, chain = g
        own = prots[sites[j][0]][0]
        kind, pa, pg = partner(k, own)
        bonds = ', '.join(n for bit, n in ((G_SB, 'salt bridge'), (G_HB, 'hydrogen bond'), (G_CATPI, 'cation-pi'), (G_CLASH, 'clash')) if flags & bit)
        confident = not flags & G_CLASH and (src == 0 or (score or 0) >= 80) and (pl_r or 0) >= 80 and (pl_p or 0) >= 80 and bool(flags & (G_SB | G_HB))
        model = (f'{base_acc(own)}_{pa}' if ref == 1 else f'{pa}_{base_acc(own)}') if src == 0 else ref
        rows.append([sid(j), 'predicted_contact', 'at_residue', kind, pa, pg, 'ProtVar' if src == 0 else 'AFDB', model,
                     f'{(score or 0) / 100:.2f}', f'{dmin / 10:.1f}',
                     '. '.join(x for x in (bonds.capitalize(), f'pLDDT {pl_r} (arginine), {pl_p} (partner)', 'Partner residues ' + ', '.join(labels)) if x), confident])
    for j, recs in d['intact_rows']:
        own = prots[sites[j][0]][0]
        for kind_, code, a, b, res, pk, pmids, acs, n in recs:
            ev = ('intact_mutation', 'intact_binding_region', 'intact_methylated')[kind_]
            what = (INTACT_EFFECTS[code] if kind_ == 0 else INTACT_REGIONS[code] if kind_ == 1 else 'methylated arginine')
            detail = f'{what}, residues {a}-{b}' + (f' to {res}' if res else '') + (f', PubMed {", ".join(map(str, pmids))}' if pmids else '') + f', IntAct {", ".join(acs)}'
            counted = (1 <= code <= 4) if kind_ == 0 else True
            for k in (pk or [None]):
                kind, pa, pg = partner(k, own) if k is not None else ('', '', '')
                rows.append([sid(j), ev, 'in_region' if a != b else 'at_residue', kind, pa, pg, 'IntAct', '', '', '', detail, counted])
    order = {e: k for k, e in enumerate(SITE_EVIDENCE)}
    site_order = {sid(j): j for j in range(len(sites))}
    rows.sort(key=lambda r: (site_order[r[0]], order[r[1]], r[2], r[6], r[4], r[7]))
    write('site_interactions.tsv', ['site_id', 'evidence_type', 'localisation', 'partner_kind', 'partner_accession', 'partner_gene',
                                    'source', 'structures', 'score', 'distance_a', 'detail', 'counts_as_evidence'], rows)

    # Evidence along each protein, as runs of consecutive residues
    rows = []
    for acc, res in d['iface_residues'].items():
        by_type = defaultdict(dict)
        for pos, keys in res.items():
            for key in keys:
                t = 'pdb_rna' if key == 'RNA' else 'pdb_dna' if key == 'DNA' else 'pdb_self' if base_acc(key) == base_acc(acc) else 'pdb_protein'
                by_type[t].setdefault(pos, set()).add(key)
        for t, at in by_type.items():
            for a, b in runs(at):
                rows.append([acc, a, b, t, 'at_residue', max(len(at[q]) for q in range(a, b + 1)) if t == 'pdb_protein' else '', '', 'PDBe-KB'])
    for acc, rs in d['pred_iface'].items():
        for a, b, n, q in rs:
            rows.append([acc, a, b, 'predicted_interface', 'at_residue', n, f'{q / 100:.2f}', 'ProtVar'])
    for i, r_pos, r_bits, d_xl, d_peps, d_bind, _in_int, frags, d_cand in d['xl_prot']:
        acc = prots[i][0]
        for q, bits in zip(r_pos, r_bits):
            rows.append([acc, q, q, 'rna_crosslink', 'at_residue', '', '', '; '.join(src for k, src in enumerate(XL_SOURCE) if bits >> k & 1)])
        for q in d_xl:
            rows.append([acc, q, q, 'dna_crosslink', 'at_residue', '', '', 'Trendel2025'])
        for q in d_cand:
            rows.append([acc, q, q, 'dna_crosslink', 'candidate_residue', '', '', 'Trendel2025'])
        for a, b in d_peps:
            rows.append([acc, a, b, 'dna_crosslink', 'in_peptide', '', '', 'Trendel2025'])
        for a, b in d_bind:
            rows.append([acc, a, b, 'dna_binding_region', 'in_region', '', '', 'UniProt'])
        for a, b in frags:
            rows.append([acc, a, b, 'rbdmap_fragment', 'in_fragment', '', '', 'RBDmap'])
    for i, regions, muts in d['intact_prot']:
        acc = prots[i][0]
        for a, b, bits, _ps, n, _rec in regions:
            rows.append([acc, a, b, 'intact_binding_region', 'in_region', n, '', 'IntAct'])
        for a, b, code, _ps, n, _rec in muts:
            rows.append([acc, a, b, 'intact_mutation', 'in_region' if a != b else 'at_residue', n, '', 'IntAct'])
    rows.sort(key=lambda r: (r[0], r[1], r[2], order[r[3]], r[4]))
    write('residue_evidence.tsv', ['uniprot_accession', 'start', 'end', 'evidence_type', 'localisation', 'n_partners', 'best_score', 'sources'], rows)

    # STRING partners of each methylated protein
    rows = []
    for i, lst in sorted(d['partners'].items()):
        own = prots[i][0]
        for k, sc in sorted(lst, key=lambda x: (-x[1], vocab[x[0]])):
            a = base_acc(vocab[k])
            classes = [c for c, ok in (('prmt', gene(a) in PRMT_GENES), ('tudor', a in d['domain_classes']['tudor']),
                                       ('wd40', a in d['domain_classes']['wd40']), ('phd', a in d['domain_classes']['phd'])) if ok]
            rows.append([own, a, gene(a), f'{sc / 1000:.3f}', a in d['meth_base'], '; '.join(classes)])
    write('protein_partners.tsv', ['uniprot_accession', 'partner_accession', 'partner_gene', 'string_score', 'partner_methylated', 'partner_classes'], rows)

    # Complexes with a methylated member (Complex Portal curated, hu.MAP predicted)
    rows, n_cx = [], Counter()
    for cid, name, predicted, slots, n, k, fold, p, q in d['complexes']:
        methylated = sorted({a for s in slots for a in s if a in d['meth_base']})
        n_cx.update(methylated)
        rows.append([cid, '' if name in ('', '-') else name, 'predicted' if predicted else 'curated', n, k, '; '.join(methylated)]
                    + ['' if x is None else f'{x:.4g}' for x in (fold, p, q)])
    write('complexes.tsv', ['complex_id', 'complex_name', 'curation', 'n_members', 'n_members_methylated', 'methylated_members',
                            'fold_enrichment', 'p', 'q'], rows)

    # Statistics, in the columns of the release's statistics table. Residue tests: pooled Fisher (Haldane-Anscombe odds
    # ratio) and, for the within-protein tests, Mantel-Haenszel with each protein as a stratum; q within each analysis
    # and method. Protein-set tests: clusterProfiler enricher, fold enrichment and its BH p.adjust.
    st = []
    for an, cat, stratum, tabs, background in d['tests']:
        a, b, c, dd = (sum(t[k] for t in tabs) for k in range(4))
        if not a + b + c + dd:
            continue
        p, odds = fisher_2x2(a, b, c, dd)
        st.append([an, cat, stratum, a, b, c, dd, odds, math.log2(odds), None, None, p, None, 'fisher', background, None])
        if len(tabs) > 1:
            mh_or, mh_p, used = mantel_haenszel(tabs)
            ok = mh_or is not None and math.isfinite(mh_or) and mh_or > 0
            st.append([an, cat, stratum, a, b, c, dd, mh_or if ok else None, math.log2(mh_or) if ok else None, None, None, mh_p,
                       None, 'mantel_haenszel', background + ', each protein a stratum', used])
    for key in dict.fromkeys((r[0], r[13]) for r in st):
        idx = [k for k, r in enumerate(st) if (r[0], r[13]) == key and r[11] is not None]
        for k, q in zip(idx, bh([st[k][11] for k in idx])):
            st[k][12] = q
    for an, cat, stratum, a, b, c, dd, fold, p, q, background in d['set_tests']:
        st.append([an, cat, stratum, a, b, c, dd, None, None, fold, math.log2(fold), p, q, 'enricher', background, None])
    fmt = lambda x: '' if x is None else f'{x:.6g}' if isinstance(x, float) else x  # noqa: E731
    write('interaction_statistics.tsv', ['analysis', 'category', 'stratum', 'a', 'b', 'c', 'd', 'odds_ratio', 'log2_odds_ratio',
                                         'fold_enrichment', 'log2_fold_enrichment', 'p', 'q', 'method', 'background', 'n_strata'],
               [[fmt(x) for x in r] for r in st])

    # One row per methylated protein
    rows = []
    for i, p in enumerate(prots):
        a = base_acc(p[0])
        go_ = d['go'].get(a)
        lvl = lambda key: ('experimental' if go_[key] & GO_EXPERIMENTAL else 'other_evidence' if go_[key] else 'none') if go_ else 'none'  # noqa: E731
        lst = d['partners'].get(i, [])
        rows.append([p[0], lvl('rna'), lvl('dna'), len(lst), sum(1 for k, _ in lst if base_acc(vocab[k]) in d['meth_base']), n_cx.get(a, 0)])
    write('protein_interaction_summary.tsv', ['uniprot_accession', 'go_rna_binding', 'go_dna_binding', 'n_string_partners',
                                              'n_partners_methylated', 'n_complexes'], rows)
    log(f'wrote the release interaction tables to {outdir}')


SITE_EVIDENCE = ['rna_crosslink', 'rbdmap_fragment', 'dna_crosslink', 'dna_binding_region', 'pdb_protein', 'pdb_rna', 'pdb_dna',
                 'pdb_self', 'predicted_interface', 'predicted_contact', 'intact_mutation', 'intact_binding_region', 'intact_methylated']


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--sites-json', default='docs/data/methylarg_sites.json')
    ap.add_argument('--out', default='docs/data/methylarg_interactions.json')
    ap.add_argument('--pred-iface-out', default='docs/data/methylarg_pred_iface.json',
                    help='every predicted interface residue of the explorer proteins, loaded by the explorer on demand')
    ap.add_argument('--cache', default='PTM_data/interactions_cache')
    ap.add_argument('--string-min-score', type=int, default=700, help='STRING physical score floor (0-1000); 700 = high confidence. Protein pages fetch lower-confidence links live.')
    ap.add_argument('--workers', type=int, default=6, help='parallel PDBe requests')
    ap.add_argument('--refresh', action='store_true', help='download every source again')
    ap.add_argument('--reference-fasta', default='',
                    help='restrict the proteome to the entries of this UniProt FASTA (the release reference, e.g. '
                         'inputs/reference/uniprot_human_reviewed_2026_03.fasta); stops if a downloaded sequence differs from it')
    ap.add_argument('--fetch-interfaces-only', action='store_true', help='only fill the PDBe interface cache')
    ap.add_argument('--fetch-contacts-only', action='store_true', help='only fill the structure-contact cache (needs the interface cache)')
    ap.add_argument('--contact-entries', type=int, default=2, help='PDB entries examined per site and partner type')
    ap.add_argument('--fetch-environments-only', action='store_true', help='only fill the 10 Å residue-environment cache (interface geometry)')
    ap.add_argument('--controls-per-protein', type=int, default=12, help='non-methylated interface arginines sampled per protein for the geometry comparison')
    ap.add_argument('--fetch-models-only', action='store_true', help='only measure methylarginines in predicted complex models (ProtVar, AlphaFold DB)')
    ap.add_argument('--model-workers', type=int, default=8, help='parallel model downloads for the predicted-contact gallery')
    ap.add_argument('--tables', default='', help='read the frozen tables in this folder (inputs/interactions) instead of --cache')
    ap.add_argument('--topup-cache', default='PTM_data/interactions_topup',
                    help='with --tables: where to save what the tables do not cover (scripts/37a_stage_interactions.py --topup stages it)')
    ap.add_argument('--measured-sites', default='releases/v1.0.0/human_methylarginine_sites_v1.0.0.tsv',
                    help='with --tables: the sites the frozen model measurements were made at (the published union)')
    ap.add_argument('--tables-out', default='', help='also write the interaction tables of a release into this folder')
    ap.add_argument('--enrichment-workdir', default='',
                    help='where the clusterProfiler runs keep their inputs and results (default: <tables-out>/enricher, '
                         'or PTM_data/interactions_enricher)')
    args = ap.parse_args()
    cache = Path(args.cache)
    enrich_dir = Path(args.enrichment_workdir or (Path(args.tables_out) / 'enricher' if args.tables_out else 'PTM_data/interactions_enricher'))
    if args.tables:
        measured_sites = defaultdict(set)
        with open(args.measured_sites, encoding='utf-8') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                measured_sites[r['uniprot_accession']].add(int(r['position']))
        use_frozen_tables(Path(args.tables), Path(args.topup_cache), measured_sites)

    explorer = json.loads(Path(args.sites_json).read_text())
    prots, sites = explorer['proteins'], explorer['sites']
    nP = len(prots)
    meth_accs = [p[0] for p in prots]
    canonical = [p[0] for p in prots if '-' not in p[0]]

    if args.fetch_interfaces_only:
        fetch_interfaces(cache, canonical, args.refresh, args.workers)
        return

    proteome, domain_classes = load_uniprot(cache, args.refresh)
    reference = None
    if args.reference_fasta:
        ref = read_fasta(Path(args.reference_fasta))
        differ = sorted(a for a, s in ref.items() if a in proteome and proteome[a][2] != s)
        missing = sorted(a for a in ref if a not in proteome)
        if differ or missing:
            raise SystemExit(f'{args.reference_fasta} is not the UniProt release of the download: {len(differ)} sequences differ '
                             f'({", ".join(differ[:5])}), {len(missing)} entries missing ({", ".join(missing[:5])})')
        log(f'Reference FASTA: {len(ref):,} entries; {len(proteome) - len(ref):,} downloaded entries outside it left out')
        proteome = {a: v for a, v in proteome.items() if a in ref}
        domain_classes = {k: {a for a in v if a in ref} for k, v in domain_classes.items()}
        reference = {'fasta': Path(args.reference_fasta).name, 'entries': len(ref),
                     'md5': hashlib.md5(Path(args.reference_fasta).read_bytes()).hexdigest()}
    tudor = domain_classes['tudor']
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

    # Complexes: membership, methylated share and enrichment against the reviewed proteome (clusterProfiler enricher,
    # curated and predicted complexes as separate jobs, complexes of at least two reviewed proteins). A complex's set is
    # its reviewed member proteins; a slot that any of several proteins can fill contributes each of them.
    N = len(proteome)
    K = sum(1 for a in proteome if a in meth_base)
    jobs = [Job(f'complexes_{kind}', [a for a in proteome if a in meth_base], proteome, whole_universe=True, min_size=2,
                sets=pd.DataFrame([(c['id'], a) for c in complexes if c['predicted'] == (kind == 'predicted') for s in c['slots']
                                   for a in s if a in proteome], columns=['term', 'gene']))
            for kind in ('curated', 'predicted')]
    res = run_enrichment(jobs, enrich_dir)
    cx_test = {x.ID: x for j in jobs for _, x in res[j.job].iterrows()}
    n_curated, n_predicted = (len(res[j.job]) for j in jobs)
    cx_rows, cx_release = [], []
    for c in complexes:
        slots = [[a for a in s if a in proteome] or s for s in c['slots']]
        reviewed = {a for s in slots for a in s if a in proteome}
        methylated = {a for a in reviewed if a in meth_base}
        if not methylated:
            continue
        members = [vid(s[0]) if len(s) == 1 else [vid(a) for a in s] for s in slots]
        t = cx_test.get(c['id'])
        fold, p, q = (float(t.FoldEnrichment), float(t.pvalue), float(t['p.adjust'])) if t is not None else (None, None, None)
        cx_rows.append([c['id'], c['name'], int(c['predicted']), members, len(reviewed), len(methylated), sig(fold), sig(p), sig(q)])
        cx_release.append([c['id'], c['name'], bool(c['predicted']), slots, len(reviewed), len(methylated), fold, p, q])

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
    iface_detail = []      # (site j, partner key, PDB entries) for the release tables
    release_tests = []     # (analysis, category, stratum, per-protein 2x2 tables, background) for the release statistics
    iface_residues = {}    # accession -> {position: {partner key: PDB entries}}
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
        iface_residues[acc] = res_partners
        for pos, j in site_pos.items():
            hit = res_partners.get(pos)
            if not hit:
                continue
            for key, ids in hit.items():
                iface_detail.append((j, key, sorted(ids)))
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
            release_tests.append(('pdb_interface', key, f, tabs, 'other arginines of the same proteins observed in a PDB structure'))
            a, b, c, d = (sum(t[k] for t in tabs) for k in range(4))
            p, odds = fisher_2x2(a, b, c, d) if a + b + c + d else (None, None)
            mh_or, mh_p, used = mantel_haenszel(tabs)
            enrich_strata[key][f] = {
                'pooled': [a, b, c, d, sig(math.log2(odds)) if odds else None, sig(p)],
                'mh': [sig(math.log2(mh_or)) if mh_or and math.isfinite(mh_or) else None, sig(mh_p), used],
            }

    # How many methylated arginines in each RG context are structurally observed at all (all canonical proteins).
    for i, acc in enumerate(meth_accs):
        if '-' in acc:
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

    # Predicted interfaces (ProtVar, AlphaFold-Multimer): flags for methylated sites and a within-protein test
    # over every arginine of the sequence (predicted models are full length, so disordered regions count too).
    pred = load_protvar(cache, args.refresh, set(meth_base))
    served = load_protvar_served(cache, args.refresh)

    def model_order(acc, other, first):
        """1 if the API serves the pair in the 2024 table's order, 2 if only reversed, 0 if not at all."""
        a, b = (acc, other) if first else (other, acc)
        return 1 if (a, b) in served else 2 if (b, a) in served else 0

    pred_rows, pred_all = [], []
    pred_strata = {(t, f): defaultdict(lambda: [0, 0, 0, 0]) for t in PRED_TIERS for f in ('all', 'ordered', 'disordered')}
    for i, acc in enumerate(meth_accs):
        if '-' in acc or acc not in pred:
            continue
        pairs = pred[acc]
        for j in by_prot_sites[i]:
            pos = sites[j][1]
            hits = sorted(((vid(o), round(sc * 100), int(first), model_order(acc, o, first)) for o, sc, pos_r, _, first in pairs if pos in pos_r), key=lambda h: -h[1])
            if hits:
                pred_rows.append([j, [list(h) for h in hits[:12]], len(hits)])
                pred_all.append((j, hits))
        seq = proteome.get(acc, ('', '', ''))[2]
        if not seq:
            continue
        dis = disordered.get(acc, [])
        site_pos = {sites[j][1] for j in by_prot_sites[i]}
        for t in PRED_TIERS:
            iface_r = set()
            for o, sc, pos_r, _, first in pairs:
                if sc >= t:
                    iface_r |= pos_r
            if not iface_r:
                continue
            for pos in range(1, len(seq) + 1):
                if seq[pos - 1] != 'R':
                    continue
                cell = (0 if pos in site_pos else 2) + (0 if pos in iface_r else 1)
                pred_strata[(t, 'all')][acc][cell] += 1
                pred_strata[(t, 'disordered' if any(a <= pos <= b for a, b in dis) else 'ordered')][acc][cell] += 1
    pred_enrich = {}
    for t in PRED_TIERS:
        pred_enrich[str(t)] = {}
        for f in ('all', 'ordered', 'disordered'):
            tabs = list(pred_strata[(t, f)].values())
            release_tests.append(('predicted_interface', f'pdockq_{t}', f, tabs, 'other arginines of the same proteins'))
            a, b, c, d = (sum(x[k] for x in tabs) for k in range(4))
            p, odds = fisher_2x2(a, b, c, d) if a + b + c + d else (None, None)
            mh_or, mh_p, used = mantel_haenszel(tabs)
            pred_enrich[str(t)][f] = {'pooled': [a, b, c, d, sig(math.log2(odds)) if odds else None, sig(p)],
                                       'mh': [sig(math.log2(mh_or)) if mh_or and math.isfinite(mh_or) else None, sig(mh_p), used]}
    log(f'Predicted interfaces: {len(pred_rows):,} methylated sites at a predicted heterodimer interface')
    pred_iface = predicted_interface_runs({acc: pred[acc] for acc in meth_accs if '-' not in acc and acc in pred})

    # Predicted-contact gallery: every methylarginine within 5 Å of the partner in a served ProtVar model
    # (pDockQ >= 0.5) or an AlphaFold DB human heterodimer, with salt bridges, H-bonds, cation-pi and pLDDT.
    sites_by_acc = defaultdict(set)
    for s in sites:
        acc = meth_accs[s[0]]
        if '-' not in acc and not prots[s[0]][10]:
            sites_by_acc[acc].add(s[1])
    measured = measure_protvar_models(cache, args.refresh, args.model_workers, sites_by_acc)
    measured += measure_afdb_complexes(cache, args.refresh, max(2, args.model_workers // 2), sites_by_acc)
    if args.fetch_models_only:
        return
    gallery = []
    for r in measured:
        j = site_index.get((r['acc'], r['pos']))
        if j is None:
            continue
        dmin, flags, pl_r, pl_p, labels = r['m']
        pct = lambda v: round(100 * v) if v is not None else None  # noqa: E731
        gallery.append([j, vid(r['partner']), r['src'], pct(r['score']), pct(r['pdockq']), pl_r, pl_p, round(dmin * 10), flags, [lab for _, lab in labels], r['ref'], r['chain']])
    gallery.sort(key=lambda g: (g[0], -(g[3] or 0)))
    log(f'Predicted-contact gallery: {len(gallery):,} methylarginine-partner contacts (<= {CONTACT_KEEP} Å) on {len({g[0] for g in gallery}):,} sites; '
        f'{sum(1 for g in gallery if g[8] & (G_SB | G_HB)):,} with a salt bridge or H-bond')

    # IntAct: curated mutations, short binding regions and methylated residues at the sites; per-protein tracks;
    # within-protein tests of methylated against other arginines.
    intact, intact_counts, intact_release = load_intact(cache, args.refresh, {a: v[2] for a, v in proteome.items()})

    def intact_partner(p):
        db, ident, ptype = p
        t = ptype.lower()
        if 'deoxyribonucleic' in t or t.endswith('dna'):
            return -2
        if 'ribonucleic' in t or t.endswith('rna'):
            return -1
        if db == 'uniprotkb' and t in ('protein', 'peptide') and ACC_RE.match(ident.split('-')[0]):
            return vid(ident)
        return -3 if 'small molecule' in t else -4

    intact_rows, intact_prot = [], []
    intact_strata = {k: defaultdict(lambda: [0, 0, 0, 0]) for k in ('mutation', 'peptide', 'peptide_disordered')}
    for i, acc in enumerate(meth_accs):
        if '-' in acc or acc not in intact:
            continue
        site_pos = {sites[j][1]: j for j in by_prot_sites[i]}
        at_site, regions, muts = defaultdict(dict), {}, {}
        r_effect = defaultdict(set)       # arginine -> effects of single-residue mutations there
        peptide_spans = []
        for kind, code, a, b, res, participants, pmid, ac in intact[acc]:
            pk = tuple(sorted({intact_partner(x) for x in participants}))
            span = b - a + 1
            if kind == 1:
                r = regions.setdefault((a, b), [0, set(), 0])
                r[0] |= 1 << code
                r[1].update(pk)
                r[2] += 1
                if span <= INTACT_PEPTIDE_MAX:
                    peptide_spans.append((a, b))
            elif kind == 0 and span <= INTACT_MUT_MAX:
                m = muts.setdefault((a, b, code), [set(), 0])
                m[0].update(pk)
                m[1] += 1
                if span == 1 and code != 5:
                    r_effect[a].add(code)
            if (kind == 1 and span > INTACT_PEPTIDE_MAX) or (kind == 0 and span > INTACT_MUT_MAX):
                continue
            for q in range(a, b + 1):
                if q in site_pos:
                    rec = at_site[site_pos[q]].setdefault((kind, code, a, b, res, pk), [set(), []])
                    if pmid:
                        rec[0].add(pmid)
                    rec[1].append(ac)
        for j, recs in sorted(at_site.items()):
            intact_rows.append([j, [[kind, code, a, b, res, list(pk), sorted(pm), sorted(set(acs))[:5], len(acs)]
                                    for (kind, code, a, b, res, pk), (pm, acs) in sorted(recs.items())]])
        intact_prot.append([i, [[a, b, bits, sorted(ps)[:8], len(ps), n] for (a, b), (bits, ps, n) in sorted(regions.items())],
                            [[a, b, code, sorted(ps)[:8], len(ps), n] for (a, b, code), (ps, n) in sorted(muts.items())]])
        seq_ = proteome.get(acc, ('', '', ''))[2]
        dis = disordered.get(acc, [])
        for q in range(1, len(seq_) + 1):
            if seq_[q - 1] != 'R':
                continue
            row = 0 if q in site_pos else 2
            effects = r_effect.get(q, set())
            if effects:
                intact_strata['mutation'][acc][row + (0 if effects - {0} else 1)] += 1
            if peptide_spans:
                inside = 0 if any(a <= q <= b for a, b in peptide_spans) else 1
                intact_strata['peptide'][acc][row + inside] += 1
                if any(a <= q <= b for a, b in dis):
                    intact_strata['peptide_disordered'][acc][row + inside] += 1
    intact_enrich = {k: or_summary(v.values()) for k, v in intact_strata.items()}
    for k, v in intact_strata.items():
        release_tests.append(('intact', INTACT_TEST[k], 'all', list(v.values()),
                              'other arginines of the same proteins' + (' with an IntAct mutation' if k == 'mutation' else ' with an IntAct binding region')))
    log(f'IntAct: {len(intact_rows):,} methylated sites with residue-level features on {len(intact_prot):,} proteins; '
        + ', '.join(f'{k} MH OR {v["mh"][0]} (p {v["mh"][3]})' for k, v in intact_enrich.items()))

    # Residue-level RNA and DNA crosslinks: site flags, per-protein tracks and within-protein tests.
    rna_xl, rna_count, rbd_frag, dna_xl, dna_peps, dna_int, dna_bind, rna_single, xl_check = load_crosslinks(cache, args.refresh, {a: v[2] for a, v in proteome.items()}, proteome)
    xl_seq_all = load_crosslink_sequences({a: v[2] for a, v in proteome.items()}, proteome)
    xl_rows, xl_prot, xl_seq = [], [], []
    xl_strata = {k: defaultdict(lambda: [0, 0, 0, 0]) for k in ('rna_xl', 'rna_xl_single', 'rna_near', 'rna_rbdmap', 'dna_xl', 'dna_cand', 'dna_pep', 'dna_bind', 'rna_xl_nonrg', 'rna_rbdmap_nonrg')}
    for i, acc in enumerate(meth_accs):
        if '-' in acc:
            continue
        r_map, d_map, peps, binds, frags = rna_xl.get(acc, {}), dna_xl.get(acc, {}), dna_peps.get(acc, []), dna_bind.get(acc, []), rbd_frag.get(acc, [])
        r_single = rna_single.get(acc, set())
        site_pos = {sites[j][1] for j in by_prot_sites[i]}
        def near_bits(pos):
            b = 0
            for q in range(pos - XL_NEAR, pos + XL_NEAR + 1):
                b |= r_map.get(q, 0)
            return b
        def dna_bits(pos):
            return (d_map.get(pos, 0)) | (D_PEP if any(a <= pos <= b for a, b in peps) else 0) | (D_BIND if any(a <= pos <= b for a, b in binds) else 0)
        for j in by_prot_sites[i]:
            pos = sites[j][1]
            rb, nb, db = r_map.get(pos, 0), near_bits(pos), dna_bits(pos)
            fb = 1 if any(a <= pos <= b for a, b in frags) else 0
            if rb or nb or db or fb:
                xl_rows.append([j, rb, nb, db, fb, 1 if pos in r_single else 0])
        # The peptides and fragments behind the evidence near each site (what the evidence pop-ups show).
        near = lambda a, b: any(a - XL_NEAR <= q <= b + XL_NEAR for q in site_pos)  # noqa: E731
        seq_rna = [[q, bit, pep, k, masses] for q, hits in sorted(xl_seq_all['rna'].get(acc, {}).items()) if near(q, q) for bit, pep, k, masses in hits]
        seq_dna = [[a, pep] for a, pep in xl_seq_all['dna'].get(acc, []) if near(a, a + len(pep) - 1)]
        seq_rbd = [[a, b, frag] for a, b, frag in xl_seq_all['rbd'].get(acc, []) if near(a, b)]
        if seq_rna or seq_dna or seq_rbd:
            xl_seq.append([i, seq_rna, seq_dna, seq_rbd])
        if r_map or d_map or peps or binds or frags or acc in dna_int:
            xl_prot.append([i, sorted(r_map), [r_map[q] for q in sorted(r_map)], sorted(q for q, v in d_map.items() if v & D_XL), sorted(set(peps)), sorted(set(binds)),
                            1 if acc in dna_int else 0, sorted(set(frags)), sorted(q for q, v in d_map.items() if v == D_CAND)])
        seq_ = proteome.get(acc, ('', '', ''))[2]
        if not seq_:
            continue
        for key, test, have in (('rna_xl', lambda q: bool(r_map.get(q)), bool(r_map)), ('rna_xl_single', lambda q: q in r_single, bool(r_single)), ('rna_near', lambda q: bool(near_bits(q)), bool(r_map)),
                                ('rna_rbdmap', lambda q: any(a <= q <= b for a, b in frags), bool(frags)),
                                ('dna_xl', lambda q: bool(d_map.get(q, 0) & D_XL), any(v & D_XL for v in d_map.values())),
                                ('dna_cand', lambda q: bool(d_map.get(q, 0)), bool(d_map)), ('dna_pep', lambda q: any(a <= q <= b for a, b in peps), bool(peps)),
                                ('dna_bind', lambda q: any(a <= q <= b for a, b in binds), bool(binds))):
            if not have:
                continue
            for q in range(1, len(seq_) + 1):
                if seq_[q - 1] == 'R':
                    cell = (0 if q in site_pos else 2) + (0 if test(q) else 1)
                    xl_strata[key][acc][cell] += 1
                    # RG/RGG regions both bind RNA and carry most methylarginines: repeat the RNA tests without
                    # arginines in an RG or GR dipeptide.
                    if key in ('rna_xl', 'rna_rbdmap') and not ((q < len(seq_) and seq_[q] == 'G') or (q > 1 and seq_[q - 2] == 'G')):
                        xl_strata[key + '_nonrg'][acc][cell] += 1
    xl_enrich = {}
    for key, tabs in xl_strata.items():
        tabs = list(tabs.values())
        a, b, c, d = (sum(x[k] for x in tabs) for k in range(4))
        pv, odds = fisher_2x2(a, b, c, d) if a + b + c + d else (None, None)
        mh_or, mh_p, used = mantel_haenszel(tabs)
        xl_enrich[key] = {'pooled': [a, b, c, d, sig(math.log2(odds)) if odds else None, sig(pv)],
                          'mh': [sig(math.log2(mh_or)) if mh_or and math.isfinite(mh_or) else None, sig(mh_p), used], 'proteins': len(tabs)}
    for key, tabs in xl_strata.items():
        xl_enrich[key]['or'] = or_summary(tabs.values())
        release_tests.append(('crosslink', XL_TEST[key], 'all', list(tabs.values()), 'other arginines of the same proteins (proteins with that evidence)'))

    # By methylation type: each class of methylated arginine vs the unmethylated arginines (not in the union) of the
    # same proteins. Maron IP classes come from methyl-specific immunoprecipitation; state calls from the sources.
    type_classes = [('rme2s', lambda r: r[17] & 4), ('rme2a', lambda r: (r[17] & 2) and not (r[17] & 4)), ('rme1ip', lambda r: r[17] == 1),
                    ('rme2', lambda r: r[2] in (2, 3)), ('rme1', lambda r: r[2] == 1), ('nostate', lambda r: r[2] == 0 and not r[17]), ('any', lambda r: True)]
    xl_types = {}
    xl_freq = {}
    for out_key, have, test in (('rna_xl', lambda acc: bool(rna_xl.get(acc)), lambda acc, q: bool(rna_xl.get(acc, {}).get(q))),
                                ('rna_xl_single', lambda acc: bool(rna_single.get(acc)), lambda acc, q: q in rna_single.get(acc, set())),
                                ('rna_rbdmap', lambda acc: bool(rbd_frag.get(acc)), lambda acc, q: any(a <= q <= b for a, b in rbd_frag.get(acc, []))),
                                ('dna_xl', lambda acc: any(v & D_XL for v in dna_xl.get(acc, {}).values()), lambda acc, q: bool(dna_xl.get(acc, {}).get(q, 0) & D_XL))):
        tabs = {name: defaultdict(lambda: [0, 0, 0, 0]) for name, _ in type_classes}
        for i, acc in enumerate(meth_accs):
            if '-' in acc or not have(acc):
                continue
            seq_ = proteome.get(acc, ('', '', ''))[2]
            if not seq_:
                continue
            row_at = {sites[j][1]: sites[j] for j in by_prot_sites[i]}
            for q in range(1, len(seq_) + 1):
                if seq_[q - 1] != 'R':
                    continue
                hit = test(acc, q)
                r = row_at.get(q)
                for name, member in type_classes:
                    if r is None:
                        tabs[name][acc][2 if hit else 3] += 1
                    elif member(r):
                        tabs[name][acc][0 if hit else 1] += 1
        xl_types[out_key] = {name: or_summary(tabs[name].values()) for name, _ in type_classes}

    # Contact frequency: among RNA-crosslinked arginines, are methylated ones crosslinked in more datasets or
    # more spectra than unmethylated ones of the same proteins?
    freq_meth, freq_un, ds_meth, ds_un = [], [], [], []
    multi = defaultdict(lambda: [0, 0, 0, 0])
    for i, acc in enumerate(meth_accs):
        if '-' in acc or not rna_xl.get(acc):
            continue
        site_pos = {sites[j][1] for j in by_prot_sites[i]}
        seq_ = proteome.get(acc, ('', '', ''))[2]
        for q, bits in rna_xl[acc].items():
            if not seq_ or q > len(seq_) or seq_[q - 1] != 'R':
                continue
            nds = bin(bits & 31).count('1')
            if not nds:
                continue    # only in Table S3's focused sites: no dataset or spectrum count to compare
            nsp = rna_count[acc].get(q, 0)
            if q in site_pos:
                freq_meth.append(nsp); ds_meth.append(nds)
            else:
                freq_un.append(nsp); ds_un.append(nds)
            multi[acc][(0 if q in site_pos else 2) + (0 if nds >= 2 else 1)] += 1
    med = lambda v: sorted(v)[len(v) // 2] if v else None  # noqa: E731
    xl_freq = {
        'spectra': {'meth': [len(freq_meth), sig(sum(freq_meth) / len(freq_meth)) if freq_meth else None, med(freq_meth)],
                    'other': [len(freq_un), sig(sum(freq_un) / len(freq_un)) if freq_un else None, med(freq_un)], 'p': sig(mann_whitney(freq_meth, freq_un)[1])},
        'datasets': {'meth': [len(ds_meth), sig(sum(ds_meth) / len(ds_meth)) if ds_meth else None, med(ds_meth)],
                     'other': [len(ds_un), sig(sum(ds_un) / len(ds_un)) if ds_un else None, med(ds_un)], 'p': sig(mann_whitney(ds_meth, ds_un)[1])},
        'multi': or_summary(multi.values()),
    }
    log(f'Crosslinks by type: ' + ', '.join(f"{k} MH OR {v['mh'][0]} (p {v['mh'][3]})" for k, v in xl_types['rna_xl'].items()))
    log(f"Crosslink frequency: methylated R mean spectra {xl_freq['spectra']['meth'][1]} vs {xl_freq['spectra']['other'][1]} (p {xl_freq['spectra']['p']}); "
        f"datasets {xl_freq['datasets']['meth'][1]} vs {xl_freq['datasets']['other'][1]} (p {xl_freq['datasets']['p']})")
    log(f'Crosslinks: {sum(1 for r in xl_rows if r[1]):,} methylarginines are RNA-crosslinked residues, {sum(1 for r in xl_rows if r[2]):,} within {XL_NEAR} residues; '
        f'{sum(1 for r in xl_rows if r[3] & D_XL):,} DNA-crosslinked (single localisation), {sum(1 for r in xl_rows if r[3] & D_CAND):,} DNA candidate residues, '
        f'{sum(1 for r in xl_rows if r[3] & D_BIND):,} in UniProt DNA-binding regions')

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
        chain_map = defaultdict(list)   # SIFTS segments of a chain; a fusion chain has segments of several proteins
        for pdb, rows in chains.items():
            for ch, sp, beg, sp_beg, sp_end in rows:
                chain_map[(pdb, ch)].append((beg, sp_beg, sp_end, base_acc(sp)))
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
                            for res_beg, sp_beg, sp_end, seg_acc in chain_map.get((pdb, hch), []):
                                u = sp_beg + (int(hlab) - res_beg)
                                if sp_beg <= u <= sp_end:
                                    unp, pacc = u, seg_acc   # the protein of the segment that holds the residue
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
    cond_conf = cdcode_threshold_enrichment(cd, proteome, meth_base, cd_uids, enrich_dir)
    # How CD-CODE confidence scores are distributed over human memberships (scores cluster at 1-2 and 4-5).
    score_hist = Counter(m[0] for uid in cd_uids for a, m in cd.get(uid, {}).items() if a in proteome)
    cond_score_hist = {str(k): v for k, v in sorted(score_hist.items())}
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
        for key, bit in F_DOMAIN.items():
            if b in domain_classes[key]:
                f |= bit
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

    # GO RNA/DNA binding among methylated proteins against all reviewed human proteins (clusterProfiler enricher, the
    # four annotations as one job): {key: [methylated with, methylated, other with, other, log2 fold enrichment, p, q]}.
    go_sets = []
    for key in ('rna', 'rna_exp', 'dna', 'dna_exp'):
        for acc in proteome:
            g = go.get(acc)
            codes = (g['rna'] if key.startswith('rna') else g['dna']) if g else set()
            if codes & GO_EXPERIMENTAL if key.endswith('_exp') else codes:
                go_sets.append((key, acc))
    go_test = run_enrichment([Job('go_binding', [a for a in proteome if a in meth_base], proteome, whole_universe=True,
                                  sets=pd.DataFrame(go_sets, columns=['term', 'gene']))], enrich_dir)['go_binding']
    go_enrich, set_tests = {}, []
    for _, x in go_test.iterrows():
        go_enrich[x.ID] = [int(x.a), int(x.n), int(x.c), int(x.N - x.n), sig(math.log2(x.FoldEnrichment)), sig(x.pvalue), sig(x['p.adjust'])]
        set_tests.append(('go_binding', x.ID.replace('_exp', '') + '_binding' + ('_experimental' if x.ID.endswith('_exp') else ''), 'all',
                          int(x.a), int(x.b), int(x.c), int(x.d), float(x.FoldEnrichment), float(x.pvalue), float(x['p.adjust']),
                          'all reviewed human proteins'))

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
            **({'reference': reference} if reference else {}),
            'complexes_tested': {'curated': n_curated, 'predicted': n_predicted},
            # protein-set tests (GO binding, complexes, condensates by confidence): clusterProfiler enricher, log2 fold enrichment
            'set_effect': 'log2_fold_enrichment',
            'methylated_in_proteome': K,
            'prmt_genes': PRMT_GENES,
            'flag_bits': {'prmt': F_PRMT, 'tudor': F_TUDOR, 'rna': F_RNA, 'dna': F_DNA, 'rna_exp': F_RNA_EXP, 'dna_exp': F_DNA_EXP, 'wd40': F_WD40, 'phd': F_PHD},
            'domain_classes': {k: list(v) for k, v in DOMAIN_CLASSES.items()},
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
                'uniprot': 'UniProt reviewed human proteome; InterPro Tudor (IPR002999), WD40 repeat (IPR001680) and PHD finger (IPR001965, IPR019787) members (CC BY 4.0)',
                'go': 'QuickGO/GOA annotations to GO:0003723 and GO:0003677 incl. descendants (CC BY 4.0)',
                'pdbe': 'PDBe-KB interface residues and SIFTS observed segments and chain mapping (EMBL-EBI)',
                'cdcode': f'CD-CODE {cd_version} condensate membership via cd-code.org API',
                'contacts': 'RCSB ModelServer residue surroundings (wwPDB data, CC0)',
                'regions': 'UniProt disordered regions (MobiDB-lite) and domain, repeat and zinc-finger features (CC BY 4.0)',
                'protvar': 'ProtVar predicted interfaces 2024.05.28, AlphaFold-Multimer, residues within 5 Å (EMBL-EBI, CC BY 4.0); served-model list from af2complexes_interfaces_v1',
                'afdb_complexes': 'AlphaFold DB predicted complexes (human heterodimers), API and model files (EMBL-EBI, CC BY 4.0)',
                'rbdmap': 'RBDmap (Castello et al. 2016 Mol Cell, HeLa) Table S1 RBDpep: RNA-bound proteolytic fragments (Oxford Research Archive)',
                'rna_crosslinks': 'pepR-MS (Sha, Kuster & Trendel, bioRxiv 2025, CC BY-NC) Tables S2-S3: RNA-crosslinked residues re-analysed from Kramer 2014, XRNAX 2019, RBS-ID 2020, pRBS-ID 2021 and pepR-MS',
                'dna_crosslinks': 'Trendel et al. 2025 Cell Tables S1-S2: DNA interactome and DNA-crosslinked peptides with localised residues; UniProt DNA-binding regions (CC BY 4.0)',
                'intact': f'IntAct curated features (mutations, binding regions, modified residues; EMBL-EBI, CC BY 4.0). {intact_release}',
            },
            'intact': {'effects': INTACT_EFFECTS, 'regions': INTACT_REGIONS, 'peptide_max': INTACT_PEPTIDE_MAX,
                       'mutation_max': INTACT_MUT_MAX, 'release': intact_release, 'counts': intact_counts},
            'xl_check': xl_check,
            'pred_iface': {'file': Path(args.pred_iface_out).name, 'proteins': len(pred_iface), 'pdockq_min': PRED_TIERS[0]},
            'fields': {
                'proteins': ['accession', 'gene', 'flags', 'ensp'],
                'string': 'per explorer protein: [delta partner index, score, ...] with partner indices ascending',
                'complexes': ['id', 'name', 'predicted', 'members (index or list of alternatives)', 'n', 'methylated', 'fold', 'p', 'q'],
                'iface': ['site', 'mask', 'partners (-1 RNA, -2 DNA, else protein index)', 'pdb ids (first 8)', 'pdb count'],
                'enrich': '[methylated at, methylated not at, other R at, other R not at, log2 OR, p] over structurally observed arginines',
                'go_enrich': '[methylated with term, methylated reviewed, other with term, other reviewed, log2 fold enrichment, p, BH q]',
                'enrich_strata': '{type: {all|ordered|domain: {pooled: [a, b, c, d, log2 OR, p], mh: [log2 MH OR, CMH p, informative proteins]}}}',
                'contacts': ['site', '[partner (-1 RNA, -2 DNA, -3 ligand, -4 unmapped chain, else protein index), pdb, residue labels]', 'methylated residue modelled (2MR/DA2/NMM)'],
                'cd_members': ['protein', 'condensate index (meta.condensates of the sites file)', 'CD-CODE confidence score (0-5)', 'role bits', 'evidence bits', 'publications', 'all PMIDs'],
                'wr': ['PRMT or Tudor protein index', 'its STRING physical partners (indices)'],
                'rg_iface': '{segment|dipeptide|other: {cover: [methylated sites, structurally observed], rna|any: {pooled, mh}}}',
                'geometry': ['site', 'pdb', 'chain', 'Levy class (0 none, 1 support, 2 rim, 3 core)', 'rASA alone', 'rASA in complex', 'nucleic-acid contact (base/backbone/contact/none)', 'bases hydrogen-bonded', 'ASA buried by nucleic acid', 'ASA buried by protein', 'hydrogen bonds to protein'],
                'geometry_summary': '{levy_counts, na_counts: {1 methylated, 0 control}, core|base: {pooled, mh}}',
                'modelled': ['site', 'PDB entries modelling a methylated arginine there', 'residue names'],
                'modelled_gaps': ['protein index', 'UniProt positions modelled as methylarginine but not in the union', 'PDB entries'],
                'cond_conf': '{min stars: [[condensate index, members, methylated, log2 fold enrichment, p, q], ...]}',
                'cond_score_hist': '{CD-CODE confidence score: memberships of reviewed human proteins}',
                'pred': ['site', '[[partner index, pDockQ x 100, 1 if this protein is ProtVar uniprot_id1, served model: 1 same order, 2 reversed, 0 none], ...] (top 12)', 'predicted partners'],
                'gallery': ['site', 'partner index', 'source (0 ProtVar served model, 1 AlphaFold DB complex)', 'score x 100 (ProtVar pDockQ; AlphaFold DB ipTM)', 'pDockQ x 100', 'pLDDT of the arginine', 'lowest pLDDT of the partner residues named', 'minimum heavy-atom distance x 10 (Å)', 'flags (1 salt bridge <= 4 Å, 2 H-bond <= 3.5 Å, 4 cation-pi <= 5 Å, 8 clash < 2.2 Å)', 'partner residues with bond codes (sb, hb, pi)', 'model reference (ProtVar: 1 this protein first, 2 partner first; AlphaFold DB: model id)', 'this protein\'s chain'],
                'xl': ['site', 'RNA crosslink bits at the residue (1 Kramer 2014, 2 XRNAX 2019, 4 RBS-ID 2020, 8 pRBS-ID 2021, 16 pepR-MS 2025, 32 pepR-MS focused site)', 'same bits within 3 residues', 'DNA bits (1 crosslink localised to this residue alone, 2 inside a crosslinked peptide, 4 UniProt DNA-binding region, 8 one of several candidate residues of a crosslink)', '1 if inside an RBDmap RNA-bound fragment', '1 if an RNA crosslink names this residue alone (single localisation)'],
                'xl_prot': ['protein index', 'RNA-crosslinked positions', 'their dataset bits', 'DNA-crosslinked positions (single localisation)', 'DNA-crosslinked peptide spans', 'UniProt DNA-binding regions', '1 if in the Trendel 2025 DNA interactome', 'RBDmap RNA-bound fragments', 'DNA candidate positions (crosslink localised to several residues)'],
                'xl_enrich': '{rna_xl|rna_near|rna_rbdmap|dna_xl|dna_cand|dna_pep|dna_bind (+ _nonrg RNA variants): {pooled, mh, proteins, or: {pooled: [a, b, c, d, OR, lo, hi, p], mh: [OR, lo, hi, p, strata]}}} methylated vs other arginines of proteins with that evidence',
                'xl_types': '{rna_xl|rna_xl_single|rna_rbdmap|dna_xl: {rme2s|rme2a|rme1ip|rme2|rme1|nostate|any: or_summary}} each methylation class vs unmethylated arginines of the same proteins',
                'xl_freq': '{spectra|datasets: {meth: [n, mean, median], other: [...], p (Mann-Whitney)}, multi: or_summary (crosslinked in >= 2 datasets)} among RNA-crosslinked arginines',
                'pred_enrich': '{pDockQ floor: {all|ordered|disordered: {pooled, mh}}} methylated vs other arginines of the same protein at a predicted interface',
                'intact': ['site', '[[kind (0 mutation, 1 binding region <= peptide_max aa, 2 methylated arginine), code (meta.intact effects or regions; 0 for kind 2), start, end, resulting sequence, partners (-1 RNA, -2 DNA, -3 small molecule, -4 other, else protein index), PubMed ids, IntAct interaction ACs (first 5), records], ...]'],
                'intact_prot': ['protein index', '[[start, end, region type bits (1 << meta.intact regions index), partners (first 8), partner count, records], ...] binding regions', '[[start, end, effect code, partners (first 8), partner count, records], ...] mutations of <= mutation_max residues'],
                'intact_enrich': '{mutation: or_summary of methylated vs other arginines, a single-residue mutation changing an interaction vs only no-effect mutations; peptide|peptide_disordered: methylated vs other arginines inside a binding region of <= peptide_max aa (all, or UniProt-disordered arginines only)}',
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
        'cond_score_hist': cond_score_hist,
        'pred': pred_rows,
        'pred_enrich': pred_enrich,
        'gallery': gallery,
        'xl': xl_rows,
        'xl_prot': xl_prot,
        'xl_seq': xl_seq,
        'xl_enrich': xl_enrich,
        'xl_sets': [label for _, label in RNA_XL_SETS],
        'xl_types': xl_types,
        'xl_freq': xl_freq,
        'intact': intact_rows,
        'intact_prot': intact_prot,
        'intact_enrich': intact_enrich,
    }
    # Explorer proteins keep their own gene labels.
    for k, p in enumerate(prots):
        out['proteins'][k][1] = p[1] or out['proteins'][k][1]
    Path(args.out).write_text(json.dumps(out, separators=(',', ':'), ensure_ascii=False))
    Path(args.pred_iface_out).write_text(json.dumps({
        'meta': {'generated': out['meta']['generated'], 'source': out['meta']['sources']['protvar'], 'pdockq_min': PRED_TIERS[0],
                 'fields': ['start', 'end', 'partners (most along the run)', 'best pDockQ x100']},
        'prot': pred_iface}, separators=(',', ':')))
    log(f'wrote {args.pred_iface_out}: predicted interface runs for {len(pred_iface):,} proteins '
        f'({Path(args.pred_iface_out).stat().st_size / 1e6:.2f} MB)')
    if args.tables_out:
        write_release_tables(Path(args.tables_out), {
            'prots': prots, 'sites': sites, 'vocab': vocab, 'gene_of': gene_of, 'iface_detail': iface_detail, 'xl_rows': xl_rows,
            'pred_all': pred_all, 'gallery': gallery, 'intact_rows': intact_rows, 'partners': partners, 'domain_classes': domain_classes,
            'meth_base': meth_base, 'cx_rows': cx_rows, 'iface_residues': iface_residues, 'pred_iface': pred_iface, 'xl_prot': xl_prot,
            'intact_prot': intact_prot, 'tests': release_tests, 'set_tests': set_tests, 'complexes': cx_release, 'go': go, 'xl_seq': xl_seq_all,
            'proteome': proteome})
    log(f'wrote {args.out}: {len(vocab):,} proteins, {sum(len(r) // 2 for r in string_rows):,} partner links, '
        f'{len(cx_rows):,} complexes with methylated members, {len(site_rows):,} sites at PDB interfaces '
        f'({Path(args.out).stat().st_size / 1e6:.2f} MB)')


if __name__ == '__main__':
    main()

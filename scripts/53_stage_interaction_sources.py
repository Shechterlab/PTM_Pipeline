#!/usr/bin/env python3
"""Stage the provenance of the interaction records at methylarginines: which experiment, which method, which
species, which paper, and the state of the site in the experiment.

The interaction inputs (scripts/37a) say that a partner contacts a residue in some PDB entries (PDBe-KB) or that an
IntAct interaction maps a binding region, a mutation or a methylated residue there; they do not say how the
interaction was found, whether the site was methylated in the experiment, what else was modified, in which species,
or which paper reported it. This script looks that up for every record at a site of a release (its sites and its
lower-confidence sites), so the build can say it:

  RCSB PDB Data API (GraphQL, data.rcsb.org; PDB data CC0 1.0), for every PDB entry with a PDBe-KB interface at a
      site or an arginine contact in pdb_arginine_contacts_*: title, method, resolution, release date, primary
      citation; every polymer entity (description, type, source organisms, UniProt alignment); non-polymer ligands;
      and, for each entity aligned to a site's protein, the residue modelled at the site and the other modified
      residues (UniProt numbering where the entity aligns, the chemical component otherwise).
  The same API and the RCSB ModelServer (models.rcsb.org), for every PDB entry modelling a methylarginine (2MR, DA2,
      NMM; pdb_methylarginine_residues_*): which sites of the release each one is, by the entity's UniProt alignment or
      else by sequence (orthologues, synthetic peptides and the other histone variants sharing a tail), and what the
      methylated side chain touches (residues of other chains within 4 Å). scripts/37 finds these only through SIFTS,
      so it misses peptides without a human UniProt mapping, such as reader domains bound to methylated histone tails
      of mouse or frog sequence.
  IntAct web service (www.ebi.ac.uk/intact/ws, CC BY 4.0), for every interaction scripts/37 places at a site (binding
      regions of at most 20 residues, mutations of at most 3, methylated residues; the staged feature tables): the
      binary records with both participants (identifier, type, species), their features, detection method,
      interaction type, host organism, negative flag, MI score and publication (PubMed, DOI, authors, journal,
      title, year).
  Ontology Lookup Service (www.ebi.ac.uk/ols4; PSI-MI and PSI-MOD, CC BY 4.0), for the identifier and name of every
      detection method, interaction type and feature label the IntAct answers use.
  NCBI E-utilities (esummary), for the title, journal, year and DOI of the papers of those IntAct records.
  UniProt REST, for the gene and organism of the partners that are not reviewed human entries (other species).

Every answer is cached as JSON under --cache (not in Git; its files are listed in raw_cache_files.tsv). Outputs go
to --outdir; --install copies them into inputs/interactions/ and puts their rows in inputs/MANIFEST.tsv and
inputs/interactions/raw_cache_files.tsv.

    python scripts/53_stage_interaction_sources.py --release releases/v1.4.2 --cache PTM_data/interactions_cache \
        --outdir /tmp/sources [--install]
"""
from __future__ import annotations

import argparse
import bisect
import datetime as dt
import glob
import hashlib
import json
import math
import re
import shutil
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import md5_file, refresh_manifest_row, update_manifest, upsert_tsv_rows  # noqa: E402

RCSB = 'https://data.rcsb.org/graphql'
MODELSERVER = 'https://models.rcsb.org/v1/{pdb}/residueSurroundings'
INTACT = 'https://www.ebi.ac.uk/intact/ws/interaction/findInteractionWithFacet'
OLS = 'https://www.ebi.ac.uk/ols4/api'
ESUMMARY = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
EFETCH = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
# an arginine the abstract names by number: R3, Arg-3, arginine 3, H4R3 (histone marks), R3me2s
ARG_NUMBER = re.compile(r'(?:(?<=H[1-4])|(?<=H2[AB])|(?<![A-Za-z0-9]))(?:R|Arg)\s?-?(\d{1,4})(?![\d.]*\s?(?:nM|µM|uM|mM|kDa))(?=me|[ACDEFGHIKLMNPQSTVWY]?(?:[^\w]|$))|\barginine[\s-]+(\d{1,4})\b')
UNIPROT = 'https://rest.uniprot.org/uniprotkb/accessions'
PEPTIDE_MAX, MUTATION_MAX = 20, 3      # as scripts/37 (INTACT_PEPTIDE_MAX, INTACT_MUT_MAX)
ENTRY_QUERY = '''{entries(entry_ids:%s){rcsb_id struct{title} exptl{method} rcsb_entry_info{resolution_combined}
 rcsb_accession_info{initial_release_date}
 rcsb_primary_citation{pdbx_database_id_PubMed pdbx_database_id_DOI year rcsb_authors journal_abbrev title}
 nonpolymer_entities{pdbx_entity_nonpoly{comp_id name}}
 polymer_entities{rcsb_id rcsb_polymer_entity{pdbx_description} entity_poly{type pdbx_seq_one_letter_code}
  rcsb_entity_source_organism{ncbi_scientific_name ncbi_taxonomy_id}
  polymer_entity_instances{rcsb_polymer_entity_instance_container_identifiers{asym_id auth_asym_id}}
  rcsb_polymer_entity_align{reference_database_name reference_database_accession aligned_regions{entity_beg_seq_id ref_beg_seq_id length}}}}}'''
CHEM_QUERY = '{chem_comps(comp_ids:%s){chem_comp{id name type one_letter_code mon_nstd_parent_comp_id}}}'
STANDARD = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
            'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
            'U': 'SEC', 'O': 'PYL'}
METHYL_R = {'2MR', 'DA2', 'NMM'}       # the omega-methylarginines of the PDB, as scripts/37 (METHYL_R)
WINDOW, MIN_MATCH = 10, 8              # residues each side of a methylarginine matched to a human sequence; shortest match
MIN_DISTINCT = 5                       # residue types a match needs: low-complexity stretches (RGRGRGRG) recur by chance
CONTACT_CUTOFF = 4.0                   # Å from the side chain, as scripts/37 (CONTACT_CUTOFF)
FLANK = 5                              # residues each side of a site checked against the human sequence
BACKBONE = {'N', 'CA', 'C', 'O', 'OXT'}
WATERS = {'HOH', 'WAT', 'DOD'}


def rd(path) -> pd.DataFrame:
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


def http(url: str, data: bytes | None = None, headers: dict | None = None, tries: int = 4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers=headers or {}, method='POST' if data is not None else 'GET')
            with urllib.request.urlopen(req, timeout=180) as r:
                return json.load(r)
        except Exception as e:  # noqa: BLE001
            if k == tries - 1:
                raise
            print(f'  retry {k + 1} after {e}', flush=True)
            time.sleep(3 * (k + 1))


def cached(cache: Path, name: str, fetch):
    """The JSON answer cached at cache/name, fetched once."""
    path = cache / name
    if path.exists():
        return json.loads(path.read_text(encoding='utf-8'))
    ans = fetch()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(ans, sort_keys=True), encoding='utf-8')
    return ans


def site_positions(release: Path) -> dict[str, set[int]]:
    """accession -> positions of the release's sites and lower-confidence sites."""
    out = defaultdict(set)
    for f in sorted(release.glob('human_methylarginine_sites_v*.tsv')) + sorted(release.glob('human_methylarginine_lower_confidence_sites_v*.tsv')):
        t = rd(f)
        for a, p in zip(t['uniprot_accession'].str.split('-').str[0], t['position']):
            out[a].add(int(p))
    return out


def spans(ranges: str) -> list[tuple[int, int]]:
    out = []
    for part in ranges.split(','):
        a, _, b = part.strip().partition('-')
        try:
            out.append((int(a), int(b or a)))
        except ValueError:
            pass
    return out


# ---- PDB ---------------------------------------------------------------------------------------------------------

def pdb_entries(inputs: Path, S) -> list[str]:
    ent = set()
    for f in sorted(glob.glob(str(inputs / 'interactions' / 'pdbe_kb_interface_residues_*.tsv'))):
        t = rd(f)
        t = t[t['uniprot_accession'].isin(S)]
        for a, s, e, pe in zip(t['uniprot_accession'], t['start'].astype(int), t['end'].astype(int), t['pdb_entries']):
            if any(s <= p <= e for p in S[a]):
                ent |= {x.strip().upper() for x in pe.replace(';', ',').split(',') if x.strip()}
    for f in sorted(glob.glob(str(inputs / 'interactions' / 'pdb_arginine_contacts_*.tsv'))) + \
            sorted(glob.glob(str(inputs / 'interactions' / 'pdb_methylarginine_residues_*.tsv'))):
        ent |= set(rd(f)['pdb_id'].str.upper())
    return sorted(ent)


def read_fasta(path: Path) -> dict[str, str]:
    out, acc = {}, None
    for line in path.read_text(encoding='utf-8').splitlines():
        if line.startswith('>'):
            acc = line[1:].split()[0]
            out[acc] = ''
        elif acc:
            out[acc] += line.strip()
    return out


def one_letter(tok: str, chem: dict) -> str:
    if len(tok) == 1:
        return tok.upper()
    one = (chem.get(tok) or ('', '', '', ''))[2]
    return one if len(one) == 1 and one.isalpha() and one != 'X' else 'X'


class Proteome:
    """The sequences of the proteins with sites, joined so one search covers them all."""
    def __init__(self, seqs: dict[str, str]):
        self.accs = sorted(seqs)
        self.text = '|'.join(seqs[a] for a in self.accs)
        self.starts, at = [], 0
        for a in self.accs:
            self.starts.append(at)
            at += len(seqs[a]) + 1

    def find(self, word: str) -> list[tuple[str, int]]:
        """(accession, 0-based offset) of every exact occurrence."""
        out, i = [], self.text.find(word)
        while i >= 0:
            j = bisect.bisect_right(self.starts, i) - 1
            out.append((self.accs[j], i - self.starts[j]))
            i = self.text.find(word, i + 1)
        return out


def sequence_matches(tok: list[str], k: int, chem: dict, proteome: Proteome, S) -> list[tuple[str, int, int, str]]:
    """Human site positions a methylarginine at entity position k (1-based) falls on by sequence: the longest stretch of
    the residues within WINDOW of it (modified residues read as their parent, none without one) that holds it and is
    found exactly in the sequence of a protein with sites, at least MIN_MATCH long and of at least MIN_DISTINCT residue
    types. -> [(accession, position, match length, matched stretch)]"""
    letters = [one_letter(t, chem) for t in tok]
    lo, hi = k - 1, k - 1
    while lo > max(0, k - 1 - WINDOW) and letters[lo - 1] != 'X':
        lo -= 1
    while hi < min(len(letters) - 1, k - 1 + WINDOW) and letters[hi + 1] != 'X':
        hi += 1
    for n in range(hi - lo + 1, MIN_MATCH - 1, -1):
        out = []
        for a in range(max(lo, k - n), min(k - 1, hi - n + 1) + 1):   # stretches [a, a + n) that hold index k - 1
            word = ''.join(letters[a:a + n])
            if len(set(word)) < MIN_DISTINCT:
                continue
            for acc, off in proteome.find(word):
                pos = off + (k - 1 - a) + 1
                if pos in S.get(acc, ()):
                    out.append((acc, pos, n, word))
        if out:
            return sorted(set(out))
    return []


def site_context(letters: list[str], k: int, useq: str, p: int, flank: int = FLANK) -> tuple[int, str]:
    """Check the entity position k the RCSB alignment gives for UniProt position p against the residues around it.
    -> (entity position, matches | realigned | differs): matches when the FLANK residues on each side (those both
    sequences have) are identical; realigned when they are not but the human residues around p occur once in the
    entity, the site at another position (an isoform the alignment does not follow); differs otherwise."""
    def same(kk):
        n = 0
        for d in range(-flank, flank + 1):
            i, j = kk - 1 + d, p - 1 + d
            if d == 0 or not (0 <= i < len(letters) and 0 <= j < len(useq)):
                continue
            if letters[i] != useq[j] and letters[i] != 'X':
                return False
            n += 1
        return n >= flank
    if not useq or same(k):
        return k, 'matches' if useq else 'differs'
    hits = [kk for kk in range(1, len(letters) + 1) if same(kk)]
    return (hits[0], 'realigned') if len(hits) == 1 else (k, 'differs')


def parse_cif_atoms(text: str) -> list[dict]:
    """The atom_site loop of an mmCIF answer, as scripts/37 (parse_cif_atoms)."""
    cols, atoms, in_loop = [], [], False
    for line in text.splitlines():
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


def residue_contacts(cache: Path, pdb: str, chain: str, label: int) -> dict:
    """Residues of other chains within CONTACT_CUTOFF of the side chain of residue label (label_seq_id) of an author
    chain, first model of the asymmetric unit (RCSB ModelServer residueSurroundings), cached as JSON:
    {resn, hits: [[author chain, label_seq_id, author residue number, residue name, distance]]}."""
    def fetch():
        url = MODELSERVER.format(pdb=pdb.lower()) + '?' + urllib.parse.urlencode(
            {'auth_asym_id': chain, 'label_seq_id': label, 'radius': CONTACT_CUTOFF + 1.5, 'model_nums': 1, 'encoding': 'cif'})
        text = ''
        for k in range(4):
            try:
                with urllib.request.urlopen(urllib.request.Request(url), timeout=120) as r:
                    text = r.read().decode('utf-8', 'replace')
                break
            except urllib.error.HTTPError as e:
                if e.code in (400, 404):
                    break
                time.sleep(3 * (k + 1))
            except Exception:  # noqa: BLE001 - retried
                if k == 3:
                    raise
                time.sleep(3 * (k + 1))
        atoms = parse_cif_atoms(text)
        me = [a for a in atoms if a.get('auth_asym_id') == chain and a.get('label_seq_id') == str(label)]
        side = [(float(a['Cartn_x']), float(a['Cartn_y']), float(a['Cartn_z'])) for a in me if a.get('label_atom_id') not in BACKBONE]
        seen = {}
        for a in atoms if side else []:
            ch, comp = a.get('auth_asym_id'), a.get('label_comp_id')
            if ch == chain or comp in WATERS:
                continue
            d = min(math.dist((float(a['Cartn_x']), float(a['Cartn_y']), float(a['Cartn_z'])), q) for q in side)
            if d <= CONTACT_CUTOFF:
                key = (ch, a.get('label_seq_id'), a.get('auth_seq_id'), comp)
                seen[key] = min(d, seen.get(key, d))
        return {'resn': me[0].get('label_comp_id') if me else '',
                'hits': [[*key, round(d, 2)] for key, d in sorted(seen.items(), key=lambda kv: kv[1])]}
    return cached(cache, f'modelserver_methyl/{pdb.lower()}_{chain}_{label}.json', fetch)


def tokens(seq: str) -> list[str]:
    """One-letter code with (XYZ) modified residues -> one token per residue."""
    return [a or b for a, b in re.findall(r'\(([A-Za-z0-9]+)\)|([A-Za-z])', seq.replace('\n', ''))]


def stage_pdb(inputs: Path, cache: Path, S, seqs: dict[str, str]) -> tuple[dict[str, pd.DataFrame], dict]:
    ids = pdb_entries(inputs, S)
    proteome = Proteome(seqs)
    print(f'PDB: {len(ids):,} entries', flush=True)
    entries = []
    qkey = hashlib.md5(ENTRY_QUERY.encode()).hexdigest()[:6]
    for i in range(0, len(ids), 50):
        batch = ids[i:i + 50]
        ans = cached(cache, f'rcsb/entries_{qkey}_{hashlib.md5(",".join(batch).encode()).hexdigest()[:12]}.json',
                     lambda b=batch: http(RCSB, json.dumps({'query': ENTRY_QUERY % json.dumps(b)}).encode(), {'Content-Type': 'application/json'}))
        entries += [e for e in (ans.get('data') or {}).get('entries') or [] if e]
        if i % 1000 == 0:
            print(f'  {i + len(batch):,} of {len(ids):,}', flush=True)
    codes = sorted({t for e in entries for pe in e.get('polymer_entities') or [] for t in tokens((pe.get('entity_poly') or {}).get('pdbx_seq_one_letter_code') or '') if len(t) > 1})
    chem = {}
    for i in range(0, len(codes), 200):
        b = codes[i:i + 200]
        ans = cached(cache, f'rcsb/chem_{hashlib.md5(",".join(b).encode()).hexdigest()[:12]}.json',
                     lambda b=b: http(RCSB, json.dumps({'query': CHEM_QUERY % json.dumps(b)}).encode(), {'Content-Type': 'application/json'}))
        for c in (ans.get('data') or {}).get('chem_comps') or []:
            if c and c.get('chem_comp'):
                cc = c['chem_comp']
                par = cc.get('mon_nstd_parent_comp_id') or ''
                chem[cc['id']] = (cc.get('name') or '', ','.join(par) if isinstance(par, list) else par, cc.get('one_letter_code') or '', cc.get('type') or '')
    rows, ents, sites, methyl, jobs = [], [], [], [], []
    for e in entries:
        pid = e['rcsb_id'].lower()
        chain_entity = {(i.get('rcsb_polymer_entity_instance_container_identifiers') or {}).get('auth_asym_id'): pe['rcsb_id'].split('_')[-1]
                        for pe in e.get('polymer_entities') or [] for i in pe.get('polymer_entity_instances') or [] if i}
        c = e.get('rcsb_primary_citation') or {}
        authors = c.get('rcsb_authors') or []
        res = (e.get('rcsb_entry_info') or {}).get('resolution_combined') or []
        lig = sorted({(n['pdbx_entity_nonpoly']['comp_id'], n['pdbx_entity_nonpoly']['name']) for n in e.get('nonpolymer_entities') or [] if n and n.get('pdbx_entity_nonpoly')})
        rows.append({'pdb_id': e['rcsb_id'].lower(), 'title': (e.get('struct') or {}).get('title') or '',
                     'method': '; '.join(x['method'] for x in e.get('exptl') or []), 'resolution_a': f'{res[0]:.2f}' if res else '',
                     'release_date': ((e.get('rcsb_accession_info') or {}).get('initial_release_date') or '')[:10],
                     'pubmed': str(c.get('pdbx_database_id_PubMed') or ''), 'doi': c.get('pdbx_database_id_DOI') or '',
                     'first_author': authors[0] if authors else '', 'authors': str(len(authors)) if authors else '', 'year': str(c.get('year') or ''),
                     'journal': c.get('journal_abbrev') or '', 'citation_title': c.get('title') or '',
                     'ligands': '; '.join(f'{k} {v}' for k, v in lig)})
        for pe in e.get('polymer_entities') or []:
            orgs = pe.get('rcsb_entity_source_organism') or []
            al = [(a['reference_database_accession'], r['entity_beg_seq_id'], r['ref_beg_seq_id'], r['length'])
                  for a in pe.get('rcsb_polymer_entity_align') or [] if a and a.get('reference_database_name') == 'UniProt'
                  for r in a.get('aligned_regions') or []]
            ent_id = pe['rcsb_id'].split('_')[-1]
            ents.append({'pdb_id': e['rcsb_id'].lower(), 'entity_id': ent_id, 'type': (pe.get('entity_poly') or {}).get('type') or '',
                         'description': (pe.get('rcsb_polymer_entity') or {}).get('pdbx_description') or '',
                         'uniprot_accessions': '; '.join(sorted({a for a, *_ in al})),
                         'organisms': '; '.join(sorted({o['ncbi_scientific_name'] for o in orgs if o.get('ncbi_scientific_name')})),
                         'taxids': '; '.join(sorted({str(o['ncbi_taxonomy_id']) for o in orgs if o.get('ncbi_taxonomy_id')}))})
            tok = tokens((pe.get('entity_poly') or {}).get('pdbx_seq_one_letter_code') or '')
            if not tok:
                continue
            # entity position (1-based) -> (accession, UniProt position)
            ref = {}
            for acc, eb, rb, n in al:
                for k in range(n):
                    ref[eb + k] = (acc, rb + k)
            mods = []
            for k, t in enumerate(tok, 1):
                if len(t) > 1 and t != 'UNK':
                    name, parent, one, _ = chem.get(t, ('', '', '', ''))
                    at = ref.get(k)
                    mods.append((t, parent, at, name))
            letters = [one_letter(t, chem) for t in tok]
            for acc in sorted({a for a, *_ in al} & set(S)):
                for p in sorted(S[acc]):
                    ks = [k for k, (a, q) in ref.items() if a == acc and q == p]
                    if not ks:
                        continue
                    k, check = site_context(letters, ks[0], seqs.get(acc, ''), p)
                    t = tok[k - 1] if k - 1 < len(tok) else ''
                    name, parent, one, _ = chem.get(t, ('', '', '', '')) if len(t) > 1 else ('', '', '', '')
                    code = STANDARD.get(t, t) if len(t) == 1 else t
                    others = [f'{(chem.get(m, ("", "", "", ""))[2] or "?")}{at[1]} {m}' if at and at[0] == acc else f'{m}'
                              for m, par, at, nm in mods if not (at and at[0] == acc and at[1] == p)]
                    sites.append({'pdb_id': e['rcsb_id'].lower(), 'entity_id': ent_id, 'uniprot_accession': acc, 'position': str(p),
                                  'residue': code, 'residue_name': name, 'residue_parent': parent or (code if len(t) == 1 else ''),
                                  'entity_position': str(k), 'context_check': check, 'other_modified_residues': '; '.join(others)})
            # A methylarginine the entity models: the human sites it is, by the UniProt alignment of the entity or else
            # by sequence (orthologues, synthetic peptides and peptides without an alignment)
            own = sorted({a for a, *_ in al})
            chains = sorted(c for c, x in chain_entity.items() if x == ent_id and c)
            for k, t in enumerate(tok, 1):
                if t not in METHYL_R:
                    continue
                at = ref.get(k)
                hits = [(at[0], at[1], 0, '')] if at and at[1] in S.get(at[0], ()) and site_context(letters, k, seqs.get(at[0], ''), at[1]) == (k, 'matches') else []
                # the same stretch in other proteins with sites (histone variants share their tails), or the human
                # site of an orthologue, a synthetic peptide or an unaligned one
                hits += [h for h in sequence_matches(tok, k, chem, proteome, S) if not hits or h[:2] != hits[0][:2]]
                for acc, p, n, win in hits:
                    methyl.append({'pdb_id': pid, 'entity_id': ent_id, 'entity_position': str(k), 'residue': t,
                                   'uniprot_accession': acc, 'position': str(p),
                                   'matched_by': 'UniProt alignment' if not n else f'sequence identity over {n} residues',
                                   'matched_window': win, 'entity_uniprot': '; '.join(own),
                                   'entity_organisms': '; '.join(sorted({o['ncbi_scientific_name'] for o in orgs if o.get('ncbi_scientific_name')})),
                                   'chains': '; '.join(chains)})
                    jobs += [(pid, c, k) for c in chains[:2]]
    # what each methylarginine at a site touches: the residues of other chains within CONTACT_CUTOFF of its side chain
    ent_of = {(e['rcsb_id'].lower(), (i.get('rcsb_polymer_entity_instance_container_identifiers') or {}).get('auth_asym_id')): pe['rcsb_id'].split('_')[-1]
              for e in entries for pe in e.get('polymer_entities') or [] for i in pe.get('polymer_entity_instances') or [] if i}
    contacts = []
    for pid, ch, k in sorted(set(jobs)):
        res = residue_contacts(cache, pid, ch, k)
        by = defaultdict(list)
        for hch, hlab, hseq, comp, d in res['hits']:
            by[(hch, ent_of.get((pid, hch), ''))].append((comp, hseq, d))
        for (hch, hent), hs in sorted(by.items()) or [(('', ''), [])]:
            contacts.append({'pdb_id': pid, 'chain': ch, 'entity_position': str(k), 'residue_modelled': res['resn'],
                             'partner_chain': hch, 'partner_entity_id': hent,
                             'partner_residues': '; '.join(sorted({f'{c}{n}' for c, n, _ in hs}, key=lambda x: (re.sub(r'-?\d+$', '', x), int(m.group()) if (m := re.search(r'-?\d+$', x)) else 0))),
                             'closest_distance_a': f'{min(d for *_, d in hs):.2f}' if hs else ''})
    pdb = pd.DataFrame(rows).sort_values('pdb_id')
    ent = pd.DataFrame(ents).sort_values(['pdb_id', 'entity_id'])
    me = pd.DataFrame(methyl).drop_duplicates().sort_values(['uniprot_accession', 'position', 'pdb_id', 'entity_id', 'entity_position'],
                                                           key=lambda s: s.astype(int) if s.name in ('position', 'entity_position') else s)
    mc = pd.DataFrame(contacts).sort_values(['pdb_id', 'chain', 'entity_position', 'partner_chain'], key=lambda s: s.astype(int) if s.name == 'entity_position' else s)
    st = pd.DataFrame(sites).drop_duplicates().sort_values(['uniprot_accession', 'position', 'pdb_id', 'entity_id'], key=lambda s: s.map(lambda x: int(x) if str(x).isdigit() else x) if s.name == 'position' else s)
    info = {'entries_asked': len(ids), 'entries_answered': len(rows), 'entities': len(ent), 'site_residues': len(st),
            'chemical_components': len(chem), 'site_residue_not_arginine': int((~st['residue_parent'].isin(['ARG', 'R'])).sum() if len(st) else 0),
            'methylated_site_rows': len(me), 'methylated_site_entries': int(me['pdb_id'].nunique()) if len(me) else 0,
            'methylated_sites': int(me[['uniprot_accession', 'position']].drop_duplicates().shape[0]) if len(me) else 0,
            'methylated_by_sequence': int((me['matched_by'] != 'UniProt alignment').sum()) if len(me) else 0,
            'methylated_residues_read': len(set(jobs)), 'methylated_residues_with_contacts': int(mc.loc[mc['partner_chain'] != '', ['pdb_id', 'chain', 'entity_position']].drop_duplicates().shape[0]) if len(mc) else 0}
    return {'pdb_entries': pdb, 'pdb_entities': ent, 'pdb_site_residues': st, 'pdb_methylated_sites': me, 'pdb_methylated_site_contacts': mc}, info


# ---- IntAct --------------------------------------------------------------------------------------------------------

def intact_acs(inputs: Path, S) -> dict[str, set[str]]:
    """interaction ac -> the kinds of feature that place it at a site, by the rules of scripts/37."""
    acs = defaultdict(set)
    for name, kind, mx in (('binding_regions', 'binding_region', PEPTIDE_MAX), ('mutations', 'mutation', MUTATION_MAX), ('methylarginines', 'methylated', 1)):
        for f in sorted(glob.glob(str(inputs / 'interactions' / f'intact_{name}_human_*.tsv'))):
            t = rd(f)
            t = t[t['uniprot_accession'].isin(S) & t['sequence_check'].isin(['matches', 'no curated sequence'])]
            for a, rg, ac in zip(t['uniprot_accession'], t['ranges'], t['interaction_ac']):
                if any(y - x + 1 <= mx and any(x <= p <= y for p in S[a]) for x, y in spans(rg)):
                    acs[ac].add(kind)
    return acs


ANN = re.compile(r'^([a-z -]+) \((.*)\)$', re.S)


def stage_intact(inputs: Path, cache: Path, S) -> tuple[pd.DataFrame, set, dict]:
    acs = intact_acs(inputs, S)
    ids = sorted(acs)
    print(f'IntAct: {len(ids):,} interactions', flush=True)
    recs = []
    for i in range(0, len(ids), 40):           # the service answers 500 to longer queries
        b = ids[i:i + 40]
        url = INTACT + '?' + urllib.parse.urlencode({'query': ' '.join(b), 'page': 0, 'pageSize': 200})
        ans = cached(cache, f'intact_ws/batch_{hashlib.md5(",".join(b).encode()).hexdigest()[:12]}.json', lambda u=url: http(u, b''))
        recs += [c for c in ans['data']['content'] if c['ac'] in acs]
    labels, rows = set(), []
    for c in recs:
        ann = {}
        for a in c.get('allAnnotations') or []:
            m = ANN.match(a.strip())
            if m:
                ann.setdefault(m.group(1), m.group(2))
        fa = (c.get('firstAuthor') or '').strip()
        year = (re.search(r'\((\d{4})\)', fa) or [None, ''])[1]
        doi = next((p.split(' ')[0] for p in c.get('publicationIdentifiers') or [] if p.endswith('(doi)')), '')
        for side in 'AB':
            for f in c.get('features' + side) or []:
                labels.add(f.split(':')[0].strip())
        rows.append({'interaction_ac': c['ac'], 'binary_id': str(c.get('binaryInteractionId') or ''), 'placed_by': '; '.join(sorted(acs[c['ac']])),
                     'detection_method_mi': c.get('detectionMethodMIIdentifier') or '', 'interaction_type_mi': c.get('typeMIIdentifier') or '',
                     'host_organism': c.get('hostOrganism') or '', 'host_taxid': str(c.get('hostOrganismTaxId') or ''),
                     'negative': 'true' if c.get('negative') else 'false', 'mi_score': str(c.get('intactMiscore') or ''),
                     'source_database': c.get('sourceDatabase') or '',
                     **{f'{k}_{s.lower()}': v for s in 'AB' for k, v in (
                         ('participant', c.get('uniqueId' + s) or ''), ('participant_type', c.get('type' + s) or ''),
                         ('participant_name', c.get('molecule' + s) or ''), ('participant_description', c.get('description' + s) or ''),
                         ('participant_taxid', str(c.get('taxId' + s) or '')), ('participant_species', c.get('species' + s) or ''),
                         ('experimental_role', c.get('experimentalRole' + s) or ''), ('biological_role', c.get('biologicalRole' + s) or ''),
                         ('features', ' | '.join(c.get('features' + s) or [])))},
                     'pubmed': c.get('publicationPubmedIdentifier') or '', 'doi': doi, 'first_author': re.sub(r'\s*\(\d{4}\)\s*$', '', fa).strip(),
                     'year': year, 'journal': ann.get('journal', ''), 'paper_title': ann.get('title', ''),
                     'authors': (ann.get('author-list', '') or '').strip()})
    t = pd.DataFrame(rows)
    missing = sorted(set(ids) - set(t['interaction_ac'])) if len(t) else ids
    info = {'interactions_asked': len(ids), 'interactions_answered': int(t['interaction_ac'].nunique()) if len(t) else 0,
            'binary_records': len(t), 'not_answered': missing[:50], 'negative': int((t['negative'] == 'true').sum()) if len(t) else 0}
    return t, labels, info


# ---- Ontology terms -----------------------------------------------------------------------------------------------

def ols_term(cache: Path, obo: str):
    ont = obo.split(':')[0].lower()
    ans = cached(cache, f'ols/{obo.replace(":", "_")}.json',
                 lambda: http(f'{OLS}/ontologies/{ont}/terms?' + urllib.parse.urlencode({'obo_id': obo})))
    t = ((ans.get('_embedded') or {}).get('terms') or [{}])[0]
    return t.get('label') or ''


def ols_label(cache: Path, label: str):
    """A feature's short label -> (obo id, name, match): the PSI-MI or PSI-MOD term with that name or synonym (IntAct
    labels are the terms' short labels, e.g. nonome2arg for MOD:00077), else the best text match."""
    slug = re.sub(r'[^A-Za-z0-9]+', '_', label)[:60]
    docs = []
    for field in ('synonym', 'label'):
        ans = cached(cache, f'ols/exact_{field}_{slug}.json',
                     lambda f=field: http(f'{OLS}/search?' + urllib.parse.urlencode({'q': label, 'ontology': 'mod,mi', 'exact': 'true', 'queryFields': f,
                                                                                     'fieldList': 'obo_id,label,synonym', 'rows': 5})))
        docs = [d for d in (ans.get('response') or {}).get('docs') or []
                if label.lower() in [s.lower() for s in (d.get('synonym') or [])] + [(d.get('label') or '').lower()]]
        if docs:
            return docs[0].get('obo_id') or '', docs[0].get('label') or '', 'exact'
    ans = cached(cache, f'ols/search_{slug}.json',
                 lambda: http(f'{OLS}/search?' + urllib.parse.urlencode({'q': label, 'ontology': 'mod,mi', 'rows': 5})))
    docs = (ans.get('response') or {}).get('docs') or []
    d = (docs or [{}])[0]
    return d.get('obo_id') or '', d.get('label') or '', 'nearest' if docs else 'none'


def stage_terms(cache: Path, ia: pd.DataFrame, labels: set) -> pd.DataFrame:
    rows = []
    for mi in sorted(set(ia['detection_method_mi']) | set(ia['interaction_type_mi']) - {''}):
        if mi:
            rows.append({'kind': 'detection method' if mi in set(ia['detection_method_mi']) else 'interaction type', 'label': '', 'obo_id': mi,
                         'name': ols_term(cache, mi), 'match': 'identifier'})
    for lab in sorted(labels - {''}):
        obo, name, how = ols_label(cache, lab)
        rows.append({'kind': 'feature', 'label': lab, 'obo_id': obo, 'name': name, 'match': how})
    return pd.DataFrame(rows)


# ---- Papers ------------------------------------------------------------------------------------------------------

def abstract_arginines(cache: Path, pmids: list[str]) -> dict[str, str]:
    """PubMed id -> the arginine numbers its abstract names ('3; 17'), from NCBI efetch. Only the numbers are kept
    (the abstracts are the publishers'); a paper whose abstract names none gets ''."""
    out = {}
    for i in range(0, len(pmids), 200):
        b = pmids[i:i + 200]
        def fetch(b=b):
            url = EFETCH + '?' + urllib.parse.urlencode({'db': 'pubmed', 'id': ','.join(b), 'rettype': 'xml', 'retmode': 'xml'})
            with urllib.request.urlopen(url, timeout=180) as r:
                xml = r.read().decode('utf-8', 'replace')
            found = {}
            for art in re.findall(r'<PubmedArticle>.*?</PubmedArticle>', xml, re.S):
                pm = (re.search(r'<PMID[^>]*>(\d+)</PMID>', art) or [None, ''])[1]
                text = re.sub(r'<[^>]+>', ' ', ' '.join(re.findall(r'<(?:ArticleTitle|AbstractText)[^>]*>(.*?)</(?:ArticleTitle|AbstractText)>', art, re.S)))
                found[pm] = sorted({int(a or b) for a, b in ARG_NUMBER.findall(text)})
            return found
        ans = cached(cache, f'ncbi/arginines_{hashlib.md5(",".join(b).encode()).hexdigest()[:12]}.json', fetch)
        out.update({k: '; '.join(map(str, v)) for k, v in ans.items()})
        time.sleep(0.4)
    return out


def stage_papers(cache: Path, pmids: list[str]) -> pd.DataFrame:
    """Title, first author, journal, year and DOI of each PubMed id (NCBI E-utilities esummary): IntAct's answers give
    the paper's id and authors but not its title; and the arginines its title and abstract name by number (efetch)."""
    named = abstract_arginines(cache, pmids)
    rows = []
    for i in range(0, len(pmids), 200):
        b = pmids[i:i + 200]
        ans = cached(cache, f'ncbi/esummary_{hashlib.md5(",".join(b).encode()).hexdigest()[:12]}.json',
                     lambda b=b: http(ESUMMARY + '?' + urllib.parse.urlencode({'db': 'pubmed', 'retmode': 'json', 'id': ','.join(b)})))
        res = ans.get('result') or {}
        for k in res.get('uids') or []:
            d = res[k]
            rows.append({'pubmed': k, 'title': (d.get('title') or '').strip(), 'first_author': d.get('sortfirstauthor') or '',
                         'journal': d.get('source') or '', 'year': (d.get('pubdate') or '')[:4],
                         'doi': next((a['value'] for a in d.get('articleids') or [] if a.get('idtype') == 'doi'), ''), 'arginines_named': named.get(k, '')})
        time.sleep(0.4)   # NCBI asks for at most three requests a second without a key
    return pd.DataFrame(rows, columns=['pubmed', 'title', 'first_author', 'journal', 'year', 'doi', 'arginines_named']).sort_values('pubmed', key=lambda s: s.astype(int))


# ---- Names of the partners outside the human reviewed set ----------------------------------------------------------

def stage_names(cache: Path, accs: list[str]) -> pd.DataFrame:
    """Primary gene, organism and recommended name of UniProt accessions (UniProt REST): the partners of other species,
    and unreviewed entries, which the release's reference does not name."""
    rows = []
    for i in range(0, len(accs), 200):
        b = accs[i:i + 200]
        def fetch(b=b):
            url = UNIPROT + '?' + urllib.parse.urlencode({'accessions': ','.join(b), 'fields': 'accession,gene_primary,organism_name,protein_name', 'format': 'tsv'})
            with urllib.request.urlopen(url, timeout=180) as r:
                return {'tsv': r.read().decode('utf-8')}
        ans = cached(cache, f'uniprot/names_{hashlib.md5(",".join(b).encode()).hexdigest()[:12]}.json', fetch)
        lines = ans['tsv'].splitlines()
        for line in lines[1:]:
            acc, gene, org, name = (line.split('\t') + ['', '', '', ''])[:4]
            rows.append({'accession': acc, 'gene': gene.split(';')[0].strip(), 'organism': re.sub(r' \(.*$', '', org), 'protein_name': re.split(r' \(', name)[0]})
    return pd.DataFrame(rows, columns=['accession', 'gene', 'organism', 'protein_name']).drop_duplicates('accession').sort_values('accession')


# ---- Main --------------------------------------------------------------------------------------------------------

DESCR = {
    'pdb_entries': ('PDB entries with an interface or arginine contact at the sites of the release, or modelling a methylarginine (pdb_methylarginine_residues): title, method, resolution, release date, primary citation and non-polymer ligands',
                    RCSB, 'RCSB PDB Data API', 'Burley et al., Nucleic Acids Res 51:D488 (2023)', 'CC0 1.0 (PDB data)'),
    'pdb_entities': ('the polymer entities of those PDB entries: type, description, UniProt accessions they align to, source organisms',
                     RCSB, 'RCSB PDB Data API', 'Burley et al., Nucleic Acids Res 51:D488 (2023)', 'CC0 1.0 (PDB data)'),
    'pdb_site_residues': ('for each entity of those entries aligned to a site\'s protein: the residue modelled at the site (a chemical component code: ARG, or a methylarginine such as 2MR, DA2 or NMM, or another residue where the construct differs), its entity position, whether the 5 residues on each side match the human sequence (matches; realigned: the alignment was shifted, as for an isoform, and the site found by its neighbours; differs), and the other modified residues of the entity (UniProt numbering where it aligns)',
                          RCSB, 'RCSB PDB Data API (entity sequences, UniProt alignments, chemical components)', 'Burley et al., Nucleic Acids Res 51:D488 (2023)', 'CC0 1.0 (PDB data)'),
    'pdb_methylated_sites': (f'every methylarginine (2MR, DA2, NMM) those entries model that is a site of the release: by the entity\'s UniProt alignment, or else by sequence (the longest stretch of the residues within {WINDOW} of it identical to the human sequence, at least {MIN_MATCH} long and of {MIN_DISTINCT} residue types: orthologues, synthetic and unaligned peptides, and the other histone variants sharing a tail), with the entity\'s own UniProt accessions, organisms and chains',
                             RCSB, 'RCSB PDB Data API (entity sequences, UniProt alignments, chemical components); UniProt reviewed human sequences of 2026-10-16', 'Burley et al., Nucleic Acids Res 51:D488 (2023)', 'CC0 1.0 (PDB data)'),
    'pdb_methylated_site_contacts': (f'what each of those methylarginines touches: the residues of other chains within {CONTACT_CUTOFF:g} Å of its side chain (first model, asymmetric unit, at most two copies of the entity), by partner chain and entity; a row with no partner chain is a residue that touches nothing',
                                     MODELSERVER.format(pdb='{pdb}'), 'RCSB ModelServer residueSurroundings', 'Sehnal et al., Nucleic Acids Res 49:W431 (2021)', 'CC0 1.0 (PDB data)'),
    'intact_interactions': ('the binary IntAct records of the interactions scripts/37 places at the sites (binding regions of at most 20 residues, mutations of at most 3, methylated residues): both participants with type and species, their features, detection method, interaction type, host organism, negative flag, MI score and publication',
                            INTACT, 'IntAct web service', 'del Toro et al., Nucleic Acids Res 50:D648 (2022)', 'CC BY 4.0'),
    'pubmed_papers': ('title, first author, journal, year and DOI of the papers of those IntAct records (IntAct gives the PubMed id and authors, not the title), and the arginines their title and abstract name by number (R3, Arg-3, arginine 3, H4R3; only the numbers are kept)',
                      ESUMMARY, 'NCBI E-utilities esummary (PubMed)', 'Sayers et al., Nucleic Acids Res 53:D20 (2025), doi:10.1093/nar/gkae979',
                      'NLM Terms and Conditions: no charge; credit "Courtesy of the U.S. National Library of Medicine"; may not reflect the current PubMed record'),
    'uniprot_partner_names': ('primary gene, organism and recommended name of the UniProt accessions those PDB entities and IntAct participants name that are not reviewed human entries of the release\'s reference (other species, unreviewed entries)',
                              UNIPROT, 'UniProt REST API', 'UniProt Consortium, Nucleic Acids Res 53:D609 (2025)', 'CC BY 4.0'),
    'intact_terms': ('names of the PSI-MI detection methods and interaction types, and the PSI-MOD or PSI-MI term of each IntAct feature label, in those records',
                     OLS, 'EMBL-EBI Ontology Lookup Service (PSI-MI, PSI-MOD)', 'McLaughlin et al., Nucleic Acids Res 53:D1318 (2025); Montecchi-Palazzi et al., Nat Biotechnol 26:864 (2008)', 'CC BY 4.0'),
}


CACHE_OF = {'pdb_entries': 'rcsb', 'pdb_entities': 'rcsb', 'pdb_site_residues': 'rcsb', 'pdb_methylated_sites': 'rcsb',
            'pdb_methylated_site_contacts': 'modelserver_methyl', 'intact_interactions': 'intact_ws', 'intact_terms': 'ols', 'pubmed_papers': 'ncbi',
            'uniprot_partner_names': 'uniprot'}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--release', default='releases/v1.4.2', help='the release whose sites and lower-confidence sites set the scope')
    ap.add_argument('--inputs', default='inputs')
    ap.add_argument('--cache', default='PTM_data/interactions_cache')
    ap.add_argument('--fasta', default='inputs/interactions/uniprot_human_reviewed_2026-10-16.fasta',
                    help='the human reviewed sequences a methylated peptide of the PDB is matched against')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--date', default=dt.date.today().isoformat())
    ap.add_argument('--install', action='store_true')
    args = ap.parse_args()
    inputs, cache, out = Path(args.inputs), Path(args.cache), Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    S = site_positions(Path(args.release))
    print(f'sites: {sum(map(len, S.values())):,} on {len(S):,} proteins ({args.release})', flush=True)
    seqs = {a: q for a, q in read_fasta(Path(args.fasta)).items() if a in S}
    tables, pinfo = stage_pdb(inputs, cache, S, seqs)
    ia, labels, iinfo = stage_intact(inputs, cache, S)
    terms = stage_terms(cache, ia, labels)
    papers = stage_papers(cache, sorted({p for p in ia['pubmed'] if p.isdigit()}, key=int))
    human = set(rd(inputs / 'reference' / 'uniprot_human_reviewed_2026_03_entries.tsv')['Entry'])
    others = {a for v in tables['pdb_entities']['uniprot_accessions'] for a in v.split('; ') if a} \
        | {re.split(r'[-_]', a)[0] for side in 'ab' for a, t in zip(ia[f'participant_{side}'], ia[f'participant_type_{side}']) if t == 'protein'}
    names = stage_names(cache, sorted(a for a in others - human if re.fullmatch(r'[A-Z0-9]{6,10}', a)))
    tables.update({'intact_interactions': ia, 'intact_terms': terms, 'pubmed_papers': papers, 'uniprot_partner_names': names})
    names = {k: f'{k}_{args.date}.tsv' for k in tables}
    for k, t in tables.items():
        # one line per row: titles carry line breaks, tabs and quotes in some answers
        t = t.apply(lambda c: c.str.replace(r'[\t\r\n]+', ' ', regex=True).str.replace('"', "'").str.strip() if c.dtype == object else c)
        tables[k] = t
        t.to_csv(out / names[k], sep='\t', index=False, lineterminator='\n')
    summary = {'release': str(args.release), 'date': args.date, 'sites': sum(map(len, S.values())), 'proteins': len(S), 'pdb': pinfo, 'intact': iinfo,
               'terms': {'rows': len(terms), 'features_not_exact': terms.loc[(terms['kind'] == 'feature') & (terms['match'] != 'exact'), 'label'].tolist()}}
    (out / f'interaction_sources_summary_{args.date}.json').write_text(json.dumps(summary, indent=1, sort_keys=True) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=1)[:3000])
    # cached answers, as raw_cache_files.tsv lists them
    raw = [{'path': str(p.relative_to(cache)), 'bytes': str(p.stat().st_size), 'md5': md5_file(p),
            'saved': dt.date.fromtimestamp(p.stat().st_mtime).isoformat()} for p in sorted(cache.rglob('*.json'))
           if p.relative_to(cache).parts[0] in ('rcsb', 'intact_ws', 'ols', 'modelserver_methyl', 'ncbi', 'uniprot')]
    mcols = (inputs / 'MANIFEST.tsv').read_text().split('\n')[0].split('\t')
    mrows = []
    for k, t in tables.items():
        what, src, ver, cite, lic = DESCR[k]
        p = out / names[k]
        mrows.append({'file': f'interactions/{names[k]}', 'bytes': str(p.stat().st_size), 'md5': md5_file(p), 'rows': str(len(t)), 'redistribution': 'git',
                      'content': what, 'source': src, 'version_or_date': f'{ver}, answers saved {args.date}', 'citation': cite, 'licence': lic,
                      'derived_from': f'{CACHE_OF[k]}/*.json in the interaction cache ({sum(1 for r in raw if r["path"].startswith(CACHE_OF[k] + "/"))} files)',
                      'how_made_or_checked': f'scripts/53_stage_interaction_sources.py --release {args.release}: the records at its sites and lower-confidence sites',
                      'original_name': ''})
    s = out / f'interaction_sources_summary_{args.date}.json'
    mrows.append({'file': f'interactions/{s.name}', 'bytes': str(s.stat().st_size), 'md5': md5_file(s), 'rows': '', 'redistribution': 'git',
                  'content': 'counts of the provenance staging: entries and interactions asked and answered, entities, site residues not modelled as arginine, methylated sites in structures, feature labels without an exact ontology term',
                  'source': 'scripts/53_stage_interaction_sources.py', 'version_or_date': args.date, 'citation': '', 'licence': 'CC0 1.0',
                  'derived_from': '', 'how_made_or_checked': 'written with the tables', 'original_name': ''})
    pd.DataFrame(mrows)[mcols].to_csv(out / 'manifest_rows.tsv', sep='\t', index=False, lineterminator='\n')
    pd.DataFrame(raw).to_csv(out / 'raw_cache_rows.tsv', sep='\t', index=False, lineterminator='\n')
    if args.install:
        dest = inputs / 'interactions'
        for n in list(names.values()) + [s.name]:
            shutil.copy2(out / n, dest / n)
        update_manifest(inputs / 'MANIFEST.tsv', mrows)
        upsert_tsv_rows(dest / 'raw_cache_files.tsv', raw)
        refresh_manifest_row(inputs / 'MANIFEST.tsv', dest / 'raw_cache_files.tsv', 'interactions/raw_cache_files.tsv')
        print(f'installed {len(names) + 1} files into {dest}')


if __name__ == '__main__':
    main()

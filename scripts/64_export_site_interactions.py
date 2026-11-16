#!/usr/bin/env python3
"""Site-level interaction tables for any set of modified residues: what scripts/37 computes at methylarginines, with
the modified residue(s) a parameter.

    python scripts/64_export_site_interactions.py --site-set phospho --out PTM_data/site_interactions/phospho
    python scripts/64_export_site_interactions.py --site-set phospho --out ... --version v1.5.1 --release releases/v1.5.1

A site set (config/site_sets.json) names the residue letters it modifies and the release tables its sites come from;
the same file gives, for each residue letter, the side-chain atoms a contact is measured from and the PDB components
that model its modified forms. Nothing here names a residue.

For the sites of the set it writes, in the columns of the methylarginine release tables:

  site_interactions.tsv         PDBe-KB interfaces at the site (pdb_protein, pdb_rna, pdb_dna, pdb_self) and the published
                                residue-level evidence from cells: RNA crosslinks (pepR-MS and the four studies it
                                re-analyses), CAPRI peptides, RBR-ID residues, DNA crosslinks (Trendel 2025) and UniProt
                                DNA-binding regions
  residue_evidence.tsv          the same evidence along each protein of the set, as runs of residues
  nucleic_acid_contacts.tsv     the RNA and DNA chains within 4 Å of the site's side chain in the PDB entries read (RCSB
                                ModelServer, SIFTS numbering), with the bound sequence around the closest nucleotide
  nucleic_acid_near_sites.tsv   the bound sequence a residue at i ± 1-5 touches where the site's own side chain touches none
  nucleic_acid_proximal_sites.tsv  RNA or DNA within 8 Å of the side chain where none is within 4 Å (the near-site table
                                leaves these sites out, so the tier is written too)

Left out, whatever the set: RBDmap fragments and pCLAP peptides (their extents are not republished), predicted
interfaces and contacts, and IntAct features (scripts/37 computes those for methylarginines only).

The steps are scripts/37's own: its functions are called for every request and measurement (fetch_interfaces,
site_contacts, chain_na_contacts, chain_na_proximity, nucleic_contacts, load_nucleic_chains, write_release_tables), and
the stages between them follow its main() line for line, so the methylarginine set gives the rows of the methylarginine
release (tests/test_site_interactions.py). Data come from the staged tables (inputs/interactions/); what they lack
(proteins never looked up, residues and chains never measured) is fetched into --topup-cache, which
scripts/37a_stage_interactions.py --topup stages as dated tables.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
_spec = importlib.util.spec_from_file_location('export37', HERE / '37_export_interaction_data.py')
x37 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(x37)
import interaction_inputs as ii  # noqa: E402

log = x37.log
base_acc = x37.base_acc
TABLES = ('site_interactions', 'residue_evidence', 'nucleic_acid_contacts', 'nucleic_acid_near_sites', 'nucleic_acid_proximal_sites')
# Evidence scripts/37 computes that this step never writes: extents that may not be republished.
NOT_REPUBLISHED = ('rbdmap_fragment', 'pclap_peptide', 'pclap_lysc_region')
# The staged tables and cache folders of the two per-residue structure steps: scripts/37's own for arginine, one shared
# set for every other residue (a position holds one residue, so the tables of several residues do not collide).
ARGININE_TABLES = {'contacts': 'pdb_arginine_contacts_', 'chains': 'pdb_chain_arginines_', 'distances': 'pdb_arginine_nucleic_acid_distances_',
                   'chain_column': 'arginines'}
RESIDUE_TABLES = {'contacts': 'pdb_residue_contacts_', 'chains': 'pdb_chain_residues_', 'distances': 'pdb_residue_nucleic_acid_distances_',
                  'chain_column': 'residues'}
RESIDUE_FOLDERS = ('residue_contacts', 'na_residues')


# ----------------------------------------------------------------------------- site sets
def load_config(path: Path = ROOT / 'config' / 'site_sets.json') -> dict:
    return json.loads(Path(path).read_text(encoding='utf-8'))


def residues_of(config: dict, letters) -> tuple:
    """The Residues (scripts/37) a set's structure steps read and the staged tables that hold them. A set of arginines
    alone is scripts/37's own ARGININE, with its tables; any other set reads every other residue of the configuration at
    once, so one request of a chain serves all of them, and keeps its own letters when it reads the answers."""
    arginine = set(x37.ARGININE.letters)
    spec_of = lambda ls: {c: config['residues'][r]['side_chain_atoms']  # noqa: E731
                          for r in ls for c in [config['residues'][r]['component'], *config['residues'][r]['modified_components']]}
    if set(letters) <= arginine:
        atoms = spec_of(sorted(arginine))
        if {c: frozenset(a) for c, a in atoms.items()} != x37.ARGININE.atoms:
            raise SystemExit('config/site_sets.json and scripts/37 disagree on the arginine side chain or its modified forms')
        return x37.ARGININE, ARGININE_TABLES
    if set(letters) & arginine:
        raise SystemExit('a site set mixes arginine with other residues: split it')
    others = sorted(set(config['residues']) - arginine)
    return x37.Residues(spec_of(others), *RESIDUE_FOLDERS, ''.join(others)), RESIDUE_TABLES


def read_sites(config: dict, name: str, root: Path = ROOT) -> list:
    """The sites of a set, [(accession, position, residue letter)] in the order of its tables, each residue once."""
    spec = config['site_sets'][name]
    cols, letters = spec['columns'], spec['residues']
    seen, out = set(), []
    for table in spec['tables']:
        with open(root / table, encoding='utf-8') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                letter = r[cols['residue']] if cols.get('residue') else (letters[0] if len(letters) == 1 else '')
                key = (r[cols['accession']], int(r[cols['position']]))
                if letter not in letters:
                    raise SystemExit(f'{table}: {key[0]} {key[1]} is residue {letter!r}, not one of {letters}')
                if key not in seen:
                    seen.add(key)
                    out.append((*key, letter))
    return out


# ----------------------------------------------------------------------------- staged tables, then the network
class Source:
    """What the steps read: the staged tables first, scripts/37's own request for what they lack (into topup)."""

    def __init__(self, folder: Path, topup: Path, residues, tables: dict, workers: int):
        self.t, self.topup, self.residues, self.workers = ii.Tables(folder), topup, residues, workers
        self.fetched = Counter()
        self.queried = ii.queried_interfaces(self.t)
        self.na_chains, self.na_looked = ii.nucleic_chains(self.t)
        self.na_res = ii.na_residue_contacts(self.t)
        self.contacts = self._contacts(tables)
        self.proximity = self._proximity(tables)

    def _rows(self, prefix: str) -> list:
        try:
            return list(self.t.rows(prefix))
        except KeyError:
            return []

    def _contacts(self, tables: dict) -> dict:
        out = {}
        for r in self._rows(tables['contacts']):
            rec = out.setdefault((r['pdb_id'], r['chain'], int(r['seqres_index'])), {'resn': r['residue'], 'hits': []})
            if r['contact_chain']:
                rec['hits'].append((r['contact_chain'], r['contact_seqres_index'], r['contact_author_number'], r['contact_residue'], float(r['distance_a'])))
        return out

    def _proximity(self, tables: dict) -> dict:
        out = {}
        for r in self._rows(tables['chains']):
            out[(r['pdb_id'], r['chain'])] = {'args': [[int(x.split(':')[0]), x.split(':')[1], x.split(':')[2]]
                                                       for x in r[tables['chain_column']].split('; ') if x], 'hits': []}
        for r in self._rows(tables['distances']):
            out[(r['pdb_id'], r['chain'])]['hits'].append(
                [int(r['seqres_index']), r['na_chain'], int(r['closest_na_seqres_index']), r['closest_nucleotide'],
                 float(r['side_chain_distance_a']), [int(x) for x in r['na_seqres_within_8a'].split('; ') if x]])
        return out

    def interfaces(self, accs: list) -> dict:
        out = ii.interfaces(self.t, [a for a in accs if a in self.queried])
        missing = [a for a in accs if a not in self.queried]
        if missing:
            log(f'PDBe interfaces: {len(missing):,} proteins not in the staged tables, read from {self.topup}')
            self.fetched['pdbe_interfaces'] += len(missing)
            out.update(x37.fetch_interfaces(self.topup, missing, False, self.workers))
        return out

    def sifts_chains(self, pdbs: set) -> dict:
        return ii.sifts_chains(self.t, pdbs)

    def nucleic_chains(self, pdbs: set) -> dict:
        out = {k: v for k, v in self.na_chains.items() if k[0] in pdbs}
        missing = set(pdbs) - self.na_looked
        if missing:
            log(f'Bound sequences: {len(missing):,} PDB entries not in the staged tables, read from {self.topup}')
            self.fetched['pdbe_molecules'] += len(missing)
            out.update(x37.load_nucleic_chains(self.topup, missing, False, self.workers))
        return out

    def site_contacts(self, pdb: str, chain: str, label: int):
        hit = self.contacts.get((pdb, chain, int(label)))
        if hit:
            return {'resn': hit['resn'], 'hits': [list(h) for h in hit['hits']]}
        self.fetched[self.residues.contacts_folder] += 1
        return x37.site_contacts(self.topup, pdb, chain, label, False, self.residues)

    def chain_na_contacts(self, pdb: str, chain: str, na_chains: set):
        hit = self.na_res.get((pdb, chain))
        if hit is not None:
            return {'hits': [list(h) for h in hit]}
        self.fetched['na_regions'] += 1
        return x37.chain_na_contacts(self.topup, pdb, chain, na_chains, False)

    def chain_na_proximity(self, pdb: str, chain: str, na_chains: set):
        hit = self.proximity.get((pdb, chain))
        if hit is not None:
            return hit
        self.fetched[self.residues.proximity_folder] += 1
        return x37.chain_na_proximity(self.topup, pdb, chain, na_chains, False, self.residues)


def reference_proteome(src: Source, fasta: Path | None) -> dict:
    """accession -> (gene, STRING id, sequence) of the staged UniProt download, kept to the release reference."""
    proteome, _ = ii.uniprot(src.t)
    if fasta:
        ref = x37.read_fasta(Path(fasta))
        differ = sorted(a for a, s in ref.items() if a in proteome and proteome[a][2] != s)
        missing = sorted(a for a in ref if a not in proteome)
        if differ or missing:
            raise SystemExit(f'{fasta} is not the UniProt release of the staged download: {len(differ)} sequences differ, {len(missing)} entries missing')
        proteome = {a: v for a, v in proteome.items() if a in ref}
    return proteome


# ----------------------------------------------------------------------------- the stages of scripts/37 main()
def compute(site_list: list, src: Source, proteome: dict, noun: str = 'site', contact_entries: int = 2, na_contact_entries: int = 0,
            fetch_only: bool = False) -> dict:
    """What scripts/37 write_release_tables needs for the tables above, for site_list: [(accession, position, letter)].
    Sites whose residue is not the letter in the staged sequence (another UniProt release) are left out and counted."""
    seqs = {a: v[2] for a, v in proteome.items()}
    letters_all = set(src.residues.letters)
    kept, dropped = [], 0
    for acc, pos, letter in site_list:
        s = seqs.get(base_acc(acc))
        if '-' not in acc and s is not None and (pos > len(s) or s[pos - 1] != letter):
            dropped += 1
            continue
        kept.append((acc, pos, letter))
    if dropped:
        log(f'{dropped:,} {noun}s left out: the residue differs in the staged UniProt sequences')
    index_of, prots, sites, site_residues = {}, [], [], []
    for acc, pos, letter in kept:
        if acc not in index_of:
            index_of[acc] = len(prots)
            prots.append([acc])
        sites.append([index_of[acc], pos])
        site_residues.append(letter)
    meth_accs = [p[0] for p in prots]
    canonical = [a for a in meth_accs if '-' not in a]
    workers = src.workers

    vocab = list(meth_accs)
    index = {base_acc(a): i for i, a in reversed(list(enumerate(vocab)))}

    def vid(acc):
        acc = base_acc(acc)
        if acc not in index:
            index[acc] = len(vocab)
            vocab.append(acc)
        return index[acc]

    gene_of = {acc: g for acc, (g, _, _) in proteome.items()}
    ifaces = src.interfaces(canonical)
    by_prot_sites = defaultdict(list)
    for j, s in enumerate(sites):
        by_prot_sites[s[0]].append(j)

    # Interfaces at the sites (scripts/37 main: "Interfaces at methylarginine sites")
    iface_detail, iface_residues, contact_jobs = [], {}, []
    for i, acc in enumerate(meth_accs):
        data = ifaces.get(acc)
        if not data:
            continue
        res_partners = defaultdict(dict)
        for part in data:
            pacc = part.get('accession') or ''
            if pacc in ('RNA', 'DNA'):
                key = pacc
            elif x37.ACC_RE.match(base_acc(pacc)):
                key = base_acc(pacc)
            else:
                continue
            for r in part.get('residues', []):
                if r.get('indexType', 'UNIPROT') != 'UNIPROT':
                    continue
                pdbs = {e.get('pdbId') for e in r.get('interactingPDBEntries') or [] if e.get('pdbId')}
                for pos in range(int(r['startIndex']), int(r['endIndex']) + 1):
                    res_partners[pos].setdefault(key, set()).update(pdbs)
        site_pos = {sites[j][1]: j for j in by_prot_sites[i]}
        iface_residues[acc] = res_partners
        for pos, j in site_pos.items():
            hit = res_partners.get(pos)
            if not hit:
                continue
            for key, ids in hit.items():
                iface_detail.append((j, key, sorted(ids)))
                cap = (na_contact_entries or None) if key in ('RNA', 'DNA') else contact_entries
                for pdb in sorted(ids, reverse=True)[:cap]:
                    contact_jobs.append((j, key, pdb, acc, pos))
    log(f'PDB interfaces: {len({j for j, _, _ in iface_detail}):,} of {len(sites):,} {noun}s at an interface of {len(ifaces):,} proteins')

    # Residue-level RNA and DNA evidence from cells (scripts/37 main: "Residue-level RNA and DNA crosslinks"), without the
    # RBDmap fragments and the pCLAP peptides
    rna_xl, _, _, dna_xl, dna_peps, dna_int, dna_bind, rna_single, _ = ii.crosslinks(src.t, seqs, proteome, x37.RNA_XL_SETS, x37.RNA_FOCUSED, x37.D_XL, x37.D_CAND)
    xl_seq_all = ii.crosslink_sequences(src.t, seqs, proteome, x37.RNA_XL_SETS)
    xl_seq_all = {'rna': xl_seq_all['rna'], 'dna': xl_seq_all['dna'], 'rbd': {}}
    capri, _ = ii.capri(src.t, seqs, proteome)
    rbrid, _ = x37.load_rbrid(src.t.folder, seqs)
    xl_rows, xl_prot = [], []
    for i, acc in enumerate(meth_accs):
        if '-' in acc:
            continue
        r_map, d_map, peps, binds = rna_xl.get(acc, {}), dna_xl.get(acc, {}), dna_peps.get(acc, []), dna_bind.get(acc, [])
        c_xl, c_adj = capri['xl'].get(acc, []), capri['adj'].get(acc, [])
        b_pos, b_spans, b_seg, b_cov = (rbrid[k].get(acc, set() if k == 'pos' else []) for k in ('pos', 'spans', 'seg', 'cov'))
        r_single = rna_single.get(acc, set())
        for j in by_prot_sites[i]:
            pos = sites[j][1]
            rb, db = r_map.get(pos, 0), d_map.get(pos, 0)
            nb = 0
            for q in range(pos - x37.XL_NEAR, pos + x37.XL_NEAR + 1):
                nb |= r_map.get(q, 0)
            db |= (x37.D_PEP if any(a <= pos <= b for a, b in peps) else 0) | (x37.D_BIND if any(a <= pos <= b for a, b in binds) else 0)
            cb = (x37.CAPRI_XL if any(a <= pos <= b for a, b, _ in c_xl) else 0) | (x37.CAPRI_ADJ if any(a <= pos <= b for a, b, _ in c_adj) else 0) | \
                (x37.RBRID_RES if pos in b_pos else 0)
            if rb or nb or db or cb:
                xl_rows.append([j, rb, nb, db, 0, 1 if pos in r_single else 0, cb])
        if r_map or d_map or peps or binds or acc in dna_int or c_xl or c_adj or b_spans or b_seg:
            xl_prot.append([i, sorted(r_map), [r_map[q] for q in sorted(r_map)], sorted(q for q, v in d_map.items() if v & x37.D_XL), sorted(set(peps)),
                            sorted(set(binds)), 1 if acc in dna_int else 0, [], sorted(q for q, v in d_map.items() if v == x37.D_CAND),
                            sorted({(a, b) for a, b, _ in c_xl}), sorted({(a, b) for a, b, _ in c_adj}), [], [], [],
                            [list(x) for x in b_spans], [list(x) for x in b_seg], sorted(b_cov)])
    log(f'Cells: {sum(1 for r in xl_rows if r[1]):,} {noun}s are RNA-crosslinked residues, {sum(1 for r in xl_rows if r[3] & (x37.D_XL | x37.D_CAND)):,} DNA-crosslinked or candidates, '
        f'{sum(1 for r in xl_rows if r[6] & (x37.CAPRI_XL | x37.CAPRI_ADJ)):,} in a CAPRI peptide, {sum(1 for r in xl_rows if r[6] & x37.RBRID_RES):,} under an RBR-ID+ peptide')

    # Every entry in which a protein of the set touches RNA or DNA (PDBe-KB)
    na_entries = defaultdict(set)
    for acc in meth_accs:
        for part in ifaces.get(acc) or []:
            if part.get('accession') in ('RNA', 'DNA'):
                for r in part.get('residues', []):
                    na_entries[acc].update(e.get('pdbId') for e in r.get('interactingPDBEntries') or [] if e.get('pdbId'))

    # Residues each interface site touches, from the structures themselves (scripts/37 main: "Structure contacts")
    na_rows, chains_na = [], {}
    if contact_jobs:
        chains = src.sifts_chains({pdb for _, _, pdb, _, _ in contact_jobs})
        requests_ = defaultdict(list)
        for j, key, pdb, acc, pos in contact_jobs:
            cands = sorted({(ch, res_beg + (pos - sp_beg)) for ch, sp, res_beg, sp_beg, sp_end in chains.get(pdb, [])
                            if base_acc(sp) == acc and sp_beg <= pos <= sp_end})
            for ch, auth in cands[:2]:
                requests_[(pdb, ch, auth)].append((j, key))
        todo = sorted(requests_)
        log(f'Structure contacts: {len(todo):,} residue neighbourhoods to read')
        results, done = {}, 0
        with ThreadPoolExecutor(max_workers=workers) as pool:
            for req, res in pool.map(lambda req: (req, src.site_contacts(*req)), todo):
                results[req] = res
                done += 1
                if done % 5000 == 0:
                    log(f'  {done:,}/{len(todo):,}')
        na_pdbs = {pdb for _, key, pdb, _, _ in contact_jobs if key in ('RNA', 'DNA')} | set().union(*na_entries.values())
        na_pdbs |= {pdb for (pdb, _, _), res in results.items() if res and any(h[3] in x37.NUC_RNA or h[3] in x37.NUC_DNA for h in res['hits'])}
        chains_na = src.nucleic_chains(na_pdbs)
        log(f'Bound sequences: {len(chains_na):,} RNA and DNA chains in {len(na_pdbs):,} PDB entries')
        na_rows = x37.nucleic_contacts(results, requests_, chains_na)
        log(f"Bound sequences: {len({r[0] for r in na_rows if r[3] != 'dna'}):,} {noun}s in contact with RNA and "
            f"{len({r[0] for r in na_rows if r[3] == 'dna'}):,} with DNA in the structures read ({len(na_rows):,} site-chain pairs)")

    # Regions, the proximal tier and the near-site tier (scripts/37 main: "Regions", "Proximal", "Bound sequences via a
    # flanking residue")
    near_rows, prox_rows = [], []
    if chains_na:
        na_by_pdb = defaultdict(list)
        for (pdb, ch), v in chains_na.items():
            na_by_pdb[pdb].append(ch)
        complexes = defaultdict(lambda: defaultdict(list))
        for acc, ents in na_entries.items():
            for pdb in ents:
                bound = tuple(sorted({(chains_na[(pdb, ch)][0], chains_na[(pdb, ch)][4]) for ch in na_by_pdb.get(pdb, [])}))
                if bound:
                    complexes[acc][bound].append(pdb)
        reps = {max(ents) for by_bound in complexes.values() for ents in by_bound.values()}
        rchains = src.sifts_chains(reps)
        seg_of = defaultdict(list)
        for pdb, rows in rchains.items():
            for ch, sp, beg, sp_beg, sp_end in rows:
                seg_of[(pdb, ch)].append((beg, sp_beg, sp_end, base_acc(sp)))
        sig_entries, jobs = [], []
        for acc, by_bound in sorted(complexes.items()):
            for _, ents in sorted(by_bound.items()):
                pdb = max(ents)
                for ch in sorted({ch for ch, sp, *_ in rchains.get(pdb, []) if base_acc(sp) == acc})[:2]:
                    jobs.append((acc, len(sig_entries), pdb, ch))
                sig_entries.append(sorted(ents, reverse=True))
        log(f'Nucleic-acid regions: {len(jobs):,} protein chains in {len(sig_entries):,} distinct complexes of {len(complexes):,} proteins')
        per_res = defaultdict(lambda: defaultdict(set))
        with ThreadPoolExecutor(max_workers=workers) as pool:
            for (acc, si, pdb, ch), res in pool.map(lambda job: (job, src.chain_na_contacts(job[2], job[3], set(na_by_pdb.get(job[2], [])))), jobs):
                for lab, _, _, nch, nlab, _, _ in res['hits']:
                    unp = next((sp_beg + (lab - beg) for beg, sp_beg, sp_end, sa in seg_of[(pdb, ch)]
                                if sa == acc and sp_beg <= sp_beg + (lab - beg) <= sp_end), None)
                    if unp is not None and (pdb, nch) in chains_na:
                        per_res[(acc, 'dna' if chains_na[(pdb, nch)][0] == 'dna' else 'rna')][unp].add((si, pdb, nch, nlab))
        site_at = {(prots[sites[j][0]][0], sites[j][1]): j for j in range(len(sites))}
        at_entry = {(r[0], r[1]) for r in na_rows}
        sites_of = defaultdict(list)
        for (acc, q), j in site_at.items():
            sites_of[acc].append((q, j))
        cx_of, si = {}, 0
        for acc, by_bound in sorted(complexes.items()):
            for _, ents in sorted(by_bound.items()):
                for e in ents:
                    cx_of[(acc, e)] = si
                si += 1
        achains = src.sifts_chains({e for _, e in cx_of})
        aseg = defaultdict(list)
        for pdb, rows in achains.items():
            for ch, sp, beg, sp_beg, sp_end in rows:
                aseg[(pdb, ch)].append((beg, sp_beg, sp_end, base_acc(sp)))
        arg_jobs = [(acc, si, pdb, ch) for (acc, pdb), si in sorted(cx_of.items())
                    for ch in sorted({ch for ch, sp, *_ in achains.get(pdb, []) if base_acc(sp) == acc})[:2]]
        log(f'Side-chain distances: {len(arg_jobs):,} protein chains in {len({p for _, p in cx_of}):,} entries read within {x37.PROXIMAL_CUTOFF:g} Å of RNA or DNA')
        site_cx, done = {}, 0
        with ThreadPoolExecutor(max_workers=workers) as pool:
            for (acc, si, pdb, ch), res in pool.map(lambda job: (job, src.chain_na_proximity(job[2], job[3], set(na_by_pdb.get(job[2], [])))), arg_jobs):
                done += 1
                if done % 5000 == 0:
                    log(f'  {done:,}/{len(arg_jobs):,}')
                ref = proteome.get(acc, ('', '', ''))[2]
                for lab, nch, k, _, dist, within in res['hits']:
                    q = next((sp_beg + (lab - beg) for beg, sp_beg, sp_end, sa in aseg[(pdb, ch)]
                              if sa == acc and sp_beg <= sp_beg + (lab - beg) <= sp_end), None)
                    if not q or not 0 < q <= len(ref) or ref[q - 1] not in letters_all or (pdb, nch) not in chains_na or not 1 <= k <= len(chains_na[(pdb, nch)][4]):
                        continue
                    kind = 'dna' if chains_na[(pdb, nch)][0] == 'dna' else 'rna'
                    if (acc, q) in site_at:
                        g = site_cx.setdefault((acc, q, si, kind), [99, None, None, None, None, set()])
                        g[5].add(pdb)
                        if dist < g[0]:
                            g[:5] = [dist, pdb, nch, k, within]
        if fetch_only:
            return {}
        for (acc, q0, si, kind), (d, pdb, nch, centre, within, ents) in sorted(site_cx.items()):
            j = site_at[(acc, q0)]
            if any((j, e) in at_entry for e in sig_entries[si]):
                continue
            seq, lo = chains_na[(pdb, nch)][4], centre - x37.NA_WINDOW
            window = ''.join(seq[k - 1] if 1 <= k <= len(seq) else '-' for k in range(lo, centre + x37.NA_WINDOW + 1))
            prox_rows.append([j, kind, pdb, nch, window, sorted(k - lo for k in within if lo <= k <= centre + x37.NA_WINDOW), round(d, 2), sorted(ents, reverse=True), centre])
        prox_cx = {(r[0], cx_of[(prots[sites[r[0]][0]][0], r[2])]) for r in prox_rows}
        for (acc, kind), res in sorted(per_res.items()):
            by_cx = defaultdict(lambda: defaultdict(set))
            for q, hits in res.items():
                for si, pdb, nch, nlab in hits:
                    by_cx[si][q].add((pdb, nch, nlab))
            for q0, j in sorted(sites_of.get(acc, [])):
                for si, at_pos in sorted(by_cx.items()):
                    near = {q: h for q, h in at_pos.items() if 0 < abs(q - q0) <= x37.NEAR_RESIDUES}
                    if not near or (j, si) in prox_cx or any((j, e) in at_entry for e in sig_entries[si]):
                        continue
                    (pdb, nch), _ = min(Counter((pd, ch) for h in near.values() for pd, ch, _ in h).items(), key=lambda kv: (-kv[1], kv[0]))
                    on = {q: sorted(n for pd, ch, n in h if (pd, ch) == (pdb, nch)) for q, h in near.items()}
                    via = min((q for q in on if on[q]), key=lambda q: (abs(q - q0), q))
                    centre = on[via][len(on[via]) // 2]
                    seq, lo = chains_na[(pdb, nch)][4], centre - x37.NA_WINDOW
                    window = ''.join(seq[k - 1] if 1 <= k <= len(seq) else '-' for k in range(lo, centre + x37.NA_WINDOW + 1))
                    offsets = sorted({n - lo for ns in on.values() for n in ns if lo <= n <= centre + x37.NA_WINDOW})
                    near_rows.append([j, kind, pdb, nch, window, offsets, via - q0, sig_entries[si], centre])
        log(f"Within {x37.PROXIMAL_CUTOFF:g} Å and none within {x37.CONTACT_CUTOFF:g} Å: {len({r[0] for r in prox_rows}):,} {noun}s; "
            f"through a residue within {x37.NEAR_RESIDUES}: {len({r[0] for r in near_rows}):,} {noun}s")
    if fetch_only:
        return {}
    empty = {'pep': {}, 'lysc': {}, 'measured': {}}
    return {'prots': prots, 'sites': sites, 'site_residues': site_residues, 'vocab': vocab, 'gene_of': gene_of, 'iface_detail': iface_detail,
            'xl_rows': xl_rows, 'pred_all': [], 'gallery': [], 'intact_rows': [], 'partners': {}, 'domain_classes': {'tudor': set(), 'wd40': set(), 'phd': set()},
            'meth_base': {}, 'iface_residues': iface_residues, 'pred_iface': {}, 'xl_prot': xl_prot, 'intact_prot': [], 'tests': [], 'set_tests': [],
            'complexes': [], 'go': {}, 'xl_seq': xl_seq_all, 'capri': capri, 'pclap': empty, 'rbrid': rbrid, 'proteome': proteome, 'na_rows': na_rows,
            'chains_na': chains_na, 'na_regions': [], 'nab': {}, 'near_rows': near_rows, 'prox_rows': prox_rows, 'arg_rows': []}


# ----------------------------------------------------------------------------- writing
def release_columns() -> dict:
    import release_codes as rc
    return {'site_interactions': rc.SITE_INTERACTION_COLUMNS, 'residue_evidence': rc.RESIDUE_EVIDENCE_COLUMNS,
            'nucleic_acid_contacts': rc.NUCLEIC_ACID_CONTACT_COLUMNS, 'nucleic_acid_near_sites': rc.NUCLEIC_ACID_NEAR_COLUMNS,
            'nucleic_acid_proximal_sites': rc.NUCLEIC_ACID_PROXIMAL_COLUMNS}


def release_text(path: Path, columns: list) -> str:
    """A table scripts/37 wrote, as scripts/build_release.py write_table writes it into a release: the declared
    columns in order, numbers in up to 6 significant digits, integers without decimals."""
    out = ['\t'.join(c for c, *_ in columns)]
    with open(path, encoding='utf-8', newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            cells = []
            for name, kind, *_ in columns:
                v = r[name]
                cells.append(v if v == '' or kind not in ('integer', 'number') else str(int(v)) if kind == 'integer' else f'{float(v):.6g}')
            out.append('\t'.join(cells))
    return '\n'.join(out) + '\n'


def write_tables(d: dict, outdir: Path, prefix: str = '', suffix: str = '') -> dict:
    """The tables, as <prefix><table><suffix>.tsv in outdir, in the release's format. Returns {table: path}."""
    outdir.mkdir(parents=True, exist_ok=True)
    raw = outdir / '.raw'
    x37.write_release_tables(raw, d, only={f'{t}.tsv' for t in TABLES})
    out = {}
    for name, columns in release_columns().items():
        text = release_text(raw / f'{name}.tsv', columns)
        if any(f'\t{e}\t' in text for e in NOT_REPUBLISHED):
            raise SystemExit(f'{name}: evidence that may not be republished')
        out[name] = outdir / f'{prefix}{name}{suffix}.tsv'
        out[name].write_text(text, encoding='utf-8', newline='')
        (raw / f'{name}.tsv').unlink()
    raw.rmdir()
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--site-set', required=True, help='a site set of --config')
    ap.add_argument('--config', default=str(ROOT / 'config' / 'site_sets.json'))
    ap.add_argument('--tables', default=str(ROOT / 'inputs' / 'interactions'), help='the staged tables')
    ap.add_argument('--topup-cache', default='PTM_data/site_interactions_topup', help='where to save what the staged tables do not cover')
    ap.add_argument('--reference-fasta', default=str(ROOT / 'inputs' / 'reference' / 'uniprot_human_reviewed_2026_03.fasta'))
    ap.add_argument('--out', required=True, help='folder for the tables')
    ap.add_argument('--version', default='', help='name the tables <site set>_<table>_<version>.tsv (a release); default <table>.tsv')
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--contact-entries', type=int, default=2, help='PDB entries examined per site and protein partner (as scripts/37)')
    ap.add_argument('--na-contact-entries', type=int, default=0, help='PDB entries examined per site with RNA or DNA (0: every entry, as scripts/37)')
    ap.add_argument('--fetch-only', action='store_true', help='only fill --topup-cache')
    args = ap.parse_args()
    config = load_config(Path(args.config))
    if args.site_set not in config['site_sets']:
        raise SystemExit(f'--site-set: one of {", ".join(config["site_sets"])}')
    spec = config['site_sets'][args.site_set]
    residues, tables = residues_of(config, spec['residues'])
    src = Source(Path(args.tables), Path(args.topup_cache), residues, tables, args.workers)
    proteome = reference_proteome(src, Path(args.reference_fasta) if args.reference_fasta else None)
    site_list = read_sites(config, args.site_set)
    log(f'{args.site_set}: {len(site_list):,} {spec["noun"]}s on {len({s[0] for s in site_list}):,} proteins; residues {", ".join(spec["residues"])}')
    d = compute(site_list, src, proteome, spec['noun'], args.contact_entries, args.na_contact_entries, args.fetch_only)
    if src.fetched:
        log('Read from the network into ' + args.topup_cache + ': ' + ', '.join(f'{n:,} {k}' for k, n in sorted(src.fetched.items()))
            + '. Stage them (scripts/37a_stage_interactions.py --topup) so the tables can be built again without it.')
    if args.fetch_only:
        return
    out = write_tables(d, Path(args.out), f'{args.site_set}_' if args.version else '', f'_{args.version}' if args.version else '')
    for name, path in out.items():
        log(f'{path}: {sum(1 for _ in open(path, encoding="utf-8")) - 1:,} rows')


if __name__ == '__main__':
    main()

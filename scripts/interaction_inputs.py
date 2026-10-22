"""The frozen interaction tables (inputs/interactions/) in the shapes scripts/37's loaders return.

scripts/37_export_interaction_data.py was written to read a cache of raw downloads. With --tables it reads these
tables instead (scripts/37b_check_staged_interactions.py shows that each table rebuilds what its loader returns from
the cache). Proteins, sites and structures the 2026-10-16 download did not cover (proteins and sites new in a release)
are fetched by scripts/37's own code into a separate top-up cache, and scripts/37a_stage_interactions.py --topup writes
what was fetched as dated tables next to the frozen ones.
"""
from __future__ import annotations

import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

csv.field_size_limit(1 << 30)


def L(value: str) -> list[str]:
    return value.split('; ') if value else []


class Tables:
    """Rows of every table of a kind in a folder: the first staging and every dated top-up, oldest first."""

    def __init__(self, folder: Path):
        self.folder = Path(folder)
        self.names = sorted(p.name for p in self.folder.glob('*.tsv'))

    def paths(self, prefix: str) -> list[Path]:
        hits = [self.folder / n for n in self.names if n.startswith(prefix)]
        if not hits:
            raise KeyError(f'no table {prefix}* in {self.folder}')
        return hits

    def rows(self, prefix: str):
        for path in self.paths(prefix):
            with open(path, encoding='utf-8') as fh:
                yield from csv.DictReader(fh, delimiter='\t')


def expand(ranges: str) -> set[int]:
    out = set()
    for part in filter(None, ranges.split(',')):
        a, _, b = part.partition('-')
        out.update(range(int(a), int(b or a) + 1))
    return out


# ---- Bulk sources ------------------------------------------------------------------------------------------------------

def uniprot(t: Tables):
    seqs, acc = {}, None
    with open(t.paths('uniprot_human_reviewed_')[0].with_suffix('.fasta'), encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('>'):
                acc = line[1:].split()[0]
                seqs[acc] = []
            else:
                seqs[acc].append(line.strip())
    seqs = {a: ''.join(v) for a, v in seqs.items()}
    proteome, classes = {}, {'tudor': set(), 'wd40': set(), 'phd': set()}
    for r in t.rows('uniprot_human_reviewed_'):
        a = r['uniprot_accession']
        proteome[a] = ('; '.join(L(r['gene_names_primary'])), r['string_id'], seqs[a])
        for key, col in (('tudor', 'tudor_domain'), ('wd40', 'wd40_repeat'), ('phd', 'phd_finger')):
            if r[col] == 'true':
                classes[key].add(a)
    return proteome, classes


def regions(t: Tables):
    dis, dom = defaultdict(list), defaultdict(list)
    for r in t.rows('uniprot_human_regions_'):
        (dis if r['feature'] == 'disordered_region' else dom)[r['uniprot_accession']].append((int(r['start']), int(r['end'])))
    return dis, dom


def go(t: Tables):
    out = defaultdict(lambda: {'rna': set(), 'dna': set()})
    for r in t.rows('quickgo_'):
        out[r['uniprot_accession']][r['binding']].add(r['evidence_code'])
    return out


def string(t: Tables, keep: set, min_score: int, sid_to_acc: dict):
    """Physical links (a, b) -> score, in the order scripts/37 meets them in the STRING file (sorted by STRING id), so
    proteins outside the union are numbered as in a run from the raw download."""
    sid = {acc: s for s, acc in sid_to_acc.items()}
    edges = {}
    for r in t.rows('string_v'):
        a, b, sc = r['protein_a'], r['protein_b'], int(r['combined_score'])
        if sc >= min_score and (a in keep or b in keep):
            key = (a, b) if a < b else (b, a)
            edges[key] = max(sc, edges.get(key, 0))
    return dict(sorted(edges.items(), key=lambda kv: tuple(sorted((sid.get(kv[0][0], ''), sid.get(kv[0][1], ''))))))


def complexes(t: Tables):
    out = {}
    for r in t.rows('complex_portal_'):
        key = (r['complex_id'], r['curation'])
        d = out.setdefault(key, {'id': r['complex_id'], 'name': r['complex_name'], 'predicted': r['curation'] == 'predicted', 'slots': []})
        k = int(r['member_slot']) - 1
        while len(d['slots']) <= k:
            d['slots'].append([])
        d['slots'][k].append(r['uniprot_accession'])
    return list(out.values())


def queried_interfaces(t: Tables) -> set:
    return {r['uniprot_accession'] for r in t.rows('pdbe_kb_interface_queries_')}


def interfaces(t: Tables, accs) -> dict:
    """PDBe-KB answers as the API gives them: acc -> [{accession, name, residues: [{startIndex, endIndex,
    interactingPDBEntries}]}]; one residue range per table row (the API gives one per residue; scripts/37 expands them)."""
    want = set(accs)
    parts = defaultdict(dict)
    for r in t.rows('pdbe_kb_interface_residues_'):
        a = r['uniprot_accession']
        if a not in want:
            continue
        part = parts[a].setdefault((r['partner'], r['partner_name']), {'accession': r['partner'], 'name': r['partner_name'], 'residues': []})
        part['residues'].append({'startIndex': int(r['start']), 'endIndex': int(r['end']), 'indexType': 'UNIPROT',
                                 'interactingPDBEntries': [{'pdbId': p} for p in L(r['pdb_entries'])]})
    return {a: list(v.values()) for a, v in parts.items()}


def observed(t: Tables, keep: set):
    obs = defaultdict(list)
    for r in t.rows('sifts_observed_ranges_'):
        if r['uniprot_accession'] in keep:
            obs[r['uniprot_accession']].append((int(r['start']), int(r['end'])))
    return obs


def sifts_chains(t: Tables, keep_pdbs: set):
    """pdb -> [(chain, accession, SEQRES start, UniProt start, UniProt end)] in the order of the SIFTS file (chain, then
    SEQRES start): scripts/37 takes the last segment of a fusion chain as the chain's protein."""
    out = defaultdict(list)
    for r in t.rows('sifts_pdb_chain_uniprot_'):
        if r['pdb_id'] in keep_pdbs:
            out[r['pdb_id']].append((r['chain'], r['uniprot_accession'], int(r['seqres_start']), int(r['uniprot_start']), int(r['uniprot_end'])))
    for rows in out.values():
        rows.sort(key=lambda x: (x[0], x[2]))
    return out


def modelled_methyl_r(t: Tables):
    return [(r['pdb_id'], r['chain'], int(r['seqres_index']), r['chem_comp']) for r in t.rows('pdb_methylarginine_')]


def labels_agree(labels, seq: str) -> bool:
    """Residue labels such as R83 agree with the sequence at their numbers (all of them, and at least one)."""
    for x in labels:
        aa, num = x[:1], x[1:]
        if not num.isdigit() or not 1 <= int(num) <= len(seq) or seq[int(num) - 1] != aa:
            return False
    return bool(labels)


def protvar(t: Tables, keep: set, floor: float, seqs: dict | None = None):
    """Predicted interface residues of each protein (ProtVar 2024, pDockQ >= floor). With seqs (the reviewed proteome),
    a pair with a partner outside it (another species, an unreviewed entry) is left out, and so is a side whose residue
    labels disagree with the protein's reference sequence (numbered on another sequence or isoform)."""
    out = defaultdict(list)
    for r in t.rows('protvar_predicted_'):
        a, b = r['uniprot_id1'].split('-')[0], r['uniprot_id2'].split('-')[0]
        score = float(r['pdockq'])
        if a == b or score < floor or (a not in keep and b not in keep):
            continue
        if seqs is not None and (a not in seqs or b not in seqs):
            continue
        for acc, other, field, first in ((a, b, 'ifresid1', True), (b, a, 'ifresid2', False)):
            if acc not in keep:
                continue
            res = [x for x in r[field].split(',') if x]
            if seqs is not None and not labels_agree(res, seqs[acc]):
                continue
            pos_all = {int(x[1:]) for x in res if x[1:].isdigit()}
            pos_r = {int(x[1:]) for x in res if x[:1] == 'R' and x[1:].isdigit()}
            if pos_all:
                out[acc].append((other, score, pos_r, pos_all, first))
    return out


def protvar_served_rows(t: Tables) -> list[dict]:
    return list(t.rows('protvar_served_'))


def intact(t: Tables, effect_mi: dict, region_mi: dict):
    """IntAct features scripts/37 uses, from the tables; counts of the features kept and left out (the tables hold
    reviewed human proteins only, so features on isoforms and other entries are not counted here)."""
    feats, counts = defaultdict(list), Counter()
    part_re = re.compile(r'([^:]+):(.*) \(([^()]*)\)$')
    for kind, label, name in ((0, 'mutations', 'mutations'), (1, 'binding_regions', 'bindings_regions'), (2, 'methylarginines', 'ptms')):
        for r in t.rows(f'intact_{label}_human_'):
            check = r['sequence_check']
            if check not in ('matches', 'no curated sequence'):
                counts[f'{name}:' + {'mismatch': 'sequence_mismatch', 'undetermined range': 'undetermined_range'}.get(check, check.replace(' ', '_'))] += 1
                continue
            a, mi = r['uniprot_accession'], r['feature_type_mi']
            code = 0 if kind == 2 else effect_mi.get(mi, 5) if kind == 0 else region_mi.get(mi)
            spans = [tuple(map(int, x.split('-'))) for x in r['ranges'].split(',')]
            parts = [part_re.match(p).groups() for p in L(r['interaction_participants'])]
            own = [p for p in parts if p[0] == 'uniprotkb' and p[1].split('-')[0] == a]
            others = [p for p in parts if not (p[0] == 'uniprotkb' and p[1].split('-')[0] == a)]
            partners = tuple(others) if others else (('uniprotkb', a, 'protein'),) if len(own) > 1 or len(parts) == 1 else ()
            pm = L(r['pubmed'])
            resulting = r['resulting_sequence'] if kind == 0 and r['resulting_sequence'] not in ('', '-') else ''
            for k, (s0, s1) in enumerate(spans):
                res = resulting.split(',')[k] if resulting and len(resulting.split(',')) == len(spans) else resulting
                feats[a].append((kind, code, s0, s1, res, partners, int(pm[0]) if pm else 0, r['interaction_ac']))
            counts[f'{name}:kept'] += 1
    return feats, dict(counts)


def crosslinks(t: Tables, seq: dict, proteome: dict, rna_sets: list, focused_bit: int, d_xl: int, d_cand: int):
    """The 9 values scripts/37 load_crosslinks returns (as rebuilt by scripts/37b)."""
    checked = Counter()

    def ok(a, pos, aa):
        s = seq.get(a)
        good = bool(s) and 0 < pos <= len(s) and s[pos - 1] == aa
        checked['ok' if good else 'mismatch'] += 1
        return good

    rna, rna_n, rbd, dna = defaultdict(dict), defaultdict(Counter), defaultdict(list), defaultdict(dict)
    rna_cands, rna_single = defaultdict(set), defaultdict(set)
    dna_peps, dna_bind = defaultdict(list), defaultdict(list)
    bit_of = {label: bit for bit, (_, label) in enumerate(rna_sets)}
    for r in t.rows('pepr_ms_v2_table_s2_'):
        if not r['position']:
            continue
        a, pos, bit = r['uniprot_accession'], int(r['position']), bit_of[r['dataset']]
        if r['crosslinked_residue'] and ok(a, pos, r['crosslinked_residue']):
            rna[a][pos] = rna[a].get(pos, 0) | (1 << bit)
            rna_n[a][pos] += 1
            rna_cands[(bit, a, r['peptide'], r['rna_mass'])].add(pos)
    for r in t.rows('pepr_ms_v2_table_s3_'):
        a = r['uniprot_accession']
        focused = list(zip(L(r['focused_residues']), L(r['focused_positions'])))
        for aa, pos in focused:
            if ok(a, int(pos), aa):
                rna[a][int(pos)] = rna[a].get(int(pos), 0) | focused_bit
                if len(focused) == 1:
                    rna_single[a].add(int(pos))
    for (_, a, _, _), positions in rna_cands.items():
        if len(positions) == 1:
            rna_single[a].update(positions)
    gene_acc = {v[0]: a for a, v in sorted(proteome.items())}
    for r in t.rows('rbdmap_'):
        frag, a = r['fragment_sequence'], r['protein_id'].split('-')[0]
        if a not in seq and r['gene']:
            a = gene_acc.get(r['gene'], a)
        s = seq.get(a, '')
        try:
            s0, s1 = int(r['fragment_start']), int(r['fragment_end'])
        except ValueError:
            continue
        if not s or not frag:
            continue
        if s[s0 - 1:s1] != frag:
            hit = s.find(frag)
            if hit < 0 or s.find(frag, hit + 1) >= 0:
                continue
            s0, s1 = hit + 1, hit + len(frag)
        rbd[a].append((s0, s1))
    for r in t.rows('trendel2025_dna_crosslink_localisation'):
        try:
            start = int(r['peptide_start'])
        except ValueError:
            continue
        pep, a = r['sequence_with_localization'], r['uniprot_accession']
        bit = d_xl if sum(c.islower() for c in pep) == 1 else d_cand
        for k, c in enumerate(pep):
            if c.islower() and ok(a, start + k, c.upper()):
                dna[a][start + k] = dna[a].get(start + k, 0) | bit
    for r in t.rows('trendel2025_dna_crosslinked_peptides'):
        try:
            dna_peps[r['uniprot_accession']].append((int(r['peptide_start']), int(r['peptide_end'])))
        except ValueError:
            continue
    dna_int = {r['uniprot_accession'] for r in t.rows('trendel2025_mcf7_dna_interactome')}
    for r in t.rows('uniprot_human_dna_binding_regions_'):
        dna_bind[r['uniprot_accession']].append((int(r['start']), int(r['end'])))
    return rna, rna_n, rbd, dna, dna_peps, dna_int, dna_bind, rna_single, {'matched': checked['ok'], 'dropped': checked['mismatch']}


def crosslink_sequences(t: Tables, seq: dict, proteome: dict, rna_sets: list) -> dict:
    """The peptides behind the crosslinks, for showing them: rna[acc][pos] -> [(dataset bit, peptide, index of the
    crosslinked residue in it, RNA masses)], dna[acc] -> [(start, peptide with the candidate residues in lower case)]
    and rbd[acc] -> [(start, end, fragment)], each only where it matches the reference sequence."""
    bit_of = {label: bit for bit, (_, label) in enumerate(rna_sets)}
    rna = defaultdict(lambda: defaultdict(dict))
    for r in t.rows('pepr_ms_v2_table_s2_'):
        if not r['position'] or not r['crosslinked_residue']:
            continue
        a, pos, pep = r['uniprot_accession'], int(r['position']), r['peptide']
        # n[42.0106] marks an N-terminal acetyl group, not a residue
        plain = re.sub(r'([A-Z])\[(57\.02\d*|15\.99\d*)\]', r'\1', re.sub(r'^n\[[^\]]*\]', '', r['modified_sequence']))
        m = re.search(r'([A-Z])\[\d+\.\d+\]', plain)
        k = len(re.sub(r'\[[^\]]*\]', '', plain[:m.start()])) if m else -1
        s = seq.get(a, '')
        if k < 0 or pep[k:k + 1] != r['crosslinked_residue'] or not s or s[pos - 1 - k:pos - 1 - k + len(pep)] != pep:
            continue
        masses = rna[a][pos].setdefault((bit_of[r['dataset']], pep, k), [])
        if r['rna_mass'] not in masses:
            masses.append(r['rna_mass'])
    rna_out = {a: {pos: [(bit, pep, k, '/'.join(m)) for (bit, pep, k), m in sorted(v.items())] for pos, v in by.items()} for a, by in rna.items()}
    dna = defaultdict(set)
    for r in t.rows('trendel2025_dna_crosslink_localisation'):
        try:
            start = int(r['peptide_start'])
        except ValueError:
            continue
        pep, a = r['sequence_with_localization'], r['uniprot_accession']
        if seq.get(a, '')[start - 1:start - 1 + len(pep)] == pep.upper():
            dna[a].add((start, pep))
    gene_acc = {v[0]: a for a, v in sorted(proteome.items())}
    rbd = defaultdict(set)
    for r in t.rows('rbdmap_'):
        frag, a = r['fragment_sequence'], r['protein_id'].split('-')[0]
        if a not in seq and r['gene']:
            a = gene_acc.get(r['gene'], a)
        s = seq.get(a, '')
        if not s or not frag:
            continue
        hit = s.find(frag)
        if hit >= 0 and s.find(frag, hit + 1) < 0:
            rbd[a].add((hit + 1, hit + len(frag), frag))
    return {'rna': rna_out, 'dna': {a: sorted(v) for a, v in dna.items()}, 'rbd': {a: sorted(v) for a, v in rbd.items()}}


def merge_member(members: dict, acc: str, value: tuple) -> None:
    """A protein listed in a condensate under its accession and under isoforms (P29590 and P29590-12): one member with the
    highest confidence, every role and evidence bit, and every PubMed ID."""
    if acc in members:
        score, roles, evid, pmids = members[acc]
        value = (max(score, value[0]), roles | value[1], evid | value[2], pmids + [p for p in value[3] if p not in pmids])
    members[acc] = value


def cdcode(t: Tables, role_bits: dict, role_other: int, evidence_bits: dict, evidence_other: int):
    out = defaultdict(dict)
    for r in t.rows('cdcode_'):
        roles = evid = 0
        for x in L(r['roles']):
            roles |= role_bits.get(x, role_other)
        for e in L(r['experimental_evidence']):
            evid |= evidence_bits.get(e, evidence_other)
        score = float(r['confidence_score'])
        merge_member(out[r['condensate_uid']], r['uniprot_accession'].split('-')[0],
                     (int(score) if score.is_integer() else score, roles, evid, L(r['pmids'])))
    version = re.search(r'cdcode_([\d.]+)_', t.paths('cdcode_')[0].name).group(1)
    return dict(out), version


# ---- Measurements (structures and predicted models) ----------------------------------------------------------------

def contacts(t: Tables) -> dict:
    """(pdb, chain, SEQRES index) -> {'resn', 'hits': [(chain, seqres index, author number, residue, distance)]}."""
    out = {}
    for r in t.rows('pdb_arginine_contacts_'):
        rec = out.setdefault((r['pdb_id'], r['chain'], int(r['seqres_index'])), {'resn': r['residue'], 'hits': []})
        if r['contact_chain']:
            rec['hits'].append((r['contact_chain'], r['contact_seqres_index'], r['contact_author_number'], r['contact_residue'], float(r['distance_a'])))
    return out


def nucleic_chains(t: Tables) -> tuple[dict, set]:
    """(pdb, author chain) -> (kind, entity, molecule name, organism, SEQRES sequence, sample preparation, entry has core
    histones) of the RNA and DNA chains staged, and the PDB entries looked up (an entry with none has a row with an empty
    chain). Empty before v1.3.0's staging."""
    out, looked = {}, set()
    try:
        rows = list(t.rows('pdb_nucleic_acid_chains_'))
    except KeyError:
        return out, looked
    for r in rows:
        looked.add(r['pdb_id'])
        if r['chain']:
            out[(r['pdb_id'], r['chain'])] = (r['nucleic_acid'], int(r['entity_id']), r['molecule_name'], r['source_organism'], r['sequence'],
                                              r.get('sample_preparation', ''), r.get('entry_core_histones') == 'true')
    return out, looked


def na_arginines(t: Tables) -> dict:
    """(pdb, protein chain) -> {'args': [[SEQRES index, author number, residue]], 'hits': [[SEQRES index, nucleic-acid
    chain, closest nucleotide's SEQRES index, nucleotide, side-chain distance, [SEQRES indices within 8 Å]]]}, as
    scripts/37 chain_na_proximity returns it, for the chains read. Empty before these tables were staged."""
    out = {}
    try:
        chains = list(t.rows('pdb_chain_arginines_'))
        dists = list(t.rows('pdb_arginine_nucleic_acid_distances_'))
    except KeyError:
        return out
    for r in chains:
        out[(r['pdb_id'], r['chain'])] = {'args': [[int(x.split(':')[0]), x.split(':')[1], x.split(':')[2]] for x in r['arginines'].split('; ') if x], 'hits': []}
    for r in dists:
        out[(r['pdb_id'], r['chain'])]['hits'].append([int(r['seqres_index']), r['na_chain'], int(r['closest_na_seqres_index']), r['closest_nucleotide'],
                                                     float(r['side_chain_distance_a']), [int(x) for x in r['na_seqres_within_8a'].split('; ') if x]])
    return out


def na_residue_contacts(t: Tables) -> dict:
    """(pdb, protein chain) -> [(SEQRES index, author number, residue, nucleic-acid chain, its SEQRES index, nucleotide,
    distance)] for the protein chains read against RNA and DNA (a chain with none has an empty list)."""
    out = {}
    try:
        rows = list(t.rows('pdb_nucleic_acid_residue_contacts_'))
    except KeyError:
        return out
    for r in rows:
        hits = out.setdefault((r['pdb_id'], r['chain']), [])
        if r['na_chain']:
            hits.append((int(r['seqres_index']), r['author_number'], r['residue'], r['na_chain'], int(r['na_seqres_index']),
                         r['nucleotide'], float(r['distance_a'])))
    return out


def environments(t: Tables) -> dict:
    """(pdb, chain, SEQRES index) -> scripts/37 analyze_environment result, or None for a residue that is not an arginine."""
    out = {}
    for r in t.rows('pdb_arginine_interface_environments_'):
        k = (r['pdb_id'], r['chain'], int(r['seqres_index']))
        if not r['levy_class']:
            out[k] = None
            continue
        out[k] = {'resn': r['residue'], 'levy': r['levy_class'], 'rasa_mono': float(r['relative_asa_alone']),
                  'rasa_cplx': float(r['relative_asa_in_complex']), 'd_na': float(r['asa_buried_by_nucleic_acid_a2']),
                  'd_prot': float(r['asa_buried_by_protein_a2']), 'na_type': r['nucleic_acid_contact'],
                  'base_hb': L(r['base_hbond_partners']), 'n_prot_hb': int(r['protein_hbonds'])}
    return out


def model_contacts(t: Tables, g_clash: int, g_sb: int, g_hb: int, g_catpi: int, seqs: dict | None = None) -> dict:
    """(source, model) -> {chain: {position: measure}} with measure = (min distance, flags, pLDDT R, lowest partner pLDDT,
    [[partner chain, residue label], ...]), as scripts/37 measure_methyl_r gives it. With seqs (the reviewed proteome), a
    measure is left out when its partner is outside it or its partner residues disagree with the partner's reference
    sequence: the chain named for the partner is another protein, or numbers another sequence. An AlphaFold DB row staged
    without its accessions (the top-ups before 2026-10-22 did not write them) takes them from the model's chains in the
    heterodimer tables."""
    out = defaultdict(lambda: defaultdict(dict))
    models = None
    for r in t.rows('predicted_model_arginine_contacts_'):
        if r['source'] == 'AlphaFold DB' and not (r['uniprot_accession'] and r['partner_accession']):
            models = afdb_heterodimers(t) if models is None else models
            x = models.get(r['model'], {})
            acc_of = dict(zip(x.get('chains', []), x.get('accessions', [])))
            if r['chain'] in acc_of and len(acc_of) == 2:
                r = dict(r, uniprot_accession=r['uniprot_accession'] or acc_of[r['chain']],
                         partner_accession=r['partner_accession'] or next(a for c, a in acc_of.items() if c != r['chain']))
        if seqs is not None:
            partner = seqs.get(r['partner_accession'], '')
            labels = [x.split(':')[1] for x in L(r['partner_residues']) if x.count(':') >= 1]   # 'A:E240:hb+sb' -> 'E240'
            if not partner or (labels and not labels_agree(labels, partner)) \
                    or seqs.get(r['uniprot_accession'], '')[int(r['position']) - 1:int(r['position'])] != 'R':
                continue
        flags = (r['clash'] == 'true') * g_clash | (r['salt_bridge'] == 'true') * g_sb | (r['hydrogen_bond'] == 'true') * g_hb | (r['cation_pi'] == 'true') * g_catpi
        labels = [list(x.split(':', 1)) for x in L(r['partner_residues'])]
        m = [float(r['min_distance_a']), flags, int(r['plddt_arginine']) if r['plddt_arginine'] else None,
             int(r['plddt_partner_min']) if r['plddt_partner_min'] else None, labels]
        out[(r['source'], r['model'])][r['chain']][str(r['position'])] = m
    return out


def model_positions(t: Tables) -> dict:
    """(source, model, accession) -> positions a top-up measured in that model (empty without top-up tables)."""
    out = defaultdict(set)
    if not any(n.startswith('predicted_model_measured_positions_') for n in t.names):
        return out
    for r in t.rows('predicted_model_measured_positions_'):
        out[(r['source'], r['model'], r['uniprot_accession'])].update(int(x) for x in L(r['positions']))
    return out


def afdb_heterodimers(t: Tables) -> dict:
    """Model id -> accessions, chains and scores, the later tables over the earlier. A top-up lists the models of a
    query without reading each, so a later row can lack the chains an earlier table gave: the model keeps those (with
    the accessions in their order), else scripts/37 would look its chains up again on every run."""
    out = {}
    for r in t.rows('alphafold_db_human_heterodimers_'):
        new = {'accessions': L(r['uniprot_accessions']), 'chains': L(r['chains']),
               'iptm': float(r['iptm']) if r['iptm'] else None, 'pdockq': float(r['pdockq']) if r['pdockq'] else None}
        old = out.get(r['model_id'])
        if old and old['chains'] and not new['chains']:
            new.update(accessions=old['accessions'], chains=old['chains'], iptm=new['iptm'] if new['iptm'] is not None else old['iptm'],
                       pdockq=new['pdockq'] if new['pdockq'] is not None else old['pdockq'])
        out[r['model_id']] = new
    return out


def afdb_queried(t: Tables) -> dict:
    """Accession -> the number of human heterodimer models its AlphaFold DB query listed."""
    return {r['uniprot_accession']: int(r['human_heterodimers'] or 0) for r in t.rows('alphafold_db_complex_queries_')}


def read_json(path: Path, default):
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return default


# ---- Staging -----------------------------------------------------------------------------------------------------------

def write_tsv(path: Path, header: list[str], rows) -> int:
    n = 0
    with open(path, 'w', encoding='utf-8', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(header)
        for r in rows:
            w.writerow(['; '.join(map(str, v)) if isinstance(v, list) else v for v in r])
            n += 1
    return n


def interface_entries(t: Tables) -> set:
    """PDB entries any interface, contact, environment or modelled-residue table names."""
    out = set()
    for r in t.rows('pdbe_kb_interface_residues_'):
        out.update(L(r['pdb_entries']))
    for prefix in ('pdb_methylarginine_', 'pdb_arginine_contacts_', 'pdb_arginine_interface_environments_'):
        out.update(r['pdb_id'] for r in t.rows(prefix))
    return out


def stage_sifts_other(t: Tables, raw_gz: Path, out: Path) -> int:
    """The SIFTS chain segments of those PDB entries that the reviewed-human table leaves out (other organisms,
    unreviewed entries), so a structure contact with such a chain names its protein as scripts/37 does from the raw file."""
    import gzip
    entries = interface_entries(t)
    have = {(r['pdb_id'], r['chain'], r['uniprot_accession'], r['seqres_start'], r['uniprot_start'], r['uniprot_end'])
            for r in t.rows('sifts_pdb_chain_uniprot_human_')}
    rows = []
    with gzip.open(raw_gz, 'rt', encoding='utf-8') as fh:
        header = None
        for r in csv.reader(fh, delimiter='\t'):
            if not r or r[0].startswith('#'):
                continue
            if header is None:
                header = {h.strip(): i for i, h in enumerate(r)}
                continue
            pdb, beg = r[header['PDB']].lower(), r[header['RES_BEG']]
            if pdb not in entries or not re.fullmatch(r'\d+', beg):
                continue
            try:
                row = (pdb, r[header['CHAIN']], r[header['SP_PRIMARY']], str(int(beg)), str(int(r[header['SP_BEG']])), str(int(r[header['SP_END']])))
            except ValueError:
                continue
            if row not in have:
                rows.append(row)
    rows.sort(key=lambda x: (x[0], x[1], int(x[4])))
    return write_tsv(out, ['pdb_id', 'chain', 'uniprot_accession', 'seqres_start', 'uniprot_start', 'uniprot_end'], rows)

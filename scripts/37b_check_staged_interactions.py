#!/usr/bin/env python3
"""Check that the tables in inputs/interactions/ hold what scripts/37 reads from its download cache.

For each source this runs the scripts/37 loader on the cache, rebuilds the same result from the tables alone and
compares the two. Network access is switched off for the run, so every loader reads only saved files. The measurement
tables (structure contacts, interface geometry, predicted-model contacts) are compared with the website data, and
every table's MD5 with MANIFEST.tsv. It needs the cache, so it runs where the downloads are kept:

    python scripts/37b_check_staged_interactions.py --cache PTM_data/interactions_cache --tables inputs/interactions
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import re
import sys
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(name, file):
    spec = importlib.util.spec_from_file_location(name, HERE / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


x37 = load('export37', '37_export_interaction_data.py')
x37a = load('stage37a', '37a_stage_interactions.py')
csv.field_size_limit(1 << 30)


def no_network(*_a, **_k):
    raise RuntimeError('network call during the check: a loader wanted a file that is not in the cache')


urllib.request.urlopen = no_network


def L(v: str) -> list[str]:
    return v.split('; ') if v else []


def expand(ranges: str) -> set[int]:
    out = set()
    for part in filter(None, ranges.split(',')):
        a, _, b = part.partition('-')
        out.update(range(int(a), int(b or a) + 1))
    return out


def merge(spans) -> list[tuple[int, int]]:
    out = []
    for a, b in sorted(spans):
        if out and a <= out[-1][1] + 1:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


class Tables:
    def __init__(self, folder: Path):
        self.folder = folder
        with open(folder / 'MANIFEST.tsv', encoding='utf-8') as fh:
            self.manifest = list(csv.DictReader(fh, delimiter='\t'))
        self.names = [r['file'].split('/', 1)[1] for r in self.manifest]

    def path(self, prefix: str, suffix: str = '.tsv') -> Path:
        hits = [n for n in self.names if n.startswith(prefix) and n.endswith(suffix)]
        if len(hits) != 1:
            raise KeyError(f'{prefix}: {hits}')
        return self.folder / hits[0]

    def rows(self, prefix: str):
        with open(self.path(prefix), encoding='utf-8') as fh:
            yield from csv.DictReader(fh, delimiter='\t')


results = []


def check(name: str, ok: bool, detail: str):
    results.append((name, ok, detail))
    print(f'{"ok  " if ok else "FAIL"} {name}: {detail}', flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cache', default='PTM_data/interactions_cache')
    ap.add_argument('--tables', default='inputs/interactions')
    ap.add_argument('--site-json', default='website/data/methylarg_interactions.json')
    args = ap.parse_args()
    c, t = Path(args.cache), Tables(Path(args.tables))
    site = json.loads(Path(args.site_json).read_text())

    # ------------------------------------------------------------------ files
    bad = [r['file'] for r in t.manifest if x37a.md5_file(t.folder / r['file'].split('/', 1)[1]) != r['md5']]
    check('MANIFEST checksums', not bad, f'{len(t.manifest)} files' + (f'; differ: {bad}' if bad else ''))

    # ------------------------------------------------------------------ UniProt
    proteome, classes = x37.load_uniprot(c, False)
    seqs, acc = {}, None
    with open(t.path('uniprot_human_reviewed_', '.fasta'), encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('>'):
                acc = line[1:].split()[0]
                seqs[acc] = []
            else:
                seqs[acc].append(line.strip())
    seqs = {a: ''.join(v) for a, v in seqs.items()}
    rows = list(t.rows('uniprot_human_reviewed_'))
    diff = [r['uniprot_accession'] for r in rows
            if r['uniprot_accession'] not in proteome
            or L(r['gene_names_primary']) != [g.strip() for g in proteome[r['uniprot_accession']][0].split(';') if g.strip()]
            or r['string_id'] != proteome[r['uniprot_accession']][1]
            or seqs.get(r['uniprot_accession']) != proteome[r['uniprot_accession']][2]
            or int(r['sequence_length']) != len(proteome[r['uniprot_accession']][2])
            or r['sequence_md5'] != hashlib.md5(proteome[r['uniprot_accession']][2].encode()).hexdigest().upper()]
    rebuilt = {k: {r['uniprot_accession'] for r in rows if r[col] == 'true'}
               for k, col in (('tudor', 'tudor_domain'), ('wd40', 'wd40_repeat'), ('phd', 'phd_finger'))}
    check('UniProt entries, genes, STRING ids, sequences', len(rows) == len(proteome) == len(seqs) and not diff,
          f'{len(rows):,} entries, {len(seqs):,} sequences' + (f'; differ: {diff[:5]}' if diff else ''))
    check('UniProt Tudor / WD40 / PHD classes', rebuilt == classes, ', '.join(f'{k} {len(v)}' for k, v in rebuilt.items()))
    seq = {a: v[2] for a, v in proteome.items()}
    reviewed = set(proteome)

    dis, dom = x37.load_uniprot_regions(c, False)
    dis2, dom2 = defaultdict(list), defaultdict(list)
    for r in t.rows('uniprot_human_regions_'):
        (dis2 if r['feature'] == 'disordered_region' else dom2)[r['uniprot_accession']].append((int(r['start']), int(r['end'])))
    same = ({k: sorted(v) for k, v in dis.items()} == {k: sorted(v) for k, v in dis2.items()}
            and {k: sorted(v) for k, v in dom.items()} == {k: sorted(v) for k, v in dom2.items()})
    check('UniProt disordered regions and domains', same,
          f'{sum(map(len, dis2.values())):,} disordered regions, {sum(map(len, dom2.values())):,} domains/repeats/zinc fingers')

    # ------------------------------------------------------------------ GO
    go = x37.load_go(c, False)
    go2 = defaultdict(lambda: {'rna': set(), 'dna': set()})
    for r in t.rows('quickgo_'):
        go2[r['uniprot_accession']][r['binding']].add(r['evidence_code'])
    check('QuickGO RNA/DNA binding evidence per protein', dict(go) == dict(go2), f'{len(go2):,} proteins')

    # ------------------------------------------------------------------ STRING
    sid_to_acc = {}
    for a, (_, sid, _) in proteome.items():
        if sid and sid not in sid_to_acc:
            sid_to_acc[sid] = a
    edges = x37.load_string(c, False, sid_to_acc, reviewed, 400)
    edges2 = {(r['protein_a'], r['protein_b']): int(r['combined_score']) for r in t.rows('string_v')}
    n700 = sum(v >= 700 for v in edges2.values())
    check('STRING physical links (score >= 400, any reviewed protein)', edges == edges2,
          f'{len(edges2):,} links, {n700:,} at >= 700')

    # ------------------------------------------------------------------ Complex Portal
    cx = [(d['id'], x37a.text(d['name']), d['predicted'], d['slots']) for d in x37.load_complexes(c, False)]
    cx2 = {}
    for r in t.rows('complex_portal_'):
        key = (r['complex_id'], r['curation'])
        d = cx2.setdefault(key, [r['complex_id'], r['complex_name'], r['curation'] == 'predicted', []])
        k = int(r['member_slot']) - 1
        while len(d[3]) <= k:
            d[3].append([])
        d[3][k].append(r['uniprot_accession'])
    cx2 = [tuple(v) for v in cx2.values()]
    check('Complex Portal complexes and member slots', cx == cx2,
          f'{sum(not x[2] for x in cx2):,} curated, {sum(x[2] for x in cx2):,} predicted')

    # ------------------------------------------------------------------ PDBe-KB and SIFTS
    q = list(t.rows('pdbe_kb_interface_queries_'))
    ifaces = x37.fetch_interfaces(c, [r['uniprot_accession'] for r in q], False, 1)
    want = set()
    for a, data in ifaces.items():
        for part in data:
            for r in part.get('residues', []):
                if r.get('indexType', 'UNIPROT') != 'UNIPROT':
                    continue
                pdbs = frozenset(e['pdbId'] for e in r.get('interactingPDBEntries') or [] if e.get('pdbId'))
                for pos in range(int(r['startIndex']), int(r['endIndex']) + 1):
                    want.add((a, part.get('accession') or '', x37a.text(part.get('name')), pos, pdbs))
    got = set()
    for r in t.rows('pdbe_kb_interface_residues_'):
        for pos in range(int(r['start']), int(r['end']) + 1):
            got.add((r['uniprot_accession'], r['partner'], r['partner_name'], pos, frozenset(L(r['pdb_entries']))))
    check('PDBe-KB interface residues', want == got,
          f'{len(q):,} proteins looked up, {len(ifaces):,} with interfaces, {len(got):,} residue-partner pairs')

    obs = x37.load_observed(c, False, reviewed)
    obs2 = defaultdict(list)
    for r in t.rows('sifts_observed_ranges_'):
        obs2[r['uniprot_accession']].append((int(r['start']), int(r['end'])))
    check('SIFTS observed residues', {k: merge(v) for k, v in obs.items()} == {k: merge(v) for k, v in obs2.items()},
          f'{len(obs2):,} proteins')

    chains_rows = list(t.rows('sifts_pdb_chain_uniprot_'))
    ch = x37.load_sifts_chains(c, False, {r['pdb_id'] for r in chains_rows})
    want = Counter((p, *x) for p, xs in ch.items() for x in xs if x37.base_acc(x[1]) in reviewed)
    got = Counter((r['pdb_id'], r['chain'], r['uniprot_accession'], int(r['seqres_start']), int(r['uniprot_start']),
                   int(r['uniprot_end'])) for r in chains_rows)
    check('SIFTS chain segments (reviewed human chains)', want == got, f'{len(got):,} segments in {len(ch):,} entries')

    modelled = x37.load_modelled_methyl_r(c, False, 1)
    got = [(r['pdb_id'], r['chain'], int(r['seqres_index']), r['chem_comp']) for r in t.rows('pdb_methylarginine_')]
    check('PDB-modelled methylarginines', [(p, ch_ or '', i, cc) for p, ch_, i, cc in modelled] == got,
          f'{len(got)} residues in {len({g[0] for g in got})} entries')

    # ------------------------------------------------------------------ ProtVar
    pv = x37.load_protvar(c, False, reviewed)
    pv2 = defaultdict(list)
    for r in t.rows('protvar_predicted_'):
        a, b = x37.base_acc(r['uniprot_id1']), x37.base_acc(r['uniprot_id2'])
        score = float(r['pdockq'])
        if a == b or score < x37.PRED_TIERS[0] or (a not in reviewed and b not in reviewed):
            continue
        for acc_, other, field, first in ((a, b, 'ifresid1', True), (b, a, 'ifresid2', False)):
            if acc_ not in reviewed:
                continue
            res = [x for x in r[field].split(',') if x]
            pos_all = {int(x[1:]) for x in res if x[1:].isdigit()}
            pos_r = {int(x[1:]) for x in res if x[:1] == 'R' and x[1:].isdigit()}
            if pos_all:
                pv2[acc_].append((other, score, pos_r, pos_all, first))
    check('ProtVar predicted interfaces', dict(pv) == dict(pv2),
          f'{sum(map(len, pv2.values())):,} interface sides on {len(pv2):,} reviewed proteins')

    served = x37.load_protvar_served(c, False)
    raw = {}
    with x37.open_text(c / Path(x37.PROTVAR_SERVED).name) as fh:
        rd = csv.reader(fh, delimiter='\t')
        ix = {h: k for k, h in enumerate(next(rd))}
        for row in rd:
            raw[row[ix['interaction_id']]] = ({int(v) for v in row[ix['ifresid1']].split(',') if v.isdigit()},
                                              {int(v) for v in row[ix['ifresid2']].split(',') if v.isdigit()}, row[ix['pdockq']])
    rows = list(t.rows('protvar_served_'))
    diff = [r['interaction_id'] for r in rows
            if raw.get(r['interaction_id']) != (expand(r['ifresid1_ranges']), expand(r['ifresid2_ranges']), r['pdockq'])]
    check('ProtVar served models and interface residues', {(r['uniprot_id1'], r['uniprot_id2']) for r in rows} == served
          and len(rows) == len(raw) and not diff, f'{len(rows):,} pairs' + (f'; differ: {diff[:5]}' if diff else ''))

    # ------------------------------------------------------------------ IntAct
    feats, counts, _ = x37.load_intact(c, False, seq)
    feats2 = defaultdict(list)
    part_re = re.compile(r'([^:]+):(.*) \(([^()]*)\)$')
    for kind, label in ((0, 'mutations'), (1, 'binding_regions'), (2, 'methylarginines')):
        for r in t.rows(f'intact_{label}_human_'):
            if r['sequence_check'] not in ('matches', 'no curated sequence'):
                continue
            a = r['uniprot_accession']
            mi = r['feature_type_mi']
            code = 0 if kind == 2 else x37.INTACT_EFFECT_MI.get(mi, 5) if kind == 0 else x37.INTACT_REGION_MI.get(mi)
            spans = [tuple(map(int, x.split('-'))) for x in r['ranges'].split(',')]
            parts = [part_re.match(p).groups() for p in L(r['interaction_participants'])]
            own = [p for p in parts if p[0] == 'uniprotkb' and p[1].split('-')[0] == a]
            others = [p for p in parts if not (p[0] == 'uniprotkb' and p[1].split('-')[0] == a)]
            partners = tuple(others) if others else (('uniprotkb', a, 'protein'),) if len(own) > 1 or len(parts) == 1 else ()
            pm = L(r['pubmed'])
            resulting = r['resulting_sequence'] if kind == 0 and r['resulting_sequence'] not in ('', '-') else ''
            for k, (s0, s1) in enumerate(spans):
                res = resulting.split(',')[k] if resulting and len(resulting.split(',')) == len(spans) else resulting
                feats2[a].append((kind, code, s0, s1, res, partners, int(pm[0]) if pm else 0, r['interaction_ac']))
    check('IntAct features scripts/37 uses', dict(feats) == dict(feats2), f'{sum(map(len, feats2.values())):,} features on {len(feats2):,} proteins')
    check('IntAct counts equal the website build', counts == site['meta']['intact']['counts'],
          ', '.join(f'{k} {v:,}' for k, v in sorted(counts.items())))

    # ------------------------------------------------------------------ crosslinks
    xl = x37.load_crosslinks(c, False, seq, proteome)
    checked = Counter()

    def ok(a, pos, aa):
        s = seq.get(a)
        good = bool(s) and 0 < pos <= len(s) and s[pos - 1] == aa
        checked['ok' if good else 'mismatch'] += 1
        return good

    rna, rna_n, rbd, dna = defaultdict(dict), defaultdict(Counter), defaultdict(list), defaultdict(dict)
    rna_cands, rna_single = defaultdict(set), defaultdict(set)
    dna_peps, dna_bind, dna_int = defaultdict(list), defaultdict(list), set()
    bit_of = {label: bit for bit, (_, label) in enumerate(x37.RNA_XL_SETS)}
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
                rna[a][int(pos)] = rna[a].get(int(pos), 0) | x37.RNA_FOCUSED
                if len(focused) == 1:
                    rna_single[a].add(int(pos))
    for (_, a, _, _), positions in rna_cands.items():
        if len(positions) == 1:
            rna_single[a].update(positions)
    gene_acc = {v[0]: a for a, v in sorted(proteome.items())}
    for r in t.rows('rbdmap_'):
        frag, a = r['fragment_sequence'], x37.base_acc(r['protein_id'])
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
        bit = x37.D_XL if sum(ch_.islower() for ch_ in pep) == 1 else x37.D_CAND
        for k, ch_ in enumerate(pep):
            if ch_.islower() and ok(a, start + k, ch_.upper()):
                dna[a][start + k] = dna[a].get(start + k, 0) | bit
    for r in t.rows('trendel2025_dna_crosslinked_peptides'):
        try:
            dna_peps[r['uniprot_accession']].append((int(r['peptide_start']), int(r['peptide_end'])))
        except ValueError:
            continue
    dna_int = {r['uniprot_accession'] for r in t.rows('trendel2025_mcf7_dna_interactome')}
    for r in t.rows('uniprot_human_dna_binding_regions_'):
        dna_bind[r['uniprot_accession']].append((int(r['start']), int(r['end'])))
    rebuilt = (rna, rna_n, rbd, dna, dna_peps, dna_int, dna_bind, rna_single, {'matched': checked['ok'], 'dropped': checked['mismatch']})
    names = ('RNA crosslinked residues', 'RNA crosslink counts', 'RBDmap fragments', 'DNA crosslinked residues',
             'DNA crosslinked peptides', 'DNA interactome', 'UniProt DNA-binding regions', 'single-localisation RNA residues',
             'sequence check')
    plain = lambda v: {k: dict(x) if isinstance(x, (dict, Counter)) else x for k, x in v.items()} if isinstance(v, dict) else v
    for name, a, b in zip(names, xl, rebuilt):
        size = sum(len(v) for v in b.values()) if isinstance(b, dict) and name != 'sequence check' else len(b) if not isinstance(b, dict) else b
        check(f'Crosslinks: {name}', plain(a) == plain(b), f'{size:,}' if isinstance(size, int) else str(size))
    check('Crosslink sequence check equals the website build', xl[8] == site['meta']['xl_check'], str(site['meta']['xl_check']))

    # ------------------------------------------------------------------ CD-CODE
    uids = sorted(p.stem for p in (c / 'cdcode').glob('*.json') if p.name != 'release.json')
    cd, version = x37.load_cdcode(c, False, uids, 1)
    cd2 = defaultdict(dict)
    for r in t.rows('cdcode_'):
        roles = evid = 0
        for x in L(r['roles']):
            roles |= x37.ROLE_BITS.get(x, x37.ROLE_OTHER)
        for e in L(r['experimental_evidence']):
            evid |= x37.EVIDENCE_BITS.get(e, x37.EVIDENCE_OTHER)
        cd2[r['condensate_uid']][r['uniprot_accession']] = (float(r['confidence_score']), roles, evid, L(r['pmids']))
    want = {u: {a: (float(m[0]), m[1], m[2], m[3]) for a, m in v.items()} for u, v in cd.items() if v}
    check(f'CD-CODE {version} members', want == dict(cd2), f'{len(cd2)} condensates with members ({len(uids)} saved)')

    # ------------------------------------------------------------------ measurements against the website data
    rows = list(t.rows('pdb_arginine_contacts_'))
    hoods = {(r['pdb_id'], r['chain'], r['seqres_index']) for r in rows}
    check('PDB arginine contacts: neighbourhoods read', len(hoods) > 0, f'{len(hoods):,} neighbourhoods, {sum(bool(r["contact_chain"]) for r in rows):,} contacts '
          f'(website: {len(site["contacts"]):,} sites with residue-level contacts)')

    env = list(t.rows('pdb_arginine_interface_environments_'))
    levy = {'none': 0, 'support': 1, 'rim': 2, 'core': 3}
    pool = defaultdict(list)
    for r in env:
        if r['levy_class']:
            pool[(r['pdb_id'], r['chain'], levy[r['levy_class']])].append((float(r['relative_asa_alone']), float(r['relative_asa_in_complex'])))
    missing = [g for g in site['geometry']
               if not any(abs(x - g[4]) < 0.002 and abs(y - g[5]) < 0.002 for x, y in pool.get((g[1], g[2], g[3]), []))]
    check('Interface geometry: every website site has its table row', not missing,
          f'{len(env):,} environments; {len(site["geometry"]):,} website sites' + (f'; unmatched: {missing[:3]}' if missing else ''))

    flags = lambda r: (r['salt_bridge'] == 'true') * x37.G_SB | (r['hydrogen_bond'] == 'true') * x37.G_HB | \
        (r['cation_pi'] == 'true') * x37.G_CATPI | (r['clash'] == 'true') * x37.G_CLASH
    got = Counter((0 if r['source'] == 'ProtVar' else 1, r['chain'], round(float(r['min_distance_a']) * 10), flags(r),
                   int(r['plddt_arginine']) if r['plddt_arginine'] else None,
                   int(r['plddt_partner_min']) if r['plddt_partner_min'] else None)
                  for r in t.rows('predicted_model_arginine_contacts_'))
    want = Counter((g[2], g[11], g[7], g[8], g[5], g[6]) for g in site['gallery'])
    check('Predicted-model contacts equal the website gallery', got == want, f'{sum(got.values()):,} contacts')

    afdb = list(t.rows('alphafold_db_human_heterodimers_'))
    check('AlphaFold DB human heterodimer models', len(afdb) > 0, f'{len(afdb):,} models')

    failed = [n for n, ok_, _ in results if not ok_]
    print(f'\n{len(results) - len(failed)} of {len(results)} checks passed' + (f'; failed: {failed}' if failed else ''))
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()

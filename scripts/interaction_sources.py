"""The provenance layer of the interaction records (scripts/release_layers.py SOURCES): for every record at a site, the
experiment behind it (PDB entry, IntAct interaction, crosslinking study or predicted model), its method, species and
paper, what the record says the partner does there, and the state of the site in that experiment.

The records are those of supporting/site_interactions (scripts/37: PDBe-KB interfaces, IntAct features, crosslinks and
RNA-bound peptides in cells, predicted complexes). Each PDB and IntAct row is split into its PDB entries or IntAct
interactions and joined to what scripts/53_stage_interaction_sources.py staged in inputs/interactions/:

  pdb_entries, pdb_entities         the entry's title, method, resolution and primary citation; each polymer entity's
                                    description, UniProt accessions and source organisms
  pdb_site_residues                 the residue each entity aligned to the site's protein models at the site (ARG, a
                                    methylarginine 2MR/DA2/NMM, or another residue in a mutated or variant construct),
                                    and the other modified residues of that entity
  pdb_methylated_sites, *_contacts  the entries modelling the site as a methylarginine, aligned or matched by sequence
                                    (orthologues, peptides, histone variants sharing a tail), and the residues of other
                                    chains within 4 Å of its side chain. An entry PDBe-KB does not list for the site
                                    adds records of its own (evidence pdb_methylated_contact)
  intact_interactions, intact_terms the binary IntAct records: both participants with species, their features,
                                    detection method, interaction type, host organism and publication

A record's class says what it shows (CLASSES): the partner methylates or demethylates the site (an IntAct methylation
or demethylation reaction); binds it methylated, unmodified or with another residue there (the residue a structure
models; IntAct curating the residue methylated); mutating the site changes the binding or does not; binds a region
holding it with the state not recorded; another reaction; RNA or DNA; a partner that is not a protein (a gene locus
in chromatin immunoprecipitation, a small molecule); another copy of the protein; a predicted model. Flags (FLAGS) mark
what a reader should check: a partner or site protein from another species, a structure matched by sequence, a
methylation state the release does not report at the site, a paper whose title names another member of the
partner's family and not the partner, no publication.

Outputs: supporting/site_interaction_records (one row per site and record), five site columns, and the website's
data/site_partners.json (and site_partners_all.json from the view with lower-confidence sites).
"""
from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

DATE = '2026-10-30'
FILES = [f'{k}_{DATE}.tsv' for k in ('pdb_entries', 'pdb_entities', 'pdb_site_residues', 'pdb_methylated_sites', 'pdb_methylated_site_contacts',
                                    'intact_interactions', 'intact_terms', 'pubmed_papers', 'uniprot_partner_names')] + [f'interaction_sources_summary_{DATE}.json']
INTACT_FEATURES = '2026-10-19'   # the feature tables of scripts/37a (intact_methylarginines_human_*, intact_mutations_human_*)
HUMAN = 'Homo sapiens'

CLASSES = {
    'writes': 'methylates the site: an IntAct methylation reaction with the partner as the enzyme and the methylarginine curated at the site, '
              'or mutating the site changes the reaction',
    'removes': 'demethylates the site: an IntAct demethylation reaction with the partner as the enzyme and the methylarginine curated at the site',
    'binds_methylated': 'binds with the site methylated in the experiment: a PDB structure modelling a methylarginine there (2MR, DA2 or NMM), '
                        'or an IntAct binding record that curates the residue as methylated',
    'binds_unmodified': 'binds with an unmodified arginine at the site: a PDB structure modelling ARG there',
    'binds_other_residue': 'binds with another residue at the site: a PDB structure of a mutated or variant construct (state_detail names it)',
    'binds_state_not_read': 'a PDB interface at the site whose residue there could not be read (no entity of the entry aligns to the site)',
    'mutation_changes_binding': 'mutating the site decreases, disrupts or increases the binding (IntAct)',
    'mutation_no_effect': 'mutating the site leaves the binding unchanged (IntAct)',
    'mutation_unspecified': 'a mutation of the site was tested, its effect on the binding not stated (IntAct feature type mutation)',
    'binds_region': 'binds a region holding the site, the state of the site not recorded (IntAct binding region)',
    'other_reaction': 'an enzymatic reaction on a region holding the site that IntAct does not place at the site: acetylation, phosphorylation '
                      'or cleavage, a methylation or demethylation of another residue (a lysine methyltransferase on a histone tail), one where '
                      'mutating the site has no effect, or one in which the site\'s protein is the enzyme',
    'nucleic_acid': 'RNA or DNA: contact in a PDB structure, a crosslink or RNA-bound peptide in cells, an IntAct nucleic-acid participant, '
                    'or a UniProt DNA-binding region (annotation)',
    'not_a_protein': 'the IntAct partner is not a protein: a gene locus (chromatin immunoprecipitation), a small molecule or another participant',
    'self': 'another copy of the same protein in a PDB structure, or the protein with itself in IntAct',
    'predicted': 'a predicted complex model only (ProtVar, AlphaFold DB), not an observation',
}
# the classes in a sentence (CHANGES.md)
CLASS_WORDS = {'writes': 'methylates the site', 'removes': 'demethylates it', 'binds_methylated': 'binds it methylated',
               'binds_unmodified': 'binds it unmodified', 'binds_other_residue': 'binds a mutant or variant there',
               'binds_state_not_read': 'at an interface whose residue could not be read', 'mutation_changes_binding': 'mutating it changes the binding',
               'mutation_no_effect': 'mutating it leaves the binding unchanged', 'mutation_unspecified': 'mutated, effect not stated',
               'binds_region': 'binds a region holding it', 'other_reaction': 'another reaction on the region', 'nucleic_acid': 'RNA or DNA',
               'not_a_protein': 'not a protein partner', 'self': 'another copy of the protein', 'predicted': 'predicted only'}
BINDS = ('binds_methylated', 'binds_unmodified', 'binds_other_residue', 'binds_state_not_read', 'binds_region')
STATES = {
    'arginine': 'an unmodified arginine (PDB residue ARG)',
    'Rme1': 'monomethylarginine (PDB NMM, PSI-MOD MOD:00414)',
    'Rme2a': 'asymmetric dimethylarginine (PDB DA2, PSI-MOD MOD:00077)',
    'Rme2s': 'symmetric dimethylarginine (PDB 2MR, PSI-MOD MOD:00076)',
    'Rme2': 'dimethylarginine, asymmetric or symmetric not stated (PSI-MOD MOD:00783)',
    'Rme': 'methylarginine, degree not stated (PSI-MOD MOD:00658)',
    'other_residue': 'another residue in the construct (state_detail names it)',
    'mutated': 'mutated in the experiment (state_detail gives the change)',
    'not_read': 'a structure whose residue at the site could not be read',
    'not_recorded': 'not recorded: an IntAct region, or a crosslink, prediction or annotation',
}
FLAGS = {
    'partner_other_species': 'the partner in the experiment is from another species (its PDB entity or IntAct participant)',
    'site_protein_other_species': 'the site\'s protein in the experiment is from another species: an orthologue identical around the site',
    'matched_by_sequence': 'the PDB entity is not aligned to this protein: its sequence around the site is identical (a peptide, an '
                           'orthologue, or another histone variant sharing the tail)',
    'state_not_reported': 'the methylation the experiment records (mono- or dimethyl, asymmetric or symmetric) is not among the states the '
                          'release reports at the site',
    'paper_names_other_arginine': 'the paper\'s title and abstract name other arginines by number and not this one (in UniProt or histone numbering): '
                                  'check that IntAct placed the feature on the residue the paper studied',
    'title_names_paralogue': 'the paper\'s title names another member of the partner\'s gene family and not the partner (IntAct records): check '
                             'which protein the paper studied',
    'state_unexpected_for_enzyme': 'a methylation reaction whose product is not the kind the enzyme\'s PRMT type makes (type I asymmetric, type II '
                                   'symmetric, type III monomethyl only)',
    'no_publication': 'the record gives no PubMed or DOI',
    'negative': 'IntAct records the interaction as not observed',
}
PARTNER_KINDS = {'protein': 'a protein', 'self': 'another copy of the site\'s protein', 'rna': 'RNA', 'dna': 'DNA',
                 'peptide': 'a peptide or complex without a UniProt entry (IntAct: a synthetic peptide, a Complex Portal complex)',
                 'gene': 'a gene locus (IntAct participant type gene, from chromatin immunoprecipitation)',
                 'small_molecule': 'a small molecule', 'other': 'another participant type'}

PDB_STATE = {'ARG': 'arginine', 'NMM': 'Rme1', 'DA2': 'Rme2a', '2MR': 'Rme2s'}
METHYL_LABEL = {'me1arg': 'Rme1', 'nonome2arg': 'Rme2a', "nono'me2arg": 'Rme2s', 'nnme2arg': 'Rme2', 'mearg': 'Rme'}
# residue marks shown beside a record: PDB chemical components (residue, suffix) and IntAct feature labels (PSI-MOD).
# Others the entries model (selenomethionine, oxidised or alkylated cysteines, caps, cofactors) are not marks.
MARKS = {'SEP': 'Sph', 'TPO': 'Tph', 'PTR': 'Yph', 'NEP': 'Hph', 'HIP': 'Hph', 'ALY': 'Kac', 'MLZ': 'Kme1', 'MLY': 'Kme2', 'M3L': 'Kme3',
         'KCR': 'Kcr', 'BTK': 'Kbu', 'SLL': 'Ksucc', 'KHB': 'Kbhb', 'MYK': 'Kmyr', 'NMM': 'Rme1', 'DA2': 'Rme2a', '2MR': 'Rme2s', 'CIR': 'Rcit',
         'HIC': 'Hme', '5CT': 'Khyp', 'HY3': 'Poh'}
LABEL_MARKS = {'n6aclys': 'ac', 'naclys': 'ac', 'ac1res': 'ac', 'n6me1lys': 'me1', 'n6me2lys': 'me2', 'n6me3+lys': 'me3', 'melys': 'me',
               'opser': 'ph', 'opthr': 'ph', 'optyr': 'ph', 'phosres': 'ph', 'meres': 'me', **{k: v.replace('R', '', 1) for k, v in METHYL_LABEL.items()}}
THREE = {'ala': 'A', 'arg': 'R', 'asn': 'N', 'asp': 'D', 'cys': 'C', 'gln': 'Q', 'glu': 'E', 'gly': 'G', 'his': 'H', 'ile': 'I', 'leu': 'L',
         'lys': 'K', 'met': 'M', 'phe': 'F', 'pro': 'P', 'ser': 'S', 'thr': 'T', 'trp': 'W', 'tyr': 'Y', 'val': 'V'}
REACTIONS = {'MI:0213': 'writes', 'MI:0871': 'removes'}
# what each PRMT type makes (Bedford & Clarke, Mol Cell 33:1 (2009)), by gene as IntAct names it; other enzymes unchecked
PRMT_TYPE = {'PRMT1': 'I', 'PRMT2': 'I', 'PRMT3': 'I', 'CARM1': 'I', 'PRMT4': 'I', 'PRMT6': 'I', 'PRMT8': 'I', 'HMT1': 'I',
             'PRMT5': 'II', 'PRMT9': 'II', 'PRMT7': 'III'}
PRMT_MAKES = {'I': {'Rme', 'Rme1', 'Rme2', 'Rme2a'}, 'II': {'Rme', 'Rme1', 'Rme2', 'Rme2s'}, 'III': {'Rme', 'Rme1'}}
NOT_REACTION = {'MI:0407', 'MI:0915', 'MI:0914', 'MI:2364', 'MI:0403', 'MI:0408'}   # direct, physical, association, proximity, colocalisation, covalent
IA_KIND = {'protein': 'protein', 'peptide': 'protein', 'stable complex': 'protein', 'gene': 'gene', 'small molecule': 'small_molecule',
           'ds dna': 'dna', 'ss dna': 'dna', 'dna': 'dna', 'snrna': 'rna', 'rna': 'rna', 'mrna': 'rna', 'nucleic acid': 'rna'}
# the papers and methods behind the records that come from one study or model each
DATASETS = {
    'Kramer2014': ('Kramer et al., Nat Methods 11:1064 (2014); re-analysed by Sha, Kuster & Trendel, bioRxiv (2025)', '10.1038/nmeth.3092', 'UV crosslinking and mass spectrometry (RNA adducts)'),
    'XRNAX': ('Trendel et al., Cell 176:391 (2019); re-analysed by Sha, Kuster & Trendel, bioRxiv (2025)', '10.1016/j.cell.2018.11.004', 'UV crosslinking, XRNAX and mass spectrometry'),
    'RBS-ID': ('Bae et al., Nat Struct Mol Biol 27:678 (2020); re-analysed by Sha, Kuster & Trendel, bioRxiv (2025)', '10.1038/s41594-020-0436-2', 'UV crosslinking, chemical RNA digestion and mass spectrometry (RBS-ID)'),
    'pRBS-ID': ('Bae et al., Nat Commun 12:6026 (2021); re-analysed by Sha, Kuster & Trendel, bioRxiv (2025)', '10.1038/s41467-021-26317-5', 'photoactivatable ribonucleoside crosslinking and mass spectrometry (pRBS-ID)'),
    'pepR-MS': ('Sha, Kuster & Trendel, bioRxiv (2025)', '10.1101/2025.09.05.674461', 'UV crosslinking and mass spectrometry (pepR-MS compendium)'),
    'pepR-MS_focused': ('Sha, Kuster & Trendel, bioRxiv (2025)', '10.1101/2025.09.05.674461', 'UV crosslinking and mass spectrometry (pepR-MS focused sites)'),
    'CAPRI': ('Panhale et al., Nat Commun 10:2682 (2019)', '10.1038/s41467-019-10585-3', 'UV crosslinking and mass spectrometry (CAPRI), HEK293'),
    'pCLAP': ('Mullari et al., J Proteome Res 16:2762 (2017)', '10.1021/acs.jproteome.7b00042', 'UV crosslinking and peptide capture (pCLAP), HEK293'),
    'RBDmap': ('Castello et al., Mol Cell 63:696 (2016)', '10.1016/j.molcel.2016.06.029', 'UV crosslinking and proteolytic mapping (RBDmap)'),
    'RBR-ID': ('Oksuz et al., Mol Cell 83:2449 (2023)', '10.1016/j.molcel.2023.06.012', '4-thiouridine crosslinking and mass spectrometry (RBR-ID), K562 nuclei'),
    'Trendel2025': ('Trendel et al., Cell 188:4424 (2025)', '10.1016/j.cell.2025.04.037', 'DNA-protein crosslinking and mass spectrometry, MCF7'),
    'UniProt': ('UniProt Consortium, Nucleic Acids Res 53:D609 (2025)', '10.1093/nar/gkae1010', 'UniProt DNA-binding region (annotation, not a measured contact)'),
    'ProtVar': ('Burke et al., Nat Struct Mol Biol 30:216 (2023); Stephenson et al., Nucleic Acids Res 52:W140 (2024)', '10.1038/s41594-022-00910-8', 'AlphaFold-Multimer model of the pair (ProtVar)'),
    'AFDB': ('Han et al., bioRxiv (2026)', '10.64898/2026.03.27.714458', 'AlphaFold DB complex model'),
}

RECORD_COLUMNS = [
    ('site_id', 'string', '', 'the site (site table)'),
    ('record_class', 'code', 'interaction_class', 'what the record shows the partner doing at the site'),
    ('evidence_types', 'list<code>', 'evidence_type', 'the rows of supporting/site_interactions it comes from (pdb_methylated_contact: an entry modelling the site '
                                                      'methylated that PDBe-KB does not list for it), in code order'),
    ('partner_kind', 'code', 'record_partner', 'what the partner is'),
    ('partner_accession', 'string', '', 'its UniProtKB accession, of any species, empty when it has none'),
    ('partner_gene', 'string', '', 'its gene (human reviewed entries), or the name IntAct gives'),
    ('partner_name', 'string', '', 'its name in the experiment: the PDB entity description or the IntAct participant description'),
    ('partner_species', 'string', '', 'the partner\'s organism in the experiment'),
    ('source', 'string', '', 'PDB, IntAct, or the study or model (CAPRI, pCLAP, RBDmap, RBR-ID, Trendel2025, the pepR-MS datasets, UniProt, ProtVar, AFDB)'),
    ('source_id', 'string', '', 'the PDB entry, the IntAct interaction (EBI-), or the predicted model'),
    ('method', 'string', '', 'how it was found: the PDB experimental method and resolution, the IntAct detection method (PSI-MI), or the study\'s method'),
    ('context', 'string', '', 'the PDB entry title; for IntAct the interaction type and host organism; the localisation of a crosslink'),
    ('site_state', 'code', 'site_state', 'the site in that experiment'),
    ('state_detail', 'string', '', 'the residue modelled (chemical component, how the entity was matched, the copies), the PSI-MOD term or the mutation'),
    ('other_marks', 'list<string>', '', 'other modified residues of the site\'s protein in the experiment (K9me3, S10ph; UniProt numbering where the entity aligns)'),
    ('site_species', 'string', '', 'the organism of the site\'s protein in the experiment'),
    ('detail', 'string', '', 'the partner residues within 4 Å of a methylated site, the IntAct features of both participants, or the evidence detail of '
                             'supporting/site_interactions'),
    ('pubmed', 'string', '', 'PubMed id of the paper'),
    ('doi', 'string', '', 'its DOI'),
    ('citation', 'string', '', 'first author, journal and year (the primary citation of a PDB entry)'),
    ('paper_title', 'string', '', 'the title of the paper'),
    ('flags', 'list<code>', 'record_flag', 'what to check, in code order'),
]

SITE_COLUMNS = [
    ('interaction_classes', 'list<code>', 'interaction_class', 'the classes of the interaction records at the site (supporting/site_interaction_records), in code order'),
    ('partners_binding_methylated', 'list<string>', '', 'partners recorded binding the site methylated, by gene (the species in brackets when not human; '
                                                        'the name when there is no gene), sorted'),
    ('methylating_enzymes', 'list<string>', '', 'enzymes IntAct records methylating the site (methylation reaction), the same way'),
    ('demethylating_enzymes', 'list<string>', '', 'enzymes IntAct records demethylating the site (demethylation reaction), the same way'),
    ('interaction_papers', 'integer', '', 'distinct papers (PubMed or DOI) behind the experimental records at the site (predicted models left out)'),
]


def read(inputs, name: str) -> pd.DataFrame:
    return pd.read_csv(Path(inputs) / 'interactions' / f'{name}_{DATE}.tsv', sep='\t', dtype=str, keep_default_na=False)


def ordered(values, order) -> list:
    return sorted(set(values), key=lambda v: list(order).index(v))


def species(name: str) -> str:
    """The binomial of an organism name ('Homo sapiens HeLa ...' -> 'Homo sapiens'), as written otherwise."""
    m = re.match(r'([A-Z][a-z]+ [a-z]+)', name or '')
    return m.group(1) if m else (name or '')


def spans(ranges: str) -> list[tuple[int, int]]:
    out = []
    for part in ranges.split(','):
        a, _, b = part.strip().partition('-')
        if a.isdigit() and (b.isdigit() or not b):
            out.append((int(a), int(b or a)))
    return out


UNIPROT_AC = re.compile(r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}')


def base(acc: str) -> str:
    """A UniProt accession without its isoform or chain suffix (P62805-1 -> P62805); other identifiers (EBI-, CPX-,
    ENSG) as they are."""
    head = re.split(r'[-_ ]', acc or '')[0]
    return head if UNIPROT_AC.fullmatch(head) else (acc or '')


def gene_family(gene: str) -> str:
    """The letters of a gene symbol before its last number (JMJD6 -> JMJD), empty when it has no number."""
    m = re.fullmatch(r'([A-Z][A-Z0-9]*?[A-Z])(\d+)[A-Z]?', gene or '')
    return m.group(1) if m else ''


def title_names_paralogue(title: str, gene: str, genes: set) -> bool:
    """The title names another gene of the partner's family (JMJD6 for a JMJD4 record) and neither the partner, its
    number written out ('DNA methyltransferase 3A'), nor a range of the family holding it (MCM2-7)."""
    fam = gene_family(gene)
    if not fam or not title:
        return False
    up = title.upper()
    words = set(re.findall(r'[A-Z0-9]+', up))
    number = gene.upper()[len(fam):]
    if gene.upper() in words or number in words:
        return False
    for a, b in re.findall(re.escape(fam) + r'(\d+)-(\d+)', up):
        if int(a) <= int(re.match(r'\d+', number).group()) <= int(b):
            return False
    return any(w != gene.upper() and gene_family(w) == fam and w in genes for w in words)


def mark_text(item: str) -> str:
    """'K9 M3L' -> 'K9me3' for a mark of MARKS ('Kme3' where the entity does not align there); '' for anything else."""
    parts = item.split(' ')
    code = parts[-1]
    if code not in MARKS:
        return ''
    aa, suffix = MARKS[code][0], MARKS[code][1:]
    pos = re.sub(r'\D', '', parts[0]) if len(parts) == 2 else ''
    return f'{aa}{pos}{suffix}' if pos else f'{aa}{suffix} (not aligned)'


def mark_key(m: str):
    d = re.search(r'\d+', m)
    return int(d.group()) if d else 10 ** 6, m


def feature_items(features: str) -> list[tuple[str, str]]:
    """'nnme2arg:(arg-4) | sufficient to bind:(Region 2-22)' -> [(label, range text)]."""
    out = []
    for f in (features or '').split(' | '):
        lab, _, rng = f.partition(':')
        # 'lys-21' of 'n6me1lys:(lys-21)', and of the longer form 'n6me1lys:MI:0335 (certain):21..21(lys-21)'
        m = re.search(r'\(([^()]*)\)\s*$', rng)
        out.append((lab.strip(), (m.group(1) if m else rng).strip().strip(',')))
    return [x for x in out if x[0]]


def methyl_state(feature_type: str) -> str:
    """The state of a methylarginine feature type of scripts/37a (a PSI-MOD or PSI-MI name)."""
    f = feature_type.lower()
    for word, state in (('asymmetric', 'Rme2a'), ('symmetric', 'Rme2s'), ('dimethyl', 'Rme2'), ('monomethyl', 'Rme1')):
        if word in f:
            return state
    return 'Rme'


DEGREE = ['Rme', 'Rme1', 'Rme2', 'Rme2s', 'Rme2a']   # the most specific of several states a record gives


def residue_key(x: str):
    m = re.search(r'-?\d+$', x)
    return re.sub(r'-?\d+$', '', x), int(m.group()) if m else 0


def degree_flag(state: str, site_states: set) -> bool:
    """The experiment's methylation is not among the release's states at the site (site_states: the codes of
    methyl_states_reported, maron2021_s2_states, uniprot_states, states_by_study_design); False when either is unknown."""
    if state not in ('Rme1', 'Rme2', 'Rme2a', 'Rme2s') or not site_states:
        return False
    degrees = {1 if s == 'Rme1' else 2 for s in site_states if s in ('Rme1', 'Rme2', 'Rme2a', 'Rme2s')}
    if degrees and (1 if state == 'Rme1' else 2) not in degrees:
        return True
    types = {s for s in site_states if s in ('Rme2a', 'Rme2s')}
    return state in ('Rme2a', 'Rme2s') and bool(types) and state not in types


class Sources:
    """The staged provenance tables, indexed."""

    def __init__(self, inputs, genes: dict[str, str]):
        # the gene of a partner: the release's reference for human reviewed entries, else UniProt (other species)
        names = read(inputs, 'uniprot_partner_names')
        self.genes = {**{a: g or n for a, g, n in zip(names['accession'], names['gene'], names['protein_name'])}, **genes}
        self.entries = read(inputs, 'pdb_entries').set_index('pdb_id')
        ent = read(inputs, 'pdb_entities')
        self.entities = defaultdict(list)          # pdb -> rows
        for r in ent.itertuples():
            self.entities[r.pdb_id].append(r)
        st = read(inputs, 'pdb_site_residues')
        self.site_res = defaultdict(list)          # (pdb, acc, pos) -> rows
        for r in st.itertuples():
            self.site_res[(r.pdb_id, r.uniprot_accession, int(r.position))].append(r)
        self.methyl = read(inputs, 'pdb_methylated_sites')
        mc = read(inputs, 'pdb_methylated_site_contacts')
        self.contacts = defaultdict(list)          # (pdb, chain, entity position) -> rows
        for r in mc.itertuples():
            self.contacts[(r.pdb_id, r.chain, r.entity_position)].append(r)
        ia = read(inputs, 'intact_interactions')
        self.intact = defaultdict(list)
        for r in ia.to_dict('records'):
            self.intact[r['interaction_ac']].append(r)
        terms = read(inputs, 'intact_terms')
        self.mi = dict(zip(terms['obo_id'], terms['name']))
        self.label = dict(zip(terms['label'], terms['name']))
        pp = read(inputs, 'pubmed_papers')
        self.pubmed = {r.pubmed: r for r in pp.itertuples()}
        feat = lambda name: pd.read_csv(Path(inputs) / 'interactions' / f'intact_{name}_human_{INTACT_FEATURES}.tsv', sep='\t', dtype=str, keep_default_na=False)  # noqa: E731
        self.ia_methyl = defaultdict(list)          # (interaction ac, accession) -> [(start, end, feature type)]
        for r in feat('methylarginines').itertuples():
            for a, b in spans(r.ranges):
                self.ia_methyl[(r.interaction_ac, r.uniprot_accession)].append((a, b, r.feature_type))
        self.ia_mut = defaultdict(list)             # (interaction ac, accession) -> [(start, end, feature type, resulting residues)]
        for r in feat('mutations').itertuples():
            sp, res = spans(r.ranges), r.resulting_sequence.split(',')
            for k, (a, b) in enumerate(sp):
                # one resulting sequence per range where IntAct gives one each, else the whole
                self.ia_mut[(r.interaction_ac, r.uniprot_accession)].append((a, b, r.feature_type, res[k] if len(res) == len(sp) else r.resulting_sequence))

    def entity(self, pdb: str, ent_id: str):
        return next((e for e in self.entities.get(pdb, []) if e.entity_id == ent_id), None)

    def paper(self, pdb: str) -> dict:
        e = self.entries.loc[pdb] if pdb in self.entries.index else None
        if e is None:
            return {'pubmed': '', 'doi': '', 'citation': '', 'paper_title': ''}
        cite = ', '.join(x for x in (f'{e["first_author"]} et al.' if e['first_author'] else '', e['journal'], f'({e["year"]})' if e['year'] else '') if x)
        return {'pubmed': e['pubmed'], 'doi': e['doi'], 'citation': cite.replace(', (', ' ('), 'paper_title': e['citation_title']}

    def method(self, pdb: str) -> tuple[str, str]:
        if pdb not in self.entries.index:
            return '', ''
        e = self.entries.loc[pdb]
        meth = e['method'].lower().replace('x-ray diffraction', 'X-ray crystallography').replace('electron microscopy', 'cryo-EM').replace('solution nmr', 'solution NMR')
        return meth + (f', {e["resolution_a"]} Å' if e['resolution_a'] else ''), e['title']

    def site_state(self, pdb: str, acc: str, pos: int):
        """(state code, detail, other marks, species) of the site in an entry, from the entities aligned to it."""
        found = self.site_res.get((pdb, acc, pos), [])
        # an entity whose residues around the site differ from the human sequence places nothing there but an arginine
        rows = [r for r in found if r.context_check != 'differs' or r.residue in PDB_STATE]
        if not rows:
            return 'not_read', 'the residues around the site differ from the human sequence' if found else '', [], ''
        codes = Counter(r.residue for r in rows)
        meth = [c for c in codes if c in ('NMM', 'DA2', '2MR')]
        code = meth[0] if meth else 'ARG' if 'ARG' in codes else sorted(codes)[0]
        state = PDB_STATE.get(code, 'other_residue')
        detail = f'{code} modelled' + (f' (other copies: {", ".join(f"{c}" for c in sorted(codes) if c != code)})' if len(codes) > 1 else '') \
            + ('; the site found by its neighbours (the RCSB alignment is shifted, as for another isoform)' if any(r.context_check == 'realigned' and r.residue == code for r in rows) else '')
        marks = sorted({m for r in rows for x in r.other_modified_residues.split('; ') if x for m in [mark_text(x)] if m},
                       key=mark_key)
        orgs = {o for r in rows for e in [self.entity(pdb, r.entity_id)] if e is not None for o in e.organisms.split('; ') if o}
        return state, detail, marks, '; '.join(sorted(orgs))

    def partner_species(self, pdb: str, acc: str) -> tuple[str, str]:
        es = [e for e in self.entities.get(pdb, []) if acc and acc in e.uniprot_accessions.split('; ')]
        return '; '.join(sorted({o for e in es for o in e.organisms.split('; ') if o})), '; '.join(sorted({e.description for e in es if e.description}))

    def nucleic_names(self, pdb: str, kind: str) -> str:
        want = 'polyribonucleotide' if kind == 'rna' else 'polydeoxyribonucleotide'
        names = sorted({e.description for e in self.entities.get(pdb, []) if e.type == want and e.description})
        return '; '.join(names[:3]) + ('; …' if len(names) > 3 else '')


def pdb_records(ev: pd.DataFrame, src: Sources, site_states: dict) -> list[dict]:
    out = []
    for r in ev[ev['evidence_type'].isin(['pdb_protein', 'pdb_rna', 'pdb_dna', 'pdb_self'])].itertuples():
        acc, pos = r.site_id.split(':R')[0], int(r.site_id.split(':R')[1])
        for pdb in [x for x in r.structures.split('; ') if x]:
            state, detail, marks, site_org = src.site_state(pdb, acc, pos)
            meth, title = src.method(pdb)
            rec = {'site_id': r.site_id, 'evidence_types': r.evidence_type, 'source': 'PDB', 'source_id': pdb, 'method': meth, 'context': title,
                   'site_state': state, 'state_detail': detail, 'other_marks': '; '.join(marks), 'site_species': site_org, 'detail': '',
                   **src.paper(pdb)}
            if r.evidence_type == 'pdb_protein':
                org, name = src.partner_species(pdb, r.partner_accession)
                rec.update(partner_kind='protein', partner_accession=r.partner_accession, partner_gene=r.partner_gene or src.genes.get(r.partner_accession, ''), partner_name=name, partner_species=org,
                           record_class={'arginine': 'binds_unmodified', 'other_residue': 'binds_other_residue', 'not_read': 'binds_state_not_read'}.get(state, 'binds_methylated'))
            elif r.evidence_type == 'pdb_self':
                rec.update(partner_kind='self', partner_accession=acc, partner_gene='', partner_name='', partner_species=site_org, record_class='self')
            else:
                kind = 'rna' if r.evidence_type == 'pdb_rna' else 'dna'
                rec.update(partner_kind=kind, partner_accession='', partner_gene='', partner_name=src.nucleic_names(pdb, kind), partner_species='', record_class='nucleic_acid')
            out.append(rec)
    return out


def methylated_records(src: Sources, sites: set, have: set) -> tuple[list[dict], dict]:
    """Records of the entries modelling a site methylated: the contact rows join the PDBe-KB record of the same entry and
    partner (detail), and the others become records of their own (pdb_methylated_contact)."""
    extra, contact_detail = [], {}
    ms = src.methyl[(src.methyl['uniprot_accession'] + ':R' + src.methyl['position']).isin(sites)]
    for r in ms.itertuples():
        site = f'{r.uniprot_accession}:R{r.position}'
        meth, title = src.method(r.pdb_id)
        by_partner = defaultdict(lambda: {'chains': set(), 'residues': set(), 'd': None})
        for ch in [c for c in r.chains.split('; ') if c][:2]:
            for c in src.contacts.get((r.pdb_id, ch, r.entity_position), []):
                if not c.partner_chain or not c.partner_entity_id:
                    continue
                g = by_partner[c.partner_entity_id]
                g['chains'].add(c.partner_chain)
                g['residues'] |= {x for x in c.partner_residues.split('; ') if x}
                d = float(c.closest_distance_a) if c.closest_distance_a else None
                g['d'] = d if g['d'] is None else min(g['d'], d)
        site_entity = src.entity(r.pdb_id, r.entity_id)
        flags = ['matched_by_sequence'] if r.matched_by != 'UniProt alignment' else []
        if r.entity_organisms and HUMAN not in r.entity_organisms.split('; '):
            flags.append('site_protein_other_species')
        for ent_id, g in sorted(by_partner.items()):
            e = src.entity(r.pdb_id, ent_id)
            if e is None:
                continue
            accs = [a for a in e.uniprot_accessions.split('; ') if a]
            kind = 'rna' if e.type == 'polyribonucleotide' else 'dna' if e.type == 'polydeoxyribonucleotide' else 'protein'
            if kind == 'protein' and (r.uniprot_accession in accs or (site_entity is not None and e.description == site_entity.description)):
                kind = 'self'
            residues = sorted(g['residues'], key=residue_key)
            text = f'within 4 Å of the {r.residue} side chain ({g["d"]:.1f} Å): ' + ', '.join(residues[:12]) + (' …' if len(residues) > 12 else '')
            acc = accs[0] if accs else ''
            key = (site, r.pdb_id, acc if kind == 'protein' else kind)
            contact_detail[key] = text
            if (site, r.pdb_id, acc if kind == 'protein' else kind) in have:
                continue
            have.add(key)
            orgs = '; '.join(sorted({o for o in e.organisms.split('; ') if o}))
            fl = list(flags) + (['partner_other_species'] if kind == 'protein' and orgs and HUMAN not in orgs.split('; ') else [])
            st, _, marks, _ = src.site_state(r.pdb_id, r.uniprot_accession, int(r.position))
            extra.append({'site_id': site, 'record_class': {'protein': 'binds_methylated', 'self': 'self'}.get(kind, 'nucleic_acid'),
                          'evidence_types': 'pdb_methylated_contact', 'partner_kind': kind, 'partner_accession': acc,
                          'partner_gene': src.genes.get(acc, '') if kind == 'protein' else '', 'partner_name': e.description, 'partner_species': orgs,
                          'source': 'PDB', 'source_id': r.pdb_id, 'method': meth, 'context': title, 'site_state': PDB_STATE[r.residue],
                          'state_detail': f'{r.residue} modelled; ' + ('entity aligned to this protein' if r.matched_by == 'UniProt alignment'
                                                                       else f'{r.matched_by} ({r.matched_window})'),
                          'other_marks': '; '.join(marks), 'site_species': r.entity_organisms, 'detail': text, '_flags': fl, **src.paper(r.pdb_id)})
    return extra, contact_detail


def host_text(host: str) -> str:
    """IntAct's host organism as a phrase: 'in vitro', 'in Homo sapiens HeLa ... cell', '' when unknown."""
    h = (host or '').strip()
    if h.lower() in ('in vitro', 'in vivo'):
        return h.lower()
    return '' if h.lower() in ('', 'unknown', '-') else f'in {h}'


CHANGES_BINDING = re.compile(r'decreas|disrupt|increas|caus')   # the IntAct mutation effects that change an interaction


def ia_kind(b: dict, side: str) -> str:
    """What an IntAct participant is: a protein, a peptide or complex without a UniProt entry, a gene locus (Ensembl gene,
    or type gene), RNA, DNA, a small molecule."""
    pid, t = b[f'participant_{side}'], b[f'participant_type_{side}']
    if pid.startswith('ENSG'):
        return 'gene'
    k = IA_KIND.get(t, 'other')
    return 'peptide' if k == 'protein' and not UNIPROT_AC.fullmatch(base(pid)) else k


def row_partner_matches(r, b: dict, other: str, acc: str) -> bool:
    """The binary's other participant is the partner of the site_interactions row r."""
    o = base(b[f'participant_{other}'])
    if r.partner_kind == 'protein':
        return o == r.partner_accession
    if r.partner_kind == 'self':
        return o == acc
    # scripts/37 files every participant that is not a protein as RNA, DNA or other; any such participant matches,
    # the one of the same kind first (intact_records sorts them)
    return ia_kind(b, other) != 'protein'


def intact_records(ev: pd.DataFrame, src: Sources) -> list[dict]:
    """One record per site, IntAct interaction and partner of site_interactions: the binary record of the site's protein
    and that partner, or, where IntAct expanded a complex so that no binary names both (spoke expansion), the partner of
    the row with the interaction's method and paper."""
    by = {}   # (site, ac, partner) -> record; the strongest class when an interaction gives several rows
    order = list(CLASSES)
    for r in ev[ev['evidence_type'].str.startswith('intact_')].itertuples():
        acc, pos = r.site_id.split(':R')[0], int(r.site_id.split(':R')[1])
        pkey = r.partner_accession if r.partner_kind == 'protein' else r.partner_kind
        for ac in re.findall(r'EBI-\d+', r.detail):
            binaries = src.intact.get(ac, [])
            if not binaries:
                continue
            sides = [(me, other, b) for b in binaries for me, other in (('a', 'b'), ('b', 'a')) if base(b[f'participant_{me}']) == acc]
            exact = sorted([(me, other, b) for me, other, b in sides if row_partner_matches(r, b, other, acc)],
                           key=lambda x: ia_kind(x[2], x[1]) != r.partner_kind)
            if not exact and r.partner_kind == 'self':
                # scripts/37 files a partner it cannot name (a Complex Portal complex) as the protein itself: the
                # interaction's own other participant
                exact = [x for x in sides if base(x[2][f'participant_{x[1]}']) != acc][:1]
            if exact:
                me, other, b = exact[0]
            else:   # no binary names the site's protein and this partner: the interaction's method and paper only
                me, other, b = sides[0] if sides else (None, None, binaries[0])
            matched = bool(exact)
            kind = ia_kind(b, other) if matched else (r.partner_kind or 'other')
            if matched and base(b[f'participant_{other}']) == acc:
                kind = 'self'
            if kind == 'other' and matched:
                kind = ia_kind(b, other)
            mine = feature_items(b[f'features_{me}']) if me else []
            # the state of the site: the methylarginine or mutation features scripts/37a staged with their ranges
            state, sdetail = 'not_recorded', ''
            here = [t for a, z, t in src.ia_methyl.get((ac, acc), []) if a <= pos <= z]
            if here:
                found = sorted({methyl_state(t) for t in here}, key=DEGREE.index)
                state, sdetail = found[-1], '; '.join(sorted({re.sub(r'^psi-mod:"(MOD:\d+)"\((.*)\)$', r'\2 (\1)', t) for t in here}))
            else:
                muts = [(a, z, t, seq) for a, z, t, seq in src.ia_mut.get((ac, acc), []) if a <= pos <= z]
                if muts:
                    state, sdetail = 'mutated', '; '.join(sorted({f'{t}: ' + (f'R{pos} to {seq}' if a == z else f'residues {a}-{z} to {seq}')
                                                                  if seq and seq != '-' else t for a, z, t, seq in muts}))
            # the other marks of the site's participant, by IntAct's feature names (the curator's, whose numbering may be
            # the mature protein's); its methylarginines elsewhere by their staged ranges
            others = sorted({f'R{a}{methyl_state(t).replace("R", "", 1)}' for a, z, t in src.ia_methyl.get((ac, acc), []) if a == z and a != pos}, key=mark_key)
            marks = others + sorted({f'{rng} {LABEL_MARKS[lab]}' if re.search(r'\d', rng) else src.label.get(lab, lab)
                                     for lab, rng in mine if lab in LABEL_MARKS and lab not in METHYL_LABEL})
            itype = b['interaction_type_mi']
            role_me = b[f'biological_role_{me}'] if me else ''
            role_other = b[f'biological_role_{other}'] if matched else ''
            # the partner writes or removes the mark only where IntAct places the methylarginine at the site, or mutating
            # the site changes the reaction; a reaction elsewhere in the region (SETD7 on H3K4) is another reaction
            at_site = state.startswith('Rme') or (state == 'mutated' and bool(CHANGES_BINDING.search(sdetail)))
            effect = (re.match(r'^([a-z ]+),', r.detail) or [None, ''])[1] if r.evidence_type == 'intact_mutation' else ''
            if itype in REACTIONS and not (role_me == 'enzyme' and role_other != 'enzyme') and at_site:
                cls = REACTIONS[itype]
            elif itype and itype not in NOT_REACTION:
                cls = 'other_reaction'
            elif kind in ('rna', 'dna'):
                cls = 'nucleic_acid'
            elif kind in ('gene', 'small_molecule', 'other'):
                cls = 'not_a_protein'
            elif kind == 'self':
                cls = 'self'
            elif r.evidence_type == 'intact_methylated' or state.startswith('Rme'):
                cls = 'binds_methylated'
            elif r.evidence_type == 'intact_mutation':
                cls = 'mutation_no_effect' if effect == 'no effect' else 'mutation_unspecified' if effect in ('unspecified', '') else 'mutation_changes_binding'
            else:
                cls = 'binds_region'
            host = b['host_organism']
            taxid = b[f'participant_taxid_{other}'] if matched else ''
            psp = species(b[f'participant_species_{other}']) if matched and kind in ('protein', 'self') and not taxid.startswith('-') else ''
            fl = []
            if kind == 'protein' and taxid not in ('9606', '') and not taxid.startswith('-'):   # -1 in vitro, -2 chemical synthesis
                fl.append('partner_other_species')
            if me and b[f'participant_taxid_{me}'] not in ('9606', '') and not b[f'participant_taxid_{me}'].startswith('-'):
                fl.append('site_protein_other_species')
            if b['negative'] == 'true':
                fl.append('negative')
            cite = ', '.join(x for x in (b['first_author'], b['journal'].split(' (')[0], f'({b["year"]})' if b['year'] else '') if x).replace(', (', ' (')
            paper = src.pubmed.get(b['pubmed'])
            if matched:
                pacc = base(b[f'participant_{other}']) if kind in ('protein', 'self') else b[f'participant_{other}']
                pname = b[f'participant_description_{other}'] or b[f'participant_name_{other}']
                pgene = (r.partner_gene if r.partner_accession and r.partner_accession == pacc else '') or src.genes.get(pacc, '') \
                    or (b[f'participant_name_{other}'] if kind in ('protein', 'self') else '')
            else:
                pacc, pname, pgene = r.partner_accession, '', r.partner_gene or src.genes.get(r.partner_accession, '')
            rec = {'site_id': r.site_id, 'record_class': cls, 'evidence_types': r.evidence_type, 'partner_kind': kind,
                   'partner_accession': pacc, 'partner_gene': pgene, 'partner_name': pname, 'partner_species': psp,
                   'source': 'IntAct', 'source_id': ac, 'method': src.mi.get(b['detection_method_mi'], b['detection_method_mi']),
                   'context': ', '.join(x for x in (src.mi.get(itype, itype), host_text(host)) if x),
                   'site_state': state, 'state_detail': sdetail, 'other_marks': '; '.join(marks),
                   'site_species': species(b[f'participant_species_{me}']) if me else '',
                   'detail': ' / '.join(x for x in ((b[f'features_{me}'] if me else ''), (b[f'features_{other}'] if matched else '')) if x)
                             + ('' if matched else 'IntAct expanded a complex: no binary record names this partner with the site\'s protein'),
                   'pubmed': b['pubmed'], 'doi': b['doi'] or (paper.doi if paper else ''), 'citation': cite,
                   'paper_title': paper.title if paper else b['paper_title'], '_flags': fl}
            named = {int(x) for x in (paper.arginines_named.split('; ') if paper and paper.arginines_named else [])}
            if named and pos not in named and pos - 1 not in named and cls in ('writes', 'removes', 'binds_methylated', 'mutation_changes_binding',
                                                                              'mutation_no_effect', 'mutation_unspecified'):
                rec['_flags'].append('paper_names_other_arginine')
                rec['detail'] = (rec['detail'] + ' / ' if rec['detail'] else '') + 'the paper names arginine ' + ', '.join(map(str, sorted(named)))
            k = (r.site_id, ac, pkey)
            if k in by:
                old = by[k]
                old['evidence_types'] = '; '.join(sorted(set(old['evidence_types'].split('; ')) | {r.evidence_type}))
                if order.index(cls) < order.index(old['record_class']):
                    old['record_class'] = cls
                if old['site_state'] == 'not_recorded':
                    old['site_state'], old['state_detail'] = state, sdetail
            else:
                by[k] = rec
    return list(by.values())


def other_records(ev: pd.DataFrame) -> list[dict]:
    out = []
    for r in ev[~ev['evidence_type'].str.startswith(('pdb_', 'intact_'))].itertuples():
        cite, doi, method = DATASETS.get(r.source, ('', '', ''))
        pred = r.evidence_type.startswith('predicted')
        if pred:
            method += f', {"pDockQ" if r.source == "ProtVar" else "ipTM"} {r.score}' if r.score else ''
        out.append({'site_id': r.site_id, 'record_class': 'predicted' if pred else 'nucleic_acid', 'evidence_types': r.evidence_type,
                    'partner_kind': r.partner_kind or ('rna' if 'rna' in r.evidence_type else 'dna'), 'partner_accession': r.partner_accession,
                    'partner_gene': r.partner_gene, 'partner_name': '', 'partner_species': HUMAN if r.partner_accession else '', 'source': r.source,
                    'source_id': r.structures if pred else '', 'method': method, 'context': r.localisation.replace('_', ' '),
                    'site_state': 'not_recorded', 'state_detail': '', 'other_marks': '', 'site_species': HUMAN, 'detail': r.detail,
                    'pubmed': '', 'doi': doi, 'citation': cite, 'paper_title': '', '_flags': []})
    return out


def records(ev: pd.DataFrame, sites: pd.DataFrame, inputs, genes: dict[str, str]) -> pd.DataFrame:
    """supporting/site_interaction_records from supporting/site_interactions (ev) and the staged provenance."""
    src = Sources(inputs, genes)
    sites = sites[~sites['flags'].astype(str).str.contains('one_sided_window', regex=False)]   # the sites the analyses count
    states = {}
    for r in sites.itertuples():
        states[r.site_id] = {x for c in ('methyl_states_reported', 'maron2021_s2_states', 'uniprot_states', 'states_by_study_design')
                             for x in str(getattr(r, c, '') or '').split('; ') if x}
    rows = pdb_records(ev, src, states)
    for x in rows:
        x['_flags'] = []
    have = {(x['site_id'], x['source_id'], x['partner_accession'] if x['partner_kind'] == 'protein' else x['partner_kind']) for x in rows}
    extra, contact_detail = methylated_records(src, set(ev['site_id']) | set(sites['site_id']), have)
    for x in rows:
        k = (x['site_id'], x['source_id'], x['partner_accession'] if x['partner_kind'] == 'protein' else x['partner_kind'])
        if k in contact_detail and not x['detail']:
            x['detail'] = contact_detail[k]
        if x['partner_kind'] == 'protein' and x['partner_species'] and HUMAN not in x['partner_species'].split('; '):
            x['_flags'].append('partner_other_species')
        if x['site_species'] and HUMAN not in x['site_species'].split('; '):
            x['_flags'].append('site_protein_other_species')
    rows += extra + intact_records(ev, src) + other_records(ev)
    all_genes = set(genes.values())
    for x in rows:
        fl = x.pop('_flags')
        if degree_flag(x['site_state'], states.get(x['site_id'], set())):
            fl.append('state_not_reported')
        if x['source'] == 'IntAct' and x['partner_kind'] == 'protein' and title_names_paralogue(x['paper_title'], x['partner_gene'], all_genes):
            fl.append('title_names_paralogue')
        kind = PRMT_TYPE.get(x['partner_gene'].upper().replace('-', ''))
        if x['record_class'] == 'writes' and kind and x['site_state'].startswith('Rme') and x['site_state'] not in PRMT_MAKES[kind]:
            fl.append('state_unexpected_for_enzyme')
        if x['record_class'] != 'predicted' and not x['pubmed'] and not x['doi']:
            fl.append('no_publication')
        x['flags'] = '; '.join(ordered(fl, FLAGS))
    t = pd.DataFrame(rows, columns=[c for c, *_ in RECORD_COLUMNS])
    # the release's format: '; ' separates the items of a list column and appears in no other value
    for c, kind, *_ in RECORD_COLUMNS:
        clean = lambda v: re.sub(r'\s*;\s*', ', ', re.sub(r'[\t\r\n]+', ' ', str(v))).strip()  # noqa: E731
        t[c] = t[c].map(lambda v: '; '.join(x for x in (clean(i) for i in str(v).split('; ')) if x) if kind.startswith('list') else clean(v))
    t['_c'] = t['record_class'].map(list(CLASSES).index)
    t = t.sort_values(['site_id', '_c', 'partner_gene', 'partner_name', 'source', 'source_id'], kind='stable').drop(columns='_c').reset_index(drop=True)
    return t


def label(r) -> str:
    """A partner as the site columns name it: the gene, else the name; the species when not human."""
    name = r.partner_gene or r.partner_name or r.partner_accession
    names = re.split(r'[;,]\s*', r.partner_species) if r.partner_species else []
    sp = species(names[0]) if names else ''
    return name + (f' ({sp})' if sp and sp != HUMAN and HUMAN not in names else '')


def site_values(t: pd.DataFrame, site_ids) -> pd.DataFrame:
    """The five site columns from the records."""
    order = list(CLASSES)
    g = t.groupby('site_id')
    classes = g['record_class'].agg(lambda s: '; '.join(ordered(s, order)))

    def named(cls):
        sub = t[t['record_class'] == cls]
        return sub.groupby('site_id').apply(lambda d: '; '.join(sorted({label(r) for r in d.itertuples()})), include_groups=False) if len(sub) else pd.Series(dtype=str)
    meth, writes, removes = named('binds_methylated'), named('writes'), named('removes')
    exp = t[t['record_class'] != 'predicted']
    papers = exp.assign(p=exp['pubmed'].where(exp['pubmed'] != '', exp['doi'])).query('p != ""').groupby('site_id')['p'].nunique()
    idx = pd.Index(site_ids)
    return pd.DataFrame({'interaction_classes': idx.map(classes).fillna(''), 'partners_binding_methylated': idx.map(meth).fillna(''),
                         'methylating_enzymes': idx.map(writes).fillna(''), 'demethylating_enzymes': idx.map(removes).fillna(''),
                         'interaction_papers': idx.map(papers).fillna(0).astype(int)}, index=idx)


def changes(files: dict) -> list[str]:
    if 'site_interaction_records' not in files:
        return []
    t = files['site_interaction_records'][1]
    sites = files['sites'][1]
    n = Counter(t['record_class'])
    exp = t[t['record_class'] != 'predicted']
    fl = Counter(f for v in t['flags'] for f in v.split('; ') if f)
    meth = t[t['record_class'] == 'binds_methylated']
    new = t[t['evidence_types'] == 'pdb_methylated_contact']
    papers = exp['pubmed'].where(exp['pubmed'] != '', exp['doi'])
    out = [f'Every interaction record at a site now says where it comes from (supporting/site_interaction_records, {len(t):,} rows for '
           f'{t["site_id"].nunique():,} sites): the PDB entry, IntAct interaction, crosslinking study or predicted model, its method, the '
           f'species of both proteins, the paper ({papers[papers != ""].nunique():,} '
           f'papers behind the experimental records), the state of the site in the experiment and the other marks of the protein.',
           'What the records show: ' + '; '.join(f'{CLASS_WORDS[k]} {n[k]:,}' for k in CLASSES if n[k]) + '.',
           f'{meth["site_id"].nunique()} sites have a partner recorded binding them methylated ({len(meth):,} records from '
           f'{meth.loc[meth["source"] == "PDB", "source_id"].nunique()} PDB entries and {meth.loc[meth["source"] == "IntAct", "source_id"].nunique()} IntAct '
           f'interactions). {new["source_id"].nunique()} PDB entries model a site methylated and are not among those PDBe-KB lists for it: peptides '
           f'without a human UniProt mapping, orthologues identical around the site, and the other histone variants sharing a tail. Their '
           f'{len(new):,} records are new: ' + ', '.join(f'{CLASS_WORDS[k]} {v:,}' for k, v in Counter(new['record_class']).most_common()) + '.',
           'Flags for the reader to check (column flags): ' + '; '.join(f'{k} {fl[k]:,}' for k in FLAGS if fl[k]) + '.']
    if 'interaction_classes' in sites:
        w = sites['methylating_enzymes'].astype(str).str.len().gt(0).sum()
        out.append(f'Site columns interaction_classes, partners_binding_methylated, methylating_enzymes ({w} sites), demethylating_enzymes and '
                   f'interaction_papers summarise them.')
    return out


# ---- website: SITE_OUT/site_partners.json (and site_partners_all.json, the view with lower-confidence sites) -------------

def website_payload(t: pd.DataFrame, version: str) -> dict:
    """Records grouped per site: one row per class, partner, state, source and paper set, with its PDB entries or IntAct
    interactions; papers and methods by index."""
    papers, paper_ix, texts, text_ix = [], {}, [], {}

    def pid(r):
        if not (r.pubmed or r.doi or r.citation):
            return -1
        k = (r.pubmed, r.doi, r.citation)
        if k not in paper_ix:
            paper_ix[k] = len(papers)
            papers.append([r.pubmed, r.doi, r.citation, r.paper_title])
        return paper_ix[k]

    def tid(s):
        if s not in text_ix:
            text_ix[s] = len(texts)
            texts.append(s)
        return text_ix[s]
    classes, kinds, states, flags = list(CLASSES), list(PARTNER_KINDS), list(STATES), list(FLAGS)
    out = defaultdict(dict)
    sources = {}   # PDB entry or IntAct interaction -> [method, paper, context (title or type and host)], by index
    for r in t.itertuples():
        if r.source in ('PDB', 'IntAct') and r.source_id and r.source_id not in sources:
            sources[r.source_id] = [tid(r.method) if r.method else -1, pid(r), tid(r.context) if r.context else -1]
        # a record without an id (a crosslink, an annotation) is a row of its own
        key = (r.record_class, r.partner_kind, r.partner_accession, r.partner_gene, r.partner_species, r.site_state, r.source, r.flags,
               r.evidence_types, '' if r.source_id else r.Index)
        g = out[r.site_id].setdefault(key, {'ids': [], 'n': 0, 'papers': set(), 'methods': set(), 'marks': set(), 'details': [], 'names': set(),
                                            'state_detail': set(), 'context': set()})
        g['n'] += 1
        if r.source_id:
            g['ids'].append(r.source_id)
        p = pid(r)
        if p >= 0:
            g['papers'].add(p)
        if r.method:
            g['methods'].add(tid(r.method))
        g['marks'] |= {m for m in r.other_marks.split('; ') if m}
        if r.detail and len(g['details']) < 3:
            g['details'].append(tid(r.detail))
        if r.partner_name:
            g['names'].add(r.partner_name)
        if r.state_detail:
            g['state_detail'].add(r.state_detail)
        if r.context and r.source != 'PDB':
            g['context'].add(tid(r.context))
    sites = {}
    for s, groups in out.items():
        rows = []
        for (cls, kind, acc, gene, sp, state, source, fl, _, _), g in groups.items():
            rows.append([classes.index(cls), kinds.index(kind) if kind in kinds else len(kinds) - 1, acc, gene, '; '.join(sorted(g['names']))[:160], sp,
                         states.index(state), '; '.join(sorted(g['state_detail']))[:160], source, sorted(set(g['ids'])), sorted(g['methods']),
                         sorted(g['papers']), sorted(g['marks'], key=mark_key), g['details'], sorted(g['context']),
                         [flags.index(f) for f in fl.split('; ') if f], g['n']])
        sites[s] = rows
    return {'version': version, 'classes': [[k, v] for k, v in CLASSES.items()], 'kinds': kinds, 'states': [[k, v] for k, v in STATES.items()],
            'flags': [[k, v] for k, v in FLAGS.items()], 'papers': papers, 'texts': texts, 'sources': sources,
            'columns': ['class', 'partner_kind', 'accession', 'gene', 'name', 'species', 'state', 'state_detail', 'source', 'ids', 'methods',
                        'papers', 'other_marks', 'details', 'context', 'flags', 'n'], 'sites': sites}


def website(release: Path, version: str, inputs: Path, site_out: Path, name: str = 'site_partners.json') -> None:
    from release_layers import release_table, write_json
    path = Path(release) / 'supporting' / f'site_interaction_records_{version}.tsv'
    if not path.exists():
        return
    t = release_table(release, f'supporting/site_interaction_records_{version}.tsv')
    write_json(Path(site_out) / name, website_payload(t, version))


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(description='the website file of the interaction records of a release folder (or of its view with lower-confidence sites)')
    ap.add_argument('--release', required=True)
    ap.add_argument('--site-out', default='docs/data')
    ap.add_argument('--name', default='site_partners.json')
    args = ap.parse_args()
    version = json.loads((Path(args.release) / 'release.json').read_text())['version']
    website(Path(args.release), version, Path('inputs'), Path(args.site_out), args.name)
    print(f'{args.name} written')


if __name__ == '__main__':
    main()

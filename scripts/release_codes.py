"""Codes and column definitions of a release (used by scripts/build_release.py and scripts/release_docs.py).

Every coded column of a release table holds codes only, what each code means is here, and the release's
data_dictionary.tsv is written from this file. A column is (name, type, codes, description): type is one of integer,
number, string, boolean, code, list<code>, list<integer>, list<string>, codes names an entry of CODES (or is empty).
"""
from __future__ import annotations

import re

CONTEXT_WINDOW = 20  # residues: the boundary window of scripts/05 and scripts/06
PADJ_CUTOFF = 0.05   # Maron 2021 Table S3: a peptide changes when its adjusted p value is at most this (scripts/16)

# Studies and databases behind the union, in the order of the studies file: those of v1.1.0 (most sites first), then
# the studies added in v1.2.0. Citations and DOIs checked against Crossref.
STUDIES = [
    ('Larsen2016', 'Larsen 2016', 'Larsen SC et al., Sci Signal 9:rs9 (2016)', '10.1126/scisignal.aaf7329'),
    ('PSP', 'PhosphoSitePlus (PSP)', 'Hornbeck PV et al., Nucleic Acids Res 43:D512-D520 (2015)', '10.1093/nar/gku1267'),
    ('Maron2021', 'Maron 2021', 'Maron MI et al., iScience 24:102971 (2021)', '10.1016/j.isci.2021.102971'),
    ('Guo2014', 'Guo 2014', 'Guo A et al., Mol Cell Proteomics 13:372-387 (2014)', '10.1074/mcp.O113.027870'),
    ('Musiani2019', 'Musiani 2019', 'Musiani D et al., Sci Signal 12:eaat8388 (2019)', '10.1126/scisignal.aat8388'),
    ('Fong2019', 'Fong 2019', 'Fong JY et al., Cancer Cell 36:194-209 (2019)', '10.1016/j.ccell.2019.07.003'),
    ('ProMetheusDB', 'ProMetheusDB', 'Massignani E et al., Mol Cell Proteomics 21:100243 (2022)', '10.1016/j.mcpro.2022.100243'),
    ('Li2021', 'Li 2021', 'Li WJ et al., Nat Commun 12:1946 (2021)', '10.1038/s41467-021-21963-1'),
    ('Fedoriw2019', 'Fedoriw 2019', 'Fedoriw A et al., Cancer Cell 36:100-114 (2019)', '10.1016/j.ccell.2019.05.014'),
    ('UniProt', 'UniProt', 'UniProt Consortium, Nucleic Acids Res 53:D609-D617 (2025)', '10.1093/nar/gkae1010'),
    ('Lim2020', 'Lim 2020', 'Lim Y et al., Proteome Sci 18:6 (2020)', '10.1186/s12953-020-00162-8'),
    ('Wei2020', 'Wei 2020', 'Wei M et al., Front Cell Dev Biol 8:545934 (2020)', '10.3389/fcell.2020.545934'),
    ('IEDB', 'IEDB', 'Vita R et al., Nucleic Acids Res 47:D339-D343 (2019)', '10.1093/nar/gky1006'),
    ('Hartel2019', 'Hartel 2019', 'Hartel NG et al., Mol Cell Proteomics 18:2149-2164 (2019)', '10.1074/mcp.RA119.001625'),
    ('Geoghegan2015', 'Geoghegan 2015', 'Geoghegan V et al., Nat Commun 6:6758 (2015)', '10.1038/ncomms7758'),
    ('Sylvestersen2014', 'Sylvestersen 2014', 'Sylvestersen KB et al., Mol Cell Proteomics 13:2072-2088 (2014)', '10.1074/mcp.O113.032748'),
    ('Lu2025', 'Lu 2025', 'Lu L et al., EMBO Rep 26:5649-5672 (2025)', '10.1038/s44319-025-00599-y'),
    ('Ma2023', 'Ma 2023', 'Ma M et al., J Am Soc Mass Spectrom 34:1692-1700 (2023)', '10.1021/jasms.3c00154'),
]
# The studies read from their own tables (scripts/00_stage_study_sites.py, inputs/sources/<table>_methylarginine_records.tsv),
# in the order the build reads them (after the three compiled tables, so that record_ids stay those of v1.1.0).
STUDY_TABLES = {
    'larsen2016': 'Larsen2016', 'maron2021': 'Maron2021', 'guo2014': 'Guo2014', 'musiani2019': 'Musiani2019',
    'fong2019': 'Fong2019', 'li2021': 'Li2021', 'fedoriw2019': 'Fedoriw2019', 'uniprot': 'UniProt', 'lim2020': 'Lim2020',
    'wei2020': 'Wei2020', 'iedb': 'IEDB', 'hartel2019': 'Hartel2019', 'geoghegan2015': 'Geoghegan2015',
    'sylvestersen2014': 'Sylvestersen2014', 'lu2025': 'Lu2025', 'ma2023': 'Ma2023',
}
# The studies new in v1.2.0 (not in the compiled tables v1.1.0 read).
STUDIES_ADDED = ['Hartel2019', 'Geoghegan2015', 'Sylvestersen2014', 'Lu2025', 'Ma2023']


def study_of(source_table: str, source, reference) -> str:
    """study_id of a source record, from its source table and, in the compiled tables, its source and reference."""
    if source_table in STUDY_TABLES:
        return STUDY_TABLES[source_table]
    if source_table == 'human_ptm_master':
        return {'unip': 'UniProt', 'iedb': 'IEDB'}[source]
    if source_table == 'prometheus_mmc2':
        return 'ProMetheusDB'
    ref = '' if reference is None or reference != reference else str(reference).strip()
    if ref.lower().startswith('this study'):
        return 'Maron2021'
    if 'PhosphoSitePlus' in ref:
        return 'PSP'
    author = ref.split(' et al.')[0].split()[-1]
    return author + re.findall(r'\((\d{4})\)', ref)[-1]

STUDY_ORDER = [s[0] for s in STUDIES]
MARON_S5 = 'Maron 2021 Table S5 (CC BY-NC-ND 4.0, open-access supplement)'
PTM_MASTER = 'the iPTMnet-format PTM master table (iPTMnet: CC BY-NC-SA 4.0)'
STUDY_LICENCE = {
    'PSP': 'CC BY-NC-SA 3.0 Unported (PhosphoSitePlus terms of 2021, inputs/licences/), compiled in ' + MARON_S5,
    'ProMetheusDB': 'CC BY-NC-ND 4.0 (open-access supplement mmc2)',
    'UniProt': 'CC BY 4.0 (UniProt release 2026_03 MOD_RES)',
    'IEDB': 'IEDB query API, terms at www.iedb.org/terms_of_use_v3.php (epitopes the iPTMnet-format master table cites)',
    'Larsen2016': 'AAAS (Sci Signal supplementary Table S1: derived records only, the table is not redistributed)',
    'Musiani2019': 'PRIDE PXD009070 MaxQuant output (CC0), with the site list of ' + MARON_S5,
    'Maron2021': 'CC BY-NC-ND 4.0 (open-access supplement, Tables S2 and S3)',
    'Guo2014': 'CC BY 4.0 (supplementary Tables S2 and S3)',
    'Fong2019': 'Elsevier user licence (supplementary Table S4: derived records only, the table is not redistributed), PRIDE PXD012007 CC0',
    'Li2021': 'CC BY 4.0 (Supplementary Data 1, 3, 5 and 6)',
    'Fedoriw2019': 'Elsevier user licence (supplementary Tables S5 and S7: derived records only, the tables are not redistributed)',
    'Lim2020': 'CC BY 4.0 (Additional files 1 and 2)',
    'Wei2020': 'CC BY 4.0 (Table 1)',
    'Hartel2019': 'CC BY 4.0 (the article and its supplement, from 2020-09-15, and PRIDE PXD012357, CC0), records staged from supplementary Tables S1-S15',
    'Geoghegan2015': 'CC BY 4.0 (the article and its Supplementary Table 1)',
    'Sylvestersen2014': 'CC BY 4.0 (the article and its supplement, from 2020-10-12, and PRIDE PXD000559, CC0)',
    'Lu2025': 'CC BY 4.0 (the article and Dataset EV1)',
    'Ma2023': 'CC0 (PRIDE PXD038203 Table S1), the article is not open access',
}
STUDY_NOTE = {
    'Hartel2019': 'HEK293T and LN229 cells. ADMA and SDMA called from characteristic neutral losses (dimethylamine, methylamine) of dimethylarginine peptides (reported_methyl_state Rme2a, Rme2s), the other dimethylarginines as Rme2',
    'Geoghegan2015': 'Jurkat and primary human T cells, heavy-methyl SILAC with methylarginine antibodies',
    'Sylvestersen2014': 'HEK293T cells, monomethylarginine antibody enrichment',
    'Lu2025': 'HeLa cells, enrichment with the SMN Tudor domain',
    'Ma2023': 'breast cancer patient-derived xenografts grown in mice, ADMA antibody enrichment. A peptide assigned to a human protein can come from the mouse stroma where the two sequences agree',
    'PSP': 'downloaded 22 May 2021 (date from the Shechter lab), as compiled in Maron et al. 2021 Table S5',
    'Maron2021': "the study's own data: Tables S2 and S3, localised with the lab's Proteome Discoverer exports",
}

# What a record holds, how the study localised the methyl group and what was kept (scripts/00_stage_study_sites.py),
# and what to keep in mind. A record is then placed only where its peptide or window names exactly one arginine in the
# reviewed proteome (scripts/02 place_localised_records).
STUDY_QUALITY = {
    'Hartel2019': 'Peptide-spectrum matches and peptide groups (Proteome Discoverer) with the methylated arginine marked, '
                  'Percolator q <= 0.01 and ptmRS localisation >= 75 %. LN229 Tables S1-S2 give no positions: their rows '
                  'take the position Tables S3-S15 give for the same protein, peptide and residue, and are dropped without one.',
    'Geoghegan2015': 'Heavy-methyl SILAC peptide pairs with the methylated arginine marked, q <= 0.05 (0.93 % FDR for the paired '
                     'set), localised at ModLS > 8 or unambiguously (rows at ModLS 8 or below are left out). Positions number '
                     'the first accession of the protein group, often an unreviewed entry, so the site is placed from the peptide.',
    'Sylvestersen2014': 'MaxQuant site tables with a 13-residue window and the peptide, localisation probability >= 0.75. '
                        'The export writes every leucine as isoleucine: L is written back where the reported entry, or else '
                        'one site of the reviewed proteome, matches with I and L taken as one.',
    'Lu2025': 'The 743 dimethylarginine sites of Dataset EV1 come without peptides or localisation scores (SEQUEST HT, '
              'Percolator 1 % FDR). A site is kept only where a PSM localises it by mass, every arginine of its peptide '
              'dimethylated (155 sites), and the other 588 are left out.',
    'Ma2023': 'MaxQuant ADMA-antibody sites with a 31-residue window, localisation probability >= 0.75, 1 % FDR. Tumours '
              'grown in mice: a peptide identical in mouse and human can come from mouse stroma (flag xenograft_only).',
    # The studies compiled in Maron et al. 2021 Table S5, and the databases, now read from their own tables: the
    # study's methods (from the article, or for Musiani 2019, whose article could not be retrieved, from the PRIDE record
    # and deposited search), and what was kept. The compiled rows of a study read here add no site of their own.
    'Larsen2016': 'HEK293 cells, peptides fractionated at high pH and enriched with methylarginine antibodies, MaxQuant with '
                  'peptide, protein and site FDR 0.01, monomethylation of R and K searched. Read from the article\'s Table S1 '
                  '(8,900 sites, each with the leading protein, the position, the localisation probability and a 31-residue '
                  'window): arginine sites at localisation probability >= 0.75. The leading protein is often an isoform or '
                  'an unreviewed entry, so the site is placed from its window.',
    'PSP': 'Curated database, trusted as curated. Through Maron Table S5: the sites PhosphoSitePlus listed on 22 May 2021, '
           'each with a 15-residue window and position, without the study behind it, a localisation score or a degree. A '
           'record is placed where its window names exactly one arginine in the reviewed proteome, in an entry the record '
           'names.',
    'Maron2021': 'A549 cells treated with GSK591 or MS023, sequential Rme1, Rme2a and Rme2s antibody immunoprecipitation, '
                 'Byonic in Proteome Discoverer 2.4 with mono- and dimethylation of R and K searched, 1 % FDR. Tables S2 and '
                 'S3 give no localisation score; the lab\'s Proteome Discoverer peptide-group exports of the same '
                 'immunoprecipitations do (ptmRS site probability, a residue named only at >= 75 %). A mark of the tables is '
                 'kept where an export localises a methyl group to that arginine of the same peptide (with the export\'s '
                 'degree when it differs, noted), or where its peptide holds no other arginine or lysine. GluC peptides and '
                 'runs the exports do not hold are kept by mass only.',
    'Guo2014': 'HCT116 cells, monomethylarginine (Me-R4-100, D5A12) and ADMA antibodies, SEQUEST with a 1 % PSM FDR, mono- '
               'and dimethylation of R and K searched. Read from Tables S2 (peptide-spectrum matches with their ModScore) '
               'and S3 (sites with a 13-residue window): a site is kept where a PSM matching its window localises it at '
               'ModScore > 13, the paper\'s threshold, the degree read from the precursor mass. The positions number RefSeq '
               'proteins, so the record is placed from its peptide.',
    'Musiani2019': 'HeLa cells (heavy-methyl SILAC), SDMA then monomethylarginine antibody enrichment, MaxQuant 1.6.1.0 with '
                   'a 1 % FDR. Read from the MaxQuant site tables deposited in PRIDE PXD009070: arginine sites at '
                   'localisation probability >= 0.75 that the study reports (listed in Maron Table S5, the article\'s own '
                   'table could not be retrieved). The 268 Table S5 windows centred on a lysine are not among them.',
    'Fong2019': 'NB4 acute promyelocytic leukaemia cells, SDMA and monomethylarginine (and, with MS023, ADMA) antibody '
                'enrichment, MaxQuant with a 1 % FDR. Read from Table S4: arginine sites at localisation probability >= 0.75 '
                'and Andromeda score >= 25 (the paper accepted localisation >= 0.5), with the marked peptide and degree.',
    'Li2021': 'HEK293 cells (SILAC) after PRMT4, PRMT5 or PRMT7 knockdown or SGC3027 treatment, monomethylarginine then ADMA '
              'or SDMA antibody enrichment, Proteome Discoverer 2.2 with FDR < 1 %, methylation of R and K searched, no '
              'localisation score given. Read from Supplementary Data 1, 3, 5 and 6, the sites whose methylation fell at '
              'least two-fold or vanished: a site is kept where its peptide holds no other arginine or lysine.',
    'Fedoriw2019': 'DLBCL and pancreatic cancer cell lines with a type I PRMT inhibitor, monomethylarginine, ADMA and SDMA '
                   'antibody enrichment, SEQUEST with a 2.5 % FDR, no localisation score given. Read from Tables S5 and S7 '
                   '(gene symbol and marked peptide, the peptides whose methylation changed): a site is kept where its '
                   'peptide holds one arginine, or every arginine dimethylated, and occurs once in the entries of its gene, '
                   'which gives the position.',
    'Lim2020': 'Pooled colorectal cancer tissue from 10 patients, monomethylarginine and ADMA antibody enrichment, SEQUEST '
               'with a 2.5 % FDR. Read from Tables S1 and S2: sites with ModScore > 13 (the paper\'s threshold), each with '
               'peptide, window and degree. The tables write some leucines as I: L is written back where the reported '
               'entry, or one site of the reviewed proteome, matches with I and L taken as one residue.',
    'Wei2020': 'HPDE6c7 and PANC-1 cells, ADMA antibody enrichment, MaxQuant 1.6.0.16 with FDR < 0.01, Andromeda score > 40 '
               'and localisation probability > 0.75. Read from Table 1: the 289 sites with peptide and 31-residue window, '
               '137 of them numbered on an isoform, placed from the peptide.',
    'ProMetheusDB': 'Methyl-proteome table mmc2 of the heavy-methyl SILAC and SILAC experiments behind ProMetheusDB (SK-OV-3, '
                    'HeLa, NB4 and U2OS cells), MaxQuant with mono- and dimethylation of R and K and trimethylation of K '
                    'searched: marked peptides, 15-residue window and degree (me or di) of each site, no localisation score. '
                    'A record is kept where one of its peptides holds no other arginine or lysine, or holds no lysine and '
                    'has every arginine dimethylated, and is placed from those peptides.',
    'UniProt': 'UniProt release 2026_03 annotations of methylated arginines (MOD_RES) with experimental evidence from the '
               'literature (ECO:0000269), curated residue by residue. Annotations from large-scale proteomics '
               '(ECO:0007744), which repeat the studies read here, and those by similarity are left out. Position only, '
               'placed where the residue is an arginine of the entry.',
    'IEDB': 'Immune Epitope Database epitopes (HLA-presented peptides) with a methylated arginine, the 85 the '
            'iPTMnet-format master table cites, read from the IEDB query API: a site is kept where the epitope holds no '
            'other arginine or lysine. The degree is not given.',
}
# Placement codes of records not placed because their evidence names several sites, or one site in a protein the
# record does not name (the studies table counts each; the others not placed name none).
PLACEMENT_SEVERAL = {'unresolved_study_not_unique', 'unresolved_ambiguous_hit'}
PLACEMENT_OTHER_PROTEIN = {'unresolved_study_other_entry'}
# Placement codes of records not placed because their source gives no localisation of the methyl group.
PLACEMENT_NOT_LOCALISED = {'superseded_by_study_table', 'unresolved_not_localised'}
# Databases whose sites count as localised by curation (localised_by curated_database), placed from their window.
CURATED_DATABASES = {'PSP'}

SOURCE_TABLES = {
    'maron2021_s5': ('maron_s5', 'Maron et al. 2021 Table S5: methylarginine sites compiled from ten studies'),
    'prometheusdb_mmc2': ('prometheus_mmc2', 'ProMetheusDB (Massignani et al. 2022) supplementary table mmc2'),
    'ptm_master': ('human_ptm_master', 'the arginine methylation rows of the iPTMnet-format PTM master table (UniProt and IEDB records)'),
    'hartel2019': ('hartel2019', 'Hartel et al. 2019 supplementary Tables S1-S15, staged (scripts/00_stage_study_sites.py)'),
    'geoghegan2015': ('geoghegan2015', 'Geoghegan et al. 2015 Supplementary Table 1, staged (scripts/00_stage_study_sites.py)'),
    'sylvestersen2014': ('sylvestersen2014', 'Sylvestersen et al. 2014 supplementary tables, staged (scripts/00_stage_study_sites.py)'),
    'lu2025': ('lu2025', 'Lu et al. 2025 Dataset EV1, staged (scripts/00_stage_study_sites.py)'),
    'ma2023': ('ma2023', 'Ma et al. 2023 Table S1 (PRIDE PXD038203), staged (scripts/00_stage_study_sites.py)'),
    'larsen2016': ('larsen2016', 'Larsen et al. 2016 Table S1, staged (scripts/00_stage_study_sites.py)'),
    'maron2021': ('maron2021', 'Maron et al. 2021 Tables S2 and S3 with the lab\'s Proteome Discoverer exports, staged (scripts/00_stage_study_sites.py)'),
    'guo2014': ('guo2014', 'Guo et al. 2014 supplementary Tables S2 and S3, staged (scripts/00_stage_study_sites.py)'),
    'musiani2019': ('musiani2019', 'Musiani et al. 2019 MaxQuant site tables (PRIDE PXD009070), staged (scripts/00_stage_study_sites.py)'),
    'fong2019': ('fong2019', 'Fong et al. 2019 supplementary Table S4, staged (scripts/00_stage_study_sites.py)'),
    'li2021': ('li2021', 'Li et al. 2021 Supplementary Data 1, 3, 5 and 6, staged (scripts/00_stage_study_sites.py)'),
    'fedoriw2019': ('fedoriw2019', 'Fedoriw et al. 2019 supplementary Tables S5 and S7, staged (scripts/00_stage_study_sites.py)'),
    'uniprot': ('uniprot', 'UniProt release 2026_03 MOD_RES with experimental evidence, staged (scripts/00_stage_study_sites.py)'),
    'lim2020': ('lim2020', 'Lim et al. 2020 Additional files 1 and 2, staged (scripts/00_stage_study_sites.py)'),
    'wei2020': ('wei2020', 'Wei et al. 2020 Table 1, staged (scripts/00_stage_study_sites.py)'),
    'iedb': ('iedb', 'IEDB epitopes with a methylated arginine (query API), staged (scripts/00_stage_study_sites.py)'),
}
SOURCE_TABLE_CODE = {raw: code for code, (raw, _) in SOURCE_TABLES.items()}

CODES: dict[str, dict[str, str]] = {
    'methyl_state': {
        'Rme1': 'monomethylarginine',
        'Rme2': 'dimethylarginine, asymmetric or symmetric not stated',
        'Rme2a': 'asymmetric dimethylarginine',
        'Rme2s': 'symmetric dimethylarginine',
        'Rme': 'methylarginine, degree not stated (UniProt "Omega-N-methylated arginine")',
    },
    'state_summary': {
        'Rme1': 'monomethylarginine only',
        'Rme2': 'dimethylarginine only (asymmetric, symmetric or not stated)',
        'Rme1_Rme2': 'both mono- and dimethylarginine reported at the site',
    },
    'state_evidence': {
        'source_records': 'the source records (ProMetheusDB me/di, methyl_states_reported)',
        'maron2021_s2': 'Maron et al. 2021 Table S2 (maron2021_s2_states)',
        'uniprot': 'UniProt MOD_RES with experimental evidence (uniprot_states)',
        'study_design': 'a study that looked for one degree only (states_by_study_design)',
    },
    'dimethyl_type': {
        'ADMA': 'asymmetric dimethylarginine: the neutral loss of dimethylamine (45.058 Da)',
        'SDMA': 'symmetric dimethylarginine: the neutral loss of methylamine (31.042 Da)',
    },
    'antibody': {
        'MMA': 'anti-monomethylarginine immunoprecipitation',
        'ADMA': 'anti-asymmetric dimethylarginine immunoprecipitation',
        'SDMA': 'anti-symmetric dimethylarginine immunoprecipitation',
    },
    'inhibitor': {
        'GSK591': 'PRMT5 (type II) inhibitor',
        'MS023': 'type I PRMT inhibitor',
    },
    'maron_table': {
        'S2': 'Maron et al. 2021 Table S2: methylarginine residues from PTMScan immunoprecipitation, with the methylation degree',
        'S3': 'Maron et al. 2021 Table S3: methylpeptides from MMA, ADMA and SDMA immunoprecipitation after GSK591, MS023 or DMSO',
    },
    'call': {
        'decreased': 'lower with the treatment: at least one decreased peptide (or row) and none increased',
        'increased': 'higher with the treatment: at least one increased peptide (or row) and none decreased',
        'mixed': 'peptides (or rows) disagree',
        'no_change': 'no significant change',
    },
    'cell_line': {'MCF7': 'MCF7 breast cancer cells', 'MDA-MB-231': 'MDA-MB-231 breast cancer cells'},
    'motif': {
        'RG': 'G at +1',
        'RGG': 'G at +1 and +2',
        'GRG': 'G at -1 and +1',
        'GAR': 'G at -2 and A at -1',
        'RXR': 'R at -2 or +2',
        'CARM1_proline_rich': 'P within -5..+5',
        'CARM1_hydrophobic_proline': 'P within -5..+5 and a hydrophobic residue (A, I, L, M, F, W, V, Y) within -3..+3',
        'DR': 'D at -1',
        'RD': 'D at +1',
        'ER': 'E at -1',
        'RE': 'E at +1',
    },
    'disorder': {
        'disordered': 'inside a MobiDB predicted disordered region',
        'disorder_boundary': f'outside, within {CONTEXT_WINDOW} residues of a disordered region',
        'ordered': f'more than {CONTEXT_WINDOW} residues from any disordered region',
        'no_disorder_annotation': 'MobiDB has no disordered region on the protein',
    },
    'domain': {
        'in_domain': 'inside an InterPro domain or repeat',
        'boundary': f'outside, within {CONTEXT_WINDOW} residues of a domain or repeat',
        'inter_domain_linker': f'between two domains or repeats, more than {CONTEXT_WINDOW} residues from both',
        'distal': f'more than {CONTEXT_WINDOW} residues from any domain or repeat, not between two',
        'no_domain_annotation': 'InterPro has no domain or repeat on the protein',
    },
    'domain_side': {
        'inside': 'inside the nearest domain',
        'N_terminal_to_domain': 'before (N-terminal to) the nearest domain',
        'C_terminal_to_domain': 'after (C-terminal to) the nearest domain',
    },
    'interaction_status': {
        'experimental': 'the site has experimental interaction evidence: an RNA or DNA crosslink at the residue or as a candidate residue, an RBDmap fragment, or a PDB interface',
        'predicted_only': 'only predicted evidence: a predicted interface (ProtVar, pDockQ >= 0.23) or a confident predicted contact (a salt bridge or hydrogen bond, no clash, both residues at pLDDT >= 80, AlphaFold DB models at ipTM >= 0.8)',
        'none': 'none recorded. Most proteins were never crystallised or crosslinked, so this is not evidence of no interaction. IntAct records do not set the status',
    },
    'evidence_type': {
        'rna_crosslink': 'RNA-protein crosslink (pepR-MS re-analysis of Kramer 2014, XRNAX, RBS-ID, pRBS-ID, pepR-MS focused sites)',
        'rbdmap_fragment': 'inside an RNA-bound proteolytic fragment (RBDmap, Castello et al. 2016)',
        'dna_crosslink': 'DNA-protein crosslink (Trendel et al. 2025)',
        'dna_binding_region': 'inside a UniProt DNA-binding region (annotation, not evidence of contact)',
        'pdb_protein': 'at a protein-protein interface in the PDB (PDBe-KB)',
        'pdb_rna': 'at a protein-RNA interface in the PDB',
        'pdb_dna': 'at a protein-DNA interface in the PDB',
        'pdb_self': 'at an interface between copies of the same protein in the PDB',
        'predicted_interface': 'at a predicted interface residue (ProtVar, AlphaFold-Multimer)',
        'predicted_contact': 'a measured contact in a predicted complex model (ProtVar served model or AlphaFold DB complex)',
        'intact_mutation': 'an IntAct mutation feature at the residue that changes an interaction',
        'intact_binding_region': 'inside an IntAct binding region',
        'intact_methylated': 'IntAct records the residue as methylated',
    },
    'localisation': {
        'at_residue': 'at this residue',
        'candidate_residue': 'a crosslink whose candidate residues include this one',
        'within_3': 'a crosslinked residue within 3 residues (not this one)',
        'in_fragment': 'inside the fragment',
        'in_peptide': 'inside a crosslinked peptide, residue not localised',
        'in_region': 'inside the region',
    },
    'partner_kind': {
        'protein': 'another protein (partner_accession)',
        'self': 'another copy of the same protein',
        'rna': 'RNA',
        'dna': 'DNA',
        'other': 'a small molecule or another participant type (IntAct)',
    },
    'go_binding': {
        'experimental': 'a Gene Ontology annotation with an experimental evidence code',
        'other_evidence': 'annotated with other evidence codes only',
        'none': 'no annotation',
    },
    'curation': {
        'curated': 'Complex Portal curated complex',
        'predicted': 'hu.MAP 3.0 predicted complex (via Complex Portal)',
    },
    'partner_class': {
        'prmt': 'a protein arginine methyltransferase (PRMT1-9)',
        'tudor': 'has a Tudor domain (InterPro IPR002999)',
        'wd40': 'has WD40 repeats',
        'phd': 'has a PHD finger',
    },
    'flag': {
        'one_sided_window': "the site's only records are windows that match the entry exactly on one side of the reported arginine (at least 6 residues, the arginine included) and not on the other, and no other reviewed entry holds the window. Some come from another isoform of the protein. The pipeline leaves these records unplaced: the site is kept, and left out of the analyses, until the policy is decided",
        'paralog_window': 'a supporting window or peptide matches another reviewed entry exactly and this entry within 2 residues (open, accepted in config/remap_evidence_allowlist.tsv). The other entry is a paralog, except at RBM39 R109, where a 6-residue peptide differs from RBM39 only at an isobaric I/L and also occurs in CCDC13',
        'near_match_only': 'every supporting record matches the entry within 2 residues, not exactly (a variant or sequence conflict)',
        'arginine_from_position': 'the evidence of every supporting record covers other arginines of the entry too (a peptide with several arginines, or a window or peptide found more than once): the positions the sources report pick this one (evidence_check agrees_reported_arginine)',
        'shared_sequence': 'residues -15 to +15 occur identically in another reviewed entry (same_window_accessions): histone variants, identical paralogs. The evidence cannot tell them apart: the protein is the one the source named',
        'single_study': 'one study or database reports the site',
        'xenograft_only': 'only Ma et al. 2023 reports the site, from breast cancer patient-derived xenografts grown in mice: where the human and mouse sequences agree, the peptide can come from the mouse stroma',
    },
    'residue': {'R': 'arginine'},
    'records_change': {
        'same': 'the same records support the site in both versions',
        'records_added': 'every record that supported the site before supports it now, and records this version reads for the first time support it too',
        'changed': 'some records moved onto or off the site, others support it in both versions',
        'replaced': 'none of the records that supported the site before support it now, other records do',
    },
    'site_change': {
        'kept': 'a site of the previous version too',
        'added': 'not a site of the previous version',
        'removed': 'a site of the previous version only (site changes table)',
    },
    'record_change': {
        'same': 'placed on the same site in both versions',
        'moved': 'placed on another site than in the previous version',
        'newly_placed': 'placed now, not placed in the previous version',
        'no_longer_placed': 'placed in the previous version, not placed now',
        'not_placed': 'placed in neither version',
        'new_record': 'a source record the previous version did not read (source_correction accession_na)',
    },
    'reason': {
        'peptide_recentred': 'the previous version moved the peptide to the arginine nearest its middle. This version keeps the reported arginine inside the peptide',
        'unreviewed_accession': "the previous version filed reviewed-entry numbering under the source's unreviewed or obsolete accession. This version files it under the reviewed entry",
        'lysine_window': 'the source window is centred on a lysine. The previous version turned it into a nearby arginine site',
        'truncated_window': 'the window is cut short at a protein end. The previous version read it as if centred, taking the arginine at or nearest the middle of the short window. This version counts from its intact end',
        'multi_arginine_peptide': 'a peptide with several arginines that does not cover the reported position. The previous version guessed one',
        'peptide_start': "a peptide with several arginines, numbered differently from the entry. The previous version took the arginine nearest its middle. This version counts the reported position from the peptide's start in Maron 2021 Tables S2/S3",
        'wrong_protein': 'the previous version filed the position under another reviewed protein than the one whose sequence matched',
        'outdated_accession_or_symbol': 'the source accession is merged or obsolete, or its gene symbol outdated. This version reaches the current entry',
        'entry_not_in_reference': 'no reviewed entry has the source accession or gene symbol (an entry deleted from Swiss-Prot, a readthrough or an unreviewed entry). The evidence occurs once in the reviewed proteome, and this version places it there. For the deleted duplicate genes CBSL, OVOS2 and SIK1B that is the paralog holding the same sequence (CBS, OVOS1, SIK1)',
        'evidence_mismatch': 'the previous version kept the reported arginine although the source peptide or window does not match the entry there',
        'window_padding': 'the window is padded where it runs off a protein end. The padding now gives the cut side',
        'ambiguous_evidence': 'the evidence supports several sites equally. This version leaves the record unplaced',
        'one_sided_window': 'the site is kept pending a decision (flag one_sided_window)',
        'reference_renumbered': 'the same residue under a new number: UniProt 2026_03 extends the sequence at the N-terminus (RAB12, +96)',
        'new_record': 'a source record this version reads for the first time (a study table read for the first time) supports the site',
        'not_localised': 'the source gives nothing that localises the methyl group to one residue (a compiled record, or a study without localisation scores whose peptide leaves several residues): this version places only localised records',
        'not_unique_in_proteome': 'a localised record whose evidence names several sites in the reviewed proteome, one site in a protein the record does not name, or none: this version places such records only where the evidence is unequivocal',
        'other': 'placement changed for a reason not covered by the codes above',
    },
    'placement': {
        'sequence_unique_canonical': 'evidence sequence found once in the reported entry',
        'sequence_unique_candidate': 'evidence sequence found once, in another entry of the same gene or the current entry of a merged accession',
        'sequence_closest_supported': 'evidence sequence found more than once: the occurrence nearest the reported position',
        'sequence_position_verified': 'the reported arginine lies inside the peptide in the reference sequence',
        'sequence_unique_proteome': 'evidence sequence found once in the whole reference proteome (no candidate entry held it)',
        'sequence_peptide_start': "peptide with several arginines: the reported position counted from the peptide's start in Maron 2021 Tables S2/S3",
        'same_accession_position_verified': 'position-only record: the reported residue is an arginine of the entry',
        'canonical_position_fallback': 'position-only record: the reported residue is an arginine of the canonical entry',
        'same_accession_evidence_near_match': 'evidence not in the entry, but matches it at the reported arginine within 2 residues (variant or conflict)',
        'canonical_evidence_near_match': 'evidence not in the canonical entry, but matches it at the reported arginine within 2 residues',
        'unresolved_no_anchor': 'not placed: the evidence, if any, is in no candidate entry, and the reported residue is not an arginine',
        'unresolved_no_hit': 'not placed: evidence not found in any candidate entry',
        'unresolved_ambiguous_hit': 'not placed: evidence supports several sites equally',
        'unresolved_evidence_elsewhere': 'not placed: the evidence lies elsewhere in the entry, not at the reported position',
        'unresolved_evidence_mismatch': 'not placed: the evidence does not match the entry at the reported position',
        'unresolved_study_not_unique': 'not placed: a localised record whose evidence names more than one site in the reference proteome',
        'unresolved_study_not_found': 'not placed: a localised record whose evidence names no arginine anywhere in the reference proteome',
        'unresolved_study_other_entry': 'not placed: a localised record whose evidence names one site, in a protein the record does not name (the study may have seen the sequence in an isoform or fusion the reference proteome does not hold)',
        'superseded_by_study_table': "not placed: a compiled record (Maron 2021 Table S5, the iPTMnet-format master table) of a study read from its own table, whose records carry the study's localisation",
        'unresolved_not_localised': 'not placed: the source gives nothing that localises the methyl group to one residue (no localisation score, and a peptide that leaves several residues it could be on, or no peptide)',
    },
    'localised_by': {
        'localisation_score': "the study's localisation score passes its threshold (MaxQuant localisation probability >= 0.75, ptmRS >= 75 %, Ascore or ModScore > 13, ModLS > 8)",
        'one_methylatable_residue': 'the peptide holds one residue the search could methylate (no other arginine, and no lysine where lysine methylation was searched)',
        'every_arginine_dimethylated': 'every residue of the peptide the search could methylate is an arginine and every one is dimethylated, so the mass leaves no other arrangement',
        'curated_experimental': 'UniProt annotates the residue with experimental evidence from the literature (ECO:0000269)',
        'curated_database': 'PhosphoSitePlus lists the residue as a methylarginine site from its curation (trusted as curated, it gives no localisation score)',
    },
    'evidence_kind': {
        'window_13': 'a 13-residue window centred on the site',
        'window_15': 'a 15-residue window centred on the site',
        'window_31': 'a 31-residue window centred on the site',
        'peptide': 'the identified peptide',
        'marked_peptide': 'the peptide with the methylated residue marked',
        'position_only': 'no sequence: the reported position only',
    },
    'source_correction': {
        'nan_restored': 'Maron 2021 Table S5 gives this Larsen 2016 window without the tripeptide NAN (28 residues, where every other Larsen window has 31, no Larsen window holds NAN, so the letters were read as a missing value when the table was made). Restored where inserting NAN at exactly one place gives the window centred on the reported position of the reported entry',
        'accession_na': 'Table S5 gives no accession (NA), only a gene symbol and a window. The previous version left the record out. It is placed from its gene symbol and window, and numbered after the records of the previous version',
        'leucine_restored': "Sylvestersen et al. 2014's export writes every leucine as isoleucine (the two have the same mass). L is written back in the window and peptide where the reported entry, read at the reported position, agrees with them once I and L are taken as one residue. The evidence as given is in evidence_as_reported",
    },
    'evidence_check': {
        'agrees': 'placed on the one site its evidence names',
        'agrees_reported_arginine': 'placed at the reported arginine, which its peptide or window covers, the evidence covers other arginines of the entry too (a peptide with several arginines, or a window or peptide found more than once), so the reported position picks this one',
        'agrees_position_only': 'a record with a position and no sequence, placed at that arginine',
        'agrees_no_site': 'not placed, and its evidence names no site on a reviewed entry',
        'one_of_several': 'placed on one of several sites its evidence names',
        'ambiguous': 'not placed: its evidence names several sites',
        'dropped_evidence_has_site': 'not placed, although its evidence names one site',
        'contradicts': 'placed on a site its evidence does not name (accepted only for the rows of config/remap_evidence_allowlist.tsv, flag paralog_window)',
    },
}
CODES['source_table'] = {code: text for code, (_, text) in SOURCE_TABLES.items()}
CODES['study'] = {sid: f'{label}: {citation}' for sid, label, citation, _ in STUDIES}
VERDICT_CODE = {
    'agrees': 'agrees', 'agrees: reported arginine': 'agrees_reported_arginine', 'agrees: position only': 'agrees_position_only',
    'agrees: no site': 'agrees_no_site', 'one of several evidence sites': 'one_of_several',
    'evidence ambiguous': 'ambiguous', 'dropped, evidence has a site': 'dropped_evidence_has_site',
    'contradicts evidence': 'contradicts',
}
# Studies that looked for one degree of methylation only: their records report that degree by design.
STUDY_DESIGN_STATE = {'Larsen2016': 'Rme1', 'Wei2020': 'Rme2'}
FLAG_ORDER = list(CODES['flag'])
MOTIF_ORDER = list(CODES['motif'])
ANTIBODY_ORDER = list(CODES['antibody'])
EVIDENCE_ORDER = list(CODES['evidence_type'])
CALL_ORDER = list(CODES['call'])

# ---- Columns of each file -------------------------------------------------------------------------------------------

SITE_COLUMNS = [
    ('site_id', 'string', '', 'UniProt accession and arginine position, e.g. P35637:R216'),
    ('uniprot_accession', 'string', '', 'reviewed UniProtKB accession of the canonical sequence'),
    ('uniprot_entry_name', 'string', '', 'UniProtKB entry name'),
    ('gene', 'list<string>', '', 'primary gene name of the UniProt entry (several for a few entries, e.g. histone genes), empty when UniProt gives none'),
    ('protein_name', 'string', '', 'recommended protein name of the UniProt entry'),
    ('position', 'integer', '', 'position of the arginine in the reference sequence (1-based, initiator methionine = 1)'),
    ('residue', 'code', 'residue', 'residue at that position in the reference sequence: always R'),
    ('window_15', 'string', '', 'residues -7 to +7 around the site, _ beyond a protein end'),
    ('window_31', 'string', '', 'residues -15 to +15 around the site, _ beyond a protein end'),
    ('protein_length', 'integer', '', 'length of the reference sequence'),
    ('reference_sequence_md5', 'string', '', 'MD5 of the reference sequence, uppercase hex as UniProt gives it (a later UniProt entry with the same MD5 has the same numbering)'),
    ('methyl_states_reported', 'list<code>', 'methyl_state', 'methylation degrees given by the source records, in code order, empty where they give none'),
    ('maron2021_s2_states', 'list<code>', 'methyl_state', 'degrees called at the site in Maron et al. 2021 Table S2'),
    ('uniprot_states', 'list<code>', 'methyl_state', 'degrees UniProt release 2026_03 annotates at the site (MOD_RES) with experimental evidence: ECO:0000269, or ECO:0007744 (large-scale proteomics). ADMA and SDMA calls in UniProt mostly follow the antibody of the experiment'),
    ('uniprot_states_inferred', 'list<code>', 'methyl_state', 'degrees UniProt annotates at the site only by similarity to another species (ECO:0000250) or without experimental evidence; not counted in methyl_state'),
    ('states_by_study_design', 'list<code>', 'methyl_state', 'the one degree a study looked for: Rme1 for Larsen 2016 (monomethylarginine antibodies; Table S1 lists monomethylarginine sites, and the deposited search, PRIDE PXD003700, has monomethylation only), Rme2 for Wei 2020 (ADMA antibody, PTMScan #13474; only dimethylation searched)'),
    ('methyl_state', 'code', 'state_summary', 'the degrees reported at the site by any of methyl_states_reported, maron2021_s2_states, uniprot_states and states_by_study_design; empty when none gives one. Mass spectrometry tells mono- from dimethylation, but asymmetric from symmetric mostly by the antibody used'),
    ('methyl_state_evidence', 'list<code>', 'state_evidence', 'which of these give a degree, in code order'),
    ('maron2021_s3_antibodies', 'list<code>', 'antibody', 'antibodies whose immunoprecipitation recovered the site in Maron et al. 2021 Table S3 (supporting/maron2021_detection)'),
    ('dimethyl_type_neutral_loss', 'list<code>', 'dimethyl_type', 'the dimethylarginine type Hartel et al. 2019 called at the site from neutral losses in the fragment spectra, independent of the antibody used (both when both were called). Shown beside maron2021_s3_antibodies, which says which antibodies recovered the site; the two can differ'),
    ('maron2021_GSK591_call', 'code', 'call', f'response to GSK591 in Maron et al. 2021 Table S3, combined over antibodies (supporting/maron2021_inhibitors: a peptide changes at adjusted p <= {PADJ_CUTOFF}), empty when not measured'),
    ('maron2021_MS023_call', 'code', 'call', 'response to MS023, as for GSK591'),
    ('carm1_ko_call', 'code', 'call', 'response to CARM1 knockout (Shishkova et al. 2017), combined over cell lines (supporting/carm1_ko), empty when not measured'),
    ('motifs', 'list<code>', 'motif', 'sequence motifs around the site (scripts/08), in code order'),
    ('disorder', 'code', 'disorder', 'MobiDB predicted disorder at the site (scripts/06)'),
    ('disorder_edge_distance', 'integer', '', 'residues to the nearest edge of a disordered region (a protein end is not an edge), empty with no annotation or when the region spans the whole protein'),
    ('domain', 'code', 'domain', 'InterPro 108.0 domain or repeat context (scripts/05)'),
    ('nearest_domain_interpro', 'string', '', 'InterPro accession of the nearest domain or repeat (the innermost when inside several)'),
    ('nearest_domain_name', 'string', '', 'its InterPro name'),
    ('nearest_domain_start', 'integer', '', 'its first residue'),
    ('nearest_domain_end', 'integer', '', 'its last residue'),
    ('nearest_domain_class', 'string', '', 'its class in config/domain_class_ontology.json, as scripts/07 assigns it, No domain annotation when none'),
    ('domain_side', 'code', 'domain_side', 'where the site lies relative to the nearest domain, empty when none'),
    ('domain_edge_distance', 'integer', '', 'residues to the nearest domain edge (0 inside a domain), empty when none'),
    ('interaction_status', 'code', 'interaction_status', 'interaction evidence at the site (supporting/site_interactions), empty for the sites kept pending'),
    ('interaction_evidence', 'list<code>', 'evidence_type', 'the evidence types at the site, in code order'),
    ('n_studies', 'integer', '', 'number of studies or databases reporting the site'),
    ('studies', 'list<code>', 'study', 'the studies or databases reporting the site, in the order of the studies file, PhosphoSitePlus is PSP, as its terms ask'),
    ('study_dois', 'list<string>', '', 'their DOIs, in the same order'),
    ('source_tables', 'list<code>', 'source_table', 'the source tables of its records, in code order'),
    ('n_records', 'integer', '', 'number of source records placed on the site'),
    ('record_ids', 'list<integer>', '', 'the records behind the site (record_id in supporting/source_records), ascending'),
    ('reported_positions', 'list<integer>', '', 'positions the records report, ascending, only when one differs from position'),
    ('flags', 'list<code>', 'flag', 'cautions, in code order'),
    ('same_window_accessions', 'list<string>', '', 'every other reviewed entry holding the 31 residues of window_31, sorted'),
    ('same_window_genes', 'list<string>', '', 'the first primary gene of each (the accession when UniProt gives no gene), same order'),
    ('in_previous_version', 'boolean', '', 'whether the previous version had this site'),
    ('change_since_previous_version', 'code', 'site_change', 'kept or added (removed sites are in supporting/site_changes)'),
    ('change_reasons', 'list<code>', 'reason', 'why the site was added, or renumbered (CHANGES.md)'),
]

STUDY_COLUMNS = [
    ('study_id', 'code', 'study', 'study or database'),
    ('label', 'string', '', 'name to show'),
    ('source_table', 'list<code>', 'source_table', 'the source tables its records come from'),
    ('citation', 'string', '', 'citation, checked against Crossref'),
    ('doi', 'string', '', 'DOI'),
    ('licence', 'string', '', 'terms of the records as this release received them'),
    ('n_records', 'integer', '', 'source records'),
    ('n_records_placed', 'integer', '', 'of these, placed on a site'),
    ('n_sites', 'integer', '', 'sites (status site) it reports'),
    ('n_sites_only_here', 'integer', '', 'sites no other study or database reports'),
    ('n_sites_in_previous_version', 'integer', '', 'of its sites, those the previous version had'),
    ('n_records_several_sites', 'integer', '', 'records not placed because their evidence names several sites'),
    ('n_records_other_protein', 'integer', '', 'records not placed because their evidence names one site, in a protein the record does not name'),
    ('n_records_not_localised', 'integer', '', 'records not placed because the source gives no localisation of the methyl group (placement superseded_by_study_table or unresolved_not_localised)'),
    ('n_records_no_site', 'integer', '', 'records not placed for another reason (placement_code in supporting/source_records)'),
    ('note', 'string', '', ''),
    ('quality', 'string', '', 'what a record holds, how the study localised the methyl group, and what was kept'),
]

RECORD_COLUMNS = [
    ('record_id', 'integer', '', 'row number of the record in the source tables as the pipeline reads them, the same in every release'),
    ('source_table', 'code', 'source_table', 'the source table the record comes from'),
    ('study', 'code', 'study', 'the study or database that reported it'),
    ('study_doi', 'string', '', 'its DOI'),
    ('reference_text', 'string', '', 'the reference as the source table gives it'),
    ('pmid', 'string', '', 'PubMed ID, where the source gives one'),
    ('reported_accession', 'list<string>', '', 'accession as reported (a few records give two)'),
    ('reported_gene', 'string', '', 'gene as reported'),
    ('reported_position', 'string', '', 'position as reported'),
    ('reported_methyl_state', 'string', '', 'methylation state as reported, in the codes of the site table (ProMetheusDB me = Rme1, di = Rme2)'),
    ('reported_enzyme', 'string', '', 'methyltransferase as reported (PTM master table only)'),
    ('reported_sample', 'list<string>', '', 'cells or tissue the study measured (records read from a study table)'),
    ('reported_enrichment', 'list<string>', '', 'how the study enriched the peptides, e.g. ADMA IAP (antibody), SCX fractionation, SMN Tudor domain (records read from a study table)'),
    ('reported_detail', 'string', '', 'what else the study gives for the record: the neutral-loss call, the localisation probability or score, the source rows merged into the record and their tables, the PRMT1 knockdown result (records read from a study table)'),
    ('evidence_kind', 'code', 'evidence_kind', 'kind of sequence evidence'),
    ('evidence_sequence', 'list<string>', '', 'the window, peptide or marked peptide (a few records give several), as the pipeline reads it: modifications other than methylation, e.g. (ox), are left out, and source_correction says when it was corrected'),
    ('evidence_window', 'list<string>', '', 'the sequence window a marked-peptide record also carries (ProMetheusDB SeqWindow), which the placement uses as well'),
    ('source_correction', 'code', 'source_correction', 'a fault of the source row this version corrects before placing the record'),
    ('evidence_as_reported', 'string', '', 'the evidence as the source gives it, when source_correction changed it'),
    ('localised_by', 'code', 'localised_by', 'what localises the methyl group to one residue, empty when nothing does (the record is then not placed)'),
    ('placed_site', 'string', '', 'the site this version places the record on, empty when not placed'),
    ('placement_code', 'code', 'placement', 'how the record was placed (scripts/02 remap_status)'),
    ('evidence_check', 'code', 'evidence_check', "independent check of the placement against the record's own evidence (scripts/qc_remap_evidence_check.py)"),
    ('accepted_exception', 'string', '', 'why a contradiction is accepted for now (config/remap_evidence_allowlist.tsv)'),
    ('pending_site', 'string', '', 'the site a one-sided window record supports (flag one_sided_window), empty otherwise'),
    ('previous_site', 'string', '', 'the site the previous version placed the record on, in its numbering'),
    ('change_since_previous_version', 'code', 'record_change', 'how the placement changed'),
    ('reason_code', 'code', 'reason', 'why it changed'),
]

CHANGE_COLUMNS = [
    ('site_id', 'string', '', 'the site: in this version\'s numbering when the site is in this version'),
    ('previous_site_id', 'string', '', 'the same site in the previous version\'s numbering, empty for an added site'),
    ('gene', 'list<string>', '', 'primary gene name (for an entry not in the reference, the gene its records report)'),
    ('uniprot_accession', 'string', '', 'accession'),
    ('position', 'integer', '', 'arginine position, in this version\'s numbering when the site is in this version'),
    ('in_previous_version', 'boolean', '', 'a site of the previous version'),
    ('in_this_version', 'boolean', '', 'a site of this version'),
    ('change', 'code', 'site_change', 'kept, removed or added'),
    ('reason_codes', 'list<code>', 'reason', 'why it was removed, added or renumbered, or, for a kept site, why records moved onto or off it'),
    ('records_change', 'code', 'records_change', 'for a kept site, whether its records changed'),
    ('records_now_at_or_previously_at', 'list<string>', '', 'for a removed site, where its records are placed now; for an added site, where they were placed before; for a kept site, both, for the records that moved'),
    ('record_ids', 'list<integer>', '', 'the records involved (supporting/source_records): those placed on a removed site before, or on an added site now, and for a kept site those that moved onto or off it'),
]

MARON_DETECTION_COLUMNS = [
    ('site_id', 'string', '', 'the site'),
    ('maron_table', 'code', 'maron_table', 'the table'),
    ('antibody', 'code', 'antibody', 'the immunoprecipitation (Table S3 only)'),
    ('methyl_states', 'list<code>', 'methyl_state', 'the degrees called (Table S2 only)'),
    ('n_peptides', 'integer', '', 'distinct peptides (sequence and modifications) carrying the site'),
]

MARON_INHIBITOR_COLUMNS = [
    ('site_id', 'string', '', 'the site'),
    ('antibody', 'code', 'antibody', 'the immunoprecipitation'),
    ('inhibitor', 'code', 'inhibitor', 'the treatment, compared with DMSO'),
    ('n_peptides', 'integer', '', 'peptide rows carrying the site with a measured difference'),
    ('n_decreased', 'integer', '', f'of these, lower with the inhibitor at adjusted p <= {PADJ_CUTOFF}'),
    ('n_increased', 'integer', '', f'higher with the inhibitor at adjusted p <= {PADJ_CUTOFF}'),
    ('log2_diff_median', 'number', '', 'median log2 intensity difference, inhibitor minus DMSO'),
    ('log2_diff_strongest', 'number', '', 'the difference of the peptide with the lowest adjusted p'),
    ('adj_p_min', 'number', '', 'the lowest adjusted p value (Maron et al. 2021)'),
    ('call', 'code', 'call', 'decreased, increased, mixed or no_change over the peptides, peptides with several arginines give each the same values'),
]

CARM1_COLUMNS = [
    ('site_id', 'string', '', 'the site'),
    ('cell_line', 'code', 'cell_line', 'the cell line'),
    ('n_peptides', 'integer', '', 'peptides carrying the site'),
    ('log2_ko_over_wt', 'number', '', 'median log2 knockout / parental intensity over the peptides'),
    ('p', 'number', '', 'the lowest p value over the peptides (Shishkova et al. 2017)'),
    ('n_decreased', 'integer', '', 'peptides at log2 <= -1 and p < 0.05'),
    ('n_increased', 'integer', '', 'peptides at log2 >= 1 and p < 0.05'),
    ('call', 'code', 'call', 'a peptide decreases at log2 <= -1 and p < 0.05, increases at log2 >= 1 and p < 0.05 (scripts/41), combined over peptides'),
]

CODES['analysis'] = {
    'panther_group': 'methylated proteins in a group of PANTHER 19.0 top-level protein classes (config/panther_class_groups.json, scripts/03b), against the other reviewed human proteins. The groups and No PANTHER protein class partition the reviewed proteome: a + c sums to every reviewed entry',
    'panther_class': 'the same for the 24 PANTHER top-level protein classes and No PANTHER protein class (again a partition)',
    'panther_subclass': 'the same for every PANTHER class below the top level with at least 10 reviewed proteins (a protein counts in a class when PANTHER gives it that class or one under it)',
    'functional_class': 'methylated proteins in a broad functional class (config/functional_ontology.json, scripts/03: the paper\'s keyword classes, which overlap), against the other reviewed human proteins',
    'functional_subclass': 'the same for functional subclasses',
    'motif': 'methylarginines with a sequence motif (scripts/08), against the other arginines of the methylated proteins',
    'disorder': 'methylarginines in a disorder class (scripts/06), against the other arginines of the methylated proteins',
    'domain_context': 'methylarginines in a domain context (scripts/05, 07), against the other arginines of the methylated proteins. A protein whose every arginine is a site has no background and is left out (MIF, 5 sites)',
    'domain_class': 'methylarginines whose nearest domain is of a class (scripts/07), against the other arginines of the methylated proteins (MIF left out, as for domain_context)',
    'disease_family': 'methylarginine-containing proteins (counted by gene) in one of the eight audited disease families (scripts/14), against the other reviewed human genes with a DisGeNET term (DOSE\'s background for the terms), q over the eight',
    'disease_family_additional': 'the same for the additional disease families (scripts/disease_families.py ADDITIONAL_FAMILY_RULES: rule-based, not audited, applied only to the terms no audited family takes), q over these families',
    'condensate': 'methylated proteins among the members of a CD-CODE condensate (scripts/12), against the other condensate proteins',
    'hpa_location': 'methylated proteins at a Human Protein Atlas location, against the other proteins HPA locates (Approved or better)',
    'pdb_interface': 'methylarginines at a PDB interface of a kind (category: rna, dna, protein including self, any), against the other arginines of the same proteins that a PDB structure covers (scripts/37). Strata: all, ordered (outside MobiDB disorder), domain (inside an InterPro domain)',
    'predicted_interface': 'methylarginines at a predicted interface residue (ProtVar, category the pDockQ floor), against the other arginines of the proteins with a predicted interface. Strata: all, ordered, disordered',
    'crosslink': 'methylarginines with RNA or DNA crosslink evidence of a kind (category), against the other arginines of the proteins with that evidence, the _outside_rg categories leave out arginines in an RG or GR dipeptide',
    'intact': 'IntAct: methylarginines whose mutation changes an interaction (against other mutated arginines), and methylarginines inside a binding region of at most 20 residues (against the other arginines of proteins with one)',
    'go_binding': 'methylated proteins with a Gene Ontology RNA- or DNA-binding annotation (category, _experimental: experimental evidence codes only), against all reviewed human proteins (clusterProfiler enricher, the four annotations as one job)',
}
CODES['method'] = {
    'enricher': 'over-representation analysis with clusterProfiler::enricher (scripts/enrich_sets.R, versions in release.json): for a set with M members in a background of N, of which n are methylated (the targets), p is the one-sided hypergeometric probability of k or more methylated members, and fold_enrichment = (k/n)/(M/N). a = k, b = n - k (methylated proteins outside the set), c = M - k (other members of the set), d = N - M - n + k. q = Benjamini-Hochberg over the sets of the analysis (clusterProfiler p.adjust). A set is tested when it has at least one methylated member and its size is within the analysis\'s range',
    'fisher': 'two-sided Fisher exact test (scipy), q = Benjamini-Hochberg within the analysis and method. The odds ratio adds 0.5 to each cell (Haldane-Anscombe)',
    'mantel_haenszel': 'Mantel-Haenszel common odds ratio with each protein a stratum, and the continuity-corrected Cochran-Mantel-Haenszel test (statsmodels), a b c d are the pooled counts. q = Benjamini-Hochberg within the analysis and method',
}
CODES['ptm_type'] = {'phospho': 'phosphorylation (S, T, Y)', 'lys_acyl': 'lysine acylation', 'ub_sumo': 'ubiquitination or SUMOylation (K)', 'lys_methyl': 'lysine methylation'}

PROTEIN_COLUMNS = [
    ('uniprot_accession', 'string', '', 'reviewed UniProtKB accession'),
    ('uniprot_entry_name', 'string', '', 'UniProtKB entry name'),
    ('gene', 'list<string>', '', 'primary gene name of the UniProt entry (several for a few entries)'),
    ('gene_synonyms', 'list<string>', '', 'other gene names UniProt gives the entry (release 2026_03)'),
    ('reported_genes', 'list<string>', '', 'other gene symbols the source records used, most frequent first'),
    ('protein_name', 'string', '', 'recommended protein name'),
    ('protein_length', 'integer', '', 'length of the reference sequence'),
    ('reference_sequence_md5', 'string', '', 'MD5 of the reference sequence, uppercase hex as UniProt gives it'),
    ('n_sites', 'integer', '', 'methylarginine sites, not counting the sites kept pending'),
    ('n_pending', 'integer', '', 'sites kept pending (flag one_sided_window)'),
    ('studies', 'list<code>', 'study', 'studies or databases reporting its sites, in the order of the studies file'),
    ('panther_group', 'string', '', 'group of its PANTHER 19.0 top-level protein class (config/panther_class_groups.json, scripts/03b), or No PANTHER protein class'),
    ('panther_class', 'string', '', 'its PANTHER top-level protein class, empty when PANTHER gives none'),
    ('panther_classes', 'list<string>', '', 'its most specific PANTHER protein classes'),
    ('panther_family', 'string', '', 'its PANTHER family (PTHR id)'),
    ('functional_classes', 'list<string>', '', 'broad functional classes (config/functional_ontology.json, scripts/03)'),
    ('functional_subclasses', 'list<string>', '', 'functional subclasses'),
    ('arginine_count', 'integer', '', 'arginines in the reference sequence'),
    ('arg_enrichment_log2', 'number', '', 'shrinkage estimate of the log2 enrichment of methylarginines over the proteome-wide rate (scripts/08)'),
    ('arg_enrichment_q', 'number', '', 'q value of the per-protein test (scripts/08), empty when not tested'),
    ('nearest_spacing_min', 'number', '', 'smallest distance in residues from a site to the nearest other site (scripts/10)'),
    ('nearest_spacing_median', 'number', '', 'median of these distances'),
    ('max_sites_3aa', 'integer', '', 'most sites in any 3-residue window (scripts/10)'),
    ('max_sites_5aa', 'integer', '', 'most sites in any 5-residue window'),
    ('max_sites_10aa', 'integer', '', 'most sites in any 10-residue window'),
    ('max_sites_20aa', 'integer', '', 'most sites in any 20-residue window'),
    ('max_sites_50aa', 'integer', '', 'most sites in any 50-residue window'),
    ('best_window_aa', 'integer', '', 'width of the protein\'s most significant density window: scripts/28 compares the most sites in any window of 5 to 100 residues with a null that places the sites on the protein\'s own arginines at random (1,000 permutations, proteins with at least 3 sites), and this is the window with the lowest q value (then the largest excess)'),
    ('best_window_start', 'integer', '', 'its first residue'),
    ('best_window_end', 'integer', '', 'its last residue'),
    ('best_window_sites', 'integer', '', 'sites in it'),
    ('best_window_log2_obs_null', 'number', '', 'log2 of observed over null mean sites'),
    ('best_window_q', 'number', '', 'its q value'),
    ('n_dense_windows', 'integer', '', 'windows with q <= 0.05 (supporting/density_windows)'),
    ('architecture_class', 'string', '', 'domain-cluster architecture of the clustered disordered sites (scripts/25: at least 3 sites within 20 residues), empty when none'),
    ('architecture_clusters', 'list<string>', '', 'the site positions of each cluster, comma-separated'),
    ('condensates', 'list<string>', '', 'CD-CODE condensates the protein belongs to (CD-CODE ids, export of 2026-03-18, scripts/12)'),
    ('hpa_locations', 'list<string>', '', 'Human Protein Atlas 25.1 subcellular locations, main locations first (supporting/hpa_locations)'),
    ('disease_families', 'list<string>', '', 'disease families whose genes include the protein\'s gene (scripts/14)'),
    ('disease_family_terms', 'list<integer>', '', 'for each of these families, how many of its DisGeNET terms name the gene'),
    ('als_ftd_family_gene', 'boolean', '', 'the gene is in the ALS / motor neuron / FTD family'),
    ('n_phospho', 'integer', '', 'phosphorylation sites in the PTM master table (supporting/other_ptm_sites)'),
    ('n_lys_acyl', 'integer', '', 'lysine acylation sites'),
    ('n_ub_sumo', 'integer', '', 'ubiquitination or SUMOylation sites'),
    ('n_lys_methyl', 'integer', '', 'lysine methylation sites'),
    ('n_sites_experimental', 'integer', '', 'sites with interaction_status experimental'),
    ('n_sites_predicted_only', 'integer', '', 'sites with interaction_status predicted_only'),
    ('go_rna_binding', 'code', 'go_binding', 'Gene Ontology RNA binding (GO:0003723 or a child term)'),
    ('go_dna_binding', 'code', 'go_binding', 'Gene Ontology DNA binding (GO:0003677 or a child term)'),
    ('n_string_partners', 'integer', '', 'STRING 12.5 physical partners at score >= 0.7 (supporting/protein_partners)'),
    ('n_partners_methylated', 'integer', '', 'of these, proteins of the site table'),
    ('n_complexes', 'integer', '', 'Complex Portal and hu.MAP complexes the protein is in (supporting/complexes)'),
]

STATISTICS_COLUMNS = [
    ('analysis', 'code', 'analysis', 'what was tested'),
    ('category', 'string', '', 'the class, motif, family, condensate (CD-CODE id) or location tested'),
    ('stratum', 'string', '', 'all, or the subset the test is restricted to'),
    ('a', 'integer', '', 'targets in the category'),
    ('b', 'integer', '', 'targets not in the category'),
    ('c', 'integer', '', 'background in the category'),
    ('d', 'integer', '', 'background not in the category'),
    ('odds_ratio', 'number', '', 'odds ratio, for the 2x2 tests (methods fisher, mantel_haenszel); empty for method enricher'),
    ('log2_odds_ratio', 'number', '', 'its log2'),
    ('fold_enrichment', 'number', '', 'fold enrichment (k/n over M/N) as clusterProfiler reports it, for method enricher; empty otherwise'),
    ('log2_fold_enrichment', 'number', '', 'its log2'),
    ('p', 'number', '', 'p value'),
    ('q', 'number', '', 'Benjamini-Hochberg adjusted p value within the analysis and method, over all its strata and categories'),
    ('method', 'code', 'method', 'the test'),
    ('background', 'string', '', 'what the targets are compared with'),
    ('n_strata', 'integer', '', 'strata of a stratified test, empty otherwise'),
]

DISEASE_TERM_COLUMNS = [
    ('term_id', 'string', '', 'DisGeNET term (UMLS CUI)'),
    ('term_name', 'string', '', 'its name'),
    ('family', 'string', '', 'its disease family (scripts/disease_families.py): audited, or else additional; empty when none'),
    ('n_genes_term', 'integer', '', 'genes of the term in the universe'),
    ('n_methylated', 'integer', '', 'of these, the genes of methylarginine-containing proteins'),
    ('fold_enrichment', 'number', '', 'observed over expected'),
    ('p', 'number', '', 'hypergeometric p value (DOSE enrichDGN)'),
    ('p_adjust', 'number', '', 'Benjamini-Hochberg adjusted p value over the tested terms'),
    ('q', 'number', '', 'q value (qvalue, lambda 0.05, as DOSE)'),
    ('genes', 'list<string>', '', 'the genes of the methylarginine-containing proteins in the term (symbols, org.Hs.eg.db)'),
    ('uniprot_accessions', 'list<string>', '', 'the proteins of the site table those genes stand for'),
]

DENSITY_COLUMNS = [
    ('uniprot_accession', 'string', '', 'protein'),
    ('window_aa', 'integer', '', 'window width in residues'),
    ('start', 'integer', '', 'first residue of the densest window of that width'),
    ('end', 'integer', '', 'last residue'),
    ('sites', 'integer', '', 'methylarginine sites in it'),
    ('arginines', 'integer', '', 'arginines in it'),
    ('null_mean_sites', 'number', '', 'mean of the most sites in a window under the null (sites placed on the protein\'s arginines at random)'),
    ('log2_obs_over_null', 'number', '', 'log2 of observed over null mean'),
    ('p', 'number', '', 'empirical p value (1,000 permutations)'),
    ('q', 'number', '', 'Benjamini-Hochberg q value over all protein-window tests. Only windows with q <= 0.05 are listed'),
    ('site_positions', 'list<integer>', '', 'the sites in the window'),
]

HPA_COLUMNS = [
    ('uniprot_accession', 'string', '', 'protein'),
    ('gene', 'string', '', 'gene as HPA gives it'),
    ('ensembl_gene', 'string', '', 'Ensembl gene HPA maps to the protein'),
    ('location', 'string', '', 'subcellular location'),
    ('main', 'boolean', '', 'a main location (otherwise additional)'),
    ('reliability', 'string', '', 'HPA reliability: Enhanced, Supported, Approved or Uncertain'),
]

OTHER_PTM_COLUMNS = [
    ('uniprot_accession', 'string', '', 'protein'),
    ('position', 'integer', '', 'position on the reference sequence'),
    ('residue', 'string', '', 'the modified residue'),
    ('ptm_type', 'code', 'ptm_type', 'the modification'),
    ('source', 'string', '', 'ptm_master: the iPTMnet-format PTM master table (rows on the canonical accession, residue checked against the reference sequence)'),
    ('evidence', 'list<string>', '', 'the databases the rows come from, as the master table names them'),
    ('n_records', 'integer', '', 'rows of the master table'),
]

CONDENSATE_COLUMNS = [
    ('condensate_id', 'string', '', 'CD-CODE condensate id'),
    ('condensate_name', 'string', '', 'its name'),
    ('condensate_type', 'string', '', 'its type in CD-CODE'),
    ('n_members', 'integer', '', 'human proteins CD-CODE lists as members'),
    ('n_members_methylated', 'integer', '', 'of these, proteins of the site table'),
]

DISEASE_TERM_TESTED_COLUMNS = [
    ('term_id', 'string', '', 'DisGeNET term (UMLS CUI)'),
    ('term_name', 'string', '', 'its name (a ";" inside a name, the list separator, is written ",")'),
    ('family', 'string', '', 'its disease family (scripts/disease_families.py); empty when none'),
    ('n_genes_term', 'integer', '', "genes of the term in DOSE's background (reviewed human genes with a DisGeNET term)"),
    ('n_methylated', 'integer', '', 'of these, the genes of methylarginine-containing proteins'),
    ('fold_enrichment', 'number', '', 'fold enrichment as DOSE reports it'),
    ('p', 'number', '', 'one-sided hypergeometric p value (DOSE enrichDGN)'),
    ('p_adjust', 'number', '', 'Benjamini-Hochberg adjusted p value over the tested terms'),
    ('q', 'number', '', 'q value (DOSE)'),
    ('enriched', 'boolean', '', 'p, p_adjust and q at most 0.05, the cut of the published analysis (supporting/disease_terms)'),
    ('genes', 'list<string>', '', 'the genes of the methylarginine-containing proteins in the term (symbols)'),
    ('uniprot_accessions', 'list<string>', '', 'the proteins of the site table those genes stand for'),
]

DISEASE_ONTOLOGY_COLUMNS = [
    ('term_id', 'string', '', 'Human Disease Ontology term (DOID)'),
    ('term_name', 'string', '', 'its name (a ";" inside a name is written ",")'),
    ('parents', 'list<string>', '', 'its parent terms (is_a)'),
    ('n_genes_term', 'integer', '', 'genes of the term and of the terms under it, among the reviewed human genes with a Disease Ontology term (the gene sets of DOSE enrichDO)'),
    ('n_methylated', 'integer', '', 'of these, the genes of methylarginine-containing proteins'),
    ('fold_enrichment', 'number', '', 'fold enrichment (clusterProfiler)'),
    ('p', 'number', '', 'one-sided hypergeometric p value (clusterProfiler::enricher on the Disease Ontology gene sets; the same as DOSE enrichDO)'),
    ('q', 'number', '', 'Benjamini-Hochberg adjusted p value over the tested terms (at least 10 genes)'),
    ('genes', 'list<string>', '', 'the genes of the methylarginine-containing proteins in the term (symbols)'),
    ('uniprot_accessions', 'list<string>', '', 'the proteins of the site table those genes stand for'),
]

DISEASE_FAMILY_COLUMNS = [
    ('family', 'string', '', 'disease family (scripts/disease_families.py): the eight audited families, then the additional ones, each in order of significance'),
    ('audited', 'boolean', '', 'one of the eight audited families (scripts/disease_families.py FAMILY_RULES, the original pipeline\'s hand-checked rules); otherwise an additional, rule-based family'),
    ('n_terms', 'integer', '', 'DisGeNET terms assigned to it'),
    ('n_genes', 'integer', '', 'genes of the universe in at least one of its terms'),
    ('n_methylated_genes', 'integer', '', 'of these, the genes of methylarginine-containing proteins'),
    ('fold_enrichment', 'number', '', 'fold enrichment of the family (clusterProfiler::enricher against the reviewed human genes with a DisGeNET term; statistics: disease_family)'),
    ('q', 'number', '', 'its Benjamini-Hochberg adjusted p value (over the audited families, or over the additional ones)'),
]

SITE_INTERACTION_COLUMNS = [
    ('site_id', 'string', '', 'the site (site table)'),
    ('evidence_type', 'code', 'evidence_type', 'the kind of evidence'),
    ('localisation', 'code', 'localisation', 'where the evidence falls relative to the site'),
    ('partner_kind', 'code', 'partner_kind', 'what the arginine contacts, empty when the record names no partner'),
    ('partner_accession', 'string', '', 'the partner protein (reviewed UniProtKB accession)'),
    ('partner_gene', 'string', '', 'its primary gene name'),
    ('source', 'string', '', 'PDBe-KB, ProtVar, AFDB (AlphaFold DB), IntAct, RBDmap, Trendel2025, UniProt, or the RNA crosslink dataset (Kramer2014, XRNAX, RBS-ID, pRBS-ID, pepR-MS, pepR-MS_focused)'),
    ('structures', 'list<string>', '', 'PDB entries (PDB rows), or the predicted model (ProtVar pair or AlphaFold DB model id)'),
    ('score', 'number', '', 'pDockQ (ProtVar) or ipTM (AlphaFold DB) of the predicted complex'),
    ('distance_a', 'number', '', 'closest distance in the predicted model between the arginine and the partner, in angstroms'),
    ('detail', 'string', '', 'for an RNA crosslink, the crosslinked peptides of that dataset with the crosslinked residue in brackets, its position and the RNA adduct masses (pepR-MS compendium Table S2); for an RBDmap fragment, its range and sequence; for a DNA crosslink, the crosslinked peptide with the residues the crosslink may sit on in lower case (Trendel et al. 2025); bonds, pLDDT and partner residues of a predicted contact; effect, residues, PubMed and IntAct ids of an IntAct record'),
    ('counts_as_evidence', 'boolean', '', 'the record is site-level evidence: crosslinks at or as a candidate residue, RBDmap fragments, PDB interfaces, predicted interfaces, confident predicted contacts, and IntAct mutations that change an interaction, short binding regions and methylated residues. interaction_status counts all of these except IntAct'),
]

RESIDUE_EVIDENCE_COLUMNS = [
    ('uniprot_accession', 'string', '', 'protein of the site table'),
    ('start', 'integer', '', 'first residue of a run of consecutive residues with the same evidence'),
    ('end', 'integer', '', 'last residue'),
    ('evidence_type', 'code', 'evidence_type', 'the kind of evidence'),
    ('localisation', 'code', 'localisation', 'at the residues, or a region or peptide that holds them'),
    ('n_partners', 'integer', '', 'partners: most along the run (PDB protein, predicted interface) or records (IntAct)'),
    ('best_score', 'number', '', 'best pDockQ along the run (predicted interfaces)'),
    ('sources', 'list<string>', '', 'sources of the evidence'),
]

PARTNER_COLUMNS = [
    ('uniprot_accession', 'string', '', 'protein of the site table'),
    ('partner_accession', 'string', '', 'STRING 12.5 physical partner at score >= 0.7 (reviewed UniProtKB accession)'),
    ('partner_gene', 'string', '', 'its primary gene name'),
    ('string_score', 'number', '', 'STRING physical interaction score (0-1)'),
    ('partner_methylated', 'boolean', '', 'the partner is a protein of the site table'),
    ('partner_classes', 'list<code>', 'partner_class', 'reader or writer classes of the partner'),
]

COMPLEX_COLUMNS = [
    ('complex_id', 'string', '', 'Complex Portal id'),
    ('complex_name', 'string', '', 'its name (empty for most predicted complexes)'),
    ('curation', 'code', 'curation', 'curated or predicted'),
    ('n_members', 'integer', '', 'reviewed member proteins (each protein that can fill a member slot counts)'),
    ('n_members_methylated', 'integer', '', 'members that are proteins of the site table'),
    ('methylated_members', 'list<string>', '', 'those proteins'),
    ('fold_enrichment', 'number', '', 'clusterProfiler FoldEnrichment: methylated share of the members over the methylated share of reviewed human proteins'),
    ('p', 'number', '', 'clusterProfiler enricher p value (one-sided hypergeometric, P(X >= n_members_methylated)), against all reviewed human proteins'),
    ('q', 'number', '', 'Benjamini-Hochberg over the tested complexes of the same curation (complexes of at least two reviewed proteins with a methylated member). Only complexes with a methylated member are listed; those of a single reviewed protein are listed without a test'),
]

"""Codes and column definitions of a release (used by scripts/build_release.py and scripts/release_docs.py).

Every coded column of a release table holds codes only, what each code means is here, and the release's
data_dictionary.tsv is written from this file. A column is (name, type, codes, description): type is one of integer,
number, string, boolean, code, list<code>, list<integer>, list<string>, codes names an entry of CODES (or is empty).
"""
from __future__ import annotations

CONTEXT_WINDOW = 20  # residues: the boundary window of scripts/05 and scripts/06
PADJ_CUTOFF = 0.05   # Maron 2021 Table S3: a peptide changes when its adjusted p value is at most this (scripts/16)

# Studies and databases behind the union, in the order of the studies file (most sites first, as in v1.1.0).
# Citations and DOIs checked against Crossref.
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
]
STUDY_ORDER = [s[0] for s in STUDIES]
MARON_S5 = 'Maron 2021 Table S5 (CC BY-NC-ND 4.0, open-access supplement)'
PTM_MASTER = 'the iPTMnet-format PTM master table (iPTMnet: CC BY-NC-SA 4.0)'
STUDY_LICENCE = {
    'PSP': 'CC BY-NC-SA 3.0 Unported (PhosphoSitePlus terms of 2021, inputs/licences/), compiled in ' + MARON_S5,
    'ProMetheusDB': 'CC BY-NC-ND 4.0 (open-access supplement mmc2)',
    'UniProt': 'CC BY 4.0, via ' + PTM_MASTER,
    'IEDB': 'via ' + PTM_MASTER,
}
STUDY_NOTE = {
    'PSP': 'downloaded 22 May 2021 (date from the Shechter lab), as compiled in Maron et al. 2021 Table S5',
    'Maron2021': "the study's own sites in Table S5 (reference \"This study\")",
}

SOURCE_TABLES = {
    'maron2021_s5': ('maron_s5', 'Maron et al. 2021 Table S5: methylarginine sites compiled from ten studies'),
    'prometheusdb_mmc2': ('prometheus_mmc2', 'ProMetheusDB (Massignani et al. 2022) supplementary table mmc2'),
    'ptm_master': ('human_ptm_master', 'the arginine methylation rows of the iPTMnet-format PTM master table (UniProt and IEDB records)'),
}
SOURCE_TABLE_CODE = {raw: code for code, (raw, _) in SOURCE_TABLES.items()}

CODES: dict[str, dict[str, str]] = {
    'status': {
        'site': 'a site of the union: at least one source record is placed on it',
        'pending': 'a site whose only records are windows crossing a splice junction (flag pending_splice_junction). Kept until the policy is decided, and left out of every analysis',
    },
    'methyl_state': {
        'Rme1': 'monomethylarginine',
        'Rme2': 'dimethylarginine, asymmetric or symmetric not stated',
        'Rme2a': 'asymmetric dimethylarginine',
        'Rme2s': 'symmetric dimethylarginine',
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
        'predicted_only': 'only predicted interface or contact evidence',
        'none': 'no interaction evidence',
    },
    'evidence_type': {
        'rna_crosslink': 'RNA-protein crosslink (pepR-MS re-analysis of Kramer 2014, XRNAX, RBS-ID, pRBS-ID, pepR-MS focused sites)',
        'rbdmap_fragment': 'inside an RNA-bound proteolytic fragment (RBDmap, Castello et al. 2016)',
        'dna_crosslink': 'DNA-protein crosslink (Trendel et al. 2025)',
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
    'flag': {
        'pending_splice_junction': 'a window crossing a splice junction: one full flank (at least 6 residues) matches the entry at the reported arginine, the other does not, and no other entry holds the window. The pipeline leaves it unplaced and the release keeps it as a pending site',
        'paralog_window': 'a supporting window matches a paralog exactly and this entry within 2 residues (open, accepted in config/remap_evidence_allowlist.tsv)',
        'near_match_only': 'every supporting record matches the entry within 2 residues, not exactly (a variant or sequence conflict)',
        'shared_sequence': 'residues -15 to +15 occur identically in another reviewed entry (same_window_accessions): histone variants, identical paralogs. The evidence cannot tell them apart: the protein is the one the source named',
        'single_study': 'one study or database reports the site',
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
    },
    'reason': {
        'peptide_recentred': 'the previous version moved the peptide to the arginine nearest its middle. This version keeps the reported arginine inside the peptide',
        'unreviewed_accession': "the previous version filed reviewed-entry numbering under the source's unreviewed or obsolete accession. This version files it under the reviewed entry",
        'lysine_window': 'the source window is centred on a lysine. The previous version turned it into a nearby arginine site',
        'truncated_window': 'the window is cut short at a protein end. The previous version read it from the wrong end',
        'multi_arginine_peptide': 'a peptide with several arginines that does not cover the reported position. The previous version guessed one',
        'peptide_start': "a peptide with several arginines, numbered differently from the entry. The previous version took the arginine nearest its middle. This version counts the reported position from the peptide's start in Maron 2021 Tables S2/S3",
        'wrong_protein': 'the previous version filed the position under another reviewed protein than the one whose sequence matched',
        'outdated_accession_or_symbol': 'the source accession is merged or obsolete, or its gene symbol outdated. This version reaches the current entry',
        'evidence_mismatch': 'the previous version kept the reported arginine although the source peptide or window does not match the entry there',
        'window_padding': 'the window is padded where it runs off a protein end. The padding now gives the cut side',
        'ambiguous_evidence': 'the evidence supports several sites equally. This version leaves the record unplaced',
        'pending_splice_junction': 'the site is kept as pending (flag pending_splice_junction)',
        'reference_renumbered': 'the same residue under a new number: UniProt 2026_03 extends the sequence at the N-terminus (RAB12, +96)',
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
        'unresolved_no_anchor': 'not placed: no evidence names an arginine and the reported residue is not one',
        'unresolved_no_hit': 'not placed: evidence not found in any candidate entry',
        'unresolved_ambiguous_hit': 'not placed: evidence supports several sites equally',
        'unresolved_evidence_elsewhere': 'not placed: the evidence lies elsewhere in the entry, not at the reported position',
        'unresolved_evidence_mismatch': 'not placed: the evidence does not match the entry at the reported position',
    },
    'evidence_kind': {
        'window_13': 'a 13-residue window centred on the site',
        'window_15': 'a 15-residue window centred on the site',
        'window_31': 'a 31-residue window centred on the site',
        'peptide': 'the identified peptide',
        'marked_peptide': 'the peptide with the methylated residue marked',
        'position_only': 'no sequence: the reported position only',
    },
    'evidence_check': {
        'agrees': 'placed on the one site its evidence names',
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
    'agrees': 'agrees', 'agrees: no site': 'agrees_no_site', 'one of several evidence sites': 'one_of_several',
    'evidence ambiguous': 'ambiguous', 'dropped, evidence has a site': 'dropped_evidence_has_site',
    'contradicts evidence': 'contradicts',
}
FLAG_ORDER = list(CODES['flag'])
MOTIF_ORDER = list(CODES['motif'])
ANTIBODY_ORDER = list(CODES['antibody'])
EVIDENCE_ORDER = list(CODES['evidence_type'])
CALL_ORDER = list(CODES['call'])

# ---- Columns of each file -------------------------------------------------------------------------------------------

SITE_COLUMNS = [
    ('site_id', 'string', '', 'UniProt accession and arginine position, e.g. P35637:R216'),
    ('status', 'code', 'status', 'site, or pending'),
    ('uniprot_accession', 'string', '', 'reviewed UniProtKB accession of the canonical sequence'),
    ('uniprot_entry_name', 'string', '', 'UniProtKB entry name'),
    ('gene', 'list<string>', '', 'primary gene name of the UniProt entry (several for a few entries, e.g. histone genes), empty when UniProt gives none'),
    ('protein_name', 'string', '', 'recommended protein name of the UniProt entry'),
    ('position', 'integer', '', 'position of the arginine in the reference sequence (1-based, initiator methionine = 1)'),
    ('residue', 'code', '', 'residue at that position in the reference sequence: always R'),
    ('window_15', 'string', '', 'residues -7 to +7 around the site, _ beyond a protein end'),
    ('window_31', 'string', '', 'residues -15 to +15 around the site, _ beyond a protein end'),
    ('protein_length', 'integer', '', 'length of the reference sequence'),
    ('reference_sequence_md5', 'string', '', 'MD5 of the reference sequence, uppercase hex as UniProt gives it (a later UniProt entry with the same MD5 has the same numbering)'),
    ('methyl_states_reported', 'list<code>', 'methyl_state', 'methylation degrees given by the source records, in code order, empty where they give none'),
    ('maron2021_s2_states', 'list<code>', 'methyl_state', 'degrees called at the site in Maron et al. 2021 Table S2'),
    ('maron2021_s3_antibodies', 'list<code>', 'antibody', 'antibodies whose immunoprecipitation recovered the site in Maron et al. 2021 Table S3 (supporting/maron2021_detection)'),
    ('maron2021_GSK591_call', 'code', 'call', f'response to GSK591 in Maron et al. 2021 Table S3, combined over antibodies (supporting/maron2021_inhibitors: a peptide changes at adjusted p <= {PADJ_CUTOFF}), empty when not measured'),
    ('maron2021_MS023_call', 'code', 'call', 'response to MS023, as for GSK591'),
    ('carm1_ko_call', 'code', 'call', 'response to CARM1 knockout (Shishkova et al. 2017), combined over cell lines (supporting/carm1_ko), empty when not measured'),
    ('motifs', 'list<code>', 'motif', 'sequence motifs around the site (scripts/08), in code order'),
    ('disorder', 'code', 'disorder', 'MobiDB predicted disorder at the site (scripts/06)'),
    ('disorder_edge_distance', 'integer', '', 'residues to the nearest edge of a disordered region, empty with no annotation'),
    ('domain', 'code', 'domain', 'InterPro 108.0 domain or repeat context (scripts/05)'),
    ('nearest_domain_interpro', 'string', '', 'InterPro accession of the nearest domain or repeat (the innermost when inside several)'),
    ('nearest_domain_name', 'string', '', 'its InterPro name'),
    ('nearest_domain_start', 'integer', '', 'its first residue'),
    ('nearest_domain_end', 'integer', '', 'its last residue'),
    ('nearest_domain_class', 'string', '', 'its class in config/domain_class_ontology.json, as scripts/07 assigns it, No domain annotation when none'),
    ('domain_side', 'code', 'domain_side', 'where the site lies relative to the nearest domain, empty when none'),
    ('domain_edge_distance', 'integer', '', 'residues to the nearest domain edge (0 inside a domain), empty when none'),
    ('interaction_status', 'code', 'interaction_status', 'interaction evidence at the site (site_interactions)'),
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
    ('note', 'string', '', ''),
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
    ('reported_methyl_state', 'string', '', 'methylation state as reported, verbatim'),
    ('evidence_kind', 'code', 'evidence_kind', 'kind of sequence evidence'),
    ('evidence_sequence', 'list<string>', '', 'the window, peptide or marked peptide, verbatim (a few records give several)'),
    ('placed_site', 'string', '', 'the site this version places the record on, empty when not placed'),
    ('placement_code', 'code', 'placement', 'how the record was placed (scripts/02 remap_status)'),
    ('evidence_check', 'code', 'evidence_check', "independent check of the placement against the record's own evidence (scripts/qc_remap_evidence_check.py)"),
    ('accepted_exception', 'string', '', 'why a contradiction is accepted for now (config/remap_evidence_allowlist.tsv)'),
    ('pending_site', 'string', '', 'the pending site the record supports (flag pending_splice_junction), empty otherwise'),
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
    ('in_this_version', 'boolean', '', 'a site of this version (status site or pending)'),
    ('change', 'code', 'site_change', 'kept, removed or added'),
    ('reason_codes', 'list<code>', 'reason', 'why it was removed, added or renumbered'),
    ('records_now_at_or_previously_at', 'list<string>', '', 'for a removed site, where its records are placed now, for an added site, where they were placed before'),
    ('record_ids', 'list<integer>', '', 'the records involved (supporting/source_records)'),
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
    ('call', 'code', 'call', 'a peptide decreases at log2 <= -1 and p < 0.05, increases at log2 >= 1 and p < 0.05 (scripts/41), combined over peptides'),
]

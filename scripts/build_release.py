#!/usr/bin/env python3
"""Build a release of the corrected human methylarginine sites, with the provenance of every site and record.

Reads a corrected run made from the frozen inputs (remap rows, union, evidence check; see hpc/build_release.sh), the
previous version (the published union, PTM_results/) and the frozen inputs (inputs/), and writes to --outdir, for
version V and previous version P:

  human_methylarginine_sites_V.tsv            one row per site: identity, sequence, methylation evidence, structural
                                              context, provenance, flags and change since P
  human_methylarginine_sites_V.xlsx           the same and the tables below, as one annotated workbook
  human_methylarginine_source_records_V.tsv   one row per source record: what the study reported, where it was placed
  human_methylarginine_changes_P_to_V.tsv     every site of P or V: kept, removed or added, and why
  human_methylarginine_studies_V.tsv          one row per study or database: citation, DOI, records and sites
  README.md                                   what the files hold and how they were built
  MD5SUMS                                     checksums of the files above (check with: md5sum -c MD5SUMS)

Structural context and motifs reuse the pipeline's own code (scripts/05, 06, 08), so they match the published analyses.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))


def load_script(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


remap = load_script('remap', '02_remap_arg_methyl_to_canonical.py')
domain_context = load_script('domain_context', '05_annotate_domain_context.py')
disorder_context = load_script('disorder_context', '06_annotate_disorder_context.py')
motif_odds = load_script('motif_odds', '08_motif_and_arg_odds.py')

# Frozen inputs read here, relative to --inputs (checksums in inputs/MANIFEST.tsv).
INPUT = {
    'fasta': 'reference/uniprot_human_reviewed_2026-04-07.fasta',
    'entries': 'reference/uniprot_human_reviewed_2026-04-07_entries.tsv',
    'disorder': 'annotation/mobidb_disorder_intervals_2026-03-30.tsv',
    'domains': 'annotation/interpro_domain_intervals_2026-04-07.tsv',
    'maron_s2': 'sources/maron2021_tableS2_methylarginine_ip_sites.xlsx',
    'maron_s3': 'sources/maron2021_tableS3_methylpeptides_prmt_inhibitors.xlsx',
}
CONTEXT_WINDOW = 20  # residues: the boundary window of scripts/05 and scripts/06
SIGNIFICANT = 0.05   # adjusted p value for an inhibitor response call

# Studies and databases behind the union: DOI (as on the website, docs/assets/explorer.js STUDY_DOI) and citation,
# both checked against Crossref.
STUDIES = {
    'Maron 2021': ('10.1016/j.isci.2021.102971', 'Maron et al., iScience 24:102971 (2021)'),
    'ProMetheusDB': ('10.1016/j.mcpro.2022.100243', 'Massignani et al., Mol Cell Proteomics 21:100243 (2022)'),
    'PhosphoSitePlus': ('10.1093/nar/gku1267', 'Hornbeck et al., Nucleic Acids Res 43:D512-D520 (2015)'),
    'Larsen 2016': ('10.1126/scisignal.aaf7329', 'Larsen et al., Sci Signal 9 (2016)'),
    'Musiani 2019': ('10.1126/scisignal.aat8388', 'Musiani et al., Sci Signal 12:eaat8388 (2019)'),
    'Fong 2019': ('10.1016/j.ccell.2019.07.003', 'Fong et al., Cancer Cell 36:194-209 (2019)'),
    'Guo 2014': ('10.1074/mcp.O113.027870', 'Guo et al., Mol Cell Proteomics 13:372-387 (2014)'),
    'Fedoriw 2019': ('10.1016/j.ccell.2019.05.014', 'Fedoriw et al., Cancer Cell 36:100-114 (2019)'),
    'Wei 2020': ('10.3389/fcell.2020.545934', 'Wei et al., Front Cell Dev Biol 8:545934 (2020)'),
    'Lim 2020': ('10.1186/s12953-020-00162-8', 'Lim et al., Proteome Sci 18:6 (2020)'),
    'Li 2021': ('10.1038/s41467-021-21963-1', 'Li et al., Nat Commun 12:1946 (2021)'),
    'IEDB': ('10.1093/nar/gky1006', 'Vita et al., Nucleic Acids Res 47:D339-D343 (2019)'),
    'UniProt': ('10.1093/nar/gkae1010', 'UniProt Consortium, Nucleic Acids Res 53:D609-D617 (2025)'),
}
SOURCE_TABLE = {'maron_s5': 'Maron 2021 Table S5', 'prometheus_mmc2': 'ProMetheusDB mmc2', 'human_ptm_master': 'PTM master table'}
SOURCE_TABLE_ABOUT = {
    'Maron 2021 Table S5': 'compilation of ten studies, with Maron 2021\'s own sites (inputs/sources/maron2021_tableS5_*)',
    'ProMetheusDB mmc2': 'ProMetheusDB methyl-site table (inputs/sources/massignani2022_prometheusdb_mmc2_*)',
    'PTM master table': 'UniProt and IEDB arginine methylation rows of the iPTMnet-format master table (inputs/sources/ptm_master_iptmnet_human.tsv)',
}
PSP_TERMS = 'PhosphoSitePlus data: CC BY-NC-SA 3.0 (non-commercial, share alike)'

# How each record was placed, in words (scripts/02 remap_status).
PLACEMENT = {
    'sequence_unique_canonical': 'evidence sequence found once in the reported entry',
    'sequence_unique_candidate': 'evidence sequence found once, in another entry of the same gene or the current entry of a merged accession',
    'sequence_closest_supported': 'evidence sequence found more than once; the occurrence nearest the reported position',
    'sequence_position_verified': 'reported arginine lies inside the peptide in the reference sequence',
    'sequence_unique_proteome': 'evidence sequence found once in the whole reference proteome (no candidate entry held it)',
    'same_accession_position_verified': 'position-only record; the reported residue is an arginine of the entry',
    'canonical_position_fallback': 'position-only record; the reported residue is an arginine of the canonical entry',
    'same_accession_evidence_near_match': 'evidence not in the entry, but matches it at the reported arginine within 2 residues (variant or conflict)',
    'canonical_evidence_near_match': 'evidence not in the canonical entry, but matches it at the reported arginine within 2 residues',
    'unresolved_no_anchor': 'not placed: no evidence names an arginine and the reported residue is not one',
    'unresolved_no_hit': 'not placed: evidence not found in any candidate entry',
    'unresolved_ambiguous_hit': 'not placed: evidence supports several sites equally',
    'unresolved_evidence_elsewhere': 'not placed: the evidence lies elsewhere in the entry, not at the reported position',
    'unresolved_evidence_mismatch': 'not placed: the evidence does not match the entry at the reported position',
}

# Why a record's placement differs from the previous version; numbers are the issues of DATA_PROVENANCE_REVIEW.md.
REASONS = {
    'peptide_recentred': ('1', 'the published run moved the peptide to the arginine nearest its middle; the release keeps the reported arginine inside the peptide'),
    'unreviewed_accession': ('2', 'the published run filed reviewed-entry numbering under the source\'s unreviewed or obsolete accession; the release files it under the reviewed entry'),
    'lysine_window': ('3', 'the source window is centred on a lysine; the published run turned it into a nearby arginine site'),
    'truncated_window': ('4', 'the window is cut short at a protein end; the published run read it from the wrong end'),
    'multi_arginine_peptide': ('5', 'a peptide with several arginines that does not cover the reported position; the published run guessed one'),
    'wrong_protein': ('6', 'the published run filed the position under another reviewed protein than the one whose sequence matched'),
    'outdated_accession_or_symbol': ('7', 'the source accession is merged or obsolete, or its gene symbol outdated; the release reaches the current entry'),
    'evidence_mismatch': ('8', 'the published run kept the reported arginine although the source peptide or window does not match the entry there'),
    'window_padding': ('13', 'the window is padded where it runs off a protein end; the padding now gives the cut side'),
    'ambiguous_evidence': ('14', 'the evidence supports several sites equally; the release leaves the record unplaced'),
    'other': ('', 'placement changed for a reason not covered by the codes above'),
}
PENDING_JUNCTION = 'pending_splice_junction'
FLAGS = {
    PENDING_JUNCTION: 'a window crossing a splice junction: one full flank (at least 6 residues) matches the entry exactly at the reported arginine, the other does not, and no other entry holds the window. The pipeline leaves it unplaced; the release includes it until the policy is decided',
    'paralog_window': 'a supporting window matches a paralog exactly and this entry within 2 residues (issue 9, open; accepted in config/remap_evidence_allowlist.tsv)',
    'near_match_only': 'every supporting record matches the entry within 2 residues, not exactly (a variant or sequence conflict)',
    'shared_sequence': 'residues -15 to +15 occur identically in another reviewed entry (listed in same_window_in): histone variants, identical paralogs. The evidence cannot tell the proteins apart; the gene is the one the source named',
    'single_study': 'one study or database reports the site',
}

# Column dictionary: (column, group, description). The site table's groups head the workbook's Sites sheet.
SITE_COLUMNS = [
    ('site_id', 'Site', 'UniProt accession and arginine position, e.g. P35637:R216'),
    ('gene', 'Site', 'primary gene name of the UniProt entry'),
    ('protein_name', 'Site', 'recommended protein name of the UniProt entry'),
    ('uniprot_accession', 'Site', 'reviewed canonical UniProtKB accession (linked to its UniProt page in the workbook)'),
    ('uniprot_entry_name', 'Site', 'UniProtKB entry name'),
    ('position', 'Site', 'position of the arginine in the reference sequence (1-based)'),
    ('residue', 'Site', 'residue at that position in the reference sequence (always R)'),
    ('window_15', 'Sequence', 'residues -7 to +7 around the site; _ beyond a protein end'),
    ('window_31', 'Sequence', 'residues -15 to +15 around the site; _ beyond a protein end'),
    ('protein_length', 'Sequence', 'length of the reference sequence'),
    ('reference_sequence_md5', 'Sequence', 'MD5 of the reference sequence (UniProt sequence MD5 format: uppercase hex); compare with a later UniProt entry to see whether the numbering still holds'),
    ('methyl_states_reported', 'Methylation evidence', 'methylation degree given by the source records: me1 (monomethyl), me2 (dimethyl, asymmetric or symmetric not stated); blank where the sources give none'),
    ('maron2021_ip_states', 'Methylation evidence', 'degree called in Maron 2021 Table S2 (PTMScan methyl-arginine immunoprecipitation): me1, me2'),
    ('maron2021_ip_antibodies', 'Methylation evidence', 'antibodies whose immunoprecipitation recovered the site in Maron 2021 Table S3: ADMA (asymmetric dimethyl), SDMA (symmetric dimethyl), MMA (monomethyl)'),
    ('maron2021_GSK591_PRMT5i', 'Methylation evidence', f'response to the type II (PRMT5) inhibitor GSK591 in Maron 2021 Table S3: decreased or increased (adjusted p < {SIGNIFICANT} in at least one immunoprecipitation), mixed, or no significant change'),
    ('maron2021_MS023_typeIi', 'Methylation evidence', 'response to the type I PRMT inhibitor MS023 in Maron 2021 Table S3, as for GSK591'),
    ('maron2021_inhibitor_detail', 'Methylation evidence', 'per inhibitor and immunoprecipitation: log2 intensity difference vs DMSO and adjusted p value of the strongest peptide (lowest adjusted p). Measured on peptides: arginines on the same peptide share its values'),
    ('motifs', 'Context', 'sequence motifs around the site, as defined by scripts/08: RG, RGG, GRG (G at -1 and +1), GAR (G at -2, A at -1), RXR (R at -2 or +2), CARM1_proline_rich (P within -5..+5), CARM1_hydrophobic_proline (P within -5..+5 and a hydrophobic residue within -3..+3), DR, ER (D or E at -1), RD, RE (D or E at +1)'),
    ('disorder', 'Context', f'MobiDB predicted disorder at the site (scripts/06): disordered; disorder boundary (outside, within {CONTEXT_WINDOW} residues of a disordered region); ordered; no disorder annotation'),
    ('disorder_edge_distance', 'Context', 'residues to the nearest edge of a predicted disordered region'),
    ('domain', 'Context', f'InterPro domain or repeat context (scripts/05): in domain; boundary (within {CONTEXT_WINDOW} residues of one); inter-domain linker; distal; no domain annotation; not staged (protein outside the InterPro intervals staged for the published union)'),
    ('nearest_domain', 'Context', 'InterPro entry nearest the site (the innermost when inside several), with its accession and span'),
    ('domain_edge_distance', 'Context', 'residues to the nearest domain edge (0 inside a domain)'),
    ('n_studies', 'Provenance', 'number of studies or databases reporting the site'),
    ('studies', 'Provenance', 'the studies or databases (one row each in the Studies table)'),
    ('study_dois', 'Provenance', 'their DOIs, in the same order'),
    ('n_records', 'Provenance', 'number of source records placed on the site'),
    ('n_records_at_reported_position', 'Provenance', 'records whose own accession and position name this site unchanged'),
    ('source_tables', 'Provenance', 'source tables of the records: Maron 2021 Table S5, ProMetheusDB mmc2, PTM master table'),
    ('record_ids', 'Provenance', 'the supporting rows of the source records table (record_id)'),
    ('release_status', 'Status', 'placed, or included pending a policy decision (see flags)'),
    ('flags', 'Status', 'cautions, separated by ";": ' + '; '.join(FLAGS)),
    ('same_window_in', 'Status', 'other reviewed entries holding the same 31 residues (see the shared_sequence flag)'),
    ('in_previous_version', 'Change', 'whether the previous version (the published union) had this site'),
    ('change_since_previous_version', 'Change', 'kept or added; removed sites are listed in the changes table'),
    ('change_reasons', 'Change', 'for an added site, why (reason codes, see the Reason codes table)'),
]
RECORD_COLUMNS = [
    ('record_id', 'row number of the record, cited by record_ids in the other tables'),
    ('source_table', 'the source table the record comes from'),
    ('study', 'the study or database that reported it'),
    ('study_doi', 'its DOI'),
    ('reference_text', 'the reference as the source table gives it'),
    ('pmid', 'PubMed ID, where the source gives one'),
    ('reported_accession', 'accession as reported'),
    ('reported_gene', 'gene as reported'),
    ('reported_position', 'position as reported'),
    ('reported_methyl_state', 'methylation state as reported'),
    ('evidence_kind', 'kind of sequence evidence: centred window (length), peptide, marked peptide, or position only'),
    ('evidence_sequence', 'the window, peptide or marked peptide, verbatim'),
    ('placed_site', 'where this version places the record (blank: not placed)'),
    ('placed_accession', 'accession of the placed site'),
    ('placed_position', 'position of the placed site'),
    ('placement', 'how the record was placed, in words'),
    ('placement_code', 'the same, as scripts/02 remap_status'),
    ('evidence_check', 'independent check of the placement against the record\'s own evidence (scripts/qc_remap_evidence_check.py)'),
    ('accepted_exception', 'why a contradiction is accepted for now (open issue 9)'),
    ('previous_site', 'where the previous version placed the record'),
    ('change_since_previous_version', 'same, moved, newly placed, no longer placed, or not placed in either'),
    ('reason_code', 'why the placement changed (Reason codes table)'),
    ('issue', 'the issue number of DATA_PROVENANCE_REVIEW.md behind the reason'),
]
CHANGE_COLUMNS = [
    ('site_id', 'the site'),
    ('gene', 'primary gene name'),
    ('uniprot_accession', 'reviewed canonical accession'),
    ('position', 'arginine position'),
    ('in_previous_version', 'site of the previous version'),
    ('in_this_version', 'site of this version'),
    ('change', 'kept, removed or added'),
    ('reason_codes', 'why it was removed or added (Reason codes table)'),
    ('issues', 'issue numbers of DATA_PROVENANCE_REVIEW.md'),
    ('records_now_at_or_previously_at', 'for a removed site, where its records are placed now; for an added site, where they were placed before'),
    ('record_ids', 'the records involved (source records table)'),
]
STUDY_COLUMNS = [
    ('study', 'study or database, as named in the other tables'),
    ('citation', 'short citation'),
    ('doi', 'DOI (linked in the workbook)'),
    ('source_tables', 'source tables its records come from'),
    ('n_records', 'source records'),
    ('n_records_placed', 'source records placed on a site'),
    ('n_sites', 'sites of this version it supports'),
    ('n_sites_only_this_study', 'sites no other study or database reports'),
    ('terms', 'terms beyond citation, where they apply'),
]
REASON_COLUMNS = [
    ('reason_code', 'code used in the records and changes tables'),
    ('issue', 'issue number in DATA_PROVENANCE_REVIEW.md'),
    ('records', 'records with this reason'),
    ('meaning', 'what the previous version did and what this version does'),
]


def md5_bytes(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, low_memory=False)


def clean(value) -> str:
    value = str(value if value is not None else '').strip()
    return '' if value.lower() in ('nan', 'none', '<na>') else value


def study_of(row) -> str:
    if row['source_table'] == 'human_ptm_master':
        return {'unip': 'UniProt', 'iedb': 'IEDB'}.get(row['source'], row['source'] or 'PTM master')
    if row['source_table'] == 'prometheus_mmc2':
        return 'ProMetheusDB'
    ref = row['reference'].strip()
    if ref.lower().startswith('this study'):
        return 'Maron 2021'
    if 'PhosphoSitePlus' in ref:
        return 'PhosphoSitePlus'
    author = ref.split(' et al.')[0].split()[-1] if ' et al.' in ref else ref[:24]
    years = re.findall(r'\((\d{4})\)', ref)
    return f'{author} {years[-1]}' if years else author


def placed_key(row) -> str:
    return f"{row['corrected_canonical_UniProtAC']}:{int(float(row['corrected_position']))}" if row['is_position_resolved'] == 'True' else ''


def r_site(key: str) -> str:
    return key.replace(':', ':R', 1) if key else ''


def evidence_summary(row: pd.Series) -> tuple[str, str]:
    """Kind and text of the record's sequence evidence, as given by the source."""
    for field, kind in (('peptides', 'marked peptide'), ('seq_window', 'window'), ('peptide_sequence', 'peptide')):
        text = clean(row.get(field, ''))
        if text:
            length = int(float(row.get(f'{field}_window_length', 0) or 0)) if field != 'peptides' else 0
            return (f'centred window ({length} aa)' if length else kind), text
    return 'position only', ''


def record_reason(row: pd.Series, old: str, new: str, reference: dict, typed: pd.Series) -> str:
    """Reason code for a record whose placement differs between the previous version and this one."""
    status = row['remap_status']
    evidence = remap.row_evidence(typed)
    windows = [e for e in evidence if e['anchor_type'] in ('centered_window', 'truncated_window')]
    if new and not old:
        if any(re.search(r'^[_X-]|[_X-]$', t.strip()) for f in ('seq_window', 'peptide_sequence') for t in clean(row.get(f, '')).split(';') if t.strip()):
            return 'window_padding'
        # The published code searched the row's accession and gene symbol; a record it could not place that now lands
        # on another entry was reached through the accession history, a gene synonym or the proteome-wide search.
        reported = remap.canonicalize_uniprot_accession(row['substrate_UniProtAC_original'])
        if status == 'sequence_unique_proteome' or clean(row.get('source_accession_history', '')) or new.split(':')[0] != reported:
            return 'outdated_accession_or_symbol'
        return 'other'
    if old and not new:
        if any(e['offsets'] and len(e['sequence']) > e['offsets'][0] >= 0 and e['sequence'][e['offsets'][0]] == 'K' and e['anchor_type'] == 'centered_window' for e in windows):
            return 'lysine_window'
        peptides = [e for e in evidence if e['anchor_type'] == 'peptide']
        if peptides and not windows and all(len(e['offsets']) > 1 for e in peptides) and status in ('unresolved_no_anchor', 'unresolved_evidence_elsewhere'):
            return 'multi_arginine_peptide'
        if status in ('unresolved_evidence_mismatch', 'unresolved_evidence_elsewhere'):
            return 'evidence_mismatch'
        if status == 'unresolved_ambiguous_hit':
            return 'ambiguous_evidence'
        return 'other'
    old_acc, new_acc = old.split(':')[0], new.split(':')[0]
    if old_acc != new_acc:
        return 'unreviewed_accession' if old_acc not in reference else 'wrong_protein'
    if status == 'sequence_position_verified':
        return 'peptide_recentred'
    if row['remap_anchor_type'] == 'truncated_window':
        return 'truncated_window'
    if not windows and any(e['anchor_type'] == 'peptide' for e in evidence):
        return 'peptide_recentred'
    return 'other'


def junction_sites(rows: pd.DataFrame, typed: pd.DataFrame, reference: dict, union_keys: set) -> dict[str, list[int]]:
    """Sites the release includes pending a policy decision: unplaced window records whose window matches one full
    flank (at least 6 residues) and the site exactly at the reported arginine of the reported entry, and matches no other
    reviewed entry exactly (a window crossing a splice junction; a window found in a paralog is issue 9 instead)."""
    proteome = remap.ProteomeIndex(reference)
    out: dict[str, list[int]] = defaultdict(list)
    for i in rows.index[rows['remap_status'] == 'unresolved_evidence_mismatch']:
        r = typed.loc[i]
        acc = remap.canonicalize_uniprot_accession(r['canonical_UniProtAC'])
        seq, pos = reference.get(acc, ''), int(r['position'])
        if not (seq and 0 < pos <= len(seq) and seq[pos - 1] == 'R') or f'{acc}:{pos}' in union_keys:
            continue
        for item in remap.row_evidence(r):
            if item['anchor_type'] not in ('centered_window', 'truncated_window'):
                continue
            t = item['sequence']
            if proteome.find(t):
                continue
            offsets = item['offsets'] if item['offsets'] is not None else [len(t) - item['full_length'] // 2 - 1, item['full_length'] // 2]
            for o in offsets:
                if not (0 <= o < len(t) and t[o] == 'R'):
                    continue
                start = pos - 1 - o
                ok = [0 <= start + k < len(seq) and seq[start + k] == aa for k, aa in enumerate(t)]
                if min(o, len(t) - 1 - o) >= 6 and (all(ok[:o + 1]) or all(ok[o:])):
                    out[f'{acc}:{pos}'].append(i)
                    break
            else:
                continue
            break
    return dict(out)


# ---- Maron 2021 Tables S2 and S3 ---------------------------------------------------------------------------------

def place_maron_peptide(reference: dict, accession: str, peptide: str, begin, position: int) -> str:
    """Site key for a Maron S2/S3 site, or '' when its peptide does not support it on the reference entry: the peptide
    must sit at the reported start, or occur once in the entry (the offset then carries over), and cover an arginine."""
    acc = remap.canonicalize_uniprot_accession(remap.normalize_accession(str(accession).strip()))
    seq, peptide = reference.get(acc, ''), re.sub(r'[^A-Z]', '', str(peptide).upper())
    begin = pd.to_numeric(begin, errors='coerce')
    if not seq or not peptide or pd.isna(begin):
        return ''
    begin, offset = int(begin), position - int(begin)
    if not 0 <= offset < len(peptide) or peptide[offset] != 'R':
        return ''
    if seq[begin - 1:begin - 1 + len(peptide)] != peptide:
        hits = [m.start() for m in re.finditer(f'(?={peptide})', seq)]
        if len(hits) != 1:
            return ''
        begin = hits[0] + 1
    return f'{acc}:{begin + offset}'


def maron_annotations(inputs: Path, reference: dict) -> tuple[dict, dict]:
    """Per site: Table S2 methylation calls, Table S3 antibodies and inhibitor responses. Also returns row counts."""
    calls: dict[str, set] = defaultdict(set)
    s2 = pd.read_excel(inputs / INPUT['maron_s2'], sheet_name=0)
    s2_placed = 0
    for _, r in s2.iterrows():
        mod = re.search(r'me([12])', str(r['Modification']))
        pos = pd.to_numeric(r['Position.In.Protein'], errors='coerce')
        key = place_maron_peptide(reference, r['Accession'], r['Sequence'], r['Begin.Pos'], int(pos)) if mod and pd.notna(pos) else ''
        if key:
            calls[key].add(f'me{mod.group(1)}')
            s2_placed += 1
    antibodies: dict[str, set] = defaultdict(set)
    response: dict[tuple, list] = defaultdict(list)
    s3_rows = s3_placed = 0
    for sheet, drug in (('Control_vs_GSK591', 'GSK591'), ('Control_vs_MS023', 'MS023')):
        df = pd.read_excel(inputs / INPUT['maron_s3'], sheet_name=sheet)
        for _, r in df.iterrows():
            for pos, degree in re.findall(r'R(\d+)me([12])', str(r['Modifications'])):
                s3_rows += 1
                key = place_maron_peptide(reference, r['Accession'], r['Sequence'], r['Begin.Pos'], int(pos))
                if not key:
                    continue
                s3_placed += 1
                ip = clean(r['IP'])
                antibodies[key].add(ip)
                diff, padj = pd.to_numeric(r[f'{drug}_vs_DMSO_diff'], errors='coerce'), pd.to_numeric(r[f'{drug}_vs_DMSO_p.adj'], errors='coerce')
                if pd.notna(diff) and pd.notna(padj):
                    response[(key, drug)].append((float(padj), float(diff), ip))
    out: dict[str, dict] = defaultdict(dict)
    for key, states in calls.items():
        out[key]['maron2021_ip_states'] = '; '.join(sorted(states))
    for key, ips in antibodies.items():
        out[key]['maron2021_ip_antibodies'] = '; '.join(sorted(ips))
    details: dict[str, list] = defaultdict(list)
    for (key, drug), values in sorted(response.items()):
        down = any(p < SIGNIFICANT and d < 0 for p, d, _ in values)
        up = any(p < SIGNIFICANT and d > 0 for p, d, _ in values)
        column = 'maron2021_GSK591_PRMT5i' if drug == 'GSK591' else 'maron2021_MS023_typeIi'
        out[key][column] = 'mixed' if down and up else 'decreased' if down else 'increased' if up else 'no significant change'
        best = {}
        for p, d, ip in sorted(values):
            best.setdefault(ip, (d, p))
        details[key].append(f'{drug}: ' + ', '.join(f'{ip} IP {d:+.2f} (adj. p {p:.2g})' for ip, (d, p) in sorted(best.items())))
    for key, parts in details.items():
        out[key]['maron2021_inhibitor_detail'] = ' | '.join(parts)
    counts = {'maron_s2_rows': len(s2), 'maron_s2_rows_placed': s2_placed, 'maron_s3_site_rows': s3_rows, 'maron_s3_site_rows_placed': s3_placed}
    return dict(out), counts


# ---- Structural context and motifs (the pipeline's own classifiers) ----------------------------------------------

def context_annotations(sites: pd.DataFrame, inputs: Path, staged_for_interpro: set) -> pd.DataFrame:
    frame = sites[['uniprot_accession', 'position']].rename(columns={'uniprot_accession': 'canonical_UniProtAC', 'position': 'corrected_position'})
    disorder = disorder_context.normalize_disorder_intervals(pd.read_csv(inputs / INPUT['disorder'], sep='\t', low_memory=False))
    by_acc = {a: g.reset_index(drop=True) for a, g in disorder.groupby('canonical_UniProtAC')}
    dis = disorder_context.annotate_frame(frame, by_acc, 'canonical_UniProtAC', 'corrected_position', CONTEXT_WINDOW)

    intervals = pd.read_csv(inputs / INPUT['domains'], sep='\t', low_memory=False)
    intervals = intervals[intervals['interpro_type'].astype(str).isin({'domain', 'repeat'})].copy()  # scripts/05 default types
    for col in ('fragment_start', 'fragment_end'):
        intervals[col] = pd.to_numeric(intervals[col], errors='coerce')
    intervals = intervals.dropna(subset=['fragment_start', 'fragment_end'])
    intervals[['fragment_start', 'fragment_end']] = intervals[['fragment_start', 'fragment_end']].astype(int)
    domains_by_acc = {a: g.sort_values(['fragment_start', 'fragment_end', 'interpro_accession']).reset_index(drop=True)
                      for a, g in intervals.groupby('canonical_UniProtAC')}
    dom = pd.DataFrame([domain_context.classify_row(r, domains_by_acc, 'canonical_UniProtAC', 'corrected_position', CONTEXT_WINDOW)
                        for _, r in frame.iterrows()], index=frame.index)

    words = {'disordered': 'disordered', 'disorder_boundary': 'disorder boundary', 'ordered': 'ordered',
             'no_disorder_annotation': 'no disorder annotation'}
    domain_words = {'in_domain': 'in domain', 'boundary': 'boundary', 'inter_domain_linker': 'inter-domain linker',
                    'distal': 'distal', 'no_domain_annotation': 'no domain annotation'}
    out = pd.DataFrame(index=sites.index)
    out['motifs'] = [', '.join(k.removeprefix('motif_') for k, v in motif_odds.family_flags(w).items() if v) for w in sites['window_15']]
    out['disorder'] = dis['disorder_context_class'].map(words)
    out['disorder_edge_distance'] = pd.to_numeric(dis['nearest_disorder_edge_distance'], errors='coerce').astype('Int64')
    staged = sites['uniprot_accession'].isin(staged_for_interpro)
    out['domain'] = dom['domain_context_class'].map(domain_words).where(staged, 'not staged')
    out['nearest_domain'] = [f"{n} ({a}, {s}-{e})" if st and a else '' for n, a, s, e, st in
                             zip(dom['nearest_interpro_name'], dom['nearest_interpro_accession'], dom['nearest_domain_start'],
                                 dom['nearest_domain_end'], staged)]
    out['domain_edge_distance'] = pd.to_numeric(dom['nearest_domain_edge_distance'], errors='coerce').astype('Int64').where(staged)
    return out


# ---- Build -------------------------------------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--release-root', required=True, help='results root of the corrected run (remapped/, integrated/, qc/)')
    ap.add_argument('--published-root', default='PTM_results', help='results root of the previous version')
    ap.add_argument('--inputs', default='inputs', help='frozen inputs with MANIFEST.tsv')
    ap.add_argument('--allow', default='config/remap_evidence_allowlist.tsv')
    ap.add_argument('--version', required=True, help='version of this release, e.g. v2.0')
    ap.add_argument('--previous-version', default='v1.0', help='version of the published union')
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    root, pub_root, inputs, outdir = Path(args.release_root), Path(args.published_root), Path(args.inputs), Path(args.outdir)
    version, previous = args.version, args.previous_version
    outdir.mkdir(parents=True, exist_ok=True)

    reference = remap.parse_fasta(inputs / INPUT['fasta'])
    meta = read_tsv(inputs / INPUT['entries']).set_index('Entry')
    rows = read_tsv(root / 'remapped' / 'human_arg_methyl_source_rows_with_remap.tsv')
    typed = pd.read_csv(root / 'remapped' / 'human_arg_methyl_source_rows_with_remap.tsv', sep='\t', low_memory=False)
    typed = remap.annotate_window_lengths(typed) if 'seq_window_window_length' not in typed else typed
    pub = read_tsv(pub_root / 'remapped' / 'human_arg_methyl_source_rows_with_remap.tsv')
    check = read_tsv(root / 'qc' / 'remap_evidence_check' / 'remap_evidence_check_rows.tsv')
    union = read_tsv(root / 'integrated' / 'human_arg_methyl_union_dedup_by_site.tsv')
    pub_union = read_tsv(pub_root / 'integrated' / 'human_arg_methyl_union_dedup_by_site.tsv')
    allow = (pd.read_csv(args.allow, sep='\t', dtype=str, comment='#').fillna('') if Path(args.allow).exists()
             else pd.DataFrame(columns=['source_table', 'reported_accession', 'reported_position', 'placed_site', 'reason']))
    allowed = {(r.source_table, r.reported_accession, r.reported_position, r.placed_site): r.reason for r in allow.itertuples()}

    # This version and the previous one read the same source records in the same order.
    key_cols = ['source_table', 'substrate_UniProtAC_original', 'position', 'peptide_sequence', 'seq_window', 'peptides']
    if len(rows) != len(pub) or not (rows[key_cols].values == pub[key_cols].values).all():
        raise SystemExit('the release run and the published run do not hold the same source records in the same order')
    if len(check) != len(rows):
        raise SystemExit('the evidence check does not cover every record of the release run')

    def entry(acc: str, field: str) -> str:
        return clean(meta[field].get(acc, '')) if field in meta else ''

    # ---- Records -----------------------------------------------------------------------------------------------
    rec = []
    for i, r in rows.iterrows():
        study = study_of(r)
        kind, text = evidence_summary(typed.loc[i])
        new, old = placed_key(r), placed_key(pub.loc[i])
        exception = allowed.get((r['source_table'], check.loc[i, 'reported_accession'], check.loc[i, 'reported_position'], check.loc[i, 'placed_site']), '')
        change = 'same' if new == old and new else 'not placed in either' if not new and not old else 'moved' if new and old else 'newly placed' if new else 'no longer placed'
        reason = '' if change in ('same', 'not placed in either') else record_reason(r, old, new, reference, typed.loc[i])
        rec.append({
            'record_id': i + 1,
            'source_table': SOURCE_TABLE.get(r['source_table'], r['source_table']),
            'study': study,
            'study_doi': STUDIES.get(study, ('', ''))[0],
            'reference_text': clean(r['reference']),
            'pmid': clean(r['pmid']),
            'reported_accession': clean(r['substrate_UniProtAC_original']),
            'reported_gene': clean(r['substrate_genename']),
            'reported_position': clean(r['position']),
            'reported_methyl_state': clean(r['arg_methyl_state']),
            'evidence_kind': kind,
            'evidence_sequence': text,
            'placed_site': r_site(new),
            'placed_accession': new.split(':')[0] if new else '',
            'placed_position': new.split(':')[1] if new else '',
            'placement': PLACEMENT.get(r['remap_status'], r['remap_status']),
            'placement_code': r['remap_status'],
            'evidence_check': check.loc[i, 'verdict'],
            'accepted_exception': exception,
            'previous_site': r_site(old),
            'change_since_previous_version': change,
            'reason_code': reason,
            'issue': REASONS[reason][0] if reason else '',
        })
    records = pd.DataFrame(rec)

    # ---- Sites -------------------------------------------------------------------------------------------------
    union_keys = {f'{a}:{int(float(p))}' for a, p in zip(union['canonical_UniProtAC'], union['corrected_position'])}
    pending = junction_sites(rows, typed, reference, union_keys)
    site_of = pd.Series([s.replace(':R', ':') for s in records['placed_site']], index=records.index)
    for site, ids in pending.items():
        site_of.loc[ids] = site
    pub_keys = {f'{a}:{int(float(p))}' for a, p in zip(pub_union['canonical_UniProtAC'], pub_union['corrected_position'])}
    observed = {f'{a}:{int(float(p))}': o for a, p, o in zip(union['canonical_UniProtAC'], union['corrected_position'], union['arg_methyl_states_observed'])}

    # A site whose surrounding sequence occurs identically in another reviewed entry (histone variants, identical
    # paralogs such as HSPA1A/HSPA1B) cannot be told apart by its peptide or window: the gene is the one the source named.
    proteome = remap.ProteomeIndex(reference)

    def same_window_in(acc: str, pos: int) -> list[str]:
        window = reference.get(acc, '')[max(0, pos - 16):pos + 15]
        if len(window) < 15:
            return []
        return [f"{entry(a, 'Gene Names (primary)') or a} ({a})" for a in sorted({a for a, _ in proteome.find(window) if a != acc})]

    sites = []
    for site, group in records.groupby(site_of, sort=False):
        if not site:
            continue
        acc, pos = site.split(':')[0], int(site.split(':')[1])
        seq = reference.get(acc, '')
        studies = sorted(set(group['study']))
        flags = []
        if site in pending:
            flags.append(PENDING_JUNCTION)
        if (group['accepted_exception'] != '').any():
            flags.append('paralog_window')
        if site not in pending and group['placement_code'].str.contains('near_match').all():
            flags.append('near_match_only')
        shared = same_window_in(acc, pos)
        if shared:
            flags.append('shared_sequence')
        if len(studies) == 1:
            flags.append('single_study')
        window_31 = ''.join(seq[k] if 0 <= k < len(seq) else '_' for k in range(pos - 16, pos + 15))
        added = site not in pub_keys
        sites.append({
            'site_id': r_site(site),
            'gene': entry(acc, 'Gene Names (primary)'),
            'protein_name': entry(acc, 'Protein names').split(' (')[0],
            'uniprot_accession': acc,
            'uniprot_entry_name': entry(acc, 'Entry Name'),
            'position': pos,
            'residue': seq[pos - 1] if seq and pos <= len(seq) else '',
            'window_15': window_31[8:23],
            'window_31': window_31,
            'protein_length': len(seq),
            'reference_sequence_md5': hashlib.md5(seq.encode()).hexdigest().upper() if seq else '',
            'methyl_states_reported': '; '.join(s.removeprefix('R') for s in clean(observed.get(site, '')).split(';') if s),
            'n_studies': len(studies),
            'studies': '; '.join(studies),
            'study_dois': '; '.join(STUDIES[s][0] for s in studies if s in STUDIES),
            'n_records': len(group),
            'n_records_at_reported_position': int(((group['reported_accession'].map(remap.canonicalize_uniprot_accession) == acc) & (group['reported_position'] == str(pos))).sum()),
            'source_tables': '; '.join(sorted(set(group['source_table']))),
            'record_ids': ';'.join(map(str, group['record_id'])),
            'release_status': 'pending policy decision: not placed by the pipeline' if site in pending else 'placed',
            'flags': ';'.join(flags),
            'same_window_in': '; '.join(shared[:10]) + (f'; and {len(shared) - 10} more' if len(shared) > 10 else ''),
            'in_previous_version': not added,
            'change_since_previous_version': 'added' if added else 'kept',
            'change_reasons': ';'.join(([PENDING_JUNCTION] if site in pending else []) + sorted(set(group.loc[group['previous_site'] != r_site(site), 'reason_code']) - {''})) if added else '',
        })
    sites = pd.DataFrame(sites)
    sites = sites.assign(_no_gene=sites['gene'] == '').sort_values(['_no_gene', 'gene', 'uniprot_accession', 'position'])
    sites = sites.drop(columns='_no_gene').reset_index(drop=True)

    maron, maron_counts = maron_annotations(inputs, reference)
    for column in ('maron2021_ip_states', 'maron2021_ip_antibodies', 'maron2021_GSK591_PRMT5i', 'maron2021_MS023_typeIi', 'maron2021_inhibitor_detail'):
        sites[column] = [maron.get(s.replace(':R', ':'), {}).get(column, '') for s in sites['site_id']]
    staged_for_interpro = set(pub_union['canonical_UniProtAC'])
    sites = pd.concat([sites, context_annotations(sites, inputs, staged_for_interpro)], axis=1)
    sites = sites[[c for c, _, _ in SITE_COLUMNS]]

    # ---- Changes since the previous version --------------------------------------------------------------------
    release_keys = {s.replace(':R', ':') for s in sites['site_id']}
    changes = []
    for site in sorted(pub_keys | release_keys, key=lambda k: (k.split(':')[0], int(k.split(':')[1]))):
        in_pub, in_rel = site in pub_keys, site in release_keys
        acc, pos = site.split(':')
        if in_pub and in_rel:
            change, involved = 'kept', records.iloc[0:0]
        elif in_pub:
            change, involved = 'removed', records[(records['previous_site'] == r_site(site)) & (site_of != site)]
        else:
            change, involved = 'added', records[(site_of == site) & (records['previous_site'] != r_site(site))]
        codes = sorted(set(involved['reason_code']) - {''})
        if change == 'added' and site in pending:
            codes = [PENDING_JUNCTION] + codes
        elsewhere = sorted(set((site_of[involved.index].map(r_site) if change == 'removed' else involved['previous_site'])) - {'', r_site(site)})
        changes.append({
            'site_id': r_site(site),
            'gene': entry(acc, 'Gene Names (primary)'),
            'uniprot_accession': acc,
            'position': int(pos),
            'in_previous_version': in_pub,
            'in_this_version': in_rel,
            'change': change,
            'reason_codes': ';'.join(codes),
            'issues': ';'.join(sorted({REASONS[c][0] for c in codes if c in REASONS and REASONS[c][0]}, key=int)),
            'records_now_at_or_previously_at': '; '.join(elsewhere),
            'record_ids': ';'.join(map(str, involved['record_id'])),
        })
    changes = pd.DataFrame(changes)

    # ---- Studies -----------------------------------------------------------------------------------------------
    site_studies = sites.set_index('site_id')['studies'].str.split('; ')
    study_rows = []
    for study, group in records.groupby('study'):
        supported = site_studies[site_studies.map(lambda s: study in s)]
        doi, citation = STUDIES.get(study, ('', ''))
        texts = group['reference_text'].replace('', pd.NA).dropna()
        study_rows.append({
            'study': study,
            'citation': citation or (texts.mode().iloc[0] if len(texts) else ''),
            'doi': doi,
            'source_tables': '; '.join(sorted(set(group['source_table']))),
            'n_records': len(group),
            'n_records_placed': int((group['placed_site'] != '').sum()),
            'n_sites': len(supported),
            'n_sites_only_this_study': int((supported.map(len) == 1).sum()),
            'terms': PSP_TERMS if study == 'PhosphoSitePlus' else '',
        })
    studies = pd.DataFrame(study_rows).sort_values('n_sites', ascending=False).reset_index(drop=True)
    reason_counts = records['reason_code'].replace('', pd.NA).dropna().value_counts()
    reasons = pd.DataFrame([{'reason_code': c, 'issue': i, 'records': int(reason_counts.get(c, 0)), 'meaning': t}
                            for c, (i, t) in REASONS.items()])
    records = records[[c for c, _ in RECORD_COLUMNS]]

    # ---- Write -------------------------------------------------------------------------------------------------
    stem = 'human_methylarginine'
    files = {
        'sites': outdir / f'{stem}_sites_{version}.tsv',
        'xlsx': outdir / f'{stem}_sites_{version}.xlsx',
        'records': outdir / f'{stem}_source_records_{version}.tsv',
        'changes': outdir / f'{stem}_changes_{previous}_to_{version}.tsv',
        'studies': outdir / f'{stem}_studies_{version}.tsv',
    }
    for p in outdir.iterdir():
        if p.is_file():
            p.unlink()
    sites.to_csv(files['sites'], sep='\t', index=False)
    records.to_csv(files['records'], sep='\t', index=False)
    changes.to_csv(files['changes'], sep='\t', index=False)
    studies.to_csv(files['studies'], sep='\t', index=False)
    counts = summary_counts(sites, records, changes, pending)
    counts.update(maron_counts)
    build = build_record(root, inputs)
    write_readme(outdir, version, previous, files, counts, sites, reasons, build)
    import release_workbook
    release_workbook.write(files['xlsx'], version, previous, counts, build, sites, records, changes, studies, reasons,
                           {'Sites': SITE_COLUMNS, 'Source records': RECORD_COLUMNS, 'Changes': CHANGE_COLUMNS,
                            'Studies': STUDY_COLUMNS, 'Reason codes': REASON_COLUMNS}, FLAGS, SOURCE_TABLE_ABOUT)
    with open(outdir / 'MD5SUMS', 'w') as fh:
        for p in sorted(outdir.iterdir()):
            if p.is_file() and p.name != 'MD5SUMS':
                fh.write(f'{md5_bytes(p)}  {p.name}\n')
    print(counts)


def summary_counts(sites, records, changes, pending) -> dict:
    placed = sites[sites['release_status'] == 'placed']
    return {
        'sites_placed': len(placed), 'proteins_placed': placed['uniprot_accession'].nunique(),
        'sites_pending': len(pending), 'records': len(records),
        'records_placed': int((records['placed_site'] != '').sum()),
        'kept': int((changes['change'] == 'kept').sum()), 'removed': int((changes['change'] == 'removed').sum()),
        'added': int((changes['change'] == 'added').sum()),
        'flags': {k: int(v) for k, v in sites['flags'].str.split(';').explode().value_counts().items() if k},
        'with_maron_ip': int((sites['maron2021_ip_antibodies'] != '').sum()),
        'with_state': int(((sites['methyl_states_reported'] != '') | (sites['maron2021_ip_states'] != '')).sum()),
        'not_staged_interpro': int((sites['domain'] == 'not staged').sum()),
        'not_staged_proteins': int(sites.loc[sites['domain'] == 'not staged', 'uniprot_accession'].nunique()),
        'studies': int(records['study'].nunique()),
    }


def build_record(root: Path, inputs: Path) -> dict:
    repo = SCRIPTS.parent
    try:
        commit = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True, text=True, cwd=repo).stdout.strip()
    except OSError:
        commit = ''
    run_info = (root / 'RUN_INFO.txt').read_text().replace(str(repo) + '/', '') if (root / 'RUN_INFO.txt').exists() else ''
    shown = Path(os.path.relpath(inputs.resolve(), repo)) if inputs.resolve().is_relative_to(repo) else inputs
    return {'commit': commit or 'unknown', 'inputs': str(shown), 'run_info': run_info.strip()}


def write_readme(outdir: Path, version: str, previous: str, files: dict, counts: dict, sites: pd.DataFrame,
                 reasons: pd.DataFrame, build: dict) -> None:
    n = {k: f'{v:,}' if isinstance(v, int) else v for k, v in counts.items()}
    lines = [
        f'# Human methylarginine sites, {version}',
        '',
        'The corrected union of human protein arginine methylation sites behind DeAngelo, Silverstein & Shechter, '
        'Cell Chem Biol 33:1404-1421 (2026), doi:10.1016/j.chembiol.2026.09.011, with the provenance of every site and of '
        f'every source record. It replaces {previous}, the published union (PTM_results/, tag v1.0.0: 13,719 sites on '
        '5,090 accessions), after the remapping corrections described in DATA_PROVENANCE_REVIEW.md.',
        '',
        f"- {n['sites_placed']} sites on {n['proteins_placed']} reviewed human proteins, plus {n['sites_pending']} sites "
        f'included pending a policy decision (flag {PENDING_JUNCTION}).',
        f"- {n['records']} source records from {n['studies']} studies and databases, {n['records_placed']} of them placed.",
        f"- Since {previous}: {n['kept']} sites kept, {n['removed']} removed, {n['added']} added.",
        '',
        '## Files',
        '',
        '| File | One row per | Use it for |',
        '|---|---|---|',
        f"| {files['sites'].name} | site | the site list: sequence, methylation evidence, context, provenance |",
        f"| {files['xlsx'].name} | (workbook) | all tables, with column notes, links and highlighting |",
        f"| {files['records'].name} | source record (one study's report of one site) | what each study reported and where it was placed |",
        f"| {files['changes'].name} | site of {previous} or {version} | what changed and why |",
        f"| {files['studies'].name} | study or database | citations, DOIs and how much each contributes |",
        '| MD5SUMS | file | checksums (`md5sum -c MD5SUMS`) |',
        '',
        '## Where the sites come from',
        '',
        'Three source tables: Maron et al. 2021 Table S5 (a compilation of ten studies: Larsen 2016, PhosphoSitePlus, '
        'Musiani 2019, Guo 2014, Fong 2019, Maron 2021\'s own data, Li 2021, Fedoriw 2019, Lim 2020, Wei 2020), ProMetheusDB '
        '(Massignani et al. 2022) and the arginine methylation rows of the PTM master table (UniProt and IEDB records). '
        'Each record\'s peptide, window or position was placed on the reviewed canonical UniProt sequence by '
        'scripts/02_remap_arg_methyl_to_canonical.py and checked against its own evidence by scripts/qc_remap_evidence_check.py. '
        'Maron 2021 Tables S2 and S3 add methylation-degree calls, the antibodies that recovered each site and its response '
        'to PRMT inhibitors; they annotate sites but add none. A Table S2 or S3 site is attached only when its peptide sits '
        'on the reference entry at the reported start, or occurs once in the entry '
        f"({n['maron_s2_rows_placed']} of {n['maron_s2_rows']} Table S2 rows, {n['maron_s3_site_rows_placed']} of "
        f"{n['maron_s3_site_rows']} Table S3 site rows).",
        '',
        '## Site table columns',
        '',
    ]
    group = ''
    for column, g, text in SITE_COLUMNS:
        if g != group:
            lines += ['', f'**{g}**', ''] if group else [f'**{g}**', '']
            group = g
        lines.append(f'- `{column}`: {text}.')
    lines += ['', '## Flags', '']
    lines += [f"- `{k}` ({counts['flags'].get(k, 0):,} sites): {v}." for k, v in FLAGS.items()]
    lines += [
        '',
        f"Methylation degree is known for {n['with_state']} sites (sources or Maron Table S2); {n['with_maron_ip']} sites "
        'have Maron 2021 antibody and inhibitor data. Domain context is "not staged" for '
        f"{n['not_staged_interpro']} sites on {n['not_staged_proteins']} proteins that were not in the published union, "
        'because the frozen InterPro intervals were staged for that union only.',
        '',
        '## Source records, changes and studies',
        '',
        'The columns of the other tables are listed on the workbook\'s Columns sheet. Every record carries the code of '
        'the reason its placement changed; every removed or added site carries the codes of its records:',
        '',
        '| reason_code | issue | records | meaning |',
        '|---|---|---|---|',
    ]
    lines += [f'| {r.reason_code} | {r.issue} | {r.records:,} | {r.meaning} |' for r in reasons.itertuples()]
    lines += [
        '',
        '## How it was built',
        '',
        f"- Code: PTM_Pipeline commit {build['commit']}; driver hpc/build_release.sh (parse, place, union, the evidence "
        'check and the residue audit, which stop the build on any failure), then scripts/build_release.py.',
        f"- Inputs: {build['inputs']}/, checked against {build['inputs']}/MANIFEST.tsv. Reference sequences: UniProtKB "
        'reviewed human (proteome UP000005640), staged 2026-04-07; release not recorded, earlier than 2026_03. The '
        'published union was mapped to the same sequences. Accession history: UniProt 2026_03.',
        '- Context: MobiDB disorder (export of 2026-03-30), InterPro domains and repeats (staged 2026-04-07), motifs, '
        'with the classifiers of scripts/06, scripts/05 and scripts/08.',
        '',
        '## Licences',
        '',
        'PhosphoSitePlus® data (www.phosphosite.org) are licensed under CC BY-NC-SA 3.0: non-commercial use, with '
        'attribution, adaptations shared under the same licence. The release as a whole is CC BY-NC-SA 4.0 (LICENSE.md). '
        'Cite the studies in study_dois for the sites you use, and the paper above for the union.',
    ]
    if build['run_info']:
        lines += ['', '## Run record', '', '```', build['run_info'], '```']
    (outdir / 'README.md').write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()

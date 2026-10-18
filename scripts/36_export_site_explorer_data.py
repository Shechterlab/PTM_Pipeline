from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from common import CPLM4_SOURCE_BITS, cplm4_evidence_mask


SOURCES = ['Maron2021', 'ProMetheusDB', 'human_ptm_master']
STATES = ['', 'Rme1', 'Rme2', 'mixed']
MOTIFS = [
    'RG', 'RGG', 'GRG', 'GAR', 'RXR',
    'CARM1_proline_rich', 'CARM1_hydrophobic_proline',
    'DR', 'RD', 'ER', 'RE',
]
MARON_IPS = ['MMA', 'ADMA', 'SDMA']
MARON_TREATMENTS = ['GSK591', 'MS023']
OTHER_PTM_GROUPS = ['Phosphorylation', 'Lys acylation', 'Ubiquitin/SUMO', 'Lys methylation']
# Bit k of an other-PTM entry's source mask (its fourth element) means OTHER_PTM_SOURCES[k].
OTHER_PTM_SOURCES = ['Cross-PTM master table', 'CPLM 4.0 (Exp.)', 'CPLM 4.0 (Dat.)']
ALS_FAMILY = 'ALS / motor neuron / FTD'
WINDOW_FLANK = 15


def read_tsv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, **kwargs)


def read_optional(path: Path, **kwargs) -> pd.DataFrame:
    if not path.exists():
        print(f'[site explorer] optional input missing, skipped: {path}')
        return pd.DataFrame()
    return read_tsv(path, **kwargs)


def to_float(value, digits: int = 3) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if pd.isna(number):
        return None
    return round(number, digits)


def to_sig(value, digits: int = 3) -> float | None:
    """Round to significant digits; used for p and q values."""
    number = to_float(value, digits=300)
    if number is None:
        return None
    return float(f'{number:.{digits}g}')


def to_int(value) -> int | None:
    number = to_float(value, digits=6)
    return None if number is None else int(number)


def split_multi(value: str, sep: str = ';') -> list[str]:
    return [part.strip() for part in str(value or '').split(sep) if part.strip()]


class Vocab:
    def __init__(self, values: list[str] | None = None):
        self.values: list[str] = []
        self.lookup: dict[str, int] = {}
        for value in values or []:
            self.index(value)

    def index(self, value: str) -> int:
        if value not in self.lookup:
            self.lookup[value] = len(self.values)
            self.values.append(value)
        return self.lookup[value]

    def mask(self, labels: list[str]) -> int:
        mask = 0
        for label in labels:
            mask |= 1 << self.index(label)
        return mask


def recommended_name(value: str) -> str:
    """First name from a UniProt-style 'Name (Alt) (Alt) [Cleaved into: ...]' or 'Name; ; Name (...)' string."""
    name = split_multi(value)[0] if split_multi(value) else ''
    name = name.split(' [Cleaved into')[0].split(' [Includes')[0]
    depth = 0
    for i, char in enumerate(name):
        if char == '(' and depth == 0 and i > 0 and name[i - 1] == ' ':
            return name[:i - 1].strip()
        depth += char == '('
        depth -= char == ')'
    return name.strip()


# Typos carried over from the Maron et al. 2021 Table S5 reference strings.
REFERENCE_FIXES = {'Cancer Cell 36, 194-209.e199 (2019)': 'Cancer Cell 36, 194-209.e9 (2019)'}


def clean_reference(reference: str) -> str:
    for bad, good in REFERENCE_FIXES.items():
        reference = reference.replace(bad, good)
    return reference


def study_label(reference: str) -> str:
    """Short 'Author Year' label for a Maron Table S5 reference string."""
    ref = reference.strip()
    if ref.lower().startswith('this study'):
        return 'Maron 2021'
    if 'PhosphoSitePlus' in ref:
        return 'PhosphoSitePlus'
    author = ref.split(' et al.')[0].split()[-1] if ' et al.' in ref else ref[:24]
    years = re.findall(r'\((\d{4})\)', ref)
    return f'{author} {years[-1]}' if years else author


PTM_MASTER_SOURCES = {'unip': 'UniProt', 'iedb': 'IEDB'}


def gene_symbols(value: str) -> list[str]:
    """Symbols from a plain symbol or a UniProt 'Name: X; Synonyms: A, B; ORFNames: ...;' gene string."""
    value = str(value or '').strip()
    if 'Name=' in value or 'Name:' in value or 'Synonyms' in value:
        out = []
        for part in value.split(';'):
            key, _, rest = part.partition(':') if ':' in part else part.partition('=')
            if key.strip() in ('Name', 'Synonyms'):
                out += [g.strip() for g in rest.split(',') if g.strip()]
        return out
    return [g for g in value.replace(';', ' ').split() if g]


def load_uniprot(path: Path | None) -> dict[str, dict]:
    """Per-accession length, names and gene symbols from a UniProt metadata TSV.

    Uses the pipeline's PTM_UNIPROT_METADATA columns (Entry, Entry Name, Gene Names (primary),
    Protein names, Length); an optional 'Gene Names' column adds synonyms for search.
    """
    if path is None or not path.exists():
        return {}
    meta = read_tsv(path)
    out = {}
    for row in meta.to_dict('records'):
        primary = split_multi(row.get('Gene Names (primary)', ''))
        out[row['Entry']] = {
            'length': to_int(row.get('Length')),
            'name': row.get('Protein names', ''),
            'entry_name': row.get('Entry Name', ''),
            'gene': primary[0] if primary else '',
            'synonyms': [g for g in str(row.get('Gene Names', '')).replace(';', ' ').split() if g],
        }
    return out


def enrichment_rows(df: pd.DataFrame, label_col: str = 'category') -> list[list]:
    rows = []
    for _, row in df.iterrows():
        rows.append([
            row[label_col],
            to_int(row.get('target_count')),
            to_float(row.get('target_fraction'), 4),
            to_float(row.get('background_fraction'), 4),
            to_float(row.get('log2_odds_ratio')),
            to_sig(row.get('q_value')),
        ])
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(
        description='Export the integrated methylarginine union and its site/protein annotations as compact JSON '
                    'for the static site explorer in docs/.',
    )
    ap.add_argument('--results-root', default=os.environ.get('PTM_RESULTS_ROOT', 'PTM_results'))
    ap.add_argument('--uniprot-metadata', default=os.environ.get('PTM_UNIPROT_METADATA', ''),
                    help='UniProt metadata TSV with Entry, Entry Name, Protein names and Length (PTM_UNIPROT_METADATA).')
    ap.add_argument('--cplm4-sites', default=os.environ.get('PTM_CPLM4_SITES', ''),
                    help='Validated CPLM 4.0 lysine methylation table from 39_stage_cplm4_lys_methylation.py '
                         '(default <results-root>/cplm4/cplm4_human_lys_methyl_validated.tsv; skipped if missing).')
    ap.add_argument('--outdir', default='docs/data')
    args = ap.parse_args()

    root = Path(args.results_root)
    meta_path = Path(args.uniprot_metadata) if args.uniprot_metadata else None
    uniprot = load_uniprot(meta_path)

    # ---- Site-level tables -------------------------------------------------------------------
    sites = read_tsv(root / 'integrated' / 'human_arg_methyl_union_dedup_by_site.tsv')
    disorder = read_tsv(root / 'disorder_context' / 'sites_with_disorder_context.tsv')
    domain = read_tsv(root / 'domain_context' / 'sites_with_domain_context.tsv')
    windows = read_tsv(root / 'motif_arg_odds' / 'sequence_windows.tsv')

    merged = (
        sites
        .merge(disorder[['site_key', 'disorder_context_class', 'nearest_disorder_edge_distance']], on='site_key', how='left')
        .merge(domain[['site_key', 'domain_context_class', 'nearest_domain_edge_distance', 'nearest_domain_side',
                       'nearest_interpro_accession', 'nearest_interpro_name', 'nearest_domain_class']],
               on='site_key', how='left')
        .merge(windows[['site_key', 'window_31'] + [f'motif_{m}' for m in MOTIFS]], on='site_key', how='left')
        .fillna('')
    )
    merged['position_int'] = pd.to_numeric(merged['corrected_position'], errors='coerce').astype(int)
    merged = merged.sort_values(['canonical_UniProtAC', 'position_int']).reset_index(drop=True)
    site_index = {(acc, pos): i for i, (acc, pos) in enumerate(zip(merged['canonical_UniProtAC'], merged['position_int']))}

    def site_lookup(acc: str, pos) -> int | None:
        pos_int = to_int(pos)
        return None if pos_int is None else site_index.get((acc, pos_int))

    # Maron Table S2 (me1/me2 calls) and S3 (methyl-IP detection and PRMT inhibitor response).
    maron_root = root / 'maron_state_treatment'
    s2_state = defaultdict(int)
    for _, row in read_optional(maron_root / 'maron_table_s2_site_state_collapsed.tsv').iterrows():
        idx = site_lookup(row['canonical_UniProtAC'], row['position'])
        if idx is not None and row['methyl_state'] in ('me1', 'me2'):
            s2_state[idx] |= 1 if row['methyl_state'] == 'me1' else 2
    ip_mask = defaultdict(int)
    for _, row in read_optional(maron_root / 'maron_table_s3_site_ip_collapsed.tsv').iterrows():
        idx = site_lookup(row['canonical_UniProtAC'], row['position'])
        if idx is not None and row['ip'] in MARON_IPS:
            ip_mask[idx] |= 1 << MARON_IPS.index(row['ip'])
    treatment = defaultdict(list)
    maron_unmatched = 0
    for _, row in read_optional(maron_root / 'maron_table_s3_site_treatment_collapsed.tsv').iterrows():
        idx = site_lookup(row['canonical_UniProtAC'], row['position'])
        if idx is None:
            maron_unmatched += 1
            continue
        if row['ip'] in MARON_IPS and row['treatment'] in MARON_TREATMENTS:
            treatment[idx].append([
                MARON_IPS.index(row['ip']),
                MARON_TREATMENTS.index(row['treatment']),
                to_float(row['treatment_diff_median']),
                to_sig(row['treatment_padj_min']),
                to_int(row['peptide_rows']),
            ])

    # ---- Protein-level tables ----------------------------------------------------------------
    functional = read_tsv(root / 'functional_class_union' / 'protein_functional_labels.tsv')
    functional = functional.drop_duplicates('canonical_UniProtAC').set_index('canonical_UniProtAC')
    arg_odds = read_tsv(root / 'motif_arg_odds' / 'per_protein_arg_odds.tsv').drop_duplicates('canonical_UniProtAC')
    arg_odds = arg_odds.set_index('canonical_UniProtAC')
    clustering_root = root / 'methyl_arg_clustering'
    local_density = read_optional(clustering_root / 'protein_local_density_summary.tsv')
    local_density = local_density.drop_duplicates('canonical_UniProtAC').set_index('canonical_UniProtAC') \
        if not local_density.empty else local_density
    best_region = read_optional(clustering_root / 'arg_methyl_multiscale_density_best_region_per_protein.tsv')
    best_region = best_region.drop_duplicates('canonical_UniProtAC').set_index('canonical_UniProtAC') \
        if not best_region.empty else best_region
    sig_regions = defaultdict(list)
    for _, row in read_optional(clustering_root / 'arg_methyl_multiscale_density_significant_regions.tsv').iterrows():
        sig_regions[row['canonical_UniProtAC']].append([
            to_int(row['window_aa']), to_int(row['best_window_start']), to_int(row['best_window_end']),
            to_int(row['observed_max_methyl_sites']), to_float(row['log2_observed_over_null_mean']),
            to_sig(row['q_value']),
        ])
    architecture = read_optional(root / 'domain_context' / 'domain_cluster_architecture_classes_ge3.tsv')
    architecture = architecture.drop_duplicates('canonical_UniProtAC').set_index('canonical_UniProtAC') \
        if not architecture.empty else architecture

    condensate_membership = read_optional(root / 'condensates' / 'condensate_membership.tsv')
    condensate_enrichment = read_optional(root / 'condensates' / 'condensate_enrichment.tsv')
    # CD-CODE condensates: id, type, members in the CD-CODE universe and among methylarginine proteins,
    # plus condensate-level enrichment where the condensate passed the module's size filter.
    union_accessions = set(merged['canonical_UniProtAC'])
    condensates = Vocab()
    condensate_info: dict[int, dict] = {}
    protein_condensates = defaultdict(set)
    for row in condensate_membership.to_dict('records'):
        idx = condensates.index(row['condensate_id'])
        info = condensate_info.setdefault(idx, {'name': row['condensate_name'], 'type': row['condensate_type'],
                                               'members': set(), 'methyl': set()})
        info['members'].add(row['canonical_UniProtAC'])
        if row['canonical_UniProtAC'] in union_accessions:
            info['methyl'].add(row['canonical_UniProtAC'])
            protein_condensates[row['canonical_UniProtAC']].add(idx)
    condensate_stats = {}
    for row in condensate_enrichment.to_dict('records'):
        if row['condensate_id'] in condensates.lookup:
            condensate_stats[condensates.lookup[row['condensate_id']]] = row
    condensate_rows = []
    for idx, cid in enumerate(condensates.values):
        info = condensate_info[idx]
        stat = condensate_stats.get(idx)
        condensate_rows.append([
            info['name'], cid, info['type'],
            to_int(stat['member_count']) if stat else len(info['members']),
            to_int(stat['methyl_protein_count']) if stat else len(info['methyl']),
            to_float(stat['log2_odds_ratio']) if stat else None,
            to_sig(stat['q_value']) if stat else None,
        ])

    # Disease: enriched DisGeNET terms (clusterProfiler) mapped to the audited primary disease families.
    disease_root = root / 'functional_class_union' / 'disease_enrichment'
    audit = read_optional(disease_root / 'all_methyl_proteins_disease_term_family_audit_audited_v1.tsv')
    term_family = dict(zip(audit.get('term_id', []), audit.get('primary_family', [])))
    overview = read_optional(disease_root / 'all_methyl_proteins_disease_overview_audited_v3_equal_size.tsv')
    disease_families = Vocab(list(overview.get('disease_family', [])))
    disease_terms = []
    gene_terms = defaultdict(list)
    enriched = read_optional(disease_root / 'all_methyl_proteins_disease_enrichment.tsv')
    for _, row in enriched.iterrows():
        family = term_family.get(row['ID'], 'Other / unassigned')
        family_idx = disease_families.lookup.get(family, -1)
        term_idx = len(disease_terms)
        disease_terms.append([row['ID'], row['Description'], family_idx, to_int(row['Count']), to_sig(row['p.adjust'])])
        for gene in split_multi(row['geneID'], sep='/'):
            gene_terms[gene].append(term_idx)
    als_genes = set()
    als = read_optional(disease_root / 'all_methyl_proteins_als_ftd_contributing_genes_audited_v1.tsv')
    for _, row in als.iterrows():
        als_genes.add(row['gene_symbol'])

    # Other PTM groups on the same canonical proteins: the cross-PTM master table (source bit 1) and, for lysine
    # methylation, the sequence-validated CPLM 4.0 sites from 39_stage_cplm4_lys_methylation.py (bits 2 and 4).
    other_ptms: dict[str, dict[tuple[int, int, str], int]] = defaultdict(dict)
    ptm_path = root / 'ptm_compare' / 'ptm_sites_with_contexts.tsv'
    if ptm_path.exists():
        ptm = read_tsv(ptm_path, usecols=['ptm_group', 'substrate_UniProtAC', 'canonical_UniProtAC', 'residue', 'position'])
        ptm = ptm[ptm['ptm_group'].isin(OTHER_PTM_GROUPS) & (ptm['substrate_UniProtAC'] == ptm['canonical_UniProtAC'])]
        for row in ptm.itertuples(index=False):
            pos = to_int(row.position)
            if pos is not None:
                key = (pos, OTHER_PTM_GROUPS.index(row.ptm_group), row.residue[:1])
                other_ptms[row.canonical_UniProtAC][key] = other_ptms[row.canonical_UniProtAC].get(key, 0) | 1
    cplm_path = Path(args.cplm4_sites) if args.cplm4_sites else root / 'cplm4' / 'cplm4_human_lys_methyl_validated.tsv'
    cplm_meta: dict = {}
    cplm_bits = sum(CPLM4_SOURCE_BITS.values())
    if cplm_path.exists():
        cplm = read_tsv(cplm_path, usecols=['canonical_UniProtAC', 'position', 'evidence'])
        summary_path = cplm_path.parent / 'cplm4_summary.json'
        if summary_path.exists():
            staged = json.loads(summary_path.read_text(encoding='utf-8'))
            cplm_meta = {k: staged.get(k) for k in ('input_sha256', 'input_member')}
        cplm_meta.update({'validated_sites': len(cplm), 'validated_proteins': int(cplm['canonical_UniProtAC'].nunique())})
        for row in cplm.itertuples(index=False):
            pos = to_int(row.position)
            if pos is not None:
                key = (pos, OTHER_PTM_GROUPS.index('Lys methylation'), 'K')
                other_ptms[row.canonical_UniProtAC][key] = other_ptms[row.canonical_UniProtAC].get(key, 0) | cplm4_evidence_mask(row.evidence)
    else:
        print(f'no CPLM 4.0 table at {cplm_path}; Lys methylation from the master table only')

    func_classes = Vocab()
    func_subclasses = Vocab()
    proteins: list[list] = []
    protein_index: dict[str, int] = {}
    protein_disease_terms: list[list[int]] = []
    for accession, group in merged.groupby('canonical_UniProtAC', sort=True):
        first = group.iloc[0]
        func = functional.loc[accession] if accession in functional.index else None
        odds = arg_odds.loc[accession] if accession in arg_odds.index else None
        up = uniprot.get(accession, {})
        # Current UniProt symbol first; source tables carry outdated or paralog symbols on some rows.
        row_genes = Counter(g for g in group['substrate_genename'] if g)
        func_genes = gene_symbols(func['genename']) if func is not None else []
        gene = up.get('gene') or (row_genes.most_common(1)[0][0] if row_genes else '') \
            or (func_genes[0] if func_genes else '')
        aliases = set(row_genes) | set(up.get('synonyms', []))
        if func is not None:
            aliases.update(func_genes)
        aliases = {a for a in aliases if a != gene and ' ' not in a and ':' not in a}
        name = recommended_name(up.get('name') or first['protein_name']
                                or (func['protein_name'] if func is not None else ''))
        uniprot_id = up.get('entry_name') or first['UniProtID'] or (func['UniProtID'] if func is not None else '')
        length = up.get('length')
        if length is not None and length < int(group['position_int'].max()):
            length = None
        if length is None and odds is not None:
            length = to_int(odds['protein_length'])

        density = None
        if not local_density.empty and accession in local_density.index:
            ld = local_density.loc[accession]
            density = [to_float(ld['nearest_distance_min'], 1), to_float(ld['nearest_distance_median'], 1)] + [
                to_int(ld[f'max_methyl_sites_in_{w}aa_window']) for w in (3, 5, 10, 20, 50)
            ]
        best = None
        if not best_region.empty and accession in best_region.index:
            br = best_region.loc[accession]
            best = [to_int(br['window_aa']), to_int(br['best_window_start']), to_int(br['best_window_end']),
                    to_int(br['observed_max_methyl_sites']), to_float(br['log2_observed_over_null_mean']),
                    to_sig(br['q_value'])]
        arch = None
        if not architecture.empty and accession in architecture.index:
            ar = architecture.loc[accession]
            clusters = []
            for chunk in split_multi(ar['best_cluster_positions']):
                positions = [to_int(p) for p in chunk.split(',') if to_int(p) is not None]
                if positions:
                    clusters.append([min(positions), max(positions), len(positions)])
            arch = [ar['architecture_class'], clusters]

        symbols = {gene} | aliases
        terms = sorted({t for symbol in symbols for t in gene_terms.get(symbol, [])})
        family_mask = 0
        for term_idx in terms:
            family_idx = disease_terms[term_idx][2]
            if family_idx >= 0:
                family_mask |= 1 << family_idx

        protein_index[accession] = len(proteins)
        protein_disease_terms.append(terms)
        proteins.append([
            accession,
            gene,
            name,
            uniprot_id,
            func_classes.mask(split_multi(func['functional_class_broad'])) if func is not None else 0,
            func_subclasses.mask(split_multi(func['functional_class_subclass'])) if func is not None else 0,
            length,
            to_int(odds['arginine_count']) if odds is not None else None,
            to_float(odds['shrinkage_log2_enrichment']) if odds is not None else None,
            to_sig(odds['q_value']) if odds is not None else None,
            int(first['accession_is_isoform'] == 'True'),
            sorted(protein_condensates.get(accession, set())),
            family_mask,
            int(bool(symbols & als_genes)),
            density,
            best,
            sig_regions.get(accession, []),
            arch,
            [[*key, mask] for key, mask in sorted(other_ptms.get(accession, {}).items())],
            sorted(aliases),
        ])

    # ---- Sites --------------------------------------------------------------------------------
    source_vocab = Vocab(SOURCES)
    refs = Vocab()
    interpro: list[list[str]] = []
    interpro_lookup: dict[str, int] = {}
    disorder_classes = Vocab()
    domain_contexts = Vocab()
    domain_classes = Vocab()
    window_flags = {'missing': 0, 'center_not_r': 0}
    site_rows: list[list] = []
    for i, row in enumerate(merged.to_dict('records')):
        state = row['arg_methyl_state_consensus']
        window = row['window_31']
        center_mismatch = 0
        if len(window) != 2 * WINDOW_FLANK + 1 or set(window) == {'_'}:
            window = ''
            window_flags['missing'] += 1
        elif window[WINDOW_FLANK] != 'R':
            center_mismatch = 1
            window_flags['center_not_r'] += 1
        motif_mask = 0
        if window:
            for bit, motif in enumerate(MOTIFS):
                if row[f'motif_{motif}'] == 'True':
                    motif_mask |= 1 << bit
        ipr_idx = -1
        if row['nearest_interpro_accession']:
            acc = row['nearest_interpro_accession']
            if acc not in interpro_lookup:
                interpro_lookup[acc] = len(interpro)
                interpro.append([acc, row['nearest_interpro_name']])
            ipr_idx = interpro_lookup[acc]
        original_positions = row['source_positions_original']
        site_rows.append([
            protein_index[row['canonical_UniProtAC']],
            int(row['position_int']),
            STATES.index(state) if state in STATES else 0,
            source_vocab.mask(split_multi(row['source_families_exact'])),
            [refs.index(clean_reference(ref)) for ref in split_multi(row['reference_any'])],
            disorder_classes.index(row['disorder_context_class']),
            to_int(row['nearest_disorder_edge_distance']),
            domain_contexts.index(row['domain_context_class']),
            ipr_idx,
            domain_classes.index(row['nearest_domain_class']),
            to_int(row['nearest_domain_edge_distance']),
            row['nearest_domain_side'],
            window,
            motif_mask,
            center_mismatch,
            '' if original_positions == str(row['position_int']) else original_positions,
            s2_state.get(i, 0),
            ip_mask.get(i, 0),
            treatment.get(i, []),
        ])

    # ---- Evidence (lazy-loaded): per-site source rows and per-protein disease terms --------------
    source_tables = Vocab()
    remap_statuses = Vocab()
    studies = Vocab()
    study_citations: dict[str, str] = {}
    study_masks = [0] * len(site_rows)
    evidence: list[list] = [[] for _ in site_rows]
    rows = read_optional(root / 'remapped' / 'human_arg_methyl_source_rows_remapped_resolved.tsv')
    evidence_unmatched = 0
    for row in rows.to_dict('records'):
        idx = site_lookup(row['corrected_canonical_UniProtAC'] or row['canonical_UniProtAC'], row['corrected_position'])
        if idx is None:
            evidence_unmatched += 1
            continue
        reference = clean_reference(row['reference'].strip())
        if row['source_table'] == 'maron_s5' and reference:
            study = study_label(reference)
            study_citations[study] = reference
        elif row['source_table'] == 'prometheus_mmc2':
            study = 'ProMetheusDB'
        else:
            study = PTM_MASTER_SOURCES.get(row['source'], row['source'] or row['source_table'])
        study_masks[idx] |= 1 << studies.index(study)
        evidence[idx].append([
            source_tables.index(row['source_table']),
            refs.index(reference) if reference else -1,
            row['pmid'],
            row['peptide_sequence'] or row['peptides'],
            row['enzyme_genename'],
            row['arg_methyl_state'],
            row['substrate_UniProtAC_original'],
            row['position'],
            remap_statuses.index(row['remap_status']),
            studies.index(study),
        ])

    for idx, mask in enumerate(study_masks):
        site_rows[idx].append(mask)
    study_counts = [0] * len(studies.values)
    study_unique = [0] * len(studies.values)
    for mask in study_masks:
        for bit in range(len(studies.values)):
            if mask & (1 << bit):
                study_counts[bit] += 1
                if mask == 1 << bit:
                    study_unique[bit] += 1

    # ---- Overview tables -------------------------------------------------------------------------
    def table(path: Path, label_col: str = 'category') -> list[list]:
        df = read_optional(path)
        return enrichment_rows(df, label_col) if not df.empty else []

    overview_tables = {
        'functional_classes': table(root / 'functional_class_union' / 'arg_methyl_functional_class_broad_enrichment.tsv'),
        'motifs': table(root / 'motif_arg_odds' / 'motif_family_enrichment_vs_methylated_protein_arginines.tsv', 'motif_family'),
        'disorder': table(root / 'disorder_context' / 'idr_enrichment_vs_arginines_in_methylated_proteins.tsv'),
        'domain_context': table(root / 'domain_context' / 'domain_context_enrichment.tsv'),
        'domain_classes': table(root / 'domain_context' / 'domain_class_enrichment.tsv'),
        'disease_families': [
            [row['disease_family'], to_int(row['methylarginine_protein_count']), to_int(row['family_gene_count']),
             to_float(row['odds_ratio']), to_sig(row['BH_adjusted_p'])]
            for _, row in overview.iterrows()
        ],
        'treatments': [
            [row['treatment'], row['ip'], to_int(row['unique_canonical_sites']), to_int(row['sig_down_sites']),
             to_int(row['sig_up_sites']), to_float(row['median_diff'])]
            for _, row in read_optional(maron_root / 'maron_s3_treatment_ip_summary.tsv').iterrows()
        ],
    }
    parent_counts = defaultdict(Counter)
    for row in functional.reset_index().to_dict('records'):
        for sub in split_multi(row['functional_class_subclass']):
            for broad in split_multi(row['functional_class_broad']):
                parent_counts[sub][broad] += 1
    overview_tables['functional_subclasses'] = [
        r + [parent_counts[r[0]].most_common(1)[0][0] if parent_counts[r[0]] else '']
        for r in table(root / 'functional_class_union' / 'arg_methyl_functional_class_subclass_enrichment.tsv')
    ]
    condensate_summary = read_optional(root / 'condensates' / 'condensate_summary.tsv')

    payload = {
        'meta': {
            'title': 'Human methylarginine site union',
            'generated': dt.date.today().isoformat(),
            'n_sites': len(site_rows),
            'n_proteins': len(proteins),
            'window_flank': WINDOW_FLANK,
            'window_flags': window_flags,
            'sources': SOURCES,
            'studies': [
                [name, study_citations.get(name, ''), study_counts[k], study_unique[k]]
                for k, name in enumerate(studies.values)
            ],
            'states': STATES,
            'motifs': MOTIFS,
            'maron_ips': MARON_IPS,
            'maron_treatments': MARON_TREATMENTS,
            'other_ptm_groups': OTHER_PTM_GROUPS,
            'other_ptm_sources': OTHER_PTM_SOURCES,
            'cplm4': {**cplm_meta,
                      'explorer_sites': sum(1 for p in proteins for item in p[18] if item[3] & cplm_bits),
                      'explorer_proteins': sum(1 for p in proteins if any(item[3] & cplm_bits for item in p[18]))} if cplm_meta else {},
            'func_classes': func_classes.values,
            'func_subclasses': func_subclasses.values,
            'disorder_classes': disorder_classes.values,
            'domain_contexts': domain_contexts.values,
            'domain_classes': domain_classes.values,
            'disease_families': disease_families.values,
            'condensates': condensate_rows,
            'condensate_fields': ['name', 'cdcode_id', 'type', 'members', 'methyl_members', 'log2_odds_ratio', 'q_value'],
            'condensate_summary': condensate_summary.iloc[0].to_dict() if not condensate_summary.empty else {},
            'refs': refs.values,
            'interpro': interpro,
            'overview': overview_tables,
            'protein_fields': [
                'accession', 'gene', 'protein_name', 'uniprot_id', 'func_class_mask', 'func_subclass_mask', 'length',
                'arginine_count', 'shrinkage_log2_enrichment', 'arg_odds_q', 'is_isoform', 'condensates',
                'disease_family_mask', 'als_ftd_audited_gene', 'local_density', 'best_density_window',
                'significant_density_windows', 'cluster_architecture', 'other_ptms', 'gene_aliases',
            ],
            'site_fields': [
                'protein', 'position', 'state', 'source_mask', 'refs', 'disorder_class', 'disorder_edge_distance',
                'domain_context', 'interpro', 'domain_class', 'domain_edge_distance', 'domain_side', 'window_31',
                'motif_mask', 'window_center_mismatch', 'source_positions_original', 'maron_s2_state_mask',
                'maron_ip_mask', 'maron_treatment', 'study_mask',
            ],
        },
        'proteins': proteins,
        'sites': site_rows,
    }
    evidence_payload = {
        'source_tables': source_tables.values,
        'studies': studies.values,
        'remap_statuses': remap_statuses.values,
        'evidence_fields': ['source_table', 'reference', 'pmid', 'peptide', 'enzyme', 'state',
                            'source_accession', 'source_position', 'remap_status', 'study'],
        'refs': refs.values,
        'sites': evidence,
        'disease_terms': disease_terms,
        'protein_disease_terms': protein_disease_terms,
    }
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for name, data in [('methylarg_sites.json', payload), ('methylarg_evidence.json', evidence_payload)]:
        path = outdir / name
        path.write_text(json.dumps(data, separators=(',', ':'), ensure_ascii=False), encoding='utf-8')
        print(f'wrote {path} ({path.stat().st_size / 1e6:.2f} MB)')
    print(f'{len(site_rows)} sites, {len(proteins)} proteins ({sum(p[6] is not None for p in proteins)} with lengths, '
          f'{sum(1 for p in proteins if p[19])} with gene aliases); '
          f'window flags {window_flags}; Maron S3 rows unmatched {maron_unmatched}; '
          f'evidence rows unmatched {evidence_unmatched}')


if __name__ == '__main__':
    main()

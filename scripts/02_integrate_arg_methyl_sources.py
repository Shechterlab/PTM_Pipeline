from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import pandas as pd

from common import (
    accession_is_isoform,
    canonicalize_uniprot_accession,
    confidence_tier,
    consensus_or_mixed,
    normalize_accession,
    normalize_arg_methyl_state,
    save_table,
    source_family_from_base_source,
    summarize_sources,
)


def first_nonempty(values: Iterable[object]) -> str:
    for value in values:
        text = str(value).strip()
        if text not in {'', 'nan', 'None'}:
            return text
    return ''


def any_true(values: Iterable[object]) -> bool:
    for value in values:
        if isinstance(value, bool):
            if value:
                return True
            continue
        text = str(value).strip().lower()
        if text in {'true', '1', 'yes'}:
            return True
    return False


def canonical_counting_status(is_isoform: bool) -> str:
    return 'isoform_needs_sequence_remap' if is_isoform else 'canonical_or_unversioned_accession'


def prepare_source_frame(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    if 'substrate_UniProtAC_original' not in out.columns:
        out['substrate_UniProtAC_original'] = out['substrate_UniProtAC']
    out['substrate_UniProtAC_original'] = out['substrate_UniProtAC_original'].astype(str).str.strip()
    out['substrate_UniProtAC'] = out['substrate_UniProtAC'].astype(str).map(normalize_accession)

    if 'canonical_UniProtAC' not in out.columns:
        out['canonical_UniProtAC'] = out['substrate_UniProtAC'].map(canonicalize_uniprot_accession)
    else:
        out['canonical_UniProtAC'] = out['canonical_UniProtAC'].astype(str).str.strip()
        out['canonical_UniProtAC'] = out['canonical_UniProtAC'].where(
            out['canonical_UniProtAC'].ne(''),
            out['substrate_UniProtAC'].map(canonicalize_uniprot_accession),
        )

    if 'accession_is_isoform' not in out.columns:
        out['accession_is_isoform'] = out['substrate_UniProtAC'].map(accession_is_isoform)
    else:
        out['accession_is_isoform'] = out['accession_is_isoform'].map(
            lambda x: str(x).strip().lower() in {'true', '1', 'yes'}
        )

    if 'mod_state' not in out.columns:
        out['mod_state'] = ''
    out['mod_state'] = out['mod_state'].astype(str).str.strip()

    if 'arg_methyl_state' not in out.columns:
        out['arg_methyl_state'] = out['mod_state'].map(normalize_arg_methyl_state)
    else:
        out['arg_methyl_state'] = out['arg_methyl_state'].map(normalize_arg_methyl_state)
    out['arg_methyl_state_known'] = out['arg_methyl_state'].astype(str).ne('')

    if 'data_source_label' not in out.columns:
        out['data_source_label'] = ''
    if 'reference' not in out.columns:
        out['reference'] = ''
    if 'protein_name' not in out.columns:
        out['protein_name'] = ''
    if 'UniProtID' not in out.columns:
        out['UniProtID'] = ''
    if 'peptide_sequence' not in out.columns:
        out['peptide_sequence'] = ''
    if 'seq_window' not in out.columns:
        out['seq_window'] = ''
    if 'peptides' not in out.columns:
        out['peptides'] = ''

    out['position'] = pd.to_numeric(out['position'], errors='coerce')
    out = out[out['position'].notna()].copy()
    out['position'] = out['position'].astype(int)
    out['residue'] = 'R'
    out['site'] = 'R' + out['position'].astype(str)
    out['site_key'] = out['substrate_UniProtAC'] + ':' + out['site']
    out['site_key_canonical_naive'] = out['canonical_UniProtAC'] + ':' + out['site']
    out['needs_sequence_remap'] = out['accession_is_isoform']
    out['canonical_counting_status'] = out['accession_is_isoform'].map(canonical_counting_status)
    out['is_arg_methyl'] = True
    out['ptm_group'] = 'Arg methylation'
    out['ptm_type'] = 'METHYLATION'
    out['residue'] = 'R'
    out['data_source_label'] = out['data_source_label'].where(
        out['data_source_label'].astype(str).ne(''),
        out['source_family'].astype(str),
    )
    return out


def fuzzy_support_for_sites(all_sites: pd.DataFrame, exact: pd.DataFrame, tolerance: int) -> pd.DataFrame:
    pos_by_acc_family: dict[tuple[str, str], list[int]] = {}
    for (acc, fam), group in all_sites.groupby(['substrate_UniProtAC', 'source_family']):
        pos_by_acc_family[(acc, fam)] = sorted(set(group['position'].astype(int)))

    out_rows = []
    fams_by_acc = all_sites.groupby('substrate_UniProtAC')['source_family'].unique().to_dict()
    for row in exact.itertuples(index=False):
        acc = row.substrate_UniProtAC
        pos = int(row.position)
        matched = []
        for fam in fams_by_acc.get(acc, []):
            fam_positions = pos_by_acc_family.get((acc, fam), [])
            hit = False
            for candidate_pos in fam_positions:
                if candidate_pos < pos - tolerance:
                    continue
                if candidate_pos > pos + tolerance:
                    break
                hit = True
                break
            if hit:
                matched.append(str(fam))
        out_rows.append({
            'substrate_UniProtAC': acc,
            'site': row.site,
            'source_family_count_fuzzy': len(sorted(set(matched))),
            'source_families_fuzzy': ';'.join(sorted(set(matched))),
        })
    return pd.DataFrame(out_rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--base-master')
    ap.add_argument('--maron')
    ap.add_argument('--prometheus')
    ap.add_argument('--premerged-remapped', default='')
    ap.add_argument('--outdir', default='results/integrated')
    ap.add_argument('--fuzzy-tolerance', type=int, default=2)
    args = ap.parse_args()

    if args.premerged_remapped:
        all_sites = pd.read_csv(args.premerged_remapped, sep='\t', low_memory=False)
        all_sites = all_sites[all_sites.get('is_position_resolved', False).fillna(False)].copy()
        all_sites = prepare_source_frame(all_sites)
        all_sites['source_position_original'] = pd.to_numeric(all_sites['position'], errors='coerce')
        all_sites['canonical_UniProtAC'] = all_sites.get('corrected_canonical_UniProtAC', all_sites['canonical_UniProtAC']).astype(str)
        all_sites['corrected_position'] = pd.to_numeric(all_sites.get('corrected_position', all_sites['position']), errors='coerce')
        all_sites.dropna(subset=['corrected_position'], inplace=True)
        all_sites['corrected_position'] = all_sites['corrected_position'].astype(int)
        all_sites['position'] = all_sites['corrected_position']
        all_sites['site'] = 'R' + all_sites['position'].astype(str)
        all_sites['site_key'] = all_sites['canonical_UniProtAC'] + ':' + all_sites['site']
        all_sites['integration_accession'] = all_sites['canonical_UniProtAC']
        all_sites['source_accessions_original'] = all_sites['substrate_UniProtAC'].astype(str)
    else:
        if not all([args.base_master, args.maron, args.prometheus]):
            raise ValueError('either --premerged-remapped or all of --base-master/--maron/--prometheus are required')
        base = pd.read_csv(args.base_master, sep='\t', low_memory=False)
        base = base[base['is_arg_methyl'].fillna(False)].copy()
        base['source_family'] = base['source'].map(source_family_from_base_source)
        base['data_source_label'] = 'human_ptm_master:' + base['source'].astype(str)

        maron = pd.read_csv(args.maron, sep='\t', low_memory=False)
        prom = pd.read_csv(args.prometheus, sep='\t', low_memory=False)

        frames = [prepare_source_frame(df) for df in [base, maron, prom]]
        keep_cols = sorted(set().union(*(frame.columns for frame in frames)))
        all_sites = pd.concat([frame.reindex(columns=keep_cols) for frame in frames], ignore_index=True, sort=False)
        all_sites['integration_accession'] = all_sites['substrate_UniProtAC']
        all_sites['source_accessions_original'] = all_sites['substrate_UniProtAC'].astype(str)

    all_sites = all_sites[all_sites['substrate_UniProtAC'].astype(str).ne('')].copy()
    all_sites['source_family'] = all_sites['source_family'].astype(str)
    if 'integration_accession' not in all_sites.columns:
        all_sites['integration_accession'] = all_sites['substrate_UniProtAC']

    exact = all_sites.groupby(['integration_accession', 'site'], as_index=False).agg(
        substrate_UniProtAC_original=('substrate_UniProtAC_original', summarize_sources),
        canonical_UniProtAC=('canonical_UniProtAC', first_nonempty),
        accession_is_isoform=('accession_is_isoform', any_true),
        residue=('residue', 'first'),
        position=('position', 'median'),
        corrected_position=('position', 'median'),
        source_positions_original=('source_position_original', summarize_sources),
        substrate_genename=('substrate_genename', first_nonempty),
        source_accessions_original=('source_accessions_original', summarize_sources),
        source_family_count_exact=('source_family', lambda s: len(set(s.astype(str)))),
        source_families_exact=('source_family', summarize_sources),
        source_labels_exact=('data_source_label', summarize_sources),
        reference_any=('reference', summarize_sources),
        protein_name=('protein_name', first_nonempty),
        UniProtID=('UniProtID', first_nonempty),
        arg_methyl_states_observed=('arg_methyl_state', summarize_sources),
        arg_methyl_state_consensus=('arg_methyl_state', consensus_or_mixed),
        n_rows_union=('site_key', 'size'),
        n_state_annotated_rows=('arg_methyl_state', lambda s: sum(str(v) not in {'', 'nan', 'None'} for v in s)),
    )
    exact = exact.rename(columns={'integration_accession': 'substrate_UniProtAC'})
    exact['position'] = exact['position'].round().astype(int)
    exact['corrected_position'] = exact['corrected_position'].round().astype(int)
    exact['site_key'] = exact['substrate_UniProtAC'] + ':' + exact['site']
    exact['site_key_canonical_naive'] = exact['canonical_UniProtAC'] + ':' + exact['site']
    exact['needs_sequence_remap'] = exact['accession_is_isoform']
    exact['canonical_counting_status'] = exact['accession_is_isoform'].map(canonical_counting_status)
    exact['arg_methyl_state_known'] = exact['arg_methyl_states_observed'].astype(str).ne('')

    state_support = all_sites[all_sites['arg_methyl_state'].astype(str).ne('')].groupby(
        ['integration_accession', 'site'],
        as_index=False,
    ).agg(
        state_annotated_source_family_count=('source_family', lambda s: len(set(s.astype(str)))),
        state_annotated_source_families=('source_family', summarize_sources),
    )
    state_support = state_support.rename(columns={'integration_accession': 'substrate_UniProtAC'})
    exact = exact.merge(state_support, how='left', on=['substrate_UniProtAC', 'site'])
    exact['state_annotated_source_family_count'] = exact['state_annotated_source_family_count'].fillna(0).astype(int)
    exact['state_annotated_source_families'] = exact['state_annotated_source_families'].fillna('')

    fuzzy = fuzzy_support_for_sites(
        all_sites[['integration_accession', 'position', 'source_family']].rename(columns={'integration_accession': 'substrate_UniProtAC'}).copy(),
        exact[['substrate_UniProtAC', 'site', 'position']].copy(),
        tolerance=args.fuzzy_tolerance,
    )
    exact = exact.merge(fuzzy, how='left', on=['substrate_UniProtAC', 'site'])
    exact['source_family_count_fuzzy'] = exact['source_family_count_fuzzy'].fillna(
        exact['source_family_count_exact']
    ).astype(int)
    exact['source_families_fuzzy'] = exact['source_families_fuzzy'].fillna(exact['source_families_exact'])
    exact['confidence_tier'] = exact.apply(confidence_tier, axis=1)
    exact['support_any_ge_2'] = exact['source_family_count_fuzzy'] >= 2
    exact['support_exact_ge_2'] = exact['source_family_count_exact'] >= 2

    accession_audit = all_sites.groupby('integration_accession', as_index=False).agg(
        substrate_UniProtAC_original=('substrate_UniProtAC_original', summarize_sources),
        canonical_UniProtAC=('canonical_UniProtAC', first_nonempty),
        accession_is_isoform=('accession_is_isoform', any_true),
        canonical_counting_status=('canonical_counting_status', first_nonempty),
        substrate_genename=('substrate_genename', first_nonempty),
        unique_sites_accession_level=('site', 'nunique'),
        source_accessions_original=('source_accessions_original', summarize_sources),
        source_families=('source_family', summarize_sources),
    )
    accession_audit = accession_audit.rename(columns={'integration_accession': 'substrate_UniProtAC'})
    accession_audit['needs_sequence_remap'] = accession_audit['accession_is_isoform']

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(all_sites, outdir / 'human_arg_methyl_union_all_rows.tsv')
    save_table(exact, outdir / 'human_arg_methyl_union_dedup_by_site.tsv')
    save_table(exact[exact['support_any_ge_2']].copy(), outdir / 'human_arg_methyl_support_ge_2_fuzzy.tsv')
    save_table(exact[exact['support_exact_ge_2']].copy(), outdir / 'human_arg_methyl_support_ge_2_exact.tsv')
    save_table(accession_audit, outdir / 'accession_canonicalization_audit.tsv')

    support = exact['confidence_tier'].value_counts().rename_axis('confidence_tier').reset_index(name='site_count')
    save_table(support, outdir / 'confidence_tier_counts.tsv')
    overlap = exact[['source_family_count_exact', 'source_family_count_fuzzy']].value_counts().rename_axis(
        ['source_family_count_exact', 'source_family_count_fuzzy']
    ).reset_index(name='site_count')
    save_table(overlap, outdir / 'support_exact_vs_fuzzy.tsv')

    summary = pd.DataFrame([{
        'rows_union_all': len(all_sites),
        'sites_union_dedup': len(exact),
        'proteins_union_dedup': exact['substrate_UniProtAC'].nunique(),
        'canonical_accessions_union_naive': exact['canonical_UniProtAC'].nunique(),
        'sites_support_ge_2_exact': int(exact['support_exact_ge_2'].sum()),
        'sites_support_ge_2_fuzzy': int(exact['support_any_ge_2'].sum()),
        'sites_with_state_annotation': int(exact['arg_methyl_state_known'].sum()),
        'sites_needing_sequence_remap': int(exact['needs_sequence_remap'].sum()),
    }])
    save_table(summary, outdir / 'summary.tsv')
    print(summary.to_dict(orient='records')[0])


if __name__ == '__main__':
    main()

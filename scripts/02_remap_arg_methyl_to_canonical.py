from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd

from common import (
    accession_is_isoform,
    canonical_gene_name,
    canonicalize_uniprot_accession,
    centered_anchor_indices,
    clean_aa_sequence,
    normalize_accession,
    parse_fasta,
    parse_modified_peptide,
    save_table,
    split_multi_value_field,
)


def load_uniprot_metadata(path: str | Path) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    meta = pd.read_csv(path, sep='\t', low_memory=False)
    meta = meta.rename(columns={
        'Entry': 'canonical_UniProtAC',
        'Gene Names (primary)': 'gene_primary',
    })
    meta['canonical_UniProtAC'] = meta['canonical_UniProtAC'].astype(str).str.strip()
    meta['gene_primary'] = meta.get('gene_primary', '').astype(str).map(canonical_gene_name)
    gene_to_accessions: dict[str, list[str]] = defaultdict(list)
    for row in meta.itertuples(index=False):
        gene = canonical_gene_name(getattr(row, 'gene_primary', ''))
        acc = str(getattr(row, 'canonical_UniProtAC', '')).strip()
        if gene and acc:
            gene_to_accessions[gene].append(acc)
    gene_to_accessions = {gene: sorted(set(accs)) for gene, accs in gene_to_accessions.items()}
    return meta, gene_to_accessions


def prepare_frame(df: pd.DataFrame, source_name: str) -> pd.DataFrame:
    out = df.copy()
    if 'substrate_UniProtAC_original' not in out.columns:
        out['substrate_UniProtAC_original'] = out['substrate_UniProtAC']
    out['substrate_UniProtAC_original'] = out['substrate_UniProtAC_original'].astype(str).str.strip()
    out['substrate_UniProtAC'] = out['substrate_UniProtAC'].astype(str).map(normalize_accession)
    out['canonical_UniProtAC'] = out.get('canonical_UniProtAC', out['substrate_UniProtAC']).astype(str).map(canonicalize_uniprot_accession)
    out['accession_is_isoform'] = out.get('accession_is_isoform', out['substrate_UniProtAC'].map(accession_is_isoform))
    out['substrate_genename'] = out.get('substrate_genename', '').astype(str).map(canonical_gene_name)
    out['position'] = pd.to_numeric(out['position'], errors='coerce')
    out = out[out['position'].notna()].copy()
    out['position'] = out['position'].astype(int)
    out['source_table'] = source_name
    if 'source_family' not in out.columns:
        out['source_family'] = source_name
    return out


def find_all_occurrences(sequence: str, anchor: str) -> list[int]:
    if not sequence or not anchor:
        return []
    starts = []
    start = 0
    while True:
        idx = sequence.find(anchor, start)
        if idx == -1:
            break
        starts.append(idx + 1)
        start = idx + 1
    return starts


def candidate_accessions(row: pd.Series, canonical_sequences: dict[str, str], gene_map: dict[str, list[str]]) -> list[str]:
    candidates: list[str] = []
    for field in ['canonical_UniProtAC', 'substrate_UniProtAC', 'UniProtAC', 'substrate_UniProtAC_original']:
        if field not in row.index:
            continue
        acc = canonicalize_uniprot_accession(row.get(field, ''))
        if acc and acc in canonical_sequences and acc not in candidates:
            candidates.append(acc)
    gene = canonical_gene_name(row.get('substrate_genename', ''))
    for acc in gene_map.get(gene, []):
        if acc in canonical_sequences and acc not in candidates:
            candidates.append(acc)
    return candidates


def build_anchors(row: pd.Series) -> list[dict]:
    anchors: list[dict] = []
    seen: set[tuple[str, tuple[int, ...], str]] = set()

    def add_anchor(sequence: str, indices: list[int], anchor_type: str) -> None:
        seq = clean_aa_sequence(sequence)
        clean_indices = sorted(set(idx for idx in indices if 1 <= idx <= len(seq)))
        if not seq or not clean_indices:
            return
        key = (seq, tuple(clean_indices), anchor_type)
        if key in seen:
            return
        seen.add(key)
        anchors.append({
            'sequence': seq,
            'indices': clean_indices,
            'anchor_type': anchor_type,
        })

    for peptide in split_multi_value_field(row.get('peptides', ''), separators=','):
        clean, modified = parse_modified_peptide(peptide)
        add_anchor(clean, modified, 'annotated_peptide')

    for field, anchor_type in [('seq_window', 'centered_window'), ('peptide_sequence', 'centered_window')]:
        for token in split_multi_value_field(row.get(field, ''), separators=';'):
            clean = clean_aa_sequence(token)
            add_anchor(clean, centered_anchor_indices(clean), anchor_type)

    return anchors


def fallback_position(row: pd.Series, canonical_sequences: dict[str, str]) -> tuple[int | None, str, str]:
    current_acc = normalize_accession(row.get('substrate_UniProtAC', ''))
    canonical_acc = canonicalize_uniprot_accession(row.get('canonical_UniProtAC', current_acc))
    pos = int(row.get('position', 0))
    if current_acc in canonical_sequences:
        seq = canonical_sequences[current_acc]
        if 1 <= pos <= len(seq) and seq[pos - 1] == 'R':
            return pos, 'same_accession_position_verified', current_acc
    if canonical_acc in canonical_sequences and not bool(row.get('accession_is_isoform', False)):
        seq = canonical_sequences[canonical_acc]
        if 1 <= pos <= len(seq) and seq[pos - 1] == 'R':
            return pos, 'canonical_position_fallback', canonical_acc
    return None, 'unresolved_no_anchor', ''


def resolve_hit(row: pd.Series, hits: list[dict], canonical_acc: str) -> tuple[int | None, str, str, str, int, str]:
    if not hits:
        return None, 'unresolved_no_hit', '', '', 0, ''

    positions_by_acc: dict[str, set[int]] = defaultdict(set)
    for hit in hits:
        positions_by_acc[hit['matched_accession']].add(hit['matched_position'])

    if canonical_acc in positions_by_acc and len(positions_by_acc[canonical_acc]) == 1:
        position = sorted(positions_by_acc[canonical_acc])[0]
        supporting = [hit for hit in hits if hit['matched_accession'] == canonical_acc and hit['matched_position'] == position]
        best = sorted(supporting, key=lambda h: (h['priority'], len(h['anchor_sequence'])))[0]
        return position, 'sequence_unique_canonical', best['anchor_type'], best['anchor_sequence'], len(hits), canonical_acc

    unique_positions = sorted(set(hit['matched_position'] for hit in hits))
    if len(unique_positions) == 1 and len(set(hit['matched_accession'] for hit in hits)) == 1:
        best = sorted(hits, key=lambda h: (h['priority'], len(h['anchor_sequence'])))[0]
        return unique_positions[0], 'sequence_unique_candidate', best['anchor_type'], best['anchor_sequence'], len(hits), best['matched_accession']

    current_position = int(row.get('position', 0))
    closest = sorted(
        hits,
        key=lambda h: (
            abs(h['matched_position'] - current_position),
            h['matched_accession'] != canonical_acc,
            h['priority'],
            -len(h['anchor_sequence']),
        ),
    )
    if closest:
        best = closest[0]
        tied = [
            hit for hit in closest
            if abs(hit['matched_position'] - current_position) == abs(best['matched_position'] - current_position)
            and hit['matched_accession'] == best['matched_accession']
        ]
        if len(set(hit['matched_position'] for hit in tied)) == 1 and abs(best['matched_position'] - current_position) <= 25:
            return best['matched_position'], 'sequence_closest_supported', best['anchor_type'], best['anchor_sequence'], len(hits), best['matched_accession']

    return None, 'unresolved_ambiguous_hit', '', '', len(hits), ''


def remap_row(row: pd.Series, canonical_sequences: dict[str, str], gene_map: dict[str, list[str]]) -> dict:
    canonical_acc = canonicalize_uniprot_accession(row.get('canonical_UniProtAC', row.get('substrate_UniProtAC', '')))
    candidates = candidate_accessions(row, canonical_sequences, gene_map)
    anchors = build_anchors(row)
    hits: list[dict] = []
    for anchor in anchors:
        for acc in candidates:
            starts = find_all_occurrences(canonical_sequences.get(acc, ''), anchor['sequence'])
            for start in starts:
                for idx in anchor['indices']:
                    position = start + idx - 1
                    hits.append({
                        'matched_accession': acc,
                        'matched_position': position,
                        'anchor_type': anchor['anchor_type'],
                        'anchor_sequence': anchor['sequence'],
                        'priority': 0 if anchor['anchor_type'] == 'annotated_peptide' else 1,
                    })

    corrected_position, remap_status, anchor_type, anchor_sequence, hit_count, corrected_accession = resolve_hit(row, hits, canonical_acc)
    if corrected_position is None:
        corrected_position, fallback_status, fallback_accession = fallback_position(row, canonical_sequences)
        if corrected_position is not None:
            remap_status = fallback_status
            anchor_type = 'position_only'
            anchor_sequence = ''
            corrected_accession = fallback_accession
        elif remap_status == 'unresolved_no_hit':
            remap_status = fallback_status

    corrected_accession = corrected_accession if corrected_position is not None else ''
    corrected_site = f'R{corrected_position}' if corrected_position is not None else ''
    remapped = corrected_position is not None
    return {
        'candidate_canonical_accessions': ';'.join(candidates),
        'candidate_canonical_accession_count': len(candidates),
        'anchor_count': len(anchors),
        'anchor_types_available': ';'.join(sorted(set(anchor['anchor_type'] for anchor in anchors))),
        'remap_hit_count': hit_count,
        'remap_status': remap_status,
        'remap_anchor_type': anchor_type,
        'remap_anchor_sequence': anchor_sequence,
        'corrected_position': corrected_position if corrected_position is not None else pd.NA,
        'corrected_site': corrected_site,
        'corrected_canonical_UniProtAC': corrected_accession,
        'position_delta': (corrected_position - int(row.get('position', 0))) if corrected_position is not None else pd.NA,
        'is_position_remapped': remapped and corrected_position != int(row.get('position', 0)),
        'is_position_resolved': remapped,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--base-master', required=True)
    ap.add_argument('--maron', required=True)
    ap.add_argument('--prometheus', required=True)
    ap.add_argument('--dbptm', default='')
    ap.add_argument('--canonical-fasta', required=True)
    ap.add_argument('--uniprot-metadata', required=True)
    ap.add_argument('--outdir', default='results/remapped')
    args = ap.parse_args()

    canonical_sequences = parse_fasta(args.canonical_fasta)
    _, gene_map = load_uniprot_metadata(args.uniprot_metadata)

    base = pd.read_csv(args.base_master, sep='\t', low_memory=False)
    base = base[base['is_arg_methyl'].fillna(False)].copy()
    maron = pd.read_csv(args.maron, sep='\t', low_memory=False)
    prom = pd.read_csv(args.prometheus, sep='\t', low_memory=False)

    frames = [
        prepare_frame(base, 'human_ptm_master'),
        prepare_frame(maron, 'maron_s5'),
        prepare_frame(prom, 'prometheus_mmc2'),
    ]
    if args.dbptm:
        dbptm = pd.read_csv(args.dbptm, sep='	', low_memory=False)
        frames.append(prepare_frame(dbptm, 'dbptm'))
    all_rows = pd.concat(frames, ignore_index=True, sort=False)
    all_rows['canonical_UniProtAC'] = all_rows['canonical_UniProtAC'].astype(str).map(canonicalize_uniprot_accession)

    remap_annotations = all_rows.apply(
        lambda row: pd.Series(remap_row(row, canonical_sequences=canonical_sequences, gene_map=gene_map)),
        axis=1,
    )
    remapped = pd.concat([all_rows, remap_annotations], axis=1)

    resolved = remapped[remapped['is_position_resolved'].fillna(False)].copy()
    resolved['canonical_UniProtAC'] = resolved['corrected_canonical_UniProtAC']
    resolved['corrected_position'] = pd.to_numeric(resolved['corrected_position'], errors='coerce').astype('Int64')
    resolved['corrected_site'] = resolved['corrected_site'].astype(str)
    resolved['corrected_site_key'] = resolved['canonical_UniProtAC'].astype(str) + ':' + resolved['corrected_site'].astype(str)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    save_table(remapped, outdir / 'human_arg_methyl_source_rows_with_remap.tsv')
    save_table(resolved, outdir / 'human_arg_methyl_source_rows_remapped_resolved.tsv')

    status = remapped['remap_status'].value_counts(dropna=False).rename_axis('remap_status').reset_index(name='row_count')
    save_table(status, outdir / 'remap_status_summary.tsv')

    summary = pd.DataFrame([{
        'rows_total': len(remapped),
        'rows_resolved': int(remapped['is_position_resolved'].fillna(False).sum()),
        'rows_with_position_shift': int(remapped['is_position_remapped'].fillna(False).sum()),
        'resolved_unique_sites': int(resolved['corrected_site_key'].nunique()),
        'resolved_canonical_accessions': int(resolved['canonical_UniProtAC'].nunique()),
    }])
    save_table(summary, outdir / 'summary.tsv')
    print(summary.to_dict(orient='records')[0])


if __name__ == '__main__':
    main()

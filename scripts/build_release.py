#!/usr/bin/env python3
"""Build the release tables of the corrected methylarginine union, with the provenance of every site and record.

Reads a corrected run made from the frozen inputs (remap rows, union, evidence check; see hpc/build_release.sh), the
published run (PTM_results/), the reference sequences and the input manifest, and writes to --outdir:

  sites.tsv                    one row per site of the release
  evidence.tsv                 one row per source record: what the study reported, where it was placed, and why
  changes_from_published.tsv   every site of the published union or the release: kept, removed or added, and why
  methylarginine_sites.xlsx    the three tables and a README sheet
  README.md                    what the files hold (columns, flags, reason codes) and how they were built
  release_manifest.tsv         size and MD5 of each file above
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
_spec = importlib.util.spec_from_file_location('remap', SCRIPTS / '02_remap_arg_methyl_to_canonical.py')
remap = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(remap)

# Studies and databases behind the union, with the DOIs the website lists (docs/assets/explorer.js STUDY_DOI).
STUDY_DOI = {
    'Maron 2021': '10.1016/j.isci.2021.102971', 'ProMetheusDB': '10.1016/j.mcpro.2022.100243',
    'PhosphoSitePlus': '10.1093/nar/gku1267', 'Larsen 2016': '10.1126/scisignal.aaf7329',
    'Musiani 2019': '10.1126/scisignal.aat8388', 'Fong 2019': '10.1016/j.ccell.2019.07.003',
    'Guo 2014': '10.1074/mcp.O113.027870', 'Fedoriw 2019': '10.1016/j.ccell.2019.05.014',
    'Wei 2020': '10.3389/fcell.2020.545934', 'Lim 2020': '10.1186/s12953-020-00162-8',
    'Li 2021': '10.1038/s41467-021-21963-1', 'IEDB': '10.1093/nar/gky1006', 'UniProt': '10.1093/nar/gkae1010',
}
SOURCE_TABLE = {'maron_s5': 'Maron 2021 Table S5', 'prometheus_mmc2': 'ProMetheusDB mmc2', 'human_ptm_master': 'PTM master table'}

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

# Why a record's placement differs from the published run; numbers are the issues of DATA_PROVENANCE_REVIEW.md.
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


def evidence_summary(row: pd.Series) -> tuple[str, str]:
    """Kind and text of the record's sequence evidence, as given by the source."""
    for field, kind in (('peptides', 'marked peptide'), ('seq_window', 'window'), ('peptide_sequence', 'peptide')):
        text = clean(row.get(field, ''))
        if text:
            length = int(float(row.get(f'{field}_window_length', 0) or 0)) if field != 'peptides' else 0
            return (f'centred window ({length} aa)' if length else kind), text
    return 'position only', ''


def record_reason(row: pd.Series, old: str, new: str, reference: dict, typed: pd.Series) -> str:
    """Reason code for a record whose placement differs between the published run and the release."""
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


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--release-root', required=True, help='results root of the corrected run (remapped/, integrated/, qc/)')
    ap.add_argument('--published-root', default='PTM_results')
    ap.add_argument('--inputs', default='inputs', help='frozen inputs with MANIFEST.tsv')
    ap.add_argument('--allow', default='config/remap_evidence_allowlist.tsv')
    ap.add_argument('--release-name', required=True)
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    root, pub_root, inputs, outdir = Path(args.release_root), Path(args.published_root), Path(args.inputs), Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    reference = remap.parse_fasta(inputs / 'uniprot_human_reviewed_canonical.fasta')
    meta = read_tsv(inputs / 'uniprot_human_reviewed_canonical_metadata.tsv').set_index('Entry')
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

    # The release and the published run read the same source records in the same order.
    key_cols = ['source_table', 'substrate_UniProtAC_original', 'position', 'peptide_sequence', 'seq_window', 'peptides']
    if len(rows) != len(pub) or not (rows[key_cols].values == pub[key_cols].values).all():
        raise SystemExit('the release run and the published run do not hold the same source records in the same order')
    if len(check) != len(rows):
        raise SystemExit('the evidence check does not cover every record of the release run')

    # ---- Records -----------------------------------------------------------------------------------------------
    rec = []
    for i, r in rows.iterrows():
        study = study_of(r)
        kind, text = evidence_summary(typed.loc[i])
        new, old = placed_key(r), placed_key(pub.loc[i])
        verdict = check.loc[i, 'verdict']
        exception = allowed.get((r['source_table'], check.loc[i, 'reported_accession'], check.loc[i, 'reported_position'], check.loc[i, 'placed_site']), '')
        change = 'same' if new == old and new else 'not placed in either' if not new and not old else 'moved' if new and old else 'newly placed' if new else 'no longer placed'
        reason = '' if change in ('same', 'not placed in either') else record_reason(r, old, new, reference, typed.loc[i])
        rec.append({
            'record_id': i + 1,
            'source_table': SOURCE_TABLE.get(r['source_table'], r['source_table']),
            'study': study,
            'study_doi': STUDY_DOI.get(study, ''),
            'reference_text': clean(r['reference']),
            'pmid': clean(r['pmid']),
            'reported_accession': clean(r['substrate_UniProtAC_original']),
            'reported_gene': clean(r['substrate_genename']),
            'reported_position': clean(r['position']),
            'reported_methyl_state': clean(r['arg_methyl_state']),
            'evidence_kind': kind,
            'evidence_sequence': text,
            'placed_accession': new.split(':')[0] if new else '',
            'placed_position': new.split(':')[1] if new else '',
            'placement': PLACEMENT.get(r['remap_status'], r['remap_status']),
            'placement_code': r['remap_status'],
            'evidence_check': verdict,
            'accepted_exception': exception,
            'published_site': old,
            'change_from_published': change,
            'reason_code': reason,
            'issue': REASONS[reason][0] if reason else '',
        })
    records = pd.DataFrame(rec)

    # ---- Sites -------------------------------------------------------------------------------------------------
    union_keys = {f'{a}:{int(float(p))}' for a, p in zip(union['canonical_UniProtAC'], union['corrected_position'])}
    pending = junction_sites(rows, typed, reference, union_keys)
    records['site'] = [f'{a}:{p}' if a else '' for a, p in zip(records['placed_accession'], records['placed_position'])]
    for site, ids in pending.items():
        records.loc[ids, 'site'] = site
    pub_keys = {f'{a}:{int(float(p))}' for a, p in zip(pub_union['canonical_UniProtAC'], pub_union['corrected_position'])}
    state = {f'{a}:{int(float(p))}': (c, o) for a, p, c, o in zip(union['canonical_UniProtAC'], union['corrected_position'], union['arg_methyl_state_consensus'], union['arg_methyl_states_observed'])}

    sites = []
    for site, group in records[records['site'] != ''].groupby('site', sort=False):
        acc, pos = site.split(':')[0], int(site.split(':')[1])
        seq = reference.get(acc, '')
        studies = sorted(set(group['study']))
        flags = []
        if site in pending:
            flags.append(PENDING_JUNCTION)
        if (group['accepted_exception'] != '').any():
            flags.append('paralog_window')
        if len(studies) == 1:
            flags.append('single_study')
        if site not in pending and group['placement_code'].str.contains('near_match').all():
            flags.append('near_match_only')
        consensus, observed = state.get(site, ('', ''))
        sites.append({
            'site': f'{acc}:R{pos}',
            'accession': acc,
            'entry_name': clean(meta['Entry Name'].get(acc, '')) if 'Entry Name' in meta else '',
            'gene': clean(meta['Gene Names (primary)'].get(acc, '')),
            'protein_name': clean(meta['Protein names'].get(acc, '')).split(' (')[0] if 'Protein names' in meta else '',
            'position': pos,
            'residue': seq[pos - 1] if seq and pos <= len(seq) else '',
            'window_31': ''.join(seq[k] if 0 <= k < len(seq) else '_' for k in range(pos - 16, pos + 15)),
            'protein_length': len(seq),
            'reference_sequence_md5': hashlib.md5(seq.encode()).hexdigest().upper() if seq else '',
            'release_status': 'pending policy decision: not placed by the pipeline' if site in pending else 'placed',
            'flags': ';'.join(flags),
            'methyl_state': clean(consensus),
            'methyl_states_reported': clean(observed),
            'n_records': len(group),
            'n_studies': len(studies),
            'studies': '; '.join(studies),
            'study_dois': '; '.join(STUDY_DOI[s] for s in studies if s in STUDY_DOI),
            'source_tables': '; '.join(sorted(set(group['source_table']))),
            'n_records_at_reported_position': int(((group['reported_accession'].map(remap.canonicalize_uniprot_accession) == acc) & (group['reported_position'] == str(pos))).sum()),
            'in_published_union': site in pub_keys,
            'record_ids': ';'.join(map(str, group['record_id'])),
        })
    sites = pd.DataFrame(sites).sort_values(['gene', 'accession', 'position']).reset_index(drop=True)

    # ---- Changes from the published union ----------------------------------------------------------------------
    release_keys = {s.replace(':R', ':') for s in sites['site']}
    changes = []
    for site in sorted(pub_keys | release_keys):
        in_pub, in_rel = site in pub_keys, site in release_keys
        acc, pos = site.split(':')
        if in_pub and in_rel:
            change, involved = 'kept', records.iloc[0:0]
        elif in_pub:
            change, involved = 'removed', records[(records['published_site'] == site) & (records['site'] != site)]
        else:
            change, involved = 'added', records[(records['site'] == site) & (records['published_site'] != site)]
        codes = sorted(set(involved['reason_code']) - {''})
        if change == 'added' and site in pending:
            codes = [PENDING_JUNCTION] + codes
        elsewhere = sorted(set(involved['site'] if change == 'removed' else involved['published_site']) - {'', site})
        changes.append({
            'site': f'{acc}:R{pos}',
            'accession': acc,
            'gene': clean(meta['Gene Names (primary)'].get(acc, '')) if acc in meta.index else '',
            'position': int(pos),
            'in_published_union': in_pub,
            'in_release': in_rel,
            'change': change,
            'reason_codes': ';'.join(codes),
            'issues': ';'.join(sorted({REASONS[c][0] for c in codes if c in REASONS and REASONS[c][0]}, key=int)),
            'other_sites': ';'.join(s.replace(':', ':R') for s in elsewhere),
            'record_ids': ';'.join(map(str, involved['record_id'])),
        })
    changes = pd.DataFrame(changes)

    records = records.drop(columns=['site'])
    records['placed_site'] = [f'{a}:R{p}' if a else '' for a, p in zip(records['placed_accession'], records['placed_position'])]
    records['published_site'] = records['published_site'].str.replace(':', ':R', regex=False)

    # ---- Write -------------------------------------------------------------------------------------------------
    sites.to_csv(outdir / 'sites.tsv', sep='\t', index=False)
    records.to_csv(outdir / 'evidence.tsv', sep='\t', index=False)
    changes.to_csv(outdir / 'changes_from_published.tsv', sep='\t', index=False)
    counts = write_readme(outdir, args.release_name, sites, records, changes, pending, inputs, root)
    with pd.ExcelWriter(outdir / 'methylarginine_sites.xlsx', engine='openpyxl') as xl:
        pd.DataFrame({'README': (outdir / 'README.md').read_text().splitlines()}).to_excel(xl, sheet_name='README', index=False)
        for name, table in (('sites', sites), ('evidence', records), ('changes_from_published', changes)):
            table.to_excel(xl, sheet_name=name, index=False, freeze_panes=(1, 0))
            sheet = xl.sheets[name]
            for k, col in enumerate(table.columns, start=1):
                width = min(60, max(len(str(col)), int(table[col].astype(str).str.len().quantile(0.9)) if len(table) else 0) + 2)
                sheet.column_dimensions[sheet.cell(row=1, column=k).column_letter].width = width
        xl.sheets['README'].column_dimensions['A'].width = 140
    manifest = pd.DataFrame([{'file': p.name, 'bytes': p.stat().st_size, 'md5': md5_bytes(p)} for p in sorted(outdir.iterdir())
                             if p.is_file() and p.name != 'release_manifest.tsv'])
    manifest.to_csv(outdir / 'release_manifest.tsv', sep='\t', index=False)
    print(counts)


def write_readme(outdir: Path, name: str, sites: pd.DataFrame, records: pd.DataFrame, changes: pd.DataFrame,
                 pending: dict, inputs: Path, root: Path) -> dict:
    placed = sites[sites['release_status'] == 'placed']
    flags = sites['flags'].str.split(';').explode().value_counts()
    counts = {
        'sites_placed': len(placed), 'proteins_placed': placed['accession'].nunique(),
        'sites_pending': len(pending), 'records': len(records),
        'records_placed': int((records['placed_accession'] != '').sum()),
        'kept': int((changes['change'] == 'kept').sum()), 'removed': int((changes['change'] == 'removed').sum()),
        'added': int((changes['change'] == 'added').sum()),
    }
    reason_counts = records['reason_code'].replace('', pd.NA).dropna().value_counts()
    try:
        commit = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True, text=True, cwd=SCRIPTS).stdout.strip()
    except OSError:
        commit = 'unknown'
    run_info = (root / 'RUN_INFO.txt').read_text() if (root / 'RUN_INFO.txt').exists() else ''
    repo = SCRIPTS.parent
    inputs = Path(os.path.relpath(inputs.resolve(), repo)) if inputs.resolve().is_relative_to(repo) else inputs
    run_info = run_info.replace(str(repo) + '/', '')
    lines = [
        f'# Human methylarginine sites: release {name}',
        '',
        'The corrected union of human protein arginine methylation sites behind DeAngelo, Silverstein & Shechter, '
        'Cell Chem Biol 33:1404-1421 (2026), doi:10.1016/j.chembiol.2026.09.011, with the provenance of every site and of '
        'every source record. It replaces the published union (PTM_results/, 13,719 sites on 5,090 accessions) after the '
        'remapping corrections described in DATA_PROVENANCE_REVIEW.md; changes_from_published.tsv lists every difference.',
        '',
        '## Contents',
        '',
        f"- {counts['sites_placed']:,} sites on {counts['proteins_placed']:,} reviewed human proteins, placed on the reference "
        'sequences below; plus ' + f"{counts['sites_pending']} sites included pending a policy decision (flag {PENDING_JUNCTION}).",
        f"- {counts['records']:,} source records, {counts['records_placed']:,} of them placed.",
        f"- Against the published union: {counts['kept']:,} sites kept, {counts['removed']:,} removed, {counts['added']:,} added.",
        '',
        '## Files',
        '',
        '| File | One row per | Use it for |',
        '|---|---|---|',
        '| sites.tsv | site | the site list, with where each site comes from |',
        '| evidence.tsv | source record (a study\'s report of a site) | what each study reported and where it was placed |',
        '| changes_from_published.tsv | site of the published union or of this release | what changed and why |',
        '| methylarginine_sites.xlsx | the three tables and this README | the same, as one workbook |',
        '| release_manifest.tsv | file | size and MD5 of each file |',
        '',
        '## Where the sites come from',
        '',
        'Three source tables: Maron et al. 2021 Table S5 (a compilation of ten studies: Larsen 2016, PhosphoSitePlus, '
        'Musiani 2019, Guo 2014, Fong 2019, Maron 2021\'s own data, Li 2021, Fedoriw 2019, Lim 2020, Wei 2020), ProMetheusDB '
        '(Massignani et al. 2022) and the arginine methylation rows of the PTM master table (UniProt and IEDB records). '
        'Each record\'s peptide, window or position was placed on the reviewed canonical UniProt sequence by '
        'scripts/02_remap_arg_methyl_to_canonical.py and checked by scripts/qc_remap_evidence_check.py.',
        '',
        '## sites.tsv columns',
        '',
        '- site, accession, entry_name, gene, protein_name, position, residue: the site on the reviewed canonical entry.',
        '- window_31: residues -15 to +15 of the reference sequence (_ beyond the ends).',
        '- protein_length, reference_sequence_md5: the reference sequence the position refers to (MD5 in UniProt\'s '
        'sequence.md5 format, to compare with the current UniProt entry).',
        '- release_status: placed, or pending a policy decision (see flags).',
        '- flags: single_study (one study reports the site); paralog_window (a supporting window matches a paralog '
        'exactly and this entry within 2 residues; issue 9, open); near_match_only (every supporting record matches the '
        f'entry within 2 residues, not exactly); {PENDING_JUNCTION} (a window crossing a splice junction: one full flank '
        'matches exactly at the reported arginine, the other does not; the pipeline leaves it unplaced, the release '
        'includes it until the policy is decided).',
        '- methyl_state, methyl_states_reported: methylation state calls, where a source gives them.',
        '- n_records, n_studies, studies, study_dois, source_tables: the support of the site.',
        '- n_records_at_reported_position: records whose own accession and position name this site unchanged.',
        '- in_published_union: whether the published union had this site.',
        '- record_ids: the supporting rows of evidence.tsv.',
        '',
        '## evidence.tsv columns',
        '',
        '- record_id; source_table; study, study_doi, reference_text, pmid: who reported it.',
        '- reported_accession, reported_gene, reported_position, reported_methyl_state: what the source reported.',
        '- evidence_kind, evidence_sequence: the window, peptide or marked peptide the source gave, verbatim.',
        '- placed_accession, placed_position, placed_site, placement, placement_code: where the release puts it and how.',
        '- evidence_check: the independent check (agrees; contradicts evidence; one of several evidence sites; ...).',
        '- accepted_exception: why a contradiction is accepted for now (open issue 9).',
        '- published_site, change_from_published, reason_code, issue: where the published run put it and why it differs.',
        '',
        '## Reason codes',
        '',
        '| reason_code | issue | records | meaning |',
        '|---|---|---|---|',
    ]
    for code, (issue, text) in REASONS.items():
        lines.append(f'| {code} | {issue} | {int(reason_counts.get(code, 0)):,} | {text} |')
    lines += [
        '',
        f'Flags in sites.tsv: ' + ', '.join(f'{k} {v:,}' for k, v in flags.items() if k) + '.',
        '',
        '## How it was built',
        '',
        f'- Code: PTM_Pipeline commit {commit} (scripts/02, scripts/qc_remap_evidence_check.py, scripts/build_release.py; '
        'driver hpc/build_release.sh).',
        f'- Inputs: {inputs}/ (checksums in {inputs}/MANIFEST.tsv). Reference sequences: UniProtKB reviewed human '
        '(proteome UP000005640), staged 2026-04-07; release not recorded, earlier than 2026_03. The published union was '
        'mapped to the same sequences. Accession history: UniProt 2026_03.',
        '- Checks passed: every record agrees with its own evidence except the accepted issue-9 rows in '
        'config/remap_evidence_allowlist.tsv; every placed site is an arginine of the reference sequence.',
        '',
        '## Licences',
        '',
        'PhosphoSitePlus® data (www.phosphosite.org) are licensed under CC BY-NC-SA 3.0: non-commercial use, with '
        'attribution, adaptations shared under the same licence. The release as a whole is CC BY-NC-SA 4.0 (LICENSE.md). '
        'Cite the studies in study_dois for the sites you use.',
    ]
    if run_info:
        lines += ['', '## Run record', '', '```', run_info.strip(), '```']
    (outdir / 'README.md').write_text('\n'.join(lines) + '\n')
    return counts


if __name__ == '__main__':
    main()

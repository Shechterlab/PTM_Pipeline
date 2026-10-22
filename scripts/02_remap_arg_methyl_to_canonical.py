from __future__ import annotations

import argparse
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

import release_codes as rc
from common import (
    accession_is_isoform,
    canonical_gene_name,
    canonicalize_uniprot_accession,
    clean_aa_sequence,
    mass_localised,
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
    # Every primary gene of an entry: UniProt lists several, separated by ';', for some (histones, SLX1A; SLX1B).
    genes_of = meta.get('gene_primary', pd.Series('', index=meta.index)).fillna('').astype(str).str.split(';')
    meta['gene_primary'] = genes_of.str[0].map(canonical_gene_name)
    gene_to_accessions: dict[str, list[str]] = defaultdict(list)
    for acc, genes in zip(meta['canonical_UniProtAC'], genes_of):
        for gene in {canonical_gene_name(g) for g in genes}:
            if gene and acc:
                gene_to_accessions[gene].append(acc)
    gene_to_accessions = {gene: sorted(set(accs)) for gene, accs in gene_to_accessions.items()}
    return meta, gene_to_accessions


def load_accession_history(path: str | Path | None) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    """(secondary accession -> current entries, gene synonym -> entries) from uniprot_accession_history.tsv (scripts/00).
    A merged accession maps to one entry; an accession split into several maps to each of them."""
    secondary: dict[str, list[str]] = defaultdict(list)
    synonyms: dict[str, list[str]] = defaultdict(list)
    if not path:
        return {}, {}
    table = pd.read_csv(path, sep='\t', dtype=str).fillna('')
    for kind, key, acc in zip(table['kind'], table['key'], table['accession']):
        target = secondary if kind == 'secondary' else synonyms if kind == 'synonym' else None
        key = key.strip() if kind == 'secondary' else canonical_gene_name(key)
        if target is not None and key and acc and acc not in target[key]:
            target[key].append(acc)
    return {k: sorted(v) for k, v in secondary.items()}, {k: sorted(v) for k, v in synonyms.items()}


class ProteomeIndex:
    """Substring search over every reference sequence at once (for rows whose candidates hold none of their evidence)."""

    def __init__(self, sequences: dict[str, str]):
        self.sequences = sequences
        self.order = sorted(sequences)
        self.offsets: list[int] = []
        parts, off = [], 0
        for acc in self.order:
            self.offsets.append(off)
            parts.append(sequences[acc] + '#')
            off += len(sequences[acc]) + 1
        self.text = ''.join(parts)

    def find(self, token: str) -> list[tuple[str, int]]:
        """(accession, 1-based start) of every occurrence."""
        import bisect
        out, i = [], self.text.find(token)
        while i != -1:
            k = bisect.bisect_right(self.offsets, i) - 1
            out.append((self.order[k], i - self.offsets[k] + 1))
            i = self.text.find(token, i + 1)
        return out


# Evidence searched across the whole reference set must be long enough to name one place.
PROTEOME_MIN_ANCHOR = 7


def own_accession(row: pd.Series, canonical_sequences: dict[str, str], secondary: dict[str, list[str]],
                  gene_map: dict[str, list[str]] | None = None) -> str:
    """The row's own entry: its accession when that is in the reference set, else the current entry UniProt merged it
    into. A split accession (HSPA1A P08107 -> P0DMV8, P0DMV9) has several; the one named by the row's gene symbol is its
    own, and with none or several the row keeps no own entry."""
    acc = canonicalize_uniprot_accession(row.get('canonical_UniProtAC', row.get('substrate_UniProtAC', '')))
    if acc in canonical_sequences:
        return acc
    current = [a for a in secondary.get(acc, []) if a in canonical_sequences]
    if len(current) > 1:
        named = [a for a in current if a in (gene_map or {}).get(canonical_gene_name(row.get('substrate_genename', '')), [])]
        current = named
    return current[0] if len(current) == 1 else acc


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


def study_frame(path: str | Path) -> pd.DataFrame:
    """A study's records as scripts/00_stage_study_sites.py stages them (inputs/sources/<study>_methylarginine_records.tsv),
    in the columns of the other source tables: the marked peptide is read like ProMetheusDB's, the window like Table
    S5's. The source table is named after the study in lower case (hartel2019)."""
    t = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    study = t['study_id'].iloc[0]
    if (t['study_id'] != study).any():
        raise SystemExit(f'{path}: one study per file')
    out = pd.DataFrame({
        'substrate_UniProtAC_original': t['reported_accession'],
        'substrate_UniProtAC': t['reported_accession'].str.split('; ').str[0],
        'substrate_genename': t['reported_gene'],
        'position': t['reported_position'],
        'residue': 'R',
        # marked as ProMetheusDB marks them, R(me) and R(di), which the evidence check reads too
        'peptides': t['peptide'].str.replace('R(me1)', 'R(me)', regex=False).str.replace('R(me2)', 'R(di)', regex=False),
        'seq_window': t['window'],
        # The degree, or for a dimethylarginine the one type a study called from neutral losses (Hartel 2019); a
        # residue with both types called stays Rme2 here, and study_rme2_type keeps the call.
        'arg_methyl_state': [{'ADMA': 'Rme2a', 'SDMA': 'Rme2s'}.get(x, d) for d, x in zip(t['methyl_degree'], t['rme2_type'])],
        'study_rme2_type': t['rme2_type'],
        'arg_methyl_state_known': t['methyl_degree'].ne(''),
        'site_localisation': t['localised_by'],
        'study_not_localised': t['not_localised'] if 'not_localised' in t.columns else '',
        'source': study,
        'source_family': study,
        'reference': study,
        'ptm_group': 'Arg methylation',
        'data_layer': 'external_methyl_table',
        'study_sample': t['cell_type'],
        'study_enrichment': t['enrichment'],
        'study_localisation': t['localisation_probability'],
        'study_rows': t['n_source_rows'],
        'study_tables': t['source_tables'],
        'study_note': t['note'],
        'staged_study': True,
    })
    return prepare_frame(out, study.lower())


# Sources whose export writes every leucine as isoleucine (isobaric in mass spectrometry).
LEUCINE_AS_ISOLEUCINE = {'sylvestersen2014', 'lim2020'}


def restore_leucines(rows: pd.DataFrame, canonical_sequences: dict[str, str]) -> pd.DataFrame:
    """For a source of LEUCINE_AS_ISOLEUCINE, write back L where the export wrote I: in the window and the marked peptide
    of a row whose reported entry, read at the reported position, agrees with them once I and L are taken as one
    residue, or else, for a row whose reported entry is not in the reference set or does not agree, from the one site
    of the whole reference set its peptide (or window) names when I and L are taken as one. The row's evidence as
    reported goes to peptide_sequence_reported and source_correction says leucine_restored; other rows are left as
    given."""
    rows = rows.copy()
    for col in ('source_correction', 'peptide_sequence_reported'):
        if col not in rows.columns:
            rows[col] = ''
    targets = rows.index[rows['source_table'].isin(LEUCINE_AS_ISOLEUCINE)]
    collapsed = ProteomeIndex({acc: seq.replace('I', 'L') for acc, seq in canonical_sequences.items()}) if len(targets) else None
    for i in targets:
        r = rows.loc[i]
        seq = canonical_sequences.get(canonicalize_uniprot_accession(r['substrate_UniProtAC']), '')
        pos = int(r['position'])
        fixed = leucines_from_entry(r, seq, pos) if seq else {}
        if not fixed:
            fixed = leucines_from_proteome(r, canonical_sequences, collapsed)
        if fixed:
            peptide, window = str(r.get('peptides', '') or ''), str(r.get('seq_window', '') or '')
            rows.at[i, 'peptide_sequence_reported'] = ', '.join(x for x in (peptide, window) if x)
            rows.at[i, 'source_correction'] = 'leucine_restored'
            for col, value in fixed.items():
                rows.at[i, col] = value
    return rows


def with_leucines(marked: str, ref: str) -> str:
    """marked (a peptide, its marks kept, or a window, its padding kept) with L written where ref, aligned residue by
    residue, has L and marked has I."""
    out, k = [], 0
    for ch in marked:
        if 'A' <= ch <= 'Z' and ch != 'X' and k < len(ref):
            out.append('L' if ch == 'I' and ref[k] == 'L' else ch)
            k += 1
        else:
            out.append(ch)
    return ''.join(out)


def leucines_from_proteome(r: pd.Series, canonical_sequences: dict[str, str], collapsed: ProteomeIndex) -> dict[str, str]:
    """The row's peptide and window with L written back from the one site of the reference set its peptide (or, without
    one, its window) names when I and L are taken as one residue: each where it matches that site so; empty when the
    sequence names none or several sites, or nothing changes."""
    peptide, window = str(r.get('peptides', '') or ''), str(r.get('seq_window', '') or '')
    pieces = []
    if peptide:
        plain, marked = parse_modified_peptide(peptide)
        if len(marked) == 1:
            pieces.append(('peptides', peptide, plain, marked[0] - 1))
    if window:
        core, lead, _ = split_window_padding(window)
        if core:
            pieces.append(('seq_window', window, core, len(window) // 2 - lead))
    pieces = [x for x in pieces if len(x[2]) >= PROTEOME_MIN_ANCHOR]
    if not pieces:
        return {}
    _, _, plain, centre = pieces[0]
    sites = {(acc, start + centre) for acc, start in collapsed.find(plain.replace('I', 'L'))}
    if len(sites) != 1:
        return {}
    acc, pos = next(iter(sites))
    seq, fixed = canonical_sequences[acc], {}
    for col, given, plain, centre in pieces:
        start = pos - centre - 1
        ref = seq[start:start + len(plain)] if start >= 0 else ''
        if ref.replace('I', 'L') == plain.replace('I', 'L'):
            restored = with_leucines(given, ref)
            if restored != given:
                fixed[col] = restored
    return fixed


def leucines_from_entry(r: pd.Series, seq: str, pos: int) -> dict[str, str]:
    """The row's window and peptide with L written back from its reported entry, read at the reported position, where
    they agree with it once I and L are taken as one residue."""
    fixed = {}
    window = str(r.get('seq_window', '') or '')
    if window:
        half = len(window) // 2
        ref = ''.join(seq[k] if 0 <= k < len(seq) else '_' for k in range(pos - 1 - half, pos - 1 - half + len(window)))
        if ref.replace('I', 'L') == window.replace('I', 'L') and ref != window:
            fixed['seq_window'] = ''.join('L' if w == 'I' and c == 'L' else w for w, c in zip(window, ref))
    peptide = str(r.get('peptides', '') or '')
    if peptide:
        plain, marked = parse_modified_peptide(peptide)
        start = pos - marked[0] if len(marked) == 1 else -1
        ref = seq[start:start + len(plain)] if start >= 0 else ''
        if ref and ref.replace('I', 'L') == plain.replace('I', 'L') and ref != plain:
            fixed['peptides'] = with_leucines(peptide, ref)
    return fixed


def lead_reviewed_accession(rows: pd.DataFrame, canonical_sequences: dict[str, str], gene_map: dict[str, list[str]]) -> pd.DataFrame:
    """A staged study's record lists every accession of its protein group, led by the one its position numbers, which can
    be an unreviewed entry (Geoghegan 2015 numbers the alphabetically first). When that one is not a reviewed entry, the
    record is searched from a listed reviewed accession instead: the one of the record's gene, or else the first listed;
    the list as reported stays in substrate_UniProtAC_original."""
    rows = rows.copy()
    staged = rows['staged_study'].astype('boolean').fillna(False).astype(bool) if 'staged_study' in rows.columns else pd.Series(False, index=rows.index)
    for i in rows.index[staged]:
        listed = [canonicalize_uniprot_accession(a) for a in str(rows.at[i, 'substrate_UniProtAC_original']).split('; ')]
        if not listed or listed[0] in canonical_sequences:
            continue
        reviewed = [a for a in listed if a in canonical_sequences]
        if reviewed:
            own = set(gene_map.get(canonical_gene_name(rows.at[i, 'substrate_genename']), []))
            lead = next((a for a in reviewed if a in own), reviewed[0])
            rows.at[i, 'substrate_UniProtAC'] = lead
            rows.at[i, 'canonical_UniProtAC'] = lead
    return rows


def proteome_sites(row: pd.Series, proteome: ProteomeIndex, peptides: list[str] | None = None) -> set[tuple[str, int]]:
    """Every site the row's evidence names anywhere in the reference set: each occurrence of its marked peptide (of
    PROTEOME_MIN_ANCHOR residues or more), at the marked arginine. A row without one is read from its centred window
    instead: the peptide is what the study measured, the window only the sequence of the entry it assigned. peptides:
    read these marked peptides instead of the row's (the ones that localise the site)."""
    if peptides is not None:
        row = row.copy()
        row['peptides'], row['seq_window'], row['peptide_sequence'] = ', '.join(peptides), '', ''
    anchors = [a for a in build_anchors(row) if len(a['sequence']) >= PROTEOME_MIN_ANCHOR]
    marked = [a for a in anchors if a['anchor_type'] == 'annotated_peptide']
    sites = set()
    for anchor in marked or anchors:
        for acc, start in proteome.find(anchor['sequence']):
            sites |= {(acc, start + k - 1) for k in anchor['indices']}
    return sites


# The residues each search could methylate (variable modifications), for the sources whose records are localised here
# from the peptide's mass: ProMetheusDB (mmc2 header: mono- and dimethylation of R and K, trimethylation of K).
SEARCH_METHYLATABLE = {'prometheus_mmc2': 'RK'}


def marked_arginines(peptide: str) -> tuple[str, dict[int, str]]:
    """'K(di)NNNLYSWDDFM(ox)ELGR(me)' -> ('KNNNLYSWDDFMELGR', {16: 'Rme1'}): the plain peptide and its methylated
    arginines, R(me) mono and R(di) dimethyl."""
    plain, marks = [], {}
    for m in re.finditer(r'([A-Za-z])(?:\(([^)]*)\))?', peptide):
        plain.append(m.group(1).upper())
        if plain[-1] == 'R' and m.group(2):
            marks[len(plain)] = {'me': 'Rme1', 'di': 'Rme2'}.get(m.group(2).strip().lower(), '')
    return ''.join(plain), marks


def localising_peptides(peptides: str, methylatable: str) -> list[tuple[str, str]]:
    """The listed marked peptides whose mass alone localises their methyl groups, each with how (common.mass_localised)."""
    out = []
    for peptide in split_multi_value_field(peptides, separators=','):
        how = mass_localised(*marked_arginines(peptide), methylatable)
        if how:
            out.append((peptide, how))
    return out


# A record that adds no site lists at most this many candidate sites (lower-confidence evidence).
MAX_CANDIDATES = 12


def candidate_sites(r: pd.Series, proteome: ProteomeIndex, peptides: list[str] | None = None, residues: bool = False) -> list[str]:
    """Where a record that adds no site could put its methyl group, as site keys (lower-confidence evidence). residues
    False: the sites its evidence names in the reviewed proteome (the marked or centred arginine wherever its peptide or
    window occurs). residues True, for a peptide that leaves several residues: every arginine of the peptide, in the
    entry remap_row placed it in when it did, else wherever the peptide occurs. Without either, the site remap_row placed
    it on. At most MAX_CANDIDATES, in order."""
    placed, acc, pos = None, r['corrected_canonical_UniProtAC'], r['corrected_position']
    if str(r['is_position_resolved']) == 'True' and isinstance(acc, str) and acc and pd.notna(pos):
        placed = (acc, int(float(pos)))
    sites = set()
    if residues:
        plain = next((seq for seq, _ in (marked_arginines(t.strip()) for t in str(r.get('peptides', '') or '').split(','))
                      if len(seq) >= PROTEOME_MIN_ANCHOR), '')
        if plain:
            hits = proteome.find(plain)
            if placed:
                hits = [(a, st) for a, st in hits if a == placed[0] and st <= placed[1] < st + len(plain)] or hits
            sites = {(a, st + k) for a, st in hits for k, aa in enumerate(plain) if aa == 'R'}
    if not sites:
        sites = proteome_sites(r, proteome, peptides)
    if not sites and placed:
        sites = {placed}
    return [f'{a}:R{pos}' for a, pos in sorted(sites)[:MAX_CANDIDATES]]


def named_paralog(r: pd.Series, sites: set) -> bool:
    """The paralog rule: the evidence names one site in each of several reviewed entries (paralogs or identical copies
    share the sequence, none holds it twice), and remap_row placed the record on one of them, in an entry the record
    names (its accession, the current entry of a merged accession, or an entry of its gene). Mass spectrometry cannot
    tell the copies apart; the site is credited to the protein the study assigned the peptide to (MaxQuant's razor
    protein, for most studies) and flagged paralog_assigned."""
    accs = [a for a, _ in sites]
    if len(sites) < 2 or len(set(accs)) != len(accs):
        return False
    acc, pos = r['corrected_canonical_UniProtAC'], r['corrected_position']
    if str(r['is_position_resolved']) != 'True' or r['remap_status'] == 'sequence_unique_proteome' or not isinstance(acc, str) or pd.isna(pos):
        return False
    return (acc, int(float(pos))) in sites


def record_tests(r: pd.Series, proteome: ProteomeIndex, how: str, peptides: list[str] | None = None,
                 residues: bool = False, sites: set | None = None) -> dict:
    """The site tests of a record (release_codes.CODES['site_test']), at the site its study reports.

    reported_site: where remap_row placed the record from what the study reports (its accession or an entry it names,
    the reported position, the evidence), an arginine of the reference; empty when it could not. evidence_sites: every
    arginine its evidence names in the reviewed proteome (the marked or centred arginine wherever the peptide or window
    occurs; residues True, a peptide that leaves several residues: every arginine of every occurrence). Each test is
    pass, fail or untested:
      evidence_at_site       the peptide or window occurs in the reported entry and covers the reported arginine
      named_entry            the reported entry is one the record names (remap_row placed it there, not in another
                             entry the whole-proteome search found: sequence_unique_proteome)
      unique_in_proteome     the evidence names one site in the reviewed proteome (residues True: the peptide occurs
                             once); paralogs: one site in each of several entries, one of them the reported entry
                             the record names (named_paralog)
      methyl_on_arginine     the methyl group can only be on an arginine: localised, or the peptide holds no lysine
                             (fail: it holds one, which could carry the methyl group); untested without a peptide
      localised              the study localises the methyl group to the arginine (site_localisation)
    sites: proteome_sites(r, proteome, peptides) when the caller has it."""
    reported = ''
    acc, pos = r['corrected_canonical_UniProtAC'], r['corrected_position']
    if str(r['is_position_resolved']) == 'True' and isinstance(acc, str) and acc and pd.notna(pos):
        seq = proteome.sequences.get(acc, '')
        if 0 < int(float(pos)) <= len(seq) and seq[int(float(pos)) - 1] == 'R':
            reported = (acc, int(float(pos)))
    plains = [plain for plain, _ in (marked_arginines(t.strip()) for t in str(r.get('peptides', '') or '').split(',')) if plain]
    if not plains and str(r.get('peptide_sequence', '') or '') not in ('', 'nan'):
        plains = [re.sub(r'[^A-Z]', '', str(r['peptide_sequence']).upper())]
    if residues:
        plain = next((x for x in plains if len(x) >= PROTEOME_MIN_ANCHOR), '')
        hits = proteome.find(plain) if plain else []
        sites = {(a, st + k) for a, st in hits for k, aa in enumerate(plain) if aa == 'R'}
        unique = len(hits) == 1
    else:
        sites = proteome_sites(r, proteome, peptides) if sites is None else sites
        unique = len(sites) == 1
    curated = how in ('curated_experimental',) and not sites
    test = lambda ok: 'pass' if ok else 'fail'  # noqa: E731
    return {
        'reported_site': f'{reported[0]}:R{reported[1]}' if reported else '',
        'evidence_sites': '; '.join(f'{a}:R{p}' for a, p in sorted(sites)),
        'test_evidence_at_site': 'curated' if curated else test(bool(reported) and reported in sites),
        'test_named_entry': test(bool(reported) and r['remap_status'] != 'sequence_unique_proteome'),
        'test_unique_in_proteome': 'curated' if curated else 'paralogs' if not residues and not unique and named_paralog(r, sites) else test(unique),
        'test_methyl_on_arginine': 'pass' if how else (test('K' not in ''.join(plains)) if plains else 'untested'),
        'test_localised': test(bool(how)),
    }


TEST_COLUMNS = ['reported_site', 'evidence_sites', 'test_evidence_at_site', 'test_named_entry', 'test_unique_in_proteome',
                'test_methyl_on_arginine', 'test_localised']


def unplace(remapped: pd.DataFrame, i, status: str) -> None:
    remapped.loc[i, ['is_position_resolved', 'corrected_canonical_UniProtAC', 'corrected_position', 'corrected_site']] = [False, '', None, '']
    remapped.at[i, 'is_position_remapped'] = False
    remapped.at[i, 'remap_status'] = status


def place_localised_records(remapped: pd.DataFrame, proteome: ProteomeIndex) -> pd.DataFrame:
    """A record adds a site only where its source localises the methyl group to one residue and its evidence names
    that residue unequivocally.

    Localised (site_localisation, codes in release_codes.CODES['localised_by']): a staged study's record as its staging
    says (a localisation score over the study's threshold, a peptide whose mass leaves one arrangement, or a curated
    experimental annotation); a ProMetheusDB record when one of its marked peptides localises it by mass; a record of a
    curated database trusted as curated (release_codes.CURATED_DATABASES: PhosphoSitePlus, curated_database). Any other
    record of a compiled table (Maron 2021 Table S5, the iPTMnet-format master table) carries no localisation: it is not
    placed, superseded_by_study_table when the study's own table is read in this run, else unresolved_not_localised.

    Placed: a curated record (UniProt) where remap_row placed it (its accession, at an arginine of that position). Any
    other localised record stays placed when its localising peptides, or without a peptide its window, searched in the
    whole reviewed proteome, name exactly one site, and remap_row placed it there in an entry the record names (its
    accession, the current entry of a merged accession, or an entry of its gene). A sequence that names one site in each
    of several entries (paralogs, identical copies) stays placed on the one remap_row chose, the entry the record names:
    sequence_paralog_named (named_paralog), then all records of that sequence on one copy (paralog_consensus; the
    release flags the site paralog_assigned). Repeats within a protein,
    isoform-specific or changed sequences, a sequence found only in another protein (which the study may have seen in an
    isoform or fusion the reviewed proteome does not hold), near matches and position-only records stay unplaced."""
    remapped = remapped.copy()
    if 'site_localisation' not in remapped.columns:
        remapped['site_localisation'] = ''
    remapped['site_localisation'] = remapped['site_localisation'].fillna('').astype(str)
    remapped['candidate_sites'] = ''
    for column in TEST_COLUMNS:
        remapped[column] = ''
    staged = remapped['staged_study'].astype('boolean').fillna(False).astype(bool) if 'staged_study' in remapped.columns else pd.Series(False, index=remapped.index)
    read = {rc.STUDY_TABLES[t] for t in remapped.loc[staged, 'source_table'].unique()}
    for i in remapped.index:
        r = remapped.loc[i]
        peptides = None
        if staged[i]:
            how = r['site_localisation']
        elif r['source_table'] in SEARCH_METHYLATABLE:
            found = localising_peptides(str(r.get('peptides', '') or ''), SEARCH_METHYLATABLE[r['source_table']])
            how = found[0][1] if found else ''
            peptides = [p for p, _ in found]
            remapped.at[i, 'site_localisation'] = how
        elif rc.study_of(r['source_table'], r.get('source', ''), r.get('reference', '')) in rc.CURATED_DATABASES:
            how = 'curated_database'
            remapped.at[i, 'site_localisation'] = how
        else:
            how = ''
        if not how:
            superseded = not staged[i] and rc.study_of(r['source_table'], r.get('source', ''), r.get('reference', '')) in read
            why = str(r.get('study_not_localised', '') or '') if staged[i] else 'several_residues'
            residues = why in ('several_residues', 'search_disagrees')
            remapped.loc[i, TEST_COLUMNS] = list(record_tests(r, proteome, how, residues=residues).values())
            if not superseded:   # lower-confidence evidence: where the methyl group could be
                remapped.at[i, 'candidate_sites'] = '; '.join(candidate_sites(r, proteome, residues=residues))
            unplace(remapped, i, 'superseded_by_study_table' if superseded else 'unresolved_not_localised')
            continue
        sites = proteome_sites(r, proteome, peptides)
        remapped.loc[i, TEST_COLUMNS] = list(record_tests(r, proteome, how, peptides, sites=sites).values())
        if how == 'curated_experimental':
            continue
        if len(sites) == 1:
            acc, pos = next(iter(sites))
            placed = r['is_position_resolved'] and r['remap_status'] != 'sequence_unique_proteome'
            if placed and (r['corrected_canonical_UniProtAC'], int(float(r['corrected_position']))) == (acc, pos):
                # The code says what placed it: the sequence found once, in the reported entry or another the record names
                # (remap_row may have reached the same site from a window that occurs more than once).
                own = acc == canonicalize_uniprot_accession(r['substrate_UniProtAC'])
                remapped.at[i, 'remap_status'] = 'sequence_unique_canonical' if own else 'sequence_unique_candidate'
                continue
            status = 'unresolved_study_other_entry'
        elif named_paralog(r, sites):
            remapped.at[i, 'remap_status'] = 'sequence_paralog_named'
            continue
        else:
            status = 'unresolved_study_not_unique' if sites else 'unresolved_study_not_found'
        remapped.at[i, 'candidate_sites'] = '; '.join(candidate_sites(r, proteome, peptides))
        unplace(remapped, i, status)
    return paralog_consensus(remapped)


def paralog_consensus(remapped: pd.DataFrame) -> pd.DataFrame:
    """The paralog rule gives one site per shared sequence: the records whose evidence names the same sites (one in each
    paralog) go on the copy most of them name, then the copy most studies name, then the first accession. Studies name
    different copies of the same sequence (their razor proteins differ); the consensus is the copy the evidence favours,
    and the site is flagged paralog_assigned with the others."""
    idx = remapped.index[remapped['remap_status'] == 'sequence_paralog_named']
    groups = defaultdict(list)
    for i in idx:
        groups[remapped.at[i, 'evidence_sites']].append(i)
    for members in groups.values():
        named, studies = Counter(), defaultdict(set)
        for i in members:
            r = remapped.loc[i]
            site = (r['corrected_canonical_UniProtAC'], int(float(r['corrected_position'])))
            named[site] += 1
            studies[site].add(rc.study_of(r['source_table'], r.get('source', ''), r.get('reference', '')))
        acc, pos = sorted(named, key=lambda s: (-named[s], -len(studies[s]), s))[0]
        for i in members:
            remapped.loc[i, ['corrected_canonical_UniProtAC', 'corrected_position', 'corrected_site']] = [acc, pos, f'R{pos}']
    return remapped


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


def candidate_accessions(row: pd.Series, canonical_sequences: dict[str, str], gene_map: dict[str, list[str]],
                         secondary: dict[str, list[str]] | None = None, synonyms: dict[str, list[str]] | None = None) -> list[str]:
    """Entries to search for the row's evidence: its own accession(s); for an accession UniProt has since merged or split,
    the current entries (secondary); entries sharing its gene symbol; and, when the symbol names no entry, entries that
    list it as a synonym (outdated symbols such as MST4 for STK26)."""
    candidates: list[str] = []
    secondary, synonyms = secondary or {}, synonyms or {}

    def add(acc: str) -> None:
        if acc and acc in canonical_sequences and acc not in candidates:
            candidates.append(acc)

    for field in ['canonical_UniProtAC', 'substrate_UniProtAC', 'UniProtAC', 'substrate_UniProtAC_original']:
        if field not in row.index:
            continue
        acc = canonicalize_uniprot_accession(row.get(field, ''))
        add(acc)
        if acc and acc not in canonical_sequences:
            for current in secondary.get(acc, []):
                add(current)
    gene = canonical_gene_name(row.get('substrate_genename', ''))
    for acc in gene_map.get(gene, []):
        add(acc)
    if gene and not any(acc in canonical_sequences for acc in gene_map.get(gene, [])):
        for acc in synonyms.get(gene, []):
            add(acc)
    return candidates


# A source record gives either a centred site window (13, 15 or 31 residues: Guo; PhosphoSitePlus, Musiani, Fong,
# ProMetheusDB; Larsen, Wei) or a peptide that contains the site somewhere (Maron, Fedoriw, Li, Lim). The token length
# alone does not tell them apart (a tryptic peptide can be 15 residues long), so the kind is read per source: a source
# whose tokens mostly share one window length gives windows. Its windows are shorter where they run off a protein end.
WINDOW_LENGTHS = (13, 15, 31)
WINDOW_SOURCE_MIN_SHARE = 0.5


def annotate_window_lengths(rows: pd.DataFrame) -> pd.DataFrame:
    """Add seq_window_window_length and peptide_sequence_window_length: the full window length of the row's source for
    that field, or 0 when the source gives peptides. A source is a source table plus, for compilations, the reference."""
    rows = rows.copy()
    reference = rows['reference'] if 'reference' in rows.columns else pd.Series('', index=rows.index)
    keys = rows['source_table'].astype(str) + '|' + reference.fillna('').astype(str)
    for field in ('seq_window', 'peptide_sequence'):
        values = rows[field] if field in rows.columns else pd.Series('', index=rows.index)
        lengths = values.map(lambda v: [len(clean_aa_sequence(t)) for t in split_multi_value_field(v, separators=';')])
        per_source = {}
        for key, group in lengths.groupby(keys):
            flat = [n for row_lengths in group for n in row_lengths]
            if flat:
                top, count = Counter(flat).most_common(1)[0]
                per_source[key] = top if top in WINDOW_LENGTHS and count / len(flat) >= WINDOW_SOURCE_MIN_SHARE else 0
        rows[f'{field}_window_length'] = keys.map(per_source).fillna(0).astype(int)
    return rows


# Guo pads its windows with X and ProMetheusDB with _ where they run off a protein end; the padding says which end.
WINDOW_PADDING = '_-X'


def split_window_padding(token: str) -> tuple[str, int, int]:
    """Residues of a window token and the number of padding characters before and after them."""
    raw = ''.join(ch for ch in str(token).upper() if 'A' <= ch <= 'Z' or ch in WINDOW_PADDING)
    core = raw.strip(WINDOW_PADDING)
    if not core:
        return '', 0, 0
    return clean_aa_sequence(core), len(raw) - len(raw.lstrip(WINDOW_PADDING)), len(raw) - len(raw.rstrip(WINDOW_PADDING))


def truncated_window_indices(window: str, full_length: int, start: int, sequence_length: int, reported: int) -> list[int]:
    """1-based index of the site in a window cut short without padding, for its occurrence at start (1-based) in a
    sequence. The site is half a window from the intact end, so the protein end the occurrence touches decides the
    reading. An occurrence inside the sequence (the source protein began or ended there: another start codon, a mature
    chain, an isoform) allows both; when both give an arginine, only one that the reported position names is taken,
    counted either in this sequence or, for a window cut at its N-terminal side, from the window's first residue."""
    half = full_length // 2
    n_cut, c_cut = len(window) - half, half + 1
    end = start + len(window) - 1
    at_n, at_c = start == 1, end == sequence_length
    if at_n and at_c:
        readings = range(max(1, n_cut), min(len(window), c_cut) + 1)
    elif at_n or at_c:
        readings = [n_cut] if at_n else [c_cut]
    else:
        readings = [n_cut, c_cut]
    readings = [i for i in readings if 1 <= i <= len(window) and window[i - 1] == 'R']
    if len(readings) > 1:
        readings = [i for i in readings if start + i - 1 == reported or (i == n_cut and i == reported)]
    return readings if len(readings) == 1 else []


def row_evidence(row: pd.Series) -> list[dict]:
    """Every sequence the source gives for a row, with the 0-based offsets it names as the site: the marked arginines
    of an annotated peptide, the centre of a window (whatever residue it is), the arginines of a plain peptide. A window
    cut short without padding has no fixed offset ('full_length' set instead): its reading depends on where it lies."""
    evidence: list[dict] = []
    for peptide in split_multi_value_field(row.get('peptides', ''), separators=','):
        clean, modified = parse_modified_peptide(peptide)
        if clean:
            evidence.append({'sequence': clean, 'offsets': [i - 1 for i in modified], 'full_length': 0, 'anchor_type': 'annotated_peptide'})
    for field in ('seq_window', 'peptide_sequence'):
        full_length = int(row.get(f'{field}_window_length', 0) or 0)
        half = full_length // 2
        for token in split_multi_value_field(row.get(field, ''), separators=';'):
            core, lead, trail = split_window_padding(token) if full_length else (clean_aa_sequence(token), 0, 0)
            if not core:
                continue
            padded = (lead or trail) and lead + len(core) + trail == full_length
            if full_length and (len(core) == full_length or padded):
                kind = 'centered_window' if len(core) == full_length else 'truncated_window'
                evidence.append({'sequence': core, 'offsets': [half - lead], 'full_length': 0, 'anchor_type': kind})
            elif full_length and len(core) < full_length:
                evidence.append({'sequence': core, 'offsets': None, 'full_length': full_length, 'anchor_type': 'truncated_window'})
            else:
                evidence.append({'sequence': core, 'offsets': [i for i, aa in enumerate(core) if aa == 'R'], 'full_length': 0, 'anchor_type': 'peptide'})
    return evidence


def build_anchors(row: pd.Series) -> list[dict]:
    """The row's evidence that can place a site. A window centred on another residue (for example a methyllysine
    window filed as an arginine site) gives no anchor: no arginine is guessed from it. A peptide names its modified
    arginine only when it has a single one; with several it is kept without indices, to verify the reported position."""
    anchors: list[dict] = []
    seen: set[tuple] = set()
    for item in row_evidence(row):
        seq, kind = item['sequence'], item['anchor_type']
        if item['full_length']:
            indices = []
        elif kind == 'peptide':
            indices = [item['offsets'][0] + 1] if len(item['offsets']) == 1 else []
        else:
            indices = sorted({i + 1 for i in item['offsets'] if 0 <= i < len(seq) and seq[i] == 'R'})
            if not indices:
                continue
        key = (seq, tuple(indices), item['full_length'], kind)
        if key not in seen:
            seen.add(key)
            anchors.append({'sequence': seq, 'indices': indices, 'full_length': item['full_length'], 'anchor_type': kind})
    return anchors


# A source sequence that is not in the entry still supports its reported arginine when, laid with its site on that
# arginine, it differs at no more than 2 residues and at no more than a fifth of them (variants, sequence conflicts).
# Residues that fall off the entry count as differences.
NEAR_MATCH_MAX_MISMATCHES = 2
NEAR_MATCH_MAX_FRACTION = 0.2
NEAR_MATCH_STATUS = {'same_accession_position_verified': 'same_accession_evidence_near_match',
                     'canonical_position_fallback': 'canonical_evidence_near_match'}


def near_match_mismatches(evidence: list[dict], sequence: str, position: int) -> int | None:
    """Fewest differences with which one of the row's sequences matches the entry with its site on position (1-based),
    or None when none matches within the limits above."""
    best = None
    for item in evidence:
        seq = item['sequence']
        if item['full_length']:
            half = item['full_length'] // 2
            offsets = [len(seq) - half - 1, half]
        else:
            offsets = item['offsets']
        for offset in offsets:
            if not 0 <= offset < len(seq) or seq[offset] != 'R':
                continue
            start = position - 1 - offset
            mismatches = sum(1 for i, aa in enumerate(seq) if not 0 <= start + i < len(sequence) or sequence[start + i] != aa)
            if mismatches <= NEAR_MATCH_MAX_MISMATCHES and mismatches <= NEAR_MATCH_MAX_FRACTION * len(seq):
                best = mismatches if best is None else min(best, mismatches)
    return best


def fallback_position(row: pd.Series, canonical_sequences: dict[str, str]) -> tuple[int | None, str, str]:
    """The reported position when it is an arginine on the row's own entry: (position, status, accession)."""
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


def load_peptide_coordinates(path: str | Path | None, canonical_sequences: dict[str, str]) -> dict[str, set[tuple[str, int, int]]]:
    """Per entry: (peptide, start in the source's numbering, start in the reference entry), from
    scripts/01_parse_maron_peptide_coordinates.py. A peptide counts when it lies at its stated start in the entry or
    occurs there once."""
    coordinates: dict[str, set[tuple[str, int, int]]] = defaultdict(set)
    if not path:
        return coordinates
    for r in pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False).itertuples():
        sequence, peptide, start = canonical_sequences.get(r.canonical_UniProtAC, ''), r.peptide, int(r.peptide_start)
        if not sequence or not peptide:
            continue
        if sequence[start - 1:start - 1 + len(peptide)] != peptide:
            found = find_all_occurrences(sequence, peptide)
            if len(found) != 1:
                continue
            coordinates[r.canonical_UniProtAC].add((peptide, start, found[0]))
        else:
            coordinates[r.canonical_UniProtAC].add((peptide, start, start))
    return coordinates


def peptide_start_position(anchors: list[dict], acc: str, reported: int, sequence: str,
                           coordinates: dict[str, set[tuple[str, int, int]]]) -> tuple[int | None, str]:
    """The arginine a plain peptide reports, counted from the peptide's start where a coordinate table gives it: the
    row's peptide equals, contains or lies inside a peptide with a known start on the same entry. Windows and annotated
    peptides mark their arginine themselves and are not used. Returns (position, peptide) when exactly one arginine
    results."""
    found: dict[int, str] = {}
    for anchor in anchors:
        peptide = anchor['sequence']
        if anchor['anchor_type'] != 'peptide':
            continue
        for known, source_start, reference_start in coordinates.get(acc, ()):
            if peptide in known:
                shift = known.index(peptide)
                source, reference = source_start + shift, reference_start + shift
            elif known in peptide:
                shift = peptide.index(known)
                source, reference = source_start - shift, reference_start - shift
            else:
                continue
            offset = reported - source
            if 0 <= offset < len(peptide) and peptide[offset] == 'R' and sequence[reference - 1:reference - 1 + len(peptide)] == peptide:
                found.setdefault(reference + offset, peptide)
    return next(iter(found.items())) if len(found) == 1 else (None, '')


def resolve_hit(row: pd.Series, hits: list[dict], canonical_acc: str) -> tuple[int | None, str, str, str, str, int]:
    """(position, status, accession, anchor type, anchor sequence, hit count) from the hits of a row's evidence."""
    if not hits:
        return None, 'unresolved_no_hit', '', '', '', 0

    def best_of(candidates: list[dict], longest: bool = False) -> dict:
        return sorted(candidates, key=lambda h: (h['priority'], -len(h['anchor_sequence']) if longest else len(h['anchor_sequence'])))[0]

    positions_by_acc: dict[str, set[int]] = defaultdict(set)
    for hit in hits:
        positions_by_acc[hit['matched_accession']].add(hit['matched_position'])

    if canonical_acc in positions_by_acc and len(positions_by_acc[canonical_acc]) == 1:
        position = sorted(positions_by_acc[canonical_acc])[0]
        best = best_of([hit for hit in hits if hit['matched_accession'] == canonical_acc and hit['matched_position'] == position])
        return position, 'sequence_unique_canonical', canonical_acc, best['anchor_type'], best['anchor_sequence'], len(hits)

    sites = {(hit['matched_accession'], hit['matched_position']) for hit in hits}
    if len(sites) == 1:
        acc, position = next(iter(sites))
        best = best_of(hits)
        return position, 'sequence_unique_candidate', acc, best['anchor_type'], best['anchor_sequence'], len(hits)

    # Several sites: take the one nearest the reported position. A tie is broken only by the row's own accession; two
    # other entries supported equally (for example paralogs behind an unreviewed accession) leave the row ambiguous.
    current_position = int(row.get('position', 0))
    distance = min(abs(position - current_position) for _, position in sites)
    nearest = {site for site in sites if abs(site[1] - current_position) == distance}
    own = {site for site in nearest if site[0] == canonical_acc}
    if len(own or nearest) == 1 and distance <= 25:
        acc, position = next(iter(own or nearest))
        best = best_of([hit for hit in hits if (hit['matched_accession'], hit['matched_position']) == (acc, position)], longest=True)
        return position, 'sequence_closest_supported', acc, best['anchor_type'], best['anchor_sequence'], len(hits)

    return None, 'unresolved_ambiguous_hit', '', '', '', len(hits)


def anchor_hits(anchor: dict, acc: str, start: int, sequence: str, reported: int) -> list[dict]:
    indices = anchor['indices']
    if anchor['full_length']:
        indices = truncated_window_indices(anchor['sequence'], anchor['full_length'], start, len(sequence), reported)
    return [{
        'matched_accession': acc,
        'matched_position': start + idx - 1,
        'anchor_type': anchor['anchor_type'],
        'anchor_sequence': anchor['sequence'],
        'priority': 0 if anchor['anchor_type'] == 'annotated_peptide' else 1,
    } for idx in indices]


def remap_row(row: pd.Series, canonical_sequences: dict[str, str], gene_map: dict[str, list[str]],
              secondary: dict[str, list[str]] | None = None, synonyms: dict[str, list[str]] | None = None,
              proteome: ProteomeIndex | None = None, coordinates: dict[str, set[tuple[str, int, int]]] | None = None) -> dict:
    source_acc = canonicalize_uniprot_accession(row.get('canonical_UniProtAC', row.get('substrate_UniProtAC', '')))
    # A merged accession's current entry counts as the row's own (ties, and peptides keeping their reported arginine).
    canonical_acc = own_accession(row, canonical_sequences, secondary or {}, gene_map)
    candidates = candidate_accessions(row, canonical_sequences, gene_map, secondary, synonyms)
    anchors = build_anchors(row)
    reported = int(row.get('position', 0))
    hits: list[dict] = []
    for anchor in anchors:
        for acc in candidates:
            sequence = canonical_sequences.get(acc, '')
            for start in find_all_occurrences(sequence, anchor['sequence']):
                hits.extend(anchor_hits(anchor, acc, start, sequence, reported))
    # No candidate at all (a deleted accession whose symbol names no current entry, such as a readthrough or an old clone
    # name): search the whole reference set and accept the evidence only if it names exactly one site there. A row with
    # candidates whose evidence lies elsewhere is left alone: placing it on a paralog is the open policy of issue 9.
    proteome_search = False
    if not hits and not candidates and proteome is not None and canonical_acc not in canonical_sequences:
        found = []
        for anchor in anchors:
            if len(anchor['sequence']) < PROTEOME_MIN_ANCHOR:
                continue
            for acc, start in proteome.find(anchor['sequence']):
                found.extend(anchor_hits(anchor, acc, start, canonical_sequences[acc], reported))
        if len({(h['matched_accession'], h['matched_position']) for h in found}) == 1:
            hits, proteome_search = found, True

    # The position belongs to the accession whose sequence matched. For a candidate found through the gene name (often
    # the reviewed entry behind an unreviewed or mislabelled source accession) that is not the row's own accession.
    corrected_position, remap_status, matched_accession, anchor_type, anchor_sequence, hit_count = resolve_hit(row, hits, canonical_acc)
    if proteome_search and remap_status == 'sequence_unique_candidate':
        remap_status = 'sequence_unique_proteome'

    # A peptide does not say which of its arginines was modified; the source position does. Keep the reported position
    # when it is an arginine inside an occurrence of one of the row's peptides on the canonical accession. Windows and
    # annotated peptides mark the modified residue themselves.
    canonical_seq = canonical_sequences.get(canonical_acc, '')
    if corrected_position != reported and canonical_seq and 1 <= reported <= len(canonical_seq) and canonical_seq[reported - 1] == 'R':
        for anchor in anchors:
            if anchor['anchor_type'] != 'peptide':
                continue
            length = len(anchor['sequence'])
            if any(start <= reported < start + length for start in find_all_occurrences(canonical_seq, anchor['sequence'])):
                corrected_position, remap_status, matched_accession = reported, 'sequence_position_verified', canonical_acc
                anchor_type, anchor_sequence = anchor['anchor_type'], anchor['sequence']
                hit_count = max(hit_count, 1)
                break
    if corrected_position is None:
        # The reported arginine stands on its own only for a row without sequence evidence. A row whose peptide or
        # window lies elsewhere in the entry is not placed; one whose evidence is not in the entry must match it near
        # the reported position (near_match_mismatches).
        position, fallback_status, fallback_acc = fallback_position(row, canonical_sequences)
        evidence = row_evidence(row)
        if position is not None and evidence:
            fallback_seq = canonical_sequences[fallback_acc]
            if any(find_all_occurrences(fallback_seq, item['sequence']) for item in evidence):
                position, fallback_status = None, 'unresolved_evidence_elsewhere'
            elif near_match_mismatches(evidence, fallback_seq, position) is None:
                position, fallback_status = None, 'unresolved_evidence_mismatch'
            else:
                fallback_status = NEAR_MATCH_STATUS[fallback_status]
        if position is not None:
            corrected_position, remap_status = position, fallback_status
            anchor_type = 'near_match' if evidence else 'position_only'
            anchor_sequence = ''
        elif remap_status == 'unresolved_no_hit':
            remap_status = fallback_status
    if corrected_position is None and coordinates:
        # A peptide with several arginines, reported in another numbering than the entry's: count the reported position
        # from the peptide's start, where the source's own coordinate table gives it.
        position, peptide = peptide_start_position(anchors, canonical_acc, reported, canonical_seq, coordinates)
        if position is not None:
            corrected_position, remap_status, matched_accession = position, 'sequence_peptide_start', canonical_acc
            anchor_type, anchor_sequence, hit_count = 'peptide', peptide, max(hit_count, 1)

    # A row without a source accession (Maron S5 rows whose accession reads NA) is placed only through its gene symbol's
    # candidates, so its accession is the matched one.
    corrected_accession = (matched_accession if remap_status.startswith('sequence_') else canonical_acc) if corrected_position is not None else ''
    n_current = len([a for a in (secondary or {}).get(source_acc, []) if a in canonical_sequences])
    history = ('split' if n_current > 1 else 'merged' if n_current == 1 else '') if source_acc not in canonical_sequences else ''
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
        'position_delta': (corrected_position - reported) if corrected_position is not None else pd.NA,
        'is_position_remapped': remapped and corrected_position != reported,
        'is_accession_remapped': remapped and corrected_accession != source_acc,
        'is_position_resolved': remapped,
        'source_accession_history': history,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--base-master', required=True)
    ap.add_argument('--maron', required=True)
    ap.add_argument('--prometheus', required=True)
    ap.add_argument('--maron-extra', default='', help='Table S5 rows read for the first time in this version (maron_s5_rows_accession_na.tsv); they come after all other rows')
    ap.add_argument('--study-records', nargs='*', default=[],
                    help='studies staged by scripts/00_stage_study_sites.py, read after the source tables above, in this order')
    ap.add_argument('--localised-only', action='store_true',
                    help='place only records whose source localises the methyl group to one residue (place_localised_records); v1.1.0 was built without it')
    ap.add_argument('--canonical-fasta', required=True)
    ap.add_argument('--uniprot-metadata', required=True)
    ap.add_argument('--accession-history', default='',
                    help='uniprot_accession_history.tsv from scripts/00: secondary accessions and gene synonyms of the reference set')
    ap.add_argument('--peptide-coordinates', default='',
                    help='maron_peptide_coordinates.tsv from scripts/01_parse_maron_peptide_coordinates.py: peptide starts in the source numbering')
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()

    canonical_sequences = parse_fasta(args.canonical_fasta)
    proteome = ProteomeIndex(canonical_sequences)
    coordinates = load_peptide_coordinates(args.peptide_coordinates, canonical_sequences)
    if args.peptide_coordinates:
        print(f'peptide coordinates: {sum(map(len, coordinates.values()))} peptides on {len(coordinates)} entries')
    _, gene_map = load_uniprot_metadata(args.uniprot_metadata)
    secondary, synonyms = load_accession_history(args.accession_history)
    if args.accession_history:
        print(f'accession history: {len(secondary)} secondary accessions, {len(synonyms)} gene synonyms')
    else:
        print('no --accession-history: merged or split accessions and outdated gene symbols find no candidate')

    base = pd.read_csv(args.base_master, sep='\t', low_memory=False)
    base = base[base['is_arg_methyl'].fillna(False)].copy()
    maron = pd.read_csv(args.maron, sep='\t', low_memory=False)
    prom = pd.read_csv(args.prometheus, sep='\t', low_memory=False)

    frames = [
        prepare_frame(base, 'human_ptm_master'),
        prepare_frame(maron, 'maron_s5'),
        prepare_frame(prom, 'prometheus_mmc2'),
    ]
    if args.maron_extra and Path(args.maron_extra).exists():
        frames.append(prepare_frame(pd.read_csv(args.maron_extra, sep='\t', low_memory=False, keep_default_na=False, na_values=['']), 'maron_s5'))
    frames += [study_frame(path) for path in args.study_records]
    all_rows = pd.concat(frames, ignore_index=True, sort=False)
    all_rows = restore_leucines(all_rows, canonical_sequences)
    all_rows = lead_reviewed_accession(all_rows, canonical_sequences, gene_map)
    all_rows['canonical_UniProtAC'] = all_rows['canonical_UniProtAC'].astype(str).map(canonicalize_uniprot_accession)
    all_rows = annotate_window_lengths(all_rows)
    kinds = all_rows.assign(reference=all_rows.get('reference', pd.Series('', index=all_rows.index)).fillna('')).groupby(
        ['source_table', 'reference'])[['seq_window_window_length', 'peptide_sequence_window_length']].first()
    print('window length per source (0 = peptides):')
    print(kinds.to_string())

    remap_annotations = all_rows.apply(
        lambda row: pd.Series(remap_row(row, canonical_sequences=canonical_sequences, gene_map=gene_map, secondary=secondary, synonyms=synonyms,
                                        proteome=proteome, coordinates=coordinates)),
        axis=1,
    )
    remapped = pd.concat([all_rows, remap_annotations], axis=1)
    if args.localised_only:
        remapped = place_localised_records(remapped, proteome)

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

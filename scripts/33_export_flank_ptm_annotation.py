#!/usr/bin/env python3
"""Export a compact non-methyl PTM annotation for an external flank table.

The methylarginine pipeline already assembles a cross-PTM site union with
enzyme attribution (``ptm_compare/ptm_sites_with_contexts.tsv``). Downstream
projects that work on a different coordinate universe -- domain-adjacent
flanks, in the MicroDomains repository -- need only the small slice of that
union that lands inside their windows, keyed by canonical accession and
position.

Exporting the slice rather than the whole union keeps the consumer
reproducible without a copy of this repository's staged inputs, and keeps the
provenance of every exported site in one place: the ``source`` and
``enzyme_genename`` columns come straight from the union and are not
re-derived.

Kinase families are labelled from the attributed enzyme, not from sequence.
CK2 is called out because it is the acidophilic kinase and is the one family a
consumer analysing acidic segments will want to separate from bulk
phosphorylation.

Usage:
    python scripts/33_export_flank_ptm_annotation.py \
        --flanks /path/to/continuous_microdomain_scores_public.tsv \
        --output PTM_results/exports/ptm_sites_in_domain_flanks.tsv
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd

from common import save_table


DEFAULT_SITES = 'PTM_results/ptm_compare/ptm_sites_with_contexts.tsv'
DEFAULT_OUTPUT = 'PTM_results/exports/ptm_sites_in_domain_flanks.tsv'

# Exported groups. Arg methylation is deliberately excluded: the consumer
# already carries its own methylarginine annotation and re-exporting it here
# would create two copies of the same evidence under different provenance.
EXPORT_GROUPS = ['Phosphorylation', 'Lys acylation', 'Lys methylation', 'Ubiquitin/SUMO']

# Catalytic and regulatory subunits of protein kinase CK2. The regulatory
# subunit CSNK2B has no catalytic activity of its own, but attribution in the
# source databases is to the holoenzyme, so a CSNK2B-attributed site is still
# a CK2 site.
CK2_GENES = {'CSNK2A1', 'CSNK2A2', 'CSNK2A3', 'CSNK2B'}

# Proline-directed kinases, kept as a labelled comparator: they are the other
# large family acting in disordered segments, and they are not acidophilic.
PRO_DIRECTED_GENES = {
    'CDK1', 'CDK2', 'CDK4', 'CDK5', 'CDK6', 'CDK7', 'CDK9',
    'MAPK1', 'MAPK3', 'MAPK8', 'MAPK9', 'MAPK14', 'GSK3A', 'GSK3B',
}

KEEP_COLUMNS = [
    'canonical_UniProtAC', 'position', 'residue', 'ptm_group', 'ptm_type',
    'substrate_genename', 'enzyme_genename', 'source', 'pmid',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--sites', default=os.environ.get('PTM_SITES_WITH_CONTEXTS', DEFAULT_SITES),
                        help='cross-PTM site table with contexts')
    parser.add_argument('--flanks', required=True,
                        help='flank table with canonical_UniProtAC, flank_start, flank_end')
    parser.add_argument('--output', default=DEFAULT_OUTPUT, help='destination TSV')
    parser.add_argument('--groups', nargs='*', default=EXPORT_GROUPS,
                        help='ptm_group values to export')
    return parser.parse_args()


def load_sites(path: Path, groups: list[str]) -> pd.DataFrame:
    sites = pd.read_csv(path, sep='\t', low_memory=False)
    missing = [c for c in KEEP_COLUMNS if c not in sites.columns]
    if missing:
        raise SystemExit(f'site table is missing required columns: {missing}')
    sites = sites.loc[sites['ptm_group'].isin(groups), KEEP_COLUMNS].copy()
    sites['position'] = pd.to_numeric(sites['position'], errors='coerce')
    sites = sites.dropna(subset=['canonical_UniProtAC', 'position'])
    sites['position'] = sites['position'].astype(int)
    return sites


def load_windows(path: Path) -> pd.DataFrame:
    flanks = pd.read_csv(path, sep='\t', low_memory=False,
                         usecols=['canonical_UniProtAC', 'flank_start', 'flank_end'])
    flanks = flanks.dropna()
    flanks['flank_start'] = flanks['flank_start'].astype(int)
    flanks['flank_end'] = flanks['flank_end'].astype(int)
    return flanks.drop_duplicates()


def in_any_window(sites: pd.DataFrame, windows: pd.DataFrame) -> np.ndarray:
    """True where a site falls inside at least one flank of the same protein.

    Grouped by accession rather than merged, because a protein contributes many
    overlapping windows and a merge would multiply each site by the number of
    windows containing it.
    """
    by_ac: dict[str, tuple[np.ndarray, np.ndarray]] = {
        ac: (sub['flank_start'].to_numpy(), sub['flank_end'].to_numpy())
        for ac, sub in windows.groupby('canonical_UniProtAC')
    }
    keep = np.zeros(len(sites), dtype=bool)
    for i, (ac, pos) in enumerate(zip(sites['canonical_UniProtAC'], sites['position'])):
        bounds = by_ac.get(ac)
        if bounds is None:
            continue
        starts, ends = bounds
        keep[i] = bool(np.any((pos >= starts) & (pos <= ends)))
    return keep


def main() -> None:
    args = parse_args()
    sites = load_sites(Path(args.sites), list(args.groups))
    windows = load_windows(Path(args.flanks))

    kept = sites.loc[in_any_window(sites, windows)].copy()
    enzyme = kept['enzyme_genename'].fillna('')
    kept['is_ck2'] = enzyme.isin(CK2_GENES)
    kept['is_proline_directed'] = enzyme.isin(PRO_DIRECTED_GENES)
    kept['has_enzyme_attribution'] = enzyme.ne('')
    kept = kept.sort_values(['canonical_UniProtAC', 'position', 'ptm_group'])
    kept = kept.drop_duplicates(
        subset=['canonical_UniProtAC', 'position', 'ptm_group', 'enzyme_genename'])

    save_table(kept, args.output)

    counts = kept.groupby('ptm_group').size().sort_values(ascending=False)
    print(f'sites read:        {len(sites)}')
    print(f'flank windows:     {len(windows)} over {windows["canonical_UniProtAC"].nunique()} proteins')
    print(f'sites in a flank:  {len(kept)} over {kept["canonical_UniProtAC"].nunique()} proteins')
    print(f'  CK2-attributed:  {int(kept["is_ck2"].sum())}')
    print(f'  Pro-directed:    {int(kept["is_proline_directed"].sum())}')
    print(f'  any attribution: {int(kept["has_enzyme_attribution"].sum())}')
    print('by group:')
    for group, n in counts.items():
        print(f'  {group:<18} {n}')
    print(f'written: {args.output}')


if __name__ == '__main__':
    main()

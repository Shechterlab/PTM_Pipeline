from __future__ import annotations

import argparse
import math
from collections import defaultdict
from pathlib import Path

import pandas as pd
from scipy.stats import binomtest


BIN_ORDER = ["Domain edge <=20 aa", "Domain edge 21-40 aa", "Domain-distal >40 aa"]


def compact_bin(distance: float | int | None) -> str | pd.NA:
    if pd.isna(distance):
        return pd.NA
    distance = float(distance)
    if distance <= 20:
        return BIN_ORDER[0]
    if distance <= 40:
        return BIN_ORDER[1]
    return BIN_ORDER[2]


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    ordered = sorted((int(start), int(end)) for start, end in intervals if pd.notna(start) and pd.notna(end) and start <= end)
    if not ordered:
        return []
    merged = [list(ordered[0])]
    for start, end in ordered[1:]:
        if start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def subtract_intervals(base: list[tuple[int, int]], blockers: list[tuple[int, int]]) -> list[tuple[int, int]]:
    blockers = merge_intervals(blockers)
    out: list[tuple[int, int]] = []
    for start, end in merge_intervals(base):
        pieces = [(start, end)]
        for block_start, block_end in blockers:
            next_pieces: list[tuple[int, int]] = []
            for piece_start, piece_end in pieces:
                if block_end < piece_start or block_start > piece_end:
                    next_pieces.append((piece_start, piece_end))
                    continue
                if piece_start < block_start:
                    next_pieces.append((piece_start, block_start - 1))
                if block_end < piece_end:
                    next_pieces.append((block_end + 1, piece_end))
            pieces = next_pieces
            if not pieces:
                break
        out.extend(pieces)
    return merge_intervals(out)


def bh_qvalues(p_values: list[float]) -> list[float]:
    n = len(p_values)
    indexed = sorted(enumerate(p_values), key=lambda item: item[1])
    q = [math.nan] * n
    running = 1.0
    for rank, (idx, p_value) in reversed(list(enumerate(indexed, start=1))):
        running = min(running, p_value * n / rank)
        q[idx] = running
    return q


def load_segments(disorder_path: Path, interpro_path: Path, accessions: set[str]) -> pd.DataFrame:
    disorder = pd.read_csv(disorder_path, sep="\t", low_memory=False)
    domains = pd.read_csv(interpro_path, sep="\t", low_memory=False)
    disorder = disorder[disorder["canonical_UniProtAC"].astype(str).isin(accessions)].copy()
    domains = domains[domains["canonical_UniProtAC"].astype(str).isin(accessions)].copy()
    domains = domains[domains["interpro_type"].astype(str).isin({"domain", "repeat"})].copy()

    disorder_by_acc = defaultdict(list)
    for row in disorder.itertuples(index=False):
        disorder_by_acc[str(row.canonical_UniProtAC)].append((int(row.fragment_start), int(row.fragment_end)))

    domain_by_acc = defaultdict(list)
    for row in domains.itertuples(index=False):
        domain_by_acc[str(row.canonical_UniProtAC)].append((int(row.fragment_start), int(row.fragment_end)))

    rows = []
    for accession in sorted(accessions):
        outside_domain_idrs = subtract_intervals(disorder_by_acc.get(accession, []), domain_by_acc.get(accession, []))
        for idx, (start, end) in enumerate(outside_domain_idrs, start=1):
            rows.append(
                {
                    "canonical_UniProtAC": accession,
                    "outside_domain_idr_segment_id": f"{accession}:{idx}:{start}-{end}",
                    "outside_domain_idr_start": start,
                    "outside_domain_idr_end": end,
                    "outside_domain_idr_length": end - start + 1,
                }
            )
    return pd.DataFrame(rows)


def assign_segments(rows: pd.DataFrame, segments: pd.DataFrame) -> pd.DataFrame:
    out = rows.copy()
    for col in ["outside_domain_idr_segment_id", "outside_domain_idr_start", "outside_domain_idr_end", "outside_domain_idr_length"]:
        out[col] = pd.NA

    for accession, group in out.groupby("canonical_UniProtAC", sort=False):
        segs = segments[segments["canonical_UniProtAC"] == accession]
        if segs.empty:
            continue
        positions = pd.to_numeric(group["position"], errors="coerce")
        for seg in segs.itertuples(index=False):
            mask = positions.between(int(seg.outside_domain_idr_start), int(seg.outside_domain_idr_end), inclusive="both")
            idx = group.index[mask]
            if idx.empty:
                continue
            out.loc[idx, "outside_domain_idr_segment_id"] = seg.outside_domain_idr_segment_id
            out.loc[idx, "outside_domain_idr_start"] = int(seg.outside_domain_idr_start)
            out.loc[idx, "outside_domain_idr_end"] = int(seg.outside_domain_idr_end)
            out.loc[idx, "outside_domain_idr_length"] = int(seg.outside_domain_idr_length)
    return out


def build_length_conditioned_table(observed: pd.DataFrame, arginine_universe: pd.DataFrame) -> pd.DataFrame:
    segment_bin_counts = (
        arginine_universe.groupby(["outside_domain_idr_segment_id", "domain_edge_bin_compact"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    for category in BIN_ORDER:
        if category not in segment_bin_counts.columns:
            segment_bin_counts[category] = 0
    segment_bin_counts = segment_bin_counts[BIN_ORDER]
    segment_totals = segment_bin_counts.sum(axis=1)
    segment_bin_fractions = segment_bin_counts.div(segment_totals, axis=0)

    observed_segments = observed["outside_domain_idr_segment_id"].dropna().astype(str)
    expected_counts = segment_bin_fractions.loc[observed_segments].sum(axis=0)
    observed_counts = observed["domain_edge_bin_compact"].value_counts().reindex(BIN_ORDER, fill_value=0)
    total_observed = int(observed_counts.sum())

    rows = []
    p_values = []
    for category in BIN_ORDER:
        expected = float(expected_counts[category])
        expected_fraction = expected / total_observed if total_observed else math.nan
        observed_count = int(observed_counts[category])
        observed_fraction = observed_count / total_observed if total_observed else math.nan
        p_value = (
            float(binomtest(observed_count, total_observed, expected_fraction, alternative="two-sided").pvalue)
            if total_observed and expected_fraction > 0
            else math.nan
        )
        p_values.append(p_value if not math.isnan(p_value) else 1.0)
        ratio = observed_fraction / expected_fraction if expected_fraction else math.nan
        rows.append(
            {
                "category": category,
                "observed_count": observed_count,
                "observed_fraction": observed_fraction,
                "expected_count_length_conditioned": expected,
                "expected_fraction_length_conditioned": expected_fraction,
                "enrichment_ratio_vs_length_null": ratio,
                "log2_enrichment_vs_length_null": math.log2(ratio) if ratio and ratio > 0 else math.nan,
                "binomial_p_value": p_value,
            }
        )
    q_values = bh_qvalues(p_values)
    for row, q_value in zip(rows, q_values):
        row["q_value"] = q_value
    return pd.DataFrame(rows)


def add_length_bin(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["outside_domain_idr_length_bin"] = pd.cut(
        pd.to_numeric(out["outside_domain_idr_length"], errors="coerce"),
        bins=[0, 40, 80, 160, 10**9],
        labels=["<=40 aa", "41-80 aa", "81-160 aa", ">160 aa"],
        include_lowest=True,
    )
    return out


def build_length_stratum_table(observed: pd.DataFrame, arginine_universe: pd.DataFrame) -> pd.DataFrame:
    rows = []
    observed = add_length_bin(observed)
    arginine_universe = add_length_bin(arginine_universe)
    for length_bin in ["<=40 aa", "41-80 aa", "81-160 aa", ">160 aa"]:
        obs_subset = observed[observed["outside_domain_idr_length_bin"].astype(str) == length_bin].copy()
        if obs_subset.empty:
            continue
        segment_ids = set(obs_subset["outside_domain_idr_segment_id"].dropna().astype(str))
        arg_subset = arginine_universe[arginine_universe["outside_domain_idr_segment_id"].astype(str).isin(segment_ids)].copy()
        table = build_length_conditioned_table(obs_subset, arg_subset)
        table.insert(0, "outside_domain_idr_length_bin", length_bin)
        table.insert(1, "observed_sites_in_length_bin", int(len(obs_subset)))
        rows.append(table)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def main() -> None:
    parser = argparse.ArgumentParser(description="Recalculate domain-edge bins with an IDR-segment-length-conditioned null.")
    parser.add_argument("--sites", default="PTM_results/domain_context/sites_with_domain_context.tsv")
    parser.add_argument("--background", default="PTM_results/domain_context/same_protein_nonmethyl_arginine_domain_background.tsv")
    parser.add_argument("--disorder", default="PTM_data/context/mobidb_human_reviewed_disorder_intervals.tsv")
    parser.add_argument("--interpro", default="PTM_data/context/interpro_human_reviewed_domain_like_intervals.tsv")
    parser.add_argument("--outdir", default="PTM_results/domain_context")
    args = parser.parse_args()

    sites = pd.read_csv(args.sites, sep="\t", low_memory=False)
    background = pd.read_csv(args.background, sep="\t", low_memory=False)

    site_rows = sites.copy()
    if "corrected_position" in site_rows.columns:
        site_rows["position"] = site_rows["corrected_position"]
    site_rows["is_methyl_site"] = True
    background["is_methyl_site"] = False
    common_cols = [
        "canonical_UniProtAC",
        "position",
        "site_key",
        "domain_context_class",
        "nearest_domain_edge_distance",
        "domain_edge_bin_compact",
        "disorder_context_class",
        "is_methyl_site",
    ]
    site_rows = site_rows[common_cols]
    background = background[common_cols]
    arginine_universe = pd.concat([site_rows, background], ignore_index=True)

    methylated_accessions = set(sites["canonical_UniProtAC"].dropna().astype(str))
    segments = load_segments(Path(args.disorder), Path(args.interpro), methylated_accessions)
    arginine_universe = assign_segments(arginine_universe, segments)

    outside_disordered_mask = (
        arginine_universe["outside_domain_idr_segment_id"].notna()
        & arginine_universe["nearest_domain_edge_distance"].notna()
        & arginine_universe["domain_context_class"].isin(["boundary", "inter_domain_linker", "distal"])
    )
    arginine_universe = arginine_universe[outside_disordered_mask].copy()
    arginine_universe["domain_edge_bin_compact"] = arginine_universe["nearest_domain_edge_distance"].map(compact_bin)

    observed = arginine_universe[arginine_universe["is_methyl_site"]].copy()
    length_null = build_length_conditioned_table(observed, arginine_universe)
    length_strata = build_length_stratum_table(observed, arginine_universe)

    observed_segment_ids = set(observed["outside_domain_idr_segment_id"].dropna().astype(str))
    observed_segments = segments[segments["outside_domain_idr_segment_id"].isin(observed_segment_ids)].copy()
    segment_arg_counts = (
        arginine_universe.groupby("outside_domain_idr_segment_id").size().rename("arginines_in_segment").reset_index()
    )
    observed_segments = observed_segments.merge(segment_arg_counts, on="outside_domain_idr_segment_id", how="left")
    observed_annotated = observed[
        [
            "canonical_UniProtAC",
            "position",
            "site_key",
            "nearest_domain_edge_distance",
            "domain_edge_bin_compact",
            "outside_domain_idr_segment_id",
            "outside_domain_idr_start",
            "outside_domain_idr_end",
            "outside_domain_idr_length",
        ]
    ].sort_values(["canonical_UniProtAC", "position"])

    qc = pd.DataFrame(
        [
            {
                "methylated_disordered_outside_domain_sites": int(len(observed)),
                "arginines_in_same_outside_domain_idr_segments": int(len(arginine_universe)),
                "outside_domain_idr_segments_in_methylated_proteins": int(len(segments)),
                "outside_domain_idr_segments_with_methyl_site": int(len(observed_segment_ids)),
                "median_segment_length_all_segments": float(segments["outside_domain_idr_length"].median()) if not segments.empty else math.nan,
                "median_segment_length_all_arginines_arg_weighted": float(arginine_universe["outside_domain_idr_length"].median()) if not arginine_universe.empty else math.nan,
                "median_segment_length_methyl_sites_site_weighted": float(observed["outside_domain_idr_length"].median()) if not observed.empty else math.nan,
                "median_domain_edge_distance_methyl_sites": float(observed["nearest_domain_edge_distance"].median()) if not observed.empty else math.nan,
            }
        ]
    )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    length_null.to_csv(outdir / "domain_edge_bin_enrichment_disordered_only_compact_idr_length_null.tsv", sep="\t", index=False)
    length_strata.to_csv(outdir / "domain_edge_bin_enrichment_disordered_only_compact_idr_length_null_by_segment_length.tsv", sep="\t", index=False)
    qc.to_csv(outdir / "domain_edge_idr_length_null_qc.tsv", sep="\t", index=False)
    observed_segments.sort_values(["canonical_UniProtAC", "outside_domain_idr_start"]).to_csv(
        outdir / "domain_edge_idr_length_null_segments_with_methyl_sites.tsv", sep="\t", index=False
    )
    observed_annotated.to_csv(outdir / "domain_edge_idr_length_null_sites.tsv", sep="\t", index=False)

    print(outdir / "domain_edge_bin_enrichment_disordered_only_compact_idr_length_null.tsv")
    print(outdir / "domain_edge_idr_length_null_qc.tsv")


if __name__ == "__main__":
    main()

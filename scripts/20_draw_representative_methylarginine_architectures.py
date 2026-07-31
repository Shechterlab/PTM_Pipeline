from __future__ import annotations

import argparse
from dataclasses import dataclass
from html import escape
from pathlib import Path

import pandas as pd

from common import parse_fasta


FIGURE_WIDTH = 1600
FIGURE_HEIGHT = 980
ROW_CARD_X = 42
ROW_CARD_W = 1516
ROW_CARD_H = 154
TRACK_X0 = 370
TRACK_X1 = 1098
INTERACTOR_X = 1140

CARD_WIDTH = 1520
CARD_HEIGHT = 230
SINGLE_ROW_Y = 54

COLORS = {
    "bg": "#FFFFFF",
    "ink": "#1A2433",
    "muted": "#637083",
    "card": "#FCFDFE",
    "card_edge": "#E3E8F0",
    "divider": "#ECF0F4",
    "domain_fill": "#E9DFC9",
    "domain_edge": "#736858",
    "idr": "#8C96A3",
    "rg_fill": "#F5D6DE",
    "prld_fill": "#F5E9C8",
    "tail_fill": "#E6EEF7",
    "site_fill": "#CBB2B9",
    "site_edge": "#8E747C",
    "site_halo": "#E8DDE1",
    "reader_fill": "#DCEBFA",
    "reader_edge": "#7B97B6",
    "target_fill": "#E0ECFA",
    "target_edge": "#4F79A8",
    "tag_fill": "#EEF4FA",
    "tag_edge": "#CAD7E7",
    "chip_fill": "#F5F8FB",
    "chip_edge": "#D9E1EA",
    "shadow": "#9FB0C8",
    "adj_fill_near": "#EAF2FB",
    "adj_fill_mid": "#F4F8FC",
}


@dataclass(frozen=True)
class Feature:
    start: int
    end: int
    label: str
    kind: str


@dataclass(frozen=True)
class TrackRow:
    panel_id: str
    display_name: str
    descriptor: str
    accession: str
    length: int
    domains: list[Feature]
    regions: list[Feature]
    site_positions: list[int]
    tag: str
    tag_detail: str
    interactor_kind: str
    interactor_label: str
    peptide_inset: str = ""


def fmt(value: float | int) -> str:
    if isinstance(value, int):
        return str(value)
    if abs(value - round(value)) < 1e-6:
        return str(int(round(value)))
    return f"{value:.2f}".rstrip("0").rstrip(".")


def attrs(**kwargs: object) -> str:
    out: list[str] = []
    for key, value in kwargs.items():
        if value is None:
            continue
        out.append(f'{key.rstrip("_").replace("_", "-")}="{escape(str(value), quote=True)}"')
    return " ".join(out)


def rect(x: float, y: float, width: float, height: float, **kwargs: object) -> str:
    return f'<rect {attrs(x=fmt(x), y=fmt(y), width=fmt(width), height=fmt(height), **kwargs)} />'


def line(x1: float, y1: float, x2: float, y2: float, **kwargs: object) -> str:
    return f'<line {attrs(x1=fmt(x1), y1=fmt(y1), x2=fmt(x2), y2=fmt(y2), **kwargs)} />'


def circle(cx: float, cy: float, r: float, **kwargs: object) -> str:
    return f'<circle {attrs(cx=fmt(cx), cy=fmt(cy), r=fmt(r), **kwargs)} />'


def path(d: str, **kwargs: object) -> str:
    return f'<path {attrs(d=d, **kwargs)} />'


def text(
    x: float,
    y: float,
    content: str,
    font_size: float,
    *,
    anchor: str = "start",
    weight: str = "400",
    fill: str = COLORS["ink"],
    italic: bool = False,
) -> str:
    text_attrs = attrs(
        x=fmt(x),
        y=fmt(y),
        fill=fill,
        font_size=fmt(font_size),
        text_anchor=anchor,
        font_weight=weight,
        font_style=("italic" if italic else "normal"),
    )
    return f'<text {text_attrs}>{escape(content)}</text>'


def multiline_text(
    x: float,
    y: float,
    lines: list[str],
    font_size: float,
    *,
    anchor: str = "start",
    fill: str = COLORS["ink"],
    weight: str = "400",
    line_height: float | None = None,
) -> str:
    line_height = line_height or font_size * 1.24
    body = [f'<text {attrs(x=fmt(x), y=fmt(y), fill=fill, font_size=fmt(font_size), text_anchor=anchor, font_weight=weight)}>' ]
    for idx, item in enumerate(lines):
        dy = "0" if idx == 0 else fmt(line_height)
        body.append(f'<tspan x="{fmt(x)}" dy="{dy}">{escape(item)}</tspan>')
    body.append("</text>")
    return "".join(body)


def wrap_text(content: str, width_px: float, font_size: float) -> list[str]:
    max_chars = max(18, int(width_px / max(font_size * 0.56, 1.0)))
    words = content.split()
    if not words:
        return [""]
    lines: list[str] = []
    current = words[0]
    for word in words[1:]:
        candidate = f"{current} {word}"
        if len(candidate) <= max_chars:
            current = candidate
        else:
            lines.append(current)
            current = word
    lines.append(current)
    return lines


def merge_intervals(intervals: list[tuple[int, int]], max_gap: int = 0) -> list[tuple[int, int]]:
    ordered = sorted((int(start), int(end)) for start, end in intervals if start <= end)
    if not ordered:
        return []
    merged = [list(ordered[0])]
    for start, end in ordered[1:]:
        if start <= merged[-1][1] + max_gap + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def intersect_intervals(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for start_a, end_a in a:
        for start_b, end_b in b:
            start = max(start_a, start_b)
            end = min(end_a, end_b)
            if start <= end:
                out.append((start, end))
    return merge_intervals(out, max_gap=0)


def rg_rich_windows(sequence: str, *, window: int = 21, rg_fraction: float = 0.45, min_motif_hits: int = 2) -> list[tuple[int, int]]:
    hits: list[tuple[int, int]] = []
    if len(sequence) < window:
        return hits
    for idx in range(len(sequence) - window + 1):
        segment = sequence[idx : idx + window]
        fraction = sum(aa in {"R", "G"} for aa in segment) / window
        motif_hits = sum(segment[pos : pos + 2] in {"RG", "GR"} for pos in range(len(segment) - 1)) + segment.count("RGG")
        if fraction >= rg_fraction and motif_hits >= min_motif_hits:
            hits.append((idx + 1, idx + window))
    return merge_intervals(hits, max_gap=5)


def pos_to_x(position: int, scale_length: int, x0: float = TRACK_X0, x1: float = TRACK_X1) -> float:
    if scale_length <= 1:
        return x0
    return x0 + (x1 - x0) * ((position - 1) / (scale_length - 1))


def wavy_segment(x0: float, x1: float, y: float, amplitude: float = 4.0, period: float = 26.0) -> str:
    if x1 <= x0:
        return ""
    steps = max(1, int((x1 - x0) / period))
    step = (x1 - x0) / steps
    parts = [f"M {fmt(x0)} {fmt(y)}"]
    cursor = x0
    for idx in range(steps):
        nxt = x1 if idx == steps - 1 else cursor + step
        mid = (cursor + nxt) / 2
        peak = y - amplitude if idx % 2 == 0 else y + amplitude
        parts.append(f"Q {fmt(mid)} {fmt(peak)} {fmt(nxt)} {fmt(y)}")
        cursor = nxt
    return " ".join(parts)


def defs_block() -> str:
    return """
<defs>
  <filter id="softShadow" x="-20%" y="-20%" width="140%" height="160%">
    <feDropShadow dx="0" dy="4" stdDeviation="5" flood-color="#9FB0C8" flood-opacity="0.18" />
  </filter>
</defs>
""".strip()


def style_block() -> str:
    return """
<style>
  svg { background: #FFFFFF; }
  text {
    font-family: Arial, "Liberation Sans", "DejaVu Sans", Helvetica, sans-serif;
    letter-spacing: 0.01em;
  }
</style>
""".strip()


def load_disorder_regions(path: Path) -> dict[str, list[tuple[int, int]]]:
    frame = pd.read_csv(path, sep="\t", low_memory=False)
    out: dict[str, list[tuple[int, int]]] = {}
    for accession, group in frame.groupby("canonical_UniProtAC", sort=False):
        pairs = [(int(row.fragment_start), int(row.fragment_end)) for row in group.itertuples(index=False)]
        out[str(accession)] = merge_intervals(pairs)
    return out


def domain_union(frame: pd.DataFrame, accession: str, query: str, *, interpro_type: str = "domain") -> tuple[int, int]:
    subset = frame[
        (frame["canonical_UniProtAC"] == accession)
        & (frame["interpro_type"] == interpro_type)
        & (frame["interpro_name"].astype(str).str.contains(query, case=False, na=False, regex=False))
    ]
    if subset.empty:
        raise ValueError(f"No interval for {accession}: {query}")
    return int(subset["fragment_start"].min()), int(subset["fragment_end"].max())


def build_rows(
    site_table: pd.DataFrame,
    remapped_rows: pd.DataFrame,
    metadata: pd.DataFrame,
    interpro: pd.DataFrame,
    disorder_by_acc: dict[str, list[tuple[int, int]]],
    sequences: dict[str, str],
) -> list[TrackRow]:
    lengths = {
        str(row["Entry"]): int(row["Length"])
        for _, row in metadata[["Entry", "Length"]].dropna().drop_duplicates(subset=["Entry"]).iterrows()
    }

    def positions_for(accession: str, gene: str, *, start: int | None = None, end: int | None = None, consensus: str | None = None) -> list[int]:
        subset = site_table[
            (site_table["canonical_UniProtAC"] == accession)
            & (site_table["substrate_genename"] == gene)
        ].copy()
        if start is not None:
            subset = subset[subset["position"] >= start]
        if end is not None:
            subset = subset[subset["position"] <= end]
        if consensus is not None:
            subset = subset[subset["arg_methyl_state_consensus"] == consensus]
        subset = subset.dropna(subset=["position"]).drop_duplicates(subset=["position"]).sort_values("position")
        return [int(value) for value in subset["position"].tolist()]

    fus_acc = "P35637"
    fus_rg = intersect_intervals(rg_rich_windows(sequences[fus_acc]), disorder_by_acc[fus_acc])
    fus_prld_end = fus_rg[0][0] - 1
    fus_sites = positions_for(fus_acc, "FUS")
    fus_py_start = max(1, sequences[fus_acc].rfind("PY") + 1 - 12)

    snrpb_acc = "P14678"
    snrpb_tail = (90, lengths[snrpb_acc])
    snrpb_subset = site_table[
        (site_table["canonical_UniProtAC"] == snrpb_acc)
        & (site_table["substrate_genename"] == "SNRPB")
        & (site_table["position"].between(snrpb_tail[0], snrpb_tail[1]))
        & (site_table["arg_methyl_state_consensus"].fillna("").astype(str) != "")
    ].copy()
    snrpb_subset = snrpb_subset.drop_duplicates(subset=["position"]).sort_values("position")
    snrpb_sites = [int(value) for value in snrpb_subset["position"].tolist()]

    ews_acc = "Q01844"
    ews_rg = intersect_intervals(rg_rich_windows(sequences[ews_acc]), disorder_by_acc[ews_acc])
    ews_sites = positions_for(ews_acc, "EWSR1")

    h3_acc = "P68431"
    h3_site_rows = remapped_rows[
        (remapped_rows["canonical_UniProtAC"] == h3_acc)
        & (remapped_rows["substrate_genename"] == "H3C1")
        & (remapped_rows["site"] == "R2")
    ].copy()
    h3_site = 3
    if not h3_site_rows.empty and h3_site_rows["corrected_position"].notna().any():
        h3_site = int(h3_site_rows["corrected_position"].mode().iloc[0])

    return [
        TrackRow(
            panel_id="fus",
            display_name="FUS",
            descriptor="clustered RGG/RG-rich IDRs",
            accession=fus_acc,
            length=lengths[fus_acc],
            domains=[
                Feature(*domain_union(interpro, fus_acc, "RNA recognition motif domain"), "RRM", "domain"),
                Feature(*domain_union(interpro, fus_acc, "Zinc finger, RanBP2-type"), "ZnF", "domain"),
            ],
            regions=[
                Feature(1, fus_prld_end, "PrLD / low complexity", "prld"),
                Feature(fus_rg[0][0], fus_rg[-1][1], "RGG/RG-rich IDRs", "rg_span"),
                Feature(fus_rg[0][0], fus_rg[0][1], "", "rg"),
                Feature(fus_rg[1][0], fus_rg[1][1], "", "rg"),
                Feature(fus_rg[2][0], fus_rg[2][1], "", "rg"),
                Feature(fus_py_start, lengths[fus_acc], "PY-NLS", "tail"),
            ],
            site_positions=fus_sites,
            tag="multivalent IDR network tuning",
            tag_detail="LLPS / transport",
            interactor_kind="fus_network",
            interactor_label="LLPS / transport",
        ),
        TrackRow(
            panel_id="snrpb",
            display_name="SmB/B′ / SNRPB",
            descriptor="methyl-dependent SMN-binding C-terminal RG tail",
            accession=snrpb_acc,
            length=lengths[snrpb_acc],
            domains=[Feature(*domain_union(interpro, snrpb_acc, "Sm domain"), "Sm core", "domain")],
            regions=[Feature(snrpb_tail[0], snrpb_tail[1], "C-terminal RG tail", "tail")],
            site_positions=snrpb_sites,
            tag="reader-routed RNP assembly",
            tag_detail="methylated RG tail -> SMN Tudor",
            interactor_kind="smn",
            interactor_label="SMN Tudor",
        ),
        TrackRow(
            panel_id="ewsr1",
            display_name="EWS / EWSR1",
            descriptor="dense disordered RGG modules adjacent to RNA-binding domains",
            accession=ews_acc,
            length=lengths[ews_acc],
            domains=[
                Feature(*domain_union(interpro, ews_acc, "RNA recognition motif domain"), "RRM", "domain"),
                Feature(*domain_union(interpro, ews_acc, "Zinc finger, RanBP2-type"), "ZnF", "domain"),
            ],
            regions=[
                Feature(ews_rg[0][0], ews_rg[-1][1], "RGG/RG-rich modules", "rg_span"),
                Feature(ews_rg[0][0], ews_rg[0][1], "", "rg"),
                Feature(ews_rg[1][0], ews_rg[1][1], "", "rg"),
                Feature(ews_rg[2][0], ews_rg[2][1], "", "rg"),
            ],
            site_positions=ews_sites,
            tag="direct G4 DNA / RNA contact tuning",
            tag_detail="EWS RGG3 -> G4 DNA/RNA",
            interactor_kind="g4",
            interactor_label="G4 DNA/RNA",
        ),
        TrackRow(
            panel_id="h3",
            display_name="H3.1 / H3C1",
            descriptor="single chromatin reader-gating site",
            accession=h3_acc,
            length=lengths[h3_acc],
            domains=[Feature(*domain_union(interpro, h3_acc, "Core Histone H2A/H2B/H3 domain"), "histone fold", "domain")],
            regions=[Feature(1, 43, "H3 N-tail", "tail")],
            site_positions=[h3_site],
            tag="single-site chromatin reader gating",
            tag_detail="H3R2 -> WDR5 pocket",
            interactor_kind="wdr5",
            interactor_label="WDR5 WD-repeat pocket",
            peptide_inset="ARTKQTARK",
        ),
    ]


def draw_chip(elements: list[str], x: float, y: float, label: str) -> None:
    width = max(158.0, 24.0 + len(label) * 7.0)
    elements.append(rect(x, y, width, 32, rx=16, fill=COLORS["chip_fill"], stroke=COLORS["chip_edge"], stroke_width="1"))
    elements.append(text(x + width / 2, y + 21, label, 11.4, anchor="middle", fill=COLORS["muted"], weight="600"))


def draw_track_regions(elements: list[str], row: TrackRow, track_y: float, scale_length: int, track_end_x: float) -> None:
    for feature in row.regions:
        x0 = pos_to_x(feature.start, scale_length)
        x1 = pos_to_x(feature.end, scale_length)
        width = max(8.0, x1 - x0)
        if feature.kind == "prld":
            elements.append(rect(x0, track_y - 10, width, 20, rx=10, fill=COLORS["prld_fill"], stroke="none"))
            elements.append(text((x0 + x1) / 2, track_y - 16, feature.label, 11.2, anchor="middle", fill=COLORS["muted"], weight="600"))
        elif feature.kind == "rg_span":
            elements.append(text((x0 + x1) / 2, track_y - 16, feature.label, 11.3, anchor="middle", fill=COLORS["muted"], weight="600"))
        elif feature.kind == "rg":
            elements.append(rect(x0, track_y - 9, width, 18, rx=9, fill=COLORS["rg_fill"], stroke="none"))
        elif feature.kind == "tail":
            elements.append(rect(x0, track_y - 10, width, 20, rx=10, fill=COLORS["tail_fill"], stroke="none"))
            label_y = track_y - 16 if width >= 72 else track_y + 34
            label_anchor = "middle" if width >= 72 else "start"
            label_x = (x0 + x1) / 2 if width >= 72 else min(track_end_x - 10, x1 + 10)
            elements.append(text(label_x, label_y, feature.label, 11.0, anchor=label_anchor, fill=COLORS["muted"], weight="600"))


def draw_domains(elements: list[str], row: TrackRow, track_y: float, scale_length: int) -> None:
    for feature in row.domains:
        x0 = pos_to_x(feature.start, scale_length)
        x1 = pos_to_x(feature.end, scale_length)
        width = max(12.0, x1 - x0)
        elements.append(
            rect(
                x0,
                track_y - 14,
                width,
                28,
                rx=14,
                fill=COLORS["domain_fill"],
                stroke=COLORS["domain_edge"],
                stroke_width="1.25",
            )
        )
        font = 10.5 if width >= 68 else 9.6
        elements.append(text((x0 + x1) / 2, track_y + 4, feature.label, font, anchor="middle", fill=COLORS["ink"], weight="600"))


def draw_sites(elements: list[str], row: TrackRow, track_y: float, scale_length: int) -> None:
    if not row.site_positions:
        return
    placements: list[tuple[float, float]] = []
    offsets = [-18, -10, -2, -26, 6, 14]
    for position in row.site_positions:
        x = pos_to_x(position, scale_length)
        y = track_y - 18
        for candidate in offsets:
            proposed_y = y + candidate
            if all(abs(x - prior_x) > 11 or abs(proposed_y - prior_y) > 9 for prior_x, prior_y in placements):
                y = proposed_y
                break
        placements.append((x, y))

    for feature in row.regions:
        if feature.kind not in {"rg", "tail"}:
            continue
        region_sites = [pos for pos in row.site_positions if feature.start <= pos <= feature.end]
        if len(region_sites) >= 4:
            hx0 = pos_to_x(max(1, min(region_sites) - 2), scale_length)
            hx1 = pos_to_x(min(row.length, max(region_sites) + 2), scale_length)
            elements.append(rect(hx0, track_y - 42, max(18.0, hx1 - hx0), 30, rx=15, fill=COLORS["site_halo"], stroke="none", opacity="0.55"))

    for x, y in placements:
        elements.append(circle(x, y, 4.5, fill=COLORS["site_fill"], stroke=COLORS["site_edge"], stroke_width="1.1"))


def draw_g4_icon(elements: list[str], x: float, y: float) -> None:
    points = [
        (x + 16, y + 16),
        (x + 44, y + 16),
        (x + 16, y + 44),
        (x + 44, y + 44),
    ]
    for px, py in points:
        elements.append(rect(px - 8, py - 8, 16, 16, rx=4, fill="none", stroke=COLORS["target_edge"], stroke_width="2.0"))
    elements.append(line(points[0][0], points[0][1], points[1][0], points[1][1], stroke=COLORS["target_edge"], stroke_width="1.6"))
    elements.append(line(points[2][0], points[2][1], points[3][0], points[3][1], stroke=COLORS["target_edge"], stroke_width="1.6"))
    elements.append(line(points[0][0], points[0][1], points[2][0], points[2][1], stroke=COLORS["target_edge"], stroke_width="1.6"))
    elements.append(line(points[1][0], points[1][1], points[3][0], points[3][1], stroke=COLORS["target_edge"], stroke_width="1.6"))
    elements.append(path(wavy_segment(x - 4, x + 62, y + 68, amplitude=3.2, period=14.0), fill="none", stroke=COLORS["target_edge"], stroke_width="2.0"))


def draw_smn_reader(elements: list[str], x: float, y: float) -> None:
    elements.append(rect(x, y, 250, 70, rx=18, fill=COLORS["reader_fill"], stroke=COLORS["reader_edge"], stroke_width="1.2"))
    for px in (x + 58, x + 90, x + 122):
        elements.append(circle(px, y + 24, 4.8, fill=COLORS["reader_edge"], stroke="none"))
    elements.append(path(f"M {fmt(x + 50)} {fmt(y + 40)} Q {fmt(x + 92)} {fmt(y + 54)} {fmt(x + 134)} {fmt(y + 40)}", fill="none", stroke=COLORS["reader_edge"], stroke_width="2.0"))
    elements.append(multiline_text(x + 170, y + 28, ["SMN Tudor", "aromatic cage"], 11.3, anchor="middle", fill=COLORS["ink"], weight="600"))


def draw_fus_network(elements: list[str], x: float, y: float) -> None:
    elements.append(rect(x, y, 260, 78, rx=22, fill=COLORS["reader_fill"], stroke=COLORS["reader_edge"], stroke_width="1.2"))
    droplet = (
        f"M {fmt(x + 28)} {fmt(y + 46)} "
        f"C {fmt(x + 18)} {fmt(y + 22)} {fmt(x + 44)} {fmt(y + 8)} {fmt(x + 64)} {fmt(y + 18)} "
        f"C {fmt(x + 92)} {fmt(y + 8)} {fmt(x + 108)} {fmt(y + 38)} {fmt(x + 88)} {fmt(y + 56)} "
        f"C {fmt(x + 82)} {fmt(y + 78)} {fmt(x + 42)} {fmt(y + 74)} {fmt(x + 28)} {fmt(y + 46)} Z"
    )
    elements.append(path(droplet, fill="#C9DAF0", stroke="#89A6C7", stroke_width="1.1"))
    elements.append(path(wavy_segment(x + 32, x + 92, y + 46, amplitude=5.0, period=16.0), fill="none", stroke=COLORS["site_edge"], stroke_width="2.0"))
    for px in (x + 46, x + 60, x + 76):
        elements.append(circle(px, y + 46, 3.7, fill=COLORS["site_fill"], stroke="white", stroke_width="0.8"))
    elements.append(rect(x + 142, y + 24, 88, 28, rx=14, fill="white", stroke=COLORS["reader_edge"], stroke_width="1.1"))
    elements.append(text(x + 186, y + 42, "Kapβ2", 10.2, anchor="middle", fill=COLORS["muted"], weight="600"))
    elements.append(text(x + 154, y + 69, "LLPS / transport outputs", 10.6, fill=COLORS["muted"], weight="600"))


def draw_wdr5(elements: list[str], x: float, y: float, peptide: str) -> None:
    elements.append(rect(x, y, 250, 76, rx=20, fill=COLORS["reader_fill"], stroke=COLORS["reader_edge"], stroke_width="1.2"))
    elements.append(path(f"M {fmt(x + 34)} {fmt(y + 30)} Q {fmt(x + 78)} {fmt(y + 52)} {fmt(x + 122)} {fmt(y + 30)}", fill="none", stroke=COLORS["reader_edge"], stroke_width="2.4"))
    elements.append(text(x + 174, y + 30, "WDR5", 12.0, anchor="middle", fill=COLORS["ink"], weight="700"))
    elements.append(text(x + 174, y + 48, "WD-repeat pocket", 10.8, anchor="middle", fill=COLORS["muted"], weight="600"))
    inset_x = x - 4
    inset_y = y - 42
    elements.append(rect(inset_x, inset_y, 192, 30, rx=15, fill="white", stroke=COLORS["card_edge"], stroke_width="1"))
    chars = list(peptide)
    cursor = inset_x + 16
    for idx, char in enumerate(chars):
        fill = COLORS["site_fill"] if idx == 1 else COLORS["ink"]
        weight = "700" if idx == 1 else "400"
        elements.append(text(cursor, inset_y + 20, char, 11.6, fill=fill, weight=weight))
        cursor += 12.0
    elements.append(text(inset_x + 146, inset_y + 20, "H3R2", 10.6, fill=COLORS["muted"], weight="600"))


def draw_interactor(elements: list[str], row: TrackRow, row_y: float) -> None:
    if row.interactor_kind == "fus_network":
        draw_fus_network(elements, INTERACTOR_X, row_y + 20)
    elif row.interactor_kind == "smn":
        draw_smn_reader(elements, INTERACTOR_X, row_y + 22)
    elif row.interactor_kind == "g4":
        elements.append(rect(INTERACTOR_X, row_y + 18, 250, 84, rx=20, fill=COLORS["target_fill"], stroke="#BFD0E6", stroke_width="1.2"))
        draw_g4_icon(elements, INTERACTOR_X + 18, row_y + 28)
        label_lines = ["G-quadruplex", "RNA target"] if row.interactor_label == "G4 RNA" else ["G-quadruplex", "DNA / RNA target"]
        elements.append(multiline_text(INTERACTOR_X + 172, row_y + 47, label_lines, 11.5, anchor="middle", fill=COLORS["ink"], weight="600"))
    elif row.interactor_kind == "wdr5":
        draw_wdr5(elements, INTERACTOR_X, row_y + 24, row.peptide_inset)


def draw_tag(elements: list[str], row: TrackRow, row_y: float) -> None:
    tag_y = row_y + 108
    elements.append(rect(INTERACTOR_X, tag_y, 250, 28, rx=14, fill=COLORS["tag_fill"], stroke=COLORS["tag_edge"], stroke_width="1"))
    elements.append(text(INTERACTOR_X + 125, tag_y + 19, row.tag, 11.2, anchor="middle", fill=COLORS["ink"], weight="700"))
    elements.append(text(INTERACTOR_X + 125, tag_y + 46, row.tag_detail, 10.6, anchor="middle", fill=COLORS["muted"], weight="600"))


def draw_row(elements: list[str], row: TrackRow, row_y: float, scale_length: int, *, standalone: bool = False) -> None:
    card_x = 20 if standalone else ROW_CARD_X
    card_w = CARD_WIDTH if standalone else ROW_CARD_W
    label_x = card_x + 30
    track_x_shift = TRACK_X0 - 370
    interactor_shift = INTERACTOR_X - 1140

    if standalone:
        local_track_x0 = TRACK_X0 + card_x - ROW_CARD_X
        local_track_x1 = TRACK_X1 + card_x - ROW_CARD_X
        local_interactor_x = INTERACTOR_X + card_x - ROW_CARD_X
    else:
        local_track_x0 = TRACK_X0
        local_track_x1 = TRACK_X1
        local_interactor_x = INTERACTOR_X

    elements.append(rect(card_x, row_y, card_w, ROW_CARD_H, rx=26, fill=COLORS["card"], stroke=COLORS["card_edge"], stroke_width="1", filter="url(#softShadow)"))
    elements.append(text(label_x, row_y + 34, row.display_name, 21.0, weight="700"))
    elements.append(text(label_x, row_y + 58, row.descriptor, 12.2, fill=COLORS["muted"], weight="600"))

    baseline_y = row_y + 78
    row_track_end_x = pos_to_x(row.length, scale_length, local_track_x0, local_track_x1)
    elements.append(line(local_track_x0, baseline_y, row_track_end_x, baseline_y, stroke=COLORS["divider"], stroke_width="2.0"))

    for feature in row.regions:
        if feature.kind == "rg_span":
            continue
        x0 = pos_to_x(feature.start, scale_length, local_track_x0, local_track_x1)
        x1 = pos_to_x(feature.end, scale_length, local_track_x0, local_track_x1)
        if feature.kind in {"prld", "rg", "tail"}:
            elements.append(path(wavy_segment(x0, x1, baseline_y, amplitude=3.4, period=20.0), fill="none", stroke=COLORS["idr"], stroke_width="2.0"))

    original_track = (TRACK_X0, TRACK_X1)
    original_interactor = INTERACTOR_X
    globals()["TRACK_X0"], globals()["TRACK_X1"] = local_track_x0, local_track_x1
    globals()["INTERACTOR_X"] = local_interactor_x
    try:
        draw_track_regions(elements, row, baseline_y, scale_length, row_track_end_x)
        draw_domains(elements, row, baseline_y, scale_length)
        draw_sites(elements, row, baseline_y, scale_length)
        draw_interactor(elements, row, row_y)
        draw_tag(elements, row, row_y)
    finally:
        globals()["TRACK_X0"], globals()["TRACK_X1"] = original_track
        globals()["INTERACTOR_X"] = original_interactor

    elements.append(line(local_track_x0, baseline_y + 18, local_track_x0, baseline_y + 24, stroke=COLORS["divider"], stroke_width="1.4"))
    elements.append(line(row_track_end_x, baseline_y + 18, row_track_end_x, baseline_y + 24, stroke=COLORS["divider"], stroke_width="1.4"))
    elements.append(text(local_track_x0, baseline_y + 38, "1", 10.4, anchor="middle", fill=COLORS["muted"]))
    elements.append(text(row_track_end_x, baseline_y + 38, f"{row.length} aa", 10.4, anchor="middle", fill=COLORS["muted"]))


def draw_legend(elements: list[str], y: float) -> None:
    elements.append(line(52, y - 16, 1548, y - 16, stroke=COLORS["divider"], stroke_width="1.2"))

    legend_y = y + 18
    items = [
        ("site", "methyl-site"),
        ("domain", "ordered domain"),
        ("idr", "RG-rich IDR"),
        ("reader", "reader"),
        ("target", "G4 target"),
    ]
    cursor = 56.0
    gaps = {"site": 18.0, "domain": 24.0, "idr": 24.0, "reader": 24.0, "target": 0.0}
    for kind, label in items:
        if kind == "site":
            elements.append(circle(cursor + 5, legend_y + 2, 4.5, fill=COLORS["site_fill"], stroke=COLORS["site_edge"], stroke_width="1"))
            icon_end = cursor + 18
        elif kind == "domain":
            elements.append(rect(cursor, legend_y - 5, 34, 16, rx=8, fill=COLORS["domain_fill"], stroke=COLORS["domain_edge"], stroke_width="1"))
            icon_end = cursor + 42
        elif kind == "idr":
            elements.append(path(wavy_segment(cursor, cursor + 38, legend_y + 2, amplitude=3.0, period=14.0), fill="none", stroke=COLORS["idr"], stroke_width="2"))
            elements.append(rect(cursor + 6, legend_y - 4, 26, 7, rx=3.5, fill=COLORS["rg_fill"], stroke="none"))
            icon_end = cursor + 46
        elif kind == "reader":
            elements.append(rect(cursor, legend_y - 5, 34, 16, rx=8, fill=COLORS["reader_fill"], stroke=COLORS["reader_edge"], stroke_width="1"))
            icon_end = cursor + 42
        else:
            base_x = cursor + 2
            base_y = legend_y - 2
            for px, py in ((base_x + 8, base_y), (base_x + 24, base_y), (base_x + 8, base_y + 16), (base_x + 24, base_y + 16)):
                elements.append(rect(px - 4.5, py - 4.5, 9, 9, rx=2.5, fill="none", stroke=COLORS["target_edge"], stroke_width="1.3"))
            elements.append(line(base_x + 8, base_y, base_x + 24, base_y, stroke=COLORS["target_edge"], stroke_width="1.1"))
            elements.append(line(base_x + 8, base_y + 16, base_x + 24, base_y + 16, stroke=COLORS["target_edge"], stroke_width="1.1"))
            elements.append(line(base_x + 8, base_y, base_x + 8, base_y + 16, stroke=COLORS["target_edge"], stroke_width="1.1"))
            elements.append(line(base_x + 24, base_y, base_x + 24, base_y + 16, stroke=COLORS["target_edge"], stroke_width="1.1"))
            icon_end = cursor + 38
        elements.append(text(icon_end, legend_y + 6, label, 10.3, fill=COLORS["muted"], weight="600"))
        cursor = icon_end + 10 + len(label) * 5.55 + gaps[kind]

    adj_title_x = 1010
    adj_title_y = y + 8
    elements.append(text(adj_title_x, adj_title_y, "domain-edge adjacency", 10.8, fill=COLORS["muted"], weight="700"))

    domain_y = y + 30
    domain_x = 1060
    domain_w = 86
    elements.append(rect(domain_x, domain_y - 12, domain_w, 24, rx=12, fill=COLORS["domain_fill"], stroke=COLORS["domain_edge"], stroke_width="1.05"))

    idr_x0 = domain_x + domain_w + 8
    idr_x1 = 1498

    def scale(position: int) -> float:
        return idr_x0 + (idr_x1 - idr_x0) * ((position - 1) / 59.0)

    band20_x1 = scale(20)
    band40_x1 = scale(40)
    elements.append(rect(idr_x0, domain_y - 16, band20_x1 - idr_x0, 32, rx=16, fill=COLORS["adj_fill_near"], stroke="none"))
    elements.append(rect(band20_x1, domain_y - 16, band40_x1 - band20_x1, 32, rx=0, fill=COLORS["adj_fill_mid"], stroke="none"))
    elements.append(path(wavy_segment(idr_x0, idr_x1, domain_y, amplitude=4.0, period=18.0), fill="none", stroke=COLORS["idr"], stroke_width="2.0"))
    elements.append(line(idr_x0, domain_y - 21, idr_x0, domain_y + 21, stroke=COLORS["reader_edge"], stroke_width="1.2"))
    elements.append(line(band20_x1, domain_y - 18, band20_x1, domain_y + 18, stroke=COLORS["reader_edge"], stroke_width="1.0", stroke_dasharray="3 4", opacity="0.7"))
    elements.append(line(band40_x1, domain_y - 18, band40_x1, domain_y + 18, stroke=COLORS["reader_edge"], stroke_width="1.0", stroke_dasharray="3 4", opacity="0.7"))
    elements.append(text(idr_x0 + 24, domain_y - 23, "<=20", 9.5, anchor="middle", fill=COLORS["reader_edge"], weight="700"))
    elements.append(text((band20_x1 + band40_x1) / 2, domain_y - 23, "21-40", 9.5, anchor="middle", fill=COLORS["reader_edge"], weight="700"))
    elements.append(text((band40_x1 + idr_x1) / 2, domain_y - 23, ">40", 9.5, anchor="middle", fill=COLORS["reader_edge"], weight="700"))

    for idx, pos in enumerate((7, 13, 28, 34, 49)):
        site_y = domain_y + (-5 if idx % 2 == 0 else 6)
        elements.append(circle(scale(pos), site_y, 4.0, fill=COLORS["site_fill"], stroke=COLORS["site_edge"], stroke_width="0.95"))


def build_figure_svg(rows: list[TrackRow]) -> str:
    scale_length = max(row.length for row in rows)
    elements: list[str] = [
        rect(0, 0, FIGURE_WIDTH, FIGURE_HEIGHT, fill=COLORS["bg"], stroke="none"),
        text(52, 58, "d.", 24.0, weight="700"),
        text(98, 58, "Representative architectures for methylarginine-dependent interaction redistribution", 22.0, weight="700"),
        text(98, 84, "Evidence-guided deployment across clustered IDRs, boundary-proximal RG tails, domain-adjacent RGG modules, and single-site chromatin marks.", 12.4, fill=COLORS["muted"], weight="600"),
    ]
    chips = [
        "RG/RGG/GAR motifs enriched",
        "IDR and domain-adjacent deployment",
        "local site clustering observed",
        "RNA / RNP classes enriched",
    ]
    x = 98.0
    for label in chips:
        draw_chip(elements, x, 104, label)
        x += max(158.0, 24.0 + len(label) * 7.0) + 12.0

    for idx, row in enumerate(rows):
        draw_row(elements, row, 164 + idx * 176, scale_length)

    draw_legend(elements, 934)
    body = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{FIGURE_WIDTH}" height="{FIGURE_HEIGHT}" viewBox="0 0 {FIGURE_WIDTH} {FIGURE_HEIGHT}">',
        defs_block(),
        style_block(),
        *elements,
        "</svg>",
    ]
    return "\n".join(body)


def build_row_card_svg(row: TrackRow, scale_length: int) -> str:
    width = 1560
    height = CARD_HEIGHT
    elements: list[str] = [rect(0, 0, width, height, fill=COLORS["bg"], stroke="none")]
    draw_row(elements, row, SINGLE_ROW_Y, scale_length, standalone=True)
    body = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        defs_block(),
        style_block(),
        *elements,
        "</svg>",
    ]
    return "\n".join(body)


def main() -> None:
    parser = argparse.ArgumentParser(description="Create representative methylarginine substrate architecture SVGs.")
    parser.add_argument("--site-table", default="PTM_results/domain_context/sites_with_domain_context.tsv")
    parser.add_argument("--remapped-rows", default="PTM_results/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv")
    parser.add_argument("--interpro", default="PTM_data/context/interpro_human_reviewed_domain_like_intervals.tsv")
    parser.add_argument("--disorder", default="PTM_data/context/mobidb_human_reviewed_disorder_intervals.tsv")
    parser.add_argument("--metadata", default="PTM_data/context/uniprot_human_reviewed_canonical_metadata.tsv")
    parser.add_argument("--fasta", default="PTM_data/context/uniprot_human_reviewed_canonical.fasta")
    parser.add_argument("--outdir", default="PTM_results/summary_figures")
    args = parser.parse_args()

    site_table = pd.read_csv(args.site_table, sep="\t", low_memory=False)
    remapped_rows = pd.read_csv(args.remapped_rows, sep="\t", low_memory=False)
    interpro = pd.read_csv(args.interpro, sep="\t", low_memory=False)
    metadata = pd.read_csv(args.metadata, sep="\t", low_memory=False)
    disorder_by_acc = load_disorder_regions(Path(args.disorder))
    sequences = parse_fasta(args.fasta)

    rows = build_rows(site_table, remapped_rows, metadata, interpro, disorder_by_acc, sequences)
    scale_length = max(row.length for row in rows)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    combined_path = outdir / "representative_methylarginine_architectures.svg"
    combined_path.write_text(build_figure_svg(rows))
    print(combined_path)

    for row in rows:
        row_path = outdir / f"{row.panel_id}_architecture_card.svg"
        row_path.write_text(build_row_card_svg(row, scale_length))
        print(row_path)


if __name__ == "__main__":
    main()

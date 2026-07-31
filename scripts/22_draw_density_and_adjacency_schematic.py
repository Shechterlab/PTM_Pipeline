from __future__ import annotations

import argparse
from html import escape
from pathlib import Path


WIDTH = 1360
HEIGHT = 920

COLORS = {
    "bg": "#FFFFFF",
    "ink": "#1A2433",
    "muted": "#637083",
    "divider": "#E6ECF2",
    "idr": "#8C96A3",
    "window_edge": "#5D82AF",
    "window_dash": "#9EB6D2",
    "site_fill": "#CBB2B9",
    "site_edge": "#8E747C",
    "domain_fill": "#E9DFC9",
    "domain_edge": "#736858",
    "accent": "#4F79A8",
    "accent_light": "#ECF3FA",
    "accent_mid": "#F4F8FC",
    "chip_fill": "#F5F8FB",
    "chip_edge": "#D9E1EA",
    "zone_20": "#EAF2FB",
    "zone_40": "#F4F8FC",
}


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


def circle(cx: float, cy: float, r: float, **kwargs: object) -> str:
    return f'<circle {attrs(cx=fmt(cx), cy=fmt(cy), r=fmt(r), **kwargs)} />'


def line(x1: float, y1: float, x2: float, y2: float, **kwargs: object) -> str:
    return f'<line {attrs(x1=fmt(x1), y1=fmt(y1), x2=fmt(x2), y2=fmt(y2), **kwargs)} />'


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
) -> str:
    return f'<text {attrs(x=fmt(x), y=fmt(y), fill=fill, font_size=fmt(font_size), text_anchor=anchor, font_weight=weight)}>{escape(content)}</text>'


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


def wavy_segment(x0: float, x1: float, y: float, amplitude: float = 4.0, period: float = 22.0) -> str:
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


def chip(x: float, y: float, label: str) -> str:
    width = max(168.0, 24.0 + len(label) * 7.0)
    return "\n".join([
        rect(x, y, width, 32, rx=16, fill=COLORS["chip_fill"], stroke=COLORS["chip_edge"], stroke_width="1"),
        text(x + width / 2, y + 21, label, 11.4, anchor="middle", fill=COLORS["muted"], weight="600"),
    ])


def draw_window_bracket(elements: list[str], x0: float, x1: float, y: float, count_label: str) -> None:
    elements.append(line(x0, y, x1, y, stroke=COLORS["window_edge"], stroke_width="2.0"))
    elements.append(line(x0, y - 10, x0, y + 10, stroke=COLORS["window_edge"], stroke_width="2.0"))
    elements.append(line(x1, y - 10, x1, y + 10, stroke=COLORS["window_edge"], stroke_width="2.0"))
    mid = (x0 + x1) / 2
    elements.append(text(mid, y - 16, "rolling 20 aa span", 10.8, anchor="middle", fill=COLORS["accent"], weight="700"))
    elements.append(rect(x1 + 22, y - 16, 62, 34, rx=17, fill="white", stroke=COLORS["window_edge"], stroke_width="1.2"))
    elements.append(text(x1 + 53, y + 7, count_label, 14.6, anchor="middle", fill=COLORS["ink"], weight="700"))
    elements.append(text(x1 + 53, y + 22, "sites", 9.5, anchor="middle", fill=COLORS["muted"], weight="600"))


def draw_density_row(
    elements: list[str],
    y: float,
    title: str,
    subtitle: str,
    positions: list[int],
    *,
    highlight_start: int,
    highlight_end: int,
    count_label: str,
) -> None:
    x0 = 310
    x1 = 1110
    scale = lambda position: x0 + (x1 - x0) * ((position - 1) / 59.0)
    elements.append(text(72, y - 12, title, 18.0, weight="700"))
    elements.append(text(72, y + 11, subtitle, 11.7, fill=COLORS["muted"], weight="600"))

    guide_y = y + 4
    elements.append(path(wavy_segment(x0, x1, guide_y, amplitude=5.0, period=24.0), fill="none", stroke=COLORS["idr"], stroke_width="2.5"))
    window_x0 = scale(highlight_start)
    window_x1 = scale(highlight_end)
    elements.append(line(window_x0, guide_y - 28, window_x0, guide_y + 34, stroke=COLORS["window_dash"], stroke_width="1.4", stroke_dasharray="4 4"))
    elements.append(line(window_x1, guide_y - 28, window_x1, guide_y + 34, stroke=COLORS["window_dash"], stroke_width="1.4", stroke_dasharray="4 4"))
    draw_window_bracket(elements, window_x0, window_x1, guide_y + 42, count_label)

    for idx, pos in enumerate(positions):
        y_offset = -5 if idx % 2 == 0 else 6
        elements.append(circle(scale(pos), guide_y + y_offset, 5.0, fill=COLORS["site_fill"], stroke=COLORS["site_edge"], stroke_width="1.0"))

    elements.append(line(x0, guide_y + 70, x0, guide_y + 76, stroke=COLORS["divider"], stroke_width="1.2"))
    elements.append(line(x1, guide_y + 70, x1, guide_y + 76, stroke=COLORS["divider"], stroke_width="1.2"))
    elements.append(text(x0, guide_y + 94, "same 60 aa IDR segment", 10.2, anchor="start", fill=COLORS["muted"], weight="600"))


def draw_adjacency_card(
    elements: list[str],
    x: float,
    y: float,
    title: str,
    zone_label: str,
    positions: list[int],
    *,
    highlight: str,
) -> None:
    card_w = 396
    card_h = 214
    elements.append(rect(x, y, card_w, card_h, rx=24, fill="white", stroke=COLORS["chip_edge"], stroke_width="1"))
    elements.append(text(x + 22, y + 30, title, 16.0, weight="700"))
    elements.append(text(x + 22, y + 52, zone_label, 11.0, fill=COLORS["muted"], weight="600"))

    domain_x = x + 26
    domain_y = y + 102
    elements.append(rect(domain_x, domain_y - 18, 88, 36, rx=18, fill=COLORS["domain_fill"], stroke=COLORS["domain_edge"], stroke_width="1.2"))
    elements.append(text(domain_x + 44, domain_y + 5, "ordered", 10.4, anchor="middle", fill=COLORS["muted"], weight="600"))
    elements.append(text(domain_x + 44, domain_y + 18, "domain", 10.4, anchor="middle", fill=COLORS["muted"], weight="600"))

    idr_x0 = domain_x + 92
    idr_x1 = x + 356
    scale = lambda position: idr_x0 + (idr_x1 - idr_x0) * ((position - 1) / 59.0)
    band20_x1 = scale(20)
    band40_x1 = scale(40)
    elements.append(rect(idr_x0, domain_y - 28, band20_x1 - idr_x0, 58, rx=16, fill=COLORS["zone_20"], stroke="none"))
    elements.append(rect(band20_x1, domain_y - 28, band40_x1 - band20_x1, 58, rx=0, fill=COLORS["zone_40"], stroke="none"))
    elements.append(path(wavy_segment(idr_x0, idr_x1, domain_y, amplitude=5.0, period=20.0), fill="none", stroke=COLORS["idr"], stroke_width="2.5"))
    elements.append(line(idr_x0, domain_y - 34, idr_x0, domain_y + 34, stroke=COLORS["accent"], stroke_width="1.3"))
    elements.append(line(band20_x1, domain_y - 30, band20_x1, domain_y + 30, stroke=COLORS["window_dash"], stroke_width="1.1", stroke_dasharray="4 4"))
    elements.append(line(band40_x1, domain_y - 30, band40_x1, domain_y + 30, stroke=COLORS["window_dash"], stroke_width="1.1", stroke_dasharray="4 4"))
    elements.append(text(idr_x0 + 36, domain_y - 38, "<=20 aa", 9.8, anchor="middle", fill=COLORS["accent"], weight="700"))
    elements.append(text((band20_x1 + band40_x1) / 2, domain_y - 38, "21-40 aa", 9.8, anchor="middle", fill=COLORS["accent"], weight="700"))
    elements.append(text((band40_x1 + idr_x1) / 2, domain_y - 38, ">40 aa", 9.8, anchor="middle", fill=COLORS["accent"], weight="700"))

    for idx, pos in enumerate(positions):
        y_offset = -6 if idx % 2 == 0 else 7
        elements.append(circle(scale(pos), domain_y + y_offset, 5.0, fill=COLORS["site_fill"], stroke=COLORS["site_edge"], stroke_width="1.0"))

    if highlight == "edge":
        hx0 = scale(min(positions) - 1)
        hx1 = scale(max(positions) + 1)
        elements.append(path(f"M {fmt(hx0)} {fmt(domain_y + 34)} Q {fmt((hx0 + hx1) / 2)} {fmt(domain_y + 56)} {fmt(hx1)} {fmt(domain_y + 34)}", fill="none", stroke=COLORS["accent"], stroke_width="1.8"))
        elements.append(text((hx0 + hx1) / 2, domain_y + 74, "edge-proximal cluster", 10.2, anchor="middle", fill=COLORS["muted"], weight="600"))
    elif highlight == "mid":
        hx0 = scale(min(positions) - 1)
        hx1 = scale(max(positions) + 1)
        elements.append(path(f"M {fmt(hx0)} {fmt(domain_y + 34)} Q {fmt((hx0 + hx1) / 2)} {fmt(domain_y + 56)} {fmt(hx1)} {fmt(domain_y + 34)}", fill="none", stroke=COLORS["accent"], stroke_width="1.8"))
        elements.append(text((hx0 + hx1) / 2, domain_y + 74, "intermediate-distance cluster", 10.2, anchor="middle", fill=COLORS["muted"], weight="600"))
    else:
        hx0 = scale(min(positions) - 1)
        hx1 = scale(max(positions) + 1)
        elements.append(path(f"M {fmt(hx0)} {fmt(domain_y + 34)} Q {fmt((hx0 + hx1) / 2)} {fmt(domain_y + 56)} {fmt(hx1)} {fmt(domain_y + 34)}", fill="none", stroke=COLORS["accent"], stroke_width="1.8"))
        elements.append(text((hx0 + hx1) / 2, domain_y + 74, "domain-distal cluster", 10.2, anchor="middle", fill=COLORS["muted"], weight="600"))


def build_svg() -> str:
    elements: list[str] = [
        rect(0, 0, WIDTH, HEIGHT, fill=COLORS["bg"], stroke="none"),
        text(52, 52, "Local methylarginine density and structured-domain adjacency", 24.0, weight="700"),
        text(52, 80, "Definitions for interpreting methylarginine spacing and distance from the nearest ordered-domain edge.", 13.0, fill=COLORS["muted"], weight="600"),
        chip(52, 104, "same site chemistry"),
        chip(242, 104, "same total site count"),
        chip(446, 104, "spacing defines density"),
        chip(646, 104, "domain-edge distance defines adjacency"),
    ]

    elements.append(text(52, 162, "Local density cartoon", 18.0, weight="700"))
    elements.append(text(52, 183, "Keep the IDR loose and schematic; use a bracketed span, not a box, to avoid looking like an ordered domain.", 11.8, fill=COLORS["muted"], weight="600"))

    draw_density_row(
        elements,
        228,
        "Diffuse spacing",
        "Same total number of methyl-sites, but spread across the IDR.",
        [5, 14, 24, 35, 46, 57],
        highlight_start=20,
        highlight_end=39,
        count_label="2",
    )
    draw_density_row(
        elements,
        368,
        "Intermediate clustering",
        "Several methyl-sites fall into one span without forming a single tight block.",
        [7, 18, 26, 30, 33, 46],
        highlight_start=20,
        highlight_end=39,
        count_label="4",
    )
    draw_density_row(
        elements,
        508,
        "High local density",
        "Clustered sites create a compact multivalent patch even when the total site count is unchanged.",
        [6, 24, 27, 29, 31, 33],
        highlight_start=20,
        highlight_end=39,
        count_label="5",
    )

    elements.append(line(52, 620, 1308, 620, stroke=COLORS["divider"], stroke_width="1.2"))

    elements.append(text(52, 656, "Structured-domain adjacency cartoon", 18.0, weight="700"))
    elements.append(text(52, 677, "A companion schematic can show how disordered methyl-sites are positioned relative to the nearest ordered-domain edge.", 11.8, fill=COLORS["muted"], weight="600"))

    draw_adjacency_card(
        elements,
        52,
        700,
        "Edge-proximal",
        "disordered sites within <=20 aa of the domain edge",
        [6, 10, 14],
        highlight="edge",
    )
    draw_adjacency_card(
        elements,
        482,
        700,
        "Intermediate distance",
        "disordered sites in the 21-40 aa shell",
        [24, 28, 33],
        highlight="mid",
    )
    draw_adjacency_card(
        elements,
        912,
        700,
        "Domain-distal",
        "disordered sites >40 aa from the nearest domain edge",
        [45, 50, 55],
        highlight="distal",
    )

    elements.append(rect(52, 890, 1256, 24, rx=12, fill="#F7F9FC", stroke=COLORS["chip_edge"], stroke_width="1"))
    elements.append(
        text(
            680,
            907,
            "Same chemistry and same total site count can yield different local-density and domain-adjacency architectures.",
            11.6,
            anchor="middle",
            fill=COLORS["muted"],
            weight="700",
        )
    )

    return "\n".join(
        [
            '<?xml version="1.0" encoding="UTF-8"?>',
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}">',
            """
<style>
  svg { background: #FFFFFF; }
  text {
    font-family: Arial, "Liberation Sans", "DejaVu Sans", Helvetica, sans-serif;
    letter-spacing: 0.01em;
  }
</style>
""".strip(),
            *elements,
            "</svg>",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Create density and domain-adjacency schematic SVGs.")
    parser.add_argument("--outdir", default="PTM_results/summary_figures")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / "density_and_domain_adjacency_schematic.svg"
    outpath.write_text(build_svg())
    print(outpath)


if __name__ == "__main__":
    main()

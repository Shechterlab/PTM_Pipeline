from __future__ import annotations

import argparse
from html import escape
from pathlib import Path


WIDTH = 1600
HEIGHT = 1040

COLORS = {
    "bg": "#FFFFFF",
    "ink": "#162131",
    "muted": "#5E6B7D",
    "hub_fill": "#EEF3F8",
    "hub_edge": "#D3DDE9",
    "arg_fill": "#6F8FBF",
    "arg_edge": "#4E6F9F",
    "me_fill": "#CF6B86",
    "me_edge": "#944960",
    "reader_fill": "#DCEAF7",
    "reader_edge": "#7E9BB9",
    "target_fill": "#E4EEFB",
    "target_edge": "#4F79A8",
    "cond_fill": "#E9EDF2",
    "cond_edge": "#8B9BAE",
    "lavender_fill": "#F4E9F9",
    "lavender_edge": "#B188C4",
    "card_edge": "#E1E7EF",
    "shadow": "#A8B7CA",
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
    return f'<text {attrs(x=fmt(x), y=fmt(y), fill=fill, font_size=fmt(font_size), text_anchor=anchor, font_weight=weight, font_style=("italic" if italic else "normal"))}>{escape(content)}</text>'


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


def wavy_segment(x0: float, x1: float, y: float, amplitude: float = 4.0, period: float = 24.0) -> str:
    if x1 <= x0:
        return ""
    segments = max(1, int((x1 - x0) / period))
    step = (x1 - x0) / segments
    parts = [f"M {fmt(x0)} {fmt(y)}"]
    cursor = x0
    for idx in range(segments):
        nxt = x1 if idx == segments - 1 else cursor + step
        mid = (cursor + nxt) / 2
        peak = y - amplitude if idx % 2 == 0 else y + amplitude
        parts.append(f"Q {fmt(mid)} {fmt(peak)} {fmt(nxt)} {fmt(y)}")
        cursor = nxt
    return " ".join(parts)


def defs_block() -> str:
    return """
<defs>
  <filter id="softShadow" x="-20%" y="-20%" width="140%" height="160%">
    <feDropShadow dx="0" dy="4" stdDeviation="5" flood-color="#A8B7CA" flood-opacity="0.18" />
  </filter>
  <marker id="arrowHead" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="8" markerHeight="8" orient="auto-start-reverse">
    <path d="M 0 0 L 10 5 L 0 10 z" fill="#5E6B7D" />
  </marker>
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


def draw_g4_icon(elements: list[str], x: float, y: float) -> None:
    points = [
        (x + 20, y + 18),
        (x + 50, y + 18),
        (x + 20, y + 48),
        (x + 50, y + 48),
    ]
    for px, py in points:
        elements.append(rect(px - 9, py - 9, 18, 18, rx=4, fill="white", stroke=COLORS["target_edge"], stroke_width="2"))
    elements.append(line(points[0][0], points[0][1], points[1][0], points[1][1], stroke=COLORS["target_edge"], stroke_width="1.8"))
    elements.append(line(points[2][0], points[2][1], points[3][0], points[3][1], stroke=COLORS["target_edge"], stroke_width="1.8"))
    elements.append(line(points[0][0], points[0][1], points[2][0], points[2][1], stroke=COLORS["target_edge"], stroke_width="1.8"))
    elements.append(line(points[1][0], points[1][1], points[3][0], points[3][1], stroke=COLORS["target_edge"], stroke_width="1.8"))
    elements.append(path(wavy_segment(x - 8, x + 72, y + 74, amplitude=3.0, period=16), fill="none", stroke=COLORS["target_edge"], stroke_width="2.2"))


def draw_reader_icon(elements: list[str], x: float, y: float) -> None:
    elements.append(rect(x, y + 8, 76, 44, rx=16, fill="white", stroke=COLORS["reader_edge"], stroke_width="1.4"))
    for px in (x + 22, x + 36, x + 50):
        elements.append(circle(px, y + 24, 4.8, fill=COLORS["reader_edge"], stroke="none"))
    elements.append(path(f"M {fmt(x + 18)} {fmt(y + 38)} Q {fmt(x + 36)} {fmt(y + 48)} {fmt(x + 54)} {fmt(y + 38)}", fill="none", stroke=COLORS["reader_edge"], stroke_width="2.0"))
    elements.append(rect(x + 88, y + 14, 80, 30, rx=15, fill="white", stroke=COLORS["reader_edge"], stroke_width="1.2"))
    elements.append(path(f"M {fmt(x + 102)} {fmt(y + 29)} Q {fmt(x + 128)} {fmt(y + 42)} {fmt(x + 154)} {fmt(y + 29)}", fill="none", stroke=COLORS["reader_edge"], stroke_width="2.0"))


def draw_condensate_icon(elements: list[str], x: float, y: float) -> None:
    droplet = (
        f"M {fmt(x + 18)} {fmt(y + 54)} "
        f"C {fmt(x + 8)} {fmt(y + 20)} {fmt(x + 42)} {fmt(y + 2)} {fmt(x + 70)} {fmt(y + 18)} "
        f"C {fmt(x + 112)} {fmt(y + 8)} {fmt(x + 132)} {fmt(y + 48)} {fmt(x + 108)} {fmt(y + 74)} "
        f"C {fmt(x + 98)} {fmt(y + 102)} {fmt(x + 46)} {fmt(y + 98)} {fmt(x + 18)} {fmt(y + 54)} Z"
    )
    elements.append(path(droplet, fill=COLORS["cond_fill"], stroke=COLORS["cond_edge"], stroke_width="1.4"))
    elements.append(path(wavy_segment(x + 24, x + 104, y + 54, amplitude=5.0, period=18), fill="none", stroke=COLORS["arg_edge"], stroke_width="2.0"))
    for px in (x + 38, x + 52, x + 68, x + 86):
        elements.append(circle(px, y + 54, 4.0, fill=COLORS["me_fill"], stroke="white", stroke_width="0.8"))
    for px in (x + 44, x + 78):
        elements.append(circle(px, y + 34, 4.2, fill=COLORS["arg_fill"], stroke="none"))


def draw_hub(elements: list[str], cx: float, cy: float) -> None:
    elements.append(circle(cx, cy, 128, fill=COLORS["hub_fill"], stroke=COLORS["hub_edge"], stroke_width="1.2"))
    blob = (
        f"M {fmt(cx - 72)} {fmt(cy + 18)} "
        f"C {fmt(cx - 102)} {fmt(cy - 42)} {fmt(cx - 38)} {fmt(cy - 108)} {fmt(cx + 28)} {fmt(cy - 82)} "
        f"C {fmt(cx + 92)} {fmt(cy - 98)} {fmt(cx + 122)} {fmt(cy - 22)} {fmt(cx + 84)} {fmt(cy + 40)} "
        f"C {fmt(cx + 70)} {fmt(cy + 94)} {fmt(cx - 8)} {fmt(cy + 110)} {fmt(cx - 72)} {fmt(cy + 18)} Z"
    )
    elements.append(path(blob, fill="#E2E7ED", stroke="#BCC7D3", stroke_width="1.2"))
    elements.append(path(wavy_segment(cx - 46, cx + 66, cy + 8, amplitude=7.0, period=20), fill="none", stroke=COLORS["arg_edge"], stroke_width="3.0"))
    for px in (cx - 16, cx + 8, cx + 34, cx + 58):
        elements.append(circle(px, cy + 8, 5.2, fill=COLORS["arg_fill"], stroke="white", stroke_width="0.8"))
    for px in (cx - 2, cx + 24, cx + 46):
        elements.append(circle(px, cy + 8, 5.2, fill=COLORS["me_fill"], stroke="white", stroke_width="0.8"))
    elements.append(rect(cx - 166, cy - 168, 332, 34, rx=17, fill="#F7F9FC", stroke=COLORS["card_edge"], stroke_width="1"))
    elements.append(text(cx, cy - 146, "Arg -> Rme1 / Rme2a / Rme2s", 14.0, anchor="middle", fill=COLORS["ink"], weight="700"))
    elements.append(text(cx, cy + 150, "RGG/RG-rich IDR hub", 16.0, anchor="middle", fill=COLORS["ink"], weight="700"))
    elements.append(text(cx, cy + 172, "Interaction redistribution", 12.8, anchor="middle", fill=COLORS["muted"], weight="700"))


def draw_card(elements: list[str], x: float, y: float, width: float, height: float, fill: str, edge: str, title: str, body: list[str]) -> None:
    elements.append(rect(x, y, width, height, rx=26, fill=fill, stroke=edge, stroke_width="1.3", filter="url(#softShadow)"))
    elements.append(text(x + 24, y + 36, title, 20.0, fill=COLORS["ink"], weight="700"))
    elements.append(multiline_text(x + 24, y + 76, body, 13.0, fill=COLORS["muted"], weight="600", line_height=18.0))


def build_svg() -> str:
    elements: list[str] = [
        rect(0, 0, WIDTH, HEIGHT, fill=COLORS["bg"], stroke="none"),
        text(60, 60, "Methylarginine redistributes interaction capacity across nucleic-acid, reader, and IDR networks", 25.0, weight="700"),
        text(60, 86, "Charge is retained, but donor topology and methyl-face geometry are reweighted across competing interaction outputs.", 13.0, fill=COLORS["muted"], weight="600"),
    ]

    left = (92, 182, 396, 246)
    right = (1112, 182, 396, 246)
    bottom = (602, 692, 396, 232)

    draw_card(
        elements,
        *left,
        COLORS["lavender_fill"],
        COLORS["lavender_edge"],
        "Nucleic-acid binding",
        [
            "RGG/RG multivalency",
            "phosphate and base-edge contacts",
            "Arg/Rme-base or aromatic contacts",
        ],
    )
    draw_g4_icon(elements, 338, 268)
    elements.append(text(374, 368, "RNA / DNA targets", 11.5, anchor="middle", fill=COLORS["muted"], weight="600"))

    draw_card(
        elements,
        *right,
        COLORS["reader_fill"],
        COLORS["reader_edge"],
        "Reader / protein recruitment",
        [
            "methyl-compatible binding surfaces",
            "Tudor domains and WD-repeat pockets",
            "transport or other methyl-sensitive interfaces",
        ],
    )
    draw_reader_icon(elements, 1324, 258)
    elements.append(text(1394, 368, "state-selective compatibility", 11.5, anchor="middle", fill=COLORS["muted"], weight="600"))

    draw_card(
        elements,
        *bottom,
        "#EEF2F6",
        COLORS["cond_edge"],
        "IDR self-association / condensates",
        [
            "Arg-Tyr cation-pi interactions",
            "Arg/Rme-aromatic contacts",
            "RGG/RG network tuning and LLPS",
        ],
    )
    draw_condensate_icon(elements, 842, 742)
    elements.append(text(876, 866, "material-state tuning", 11.5, anchor="middle", fill=COLORS["muted"], weight="600"))

    draw_hub(elements, 800, 446)

    elements.append(path(f"M 690 416 C 618 366 560 340 492 322", fill="none", stroke=COLORS["muted"], stroke_width="2.4", marker_end="url(#arrowHead)"))
    elements.append(path(f"M 910 416 C 982 366 1040 340 1108 322", fill="none", stroke=COLORS["muted"], stroke_width="2.4", marker_end="url(#arrowHead)"))
    elements.append(path(f"M 800 574 C 800 624 800 652 800 690", fill="none", stroke=COLORS["muted"], stroke_width="2.4", marker_end="url(#arrowHead)"))

    elements.append(rect(194, 944, 1212, 38, rx=19, fill="#F7F9FC", stroke=COLORS["card_edge"], stroke_width="1"))
    elements.append(
        text(
            800,
            968,
            "Same Arg-rich surface, different dominant outputs: RNA engagement, reader recruitment, or network-driven condensation.",
            13.0,
            anchor="middle",
            fill=COLORS["muted"],
            weight="700",
        )
    )

    body = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}">',
        defs_block(),
        style_block(),
        *elements,
        "</svg>",
    ]
    return "\n".join(body)


def main() -> None:
    parser = argparse.ArgumentParser(description="Create an interaction-redistribution model SVG.")
    parser.add_argument("--outdir", default="PTM_results/summary_figures")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / "interaction_redistribution_model.svg"
    outpath.write_text(build_svg())
    print(outpath)


if __name__ == "__main__":
    main()

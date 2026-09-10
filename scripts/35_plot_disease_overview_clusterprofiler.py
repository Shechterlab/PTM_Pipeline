#!/usr/bin/env python3
"""ClusterProfiler-style disease-family dot plot from the audited family table."""

import csv
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from reportlab.pdfgen import canvas

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "PTM_results/functional_class_union/disease_enrichment"
SOURCE = OUT / "all_methyl_proteins_disease_overview_audited_v1.tsv"
PREFIX = None
SIZE_MODE = "equal"

WIDTH, HEIGHT = 1280, 760
LEFT, RIGHT, TOP, BOTTOM = 430, 260, 78, 88
XMIN, XMAX = 1.0, 2.5
LIGHT = (215, 227, 234)
DARK = (31, 78, 121)


def load_rows():
    with SOURCE.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        row["odds_ratio"] = float(row["odds_ratio"])
        row["BH_adjusted_p"] = float(row["BH_adjusted_p"])
        row["methylarginine_protein_count"] = int(row["methylarginine_protein_count"])
        row["minus_log10_BH_adjusted_p"] = -math.log10(row["BH_adjusted_p"])
    return sorted(rows, key=lambda row: row["odds_ratio"], reverse=True)


def x_pos(value):
    return LEFT + (value - XMIN) / (XMAX - XMIN) * (WIDTH - LEFT - RIGHT)


def radius(count, counts):
    if SIZE_MODE == "equal":
        return 12
    lo, hi = math.sqrt(min(counts)), math.sqrt(max(counts))
    return 10 + (math.sqrt(count) - lo) / (hi - lo) * 6


def color(q):
    score = min(60.0, max(10.0, -math.log10(q)))
    t = (score - 10.0) / 50.0
    return tuple(round(LIGHT[i] + t * (DARK[i] - LIGHT[i])) for i in range(3))


def rgb_hex(rgb):
    return "#%02x%02x%02x" % rgb


def y_positions(rows):
    gap = (HEIGHT - TOP - BOTTOM) / len(rows)
    return [TOP + (i + 0.48) * gap for i in range(len(rows))]


def write_table(rows):
    fields = list(rows[0])
    with PREFIX.with_suffix(".tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def draw_png(rows):
    image = Image.new("RGB", (WIDTH, HEIGHT), "white")
    draw = ImageDraw.Draw(image)
    regular = lambda size: ImageFont.truetype("arial.ttf", size)
    bold = lambda size: ImageFont.truetype("arialbd.ttf", size)
    draw.text((LEFT, 18), "Disease ontology enrichment among methylarginine proteins", font=bold(25), fill="#1d1d1d")
    for tick in (1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5):
        x = x_pos(tick)
        draw.line((x, TOP - 12, x, HEIGHT - BOTTOM), fill="#e8e8e8", width=1)
        label = f"{tick:g}"
        draw.text((x - draw.textlength(label, font=regular(14)) / 2, HEIGHT - BOTTOM + 12), label, font=regular(14), fill="#222222")
    draw.line((x_pos(1), TOP - 12, x_pos(1), HEIGHT - BOTTOM), fill="#555555", width=2)
    counts = [row["methylarginine_protein_count"] for row in rows]
    for row, y in zip(rows, y_positions(rows)):
        label = row["disease_family"]
        draw.text((LEFT - 18 - draw.textlength(label, font=regular(19)), y - 11), label, font=regular(19), fill="#1d1d1d")
        r = radius(row["methylarginine_protein_count"], counts)
        x = x_pos(row["odds_ratio"])
        draw.ellipse((x-r, y-r, x+r, y+r), fill=color(row["BH_adjusted_p"]), outline="white", width=2)
    axis = "Odds ratio vs reviewed-human background"
    draw.text((LEFT + (WIDTH-LEFT-RIGHT)/2 - draw.textlength(axis, font=regular(19))/2, HEIGHT-32), axis, font=regular(19), fill="#1d1d1d")
    draw_legends_png(draw, rows, regular, bold)
    image.save(PREFIX.with_suffix(".png"), dpi=(300, 300))


def draw_legends_png(draw, rows, regular, bold):
    lx = WIDTH - RIGHT + 48
    draw.text((lx, 105), "BH-adjusted p", font=bold(16), fill="#1d1d1d")
    for i, q in enumerate((1e-10, 1e-25, 1e-50)):
        y = 143 + i * 36
        draw.ellipse((lx, y-8, lx+16, y+8), fill=color(q))
        draw.text((lx+27, y-10), f"1×10^{int(math.log10(q))}", font=regular(14), fill="#1d1d1d")
    if SIZE_MODE == "equal":
        return
    draw.text((lx, 285), "Methylarginine proteins", font=bold(16), fill="#1d1d1d")
    counts = [row["methylarginine_protein_count"] for row in rows]
    for i, count in enumerate((500, 1500, 3500)):
        y = 330 + i * 53
        r = radius(count, counts)
        draw.ellipse((lx+18-r, y-r, lx+18+r, y+r), fill="#416f96", outline="white", width=2)
        draw.text((lx+50, y-9), f"{count:,}", font=regular(14), fill="#1d1d1d")


def draw_svg(rows):
    counts = [row["methylarginine_protein_count"] for row in rows]
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}">', '<rect width="100%" height="100%" fill="white"/>', '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#1d1d1d}.title{font-size:25px;font-weight:700}.label{font-size:19px}.tick{font-size:14px}.legend-title{font-size:16px;font-weight:700}.legend{font-size:14px}</style>', f'<text x="{LEFT}" y="38" class="title">Disease ontology enrichment among methylarginine proteins</text>']
    for tick in (1.0,1.25,1.5,1.75,2.0,2.25,2.5):
        x=x_pos(tick); parts += [f'<line x1="{x:.1f}" y1="{TOP-12}" x2="{x:.1f}" y2="{HEIGHT-BOTTOM}" stroke="#e8e8e8"/>', f'<text x="{x:.1f}" y="{HEIGHT-BOTTOM+27}" text-anchor="middle" class="tick">{tick:g}</text>']
    parts.append(f'<line x1="{x_pos(1):.1f}" y1="{TOP-12}" x2="{x_pos(1):.1f}" y2="{HEIGHT-BOTTOM}" stroke="#555" stroke-width="1.5"/>')
    for row,y in zip(rows,y_positions(rows)):
        r=radius(row["methylarginine_protein_count"],counts); x=x_pos(row["odds_ratio"]); parts += [f'<text x="{LEFT-18}" y="{y+6:.1f}" text-anchor="end" class="label">{row["disease_family"]}</text>', f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}" fill="{rgb_hex(color(row["BH_adjusted_p"]))}" stroke="white" stroke-width="2"/>']
    parts.append(f'<text x="{LEFT+(WIDTH-LEFT-RIGHT)/2:.1f}" y="{HEIGHT-20}" text-anchor="middle" class="label">Odds ratio vs reviewed-human background</text>')
    lx=WIDTH-RIGHT+48; parts.append(f'<text x="{lx}" y="121" class="legend-title">BH-adjusted p</text>')
    for i,q in enumerate((1e-10,1e-25,1e-50)):
        y=151+i*36; parts += [f'<circle cx="{lx+8}" cy="{y}" r="8" fill="{rgb_hex(color(q))}"/>', f'<text x="{lx+27}" y="{y+5}" class="legend">1×10^{int(math.log10(q))}</text>']
    if SIZE_MODE != "equal":
        parts.append(f'<text x="{lx}" y="301" class="legend-title">Methylarginine proteins</text>')
        for i,count in enumerate((500,1500,3500)):
            y=340+i*53; r=radius(count,counts); parts += [f'<circle cx="{lx+18}" cy="{y}" r="{r:.1f}" fill="#416f96" stroke="white" stroke-width="2"/>', f'<text x="{lx+50}" y="{y+5}" class="legend">{count:,}</text>']
    parts.append('</svg>')
    PREFIX.with_suffix(".svg").write_text("\n".join(parts),encoding="utf-8")


def draw_pdf(rows):
    c=canvas.Canvas(str(PREFIX.with_suffix(".pdf")),pagesize=(WIDTH,HEIGHT)); counts=[r["methylarginine_protein_count"] for r in rows]
    c.setFillColor("#1d1d1d"); c.setFont("Helvetica-Bold",25); c.drawString(LEFT,HEIGHT-38,"Disease ontology enrichment among methylarginine proteins")
    for tick in (1.0,1.25,1.5,1.75,2.0,2.25,2.5):
        x=x_pos(tick); c.setStrokeColor("#e8e8e8"); c.line(x,BOTTOM,x,HEIGHT-TOP+12); c.setFont("Helvetica",14); c.setFillColor("#222222"); c.drawCentredString(x,BOTTOM-27,f"{tick:g}")
    c.setStrokeColor("#555555"); c.line(x_pos(1),BOTTOM,x_pos(1),HEIGHT-TOP+12)
    for row,y_top in zip(rows,y_positions(rows)):
        y=HEIGHT-y_top; c.setFont("Helvetica",19); c.setFillColor("#1d1d1d"); c.drawRightString(LEFT-18,y-6,row["disease_family"]); r=radius(row["methylarginine_protein_count"],counts); c.setFillColor(rgb_hex(color(row["BH_adjusted_p"]))); c.circle(x_pos(row["odds_ratio"]),y,r,fill=1,stroke=0)
    c.setFont("Helvetica",19); c.setFillColor("#1d1d1d"); c.drawCentredString(LEFT+(WIDTH-LEFT-RIGHT)/2,20,"Odds ratio vs reviewed-human background")
    lx=WIDTH-RIGHT+48; c.setFont("Helvetica-Bold",16); c.drawString(lx,HEIGHT-121,"BH-adjusted p")
    for i,q in enumerate((1e-10,1e-25,1e-50)):
        y=HEIGHT-(151+i*36); c.setFillColor(rgb_hex(color(q))); c.circle(lx+8,y,8,fill=1,stroke=0); c.setFillColor("#1d1d1d"); c.setFont("Helvetica",14); c.drawString(lx+27,y-5,f"1e{int(math.log10(q))}")
    if SIZE_MODE != "equal":
        c.setFont("Helvetica-Bold",16); c.drawString(lx,HEIGHT-301,"Methylarginine proteins")
        for i,count in enumerate((500,1500,3500)):
            y=HEIGHT-(340+i*53); r=radius(count,counts); c.setFillColor("#416f96"); c.circle(lx+18,y,r,fill=1,stroke=0); c.setFillColor("#1d1d1d"); c.setFont("Helvetica",14); c.drawString(lx+50,y-5,f"{count:,}")
    c.save()


def main():
    global PREFIX, SIZE_MODE
    rows=load_rows()
    for SIZE_MODE, suffix in (("equal", "audited_v3_equal_size"), ("compressed", "audited_v3_compressed_count")):
        PREFIX = OUT / f"all_methyl_proteins_disease_overview_{suffix}"
        write_table(rows); draw_svg(rows); draw_png(rows); draw_pdf(rows)
        print(f"wrote {PREFIX.name} (top={rows[0]['disease_family']})")


if __name__ == "__main__":
    main()

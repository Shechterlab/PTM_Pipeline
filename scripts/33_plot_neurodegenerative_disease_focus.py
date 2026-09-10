import csv
import html
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "PTM_results/functional_class_union/disease_enrichment/all_methyl_proteins_neurodegenerative_focus.tsv"
PREFIX = ROOT / "PTM_results/functional_class_union/disease_enrichment/all_methyl_proteins_als_ftd_neurodegeneration_compact"
DISPLAY = [
    ("Motor neuron atrophy", "Motor neuron atrophy"),
    ("Motor Neuron Disease", "Motor neuron disease"),
    ("Frontotemporal dementia", "Frontotemporal dementia (FTD)"),
    ("Tauopathies", "Tauopathies"),
    ("Amyotrophic Lateral Sclerosis, Familial", "Familial ALS"),
]


def q_label(value: float) -> str:
    if value < 0.001:
        coefficient, exponent = f"{value:.1e}".split("e")
        return f"qBH={float(coefficient):g}×10^{int(exponent)}"
    return f"qBH={value:.3g}"


def main() -> None:
    with INPUT.open(encoding="utf-8", newline="") as handle:
        source = {row["Description"]: row for row in csv.DictReader(handle, delimiter="\t")}
    rows = [(label, source[name]) for name, label in DISPLAY]
    width, height = 1050, 480
    left, right, top, bottom = 335, 275, 72, 68
    plot_width = width - left - right
    x_min, x_max = 0.9, 2.05
    row_gap = 67

    def x(value: float) -> float:
        return left + (value - x_min) / (x_max - x_min) * plot_width

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#1d1d1d}.label{font-size:20px}.value{font-size:17px}.tick{font-size:15px}.note{font-size:14px;fill:#555}.title{font-size:24px;font-weight:700}</style>',
        f'<text x="{left}" y="34" class="title">Neurodegenerative disease ontology enrichment</text>',
    ]
    for tick in (1.0, 1.25, 1.5, 1.75, 2.0):
        xt = x(tick)
        parts.append(f'<line x1="{xt:.1f}" y1="{top-12}" x2="{xt:.1f}" y2="{height-bottom}" stroke="#e7e7e7" stroke-width="1"/>')
        parts.append(f'<text x="{xt:.1f}" y="{height-bottom+25}" text-anchor="middle" class="tick">{tick:g}</text>')
    parts.append(f'<line x1="{x(1):.1f}" y1="{top-12}" x2="{x(1):.1f}" y2="{height-bottom}" stroke="#4b4b4b" stroke-width="1.5"/>')
    for index, (label, row) in enumerate(rows):
        y = top + index * row_gap
        fold = float(row["fold_enrichment"])
        q = float(row["p.adjust"])
        count = int(row["Count"])
        parts.extend([
            f'<text x="{left-18}" y="{y+6}" text-anchor="end" class="label">{html.escape(label)}</text>',
            f'<line x1="{x(1):.1f}" y1="{y}" x2="{x(fold):.1f}" y2="{y}" stroke="#d8c5bc" stroke-width="7" stroke-linecap="round"/>',
            f'<circle cx="{x(fold):.1f}" cy="{y}" r="10" fill="#a85b45" stroke="white" stroke-width="2"/>',
            f'<text x="{x(fold)+18:.1f}" y="{y+6}" class="value">{html.escape(q_label(q))}; n={count}</text>',
        ])
    parts.extend([
        f'<text x="{left+plot_width/2:.1f}" y="{height-20}" text-anchor="middle" class="label">Fold enrichment vs reviewed-human background</text>',
        f'<text x="{width-12}" y="{height-8}" text-anchor="end" class="note">Two-sided enrichment test; BH adjustment across ontology terms</text>',
        '</svg>',
    ])
    PREFIX.with_suffix(".svg").write_text("\n".join(parts), encoding="utf-8")
    fields = ["Description", "display_label", "Count", "GeneRatio", "BgRatio", "fold_enrichment", "pvalue", "p.adjust", "qvalue"]
    with PREFIX.with_suffix(".tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for (name, label), (_, row) in zip(DISPLAY, rows):
            writer.writerow({**{key: row[key] for key in fields if key not in {"Description", "display_label"}}, "Description": name, "display_label": label})


if __name__ == "__main__":
    main()

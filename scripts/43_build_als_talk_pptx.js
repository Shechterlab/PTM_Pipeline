#!/usr/bin/env node
// Editable PowerPoint version of the ALS talk slides (script 41).
// Input : PTM_results/disease_gene_tests/als_talk_figures/als_talk_slide_data.json (script 42)
// Output: PTM_results/disease_gene_tests/als_talk_figures/methylarginine_ALS_talk.pptx
// Charts are native PowerPoint charts (right-click > Edit Data); protein maps are plain shapes.
// Requires pptxgenjs (npm). Usage: node scripts/43_build_als_talk_pptx.js [path/to/apply_theme.js]

const path = require("path");
const fs = require("fs");
const pptxgen = require("pptxgenjs");

const ROOT = path.resolve(__dirname, "..");
const RESULTS = process.env.PTM_RESULTS_ROOT || path.join(ROOT, "PTM_results");
const DIR = path.join(RESULTS, "disease_gene_tests", "als_talk_figures");
const D = JSON.parse(fs.readFileSync(path.join(DIR, "als_talk_slide_data.json"), "utf8"));
const OUT = path.join(DIR, "methylarginine_ALS_talk.pptx");

const HEX = {
  ink: "1F2937", navy: "1F3F7A", white: "FFFFFF", paper: "F3F4F6",
  red: "B5332A", slate: "687388", light: "C9CED8", domain: "3F5E8C", lcd: "EFE6D2", muted: "6B7280",
  grid: "E5E7EB", backbone: "B8BEC9",
};
const THEME = {
  name: "Methylarginine ALS",
  headFontFace: "Arial",
  bodyFontFace: "Arial",
  colors: {
    dk1: HEX.ink, lt1: HEX.white, dk2: HEX.navy, lt2: HEX.paper,
    accent1: HEX.red, accent2: HEX.slate, accent3: HEX.light, accent4: HEX.domain, accent5: HEX.lcd, accent6: HEX.muted,
    hlink: HEX.domain, folHlink: HEX.slate,
  },
};

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE"; // 13.333 x 7.5 in
pres.title = "Arginine methylation marks the ALS proteome";
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
const C = pres.SchemeColor;

const fmtE = (x) => { const [m, e] = x.toExponential(0).split("e"); return `${m}e${e.replace("+", "")}`; };

// ---------- layouts ----------
pres.defineSlideMaster({
  title: "TITLE",
  background: { color: HEX.navy },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.8, y: 2.3, w: 11.7, h: 1.6, fontSize: 40, bold: true,
        color: C.background1, valign: "bottom", align: "left", margin: 0 }, text: "" } },
    { placeholder: { options: { name: "body", type: "body", x: 0.8, y: 4.1, w: 11.7, h: 1.0, fontSize: 20,
        color: C.accent3, valign: "top", margin: 0 }, text: "" } },
  ],
});
pres.defineSlideMaster({
  title: "CONTENT",
  background: { color: HEX.white },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.6, y: 0.35, w: 12.1, h: 0.75, fontSize: 28, bold: true,
        color: C.text2, valign: "middle", align: "left", margin: 0 }, text: "" } },
    { placeholder: { options: { name: "body", type: "body", x: 0.6, y: 1.1, w: 12.1, h: 0.45, fontSize: 16,
        color: C.accent6, valign: "top", align: "left", margin: 0 }, text: "" } },
    { placeholder: { options: { name: "source", type: "body", x: 0.6, y: 6.75, w: 11.4, h: 0.5, fontSize: 10,
        italic: true, color: C.accent6, valign: "bottom", align: "left", margin: 0 }, text: "" } },
  ],
  slideNumber: { x: 12.3, y: 6.95, w: 0.5, h: 0.3, fontSize: 10, color: HEX.muted, align: "right" },
});

function content(section, title, subtitle, source, notes) {
  const s = pres.addSlide({ masterName: "CONTENT", sectionTitle: section });
  s.addText(title, { placeholder: "title" });
  s.addText(subtitle, { placeholder: "body" });
  s.addText(source, { placeholder: "source" });
  if (notes) s.addNotes(notes);
  return s;
}

function sidePanel(slide, x, heading, genes, more, stat) {
  slide.addText(heading, { x, y: 1.85, w: 3.0, h: 0.8, fontSize: 18, bold: true, color: C.accent1, margin: 0,
    valign: "top", isTextBox: true, objectName: "Side panel heading" });
  slide.addText([
    { text: genes.join("  "), options: { bold: true, color: C.text1, fontSize: 15, breakLine: true } },
    { text: more, options: { color: C.accent6, fontSize: 14 } },
  ], { x, y: 2.75, w: 3.0, h: 2.4, margin: 0, valign: "top", align: "left", isTextBox: true, paraSpaceBefore: 6,
    objectName: "Side panel proteins" });
  if (stat) slide.addText(stat, { x, y: 4.75, w: 3.0, h: 1.2, fontSize: 14, color: C.text1, margin: 0, valign: "top",
    isTextBox: true, objectName: "Side panel stat" });
}

function hbar(slide, rows, opts) {
  // Horizontal bar, top row first. Bar charts plot the first category at the bottom, so reverse.
  const r = [...rows].reverse();
  slide.addChart(pres.charts.BAR, [{ name: opts.series, labels: r.map((x) => x.label), values: r.map((x) => x.value) }], {
    x: opts.x, y: opts.y, w: opts.w, h: opts.h, barDir: "bar", barGapWidthPct: 45,
    chartColors: r.map((x) => (x.highlight ? HEX.red : HEX.slate)),
    catAxisLabelColor: HEX.ink, catAxisLabelFontSize: opts.catFont || 13, catAxisLabelFontFace: "+mn-lt",
    catAxisLineShow: false, catAxisMajorTickMark: "none",
    valAxisLabelColor: HEX.muted, valAxisLabelFontSize: 11, valAxisLabelFontFace: "+mn-lt",
    valAxisTitle: opts.valTitle, showValAxisTitle: true, valAxisTitleColor: HEX.muted, valAxisTitleFontSize: 12,
    valAxisTitleFontFace: "+mn-lt", valAxisLineShow: false, valAxisMinVal: 0,
    valGridLine: { color: HEX.grid, size: 0.75 }, catGridLine: { style: "none" },
    showLegend: false, showValue: false, objectName: opts.name,
  });
}

// ---------- slide 1: title ----------
pres.addSection({ title: "Title" });
const t = pres.addSlide({ masterName: "TITLE", sectionTitle: "Title" });
t.addText("Arginine methylation marks the ALS proteome", { placeholder: "title" });
t.addText("FET proteins, TDP-43, hnRNPs and most ALS-causing proteins carry methylarginine", { placeholder: "body" });

pres.addSection({ title: "ALS enrichment" });

// ---------- slide 2: ClinGen ----------
{
  const cg = D.clingen;
  const s = content("ALS enrichment", `${Math.round(cg.bars[0].pct)}% of ALS-causing proteins are arginine-methylated`,
    `The most of any disease group: ALS ranks 1st of ${cg.panels.length} ClinGen disease panels`,
    `ClinGen gene–disease validity (Definitive/Strong/Moderate; panels with ≥10 proteins). Methylated = any arginine-methylation site in the PTM_pipeline compendium. Adjusted for RNA binding, length and abundance: ALS odds ratio ${cg.als_adj_or[0].toFixed(1)} (95% CI ${cg.als_adj_or[1].toFixed(1)}–${cg.als_adj_or[2].toFixed(1)}).`,
    `Methylated ALS-causing proteins (${cg.bars[0].n}): ${cg.als_methyl_genes.join(", ")}.\n` +
    `Rank 1 of ${cg.panels.length} ClinGen panels by fraction methylated. Adjusted odds ratio ${cg.als_adj_or[0]} (95% CI ${cg.als_adj_or[1]}–${cg.als_adj_or[2]}), p = ${cg.als_adj_p.toFixed(3)}; this p-value does not survive correction across 42 panels, so quote the rank.\n` +
    `Not enriched for heavy methylation (>=3 sites) after adjustment.`);
  s.addChart(pres.charts.BAR, [{ name: "% methylated", labels: cg.bars.map((b) => `${b.label}\n(${b.n})`), values: cg.bars.map((b) => b.pct) }], {
    x: 0.6, y: 1.75, w: 6.6, h: 4.85, barDir: "col", barGapWidthPct: 55,
    chartColors: [HEX.red, HEX.slate, HEX.light],
    showValue: true, dataLabelPosition: "outEnd", dataLabelFormatCode: '0"%"', dataLabelFontSize: 28, dataLabelFontBold: true,
    dataLabelColor: HEX.ink, dataLabelFontFace: "+mn-lt",
    catAxisLabelColor: HEX.ink, catAxisLabelFontSize: 13, catAxisLabelFontFace: "+mn-lt", catAxisMajorTickMark: "none",
    valAxisHidden: true, valAxisMinVal: 0, valAxisMaxVal: 85, valGridLine: { style: "none" }, catGridLine: { style: "none" },
    showLegend: false, objectName: "ALS vs disease vs proteome",
  });
  const panels = [...cg.panels].reverse();
  s.addChart(pres.charts.BAR, [{ name: "% methylated", labels: panels.map((p) => p.panel), values: panels.map((p) => p.pct) }], {
    x: 7.7, y: 2.05, w: 5.0, h: 4.55, barDir: "bar", barGapWidthPct: 25,
    chartColors: panels.map((p) => (p.als ? HEX.red : HEX.light)),
    catAxisHidden: true, catAxisLineShow: false,
    valAxisMinVal: 0, valAxisMaxVal: 100, valAxisMajorUnit: 20, valAxisLabelColor: HEX.muted, valAxisLabelFontSize: 11,
    valAxisLabelFontFace: "+mn-lt", valAxisTitle: "% of panel proteins methylated", showValAxisTitle: true,
    valAxisTitleColor: HEX.muted, valAxisTitleFontSize: 12, valAxisTitleFontFace: "+mn-lt", valAxisLineShow: false,
    valGridLine: { color: HEX.grid, size: 0.75 }, catGridLine: { style: "none" }, showLegend: false,
    objectName: "All ClinGen panels",
  });
  s.addText([{ text: `All ${cg.panels.length} ClinGen disease panels: `, options: { color: C.text1 } }, { text: "ALS in red, at the top", options: { color: C.accent1 } }], { x: 7.7, y: 1.7, w: 5.0, h: 0.35, fontSize: 14, bold: true, margin: 0, isTextBox: true, objectName: "Panel chart heading" });
}

// ---------- slide 3: KEGG ----------
{
  const k = D.kegg;
  const s = content("ALS enrichment", `ALS is the #${k.als.rank} KEGG pathway for methylated proteins`,
    `${k.n_query.toLocaleString("en-US")} KEGG-annotated methylated proteins; top 12 of ${k.n_significant} significant pathways (n = methylated proteins)`,
    "clusterProfiler enrichKEGG (organism hsa, default KEGG background, BH), as in the 2022 analysis of Maron et al. 2021 Table S5, re-run on the PTM_pipeline methylarginine compendium.",
    `ALS pathway: ${k.als.count} of ${k.als.size} pathway proteins methylated, ${k.als.fold}-fold, FDR ${fmtE(k.als.fdr)}.\n` +
    "Caveats: whole-genome KEGG background, as in 2022. Against an MS-detectable background the ALS pathway ranks 11 (FDR 0.001). " +
    "The KEGG ALS map (hsa05014) also includes proteasome, mitochondrial and cytoskeletal proteins.");
  hbar(s, k.top_by_p.map((r) => ({ label: `${r.name}  (${r.count})`, value: r.neglog10fdr, highlight: r.highlight })),
    { x: 0.6, y: 1.75, w: 8.6, h: 4.9, series: "-log10 FDR", valTitle: "Enrichment, −log10 FDR", name: "KEGG top pathways" });
  sidePanel(s, 9.7, "Methylated proteins in the KEGG ALS pathway", k.showcase, `+ ${k.als.count - k.showcase.length} more`,
    `${k.als.count} of ${k.als.size} pathway proteins are methylated (${k.als.fold.toFixed(1)}-fold over expected; FDR ${fmtE(k.als.fdr)}).`);
}

// ---------- slide 4: protein maps ----------
{
  const P = D.proteins;
  const s = content("ALS enrichment", "Methylarginine sits in the RGG regions of ALS proteins",
    "Arginine-methylation sites on FET proteins, TDP-43, hnRNPA1/A2B1 and MATR3",
    "Sites: PTM_pipeline methylarginine compendium, canonical UniProt positions. Domains and Gly-rich / low-complexity regions: UniProt. TDP-43 R293 (C-terminal Gly-rich region) reported monomethylated in 7 proteomic studies; R275 in 3.",
    P.map((p) => `${p.label} (${p.uniprot}): ${p.sites.length} sites at R${p.sites.join(", R")}`).join("\n") +
    "\nSite counts use canonical UniProt positions only; rows mapped to isoform/TrEMBL accessions are excluded.");
  const x0 = 2.1, x1 = 11.7, maxLen = 860, sc = (x1 - x0) / maxLen;
  const X = (res) => x0 + res * sc;
  const top = 2.2, gap = 0.58;
  // legend
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 8.0, y: 1.62, w: 0.3, h: 0.2, rectRadius: 0.04, fill: { color: C.accent4 }, line: { type: "none" }, objectName: "Legend domain" });
  s.addText("Folded domain (RRM, ZnF)", { x: 8.38, y: 1.57, w: 2.3, h: 0.3, fontSize: 12, color: C.text1, margin: 0, isTextBox: true });
  s.addShape(pres.shapes.RECTANGLE, { x: 10.7, y: 1.62, w: 0.3, h: 0.2, fill: { color: C.accent5 }, line: { type: "none" }, objectName: "Legend LCD" });
  s.addText("Gly-rich / RGG / LCD", { x: 11.08, y: 1.57, w: 1.7, h: 0.3, fontSize: 12, color: C.text1, margin: 0, isTextBox: true });
  P.forEach((p, i) => {
    const y = top + i * gap;
    s.addText(p.label, { x: 0.5, y: y - 0.2, w: 1.45, h: 0.4, fontSize: 16, bold: true, color: C.text1, align: "right",
      margin: 0, isTextBox: true, objectName: `${p.label} name` });
    s.addShape(pres.shapes.RECTANGLE, { x: X(0), y: y - 0.04, w: p.length * sc, h: 0.08, fill: { color: HEX.backbone },
      line: { type: "none" }, objectName: `${p.label} backbone` });
    p.features.filter((f) => f.kind === "lcd").forEach((f) => s.addShape(pres.shapes.RECTANGLE, {
      x: X(f.start), y: y - 0.12, w: Math.max((f.end - f.start) * sc, 0.03), h: 0.24, fill: { color: C.accent5 },
      line: { type: "none" }, objectName: `${p.label} ${f.label} ${f.start}-${f.end}` }));
    p.features.filter((f) => f.kind === "domain").forEach((f) => s.addText(f.label, {
      shape: pres.shapes.ROUNDED_RECTANGLE, x: X(f.start), y: y - 0.15, w: (f.end - f.start) * sc, h: 0.3, rectRadius: 0.06,
      fill: { color: C.accent4 }, line: { type: "none" }, fontSize: 10, bold: true, color: C.background1, align: "center",
      valign: "middle", margin: 0, objectName: `${p.label} ${f.label} ${f.start}-${f.end}` }));
    p.sites.forEach((r) => {
      s.addShape(pres.shapes.LINE, { x: X(r), y: y - 0.3, w: 0, h: 0.14, line: { color: C.accent1, width: 1 },
        objectName: `${p.label} R${r} stem` });
      s.addShape(pres.shapes.OVAL, { x: X(r) - 0.03, y: y - 0.34, w: 0.06, h: 0.06, fill: { color: C.accent1 },
        line: { type: "none" }, objectName: `${p.label} R${r}` });
    });
    const countText = p.gene === "TARDBP" ? `${p.sites.length} sites (R${p.sites.join(", R")})` : `${p.sites.length} sites`;
    s.addText(countText, { x: X(p.length) + 0.12, y: y - 0.16, w: 2.6, h: 0.32, fontSize: 14, bold: true,
      color: C.accent1, margin: 0, isTextBox: true, objectName: `${p.label} site count` });
  });
  // residue axis
  const ay = top + P.length * gap - 0.15;
  s.addShape(pres.shapes.LINE, { x: X(0), y: ay, w: X(800) - X(0), h: 0, line: { color: HEX.muted, width: 0.75 }, objectName: "Residue axis" });
  [0, 200, 400, 600, 800].forEach((v) => {
    s.addShape(pres.shapes.LINE, { x: X(v), y: ay, w: 0, h: 0.06, line: { color: HEX.muted, width: 0.75 } });
    s.addText(String(v), { x: X(v) - 0.3, y: ay + 0.07, w: 0.6, h: 0.25, fontSize: 11, color: C.accent6, align: "center", margin: 0, isTextBox: true });
  });
  s.addText("Residue", { x: X(400) - 0.6, y: ay + 0.3, w: 1.2, h: 0.25, fontSize: 12, color: C.accent6, align: "center", margin: 0, isTextBox: true });
}

// ---------- slide 5: Disease Ontology ----------
{
  const d = D.do;
  const s = content("ALS enrichment", `Motor neuron disease ranks #${d.mnd.rank} of ${d.n_terms.toLocaleString("en-US")} disease terms`,
    `Disease Ontology enrichment of ${d.n_query.toLocaleString("en-US")} DO-annotated arginine-methylated proteins. Top 10 terms, plus ALS`,
    `DOSE enrichDO, reviewed human proteome background (BH). Motor neuron disease FDR ${fmtE(d.mnd.fdr)}; against a stricter MS-detectable background it ranks ${d.mnd_ms_detectable.rank} (FDR ${d.mnd_ms_detectable.fdr.toFixed(2)}).`,
    "Top Disease Ontology terms are intellectual-disability and developmental terms; motor neuron disease is the top neurodegenerative term. " +
    `Motor neuron disease: ${d.mnd.count} methylated proteins, FDR ${fmtE(d.mnd.fdr)}. ALS: rank 20, FDR 0.014.`);
  hbar(s, d.rows.map((r) => ({ label: `${r.rank > 10 ? `(rank ${r.rank})  ` : ""}${r.name}${r.highlight ? `  (${r.count} of ${r.size})` : ""}`,
    value: r.neglog10fdr, highlight: r.highlight })),
    { x: 0.6, y: 1.75, w: 8.8, h: 4.9, series: "-log10 FDR", valTitle: "Enrichment, −log10 FDR", name: "Disease Ontology top terms", catFont: 12 });
  sidePanel(s, 9.9, "Methylated motor neuron disease proteins", d.showcase, `+ ${d.mnd.count - d.showcase.length} more`, null);
}

// ---------- backup: KEGG ordered by gene ratio (2022 dot-plot order) ----------
pres.addSection({ title: "Backup" });
{
  const k = D.kegg;
  const s = content("Backup", "2022 KEGG re-run: neurodegeneration and ALS lead",
    "Top 20 pathways by gene ratio, as in the 2022 dot plot (methylated proteins; FDR)",
    "clusterProfiler enrichKEGG(organism = \"hsa\", pvalueCutoff = 0.05), default KEGG background, on the PTM_pipeline methylarginine compendium. Gene ratio = methylated proteins in pathway / KEGG-annotated methylated proteins.",
    "Same call as the 2022 enrichKEGG analysis. With Maron 2021 Table S5 alone, Pathways of neurodegeneration (173) and ALS (162, FDR 7e-15) are also the top two by gene ratio.");
  hbar(s, k.top_by_gene_ratio.map((r) => ({ label: `${r.name}  (${r.count}; ${fmtE(r.fdr)})`, value: +(100 * r.gene_ratio).toFixed(2), highlight: r.highlight })),
    { x: 0.6, y: 1.75, w: 12.1, h: 4.9, series: "Gene ratio (%)", valTitle: "Gene ratio, % of methylated proteins", name: "KEGG by gene ratio", catFont: 11 });
}

(async () => {
  await pres.writeFile({ fileName: OUT });
  const themeJs = process.argv[2];
  if (themeJs) {
    const { applyTheme } = require(themeJs);
    await applyTheme(OUT, THEME);
  }
  console.log("wrote", OUT);
})();

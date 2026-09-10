#!/usr/bin/env python3
"""Audited, mutually exclusive disease-family enrichment and figure export."""

from __future__ import annotations

import argparse
import csv
import html
import math
import re
import sqlite3
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "PTM_results/functional_class_union/disease_enrichment"
VERSION = "audited_v1"

# First match wins. This precedence is the reproducible primary assignment rule.
FAMILY_RULES = [
    ("ALS / motor neuron / FTD", r"amyotrophic|\bals\b|motor neuron|frontotemporal|\bftd\b|bulbar palsy"),
    ("Alzheimer / tau / dementia", r"alzheimer|tauopath|\bdementia\b|cognitive disorder"),
    ("Parkinsonism", r"parkinson"),
    ("Ataxia / cerebellar degeneration", r"ataxia|spinocerebellar|friedreich|cerebellar"),
    ("Myopathy / muscle / cardiomyopathy", r"myopath|muscular dystrophy|cardiomyopath|rhabdomyolysis|myositis"),
    ("Cancer / neoplasia", r"cancer|carcinoma|sarcoma|glioma|leukemia|lymphoma|melanoma|tumou?r|adenoma|neoplas|myeloma"),
    ("Immune / inflammatory", r"autoimmune|inflamm|lupus|arthritis|colitis|crohn|asthma|dermatitis|psoriasis"),
    ("Developmental / congenital syndromes", r"congenital|dysmorph|cranio|palate|facies|syndrom|hypoplasia|epicanthus|forehead|midface|feeding difficulties"),
]

# Exact-name overrides are deliberately separate from label matching.
CURATED_OVERRIDES = {
    "Ataxia Telangiectasia": (
        "Excluded: genome-instability syndrome",
        "excluded_from_overview",
        "Curated override: excluded rather than counted as general ataxia; genome-instability family is outside the 8-family overview.",
    ),
    "INCLUSION BODY MYOPATHY WITH EARLY-ONSET PAGET DISEASE AND FRONTOTEMPORAL DEMENTIA": (
        "ALS / motor neuron / FTD",
        "included",
        "Curated override: FTD takes precedence over the myopathy label.",
    ),
    "Inclusion Body Myopathy with Early-Onset Paget Disease with or without Frontotemporal Dementia 1": (
        "ALS / motor neuron / FTD",
        "included",
        "Curated override: FTD takes precedence over the myopathy label.",
    ),
    "INCLUSION BODY MYOPATHY WITH EARLY-ONSET PAGET DISEASE WITH OR WITHOUT FRONTOTEMPORAL DEMENTIA 2": (
        "ALS / motor neuron / FTD",
        "included",
        "Curated override: FTD takes precedence over the myopathy label.",
    ),
    "INCLUSION BODY MYOPATHY WITH EARLY-ONSET PAGET DISEASE WITH OR WITHOUT FRONTOTEMPORAL DEMENTIA 3": (
        "ALS / motor neuron / FTD",
        "included",
        "Curated override: FTD takes precedence over the myopathy label.",
    ),
}

FOCUS_TERMS = [
    ("C4024896", "Motor neuron atrophy"),
    ("C0085084", "Motor neuron disease"),
    ("C0338451", "Frontotemporal dementia (FTD)"),
    ("C4551993", "Familial ALS"),
    ("C0949664", "Tauopathies"),
    ("C0751072", "Frontotemporal lobar degeneration"),
]

PALETTE = {"line": "#D8C5BC", "dot": "#A85B45", "text": "#1D1D1D", "grid": "#E7E7E7", "muted": "#666666"}


def args_parser():
    p = argparse.ArgumentParser()
    p.add_argument("--labels", type=Path, default=ROOT / "PTM_results/functional_class_union/protein_functional_labels.tsv")
    p.add_argument("--term-gene", type=Path, required=True, help="TSV: term_id, term_name, entrez_id from DOSE DGN_PATHID2EXTID/PATHID2NAME")
    p.add_argument("--org-sqlite", type=Path, required=True, help="org.Hs.eg.db SQLite file used for SYMBOL to ENTREZ mapping")
    p.add_argument("--term-results", type=Path, default=DEFAULT_OUT / "all_methyl_proteins_disease_enrichment.tsv")
    p.add_argument("--input-summary", type=Path, default=DEFAULT_OUT / "disease_input_summary.tsv")
    p.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    return p.parse_args()


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def write_tsv(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader(); w.writerows(rows)


def primary_gene(text):
    return (text or "").split(";", 1)[0].strip()


def symbol_entrez_map(sqlite_path, symbols):
    con = sqlite3.connect(sqlite_path)
    tables = {r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    if not {"genes", "gene_info"}.issubset(tables):
        raise RuntimeError("Unexpected org.Hs.eg.db schema")
    # genes(_id, gene_id); gene_info(_id, symbol, gene_name)
    result = defaultdict(set)
    batch = 800
    for i in range(0, len(symbols), batch):
        chunk = symbols[i:i+batch]
        q = "SELECT gi.symbol, g.gene_id FROM gene_info gi JOIN genes g USING (_id) WHERE gi.symbol IN (%s)" % ",".join("?" * len(chunk))
        for symbol, entrez in con.execute(q, chunk): result[symbol].add(str(entrez))
    con.close()
    return result


def classify(name):
    if name in CURATED_OVERRIDES:
        fam, status, reason = CURATED_OVERRIDES[name]
        return fam, status, "curated_override", reason, ""
    matches = [fam for fam, pat in FAMILY_RULES if re.search(pat, name, re.I)]
    if not matches:
        return "Other / unassigned", "unassigned", "no_match", "No overview-family label rule matched.", ""
    return matches[0], "included", "precedence_rule", f"Primary assignment is the first matching family in precedence order ({matches[0]}).", " | ".join(matches[1:])


def logchoose(n, k):
    if k < 0 or k > n: return float("-inf")
    return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)


def fisher_two_sided(a, b, c, d):
    n1, n2, m1 = a+b, c+d, a+c
    lo, hi = max(0, m1-n2), min(n1, m1)
    def lp(x): return logchoose(m1, x)+logchoose(n1+n2-m1, n1-x)-logchoose(n1+n2, n1)
    observed = lp(a)
    vals = [lp(x) for x in range(lo, hi+1) if lp(x) <= observed + 1e-12]
    mx = max(vals)
    return min(1.0, math.exp(mx) * sum(math.exp(v-mx) for v in vals))


def bh(values):
    n=len(values); order=sorted(range(n), key=values.__getitem__); out=[0.0]*n; running=1.0
    for rank0 in range(n-1, -1, -1):
        idx=order[rank0]; running=min(running, values[idx]*n/(rank0+1)); out[idx]=running
    return out


def odds(a,b,c,d):
    if b*c == 0: return ((a+0.5)*(d+0.5))/((b+0.5)*(c+0.5))
    return (a*d)/(b*c)


def qfmt(x):
    return f"{x:.1e}" if x < 0.001 else f"{x:.3g}"


def svg_plot(rows, path, title, label_key, count_key, q_key, width=1280, height=None):
    height = height or (150 + 72*len(rows))
    left, right, top, bottom = 420, 300, 82, 72
    pw = width-left-right; vals=[float(r["odds_ratio"] if "odds_ratio" in r else r["fold_enrichment"]) for r in rows]
    xmin=.85; xmax=max(2.05, max(vals)*1.12)
    xp=lambda v: left+(v-xmin)/(xmax-xmin)*pw
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">', '<rect width="100%" height="100%" fill="white"/>', '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#1d1d1d}.title{font-size:25px;font-weight:700}.label{font-size:19px}.value{font-size:16px}.tick{font-size:14px}.note{font-size:14px;fill:#666}</style>', f'<text x="{left}" y="38" class="title">{html.escape(title)}</text>']
    ticks=[1,1.25,1.5,1.75,2] + ([2.25] if xmax>2.15 else [])
    for t in ticks:
        if t>xmax: continue
        xx=xp(t); parts += [f'<line x1="{xx:.1f}" y1="{top-18}" x2="{xx:.1f}" y2="{height-bottom}" stroke="{PALETTE["grid"]}"/>', f'<text x="{xx:.1f}" y="{height-bottom+25}" text-anchor="middle" class="tick">{t:g}</text>']
    parts.append(f'<line x1="{xp(1):.1f}" y1="{top-18}" x2="{xp(1):.1f}" y2="{height-bottom}" stroke="#444" stroke-width="1.5"/>')
    gap=(height-top-bottom)/max(1,len(rows))
    for i,r in enumerate(rows):
        y=top+(i+.45)*gap; v=float(r["odds_ratio"] if "odds_ratio" in r else r["fold_enrichment"]); q=float(r[q_key]); count=r[count_key]
        parts += [f'<text x="{left-18}" y="{y+6:.1f}" text-anchor="end" class="label">{html.escape(r[label_key])}</text>', f'<line x1="{xp(1):.1f}" y1="{y:.1f}" x2="{xp(v):.1f}" y2="{y:.1f}" stroke="{PALETTE["line"]}" stroke-width="7" stroke-linecap="round"/>', f'<circle cx="{xp(v):.1f}" cy="{y:.1f}" r="10" fill="{PALETTE["dot"]}" stroke="white" stroke-width="2"/>', f'<text x="{xp(v)+18:.1f}" y="{y+6:.1f}" class="value">BH-adjusted p={qfmt(q)}; n={count}</text>']
    axis="Odds ratio vs reviewed-human background" if "odds_ratio" in rows[0] else "Fold enrichment vs reviewed-human background"
    parts += [f'<text x="{left+pw/2:.1f}" y="{height-22}" text-anchor="middle" class="label">{axis}</text>', '</svg>']
    path.write_text("\n".join(parts), encoding="utf-8")


def raster_plot(rows, path, title, label_key, count_key, q_key, width=1280, height=None):
    height=height or (150+72*len(rows)); scale=1
    im=Image.new("RGB",(width,height),"white"); d=ImageDraw.Draw(im)
    font=lambda n,b=False: ImageFont.truetype("arialbd.ttf" if b else "arial.ttf",n)
    left,right,top,bottom=420,300,82,72; pw=width-left-right
    vals=[float(r.get("odds_ratio",r.get("fold_enrichment"))) for r in rows]; xmin=.85; xmax=max(2.05,max(vals)*1.12); xp=lambda v:left+(v-xmin)/(xmax-xmin)*pw
    d.text((left,18),title,font=font(25,True),fill=PALETTE["text"])
    for t in [1,1.25,1.5,1.75,2,2.25]:
        if t>xmax: continue
        x=xp(t); d.line((x,top-18,x,height-bottom),fill=PALETTE["grid"],width=1); txt=f"{t:g}"; d.text((x-d.textlength(txt,font=font(14))/2,height-bottom+8),txt,font=font(14),fill=PALETTE["text"])
    d.line((xp(1),top-18,xp(1),height-bottom),fill="#444444",width=2); gap=(height-top-bottom)/len(rows)
    for i,r in enumerate(rows):
        y=top+(i+.45)*gap; v=float(r.get("odds_ratio",r.get("fold_enrichment"))); label=r[label_key]; lx=left-18-d.textlength(label,font=font(19)); d.text((lx,y-11),label,font=font(19),fill=PALETTE["text"]); d.line((xp(1),y,xp(v),y),fill=PALETTE["line"],width=7); d.ellipse((xp(v)-10,y-10,xp(v)+10,y+10),fill=PALETTE["dot"],outline="white",width=2); ann=f"BH-adjusted p={qfmt(float(r[q_key]))}; n={r[count_key]}"; d.text((xp(v)+18,y-10),ann,font=font(16),fill=PALETTE["text"])
    axis="Odds ratio vs reviewed-human background" if "odds_ratio" in rows[0] else "Fold enrichment vs reviewed-human background"; d.text((left+pw/2-d.textlength(axis,font=font(19))/2,height-30),axis,font=font(19),fill=PALETTE["text"])
    im.save(path,dpi=(300,300))


def pdf_plot(rows, path, title, label_key, count_key, q_key, width=1280, height=None):
    height=height or (150+72*len(rows)); c=canvas.Canvas(str(path),pagesize=(width,height)); left,right,top,bottom=420,300,82,72; pw=width-left-right
    vals=[float(r.get("odds_ratio",r.get("fold_enrichment"))) for r in rows]; xmin=.85; xmax=max(2.05,max(vals)*1.12); xp=lambda v:left+(v-xmin)/(xmax-xmin)*pw
    c.setFillColor(PALETTE["text"]); c.setFont("Helvetica-Bold",25); c.drawString(left,height-38,title)
    for t in [1,1.25,1.5,1.75,2,2.25]:
        if t>xmax: continue
        x=xp(t); c.setStrokeColor(PALETTE["grid"]); c.line(x,bottom,x,height-top+18); c.setFont("Helvetica",14); c.setFillColor(PALETTE["text"]); c.drawCentredString(x,bottom-25,f"{t:g}")
    c.setStrokeColor("#444444"); c.line(xp(1),bottom,xp(1),height-top+18); gap=(height-top-bottom)/len(rows)
    for i,r in enumerate(rows):
        y=height-(top+(i+.45)*gap); v=float(r.get("odds_ratio",r.get("fold_enrichment"))); c.setFont("Helvetica",19); c.setFillColor(PALETTE["text"]); c.drawRightString(left-18,y-6,r[label_key]); c.setStrokeColor(PALETTE["line"]); c.setLineWidth(7); c.line(xp(1),y,xp(v),y); c.setFillColor(PALETTE["dot"]); c.circle(xp(v),y,10,fill=1,stroke=0); c.setFont("Helvetica",16); c.setFillColor(PALETTE["text"]); c.drawString(xp(v)+18,y-6,f"BH-adjusted p={qfmt(float(r[q_key]))}; n={r[count_key]}")
    axis="Odds ratio vs reviewed-human background" if "odds_ratio" in rows[0] else "Fold enrichment vs reviewed-human background"; c.setFont("Helvetica",19); c.drawCentredString(left+pw/2,22,axis); c.save()


def main():
    a=args_parser(); a.outdir.mkdir(parents=True,exist_ok=True); prefix=a.outdir/f"all_methyl_proteins_disease_overview_{VERSION}"
    labels=read_tsv(a.labels)
    symbols=[]; target_symbols=set()
    for r in labels:
        s=primary_gene(r.get("substrate_genename")) or primary_gene(r.get("genename"));
        if not s: continue
        symbols.append(s)
        if str(r.get("is_arg_methyl_protein","")).strip().lower() in {"1","true","t","yes","y"}: target_symbols.add(s)
    symbols=sorted(set(symbols)); smap=symbol_entrez_map(a.org_sqlite,symbols)
    universe={e for s in symbols for e in smap.get(s,set())}; target={e for s in target_symbols for e in smap.get(s,set())}; background=universe-target
    expected={r["target"]:int(r["entrez_count"]) for r in read_tsv(a.input_summary) if r.get("entrez_count")}
    if len(universe) != expected.get("Universe") or len(target) != expected.get("All methyl proteins"):
        raise RuntimeError(f"Identifier-map drift: observed universe/target {len(universe)}/{len(target)}, expected {expected.get('Universe')}/{expected.get('All methyl proteins')}")
    entrez_to_symbol={e:s for s in symbols for e in smap.get(s,set())}

    term_genes=defaultdict(set); term_names={}
    with a.term_gene.open(encoding="utf-8",newline="") as h:
        for r in csv.DictReader(h,delimiter="\t"):
            term_names[r["term_id"]]=r["term_name"]
            if r["entrez_id"] in universe: term_genes[r["term_id"]].add(r["entrez_id"])
    audit=[]; family_terms=defaultdict(list)
    for tid,name in sorted(term_names.items()):
        fam,status,method,reason,secondary=classify(name)
        audit.append({"term_id":tid,"term_name":name,"primary_family":fam,"audit_status":status,"assignment_method":method,"secondary_rule_matches":secondary,"override_reason":reason,"genes_in_reviewed_universe":len(term_genes.get(tid,set()))})
        if status=="included": family_terms[fam].append(tid)
    write_tsv(a.outdir/f"all_methyl_proteins_disease_term_family_audit_{VERSION}.tsv",audit,list(audit[0]))

    fam_rows=[]; fam_gene_sets={}
    for fam,_ in FAMILY_RULES:
        genes=set().union(*(term_genes[t] for t in family_terms[fam])) if family_terms[fam] else set(); fam_gene_sets[fam]=genes
        aa=len(target&genes); cc=len(background&genes); bb=len(target)-aa; dd=len(background)-cc
        fam_rows.append({"disease_family":fam,"contributing_ontology_terms":len(family_terms[fam]),"methylarginine_protein_count":aa,"background_count":cc,"family_gene_count":aa+cc,"methylarginine_denominator":len(target),"background_denominator":len(background),"odds_ratio":odds(aa,bb,cc,dd),"p_value":fisher_two_sided(aa,bb,cc,dd)})
    qs=bh([r["p_value"] for r in fam_rows])
    for r,q in zip(fam_rows,qs): r["BH_adjusted_p"]=q
    fam_rows.sort(key=lambda r:r["odds_ratio"])
    write_tsv(prefix.with_suffix(".tsv"),fam_rows,list(fam_rows[0]))

    # Pairwise gene overlap, including Jaccard and overlap coefficient.
    overlap=[]
    fams=[f for f,_ in FAMILY_RULES]
    for f1 in fams:
        for f2 in fams:
            inter=len(fam_gene_sets[f1]&fam_gene_sets[f2]); union=len(fam_gene_sets[f1]|fam_gene_sets[f2]); smaller=min(len(fam_gene_sets[f1]),len(fam_gene_sets[f2]))
            overlap.append({"family_1":f1,"family_2":f2,"family_1_genes":len(fam_gene_sets[f1]),"family_2_genes":len(fam_gene_sets[f2]),"overlap_genes":inter,"jaccard":inter/union if union else 0,"overlap_coefficient":inter/smaller if smaller else 0,"substantial_overlap_flag":inter/smaller>=0.50 if smaller else False})
    write_tsv(a.outdir/f"all_methyl_proteins_disease_family_gene_overlap_{VERSION}.tsv",overlap,list(overlap[0]))

    als="ALS / motor neuron / FTD"; als_terms=[]
    for tid in family_terms[als]:
        genes=term_genes[tid]; als_terms.append({"term_id":tid,"term_name":term_names[tid],"methylarginine_genes":len(target&genes),"background_genes":len(background&genes),"total_genes_in_universe":len(genes)})
    als_terms.sort(key=lambda r:(-r["methylarginine_genes"],r["term_name"]))
    write_tsv(a.outdir/f"all_methyl_proteins_als_ftd_contributing_terms_{VERSION}.tsv",als_terms,list(als_terms[0]))
    gene_rows=[]
    for e in sorted(fam_gene_sets[als],key=lambda x:(entrez_to_symbol.get(x,x),x)):
        tids=[t for t in family_terms[als] if e in term_genes[t]]
        gene_rows.append({"entrez_id":e,"gene_symbol":entrez_to_symbol.get(e,""),"is_methylarginine_protein":e in target,"contributing_term_count":len(tids),"contributing_term_ids":";".join(tids),"contributing_term_names":";".join(term_names[t] for t in tids)})
    write_tsv(a.outdir/f"all_methyl_proteins_als_ftd_contributing_genes_{VERSION}.tsv",gene_rows,list(gene_rows[0]))

    term_results={r["ID"]:r for r in read_tsv(a.term_results)}; focus=[]
    for tid,label in FOCUS_TERMS:
        if tid not in term_results: continue
        r=term_results[tid].copy(); r["source_term_id"]=tid; r["source_term_name"]=r["Description"]; r["display_label"]=label; r["BH_adjusted_p"]=r["p.adjust"]; focus.append(r)
    focus.sort(key=lambda r:float(r["fold_enrichment"]))
    focus_fields=["source_term_id","source_term_name","display_label","Count","GeneRatio","BgRatio","fold_enrichment","pvalue","p.adjust","BH_adjusted_p","qvalue","geneID"]
    focus_prefix=a.outdir/f"all_methyl_proteins_neurodegeneration_focus_{VERSION}"
    write_tsv(focus_prefix.with_suffix(".tsv"),focus,focus_fields)

    for ext,fun in [(".svg",svg_plot),(".png",raster_plot),(".pdf",pdf_plot)]:
        fun(fam_rows,prefix.with_suffix(ext),"Disease ontology enrichment among methylarginine proteins","disease_family","methylarginine_protein_count","BH_adjusted_p")
        fun(focus,focus_prefix.with_suffix(ext),"Neurodegenerative disease ontology enrichment","display_label","Count","BH_adjusted_p")

    old=next((r for r in read_tsv(a.outdir/"all_methyl_proteins_disease_theme_summary.tsv") if r["disease_theme"]==als),None); new=next(r for r in fam_rows if r["disease_family"]==als)
    max_pairs=sorted((r for r in overlap if r["family_1"]<r["family_2"]),key=lambda r:r["overlap_coefficient"],reverse=True)[:5]
    note=a.outdir/f"all_methyl_proteins_disease_ontology_analysis_{VERSION}.md"
    note.write_text(f"""# Audited disease-ontology consolidation\n\n## Methods and denominators\n\nThe analysis used the reviewed-human universe and identifier mapping of the existing pipeline. After SYMBOL-to-Entrez mapping, the methylarginine set contained **{len(target):,} genes** and the non-methylarginine background contained **{len(background):,} genes** ({len(universe):,} total). Each ontology term was assigned to at most one primary family using first-match precedence. Exact curated overrides were applied afterward and are visible in the audit table. Family gene sets are unions of unique reviewed-universe genes across assigned terms. Each gene is counted once per family. Two-sided Fisher exact tests were recalculated from each family union, followed by BH correction across the eight final family tests.\n\n## Consolidation decisions\n\nPrecedence is ALS/motor-neuron/FTD, Alzheimer/tau/dementia, Parkinsonism, ataxia/cerebellar degeneration, myopathy/muscle/cardiomyopathy, cancer/neoplasia, immune/inflammatory, then developmental/congenital. This prevents FTD labels from entering general dementia. The ALS abbreviation is matched as a standalone token rather than as the substring `als`; this removes false matches inside unrelated words. Ataxia telangiectasia is explicitly excluded from the plotted ataxia family because its defining placement here is genome instability; the audit retains that decision. Four FTD-plus-inclusion-body-myopathy labels are explicitly assigned to ALS/motor-neuron/FTD. Unmatched terms remain `Other / unassigned` and are not plotted. Term assignment never used enrichment magnitude or significance.\n\n## ALS / motor neuron / FTD comparison\n\nThe audited family contains **{new['methylarginine_protein_count']:,} methylarginine genes among {new['family_gene_count']:,} family genes**, OR **{new['odds_ratio']:.3f}**, two-sided p **{new['p_value']:.3g}**, BH-adjusted p **{new['BH_adjusted_p']:.3g}**. The prior regex result was {old['target_count']}/{old['total_theme_gene_count']} genes, OR {float(old['odds_ratio']):.3f}, BH q {float(old['q_value']):.3g}. The reduction is driven mainly by correcting the unbounded `als` substring rule, which had admitted unrelated labels, plus mutually exclusive precedence and the documented exact-name overrides. The tested universe is unchanged.\n\nThe contributing-term table is ranked by methylarginine-gene count. The contributing-gene table gives every family gene, methylarginine membership, and source terms.\n\n## Family overlap sensitivity\n\nDisease families are not biologically independent. The five largest overlap coefficients are:\n\n""" + "\n".join(f"- {r['family_1']} vs {r['family_2']}: {r['overlap_genes']} shared genes; Jaccard {r['jaccard']:.3f}; overlap coefficient {r['overlap_coefficient']:.3f}{' (flagged)' if r['substantial_overlap_flag'] else ''}." for r in max_pairs) + """\n\n## Plot-table agreement and interpretation\n\nBoth figures are generated directly from their exact exported TSV tables in the same script. Overview annotations use family-union counts and family-level BH-adjusted p-values. The focused panel retains source term identifiers and names and displays the existing clusterProfiler `p.adjust` values; Storey `qvalue` is retained only for provenance. These results show enrichment of disease-associated gene sets among methylarginine proteins. They do not establish that methylarginine causes disease or that the ontology families are independent mechanisms.\n\n## Recommended legend language\n\n**Disease ontology enrichment among methylarginine proteins.** Disease terms were assigned to mutually exclusive primary families by prespecified label precedence plus documented exact-name overrides. Genes were unioned within each family and counted once. Odds ratios compare mapped methylarginine proteins with the remaining reviewed-human background. P values are from two-sided Fisher exact tests, with Benjamini-Hochberg adjustment across the eight plotted families. The focused panel shows selected direct neurodegenerative terms and uses the original clusterProfiler BH-adjusted `p.adjust` values.\n""",encoding="utf-8")
    top_gene_rows=sorted((r for r in gene_rows if r["is_methylarginine_protein"]),key=lambda r:(-r["contributing_term_count"],r["gene_symbol"]))[:15]
    with note.open("a",encoding="utf-8") as h:
        h.write("\n## Top ALS / motor neuron / FTD contributors\n\n### Ontology terms\n\n")
        h.write("\n".join(f"- {r['term_name']} ({r['term_id']}): {r['methylarginine_genes']} methylarginine genes among {r['total_genes_in_universe']} reviewed-universe genes." for r in als_terms[:10]))
        h.write("\n\n### Methylarginine genes\n\n")
        h.write("\n".join(f"- {r['gene_symbol']} ({r['entrez_id']}): present in {r['contributing_term_count']} contributing terms." for r in top_gene_rows))
        h.write("\n")
    print(f"mapped_target={len(target)} background={len(background)} universe={len(universe)}")
    print(f"ALS={new['methylarginine_protein_count']}/{new['family_gene_count']} OR={new['odds_ratio']:.6g} BH={new['BH_adjusted_p']:.6g}")


if __name__ == "__main__": main()

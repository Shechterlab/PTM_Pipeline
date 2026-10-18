/* Methylarginine Site Explorer (Shechter Lab).
   Static, client-side app. Data: docs/data/*.json written by
   scripts/36_export_site_explorer_data.py from PTM_results/. Protein
   features and sequences are fetched live from the UniProt REST API. */
(function () {
  'use strict';

  // ---------------------------------------------------------------- helpers
  const $ = (sel, root = document) => root.querySelector(sel);
  const $$ = (sel, root = document) => Array.from(root.querySelectorAll(sel));
  const ESC = { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' };
  const esc = (s) => String(s == null ? '' : s).replace(/[&<>"']/g, (c) => ESC[c]);
  const fmt = (n) => (n == null ? '—' : Number(n).toLocaleString('en-US'));
  const compact = (n) => (n >= 1000 ? (n / 1000).toFixed(n >= 10000 ? 0 : 1).replace(/\.0$/, '') + 'K' : String(n));
  const fmtP = (v) => (v == null ? '—' : v === 0 ? '< 1e-300' : v < 0.001 ? v.toExponential(1) : v < 0.1 ? v.toPrecision(2) : v.toFixed(2));
  const fmtL2 = (v) => (v == null ? '—' : (v > 0 ? '+' : '') + v.toFixed(2));
  const pct = (a, b) => (b ? Math.round((100 * a) / b) : 0);
  const pctText = (a, b) => (a > 0 && pct(a, b) === 0 ? '<1' : String(pct(a, b)));
  const popcount = (n) => { let c = 0; while (n) { c += n & 1; n >>= 1; } return c; };
  const bits = (mask) => { const out = []; for (let i = 0; mask; i++, mask >>= 1) if (mask & 1) out.push(i); return out; };
  const debounce = (fn, ms) => { let t; return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); }; };
  const r1 = (v) => Math.round(v * 10) / 10;
  // Residue numbers as shown. Histones follow the convention without the initiator methionine (H3R2 is UniProt R3 of
  // P68431); every other protein, and every stored position, uses UniProt numbering with Met as 1. p defaults to the
  // protein page that is open.
  const HISTONE_GENE = /^(H1-\d+|H2A[A-Z0-9-]*|H2B[A-Z0-9-]*|H3C\d+|H3-\d+[A-Z]?|H3Y\d|H4C\d+|H4-\d+|MACROH2A\d|CENPA)$/;
  const HISTONE_NAME = /^(putative |core )?histone (h1|h2a|h2b|h3|h4|macro-h2a)\b|^testis-specific h1 histone|^histone h3-like centromeric/i;
  const rn = (pos, p) => { p = p || (UI.detail != null ? P[UI.detail] : null); return pos - (p && p.hist ? 1 : 0); };
  const rl = (pos, p) => 'R' + rn(pos, p);
  const store = {
    get(k) { try { return localStorage.getItem(k); } catch (e) { return null; } },
    set(k, v) { try { if (v == null) localStorage.removeItem(k); else localStorage.setItem(k, v); } catch (e) { /* storage unavailable */ } },
  };

  // -------------------------------------------------------------- vocabulary
  const STATE_LABEL = ['No state call', 'Rme1', 'Rme2', 'Rme1 + Rme2'];
  const SOURCE_LABEL = { Maron2021: 'Maron 2021 Table S5', ProMetheusDB: 'ProMetheusDB', human_ptm_master: 'PTM master table' };
  const SOURCE_HELP = {
    Maron2021: 'Maron et al. 2021 (iScience) Table S5: compiled human methylarginine sites from published proteomics',
    ProMetheusDB: 'ProMetheusDB methyl-proteome (mmc2), with Rme1/Rme2 calls',
    human_ptm_master: 'Cross-PTM human master table (UniProt and IEDB methylarginine annotations)',
  };
  const STUDY_HELP = {
    ProMetheusDB: 'ProMetheusDB methyl-proteome table (mmc2)',
    UniProt: 'UniProt methylarginine annotation, via the cross-PTM master table',
    IEDB: 'Immune Epitope Database record with a methylated arginine, via the cross-PTM master table',
    'Maron 2021': 'Sites first reported in Maron et al. 2021 (Table S5 "This study")',
  };
  // DOIs of the studies and databases behind the union (each checked against Crossref).
  const STUDY_DOI = {
    'Maron 2021': '10.1016/j.isci.2021.102971', ProMetheusDB: '10.1016/j.mcpro.2022.100243', PhosphoSitePlus: '10.1093/nar/gku1267',
    'Larsen 2016': '10.1126/scisignal.aaf7329', 'Musiani 2019': '10.1126/scisignal.aat8388', 'Fong 2019': '10.1016/j.ccell.2019.07.003',
    'Guo 2014': '10.1074/mcp.O113.027870', 'Fedoriw 2019': '10.1016/j.ccell.2019.05.014', 'Wei 2020': '10.3389/fcell.2020.545934',
    'Lim 2020': '10.1186/s12953-020-00162-8', 'Li 2021': '10.1038/s41467-021-21963-1', IEDB: '10.1093/nar/gky1006', UniProt: '10.1093/nar/gkae1010',
  };
  const doiLink = (doi, text) => '<a href="https://doi.org/' + doi + '" target="_blank" rel="noopener">' + (text == null ? 'doi:' + doi : text) + '</a>';
  const PSP_NOTICE = 'PhosphoSitePlus® data (www.phosphosite.org) are licensed under CC BY-NC-SA 3.0: non-commercial use, with attribution, and adaptations shared under the same licence.';
  // Data sources, versions and licences shown on the Documentation page. [name, used for, version, licence, citation, doi]
  const DATA_SOURCES = [
    ['Site union', [
      ['Maron et al. 2021, Table S5', 'Compiled methylarginine sites; Table S2/S3 state and PRMT inhibitor response', 'Supplementary tables', 'Journal article supplement', 'iScience 24:102971', '10.1016/j.isci.2021.102971'],
      ['PhosphoSitePlus®', 'Sites compiled in Maron Table S5', 'As compiled in Maron et al. 2021', 'CC BY-NC-SA 3.0 (non-commercial)', 'Hornbeck et al., Nucleic Acids Res 43:D512 (2015)', '10.1093/nar/gku1267'],
      ['ProMetheusDB', 'Methyl-proteome table mmc2 (Rme1/Rme2 calls)', 'mmc2, sheet RK-methylsites', 'Journal article supplement', 'Massignani et al., Mol Cell Proteomics 21:100243 (2022)', '10.1016/j.mcpro.2022.100243'],
      ['iPTMnet (PTM master table)', 'Other-PTM sites shown as crosstalk; methylarginine rows', 'Lab export', 'See iPTMnet terms', 'Huang et al., Nucleic Acids Res 46:D542 (2018)', '10.1093/nar/gkx1104'],
      ['IEDB (PTM master table)', 'Epitope records with a methylated arginine', 'Lab export', 'IEDB terms of use', 'Vita et al., Nucleic Acids Res 47:D339 (2019)', '10.1093/nar/gky1006'],
    ]],
    ['Annotations in the pipeline', [
      ['UniProt', 'Canonical sequences, remapping, gene symbols and synonyms', 'Reviewed human proteome at staging', 'CC BY 4.0', 'UniProt Consortium, Nucleic Acids Res 53:D609 (2025)', '10.1093/nar/gkae1010'],
      ['InterPro', 'Domain context and nearest domain; Tudor, WD40 and PHD classes of interaction partners', 'Current release at staging', 'CC0', 'Blum et al., Nucleic Acids Res 53:D444 (2025)', '10.1093/nar/gkae1082'],
      ['MobiDB', 'Disorder context (predicted consensus disorder, prediction-disorder-th_50)', 'Human export of 30 Mar 2026', 'Not stated on mobidb.org (CC BY 4.0 reported elsewhere)', 'Piovesan et al., Nucleic Acids Res 53:D495 (2025)', '10.1093/nar/gkae969'],
      ['CD-CODE', 'Condensate membership (filters and enrichment); per-protein confidence on protein pages', 'Export of 18 Mar 2026; API release 2.3', 'Not stated by CD-CODE', 'Rostam et al., Nat Methods 20:673 (2023)', '10.1038/s41592-023-01831-0'],
      ['DisGeNET (via DOSE)', 'Disease enrichment terms and families', 'As bundled with DOSE / clusterProfiler 4.18 (Bioconductor, R 4.5)', 'CC BY-NC-SA 4.0 (non-commercial)', 'Piñero et al., Nucleic Acids Res (2019 update)', '10.1093/nar/gkz1021'],
      ['clusterProfiler / DOSE', 'Over-representation analysis', 'clusterProfiler 4.18', 'Artistic-2.0 (software)', 'Wu et al., The Innovation 2:100141 (2021); Yu et al., Bioinformatics 31:608 (2015)', '10.1016/j.xinn.2021.100141'],
    ]],
    ['Interactions and structures (exported with scripts/37)', [
      ['STRING', 'Physical partners (score ≥ 0.7 exported; lower scores live)', 'v12.5', 'CC BY 4.0', 'Szklarczyk et al., Nucleic Acids Res 53:D730 (2025)', '10.1093/nar/gkae1113'],
      ['Complex Portal', 'Curated human complexes', 'Current release at export', 'CC0', 'Meldal et al., Nucleic Acids Res 50:D578 (2022)', '10.1093/nar/gkab991'],
      ['hu.MAP 3.0', 'Predicted complexes (via Complex Portal)', 'hu.MAP 3.0', 'CC0', 'Fischer et al., Mol Syst Biol 21 (2025)', '10.1038/s44320-025-00121-5'],
      ['Gene Ontology / QuickGO', 'RNA- and DNA-binding annotations', 'Current release at export', 'CC BY 4.0', 'Gene Ontology Consortium, Genetics 224:iyad031 (2023); Binns et al., Bioinformatics 25:3045 (2009)', '10.1093/genetics/iyad031'],
      ['PDBe-KB and SIFTS', 'Interface residues; PDB–UniProt residue mapping', 'Current release at export', 'EMBL-EBI terms of use', 'PDBe-KB consortium, Nucleic Acids Res 50:D534 (2022); Dana et al., Nucleic Acids Res 47:D482 (2019)', '10.1093/nar/gkab988'],
      ['RCSB PDB', 'Structure files, ModelServer contacts, search, entity metadata', 'Live / at export', 'CC0 (PDB data)', 'Burley et al., Nucleic Acids Res 51:D488 (2023)', '10.1093/nar/gkac1077'],
      ['ProtVar', 'Predicted interfaces (5 Å) and served AlphaFold-Multimer models', 'Interface table 2024.05.28; served-model list af2complexes_interfaces_v1', 'CC BY 4.0', 'Stephenson et al., Nucleic Acids Res 52:W140 (2024)', '10.1093/nar/gkae413'],
      ['AlphaFold DB complexes', 'Predicted heterodimers (gallery and protein pages)', 'Live API', 'CC BY 4.0', 'Han et al., preprint (2026): AlphaFold Database expands to proteome-scale quaternary structures', '10.64898/2026.03.27.714458'],
    ]],
    ['RNA and DNA contact in living cells (exported with scripts/37)', [
      ['pepR-MS compendium', 'RNA-crosslinked residues: pepR-MS and re-analysed Kramer 2014, XRNAX, RBS-ID, pRBS-ID (Tables S2–S3)', 'bioRxiv v2 (2026) supplement', 'CC BY-NC', 'Sha, Kuster & Trendel, bioRxiv (2025)', '10.1101/2025.09.05.674461'],
      ['Kramer et al. 2014', 'RNA crosslinks (via pepR-MS re-analysis)', '', 'Journal article', 'Nat Methods 11:1064 (2014)', '10.1038/nmeth.3092'],
      ['XRNAX', 'RNA crosslinks (via pepR-MS re-analysis)', '', 'Journal article', 'Trendel et al., Cell 176:391 (2019)', '10.1016/j.cell.2018.11.004'],
      ['RBS-ID', 'RNA crosslinks at single-residue resolution (via pepR-MS re-analysis)', '', 'Journal article', 'Bae et al., Nat Struct Mol Biol 27:678 (2020)', '10.1038/s41594-020-0436-2'],
      ['pRBS-ID', 'Base-specific RNA crosslinks (via pepR-MS re-analysis)', '', 'Journal article', 'Bae et al., Nat Commun 12:6026 (2021)', '10.1038/s41467-021-26317-5'],
      ['RBDmap', 'RNA-bound fragments (RBDpep, HeLa)', 'Table S1 (Oxford Research Archive)', 'ORA terms; article open access', 'Castello et al., Mol Cell 63:696 (2016)', '10.1016/j.molcel.2016.06.029'],
      ['DNA interactome', 'Proteins in direct contact with DNA; DNA-crosslinked peptides and residues (MCF7)', 'Tables S1–S2', 'Journal supplement (Elsevier)', 'Trendel et al., Cell 188:4424 (2025)', '10.1016/j.cell.2025.04.037'],
      ['UniProt DNA-binding regions', 'Curated DNA-binding regions (DNA_BIND)', 'At export', 'CC BY 4.0', 'UniProt Consortium, Nucleic Acids Res 53:D609 (2025)', '10.1093/nar/gkae1010'],
    ]],
    ['Live on protein pages', [
      ['AlphaFold DB', 'Predicted models and pLDDT', 'Live API', 'CC BY 4.0', 'Jumper et al., Nature 596:583 (2021); Varadi et al., Nucleic Acids Res 52:D368 (2024)', '10.1093/nar/gkad1011'],
      ['AlphaMissense', 'R→X pathogenicity at each site', 'Via AlphaFold DB', 'CC BY 4.0', 'Cheng et al., Science 381:eadg7492 (2023)', '10.1126/science.adg7492'],
      ['ENCODE', 'eCLIP and TF ChIP-seq experiments', 'Live API', 'Unrestricted use (ENCODE data use policy)', 'ENCODE Project Consortium, Nature 583:699 (2020)', '10.1038/s41586-020-2493-4'],
      ['ChIP-Atlas', 'TF target genes', 'Live', 'CC BY 4.0', 'Zou et al., Nucleic Acids Res 52:W45 (2024)', '10.1093/nar/gkae358'],
      ['Open Targets Platform', 'Gene–disease association scores', 'Live GraphQL API', 'CC0', 'Buniello et al., Nucleic Acids Res 53:D1467 (2025)', '10.1093/nar/gkae1128'],
      ['RCSB sequence search', 'Homologous structures (MMseqs2)', 'Live API', 'CC0 (PDB data)', 'Steinegger & Söding, Nat Biotechnol 35:1026 (2017)', '10.1038/nbt.3988'],
      ['NCBI E-utilities', 'Reference titles for CD-CODE PubMed IDs', 'Live API', 'NCBI policies', '', ''],
      ['Ensembl Compara', 'Orthologue alignments for site conservation', 'Live REST API (rest.ensembl.org)', 'No restrictions (Ensembl data policy)', 'Herrero et al., Database 2016:bav096; Yates et al., Bioinformatics 31:143 (2015)', '10.1093/database/bav096'],
    ]],
    ['Methods and software', [
      ['3Dmol.js 2.5.3', '3D viewer', '2.5.3', 'BSD-3-Clause', 'Rego & Koes, Bioinformatics 31:1322 (2015)', '10.1093/bioinformatics/btu829'],
      ['Sequence logos', 'Information content with small-sample correction; WebLogo chemistry colours', '', '', 'Schneider et al., J Mol Biol 188:415 (1986); Crooks et al., Genome Res 14:1188 (2004)', '10.1101/gr.849004'],
      ['Levy interface classes', 'Core, rim and support', '', '', 'Levy, J Mol Biol 403:660 (2010)', '10.1016/j.jmb.2010.09.028'],
      ['EMBOSS charge', 'Sliding 5-residue charge on protein maps (K/R +1, H +0.5, D/E −1)', 'EMBOSS charge scale (Eamino.dat)', 'GPL (software)', 'Rice, Longden & Bleasby, Trends Genet 16:276 (2000)', '10.1016/S0168-9525(00)02024-2'],
      ['Solvent accessibility', 'Shrake–Rupley; Arg maximum 274 Å²', '', '', 'Shrake & Rupley, J Mol Biol 79:351 (1973); Tien et al., PLoS ONE 8:e80635 (2013)', '10.1016/0022-2836(73)90011-9'],
      ['AlphaFold-Multimer and pDockQ', 'Predicted complexes and their interface score', '', '', 'Evans et al., bioRxiv (2021); Bryant et al., Nat Commun 13:1265 (2022)', '10.1038/s41467-022-28865-w'],
      ['pDockQ thresholds', 'pDockQ ≥ 0.23 (acceptable) and ≥ 0.5 (confident)', '', '', 'Burke et al., Nat Struct Mol Biol 30:216 (2023)', '10.1038/s41594-022-00910-8'],
    ]],
  ];

  const PAPER = {
    cite: 'DeAngelo JD, Silverstein AM, Shechter D. Protein arginine methylation rewires molecular interactions. Cell Chem Biol 33(10):1404–1421 (2026).',
    short: 'DeAngelo, Silverstein & Shechter, Cell Chem Biol 2026',
    doi: '10.1016/j.chembiol.2026.09.011',
  };
  // The site's name and address come from the page itself (<meta name="application-name"> in index.html and the
  // current URL), so renaming or moving the site updates every figure credit.
  const siteName = () => { const m = document.querySelector('meta[name="application-name"]'); return (m && m.content) || document.title.split(' · ')[0]; };
  const siteUrl = () => (location.protocol === 'file:' ? location.pathname : location.origin + location.pathname).replace(/index\.html$/, '');
  const creditLines = () => [siteName() + ' · ' + siteUrl().replace(/^https?:\/\//, ''), 'Please cite ' + PAPER.short + ', doi:' + PAPER.doi + ' · and the original data sources'];

  // ----- citing the original data. Every list of sites can say which studies reported them and which resources
  // supplied the evidence shown, so users cite the primary work as well as the compilation.
  const SOURCE_CITE = {
    'Maron 2021': ['Maron MI et al., iScience 24:102971 (2021)', '10.1016/j.isci.2021.102971'],
    ProMetheusDB: ['Massignani E et al., Mol Cell Proteomics 21:100243 (2022)', '10.1016/j.mcpro.2022.100243'],
    PhosphoSitePlus: ['Hornbeck PV et al., Nucleic Acids Res 43:D512 (2015)', '10.1093/nar/gku1267'],
    UniProt: ['UniProt Consortium, Nucleic Acids Res 53:D609 (2025)', '10.1093/nar/gkae1010'],
    IEDB: ['Vita R et al., Nucleic Acids Res 47:D339 (2019)', '10.1093/nar/gky1006'],
  };
  // Full reference of a study or database behind the union: [text, doi].
  function studyRef(k) {
    const name = M.studies[k][0], ref = M.studies[k][1];
    if (SOURCE_CITE[name]) return SOURCE_CITE[name];
    return [ref && ref !== 'This study.' ? ref.replace(/\s+/g, ' ').trim() : name, STUDY_DOI[name] || ''];
  }
  const dsCite = (name) => { for (const [, rows] of DATA_SOURCES) { const r = rows.find((x) => x[0] === name); if (r) return [r[4], r[5]]; } return [name, '']; };
  // The studies reporting these sites, most sites first, and the evidence resources their annotations come from.
  function citeRefs(sites) {
    const n = new Array(M.studies.length).fill(0);
    for (const s of sites) for (const b of bits(s.study)) n[b]++;
    const studies = n.map((c, k) => ({ k, name: M.studies[k][0], n: c, ref: studyRef(k) })).filter((x) => x.n).sort((a, b) => b.n - a.n || a.name.localeCompare(b.name));
    const ev = [];
    const use = (label, names, count) => { if (count) ev.push({ label, n: count, refs: names.map((x) => (Array.isArray(x) ? x : dsCite(x))) }); };
    const cnt = (f) => sites.reduce((a, s) => a + (f(s) ? 1 : 0), 0);
    use('PRMT inhibitor response and methyl-IP calls (Maron Tables S2–S3)', [SOURCE_CITE['Maron 2021']], cnt((s) => s.tr.length || s.ip || s.s2));
    if (ixReady()) {
      const sets = new Set();
      for (const s of sites) if (s.xl && (s.xl.r || s.xl.n)) xlSets((s.xl.r || 0) | (s.xl.n || 0)).forEach((x) => sets.add(x));
      const setNames = { 'Kramer 2014': 'Kramer et al. 2014', 'XRNAX 2019': 'XRNAX', 'RBS-ID 2020': 'RBS-ID', 'pRBS-ID 2021': 'pRBS-ID' };
      use('RNA crosslinks in cells', ['pepR-MS compendium'].concat(Array.from(sets).map((x) => setNames[x]).filter(Boolean)), cnt((s) => s.xl && (s.xl.r || s.xl.n)));
      use('RBDmap RNA-bound fragments', ['RBDmap'], cnt((s) => s.xl && s.xl.f));
      use('DNA crosslinks in cells', ['DNA interactome'], cnt((s) => s.xl && s.xl.d & 11));
      use('PDB interfaces and contacts', ['PDBe-KB and SIFTS', 'RCSB PDB'], cnt((s) => s.iface || s.contacts));
      use('Predicted complexes', ['ProtVar', 'AlphaFold-Multimer and pDockQ'], cnt((s) => s.pred || (s.gal || []).length));
    }
    return { studies, ev };
  }
  const refLine = ([text, doi]) => esc(text) + (doi ? ' ' + doiLink(doi) : '');
  const refText = ([text, doi]) => text + (doi ? ' doi:' + doi : '');
  // "Cite the sources": the compilation, the reporting studies (with how many of these sites each reported) and the
  // evidence resources. open: expanded (reports); what: the noun for the sites ("these 209 sites", "this site").
  function citeBoxHtml(sites, what, open) {
    if (!sites.length) return '';
    const { studies, ev } = citeRefs(sites);
    const text = ['Compilation and explorer: ' + PAPER.cite + ' doi:' + PAPER.doi, '', 'Studies and databases that reported ' + what + ':']
      .concat(studies.map((x) => '- ' + refText(x.ref) + ' [' + x.n + ' site' + (x.n === 1 ? '' : 's') + ']'))
      .concat(ev.length ? ['', 'Evidence shown for ' + what + ':'].concat(ev.map((e) => '- ' + e.label + ': ' + e.refs.map(refText).join('; '))) : []).join('\n');
    return '<details class="cite-box"' + (open ? ' open' : '') + '><summary>How to cite ' + esc(what) + '</summary>' +
      '<p class="small">Please cite the original studies that reported the sites you use, the resources behind any evidence you rely on, and ' + esc(PAPER.short) + ' (' + doiLink(PAPER.doi) + ') for the compilation and this explorer. Full list: ' + docLink('cite').replace('>Method<', '>How to cite<') + '</p>' +
      '<h4>Reported by</h4><ol class="cite-list">' + studies.map((x) => '<li>' + refLine(x.ref) + ' <span class="muted">· ' + fmt(x.n) + ' site' + (x.n === 1 ? '' : 's') + '</span></li>').join('') + '</ol>' +
      (ev.length ? '<h4>Evidence shown</h4><ul class="cite-list">' + ev.map((e) => '<li><b>' + esc(e.label) + '</b> <span class="muted">(' + fmt(e.n) + ' site' + (e.n === 1 ? '' : 's') + ')</span>: ' + e.refs.map(refLine).join('; ') + '</li>').join('') + '</ul>' : '') +
      (studies.some((x) => x.name === 'PhosphoSitePlus') ? '<p class="note">' + esc(PSP_NOTICE) + '</p>' : '') +
      '<p class="screen-only"><button class="btn quiet small" type="button" data-copy-refs="' + esc(text) + '">Copy references</button></p></details>';
  }
  function bindCiteBox(root) {
    $$('[data-copy-refs]', root).forEach((b) => b.addEventListener('click', (e) => copyText(b.dataset.copyRefs, e.currentTarget)));
    bindDocLinks(root);
  }
  const MOTIF_LABEL = {
    RG: 'RG', RGG: 'RGG', GRG: 'GRG', GAR: 'GAR', RXR: 'RXR',
    CARM1_proline_rich: 'Pro ±5', CARM1_hydrophobic_proline: 'Pro ±5 + Φ ±3',
    DR: 'DR', RD: 'RD', ER: 'ER', RE: 'RE',
  };
  const MOTIF_HELP = {
    RG: 'Gly at +1', RGG: 'Gly at +1 and +2', GRG: 'Gly at −1 and +1', GAR: 'Gly at −2 and Ala at −1',
    RXR: 'Another Arg at −2 or +2', CARM1_proline_rich: 'Pro within ±5 (CARM1-like)',
    CARM1_hydrophobic_proline: 'Pro within ±5 and a hydrophobic residue (AILMFWVY) within ±3 (CARM1-like)',
    DR: 'Asp at −1', RD: 'Asp at +1', ER: 'Glu at −1', RE: 'Glu at +1',
  };
  const DISORDER_LABEL = { disordered: 'IDR', disorder_boundary: 'IDR edge ≤20 aa', ordered: 'Ordered', no_disorder_annotation: 'No disorder data' };
  const DISORDER_HELP = {
    disordered: 'Inside a MobiDB disorder interval',
    disorder_boundary: 'Outside, but within 20 aa of, a disorder interval',
    ordered: 'More than 20 aa from any disorder interval',
    no_disorder_annotation: 'Protein has no MobiDB disorder annotation',
  };
  const DOMAIN_LABEL = { in_domain: 'In domain', boundary: 'Domain edge ≤20 aa', inter_domain_linker: 'Inter-domain linker', distal: 'Domain-distal', no_domain_annotation: 'No domain data' };
  const DOMAIN_HELP = {
    in_domain: 'Inside an InterPro domain or repeat',
    boundary: 'Outside, but within 20 aa of, a domain edge',
    inter_domain_linker: 'Between two domains, more than 20 aa from both',
    distal: 'More than 20 aa from any domain and not between two domains',
    no_domain_annotation: 'Protein has no InterPro domain or repeat annotation',
  };
  // How each evidence record was placed on the canonical sequence (PTM_Pipeline remapping).
  const REMAP_HELP = {
    same_accession_position_verified: 'Same accession; position checked',
    sequence_unique_canonical: 'Peptide fits one place in the canonical sequence',
    sequence_unique_candidate: 'Peptide fits one place (from another isoform)',
    sequence_closest_supported: 'Ambiguous: peptide fits several places; nearest supported position used',
    sequence_position_verified: 'Reported position is an arginine inside the matching peptide',
  };
  const TREAT_HELP = { GSK591: 'PRMT5 inhibitor', MS023: 'type I PRMT inhibitor' };
  const PTM_LABEL = ['Phospho', 'Lys acyl', 'Ub/SUMO', 'Lys methyl'];
  const PTM_PREFIX = ['p', 'ac', 'ub', 'me'];
  const PADJ = 0.05;
  // Maron 2021 methyl-IP antibodies named by the residue state they recognise.
  const IP_LABEL = { MMA: 'Rme1', ADMA: 'Rme2a', SDMA: 'Rme2s' };
  const IP_HELP = { MMA: 'monomethylarginine (MMA) antibody', ADMA: 'asymmetric dimethylarginine (ADMA) antibody', SDMA: 'symmetric dimethylarginine (SDMA) antibody' };
  const ipLabel = (name) => IP_LABEL[name] || name;
  const TREAT_OPTIONS = [
    ['', 'Any site'],
    ['tested', 'Measured in Maron Table S3'],
    ['0:down', 'Decreased with GSK591 (PRMT5i)'],
    ['1:down', 'Decreased with MS023 (type I PRMTi)'],
    ['0:up', 'Increased with GSK591'],
    ['1:up', 'Increased with MS023'],
  ];
  const RG_HELP = 'At least 3 RG dipeptides with ≤12 aa between consecutive RG arginines, found in the UniProt sequence; RGG-rich when ≥2 are RGG';

  // Site colouring for the protein figure. Category order is fixed so a colour always means the same thing.
  const COLOR_MODES = {
    state: {
      name: 'Methylation state', cats: ['Rme1', 'Rme2', 'Rme1 + Rme2', 'No state call'],
      cat: (s) => [3, 0, 1, 2][s.st], colors: (pal) => [pal.cat[0], pal.cat[1], pal.cat[2], pal.none],
    },
    ip: {
      name: 'Maron methyl-IP', cats: ['Rme2s IP', 'Rme2a IP (no Rme2s)', 'Rme1 IP only', 'Not in Maron IP'],
      cat: (s) => (s.ip & 4 ? 0 : s.ip & 2 ? 1 : s.ip & 1 ? 2 : 3), colors: (pal) => [pal.cat[0], pal.cat[1], pal.cat[2], pal.none],
    },
    prmti: {
      name: 'PRMT inhibitor response', cats: ['↓ GSK591 (PRMT5i)', '↓ MS023 (type I PRMTi)', '↓ both', 'Measured, not decreased', 'Not measured'],
      cat: (s) => { if (!s.tr.length) return 4; const g = s.trf & 1, m = s.trf & 4; return g && m ? 2 : g ? 0 : m ? 1 : 3; },
      colors: (pal) => [pal.cat[0], pal.cat[1], pal.cat[2], pal.dark, pal.none],
    },
    reports: {
      name: 'Independent reports', cats: ['1 report', '2 reports', '3–4 reports', '5 or more'],
      cat: (s) => { const n = popcount(s.study); return n >= 5 ? 3 : n >= 3 ? 2 : n === 2 ? 1 : 0; }, colors: (pal) => pal.seq,
    },
    ixstat: {
      name: 'Interaction evidence', cats: ['Experimental contact', 'Predicted only', 'None recorded'],
      cat: (s) => { const st = ixStatus(s); return st === 'exp' ? 0 : st === 'pred' ? 1 : 2; }, colors: (pal) => [pal.ixExp || '#0E7490', pal.ixPred || '#C2410C', pal.none],
    },
    single: { name: 'Single colour', cats: ['Methylarginine'], cat: () => 0, colors: (pal) => [pal.site] },
  };
  // Protein map: overlays on the backbone, then rows the reader switches on and reorders in the Rows menu (remembered
  // in this browser). Lysine methylation and the other PTMs come from the pipeline's cross-PTM master table.
  const MAP_OVERLAYS = [['domains', 'Domains'], ['dense', 'Dense methyl-R windows']];
  const MAP_ROWS = [['regions', 'Regions'], ['disorder', 'Disorder'], ['charge', 'Charge (EMBOSS, 5-aa mean)'], ['bias', 'Composition'], ['rg', 'RG/RGG-rich'],
    ['args', 'All arginines'], ['variants', 'Variants and mutagenesis'], ['kme', 'Lys methylation'], ['phospho', 'Phosphorylation'], ['kacyl', 'Lys acylation'],
    ['ub', 'Ubiquitin / SUMO'], ['iface', 'PDB interfaces'], ['xl', 'RNA/DNA crosslinks']];
  const MAP_PTM_KEYS = ['phospho', 'kacyl', 'ub', 'kme'];
  const MAP_DEFAULT_ON = ['domains', 'dense', 'regions', 'disorder', 'charge', 'rg', 'args', 'variants', 'kme', 'phospho', 'kacyl', 'ub', 'iface', 'xl'];
  const KME_COLOR = '#A21CAF';
  const mapPrefs = (() => { try { return JSON.parse(store.get('mre-map-rows') || 'null') || {}; } catch (e) { return {}; } })();
  const mapRowOrder = () => {
    const keys = MAP_ROWS.map(([k]) => k), order = (UI.fig.order || []).filter((k) => keys.includes(k));
    return order.concat(keys.filter((k) => !order.includes(k)));
  };
  const saveMapPrefs = () => store.set('mre-map-rows', JSON.stringify({ on: Array.from(UI.fig.tracks), order: mapRowOrder() }));
  const CHARGE_AA = { K: 1, R: 1, H: 0.5, D: -1, E: -1 };
  // Mean charge of the residues centred on each position, i-2 to i+2, on the EMBOSS charge scale. Near either
  // terminus the window is cut short and the mean is taken over the residues actually present (3 or 4), so every
  // residue has a value and the ends are not diluted by imaginary neutral residues.
  function chargeProfile(u) {
    if (!u || !u.seq) return null;
    if (!u.charge5) {
      const q = u.seq, L = q.length, out = new Array(L);
      for (let i = 0; i < L; i++) {
        let t = 0, n = 0;
        for (let k = Math.max(0, i - 2); k <= Math.min(L - 1, i + 2); k++) { t += CHARGE_AA[q[k]] || 0; n++; }
        out[i] = t / n;
      }
      u.charge5 = out;
    }
    return u.charge5;
  }
  // AlphaFold DB confidence bands and colours.
  const PLDDT_BANDS = [[90, '#0053D6', 'Very high (> 90)'], [70, '#65CBF3', 'Confident (70–90)'], [50, '#FFDB13', 'Low (50–70)'], [-Infinity, '#FF7D45', 'Very low (< 50)']];
  const plddtBand = (v) => PLDDT_BANDS.findIndex((b) => v > b[0]);
  // AlphaMissense thresholds (mean pathogenicity of the 19 R->X substitutions at a position).
  const amClass = (v) => (v >= 0.564 ? 'likely pathogenic' : v >= 0.34 ? 'ambiguous' : 'likely benign');
  // Light, self-contained palette for exported figures.
  const EXPORT_PAL = {
    bg: '#FFFFFF', ink: '#111827', muted: '#5F6F78', grid: '#E0E8EB', backbone: '#9AA6B2', stem: '#C7D4D9',
    backboneFill: '#EEF5F7', backboneEdge: '#C7D4D9', domainFill: '#B9DDE6', domainInk: '#0F3B47',
    disFill: '#D2DBDF', disEdge: '#8E9DA4', regionFill: '#E3ECF6', biasFill: '#FDE3D4', rgFill: '#F5D6DE', rgInk: '#7A2E45',
    dense: 'rgba(194,65,12,0.12)', denseInk: '#9A3412', ptm: '#6B7F88', rtick: '#C7D4D9', ring: '#FFFFFF',
    site: '#B4485F', brush: 'rgba(194,65,12,0.16)', accent: '#C2410C', variant: '#9AA6B2',
    cat: ['#2A78D6', '#EB6834', '#1BAF7A'], none: '#A3AEB3', dark: '#5F6F78', seq: ['#86B6EF', '#3987E5', '#1C5CAB', '#0D366B'],
    bar: '#0E7490', barNeg: '#5F6F78', font: 'Arial, Helvetica, "Liberation Sans", sans-serif',
  };

  // ------------------------------------------------------------------ model
  let M = null;      // meta
  let P = [];        // proteins
  let S = [];        // sites
  let IDX = null;    // lookup maps
  let R = null;      // current filtered result
  let evidencePromise = null;
  const uniprotCache = new Map();
  const UNIDATA = new Map();   // accession -> resolved UniProt features (or null)
  const afCache = new Map();
  const AFDATA = new Map();    // accession -> resolved AlphaFold data (or null)
  let mol = null;              // live 3Dmol viewer for the open protein

  const DEFAULT_FILTERS = () => ({
    // Motif, study, domain class, functional class, disease, condensate, complex, GO and partner filters are search
    // terms (motif:RGG …) held in Q, not here.
    states: new Set(), ips: new Set(), sources: new Set(), s2: new Set(), multi: false, reviewed: false, minRep: 0,
    treat: '', disorder: new Set(), domctx: new Set(),
    nn: 0, ptmGroup: '', ptmDist: 0,
    anyCond: false, minSites: 1, dense: false, arch: '', argEnriched: false,
    ixScore: 700, goExp: false, iface: new Set(), condConf: 0, modelledPdb: false, pred: 0, galOnly: false, xl: new Set(),
    evAll: 0, evExact: false, pevAll: 0, pevExact: false,
  });
  const cloneFilters = (f) => {
    const out = {};
    for (const k of Object.keys(f)) out[k] = f[k] instanceof Set ? new Set(f[k]) : f[k];
    return out;
  };
  let F = DEFAULT_FILTERS();
  let Q = { raw: '', plain: '', proteins: null, positions: null, note: '', tokens: [], errors: [], siteT: [], protT: [] };
  const UI = {
    view: 'overview', detail: null, highlight: null,
    sortP: { key: 'sites', dir: -1 }, sortS: { key: 'gene', dir: 1 },
    pageP: 50, pageS: 100, zoom: null, openEvidence: new Set(),
    allCond: false, funcExpanded: false, charts: {},
    fig: { color: 'state', labels: 'auto', tracks: new Set(Array.isArray(mapPrefs.on) ? mapPrefs.on : MAP_DEFAULT_ON), order: Array.isArray(mapPrefs.order) ? mapPrefs.order : null },
    ixTab: 'proteins', ixTop: 20, ixMin: 700, ixColor: 'meth', ixAll: false, ixNet: null,
    cxPred: false, cxMin: 3, cxSort: 'q', cxAll: false, refsOpen: new Set(), ctx: 'protein', condThr: 3, pendingContact: null, pendingModel: null, ixOpen: new Set(), ixJump: null,
    gal: { conf: 'confident', pl: 80, kind: 'bond', noPdb: false, xl: false, all: false },
    figs: { flank: 7, rows: 50, order: 'protein', logo: 'freq', nProt: 8, scale: 'common', color: 'state', lanes: true, sortP: 'default' }, consPos: null, consFlank: 10, docJump: null,
    evTop: 1, evAllRows: false, evMode: 'sites', cxOpen: new Set(), ixScrollMap: false, reportPrint: false,
  };

  // ------------------------------------------------------------------- load
  fetch('data/methylarg_sites.json')
    .then((r) => { if (!r.ok) throw new Error(r.status + ' ' + r.statusText); return r.json(); })
    .then(init)
    .catch((err) => {
      $('#loading').innerHTML = '<div class="empty"><b>The site data did not load</b>' +
        'Could not read <code>data/methylarg_sites.json</code> (' + esc(err.message) + '). ' +
        'If you opened index.html from disk, serve the folder instead, for example <code>python -m http.server</code> inside <code>docs/</code>.</div>';
    });

  function loadEvidence() {
    if (!evidencePromise) {
      evidencePromise = fetch('data/methylarg_evidence.json')
        .then((r) => { if (!r.ok) throw new Error(r.status); return r.json(); })
        .catch(() => null);
    }
    return evidencePromise;
  }

  function loadUniProt(acc) {
    if (uniprotCache.has(acc)) return uniprotCache.get(acc);
    const fields = 'accession,length,sequence,protein_name,gene_primary,ft_domain,ft_repeat,ft_zn_fing,ft_region,ft_compbias,ft_motif,ft_variant,ft_mutagen,ft_mod_res,cc_disease,xref_ensembl';
    const ctrl = typeof AbortController === 'function' ? new AbortController() : null;
    const timer = ctrl ? setTimeout(() => ctrl.abort(), 20000) : null;
    const promise = fetch('https://rest.uniprot.org/uniprotkb/' + encodeURIComponent(acc) + '.json?fields=' + fields,
      ctrl ? { signal: ctrl.signal } : undefined)
      .then((r) => { if (!r.ok) throw new Error(r.status); return r.json(); })
      .then((d) => {
        const desc = d.proteinDescription || {};
        const nameObj = desc.recommendedName || (desc.submissionNames || [])[0] || {};
        const seq = (d.sequence && d.sequence.value) || '';
        const all = (d.features || []).map((f) => {
          const alt = f.alternativeSequence || {};
          return {
            type: f.type, desc: f.description || '',
            start: f.location && f.location.start && f.location.start.value,
            end: f.location && f.location.end && f.location.end.value,
            from: alt.originalSequence || '', to: (alt.alternativeSequences || []).join('/'), id: f.featureId || '',
          };
        }).filter((f) => f.start && f.end);
        const special = ['Natural variant', 'Mutagenesis', 'Modified residue'];
        return {
          seq,
          len: (d.sequence && d.sequence.length) || null,
          name: (nameObj.fullName && nameObj.fullName.value) || '',
          gene: (d.genes && d.genes[0] && d.genes[0].geneName && d.genes[0].geneName.value) || '',
          feats: all.filter((f) => !special.includes(f.type)),
          variants: all.filter((f) => f.type === 'Natural variant'),
          mutagen: all.filter((f) => f.type === 'Mutagenesis'),
          methylR: new Map(all.filter((f) => f.type === 'Modified residue' && /arginine/i.test(f.desc)).map((f) => [f.start, f.desc])),
          rg: rgBoxes(seq),
          diseases: (d.comments || []).filter((c) => c.commentType === 'DISEASE').map((c) => ({
            name: (c.disease && c.disease.diseaseId) || '', acronym: (c.disease && c.disease.acronym) || '', id: (c.disease && c.disease.diseaseAccession) || '',
            mim: ((c.disease && c.disease.diseaseCrossReference) || {}).id || '', desc: (c.disease && c.disease.description) || '',
            note: ((c.note && c.note.texts) || []).map((t) => t.value).join(' '),
          })).filter((x) => x.name),
          ensg: Array.from(new Set((d.uniProtKBCrossReferences || []).filter((x) => x.database === 'Ensembl')
            .map((x) => ((x.properties || []).find((pp) => pp.key === 'GeneId') || {}).value || '').filter(Boolean).map((g) => g.split('.')[0]))),
        };
      })
      .catch(() => null)
      .then((d) => { UNIDATA.set(acc, d); return d; })
      .finally(() => { if (timer) clearTimeout(timer); });
    uniprotCache.set(acc, promise);
    return promise;
  }

  // UniProt's own modified-residue term for a methylarginine, shortened.
  function uniprotMethylLabel(desc) {
    const d = desc.toLowerCase();
    const t = d.includes('asymmetric') ? 'Rme2a' : d.includes('symmetric') ? 'Rme2s' : d.includes('dimethyl') ? 'Rme2' : d.includes('methyl') ? 'Rme1' : 'modified';
    return t + (d.includes('alternate') ? ' (alt.)' : '');
  }

  // ------------------------------------------------------------- AlphaFold
  const fetchJson = (u) => fetch(u).then((r) => { if (!r.ok) throw new Error(r.status); return r.json(); });
  const fetchText = (u) => fetch(u).then((r) => { if (!r.ok) throw new Error(r.status); return r.text(); });

  // AlphaFold DB entry for an accession: per-residue pLDDT and, for arginines, AlphaMissense R->X scores.
  function loadAF(acc) {
    if (afCache.has(acc)) return afCache.get(acc);
    const pr = fetchJson('https://alphafold.ebi.ac.uk/api/prediction/' + encodeURIComponent(acc))
      .then((list) => {
        const entries = Array.isArray(list) ? list : [];
        const e = entries.find((x) => x.uniprotAccession === acc) || entries[0];
        if (!e) return null;
        return Promise.all([
          e.plddtDocUrl ? fetchJson(e.plddtDocUrl).catch(() => null) : null,
          e.amAnnotationsUrl ? fetchText(e.amAnnotationsUrl).catch(() => null) : null,
        ]).then(([conf, am]) => {
          const plddt = new Map();
          if (conf && conf.residueNumber) conf.residueNumber.forEach((r, k) => plddt.set(r, conf.confidenceScore[k]));
          const amR = new Map();
          if (am) {
            for (const line of am.split('\n')) {
              if (line.charCodeAt(0) !== 82) continue; // rows starting with 'R' (arginine substitutions)
              const [v, sc, cls] = line.split(',');
              const m = /^R(\d+)[A-Z]$/.exec(v);
              if (!m) continue;
              const pos = +m[1];
              let a = amR.get(pos);
              if (!a) amR.set(pos, (a = { sum: 0, n: 0, lp: 0, max: 0 }));
              const s = +sc;
              a.sum += s; a.n++; a.max = Math.max(a.max, s);
              if (/path/i.test(cls || '')) a.lp++;
            }
          }
          return { entry: e, plddt, am: amR, seq: e.uniprotSequence || e.sequence || '' };
        });
      })
      .catch(() => null)
      .then((d) => { AFDATA.set(acc, d); return d; });
    afCache.set(acc, pr);
    return pr;
  }

  function load3Dmol() {
    if (window.$3Dmol) return Promise.resolve(window.$3Dmol);
    if (!load3Dmol.p) {
      load3Dmol.p = new Promise((resolve, reject) => {
        const s = document.createElement('script');
        s.src = 'https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.5.3/3Dmol-min.js';
        s.async = true;
        s.onload = () => (window.$3Dmol ? resolve(window.$3Dmol) : reject(new Error('3Dmol missing')));
        s.onerror = () => { load3Dmol.p = null; reject(new Error('3Dmol failed to load')); };
        document.head.appendChild(s);
      });
    }
    return load3Dmol.p;
  }

  // RG/RGG-rich segments: runs of >=3 RG dipeptides with <=12 aa between consecutive RG arginines.
  function rgBoxes(seq) {
    const out = [];
    let cur = null;
    const close = () => { if (cur && cur.n >= 3) out.push(cur); };
    for (let k = 0; k < seq.length - 1; k++) {
      if (seq[k] !== 'R' || seq[k + 1] !== 'G') continue;
      const pos = k + 1;
      const rgg = seq[k + 2] === 'G' ? 1 : 0;
      if (cur && pos - cur.last <= 12) { cur.last = pos; cur.n++; cur.rgg += rgg; cur.end = pos + 1 + rgg; }
      else { close(); cur = { start: pos, last: pos, end: pos + 1 + rgg, n: 1, rgg }; }
    }
    close();
    return out.map((c) => ({ start: c.start, end: c.end, n: c.n, rgg: c.rgg, label: c.rgg >= 2 ? 'RGG-rich' : 'RG-rich' }));
  }

  function init(data) {
    M = data.meta;
    M.vi = { disordered: M.disorder_classes.indexOf('disordered') };
    M.studyOrder = M.studies.map((st, k) => k).sort((a, b) => M.studies[b][2] - M.studies[a][2]);
    M.subParent = {};
    for (const r of M.overview.functional_subclasses || []) M.subParent[r[0]] = r[6];
    P = data.proteins.map((r, i) => ({
      i, acc: r[0], gene: r[1] || '', name: r[2] || '', uid: r[3] || '', func: r[4], sub: r[5], len: r[6],
      argN: r[7], shrink: r[8], argQ: r[9], iso: !!r[10], cond: r[11], dis: r[12], als: !!r[13],
      dens: r[14], best: r[15], sig: r[16], arch: r[17], ptms: r[18], aliases: r[19] || [], sites: [],
    }));
    S = data.sites.map((r, j) => ({
      j, p: r[0], pos: r[1], st: r[2], src: r[3], refs: r[4], dc: r[5], dEdge: r[6], dctx: r[7], ipr: r[8],
      dcls: r[9], domDist: r[10], domSide: r[11], win: r[12], mot: r[13], mis: !!r[14], orig: r[15],
      s2: r[16], ip: r[17], tr: r[18], study: r[19] || 0, nn: null, ptmD: null, trf: 0,
      iface: 0, ifP: null, ifPdb: null, ifN: 0, contacts: null, modelledMe: '', geo: null, modelled: null,
    }));
    for (const s of S) P[s.p].sites.push(s.j);
    // Unreviewed (TrEMBL) entries have no Swiss-Prot entry name. Their sites often repeat a site of the reviewed entry
    // of the same gene at the same position; keep both (the union as published) but flag the repeat.
    const reviewedByGene = new Map();
    for (const p of P) {
      p.unrev = !p.uid || p.uid === p.acc.split('-')[0] + '_HUMAN';
      if (!p.unrev && p.gene) { if (!reviewedByGene.has(p.gene)) reviewedByGene.set(p.gene, []); reviewedByGene.get(p.gene).push(p); }
    }
    for (const p of P) {
      if (!p.unrev) continue;
      p.reviewedTwin = (reviewedByGene.get(p.gene) || [])[0] || null;
      for (const j of p.sites) {
        for (const q of reviewedByGene.get(p.gene) || []) {
          const k = q.sites.find((x) => S[x].pos === S[j].pos);
          if (k != null) { S[j].dup = k; break; }
        }
      }
    }

    for (const p of P) {
      p.hist = HISTONE_GENE.test(p.gene) || (!p.gene && HISTONE_NAME.test(p.name));
      p.label = p.gene || p.acc;
      p.maxPos = p.sites.length ? S[p.sites[p.sites.length - 1]].pos : 0;
      p.ptmPos = [[], [], [], []];
      for (const [pos, g] of p.ptms) p.ptmPos[g].push(pos);
      p.ptmPos.forEach((a) => a.sort((x, y) => x - y));
      const ss = p.sites;
      for (let k = 0; k < ss.length; k++) {
        const s = S[ss[k]];
        const prev = k > 0 ? s.pos - S[ss[k - 1]].pos : Infinity;
        const next = k < ss.length - 1 ? S[ss[k + 1]].pos - s.pos : Infinity;
        const d = Math.min(prev, next);
        s.nn = Number.isFinite(d) ? d : null;
        s.ptmD = p.ptmPos.map((arr) => nearestDistance(arr, s.pos));
      }
      p.nIdr = ss.reduce((n, j) => n + (S[j].dc === M.vi.disordered ? 1 : 0), 0);
    }
    for (const s of S) {
      for (const [, t, diff, padj] of s.tr) {
        if (padj != null && padj <= PADJ && diff != null) s.trf |= 1 << (t * 2 + (diff < 0 ? 0 : 1));
      }
    }
    IDX = { gene: new Map(), acc: new Map(), uid: new Map(), alias: new Map() };
    const add = (map, key, i) => { if (!key) return; key = key.toUpperCase(); if (!map.has(key)) map.set(key, []); if (!map.get(key).includes(i)) map.get(key).push(i); };
    for (const p of P) {
      add(IDX.gene, p.gene, p.i);
      add(IDX.acc, p.acc, p.i);
      if (p.acc.includes('-')) add(IDX.acc, p.acc.split('-')[0], p.i);
      add(IDX.uid, p.uid, p.i);
      for (const a of p.aliases) add(IDX.alias, a, p.i);
    }

    $('#data-version').textContent = 'data exported ' + M.generated;
    if (window.matchMedia && window.matchMedia('(max-width: 560px)').matches) $('#q').placeholder = 'Gene, site, gene list or motif:RGG';
    $('#year').textContent = new Date().getFullYear();
    renderHero();
    const params = new URLSearchParams(location.search);
    F = filtersFromParams(params);
    const q0 = [params.get('q') || ''].concat(legacyTerms(params)).join(' ').trim();
    if (params.get('A') != null) UI.cmpA = { qs: params.get('A'), label: '' };
    if (q0) { $('#q').value = q0; $('#q-clear').hidden = false; Q = parseQuery(q0); }
    renderRail();
    bindGlobal();
    bindDocLinks($('footer.site'));
    compute();
    route();
    loadInteractions();
  }

  function nearestDistance(sorted, pos) {
    if (!sorted.length) return Infinity;
    let lo = 0, hi = sorted.length - 1;
    while (lo < hi) { const mid = (lo + hi) >> 1; if (sorted[mid] < pos) lo = mid + 1; else hi = mid; }
    let d = Math.abs(sorted[lo] - pos);
    if (lo > 0) d = Math.min(d, Math.abs(sorted[lo - 1] - pos));
    return d;
  }

  // ----------------------------------------------------------------- search
  // Exact lookup: accession, current symbol, UniProt entry name, then synonyms and source-reported symbols.
  function lookup(token) {
    const key = token.toUpperCase();
    const hit = IDX.acc.get(key) || IDX.gene.get(key) || IDX.uid.get(key);
    if (hit) return { hit, alias: false };
    const al = IDX.alias.get(key);
    return al ? { hit: al, alias: true } : null;
  }
  const lookupExact = (token) => { const l = lookup(token); return l ? l.hit : null; };

  function parseQuery(text) {
    const raw = text.trim();
    const { plain, tokens, errors } = splitTerms(raw);
    const q = parsePlain(plain);
    Object.assign(q, { raw, plain, tokens, errors, siteT: tokens.filter((t) => t.vals[0].kind === 'site'), protT: tokens.filter((t) => t.vals[0].kind === 'prot') });
    if (errors.length) q.note = (q.note ? q.note + ' ' : '') + '<span class="warn-text">' + errors.map(esc).join(' ') + '</span>';
    return q;
  }
  function parsePlain(text) {
    const raw = text.trim();
    if (!raw) return { raw: '', proteins: null, positions: null, note: '' };
    const tokens = raw.split(/[\s,;]+/).filter(Boolean);
    const proteins = new Set();
    const positions = new Map();
    const missing = [];
    const missingSites = [];
    const viaAlias = [];
    let last = null;
    let histoneNote = false;
    const addPos = (hit, typed) => {
      // Histone positions are typed in histone numbering (no initiator Met): H3C1 R2 is UniProt R3.
      const uni = (i) => typed + (P[i].hist ? 1 : 0);
      for (const i of hit) {
        if (!positions.has(i)) positions.set(i, new Set());
        positions.get(i).add(uni(i));
        proteins.add(i);
        if (P[i].hist) histoneNote = true;
      }
      // A symbol can name several entries (isoforms, unreviewed twins): warn only when none of them has the site.
      if (!hit.some((i) => P[i].sites.some((j) => S[j].pos === uni(i)))) missingSites.push(P[hit[0]].label + ' R' + typed);
    };
    const noteAlias = (tok, l) => { if (l.alias) viaAlias.push(tok + ' → ' + l.hit.map((i) => P[i].label).join('/')); };
    for (const tok of tokens) {
      const l = lookup(tok);
      if (l) { l.hit.forEach((i) => proteins.add(i)); last = l.hit; noteAlias(tok, l); continue; }
      let m = tok.match(/^(.+?)[:_/-]R?(\d+)$/i);
      if (m) {
        const lm = lookup(m[1]);
        if (lm) { addPos(lm.hit, +m[2]); last = lm.hit; noteAlias(m[1], lm); continue; }
      }
      m = tok.match(/^R?(\d+)$/i);
      if (m && last) { addPos(last, +m[1]); continue; }
      missing.push(tok);
    }
    if (!proteins.size) {
      const t = raw.toLowerCase();
      const hits = P.filter((p) => p.gene.toLowerCase().startsWith(t) || p.acc.toLowerCase().startsWith(t) ||
        p.aliases.some((a) => a.toLowerCase().startsWith(t)) || (t.length >= 3 && p.name.toLowerCase().includes(t)));
      if (hits.length) {
        return { raw, proteins: new Set(hits.map((p) => p.i)), positions: null,
          note: '<b>' + fmt(hits.length) + '</b> protein' + (hits.length === 1 ? '' : 's') + ' match “' + esc(raw) + '” by symbol prefix, synonym or protein name.' };
      }
      return { raw, proteins: new Set(), positions: null,
        note: 'No protein in the union matches “' + esc(raw) + '”. Try a gene symbol (PRMT5), a UniProt accession (O14744) or part of a protein name.' };
    }
    const notes = [];
    const matched = tokens.length - missing.length;
    // A gene list gets a count; a single protein with its sites needs none.
    if (tokens.length > 1 && (proteins.size > 1 || missing.length)) notes.push('<b>' + fmt(proteins.size) + '</b> protein' + (proteins.size === 1 ? '' : 's') + ' from ' + fmt(matched) + ' of ' + fmt(tokens.length) + ' search words.');
    if (viaAlias.length) notes.push('Matched by synonym: ' + viaAlias.slice(0, 10).map(esc).join(', ') + '.');
    if (missing.length) notes.push('Not in the methylarginine union: ' + missing.slice(0, 25).map(esc).join(', ') + (missing.length > 25 ? ' and ' + (missing.length - 25) + ' more' : '') + '.');
    if (missingSites.length) notes.push('No methylation reported at ' + missingSites.slice(0, 10).map(esc).join(', ') + '.');
    if (histoneNote) notes.push('Histone residues are numbered without the initiator methionine (H3 R2 is UniProt R3).');
    return { raw, proteins, positions: positions.size ? positions : null, note: notes.join(' ') };
  }

  function suggestions(text) {
    const tok = text.split(/[\s,;]+/).pop() || '';
    if (!tok) return [];
    const t = tok.toLowerCase();
    const exact = [], prefix = [], alias = [], named = [];
    for (const p of P) {
      const g = p.gene.toLowerCase();
      if (g === t || p.acc.toLowerCase() === t) exact.push([p, '']);
      else if (g.startsWith(t) || p.acc.toLowerCase().startsWith(t)) prefix.push([p, '']);
      else {
        const a = p.aliases.find((x) => x.toLowerCase().startsWith(t));
        if (a) alias.push([p, a]);
        else if (t.length >= 3 && p.name.toLowerCase().includes(t)) named.push([p, '']);
      }
    }
    const bySites = (a, b) => b[0].sites.length - a[0].sites.length;
    return exact.concat(prefix.sort(bySites), alias.sort(bySites), named.sort(bySites)).slice(0, 8);
  }

  // ----------------------------------------------------------- search terms
  // Readable filters inside the search, e.g. motif:RGG window:PRPR class:rna-binding disease:als go:rna binds:PRMT5
  // interaction:predicted. A comma inside a term means any of its values (motif:RGG,GAR); separate terms must all
  // match. Every other word still searches genes, accessions and sites, so the search box and the link say exactly
  // what is selected.
  const TERM_KEYS = ['motif', 'window', 'study', 'domain', 'class', 'subclass', 'disease', 'condensate', 'complex', 'go', 'binds', 'interaction'];
  const TERM_GROUP = { motif: 'motif', window: 'motif', study: 'prov', domain: 'struct', class: 'protein', subclass: 'protein', disease: 'protein', condensate: 'protein', complex: 'ix', go: 'ix', binds: 'ix', interaction: 'ix' };
  const TERM_HELP = {
    motif: 'Motif covering the methyl-R: a family (RG, RGG, GRG, GAR, RXR, DR, RD, ER, RE, CARM1-Pro) or a pattern such as GXRG or [RK]XR',
    window: 'Motif anywhere within ±15 aa of the site, e.g. PRPR',
    study: 'Reported by a study or database, e.g. maron-2021, phosphositeplus',
    domain: 'Nearest InterPro domain class (lab domain ontology), e.g. rrm, rna-helicase',
    class: 'Functional class (lab protein ontology), e.g. rna-binding, chromatin',
    subclass: 'Functional subclass (lab protein ontology), e.g. hnrnp, chromatin-readers',
    disease: 'Disease family (DisGeNET), e.g. als, cancer; als-all is the full audited ALS/MN/FTD gene set',
    condensate: 'CD-CODE condensate, e.g. stress-granule, nuclear-speckle',
    complex: 'Complex Portal complex, by ID (CPX-1234) or name',
    go: 'Gene Ontology binding: rna or dna',
    binds: 'STRING physical partner: prmt, tudor, wd40, phd, methylated, other, or a gene symbol',
    interaction: 'Contact evidence at the site: experimental, predicted, predicted-only, any or none',
  };
  const TERM_NOUN = { study: 'study or database', domain: 'domain class', class: 'functional class', subclass: 'functional subclass', disease: 'disease family', condensate: 'condensate', complex: 'complex' };
  const MOTIF_TERM = { CARM1_proline_rich: 'CARM1-Pro', CARM1_hydrophobic_proline: 'CARM1-Pro-hydrophobic' };
  const motifTerm = (m) => MOTIF_TERM[m] || m;
  const motifFamily = (v) => M.motifs.findIndex((m) => motifTerm(m).toLowerCase() === String(v).toLowerCase());
  const BIND_TERMS = [['prmt', 1], ['tudor', 2], ['wd40', 4], ['phd', 5], ['methylated', 3], ['other', 6]];
  const IX_STATUS_TERMS = [
    ['experimental', 'Experimental contact evidence', (m) => m & EV_EXP],
    ['predicted', 'Predicted interface or contact', (m) => m & EV_PRED],
    ['predicted-only', 'Predicted only (no experimental evidence)', (m) => (m & EV_PRED) && !(m & EV_EXP)],
    ['any', 'Any contact evidence', (m) => m],
    ['none', 'No contact evidence recorded', (m) => !m],
  ];
  const slugify = (s) => String(s).toLowerCase().replace(/['’]/g, '').replace(/[^a-z0-9+]+/g, '-').replace(/^-+|-+$/g, '');
  const termNameCache = {};
  // Names a term can take, each with a short slug (the part before the first " / ", made unique within its list).
  function termNames(key) {
    if (termNameCache[key]) return termNameCache[key];
    if (key === 'complex' && !ixReady()) return [];
    const pairs = key === 'class' ? M.func_classes.map((n, k) => [k, n])
      : key === 'subclass' ? M.func_subclasses.map((n, k) => [k, n])
        : key === 'disease' ? M.disease_families.map((n, k) => [k, n])
          : key === 'condensate' ? M.condensates.map((c, k) => [k, c[0]])
            : key === 'study' ? M.studies.map((st, k) => [k, st[0]])
              : key === 'domain' ? M.domain_classes.map((n, k) => [k, n])
                : key === 'complex' ? IX.cx.map((c) => [c.k, c.name]) : [];
    let out;
    if (key === 'complex') out = pairs.map(([k, n]) => ({ k, name: n, slug: IX.cx[k].id }));
    else {
      const short = pairs.map(([, n]) => slugify(n.split(' / ')[0]));
      const seen = {};
      short.forEach((s) => { seen[s] = (seen[s] || 0) + 1; });
      out = pairs.map(([k, n], i) => ({ k, name: n, slug: seen[short[i]] > 1 ? slugify(n) : short[i] }));
      const again = {};
      out.forEach((x) => { again[x.slug] = (again[x.slug] || 0) + 1; });
      out.forEach((x) => { if (again[x.slug] > 1) x.slug += '-' + (key === 'condensate' ? slugify(M.condensates[x.k][1]) : x.k); });
    }
    if (key === 'disease') out.push({ k: 'als', name: 'ALS / MN / FTD, full audited gene set', slug: 'als-all' });
    return (termNameCache[key] = out);
  }
  const termSlug = (key, k) => { const x = termNames(key).find((y) => y.k === k); return x ? x.slug : ''; };
  // Exact slug or name first, then slug or name prefix, then any substring; the shortest name wins a tie.
  function pickName(list, v) {
    const t = v.toLowerCase(), sv = slugify(v);
    if (!sv) return null;
    const tiers = [(x) => x.slug.toLowerCase() === sv || x.name.toLowerCase() === t, (x) => x.slug.toLowerCase().startsWith(sv) || x.name.toLowerCase().startsWith(t),
      (x) => x.name.toLowerCase().includes(t) || x.slug.toLowerCase().includes(sv)];
    for (const tier of tiers) {
      const hits = list.filter(tier);
      if (hits.length) return hits.sort((a, b) => a.name.length - b.name.length)[0];
    }
    return null;
  }
  // One value of a term: {text, label, kind: 'site' | 'prot', test(item, filters)} or {error}.
  function resolveValue(key, v) {
    if (key === 'motif' || key === 'window') {
      const fam = key === 'motif' ? motifFamily(v) : -1;
      if (fam >= 0) { const bit = 1 << fam; return { text: motifTerm(M.motifs[fam]), label: MOTIF_LABEL[M.motifs[fam]], kind: 'site', family: true, test: (s) => s.mot & bit }; }
      const m = parseMotif(v);
      if (!m || m.error) return { error: m ? m.error : 'Empty motif.' };
      const where = key === 'motif' ? 'site' : 'window';
      return { text: m.text, label: m.text, kind: 'site', test: (s) => motifMatch(s.win, m, where) };
    }
    if (TERM_NOUN[key]) {
      if (key === 'complex' && !IX) return { text: v, label: v, kind: 'prot', pending: true, test: () => false };
      if (key === 'complex' && IX.failed) return { error: 'The interaction data did not load, so complex: is unavailable.' };
      const x = pickName(termNames(key), v);
      if (!x) return { error: 'No ' + TERM_NOUN[key] + ' matches “' + v + '”.' };
      const k = x.k, bit = typeof k === 'number' ? 2 ** k : 0;
      const test = key === 'study' ? (s) => s.study & bit
        : key === 'domain' ? (s) => s.dcls === k
          : key === 'class' ? (p) => p.func & bit
            : key === 'subclass' ? (p) => p.sub & bit
              : key === 'disease' ? (k === 'als' ? (p) => p.als : (p) => p.dis & bit)
                : key === 'complex' ? (p) => (p.cx || []).includes(k)
                  : (p, f) => (f.condConf ? (cdMemberships(p) || []).some((c) => c.ci === k && c.score >= f.condConf) : p.cond.includes(k));
      return { text: x.slug, label: key === 'complex' ? x.name + ' (' + x.slug + ')' : x.name, kind: key === 'study' || key === 'domain' ? 'site' : 'prot', k, test };
    }
    if (key === 'go') {
      const w = v.toLowerCase();
      if (w !== 'rna' && w !== 'dna') return { error: 'go: takes rna or dna.' };
      return { text: w, label: w.toUpperCase() + ' binding', kind: 'prot',
        test: (p, f) => ixReady() && !!(p.ixFlags & (w === 'rna' ? (f.goExp ? IXF.rnaExp : IXF.rna) : (f.goExp ? IXF.dnaExp : IXF.dna))) };
    }
    if (key === 'binds') {
      const cls = BIND_TERMS.find(([t]) => t === v.toLowerCase());
      if (cls) return { text: cls[0], label: BIND_OPTS.find((o) => o[0] === cls[1])[1], kind: 'prot', bindClass: true, test: (p, f) => ixReady() && partnersAt(p, f.ixScore).some(([kk]) => bindMatch(cls[1], kk)) };
      const g = v.toUpperCase();
      if (!IX) return { text: g, label: g, kind: 'prot', pending: true, test: () => false };
      if (IX.failed) return { error: 'The interaction data did not load, so binds: is unavailable.' };
      let k = partnerByGene(g);
      if (k < 0) { const hit = lookupExact(g); if (hit && hit[0] < IX.nP) k = hit[0]; }
      if (k < 0) return { error: 'No protein with the gene symbol ' + g + ' in the interaction data.' };
      return { text: IX.prot[k].gene, label: IX.prot[k].gene, kind: 'prot', partner: k, test: (p, f) => partnersAt(p, f.ixScore).some(([kk]) => kk === k) };
    }
    if (key === 'interaction') {
      const w = v.toLowerCase();
      const st = IX_STATUS_TERMS.find(([t]) => t === w) || IX_STATUS_TERMS.find(([t]) => t.startsWith(w));
      if (!st) return { error: 'interaction: takes experimental, predicted, predicted-only, any or none.' };
      return { text: st[0], label: st[1], kind: 'site', test: (s) => ixReady() && !!st[2](s.ev || 0) };
    }
    return { error: 'Unknown term ' + key + ':' };
  }
  // key:value, key:"two words" or key:a,b at the start of the search or after a space, comma or semicolon.
  const TERM_RE = /(^|[\s,;])([A-Za-z]+):(?:"([^"]*)"?|([^\s"]*))/g;
  function splitTerms(raw) {
    const tokens = [], errors = [];
    const plain = raw.replace(TERM_RE, (all, pre, key, quoted, bare) => {
      const k = key.toLowerCase();
      if (!TERM_KEYS.includes(k)) return all;
      const vals = [];
      for (const v of (quoted != null ? [quoted] : bare.split(',')).map((x) => x.trim()).filter(Boolean)) {
        const r = resolveValue(k, v);
        if (r.error) errors.push(k + ':' + v + ': ' + r.error);
        else if (!vals.some((x) => x.text === r.text)) vals.push(r);
      }
      if (vals.length) tokens.push({ key: k, vals });
      return pre + ' ';
    });
    return { plain: plain.replace(/\s+/g, ' ').trim(), tokens, errors };
  }
  const quoteTerm = (t) => (/[\s",]/.test(t) ? '"' + t.replace(/"/g, '') + '"' : t);
  const tokenText = (t) => t.key + ':' + t.vals.map((v) => quoteTerm(v.text)).join(',');
  const queryText = (plain, tokens) => [plain || ''].concat(tokens.filter((t) => t.vals.length).map(tokenText)).join(' ').replace(/\s+/g, ' ').trim();
  const sameToken = (a, b) => a.key === b.key && tokenText(a) === tokenText(b);
  function tokenLabel(t) {
    const any = (vals) => vals.map((v) => v.label).join(' or ');
    switch (t.key) {
      case 'motif': return 'Motif ' + any(t.vals) + (t.vals.every((v) => v.family) ? '' : ' at the site');
      case 'window': return 'Motif ' + any(t.vals) + ' within ±15 aa';
      case 'study': return 'Reported by ' + any(t.vals);
      case 'domain': return 'Nearest domain: ' + any(t.vals);
      case 'class': return 'Class: ' + any(t.vals);
      case 'subclass': return 'Subclass: ' + any(t.vals);
      case 'disease': return 'Disease: ' + any(t.vals);
      case 'condensate': return 'Condensate: ' + any(t.vals) + (F.condConf ? ' (CD-CODE ≥ ' + F.condConf + '/5)' : '');
      case 'complex': return 'Complex: ' + any(t.vals);
      case 'go': return any(t.vals) + ' (GO' + (F.goExp ? ', experimental' : '') + ')';
      case 'binds': return 'Binds ' + any(t.vals) + ' (STRING ≥ ' + (F.ixScore / 1000).toFixed(1) + ')';
      case 'interaction': return 'Site: ' + any(t.vals).replace(/^./, (c) => c.toLowerCase());
      default: return tokenText(t);
    }
  }
  // Change the terms of the current search, keep its plain words, and refresh. mode: 'set' replaces every term of the
  // key, 'add' adds a separate term (narrows), 'toggle' adds or removes one value inside the key's chip term.
  function editTerms(key, text, mode, view, chipOf) {
    let toks = Q.tokens.map((t) => ({ key: t.key, vals: t.vals.slice() }));
    const same = (v) => v.text.toLowerCase() === String(text).toLowerCase();
    if (mode === 'set') toks = toks.filter((t) => t.key !== key);
    if (mode === 'toggle') {
      const own = toks.find((t) => t.key === key && t.vals.every((v) => !chipOf || chipOf(v)));
      if (own && own.vals.some(same)) own.vals = own.vals.filter((v) => !same(v));
      else if (own) own.vals.push({ text });
      else toks.push({ key, vals: [{ text }] });
    } else if (text !== '' && text != null) toks.push({ key, vals: [{ text }] });
    applyQueryText(queryText(Q.plain, toks), view);
  }
  const setTerm = (key, text, view) => editTerms(key, text, 'set', view);
  const addTerm = (key, text, view) => editTerms(key, text, 'add', view);
  function applyQueryText(text, view) {
    setQuery(text);
    UI.pageP = 50; UI.pageS = 100;
    compute();
    if (view && view !== UI.view) go(view); else renderView();
  }
  const termValues = (key) => Q.tokens.filter((t) => t.key === key).flatMap((t) => t.vals);
  // The single value of a key when the search holds exactly one, for context cards and rail selects.
  const singleTerm = (key) => { const t = Q.tokens.filter((x) => x.key === key); return t.length === 1 && t[0].vals.length === 1 ? t[0].vals[0] : null; };
  // Where a search made only of terms should land: sites when any term is about sites.
  const termView = () => (Q.siteT.length ? 'sites' : 'proteins');

  // Values offered while typing key:…, and keys offered while typing their first letters.
  function termSuggestions(tok) {
    const m = /^([A-Za-z]+):"?([^"]*)$/.exec(tok);
    if (!m) {
      const t = tok.toLowerCase();
      return t.length >= 2 ? TERM_KEYS.filter((k) => k.startsWith(t)).map((k) => ({ term: k + ':', label: k + ':', sub: TERM_HELP[k], open: true })) : [];
    }
    const key = m[1].toLowerCase(), part = m[2].split(',').pop().toLowerCase();
    if (!TERM_KEYS.includes(key)) return [];
    const head = m[2].includes(',') ? m[2].slice(0, m[2].lastIndexOf(',') + 1) : '';
    let list;
    if (key === 'motif') list = M.motifs.map((mm) => [motifTerm(mm), MOTIF_LABEL[mm] + ' · ' + MOTIF_HELP[mm]]);
    else if (key === 'window') list = [['PRPR', 'Pro-Arg repeats'], ['RGRG', 'RG repeat'], ['YGG', 'Tyr-Gly-Gly']];
    else if (key === 'go') list = [['rna', 'RNA binding (GO:0003723)'], ['dna', 'DNA binding (GO:0003677)']];
    else if (key === 'interaction') list = IX_STATUS_TERMS.map(([t, l]) => [t, l]);
    else if (key === 'binds') {
      list = BIND_TERMS.map(([t, v]) => [t, BIND_OPTS.find((o) => o[0] === v)[1]]);
      if (ixReady() && part.length >= 2) {
        const bound = boundBy(F.ixScore);
        list = list.concat(IX.prot.filter((x) => bound.get(x.k) && x.gene.toLowerCase().startsWith(part)).sort((a, b) => bound.get(b.k) - bound.get(a.k)).slice(0, 6)
          .map((x) => [x.gene, (partnerClassText(x) ? partnerClassText(x) + ' · ' : '') + fmt(bound.get(x.k)) + ' methylated partners']));
      }
    } else {
      const names = termNames(key).filter((x) => key !== 'condensate' || M.condensates[x.k][4] > 0);
      list = names.filter((x) => !part || x.slug.toLowerCase().startsWith(part) || x.name.toLowerCase().includes(part)).slice(0, 8).map((x) => [x.slug, x.name]);
      return list.map(([v, l]) => ({ term: key + ':' + head + v, label: key + ':' + v, sub: l }));
    }
    return list.filter(([v, l]) => !part || v.toLowerCase().startsWith(part) || l.toLowerCase().includes(part)).slice(0, 8).map(([v, l]) => ({ term: key + ':' + head + v, label: key + ':' + v, sub: l }));
  }

  // Links made before terms existed carried numeric filters (motifs=2, disease=3, func=1 …); read them as terms.
  function legacyTerms(params) {
    const out = [];
    const nums = (k) => (params.get(k) || '').split(',').filter((x) => x !== '' && Number.isFinite(+x)).map(Number);
    const one = (k) => (params.has(k) && params.get(k) !== '' ? params.get(k) : null);
    const fam = nums('motifs').filter((b) => M.motifs[b]).map((b) => motifTerm(M.motifs[b]));
    if (fam.length) out.push('motif:' + fam.join(','));
    if (one('motifQ')) out.push((params.get('motifWhere') === 'window' ? 'window:' : 'motif:') + one('motifQ'));
    const st = nums('studies').map((k) => termSlug('study', k)).filter(Boolean);
    if (st.length) out.push('study:' + st.join(','));
    [['domclass', 'domain'], ['func', 'class'], ['sub', 'subclass'], ['cond', 'condensate']].forEach(([p, key]) => {
      const v = one(p);
      if (v != null && Number.isFinite(+v) && termSlug(key, +v)) out.push(key + ':' + termSlug(key, +v));
    });
    const dis = one('disease');
    if (dis === 'als') out.push('disease:als-all');
    else if (dis != null && termSlug('disease', +dis)) out.push('disease:' + termSlug('disease', +dis));
    const binds = nums('binds').map((v) => (BIND_TERMS.find((b) => b[1] === v) || [])[0]).filter(Boolean);
    if (binds.length) out.push('binds:' + binds.join(','));
    const go = nums('go').map((v) => (v === 1 ? 'rna' : v === 2 ? 'dna' : '')).filter(Boolean);
    if (go.length) out.push('go:' + go.join(','));
    // Partners and complexes are indices into the interaction data, resolved once it loads.
    UI.legacyIx = { partner: one('partner'), cpx: one('cpx') };
    return out;
  }

  // --------------------------------------------------------------- filtering
  function sitePass(s, f) {
    if (f.states.size && !f.states.has(s.st)) return false;
    if (f.ips.size && !bits(s.ip).some((b) => f.ips.has(b))) return false;
    if (f.s2.size && !bits(s.s2).some((b) => f.s2.has(b))) return false;
    if (f.sources.size && !bits(s.src).some((b) => f.sources.has(b))) return false;
    if (f.multi && popcount(s.src) < 2) return false;
    if (f.reviewed && P[s.p].unrev) return false;
    if (f.minRep && popcount(s.study) < f.minRep) return false;
    if (f.treat) {
      if (f.treat === 'tested') { if (!s.tr.length) return false; }
      else {
        const [t, dir] = f.treat.split(':');
        if (!(s.trf & (1 << (+t * 2 + (dir === 'down' ? 0 : 1))))) return false;
      }
    }
    if (f.disorder.size && !f.disorder.has(s.dc)) return false;
    if (f.domctx.size && !f.domctx.has(s.dctx)) return false;
    if (f.nn && !(s.nn != null && s.nn <= f.nn)) return false;
    if (f.ptmDist) {
      const d = f.ptmGroup === '' ? Math.min.apply(null, s.ptmD) : s.ptmD[+f.ptmGroup];
      if (!(d <= f.ptmDist)) return false;
    }
    if (!ixSitePass(s, f)) return false;
    return true;
  }

  // condTerm: the search names a condensate, which then carries the CD-CODE confidence threshold itself.
  function proteinPass(p, f, condTerm) {
    if (f.condConf && !condTerm) { const mem = cdMemberships(p); if (!mem || !mem.some((x) => x.score >= f.condConf)) return false; }
    if (f.anyCond && !p.cond.length) return false;
    if (f.minSites > 1 && p.sites.length < f.minSites) return false;
    if (f.dense && !p.sig.length) return false;
    if (f.arch !== '' && !(p.arch && p.arch[0] === f.arch)) return false;
    if (f.argEnriched && !(p.argQ != null && p.argQ <= PADJ)) return false;
    if (!ixProteinPass(p, f)) return false;
    return true;
  }

  function filterSites(f, q) {
    const protT = q.protT || [], siteT = q.siteT || [];
    const condTerm = protT.some((t) => t.key === 'condensate');
    const okP = new Uint8Array(P.length);
    for (const p of P) okP[p.i] = (!q.proteins || q.proteins.has(p.i)) && proteinPass(p, f, condTerm) && protT.every((t) => t.vals.some((v) => v.test(p, f))) ? 1 : 0;
    const sites = [];
    const counts = new Uint32Array(P.length);
    for (const s of S) {
      if (!okP[s.p]) continue;
      if (q.positions && q.positions.has(s.p) && !q.positions.get(s.p).has(s.pos)) continue;
      if (!sitePass(s, f)) continue;
      if (siteT.length && !siteT.every((t) => t.vals.some((v) => v.test(s, f)))) continue;
      sites.push(s);
      counts[s.p]++;
    }
    return { sites, counts };
  }

  function compute() {
    const { sites, counts } = filterSites(F, Q);
    const prots = P.filter((p) => counts[p.i] > 0);
    R = { sites, prots, counts, shown: new Set(sites.map((s) => s.j)) };
    $('#count-proteins').textContent = fmt(prots.length);
    $('#count-sites').textContent = fmt(sites.length);
    $('#query-note').innerHTML = Q.note || '';
    renderScope();
    renderRailState();
    syncUrl();
  }

  // Shareable state: the search and every non-default filter live in the query string, e.g.
  // ?q=FUS&motifQ=PRPR&motifWhere=site#sites. The hash keeps the view or protein, so links reopen the same page.
  function filterParams(f) {
    const d = DEFAULT_FILTERS(), out = [];
    for (const k of Object.keys(d)) {
      const v = f[k], dv = d[k];
      if (dv instanceof Set) { if (v.size) out.push([k, Array.from(v).sort((a, b) => a - b).join(',')]); }
      else if (typeof dv === 'boolean') { if (v) out.push([k, '1']); }
      else if (String(v) !== String(dv)) out.push([k, String(v)]);
    }
    return out;
  }
  function filtersFromParams(params) {
    const f = DEFAULT_FILTERS();
    for (const k of Object.keys(f)) {
      if (!params.has(k)) continue;
      const raw = params.get(k), dv = f[k];
      if (dv instanceof Set) raw.split(',').filter((x) => x !== '' && Number.isFinite(+x)).forEach((x) => dv.add(+x));
      else if (typeof dv === 'boolean') f[k] = raw === '1' || raw === 'true';
      else if (typeof dv === 'number') { const n = +raw; if (Number.isFinite(n)) f[k] = n; }
      else f[k] = raw;
    }
    return f;
  }
  function syncUrl() {
    const qs = new URLSearchParams();
    if (Q.raw) qs.set('q', Q.raw);
    for (const [k, v] of filterParams(F)) qs.set(k, v);
    if (UI.cmpA) qs.set('A', UI.cmpA.qs);
    const s = qs.toString();
    const next = location.pathname + (s ? '?' + s : '') + location.hash;
    if (next !== location.pathname + location.search + location.hash) history.replaceState(null, '', next);
  }
  function copyLink(btn) {
    syncUrl();
    const done = (ok) => { if (!btn) return; const old = btn.textContent; btn.textContent = ok ? 'Link copied' : 'Copy failed'; setTimeout(() => { btn.textContent = old; }, 1600); };
    if (navigator.clipboard && navigator.clipboard.writeText) navigator.clipboard.writeText(location.href).then(() => done(true), () => done(false));
    else { window.prompt('Copy this link:', location.href); done(true); }
  }

  // Every active constraint, search and search terms included, as {group, label, clear, dryQ}. dryQ is the search
  // without this constraint, for the dry runs of the empty state.
  function constraints() {
    const out = [];
    const add = (group, label, clear, dryQ) => out.push({ group, label, clear, dryQ });
    if (Q.plain) add('search', 'Search: ' + (Q.plain.length > 40 ? Q.plain.slice(0, 38) + '…' : Q.plain), () => setQuery(queryText('', Q.tokens)), () => parseQuery(queryText('', Q.tokens)));
    Q.tokens.forEach((t) => {
      const rest = () => queryText(Q.plain, Q.tokens.filter((x) => !sameToken(x, t)));
      add(TERM_GROUP[t.key], tokenLabel(t), () => setQuery(rest()), () => parseQuery(rest()));
    });
    F.states.forEach((v) => add('methyl', 'State: ' + STATE_LABEL[v], () => F.states.delete(v)));
    F.ips.forEach((v) => add('methyl', ipLabel(M.maron_ips[v]) + ' IP', () => F.ips.delete(v)));
    F.s2.forEach((v) => add('methyl', 'Table S2 ' + ['Rme1', 'Rme2'][v], () => F.s2.delete(v)));
    if (F.minRep) add('prov', '≥' + F.minRep + ' independent reports', () => { F.minRep = 0; });
    F.sources.forEach((v) => add('prov', SOURCE_LABEL[M.sources[v]], () => F.sources.delete(v)));
    if (F.multi) add('prov', '≥2 source tables', () => { F.multi = false; });
    if (F.reviewed) add('prov', 'Reviewed UniProt entries only', () => { F.reviewed = false; });
    if (F.treat) add('treat', TREAT_OPTIONS.find((o) => o[0] === F.treat)[1], () => { F.treat = ''; });
    F.disorder.forEach((v) => add('struct', DISORDER_LABEL[M.disorder_classes[v]], () => F.disorder.delete(v)));
    F.domctx.forEach((v) => add('struct', DOMAIN_LABEL[M.domain_contexts[v]], () => F.domctx.delete(v)));
    if (F.nn) add('cluster', 'Methyl-R within ' + F.nn + ' aa', () => { F.nn = 0; });
    if (F.ptmDist) add('cluster', (F.ptmGroup === '' ? 'Any PTM' : PTM_LABEL[+F.ptmGroup]) + ' within ±' + F.ptmDist + ' aa', () => { F.ptmDist = 0; });
    if (F.dense) add('cluster', 'Dense methyl-R region', () => { F.dense = false; });
    if (F.arch !== '') add('cluster', F.arch, () => { F.arch = ''; });
    if (F.anyCond) add('protein', 'In a CD-CODE condensate', () => { F.anyCond = false; });
    if (F.condConf && !termValues('condensate').length) add('protein', 'CD-CODE confidence ≥ ' + F.condConf + '/5', () => { F.condConf = 0; });
    if (F.minSites > 1) add('protein', '≥' + F.minSites + ' sites per protein', () => { F.minSites = 1; });
    if (F.argEnriched) add('protein', 'Arg-methylation enriched', () => { F.argEnriched = false; });
    ixConstraints(add);
    return out;
  }

  function clearAll() {
    F = DEFAULT_FILTERS();
    setQuery('');
    afterFilterChange();
  }
  // Filters and search terms go; the plain words of the search stay.
  function clearFilters() {
    F = DEFAULT_FILTERS();
    setQuery(Q.plain || '');
    afterFilterChange();
  }

  function afterFilterChange() {
    UI.pageP = 50; UI.pageS = 100;
    renderRail();
    compute();
    renderView();
  }

  function setFilters(mutator, view) {
    mutator(F);
    UI.pageP = 50; UI.pageS = 100;
    compute();
    renderRail();
    if (view && view !== UI.view) go(view); else renderView();
  }

  // The scope bar sits under the search box on every view and lists everything that narrows the results.
  function renderScope() {
    const box = $('#scope');
    if (!box) return;
    const cons = constraints();
    box.classList.toggle('on', cons.length > 0);
    const cmp = UI.cmpA ? '<button class="btn small quiet scope-cmp" type="button" id="scope-cmp" title="A: ' + esc(UI.cmpA.label || '') + '">Compare with A</button>' : '';
    const bindCmp = () => { const b = $('#scope-cmp', box); if (b) b.addEventListener('click', () => { go('compare'); window.scrollTo({ top: 0 }); }); };
    if (!cons.length) {
      box.innerHTML = '<span class="scope-count">Showing all <b>' + fmt(S.length) + '</b> sites on <b>' + fmt(P.length) + '</b> proteins</span><span class="scope-hint">Search or use the filters to narrow this down</span>' + cmp;
      bindCmp();
      return;
    }
    box.innerHTML = '<span class="scope-count"><span class="scope-flag">Filtered</span><b>' + fmt(R.sites.length) + '</b> of ' + fmt(S.length) + ' sites · <b>' + fmt(R.prots.length) + '</b> proteins</span>' +
      '<span class="scope-chips">' + cons.map((c, k) => '<span class="af' + (c.group === 'search' ? ' search' : '') + '">' + esc(c.label) + '<button type="button" data-k="' + k + '" aria-label="Remove ' + esc(c.label) + '">×</button></span>').join('') + '</span>' +
      cmp + '<button class="btn small quiet scope-link" type="button" id="scope-link" title="Copy a link that reopens this page with the same search and filters">Copy link</button>' +
      '<button class="btn small scope-clear" type="button" id="scope-clear">Clear all</button>';
    $$('.af button', box).forEach((b) => b.addEventListener('click', () => { cons[+b.dataset.k].clear(); afterFilterChange(); }));
    $('#scope-link', box).addEventListener('click', (e) => copyLink(e.currentTarget));
    bindCmp();
    $('#scope-clear', box).addEventListener('click', clearAll);
    syncToolbar();
  }

  // --------------------------------------------------------------- routing
  function go(view, acc) {
    const hash = view === 'protein' ? acc : view;
    if (location.hash.slice(1) === hash) route();
    else location.hash = hash;
  }

  function route() {
    const h = decodeURIComponent(location.hash.slice(1));
    let view = 'overview';
    if (['overview', 'proteins', 'sites', 'interactions', 'figures', 'about', 'report', 'compare'].includes(h)) view = h;
    else if (h === 'docs') view = 'about';
    else if (h) {
      // #P35637 opens a protein; #P35637/R216 opens it with that site selected (shown numbering, so histones are
      // numbered without the initiator Met).
      const [accPart, sitePart] = h.split('/');
      const hit = lookupExact(accPart);
      if (hit && hit.length) {
        view = 'protein';
        const site = /^R?(\d+)$/i.exec(sitePart || '');
        if (UI.detail !== hit[0]) {
          UI.detail = hit[0]; UI.zoom = null; UI.openEvidence.clear(); UI.pIx = '';
          // A link such as ?q=FUS R216#P35637 opens with that site selected.
          const want = Q.positions && Q.positions.get(hit[0]);
          if (UI.highlight == null && want && want.size === 1) UI.highlight = want.values().next().value;
        }
        if (site) UI.highlight = +site[1] + (P[hit[0]].hist ? 1 : 0);
        rememberProtein(hit[0]);
      }
    }
    UI.view = view;
    if (!['protein', 'about'].includes(view)) UI.lastView = view;
    setDocTitle();
    if (view !== 'report' && view !== 'compare') setPrintPageStyle(false);
    $$('.tab').forEach((t) => t.setAttribute('aria-selected', String(t.dataset.view === view)));
    const ptab = $('#tab-protein');
    ptab.hidden = UI.detail == null;
    if (UI.detail != null) ptab.textContent = P[UI.detail].label;
    tabsInView();
    $('#layout').classList.toggle('no-rail', view === 'about' || view === 'report' || view === 'compare');
    renderView();
    if (view !== 'protein') UI.highlight = null;
  }

  // Narrow screens scroll the tab bar: keep the selected tab visible and fade the edge while tabs are hidden.
  function tabsInView() {
    const bar = $('.tabs');
    if (!bar) return;
    const sel = $('.tab[aria-selected="true"]', bar);
    if (sel && (sel.offsetLeft < bar.scrollLeft || sel.offsetLeft + sel.offsetWidth > bar.scrollLeft + bar.clientWidth)) bar.scrollLeft = Math.max(0, sel.offsetLeft - 16);
    tabsFade();
  }
  function tabsFade() {
    const bar = $('.tabs');
    if (bar) bar.classList.toggle('more-right', bar.scrollWidth - bar.clientWidth - bar.scrollLeft > 4);
  }

  function renderView() {
    const el = $('#content');
    hideTip();
    if (UI.view === 'overview') renderOverview(el);
    else if (UI.view === 'proteins') renderProteins(el);
    else if (UI.view === 'sites') renderSites(el);
    else if (UI.view === 'interactions') renderInteractions(el);
    else if (UI.view === 'figures') renderFigures(el);
    else if (UI.view === 'report') renderReport(el);
    else if (UI.view === 'compare') renderCompare(el);
    else if (UI.view === 'about') renderAbout(el);
    else if (UI.view === 'protein') renderProtein(el, P[UI.detail]);
  }

  function openProtein(i, pos) {
    UI.highlight = pos == null ? null : pos;
    if (UI.detail !== i) { UI.detail = i; UI.zoom = null; UI.openEvidence.clear(); UI.pIx = ''; }
    go('protein', P[i].acc + (pos == null ? '' : '/' + rl(pos, P[i])));
    window.scrollTo({ top: 0 });
  }
  // The open site is part of the address (#ACC/R123) without adding a history step for every click.
  function setSiteHash(p, pos) {
    const want = '#' + p.acc + (pos == null ? '' : '/' + rl(pos, p));
    if (location.hash !== want) history.replaceState(null, '', location.pathname + location.search + want);
    setDocTitle();
  }
  const VIEW_NAMES = { overview: 'Overview', proteins: 'Proteins', sites: 'Sites', interactions: 'Interactions', figures: 'Figures', about: 'Documentation', report: 'Insight report', compare: 'Comparison' };
  let baseTitle = '';
  function setDocTitle() {
    baseTitle = baseTitle || document.title;
    const p = UI.view === 'protein' && UI.detail != null ? P[UI.detail] : null;
    document.title = (p ? p.label + (UI.highlight != null ? ' ' + rl(UI.highlight, p) : '') : VIEW_NAMES[UI.view] || '') + ' · ' + baseTitle.split(' · ')[0];
  }
  // Proteins opened in this session, newest first, for the breadcrumb bar.
  function rememberProtein(i) {
    UI.recent = [i].concat((UI.recent || []).filter((x) => x !== i)).slice(0, 6);
  }
  // Breadcrumbs on protein pages: where the reader came from (view, with the selection it shows), this protein and the
  // open site, then the proteins seen recently.
  function crumbsHtml(p) {
    const from = UI.lastView || 'proteins', cons = constraints();
    const sel = cons.length ? fmt(R.sites.length) + ' site' + (R.sites.length === 1 ? '' : 's') + ' · ' + (selectionText().length > 60 ? selectionText().slice(0, 58) + '…' : selectionText()) : 'all sites';
    const recent = (UI.recent || []).filter((i) => i !== p.i).slice(0, 5);
    return '<nav class="crumbs" aria-label="Where you are"><ol>' +
      (from !== 'overview' ? '<li><button class="linkish" type="button" data-crumb="overview">Overview</button></li>' : '') +
      '<li><button class="linkish" type="button" data-crumb="' + from + '">' + esc(VIEW_NAMES[from] || 'Results') + '</button> <span class="muted small">(' + esc(sel) + ')</span></li>' +
      '<li aria-current="page"><b>' + esc(p.label) + '</b></li>' +
      '<li id="crumb-site"' + (UI.highlight == null ? ' hidden' : '') + '>' + (UI.highlight == null ? '' : esc(rl(UI.highlight, p))) + '</li></ol>' +
      (recent.length ? '<span class="crumb-recent"><span class="muted small">Recent</span> ' + recent.map((i) => '<button class="linkish" type="button" data-open-p="' + i + '">' + esc(P[i].label) + '</button>').join(' · ') + '</span>' : '') + '</nav>';
  }

  // ---------------------------------------------------------------- hero
  function renderHero() {
    const multi = S.filter((s) => popcount(s.study) >= 2).length;
    const inCond = P.filter((p) => p.cond.length).length;
    const measured = S.filter((s) => s.tr.length).length;
    const dense = P.filter((p) => p.sig.length).length;
    const stats = [
      [fmt(M.n_sites), 'methylarginine sites', () => { F = DEFAULT_FILTERS(); setQuery(''); setFilters(() => {}, 'sites'); }],
      [fmt(M.n_proteins), 'proteins', () => { F = DEFAULT_FILTERS(); setQuery(''); setFilters(() => {}, 'proteins'); }],
      [fmt(multi), 'sites with ≥2 independent reports', () => setFilters((f) => { f.minRep = 2; }, 'sites')],
      [fmt(measured), 'sites with PRMTi data', () => setFilters((f) => { f.treat = 'tested'; }, 'sites')],
      [fmt(inCond), 'proteins in condensates', () => setFilters((f) => { f.anyCond = true; }, 'proteins')],
      [fmt(dense), 'proteins with dense clusters', () => setFilters((f) => { f.dense = true; }, 'proteins')],
    ];
    const box = $('#hero-stats');
    box.innerHTML = stats.map((s, k) => '<button class="stat" type="button" data-k="' + k + '"><b>' + s[0] + '</b><span>' + s[1] + '</span></button>').join('');
    $$('.stat', box).forEach((b) => b.addEventListener('click', () => { stats[+b.dataset.k][2](); document.getElementById('toolbar').scrollIntoView({ behavior: 'smooth' }); }));

    // Decorative lollipop silhouette of a real substrate (FUS when present).
    const hit = lookupExact('FUS');
    const p = hit ? P[hit[0]] : P.slice().sort((a, b) => b.sites.length - a.sites.length)[0];
    const len = p.len || p.maxPos * 1.05;
    const W = 600, H = 160, base = 150;
    const x = (pos) => 6 + (pos / len) * (W - 12);
    const levels = [];
    let marks = '';
    for (const j of p.sites) {
      const cx = x(S[j].pos);
      let lv = 0;
      while (levels[lv] != null && cx - levels[lv] < 11 && lv < 8) lv++;
      levels[lv] = cx;
      const cy = base - 24 - lv * 12;
      marks += '<line x1="' + cx.toFixed(1) + '" y1="' + base + '" x2="' + cx.toFixed(1) + '" y2="' + cy + '" stroke="#fff" stroke-width="1.2"/>' +
        '<circle cx="' + cx.toFixed(1) + '" cy="' + cy + '" r="4.5" fill="#fff"/>';
    }
    $('#hero-art').innerHTML = '<svg viewBox="0 0 ' + W + ' ' + H + '" preserveAspectRatio="xMaxYMax meet"><rect x="6" y="' + (base - 3) + '" width="' + (W - 12) + '" height="6" rx="3" fill="#fff"/>' + marks + '</svg>';
    $('#hero-art-caption').textContent = p.label + ' · ' + p.sites.length + ' methylated arginines';
  }

  // ------------------------------------------------------------------ rail
  const SECTIONS = [
    ['motif', 'Sequence motif', true], ['methyl', 'Methylation evidence', true], ['prov', 'Provenance', false], ['treat', 'PRMT inhibitors', false], ['struct', 'Disorder & domains', false], ['cluster', 'Clustering & crosstalk', false],
    ['protein', 'Protein features', false], ['ix', 'Interactions & binding', false],
  ];

  function chipRow(name, items, selected, dots) {
    return '<div class="chips" data-set="' + name + '">' + items.map(([v, label, help]) =>
      '<button class="chip" type="button" data-v="' + v + '" aria-pressed="' + selected.has(v) + '"' + (help ? ' title="' + esc(help) + '"' : '') + '>' +
      (dots ? '<i class="dot st-' + v + '"></i>' : '') + esc(label) + '</button>').join('') + '</div>';
  }

  function options(list, current, first) {
    return (first ? '<option value="">' + esc(first) + '</option>' : '') +
      list.map(([v, label]) => '<option value="' + esc(v) + '"' + (String(current) === String(v) ? ' selected' : '') + '>' + esc(label) + '</option>').join('');
  }

  function subclassOptions(current) {
    const groups = {};
    M.func_subclasses.forEach((name, k) => { const parent = M.subParent[name] || 'Other'; (groups[parent] = groups[parent] || []).push([termSlug('subclass', k), name]); });
    return Object.keys(groups).sort().map((g) => '<optgroup label="' + esc(g) + '">' +
      groups[g].sort((a, b) => a[1].localeCompare(b[1])).map(([v, n]) => '<option value="' + esc(v) + '"' + (current === v ? ' selected' : '') + '>' + esc(n) + '</option>').join('') + '</optgroup>').join('');
  }

  // Rail controls that write search terms: chips toggle one value of a key (any of them), selects set the key.
  function termChips(key, items) {
    const on = new Set(termValues(key).map((v) => v.text.toLowerCase()));
    return '<div class="chips" data-term="' + key + '">' + items.map(([v, label, help]) =>
      '<button class="chip" type="button" data-v="' + esc(v) + '" aria-pressed="' + on.has(String(v).toLowerCase()) + '" title="' + esc((help ? help + ' · ' : '') + 'search ' + key + ':' + v) + '">' + esc(label) + '</button>').join('') + '</div>';
  }
  const CHIP_OF = { motif: (v) => v.family, binds: (v) => v.bindClass };
  function termSelect(id, key, html, label) {
    const n = termValues(key).length, one = singleTerm(key);
    return '<select id="' + id + '" data-termsel="' + key + '"' + (n ? ' class="set"' : '') + (label ? ' aria-label="' + esc(label) + '"' : '') + '>' +
      (n > 1 ? '<option value="__multi" selected>Several (see the search)</option>' : '') + '<option value="">Any</option>' + html.replace(one ? 'value="' + esc(one.text) + '"' : '\u0000', (m) => m + ' selected') + '</select>';
  }
  const termOptions = (list) => list.map(([v, label]) => '<option value="' + esc(v) + '">' + esc(label) + '</option>').join('');
  // The custom-pattern term the motif box edits: a window: term, or a motif: term that names no family.
  const patternToken = () => Q.tokens.find((t) => t.key === 'window' || (t.key === 'motif' && !t.vals.some((v) => v.family)));
  const partnerToken = () => Q.tokens.find((t) => t.key === 'binds' && !t.vals.some((v) => v.bindClass));
  // Replace one term of the search (or drop it when text is empty), keeping everything else.
  function replaceToken(old, key, text, view) {
    const toks = Q.tokens.filter((t) => !old || !sameToken(t, old));
    if (text) toks.push({ key, vals: text.split(',').filter(Boolean).map((x) => ({ text: x })) });
    applyQueryText(queryText(Q.plain, toks), view);
  }

  function sectionBody(key) {
    const condOpts = M.condensates.map((c, k) => [k, c[0] + ' (' + c[4] + ')']).filter(([k]) => M.condensates[k][4] > 0).sort((a, b) => a[1].localeCompare(b[1]));
    const archClasses = Array.from(new Set(P.filter((p) => p.arch).map((p) => p.arch[0]))).sort();
    const sel = (id, html, isSet, label) => '<select id="' + id + '"' + (isSet ? ' class="set"' : '') + (label ? ' aria-label="' + esc(label) + '"' : '') + '>' + html + '</select>';
    switch (key) {
      case 'methyl': return '<div class="fgroup"><span class="flabel" title="Consensus Rme1/Rme2 calls across sources; most compiled sites have no state call">Methylation state</span>' +
        chipRow('states', [[1, 'Rme1'], [2, 'Rme2'], [3, 'Rme1 + Rme2'], [0, 'No call']], F.states, true) + '</div>' +
        '<div class="fgroup"><span class="flabel" title="Maron et al. 2021 Table S3: methyl-specific immunoprecipitation that recovered the peptide">Maron methyl-IP</span>' +
        chipRow('ips', M.maron_ips.map((n, k) => [k, ipLabel(n), 'Recovered by the ' + (IP_HELP[n] || n) + ' IP']), F.ips) + '</div>' +
        '<div class="fgroup"><span class="flabel" title="Maron et al. 2021 Table S2 Rme1/Rme2 site calls">Maron Table S2 call</span>' + chipRow('s2', [[0, 'Rme1'], [1, 'Rme2']], F.s2) + '</div>';
      case 'prov': return '<div class="fgroup"><span class="flabel">Reported by study or database (match any)</span>' +
        termChips('study', M.studyOrder.map((k) => [termSlug('study', k), M.studies[k][0] + ' · ' + compact(M.studies[k][2]), M.studies[k][1] || STUDY_HELP[M.studies[k][0]] || ''])) + '</div>' +
        '<div class="fgroup"><label class="flabel" for="f-minrep">Independent reports per site</label>' + sel('f-minrep', options([[2, '≥2 studies or databases'], [3, '≥3'], [5, '≥5']], F.minRep || '', 'Any'), !!F.minRep) + '</div>' +
        '<div class="fgroup"><span class="flabel">Pipeline source table</span>' + chipRow('sources', M.sources.map((n, k) => [k, SOURCE_LABEL[n], SOURCE_HELP[n]]), F.sources) +
        '<label class="check"><input type="checkbox" id="f-multi"' + (F.multi ? ' checked' : '') + '> In ≥2 source tables at the same position</label>' +
        '<label class="check" title="' + esc(unrevHelp()) + '"><input type="checkbox" id="f-reviewed"' + (F.reviewed ? ' checked' : '') + '> Reviewed UniProt entries only</label></div>';
      case 'treat': return '<div class="fgroup"><label class="flabel" for="f-treat">Maron Table S3 response (adj. p ≤ 0.05)</label>' + sel('f-treat', options(TREAT_OPTIONS, F.treat), !!F.treat) + '</div>';
      case 'motif': {
        const pt = patternToken();
        return '<div class="fgroup"><label class="flabel" for="f-motifq">Find a motif (2–10 residues)</label>' +
          '<div class="motif-search"><input id="f-motifq" type="text" maxlength="60" autocomplete="off" spellcheck="false" placeholder="e.g. GXRG, [RK]XR, PRPR" value="' + esc(pt ? pt.vals.map((v) => v.text).join(',') : '') + '"' + (pt ? ' class="set"' : '') + ' aria-describedby="f-motifmsg">' +
          '<button class="btn small" type="button" id="f-motifgo">Search</button></div>' +
          sel('f-motifwhere', options([['motif', 'Covering the methyl-R'], ['window', 'Anywhere within ±15 aa']], pt ? pt.key : 'motif'), !!pt && pt.key === 'window', 'Where the motif may sit') +
          '<p class="note" id="f-motifmsg">' + motifHint() + '</p></div>' +
          '<div class="fgroup"><span class="flabel">Motif family at the methyl-R (match any)</span>' + termChips('motif', M.motifs.map((m) => [motifTerm(m), MOTIF_LABEL[m], MOTIF_HELP[m]])) + '</div>';
      }
      case 'struct': return '<div class="fgroup"><span class="flabel">Disorder (MobiDB)</span>' + chipRow('disorder', M.disorder_classes.map((c, k) => [k, DISORDER_LABEL[c] || c, DISORDER_HELP[c]]), F.disorder) + '</div>' +
        '<div class="fgroup"><span class="flabel">Domain context (InterPro)</span>' + chipRow('domctx', M.domain_contexts.map((c, k) => [k, DOMAIN_LABEL[c] || c, DOMAIN_HELP[c]]), F.domctx) + '</div>' +
        '<div class="fgroup"><label class="flabel" for="f-domclass">Nearest domain class (lab domain ontology)</label>' + termSelect('f-domclass', 'domain', termOptions(M.domain_classes.map((c, k) => [termSlug('domain', k), c]).sort((a, b) => a[1].localeCompare(b[1])))) + '</div>';
      case 'cluster': return '<div class="fgroup"><label class="flabel" for="f-nn">Another methyl-R within</label>' + sel('f-nn', options([[2, '2 aa'], [3, '3 aa'], [5, '5 aa'], [10, '10 aa'], [20, '20 aa']], F.nn || '', 'Any spacing'), !!F.nn) + '</div>' +
        '<div class="fgroup"><label class="flabel" for="f-ptm" title="Sites from the pipeline\'s cross-PTM master table (iPTMnet sources and UniProt). Not exhaustive: sites reported elsewhere may be missing.">Other PTM nearby (PTM master table; not exhaustive)</label>' + sel('f-ptm', options(PTM_LABEL.map((l, k) => [k, l]), F.ptmGroup, 'Any PTM type'), !!F.ptmDist) +
        sel('f-ptmdist', options([[2, 'within ±2 aa'], [5, 'within ±5 aa'], [10, 'within ±10 aa'], [20, 'within ±20 aa']], F.ptmDist || '', 'Off'), !!F.ptmDist, 'PTM distance') + '</div>' +
        '<div class="fgroup"><label class="check"><input type="checkbox" id="f-dense"' + (F.dense ? ' checked' : '') + '> Protein has a significant dense methyl-R region (multiscale null, q ≤ 0.05)</label></div>' +
        '<div class="fgroup"><label class="flabel" for="f-arch">Clustered architecture (≥3 sites in 20 aa)</label>' + sel('f-arch', options(archClasses.map((a) => [a, a]), F.arch, 'Any'), F.arch !== '') + '</div>';
      case 'protein': return '<div class="fgroup"><label class="flabel" for="f-func">Functional class (lab protein ontology)</label>' + termSelect('f-func', 'class', termOptions(M.func_classes.map((c, k) => [termSlug('class', k), c]).sort((a, b) => a[1].localeCompare(b[1])))) + '</div>' +
        '<div class="fgroup"><label class="flabel" for="f-sub">Functional subclass</label>' + termSelect('f-sub', 'subclass', subclassOptions('')) + '</div>' +
        '<div class="fgroup"><label class="flabel" for="f-disease">Disease family (DisGeNET)</label>' + termSelect('f-disease', 'disease', termOptions(termNames('disease').map((x) => [x.slug, x.name]))) + '</div>' +
        '<div class="fgroup"><label class="flabel" for="f-cond">CD-CODE condensate (methylated members)</label>' + termSelect('f-cond', 'condensate', termOptions(condOpts.map(([k, l]) => [termSlug('condensate', k), l]))) +
        '<label class="check"><input type="checkbox" id="f-anycond"' + (F.anyCond ? ' checked' : '') + '> In any CD-CODE condensate</label>' +
        '<label class="flabel" for="f-condconf" title="CD-CODE\'s 1–5 confidence that the protein belongs to the condensate">Membership confidence (CD-CODE)</label>' +
        sel('f-condconf', options([[2, '≥ 2 of 5'], [3, '≥ 3 of 5'], [4, '≥ 4 of 5'], [5, '5 of 5']], F.condConf || '', 'Any'), !!F.condConf) + '</div>' +
        '<div class="fgroup"><label class="flabel" for="f-min">Minimum methyl-R sites per protein</label>' +
        '<input id="f-min" type="number" min="1" max="100" step="1" value="' + F.minSites + '"' + (F.minSites > 1 ? ' class="set"' : '') + '></div>' +
        '<div class="fgroup"><label class="check"><input type="checkbox" id="f-arg"' + (F.argEnriched ? ' checked' : '') + '> Arginines methylated more than the proteome rate (Fisher q ≤ 0.05)</label></div>';
      case 'ix': return ixRailBody(sel);
      default: return '';
    }
  }

  function renderRail() {
    const rail = $('#rail');
    const open = new Set($$('details[open]', rail).map((d) => d.dataset.sec));
    const firstRender = !rail.children.length;
    const n = constraints().filter((c) => c.group !== 'search').length;
    rail.innerHTML =
      '<div class="rail-head"><h2>Filters</h2><div class="actions">' +
      '<button class="btn small" type="button" id="reset-filters"' + (n ? '' : ' disabled') + '>Clear filters' + (n ? ' (' + n + ')' : '') + '</button>' +
      '<button class="btn quiet small rail-toggle" type="button" id="rail-toggle" aria-expanded="' + (rail.dataset.collapsed === 'false') + '">' + (rail.dataset.collapsed === 'false' ? 'Hide' : 'Show') + '</button></div></div>' +
      SECTIONS.map(([key, title, dflt]) => '<details data-sec="' + key + '"' + ((firstRender ? dflt : open.has(key)) ? ' open' : '') + '><summary><span>' + esc(title) + '</span>' +
        '<span class="sec-count" data-count="' + key + '" hidden></span><span class="sec-clear" role="button" tabindex="0" data-clear="' + key + '" hidden>Clear</span></summary>' +
        sectionBody(key) + '</details>').join('');
    bindRail();
    renderRailState();
  }

  function renderRailState() {
    // Rail controls mirror the search terms (typed, removed from a chip, or restored from a link) unless being edited.
    $$('#rail .chips[data-term]').forEach((row) => {
      const on = new Set(termValues(row.dataset.term).map((v) => v.text.toLowerCase()));
      $$('.chip', row).forEach((b) => b.setAttribute('aria-pressed', String(on.has(b.dataset.v.toLowerCase()))));
    });
    $$('#rail select[data-termsel]').forEach((el) => {
      if (document.activeElement === el) return;
      const n = termValues(el.dataset.termsel).length, one = singleTerm(el.dataset.termsel);
      if (n > 1) return;
      const want = one ? one.text : '';
      if (el.value !== want && Array.from(el.options).some((o) => o.value === want)) el.value = want;
      el.classList.toggle('set', !!n);
    });
    const mq = $('#f-motifq'), pt = patternToken();
    if (mq && document.activeElement !== mq) {
      const want = pt ? pt.vals.map((v) => v.text).join(',') : '';
      if (mq.value.toUpperCase().replace(/\s+/g, '') !== want) mq.value = want;
      mq.classList.toggle('set', !!pt);
      const mw = $('#f-motifwhere');
      if (mw && pt) mw.value = pt.key;
    }
    const pq = $('#f-partner');
    if (pq && document.activeElement !== pq && ixReady()) {
      const tk = partnerToken(), want = tk ? tk.vals.map((v) => v.text).join(',') : '';
      if (pq.value !== want) { pq.value = want; pq.classList.toggle('set', !!want); const pm = $('#f-partnermsg'); if (pm) pm.textContent = partnerHint(); }
    }
    const counts = {};
    for (const c of constraints()) counts[c.group] = (counts[c.group] || 0) + 1;
    $$('[data-count]').forEach((el) => { const n = counts[el.dataset.count] || 0; el.hidden = !n; el.textContent = n; });
    $$('[data-clear]').forEach((el) => { el.hidden = !counts[el.dataset.clear]; });
    $$('#rail details').forEach((d) => d.classList.toggle('active', !!counts[d.dataset.sec]));
    const reset = $('#reset-filters');
    if (reset) {
      const n = Object.keys(counts).filter((g) => g !== 'search').reduce((a, g) => a + counts[g], 0);
      reset.disabled = !n;
      reset.textContent = 'Clear filters' + (n ? ' (' + n + ')' : '');
    }
  }

  function bindRail() {
    const rail = $('#rail');
    $$('.chips[data-set]', rail).forEach((row) => {
      row.addEventListener('click', (e) => {
        const b = e.target.closest('.chip');
        if (!b) return;
        const set = F[row.dataset.set];
        const v = +b.dataset.v;
        if (set.has(v)) set.delete(v); else set.add(v);
        b.setAttribute('aria-pressed', String(set.has(v)));
        UI.pageP = 50; UI.pageS = 100;
        compute(); renderView();
      });
    });
    const bind = (id, fn, ev = 'change') => {
      const el = $('#' + id, rail);
      if (el) el.addEventListener(ev, () => { fn(el); UI.pageP = 50; UI.pageS = 100; compute(); renderView(); if (el.tagName === 'SELECT' || el.type === 'number') el.classList.toggle('set', !!el.value && el.value !== '1'); });
    };
    bind('f-multi', (el) => { F.multi = el.checked; });
    bind('f-reviewed', (el) => { F.reviewed = el.checked; });
    bind('f-minrep', (el) => { F.minRep = +el.value || 0; });
    bind('f-treat', (el) => { F.treat = el.value; });
    bind('f-nn', (el) => { F.nn = +el.value || 0; });
    bind('f-ptm', (el) => { F.ptmGroup = el.value; if (!F.ptmDist) { F.ptmDist = 5; $('#f-ptmdist').value = '5'; $('#f-ptmdist').classList.add('set'); } });
    bind('f-ptmdist', (el) => { F.ptmDist = +el.value || 0; });
    bind('f-dense', (el) => { F.dense = el.checked; });
    bind('f-arch', (el) => { F.arch = el.value; });
    bind('f-anycond', (el) => { F.anyCond = el.checked; });
    bind('f-condconf', (el) => { F.condConf = +el.value || 0; });
    bind('f-min', (el) => { F.minSites = Math.max(1, Math.min(100, +el.value || 1)); }, 'input');
    bind('f-arg', (el) => { F.argEnriched = el.checked; });
    $$('.chips[data-term]', rail).forEach((row) => row.addEventListener('click', (e) => {
      const b = e.target.closest('.chip');
      if (b) editTerms(row.dataset.term, b.dataset.v, 'toggle', null, CHIP_OF[row.dataset.term]);
    }));
    $$('select[data-termsel]', rail).forEach((el) => el.addEventListener('change', () => { if (el.value !== '__multi') setTerm(el.dataset.termsel, el.value); }));
    const mq = $('#f-motifq', rail);
    if (mq) bindMotifBox(mq, $('#f-motifwhere', rail), $('#f-motifgo', rail), $('#f-motifmsg', rail));
    bindIxRail(rail, () => { UI.pageP = 50; UI.pageS = 100; compute(); renderView(); });
    $('#reset-filters').addEventListener('click', clearFilters);
    $$('[data-clear]', rail).forEach((el) => {
      const run = (e) => {
        e.preventDefault(); e.stopPropagation();
        for (const c of constraints()) if (c.group === el.dataset.clear) c.clear();
        afterFilterChange();
      };
      el.addEventListener('click', run);
      el.addEventListener('keydown', (e) => { if (e.key === 'Enter' || e.key === ' ') run(e); });
    });
    $('#rail-toggle').addEventListener('click', () => {
      const collapsed = rail.dataset.collapsed !== 'false';
      rail.dataset.collapsed = collapsed ? 'false' : 'true';
      const t = $('#rail-toggle');
      t.textContent = collapsed ? 'Hide' : 'Show';
      t.setAttribute('aria-expanded', String(collapsed));
    });
  }

  // --------------------------------------------------------- global events
  // Replace the plain words of the search and keep its terms.
  const setPlain = (text) => setQuery(queryText(text, Q.tokens));
  function setQuery(text) {
    $('#q').value = text;
    Q = parseQuery(text);
    $('#q-clear').hidden = !text;
    const url = new URL(location.href);
    if (text.trim()) url.searchParams.set('q', text.trim()); else url.searchParams.delete('q');
    history.replaceState(null, '', url.pathname + url.search + url.hash);
  }

  function syncToolbar() {
    const tb = $('#toolbar');
    if (tb) document.documentElement.style.setProperty('--toolbar-h', tb.offsetHeight + 'px');
  }

  function bindGlobal() {
    const input = $('#q');
    const list = $('#suggest');
    let sel = -1;
    let items = [];
    const closeList = () => { list.hidden = true; input.setAttribute('aria-expanded', 'false'); sel = -1; };
    // The word being typed: key:value terms get value suggestions, a lone plain word gets protein suggestions.
    const lastWord = () => { const m = /(?:^|\s)((?:[A-Za-z]+:"[^"]*"?)|[^\s]*)$/.exec(input.value); return m ? m[1] : ''; };
    const drawList = () => {
      const w = lastWord();
      const terms = termSuggestions(w);
      const plainWords = splitTerms(input.value).plain.split(/[\s,;]+/).filter(Boolean);
      const prots = !/:/.test(w) && plainWords.length <= 1 ? suggestions(w) : [];
      items = terms.map((t) => ({ term: t })).concat(prots.map(([p, alias]) => ({ p, alias }))).slice(0, 10);
      if (!items.length || !w) { closeList(); return; }
      list.innerHTML = items.map((it, k) => '<li role="option" id="sg-' + k + '" aria-selected="' + (k === sel) + '" data-k="' + k + '">' + (it.term
        ? '<b class="sg-term">' + esc(it.term.label) + '</b><span>' + esc(it.term.sub || '') + '</span><em>filter</em>'
        : '<b>' + esc(it.p.label) + '</b><span>' + (it.alias ? 'synonym ' + esc(it.alias) + ' · ' : '') + esc(it.p.name || it.p.acc) + '</span><em>' + it.p.sites.length + ' site' + (it.p.sites.length === 1 ? '' : 's') + '</em>') + '</li>').join('');
      list.hidden = false;
      input.setAttribute('aria-expanded', 'true');
    };
    // A search made only of terms lists sites (or proteins for protein terms); plain words list proteins.
    const landing = () => (!Q.plain && Q.tokens.length ? termView() : 'proteins');
    const apply = debounce(() => {
      setQuery(input.value);
      UI.pageP = 50; UI.pageS = 100;
      compute();
      if (UI.view === 'overview' && Q.raw && !lastWord().endsWith(':')) go(landing());
      else renderView();
    }, 220);
    input.addEventListener('input', () => { $('#q-clear').hidden = !input.value; sel = -1; drawList(); apply(); });
    input.addEventListener('keydown', (e) => {
      if (e.key === 'ArrowDown' && !list.hidden) { sel = Math.min(items.length - 1, sel + 1); drawList(); e.preventDefault(); }
      else if (e.key === 'ArrowUp' && !list.hidden) { sel = Math.max(-1, sel - 1); drawList(); e.preventDefault(); }
      else if (e.key === 'Escape') { if (!list.hidden) closeList(); else if (input.value) { setQuery(''); compute(); renderView(); } }
      else if (e.key === 'Enter') {
        e.preventDefault();
        if (sel >= 0 && items[sel]) { pick(items[sel]); return; }
        closeList();
        setQuery(input.value);
        compute();
        if (R.prots.length === 1) openProtein(R.prots[0].i, Q.positions && Q.positions.size ? Array.from(Q.positions.values())[0].values().next().value : null);
        // One protein named with a site it lacks: open the protein; the note under the search says the site is not in the set.
        else if (!R.prots.length && Q.proteins && Q.proteins.size === 1 && Q.positions) openProtein(Q.proteins.values().next().value);
        else if (Q.positions) go('sites');
        else go(landing());
      }
    });
    input.addEventListener('blur', () => setTimeout(closeList, 150));
    list.addEventListener('mousedown', (e) => {
      const li = e.target.closest('li');
      if (li) { e.preventDefault(); pick(items[+li.dataset.k]); }
    });
    const pick = (it) => {
      const w = lastWord(), head = input.value.slice(0, input.value.length - w.length);
      if (it.term) {
        // A key keeps the list open for its values; a value completes the term and applies it.
        input.value = head + it.term.term + (it.term.open ? '' : ' ');
        sel = -1;
        input.focus();
        if (it.term.open) { drawList(); return; }
        closeList();
        setQuery(input.value);
        UI.pageP = 50; UI.pageS = 100;
        compute();
        if (UI.view === 'overview') go(landing()); else renderView();
        return;
      }
      closeList();
      setQuery((head + (it.p.gene || it.p.acc)).trim());
      compute();
      openProtein(it.p.i);
    };
    $('#q-clear').addEventListener('click', () => { setQuery(''); compute(); renderView(); input.focus(); });
    $$('.tab').forEach((t) => t.addEventListener('click', () => {
      if (t.dataset.view === 'protein') { if (UI.detail != null) openProtein(UI.detail); }
      else go(t.dataset.view);
    }));
    window.addEventListener('hashchange', route);
    const home = $('.home-link');
    if (home) home.addEventListener('click', (e) => { e.preventDefault(); go('overview'); window.scrollTo({ top: 0, behavior: 'smooth' }); });
    $$('[data-theme-set]').forEach((b) => b.addEventListener('click', () => { setTheme(b.dataset.themeSet); if (UI.view === 'protein') renderView(); }));
    setTheme(store.get('mre-theme') || 'system', true);
    syncToolbar();
    window.addEventListener('resize', debounce(syncToolbar, 150));
    window.addEventListener('resize', debounce(tabsFade, 150));
    $('.tabs').addEventListener('scroll', debounce(tabsFade, 60), { passive: true });
  }

  function setTheme(mode, initial) {
    if (mode === 'light' || mode === 'dark') document.documentElement.setAttribute('data-theme', mode);
    else document.documentElement.removeAttribute('data-theme');
    if (!initial) store.set('mre-theme', mode === 'system' ? null : mode);
    $$('[data-theme-set]').forEach((b) => b.setAttribute('aria-pressed', String(b.dataset.themeSet === mode)));
  }

  // ------------------------------------------------------------- rendering
  // Custom motif search on the ±15 aa site windows: residues, X . * ? for any residue, [ST] for a choice.
  const AA20 = 'ACDEFGHIKLMNPQRSTVWY';
  const motifCache = new Map();
  function parseMotif(raw) {
    const s = String(raw || '').toUpperCase().replace(/\s+/g, '');
    if (!s) return null;
    if (motifCache.has(s)) return motifCache.get(s);
    let out;
    const toks = [];
    for (let k = 0; k < s.length && !out; k++) {
      const c = s[k];
      if (c === '[') {
        const e = s.indexOf(']', k);
        const set = e < 0 ? '' : s.slice(k + 1, e);
        if (e < 0) out = { error: 'Close the [ … ] group.' };
        else if (!set || /[^ACDEFGHIKLMNPQRSTVWY]/.test(set)) out = { error: 'Inside [ ] list amino-acid letters only, e.g. [ST].' };
        else { toks.push(set); k = e; }
      } else if ('X.*?'.includes(c)) toks.push(null);
      else if (AA20.includes(c)) toks.push(c);
      else out = { error: '"' + c + '" is not an amino-acid letter or wildcard.' };
    }
    if (!out) out = toks.length < 2 || toks.length > 10 ? { error: 'Use 2 to 10 residues (each [ ] group counts as one).' } : { toks, text: s };
    motifCache.set(s, out);
    return out;
  }
  // "site": some position of the motif that allows R sits on the methylarginine; "window": anywhere in ±15 aa.
  function motifMatch(win, m, where) {
    if (!win) return false;
    const ok = (start) => start >= 0 && m.toks.every((tk, i) => { const ch = win[start + i]; return ch && ch !== '_' && (tk === null || tk.includes(ch)); });
    if (where === 'window') { for (let st = 0; st + m.toks.length <= win.length; st++) if (ok(st)) return true; return false; }
    return m.toks.some((tk, i) => (tk === null || tk.includes('R')) && ok(15 - i));
  }
  function motifHint() {
    const n = S.filter((s) => !s.win).length;
    return 'X, . or * = any residue; [ST] = S or T; RGG,GAR = either. Adds <code>motif:</code> or <code>window:</code> to the search.' + (n ? ' ' + fmt(n) + ' sites without a sequence window are left out.' : '');
  }
  // A motif box (rail or overview) writes one motif: or window: term, replacing the pattern term it shows.
  function bindMotifBox(input, where, btn, msg, view) {
    const run = () => {
      const raw = input.value.toUpperCase().replace(/\s+/g, '');
      const key = where ? where.value : 'motif';
      const bad = raw.split(',').filter(Boolean).map((v) => resolveValue(key, v)).find((r) => r.error);
      input.setAttribute('aria-invalid', String(!!bad));
      if (msg) msg.innerHTML = bad ? '<span class="warn-text">' + esc(bad.error) + '</span>' : motifHint();
      if (bad || (!raw && !patternToken())) return;
      replaceToken(patternToken(), key, raw, view);
    };
    input.addEventListener('keydown', (e) => { if (e.key === 'Enter') { e.preventDefault(); run(); } });
    if (btn) btn.addEventListener('click', run);
    if (where) where.addEventListener('change', () => { if (input.value.trim()) run(); });
  }

  function windowHtml(win, flank) {
    if (!win) return '<span class="muted">—</span>';
    const c = 15;
    let out = '';
    for (let k = c - flank; k <= c + flank; k++) {
      const ch = win[k];
      if (k === c) out += '<span class="c">' + ch + '</span>';
      else if (ch === '_') out += '<span class="pad">·</span>';
      else if (ch === 'R') out += '<span class="r">R</span>';
      else if (ch === 'G') out += '<span class="g">G</span>';
      else out += ch;
    }
    return '<span class="win">' + out + '</span>';
  }

  const pill = (st) => '<span class="pill"><i class="st-' + st + '"></i>' + STATE_LABEL[st] + '</span>';
  const motifTags = (mask) => bits(mask).filter((b) => !M.motifs[b].startsWith('CARM1')).map((b) => '<span class="tag">' + MOTIF_LABEL[M.motifs[b]] + '</span>').join(' ') || '<span class="muted">—</span>';
  const ipTags = (mask) => bits(mask).map((b) => '<span class="tag ip" title="Recovered by the Maron 2021 ' + esc(IP_HELP[M.maron_ips[b]] || M.maron_ips[b]) + ' IP">' + ipLabel(M.maron_ips[b]) + '</span>').join(' ');
  function studyTags(mask, max) {
    const list = M.studyOrder.filter((k) => mask & (1 << k));
    const shown = list.slice(0, max || 3);
    const rest = list.slice(shown.length);
    return shown.map((k) => '<span class="tag src" title="' + esc(M.studies[k][1] || STUDY_HELP[M.studies[k][0]] || '') + '">' + esc(M.studies[k][0]) + '</span>').join(' ') +
      (rest.length ? ' <span class="tag" title="' + esc(rest.map((k) => M.studies[k][0]).join(', ')) + '">+' + rest.length + '</span>' : '');
  }

  function treatHtml(s) {
    if (!s.tr.length) return '<span class="muted">—</span>';
    const parts = [];
    for (let t = 0; t < 2; t++) {
      if (s.trf & (1 << (t * 2))) parts.push('<span class="trend down" title="Decreased with ' + M.maron_treatments[t] + ' (' + TREAT_HELP[M.maron_treatments[t]] + '), adj. p ≤ 0.05"><b>▼</b>' + M.maron_treatments[t] + '</span>');
      if (s.trf & (1 << (t * 2 + 1))) parts.push('<span class="trend up" title="Increased with ' + M.maron_treatments[t] + ', adj. p ≤ 0.05"><b>▲</b>' + M.maron_treatments[t] + '</span>');
    }
    return parts.join(' ') || '<span class="muted" title="Measured; no significant change">n.s.</span>';
  }

  function contextHtml(s) {
    const d = M.disorder_classes[s.dc];
    const c = M.domain_contexts[s.dctx];
    return '<span title="' + esc(DISORDER_HELP[d] || '') + '">' + esc(DISORDER_LABEL[d] || d) + '</span><br><span class="muted" title="' + esc(DOMAIN_HELP[c] || '') + '">' + esc(DOMAIN_LABEL[c] || c) + '</span>';
  }

  function miniMap(p, shown) {
    const W = 150, H = 20, len = p.len || p.maxPos * 1.05 || 1;
    let marks = '';
    for (const j of p.sites) {
      const s = S[j];
      const x = 1 + (s.pos / len) * (W - 2);
      const on = shown.has(j);
      marks += '<line x1="' + x.toFixed(1) + '" x2="' + x.toFixed(1) + '" y1="' + (on ? 2 : 7) + '" y2="' + (on ? 18 : 13) + '" stroke="' + (on ? 'var(--st-' + ['none', '1', '2', 'mix'][s.st] + ')' : 'var(--line-strong)') + '" stroke-width="' + (on ? 1.6 : 1) + '"/>';
    }
    return '<svg class="minimap" width="' + W + '" height="' + H + '" viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="Site positions along ' + esc(p.label) + '"><line x1="1" x2="' + (W - 1) + '" y1="10" y2="10" stroke="var(--line-strong)" stroke-width="2" stroke-linecap="round"/>' + marks + '</svg>';
  }

  // ----- result tables. Core columns always show; optional ones are switched in the Columns menu, remembered per table.
  // A column: {k, label, sort?, cls?, opt?, td(row)}. cls also hides a column on narrow screens (hide-sm, hide-xs).
  const COL_DEFAULT = { proteins: ['ix', 'idr', 'cond', 'density'], sites: ['ix', 'prmti', 'src'], psites: ['ix', 'prmti', 'src'] };
  function colsOn(table) {
    UI.cols = UI.cols || {};
    if (!UI.cols[table]) {
      let saved = null;
      try { saved = JSON.parse(store.get('mre-cols-' + table) || 'null'); } catch (e) { saved = null; }
      UI.cols[table] = new Set(Array.isArray(saved) ? saved : COL_DEFAULT[table]);
    }
    return UI.cols[table];
  }
  const colsShown = (table, cols) => { const on = colsOn(table); return cols.filter((c) => !c.opt || on.has(c.k)); };
  function colMenu(table, cols) {
    const on = colsOn(table);
    return '<details class="col-menu"><summary class="btn quiet small">Columns</summary><div class="col-pop" role="group" aria-label="Optional columns">' +
      cols.filter((c) => c.opt).map((c) => '<label class="check"><input type="checkbox" data-col="' + c.k + '"' + (on.has(c.k) ? ' checked' : '') + '> ' + esc(c.menu || c.label) + '</label>').join('') +
      '<button class="linkish" type="button" data-col-reset>Default columns</button></div></details>';
  }
  function bindColMenu(el, table, redraw) {
    const box = $('.col-menu', el);
    if (!box) return;
    const save = () => store.set('mre-cols-' + table, JSON.stringify(Array.from(colsOn(table))));
    $$('[data-col]', box).forEach((cb) => cb.addEventListener('change', () => {
      const on = colsOn(table);
      if (cb.checked) on.add(cb.dataset.col); else on.delete(cb.dataset.col);
      save(); redraw(true);
    }));
    $('[data-col-reset]', box).addEventListener('click', () => { UI.cols[table] = new Set(COL_DEFAULT[table]); store.set('mre-cols-' + table, null); redraw(true); });
  }
  function colHead(cols, sortState, key) {
    return '<thead><tr>' + cols.map((c) => '<th' + (c.cls ? ' class="' + c.cls + '"' : '') + (c.tip ? ' title="' + esc(c.tip) + '"' : '') + (c.sort ? ' data-sort="' + c.sort + '"' + (key === c.sort ? ' aria-sort="' + (sortState.dir > 0 ? 'ascending' : 'descending') + '"' : '') : '') + '>' + c.label + '</th>').join('') + '</tr></thead>';
  }
  const colLabel = (c) => c.label.replace(/<[^>]+>/g, '').replace(/^Actions$/, '');
  const colRow = (cols, row) => cols.map((c) => '<td' + (c.tdCls || c.cls ? ' class="' + ((c.tdCls || '') + ' ' + (c.cls || '')).trim() + '"' : '') + ' data-label="' + esc(colLabel(c)) + '"><div class="cv">' + c.td(row) + '</div></td>').join('');

  function resultsHead(title, sub, actions) {
    return '<div class="results-head"><div><h2>' + title + '</h2><p>' + sub + '</p></div><div class="actions">' + (actions || '') + '</div></div>';
  }

  // Summary card for a protein-level selection (condensate, class, disease or study) above the results.
  function contextCards() {
    const ov = M.overview;
    const cards = [];
    const card = (eyebrow, title, body, links) => cards.push('<div class="ctx-card"><div><span class="eyebrow">' + esc(eyebrow) + '</span><h3>' + esc(title) + '</h3><p>' + body + '</p></div>' + (links ? '<div class="actions">' + links + '</div>' : '') + '</div>');
    const one = (key) => { const v = singleTerm(key); return v ? v.k : null; };
    if (one('condensate') != null) {
      const c = M.condensates[one('condensate')];
      card('CD-CODE condensate · ' + c[2], c[0], fmt(c[4]) + ' of ' + fmt(c[3]) + ' members (' + pct(c[4], c[3]) + '%) carry a methylarginine. ' +
        (c[5] != null ? 'Methylated fraction vs the global rate: log2 OR <b>' + fmtL2(c[5]) + '</b>, q ' + fmtP(c[6]) + '.' : 'Too few members for the condensate enrichment test.'),
        '<a class="btn quiet small" href="https://cd-code.org/condensates/' + encodeURIComponent(c[1]) + '" target="_blank" rel="noopener">CD-CODE ' + esc(c[1]) + ' ↗</a>');
    }
    if (one('class') != null) {
      const name = M.func_classes[one('class')];
      const r = ov.functional_classes.find((x) => x[0] === name);
      if (r) card('Functional class', name, fmt(r[1]) + ' methylarginine proteins (' + (100 * r[2]).toFixed(1) + '% of the set vs ' + (100 * r[3]).toFixed(1) + '% of the reviewed background); log2 OR <b>' + fmtL2(r[4]) + '</b>, q ' + fmtP(r[5]) + '.');
    }
    if (one('subclass') != null) {
      const name = M.func_subclasses[one('subclass')];
      const r = (ov.functional_subclasses || []).find((x) => x[0] === name);
      if (r) card('Functional subclass of ' + (r[6] || '—'), name, fmt(r[1]) + ' methylarginine proteins; log2 OR <b>' + fmtL2(r[4]) + '</b>, q ' + fmtP(r[5]) + ' vs the reviewed background.');
    }
    if (one('disease') != null) {
      if (one('disease') === 'als') card('Audited disease family', 'ALS / motor neuron / FTD (full gene set)', 'Every gene in the audited ALS/MN/FTD family, not only genes in enriched terms.');
      else {
        const name = M.disease_families[one('disease')];
        const r = ov.disease_families.find((x) => x[0] === name);
        if (r) card('Audited disease family (DisGeNET)', name, fmt(r[1]) + ' methylarginine genes among ' + fmt(r[2]) + ' family genes; OR <b>' + (r[3] || 0).toFixed(2) + '</b>, BH-adjusted p ' + fmtP(r[4]) + '. Proteins listed here belong to an enriched term of this family.');
      }
    }
    if (one('study') != null) {
      const k = one('study');
      const st = M.studies[k];
      card('Provenance', st[0], esc(st[1] && st[1] !== 'This study.' ? st[1] : STUDY_HELP[st[0]] || '') + ' · ' + fmt(st[2]) + ' sites, ' + fmt(st[3]) + ' reported only here.');
    }
    return cards.length ? '<div class="ctx-cards">' + cards.join('') + '</div>' : '';
  }

  // When nothing matches, show what each constraint would give back on its own.
  function unrevHelp() {
    const prots = P.filter((p) => p.unrev);
    const sites = prots.reduce((n, p) => n + p.sites.length, 0);
    const dups = S.filter((s) => s.dup != null).length;
    return 'Drops ' + fmt(sites) + ' sites on ' + fmt(prots.length) + ' unreviewed (TrEMBL) entries. ' + fmt(dups) + ' of them repeat a site at the same position on the reviewed entry of the same gene, so the union counts them twice.';
  }

  function emptyState() {
    const cons = constraints();
    if (!cons.length) return '<div class="card empty"><b>No sites match</b>Broaden the search.</div>';
    const saveF = cloneFilters(F), saveQ = Q;
    const tries = cons.map((c, k) => {
      // Dry run: clearing the search or a term also rewrites the box and the URL, so swap Q directly.
      if (c.dryQ) Q = c.dryQ(); else c.clear();
      const n = filterSites(F, Q).sites.length;
      F = cloneFilters(saveF); Q = saveQ;
      return [k, c.label, n];
    }).filter((t) => t[2] > 0).sort((a, b) => b[2] - a[2]).slice(0, 5);
    return '<div class="card empty"><b>' + (cons.length === 1 ? 'No sites match' : 'No sites match all ' + cons.length + ' constraints') + '</b>' +
      (tries.length ? 'Removing one of these brings sites back:<div class="empty-tries">' + tries.map((t) => '<button class="chip" type="button" data-k="' + t[0] + '">Remove “' + esc(t[1]) + '” · ' + fmt(t[2]) + ' sites</button>').join('') + '</div>' : 'Remove a filter or broaden the search.') +
      '<p><button class="btn small" type="button" id="empty-reset">Clear all</button></p></div>';
  }

  function bindEmpty(root) {
    const b = $('#empty-reset', root);
    if (b) b.addEventListener('click', clearAll);
    const cons = constraints();
    $$('.empty-tries .chip', root).forEach((c) => c.addEventListener('click', () => { cons[+c.dataset.k].clear(); afterFilterChange(); }));
  }

  // -------------------------------------------------------------- proteins
  // Interaction evidence of a protein's sites in the current selection: experimental, predicted only.
  function protIx(p) {
    let e = 0, pr = 0;
    for (const j of p.sites) { if (!R.shown.has(j)) continue; const st = ixStatus(S[j]); if (st === 'exp') e++; else if (st === 'pred') pr++; }
    return { e, pr };
  }
  const PROTEIN_COLS = [
    { k: 'protein', label: 'Protein', sort: 'gene', td: (p) => '<span class="gene">' + esc(p.label) + '</span> <span class="acc-in">' + esc(p.acc) + '</span><div class="pname" title="' + esc(p.name) + '">' + esc(p.name || '—') + '</div>' },
    { k: 'sites', label: 'Sites', sort: 'sites', cls: 'r', tdCls: 'num', td: (p) => { const n = R.counts[p.i]; return n === p.sites.length ? fmt(n) : fmt(n) + '<span class="muted"> / ' + fmt(p.sites.length) + '</span>'; } },
    { k: 'map', label: 'Site map', cls: 'hide-xs', td: (p) => miniMap(p, R.shown) },
    { k: 'ix', label: 'Interaction evidence', menu: 'Interaction evidence (sites)', sort: 'ix', opt: true, tip: 'Sites in this selection with experimental contact evidence, and with predicted evidence only',
      td: (p) => { if (!ixReady()) return '<span class="muted">…</span>'; const c = protIx(p); return c.e || c.pr ? (c.e ? '<span class="ixs ixs-exp" title="Sites with experimental contact evidence"><i></i>' + c.e + ' exp.</span> ' : '') + (c.pr ? '<span class="ixs ixs-pred" title="Sites with predicted evidence only"><i></i>' + c.pr + ' pred.</span>' : '') : '<span class="muted">—</span>'; } },
    { k: 'idr', label: 'In IDR', sort: 'idr', opt: true, cls: 'r hide-xs', tdCls: 'num', td: (p) => pct(p.nIdr, p.sites.length) + '%' },
    { k: 'cond', label: 'Condensates', sort: 'cond', opt: true, cls: 'r hide-sm', tdCls: 'num', td: (p) => (p.cond.length || '<span class="muted">0</span>') },
    { k: 'density', label: 'Densest window', sort: 'density', opt: true, cls: 'hide-sm', tdCls: 'nw', td: (p) => (p.best && p.best[5] != null && p.best[5] <= PADJ ? p.best[3] + ' in ' + p.best[0] + ' aa' : '<span class="muted">—</span>') },
    { k: 'arg', label: 'Arg enrichment', sort: 'arg', opt: true, cls: 'r hide-sm', tdCls: 'num', tip: 'Shrunken log2 odds that this protein\'s arginines are methylated more than the proteome rate', td: (p) => (p.shrink == null ? '<span class="muted">—</span>' : fmtL2(p.shrink)) },
    { k: 'class', label: 'Functional class', opt: true, cls: 'hide-sm', td: (p) => '<div class="tags">' + (bits(p.func).map((b) => '<span class="tag">' + esc(M.func_classes[b].split(' / ')[0]) + '</span>').join('') || '<span class="muted">—</span>') + '</div>' },
  ];
  const P_SORT = {
    gene: (p) => p.label.toLowerCase(),
    ix: (p) => { if (!ixReady()) return 0; const c = protIx(p); return c.e * 1000 + c.pr; },
    sites: (p) => R.counts[p.i],
    total: (p) => p.sites.length,
    idr: (p) => (p.sites.length ? p.nIdr / p.sites.length : 0),
    cond: (p) => p.cond.length,
    density: (p) => (p.best && p.best[5] != null && p.best[5] <= PADJ ? p.best[4] : -99),
    arg: (p) => (p.shrink == null ? -99 : p.shrink),
  };

  function renderProteins(el) {
    const key = P_SORT[UI.sortP.key] ? UI.sortP.key : 'sites';
    const fn = P_SORT[key];
    const dir = UI.sortP.dir;
    const rows = R.prots.slice().sort((a, b) => {
      const x = fn(a), y = fn(b);
      if (x < y) return -dir;
      if (x > y) return dir;
      return b.sites.length - a.sites.length;
    });
    const shownRows = rows.slice(0, UI.pageP);
    const cols = colsShown('proteins', PROTEIN_COLS);
    el.innerHTML = resultsHead('Proteins', fmt(rows.length) + ' of ' + fmt(P.length) + ' proteins · ' + fmt(R.sites.length) + ' sites',
      '<button class="btn quiet small" type="button" id="copy-genes">Copy gene list</button><button class="btn quiet small" type="button" data-go-figs>Make figures</button><button class="btn small" type="button" id="dl-proteins">Download TSV</button>') +
      contextCards() + profileCardHtml() +
      (rows.length ? '<div class="card tbl-card"><div class="tbl-tools"><span class="muted small">Click a protein to open it.</span>' + colMenu('proteins', PROTEIN_COLS) + '</div><div class="tbl-wrap"><table class="tbl stackable">' + colHead(cols, UI.sortP, key) + '<tbody>' +
        shownRows.map((p) => '<tr class="row" tabindex="0" data-i="' + p.i + '">' + colRow(cols, p) + '</tr>').join('') +
        '</tbody></table></div>' +
        '<div class="tbl-foot"><span>Showing ' + fmt(shownRows.length) + ' of ' + fmt(rows.length) + '</span>' +
        (rows.length > shownRows.length ? '<button class="btn quiet small" type="button" id="more">Show 50 more</button>' : '') + '</div></div>' +
        citeBoxHtml(R.sites, R.sites.length === S.length ? 'the sites in the union' : 'the ' + fmt(R.sites.length) + ' sites on these proteins') : emptyState());
    bindEmpty(el);
    bindRows(el, (tr) => { if (tr.dataset.i != null) openProtein(+tr.dataset.i); });
    bindSort(el, UI.sortP, () => renderProteins(el));
    bindColMenu(el, 'proteins', () => { renderProteins(el); const m = $('.col-menu', el); if (m) m.open = true; });
    bindCiteBox(el);
    const more = $('#more', el);
    if (more) more.addEventListener('click', () => { UI.pageP += 50; renderProteins(el); });
    $('#dl-proteins', el).addEventListener('click', () => download('methylarginine_proteins.tsv', proteinTsv(rows)));
    bindProfile(el);
    $$('[data-go-figs]', el).forEach((b) => b.addEventListener('click', () => { go('figures'); window.scrollTo({ top: 0 }); }));
    $('#copy-genes', el).addEventListener('click', (e) => copyText(rows.map((p) => p.gene || p.acc).join('\n'), e.currentTarget));
  }

  function bindRows(el, fn) {
    $$('tr.row', el).forEach((tr) => {
      tr.addEventListener('click', (e) => { if (!e.target.closest('a,button')) fn(tr); });
      tr.addEventListener('keydown', (e) => { if (e.key === 'Enter') fn(tr); });
    });
  }

  function bindSort(el, sort, redraw) {
    $$('th[data-sort]', el).forEach((th) => th.addEventListener('click', () => {
      const k = th.dataset.sort;
      if (sort.key === k) sort.dir = -sort.dir; else { sort.key = k; sort.dir = ['gene', 'pos', 'win', 'mot', 'ctx'].includes(k) ? 1 : -1; }
      redraw();
    }));
  }

  // ----------------------------------------------------------------- sites
  // Quick filter above the sites table: the selection split by interaction status, counted without any interaction: term.
  function ixSegHtml() {
    if (!ixReady()) return '';
    const cur = Q.tokens.filter((t) => t.key === 'interaction');
    const base = cur.length ? filterSites(F, parseQuery(queryText(Q.plain, Q.tokens.filter((t) => t.key !== 'interaction')))).sites : R.sites;
    const n = { exp: 0, pred: 0, none: 0 };
    for (const s of base) n[ixStatus(s)]++;
    const active = !cur.length ? '' : cur.length === 1 && cur[0].vals.length === 1 ? cur[0].vals[0].text : 'other';
    const b = (v, label, cnt, dot) => '<button type="button" class="seg-btn" data-ixseg="' + v + '" aria-pressed="' + (active === v) + '"' + (cnt ? '' : ' disabled') + '>' + (dot ? '<i class="ixd ixd-' + dot + '"></i>' : '') + label + ' <b>' + fmt(cnt) + '</b></button>';
    return '<div class="ix-seg card-flat" role="group" aria-label="Interaction evidence at the site"><span class="flabel">Interaction evidence at the site' + docLink('ix-status') + '</span><div class="ixseg-btns">' +
      b('', 'All', base.length) + b('experimental', 'Experimental', n.exp, 'exp') + b('predicted-only', 'Predicted only', n.pred, 'pred') + b('none', 'None recorded', n.none, 'none') + '</div></div>';
  }
  function bindIxSeg(el) {
    $$('[data-ixseg]', el).forEach((b) => b.addEventListener('click', () => setTerm('interaction', b.dataset.ixseg)));
  }
  // Window and motif share a cell: the motif tags sit under the sequence they describe.
  const winCell = (s) => windowHtml(s.win, 7) + (s.mis ? ' <span class="warn" title="Canonical residue at this position is not R">!</span>' : '') +
    (bits(s.mot).some((b) => !M.motifs[b].startsWith('CARM1')) ? '<div class="tags win-mot">' + motifTags(s.mot) + '</div>' : '');
  const SITE_COLS = [
    { k: 'protein', label: 'Protein', sort: 'gene', td: (s) => '<span class="gene">' + esc(P[s.p].label) + '</span><span class="acc">' + esc(P[s.p].acc) + '</span>' },
    { k: 'site', label: 'Site', sort: 'pos', tdCls: 'num', td: (s) => '<b>' + rl(s.pos, P[s.p]) + '</b>' + (P[s.p].hist ? '<span class="acc">UniProt R' + s.pos + '</span>' : '') },
    { k: 'state', label: 'State', sort: 'state', td: (s) => pill(s.st) },
    { k: 'win', label: 'Window ±7 · motifs', sort: 'win', td: winCell },
    { k: 'ix', label: 'Interaction evidence', sort: 'ix', opt: true, tdCls: 'ix-col', tip: 'Experimental: RNA/DNA crosslink, RBDmap or a PDB interface at this residue. Predicted: at the interface of a predicted complex.', td: (s) => ixCell(s, 2) },
    { k: 'prmti', label: 'PRMTi', sort: 'prmti', opt: true, menu: 'PRMT inhibitor response', tip: 'Maron 2021 Table S3: change after GSK591 (PRMT5i) or MS023 (type I PRMTi), adj. p ≤ 0.05', td: treatHtml },
    { k: 'src', label: 'Reported by', sort: 'src', opt: true, cls: 'hide-sm', td: (s) => '<div class="tags">' + studyTags(s.study, 2) + '</div>' },
    { k: 'ctx', label: 'Context', sort: 'ctx', opt: true, menu: 'Disorder and domain context', cls: 'hide-xs', tdCls: 'nw', td: contextHtml },
    { k: 'nn', label: 'Next methyl-R', sort: 'nn', opt: true, cls: 'r hide-sm', tdCls: 'num', td: (s) => (s.nn == null ? '<span class="muted">—</span>' : s.nn + ' aa') },
    { k: 'ip', label: 'Maron IP', sort: 'ip', opt: true, menu: 'Maron methyl-IP', cls: 'hide-sm', td: (s) => '<div class="tags">' + (ipTags(s.ip) || '<span class="muted">—</span>') + '</div>' },
  ];
  const S_SORT = {
    ix: (s) => ixRank(s),
    gene: (s) => P[s.p].label.toLowerCase(),
    pos: (s) => s.pos,
    state: (s) => s.st,
    nn: (s) => s.nn,
    src: (s) => popcount(s.study),
    win: (s) => (s.win ? s.win.slice(8, 23) : null),
    mot: (s) => bits(s.mot).filter((b) => !M.motifs[b].startsWith('CARM1')).map((b) => MOTIF_LABEL[M.motifs[b]]).join(' ') || null,
    ctx: (s) => (DISORDER_LABEL[M.disorder_classes[s.dc]] || M.disorder_classes[s.dc] || '') + ' ' + (DOMAIN_LABEL[M.domain_contexts[s.dctx]] || M.domain_contexts[s.dctx] || ''),
    ip: (s) => s.ip || null,
    iface: (s) => s.iface || null,
    // Most informative first when descending: decreased by both inhibitors, GSK591, MS023, increased, n.s., not measured.
    prmti: (s) => (!s.tr.length ? null : (s.trf & 1) && (s.trf & 4) ? 4 : s.trf & 1 ? 3 : s.trf & 4 ? 2 : s.trf & 10 ? 1 : 0),
  };

  function renderSites(el) {
    const key = S_SORT[UI.sortS.key] ? UI.sortS.key : 'gene';
    const fn = S_SORT[key];
    const dir = UI.sortS.dir;
    const rows = R.sites.slice().sort((a, b) => {
      const x = fn(a), y = fn(b);
      // Missing values sort last in either direction.
      if (x == null || y == null) { if (x != null) return -1; if (y != null) return 1; }
      else if (x < y) return -dir;
      else if (x > y) return dir;
      return a.p - b.p || a.pos - b.pos;
    });
    const shownRows = rows.slice(0, UI.pageS);
    const cols = colsShown('sites', SITE_COLS);
    el.innerHTML = resultsHead('Sites', fmt(rows.length) + ' of ' + fmt(S.length) + ' methylarginine sites on ' + fmt(R.prots.length) + ' proteins',
      '<button class="btn quiet small" type="button" data-go-figs>Make figures</button><button class="btn small" type="button" id="dl-sites">Download TSV</button>') +
      contextCards() + profileCardHtml() + ixSegHtml() +
      (rows.length ? '<div class="card tbl-card"><div class="tbl-tools"><span class="muted small">Click a site to open it on its protein.</span>' + colMenu('sites', SITE_COLS) + '</div><div class="tbl-wrap"><table class="tbl stackable">' + colHead(cols, UI.sortS, key) + '<tbody>' +
        shownRows.map((s) => '<tr class="row" tabindex="0" data-j="' + s.j + '">' + colRow(cols, s) + '</tr>').join('') +
        '</tbody></table></div><div class="tbl-foot"><span>Showing ' + fmt(shownRows.length) + ' of ' + fmt(rows.length) + '</span>' +
        (rows.length > shownRows.length ? '<button class="btn quiet small" type="button" id="more">Show 100 more</button>' : '') + '</div></div>' +
        citeBoxHtml(R.sites, R.sites.length === S.length ? 'the sites in the union' : 'these ' + fmt(R.sites.length) + ' sites') : emptyState());
    bindEmpty(el);
    bindRows(el, (tr) => { if (tr.dataset.j == null) return; const s = S[+tr.dataset.j]; openProtein(s.p, s.pos); });
    bindContactLinks(el);
    bindSort(el, UI.sortS, () => renderSites(el));
    bindColMenu(el, 'sites', () => { renderSites(el); const m = $('.col-menu', el); if (m) m.open = true; });
    bindCiteBox(el);
    bindIxSeg(el);
    const more = $('#more', el);
    if (more) more.addEventListener('click', () => { UI.pageS += 100; renderSites(el); });
    $('#dl-sites', el).addEventListener('click', () => download('methylarginine_sites.tsv', siteTsv(rows)));
    bindProfile(el);
    $$('[data-go-figs]', el).forEach((b) => b.addEventListener('click', () => { go('figures'); window.scrollTo({ top: 0 }); }));
  }

  // -------------------------------------------------------------- charts
  // A chart spec renders as interactive HTML bars and exports as a standalone SVG/PNG/TSV.
  function niceScale(rows) {
    const vals = rows.map((r) => r.value).filter((v) => v != null);
    let lo = Math.min(0, ...vals), hi = Math.max(0, ...vals);
    const raw = (hi - lo || 1) / 5;
    const mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const step = [1, 2, 5, 10].map((m) => m * mag).find((st) => st >= raw);
    lo = Math.floor(lo / step) * step; hi = Math.ceil(hi / step) * step;
    const ticks = [];
    for (let t = lo; t <= hi + step / 2; t += step) ticks.push(+t.toFixed(6));
    return { lo, hi, ticks, xp: (v) => ((v - lo) / (hi - lo)) * 100 };
  }

  function hbars(spec) {
    const rows = spec.rows;
    const sc = niceScale(rows);
    const vfmt = spec.format || fmtL2;
    const tfmt = spec.tick || ((t) => t);
    return '<div class="hbars' + (spec.extraCol ? ' has-extra' : '') + '" role="list">' + rows.map((r, k) => {
      const v = r.value == null ? 0 : r.value;
      const a = sc.xp(Math.min(0, v)), b = sc.xp(Math.max(0, v));
      return '<button class="hb' + (r.sub ? ' sub' : '') + '" type="button" role="listitem" data-k="' + k + '" aria-label="' + esc(r.label + ': ' + (spec.valueLabel || 'log2 OR') + ' ' + vfmt(r.value) + (spec.noSig ? '' : r.sig ? ', significant' : ', not significant')) + '">' +
        '<span class="hb-label">' + esc(r.label) + '</span>' + (spec.extraCol ? '<span class="hb-extra">' + esc(r.extra || '') + '</span>' : '') +
        '<span class="hb-track"><span class="hb-zero" style="left:' + sc.xp(0).toFixed(2) + '%"></span>' +
        '<span class="hb-bar' + (r.sig || spec.noSig ? '' : ' ns') + (v < 0 ? ' neg' : '') + '" style="left:' + a.toFixed(2) + '%;width:' + Math.max(0.6, b - a).toFixed(2) + '%"></span></span>' +
        '<span class="hb-val">' + vfmt(r.value) + '</span></button>';
    }).join('') + '</div>' +
      '<div class="hb-axis' + (spec.extraCol ? ' has-extra' : '') + '" aria-hidden="true"><span></span>' + (spec.extraCol ? '<span></span>' : '') + '<span class="ticks">' + sc.ticks.map((t) => '<span style="left:' + sc.xp(t).toFixed(2) + '%">' + tfmt(t) + '</span>').join('') + '</span><span></span></div>' +
      '<div class="chart-key">' + (spec.noSig ? '' : '<span><i></i>q ≤ 0.05</span><span><i class="ns"></i>not significant</span>') + '<span>' + esc(spec.axis || 'log2 odds ratio vs background') + '</span></div>';
  }

  const expButtons = (id) => '<div class="exp" role="group" aria-label="Export chart"><button type="button" data-exp="svg" data-chart="' + id + '">SVG</button><button type="button" data-exp="png" data-chart="' + id + '">PNG</button><button type="button" data-exp="tsv" data-chart="' + id + '">TSV</button></div>';

  function chartCard(spec, cls) {
    UI.charts[spec.id] = spec;
    return '<div class="card chart-card' + (cls ? ' ' + cls : '') + '" id="' + spec.id + '"><div class="card-head"><div><h3>' + esc(spec.title) + '</h3><p class="sub">' + spec.subHtml + '</p></div>' + expButtons(spec.id) + '</div>' +
      hbars(spec) + (spec.footHtml || '') + '</div>';
  }

  function chartSvg(spec) {
    const pal = EXPORT_PAL;
    const rows = spec.rows;
    const W = 760, labelW = 270, extraW = spec.extraCol ? 64 : 0, valW = 64;
    const x0 = 24 + labelW + extraW, x1 = W - valW - 16;
    const rowH = 22, top = 74;
    const sc = niceScale(rows);
    const X = (v) => x0 + (sc.xp(v) / 100) * (x1 - x0);
    const vfmt = spec.format || fmtL2;
    const tfmt = spec.tick || ((t) => t);
    const H = top + rows.length * rowH + 46;
    let g = '<rect width="' + W + '" height="' + H + '" fill="' + pal.bg + '"/>' +
      '<text x="24" y="32" font-size="17" font-weight="700" fill="' + pal.ink + '">' + esc(spec.title) + '</text>' +
      '<text x="24" y="52" font-size="11" fill="' + pal.muted + '">' + esc(spec.subText || '') + '</text>';
    for (const t of sc.ticks) g += '<line x1="' + r1(X(t)) + '" x2="' + r1(X(t)) + '" y1="' + (top - 6) + '" y2="' + (top + rows.length * rowH) + '" stroke="' + pal.grid + '"/>';
    g += '<line x1="' + r1(X(0)) + '" x2="' + r1(X(0)) + '" y1="' + (top - 8) + '" y2="' + (top + rows.length * rowH + 2) + '" stroke="' + pal.backbone + '"/>';
    rows.forEach((r, k) => {
      const y = top + k * rowH;
      const v = r.value == null ? 0 : r.value;
      const a = X(Math.min(0, v)), b = X(Math.max(0, v));
      const fill = v < 0 ? pal.barNeg : pal.bar;
      const ns = !spec.noSig && !r.sig;
      g += '<text x="' + (r.sub ? 40 : 24) + '" y="' + (y + 14) + '" font-size="' + (r.sub ? 10.5 : 11.5) + '" fill="' + (r.sub ? pal.muted : pal.ink) + '">' + esc(r.label.length > 46 ? r.label.slice(0, 44) + '…' : r.label) + '</text>';
      if (spec.extraCol) g += '<text x="' + (24 + labelW + extraW - 6) + '" y="' + (y + 14) + '" font-size="10.5" fill="' + pal.muted + '" text-anchor="end">' + esc(r.extra || '') + '</text>';
      g += '<rect x="' + r1(a) + '" y="' + (y + 6) + '" width="' + r1(Math.max(1, b - a)) + '" height="10" rx="2" fill="' + (ns ? pal.bg : fill) + '" stroke="' + fill + '" stroke-width="' + (ns ? 1.4 : 0) + '"/>';
      g += '<text x="' + (W - 16) + '" y="' + (y + 14) + '" font-size="10.5" fill="' + pal.muted + '" text-anchor="end">' + esc(vfmt(r.value)) + '</text>';
    });
    const ay = top + rows.length * rowH + 16;
    for (const t of sc.ticks) g += '<text x="' + r1(X(t)) + '" y="' + ay + '" font-size="10" fill="' + pal.muted + '" text-anchor="middle">' + esc(String(tfmt(t))) + '</text>';
    g += '<text x="' + r1((x0 + x1) / 2) + '" y="' + (ay + 16) + '" font-size="10.5" fill="' + pal.muted + '" text-anchor="middle">' + esc(spec.axis || 'log2 odds ratio vs background') + '</text>';
    if (!spec.noSig) {
      g += '<rect x="24" y="' + (ay + 6) + '" width="14" height="9" rx="2" fill="' + pal.bar + '"/><text x="42" y="' + (ay + 14) + '" font-size="10" fill="' + pal.muted + '">q ≤ 0.05</text>' +
        '<rect x="96" y="' + (ay + 6) + '" width="14" height="9" rx="2" fill="' + pal.bg + '" stroke="' + pal.bar + '" stroke-width="1.4"/><text x="114" y="' + (ay + 14) + '" font-size="10" fill="' + pal.muted + '">not significant</text>';
    }
    return { svg: '<svg xmlns="http://www.w3.org/2000/svg" width="' + W + '" height="' + H + '" viewBox="0 0 ' + W + ' ' + H + '" font-family=\'' + pal.font + '\'>' + g + '</svg>', W, H };
  }

  function chartTsv(spec) {
    return ['label\t' + (spec.valueLabel || 'log2_odds_ratio') + '\tq_value\t' + (spec.extraName || 'detail'),
      ...spec.rows.map((r) => [r.label, r.value, r.q == null ? '' : r.q, r.extra || ''].map(cell).join('\t'))].join('\n') + '\n';
  }

  function bindChartExports(root) {
    $$('[data-exp][data-chart]', root).forEach((b) => b.addEventListener('click', () => {
      const spec = UI.charts[b.dataset.chart];
      const name = 'methylarginine_' + spec.id.replace(/^ch-/, '');
      if (b.dataset.exp === 'tsv') { download(name + '.tsv', chartTsv(spec)); return; }
      exportFigure(b.dataset.exp, chartSvg(spec), name);
    }));
  }

  function bindBars(container, rows, onClick) {
    $$('.hb[data-k]', container).forEach((b) => {
      const r = rows[+b.dataset.k];
      b.addEventListener('click', () => onClick(r));
      b.addEventListener('mouseenter', (e) => showTip(e, r.title));
      b.addEventListener('mousemove', moveTip);
      b.addEventListener('mouseleave', hideTip);
    });
  }

  // -------------------------------------------------------------- overview
  function renderOverview(el) {
    const ov = M.overview;
    UI.charts = {};
    const topSites = P.slice().sort((a, b) => b.sites.length - a.sites.length).slice(0, 14);
    const topArg = P.filter((p) => p.argQ != null && p.argQ <= PADJ && p.shrink != null).sort((a, b) => b.shrink - a.shrink).slice(0, 14);
    const topDense = P.filter((p) => p.sig.length).sort((a, b) => Math.max(...b.sig.map((w) => w[4])) - Math.max(...a.sig.map((w) => w[4]))).slice(0, 14);
    const examples = [['FUS', 'One protein'], ['FUS R216', 'One site'], ['TAF15 EWSR1 FUS', 'A gene list (FET family)'], ['SNRPB SNRPD1 SNRPD3', 'Sm proteins'], ['TDP43', 'A synonym'], ['splicing factor', 'Protein name']];
    const geneChips = (list, sub) => '<div class="gene-chips">' + list.map((p) => '<button class="gene-chip" type="button" data-i="' + p.i + '">' + esc(p.label) + '<span>' + sub(p) + '</span></button>').join('') + '</div>';
    const sigq = (q) => q != null && q <= PADJ;

    const study = {
      id: 'ch-study', title: 'Provenance', subHtml: 'Sites reported by each study or database; hover for sites reported only there.' + docLink('ov-provenance'),
      subText: 'Methylarginine sites reported by each study or database', noSig: true, format: (v) => fmt(v), tick: compact, axis: 'methylarginine sites reported', valueLabel: 'sites',
      extraCol: true, extraName: 'reported_only_here',
      rows: M.studyOrder.map((k) => ({ label: M.studies[k][0], value: M.studies[k][2], sig: true, key: k, extra: fmt(M.studies[k][3]) + ' only',
        title: '<b>' + esc(M.studies[k][0]) + '</b><br>' + esc(M.studies[k][1] && M.studies[k][1] !== 'This study.' ? M.studies[k][1] : STUDY_HELP[M.studies[k][0]] || '') + '<br>' + fmt(M.studies[k][2]) + ' sites, ' + fmt(M.studies[k][3]) + ' reported only here' })),
    };
    const repDist = [0, 0, 0, 0, 0];
    for (const s of S) repDist[Math.min(5, popcount(s.study)) - 1]++;
    const motif = {
      id: 'ch-motif', title: 'Sequence motifs', subHtml: 'Methylated vs other arginines of the same proteins (±7 aa).' + docLink('ov-motifs'),
      subText: 'Methylated vs unmodified arginines in the same proteins (±7 window)',
      rows: ov.motifs.map((r) => ({ label: MOTIF_LABEL[r[0]] + ' · ' + MOTIF_HELP[r[0]], value: r[4], q: r[5], sig: sigq(r[5]), key: r[0], extra: fmt(r[1]),
        title: '<b>' + esc(MOTIF_LABEL[r[0]]) + '</b> (' + esc(MOTIF_HELP[r[0]]) + ')<br>' + fmt(r[1]) + ' methyl-R sites (' + (100 * r[2]).toFixed(1) + '%) vs ' + (100 * r[3]).toFixed(1) + '% of unmodified R in the same proteins<br>log2 OR ' + fmtL2(r[4]) + ', q ' + fmtP(r[5]) })),
      extraName: 'sites',
      footHtml: '<form class="ov-motif" id="ov-motif" autocomplete="off"><label class="flabel" for="ov-motifq">Find any motif</label><div class="motif-search">' +
        '<input id="ov-motifq" type="text" maxlength="60" spellcheck="false" placeholder="e.g. GXRG, [RK]XR, PRPR" aria-describedby="ov-motifmsg">' +
        '<select id="ov-motifwhere" aria-label="Where the motif may sit"><option value="motif">covering the methyl-R</option><option value="window">within ±15 aa</option></select>' +
        '<button class="btn small" type="submit">Search</button></div><p class="note" id="ov-motifmsg">Click a bar, or type a pattern: X = any residue, [ST] = S or T. Both add <code>motif:</code> to the search.</p></form>',
    };
    // Members counted: the pipeline's CD-CODE table (all memberships) or only CD-CODE confidence >= threshold.
    const thr = ixReady() && IX.condConf[String(UI.condThr)] ? UI.condThr : 0;
    const condAll = (thr ? IX.condConf[String(thr)].map(([k, n, m, l2, , q]) => ({ k, c: M.condensates[k], n, m, l2, q })) :
      M.condensates.map((c, k) => ({ k, c, n: c[3], m: c[4], l2: c[5], q: c[6] }))).filter((x) => x.l2 != null).sort((a, b) => b.l2 - a.l2);
    const condShown = UI.allCond ? condAll : condAll.slice(0, 14);
    const thrSel = '<label class="fig-ctl cond-thr">Count members <select id="cond-thr">' + [[0, 'all memberships'], [2, 'CD-CODE confidence ≥ 2'], [3, 'CD-CODE confidence ≥ 3'], [4, 'CD-CODE confidence ≥ 4']]
      .map(([v, l]) => '<option value="' + v + '"' + (v === thr ? ' selected' : '') + (v && !ixReady() ? ' disabled' : '') + '>' + l + '</option>').join('') + '</select></label>';
    const cond = {
      id: 'ch-cond', title: 'Biomolecular condensates',
      subHtml: 'CD-CODE condensates: methylated share of members vs the global rate.' + docLink('ov-condensates') + '<br>' + thrSel,
      subText: 'CD-CODE condensates: methylated members vs the global rate' + (thr ? ' (confidence ≥ ' + thr + ')' : ' (all memberships)'), extraCol: true, extraName: 'methylated_of_members',
      rows: condShown.map(({ c, k, n, m, l2, q }) => ({ label: c[0], value: l2, q, sig: sigq(q), key: k, extra: m + '/' + n,
        title: '<b>' + esc(c[0]) + '</b> · CD-CODE ' + esc(c[1]) + (c[2] !== 'biomolecular' ? ' · ' + esc(c[2]) : '') + '<br>' + m + ' of ' + n + ' counted members methylated' + (thr ? ' (confidence ≥ ' + thr + ')' : '') + '<br>log2 OR ' + fmtL2(l2) + ', q ' + fmtP(q) })),
      footHtml: '<button class="linkish" type="button" id="cond-more" style="margin-top:.5rem">' + (UI.allCond ? 'Show the top 14' : 'Show all ' + condAll.length + ' tested condensates') + '</button>',
    };
    // Functional ontology: broad classes, each optionally expanded into its subclasses.
    const subs = ov.functional_subclasses || [];
    const funcRows = [];
    ov.functional_classes.slice().sort((a, b) => b[4] - a[4]).forEach((r) => {
      funcRows.push({ label: r[0], value: r[4], q: r[5], sig: sigq(r[5]), key: r[0], kind: 'broad', extra: fmt(r[1]),
        title: '<b>' + esc(r[0]) + '</b><br>' + fmt(r[1]) + ' methylarginine proteins (' + (100 * r[2]).toFixed(1) + '%) vs ' + (100 * r[3]).toFixed(1) + '% of reviewed background<br>log2 OR ' + fmtL2(r[4]) + ', q ' + fmtP(r[5]) });
      if (UI.funcExpanded) {
        subs.filter((x) => x[6] === r[0] && x[0] !== r[0]).sort((a, b) => b[4] - a[4]).forEach((x) => funcRows.push({ label: x[0], value: x[4], q: x[5], sig: sigq(x[5]), key: x[0], kind: 'sub', sub: true, extra: fmt(x[1]),
          title: '<b>' + esc(x[0]) + '</b> (subclass of ' + esc(r[0]) + ')<br>' + fmt(x[1]) + ' methylarginine proteins<br>log2 OR ' + fmtL2(x[4]) + ', q ' + fmtP(x[5]) }));
      }
    });
    const func = {
      id: 'ch-func', title: 'Functional ontology', subHtml: 'Functional classes of methylated proteins vs the reviewed human proteome.' + docLink('ov-ontology'),
      subText: 'Functional classes' + (UI.funcExpanded ? ' and subclasses' : '') + ' vs the reviewed human background', rows: funcRows, extraCol: true, extraName: 'proteins',
      footHtml: '<button class="linkish" type="button" id="func-more" style="margin-top:.5rem">' + (UI.funcExpanded ? 'Hide subclasses' : 'Show ' + subs.length + ' subclasses') + '</button>',
    };
    const dis = {
      id: 'ch-dis', title: 'Disease families', subHtml: 'DisGeNET disease families, methylated proteins vs the reviewed human proteome. Association, not causation.' + docLink('ov-disease'),
      subText: 'Audited DisGeNET disease families vs the reviewed background', axis: 'log2 odds ratio vs reviewed background', extraCol: true, extraName: 'methylarginine_genes_of_family',
      rows: ov.disease_families.map((r) => ({ label: r[0], value: r[3] ? Math.log2(r[3]) : null, q: r[4], sig: sigq(r[4]), key: r[0], extra: fmt(r[1]) + '/' + compact(r[2]),
        title: '<b>' + esc(r[0]) + '</b><br>' + fmt(r[1]) + ' methylarginine genes of ' + fmt(r[2]) + ' family genes<br>OR ' + (r[3] || 0).toFixed(2) + ', BH-adjusted p ' + fmtP(r[4]) })).sort((a, b) => b.value - a.value),
    };
    const ctxRows = (table, labels, help) => table.map((r) => ({ label: labels[r[0]] || r[0], value: r[4], q: r[5], sig: sigq(r[5]), key: r[0], extra: fmt(r[1]),
      title: '<b>' + esc(labels[r[0]] || r[0]) + '</b>' + (help[r[0]] ? ' · ' + esc(help[r[0]]) : '') + '<br>' + fmt(r[1]) + ' sites (' + (100 * r[2]).toFixed(1) + '%) vs ' + (100 * r[3]).toFixed(1) + '% of unmodified R in the same proteins<br>log2 OR ' + fmtL2(r[4]) + ', q ' + fmtP(r[5]) }));
    const disorder = {
      id: 'ch-disorder', title: 'Disorder context', subHtml: 'MobiDB disorder, methylated vs other arginines of the same proteins.' + docLink('ov-disorder'),
      subText: 'Methylated vs unmodified arginines in the same proteins (MobiDB)', extraCol: true, extraName: 'sites',
      rows: ctxRows(ov.disorder.filter((r) => r[0] !== 'no_disorder_annotation'), DISORDER_LABEL, DISORDER_HELP),
    };
    const domctx = {
      id: 'ch-domctx', title: 'Domain context', subHtml: 'Position relative to InterPro domains, methylated vs other arginines of the same proteins.' + docLink('ov-domains'),
      subText: 'Position relative to InterPro domains, methylated vs unmodified arginines', extraCol: true, extraName: 'sites',
      rows: ctxRows(ov.domain_context, DOMAIN_LABEL, DOMAIN_HELP),
    };
    const domcls = {
      id: 'ch-domcls', title: 'Nearest domain class', subHtml: 'Class of the nearest InterPro domain.',
      subText: 'Class of the nearest InterPro domain, methylated vs unmodified arginines', extraCol: true, extraName: 'sites',
      rows: ctxRows(ov.domain_classes || [], {}, {}).sort((a, b) => b.value - a.value),
    };

    const multi = M.n_sites - repDist[0];
    el.innerHTML =
      '<p class="ov-hint">Click any bar to filter the Proteins and Sites views; the charts always show the whole dataset.</p>' +
      '<div class="ov-grid">' +
      chartCard(func) + chartCard(dis) +
      ixOverviewHtml() +
      chartCard(motif) + chartCard(cond) + chartCard(disorder) +
      '<div class="card chart-card" id="ch-dom"><div class="card-head"><div><h3>Domains</h3><p class="sub">' + domctx.subHtml + '</p></div>' + expButtons('ch-domctx') + '</div>' + hbars(domctx) +
      '<div class="card-head" style="margin-top:1.1rem"><div><h3 class="h4">Nearest domain class</h3><p class="sub">' + domcls.subHtml + '</p></div>' + expButtons('ch-domcls') + '</div><div id="ch-domcls-bars">' + hbars(domcls) + '</div></div>' +
      '<div class="card wide chart-card" id="ch-study"><div class="card-head"><div><h3>Provenance</h3><p class="sub">' + study.subHtml + '</p></div>' + expButtons('ch-study') + '</div>' +
      '<div class="two-col"><div>' + hbars(study) + '</div>' +
      '<div><p class="flabel" style="margin:0 0 .5rem">Independent reports per site</p><div class="hbars rep-dist">' +
      repDist.map((n, k) => '<button class="hb" type="button" data-rep="' + (k + 1) + '"><span class="hb-label">' + (k === 4 ? '5 or more' : (k + 1) + (k ? ' reports' : ' report')) + '</span><span class="hb-track"><span class="hb-bar" style="left:0;width:' + (100 * n / repDist[0]).toFixed(1) + '%"></span></span><span class="hb-val">' + fmt(n) + '</span></button>').join('') +
      '</div><p class="note">' + fmt(multi) + ' sites (' + pct(multi, M.n_sites) + '%) have two or more. Click a row to keep sites with at least that many.</p></div></div></div>' +
      '<div class="card wide" id="ch-treat"><h3>PRMT inhibitor response</h3><p class="sub">Maron et al. 2021: sites that changed after GSK591 (PRMT5) or MS023 (type I PRMTs) vs DMSO, adj. p ≤ 0.05.' + docLink('ov-inhibitors') + '</p>' + treatTable(ov.treatments) + '</div>' +
      '<div class="card wide starts" id="ov-start"><h3>Starting points</h3><div class="starts-grid">' +
      '<div><h4>Try a search</h4><div class="gene-chips">' + examples.map((e) => '<button class="gene-chip" type="button" data-q="' + esc(e[0]) + '">' + esc(e[0]) + '<span>' + esc(e[1]) + '</span></button>').join('') + '</div></div>' +
      '<div><h4>Most methyl-R sites</h4>' + geneChips(topSites, (p) => p.sites.length) + '</div>' +
      '<div><h4>Strongest arginine-methylation enrichment</h4>' + geneChips(topArg, (p) => fmtL2(p.shrink)) + '</div>' +
      '<div><h4>Densest methyl-R clusters</h4>' + geneChips(topDense, (p) => { const w = p.sig.slice().sort((a, b) => b[4] - a[4])[0]; return w[3] + '/' + w[0] + ' aa'; }) + '</div>' +
      '</div></div>' +
      '</div>';
    bindDocLinks(el);
    UI.charts['ch-study'] = study; UI.charts['ch-domctx'] = domctx; UI.charts['ch-domcls'] = domcls;

    $$('.gene-chip[data-i]', el).forEach((b) => b.addEventListener('click', () => openProtein(+b.dataset.i)));
    $$('.gene-chip[data-q]', el).forEach((b) => b.addEventListener('click', () => {
      setPlain(b.dataset.q); compute();
      if (R.prots.length === 1) openProtein(R.prots[0].i, Q.positions ? Array.from(Q.positions.values())[0].values().next().value : null);
      else go('proteins');
    }));
    bindBars($('#ch-study .hbars:not(.rep-dist)', el), study.rows, (r) => setTerm('study', termSlug('study', r.key), 'sites'));
    $$('.rep-dist .hb', el).forEach((b) => b.addEventListener('click', () => setFilters((f) => { f.minRep = +b.dataset.rep > 1 ? +b.dataset.rep : 0; }, 'sites')));
    bindBars($('#ch-motif', el), motif.rows, (r) => setTerm('motif', motifTerm(r.key), 'sites'));
    const omf = $('#ov-motif', el);
    if (omf) {
      omf.addEventListener('submit', (e) => e.preventDefault());
      bindMotifBox($('#ov-motifq', el), $('#ov-motifwhere', el), $('button[type=submit]', omf), $('#ov-motifmsg', el), 'sites');
    }
    bindBars($('#ch-cond', el), cond.rows, (r) => { F.condConf = thr; setTerm('condensate', termSlug('condensate', +r.key), 'proteins'); });
    const ct = $('#cond-thr', el);
    if (ct) ct.addEventListener('change', () => { UI.condThr = +ct.value; UI.allCond = false; renderOverview(el); });
    $('#cond-more', el).addEventListener('click', () => { UI.allCond = !UI.allCond; renderOverview(el); });
    bindBars($('#ch-func', el), func.rows, (r) => (r.kind === 'sub' ? setTerm('subclass', termSlug('subclass', M.func_subclasses.indexOf(r.key)), 'proteins')
      : setTerm('class', termSlug('class', M.func_classes.indexOf(r.key)), 'proteins')));
    $('#func-more', el).addEventListener('click', () => { UI.funcExpanded = !UI.funcExpanded; renderOverview(el); });
    bindBars($('#ch-dis', el), dis.rows, (r) => setTerm('disease', termSlug('disease', M.disease_families.indexOf(r.key)), 'proteins'));
    bindBars($('#ch-disorder', el), disorder.rows, (r) => setFilters((f) => { const k = M.disorder_classes.indexOf(r.key); f.disorder = new Set(k >= 0 ? [k] : []); }, 'sites'));
    bindBars($('#ch-dom > .hbars', el), domctx.rows, (r) => setFilters((f) => { const k = M.domain_contexts.indexOf(r.key); f.domctx = new Set(k >= 0 ? [k] : []); }, 'sites'));
    bindBars($('#ch-domcls-bars', el), domcls.rows, (r) => setTerm('domain', termSlug('domain', M.domain_classes.indexOf(r.key)), 'sites'));
    $$('#ch-treat tr.row', el).forEach((tr) => {
      const go2 = () => setFilters((f) => { f.treat = tr.dataset.t + ':' + tr.dataset.dir; f.ips = new Set([+tr.dataset.ip]); }, 'sites');
      tr.addEventListener('click', go2);
      tr.addEventListener('keydown', (e) => { if (e.key === 'Enter') go2(); });
    });
    bindIxGo(el);
    bindCheckSite(el);
    bindChartExports(el);
  }

  function treatTable(rows) {
    const max = Math.max(...rows.map((r) => (r[2] ? Math.max(r[3], r[4]) / r[2] : 0)), 0.01);
    const order = rows.slice().sort((a, b) => a[0].localeCompare(b[0]) || M.maron_ips.indexOf(a[1]) - M.maron_ips.indexOf(b[1]));
    return '<div class="tbl-wrap"><table class="tbl treat-tbl"><thead><tr><th>Inhibitor</th><th>Methyl-IP</th><th class="r">Sites</th><th>Decreased</th><th class="hide-xs">Increased</th><th class="r hide-sm">Median change</th></tr></thead><tbody>' +
      order.map((r) => {
        const t = M.maron_treatments.indexOf(r[0]);
        const ip = M.maron_ips.indexOf(r[1]);
        const down = r[2] ? r[3] / r[2] : 0, up = r[2] ? r[4] / r[2] : 0;
        return '<tr class="row" tabindex="0" data-t="' + t + '" data-ip="' + ip + '" data-dir="down"><td><b>' + esc(r[0]) + '</b> <span class="muted">' + esc(TREAT_HELP[r[0]] || '') + '</span></td><td><span class="tag ip" title="' + esc(IP_HELP[r[1]] || '') + '">' + esc(ipLabel(r[1])) + '</span></td><td class="r num">' + fmt(r[2]) + '</td>' +
          '<td><div class="meter-row"><span class="meter"><i style="width:' + (100 * down / max).toFixed(1) + '%"></i></span><span class="num">' + r[3] + ' (' + (100 * down).toFixed(0) + '%)</span></div></td>' +
          '<td class="hide-xs"><div class="meter-row"><span class="meter up"><i style="width:' + (100 * up / max).toFixed(1) + '%"></i></span><span class="num">' + r[4] + ' (' + (100 * up).toFixed(0) + '%)</span></div></td>' +
          '<td class="r num hide-sm">' + fmtL2(r[5]) + '</td></tr>';
      }).join('') + '</tbody></table></div><p class="note">Median change: Table S3 <code>&lt;inhibitor&gt;_vs_DMSO_diff</code> across sites. Click a row for the sites.</p>';
  }

  // --------------------------------------------------------------- protein
  function renderProtein(el, p) {
    if (!p) { go('overview'); return; }
    const shownN = R.counts[p.i];
    const argN = p.argN;
    const best = p.best;
    const nnMed = p.dens ? p.dens[1] : null;
    const ext = (href, label) => '<a class="btn quiet small" href="' + href + '" target="_blank" rel="noopener">' + label + ' ↗</a>';
    const base = p.acc.split('-')[0];
    const filtered = constraints().length > 0;
    el.innerHTML =
      '<div class="print-only print-head" id="print-head"></div>' +
      crumbsHtml(p) +
      '<div class="detail-head">' +
      '<div class="detail-title"><h2>' + esc(p.label) + '</h2><span class="pn" id="pname">' + esc(p.name) + '</span></div>' +
      (p.aliases.length ? '<p class="aliases">Also reported as ' + p.aliases.map((a) => '<span>' + esc(a) + '</span>').join(', ') + '</p>' : '') +
      '<div class="detail-links"><span class="tag mono">' + esc(p.acc) + '</span>' + (p.uid ? '<span class="tag mono">' + esc(p.uid) + '</span>' : '') +
      (p.iso ? '<span class="tag" title="Sites were reported on this isoform accession">isoform</span>' : '') +
      (p.unrev ? '<span class="tag warn-tag" title="Unreviewed UniProt entry (TrEMBL)">unreviewed entry</span>' : '') +
      ext('https://www.uniprot.org/uniprotkb/' + encodeURIComponent(p.acc) + '/entry', 'UniProt') +
      ext('https://alphafold.ebi.ac.uk/entry/' + encodeURIComponent(base), 'AlphaFold') +
      ext('https://www.ebi.ac.uk/interpro/protein/UniProt/' + encodeURIComponent(base) + '/', 'InterPro') +
      '<button class="btn small" type="button" id="save-pdf" title="Print this protein page, or save it as a PDF from the print dialog">Save PDF</button></div></div>' +
      ctxToggleHtml(p) +
      (filtered && shownN < p.sites.length ? '<div class="filter-note">' + fmt(shownN) + ' of ' + fmt(p.sites.length) + ' sites match your filters; the others are faded. <button class="linkish" type="button" id="pf-clear">Clear filters</button></div>' : '') +
      (p.unrev ? (() => {
        const twin = p.reviewedTwin, nDup = p.sites.filter((j) => S[j].dup != null).length;
        return '<div class="filter-note">This is an unreviewed UniProt entry (TrEMBL), so sequence-based context (motifs, disorder, domains) may be missing. ' +
          (twin ? 'The reviewed entry for ' + esc(p.gene) + ' is <button class="linkish" type="button" data-open-p="' + twin.i + '">' + esc(twin.label) + ' (' + esc(twin.acc) + ')</button>; ' + (p.sites.length === 1 ? (nDup ? 'the site here is also a site there.' : 'the site here is not among its sites.') : (nDup === p.sites.length ? 'all ' + fmt(nDup) + ' sites here are also sites there.' : fmt(nDup) + ' of the ' + fmt(p.sites.length) + ' sites here are also sites there, at the same position.')) : 'No reviewed entry for this gene is in the set.') + '</div>';
      })() : '') +
      '<div class="facts">' +
      fact(fmt(p.sites.length), 'Methylarginine sites', shownN === p.sites.length ? (filtered ? 'all match your filters' : 'all reported sites') : fmt(shownN) + ' match current filters') +
      fact('<span id="arg-count">' + (argN == null ? '—' : fmt(argN)) + '</span>', 'Arginines', '<span id="arg-frac">' + (argN ? pctText(p.sites.length, argN) + '% methylated' : 'loading from UniProt') + '</span>') +
      fact(p.shrink == null ? '—' : fmtL2(p.shrink), 'Arg enrichment (log2)', p.argQ == null ? 'not tested (needs sequence-confirmed sites)' : 'shrinkage estimate; Fisher q ' + fmtP(p.argQ)) +
      fact(best ? best[3] + ' / ' + best[0] + ' aa' : '—', 'Densest window', best ? rl(best[1], p) + '–' + rn(best[2], p) + ', q ' + fmtP(best[5]) + (best[5] <= PADJ ? '' : ' (n.s.)') : 'needs ≥3 sequence-confirmed sites') +
      fact('<span id="af-low">…</span>', 'Sites in low-confidence model', '<span id="af-low-sub">AlphaFold pLDDT &lt; 70</span>') +
      fact(fmt(p.cond.length), 'CD-CODE condensates', '<span id="cond-fact-sub">' + (p.cond.length ? '' : 'none annotated') + '</span>') +
      '</div>' +
      figureCard(p) +
      '<div class="card" id="site-card" hidden></div>' +
      contactCardHtml() +
      '<div class="card tbl-card" id="card-sites"><div class="card-pad"><div class="card-head"><div><h3>Sites</h3><p class="sub">Every reported methylarginine on ' + esc(p.label) + '. Click a row for a summary of the site; faded rows do not match the current filters.</p></div><div class="actions" id="psites-tools"></div></div><div id="psites-ix"></div></div>' +
      '<div class="tbl-wrap" id="psites"></div><p class="note uni-extra" id="uni-extra" hidden></p><div class="card-pad" id="psites-cite"></div></div>' +
      '<div class="card" id="card-ix"></div>' +
      '<div class="card" id="card-cons"></div>' +
      structureCard(p) +
      '<div class="two-col"><div class="card" id="card-dom"></div><div class="card" id="card-cond"></div></div>' +
      '<div class="two-col"><div class="card" id="card-func"></div><div class="card" id="card-dis"></div></div>' +
      '<div class="two-col"><div class="card" id="card-cluster"></div><div class="card" id="card-ptm"></div></div>' +
      '<div class="card"><div class="card-head"><div><h3>Sequence</h3><p class="sub" id="seq-sub">Loading sequence from UniProt…</p></div><div class="exp" role="group" aria-label="Sequence export"><button type="button" id="copy-fasta" disabled>Copy FASTA</button><button type="button" id="seq-svg" disabled title="Sequence figure with gene, protein name, UniProt accession and legend">SVG</button><button type="button" id="seq-png" disabled title="Sequence figure with gene, protein name, UniProt accession and legend">PNG</button></div></div><div id="seq"></div></div>' +
      '<div class="print-only print-foot"><p>' + esc(siteName()) + ' · ' + esc(siteUrl()) + '</p><p>Cite: ' + esc(PAPER.cite) + ' https://doi.org/' + esc(PAPER.doi) + '</p><p>Sites: methylarginine union from PTM_Pipeline. Features, variants and sequence: UniProt. Predicted structure, pLDDT and AlphaMissense: AlphaFold DB. Experimental structures: RCSB PDB and PDBe (SIFTS). Condensates: CD-CODE.</p></div>';

    $$('[data-crumb]', el).forEach((b) => b.addEventListener('click', () => { go(b.dataset.crumb); window.scrollTo({ top: 0 }); }));
    $('#save-pdf', el).addEventListener('click', (e) => savePdf(p, e.currentTarget));
    bindCtx(p);
    $$('[data-open-p]', el).forEach((b) => b.addEventListener('click', () => openProtein(+b.dataset.openP)));
    const pfc = $('#pf-clear', el);
    if (pfc) pfc.addEventListener('click', clearFilters);
    renderSiteTable(p);
    renderDomainCard(p, null);
    renderCondCard(p);
    renderFuncCard(p);
    renderDiseaseCard(p, null);
    renderClusterCard(p);
    renderPtmCard(p);
    renderIxCard(p);
    renderConsCard(p, UNIDATA.has(p.acc) ? UNIDATA.get(p.acc) : undefined);
    bindFigure(p);
    drawMap(p, UNIDATA.get(p.acc));
    initStructure(p);
    const current = () => UI.detail === p.i && UI.view === 'protein';
    loadEvidence().then((ev) => { if (current()) renderDiseaseCard(p, ev); });
    loadUniProt(p.acc).then((u) => {
      if (!current()) return;
      drawMap(p, u);
      loadEvidence().then((ev) => { if (current()) renderDiseaseCard(p, ev); });
      renderSequence(p, u);
      renderConsCard(p, u);
      renderDomainCard(p, u);
      renderSiteTable(p);
      renderContactCard(p);
      if (u) {
        if (!p.name && u.name) $('#pname').textContent = u.name;
        if (argN == null && u.seq) {
          const n = (u.seq.match(/R/g) || []).length;
          $('#arg-count').textContent = fmt(n);
          $('#arg-frac').textContent = pctText(p.sites.length, n) + '% methylated (UniProt sequence)';
        }
      } else if (argN == null) $('#arg-frac').textContent = 'UniProt unavailable';
    });
    loadIface(p.acc).then(() => { if (current()) drawMap(p, UNIDATA.get(p.acc)); });
    loadAF(p.acc).then((af) => {
      if (!current()) return;
      const lowEl = $('#af-low'), subEl = $('#af-low-sub');
      if (af && af.plddt.size) {
        const vals = p.sites.map((j) => af.plddt.get(S[j].pos)).filter((v) => v != null);
        const low = vals.filter((v) => v < 70).length;
        lowEl.textContent = pct(low, vals.length) + '%';
        subEl.textContent = low + ' of ' + vals.length + ' sites at pLDDT < 70 (model mean ' + (af.entry.globalMetricValue != null ? af.entry.globalMetricValue.toFixed(1) : '—') + ')';
      } else { lowEl.textContent = '—'; subEl.textContent = 'no AlphaFold DB model'; }
      renderSiteTable(p);
    });
    if (UI.highlight != null) setTimeout(() => flashSite(UI.highlight, true), 60);
    if (UI.pendingContact) { const pc = UI.pendingContact; UI.pendingContact = null; setTimeout(() => openContactStructure(p, pc.pdb, pc.pos), 400); }
    if (UI.pendingModel) { const g = UI.pendingModel; UI.pendingModel = null; setTimeout(() => viewGalleryModel(p, g), 400); }
  }

  function fact(value, label, sub) {
    return '<div class="fact"><b>' + value + '</b><span>' + label + '</span>' + (sub ? '<small>' + sub + '</small>' : '') + '</div>';
  }

  const shortDisease = (desc) => { const m = /^in ([^;]+)/.exec(desc || ''); return m ? m[1].replace(/\bdbSNP.*$/, '').trim() : ''; };

  function afCell(p, s) {
    const af = AFDATA.get(p.acc);
    if (af === undefined) return '<span class="muted">…</span>';
    if (!af) return '<span class="muted" title="No AlphaFold DB model">—</span>';
    const v = af.plddt.get(s.pos);
    const a = af.am.get(s.pos);
    const band = v == null ? -1 : plddtBand(v);
    return (v == null ? '<span class="muted">—</span>' : '<span class="plddt" title="AlphaFold pLDDT ' + v.toFixed(1) + ' · ' + esc(PLDDT_BANDS[band][2]) + '"><i style="background:' + PLDDT_BANDS[band][1] + '"></i>' + Math.round(v) + '</span>') +
      (a && a.n ? ' <span class="am am-' + amClass(a.sum / a.n).split(' ')[1] + '" title="Mean AlphaMissense pathogenicity of ' + a.n + ' R→X substitutions: ' + (a.sum / a.n).toFixed(2) + ' (' + amClass(a.sum / a.n) + '); ' + a.lp + ' of ' + a.n + ' likely pathogenic">AM ' + (a.sum / a.n).toFixed(2) + '</span>' : '');
  }

  // Sites table of a protein page. Row actions: 3D zooms the structure to the site, Evidence expands the source rows.
  function psiteCols(p, u, varsAt) {
    return [
      { k: 'site', label: 'Site', tdCls: 'num site-col', td: (s) => {
        const vs = (varsAt.get(s.pos) || []).filter((v) => v.from === 'R');
        return '<b>' + rl(s.pos, P[s.p]) + '</b>' + (p.hist ? '<br><span class="muted small">UniProt R' + s.pos + '</span>' : '') +
          (s.dup != null ? '<br><span class="tag" title="The same gene and position is also a site on the reviewed UniProt entry ' + esc(P[S[s.dup].p].acc) + '">also on ' + esc(P[S[s.dup].p].acc) + '</span>' : '') +
          (s.orig ? '<br><span class="muted small" title="Positions given by the source studies before mapping to the current canonical UniProt sequence. Several can land on one site when a source used another isoform or its peptide fits more than one arginine; expand the row to see how each record was mapped.">source position' + (new Set(s.orig.split(';')).size > 1 ? 's ' : ' ') + esc(Array.from(new Set(s.orig.split(';'))).join(', ')) + '</span>' : '') +
          vs.map((v) => '<br><span class="tag var" title="UniProt natural variant ' + esc(v.id) + ': ' + rl(s.pos, P[s.p]) + v.to + ' · ' + esc(v.desc) + '">' + rl(s.pos, P[s.p]) + esc(v.to) + (shortDisease(v.desc) ? ' · ' + esc(shortDisease(v.desc).slice(0, 18)) : '') + '</span>').join('');
      } },
      { k: 'state', label: 'State', td: (s) => pill(s.st) },
      { k: 'win', label: 'Window ±7 · motifs', td: winCell },
      { k: 'ix', label: 'Interaction evidence', opt: true, tdCls: 'ix-col', tip: 'Experimental: RNA/DNA crosslink, RBDmap or a PDB interface at this residue. Predicted: at the interface of a predicted complex.', td: (s) => ixCell(s, 3) },
      { k: 'prmti', label: 'PRMTi', opt: true, menu: 'PRMT inhibitor response', td: treatHtml },
      { k: 'src', label: 'Reported by', opt: true, cls: 'hide-sm', td: (s) => '<div class="tags">' + studyTags(s.study, 3) + '</div>' },
      { k: 'af', label: 'AlphaFold', opt: true, menu: 'AlphaFold pLDDT and AlphaMissense', tdCls: 'nw', tip: 'AlphaFold pLDDT at the site and mean AlphaMissense score for R→X', td: (s) => afCell(p, s) },
      { k: 'ctx', label: 'Context', opt: true, menu: 'Disorder and domain context', tdCls: 'nw', td: contextHtml },
      { k: 'dom', label: 'Nearest domain', opt: true, td: (s) => {
        const ipr = s.ipr >= 0 ? M.interpro[s.ipr] : null;
        return ipr ? '<a href="https://www.ebi.ac.uk/interpro/entry/InterPro/' + encodeURIComponent(ipr[0]) + '/" target="_blank" rel="noopener" title="' + esc(ipr[0]) + ' · ' + esc(M.domain_classes[s.dcls]) + '">' + esc(ipr[1]) + '</a><br><span class="muted">' + (s.domDist === 0 ? 'inside' : s.domDist + ' aa ' + (s.domSide === 'N_terminal_to_domain' ? 'N-term' : 'C-term')) + '</span>' : '<span class="muted">—</span>';
      } },
      { k: 'nn', label: 'Next methyl-R', opt: true, cls: 'r', tdCls: 'num', td: (s) => (s.nn == null ? '—' : s.nn + ' aa') },
      { k: 'stev', label: 'State evidence', opt: true, menu: 'State evidence (Maron IP, Table S2, UniProt)', tip: 'Maron methyl-IP and Table S2 calls, and UniProt\'s curated modification', td: (s) => {
        const s2 = bits(s.s2).map((b) => ['Rme1', 'Rme2'][b]).join('/'), uniType = u && u.methylR.get(s.pos);
        return '<div class="tags">' + ([ipTags(s.ip), s2 ? '<span class="tag">S2 ' + s2 + '</span>' : '', uniType ? '<span class="tag uni" title="UniProt: ' + esc(uniType) + '">UniProt ' + esc(uniprotMethylLabel(uniType)) + '</span>' : ''].filter(Boolean).join(' ') || '<span class="muted">—</span>') + '</div>';
      } },
      { k: 'act', label: '<span class="sr-only">Actions</span>', tdCls: 'nw', td: (s) => { const open = UI.openEvidence.has(s.j); return '<button class="btn quiet small view3d" type="button" data-pos="' + s.pos + '" title="Zoom to ' + rl(s.pos, P[s.p]) + ' in the 3D structure">3D</button> <button class="btn quiet small ev-toggle" type="button" aria-expanded="' + open + '" data-j="' + s.j + '" title="Source records, mapping, variants and AlphaMissense">' + (open ? 'Hide' : 'Evidence') + '</button>'; } },
    ];
  }
  // Interaction status counts for one protein's sites; a click shows only those sites and colours the map by status.
  function pIxSegHtml(p) {
    if (!ixReady()) return '<p class="note">Loading interaction evidence…</p>';
    const n = { exp: 0, pred: 0, none: 0 };
    for (const j of p.sites) n[ixStatus(S[j])]++;
    const b = (k) => '<button type="button" class="seg-btn" data-pix="' + k + '" aria-pressed="' + (UI.pIx === k) + '"' + (n[k] ? '' : ' disabled') + '><i class="ixd ixd-' + k + '"></i>' + IXS[k][0] + ' <b>' + n[k] + '</b></button>';
    return '<div class="ix-seg"><span class="flabel">Interaction evidence at these sites' + docLink('ix-status') + '</span><div class="ixseg-btns">' + b('exp') + b('pred') + b('none') + '</div></div>';
  }
  // Zoom the 3D structure to a site.
  function show3d(pos) {
    UI.highlight = pos;
    const card = $('#card-af');
    if (card) card.scrollIntoView({ block: 'start', behavior: 'smooth' });
    if (UI.startStructure) UI.startStructure();
    flashSite(pos, false);
  }

  function renderSiteTable(p) {
    const box = $('#psites');
    if (!box) return;
    const u = UNIDATA.get(p.acc);
    const rows = p.sites.map((j) => S[j]);
    const varsAt = new Map();
    if (u) for (const v of u.variants) if (v.start === v.end) { if (!varsAt.has(v.start)) varsAt.set(v.start, []); varsAt.get(v.start).push(v); }
    const cols = colsShown('psites', psiteCols(p, u, varsAt));
    const tools = $('#psites-tools');
    if (tools) tools.innerHTML = colMenu('psites', psiteCols(p, u, varsAt));
    const seg = $('#psites-ix');
    if (seg) seg.innerHTML = pIxSegHtml(p);
    const ixF = UI.pIx || '';
    box.innerHTML = '<table class="tbl stackable">' + colHead(cols, { dir: 1 }, null) + '<tbody>' +
      rows.map((s) => {
        const dim = !R.shown.has(s.j), open = UI.openEvidence.has(s.j);
        if (ixF && ixStatus(s) !== ixF) return '';
        return '<tr class="row' + (dim ? ' dim' : '') + (UI.highlight === s.pos ? ' sel' : '') + '" tabindex="0" data-j="' + s.j + '" data-pos="' + s.pos + '">' + colRow(cols, s) + '</tr>' +
          (open ? '<tr class="ev-row" data-for="' + s.j + '"><td colspan="' + cols.length + '" class="ev" id="ev-' + s.j + '">Loading evidence…</td></tr>' : '');
      }).join('') + '</tbody></table>' +
      (ixF ? '<p class="note">Showing sites with ' + esc(IXS[ixF][0].toLowerCase()) + ' interaction evidence. <button class="linkish" type="button" data-pix="">Show all sites</button></p>' : '');
    bindRows(box, (tr) => { UI.highlight = +tr.dataset.pos; flashSite(+tr.dataset.pos, false); });
    bindContactLinks(box);
    $$('.ev-toggle', box).forEach((b) => b.addEventListener('click', () => {
      const j = +b.dataset.j;
      if (UI.openEvidence.has(j)) UI.openEvidence.delete(j); else UI.openEvidence.add(j);
      renderSiteTable(p);
    }));
    $$('.view3d', box).forEach((b) => b.addEventListener('click', () => show3d(+b.dataset.pos)));
    const card = box.closest('.card');
    if (card) {
      bindColMenu(card, 'psites', () => { renderSiteTable(p); const m = $('.col-menu', card); if (m) m.open = true; });
      $$('[data-pix]', card).forEach((b) => b.addEventListener('click', () => {
        UI.pIx = UI.pIx === b.dataset.pix ? '' : b.dataset.pix;
        if (UI.pIx && UI.fig.color !== 'ixstat') { UI.fig.color = 'ixstat'; const fc = $('#fig-color'); if (fc) fc.value = 'ixstat'; drawMap(p, UNIDATA.get(p.acc)); }
        renderSiteTable(p);
      }));
    }
    if (UI.openEvidence.size) loadEvidence().then((ev) => UI.openEvidence.forEach((j) => fillEvidence(j, ev)));
    // Sources to cite for this protein's sites; the evidence part fills in once the interaction data is here.
    const pc = $('#psites-cite');
    if (pc) {
      const wasOpen = !!$('details[open]', pc);
      pc.innerHTML = citeBoxHtml(rows, 'the ' + rows.length + ' site' + (rows.length === 1 ? '' : 's') + ' on ' + p.label, wasOpen);
      bindCiteBox(pc);
    }
    // The site card shows the same per-site data (AlphaFold, variants), so it follows each refresh.
    const sc = $('#site-card');
    if (sc && !sc.hidden && UI.highlight != null) renderSiteCard(UI.highlight);
    // UniProt-annotated methylarginines that the union does not contain: a pointer to coverage gaps.
    const extra = $('#uni-extra');
    if (extra && u) {
      const have = new Set(p.sites.map((j) => S[j].pos));
      const missing = Array.from(u.methylR.entries()).filter(([pos]) => !have.has(pos)).sort((a, b) => a[0] - b[0]);
      extra.hidden = !missing.length;
      extra.innerHTML = missing.length ? 'UniProt also annotates ' + missing.length + ' methylarginine' + (missing.length === 1 ? '' : 's') + ' that ' + (missing.length === 1 ? 'is' : 'are') + ' not in the union: ' +
        missing.map(([pos, d]) => '<span class="tag uni" title="' + esc(d) + '">' + rl(pos) + ' ' + esc(uniprotMethylLabel(d)) + '</span>').join(' ') + '. Positions follow the current UniProt sequence.' : '';
    }
    const gap = ixReady() ? IX.modelledGaps.get(p.i) : null;
    if (extra && gap) {
      extra.hidden = false;
      extra.innerHTML += (extra.innerHTML ? '<br>' : '') + 'PDB entries model a methylated arginine at ' + gap.pos.map((x) => '<span class="tag var">' + rl(x) + '</span>').join(' ') + ', which the union does not contain (' +
        gap.pdbs.map((x) => '<a href="https://www.rcsb.org/structure/' + x.toUpperCase() + '" target="_blank" rel="noopener">' + x.toUpperCase() + '</a>').join(', ') + ').';
    }
  }

  function fillEvidence(j, ev) {
    const cellEl = document.getElementById('ev-' + j);
    if (!cellEl) return;
    const s = S[j];
    const p = P[s.p];
    const rows = ev ? ev.sites[j] : null;
    const u = UNIDATA.get(p.acc);
    const af = AFDATA.get(p.acc);
    let html = '';
    if (s.tr.length) {
      html += '<p><b>PRMT inhibitor response</b> (Maron et al. 2021 Table S3, change vs DMSO, adj. p)</p><table><thead><tr><th>Methyl-IP</th><th>Inhibitor</th><th>Change</th><th>adj. p</th><th>Peptides</th></tr></thead><tbody>' +
        s.tr.slice().sort((a, b) => a[0] - b[0] || a[1] - b[1]).map(([ip, t, d, q, n]) => '<tr><td>' + ipLabel(M.maron_ips[ip]) + '</td><td>' + M.maron_treatments[t] + ' <span class="muted">' + TREAT_HELP[M.maron_treatments[t]] + '</span></td><td class="num">' + fmtL2(d) + '</td><td class="num">' + (q != null && q <= PADJ ? '<b>' + fmtP(q) + '</b>' : fmtP(q)) + '</td><td class="num">' + (n == null ? '—' : n) + '</td></tr>').join('') + '</tbody></table>';
    }
    if (u) {
      const at = (list) => list.filter((f) => f.start <= s.pos && f.end >= s.pos);
      const vs = at(u.variants), mg = at(u.mutagen);
      const uni = u.methylR.get(s.pos);
      if (uni || vs.length || mg.length) {
        html += '<p><b>UniProt at ' + rl(s.pos, P[s.p]) + '</b></p><ul class="term-list">' +
          (uni ? '<li>Modified residue: ' + esc(uni) + '</li>' : '') +
          vs.map((v) => '<li>Natural variant ' + esc(v.from) + '→' + esc(v.to) + ': ' + esc(v.desc || '—') + (v.id ? ' <span class="muted">(' + esc(v.id) + ')</span>' : '') + '</li>').join('') +
          mg.map((v) => '<li>Mutagenesis ' + esc(v.from) + '→' + esc(v.to) + (v.start !== v.end ? ' (' + v.start + '–' + v.end + ')' : '') + ': ' + esc(v.desc || '—') + '</li>').join('') + '</ul>';
      }
    }
    if (af) {
      const v = af.plddt.get(s.pos), a = af.am.get(s.pos);
      if (v != null || (a && a.n)) {
        html += '<p><b>AlphaFold</b> (' + esc(af.entry.entryId) + '): ' + (v != null ? 'pLDDT ' + v.toFixed(1) + ' (' + PLDDT_BANDS[plddtBand(v)][2].toLowerCase() + ')' : '') +
          (a && a.n ? '; AlphaMissense mean ' + (a.sum / a.n).toFixed(3) + ' over ' + a.n + ' R→X substitutions (' + amClass(a.sum / a.n) + '), max ' + a.max.toFixed(3) + ', ' + a.lp + ' likely pathogenic' : '') + '.</p>';
      }
    }
    const near = [];
    for (const [pos, g, res] of p.ptms) { const d = pos - s.pos; if (Math.abs(d) <= 10) near.push(PTM_PREFIX[g] + res + pos + ' (' + (d > 0 ? '+' : '') + d + ')'); }
    if (near.length) html += '<p><b>Other PTMs within ±10 aa</b> <span class="muted">(PTM master table; not exhaustive)</span>: ' + esc(near.join(', ')) + '</p>';
    if (s.win) html += '<p><b>Window ±15:</b> ' + windowHtml(s.win, 15) + '</p>';
    if (rows && rows.length) {
      html += '<p><b>Source rows</b> (' + rows.length + ')</p><table><thead><tr><th>Study / database</th><th>Source table</th><th>Reference</th><th>PMID</th><th>Peptide</th><th>Enzyme</th><th>State</th><th>Reported as</th></tr></thead><tbody>' +
        rows.map((r) => { const st = r[9] != null ? ev.studies[r[9]] : null; const ref = r[1] >= 0 ? ev.refs[r[1]] : '';
          return '<tr><td><b>' + (st && STUDY_DOI[st] ? doiLink(STUDY_DOI[st], esc(st)) : esc(st || '—')) + '</b></td><td>' + esc(ev.source_tables[r[0]]) + '</td><td>' + esc(ref === 'This study.' ? 'Maron et al. 2021 (first reported in that study)' : ref || '—') + '</td>' +
          '<td>' + (r[2] ? r[2].split(/[;,\s]+/).filter(Boolean).map((id) => '<a href="https://pubmed.ncbi.nlm.nih.gov/' + encodeURIComponent(id) + '/" target="_blank" rel="noopener">' + esc(id) + '</a>').join(' ') : '—') + '</td>' +
          '<td class="pep">' + esc(r[3] || '—') + '</td><td>' + esc(r[4] || '—') + '</td><td>' + esc(r[5] || '—') + '</td>' +
          '<td class="mono">' + esc(r[6]) + ':' + esc(r[7]) + (r[7] && String(r[7]) !== String(s.pos) ? ' <span class="muted">→ ' + rl(s.pos, P[s.p]) + '</span>' : '') + '<br><span class="muted"' + (REMAP_HELP[ev.remap_statuses[r[8]]] ? ' title="' + esc(ev.remap_statuses[r[8]]) + '"' : '') + '>' + esc((REMAP_HELP[ev.remap_statuses[r[8]]] || ev.remap_statuses[r[8]].replace(/_/g, ' '))) + '</span></td></tr>'; }).join('') +
        '</tbody></table>' +
        (rows.some((r) => r[9] != null && ev.studies[r[9]] === 'PhosphoSitePlus') ? '<p class="note">' + esc(PSP_NOTICE) + '</p>' : '');
    } else if (!ev) html += '<p class="note">The evidence file did not load.</p>';
    if (s.refs.length && !(rows && rows.length)) html += '<p><b>References:</b> ' + s.refs.map((k) => esc(M.refs[k])).join('; ') + '</p>';
    cellEl.innerHTML = html || '<span class="muted">No additional evidence.</span>';
  }

  // Site card: one methylarginine at a glance, interaction evidence first. Opens when a site is searched or clicked.
  function renderSiteCard(pos, scroll) {
    const box = $('#site-card');
    if (!box || UI.detail == null) return;
    const p = P[UI.detail];
    const j = pos == null ? null : p.sites.find((x) => S[x].pos === pos);
    const s = j == null ? null : S[j];
    box.hidden = !s;
    const crumb = $('#crumb-site');
    if (crumb) { crumb.hidden = !s; crumb.textContent = s ? rl(s.pos, p) : ''; }
    if (UI.view === 'protein') setSiteHash(p, s ? s.pos : null);
    if (!s) { box.innerHTML = ''; return; }
    const st = ixStatus(s);
    let ix;
    if (!st) ix = '<p class="note">Loading interaction evidence…</p>';
    else {
      const m = s.ev || 0, pt = ixPartners(s);
      const lines = EV_SETS.filter(([b]) => m & b).map(([b, l, , help]) => '<li><i class="ixd ixd-' + (b & EV_EXP ? 'exp' : 'pred') + '"></i><b>' + esc(l) + '</b> <span class="muted">' + esc(help) + '</span></li>').join('');
      const tags = ifaceTags(s, 8) + modelledTag(s) + predTag(s) + xlTags(s);
      ix = '<p class="sc-status">' + ixBadge(s) + ' ' + (st === 'exp' ? 'This arginine has experimental contact evidence' + (m & EV_PRED ? ' and is also at a predicted interface' : '') + '.'
        : st === 'pred' ? 'Predicted to sit at a protein interface; no experimental contact recorded yet.' : 'No contact evidence is recorded at this residue. Most proteins were never crystallised or crosslinked, so this is not evidence of no contact.') + '</p>' +
        (lines ? '<ul class="sc-ev">' + lines + '</ul>' : '') +
        (pt.exp.length ? '<p class="sc-who"><span class="flabel">Contacts (experimental)</span> ' + pt.exp.map(esc).join(', ') + '</p>' : '') +
        (pt.pred.length ? '<p class="sc-who"><span class="flabel">Partners (predicted)</span> ' + pt.pred.map(esc).join(', ') + '</p>' : '') +
        (tags ? '<div class="tags">' + tags + '</div>' : '') + contactHtml(s, 6);
      // What is known about the whole protein, labelled as such so it is not read as site evidence.
      const part = partnersAt(p, 700), cls = (v) => part.filter(([k]) => bindMatch(v, k)).map(([k]) => ixGene(k));
      const pr = cls(1), tu = cls(2), wd = cls(4), cxN = (p.cx || []).filter((k) => !IX.cx[k].pred).length;
      const go = [p.ixFlags & IXF.rna ? 'RNA' : '', p.ixFlags & IXF.dna ? 'DNA' : ''].filter(Boolean);
      const prot = [part.length ? fmt(part.length) + ' STRING physical partner' + (part.length === 1 ? '' : 's') + ' (≥ 0.7)' : 'no STRING physical partner ≥ 0.7',
        pr.length ? 'PRMTs: ' + pr.join(', ') : '', tu.length ? 'Tudor readers: ' + tu.slice(0, 6).join(', ') + (tu.length > 6 ? ' +' + (tu.length - 6) : '') : '', wd.length ? 'WD40: ' + wd.slice(0, 4).join(', ') : '',
        cxN ? cxN + ' curated complex' + (cxN === 1 ? '' : 'es') : '', go.length ? 'GO ' + go.join(' and ') + ' binding' : ''].filter(Boolean);
      ix += '<p class="sc-prot"><span class="flabel">Whole protein (not specific to ' + rl(s.pos, P[s.p]) + ')</span> ' + esc(prot.join(' · ')) + ' <button class="linkish" type="button" data-sc="ixcard">Interactions card</button></p>';
    }
    const pub = (k, v) => '<div><dt>' + k + '</dt><dd>' + v + '</dd></div>';
    box.innerHTML = '<div class="card-head"><div><span class="eyebrow">Selected site</span><h3>' + esc(p.label) + ' ' + rl(s.pos, P[s.p]) + ' ' + pill(s.st) + '</h3>' +
      (p.hist ? '<p class="muted small">Histone numbering, without the initiator methionine · UniProt R' + s.pos + '</p>' : '') + '</div>' +
      '<div class="actions"><button class="btn quiet small" type="button" data-sc="3d">Show in 3D</button><button class="btn quiet small" type="button" data-sc="ev">Source evidence</button>' +
      '<button class="btn quiet small" type="button" data-sc="link" title="Copy a link that opens this page with ' + rl(s.pos, P[s.p]) + ' selected">Copy link</button><button class="btn quiet small" type="button" data-sc="close" aria-label="Close the site summary">×</button></div></div>' +
      '<div class="sc-grid"><section class="sc-ix"><h4>Interaction evidence' + docLink('ix-status') + '</h4>' + ix + '</section>' +
      '<section><h4>Site</h4><dl class="sc-kv">' + pub('Window ±7', windowHtml(s.win, 7)) + pub('Motifs', motifTags(s.mot)) + pub('Context', contextHtml(s).replace('<br>', ' · ')) +
      pub('Next methyl-R', s.nn == null ? '—' : s.nn + ' aa') + pub('PRMT inhibitors', treatHtml(s)) + pub('Maron IP', ipTags(s.ip) || '—') + pub('AlphaFold', afCell(p, s)) +
      pub('Reported by', '<ul class="sc-refs">' + M.studyOrder.filter((k) => s.study & (1 << k)).map((k) => { const r = studyRef(k); return '<li>' + (r[1] ? '<a href="https://doi.org/' + r[1] + '" target="_blank" rel="noopener" title="' + esc(r[0]) + '">' + esc(M.studies[k][0]) + '</a>' : '<span title="' + esc(r[0]) + '">' + esc(M.studies[k][0]) + '</span>') + '</li>'; }).join('') + '</ul>' +
        '<p class="muted small">Cite these for ' + rl(s.pos, P[s.p]) + ', with ' + esc(PAPER.short) + ' for the compilation.</p>') + '</dl></section></div>';
    bindDocLinks(box);
    bindContactLinks(box);
    $$('[data-sc]', box).forEach((b) => b.addEventListener('click', (e) => {
      const act = b.dataset.sc;
      if (act === '3d') show3d(s.pos);
      else if (act === 'ev') {
        UI.openEvidence.add(s.j);
        renderSiteTable(p);
        const tr = $('#psites tr.row[data-pos="' + s.pos + '"]');
        if (tr) tr.scrollIntoView({ block: 'start', behavior: 'smooth' });
      } else if (act === 'link') {
        const url = location.origin + location.pathname + '?q=' + encodeURIComponent(p.label + ' ' + rl(s.pos, P[s.p])) + '#' + encodeURIComponent(p.acc);
        copyText(url, e.currentTarget);
      } else if (act === 'ixcard') { const c = $('#card-ix'); if (c) c.scrollIntoView({ block: 'start', behavior: 'smooth' }); }
      else if (act === 'close') { UI.highlight = null; box.hidden = true; box.innerHTML = ''; setSiteHash(p, null); const cs = $('#crumb-site'); if (cs) cs.hidden = true; $$('#psites tr.row.sel').forEach((r) => r.classList.remove('sel')); drawMap(p, UNIDATA.get(p.acc)); }
    }));
    if (scroll) box.scrollIntoView({ block: 'nearest', behavior: 'smooth' });
  }

  function flashSite(pos, scroll) {
    renderSiteCard(pos, scroll);
    $$('#psites tr.row.sel').forEach((r) => r.classList.remove('sel'));
    const tr = $('#psites tr.row[data-pos="' + pos + '"]');
    if (tr) tr.classList.add('sel');
    $$('.pmap .head').forEach((h) => {
      const sel = +h.dataset.pos === pos;
      h.setAttribute('stroke-width', sel ? 2.6 : 2);
      h.setAttribute('stroke', sel ? getComputedStyle(document.documentElement).getPropertyValue('--ink').trim() : getComputedStyle(document.documentElement).getPropertyValue('--surface').trim());
    });
    focusStructure(pos);
    if (!tr) return;
    if (scroll && !$('#site-card:not([hidden])')) tr.scrollIntoView({ block: 'center', behavior: 'smooth' });
    tr.classList.remove('flash'); void tr.offsetWidth; tr.classList.add('flash');
  }

  // Domains & regions: UniProt features with the methyl-R sites they hold, plus the pipeline's InterPro classes.
  function renderDomainCard(p, u) {
    const box = $('#card-dom');
    if (!box) return;
    const pos = p.sites.map((j) => S[j].pos);
    const inside = (f) => pos.filter((x) => x >= f.start && x <= f.end).length;
    const near = (f) => pos.filter((x) => (x < f.start && f.start - x <= 20) || (x > f.end && x - f.end <= 20)).length;
    const cls = {};
    for (const j of p.sites) { const c = M.domain_classes[S[j].dcls]; cls[c] = (cls[c] || 0) + 1; }
    const ctx = {};
    for (const j of p.sites) { const c = DOMAIN_LABEL[M.domain_contexts[S[j].dctx]]; ctx[c] = (ctx[c] || 0) + 1; }
    let table = '<p class="note">' + (uniprotCache.has(p.acc) && u === null ? 'UniProt features could not be loaded.' : 'Loading UniProt features…') + '</p>';
    if (u) {
      const kinds = { Domain: 'Domain', Repeat: 'Repeat', 'Zinc finger': 'Zinc finger', Region: 'Region', Motif: 'Motif', 'Compositional bias': 'Composition' };
      const feats = u.feats.filter((f) => kinds[f.type]).map((f) => ({ name: f.desc || f.type, kind: kinds[f.type], start: f.start, end: f.end }))
        .concat(u.rg.map((b) => ({ name: b.label + ' (' + b.n + ' RG)', kind: 'RG/RGG', start: b.start, end: b.end, rg: true })))
        .sort((a, b) => a.start - b.start || b.end - a.end);
      table = feats.length ? '<div class="tbl-wrap"><table class="tbl compact"><thead><tr><th>Feature</th><th>Type</th><th class="r">Range</th><th class="r" title="Methylarginine sites inside the feature">Sites in</th><th class="r hide-xs" title="Sites within 20 aa outside the feature">±20 aa</th></tr></thead><tbody>' +
        feats.map((f) => { const n = inside(f); return '<tr' + (n ? ' class="has"' : '') + '><td>' + esc(f.name) + '</td><td><span class="tag kind-' + f.kind.replace(/[^A-Za-z]/g, '') + '">' + esc(f.kind) + '</span></td><td class="r num">' + f.start + '–' + f.end + '</td><td class="r num">' + (n || '<span class="muted">0</span>') + '</td><td class="r num hide-xs">' + (near(f) || '<span class="muted">0</span>') + '</td></tr>'; }).join('') +
        '</tbody></table></div>' : '<p class="muted">UniProt annotates no domains or regions for this entry.</p>';
    }
    box.innerHTML = '<h3>Domains &amp; regions</h3><p class="sub">UniProt features (live) with the methylarginines they contain. RG/RGG-rich segments are found in the sequence: ' + esc(RG_HELP.charAt(0).toLowerCase() + RG_HELP.slice(1)) + '.</p>' + table +
      '<p class="flabel" style="margin:1rem 0 .4rem">Pipeline context (InterPro, MobiDB)</p><div class="tags">' +
      Object.entries(ctx).sort((a, b) => b[1] - a[1]).map(([k, n]) => '<span class="tag">' + esc(k) + ' · ' + n + '</span>').join('') + '</div>' +
      '<div class="tags" style="margin-top:.35rem">' + Object.entries(cls).filter(([k]) => k && k !== 'No domain annotation').sort((a, b) => b[1] - a[1]).map(([k, n]) => '<button class="chip" type="button" data-cls="' + esc(k) + '" title="Sites whose nearest InterPro domain is in this class; click to filter all sites by it">' + esc(k) + ' · ' + n + '</button>').join('') + '</div>';
    $$('[data-cls]', box).forEach((b) => b.addEventListener('click', () => setTerm('domain', termSlug('domain', M.domain_classes.indexOf(b.dataset.cls)), 'sites')));
  }

  // PubMed summaries (first author, year, journal, title) for reference lists, from NCBI E-utilities.
  const PUBMED = new Map();
  function loadPubmed(ids) {
    const need = Array.from(new Set(ids.map(String))).filter((id) => /^\d+$/.test(id) && !PUBMED.has(id));
    if (!need.length) return Promise.resolve(false);
    need.forEach((id) => PUBMED.set(id, null));
    const batches = [];
    for (let k = 0; k < need.length; k += 200) batches.push(need.slice(k, k + 200));
    return Promise.all(batches.map((b) => fetchJson('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&tool=methylarginine_explorer&id=' + b.join(','))
      .then((d) => {
        const r = d.result || {};
        (r.uids || []).forEach((u) => { const x = r[u] || {}; PUBMED.set(u, { au: x.sortfirstauthor || '', year: (x.pubdate || '').slice(0, 4), journal: x.source || '', title: (x.title || '').replace(/\.$/, '') }); });
      }).catch(() => null))).then(() => true);
  }
  const pubmedUrl = (pm) => 'https://pubmed.ncbi.nlm.nih.gov/' + encodeURIComponent(pm) + '/';
  const pubmedAll = (ids) => 'https://pubmed.ncbi.nlm.nih.gov/?term=' + encodeURIComponent(ids.map((x) => x + '[pmid]').join(' OR '));
  const refShort = (pm) => { const x = PUBMED.get(String(pm)); return x && x.au ? x.au.split(' ')[0] + ' ' + x.year : 'PMID ' + pm; };
  const refFull = (pm) => {
    const x = PUBMED.get(String(pm));
    return x ? esc(x.au) + (x.au ? ' et al. ' : '') + (x.year ? '(' + esc(x.year) + ') ' : '') + esc(x.title) + (x.journal ? '. <i>' + esc(x.journal) + '</i>.' : '.') : 'PMID ' + esc(pm);
  };
  // Up to three short links, then a "+N more" expander with the full list and a PubMed search for all of them.
  function refsBlock(pmids, key) {
    if (!pmids || !pmids.length) return '';
    const link = (pm) => { const x = PUBMED.get(String(pm)); return '<a href="' + pubmedUrl(pm) + '" target="_blank" rel="noopener"' + (x ? ' title="' + esc(x.title) + '"' : '') + '>' + esc(refShort(pm)) + '</a>'; };
    const list = '<ol class="ref-list">' + pmids.map((pm) => '<li>' + refFull(pm) + ' <a href="' + pubmedUrl(pm) + '" target="_blank" rel="noopener">PMID ' + esc(pm) + ' ↗</a></li>').join('') + '</ol>';
    return '<div class="refs small">' + pmids.slice(0, 3).map(link).join(', ') +
      (pmids.length > 1 ? ' · <a href="' + pubmedAll(pmids) + '" target="_blank" rel="noopener">All ' + pmids.length + ' in PubMed ↗</a>' : '') +
      '<details class="refs-more" data-refkey="' + esc(key) + '"' + (UI.refsOpen.has(key) ? ' open' : '') + '><summary>' + (pmids.length > 3 ? '+' + (pmids.length - 3) + ' more · full list' : 'Full reference' + (pmids.length === 1 ? '' : 's')) + '</summary>' + list + '</details></div>';
  }
  function bindRefs(box) {
    $$('details[data-refkey]', box).forEach((d) => d.addEventListener('toggle', () => { if (d.open) UI.refsOpen.add(d.dataset.refkey); else UI.refsOpen.delete(d.dataset.refkey); }));
  }

  // One switch for every context card: what is known about this protein, or how its classes, condensates,
  // complexes and diseases behave across all methylarginine proteins.
  function ctxToggleHtml(p) {
    return '<div class="ctx-bar"><span class="flabel">Context cards</span><div class="seg" role="group" aria-label="Context cards show">' +
      '<button type="button" data-ctx="protein" aria-pressed="' + (UI.ctx !== 'dataset') + '">' + esc(p.label) + ' evidence</button>' +
      '<button type="button" data-ctx="dataset" aria-pressed="' + (UI.ctx === 'dataset') + '">Dataset context</button></div>' +
      '<span class="muted small">' + (UI.ctx === 'dataset' ? 'Condensates, ontology, disease and complexes show statistics across all methylarginine proteins.' : 'Condensates, ontology, disease and complexes show what is recorded for ' + esc(p.label) + ' itself.') + '</span></div>';
  }
  function bindCtx(p) {
    $$('.ctx-bar [data-ctx]').forEach((b) => b.addEventListener('click', () => setCtx(p, b.dataset.ctx)));
  }
  function setCtx(p, v) {
    UI.ctx = v === 'dataset' ? 'dataset' : 'protein';
    const bar = $('.ctx-bar');
    if (bar) { bar.outerHTML = ctxToggleHtml(p); bindCtx(p); }
    renderCondCard(p);
    renderFuncCard(p);
    loadEvidence().then((ev) => { if (UI.detail === p.i) renderDiseaseCard(p, ev); });
    if (UI.ixTab === 'complexes') renderIxCard(p);
  }

  const ROLE_LABEL = [[1, 'member'], [2, 'driver'], [4, 'client'], [8, 'regulator'], [16, 'other role']];
  const EVID_LABEL = [[1, 'in vivo'], [2, 'in cellulo'], [4, 'in vitro'], [8, 'FRAP'], [16, 'colocalisation'], [32, 'mass spectrometry'], [64, 'emergent properties'], [128, 'other']];
  const stars = (n) => '<span class="stars" role="img" aria-label="' + n + ' of 5">' + [1, 2, 3, 4, 5].map((k) => '<i class="' + (k <= n ? 'on' : '') + '">★</i>').join('') + '</span>';
  const cdMemberships = (p) => (ixReady() ? IX.cd.get(p.i) || [] : null);
  const condConfidence = (p, ci) => { const m = cdMemberships(p); if (!m) return null; const r = m.find((x) => x.ci === ci); return r ? r.score : 0; };

  // Condensates: this protein's own CD-CODE membership first; condensate-wide statistics in a separate, labelled view.
  function renderCondCard(p) {
    const box = $('#card-cond');
    if (!box) return;
    const view = UI.ctx === 'dataset' ? 'context' : 'protein';
    const mem = cdMemberships(p);
    const tabs = '<div class="seg cond-tabs" role="tablist"><button type="button" data-cv="protein" aria-pressed="' + (view === 'protein') + '">' + esc(p.label) + '</button><button type="button" data-cv="context" aria-pressed="' + (view === 'context') + '">Condensate context</button></div>';
    let body;
    if (view === 'context') {
      const ids = new Set(p.cond);
      if (mem) mem.forEach((m) => ids.add(m.ci));
      const rows = Array.from(ids).map((k) => ({ k, c: M.condensates[k] })).sort((x, y) => (y.c[5] == null ? -99 : y.c[5]) - (x.c[5] == null ? -99 : x.c[5]) || x.c[0].localeCompare(y.c[0]));
      const max = Math.max(1, ...rows.map((r) => Math.abs(r.c[5] || 0)));
      body = '<p class="sub">About the condensates, not about ' + esc(p.label) + ': how many members of each condensate are arginine-methylated in this dataset (methylated/all members) and the enrichment over the proteome (log2 OR). Click a name to list its methylated members.</p>' +
        (rows.length ? '<ul class="cond-list">' + rows.map(({ k, c }) => '<li><div class="cond-name"><button class="linkish strong" type="button" data-k="' + k + '" title="List all methylated members">' + esc(c[0]) + '</button></div>' +
          '<div class="cond-stat"><span class="num muted">' + c[4] + '/' + c[3] + '</span>' +
          (c[5] != null ? '<span class="mini-bar" title="log2 OR ' + fmtL2(c[5]) + ', q ' + fmtP(c[6]) + '"><i class="' + (c[6] != null && c[6] <= PADJ ? '' : 'ns') + (c[5] < 0 ? ' neg' : '') + '" style="width:' + (100 * Math.abs(c[5]) / max).toFixed(0) + '%"></i></span><span class="num">' + fmtL2(c[5]) + '</span>' : '<span class="muted small">not tested (&lt;5 members)</span>') +
          '</div></li>').join('') + '</ul>' : '<p class="muted">No CD-CODE condensate membership recorded.</p>');
    } else if (!mem) {
      body = '<p class="note">Loading CD-CODE membership evidence…</p>';
    } else {
      const byCi = new Map(mem.map((m) => [m.ci, m]));
      const ids = Array.from(new Set(p.cond.concat(mem.map((m) => m.ci))));
      const rows = ids.map((ci) => ({ ci, c: M.condensates[ci], m: byCi.get(ci) || null }))
        .sort((x, y) => ((y.m ? y.m.score : -1) - (x.m ? x.m.score : -1)) || ((y.m ? y.m.nPub : 0) - (x.m ? x.m.nPub : 0)) || x.c[0].localeCompare(y.c[0]));
      body = '<p class="sub">' + esc(p.label) + '\'s membership in each condensate as curated in CD-CODE ' + esc(IX.meta.cdcode_version || '') + '. Stars are CD-CODE\'s confidence that this protein belongs to the condensate (1–5); role, evidence types and publications are for this protein.</p>' +
        (rows.length ? '<ul class="cond-list mine">' + rows.map(({ ci, c, m }) => {
          const roles = m ? ROLE_LABEL.filter(([b2]) => m.roles & b2).map(([, l]) => l) : [];
          const ev = m ? EVID_LABEL.filter(([b2]) => m.evid & b2).map(([, l]) => l) : [];
          return '<li><div class="cond-name"><b>' + esc(c[0]) + '</b>' + (c[2] !== 'biomolecular' ? ' <span class="tag">' + esc(c[2]) + '</span>' : '') +
            '<a class="ext" href="https://cd-code.org/condensates/' + encodeURIComponent(c[1]) + '" target="_blank" rel="noopener" title="Open in CD-CODE">CD-CODE ↗</a></div>' +
            (m ? '<div class="cond-conf">' + stars(m.score) + '<span class="muted small">' + esc(roles.join(', ') || 'role not given') + '</span></div>' +
              '<div class="cond-ev small">' + (ev.length ? esc(ev.join(' · ')) : '<span class="muted">no evidence type recorded</span>') + ' · ' +
              (m.nPub ? m.nPub + ' publication' + (m.nPub === 1 ? '' : 's') : '<span class="muted">no publication recorded</span>') + '</div>' + refsBlock(m.pmids, 'cd:' + p.i + ':' + ci) :
              '<div class="cond-ev small muted">In the pipeline\'s membership table but not in the current CD-CODE release.</div>') + '</li>';
        }).join('') + '</ul>' : '<p class="muted">No CD-CODE condensate membership recorded.</p>') +
        '<p class="note">Low-star memberships usually rest on a single proteomics list. Condensate-wide methylation statistics are under “Condensate context”.</p>';
    }
    box.innerHTML = '<div class="card-head"><div><h3>Condensates</h3></div></div>' + tabs + body;
    $$('[data-cv]', box).forEach((b2) => b2.addEventListener('click', () => setCtx(p, b2.dataset.cv === 'context' ? 'dataset' : 'protein')));
    bindRefs(box);
    if (mem && view === 'protein') loadPubmed(mem.flatMap((m) => m.pmids)).then((fresh) => { if (fresh && UI.detail === p.i && UI.ctx !== 'dataset') renderCondCard(p); });
    $$('button[data-k]', box).forEach((b2) => b2.addEventListener('click', () => setTerm('condensate', termSlug('condensate', +b2.dataset.k), 'proteins')));
    const factSub = $('#cond-fact-sub');
    if (factSub && mem) {
      const best = mem.slice().sort((x, y) => y.score - x.score)[0];
      factSub.textContent = best ? 'highest confidence: ' + M.condensates[best.ci][0] + ' (' + best.score + '/5)' : (p.cond.length ? 'not in current CD-CODE' : 'none annotated');
    }
  }

  // Functional ontology tree: broad classes with this protein's subclasses nested under their parent class.
  function renderFuncCard(p) {
    const box = $('#card-func');
    const ov = M.overview;
    const broad = bits(p.func).map((b) => M.func_classes[b]);
    const subs = bits(p.sub).map((b) => M.func_subclasses[b]);
    // Dataset-level class enrichment stays in the tooltip; the card itself lists only this protein's classes.
    const stat = (name, table) => {
      if (UI.ctx !== 'dataset' || !name) return '';
      const r = (table || []).find((x) => x[0] === name);
      return r ? ' <span class="tag' + (r[5] != null && r[5] <= PADJ ? '' : ' src') + '" title="Class enrichment among methylarginine proteins: log2 OR ' + fmtL2(r[4]) + ', q ' + fmtP(r[5]) + '">' + fmtL2(r[4]) + '</span>' : '';
    };
    const tipOf = (name, table) => { const r = (table || []).find((x) => x[0] === name); return r ? ' title="' + esc('Click to list the methylated proteins in this class. Across the dataset: log2 OR ' + fmtL2(r[4]) + ', q ' + fmtP(r[5])) + '"' : ''; };
    const used = new Set();
    const tree = broad.map((b) => {
      const kids = subs.filter((s) => (M.subParent[s] || '') === b && s !== b);
      kids.forEach((k) => used.add(k));
      return '<li><button class="linkish strong" type="button" data-func="' + esc(b) + '"' + tipOf(b, ov.functional_classes) + '>' + esc(b) + '</button>' + stat(b, ov.functional_classes) +
        (kids.length ? '<ul>' + kids.map((k) => '<li><button class="linkish" type="button" data-sub="' + esc(k) + '"' + tipOf(k, ov.functional_subclasses) + '>' + esc(k) + '</button>' + stat(k, ov.functional_subclasses) + '</li>').join('') + '</ul>' : '') + '</li>';
    });
    const rest = subs.filter((s) => !used.has(s) && !broad.includes(s));
    box.innerHTML = '<h3>Functional ontology</h3><p class="sub">The lab functional classes (config/functional_ontology.json) assigned to ' + esc(p.label) + ', with subclasses nested under their usual parent. Click a class to list the other methylated proteins in it. ' + (UI.ctx === 'dataset' ? 'Tags: each class\'s enrichment among all methylarginine proteins (log2 OR; outlined = not significant).' : 'Switch to Dataset context to see how strongly each class is enriched across the dataset.') + '</p>' +
      (broad.length || subs.length ? '<ul class="tree">' + tree.join('') + (rest.length ? '<li><span class="muted">Other subclasses</span><ul>' + rest.map((k) => '<li><button class="linkish" type="button" data-sub="' + esc(k) + '"' + tipOf(k, ov.functional_subclasses) + '>' + esc(k) + '</button>' + stat(k, ov.functional_subclasses) + '</li>').join('') + '</ul></li>' : '') + '</ul>' : '<p class="muted">No functional class assigned.</p>');
    $$('[data-func]', box).forEach((b) => b.addEventListener('click', () => setTerm('class', termSlug('class', M.func_classes.indexOf(b.dataset.func)), 'proteins')));
    $$('[data-sub]', box).forEach((b) => b.addEventListener('click', () => setTerm('subclass', termSlug('subclass', M.func_subclasses.indexOf(b.dataset.sub)), 'proteins')));
  }

  // Open Targets Platform: this gene's disease associations and their overall scores (CC0).
  function loadOpenTargets(ensg) {
    return live('ot:' + ensg, () => fetch('https://api.platform.opentargets.org/api/v4/graphql', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ query: 'query($id:String!){target(ensemblId:$id){approvedSymbol associatedDiseases(page:{index:0,size:12}){count rows{score disease{id name} datatypeScores{id score}}}}}', variables: { id: ensg } }),
    }).then((r) => { if (!r.ok) throw new Error(r.status); return r.json(); })
      .then((d) => (d.data && d.data.target ? d.data.target.associatedDiseases : { count: 0, rows: [] })));
  }
  const OT_TYPES = { genetic_association: 'genetics', somatic_mutation: 'somatic', known_drug: 'drugs', affected_pathway: 'pathways', rna_expression: 'expression', literature: 'literature', animal_model: 'animal models' };

  // Disease: what is known about this gene first (UniProt curation, Open Targets scores); the dataset-level
  // DisGeNET enrichment terms that include it come last and are labelled as such.
  function renderDiseaseCard(p, ev) {
    const box = $('#card-dis');
    if (!box) return;
    const u = UNIDATA.get(p.acc);
    const uni = u === undefined ? '<p class="note">Loading UniProt…</p>' : !u ? '<p class="muted">UniProt unavailable.</p>' :
      u.diseases.length ? '<ul class="term-list uni-dis">' + u.diseases.map((d) => '<li><b>' + esc(d.name) + '</b>' + (d.acronym ? ' (' + esc(d.acronym) + ')' : '') +
        (d.mim ? ' <a href="https://omim.org/entry/' + encodeURIComponent(d.mim) + '" target="_blank" rel="noopener">MIM ' + esc(d.mim) + ' ↗</a>' : '') +
        (d.note ? '<br><span class="muted small">' + esc(d.note.length > 220 ? d.note.slice(0, 218) + '…' : d.note) + '</span>' : '') + '</li>').join('') + '</ul>' :
        '<p class="muted">No disease involvement curated in UniProt.</p>';
    let dg = '';
    if (ev) {
      const list = (ev.protein_disease_terms[p.i] || []).map((k) => ev.disease_terms[k]);
      const groups = new Map();
      for (const t of list) { const g = t[2] >= 0 ? M.disease_families[t[2]] : 'Other / unassigned'; if (!groups.has(g)) groups.set(g, []); groups.get(g).push(t); }
      const order = Array.from(groups.keys()).sort((x, y) => (x === 'Other / unassigned') - (y === 'Other / unassigned') || groups.get(y).length - groups.get(x).length);
      dg = order.length ? '<div class="dis-groups">' + order.map((g) => {
        const terms = groups.get(g).sort((x, y) => (x[4] || 1) - (y[4] || 1));
        const fi = M.disease_families.indexOf(g);
        return '<details><summary>' + (fi >= 0 ? '<button class="linkish strong" type="button" data-fam="' + fi + '">' + esc(g) + '</button>' : '<span class="strong">' + esc(g) + '</span>') +
          ' <span class="muted">' + terms.length + ' term' + (terms.length === 1 ? '' : 's') + '</span></summary><ul class="term-list">' +
          terms.map((t) => '<li title="' + esc('Term enrichment across all methylarginine proteins: p.adjust ' + fmtP(t[4])) + '">' + esc(t[1]) + ' <span class="muted">· ' + esc(t[0]) + (UI.ctx === 'dataset' ? ' · p.adj ' + fmtP(t[4]) : '') + '</span></li>').join('') + '</ul></details>';
      }).join('') + '</div>' : '<p class="muted">' + esc(p.label) + ' is in none of the enriched DisGeNET terms.</p>';
    } else dg = '<p class="note">Loading terms…</p>';
    box.innerHTML = '<h3>Disease associations</h3>' +
      (p.als ? '<p><span class="tag ip" title="Member of the full audited ALS / motor neuron / FTD family gene set">In the audited ALS / MN / FTD gene set</span></p>' : '') +
      '<p class="flabel">Curated in UniProt</p>' + uni +
      '<p class="flabel">Open Targets association scores</p><div id="ot-box">' + (u && u.ensg && u.ensg.length ? '<p class="note">Loading Open Targets…</p>' : '<p class="muted">' + (u === undefined ? 'Waiting for UniProt…' : 'No Ensembl gene for this entry.') + '</p>') + '</div>' +
      '<details class="dataset-ctx"' + (UI.ctx === 'dataset' ? ' open' : '') + '><summary class="flabel">Dataset context: DisGeNET terms enriched among all methylarginine proteins that include ' + esc(p.label) + '</summary>' +
      '<p class="note">These terms were tested for enrichment across the whole methylarginine set (clusterProfiler, p.adjust ≤ 0.05); the gene–term link is DisGeNET\'s, the statistics are about the set. Association, not causation.</p>' + dg + '</details>';
    $$('[data-fam]', box).forEach((b2) => b2.addEventListener('click', (e) => { e.preventDefault(); setTerm('disease', termSlug('disease', b2.dataset.fam === 'als' ? 'als' : +b2.dataset.fam), 'proteins'); }));
    if (u && u.ensg && u.ensg.length) {
      loadOpenTargets(u.ensg[0]).then((d) => {
        const el = $('#ot-box', box);
        if (!el || UI.detail !== p.i) return;
        if (!d) { el.innerHTML = '<p class="muted">Open Targets did not answer.</p>'; return; }
        if (!d.rows.length) { el.innerHTML = '<p class="muted">No Open Targets disease associations.</p>'; return; }
        el.innerHTML = '<ul class="ot-list">' + d.rows.map((r) => {
          const types = (r.datatypeScores || []).filter((t) => t.score > 0).sort((x, y) => y.score - x.score).slice(0, 3).map((t) => OT_TYPES[t.id] || t.id.replace(/_/g, ' '));
          return '<li><span class="ot-name"><a href="https://platform.opentargets.org/evidence/' + encodeURIComponent(u.ensg[0]) + '/' + encodeURIComponent(r.disease.id) + '" target="_blank" rel="noopener">' + esc(r.disease.name) + '</a>' +
            (types.length ? ' <span class="muted small">' + esc(types.join(', ')) + '</span>' : '') + '</span><span class="mini-bar" title="Open Targets overall association score ' + r.score.toFixed(2) + '"><i style="width:' + (100 * r.score).toFixed(0) + '%"></i></span><span class="num">' + r.score.toFixed(2) + '</span></li>';
        }).join('') + '</ul><p class="note">Top ' + d.rows.length + ' of ' + fmt(d.count) + ' diseases by overall association score (0–1) for this gene. <a href="https://platform.opentargets.org/target/' + encodeURIComponent(u.ensg[0]) + '/associations" target="_blank" rel="noopener">All in Open Targets ↗</a> · CC0</p>';
      });
    }
  }

  function renderClusterCard(p) {
    const box = $('#card-cluster');
    const d = p.dens;
    const sig = p.sig.slice().sort((a, b) => a[0] - b[0]);
    box.innerHTML = '<h3>Clustering</h3><p class="sub">Local methyl-R density from the multiscale rolling-window scan against a same-protein arginine null (1,000 permutations).</p>' +
      (d ? '<dl class="kv"><dt>Closest / median spacing</dt><dd>' + d[0] + ' aa / ' + d[1] + ' aa</dd><dt>Max sites within</dt><dd>' +
        [3, 5, 10, 20, 50].map((w, k) => w + ' aa: <b>' + d[2 + k] + '</b>').join(' · ') + '</dd>' +
        (p.best ? '<dt>Densest window</dt><dd>' + p.best[3] + ' sites in ' + p.best[0] + ' aa (' + rl(p.best[1], p) + '–' + rn(p.best[2], p) + '), log2 obs/null ' + fmtL2(p.best[4]) + ', q ' + fmtP(p.best[5]) + '</dd>' : '') +
        (p.arch ? '<dt>Architecture</dt><dd>' + esc(p.arch[0]) + (p.arch[1].length ? '<br><span class="muted">clusters ' + p.arch[1].map((c) => c[0] + '–' + c[1] + ' (' + c[2] + ')').join(', ') + '</span>' : '') + '</dd>' : '') +
        '</dl>' : '<p class="muted">Clustering statistics need sequence-confirmed sites.</p>') +
      (sig.length ? '<p style="margin:.9rem 0 .4rem"><b>Significant dense windows</b> (q ≤ 0.05)</p><div class="chips">' +
        sig.map((w, k) => '<button class="chip" type="button" data-k="' + k + '" title="log2 obs/null ' + fmtL2(w[4]) + ', q ' + fmtP(w[5]) + '">' + w[3] + ' sites / ' + w[0] + ' aa · ' + w[1] + '–' + w[2] + '</button>').join('') + '</div><p class="note">Click a window to zoom the figure.</p>' : '');
    $$('.chip', box).forEach((b) => b.addEventListener('click', () => {
      const w = sig[+b.dataset.k];
      const pad = Math.max(10, Math.round((w[2] - w[1]) * 0.6));
      UI.zoom = [Math.max(1, w[1] - pad), w[2] + pad];
      loadUniProt(p.acc).then((u) => drawMap(p, u));
      $('#map-box').scrollIntoView({ block: 'center', behavior: 'smooth' });
    }));
  }

  function renderPtmCard(p) {
    const box = $('#card-ptm');
    const counts = [0, 0, 0, 0];
    for (const [, g] of p.ptms) counts[g]++;
    const near = [];
    for (const [pos, g, res] of p.ptms) {
      let best = null;
      for (const j of p.sites) { const d = Math.abs(S[j].pos - pos); if (d <= 10 && (!best || d < best[1])) best = [S[j].pos, d]; }
      if (best) near.push([pos, g, res, best[0], best[1]]);
    }
    near.sort((a, b) => a[4] - b[4] || a[0] - b[0]);
    box.innerHTML = '<h3>Other PTMs</h3><p class="sub">Other modification types recorded for this protein in the pipeline\'s cross-PTM master table.</p>' +
      '<p class="note ptm-caveat"><b>Not an exhaustive list:</b> it holds only what the pipeline\'s cross-PTM master table contains (an iPTMnet export of HPRD, literature mining and other curated sources, plus UniProt), so modifications reported elsewhere, for example in PhosphoSitePlus or recent proteomics studies, may be missing.</p>' +
      (p.ptms.length ? '<dl class="kv">' + counts.map((c, g) => (c ? '<dt>' + PTM_LABEL[g] + '</dt><dd>' + c + '</dd>' : '')).join('') + '</dl>' +
        (near.length ? '<p style="margin:.9rem 0 .4rem"><b>Within 10 aa of a methyl-R</b></p><div class="tags">' + near.slice(0, 40).map((n) => '<span class="tag" title="' + n[4] + ' aa from ' + rl(n[3], p) + '">' + PTM_PREFIX[n[1]] + esc(n[2]) + n[0] + ' · ' + rl(n[3], p) + '</span>').join('') + '</div>' + (near.length > 40 ? '<p class="note">' + (near.length - 40) + ' more not shown.</p>' : '') : '<p class="note">None within 10 aa of a methylarginine.</p>')
        : '<p class="muted">None in the master table for this accession, which does not mean the protein has no other modifications.</p>');
  }

  function renderSequence(p, u) {
    const box = $('#seq');
    const sub = $('#seq-sub');
    if (!box) return;
    if (!u || !u.seq) { sub.textContent = 'The UniProt sequence could not be loaded. Open the UniProt entry for the full record.'; return; }
    const seq = u.seq;
    const siteAt = new Map(p.sites.map((j) => [S[j].pos, S[j]]));
    let mismatch = 0;
    siteAt.forEach((s, pos) => { if (seq[pos - 1] !== 'R') mismatch++; });
    const inRg = new Uint8Array(seq.length + 2);
    for (const b of u.rg) for (let k = b.start; k <= b.end; k++) inRg[k] = 1;
    sub.innerHTML = 'Current UniProt sequence, ' + fmt(seq.length) + ' aa. Methylated arginines are filled by state, other arginines are bold, glycines are muted, and RG/RGG-rich segments are underlined. Click a filled R to find its row.' +
      (mismatch ? ' <span class="warn">' + mismatch + ' reported site' + (mismatch === 1 ? ' does' : 's do') + ' not fall on R in this sequence version; positions follow the pipeline\'s staged canonical FASTA.</span>' : '');
    let lines = '';
    for (let start = 0; start < seq.length; start += 60) {
      let line = '';
      for (let k = start; k < Math.min(seq.length, start + 60); k++) {
        if (k > start && (k - start) % 10 === 0) line += ' ';
        const ch = seq[k];
        const s = siteAt.get(k + 1);
        const rg = inRg[k + 1] ? ' rg' : '';
        if (s && ch === 'R') line += '<span class="m st' + s.st + rg + '" data-pos="' + (k + 1) + '" title="' + rl(k + 1, p) + ' · ' + STATE_LABEL[s.st] + '">R</span>';
        else if (ch === 'R') line += '<span class="R' + rg + '">R</span>';
        else if (ch === 'G') line += '<span class="G' + rg + '">G</span>';
        else line += rg ? '<span class="rg">' + ch + '</span>' : ch;
      }
      lines += '<span class="ln">' + (start + 1) + '</span><span class="line">' + line + '</span>';
    }
    box.innerHTML = '<div class="seq">' + lines + '</div>';
    $$('.seq span.m', box).forEach((m) => m.addEventListener('click', () => flashSite(+m.dataset.pos, true)));
    const btn = $('#copy-fasta');
    btn.disabled = false;
    btn.onclick = (e) => {
      const methyl = new Set(p.sites.map((j) => S[j].pos));
      const body = seq.split('').map((c, k) => (methyl.has(k + 1) && c === 'R' ? 'r' : c)).join('').replace(/(.{60})/g, '$1\n');
      copyText('>' + p.acc + ' ' + p.label + ' methylarginines in lowercase (' + p.sites.length + ' sites)\n' + body.trim() + '\n', e.currentTarget);
    };
    const exportSeq = (kind) => {
      exportFigure(kind, seqSvg(p, u), (p.gene || p.label || p.acc) + '_' + p.acc + '_sequence');
    };
    for (const kind of ['svg', 'png']) {
      const b = $('#seq-' + kind);
      b.disabled = false;
      b.onclick = () => exportSeq(kind);
    }
  }

  // The sequence card as a stand-alone figure: title (gene, protein name), UniProt accession and length, legend,
  // 60 residues per line in blocks of 10. Every character gets an explicit x so columns line up in any font.
  function seqSvg(p, u) {
    const pal = EXPORT_PAL, seq = u.seq;
    const stColor = [pal.none, pal.cat[0], pal.cat[1], pal.cat[2]];
    const siteAt = new Map(p.sites.map((j) => [S[j].pos, S[j]]));
    const inRg = new Uint8Array(seq.length + 2);
    for (const b of u.rg || []) for (let k = b.start; k <= b.end; k++) inRg[k] = 1;
    const monoAttr = ' font-family=\'"DejaVu Sans Mono", Menlo, Consolas, "Courier New", monospace\'';
    const cw = 8.4, lh = 19, fs = 13, numW = 46, x0 = 16 + numW + 12;
    const W = Math.ceil(x0 + 65 * cw + 16);
    const colX = (c) => x0 + (c + Math.floor(c / 10)) * cw;
    let g = '';
    const states = Array.from(new Set(p.sites.map((j) => S[j].st))).sort();
    g += '<text x="16" y="28" font-size="18" font-weight="700" fill="' + pal.ink + '">' + esc(p.label) + (p.name ? ' <tspan font-size="13" font-weight="400" fill="' + pal.muted + '">' + esc(p.name) + '</tspan>' : '') + '</text>';
    g += '<text x="16" y="47" font-size="11.5" fill="' + pal.muted + '">UniProt ' + esc(p.acc) + ' · ' + fmt(seq.length) + ' aa · ' + p.sites.length + ' methylarginine site' + (p.sites.length === 1 ? '' : 's') + ' · Homo sapiens</text>';
    // Legend
    let lx = 16, ly = 70;
    // Legend items wrap to a new row rather than run past the right edge.
    const item = (draw, label) => {
      const w = 30 + label.length * 6.6;
      if (lx > 16 && lx + w > W - 16) { lx = 16; ly += 22; }
      g += draw(lx, ly) + '<text x="' + (lx + 18) + '" y="' + (ly + 4) + '" font-size="11" fill="' + pal.ink + '">' + esc(label) + '</text>';
      lx += w;
    };
    states.sort((a, b) => (a === 0) - (b === 0) || a - b).forEach((st) => item((x, ly) => '<rect x="' + x + '" y="' + (ly - 7) + '" width="12" height="14" rx="2" fill="' + stColor[st] + '"/><text x="' + (x + 6) + '" y="' + (ly + 4) + '" font-size="10" font-weight="700" text-anchor="middle" fill="#FFFFFF"' + monoAttr + '>R</text>', 'Methylated R, ' + STATE_LABEL[st]));
    item((x, ly) => '<text x="' + (x + 6) + '" y="' + (ly + 4) + '" font-size="12" font-weight="700" text-anchor="middle" fill="#155E75"' + monoAttr + '>R</text>', 'Other arginine');
    if (u.rg && u.rg.length) item((x, ly) => '<line x1="' + x + '" x2="' + (x + 12) + '" y1="' + (ly + 7) + '" y2="' + (ly + 7) + '" stroke="' + pal.site + '" stroke-width="2"/>', 'RG/RGG-rich segment');
    // Sequence
    let y = ly + 34;
    for (let start = 0; start < seq.length; start += 60) {
      g += '<text x="' + (16 + numW) + '" y="' + y + '" font-size="11" text-anchor="end" fill="' + pal.muted + '">' + (start + 1) + '</text>';
      const runs = [];
      for (let k = start; k < Math.min(seq.length, start + 60); k++) {
        const ch = seq[k], c = k - start, s = siteAt.get(k + 1);
        const cls = s && ch === 'R' ? 'm' : ch === 'R' ? 'R' : ch === 'G' ? 'G' : 'n';
        if (cls === 'm') g += '<rect x="' + r1(colX(c) - 0.6) + '" y="' + (y - fs + 1) + '" width="' + r1(cw + 1.2) + '" height="' + (fs + 4) + '" rx="2" fill="' + stColor[s.st] + '"/>';
        if (inRg[k + 1]) g += '<line x1="' + r1(colX(c)) + '" x2="' + r1(colX(c) + cw) + '" y1="' + (y + 4.5) + '" y2="' + (y + 4.5) + '" stroke="' + pal.site + '" stroke-width="2"/>';
        const last = runs[runs.length - 1];
        if (last && last.cls === cls) { last.xs.push(colX(c) + cw / 2); last.t += ch; } else runs.push({ cls, xs: [colX(c) + cw / 2], t: ch });
      }
      for (const r of runs) {
        const style = r.cls === 'm' ? ' font-weight="700" fill="#FFFFFF"' : r.cls === 'R' ? ' font-weight="700" fill="#155E75"' : r.cls === 'G' ? ' fill="' + pal.muted + '"' : ' fill="' + pal.ink + '"';
        g += '<text x="' + r.xs.map(r1).join(' ') + '" y="' + y + '" font-size="' + fs + '" text-anchor="middle"' + monoAttr + style + '>' + esc(r.t) + '</text>';
      }
      y += lh;
    }
    y += 10;
    g += '<text x="16" y="' + y + '" font-size="9" fill="' + pal.muted + '">Sequence: UniProt ' + esc(p.acc) + ' (CC BY 4.0) · sites: methylarginine union, PTM_Pipeline</text>';
    const H = Math.ceil(y + 12);
    const svg = '<svg xmlns="http://www.w3.org/2000/svg" width="' + W + '" height="' + H + '" viewBox="0 0 ' + W + ' ' + H + '" font-family=\'' + pal.font.replace(/'/g, '"') + '\'>' +
      '<rect width="' + W + '" height="' + H + '" fill="' + pal.bg + '"/>' + g + '</svg>';
    return { svg, W, H };
  }

  // ------------------------------------------------------------- figure
  // ----- contacts at a glance: one protein as a linear cartoon, its methylarginines coloured by interaction status
  // and, underneath on the same residue scale, the contact evidence: RNA and DNA crosslinks in cells, PDB interfaces
  // (PDBe-KB, live) and predicted-complex interfaces (assessed at methylarginines only). A dashed guide drops from
  // every methylarginine that has contact evidence at the residue itself, so coincidences read at a glance.
  const IXS_COLOR = { exp: '#0E7490', pred: '#C2410C' };
  function contactLanes(p, ifd) {
    const xp = ixReady() ? IX.xlProt.get(p.i) : null;
    const sites = p.sites.map((j) => S[j]);
    return [
      { key: 'rna', short: 'RNA xl', label: 'RNA crosslink', sub: 'in living cells', color: IFACE_COLOR[4],
        ticks: xp ? xp.rnaPos : [], boxes: xp ? xp.rbd.map(([a, b]) => [a, b, 0.25, 'RBDmap RNA-bound fragment']) : [],
        hit: (s) => !!(s.xl && s.xl.r), weak: (s) => !!(s.xl && (s.xl.f || s.xl.n)) },
      { key: 'dna', short: 'DNA xl', label: 'DNA crosslink', sub: 'in living cells', color: IFACE_COLOR[8],
        ticks: xp ? xp.dnaPos : [], boxes: xp ? xp.dnaPeps.map(([a, b]) => [a, b, 0.35, 'DNA-crosslinked peptide']).concat(xp.dnaCand.map((c) => [c, c, 0.7, 'DNA crosslink candidate residue'])) : [],
        hit: (s) => !!(s.xl && s.xl.d & 9), weak: (s) => !!(s.xl && s.xl.d & 2) },
      { key: 'pdbna', short: 'PDB RNA/DNA', label: 'PDB: RNA/DNA', sub: 'interface residues', color: IFACE_COLOR[4],
        ranges: ifd ? ifd.rna.map((r) => r.concat(IFACE_COLOR[4])).concat(ifd.dna.map((r) => r.concat(IFACE_COLOR[8]))) : null,
        hit: (s) => !!(s.iface & 12) },
      { key: 'pdbprot', short: 'PDB protein', label: 'PDB: protein', sub: 'interface residues', color: IFACE_COLOR[1],
        ranges: ifd ? ifd.protein.map((r) => r.concat(IFACE_COLOR[1])).concat(ifd.self.map((r) => r.concat(IFACE_COLOR[2]))) : null,
        hit: (s) => !!(s.iface & 3) },
      { key: 'pred', short: 'Predicted', label: 'Predicted interface', sub: 'checked at methyl-R only', color: IXS_COLOR.pred,
        ticks: sites.filter((s) => (s.ev || 0) & EV_PRED).map((s) => s.pos), hit: (s) => !!((s.ev || 0) & EV_PRED) },
    ];
  }
  function contactCartoonSvg(p, u, ifd, o) {
    const pal = o.pal || EXPORT_PAL, exp = !!o.forExport, panel = !!o.panel;
    const W = o.W || 1000, m = panel ? 4 : 16, labelW = panel ? 64 : 150, x0 = m + labelW, x1 = W - (panel ? 6 : 96);
    const len = (u && u.len) || p.len || Math.ceil(p.maxPos * 1.05) || 100;
    const X = (pos) => x0 + ((pos - 0.5) / len) * (x1 - x0);
    const sites = p.sites.map((j) => S[j]);
    const lanes = contactLanes(p, ifd);
    const stColor = (s) => { const st = ixStatus(s); return st === 'exp' ? IXS_COLOR.exp : st === 'pred' ? IXS_COLOR.pred : pal.none; };
    let g = '', y = 8;
    if (exp && !panel) {
      g += '<text x="16" y="26" font-size="17" font-weight="700" fill="' + pal.ink + '">' + esc(o.title || p.label + ': contact evidence along the protein') + '</text>' +
        '<text x="16" y="44" font-size="11" fill="' + pal.muted + '">' + esc(o.sub || p.acc + ' · ' + fmt(len) + ' aa · ' + sites.length + ' methylarginines') + '</text>';
      y = 58;
    }
    // Protein row: disorder shading, backbone, domains, then lollipops stacked to avoid overlap.
    const headTop = y + 17, base = y + 76;
    const feats = (u && u.feats) || [];
    g += '<text x="' + m + '" y="' + (base + 4) + '" font-size="' + r1(Math.min(11.5, (labelW - 4) / (p.label.length * 0.62))) + '" font-weight="700" fill="' + pal.ink + '">' + esc(p.label) + '</text>' +
      '<text x="' + m + '" y="' + (base + 16) + '" font-size="' + (panel ? 8.5 : 9.5) + '" fill="' + pal.muted + '">' + (panel ? '' : esc(p.acc) + ' · ') + fmt(len) + ' aa</text>';
    for (const f of feats.filter((f) => f.type === 'Region' && /disordered/i.test(f.desc))) {
      g += '<rect x="' + r1(X(f.start)) + '" y="' + (base - 10) + '" width="' + r1(Math.max(1, X(f.end + 1) - X(f.start))) + '" height="20" fill="' + pal.disFill + '" opacity="0.6"><title>Disordered ' + f.start + '–' + f.end + ' (UniProt)</title></rect>';
    }
    g += '<rect x="' + r1(X(1)) + '" y="' + (base - 3) + '" width="' + r1(X(len + 1) - X(1)) + '" height="6" rx="3" fill="' + pal.backboneFill + '" stroke="' + pal.backboneEdge + '"/>';
    for (const f of feats.filter((f) => ['Domain', 'Repeat', 'Zinc finger'].includes(f.type))) {
      const a = X(f.start), wd = Math.max(2, X(f.end + 1) - a), rep = f.type === 'Repeat';
      g += '<rect x="' + r1(a) + '" y="' + (base - (rep ? 5 : 8)) + '" width="' + r1(wd) + '" height="' + (rep ? 10 : 16) + '" rx="3" fill="' + pal.domainFill + '"' + (rep ? ' opacity="0.6"' : '') + '><title>' + esc(f.desc + ' ' + f.start + '–' + f.end) + '</title></rect>';
      const lab = rep ? '' : (f.desc || '').replace(/ domain$/i, '');
      if (lab && lab.length * 5.4 < wd - 4) g += '<text x="' + r1(a + wd / 2) + '" y="' + (base + 3.5) + '" font-size="9" text-anchor="middle" fill="' + pal.domainInk + '">' + esc(lab) + '</text>';
    }
    // Lysine methylation (PTM master table) as diamonds on the backbone, for methyl-R/K-methyl crosstalk.
    const kme = p.ptms.filter((x) => x[1] === 3);
    for (const [pos] of kme) g += '<rect x="' + r1(X(pos) - 3) + '" y="' + (base - 3) + '" width="6" height="6" transform="rotate(45 ' + r1(X(pos)) + ' ' + base + ')" fill="' + KME_COLOR + '" stroke="' + pal.ring + '" stroke-width="0.6"><title>' + esc(p.label + ' K' + pos + ' methylation (PTM master table)') + '</title></rect>';
    // Lanes start under the backbone; guides need the lane block's extent first.
    const laneTop = base + 22, laneH = 22;
    const laneY = (k) => laneTop + k * laneH + laneH / 2;
    const bottom = laneTop + lanes.length * laneH;
    // Guides from methylarginines with evidence, drawn first so marks sit on top.
    const levels = [];
    const order = sites.slice().sort((a, b) => (ixStatus(a) === 'none') - (ixStatus(b) === 'none') || a.pos - b.pos);
    const heads = [];
    for (const s of sites.slice().sort((a, b) => a.pos - b.pos)) {
      const hx = X(s.pos);
      let lv = 0;
      while (levels[lv] != null && hx - levels[lv] < 7 && lv < 4) lv++;
      levels[lv] = hx;
      heads.push({ s, hx, hy: base - 16 - lv * 9 });
    }
    for (const h of heads) {
      const st = ixStatus(h.s);
      if (st && st !== 'none') g += '<line x1="' + r1(h.hx) + '" x2="' + r1(h.hx) + '" y1="' + (base + 8) + '" y2="' + bottom + '" stroke="' + IXS_COLOR[st] + '" stroke-width="0.9" stroke-dasharray="2 2" opacity="0.75"/>';
    }
    // Lanes.
    lanes.forEach((ln, k) => {
      const cy = laneY(k);
      g += '<line x1="' + x0 + '" x2="' + x1 + '" y1="' + cy + '" y2="' + cy + '" stroke="' + pal.grid + '"/>' +
        (panel ? '<text x="' + (x0 - 4) + '" y="' + (cy + 3) + '" font-size="8.5" text-anchor="end" fill="' + pal.ink + '">' + esc(ln.short) + '</text>'
          : '<text x="' + (x0 - 8) + '" y="' + (cy - 1) + '" font-size="10.5" text-anchor="end" fill="' + pal.ink + '">' + esc(ln.label) + '</text>' +
            '<text x="' + (x0 - 8) + '" y="' + (cy + 9) + '" font-size="8.5" text-anchor="end" fill="' + pal.muted + '">' + esc(ln.sub) + '</text>');
      for (const [a, b, op, what] of ln.boxes || []) g += '<rect x="' + r1(X(a)) + '" y="' + (cy - 4) + '" width="' + r1(Math.max(1.6, X(b + 1) - X(a))) + '" height="8" rx="1" fill="' + ln.color + '" fill-opacity="' + op + '"><title>' + esc(what + ' ' + (a === b ? a : a + '–' + b)) + '</title></rect>';
      for (const [a, b, name, color] of ln.ranges || []) g += '<rect x="' + r1(X(a)) + '" y="' + (cy - 5) + '" width="' + r1(Math.max(1.6, X(b + 1) - X(a))) + '" height="10" rx="1" fill="' + color + '" fill-opacity="0.85"><title>' + esc(ln.label + ' ' + (a === b ? a : a + '–' + b) + (name ? ' · ' + name : '')) + '</title></rect>';
      for (const t of ln.ticks || []) g += '<rect x="' + r1(X(t) - 0.9) + '" y="' + (cy - 6) + '" width="' + r1(Math.max(1.8, X(t + 1) - X(t))) + '" height="12" fill="' + ln.color + '"><title>' + esc(ln.label + ' at residue ' + t) + '</title></rect>';
      if (ln.ranges === null || (ln.key.startsWith('pdb') && ifd && ifd.failed)) g += '<text x="' + ((x0 + x1) / 2) + '" y="' + (cy + 3) + '" font-size="9.5" text-anchor="middle" fill="' + pal.muted + '">' + (ln.ranges === null ? 'loading PDBe-KB…' : 'PDBe-KB did not answer; site-level PDB evidence only') + '</text>';
      // Rings where a methylarginine itself carries this evidence; a count on the right.
      const hits = sites.filter(ln.hit), weak = ln.weak ? sites.filter((s) => !ln.hit(s) && ln.weak(s)) : [];
      for (const s of weak) g += '<circle cx="' + r1(X(s.pos)) + '" cy="' + cy + '" r="3.6" fill="none" stroke="' + ln.color + '" stroke-width="1" stroke-dasharray="1.5 1.2"/>';
      for (const s of hits) g += '<circle cx="' + r1(X(s.pos)) + '" cy="' + cy + '" r="4.2" fill="' + pal.bg + '" stroke="' + ln.color + '" stroke-width="1.8"><title>' + esc(p.label + ' ' + rl(s.pos, P[s.p]) + ': ' + ln.label + ' at this methylarginine') + '</title></circle>';
      if (!panel) g += '<text x="' + (x1 + 10) + '" y="' + (cy + 3.5) + '" font-size="10" fill="' + (hits.length ? pal.ink : pal.muted) + '">' + hits.length + ' methyl-R</text>';
    });
    // Lollipops last: stems, heads coloured by status, labels for sites with evidence when there is room.
    for (const h of heads) {
      const st = ixStatus(h.s);
      g += '<line x1="' + r1(h.hx) + '" x2="' + r1(h.hx) + '" y1="' + (base - 8) + '" y2="' + h.hy + '" stroke="' + pal.stem + '" stroke-width="1"/>' +
        '<circle class="cc-site" data-pos="' + h.s.pos + '" cx="' + r1(h.hx) + '" cy="' + h.hy + '" r="' + (st && st !== 'none' ? 4.2 : 3.2) + '" fill="' + stColor(h.s) + '" stroke="' + pal.ring + '" stroke-width="0.8"' + (exp ? '' : ' style="cursor:pointer"') + '>' +
        '<title>' + esc(p.label + ' ' + rl(h.s.pos, P[h.s.p]) + ' · ' + (st ? IXS[st][0] : 'loading') + ((h.s.ev || 0) ? ': ' + evLabel(h.s.ev) : '')) + '</title></circle>';
    }
    // Label sites with evidence; neighbours closer than a label's width share one (R179/R181/R183, or R12 +4).
    const groups = [];
    for (const h of heads.filter((x) => (x.s.ev || 0)).sort((a, b) => a.hx - b.hx)) {
      const last = groups[groups.length - 1];
      if (last && h.hx - last[0].hx < 34) last.push(h); else groups.push([h]);
    }
    // Two label lines: a label that would touch the previous one on the first line drops to the second.
    const lineEnd = [-1e9, -1e9];
    for (const grp of groups) {
      const text = grp.length <= 3 ? grp.map((h) => rl(h.s.pos, p)).join('/') : '' + rl(grp[0].s.pos, p) + ' +' + (grp.length - 1);
      const cx = (grp[0].hx + grp[grp.length - 1].hx) / 2, half = text.length * 2.6;
      const line = cx - half > lineEnd[0] + 6 ? 0 : cx - half > lineEnd[1] + 6 ? 1 : -1;
      if (line < 0) continue;
      lineEnd[line] = cx + half;
      const col = grp.some((h) => ixStatus(h.s) === 'exp') ? IXS_COLOR.exp : IXS_COLOR.pred;
      g += '<text x="' + r1(cx) + '" y="' + (headTop - 2 - line * 11) + '" font-size="9" text-anchor="middle" fill="' + col + '" font-weight="700">' + text + '</text>';
    }
    // Axis and legend.
    let yy = bottom + 6;
    const step = niceStep(len) * (panel && W < 500 ? 2 : 1);
    g += '<line x1="' + x0 + '" x2="' + x1 + '" y1="' + yy + '" y2="' + yy + '" stroke="' + pal.backbone + '"/>';
    const hOff = p.hist ? 1 : 0;
    for (let t = 0; t + hOff <= len; t += step) { const tx = r1(X(Math.max(t + hOff, 0.5))); g += '<line x1="' + tx + '" x2="' + tx + '" y1="' + yy + '" y2="' + (yy + 4) + '" stroke="' + pal.backbone + '"/><text x="' + tx + '" y="' + (yy + 15) + '" font-size="9.5" text-anchor="middle" fill="' + pal.muted + '">' + fmt(t) + '</text>'; }
    yy += 32;
    let lx = m;
    const lfs = panel ? 8.5 : 10, lcw = panel ? 4.7 : 5.3;
    const item = (draw, label) => {
      const w = 16 + label.length * lcw;
      if (lx > m && lx + w > W - m) { lx = m; yy += 15; }
      g += draw(lx, yy) + '<text x="' + (lx + 16) + '" y="' + (yy + 3.5) + '" font-size="' + lfs + '" fill="' + pal.muted + '">' + esc(label) + '</text>';
      lx += w + 12;
    };
    item((a, b) => '<circle cx="' + (a + 5) + '" cy="' + b + '" r="4.2" fill="' + IXS_COLOR.exp + '"/>', 'Methyl-R with experimental contact');
    item((a, b) => '<circle cx="' + (a + 5) + '" cy="' + b + '" r="4.2" fill="' + IXS_COLOR.pred + '"/>', 'Predicted only');
    item((a, b) => '<circle cx="' + (a + 5) + '" cy="' + b + '" r="3.2" fill="' + pal.none + '"/>', 'None recorded');
    item((a, b) => '<circle cx="' + (a + 5) + '" cy="' + b + '" r="4.2" fill="' + pal.bg + '" stroke="' + pal.ink + '" stroke-width="1.6"/>', 'Evidence at that methyl-R');
    item((a, b) => '<rect x="' + (a + 4) + '" y="' + (b - 6) + '" width="2" height="12" fill="' + pal.muted + '"/>', 'Crosslinked residue');
    item((a, b) => '<rect x="' + a + '" y="' + (b - 5) + '" width="12" height="10" fill="' + pal.disFill + '" opacity="0.6"/>', 'Disordered');
    if (kme.length) item((a, b) => '<rect x="' + (a + 2) + '" y="' + (b - 3) + '" width="6" height="6" transform="rotate(45 ' + (a + 5) + ' ' + b + ')" fill="' + KME_COLOR + '"/>', 'Lys methylation');
    yy += 12;
    const H = yy + 6;
    return { W, H, svg: '<svg xmlns="http://www.w3.org/2000/svg" class="ccart" width="' + W + '" height="' + H + '" viewBox="0 0 ' + W + ' ' + H + '" font-family=\'' + (pal.font || EXPORT_PAL.font).replace(/'/g, '"') + '\' role="img" aria-label="' + esc(p.label + ': methylarginines and contact evidence along the sequence') + '">' +
      (exp ? '<rect width="' + W + '" height="' + H + '" fill="' + pal.bg + '"/>' : '') + g + '</svg>' };
  }
  // The card on protein pages, redrawn as UniProt features, interaction data and PDBe-KB interfaces arrive.
  function contactCardHtml() {
    return '<div class="card" id="card-contacts"><div class="card-head"><div><h3>Contacts at a glance</h3><p class="sub">Each methylarginine, coloured by its interaction status, above the contact evidence on the same residue scale: RNA and DNA crosslinks in living cells, residues at PDB interfaces, and predicted-complex interfaces. A dashed line drops from every methylarginine with evidence at the residue itself; a ring marks which kind. Click a site for its summary.' + docLink('ix-status') + '</p></div>' +
      '<div class="exp" role="group" aria-label="Export"><button type="button" data-cc="svg">SVG</button><button type="button" data-cc="png">PNG</button></div></div>' +
      '<div class="fig-controls">' + panelControlsHtml() + '</div><div class="cc-box" id="cc-box"><p class="note">Loading…</p></div></div>';
  }
  function renderContactCard(p) {
    const box = $('#cc-box');
    if (!box) return;
    const ifd = p.iso ? { rna: [], dna: [], protein: [], self: [] } : IFDATA.get(p.acc.split('-')[0]) || null;
    box.innerHTML = '<div class="panel-prev"' + panelPreviewStyle() + '>' + contactCartoonSvg(p, UNIDATA.get(p.acc), ifd, Object.assign({ pal: Object.assign(screenPalette(), { bg: cssVar('--surface') }) }, panelOpts())).svg + '</div>';
    $$('.cc-site', box).forEach((c) => c.addEventListener('click', () => { UI.highlight = +c.dataset.pos; flashSite(+c.dataset.pos, true); drawMap(p, UNIDATA.get(p.acc)); }));
    const card = $('#card-contacts');
    if (card && !card.dataset.bound) {
      card.dataset.bound = '1';
      $$('[data-cc]', card).forEach((b) => b.addEventListener('click', () => {
        const f = IFDATA.get(p.acc.split('-')[0]) || null;
        exportPanel(b.dataset.cc, contactCartoonSvg(p, UNIDATA.get(p.acc), f, Object.assign({ forExport: true }, panelOpts())), p.label + '_contacts');
      }));
      bindPanelControls(card, () => renderContactCard(p));
      bindDocLinks(card);
    }
    if (!p.iso && !IFDATA.has(p.acc.split('-')[0])) loadIface(p.acc).then(() => { if (UI.detail === p.i) renderContactCard(p); }, () => {});
  }

  // Rows menu: overlays and rows with a checkbox each; rows also move up and down. Choices are remembered.
  function rowsMenuHtml() {
    const T = UI.fig.tracks, name = new Map(MAP_ROWS);
    return '<p class="flabel">On the backbone</p>' + MAP_OVERLAYS.map(([k, l]) => '<label class="check"><input type="checkbox" data-row="' + k + '"' + (T.has(k) ? ' checked' : '') + '> ' + esc(l) + '</label>').join('') +
      '<p class="flabel">Rows, top to bottom</p><ol class="rows-list">' + mapRowOrder().map((k, i, all) => '<li><label class="check"><input type="checkbox" data-row="' + k + '"' + (T.has(k) ? ' checked' : '') + '> ' + esc(name.get(k)) + '</label>' +
        '<span class="row-move"><button type="button" class="linkish" data-move="' + k + ':-1"' + (i ? '' : ' disabled') + ' aria-label="Move ' + esc(name.get(k)) + ' up">↑</button><button type="button" class="linkish" data-move="' + k + ':1"' + (i < all.length - 1 ? '' : ' disabled') + ' aria-label="Move ' + esc(name.get(k)) + ' down">↓</button></span></li>').join('') + '</ol>' +
      '<button class="linkish" type="button" data-rows-reset>Default rows</button>';
  }
  function bindRowsMenu(p) {
    const pop = $('#rows-pop');
    if (!pop) return;
    const redraw = () => { saveMapPrefs(); pop.innerHTML = rowsMenuHtml(); bindRowsMenu(p); drawMap(p, UNIDATA.get(p.acc)); renderContactCard(p); };
    $$('[data-row]', pop).forEach((cb) => cb.addEventListener('change', () => { if (cb.checked) UI.fig.tracks.add(cb.dataset.row); else UI.fig.tracks.delete(cb.dataset.row); redraw(); }));
    $$('[data-move]', pop).forEach((b) => b.addEventListener('click', () => {
      const [k, d] = b.dataset.move.split(':'), order = mapRowOrder(), i = order.indexOf(k), j = i + +d;
      if (j < 0 || j >= order.length) return;
      [order[i], order[j]] = [order[j], order[i]];
      UI.fig.order = order;
      redraw();
      const again = $('[data-move="' + b.dataset.move + '"]', pop);
      if (again && !again.disabled) again.focus();
    }));
    $('[data-rows-reset]', pop).addEventListener('click', () => { UI.fig.tracks = new Set(MAP_DEFAULT_ON); UI.fig.order = null; redraw(); });
  }
  function figureCard(p) {
    const f = UI.fig;
    return '<div class="card fig-card"><div class="card-head"><div><h3>' + esc(p.label) + ' protein map</h3><p class="sub" id="map-sub">' + esc((p.name ? p.name + ' · ' : '') + p.acc) + ' · loading UniProt features…</p></div>' +
      '<div class="exp" role="group" aria-label="Figure export"><button type="button" id="zoom-reset" hidden>Reset zoom</button><button type="button" id="fig-svg">SVG</button><button type="button" id="fig-png">PNG</button></div></div>' +
      '<div class="fig-controls"><label class="fig-ctl">Colour sites by <select id="fig-color">' + Object.keys(COLOR_MODES).map((k) => '<option value="' + k + '"' + (f.color === k ? ' selected' : '') + '>' + esc(COLOR_MODES[k].name) + '</option>').join('') + '</select></label>' +
      '<div class="seg" role="group" aria-label="Site labels"><span>Site labels</span>' + [['auto', 'Auto'], ['all', 'All'], ['none', 'None']].map(([k, l]) => '<button type="button" data-lab="' + k + '" aria-pressed="' + (f.labels === k) + '">' + l + '</button>').join('') + '</div>' +
      '<details class="col-menu rows-menu"><summary class="btn quiet small">Rows</summary><div class="col-pop" id="rows-pop">' + rowsMenuHtml() + '</div></details>' + '</div>' +
      '<div class="legend" id="fig-legend"></div><div class="map-box" id="map-box"></div>' +
      '<p class="note screen-only">Drag across the map to zoom. SVG and PNG exports reproduce this view on a white background with a title and legend.</p></div>';
  }

  function bindFigure(p) {
    const rerender = () => drawMap(p, UNIDATA.get(p.acc));
    $('#fig-color').addEventListener('change', (e) => { UI.fig.color = e.target.value; rerender(); styleStructure(); renderAfSide(p, AFDATA.get(p.acc)); });
    $$('[data-lab]').forEach((b) => b.addEventListener('click', () => { UI.fig.labels = b.dataset.lab; $$('[data-lab]').forEach((x) => x.setAttribute('aria-pressed', String(x === b))); rerender(); }));
    bindRowsMenu(p);
    const exportFig = (kind) => {
      const out = mapSvg(p, UNIDATA.get(p.acc), { pal: EXPORT_PAL, forExport: true });
      exportFigure(kind, out, (p.gene || p.acc) + '_methylarginine_map' + (UI.zoom ? '_' + out.v0 + '-' + out.v1 : ''));
    };
    $('#fig-svg').addEventListener('click', () => exportFig('svg'));
    $('#fig-png').addEventListener('click', () => exportFig('png'));
    $('#zoom-reset').addEventListener('click', () => { UI.zoom = null; drawMap(p, UNIDATA.get(p.acc)); });
  }

  const cssVar = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();

  function screenPalette() {
    const v = cssVar;
    return {
      bg: v('--surface'), ink: v('--ink'), muted: v('--muted'), grid: v('--line'), backbone: v('--line-strong'), stem: v('--line-strong'),
      backboneFill: v('--surface-2'), backboneEdge: v('--line-strong'), domainFill: v('--fig-domain'), domainInk: v('--fig-domain-ink'),
      disFill: v('--fig-dis'), disEdge: v('--line-strong'), regionFill: v('--fig-region'), biasFill: v('--fig-bias'), rgFill: v('--fig-rg'), rgInk: v('--fig-rg-ink'),
      dense: v('--dense'), denseInk: v('--accent-dark'), ptm: v('--ptm'), rtick: v('--line-strong'), ring: v('--surface'),
      site: v('--fig-site'), brush: v('--fig-brush'), accent: v('--accent'), variant: v('--muted'),
      cat: [v('--st-1'), v('--st-2'), v('--st-mix')], none: v('--st-none'), dark: v('--muted'),
      seq: [v('--seq-1'), v('--seq-2'), v('--seq-3'), v('--seq-4')], font: v('--font-body'),
    };
  }

  // Protein map: backbone with domain blocks, methylarginine lollipops and context tracks underneath.
  function mapSvg(p, u, o) {
    const pal = o.pal;
    const fig = UI.fig;
    const T = fig.tracks;
    const mode = COLOR_MODES[fig.color] || COLOR_MODES.state;
    const colors = mode.colors(pal);
    const len = (u && u.len) || p.len || Math.ceil(p.maxPos * 1.05) || 100;
    const feats = u ? u.feats : [];
    const z = UI.zoom || [1, len];
    const v0 = Math.max(1, z[0]), v1 = Math.min(len, Math.max(z[1], v0 + 10));
    const W = 1000, ML = 96, MR = 18;
    const x = (pos) => ML + ((pos - v0) / (v1 - v0 + 1)) * (W - ML - MR);
    const visible = (a, b) => b >= v0 && a <= v1;
    const cx = (pos) => Math.max(ML, Math.min(W - MR, x(pos)));
    const sites = p.sites.map((j) => S[j]).filter((s) => s.pos >= v0 && s.pos <= v1);
    const LABEL_GAP = 9.5;
    const labelsOn = fig.labels === 'all' || (fig.labels === 'auto' && sites.length <= 60 && sites.length * LABEL_GAP <= (W - ML - MR) * 0.9);
    const R0 = 5;
    const levels = [];
    const heads = sites.map((s) => {
      const hx = x(s.pos);
      let lv = 0;
      while (levels[lv] != null && hx - levels[lv] < 2 * R0 + 1 && lv < 9) lv++;
      levels[lv] = hx;
      return { s, hx, lv };
    });
    const nLv = Math.max(1, levels.length);
    const y0 = o.forExport ? 70 : 8;
    const lift = labelsOn ? 8 : 0;
    const headTop = y0 + (labelsOn ? 46 : 0);
    const base = headTop + 18 + lift + nLv * 12 + 12;
    let g = '';
    let y;
    const tlabel = (txt, yy) => { g += '<text x="0" y="' + (yy + 4) + '" font-size="10.5" font-weight="600" letter-spacing=".04em" fill="' + pal.muted + '">' + esc(txt.toUpperCase()) + '</text>'; };
    const tip = (html) => ' data-tip="' + esc(html) + '"';
    const gridLine = (yy) => { g += '<line x1="' + ML + '" x2="' + (W - MR) + '" y1="' + yy + '" y2="' + yy + '" stroke="' + pal.grid + '"/>'; };

    if (o.forExport) {
      g += '<text x="0" y="26" font-size="19" font-weight="700" fill="' + pal.ink + '">' + esc(p.label) + (p.name ? ' <tspan font-weight="400" fill="' + pal.muted + '" font-size="14">' + esc(p.name) + '</tspan>' : '') + '</text>' +
        '<text x="0" y="46" font-size="11.5" fill="' + pal.muted + '">' + esc(p.acc + ' · ' + fmt(len) + ' aa · ' + p.sites.length + ' methylarginine site' + (p.sites.length === 1 ? '' : 's') + (UI.zoom ? ' · showing residues ' + v0 + '–' + v1 : '')) + '</text>';
    }
    if (T.has('dense')) {
      for (const w of p.sig) {
        if (!visible(w[1], w[2])) continue;
        const a = cx(w[1]), b = cx(w[2]);
        g += '<rect class="feat" x="' + r1(a) + '" y="' + headTop + '" width="' + r1(Math.max(2, b - a)) + '" height="' + (base - headTop + 12) + '" rx="3" fill="' + pal.dense + '"' + tip('<b>Dense methyl-R window</b><br>' + w[3] + ' sites in ' + w[0] + ' aa (' + rl(w[1], p) + '–' + rn(w[2], p) + '), log2 obs/null ' + fmtL2(w[4]) + ', q ' + fmtP(w[5])) + '/>';
      }
    }
    tlabel('Sites', headTop + 10);
    g += '<rect x="' + r1(cx(1)) + '" y="' + (base - 4) + '" width="' + r1(cx(len) - cx(1)) + '" height="8" rx="4" fill="' + pal.backboneFill + '" stroke="' + pal.backboneEdge + '"/>';
    const below = [];
    if (T.has('domains')) {
      const doms = feats.filter((f) => ['Domain', 'Repeat', 'Zinc finger'].includes(f.type) && visible(f.start, f.end)).sort((a, b) => (a.type === 'Repeat') - (b.type === 'Repeat') || a.start - b.start);
      for (const f of doms) {
        const a = cx(f.start), b = cx(f.end), w = Math.max(2, b - a);
        const rep = f.type === 'Repeat';
        const h = rep ? 12 : 18;
        g += '<rect class="feat" x="' + r1(a) + '" y="' + (base - h / 2) + '" width="' + r1(w) + '" height="' + h + '" rx="4" fill="' + pal.domainFill + '"' + (rep ? ' opacity=".7"' : '') + tip('<b>' + esc(f.type) + '</b> ' + esc(f.desc) + '<br>' + f.start + '–' + f.end) + '/>';
        if (!rep) {
          const text = f.desc || f.type;
          if (w > text.length * 6 + 8) g += '<text x="' + r1(a + w / 2) + '" y="' + (base + 4) + '" font-size="10.5" font-weight="600" text-anchor="middle" fill="' + pal.domainInk + '" pointer-events="none">' + esc(text) + '</text>';
          else below.push([a + w / 2, text]);
        }
      }
    }
    // Site labels sit in a row above the heads, spread apart with leader lines so clusters stay readable.
    const labelX = new Map();
    if (labelsOn && heads.length) {
      const order = heads.slice().sort((a, b) => a.hx - b.hx);
      const gap = Math.min(LABEL_GAP, (W - ML - MR) / order.length);
      const L = order.map((h) => h.hx);
      for (let k = 1; k < L.length; k++) L[k] = Math.max(L[k], L[k - 1] + gap);
      if (L[L.length - 1] > W - MR - 4) {
        L[L.length - 1] = W - MR - 4;
        for (let k = L.length - 2; k >= 0; k--) L[k] = Math.min(L[k], L[k + 1] - gap);
      }
      order.forEach((h, k) => labelX.set(h, Math.max(ML, L[k])));
    }
    const labelY = headTop + 4;
    for (const h of heads) {
      const hy = headTop + 18 + lift + (nLv - 1 - h.lv) * 12;
      const col = colors[mode.cat(h.s)] || pal.none;
      const dim = !R.shown.has(h.s.j);
      const sel = UI.highlight === h.s.pos && !o.forExport;
      g += '<line x1="' + r1(h.hx) + '" x2="' + r1(h.hx) + '" y1="' + (base - 4) + '" y2="' + hy + '" stroke="' + pal.stem + '" stroke-width="1"' + (dim ? ' opacity=".4"' : '') + '/>';
      g += '<circle class="head" data-j="' + h.s.j + '" data-pos="' + h.s.pos + '" cx="' + r1(h.hx) + '" cy="' + hy + '" r="' + R0 + '" fill="' + col + '" stroke="' + (sel ? pal.ink : pal.ring) + '" stroke-width="' + (sel ? 2.6 : 2) + '"' + (dim ? ' opacity=".28"' : '') + '/>';
      if (labelsOn) {
        const lx = labelX.get(h);
        g += '<line x1="' + r1(h.hx) + '" y1="' + (hy - R0 - 1) + '" x2="' + r1(lx) + '" y2="' + (labelY + 3) + '" stroke="' + pal.stem + '" stroke-width=".7"' + (dim ? ' opacity=".4"' : '') + '/>' +
          '<text x="' + r1(lx + 1) + '" y="' + labelY + '" font-size="9.5" fill="' + pal.ink + '" transform="rotate(-55 ' + r1(lx + 1) + ' ' + labelY + ')"' + (dim ? ' opacity=".4"' : '') + '>' + rl(h.s.pos, P[h.s.p]) + '</text>';
      }
    }
    y = base + 22;
    if (below.length) {
      let right = -1e9;
      for (const [mx, text] of below.sort((a, b) => a[0] - b[0])) {
        const w = text.length * 5.4;
        if (mx - w / 2 > right + 6) { g += '<text x="' + r1(mx) + '" y="' + (y - 4) + '" font-size="9.5" text-anchor="middle" fill="' + pal.muted + '">' + esc(text) + '</text>'; right = mx + w / 2; }
      }
      y += 12;
    }
    // A row of interval features. Overlapping features (UniProt regions often nest, e.g. 1-186 and 1-117) are packed
    // onto separate lines so no box or label hides another; a label sits inside its box, else just after it when there
    // is room on that line, else only in the tooltip. Motifs are outlined to tell them from regions.
    const boxTrack = (name, list, fill, ink, labelOf, tipOf, dashed) => {
      if (!list.length) return;
      const items = list.slice().sort((p1, p2) => p1.start - p2.start || (p2.end - p2.start) - (p1.end - p1.start)).map((f) => {
        const a = cx(f.start), b = Math.max(cx(f.end + 1), a + 2), text = labelOf ? labelOf(f) : '', tw = text.length * 5.6;
        return { f, a, w: b - a, text, tw };
      });
      const lanes = [], MAXL = 5;
      for (const it of items) {
        const inside = it.text && it.w > it.tw + 6;
        const after = it.text && !inside && it.a + it.w + 8 + it.tw < W - MR;
        const reach = after ? it.a + it.w + 6 + it.tw : it.a + it.w;
        let k = lanes.findIndex((r) => r < it.a - 1);
        if (k < 0 && lanes.length < MAXL) { lanes.push(-1e9); k = lanes.length - 1; }
        if (k < 0) { k = lanes.indexOf(Math.min(...lanes)); it.text = ''; }
        else it.label = inside ? 'in' : after ? 'after' : it.text && it.w > 40 ? 'cut' : '';
        if (it.label === 'cut') it.text = it.text.slice(0, Math.max(3, Math.floor((it.w - 10) / 5.6))) + '…';
        // A label placed after a box may not run into the next box on its line.
        it.lane = k;
        lanes[k] = Math.max(lanes[k], it.label === 'after' ? reach : it.a + it.w);
      }
      tlabel(name, y);
      if (dashed) gridLine(y);
      for (const it of items) {
        const yy = y + it.lane * 14, motif = it.f.type === 'Motif';
        g += '<rect class="feat" x="' + r1(it.a) + '" y="' + (yy - 5) + '" width="' + r1(it.w) + '" height="10" rx="3" fill="' + (motif ? pal.bg : fill) + '"' +
          (dashed ? ' stroke="' + pal.disEdge + '" stroke-dasharray="3 2"' : motif ? ' stroke="' + pal.muted + '" stroke-width="1"' : '') + tip(tipOf(it.f)) + '/>';
        if (it.label === 'in' || it.label === 'cut') g += '<text x="' + r1(it.a + 3) + '" y="' + (yy + 3.5) + '" font-size="9.5" fill="' + ink + '" pointer-events="none">' + esc(it.text) + '</text>';
        else if (it.label === 'after') g += '<text x="' + r1(it.a + it.w + 4) + '" y="' + (yy + 3.5) + '" font-size="9.5" fill="' + ink + '">' + esc(it.text) + '</text>';
      }
      y += 20 + (lanes.length - 1) * 14;
    };
    // Rows under the backbone, drawn in the order the reader chose (Rows menu); each row skips itself when its data
    // is missing or still loading.
    const rows = {};
    rows.regions = () => boxTrack('Regions', feats.filter((f) => ((f.type === 'Region' && f.desc !== 'Disordered') || f.type === 'Motif') && visible(f.start, f.end)), pal.regionFill, pal.ink, (f) => f.desc, (f) => '<b>' + esc(f.type) + '</b> ' + esc(f.desc) + '<br>' + f.start + '–' + f.end);
    rows.disorder = () => {
      if (!u) return;
      const dis = feats.filter((f) => f.type === 'Region' && f.desc === 'Disordered' && visible(f.start, f.end));
      if (dis.length) boxTrack('Disorder', dis, pal.disFill, pal.muted, null, (f) => '<b>Disordered</b> (UniProt)<br>' + f.start + '–' + f.end, true);
      else { tlabel('Disorder', y); gridLine(y); y += 20; }
    };
    rows.bias = () => boxTrack('Composition', feats.filter((f) => f.type === 'Compositional bias' && visible(f.start, f.end)), pal.biasFill, pal.ink, (f) => f.desc.replace(/ residues$/, '-rich').replace('Basic and acidic', 'Charged'), (f) => '<b>Compositional bias</b> ' + esc(f.desc) + '<br>' + f.start + '–' + f.end);
    rows.rg = () => u && boxTrack('RG/RGG', u.rg.filter((b) => visible(b.start, b.end)).map((b) => Object.assign({}, b)), pal.rgFill, pal.rgInk, (b) => b.label, (b) => '<b>' + b.label + '</b> ' + b.start + '–' + b.end + '<br>' + b.n + ' RG dipeptides, ' + b.rgg + ' RGG<br>' + RG_HELP);
    rows.args = () => {
      if (!(u && u.seq)) return;
      tlabel('Arginines', y);
      const methyl = new Set(p.sites.map((j) => S[j].pos));
      for (let k = Math.max(0, v0 - 1); k < Math.min(u.seq.length, v1); k++) {
        if (u.seq[k] !== 'R' || methyl.has(k + 1)) continue;
        g += '<line x1="' + r1(x(k + 1)) + '" x2="' + r1(x(k + 1)) + '" y1="' + (y - 5) + '" y2="' + (y + 5) + '" stroke="' + pal.rtick + '"/>';
      }
      for (const s of sites) g += '<line x1="' + r1(x(s.pos)) + '" x2="' + r1(x(s.pos)) + '" y1="' + (y - 7) + '" y2="' + (y + 7) + '" stroke="' + (colors[mode.cat(s)] || pal.none) + '" stroke-width="2"/>';
      y += 20;
    };
    rows.variants = () => {
      if (!(u && (u.variants.length || u.mutagen.length))) return;
      const methyl = new Set(p.sites.map((j) => S[j].pos));
      const track = (name, list, isMut) => {
        list = list.filter((f) => visible(f.start, f.end));
        if (!list.length) return;
        tlabel(name, y);
        gridLine(y);
        let right = -1e9;
        const labelled = new Set();
        for (const f of list.sort((a, b) => a.start - b.start)) {
          const hit = f.from === 'R' && f.start === f.end && methyl.has(f.start);
          const xx = r1(x(f.start));
          const label = f.from + f.start + (f.to || '');
          const tipHtml = '<b>' + (isMut ? 'Mutagenesis' : 'Natural variant') + ' ' + esc(label) + '</b>' + (hit ? ' · at a methylarginine' : '') + '<br>' + esc(f.desc || '');
          if (isMut) g += '<rect class="feat" x="' + (xx - 3) + '" y="' + (y - 3) + '" width="6" height="6" transform="rotate(45 ' + xx + ' ' + y + ')" fill="' + (hit ? pal.accent : pal.bg) + '" stroke="' + (hit ? pal.accent : pal.variant) + '" stroke-width="1.2"' + tip(tipHtml) + '/>';
          else g += '<rect class="feat" x="' + (xx - (hit ? 1.5 : 1)) + '" y="' + (y - (hit ? 7 : 4)) + '" width="' + (hit ? 3 : 2) + '" height="' + (hit ? 14 : 8) + '" rx="1" fill="' + (hit ? pal.accent : pal.variant) + '"' + tip(tipHtml) + '/>';
          if (hit && +xx + 4 > right && !labelled.has(f.start)) {
            // Several substitutions at one site share a label (R514G/S).
            labelled.add(f.start);
            const tos = Array.from(new Set(list.filter((v) => v.start === f.start && v.end === f.end && v.from === 'R').map((v) => v.to).filter(Boolean)));
            const t = f.from + f.start + tos.join('/') + (shortDisease(f.desc) ? ' ' + shortDisease(f.desc).split(' and ')[0] : '');
            const tw = t.length * 5.2;
            const flip = +xx + 4 + tw > W - 2;
            g += '<text x="' + (flip ? +xx - 4 : +xx + 4) + '" y="' + (y - 5) + '" font-size="9" fill="' + pal.accent + '"' + (flip ? ' text-anchor="end"' : '') + '>' + esc(t) + '</text>';
            right = flip ? +xx : +xx + 4 + tw;
          }
        }
        y += 20;
      };
      track('Variants', u.variants, false);
      track('Mutagenesis', u.mutagen.filter((f) => f.from === 'R' || f.from.includes('R')), true);
    };
    // Other PTMs from the cross-PTM master table, one row per group; lysine methylation in its own colour.
    const ptmRow = (gi) => {
        const list = p.ptms.filter((t) => t[1] === gi && t[0] >= v0 && t[0] <= v1);
        if (!list.length) return;
        tlabel(PTM_LABEL[gi], y);
        gridLine(y);
        for (const t of list) g += '<rect class="feat" x="' + r1(x(t[0]) - 1.5) + '" y="' + (y - 5) + '" width="3" height="10" rx="1" fill="' + (gi === 3 ? KME_COLOR : pal.ptm) + '"' + tip('<b>' + PTM_PREFIX[gi] + t[2] + t[0] + '</b> · ' + PTM_LABEL[gi]) + '/>';
        y += 18;
    };
    // Residues at interfaces in PDB structures (PDBe-KB), one row per partner type.
    rows.iface = () => {
      const ifd = !p.iso ? IFDATA.get(p.acc.split('-')[0]) : null;
      if (!ifd) return;
      const methylPos = new Set(p.sites.map((j) => S[j].pos));
      for (const [key, label, color] of [['rna', 'RNA contact', IFACE_COLOR[4]], ['dna', 'DNA contact', IFACE_COLOR[8]], ['protein', 'Other protein', IFACE_COLOR[1]], ['self', 'Self-contact', IFACE_COLOR[2]]]) {
        const byPos = new Map();
        for (const [a, b, name] of ifd[key]) for (let r = Math.max(a, v0); r <= Math.min(b, v1); r++) { if (!byPos.has(r)) byPos.set(r, new Set()); byPos.get(r).add(name); }
        if (!byPos.size) continue;
        tlabel(label, y);
        gridLine(y);
        const positions = Array.from(byPos.keys()).sort((m, n) => m - n);
        let start = positions[0], prev = positions[0];
        const flush = (a, b) => {
          const names = new Set();
          let nm = 0;
          for (let r = a; r <= b; r++) { byPos.get(r).forEach((nme) => names.add(nme)); if (methylPos.has(r)) nm++; }
          const xa = x(a) - 0.8, w = Math.max(1.6, x(b + 1) - x(a));
          g += '<rect class="feat" x="' + r1(xa) + '" y="' + (y - 5) + '" width="' + r1(w) + '" height="10" rx="1" fill="' + color + '"' +
            tip('<b>' + label + '</b> ' + (a === b ? a : a + '–' + b) + (key === 'self' ? '<br>Touches another copy of this same protein (a homodimer or larger homo-oligomer)' : key === 'protein' ? '<br>Touches a different protein' : '') + '<br>' + esc(Array.from(names).slice(0, 6).join(', ')) + (names.size > 6 ? ' and ' + (names.size - 6) + ' more' : '') + (nm ? '<br>' + nm + ' methylarginine' + (nm === 1 ? '' : 's') + ' here' : '') + '<br>PDBe-KB') + '/>';
        };
        for (let k = 1; k <= positions.length; k++) {
          if (k < positions.length && positions[k] === prev + 1) { prev = positions[k]; continue; }
          flush(start, prev);
          if (k < positions.length) { start = prev = positions[k]; }
        }
        y += 18;
      }
    };
    rows.xl = () => {
      const xp = !p.iso && ixReady() ? IX.xlProt.get(p.i) : null;
      if (!xp) return;
      const tick = (pos, color, label) => '<rect class="feat" x="' + r1(x(pos) - 0.9) + '" y="' + (y - 6) + '" width="' + r1(Math.max(1.8, x(pos + 1) - x(pos))) + '" height="12" fill="' + color + '"' + tip(label) + '/>';
      const box = (a, b, color, op, label) => '<rect class="feat" x="' + r1(x(Math.max(a, v0)) - 0.5) + '" y="' + (y - 4) + '" width="' + r1(Math.max(1.5, x(Math.min(b, v1) + 1) - x(Math.max(a, v0)))) + '" height="8" rx="1" fill="' + color + '" fill-opacity="' + op + '"' + tip(label) + '/>';
      const inView = (a, b) => b >= v0 && a <= v1;
      if (xp.rnaPos.length || xp.rbd.length) {
        tlabel('RNA crosslinks', y);
        gridLine(y);
        for (const [a, b] of xp.rbd) if (inView(a, b)) g += box(a, b, IFACE_COLOR[4], 0.28, '<b>RBDmap RNA-bound fragment</b> ' + a + '–' + b + '<br>Castello et al. 2016, HeLa');
        xp.rnaPos.forEach((pos, k) => { if (pos >= v0 && pos <= v1) g += tick(pos, IFACE_COLOR[4], '<b>RNA-crosslinked residue</b> ' + pos + '<br>' + esc(xlSets(xp.rnaBits[k]).join(', '))); });
        y += 18;
      }
      if (xp.dnaPos.length || xp.dnaCand.length || xp.dnaPeps.length || xp.dnaBind.length) {
        tlabel('DNA crosslinks', y);
        gridLine(y);
        for (const [a, b] of xp.dnaBind) if (inView(a, b)) g += box(a, b, IFACE_COLOR[8], 0.22, '<b>UniProt DNA-binding region</b> ' + a + '–' + b);
        for (const [a, b] of xp.dnaPeps) if (inView(a, b)) g += box(a, b, IFACE_COLOR[8], 0.45, '<b>DNA-crosslinked peptide</b> ' + a + '–' + b + '<br>Trendel et al. 2025, MCF7');
        for (const pos of xp.dnaCand) if (pos >= v0 && pos <= v1) g += box(pos, pos, IFACE_COLOR[8], 0.75, '<b>DNA crosslink candidate</b> ' + pos + '<br>One of several residues the crosslink may sit on · Trendel et al. 2025, MCF7');
        for (const pos of xp.dnaPos) if (pos >= v0 && pos <= v1) g += tick(pos, IFACE_COLOR[8], '<b>DNA crosslink at this residue</b> ' + pos + '<br>Trendel et al. 2025, MCF7');
        y += 18;
      }
    };
    rows.phospho = () => ptmRow(0);
    rows.kacyl = () => ptmRow(1);
    rows.ub = () => ptmRow(2);
    rows.kme = () => ptmRow(3);
    // Mean residue charge of residues i-2 to i+2 (EMBOSS charge scale: K/R +1, H +0.5, D/E -1), drawn at residue i.
    rows.charge = () => {
      const ch = u && chargeProfile(u);
      if (!ch) return;
      tlabel('Charge', y + 2);
      const amp = 11, mid = y + 2;
      gridLine(mid);
      const n = v1 - v0 + 1, stepR = Math.max(1, Math.floor(n / 700));
      let pos = '', neg = '';
      for (let r = v0; r <= v1; r += stepR) {
        const v = ch[r - 1] == null ? 0 : ch[r - 1], xx = r1(x(r));
        pos += ' ' + xx + ',' + r1(mid - Math.max(0, v) * amp);
        neg += ' ' + xx + ',' + r1(mid - Math.min(0, v) * amp);
      }
      const ends = (pts) => r1(x(v0)) + ',' + mid + pts + ' ' + r1(x(v1)) + ',' + mid;
      g += '<polygon points="' + ends(pos) + '" fill="' + CHARGE_COLORS.basic + '" fill-opacity="0.75"/>' +
        '<polygon points="' + ends(neg) + '" fill="' + CHARGE_COLORS.acidic + '" fill-opacity="0.75"/>' +
        '<rect class="feat" x="' + ML + '" y="' + (mid - amp) + '" width="' + (W - ML - MR) + '" height="' + (2 * amp) + '" fill="transparent"' +
        tip('<b>Charge, sliding 5-residue mean</b> (EMBOSS charge: K and R +1, H +0.5, D and E −1)<br>Residues i−2 to i+2, drawn at residue i; near either end the mean is over the residues present. Blue above the line is basic, red below is acidic; full height = ±1') + '/>' +
        '<text x="' + (W - MR + 3) + '" y="' + (mid - amp + 6) + '" font-size="8" fill="' + pal.muted + '">+1</text>' +
        '<text x="' + (W - MR + 3) + '" y="' + (mid + amp) + '" font-size="8" fill="' + pal.muted + '">−1</text>';
      y += 30;
    };
    for (const k of mapRowOrder()) if (T.has(k) && rows[k]) rows[k]();
    // Axis
    y += 4;
    gridLine(y);
    const span = v1 - v0;
    const raw = span / 8;
    const mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const step = [1, 2, 5, 10].map((m) => m * mag).find((s) => s >= raw) || mag * 10;
    const hOff = p.hist ? 1 : 0;
    for (let t = Math.ceil((v0 - hOff) / step) * step + hOff; t <= v1; t += step) {
      g += '<line x1="' + r1(x(t)) + '" x2="' + r1(x(t)) + '" y1="' + y + '" y2="' + (y + 4) + '" stroke="' + pal.grid + '"/><text x="' + r1(x(t)) + '" y="' + (y + 15) + '" font-size="10" text-anchor="middle" fill="' + pal.muted + '">' + fmt(t - hOff) + '</text>';
    }
    tlabel(p.hist ? 'aa (histone)' : 'aa', y + 11);
    y += 24;
    if (o.forExport) {
      const ly = y + 8;
      g += '<text x="0" y="' + (ly + 4) + '" font-size="10.5" font-weight="700" fill="' + pal.muted + '">' + esc(mode.name) + '</text>';
      let lx = 8 + mode.name.length * 6.2;
      mode.cats.forEach((c, k) => { g += '<circle cx="' + r1(lx + 6) + '" cy="' + ly + '" r="5" fill="' + colors[k] + '"/><text x="' + r1(lx + 15) + '" y="' + (ly + 4) + '" font-size="10.5" fill="' + pal.ink + '">' + esc(c) + '</text>'; lx += 26 + c.length * 5.9; });
      let ex = 0;
      const ey = ly + 22;
      const item = (draw, lab) => { g += draw(ex, ey) + '<text x="' + (ex + 26) + '" y="' + (ey + 4) + '" font-size="10.5" fill="' + pal.ink + '">' + esc(lab) + '</text>'; ex += 38 + lab.length * 5.9; };
      const shown = (list) => list.some((f) => visible(f.start, f.end));
      if (T.has('domains') && shown(feats.filter((f) => ['Domain', 'Zinc finger', 'Repeat'].includes(f.type)))) item((a, b) => '<rect x="' + a + '" y="' + (b - 6) + '" width="20" height="12" rx="3" fill="' + pal.domainFill + '"/>', 'Domain / repeat');
      if (T.has('charge') && u && u.seq) item((a, b) => '<rect x="' + a + '" y="' + (b - 6) + '" width="20" height="6" fill="' + CHARGE_COLORS.basic + '"/><rect x="' + a + '" y="' + b + '" width="20" height="6" fill="' + CHARGE_COLORS.acidic + '"/>', 'Charge, 5-aa mean (+ / −)');
      if (T.has('kme') && p.ptms.some((t) => t[1] === 3 && visible(t[0], t[0]))) item((a, b) => '<rect x="' + (a + 8) + '" y="' + (b - 5) + '" width="3" height="10" rx="1" fill="' + KME_COLOR + '"/>', 'Lys methylation');
      if (T.has('xl') && ixReady() && IX.xlProt.get(p.i)) {
        item((a, b) => '<rect x="' + (a + 8) + '" y="' + (b - 6) + '" width="3" height="12" fill="' + IFACE_COLOR[4] + '"/>', 'RNA-crosslinked residue');
        item((a, b) => '<rect x="' + (a + 8) + '" y="' + (b - 6) + '" width="3" height="12" fill="' + IFACE_COLOR[8] + '"/>', 'DNA-crosslinked residue');
      }
      if (T.has('disorder') && shown(feats.filter((f) => f.type === 'Region' && f.desc === 'Disordered'))) item((a, b) => '<rect x="' + a + '" y="' + (b - 5) + '" width="20" height="10" rx="3" fill="' + pal.disFill + '" stroke="' + pal.disEdge + '" stroke-dasharray="3 2"/>', 'Disordered (UniProt)');
      if (T.has('rg') && u && shown(u.rg)) item((a, b) => '<rect x="' + a + '" y="' + (b - 5) + '" width="20" height="10" rx="3" fill="' + pal.rgFill + '"/>', 'RG/RGG-rich');
      if (T.has('dense') && p.sig.some((w) => visible(w[1], w[2]))) item((a, b) => '<rect x="' + a + '" y="' + (b - 5) + '" width="20" height="10" rx="3" fill="' + pal.dense + '"/>', 'Dense methyl-R window');
      if (T.has('variants') && u && u.variants.some((f) => f.from === 'R' && p.sites.some((j) => S[j].pos === f.start))) item((a, b) => '<rect x="' + (a + 9) + '" y="' + (b - 7) + '" width="3" height="14" rx="1" fill="' + pal.accent + '"/>', 'Variant at a methylarginine');
      y = ey;
      y += 18;
      g += '<text x="' + (W - MR) + '" y="' + (y + 4) + '" font-size="9" text-anchor="end" fill="' + pal.muted + '">Sites: methylarginine union, PTM_Pipeline · features and variants: UniProt' + (MAP_PTM_KEYS.some((k) => T.has(k)) && p.ptms.length ? ' · other PTMs incl. Lys methylation: PTM master table (not exhaustive)' : '') + (T.has('iface') && !p.iso && IFDATA.get(p.acc.split('-')[0]) ? ' · interfaces: PDBe-KB' : '') + (T.has('charge') && u && u.seq ? ' · charge: EMBOSS scale, 5-aa mean' : '') + '</text>';
      y += 12;
    }
    const H = Math.ceil(y + 4);
    const svg = '<svg xmlns="http://www.w3.org/2000/svg" class="pmap" width="' + W + '" height="' + H + '" viewBox="0 0 ' + W + ' ' + H + '" font-family=\'' + (pal.font || EXPORT_PAL.font).replace(/'/g, '"') + '\'' +
      ' role="img" aria-label="Map of ' + esc(p.label) + ' with ' + sites.length + ' methylated arginines">' +
      (o.forExport ? '<rect width="' + W + '" height="' + H + '" fill="' + pal.bg + '"/>' : '') + g +
      (o.forExport ? '' : '<rect id="brush" x="0" y="0" width="0" height="' + H + '" fill="' + pal.brush + '" stroke="' + pal.accent + '" hidden/>') + '</svg>';
    return { svg, W, H, ML, MR, v0, v1, len, mode, colors };
  }

  function drawMap(p, u) {
    const box = $('#map-box');
    if (!box) return;
    const sub = $('#map-sub');
    const out = mapSvg(p, u, { pal: screenPalette(), forExport: false });
    if (sub) {
      const who = (p.name || (u && u.name) || '') + (p.name || (u && u.name) ? ' · ' : '') + p.acc + ' · ' + (p.hist ? 'histone numbering (no initiator Met; R2 = UniProt R3) · ' : '');
      sub.textContent = u ? who + fmt(out.len) + ' aa · domains, regions, disorder and variants from UniProt' + (MAP_PTM_KEYS.some((k) => UI.fig.tracks.has(k)) && p.ptms.length ? ' · other PTMs from the PTM master table (not exhaustive)' : '') + (UI.fig.tracks.has('iface') && IFDATA.get(p.acc.split('-')[0]) ? ' · PDB interfaces from PDBe-KB' : '') + (UI.zoom ? ' · showing ' + out.v0 + '–' + out.v1 : '') :
        (uniprotCache.has(p.acc) && UNIDATA.has(p.acc) ? 'UniProt features unavailable; showing sites on a ' + (p.len ? fmt(out.len) + ' aa' : 'length-estimated') + ' backbone.' : 'Loading UniProt features…');
    }
    $('#zoom-reset').hidden = !UI.zoom;
    $('#fig-legend').innerHTML = '<span class="legend-title">' + esc(out.mode.name) + '</span>' + out.mode.cats.map((c, k) => '<span><i style="background:' + out.colors[k] + '"></i>' + esc(c) + '</span>').join('') +
      '<span><i class="sq" style="background:var(--fig-domain)"></i>Domain</span>' +
      (UI.fig.tracks.has('iface') && out.svg.includes('>SELF-CONTACT<') ? '<span class="legend-note">Self-contact: the residue touches another copy of this protein (homodimer or homo-oligomer)</span>' : '');
    box.innerHTML = out.svg;
    bindMap(box, p, u, out);
  }

  // ------------------------------------------------------------- 3D structure
  function structureCard(p) {
    return '<div class="card" id="card-af"><div class="card-head"><div><h3>3D structure</h3><p class="sub" id="af-sub">The AlphaFold DB model, or any experimental structure listed below. Methylarginine side chains use the protein-map colours. Click any residue, or a row in the site table, to zoom to it and show its contacts; the options change the style, colouring, surface and lighting.</p></div>' +
      '<div class="exp" role="group" aria-label="Structure controls"><button type="button" id="af-reset" disabled>Reset view</button><button type="button" id="af-spin" disabled aria-pressed="false">Spin</button><button type="button" id="af-png" disabled>PNG</button></div></div>' +
      '<div class="af-controls" id="af-controls"></div>' +
      '<div class="af-wrap"><div><div class="af-viewer" id="af-viewer"><div class="af-placeholder">Loading the AlphaFold model…</div></div><p class="note af-msg" id="af-msg" aria-live="polite"></p></div><div class="af-side" id="af-side"></div></div>' +
      '<div class="pdb-block"><div class="card-head"><div><h3 class="h4">Experimental structures (PDB)</h3><p class="sub" id="pdb-sub">Looking up PDB entries…</p></div>' +
      '<div class="exp"><a class="exp-link" href="https://www.ebi.ac.uk/pdbe/pdbe-kb/proteins/' + encodeURIComponent(p.acc.split('-')[0]) + '" target="_blank" rel="noopener">PDBe-KB ↗</a></div></div><div id="pdb-list"></div></div>' +
      '<div class="pdb-block" id="hom-block"><div class="card-head"><div><h3 class="h4">Homologous experimental structures</h3><p class="sub" id="hom-sub">Searching the PDB for related proteins…</p></div></div><div id="hom-list"></div></div>' +
      '<div class="pdb-block" id="pred-block"><div class="card-head"><div><h3 class="h4">Predicted complexes</h3><p class="sub" id="pred-sub">Looking up predicted complexes…</p></div></div><div id="pred-list"></div></div></div>';
  }

  function disposeStructure() {
    if (!mol) return;
    try {
      const c = mol.el && mol.el.querySelector('canvas');
      const gl = c && (c.getContext('webgl2') || c.getContext('webgl'));
      const ext = gl && gl.getExtension('WEBGL_lose_context');
      if (ext) ext.loseContext();
    } catch (e) { /* ignore */ }
    mol = null;
  }

  function renderAfSide(p, af) {
    const side = $('#af-side');
    if (!side) return;
    if (af === undefined) { side.innerHTML = '<p class="note">Looking up AlphaFold DB…</p>'; return; }
    const legend = (counts) => '<ul class="af-legend">' + PLDDT_BANDS.map((b, k) => '<li><i style="background:' + b[1] + '"></i>' + esc(b[2]) + (counts ? '<b class="num">' + counts[k] + '</b>' : '') + '</li>').join('') + '</ul>';
    if (mol && mol.mode === 'pdb' && mol.p.i === p.i) { renderPdbSide(p); return; }
    if (!af) { side.innerHTML = '<p class="note">No AlphaFold DB model for ' + esc(p.acc) + '. Experimental structures, if any, are listed below.</p>'; return; }
    const counts = [0, 0, 0, 0];
    let n = 0, amSum = 0, amN = 0, amLp = 0;
    for (const j of p.sites) {
      const v = af.plddt.get(S[j].pos);
      if (v != null) { counts[plddtBand(v)]++; n++; }
      const a = af.am.get(S[j].pos);
      if (a && a.n) { const m = a.sum / a.n; amSum += m; amN++; if (m >= 0.564) amLp++; }
    }
    const pal = screenPalette();
    const mode = COLOR_MODES[UI.fig.color] || COLOR_MODES.state;
    const colors = mode.colors(pal);
    side.innerHTML = '<dl class="kv"><dt>Model</dt><dd><a href="https://alphafold.ebi.ac.uk/entry/' + encodeURIComponent(p.acc.split('-')[0]) + '" target="_blank" rel="noopener">' + esc(af.entry.entryId) + ' v' + esc(af.entry.latestVersion || '') + ' ↗</a></dd>' +
      '<dt>Mean pLDDT</dt><dd>' + (af.entry.globalMetricValue != null ? af.entry.globalMetricValue.toFixed(1) : '—') + '</dd>' +
      (af.entry.paeImageUrl ? '<dt>PAE</dt><dd><a href="' + esc(af.entry.paeImageUrl) + '" target="_blank" rel="noopener">Predicted aligned error ↗</a></dd>' : '') + '</dl>' +
      '<p class="flabel">Methylarginines by model confidence (pLDDT)</p>' +
      (n ? '<div class="af-stack">' + counts.map((c, k) => (c ? '<i style="flex:' + c + ';background:' + PLDDT_BANDS[k][1] + '" title="' + esc(PLDDT_BANDS[k][2]) + ': ' + c + '"></i>' : '')).join('') + '</div>' : '') +
      legend(n ? counts : null) +
      (amN ? '<p class="note">AlphaMissense: mean R→X pathogenicity ' + (amSum / amN).toFixed(2) + ' across sites; ' + amLp + ' of ' + amN + ' sites average “likely pathogenic”.</p>' : '') +
      '<p class="flabel screen-only">Focus a site</p><div class="chips af-sites screen-only">' + p.sites.map((j) => S[j]).map((s) => '<button class="chip" type="button" data-pos="' + s.pos + '"><i class="dot" style="background:' + (colors[mode.cat(s)] || pal.none) + '"></i>' + rl(s.pos, P[s.p]) + '</button>').join('') + '</div>';
    $$('.af-sites .chip', side).forEach((b) => b.addEventListener('click', () => { UI.highlight = +b.dataset.pos; flashSite(+b.dataset.pos, false); }));
  }

  function initStructure(p) {
    const holder = $('#af-viewer');
    if (!holder) return;
    renderAfSide(p, AFDATA.get(p.acc));
    const enable = () => ['#af-reset', '#af-spin', '#af-png'].forEach((id) => { const b = $(id); if (b) b.disabled = false; });
    const bindButtons = () => {
      $('#af-reset').onclick = () => {
        if (!mol) return;
        if (mol.label) { mol.viewer.removeLabel(mol.label); mol.label = null; }
        mol.focus = null;
        styleStructure(false);
        mol.viewer.zoomTo();
        mol.viewer.render();
        const msg = $('#af-msg');
        if (msg) msg.textContent = '';
      };
      $('#af-spin').onclick = (e) => { if (!mol) return; mol.spin = !mol.spin; mol.viewer.spin(mol.spin ? 'y' : false); e.currentTarget.setAttribute('aria-pressed', String(mol.spin)); };
      $('#af-png').onclick = () => {
        if (!mol) return;
        const src = mol.mode === 'pdb' ? 'PDB ' + mol.pdb.e.id.toUpperCase() : 'AlphaFold DB model';
        pngWithCredit(mol.viewer.pngURI(), (p.gene || p.acc) + '_' + (mol.mode === 'pdb' ? mol.pdb.e.id.toUpperCase() : 'alphafold') + '_methylarginines.png', p.label + ' (' + p.acc + ') · ' + src + ' · methylarginines highlighted');
      };
    };
    // Re-renders of the same protein keep the live viewer instead of rebuilding it.
    if (mol && mol.p.i === p.i && mol.el) {
      holder.replaceWith(mol.el);
      mol.el.id = 'af-viewer';
      mol.viewer.resize();
      renderViewControls();
      styleStructure(false);
      renderAfSide(p, AFDATA.get(p.acc));
      renderPdbList(p);
      renderHomologList(p);
      renderPredList(p);
      enable(); bindButtons();
      if (mol.spin) $('#af-spin').setAttribute('aria-pressed', 'true');
      UI.startStructure = () => {};
      return;
    }
    disposeStructure();
    renderViewControls();
    const start = () => { renderPdbList(p); renderHomologList(p); renderPredList(p); return startAf(); };
    const startAf = () => loadAF(p.acc).then((af) => {
      if (UI.detail !== p.i || !document.body.contains(holder)) return;
      renderAfSide(p, af);
      if (!af) {
        if (!mol || mol.p.i !== p.i) holder.innerHTML = '<div class="af-placeholder">No AlphaFold DB model for ' + esc(p.acc) + '. AlphaFold DB does not cover every isoform, and for human proteins longer than 2,700 residues it offers only fragment downloads. Experimental structures, if any, can be viewed from the list below.</div>';
        return;
      }
      if (mol && mol.p.i === p.i) return;
      holder.innerHTML = '<div class="af-placeholder">Loading the 3D viewer…</div>';
      return Promise.all([load3Dmol(), fetchText(af.entry.pdbUrl)]).then(([lib, pdb]) => {
        if (UI.detail !== p.i || !document.body.contains(holder)) return;
        holder.innerHTML = '';
        const viewer = lib.createViewer(holder, { backgroundColor: cssVar('--surface-2') || '#eef5f7', antialias: true });
        viewer.addModel(pdb, 'pdb');
        mol = { viewer, p, el: holder, label: null, spin: false, mode: 'af', afPdb: pdb };
        renderViewControls();
        styleStructure();
        viewer.zoomTo();
        viewer.render();
        enable(); bindButtons();
        if (UI.highlight != null) focusStructure(UI.highlight);
      });
    }).catch(() => {
      holder.innerHTML = '<div class="af-placeholder">The 3D viewer could not load. <a href="https://alphafold.ebi.ac.uk/entry/' + encodeURIComponent(p.acc.split('-')[0]) + '" target="_blank" rel="noopener">Open the model in AlphaFold DB ↗</a></div>';
    });
    let started = false;
    const startOnce = () => { if (!started) { started = true; start(); } };
    UI.startStructure = startOnce;
    if ('IntersectionObserver' in window) {
      const io = new IntersectionObserver((entries) => { if (entries.some((e) => e.isIntersecting)) { io.disconnect(); startOnce(); } }, { rootMargin: '300px' });
      io.observe(holder);
    } else startOnce();
  }

  // 3D view options shared by the AlphaFold model and PDB entries; remembered between visits.
  const VIEW_STYLES = [['cartoon', 'Cartoon'], ['trace', 'Trace'], ['sticks', 'Sticks'], ['spheres', 'Spheres']];
  const VIEW_COLORS = [
    ['plddt', 'Confidence (pLDDT)', 'af'], ['focus', 'This protein\'s chain', 'pdb'], ['chain', 'By chain', 'pdb'],
    ['rainbow', 'Rainbow (N→C)'], ['ss', 'Secondary structure'], ['charge', 'Charge'], ['hydro', 'Hydrophobicity'], ['xl', 'RNA/DNA crosslinks (cells)'], ['uniform', 'Uniform grey'],
  ];
  const XL_COLORS = { rna: '#7B4FBF', dna: '#B7791F', both: '#C2410C', rbd: '#CDB8EE', dnaReg: '#EBD3A6', none: '#C9D2D8', other: '#E6EAED' };
  const VIEW_SURFACES = [['off', 'Off'], ['glass', 'Transparent'], ['solid', 'Opaque']];
  const VIEW_BGS = [['page', 'Page'], ['white', 'White'], ['black', 'Black']];
  const VIEW_LIGHTS = [['ao', 'Ambient occlusion'], ['aoline', 'Occlusion + outline'], ['outline', 'Outline'], ['plain', 'Standard']];
  const VIEW_TOGGLES = [
    ['contacts', 'Contacts (5 Å)', 'Show side chains within 5 Å of the focused residue'],
    ['labels', 'Site labels', 'Label every resolved methylarginine'],
    ['others', 'Other chains', 'Show the other chains of the PDB entry', 'pdb'],
    ['ligands', 'Ligands', 'Show ligands and other hetero groups (waters are always hidden)', 'pdb'],
  ];
  const VIEW_DEFAULT = { style: 'cartoon', color: 'auto', surface: 'off', light: 'ao', bg: 'page', contacts: true, labels: false, others: true, ligands: true };
  UI.mv = (() => {
    const v = Object.assign({}, VIEW_DEFAULT);
    let saved = {};
    try { saved = JSON.parse(store.get('mx-view') || '{}') || {}; } catch (e) { /* ignore */ }
    const ok = { style: VIEW_STYLES, color: VIEW_COLORS, surface: VIEW_SURFACES, light: VIEW_LIGHTS, bg: VIEW_BGS };
    for (const k of Object.keys(v)) {
      if (!(k in saved)) continue;
      if (ok[k] ? ok[k].some((o) => o[0] === saved[k]) : typeof saved[k] === 'boolean') v[k] = saved[k];
    }
    return v;
  })();
  // Kyte & Doolittle (1982) hydropathy.
  const KD = { ILE: 4.5, VAL: 4.2, LEU: 3.8, PHE: 2.8, CYS: 2.5, MET: 1.9, ALA: 1.8, GLY: -0.4, THR: -0.7, SER: -0.8, TRP: -0.9, TYR: -1.3, PRO: -1.6, HIS: -3.2, GLU: -3.5, GLN: -3.5, ASP: -3.5, ASN: -3.5, LYS: -3.9, ARG: -4.5 };
  const SS_COLORS = { h: '#D9466F', s: '#E8B730', c: '#B7C2CA' };
  const CHARGE_COLORS = { basic: '#2F6FD6', his: '#8DA9E8', acidic: '#D64545', other: '#CDD5DB' };
  const RAINBOW = ['#3B4CC0', '#2C9FD9', '#3DBE7A', '#F2C230', '#E8543B'];
  const HYDRO = ['#2F6FA7', '#F1F3F5', '#C8742B'];
  const CHAIN_COLORS = ['#0E7490', '#E07A2E', '#6C5CE7', '#2E9E5B', '#C2417A', '#8C6D1F', '#3D7DD8', '#B33939', '#5B8C85', '#9B59B6'];
  const ELEM_GREY = { prop: 'elem', map: { C: '#A3AEB8', N: '#3B5BDB', O: '#E03131', S: '#E6B800', P: '#F08C00', SE: '#F08C00' } };
  const ELEM_FOCUS = { prop: 'elem', map: { C: '#F59E0B', N: '#3B5BDB', O: '#E03131', S: '#E6B800' } };
  const WATER = ['HOH', 'WAT', 'DOD'];
  const mixHex = (a, b, t) => {
    const pa = parseInt(a.slice(1), 16), pb = parseInt(b.slice(1), 16);
    const ch = (s) => Math.round(((pa >> s) & 255) * (1 - t) + ((pb >> s) & 255) * t);
    return '#' + ((1 << 24) | (ch(16) << 16) | (ch(8) << 8) | ch(0)).toString(16).slice(1);
  };
  const ramp = (stops, t) => {
    const n = stops.length - 1;
    const u = Math.max(0, Math.min(1, t)) * n;
    const k = Math.min(n - 1, Math.floor(u));
    return mixHex(stops[k], stops[k + 1], u - k);
  };
  const baseResn = (r) => (METHYL_R_CODES.includes(r) ? 'ARG' : r);
  const isPolymer = (a) => !a.hetflag || METHYL_R_CODES.includes(a.resn);

  // Colour options that apply to the loaded model: pLDDT for AlphaFold models (B-factors in predicted complexes).
  function viewColorOk(o) {
    const pdbMode = !!(mol && mol.mode === 'pdb');
    return !o[2] || o[2] === (pdbMode ? 'pdb' : 'af') || (o[0] === 'plddt' && pdbMode && !!mol.pdb.pred);
  }
  function viewColorKey() {
    const pdbMode = !!(mol && mol.mode === 'pdb');
    const opt = VIEW_COLORS.find((o) => o[0] === UI.mv.color);
    if (!opt || !viewColorOk(opt)) return pdbMode ? 'focus' : 'plddt';
    return opt[0];
  }

  // Per-chain residue ranges and chain order of the loaded model, for rainbow and chain colouring.
  function modelChains() {
    if (mol.chains) return mol.chains;
    const ranges = new Map();
    for (const a of mol.viewer.selectedAtoms({})) {
      if (!isPolymer(a)) continue;
      const r = ranges.get(a.chain);
      if (!r) ranges.set(a.chain, [a.resi, a.resi]);
      else { if (a.resi < r[0]) r[0] = a.resi; if (a.resi > r[1]) r[1] = a.resi; }
    }
    mol.chains = { ranges, order: new Map(Array.from(ranges.keys()).map((c, k) => [c, k])) };
    return mol.chains;
  }

  function viewColorFn(key) {
    const grey = '#9AA6B2';
    if (key === 'plddt' && mol.mode === 'pdb') return (a) => PLDDT_BANDS[plddtBand(a.b)][1];
    if (key === 'plddt') {
      const af = AFDATA.get(mol.p.acc);
      return (a) => { const v = af && af.plddt.has(a.resi) ? af.plddt.get(a.resi) : a.b; return v == null ? grey : PLDDT_BANDS[plddtBand(v)][1]; };
    }
    if (key === 'focus') { const c = mol.pdb.chain; return (a) => (a.chain === c ? '#0E7490' : '#C7D4D9'); }
    if (key === 'chain') { const { order } = modelChains(); return (a) => CHAIN_COLORS[(order.get(a.chain) || 0) % CHAIN_COLORS.length]; }
    if (key === 'rainbow') {
      const { ranges } = modelChains();
      return (a) => { const r = ranges.get(a.chain); return r ? ramp(RAINBOW, (a.resi - r[0]) / Math.max(1, r[1] - r[0])) : grey; };
    }
    if (key === 'ss') return (a) => SS_COLORS[a.ss] || SS_COLORS.c;
    if (key === 'charge') {
      return (a) => { const r = baseResn(a.resn); return r === 'ARG' || r === 'LYS' ? CHARGE_COLORS.basic : r === 'HIS' ? CHARGE_COLORS.his : r === 'ASP' || r === 'GLU' ? CHARGE_COLORS.acidic : CHARGE_COLORS.other; };
    }
    if (key === 'hydro') return (a) => { const v = KD[baseResn(a.resn)]; return v == null ? CHARGE_COLORS.other : ramp(HYDRO, (v + 4.5) / 9); };
    if (key === 'xl') {
      // This protein's residues in UniProt numbering (AlphaFold and predicted models directly, PDB entries through SIFTS).
      const xp = ixReady() ? IX.xlProt.get(mol.p.i) : null;
      const rna = new Set(xp ? xp.rnaPos : []), dna = new Set(xp ? xp.dnaPos : []);
      const inR = (pos, list) => list.some(([a, b]) => pos >= a && pos <= b);
      const pdbMode = mol.mode === 'pdb';
      const cache = new Map();
      const unpOf = (atom) => {
        if (!pdbMode) return atom.resi;
        if (atom.chain !== mol.pdb.chain) return null;
        if (mol.pdb.pred) return atom.resi;
        if (!cache.has(atom.resi)) cache.set(atom.resi, authToUnp(atom.resi));
        return cache.get(atom.resi);
      };
      return (a) => {
        const pos = unpOf(a);
        if (pos == null) return XL_COLORS.other;
        if (!xp) return XL_COLORS.none;
        const r = rna.has(pos), d = dna.has(pos);
        if (r && d) return XL_COLORS.both;
        if (r) return XL_COLORS.rna;
        if (d) return XL_COLORS.dna;
        if (inR(pos, xp.rbd)) return XL_COLORS.rbd;
        if (inR(pos, xp.dnaPeps) || inR(pos, xp.dnaBind)) return XL_COLORS.dnaReg;
        return XL_COLORS.none;
      };
    }
    return () => grey;
  }

  function viewColorKeyHtml(key) {
    const sw = (c, l) => '<span><i class="sq" style="background:' + c + '"></i>' + esc(l) + '</span>';
    const grad = (stops, a, b) => '<span>' + esc(a) + '<i class="grad" style="background:linear-gradient(90deg,' + stops.join(',') + ')"></i>' + esc(b) + '</span>';
    const parts = {
      plddt: () => PLDDT_BANDS.map((b) => sw(b[1], b[2])).join(''),
      focus: () => sw('#0E7490', 'Chain ' + mol.pdb.chain + (mol.pdb.homolog ? ' (homolog)' : ' (this protein)')) + sw('#C7D4D9', mol.pdb.pred ? 'Partner: ' + mol.pdb.pred.partnerLabel : 'Other chains'),
      chain: () => Array.from(modelChains().order.keys()).slice(0, 8).map((c, k) => sw(CHAIN_COLORS[k % CHAIN_COLORS.length], 'Chain ' + c + (mol.pdb && c === mol.pdb.chain ? ' (this protein)' : mol.pdb && mol.pdb.pred && mol.pdb.pred.names[c] ? ' (' + mol.pdb.pred.names[c] + ')' : ''))).join('') + (modelChains().order.size > 8 ? '<span class="muted">+' + (modelChains().order.size - 8) + ' more</span>' : ''),
      rainbow: () => grad(RAINBOW, 'N-terminus', 'C-terminus'),
      ss: () => sw(SS_COLORS.h, 'Helix') + sw(SS_COLORS.s, 'Strand') + sw(SS_COLORS.c, 'Coil'),
      charge: () => sw(CHARGE_COLORS.basic, 'Basic (R, K)') + sw(CHARGE_COLORS.his, 'His') + sw(CHARGE_COLORS.acidic, 'Acidic (D, E)') + sw(CHARGE_COLORS.other, 'Other'),
      hydro: () => grad(HYDRO, 'Hydrophilic', 'Hydrophobic') + '<span class="muted">Kyte–Doolittle</span>',
      xl: () => (ixReady() && IX.xlProt.get(mol.p.i) ? sw(XL_COLORS.rna, 'RNA-crosslinked residue') + sw(XL_COLORS.dna, 'DNA-crosslinked residue') + sw(XL_COLORS.both, 'Both') + sw(XL_COLORS.rbd, 'RBDmap RNA-bound fragment') + sw(XL_COLORS.dnaReg, 'DNA crosslink candidate, crosslinked peptide or DNA-binding region') + sw(XL_COLORS.none, 'No crosslink') + (mol.mode === 'pdb' && !mol.pdb.pred ? sw(XL_COLORS.other, 'Other chains') : '') :
        '<span class="muted">No RNA or DNA crosslinks recorded for this protein</span>'),
      uniform: () => '',
    };
    return '<span class="legend-title">' + esc((VIEW_COLORS.find((o) => o[0] === key) || ['', ''])[1]) + '</span>' + (parts[key] || parts.uniform)() +
      '<span><i style="background:var(--fig-site)"></i>Methylarginines: map colours (this colour where the map shows grey)</span>' +
      (mol && mol.mode === 'pdb' && mol.naKinds ? Array.from(mol.naKinds).map((k) => '<span><i class="sq" style="background:' + NA_COLOR[k] + '"></i>' + k + '</span>').join('') : '') +
      (UI.mv.contacts && mol && mol.focus ? '<span><i class="sq" style="background:' + ELEM_GREY.map.C + '"></i>Contacts within 5 Å</span>' : '');
  }

  function renderViewControls() {
    const box = $('#af-controls');
    if (!box) return;
    const V = UI.mv;
    const pdbMode = !!(mol && mol.mode === 'pdb');
    const key = viewColorKey();
    const sel = (id, label, list, cur) => '<label class="fig-ctl">' + label + ' <select id="' + id + '">' + list.map((o) => '<option value="' + o[0] + '"' + (o[0] === cur ? ' selected' : '') + '>' + esc(o[1]) + '</option>').join('') + '</select></label>';
    box.innerHTML = '<div class="fig-controls">' +
      sel('vw-style', 'Style', VIEW_STYLES, V.style) +
      sel('vw-color', 'Colour', VIEW_COLORS.filter(viewColorOk), key) +
      sel('vw-surface', 'Surface', VIEW_SURFACES, V.surface) +
      sel('vw-light', 'Lighting', VIEW_LIGHTS, V.light) +
      sel('vw-bg', 'Background', VIEW_BGS, V.bg) + '</div>' +
      '<div class="chips af-toggles" aria-label="Show">' + VIEW_TOGGLES.filter((t) => !t[3] || pdbMode).map((t) => '<button class="chip" type="button" data-vw="' + t[0] + '" aria-pressed="' + !!V[t[0]] + '" title="' + esc(t[2]) + '">' + esc(t[1]) + '</button>').join('') +
      '<button class="linkish" type="button" id="vw-default">Default view</button></div>' +
      '<div class="legend af-key" id="af-key">' + (mol ? viewColorKeyHtml(key) : '') + '</div>';
    const save = () => store.set('mx-view', JSON.stringify(V));
    const change = (k, val, reSurface) => { V[k] = val; save(); if (mol) styleStructure(reSurface); renderViewControls(); };
    $('#vw-style').onchange = (e) => change('style', e.target.value, false);
    $('#vw-color').onchange = (e) => change('color', e.target.value, true);
    $('#vw-surface').onchange = (e) => change('surface', e.target.value, true);
    $('#vw-light').onchange = (e) => change('light', e.target.value, false);
    $('#vw-bg').onchange = (e) => change('bg', e.target.value, false);
    $$('[data-vw]', box).forEach((b) => { b.onclick = () => change(b.dataset.vw, !V[b.dataset.vw], false); });
    $('#vw-default').onclick = () => { Object.assign(V, VIEW_DEFAULT); store.set('mx-view', null); if (mol) styleStructure(true); renderViewControls(); };
  }

  function viewBackground() {
    return UI.mv.bg === 'white' ? '#FFFFFF' : UI.mv.bg === 'black' ? '#000000' : (cssVar('--surface-2') || '#eef5f7');
  }

  function styleStructure(reSurface) {
    if (!mol || mol.p.i !== UI.detail) return;
    const { viewer, p } = mol;
    const V = UI.mv;
    const pdbMode = mol.mode === 'pdb';
    const pdb = mol.pdb;
    const af = AFDATA.get(p.acc);
    const inChain = (a) => !pdbMode || a.chain === pdb.chain;
    viewer.setBackgroundColor(viewBackground());
    // Screen-space ambient occlusion and/or silhouette outlines (WebGL2; ignored by WebGL1 browsers).
    const light = { ao: 'ambientOcclusion', aoline: 'ambientOcclusion outline', outline: 'outline' }[V.light] || '';
    viewer.setViewStyle({ style: light, strength: V.style === 'spheres' || V.surface === 'solid' ? 1.6 : 1.9, color: UI.mv.bg === 'black' ? '#FFFFFF' : '#000000', width: 0.08 });

    // Methylarginine positions in model numbering with their map colours.
    const pal = screenPalette();
    const mode = COLOR_MODES[UI.fig.color] || COLOR_MODES.state;
    const colors = mode.colors(pal);
    const siteColor = new Map();
    for (const j of p.sites) {
      const r = pdbMode ? pdb.toAuth.get(S[j].pos) : S[j].pos;
      // Uncategorised sites take the accent colour in 3D, where the map's neutral grey would vanish against grey chains.
      const c = colors[mode.cat(S[j])] || pal.none;
      if (r != null) siteColor.set(r, c === pal.none ? pal.site : c);
    }
    if (pdbMode) for (const r of pdb.det.modres) if (r.chain_id === pdb.chain && !siteColor.has(+r.author_residue_number)) siteColor.set(+r.author_residue_number, '#C2410C');
    const siteAt = (a) => (inChain(a) && siteColor.has(a.resi) && baseResn(a.resn) === 'ARG' ? siteColor.get(a.resi) : null);

    // Polymer in the chosen representation and colouring; sites tint spheres and surfaces so they stay visible.
    const key = viewColorKey();
    const base = viewColorFn(key);
    const tinted = (a) => siteAt(a) || base(a);
    const rep = {
      cartoon: () => ({ cartoon: { colorfunc: base } }),
      trace: () => ({ cartoon: { style: 'trace', thickness: 0.35, colorfunc: base } }),
      sticks: () => ({ stick: { radius: 0.16, colorfunc: base } }),
      spheres: () => ({ sphere: { colorfunc: tinted } }),
    }[V.style] || (() => ({ cartoon: { colorfunc: base } }));
    viewer.setStyle({}, rep());
    if (pdbMode) {
      viewer.setStyle({ predicate: (a) => a.hetflag && !isPolymer(a) }, V.ligands ? { stick: { radius: 0.18, colorscheme: 'greenCarbon' } } : {});
      if (!V.others) viewer.setStyle({ predicate: (a) => a.chain !== pdb.chain }, {});
      // Nucleic acids: cartoon backbone with thin base sticks, purple for RNA and ochre for DNA.
      if (mol.naKinds == null) mol.naKinds = new Set(viewer.selectedAtoms({ predicate: isNuc }).map(nucKind));
      for (const kind of mol.naKinds) {
        viewer.setStyle({ predicate: (a) => isNuc(a) && nucKind(a) === kind && (V.others || a.chain === pdb.chain) },
          { cartoon: { color: NA_COLOR[kind] }, stick: { radius: 0.13, color: NA_COLOR[kind] } });
      }
    }
    viewer.setStyle({ resn: WATER }, {});

    // Methylarginine side chains (union sites in map colours; methylated residues modelled in the entry in orange).
    const byColor = new Map();
    siteColor.forEach((c, r) => { if (!byColor.has(c)) byColor.set(c, []); byColor.get(c).push(r); });
    const chainSel = pdbMode ? { chain: pdb.chain } : {};
    byColor.forEach((list, c) => {
      viewer.addStyle(Object.assign({ resi: list, atom: ['CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'] }, chainSel), { stick: { radius: 0.28, color: c } });
      viewer.addStyle(Object.assign({ resi: list, atom: 'CZ' }, chainSel), { sphere: { radius: 1.15, color: c } });
    });
    if (pdbMode) for (const r of pdb.det.modres) viewer.addStyle({ chain: r.chain_id, resi: r.author_residue_number }, { stick: { radius: 0.3, color: '#C2410C' } });

    // Side chains around the focused residue.
    if (V.contacts && mol.focus) {
      const f = mol.focus;
      const near = viewer.selectedAtoms({ within: { distance: 5, sel: { chain: f.chain, resi: f.resi } } });
      const keys = new Set(near.filter((a) => !WATER.includes(a.resn) && (V.others || !pdbMode || a.chain === pdb.chain) && (V.ligands || isPolymer(a))).map((a) => a.chain + ':' + a.resi));
      keys.delete(f.chain + ':' + f.resi);
      if (keys.size) viewer.addStyle({ predicate: (a) => keys.has(a.chain + ':' + a.resi) && (a.hetflag || !['N', 'C', 'O'].includes(a.atom)) }, { stick: { radius: 0.15, colorscheme: ELEM_GREY } });
      if (!siteColor.has(f.resi) || (pdbMode && f.chain !== pdb.chain)) viewer.addStyle({ chain: f.chain, resi: f.resi }, { stick: { radius: 0.24, colorscheme: ELEM_FOCUS } });
    }

    // Persistent site labels.
    (mol.siteLabels || []).forEach((l) => viewer.removeLabel(l));
    mol.siteLabels = [];
    if (V.labels) {
      for (const j of p.sites) {
        const r = pdbMode ? pdb.toAuth.get(S[j].pos) : S[j].pos;
        if (r == null) continue;
        const at = viewer.selectedAtoms(Object.assign({ resi: r, atom: 'CZ' }, chainSel))[0];
        if (at) mol.siteLabels.push(viewer.addLabel('' + rl(S[j].pos, P[S[j].p]), { position: { x: at.x, y: at.y, z: at.z }, backgroundColor: siteColor.get(r), backgroundOpacity: 0.85, fontColor: '#FFFFFF', fontSize: 11, inFront: true }));
      }
    }

    bindResidueClicks();
    const siteResi = Array.from(siteColor.keys());
    const hoverSel = Object.assign({ resi: siteResi, atom: 'CZ' }, chainSel);
    const unpOf = (atom) => (pdbMode ? pdb.toUnp.get(atom.resi) : atom.resi);
    viewer.setHoverable(hoverSel, true,
      (atom, v) => {
        if (atom.hoverLabel || UI.mv.labels) return;
        const pos = unpOf(atom);
        const pl = !pdbMode && af ? af.plddt.get(atom.resi) : pdbMode && mol.pdb.pred ? atom.b : null;
        const text = pos != null ? '' + rl(pos) + (pdbMode ? ' (' + atom.chain + atom.resi + ')' : '') + (pl != null ? ' · pLDDT ' + Math.round(pl) : '') : atom.resn + ' ' + atom.chain + atom.resi;
        atom.hoverLabel = v.addLabel(text, { position: { x: atom.x, y: atom.y, z: atom.z }, backgroundColor: '#12202C', backgroundOpacity: 0.85, fontColor: '#FFFFFF', fontSize: 12, inFront: true });
      },
      (atom, v) => { if (atom.hoverLabel) { v.removeLabel(atom.hoverLabel); delete atom.hoverLabel; } });
    viewer.setClickable(hoverSel, true, (atom) => {
      const pos = unpOf(atom);
      if (pos == null) { focusResidue(atom); return; }
      UI.highlight = pos;
      flashSite(pos, true);
    });

    if (reSurface !== false) updateSurface(tinted);
    viewer.render();
    const keyBox = $('#af-key');
    if (keyBox) keyBox.innerHTML = viewColorKeyHtml(key);
  }

  function updateSurface(colorfunc) {
    const { viewer } = mol;
    mol.surfGen = (mol.surfGen || 0) + 1;
    const gen = mol.surfGen;
    viewer.removeAllSurfaces();
    if (UI.mv.surface === 'off') return;
    const chain = mol.mode === 'pdb' ? mol.pdb.chain : null;
    const sel = { predicate: (a) => isPolymer(a) && (!chain || a.chain === chain) };
    const owner = mol;
    viewer.addSurface('MS', { opacity: UI.mv.surface === 'glass' ? 0.55 : 1, colorfunc }, sel, sel, undefined, (id) => {
      if (mol !== owner || owner.surfGen !== gen) { try { viewer.removeSurface(id); } catch (e) { /* viewer gone */ } }
      viewer.render();
    });
  }

  // Click any residue in the 3D view to zoom to it, label it with UniProt numbering where known and show its contacts.
  function bindResidueClicks() {
    mol.viewer.setClickable({}, true, (atom) => focusResidue(atom));
  }

  function focusResidue(atom) {
    const { viewer } = mol;
    const pdbMode = mol.mode === 'pdb';
    const pred = pdbMode && mol.pdb.pred;
    const unp = pdbMode ? (atom.chain === mol.pdb.chain ? authToUnp(atom.resi) : null) : atom.resi;
    mol.focus = { chain: atom.chain, resi: atom.resi };
    styleStructure(false);
    if (mol.label) { viewer.removeLabel(mol.label); mol.label = null; }
    const ca = viewer.selectedAtoms({ chain: atom.chain, resi: atom.resi, atom: 'CA' })[0] || atom;
    const hom = pdbMode && mol.pdb.homolog;
    mol.label = viewer.addLabel(pred ? atom.resn + ' ' + atom.resi + ' · ' + (pred.names[atom.chain] || 'chain ' + atom.chain) + ' · pLDDT ' + Math.round(atom.b) :
      atom.resn + ' ' + (pdbMode ? atom.chain : '') + atom.resi + (pdbMode && unp != null ? (hom ? ' (≙ ' + mol.p.label + ' ' + unp + ')' : ' (UniProt ' + unp + ')') : ''),
      { position: { x: ca.x, y: ca.y, z: ca.z }, backgroundColor: '#12202C', backgroundOpacity: 0.9, fontColor: '#FFFFFF', fontSize: 12, inFront: true });
    const around = [];
    for (let r = atom.resi - 6; r <= atom.resi + 6; r++) around.push(r);
    viewer.zoomTo({ chain: atom.chain, resi: around }, 600);
    viewer.render();
    const msg = $('#af-msg');
    if (msg) msg.textContent = '';
  }

  // Author residue number in the loaded PDB chain -> UniProt position, via the SIFTS segment and observed ranges.
  function authToUnp(auth) {
    if (mol.pdb.pred) return auth;
    if (mol.pdb.homolog) { const q = mol.pdb.homolog.authToQuery.get(auth); return q == null ? null : q; }
    const { e, det } = mol.pdb;
    const seg = (det.coverage || []).find(([a, b, au]) => au != null && auth >= au && auth <= au + (b - a));
    if (!seg) return null;
    const rn = seg[0] + (auth - seg[2]);
    const m = e.maps.find((x) => x.chain === det.chain && rn >= x.start && rn <= x.end);
    return m ? m.unpStart + (rn - m.start) : null;
  }

  function focusStructure(pos) {
    if (!mol || mol.p.i !== UI.detail) return;
    const { viewer } = mol;
    if (mol.label) { viewer.removeLabel(mol.label); mol.label = null; }
    const pdbMode = mol.mode === 'pdb';
    const auth = pdbMode ? mol.pdb.toAuth.get(pos) : pos;
    const chainSel = pdbMode ? { chain: mol.pdb.chain } : {};
    const atoms = auth == null ? [] : viewer.selectedAtoms(Object.assign({ resi: auth, atom: 'CZ' }, chainSel));
    const msg = $('#af-msg');
    if (msg) msg.textContent = atoms.length ? '' : '' + rl(pos) + ' is not ' + (pdbMode && mol.pdb.pred ? 'in the ' + mol.pdb.pred.label + ' model' : 'resolved in ' + (pdbMode ? mol.pdb.e.id.toUpperCase() + ' chain ' + mol.pdb.chain : 'the AlphaFold model')) + '.';
    if (!atoms.length) return;
    const a = atoms[0];
    mol.focus = { chain: a.chain, resi: a.resi };
    styleStructure(false);
    mol.label = viewer.addLabel('' + rl(pos) + (pdbMode ? ' (' + mol.pdb.chain + auth + ')' : ''), { position: { x: a.x, y: a.y, z: a.z }, backgroundColor: '#C2410C', backgroundOpacity: 0.95, fontColor: '#FFFFFF', fontSize: 13, inFront: true });
    const around = [];
    for (let r = auth - 6; r <= auth + 6; r++) around.push(r);
    viewer.zoomTo(Object.assign({ resi: around }, chainSel), 600);
    viewer.render();
  }
  // ------------------------------------------------------------- PDB structures
  const pdbCache = new Map();
  const pdbDetail = new Map();
  const METHYL_R_CODES = ['2MR', 'DA2', 'NMM', 'AGM'];
  const withTimeout = (promise, ms) => Promise.race([promise, new Promise((resolve, reject) => setTimeout(() => reject(new Error('timeout')), ms))]);

  // PDB entries for an accession (PDBe SIFTS best structures) plus the entries whose chain for this protein
  // contains a methylated arginine residue (RCSB search).
  function loadPdb(acc) {
    const base = acc.split('-')[0];
    if (pdbCache.has(base)) return pdbCache.get(base);
    const q = { query: { type: 'group', logical_operator: 'and', nodes: [
      { type: 'terminal', service: 'text', parameters: { attribute: 'rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession', operator: 'exact_match', value: base } },
      { type: 'terminal', service: 'text', parameters: { attribute: 'rcsb_polymer_entity_container_identifiers.chem_comp_monomers', operator: 'in', value: METHYL_R_CODES } },
    ] }, return_type: 'polymer_entity', request_options: { paginate: { start: 0, rows: 100 } } };
    const methyl = withTimeout(fetch('https://search.rcsb.org/rcsbsearch/v2/query?json=' + encodeURIComponent(JSON.stringify(q)))
      .then((r) => (r.status === 200 ? r.json() : { result_set: [] })), 20000)
      .then((d) => new Set((d.result_set || []).map((x) => x.identifier.split('_')[0].toLowerCase())))
      .catch(() => new Set());
    const best = withTimeout(fetchJson('https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/' + encodeURIComponent(base)), 25000)
      .then((d) => d[base] || [])
      .catch((e) => (String(e.message) === '404' ? [] : null));
    // Entries of this protein that also contain RNA or DNA polymers.
    const naSearch = (attr) => withTimeout(fetch('https://search.rcsb.org/rcsbsearch/v2/query?json=' + encodeURIComponent(JSON.stringify({ query: { type: 'group', logical_operator: 'and', nodes: [
      { type: 'terminal', service: 'text', parameters: { attribute: 'rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession', operator: 'exact_match', value: base } },
      { type: 'terminal', service: 'text', parameters: { attribute: attr, operator: 'greater_or_equal', value: 1 } },
    ] }, return_type: 'entry', request_options: { return_all_hits: true } })))
      .then((r) => (r.status === 200 ? r.json() : r.status === 204 ? { result_set: [] } : Promise.reject(new Error(r.status)))), 20000)
      .then((d) => new Set((d.result_set || []).map((x) => x.identifier.toLowerCase())));
    const na = Promise.all([naSearch('rcsb_entry_info.polymer_entity_count_RNA'), naSearch('rcsb_entry_info.polymer_entity_count_DNA')]).catch(() => null);
    const pr = Promise.all([best, methyl, na]).then(([list, meth, naSets]) => {
      const byId = new Map();
      for (const m of list || []) {
        let e = byId.get(m.pdb_id);
        if (!e) byId.set(m.pdb_id, (e = { id: m.pdb_id, method: m.experimental_method, resolution: m.resolution, chains: [], maps: [] }));
        if (!e.chains.includes(m.chain_id)) e.chains.push(m.chain_id);
        if (!e.maps.some((x) => x.chain === m.chain_id && x.unpStart === m.unp_start && x.start === m.start)) e.maps.push({ chain: m.chain_id, unpStart: m.unp_start, unpEnd: m.unp_end, start: m.start, end: m.end });
      }
      for (const id of meth) if (!byId.has(id)) byId.set(id, { id, method: '', resolution: null, chains: [], maps: [] });
      const naMap = new Map();
      if (naSets) {
        naSets[0].forEach((id) => naMap.set(id, { rna: true, dna: false }));
        naSets[1].forEach((id) => { const v = naMap.get(id) || { rna: false, dna: false }; v.dna = true; naMap.set(id, v); });
      }
      return { entries: Array.from(byId.values()), methyl: meth, failed: list === null, na: naMap, naFailed: !naSets };
    });
    pdbCache.set(base, pr);
    return pr;
  }

  // Observed residues for the mapped chain and any methylarginine residues modelled in the entry.
  function loadPdbDetail(e, withMods) {
    if (pdbDetail.has(e.id)) return pdbDetail.get(e.id);
    const chain = e.maps.length ? e.maps[0].chain : null;
    const cov = chain ? fetchJson('https://www.ebi.ac.uk/pdbe/api/pdb/entry/polymer_coverage/' + e.id + '/chain/' + encodeURIComponent(chain))
      .then((d) => {
        const out = [];
        for (const m of (d[e.id] || {}).molecules || []) for (const c of m.chains || []) for (const o of c.observed || []) out.push([o.start.residue_number, o.end.residue_number, o.start.author_residue_number]);
        return out;
      }).catch(() => null) : Promise.resolve(null);
    const mod = !withMods ? Promise.resolve([]) : fetch('https://www.ebi.ac.uk/pdbe/api/pdb/entry/modified_AA_or_NA/' + e.id)
      .then((r) => (r.ok ? r.json() : null))
      .then((d) => (d ? (d[e.id] || []).filter((r) => METHYL_R_CODES.includes(r.chem_comp_id) || /arginine/i.test(r.chem_comp_name || '')) : []))
      .catch(() => []);
    const pr = Promise.all([cov, mod]).then(([coverage, modres]) => ({ chain, coverage, modres }));
    pdbDetail.set(e.id, pr);
    return pr;
  }

  // Methylarginine sites of protein p that fall in entry e's mapped chain, as [uniprotPos, authorResidue|null, observed].
  function pdbSites(p, e, det) {
    const chain = det && det.chain;
    const maps = e.maps.filter((m) => m.chain === chain);
    const out = [];
    for (const j of p.sites) {
      const pos = S[j].pos;
      const m = maps.find((x) => pos >= x.unpStart && pos <= x.unpEnd);
      if (!m) continue;
      const rn = m.start + (pos - m.unpStart);
      const seg = det && det.coverage ? det.coverage.find(([a, b]) => rn >= a && rn <= b) : null;
      const auth = seg && seg[2] != null ? seg[2] + (rn - seg[0]) : null;
      out.push([pos, auth, !!seg]);
    }
    return out;
  }

  const methodShort = (m) => ({ 'X-ray diffraction': 'X-ray', 'Electron Microscopy': 'Cryo-EM', 'Solution NMR': 'NMR', 'Solid-state NMR': 'ssNMR', 'Electron crystallography': 'EC' }[m] || m || '—');

  function renderPdbList(p) {
    const box = $('#pdb-list');
    const sub = $('#pdb-sub');
    if (!box) return;
    loadPdb(p.acc).then((res) => {
      if (UI.detail !== p.i || !document.body.contains(box)) return;
      const inRange = (e) => e.maps.length ? p.sites.filter((j) => e.maps.some((m) => S[j].pos >= m.unpStart && S[j].pos <= m.unpEnd)).length : 0;
      const entries = res.entries.map((e) => ({ e, n: inRange(e), meth: res.methyl.has(e.id) }))
        .sort((a, b) => (b.meth - a.meth) || (b.n - a.n) || ((a.e.resolution || 99) - (b.e.resolution || 99)));
      if (!entries.length) {
        sub.textContent = res.failed ? 'The PDBe lookup did not answer in time; use PDBe-KB for the full list.' : 'No experimental structures in the PDB map to ' + p.acc + '.';
        box.innerHTML = '';
        return;
      }
      const shown = UI.pdbAll ? entries : entries.slice(0, 10);
      sub.innerHTML = fmt(entries.length) + ' PDB entr' + (entries.length === 1 ? 'y maps' : 'ies map') + ' to ' + esc(p.acc.split('-')[0]) + (res.methyl.size ? '; <b>' + res.methyl.size + '</b> model a methylated arginine in this protein' : '') +
        '. Counts are methylarginines resolved in the mapped chain (inside the construct if the coverage lookup fails).' + (res.failed ? ' The PDBe list timed out, so only methylarginine-containing entries are shown.' : '');
      box.innerHTML = '<div class="tbl-wrap"><table class="tbl compact"><thead><tr><th>PDB</th><th>Method</th><th class="r hide-xs">Resolution</th><th class="hide-xs">Chain · UniProt range</th><th>Methylarginines</th><th><span class="sr-only">View</span></th></tr></thead><tbody>' +
        shown.map(({ e, n, meth }) => '<tr data-id="' + e.id + '"><td><a href="https://www.rcsb.org/structure/' + e.id.toUpperCase() + '" target="_blank" rel="noopener"><b>' + e.id.toUpperCase() + '</b></a>' + (meth ? ' <span class="tag var" title="A methylated arginine residue is modelled in this protein\'s chain">methyl-R in model</span>' : '') +
          (res.na.get(e.id) && res.na.get(e.id).rna ? ' <span class="tag na-rna" title="The entry contains RNA">RNA</span>' : '') + (res.na.get(e.id) && res.na.get(e.id).dna ? ' <span class="tag na-dna" title="The entry contains DNA">DNA</span>' : '') + '</td>' +
          '<td>' + esc(methodShort(e.method)) + '</td><td class="r num hide-xs">' + (e.resolution ? e.resolution.toFixed(1) + ' Å' : '—') + '</td>' +
          '<td class="hide-xs">' + (e.maps.length ? esc(e.maps[0].chain) + (e.chains.length > 1 ? ' (+' + (e.chains.length - 1) + ')' : '') + ' · ' + e.maps.filter((m) => m.chain === e.maps[0].chain).map((m) => m.unpStart + '–' + m.unpEnd).join(', ') : '—') + '</td>' +
          '<td class="pdb-sites" data-n="' + n + '">' + (n ? n + ' in construct…' : '<span class="muted">none in construct</span>') + '</td>' +
          '<td><button class="btn quiet small" type="button" data-view="' + e.id + '">View</button></td></tr>').join('') +
        '</tbody></table></div>' + (entries.length > 10 ? '<button class="linkish" type="button" id="pdb-more" style="margin-top:.5rem">' + (UI.pdbAll ? 'Show the first 10' : 'Show all ' + fmt(entries.length) + ' entries') + '</button>' : '');
      const more = $('#pdb-more', box);
      if (more) more.addEventListener('click', () => { UI.pdbAll = !UI.pdbAll; renderPdbList(p); });
      shown.forEach(({ e, n, meth }) => {
        if (!n && !meth) return;
        loadPdbDetail(e, meth).then((det) => {
          const cell = $('tr[data-id="' + e.id + '"] .pdb-sites', box);
          if (!cell) return;
          const ss = pdbSites(p, e, det);
          const obs = det.coverage ? ss.filter((x) => x[2]) : ss;
          const mods = det.modres.map((r) => '<span class="tag uni" title="' + esc(r.chem_comp_name) + ' (' + esc(r.chem_comp_id) + ')">' + esc(r.chem_comp_id) + ' ' + esc(r.chain_id) + esc(r.author_residue_number) + '</span>').join(' ');
          cell.innerHTML = (obs.length ? '<b>' + obs.length + '</b> resolved: ' + obs.slice(0, 8).map((x) => '' + rl(x[0])).join(', ') + (obs.length > 8 ? ' …' : '') : '<span class="muted">' + (n ? n + ' in construct, none resolved' : '—') + '</span>') + (mods ? '<br>' + mods : '');
        });
      });
      $$('[data-view]', box).forEach((b) => b.addEventListener('click', () => {
        const e = entries.find((x) => x.e.id === b.dataset.view).e;
        showPdb(p, e, b);
      }));
    });
  }

  function showPdb(p, e, btn) {
    const holder = mol && mol.p.i === p.i ? mol.el : $('#af-viewer');
    if (!holder) return;
    if (btn) { btn.disabled = true; btn.textContent = 'Loading…'; }
    const nmr = /NMR/i.test(e.method || '');
    const file = nmr
      ? fetchText('https://files.rcsb.org/download/' + e.id.toUpperCase() + '.pdb').then((t) => [t, 'pdb']).catch(() => fetchText('https://files.rcsb.org/download/' + e.id.toUpperCase() + '.cif').then((t) => [t, 'cif']))
      : fetchText('https://files.rcsb.org/download/' + e.id.toUpperCase() + '.cif').then((t) => [t, 'cif']);
    return Promise.all([load3Dmol(), file, loadPdb(p.acc).then((res) => loadPdbDetail(e, res.methyl.has(e.id)))]).then(([lib, [cif, fmtName], det]) => {
      if (UI.detail !== p.i) return;
      if (!mol || mol.p.i !== p.i) {
        disposeStructure();
        holder.innerHTML = '';
        mol = { viewer: lib.createViewer(holder, { backgroundColor: cssVar('--surface-2') || '#eef5f7', antialias: true }), p, el: holder, label: null, spin: false, mode: 'af', afPdb: null };
      }
      const viewer = mol.viewer;
      viewer.removeAllLabels();
      mol.label = null;
      viewer.clear();
      viewer.removeAllSurfaces();
      viewer.addModel(cif, fmtName);
      Object.assign(mol, { chains: null, focus: null, siteLabels: [], naKinds: null });
      const sites = pdbSites(p, e, det).filter((x) => x[1] != null && x[2]);
      mol.mode = 'pdb';
      mol.pdb = { e, det, chain: det.chain, toAuth: new Map(sites.map((x) => [x[0], x[1]])), toUnp: new Map(sites.map((x) => [x[1], x[0]])) };
      renderViewControls();
      styleStructure();
      viewer.zoomTo();
      viewer.render();
      ['#af-reset', '#af-spin', '#af-png'].forEach((id) => { const b = $(id); if (b) b.disabled = false; });
      renderPdbSide(p);
      $('#af-viewer').scrollIntoView({ block: 'center', behavior: 'smooth' });
    }).catch(() => { if (btn) btn.textContent = 'Failed'; })
      .finally(() => { if (btn) { btn.disabled = false; if (btn.textContent === 'Loading…') btn.textContent = 'View'; } });
  }

  function renderPdbSide(p) {
    const side = $('#af-side');
    if (!side || !mol || mol.mode !== 'pdb') return;
    const { e, det, toAuth, pred } = mol.pdb;
    const pal = screenPalette();
    const mode = COLOR_MODES[UI.fig.color] || COLOR_MODES.state;
    const colors = mode.colors(pal);
    const sites = p.sites.map((j) => S[j]).filter((s) => toAuth.has(s.pos));
    const chips = () => '<div class="chips af-sites">' + sites.map((s) => '<button class="chip" type="button" data-pos="' + s.pos + '"><i class="dot" style="background:' + (colors[mode.cat(s)] || pal.none) + '"></i>' + rl(s.pos, P[s.p]) + '</button>').join('') + '</div>';
    if (pred) {
      side.innerHTML = '<p class="flabel">Showing a predicted complex</p><dl class="kv"><dt>Model</dt><dd><a href="' + esc(pred.link) + '" target="_blank" rel="noopener">' + esc(pred.label) + ' ↗</a></dd>' +
        '<dt>Source</dt><dd>' + esc(pred.method) + '</dd><dt>Scores</dt><dd>' + esc(pred.scores) + '</dd>' +
        '<dt>Chains</dt><dd>' + Object.keys(pred.names).sort().map((c) => esc(c) + ': ' + esc(pred.names[c]) + (c === pred.chain ? ' (teal)' : '')).join(' · ') + '</dd></dl>' +
        '<p class="flabel">Methylarginines in the model</p>' + (sites.length ? chips() : '<p class="note">None of this protein\'s methylarginines are in the modelled range.</p>') +
        predDistHtml() +
        '<p class="note">A predicted model, not an experimental structure: weigh each contact by the pLDDT of both residues and by the model\'s predicted aligned error' + (pred.pae ? ' (<a href="' + esc(pred.pae) + '" target="_blank" rel="noopener">PAE plot ↗</a>)' : '') + '. Colour by pLDDT under Colour.</p>' +
        (mol.afPdb ? '<p><button class="btn quiet small" type="button" id="back-af">Back to AlphaFold model</button></p>' : '');
    } else if (mol.pdb.homolog) {
      const h = mol.pdb.homolog;
      side.innerHTML = '<p class="flabel">Showing a homologous structure</p><dl class="kv"><dt>Entry</dt><dd><a href="https://www.rcsb.org/structure/' + e.id.toUpperCase() + '" target="_blank" rel="noopener">' + e.id.toUpperCase() + ' ↗</a></dd>' +
        '<dt>Method</dt><dd>' + esc(methodShort(e.method)) + (e.resolution ? ', ' + e.resolution.toFixed(1) + ' Å' : '') + '</dd>' +
        '<dt>Homolog</dt><dd>' + esc(h.desc) + (h.org ? ' <i>(' + esc(h.org) + ')</i>' : '') + (h.uniprot.length ? ' · ' + h.uniprot.map((a) => '<a href="https://www.uniprot.org/uniprotkb/' + encodeURIComponent(a) + '/entry" target="_blank" rel="noopener">' + esc(a) + '</a>').join(', ') : '') + '</dd>' +
        '<dt>Alignment</dt><dd>' + Math.round(100 * h.identity) + '% identity over ' + esc(p.label) + ' ' + h.qBeg + '–' + h.qEnd + '</dd><dt>Chain</dt><dd>' + esc(det.chain) + ' (teal; other chains grey)</dd>' +
        (h.partners.length ? '<dt>With</dt><dd>' + h.partners.map((x) => esc(x)).join('; ') + '</dd>' : '') + '</dl>' +
        '<p class="flabel">Aligned methylarginines</p>' + (sites.length ? '<div class="chips af-sites">' + sites.map((s) => '<button class="chip" type="button" data-pos="' + s.pos + '"><i class="dot" style="background:' + (colors[mode.cat(s)] || pal.none) + '"></i>' + rl(s.pos, P[s.p]) + ' ≙ ' + esc(det.chain) + toAuth.get(s.pos) + '</button>').join('') + '</div>' :
          '<p class="note">No methylarginine of ' + esc(p.label) + ' aligns to a resolved arginine in this chain.</p>') +
        (det.modres.length ? '<p class="flabel">Methylated arginines in the model</p><div class="tags">' + det.modres.map((r) => '<span class="tag var" title="' + esc(r.chem_comp_name) + '">' + esc(r.chem_comp_id) + ' ' + esc(r.chain_id) + esc(r.author_residue_number) + '</span>').join('') + '</div>' : '') +
        contactsHtml() +
        '<p class="note">A related protein, not ' + esc(p.label) + ' itself: positions come from the RCSB sequence alignment (MMseqs2), and the contacts are those of the homologous arginine.</p>' +
        (mol.afPdb ? '<p><button class="btn quiet small" type="button" id="back-af">Back to AlphaFold model</button></p>' : '');
    } else side.innerHTML = '<p class="flabel">Showing PDB entry</p><dl class="kv"><dt>Entry</dt><dd><a href="https://www.rcsb.org/structure/' + e.id.toUpperCase() + '" target="_blank" rel="noopener">' + e.id.toUpperCase() + ' ↗</a></dd>' +
      '<dt>Method</dt><dd>' + esc(methodShort(e.method)) + (e.resolution ? ', ' + e.resolution.toFixed(1) + ' Å' : '') + (/NMR/i.test(e.method || '') ? ' (first model of the ensemble)' : '') + '</dd><dt>Chain</dt><dd>' + esc(det.chain || '—') + ' (teal; other chains grey)</dd></dl>' +
      (det.modres.length ? '<p class="flabel">Methylated arginines in the model</p><div class="tags">' + det.modres.map((r) => '<span class="tag var" title="' + esc(r.chem_comp_name) + '">' + esc(r.chem_comp_id) + ' ' + esc(r.chain_id) + esc(r.author_residue_number) + '</span>').join('') + '</div><p class="note">Drawn in orange.</p>' : '') +
      '<p class="flabel">Resolved methylarginines from the union</p>' +
      (sites.length ? '<div class="chips af-sites">' + sites.map((s) => '<button class="chip" type="button" data-pos="' + s.pos + '"><i class="dot" style="background:' + (colors[mode.cat(s)] || pal.none) + '"></i>' + rl(s.pos, P[s.p]) + ' · ' + esc(det.chain) + toAuth.get(s.pos) + '</button>').join('') + '</div>' : '<p class="note">None of this protein\'s methylarginines are resolved in this chain.</p>') +
      contactsHtml() +
      (mol.afPdb ? '<p><button class="btn quiet small" type="button" id="back-af">Back to AlphaFold model</button></p>' : '');
    $$('.contact-list [data-pos]', side).forEach((b) => b.addEventListener('click', () => { UI.highlight = +b.dataset.pos; flashSite(+b.dataset.pos, false); }));
    $$('.af-sites .chip', side).forEach((b) => b.addEventListener('click', () => { UI.highlight = +b.dataset.pos; flashSite(+b.dataset.pos, false); }));
    const back = $('#back-af', side);
    if (back) back.addEventListener('click', () => {
      const v = mol.viewer;
      v.removeAllLabels(); mol.label = null; v.clear(); v.removeAllSurfaces();
      v.addModel(mol.afPdb, 'pdb');
      Object.assign(mol, { mode: 'af', pdb: null, chains: null, focus: null, siteLabels: [], naKinds: null });
      renderViewControls();
      styleStructure(); v.zoomTo(); v.render();
      renderAfSide(p, AFDATA.get(p.acc));
    });
  }



  // ------------------------------------------------------------- homologous structures
  // PDB entities related to this protein by sequence (RCSB MMseqs2 search, experimental entries), excluding entities
  // mapped to this accession. Each methylarginine is carried through the alignment to the homolog's residue.
  function loadHomologs(p) {
    const base = p.acc.split('-')[0];
    return live('homologs:' + p.acc, () => loadUniProt(p.acc).then((u) => {
      if (!u || !u.seq) return null;
      const q = { query: { type: 'terminal', service: 'sequence', parameters: { evalue_cutoff: 1e-3, identity_cutoff: 0.3, sequence_type: 'protein', value: u.seq } },
        return_type: 'polymer_entity', request_options: { paginate: { start: 0, rows: 500 }, results_content_type: ['experimental'], scoring_strategy: 'sequence', results_verbosity: 'verbose' } };
      return withTimeout(fetch('https://search.rcsb.org/rcsbsearch/v2/query', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(q) }), 45000)
        .then((r) => (r.status === 204 ? { result_set: [] } : r.ok ? r.json() : Promise.reject(new Error(r.status))));
    }).then((d) => {
      if (!d) return null;
      const sitePos = p.sites.map((j) => S[j].pos);
      const hits = [];
      for (const r of d.result_set || []) {
        const m = (((r.services || [])[0] || {}).nodes || [])[0];
        const ctx = m && (m.match_context || [])[0];
        if (!ctx || !ctx.query_aligned_seq) continue;
        // Walk the alignment: query position -> subject (entity sequence) index and residue.
        const qa = ctx.query_aligned_seq, sa = ctx.subject_aligned_seq;
        const align = new Map();
        let qi = ctx.query_beg - 1, si = ctx.subject_beg - 1;
        for (let k = 0; k < qa.length; k++) {
          if (qa[k] !== '-') qi++;
          if (sa[k] !== '-') si++;
          if (qa[k] !== '-' && sa[k] !== '-') align.set(qi, [si, sa[k]]);
        }
        const sites = sitePos.filter((pos) => align.has(pos)).map((pos) => ({ pos, sub: align.get(pos)[0], aa: align.get(pos)[1] }));
        if (!sites.length) continue;
        hits.push({ id: r.identifier, entry: r.identifier.split('_')[0], identity: ctx.sequence_identity, qBeg: ctx.query_beg, qEnd: ctx.query_end, align, sites });
      }
      const entries = Array.from(new Set(hits.map((h) => h.entry)));
      const chunks = [];
      for (let k = 0; k < entries.length; k += 80) chunks.push(entries.slice(k, k + 80));
      const gql = (ids) => fetch('https://data.rcsb.org/graphql', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ query:
        '{ entries(entry_ids: ' + JSON.stringify(ids) + ') { rcsb_id struct { title } rcsb_entry_info { resolution_combined experimental_method } polymer_entities { rcsb_id rcsb_polymer_entity { pdbx_description } ' +
        'rcsb_polymer_entity_container_identifiers { asym_ids auth_asym_ids uniprot_ids } entity_poly { rcsb_entity_polymer_type } rcsb_entity_source_organism { ncbi_scientific_name } } } }' }) })
        .then((r) => (r.ok ? r.json() : Promise.reject(new Error(r.status))));
      return Promise.all(chunks.map(gql)).then((parts) => {
        const meta = new Map();
        for (const part of parts) for (const en of ((part.data || {}).entries || [])) meta.set(en.rcsb_id, en);
        const out = [];
        for (const h of hits) {
          const en = meta.get(h.entry);
          if (!en) continue;
          const ent = (en.polymer_entities || []).find((x) => x.rcsb_id === h.id);
          if (!ent) continue;
          const ids = ent.rcsb_polymer_entity_container_identifiers || {};
          const uniprot = ids.uniprot_ids || [];
          if (uniprot.includes(base)) continue; // this protein's own entries are listed above
          const orgs = (ent.rcsb_entity_source_organism || []).map((o) => o.ncbi_scientific_name).filter(Boolean);
          const partners = (en.polymer_entities || []).filter((x) => x.rcsb_id !== h.id).map((x) => {
            const type = (x.entity_poly || {}).rcsb_entity_polymer_type || '';
            const desc = pdbDesc((x.rcsb_polymer_entity || {}).pdbx_description || type);
            return /RNA/.test(type) ? 'RNA: ' + desc : /DNA/.test(type) ? 'DNA: ' + desc : desc;
          });
          const info = en.rcsb_entry_info || {};
          out.push(Object.assign(h, {
            desc: pdbDesc((ent.rcsb_polymer_entity || {}).pdbx_description || ''), org: orgs[0] || '', uniprot, chain: (ids.auth_asym_ids || [])[0] || null, asym: (ids.asym_ids || [])[0] || null,
            partners: Array.from(new Set(partners)), title: (en.struct || {}).title || '', method: info.experimental_method || '', resolution: (info.resolution_combined || [])[0] || null,
            nR: h.sites.filter((x) => x.aa === 'R').length,
          }));
        }
        return out.filter((h) => h.chain && h.asym);
      }).then((out) => {
        // Which aligned residues are modelled: unobserved ranges of the first chain of each entity (label numbering).
        const inst = out.map((h) => h.entry + '.' + h.asym);
        const chunks2 = [];
        for (let k = 0; k < inst.length; k += 100) chunks2.push(inst.slice(k, k + 100));
        const gql2 = (ids) => fetch('https://data.rcsb.org/graphql', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ query:
          '{ polymer_entity_instances(instance_ids: ' + JSON.stringify(ids) + ') { rcsb_id rcsb_polymer_entity_instance_container_identifiers { auth_asym_id } rcsb_polymer_instance_feature { type feature_positions { beg_seq_id end_seq_id } } } }' }) })
          .then((r) => (r.ok ? r.json() : Promise.reject(new Error(r.status))));
        return Promise.all(chunks2.map(gql2)).then((parts) => {
          const gaps = new Map();
          for (const part of parts) for (const x of ((part.data || {}).polymer_entity_instances || [])) {
            const ranges = [];
            for (const f of x.rcsb_polymer_instance_feature || []) if (f.type === 'UNOBSERVED_RESIDUE_XYZ') for (const q of f.feature_positions || []) ranges.push([q.beg_seq_id, q.end_seq_id == null ? q.beg_seq_id : q.end_seq_id]);
            gaps.set(x.rcsb_id, { ranges, auth: (x.rcsb_polymer_entity_instance_container_identifiers || {}).auth_asym_id });
          }
          for (const h of out) {
            const g = gaps.get(h.entry + '.' + h.asym);
            if (g && g.auth) h.chain = g.auth;
            for (const x of h.sites) x.obs = !g || !g.ranges.some(([a, b]) => x.sub >= a && x.sub <= b);
            h.nR = h.sites.filter((x) => x.aa === 'R' && x.obs).length;
            h.nObs = h.sites.filter((x) => x.obs).length;
          }
          // Modelled, conserved arginines first; then entries with a partner; then identity.
          return out.filter((h) => h.nObs).sort((a, b) => (b.nR > 0) - (a.nR > 0) || (b.partners.length > 0) - (a.partners.length > 0) || b.nR - a.nR || b.identity - a.identity);
        });
      });
    }));
  }

  // PDB entity descriptions are often upper case ("PROTEIN (IMPORTIN BETA SUBUNIT)"); show them in sentence case.
  const pdbDesc = (s) => {
    const x = String(s || '').replace(/^PROTEIN \((.*)\)$/i, '$1');
    return x && x === x.toUpperCase() ? x.charAt(0) + x.slice(1).toLowerCase() : x;
  };
  const shorten = (s, n) => (s.length > n ? s.slice(0, n - 1) + '…' : s);

  function renderHomologList(p) {
    const box = $('#hom-list');
    const sub = $('#hom-sub');
    if (!box) return;
    loadHomologs(p).then((list) => {
      if (UI.detail !== p.i || !document.body.contains(box)) return;
      if (list == null) { sub.textContent = 'The RCSB sequence search did not answer.'; box.innerHTML = ''; return; }
      // Per methylarginine: entries where the aligned residue is a modelled arginine. Entries with a partner first,
      // then higher identity, then better resolution.
      const bySite = new Map();
      for (const h of list) for (const x of h.sites) {
        if (x.aa !== 'R' || !x.obs) continue;
        if (!bySite.has(x.pos)) bySite.set(x.pos, []);
        bySite.get(x.pos).push(h);
      }
      const rank = (a, b) => (b.partners.length > 0) - (a.partners.length > 0) || b.identity - a.identity || (a.resolution || 99) - (b.resolution || 99);
      const sites = Array.from(bySite.keys()).sort((a, b) => a - b);
      if (!sites.length) { sub.textContent = 'No experimental structure of a related protein (≥ 30% identity) models an arginine at these methylarginine positions.'; box.innerHTML = ''; return; }
      const nEntries = new Set(sites.flatMap((pos) => bySite.get(pos).map((h) => h.entry))).size;
      sub.innerHTML = fmt(sites.length) + ' of ' + fmt(p.sites.length) + ' methylarginines align to a modelled arginine in ' + fmt(nEntries) + ' PDB entr' + (nEntries === 1 ? 'y' : 'ies') + ' of related proteins (≥ 30% identity, RCSB sequence search). ' +
        'Where ' + esc(p.label) + ' has no structure of its own, the homolog\'s arginine and its contacts are the closest experimental evidence; <b>View</b> opens the entry at that residue.';
      UI.homOpen = UI.homOpen || new Set();
      const line = (h, pos) => {
        const x = h.sites.find((y) => y.pos === pos);
        return '<li title="' + esc('' + rl(pos) + ' aligns to the arginine at position ' + x.sub + ' of the ' + h.entry + ' chain ' + h.chain + ' sequence') + '"><button class="btn quiet small" type="button" data-hom="' + esc(h.id) + '" data-hpos="' + pos + '">View</button> ' +
          '<a href="https://www.rcsb.org/structure/' + h.entry + '" target="_blank" rel="noopener"><b>' + h.entry + '</b></a> <span class="muted">' + esc(methodShort(h.method)) + (h.resolution ? ' ' + (+h.resolution).toFixed(1) + ' Å' : '') + '</span> · ' +
          esc(h.desc) + (h.org ? ' <i class="muted">' + esc(h.org) + '</i>' : '') + ' · <span class="num">' + Math.round(100 * h.identity) + '%</span> identity' +
          (h.partners.length ? ' · with ' + h.partners.slice(0, 2).map((y) => '<span class="' + (/^RNA:/.test(y) ? 'tag na-rna' : /^DNA:/.test(y) ? 'tag na-dna' : 'tag') + '" title="' + esc(y) + '">' + esc(shorten(y.replace(/^(RNA|DNA): /, ''), 34)) + '</span>').join(' ') + (h.partners.length > 2 ? ' +' + (h.partners.length - 2) : '') : ' · <span class="muted">no other chain type</span>') +
          '</li>';
      };
      const shownSites = UI.homAll ? sites : sites.slice(0, 12);
      box.innerHTML = '<div class="tbl-wrap"><table class="tbl compact hom-tbl"><thead><tr><th>Site</th><th>Homologous entries with a modelled arginine here</th></tr></thead><tbody>' +
        shownSites.map((pos) => {
          const hs = bySite.get(pos).slice().sort(rank);
          const open = UI.homOpen.has(pos);
          return '<tr><td class="nw"><b>' + rl(pos) + '</b><br><span class="muted">' + fmt(hs.length) + ' entr' + (hs.length === 1 ? 'y' : 'ies') + '</span></td><td><ul class="hom-list">' +
            hs.slice(0, open ? hs.length : 3).map((h) => line(h, pos)).join('') + '</ul>' +
            (hs.length > 3 ? '<button class="linkish" type="button" data-homopen="' + pos + '">' + (open ? 'Show the top 3' : 'Show all ' + fmt(hs.length)) + '</button>' : '') + '</td></tr>';
        }).join('') + '</tbody></table></div>' +
        (sites.length > 12 ? '<button class="linkish" type="button" id="hom-more" style="margin-top:.5rem">' + (UI.homAll ? 'Show the first 12 sites' : 'Show all ' + fmt(sites.length) + ' sites') + '</button>' : '') +
        '<p class="note">Entries with another protein, RNA or DNA come first, then higher identity. Positions follow the RCSB MMseqs2 alignment; only arginines modelled in the deposited chain are listed.</p>';
      const byId = new Map(list.map((h) => [h.id, h]));
      $$('[data-hom]', box).forEach((b) => b.addEventListener('click', () => showHomolog(p, byId.get(b.dataset.hom), b, +b.dataset.hpos)));
      $$('[data-homopen]', box).forEach((b) => b.addEventListener('click', () => { const k = +b.dataset.homopen; if (UI.homOpen.has(k)) UI.homOpen.delete(k); else UI.homOpen.add(k); renderHomologList(p); }));
      const more = $('#hom-more', box);
      if (more) more.addEventListener('click', () => { UI.homAll = !UI.homAll; renderHomologList(p); });
    });
  }

  // Load a homologous PDB entry; residues map through the alignment (entity index -> author number via PDBe coverage).
  function showHomolog(p, h, btn, focusPos) {
    const holder = mol && mol.p.i === p.i ? mol.el : $('#af-viewer');
    if (!holder) return Promise.resolve();
    if (btn) { btn.disabled = true; btn.textContent = 'Loading…'; }
    const id = h.entry.toLowerCase();
    const nmr = /NMR/i.test(h.method || '');
    const file = nmr ? fetchText('https://files.rcsb.org/download/' + h.entry + '.pdb').then((t) => [t, 'pdb']).catch(() => fetchText('https://files.rcsb.org/download/' + h.entry + '.cif').then((t) => [t, 'cif']))
      : fetchText('https://files.rcsb.org/download/' + h.entry + '.cif').then((t) => [t, 'cif']);
    const cov = fetchJson('https://www.ebi.ac.uk/pdbe/api/pdb/entry/polymer_coverage/' + id + '/chain/' + encodeURIComponent(h.chain))
      .then((d) => {
        const segs = [];
        for (const m of (d[id] || {}).molecules || []) for (const c of m.chains || []) for (const o of c.observed || []) segs.push([o.start.residue_number, o.end.residue_number, o.start.author_residue_number]);
        return segs;
      }).catch(() => []);
    const mods = fetch('https://www.ebi.ac.uk/pdbe/api/pdb/entry/modified_AA_or_NA/' + id).then((r) => (r.ok ? r.json() : null))
      .then((d) => (d ? (d[id] || []).filter((r) => METHYL_R_CODES.includes(r.chem_comp_id)) : [])).catch(() => []);
    const card = $('#card-af');
    if (card) card.scrollIntoView({ block: 'start', behavior: 'smooth' });
    return Promise.all([load3Dmol(), file, cov, mods]).then(([lib, [text, fmtName], segs, modres]) => {
      if (UI.detail !== p.i) return;
      if (!mol || mol.p.i !== p.i) {
        disposeStructure();
        holder.innerHTML = '';
        mol = { viewer: lib.createViewer(holder, { backgroundColor: cssVar('--surface-2') || '#eef5f7', antialias: true }), p, el: holder, label: null, spin: false, mode: 'af', afPdb: null };
      }
      const viewer = mol.viewer;
      viewer.removeAllLabels();
      mol.label = null;
      viewer.clear();
      viewer.removeAllSurfaces();
      viewer.addModel(text, fmtName);
      Object.assign(mol, { chains: null, focus: null, siteLabels: [], naKinds: null });
      // Entity index -> author residue number from the observed segments (no mapping without coverage).
      const toAuthNum = (sub) => { const s = segs.find(([a, b]) => sub >= a && sub <= b); return s && s[2] != null ? s[2] + (sub - s[0]) : null; };
      const authToQuery = new Map();
      h.align.forEach(([sub], q) => { const a = toAuthNum(sub); if (a != null) authToQuery.set(a, q); });
      const pairs = h.sites.filter((x) => x.aa === 'R').map((x) => [x.pos, toAuthNum(x.sub)]).filter((x) => x[1] != null);
      mol.mode = 'pdb';
      mol.pdb = { e: { id, method: h.method, resolution: h.resolution ? +h.resolution : null, maps: [] }, det: { chain: h.chain, coverage: segs, modres }, chain: h.chain,
        toAuth: new Map(pairs), toUnp: new Map(pairs.map(([a, b]) => [b, a])), homolog: Object.assign({}, h, { authToQuery }) };
      renderViewControls();
      styleStructure();
      viewer.zoomTo();
      viewer.render();
      ['#af-reset', '#af-spin', '#af-png'].forEach((k) => { const b = $(k); if (b) b.disabled = false; });
      renderPdbSide(p);
      const first = focusPos != null ? focusPos : (pairs[0] || [])[0];
      if (first != null && mol.pdb.toAuth.has(first)) { UI.highlight = first; focusStructure(first); }
    }).catch(() => { if (btn) btn.textContent = 'Failed'; })
      .finally(() => { if (btn) { btn.disabled = false; if (btn.textContent === 'Loading…') btn.textContent = 'View'; } });
  }

  // ------------------------------------------------------------- predicted complexes
  // ProtVar (static export: partners whose AlphaFold-Multimer interface includes a methylarginine) and AlphaFold DB
  // complexes (live API). Both use UniProt numbering in every chain, so sites map onto the model directly.
  const PROTVAR_MODEL = (a, b) => 'https://www.ebi.ac.uk/ProtVar/api/prediction/interaction/' + encodeURIComponent(a) + '/' + encodeURIComponent(b) + '/model';

  // This protein's ProtVar partners with the methylarginines at each predicted interface, best pDockQ first.
  function predPartners(p) {
    const out = new Map();
    for (const j of p.sites) {
      const s = S[j];
      if (!s.pred) continue;
      for (const [k, sc, first, model] of s.pred.hits) {
        if (!out.has(k)) out.set(k, { k, sc, first, model, sites: [] });
        out.get(k).sites.push(s.pos);
      }
    }
    return Array.from(out.values()).sort((a, b) => b.sc - a.sc || b.sites.length - a.sites.length);
  }

  function loadAfComplexes(acc) {
    return live('afcx:' + acc, () => fetch('https://alphafold.ebi.ac.uk/api/complex/' + encodeURIComponent(acc))
      .then((r) => { if (r.status === 404) return []; if (!r.ok) throw new Error(r.status); return r.json(); })
      .then((d) => (Array.isArray(d) ? d : []).filter((x) => x.isComplex !== false && (x.uniprotAccession || []).includes(acc))
        .sort((a, b) => (b.complexPredictionAccuracy_ipTM || 0) - (a.complexPredictionAccuracy_ipTM || 0))));
  }

  function renderPredList(p) {
    const box = $('#pred-list');
    const sub = $('#pred-sub');
    if (!box) return;
    const base = p.acc.split('-')[0];
    const isoform = p.acc.includes('-');
    const pv = ixReady() && !isoform ? predPartners(p) : [];
    const geneLink = (k) => {
      const x = IX.prot[k];
      return x.p >= 0 && x.p !== p.i ? '<button class="linkish strong" type="button" data-open="' + x.p + '">' + esc(x.gene) + '</button>' : '<b>' + esc(x.gene) + '</b>';
    };
    const pvRows = pv.slice(0, UI.predAll ? pv.length : 12).map((r) => '<tr><td>' + geneLink(r.k) + ' <span class="muted mono">' + esc(IX.prot[r.k].acc) + '</span></td>' +
      '<td class="r num">' + (r.sc / 100).toFixed(2) + '</td><td>' + r.sites.slice(0, 10).map((x) => '' + rl(x)).join(', ') + (r.sites.length > 10 ? ' +' + (r.sites.length - 10) : '') + '</td>' +
      '<td>' + (r.model ? '<button class="btn quiet small" type="button" data-pv="' + r.k + '">View</button>' : '<span class="muted" title="The ProtVar API serves no model for this pair; its interface comes from an earlier model set">no model</span>') + '</td></tr>').join('');
    const pvHtml = !ixReady() ? '<p class="note">Loading interaction data…</p>' :
      isoform ? '<p class="note">Predicted complexes use canonical sequences; open the canonical entry for ' + esc(base) + '.</p>' :
      pv.length ? '<p class="flabel">ProtVar: partners with a methylarginine at the predicted interface</p><div class="tbl-wrap"><table class="tbl compact"><thead><tr><th>Partner</th><th class="r">pDockQ</th><th title="Methylarginines within 5 Å of the partner in the ProtVar 2024 interface table">Methyl-R ≤ 5 Å</th><th><span class="sr-only">View</span></th></tr></thead><tbody>' + pvRows + '</tbody></table></div>' +
        (pv.length > 12 ? '<button class="linkish" type="button" id="pred-more">' + (UI.predAll ? 'Show the top 12' : 'Show all ' + fmt(pv.length) + ' partners') + '</button>' : '') :
      '<p class="note">No ProtVar model (pDockQ ≥ 0.23) places one of these methylarginines at an interface.</p>';
    box.innerHTML = pvHtml + '<div id="afcx-list"><p class="note">Looking up AlphaFold DB complexes…</p></div>';
    if (sub) sub.innerHTML = 'Putative contacts from AlphaFold-Multimer models of protein pairs. pDockQ and ipTM score the predicted interface (pDockQ ≥ 0.23 acceptable, ≥ 0.5 confident). These are predictions of unmodified proteins; <b>View</b> opens a model above and measures each methylarginine\'s distance to the partner.';
    $$('[data-open]', box).forEach((b) => b.addEventListener('click', () => openProtein(+b.dataset.open)));
    const more = $('#pred-more', box);
    if (more) more.addEventListener('click', () => { UI.predAll = !UI.predAll; renderPredList(p); });
    $$('[data-pv]', box).forEach((b) => b.addEventListener('click', () => viewProtvar(p, pv.find((x) => x.k === +b.dataset.pv), b)));
    if (isoform) { $('#afcx-list', box).innerHTML = ''; return; }
    loadAfComplexes(base).then((list) => {
      const hold = $('#afcx-list', box);
      if (!hold || UI.detail !== p.i) return;
      if (list == null) { hold.innerHTML = '<p class="note">AlphaFold DB did not answer.</p>'; return; }
      if (!list.length) { hold.innerHTML = '<p class="note">AlphaFold DB has no predicted complex for ' + esc(base) + '.</p>'; return; }
      const partners = (x) => {
        const accs = x.uniprotAccession || [];
        const others = accs.map((acc, k) => [acc, k]).filter(([acc]) => acc !== base);
        if (!others.length) return '<b>' + esc(p.label) + '</b> <span class="tag">homodimer</span>';
        return others.map(([acc, k]) => {
          const ki = IX.byAcc.get(acc);
          const name = (x.gene || [])[k] || acc;
          const org = (x.taxId || [])[k] !== 9606 ? ' <span class="tag" title="' + esc((x.organismScientificName || [])[k] || '') + '">' + esc((x.organismScientificName || [])[k] || 'non-human') + '</span>' : '';
          return (ki != null ? geneLink(ki) : '<b>' + esc(name) + '</b>') + ' <span class="muted mono">' + esc(acc) + '</span>' + org;
        }).join(', ');
      };
      const num = (v) => (v == null ? '—' : (+v).toFixed(2));
      hold.innerHTML = '<p class="flabel" style="margin-top:.8rem">AlphaFold DB complexes</p><div class="tbl-wrap"><table class="tbl compact"><thead><tr><th>Partner</th><th class="r">ipTM</th><th class="r">pDockQ</th><th class="r hide-xs">Mean pLDDT</th><th><span class="sr-only">View</span></th></tr></thead><tbody>' +
        list.map((x) => '<tr><td>' + partners(x) + '</td><td class="r num">' + num(x.complexPredictionAccuracy_ipTM) + '</td><td class="r num">' + num(x.complexPredictionAccuracy_pDockQ) + '</td><td class="r num hide-xs">' + (x.globalMetricValue != null ? (+x.globalMetricValue).toFixed(1) : '—') + '</td>' +
          '<td><button class="btn quiet small" type="button" data-afcx="' + esc(x.modelEntityId) + '">View</button></td></tr>').join('') + '</tbody></table></div>' +
        '<p class="note">AlphaFold DB predicted complexes (CC BY 4.0). EMBL-EBI guidance: ipTM above 0.8 is a confident prediction, below 0.6 likely a failed one, and 0.6–0.8 a grey zone; low pLDDT in disordered regions can lower ipTM even when the complex is right (<a href="https://www.ebi.ac.uk/training/online/courses/alphafold/inputs-and-outputs/evaluating-alphafolds-predicted-structures-using-confidence-scores/confidence-scores-in-alphafold-multimer/" target="_blank" rel="noopener">source ↗</a>).</p>';
      $$('[data-open]', hold).forEach((b) => b.addEventListener('click', () => openProtein(+b.dataset.open)));
      $$('[data-afcx]', hold).forEach((b) => b.addEventListener('click', () => {
        const x = list.find((y) => y.modelEntityId === b.dataset.afcx);
        viewAfdb(p, x.modelEntityId, 'ipTM ' + num(x.complexPredictionAccuracy_ipTM) + ' · pDockQ ' + num(x.complexPredictionAccuracy_pDockQ) + ' · mean pLDDT ' + (x.globalMetricValue != null ? (+x.globalMetricValue).toFixed(1) : '—'), b, null);
      }));
    });
  }

  // ProtVar model of this protein with partner r.k. The API serves a later model run than the 2024 interface table,
  // sometimes with the accessions the other way round; its first accession is chain A.
  function viewProtvar(p, r, btn) {
    if (!r || !r.model) return Promise.resolve();
    const base = p.acc.split('-')[0];
    const other = IX.prot[r.k];
    const me = IX.prot[p.i];
    const first = r.model === 2 ? !r.first : !!r.first;
    const [a, b] = first ? [base, other.acc] : [other.acc, base];
    const url = PROTVAR_MODEL(a, b);
    if (UI.startStructure) UI.startStructure();
    return showModel(p, () => Promise.resolve({
      id: 'protvar_' + a + '_' + b, url, chain: first ? 'A' : 'B', label: 'ProtVar ' + me.gene + '–' + other.gene,
      method: 'AlphaFold-Multimer via ProtVar (EMBL-EBI)', scores: r.scoreLabel || 'pDockQ ' + (r.sc / 100).toFixed(2) + ' (2024 interface table)', link: url, pae: null,
      names: first ? { A: me.gene, B: other.gene } : { A: other.gene, B: me.gene }, partnerLabel: other.gene, tableSites: r.noTable ? null : r.sites,
    }), btn, r.sites[0]);
  }

  // AlphaFold DB complex model: chain ids and file URL from the prediction API.
  function viewAfdb(p, id, scores, btn, focusPos) {
    const base = p.acc.split('-')[0];
    if (UI.startStructure) UI.startStructure();
    return showModel(p, () => fetchJson('https://alphafold.ebi.ac.uk/api/prediction/' + encodeURIComponent(id)).then((chains) => {
      const mine = chains.find((c) => c.uniprotAccession === base) || chains[0];
      const names = {};
      chains.forEach((c) => { names[c.chainId] = c.gene || c.uniprotAccession; });
      const partnerLabel = chains.filter((c) => c.chainId !== mine.chainId).map((c) => c.gene || c.uniprotAccession).join(', ') || p.label;
      return {
        id, url: mine.pdbUrl, chain: mine.chainId, label: id, names, partnerLabel, scores,
        method: (mine.toolUsed || 'AlphaFold-Multimer') + ' via AlphaFold DB', link: 'https://alphafold.ebi.ac.uk/entry/' + encodeURIComponent(id), pae: mine.paeImageUrl || null,
      };
    }), btn, focusPos);
  }

  // Load a predicted complex into the 3D panel. spec() resolves to {url, chain, names, ...}; numbering is UniProt.
  function showModel(p, spec, btn, focusPos) {
    const holder = mol && mol.p.i === p.i ? mol.el : $('#af-viewer');
    if (!holder) return Promise.resolve();
    if (btn) { btn.disabled = true; btn.textContent = 'Loading…'; }
    const card = $('#card-af');
    if (card) card.scrollIntoView({ block: 'start', behavior: 'smooth' });
    return Promise.all([load3Dmol(), spec()]).then(([lib, m]) => fetchText(m.url).then((text) => [lib, m, text])).then(([lib, m, text]) => {
      if (UI.detail !== p.i) return;
      if (!/^ATOM/m.test(text)) throw new Error('empty model');
      if (!mol || mol.p.i !== p.i) {
        disposeStructure();
        holder.innerHTML = '';
        mol = { viewer: lib.createViewer(holder, { backgroundColor: cssVar('--surface-2') || '#eef5f7', antialias: true }), p, el: holder, label: null, spin: false, mode: 'af', afPdb: null };
      }
      const viewer = mol.viewer;
      viewer.removeAllLabels();
      mol.label = null;
      viewer.clear();
      viewer.removeAllSurfaces();
      viewer.addModel(text, 'pdb');
      Object.assign(mol, { chains: null, focus: null, siteLabels: [], naKinds: null });
      const present = new Set(viewer.selectedAtoms({ chain: m.chain, atom: 'CA' }).map((a) => a.resi));
      const sites = p.sites.map((j) => S[j].pos).filter((pos) => present.has(pos));
      mol.mode = 'pdb';
      mol.pdb = { e: { id: m.id, method: m.method, resolution: null, maps: [] }, det: { chain: m.chain, coverage: null, modres: [] }, chain: m.chain,
        toAuth: new Map(sites.map((x) => [x, x])), toUnp: new Map(sites.map((x) => [x, x])), pred: m };
      renderViewControls();
      styleStructure();
      viewer.zoomTo();
      viewer.render();
      ['#af-reset', '#af-spin', '#af-png'].forEach((id) => { const b = $(id); if (b) b.disabled = false; });
      renderPdbSide(p);
      if (focusPos != null && present.has(focusPos)) { UI.highlight = focusPos; focusStructure(focusPos); }
    }).catch(() => { if (btn) btn.textContent = 'Failed'; })
      .finally(() => { if (btn) { btn.disabled = false; if (btn.textContent === 'Loading…') btn.textContent = 'View'; } });
  }

  // ------------------------------------------------------------- interactions
  // Static context from scripts/37_export_interaction_data.py (STRING physical links, Complex Portal, GO,
  // PDBe-KB interfaces) plus live STRING, QuickGO, ENCODE and ChIP-Atlas calls on protein pages.
  const IXF = { prmt: 1, tudor: 2, rna: 4, dna: 8, rnaExp: 16, dnaExp: 32, wd40: 64, phd: 128 };
  // Partner classes named in the filters; "any other protein" is a partner in none of them.
  const IXF_CLASS = IXF.prmt | IXF.tudor | IXF.wd40 | IXF.phd;
  const IFACE_OPTS = [[4, 'RNA'], [8, 'DNA'], [1, 'Another protein'], [2, 'A copy of itself', 'Self-contact: the arginine touches another copy of the same protein in a homodimer or larger homo-oligomer']];
  const IFACE_NAME = { 4: 'RNA', 8: 'DNA', 1: 'another protein', 2: 'a copy of itself (homo-oligomer)' };
  const IFACE_COLOR = { 4: '#7B4FBF', 8: '#B7791F', 1: '#0E7490', 2: '#8A99A3' };
  const NA_COLOR = { RNA: '#7B4FBF', DNA: '#B7791F' };
  const BIND_OPTS = [
    [1, 'a PRMT', 'PRMT1, PRMT2, PRMT3, CARM1 (PRMT4), PRMT5, PRMT6, PRMT7, PRMT8 or PRMT9'],
    [2, 'a Tudor-domain protein', 'InterPro IPR002999: methylarginine readers such as SMN1/2, SMNDC1, TDRD3 and SND1 (TDRKH binds methylated PIWI proteins but can also bind PIWIL1 independently of methylation), and methyllysine readers such as PHF1, SETDB1 and KDM4A'],
    [4, 'a WD40-repeat protein', 'InterPro IPR001680, e.g. WDR77 (MEP50), the PRMT5 cofactor, and WDR5, which reads symmetric dimethyl-H3R2'],
    [5, 'a PHD-finger protein', 'InterPro IPR001965 or IPR019787: zinc fingers, many of which read the histone H3 N-terminus, where H3R2 methylation can block binding (e.g. UHRF1)'],
    [3, 'another methylated protein', 'A partner that is itself in this methylarginine set'],
    [6, 'any other protein', 'A partner that is not a PRMT and has no Tudor, WD40 or PHD domain'],
  ];
  function bindMatch(v, k) {
    const fl = IX.prot[k].flags;
    switch (v) {
      case 1: return !!(fl & IXF.prmt);
      case 2: return !!(fl & IXF.tudor);
      case 3: return k < IX.nP;
      case 4: return !!(fl & IXF.wd40);
      case 5: return !!(fl & IXF.phd);
      case 6: return !(fl & IXF_CLASS);
      default: return false;
    }
  }
  const IX_SCORES = [[400, 'Medium (≥ 0.4, live from STRING)'], [700, 'High (≥ 0.7)'], [900, 'Highest (≥ 0.9)']];
  const RAIL_SCORES = [[700, 'High (≥ 0.7)'], [900, 'Highest (≥ 0.9)']];
  const PRED_OPTS = [[0, 'Any'], [23, 'pDockQ ≥ 0.23 (acceptable)'], [50, 'pDockQ ≥ 0.5 (confident)']];
  const STRING_CALLER = 'shechterlab_methylarginine_explorer';
  const IX_WD_COLOR = '#0E7490', IX_PHD_COLOR = '#7B4FBF';
  let IX = null;
  let ixPromise = null;

  function loadInteractions() {
    if (!ixPromise) {
      ixPromise = fetch('data/methylarg_interactions.json')
        .then((r) => { if (!r.ok) throw new Error(r.status); return r.json(); })
        .then(initInteractions)
        .catch(() => { IX = { failed: true }; return null; })
        .then((ix) => { afterInteractionsLoaded(); return ix; });
    }
    return ixPromise;
  }

  function initInteractions(d) {
    const nP = P.length;
    const prot = d.proteins.map((r, k) => ({ k, acc: r[0], gene: (r[1] || r[0]).replace(/; /g, '/'), flags: r[2], ensp: r[3], p: k < nP ? k : -1 }));
    const byEnsp = new Map();
    prot.forEach((x) => { if (x.ensp && !byEnsp.has(x.ensp)) byEnsp.set(x.ensp, x.k); });
    const geneOf = (k) => (prot[k] ? prot[k].gene : '?');
    const cx = d.complexes.map((r, k) => {
      const c = { k, id: r[0], name: r[1], pred: !!r[2], members: r[3], n: r[4], m: r[5], fold: r[6], pv: r[7], q: r[8] };
      if (!c.name || c.name === '-') {
        // hu.MAP predictions have no names: label them by their members, methylated ones first.
        const ids = c.members.map((m) => (Array.isArray(m) ? m[0] : m)).sort((a, b) => (b < nP) - (a < nP));
        c.name = ids.slice(0, 4).map(geneOf).join(', ') + (ids.length > 4 ? ' +' + (ids.length - 4) : '');
      }
      return c;
    });
    const cxOf = new Map();
    for (const c of cx) {
      for (const slot of c.members) for (const id of (Array.isArray(slot) ? slot : [slot])) {
        if (!cxOf.has(id)) cxOf.set(id, []);
        cxOf.get(id).push(c.k);
      }
    }
    for (const p of P) {
      const flat = d.string[p.i] || [];
      p.partners = [];
      let last = 0;
      for (let k = 0; k < flat.length; k += 2) { last += flat[k]; p.partners.push([last, flat[k + 1]]); }
      p.partners.sort((a, b) => b[1] - a[1]);
      p.ixFlags = prot[p.i] ? prot[p.i].flags : 0;
      p.cx = cxOf.get(p.i) || [];
      // Summaries used to colour this protein when it appears as someone's partner.
      let has1 = false, has2 = false, gsk = false, ms = false, measured = false;
      for (const j of p.sites) {
        const s = S[j];
        if (s.st === 1 || s.st === 3) has1 = true;
        if (s.st === 2 || s.st === 3) has2 = true;
        if (s.trf & 1) gsk = true;
        if (s.trf & 4) ms = true;
        if (s.tr.length) measured = true;
      }
      p.stCat = has1 && has2 ? 2 : has2 ? 1 : has1 ? 0 : 3;
      p.trCat = gsk && ms ? 2 : gsk ? 0 : ms ? 1 : measured ? 3 : 4;
    }
    for (const [j, mask, partners, pdbs, nPdb] of d.iface) {
      const s = S[j];
      if (!s) continue;
      s.iface = mask; s.ifP = partners; s.ifPdb = pdbs; s.ifN = nPdb;
    }
    // CD-CODE: this protein's own confidence, role and evidence in each condensate.
    const cd = new Map();
    for (const [i, ci, score, roles, evid, nPub, pmids] of d.cd_members || []) {
      if (!cd.has(i)) cd.set(i, []);
      cd.get(i).push({ ci, score, roles, evid, nPub, pmids });
    }
    for (const [j, rows, modelled] of d.contacts || []) {
      const s = S[j];
      if (!s) continue;
      s.contacts = rows;
      s.modelledMe = modelled;
    }
    const wr = new Map((d.wr || []).map(([k, list]) => [k, new Set(list)]));
    for (const g of d.geometry || []) {
      const s = S[g[0]];
      if (s) s.geo = { pdb: g[1], chain: g[2], levy: g[3], rMono: g[4], rCplx: g[5], na: g[6], baseHb: g[7] || [], dNa: g[8], dProt: g[9], protHb: g[10] };
    }
    for (const [j, pdbs, comps] of d.modelled || []) { if (S[j]) S[j].modelled = { pdbs, comps }; }
    // ProtVar predicted heterodimer interfaces: [partner, pDockQ x 100, 1 if this protein is ProtVar's first chain].
    for (const [j, hits, n] of d.pred || []) { if (S[j]) S[j].pred = { hits, n }; }
    // Methylarginines measured in served ProtVar models and AlphaFold DB complexes (<= 5 Å of the partner).
    const gallery = (d.gallery || []).map(([j, k, src, score, pdq, plR, plP, d10, flags, res, ref, chain]) => ({ j, k, src, score, pdq, plR, plP, d: d10 / 10, flags, res, ref, chain }));
    for (const g of gallery) { const s = S[g.j]; if (s) (s.gal = s.gal || []).push(g); }
    // Residue-level RNA/DNA crosslinks (pepR-MS compendium, RBDmap, Trendel 2025, UniProt DNA-binding regions).
    for (const [j, r, n, dna, f, one] of d.xl || []) { if (S[j]) S[j].xl = { r, n, d: dna, f, one: !!one }; }
    const xlProt = new Map();
    for (const [i, rnaPos, rnaBits, dnaPos, dnaPeps, dnaBind, dnaInt, rbd, dnaCand] of d.xl_prot || []) xlProt.set(i, { rnaPos, rnaBits, dnaPos, dnaPeps, dnaBind, dnaInt: !!dnaInt, rbd: rbd || [], dnaCand: dnaCand || [] });
    const modelledGaps = new Map((d.modelled_gaps || []).map(([k, pos, pdbs]) => [k, { pos, pdbs }]));
    IX = { meta: d.meta, prot, byEnsp, cx, cxOf, enrich: d.enrich || {}, goEnrich: d.go_enrich || {}, strata: d.enrich_strata || {}, cd, wr, nP,
      geoSummary: d.geometry_summary || {}, rg: d.rg_iface || {}, modelledGaps, condConf: d.cond_conf || {}, nModelled: (d.modelled || []).length,
      predEnrich: d.pred_enrich || {}, nPred: (d.pred || []).length, condHist: d.cond_score_hist || {}, gallery,
      xlProt, xlEnrich: d.xl_enrich || {}, xlSets: d.xl_sets || [], xlTypes: d.xl_types || {}, xlFreq: d.xl_freq || null };
    IX.byAcc = new Map(prot.map((x) => [x.acc, x.k]));
    IX.prmts = prot.filter((x) => x.flags & IXF.prmt).sort((a, b) => a.gene.localeCompare(b.gene, 'en', { numeric: true }));
    IX.tudors = prot.filter((x) => x.flags & IXF.tudor).sort((a, b) => a.gene.localeCompare(b.gene));
    IX.wd40s = prot.filter((x) => x.flags & IXF.wd40).sort((a, b) => a.gene.localeCompare(b.gene));
    IX.phds = prot.filter((x) => x.flags & IXF.phd).sort((a, b) => a.gene.localeCompare(b.gene));
    return IX;
  }

  function afterInteractionsLoaded() {
    for (const s of S) s.ev = evMask(s);
    for (const p of P) p.pev = pevMask(p);
    // Partner and complex terms resolve against the interaction data: re-read the search now that it is here.
    const lx = UI.legacyIx, extra = [];
    UI.legacyIx = null;
    if (lx && ixReady()) {
      if (lx.partner != null && IX.prot[+lx.partner]) extra.push('binds:' + quoteTerm(IX.prot[+lx.partner].gene));
      if (lx.cpx != null && IX.cx[+lx.cpx]) extra.push('complex:' + IX.cx[+lx.cpx].id);
    }
    if (extra.length) setQuery([Q.raw].concat(extra).join(' ').trim());
    else if (Q.raw) Q = parseQuery(Q.raw);
    renderRail();
    compute();
    if (['overview', 'proteins', 'sites', 'interactions', 'figures', 'report', 'compare'].includes(UI.view)) renderView();
    else if (UI.view === 'protein' && UI.detail != null) {
      const p = P[UI.detail];
      renderIxCard(p);
      renderCondCard(p);
      renderSiteTable(p);
      renderContactCard(p);
      if (UI.highlight != null) renderSiteCard(UI.highlight);
      drawMap(p, UNIDATA.get(p.acc));
      const pl = $('#pred-list');
      if (pl && pl.childElementCount) renderPredList(p);
    }
  }

  const ixReady = () => IX && !IX.failed;
  const ixGene = (k) => (k === -1 ? 'RNA' : k === -2 ? 'DNA' : IX.prot[k] ? IX.prot[k].gene : '?');
  const partnersAt = (p, min) => (p.partners || []).filter(([, sc]) => sc >= min);

  function ifaceTags(s, max) {
    if (!s.iface) return '';
    const title = esc('PDB entries: ' + (s.ifPdb || []).map((x) => x.toUpperCase()).join(', ') + (s.ifN > (s.ifPdb || []).length ? ' and ' + (s.ifN - s.ifPdb.length) + ' more' : ''));
    const tags = [];
    if (s.iface & 4) tags.push('<span class="tag na-rna" title="' + title + '">RNA</span>');
    if (s.iface & 8) tags.push('<span class="tag na-dna" title="' + title + '">DNA</span>');
    const prots = (s.ifP || []).filter((k) => k >= 0 && k !== P[s.p].i);
    prots.slice(0, max || 3).forEach((k) => tags.push('<span class="tag iface-p" title="' + title + '">' + esc(ixGene(k)) + '</span>'));
    if (prots.length > (max || 3)) tags.push('<span class="tag" title="' + esc(prots.slice(max || 3).map(ixGene).join(', ')) + '">+' + (prots.length - (max || 3)) + '</span>');
    if (s.iface & 2) tags.push('<span class="tag" title="' + esc('Self-contact: touches another copy of this same protein (homo-oligomer). ') + title + '">self-contact</span>');
    return tags.join(' ');
  }

  function ifaceText(s) {
    if (!s.iface) return '';
    return (s.ifP || []).map((k) => (k === P[s.p].i ? 'self (homo-oligomer)' : ixGene(k))).join(';');
  }

  // Residues each methylarginine touches in PDB structures (precomputed side-chain contacts, <= 4 Å).
  const CONTACT_KIND = { '-1': 'RNA', '-2': 'DNA', '-3': 'Ligand', '-4': 'Protein' };
  function contactGroups(s) {
    if (!s.contacts || !s.contacts.length) return [];
    const groups = new Map();
    for (const [k, pdb, labels] of s.contacts) {
      if (!groups.has(k)) groups.set(k, { k, pdbs: [], labels: [] });
      const g = groups.get(k);
      if (!g.pdbs.includes(pdb)) g.pdbs.push(pdb);
      labels.forEach((l) => { if (!g.labels.includes(l)) g.labels.push(l); });
    }
    return Array.from(groups.values()).map((g) => ({
      ...g, name: g.k >= 0 ? (g.k === P[s.p].i ? P[s.p].label + ' (other copy)' : ixGene(g.k)) : CONTACT_KIND[g.k],
      cls: g.k === -1 ? 'na-rna' : g.k === -2 ? 'na-dna' : g.k >= 0 ? 'iface-p' : '',
    })).sort((x, y) => ((x.k === -1 || x.k === -2) ? 0 : 1) - ((y.k === -1 || y.k === -2) ? 0 : 1));
  }
  const LEVY_LABEL = ['not buried in this entry', 'interface support', 'interface rim', 'interface core'];
  const pdbButtons = (s, pdbs) => pdbs.slice(0, 3).map((x) => '<button class="linkish pdb-open" type="button" data-cpdb="' + esc(x) + '" data-cpos="' + s.pos + '" data-cp="' + s.p + '" title="Open ' + x.toUpperCase() + ' in the 3D viewer at ' + rl(s.pos, P[s.p]) + '">' + x.toUpperCase() + '</button>').join(', ') +
    (pdbs.length > 3 ? ' <span title="' + esc(pdbs.map((x) => x.toUpperCase()).join(', ')) + '">+' + (pdbs.length - 3) + '</span>' : '');
  function geoLine(s) {
    const g = s.geo;
    if (!g) return '';
    const na = g.na === 'base' ? 'H-bonds to ' + (g.baseHb.length ? 'bases ' + g.baseHb.slice(0, 4).join(' ') : 'nucleobases') : g.na === 'backbone' ? 'H-bonds to the nucleic-acid backbone' : g.na === 'contact' ? 'touches nucleic acid without an H-bond' : '';
    return '<span class="geo" title="' + esc('Relative accessibility of the arginine alone ' + Math.round(100 * g.rMono) + '%, in the complex ' + Math.round(100 * g.rCplx) + '% (' + g.pdb.toUpperCase() + ' chain ' + g.chain + '; Levy 2010 classes). ' + (g.protHb ? g.protHb + ' H-bond(s) to protein.' : '')) + '">' +
      esc(LEVY_LABEL[g.levy]) + ' <span class="muted">' + g.pdb.toUpperCase() + '</span>' + (na ? ' · ' + esc(na) : '') + '</span>';
  }
  // ProtVar predicted interface: the best-scoring partner, the rest in the tooltip.
  const predTag = (s) => {
    const best = (s.gal || []).filter((g) => galPass(g, GAL_DEFAULT)).sort(galRank)[0];
    if (!s.pred && !best) return '';
    const h = s.pred ? s.pred.hits : [];
    const tip = (h.length ? 'Predicted heterodimer interface (ProtVar, AlphaFold-Multimer, ≤ 5 Å): ' + h.map(([k, sc]) => ixGene(k) + ' pDockQ ' + (sc / 100).toFixed(2)).join(', ') + (s.pred.n > h.length ? ' and ' + (s.pred.n - h.length) + ' more' : '') + '. ' : '') +
      (best ? 'Best-supported predicted contact: ' + galBondText(best) + ' (' + galModelText(best) + ', pLDDT ' + best.plR + '/' + best.plP + ').' : '');
    const name = best ? ixGene(best.k) : ixGene(h[0][0]);
    const more = new Set(h.map((x) => x[0]).concat((s.gal || []).map((g) => g.k)));
    return '<span class="tag pred' + (best ? ' pred-hi' : '') + '" title="' + esc(tip) + '">pred. ' + esc(name) + (more.size > 1 ? ' +' + (more.size - 1) : '') + '</span>';
  };
  // RNA/DNA crosslinking mass spectrometry in living human cells.
  const XL_OPTS = [[1, 'RNA-crosslinked residue'], [2, 'Within 3 aa of an RNA crosslink'], [4, 'RBDmap RNA-bound fragment'],
    [8, 'DNA crosslink at this residue', 'The DNA crosslink is localised to this residue alone (Trendel et al. 2025)'],
    [32, 'DNA crosslink candidate', 'One of the residues a DNA crosslink may sit on; most DNA crosslinks are not localised to a single residue'], [16, 'UniProt DNA-binding region']];
  const xlSets = (bits) => (IX.xlSets || []).filter((_, k) => bits & (1 << k)).concat(bits & 32 ? ['pepR-MS focused site'] : []);
  function xlMask(s) {
    const x = s.xl;
    if (!x) return 0;
    return (x.r ? 1 : 0) | (x.n ? 2 : 0) | (x.f ? 4 : 0) | (x.d & 1 ? 8 : 0) | (x.d & 4 ? 16 : 0) | (x.d & 9 ? 32 : 0);
  }
  // Interaction evidence per site as bits, for the evidence map and its filter. No bit means nothing recorded: most
  // proteins were never crystallised or crosslinked, so it is not evidence of no contact.
  const EV_SETS = [
    [1, 'RNA crosslink', 'RNA XL', 'Crosslinked to RNA at this residue in living cells (pepR-MS compendium)', (s) => xlMask(s) & 1],
    [2, 'RNA-bound fragment', 'RBDmap', 'Inside an RBDmap RNA-bound fragment (Castello et al. 2016)', (s) => xlMask(s) & 4],
    [4, 'DNA crosslink', 'DNA XL', 'A DNA crosslink localised to, or with a candidate at, this residue (Trendel et al. 2025)', (s) => xlMask(s) & 40],
    [8, 'PDB protein interface', 'PDB protein', 'Contacts another protein, or another copy of itself, in a PDB structure (PDBe-KB)', (s) => s.iface & 3],
    [16, 'PDB RNA/DNA interface', 'PDB RNA/DNA', 'Contacts RNA or DNA in a PDB structure (PDBe-KB)', (s) => s.iface & 12],
    [32, 'Predicted interface', 'Predicted', 'At the partner interface of an AlphaFold-Multimer model (ProtVar, pDockQ ≥ 0.23)', (s) => s.pred && s.pred.hits[0][1] >= 23],
    [64, 'Confident predicted contact', 'Best contact', 'Salt bridge or hydrogen bond to a partner in a confident predicted complex, both residues at pLDDT ≥ 80', (s) => (s.gal || []).some((g) => galPass(g, GAL_DEFAULT))],
  ];
  const evMask = (s) => EV_SETS.reduce((m, [b, , , , t]) => (t(s) ? m | b : m), 0);
  // Experimental evidence (crosslinks, RBDmap, PDB interfaces) and predicted evidence (predicted-complex interfaces and
  // confident predicted contacts), the two levels of a site's interaction status.
  const EV_EXP = 1 | 2 | 4 | 8 | 16, EV_PRED = 32 | 64;
  // A site's interaction status: 'exp' (any experimental evidence), 'pred' (predicted only) or 'none'; null until loaded.
  const ixStatus = (s) => (!ixReady() ? null : s.ev & EV_EXP ? 'exp' : s.ev & EV_PRED ? 'pred' : 'none');
  const IXS = { exp: ['Experimental', 'experimental'], pred: ['Predicted only', 'predicted-only'], none: ['None recorded', 'none'] };
  const ixRank = (s) => { const st = ixStatus(s); return st == null ? null : (st === 'exp' ? 20 : st === 'pred' ? 10 : 0) + popcount(s.ev || 0); };
  // Partners named by the evidence at a site: RNA/DNA, PDB interface partners and predicted-complex partners.
  function ixPartners(s) {
    const exp = new Set(), pred = new Set();
    if (s.xl && (s.xl.r || s.xl.f)) exp.add('RNA');
    if (s.xl && s.xl.d & 9) exp.add('DNA');
    if (s.iface & 4) exp.add('RNA');
    if (s.iface & 8) exp.add('DNA');
    (s.ifP || []).filter((k) => k >= 0 && k !== s.p).forEach((k) => exp.add(ixGene(k)));
    if (s.iface & 2) exp.add(P[s.p].label + ' (self)');
    if (s.pred && s.pred.hits[0][1] >= 23) s.pred.hits.filter((h) => h[1] >= 23).forEach(([k]) => pred.add(ixGene(k)));
    (s.gal || []).filter((g) => galPass(g, GAL_DEFAULT)).forEach((g) => pred.add(ixGene(g.k)));
    return { exp: Array.from(exp), pred: Array.from(pred).filter((g) => !exp.has(g)) };
  }
  function ixBadge(s, withPartners) {
    const st = ixStatus(s);
    if (!st) return '<span class="ixs ixs-wait">…</span>';
    const m = s.ev || 0, pt = ixPartners(s);
    const tip = st === 'none' ? 'No contact evidence recorded at this residue. Most proteins were never crystallised or crosslinked, so this is not evidence of no contact.'
      : EV_SETS.filter(([b]) => m & b).map(([, l]) => l).join(', ') + (pt.exp.length ? '. Experimental partners: ' + pt.exp.join(', ') : '') + (pt.pred.length ? '. Predicted partners: ' + pt.pred.join(', ') : '');
    const label = st === 'exp' ? 'Experimental' + (m & EV_PRED ? ' + pred.' : '') : st === 'pred' ? 'Predicted' : 'None';
    const who = withPartners && st !== 'none' ? '<span class="ixs-who">' + esc(pt.exp.concat(pt.pred).slice(0, 3).join(', ') + (pt.exp.length + pt.pred.length > 3 ? ' +' + (pt.exp.length + pt.pred.length - 3) : '')) + '</span>' : '';
    return '<span class="ixs ixs-' + st + '" title="' + esc(tip) + '"><i></i>' + label + '</span>' + who;
  }
  // Table cell: the status, then the evidence tags and residue contacts.
  function ixCell(s, max) {
    if (!ixReady()) return '<span class="muted">' + (IX ? '—' : '…') + '</span>';
    const tags = ifaceTags(s, max) + modelledTag(s) + predTag(s) + xlTags(s);
    return ixBadge(s) + (tags ? '<div class="tags ixs-tags">' + tags + '</div>' : '') + contactHtml(s, Math.min(2, max));
  }
  const evLabel = (mask) => EV_SETS.filter(([b]) => mask & b).map(([, l]) => l).join(' + ');
  // Protein-level evidence, for the protein view of the evidence map and its filter (STRING physical links ≥ 0.7).
  const PEV_SETS = [
    [1, 'GO RNA binding', 'RNA-binding', 'Gene Ontology RNA binding (GO:0003723 or a child term)', (p) => p.ixFlags & IXF.rna],
    [2, 'GO DNA binding', 'DNA-binding', 'Gene Ontology DNA binding (GO:0003677 or a child term)', (p) => p.ixFlags & IXF.dna],
    [4, 'In a condensate', 'Condensate', 'Member of a CD-CODE condensate, any confidence', (p) => p.cond.length > 0],
    [8, 'Binds a PRMT', 'PRMT partner', 'STRING physical link ≥ 0.7 to a PRMT', (p) => partnersAt(p, 700).some(([k]) => bindMatch(1, k))],
    [16, 'Binds a Tudor protein', 'Tudor partner', 'STRING physical link ≥ 0.7 to a Tudor-domain protein (methylarginine readers such as SMN1, SMNDC1, TDRD3 and SND1)', (p) => partnersAt(p, 700).some(([k]) => bindMatch(2, k))],
    [32, 'In a curated complex', 'Complex', 'Member of a Complex Portal curated complex', (p) => p.cx.some((k) => !IX.cx[k].pred)],
    [64, 'Site contact evidence', 'Site contact', 'At least one methylarginine with site-level contact evidence (the site view of this map)', (p) => p.sites.some((j) => S[j].ev)],
  ];
  const pevMask = (p) => PEV_SETS.reduce((m, [b, , , , t]) => (t(p) ? m | b : m), 0);
  const pevLabel = (mask) => PEV_SETS.filter(([b]) => mask & b).map(([, l]) => l).join(' + ');

  function xlTags(s) {
    const x = s.xl;
    if (!x) return '';
    const t = [];
    if (x.r) t.push('<span class="tag na-rna' + (x.one ? '' : ' xl-near') + '" title="' + esc('RNA-crosslinked residue in living human cells (UV crosslinking mass spectrometry): ' + xlSets(x.r).join(', ') + (x.one ? '. At least one crosslink is localised to this residue alone.' : '. Only as one of several candidate residues of a crosslink.')) + '">RNA xl</span>');
    else if (x.n) t.push('<span class="tag na-rna xl-near" title="' + esc('Within 3 residues of an RNA-crosslinked residue: ' + xlSets(x.n).join(', ')) + '">near RNA xl</span>');
    if (x.f) t.push('<span class="tag na-rna xl-near" title="Inside an RNA-bound fragment (RBDmap, HeLa; Castello et al. 2016)">RBDmap</span>');
    if (x.d & 1) t.push('<span class="tag na-dna" title="DNA crosslink localised to this residue in living MCF7 cells (Trendel et al. 2025)">DNA xl</span>');
    else if (x.d & 8) t.push('<span class="tag na-dna xl-near" title="One of several candidate residues of a DNA crosslink in living MCF7 cells (Trendel et al. 2025)">DNA xl candidate</span>');
    else if (x.d & 2) t.push('<span class="tag na-dna xl-near" title="Inside a DNA-crosslinked peptide (Trendel et al. 2025)">DNA xl peptide</span>');
    if (x.d & 4) t.push('<span class="tag na-dna xl-near" title="Inside a UniProt DNA-binding region">DNA-binding domain</span>');
    return t.join(' ');
  }
  const xlText = (s) => (s.xl ? [s.xl.r ? 'RNA:' + xlSets(s.xl.r).join('/') : '', s.xl.n && !s.xl.r ? 'RNA<=3aa' : '', s.xl.f ? 'RBDmap' : '', s.xl.d & 1 ? 'DNA' : s.xl.d & 8 ? 'DNA candidate' : '', s.xl.d & 2 && !(s.xl.d & 9) ? 'DNA peptide' : '', s.xl.d & 4 ? 'UniProt DNA-binding' : ''].filter(Boolean).join(';') : '');
  const predText = (s) => (s.pred ? s.pred.hits.map(([k, sc]) => ixGene(k) + ':' + (sc / 100).toFixed(2)).join(';') : '');
  const modelledTag = (s) => (s.modelled ? '<span class="tag var" title="' + esc('Modelled as ' + s.modelled.comps.join('/') + ' in ' + s.modelled.pdbs.map((x) => x.toUpperCase()).join(', ')) + '">PDB ' + esc(s.modelled.comps[0]) + '</span>' : '');
  function contactHtml(s, max) {
    const gs = contactGroups(s);
    if (!gs.length) return s.geo ? '<div class="contacts">' + geoLine(s) + '</div>' : '';
    return '<div class="contacts">' + gs.slice(0, max || 4).map((g) => '<span title="' + esc(g.name + ': ' + g.labels.join(' ')) + '"><b class="' + g.cls + '">' + esc(g.name) + '</b> ' + esc(g.labels.slice(0, 6).join(' ')) + (g.labels.length > 6 ? ' +' + (g.labels.length - 6) : '') +
      ' <span class="muted">' + pdbButtons(s, g.pdbs) + '</span></span>').join('') + (gs.length > (max || 4) ? '<span class="muted" title="' + esc(gs.slice(max || 4).map((g) => g.name + ': ' + g.labels.join(' ')).join(' · ')) + '">+' + (gs.length - (max || 4)) + ' more partners</span>' : '') +
      geoLine(s) + '</div>';
  }
  function contactText(s) {
    return contactGroups(s).map((g) => g.name + ': ' + g.labels.join(' ') + ' (' + g.pdbs.join(',') + ')').join('; ');
  }

  // ----- live sources (cached per protein)
  const liveCache = new Map();
  function live(key, make) {
    if (!liveCache.has(key)) liveCache.set(key, make().catch(() => null));
    return liveCache.get(key);
  }
  const enspId = (k) => '9606.ENSP' + String(IX.prot[k].ensp).padStart(11, '0');

  // Edges among a set of proteins from the STRING API (physical subnetwork).
  function loadStringNet(ids, minScore) {
    const usable = ids.filter((k) => IX.prot[k] && IX.prot[k].ensp);
    if (usable.length < 2) return Promise.resolve([]);
    return live('string:' + usable.join(',') + ':' + minScore, () =>
      fetch('https://string-db.org/api/json/network', {
        method: 'POST',
        body: new URLSearchParams({ identifiers: usable.map(enspId).join('\r'), species: '9606', network_type: 'physical', required_score: String(minScore), caller_identity: STRING_CALLER }),
      }).then((r) => { if (!r.ok) throw new Error(r.status); return r.json(); })
        .then((rows) => {
          const out = new Map();
          for (const e of rows) {
            const a = IX.byEnsp.get(+String(e.stringId_A).replace(/^9606\.ENSP/, ''));
            const b = IX.byEnsp.get(+String(e.stringId_B).replace(/^9606\.ENSP/, ''));
            if (a == null || b == null || a === b) continue;
            const key = a < b ? a + ',' + b : b + ',' + a;
            out.set(key, Math.max(out.get(key) || 0, Math.round(e.score * 1000)));
          }
          return Array.from(out.entries()).map(([key, sc]) => key.split(',').map(Number).concat([sc]));
        }));
  }

  function loadQuickGO(acc, term) {
    return live('go:' + acc + ':' + term, () => fetchJson('https://www.ebi.ac.uk/QuickGO/services/annotation/search?geneProductId=' + encodeURIComponent(acc) +
      '&goId=' + term + '&goUsage=descendants&goUsageRelationships=is_a,part_of,occurs_in&limit=100&includeFields=goName')
      .then((d) => d.results || []));
  }

  function loadEncode(gene, assay) {
    return live('encode:' + gene + ':' + assay, () => fetchJson('https://www.encodeproject.org/search/?type=Experiment&status=released&format=json&limit=50' +
      '&assay_title=' + encodeURIComponent(assay) + '&target.label=' + encodeURIComponent(gene) +
      '&field=accession&field=biosample_summary&field=date_released&field=assay_title')
      .then((d) => d['@graph'] || []).catch((e) => (String(e.message) === '404' ? [] : Promise.reject(e))));
  }

  function loadChipAtlas(gene) {
    return live('chipatlas:' + gene, () => fetch('https://chip-atlas.dbcls.jp/data/hg38/target/' + encodeURIComponent(gene) + '.1.tsv')
      .then((r) => { if (r.status === 404) return null; if (!r.ok) throw new Error(r.status); return r.text(); })
      .then((t) => {
        if (!t) return { genes: [], experiments: 0 };
        const lines = t.trim().split('\n');
        const head = lines[0].split('\t');
        const genes = lines.slice(1).map((l) => { const c = l.split('\t'); return [c[0], +c[1]]; }).filter((g) => g[1] > 0);
        const cells = Array.from(new Set(head.slice(2).filter((h) => h.includes('|')).map((h) => h.split('|')[1])));
        return { genes, experiments: head.filter((h) => /^SRX|^ERX|^DRX/.test(h)).length, cells };
      }));
  }

  // PDBe-KB residues at interfaces with proteins, RNA and DNA, in UniProt numbering (all residues, not only sites).
  const IFDATA = new Map();
  function loadIface(acc) {
    const base = acc.split('-')[0];
    // PDBe can be slow; give up after 20 s and record the entry as unavailable rather than waiting forever.
    const ac = window.AbortController ? new AbortController() : null;
    if (ac) setTimeout(() => ac.abort(), 20000);
    return live('iface:' + base, () => fetch('https://www.ebi.ac.uk/pdbe/graph-api/uniprot/interface_residues/' + encodeURIComponent(base), ac ? { signal: ac.signal } : {})
      .then((r) => { if (r.status === 404) return {}; if (!r.ok) throw new Error(r.status); return r.json(); })
      .then((d) => {
        const rec = d[base];
        const out = { rna: [], dna: [], protein: [], self: [] };
        for (const part of (rec && rec.data) || []) {
          const key = part.accession === 'RNA' ? 'rna' : part.accession === 'DNA' ? 'dna' : part.accession === base ? 'self' : 'protein';
          for (const r of part.residues || []) out[key].push([r.startIndex, r.endIndex, part.name || part.accession]);
        }
        IFDATA.set(base, out);
        return out;
      })
      .catch(() => { const out = { rna: [], dna: [], protein: [], self: [], failed: true }; IFDATA.set(base, out); return out; }));
  }

  // ----- filters
  function ixProteinPass(p, f) {
    if (f.pevAll) { const m = p.pev || 0; if ((m & f.pevAll) !== f.pevAll || (f.pevExact && m !== f.pevAll)) return false; }
    return true;
  }

  function ixSitePass(s, f) {
    if (f.evAll) { const m = s.ev || 0; if ((m & f.evAll) !== f.evAll || (f.evExact && m !== f.evAll)) return false; }
    if (f.modelledPdb && !s.modelled) return false;
    if (f.pred && !(s.pred && s.pred.hits[0][1] >= f.pred)) return false;
    if (f.galOnly && !(s.gal || []).some((g) => galPass(g, UI.gal))) return false;
    if (f.xl.size) { const m = xlMask(s); let want = 0; f.xl.forEach((b) => { want |= b; }); if (!(m & want)) return false; }
    if (!f.iface.size) return true;
    let m = 0;
    f.iface.forEach((b) => { m |= b; });
    return !!(s.iface & m);
  }

  function ixConstraints(add) {
    if (F.ixScore !== 700 && !termValues('binds').length) add('ix', 'STRING confidence ≥ ' + (F.ixScore / 1000).toFixed(1), () => { F.ixScore = 700; });
    if (F.goExp && !termValues('go').length) add('ix', 'GO experimental evidence only', () => { F.goExp = false; });
    F.iface.forEach((v) => add('ix', 'Site at PDB interface with ' + IFACE_NAME[v], () => F.iface.delete(v)));
    if (F.modelledPdb) add('ix', 'Modelled as methyl-R in a PDB entry', () => { F.modelledPdb = false; });
    if (F.pred) add('ix', 'Predicted interface (pDockQ ≥ ' + (F.pred / 100).toFixed(2) + ')', () => { F.pred = 0; });
    if (F.galOnly) add('ix', 'Best-supported predicted contact', () => { F.galOnly = false; });
    F.xl.forEach((v) => add('ix', (XL_OPTS.find((o) => o[0] === v) || ['', 'Crosslink'])[1], () => F.xl.delete(v)));
    if (F.evAll) add('ix', 'Evidence: ' + evLabel(F.evAll) + (F.evExact ? ' (only these)' : ''), () => { F.evAll = 0; F.evExact = false; });
    if (F.pevAll) add('ix', 'Protein: ' + pevLabel(F.pevAll) + (F.pevExact ? ' (only these)' : ''), () => { F.pevAll = 0; F.pevExact = false; });
  }

  function ixRailBody(sel) {
    if (!IX) return '<p class="note">Loading interaction data…</p>';
    if (IX.failed) return '<p class="note">Interaction data (<code>data/methylarg_interactions.json</code>) is not available.</p>';
    // Less-used controls sit behind "+ More", opened whenever one of them is set.
    const more = (label, isSet, html) => '<details class="more-opts"' + (isSet ? ' open' : '') + '><summary>' + esc(label) + '</summary>' + html + '</details>';
    const pt = partnerToken(), ptText = pt ? pt.vals.map((v) => v.text).join(',') : '';
    return '<p class="note" style="margin-top:0">Site evidence and crosslink and interface filters keep individual sites; partner and Gene Ontology filters keep whole proteins. <button class="linkish" type="button" data-ixgo="">See the analyses</button></p>' +
      '<div class="fgroup"><span class="flabel" title="Experimental: RNA or DNA crosslink, RBDmap fragment, or a PDB interface at this residue. Predicted: at the interface of an AlphaFold-Multimer model (pDockQ ≥ 0.23) or a confident predicted contact.">Contact evidence at the site</span>' +
      termChips('interaction', [['experimental', 'Experimental', 'Crosslink, RBDmap or PDB interface at this residue'], ['predicted', 'Predicted', 'Predicted interface or confident predicted contact, with or without experimental evidence'],
        ['predicted-only', 'Predicted only', 'Predicted, with no experimental evidence yet: candidates to test'], ['none', 'None recorded', 'No contact evidence recorded; most proteins were never crystallised or crosslinked']]) + '</div>' +
      '<div class="fgroup"><span class="flabel" title="STRING v' + esc(IX.meta.string_version) + ' physical subnetwork">Has a physical partner that is (match any)</span>' + termChips('binds', BIND_TERMS.map(([t, v]) => { const o = BIND_OPTS.find((x) => x[0] === v); return [t, o[1], o[2]]; })) +
      more('A specific partner, confidence', !!pt || F.ixScore !== 700,
        '<label class="flabel" for="f-partner">Binds a specific protein (any gene)</label>' +
        '<div class="motif-search"><input id="f-partner" type="text" list="f-partner-list" autocomplete="off" spellcheck="false" placeholder="e.g. PRMT5, SMN1, WDR5" value="' + esc(ptText) + '"' + (pt ? ' class="set"' : '') + ' aria-describedby="f-partnermsg">' +
        '<button class="btn small" type="button" id="f-partnergo">Apply</button></div><datalist id="f-partner-list"></datalist>' +
        '<p class="note" id="f-partnermsg">' + esc(partnerHint()) + '</p>' +
        '<label class="flabel" for="f-ixscore">STRING physical confidence</label>' + sel('f-ixscore', options(RAIL_SCORES, F.ixScore), false)) + '</div>' +
      '<div class="fgroup"><span class="flabel" title="UV crosslinking mass spectrometry in living human cells: pepR-MS compendium (Kramer 2014, XRNAX, RBS-ID, pRBS-ID, pepR-MS), RBDmap, Trendel 2025 DNA crosslinks; and UniProt DNA-binding regions">Site contacts RNA/DNA in cells (match any)</span>' + chipRow('xl', XL_OPTS, F.xl) + '</div>' +
      '<div class="fgroup"><span class="flabel" title="PDBe-KB: in at least one PDB structure the residue contacts this kind of partner">Site at a PDB interface with (match any)</span>' + chipRow('iface', IFACE_OPTS, F.iface) +
      more('Predicted complexes, modelled methyl-R', !!(F.pred || F.galOnly || F.modelledPdb),
        '<label class="flabel" for="f-pred" title="ProtVar: residues within 5 Å of the partner in AlphaFold-Multimer models of human protein pairs. pDockQ estimates model quality; 0.23 is the usual floor for an acceptable interface, 0.5 for a confident one.">Site at a predicted interface (ProtVar)</label>' +
        sel('f-pred', options(PRED_OPTS, F.pred), !!F.pred) +
        '<label class="check"><input type="checkbox" id="f-gal"' + (F.galOnly ? ' checked' : '') + '> Best-supported predicted contact (salt bridge or H-bond, confident model, pLDDT ≥ ' + UI.gal.pl + ')</label>' +
        '<label class="check"><input type="checkbox" id="f-modelled"' + (F.modelledPdb ? ' checked' : '') + '> Modelled as a methylated arginine (2MR, DA2, NMM, AGM) in a PDB entry</label>') + '</div>' +
      '<div class="fgroup"><span class="flabel" title="Gene Ontology GO:0003723 (RNA binding) or GO:0003677 (DNA binding), including child terms">Protein binds (Gene Ontology, match any)</span>' + termChips('go', [['rna', 'RNA', 'GO:0003723 or a child term'], ['dna', 'DNA', 'GO:0003677 or a child term']]) +
      more('Evidence', F.goExp, '<label class="check"><input type="checkbox" id="f-goexp"' + (F.goExp ? ' checked' : '') + '> Experimental evidence only (not electronic or inferred)</label>') + '</div>';
  }

  // Open a PDB entry from a contact list in the 3D panel, focused on the site.
  function openContactStructure(p, pdb, pos) {
    const card = $('#card-af');
    if (card) card.scrollIntoView({ block: 'start', behavior: 'smooth' });
    if (UI.startStructure) UI.startStructure();
    loadPdb(p.acc).then((res) => {
      const e = res.entries.find((x) => x.id === pdb);
      if (!e) { window.open('https://www.rcsb.org/structure/' + pdb.toUpperCase(), '_blank', 'noopener'); return null; }
      return showPdb(p, e, null);
    }).then(() => { if (UI.detail === p.i && mol && mol.mode === 'pdb') { UI.highlight = pos; focusStructure(pos); } });
  }
  function bindContactLinks(box) {
    $$('[data-cpdb]', box).forEach((b) => b.addEventListener('click', (e) => {
      e.stopPropagation();
      const pi = +b.dataset.cp, pos = +b.dataset.cpos, pdb = b.dataset.cpdb;
      if (UI.view === 'protein' && UI.detail === pi) openContactStructure(P[pi], pdb, pos);
      else { UI.pendingContact = { pdb, pos }; openProtein(pi, pos); }
    }));
  }

  // Methylated proteins each protein binds at a STRING score (for the partner picker and the overview chart).
  const boundByCache = new Map();
  function boundBy(min) {
    if (!boundByCache.has(min)) {
      const m = new Map();
      for (const p of P) for (const [k, sc] of p.partners || []) if (sc >= min) m.set(k, (m.get(k) || 0) + 1);
      boundByCache.set(min, m);
    }
    return boundByCache.get(min);
  }
  function partnerClassText(x) {
    const out = [];
    if (x.flags & IXF.prmt) out.push('PRMT');
    if (x.flags & IXF.tudor) out.push('Tudor');
    if (x.flags & IXF.wd40) out.push('WD40');
    if (x.flags & IXF.phd) out.push('PHD');
    return out.join(', ');
  }
  function partnerHint() {
    if (!ixReady()) return '';
    const tk = partnerToken(), v = tk && tk.vals.find((x) => x.partner != null);
    if (!v) return 'Type a gene symbol; adds binds:GENE to the search. Suggestions show how many methylated proteins it binds.';
    const n = boundBy(F.ixScore).get(v.partner) || 0;
    return ixGene(v.partner) + ' binds ' + fmt(n) + ' methylated protein' + (n === 1 ? '' : 's') + ' at STRING ≥ ' + (F.ixScore / 1000).toFixed(1) + '.';
  }
  function partnerByGene(q) {
    q = q.trim().toUpperCase();
    if (!q) return -1;
    const bound = boundBy(F.ixScore);
    let best = -1;
    for (const x of IX.prot) {
      if (x.gene.toUpperCase() !== q) continue;
      if (best < 0 || (bound.get(x.k) || 0) > (bound.get(best) || 0)) best = x.k;
    }
    return best;
  }
  function bindIxRail(rail, refresh) {
    bindIxGo(rail);
    const part = $('#f-partner', rail);
    if (part) {
      const list = $('#f-partner-list', rail), msg = $('#f-partnermsg', rail);
      const suggest = () => {
        const q = part.value.trim().toUpperCase();
        if (!q) { list.innerHTML = ''; return; }
        const bound = boundBy(F.ixScore);
        const hits = IX.prot.filter((x) => bound.get(x.k) && x.gene.toUpperCase().startsWith(q))
          .sort((a, b) => (a.gene.length - b.gene.length) || (bound.get(b.k) - bound.get(a.k)) || a.gene.localeCompare(b.gene)).slice(0, 25);
        list.innerHTML = hits.map((x) => { const c = partnerClassText(x); return '<option value="' + esc(x.gene) + '" label="' + esc((c ? c + ' · ' : '') + fmt(bound.get(x.k)) + ' methylated partners') + '"></option>'; }).join('');
      };
      const apply = () => {
        const q = part.value.trim().replace(/\s+/g, '');
        const bad = q.split(',').filter(Boolean).map((g) => resolveValue('binds', g)).find((r) => r.error);
        if (bad) { part.setAttribute('aria-invalid', 'true'); msg.textContent = bad.error; return; }
        part.removeAttribute('aria-invalid');
        replaceToken(partnerToken(), 'binds', q.toUpperCase());
        msg.textContent = partnerHint();
      };
      part.addEventListener('input', (e) => {
        suggest();
        // Picking a suggestion (no typed input type) applies it at once.
        if (!e.inputType || e.inputType === 'insertReplacementText') apply();
      });
      part.addEventListener('keydown', (e) => { if (e.key === 'Enter') { e.preventDefault(); apply(); } });
      $('#f-partnergo', rail).addEventListener('click', apply);
    }
    const sc = $('#f-ixscore', rail);
    if (sc) sc.addEventListener('change', () => { F.ixScore = +sc.value || 700; refresh(); });
    const ge = $('#f-goexp', rail);
    if (ge) ge.addEventListener('change', () => { F.goExp = ge.checked; refresh(); });
    const mo = $('#f-modelled', rail);
    if (mo) mo.addEventListener('change', () => { F.modelledPdb = mo.checked; refresh(); });
    const pr = $('#f-pred', rail);
    if (pr) pr.addEventListener('change', () => { F.pred = +pr.value || 0; refresh(); });
    const gl = $('#f-gal', rail);
    if (gl) gl.addEventListener('change', () => { F.galOnly = gl.checked; refresh(); });
  }

  // ----- protein network
  const NODE_MODES = {
    meth: { name: 'Methylated', cats: ['In this methylarginine set'], cat: () => 0, colors: (pal) => [pal.site] },
    sites: {
      name: 'Methyl-R sites', cats: ['1 site', '2–4 sites', '5–9 sites', '10 or more'],
      cat: (p) => { const n = p.sites.length; return n >= 10 ? 3 : n >= 5 ? 2 : n >= 2 ? 1 : 0; }, colors: (pal) => pal.seq,
    },
    state: { name: 'Methylation state', cats: ['Rme1', 'Rme2', 'Rme1 + Rme2', 'No state call'], cat: (p) => p.stCat, colors: (pal) => [pal.cat[0], pal.cat[1], pal.cat[2], pal.none] },
    prmti: { name: 'PRMT inhibitor response', cats: ['↓ GSK591', '↓ MS023', '↓ both', 'Measured, n.s.', 'Not measured'], cat: (p) => p.trCat, colors: (pal) => [pal.cat[0], pal.cat[1], pal.cat[2], pal.dark, pal.none] },
  };

  // Radial force layout: the protein is pinned at the centre, partners sit at a distance that shrinks with
  // their STRING score, links among partners pull gently, nodes and labels repel, then the layout is
  // stretched to fill the panel and residual overlaps are pushed apart. Deterministic for a given input.
  function layoutNetwork(nodes, edges, W, H) {
    const n = nodes.length;
    const cx = W / 2, cy = H / 2;
    const R = Math.min(W, H) * 0.42;
    const pos = nodes.map((nd, k) => {
      if (!k) return { x: cx, y: cy };
      const a = ((k - 1) / Math.max(1, n - 1)) * 2 * Math.PI - Math.PI / 2;
      return { x: cx + R * Math.cos(a), y: cy + R * Math.sin(a) };
    });
    const star = new Map(edges.filter((e) => e[0] === 0).map((e) => [e[1], e[2]]));
    const other = edges.filter((e) => e[0] !== 0);
    let temp = R / 3;
    for (let it = 0; it < 360; it++) {
      const fx = new Float64Array(n), fy = new Float64Array(n);
      for (let i = 1; i < n; i++) {
        for (let j = i + 1; j < n; j++) {
          const dx = pos[i].x - pos[j].x, dy = pos[i].y - pos[j].y;
          const d2 = Math.max(25, dx * dx + dy * dy);
          const f = (R * R * 0.09) / d2;
          fx[i] += dx * f; fy[i] += dy * f; fx[j] -= dx * f; fy[j] -= dy * f;
        }
        // Spring to the centre with a rest length set by the score (stronger link = closer).
        const w = star.has(i) ? star.get(i) : 0.4;
        const dx = pos[i].x - cx, dy = pos[i].y - cy;
        const d = Math.max(1, Math.hypot(dx, dy));
        const L = R * (1.15 - 0.55 * w);
        fx[i] -= (dx / d) * (d - L) * 0.12; fy[i] -= (dy / d) * (d - L) * 0.12;
      }
      for (const [a, b, w] of other) {
        const dx = pos[a].x - pos[b].x, dy = pos[a].y - pos[b].y;
        const d = Math.max(1, Math.hypot(dx, dy));
        const f = (d - R * 0.5) * 0.012 * w;
        fx[a] -= (dx / d) * f; fy[a] -= (dy / d) * f; fx[b] += (dx / d) * f; fy[b] += (dy / d) * f;
      }
      for (let i = 1; i < n; i++) {
        const len = Math.hypot(fx[i], fy[i]);
        if (len > 1e-9) { const st = Math.min(len, temp); pos[i].x += (fx[i] / len) * st; pos[i].y += (fy[i] / len) * st; }
      }
      temp = Math.max(0.5, temp * 0.986);
    }
    // Stretch to the panel (the centre stays put), keeping room for labels.
    const mx = 46, myTop = 30, myBot = 34;
    let ex = 1, ey = 1;
    for (let i = 1; i < n; i++) { ex = Math.max(ex, Math.abs(pos[i].x - cx)); ey = Math.max(ey, Math.abs(pos[i].y - cy)); }
    const sx = (W / 2 - mx) / ex, sy = (H / 2 - Math.max(myTop, myBot)) / ey;
    for (let i = 1; i < n; i++) { pos[i].x = cx + (pos[i].x - cx) * sx; pos[i].y = cy + (pos[i].y - cy) * sy; }
    // Push apart anything still overlapping (node discs plus their labels).
    const box = (k) => ({ hw: Math.max(nodes[k].r, nodes[k].label.length * 3.4) + 4, hh: nodes[k].r + 9 });
    for (let it = 0; it < 120; it++) {
      let moved = false;
      for (let i = 0; i < n; i++) {
        for (let j = i + 1; j < n; j++) {
          const bi = box(i), bj = box(j);
          const dx = pos[j].x - pos[i].x, dy = pos[j].y - pos[i].y;
          const ox = bi.hw + bj.hw - Math.abs(dx), oy = bi.hh + bj.hh - Math.abs(dy);
          if (ox <= 0 || oy <= 0) continue;
          moved = true;
          const alongX = ox < oy;
          const push = (alongX ? ox : oy) / 2 + 0.5;
          const sgn = alongX ? (dx >= 0 ? 1 : -1) : (dy >= 0 ? 1 : -1);
          if (i === 0) { if (alongX) pos[j].x += sgn * push * 2; else pos[j].y += sgn * push * 2; }
          else if (alongX) { pos[i].x -= sgn * push; pos[j].x += sgn * push; }
          else { pos[i].y -= sgn * push; pos[j].y += sgn * push; }
        }
      }
      for (let i = 1; i < n; i++) {
        pos[i].x = Math.min(W - mx, Math.max(mx, pos[i].x));
        pos[i].y = Math.min(H - myBot, Math.max(myTop, pos[i].y));
      }
      if (!moved) break;
    }
    return pos;
  }

  function nodeRadius(k, focal) {
    if (focal) return 17;
    const x = IX.prot[k];
    return x.p >= 0 ? Math.min(15, 6 + 2.2 * Math.sqrt(P[x.p].sites.length)) : 6;
  }

  function ixNetModel(p, all, min, extra) {
    const list = all.slice(0, UI.ixTop);
    const ids = [p.i].concat(list.map(([k]) => k));
    const idx = new Map(ids.map((k, n) => [k, n]));
    const edges = list.map(([k, sc]) => [0, idx.get(k), sc]);
    if (extra) for (const [a, b, sc] of extra) if (idx.has(a) && idx.has(b) && a !== p.i && b !== p.i) edges.push([idx.get(a), idx.get(b), sc]);
    // Second ring: PRMTs and Tudor-domain proteins that physically bind one or more of the partners
    // (writer-substrate-reader paths through this protein's neighbourhood).
    const ring = [];
    if (UI.ixRing && IX.wr) {
      const inNet = new Set(ids);
      IX.wr.forEach((set, w) => {
        if (inNet.has(w)) return;
        const hits = ids.slice(1).map((k, n) => (set.has(k) ? n + 1 : -1)).filter((n) => n > 0);
        if (hits.length) ring.push({ k: w, hits });
      });
      ring.sort((a, b) => b.hits.length - a.hits.length || ((IX.prot[b.k].flags & IXF.prmt) - (IX.prot[a.k].flags & IXF.prmt)) || IX.prot[a.k].gene.localeCompare(IX.prot[b.k].gene));
      ring.splice(14);
    }
    const W = 900, H = ring.length ? 660 : 540;
    const nodes = ids.map((k, n) => ({ r: nodeRadius(k, n === 0), label: IX.prot[k].gene }));
    const pos = layoutNetwork(nodes, edges.map(([a, b, sc]) => [a, b, sc / 1000]), W, H);
    const ringEdges = [];
    if (ring.length) {
      // Pull the core in, then seat each ring protein on an outer ellipse near the partners it binds.
      const cx = W / 2, cy = H / 2;
      pos.forEach((q) => { q.x = cx + (q.x - cx) * 0.66; q.y = cy + (q.y - cy) * 0.64; });
      const rx = W / 2 - 58, ry = H / 2 - 34;
      const items = ring.map((r) => {
        let sx = 0, sy = 0;
        r.hits.forEach((n) => { const a = Math.atan2(pos[n].y - cy, pos[n].x - cx); sx += Math.cos(a); sy += Math.sin(a); });
        return { ...r, a: Math.atan2(sy, sx), need: (IX.prot[r.k].gene.length * 6.6 + 26) / ((rx + ry) / 2) };
      }).sort((x, y) => x.a - y.a);
      for (let it = 0; it < 80; it++) {
        for (let i = 0; i < items.length && items.length > 1; i++) {
          const j = (i + 1) % items.length;
          let gap = items[j].a - items[i].a;
          if (j === 0) gap += 2 * Math.PI;
          const req = (items[i].need + items[j].need) / 2;
          if (gap < req) { const d = (req - gap) / 2; items[i].a -= d; items[j].a += d; }
        }
      }
      items.forEach((it) => {
        const n = ids.length;
        ids.push(it.k);
        pos.push({ x: cx + rx * Math.cos(it.a), y: cy + ry * Math.sin(it.a) });
        it.hits.forEach((h) => ringEdges.push([n, h]));
      });
    }
    return { p, ids, edges, ringEdges, nCore: ids.length - ring.length, ring, pos, W, H, live: !!extra, total: all.length, min };
  }

  function ixNetSvg(net, pal, forExport) {
    const { p, ids, edges, pos } = net;
    const mode = NODE_MODES[UI.ixColor] || NODE_MODES.meth;
    const colors = mode.colors(pal);
    const W = net.W, top = forExport ? 54 : 0, H = net.H + top + (forExport ? 70 : 0);
    let g = '';
    if (forExport) {
      g += '<text x="16" y="24" font-size="17" font-weight="700" fill="' + pal.ink + '">' + esc(p.label) + ' physical interaction partners</text>' +
        '<text x="16" y="42" font-size="11" fill="' + pal.muted + '">STRING v' + esc(IX.meta.string_version) + ' physical subnetwork, score ≥ ' + (net.min / 1000).toFixed(1) + ' · top ' + ((net.nCore || ids.length) - 1) + ' of ' + net.total + ' partners · nodes coloured by ' + esc(mode.name.toLowerCase()) + ' in this methylarginine set</text>';
    }
    g += '<g transform="translate(0,' + top + ')">';
    for (const [a, b] of net.ringEdges || []) {
      g += '<line x1="' + r1(pos[a].x) + '" y1="' + r1(pos[a].y) + '" x2="' + r1(pos[b].x) + '" y2="' + r1(pos[b].y) + '" stroke="' + pal.accent + '" stroke-width="1" stroke-dasharray="4 3" stroke-opacity="0.55"/>';
    }
    for (const [a, b, sc] of edges) {
      const w = 0.6 + 2.6 * Math.max(0, (sc - 400) / 600);
      g += '<line x1="' + r1(pos[a].x) + '" y1="' + r1(pos[a].y) + '" x2="' + r1(pos[b].x) + '" y2="' + r1(pos[b].y) + '" stroke="' + (a === 0 ? pal.backbone : pal.grid) + '" stroke-width="' + r1(w) + '"' + (a === 0 ? '' : ' stroke-opacity="0.9"') + '/>';
    }
    ids.forEach((k, n) => {
      const x = IX.prot[k];
      const q = x.p >= 0 ? P[x.p] : null;
      const outer = n >= (net.nCore || ids.length);
      const r = outer ? Math.min(9, nodeRadius(k, false)) : nodeRadius(k, n === 0);
      const fill = q ? colors[mode.cat(q)] || pal.none : pal.bg;
      const stroke = q ? pal.ring : pal.backbone;
      const cxp = r1(pos[n].x), cyp = r1(pos[n].y);
      const binds = outer ? (net.ringEdges || []).filter((e) => e[0] === n).map((e) => IX.prot[ids[e[1]]].gene) : [];
      const tipHtml = '<b>' + esc(x.gene) + '</b> ' + esc(x.acc) + (outer ? '<br>Second ring: binds ' + esc(binds.join(', ')) + ' (STRING ≥ 0.7)' : n ? '<br>STRING physical score ' + (net.edges.find((e) => e[0] === 0 && e[1] === n) || [0, 0, 0])[2] / 1000 : '') +
        '<br>' + (q ? q.sites.length + ' methylarginine site' + (q.sites.length === 1 ? '' : 's') + ' in this set' : 'Not in this methylarginine set') +
        (x.flags & IXF.prmt ? '<br>Protein arginine methyltransferase' : '') + (x.flags & IXF.tudor ? '<br>Tudor domain' : '') +
        (x.flags & IXF.wd40 ? '<br>WD40 repeat' : '') + (x.flags & IXF.phd ? '<br>PHD finger' : '') +
        (x.flags & IXF.rna ? '<br>RNA binding (GO)' : '') + (x.flags & IXF.dna ? '<br>DNA binding (GO)' : '');
      g += '<g class="nd" data-k="' + k + '"' + (forExport ? '' : ' tabindex="0" role="button" data-tip="' + esc(tipHtml) + '"') + ' style="cursor:pointer">';
      let ringR = r + 4;
      if (x.flags & IXF.tudor) { g += '<circle cx="' + cxp + '" cy="' + cyp + '" r="' + ringR + '" fill="none" stroke="' + pal.accent + '" stroke-width="1.4" stroke-dasharray="3 2"/>'; ringR += 3; }
      if (x.flags & IXF.wd40) { g += '<circle cx="' + cxp + '" cy="' + cyp + '" r="' + ringR + '" fill="none" stroke="' + IX_WD_COLOR + '" stroke-width="1.4"/>'; ringR += 3; }
      if (x.flags & IXF.phd) g += '<circle cx="' + cxp + '" cy="' + cyp + '" r="' + ringR + '" fill="none" stroke="' + IX_PHD_COLOR + '" stroke-width="1.6" stroke-dasharray="1 2" stroke-linecap="round"/>';
      if (x.flags & IXF.prmt) g += '<rect x="' + r1(pos[n].x - r) + '" y="' + r1(pos[n].y - r) + '" width="' + r1(2 * r) + '" height="' + r1(2 * r) + '" rx="3" transform="rotate(45 ' + cxp + ' ' + cyp + ')" fill="' + fill + '" stroke="' + pal.ink + '" stroke-width="1.6"/>';
      else g += '<circle cx="' + cxp + '" cy="' + cyp + '" r="' + r1(r) + '" fill="' + fill + '" stroke="' + (n === 0 ? pal.ink : stroke) + '" stroke-width="' + (n === 0 ? 2 : q ? 1.2 : 1.4) + '"/>';
      g += '<text x="' + cxp + '" y="' + r1(pos[n].y + r + 12) + '" font-size="' + (n === 0 ? 12.5 : outer ? 9.5 : 10.5) + '" font-weight="' + (n === 0 || (q && !outer) ? 700 : 400) + '" text-anchor="middle" fill="' + (q && !outer ? pal.ink : pal.muted) + '"' + (outer ? ' font-style="italic"' : '') + '>' + esc(x.gene) + '</text></g>';
    });
    g += '</g>';
    if (forExport) {
      let lx = 16;
      const ly = top + net.H + 22;
      g += '<text x="16" y="' + ly + '" font-size="10.5" font-weight="700" fill="' + pal.muted + '">' + esc(mode.name) + '</text>';
      lx = 16 + mode.name.length * 6.2 + 12;
      mode.cats.forEach((c, k) => { g += '<circle cx="' + (lx + 5) + '" cy="' + (ly - 4) + '" r="5" fill="' + colors[k] + '"/><text x="' + (lx + 14) + '" y="' + ly + '" font-size="10.5" fill="' + pal.ink + '">' + esc(c) + '</text>'; lx += 26 + c.length * 5.9; });
      g += '<circle cx="' + (lx + 5) + '" cy="' + (ly - 4) + '" r="5" fill="' + pal.bg + '" stroke="' + pal.backbone + '"/><text x="' + (lx + 14) + '" y="' + ly + '" font-size="10.5" fill="' + pal.ink + '">Not in this set</text>';
      const ly2 = ly + 20;
      g += '<rect x="18" y="' + (ly2 - 10) + '" width="10" height="10" transform="rotate(45 23 ' + (ly2 - 5) + ')" fill="none" stroke="' + pal.ink + '"/><text x="34" y="' + ly2 + '" font-size="10.5" fill="' + pal.ink + '">PRMT</text>' +
        '<circle cx="90" cy="' + (ly2 - 4) + '" r="7" fill="none" stroke="' + pal.accent + '" stroke-dasharray="3 2"/><text x="102" y="' + ly2 + '" font-size="10.5" fill="' + pal.ink + '">Tudor domain</text>' +
        '<circle cx="190" cy="' + (ly2 - 4) + '" r="7" fill="none" stroke="' + IX_WD_COLOR + '"/><text x="202" y="' + ly2 + '" font-size="10.5" fill="' + pal.ink + '">WD40 repeat</text>' +
        '<circle cx="290" cy="' + (ly2 - 4) + '" r="7" fill="none" stroke="' + IX_PHD_COLOR + '" stroke-width="1.6" stroke-dasharray="1 2" stroke-linecap="round"/><text x="302" y="' + ly2 + '" font-size="10.5" fill="' + pal.ink + '">PHD finger</text>' +
        (net.ring && net.ring.length ? '<line x1="380" y1="' + (ly2 - 4) + '" x2="404" y2="' + (ly2 - 4) + '" stroke="' + pal.accent + '" stroke-dasharray="4 3"/><text x="410" y="' + ly2 + '" font-size="10.5" font-style="italic" fill="' + pal.ink + '">Second ring: PRMTs and Tudor proteins that bind the partners</text>' : '') +
        '<text x="' + (W - 16) + '" y="' + (ly2 + 18) + '" font-size="9" text-anchor="end" fill="' + pal.muted + '">Edges: STRING physical links (CC BY 4.0) · methylarginine set: PTM_Pipeline</text>';
    }
    return '<svg xmlns="http://www.w3.org/2000/svg" class="ix-net" width="' + W + '" height="' + H + '" viewBox="0 0 ' + W + ' ' + H + '" font-family=\'' + (pal.font || EXPORT_PAL.font).replace(/'/g, '"') + '\'' +
      ' role="img" aria-label="Physical interaction partners of ' + esc(p.label) + '">' + (forExport ? '<rect width="' + W + '" height="' + H + '" fill="' + pal.bg + '"/>' : '') + g + '</svg>';
  }

  // ----- protein page card
  function renderIxCard(p) {
    const box = $('#card-ix');
    if (!box) return;
    const head = '<div class="card-head"><div><h3>Interactions</h3><p class="sub">Physical partners (STRING), protein complexes (Complex Portal), RNA and DNA binding (Gene Ontology, ENCODE, ChIP-Atlas) and contacts in PDB structures (PDBe-KB). Partners are coloured by their own arginine methylation in this dataset.</p></div><div class="exp" id="ix-exp"></div></div>';
    if (!IX) { box.innerHTML = head + '<p class="note">Loading interaction data…</p>'; return; }
    if (IX.failed) { box.innerHTML = head + '<p class="note">Interaction data is not available.</p>'; return; }
    const tabs = [['proteins', 'Protein partners', partnersAt(p, 700).length], ['complexes', 'Complexes', p.cx.length], ['rna', 'RNA', null], ['dna', 'DNA', null]];
    const ev = EV_SETS.map(([b, l]) => [b, l, p.sites.filter((j) => (S[j].ev || 0) & b).length]).filter((x) => x[2]);
    const strip = ev.length ? '<p class="ev-strip"><span class="muted small">Evidence at its ' + fmt(p.sites.length) + ' methylarginines:</span> ' +
      ev.map(([b, l, n]) => '<button class="tag ev-tag" type="button" data-evp="' + b + '" title="List the ' + n + ' sites of ' + esc(p.label) + ' with: ' + esc(l) + '">' + esc(l) + ' <b>' + n + '</b></button>').join(' ') + '</p>' : '';
    box.innerHTML = head + strip + '<div class="seg ix-tabs" role="tablist">' + tabs.map(([k, l, n]) => '<button type="button" role="tab" data-ixtab="' + k + '" aria-pressed="' + (UI.ixTab === k) + '">' + esc(l) + (n != null ? ' <span class="n">' + fmt(n) + '</span>' : '') + '</button>').join('') + '</div><div id="ix-body"></div>';
    $$('[data-ixtab]', box).forEach((b) => b.addEventListener('click', () => { UI.ixTab = b.dataset.ixtab; renderIxCard(p); }));
    $$('[data-evp]', box).forEach((b) => b.addEventListener('click', () => { setPlain(p.acc); setFilters((f) => { f.evAll = +b.dataset.evp; f.evExact = false; }, 'sites'); window.scrollTo({ top: 0 }); }));
    const body = $('#ix-body', box);
    if (UI.ixTab === 'complexes') renderIxComplexes(p, body);
    else if (UI.ixTab === 'rna' || UI.ixTab === 'dna') renderIxNucleic(p, body, UI.ixTab);
    else renderIxProteins(p, body);
  }

  function partnerTags(x) {
    return (x.flags & IXF.prmt ? '<span class="tag ix-prmt">PRMT</span> ' : '') + (x.flags & IXF.tudor ? '<span class="tag ix-tudor">Tudor</span> ' : '') +
      (x.flags & IXF.wd40 ? '<span class="tag ix-dom" title="WD40 repeat (InterPro IPR001680)">WD40</span> ' : '') + (x.flags & IXF.phd ? '<span class="tag ix-dom" title="PHD finger (InterPro IPR001965/IPR019787)">PHD</span> ' : '') +
      (x.flags & IXF.rna ? '<span class="tag na-rna" title="GO RNA binding' + (x.flags & IXF.rnaExp ? ', experimental evidence' : '') + '">RNA' + (x.flags & IXF.rnaExp ? '' : '*') + '</span> ' : '') +
      (x.flags & IXF.dna ? '<span class="tag na-dna" title="GO DNA binding' + (x.flags & IXF.dnaExp ? ', experimental evidence' : '') + '">DNA' + (x.flags & IXF.dnaExp ? '' : '*') + '</span>' : '');
  }

  // Partners at a STRING score: from the static export (>= its floor) or live from the STRING API below it.
  function partnersFor(p, min) {
    const floor = IX.meta.string_min_score;
    const fallback = () => ({ list: partnersAt(p, floor), live: false, failed: true });
    if (min >= floor) return Promise.resolve({ list: partnersAt(p, min), live: false });
    const me = IX.prot[p.i];
    if (!me || !me.ensp) return Promise.resolve(fallback());
    return live('partners:' + p.i + ':' + min, () => fetch('https://string-db.org/api/json/interaction_partners', {
      method: 'POST',
      body: new URLSearchParams({ identifiers: enspId(p.i), species: '9606', network_type: 'physical', required_score: String(min), limit: '2000', caller_identity: STRING_CALLER }),
    }).then((r) => { if (!r.ok) throw new Error(r.status); return r.json(); })
      .then((rows) => {
        const seen = new Map();
        for (const e of rows) {
          const k = IX.byEnsp.get(+String(e.stringId_B).replace(/^9606\.ENSP/, ''));
          if (k == null || k === p.i) continue;
          seen.set(k, Math.max(seen.get(k) || 0, Math.round(e.score * 1000)));
        }
        return Array.from(seen.entries()).sort((x, y) => y[1] - x[1]);
      })).then((list) => (list ? { list, live: true } : fallback()));
  }

  function renderIxProteins(p, body) {
    const mode = NODE_MODES[UI.ixColor] || NODE_MODES.meth;
    const pal = screenPalette();
    const colors = mode.colors(pal);
    const exp = $('#ix-exp');
    if (exp) exp.innerHTML = '<button type="button" id="ix-svg">SVG</button><button type="button" id="ix-png">PNG</button><button type="button" id="ix-tsv">TSV</button>';
    body.innerHTML = '<div class="fig-controls">' +
      '<label class="fig-ctl">Show top <select id="ix-top">' + [10, 20, 30, 50].map((n) => '<option value="' + n + '"' + (UI.ixTop === n ? ' selected' : '') + '>' + n + '</option>').join('') + '</select></label>' +
      '<label class="fig-ctl">Confidence <select id="ix-min">' + IX_SCORES.map(([v, l]) => '<option value="' + v + '"' + (UI.ixMin === v ? ' selected' : '') + '>' + esc(l) + '</option>').join('') + '</select></label>' +
      '<label class="fig-ctl">Colour partners by <select id="ix-color">' + Object.keys(NODE_MODES).map((k) => '<option value="' + k + '"' + (UI.ixColor === k ? ' selected' : '') + '>' + esc(NODE_MODES[k].name) + '</option>').join('') + '</select></label>' +
      '<label class="check fig-ctl" title="Add the PRMTs and Tudor-domain proteins that bind these partners (STRING ≥ 0.7), to trace writer–substrate–reader paths"><input type="checkbox" id="ix-ring"' + (UI.ixRing ? ' checked' : '') + '> Second ring: PRMTs and Tudor proteins</label></div>' +
      '<div class="legend"><span class="legend-title">' + esc(mode.name) + '</span>' + mode.cats.map((c, k) => '<span><i style="background:' + colors[k] + '"></i>' + esc(c) + '</span>').join('') +
      '<span><i class="hollow"></i>Not in this set</span><span><i class="diamond"></i>PRMT</span><span><i class="ring"></i>Tudor domain</span><span><i class="ring wd"></i>WD40 repeat</span><span><i class="ring phd"></i>PHD finger</span></div>' +
      '<div id="ix-partners"><p class="note">Loading STRING partners…</p></div>';
    $('#ix-top', body).addEventListener('change', (e) => { UI.ixTop = +e.target.value; renderIxProteins(p, body); });
    $('#ix-min', body).addEventListener('change', (e) => { UI.ixMin = +e.target.value; renderIxProteins(p, body); });
    $('#ix-color', body).addEventListener('change', (e) => { UI.ixColor = e.target.value; renderIxProteins(p, body); });
    $('#ix-ring', body).addEventListener('change', (e) => { UI.ixRing = e.target.checked; renderIxProteins(p, body); });
    UI.ixNet = null;
    const minAsked = UI.ixMin;
    partnersFor(p, minAsked).then((res) => {
      const box = $('#ix-partners', body);
      if (!box || UI.detail !== p.i || UI.ixTab !== 'proteins' || UI.ixMin !== minAsked) return;
      const all = res.list;
      const min = res.failed ? IX.meta.string_min_score : minAsked;
      const nMeth = all.filter(([k]) => k < IX.nP).length;
      const predMap = new Map(predPartners(p).map((r) => [r.k, r]));
      box.innerHTML = (res.failed && minAsked < IX.meta.string_min_score ? '<p class="note">STRING did not answer; showing the high-confidence (≥ 0.7) links from the export instead.</p>' : '') +
        (all.length ? '<div class="ix-net-box" id="ix-net"></div><p class="note" id="ix-net-note"></p><div id="ix-paths"></div>' : '<p class="note">No STRING physical partners at this confidence.</p>') +
        '<p class="sub">' + fmt(all.length) + ' physical partners at score ≥ ' + (min / 1000).toFixed(1) + (res.live ? ' (live from STRING)' : '') + '; <b>' + fmt(nMeth) + '</b> (' + pctText(nMeth, all.length || 1) + '%) are themselves in this methylarginine set, against ' +
        pctText(IX.meta.methylated_in_proteome, IX.meta.proteome_size) + '% of the reviewed proteome. ' +
        (nMeth ? '<button class="linkish" type="button" id="ix-list">List the methylated partners' + (min < 700 ? ' (≥ 0.7)' : '') + '</button>' : '') + '</p>' +
        (all.length ? '<div class="tbl-wrap ix-tbl"><table class="tbl compact"><thead><tr><th>Partner</th><th class="r">STRING</th><th>Methyl-R sites</th><th>Annotations</th><th class="hide-sm" title="ProtVar AlphaFold-Multimer model of this pair places a methylarginine of ' + esc(p.label) + ' within 5 Å of the partner">Predicted interface</th><th class="hide-sm">Shared complexes</th></tr></thead><tbody>' +
          all.slice(0, UI.ixAll ? all.length : 40).map(([k, sc]) => {
            const x = IX.prot[k];
            const q = x.p >= 0 ? P[x.p] : null;
            const shared = (IX.cxOf.get(k) || []).filter((c) => p.cx.includes(c));
            return '<tr' + (q ? ' class="row" tabindex="0" data-i="' + q.i + '"' : '') + '><td><b>' + esc(x.gene) + '</b> <span class="muted mono">' + esc(x.acc) + '</span></td><td class="r num">' + (sc / 1000).toFixed(3) + '</td>' +
              '<td>' + (q ? '<span class="pill"><i style="background:' + (colors[mode.cat(q)] || pal.none) + '"></i>' + q.sites.length + '</span>' : '<a href="https://www.uniprot.org/uniprotkb/' + encodeURIComponent(x.acc) + '/entry" target="_blank" rel="noopener" class="muted">not in set ↗</a>') + '</td>' +
              '<td class="nw">' + (partnerTags(x) || '<span class="muted">—</span>') + '</td>' +
              '<td class="hide-sm">' + (predMap.has(k) ? (predMap.get(k).model ? '<button class="linkish" type="button" data-pvk="' + k + '" title="Open the ProtVar model">pDockQ ' + (predMap.get(k).sc / 100).toFixed(2) + '</button>' : 'pDockQ ' + (predMap.get(k).sc / 100).toFixed(2)) + ' <span class="muted">' + predMap.get(k).sites.slice(0, 4).map((x) => '' + rl(x)).join(', ') + (predMap.get(k).sites.length > 4 ? ' +' + (predMap.get(k).sites.length - 4) : '') + '</span>' : '<span class="muted">—</span>') + '</td>' +
              '<td class="hide-sm">' + (shared.length ? shared.slice(0, 2).map((c) => esc(IX.cx[c].name)).join('; ') + (shared.length > 2 ? ' +' + (shared.length - 2) : '') : '<span class="muted">—</span>') + '</td></tr>';
          }).join('') + '</tbody></table></div>' + (all.length > 40 ? '<button class="linkish" type="button" id="ix-all">' + (UI.ixAll ? 'Show the first 40' : 'Show all ' + fmt(all.length)) + '</button>' : '') : '') +
        '<p class="note">STRING v' + esc(IX.meta.string_version) + ' physical subnetwork (experiments, curated databases and text mining of physical binding). RNA*/DNA* mark GO annotations without experimental evidence. Click a methylated partner to open it.</p>';
      const lst = $('#ix-list', box);
      if (lst) lst.addEventListener('click', () => { F.ixScore = Math.max(700, min); replaceToken(partnerToken(), 'binds', IX.prot[p.i].gene, 'proteins'); });
      const allBtn = $('#ix-all', box);
      if (allBtn) allBtn.addEventListener('click', () => { UI.ixAll = !UI.ixAll; renderIxProteins(p, body); });
      $$('[data-pvk]', box).forEach((b) => b.addEventListener('click', (e) => { e.stopPropagation(); viewProtvar(p, predMap.get(+b.dataset.pvk), null); }));
      bindRows(box, (tr) => openProtein(+tr.dataset.i));
      if (!all.length) return;
      const draw = (net) => {
        const holder = $('#ix-net', box);
        if (!holder) return;
        holder.innerHTML = ixNetSvg(net, screenPalette(), false);
        $$('.nd', holder).forEach((nd) => {
          const k = +nd.dataset.k;
          nd.addEventListener('mouseenter', (e) => showTip(e, nd.getAttribute('data-tip')));
          nd.addEventListener('mousemove', moveTip);
          nd.addEventListener('mouseleave', hideTip);
          const open = () => { hideTip(); const x = IX.prot[k]; if (x.p >= 0) { if (x.p !== p.i) openProtein(x.p); } else window.open('https://www.uniprot.org/uniprotkb/' + encodeURIComponent(x.acc) + '/entry', '_blank', 'noopener'); };
          nd.addEventListener('click', open);
          nd.addEventListener('keydown', (e) => { if (e.key === 'Enter') open(); });
        });
        const note = $('#ix-net-note', box);
        if (note) note.textContent = 'Top ' + (net.nCore - 1) + ' partners by STRING score. ' + (net.live ? 'Thin grey lines join partners that are STRING physical partners of each other.' : 'Links among partners come from the STRING API; showing direct links only.') +
          (net.ring.length ? ' Dashed lines: PRMTs and Tudor-domain proteins (outer ring, italic) that bind those partners.' : '');
        const paths = $('#ix-paths', box);
        if (paths) {
          paths.innerHTML = net.ring.length ? '<p class="flabel">Writer and reader paths</p><ul class="paths">' + net.ring.map((r) => '<li><b>' + esc(IX.prot[r.k].gene) + '</b> <span class="tag ' + (IX.prot[r.k].flags & IXF.prmt ? 'ix-prmt">PRMT' : 'ix-tudor">Tudor') + '</span> binds ' +
            r.hits.map((h) => esc(IX.prot[net.ids[h]].gene)).join(', ') + ' <span class="muted">(partners of ' + esc(p.label) + ')</span></li>').join('') + '</ul>' :
            (UI.ixRing ? '<p class="note">No PRMT or Tudor-domain protein outside this network binds these partners at STRING ≥ 0.7.</p>' : '');
        }
        UI.ixNet = net;
      };
      draw(ixNetModel(p, all, min, null));
      const ids = [p.i].concat(all.slice(0, UI.ixTop).map(([k]) => k));
      loadStringNet(ids, min).then((extra) => { if (extra && UI.detail === p.i && UI.ixTab === 'proteins' && UI.ixMin === minAsked && document.body.contains(box)) draw(ixNetModel(p, all, min, extra)); });
      const name = (p.gene || p.acc) + '_physical_partners';
      $('#ix-svg').onclick = () => { if (UI.ixNet) exportFigure('svg', ixNetSvg(UI.ixNet, EXPORT_PAL, true), name); };
      $('#ix-png').onclick = () => { if (UI.ixNet) exportFigure('png', ixNetSvg(UI.ixNet, EXPORT_PAL, true), name); };
      $('#ix-tsv').onclick = () => download(name + '.tsv', ['partner_gene\tpartner_accession\tstring_physical_score\tmethylarginine_sites\tprmt\ttudor_domain\tgo_rna_binding\tgo_dna_binding\tprotvar_pdockq\tmethylarginines_at_predicted_interface']
        .concat(all.map(([k, sc]) => { const x = IX.prot[k]; const r = predMap.get(k); return [x.gene, x.acc, (sc / 1000).toFixed(3), x.p >= 0 ? P[x.p].sites.length : 0, x.flags & IXF.prmt ? 1 : 0, x.flags & IXF.tudor ? 1 : 0, x.flags & IXF.rna ? (x.flags & IXF.rnaExp ? 'experimental' : 'yes') : '', x.flags & IXF.dna ? (x.flags & IXF.dnaExp ? 'experimental' : 'yes') : '', r ? (r.sc / 100).toFixed(2) : '', r ? r.sites.map((v) => rl(v, p)).join(';') : ''].join('\t'); })).join('\n') + '\n');
    });
  }

  function memberChips(c, focal) {
    const slots = c.members.map((m) => (Array.isArray(m) ? m : [m]));
    const shown = slots.slice(0, 48);
    const chip = (ids) => {
      const meth = ids.find((k) => k < IX.nP);
      const k = meth != null ? meth : ids[0];
      const label = ids.map(ixGene).join('/');
      return meth != null ? '<button class="tag cx-m' + (meth === focal ? ' self' : '') + '" type="button" data-i="' + meth + '" title="' + esc(P[meth].sites.length + ' methylarginine sites') + '">' + esc(label) + '</button>' :
        '<span class="tag cx-u" title="' + esc(IX.prot[k] ? IX.prot[k].acc : '') + ' · not in this set">' + esc(label) + '</span>';
    };
    return '<div class="tags cx-members">' + shown.map(chip).join('') +
      (slots.length > shown.length ? '<details class="inline-more"><summary>+' + (slots.length - shown.length) + ' more</summary><div class="tags">' + slots.slice(shown.length).map(chip).join('') + '</div></details>' : '') + '</div>';
  }

  function cxLink(c) {
    return '<a href="https://www.ebi.ac.uk/complexportal/complex/' + encodeURIComponent(c.id) + '" target="_blank" rel="noopener">' + esc(c.id) + ' ↗</a>';
  }

  function renderIxComplexes(p, body) {
    const exp = $('#ix-exp');
    if (exp) exp.innerHTML = '';
    const list = p.cx.map((k) => IX.cx[k]).sort((a, b) => (a.pred - b.pred) || (a.q || 1) - (b.q || 1));
    if (!list.length) { body.innerHTML = '<p class="note">' + esc(p.label) + ' is not a member of any curated or predicted Complex Portal complex.</p>'; return; }
    const cur = list.filter((c) => !c.pred), pred = list.filter((c) => c.pred);
    const block = (c) => '<li><div class="cx-head"><b>' + esc(c.name) + '</b> <span class="tag' + (c.pred ? '' : ' cx-cur') + '">' + (c.pred ? 'predicted (hu.MAP 3.0)' : 'curated') + '</span> ' + cxLink(c) +
      '<span class="cx-stat" title="About the complex: members that are methylated in this dataset">' + c.m + ' of ' + c.n + ' members methylated' +
      (UI.ctx === 'dataset' && c.q != null ? ' · ' + (c.fold != null ? c.fold.toFixed(1) + '× the proteome rate, ' : '') + 'q ' + fmtP(c.q) : '') + '</span>' +
      '<button class="linkish" type="button" data-cx="' + c.k + '">Filter to members</button></div>' + memberChips(c, p.i) + '</li>';
    body.innerHTML = '<p class="sub">Complex Portal complexes that contain ' + esc(p.label) + ' (curated, or predicted by hu.MAP 3.0). Members that are methylated in this dataset are filled; click one to open it. Complex-level enrichment is ranked on the Overview.</p>' +
      (cur.length ? '<p class="flabel">Curated complexes (' + cur.length + ')</p><ul class="cx-list">' + cur.map(block).join('') + '</ul>' : '') +
      (pred.length ? '<details' + (cur.length ? '' : ' open') + '><summary class="flabel">Predicted complexes, hu.MAP 3.0 via Complex Portal (' + pred.length + ')</summary><ul class="cx-list">' + pred.map(block).join('') + '</ul></details>' : '');
    $$('.cx-m', body).forEach((b) => b.addEventListener('click', () => { if (+b.dataset.i !== p.i) openProtein(+b.dataset.i); }));
    $$('[data-cx]', body).forEach((b) => b.addEventListener('click', () => setTerm('complex', IX.cx[+b.dataset.cx].id, 'proteins')));
  }

  // Crosslinking evidence for one protein: crosslinked residues, methylarginines at or near them, regions.
  function xlSectionHtml(p, kind) {
    const xp = IX.xlProt.get(p.i);
    const sites = p.sites.map((j) => S[j]);
    const list = (arr) => arr.map((s) => '<button class="linkish strong" type="button" data-pos="' + s.pos + '">' + rl(s.pos, P[s.p]) + '</button>').join(', ');
    const ranges = (arr) => arr.slice(0, 8).map(([a, b]) => a + '–' + b).join(', ') + (arr.length > 8 ? ' +' + (arr.length - 8) : '');
    if (kind === 'rna') {
      const at = sites.filter((s) => s.xl && s.xl.r), near = sites.filter((s) => s.xl && !s.xl.r && s.xl.n), frag = sites.filter((s) => s.xl && s.xl.f);
      const sets = xp ? xlSets(xp.rnaBits.reduce((a, b) => a | b, 0)) : [];
      return '<section><h4>RNA crosslinks in living cells</h4>' + (xp && (xp.rnaPos.length || xp.rbd.length) ?
        '<p>' + (xp.rnaPos.length ? fmt(xp.rnaPos.length) + ' RNA-crosslinked residues (' + esc(sets.join(', ')) + ')' : 'No residue-level crosslinks') + (xp.rbd.length ? '; RBDmap RNA-bound fragments ' + ranges(xp.rbd) : '') + '.</p>' +
        '<ul class="ix-sites">' + (at.length ? '<li><b>At a crosslinked residue:</b> ' + list(at) + '</li>' : '') + (near.length ? '<li><b>Within 3 residues:</b> ' + list(near) + '</li>' : '') + (frag.length ? '<li><b>Inside an RBDmap fragment:</b> ' + list(frag) + '</li>' : '') +
        (!at.length && !near.length && !frag.length ? '<li class="muted">No methylarginine at or near a crosslinked residue.</li>' : '') + '</ul>' :
        '<p class="muted">No RNA crosslinks recorded for this protein in the pepR-MS compendium or RBDmap.</p>') +
        '<p class="note">UV crosslinking mass spectrometry localises the amino acid that touched RNA in living cells (pepR-MS re-analysis of Kramer 2014, XRNAX, RBS-ID, pRBS-ID and pepR-MS; RBDmap gives RNA-bound fragments). None of the pepR-MS crosslinked peptides carries a methylation mass shift, so a methylated arginine is seen there through its unmethylated fraction.</p></section>';
    }
    const at = sites.filter((s) => s.xl && s.xl.d & 1), cand = sites.filter((s) => s.xl && !(s.xl.d & 1) && s.xl.d & 8), pep = sites.filter((s) => s.xl && !(s.xl.d & 9) && s.xl.d & 2), bind = sites.filter((s) => s.xl && s.xl.d & 4);
    return '<section><h4>DNA contact in living cells</h4>' +
      '<p>' + (xp && xp.dnaInt ? '<b>In the DNA interactome</b> of MCF7 cells (direct protein–DNA photo-crosslinking; Trendel et al. 2025).' : 'Not in the Trendel et al. 2025 DNA interactome (MCF7).') +
      (xp && (xp.dnaPos.length || xp.dnaCand.length) ? ' Crosslinks localised to ' + fmt(xp.dnaPos.length) + ' residue' + (xp.dnaPos.length === 1 ? '' : 's') + ' alone and to ' + fmt(xp.dnaCand.length) + ' further candidate residue' + (xp.dnaCand.length === 1 ? '' : 's') + '.' : '') + (xp && xp.dnaBind.length ? ' UniProt DNA-binding regions: ' + ranges(xp.dnaBind) + '.' : '') + '</p>' +
      '<ul class="ix-sites">' + (at.length ? '<li><b>At a DNA-crosslinked residue:</b> ' + list(at) + '</li>' : '') + (cand.length ? '<li><b>A candidate residue of a DNA crosslink:</b> ' + list(cand) + '</li>' : '') + (pep.length ? '<li><b>Inside a DNA-crosslinked peptide:</b> ' + list(pep) + '</li>' : '') + (bind.length ? '<li><b>Inside a UniProt DNA-binding region:</b> ' + list(bind) + '</li>' : '') +
      (!at.length && !cand.length && !pep.length && !bind.length ? '<li class="muted">No methylarginine at a DNA-crosslinked residue or in a DNA-binding region.</li>' : '') + '</ul>' +
      '<p class="note">Most DNA crosslinks in Trendel et al. 2025 are localised only to a short run of candidate residues, so "candidate" marks residues the crosslink may sit on.</p></section>';
  }

  function renderIxNucleic(p, body, kind) {
    const exp = $('#ix-exp');
    if (exp) exp.innerHTML = '';
    const KIND = kind.toUpperCase();
    const bit = kind === 'rna' ? 4 : 8;
    const flag = kind === 'rna' ? IXF.rna : IXF.dna, flagExp = kind === 'rna' ? IXF.rnaExp : IXF.dnaExp;
    const term = kind === 'rna' ? 'GO:0003723' : 'GO:0003677';
    const sites = p.sites.map((j) => S[j]).filter((s) => s.iface & bit);
    const gene = p.gene;
    body.innerHTML =
      '<div class="ix-grid">' +
      '<section><h4>Gene Ontology</h4><p>' + (p.ixFlags & flag ? '<b>' + KIND + ' binding</b> (' + term + ' or a child term), ' + (p.ixFlags & flagExp ? 'with experimental evidence.' : 'from electronic or inferred annotation only.') : 'No ' + KIND + '-binding annotation (' + term + ').') + '</p><div id="ix-go" class="note">Loading QuickGO evidence…</div></section>' +
      '<section><h4>Methylarginines at ' + KIND + ' interfaces</h4>' + (sites.length ?
        '<ul class="ix-sites">' + sites.map((s) => { const g = contactGroups(s).find((x) => x.k === (kind === 'rna' ? -1 : -2)); return '<li><button class="linkish strong" type="button" data-pos="' + s.pos + '">' + rl(s.pos, P[s.p]) + '</button> ' + pill(s.st) +
          (g ? ' contacts <b>' + esc(g.labels.slice(0, 8).join(' ')) + '</b> <span class="muted">(' + g.pdbs.map((x) => x.toUpperCase()).join(', ') + ')</span>' : '') +
          ' <span class="muted small">at the interface in ' + (s.ifN || (s.ifPdb || []).length) + ' entr' + ((s.ifN || 0) === 1 ? 'y' : 'ies') + '</span></li>'; }).join('') + '</ul>' :
        '<p class="muted">None of the methylated arginines contacts ' + KIND + ' in a PDB structure (PDBe-KB).</p>') +
      '<p class="note">PDBe-KB interface residues across all PDB entries of this protein; methylation of these arginines would change the contact.</p></section>' +
      xlSectionHtml(p, kind) +
      '<section><h4>PDB structures with ' + KIND + '</h4><div id="ix-pdb" class="note">Looking up PDB entries…</div></section>' +
      (kind === 'rna' ? '<section><h4>ENCODE eCLIP</h4><div id="ix-encode" class="note">Searching ENCODE…</div></section>' :
        '<section><h4>ChIP-Atlas target genes</h4><div id="ix-chip" class="note">Loading ChIP-Atlas…</div></section><section><h4>ENCODE TF ChIP-seq</h4><div id="ix-encode" class="note">Searching ENCODE…</div></section>') +
      '</div>';
    $$('.ix-sites [data-pos]', body).forEach((b) => b.addEventListener('click', () => { UI.highlight = +b.dataset.pos; flashSite(+b.dataset.pos, false); }));
    const alive = () => UI.detail === p.i && UI.ixTab === kind && document.body.contains(body);
    loadQuickGO(p.acc.split('-')[0], term).then((rows) => {
      const el = $('#ix-go', body);
      if (!el || !alive()) return;
      if (!rows) { el.textContent = 'QuickGO did not answer.'; return; }
      if (!rows.length) { el.textContent = 'No QuickGO annotations.'; return; }
      const byCode = new Map();
      for (const r of rows) {
        const key = r.goEvidence;
        if (!byCode.has(key)) byCode.set(key, { terms: new Set(), refs: new Set() });
        byCode.get(key).terms.add(r.goName || r.goId);
        if (r.reference) byCode.get(key).refs.add(r.reference);
      }
      el.className = '';
      el.innerHTML = '<ul class="go-list">' + Array.from(byCode.entries()).sort((a, b) => (IX.meta.go_experimental_codes.includes(b[0]) - IX.meta.go_experimental_codes.includes(a[0]))).map(([code, v]) =>
        '<li><span class="tag' + (IX.meta.go_experimental_codes.includes(code) ? ' go-exp' : '') + '">' + esc(code) + '</span> ' + esc(Array.from(v.terms).join(', ')) +
        ' <span class="muted">' + (() => {
          const refs = Array.from(v.refs);
          const fmtRef = (ref) => (/^PMID:/.test(ref) ? '<a href="' + pubmedUrl(ref.slice(5)) + '" target="_blank" rel="noopener">' + esc(ref) + '</a>' : esc(ref));
          return refs.slice(0, 4).map(fmtRef).join(', ') + (refs.length > 4 ? ' <details class="inline-more"><summary>+' + (refs.length - 4) + ' more</summary>' + refs.slice(4).map(fmtRef).join(', ') + '</details>' : '');
        })() + '</span></li>').join('') + '</ul>' +
        '<p class="note"><a href="https://www.ebi.ac.uk/QuickGO/annotations?geneProductId=' + encodeURIComponent(p.acc.split('-')[0]) + '&goId=' + term + '&goUsage=descendants" target="_blank" rel="noopener">All annotations in QuickGO ↗</a></p>';
    });
    loadPdb(p.acc).then((res) => {
      const el = $('#ix-pdb', body);
      if (!el || !alive()) return;
      const hits = res.entries.filter((e) => res.na && res.na.get(e.id) && res.na.get(e.id)[kind]);
      if (res.naFailed) { el.textContent = 'The RCSB search did not answer.'; return; }
      if (!hits.length) { el.textContent = 'No PDB entry of this protein contains ' + KIND + '.'; return; }
      el.className = '';
      el.innerHTML = '<div class="tags">' + hits.slice(0, 40).map((e) => '<button class="tag pdb-btn" type="button" data-id="' + e.id + '" title="' + esc(methodShort(e.method)) + (e.resolution ? ', ' + e.resolution.toFixed(1) + ' Å' : '') + '">' + e.id.toUpperCase() + '</button>').join('') + (hits.length > 40 ? '<span class="tag">+' + (hits.length - 40) + '</span>' : '') + '</div>' +
        '<p class="note">' + hits.length + ' entr' + (hits.length === 1 ? 'y' : 'ies') + '. Click one to view it in the 3D panel, with ' + KIND + ' coloured and methylarginine contacts listed.</p>';
      $$('.pdb-btn', el).forEach((b) => b.addEventListener('click', () => {
        const e = res.entries.find((x) => x.id === b.dataset.id);
        const card = $('#card-af');
        if (card) card.scrollIntoView({ block: 'start', behavior: 'smooth' });
        if (UI.startStructure) UI.startStructure();
        showPdb(p, e, null);
      }));
    });
    loadEncode(gene, kind === 'rna' ? 'eCLIP' : 'TF ChIP-seq').then((rows) => {
      const el = $('#ix-encode', body);
      if (!el || !alive()) return;
      if (!rows) { el.textContent = 'ENCODE did not answer.'; return; }
      if (!rows.length) { el.textContent = 'No released ENCODE ' + (kind === 'rna' ? 'eCLIP' : 'TF ChIP-seq') + ' experiment targets ' + gene + '.'; return; }
      el.className = '';
      el.innerHTML = '<ul class="enc-list">' + rows.map((r) => '<li><a href="https://www.encodeproject.org/experiments/' + encodeURIComponent(r.accession) + '/" target="_blank" rel="noopener">' + esc(r.accession) + ' ↗</a> ' + esc((r.biosample_summary || '').replace(/^Homo sapiens /, '')) + ' <span class="muted">' + esc(r.date_released || '') + '</span></li>').join('') + '</ul>' +
        '<p class="note">Released ENCODE experiments; peaks and signal tracks are on each experiment page.</p>';
    });
    if (kind === 'dna') {
      loadChipAtlas(gene).then((d) => {
        const el = $('#ix-chip', body);
        if (!el || !alive()) return;
        if (!d) { el.textContent = 'ChIP-Atlas did not answer.'; return; }
        if (!d.genes.length) { el.textContent = 'ChIP-Atlas has no ChIP-seq for ' + gene + ' (hg38).'; return; }
        el.className = '';
        el.innerHTML = '<p><b>' + fmt(d.genes.length) + '</b> genes with a ' + esc(gene) + ' peak within 1 kb of the TSS, from ' + d.experiments + ' ChIP-seq experiment' + (d.experiments === 1 ? '' : 's') + (d.cells && d.cells.length ? ' (' + esc(d.cells.slice(0, 6).join(', ').replace(/_/g, ' ')) + (d.cells.length > 6 ? ', …' : '') + ')' : '') + '. Top genes by mean MACS2 score:</p>' +
          '<div class="tags">' + d.genes.slice(0, 20).map((g) => '<span class="tag" title="mean score ' + g[1].toFixed(0) + '">' + esc(g[0]) + '</span>').join('') + '</div>' +
          '<p class="note"><a href="https://chip-atlas.dbcls.jp/data/hg38/target/' + encodeURIComponent(gene) + '.1.tsv" target="_blank" rel="noopener">Full target-gene table ↗</a> · ChIP-Atlas, CC BY 4.0</p>';
      });
    }
  }

  // ----- overview cards
  // ----- Interactions view: key findings first, then one collapsible section per analysis.
  const IX_SECTIONS = [
    ['cells', 'RNA and DNA contact in living cells', 'Crosslinking mass spectrometry: is the methylated arginine itself the residue that touches RNA or DNA?'],
    ['pdb', 'Interfaces in PDB structures', 'Protein, RNA and DNA interfaces in experimental structures, within-protein tests and interface geometry'],
    ['predicted', 'Predicted complexes', 'AlphaFold-Multimer interfaces and the best-supported predicted contacts'],
    ['partners', 'Writers, readers and other partners', 'PRMTs and Tudor, WD40 and PHD proteins that physically bind methylated proteins'],
    ['complexes', 'Protein complexes', 'Complex Portal complexes enriched in methylated members'],
    ['go', 'RNA- and DNA-binding proteins (Gene Ontology)', 'Whole proteins: methylated vs other reviewed human proteins'],
  ];

  // Headline statistics, worded from the numbers so the text stays true if the data change.
  function ixFindings() {
    if (!ixReady()) return [];
    const out = [];
    const or2 = (l2) => Math.pow(2, l2);
    const orNum = (v) => (v >= 10 ? v.toFixed(0) : v.toFixed(2));
    const dirCi = (lo, hi) => (lo > 1 ? 'more often' : hi < 1 ? 'less often' : 'not significantly more or less often');
    const dirP = (l2, p) => (p <= 0.05 ? (l2 > 0 ? 'more often' : 'less often') : 'not significantly more or less often');
    const orCi = (o) => 'OR ' + fmtOr(o[0], o[1], o[2]) + ', p ' + fmtP(o[3]);
    const xe = IX.xlEnrich || {}, T = IX.xlTypes || {};
    if (xe.rna_xl && xe.rna_xl.or) {
      const [a, b, c, d] = xe.rna_xl.or.pooled, m = xe.rna_xl.or.mh;
      const one = xe.rna_xl_single && xe.rna_xl_single.or ? xe.rna_xl_single.or.mh : null;
      out.push({ sec: 'cells', html: 'In living cells, methylated arginines are <b>' + dirCi(m[1], m[2]) + '</b> the residue crosslinked to RNA than the other arginines of the same proteins: ' +
        pctText(a, a + b || 1) + '% vs ' + pctText(c, c + d || 1) + '% (within-protein ' + orCi(m) + ', ' + fmt(m[4]) + ' proteins' + (one ? '; crosslinks localised to a single residue only: ' + orCi(one) : '') + '). ' +
        'The crosslink searches carry no methylation mass, so these crosslinks were caught on unmethylated copies of the peptides.' });
    }
    const t = T.rna_xl || {}, t1 = T.rna_xl_single || {};
    if (t.rme2s && t.rme2a && t.nostate) {
      const part = (o, label) => label + ', <b>' + (o.mh[1] > 1 ? 'more often' : o.mh[2] < 1 ? 'less often' : 'no significant difference') + '</b> (' + orCi(o.mh) + ', ' + fmt(o.pooled[0]) + ' crosslinked)';
      const s1 = t1.rme2s && t1.rme2s.mh[0] != null ? ' With single-localisation crosslinks only, SDMA-IP sites: ' + orCi(t1.rme2s.mh) + '.' : '';
      out.push({ sec: 'cells', html: 'RNA crosslinking by methylation class, against unmethylated arginines of the same proteins: ' + part(t.rme2s, 'sites recovered by the SDMA IP (Rme2s)') + '; ' + part(t.rme2a, 'ADMA IP (Rme2a)') + '; ' + part(t.nostate, 'sites without a state call') + '.' + s1 +
        ' Few residues carry these classes, so treat the class differences as suggestive.' });
    }
    if (xe.dna_cand && xe.dna_cand.or) {
      const one = xe.dna_xl && xe.dna_xl.or ? xe.dna_xl.or.pooled : null, m = xe.dna_cand.or.mh, bind = xe.dna_bind && xe.dna_bind.or ? xe.dna_bind.or.mh : null;
      out.push({ sec: 'cells', html: 'DNA crosslinks (Trendel et al. 2025, MCF7) are rarely localised to one residue: ' + (one ? fmt(one[0]) + ' methylarginine' + (one[0] === 1 ? ' carries one' : 's carry one') + ' localised to it alone. ' : '') +
        'Counting every candidate residue, methylated arginines are <b>' + dirCi(m[1], m[2]) + '</b> a candidate than other arginines of the same proteins (' + orCi(m) + ')' +
        (bind && bind[0] != null ? ', and they sit <b>' + dirCi(bind[1], bind[2]) + '</b> inside UniProt DNA-binding regions (' + orCi(bind) + ')' : '') + '.' });
    }
    const st = IX.strata || {};
    if (st.rna && st.protein) {
      const r = st.rna.all, pr = st.protein.all;
      out.push({ sec: 'pdb', html: 'In PDB structures, methylated arginines sit at RNA interfaces <b>' + dirP(r.mh[0], r.mh[1]) + '</b> than other arginines of the same protein (OR ' + orNum(or2(r.mh[0])) + ', p ' + fmtP(r.mh[1]) +
        '; pooled across proteins OR ' + orNum(or2(r.pooled[4])) + ') and at protein interfaces <b>' + dirP(pr.mh[0], pr.mh[1]) + '</b> (OR ' + orNum(or2(pr.mh[0])) + ', p ' + fmtP(pr.mh[1]) + ').' });
    }
    const pe = IX.predEnrich && IX.predEnrich['0.23'];
    if (pe && pe.all && pe.all.mh[0] != null) {
      const m = pe.all.mh;
      out.push({ sec: 'predicted', html: 'In AlphaFold-Multimer models of human protein pairs (ProtVar, pDockQ ≥ 0.23), methylated arginines sit at the predicted interface <b>' + dirP(m[0], m[1]) + '</b> than other arginines of the same proteins (OR ' + orNum(or2(m[0])) + ', p ' + fmtP(m[1]) + ').' });
    }
    if (IX.gallery && IX.gallery.length) {
      const rows = galRows(Object.assign({}, GAL_DEFAULT));
      out.push({ sec: 'predicted', html: '<b>' + fmt(new Set(rows.map((r) => r.j)).size) + '</b> methylarginines form a salt bridge or hydrogen bond with a partner in a confident predicted complex, with both residues at pLDDT ≥ ' + GAL_DEFAULT.pl + '.' });
    }
    let edgesMeth = 0, edgesAll = 0;
    for (const p of P) for (const [k, sc] of p.partners || []) if (sc >= 700) { edgesAll++; if (k < IX.nP) edgesMeth++; }
    if (edgesAll) {
      out.push({ sec: 'partners', html: 'Of the high-confidence physical partners of methylated proteins, <b>' + pctText(edgesMeth, edgesAll) + '%</b> are themselves methylated, against ' + pctText(IX.meta.methylated_in_proteome, IX.meta.proteome_size) + '% of human proteins (STRING ≥ 0.7).' });
    }
    const g = IX.goEnrich && IX.goEnrich.rna;
    if (g) out.push({ sec: 'go', html: '<b>' + pctText(g[0], g[1] || 1) + '%</b> of methylated proteins are annotated as RNA-binding, against ' + pctText(g[2], g[3] || 1) + '% of other reviewed human proteins (Gene Ontology; OR ' + orNum(or2(g[4])) + ', p ' + fmtP(g[5]) + ').' });
    return out;
  }

  function findingsList(list) {
    return '<ul class="findings">' + list.map((f) => '<li><span>' + f.html + '</span> <button class="linkish more-link" type="button" data-ixgo="' + f.sec + '" title="Open the analysis">Details ＋</button></li>').join('') + '</ul>';
  }

  // Short version on the overview, pointing to the Interactions tab.
  function ixOverviewHtml() {
    if (!IX) return '<div class="card wide ix-feature" id="ch-ix"><h3>Interactions and binding</h3><p class="note">Loading interaction data…</p></div>';
    if (IX.failed) return '';
    const f = ixFindings();
    const pick = ['cells', 'pdb', 'predicted', 'partners'].map((k) => f.find((x) => x.sec === k)).filter(Boolean);
    return '<div class="card wide ix-feature" id="ch-ix"><div class="card-head"><div><h3>Interactions and binding</h3><p class="sub">What methylated arginines touch: protein partners, complexes, PDB and predicted structures, and RNA or DNA crosslinks in living cells.</p></div>' +
      '<div class="exp"><button type="button" class="btn" data-ixgo="">Explore interactions →</button></div></div>' + checkSiteHtml('ov-check') + findingsList(pick) + '</div>';
  }
  // One protein drawn as the contacts cartoon, to show what each kind of evidence looks like along a sequence.
  // Default: the protein of the current selection with most experimental-contact sites, else HNRNPA1.
  // Of the entries a symbol names, the reviewed one with most sites (unreviewed TrEMBL twins come last).
  const pickReviewed = (hit) => (hit && hit.length ? hit.map((i) => P[i]).sort((a, b) => a.unrev - b.unrev || b.sites.length - a.sites.length)[0] : null);
  function ixExampleProtein() {
    if (UI.ixEx != null && P[UI.ixEx]) return P[UI.ixEx];
    let best = null, bn = 0;
    for (const p of R.prots) { const n = p.sites.filter((j) => R.shown.has(j) && ixStatus(S[j]) === 'exp').length; if (n > bn) { bn = n; best = p; } }
    return constraints().length && best ? best : pickReviewed(lookupExact('HNRNPA1')) || best;
  }
  function ixExampleHtml() {
    if (!ixReady()) return '';
    const p = ixExampleProtein();
    if (!p) return '';
    return '<div class="card" id="ix-example"><div class="card-head"><div><h3>What the evidence looks like on one protein</h3><p class="sub">Methylarginines coloured by interaction status above, on the same residue scale, RNA and DNA crosslinks in living cells, PDB interface residues and predicted-complex interfaces. Dashed lines and rings mark evidence at the methylarginine itself. Every protein page has this view.' + docLink('ix-status') + '</p></div>' +
      '<div class="exp" role="group" aria-label="Export"><button type="button" data-ixex="svg">SVG</button><button type="button" data-ixex="png">PNG</button></div></div>' +
      '<form class="ix-ex-form" autocomplete="off"><div class="motif-search"><input id="ix-ex-q" type="text" spellcheck="false" placeholder="Gene, e.g. HNRNPA1, FUS, BRD4" value="' + esc(p.label) + '" aria-label="Protein to draw"><button class="btn small" type="submit">Show</button>' +
      '<button class="btn quiet small" type="button" id="ix-ex-open">Open ' + esc(p.label) + ' →</button></div><p class="note" id="ix-ex-msg"></p></form><div class="cc-box" id="ix-ex-box"><p class="note">Loading ' + esc(p.label) + '…</p></div></div>';
  }
  function bindIxExample(el) {
    const card = $('#ix-example', el);
    if (!card) return;
    const p = ixExampleProtein(), base = p.acc.split('-')[0];
    const draw = () => {
      const box = $('#ix-ex-box', card);
      if (!box || UI.view !== 'interactions') return;
      box.innerHTML = contactCartoonSvg(p, UNIDATA.get(p.acc), p.iso ? null : IFDATA.get(base) || null, { pal: Object.assign(screenPalette(), { bg: cssVar('--surface') }) }).svg;
      $$('.cc-site', box).forEach((c) => c.addEventListener('click', () => openProtein(p.i, +c.dataset.pos)));
    };
    Promise.all([loadUniProt(p.acc), p.iso ? null : loadIface(p.acc).catch(() => null)]).then(draw, draw);
    $('form', card).addEventListener('submit', (e) => {
      e.preventDefault();
      const v = $('#ix-ex-q', card).value.trim(), hit = v ? pickReviewed(lookupExact(v)) : null;
      if (!hit) { $('#ix-ex-msg', card).textContent = 'No protein in the set matches “' + v + '”.'; return; }
      UI.ixEx = hit.i;
      card.outerHTML = ixExampleHtml();
      bindIxExample(el);
    });
    $('#ix-ex-open', card).addEventListener('click', () => openProtein(p.i));
    $$('[data-ixex]', card).forEach((b) => b.addEventListener('click', () => exportFigure(b.dataset.ixex, contactCartoonSvg(p, UNIDATA.get(p.acc), IFDATA.get(base) || null, { forExport: true }), p.label + '_contacts')));
    bindDocLinks(card);
  }
  // "Is this site in an interaction?" answered in place: a gene and site (FUS R216), or a gene for its counts.
  const checkSiteHtml = (id) => '<form class="check-site" data-check-site autocomplete="off"><label class="flabel" for="' + id + '">Is a site in an interaction? Check one</label>' +
    '<div class="motif-search"><input id="' + id + '" type="text" spellcheck="false" placeholder="e.g. FUS R216, HNRNPA1 R194, or a gene"><button class="btn small" type="submit">Check</button></div><div class="check-out" data-check-out aria-live="polite"></div></form>';
  function bindCheckSite(el) {
    $$('[data-check-site]', el).forEach((f) => {
      const input = $('input', f), out = $('[data-check-out]', f);
      f.addEventListener('submit', (e) => {
        e.preventDefault();
        const v = input.value.trim();
        if (!v) return;
        const q = parsePlain(v);
        if (!q.proteins || !q.proteins.size) { out.innerHTML = '<p class="note">No protein in the set matches “' + esc(v) + '”. Try a gene and site, e.g. FUS R216.</p>'; return; }
        const i = q.positions ? q.positions.keys().next().value : q.proteins.values().next().value, p = P[i];
        const pos = q.positions && q.positions.get(i) ? q.positions.get(i).values().next().value : null;
        const s = pos == null ? null : p.sites.map((j) => S[j]).find((x) => x.pos === pos);
        if (!ixReady()) { out.innerHTML = '<p class="note">Interaction data is still loading.</p>'; return; }
        if (pos != null && !s) {
          const near = p.sites.map((j) => S[j].pos).sort((a, b) => Math.abs(a - pos) - Math.abs(b - pos)).slice(0, 3).sort((a, b) => a - b);
          out.innerHTML = '<p class="note">' + esc(p.label) + ' ' + rl(pos, p) + ' is not a reported methylarginine. Nearest reported sites: ' + near.map((x) => '<button class="linkish" type="button" data-ck="' + x + '">' + rl(x) + '</button>').join(', ') + '.</p>';
        } else if (s) {
          const m = s.ev || 0;
          out.innerHTML = '<div class="check-res"><b>' + esc(p.label) + ' ' + rl(s.pos, P[s.p]) + '</b> ' + ixBadge(s, true) +
            '<span class="muted small">' + (m ? esc(EV_SETS.filter(([b]) => m & b).map(([, l]) => l).join(' · ')) : 'no contact evidence recorded') + '</span>' +
            '<button class="btn quiet small" type="button" data-ck-open="' + s.pos + '">Open site →</button></div>';
        } else {
          const n = { exp: 0, pred: 0, none: 0 };
          for (const j of p.sites) n[ixStatus(S[j])]++;
          out.innerHTML = '<div class="check-res"><b>' + esc(p.label) + '</b> <span class="ixs ixs-exp"><i></i>' + n.exp + ' experimental</span> <span class="ixs ixs-pred"><i></i>' + n.pred + ' predicted only</span> <span class="ixs ixs-none"><i></i>' + n.none + ' none</span>' +
            '<button class="btn quiet small" type="button" data-ck-open="">Open protein →</button></div>';
        }
        $$('[data-ck]', out).forEach((b) => b.addEventListener('click', () => { input.value = p.label + ' ' + rl(+b.dataset.ck, p); f.requestSubmit(); }));
        $$('[data-ck-open]', out).forEach((b) => b.addEventListener('click', () => {
          if (b.dataset.ckOpen === '') { UI.pIx = ''; openProtein(p.i); }
          else openProtein(p.i, +b.dataset.ckOpen);
        }));
      });
    });
  }

  // ------------------------------------------------- evidence map (UpSet)
  // How the interaction evidence combines on the current selection, counted per site or per protein, and the protein
  // context of each kind of site evidence. Totals, rows and cells set the evidence filters; the actions hand the
  // selection to the Figures tab or the insight report.
  function upsetData(items, maskOf, sets) {
    const combos = new Map(), tot = sets.map(() => 0);
    let none = 0;
    for (const it of items) {
      const m = maskOf(it) || 0;
      if (!m) { none++; continue; }
      combos.set(m, (combos.get(m) || 0) + 1);
      sets.forEach(([b], k) => { if (m & b) tot[k]++; });
    }
    return { rows: Array.from(combos.entries()).sort((a, b) => b[1] - a[1] || a[0] - b[0]), tot, none, n: items.length };
  }
  const evCounts = (sites) => upsetData(sites, (s) => s.ev, EV_SETS);

  function upsetTableHtml(sets, d, attr, noun, all) {
    const shown = all ? d.rows : d.rows.slice(0, 12);
    const max = Math.max(1, ...d.rows.map((r) => r[1]));
    const label = (m) => sets.filter(([b]) => m & b).map(([, l]) => l).join(' + ');
    return '<div class="tbl-wrap"><table class="tbl compact upset"><thead><tr>' + sets.map(([, l, sh, help]) => '<th class="ev-col" title="' + esc(l + ': ' + help) + '">' + esc(sh) + '</th>').join('') + '<th>' + esc(noun) + '</th></tr>' +
      '<tr class="ev-tot">' + sets.map(([b, l], k) => '<td><button type="button" class="linkish" data-' + attr + 'set="' + b + '" title="List the ' + fmt(d.tot[k]) + ' ' + noun.toLowerCase() + ' with: ' + esc(l) + '">' + fmt(d.tot[k]) + '</button></td>').join('') + '<td class="muted small">' + noun.toLowerCase() + ' with each type</td></tr></thead><tbody>' +
      shown.map(([m, n]) => '<tr class="row" tabindex="0" data-' + attr + 'combo="' + m + '" title="' + esc(label(m)) + ' only">' + sets.map(([b]) => '<td class="dot">' + (m & b ? '<i class="on"></i>' : '<i></i>') + '</td>').join('') +
        '<td><div class="meter-row"><span class="meter"><i style="width:' + (100 * n / max).toFixed(1) + '%"></i></span><span class="num">' + fmt(n) + '</span></div></td></tr>').join('') + '</tbody></table></div>';
  }

  // Share of the sites with each kind of site evidence that sit on proteins with each protein property.
  const CTX_SETS = () => PEV_SETS.filter(([b]) => b !== 64);
  function contextRows(sites) {
    const groups = [[0, 'All selected sites', (s) => true]].concat(EV_SETS.map(([b, l]) => [b, l, (s) => (s.ev || 0) & b]));
    return groups.map(([b, l, test]) => {
      const list = sites.filter(test);
      return { b, l, n: list.length, cells: CTX_SETS().map(([pb]) => list.filter((s) => (P[s.p].pev || 0) & pb).length) };
    }).filter((r) => r.n);
  }
  function contextTableHtml(sites, clickable) {
    const rows = contextRows(sites), ref = rows[0];
    const cell = (r, k) => {
      const v = r.cells[k], f = r.n ? v / r.n : 0, rf = ref && ref.n ? ref.cells[k] / ref.n : 0;
      const ratio = rf ? f / rf : null;
      const shade = Math.min(0.55, f * 0.6);
      const title = fmt(v) + ' of ' + fmt(r.n) + ' sites' + (r.b ? '; ' + (ratio != null ? ratio.toFixed(1) + '× the share among all selected sites' : '') : '');
      const inner = Math.round(100 * f) + '%' + (r.b && ratio != null && v >= 3 && (ratio >= 1.25 || ratio <= 0.8) ? '<span class="ctx-ratio ' + (ratio > 1 ? 'up' : 'down') + '">' + (ratio > 1 ? '▲' : '▼') + '</span>' : '');
      return '<td class="r num ctx-cell" style="background:rgba(14,116,144,' + shade.toFixed(2) + ')" title="' + esc(title) + '">' +
        (clickable && v ? '<button type="button" class="linkish" data-ctx="' + r.b + ':' + CTX_SETS()[k][0] + '">' + inner + '</button>' : inner) + '</td>';
    };
    return '<div class="tbl-wrap"><table class="tbl compact ctx-tbl"><thead><tr><th>Sites</th><th class="r">n</th>' + CTX_SETS().map(([, l, sh, help]) => '<th class="r" title="' + esc(l + ': ' + help) + '">' + esc(sh) + '</th>').join('') + '</tr></thead><tbody>' +
      rows.map((r) => '<tr' + (r.b ? '' : ' class="ctx-ref"') + '><td>' + esc(r.l) + '</td><td class="r num">' + fmt(r.n) + '</td>' + r.cells.map((_, k) => cell(r, k)).join('') + '</tr>').join('') + '</tbody></table></div>';
  }

  function ixEvidenceMapHtml() {
    const scoped = constraints().length > 0, mode = UI.evMode;
    const seg = '<div class="seg" role="group" aria-label="Count"><span>Count</span><button type="button" data-evmode="sites" aria-pressed="' + (mode === 'sites') + '">Sites</button><button type="button" data-evmode="proteins" aria-pressed="' + (mode === 'proteins') + '">Proteins</button></div>';
    const head = (sub) => '<div class="card ix-evmap" id="ix-evmap"><div class="card-head"><div><h3>How the evidence fits together</h3><p class="sub">' + sub + docLink('ix-evmap') + '</p></div>' +
      '<div class="exp"><button type="button" id="ev-report" title="Profile, evidence map, logo and protein comparison of this selection as one PDF">Insight report</button><button type="button" id="ev-tsv">TSV</button></div></div><div class="fig-controls">' + seg + '</div>';
    if (mode === 'proteins') {
      const prots = R.prots, d = upsetData(prots, (p) => p.pev, PEV_SETS);
      return head((scoped ? 'The ' + fmt(prots.length) + ' proteins in your selection' : 'All ' + fmt(prots.length) + ' methylated proteins') + ', by protein-level evidence. Click a column total for every protein with that property, or a row for the proteins with exactly that combination.') +
        upsetTableHtml(PEV_SETS, d, 'pev', 'Proteins', UI.evAllRows) +
        '<p class="note">' + fmt(d.none) + ' proteins (' + pctText(d.none, d.n || 1) + '%) have none of these. Partners are STRING physical links ≥ 0.7; condensates are CD-CODE memberships at any confidence.' +
        (d.rows.length > 12 ? ' <button class="linkish" type="button" id="ev-more">' + (UI.evAllRows ? 'Show the top 12 combinations' : 'Show all ' + d.rows.length + ' combinations') + '</button>' : '') + '</p></div>';
    }
    const sites = R.sites, d = evCounts(sites);
    const key = UI.evTop;
    const per = new Map();
    for (const s of sites) if ((s.ev || 0) & key) per.set(s.p, (per.get(s.p) || 0) + 1);
    const top = Array.from(per.entries()).sort((a, b) => b[1] - a[1] || P[a[0]].label.localeCompare(P[b[0]].label)).slice(0, 15);
    const keyLabel = EV_SETS.find(([b]) => b === key)[1];
    return head((scoped ? 'The ' + fmt(sites.length) + ' sites in your selection' : 'All ' + fmt(sites.length) + ' methylarginines') + ', by the interaction evidence each carries. Click a column total for every site with that evidence, or a row for the sites with exactly that combination.') +
      upsetTableHtml(EV_SETS, d, 'ev', 'Sites', UI.evAllRows) +
      '<p class="note">' + fmt(d.none) + ' sites (' + pctText(d.none, d.n || 1) + '%) carry none of these; most sit in proteins that were never crystallised, modelled with a partner or crosslinked, so no evidence is not evidence of no contact.' +
      (d.rows.length > 12 ? ' <button class="linkish" type="button" id="ev-more">' + (UI.evAllRows ? 'Show the top 12 combinations' : 'Show all ' + d.rows.length + ' combinations') + '</button>' : '') + '</p>' +
      '<h4 class="ev-h">Protein context of these sites</h4><p class="sub">Share of the sites in each row that sit on proteins with each property; ▲▼ mark rows at least 1.25× above or below the share among all selected sites. Click a cell to list those sites.</p>' +
      contextTableHtml(sites, true) +
      '<div class="ev-top"><div class="fig-controls"><label class="fig-ctl">Proteins with most sites carrying <select id="ev-top">' + EV_SETS.map(([b, l]) => '<option value="' + b + '"' + (b === key ? ' selected' : '') + '>' + esc(l) + '</option>').join('') + '</select></label></div>' +
      (top.length ? '<div class="gene-chips">' + top.map(([i, n]) => '<button class="gene-chip" type="button" data-i="' + i + '" title="' + n + ' of ' + P[i].sites.length + ' methylarginines on ' + esc(P[i].label) + ' carry: ' + esc(keyLabel) + '">' + esc(P[i].label) + '<span>' + n + '/' + P[i].sites.length + '</span></button>').join('') + '</div>' +
        '<div class="ev-actions"><button class="btn quiet small" type="button" data-evact="sites">List these ' + fmt(d.tot[EV_SETS.findIndex(([b]) => b === key)]) + ' sites</button><button class="btn quiet small" type="button" data-evact="logo">Logo and alignment of these sites</button><button class="btn quiet small" type="button" data-evact="compare">Compare the top ' + Math.min(8, top.length) + ' proteins</button></div>' :
        '<p class="note">No site in the selection carries this evidence.</p>') + '</div></div>';
  }

  function bindEvidenceMap(el) {
    const box = $('#ix-evmap', el);
    if (!box) return;
    const redraw = () => { box.outerHTML = ixEvidenceMapHtml(); bindEvidenceMap(el); };
    $$('[data-evmode]', box).forEach((b) => b.addEventListener('click', () => { UI.evMode = b.dataset.evmode; UI.evAllRows = false; redraw(); }));
    $$('[data-evset]', box).forEach((b) => b.addEventListener('click', () => setFilters((f) => { f.evAll = +b.dataset.evset; f.evExact = false; }, 'sites')));
    $$('[data-pevset]', box).forEach((b) => b.addEventListener('click', () => setFilters((f) => { f.pevAll = +b.dataset.pevset; f.pevExact = false; }, 'proteins')));
    bindRows(box, (tr) => {
      if (tr.dataset.evcombo) setFilters((f) => { f.evAll = +tr.dataset.evcombo; f.evExact = true; }, 'sites');
      else if (tr.dataset.pevcombo) setFilters((f) => { f.pevAll = +tr.dataset.pevcombo; f.pevExact = true; }, 'proteins');
    });
    $$('[data-ctx]', box).forEach((b) => b.addEventListener('click', () => {
      const [eb, pb] = b.dataset.ctx.split(':').map(Number);
      setFilters((f) => { f.evAll = eb; f.evExact = false; f.pevAll = pb; f.pevExact = false; }, 'sites');
      window.scrollTo({ top: 0 });
    }));
    const more = $('#ev-more', box);
    if (more) more.addEventListener('click', () => { UI.evAllRows = !UI.evAllRows; redraw(); });
    const sel = $('#ev-top', box);
    if (sel) sel.addEventListener('change', (e) => { UI.evTop = +e.target.value; redraw(); });
    $$('.gene-chip[data-i]', box).forEach((b) => b.addEventListener('click', () => openProtein(+b.dataset.i)));
    $$('[data-evact]', box).forEach((b) => b.addEventListener('click', () => {
      const key = UI.evTop, act = b.dataset.evact;
      if (act === 'compare') {
        const per = new Map();
        for (const s of R.sites) if ((s.ev || 0) & key) per.set(s.p, (per.get(s.p) || 0) + 1);
        const genes = Array.from(per.entries()).sort((a, c) => c[1] - a[1]).slice(0, 8).map(([i]) => P[i].label);
        setPlain(genes.join(' '));
      }
      setFilters((f) => { f.evAll = key; f.evExact = false; }, act === 'sites' ? 'sites' : 'figures');
      window.scrollTo({ top: 0 });
    }));
    $('#ev-report', box).addEventListener('click', () => openReport());
    $('#ev-tsv', box).addEventListener('click', () => {
      const prot = UI.evMode === 'proteins', sets = prot ? PEV_SETS : EV_SETS;
      const d = prot ? upsetData(R.prots, (p) => p.pev, PEV_SETS) : evCounts(R.sites);
      const lines = [sets.map(([, l]) => l.toLowerCase().replace(/[^a-z0-9]+/g, '_')).join('\t') + '\t' + (prot ? 'proteins' : 'sites')]
        .concat(d.rows.map(([m, n]) => sets.map(([b]) => (m & b ? 1 : 0)).join('\t') + '\t' + n));
      if (!prot) {
        lines.push('', '# protein context: share of the sites in each group on proteins with each property', 'site_group\tsites\t' + CTX_SETS().map(([, l]) => l.toLowerCase().replace(/[^a-z0-9]+/g, '_')).join('\t'));
        for (const r of contextRows(R.sites)) lines.push(r.l + '\t' + r.n + '\t' + r.cells.join('\t'));
      }
      download('methylarginine_evidence_' + (prot ? 'proteins' : 'sites') + '.tsv', lines.join('\n') + '\n');
    });
  }

  // ------------------------------------------------------------ insight report
  // One printable page for the current selection: what it is, what sets it apart, how its interaction evidence combines
  // (sites, their protein context, proteins), its sequence logo and alignment, its proteins on one scale, the proteins
  // in it and the methods. The browser's print dialog saves it as a PDF; the URL reproduces the selection.
  function openReport() { UI.reportPrint = true; go('report'); window.scrollTo({ top: 0 }); }

  // Running header and footer on every printed page of the report and of protein pages: the site and its address, the paper, and the page
  // number. CSS page-margin boxes (Chromium-based browsers) carry them; elsewhere a repeated fixed footer carries the
  // citation without page numbers. The text comes from the page itself, so a renamed or moved site prints correctly.
  const cssStr = (t) => '"' + String(t).replace(/["\\]/g, '\\$&').replace(/\n/g, ' ') + '"';
  const reportCite = () => 'Please cite DeAngelo, Silverstein & Shechter, Cell Chem Biol 33:1404–1421 (2026), doi:' + PAPER.doi;
  function setPrintPageStyle(on, label) {
    let st = document.getElementById('report-page-style');
    if (!on) { if (st) st.remove(); return; }
    if (!st) { st = document.createElement('style'); st.id = 'report-page-style'; document.head.appendChild(st); }
    const box = 'font-family: Arial, Helvetica, sans-serif; font-size: 7pt; color: #5f6f78;';
    st.textContent = '@media print { @page { margin: 16mm 11mm 15mm; ' +
      '@top-left { content: ' + cssStr(siteName() + ' · ' + siteUrl().replace(/^https?:\/\//, '')) + '; ' + box + ' } ' +
      '@top-right { content: ' + cssStr(label) + '; ' + box + ' } ' +
      '@bottom-left { content: ' + cssStr(reportCite()) + '; ' + box + ' } ' +
      '@bottom-right { content: "Page " counter(page) " of " counter(pages); ' + box + ' } } }';
  }

  function upsetSvg(sets, d, title, noun) {
    const pal = EXPORT_PAL, rows = d.rows.slice(0, 12), max = Math.max(1, ...rows.map((r) => r[1]));
    const colW = 78, x0 = 16, barX = x0 + sets.length * colW + 14, barW = 230, W = barX + barW + 64, top = 92, rh = 18;
    let g = '<text x="16" y="26" font-size="15" font-weight="700" fill="' + pal.ink + '">' + esc(title) + '</text>' +
      '<text x="16" y="44" font-size="10.5" fill="' + pal.muted + '">' + esc(fmt(d.n) + ' ' + noun + '; the ' + rows.length + ' most frequent combinations of evidence (dots), with the number of ' + noun + ' carrying exactly that combination') + '</text>';
    sets.forEach(([, , sh], k) => {
      const cx = x0 + k * colW + colW / 2;
      g += '<text x="' + cx + '" y="' + (top - 26) + '" font-size="9.5" font-weight="700" text-anchor="middle" fill="' + pal.ink + '">' + esc(sh) + '</text>' +
        '<text x="' + cx + '" y="' + (top - 12) + '" font-size="9" text-anchor="middle" fill="' + pal.muted + '">' + fmt(d.tot[k]) + '</text>';
    });
    rows.forEach(([m, n], r) => {
      const y = top + r * rh + rh / 2;
      if (r % 2 === 0) g += '<rect x="' + x0 + '" y="' + (y - rh / 2) + '" width="' + (W - 2 * x0) + '" height="' + rh + '" fill="' + pal.backboneFill + '"/>';
      const on = sets.map(([b], k) => (m & b ? k : -1)).filter((k) => k >= 0);
      if (on.length > 1) g += '<line x1="' + (x0 + on[0] * colW + colW / 2) + '" x2="' + (x0 + on[on.length - 1] * colW + colW / 2) + '" y1="' + y + '" y2="' + y + '" stroke="' + pal.ink + '" stroke-width="1.5"/>';
      sets.forEach(([b], k) => { g += '<circle cx="' + (x0 + k * colW + colW / 2) + '" cy="' + y + '" r="4.2" fill="' + (m & b ? pal.ink : pal.grid) + '"/>'; });
      g += '<rect x="' + barX + '" y="' + (y - 5) + '" width="' + r1(Math.max(1, (n / max) * barW)) + '" height="10" rx="2" fill="' + pal.bar + '"/><text x="' + (barX + (n / max) * barW + 5) + '" y="' + (y + 3.5) + '" font-size="9.5" fill="' + pal.ink + '">' + fmt(n) + '</text>';
    });
    const H = top + rows.length * rh + 24;
    g += '<text x="16" y="' + (H - 8) + '" font-size="9.5" fill="' + pal.muted + '">' + fmt(d.none) + ' ' + noun + ' (' + pctText(d.none, d.n || 1) + '%) carry none of these; absence of evidence is not evidence of no contact.</text>';
    return { W, H, svg: '<svg xmlns="http://www.w3.org/2000/svg" width="' + W + '" height="' + H + '" viewBox="0 0 ' + W + ' ' + H + '" font-family=\'' + pal.font.replace(/'/g, '"') + '\'><rect width="' + W + '" height="' + H + '" fill="' + pal.bg + '"/>' + g + '</svg>' };
  }

  function renderReport(el) {
    const sites = R.sites, prots = R.prots, scoped = constraints().length > 0;
    const date = new Date().toISOString().slice(0, 10);
    const link = siteUrl() + location.search + '#sites';
    const share = (n) => pctText(n, sites.length || 1) + '%';
    const nIdr = sites.filter((s) => s.dc === M.vi.disordered).length;
    const nEv = sites.filter((s) => s.ev).length;
    const nRnaProt = sites.filter((s) => (P[s.p].pev || 0) & 1).length;
    const nRep = sites.filter((s) => popcount(s.study) >= 2).length;
    const facts = [[fmt(sites.length), 'methylarginine sites'], [fmt(prots.length), 'proteins'], [share(nIdr), 'inside an IDR'],
      [ixReady() ? share(nEv) : '—', 'with interaction evidence at the site'], [ixReady() ? share(nRnaProt) : '—', 'on GO RNA-binding proteins'], [share(nRep), 'reported by ≥2 studies']];
    const prof = scoped ? profileRows() : null;
    const sig = prof ? prof.filter((r) => r.q <= 0.05) : [];
    const profUp = sig.filter((r) => r.fold > 1).slice(0, 12), profDown = sig.filter((r) => r.fold < 1).slice(0, 6);
    const pc = (x, n) => (100 * x / n).toFixed(x / n < 0.1 ? 1 : 0) + '%';
    const alignList = figSites();
    const cmpProts = figProteins().slice(0, 8);
    const topProts = prots.slice().sort((a, b) => R.counts[b.i] - R.counts[a.i] || a.label.localeCompare(b.label)).slice(0, 25);
    const evCount = (p, b) => p.sites.filter((j) => R.shown.has(j) && (S[j].ev || 0) & b).length;
    const sec = (id, title, body, cls) => '<section class="card report-sec' + (cls ? ' ' + cls : '') + '" id="rep-' + id + '"><h3>' + esc(title) + '</h3>' + body + '</section>';
    const marginBoxes = /Chrome\/|Chromium\/|Edg\//.test(navigator.userAgent);
    setPrintPageStyle(true, 'Insight report · ' + date);
    el.innerHTML = '<div class="report' + (marginBoxes ? '' : ' runfoot') + '">' +
      '<div class="print-only report-runfoot">' + esc(siteName() + ' · ' + siteUrl() + ' · ' + reportCite()) + '</div>' +
      '<div class="report-bar screen-only"><button class="linkish" type="button" id="rep-back">← Back</button><button class="btn" type="button" id="rep-print">Save as PDF</button>' + cmpButtons() + '<span class="muted small">Opens the print dialog; choose “Save as PDF”. The link below reopens this selection.</span></div>' +
      '<header class="report-head"><p class="report-kicker">' + esc(siteName()) + ' · insight report · ' + date + '</p><h2>' + esc(scoped ? selectionText() : 'The whole methylarginine union') + '</h2>' +
      '<p class="muted small">Reopen: <span class="mono">' + esc(link) + '</span><br>Data exported ' + esc(M.generated) + (IX && IX.meta ? '; interaction data ' + esc(IX.meta.generated || '') : '') + '. Please cite ' + esc(PAPER.cite) + ' doi:' + esc(PAPER.doi) + '</p></header>' +
      '<div class="facts report-facts">' + facts.map(([v, l]) => '<div class="fact"><b>' + v + '</b><span>' + esc(l) + '</span></div>').join('') + '</div>' +
      sec('summary', 'In words', '<p class="report-summary">' + esc(summaryText(currentSelection(), prof)) + '</p><p class="note">Written from the numbers in this report; it states differences, not causes.</p>') +
      (scoped ? sec('profile', 'What sets this selection apart', sig.length ?
        '<div class="tbl-wrap"><table class="tbl compact prof-tbl"><thead><tr><th>Feature</th><th class="r">Here</th><th class="r hide-xs">Other sites</th><th class="r">Ratio</th><th class="r">q</th></tr></thead><tbody>' +
        profUp.concat(profDown).map((r) => '<tr' + (r.fold < 1 ? ' class="depleted"' : '') + '><td><span class="muted small">' + esc(r.ft.group) + '</span> ' + esc(r.ft.label) + '</td><td class="r num">' + fmt(r.a) + ' (' + pc(r.a, r.n1) + ')</td><td class="r num hide-xs">' + pc(r.c, r.n2) + '</td><td class="r num"><b>' + (r.fold >= 10 ? r.fold.toFixed(0) : r.fold.toFixed(1)) + '×</b></td><td class="r num">' + fmtP(r.q) + '</td></tr>').join('') +
        '</tbody></table></div><p class="note">The ' + fmt(sites.length) + ' selected sites against the other ' + fmt(S.length - sites.length) + ': Fisher\'s exact test per feature, Benjamini–Hochberg across the ' + (prof ? prof.length : 0) + ' features tested; q ≤ 0.05 shown. Sites in one protein are not independent, so read this as a pointer.</p>' :
        '<p class="note">No site feature differs significantly (BH q ≤ 0.05) between this selection and the other sites.</p>') : '') +
      (ixReady() ? sec('evidence', 'Interaction evidence', '<div class="fig-preview">' + upsetSvg(EV_SETS, evCounts(sites), 'Site-level evidence', 'sites').svg + '</div>' +
        '<h4 class="ev-h">Protein context of these sites</h4><p class="sub">Share of the sites in each row on proteins with each property; ▲▼ mark rows at least 1.25× above or below the share among all selected sites.</p>' + contextTableHtml(sites, false) +
        '<div class="fig-preview">' + upsetSvg(PEV_SETS, upsetData(prots, (p) => p.pev, PEV_SETS), 'Protein-level evidence', 'proteins').svg + '</div>' +
        '', 'report-break') : '') +
      sec('logo', 'Sequence logo and alignment', alignList.length ? '<div class="fig-preview">' + alignSvg(alignList, { flank: 7, rows: 30, logo: 'freq', title: 'Selected sites aligned at the methylarginine', sub: fmt(alignList.length) + ' sites with a sequence window; the logo uses all of them' }).svg + '</div>' : '<p class="note">No selected site has a sequence window.</p>', 'report-break') +
      sec('compare', 'Proteins on one scale', '<div class="fig-preview" id="rep-cmp"><p class="note">Loading UniProt features for ' + fmt(cmpProts.length) + ' proteins…</p></div>', 'report-break') +
      sec('proteins', 'Proteins in the selection', '<div class="tbl-wrap"><table class="tbl compact rep-prot"><thead><tr><th>Protein</th><th class="r">Sites here / all</th>' + (ixReady() ? EV_SETS.slice(0, 6).map(([, , sh]) => '<th class="r">' + esc(sh) + '</th>').join('') + '<th>Protein evidence</th>' : '') + '</tr></thead><tbody>' +
        topProts.map((p) => '<tr><td><b>' + esc(p.label) + '</b> <span class="muted mono">' + esc(p.acc) + '</span><br><span class="muted small">' + esc((p.name || '').slice(0, 60)) + '</span></td><td class="r num">' + R.counts[p.i] + ' / ' + p.sites.length + '</td>' +
          (ixReady() ? EV_SETS.slice(0, 6).map(([b]) => { const n = evCount(p, b); return '<td class="r num">' + (n || '<span class="muted">·</span>') + '</td>'; }).join('') + '<td class="small">' + esc(PEV_SETS.filter(([b]) => b !== 64 && (p.pev || 0) & b).map(([, , sh]) => sh).join(', ')) + '</td>' : '') + '</tr>').join('') +
        '</tbody></table></div>' + (prots.length > topProts.length ? '<p class="note">The ' + topProts.length + ' proteins with most selected sites of ' + fmt(prots.length) + '.</p>' : '')) +
      sec('cite', 'Sources to cite', citeBoxHtml(sites, scoped ? 'these ' + fmt(sites.length) + ' sites' : 'the sites in the union', true)) +
      sec('methods', 'Methods and sources', '<ul class="small">' +
        '<li>Sites: the integrated human methylarginine union (Maron et al. 2021 Table S5, ProMetheusDB, UniProt and IEDB) on canonical UniProt sequences, as described in the Documentation of ' + esc(siteName()) + '.</li>' +
        '<li>Interaction evidence: pepR-MS RNA crosslinks and RBDmap fragments, Trendel et al. 2025 DNA crosslinks, PDBe-KB interfaces, ProtVar AlphaFold-Multimer interfaces (pDockQ ≥ 0.23) and confident predicted contacts; protein evidence from Gene Ontology, CD-CODE, STRING physical links (≥ 0.7) and Complex Portal.</li>' +
        '<li>Logo: residue frequencies at each position of the ±7 windows; colours follow WebLogo\'s chemistry scheme. Protein features and disordered regions: UniProt (live).</li>' +
        '<li>Feature comparison: Fisher\'s exact test, Benjamini–Hochberg; exploratory, since sites of one protein are not independent.</li>' +
        (ixReady() ? '<li>Site evidence: ' + esc(EV_SETS.map(([, , sh, h]) => sh + ' = ' + h).join('; ')) + '. Protein evidence: ' + esc(PEV_SETS.map(([, , sh, h]) => sh + ' = ' + h).join('; ')) + '.</li>' : '') + '</ul>' +
        '<p class="small">' + esc(creditLines().join(' · ')) + '</p>') +
      '</div>';
    $('#rep-back', el).addEventListener('click', () => { if (history.length > 1) history.back(); else go('sites'); });
    const print = () => {
      const t = document.title;
      document.title = siteName() + ' insight report ' + date;
      window.print();
      document.title = t;
    };
    $('#rep-print', el).addEventListener('click', print);
    bindCmpButtons(el);
    bindCiteBox(el);
    let printed = false;
    const ready = () => { if (UI.reportPrint && !printed && UI.view === 'report') { printed = true; UI.reportPrint = false; setTimeout(print, 300); } };
    const timer = setTimeout(ready, 12000);
    Promise.all(cmpProts.map((p) => loadUniProt(p.acc))).then((us) => {
      const box = $('#rep-cmp');
      if (!box || UI.view !== 'report') return;
      const items = cmpProts.map((p, k) => ({ p, u: us[k] }));
      box.innerHTML = items.length ? compareSvg(items, { scale: 'common', color: 'state', title: 'Methylarginine sites on ' + fmt(items.length) + ' protein' + (items.length === 1 ? '' : 's'), sub: (Q.proteins && Q.proteins.size ? 'In search order' : 'The proteins with most selected sites') + ' · common residue scale · UniProt domains and disordered regions' }).svg : '<p class="note">No protein selected.</p>';
      clearTimeout(timer);
      ready();
    });
  }

  // ------------------------------------------------ summaries and comparison
  // A selection as plain data: its sites, proteins and per-protein counts, from saved filters and search text.
  function selectionOf(fl, qText) {
    const q = parseQuery(qText || '');
    const { sites, counts } = filterSites(fl, q);
    return { sites, counts, prots: P.filter((p) => counts[p.i] > 0), shown: new Set(sites.map((s) => s.j)) };
  }
  const currentSelection = () => ({ sites: R.sites, counts: R.counts, prots: R.prots, shown: R.shown });
  const listAnd = (a) => (a.length < 2 ? a.join('') : a.slice(0, -1).join(', ') + ' and ' + a[a.length - 1]);
  // A feature label inside a sentence: ordinary leading words lower-cased, motifs named as motifs.
  const lcLabel = (ft) => (ft.group === 'Motif' ? 'the ' + ft.label.split(' · ')[0] + ' motif'
    : /^(Decreased|Increased|Recovered|Inside|In|Another|Reported|Predicted|Experimental|Confident)\b/.test(ft.label) ? lc1(ft.label) : ft.label);
  // A paragraph worded from the numbers of a selection only (no interpretation beyond them).
  function summaryText(sel, prof) {
    const nAll = S.length, n = sel.sites.length, whole = n === nAll;
    const pc = (a, b) => pctText(a, b || 1) + '%';
    const share = (test) => { let a = 0, b = 0; for (const s of sel.sites) if (test(s)) a++; for (const s of S) if (test(s)) b++; return [a, b]; };
    if (!n) return 'No sites match this selection.';
    const out = [whole ? 'The union holds ' + fmt(nAll) + ' methylarginine sites on ' + fmt(P.length) + ' proteins.'
      : 'This selection holds ' + fmt(n) + ' methylarginine site' + (n === 1 ? '' : 's') + ' on ' + fmt(sel.prots.length) + ' protein' + (sel.prots.length === 1 ? '' : 's') + ', ' + pc(n, nAll) + ' of the ' + fmt(nAll) + ' sites in the union.'];
    if (ixReady()) {
      const [e, e0] = share((s) => ixStatus(s) === 'exp'), [pr, p0] = share((s) => ixStatus(s) === 'pred');
      out.push(fmt(e) + ' (' + pc(e, n) + ') have experimental contact evidence at the residue and ' + fmt(pr) + ' (' + pc(pr, n) + ') sit at a predicted interface only' +
        (whole ? '' : ', against ' + pc(e0, nAll) + ' and ' + pc(p0, nAll) + ' of all sites') + '.');
    }
    if (prof) {
      const sig = prof.filter((r) => r.q <= 0.05);
      const f = (r) => lcLabel(r.ft) + ' (' + pc(r.a, r.n1) + ' here vs ' + pc(r.c, r.n2) + ', q ' + fmtP(r.q) + ')';
      const up = sig.filter((r) => r.fold > 1).slice(0, 3), down = sig.filter((r) => r.fold < 1).slice(0, 2);
      if (up.length) out.push('Against the other sites, the clearest differences are ' + listAnd(up.map(f)) + '.');
      if (down.length) out.push('Less common here: ' + listAnd(down.map(f)) + '.');
      if (!sig.length) out.push('No site feature tested differs from the other sites at BH q ≤ 0.05.');
    }
    const [idr, idr0] = share((s) => s.dc === M.vi.disordered);
    out.push(pc(idr, n) + ' lie in intrinsically disordered regions' + (whole ? '' : ' (' + pc(idr0, nAll) + ' overall)') + '.');
    if (ixReady()) {
      const [rb, rb0] = share((s) => (P[s.p].pev || 0) & 1);
      out.push(pc(rb, n) + ' are on proteins with Gene Ontology RNA binding' + (whole ? '' : ' (' + pc(rb0, nAll) + ' overall)') + '.');
    }
    if (sel.prots.length > 1) {
      const top = sel.prots.slice().sort((a, b) => sel.counts[b.i] - sel.counts[a.i] || a.label.localeCompare(b.label)).slice(0, 3);
      out.push('The proteins with most sites are ' + listAnd(top.map((p) => p.label + ' (' + sel.counts[p.i] + ')')) + '.');
    }
    return out.join(' ');
  }

  // Selection A for the comparison view: the search and filters saved from an earlier selection.
  function saveSelectionA() {
    const qs = new URLSearchParams();
    if (Q.raw) qs.set('q', Q.raw);
    for (const [k, v] of filterParams(F)) qs.set(k, v);
    UI.cmpA = { qs: qs.toString(), label: selectionText() };
    syncUrl();
    renderScope();
  }
  function selectionA() {
    if (!UI.cmpA) return null;
    const params = new URLSearchParams(UI.cmpA.qs);
    const fl = filtersFromParams(params), qText = params.get('q') || '';
    if (!UI.cmpA.label) {
      // A link carries A's search and filters only: describe them the way the scope bar would.
      const keepF = F, keepQ = Q;
      F = fl; Q = parseQuery(qText);
      UI.cmpA.label = selectionText();
      F = keepF; Q = keepQ;
    }
    return Object.assign(selectionOf(fl, qText), { label: UI.cmpA.label });
  }
  // Features compared between A and B: share of each selection's sites, ratio, Fisher's exact test, BH across features.
  function compareRows(A, B) {
    const feats = profileFeatures();
    if (ixReady()) {
      feats.unshift({ key: 'ixp', group: 'Interactions', label: 'Predicted only (no experimental evidence)', test: (s) => ixStatus(s) === 'pred' });
      feats.unshift({ key: 'ixe', group: 'Interactions', label: 'Experimental contact evidence', test: (s) => ixStatus(s) === 'exp' });
    }
    const nA = A.sites.length, nB = B.sites.length;
    const rows = feats.map((ft) => {
      let a = 0, b = 0;
      for (const s of A.sites) if (ft.test(s)) a++;
      for (const s of B.sites) if (ft.test(s)) b++;
      return { ft, a, b, nA, nB, fold: (a / (nA || 1)) / ((b || 0.5) / (nB || 1)), p: fisherP(a, nA - a, b, nB - b) };
    }).filter((r) => r.a + r.b > 0);
    const order = rows.slice().sort((x, y) => x.p - y.p);
    let prev = 1;
    for (let k = order.length - 1; k >= 0; k--) { prev = Math.min(prev, (order[k].p * order.length) / (k + 1)); order[k].q = prev; }
    return order;
  }

  // Save as A / Compare with A, offered wherever a selection is summarised.
  function cmpButtons() {
    return '<span class="cmp-btns"><button class="btn quiet small" type="button" data-cmp="save" title="Keep this selection as A, then change the search or filters to compare another selection with it">' + (UI.cmpA ? 'Replace A with this' : 'Save as A for a comparison') + '</button>' +
      (UI.cmpA ? '<button class="btn small" type="button" data-cmp="open" title="A: ' + esc(UI.cmpA.label || '') + '">Compare with A</button>' : '') + '</span>';
  }
  function bindCmpButtons(root) {
    $$('[data-cmp]', root).forEach((b) => b.addEventListener('click', () => {
      if (b.dataset.cmp === 'save') { saveSelectionA(); b.textContent = 'Saved as A'; setTimeout(() => renderView(), 700); }
      else { go('compare'); window.scrollTo({ top: 0 }); }
    }));
  }
  function renderCompare(el) {
    const A = selectionA(), B = Object.assign(currentSelection(), { label: selectionText() });
    const date = new Date().toISOString().slice(0, 10);
    if (!A) {
      setPrintPageStyle(false);
      el.innerHTML = '<div class="card"><h3>Compare two selections</h3><p class="sub">Compare one selection (A) with another (B): their site features, interaction evidence, sequence logos and proteins, side by side, as a page you can save as PDF.</p>' +
        '<ol class="steps"><li>Search or filter to the first selection, e.g. <code>binds:PRMT1</code>, and press <b>Save as A</b>.</li><li>Change the search or filters to the second selection, e.g. <code>binds:PRMT5</code>.</li><li>Open <b>Compare</b>; the current selection is B.</li></ol>' +
        '<p><button class="btn" type="button" id="cmp-save">Save the current selection as A</button> <span class="muted small">Current: ' + esc(B.label) + '</span></p></div>';
      $('#cmp-save', el).addEventListener('click', () => { saveSelectionA(); renderCompare(el); });
      return;
    }
    setPrintPageStyle(true, 'Comparison report · ' + date);
    const both = A.sites.filter((s) => B.shown.has(s.j)).length;
    const protBoth = A.prots.filter((p) => B.counts[p.i] > 0).length;
    const rows = A.sites.length && B.sites.length ? compareRows(A, B) : [];
    const sig = rows.filter((r) => r.q <= 0.05);
    const show = (sig.length ? sig : rows).slice(0, 16);
    const pc = (x, n) => (n ? (100 * x / n).toFixed(x / n < 0.1 ? 1 : 0) : '0') + '%';
    const status = (sel) => { const n = { exp: 0, pred: 0, none: 0 }; if (ixReady()) for (const s of sel.sites) n[ixStatus(s)]++; return n; };
    const stA = status(A), stB = status(B);
    const top = (sel) => sel.prots.slice().sort((a, b) => sel.counts[b.i] - sel.counts[a.i] || a.label.localeCompare(b.label)).slice(0, 12);
    const logo = (sel, name) => { const list = sel.sites.filter((s) => s.win && s.win.length === 31); return list.length ? alignSvg(list, { flank: 7, rows: 0, logo: 'freq', title: name, sub: fmt(list.length) + ' sites with a sequence window' }).svg : '<p class="note">No sequence windows.</p>'; };
    const side = (sel, tag) => '<div class="cmp-side"><span class="cmp-tag cmp-' + tag.toLowerCase() + '">' + tag + '</span><p class="cmp-label">' + esc(sel.label) + '</p>' +
      '<p class="cmp-nums"><b>' + fmt(sel.sites.length) + '</b> sites · <b>' + fmt(sel.prots.length) + '</b> proteins</p></div>';
    const marginBoxes = /Chrome\/|Chromium\/|Edg\//.test(navigator.userAgent);
    el.innerHTML = '<div class="report cmp' + (marginBoxes ? '' : ' runfoot') + '">' +
      '<div class="print-only report-runfoot">' + esc(siteName() + ' · ' + siteUrl() + ' · ' + reportCite()) + '</div>' +
      '<div class="report-bar screen-only"><button class="linkish" type="button" id="cmp-back">← Back</button><button class="btn" type="button" id="cmp-print">Save as PDF</button>' +
      '<button class="btn quiet small" type="button" id="cmp-swap" title="Make the current selection A">Save current as A</button><button class="btn quiet small" type="button" id="cmp-clear">Clear A</button>' +
      '<span class="muted small">B is the current search and filters; change them to compare something else with A.</span></div>' +
      '<header class="report-head"><p class="report-kicker">' + esc(siteName()) + ' · comparison report · ' + date + '</p><h2>A vs B</h2>' +
      '<div class="cmp-sides">' + side(A, 'A') + side(B, 'B') + '</div>' +
      '<p class="muted small">Overlap: ' + fmt(both) + ' site' + (both === 1 ? '' : 's') + ' and ' + fmt(protBoth) + ' protein' + (protBoth === 1 ? '' : 's') + ' are in both. Reopen: <span class="mono">' + esc(siteUrl() + location.search + '#compare') + '</span></p></header>' +
      '<section class="card report-sec"><h3>In words</h3><div class="cmp-grid"><p><span class="cmp-tag cmp-a">A</span> ' + esc(summaryText(A, null)) + '</p><p><span class="cmp-tag cmp-b">B</span> ' + esc(summaryText(B, null)) + '</p></div></section>' +
      (ixReady() ? '<section class="card report-sec"><h3>Interaction evidence at the sites</h3><div class="tbl-wrap"><table class="tbl compact cmp-tbl"><thead><tr><th>Status</th><th class="r">A</th><th class="r">B</th></tr></thead><tbody>' +
        ['exp', 'pred', 'none'].map((k) => '<tr><td><span class="ixs ixs-' + k + '"><i></i>' + IXS[k][0] + '</span></td><td class="r num">' + fmt(stA[k]) + ' <span class="muted">(' + pc(stA[k], A.sites.length) + ')</span></td><td class="r num">' + fmt(stB[k]) + ' <span class="muted">(' + pc(stB[k], B.sites.length) + ')</span></td></tr>').join('') +
        '</tbody></table></div></section>' : '') +
      '<section class="card report-sec"><h3>Site features, A against B</h3>' + (rows.length ?
        '<div class="tbl-wrap"><table class="tbl compact cmp-tbl"><thead><tr><th>Feature</th><th class="r">A</th><th class="r">B</th><th class="r">A / B</th><th class="r">q</th></tr></thead><tbody>' +
        show.map((r) => '<tr' + (r.q > 0.05 ? ' class="dim"' : '') + '><td><span class="muted small">' + esc(r.ft.group) + '</span> ' + esc(r.ft.label) + '</td><td class="r num">' + pc(r.a, r.nA) + '</td><td class="r num">' + pc(r.b, r.nB) + '</td><td class="r num"><b>' + (r.fold >= 10 ? r.fold.toFixed(0) : r.fold.toFixed(1)) + '×</b></td><td class="r num">' + fmtP(r.q) + '</td></tr>').join('') +
        '</tbody></table></div><p class="note">' + (sig.length ? 'Features that differ at BH q ≤ 0.05' : 'No feature differs at BH q ≤ 0.05; the closest are shown, faded') + '. Fisher\'s exact test per feature on the two sets of sites' + (both ? '; the ' + fmt(both) + ' sites in both selections count on both sides, which makes the test conservative' : '') + '. Sites of one protein are not independent, so read this as a pointer.</p>'
        : '<p class="note">One of the selections has no sites.</p>') + '</section>' +
      '<section class="card report-sec"><h3>Sequence logos</h3><div class="cmp-grid"><div class="fig-preview">' + logo(A, 'A · ' + A.label.slice(0, 60)) + '</div><div class="fig-preview">' + logo(B, 'B · ' + B.label.slice(0, 60)) + '</div></div></section>' +
      '<section class="card report-sec"><h3>Proteins with most sites</h3><div class="cmp-grid">' + [[A, 'A'], [B, 'B']].map(([sel, t]) => '<div><div class="tbl-wrap"><table class="tbl compact"><thead><tr><th>' + t + '</th><th class="r">Sites</th><th class="r">In ' + (t === 'A' ? 'B' : 'A') + '</th></tr></thead><tbody>' +
        top(sel).map((p) => '<tr><td><button class="linkish" type="button" data-open-p="' + p.i + '">' + esc(p.label) + '</button></td><td class="r num">' + sel.counts[p.i] + '</td><td class="r num">' + ((t === 'A' ? B : A).counts[p.i] || '<span class="muted">·</span>') + '</td></tr>').join('') + '</tbody></table></div></div>').join('') + '</div></section>' +
      '<section class="card report-sec"><h3>Sources to cite</h3>' + citeBoxHtml(Array.from(new Map(A.sites.concat(B.sites).map((x) => [x.j, x])).values()), 'the sites in A and B', true) + '</section>' +
      '<section class="card report-sec"><h3>Methods</h3><p class="small">Selections are the search terms and filters shown above, applied to the integrated human methylarginine union. Interaction status: experimental = RNA or DNA crosslink, RBDmap fragment or PDB interface at the residue; predicted only = predicted-complex interface or confident predicted contact without experimental evidence. Logos show residue frequencies in the ±7 windows. ' + esc(creditLines().join(' · ')) + '</p></section>' +
      '</div>';
    $('#cmp-back', el).addEventListener('click', () => { if (history.length > 1) history.back(); else go('sites'); });
    $('#cmp-print', el).addEventListener('click', () => { const t = document.title; document.title = siteName() + ' comparison ' + date; window.print(); document.title = t; });
    $('#cmp-swap', el).addEventListener('click', () => { saveSelectionA(); renderCompare(el); });
    $('#cmp-clear', el).addEventListener('click', () => { UI.cmpA = null; syncUrl(); renderScope(); renderCompare(el); });
    $$('[data-open-p]', el).forEach((b) => b.addEventListener('click', () => openProtein(+b.dataset.openP)));
    bindCiteBox(el);
  }

  // ------------------------------------------- what sets this selection apart
  // Site features compared between the current selection and the remaining sites (Fisher's exact test, BH across the
  // features tested). Sites of one protein are not independent, so this points to differences rather than testing them.
  function profileFeatures() {
    const out = [];
    const add = (key, group, label, test, apply, term) => out.push({ key, group, label, test, apply, term });
    M.motifs.forEach((m, b) => add('mot' + b, 'Motif', MOTIF_LABEL[m] + ' · ' + MOTIF_HELP[m], (s) => s.mot & (1 << b), null, ['motif', motifTerm(m)]));
    const dIdx = (k) => M.disorder_classes.indexOf(k);
    [['disordered', 'Inside an IDR'], ['ordered', 'In an ordered region']].forEach(([k, l]) => { const i = dIdx(k); if (i >= 0) add('dis' + i, 'Structure', l, (s) => s.dc === i, (f) => { f.disorder = new Set([i]); }); });
    [[1, 'Rme1 (source call)'], [2, 'Rme2 (source call)']].forEach(([v, l]) => add('st' + v, 'Methylation', l, (s) => s.st === v || s.st === 3, (f) => { f.states.add(v); }));
    M.maron_ips.forEach((n, b) => add('ip' + b, 'Methylation', 'Recovered by the ' + n + ' IP (' + ipLabel(n) + ')', (s) => s.ip & (1 << b), (f) => { f.ips.add(b); }));
    add('gsk', 'PRMT inhibitors', 'Decreased with GSK591 (PRMT5i)', (s) => s.trf & 1, (f) => { f.treat = '0:down'; });
    add('ms', 'PRMT inhibitors', 'Decreased with MS023 (type I PRMTi)', (s) => s.trf & 4, (f) => { f.treat = '1:down'; });
    add('nn5', 'Clustering', 'Another methyl-R within 5 aa', (s) => s.nn != null && s.nn <= 5, (f) => { f.nn = 5; });
    add('rep2', 'Provenance', 'Reported by ≥2 studies', (s) => popcount(s.study) >= 2, (f) => { f.minRep = 2; });
    if (ixReady()) EV_SETS.forEach(([b, l]) => add('ev' + b, 'Interactions', l, (s) => (s.ev || 0) & b, (f) => { f.evAll |= b; f.evExact = false; }));
    return out;
  }

  function profileRows() {
    const inSel = R.shown, n1 = R.sites.length, n2 = S.length - n1;
    if (n1 < 5 || n2 < 5) return null;
    const rows = profileFeatures().map((ft) => {
      let a = 0, c = 0;
      for (const s of S) if (ft.test(s)) { if (inSel.has(s.j)) a++; else c++; }
      return { ft, a, c, n1, n2 };
    }).filter((r) => r.a + r.c >= 10 && r.a < r.n1 && !(r.a === 0 && r.c === 0));
    for (const r of rows) {
      r.p = fisherP(r.a, r.n1 - r.a, r.c, r.n2 - r.c);
      r.fold = (r.a / r.n1) / ((r.c || 0.5) / r.n2);
    }
    const order = rows.slice().sort((x, y) => x.p - y.p);
    let prev = 1;
    for (let k = order.length - 1; k >= 0; k--) { prev = Math.min(prev, (order[k].p * order.length) / (k + 1)); order[k].q = prev; }
    return order;
  }

  function profileCardHtml() {
    if (!constraints().length) return '';
    const rows = profileRows();
    if (!rows) return '';
    const sig = rows.filter((r) => r.q <= 0.05);
    const up = sig.filter((r) => r.fold > 1).slice(0, 8), down = sig.filter((r) => r.fold < 1).slice(0, 4);
    const show = up.concat(down);
    const pc = (x, n) => (100 * x / n).toFixed(x / n < 0.1 ? 1 : 0) + '%';
    return '<div class="card profile-card" id="profile"><div class="card-head"><div><h3>What sets this selection apart</h3><p class="sub">Site features in these ' + fmt(R.sites.length) + ' sites against the other ' + fmt(S.length - R.sites.length) + '. Click a feature to narrow the selection to it.' + docLink('profile') + '</p></div>' +
      '<div class="exp"><button type="button" data-prof="figures">Logo and comparison</button><button type="button" data-prof="evidence">Evidence map</button><button type="button" data-prof="report" title="Profile, evidence map, logo and protein comparison of this selection as one PDF">Insight report</button></div></div>' +
      '<div class="cmp-row">' + cmpButtons() + '</div>' +
      (show.length ? '<div class="tbl-wrap"><table class="tbl compact prof-tbl"><thead><tr><th>Feature</th><th class="r">Here</th><th class="r hide-xs">Other sites</th><th class="r">Ratio</th><th class="r">q</th></tr></thead><tbody>' +
        show.map((r) => '<tr class="row' + (r.fold < 1 ? ' depleted' : '') + '" tabindex="0" data-prof-key="' + r.ft.key + '" title="' + (r.fold > 1 ? 'Narrow the selection to sites with this feature' : 'Less common here than in the other sites') + '"><td><span class="muted small">' + esc(r.ft.group) + '</span> ' + esc(r.ft.label) + '</td>' +
          '<td class="r num">' + fmt(r.a) + ' <span class="muted">(' + pc(r.a, r.n1) + ')</span></td><td class="r num hide-xs">' + pc(r.c, r.n2) + '</td><td class="r num"><b>' + (r.fold >= 10 ? r.fold.toFixed(0) : r.fold.toFixed(1)) + '×</b></td><td class="r num">' + fmtP(r.q) + '</td></tr>').join('') + '</tbody></table></div>' :
        '<p class="note">No site feature differs significantly (BH q ≤ 0.05) between this selection and the other sites.</p>') +
      '<p class="note">Fisher\'s exact test per feature, Benjamini–Hochberg across the ' + rows.length + ' features tested. Sites in one protein are not independent, so read this as a pointer.</p></div>';
  }

  function bindProfile(el) {
    const box = $('#profile', el);
    if (!box) return;
    const feats = new Map(profileFeatures().map((f) => [f.key, f]));
    bindRows(box, (tr) => { const ft = feats.get(tr.dataset.profKey); if (!ft || tr.classList.contains('depleted')) return; if (ft.term) addTerm(ft.term[0], ft.term[1], UI.view); else setFilters((f) => ft.apply(f), UI.view); });
    bindCmpButtons(box);
    $$('[data-prof]', box).forEach((b) => b.addEventListener('click', () => {
      if (b.dataset.prof === 'figures') { go('figures'); window.scrollTo({ top: 0 }); }
      else if (b.dataset.prof === 'report') openReport();
      else { UI.ixScrollMap = true; go('interactions'); }
    }));
    bindDocLinks(box);
  }

  // ----- complex members, expanded in place from the complexes table
  function cxMembersHtml(c) {
    const ids = c.members.map((m) => (Array.isArray(m) ? m : [m]));
    const meth = [], other = [];
    for (const alts of ids) {
      const hit = alts.find((k) => k < IX.nP);
      if (hit != null) meth.push(hit); else other.push(alts[0]);
    }
    meth.sort((a, b) => P[b].sites.length - P[a].sites.length);
    return '<div class="cx-members"><div class="gene-chips">' + meth.map((i) => '<button class="gene-chip" type="button" data-i="' + i + '" title="Open ' + esc(P[i].label) + '">' + esc(P[i].label) + '<span>' + P[i].sites.length + ' site' + (P[i].sites.length === 1 ? '' : 's') + '</span></button>').join('') +
      other.map((k) => '<span class="gene-chip off" title="No methylarginine in the set">' + esc(ixGene(k)) + '</span>').join('') + '</div>' +
      '<div class="ev-actions"><button class="btn quiet small" type="button" data-cxact="list" data-k="' + c.k + '">List methylated members</button><button class="btn quiet small" type="button" data-cxact="compare" data-k="' + c.k + '">Compare members (figure)</button><button class="btn quiet small" type="button" data-cxact="evidence" data-k="' + c.k + '">Evidence map of their sites</button></div></div>';
  }

  function bindComplexMembers(el) {
    const tbl = $('#ch-cx', el);
    if (!tbl) return;
    $$('[data-cxmem]', tbl).forEach((b) => b.addEventListener('click', (e) => {
      e.stopPropagation();
      const k = +b.dataset.cxmem, tr = b.closest('tr');
      const open = tr.nextElementSibling && tr.nextElementSibling.classList.contains('cx-mem-row');
      if (open) { tr.nextElementSibling.remove(); UI.cxOpen.delete(k); b.setAttribute('aria-expanded', 'false'); return; }
      UI.cxOpen.add(k);
      b.setAttribute('aria-expanded', 'true');
      tr.insertAdjacentHTML('afterend', '<tr class="cx-mem-row"><td colspan="5">' + cxMembersHtml(IX.cx[k]) + '</td></tr>');
      bindMemberRow(tr.nextElementSibling);
    }));
    $$('tr.cx-mem-row', tbl).forEach(bindMemberRow);
  }
  function bindMemberRow(row) {
    $$('.gene-chip[data-i]', row).forEach((b) => b.addEventListener('click', () => openProtein(+b.dataset.i)));
    $$('[data-cxact]', row).forEach((b) => b.addEventListener('click', () => {
      const k = b.dataset.k, act = b.dataset.cxact;
      if (act === 'compare') {
        const c = IX.cx[+k];
        const genes = c.members.map((m) => (Array.isArray(m) ? m : [m]).find((x) => x < IX.nP)).filter((x) => x != null).map((i) => P[i].label);
        applyQueryText(queryText(Array.from(new Set(genes)).join(' '), Q.tokens.filter((t) => t.key !== 'complex').concat([{ key: 'complex', vals: [{ text: c.id }] }])), 'figures');
      } else if (act === 'evidence') { UI.ixScrollMap = true; setTerm('complex', IX.cx[+k].id, 'interactions'); }
      else setTerm('complex', IX.cx[+k].id, 'proteins');
      window.scrollTo({ top: 0 });
    }));
  }

  function renderInteractions(el) {
    if (!IX) { el.innerHTML = '<div class="card"><h3>Interactions and binding</h3><p class="note">Loading interaction data…</p></div>'; return; }
    if (IX.failed) { el.innerHTML = '<div class="card"><h3>Interactions and binding</h3><p class="note">Interaction data (<code>data/methylarg_interactions.json</code>) is not available.</p></div>'; return; }
    const body = { cells: ixSecCells, pdb: ixSecPdb, predicted: ixSecPredicted, partners: ixSecPartners, complexes: ixSecComplexes, go: ixSecGo };
    el.innerHTML = '<div class="ix-page"><div class="card ix-intro"><h3>Interactions and binding</h3>' +
      '<p class="sub">What methylated arginines and their proteins touch, from five kinds of evidence: physical partners (STRING), complexes (Complex Portal), experimental structures (PDB), predicted complexes (ProtVar, AlphaFold DB) and RNA or DNA crosslinking in living human cells. ' +
      'Most tests compare methylated arginines with the other arginines of the same proteins (Mantel–Haenszel, each protein a stratum), so a protein with many mapped sites does not dominate. Click any row in a section to filter the sites.</p>' +
      checkSiteHtml('ix-check') +
      '<h4>Key findings</h4>' + findingsList(ixFindings()) +
      '<div class="ix-toc"><span class="flabel">Sections</span>' + IX_SECTIONS.map(([k, title]) => '<button class="chip" type="button" data-ixgo="' + k + '">' + esc(title) + '</button>').join('') +
      '<button class="linkish" type="button" id="ix-all">' + (UI.ixOpen.size === IX_SECTIONS.length ? 'Collapse all' : 'Expand all') + '</button></div></div>' +
      ixExampleHtml() +
      ixEvidenceMapHtml() +
      IX_SECTIONS.map(([k, title, tag]) => '<details class="card ix-sec" id="ixs-' + k + '" data-ixsec="' + k + '"' + (UI.ixOpen.has(k) ? ' open' : '') + '><summary><span class="pm" aria-hidden="true"></span><span class="ix-sec-t"><b>' + esc(title) + '</b><span class="sub">' + esc(tag) + '</span></span></summary>' +
        '<div class="ix-sec-body">' + body[k]() + '</div></details>').join('') + '</div>';
    $$('details.ix-sec', el).forEach((d) => d.addEventListener('toggle', () => {
      if (d.open) UI.ixOpen.add(d.dataset.ixsec); else UI.ixOpen.delete(d.dataset.ixsec);
      const all = $('#ix-all', el);
      if (all) all.textContent = UI.ixOpen.size === IX_SECTIONS.length ? 'Collapse all' : 'Expand all';
    }));
    $('#ix-all', el).addEventListener('click', () => {
      const open = UI.ixOpen.size !== IX_SECTIONS.length;
      $$('details.ix-sec', el).forEach((d) => { d.open = open; });
    });
    bindIxGo(el);
    bindIxOverview(el);
    bindCheckSite(el);
    bindIxExample(el);
    bindChartExports(el);
    bindEvidenceMap(el);
    bindComplexMembers(el);
    bindDocLinks(el);
    if (UI.ixJump) { const k = UI.ixJump; UI.ixJump = null; openIxSection(k); }
    if (UI.ixScrollMap) { UI.ixScrollMap = false; const m = $('#ix-evmap', el); if (m) setTimeout(() => m.scrollIntoView({ block: 'start' }), 0); }
  }

  function openIxSection(k) {
    const d = $('#ixs-' + k);
    if (!d) return;
    d.open = true;
    UI.ixOpen.add(k);
    d.scrollIntoView({ block: 'start', behavior: 'smooth' });
  }

  // "Details" links: open a section here, or go to the Interactions tab first.
  function bindIxGo(el) {
    $$('[data-ixgo]', el).forEach((b) => b.addEventListener('click', () => {
      const k = b.dataset.ixgo;
      if (UI.view === 'interactions') { if (k) openIxSection(k); return; }
      UI.ixJump = k || null;
      go('interactions');
      window.scrollTo({ top: 0 });
    }));
  }

  function ixSecPartners() {
    const bound = boundBy(700);
    const readerRows = (list) => list.map((x) => ({ x, n: bound.get(x.k) || 0 })).filter((r) => r.n).sort((a, b) => b.n - a.n || a.x.gene.localeCompare(b.x.gene));
    const prmtRows = readerRows(IX.prmts), tudorAll = readerRows(IX.tudors), wdAll = readerRows(IX.wd40s), phdAll = readerRows(IX.phds);
    const groups = [['Protein arginine methyltransferases', prmtRows, prmtRows.length], ['Tudor-domain proteins', tudorAll.slice(0, 12), tudorAll.length],
      ['WD40-repeat proteins', wdAll.slice(0, 12), wdAll.length], ['PHD-finger proteins', phdAll.slice(0, 12), phdAll.length]];
    const maxN = Math.max(1, ...groups.flatMap((gr) => gr[1].map((r) => r.n)));
    const barList = (list) => '<div class="hbars">' + list.map(({ x, n }) => '<button class="hb" type="button" data-partner="' + x.k + '"><span class="hb-label">' + esc(x.gene) + (x.p >= 0 ? ' <span class="muted">· methylated</span>' : '') + '</span><span class="hb-track"><span class="hb-bar" style="left:0;width:' + (100 * n / maxN).toFixed(1) + '%"></span></span><span class="hb-val">' + fmt(n) + '</span></button>').join('') + '</div>';
    let edgesMeth = 0, edgesAll = 0;
    for (const p of P) for (const [k, sc] of p.partners || []) if (sc >= 700) { edgesAll++; if (k < IX.nP) edgesMeth++; }
    return '<div id="ch-readers"><p class="sub">Methylated proteins with a high-confidence (≥ 0.7) STRING physical link to each protein. Tudor domains (InterPro IPR002999) include the methylarginine readers SMN1, SMNDC1, TDRD3 and SND1; WD40 repeats (IPR001680) include WDR77 (MEP50), the PRMT5 cofactor, and WDR5; PHD fingers (IPR001965, IPR019787), many of which read the histone H3 N-terminus. Click a bar to list the methylated proteins it binds; the filter panel also takes any other partner by gene symbol.</p>' +
      '<div class="two-col">' + groups.map(([label, list, n]) => '<div><p class="flabel">' + esc(label) + (n > list.length ? ' <span class="muted">(top ' + list.length + ' of ' + n + ')</span>' : '') + '</p>' + (list.length ? barList(list) : '<p class="note">None at this confidence.</p>') + '</div>').join('') + '</div>' +
      '<p class="note">Summed over methylated proteins, ' + fmt(edgesMeth) + ' of their ' + fmt(edgesAll) + ' high-confidence physical partners (' + pctText(edgesMeth, edgesAll || 1) + '%) are themselves methylated, against ' + pctText(IX.meta.methylated_in_proteome, IX.meta.proteome_size) + '% of reviewed human proteins. Well-studied proteins carry both more mapped sites and more mapped interactions, so read this as association, not mechanism.</p></div>';
  }

  function ixSecComplexes() {
    const min = UI.cxMin, curated = !UI.cxPred;
    const rows = IX.cx.filter((c) => c.pred === !curated && c.n >= min && c.m > 0)
      .sort(UI.cxSort === 'share' ? (a, b) => b.m / b.n - a.m / a.n || b.n - a.n : (a, b) => (a.q - b.q) || (b.m / b.n - a.m / a.n));
    const shown = rows.slice(0, UI.cxAll ? rows.length : 15);
    return '<div id="ch-cx"><div class="card-head"><div><p class="sub">Complex Portal complexes ranked by their arginine-methylated members (hypergeometric test against the ' + fmt(IX.meta.proteome_size) + ' reviewed human proteins, ' + pctText(IX.meta.methylated_in_proteome, IX.meta.proteome_size) + '% of which are in this set; BH-adjusted). Click a complex to list its methylated members, or Members to see them all here.</p></div>' +
      '<div class="exp"><button type="button" id="cx-tsv">TSV</button></div></div>' +
      '<div class="fig-controls"><div class="seg" role="group" aria-label="Complex set"><span>Complexes</span><button type="button" data-cxset="cur" aria-pressed="' + curated + '">Curated</button><button type="button" data-cxset="pred" aria-pressed="' + !curated + '">Predicted (hu.MAP 3.0)</button></div>' +
      '<label class="fig-ctl">Min. members <select id="cx-min">' + [2, 3, 5, 10, 20].map((n) => '<option' + (n === min ? ' selected' : '') + '>' + n + '</option>').join('') + '</select></label>' +
      '<label class="fig-ctl">Rank by <select id="cx-sort"><option value="q"' + (UI.cxSort !== 'share' ? ' selected' : '') + '>Enrichment q</option><option value="share"' + (UI.cxSort === 'share' ? ' selected' : '') + '>Share methylated</option></select></label></div>' +
      '<div class="tbl-wrap"><table class="tbl compact cx-tbl"><thead><tr><th>Complex</th><th class="r">Members</th><th>Methylated</th><th class="r hide-xs">Fold</th><th class="r">q</th></tr></thead><tbody>' +
      shown.map((c) => '<tr class="row" tabindex="0" data-cx="' + c.k + '"><td><button class="linkish cx-mem-btn" type="button" data-cxmem="' + c.k + '" aria-expanded="' + UI.cxOpen.has(c.k) + '" title="Show the members">Members</button> <b>' + esc(c.name) + '</b> <span class="muted mono">' + esc(c.id) + '</span></td><td class="r num">' + c.n + '</td>' +
        '<td><div class="meter-row"><span class="meter"><i style="width:' + (100 * c.m / c.n).toFixed(1) + '%"></i></span><span class="num">' + c.m + ' (' + Math.round(100 * c.m / c.n) + '%)</span></div></td>' +
        '<td class="r num hide-xs">' + (c.fold != null ? c.fold.toFixed(1) + '×' : '—') + '</td><td class="r num">' + fmtP(c.q) + '</td></tr>' +
        (UI.cxOpen.has(c.k) ? '<tr class="cx-mem-row"><td colspan="5">' + cxMembersHtml(c) + '</td></tr>' : '')).join('') +
      '</tbody></table></div>' + (rows.length > 15 ? '<button class="linkish" type="button" id="cx-more" style="margin-top:.5rem">' + (UI.cxAll ? 'Show the top 15' : 'Show all ' + fmt(rows.length) + ' complexes') + '</button>' : '') + '</div>';
  }

  function ixSecPdb() {
    const e = IX.enrich;
    const ifRow = (key, label, chip) => {
      const r = e[key];
      if (!r) return '';
      const [a, b, c, d, l2, pv] = r;
      return '<tr class="row" tabindex="0" data-iface="' + chip + '"><td><b>' + esc(label) + '</b></td><td class="r num">' + fmt(a) + ' / ' + fmt(a + b) + ' (' + pctText(a, a + b || 1) + '%)</td><td class="r num">' + fmt(c) + ' / ' + fmt(c + d) + ' (' + pctText(c, c + d || 1) + '%)</td><td class="r num">' + fmtL2(l2) + '</td><td class="r num">' + fmtP(pv) + '</td></tr>';
    };
    return '<p class="sub">Do methylated arginines sit at protein, RNA or DNA interfaces in PDB structures (PDBe-KB) more often than the other arginines of the same proteins? Only structurally observed arginines (SIFTS) count, and most methylarginines sit in disordered regions that are rarely crystallised.</p>' +
      strataTable() +
      '<details class="dataset-ctx"><summary class="flabel">Pooled across proteins (not controlled; click a row to filter sites)</summary><div class="tbl-wrap"><table class="tbl compact"><thead><tr><th>Arginines at an interface with</th><th class="r">Methylated R</th><th class="r">Other R</th><th class="r">log2 OR</th><th class="r">p</th></tr></thead><tbody>' +
      ifRow('rna', 'RNA', 4) + ifRow('dna', 'DNA', 8) + ifRow('protein', 'Protein', 1) + '</tbody></table></div>' +
      '<p class="note">Pooling mixes proteins: proteins solved bound to RNA or DNA also tend to carry more mapped methylarginines, which inflates these odds ratios. The within-protein test above removes that.</p></details>' +
      geometryHtml();
  }

  function ixSecPredicted() {
    return predEnrichHtml().replace('<h4 style="margin-top:1rem">Predicted interfaces (AlphaFold-Multimer)</h4>', '<h4>Predicted interfaces (AlphaFold-Multimer)</h4>') +
      '<h4 style="margin-top:1.2rem">Best-supported predicted contacts</h4>' + galleryHtml();
  }

  function ixSecCells() {
    return xlEnrichHtml().replace('<h4 style="margin-top:1rem">RNA and DNA contact in living cells (crosslinking mass spectrometry)</h4>', '');
  }

  function ixSecGo() {
    const g = IX.goEnrich;
    const goRow = (key, label, chip, exp) => {
      const r = g[key];
      if (!r) return '';
      const [a, n1, c, n2, l2, pv] = r;
      return '<tr class="row" tabindex="0" data-go="' + chip + '" data-goexp="' + (exp ? 1 : 0) + '"><td><b>' + esc(label) + '</b></td><td class="r num">' + fmt(a) + ' / ' + fmt(n1) + ' (' + pctText(a, n1 || 1) + '%)</td><td class="r num">' + fmt(c) + ' / ' + fmt(n2) + ' (' + pctText(c, n2 || 1) + '%)</td><td class="r num">' + fmtL2(l2) + '</td><td class="r num">' + fmtP(pv) + '</td></tr>';
    };
    return '<div class="tbl-wrap"><table class="tbl compact"><thead><tr><th>Gene Ontology</th><th class="r">Methylated proteins</th><th class="r">Other proteins</th><th class="r">log2 OR</th><th class="r">p</th></tr></thead><tbody>' +
      goRow('rna', 'RNA binding', 1, false) + goRow('rna_exp', 'RNA binding, experimental', 1, true) + goRow('dna', 'DNA binding', 2, false) + goRow('dna_exp', 'DNA binding, experimental', 2, true) + '</tbody></table></div>' +
      '<p class="note">Methylated vs other reviewed human proteins (Fisher exact test), GO:0003723 and GO:0003677 including child terms. Click a row to list the proteins.</p>';
  }

  // The interface test repeated within proteins (Mantel–Haenszel, each protein a stratum) and restricted to
  // ordered residues (outside UniProt/MobiDB-lite disordered regions) or to residues inside annotated domains.
  function strataTable() {
    const st = IX.strata || {};
    if (!st.rna) return '';
    const cell = (key, f) => {
      const x = st[key] && st[key][f];
      if (!x) return '<td class="r">—</td>';
      const [a, b, c, d] = x.pooled;
      const [l2, pv, nStr] = x.mh;
      return '<td class="r num" title="' + esc('Methylated R at the interface ' + a + '/' + (a + b) + ' · other R ' + c + '/' + (c + d) + ' · pooled log2 OR ' + fmtL2(x.pooled[4]) + ' · Mantel–Haenszel over ' + nStr + ' informative proteins') + '">' +
        (l2 == null ? '—' : '<b>' + fmtL2(l2) + '</b> <span class="muted">p ' + fmtP(pv) + '</span>') + '</td>';
    };
    return '<div class="tbl-wrap"><table class="tbl compact strata-tbl"><thead><tr><th>Within-protein test (log2 MH OR)</th><th class="r">All observed R</th><th class="r">Ordered regions only</th><th class="r">Inside domains only</th></tr></thead><tbody>' +
      [['rna', 'RNA'], ['dna', 'DNA'], ['protein', 'Protein']].map(([k, l]) => '<tr><td><b>' + l + '</b></td>' + cell(k, 'all') + cell(k, 'ordered') + cell(k, 'domain') + '</tr>').join('') +
      '</tbody></table></div><p class="note">Each protein is its own stratum (Cochran–Mantel–Haenszel), so methylated arginines are only compared with other arginines of the same protein. "Ordered" excludes UniProt disordered regions (MobiDB-lite); "inside domains" keeps residues in UniProt domain, repeat or zinc-finger features. Values near 0 mean no difference; hover a cell for the counts.</p>';
  }

  // Where in an interface methylated arginines sit (Levy core/rim/support, nucleic-acid hydrogen bonds) and the
  // RG/RGG context comparison, from the export's residue-environment analysis.
  function geometryHtml() {
    const gs = IX.geoSummary || {};
    const rg = IX.rg || {};
    const mh = (x) => (x && x.mh && x.mh[0] != null ? '<b>' + fmtL2(x.mh[0]) + '</b> <span class="muted">p ' + fmtP(x.mh[1]) + '</span>' : '—');
    const pooledTip = (x) => (x && x.pooled ? ' title="' + esc('Methylated ' + x.pooled[0] + '/' + (x.pooled[0] + x.pooled[1]) + ' · other ' + x.pooled[2] + '/' + (x.pooled[2] + x.pooled[3]) + ' · Mantel–Haenszel over ' + x.mh[2] + ' informative proteins') + '"' : '');
    let html = '';
    if (gs.levy_counts) {
      const m = gs.levy_counts['1'] || {}, c = gs.levy_counts['0'] || {};
      const buried = (x) => (x.core || 0) + (x.rim || 0) + (x.support || 0);
      const row = (k, label) => '<tr><td>' + label + '</td><td class="r num">' + fmt(m[k] || 0) + ' <span class="muted">(' + pctText(m[k] || 0, buried(m) || 1) + '%)</span></td><td class="r num">' + fmt(c[k] || 0) + ' <span class="muted">(' + pctText(c[k] || 0, buried(c) || 1) + '%)</span></td></tr>';
      const nm = gs.na_counts['1'] || {}, nc = gs.na_counts['0'] || {};
      const naTot = (x) => (x.base || 0) + (x.backbone || 0) + (x.contact || 0);
      const nrow = (k, label) => '<tr><td>' + label + '</td><td class="r num">' + fmt(nm[k] || 0) + ' <span class="muted">(' + pctText(nm[k] || 0, naTot(nm) || 1) + '%)</span></td><td class="r num">' + fmt(nc[k] || 0) + ' <span class="muted">(' + pctText(nc[k] || 0, naTot(nc) || 1) + '%)</span></td></tr>';
      html += '<h4 style="margin-top:1rem">Where in the interface?</h4><div class="two-col geo-cols"><div class="tbl-wrap"><table class="tbl compact"><thead><tr><th>Position (Levy classes)</th><th class="r">Methylated R</th><th class="r">Other interface R</th></tr></thead><tbody>' +
        row('core', 'Interface core') + row('rim', 'Interface rim') + row('support', 'Interface support') +
        '<tr class="muted"><td>Not buried in the chosen entry</td><td class="r num">' + fmt(m.none || 0) + '</td><td class="r num">' + fmt(c.none || 0) + '</td></tr></tbody></table>' +
        '<p class="note"' + pooledTip(gs.core) + '>Core vs rim or support, within protein: ' + mh(gs.core) + '</p></div>' +
        '<div class="tbl-wrap"><table class="tbl compact"><thead><tr><th>Contact with RNA or DNA</th><th class="r">Methylated R</th><th class="r">Other interface R</th></tr></thead><tbody>' +
        nrow('base', 'Hydrogen bond to a base') + nrow('backbone', 'Hydrogen bond to the backbone only') + nrow('contact', 'Contact without hydrogen bond') + '</tbody></table>' +
        '<p class="note"' + pooledTip(gs.base) + '>Base H-bond vs other nucleic-acid contacts, within protein: ' + mh(gs.base) + '</p></div></div>' +
        '<p class="note">One PDB entry per residue (the most recent showing the interface); controls are up to ' + (gs.controls_per_protein || 12) + ' non-methylated interface arginines per protein. Relative accessibility by Shrake–Rupley over a 10 Å environment (RCSB ModelServer), Arg maximum 274 Å² (Tien et al. 2013); Levy (2010) classes: core = exposed alone and buried (&lt; 25%) in the complex. Hydrogen bonds: guanidino N to partner N/O ≤ 3.5 Å.</p>';
    }
    if (rg.segment) {
      const ctxRow = (k, label) => {
        const x = rg[k] || {};
        const cov = x.cover || [0, 0];
        return '<tr><td>' + label + '</td><td class="r num">' + fmt(cov[0]) + '</td><td class="r num">' + fmt(cov[1]) + ' <span class="muted">(' + pctText(cov[1], cov[0] || 1) + '%)</span></td><td class="r num"' + pooledTip(x.rna) + '>' + mh(x.rna) + '</td><td class="r num"' + pooledTip(x.any) + '>' + mh(x.any) + '</td></tr>';
      };
      html += '<h4 style="margin-top:1rem">RG/RGG context</h4><div class="tbl-wrap"><table class="tbl compact"><thead><tr><th>Arginine context</th><th class="r">Methylated sites</th><th class="r">Seen in a PDB structure</th><th class="r">RNA interface (log2 MH OR)</th><th class="r">Any interface</th></tr></thead><tbody>' +
        ctxRow('segment', 'In an RG/RGG-rich segment') + ctxRow('dipeptide', 'RG dipeptide elsewhere') + ctxRow('other', 'Other arginines') + '</tbody></table></div>' +
        '<p class="note">Within each context, methylated arginines are compared with the other arginines of the same context and protein. RG/RGG-rich segments: ≥ 3 RG dipeptides with ≤ 12 aa between consecutive RG arginines.' +
        (rg.segment.rna && rg.segment.rna.mh && rg.segment.rna.mh[0] == null ? ' Only ' + fmt((rg.segment.cover || [0, 0])[1]) + ' methylated arginines in RG/RGG-rich segments are resolved in any PDB structure, too few for a within-protein test: these low-complexity regions are rarely crystallised or modelled.' : '') + '</p>';
    }
    if (IX.nModelled) html += '<p class="note"><button class="linkish" type="button" id="geo-modelled">' + fmt(IX.nModelled) + ' union sites are modelled as a methylated arginine in at least one PDB entry</button>; PDB entries model a further ' + fmt(Array.from(IX.modelledGaps.values()).reduce((n, g) => n + g.pos.length, 0)) + ' methylated positions that the union lacks (listed on each protein page).</p>';
    return html;
  }

  // Crosslinking MS in living cells: methylated vs other arginines of the same proteins.
  function xlEnrichHtml() {
    const e = IX.xlEnrich || {};
    if (!e.rna_xl) return '';
    const cellOr = (o) => (o && o.mh && o.mh[0] != null ? fmtOr(o.mh[0], o.mh[1], o.mh[2]) + ' <span class="muted">p ' + fmtP(o.mh[3]) + '</span>' : '—');
    const row = (key, label, filt) => {
      const x = e[key];
      if (!x || !x.or) return '';
      const [a, b, c, d] = x.or.pooled;
      return '<tr' + (filt ? ' class="row" tabindex="0" data-xl="' + filt + '"' : '') + '><td><b>' + label + '</b></td><td class="r num">' + cellOr(x.or) + '</td><td class="r num">' + fmt(a) + ' / ' + fmt(a + b) + ' <span class="muted">(' + pctText(a, a + b || 1) + '%)</span></td>' +
        '<td class="r num">' + fmt(c) + ' / ' + fmt(c + d) + ' <span class="muted">(' + pctText(c, c + d || 1) + '%)</span></td><td class="r num hide-xs">' + fmt(x.or.mh[4]) + ' / ' + fmt(x.proteins) + '</td></tr>';
    };
    const T = IX.xlTypes || {};
    const typeRow = (key, label) => {
      const cell = (out) => {
        const o = T[out] && T[out][key];
        if (!o) return '<td class="r">—</td><td class="r">—</td>';
        const [a, b] = o.pooled;
        return '<td class="r num">' + fmt(a) + ' / ' + fmt(a + b) + ' <span class="muted">(' + pctText(a, a + b || 1) + '%)</span></td><td class="r num">' + cellOr(o) + '</td>';
      };
      const click = key === 'rme2s' || key === 'any' ? ' class="row" tabindex="0" data-xlcls="' + key + '" title="List these RNA-crosslinked sites"' : '';
      return '<tr' + click + '><td><b>' + label + '</b></td>' + cell('rna_xl') + cell('rna_xl_single') + cell('rna_rbdmap') + '</tr>';
    };
    const ref = T.rna_xl && T.rna_xl.any ? T.rna_xl.any.pooled : null, refF = T.rna_rbdmap && T.rna_rbdmap.any ? T.rna_rbdmap.any.pooled : null;
    const refS = T.rna_xl_single && T.rna_xl_single.any ? T.rna_xl_single.any.pooled : null;
    const refCell = (r) => (r ? '<td class="r num">' + fmt(r[2]) + ' / ' + fmt(r[2] + r[3]) + ' (' + pctText(r[2], r[2] + r[3] || 1) + '%)</td><td class="r">1</td>' : '<td class="r">—</td><td class="r">1</td>');
    const sd = (k) => (T[k] && T[k].rme2s ? T[k].rme2s : null);
    const orTxt = (o, noun) => (o && o.mh && o.mh[0] != null ? 'OR ' + fmtOr(o.mh[0], o.mh[1], o.mh[2]) + ', p ' + fmtP(o.mh[3]) + ', ' + fmt(o.pooled[0]) + ' ' + (noun || 'crosslinked') : '—');
    const F_ = IX.xlFreq;
    return '<h4 style="margin-top:1rem">RNA and DNA contact in living cells (crosslinking mass spectrometry)</h4>' +
      '<p class="sub">UV crosslinking locates the residue that touches RNA or DNA in human cells, including in disordered regions that structures miss. Each row compares methylated arginines with the other arginines of the proteins that have that evidence: odds ratio within proteins (Mantel–Haenszel, each protein a stratum) with its 95% confidence interval (Robins–Breslow–Greenland) and CMH p value.</p>' +
      '<div class="tbl-wrap"><table class="tbl compact strata-tbl"><thead><tr><th>Arginine is…</th><th class="r">OR within protein [95% CI]</th><th class="r">Methylated R</th><th class="r">Other R</th><th class="r hide-xs">Informative / proteins</th></tr></thead><tbody>' +
      row('rna_xl', 'An RNA-crosslinked residue', 1) + row('rna_xl_single', '… named alone by a crosslink (single localisation)') + row('rna_xl_nonrg', '… excluding RG/GR dipeptides') + row('rna_near', 'Within 3 aa of an RNA crosslink', 2) + row('rna_rbdmap', 'Inside an RBDmap RNA-bound fragment', 4) + row('rna_rbdmap_nonrg', '… excluding RG/GR dipeptides') +
      row('dna_xl', 'A DNA crosslink localised to this residue', 8) + row('dna_cand', '… or one of its candidate residues', 32) + row('dna_pep', 'Inside a DNA-crosslinked peptide') + row('dna_bind', 'Inside a UniProt DNA-binding region', 16) + '</tbody></table></div>' +
      '<p class="note">About half of the crosslinks name several candidate residues; "single localisation" keeps only residues that at least one crosslink names alone. Few DNA crosslinks are localised to one residue, so the DNA test rests on candidates.</p>' +
      '<h4 style="margin-top:1rem">By methylation type</h4>' +
      '<p class="sub">Each class of methylated arginine against the unmethylated arginines (not in the union) of the same proteins. IP classes come from methyl-specific immunoprecipitation (Maron et al. 2021: SDMA, ADMA, MMA antibodies); state calls come from the source studies.</p>' +
      '<div class="tbl-wrap"><table class="tbl compact strata-tbl"><thead><tr><th>Methylarginines</th><th class="r">At the crosslinked residue</th><th class="r">OR [95% CI]</th><th class="r">Single localisation</th><th class="r">OR [95% CI]</th><th class="r">Inside an RBDmap fragment</th><th class="r">OR [95% CI]</th></tr></thead><tbody>' +
      typeRow('rme2s', 'Recovered by the SDMA IP (Rme2s)') + typeRow('rme2a', 'ADMA IP, not SDMA (Rme2a)') + typeRow('rme1ip', 'MMA IP only (Rme1)') + typeRow('rme2', 'Rme2 (source call)') + typeRow('rme1', 'Rme1 only (source call)') + typeRow('nostate', 'No state call or IP') + typeRow('any', 'All methylarginines') +
      (ref ? '<tr class="muted"><td>Unmethylated arginines (reference)</td>' + refCell(ref) + refCell(refS) + refCell(refF) + '</tr>' : '') +
      '</tbody></table></div>' +
      '<p class="note">Sites recovered by the SDMA IP sit inside RNA-bound RBDmap fragments more often than other arginines of the same proteins (' + orTxt(sd('rna_rbdmap'), 'inside fragments') + '), but are less often the exact residue caught crosslinked to RNA (' + orTxt(sd('rna_xl')) + '; single localisation only: ' + orTxt(sd('rna_xl_single')) + '). ' +
        'Treat the second result as suggestive: it rests on few crosslinked residues, and the pepR-MS searches identify only unmethylated crosslinked peptides (none carries a methylation mass shift), so a site that is almost fully methylated can look depleted even if its methylated copies touch RNA. Methylation weakening that residue\'s own contact would also give this pattern. The SDMA IP class means the peptide was recovered by the SDMA antibody, not that the site was shown to be symmetric dimethylated.</p>' +
      (F_ ? '<p class="note"><b>Contact frequency.</b> Among RNA-crosslinked arginines, methylated ones were seen in ' + F_.spectra.meth[1] + ' crosslinked spectra on average (median ' + F_.spectra.meth[2] + ', n = ' + fmt(F_.spectra.meth[0]) + ') against ' + F_.spectra.other[1] + ' (median ' + F_.spectra.other[2] + ', n = ' + fmt(F_.spectra.other[0]) + ') for unmethylated ones (Mann–Whitney p ' + fmtP(F_.spectra.p) + '), and in ' + F_.datasets.meth[1] + ' vs ' + F_.datasets.other[1] + ' datasets (p ' + fmtP(F_.datasets.p) + '); crosslinked in two or more datasets, within protein: OR ' + fmtOr(F_.multi.mh[0], F_.multi.mh[1], F_.multi.mh[2]) + ', p ' + fmtP(F_.multi.mh[3]) + '.</p>' : '') +
      '<p class="note">RNA: pepR-MS compendium (Sha, Kuster & Trendel 2025; with Kramer 2014, XRNAX 2019, RBS-ID 2020 and pRBS-ID 2021 re-analysed) and RBDmap (Castello et al. 2016). DNA: Trendel et al. 2025 (MCF7). Click a row of the first table to list the sites.</p>';
  }

  // CD-CODE scores cluster at 1-2 and 4-5, so a threshold of 3 or 4 keeps nearly the same members.
  function condHistNote() {
    const h = (IX && IX.condHist) || {};
    const n = (k) => h[String(k)] || 0;
    const tot = Object.values(h).reduce((a, b) => a + b, 0);
    if (!tot) return '';
    return 'Of ' + fmt(tot) + ' human memberships, ' + fmt(n(0) + n(1) + n(2)) + ' score 0–2, ' + fmt(n(3)) + ' score 3 and ' + fmt(n(4) + n(5)) + ' score 4–5, so from 3 upward the chart rests on the small high-confidence set. ';
  }

  // ----- small statistics helpers (same methods as scripts/37)
  const lfactCache = [0];
  const lfact = (n) => { for (let k = lfactCache.length; k <= n; k++) lfactCache[k] = lfactCache[k - 1] + Math.log(k); return lfactCache[n]; };
  function fisherP(a, b, c, d) {
    const r1 = a + b, r2 = c + d, c1 = a + c, n = r1 + r2;
    const lp = (x) => lfact(r1) + lfact(r2) + lfact(c1) + lfact(n - c1) - lfact(n) - lfact(x) - lfact(r1 - x) - lfact(c1 - x) - lfact(r2 - c1 + x);
    const obs = lp(a);
    let pv = 0;
    for (let x = Math.max(0, c1 - r2); x <= Math.min(r1, c1); x++) { const v = lp(x); if (v <= obs + 1e-7) pv += Math.exp(v); }
    return Math.min(1, pv);
  }
  function erfc(x) { // Numerical Recipes erfcc, fractional error < 1.2e-7
    const z = Math.abs(x), t = 1 / (1 + 0.5 * z);
    const r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
    return x >= 0 ? r : 2 - r;
  }
  function oddsCi(a, b, c, d) {
    const or = b * c ? (a * d) / (b * c) : null;
    if (min4(a, b, c, d) === 0) { a += 0.5; b += 0.5; c += 0.5; d += 0.5; }
    const l = Math.log((a * d) / (b * c)), se = Math.sqrt(1 / a + 1 / b + 1 / c + 1 / d);
    return { or, lo: Math.exp(l - 1.96 * se), hi: Math.exp(l + 1.96 * se) };
  }
  const min4 = (a, b, c, d) => Math.min(a, b, c, d);
  function mhTest(tables) {
    let num = 0, den = 0, sa = 0, sea = 0, sva = 0, used = 0, R = 0, Sx = 0, PR = 0, PSQR = 0, QS = 0;
    for (const [a, b, c, d] of tables) {
      const n = a + b + c + d;
      if (n < 2 || !(a + b) || !(c + d) || !(a + c) || !(b + d)) continue;
      used++;
      num += (a * d) / n; den += (b * c) / n;
      sa += a; sea += ((a + b) * (a + c)) / n; sva += ((a + b) * (c + d) * (a + c) * (b + d)) / (n * n * (n - 1));
      const P_ = (a + d) / n, Q_ = (b + c) / n, r = (a * d) / n, s = (b * c) / n;
      R += r; Sx += s; PR += P_ * r; PSQR += P_ * s + Q_ * r; QS += Q_ * s;
    }
    if (!used || sva <= 0) return { or: null, lo: null, hi: null, p: null, used };
    const or = den > 0 ? num / den : Infinity;
    const stat = (Math.abs(sa - sea) - 0.5) ** 2 / sva;
    const pv = erfc(Math.sqrt(Math.max(stat, 0) / 2));
    let lo = null, hi = null;
    if (R > 0 && Sx > 0 && Number.isFinite(or) && or > 0) {
      const se = Math.sqrt(PR / (2 * R * R) + PSQR / (2 * R * Sx) + QS / (2 * Sx * Sx));
      lo = Math.exp(Math.log(or) - 1.96 * se); hi = Math.exp(Math.log(or) + 1.96 * se);
    }
    return { or, lo, hi, p: pv, used };
  }
  const fmtOr = (or, lo, hi) => (or == null || !Number.isFinite(or) ? '—' : '<b>' + or.toFixed(2) + '</b>' + (lo != null || hi != null ? ' <span class="muted">[' + (lo != null ? lo.toFixed(2) : '0') + '–' + (hi != null && Number.isFinite(hi) ? hi.toFixed(2) : '∞') + ']</span>' : ''));

  // ----- predicted-contact gallery
  const GAL_DEFAULT = { conf: 'confident', pl: 80, kind: 'bond', noPdb: false, xl: false };
  const GAL_BOND = { sb: 'salt bridge', hb: 'H-bond', pi: 'cation–π' };
  // Confident model: ProtVar pDockQ >= 0.5 (all measured models) or AlphaFold DB ipTM >= 0.8 (0.6 when "acceptable").
  function galPass(g, o) {
    if (g.flags & 8) return false;
    if (g.src === 1 && (g.score == null || g.score < (o.conf === 'acceptable' ? 60 : 80))) return false;
    if (g.plR == null || g.plR < o.pl || g.plP == null || g.plP < o.pl) return false;
    if (o.kind === 'bond' && !(g.flags & 3)) return false;
    if (o.kind === 'sb' && !(g.flags & 1)) return false;
    if (o.kind === 'pi' && !(g.flags & 4)) return false;
    if (o.kind === 'any' && g.d > 4) return false;
    if (o.noPdb && S[g.j].iface) return false;
    if (o.xl && !(xlMask(S[g.j]) & (1 | 4 | 8))) return false;
    return true;
  }

  // Do sites with a best-supported predicted protein contact also contact RNA in cells (crosslinked residue or
  // RBDmap fragment) more or less often than the other methylarginines of the same proteins?
  function galOverlap(o) {
    const opt = Object.assign({}, o, { xl: false });
    const hit = new Set((IX.gallery || []).filter((g) => galPass(g, opt)).map((g) => g.j));
    const measured = new Set((IX.gallery || []).map((g) => P[S[g.j].p].i));
    const tabs = [];
    let a = 0, b = 0, c = 0, d = 0;
    for (const pi of measured) {
      const xp = IX.xlProt.get(pi);
      if (!xp || !(xp.rnaPos.length || xp.rbd.length)) continue;
      const tab = [0, 0, 0, 0];
      for (const j of P[pi].sites) {
        const s = S[j], rna = !!(xlMask(s) & (1 | 4));
        tab[(hit.has(j) ? 0 : 2) + (rna ? 0 : 1)]++;
      }
      tabs.push(tab);
      a += tab[0]; b += tab[1]; c += tab[2]; d += tab[3];
    }
    if (!(a + b) || !(c + d)) return null;
    return { a, b, c, d, pooled: oddsCi(a, b, c, d), p: fisherP(a, b, c, d), mh: mhTest(tabs), proteins: tabs.length };
  }
  const galRank = (a, b) => (b.flags & 1) - (a.flags & 1) || Math.min(b.plR, b.plP) - Math.min(a.plR, a.plP) || (b.score || 0) - (a.score || 0);
  const galBondText = (g) => g.res.filter((r) => r.includes(':')).map((r) => { const [res, codes] = r.split(':'); return ixGene(g.k) + ' ' + res + ' ' + codes.split('+').map((c) => GAL_BOND[c] || c).join(' + '); }).join(', ') || 'nearest ' + ixGene(g.k) + ' ' + g.res[0] + ' at ' + g.d.toFixed(1) + ' Å';
  const galModelText = (g) => (g.src === 0 ? 'ProtVar pDockQ ' + (g.score / 100).toFixed(2) : 'AlphaFold DB ipTM ' + (g.score / 100).toFixed(2));

  function galRows(o) {
    const groups = new Map();
    for (const g of IX.gallery || []) {
      if (!galPass(g, o)) continue;
      const key = g.j + ':' + g.k;
      if (!groups.has(key)) groups.set(key, { j: g.j, k: g.k, models: [] });
      groups.get(key).models.push(g);
    }
    const rows = Array.from(groups.values());
    rows.forEach((r) => r.models.sort(galRank));
    return rows.sort((a, b) => b.models.length - a.models.length || galRank(a.models[0], b.models[0]));
  }

  function galleryHtml() {
    if (!IX.gallery || !IX.gallery.length) return '';
    const o = UI.gal;
    const rows = galRows(o);
    const shown = o.all ? rows : rows.slice(0, 25);
    const nSites = new Set(rows.map((r) => r.j)).size;
    const clashes = IX.gallery.filter((g) => g.flags & 8).length;
    const sel = (id, label, list, cur) => '<label class="fig-ctl">' + label + ' <select id="' + id + '">' + list.map(([v, l]) => '<option value="' + v + '"' + (String(v) === String(cur) ? ' selected' : '') + '>' + esc(l) + '</option>').join('') + '</select></label>';
    const bonds = (g) => g.res.map((r) => { const [res, codes] = r.split(':'); return codes ? '<b>' + esc(res) + '</b> <span class="muted">' + codes.split('+').map((c) => GAL_BOND[c] || c).join(' + ') + '</span>' : '<span class="muted">' + esc(res) + '</span>'; }).slice(0, 3).join('<br>');
    return '<div id="ch-gal"><div class="card-head"><div><p class="sub">Methylarginines that a predicted complex places in a specific contact with the partner: a salt bridge (guanidino N to Asp/Glu carboxylate O ≤ 4 Å) or hydrogen bond (guanidino N to N/O ≤ 3.5 Å), in a confident model, with both residues at pLDDT ≥ the floor and no steric clash. Models: ProtVar served AlphaFold-Multimer models with pDockQ ≥ 0.5, and AlphaFold DB human heterodimers (ipTM ≥ 0.8 when "confident"). Rows where two models agree come first.</p></div>' +
      '<div class="exp"><button type="button" id="gal-tsv">TSV</button></div></div>' +
      '<div class="fig-controls">' + sel('gal-kind', 'Contact', [['bond', 'Salt bridge or H-bond'], ['sb', 'Salt bridge'], ['pi', 'Cation–π (≤ 5 Å)'], ['any', 'Any atom ≤ 4 Å']], o.kind) +
      sel('gal-pl', 'pLDDT ≥', [[70, '70'], [80, '80'], [90, '90']], o.pl) + sel('gal-conf', 'Model', [['confident', 'Confident (pDockQ ≥ 0.5 / ipTM ≥ 0.8)'], ['acceptable', 'Also ipTM 0.6–0.8']], o.conf) +
      '<label class="check fig-ctl"><input type="checkbox" id="gal-nopdb"' + (o.noPdb ? ' checked' : '') + '> Only sites without a PDB interface</label>' +
      '<label class="check fig-ctl" title="RNA-crosslinked residue, RBDmap RNA-bound fragment or DNA-crosslinked residue in living cells"><input type="checkbox" id="gal-xl"' + (o.xl ? ' checked' : '') + '> Only sites also contacting RNA/DNA in cells</label></div>' +
      '<p class="note"><b>' + fmt(rows.length) + '</b> methylarginine–partner contacts on <b>' + fmt(nSites) + '</b> sites pass. <button class="linkish" type="button" id="gal-sites">List these sites</button>. ' + fmt(clashes) + ' measured contacts with a steric clash (< 2.2 Å) are left out.</p>' + (() => {
        const ov = galOverlap(o);
        if (!ov) return '';
        return '<p class="note"><b>Protein and RNA contact at the same arginine.</b> In the ' + fmt(ov.proteins) + ' proteins with both a measured model and RNA crosslinking data, ' + fmt(ov.a) + ' of ' + fmt(ov.a + ov.b) + ' methylarginines with a contact passing these settings (' + pctText(ov.a, ov.a + ov.b || 1) + '%) also contact RNA in cells (crosslinked residue or RBDmap fragment), against ' + fmt(ov.c) + ' of ' + fmt(ov.c + ov.d) + ' (' + pctText(ov.c, ov.c + ov.d || 1) + '%) of their other methylarginines: pooled OR ' + fmtOr(ov.pooled.or, ov.pooled.lo, ov.pooled.hi) + ', Fisher p ' + fmtP(ov.p) + '; within protein OR ' + fmtOr(ov.mh.or, ov.mh.lo, ov.mh.hi) + ', CMH p ' + fmtP(ov.mh.p) + ' (' + fmt(ov.mh.used) + ' informative proteins).</p>';
      })() +
      (rows.length ? '<div class="tbl-wrap"><table class="tbl compact gal-tbl"><thead><tr><th>Site</th><th>Partner</th><th>Contact in the model</th><th>Models</th><th class="r">pLDDT R / partner</th><th class="r hide-xs">Closest</th><th class="hide-sm">PDB</th><th class="hide-sm" title="RNA or DNA contact of this arginine in living cells (crosslinking mass spectrometry)">In cells</th></tr></thead><tbody>' +
        shown.map((r, ri) => {
          const s = S[r.j], p = P[s.p], g = r.models[0], x = IX.prot[r.k];
          return '<tr><td class="nw"><button class="linkish strong" type="button" data-galsite="' + ri + '">' + esc(p.label) + ' ' + rl(s.pos, P[s.p]) + '</button></td>' +
            '<td>' + (x.p >= 0 ? '<button class="linkish" type="button" data-open="' + x.p + '">' + esc(x.gene) + '</button>' : '<b>' + esc(x.gene) + '</b>') + ' ' + partnerTags(x) + '</td>' +
            '<td>' + bonds(g) + '</td>' +
            '<td>' + r.models.map((m, mi) => '<button class="tag pred' + (m.src === 0 ? '' : ' pred-hi') + '" type="button" data-galmodel="' + ri + ':' + mi + '" title="Open this model at ' + rl(s.pos, P[s.p]) + '">' + esc(galModelText(m)) + '</button>').join(' ') + '</td>' +
            '<td class="r num">' + g.plR + ' / ' + g.plP + '</td><td class="r num hide-xs">' + g.d.toFixed(1) + ' Å</td>' +
            '<td class="hide-sm">' + (s.iface ? '<span class="tag iface-p" title="This site is at an interface in a PDB structure">PDB</span>' : '<span class="muted">—</span>') + '</td>' +
            '<td class="hide-sm nw">' + (xlTags(s).replace(/<span class="tag[^"]*xl-near[^"]*"[^>]*>DNA-binding domain<\/span>/, '') || '<span class="muted">—</span>') + '</td></tr>';
        }).join('') + '</tbody></table></div>' + (rows.length > 25 ? '<button class="linkish" type="button" id="gal-more" style="margin-top:.5rem">' + (o.all ? 'Show the top 25' : 'Show all ' + fmt(rows.length)) + '</button>' : '') :
        '<p class="note">No measured contact passes these settings.</p>') +
      '<p class="note">Predicted models of unmodified proteins: a contact suggests where methylation could change binding, not that it does. Distances come from the coordinates of each model; open a model to inspect it in 3D.</p></div>';
  }

  function bindGallery(el) {
    const card = $('#ch-gal', el);
    if (!card) return;
    const o = UI.gal;
    const rows = galRows(o);
    const shown = o.all ? rows : rows.slice(0, 25);
    const redo = () => { const fresh = document.createElement('div'); fresh.innerHTML = galleryHtml(); card.replaceWith(fresh.firstChild); bindGallery(el); };
    $('#gal-kind', card).addEventListener('change', (e) => { o.kind = e.target.value; o.all = false; redo(); });
    $('#gal-pl', card).addEventListener('change', (e) => { o.pl = +e.target.value; redo(); });
    $('#gal-conf', card).addEventListener('change', (e) => { o.conf = e.target.value; redo(); });
    $('#gal-nopdb', card).addEventListener('change', (e) => { o.noPdb = e.target.checked; redo(); });
    $('#gal-xl', card).addEventListener('change', (e) => { o.xl = e.target.checked; o.all = false; redo(); });
    const more = $('#gal-more', card);
    if (more) more.addEventListener('click', () => { o.all = !o.all; redo(); });
    $('#gal-sites', card).addEventListener('click', () => setFilters((f) => { f.galOnly = true; }, 'sites'));
    $$('[data-open]', card).forEach((b) => b.addEventListener('click', () => openProtein(+b.dataset.open)));
    $$('[data-galsite]', card).forEach((b) => b.addEventListener('click', () => { const r = shown[+b.dataset.galsite]; openGalleryModel(r.models[0]); }));
    $$('[data-galmodel]', card).forEach((b) => b.addEventListener('click', () => { const [ri, mi] = b.dataset.galmodel.split(':').map(Number); openGalleryModel(shown[ri].models[mi]); }));
    $('#gal-tsv', card).addEventListener('click', () => download('methylarginine_predicted_contacts.tsv', ['gene\taccession\tsite\tpartner_gene\tpartner_accession\tsource\tscore\tpdockq\tplddt_arginine\tplddt_partner_min\tmin_distance_A\tsalt_bridge\th_bond\tcation_pi\tpartner_residues\tmodel\trna_dna_crosslinks']
      .concat(rows.flatMap((r) => r.models.map((g) => { const s = S[g.j], p = P[s.p], x = IX.prot[g.k];
        return [p.label, p.acc, '' + rl(s.pos, P[s.p]), x.gene, x.acc, g.src === 0 ? 'ProtVar' : 'AlphaFold DB', g.src === 0 ? 'pDockQ ' + (g.score / 100).toFixed(2) : 'ipTM ' + (g.score / 100).toFixed(2), g.pdq == null ? '' : (g.pdq / 100).toFixed(2), g.plR, g.plP, g.d.toFixed(1),
          g.flags & 1 ? 1 : 0, g.flags & 2 ? 1 : 0, g.flags & 4 ? 1 : 0, g.res.join(';'), g.src === 0 ? PROTVAR_MODEL.apply(null, g.ref === 1 ? [p.acc.split('-')[0], x.acc] : [x.acc, p.acc.split('-')[0]]) : 'https://alphafold.ebi.ac.uk/entry/' + g.ref, xlText(s)].join('\t'); }))).join('\n') + '\n'));
  }

  // Open a gallery model on its protein page, focused on the site.
  function openGalleryModel(g) {
    const s = S[g.j], p = P[s.p];
    if (UI.view === 'protein' && UI.detail === p.i) { viewGalleryModel(p, g); return; }
    UI.pendingModel = g;
    openProtein(p.i, s.pos);
  }
  function viewGalleryModel(p, g) {
    const pos = S[g.j].pos;
    if (g.src === 0) return viewProtvar(p, { k: g.k, sc: g.score, first: g.ref === 1 ? 1 : 0, model: 1, sites: [pos], scoreLabel: 'pDockQ ' + (g.score / 100).toFixed(2) + ' (served model)', noTable: true }, null);
    return viewAfdb(p, g.ref, 'ipTM ' + (g.score / 100).toFixed(2) + (g.pdq != null ? ' · pDockQ ' + (g.pdq / 100).toFixed(2) : ''), null, pos);
  }

  // Methylated vs other arginines of the same protein at ProtVar predicted heterodimer interfaces.
  function predEnrichHtml() {
    const pe = IX.predEnrich || {};
    if (!pe['0.23']) return '';
    const cell = (t, f) => {
      const x = pe[t] && pe[t][f];
      if (!x) return '<td class="r">—</td>';
      const [a, b, c, d] = x.pooled;
      return '<td class="r num" title="' + esc('Methylated R at a predicted interface ' + fmt(a) + '/' + fmt(a + b) + ' (' + pctText(a, a + b || 1) + '%) · other R ' + fmt(c) + '/' + fmt(c + d) + ' (' + pctText(c, c + d || 1) + '%) · pooled log2 OR ' + fmtL2(x.pooled[4]) + ' · Mantel–Haenszel over ' + fmt(x.mh[2]) + ' informative proteins') + '">' +
        (x.mh[0] == null ? '—' : '<b>' + fmtL2(x.mh[0]) + '</b> <span class="muted">p ' + fmtP(x.mh[1]) + '</span>') + '</td>';
    };
    const all = pe['0.23'].all.pooled, ord = pe['0.23'].ordered.pooled, dis = pe['0.23'].disordered.pooled;
    return '<h4 style="margin-top:1rem">Predicted interfaces (AlphaFold-Multimer)</h4>' +
      '<p class="sub">PDB structures cover few methylarginines, so the same within-protein question is asked of predicted complexes: ProtVar residues within 5 Å of the partner in AlphaFold-Multimer models of human protein pairs. Models are full length, so disordered regions count, and every arginine of the sequence is compared.</p>' +
      '<div class="tbl-wrap"><table class="tbl compact strata-tbl"><thead><tr><th>Model quality (log2 MH OR)</th><th class="r">All R</th><th class="r">Ordered regions</th><th class="r">Disordered regions</th></tr></thead><tbody>' +
      [['0.23', 'pDockQ ≥ 0.23 (acceptable)'], ['0.5', 'pDockQ ≥ 0.5 (confident)']].map(([t, l]) => '<tr><td><b>' + l + '</b></td>' + cell(t, 'all') + cell(t, 'ordered') + cell(t, 'disordered') + '</tr>').join('') +
      '</tbody></table></div>' +
      '<p class="note">' + fmt(all[0]) + ' of ' + fmt(all[0] + all[1]) + ' methylarginines (' + pctText(all[0], all[0] + all[1] || 1) + '%) in proteins with a modelled partner sit at a predicted interface, against ' + pctText(all[2], all[2] + all[3] || 1) + '% of their other arginines; pooled, that looks like enrichment, but within each protein methylarginines are less often at the predicted interface. ' +
      'In these models ' + pctText(dis[0] + dis[2], dis[0] + dis[1] + dis[2] + dis[3] || 1) + '% of arginines in disordered regions sit at an interface, against ' + pctText(ord[0] + ord[2], ord[0] + ord[1] + ord[2] + ord[3] || 1) + '% in ordered regions. The models are of unmodified sequences and of pairs chosen from interaction data, so this describes the predicted structures, not binding in cells. ' +
      '<button class="linkish" type="button" id="geo-pred">List the ' + fmt(IX.nPred) + ' sites at a predicted interface</button>.</p>';
  }

  function bindIxOverview(el) {
    if (!ixReady()) return;
    $$('[data-cxset]', el).forEach((b) => b.addEventListener('click', () => { UI.cxPred = b.dataset.cxset === 'pred'; UI.cxAll = false; renderView(); }));
    const mn = $('#cx-min', el);
    if (mn) mn.addEventListener('change', () => { UI.cxMin = +mn.value; renderView(); });
    const so = $('#cx-sort', el);
    if (so) so.addEventListener('change', () => { UI.cxSort = so.value; renderView(); });
    const more = $('#cx-more', el);
    if (more) more.addEventListener('click', () => { UI.cxAll = !UI.cxAll; renderView(); });
    bindRows($('#ch-cx', el), (tr) => setTerm('complex', IX.cx[+tr.dataset.cx].id, 'proteins'));
    $$('#ch-readers [data-partner]', el).forEach((b) => b.addEventListener('click', () => { F.ixScore = 700; replaceToken(partnerToken(), 'binds', ixGene(+b.dataset.partner), 'proteins'); }));
    $$('tr[data-iface]', el).forEach((tr) => { const run = () => setFilters((f) => { f.iface = new Set([+tr.dataset.iface]); }, 'sites'); tr.addEventListener('click', run); tr.addEventListener('keydown', (e) => { if (e.key === 'Enter') run(); }); });
    const gm = $('#geo-modelled', el);
    if (gm) gm.addEventListener('click', () => setFilters((f) => { f.modelledPdb = true; }, 'sites'));
    const gp = $('#geo-pred', el);
    if (gp) gp.addEventListener('click', () => setFilters((f) => { f.pred = 23; }, 'sites'));
    bindGallery(el);
    $$('tr[data-xl]', el).forEach((tr) => { const run = () => setFilters((f) => { f.xl = new Set([+tr.dataset.xl]); }, 'sites'); tr.addEventListener('click', run); tr.addEventListener('keydown', (e) => { if (e.key === 'Enter') run(); }); });
    $$('tr[data-xlcls]', el).forEach((tr) => {
      const run = () => setFilters((f) => { f.xl = new Set([1]); if (tr.dataset.xlcls === 'rme2s') f.ips = new Set([M.maron_ips.indexOf('SDMA')]); }, 'sites');
      tr.addEventListener('click', run); tr.addEventListener('keydown', (e) => { if (e.key === 'Enter') run(); });
    });
    $$('tr[data-go]', el).forEach((tr) => { const run = () => { F.goExp = tr.dataset.goexp === '1'; setTerm('go', tr.dataset.go === '1' ? 'rna' : 'dna', 'proteins'); }; tr.addEventListener('click', run); tr.addEventListener('keydown', (e) => { if (e.key === 'Enter') run(); }); });
    const tsv = $('#cx-tsv', el);
    if (tsv) tsv.addEventListener('click', () => download('complexes_methylarginine_members.tsv', ['complex_ac\tname\tset\tmembers\tmethylated_members\tfold\tp\tq\tmethylated_genes']
      .concat(IX.cx.filter((c) => c.m > 0).map((c) => [c.id, c.name, c.pred ? 'predicted' : 'curated', c.n, c.m, c.fold, c.pv, c.q,
        c.members.map((m) => (Array.isArray(m) ? m : [m])).filter((ids) => ids.some((k) => k < IX.nP)).map((ids) => ids.map(ixGene).join('/')).join(';')].join('\t'))).join('\n') + '\n'));
  }

  // ----- nucleic acids in the 3D view
  const NUC_RNA = ['A', 'C', 'G', 'U', 'I', 'PSU', '5MC', '7MG', 'OMG', 'OMC', 'OMU', '2MG', 'M2G', 'H2U', '1MA', '5MU', '4SU'];
  const NUC_DNA = ['DA', 'DC', 'DG', 'DT', 'DI', 'DU', '5CM'];
  const NUC_ALL = new Set(NUC_RNA.concat(NUC_DNA));
  const isNuc = (a) => NUC_ALL.has(a.resn);
  const nucKind = (a) => (NUC_DNA.includes(a.resn) ? 'DNA' : 'RNA');

  // Methylarginines (author numbering in the mapped chain) within 4 Å of nucleic acid or of another protein chain.
  function pdbContacts() {
    const { viewer, pdb } = mol;
    const out = [];
    const sideChain = ['NE', 'CZ', 'NH1', 'NH2', 'CD'];
    for (const [pos, auth] of pdb.toAuth) {
      const atoms = viewer.selectedAtoms({ chain: pdb.chain, resi: auth, atom: sideChain });
      if (!atoms.length) continue;
      const near = viewer.selectedAtoms({ within: { distance: 4, sel: { chain: pdb.chain, resi: auth, atom: sideChain } }, predicate: (a) => a.chain !== pdb.chain || isNuc(a) });
      const hits = new Map();
      for (const a of near) {
        if (WATER.includes(a.resn)) continue;
        const kind = isNuc(a) ? nucKind(a) : a.hetflag && !isPolymer(a) ? null : 'protein';
        if (!kind) continue;
        const key = kind + ':' + a.chain + ':' + a.resi;
        if (!hits.has(key)) hits.set(key, { kind, chain: a.chain, resi: a.resi, resn: a.resn });
      }
      if (hits.size) out.push({ pos, auth, hits: Array.from(hits.values()), b: atoms[0].b });
    }
    return out;
  }

  // Predicted complexes: each methylarginine's nearest partner heavy atom in the loaded model.
  function predDistances() {
    const { viewer, pdb } = mol;
    const out = [];
    for (const [pos] of pdb.toAuth) {
      const own = viewer.selectedAtoms({ chain: pdb.chain, resi: pos });
      if (!own.length) continue;
      const near = viewer.selectedAtoms({ within: { distance: 8, sel: { chain: pdb.chain, resi: pos } }, predicate: (a) => a.chain !== pdb.chain });
      let best = null;
      for (const a of near) for (const o of own) {
        const d = Math.hypot(a.x - o.x, a.y - o.y, a.z - o.z);
        if (!best || d < best.d) best = { d, a };
      }
      const ca = own.find((o) => o.atom === 'CA') || own[0];
      out.push({ pos, d: best ? best.d : null, partner: best ? best.a : null, b: ca.b });
    }
    return out.sort((x, y) => (x.d == null ? 99 : x.d) - (y.d == null ? 99 : y.d));
  }

  function predDistHtml() {
    const pred = mol.pdb.pred;
    if (!mol.pdb.dist) mol.pdb.dist = predDistances();
    const list = mol.pdb.dist;
    if (!list.length) return '';
    const table = new Set(pred.tableSites || []);
    return '<p class="flabel">Distance to the partner in this model</p><ul class="contact-list">' + list.map((x) => '<li><button class="linkish strong" type="button" data-pos="' + x.pos + '">' + rl(x.pos) + '</button> ' +
      (x.d == null ? '<span class="muted">more than 8 Å from ' + esc(pred.partnerLabel) + '</span>' :
        '<b>' + x.d.toFixed(1) + ' Å</b> to ' + esc((pred.names[x.partner.chain] || 'chain ' + x.partner.chain) + ' ' + x.partner.resn + x.partner.resi) + (x.d < 2.2 ? ' <span class="tag clash" title="Heavy atoms closer than 2.2 Å: a steric clash in the model, not a real contact">clash</span>' : x.d <= 4 ? ' <span class="tag iface-p">contact</span>' : '')) +
      ' <span class="muted">· pLDDT ' + Math.round(x.b) + '</span>' +
      (table.has(x.pos) && (x.d == null || x.d > 5) ? '<br><span class="muted">ProtVar\'s 2024 interface table puts ' + rl(x.pos) + ' within 5 Å; this later model does not.</span>' : '') + '</li>').join('') + '</ul>' +
      '<p class="note">Nearest heavy atoms between the arginine and any partner chain; ≤ 4 Å is marked as a contact, < 2.2 Å as a steric clash (a modelling error).' + (pred.tableSites ? ' The ProtVar API serves a later prediction run than the 2024 interface table used for the site flags and statistics, so distances can differ.' : '') + '</p>';
  }

  function contactsHtml() {
    if (!mol || mol.mode !== 'pdb') return '';
    if (!mol.pdb.contacts) mol.pdb.contacts = pdbContacts();
    const list = mol.pdb.contacts;
    const hasNa = mol.viewer.selectedAtoms({ predicate: isNuc }).length > 0;
    if (!list.length) return '<p class="flabel">Contacts</p><p class="note">No resolved methylarginine side chain is within 4 Å of ' + (hasNa ? 'nucleic acid or ' : '') + 'another chain.</p>';
    return '<p class="flabel">Methylarginines in contact (≤ 4 Å)</p><ul class="contact-list">' + list.map((c) => {
      const kinds = ['RNA', 'DNA', 'protein'].map((k) => {
        const h = c.hits.filter((x) => x.kind === k);
        if (!h.length) return '';
        const pred = mol.pdb.pred;
        const lab = h.slice(0, 4).map((x) => (pred ? (pred.names[x.chain] || x.chain) + ' ' + x.resn + x.resi : x.resn + ' ' + x.chain + x.resi)).join(', ') + (h.length > 4 ? ' +' + (h.length - 4) : '');
        return '<span class="tag ' + (k === 'RNA' ? 'na-rna' : k === 'DNA' ? 'na-dna' : 'iface-p') + '" title="' + esc(lab) + '">' + k + '</span> <span class="muted">' + esc(lab) + '</span>';
      }).filter(Boolean).join('<br>');
      return '<li><button class="linkish strong" type="button" data-pos="' + c.pos + '">' + rl(c.pos) + '</button> <span class="muted">(' + (mol.pdb.pred ? 'pLDDT ' + Math.round(c.b) : esc(mol.pdb.chain) + c.auth) + ')</span><br>' + kinds + '</li>';
    }).join('') + '</ul>';
  }

  // ------------------------------------------------------------- print / PDF
  // A protein page prints as a report; the browser's "Save as PDF" keeps text, tables and the map as vectors.
  // Before printing, the map is redrawn in the light export palette and the WebGL view is replaced by a snapshot.
  let printUndo = [];
  let printReady = false;   // set when Save PDF prepared the page itself, so beforeprint does not redo it
  function restorePrint() {
    printReady = false;
    const undo = printUndo;
    printUndo = [];
    undo.reverse().forEach((f) => { try { f(); } catch (e) { /* element gone */ } });
  }
  function preparePrint() {
    restorePrint();
    const p = UI.view === 'protein' ? P[UI.detail] : null;
    if (!p) return;
    setPrintPageStyle(true, p.label + ' (' + p.acc + ') · ' + new Date().toISOString().slice(0, 10));
    printUndo.push(() => setPrintPageStyle(false));
    const head = $('#print-head');
    if (head) head.innerHTML = '<b>' + esc(siteName()) + '</b> · Shechter Lab · ' + esc(location.href) + ' · printed ' + new Date().toISOString().slice(0, 10);
    const box = $('#map-box');
    if (box) {
      const div = document.createElement('div');
      div.className = 'print-only print-map';
      div.innerHTML = mapSvg(p, UNIDATA.get(p.acc), { pal: EXPORT_PAL, forExport: true }).svg;
      box.after(div);
      printUndo.push(() => div.remove());
    }
    const netBox = $('#ix-net');
    if (netBox && UI.ixNet && UI.ixNet.p === p) {
      const div = document.createElement('div');
      div.className = 'print-only print-net';
      div.innerHTML = ixNetSvg(UI.ixNet, EXPORT_PAL, true);
      netBox.after(div);
      printUndo.push(() => div.remove());
    }
    if (mol && mol.p.i === p.i && mol.el && document.body.contains(mol.el)) {
      try {
        const v = mol.viewer;
        v.setBackgroundColor(UI.mv.bg === 'black' ? '#000000' : '#FFFFFF');
        v.render();
        const fig = document.createElement('figure');
        fig.className = 'print-only print-3d';
        fig.innerHTML = '<img alt="3D structure snapshot" src="' + v.pngURI() + '"><figcaption>' +
          esc(mol.mode === 'pdb' ? 'PDB ' + mol.pdb.e.id.toUpperCase() + ', chain ' + mol.pdb.chain : 'AlphaFold DB model ' + ((AFDATA.get(p.acc) || {}).entry || {}).entryId) +
          ' as shown on screen; methylarginine side chains in the protein-map colours.</figcaption>';
        v.setBackgroundColor(viewBackground());
        v.render();
        mol.el.after(fig);
        printUndo.push(() => fig.remove());
      } catch (e) { /* no snapshot; the viewer is hidden in print */ }
    }
    $$('#content details:not([open])').forEach((d) => { d.open = true; printUndo.push(() => { d.open = false; }); });
    const title = document.title;
    document.title = (p.gene || p.acc) + '_' + p.acc + '_methylarginines';
    printUndo.push(() => { document.title = title; });
  }
  window.addEventListener('beforeprint', () => { if (!printReady) preparePrint(); });
  window.addEventListener('afterprint', restorePrint);

  // Load what the page fetches lazily (UniProt, AlphaFold, the 3D view) before opening the print dialog.
  function savePdf(p, btn) {
    btn.disabled = true;
    btn.textContent = 'Preparing…';
    if (UI.startStructure) UI.startStructure();
    const settle = (pr, ms) => Promise.race([pr.catch(() => null), new Promise((r) => setTimeout(r, ms))]);
    const viewerReady = () => new Promise((resolve) => {
      const t0 = Date.now();
      (function poll() {
        const af = AFDATA.get(p.acc);
        if ((mol && mol.p.i === p.i) || af === null || Date.now() - t0 > 15000) resolve();
        else setTimeout(poll, 200);
      })();
    });
    settle(Promise.all([loadUniProt(p.acc), loadAF(p.acc), loadEvidence()]), 15000)
      .then(viewerReady)
      .then(() => new Promise((r) => setTimeout(r, 500)))
      .then(() => {
        if (UI.detail !== p.i) return null;
        preparePrint();
        printReady = true;
        // The snapshot is a data URI; wait for it to decode so it is in the printout.
        const img = $('.print-3d img');
        return img && img.decode ? settle(img.decode(), 5000) : null;
      })
      .then(() => {
        btn.disabled = false;
        btn.textContent = 'Save PDF';
        if (UI.detail === p.i && printReady) {
          window.print();
        }
      });
  }

  function bindMap(box, p, u, g) {
    const svg = $('svg', box);
    const brush = $('#brush', svg);
    const toX = (e) => { const r = svg.getBoundingClientRect(); return ((e.clientX - r.left) / r.width) * g.W; };
    const toPos = (vx) => g.v0 + ((vx - g.ML) / (g.W - g.ML - g.MR)) * (g.v1 - g.v0 + 1);
    let start = null;
    svg.addEventListener('pointerdown', (e) => {
      if (e.target.classList.contains('head') || e.button !== 0) return;
      start = toX(e);
      if (start < g.ML) { start = null; return; }
      svg.setPointerCapture(e.pointerId);
      brush.removeAttribute('hidden');
      brush.setAttribute('x', start); brush.setAttribute('width', 0);
    });
    svg.addEventListener('pointermove', (e) => {
      if (start == null) return;
      const cur = Math.max(g.ML, Math.min(g.W - g.MR, toX(e)));
      brush.setAttribute('x', Math.min(start, cur)); brush.setAttribute('width', Math.abs(cur - start));
    });
    svg.addEventListener('pointerup', (e) => {
      if (start == null) return;
      const cur = Math.max(g.ML, Math.min(g.W - g.MR, toX(e)));
      const a = Math.round(toPos(Math.min(start, cur))), b = Math.round(toPos(Math.max(start, cur)));
      start = null;
      brush.setAttribute('hidden', '');
      if (b - a >= 8) { UI.zoom = [a, b]; drawMap(p, u); }
    });
    $$('.head', svg).forEach((h) => {
      const s = S[+h.dataset.j];
      h.addEventListener('mouseenter', (e) => showTip(e, siteTip(s)));
      h.addEventListener('mousemove', moveTip);
      h.addEventListener('mouseleave', hideTip);
      h.addEventListener('click', () => { UI.highlight = s.pos; flashSite(s.pos, true); drawMap(p, u); });
    });
    $$('.feat', svg).forEach((f) => {
      f.addEventListener('mouseenter', (e) => showTip(e, f.getAttribute('data-tip')));
      f.addEventListener('mousemove', moveTip);
      f.addEventListener('mouseleave', hideTip);
    });
  }

  function siteTip(s) {
    const ipr = s.ipr >= 0 ? M.interpro[s.ipr][1] : null;
    return '<b>' + rl(s.pos, P[s.p]) + '</b> · ' + STATE_LABEL[s.st] + (s.ip ? ' · ' + bits(s.ip).map((b) => ipLabel(M.maron_ips[b])).join('/') + ' IP' : '') + '<br>' +
      (s.win ? windowHtml(s.win, 7) + '<br>' : '') +
      esc(DISORDER_LABEL[M.disorder_classes[s.dc]] || '') + ' · ' + esc(DOMAIN_LABEL[M.domain_contexts[s.dctx]] || '') + (ipr ? ' (' + esc(ipr) + ')' : '') + '<br>' +
      'Reported by ' + M.studyOrder.filter((k) => s.study & (1 << k)).map((k) => esc(M.studies[k][0])).join(', ') +
      (s.trf ? '<br>' + [0, 1].map((t) => (s.trf & (1 << (t * 2)) ? '▼' + M.maron_treatments[t] + ' ' : '') + (s.trf & (1 << (t * 2 + 1)) ? '▲' + M.maron_treatments[t] + ' ' : '')).join('') : '');
  }

  // ------------------------------------------------------------------ tooltip
  const tipEl = () => $('#tip');
  function showTip(e, html) { const t = tipEl(); t.innerHTML = html; t.hidden = false; moveTip(e); }
  function moveTip(e) {
    const t = tipEl();
    if (t.hidden) return;
    const pad = 14, w = t.offsetWidth, h = t.offsetHeight;
    let left = e.clientX + pad, topY = e.clientY + pad;
    if (left + w > window.innerWidth - 8) left = e.clientX - w - pad;
    if (topY + h > window.innerHeight - 8) topY = e.clientY - h - pad;
    t.style.left = Math.max(8, left) + 'px'; t.style.top = Math.max(8, topY) + 'px';
  }
  function hideTip() { const t = tipEl(); if (t) t.hidden = true; }

  // ---------------------------------------------------------------- about
  function renderAbout(el) {
    const sec = (id, title, html) => '<section class="doc-sec" id="doc-' + id + '"><h3>' + esc(title) + '</h3>' + html + '</section>';
    const sub = (id, title, html) => '<div class="doc-sub" id="doc-' + id + '"><h4>' + esc(title) + '</h4>' + html + '</div>';
    const parts = [
      ['start', 'Getting started', '<ul>' +
        '<li><b>Search</b> for a gene symbol or synonym (FUS, TDP43), a UniProt accession, part of a protein name, a site (FUS R216) or a pasted gene list. The bar under the search box lists every active filter; click one to remove it.</li>' +
        '<li><b>Filter</b> with the panel on the left, or click a bar in any chart. The URL keeps the search and filters, so a link reopens the same selection.</li>' +
        '<li><b>Views</b>: Overview (dataset charts), Proteins and Sites (tables, TSV download), Interactions (partners, complexes, structures and RNA/DNA crosslinks), Figures (alignment, logo and protein comparison from the current selection) and a page per protein (map, sites, interactions, conservation, 3D structure and sequence).</li>' +
        '<li><b>Human only.</b> Every site, protein and annotation is <i>Homo sapiens</i> (canonical UniProt sequences, human cells and tissues). Related proteins from other species appear only in the homologous-structures panel and the orthologue alignment on protein pages, labelled with their species.</li>' +
        '<li><b>Residue numbering</b> follows the canonical UniProt sequence with the initiator methionine as 1, except for histones (H1, H2A, H2B, H3, H4 and their variants, including CENP-A), which follow the histone convention without the initiator methionine: H3R2 is UniProt R3 of P68431. Histone sites are searched, shown and linked in histone numbering, with the UniProt number alongside; stored positions and the TSV <code>position</code> column stay in UniProt numbering.</li>' +
        '<li><b>Where you are.</b> The address always names the view, or the protein and the open site (<code>#P35637/R216</code>), together with the search and filters, so a copied link reopens the same place; protein pages show breadcrumbs back to the list you came from and the proteins you opened recently.</li></ul>' +
        sub('search', 'Search terms', '<p>Filters are written into the search as readable terms, so the search box and the link always say what is selected. Clicking a chart bar, a card or a filter chip adds the same term; type them directly too, with suggestions as you type. ' +
          'A comma means any of the values (<code>motif:RGG,GAR</code>); separate terms must all match (<code>class:rna-binding interaction:experimental</code>); plain words still search genes and sites (<code>FUS motif:RGG</code>). Names can be shortened to their first word (<code>disease:als</code>, <code>condensate:stress</code>) or quoted in full.</p>' +
          '<table class="tbl compact doc-terms"><thead><tr><th>Term</th><th>Selects</th></tr></thead><tbody>' + TERM_KEYS.map((k) => '<tr><td><code>' + k + ':</code></td><td>' + esc(TERM_HELP[k]) + '</td></tr>').join('') + '</tbody></table>' +
          '<p>Links made before terms existed (numeric <code>motifs=</code>, <code>disease=</code>, <code>func=</code> … parameters) still open and are rewritten as terms.</p>') +
        sub('tables', 'Tables, columns and the site card', '<p>The Proteins, Sites and protein-page tables show their core columns and a <b>Columns</b> menu for the rest (context, next methyl-R, Maron IP, AlphaFold, nearest domain …); the choice is remembered in this browser. On a phone each row becomes a labelled card. ' +
          'Clicking a site, or opening a link with a site in the search (<code>FUS R216</code>), opens the <b>site card</b> under the protein map: its interaction status and evidence with partners, what is known about the whole protein (labelled as such), sequence context, inhibitor response, AlphaFold confidence and reports.</p>')],
      ['overview', 'Overview charts', '<p>Bars show log2 odds ratios; filled bars have a Benjamini–Hochberg q ≤ 0.05, outlined bars do not. The charts always describe the whole dataset; clicking a bar filters the Proteins and Sites views.</p>' +
        sub('ov-ontology', 'Functional ontology', '<p>Classes from the lab functional ontology (config/functional_ontology.json); a protein can belong to several. Each class is compared between methylated proteins and the reviewed human proteome. "Show subclasses" lists subclasses under the broad class most of their proteins share.</p>') +
        sub('ov-disease', 'Disease families', '<p>Audited DisGeNET disease families (via DOSE/clusterProfiler), genes counted once per family, BH-adjusted across the 8 families. This is enrichment of disease gene sets among methylated proteins, not causation. The ALS/MN/FTD filter uses the full audited family gene set.</p>') +
        sub('ov-motifs', 'Sequence motifs', '<p>Motifs are read from a ±7 window on the canonical sequence and compared between methylated arginines and the other arginines of the same proteins: ' + M.motifs.map((m) => esc(MOTIF_LABEL[m]) + ' (' + esc(MOTIF_HELP[m]) + ')').join('; ') + '.</p>') +
        sub('ov-condensates', 'Condensates', '<p>CD-CODE condensates with at least five counted members; the bar is the log2 odds ratio of methylated members against the global rate, and the numbers are methylated / members. "Count members" keeps only memberships with CD-CODE confidence ≥ 2, 3 or 4 of 5, so condensates whose member lists rest on single proteomics screens do not dominate. ' + condHistNote() + '</p>') +
        sub('ov-disorder', 'Disorder context', '<p>MobiDB disorder intervals: ' + Object.keys(DISORDER_LABEL).map((k) => esc(DISORDER_LABEL[k]) + ' = ' + esc(lc1(DISORDER_HELP[k] || ''))).join('; ') + '. Methylated arginines are compared with the other arginines of the same proteins. Sites on proteins without a MobiDB annotation are left out because their background has only a few dozen arginines.</p>') +
        sub('ov-domains', 'Domains', '<p>Position relative to InterPro domains and repeats: ' + Object.keys(DOMAIN_LABEL).map((k) => esc(DOMAIN_LABEL[k]) + ' = ' + esc(lc1(DOMAIN_HELP[k]))).join('; ') + '. The nearest-domain class comes from the lab domain-class ontology (config/domain_class_ontology.json).</p>') +
        sub('ov-provenance', 'Provenance', '<p>Sites reported by each study compiled in Maron et al. 2021 Table S5, by ProMetheusDB and by curated UniProt/IEDB records; the column gives the sites reported only there. A report is one study or database.</p>') +
        sub('ov-inhibitors', 'PRMT inhibitor response', '<p>Maron et al. 2021 Table S3: sites per methyl-IP that changed significantly (adjusted p ≤ 0.05) after GSK591 (PRMT5 inhibitor) or MS023 (type I PRMT inhibitor) compared with DMSO. The median change is the Table S3 <code>&lt;inhibitor&gt;_vs_DMSO_diff</code> value across sites.</p>')],
      ['figures', 'Figures and alignments', '<p>Every exported figure (SVG or PNG, including 3D snapshots and the interaction network) ends with a credit line naming this site, its address and ' + esc(PAPER.short) + '. Printed protein pages carry the same.</p>' +
        sub('fig-align', 'Site alignment and sequence logo', '<p>The selected sites (Figures tab) are aligned at the methylated arginine using the ±15 residue windows of the canonical sequence; positions past a protein end are left blank. The logo uses every selected site with a window. "Frequency" stacks residue frequencies at each position; "Information" scales each stack by its information content, log2(20) minus the Shannon entropy and the small-sample correction 19/(2 ln2 n) (Schneider et al. 1986, ' + doiLink('10.1016/0022-2836(86)90165-8', 'J Mol Biol 188:415') + '), so the invariant central arginine reaches about 4.3 bits. Colours follow WebLogo\'s chemistry scheme (Crooks et al. 2004, ' + doiLink('10.1101/gr.849004', 'Genome Res 14:1188') + '). FASTA export writes the aligned windows with - for positions past a protein end.</p>') +
        sub('fig-compare', 'Protein comparison', '<p>The selected proteins (a searched gene list keeps its order) drawn on one residue scale or each to its own length, with UniProt domains, repeats and zinc fingers, UniProt "Disordered" regions, and the selected methylarginines coloured as on protein maps. Features are fetched live from UniProt.</p>') +
        sub('fig-ortho', 'Orthologue alignment', '<p>On protein pages, the Ensembl gene from UniProt\'s cross-references is sent to Ensembl Compara (REST, live) for orthologues in mouse, chicken, frog (<i>X. tropicalis</i>), zebrafish, fly, worm and yeast; for each species the orthologue with the highest identity is used, with Ensembl\'s pairwise protein alignment. UniProt positions are carried to the Ensembl human protein when the sequences are identical or one contains the other, otherwise through a unique 15-residue anchor around each site; positions without one are shown as not mapped. Columns follow the human sequence, so insertions in other species are not shown. R means the arginine is conserved, K a conservative change. ' + doiLink('10.1093/database/bav096', 'Herrero et al. 2016') + '; ' + doiLink('10.1093/bioinformatics/btu613', 'Yates et al. 2015') + '.</p>')],
      ['insight', 'Evidence map, selection profile and report', sub('ix-status', 'Interaction status of a site', '<p>Every site has one of three statuses, shown as a badge in the tables, on the site card, in the map colouring (<i>Colour sites by: Interaction evidence</i>) and as the quick filter above the Sites table. ' +
          '<b>Experimental</b>: an RNA crosslink at the residue (pepR-MS compendium), an RBDmap RNA-bound fragment, a DNA crosslink localised to or with a candidate at the residue (Trendel et al. 2025), or a PDB interface with protein, RNA or DNA (PDBe-KB). ' +
          '<b>Predicted only</b>: at the interface of an AlphaFold-Multimer model of a human protein pair (ProtVar, pDockQ ≥ 0.23) or a confident predicted contact (salt bridge or hydrogen bond, both residues pLDDT ≥ 80), with no experimental evidence. <b>None recorded</b>: neither. ' +
          'Dashed tags (near an RNA crosslink, inside a DNA-crosslinked peptide, a UniProt DNA-binding region) are shown for context but do not set the status. Most proteins were never crystallised or crosslinked, so "none recorded" is not evidence of no contact. ' +
          'Partner and complex information about the whole protein (STRING, Complex Portal, Gene Ontology) is never counted as site evidence. Search with <code>interaction:experimental</code>, <code>interaction:predicted-only</code>, <code>interaction:none</code>, or <code>interaction:predicted</code> for any predicted evidence. ' +
          'The <b>Check a site</b> box on the Overview and Interactions pages answers the question for one site without leaving the page.</p>') +
        sub('ix-evmap', 'How the evidence fits together', '<p>Counted per site or per protein. For every methylarginine in the current selection, which interaction evidence it carries: crosslinked to RNA at the residue (pepR-MS compendium), inside an RBDmap RNA-bound fragment, a DNA crosslink at or possibly at the residue (Trendel et al. 2025), a protein or RNA/DNA interface in a PDB structure (PDBe-KB), the interface of an AlphaFold-Multimer model (ProtVar, pDockQ ≥ 0.23), and a salt bridge or hydrogen bond in a confident predicted complex. Rows are the combinations that occur, largest first (an UpSet layout). A column total lists every site with that type; a row lists the sites with exactly that combination. Sites without any type are mostly in proteins never crystallised, modelled with a partner or crosslinked. The protein-context table shows, for the sites with each kind of evidence, the share on proteins annotated as RNA- or DNA-binding (Gene Ontology), in a CD-CODE condensate, physically linked to a PRMT or a Tudor-domain reader (STRING ≥ 0.7) or in a Complex Portal curated complex, against all selected sites; each cell lists those sites. Counted per protein, the same layout combines these protein properties and whether any of the protein\'s sites carries contact evidence. Below the site view, the proteins with most sites carrying a chosen type, with links to list the sites, draw their logo and alignment, or compare the proteins.</p>') +
        sub('profile', 'What sets this selection apart', '<p>Shown on the Proteins and Sites views whenever a search or filter is set. Each site feature (motifs, disorder, methylation calls and IPs, PRMT inhibitor response, close spacing, number of reports and the interaction evidence types) is compared between the selected sites and all other sites with Fisher\'s exact test, Benjamini–Hochberg-adjusted across the features tested; features with q ≤ 0.05 are listed, more common ones first. The ratio is the share in the selection over the share in the other sites. Clicking a more common feature narrows the selection to it. Sites of the same protein are not independent, so this points to differences rather than testing them; a selection built from a few proteins mostly reflects those proteins.</p>') +
        sub('compare', 'Comparing two selections', '<p>Press <b>Save as A for a comparison</b> on a selection (on the selection profile or the report), change the search or filters, then <b>Compare with A</b> (also in the bar under the search). The comparison shows both selections, their overlap, a paragraph on each, the interaction status split, every site feature as the share of A and of B with the ratio and a Fisher\'s exact test (Benjamini–Hochberg across features), side-by-side logos and the proteins with most sites. ' +
          'Sites in both selections count on both sides, which makes the test conservative. The link keeps both selections (<code>A=</code> holds the first), and the page prints with the same header, citation and page numbers as the report.</p>') +
        sub('report', 'Insight report', '<p>One printable page for the current selection: the selection and a link that reopens it, headline numbers, what sets it apart, the site and protein evidence maps with the protein-context table, the sequence logo and alignment (±7, first 30 sites), the proteins on one residue scale (up to 8, UniProt features fetched live), the proteins with most selected sites with their evidence, and the methods and citation. It opens the browser\'s print dialog; choose Save as PDF.</p>')],
      ['about', 'About', '<p>The explorer indexes the integrated human methylarginine union assembled by the <a href="https://github.com/Shechterlab/PTM_Pipeline">PTM_Pipeline</a> for the review below. It runs entirely in your browser; UniProt, AlphaFold DB, PDB, Ensembl and other annotations on protein pages are fetched live.</p>' +
        '<p><b>Cite:</b> ' + esc(PAPER.cite) + ' <a href="https://doi.org/' + PAPER.doi + '">doi:' + PAPER.doi + '</a>, and the original studies and resources (see How to cite).</p>'],
      ['cite', 'How to cite', '<p>The sites here were measured and first reported by other groups; this explorer compiles them. When you use sites or evidence from it, please cite:</p><ol>' +
        '<li><b>The original studies or databases that reported the sites you use.</b> Every list of sites (Proteins, Sites, protein pages, the insight report and the comparison) ends with <i>How to cite these sites</i>: the reporting studies with how many of the sites each reported, full references and DOIs, and a <i>Copy references</i> button. The site card lists the studies behind one site, and the Sites TSV has a <code>reported_by_doi</code> column.</li>' +
        '<li><b>The resources behind any evidence you rely on</b>: for example the pepR-MS compendium and the datasets it re-analyses for RNA crosslinks, Trendel et al. 2025 for DNA crosslinks, PDBe-KB for interfaces, ProtVar and AlphaFold for predicted complexes, CD-CODE for condensates and STRING or Complex Portal for partners. The same lists name them for the sites in view; every resource is in <i>Data sources, versions and licences</i> below with its citation and licence.</li>' +
        '<li><b>' + esc(PAPER.cite) + '</b> ' + doiLink(PAPER.doi) + ', for the compiled, remapped union and this explorer.</li></ol>' +
        '<p>PhosphoSitePlus® sites carry a non-commercial licence (CC BY-NC-SA 3.0); see Sources in the union.</p>' +
        '<h4>Studies and databases behind the union</h4><ol class="cite-list">' + M.studyOrder.map((k) => '<li>' + refLine(studyRef(k)) + ' <span class="muted">· ' + fmt(M.studies[k][2]) + ' sites</span></li>').join('') + '</ol>'],
      ['sources', 'Sources in the union', '' +
      '<dl>' +
      '<dt>Maron 2021</dt><dd>Table S5 of Maron et al., <a href="https://doi.org/10.1016/j.isci.2021.102971">iScience 24:102971 (2021)</a>, a compilation of human methylarginine sites from published proteomics.</dd>' +
      '<dt>ProMetheusDB</dt><dd>Methyl-proteome table (mmc2) with monomethyl (Rme1) and dimethyl (Rme2) calls: Massignani et al., ' + doiLink('10.1016/j.mcpro.2022.100243', 'Mol Cell Proteomics 21:100243 (2022)') + '.</dd>' +
      '<dt>PTM master</dt><dd>Methylarginine rows of the lab\'s cross-PTM human master table (UniProt and IEDB annotations); the same table supplies the phospho, lysine acyl, ubiquitin/SUMO and lysine methyl sites shown as crosstalk. Those other-PTM sites are not exhaustive: they come from an iPTMnet export (HPRD, literature mining and other curated sources) plus UniProt, so modifications reported elsewhere, for example in PhosphoSitePlus or recent proteomics studies, may be missing. iPTMnet: Huang et al., ' + doiLink('10.1093/nar/gkx1104', 'Nucleic Acids Res 46:D542 (2018)') + '; IEDB: Vita et al., ' + doiLink('10.1093/nar/gky1006', 'Nucleic Acids Res 47:D339 (2019)') + '.</dd></dl>' +
      '<p class="note"><b>Licence notice.</b> ' + fmt((M.studies.find((x) => x[0] === 'PhosphoSitePlus') || [0, 0, 0])[2]) + ' sites in the union are reported by PhosphoSitePlus® (via Maron et al. 2021 Table S5). ' + esc(PSP_NOTICE) + ' Please acknowledge PhosphoSitePlus® and cite Hornbeck et al. 2015 (' + doiLink('10.1093/nar/gku1267') + ') when using these sites, including in downloads.</p>' +
      '<p>Studies and databases behind the sites (sites reported · reported only there):</p><ul>' +
      M.studyOrder.map((k) => '<li><b>' + esc(M.studies[k][0]) + '</b> · ' + fmt(M.studies[k][2]) + ' · ' + fmt(M.studies[k][3]) + (STUDY_DOI[M.studies[k][0]] ? ' · ' + doiLink(STUDY_DOI[M.studies[k][0]]) : '') + '<br><span class="muted">' + esc(M.studies[k][1] && M.studies[k][1] !== 'This study.' ? M.studies[k][1] : STUDY_HELP[M.studies[k][0]] || '') + '</span></li>').join('') + '</ul>' +
      ''],
      ['processing', 'How sites were processed', '' +
      '<ul>' +
      '<li>Source rows were remapped to canonical UniProt sequences by peptide/sequence anchors and deduplicated by canonical accession and position. The original reported accession and position are kept in each site\'s evidence.</li>' +
      '<li><b>Gene symbols</b> follow current UniProt. Symbols reported by the sources (some outdated, a few belonging to paralogs with near-identical sequence, such as SNRPN on SmB/B′ rows) and UniProt synonyms are kept as searchable aliases.</li>' +
      '<li><b>State</b> is the consensus Rme1/Rme2 call where a source provides one; most compiled sites have no state call. Maron Table S2 Rme1/Rme2 calls and Table S3 methyl-IP recovery are shown separately; the IPs are named by the state their antibody recognises: Rme1 (MMA antibody), Rme2a (ADMA antibody) and Rme2s (SDMA antibody).</li>' +
      '<li><b>PRMT inhibitor response</b> uses Maron Table S3 <code>GSK591_vs_DMSO</code> and <code>MS023_vs_DMSO</code> differences (median over peptides) with adjusted p ≤ 0.05. GSK591 inhibits PRMT5; MS023 inhibits type I PRMTs.</li>' +
      '<li><b>Motifs</b> are read from a ±7 window on the canonical sequence: ' + M.motifs.map((m) => MOTIF_LABEL[m] + ' (' + MOTIF_HELP[m] + ')').join('; ') + '.</li>' +
      '<li><b>Disorder</b> uses MobiDB intervals: ' + Object.keys(DISORDER_LABEL).map((k) => DISORDER_LABEL[k] + ' = ' + lc1(DISORDER_HELP[k] || '')).join('; ') + '. The architecture figure draws UniProt "Disordered" regions as a wavy backbone.</li>' +
      '<li><b>Domain context</b> uses InterPro domains and repeats: ' + Object.keys(DOMAIN_LABEL).map((k) => DOMAIN_LABEL[k] + ' = ' + lc1(DOMAIN_HELP[k])).join('; ') + '. Nearest-domain classes come from the lab domain-class ontology (config/domain_class_ontology.json).</li>' +
      '<li><b>RG/RGG-rich segments</b> in the figure and domain table: ' + esc(RG_HELP.charAt(0).toLowerCase() + RG_HELP.slice(1)) + '.</li>' +
      '<li><b>Functional classes</b> come from the lab functional ontology (config/functional_ontology.json). Subclasses are shown under the broad class that most of their proteins share.</li>' +
      '<li><b>Clustering</b>: nearest methyl-R spacing and maximum sites in 3–50 aa windows per protein; dense windows come from a multiscale rolling-window scan against a same-protein arginine null (1,000 permutations, q ≤ 0.05). Architecture classes describe proteins with ≥3 sites in a 20 aa span.</li>' +
      '<li><b>Arg enrichment</b> is the per-protein shrinkage log2 enrichment of methylated arginines relative to the proteome-wide rate, with a Fisher test q value.</li>' +
      '<li><b>Condensates</b> are CD-CODE memberships; condensate-level enrichment uses condensates with at least five members.</li>' +
      '<li><b>Disease</b> associations list DisGeNET terms significantly enriched among methylarginine proteins (clusterProfiler) that include the gene or one of its synonyms, grouped into the audited disease families. The ALS/MN/FTD option uses the full audited family gene set.</li>' +
      '</ul>' +
      ''],
      ['protein-pages', 'Live annotations on protein pages', '' +
        '<p><b>Protein map rows.</b> The Rows menu switches each row on or off and moves it up or down; the choice and order are remembered in this browser and used by the SVG and PNG exports. ' +
        '<b>Charge</b> is the mean residue charge of a sliding 5-residue window on the EMBOSS <i>charge</i> scale (K and R +1, H +0.5, D and E −1): the value at residue <i>i</i> is the mean of residues <i>i</i>−2 to <i>i</i>+2, and within two residues of either terminus the window is cut short and the mean is taken over the residues actually present (3 or 4). Blue above the line is basic, red below acidic. ' +
        '<b>Lysine methylation</b>, phosphorylation, lysine acylation and ubiquitin/SUMO sites come from the pipeline\'s cross-PTM master table (iPTMnet sources, UniProt and IEDB) and are not exhaustive: it holds 234 lysine-methylation sites on 130 proteins. Composition (UniProt compositional bias) is off by default.</p>' +
        '<p><b>Regions</b> are UniProt Region and Motif features (motifs outlined). Overlapping features, which UniProt often nests (NPM1 1–186 and 1–117), are stacked on separate lines so none hides another; a label sits inside its box, after it, or shortened, and the tooltip always gives the full text and coordinates.</p>' +
        '<p><b>Journal panels.</b> The contacts cartoon and the Figures-tab protein comparison can be exported at a printed column width (85, 114 or 174 mm, the usual Cell Press columns): the layout is narrowed so the smallest labels print at 6 pt, the SVG carries its width in millimetres, and the PNG is rendered at 300, 600 or 1000 dpi with the resolution recorded in the file.</p>' +
        '<p><b>Contacts at a glance</b> draws the protein once more as a simple linear cartoon: the methylarginines coloured by interaction status, and under them, on the same residue scale, RNA and DNA crosslinks in living cells, PDB interface residues (PDBe-KB, live), and predicted-complex interfaces (assessed at methylarginines only). A dashed line drops from each methylarginine with evidence at the residue itself and a ring marks which kind; diamonds on the backbone are lysine methylation. The Interactions tab shows the same cartoon for any protein, and the Figures tab adds the crosslink and PDB lanes under every protein of the comparison, which can be ordered by contact evidence, length or name.</p>' +
      '<ul>' +
      '<li><b>UniProt</b>: domains, regions, disordered segments and compositional bias for the map; natural variants and mutagenesis (variants at a methylarginine are highlighted, with their disease notes); and UniProt\'s own modified-residue term for each site (asymmetric, symmetric or monomethylarginine). Methylarginines that UniProt annotates but the union lacks are listed under the site table.</li>' +
      '<li><b>AlphaFold DB</b>: the predicted model coloured by pLDDT (very high > 90, confident 70–90, low 50–70, very low < 50), the pLDDT at each site, and the mean AlphaMissense pathogenicity of the 19 R→X substitutions at each site (likely benign < 0.34, likely pathogenic ≥ 0.564). Low pLDDT usually marks disordered segments, where most methylarginines sit.</li>' +
      '<li><b>Experimental structures</b>: PDBe SIFTS mappings list the PDB entries for the protein; observed-residue ranges decide which methylarginines are actually resolved in the mapped chain. An RCSB search flags entries whose chain for this protein contains a methylated arginine residue (chemical components 2MR, DA2, NMM or AGM); those residues are drawn in orange when the entry is viewed. NMR entries show the first model of the ensemble.</li>' +
      '</ul>' +
      ''],
      ['interactions', 'Interactions', '' +
      '<ul>' +
      '<li><b>Protein partners</b>: STRING v12.5 physical subnetwork (CC BY 4.0). High-confidence links (score ≥ 0.7) are part of the exported data and drive the filters and overview; protein pages fetch medium-confidence links (≥ 0.4) and the links among a protein\'s partners live from the STRING API. Partners are coloured by their own methylation in this set. PRMTs and proteins with a Tudor domain (InterPro IPR002999), WD40 repeats (IPR001680) or a PHD finger (IPR001965, IPR019787) are marked, from UniProt InterPro cross-references of reviewed human entries; "any other protein" in the filters means a partner in none of these classes.</li>' +
      '<li><b>Complexes</b>: Complex Portal curated complexes and hu.MAP 3.0 predicted complexes (CC0). Enrichment of methylated members is a one-sided hypergeometric test against the reviewed human proteome, BH-adjusted separately for curated and predicted complexes.</li>' +
      '<li><b>RNA and DNA binding</b>: Gene Ontology GO:0003723 and GO:0003677 including child terms (QuickGO/GOA, CC BY 4.0); experimental evidence means EXP, IDA, IPI, IMP, IGI, IEP or their high-throughput equivalents. Protein pages add live ENCODE eCLIP and TF ChIP-seq experiments and ChIP-Atlas target genes (CC BY 4.0).</li>' +
      '<li><b>Interfaces</b>: PDBe-KB lists every residue that contacts another protein, RNA or DNA in at least one PDB structure, in UniProt numbering. Sites carry these flags. For the two most recent entries per site and partner type, the residues within 4 Å of the arginine side chain (NE, CZ, NH1, NH2, CD) were read from the structure itself (RCSB ModelServer, SIFTS SEQRES mapping), and partner residues are given in UniProt numbering. The overview compares methylated arginines with the other structurally observed arginines of the same protein (Cochran–Mantel–Haenszel, each protein a stratum), also restricted to ordered regions (outside UniProt/MobiDB-lite disordered regions) and to residues inside UniProt domains. In the 3D view, nucleic acids are coloured and methylarginine side chains within 4 Å of RNA, DNA or another chain are listed.</li>' +
      '<li><b>RNA and DNA contact in living cells</b>: UV crosslinking mass spectrometry locates the residue that touched RNA or DNA. RNA-crosslinked residues come from the pepR-MS supplement (Sha, Kuster & Trendel 2025), which re-analyses Kramer 2014, XRNAX 2019, RBS-ID 2020 and pRBS-ID 2021 alongside pepR-MS with one pipeline; RBDmap adds RNA-bound fragments (Castello et al. 2016). DNA-crosslinked peptides with localised residues and the DNA interactome come from Trendel et al. 2025 (MCF7); UniProt adds curated DNA-binding regions. Every position was checked against the canonical sequence (26,012 matched, 68 dropped). Sites are flagged at a crosslinked residue, within 3 residues, or inside a fragment, peptide or region, and the overview compares them with the other arginines of the same proteins, also without RG/GR dipeptides.</li>' +
      '<li><b>Writer and reader paths</b>: the protein network can add a second ring of PRMTs and Tudor-domain proteins that bind its partners (STRING ≥ 0.7).</li>' +
      '<li><b>Residue numbering</b>: positions follow the canonical UniProt sequence with the initiator methionine as 1. Papers often number differently: histone marks omit the methionine (H3R2 is UniProt R3 of P68431, H4R3 is R4 of P62805), and some receptors are numbered from the mature protein (EGFR R198/R200 in Liao et al. 2015 are UniProt R222/R224).</li>' +
      '<li><b>Unreviewed UniProt entries</b>: about 870 sites in the union sit on unreviewed (TrEMBL) entries, mostly without a staged sequence, so they lack motif, disorder and domain context. About 580 of them repeat a site at the same position on the reviewed entry of the same gene, so the union counts those residues twice. The explorer keeps the union as published, tags these entries and sites, and offers a "Reviewed UniProt entries only" filter (Provenance); the sites TSV has <code>uniprot_reviewed</code> and <code>same_site_on_reviewed_entry</code> columns.</li>' +
      '<li><b>Statistics</b> (interaction export, <code>scripts/37_export_interaction_data.py</code>): SciPy for Fisher\'s exact test, the hypergeometric test (complexes), Benjamini–Hochberg adjustment and the Mann–Whitney U test; statsmodels for the Mantel–Haenszel pooled odds ratio within proteins, its Robins–Breslow–Greenland 95% confidence interval and the continuity-corrected Cochran–Mantel–Haenszel test. Pooled odds ratios across proteins carry Woolf (log) 95% intervals; log2 odds ratios in the GO and interface tables add 0.5 to each cell (Haldane–Anscombe).</li>' +
      '<li><b>Interactions tab</b>: the key findings are written from the exported statistics, so their wording (more often, less often, not significantly different) follows the numbers. "Self-contact" marks residues that touch another copy of the same protein (a homodimer or homo-oligomer) in a PDB entry.</li>' +
      '<li><b>Interface geometry</b>: for every methylarginine at a PDB interface and up to 12 non-methylated interface arginines per protein, the 10 Å environment in one entry (the most recent showing the interface) gives relative solvent accessibility alone and in the complex (Shrake–Rupley, Arg maximum 274 Å² from Tien et al. 2013) and the Levy (' + doiLink('10.1016/j.jmb.2010.09.028', '2010') + ') class: core, rim or support. Guanidino nitrogens within 3.5 Å of a nucleobase N/O count as base hydrogen bonds, of a phosphate or sugar O as backbone hydrogen bonds. Both are compared within proteins (Mantel–Haenszel). The RG/RGG table repeats the interface test within each sequence context.</li>' +
      '<li><b>Methylarginines modelled in the PDB</b>: entries containing 2MR, DA2, NMM or AGM (RCSB search), their residues from PDBe, mapped to UniProt through SIFTS; union sites are flagged and positions the union lacks are listed on protein pages.</li>' +
      '<li><b>References</b>: CD-CODE publications for each membership are listed in full (titles from NCBI E-utilities) with a PubMed search for all of them.</li>' +
      '<li><b>Condensates</b>: protein pages show CD-CODE\'s own confidence score (1–5), role and evidence types for that protein in each condensate (CD-CODE API, release ' + esc((IX && IX.meta && IX.meta.cdcode_version) || '2.3') + '). Condensate-wide statistics (methylated share and enrichment) are shown separately and labelled as being about the condensate.</li>' +
      '<li><b>Disease</b>: protein pages show UniProt\'s curated disease involvement and Open Targets Platform association scores for the gene (live; Open Targets data is CC0). DisGeNET terms come from the dataset-level enrichment and are listed separately as context.</li>' +
      '</ul>' +
      ''],
      ['licences', 'Data sources, versions and licences', '' +
      '<p>Every external source behind the explorer, with the release used and its licence. Please cite the sources you rely on.</p>' +
      DATA_SOURCES.map(([group, rows]) => '<h4>' + esc(group) + '</h4><div class="tbl-wrap"><table class="tbl compact src-tbl"><thead><tr><th>Source</th><th>Used for</th><th>Version</th><th>Licence</th><th>Citation</th></tr></thead><tbody>' +
        rows.map(([name, use, ver, lic, cite, doi]) => '<tr><td><b>' + esc(name) + '</b></td><td>' + esc(use) + '</td><td>' + esc(ver || '—') + '</td><td>' + esc(lic || '—') + '</td><td>' + (cite ? esc(cite) + (doi ? ' · ' + doiLink(doi) : '') : '—') + '</td></tr>').join('') + '</tbody></table></div>').join('') +
      '<p><b>Licence of this site.</b> The explorer\'s data, figures and text are licensed CC BY-NC-SA 4.0 and its code MIT (see LICENSE.md in the repository); share-alike, non-commercial terms are what PhosphoSitePlus, DisGeNET and pepR-MS require of the combined data. Third-party data keep the licences listed above.</p>' +
      '<p class="note">Interaction data exported ' + esc((IX && IX.meta && IX.meta.generated) || '—') + '; site data exported ' + esc(M.generated) + '. "At staging" and "at export" mean the release current on those dates. Software licences cover the code, not the data.</p>' +
      ''],
      ['guardrails', 'Interpretive guardrails', '' +
      '<ul>' +
      '<li>The integrated union is the primary analysis universe; source support is provenance, not a confidence tier.</li>' +
      '<li>State labels are useful where a source supports them, but the union is not a clean Rme1 / Rme2a / Rme2s state-resolved dataset.</li>' +
      '<li>Domain adjacency is descriptive: an IDR-length-conditioned null shows that simple domain-edge proximity is largely explained by IDR geometry and length.</li>' +
      '<li>The stronger architectural signal is short-range methylarginine clustering relative to same-protein arginine positions.</li>' +
      '<li>Motif and clustering annotations need sequence-confirmed canonical residues; ' + fmt(M.window_flags.missing) + ' sites lack a usable sequence window and ' + fmt(M.window_flags.center_not_r) + ' have a non-R canonical residue (flagged with !).</li>' +
      '<li>Coverage follows the source compilations. Methylation reported only in low-throughput studies that none of the sources compiled is not in the union.</li>' +
      '</ul>' +
      ''],
      ['downloads', 'Downloads', '<p class="note">' + esc(PSP_NOTICE) + ' Site downloads that include PhosphoSitePlus® sites carry these terms; DisGeNET-derived disease terms are CC BY-NC-SA 4.0.</p>' +
        '<p>Use <b>Download TSV</b> on the Proteins or Sites view to export the current selection. The data files are regenerated from pipeline outputs with <code>python scripts/36_export_site_explorer_data.py</code>; see the repository README.</p>' +
        '<p class="note">Data exported ' + esc(M.generated) + ' · ' + fmt(M.n_sites) + ' sites · ' + fmt(M.n_proteins) + ' proteins.</p>'],
    ];
    el.innerHTML = '<div class="docs"><nav class="docs-toc card" aria-label="Documentation contents"><p class="flabel">Contents</p><ol>' +
      parts.map(([id, title]) => '<li><button class="linkish" type="button" data-doc="' + id + '">' + esc(title) + '</button></li>').join('') + '</ol></nav>' +
      '<div class="card prose docs-body">' + parts.map(([id, title, html]) => sec(id, title, html)).join('') + '</div></div>';
    bindDocLinks(el);
    if (UI.docJump) setTimeout(jumpDoc, 0);
  }

  // ---------------------------------------------------------- docs links
  // Short "Method" links from charts and figures to their explanation on the Documentation page.
  const lc1 = (t) => t.charAt(0).toLowerCase() + t.slice(1);
  const docLink = (key) => ' <button class="linkish doc-link" type="button" data-doc="' + key + '">Method</button>';
  function bindDocLinks(root) {
    $$('[data-doc]', root).forEach((b) => b.addEventListener('click', () => openDoc(b.dataset.doc)));
  }
  function openDoc(key) {
    UI.docJump = key;
    if (UI.view === 'about') jumpDoc(); else go('about');
  }
  function jumpDoc() {
    const t = UI.docJump && document.getElementById('doc-' + UI.docJump);
    UI.docJump = null;
    if (!t) { window.scrollTo({ top: 0 }); return; }
    t.scrollIntoView({ block: 'start' });
    t.classList.add('flash');
    setTimeout(() => t.classList.remove('flash'), 1600);
  }

  // ---------------------------------------------------------------- figures
  // Figures from the current selection (search and filters): sites aligned at the methylarginine with a sequence logo,
  // and the selected proteins on one residue scale. Exports go through exportFigure, which adds the credit band.
  // Residue colours follow WebLogo's "chemistry" scheme.
  const AA_CHEM = {
    G: 'polar', S: 'polar', T: 'polar', Y: 'polar', C: 'polar', Q: 'polar', N: 'polar', K: 'basic', R: 'basic', H: 'basic', D: 'acidic', E: 'acidic',
    A: 'hydrophobic', V: 'hydrophobic', L: 'hydrophobic', I: 'hydrophobic', P: 'hydrophobic', W: 'hydrophobic', F: 'hydrophobic', M: 'hydrophobic',
  };
  const CHEM_COLOR = { polar: '#1A9850', basic: '#2166AC', acidic: '#C0392B', hydrophobic: '#1F2933' };
  const CHEM_TINT = { polar: '#DDF0E3', basic: '#DCE8F6', acidic: '#F7DEDA', hydrophobic: '#ECEEF0' };
  const CHEM_LABEL = { basic: 'Basic (K R H)', acidic: 'Acidic (D E)', polar: 'Polar (G S T Y C Q N)', hydrophobic: 'Hydrophobic (A V L I P W F M)' };
  const MONO = 'Menlo, Consolas, &quot;DejaVu Sans Mono&quot;, &quot;Liberation Mono&quot;, monospace';
  const FIG_ORDERS = [['protein', 'Protein and position'], ['right', 'Sequence after the site'], ['left', 'Sequence before the site'], ['reports', 'Most reports first']];

  // Selected sites that have a sequence window, in the chosen order.
  function figSites() {
    const list = R.sites.filter((s) => s.win && s.win.length === 31);
    const o = UI.figs.order;
    if (o === 'reports') return list.slice().sort((a, b) => popcount(b.study) - popcount(a.study) || P[a.p].label.localeCompare(P[b.p].label) || a.pos - b.pos);
    // "Protein and position" keeps the order of a searched gene list, otherwise alphabetical.
    const rank = new Map(Q.proteins && Q.proteins.size ? Array.from(Q.proteins).map((i, k) => [i, String(k).padStart(5, '0')]) : []);
    const key = o === 'right' ? (s) => s.win.slice(16) : o === 'left' ? (s) => s.win.slice(0, 15).split('').reverse().join('') : (s) => (rank.get(s.p) || P[s.p].label) + ':' + String(s.pos).padStart(6, '0');
    return list.map((s) => [key(s), s]).sort((a, b) => (a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : 0)).map((x) => x[1]);
  }

  // Per-position residue frequencies over the windows (terminal padding left out) and the information content in bits,
  // with the small-sample correction of Schneider et al. (1986).
  function logoColumns(list, flank) {
    const cols = [];
    for (let k = 15 - flank; k <= 15 + flank; k++) {
      const c = {};
      let n = 0;
      for (const s of list) { const a = s.win[k]; if (AA_CHEM[a]) { c[a] = (c[a] || 0) + 1; n++; } }
      let H = 0;
      for (const a in c) { const f = c[a] / n; H -= f * Math.log2(f); }
      const ic = n ? Math.max(0, Math.log2(20) - H - 19 / (2 * Math.LN2 * n)) : 0;
      cols.push({ rel: k - 15, n, ic, freqs: Object.keys(c).map((a) => [a, c[a] / n]).sort((x, y) => x[1] - y[1]) });
    }
    return cols;
  }

  // One residue cell: chemistry tint and letter, monospace.
  const cellSvg = (a, x, y, cw, rh) => (AA_CHEM[a] ? '<rect x="' + x + '" y="' + y + '" width="' + cw + '" height="' + rh + '" fill="' + CHEM_TINT[AA_CHEM[a]] + '"/>' : '') +
    (a && a !== '_' ? '<text x="' + (x + cw / 2) + '" y="' + (y + rh - 4) + '" font-size="11" font-family="' + MONO + '" text-anchor="middle" fill="' + (AA_CHEM[a] ? CHEM_COLOR[AA_CHEM[a]] : '#5F6F78') + '">' + esc(a) + '</text>' : '');

  // Chemistry legend; two rows when the figure is narrower than one row of entries. Returns { svg, h }.
  function chemLegend(x, y, W) {
    const keys = Object.keys(CHEM_LABEL), widths = [110, 100, 160, 210];
    const perRow = widths.reduce((a, b) => a + b, 0) <= W - 2 * x ? 4 : 2;
    let svg = '', cx = x, cy = y;
    keys.forEach((k, i) => {
      if (i && i % perRow === 0) { cx = x; cy += 16; }
      svg += '<rect x="' + cx + '" y="' + (cy - 9) + '" width="11" height="11" fill="' + CHEM_TINT[k] + '" stroke="' + CHEM_COLOR[k] + '"/><text x="' + (cx + 16) + '" y="' + cy + '" font-size="10" fill="#5F6F78">' + esc(CHEM_LABEL[k]) + '</text>';
      cx += perRow === 4 ? widths[i] : 230;
    });
    return { svg, h: cy - y + 16 };
  }

  function alignSvg(all, o) {
    const pal = EXPORT_PAL, flank = o.flank, ncol = 2 * flank + 1;
    const rows = all.slice(0, o.rows);
    const cw = 15, rh = 15, labelW = 150, x0 = 16 + labelW;
    const W = Math.max(640, x0 + ncol * cw + 40);
    const colX = (k) => x0 + k * cw;
    let g = '<text x="16" y="28" font-size="17" font-weight="700" fill="' + pal.ink + '">' + esc(o.title) + '</text>' +
      '<text x="16" y="47" font-size="11" fill="' + pal.muted + '">' + esc(o.sub) + '</text>';
    let y = 66;
    if (o.logo !== 'none' && all.length) {
      const cols = logoColumns(all, flank);
      const LH = 100, bits = o.logo === 'bits', maxV = bits ? Math.log2(20) : 1;
      const ticks = bits ? [0, 1, 2, 3, 4] : [0, 0.5, 1];
      g += '<line x1="' + (x0 - 6) + '" x2="' + (x0 - 6) + '" y1="' + y + '" y2="' + (y + LH) + '" stroke="' + pal.backbone + '"/>';
      for (const t of ticks) {
        const ty = y + LH - (t / maxV) * LH;
        g += '<line x1="' + (x0 - 9) + '" x2="' + (x0 - 6) + '" y1="' + r1(ty) + '" y2="' + r1(ty) + '" stroke="' + pal.backbone + '"/><text x="' + (x0 - 12) + '" y="' + r1(ty + 3) + '" font-size="9.5" text-anchor="end" fill="' + pal.muted + '">' + t + '</text>';
      }
      g += '<text transform="translate(' + (x0 - 34) + ' ' + (y + LH / 2) + ') rotate(-90)" font-size="10" text-anchor="middle" fill="' + pal.muted + '">' + (bits ? 'bits' : 'frequency') + '</text>';
      cols.forEach((c, k) => {
        let base = y + LH;
        const total = bits ? c.ic : 1;
        for (const [a, f] of c.freqs) {
          const h = (f * total / maxV) * LH;
          if (h < 0.6) continue;
          // Arial capitals stand 71.6% of the font size above the baseline: scale each letter to its share of the stack.
          g += '<text x="' + colX(k) + '" y="0" font-size="100" font-weight="700" font-family="Arial, Helvetica, sans-serif" fill="' + CHEM_COLOR[AA_CHEM[a]] + '" textLength="' + (cw - 1) +
            '" lengthAdjust="spacingAndGlyphs" transform="translate(0 ' + r1(base) + ') scale(1 ' + (h / 71.6).toFixed(4) + ')">' + a + '</text>';
          base -= h;
        }
      });
      y += LH + 6;
    }
    for (let k = 0; k < ncol; k++) {
      const rel = k - flank;
      if (rel === 0 || rel % (flank > 10 ? 5 : flank > 5 ? 2 : 1) === 0) {
        g += '<text x="' + (colX(k) + cw / 2) + '" y="' + (y + 10) + '" font-size="9.5" text-anchor="middle" fill="' + (rel ? pal.muted : pal.ink) + '"' + (rel ? '' : ' font-weight="700"') + '>' + (rel > 0 ? '+' + rel : rel === 0 ? '0' : '−' + -rel) + '</text>';
      }
    }
    y += 16;
    const top = y;
    rows.forEach((s, r) => {
      const yy = top + r * rh;
      g += '<text x="16" y="' + (yy + rh - 4) + '" font-size="11" fill="' + pal.ink + '"><tspan font-weight="700">' + esc(P[s.p].label) + '</tspan> ' + rl(s.pos, P[s.p]) + '</text>';
      for (let k = 0; k < ncol; k++) g += cellSvg(s.win[15 - flank + k], colX(k), yy, cw, rh);
    });
    const bodyH = rows.length * rh;
    if (rows.length) g += '<rect x="' + colX(flank) + '" y="' + (top - 1) + '" width="' + cw + '" height="' + (bodyH + 2) + '" fill="none" stroke="' + pal.ink + '" stroke-width="1.2"/>';
    y = top + bodyH + 18;
    if (rows.length && all.length > rows.length) { g += '<text x="16" y="' + y + '" font-size="10.5" fill="' + pal.muted + '">Showing ' + fmt(rows.length) + ' of ' + fmt(all.length) + ' sites; the logo uses all ' + fmt(all.length) + '.</text>'; y += 18; }
    const leg = chemLegend(16, y + 4, W);
    g += leg.svg;
    const H = y + leg.h + 6;
    return { W, H, svg: '<svg xmlns="http://www.w3.org/2000/svg" width="' + W + '" height="' + H + '" viewBox="0 0 ' + W + ' ' + H + '" font-family=\'' + pal.font.replace(/'/g, '"') + '\'><rect width="' + W + '" height="' + H + '" fill="' + pal.bg + '"/>' + g + '</svg>' };
  }

  // Proteins for the comparison figure: a searched gene list keeps its order, otherwise the proteins with most sites.
  function figProteins() {
    const inSel = new Set(R.prots.map((p) => p.i));
    const ordered = Q.proteins && Q.proteins.size ? Array.from(Q.proteins).filter((i) => inSel.has(i)).map((i) => P[i]) : R.prots.slice().sort((a, b) => R.counts[b.i] - R.counts[a.i] || a.label.localeCompare(b.label));
    // An unreviewed entry whose reviewed entry is also selected repeats its sites in other numbering; leave it out.
    const ids = new Set(ordered.map((p) => p.i));
    return ordered.filter((p) => !(p.unrev && p.reviewedTwin && ids.has(p.reviewedTwin.i))).slice(0, UI.figs.nProt);
  }

  // Proteins on one residue scale. o.W narrows the canvas and o.panel switches to the journal-panel layout (no title,
  // a slim label column, short lane names, wrapped legend) so text keeps its printed size at a column width.
  function compareSvg(items, o) {
    const pal = EXPORT_PAL, W = o.W || 1000, panel = !!o.panel;
    const m = panel ? 4 : 16, labelW = panel ? 60 : 170, x0 = m + labelW, x1 = W - (panel ? 6 : 24);
    const fsGene = panel ? 10.5 : 12.5, fsSmall = panel ? 8.5 : 10, fsLane = 8.5;
    const mode = COLOR_MODES[o.color] || COLOR_MODES.state, colors = mode.colors(pal);
    const lenOf = ({ p, u }) => (u && u.len) || p.len || Math.ceil(p.maxPos * 1.05) || 100;
    const maxLen = Math.max(...items.map(lenOf));
    let g = '', y = 8;
    if (!panel) {
      g += '<text x="16" y="28" font-size="17" font-weight="700" fill="' + pal.ink + '">' + esc(o.title) + '</text>' +
        '<text x="16" y="47" font-size="11" fill="' + pal.muted + '">' + esc(o.sub) + '</text>';
      y = 62;
    }
    const rowH = o.lanes ? 112 : 78;
    for (const it of items) {
      const { p, u } = it;
      const len = lenOf(it), span = o.scale === 'fit' ? len : maxLen;
      const X = (pos) => x0 + ((pos - 0.5) / span) * (x1 - x0);
      const base = y + 56;
      const shown = p.sites.filter((j) => R.shown.has(j)).map((j) => S[j]);
      const gfs = Math.min(fsGene, (labelW - 4) / (p.label.length * 0.62));
      g += '<text x="' + m + '" y="' + (base - 2) + '" font-size="' + r1(gfs) + '" font-weight="700" fill="' + pal.ink + '">' + esc(p.label) + '</text>' +
        (panel ? '<text x="' + m + '" y="' + (base + 9) + '" font-size="' + fsSmall + '" fill="' + pal.muted + '">' + fmt(len) + ' aa</text>'
          : '<text x="16" y="' + (base + 12) + '" font-size="10" fill="' + pal.muted + '">' + esc(p.acc) + ' · ' + fmt(len) + ' aa · ' + shown.length + (shown.length === p.sites.length ? '' : '/' + p.sites.length) + ' sites</text>');
      const feats = (u && u.feats) || [];
      for (const f of feats.filter((f) => f.type === 'Region' && /disordered/i.test(f.desc))) {
        g += '<rect x="' + r1(X(f.start)) + '" y="' + (base - 9) + '" width="' + r1(Math.max(1, X(f.end) - X(f.start))) + '" height="18" fill="' + pal.disFill + '" opacity="0.55"/>';
      }
      g += '<rect x="' + r1(X(1)) + '" y="' + (base - 3) + '" width="' + r1(X(len) - X(1)) + '" height="6" rx="3" fill="' + pal.backboneFill + '" stroke="' + pal.backboneEdge + '"/>';
      for (const f of feats.filter((f) => ['Domain', 'Repeat', 'Zinc finger'].includes(f.type))) {
        const a = X(f.start), b = X(f.end), wd = Math.max(2, b - a);
        const rep = f.type === 'Repeat';
        g += '<rect x="' + r1(a) + '" y="' + (base - (rep ? 5 : 7)) + '" width="' + r1(wd) + '" height="' + (rep ? 10 : 14) + '" rx="' + (rep ? 2 : 3) + '" fill="' + pal.domainFill + '"' + (rep ? ' opacity="0.6"' : '') + ' stroke="' + pal.domainInk + '" stroke-width="0.6"><title>' + esc(f.desc || f.type) + ' ' + f.start + '–' + f.end + '</title></rect>';
        const lab = rep ? '' : (f.desc || '').replace(/ domain$/i, '');
        if (lab && lab.length * 5.6 < wd - 4) g += '<text x="' + r1(a + wd / 2) + '" y="' + (base + 3.5) + '" font-size="9" text-anchor="middle" fill="' + pal.domainInk + '">' + esc(lab) + '</text>';
      }
      const levels = [];
      for (const s of shown) {
        const hx = X(s.pos);
        let lv = 0;
        while (levels[lv] != null && hx - levels[lv] < 7 && lv < 5) lv++;
        levels[lv] = hx;
        const hy = base - 18 - lv * 8;
        g += '<line x1="' + r1(hx) + '" x2="' + r1(hx) + '" y1="' + (base - 7) + '" y2="' + hy + '" stroke="' + pal.stem + '" stroke-width="1"/>' +
          '<circle cx="' + r1(hx) + '" cy="' + hy + '" r="3.2" fill="' + colors[mode.cat(s)] + '" stroke="#FFFFFF" stroke-width="0.6"><title>' + esc(p.label + ' ' + rl(s.pos, P[s.p])) + '</title></circle>';
      }
      // Contact lanes under the protein: crosslinks in cells and PDB interface residues, ringed where a methylarginine
      // itself carries the evidence.
      if (o.lanes) {
        const ln = contactLanes(p, p.iso ? null : IFDATA.get(p.acc.split('-')[0]) || null);
        const pdb = { ranges: (ln[2].ranges || []).concat(ln[3].ranges || []), hit: (s) => ln[2].hit(s) || ln[3].hit(s), color: IFACE_COLOR[1] };
        [[panel ? 'RNA xl' : 'RNA crosslink', ln[0]], [panel ? 'DNA xl' : 'DNA crosslink', ln[1]], [panel ? 'PDB' : 'PDB interface', pdb]].forEach(([name, lane], k) => {
          const cy = base + 19 + k * 11;
          g += '<line x1="' + r1(X(1)) + '" x2="' + r1(X(len + 1)) + '" y1="' + cy + '" y2="' + cy + '" stroke="' + pal.grid + '"/>' +
            '<text x="' + (x0 - 4) + '" y="' + (cy + 3) + '" font-size="' + fsLane + '" text-anchor="end" fill="' + pal.muted + '">' + name + '</text>';
          for (const [a, b, op] of lane.boxes || []) g += '<rect x="' + r1(X(a)) + '" y="' + (cy - 3) + '" width="' + r1(Math.max(1.4, X(b + 1) - X(a))) + '" height="6" fill="' + lane.color + '" fill-opacity="' + op + '"/>';
          for (const [a, b, , color] of lane.ranges || []) g += '<rect x="' + r1(X(a)) + '" y="' + (cy - 3.5) + '" width="' + r1(Math.max(1.4, X(b + 1) - X(a))) + '" height="7" fill="' + color + '" fill-opacity="0.85"/>';
          for (const t of lane.ticks || []) g += '<rect x="' + r1(X(t) - 0.7) + '" y="' + (cy - 4) + '" width="' + r1(Math.max(1.4, X(t + 1) - X(t))) + '" height="8" fill="' + lane.color + '"/>';
          for (const s of shown.filter(lane.hit)) g += '<circle cx="' + r1(X(s.pos)) + '" cy="' + cy + '" r="3.1" fill="' + pal.bg + '" stroke="' + (lane.color || pal.ink) + '" stroke-width="1.4"><title>' + esc(p.label + ' ' + rl(s.pos, P[s.p]) + ': ' + name + ' at this methylarginine') + '</title></circle>';
        });
      }
      y += rowH;
    }
    if (o.scale !== 'fit') {
      const step = niceStep(maxLen) * (panel && W < 500 ? 2 : 1);
      g += '<line x1="' + x0 + '" x2="' + x1 + '" y1="' + (y + 4) + '" y2="' + (y + 4) + '" stroke="' + pal.backbone + '"/>';
      for (let t = 0; t <= maxLen; t += step) {
        const tx = x0 + (t / maxLen) * (x1 - x0);
        g += '<line x1="' + r1(tx) + '" x2="' + r1(tx) + '" y1="' + (y + 4) + '" y2="' + (y + 8) + '" stroke="' + pal.backbone + '"/><text x="' + r1(tx) + '" y="' + (y + 19) + '" font-size="9.5" text-anchor="middle" fill="' + pal.muted + '">' + fmt(t) + '</text>';
      }
      g += '<text x="' + ((x0 + x1) / 2) + '" y="' + (y + 33) + '" font-size="' + fsSmall + '" text-anchor="middle" fill="' + pal.muted + '">residue</text>';
      y += 40;
    }
    // Legend, wrapped to the canvas width.
    const cw = panel ? 4.7 : 5.5;
    let lx = m;
    y += 14;
    const leg = (draw, label) => {
      const w = 16 + label.length * cw;
      if (lx > m && lx + w > W - m) { lx = m; y += 15; }
      g += draw(lx, y - 4) + '<text x="' + (lx + 13) + '" y="' + y + '" font-size="' + fsSmall + '" fill="' + pal.muted + '">' + esc(label) + '</text>';
      lx += w + 8;
    };
    mode.cats.forEach((c, k) => leg((a, b) => '<circle cx="' + (a + 4) + '" cy="' + b + '" r="4" fill="' + colors[k] + '"/>', c));
    leg((a, b) => '<rect x="' + a + '" y="' + (b - 4) + '" width="10" height="8" rx="2" fill="' + pal.domainFill + '" stroke="' + pal.domainInk + '" stroke-width="0.6"/>', 'UniProt domain or repeat');
    leg((a, b) => '<rect x="' + a + '" y="' + (b - 4) + '" width="10" height="8" fill="' + pal.disFill + '" opacity="0.55"/>', 'Disordered (UniProt)');
    if (o.lanes) {
      leg((a, b) => '<rect x="' + (a + 4) + '" y="' + (b - 4) + '" width="2" height="8" fill="' + IFACE_COLOR[4] + '"/>', 'RNA-crosslinked residue');
      leg((a, b) => '<rect x="' + (a + 4) + '" y="' + (b - 4) + '" width="2" height="8" fill="' + IFACE_COLOR[8] + '"/>', 'DNA-crosslinked residue');
      leg((a, b) => '<rect x="' + a + '" y="' + (b - 3.5) + '" width="10" height="7" fill="' + IFACE_COLOR[1] + '"/>', 'PDB protein interface');
      leg((a, b) => '<rect x="' + a + '" y="' + (b - 3.5) + '" width="10" height="7" fill="' + IFACE_COLOR[4] + '"/>', 'PDB RNA interface');
      leg((a, b) => '<rect x="' + a + '" y="' + (b - 3.5) + '" width="10" height="7" fill="' + IFACE_COLOR[8] + '"/>', 'PDB DNA interface');
      leg((a, b) => '<rect x="' + a + '" y="' + (b - 3.5) + '" width="10" height="7" fill="' + IFACE_COLOR[2] + '"/>', 'PDB self-contact');
      leg((a, b) => '<circle cx="' + (a + 5) + '" cy="' + b + '" r="3.1" fill="' + pal.bg + '" stroke="' + pal.ink + '" stroke-width="1.4"/>', 'Evidence at the methylarginine');
    }
    const H = y + 10;
    return { W, H, svg: '<svg xmlns="http://www.w3.org/2000/svg" width="' + W + '" height="' + H + '" viewBox="0 0 ' + W + ' ' + H + '" font-family=\'' + pal.font.replace(/'/g, '"') + '\'><rect width="' + W + '" height="' + H + '" fill="' + pal.bg + '"/>' + g + '</svg>' };
  }

  function niceStep(span) {
    const raw = span / 8, mag = Math.pow(10, Math.floor(Math.log10(raw)));
    return [1, 2, 5, 10].map((m) => m * mag).find((s) => s >= raw) || 10 * mag;
  }

  // Short description of the current selection, for figure subtitles.
  function selectionText() {
    const c = constraints().map((x) => x.label);
    return c.length ? c.join('; ') : 'all sites';
  }

  function renderFigures(el) {
    const f = UI.figs;
    const sel = figSites();
    const prots = figProteins();
    const sel2 = (id, opts, val) => '<select id="' + id + '">' + opts.map(([v, l]) => '<option value="' + v + '"' + (String(v) === String(val) ? ' selected' : '') + '>' + esc(l) + '</option>').join('') + '</select>';
    const scoped = constraints().length > 0;
    el.innerHTML = resultsHead('Figures', 'From the current selection: ' + fmt(R.sites.length) + ' sites on ' + fmt(R.prots.length) + ' proteins. Exports credit this site and ' + esc(PAPER.short) + '.') +
      (scoped ? '' : '<div class="filter-note">No search or filter is set, so these figures use the whole dataset. Search for proteins or set filters to focus them.</div>') +
      '<div class="card" id="fg-align"><div class="card-head"><div><h3>Site alignment and sequence logo</h3><p class="sub">Selected sites aligned at the methylarginine; the logo uses all ' + fmt(sel.length) + ' with a sequence window.' + docLink('fig-align') + '</p></div>' +
      '<div class="exp" role="group" aria-label="Alignment export"><button type="button" data-fx="align-svg">SVG</button><button type="button" data-fx="align-png">PNG</button><button type="button" data-fx="align-fasta">FASTA</button></div></div>' +
      '<div class="fig-controls"><label class="fig-ctl">Flank ' + sel2('fg-flank', [[5, '±5'], [7, '±7'], [10, '±10'], [15, '±15']], f.flank) + '</label>' +
      '<label class="fig-ctl">Rows ' + sel2('fg-rows', [[25, '25'], [50, '50'], [100, '100'], [250, '250']], f.rows) + '</label>' +
      '<label class="fig-ctl">Order ' + sel2('fg-order', FIG_ORDERS, f.order) + '</label>' +
      '<label class="fig-ctl">Logo ' + sel2('fg-logo', [['freq', 'Frequency'], ['bits', 'Information (bits)'], ['none', 'None']], f.logo) + '</label></div>' +
      '<div class="fig-preview" id="fg-align-box"></div></div>' +
      '<div class="card" id="fg-cmp"><div class="card-head"><div><h3>Protein comparison</h3><p class="sub">' + (Q.proteins && Q.proteins.size ? 'Proteins in your search order' : 'The proteins with most selected sites') + ', with UniProt domains and disordered regions.' + docLink('fig-compare') + '</p></div>' +
      '<div class="exp" role="group" aria-label="Comparison export"><button type="button" data-fx="cmp-svg">SVG</button><button type="button" data-fx="cmp-png">PNG</button></div></div>' +
      '<div class="fig-controls"><label class="fig-ctl">Proteins ' + sel2('fg-nprot', [[4, '4'], [8, '8'], [12, '12'], [20, '20']], f.nProt) + '</label>' +
      '<label class="fig-ctl">Scale ' + sel2('fg-scale', [['common', 'Same residue scale'], ['fit', 'Fit each protein']], f.scale) + '</label>' +
      '<label class="fig-ctl">Colour sites by ' + sel2('fg-color', Object.keys(COLOR_MODES).map((k) => [k, COLOR_MODES[k].name]), f.color) + '</label>' +
      '<label class="fig-ctl">Order ' + sel2('fg-psort', [['default', Q.proteins && Q.proteins.size ? 'Search order' : 'Most selected sites'], ['evidence', 'Most sites with contact evidence'], ['length', 'Longest first'], ['name', 'Name']], f.sortP) + '</label>' +
      '<label class="check fig-ctl"><input type="checkbox" id="fg-lanes"' + (f.lanes ? ' checked' : '') + '> Contact lanes (RNA/DNA crosslinks, PDB interfaces)</label>' + panelControlsHtml() + '</div>' +
      '<div class="fig-preview" id="fg-cmp-box"><p class="note">Loading UniProt features for ' + fmt(prots.length) + ' proteins…</p></div></div>' +
      '<div class="card"><div class="card-head"><div><h3>More figures</h3></div><div class="exp"><button type="button" id="fg-report" title="Profile, evidence map, logo and protein comparison of this selection as one PDF">Insight report (PDF)</button></div></div><ul class="fig-more">' +
      '<li><b>Insight report</b>: this selection\'s profile, evidence map, logo and alignment, protein comparison and protein table on one printable page.</li>' +
      '<li><b>Protein pages</b>: protein map, sequence, 3D structure, interaction network and orthologue alignment, each with SVG or PNG export.</li>' +
      '<li><b>Overview and Interactions</b>: every chart exports as SVG, PNG or TSV.</li></ul></div>';
    const alignOpts = () => ({
      flank: f.flank, rows: f.rows, logo: f.logo,
      title: 'Methylarginine sites aligned at the modified arginine',
      sub: fmt(sel.length) + ' sites with a sequence window · selection: ' + selectionText() + ' · order: ' + FIG_ORDERS.find((x) => x[0] === f.order)[1].toLowerCase(),
    });
    $('#fg-align-box', el).innerHTML = sel.length ? alignSvg(sel, alignOpts()).svg : '<p class="note">No selected site has a sequence window.</p>';
    const cmpOpts = () => ({ scale: f.scale, color: f.color, lanes: f.lanes, title: 'Methylarginine sites on ' + fmt(prots.length) + ' protein' + (prots.length === 1 ? '' : 's'), sub: 'Selection: ' + selectionText() + (f.scale === 'fit' ? ' · each protein scaled to its own length' : ' · common residue scale') });
    let items = null;
    // Rows in the chosen order; contact evidence counts sites of the selection with experimental or predicted status.
    const sortedItems = () => {
      const ev = (p) => p.sites.filter((j) => R.shown.has(j) && ixReady() && ixStatus(S[j]) !== 'none').length;
      const len = (it) => (it.u && it.u.len) || it.p.len || 0;
      const by = { evidence: (a, b) => ev(b.p) - ev(a.p), length: (a, b) => len(b) - len(a), name: (a, b) => a.p.label.localeCompare(b.p.label) }[f.sortP];
      return by ? items.slice().sort(by) : items;
    };
    const drawCmp = () => {
      const b = $('#fg-cmp-box');
      if (b && items) b.innerHTML = items.length ? '<div class="panel-prev"' + panelPreviewStyle() + '>' + compareSvg(sortedItems(), Object.assign(cmpOpts(), panelOpts())).svg + '</div>' +
        (panelState.mm ? '<p class="note">Previewed at twice its printed width of ' + panelState.mm + ' mm; the smallest labels print at 6 pt.</p>' : '') : '<p class="note">No protein selected.</p>';
    };
    const ifaces = () => Promise.all(prots.filter((p) => !p.iso).map((p) => loadIface(p.acc).catch(() => null)));
    Promise.all([Promise.all(prots.map((p) => loadUniProt(p.acc))), f.lanes ? ifaces() : null]).then(([us]) => {
      if (UI.view !== 'figures') return;
      items = prots.map((p, k) => ({ p, u: us[k] }));
      drawCmp();
    });
    const rerender = () => renderFigures(el);
    $('#fg-flank', el).addEventListener('change', (e) => { f.flank = +e.target.value; rerender(); });
    $('#fg-rows', el).addEventListener('change', (e) => { f.rows = +e.target.value; rerender(); });
    $('#fg-order', el).addEventListener('change', (e) => { f.order = e.target.value; rerender(); });
    $('#fg-logo', el).addEventListener('change', (e) => { f.logo = e.target.value; rerender(); });
    $('#fg-nprot', el).addEventListener('change', (e) => { f.nProt = +e.target.value; rerender(); });
    $('#fg-scale', el).addEventListener('change', (e) => { f.scale = e.target.value; drawCmp(); });
    $('#fg-color', el).addEventListener('change', (e) => { f.color = e.target.value; drawCmp(); });
    $('#fg-psort', el).addEventListener('change', (e) => { f.sortP = e.target.value; drawCmp(); });
    bindPanelControls($('#fg-cmp', el), drawCmp);
    $('#fg-lanes', el).addEventListener('change', (e) => { f.lanes = e.target.checked; drawCmp(); if (f.lanes) ifaces().then(() => { if (UI.view === 'figures') drawCmp(); }); });
    $('#fg-report', el).addEventListener('click', () => openReport());
    $$('[data-fx]', el).forEach((b) => b.addEventListener('click', () => {
      const [what, kind] = b.dataset.fx.split('-');
      if (what === 'align') {
        if (!sel.length) return;
        if (kind === 'fasta') {
          download('methylarginine_site_alignment.fasta', sel.map((s) => '>' + P[s.p].label + '_' + rl(s.pos, P[s.p]) + ' ' + P[s.p].acc + '\n' + s.win.slice(15 - f.flank, 16 + f.flank).replace(/_/g, '-')).join('\n') + '\n');
        } else exportFigure(kind, alignSvg(sel, alignOpts()), 'methylarginine_site_alignment');
      } else if (items && items.length) exportPanel(kind, compareSvg(sortedItems(), Object.assign(cmpOpts(), panelOpts())), 'methylarginine_protein_comparison');
    }));
    bindDocLinks(el);
  }

  // ------------------------------------------------- orthologue conservation
  // Ensembl Compara orthologues (REST, live): for each methylarginine, the residue aligned to it in each species.
  // Pairwise alignments are anchored on the human protein, so insertions in other species are not shown.
  const ORTHO_SPECIES = [['mus_musculus', 'Mouse'], ['gallus_gallus', 'Chicken'], ['xenopus_tropicalis', 'Frog'], ['danio_rerio', 'Zebrafish'],
    ['drosophila_melanogaster', 'Fly'], ['caenorhabditis_elegans', 'Worm'], ['saccharomyces_cerevisiae', 'Yeast']];
  const orthoCache = new Map();
  function loadOrthologs(ensg) {
    if (orthoCache.has(ensg)) return orthoCache.get(ensg);
    const url = 'https://rest.ensembl.org/homology/id/human/' + encodeURIComponent(ensg) + '?type=orthologues;aligned=1;sequence=protein;cigar_line=0;' +
      ORTHO_SPECIES.map(([sp]) => 'target_species=' + sp).join(';') + ';content-type=application/json';
    const pr = fetch(url).then((r) => { if (!r.ok) throw new Error('Ensembl returned ' + r.status); return r.json(); }).then((d) => {
      const best = new Map();
      for (const h of ((d.data || [])[0] || {}).homologies || []) {
        const sp = h.target && h.target.species;
        if (sp && (!best.has(sp) || h.target.perc_id > best.get(sp).target.perc_id)) best.set(sp, h);
      }
      return ORTHO_SPECIES.filter(([sp]) => best.has(sp)).map(([sp, name]) => ({ sp, name, h: best.get(sp) }));
    });
    pr.catch(() => orthoCache.delete(ensg));
    orthoCache.set(ensg, pr);
    return pr;
  }

  // Map UniProt positions onto the Ensembl human protein (identical, contained, or by a unique 15-residue anchor), then
  // read the aligned residue of each orthologue at every column of a human-anchored window.
  function orthoModel(p, u, orthos) {
    const uni = u.seq;
    const rows = orthos.map((o) => {
      const src = o.h.source.align_seq, tgt = o.h.target.align_seq;
      const ens = src.replace(/-/g, '');
      const colOf = [];
      for (let c = 0, e = 0; c < src.length; c++) if (src[c] !== '-') colOf[++e] = c;
      const off = ens === uni ? 0 : ens.includes(uni) ? ens.indexOf(uni) : uni.includes(ens) ? -uni.indexOf(ens) : null;
      const toEns = (pos) => {
        if (off != null) { const e = pos + off; return e >= 1 && e <= ens.length ? e : null; }
        const a = Math.max(0, pos - 8), w = uni.slice(a, pos + 7);
        const i = ens.indexOf(w);
        return i >= 0 && ens.indexOf(w, i + 1) < 0 ? i + 1 + (pos - 1 - a) : null;
      };
      const at = (pos) => { const e = toEns(pos); return e == null ? null : tgt[colOf[e]]; };
      return { ...o, at, pid: o.h.target.perc_id, prot: o.h.target.protein_id };
    });
    return { rows, sites: p.sites.map((j) => S[j].pos) };
  }

  function consSvg(p, u, model, pos, flank) {
    const pal = EXPORT_PAL, cw = 15, rh = 17, labelW = 120, x0 = 16 + labelW, ncol = 2 * flank + 1;
    const W = Math.max(560, x0 + ncol * cw + 90);
    const sitesSet = new Set(model.sites);
    let g = '<text x="16" y="28" font-size="17" font-weight="700" fill="' + pal.ink + '">' + esc(p.label + ' ' + rl(pos, p) + ' across species') + '</text>' +
      '<text x="16" y="47" font-size="11" fill="' + pal.muted + '">' + esc(p.acc + ' · Ensembl Compara orthologues, columns anchored on the human sequence · residues ' + Math.max(1, pos - flank) + '–' + (pos + flank)) + '</text>';
    let y = 64;
    for (let k = 0; k < ncol; k++) {
      const q = pos - flank + k;
      if (q === pos || q % 5 === 0) g += '<text x="' + (x0 + k * cw + cw / 2) + '" y="' + (y + 9) + '" font-size="9" text-anchor="middle" fill="' + (q === pos ? pal.ink : pal.muted) + '"' + (q === pos ? ' font-weight="700"' : '') + '>' + q + '</text>';
      if (sitesSet.has(q) && q !== pos) g += '<circle cx="' + (x0 + k * cw + cw / 2) + '" cy="' + (y + 14) + '" r="2" fill="' + pal.site + '"/>';
    }
    y += 18;
    const top = y;
    const line = (label, sub, get) => {
      g += '<text x="16" y="' + (y + rh - 5) + '" font-size="11" fill="' + pal.ink + '">' + esc(label) + '</text>';
      for (let k = 0; k < ncol; k++) { const q = pos - flank + k; g += cellSvg(q >= 1 && q <= u.seq.length ? get(q) || '' : '', x0 + k * cw, y + 1, cw, rh - 2); }
      if (sub) g += '<text x="' + (x0 + ncol * cw + 8) + '" y="' + (y + rh - 5) + '" font-size="10" fill="' + pal.muted + '">' + esc(sub) + '</text>';
      y += rh;
    };
    line('Human', '', (q) => u.seq[q - 1]);
    for (const r of model.rows) line(r.name, r.pid != null ? Math.round(r.pid) + '% id.' : '', (q) => r.at(q));
    g += '<rect x="' + (x0 + flank * cw) + '" y="' + (top - 1) + '" width="' + cw + '" height="' + (y - top + 2) + '" fill="none" stroke="' + pal.ink + '" stroke-width="1.2"/>';
    y += 18;
    const leg = chemLegend(16, y, W);
    g += leg.svg;
    y += leg.h;
    g += '<text x="16" y="' + (y + 4) + '" font-size="10" fill="' + pal.muted + '">Dots: other methylarginine sites. % id.: identity of the orthologue to the human protein (Ensembl).</text>';
    const H = y + 14;
    return { W, H, svg: '<svg xmlns="http://www.w3.org/2000/svg" width="' + W + '" height="' + H + '" viewBox="0 0 ' + W + ' ' + H + '" font-family=\'' + pal.font.replace(/'/g, '"') + '\'><rect width="' + W + '" height="' + H + '" fill="' + pal.bg + '"/>' + g + '</svg>' };
  }

  function renderConsCard(p, u) {
    const box = $('#card-cons');
    if (!box) return;
    const head = '<div class="card-head"><div><h3>Conservation across species</h3><p class="sub">The residue aligned to each methylarginine in mouse, chicken, frog, zebrafish, fly, worm and yeast orthologues (Ensembl Compara).' + docLink('fig-ortho') + '</p></div>' +
      '<div class="exp" role="group" aria-label="Conservation export"><button type="button" id="cons-svg" disabled>SVG</button><button type="button" id="cons-png" disabled>PNG</button><button type="button" id="cons-tsv" disabled>TSV</button></div></div>';
    const ensg = u && u.ensg && u.ensg[0];
    if (!u || !u.seq) { box.innerHTML = head + '<p class="note">' + (u === undefined ? 'Loading UniProt…' : 'UniProt sequence unavailable.') + '</p>'; bindDocLinks(box); return; }
    if (!ensg) { box.innerHTML = head + '<p class="note">UniProt lists no Ensembl gene for ' + esc(p.acc) + '.</p>'; bindDocLinks(box); return; }
    box.innerHTML = head + '<div id="cons-body"><button class="btn small" type="button" id="cons-load">Load orthologue alignments</button> <span class="muted small">from Ensembl ' + esc(ensg) + '; usually 5–30 s</span></div>';
    bindDocLinks(box);
    const go2 = () => {
      const body = $('#cons-body', box);
      body.innerHTML = '<p class="note">Loading orthologues from Ensembl Compara…</p>';
      loadOrthologs(ensg).then((orthos) => {
        if (UI.detail !== p.i || !document.body.contains(box)) return;
        if (!orthos.length) { body.innerHTML = '<p class="note">Ensembl lists no orthologue in these species.</p>'; return; }
        const model = orthoModel(p, u, orthos);
        const sites = model.sites;
        const cls = (a) => (a === 'R' ? 'c-r' : a === 'K' ? 'c-k' : a && a !== '-' ? 'c-x' : 'c-gap');
        const per = model.rows.map((r) => { const v = sites.map((q) => r.at(q)); return { r, v, nR: v.filter((a) => a === 'R').length, nMap: v.filter((a) => a != null).length }; });
        UI.consPos = sites.includes(UI.consPos) ? UI.consPos : sites[0];
        const absent = ORTHO_SPECIES.filter(([sp]) => !orthos.some((o) => o.sp === sp)).map(([, n]) => n);
        body.innerHTML = '<div class="cons-sum">' + per.map((x) => '<span class="tag"><b>' + esc(x.r.name) + '</b> R at ' + x.nR + '/' + x.nMap + '</span>').join(' ') +
          (absent.length ? ' <span class="muted small">No orthologue in Ensembl: ' + esc(absent.join(', ')) + '.</span>' : '') + '</div>' +
          '<div class="tbl-wrap cons-wrap"><table class="cons-tbl"><thead><tr><th>Species</th>' + sites.map((q) => '<th><button type="button" class="linkish" data-cpos="' + q + '">' + q + '</button></th>').join('') + '</tr></thead><tbody>' +
          '<tr><th>Human</th>' + sites.map(() => '<td class="c-r">R</td>').join('') + '</tr>' +
          per.map((x) => '<tr><th>' + esc(x.r.name) + ' <span class="muted">' + Math.round(x.r.pid) + '%</span></th>' + x.v.map((a) => '<td class="' + cls(a) + '">' + esc(a == null ? '·' : a) + '</td>').join('') + '</tr>').join('') +
          '</tbody></table></div><p class="note">R: the arginine is conserved; K: conservative change; other letters: substituted; – gap; · position not mapped (the Ensembl protein differs here). Click a position for its alignment.</p>' +
          '<div class="fig-controls"><label class="fig-ctl">Site <select id="cons-pos">' + sites.map((q) => '<option value="' + q + '"' + (q === UI.consPos ? ' selected' : '') + '>' + rl(q) + '</option>').join('') + '</select></label>' +
          '<label class="fig-ctl">Flank <select id="cons-flank">' + [5, 10, 15, 20].map((v) => '<option value="' + v + '"' + (v === UI.consFlank ? ' selected' : '') + '>±' + v + '</option>').join('') + '</select></label></div>' +
          '<div class="fig-preview" id="cons-fig"></div>';
        const draw = () => { $('#cons-fig', body).innerHTML = consSvg(p, u, model, UI.consPos, UI.consFlank).svg; };
        draw();
        $('#cons-pos', body).addEventListener('change', (e) => { UI.consPos = +e.target.value; draw(); });
        $('#cons-flank', body).addEventListener('change', (e) => { UI.consFlank = +e.target.value; draw(); });
        $$('[data-cpos]', body).forEach((b) => b.addEventListener('click', () => { UI.consPos = +b.dataset.cpos; $('#cons-pos', body).value = String(UI.consPos); draw(); }));
        const name = (p.gene || p.acc) + '_orthologue_alignment';
        $('#cons-svg', box).disabled = false; $('#cons-png', box).disabled = false; $('#cons-tsv', box).disabled = false;
        $('#cons-svg', box).onclick = () => exportFigure('svg', consSvg(p, u, model, UI.consPos, UI.consFlank), name + '_' + rl(UI.consPos, p));
        $('#cons-png', box).onclick = () => exportFigure('png', consSvg(p, u, model, UI.consPos, UI.consFlank), name + '_' + rl(UI.consPos, p));
        $('#cons-tsv', box).onclick = () => download(name + '.tsv', ['site\t' + per.map((x) => x.r.name.toLowerCase() + '_residue').join('\t')]
          .concat(sites.map((q, k) => '' + rl(q) + '\t' + per.map((x) => x.v[k] == null ? '' : x.v[k]).join('\t'))).join('\n') + '\n# orthologues: ' + per.map((x) => x.r.name + ' ' + x.r.prot + ' (' + Math.round(x.r.pid) + '% identity)').join('; ') + '; Ensembl Compara via rest.ensembl.org\n');
      }).catch((e) => { if (document.body.contains(box)) $('#cons-body', box).innerHTML = '<p class="note">Could not load orthologues from Ensembl (' + esc(e.message || e) + '). <button class="linkish" type="button" id="cons-retry">Try again</button></p>'; const r = $('#cons-retry', box); if (r) r.addEventListener('click', go2); });
    };
    $('#cons-load', box).addEventListener('click', go2);
    if (orthoCache.has(ensg)) go2();
  }

  // ----------------------------------------------------------- export / copy
  function siteTsv(rows) {
    const head = ['gene', 'uniprot_accession', 'uniprot_reviewed', 'same_site_on_reviewed_entry', 'protein_name', 'site', 'position', 'numbering', 'state', 'maron_s2_state', 'maron_ip', 'sources', 'source_count',
      'window_15', 'window_31', 'motifs', 'disorder_context', 'disorder_edge_distance', 'domain_context', 'nearest_interpro', 'nearest_domain_class',
      'domain_edge_distance', 'domain_side', 'nearest_methyl_r_distance', 'GSK591_decreased', 'GSK591_increased', 'MS023_decreased', 'MS023_increased',
      'maron_s3_measurements', 'condensates', 'functional_classes', 'reported_by', 'independent_reports', 'references', 'reported_by_doi', 'pdb_interface_partners', 'pdb_interface_entries', 'pdb_contacts', 'predicted_interface_partners_pdockq', 'rna_dna_crosslinks'];
    const out = [head.join('\t')];
    for (const s of rows) {
      const p = P[s.p];
      const ipr = s.ipr >= 0 ? M.interpro[s.ipr] : null;
      out.push([
        p.gene, p.acc, p.unrev ? 'no' : 'yes', s.dup != null ? P[S[s.dup].p].acc : '', p.name, '' + rl(s.pos, P[s.p]), s.pos, P[s.p].hist ? 'histone (site without initiator Met; position is UniProt)' : 'UniProt', STATE_LABEL[s.st], bits(s.s2).map((b) => ['Rme1', 'Rme2'][b]).join(';'),
        bits(s.ip).map((b) => ipLabel(M.maron_ips[b])).join(';'), bits(s.src).map((b) => M.sources[b]).join(';'), popcount(s.src),
        s.win ? s.win.slice(8, 23) : '', s.win, bits(s.mot).map((b) => M.motifs[b]).join(';'),
        M.disorder_classes[s.dc], s.dEdge == null ? '' : s.dEdge, M.domain_contexts[s.dctx], ipr ? ipr[0] + ' ' + ipr[1] : '',
        M.domain_classes[s.dcls], s.domDist == null ? '' : s.domDist, s.domSide, s.nn == null ? '' : s.nn,
        s.trf & 1 ? 1 : 0, s.trf & 2 ? 1 : 0, s.trf & 4 ? 1 : 0, s.trf & 8 ? 1 : 0,
        s.tr.map(([ip, t, d, q]) => ipLabel(M.maron_ips[ip]) + ':' + M.maron_treatments[t] + ':' + d + ':' + q).join(';'),
        p.cond.map((k) => M.condensates[k][0]).join(';'), bits(p.func).map((b) => M.func_classes[b]).join(';'),
        bits(s.study).map((b) => M.studies[b][0]).join(';'), popcount(s.study),
        s.refs.map((k) => M.refs[k]).join(' | '), bits(s.study).map((b) => studyRef(b)[1]).filter(Boolean).join(';'),
        ixReady() ? ifaceText(s) : '', (s.ifPdb || []).join(';'), ixReady() ? contactText(s) : '', ixReady() ? predText(s) : '', ixReady() ? xlText(s) : '',
      ].map(cell).join('\t'));
    }
    return out.join('\n') + '\n';
  }

  function proteinTsv(rows) {
    const head = ['gene', 'uniprot_accession', 'uniprot_id', 'protein_name', 'gene_aliases', 'methyl_r_sites_shown', 'methyl_r_sites_total', 'length', 'arginines',
      'shrinkage_log2_enrichment', 'arg_odds_q', 'sites_in_idr', 'condensates', 'disease_families', 'als_ftd_audited_gene', 'functional_classes', 'functional_subclasses',
      'best_window_aa', 'best_window_sites', 'best_window_start', 'best_window_end', 'best_window_q', 'architecture_class'];
    const out = [head.join('\t')];
    for (const p of rows) {
      out.push([
        p.gene, p.acc, p.uid, p.name, p.aliases.join(';'), R.counts[p.i], p.sites.length, p.len == null ? '' : p.len, p.argN == null ? '' : p.argN,
        p.shrink == null ? '' : p.shrink, p.argQ == null ? '' : p.argQ, p.nIdr, p.cond.map((k) => M.condensates[k][0]).join(';'),
        bits(p.dis).map((b) => M.disease_families[b]).join(';'), p.als ? 1 : 0, bits(p.func).map((b) => M.func_classes[b]).join(';'),
        bits(p.sub).map((b) => M.func_subclasses[b]).join(';'),
        p.best ? p.best[0] : '', p.best ? p.best[3] : '', p.best ? p.best[1] : '', p.best ? p.best[2] : '', p.best ? p.best[5] : '',
        p.arch ? p.arch[0] : '',
      ].map(cell).join('\t'));
    }
    return out.join('\n') + '\n';
  }

  const cell = (v) => String(v == null ? '' : v).replace(/[\t\n\r]+/g, ' ');

  function downloadBlob(name, blob) {
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url; a.download = name;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 4000);
  }

  function download(name, text) {
    downloadBlob(name, new Blob([text], { type: 'text/tab-separated-values;charset=utf-8' }));
  }

  // Every exported figure ends with a credit band naming the site and the paper. Takes {svg, W, H} or an SVG string.
  function withCredit(out) {
    if (typeof out === 'string') {
      const m = /width="([\d.]+)" height="([\d.]+)"/.exec(out);
      out = { svg: out, W: +m[1], H: +m[2] };
    }
    const pal = EXPORT_PAL;
    const lines = [];
    // About 5.4 units per character at 9.5; wrap long lines at their separators, then between words (journal panels
    // can be only ~340 units wide).
    const room = Math.max(120, out.W - 32), fits = (t) => t.length * 5.4 <= room;
    for (const line of creditLines()) {
      if (fits(line)) { lines.push(line); continue; }
      let cur = '';
      for (const part of line.split(/ · |, (?=doi:)/)) {
        const next = cur ? cur + ' · ' + part : part;
        if (fits(next)) { cur = next; continue; }
        if (cur) lines.push(cur);
        cur = '';
        for (const word of part.split(' ')) { const t = cur ? cur + ' ' + word : word; if (fits(t) || !cur) cur = t; else { lines.push(cur); cur = word; } }
      }
      if (cur) lines.push(cur);
    }
    const H = out.H + 12 + lines.length * 13;
    const text = lines.map((t, k) => '<text x="16" y="' + (out.H + 14 + k * 13) + '" font-size="9.5" fill="' + pal.muted + '">' + esc(t) + '</text>').join('');
    return {
      W: out.W, H,
      svg: '<svg xmlns="http://www.w3.org/2000/svg" width="' + out.W + '" height="' + H + '" viewBox="0 0 ' + out.W + ' ' + H + '" font-family=\'' + pal.font.replace(/'/g, '"') + '\'>' +
        '<rect width="' + out.W + '" height="' + H + '" fill="' + pal.bg + '"/>' + out.svg.replace(/^<svg\b/, '<svg x="0" y="0"') + text + '</svg>',
    };
  }

  // ----- journal panels. A figure laid out for a printed column: the canvas is narrowed until its smallest labels
  // (8.5 units) print at 6 pt, the SVG states its width in millimetres, and the PNG is rendered at the chosen
  // resolution and records it, so it opens at its printed size. Widths are the usual Cell Press columns.
  const PANEL_WIDTHS = [['', 'Screen layout'], ['85', 'Single column, 85 mm'], ['114', '1.5 columns, 114 mm'], ['174', 'Full width, 174 mm']];
  const PANEL_DPI = [[300, '300 dpi'], [600, '600 dpi'], [1000, '1000 dpi']];
  const panelUnits = (mm) => Math.round(mm / (6 * 25.4 / 72 / 8.5));
  const panelPrefs = (() => { try { return JSON.parse(store.get('mre-panel') || 'null') || {}; } catch (e) { return {}; } })();
  const panelState = { mm: PANEL_WIDTHS.some(([v]) => v === panelPrefs.mm) ? panelPrefs.mm : '', dpi: PANEL_DPI.some(([v]) => v === panelPrefs.dpi) ? panelPrefs.dpi : 600 };
  const panelOpts = () => (panelState.mm ? { W: panelUnits(+panelState.mm), panel: true } : {});
  // On screen a panel is previewed at twice its printed size (96 px per inch).
  const panelPreviewStyle = () => (panelState.mm ? ' style="max-width:' + Math.round(+panelState.mm / 25.4 * 96 * 2) + 'px"' : '');
  function panelControlsHtml() {
    return '<label class="fig-ctl">Export size <select data-panel="mm">' + PANEL_WIDTHS.map(([v, l]) => '<option value="' + v + '"' + (v === panelState.mm ? ' selected' : '') + '>' + esc(l) + '</option>').join('') + '</select></label>' +
      '<label class="fig-ctl">PNG <select data-panel="dpi"' + (panelState.mm ? '' : ' disabled title="Pick a column width first"') + '>' + PANEL_DPI.map(([v, l]) => '<option value="' + v + '"' + (v === panelState.dpi ? ' selected' : '') + '>' + l + '</option>').join('') + '</select></label>';
  }
  function bindPanelControls(root, redraw) {
    $$('[data-panel]', root).forEach((el) => el.addEventListener('change', () => {
      if (el.dataset.panel === 'mm') panelState.mm = el.value; else panelState.dpi = +el.value;
      store.set('mre-panel', JSON.stringify(panelState));
      const d = $('[data-panel="dpi"]', root);
      if (d) d.disabled = !panelState.mm;
      redraw();
    }));
  }
  // Export a figure at the chosen panel width: mm-sized SVG, or a PNG at the chosen dpi with that dpi recorded.
  function exportPanel(kind, out, name) {
    const mm = +panelState.mm;
    if (!mm) { exportFigure(kind, out, name); return; }
    const fig = withCredit(out), tag = '_' + mm + 'mm';
    if (kind === 'svg') {
      const svg = fig.svg.replace(/width="[\d.]+" height="[\d.]+"/, 'width="' + mm + 'mm" height="' + (fig.H * mm / fig.W).toFixed(1) + 'mm"');
      downloadBlob(name + tag + '.svg', new Blob([svg], { type: 'image/svg+xml;charset=utf-8' }));
    } else svgToPng(fig.svg, fig.W, fig.H, (mm / 25.4 * panelState.dpi) / fig.W, name + tag + '_' + panelState.dpi + 'dpi.png', panelState.dpi);
  }
  // Record the resolution in a PNG (pHYs chunk, pixels per metre), replacing one the browser may have written.
  const CRC_TABLE = (() => { const t = new Uint32Array(256); for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xEDB88320 ^ (c >>> 1) : c >>> 1; t[n] = c >>> 0; } return t; })();
  const crc32 = (bytes) => { let c = 0xFFFFFFFF; for (let i = 0; i < bytes.length; i++) c = CRC_TABLE[(c ^ bytes[i]) & 0xFF] ^ (c >>> 8); return (c ^ 0xFFFFFFFF) >>> 0; };
  function pngWithDpi(blob, dpi) {
    return blob.arrayBuffer().then((buf) => {
      const src = new Uint8Array(buf), dv = new DataView(buf), parts = [src.subarray(0, 8)];
      const ppm = Math.round(dpi / 0.0254), chunk = new Uint8Array(21), cv = new DataView(chunk.buffer);
      cv.setUint32(0, 9); chunk.set([0x70, 0x48, 0x59, 0x73], 4); cv.setUint32(8, ppm); cv.setUint32(12, ppm); chunk[16] = 1;
      cv.setUint32(17, crc32(chunk.subarray(4, 17)));
      for (let off = 8; off < src.length;) {
        const len = dv.getUint32(off), type = String.fromCharCode(src[off + 4], src[off + 5], src[off + 6], src[off + 7]), end = off + 12 + len;
        if (type !== 'pHYs') parts.push(src.subarray(off, end));
        if (type === 'IHDR') parts.push(chunk);
        off = end;
      }
      return new Blob(parts, { type: 'image/png' });
    });
  }
  function exportFigure(kind, out, name) {
    const fig = withCredit(out);
    if (kind === 'svg') downloadBlob(name + '.svg', new Blob([fig.svg], { type: 'image/svg+xml;charset=utf-8' }));
    else svgToPng(fig.svg, fig.W, fig.H, 3, name + '.png');
  }

  // Raster figures (the 3D view) get the same credit, drawn on the canvas under the image.
  function pngWithCredit(dataUrl, name, caption) {
    const img = new Image();
    img.onload = () => {
      const lines = [caption].concat(creditLines()).filter(Boolean);
      const fs = Math.max(12, Math.round(img.width / 80)), lh = Math.round(fs * 1.45), pad = Math.round(fs * 0.9);
      const c = document.createElement('canvas');
      c.width = img.width; c.height = img.height + 2 * pad + lh * lines.length;
      const ctx = c.getContext('2d');
      ctx.fillStyle = '#FFFFFF'; ctx.fillRect(0, 0, c.width, c.height);
      ctx.drawImage(img, 0, 0);
      ctx.font = fs + 'px ' + EXPORT_PAL.font; ctx.textBaseline = 'top';
      lines.forEach((t, k) => { ctx.fillStyle = k === 0 && caption ? EXPORT_PAL.ink : EXPORT_PAL.muted; ctx.fillText(t, pad, img.height + pad + k * lh); });
      c.toBlob((b) => { if (b) downloadBlob(name, b); }, 'image/png');
    };
    img.src = dataUrl;
  }

  // Rasterise an exported SVG on a white canvas at `scale`× for PNG download.
  function svgToPng(svgText, w, h, scale, name, dpi) {
    const img = new Image();
    const url = URL.createObjectURL(new Blob([svgText], { type: 'image/svg+xml;charset=utf-8' }));
    img.onload = () => {
      const c = document.createElement('canvas');
      c.width = Math.round(w * scale); c.height = Math.round(h * scale);
      const ctx = c.getContext('2d');
      ctx.fillStyle = '#FFFFFF'; ctx.fillRect(0, 0, c.width, c.height);
      ctx.drawImage(img, 0, 0, c.width, c.height);
      URL.revokeObjectURL(url);
      c.toBlob((b) => { if (!b) return; if (dpi) pngWithDpi(b, dpi).then((b2) => downloadBlob(name, b2)); else downloadBlob(name, b); }, 'image/png');
    };
    img.src = url;
  }

  function copyText(text, btn) {
    const done = (ok) => { const old = btn.textContent; btn.textContent = ok ? 'Copied' : 'Copy failed'; setTimeout(() => { btn.textContent = old; }, 1500); };
    if (navigator.clipboard && navigator.clipboard.writeText) navigator.clipboard.writeText(text).then(() => done(true), () => fallback());
    else fallback();
    function fallback() {
      const ta = document.createElement('textarea');
      ta.value = text; document.body.appendChild(ta); ta.select();
      let ok = false;
      try { ok = document.execCommand('copy'); } catch (e) { ok = false; }
      ta.remove(); done(ok);
    }
  }
})();

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
  const fmtP = (v) => (v == null ? '—' : v === 0 ? '< 1e-300' : v < 0.001 ? v.toExponential(1) : v < 0.01 ? v.toFixed(3) : v.toFixed(2));
  const fmtL2 = (v) => (v == null ? '—' : (v > 0 ? '+' : '') + v.toFixed(2));
  const pct = (a, b) => (b ? Math.round((100 * a) / b) : 0);
  const pctText = (a, b) => (a > 0 && pct(a, b) === 0 ? '<1' : String(pct(a, b)));
  const popcount = (n) => { let c = 0; while (n) { c += n & 1; n >>= 1; } return c; };
  const bits = (mask) => { const out = []; for (let i = 0; mask; i++, mask >>= 1) if (mask & 1) out.push(i); return out; };
  const debounce = (fn, ms) => { let t; return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); }; };
  const r1 = (v) => Math.round(v * 10) / 10;
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
  // Data sources, versions and licences shown on the Data & methods page. [name, used for, version, licence, citation, doi]
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
      ['InterPro', 'Domain context and nearest domain', 'Current release at staging', 'CC0', 'Blum et al., Nucleic Acids Res 53:D444 (2025)', '10.1093/nar/gkae1082'],
      ['MobiDB', 'Disorder context (predicted consensus disorder, prediction-disorder-th_50)', 'Human export of 30 Mar 2026', 'Not stated on mobidb.org (CC BY 4.0 reported elsewhere)', 'Piovesan et al., Nucleic Acids Res 53:D495 (2025)', '10.1093/nar/gkae969'],
      ['CD-CODE', 'Condensate membership (filters and enrichment); per-protein confidence on protein pages', 'Export of 18 Mar 2026; API release 2.3', 'Not stated by CD-CODE', 'Rostam et al., Nat Methods 20:673 (2023)', '10.1038/s41592-023-01831-0'],
      ['DisGeNET (via DOSE)', 'Disease enrichment terms and families', 'As bundled with DOSE / clusterProfiler 4.18 (Bioconductor, R 4.5)', 'CC BY-NC-SA 4.0 (non-commercial)', 'Piñero et al., Nucleic Acids Res (2019 update)', '10.1093/nar/gkz1021'],
      ['clusterProfiler / DOSE', 'Over-representation analysis', 'clusterProfiler 4.18', 'Artistic-2.0 (software)', 'Wu et al., The Innovation 2:100141 (2021); Yu et al., Bioinformatics 31:608 (2015)', '10.1016/j.xinn.2021.100141'],
    ]],
    ['Interactions and structures (exported with scripts/37)', [
      ['STRING', 'Physical partners (score ≥ 0.7 exported; lower scores live)', 'v12.5', 'CC BY 4.0', 'Szklarczyk et al., Nucleic Acids Res 51:D638 (2023)', '10.1093/nar/gkac1000'],
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
    ]],
    ['Methods and software', [
      ['3Dmol.js 2.5.3', '3D viewer', '2.5.3', 'BSD-3-Clause', 'Rego & Koes, Bioinformatics 31:1322 (2015)', '10.1093/bioinformatics/btu829'],
      ['Levy interface classes', 'Core, rim and support', '', '', 'Levy, J Mol Biol 403:660 (2010)', '10.1016/j.jmb.2010.09.028'],
      ['Solvent accessibility', 'Shrake–Rupley; Arg maximum 274 Å²', '', '', 'Shrake & Rupley, J Mol Biol 79:351 (1973); Tien et al., PLoS ONE 8:e80635 (2013)', '10.1016/0022-2836(73)90011-9'],
      ['AlphaFold-Multimer and pDockQ', 'Predicted complexes and their interface score', '', '', 'Evans et al., bioRxiv (2021); Bryant et al., Nat Commun 13:1265 (2022)', '10.1038/s41467-022-28865-w'],
    ]],
  ];

  const PAPER = {
    cite: 'DeAngelo JD, Silverstein AM, Shechter D. Protein arginine methylation rewires molecular interactions. Cell Chem Biol 33(10):1404–1421 (2026).',
    short: 'DeAngelo, Silverstein & Shechter, Cell Chem Biol 2026',
    doi: '10.1016/j.chembiol.2026.09.011',
  };
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
    single: { name: 'Single colour', cats: ['Methylarginine'], cat: () => 0, colors: (pal) => [pal.site] },
  };
  const TRACKS = [
    ['dense', 'Dense windows'], ['domains', 'Domains'], ['regions', 'Regions'], ['disorder', 'Disorder'],
    ['bias', 'Composition'], ['rg', 'RG/RGG-rich'], ['args', 'All arginines'], ['variants', 'Variants'], ['ptms', 'Other PTMs'], ['iface', 'PDB interfaces'], ['xl', 'RNA/DNA crosslinks'],
  ];
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
    states: new Set(), ips: new Set(), sources: new Set(), s2: new Set(), multi: false, studies: new Set(), minRep: 0,
    treat: '', motifs: new Set(), motifQ: '', motifWhere: 'site', disorder: new Set(), domctx: new Set(), domclass: '',
    nn: 0, ptmGroup: '', ptmDist: 0,
    func: '', sub: '', cond: '', anyCond: false, disease: '', minSites: 1, dense: false, arch: '', argEnriched: false,
    binds: new Set(), ixScore: 700, partner: '', go: new Set(), goExp: false, iface: new Set(), cpx: '', condConf: 0, modelledPdb: false, pred: 0, galOnly: false, xl: new Set(),
  });
  const cloneFilters = (f) => {
    const out = {};
    for (const k of Object.keys(f)) out[k] = f[k] instanceof Set ? new Set(f[k]) : f[k];
    return out;
  };
  let F = DEFAULT_FILTERS();
  let Q = { raw: '', proteins: null, positions: null, note: '' };
  const UI = {
    view: 'overview', detail: null, highlight: null,
    sortP: { key: 'sites', dir: -1 }, sortS: { key: 'gene', dir: 1 },
    pageP: 50, pageS: 100, zoom: null, openEvidence: new Set(),
    allCond: false, funcExpanded: false, charts: {},
    fig: { color: 'state', labels: 'auto', tracks: new Set(['dense', 'domains', 'regions', 'disorder', 'bias', 'rg', 'args', 'variants', 'ptms', 'iface', 'xl']) },
    ixTab: 'proteins', ixTop: 20, ixMin: 700, ixColor: 'meth', ixAll: false, ixNet: null,
    cxPred: false, cxMin: 3, cxSort: 'q', cxAll: false, refsOpen: new Set(), ctx: 'protein', condThr: 3, pendingContact: null, pendingModel: null,
    gal: { conf: 'confident', pl: 80, kind: 'bond', noPdb: false, all: false },
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

    for (const p of P) {
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
    if (window.matchMedia && window.matchMedia('(max-width: 560px)').matches) $('#q').placeholder = 'Gene, accession, site or gene list';
    $('#year').textContent = new Date().getFullYear();
    renderHero();
    const params = new URLSearchParams(location.search);
    F = filtersFromParams(params);
    renderRail();
    bindGlobal();
    if (params.get('q')) { $('#q').value = params.get('q'); $('#q-clear').hidden = false; Q = parseQuery(params.get('q')); }
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
    if (!raw) return { raw: '', proteins: null, positions: null, note: '' };
    const tokens = raw.split(/[\s,;]+/).filter(Boolean);
    const proteins = new Set();
    const positions = new Map();
    const missing = [];
    const missingSites = [];
    const viaAlias = [];
    let last = null;
    const addPos = (hit, pos) => {
      for (const i of hit) {
        if (!positions.has(i)) positions.set(i, new Set());
        positions.get(i).add(pos);
        proteins.add(i);
        if (!P[i].sites.some((j) => S[j].pos === pos)) missingSites.push(P[i].label + ' R' + pos);
      }
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
    if (tokens.length > 1) notes.push('<b>' + fmt(proteins.size) + '</b> protein' + (proteins.size === 1 ? '' : 's') + ' from ' + fmt(matched) + ' matched term' + (matched === 1 ? '' : 's') + '.');
    if (viaAlias.length) notes.push('Matched by synonym: ' + viaAlias.slice(0, 10).map(esc).join(', ') + '.');
    if (missing.length) notes.push('Not in the methylarginine union: ' + missing.slice(0, 25).map(esc).join(', ') + (missing.length > 25 ? ' and ' + (missing.length - 25) + ' more' : '') + '.');
    if (missingSites.length) notes.push('No methylation reported at ' + missingSites.slice(0, 10).map(esc).join(', ') + '.');
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

  // --------------------------------------------------------------- filtering
  function sitePass(s, f) {
    if (f.states.size && !f.states.has(s.st)) return false;
    if (f.ips.size && !bits(s.ip).some((b) => f.ips.has(b))) return false;
    if (f.s2.size && !bits(s.s2).some((b) => f.s2.has(b))) return false;
    if (f.sources.size && !bits(s.src).some((b) => f.sources.has(b))) return false;
    if (f.multi && popcount(s.src) < 2) return false;
    if (f.studies.size && !bits(s.study).some((b) => f.studies.has(b))) return false;
    if (f.minRep && popcount(s.study) < f.minRep) return false;
    if (f.treat) {
      if (f.treat === 'tested') { if (!s.tr.length) return false; }
      else {
        const [t, dir] = f.treat.split(':');
        if (!(s.trf & (1 << (+t * 2 + (dir === 'down' ? 0 : 1))))) return false;
      }
    }
    if (f.motifs.size && !bits(s.mot).some((b) => f.motifs.has(b))) return false;
    if (f.motifQ) { const m = parseMotif(f.motifQ); if (m && m.toks && !motifMatch(s.win, m, f.motifWhere)) return false; }
    if (f.disorder.size && !f.disorder.has(s.dc)) return false;
    if (f.domctx.size && !f.domctx.has(s.dctx)) return false;
    if (f.domclass !== '' && s.dcls !== +f.domclass) return false;
    if (f.nn && !(s.nn != null && s.nn <= f.nn)) return false;
    if (f.ptmDist) {
      const d = f.ptmGroup === '' ? Math.min.apply(null, s.ptmD) : s.ptmD[+f.ptmGroup];
      if (!(d <= f.ptmDist)) return false;
    }
    if (!ixSitePass(s, f)) return false;
    return true;
  }

  function proteinPass(p, f) {
    if (f.func !== '' && !(p.func & (1 << +f.func))) return false;
    if (f.sub !== '' && !(p.sub & (1 << +f.sub))) return false;
    if (f.condConf) {
      // CD-CODE confidence threshold: for the chosen condensate, or for any condensate.
      const mem = cdMemberships(p);
      if (!mem) return false;
      if (f.cond !== '') { const c = mem.find((x) => x.ci === +f.cond); if (!c || c.score < f.condConf) return false; }
      else if (!mem.some((x) => x.score >= f.condConf)) return false;
    } else {
      if (f.cond !== '' && !p.cond.includes(+f.cond)) return false;
      if (f.anyCond && !p.cond.length) return false;
    }
    if (f.disease !== '') {
      if (f.disease === 'als') { if (!p.als) return false; }
      else if (!(p.dis & (1 << +f.disease))) return false;
    }
    if (f.minSites > 1 && p.sites.length < f.minSites) return false;
    if (f.dense && !p.sig.length) return false;
    if (f.arch !== '' && !(p.arch && p.arch[0] === f.arch)) return false;
    if (f.argEnriched && !(p.argQ != null && p.argQ <= PADJ)) return false;
    if (!ixProteinPass(p, f)) return false;
    return true;
  }

  function filterSites(f, q) {
    const okP = new Uint8Array(P.length);
    for (const p of P) okP[p.i] = (!q.proteins || q.proteins.has(p.i)) && proteinPass(p, f) ? 1 : 0;
    const sites = [];
    const counts = new Uint32Array(P.length);
    for (const s of S) {
      if (!okP[s.p]) continue;
      if (q.positions && q.positions.has(s.p) && !q.positions.get(s.p).has(s.pos)) continue;
      if (!sitePass(s, f)) continue;
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

  // Every active constraint, search included, as {label, group, clear}.
  function constraints() {
    const out = [];
    const add = (group, label, clear) => out.push({ group, label, clear });
    if (Q.raw) add('search', 'Search: ' + (Q.raw.length > 40 ? Q.raw.slice(0, 38) + '…' : Q.raw), () => setQuery(''));
    F.states.forEach((v) => add('methyl', 'State: ' + STATE_LABEL[v], () => F.states.delete(v)));
    F.ips.forEach((v) => add('methyl', ipLabel(M.maron_ips[v]) + ' IP', () => F.ips.delete(v)));
    F.s2.forEach((v) => add('methyl', 'Table S2 ' + ['Rme1', 'Rme2'][v], () => F.s2.delete(v)));
    F.studies.forEach((v) => add('prov', 'Reported by ' + M.studies[v][0], () => F.studies.delete(v)));
    if (F.minRep) add('prov', '≥' + F.minRep + ' independent reports', () => { F.minRep = 0; });
    F.sources.forEach((v) => add('prov', SOURCE_LABEL[M.sources[v]], () => F.sources.delete(v)));
    if (F.multi) add('prov', '≥2 source tables', () => { F.multi = false; });
    if (F.treat) add('treat', TREAT_OPTIONS.find((o) => o[0] === F.treat)[1], () => { F.treat = ''; });
    F.motifs.forEach((v) => add('motif', 'Motif ' + MOTIF_LABEL[M.motifs[v]], () => F.motifs.delete(v)));
    if (F.motifQ) { const m = parseMotif(F.motifQ); if (m && m.toks) add('motif', 'Motif ' + m.text + (F.motifWhere === 'window' ? ' within ±15 aa' : ' at the site'), () => { F.motifQ = ''; }); }
    F.disorder.forEach((v) => add('struct', DISORDER_LABEL[M.disorder_classes[v]], () => F.disorder.delete(v)));
    F.domctx.forEach((v) => add('struct', DOMAIN_LABEL[M.domain_contexts[v]], () => F.domctx.delete(v)));
    if (F.domclass !== '') add('struct', 'Nearest domain: ' + M.domain_classes[+F.domclass], () => { F.domclass = ''; });
    if (F.nn) add('cluster', 'Methyl-R within ' + F.nn + ' aa', () => { F.nn = 0; });
    if (F.ptmDist) add('cluster', (F.ptmGroup === '' ? 'Any PTM' : PTM_LABEL[+F.ptmGroup]) + ' within ±' + F.ptmDist + ' aa', () => { F.ptmDist = 0; });
    if (F.dense) add('cluster', 'Dense methyl-R region', () => { F.dense = false; });
    if (F.arch !== '') add('cluster', F.arch, () => { F.arch = ''; });
    if (F.func !== '') add('protein', 'Class: ' + M.func_classes[+F.func], () => { F.func = ''; });
    if (F.sub !== '') add('protein', 'Subclass: ' + M.func_subclasses[+F.sub], () => { F.sub = ''; });
    if (F.cond !== '') add('protein', 'Condensate: ' + M.condensates[+F.cond][0], () => { F.cond = ''; });
    if (F.anyCond) add('protein', 'In a CD-CODE condensate', () => { F.anyCond = false; });
    if (F.condConf) add('protein', 'CD-CODE confidence ≥ ' + F.condConf + '/5', () => { F.condConf = 0; });
    if (F.disease !== '') add('protein', F.disease === 'als' ? 'ALS/MN/FTD audited gene set' : 'Disease: ' + M.disease_families[+F.disease], () => { F.disease = ''; });
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
    if (!cons.length) {
      box.innerHTML = '<span class="scope-count">Showing all <b>' + fmt(S.length) + '</b> sites on <b>' + fmt(P.length) + '</b> proteins</span><span class="scope-hint">Search or use the filters to narrow this down</span>';
      return;
    }
    box.innerHTML = '<span class="scope-count"><span class="scope-flag">Filtered</span><b>' + fmt(R.sites.length) + '</b> of ' + fmt(S.length) + ' sites · <b>' + fmt(R.prots.length) + '</b> proteins</span>' +
      '<span class="scope-chips">' + cons.map((c, k) => '<span class="af' + (c.group === 'search' ? ' search' : '') + '">' + esc(c.label) + '<button type="button" data-k="' + k + '" aria-label="Remove ' + esc(c.label) + '">×</button></span>').join('') + '</span>' +
      '<button class="btn small quiet scope-link" type="button" id="scope-link" title="Copy a link that reopens this page with the same search and filters">Copy link</button>' +
      '<button class="btn small scope-clear" type="button" id="scope-clear">Clear all</button>';
    $$('.af button', box).forEach((b) => b.addEventListener('click', () => { cons[+b.dataset.k].clear(); afterFilterChange(); }));
    $('#scope-link', box).addEventListener('click', (e) => copyLink(e.currentTarget));
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
    if (['overview', 'proteins', 'sites', 'about'].includes(h)) view = h;
    else if (h) {
      const hit = lookupExact(h);
      if (hit && hit.length) { view = 'protein'; if (UI.detail !== hit[0]) { UI.detail = hit[0]; UI.zoom = null; UI.openEvidence.clear(); } }
    }
    UI.view = view;
    $$('.tab').forEach((t) => t.setAttribute('aria-selected', String(t.dataset.view === view)));
    const ptab = $('#tab-protein');
    ptab.hidden = UI.detail == null;
    if (UI.detail != null) ptab.textContent = P[UI.detail].label;
    $('#layout').classList.toggle('no-rail', view === 'about');
    renderView();
    if (view !== 'protein') UI.highlight = null;
  }

  function renderView() {
    const el = $('#content');
    hideTip();
    if (UI.view === 'overview') renderOverview(el);
    else if (UI.view === 'proteins') renderProteins(el);
    else if (UI.view === 'sites') renderSites(el);
    else if (UI.view === 'about') renderAbout(el);
    else if (UI.view === 'protein') renderProtein(el, P[UI.detail]);
  }

  function openProtein(i, pos) {
    UI.highlight = pos == null ? null : pos;
    if (UI.detail !== i) { UI.detail = i; UI.zoom = null; UI.openEvidence.clear(); }
    go('protein', P[i].acc);
    window.scrollTo({ top: 0 });
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
    ['methyl', 'Methylation evidence', true], ['prov', 'Provenance', true], ['treat', 'PRMT inhibitors', true],
    ['motif', 'Sequence motif', false], ['struct', 'Disorder & domains', false], ['cluster', 'Clustering & crosstalk', false],
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
    M.func_subclasses.forEach((name, k) => { const parent = M.subParent[name] || 'Other'; (groups[parent] = groups[parent] || []).push([k, name]); });
    return '<option value="">Any</option>' + Object.keys(groups).sort().map((g) => '<optgroup label="' + esc(g) + '">' +
      groups[g].sort((a, b) => a[1].localeCompare(b[1])).map(([k, n]) => '<option value="' + k + '"' + (String(current) === String(k) ? ' selected' : '') + '>' + esc(n) + '</option>').join('') + '</optgroup>').join('');
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
        chipRow('studies', M.studyOrder.map((k) => [k, M.studies[k][0] + ' · ' + compact(M.studies[k][2]), M.studies[k][1] || STUDY_HELP[M.studies[k][0]] || '']), F.studies) + '</div>' +
        '<div class="fgroup"><label class="flabel" for="f-minrep">Independent reports per site</label>' + sel('f-minrep', options([[2, '≥2 studies or databases'], [3, '≥3'], [5, '≥5']], F.minRep || '', 'Any'), !!F.minRep) + '</div>' +
        '<div class="fgroup"><span class="flabel">Pipeline source table</span>' + chipRow('sources', M.sources.map((n, k) => [k, SOURCE_LABEL[n], SOURCE_HELP[n]]), F.sources) +
        '<label class="check"><input type="checkbox" id="f-multi"' + (F.multi ? ' checked' : '') + '> In ≥2 source tables at the same position</label></div>';
      case 'treat': return '<div class="fgroup"><label class="flabel" for="f-treat">Maron Table S3 response (adj. p ≤ 0.05)</label>' + sel('f-treat', options(TREAT_OPTIONS, F.treat), !!F.treat) + '</div>';
      case 'motif': return '<div class="fgroup"><span class="flabel">Around the methyl-R (match any)</span>' + chipRow('motifs', M.motifs.map((m, k) => [k, MOTIF_LABEL[m], MOTIF_HELP[m]]), F.motifs) + '</div>' +
        '<div class="fgroup"><label class="flabel" for="f-motifq">Custom motif (2–10 residues)</label>' +
        '<div class="motif-search"><input id="f-motifq" type="text" maxlength="40" autocomplete="off" spellcheck="false" placeholder="e.g. RGG, GXRG, [RK]XR" value="' + esc(F.motifQ) + '"' + (F.motifQ ? ' class="set"' : '') + ' aria-describedby="f-motifmsg">' +
        '<button class="btn small" type="button" id="f-motifgo">Search</button></div>' +
        sel('f-motifwhere', options([['site', 'Covering the methyl-R'], ['window', 'Anywhere within ±15 aa']], F.motifWhere), F.motifWhere !== 'site', 'Where the motif may sit') +
        '<p class="note" id="f-motifmsg">' + motifHint() + '</p></div>';
      case 'struct': return '<div class="fgroup"><span class="flabel">Disorder (MobiDB)</span>' + chipRow('disorder', M.disorder_classes.map((c, k) => [k, DISORDER_LABEL[c] || c, DISORDER_HELP[c]]), F.disorder) + '</div>' +
        '<div class="fgroup"><span class="flabel">Domain context (InterPro)</span>' + chipRow('domctx', M.domain_contexts.map((c, k) => [k, DOMAIN_LABEL[c] || c, DOMAIN_HELP[c]]), F.domctx) + '</div>' +
        '<div class="fgroup"><label class="flabel" for="f-domclass">Nearest domain class</label>' + sel('f-domclass', options(M.domain_classes.map((c, k) => [k, c]).sort((a, b) => a[1].localeCompare(b[1])), F.domclass, 'Any'), F.domclass !== '') + '</div>';
      case 'cluster': return '<div class="fgroup"><label class="flabel" for="f-nn">Another methyl-R within</label>' + sel('f-nn', options([[2, '2 aa'], [3, '3 aa'], [5, '5 aa'], [10, '10 aa'], [20, '20 aa']], F.nn || '', 'Any spacing'), !!F.nn) + '</div>' +
        '<div class="fgroup"><label class="flabel" for="f-ptm" title="Sites from the pipeline\'s cross-PTM master table (iPTMnet sources and UniProt). Not exhaustive: sites reported elsewhere may be missing.">Other PTM nearby (PTM master table; not exhaustive)</label>' + sel('f-ptm', options(PTM_LABEL.map((l, k) => [k, l]), F.ptmGroup, 'Any PTM type'), !!F.ptmDist) +
        sel('f-ptmdist', options([[2, 'within ±2 aa'], [5, 'within ±5 aa'], [10, 'within ±10 aa'], [20, 'within ±20 aa']], F.ptmDist || '', 'Off'), !!F.ptmDist, 'PTM distance') + '</div>' +
        '<div class="fgroup"><label class="check"><input type="checkbox" id="f-dense"' + (F.dense ? ' checked' : '') + '> Protein has a significant dense methyl-R region (multiscale null, q ≤ 0.05)</label></div>' +
        '<div class="fgroup"><label class="flabel" for="f-arch">Clustered architecture (≥3 sites in 20 aa)</label>' + sel('f-arch', options(archClasses.map((a) => [a, a]), F.arch, 'Any'), F.arch !== '') + '</div>';
      case 'protein': return '<div class="fgroup"><label class="flabel" for="f-func">Functional class</label>' + sel('f-func', options(M.func_classes.map((c, k) => [k, c]).sort((a, b) => a[1].localeCompare(b[1])), F.func, 'Any'), F.func !== '') + '</div>' +
        '<div class="fgroup"><label class="flabel" for="f-sub">Functional subclass</label>' + sel('f-sub', subclassOptions(F.sub), F.sub !== '') + '</div>' +
        '<div class="fgroup"><label class="flabel" for="f-cond">CD-CODE condensate (methylated members)</label>' + sel('f-cond', options(condOpts, F.cond, 'Any'), F.cond !== '') +
        '<label class="check"><input type="checkbox" id="f-anycond"' + (F.anyCond ? ' checked' : '') + '> In any CD-CODE condensate</label>' +
        '<label class="flabel" for="f-condconf" title="CD-CODE\'s 1–5 confidence that the protein belongs to the condensate">Membership confidence (CD-CODE)</label>' +
        sel('f-condconf', options([[2, '≥ 2 of 5'], [3, '≥ 3 of 5'], [4, '≥ 4 of 5'], [5, '5 of 5']], F.condConf || '', 'Any'), !!F.condConf) + '</div>' +
        '<div class="fgroup"><label class="flabel" for="f-disease">Disease family (DisGeNET)</label>' + sel('f-disease', options(M.disease_families.map((c, k) => [k, c]).concat([['als', 'ALS / MN / FTD (full audited gene set)']]), F.disease, 'Any'), F.disease !== '') + '</div>' +
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
    // The motif box mirrors the filter (cleared from a chip, or restored from a shared link) unless being edited.
    const mq = $('#f-motifq');
    if (mq && document.activeElement !== mq && mq.value.toUpperCase().replace(/\s+/g, '') !== F.motifQ) { mq.value = F.motifQ; mq.classList.toggle('set', !!F.motifQ); }
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
    bind('f-minrep', (el) => { F.minRep = +el.value || 0; });
    bind('f-treat', (el) => { F.treat = el.value; });
    bind('f-domclass', (el) => { F.domclass = el.value; });
    bind('f-nn', (el) => { F.nn = +el.value || 0; });
    bind('f-ptm', (el) => { F.ptmGroup = el.value; if (!F.ptmDist) { F.ptmDist = 5; $('#f-ptmdist').value = '5'; $('#f-ptmdist').classList.add('set'); } });
    bind('f-ptmdist', (el) => { F.ptmDist = +el.value || 0; });
    bind('f-dense', (el) => { F.dense = el.checked; });
    bind('f-arch', (el) => { F.arch = el.value; });
    bind('f-func', (el) => { F.func = el.value; });
    bind('f-sub', (el) => { F.sub = el.value; });
    bind('f-cond', (el) => { F.cond = el.value; });
    bind('f-anycond', (el) => { F.anyCond = el.checked; });
    bind('f-condconf', (el) => { F.condConf = +el.value || 0; });
    bind('f-disease', (el) => { F.disease = el.value; });
    bind('f-min', (el) => { F.minSites = Math.max(1, Math.min(100, +el.value || 1)); }, 'input');
    bind('f-arg', (el) => { F.argEnriched = el.checked; });
    const mq = $('#f-motifq', rail);
    if (mq) {
      const runMotif = () => {
        const m = parseMotif(mq.value);
        F.motifQ = m && m.toks ? m.text : '';
        mq.classList.toggle('set', !!F.motifQ);
        mq.setAttribute('aria-invalid', String(!!(m && m.error)));
        UI.pageP = 50; UI.pageS = 100; compute(); renderView();
        const msg = $('#f-motifmsg', rail);
        if (msg) msg.innerHTML = m && m.error ? '<span class="warn-text">' + esc(m.error) + '</span>' : motifHint();
      };
      mq.addEventListener('keydown', (e) => { if (e.key === 'Enter') { e.preventDefault(); runMotif(); } });
      mq.addEventListener('change', runMotif);
      $('#f-motifgo', rail).addEventListener('click', runMotif);
    }
    bind('f-motifwhere', (el) => { F.motifWhere = el.value; });
    bindIxRail(rail, () => { UI.pageP = 50; UI.pageS = 100; compute(); renderView(); });
    $('#reset-filters').addEventListener('click', () => { const q = Q.raw; F = DEFAULT_FILTERS(); afterFilterChange(); if (q) $('#q').value = q; });
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
    const drawList = () => {
      items = suggestions(input.value);
      if (!items.length || input.value.split(/[\s,;]+/).filter(Boolean).length > 1) { closeList(); return; }
      list.innerHTML = items.map(([p, alias], k) => '<li role="option" id="sg-' + k + '" aria-selected="' + (k === sel) + '" data-i="' + p.i + '"><b>' + esc(p.label) + '</b><span>' +
        (alias ? 'synonym ' + esc(alias) + ' · ' : '') + esc(p.name || p.acc) + '</span><em>' + p.sites.length + ' site' + (p.sites.length === 1 ? '' : 's') + '</em></li>').join('');
      list.hidden = false;
      input.setAttribute('aria-expanded', 'true');
    };
    const apply = debounce(() => {
      setQuery(input.value);
      UI.pageP = 50; UI.pageS = 100;
      compute();
      if (UI.view === 'overview' && Q.raw) go('proteins');
      else renderView();
    }, 220);
    input.addEventListener('input', () => { $('#q-clear').hidden = !input.value; sel = -1; drawList(); apply(); });
    input.addEventListener('keydown', (e) => {
      if (e.key === 'ArrowDown' && !list.hidden) { sel = Math.min(items.length - 1, sel + 1); drawList(); e.preventDefault(); }
      else if (e.key === 'ArrowUp' && !list.hidden) { sel = Math.max(-1, sel - 1); drawList(); e.preventDefault(); }
      else if (e.key === 'Escape') { if (!list.hidden) closeList(); else if (input.value) { setQuery(''); compute(); renderView(); } }
      else if (e.key === 'Enter') {
        e.preventDefault();
        if (sel >= 0 && items[sel]) { pick(items[sel][0]); return; }
        closeList();
        setQuery(input.value);
        compute();
        if (R.prots.length === 1) openProtein(R.prots[0].i, Q.positions && Q.positions.size ? Array.from(Q.positions.values())[0].values().next().value : null);
        else if (Q.positions) go('sites');
        else go('proteins');
      }
    });
    input.addEventListener('blur', () => setTimeout(closeList, 150));
    list.addEventListener('mousedown', (e) => {
      const li = e.target.closest('li');
      if (li) { e.preventDefault(); pick(P[+li.dataset.i]); }
    });
    const pick = (p) => {
      closeList();
      const tokens = input.value.split(/[\s,;]+/).filter(Boolean);
      tokens[Math.max(0, tokens.length - 1)] = p.gene || p.acc;
      setQuery(tokens.join(' '));
      compute();
      openProtein(p.i);
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
    return 'X, . or * = any residue; [ST] = S or T. Press Enter or Search.' + (n ? ' ' + fmt(n) + ' sites without a sequence window are left out.' : '');
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

  function resultsHead(title, sub, actions) {
    return '<div class="results-head"><div><h2>' + title + '</h2><p>' + sub + '</p></div><div class="actions">' + (actions || '') + '</div></div>';
  }

  // Summary card for a protein-level selection (condensate, class, disease or study) above the results.
  function contextCards() {
    const ov = M.overview;
    const cards = [];
    const card = (eyebrow, title, body, links) => cards.push('<div class="ctx-card"><div><span class="eyebrow">' + esc(eyebrow) + '</span><h3>' + esc(title) + '</h3><p>' + body + '</p></div>' + (links ? '<div class="actions">' + links + '</div>' : '') + '</div>');
    if (F.cond !== '') {
      const c = M.condensates[+F.cond];
      card('CD-CODE condensate · ' + c[2], c[0], fmt(c[4]) + ' of ' + fmt(c[3]) + ' members (' + pct(c[4], c[3]) + '%) carry a methylarginine. ' +
        (c[5] != null ? 'Methylated fraction vs the global rate: log2 OR <b>' + fmtL2(c[5]) + '</b>, q ' + fmtP(c[6]) + '.' : 'Too few members for the condensate enrichment test.'),
        '<a class="btn quiet small" href="https://cd-code.org/condensates/' + encodeURIComponent(c[1]) + '" target="_blank" rel="noopener">CD-CODE ' + esc(c[1]) + ' ↗</a>');
    }
    if (F.func !== '') {
      const name = M.func_classes[+F.func];
      const r = ov.functional_classes.find((x) => x[0] === name);
      if (r) card('Functional class', name, fmt(r[1]) + ' methylarginine proteins (' + (100 * r[2]).toFixed(1) + '% of the set vs ' + (100 * r[3]).toFixed(1) + '% of the reviewed background); log2 OR <b>' + fmtL2(r[4]) + '</b>, q ' + fmtP(r[5]) + '.');
    }
    if (F.sub !== '') {
      const name = M.func_subclasses[+F.sub];
      const r = (ov.functional_subclasses || []).find((x) => x[0] === name);
      if (r) card('Functional subclass of ' + (r[6] || '—'), name, fmt(r[1]) + ' methylarginine proteins; log2 OR <b>' + fmtL2(r[4]) + '</b>, q ' + fmtP(r[5]) + ' vs the reviewed background.');
    }
    if (F.disease !== '') {
      if (F.disease === 'als') card('Audited disease family', 'ALS / motor neuron / FTD (full gene set)', 'Every gene in the audited ALS/MN/FTD family, not only genes in enriched terms.');
      else {
        const name = M.disease_families[+F.disease];
        const r = ov.disease_families.find((x) => x[0] === name);
        if (r) card('Audited disease family (DisGeNET)', name, fmt(r[1]) + ' methylarginine genes among ' + fmt(r[2]) + ' family genes; OR <b>' + (r[3] || 0).toFixed(2) + '</b>, BH-adjusted p ' + fmtP(r[4]) + '. Proteins listed here belong to an enriched term of this family.');
      }
    }
    if (F.studies.size === 1) {
      const k = Array.from(F.studies)[0];
      const st = M.studies[k];
      card('Provenance', st[0], esc(st[1] && st[1] !== 'This study.' ? st[1] : STUDY_HELP[st[0]] || '') + ' · ' + fmt(st[2]) + ' sites, ' + fmt(st[3]) + ' reported only here.');
    }
    return cards.length ? '<div class="ctx-cards">' + cards.join('') + '</div>' : '';
  }

  // When nothing matches, show what each constraint would give back on its own.
  function emptyState() {
    const cons = constraints();
    if (!cons.length) return '<div class="card empty"><b>No sites match</b>Broaden the search.</div>';
    const saveF = cloneFilters(F), saveQ = Q;
    const tries = cons.map((c, k) => {
      c.clear();
      const n = filterSites(F, Q).sites.length;
      F = cloneFilters(saveF); Q = saveQ;
      return [k, c.label, n];
    }).filter((t) => t[2] > 0).sort((a, b) => b[2] - a[2]).slice(0, 5);
    return '<div class="card empty"><b>No sites match all ' + cons.length + ' constraints</b>' +
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
  const P_SORT = {
    gene: (p) => p.label.toLowerCase(),
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
    const th = (k, label, cls) => '<th' + (cls ? ' class="' + cls + '"' : '') + ' data-sort="' + k + '"' + (key === k ? ' aria-sort="' + (dir > 0 ? 'ascending' : 'descending') + '"' : '') + '>' + label + '</th>';
    el.innerHTML = resultsHead('Proteins', fmt(rows.length) + ' of ' + fmt(P.length) + ' proteins · ' + fmt(R.sites.length) + ' sites',
      '<button class="btn quiet small" type="button" id="copy-genes">Copy gene list</button><button class="btn small" type="button" id="dl-proteins">Download TSV</button>') +
      contextCards() +
      (rows.length ? '<div class="card tbl-card"><div class="tbl-wrap"><table class="tbl"><thead><tr>' +
        th('gene', 'Protein') + '<th class="hide-sm">Name</th>' + th('sites', 'Sites', 'r') + '<th class="hide-xs">Site map</th>' +
        th('idr', 'In IDR', 'r hide-xs') + th('cond', 'Condensates', 'r hide-sm') + th('density', 'Densest window', 'hide-sm') + th('arg', 'Arg enrichment', 'r hide-sm') +
        '</tr></thead><tbody>' +
        shownRows.map((p) => {
          const n = R.counts[p.i];
          const best = p.best && p.best[5] != null && p.best[5] <= PADJ ? p.best[3] + ' in ' + p.best[0] + ' aa' : '<span class="muted">—</span>';
          return '<tr class="row" tabindex="0" data-i="' + p.i + '"><td><span class="gene">' + esc(p.label) + '</span><span class="acc">' + esc(p.acc) + '</span></td>' +
            '<td class="hide-sm"><div class="pname" title="' + esc(p.name) + '">' + esc(p.name || '—') + '</div></td>' +
            '<td class="r num">' + (n === p.sites.length ? fmt(n) : fmt(n) + '<span class="muted"> / ' + fmt(p.sites.length) + '</span>') + '</td>' +
            '<td class="hide-xs">' + miniMap(p, R.shown) + '</td>' +
            '<td class="r num hide-xs">' + pct(p.nIdr, p.sites.length) + '%</td>' +
            '<td class="r num hide-sm">' + (p.cond.length || '<span class="muted">0</span>') + '</td>' +
            '<td class="hide-sm">' + best + '</td>' +
            '<td class="r num hide-sm">' + (p.shrink == null ? '<span class="muted">—</span>' : fmtL2(p.shrink)) + '</td></tr>';
        }).join('') +
        '</tbody></table></div>' +
        '<div class="tbl-foot"><span>Showing ' + fmt(shownRows.length) + ' of ' + fmt(rows.length) + '</span>' +
        (rows.length > shownRows.length ? '<button class="btn quiet small" type="button" id="more">Show 50 more</button>' : '') + '</div></div>' : emptyState());
    bindEmpty(el);
    bindRows(el, (tr) => openProtein(+tr.dataset.i));
    bindSort(el, UI.sortP, () => renderProteins(el));
    const more = $('#more', el);
    if (more) more.addEventListener('click', () => { UI.pageP += 50; renderProteins(el); });
    $('#dl-proteins', el).addEventListener('click', () => download('methylarginine_proteins.tsv', proteinTsv(rows)));
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
  const S_SORT = {
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
    const th = (k, label, cls) => '<th' + (cls ? ' class="' + cls + '"' : '') + ' data-sort="' + k + '"' + (key === k ? ' aria-sort="' + (dir > 0 ? 'ascending' : 'descending') + '"' : '') + '>' + label + '</th>';
    el.innerHTML = resultsHead('Sites', fmt(rows.length) + ' of ' + fmt(S.length) + ' methylarginine sites on ' + fmt(R.prots.length) + ' proteins',
      '<button class="btn small" type="button" id="dl-sites">Download TSV</button>') +
      contextCards() +
      (rows.length ? '<div class="card tbl-card"><div class="tbl-wrap"><table class="tbl"><thead><tr>' +
        th('gene', 'Protein') + th('pos', 'Site') + th('state', 'State') + th('win', 'Window ±7') + th('mot', 'Motifs', 'hide-sm') +
        th('ctx', 'Context', 'hide-xs') + th('nn', 'Next methyl-R', 'r hide-sm') + th('ip', 'Maron IP', 'hide-sm') + th('prmti', 'PRMTi') + th('iface', 'Interfaces & binding', 'hide-sm') + th('src', 'Reported by', 'hide-sm') +
        '</tr></thead><tbody>' +
        shownRows.map((s) => {
          const p = P[s.p];
          return '<tr class="row" tabindex="0" data-j="' + s.j + '"><td><span class="gene">' + esc(p.label) + '</span><span class="acc">' + esc(p.acc) + '</span></td>' +
            '<td class="num"><b>R' + s.pos + '</b></td><td>' + pill(s.st) + '</td><td>' + windowHtml(s.win, 7) + (s.mis ? ' <span class="warn" title="Canonical residue at this position is not R">!</span>' : '') + '</td>' +
            '<td class="nw hide-sm">' + motifTags(s.mot) + '</td><td class="nw hide-xs">' + contextHtml(s) + '</td>' +
            '<td class="r num hide-sm">' + (s.nn == null ? '<span class="muted">—</span>' : s.nn + ' aa') + '</td>' +
            '<td class="nw hide-sm">' + (ipTags(s.ip) || '<span class="muted">—</span>') + '</td><td class="nw">' + treatHtml(s) + '</td>' +
            '<td class="hide-sm"><div class="tags">' + ((ifaceTags(s, 2) + modelledTag(s) + predTag(s) + (ixReady() ? xlTags(s) : '')) || '<span class="muted">—</span>') + '</div>' + contactHtml(s, 2) + '</td>' +
            '<td class="nw hide-sm"><div class="tags">' + studyTags(s.study, 2) + '</div></td></tr>';
        }).join('') +
        '</tbody></table></div><div class="tbl-foot"><span>Showing ' + fmt(shownRows.length) + ' of ' + fmt(rows.length) + '</span>' +
        (rows.length > shownRows.length ? '<button class="btn quiet small" type="button" id="more">Show 100 more</button>' : '') + '</div></div>' : emptyState());
    bindEmpty(el);
    bindRows(el, (tr) => { const s = S[+tr.dataset.j]; openProtein(s.p, s.pos); });
    bindContactLinks(el);
    bindSort(el, UI.sortS, () => renderSites(el));
    const more = $('#more', el);
    if (more) more.addEventListener('click', () => { UI.pageS += 100; renderSites(el); });
    $('#dl-sites', el).addEventListener('click', () => download('methylarginine_sites.tsv', siteTsv(rows)));
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
    const H = top + rows.length * rowH + 64;
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
    g += '<text x="' + (W - 16) + '" y="' + (H - 10) + '" font-size="9" fill="' + pal.muted + '" text-anchor="end">Methylarginine Site Explorer · ' + esc(PAPER.short) + '</text>';
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
      const out = chartSvg(spec);
      if (b.dataset.exp === 'svg') downloadBlob(name + '.svg', new Blob([out.svg], { type: 'image/svg+xml;charset=utf-8' }));
      else svgToPng(out.svg, out.W, out.H, 3, name + '.png');
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
    const examples = [['FUS', 'One protein'], ['FUS R216', 'One site'], ['TAF15 EWSR1 FUS', 'A gene list (FET family)'], ['SNRPB SNRPD1 SNRPD3', 'Sm proteins'], ['TDP43', 'A synonym'], ['heterogeneous nuclear ribonucleoprotein', 'Protein name']];
    const geneChips = (list, sub) => '<div class="gene-chips">' + list.map((p) => '<button class="gene-chip" type="button" data-i="' + p.i + '">' + esc(p.label) + '<span>' + sub(p) + '</span></button>').join('') + '</div>';
    const sigq = (q) => q != null && q <= PADJ;

    const study = {
      id: 'ch-study', title: 'Provenance', subHtml: 'Sites reported by each study compiled in Maron et al. 2021 Table S5, ProMetheusDB, and curated UniProt/IEDB records. Hover for sites reported only there.',
      subText: 'Methylarginine sites reported by each study or database', noSig: true, format: (v) => fmt(v), tick: compact, axis: 'methylarginine sites reported', valueLabel: 'sites',
      extraCol: true, extraName: 'reported_only_here',
      rows: M.studyOrder.map((k) => ({ label: M.studies[k][0], value: M.studies[k][2], sig: true, key: k, extra: fmt(M.studies[k][3]) + ' only',
        title: '<b>' + esc(M.studies[k][0]) + '</b><br>' + esc(M.studies[k][1] && M.studies[k][1] !== 'This study.' ? M.studies[k][1] : STUDY_HELP[M.studies[k][0]] || '') + '<br>' + fmt(M.studies[k][2]) + ' sites, ' + fmt(M.studies[k][3]) + ' reported only here' })),
    };
    const repDist = [0, 0, 0, 0, 0];
    for (const s of S) repDist[Math.min(5, popcount(s.study)) - 1]++;
    const motif = {
      id: 'ch-motif', title: 'Sequence motifs', subHtml: 'Motif frequency at methylated vs unmodified arginines in the same proteins. Click a motif to list its sites.',
      subText: 'Methylated vs unmodified arginines in the same proteins (±7 window)',
      rows: ov.motifs.map((r) => ({ label: MOTIF_LABEL[r[0]] + ' · ' + MOTIF_HELP[r[0]], value: r[4], q: r[5], sig: sigq(r[5]), key: r[0], extra: fmt(r[1]),
        title: '<b>' + esc(MOTIF_LABEL[r[0]]) + '</b> (' + esc(MOTIF_HELP[r[0]]) + ')<br>' + fmt(r[1]) + ' methyl-R sites (' + (100 * r[2]).toFixed(1) + '%) vs ' + (100 * r[3]).toFixed(1) + '% of unmodified R in the same proteins<br>log2 OR ' + fmtL2(r[4]) + ', q ' + fmtP(r[5]) })),
      extraName: 'sites',
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
      subHtml: 'CD-CODE condensates with ≥5 counted members: methylated fraction vs the global rate (methylated / members shown). ' +
        (thr ? 'Only memberships with CD-CODE confidence ≥ ' + thr + ' of 5 are counted (CD-CODE ' + esc(IX.meta.cdcode_version || '') + '), so condensates whose member lists rest on single proteomics screens do not dominate. ' + condHistNote() : 'All memberships are counted, including single proteomics lists (confidence 1). ') +
        'Click a bar to list its members.<br>' + thrSel,
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
      id: 'ch-func', title: 'Functional ontology', subHtml: 'Lab functional classes (multi-label) among methylarginine proteins vs the reviewed human background. Expand to see subclasses under the class most of their proteins share.',
      subText: 'Functional classes' + (UI.funcExpanded ? ' and subclasses' : '') + ' vs the reviewed human background', rows: funcRows, extraCol: true, extraName: 'proteins',
      footHtml: '<button class="linkish" type="button" id="func-more" style="margin-top:.5rem">' + (UI.funcExpanded ? 'Hide subclasses' : 'Show ' + subs.length + ' subclasses') + '</button>',
    };
    const dis = {
      id: 'ch-dis', title: 'Disease families', subHtml: 'Audited DisGeNET families (genes counted once per family; BH across 8 families). Enrichment of disease gene sets, not causation.',
      subText: 'Audited DisGeNET disease families vs the reviewed background', axis: 'log2 odds ratio vs reviewed background', extraCol: true, extraName: 'methylarginine_genes_of_family',
      rows: ov.disease_families.map((r) => ({ label: r[0], value: r[3] ? Math.log2(r[3]) : null, q: r[4], sig: sigq(r[4]), key: r[0], extra: fmt(r[1]) + '/' + compact(r[2]),
        title: '<b>' + esc(r[0]) + '</b><br>' + fmt(r[1]) + ' methylarginine genes of ' + fmt(r[2]) + ' family genes<br>OR ' + (r[3] || 0).toFixed(2) + ', BH-adjusted p ' + fmtP(r[4]) })).sort((a, b) => b.value - a.value),
    };
    const ctxRows = (table, labels, help) => table.map((r) => ({ label: labels[r[0]] || r[0], value: r[4], q: r[5], sig: sigq(r[5]), key: r[0], extra: fmt(r[1]),
      title: '<b>' + esc(labels[r[0]] || r[0]) + '</b>' + (help[r[0]] ? ' · ' + esc(help[r[0]]) : '') + '<br>' + fmt(r[1]) + ' sites (' + (100 * r[2]).toFixed(1) + '%) vs ' + (100 * r[3]).toFixed(1) + '% of unmodified R in the same proteins<br>log2 OR ' + fmtL2(r[4]) + ', q ' + fmtP(r[5]) }));
    const disorder = {
      id: 'ch-disorder', title: 'Disorder context', subHtml: 'Methylated vs unmodified arginines in the same proteins (MobiDB intervals). Sites on proteins without a MobiDB annotation are left out because their background has only a few dozen arginines.',
      subText: 'Methylated vs unmodified arginines in the same proteins (MobiDB)', extraCol: true, extraName: 'sites',
      rows: ctxRows(ov.disorder.filter((r) => r[0] !== 'no_disorder_annotation'), DISORDER_LABEL, DISORDER_HELP),
    };
    const domctx = {
      id: 'ch-domctx', title: 'Domain context', subHtml: 'Position relative to InterPro domains and repeats, methylated vs unmodified arginines in the same proteins.',
      subText: 'Position relative to InterPro domains, methylated vs unmodified arginines', extraCol: true, extraName: 'sites',
      rows: ctxRows(ov.domain_context, DOMAIN_LABEL, DOMAIN_HELP),
    };
    const domcls = {
      id: 'ch-domcls', title: 'Nearest domain class', subHtml: 'Class of the nearest InterPro domain (lab domain-class ontology), methylated vs unmodified arginines in the same proteins.',
      subText: 'Class of the nearest InterPro domain, methylated vs unmodified arginines', extraCol: true, extraName: 'sites',
      rows: ctxRows(ov.domain_classes || [], {}, {}).sort((a, b) => b.value - a.value),
    };

    el.innerHTML =
      '<div class="ov-intro"><div class="card"><h3>What you can ask</h3>' +
      '<p>This explorer indexes the integrated human methylarginine union from the <a href="https://github.com/Shechterlab/PTM_Pipeline">PTM_Pipeline</a>: ' +
      fmt(M.n_sites) + ' canonically remapped sites on ' + fmt(M.n_proteins) + ' proteins, compiled from Maron et al. 2021, ProMetheusDB and curated UniProt/IEDB annotations, with each site annotated for motif, disorder and domain context, local clustering, nearby PTMs, PRMT inhibitor response, condensate membership and disease association.</p>' +
      '<p>Search for a protein, a site or a gene list, or click a bar below to add a filter. The charts always describe the whole dataset; the bar under the search box shows what is currently filtered.</p>' +
      '<div class="quick"><div><h4>Try</h4><div class="gene-chips">' + examples.map((e) => '<button class="gene-chip" type="button" data-q="' + esc(e[0]) + '">' + esc(e[0]) + '<span>' + esc(e[1]) + '</span></button>').join('') + '</div></div></div></div>' +
      '<div class="card quick"><div><h4>Most methyl-R sites</h4>' + geneChips(topSites, (p) => p.sites.length) + '</div>' +
      '<div><h4>Strongest arginine-methylation enrichment</h4>' + geneChips(topArg, (p) => fmtL2(p.shrink)) + '</div>' +
      '<div><h4>Densest methyl-R clusters</h4>' + geneChips(topDense, (p) => { const w = p.sig.slice().sort((a, b) => b[4] - a[4])[0]; return w[3] + '/' + w[0] + ' aa'; }) + '</div></div></div>' +
      '<div class="ov-grid">' +
      '<div class="card wide chart-card" id="ch-study"><div class="card-head"><div><h3>Provenance</h3><p class="sub">' + study.subHtml + '</p></div>' + expButtons('ch-study') + '</div>' +
      '<div class="two-col"><div>' + hbars(study) + '</div>' +
      '<div><p class="flabel" style="margin:0 0 .5rem">Independent reports per site</p><div class="hbars rep-dist">' +
      repDist.map((n, k) => '<button class="hb" type="button" data-rep="' + (k + 1) + '"><span class="hb-label">' + (k === 4 ? '5 or more' : (k + 1) + (k ? ' reports' : ' report')) + '</span><span class="hb-track"><span class="hb-bar" style="left:0;width:' + (100 * n / repDist[0]).toFixed(1) + '%"></span></span><span class="hb-val">' + fmt(n) + '</span></button>').join('') +
      '</div><p class="note">A report is one study or database. ' + fmt(M.n_sites - repDist[0]) + ' sites (' + pct(M.n_sites - repDist[0], M.n_sites) + '%) have two or more. Click a row to keep sites with at least that many.</p></div></div></div>' +
      chartCard(motif) + chartCard(cond) + chartCard(func) + chartCard(dis) + chartCard(disorder) +
      '<div class="card chart-card" id="ch-dom"><div class="card-head"><div><h3>Domains</h3><p class="sub">' + domctx.subHtml + '</p></div>' + expButtons('ch-domctx') + '</div>' + hbars(domctx) +
      '<div class="card-head" style="margin-top:1.1rem"><div><h3 class="h4">Nearest domain class</h3><p class="sub">' + domcls.subHtml + '</p></div>' + expButtons('ch-domcls') + '</div><div id="ch-domcls-bars">' + hbars(domcls) + '</div></div>' +
      '<div class="card wide" id="ch-treat"><h3>PRMT inhibitor response</h3><p class="sub">Maron et al. 2021 Table S3: sites per methyl-IP that changed significantly (adj. p ≤ 0.05) after GSK591 (PRMT5 inhibitor) or MS023 (type I PRMT inhibitor) vs DMSO. Click a row to list those sites.</p>' + treatTable(ov.treatments) + '</div>' +
      ixOverviewHtml() +
      '</div>';
    UI.charts['ch-study'] = study; UI.charts['ch-domctx'] = domctx; UI.charts['ch-domcls'] = domcls;

    $$('.gene-chip[data-i]', el).forEach((b) => b.addEventListener('click', () => openProtein(+b.dataset.i)));
    $$('.gene-chip[data-q]', el).forEach((b) => b.addEventListener('click', () => {
      setQuery(b.dataset.q); compute();
      if (R.prots.length === 1) openProtein(R.prots[0].i, Q.positions ? Array.from(Q.positions.values())[0].values().next().value : null);
      else go('proteins');
    }));
    bindBars($('#ch-study .hbars:not(.rep-dist)', el), study.rows, (r) => setFilters((f) => { f.studies = new Set([r.key]); }, 'sites'));
    $$('.rep-dist .hb', el).forEach((b) => b.addEventListener('click', () => setFilters((f) => { f.minRep = +b.dataset.rep > 1 ? +b.dataset.rep : 0; }, 'sites')));
    bindBars($('#ch-motif', el), motif.rows, (r) => setFilters((f) => { f.motifs = new Set([M.motifs.indexOf(r.key)]); }, 'sites'));
    bindBars($('#ch-cond', el), cond.rows, (r) => setFilters((f) => { f.cond = String(r.key); f.condConf = thr; }, 'proteins'));
    const ct = $('#cond-thr', el);
    if (ct) ct.addEventListener('change', () => { UI.condThr = +ct.value; UI.allCond = false; renderOverview(el); });
    $('#cond-more', el).addEventListener('click', () => { UI.allCond = !UI.allCond; renderOverview(el); });
    bindBars($('#ch-func', el), func.rows, (r) => setFilters((f) => {
      if (r.kind === 'sub') { f.sub = String(M.func_subclasses.indexOf(r.key)); } else { f.func = String(M.func_classes.indexOf(r.key)); }
    }, 'proteins'));
    $('#func-more', el).addEventListener('click', () => { UI.funcExpanded = !UI.funcExpanded; renderOverview(el); });
    bindBars($('#ch-dis', el), dis.rows, (r) => setFilters((f) => { f.disease = String(M.disease_families.indexOf(r.key)); }, 'proteins'));
    bindBars($('#ch-disorder', el), disorder.rows, (r) => setFilters((f) => { const k = M.disorder_classes.indexOf(r.key); f.disorder = new Set(k >= 0 ? [k] : []); }, 'sites'));
    bindBars($('#ch-dom > .hbars', el), domctx.rows, (r) => setFilters((f) => { const k = M.domain_contexts.indexOf(r.key); f.domctx = new Set(k >= 0 ? [k] : []); }, 'sites'));
    bindBars($('#ch-domcls-bars', el), domcls.rows, (r) => setFilters((f) => { const k = M.domain_classes.indexOf(r.key); f.domclass = k >= 0 ? String(k) : ''; }, 'sites'));
    $$('#ch-treat tr.row', el).forEach((tr) => {
      const go2 = () => setFilters((f) => { f.treat = tr.dataset.t + ':' + tr.dataset.dir; f.ips = new Set([+tr.dataset.ip]); }, 'sites');
      tr.addEventListener('click', go2);
      tr.addEventListener('keydown', (e) => { if (e.key === 'Enter') go2(); });
    });
    bindIxOverview(el);
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
      }).join('') + '</tbody></table></div><p class="note">Median change is the Table S3 <code>&lt;inhibitor&gt;_vs_DMSO_diff</code> value across sites. Rows open the decreased sites.</p>';
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
      '<div class="detail-head"><button class="linkish back" type="button" id="back">← Back to results</button>' +
      '<div class="detail-title"><h2>' + esc(p.label) + '</h2><span class="pn" id="pname">' + esc(p.name) + '</span></div>' +
      (p.aliases.length ? '<p class="aliases">Also reported as ' + p.aliases.map((a) => '<span>' + esc(a) + '</span>').join(', ') + '</p>' : '') +
      '<div class="detail-links"><span class="tag mono">' + esc(p.acc) + '</span>' + (p.uid ? '<span class="tag mono">' + esc(p.uid) + '</span>' : '') +
      (p.iso ? '<span class="tag" title="Sites were reported on this isoform accession">isoform</span>' : '') +
      ext('https://www.uniprot.org/uniprotkb/' + encodeURIComponent(p.acc) + '/entry', 'UniProt') +
      ext('https://alphafold.ebi.ac.uk/entry/' + encodeURIComponent(base), 'AlphaFold') +
      ext('https://www.ebi.ac.uk/interpro/protein/UniProt/' + encodeURIComponent(base) + '/', 'InterPro') +
      '<button class="btn small" type="button" id="save-pdf" title="Print this protein page, or save it as a PDF from the print dialog">Save PDF</button></div></div>' +
      ctxToggleHtml(p) +
      (filtered && shownN < p.sites.length ? '<div class="filter-note">' + fmt(shownN) + ' of ' + fmt(p.sites.length) + ' sites match your filters; the others are faded. <button class="linkish" type="button" id="pf-clear">Clear filters</button></div>' : '') +
      '<div class="facts">' +
      fact(fmt(p.sites.length), 'Methylarginine sites', shownN === p.sites.length ? (filtered ? 'all match your filters' : 'all reported sites') : fmt(shownN) + ' match current filters') +
      fact('<span id="arg-count">' + (argN == null ? '—' : fmt(argN)) + '</span>', 'Arginines', '<span id="arg-frac">' + (argN ? pctText(p.sites.length, argN) + '% methylated' : 'loading from UniProt') + '</span>') +
      fact(p.shrink == null ? '—' : fmtL2(p.shrink), 'Arg enrichment (log2)', p.argQ == null ? 'not tested (needs sequence-confirmed sites)' : 'shrinkage estimate; Fisher q ' + fmtP(p.argQ)) +
      fact(best ? best[3] + ' / ' + best[0] + ' aa' : '—', 'Densest window', best ? 'R' + best[1] + '–' + best[2] + ', q ' + fmtP(best[5]) + (best[5] <= PADJ ? '' : ' (n.s.)') : 'needs ≥3 sequence-confirmed sites') +
      fact('<span id="af-low">…</span>', 'Sites in low-confidence model', '<span id="af-low-sub">AlphaFold pLDDT &lt; 70</span>') +
      fact(fmt(p.cond.length), 'CD-CODE condensates', '<span id="cond-fact-sub">' + (p.cond.length ? '' : 'none annotated') + '</span>') +
      '</div>' +
      figureCard(p) +
      structureCard(p) +
      '<div class="card tbl-card"><div class="card-pad"><h3>Sites</h3><p class="sub">Every reported methylarginine on ' + esc(p.label) + '. Faded rows do not match the current filters. Expand a row for source evidence, UniProt variants and AlphaMissense scores.</p></div>' +
      '<div class="tbl-wrap" id="psites"></div><p class="note uni-extra" id="uni-extra" hidden></p></div>' +
      '<div class="card" id="card-ix"></div>' +
      '<div class="two-col"><div class="card" id="card-dom"></div><div class="card" id="card-cond"></div></div>' +
      '<div class="two-col"><div class="card" id="card-func"></div><div class="card" id="card-dis"></div></div>' +
      '<div class="two-col"><div class="card" id="card-cluster"></div><div class="card" id="card-ptm"></div></div>' +
      '<div class="card"><div class="card-head"><div><h3>Sequence</h3><p class="sub" id="seq-sub">Loading sequence from UniProt…</p></div><div class="exp"><button type="button" id="copy-fasta" disabled>Copy FASTA</button></div></div><div id="seq"></div></div>' +
      '<div class="print-only print-foot"><p>Cite: ' + esc(PAPER.cite) + ' https://doi.org/' + esc(PAPER.doi) + '</p><p>Sites: methylarginine union from PTM_Pipeline. Features, variants and sequence: UniProt. Predicted structure, pLDDT and AlphaMissense: AlphaFold DB. Experimental structures: RCSB PDB and PDBe (SIFTS). Condensates: CD-CODE.</p></div>';

    $('#back', el).addEventListener('click', () => { if (history.length > 1) history.back(); else go('proteins'); });
    $('#save-pdf', el).addEventListener('click', (e) => savePdf(p, e.currentTarget));
    bindCtx(p);
    const pfc = $('#pf-clear', el);
    if (pfc) pfc.addEventListener('click', () => { const q = Q.raw; F = DEFAULT_FILTERS(); afterFilterChange(); if (q) $('#q').value = q; });
    renderSiteTable(p);
    renderDomainCard(p, null);
    renderCondCard(p);
    renderFuncCard(p);
    renderDiseaseCard(p, null);
    renderClusterCard(p);
    renderPtmCard(p);
    renderIxCard(p);
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
      renderDomainCard(p, u);
      renderSiteTable(p);
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

  function renderSiteTable(p) {
    const box = $('#psites');
    if (!box) return;
    const u = UNIDATA.get(p.acc);
    const rows = p.sites.map((j) => S[j]);
    const varsAt = new Map();
    if (u) for (const v of u.variants) if (v.start === v.end) { if (!varsAt.has(v.start)) varsAt.set(v.start, []); varsAt.get(v.start).push(v); }
    box.innerHTML = '<table class="tbl"><thead><tr><th>Site</th><th>State</th><th>Window ±7</th><th class="hide-sm">Motifs</th><th class="hide-xs">Context</th>' +
      '<th class="hide-sm">Nearest domain</th><th class="r hide-sm">Next methyl-R</th><th class="hide-sm" title="Maron methyl-IP and Table S2 calls, and UniProt\'s curated modification">State evidence</th><th>PRMTi</th>' +
      '<th class="hide-xs" title="AlphaFold pLDDT at the site and mean AlphaMissense score for R→X">AlphaFold</th><th class="hide-sm" title="PDB interfaces (PDBe-KB), predicted complexes and RNA/DNA crosslinks in living cells">Interfaces &amp; binding</th><th class="hide-sm">Reported by</th><th><span class="sr-only">Evidence</span></th></tr></thead><tbody>' +
      rows.map((s) => {
        const dim = !R.shown.has(s.j);
        const ipr = s.ipr >= 0 ? M.interpro[s.ipr] : null;
        const dom = ipr ? '<a href="https://www.ebi.ac.uk/interpro/entry/InterPro/' + encodeURIComponent(ipr[0]) + '/" target="_blank" rel="noopener" title="' + esc(ipr[0]) + ' · ' + esc(M.domain_classes[s.dcls]) + '">' + esc(ipr[1]) + '</a><br><span class="muted">' + (s.domDist === 0 ? 'inside' : s.domDist + ' aa ' + (s.domSide === 'N_terminal_to_domain' ? 'N-term' : 'C-term')) + '</span>' : '<span class="muted">—</span>';
        const s2 = bits(s.s2).map((b) => ['Rme1', 'Rme2'][b]).join('/');
        const uniType = u && u.methylR.get(s.pos);
        const vs = (varsAt.get(s.pos) || []).filter((v) => v.from === 'R');
        const open = UI.openEvidence.has(s.j);
        return '<tr class="row' + (dim ? ' dim' : '') + '" tabindex="0" data-j="' + s.j + '" data-pos="' + s.pos + '"><td class="num"><b>R' + s.pos + '</b>' +
          (s.orig ? '<br><span class="muted" title="Position as reported by the source before canonical remapping">reported ' + esc(Array.from(new Set(s.orig.split(';'))).join(', ')) + '</span>' : '') +
          vs.map((v) => '<br><span class="tag var" title="UniProt natural variant ' + esc(v.id) + ': R' + s.pos + v.to + ' · ' + esc(v.desc) + '">R' + s.pos + esc(v.to) + (shortDisease(v.desc) ? ' · ' + esc(shortDisease(v.desc).slice(0, 18)) : '') + '</span>').join('') + '</td>' +
          '<td>' + pill(s.st) + '</td><td>' + windowHtml(s.win, 7) + (s.mis ? ' <span class="warn" title="Canonical residue at this position is not R">!</span>' : '') + '</td>' +
          '<td class="nw hide-sm">' + motifTags(s.mot) + '</td><td class="nw hide-xs">' + contextHtml(s) + '</td><td class="hide-sm">' + dom + '</td>' +
          '<td class="r num hide-sm">' + (s.nn == null ? '—' : s.nn + ' aa') + '</td>' +
          '<td class="nw hide-sm">' + ([ipTags(s.ip), s2 ? '<span class="tag">S2 ' + s2 + '</span>' : '', uniType ? '<span class="tag uni" title="UniProt: ' + esc(uniType) + '">UniProt ' + esc(uniprotMethylLabel(uniType)) + '</span>' : ''].filter(Boolean).join(' ') || '<span class="muted">—</span>') + '</td>' +
          '<td class="nw">' + treatHtml(s) + '</td><td class="nw hide-xs">' + afCell(p, s) + '</td><td class="hide-sm"><div class="tags">' + ((ifaceTags(s, 3) + modelledTag(s) + predTag(s) + (ixReady() ? xlTags(s) : '')) || '<span class="muted">' + (IX ? '—' : '…') + '</span>') + '</div>' + contactHtml(s, 3) + '</td><td class="nw hide-sm"><div class="tags">' + studyTags(s.study, 3) + '</div></td>' +
          '<td class="nw"><button class="btn quiet small view3d" type="button" data-pos="' + s.pos + '" title="Zoom to R' + s.pos + ' in the 3D structure">3D</button> <button class="btn quiet small ev-toggle" type="button" aria-expanded="' + open + '" data-j="' + s.j + '">' + (open ? 'Hide' : 'Evidence') + '</button></td></tr>' +
          (open ? '<tr class="ev-row" data-for="' + s.j + '"><td colspan="13" class="ev" id="ev-' + s.j + '">Loading evidence…</td></tr>' : '');
      }).join('') + '</tbody></table>';
    bindRows(box, (tr) => { UI.highlight = +tr.dataset.pos; flashSite(+tr.dataset.pos, false); });
    bindContactLinks(box);
    $$('.ev-toggle', box).forEach((b) => b.addEventListener('click', () => {
      const j = +b.dataset.j;
      if (UI.openEvidence.has(j)) UI.openEvidence.delete(j); else UI.openEvidence.add(j);
      renderSiteTable(p);
    }));
    $$('.view3d', box).forEach((b) => b.addEventListener('click', () => {
      const pos = +b.dataset.pos;
      UI.highlight = pos;
      const card = $('#card-af');
      if (card) card.scrollIntoView({ block: 'start', behavior: 'smooth' });
      flashSite(pos, false);
    }));
    if (UI.openEvidence.size) loadEvidence().then((ev) => UI.openEvidence.forEach((j) => fillEvidence(j, ev)));
    // UniProt-annotated methylarginines that the union does not contain: a pointer to coverage gaps.
    const extra = $('#uni-extra');
    if (extra && u) {
      const have = new Set(p.sites.map((j) => S[j].pos));
      const missing = Array.from(u.methylR.entries()).filter(([pos]) => !have.has(pos)).sort((a, b) => a[0] - b[0]);
      extra.hidden = !missing.length;
      extra.innerHTML = missing.length ? 'UniProt also annotates ' + missing.length + ' methylarginine' + (missing.length === 1 ? '' : 's') + ' that ' + (missing.length === 1 ? 'is' : 'are') + ' not in the union: ' +
        missing.map(([pos, d]) => '<span class="tag uni" title="' + esc(d) + '">R' + pos + ' ' + esc(uniprotMethylLabel(d)) + '</span>').join(' ') + '. Positions follow the current UniProt sequence.' : '';
    }
    const gap = ixReady() ? IX.modelledGaps.get(p.i) : null;
    if (extra && gap) {
      extra.hidden = false;
      extra.innerHTML += (extra.innerHTML ? '<br>' : '') + 'PDB entries model a methylated arginine at ' + gap.pos.map((x) => '<span class="tag var">R' + x + '</span>').join(' ') + ', which the union does not contain (' +
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
        html += '<p><b>UniProt at R' + s.pos + '</b></p><ul class="term-list">' +
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
          '<td class="mono">' + esc(r[6]) + ':' + esc(r[7]) + (r[7] && String(r[7]) !== String(s.pos) ? ' <span class="muted">→ R' + s.pos + '</span>' : '') + '<br><span class="muted">' + esc(ev.remap_statuses[r[8]].replace(/_/g, ' ')) + '</span></td></tr>'; }).join('') +
        '</tbody></table>' +
        (rows.some((r) => r[9] != null && ev.studies[r[9]] === 'PhosphoSitePlus') ? '<p class="note">' + esc(PSP_NOTICE) + '</p>' : '');
    } else if (!ev) html += '<p class="note">The evidence file did not load.</p>';
    if (s.refs.length && !(rows && rows.length)) html += '<p><b>References:</b> ' + s.refs.map((k) => esc(M.refs[k])).join('; ') + '</p>';
    cellEl.innerHTML = html || '<span class="muted">No additional evidence.</span>';
  }

  function flashSite(pos, scroll) {
    const tr = $('#psites tr.row[data-pos="' + pos + '"]');
    $$('.pmap .head').forEach((h) => {
      const sel = +h.dataset.pos === pos;
      h.setAttribute('stroke-width', sel ? 2.6 : 2);
      h.setAttribute('stroke', sel ? getComputedStyle(document.documentElement).getPropertyValue('--ink').trim() : getComputedStyle(document.documentElement).getPropertyValue('--surface').trim());
    });
    focusStructure(pos);
    if (!tr) return;
    if (scroll) tr.scrollIntoView({ block: 'center', behavior: 'smooth' });
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
    $$('[data-cls]', box).forEach((b) => b.addEventListener('click', () => setFilters((f) => { f.domclass = String(M.domain_classes.indexOf(b.dataset.cls)); }, 'sites')));
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
    $$('button[data-k]', box).forEach((b2) => b2.addEventListener('click', () => setFilters((f) => { f.cond = b2.dataset.k; }, 'proteins')));
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
    $$('[data-func]', box).forEach((b) => b.addEventListener('click', () => setFilters((f) => { f.func = String(M.func_classes.indexOf(b.dataset.func)); }, 'proteins')));
    $$('[data-sub]', box).forEach((b) => b.addEventListener('click', () => setFilters((f) => { f.sub = String(M.func_subclasses.indexOf(b.dataset.sub)); }, 'proteins')));
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
    $$('[data-fam]', box).forEach((b2) => b2.addEventListener('click', (e) => { e.preventDefault(); setFilters((f) => { f.disease = b2.dataset.fam; }, 'proteins'); }));
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
        (p.best ? '<dt>Densest window</dt><dd>' + p.best[3] + ' sites in ' + p.best[0] + ' aa (R' + p.best[1] + '–' + p.best[2] + '), log2 obs/null ' + fmtL2(p.best[4]) + ', q ' + fmtP(p.best[5]) + '</dd>' : '') +
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
        (near.length ? '<p style="margin:.9rem 0 .4rem"><b>Within 10 aa of a methyl-R</b></p><div class="tags">' + near.slice(0, 40).map((n) => '<span class="tag" title="' + n[4] + ' aa from R' + n[3] + '">' + PTM_PREFIX[n[1]] + esc(n[2]) + n[0] + ' · R' + n[3] + '</span>').join('') + '</div>' + (near.length > 40 ? '<p class="note">' + (near.length - 40) + ' more not shown.</p>' : '') : '<p class="note">None within 10 aa of a methylarginine.</p>')
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
        if (s && ch === 'R') line += '<span class="m st' + s.st + rg + '" data-pos="' + (k + 1) + '" title="R' + (k + 1) + ' · ' + STATE_LABEL[s.st] + '">R</span>';
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
  }

  // ------------------------------------------------------------- figure
  function figureCard(p) {
    const f = UI.fig;
    return '<div class="card fig-card"><div class="card-head"><div><h3>Protein map</h3><p class="sub" id="map-sub">Loading UniProt features…</p></div>' +
      '<div class="exp" role="group" aria-label="Figure export"><button type="button" id="zoom-reset" hidden>Reset zoom</button><button type="button" id="fig-svg">SVG</button><button type="button" id="fig-png">PNG</button></div></div>' +
      '<div class="fig-controls"><label class="fig-ctl">Colour sites by <select id="fig-color">' + Object.keys(COLOR_MODES).map((k) => '<option value="' + k + '"' + (f.color === k ? ' selected' : '') + '>' + esc(COLOR_MODES[k].name) + '</option>').join('') + '</select></label>' +
      '<div class="seg" role="group" aria-label="Site labels"><span>Site labels</span>' + [['auto', 'Auto'], ['all', 'All'], ['none', 'None']].map(([k, l]) => '<button type="button" data-lab="' + k + '" aria-pressed="' + (f.labels === k) + '">' + l + '</button>').join('') + '</div></div>' +
      '<div class="chips fig-tracks" aria-label="Tracks">' + TRACKS.map(([k, l]) => '<button class="chip" type="button" data-track="' + k + '" aria-pressed="' + f.tracks.has(k) + '">' + l + '</button>').join('') + '</div>' +
      '<div class="legend" id="fig-legend"></div><div class="map-box" id="map-box"></div>' +
      '<p class="note screen-only">Drag across the map to zoom. SVG and PNG exports reproduce this view on a white background with a title and legend.</p></div>';
  }

  function bindFigure(p) {
    const rerender = () => drawMap(p, UNIDATA.get(p.acc));
    $('#fig-color').addEventListener('change', (e) => { UI.fig.color = e.target.value; rerender(); styleStructure(); renderAfSide(p, AFDATA.get(p.acc)); });
    $$('[data-lab]').forEach((b) => b.addEventListener('click', () => { UI.fig.labels = b.dataset.lab; $$('[data-lab]').forEach((x) => x.setAttribute('aria-pressed', String(x === b))); rerender(); }));
    $$('[data-track]').forEach((b) => b.addEventListener('click', () => {
      const k = b.dataset.track;
      if (UI.fig.tracks.has(k)) UI.fig.tracks.delete(k); else UI.fig.tracks.add(k);
      b.setAttribute('aria-pressed', String(UI.fig.tracks.has(k)));
      rerender();
    }));
    const exportFig = (kind) => {
      const out = mapSvg(p, UNIDATA.get(p.acc), { pal: EXPORT_PAL, forExport: true });
      const name = (p.gene || p.acc) + '_methylarginine_map' + (UI.zoom ? '_' + out.v0 + '-' + out.v1 : '');
      if (kind === 'svg') downloadBlob(name + '.svg', new Blob([out.svg], { type: 'image/svg+xml;charset=utf-8' }));
      else svgToPng(out.svg, out.W, out.H, 3, name + '.png');
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
        g += '<rect class="feat" x="' + r1(a) + '" y="' + headTop + '" width="' + r1(Math.max(2, b - a)) + '" height="' + (base - headTop + 12) + '" rx="3" fill="' + pal.dense + '"' + tip('<b>Dense methyl-R window</b><br>' + w[3] + ' sites in ' + w[0] + ' aa (R' + w[1] + '–' + w[2] + '), log2 obs/null ' + fmtL2(w[4]) + ', q ' + fmtP(w[5])) + '/>';
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
          '<text x="' + r1(lx + 1) + '" y="' + labelY + '" font-size="9.5" fill="' + pal.ink + '" transform="rotate(-55 ' + r1(lx + 1) + ' ' + labelY + ')"' + (dim ? ' opacity=".4"' : '') + '>R' + h.s.pos + '</text>';
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
    const boxTrack = (name, list, fill, ink, labelOf, tipOf, dashed) => {
      if (!list.length) return;
      tlabel(name, y);
      if (dashed) gridLine(y);
      let right = -1e9;
      list.sort((p1, p2) => p1.start - p2.start).forEach((f, k) => {
        const a = cx(f.start), b = cx(f.end), w = Math.max(2, b - a);
        const nextX = k + 1 < list.length ? cx(list[k + 1].start) : W - MR;
        g += '<rect class="feat" x="' + r1(a) + '" y="' + (y - 5) + '" width="' + r1(w) + '" height="10" rx="3" fill="' + fill + '"' +
          (dashed ? ' stroke="' + pal.disEdge + '" stroke-dasharray="3 2"' : '') + tip(tipOf(f)) + '/>';
        const text = labelOf ? labelOf(f) : '';
        const tw = text.length * 5.6;
        if (!text) return;
        if (w > tw + 6) g += '<text x="' + r1(a + 3) + '" y="' + (y + 3.5) + '" font-size="9.5" fill="' + ink + '" pointer-events="none">' + esc(text) + '</text>';
        else if (a + w + 4 > right && a + w + 8 + tw < nextX) { g += '<text x="' + r1(a + w + 4) + '" y="' + (y + 3.5) + '" font-size="9.5" fill="' + ink + '">' + esc(text) + '</text>'; right = a + w + 4 + tw; }
      });
      y += 20;
    };
    if (T.has('regions')) boxTrack('Regions', feats.filter((f) => ((f.type === 'Region' && f.desc !== 'Disordered') || f.type === 'Motif') && visible(f.start, f.end)), pal.regionFill, pal.ink, (f) => f.desc, (f) => '<b>' + esc(f.type) + '</b> ' + esc(f.desc) + '<br>' + f.start + '–' + f.end);
    if (T.has('disorder') && u) {
      const dis = feats.filter((f) => f.type === 'Region' && f.desc === 'Disordered' && visible(f.start, f.end));
      if (dis.length) boxTrack('Disorder', dis, pal.disFill, pal.muted, null, (f) => '<b>Disordered</b> (UniProt)<br>' + f.start + '–' + f.end, true);
      else { tlabel('Disorder', y); gridLine(y); y += 20; }
    }
    if (T.has('bias')) boxTrack('Composition', feats.filter((f) => f.type === 'Compositional bias' && visible(f.start, f.end)), pal.biasFill, pal.ink, (f) => f.desc.replace(/ residues$/, '-rich').replace('Basic and acidic', 'Charged'), (f) => '<b>Compositional bias</b> ' + esc(f.desc) + '<br>' + f.start + '–' + f.end);
    if (T.has('rg') && u) boxTrack('RG/RGG', u.rg.filter((b) => visible(b.start, b.end)).map((b) => Object.assign({}, b)), pal.rgFill, pal.rgInk, (b) => b.label, (b) => '<b>' + b.label + '</b> ' + b.start + '–' + b.end + '<br>' + b.n + ' RG dipeptides, ' + b.rgg + ' RGG<br>' + RG_HELP);
    if (T.has('args') && u && u.seq) {
      tlabel('Arginines', y);
      const methyl = new Set(p.sites.map((j) => S[j].pos));
      for (let k = Math.max(0, v0 - 1); k < Math.min(u.seq.length, v1); k++) {
        if (u.seq[k] !== 'R' || methyl.has(k + 1)) continue;
        g += '<line x1="' + r1(x(k + 1)) + '" x2="' + r1(x(k + 1)) + '" y1="' + (y - 5) + '" y2="' + (y + 5) + '" stroke="' + pal.rtick + '"/>';
      }
      for (const s of sites) g += '<line x1="' + r1(x(s.pos)) + '" x2="' + r1(x(s.pos)) + '" y1="' + (y - 7) + '" y2="' + (y + 7) + '" stroke="' + (colors[mode.cat(s)] || pal.none) + '" stroke-width="2"/>';
      y += 20;
    }
    if (T.has('variants') && u && (u.variants.length || u.mutagen.length)) {
      const methyl = new Set(p.sites.map((j) => S[j].pos));
      const track = (name, list, isMut) => {
        list = list.filter((f) => visible(f.start, f.end));
        if (!list.length) return;
        tlabel(name, y);
        gridLine(y);
        let right = -1e9;
        for (const f of list.sort((a, b) => a.start - b.start)) {
          const hit = f.from === 'R' && f.start === f.end && methyl.has(f.start);
          const xx = r1(x(f.start));
          const label = f.from + f.start + (f.to || '');
          const tipHtml = '<b>' + (isMut ? 'Mutagenesis' : 'Natural variant') + ' ' + esc(label) + '</b>' + (hit ? ' · at a methylarginine' : '') + '<br>' + esc(f.desc || '');
          if (isMut) g += '<rect class="feat" x="' + (xx - 3) + '" y="' + (y - 3) + '" width="6" height="6" transform="rotate(45 ' + xx + ' ' + y + ')" fill="' + (hit ? pal.accent : pal.bg) + '" stroke="' + (hit ? pal.accent : pal.variant) + '" stroke-width="1.2"' + tip(tipHtml) + '/>';
          else g += '<rect class="feat" x="' + (xx - (hit ? 1.5 : 1)) + '" y="' + (y - (hit ? 7 : 4)) + '" width="' + (hit ? 3 : 2) + '" height="' + (hit ? 14 : 8) + '" rx="1" fill="' + (hit ? pal.accent : pal.variant) + '"' + tip(tipHtml) + '/>';
          if (hit && +xx + 4 > right) {
            const t = label + (shortDisease(f.desc) ? ' ' + shortDisease(f.desc).split(' and ')[0] : '');
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
    }
    if (T.has('ptms')) {
      for (let gi = 0; gi < 4; gi++) {
        const list = p.ptms.filter((t) => t[1] === gi && t[0] >= v0 && t[0] <= v1);
        if (!list.length) continue;
        tlabel(PTM_LABEL[gi], y);
        gridLine(y);
        for (const t of list) g += '<rect class="feat" x="' + r1(x(t[0]) - 1.5) + '" y="' + (y - 5) + '" width="3" height="10" rx="1" fill="' + pal.ptm + '"' + tip('<b>' + PTM_PREFIX[gi] + t[2] + t[0] + '</b> · ' + PTM_LABEL[gi]) + '/>';
        y += 18;
      }
    }
    // Residues at interfaces in PDB structures (PDBe-KB), one row per partner type.
    const ifd = T.has('iface') && !p.iso ? IFDATA.get(p.acc.split('-')[0]) : null;
    if (ifd) {
      const methylPos = new Set(p.sites.map((j) => S[j].pos));
      for (const [key, label, color] of [['rna', 'RNA contact', IFACE_COLOR[4]], ['dna', 'DNA contact', IFACE_COLOR[8]], ['protein', 'Protein contact', IFACE_COLOR[1]], ['self', 'Homomeric', IFACE_COLOR[2]]]) {
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
            tip('<b>' + label + '</b> ' + (a === b ? a : a + '–' + b) + '<br>' + esc(Array.from(names).slice(0, 6).join(', ')) + (names.size > 6 ? ' and ' + (names.size - 6) + ' more' : '') + (nm ? '<br>' + nm + ' methylarginine' + (nm === 1 ? '' : 's') + ' here' : '') + '<br>PDBe-KB') + '/>';
        };
        for (let k = 1; k <= positions.length; k++) {
          if (k < positions.length && positions[k] === prev + 1) { prev = positions[k]; continue; }
          flush(start, prev);
          if (k < positions.length) { start = prev = positions[k]; }
        }
        y += 18;
      }
    }
    const xp = T.has('xl') && !p.iso && ixReady() ? IX.xlProt.get(p.i) : null;
    if (xp) {
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
      if (xp.dnaPos.length || xp.dnaPeps.length || xp.dnaBind.length) {
        tlabel('DNA crosslinks', y);
        gridLine(y);
        for (const [a, b] of xp.dnaBind) if (inView(a, b)) g += box(a, b, IFACE_COLOR[8], 0.22, '<b>UniProt DNA-binding region</b> ' + a + '–' + b);
        for (const [a, b] of xp.dnaPeps) if (inView(a, b)) g += box(a, b, IFACE_COLOR[8], 0.45, '<b>DNA-crosslinked peptide</b> ' + a + '–' + b + '<br>Trendel et al. 2025, MCF7');
        for (const pos of xp.dnaPos) if (pos >= v0 && pos <= v1) g += tick(pos, IFACE_COLOR[8], '<b>DNA-crosslinked residue</b> ' + pos + '<br>Trendel et al. 2025, MCF7');
        y += 18;
      }
    }
    // Axis
    y += 4;
    gridLine(y);
    const span = v1 - v0;
    const raw = span / 8;
    const mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const step = [1, 2, 5, 10].map((m) => m * mag).find((s) => s >= raw) || mag * 10;
    for (let t = Math.ceil(v0 / step) * step; t <= v1; t += step) {
      g += '<line x1="' + r1(x(t)) + '" x2="' + r1(x(t)) + '" y1="' + y + '" y2="' + (y + 4) + '" stroke="' + pal.grid + '"/><text x="' + r1(x(t)) + '" y="' + (y + 15) + '" font-size="10" text-anchor="middle" fill="' + pal.muted + '">' + fmt(t) + '</text>';
    }
    tlabel('aa', y + 11);
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
      g += '<text x="' + (W - MR) + '" y="' + (y + 4) + '" font-size="9" text-anchor="end" fill="' + pal.muted + '">Sites: methylarginine union, PTM_Pipeline · features and variants: UniProt' + (T.has('ptms') && p.ptms.length ? ' · other PTMs: PTM master table (not exhaustive)' : '') + (ifd ? ' · interfaces: PDBe-KB' : '') + ' · ' + esc(PAPER.short) + '</text>';
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
      sub.textContent = u ? fmt(out.len) + ' aa · domains, regions, disorder and variants from UniProt' + (UI.fig.tracks.has('ptms') && p.ptms.length ? ' · other PTMs from the PTM master table (not exhaustive)' : '') + (UI.fig.tracks.has('iface') && IFDATA.get(p.acc.split('-')[0]) ? ' · PDB interfaces from PDBe-KB' : '') + (UI.zoom ? ' · showing ' + out.v0 + '–' + out.v1 : '') :
        (uniprotCache.has(p.acc) && UNIDATA.has(p.acc) ? 'UniProt features unavailable; showing sites on a ' + (p.len ? fmt(out.len) + ' aa' : 'length-estimated') + ' backbone.' : 'Loading UniProt features…');
    }
    $('#zoom-reset').hidden = !UI.zoom;
    $('#fig-legend').innerHTML = '<span class="legend-title">' + esc(out.mode.name) + '</span>' + out.mode.cats.map((c, k) => '<span><i style="background:' + out.colors[k] + '"></i>' + esc(c) + '</span>').join('') +
      '<span><i class="sq" style="background:var(--fig-domain)"></i>Domain</span>';
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
      '<p class="flabel screen-only">Focus a site</p><div class="chips af-sites screen-only">' + p.sites.map((j) => S[j]).map((s) => '<button class="chip" type="button" data-pos="' + s.pos + '"><i class="dot" style="background:' + (colors[mode.cat(s)] || pal.none) + '"></i>R' + s.pos + '</button>').join('') + '</div>';
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
        fetch(mol.viewer.pngURI()).then((r) => r.blob()).then((b) => downloadBlob((p.gene || p.acc) + '_' + (mol.mode === 'pdb' ? mol.pdb.e.id.toUpperCase() : 'alphafold') + '_methylarginines.png', b));
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
    ['rainbow', 'Rainbow (N→C)'], ['ss', 'Secondary structure'], ['charge', 'Charge'], ['hydro', 'Hydrophobicity'], ['uniform', 'Uniform grey'],
  ];
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
        if (at) mol.siteLabels.push(viewer.addLabel('R' + S[j].pos, { position: { x: at.x, y: at.y, z: at.z }, backgroundColor: siteColor.get(r), backgroundOpacity: 0.85, fontColor: '#FFFFFF', fontSize: 11, inFront: true }));
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
        const text = pos != null ? 'R' + pos + (pdbMode ? ' (' + atom.chain + atom.resi + ')' : '') + (pl != null ? ' · pLDDT ' + Math.round(pl) : '') : atom.resn + ' ' + atom.chain + atom.resi;
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
    if (msg) msg.textContent = atoms.length ? '' : 'R' + pos + ' is not ' + (pdbMode && mol.pdb.pred ? 'in the ' + mol.pdb.pred.label + ' model' : 'resolved in ' + (pdbMode ? mol.pdb.e.id.toUpperCase() + ' chain ' + mol.pdb.chain : 'the AlphaFold model')) + '.';
    if (!atoms.length) return;
    const a = atoms[0];
    mol.focus = { chain: a.chain, resi: a.resi };
    styleStructure(false);
    mol.label = viewer.addLabel('R' + pos + (pdbMode ? ' (' + mol.pdb.chain + auth + ')' : ''), { position: { x: a.x, y: a.y, z: a.z }, backgroundColor: '#C2410C', backgroundOpacity: 0.95, fontColor: '#FFFFFF', fontSize: 13, inFront: true });
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
          cell.innerHTML = (obs.length ? '<b>' + obs.length + '</b> resolved: ' + obs.slice(0, 8).map((x) => 'R' + x[0]).join(', ') + (obs.length > 8 ? ' …' : '') : '<span class="muted">' + (n ? n + ' in construct, none resolved' : '—') + '</span>') + (mods ? '<br>' + mods : '');
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
    const chips = () => '<div class="chips af-sites">' + sites.map((s) => '<button class="chip" type="button" data-pos="' + s.pos + '"><i class="dot" style="background:' + (colors[mode.cat(s)] || pal.none) + '"></i>R' + s.pos + '</button>').join('') + '</div>';
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
        '<p class="flabel">Aligned methylarginines</p>' + (sites.length ? '<div class="chips af-sites">' + sites.map((s) => '<button class="chip" type="button" data-pos="' + s.pos + '"><i class="dot" style="background:' + (colors[mode.cat(s)] || pal.none) + '"></i>R' + s.pos + ' ≙ ' + esc(det.chain) + toAuth.get(s.pos) + '</button>').join('') + '</div>' :
          '<p class="note">No methylarginine of ' + esc(p.label) + ' aligns to a resolved arginine in this chain.</p>') +
        (det.modres.length ? '<p class="flabel">Methylated arginines in the model</p><div class="tags">' + det.modres.map((r) => '<span class="tag var" title="' + esc(r.chem_comp_name) + '">' + esc(r.chem_comp_id) + ' ' + esc(r.chain_id) + esc(r.author_residue_number) + '</span>').join('') + '</div>' : '') +
        contactsHtml() +
        '<p class="note">A related protein, not ' + esc(p.label) + ' itself: positions come from the RCSB sequence alignment (MMseqs2), and the contacts are those of the homologous arginine.</p>' +
        (mol.afPdb ? '<p><button class="btn quiet small" type="button" id="back-af">Back to AlphaFold model</button></p>' : '');
    } else side.innerHTML = '<p class="flabel">Showing PDB entry</p><dl class="kv"><dt>Entry</dt><dd><a href="https://www.rcsb.org/structure/' + e.id.toUpperCase() + '" target="_blank" rel="noopener">' + e.id.toUpperCase() + ' ↗</a></dd>' +
      '<dt>Method</dt><dd>' + esc(methodShort(e.method)) + (e.resolution ? ', ' + e.resolution.toFixed(1) + ' Å' : '') + (/NMR/i.test(e.method || '') ? ' (first model of the ensemble)' : '') + '</dd><dt>Chain</dt><dd>' + esc(det.chain || '—') + ' (teal; other chains grey)</dd></dl>' +
      (det.modres.length ? '<p class="flabel">Methylated arginines in the model</p><div class="tags">' + det.modres.map((r) => '<span class="tag var" title="' + esc(r.chem_comp_name) + '">' + esc(r.chem_comp_id) + ' ' + esc(r.chain_id) + esc(r.author_residue_number) + '</span>').join('') + '</div><p class="note">Drawn in orange.</p>' : '') +
      '<p class="flabel">Resolved methylarginines from the union</p>' +
      (sites.length ? '<div class="chips af-sites">' + sites.map((s) => '<button class="chip" type="button" data-pos="' + s.pos + '"><i class="dot" style="background:' + (colors[mode.cat(s)] || pal.none) + '"></i>R' + s.pos + ' · ' + esc(det.chain) + toAuth.get(s.pos) + '</button>').join('') + '</div>' : '<p class="note">None of this protein\'s methylarginines are resolved in this chain.</p>') +
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
        return '<li title="' + esc('R' + pos + ' aligns to the arginine at position ' + x.sub + ' of the ' + h.entry + ' chain ' + h.chain + ' sequence') + '"><button class="btn quiet small" type="button" data-hom="' + esc(h.id) + '" data-hpos="' + pos + '">View</button> ' +
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
          return '<tr><td class="nw"><b>R' + pos + '</b><br><span class="muted">' + fmt(hs.length) + ' entr' + (hs.length === 1 ? 'y' : 'ies') + '</span></td><td><ul class="hom-list">' +
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
      '<td class="r num">' + (r.sc / 100).toFixed(2) + '</td><td>' + r.sites.slice(0, 10).map((x) => 'R' + x).join(', ') + (r.sites.length > 10 ? ' +' + (r.sites.length - 10) : '') + '</td>' +
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
  const IXF = { prmt: 1, tudor: 2, rna: 4, dna: 8, rnaExp: 16, dnaExp: 32 };
  const IFACE_OPTS = [[4, 'RNA'], [8, 'DNA'], [1, 'Another protein'], [2, 'Same protein (homomeric)']];
  const IFACE_NAME = { 4: 'RNA', 8: 'DNA', 1: 'another protein', 2: 'the same protein' };
  const IFACE_COLOR = { 4: '#7B4FBF', 8: '#B7791F', 1: '#0E7490', 2: '#8A99A3' };
  const NA_COLOR = { RNA: '#7B4FBF', DNA: '#B7791F' };
  const BIND_OPTS = [
    [1, 'a PRMT', 'PRMT1, PRMT2, PRMT3, CARM1 (PRMT4), PRMT5, PRMT6, PRMT7, PRMT8 or PRMT9'],
    [2, 'a Tudor-domain protein', 'InterPro IPR002999: methylarginine readers such as SMN1/2, SMNDC1, TDRD3, SND1 and TDRKH, and methyllysine readers such as PHF1, SETDB1 and KDM4A'],
    [3, 'another methylated protein', 'A partner that is itself in this methylarginine set'],
  ];
  const IX_SCORES = [[400, 'Medium (≥ 0.4, live from STRING)'], [700, 'High (≥ 0.7)'], [900, 'Highest (≥ 0.9)']];
  const RAIL_SCORES = [[700, 'High (≥ 0.7)'], [900, 'Highest (≥ 0.9)']];
  const PRED_OPTS = [[0, 'Any'], [23, 'pDockQ ≥ 0.23 (acceptable)'], [50, 'pDockQ ≥ 0.5 (confident)']];
  const STRING_CALLER = 'shechterlab_methylarginine_explorer';
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
    for (const [j, r, n, dna, f] of d.xl || []) { if (S[j]) S[j].xl = { r, n, d: dna, f }; }
    const xlProt = new Map();
    for (const [i, rnaPos, rnaBits, dnaPos, dnaPeps, dnaBind, dnaInt, rbd] of d.xl_prot || []) xlProt.set(i, { rnaPos, rnaBits, dnaPos, dnaPeps, dnaBind, dnaInt: !!dnaInt, rbd: rbd || [] });
    const modelledGaps = new Map((d.modelled_gaps || []).map(([k, pos, pdbs]) => [k, { pos, pdbs }]));
    IX = { meta: d.meta, prot, byEnsp, cx, cxOf, enrich: d.enrich || {}, goEnrich: d.go_enrich || {}, strata: d.enrich_strata || {}, cd, wr, nP,
      geoSummary: d.geometry_summary || {}, rg: d.rg_iface || {}, modelledGaps, condConf: d.cond_conf || {}, nModelled: (d.modelled || []).length,
      predEnrich: d.pred_enrich || {}, nPred: (d.pred || []).length, condHist: d.cond_score_hist || {}, gallery,
      xlProt, xlEnrich: d.xl_enrich || {}, xlSets: d.xl_sets || [] };
    IX.byAcc = new Map(prot.map((x) => [x.acc, x.k]));
    IX.prmts = prot.filter((x) => x.flags & IXF.prmt).sort((a, b) => a.gene.localeCompare(b.gene, 'en', { numeric: true }));
    IX.tudors = prot.filter((x) => x.flags & IXF.tudor).sort((a, b) => a.gene.localeCompare(b.gene));
    return IX;
  }

  function afterInteractionsLoaded() {
    renderRail();
    compute();
    if (UI.view === 'overview' || UI.view === 'sites') renderView();
    else if (UI.view === 'protein' && UI.detail != null) {
      const p = P[UI.detail];
      renderIxCard(p);
      renderCondCard(p);
      renderSiteTable(p);
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
    if (s.iface & 2) tags.push('<span class="tag" title="' + title + '">homomeric</span>');
    return tags.join(' ');
  }

  function ifaceText(s) {
    if (!s.iface) return '';
    return (s.ifP || []).map((k) => (k === P[s.p].i ? 'homomeric' : ixGene(k))).join(';');
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
  const pdbButtons = (s, pdbs) => pdbs.slice(0, 3).map((x) => '<button class="linkish pdb-open" type="button" data-cpdb="' + esc(x) + '" data-cpos="' + s.pos + '" data-cp="' + s.p + '" title="Open ' + x.toUpperCase() + ' in the 3D viewer at R' + s.pos + '">' + x.toUpperCase() + '</button>').join(', ') +
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
  const XL_OPTS = [[1, 'RNA-crosslinked residue'], [2, 'Within 3 aa of an RNA crosslink'], [4, 'RBDmap RNA-bound fragment'], [8, 'DNA-crosslinked residue'], [16, 'UniProt DNA-binding region']];
  const xlSets = (bits) => (IX.xlSets || []).filter((_, k) => bits & (1 << k)).concat(bits & 32 ? ['pepR-MS focused site'] : []);
  function xlMask(s) {
    const x = s.xl;
    if (!x) return 0;
    return (x.r ? 1 : 0) | (x.n ? 2 : 0) | (x.f ? 4 : 0) | (x.d & 1 ? 8 : 0) | (x.d & 4 ? 16 : 0);
  }
  function xlTags(s) {
    const x = s.xl;
    if (!x) return '';
    const t = [];
    if (x.r) t.push('<span class="tag na-rna" title="' + esc('RNA-crosslinked residue in living human cells (UV crosslinking mass spectrometry): ' + xlSets(x.r).join(', ')) + '">RNA xl</span>');
    else if (x.n) t.push('<span class="tag na-rna xl-near" title="' + esc('Within 3 residues of an RNA-crosslinked residue: ' + xlSets(x.n).join(', ')) + '">near RNA xl</span>');
    if (x.f) t.push('<span class="tag na-rna xl-near" title="Inside an RNA-bound fragment (RBDmap, HeLa; Castello et al. 2016)">RBDmap</span>');
    if (x.d & 1) t.push('<span class="tag na-dna" title="DNA-crosslinked residue in living MCF7 cells (Trendel et al. 2025)">DNA xl</span>');
    else if (x.d & 2) t.push('<span class="tag na-dna xl-near" title="Inside a DNA-crosslinked peptide (Trendel et al. 2025)">DNA xl peptide</span>');
    if (x.d & 4) t.push('<span class="tag na-dna xl-near" title="Inside a UniProt DNA-binding region">DNA-binding domain</span>');
    return t.join(' ');
  }
  const xlText = (s) => (s.xl ? [s.xl.r ? 'RNA:' + xlSets(s.xl.r).join('/') : '', s.xl.n && !s.xl.r ? 'RNA<=3aa' : '', s.xl.f ? 'RBDmap' : '', s.xl.d & 1 ? 'DNA' : '', s.xl.d & 2 && !(s.xl.d & 1) ? 'DNA peptide' : '', s.xl.d & 4 ? 'UniProt DNA-binding' : ''].filter(Boolean).join(';') : '');
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
    return live('iface:' + base, () => fetch('https://www.ebi.ac.uk/pdbe/graph-api/uniprot/interface_residues/' + encodeURIComponent(base))
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
      }));
  }

  // ----- filters
  function ixProteinPass(p, f) {
    if (!(f.binds.size || f.partner !== '' || f.go.size || f.cpx !== '')) return true;
    if (!ixReady()) return false;
    if (f.go.size) {
      let ok = false;
      f.go.forEach((v) => {
        const bit = v === 1 ? (f.goExp ? IXF.rnaExp : IXF.rna) : (f.goExp ? IXF.dnaExp : IXF.dna);
        if (p.ixFlags & bit) ok = true;
      });
      if (!ok) return false;
    }
    if (f.binds.size || f.partner !== '') {
      const list = partnersAt(p, f.ixScore);
      if (f.partner !== '' && !list.some(([k]) => k === +f.partner)) return false;
      if (f.binds.size && !list.some(([k]) => (f.binds.has(1) && IX.prot[k].flags & IXF.prmt) || (f.binds.has(2) && IX.prot[k].flags & IXF.tudor) || (f.binds.has(3) && k < IX.nP))) return false;
    }
    if (f.cpx !== '' && !p.cx.includes(+f.cpx)) return false;
    return true;
  }

  function ixSitePass(s, f) {
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
    const conf = ' (STRING ≥ ' + (F.ixScore / 1000).toFixed(1) + ')';
    F.binds.forEach((v) => add('ix', 'Binds ' + BIND_OPTS.find((o) => o[0] === v)[1] + conf, () => F.binds.delete(v)));
    if (F.partner !== '' && ixReady()) add('ix', 'Binds ' + ixGene(+F.partner) + conf, () => { F.partner = ''; });
    F.go.forEach((v) => add('ix', (v === 1 ? 'RNA' : 'DNA') + '-binding (GO' + (F.goExp ? ', experimental' : '') + ')', () => F.go.delete(v)));
    F.iface.forEach((v) => add('ix', 'Site at PDB interface with ' + IFACE_NAME[v], () => F.iface.delete(v)));
    if (F.cpx !== '' && ixReady()) add('ix', 'Complex: ' + IX.cx[+F.cpx].name, () => { F.cpx = ''; });
    if (F.modelledPdb) add('ix', 'Modelled as methyl-R in a PDB entry', () => { F.modelledPdb = false; });
    if (F.pred) add('ix', 'Predicted interface (pDockQ ≥ ' + (F.pred / 100).toFixed(2) + ')', () => { F.pred = 0; });
    if (F.galOnly) add('ix', 'Best-supported predicted contact', () => { F.galOnly = false; });
    F.xl.forEach((v) => add('ix', (XL_OPTS.find((o) => o[0] === v) || ['', 'Crosslink'])[1], () => F.xl.delete(v)));
  }

  function ixRailBody(sel) {
    if (!IX) return '<p class="note">Loading interaction data…</p>';
    if (IX.failed) return '<p class="note">Interaction data (<code>data/methylarg_interactions.json</code>) is not available.</p>';
    const partnerOpts = [['', 'Any']].concat(
      [['__p', '— PRMTs —']], IX.prmts.map((x) => [x.k, x.gene]),
      [['__t', '— Tudor-domain proteins —']], IX.tudors.map((x) => [x.k, x.gene]));
    if (F.partner !== '' && !partnerOpts.some((o) => String(o[0]) === String(F.partner))) partnerOpts.push([F.partner, ixGene(+F.partner)]);
    const optHtml = partnerOpts.map(([v, l]) => String(v).startsWith('__') ? '<option disabled>' + esc(l) + '</option>' :
      '<option value="' + esc(v) + '"' + (String(F.partner) === String(v) ? ' selected' : '') + '>' + esc(l) + '</option>').join('');
    return '<div class="fgroup"><span class="flabel" title="STRING v' + esc(IX.meta.string_version) + ' physical subnetwork">Has a physical partner that is (match any)</span>' + chipRow('binds', BIND_OPTS, F.binds) +
      '<label class="flabel" for="f-partner">Binds a specific protein</label>' + sel('f-partner', optHtml, F.partner !== '') +
      '<label class="flabel" for="f-ixscore">STRING physical confidence</label>' + sel('f-ixscore', options(RAIL_SCORES, F.ixScore), false) + '</div>' +
      '<div class="fgroup"><span class="flabel" title="Gene Ontology GO:0003723 (RNA binding) or GO:0003677 (DNA binding), including child terms">Protein binds (Gene Ontology, match any)</span>' + chipRow('go', [[1, 'RNA'], [2, 'DNA']], F.go) +
      '<label class="check"><input type="checkbox" id="f-goexp"' + (F.goExp ? ' checked' : '') + '> Experimental evidence only (not electronic or inferred)</label></div>' +
      '<div class="fgroup"><span class="flabel" title="PDBe-KB: in at least one PDB structure the residue contacts this kind of partner">Site at a PDB interface with (match any)</span>' + chipRow('iface', IFACE_OPTS, F.iface) +
      '<label class="check"><input type="checkbox" id="f-modelled"' + (F.modelledPdb ? ' checked' : '') + '> Modelled as a methylated arginine (2MR, DA2, NMM, AGM) in a PDB entry</label>' +
      '<label class="flabel" for="f-pred" title="ProtVar: residues within 5 Å of the partner in AlphaFold-Multimer models of human protein pairs. pDockQ estimates model quality; 0.23 is the usual floor for an acceptable interface, 0.5 for a confident one.">Site at a predicted interface (ProtVar)</label>' +
      sel('f-pred', options(PRED_OPTS, F.pred), !!F.pred) +
      '<label class="check"><input type="checkbox" id="f-gal"' + (F.galOnly ? ' checked' : '') + '> Best-supported predicted contact (salt bridge or H-bond, confident model, pLDDT ≥ ' + UI.gal.pl + ')</label></div>' +
      '<div class="fgroup"><span class="flabel" title="UV crosslinking mass spectrometry in living human cells: pepR-MS compendium (Kramer 2014, XRNAX, RBS-ID, pRBS-ID, pepR-MS), RBDmap, Trendel 2025 DNA crosslinks; and UniProt DNA-binding regions">RNA/DNA binding in cells (match any)</span>' + chipRow('xl', XL_OPTS, F.xl) + '</div>';
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

  function bindIxRail(rail, refresh) {
    const part = $('#f-partner', rail);
    if (part) part.addEventListener('change', () => { F.partner = part.value; refresh(); });
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
        (x.flags & IXF.rna ? '<br>RNA binding (GO)' : '') + (x.flags & IXF.dna ? '<br>DNA binding (GO)' : '');
      g += '<g class="nd" data-k="' + k + '"' + (forExport ? '' : ' tabindex="0" role="button" data-tip="' + esc(tipHtml) + '"') + ' style="cursor:pointer">';
      if (x.flags & IXF.tudor) g += '<circle cx="' + cxp + '" cy="' + cyp + '" r="' + (r + 4) + '" fill="none" stroke="' + pal.accent + '" stroke-width="1.4" stroke-dasharray="3 2"/>';
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
        (net.ring && net.ring.length ? '<line x1="190" y1="' + (ly2 - 4) + '" x2="214" y2="' + (ly2 - 4) + '" stroke="' + pal.accent + '" stroke-dasharray="4 3"/><text x="220" y="' + ly2 + '" font-size="10.5" font-style="italic" fill="' + pal.ink + '">Second ring: PRMTs and Tudor proteins that bind the partners</text>' : '') +
        '<text x="' + (W - 16) + '" y="' + ly2 + '" font-size="9" text-anchor="end" fill="' + pal.muted + '">Edges: STRING physical links (CC BY 4.0) · methylarginine set: PTM_Pipeline · ' + esc(PAPER.short) + '</text>';
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
    box.innerHTML = head + '<div class="seg ix-tabs" role="tablist">' + tabs.map(([k, l, n]) => '<button type="button" role="tab" data-ixtab="' + k + '" aria-pressed="' + (UI.ixTab === k) + '">' + esc(l) + (n != null ? ' <span class="n">' + fmt(n) + '</span>' : '') + '</button>').join('') + '</div><div id="ix-body"></div>';
    $$('[data-ixtab]', box).forEach((b) => b.addEventListener('click', () => { UI.ixTab = b.dataset.ixtab; renderIxCard(p); }));
    const body = $('#ix-body', box);
    if (UI.ixTab === 'complexes') renderIxComplexes(p, body);
    else if (UI.ixTab === 'rna' || UI.ixTab === 'dna') renderIxNucleic(p, body, UI.ixTab);
    else renderIxProteins(p, body);
  }

  function partnerTags(x) {
    return (x.flags & IXF.prmt ? '<span class="tag ix-prmt">PRMT</span> ' : '') + (x.flags & IXF.tudor ? '<span class="tag ix-tudor">Tudor</span> ' : '') +
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
      '<span><i class="hollow"></i>Not in this set</span><span><i class="diamond"></i>PRMT</span><span><i class="ring"></i>Tudor domain</span></div>' +
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
              '<td class="hide-sm">' + (predMap.has(k) ? (predMap.get(k).model ? '<button class="linkish" type="button" data-pvk="' + k + '" title="Open the ProtVar model">pDockQ ' + (predMap.get(k).sc / 100).toFixed(2) + '</button>' : 'pDockQ ' + (predMap.get(k).sc / 100).toFixed(2)) + ' <span class="muted">' + predMap.get(k).sites.slice(0, 4).map((x) => 'R' + x).join(', ') + (predMap.get(k).sites.length > 4 ? ' +' + (predMap.get(k).sites.length - 4) : '') + '</span>' : '<span class="muted">—</span>') + '</td>' +
              '<td class="hide-sm">' + (shared.length ? shared.slice(0, 2).map((c) => esc(IX.cx[c].name)).join('; ') + (shared.length > 2 ? ' +' + (shared.length - 2) : '') : '<span class="muted">—</span>') + '</td></tr>';
          }).join('') + '</tbody></table></div>' + (all.length > 40 ? '<button class="linkish" type="button" id="ix-all">' + (UI.ixAll ? 'Show the first 40' : 'Show all ' + fmt(all.length)) + '</button>' : '') : '') +
        '<p class="note">STRING v' + esc(IX.meta.string_version) + ' physical subnetwork (experiments, curated databases and text mining of physical binding). RNA*/DNA* mark GO annotations without experimental evidence. Click a methylated partner to open it.</p>';
      const lst = $('#ix-list', box);
      if (lst) lst.addEventListener('click', () => setFilters((f) => { f.partner = String(p.i); f.ixScore = Math.max(700, min); }, 'proteins'));
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
      $('#ix-svg').onclick = () => { if (UI.ixNet) downloadBlob(name + '.svg', new Blob([ixNetSvg(UI.ixNet, EXPORT_PAL, true)], { type: 'image/svg+xml;charset=utf-8' })); };
      $('#ix-png').onclick = () => { if (!UI.ixNet) return; const svg = ixNetSvg(UI.ixNet, EXPORT_PAL, true); const m = /width="(\d+)" height="(\d+)"/.exec(svg); svgToPng(svg, +m[1], +m[2], 3, name + '.png'); };
      $('#ix-tsv').onclick = () => download(name + '.tsv', ['partner_gene\tpartner_accession\tstring_physical_score\tmethylarginine_sites\tprmt\ttudor_domain\tgo_rna_binding\tgo_dna_binding\tprotvar_pdockq\tmethylarginines_at_predicted_interface']
        .concat(all.map(([k, sc]) => { const x = IX.prot[k]; const r = predMap.get(k); return [x.gene, x.acc, (sc / 1000).toFixed(3), x.p >= 0 ? P[x.p].sites.length : 0, x.flags & IXF.prmt ? 1 : 0, x.flags & IXF.tudor ? 1 : 0, x.flags & IXF.rna ? (x.flags & IXF.rnaExp ? 'experimental' : 'yes') : '', x.flags & IXF.dna ? (x.flags & IXF.dnaExp ? 'experimental' : 'yes') : '', r ? (r.sc / 100).toFixed(2) : '', r ? r.sites.map((v) => 'R' + v).join(';') : ''].join('\t'); })).join('\n') + '\n');
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
    $$('[data-cx]', body).forEach((b) => b.addEventListener('click', () => setFilters((f) => { f.cpx = b.dataset.cx; }, 'proteins')));
  }

  // Crosslinking evidence for one protein: crosslinked residues, methylarginines at or near them, regions.
  function xlSectionHtml(p, kind) {
    const xp = IX.xlProt.get(p.i);
    const sites = p.sites.map((j) => S[j]);
    const list = (arr) => arr.map((s) => '<button class="linkish strong" type="button" data-pos="' + s.pos + '">R' + s.pos + '</button>').join(', ');
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
    const at = sites.filter((s) => s.xl && s.xl.d & 1), pep = sites.filter((s) => s.xl && !(s.xl.d & 1) && s.xl.d & 2), bind = sites.filter((s) => s.xl && s.xl.d & 4);
    return '<section><h4>DNA contact in living cells</h4>' +
      '<p>' + (xp && xp.dnaInt ? '<b>In the DNA interactome</b> of MCF7 cells (direct protein–DNA photo-crosslinking; Trendel et al. 2025).' : 'Not in the Trendel et al. 2025 DNA interactome (MCF7).') +
      (xp && xp.dnaPos.length ? ' ' + fmt(xp.dnaPos.length) + ' DNA-crosslinked residues.' : '') + (xp && xp.dnaBind.length ? ' UniProt DNA-binding regions: ' + ranges(xp.dnaBind) + '.' : '') + '</p>' +
      '<ul class="ix-sites">' + (at.length ? '<li><b>At a DNA-crosslinked residue:</b> ' + list(at) + '</li>' : '') + (pep.length ? '<li><b>Inside a DNA-crosslinked peptide:</b> ' + list(pep) + '</li>' : '') + (bind.length ? '<li><b>Inside a UniProt DNA-binding region:</b> ' + list(bind) + '</li>' : '') +
      (!at.length && !pep.length && !bind.length ? '<li class="muted">No methylarginine at a DNA-crosslinked residue or in a DNA-binding region.</li>' : '') + '</ul></section>';
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
        '<ul class="ix-sites">' + sites.map((s) => { const g = contactGroups(s).find((x) => x.k === (kind === 'rna' ? -1 : -2)); return '<li><button class="linkish strong" type="button" data-pos="' + s.pos + '">R' + s.pos + '</button> ' + pill(s.st) +
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
  function ixOverviewHtml() {
    if (!IX) return '<div class="card wide" id="ch-ix"><h3>Interactions</h3><p class="note">Loading interaction data…</p></div>';
    if (IX.failed) return '';
    const e = IX.enrich, g = IX.goEnrich;
    const min = UI.cxMin, curated = !UI.cxPred;
    const rows = IX.cx.filter((c) => c.pred === !curated && c.n >= min && c.m > 0)
      .sort(UI.cxSort === 'share' ? (a, b) => b.m / b.n - a.m / a.n || b.n - a.n : (a, b) => (a.q - b.q) || (b.m / b.n - a.m / a.n));
    const shown = rows.slice(0, UI.cxAll ? rows.length : 15);
    // Methylated proteins bound by each PRMT / Tudor protein (STRING physical, high confidence).
    const boundBy = new Map();
    for (const p of P) for (const [k, sc] of p.partners || []) if (sc >= 700 && IX.prot[k].flags & (IXF.prmt | IXF.tudor)) boundBy.set(k, (boundBy.get(k) || 0) + 1);
    const readerRows = (list) => list.map((x) => ({ x, n: boundBy.get(x.k) || 0 })).filter((r) => r.n).sort((a, b) => b.n - a.n);
    const prmtRows = readerRows(IX.prmts), tudorRows = readerRows(IX.tudors).slice(0, 12);
    const maxN = Math.max(1, ...prmtRows.map((r) => r.n), ...tudorRows.map((r) => r.n));
    const barList = (list) => '<div class="hbars">' + list.map(({ x, n }) => '<button class="hb" type="button" data-partner="' + x.k + '"><span class="hb-label">' + esc(x.gene) + (x.p >= 0 ? ' <span class="muted">· methylated</span>' : '') + '</span><span class="hb-track"><span class="hb-bar" style="left:0;width:' + (100 * n / maxN).toFixed(1) + '%"></span></span><span class="hb-val">' + fmt(n) + '</span></button>').join('') + '</div>';
    let edgesMeth = 0, edgesAll = 0;
    for (const p of P) for (const [k, sc] of p.partners || []) if (sc >= 700) { edgesAll++; if (k < IX.nP) edgesMeth++; }
    const ifRow = (key, label, chip) => {
      const r = e[key];
      if (!r) return '';
      const [a, b, c, d, l2, pv] = r;
      return '<tr class="row" tabindex="0" data-iface="' + chip + '"><td><b>' + esc(label) + '</b></td><td class="r num">' + fmt(a) + ' / ' + fmt(a + b) + ' (' + pctText(a, a + b || 1) + '%)</td><td class="r num">' + fmt(c) + ' / ' + fmt(c + d) + ' (' + pctText(c, c + d || 1) + '%)</td><td class="r num">' + fmtL2(l2) + '</td><td class="r num">' + fmtP(pv) + '</td></tr>';
    };
    const goRow = (key, label, chip, exp) => {
      const r = g[key];
      if (!r) return '';
      const [a, n1, c, n2, l2, pv] = r;
      return '<tr class="row" tabindex="0" data-go="' + chip + '" data-goexp="' + (exp ? 1 : 0) + '"><td><b>' + esc(label) + '</b></td><td class="r num">' + fmt(a) + ' / ' + fmt(n1) + ' (' + pctText(a, n1 || 1) + '%)</td><td class="r num">' + fmt(c) + ' / ' + fmt(n2) + ' (' + pctText(c, n2 || 1) + '%)</td><td class="r num">' + fmtL2(l2) + '</td><td class="r num">' + fmtP(pv) + '</td></tr>';
    };
    return '<div class="card wide" id="ch-cx"><div class="card-head"><div><h3>Protein complexes</h3><p class="sub">Complex Portal complexes ranked by their arginine-methylated members (hypergeometric test against the ' + fmt(IX.meta.proteome_size) + ' reviewed human proteins, ' + pctText(IX.meta.methylated_in_proteome, IX.meta.proteome_size) + '% of which are in this set; BH-adjusted). Click a complex to list its methylated members.</p></div>' +
      '<div class="exp"><button type="button" id="cx-tsv">TSV</button></div></div>' +
      '<div class="fig-controls"><div class="seg" role="group" aria-label="Complex set"><span>Complexes</span><button type="button" data-cxset="cur" aria-pressed="' + curated + '">Curated</button><button type="button" data-cxset="pred" aria-pressed="' + !curated + '">Predicted (hu.MAP 3.0)</button></div>' +
      '<label class="fig-ctl">Min. members <select id="cx-min">' + [2, 3, 5, 10, 20].map((n) => '<option' + (n === min ? ' selected' : '') + '>' + n + '</option>').join('') + '</select></label>' +
      '<label class="fig-ctl">Rank by <select id="cx-sort"><option value="q"' + (UI.cxSort !== 'share' ? ' selected' : '') + '>Enrichment q</option><option value="share"' + (UI.cxSort === 'share' ? ' selected' : '') + '>Share methylated</option></select></label></div>' +
      '<div class="tbl-wrap"><table class="tbl compact cx-tbl"><thead><tr><th>Complex</th><th class="r">Members</th><th>Methylated</th><th class="r hide-xs">Fold</th><th class="r">q</th></tr></thead><tbody>' +
      shown.map((c) => '<tr class="row" tabindex="0" data-cx="' + c.k + '"><td><b>' + esc(c.name) + '</b> <span class="muted mono">' + esc(c.id) + '</span></td><td class="r num">' + c.n + '</td>' +
        '<td><div class="meter-row"><span class="meter"><i style="width:' + (100 * c.m / c.n).toFixed(1) + '%"></i></span><span class="num">' + c.m + ' (' + Math.round(100 * c.m / c.n) + '%)</span></div></td>' +
        '<td class="r num hide-xs">' + (c.fold != null ? c.fold.toFixed(1) + '×' : '—') + '</td><td class="r num">' + fmtP(c.q) + '</td></tr>').join('') +
      '</tbody></table></div>' + (rows.length > 15 ? '<button class="linkish" type="button" id="cx-more" style="margin-top:.5rem">' + (UI.cxAll ? 'Show the top 15' : 'Show all ' + fmt(rows.length) + ' complexes') + '</button>' : '') + '</div>' +
      '<div class="card chart-card" id="ch-readers"><h3>Writers and readers as partners</h3><p class="sub">Methylated proteins with a high-confidence (≥ 0.7) STRING physical link to each PRMT or Tudor-domain protein. Click a bar to list them.</p>' +
      '<p class="flabel">Protein arginine methyltransferases</p>' + barList(prmtRows) + '<p class="flabel" style="margin-top:.8rem">Tudor-domain proteins</p>' + barList(tudorRows) +
      '<p class="note">Summed over methylated proteins, ' + fmt(edgesMeth) + ' of their ' + fmt(edgesAll) + ' high-confidence physical partners (' + pctText(edgesMeth, edgesAll || 1) + '%) are themselves methylated, against ' + pctText(IX.meta.methylated_in_proteome, IX.meta.proteome_size) + '% of reviewed human proteins. Well-studied proteins carry both more mapped sites and more mapped interactions, so read this as association, not mechanism.</p></div>' +
      '<div class="card wide chart-card" id="ch-iface"><h3>Interfaces and nucleic-acid binding</h3><p class="sub">Do methylated arginines sit at protein, RNA or DNA interfaces in PDB structures (PDBe-KB) more often than the other arginines of the same proteins? Only structurally observed arginines (SIFTS) count, and most methylarginines sit in disordered regions that are rarely crystallised.</p>' +
      strataTable() + xlEnrichHtml() + geometryHtml() +
      '<details class="dataset-ctx"><summary class="flabel">Pooled across proteins (not controlled; click a row to filter sites)</summary><div class="tbl-wrap"><table class="tbl compact"><thead><tr><th>Arginines at an interface with</th><th class="r">Methylated R</th><th class="r">Other R</th><th class="r">log2 OR</th><th class="r">p</th></tr></thead><tbody>' +
      ifRow('rna', 'RNA', 4) + ifRow('dna', 'DNA', 8) + ifRow('protein', 'Protein', 1) + '</tbody></table></div>' +
      '<p class="note">Pooling mixes proteins: proteins solved bound to RNA or DNA also tend to carry more mapped methylarginines, which inflates these odds ratios. The within-protein test above removes that.</p></details>' +
      '<h4 style="margin-top:1rem">Gene Ontology binding (whole proteins)</h4><div class="tbl-wrap"><table class="tbl compact"><thead><tr><th>Gene Ontology</th><th class="r">Methylated proteins</th><th class="r">Other proteins</th><th class="r">log2 OR</th><th class="r">p</th></tr></thead><tbody>' +
      goRow('rna', 'RNA binding', 1, false) + goRow('rna_exp', 'RNA binding, experimental', 1, true) + goRow('dna', 'DNA binding', 2, false) + goRow('dna_exp', 'DNA binding, experimental', 2, true) + '</tbody></table></div>' +
      '<p class="note">Methylated vs other reviewed human proteins (Fisher exact test). Click a row to list the proteins.</p></div>' +
      galleryHtml();
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
    html += predEnrichHtml();
    if (IX.nModelled) html += '<p class="note"><button class="linkish" type="button" id="geo-modelled">' + fmt(IX.nModelled) + ' union sites are modelled as a methylated arginine in at least one PDB entry</button>; PDB entries model a further ' + fmt(Array.from(IX.modelledGaps.values()).reduce((n, g) => n + g.pos.length, 0)) + ' methylated positions that the union lacks (listed on each protein page).</p>';
    return html;
  }

  // Crosslinking MS in living cells: methylated vs other arginines of the same proteins.
  function xlEnrichHtml() {
    const e = IX.xlEnrich || {};
    if (!e.rna_xl) return '';
    const row = (key, label, filt) => {
      const x = e[key];
      if (!x) return '';
      const [a, b, c, d] = x.pooled;
      return '<tr' + (filt ? ' class="row" tabindex="0" data-xl="' + filt + '"' : '') + '><td><b>' + label + '</b></td><td class="r num">' + fmt(a) + ' / ' + fmt(a + b) + ' <span class="muted">(' + pctText(a, a + b || 1) + '%)</span></td>' +
        '<td class="r num">' + pctText(c, c + d || 1) + '%</td><td class="r num">' + (x.mh[0] == null ? '—' : '<b>' + fmtL2(x.mh[0]) + '</b> <span class="muted">p ' + fmtP(x.mh[1]) + '</span>') + '</td><td class="r num hide-xs">' + fmt(x.mh[2]) + ' / ' + fmt(x.proteins) + '</td></tr>';
    };
    return '<h4 style="margin-top:1rem">RNA and DNA contact in living cells (crosslinking mass spectrometry)</h4>' +
      '<p class="sub">UV crosslinking locates the residue that touches RNA or DNA in human cells, including in disordered regions that structures miss. Each row compares methylated arginines with the other arginines of proteins that have that evidence (Mantel–Haenszel, each protein a stratum).</p>' +
      '<div class="tbl-wrap"><table class="tbl compact strata-tbl"><thead><tr><th>Arginine is…</th><th class="r">Methylated R</th><th class="r">Other R</th><th class="r">Within protein (log2 MH OR)</th><th class="r hide-xs">Informative / proteins</th></tr></thead><tbody>' +
      row('rna_xl', 'An RNA-crosslinked residue', 1) + row('rna_xl_nonrg', '… excluding RG/GR dipeptides') + row('rna_near', 'Within 3 aa of an RNA crosslink', 2) + row('rna_rbdmap', 'Inside an RBDmap RNA-bound fragment', 4) + row('rna_rbdmap_nonrg', '… excluding RG/GR dipeptides') +
      row('dna_xl', 'A DNA-crosslinked residue', 8) + row('dna_pep', 'Inside a DNA-crosslinked peptide') + row('dna_bind', 'Inside a UniProt DNA-binding region', 16) + '</tbody></table></div>' +
      '<p class="note">RNA: pepR-MS compendium (Sha, Kuster & Trendel 2025; with Kramer 2014, XRNAX 2019, RBS-ID 2020 and pRBS-ID 2021 re-analysed) and RBDmap (Castello et al. 2016). DNA: Trendel et al. 2025 (MCF7). None of the pepR-MS crosslinked peptides carries a methylation mass shift, so the RNA counts see each site through its unmethylated fraction. Click a row to list the sites.</p>';
  }

  // CD-CODE scores cluster at 1-2 and 4-5, so a threshold of 3 or 4 keeps nearly the same members.
  function condHistNote() {
    const h = (IX && IX.condHist) || {};
    const n = (k) => h[String(k)] || 0;
    const tot = Object.values(h).reduce((a, b) => a + b, 0);
    if (!tot) return '';
    return 'Of ' + fmt(tot) + ' human memberships, ' + fmt(n(0) + n(1) + n(2)) + ' score 0–2, ' + fmt(n(3)) + ' score 3 and ' + fmt(n(4) + n(5)) + ' score 4–5, so from 3 upward the chart rests on the small high-confidence set. ';
  }

  // ----- predicted-contact gallery
  const GAL_DEFAULT = { conf: 'confident', pl: 80, kind: 'bond', noPdb: false };
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
    return true;
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
    return '<div class="card wide" id="ch-gal"><div class="card-head"><div><h3>Best-supported predicted contacts</h3><p class="sub">Methylarginines that a predicted complex places in a specific contact with the partner: a salt bridge (guanidino N to Asp/Glu carboxylate O ≤ 4 Å) or hydrogen bond (guanidino N to N/O ≤ 3.5 Å), in a confident model, with both residues at pLDDT ≥ the floor and no steric clash. Models: ProtVar served AlphaFold-Multimer models with pDockQ ≥ 0.5, and AlphaFold DB human heterodimers (ipTM ≥ 0.8 when "confident"). Rows where two models agree come first.</p></div>' +
      '<div class="exp"><button type="button" id="gal-tsv">TSV</button></div></div>' +
      '<div class="fig-controls">' + sel('gal-kind', 'Contact', [['bond', 'Salt bridge or H-bond'], ['sb', 'Salt bridge'], ['pi', 'Cation–π (≤ 5 Å)'], ['any', 'Any atom ≤ 4 Å']], o.kind) +
      sel('gal-pl', 'pLDDT ≥', [[70, '70'], [80, '80'], [90, '90']], o.pl) + sel('gal-conf', 'Model', [['confident', 'Confident (pDockQ ≥ 0.5 / ipTM ≥ 0.8)'], ['acceptable', 'Also ipTM 0.6–0.8']], o.conf) +
      '<label class="check fig-ctl"><input type="checkbox" id="gal-nopdb"' + (o.noPdb ? ' checked' : '') + '> Only sites without a PDB interface</label></div>' +
      '<p class="note"><b>' + fmt(rows.length) + '</b> methylarginine–partner contacts on <b>' + fmt(nSites) + '</b> sites pass. <button class="linkish" type="button" id="gal-sites">List these sites</button>. ' + fmt(clashes) + ' measured contacts with a steric clash (< 2.2 Å) are left out.</p>' +
      (rows.length ? '<div class="tbl-wrap"><table class="tbl compact gal-tbl"><thead><tr><th>Site</th><th>Partner</th><th>Contact in the model</th><th>Models</th><th class="r">pLDDT R / partner</th><th class="r hide-xs">Closest</th><th class="hide-sm">PDB</th></tr></thead><tbody>' +
        shown.map((r, ri) => {
          const s = S[r.j], p = P[s.p], g = r.models[0], x = IX.prot[r.k];
          return '<tr><td class="nw"><button class="linkish strong" type="button" data-galsite="' + ri + '">' + esc(p.label) + ' R' + s.pos + '</button></td>' +
            '<td>' + (x.p >= 0 ? '<button class="linkish" type="button" data-open="' + x.p + '">' + esc(x.gene) + '</button>' : '<b>' + esc(x.gene) + '</b>') + ' ' + partnerTags(x) + '</td>' +
            '<td>' + bonds(g) + '</td>' +
            '<td>' + r.models.map((m, mi) => '<button class="tag pred' + (m.src === 0 ? '' : ' pred-hi') + '" type="button" data-galmodel="' + ri + ':' + mi + '" title="Open this model at R' + s.pos + '">' + esc(galModelText(m)) + '</button>').join(' ') + '</td>' +
            '<td class="r num">' + g.plR + ' / ' + g.plP + '</td><td class="r num hide-xs">' + g.d.toFixed(1) + ' Å</td>' +
            '<td class="hide-sm">' + (s.iface ? '<span class="tag iface-p" title="This site is at an interface in a PDB structure">PDB</span>' : '<span class="muted">—</span>') + '</td></tr>';
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
    const more = $('#gal-more', card);
    if (more) more.addEventListener('click', () => { o.all = !o.all; redo(); });
    $('#gal-sites', card).addEventListener('click', () => setFilters((f) => { f.galOnly = true; }, 'sites'));
    $$('[data-open]', card).forEach((b) => b.addEventListener('click', () => openProtein(+b.dataset.open)));
    $$('[data-galsite]', card).forEach((b) => b.addEventListener('click', () => { const r = shown[+b.dataset.galsite]; openGalleryModel(r.models[0]); }));
    $$('[data-galmodel]', card).forEach((b) => b.addEventListener('click', () => { const [ri, mi] = b.dataset.galmodel.split(':').map(Number); openGalleryModel(shown[ri].models[mi]); }));
    $('#gal-tsv', card).addEventListener('click', () => download('methylarginine_predicted_contacts.tsv', ['gene\taccession\tsite\tpartner_gene\tpartner_accession\tsource\tscore\tpdockq\tplddt_arginine\tplddt_partner_min\tmin_distance_A\tsalt_bridge\th_bond\tcation_pi\tpartner_residues\tmodel']
      .concat(rows.flatMap((r) => r.models.map((g) => { const s = S[g.j], p = P[s.p], x = IX.prot[g.k];
        return [p.label, p.acc, 'R' + s.pos, x.gene, x.acc, g.src === 0 ? 'ProtVar' : 'AlphaFold DB', g.src === 0 ? 'pDockQ ' + (g.score / 100).toFixed(2) : 'ipTM ' + (g.score / 100).toFixed(2), g.pdq == null ? '' : (g.pdq / 100).toFixed(2), g.plR, g.plP, g.d.toFixed(1),
          g.flags & 1 ? 1 : 0, g.flags & 2 ? 1 : 0, g.flags & 4 ? 1 : 0, g.res.join(';'), g.src === 0 ? PROTVAR_MODEL.apply(null, g.ref === 1 ? [p.acc.split('-')[0], x.acc] : [x.acc, p.acc.split('-')[0]]) : 'https://alphafold.ebi.ac.uk/entry/' + g.ref].join('\t'); }))).join('\n') + '\n'));
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
    $$('[data-cxset]', el).forEach((b) => b.addEventListener('click', () => { UI.cxPred = b.dataset.cxset === 'pred'; UI.cxAll = false; renderOverview(el); }));
    const mn = $('#cx-min', el);
    if (mn) mn.addEventListener('change', () => { UI.cxMin = +mn.value; renderOverview(el); });
    const so = $('#cx-sort', el);
    if (so) so.addEventListener('change', () => { UI.cxSort = so.value; renderOverview(el); });
    const more = $('#cx-more', el);
    if (more) more.addEventListener('click', () => { UI.cxAll = !UI.cxAll; renderOverview(el); });
    bindRows($('#ch-cx', el), (tr) => setFilters((f) => { f.cpx = tr.dataset.cx; }, 'proteins'));
    $$('#ch-readers [data-partner]', el).forEach((b) => b.addEventListener('click', () => setFilters((f) => { f.partner = b.dataset.partner; f.ixScore = 700; }, 'proteins')));
    $$('#ch-iface tr[data-iface]', el).forEach((tr) => { const run = () => setFilters((f) => { f.iface = new Set([+tr.dataset.iface]); }, 'sites'); tr.addEventListener('click', run); tr.addEventListener('keydown', (e) => { if (e.key === 'Enter') run(); }); });
    const gm = $('#geo-modelled', el);
    if (gm) gm.addEventListener('click', () => setFilters((f) => { f.modelledPdb = true; }, 'sites'));
    const gp = $('#geo-pred', el);
    if (gp) gp.addEventListener('click', () => setFilters((f) => { f.pred = 23; }, 'sites'));
    bindGallery(el);
    $$('#ch-iface tr[data-xl]', el).forEach((tr) => { const run = () => setFilters((f) => { f.xl = new Set([+tr.dataset.xl]); }, 'sites'); tr.addEventListener('click', run); tr.addEventListener('keydown', (e) => { if (e.key === 'Enter') run(); }); });
    $$('#ch-iface tr[data-go]', el).forEach((tr) => { const run = () => setFilters((f) => { f.go = new Set([+tr.dataset.go]); f.goExp = tr.dataset.goexp === '1'; }, 'proteins'); tr.addEventListener('click', run); tr.addEventListener('keydown', (e) => { if (e.key === 'Enter') run(); }); });
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
    return '<p class="flabel">Distance to the partner in this model</p><ul class="contact-list">' + list.map((x) => '<li><button class="linkish strong" type="button" data-pos="' + x.pos + '">R' + x.pos + '</button> ' +
      (x.d == null ? '<span class="muted">more than 8 Å from ' + esc(pred.partnerLabel) + '</span>' :
        '<b>' + x.d.toFixed(1) + ' Å</b> to ' + esc((pred.names[x.partner.chain] || 'chain ' + x.partner.chain) + ' ' + x.partner.resn + x.partner.resi) + (x.d < 2.2 ? ' <span class="tag clash" title="Heavy atoms closer than 2.2 Å: a steric clash in the model, not a real contact">clash</span>' : x.d <= 4 ? ' <span class="tag iface-p">contact</span>' : '')) +
      ' <span class="muted">· pLDDT ' + Math.round(x.b) + '</span>' +
      (table.has(x.pos) && (x.d == null || x.d > 5) ? '<br><span class="muted">ProtVar\'s 2024 interface table puts R' + x.pos + ' within 5 Å; this later model does not.</span>' : '') + '</li>').join('') + '</ul>' +
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
      return '<li><button class="linkish strong" type="button" data-pos="' + c.pos + '">R' + c.pos + '</button> <span class="muted">(' + (mol.pdb.pred ? 'pLDDT ' + Math.round(c.b) : esc(mol.pdb.chain) + c.auth) + ')</span><br>' + kinds + '</li>';
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
    const head = $('#print-head');
    if (head) head.innerHTML = '<b>Methylarginine Site Explorer</b> · Shechter Lab · ' + esc(location.href) + ' · printed ' + new Date().toISOString().slice(0, 10);
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
    return '<b>R' + s.pos + '</b> · ' + STATE_LABEL[s.st] + (s.ip ? ' · ' + bits(s.ip).map((b) => ipLabel(M.maron_ips[b])).join('/') + ' IP' : '') + '<br>' +
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
    el.innerHTML = '<div class="card prose">' +
      '<h3>About this resource</h3>' +
      '<p><b>Human only.</b> Every site, protein and annotation is <i>Homo sapiens</i>: canonical UniProt sequences and data from human cells and tissues. The homologous-structures panel on protein pages is the one place related proteins from other species can appear, and they are labelled with their organism.</p>' +
      '<p>The explorer indexes the integrated human methylarginine union assembled by the <a href="https://github.com/Shechterlab/PTM_Pipeline">PTM_Pipeline</a> for the review below. It runs entirely in your browser. Protein features, variants and sequences (UniProt), predicted models with per-residue confidence and AlphaMissense scores (AlphaFold DB), and experimental structures (RCSB PDB and PDBe) are fetched live when you open a protein.</p>' +
      '<h3>How to cite</h3><p>' + esc(PAPER.cite) + ' <a href="https://doi.org/' + PAPER.doi + '">doi:' + PAPER.doi + '</a></p>' +
      '<p>Please also cite the original studies for any sites you use; each site\'s evidence lists them.</p>' +
      '<h3>Sources in the union</h3><dl>' +
      '<dt>Maron 2021</dt><dd>Table S5 of Maron et al., <a href="https://doi.org/10.1016/j.isci.2021.102971">iScience 24:102971 (2021)</a>, a compilation of human methylarginine sites from published proteomics.</dd>' +
      '<dt>ProMetheusDB</dt><dd>Methyl-proteome table (mmc2) with monomethyl (Rme1) and dimethyl (Rme2) calls: Massignani et al., ' + doiLink('10.1016/j.mcpro.2022.100243', 'Mol Cell Proteomics 21:100243 (2022)') + '.</dd>' +
      '<dt>PTM master</dt><dd>Methylarginine rows of the lab\'s cross-PTM human master table (UniProt and IEDB annotations); the same table supplies the phospho, lysine acyl, ubiquitin/SUMO and lysine methyl sites shown as crosstalk. Those other-PTM sites are not exhaustive: they come from an iPTMnet export (HPRD, literature mining and other curated sources) plus UniProt, so modifications reported elsewhere, for example in PhosphoSitePlus or recent proteomics studies, may be missing. iPTMnet: Huang et al., ' + doiLink('10.1093/nar/gkx1104', 'Nucleic Acids Res 46:D542 (2018)') + '; IEDB: Vita et al., ' + doiLink('10.1093/nar/gky1006', 'Nucleic Acids Res 47:D339 (2019)') + '.</dd></dl>' +
      '<p class="note"><b>Licence notice.</b> ' + fmt((M.studies.find((x) => x[0] === 'PhosphoSitePlus') || [0, 0, 0])[2]) + ' sites in the union are reported by PhosphoSitePlus® (via Maron et al. 2021 Table S5). ' + esc(PSP_NOTICE) + ' Please acknowledge PhosphoSitePlus® and cite Hornbeck et al. 2015 (' + doiLink('10.1093/nar/gku1267') + ') when using these sites, including in downloads.</p>' +
      '<p>Studies and databases behind the sites (sites reported · reported only there):</p><ul>' +
      M.studyOrder.map((k) => '<li><b>' + esc(M.studies[k][0]) + '</b> · ' + fmt(M.studies[k][2]) + ' · ' + fmt(M.studies[k][3]) + (STUDY_DOI[M.studies[k][0]] ? ' · ' + doiLink(STUDY_DOI[M.studies[k][0]]) : '') + '<br><span class="muted">' + esc(M.studies[k][1] && M.studies[k][1] !== 'This study.' ? M.studies[k][1] : STUDY_HELP[M.studies[k][0]] || '') + '</span></li>').join('') + '</ul>' +
      '<h3>How sites were processed</h3><ul>' +
      '<li>Source rows were remapped to canonical UniProt sequences by peptide/sequence anchors and deduplicated by canonical accession and position. The original reported accession and position are kept in each site\'s evidence.</li>' +
      '<li><b>Gene symbols</b> follow current UniProt. Symbols reported by the sources (some outdated, a few belonging to paralogs with near-identical sequence, such as SNRPN on SmB/B′ rows) and UniProt synonyms are kept as searchable aliases.</li>' +
      '<li><b>State</b> is the consensus Rme1/Rme2 call where a source provides one; most compiled sites have no state call. Maron Table S2 Rme1/Rme2 calls and Table S3 methyl-IP recovery are shown separately; the IPs are named by the state their antibody recognises: Rme1 (MMA antibody), Rme2a (ADMA antibody) and Rme2s (SDMA antibody).</li>' +
      '<li><b>PRMT inhibitor response</b> uses Maron Table S3 <code>GSK591_vs_DMSO</code> and <code>MS023_vs_DMSO</code> differences (median over peptides) with adjusted p ≤ 0.05. GSK591 inhibits PRMT5; MS023 inhibits type I PRMTs.</li>' +
      '<li><b>Motifs</b> are read from a ±7 window on the canonical sequence: ' + M.motifs.map((m) => MOTIF_LABEL[m] + ' (' + MOTIF_HELP[m] + ')').join('; ') + '.</li>' +
      '<li><b>Disorder</b> uses MobiDB intervals: ' + Object.keys(DISORDER_LABEL).map((k) => DISORDER_LABEL[k] + ' = ' + DISORDER_HELP[k].toLowerCase()).join('; ') + '. The architecture figure draws UniProt "Disordered" regions as a wavy backbone.</li>' +
      '<li><b>Domain context</b> uses InterPro domains and repeats: ' + Object.keys(DOMAIN_LABEL).map((k) => DOMAIN_LABEL[k] + ' = ' + DOMAIN_HELP[k].toLowerCase()).join('; ') + '. Nearest-domain classes come from the lab domain-class ontology (config/domain_class_ontology.json).</li>' +
      '<li><b>RG/RGG-rich segments</b> in the figure and domain table: ' + esc(RG_HELP.charAt(0).toLowerCase() + RG_HELP.slice(1)) + '.</li>' +
      '<li><b>Functional classes</b> come from the lab functional ontology (config/functional_ontology.json). Subclasses are shown under the broad class that most of their proteins share.</li>' +
      '<li><b>Clustering</b>: nearest methyl-R spacing and maximum sites in 3–50 aa windows per protein; dense windows come from a multiscale rolling-window scan against a same-protein arginine null (1,000 permutations, q ≤ 0.05). Architecture classes describe proteins with ≥3 sites in a 20 aa span.</li>' +
      '<li><b>Arg enrichment</b> is the per-protein shrinkage log2 enrichment of methylated arginines relative to the proteome-wide rate, with a Fisher test q value.</li>' +
      '<li><b>Condensates</b> are CD-CODE memberships; condensate-level enrichment uses condensates with at least five members.</li>' +
      '<li><b>Disease</b> associations list DisGeNET terms significantly enriched among methylarginine proteins (clusterProfiler) that include the gene or one of its synonyms, grouped into the audited disease families. The ALS/MN/FTD option uses the full audited family gene set.</li>' +
      '</ul>' +
      '<h3>Live annotations on protein pages</h3><ul>' +
      '<li><b>UniProt</b>: domains, regions, disordered segments and compositional bias for the map; natural variants and mutagenesis (variants at a methylarginine are highlighted, with their disease notes); and UniProt\'s own modified-residue term for each site (asymmetric, symmetric or monomethylarginine). Methylarginines that UniProt annotates but the union lacks are listed under the site table.</li>' +
      '<li><b>AlphaFold DB</b>: the predicted model coloured by pLDDT (very high > 90, confident 70–90, low 50–70, very low < 50), the pLDDT at each site, and the mean AlphaMissense pathogenicity of the 19 R→X substitutions at each site (likely benign < 0.34, likely pathogenic ≥ 0.564). Low pLDDT usually marks disordered segments, where most methylarginines sit.</li>' +
      '<li><b>Experimental structures</b>: PDBe SIFTS mappings list the PDB entries for the protein; observed-residue ranges decide which methylarginines are actually resolved in the mapped chain. An RCSB search flags entries whose chain for this protein contains a methylated arginine residue (chemical components 2MR, DA2, NMM or AGM); those residues are drawn in orange when the entry is viewed. NMR entries show the first model of the ensemble.</li>' +
      '</ul>' +
      '<h3>Interactions</h3><ul>' +
      '<li><b>Protein partners</b>: STRING v12.5 physical subnetwork (CC BY 4.0). High-confidence links (score ≥ 0.7) are part of the exported data and drive the filters and overview; protein pages fetch medium-confidence links (≥ 0.4) and the links among a protein\'s partners live from the STRING API. Partners are coloured by their own methylation in this set; PRMTs and Tudor-domain proteins (InterPro IPR002999) are marked.</li>' +
      '<li><b>Complexes</b>: Complex Portal curated complexes and hu.MAP 3.0 predicted complexes (CC0). Enrichment of methylated members is a one-sided hypergeometric test against the reviewed human proteome, BH-adjusted separately for curated and predicted complexes.</li>' +
      '<li><b>RNA and DNA binding</b>: Gene Ontology GO:0003723 and GO:0003677 including child terms (QuickGO/GOA, CC BY 4.0); experimental evidence means EXP, IDA, IPI, IMP, IGI, IEP or their high-throughput equivalents. Protein pages add live ENCODE eCLIP and TF ChIP-seq experiments and ChIP-Atlas target genes (CC BY 4.0).</li>' +
      '<li><b>Interfaces</b>: PDBe-KB lists every residue that contacts another protein, RNA or DNA in at least one PDB structure, in UniProt numbering. Sites carry these flags. For the two most recent entries per site and partner type, the residues within 4 Å of the arginine side chain (NE, CZ, NH1, NH2, CD) were read from the structure itself (RCSB ModelServer, SIFTS SEQRES mapping), and partner residues are given in UniProt numbering. The overview compares methylated arginines with the other structurally observed arginines of the same protein (Cochran–Mantel–Haenszel, each protein a stratum), also restricted to ordered regions (outside UniProt/MobiDB-lite disordered regions) and to residues inside UniProt domains. In the 3D view, nucleic acids are coloured and methylarginine side chains within 4 Å of RNA, DNA or another chain are listed.</li>' +
      '<li><b>RNA and DNA contact in living cells</b>: UV crosslinking mass spectrometry locates the residue that touched RNA or DNA. RNA-crosslinked residues come from the pepR-MS supplement (Sha, Kuster & Trendel 2025), which re-analyses Kramer 2014, XRNAX 2019, RBS-ID 2020 and pRBS-ID 2021 alongside pepR-MS with one pipeline; RBDmap adds RNA-bound fragments (Castello et al. 2016). DNA-crosslinked peptides with localised residues and the DNA interactome come from Trendel et al. 2025 (MCF7); UniProt adds curated DNA-binding regions. Every position was checked against the canonical sequence (26,012 matched, 68 dropped). Sites are flagged at a crosslinked residue, within 3 residues, or inside a fragment, peptide or region, and the overview compares them with the other arginines of the same proteins, also without RG/GR dipeptides.</li>' +
      '<li><b>Writer and reader paths</b>: the protein network can add a second ring of PRMTs and Tudor-domain proteins that bind its partners (STRING ≥ 0.7).</li>' +
      '<li><b>Interface geometry</b>: for every methylarginine at a PDB interface and up to 12 non-methylated interface arginines per protein, the 10 Å environment in one entry (the most recent showing the interface) gives relative solvent accessibility alone and in the complex (Shrake–Rupley, Arg maximum 274 Å² from Tien et al. 2013) and the Levy (' + doiLink('10.1016/j.jmb.2010.09.028', '2010') + ') class: core, rim or support. Guanidino nitrogens within 3.5 Å of a nucleobase N/O count as base hydrogen bonds, of a phosphate or sugar O as backbone hydrogen bonds. Both are compared within proteins (Mantel–Haenszel). The RG/RGG table repeats the interface test within each sequence context.</li>' +
      '<li><b>Methylarginines modelled in the PDB</b>: entries containing 2MR, DA2, NMM or AGM (RCSB search), their residues from PDBe, mapped to UniProt through SIFTS; union sites are flagged and positions the union lacks are listed on protein pages.</li>' +
      '<li><b>References</b>: CD-CODE publications for each membership are listed in full (titles from NCBI E-utilities) with a PubMed search for all of them.</li>' +
      '<li><b>Condensates</b>: protein pages show CD-CODE\'s own confidence score (1–5), role and evidence types for that protein in each condensate (CD-CODE API, release ' + esc((IX && IX.meta && IX.meta.cdcode_version) || '2.3') + '). Condensate-wide statistics (methylated share and enrichment) are shown separately and labelled as being about the condensate.</li>' +
      '<li><b>Disease</b>: protein pages show UniProt\'s curated disease involvement and Open Targets Platform association scores for the gene (live; Open Targets data is CC0). DisGeNET terms come from the dataset-level enrichment and are listed separately as context.</li>' +
      '</ul>' +
      '<h3>Data sources, versions and licences</h3>' +
      '<p>Every external source behind the explorer, with the release used and its licence. Please cite the sources you rely on.</p>' +
      DATA_SOURCES.map(([group, rows]) => '<h4>' + esc(group) + '</h4><div class="tbl-wrap"><table class="tbl compact src-tbl"><thead><tr><th>Source</th><th>Used for</th><th>Version</th><th>Licence</th><th>Citation</th></tr></thead><tbody>' +
        rows.map(([name, use, ver, lic, cite, doi]) => '<tr><td><b>' + esc(name) + '</b></td><td>' + esc(use) + '</td><td>' + esc(ver || '—') + '</td><td>' + esc(lic || '—') + '</td><td>' + (cite ? esc(cite) + (doi ? ' · ' + doiLink(doi) : '') : '—') + '</td></tr>').join('') + '</tbody></table></div>').join('') +
      '<p><b>Licence of this site.</b> The explorer\'s data, figures and text are licensed CC BY-NC-SA 4.0 and its code MIT (see LICENSE.md in the repository); share-alike, non-commercial terms are what PhosphoSitePlus, DisGeNET and pepR-MS require of the combined data. Third-party data keep the licences listed above.</p>' +
      '<p class="note">Interaction data exported ' + esc((IX && IX.meta && IX.meta.generated) || '—') + '; site data exported ' + esc(M.generated) + '. "At staging" and "at export" mean the release current on those dates. Software licences cover the code, not the data.</p>' +
      '<h3>Interpretive guardrails</h3><ul>' +
      '<li>The integrated union is the primary analysis universe; source support is provenance, not a confidence tier.</li>' +
      '<li>State labels are useful where a source supports them, but the union is not a clean Rme1 / Rme2a / Rme2s state-resolved dataset.</li>' +
      '<li>Domain adjacency is descriptive: an IDR-length-conditioned null shows that simple domain-edge proximity is largely explained by IDR geometry and length.</li>' +
      '<li>The stronger architectural signal is short-range methylarginine clustering relative to same-protein arginine positions.</li>' +
      '<li>Motif and clustering annotations need sequence-confirmed canonical residues; ' + fmt(M.window_flags.missing) + ' sites lack a usable sequence window and ' + fmt(M.window_flags.center_not_r) + ' have a non-R canonical residue (flagged with !).</li>' +
      '<li>Coverage follows the source compilations. Methylation reported only in low-throughput studies that none of the sources compiled is not in the union.</li>' +
      '</ul>' +
      '<h3>Downloads and figures</h3>' +
      '<p class="note">' + esc(PSP_NOTICE) + ' Site downloads that include PhosphoSitePlus® sites carry these terms; DisGeNET-derived disease terms are CC BY-NC-SA 4.0.</p>' +
      '<p>Use <b>Download TSV</b> on the Proteins or Sites view to export the current selection. Every overview chart exports as SVG, PNG or TSV, each protein map exports as SVG or PNG at the current zoom, colouring and tracks, and the 3D view exports as PNG. The data files are regenerated from pipeline outputs with <code>python scripts/36_export_site_explorer_data.py</code>; see the repository README.</p>' +
      '<p class="note">Data exported ' + esc(M.generated) + ' · ' + fmt(M.n_sites) + ' sites · ' + fmt(M.n_proteins) + ' proteins.</p>' +
      '</div>';
  }

  // ----------------------------------------------------------- export / copy
  function siteTsv(rows) {
    const head = ['gene', 'uniprot_accession', 'protein_name', 'site', 'position', 'state', 'maron_s2_state', 'maron_ip', 'sources', 'source_count',
      'window_15', 'window_31', 'motifs', 'disorder_context', 'disorder_edge_distance', 'domain_context', 'nearest_interpro', 'nearest_domain_class',
      'domain_edge_distance', 'domain_side', 'nearest_methyl_r_distance', 'GSK591_decreased', 'GSK591_increased', 'MS023_decreased', 'MS023_increased',
      'maron_s3_measurements', 'condensates', 'functional_classes', 'reported_by', 'independent_reports', 'references', 'pdb_interface_partners', 'pdb_interface_entries', 'pdb_contacts', 'predicted_interface_partners_pdockq', 'rna_dna_crosslinks'];
    const out = [head.join('\t')];
    for (const s of rows) {
      const p = P[s.p];
      const ipr = s.ipr >= 0 ? M.interpro[s.ipr] : null;
      out.push([
        p.gene, p.acc, p.name, 'R' + s.pos, s.pos, STATE_LABEL[s.st], bits(s.s2).map((b) => ['Rme1', 'Rme2'][b]).join(';'),
        bits(s.ip).map((b) => ipLabel(M.maron_ips[b])).join(';'), bits(s.src).map((b) => M.sources[b]).join(';'), popcount(s.src),
        s.win ? s.win.slice(8, 23) : '', s.win, bits(s.mot).map((b) => M.motifs[b]).join(';'),
        M.disorder_classes[s.dc], s.dEdge == null ? '' : s.dEdge, M.domain_contexts[s.dctx], ipr ? ipr[0] + ' ' + ipr[1] : '',
        M.domain_classes[s.dcls], s.domDist == null ? '' : s.domDist, s.domSide, s.nn == null ? '' : s.nn,
        s.trf & 1 ? 1 : 0, s.trf & 2 ? 1 : 0, s.trf & 4 ? 1 : 0, s.trf & 8 ? 1 : 0,
        s.tr.map(([ip, t, d, q]) => ipLabel(M.maron_ips[ip]) + ':' + M.maron_treatments[t] + ':' + d + ':' + q).join(';'),
        p.cond.map((k) => M.condensates[k][0]).join(';'), bits(p.func).map((b) => M.func_classes[b]).join(';'),
        bits(s.study).map((b) => M.studies[b][0]).join(';'), popcount(s.study),
        s.refs.map((k) => M.refs[k]).join(' | '),
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

  // Rasterise an exported SVG on a white canvas at `scale`× for PNG download.
  function svgToPng(svgText, w, h, scale, name) {
    const img = new Image();
    const url = URL.createObjectURL(new Blob([svgText], { type: 'image/svg+xml;charset=utf-8' }));
    img.onload = () => {
      const c = document.createElement('canvas');
      c.width = Math.round(w * scale); c.height = Math.round(h * scale);
      const ctx = c.getContext('2d');
      ctx.fillStyle = '#FFFFFF'; ctx.fillRect(0, 0, c.width, c.height);
      ctx.drawImage(img, 0, 0, c.width, c.height);
      URL.revokeObjectURL(url);
      c.toBlob((b) => { if (b) downloadBlob(name, b); }, 'image/png');
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

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
  const popcount = (n) => { let c = 0; while (n) { c += n & 1; n >>= 1; } return c; };
  const bits = (mask) => { const out = []; for (let i = 0; mask; i++, mask >>= 1) if (mask & 1) out.push(i); return out; };
  const debounce = (fn, ms) => { let t; return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); }; };
  const store = {
    get(k) { try { return localStorage.getItem(k); } catch (e) { return null; } },
    set(k, v) { try { if (v == null) localStorage.removeItem(k); else localStorage.setItem(k, v); } catch (e) { /* storage unavailable */ } },
  };

  // -------------------------------------------------------------- vocabulary
  const STATE_LABEL = ['No state call', 'Rme1', 'Rme2', 'Rme1 + Rme2'];
  const SOURCE_LABEL = { Maron2021: 'Maron 2021 Table S5', ProMetheusDB: 'ProMetheusDB', human_ptm_master: 'PTM master table' };
  const STUDY_HELP = {
    ProMetheusDB: 'ProMetheusDB methyl-proteome table (mmc2)',
    UniProt: 'UniProt methylarginine annotation, via the cross-PTM master table',
    IEDB: 'Immune Epitope Database record with a methylated arginine, via the cross-PTM master table',
    'Maron 2021 (this study)': 'Sites first reported in Maron et al. 2021 (Table S5 "This study")',
  };
  const PAPER = {
    cite: 'DeAngelo JD, Silverstein AM, Shechter D. Protein arginine methylation rewires molecular interactions. Cell Chem Biol 33(10):1404–1421 (2026).',
    doi: '10.1016/j.chembiol.2026.09.011',
  };
  const SOURCE_HELP = {
    Maron2021: 'Maron et al. 2021 (iScience) Table S5: compiled human methylarginine sites from published proteomics',
    ProMetheusDB: 'ProMetheusDB methyl-proteome (mmc2), with Rme1/Rme2 calls',
    human_ptm_master: 'Cross-PTM human master table (UniProt and IEDB methylarginine annotations)',
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
  const TREAT_LABEL = { GSK591: 'GSK591', MS023: 'MS023' };
  const TREAT_HELP = { GSK591: 'PRMT5 inhibitor', MS023: 'type I PRMT inhibitor' };
  const PTM_LABEL = ['Phospho', 'Lys acyl', 'Ub/SUMO', 'Lys methyl'];
  const PTM_PREFIX = ['p', 'ac', 'ub', 'me'];
  const SIDE_LABEL = { inside: 'inside', N_terminal_to_domain: 'N-terminal to domain', C_terminal_to_domain: 'C-terminal to domain' };
  const PADJ = 0.05;
  const TREAT_OPTIONS = [
    ['', 'Any site'],
    ['tested', 'Measured in Maron Table S3'],
    ['0:down', 'Decreased with GSK591 (PRMT5i)'],
    ['1:down', 'Decreased with MS023 (type I PRMTi)'],
    ['0:up', 'Increased with GSK591'],
    ['1:up', 'Increased with MS023'],
  ];

  // ------------------------------------------------------------------ model
  let M = null;      // meta
  let P = [];        // proteins
  let S = [];        // sites
  let IDX = null;    // lookup maps
  let R = null;      // current filtered result
  let evidencePromise = null;
  const uniprotCache = new Map();

  const DEFAULT_FILTERS = () => ({
    states: new Set(), ips: new Set(), sources: new Set(), s2: new Set(), multi: false, studies: new Set(), minRep: 0,
    treat: '', motifs: new Set(), disorder: new Set(), domctx: new Set(), domclass: '',
    nn: 0, ptmGroup: '', ptmDist: 0,
    func: '', cond: '', anyCond: false, disease: '', minSites: 1, dense: false, arch: '', argEnriched: false,
  });
  let F = DEFAULT_FILTERS();
  let Q = { raw: '', proteins: null, positions: null, note: '' };
  const UI = {
    view: 'overview', detail: null, highlight: null,
    sortP: { key: 'sites', dir: -1 }, sortS: { key: 'gene', dir: 1 },
    pageP: 50, pageS: 100, zoom: null, openEvidence: new Set(),
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
    const fields = 'accession,length,sequence,protein_name,gene_primary,ft_domain,ft_repeat,ft_zn_fing,ft_region,ft_compbias,ft_motif';
    const ctrl = typeof AbortController === 'function' ? new AbortController() : null;
    const timer = ctrl ? setTimeout(() => ctrl.abort(), 15000) : null;
    const promise = fetch('https://rest.uniprot.org/uniprotkb/' + encodeURIComponent(acc) + '.json?fields=' + fields,
      ctrl ? { signal: ctrl.signal } : undefined)
      .then((r) => { if (!r.ok) throw new Error(r.status); return r.json(); })
      .then((d) => {
        const desc = d.proteinDescription || {};
        const nameObj = desc.recommendedName || (desc.submissionNames || [])[0] || {};
        return {
          seq: (d.sequence && d.sequence.value) || '',
          len: (d.sequence && d.sequence.length) || null,
          name: (nameObj.fullName && nameObj.fullName.value) || '',
          gene: (d.genes && d.genes[0] && d.genes[0].geneName && d.genes[0].geneName.value) || '',
          feats: (d.features || []).map((f) => ({
            type: f.type, desc: f.description || '',
            start: f.location && f.location.start && f.location.start.value,
            end: f.location && f.location.end && f.location.end.value,
          })).filter((f) => f.start && f.end),
        };
      })
      .catch(() => null)
      .finally(() => { if (timer) clearTimeout(timer); });
    uniprotCache.set(acc, promise);
    return promise;
  }

  function init(data) {
    M = data.meta;
    const vi = (list, key) => list.indexOf(key);
    M.vi = {
      disordered: vi(M.disorder_classes, 'disordered'),
    };
    M.studyOrder = M.studies.map((st, k) => k).sort((a, b) => M.studies[b][2] - M.studies[a][2]);
    P = data.proteins.map((r, i) => ({
      i, acc: r[0], gene: r[1] || '', name: r[2] || '', uid: r[3] || '', func: r[4], sub: r[5], len: r[6],
      argN: r[7], shrink: r[8], argQ: r[9], iso: !!r[10], cond: r[11], dis: r[12], als: !!r[13],
      dens: r[14], best: r[15], sig: r[16], arch: r[17], ptms: r[18], sites: [],
    }));
    S = data.sites.map((r, j) => ({
      j, p: r[0], pos: r[1], st: r[2], src: r[3], refs: r[4], dc: r[5], dEdge: r[6], dctx: r[7], ipr: r[8],
      dcls: r[9], domDist: r[10], domSide: r[11], win: r[12], mot: r[13], mis: !!r[14], orig: r[15],
      s2: r[16], ip: r[17], tr: r[18], study: r[19] || 0, nn: null, ptmD: null, trf: 0,
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
    IDX = { gene: new Map(), acc: new Map(), uid: new Map() };
    const add = (map, key, i) => { if (!key) return; key = key.toUpperCase(); if (!map.has(key)) map.set(key, []); map.get(key).push(i); };
    for (const p of P) {
      add(IDX.gene, p.gene, p.i);
      add(IDX.acc, p.acc, p.i);
      if (p.acc.includes('-')) add(IDX.acc, p.acc.split('-')[0], p.i);
      add(IDX.uid, p.uid, p.i);
    }

    $('#data-version').textContent = 'data exported ' + M.generated;
    if (window.matchMedia && window.matchMedia('(max-width: 560px)').matches) $('#q').placeholder = 'Gene, accession, site or gene list';
    $('#year').textContent = new Date().getFullYear();
    renderHero();
    renderRail();
    bindGlobal();
    const params = new URLSearchParams(location.search);
    if (params.get('q')) { $('#q').value = params.get('q'); Q = parseQuery(params.get('q')); }
    compute();
    route();
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
  function lookupExact(token) {
    const key = token.toUpperCase();
    return IDX.acc.get(key) || IDX.gene.get(key) || IDX.uid.get(key) || null;
  }

  function parseQuery(text) {
    const raw = text.trim();
    if (!raw) return { raw: '', proteins: null, positions: null, note: '' };
    const tokens = raw.split(/[\s,;]+/).filter(Boolean);
    const proteins = new Set();
    const positions = new Map();
    const missing = [];
    const missingSites = [];
    let last = null;
    const addPos = (hit, pos) => {
      for (const i of hit) {
        if (!positions.has(i)) positions.set(i, new Set());
        positions.get(i).add(pos);
        proteins.add(i);
        if (!P[i].sites.some((j) => S[j].pos === pos)) missingSites.push(P[i].label + ' R' + pos);
      }
    };
    for (const tok of tokens) {
      const hit = lookupExact(tok);
      if (hit) { hit.forEach((i) => proteins.add(i)); last = hit; continue; }
      let m = tok.match(/^(.+?)[:_/-]R?(\d+)$/i);
      if (m) {
        const h = lookupExact(m[1]);
        if (h) { addPos(h, +m[2]); last = h; continue; }
      }
      m = tok.match(/^R?(\d+)$/i);
      if (m && last) { addPos(last, +m[1]); continue; }
      missing.push(tok);
    }
    if (!proteins.size) {
      const t = raw.toLowerCase();
      const hits = P.filter((p) => p.gene.toLowerCase().startsWith(t) || p.acc.toLowerCase().startsWith(t) ||
        (t.length >= 3 && p.name.toLowerCase().includes(t)));
      if (hits.length) {
        return { raw, proteins: new Set(hits.map((p) => p.i)), positions: null,
          note: '<b>' + fmt(hits.length) + '</b> protein' + (hits.length === 1 ? '' : 's') + ' match “' + esc(raw) + '” by gene prefix or protein name.' };
      }
      return { raw, proteins: new Set(), positions: null,
        note: 'No protein in the union matches “' + esc(raw) + '”. Try a HGNC symbol (PRMT5), a UniProt accession (O14744) or part of a protein name.' };
    }
    const notes = [];
    const geneTokens = tokens.length - missing.length;
    if (tokens.length > 1) notes.push('<b>' + fmt(proteins.size) + '</b> protein' + (proteins.size === 1 ? '' : 's') + ' from ' + fmt(geneTokens) + ' matched term' + (geneTokens === 1 ? '' : 's') + '.');
    if (missing.length) notes.push('Not in the methylarginine union: ' + missing.slice(0, 25).map(esc).join(', ') + (missing.length > 25 ? ' and ' + (missing.length - 25) + ' more' : '') + '.');
    if (missingSites.length) notes.push('No methylation reported at ' + missingSites.slice(0, 10).map(esc).join(', ') + '.');
    return { raw, proteins, positions: positions.size ? positions : null, note: notes.join(' ') };
  }

  function suggestions(text) {
    const tok = text.split(/[\s,;]+/).pop() || '';
    if (!tok || tok.length < 1) return [];
    const t = tok.toLowerCase();
    const exact = [], prefix = [], named = [];
    for (const p of P) {
      const g = p.gene.toLowerCase();
      if (g === t || p.acc.toLowerCase() === t) exact.push(p);
      else if (g.startsWith(t) || p.acc.toLowerCase().startsWith(t)) prefix.push(p);
      else if (t.length >= 3 && p.name.toLowerCase().includes(t)) named.push(p);
    }
    const bySites = (a, b) => b.sites.length - a.sites.length;
    return exact.concat(prefix.sort(bySites), named.sort(bySites)).slice(0, 8);
  }

  // --------------------------------------------------------------- filtering
  function sitePass(s) {
    if (F.states.size && !F.states.has(s.st)) return false;
    if (F.ips.size && !bits(s.ip).some((b) => F.ips.has(b))) return false;
    if (F.s2.size && !bits(s.s2).some((b) => F.s2.has(b))) return false;
    if (F.sources.size && !bits(s.src).some((b) => F.sources.has(b))) return false;
    if (F.multi && popcount(s.src) < 2) return false;
    if (F.studies.size && !bits(s.study).some((b) => F.studies.has(b))) return false;
    if (F.minRep && popcount(s.study) < F.minRep) return false;
    if (F.treat) {
      if (F.treat === 'tested') { if (!s.tr.length) return false; }
      else {
        const [t, dir] = F.treat.split(':');
        if (!(s.trf & (1 << (+t * 2 + (dir === 'down' ? 0 : 1))))) return false;
      }
    }
    if (F.motifs.size && !bits(s.mot).some((b) => F.motifs.has(b))) return false;
    if (F.disorder.size && !F.disorder.has(s.dc)) return false;
    if (F.domctx.size && !F.domctx.has(s.dctx)) return false;
    if (F.domclass !== '' && s.dcls !== +F.domclass) return false;
    if (F.nn && !(s.nn != null && s.nn <= F.nn)) return false;
    if (F.ptmDist) {
      const d = F.ptmGroup === '' ? Math.min.apply(null, s.ptmD) : s.ptmD[+F.ptmGroup];
      if (!(d <= F.ptmDist)) return false;
    }
    return true;
  }

  function proteinPass(p) {
    if (F.func !== '' && !(p.func & (1 << +F.func))) return false;
    if (F.cond !== '' && !p.cond.includes(+F.cond)) return false;
    if (F.anyCond && !p.cond.length) return false;
    if (F.disease !== '') {
      if (F.disease === 'als') { if (!p.als) return false; }
      else if (!(p.dis & (1 << +F.disease))) return false;
    }
    if (F.minSites > 1 && p.sites.length < F.minSites) return false;
    if (F.dense && !p.sig.length) return false;
    if (F.arch !== '' && !(p.arch && p.arch[0] === F.arch)) return false;
    if (F.argEnriched && !(p.argQ != null && p.argQ <= PADJ)) return false;
    return true;
  }

  function compute() {
    const okP = new Uint8Array(P.length);
    for (const p of P) okP[p.i] = (!Q.proteins || Q.proteins.has(p.i)) && proteinPass(p) ? 1 : 0;
    const sites = [];
    const counts = new Uint32Array(P.length);
    for (const s of S) {
      if (!okP[s.p]) continue;
      if (Q.positions && Q.positions.has(s.p) && !Q.positions.get(s.p).has(s.pos)) continue;
      if (!sitePass(s)) continue;
      sites.push(s);
      counts[s.p]++;
    }
    const prots = P.filter((p) => counts[p.i] > 0);
    R = { sites, prots, counts, shown: new Set(sites.map((s) => s.j)) };
    $('#count-proteins').textContent = fmt(prots.length);
    $('#count-sites').textContent = fmt(sites.length);
    $('#query-note').innerHTML = Q.note || '';
    renderRailState();
  }

  function activeFilters() {
    const out = [];
    const chip = (label, clear) => out.push({ label, clear });
    F.states.forEach((v) => chip('State: ' + STATE_LABEL[v], () => F.states.delete(v)));
    F.s2.forEach((v) => chip('Table S2: ' + ['me1', 'me2'][v], () => F.s2.delete(v)));
    F.ips.forEach((v) => chip(M.maron_ips[v] + ' IP', () => F.ips.delete(v)));
    F.sources.forEach((v) => chip(SOURCE_LABEL[M.sources[v]], () => F.sources.delete(v)));
    if (F.multi) chip('≥2 source tables', () => { F.multi = false; });
    F.studies.forEach((v) => chip('Reported by ' + M.studies[v][0], () => F.studies.delete(v)));
    if (F.minRep) chip('≥' + F.minRep + ' independent reports', () => { F.minRep = 0; });
    if (F.treat) chip(TREAT_OPTIONS.find((o) => o[0] === F.treat)[1], () => { F.treat = ''; });
    F.motifs.forEach((v) => chip('Motif ' + MOTIF_LABEL[M.motifs[v]], () => F.motifs.delete(v)));
    F.disorder.forEach((v) => chip(DISORDER_LABEL[M.disorder_classes[v]], () => F.disorder.delete(v)));
    F.domctx.forEach((v) => chip(DOMAIN_LABEL[M.domain_contexts[v]], () => F.domctx.delete(v)));
    if (F.domclass !== '') chip('Nearest domain: ' + M.domain_classes[+F.domclass], () => { F.domclass = ''; });
    if (F.nn) chip('Methyl-R within ' + F.nn + ' aa', () => { F.nn = 0; });
    if (F.ptmDist) chip((F.ptmGroup === '' ? 'Any PTM' : PTM_LABEL[+F.ptmGroup]) + ' within ±' + F.ptmDist + ' aa', () => { F.ptmDist = 0; });
    if (F.func !== '') chip(M.func_classes[+F.func], () => { F.func = ''; });
    if (F.cond !== '') chip(M.condensates[+F.cond], () => { F.cond = ''; });
    if (F.anyCond) chip('In a CD-CODE condensate', () => { F.anyCond = false; });
    if (F.disease !== '') chip(F.disease === 'als' ? 'ALS/MN/FTD audited gene set' : M.disease_families[+F.disease], () => { F.disease = ''; });
    if (F.minSites > 1) chip('≥' + F.minSites + ' sites per protein', () => { F.minSites = 1; });
    if (F.dense) chip('Dense methyl-R region', () => { F.dense = false; });
    if (F.arch !== '') chip(F.arch, () => { F.arch = ''; });
    if (F.argEnriched) chip('Arg-methylation enriched', () => { F.argEnriched = false; });
    return out;
  }

  function setFilters(mutator, view) {
    mutator(F);
    UI.pageP = 50; UI.pageS = 100;
    compute();
    renderRail();
    if (view && view !== UI.view) go(view); else renderView();
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
  function chipRow(name, items, selected, dots) {
    return '<div class="chips" data-set="' + name + '">' + items.map(([v, label, help]) =>
      '<button class="chip" type="button" data-v="' + v + '" aria-pressed="' + selected.has(v) + '"' + (help ? ' title="' + esc(help) + '"' : '') + '>' +
      (dots ? '<i class="dot st-' + v + '"></i>' : '') + esc(label) + '</button>').join('') + '</div>';
  }

  function options(list, current, first) {
    return (first ? '<option value="">' + esc(first) + '</option>' : '') +
      list.map(([v, label]) => '<option value="' + esc(v) + '"' + (String(current) === String(v) ? ' selected' : '') + '>' + esc(label) + '</option>').join('');
  }

  function renderRail() {
    const rail = $('#rail');
    const open = new Set($$('details[open]', rail).map((d) => d.dataset.sec));
    const firstRender = !rail.children.length;
    const isOpen = (k, dflt) => (firstRender ? dflt : open.has(k)) ? ' open' : '';
    const condOpts = M.condensates.map((n, k) => [k, n]).filter(([k]) => P.some((p) => p.cond.includes(k)))
      .sort((a, b) => a[1].localeCompare(b[1]));
    const archClasses = Array.from(new Set(P.filter((p) => p.arch).map((p) => p.arch[0]))).sort();
    rail.innerHTML =
      '<div class="rail-head"><h2>Filters</h2><div class="actions">' +
      '<button class="linkish" type="button" id="reset-filters">Reset all</button>' +
      '<button class="btn quiet small rail-toggle" type="button" id="rail-toggle" aria-expanded="' + (rail.dataset.collapsed === 'false') + '">' + (rail.dataset.collapsed === 'false' ? 'Hide' : 'Show') + '</button></div></div>' +

      '<details data-sec="methyl"' + isOpen('methyl', true) + '><summary>Methylation evidence<span class="on" data-on="methyl" hidden></span></summary>' +
      '<div class="fgroup"><span class="flabel" title="Consensus Rme1/Rme2 calls across sources; most compiled sites have no state call">Methylation state</span>' +
      chipRow('states', [[1, 'Rme1'], [2, 'Rme2'], [3, 'Rme1 + Rme2'], [0, 'No call']], F.states, true) + '</div>' +
      '<div class="fgroup"><span class="flabel" title="Maron et al. 2021 Table S3: methyl-specific immunoprecipitation that recovered the peptide">Maron methyl-IP</span>' +
      chipRow('ips', M.maron_ips.map((n, k) => [k, n, n === 'SDMA' ? 'Symmetric dimethyl-Arg antibody IP' : n === 'ADMA' ? 'Asymmetric dimethyl-Arg antibody IP' : 'Monomethyl-Arg antibody IP']), F.ips) + '</div>' +
      '<div class="fgroup"><span class="flabel" title="Maron et al. 2021 Table S2 me1/me2 site calls">Maron Table S2 call</span>' +
      chipRow('s2', [[0, 'me1'], [1, 'me2']], F.s2) + '</div>' +
      '</details>' +

      '<details data-sec="prov"' + isOpen('prov', true) + '><summary>Provenance<span class="on" data-on="prov" hidden></span></summary>' +
      '<div class="fgroup"><span class="flabel">Reported by study or database (match any)</span>' +
      chipRow('studies', M.studyOrder.map((k) => [k, M.studies[k][0] + ' · ' + compact(M.studies[k][2]), M.studies[k][1] || STUDY_HELP[M.studies[k][0]] || '']), F.studies) + '</div>' +
      '<div class="fgroup"><label class="flabel" for="f-minrep">Independent reports per site</label><select id="f-minrep">' +
      options([[2, '≥2 studies or databases'], [3, '≥3'], [5, '≥5']], F.minRep || '', 'Any') + '</select></div>' +
      '<div class="fgroup"><span class="flabel">Pipeline source table</span>' +
      chipRow('sources', M.sources.map((n, k) => [k, SOURCE_LABEL[n], SOURCE_HELP[n]]), F.sources) +
      '<label class="check"><input type="checkbox" id="f-multi"' + (F.multi ? ' checked' : '') + '> In ≥2 source tables at the same position</label></div>' +
      '</details>' +

      '<details data-sec="treat"' + isOpen('treat', true) + '><summary>PRMT inhibitors<span class="on" data-on="treat" hidden></span></summary>' +
      '<div class="fgroup"><label class="flabel" for="f-treat">Maron Table S3 response (adj. p ≤ 0.05)</label>' +
      '<select id="f-treat">' + options(TREAT_OPTIONS, F.treat) + '</select></div></details>' +

      '<details data-sec="motif"' + isOpen('motif', false) + '><summary>Sequence motif<span class="on" data-on="motif" hidden></span></summary>' +
      '<div class="fgroup"><span class="flabel">Around the methyl-R (match any)</span>' +
      chipRow('motifs', M.motifs.map((m, k) => [k, MOTIF_LABEL[m], MOTIF_HELP[m]]), F.motifs) + '</div></details>' +

      '<details data-sec="struct"' + isOpen('struct', false) + '><summary>Disorder &amp; domains<span class="on" data-on="struct" hidden></span></summary>' +
      '<div class="fgroup"><span class="flabel">Disorder (MobiDB)</span>' +
      chipRow('disorder', M.disorder_classes.map((c, k) => [k, DISORDER_LABEL[c] || c, DISORDER_HELP[c]]), F.disorder) + '</div>' +
      '<div class="fgroup"><span class="flabel">Domain context (InterPro)</span>' +
      chipRow('domctx', M.domain_contexts.map((c, k) => [k, DOMAIN_LABEL[c] || c, DOMAIN_HELP[c]]), F.domctx) + '</div>' +
      '<div class="fgroup"><label class="flabel" for="f-domclass">Nearest domain class</label><select id="f-domclass">' +
      options(M.domain_classes.map((c, k) => [k, c]).sort((a, b) => a[1].localeCompare(b[1])), F.domclass, 'Any') + '</select></div></details>' +

      '<details data-sec="cluster"' + isOpen('cluster', false) + '><summary>Clustering &amp; crosstalk<span class="on" data-on="cluster" hidden></span></summary>' +
      '<div class="fgroup"><label class="flabel" for="f-nn">Another methyl-R within</label><select id="f-nn">' +
      options([[2, '2 aa'], [3, '3 aa'], [5, '5 aa'], [10, '10 aa'], [20, '20 aa']], F.nn || '', 'Any spacing') + '</select></div>' +
      '<div class="fgroup"><label class="flabel" for="f-ptm">Other PTM nearby (PTM master table)</label><select id="f-ptm">' +
      options(PTM_LABEL.map((l, k) => [k, l]), F.ptmGroup, 'Any PTM type') + '</select>' +
      '<select id="f-ptmdist" aria-label="PTM distance">' + options([[2, 'within ±2 aa'], [5, 'within ±5 aa'], [10, 'within ±10 aa'], [20, 'within ±20 aa']], F.ptmDist || '', 'Off') + '</select></div>' +
      '<div class="fgroup"><label class="check"><input type="checkbox" id="f-dense"' + (F.dense ? ' checked' : '') + '> Protein has a significant dense methyl-R region (multiscale null, q ≤ 0.05)</label></div>' +
      '<div class="fgroup"><label class="flabel" for="f-arch">Clustered architecture (≥3 sites in 20 aa)</label><select id="f-arch">' +
      options(archClasses.map((a) => [a, a]), F.arch, 'Any') + '</select></div></details>' +

      '<details data-sec="protein"' + isOpen('protein', false) + '><summary>Protein features<span class="on" data-on="protein" hidden></span></summary>' +
      '<div class="fgroup"><label class="flabel" for="f-func">Functional class</label><select id="f-func">' +
      options(M.func_classes.map((c, k) => [k, c]).sort((a, b) => a[1].localeCompare(b[1])), F.func, 'Any') + '</select></div>' +
      '<div class="fgroup"><label class="flabel" for="f-cond">CD-CODE condensate</label><select id="f-cond">' +
      options(condOpts, F.cond, 'Any') + '</select>' +
      '<label class="check"><input type="checkbox" id="f-anycond"' + (F.anyCond ? ' checked' : '') + '> In any CD-CODE condensate</label></div>' +
      '<div class="fgroup"><label class="flabel" for="f-disease">Disease family (DisGeNET)</label><select id="f-disease">' +
      options(M.disease_families.map((c, k) => [k, c]).concat([['als', 'ALS / MN / FTD (full audited gene set)']]), F.disease, 'Any') + '</select></div>' +
      '<div class="fgroup"><label class="flabel" for="f-min">Minimum methyl-R sites per protein</label>' +
      '<input id="f-min" type="number" min="1" max="100" step="1" value="' + F.minSites + '"></div>' +
      '<div class="fgroup"><label class="check"><input type="checkbox" id="f-arg"' + (F.argEnriched ? ' checked' : '') + '> Arginines methylated more than the proteome rate (Fisher q ≤ 0.05)</label></div>' +
      '</details>';
    bindRail();
    renderRailState();
  }

  function renderRailState() {
    const on = {
      methyl: F.states.size || F.ips.size || F.s2.size,
      prov: F.studies.size || F.minRep || F.sources.size || F.multi,
      treat: !!F.treat,
      motif: F.motifs.size,
      struct: F.disorder.size || F.domctx.size || F.domclass !== '',
      cluster: F.nn || F.ptmDist || F.dense || F.arch !== '',
      protein: F.func !== '' || F.cond !== '' || F.anyCond || F.disease !== '' || F.minSites > 1 || F.argEnriched,
    };
    $$('[data-on]').forEach((el) => { el.hidden = !on[el.dataset.on]; });
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
    const bind = (id, fn, ev = 'change') => { const el = $('#' + id, rail); if (el) el.addEventListener(ev, () => { fn(el); UI.pageP = 50; UI.pageS = 100; compute(); renderView(); }); };
    bind('f-multi', (el) => { F.multi = el.checked; });
    bind('f-minrep', (el) => { F.minRep = +el.value || 0; });
    bind('f-treat', (el) => { F.treat = el.value; });
    bind('f-domclass', (el) => { F.domclass = el.value; });
    bind('f-nn', (el) => { F.nn = +el.value || 0; });
    bind('f-ptm', (el) => { F.ptmGroup = el.value; if (!F.ptmDist) { F.ptmDist = 5; $('#f-ptmdist').value = '5'; } });
    bind('f-ptmdist', (el) => { F.ptmDist = +el.value || 0; });
    bind('f-dense', (el) => { F.dense = el.checked; });
    bind('f-arch', (el) => { F.arch = el.value; });
    bind('f-func', (el) => { F.func = el.value; });
    bind('f-cond', (el) => { F.cond = el.value; });
    bind('f-anycond', (el) => { F.anyCond = el.checked; });
    bind('f-disease', (el) => { F.disease = el.value; });
    bind('f-min', (el) => { F.minSites = Math.max(1, Math.min(100, +el.value || 1)); }, 'input');
    bind('f-arg', (el) => { F.argEnriched = el.checked; });
    $('#reset-filters').addEventListener('click', () => { F = DEFAULT_FILTERS(); renderRail(); compute(); renderView(); });
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

  function bindGlobal() {
    const input = $('#q');
    const list = $('#suggest');
    let sel = -1;
    let items = [];
    const closeList = () => { list.hidden = true; input.setAttribute('aria-expanded', 'false'); sel = -1; };
    const drawList = () => {
      items = suggestions(input.value);
      if (!items.length || input.value.split(/[\s,;]+/).filter(Boolean).length > 1) { closeList(); return; }
      list.innerHTML = items.map((p, k) => '<li role="option" id="sg-' + k + '" aria-selected="' + (k === sel) + '" data-i="' + p.i + '"><b>' + esc(p.label) + '</b><span>' + esc(p.name || p.acc) + '</span><em>' + p.sites.length + ' site' + (p.sites.length === 1 ? '' : 's') + '</em></li>').join('');
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
      else if (e.key === 'Escape') closeList();
      else if (e.key === 'Enter') {
        e.preventDefault();
        if (sel >= 0 && items[sel]) { pick(items[sel]); return; }
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
    $$('[data-theme-set]').forEach((b) => b.addEventListener('click', () => setTheme(b.dataset.themeSet)));
    setTheme(store.get('mre-theme') || 'system', true);
    const syncToolbar = () => document.documentElement.style.setProperty('--toolbar-h', $('#toolbar').offsetHeight + 'px');
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
  const ipTags = (mask) => bits(mask).map((b) => '<span class="tag ip">' + M.maron_ips[b] + '</span>').join(' ');
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
    const af = activeFilters();
    return '<div class="results-head"><div><h2>' + title + '</h2><p>' + sub + '</p></div><div class="actions">' + (actions || '') + '</div></div>' +
      '<div class="active-filters" id="active-filters">' + af.map((f, k) => '<span class="af">' + esc(f.label) + '<button type="button" data-k="' + k + '" aria-label="Remove filter ' + esc(f.label) + '">×</button></span>').join('') + '</div>';
  }

  function bindActiveFilters(root) {
    const af = activeFilters();
    $$('#active-filters button', root).forEach((b) => b.addEventListener('click', () => {
      af[+b.dataset.k].clear();
      renderRail(); compute(); renderView();
    }));
  }

  function emptyState() {
    return '<div class="card empty"><b>No sites match</b>Remove a filter or broaden the search. ' +
      (activeFilters().length ? '<button class="linkish" type="button" id="empty-reset">Reset all filters</button>' : '') + '</div>';
  }

  function bindEmpty(root) {
    const b = $('#empty-reset', root);
    if (b) b.addEventListener('click', () => { F = DEFAULT_FILTERS(); renderRail(); compute(); renderView(); });
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
    bindActiveFilters(el); bindEmpty(el);
    bindRows(el, (tr) => openProtein(+tr.dataset.i));
    bindSort(el, UI.sortP, () => renderProteins(el));
    const more = $('#more', el);
    if (more) more.addEventListener('click', () => { UI.pageP += 50; renderProteins(el); });
    $('#dl-proteins', el) && $('#dl-proteins', el).addEventListener('click', () => download('methylarginine_proteins.tsv', proteinTsv(rows)));
    $('#copy-genes', el) && $('#copy-genes', el).addEventListener('click', (e) => copyText(rows.map((p) => p.gene || p.acc).join('\n'), e.currentTarget));
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
      if (sort.key === k) sort.dir = -sort.dir; else { sort.key = k; sort.dir = k === 'gene' || k === 'pos' ? 1 : -1; }
      redraw();
    }));
  }

  // ----------------------------------------------------------------- sites
  const S_SORT = {
    gene: (s) => P[s.p].label.toLowerCase(),
    pos: (s) => s.pos,
    state: (s) => s.st,
    nn: (s) => (s.nn == null ? 1e9 : s.nn),
    src: (s) => popcount(s.study),
  };

  function renderSites(el) {
    const key = S_SORT[UI.sortS.key] ? UI.sortS.key : 'gene';
    const fn = S_SORT[key];
    const dir = UI.sortS.dir;
    const rows = R.sites.slice().sort((a, b) => {
      const x = fn(a), y = fn(b);
      if (x < y) return -dir;
      if (x > y) return dir;
      return a.p - b.p || a.pos - b.pos;
    });
    const shownRows = rows.slice(0, UI.pageS);
    const th = (k, label, cls) => '<th' + (cls ? ' class="' + cls + '"' : '') + ' data-sort="' + k + '"' + (key === k ? ' aria-sort="' + (dir > 0 ? 'ascending' : 'descending') + '"' : '') + '>' + label + '</th>';
    el.innerHTML = resultsHead('Sites', fmt(rows.length) + ' of ' + fmt(S.length) + ' methylarginine sites on ' + fmt(R.prots.length) + ' proteins',
      '<button class="btn small" type="button" id="dl-sites">Download TSV</button>') +
      (rows.length ? '<div class="card tbl-card"><div class="tbl-wrap"><table class="tbl"><thead><tr>' +
        th('gene', 'Protein') + th('pos', 'Site') + th('state', 'State') + '<th>Window ±7</th><th class="hide-sm">Motifs</th>' +
        '<th class="hide-xs">Context</th>' + th('nn', 'Next methyl-R', 'r hide-sm') + '<th class="hide-sm">Maron IP</th><th>PRMTi</th>' + th('src', 'Reported by', 'hide-sm') +
        '</tr></thead><tbody>' +
        shownRows.map((s) => {
          const p = P[s.p];
          return '<tr class="row" tabindex="0" data-j="' + s.j + '"><td><span class="gene">' + esc(p.label) + '</span><span class="acc">' + esc(p.acc) + '</span></td>' +
            '<td class="num"><b>R' + s.pos + '</b></td><td>' + pill(s.st) + '</td><td>' + windowHtml(s.win, 7) + (s.mis ? ' <span class="warn" title="Canonical residue at this position is not R">!</span>' : '') + '</td>' +
            '<td class="nw hide-sm">' + motifTags(s.mot) + '</td><td class="nw hide-xs">' + contextHtml(s) + '</td>' +
            '<td class="r num hide-sm">' + (s.nn == null ? '<span class="muted">—</span>' : s.nn + ' aa') + '</td>' +
            '<td class="nw hide-sm">' + (ipTags(s.ip) || '<span class="muted">—</span>') + '</td><td class="nw">' + treatHtml(s) + '</td>' +
            '<td class="nw hide-sm"><div class="tags">' + studyTags(s.study, 2) + '</div></td></tr>';
        }).join('') +
        '</tbody></table></div><div class="tbl-foot"><span>Showing ' + fmt(shownRows.length) + ' of ' + fmt(rows.length) + '</span>' +
        (rows.length > shownRows.length ? '<button class="btn quiet small" type="button" id="more">Show 100 more</button>' : '') + '</div></div>' : emptyState());
    bindActiveFilters(el); bindEmpty(el);
    bindRows(el, (tr) => { const s = S[+tr.dataset.j]; openProtein(s.p, s.pos); });
    bindSort(el, UI.sortS, () => renderSites(el));
    const more = $('#more', el);
    if (more) more.addEventListener('click', () => { UI.pageS += 100; renderSites(el); });
    $('#dl-sites', el) && $('#dl-sites', el).addEventListener('click', () => download('methylarginine_sites.tsv', siteTsv(rows)));
  }

  // -------------------------------------------------------------- overview
  function hbars(rows, opts) {
    // rows: {label, value, sig, title, id}
    const vals = rows.map((r) => r.value).filter((v) => v != null);
    let lo = Math.min(0, ...vals), hi = Math.max(0, ...vals);
    const raw = (hi - lo || 1) / 5;
    const mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const step = [1, 2, 5, 10].map((m) => m * mag).find((st) => st >= raw);
    lo = Math.floor(lo / step) * step; hi = Math.ceil(hi / step) * step;
    const span = hi - lo;
    const xp = (v) => ((v - lo) / span) * 100;
    const ticks = [];
    for (let t = lo; t <= hi + step / 2; t += step) ticks.push(+t.toFixed(6));
    const vfmt = opts.format || fmtL2;
    const tfmt = opts.tick || ((t) => t);
    return '<div class="hbars" role="list">' + rows.map((r, k) => {
      const v = r.value == null ? 0 : r.value;
      const a = xp(Math.min(0, v)), b = xp(Math.max(0, v));
      return '<button class="hb" type="button" role="listitem" data-k="' + k + '" aria-label="' + esc(r.label + ': ' + (opts.valueLabel || 'log2 OR') + ' ' + vfmt(r.value) + (opts.noSig ? '' : r.sig ? ', significant' : ', not significant')) + '">' +
        '<span class="hb-label">' + esc(r.label) + '</span><span class="hb-track"><span class="hb-zero" style="left:' + xp(0).toFixed(2) + '%"></span>' +
        '<span class="hb-bar' + (r.sig ? '' : ' ns') + (v < 0 ? ' neg' : '') + '" style="left:' + a.toFixed(2) + '%;width:' + Math.max(0.6, b - a).toFixed(2) + '%"></span></span>' +
        '<span class="hb-val">' + vfmt(r.value) + '</span></button>';
    }).join('') + '</div>' +
      '<div class="hb-axis" aria-hidden="true"><span></span><span class="ticks">' + ticks.map((t) => '<span style="left:' + xp(t).toFixed(2) + '%">' + tfmt(t) + '</span>').join('') + '</span><span></span></div>' +
      '<div class="chart-key">' + (opts.noSig ? '' : '<span><i></i>q ≤ 0.05</span><span><i class="ns"></i>not significant</span>') + '<span>' + esc(opts.axis || 'log2 odds ratio vs background') + '</span></div>';
  }

  function bindBars(container, rows, onClick) {
    $$('.hb', container).forEach((b) => {
      const r = rows[+b.dataset.k];
      b.addEventListener('click', () => onClick(r));
      b.addEventListener('mouseenter', (e) => showTip(e, r.title));
      b.addEventListener('mousemove', moveTip);
      b.addEventListener('mouseleave', hideTip);
    });
  }

  function renderOverview(el) {
    const ov = M.overview;
    const topSites = P.slice().sort((a, b) => b.sites.length - a.sites.length).slice(0, 14);
    const topArg = P.filter((p) => p.argQ != null && p.argQ <= PADJ && p.shrink != null).sort((a, b) => b.shrink - a.shrink).slice(0, 14);
    const topDense = P.filter((p) => p.sig.length).sort((a, b) => Math.max(...b.sig.map((w) => w[4])) - Math.max(...a.sig.map((w) => w[4]))).slice(0, 14);
    const examples = [['FUS', 'One protein'], ['FUS R216', 'One site'], ['TAF15 EWSR1 FUS', 'A gene list (FET family)'], ['SNRPB SNRPD1 SNRPD3', 'Sm proteins'], ['heterogeneous nuclear ribonucleoprotein', 'Protein name']];
    const geneChips = (list, sub) => '<div class="gene-chips">' + list.map((p) => '<button class="gene-chip" type="button" data-i="' + p.i + '">' + esc(p.label) + '<span>' + sub(p) + '</span></button>').join('') + '</div>';

    const motifRows = ov.motifs.map((r) => ({ label: MOTIF_LABEL[r[0]] + (MOTIF_HELP[r[0]] ? ' · ' + MOTIF_HELP[r[0]] : ''), value: r[4], sig: r[5] != null && r[5] <= PADJ, key: r[0],
      title: '<b>' + esc(MOTIF_LABEL[r[0]]) + '</b> (' + esc(MOTIF_HELP[r[0]]) + ')<br>' + fmt(r[1]) + ' methyl-R sites (' + (100 * r[2]).toFixed(1) + '%) vs ' + (100 * r[3]).toFixed(1) + '% of unmodified R in the same proteins<br>log2 OR ' + fmtL2(r[4]) + ', q ' + fmtP(r[5]) }));
    const funcRows = ov.functional_classes.map((r) => ({ label: r[0], value: r[4], sig: r[5] != null && r[5] <= PADJ, key: r[0],
      title: '<b>' + esc(r[0]) + '</b><br>' + fmt(r[1]) + ' methylarginine proteins (' + (100 * r[2]).toFixed(1) + '%) vs ' + (100 * r[3]).toFixed(1) + '% of reviewed background<br>log2 OR ' + fmtL2(r[4]) + ', q ' + fmtP(r[5]) }))
      .sort((a, b) => b.value - a.value);
    const condRows = Object.entries(M.condensate_stats).map(([k, v]) => ({ label: M.condensates[+k], value: v[2], sig: v[3] != null && v[3] <= PADJ, key: +k,
      title: '<b>' + esc(M.condensates[+k]) + '</b><br>' + v[1] + ' of ' + v[0] + ' members methylated<br>log2 OR ' + fmtL2(v[2]) + ', q ' + fmtP(v[3]) }))
      .sort((a, b) => b.value - a.value);
    const disRows = ov.disease_families.map((r) => ({ label: r[0], value: r[3] ? Math.log2(r[3]) : null, sig: r[4] != null && r[4] <= PADJ, key: r[0],
      title: '<b>' + esc(r[0]) + '</b><br>' + fmt(r[1]) + ' methylarginine genes of ' + fmt(r[2]) + ' family genes<br>OR ' + (r[3] || 0).toFixed(2) + ', BH-adjusted p ' + fmtP(r[4]) }))
      .sort((a, b) => b.value - a.value);
    const ctxRows = (table, labels, help) => table.map((r) => ({ label: labels[r[0]] || r[0], value: r[4], sig: r[5] != null && r[5] <= PADJ, key: r[0],
      title: '<b>' + esc(labels[r[0]] || r[0]) + '</b> · ' + esc(help[r[0]] || '') + '<br>' + fmt(r[1]) + ' sites (' + (100 * r[2]).toFixed(1) + '%) vs ' + (100 * r[3]).toFixed(1) + '% of unmodified R in the same proteins<br>log2 OR ' + fmtL2(r[4]) + ', q ' + fmtP(r[5]) }));
    const studyRows = M.studyOrder.map((k) => ({ label: M.studies[k][0], value: M.studies[k][2], sig: true, key: k,
      title: '<b>' + esc(M.studies[k][0]) + '</b><br>' + esc(M.studies[k][1] || STUDY_HELP[M.studies[k][0]] || '') + '<br>' + fmt(M.studies[k][2]) + ' sites, ' + fmt(M.studies[k][3]) + ' reported only here' }));
    const repDist = [0, 0, 0, 0, 0];
    for (const s of S) repDist[Math.min(5, popcount(s.study)) - 1]++;
    const disorderRows = ctxRows(ov.disorder.filter((r) => r[0] !== 'no_disorder_annotation'), DISORDER_LABEL, DISORDER_HELP);
    const domainRows = ctxRows(ov.domain_context, DOMAIN_LABEL, DOMAIN_HELP);

    el.innerHTML =
      '<div class="ov-intro"><div class="card"><h3>What you can ask</h3>' +
      '<p>This explorer indexes the integrated human methylarginine union from the <a href="https://github.com/Shechterlab/PTM_Pipeline">PTM_Pipeline</a>: ' +
      fmt(M.n_sites) + ' canonically remapped sites on ' + fmt(M.n_proteins) + ' proteins, compiled from Maron et al. 2021, ProMetheusDB and curated UniProt/IEDB annotations, with each site annotated for motif, disorder and domain context, local clustering, nearby PTMs, PRMT inhibitor response, condensate membership and disease association.</p>' +
      '<p>Search for a protein or site, paste a gene list, or click any bar below to filter. Open a protein for its site map, sequence and evidence.</p>' +
      '<div class="quick"><div><h4>Try</h4><div class="gene-chips">' + examples.map((e) => '<button class="gene-chip" type="button" data-q="' + esc(e[0]) + '">' + esc(e[0]) + '<span>' + esc(e[1]) + '</span></button>').join('') + '</div></div></div></div>' +
      '<div class="card quick"><div><h4>Most methyl-R sites</h4>' + geneChips(topSites, (p) => p.sites.length) + '</div>' +
      '<div><h4>Strongest arginine-methylation enrichment</h4>' + geneChips(topArg, (p) => fmtL2(p.shrink)) + '</div>' +
      '<div><h4>Densest methyl-R clusters</h4>' + geneChips(topDense, (p) => { const w = p.sig.slice().sort((a, b) => b[4] - a[4])[0]; return w[3] + '/' + w[0] + ' aa'; }) + '</div></div></div>' +

      '<div class="ov-grid">' +
      '<div class="card wide" id="ch-study"><h3>Provenance</h3><p class="sub">Studies compiled in Maron et al. 2021 Table S5, plus ProMetheusDB and curated UniProt/IEDB records. Bars show sites reported by each; hover for the number reported only there. Click to filter.</p>' +
      '<div class="two-col"><div>' + hbars(studyRows, { noSig: true, format: (v) => fmt(v), tick: compact, axis: 'methylarginine sites reported', valueLabel: 'sites' }) + '</div>' +
      '<div><p class="flabel" style="margin:0 0 .5rem">Independent reports per site</p><div class="hbars rep-dist">' +
      repDist.map((n, k) => '<button class="hb" type="button" data-rep="' + (k + 1) + '"><span class="hb-label">' + (k === 4 ? '5 or more' : (k + 1) + (k ? ' reports' : ' report')) + '</span><span class="hb-track"><span class="hb-bar" style="left:0;width:' + (100 * n / repDist[0]).toFixed(1) + '%"></span></span><span class="hb-val">' + fmt(n) + '</span></button>').join('') +
      '</div><p class="note">A report is one study or database. ' + fmt(M.n_sites - repDist[0]) + ' sites (' + pct(M.n_sites - repDist[0], M.n_sites) + '%) have two or more. Click a row to keep sites with at least that many.</p></div></div></div>' +
      '<div class="card" id="ch-motif"><h3>Sequence motifs</h3><p class="sub">Motif frequency at methylated vs unmodified arginines in the same proteins. Click a motif to list its sites.</p>' + hbars(motifRows, {}) + '</div>' +
      '<div class="card" id="ch-cond"><h3>Biomolecular condensates</h3><p class="sub">CD-CODE condensates with ≥5 members: methylated fraction vs the global rate. ' + fmt(M.condensate_summary.methyl_proteins_in_any_condensate) + ' methylarginine proteins sit in at least one condensate. Click to list members.</p>' +
      hbars(UI.allCond ? condRows : condRows.slice(0, 14), {}) + '<button class="linkish" type="button" id="cond-more" style="margin-top:.5rem">' + (UI.allCond ? 'Show the top 14' : 'Show all ' + condRows.length + ' condensates') + '</button></div>' +
      '<div class="card" id="ch-func"><h3>Functional classes</h3><p class="sub">Multi-label functional classes among methylarginine proteins vs the reviewed human background.</p>' + hbars(funcRows, {}) + '</div>' +
      '<div class="card" id="ch-dis"><h3>Disease families</h3><p class="sub">Audited DisGeNET families (genes counted once per family; BH across 8 families). Enrichment of disease gene sets, not causation.</p>' + hbars(disRows, { axis: 'log2 odds ratio vs reviewed background' }) + '</div>' +
      '<div class="card" id="ch-dis2"><h3>Disorder context</h3><p class="sub">Methylated vs unmodified arginines in the same proteins (MobiDB intervals). Sites on proteins without a MobiDB annotation are left out because their background has only a few dozen arginines.</p>' + hbars(disorderRows, {}) + '</div>' +
      '<div class="card" id="ch-dom"><h3>Domain context</h3><p class="sub">Methylated vs unmodified arginines in the same proteins (InterPro domains and repeats).</p>' + hbars(domainRows, {}) + '</div>' +
      '<div class="card wide" id="ch-treat"><h3>PRMT inhibitor response</h3><p class="sub">Maron et al. 2021 Table S3: sites per methyl-IP that changed significantly (adj. p ≤ 0.05) after GSK591 (PRMT5 inhibitor) or MS023 (type I PRMT inhibitor) vs DMSO. Click a row to list those sites.</p>' + treatTable(ov.treatments) + '</div>' +
      '</div>';

    $$('.gene-chip[data-i]', el).forEach((b) => b.addEventListener('click', () => openProtein(+b.dataset.i)));
    $$('.gene-chip[data-q]', el).forEach((b) => b.addEventListener('click', () => {
      setQuery(b.dataset.q); compute();
      if (R.prots.length === 1) openProtein(R.prots[0].i, Q.positions ? Array.from(Q.positions.values())[0].values().next().value : null);
      else go('proteins');
    }));
    bindBars($('#ch-study', el), studyRows, (r) => setFilters((f) => { f.studies = new Set([r.key]); }, 'sites'));
    $$('.rep-dist .hb', el).forEach((b) => b.addEventListener('click', () => setFilters((f) => { f.minRep = +b.dataset.rep > 1 ? +b.dataset.rep : 0; }, 'sites')));
    bindBars($('#ch-motif', el), motifRows, (r) => setFilters((f) => { f.motifs = new Set([M.motifs.indexOf(r.key)]); }, 'sites'));
    bindBars($('#ch-cond', el), UI.allCond ? condRows : condRows.slice(0, 14), (r) => setFilters((f) => { f.cond = String(r.key); }, 'proteins'));
    $('#cond-more', el).addEventListener('click', () => { UI.allCond = !UI.allCond; renderOverview(el); });
    bindBars($('#ch-func', el), funcRows, (r) => setFilters((f) => { f.func = String(M.func_classes.indexOf(r.key)); }, 'proteins'));
    bindBars($('#ch-dis', el), disRows, (r) => setFilters((f) => { f.disease = String(M.disease_families.indexOf(r.key)); }, 'proteins'));
    bindBars($('#ch-dis2', el), disorderRows, (r) => setFilters((f) => { const k = M.disorder_classes.indexOf(r.key); f.disorder = new Set(k >= 0 ? [k] : []); }, 'sites'));
    bindBars($('#ch-dom', el), domainRows, (r) => setFilters((f) => { const k = M.domain_contexts.indexOf(r.key); f.domctx = new Set(k >= 0 ? [k] : []); }, 'sites'));
    $$('#ch-treat tr.row', el).forEach((tr) => {
      const go2 = () => setFilters((f) => { f.treat = tr.dataset.t + ':' + tr.dataset.dir; f.ips = new Set([+tr.dataset.ip]); }, 'sites');
      tr.addEventListener('click', go2);
      tr.addEventListener('keydown', (e) => { if (e.key === 'Enter') go2(); });
    });
  }

  function treatTable(rows) {
    const max = Math.max(...rows.map((r) => r[2] ? Math.max(r[3], r[4]) / r[2] : 0), 0.01);
    const order = rows.slice().sort((a, b) => a[0].localeCompare(b[0]) || M.maron_ips.indexOf(a[1]) - M.maron_ips.indexOf(b[1]));
    return '<div class="tbl-wrap"><table class="tbl treat-tbl"><thead><tr><th>Inhibitor</th><th>Methyl-IP</th><th class="r">Sites</th><th>Decreased</th><th class="hide-xs">Increased</th><th class="r hide-sm">Median change</th></tr></thead><tbody>' +
      order.map((r) => {
        const t = M.maron_treatments.indexOf(r[0]);
        const ip = M.maron_ips.indexOf(r[1]);
        const down = r[2] ? r[3] / r[2] : 0, up = r[2] ? r[4] / r[2] : 0;
        return '<tr class="row" tabindex="0" data-t="' + t + '" data-ip="' + ip + '" data-dir="down"><td><b>' + esc(r[0]) + '</b> <span class="muted">' + esc(TREAT_HELP[r[0]] || '') + '</span></td><td><span class="tag ip">' + esc(r[1]) + '</span></td><td class="r num">' + fmt(r[2]) + '</td>' +
          '<td><div style="display:flex;gap:.5rem;align-items:center"><span class="meter" style="flex:1"><i style="width:' + (100 * down / max).toFixed(1) + '%"></i></span><span class="num">' + r[3] + ' (' + (100 * down).toFixed(0) + '%)</span></div></td>' +
          '<td class="hide-xs"><div style="display:flex;gap:.5rem;align-items:center"><span class="meter up" style="flex:1"><i style="width:' + (100 * up / max).toFixed(1) + '%"></i></span><span class="num">' + r[4] + ' (' + (100 * up).toFixed(0) + '%)</span></div></td>' +
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
    const funcs = bits(p.func).map((b) => M.func_classes[b]);
    const subs = bits(p.sub).map((b) => M.func_subclasses[b]).filter((s) => !funcs.includes(s));
    const ext = (href, label) => '<a class="btn quiet small" href="' + href + '" target="_blank" rel="noopener">' + label + ' ↗</a>';
    el.innerHTML =
      '<div class="detail-head"><button class="linkish back" type="button" id="back">← Back to results</button>' +
      '<div class="detail-title"><h2>' + esc(p.label) + '</h2><span class="pn" id="pname">' + esc(p.name) + '</span></div>' +
      '<div class="detail-links"><span class="tag mono">' + esc(p.acc) + '</span>' + (p.uid ? '<span class="tag mono">' + esc(p.uid) + '</span>' : '') +
      (p.iso ? '<span class="tag" title="Sites were reported on this isoform accession">isoform</span>' : '') +
      ext('https://www.uniprot.org/uniprotkb/' + encodeURIComponent(p.acc) + '/entry', 'UniProt') +
      ext('https://alphafold.ebi.ac.uk/entry/' + encodeURIComponent(p.acc.split('-')[0]), 'AlphaFold') +
      ext('https://www.ebi.ac.uk/interpro/protein/UniProt/' + encodeURIComponent(p.acc.split('-')[0]) + '/', 'InterPro') + '</div>' +
      (funcs.length ? '<div class="tags">' + funcs.map((f) => '<span class="tag">' + esc(f) + '</span>').join('') + subs.map((f) => '<span class="tag src">' + esc(f) + '</span>').join('') + '</div>' : '') +
      '</div>' +
      '<div class="facts">' +
      fact(fmt(p.sites.length), 'Methylarginine sites', shownN === p.sites.length ? 'all shown' : fmt(shownN) + ' match current filters') +
      fact('<span id="arg-count">' + (argN == null ? '—' : fmt(argN)) + '</span>', 'Arginines', '<span id="arg-frac">' + (argN ? pct(p.sites.length, argN) + '% methylated' : 'loading from UniProt') + '</span>') +
      fact(p.shrink == null ? '—' : fmtL2(p.shrink), 'Arg enrichment (log2)', p.argQ == null ? 'not tested (needs sequence-confirmed sites)' : 'shrinkage estimate; Fisher q ' + fmtP(p.argQ)) +
      fact(best ? best[3] + ' / ' + best[0] + ' aa' : '—', 'Densest window', best ? 'R' + best[1] + '–' + best[2] + ', q ' + fmtP(best[5]) + (best[5] <= PADJ ? '' : ' (n.s.)') : 'needs ≥3 sequence-confirmed sites') +
      fact(nnMed == null ? '—' : nnMed + ' aa', 'Median methyl-R spacing', p.dens ? 'closest pair ' + p.dens[0] + ' aa apart' : '') +
      fact(fmt(p.cond.length), 'CD-CODE condensates', p.cond.length ? '' : 'none annotated') +
      '</div>' +
      '<div class="card"><div class="map-tools"><div><h3>Protein map</h3><p class="sub" id="map-sub" style="margin:0">Loading UniProt features…</p></div>' +
      '<div class="actions"><button class="btn quiet small" type="button" id="zoom-reset" hidden>Reset zoom</button></div></div>' +
      '<div class="legend">' + [1, 2, 3, 0].map((k) => '<span><i class="st-' + k + '"></i>' + STATE_LABEL[k] + '</span>').join('') +
      '<span><i class="sq" style="background:var(--feat-domain)"></i>Domain / repeat</span><span><i class="sq" style="background:var(--dense)"></i>Dense methyl-R window</span><span class="muted">Drag across the map to zoom</span></div>' +
      '<div class="map-box" id="map-box"></div></div>' +
      '<div class="card tbl-card"><div style="padding:1.25rem 1.25rem 0.5rem"><h3>Sites</h3><p class="sub" style="margin:0">Every reported methylarginine on ' + esc(p.label) + '. Faded rows do not match the current filters. Expand a row for source evidence.</p></div>' +
      '<div class="tbl-wrap" id="psites"></div></div>' +
      '<div class="card"><h3>Sequence</h3><p class="sub" id="seq-sub">Loading sequence from UniProt…</p><div id="seq"></div></div>' +
      '<div class="two-col"><div class="card" id="card-cond"></div><div class="card" id="card-dis"></div></div>' +
      '<div class="two-col"><div class="card" id="card-cluster"></div><div class="card" id="card-ptm"></div></div>';

    $('#back', el).addEventListener('click', () => { if (history.length > 1) history.back(); else go('proteins'); });
    renderSiteTable(p);
    renderCondCard(p);
    renderDiseaseCard(p, null);
    renderClusterCard(p);
    renderPtmCard(p);
    drawMap(p, null);
    loadEvidence().then((ev) => { if (UI.detail === p.i && UI.view === 'protein') renderDiseaseCard(p, ev); });
    loadUniProt(p.acc).then((u) => {
      if (UI.detail !== p.i || UI.view !== 'protein') return;
      drawMap(p, u);
      renderSequence(p, u);
      if (u) {
        if (!p.name && u.name) $('#pname').textContent = u.name;
        if (argN == null && u.seq) {
          const n = (u.seq.match(/R/g) || []).length;
          $('#arg-count').textContent = fmt(n);
          $('#arg-frac').textContent = pct(p.sites.length, n) + '% methylated (UniProt sequence)';
        }
      } else if (argN == null) $('#arg-frac').textContent = 'UniProt unavailable';
    });
    if (UI.highlight != null) setTimeout(() => flashSite(UI.highlight, true), 60);
  }

  function fact(value, label, sub) {
    return '<div class="fact"><b>' + value + '</b><span>' + label + '</span>' + (sub ? '<small>' + sub + '</small>' : '') + '</div>';
  }

  function renderSiteTable(p) {
    const box = $('#psites');
    const rows = p.sites.map((j) => S[j]);
    box.innerHTML = '<table class="tbl"><thead><tr><th>Site</th><th>State</th><th>Window ±7</th><th class="hide-sm">Motifs</th><th class="hide-xs">Context</th>' +
      '<th class="hide-sm">Nearest domain</th><th class="r hide-sm">Next methyl-R</th><th class="hide-sm">Maron</th><th>PRMTi</th><th class="hide-sm">Reported by</th><th><span class="sr-only">Evidence</span></th></tr></thead><tbody>' +
      rows.map((s) => {
        const dim = !R.shown.has(s.j);
        const ipr = s.ipr >= 0 ? M.interpro[s.ipr] : null;
        const dom = ipr ? '<span title="' + esc(ipr[0]) + '">' + esc(ipr[1]) + '</span><br><span class="muted">' + (s.domDist === 0 ? 'inside' : s.domDist + ' aa ' + (s.domSide === 'N_terminal_to_domain' ? 'N-term' : 'C-term')) + '</span>' : '<span class="muted">—</span>';
        const s2 = bits(s.s2).map((b) => ['me1', 'me2'][b]).join('/');
        const open = UI.openEvidence.has(s.j);
        return '<tr class="row" tabindex="0" data-j="' + s.j + '" data-pos="' + s.pos + '"' + (dim ? ' style="opacity:.5"' : '') + '><td class="num"><b>R' + s.pos + '</b>' + (s.orig ? '<br><span class="muted" title="Position as reported by the source before canonical remapping">reported ' + esc(Array.from(new Set(s.orig.split(';'))).join(', ')) + '</span>' : '') + '</td>' +
          '<td>' + pill(s.st) + '</td><td>' + windowHtml(s.win, 7) + (s.mis ? ' <span class="warn" title="Canonical residue at this position is not R">!</span>' : '') + '</td>' +
          '<td class="nw hide-sm">' + motifTags(s.mot) + '</td><td class="nw hide-xs">' + contextHtml(s) + '</td><td class="hide-sm">' + dom + '</td>' +
          '<td class="r num hide-sm">' + (s.nn == null ? '—' : s.nn + ' aa') + '</td>' +
          '<td class="nw hide-sm">' + ([ipTags(s.ip), s2 ? '<span class="tag">S2 ' + s2 + '</span>' : ''].filter(Boolean).join(' ') || '<span class="muted">—</span>') + '</td>' +
          '<td class="nw">' + treatHtml(s) + '</td><td class="nw hide-sm"><div class="tags">' + studyTags(s.study, 3) + '</div></td>' +
          '<td><button class="btn quiet small ev-toggle" type="button" aria-expanded="' + open + '" data-j="' + s.j + '">' + (open ? 'Hide' : 'Evidence') + '</button></td></tr>' +
          (open ? '<tr class="ev-row" data-for="' + s.j + '"><td colspan="11" class="ev" id="ev-' + s.j + '">Loading evidence…</td></tr>' : '');
      }).join('') + '</tbody></table>';
    bindRows(box, (tr) => { flashSite(+tr.dataset.pos, false); });
    $$('.ev-toggle', box).forEach((b) => b.addEventListener('click', () => {
      const j = +b.dataset.j;
      if (UI.openEvidence.has(j)) UI.openEvidence.delete(j); else UI.openEvidence.add(j);
      renderSiteTable(p);
    }));
    if (UI.openEvidence.size) loadEvidence().then((ev) => UI.openEvidence.forEach((j) => fillEvidence(j, ev)));
  }

  function fillEvidence(j, ev) {
    const cell = document.getElementById('ev-' + j);
    if (!cell) return;
    const s = S[j];
    const p = P[s.p];
    const rows = ev ? ev.sites[j] : null;
    let html = '';
    if (s.tr.length) {
      html += '<p><b>PRMT inhibitor response</b> (Maron et al. 2021 Table S3, change vs DMSO, adj. p)</p><table><thead><tr><th>Methyl-IP</th><th>Inhibitor</th><th>Change</th><th>adj. p</th><th>Peptides</th></tr></thead><tbody>' +
        s.tr.slice().sort((a, b) => a[0] - b[0] || a[1] - b[1]).map(([ip, t, d, q, n]) => '<tr><td>' + M.maron_ips[ip] + '</td><td>' + M.maron_treatments[t] + ' <span class="muted">' + TREAT_HELP[M.maron_treatments[t]] + '</span></td><td class="num">' + fmtL2(d) + '</td><td class="num">' + (q != null && q <= PADJ ? '<b>' + fmtP(q) + '</b>' : fmtP(q)) + '</td><td class="num">' + (n == null ? '—' : n) + '</td></tr>').join('') + '</tbody></table>';
    }
    const near = [];
    for (const [pos, g, res] of p.ptms) { const d = pos - s.pos; if (Math.abs(d) <= 10) near.push(PTM_PREFIX[g] + res + pos + ' (' + (d > 0 ? '+' : '') + d + ')'); }
    if (near.length) html += '<p><b>Other PTMs within ±10 aa:</b> ' + esc(near.join(', ')) + '</p>';
    if (s.win) html += '<p><b>Window ±15:</b> ' + windowHtml(s.win, 15) + '</p>';
    if (rows && rows.length) {
      html += '<p><b>Source rows</b> (' + rows.length + ')</p><table><thead><tr><th>Study / database</th><th>Source table</th><th>Reference</th><th>PMID</th><th>Peptide</th><th>Enzyme</th><th>State</th><th>Reported as</th></tr></thead><tbody>' +
        rows.map((r) => '<tr><td><b>' + esc(r[9] != null ? ev.studies[r[9]] : '—') + '</b></td><td>' + esc(ev.source_tables[r[0]]) + '</td><td>' + esc(r[1] >= 0 ? ev.refs[r[1]] : '—') + '</td>' +
          '<td>' + (r[2] ? r[2].split(/[;,\s]+/).filter(Boolean).map((id) => '<a href="https://pubmed.ncbi.nlm.nih.gov/' + encodeURIComponent(id) + '/" target="_blank" rel="noopener">' + esc(id) + '</a>').join(' ') : '—') + '</td>' +
          '<td class="pep">' + esc(r[3] || '—') + '</td><td>' + esc(r[4] || '—') + '</td><td>' + esc(r[5] || '—') + '</td>' +
          '<td class="mono">' + esc(r[6]) + ':' + esc(r[7]) + (r[7] && String(r[7]) !== String(s.pos) ? ' <span class="muted">→ R' + s.pos + '</span>' : '') + '<br><span class="muted">' + esc(ev.remap_statuses[r[8]].replace(/_/g, ' ')) + '</span></td></tr>').join('') +
        '</tbody></table>';
    } else if (!ev) html += '<p class="note">The evidence file did not load.</p>';
    if (s.refs.length && !(rows && rows.length)) html += '<p><b>References:</b> ' + s.refs.map((k) => esc(M.refs[k])).join('; ') + '</p>';
    cell.innerHTML = html || '<span class="muted">No additional evidence.</span>';
  }

  function flashSite(pos, scroll) {
    const tr = $('#psites tr.row[data-pos="' + pos + '"]');
    $$('.pmap .head').forEach((h) => h.classList.toggle('sel', +h.dataset.pos === pos));
    if (!tr) return;
    if (scroll) tr.scrollIntoView({ block: 'center', behavior: 'smooth' });
    tr.classList.remove('flash'); void tr.offsetWidth; tr.classList.add('flash');
  }

  function renderCondCard(p) {
    const box = $('#card-cond');
    const cond = p.cond.map((k) => [k, M.condensates[k], M.condensate_stats[String(k)]]).sort((a, b) => a[1].localeCompare(b[1]));
    const funcs = bits(p.func).map((b) => M.func_classes[b]);
    box.innerHTML = '<h3>Condensates &amp; function</h3><p class="sub">CD-CODE membership; click one to list its methylated members.</p>' +
      (cond.length ? '<div class="chips">' + cond.map(([k, name, st]) => '<button class="chip" type="button" data-k="' + k + '"' + (st ? ' title="Condensate-level enrichment: log2 OR ' + fmtL2(st[2]) + ', q ' + fmtP(st[3]) + '"' : '') + '>' + esc(name) + '</button>').join('') + '</div>' : '<p class="muted">No CD-CODE condensate membership recorded.</p>') +
      '<dl class="kv" style="margin-top:1rem"><dt>Functional class</dt><dd>' + (funcs.map(esc).join('; ') || '—') + '</dd>' +
      '<dt>Subclass</dt><dd>' + (bits(p.sub).map((b) => esc(M.func_subclasses[b])).join('; ') || '—') + '</dd></dl>';
    $$('.chip', box).forEach((b) => b.addEventListener('click', () => setFilters((f) => { f.cond = b.dataset.k; }, 'proteins')));
  }

  function renderDiseaseCard(p, ev) {
    const box = $('#card-dis');
    if (!box) return;
    const fams = bits(p.dis).map((b) => M.disease_families[b]);
    let terms = '';
    if (ev) {
      const list = (ev.protein_disease_terms[p.i] || []).map((k) => ev.disease_terms[k]).sort((a, b) => (a[4] || 1) - (b[4] || 1));
      terms = list.length ? '<details' + (list.length <= 8 ? ' open' : '') + '><summary class="linkish" style="display:inline">' + list.length + ' enriched DisGeNET term' + (list.length === 1 ? '' : 's') + ' include ' + esc(p.label) + '</summary><ul class="term-list">' +
        list.map((t) => '<li>' + esc(t[1]) + ' <span class="muted">' + (t[2] >= 0 ? '· ' + esc(M.disease_families[t[2]]) + ' ' : '') + '· p.adj ' + fmtP(t[4]) + '</span></li>').join('') + '</ul></details>' : '';
    }
    box.innerHTML = '<h3>Disease associations</h3><p class="sub">DisGeNET terms enriched among methylarginine proteins (clusterProfiler, p.adjust ≤ 0.05) that contain this gene, grouped by audited family. Association, not causation.</p>' +
      (fams.length || p.als ? '<div class="tags">' + fams.map((f) => '<span class="tag">' + esc(f) + '</span>').join('') + (p.als ? '<span class="tag ip" title="Member of the full audited ALS / motor neuron / FTD family gene set">ALS/MN/FTD audited set</span>' : '') + '</div>' : '<p class="muted">Not in an enriched disease term.</p>') +
      (terms ? '<div style="margin-top:.75rem">' + terms + '</div>' : ev ? '' : '<p class="note">Loading terms…</p>');
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
        sig.map((w, k) => '<button class="chip" type="button" data-k="' + k + '" title="log2 obs/null ' + fmtL2(w[4]) + ', q ' + fmtP(w[5]) + '">' + w[3] + ' sites / ' + w[0] + ' aa · ' + w[1] + '–' + w[2] + '</button>').join('') + '</div><p class="note">Click a window to zoom the map.</p>' : '');
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
    box.innerHTML = '<h3>Other PTMs</h3><p class="sub">Sites of other modification types on this protein from the cross-PTM master table (HPRD, UniProt and other curated sources).</p>' +
      (p.ptms.length ? '<dl class="kv">' + counts.map((c, g) => c ? '<dt>' + PTM_LABEL[g] + '</dt><dd>' + c + '</dd>' : '').join('') + '</dl>' +
        (near.length ? '<p style="margin:.9rem 0 .4rem"><b>Within 10 aa of a methyl-R</b></p><div class="tags">' + near.slice(0, 40).map((n) => '<span class="tag" title="' + n[4] + ' aa from R' + n[3] + '">' + PTM_PREFIX[n[1]] + esc(n[2]) + n[0] + ' · R' + n[3] + '</span>').join('') + '</div>' + (near.length > 40 ? '<p class="note">' + (near.length - 40) + ' more not shown.</p>' : '') : '<p class="note">None within 10 aa of a methylarginine.</p>')
        : '<p class="muted">No other PTMs recorded for this accession.</p>');
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
    sub.innerHTML = 'Current UniProt sequence, ' + fmt(seq.length) + ' aa. Methylated arginines are filled by state, other arginines are bold and glycines are muted. Click a filled R to find its row.' +
      (mismatch ? ' <span class="warn">' + mismatch + ' reported site' + (mismatch === 1 ? ' does' : 's do') + ' not fall on R in this sequence version; positions follow the pipeline\'s staged canonical FASTA.</span>' : '');
    let lines = '';
    let nums = '';
    for (let start = 0; start < seq.length; start += 60) {
      nums += '<span class="ln">' + (start + 1) + '</span>';
      let line = '';
      for (let k = start; k < Math.min(seq.length, start + 60); k++) {
        if (k > start && (k - start) % 10 === 0) line += ' ';
        const ch = seq[k];
        const s = siteAt.get(k + 1);
        if (s && ch === 'R') line += '<span class="m st' + s.st + '" data-pos="' + (k + 1) + '" title="R' + (k + 1) + ' · ' + STATE_LABEL[s.st] + '">R</span>';
        else if (ch === 'R') line += '<span class="R">R</span>';
        else if (ch === 'G') line += '<span class="G">G</span>';
        else line += ch;
      }
      lines += '<span class="ln">' + (start + 1) + '</span><span class="line">' + line + '</span>';
    }
    box.innerHTML = '<div class="seq">' + lines + '</div>';
    $$('.seq span.m', box).forEach((m) => m.addEventListener('click', () => flashSite(+m.dataset.pos, true)));
  }

  // ------------------------------------------------------------- protein map
  function drawMap(p, u) {
    const box = $('#map-box');
    if (!box) return;
    const sub = $('#map-sub');
    const len = (u && u.len) || p.len || Math.ceil(p.maxPos * 1.05) || 100;
    const feats = u ? u.feats : [];
    if (sub) {
      sub.textContent = u ? fmt(len) + ' aa · features from UniProt' + (feats.length ? '' : ' (none annotated)') :
        (uniprotCache.has(p.acc) ? 'UniProt features unavailable; showing sites on a ' + (p.len ? fmt(len) + ' aa' : 'length-estimated') + ' backbone.' : 'Loading UniProt features…');
    }
    const z = UI.zoom || [1, len];
    const v0 = Math.max(1, z[0]), v1 = Math.min(len, Math.max(z[1], v0 + 10));
    $('#zoom-reset').hidden = !UI.zoom;
    const W = 1000, ML = 86, MR = 14;
    const x = (pos) => ML + ((pos - v0) / (v1 - v0 + 1)) * (W - ML - MR);
    const visible = (a, b) => b >= v0 && a <= v1;
    const clampX = (pos) => Math.max(ML, Math.min(W - MR, x(pos)));
    const sites = p.sites.map((j) => S[j]).filter((s) => s.pos >= v0 && s.pos <= v1);

    // Lay out lollipop heads on levels so dense clusters stay legible.
    const R0 = 5;
    const levels = [];
    const heads = sites.map((s) => {
      const cx = x(s.pos);
      let lv = 0;
      while (levels[lv] != null && cx - levels[lv] < 2 * R0 + 1 && lv < 9) lv++;
      levels[lv] = cx;
      return { s, cx, lv };
    });
    const nLv = Math.max(1, levels.length);
    const top = 14;
    const base = top + 26 + nLv * 12;
    let y = base;
    let svg = '';
    const label = (txt, yy) => { svg += '<text class="track-label" x="0" y="' + (yy + 4) + '">' + txt + '</text>'; };

    // Dense windows (background bands behind lollipops)
    for (const w of p.sig) if (visible(w[1], w[2])) svg += '<rect class="dense" x="' + clampX(w[1]).toFixed(1) + '" y="' + top + '" width="' + Math.max(2, clampX(w[2]) - clampX(w[1])).toFixed(1) + '" height="' + (base - top + 12) + '" rx="3"><title>Dense window: ' + w[3] + ' sites in ' + w[0] + ' aa (R' + w[1] + '–' + w[2] + '), q ' + fmtP(w[5]) + '</title></rect>';

    // Backbone + domains / repeats / zinc fingers
    label('Sites', top + 10);
    svg += '<rect class="backbone" x="' + clampX(1).toFixed(1) + '" y="' + (base - 4) + '" width="' + (clampX(len) - clampX(1)).toFixed(1) + '" height="8" rx="4"/>';
    const domainLike = feats.filter((f) => ['Domain', 'Repeat', 'Zinc finger'].includes(f.type) && visible(f.start, f.end));
    for (const f of domainLike) {
      const x0 = clampX(f.start), x1 = clampX(f.end), w = Math.max(2, x1 - x0);
      const h = f.type === 'Repeat' ? 12 : 18;
      svg += '<g class="feat" data-tip="' + esc('<b>' + f.type + '</b> ' + esc(f.desc) + '<br>' + f.start + '–' + f.end) + '"><rect class="feat-domain" x="' + x0.toFixed(1) + '" y="' + (base - h / 2) + '" width="' + w.toFixed(1) + '" height="' + h + '" rx="4"' + (f.type === 'Repeat' ? ' opacity=".7"' : '') + '/>' +
        (w > f.desc.length * 6 + 8 && f.type !== 'Repeat' ? '<text class="feat-domain-label" x="' + ((x0 + x1) / 2).toFixed(1) + '" y="' + (base + 4) + '" text-anchor="middle">' + esc(f.desc) + '</text>' : '') + '</g>';
    }
    // Lollipops
    for (const h of heads) {
      const cy = base - 20 - h.lv * 12;
      const dim = !R.shown.has(h.s.j);
      svg += '<line class="stem" x1="' + h.cx.toFixed(1) + '" x2="' + h.cx.toFixed(1) + '" y1="' + (base - 4) + '" y2="' + cy + '"/>';
      svg += '<circle class="head' + (dim ? ' dim' : '') + (UI.highlight === h.s.pos ? ' sel' : '') + '" data-j="' + h.s.j + '" data-pos="' + h.s.pos + '" cx="' + h.cx.toFixed(1) + '" cy="' + cy + '" r="' + R0 + '" fill="var(--st-' + ['none', '1', '2', 'mix'][h.s.st] + ')"/>';
    }
    y = base + 22;

    // Regions: other annotated regions and motifs (with labels when they fit)
    const regions = feats.filter((f) => (f.type === 'Region' && f.desc !== 'Disordered') || f.type === 'Motif').filter((f) => visible(f.start, f.end));
    if (regions.length) {
      label('Regions', y);
      for (const f of regions) {
        const x0 = clampX(f.start), x1 = clampX(f.end), w = Math.max(2, x1 - x0);
        svg += '<g class="feat" data-tip="' + esc('<b>' + f.type + '</b> ' + esc(f.desc) + '<br>' + f.start + '–' + f.end) + '"><rect class="feat-bias" x="' + x0.toFixed(1) + '" y="' + (y - 5) + '" width="' + w.toFixed(1) + '" height="10" rx="3"/>' +
          (w > f.desc.length * 5.6 + 6 ? '<text x="' + (x0 + 3).toFixed(1) + '" y="' + (y + 3.5) + '" font-size="9.5">' + esc(f.desc) + '</text>' : '') + '</g>';
      }
      y += 20;
    }
    // Disorder (UniProt, MobiDB-lite derived)
    const dis = feats.filter((f) => f.type === 'Region' && f.desc === 'Disordered' && visible(f.start, f.end));
    if (u) {
      label('Disorder', y);
      svg += '<line class="grid" x1="' + ML + '" x2="' + (W - MR) + '" y1="' + y + '" y2="' + y + '"/>';
      for (const f of dis) {
        const x0 = clampX(f.start), x1 = clampX(f.end);
        svg += '<g class="feat" data-tip="' + esc('<b>Disordered</b> (UniProt)<br>' + f.start + '–' + f.end) + '"><rect class="feat-region" x="' + x0.toFixed(1) + '" y="' + (y - 5) + '" width="' + Math.max(2, x1 - x0).toFixed(1) + '" height="10" rx="3"/></g>';
      }
      y += 20;
    }
    // Compositional bias
    const bias = feats.filter((f) => f.type === 'Compositional bias' && visible(f.start, f.end));
    if (bias.length) {
      label('Bias', y);
      for (const f of bias) {
        const x0 = clampX(f.start), x1 = clampX(f.end), w = Math.max(2, x1 - x0);
        const short = f.desc.replace(/ residues$/, '-rich').replace('Basic and acidic', 'Charged');
        svg += '<g class="feat" data-tip="' + esc('<b>Compositional bias</b> ' + esc(f.desc) + '<br>' + f.start + '–' + f.end) + '"><rect class="feat-bias" x="' + x0.toFixed(1) + '" y="' + (y - 5) + '" width="' + w.toFixed(1) + '" height="10" rx="3"/>' +
          (w > short.length * 5.6 + 6 ? '<text x="' + (x0 + 3).toFixed(1) + '" y="' + (y + 3.5) + '" font-size="9.5">' + esc(short) + '</text>' : '') + '</g>';
      }
      y += 20;
    }
    // All arginines (from sequence)
    if (u && u.seq) {
      label('Arginines', y);
      const methyl = new Set(p.sites.map((j) => S[j].pos));
      for (let k = Math.max(0, v0 - 1); k < Math.min(u.seq.length, v1); k++) {
        if (u.seq[k] !== 'R' || methyl.has(k + 1)) continue;
        const xx = x(k + 1).toFixed(1);
        svg += '<line class="rtick" x1="' + xx + '" x2="' + xx + '" y1="' + (y - 5) + '" y2="' + (y + 5) + '"/>';
      }
      for (const s of sites) svg += '<line x1="' + x(s.pos).toFixed(1) + '" x2="' + x(s.pos).toFixed(1) + '" y1="' + (y - 7) + '" y2="' + (y + 7) + '" stroke="var(--st-' + ['none', '1', '2', 'mix'][s.st] + ')" stroke-width="2"/>';
      y += 20;
    }
    // Other PTMs
    for (let g = 0; g < 4; g++) {
      const list = p.ptms.filter((t) => t[1] === g && t[0] >= v0 && t[0] <= v1);
      if (!list.length) continue;
      label(PTM_LABEL[g], y);
      svg += '<line class="grid" x1="' + ML + '" x2="' + (W - MR) + '" y1="' + y + '" y2="' + y + '"/>';
      for (const t of list) svg += '<g class="feat" data-tip="' + esc('<b>' + PTM_PREFIX[g] + t[2] + t[0] + '</b> · ' + PTM_LABEL[g]) + '"><rect class="ptm" x="' + (x(t[0]) - 1.5).toFixed(1) + '" y="' + (y - 5) + '" width="3" height="10" rx="1"/></g>';
      y += 18;
    }
    // Axis
    y += 4;
    svg += '<line class="grid" x1="' + ML + '" x2="' + (W - MR) + '" y1="' + y + '" y2="' + y + '"/>';
    const span = v1 - v0;
    const raw = span / 8;
    const mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const step = [1, 2, 5, 10].map((m) => m * mag).find((s) => s >= raw) || mag * 10;
    for (let t = Math.ceil(v0 / step) * step; t <= v1; t += step) {
      svg += '<line class="grid" x1="' + x(t).toFixed(1) + '" x2="' + x(t).toFixed(1) + '" y1="' + y + '" y2="' + (y + 4) + '"/><text class="tick" x="' + x(t).toFixed(1) + '" y="' + (y + 15) + '" text-anchor="middle">' + fmt(t) + '</text>';
    }
    label('aa', y + 11);
    const H = y + 22;
    box.innerHTML = '<svg class="pmap" viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="Map of ' + esc(p.label) + ' with ' + sites.length + ' methylated arginines">' + svg + '<rect class="brush" id="brush" x="0" y="0" width="0" height="' + H + '" hidden/></svg>';
    bindMap(box, p, u, { W, ML, MR, v0, v1, H });
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
      h.addEventListener('click', () => { UI.highlight = s.pos; flashSite(s.pos, true); });
    });
    $$('.feat', svg).forEach((f) => {
      f.addEventListener('mouseenter', (e) => showTip(e, f.dataset.tip));
      f.addEventListener('mousemove', moveTip);
      f.addEventListener('mouseleave', hideTip);
    });
    const reset = $('#zoom-reset');
    reset.onclick = () => { UI.zoom = null; drawMap(p, u); };
  }

  function siteTip(s) {
    const ipr = s.ipr >= 0 ? M.interpro[s.ipr][1] : null;
    return '<b>R' + s.pos + '</b> · ' + STATE_LABEL[s.st] + (s.ip ? ' · ' + bits(s.ip).map((b) => M.maron_ips[b]).join('/') + ' IP' : '') + '<br>' +
      (s.win ? windowHtml(s.win, 7) + '<br>' : '') +
      esc(DISORDER_LABEL[M.disorder_classes[s.dc]] || '') + ' · ' + esc(DOMAIN_LABEL[M.domain_contexts[s.dctx]] || '') + (ipr ? ' (' + esc(ipr) + ')' : '') + '<br>' +
      'Reported by ' + M.studyOrder.filter((k) => s.study & (1 << k)).map((k) => esc(M.studies[k][0])).join(', ') +
      (s.trf ? '<br>' + [0, 1].map((t) => (s.trf & (1 << (t * 2)) ? '▼' + M.maron_treatments[t] + ' ' : '') + (s.trf & (1 << (t * 2 + 1)) ? '▲' + M.maron_treatments[t] + ' ' : '')).join('') : '');
  }

  // ------------------------------------------------------------------ tooltip
  const tip = () => $('#tip');
  function showTip(e, html) { const t = tip(); t.innerHTML = html; t.hidden = false; moveTip(e); }
  function moveTip(e) {
    const t = tip();
    if (t.hidden) return;
    const pad = 14, w = t.offsetWidth, h = t.offsetHeight;
    let left = e.clientX + pad, topY = e.clientY + pad;
    if (left + w > window.innerWidth - 8) left = e.clientX - w - pad;
    if (topY + h > window.innerHeight - 8) topY = e.clientY - h - pad;
    t.style.left = Math.max(8, left) + 'px'; t.style.top = Math.max(8, topY) + 'px';
  }
  function hideTip() { const t = tip(); if (t) t.hidden = true; }

  // ---------------------------------------------------------------- about
  function renderAbout(el) {
    el.innerHTML = '<div class="card prose">' +
      '<h3>About this resource</h3>' +
      '<p>The explorer indexes the integrated human methylarginine union assembled by the <a href="https://github.com/Shechterlab/PTM_Pipeline">PTM_Pipeline</a> for the review below. It runs entirely in your browser; only protein features and sequences are fetched live from UniProt.</p>' +
      '<h3>How to cite</h3><p>' + esc(PAPER.cite) + ' <a href="https://doi.org/' + PAPER.doi + '">doi:' + PAPER.doi + '</a></p>' +
      '<p>Please also cite the original studies for any sites you use; each site\'s evidence lists them.</p>' +
      '<h3>Sources in the union</h3><dl>' +
      '<dt>Maron 2021</dt><dd>Table S5 of Maron et al., <a href="https://doi.org/10.1016/j.isci.2021.102971">iScience 24:102971 (2021)</a>, a compilation of human methylarginine sites from published proteomics.</dd>' +
      '<dt>ProMetheusDB</dt><dd>Methyl-proteome table (mmc2) with monomethyl (Rme1) and dimethyl (Rme2) calls.</dd>' +
      '<dt>PTM master</dt><dd>Methylarginine rows of the lab\'s cross-PTM human master table (UniProt and IEDB annotations); the same table supplies the phospho, lysine acyl, ubiquitin/SUMO and lysine methyl sites shown as crosstalk.</dd></dl>' +
      '<p>Studies and databases behind the sites (sites reported · reported only there):</p><ul>' +
      M.studyOrder.map((k) => '<li><b>' + esc(M.studies[k][0]) + '</b> · ' + fmt(M.studies[k][2]) + ' · ' + fmt(M.studies[k][3]) + '<br><span class="muted">' + esc(M.studies[k][1] && M.studies[k][1] !== 'This study.' ? M.studies[k][1] : STUDY_HELP[M.studies[k][0]] || '') + '</span></li>').join('') + '</ul>' +
      '<h3>How sites were processed</h3><ul>' +
      '<li>Source rows were remapped to canonical UniProt sequences by peptide/sequence anchors and deduplicated by canonical accession and position. The original reported accession and position are kept in each site\'s evidence.</li>' +
      '<li><b>State</b> is the consensus Rme1/Rme2 call where a source provides one; most compiled sites have no state call. Maron Table S2 me1/me2 calls and Table S3 methyl-IP (MMA, ADMA, SDMA) recovery are shown separately.</li>' +
      '<li><b>PRMT inhibitor response</b> uses Maron Table S3 <code>GSK591_vs_DMSO</code> and <code>MS023_vs_DMSO</code> differences (median over peptides) with adjusted p ≤ 0.05. GSK591 inhibits PRMT5; MS023 inhibits type I PRMTs.</li>' +
      '<li><b>Motifs</b> are read from a ±7 window on the canonical sequence: ' + M.motifs.map((m) => MOTIF_LABEL[m] + ' (' + MOTIF_HELP[m] + ')').join('; ') + '.</li>' +
      '<li><b>Disorder</b> uses MobiDB intervals: ' + Object.keys(DISORDER_LABEL).map((k) => DISORDER_LABEL[k] + ' = ' + DISORDER_HELP[k].toLowerCase()).join('; ') + '.</li>' +
      '<li><b>Domain context</b> uses InterPro domains and repeats: ' + Object.keys(DOMAIN_LABEL).map((k) => DOMAIN_LABEL[k] + ' = ' + DOMAIN_HELP[k].toLowerCase()).join('; ') + '.</li>' +
      '<li><b>Clustering</b>: nearest methyl-R spacing and maximum sites in 3–50 aa windows per protein; dense windows come from a multiscale rolling-window scan against a same-protein arginine null (1,000 permutations, q ≤ 0.05). Architecture classes describe proteins with ≥3 sites in a 20 aa span.</li>' +
      '<li><b>Arg enrichment</b> is the per-protein shrinkage log2 enrichment of methylated arginines relative to the proteome-wide rate, with a Fisher test q value.</li>' +
      '<li><b>Condensates</b> are CD-CODE memberships; condensate-level enrichment uses condensates with at least five members.</li>' +
      '<li><b>Disease</b> associations list DisGeNET terms significantly enriched among methylarginine proteins (clusterProfiler) that include the gene, grouped into the audited disease families. The ALS/MN/FTD option uses the full audited family gene set.</li>' +
      '</ul>' +
      '<h3>Interpretive guardrails</h3><ul>' +
      '<li>The integrated union is the primary analysis universe; source support is provenance, not a confidence tier.</li>' +
      '<li>State labels are useful where a source supports them, but the union is not a clean Rme1 / Rme2a / Rme2s state-resolved dataset.</li>' +
      '<li>Domain adjacency is descriptive: an IDR-length-conditioned null shows that simple domain-edge proximity is largely explained by IDR geometry and length.</li>' +
      '<li>The stronger architectural signal is short-range methylarginine clustering relative to same-protein arginine positions.</li>' +
      '<li>Motif and clustering annotations need sequence-confirmed canonical residues; ' + fmt(M.window_flags.missing) + ' sites lack a usable sequence window and ' + fmt(M.window_flags.center_not_r) + ' have a non-R canonical residue (flagged with !).</li>' +
      '</ul>' +
      '<h3>Downloads and updates</h3>' +
      '<p>Use <b>Download TSV</b> on the Proteins or Sites view to export the current selection. The data files are regenerated from pipeline outputs with <code>python scripts/36_export_site_explorer_data.py</code>; see the repository README.</p>' +
      '<p class="note">Data exported ' + esc(M.generated) + ' · ' + fmt(M.n_sites) + ' sites · ' + fmt(M.n_proteins) + ' proteins.</p>' +
      '</div>';
  }

  // ----------------------------------------------------------- export / copy
  function siteTsv(rows) {
    const head = ['gene', 'uniprot_accession', 'protein_name', 'site', 'position', 'state', 'maron_s2_state', 'maron_ip', 'sources', 'source_count',
      'window_15', 'window_31', 'motifs', 'disorder_context', 'disorder_edge_distance', 'domain_context', 'nearest_interpro', 'nearest_domain_class',
      'domain_edge_distance', 'domain_side', 'nearest_methyl_r_distance', 'GSK591_decreased', 'GSK591_increased', 'MS023_decreased', 'MS023_increased',
      'maron_s3_measurements', 'condensates', 'functional_classes', 'reported_by', 'independent_reports', 'references'];
    const out = [head.join('\t')];
    for (const s of rows) {
      const p = P[s.p];
      const ipr = s.ipr >= 0 ? M.interpro[s.ipr] : null;
      out.push([
        p.gene, p.acc, p.name, 'R' + s.pos, s.pos, STATE_LABEL[s.st], bits(s.s2).map((b) => ['me1', 'me2'][b]).join(';'),
        bits(s.ip).map((b) => M.maron_ips[b]).join(';'), bits(s.src).map((b) => M.sources[b]).join(';'), popcount(s.src),
        s.win ? s.win.slice(8, 23) : '', s.win, bits(s.mot).map((b) => M.motifs[b]).join(';'),
        M.disorder_classes[s.dc], s.dEdge == null ? '' : s.dEdge, M.domain_contexts[s.dctx], ipr ? ipr[0] + ' ' + ipr[1] : '',
        M.domain_classes[s.dcls], s.domDist == null ? '' : s.domDist, s.domSide, s.nn == null ? '' : s.nn,
        s.trf & 1 ? 1 : 0, s.trf & 2 ? 1 : 0, s.trf & 4 ? 1 : 0, s.trf & 8 ? 1 : 0,
        s.tr.map(([ip, t, d, q]) => M.maron_ips[ip] + ':' + M.maron_treatments[t] + ':' + d + ':' + q).join(';'),
        p.cond.map((k) => M.condensates[k]).join(';'), bits(p.func).map((b) => M.func_classes[b]).join(';'),
        bits(s.study).map((b) => M.studies[b][0]).join(';'), popcount(s.study),
        s.refs.map((k) => M.refs[k]).join(' | '),
      ].map(cell).join('\t'));
    }
    return out.join('\n') + '\n';
  }

  function proteinTsv(rows) {
    const head = ['gene', 'uniprot_accession', 'uniprot_id', 'protein_name', 'methyl_r_sites_shown', 'methyl_r_sites_total', 'length', 'arginines',
      'shrinkage_log2_enrichment', 'arg_odds_q', 'sites_in_idr', 'condensates', 'disease_families', 'als_ftd_audited_gene', 'functional_classes',
      'best_window_aa', 'best_window_sites', 'best_window_start', 'best_window_end', 'best_window_q', 'architecture_class'];
    const out = [head.join('\t')];
    for (const p of rows) {
      out.push([
        p.gene, p.acc, p.uid, p.name, R.counts[p.i], p.sites.length, p.len == null ? '' : p.len, p.argN == null ? '' : p.argN,
        p.shrink == null ? '' : p.shrink, p.argQ == null ? '' : p.argQ, p.nIdr, p.cond.map((k) => M.condensates[k]).join(';'),
        bits(p.dis).map((b) => M.disease_families[b]).join(';'), p.als ? 1 : 0, bits(p.func).map((b) => M.func_classes[b]).join(';'),
        p.best ? p.best[0] : '', p.best ? p.best[3] : '', p.best ? p.best[1] : '', p.best ? p.best[2] : '', p.best ? p.best[5] : '',
        p.arch ? p.arch[0] : '',
      ].map(cell).join('\t'));
    }
    return out.join('\n') + '\n';
  }

  const cell = (v) => String(v == null ? '' : v).replace(/[\t\n\r]+/g, ' ');

  function download(name, text) {
    const blob = new Blob([text], { type: 'text/tab-separated-values;charset=utf-8' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url; a.download = name;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 2000);
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

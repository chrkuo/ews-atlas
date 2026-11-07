/* Single-Cell Atlas viewer — static, no build step.
 * Reads data/datasets.json and data/<id>/... written by tools/export_h5ad.py.
 * URL hash keeps state: #dataset=tumor&color=cell_type&gene=CD68&split=timepoint
 */
(() => {
  "use strict";
  const $ = (s) => document.querySelector(s);
  const DATA = "data";

  const S = {
    registry: null, ds: null, coords: {}, obs: {}, order: null,
    emb: null, colorField: null, gene: null, geneVals: null,
    split: null, hidden: new Set(), ptSize: 2,
    view: { k: 1, cx: 0, cy: 0 }, bounds: null,
    geneCache: new Map(), panels: [],
  };

  // ---------- utils
  const css = (v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
  const hexRgb = (h) => {
    h = h.replace("#", "");
    if (h.length === 3) h = h.split("").map((c) => c + c).join("");
    const n = parseInt(h, 16);
    return [(n >> 16) & 255, (n >> 8) & 255, n & 255];
  };
  const pack = ([r, g, b]) => ((255 << 24) | (b << 16) | (g << 8) | r) >>> 0;
  const isDark = () => {
    const [r, g, b] = hexRgb(css("--surface") || "#fff");
    return r + g + b < 300;
  };
  const RAMP_LIGHT = ["#f6e8b1", "#f2b14f", "#de6a3a", "#a3254f", "#3b0f4f"];
  const RAMP_DARK = ["#33305c", "#6a3c8f", "#c34f72", "#f39654", "#fbf2a6"];
  function rampLUT() {
    const stops = (isDark() ? RAMP_DARK : RAMP_LIGHT).map(hexRgb);
    const lut = new Uint32Array(256), rgb = [];
    for (let i = 0; i < 256; i++) {
      const t = (i / 255) * (stops.length - 1), j = Math.min(Math.floor(t), stops.length - 2), f = t - j;
      const c = [0, 1, 2].map((k) => Math.round(stops[j][k] + (stops[j + 1][k] - stops[j][k]) * f));
      lut[i] = pack(c); rgb.push(c);
    }
    return { lut, rgb };
  }
  const fmt = (n) => n.toLocaleString();
  const el = (tag, attrs = {}, ...kids) => {
    const e = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs)) {
      if (k === "class") e.className = v; else if (k.startsWith("on")) e.addEventListener(k.slice(2), v);
      else e.setAttribute(k, v);
    }
    kids.flat().forEach((k) => e.append(k));
    return e;
  };
  async function getJSON(p) { const r = await fetch(p); if (!r.ok) throw new Error(p + " " + r.status); return r.json(); }
  async function getBin(p) { const r = await fetch(p); if (!r.ok) throw new Error(p + " " + r.status); return r.arrayBuffer(); }

  // ---------- hash state
  function readHash() {
    const p = new URLSearchParams(location.hash.slice(1));
    return Object.fromEntries(p.entries());
  }
  function writeHash() {
    const p = new URLSearchParams();
    p.set("dataset", S.ds.id);
    if (S.emb && S.emb !== S.ds.embeddings[0].name) p.set("emb", S.emb);
    if (S.gene) p.set("gene", S.gene); else p.set("color", S.colorField);
    if (S.split) p.set("split", S.split);
    const h = "#" + p.toString();
    if (h !== location.hash) history.replaceState(null, "", h);
  }

  // ---------- loading
  async function init() {
    let saved = null; try { saved = localStorage.getItem("atlas-theme"); } catch (_) {}
    applyTheme(saved);
    S.registry = await getJSON(`${DATA}/datasets.json`);
    if (S.registry.title) { $("#site-title").textContent = S.registry.title; document.title = S.registry.title; }
    $("#footer").innerHTML = S.registry.footer || "";
    const tabs = $("#dataset-tabs");
    S.registry.datasets.forEach((d) =>
      tabs.append(el("button", { class: "tab", "data-id": d.id, onclick: () => switchDataset(d.id) }, d.name)));
    bindUI();
    const h = readHash();
    await switchDataset(h.dataset || S.registry.datasets[0].id, h);
    window.addEventListener("hashchange", () => {
      const h2 = readHash();
      if (h2.dataset && h2.dataset !== S.ds.id) switchDataset(h2.dataset, h2); else applyState(h2);
    });
  }

  async function switchDataset(id, state = {}) {
    if (!S.registry.datasets.find((d) => d.id === id)) id = S.registry.datasets[0].id;
    $("#loading").hidden = false;
    document.querySelectorAll("#dataset-tabs .tab").forEach((t) => t.classList.toggle("active", t.dataset.id === id));
    const ds = await getJSON(`${DATA}/${id}/meta.json`);
    S.ds = ds; S.coords = {}; S.obs = {}; S.geneCache.clear(); S.hidden.clear();
    S.geneIndex = new Map(ds.genes.map((g, i) => [g.toUpperCase(), i]));
    await Promise.all(ds.embeddings.map(async (e) => { S.coords[e.name] = new Float32Array(await getBin(`${DATA}/${id}/${e.file}`)); }));
    await Promise.all(ds.fields.map(async (f) => {
      const buf = await getBin(`${DATA}/${id}/${f.file}`);
      S.obs[f.name] = f.dtype === "uint16" ? new Uint16Array(buf) : new Uint8Array(buf);
    }));
    // shuffled draw order so no category systematically sits on top
    const n = ds.n_cells, ord = new Uint32Array(n);
    for (let i = 0; i < n; i++) ord[i] = i;
    let seed = 12345;
    for (let i = n - 1; i > 0; i--) { seed = (seed * 1103515245 + 12345) & 0x7fffffff; const j = seed % (i + 1); [ord[i], ord[j]] = [ord[j], ord[i]]; }
    S.order = ord;

    $("#ds-info").innerHTML = "";
    $("#ds-info").append(el("h2", {}, ds.name), el("p", {}, ds.description || ""),
      el("div", { class: "stat" }, `${fmt(ds.n_cells)} cells · ${fmt(ds.genes.length)} genes`));
    const fillSel = (sel, withNone) => {
      sel.innerHTML = ""; if (withNone) sel.append(el("option", { value: "" }, "—"));
      ds.fields.forEach((f) => sel.append(el("option", { value: f.name }, f.name)));
    };
    fillSel($("#color-field")); fillSel($("#split-field")); fillSel($("#dot-group")); fillSel($("#violin-group"));
    const embSel = $("#embedding"); embSel.innerHTML = "";
    ds.embeddings.forEach((e) => embSel.append(el("option", { value: e.name }, e.name)));
    embSel.hidden = ds.embeddings.length < 2;
    $("#dot-group").value = ds.default_color; $("#violin-group").value = ds.default_color;
    $("#dot-plot").innerHTML = `<div class="empty">Enter genes and press Plot.</div>`;
    S.dotLast = null;
    await applyState({ color: ds.default_color, ...state }, true);
    $("#loading").hidden = true;
  }

  async function applyState(st, fresh = false) {
    const ds = S.ds;
    S.emb = ds.embeddings.find((e) => e.name === st.emb) ? st.emb : ds.embeddings[0].name;
    $("#embedding").value = S.emb;
    S.colorField = ds.fields.find((f) => f.name === st.color) ? st.color : (S.colorField && ds.fields.find((f) => f.name === S.colorField) ? S.colorField : ds.default_color);
    $("#color-field").value = S.colorField;
    S.split = ds.fields.find((f) => f.name === st.split) ? st.split : null;
    $("#split-toggle").checked = !!S.split; $("#split-field").disabled = !S.split;
    if (S.split) $("#split-field").value = S.split;
    if (fresh || st.emb) fitBounds();
    if (st.gene) await setGene(st.gene, false); else { S.gene = null; S.geneVals = null; $("#gene-input").value = ""; }
    layoutPanels(); renderAll(); renderLegend(); renderViolin(); writeHash();
  }

  async function loadGene(name) {
    const i = S.geneIndex.get(name.toUpperCase());
    if (i === undefined) return null;
    if (S.geneCache.has(i)) return S.geneCache.get(i);
    const buf = await getBin(`${DATA}/${S.ds.id}/genes/${Math.floor(i / 1000)}/${i}.bin`);
    const dv = new DataView(buf), nnz = dv.getUint32(0, true), max = dv.getFloat32(4, true);
    const idx = new Uint32Array(buf, 8, nnz), q = new Uint8Array(buf, 8 + 4 * nnz, nnz);
    const vals = new Float32Array(S.ds.n_cells);
    for (let j = 0; j < nnz; j++) vals[idx[j]] = (q[j] / 255) * max;
    const g = { name: S.ds.genes[i], vals, max, nnz };
    S.geneCache.set(i, g);
    if (S.geneCache.size > 60) S.geneCache.delete(S.geneCache.keys().next().value);
    return g;
  }

  async function setGene(name, redraw = true) {
    const g = name ? await loadGene(name) : null;
    if (name && !g) { flash(`“${name}” is not in this dataset`); return; }
    S.gene = g ? g.name : null; S.geneVals = g;
    $("#gene-input").value = S.gene || "";
    if (redraw) { renderAll(); renderLegend(); renderViolin(); writeHash(); }
  }

  function flash(msg) {
    const t = $("#tooltip"); t.textContent = msg; t.hidden = false; t.style.left = "12px"; t.style.top = "12px";
    clearTimeout(flash._t); flash._t = setTimeout(() => (t.hidden = true), 2200);
  }

  // ---------- geometry
  function fitBounds() {
    const xy = S.coords[S.emb]; let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    for (let i = 0; i < xy.length; i += 2) {
      const x = xy[i], y = xy[i + 1];
      if (x < x0) x0 = x; if (x > x1) x1 = x; if (y < y0) y0 = y; if (y > y1) y1 = y;
    }
    S.bounds = { x0, x1, y0, y1 };
    S.view = { k: 1, cx: (x0 + x1) / 2, cy: (y0 + y1) / 2 };
  }
  // data -> CSS px within a panel of size w,h
  function scaleFor(w, h) {
    const b = S.bounds, pad = 0.06;
    return Math.min(w / ((b.x1 - b.x0) * (1 + 2 * pad)), h / ((b.y1 - b.y0) * (1 + 2 * pad))) * S.view.k;
  }

  // ---------- panels
  function layoutPanels() {
    const grid = $("#grid"); grid.innerHTML = ""; S.panels = [];
    let cats = [null];
    if (S.split) { const f = S.ds.fields.find((x) => x.name === S.split); cats = f.categories.map((c, i) => ({ label: c, code: i })); }
    const n = cats.length, W = grid.clientWidth, H = grid.clientHeight;
    let cols = 1;
    if (n > 1) { let best = 0; for (let c = 1; c <= n; c++) { const r = Math.ceil(n / c); const s = Math.min(W / c, H / r); if (s > best) { best = s; cols = c; } } }
    grid.style.gridTemplateColumns = `repeat(${cols}, 1fr)`;
    grid.style.gridTemplateRows = `repeat(${Math.ceil(n / cols)}, 1fr)`;
    cats.forEach((c) => {
      const cv = el("canvas"); const p = el("div", { class: "cell-panel" }, cv);
      if (c) p.append(el("div", { class: "ptitle" }, c.label));
      grid.append(p); S.panels.push({ cv, split: c }); bindCanvas(cv);
    });
    updateTitle();
  }

  function updateTitle() {
    const t = $("#plot-title"); t.innerHTML = "";
    t.append(S.gene ? S.gene : S.colorField);
    t.append(el("span", { class: "sub" }, `${S.emb}${S.split ? " · split by " + S.split : ""}`));
  }

  let rafPending = false;
  function renderAll() {
    if (rafPending) return; rafPending = true;
    requestAnimationFrame(() => { rafPending = false; S.panels.forEach(drawPanel); updateTitle(); });
  }

  function discOffsets(r) {
    const out = [], R = Math.max(0.5, r);
    for (let dy = -Math.ceil(R); dy <= Math.ceil(R); dy++)
      for (let dx = -Math.ceil(R); dx <= Math.ceil(R); dx++)
        if (dx * dx + dy * dy <= R * R + 0.25) out.push([dx, dy]);
    return out.length ? out : [[0, 0]];
  }

  function drawPanel(p) {
    const cv = p.cv, dpr = window.devicePixelRatio || 1;
    const w = cv.clientWidth, h = cv.clientHeight; if (!w || !h) return;
    cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr);
    const ctx = cv.getContext("2d"), img = ctx.createImageData(cv.width, cv.height);
    const px = new Uint32Array(img.data.buffer), W = cv.width, H = cv.height;
    px.fill(pack(hexRgb(css("--surface"))));
    const bg = pack(hexRgb(css("--cell-bg")));
    const sc = scaleFor(w, h) * dpr, ox = W / 2, oy = H / 2, cx = S.view.cx, cy = S.view.cy;
    const xy = S.coords[S.emb], ord = S.order, n = S.ds.n_cells;
    const offs = discOffsets((S.ptSize * dpr) / 2);
    const put = (i, col) => {
      const x = Math.round(ox + (xy[2 * i] - cx) * sc), y = Math.round(oy - (xy[2 * i + 1] - cy) * sc);
      if (x < -4 || y < -4 || x > W + 4 || y > H + 4) return;
      for (const [dx, dy] of offs) { const X = x + dx, Y = y + dy; if (X >= 0 && Y >= 0 && X < W && Y < H) px[Y * W + X] = col; }
    };
    const splitCodes = p.split ? S.obs[S.split] : null, sc0 = p.split ? p.split.code : -1;
    const inPanel = (i) => !splitCodes || splitCodes[i] === sc0;

    if (S.geneVals) {
      const v = S.geneVals.vals, max = S.geneVals.max || 1, { lut } = rampLUT();
      const fg = [];
      for (let k = 0; k < n; k++) { const i = ord[k]; if (inPanel(i) && v[i] > 0) fg.push(i); else put(i, bg); }
      fg.sort((a, b) => v[a] - v[b]);
      for (const i of fg) put(i, lut[Math.min(255, Math.round((v[i] / max) * 255))]);
    } else {
      const f = S.ds.fields.find((x) => x.name === S.colorField), codes = S.obs[f.name];
      const cols = f.colors.map((c) => pack(hexRgb(c)));
      const fg = [];
      for (let k = 0; k < n; k++) {
        const i = ord[k];
        if (inPanel(i) && !S.hidden.has(codes[i]) && codes[i] < cols.length) fg.push(i); else put(i, bg);
      }
      for (const i of fg) put(i, cols[codes[i]]);
    }
    ctx.putImageData(img, 0, 0);
  }

  // ---------- interaction
  function bindCanvas(cv) {
    let drag = null;
    cv.addEventListener("pointerdown", (e) => { drag = { x: e.clientX, y: e.clientY, cx: S.view.cx, cy: S.view.cy }; cv.setPointerCapture(e.pointerId); $("#tooltip").hidden = true; });
    cv.addEventListener("pointermove", (e) => {
      if (drag) {
        const sc = scaleFor(cv.clientWidth, cv.clientHeight);
        S.view.cx = drag.cx - (e.clientX - drag.x) / sc; S.view.cy = drag.cy + (e.clientY - drag.y) / sc;
        S.panels.forEach(drawPanel);
      } else hover(e, cv);
    });
    cv.addEventListener("pointerup", () => (drag = null));
    cv.addEventListener("pointerleave", () => { $("#tooltip").hidden = true; });
    cv.addEventListener("wheel", (e) => {
      e.preventDefault();
      zoomAt(cv, e, Math.exp(-e.deltaY * 0.0015));
    }, { passive: false });
    cv.addEventListener("dblclick", (e) => zoomAt(cv, e, 2));
  }
  function zoomAt(cv, e, f) {
    const r = cv.getBoundingClientRect(), w = r.width, h = r.height;
    const sc = scaleFor(w, h), mx = e.clientX - r.left - w / 2, my = e.clientY - r.top - h / 2;
    const dx = S.view.cx + mx / sc, dy = S.view.cy - my / sc;
    S.view.k = Math.min(60, Math.max(0.4, S.view.k * f));
    const sc2 = scaleFor(w, h);
    S.view.cx = dx - mx / sc2; S.view.cy = dy + my / sc2;
    renderAll();
  }
  let hoverRaf = null;
  function hover(e, cv) {
    if (hoverRaf) return;
    hoverRaf = requestAnimationFrame(() => {
      hoverRaf = null;
      const r = cv.getBoundingClientRect(), w = r.width, h = r.height, sc = scaleFor(w, h);
      const mx = e.clientX - r.left, my = e.clientY - r.top;
      const xy = S.coords[S.emb], n = S.ds.n_cells, cx = S.view.cx, cy = S.view.cy;
      let best = -1, bd = 64;
      for (let i = 0; i < n; i++) {
        const x = w / 2 + (xy[2 * i] - cx) * sc - mx, y = h / 2 - (xy[2 * i + 1] - cy) * sc - my, d = x * x + y * y;
        if (d < bd) { bd = d; best = i; }
      }
      const tt = $("#tooltip");
      if (best < 0) { tt.hidden = true; return; }
      const lines = S.ds.fields.map((f) => `${f.name}: ${f.categories[S.obs[f.name][best]] ?? "NA"}`);
      if (S.geneVals) lines.unshift(`${S.gene}: ${S.geneVals.vals[best].toFixed(2)}`);
      tt.innerHTML = lines.map((l) => `<div>${l.replace(/</g, "&lt;")}</div>`).join("");
      const wrap = $("#canvas-wrap").getBoundingClientRect();
      tt.hidden = false;
      let L = e.clientX - wrap.left + 14, T = e.clientY - wrap.top + 14;
      if (L + tt.offsetWidth > wrap.width) L -= tt.offsetWidth + 28;
      if (T + tt.offsetHeight > wrap.height) T -= tt.offsetHeight + 28;
      tt.style.left = L + "px"; tt.style.top = T + "px";
    });
  }

  // ---------- legend
  function renderLegend() {
    const lg = $("#legend"); lg.innerHTML = "";
    if (S.geneVals) {
      $("#legend-title").textContent = S.gene; $("#legend-all").hidden = true;
      const { rgb } = rampLUT();
      const stops = [0, 64, 128, 192, 255].map((i) => `rgb(${rgb[i].join(",")}) ${(i / 255) * 100}%`).join(",");
      const pct = ((S.geneVals.nnz / S.ds.n_cells) * 100).toFixed(1);
      lg.append(el("div", { class: "bar", style: `background:linear-gradient(90deg,${stops})` }),
        el("div", { class: "scale" }, el("span", {}, "0"), el("span", {}, S.geneVals.max.toFixed(2))),
        el("p", { class: "muted", style: "font-size:12px;margin:8px 0 0" },
          `Normalized expression. Detected in ${fmt(S.geneVals.nnz)} cells (${pct}%).`),
        el("button", { class: "link-btn", style: "margin-top:6px", onclick: () => setGene(null) }, "← back to " + S.colorField));
      return;
    }
    const f = S.ds.fields.find((x) => x.name === S.colorField);
    $("#legend-title").textContent = f.name; $("#legend-all").hidden = S.hidden.size === 0;
    f.categories.forEach((c, i) => {
      const it = el("div", { class: "item" + (S.hidden.has(i) ? " off" : ""), title: "Click to hide; double-click to isolate" },
        el("span", { class: "sw", style: `background:${f.colors[i]}` }), el("span", {}, c), el("span", { class: "n" }, fmt(f.counts[i])));
      it.addEventListener("click", () => { S.hidden.has(i) ? S.hidden.delete(i) : S.hidden.add(i); renderLegend(); renderAll(); });
      it.addEventListener("dblclick", () => { S.hidden = new Set(f.categories.map((_, j) => j).filter((j) => j !== i)); renderLegend(); renderAll(); });
      lg.append(it);
    });
  }

  // ---------- gene autocomplete
  function bindGeneSearch() {
    const inp = $("#gene-input"), ul = $("#gene-suggest"); let hi = -1, items = [];
    const show = () => {
      const q = inp.value.trim().toUpperCase(); ul.innerHTML = ""; hi = -1;
      if (!q) { items = []; return; }
      const pre = [], mid = [];
      S.ds.genes.forEach((g, i) => { const G = g.toUpperCase(); if (G.startsWith(q)) pre.push(i); else if (G.includes(q)) mid.push(i); });
      items = pre.concat(mid).slice(0, 30);
      items.forEach((i) => ul.append(el("li", { onmousedown: (e) => { e.preventDefault(); pick(S.ds.genes[i]); } },
        S.ds.genes[i], el("span", {}, fmt(S.ds.gene_nnz[i]) + " cells"))));
    };
    const pick = (g) => { ul.innerHTML = ""; inp.blur(); setGene(g); };
    inp.addEventListener("input", show);
    inp.addEventListener("focus", show);
    inp.addEventListener("blur", () => setTimeout(() => (ul.innerHTML = ""), 100));
    inp.addEventListener("keydown", (e) => {
      const lis = ul.querySelectorAll("li");
      if (e.key === "ArrowDown" || e.key === "ArrowUp") {
        e.preventDefault(); hi = (hi + (e.key === "ArrowDown" ? 1 : -1) + lis.length) % lis.length;
        lis.forEach((l, j) => l.classList.toggle("hi", j === hi));
      } else if (e.key === "Enter") {
        e.preventDefault();
        const typed = inp.value.trim(), exact = S.geneIndex.has(typed.toUpperCase());
        pick(hi >= 0 ? S.ds.genes[items[hi]] : exact || !items.length ? typed : S.ds.genes[items[0]]);
      } else if (e.key === "Escape") { ul.innerHTML = ""; }
    });
  }

  // ---------- charts
  const SVGNS = "http://www.w3.org/2000/svg";
  const sv = (tag, a = {}) => { const e = document.createElementNS(SVGNS, tag); for (const [k, v] of Object.entries(a)) e.setAttribute(k, v); return e; };

  async function renderDot() {
    const box = $("#dot-plot");
    const names = $("#dot-genes").value.split(/[\s,;]+/).filter(Boolean);
    if (!names.length) { box.innerHTML = `<div class="empty">Enter genes and press Plot.</div>`; return; }
    const genes = [], missing = [];
    for (const nm of names) { const g = await loadGene(nm); g ? genes.push(g) : missing.push(nm); }
    const f = S.ds.fields.find((x) => x.name === $("#dot-group").value), codes = S.obs[f.name], K = f.categories.length;
    const stats = genes.map((g) => {
      const sum = new Float64Array(K), pos = new Float64Array(K), cnt = new Float64Array(K);
      for (let i = 0; i < codes.length; i++) { const c = codes[i]; if (c >= K) continue; cnt[c]++; sum[c] += g.vals[i]; if (g.vals[i] > 0) pos[c]++; }
      const mean = Array.from(sum, (s, c) => (cnt[c] ? s / cnt[c] : 0)), frac = Array.from(pos, (p, c) => (cnt[c] ? p / cnt[c] : 0));
      const lo = Math.min(...mean), hi = Math.max(...mean);
      return { mean, frac, scaled: mean.map((m) => (hi > lo ? (m - lo) / (hi - lo) : 0)) };
    });
    S.dotLast = true;
    const labW = Math.max(...f.categories.map((c) => c.length)) * 7 + 16, cell = 26, top = 70;
    const W = labW + genes.length * cell + 210, H = top + K * cell + 20;
    const svg = sv("svg", { width: W, height: Math.max(H, 190), viewBox: `0 0 ${W} ${Math.max(H, 190)}` });
    const { rgb } = rampLUT(), colr = (t) => `rgb(${rgb[Math.round(t * 255)].join(",")})`;
    genes.forEach((g, j) => {
      const x = labW + j * cell + cell / 2;
      const t = sv("text", { x, y: top - 8, transform: `rotate(-50 ${x} ${top - 8})`, "font-style": "italic" }); t.textContent = g.name; svg.append(t);
    });
    f.categories.forEach((c, k) => {
      const y = top + k * cell + cell / 2;
      const t = sv("text", { x: labW - 8, y: y + 4, "text-anchor": "end" }); t.textContent = c; svg.append(t);
      svg.append(sv("line", { x1: labW, x2: labW + genes.length * cell, y1: y, y2: y, class: "ax", "stroke-width": 0.5 }));
      stats.forEach((s, j) => {
        const r = Math.sqrt(s.frac[k]) * (cell / 2 - 1);
        if (r < 0.4) return;
        const dot = sv("circle", { cx: labW + j * cell + cell / 2, cy: y, r, fill: colr(s.scaled[k]), stroke: "var(--ink-2)", "stroke-width": 0.4 });
        const ti = sv("title"); ti.textContent = `${genes[j].name} in ${c}\nmean ${s.mean[k].toFixed(2)} · ${(s.frac[k] * 100).toFixed(1)}% expressing`;
        dot.append(ti); svg.append(dot);
      });
    });
    // legends
    const lx = labW + genes.length * cell + 30;
    const t1 = sv("text", { x: lx, y: top + 4, "font-size": 11 }); t1.textContent = "% expressing"; svg.append(t1);
    [0.25, 0.5, 0.75, 1].forEach((p, i) => {
      svg.append(sv("circle", { cx: lx + 8 + i * 28, cy: top + 22, r: Math.sqrt(p) * (cell / 2 - 1), fill: "var(--ink-2)" }));
      const t = sv("text", { x: lx + 8 + i * 28, y: top + 46, "text-anchor": "middle", "font-size": 10 }); t.textContent = p * 100; svg.append(t);
    });
    const t2 = sv("text", { x: lx, y: top + 76, "font-size": 11 }); t2.textContent = "Mean expr. (scaled per gene)"; svg.append(t2);
    const defs = sv("defs"), lg = sv("linearGradient", { id: "dotramp" });
    [0, 0.25, 0.5, 0.75, 1].forEach((o) => lg.append(sv("stop", { offset: o, "stop-color": colr(o) })));
    defs.append(lg); svg.append(defs);
    svg.append(sv("rect", { x: lx, y: top + 84, width: 112, height: 10, rx: 2, fill: "url(#dotramp)" }));
    const a = sv("text", { x: lx, y: top + 108, "font-size": 10 }); a.textContent = "0"; svg.append(a);
    const b = sv("text", { x: lx + 112, y: top + 108, "font-size": 10, "text-anchor": "end" }); b.textContent = "1"; svg.append(b);
    box.innerHTML = "";
    if (missing.length) box.append(el("div", { class: "empty" }, `Not found: ${missing.join(", ")}`));
    if (genes.length) box.append(svg);
  }

  function renderViolin() {
    const box = $("#violin-plot");
    if (!S.geneVals) { box.innerHTML = `<div class="empty">Search a gene to see its distribution by group.</div>`; return; }
    const f = S.ds.fields.find((x) => x.name === $("#violin-group").value), codes = S.obs[f.name], K = f.categories.length;
    const v = S.geneVals.vals, max = S.geneVals.max || 1, groups = Array.from({ length: K }, () => []);
    for (let i = 0; i < codes.length; i++) if (codes[i] < K) groups[codes[i]].push(v[i]);
    const colW = 56, left = 46, top = 14, plotH = 200, bottom = 90, W = left + K * colW + 10, H = top + plotH + bottom;
    const svg = sv("svg", { width: W, height: H, viewBox: `0 0 ${W} ${H}` });
    const yOf = (val) => top + plotH - (val / max) * plotH;
    [0, 0.5, 1].forEach((t) => {
      svg.append(sv("line", { x1: left, x2: W - 10, y1: yOf(t * max), y2: yOf(t * max), class: "ax", "stroke-width": 0.5 }));
      const tx = sv("text", { x: left - 6, y: yOf(t * max) + 4, "text-anchor": "end", "font-size": 10 }); tx.textContent = (t * max).toFixed(1); svg.append(tx);
    });
    const NB = 48, bw = max / 14;
    groups.forEach((g, k) => {
      const cx = left + k * colW + colW / 2;
      if (g.length) {
        const dens = new Float64Array(NB + 1);
        const sample = g.length > 4000 ? g.filter((_, i) => i % Math.ceil(g.length / 4000) === 0) : g;
        for (const x of sample) for (let b = 0; b <= NB; b++) { const z = (x - (b / NB) * max) / bw; dens[b] += Math.exp(-0.5 * z * z); }
        const dm = Math.max(...dens) || 1;
        let d = "";
        for (let b = 0; b <= NB; b++) d += `${b ? "L" : "M"}${cx + (dens[b] / dm) * (colW / 2 - 6)},${yOf((b / NB) * max)}`;
        for (let b = NB; b >= 0; b--) d += `L${cx - (dens[b] / dm) * (colW / 2 - 6)},${yOf((b / NB) * max)}`;
        const color = S.colorField === f.name ? f.colors[k] : "var(--accent)";
        svg.append(sv("path", { d: d + "Z", fill: color, "fill-opacity": 0.75, stroke: color }));
        const sorted = Float32Array.from(g).sort(), med = sorted[Math.floor(sorted.length / 2)];
        svg.append(sv("line", { x1: cx - 8, x2: cx + 8, y1: yOf(med), y2: yOf(med), stroke: "var(--ink)", "stroke-width": 2 }));
        const pct = (g.filter((x) => x > 0).length / g.length) * 100;
        const pt = sv("text", { x: cx, y: top + plotH + 14, "text-anchor": "middle", "font-size": 10, fill: "var(--ink-2)" }); pt.textContent = pct.toFixed(0) + "%"; svg.append(pt);
      }
      const lx = cx, ly = top + plotH + 26;
      const t = sv("text", { x: lx, y: ly, "text-anchor": "end", transform: `rotate(-40 ${lx} ${ly})`, "font-size": 11 }); t.textContent = f.categories[k]; svg.append(t);
    });
    box.innerHTML = ""; box.append(el("div", { class: "muted", style: "font-size:12px" }, `${S.gene} — bar = median, % = cells with detected expression`), svg);
  }

  // ---------- misc UI
  function applyTheme(t) {
    if (t) document.documentElement.setAttribute("data-theme", t); else document.documentElement.removeAttribute("data-theme");
  }
  function savePNG() {
    const cvs = S.panels.map((p) => p.cv); const first = cvs[0];
    const grid = $("#grid").getBoundingClientRect(), dpr = window.devicePixelRatio || 1;
    const out = document.createElement("canvas"); out.width = grid.width * dpr; out.height = grid.height * dpr;
    const ctx = out.getContext("2d"); ctx.fillStyle = css("--line"); ctx.fillRect(0, 0, out.width, out.height);
    S.panels.forEach((p) => {
      const r = p.cv.getBoundingClientRect(); const x = (r.left - grid.left) * dpr, y = (r.top - grid.top) * dpr;
      ctx.drawImage(p.cv, x, y);
      if (p.split) { ctx.fillStyle = css("--ink"); ctx.font = `${12 * dpr}px IBM Plex Sans, sans-serif`; ctx.fillText(p.split.label, x + 8 * dpr, y + 18 * dpr); }
    });
    const a = document.createElement("a");
    a.download = `${S.ds.id}_${S.gene || S.colorField}${S.split ? "_by_" + S.split : ""}.png`;
    a.href = (cvs.length > 1 ? out : first).toDataURL("image/png"); a.click();
  }

  function bindUI() {
    bindGeneSearch();
    $("#color-field").addEventListener("change", (e) => { S.colorField = e.target.value; S.hidden.clear(); S.gene = null; S.geneVals = null; $("#gene-input").value = ""; renderAll(); renderLegend(); renderViolin(); writeHash(); });
    $("#split-toggle").addEventListener("change", (e) => {
      $("#split-field").disabled = !e.target.checked;
      S.split = e.target.checked ? $("#split-field").value : null; layoutPanels(); renderAll(); writeHash();
    });
    $("#split-field").addEventListener("change", (e) => { S.split = e.target.value; layoutPanels(); renderAll(); writeHash(); });
    $("#pt-size").addEventListener("input", (e) => { S.ptSize = +e.target.value; renderAll(); });
    $("#embedding").addEventListener("change", (e) => { S.emb = e.target.value; fitBounds(); renderAll(); writeHash(); });
    $("#reset-view").addEventListener("click", () => { fitBounds(); renderAll(); });
    $("#save-png").addEventListener("click", savePNG);
    $("#legend-all").addEventListener("click", () => { S.hidden.clear(); renderLegend(); renderAll(); });
    $("#dot-run").addEventListener("click", renderDot);
    $("#dot-genes").addEventListener("keydown", (e) => { if (e.key === "Enter") renderDot(); });
    $("#dot-group").addEventListener("change", () => S.dotLast && renderDot());
    $("#violin-group").addEventListener("change", renderViolin);
    document.querySelectorAll(".lower .tab").forEach((b) => b.addEventListener("click", () => {
      document.querySelectorAll(".lower .tab").forEach((x) => x.classList.toggle("active", x === b));
      $("#tab-dot").hidden = b.dataset.tab !== "dot"; $("#tab-violin").hidden = b.dataset.tab !== "violin";
    }));
    $("#theme-toggle").addEventListener("click", () => {
      const next = isDark() ? "light" : "dark"; applyTheme(next);
      try { localStorage.setItem("atlas-theme", next); } catch (_) {}
      renderAll(); renderLegend(); renderViolin(); if (S.dotLast) renderDot();
    });
    new ResizeObserver(() => { if (S.ds) { layoutPanels(); renderAll(); } }).observe($("#canvas-wrap"));
  }

  init().catch((err) => { $("#loading").textContent = "Failed to load data: " + err.message; console.error(err); });
})();

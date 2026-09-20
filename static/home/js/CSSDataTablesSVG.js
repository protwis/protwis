/**
 * DataTablesSVG.js
 * ----------------
 * Vector export of the CURRENT table DOM, plus an optional composer that adds
 * the axis badges in SVG. Handles 3-line headers, clamps with ellipsis, and clips text.
 *
 * API:
 *   // Heatmap grid (your main matrix)
 *   DataTablesSVG.export('#ccsTable','ccs-table', { filename:'My.svg' })
 *   DataTablesSVG.exportWithBadges('.ccs-axis-wrap', '#ccsTable', 'ccs-heatmap', {
 *     padTop:60, padLeft:60, gapX:8, gapY:10, badgeTighten:40, background:null, filename:'MyWithBadges.svg'
 *   })
 *
 *   // NEW: generic/flyout table (multi-row headers, colspans/rowspans)
 *   DataTablesSVG.exportFlyout('#ccsFlyTable','ccs-fly', { filename:'Flyout.svg' })
 */

(function (global) {
  const svgNS = 'http://www.w3.org/2000/svg';

  // ---------- filename helpers ----------
  function basenameFromSelector(selector, fallback) {
    if (typeof selector === 'string') {
      const m = selector.match(/#([\w-]+)/);
      if (m && m[1]) return m[1];
    }
    return fallback || 'table';
  }
  function stamped(base, ext) {
    const d = new Date(), pad = n => String(n).padStart(2, '0');
    const stamp = `${d.getFullYear()}-${pad(d.getMonth()+1)}-${pad(d.getDate())}_` +
                  `${pad(d.getHours())}-${pad(d.getMinutes())}-${pad(d.getSeconds())}`;
    return `${base}-${stamp}.${ext}`;
  }
  function resolveFilename(selectorOrBase, ext, explicit, fallback) {
    if (explicit) return explicit;
    if (selectorOrBase && selectorOrBase[0] === '#') {
      return stamped(basenameFromSelector(selectorOrBase, fallback), ext);
    }
    const base = (selectorOrBase && String(selectorOrBase).trim()) || fallback || 'table';
    return stamped(base, ext);
  }
  function triggerDownload(text, filename, type='image/svg+xml') {
    const blob = new Blob([text], { type });
    const url  = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url; a.download = filename;
    document.body.appendChild(a); a.click(); a.remove();
    URL.revokeObjectURL(url);
  }

  // ---------- DOM-attached measuring SVG ----------
  function getMeasureSVG() {
    let m = document.getElementById('DTSVG_MEASURE');
    if (m) return m;
    m = document.createElementNS(svgNS, 'svg');
    m.setAttribute('id', 'DTSVG_MEASURE');
    m.setAttribute('xmlns', svgNS);
    m.setAttribute('width', '10');
    m.setAttribute('height', '10');
    m.style.position = 'fixed';
    m.style.left = '-10000px';
    m.style.top  = '0';
    m.style.opacity = '0';
    m.style.pointerEvents = 'none';
    document.body.appendChild(m);
    return m;
  }
  const MEASURE_SVG = getMeasureSVG();

  function measureString(str, opts = {}) {
    const t = document.createElementNS(svgNS, 'text');
    t.textContent = str || '';
    if (opts.fontFamily) t.setAttribute('font-family', opts.fontFamily);
    if (opts.fontSize)   t.setAttribute('font-size',   String(opts.fontSize));
    if (opts.fontWeight) t.setAttribute('font-weight', String(opts.fontWeight));
    MEASURE_SVG.appendChild(t);
    const w = t.getComputedTextLength() || 0;
    MEASURE_SVG.removeChild(t);
    return w;
  }
  function copyTextStyle(from, to) {
    const attrs = ['font-family','font-size','font-weight','text-anchor','dominant-baseline','letter-spacing'];
    attrs.forEach(a => {
      const v = from.getAttribute(a);
      if (v != null) to.setAttribute(a, v);
    });
  }
  function clampTextToWidth(textEl, maxWidth) {
    if (!textEl || !isFinite(maxWidth) || maxWidth <= 0) return;
    const probe = document.createElementNS(svgNS, 'text');
    copyTextStyle(textEl, probe);
    const original = textEl.textContent || '';
    probe.textContent = original;
    MEASURE_SVG.appendChild(probe);
    if (probe.getComputedTextLength() <= maxWidth) {
      MEASURE_SVG.removeChild(probe);
      return;
    }
    let lo = 0, hi = original.length, best = '…';
    while (lo <= hi) {
      const mid = (lo + hi) >> 1;
      const candidate = original.slice(0, mid) + '…';
      probe.textContent = candidate;
      if (probe.getComputedTextLength() <= maxWidth) { best = candidate; lo = mid + 1; }
      else hi = mid - 1;
    }
    textEl.textContent = best;
    MEASURE_SVG.removeChild(probe);
  }

  // ---------- parsing / colors ----------
  function parseHeaderCell(th) {
    const small = th.querySelector('small');
    const line1 = th.cloneNode(true);
    Array.from(line1.querySelectorAll('small')).forEach(el => el.remove());
    const primary   = (line1.textContent || '').trim().replace(/\s+/g, ' ');
    const secondary = small ? (small.textContent || '').trim().replace(/\s+/g, ' ') : '';
    return { primary, secondary };
  }
  function toRGB(color) {
    if (!color) return 'rgb(255,255,255)';
    if (/^rgba?\(/i.test(color)) return color;
    const tmp = document.createElement('div');
    tmp.style.color = color;
    document.body.appendChild(tmp);
    const rgb = getComputedStyle(tmp).color;
    document.body.removeChild(tmp);
    return rgb || 'rgb(255,255,255)';
  }

  // --- helpers to embed <img> icons into SVG ---
  async function urlToDataURL(url) {
    const res = await fetch(url, { credentials: 'same-origin' });
    const blob = await res.blob();
    return await new Promise(r => {
      const fr = new FileReader();
      fr.onload = () => r(fr.result);
      fr.readAsDataURL(blob);
    });
  }
  function svgImage(href, x, y, w, h) {
    const im = document.createElementNS(svgNS, 'image');
    im.setAttributeNS('http://www.w3.org/1999/xlink', 'href', href);
    im.setAttribute('x', x); im.setAttribute('y', y);
    im.setAttribute('width', w); im.setAttribute('height', h);
    return im;
  }

  // ---------- SVG prims ----------
  function createSVG(w, h, background /* null => transparent */) {
    const svg = document.createElementNS(svgNS, 'svg');
    svg.setAttribute('xmlns', svgNS);
    svg.setAttribute('width',  String(w));
    svg.setAttribute('height', String(h));
    svg.setAttribute('viewBox', `0 0 ${w} ${h}`);
    if (background && String(background).toLowerCase() !== 'transparent') {
      const bg = document.createElementNS(svgNS, 'rect');
      bg.setAttribute('x', '0'); bg.setAttribute('y', '0');
      bg.setAttribute('width', String(w)); bg.setAttribute('height', String(h));
      bg.setAttribute('fill', background);
      svg.appendChild(bg);
    }
    return svg;
  }
  function rect(x,y,w,h, fill, stroke, rx) {
    const r = document.createElementNS(svgNS, 'rect');
    r.setAttribute('x', x); r.setAttribute('y', y);
    r.setAttribute('width', w); r.setAttribute('height', h);
    if (fill   != null) r.setAttribute('fill', fill);
    if (stroke != null) r.setAttribute('stroke', stroke);
    if (rx     != null) { r.setAttribute('rx', rx); r.setAttribute('ry', rx); }
    return r;
  }
  function text(x,y, content, opts={}) {
    const t = document.createElementNS(svgNS, 'text');
    t.setAttribute('x', x); t.setAttribute('y', y);
    t.setAttribute('font-family', opts.fontFamily || 'system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial');
    t.setAttribute('font-size', (opts.fontSize || 14));
    if (opts.fontWeight) t.setAttribute('font-weight', opts.fontWeight);
    if (opts.anchor)     t.setAttribute('text-anchor', opts.anchor);
    if (opts.fill)       t.setAttribute('fill', opts.fill);
    if (opts.baseline)   t.setAttribute('dominant-baseline', opts.baseline);
    t.textContent = content;
    return t;
  }

  // ---------- wrap "(secondary text)" to two lines (with parens) ----------
  function wrapSecondary(secondary, maxW, size, weight, fontFamily) {
    if (!secondary) return [];
    const raw = secondary.replace(/^\(|\)$/g,'').trim();
    if (!raw) return [];

    const measureOpts = { fontSize:size, fontWeight:weight, fontFamily };
    const SOFT_WRAP_RATIO = 0.90;

    const oneLineWidth = measureString('(' + raw + ')', measureOpts);
    const mustWrap = oneLineWidth > maxW * SOFT_WRAP_RATIO;
    if (!mustWrap) return ['(' + raw + ')'];

    const words = raw.split(/\s+/);
    let best = null;

    for (let i = 1; i < words.length; i++) {
      const left  = words.slice(0, i).join(' ');
      const right = words.slice(i).join(' ');

      const wL = measureString('(' + left,  measureOpts);
      const wR = measureString(right + ')', measureOpts);

      const fits = (wL <= maxW && wR <= maxW);
      const max  = Math.max(wL, wR);
      if (fits) {
        if (!best || max < best.max) best = { left, right, max, fits:true };
      } else if (!best || best.fits || max < best.max) {
        best = { left, right, max, fits:false };
      }
    }
    return ['(' + best.left, best.right + ')'];
  }

  // ---------- clipping ----------
  let _clipIdCounter = 0;
  function makeClip(svgRoot, x, y, w, h) {
    const defs = svgRoot.querySelector('defs') || svgRoot.insertBefore(document.createElementNS(svgNS, 'defs'), svgRoot.firstChild);
    const id = 'cellClip' + (++_clipIdCounter);
    const cp = document.createElementNS(svgNS, 'clipPath');
    cp.setAttribute('id', id);
    cp.appendChild(rect(x, y, w, h, '#000', null, null));
    defs.appendChild(cp);
    return `url(#${id})`;
  }

  // ---------- HEATMAP table builder (your grid) ----------
  function buildTableSVG(table) {
    const thead = table.querySelector('thead');
    const tbody = table.querySelector('tbody');
    if (!thead || !tbody) throw new Error('DataTablesSVG: expected <thead> and <tbody>.');

    const firstRow = tbody.rows[0] || thead.rows[0];
    const firstDataCell = tbody.querySelector('td') || table.querySelector('td');
    const firstRowHeader = tbody.querySelector('th') || table.querySelector('th');
    if (!firstRow || !firstDataCell || !firstRowHeader) {
      throw new Error('DataTablesSVG: missing cells to infer sizes.');
    }

    const cellW = Math.round(firstDataCell.getBoundingClientRect().width);
    const cellH = Math.round(firstDataCell.getBoundingClientRect().height);
    const rowHeaderW = Math.round(firstRowHeader.getBoundingClientRect().width);
    const borderPx = 1;

    const numFontSize  = parseFloat(getComputedStyle(firstDataCell).fontSize) || 14;
    const headFontSize = parseFloat(getComputedStyle(thead.querySelector('th')).fontSize) || 16;
    const smallEl = thead.querySelector('small') || table.querySelector('small');
    const head2FontSize = smallEl ? parseFloat(getComputedStyle(smallEl).fontSize) : Math.max(11, headFontSize - 4);
    const numFontWeight = getComputedStyle(firstDataCell).fontWeight || 550;

    const cols = thead.rows[0] ? thead.rows[0].cells.length : (firstRow.cells.length);
    const width  = rowHeaderW + (cols - 1) * cellW + borderPx;
    const height = (1 * cellH) + (tbody.rows.length) * cellH + borderPx;

    const svg = createSVG(width, height, '#ffffff');

    // Corner
    svg.appendChild(rect(0, 0, rowHeaderW, cellH, '#ffffff', '#dddddd'));
    svg.appendChild(text(rowHeaderW/2, cellH/2, 'Classes', {
      fontSize: headFontSize, fontWeight: 600, anchor: 'middle', baseline: 'middle', fill: '#000'
    }));

    const innerPad = 6;
    const lineGap  = Math.max(2, headFontSize * 0.2);
    const fontFamily = 'system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial';

    // Column headers
    for (let c = 1; c < cols; c++) {
      const x = rowHeaderW + (c-1)*cellW;
      const th = thead.rows[0].cells[c];
      const { primary, secondary } = parseHeaderCell(th);

      svg.appendChild(rect(x, 0, cellW, cellH, '#ffffff', '#dddddd'));
      const clip = makeClip(svg, x+1, 1, cellW-2, cellH-2);
      const g = document.createElementNS(svgNS, 'g'); g.setAttribute('clip-path', clip);

      const lines = [{ txt: primary, size: headFontSize, weight: 600, color: '#000' }];
      if (secondary) {
        const maxW = cellW - innerPad*2;
        wrapSecondary(secondary, maxW, head2FontSize, 500, fontFamily)
          .forEach(s => lines.push({ txt: s, size: head2FontSize, weight: 500, color: '#6c757d' }));
      }
      const totalH = lines.reduce((a,ln,i)=>a+ln.size+(i?lineGap:0),0);
      let y = (cellH - totalH)/2;

      lines.forEach(ln => {
        const ty = y + ln.size/2;
        const tNode = text(x + cellW/2, ty, ln.txt, {
          fontSize: ln.size, fontWeight: ln.weight, anchor: 'middle', baseline: 'middle', fill: ln.color, fontFamily
        });
        g.appendChild(tNode);
        clampTextToWidth(tNode, cellW - innerPad*2);
        y += ln.size + lineGap;
      });

      svg.appendChild(g);
    }

    // Body rows
    for (let r = 0; r < tbody.rows.length; r++) {
      const y = (r+1)*cellH;
      const tr = tbody.rows[r];

      // Row header
      const th = tr.cells[0];
      const { primary, secondary } = parseHeaderCell(th);
      svg.appendChild(rect(0, y, rowHeaderW, cellH, '#ffffff', '#dddddd'));
      const clip = makeClip(svg, 1, y+1, rowHeaderW-2, cellH-2);
      const g = document.createElementNS(svgNS, 'g'); g.setAttribute('clip-path', clip);

      const lines = [{ txt: primary, size: headFontSize, weight: 600, color: '#000' }];
      if (secondary) {
        const maxWRow = rowHeaderW - innerPad*2;
        wrapSecondary(secondary, maxWRow, head2FontSize, 500, fontFamily)
          .forEach(s => lines.push({ txt: s, size: head2FontSize, weight: 500, color: '#6c757d' }));
      }
      const totalH = lines.reduce((a,ln,i)=>a+ln.size+(i?lineGap:0),0);
      let yy = y + (cellH - totalH)/2;

      lines.forEach(ln => {
        const ty = yy + ln.size/2;
        const tNode = text(rowHeaderW/2, ty, ln.txt, {
          fontSize: ln.size, fontWeight: ln.weight, anchor: 'middle', baseline: 'middle', fill: ln.color, fontFamily
        });
        g.appendChild(tNode);
        clampTextToWidth(tNode, rowHeaderW - innerPad*2);
        yy += ln.size + lineGap;
      });

      svg.appendChild(g);

      // Data cells
      for (let c = 1; c < tr.cells.length; c++) {
        const x = rowHeaderW + (c-1)*cellW;
        const td = tr.cells[c];
        const bg = getComputedStyle(td).backgroundColor || 'rgb(255,255,255)';
        const fill = toRGB(bg);
        const txt = (td.textContent || '').trim();
        svg.appendChild(rect(x, y, cellW, cellH, fill, '#dddddd'));
        if (txt) {
          svg.appendChild(text(x + cellW/2, y + cellH/2, txt, {
            fontSize: numFontSize, fontWeight: numFontWeight, anchor: 'middle', baseline: 'middle', fill: '#000', fontFamily
          }));
        }
      }
    }

    const rowsCount = tbody.rows.length;
    return { svg, width: width, height: height, rowsCount };
  }

  // ---------- GENERIC table builder (flyout & friends) ----------
async function buildFlyoutSVG(flyBox, table) {
  const thead = table.querySelector('thead');
  const tbody = table.querySelector('tbody');
  if (!thead || !tbody) throw new Error('Flyout: missing thead/tbody');

  // measure
  const bodyRow0 = tbody.rows[0] || table.rows[0];
  if (!bodyRow0) throw new Error('Flyout: empty tbody');

  const colW  = Array.from(bodyRow0.cells).map(td => Math.round(td.getBoundingClientRect().width));
  const rowH  = Math.round(bodyRow0.getBoundingClientRect().height);
  const headRows    = Array.from(thead.rows);
  const headHeights = headRows.map(tr => Math.round(tr.getBoundingClientRect().height));
  const headH = headHeights.reduce((a,b)=>a+b,0);

  const totalW = colW.reduce((a,b)=>a+b,0);
  const totalH = headH + rowH * tbody.rows.length;

  // palette
  const bgPanel   = '#ffffff';
  const headerBg  = '#f6f6f6';
  const headerTxt = '#222';
  const gridColor = '#e6e6e6';
  const hardLine  = '#000000';
  const radius    = 10;

  // root / panel
  const svg = createSVG(totalW, totalH, bgPanel);
  svg.setAttribute('shape-rendering', 'crispEdges');
  svg.appendChild(rect(0, 0, totalW, totalH, bgPanel, '#cfcfcf', radius));

  // cumulative Xs
  const cumX = [0];
  for (let i = 0; i < colW.length; i++) cumX.push(cumX[i] + colW[i]); // len = cols+1

  // ----- HEADER -----
  let y = 0;
  let superBottomY = null; // bottom Y of the .ccs-super row
  headRows.forEach((tr, rIdx) => {
    const h = headHeights[rIdx];
    let colIdx = 0;

    Array.from(tr.cells).forEach(th => {
      const span = th.colSpan || 1;
      const x1 = cumX[colIdx];
      const x2 = cumX[colIdx + span];
      const w  = x2 - x1;

      // header cell bg (no inner strokes -> no grey verticals in header)
      svg.appendChild(rect(x1, y, w, h, headerBg, null));

      // centered text
      const label = (th.textContent || '').replace(/\s+/g,' ').trim();
      svg.appendChild(text(x1 + w/2, y + h/2, label, {
        fontSize: 14, fontWeight: 600, anchor: 'middle', baseline: 'middle', fill: headerTxt
      }));

      colIdx += span;
    });

    y += h;
    if (tr.classList.contains('ccs-super')) superBottomY = y;
  });

  // black horizontal segments ONLY under the 2 group superheaders (“B1…”, “B2…”)
  if (superBottomY != null) {
    const yLine = Math.round(superBottomY) + 0.5;
    const seg = (x1, x2) => {
      const p = document.createElementNS(svgNS,'path');
      p.setAttribute('d', `M ${x1} ${yLine} H ${x2}`);
      p.setAttribute('stroke', hardLine);
      p.setAttribute('stroke-width', 1);
      p.setAttribute('stroke-linecap', 'square');
      svg.appendChild(p);
    };
    // groups start after the 2 metric cols: [2..6) and [6..end)
    seg(cumX[2] + 0.5, cumX[6] - 0.5);                     // left group (cols 3..6)
    seg(cumX[6] + 0.5, cumX[cumX.length - 1] - 0.5);       // right group (cols 7..end)
  }

  // ----- BODY (cells keep their computed bg incl. ID/SIM colors) -----
  for (let r = 0; r < tbody.rows.length; r++) {
    const tr = tbody.rows[r];
    const top = headH + r * rowH;

    for (let c = 0; c < tr.cells.length; c++) {
      const td = tr.cells[c];
      const x1 = cumX[c], w = colW[c];

      // background (zebra/white or metric color)
      const fill = toRGB(getComputedStyle(td).backgroundColor || '#ffffff');
      svg.appendChild(rect(x1, top, w, rowH, fill, null));

      // icon vs text
      const img = td.querySelector('img');
      if (img && img.src) {
        try {
          const href = await urlToDataURL(img.src);
          const sz = 12;
          svg.appendChild(svgImage(href, x1 + (w - sz)/2, top + (rowH - sz)/2, sz, sz));
        } catch(_) {}
      } else {
        const txt = (td.textContent || '').trim();
        if (txt) {
          const bold = (c < 2) ? 600 : 500; // ID/SIM slightly bolder
          svg.appendChild(text(x1 + w/2, top + rowH/2, txt, {
            fontSize: 14, fontWeight: bold, anchor: 'middle', baseline: 'middle', fill: '#000'
          }));
        }
      }
    }
  }

  // ----- GRID LINES -----

  // body horizontals (light)
  for (let r = 1; r < tbody.rows.length; r++) {
    const yLine = headH + r*rowH + 0.5;
    const ph = document.createElementNS(svgNS,'path');
    ph.setAttribute('d', `M 0 ${yLine} H ${totalW}`);
    ph.setAttribute('stroke', gridColor);
    ph.setAttribute('stroke-width', 1);
    svg.appendChild(ph);
  }

  // light verticals ONLY in body (no grey lines in header)
  for (let c = 1; c < colW.length; c++) {
    const xLine = cumX[c] + 0.5;
    const pv = document.createElementNS(svgNS,'path');
    pv.setAttribute('d', `M ${xLine} ${headH + 0.5} V ${totalH - 0.5}`);
    pv.setAttribute('stroke', gridColor);
    pv.setAttribute('stroke-width', 1);
    svg.appendChild(pv);
  }

  // thick vertical dividers after columns 2 and 6 (across whole panel)
  [2, 6].forEach(ix => {
    if (ix < colW.length) {
      const xDiv = cumX[ix] + 0.5;
      const pv = document.createElementNS(svgNS,'path');
      pv.setAttribute('d', `M ${xDiv} 0 V ${totalH}`);
      pv.setAttribute('stroke', hardLine);
      pv.setAttribute('stroke-width', 1);
      pv.setAttribute('stroke-linecap', 'square');
      svg.appendChild(pv);
    }
  });

  // === HEADER/BODY BLACK DIVIDER ON TOP (after grid lines) ===
  {
    const yDivider = Math.round(headH) + 0.5;
    const p = document.createElementNS(svgNS,'path');
    p.setAttribute('d', `M 0 ${yDivider} H ${totalW}`);
    p.setAttribute('stroke', hardLine);
    p.setAttribute('stroke-width', 1.5);   // your chosen thickness
    p.setAttribute('stroke-linecap', 'square');
    svg.appendChild(p);
  }

  // outer rounded border (optional – keep if you like)
  const outline = document.createElementNS(svgNS,'rect');
  outline.setAttribute('x','0.5'); outline.setAttribute('y','0.5');
  outline.setAttribute('width', String(totalW-1)); outline.setAttribute('height', String(totalH-1));
  outline.setAttribute('fill','none'); outline.setAttribute('stroke', hardLine);
  outline.setAttribute('rx', String(radius)); outline.setAttribute('ry', String(radius));
  svg.appendChild(outline);

  return { svg, width: totalW, height: totalH };
}


  // ---------- public API ----------
  function exportTable(selector, baseName, opts={}) {
    const table = (typeof selector === 'string') ? document.querySelector(selector) : selector;
    if (!table) return console.error('DataTablesSVG: table not found:', selector);
    const { svg } = buildTableSVG(table);
    const xml = new XMLSerializer().serializeToString(svg);
    const filename = resolveFilename(
      typeof selector === 'string' ? selector : (baseName || 'table'),
      'svg',
      opts.filename,
      baseName || 'table'
    );
    triggerDownload(xml, filename);
  }

  function exportWithBadges(axisWrapSelector, tableSelector, baseName, userOpts) {
    const opts = Object.assign({
      padTop: 60,
      padLeft: 60,
      gapX: 8,
      gapY: 10,
      badgeTighten: 40,
      background: null,
      filename: null
    }, userOpts || {});

    const axisWrap = (typeof axisWrapSelector === 'string') ? document.querySelector(axisWrapSelector) : axisWrapSelector;
    const table    = (typeof tableSelector === 'string')    ? document.querySelector(tableSelector)    : tableSelector;
    if (!axisWrap || !table) return console.error('DataTablesSVG.withBadges: missing axisWrap or table.');

    const { svg: tableSVG, width: tW, height: tH } = buildTableSVG(table);

    function readBadge(el) {
      if (!el) return null;
      const cs = getComputedStyle(el);
      const r  = el.getBoundingClientRect();
      return {
        text: (el.textContent || '').trim(),
        w: Math.round(r.width),
        h: Math.round(r.height),
        fill: toRGB(cs.backgroundColor || '#ffffff'),
        stroke: toRGB(cs.borderColor || '#bbbbbb'),
        color: toRGB(cs.color || '#000000'),
        radius: Math.max(0, parseFloat(cs.borderTopLeftRadius) || 12),
        fontSize: Math.max(10, parseFloat(cs.fontSize) || 16),
        fontWeight: cs.fontWeight || 600
      };
    }
    const xb = readBadge(axisWrap.querySelector('.axis-label-x .ccs-axis-badge'));
    const yb = readBadge(axisWrap.querySelector('.axis-label-y .ccs-axis-badge'));

    const yRotatedWidth = yb ? yb.h : 0;
    const padTop  = Math.max(opts.padTop,  xb ? xb.h + opts.gapX : opts.padTop);
    const padLeft = Math.max(opts.padLeft, yRotatedWidth + opts.gapY);
    const W = padLeft + tW;
    const H = padTop  + tH;

    const out = createSVG(W, H, opts.background);

    // Table group
    const gTable = document.createElementNS(svgNS, 'g');
    gTable.setAttribute('transform', `translate(${padLeft}, ${padTop})`);
    while (tableSVG.firstChild) gTable.appendChild(tableSVG.firstChild);
    out.appendChild(gTable);

    // X badge
    if (xb) {
      const bx = padLeft + Math.round((tW - xb.w) / 2);
      const by = Math.max(0, padTop - opts.gapX - xb.h);
      out.appendChild(rect(bx, by, xb.w, xb.h, xb.fill, xb.stroke, xb.radius));
      out.appendChild(text(bx + xb.w/2, by + xb.h/2, xb.text, {
        fontSize: xb.fontSize, fontWeight: xb.fontWeight, anchor: 'middle', baseline: 'middle', fill: xb.color
      }));
    }

    // Y badge (rotated)
    if (yb) {
      const Wb = yb.w, Hb = yb.h;
      const targetRight = padLeft - Math.max(0, opts.gapY) + Math.max(0, opts.badgeTighten);
      const cx = targetRight - (Hb / 2);
      const cy = padTop + (tH / 2);
      const gY = document.createElementNS(svgNS, 'g');
      gY.setAttribute('transform', `rotate(-90, ${cx}, ${cy})`);
      gY.appendChild(rect(cx - Hb/2, cy - Wb/2, Hb, Wb, yb.fill, yb.stroke, yb.radius));
      gY.appendChild(text(cx, cy, yb.text, {
        fontSize: yb.fontSize, fontWeight: yb.fontWeight, anchor: 'middle', baseline: 'middle', fill: yb.color
      }));
      out.appendChild(gY);
    }

    const xml = new XMLSerializer().serializeToString(out);
    const filename = resolveFilename(
      (typeof tableSelector === 'string') ? tableSelector : (baseName || 'ccs-heatmap'),
      'svg',
      opts.filename,
      baseName || 'ccs-heatmap'
    );
    triggerDownload(xml, filename);
  }

  // NEW: export any generic table (used for flyout)
  async function exportFlyout(flySelector, tableSelector, baseName, opts={}) {
    const fly  = typeof flySelector === 'string' ? document.querySelector(flySelector) : flySelector;
    const tbl  = typeof tableSelector === 'string' ? document.querySelector(tableSelector) : tableSelector;
    if (!fly || !tbl) return console.error('DataTablesSVG.exportFlyout: selectors not found');

    const { svg } = await buildFlyoutSVG(fly, tbl);
    const xml = new XMLSerializer().serializeToString(svg);
    const filename = resolveFilename(baseName || 'ccs-flyout', 'svg', opts.filename, 'ccs-flyout');
    triggerDownload(xml, filename);
  }


  global.DataTablesSVG = { export: exportTable, exportWithBadges, exportFlyout };
})(window);

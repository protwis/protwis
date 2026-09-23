/**
 * DataTablesImage.js (enhanced with link/selector stripping)
 * ----------------------------------------------------------
 * Tight-crop raster export of a single DOM element via html2canvas.
 *
 * Usage:
 *   // keep links (default)
 *   DataTablesImage.export('#ccsTable', 'ccs-table', { format:'png', scale:3 });
 *
 *   // strip links + clickable spans to plain black, no underline
 *   DataTablesImage.export('#gpcrTable', 'gpcr-table', {
 *     format:'jpg', scale:3, links:'strip'
 *   });
 *
 * Options:
 *   format: 'png' | 'jpg' (default 'png')
 *   quality: 0.92 (JPG)
 *   scale: number (default max(2, devicePixelRatio*2))
 *   background: '#ffffff'
 *   filename: string
 *   // link/selector handling:
 *   links: 'keep' | 'strip'         (default 'keep')
 *   linkColor: '#000'               (used when stripping)
 *   linkUnderline: false            (used when stripping)
 *   stripSelectors: string          (CSS selector list to neutralize; default 'a, .popup-cell')
 *   sandboxClone: true              (render from hidden clone; default true)
 *
 * Depends on: html2canvas
 */
(function (global) {
  if (typeof html2canvas === 'undefined') {
    console.error('DataTablesImage.js: html2canvas is required.');
    return;
  }

  /* ---------------- util ---------------- */
  function triggerDownload(dataUrl, filename) {
    const a = document.createElement('a');
    a.href = dataUrl;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
  }

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

  function resolveFilename(selector, baseName, ext, explicit) {
    if (explicit) return explicit; // honor explicit filename
    const base = (baseName && String(baseName).trim()) || basenameFromSelector(selector, 'table');
    return stamped(base, ext);
  }

  function loosenOverflow(el) {
    const parent = el && el.parentElement;
    if (!parent) return () => {};
    const prev = {
      overflow: parent.style.overflow,
      overflowX: parent.style.overflowX,
      overflowY: parent.style.overflowY,
      maxWidth: parent.style.maxWidth,
      maxHeight: parent.style.maxHeight
    };
    parent.style.overflow = 'visible';
    parent.style.overflowX = 'visible';
    parent.style.overflowY = 'visible';
    parent.style.maxWidth = 'none';
    parent.style.maxHeight = 'none';
    return () => {
      parent.style.overflow   = prev.overflow;
      parent.style.overflowX  = prev.overflowX;
      parent.style.overflowY  = prev.overflowY;
      parent.style.maxWidth   = prev.maxWidth;
      parent.style.maxHeight  = prev.maxHeight;
    };
  }

  /* --------- clone sandbox for html2canvas --------- */
  function buildSandboxClone(el, options) {
    // Create an offscreen wrapper so layout matches page CSS
    const wrap = document.createElement('div');
    wrap.style.position = 'fixed';
    wrap.style.left = '-10000px';
    wrap.style.top = '0';
    wrap.style.pointerEvents = 'none';
    wrap.style.opacity = '1';

    // Match width of original to avoid reflow differences
    const w = el.getBoundingClientRect().width;
    if (w) wrap.style.width = Math.ceil(w) + 'px';

    // Deep clone the target subtree and attach it first
    const clone = el.cloneNode(true);
    wrap.appendChild(clone);
    document.body.appendChild(wrap);

    // Optionally neutralize links / clickable UI
    if ((options.links || 'keep').toLowerCase() === 'strip') {
      const color = options.linkColor || '#000';
      const underline = !!options.linkUnderline;
      const selectorList = options.stripSelectors || 'a, .popup-cell';

      // Replace <a> with <span>, and restyle other matches in-place
      const nodes = clone.querySelectorAll(selectorList);
      nodes.forEach(node => {
        if (node.tagName && node.tagName.toLowerCase() === 'a') {
          const span = document.createElement('span');
          span.textContent = node.textContent || '';
          span.style.color = color;
          span.style.textDecoration = underline ? 'underline' : 'none';
          span.style.cursor = 'default';
          // Preserve basic typography
          span.style.font = 'inherit';
          span.style.fontWeight = 'inherit';
          span.style.fontSize = 'inherit';
          span.style.lineHeight = 'inherit';
          span.style.whiteSpace = 'inherit';
          node.parentNode.replaceChild(span, node);
        } else {
          // Any non-<a> match (e.g. .popup-cell spans)
          node.style.color = color;
          node.style.textDecoration = underline ? 'underline' : 'none';
          node.style.cursor = 'default';
        }
      });
    }

    return { wrap, clone, cleanup: () => wrap.remove() };
  }

  /* ---------------- render/export ---------------- */
  async function renderElement(selector, opts) {
    const options = Object.assign(
      {
        format: 'png',
        scale: Math.max(2, (window.devicePixelRatio || 1) * 2),
        background: '#ffffff',
        links: 'keep',                 // 'keep' | 'strip'
        linkColor: '#000',
        linkUnderline: false,
        stripSelectors: 'a, .popup-cell',
        sandboxClone: true             // render a hidden clone when true
      },
      opts || {}
    );

    const el = (typeof selector === 'string') ? document.querySelector(selector) : selector;
    if (!el) throw new Error('DataTablesImage.render: selector not found: ' + selector);

    let target = el;
    let cleanupSandbox = null;

    // Render from a hidden clone if requested or when stripping
    if (options.sandboxClone || (options.links && options.links.toLowerCase() === 'strip')) {
      const { clone, cleanup } = buildSandboxClone(el, options);
      target = clone;
      cleanupSandbox = cleanup;

      const restore = loosenOverflow(target);
      await new Promise(r => requestAnimationFrame(r));
      const canvas = await html2canvas(target, {
        backgroundColor: options.background,
        useCORS: true,
        scale: options.scale,
        scrollX: 0,
        scrollY: 0
      });
      restore();
      cleanupSandbox && cleanupSandbox();
      return canvas;
    }

    // Live element path (no clone)
    const restore = loosenOverflow(target);
    await new Promise(r => requestAnimationFrame(r));
    const canvas = await html2canvas(target, {
      backgroundColor: options.background,
      useCORS: true,
      scale: options.scale,
      scrollX: 0,
      scrollY: 0
    });
    restore();
    return canvas;
  }

  async function exportElement(selector, baseName, opts) {
    const options = Object.assign(
      {
        format: 'png',
        quality: 0.92,
        scale: Math.max(2, (window.devicePixelRatio || 1) * 2),
        background: '#ffffff',
        filename: null,
        links: 'keep',
        linkColor: '#000',
        linkUnderline: false,
        stripSelectors: 'a, .popup-cell',
        sandboxClone: true
      },
      opts || {}
    );

    const fmt = String(options.format || 'png').toLowerCase();
    const canvas = await renderElement(selector, options);
    const filename = resolveFilename(selector, baseName, fmt, options.filename);

    if (fmt === 'png') {
      triggerDownload(canvas.toDataURL('image/png'), filename);
    } else if (fmt === 'jpg' || fmt === 'jpeg') {
      triggerDownload(canvas.toDataURL('image/jpeg', options.quality), filename);
    } else {
      console.error('DataTablesImage.export: unsupported format (use png/jpg):', fmt);
    }
  }

  global.DataTablesImage = {
    render: renderElement,
    export: exportElement
  };
})(window);

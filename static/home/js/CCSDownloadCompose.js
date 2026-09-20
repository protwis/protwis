/**
 * CCSDownloadCompose.js
 * ---------------------
 * Page-specific: compose a bitmap of the table (from DataTablesImage.render)
 * with top/left padding and the axis badges, then download as PNG/JPG/PDF.
 *
 * Depends on: html2canvas, jsPDF (UMD), DataTablesImage.render
 */
(function (global) {
  if (!global.DataTablesImage || !global.DataTablesImage.render) {
    console.error('CCSDownloadCompose: DataTablesImage.render is required.');
    return;
  }
  if (typeof html2canvas === 'undefined') {
    console.error('CCSDownloadCompose: html2canvas is required.');
    return;
  }

  function basenameFromSelector(selector, fallback) {
    if (typeof selector === 'string') {
      const m = selector.match(/#([\w-]+)/);
      if (m && m[1]) return m[1];
    }
    return fallback || 'ccs-heatmap';
  }
  function stamped(base, ext) {
    const d = new Date(), pad = n => String(n).padStart(2, '0');
    const stamp = `${d.getFullYear()}-${pad(d.getMonth()+1)}-${pad(d.getDate())}_` +
                  `${pad(d.getHours())}-${pad(d.getMinutes())}-${pad(d.getSeconds())}`;
    return `${base}-${stamp}.${ext}`;
  }
  function resolveFilename(tableSelector, baseName, ext, explicit) {
    if (explicit) return explicit;
    const base = (baseName && String(baseName).trim()) || basenameFromSelector(tableSelector, 'ccs-heatmap');
    return stamped(base, ext);
  }

  function triggerDownload(dataUrl, filename) {
    const a = document.createElement('a');
    a.href = dataUrl; a.download = filename; document.body.appendChild(a); a.click(); a.remove();
  }

  async function captureBadgeCanvas(el, scale, { fixVertical=false } = {}) {
    if (!el) return null;

    const clone = el.cloneNode(true);
    clone.style.position = 'fixed';
    clone.style.left = '-10000px';
    clone.style.top = '0';
    clone.style.margin = '0';

    if (fixVertical) {
      clone.style.writingMode = 'horizontal-tb';
      clone.style.transform = 'rotate(-90deg)';
      clone.style.transformOrigin = 'left top';
      clone.style.display = 'inline-block';
    }

    document.body.appendChild(clone);
    await new Promise(r => requestAnimationFrame(r));

    const canvas = await html2canvas(clone, {
      backgroundColor: null,
      useCORS: true,
      scale: scale || Math.max(2, (window.devicePixelRatio || 1) * 2),
      scrollX: 0,
      scrollY: 0
    });

    document.body.removeChild(clone);
    return canvas;
  }

  async function composeAndDownload(format, userOpts) {
    const opts = Object.assign({
      tableSelector: '#ccsTable',
      axisWrapSelector: '.ccs-axis-wrap',
      padTop: 60,
      padLeft: 60,
      gapX: 8,
      gapY: 10,
      scale: Math.max(2, (window.devicePixelRatio || 1) * 2),
      quality: 0.95,
      baseName: null,        // if null -> derive from table ID
      filename: null,        // NEW: use this if provided (no timestamp)
      // backgrounds:
      tableBackground: '#ffffff',
      background: null,
      pdfBackground: null
    }, userOpts || {});

    const ext = String(format || 'png').toLowerCase();
    const wantsTransparentOuter =
      opts.background === null || String(opts.background).toLowerCase() === 'transparent';

    // 1) Render the table
    const tableCanvas = await window.DataTablesImage.render(
      opts.tableSelector,
      { scale: opts.scale, background: opts.tableBackground }
    );

    // 2) Badges
    const axisWrap = document.querySelector(opts.axisWrapSelector);
    const xBadgeEl = axisWrap && axisWrap.querySelector('.axis-label-x .ccs-axis-badge');
    const yBadgeEl = axisWrap && axisWrap.querySelector('.axis-label-y .ccs-axis-badge');

    const [xBadgeCanvas, yBadgeCanvas] = await Promise.all([
      captureBadgeCanvas(xBadgeEl, opts.scale),
      captureBadgeCanvas(yBadgeEl, opts.scale, { fixVertical: true })
    ]);

    // 3) Padding
    const padTop  = Math.max(opts.padTop,  (xBadgeCanvas ? xBadgeCanvas.height + opts.gapX : 0));
    const padLeft = Math.max(opts.padLeft, (yBadgeCanvas ? yBadgeCanvas.width  + opts.gapY : 0));

    // 4) Compose
    const W = padLeft + tableCanvas.width;
    const H = padTop  + tableCanvas.height;

    const out = document.createElement('canvas');
    out.width = W; out.height = H;
    const ctx = out.getContext('2d');

    if (!wantsTransparentOuter) {
      ctx.fillStyle = opts.background || '#ffffff';
      ctx.fillRect(0, 0, W, H);
    } else {
      ctx.clearRect(0, 0, W, H);
    }

    if (xBadgeCanvas) {
      const bx = padLeft + Math.round((tableCanvas.width - xBadgeCanvas.width) / 2);
      const by = Math.max(0, padTop - opts.gapX - xBadgeCanvas.height);
      ctx.drawImage(xBadgeCanvas, bx, by);
    }

    if (yBadgeCanvas) {
      const bx = Math.max(0, padLeft - opts.gapY - yBadgeCanvas.width);
      const by = padTop + Math.round((tableCanvas.height - yBadgeCanvas.height) / 2);
      ctx.drawImage(yBadgeCanvas, bx, by);
    }

    ctx.drawImage(tableCanvas, padLeft, padTop);

    // 5) Filename + save
    const filename = resolveFilename(opts.tableSelector, opts.baseName, ext, opts.filename);

    if (ext === 'png') {
      triggerDownload(out.toDataURL('image/png'), filename);
    } else if (ext === 'jpg' || ext === 'jpeg') {
      if (wantsTransparentOuter) {
        const flat = document.createElement('canvas');
        flat.width = W; flat.height = H;
        const fctx = flat.getContext('2d');
        fctx.fillStyle = '#ffffff';
        fctx.fillRect(0, 0, W, H);
        fctx.drawImage(out, 0, 0);
        triggerDownload(flat.toDataURL('image/jpeg', opts.quality), filename);
      } else {
        triggerDownload(out.toDataURL('image/jpeg', opts.quality), filename);
      }
    } else if (ext === 'pdf') {
      if (typeof window.jspdf === 'undefined') {
        console.error('CCSDownloadCompose: jsPDF is required for PDF output.');
        return;
      }
      const { jsPDF } = window.jspdf;
      const pdf = new jsPDF({ orientation: (W >= H) ? 'l' : 'p', unit: 'pt', format: [W, H] });
      if (opts.pdfBackground) {
        pdf.setFillColor(opts.pdfBackground);
        pdf.rect(0, 0, W, H, 'F');
      }
      pdf.addImage(out.toDataURL('image/png'), 'PNG', 0, 0, W, H);
      pdf.save(filename);
    } else {
      console.error('CCSDownloadCompose: unsupported format:', ext);
    }
  }

  global.CCSDownload = {
    withBadges: composeAndDownload
  };
})(window);

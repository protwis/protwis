// "Screen Snaps" -- site-wide dev tool. Floating button (bottom-right, every page) that lets
// you grab any DOM element as a downloadable file: a real, vector .svg if the element is or
// contains an <svg> (charts drawn with D3/etc.), otherwise a high-resolution .png rasterized
// via html2canvas (tables, form controls, arbitrary HTML -- there's no meaningful way to turn
// that into real vector SVG, see the "why" in the popover's own copy).
//
// Fully self-contained: builds its own button/popover DOM and injects its own <style> block at
// runtime. Disable site-wide by commenting out this file's <script> tag (and the html2canvas
// one above it) in home/templates/home/base.html.
(function () {
  "use strict";

  var QUICK_PICK_CLASS = "snap";
  var FLASH_MS = 1500;

  // ---------------------------------------------------------------------
  // SVG export -- ported from static/home/js/saveSvg.js's proven approach
  // (clone the live node so inline/D3-applied styles come along for free,
  // then inline any CSS-*class*-based rules too so the file is self-contained).
  // Renamed to avoid colliding with the separate global saveSvg()/collectSvgCss()
  // some pages already load -- this tool intentionally carries its own copy.
  // ---------------------------------------------------------------------
  function snapCollectSvgCss(svgEl) {
    var css = "";
    var sheets = document.styleSheets;
    for (var i = 0; i < sheets.length; i++) {
      var rules;
      try { rules = sheets[i].cssRules; } catch (e) { continue; }
      if (!rules) { continue; }
      for (var j = 0; j < rules.length; j++) {
        var rule = rules[j];
        try {
          if (rule instanceof CSSStyleRule && svgEl.querySelector(rule.selectorText)) {
            css += rule.selectorText + " { " + rule.style.cssText + " }\n";
          } else if (rule.cssText && rule.cssText.indexOf("@font-face") === 0) {
            css += rule.cssText + "\n";
          }
        } catch (e) { /* unparsable selector -- skip */ }
      }
    }
    return css;
  }

  function snapSaveSvg(svgEl, filename) {
    var clone = svgEl.cloneNode(true);

    var svgNS = "http://www.w3.org/2000/svg";
    var css = snapCollectSvgCss(svgEl);
    if (css) {
      var defs = clone.querySelector("defs");
      if (!defs) {
        defs = document.createElementNS(svgNS, "defs");
        clone.insertBefore(defs, clone.firstChild);
      }
      var styleEl = document.createElementNS(svgNS, "style");
      styleEl.setAttribute("type", "text/css");
      styleEl.textContent = css;
      defs.insertBefore(styleEl, defs.firstChild);
    }

    clone.setAttribute("xmlns", svgNS);
    clone.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink");

    var svgData = clone.outerHTML;
    var preface = '<?xml version="1.0" standalone="no"?>\r\n';
    var blob = new Blob([preface, svgData], { type: "image/svg+xml;charset=utf-8" });
    triggerDownload(URL.createObjectURL(blob), filename, true);
  }

  // ---------------------------------------------------------------------
  // PNG resolution metadata. Canvas-encoded PNGs carry no pHYs chunk, so
  // Illustrator/Word/journal checks read them as 72 ppi. Insert one right
  // after IHDR (8-byte signature + 25-byte IHDR = offset 33).
  // ---------------------------------------------------------------------
  var CRC_TABLE = null;

  function crc32(bytes) {
    if (!CRC_TABLE) {
      CRC_TABLE = [];
      for (var n = 0; n < 256; n++) {
        var c = n;
        for (var k = 0; k < 8; k++) { c = (c & 1) ? (0xEDB88320 ^ (c >>> 1)) : (c >>> 1); }
        CRC_TABLE[n] = c >>> 0;
      }
    }
    var crc = 0xFFFFFFFF;
    for (var i = 0; i < bytes.length; i++) { crc = CRC_TABLE[(crc ^ bytes[i]) & 0xFF] ^ (crc >>> 8); }
    return (crc ^ 0xFFFFFFFF) >>> 0;
  }

  function snapPngWithPpi(blob, ppi, cb) {
    var reader = new FileReader();
    reader.onerror = function () { cb(reader.error || new Error("could not read PNG blob")); };
    reader.onload = function () {
      var png = new Uint8Array(reader.result);
      var ppm = Math.round(ppi / 0.0254);  // pixels per metre

      var chunk = new Uint8Array(21);
      var view = new DataView(chunk.buffer);
      view.setUint32(0, 9);                                 // data length
      chunk.set([0x70, 0x48, 0x59, 0x73], 4);               // "pHYs"
      view.setUint32(8, ppm);                               // X pixels per unit
      view.setUint32(12, ppm);                              // Y pixels per unit
      chunk[16] = 1;                                        // unit: metre
      view.setUint32(17, crc32(chunk.subarray(4, 17)));     // CRC over type + data

      cb(null, new Blob([png.subarray(0, 33), chunk, png.subarray(33)], { type: "image/png" }));
    };
    reader.readAsArrayBuffer(blob);
  }

  function snapSavePng(targetEl, filename, ppi, onDone) {
    if (typeof window.html2canvas !== "function") {
      onDone(new Error("html2canvas is not loaded"));
      return;
    }

    // html2canvas rasterizes each inline <svg> by serializing it on its own, so CSS-*class*-based
    // svg styling (D3 charts) would be lost in the PNG. Tag the originals, then in the cloned
    // document inject the same collected CSS the vector export uses into each matching clone.
    var svgs = Array.prototype.slice.call(targetEl.querySelectorAll("svg"));
    if (isSvg(targetEl)) svgs.unshift(targetEl);
    svgs.forEach(function (svg, i) { svg.setAttribute("data-snap-svg-idx", i); });
    targetEl.setAttribute("data-snap-target", "1");
    function untag() {
      svgs.forEach(function (svg) { svg.removeAttribute("data-snap-svg-idx"); });
      targetEl.removeAttribute("data-snap-target");
    }

    // html2canvas doesn't paint the live page -- it re-lays-out a clone in its own iframe, which
    // has no scrollbar and so is a little wider. Anything sized in % / flex then lands somewhere
    // else than on screen, while JS-drawn svg overlays keep the pixel coordinates measured on the
    // live page (e.g. #all-gpcr's connector curves vs. its flex buttons) -- they drift apart, worst
    // when zoomed out, since the scrollbar is then more CSS px. Pinning the target and its svgs to
    // their live pixel size in the clone makes it lay out exactly as on screen.
    var frozen = [targetEl].concat(svgs).map(function (node) {
      var cs = window.getComputedStyle(node);
      return { node: node, width: cs.width, height: cs.height };
    });

    window.html2canvas(targetEl, {
      // One CSS px = 1/96 in, so this scale makes the file print at its on-screen size.
      scale: ppi / 96, backgroundColor: "#ffffff", useCORS: true,
      onclone: function (clonedDoc) {
        frozen.forEach(function (f) {
          var clone = (f.node === targetEl)
            ? clonedDoc.querySelector('[data-snap-target="1"]')
            : clonedDoc.querySelector('[data-snap-svg-idx="' + f.node.getAttribute("data-snap-svg-idx") + '"]');
          if (!clone || f.width === "auto" || f.height === "auto") return;
          clone.style.width = clone.style.minWidth = clone.style.maxWidth = f.width;
          clone.style.height = clone.style.minHeight = clone.style.maxHeight = f.height;
          clone.style.flex = "none";
        });
        svgs.forEach(function (svg, i) {
          var clone = clonedDoc.querySelector('[data-snap-svg-idx="' + i + '"]');
          var css = snapCollectSvgCss(svg);
          if (!clone || !css) return;
          var styleEl = clonedDoc.createElementNS("http://www.w3.org/2000/svg", "style");
          styleEl.textContent = css;
          clone.insertBefore(styleEl, clone.firstChild);
        });
      },
    })
      .then(function (canvas) {
        untag();
        canvas.toBlob(function (blob) {
          if (!blob) { onDone(new Error("canvas.toBlob returned nothing (image too large?)")); return; }
          snapPngWithPpi(blob, ppi, function (err, pngBlob) {
            if (err) { onDone(err); return; }
            triggerDownload(URL.createObjectURL(pngBlob), filename, true);
            onDone(null);
          });
        }, "image/png");
      })
      .catch(function (err) { untag(); onDone(err); });
  }

  function triggerDownload(url, filename, revoke) {
    var a = document.createElement("a");
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    if (revoke) { setTimeout(function () { URL.revokeObjectURL(url); }, 1000); }
  }

  // ---------------------------------------------------------------------
  // Visual "here's what I grabbed (or didn't)" outline flash on the actual
  // target element, plus the shared resolve -> export pipeline.
  // ---------------------------------------------------------------------
  function flashOutline(el, ok) {
    if (!el || !el.style) return;
    var prevOutline = el.style.outline;
    var prevOffset = el.style.outlineOffset;
    el.style.outline = (ok ? "3px solid #2ecc71" : "3px solid #e74c3c");
    el.style.outlineOffset = "2px";
    setTimeout(function () {
      el.style.outline = prevOutline;
      el.style.outlineOffset = prevOffset;
    }, FLASH_MS);
  }

  // ---------------------------------------------------------------------
  // Vector vs. raster: the snap should look like what's on screen. Only export
  // a real .svg when that one <svg> *is* everything visible in the element --
  // an element that merely contains an svg next to other HTML (e.g. the
  // segment-selection panel's connector-line overlay behind its buttons)
  // would otherwise lose all of that HTML, so it gets a PNG of the whole thing.
  // ---------------------------------------------------------------------
  var CONTENT_TAGS = { img: 1, canvas: 1, video: 1, iframe: 1, input: 1, select: 1, textarea: 1, button: 1 };

  function isRendered(el) {
    return !!el && el.getClientRects().length > 0 && window.getComputedStyle(el).visibility !== "hidden";
  }

  function isSvg(el) {
    return !!el && !!el.tagName && el.tagName.toLowerCase() === "svg";
  }

  // Returns the <svg> to export as vector, or null if the element should be rasterized.
  function svgOnlyContent(el) {
    if (isSvg(el)) return el;

    // Visible, top-level svgs only (an svg nested inside another svg is part of its parent).
    var svgs = Array.prototype.filter.call(el.querySelectorAll("svg"), function (svg) {
      return isRendered(svg) && !(svg.parentElement && svg.parentElement.closest("svg"));
    });
    if (svgs.length !== 1) return null;
    var svg = svgs[0];

    // Any visible text or content element outside that svg means it isn't the whole picture.
    var walker = document.createTreeWalker(el, NodeFilter.SHOW_ELEMENT | NodeFilter.SHOW_TEXT, {
      acceptNode: function (node) {
        return isSvg(node) ? NodeFilter.FILTER_REJECT : NodeFilter.FILTER_ACCEPT;
      },
    });
    var node;
    while ((node = walker.nextNode())) {
      if (node.nodeType === Node.TEXT_NODE) {
        if (node.nodeValue.trim() && isRendered(node.parentElement)) return null;
      } else if (CONTENT_TAGS[node.tagName.toLowerCase()] && isRendered(node)) {
        return null;
      }
    }
    return svg;
  }

  function snapElement(el, baseName, statusEl) {
    var svg = svgOnlyContent(el);
    var hasSvg = !svg && !!el.querySelector("svg");
    if (svg) {
      try {
        snapSaveSvg(svg, baseName + ".svg");
        flashOutline(el, true);
        setStatus(statusEl, "Downloaded " + baseName + ".svg", false);
      } catch (e) {
        flashOutline(el, false);
        setStatus(statusEl, "SVG export failed (see console)", true);
        window.console && console.error("[screen-snaps]", e);
      }
      return;
    }

    var ppiSelect = document.getElementById("snaptool-ppi");
    var ppi = (ppiSelect && parseInt(ppiSelect.value, 10)) || 300;

    setStatus(statusEl, hasSvg ? "Mixed HTML + SVG — rendering PNG of everything…" : "No <svg> found — rendering PNG…", false);
    snapSavePng(el, baseName + ".png", ppi, function (err) {
      if (err) {
        flashOutline(el, false);
        setStatus(statusEl, "PNG export failed (see console)", true);
        window.console && console.error("[screen-snaps]", err);
      } else {
        flashOutline(el, true);
        setStatus(statusEl, "Downloaded " + baseName + ".png (" + ppi + " ppi)" +
          (hasSvg ? " (enter the SVG's own id for a vector-only .svg)" : ""), false);
      }
    });
  }

  function setStatus(statusEl, message, isError) {
    statusEl.textContent = message;
    statusEl.style.color = isError ? "#e74c3c" : "#3f5368";
  }

  // ---------------------------------------------------------------------
  // UI -- floating pill button + lightweight (non-blocking) popover.
  // ---------------------------------------------------------------------
  function injectStyles() {
    var style = document.createElement("style");
    style.textContent =
      "#snaptool-btn{position:fixed;right:20px;bottom:20px;z-index:100000;" +
      "background:#1f78b4;color:#fff;border:none;border-radius:999px;" +
      "padding:10px 16px;font-size:13px;font-weight:600;box-shadow:0 4px 14px rgba(0,0,0,.25);" +
      "cursor:pointer;font-family:Arial,sans-serif;}" +
      "#snaptool-btn:hover{background:#1a659e;}" +
      "#snaptool-panel{position:fixed;right:20px;bottom:66px;z-index:100000;width:280px;" +
      "background:#fff;border:1px solid #d9e2ef;border-radius:10px;" +
      "box-shadow:0 8px 24px rgba(0,0,0,.2);padding:12px 14px;font-family:Arial,sans-serif;" +
      "font-size:13px;color:#3f5368;display:none;}" +
      "#snaptool-panel.open{display:block;}" +
      "#snaptool-panel h4{margin:0 0 8px;font-size:13px;font-weight:700;color:#213547;}" +
      "#snaptool-quickpicks{display:flex;flex-wrap:wrap;gap:6px;margin-bottom:8px;}" +
      "#snaptool-quickpicks button{font-size:12px;padding:4px 8px;border-radius:6px;" +
      "border:1px solid #89abd3;background:#f3f8fe;color:#1f3347;cursor:pointer;}" +
      "#snaptool-quickpicks button:hover{background:#e3edfa;}" +
      "#snaptool-hint{font-size:11px;color:#6d8194;margin-bottom:8px;line-height:1.4;}" +
      "#snaptool-input{width:100%;box-sizing:border-box;padding:6px 8px;margin-bottom:8px;" +
      "border:1px solid #d9e2ef;border-radius:6px;font-size:13px;}" +
      "#snaptool-ppi-row{display:flex;align-items:center;gap:8px;margin-bottom:8px;font-size:12px;}" +
      "#snaptool-ppi{flex:1;padding:4px 6px;border:1px solid #d9e2ef;border-radius:6px;font-size:13px;}" +
      "#snaptool-go{width:100%;padding:6px 8px;background:#1f78b4;color:#fff;border:none;" +
      "border-radius:6px;font-size:13px;font-weight:600;cursor:pointer;}" +
      "#snaptool-go:hover{background:#1a659e;}" +
      "#snaptool-status{margin-top:8px;font-size:12px;min-height:14px;word-break:break-word;}";
    document.head.appendChild(style);
  }

  function renderQuickPicks(container) {
    container.innerHTML = "";
    var targets = document.querySelectorAll("." + QUICK_PICK_CLASS);
    if (!targets.length) return;

    targets.forEach(function (el, idx) {
      var label = el.getAttribute("data-snap-label") || el.id || ("Element " + (idx + 1));
      var btn = document.createElement("button");
      btn.type = "button";
      btn.textContent = label;
      btn.addEventListener("click", function () {
        var baseName = el.id || ("snap-" + (idx + 1));
        snapElement(el, baseName, document.getElementById("snaptool-status"));
      });
      container.appendChild(btn);
    });
  }

  function buildUi() {
    var btn = document.createElement("button");
    btn.id = "snaptool-btn";
    btn.type = "button";
    btn.textContent = "📸 Screen snaps";
    document.body.appendChild(btn);

    var panel = document.createElement("div");
    panel.id = "snaptool-panel";
    panel.innerHTML =
      '<h4>Screen snaps</h4>' +
      '<div id="snaptool-quickpicks"></div>' +
      '<div id="snaptool-hint">Real .svg if the element contains a chart, otherwise a ' +
      'high-res .png at the chosen ppi (tables/UI can\'t become real vector SVG).</div>' +
      '<input type="text" id="snaptool-input" placeholder="Element id, e.g. gpcrListTable_wrapper">' +
      '<div id="snaptool-ppi-row"><label for="snaptool-ppi">PNG resolution</label>' +
      '<select id="snaptool-ppi"><option value="300" selected>300 ppi</option>' +
      '<option value="600">600 ppi</option></select></div>' +
      '<button type="button" id="snaptool-go">Snap</button>' +
      '<div id="snaptool-status"></div>';
    document.body.appendChild(panel);

    // A page's own Bootstrap dropdowns (e.g. the DataMapper's "Colors" panel)
    // auto-close on any click that bubbles to `document` -- including clicks
    // on our own button/panel, which live outside that dropdown's DOM. Same
    // stopPropagation guard already used elsewhere in this codebase for
    // clicks *inside* a dropdown-menu (classification/wheel.js:6,
    // Mapper_GPCRomeWheel.html:297-299) -- here applied to our own UI so
    // interacting with Screen Snaps never counts as an "outside click" and
    // never closes whatever dropdown the user has open.
    panel.addEventListener("click", function (e) { e.stopPropagation(); });

    var input = document.getElementById("snaptool-input");
    var statusEl = document.getElementById("snaptool-status");

    btn.addEventListener("click", function (e) {
      e.stopPropagation();
      var willOpen = !panel.classList.contains("open");
      panel.classList.toggle("open", willOpen);
      if (willOpen) {
        renderQuickPicks(document.getElementById("snaptool-quickpicks"));
        input.focus();
      }
    });

    function runIdLookup() {
      var id = input.value.trim();
      input.style.borderColor = "#d9e2ef";
      if (!id) return;

      var el = document.getElementById(id);
      if (!el) {
        input.style.borderColor = "#e74c3c";
        setStatus(statusEl, "Element not found: " + id, true);
        return;
      }
      snapElement(el, id, statusEl);
    }

    document.getElementById("snaptool-go").addEventListener("click", runIdLookup);
    input.addEventListener("keydown", function (e) {
      if (e.key === "Enter") { runIdLookup(); }
    });
  }

  function init() {
    injectStyles();
    buildUi();
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();

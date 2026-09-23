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

  function snapSavePng(targetEl, filename, onDone) {
    if (typeof window.html2canvas !== "function") {
      onDone(new Error("html2canvas is not loaded"));
      return;
    }
    window.html2canvas(targetEl, { scale: 3, backgroundColor: "#ffffff", useCORS: true })
      .then(function (canvas) {
        triggerDownload(canvas.toDataURL("image/png"), filename, false);
        onDone(null);
      })
      .catch(function (err) { onDone(err); });
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

  function snapElement(el, baseName, statusEl) {
    var svg = (el.tagName && el.tagName.toLowerCase() === "svg") ? el : el.querySelector("svg");
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

    setStatus(statusEl, "No <svg> found — rendering PNG…", false);
    snapSavePng(el, baseName + ".png", function (err) {
      if (err) {
        flashOutline(el, false);
        setStatus(statusEl, "PNG export failed (see console)", true);
        window.console && console.error("[screen-snaps]", err);
      } else {
        flashOutline(el, true);
        setStatus(statusEl, "Downloaded " + baseName + ".png", false);
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
      'high-res .png (tables/UI can\'t become real vector SVG).</div>' +
      '<input type="text" id="snaptool-input" placeholder="Element id, e.g. gpcrListTable_wrapper">' +
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

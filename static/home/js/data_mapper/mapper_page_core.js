/**
 * Data Mapper — shared core utilities used by the per-page controllers
 * (mapper_tree_page.js, mapper_cluster_page.js, mapper_list_page.js,
 * mapper_heatmap_page.js, mapper_gpcrome_page.js).
 *
 * Deliberately small: only pieces that are genuinely identical in behavior
 * across pages live here. Row-management, autocomplete and snapshot/restore
 * logic differ in real ways per page (different fields count as "blank",
 * different label-type sources, species handling on some pages only) and
 * are kept in each page's own file rather than forced into a one-size-fits
 * config to avoid behavior changes that can't be verified without a browser.
 */
(function (window) {
  'use strict';

  /**
   * FNV-1a-style 32-bit string hash (shift-based variant). This exact
   * implementation was already shared verbatim by the tree, list and
   * GPCRome-wheel pages (as mapperWheelFNV1a32) — Cluster carried a different
   * multiplicative FNV-1a variant with the same seed; it now uses this one
   * too, so a given label's default colour is consistent across all 5 pages.
   */
  function fnv1a32(str) {
    var h = 0x811c9dc5;
    var s = String(str || '');
    var i, code;
    for (i = 0; i < s.length; i++) {
      code = s.charCodeAt(i);
      h ^= code;
      h += (h << 1) + (h << 4) + (h << 7) + (h << 8) + (h << 24);
      h = h >>> 0;
    }
    return h >>> 0;
  }

  /** Deterministic default HSL colour for a categorical label, derived from fnv1a32. */
  function defaultColorForLabel(lbl) {
    var hue = fnv1a32(String(lbl || '').toLowerCase()) % 360;
    return 'hsl(' + hue + ', 68%, 48%)';
  }

  /**
   * Returns a schedule() function that debounces calls to fn by ms.
   * Each page keeps its own suppressRedraw flag/checks and calls
   * schedule() from its own scheduleRedraw() wrapper.
   */
  function debounce(fn, ms) {
    var timer = null;
    return {
      schedule: function () {
        window.clearTimeout(timer);
        timer = window.setTimeout(fn, ms);
      },
      cancel: function () {
        window.clearTimeout(timer);
      }
    };
  }

  /**
   * Classifies clipboard text pasted into an input table, using the GPCRome
   * Wheel's rules (mapper_gpcrome_page.js paste handler):
   *   'grid'   – a TAB anywhere, or every non-empty line contains ';'
   *   'column' – no column delimiter but a line break: one value per line
   *   null     – a plain single-line value; leave it to the browser
   * `lines` holds the non-empty lines.
   */
  function classifyPaste(text) {
    var s = String(text || '');
    var lines = s.split(/\r?\n/).filter(function (ln) { return ln.trim() !== ''; });
    if (!lines.length) return { kind: null, lines: lines };
    var semiGrid = lines.every(function (ln) { return ln.indexOf(';') !== -1; });
    if (s.indexOf('\t') !== -1 || semiGrid) return { kind: 'grid', lines: lines };
    if (/[\r\n]/.test(s)) return { kind: 'column', lines: lines };
    return { kind: null, lines: lines };
  }

  /** Splits one pasted grid line into trimmed cells: TAB columns, else ';' columns. */
  function splitPasteLine(line) {
    var s = String(line || '');
    var parts = s.indexOf('\t') !== -1 ? s.split('\t') : s.indexOf(';') !== -1 ? s.split(';') : [s];
    return parts.map(function (p) { return p.trim(); });
  }

  /**
   * Writes one pasted line per row down the column of the cell that received
   * the paste, starting at that row (or at findBlankRow() in the receptor
   * column when the paste didn't come from a row). Receptor cells go through
   * setReceptor($tr, raw); other cells through their first visible, enabled
   * text input, then the optional onValue($input) hook. appendRow() must return
   * the new row, or an empty jQuery object once the page's row cap is reached
   * (remaining lines are dropped).
   * Returns false, touching nothing, when the target cell has nothing writable.
   */
  function pasteColumnDown(opts) {
    var $ = window.jQuery;
    var $tbody = opts.$tbody;
    var $target = $(opts.target);
    var $tr0 = $target.closest('tr');
    var inBody = $tr0.length && $.contains($tbody[0], $tr0[0]);
    var $td = inBody ? $target.closest('td') : $();
    var isReceptor = !inBody || $td.hasClass('mapper-core-receptor-cell');
    var colIdx = $td.length ? $td.index() : -1;
    var writableSel = 'input[type="text"]:visible:not(:disabled):not([readonly])';
    if (!isReceptor && !$td.find(writableSel).length) return false;

    var $start = inBody ? $tr0 : opts.findBlankRow();
    var $tr = $start && $start.length ? $start : opts.appendRow();
    for (var i = 0; i < opts.lines.length && $tr && $tr.length; i++) {
      var val = String(opts.lines[i]).trim();
      if (isReceptor) {
        opts.setReceptor($tr, val);
      } else {
        var $inp = $tr.children().eq(colIdx).find(writableSel).first().val(val);
        if (opts.onValue) opts.onValue($inp);
      }
      if (i + 1 < opts.lines.length) {
        var $next = $tr.next('tr');
        $tr = $next.length ? $next : opts.appendRow();
      }
    }
    return true;
  }

  /**
   * Receptor text for a download. UniProt mode writes the UniProt entry name, which pastes
   * straight back: 5HT1A_HUMAN, or for a no-human-ortholog row its species entry (5HT5B_MOUSE),
   * or the ortholog picked in the row's Species dropdown. Other modes write opts.label(entryId),
   * the page's cell text for the active Receptor Names mode.
   */
  function downloadReceptorLabel(entryId, $tr, opts) {
    if (!opts.uniprotMode) return opts.label(entryId);
    if (/_human$/i.test(entryId)) return entryId.toUpperCase();
    var spec = opts.speciesSelect ? ($tr.find(opts.speciesSelect).val() || '').trim() : '';
    return (spec || entryId).toUpperCase();
  }

  /**
   * Saves a page's input table as .xlsx: one row per receptor row, unmatched receptors as typed.
   * opts: tbody (selector), headers, sheetName, fileName, values($tr) → the row's value cells,
   *       label(entryId) and uniprotMode / speciesSelect (see downloadReceptorLabel).
   */
  function downloadInputTable(opts) {
    var $ = window.jQuery;
    var meta = window.MAPPER_CORE_ENTRY_META || {};
    var rows = [opts.headers];
    $(opts.tbody + ' tr').each(function () {
      var $tr = $(this);
      var entryId = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
      var unmatched = $tr.data('mapperCoreUnmatchedRaw');
      var receptor;
      if (entryId && meta[entryId]) {
        receptor = downloadReceptorLabel(entryId, $tr, opts) || meta[entryId].name_plain || meta[entryId].text || entryId;
      } else if (unmatched != null && String(unmatched).trim()) {
        receptor = String(unmatched).trim();
      } else {
        receptor = ($tr.find('.mapper-core-in-receptor').val() || '').trim();
      }
      if (!receptor) return;
      rows.push([receptor].concat(opts.values($tr)));
    });
    var wb = window.XLSX.utils.book_new();
    window.XLSX.utils.book_append_sheet(wb, window.XLSX.utils.aoa_to_sheet(rows), opts.sheetName);
    window.XLSX.writeFile(wb, opts.fileName);
  }

  window.MapperPageCore = {
    fnv1a32: fnv1a32,
    defaultColorForLabel: defaultColorForLabel,
    debounce: debounce,
    classifyPaste: classifyPaste,
    splitPasteLine: splitPasteLine,
    pasteColumnDown: pasteColumnDown,
    downloadInputTable: downloadInputTable
  };
})(window);

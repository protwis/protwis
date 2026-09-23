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

  window.MapperPageCore = {
    fnv1a32: fnv1a32,
    defaultColorForLabel: defaultColorForLabel,
    debounce: debounce
  };
})(window);

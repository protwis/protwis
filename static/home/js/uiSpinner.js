// uiSpinner.js — tiny, framework-friendly spinner helpers
// Works with jQuery or plain DOM. Also includes a DataTables helper.

(function (global) {
  const defaultCfg = {
    spinnerSelector: '.btn-spinner',
    spinnerHTML: '<span class="btn-spinner" aria-hidden="true"></span>',
    hostClass: 'has-spinner',
    busyClass: 'is-busy',
    hideDelayMs: 150,           // give the eye a beat before hiding
    posAttr: 'data-spinner-pos' // e.g., "center" (default), "right", "inline"
  };

  function isJQ(x) { return !!(x && x.jquery); }
  function toJQ(x) { return isJQ(x) ? x : (x instanceof Element ? window.jQuery(x) : window.jQuery()); }

  function ensureHostClass($btn, cfg) {
    if (!$btn.hasClass(cfg.hostClass)) $btn.addClass(cfg.hostClass);
  }

  function ensureSpinner($btn, cfg) {
    if (!$btn.children(cfg.spinnerSelector).length) {
      $btn.append(cfg.spinnerHTML);
    }
  }

  /**
   * Show a spinner on a button-like element.
   * @param {HTMLElement|jQuery} btn
   * @param {Object} cfg
   */
  function showBtnSpinner(btn, cfg = {}) {
    const opt = { ...defaultCfg, ...cfg };
    const $btn = toJQ(btn);
    if (!$btn.length) return;

    ensureHostClass($btn, opt);
    ensureSpinner($btn, opt);

    // Accessibility
    $btn.attr('aria-busy', 'true');

    $btn.addClass(opt.busyClass);
  }

  /**
   * Hide spinner/busy state.
   * @param {HTMLElement|jQuery} btn
   * @param {Object} cfg
   */
  function hideBtnSpinner(btn, cfg = {}) {
    const opt = { ...defaultCfg, ...cfg };
    const $btn = toJQ(btn);
    if (!$btn.length) return;
    $btn.removeClass(opt.busyClass).removeAttr('aria-busy');
  }

  /**
   * Run a function while showing a spinner on a host element.
   * Prevents concurrent runs if already busy.
   * @param {HTMLElement|jQuery} btn
   * @param {() => any|Promise<any>} workFn
   * @param {Object} cfg
   * @returns {Promise<any>}
   */
  function runWithButtonSpinner(btn, workFn, cfg = {}) {
    const opt = { ...defaultCfg, ...cfg };
    const $btn = toJQ(btn);
    if (!$btn.length) return Promise.resolve();

    // If already busy, ignore this click
    if ($btn.hasClass(opt.busyClass)) return Promise.resolve();

    showBtnSpinner($btn, opt);

    // Double rAF so spinner paints before heavy work starts
    const yieldTwoFrames = () =>
      new Promise(r => requestAnimationFrame(() => requestAnimationFrame(r)));

    return yieldTwoFrames()
      .then(() => Promise.resolve().then(workFn))
      .finally(() => setTimeout(() => hideBtnSpinner($btn, opt), opt.hideDelayMs));
  }

  /**
   * DataTables convenience: find a named button on this table and run with spinner.
   * @param {DataTables.Api} dt
   * @param {string} buttonSelector e.g. 'exportExcel:name'
   * @param {() => any|Promise<any>} workFn
   * @param {Object} cfg
   */
  function runWithDtButtonSpinner(dt, buttonSelector, workFn, cfg = {}) {
    let node = null;
    try {
      const api = dt && typeof dt.button === 'function' ? dt.button(buttonSelector) : null;
      node = api ? api.node() : null;
    } catch (_) { /* noop */ }

    if (!node) return Promise.resolve().then(workFn); // Fallback: do the work, no spinner
    return runWithButtonSpinner(node, workFn, cfg);
  }

  // Public API
  const UiSpinner = {
    showBtnSpinner,
    hideBtnSpinner,
    runWithButtonSpinner,
    runWithDtButtonSpinner,
  };

  if (typeof module !== 'undefined' && module.exports) module.exports = UiSpinner;
  else global.UiSpinner = UiSpinner;

})(window);

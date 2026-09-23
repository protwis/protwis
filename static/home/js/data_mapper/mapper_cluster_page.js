/**
 * Mapper 2.0 — Cluster page.
 * Left panel: receptor + gradient (numbers) / category (categories) input table.
 * Right panel: Plotly scatter plot using tSNE positions from backend.
 *
 * IMPORTANT: The globals below (currentClusterData, labelsVisible, cluster_DataStyling,
 * colorPalette, getActiveColorOption) are declared at window scope because datamapper.js
 * references them directly from updatePlotWithAnnotations(), createTraces(), createAnnotations().
 */

// ── Globals required by datamapper.js ──────────────────────────────────────
var currentClusterData = [];
var labelsVisible = false;
// Tracks user zoom: null = use data range, [[x0,x1],[y0,y1]] = preserve zoom
var _clusterZoomRange = null;
var colorPalette = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
    '#393b79', '#e7969c', '#17becf', '#98df8a', '#ffbb78',
    '#9edae5', '#c5b0d5', '#f7b6d2', '#dbdb8d', '#c49c94'
];
var cluster_DataStyling = {
    labelFontSize: 14,
    markerSize: 10,
    strokeWidth: 0.5,
    borderOn: true,
    textColorEnabled: true
};

// ── Custom color state (read by datamapper.js) ──────────────────────────────
// Gradient: 3-stop min/mid/max colors
var CLUSTER_GRADIENT_COLORS = { min: '#3c5488', mid: '#ffffff', max: '#e64b35' };
// Gradient stop count: 1=max-only (white→max), 2=min+max, 3=min+mid+max
var CLUSTER_GRADIENT_STOPS = 3;
// Per-label colors for all categorical schemes; key = label text (e.g. 'Cluster 1', 'Class A', 'My group')
var CLUSTER_LABEL_COLORS = {};
var _lastClusterEntries = '';  // fingerprint to detect receptor set changes

var CLUSTER_COLOR_PALETTE = [
    ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f'],
    ['#3c5488','#e64b35','#00a087','#4dbbd5','#f39b7f','#b09c85','#91d1c2','#dc0000'],
    ['#000000','#ffffff','#c0c0c0','#808080','#404040','#a0522d','#556b2f','#191970']
];

function getActiveColorOption() {
    var opts = [
        ['cluster-colorByCluster',       'cluster'],
        ['cluster-colorByGradient',      'gradient'],
        ['cluster-colorByUserCategory',  'userCategory'],
        ['cluster-colorByClass',         'Class'],
        ['cluster-colorByLigandType',    'Ligand type'],
        ['cluster-colorByReceptorFamily','Receptor family']
    ];
    for (var i = 0; i < opts.length; i++) {
        var btn = document.getElementById(opts[i][0]);
        if (btn && btn.classList.contains('active')) return opts[i][1];
    }
    return 'cluster';
}
// ── Color helpers (window-scoped so datamapper.js can call them) ────────────
// FNV1a32 hash + default label colour moved to MapperPageCore (mapper_page_core.js) —
// shared verbatim with the tree/list/GPCRome-wheel pages.
function mapperClusterGetLabelColor(lbl) {
    if (lbl == null || lbl === '') return '#cccccc';
    if (!CLUSTER_LABEL_COLORS[lbl]) CLUSTER_LABEL_COLORS[lbl] = MapperPageCore.defaultColorForLabel(lbl);
    return CLUSTER_LABEL_COLORS[lbl];
}
// ───────────────────────────────────────────────────────────────────────────

(function ($) {
  'use strict';

  var MAPPER_CLUSTER_MAX_ROWS  = 500;
  var MAPPER_CLUSTER_PLACEHOLDER = 'Paste in 1–2 columns\nor type';
  var DEBOUNCE_MS = 400;

  var posCheckTimer     = null;   // debounce for position input checks
  var suppressRedraw    = false;
  var plotMode          = 'numbers';  // 'numbers' | 'categories'
  var MAPPER_CLUSTER_MODE_SNAPSHOTS = { numbers: null, categories: null };
  // Saved color button per mode — restored when returning to a mode.
  var MAPPER_CLUSTER_COLOR_STATE    = { numbers: 'cluster-colorByGradient', categories: 'cluster-colorByUserCategory' };
  var inPositionMode    = false;
  var positionProcessed = false;
  var currentPositionXhr = null;
  var processTimer      = null;
  // Per-mode tSNE result cache so mode switches can restore without re-running AJAX
  var positionResultCache = { numbers: null, categories: null };
  var _plotlyFrozen = false;    // true while position overlay is visible
  // Counters updated incrementally so mapperClusterCheckPositionMode needs no DOM scan
  var _posCtReceptors = 0;   // rows that have a resolved receptor
  var _posCtFilled    = 0;   // rows that have both a receptor AND a non-empty position
  var labelDict_IUPHAR = {};
  var labelDict_Gene   = {};

  // ── Colors panel ─────────────────────────────────────────────────────────

  function mapperClusterRebuildColorsPanel() {
    var scheme = getActiveColorOption();
    var isGrad = (scheme === 'gradient');
    var $panel = $('#cluster-colors-panel');
    $panel.toggleClass('mapper-core-colors-show-numeric', isGrad)
          .toggleClass('mapper-core-colors-show-text', !isGrad);
    $('#cluster-cp-empty').hide();
    if (isGrad) return; // gradient pickers already initialised at boot

    // Build label list for categorical schemes
    var labels = [];
    if (!currentClusterData.length) {
      $('#cluster-cp-empty').show(); return;
    }
    if (scheme === 'cluster') {
      var clusterIds = Array.from(new Set(currentClusterData.map(function(d){ return d.cluster; }))).sort(function(a,b){ return a-b; });
      labels = clusterIds.map(function(c){ return 'Cluster ' + (c+1); });
    } else {
      var prop = scheme; // 'Class','Ligand type','Receptor family','userCategory'
      labels = Array.from(new Set(currentClusterData.map(function(d){
        return d[prop];
      }))).filter(Boolean);
      labels.sort();
    }

    var titleMap = { cluster:'Cluster colors', Class:'Class colors', 'Ligand type':'Chemotype colors', 'Receptor family':'Receptor family colors', userCategory:'Category colors' };

    // Smart rebuild: if labels unchanged, just update swatch colors (keeps open Spectrum popups)
    var $list = $('#cluster-cp-label-list');
    var $existing = $list.find('.mapper-cluster-lcat-inp');
    var existingLabels = $existing.map(function(){ return $(this).data('cp-label'); }).get();
    if (JSON.stringify(existingLabels) === JSON.stringify(labels)) {
      $existing.each(function() {
        var lbl = $(this).data('cp-label');
        try { $(this).spectrum('set', mapperClusterGetLabelColor(lbl)); } catch(e) {}
      });
      return;
    }

    // Full rebuild
    $existing.each(function(){ try { $(this).spectrum('destroy'); } catch(e) {} });
    $list.empty();
    // Scheme header
    $list.append($('<div class="mapper-cluster-lcat-scheme-head">').text(titleMap[scheme] || scheme));
    labels.forEach(function(lbl) {
      var hex = mapperClusterGetLabelColor(lbl);
      var $row = $('<div class="color-item">');
      var $name = $('<label class="color-label">').text(lbl).attr('title', lbl);
      var $inp = $('<input type="text" class="mapper-cluster-lcat-inp">').data('cp-label', lbl).hide();
      $row.append($name, $inp);
      $list.append($row);
      (function(label, initHex) {
        $inp.spectrum({
          color: initHex,
          showPalette: true, showInput: true, showButtons: false, preferredFormat: 'hex',
          appendTo: '#cluster-colors-panel',
          containerClassName: 'mapper-core-lcat-sp-container',
          replacerClassName: 'mapper-core-lcat-replacer',
          palette: CLUSTER_COLOR_PALETTE,
          change: function(c) {
            CLUSTER_LABEL_COLORS[label] = c ? c.toHexString() : initHex;
            mapperClusterSyncRowSwatchesByLabel(label);
            mapperClusterScheduleRedraw();
          },
          move: function(c) {
            CLUSTER_LABEL_COLORS[label] = c ? c.toHexString() : initHex;
            mapperClusterSyncRowSwatchesByLabel(label);
            mapperClusterScheduleRedraw();
          }
        });
      }(lbl, hex));
    });
  }

  // Sync all row swatches that share the given category label
  function mapperClusterSyncRowSwatchesByLabel(lbl) {
    var hex = mapperClusterGetLabelColor(lbl);
    $('#mapper-cluster-input-tbody tr').each(function() {
      var rowCat = $(this).find('.mapper-cluster-cat').val().trim();
      if (rowCat === lbl) {
        try { $(this).find('.mapper-cluster-label-swatch').spectrum('set', hex); } catch(e) {}
      }
    });
  }

  // Demo data (subset from heatmap demo, gradient only)
  var MAPPER_CLUSTER_DEMO_ROWS = [
    { receptor: '5HT1A', gradient: '1',  category: 'Group A' },
    { receptor: '5HT1B', gradient: '2',  category: 'Group B' },
    { receptor: '5HT2A', gradient: '6',  category: 'Group C' },
    { receptor: 'ACKR1', gradient: '13', category: 'Group A' },
    { receptor: 'ACKR2', gradient: '14', category: 'Group B' },
    { receptor: 'ACM1',  gradient: '17', category: 'Group C' },
    { receptor: 'ACM2',  gradient: '18', category: 'Group A' },
    { receptor: 'ADA1A', gradient: '23', category: 'Group B' },
    { receptor: 'ADRB1', gradient: '29', category: 'Group C' },
    { receptor: 'ADRB2', gradient: '30', category: 'Group A' }
  ];

  // ── Helpers ────────────────────────────────────────────────────────────────
  function parseNumLoose(s) {
    if (s == null || String(s).trim() === '') return null;
    var x = parseFloat(String(s).trim().replace(',', '.'));
    return isFinite(x) ? x : null;
  }

  // ── Label converter ────────────────────────────────────────────────────────
  function mapperClusterBuildLabelConverter() {
    var meta = window.MAPPER_CORE_ENTRY_META || {};
    Object.keys(meta).forEach(function (entry_name) {
      var m = meta[entry_name];
      var stripped = entry_name.toLowerCase().replace('_human', '');
      if (m.name_plain) {
        labelDict_IUPHAR[entry_name.toLowerCase()] = m.name_plain;
        labelDict_IUPHAR[stripped]                 = m.name_plain;
      }
      if (m.gene) {
        labelDict_Gene[entry_name.toLowerCase()] = m.gene;
        labelDict_Gene[stripped]                 = m.gene;
      }
    });
  }

  // ── Data construction ──────────────────────────────────────────────────────
  function mapperClusterBuildData() {
    var data = {};
    $('#mapper-cluster-input-tbody tr').each(function () {
      var $tr   = $(this);
      var entry = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
      if (!entry) {
        var ta  = ($tr.find('.mapper-core-in-receptor').val() || '').trim();
        var up  = ta.toUpperCase();
        if (up && window.MAPPER_CORE_RESOLVE) entry = window.MAPPER_CORE_RESOLVE[up] || '';
      }
      if (!entry) return;

      if (plotMode === 'numbers') {
        var v1 = parseNumLoose($tr.find('.mapper-cluster-gradient').val());
        data[entry] = { Value1: v1 != null ? v1 : 0 };
      } else {
        var cat = ($tr.find('.mapper-cluster-cat').val() || '').trim();
        data[entry] = { Value1: 0, userCategory: cat || 'Unknown' };
      }
    });
    return data;
  }

  // ── Placeholder ────────────────────────────────────────────────────────────
  function showPlaceholder() {
    if (typeof Plotly !== 'undefined') Plotly.purge('plotContainer_cluster');
    else $('#plotContainer_cluster').empty();
    $('#mapper-cluster-placeholder').show();
  }
  function hidePlaceholder() {
    $('#mapper-cluster-placeholder').hide();
  }

  // ── Position lookup (built once from pre-loaded data) ─────────────────────
  var _clusterPosMap = null;
  function getClusterPosMap() {
    if (_clusterPosMap) return _clusterPosMap;
    _clusterPosMap = {};
    var allPos = window.MAPPER_CLUSTER_ALL_POSITIONS;
    if (allPos && allPos.length) {
      allPos.forEach(function (pos) {
        if (pos.label) _clusterPosMap[pos.label.toLowerCase()] = pos;
      });
    }
    return _clusterPosMap;
  }

  // ── Multi-class warning ────────────────────────────────────────────────────
  function updateMultiClassWarning(data) {
    var classMap = window.MAPPER_CLUSTER_CLASS_MAP || {};
    var seen = {};
    Object.keys(data || {}).forEach(function (entry) {
      var cls = classMap[entry];
      if (cls) seen[cls] = 1;
    });
    var isMulti = Object.keys(seen).length > 1;
    $('#cluster-multiclass-warning').toggle(isMulti);
  }

  // ── Plot render (client-side, no AJAX) ────────────────────────────────────
  function mapperClusterRenderPlot() {
    if (suppressRedraw) return;
    var data = mapperClusterBuildData();
    updateMultiClassWarning(data);
    var count = Object.keys(data).length;
    if (count < 2) { showPlaceholder(); return; }

    var posMap = getClusterPosMap();
    if (!Object.keys(posMap).length) {
      $('#mapper-cluster-messages').text('Position data unavailable. Please refresh the page.');
      showPlaceholder();
      return;
    }

    var newClusterData = [];
    Object.keys(data).forEach(function (entry_name) {
      var stripped = entry_name.toLowerCase().replace(/_human$/, '');
      var pos = posMap[stripped];
      if (!pos) return;
      var point = {
        x:                 pos.x,
        y:                 pos.y,
        label:             stripped,
        original_label:    stripped,
        'Class':           pos['Class'] || '',
        'Ligand type':     pos['Ligand type'] || '',
        'Receptor family': pos['Receptor family'] || '',
        fill:              data[entry_name].Value1 != null ? data[entry_name].Value1 : 0,
        cluster:           0
      };
      if (plotMode === 'categories') {
        point.userCategory = data[entry_name].userCategory || 'Unknown';
      }
      newClusterData.push(point);
    });

    if (newClusterData.length < 2) { showPlaceholder(); return; }

    currentClusterData = newClusterData;
    _clusterZoomRange = null;
    var _liveEl = document.getElementById('plotContainer_cluster');
    if (_liveEl) { _liveEl.style.width = ''; _liveEl.style.height = ''; }
    hidePlaceholder();
    $('#mapper-cluster-messages').text('');

    var maxClusters = Math.min(20, currentClusterData.length - 1);
    var slider = document.getElementById('mapper-cluster-slider');
    if (slider) {
      slider.max = maxClusters;
      if (parseInt(slider.value, 10) > maxClusters) {
        slider.value = Math.min(5, maxClusters);
      }
      document.getElementById('mapper-cluster-slider-val').textContent = slider.value;
    }

    var _entFP = currentClusterData.map(function(d){ return d.original_label; }).sort().join(',');
    var isNewData = _entFP !== _lastClusterEntries;
    _lastClusterEntries = _entFP;
    var numClusters = Math.min(5, currentClusterData.length - 1);
    performClustering(numClusters, isNewData);
  }

  var _clusterRedrawDebounced = MapperPageCore.debounce(function () { mapperClusterRenderPlot(); }, DEBOUNCE_MS);
  function mapperClusterScheduleRedraw() {
    if (inPositionMode && positionProcessed) {
      mapperClusterCosmeticRedraw();
      return;
    }
    if (inPositionMode && !positionProcessed) {
      return;
    }
    _clusterRedrawDebounced.schedule();
  }

  // ── Clustering ─────────────────────────────────────────────────────────────
  function setActiveColorButton(activeButtonId) {
    MAPPER_CLUSTER_COLOR_STATE[plotMode] = activeButtonId;
    var colorOptions = [
      'cluster-colorByCluster', 'cluster-colorByGradient', 'cluster-colorByUserCategory',
      'cluster-colorByClass', 'cluster-colorByLigandType', 'cluster-colorByReceptorFamily'
    ];
    colorOptions.forEach(function (option) {
      var btn = document.getElementById(option);
      if (!btn) return;
      if (option === activeButtonId) {
        btn.classList.add('active', 'btn-primary');
        btn.classList.remove('btn-outline-primary');
      } else {
        btn.classList.remove('active', 'btn-primary');
        btn.classList.add('btn-outline-primary');
      }
    });
    // Update Colors panel sub-panel if it's open
    if ($('#cluster-colors-btn').closest('.dropdown').hasClass('open')) mapperClusterRebuildColorsPanel();
  }

  function performClustering(numClusters, setDefaultColor) {
    if (!currentClusterData.length) return;
    var coordinates = currentClusterData.map(function (d) { return [d.x, d.y]; });
    if (typeof ML !== 'undefined' && typeof ML.KMeans === 'function') {
      try {
        var options = { maxIterations: 100, tolerance: 1e-6, initialization: 'kmeans++' };
        var result  = ML.KMeans(coordinates, numClusters, options);
        currentClusterData.forEach(function (point, idx) {
          point.cluster = result.clusters[idx];
        });
      } catch (e) {}
    }
    // Color default and render happen unconditionally so they are never
    // silently skipped when ML is unavailable or KMeans throws.
    if (setDefaultColor) {
      setActiveColorButton(plotMode === 'numbers' ? 'cluster-colorByGradient' : 'cluster-colorByUserCategory');
    }
    switchClusterLabels(getActiveLabelType());
  }

  // ── Label switching ────────────────────────────────────────────────────────
  function switchClusterLabels(labelType) {
    currentClusterData.forEach(function (point) {
      point.label = point.original_label.toUpperCase();
    });
    var dict = null;
    if (labelType === 'IUPHAR') dict = labelDict_IUPHAR;
    else if (labelType === 'Entrez') dict = labelDict_Gene;

    if (dict) {
      currentClusterData.forEach(function (point) {
        var key = point.original_label.toLowerCase() + '_human';
        if (dict[key]) point.label = dict[key];
        else if (dict[point.original_label.toLowerCase()]) point.label = dict[point.original_label.toLowerCase()];
      });
    }
    updatePlotWithAnnotations();
  }

  function getActiveLabelType() {
    var active = document.querySelector('.cluster-leaf-btn.btn-primary');
    return active ? active.getAttribute('data-value') : 'IUPHAR';
  }

  // ── Mode snapshot helpers ──────────────────────────────────────────────────
  function mapperClusterCaptureSnapshot() {
    var rows = [];
    $('#mapper-cluster-input-tbody tr').each(function () {
      var $tr      = $(this);
      var entry    = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
      var rawText  = ($tr.find('.mapper-core-in-receptor').val()    || '').trim();
      if (!entry && !rawText) return; // skip the trailing blank row
      var row = { entry: entry, receptorText: rawText, pos: ($tr.find('.mapper-cluster-pos').val() || '').trim() };
      if (plotMode === 'numbers') {
        row.gradient = ($tr.find('.mapper-cluster-gradient').val() || '').trim();
      } else {
        row.category = ($tr.find('.mapper-cluster-cat').val() || '').trim();
      }
      rows.push(row);
    });
    MAPPER_CLUSTER_MODE_SNAPSHOTS[plotMode] = {
      rows:     rows,
      posCache: positionResultCache[plotMode] || null
    };
  }

  function mapperClusterRestoreSnapshot(targetMode) {
    // Destroy existing AC widgets before clearing DOM
    $('#mapper-cluster-input-tbody tr').each(function () {
      var $ta = $(this).find('.mapper-core-in-receptor');
      try { if ($ta.data('ui-autocomplete')) $ta.autocomplete('destroy'); } catch (e) {}
    });
    $('#mapper-cluster-input-tbody').empty();

    var snap = MAPPER_CLUSTER_MODE_SNAPSHOTS[targetMode];
    if (snap && snap.rows && snap.rows.length) {
      snap.rows.forEach(function (r) {
        var rowOpts = { entry: r.entry };
        if (r.receptorText) rowOpts.receptorText = r.receptorText;
        if (targetMode === 'numbers') {
          rowOpts.gradient = r.gradient != null ? r.gradient : '';
        } else {
          rowOpts.category = r.category != null ? r.category : '';
        }
        if (r.pos) rowOpts.pos = r.pos;
        mapperClusterAppendRow(rowOpts);
      });
    }
    mapperClusterEnsureTrailingBlankRow();
    _posCountersReset();
  }

  // ── Input mode toggle ──────────────────────────────────────────────────────
  function setPlotMode(mode) {
    if (mode === plotMode) return;

    // Save current table state into the current mode's snapshot
    mapperClusterCaptureSnapshot();

    // Abort any in-flight position calculation and reset position state.
    // The snapshot already persisted the positions, so nothing is lost.
    clearTimeout(processTimer); processTimer = null;
    if (currentPositionXhr) { currentPositionXhr.abort(); currentPositionXhr = null; }
    inPositionMode = false;
    positionProcessed = false;
    hidePositionOverlay();

    plotMode = mode;
    var isNum = (mode === 'numbers');

    $('#mapper-cluster-numbers').toggleClass('active', isNum);
    $('#mapper-cluster-categories').toggleClass('active', !isNum);
    $('#mapper-cluster-val-header').text(isNum ? 'Gradient' : 'Category');

    // Rebuild table from the target mode's snapshot (keeps the two modes independent)
    suppressRedraw = true;
    mapperClusterRestoreSnapshot(mode);
    suppressRedraw = false;

    // Sync col visibility for any rows just restored
    $('#mapper-cluster-input-tbody tr').each(function () {
      $(this).find('.mapper-cluster-gradient').toggle(isNum);
      $(this).find('.mapper-cluster-cat').toggle(!isNum);
      // Init/destroy row swatch on mode switch
      if (!isNum) mapperClusterInitRowSwatch($(this));
      else {
        try { $(this).find('.mapper-cluster-label-swatch').spectrum('destroy'); } catch(e) {}
      }
    });

    // CSS-controlled: mapper-core-text-mode on the table shows/hides the swatch column
    $('#mapper-cluster-input-table').toggleClass('mapper-core-text-mode', !isNum);

    $('.mapper-cluster-numbers-only').toggle(isNum);
    $('.mapper-cluster-categories-only').toggle(!isNum);

    setActiveColorButton(MAPPER_CLUSTER_COLOR_STATE[mode]);

    mapperClusterSyncFirstRowPlaceholder();

    // Re-seed the live position cache from the snapshot if it was cleared since the
    // snapshot was taken.  The fingerprint check inside TryRestorePositionCache still
    // guards against stale geometry (different receptors/positions).
    var _snap = MAPPER_CLUSTER_MODE_SNAPSHOTS[mode];
    if (_snap && _snap.posCache && !positionResultCache[mode]) {
      positionResultCache[mode] = _snap.posCache;
    }

    // Try to restore a cached tSNE result for this mode (same receptors + positions).
    // If the cache is missing or stale, fall back to the normal position-check flow.
    if (!mapperClusterTryRestorePositionCache()) {
      mapperClusterCheckPositionMode();
    }
    // Schedule a redraw only when we are not blocked waiting for positions to be filled.
    if (!inPositionMode || positionProcessed) {
      mapperClusterScheduleRedraw();
    }

    mapperClusterSyncClearDropdown();
  }

  function mapperClusterSyncClearDropdown() {
    var isNumbers = (plotMode === 'numbers');
    $('#mapper-cluster-clear-mode-label').text(isNumbers ? 'Numbers' : 'Categories');
    $('.mapper-cluster-clear-col-num').toggle(isNumbers);
    $('.mapper-cluster-clear-col-text').toggle(!isNumbers);
  }

  // ── Position mode ──────────────────────────────────────────────────────────
  function _freezePlotly() {
    if (_plotlyFrozen) return;
    _plotlyFrozen = true;
    var gd = document.getElementById('plotContainer_cluster');
    if (gd && gd._fullLayout) {
      try { Plotly.relayout(gd, { hovermode: false }); } catch (e) {}
    }
  }
  function _unfreezePlotly() {
    if (!_plotlyFrozen) return;
    _plotlyFrozen = false;
    var gd = document.getElementById('plotContainer_cluster');
    if (gd && gd._fullLayout) {
      try { Plotly.relayout(gd, { hovermode: 'closest' }); } catch (e) {}
    }
  }
  function showPositionOverlayFilling() {
    $('#cluster-pos-info-box').show();
    $('#cluster-pos-fill-msg').show();
    $('#cluster-pos-processing-msg').hide();
    // No Plotly freeze during fill — plot remains interactive
  }
  function showPositionOverlayProcessing() {
    $('#cluster-pos-info-box').show();
    $('#cluster-pos-fill-msg').hide();
    $('#cluster-pos-processing-msg').show();
  }
  function hidePositionOverlay() {
    $('#cluster-pos-info-box').hide();
  }

  function mapperClusterCheckPositionMode() {
    var hasReceptor = _posCtReceptors > 0;
    var anyPos      = _posCtFilled    > 0;
    var allFilled   = hasReceptor && _posCtFilled === _posCtReceptors;

    var wasIn = inPositionMode;
    if (!hasReceptor || !anyPos) {
      clearTimeout(processTimer); processTimer = null;
      if (currentPositionXhr) { currentPositionXhr.abort(); currentPositionXhr = null; }
      inPositionMode = false;
      positionProcessed = false;
      hidePositionOverlay();
      if (wasIn) mapperClusterScheduleRedraw();
      return;
    }

    inPositionMode = true;

    if (!allFilled) {
      clearTimeout(processTimer); processTimer = null;
      if (currentPositionXhr) { currentPositionXhr.abort(); currentPositionXhr = null; }
      positionProcessed = false;
      showPositionOverlayFilling();
      return;
    }

    if (positionProcessed) {
      hidePositionOverlay();
      return;
    }

    // All filled, not yet processed — schedule auto-calculation
    mapperClusterScheduleProcess();
  }

  function mapperClusterScheduleProcess() {
    clearTimeout(processTimer);
    if (currentPositionXhr) { currentPositionXhr.abort(); currentPositionXhr = null; }
    showPositionOverlayProcessing();
    processTimer = setTimeout(mapperClusterDoPositionProcess, 500);
  }

  function getCsrfToken() {
    var el = document.querySelector('[name=csrfmiddlewaretoken]');
    if (el) return el.value;
    var m = document.cookie.match(/csrftoken=([^;]+)/);
    return m ? m[1] : '';
  }

  function mapperClusterCosmeticRedraw() {
    if (!currentClusterData.length) return;
    var inputData = mapperClusterBuildData();
    currentClusterData.forEach(function (point) {
      var key = point.label + '_human';
      var row = inputData[key] || inputData[point.label];
      if (!row) return;
      point.fill = row.Value1 != null ? row.Value1 : 0;
      if (plotMode === 'categories') point.userCategory = row.userCategory || 'Unknown';
    });
    updateMultiClassWarning(inputData);
    if (typeof updatePlotWithAnnotations === 'function') updatePlotWithAnnotations();
  }

  function mapperClusterDoPositionProcess() {
    var data = {};
    var allNumeric = true;
    $('#mapper-cluster-input-tbody tr').each(function () {
      var $tr = $(this);
      var entry = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
      if (!entry) return;
      var posStr = ($tr.find('.mapper-cluster-pos').val() || '').trim();
      var pos = parseFloat(posStr);
      if (isNaN(pos)) { allNumeric = false; return false; }
      data[entry] = { Value2: pos };
      if (plotMode === 'numbers') {
        var v1 = parseNumLoose($tr.find('.mapper-cluster-gradient').val());
        data[entry].Value1 = v1 != null ? v1 : 0;
      }
    });

    if (!allNumeric) {
      $('#mapper-cluster-messages').text('All position values must be numbers.');
      showPositionOverlayFilling();
      return;
    }
    var csrf = getCsrfToken();
    if (!csrf) {
      $('#mapper-cluster-messages').text('CSRF token missing — please refresh the page.');
      return;
    }

    currentPositionXhr = $.ajax({
      url: '/mapper/Cluster',
      type: 'POST',
      data: { Data: JSON.stringify(data), csrfmiddlewaretoken: csrf },
      headers: { 'X-Requested-With': 'XMLHttpRequest' },
      success: function (resp) {
        currentPositionXhr = null;
        var seq = resp.cluster_data_seq;
        if (!seq || seq.length < 2) {
          $('#mapper-cluster-messages').text('Not enough data returned. Check your inputs.');
          showPositionOverlayFilling();
          return;
        }
        var inputData = mapperClusterBuildData();
        var posMap = getClusterPosMap();
        var newData = [];
        seq.forEach(function (pt) {
          var key = pt.label + '_human';
          var row = inputData[key] || inputData[pt.label];
          var meta = posMap[String(pt.label || '').toLowerCase()] || {};
          var newPt = {
            x: pt.x, y: pt.y, label: pt.label, original_label: pt.label,
            'Class': meta['Class'] || '', 'Ligand type': meta['Ligand type'] || '',
            'Receptor family': meta['Receptor family'] || '',
            fill: row ? (row.Value1 != null ? row.Value1 : 0) : 0,
            cluster: 0
          };
          if (plotMode === 'categories') newPt.userCategory = row ? (row.userCategory || 'Unknown') : 'Unknown';
          newData.push(newPt);
        });
        if (newData.length < 2) {
          $('#mapper-cluster-messages').text('Not enough matched receptors.');
          showPositionOverlayFilling();
          return;
        }
        currentClusterData = newData;
        positionProcessed = true;
        _clusterZoomRange = null;
        // Cache geometry (x/y/label/metadata) keyed by current receptor+position fingerprint.
        // fill and userCategory are NOT stored — cosmetic redraw re-applies them on restore.
        positionResultCache[plotMode] = {
          points: currentClusterData.map(function (pt) {
            return {
              x: pt.x, y: pt.y, label: pt.label, original_label: pt.original_label,
              'Class': pt['Class'] || '', 'Ligand type': pt['Ligand type'] || '',
              'Receptor family': pt['Receptor family'] || '', cluster: pt.cluster || 0
            };
          }),
          fingerprint: mapperClusterGetPosFingerprint()
        };
        hidePositionOverlay();
        hidePlaceholder();
        $('#mapper-cluster-messages').text('');
        var maxClusters = Math.min(20, currentClusterData.length - 1);
        var slider = document.getElementById('mapper-cluster-slider');
        if (slider) {
          slider.max = maxClusters;
          if (parseInt(slider.value, 10) > maxClusters) slider.value = Math.min(5, maxClusters);
          document.getElementById('mapper-cluster-slider-val').textContent = slider.value;
        }
        var _posEl = document.getElementById('plotContainer_cluster');
        if (_posEl) { _posEl.style.width = ''; _posEl.style.height = ''; }
        var _posFP = currentClusterData.map(function(d){ return d.original_label; }).sort().join(',');
        var _posIsNew = _posFP !== _lastClusterEntries;
        _lastClusterEntries = _posFP;
        performClustering(Math.min(5, currentClusterData.length - 1), _posIsNew);
      },
      error: function (xhr, status) {
        currentPositionXhr = null;
        if (status === 'abort') return;
        $('#mapper-cluster-messages').text('Error communicating with server. Please try again.');
        showPositionOverlayFilling();
      }
    });
  }

  // ── Row management ─────────────────────────────────────────────────────────
  function mapperClusterCreateReceptorTd() {
    var $wrap  = $('<div class="mapper-core-receptor-input-wrap is-empty">');
    var $hidden = $('<input type="hidden" class="mapper-core-receptor-entry">');
    var $ta    = $('<textarea rows="1" class="mapper-core-in-receptor form-control input-sm">');
    var $view  = $('<div class="mapper-core-receptor-html-view" tabindex="0">');
    var $clear = $('<button type="button" class="mapper-core-receptor-clear" aria-label="Clear receptor">\xd7</button>');
    $wrap.append($hidden, $ta, $view, $clear);
    return $wrap;
  }

  function mapperClusterFindBlankRow() {
    var $blank = $();
    $('#mapper-cluster-input-tbody tr').each(function () {
      var entry = $(this).find('.mapper-core-receptor-entry').val() || '';
      var ta    = $(this).find('.mapper-core-in-receptor').val()    || '';
      if (!entry && !ta.trim()) { $blank = $(this); return false; }
    });
    return $blank;
  }

  function mapperClusterAppendRow(opts) {
    opts = opts || {};
    var $tbody = $('#mapper-cluster-input-tbody');
    if ($tbody.find('tr').length >= MAPPER_CLUSTER_MAX_ROWS) return null;

    var $tr = $('<tr>');

    // Remove cell
    var $removeCell = $('<td class="mapper-core-remove-cell">');
    var $removeBtn  = $('<button type="button" class="mapper-core-remove-row" aria-label="Remove row" title="Remove row">\xd7</button>');
    $removeCell.append($removeBtn);
    $tr.append($removeCell);

    // Receptor cell
    var $receptorCell = $('<td class="mapper-core-receptor-cell">');
    $receptorCell.append(mapperClusterCreateReceptorTd());
    $tr.append($receptorCell);

    // Gradient / category value cell
    var $valCell    = $('<td class="mapper-cluster-val-cell">');
    var $gradInput  = $('<input type="text" class="form-control input-sm mapper-cluster-gradient val-empty" maxlength="30">');
    var $catInput   = $('<input type="text" class="form-control input-sm mapper-cluster-cat val-empty" maxlength="50">').hide();
    if (plotMode === 'categories') { $gradInput.hide(); $catInput.show(); }
    $valCell.append($gradInput, $catInput);
    $tr.append($valCell);

    // Color swatch cell — visibility controlled by mapper-core-text-mode on the table
    var $colorCell  = $('<td class="mapper-cluster-swatch-cell">');
    var $swatchInp  = $('<input type="text" class="mapper-cluster-label-swatch mapper-core-row-color-picker">').hide();
    $colorCell.append($swatchInp);
    $tr.append($colorCell);

    // Position cell
    var $posCell  = $('<td class="mapper-cluster-pos-cell">');
    var $posInput = $('<input type="text" class="form-control input-sm mapper-cluster-pos pos-empty" maxlength="20">');
    $posCell.append($posInput);
    $tr.append($posCell);

    // Position input handler — validate on key release, update counter, debounce check.
    $posInput.on('keyup change', function () {
      var val = $(this).val();
      var hasV = _isValidPos(val);
      $(this).toggleClass('pos-invalid', val.trim() !== '' && !hasV);
      $(this).toggleClass('pos-empty',   val.trim() === '');
      var hasR = !!$tr.data('pos-r');
      _posRowUpdate($tr, hasR, hasV);
      // Only invalidate processed state when this row has a receptor.
      // Edits to rows without a receptor cannot affect the cluster calculation,
      // and resetting positionProcessed here was causing repeated AJAX cycles.
      if (hasR) {
        positionProcessed = false;
        positionResultCache[plotMode] = null;
      }
      clearTimeout(posCheckTimer);
      posCheckTimer = setTimeout(mapperClusterCheckPositionMode, 120);
    });

    $tbody.append($tr);

    // Remove handler
    $removeBtn.on('click', function () {
      _posRowRemove($tr);
      $tr.remove();
      positionResultCache[plotMode] = null;
      if (inPositionMode) positionProcessed = false;
      _clusterTrailingBlankRowDebounced.schedule();
      mapperClusterCheckPositionMode();
      if (!inPositionMode) mapperClusterScheduleRedraw();
    });

    // Bind autocomplete
    mapperClusterBindAc($tr.find('.mapper-core-in-receptor'));

    // Set values if provided
    if (opts.entry) {
      $tr.find('.mapper-core-receptor-entry').val(opts.entry);
      var meta = window.MAPPER_CORE_ENTRY_META && window.MAPPER_CORE_ENTRY_META[opts.entry];
      if (meta) {
        $tr.find('.mapper-core-in-receptor').hide();
        $tr.find('.mapper-core-receptor-html-view').html(mapperClusterResolvedDisplay(opts.entry)).show();
        $tr.find('.mapper-core-receptor-input-wrap').removeClass('is-empty');
      }
      $tr.addClass('has-receptor');
    }
    if (!opts.entry && opts.receptorText) {
      $tr.find('.mapper-core-in-receptor').val(opts.receptorText);
      if (opts.receptorText.trim()) $tr.find('.mapper-core-receptor-input-wrap').removeClass('is-empty');
    }
    if (opts.gradient != null) {
      var gv = String(opts.gradient);
      $tr.find('.mapper-cluster-gradient').val(gv).toggleClass('val-empty', gv.trim() === '');
    }
    if (opts.category != null) {
      var cv = String(opts.category);
      $tr.find('.mapper-cluster-cat').val(cv).toggleClass('val-empty', cv.trim() === '');
    }
    if (opts.pos != null && opts.pos !== '') {
      $posInput.val(String(opts.pos));
    }

    // Init row color swatch (categories mode)
    if (plotMode === 'categories') mapperClusterInitRowSwatch($tr);

    return $tr;
  }

  function mapperClusterInitRowSwatch($tr) {
    var $inp = $tr.find('.mapper-cluster-label-swatch');
    if (!$inp.length || !$.fn.spectrum) return;
    if (!$tr.hasClass('has-receptor')) {
      try { $inp.spectrum('destroy'); } catch(e) {}
      return;
    }
    try { $inp.spectrum('destroy'); } catch(e) {}
    var lbl = $tr.find('.mapper-cluster-cat').val().trim();
    var hex = mapperClusterGetLabelColor(lbl || '');
    $inp.spectrum({
      color: hex,
      showPalette: true, showInput: true, showButtons: false, preferredFormat: 'hex',
      appendTo: '.mapper-core-wheel-left',
      containerClassName: 'mapper-cluster-row-swatch-sp-container',
      replacerClassName: 'mapper-cluster-row-swatch-replacer',
      palette: CLUSTER_COLOR_PALETTE,
      change: function(c) {
        var newHex = c ? c.toHexString() : hex;
        var label = $tr.find('.mapper-cluster-cat').val().trim();
        if (label) {
          CLUSTER_LABEL_COLORS[label] = newHex;
          mapperClusterSyncRowSwatchesByLabel(label);
          if ($('#cluster-colors-btn').closest('.dropdown').hasClass('open')) mapperClusterRebuildColorsPanel();
          mapperClusterScheduleRedraw();
        }
      },
      move: function(c) {
        var newHex = c ? c.toHexString() : hex;
        var label = $tr.find('.mapper-cluster-cat').val().trim();
        if (label) {
          CLUSTER_LABEL_COLORS[label] = newHex;
          mapperClusterSyncRowSwatchesByLabel(label);
          mapperClusterScheduleRedraw();
        }
      }
    });
  }

  // Batches mapperClusterEnsureTrailingBlankRow's full-table scan so repeated row
  // deletions (the only hot path that hits it — value edits are already row-scoped)
  // coalesce into one pass per debounce window instead of one per click.
  var _clusterTrailingBlankRowDebounced = MapperPageCore.debounce(function () {
    mapperClusterEnsureTrailingBlankRow();
  }, 200);

  function mapperClusterEnsureTrailingBlankRow() {
    var $tbody = $('#mapper-cluster-input-tbody');
    var $rows  = $tbody.find('tr');
    var blankCount = 0;
    $rows.each(function () {
      var entry = $(this).find('.mapper-core-receptor-entry').val() || '';
      var ta    = $(this).find('.mapper-core-in-receptor').val()    || '';
      if (!entry && !ta.trim()) blankCount++;
    });
    if (blankCount < 1) mapperClusterAppendRow();
    if (blankCount > 2) {
      var removed = 0;
      $rows.get().reverse().forEach(function (tr) {
        if (removed >= blankCount - 1) return;
        var $tr   = $(tr);
        var entry = $tr.find('.mapper-core-receptor-entry').val() || '';
        var ta    = $tr.find('.mapper-core-in-receptor').val()    || '';
        if (!entry && !ta.trim()) { $tr.remove(); removed++; }
      });
    }
    mapperClusterSyncFirstRowPlaceholder();
  }

  // ── Position counter helpers ────────────────────────────────────────────────
  // Strict numeric validator — rejects "5.222a" or "abc5"; only pure float/int strings pass.
  function _isValidPos(val) {
    var s = (val || '').trim().replace(',', '.');
    return s !== '' && /^-?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/.test(s);
  }
  // Sync tbody class that enables orange-outline on empty position cells.
  function _posTableSync() {
    $('#mapper-cluster-input-tbody').toggleClass('pos-partially-filled', _posCtFilled > 0);
  }
  // Called whenever a row's receptor or position state changes.
  // Reads the previous state from jQuery .data() and updates the counters.
  function _posRowUpdate($tr, hasReceptor, hasValue) {
    var prevR = !!$tr.data('pos-r');
    var prevV = !!$tr.data('pos-v');
    $tr.data('pos-r', hasReceptor).data('pos-v', hasValue);
    if (prevR !== hasReceptor) _posCtReceptors += hasReceptor ? 1 : -1;
    var prevFilled = prevR && prevV, newFilled = hasReceptor && hasValue;
    if (prevFilled !== newFilled) _posCtFilled += newFilled ? 1 : -1;
    _posTableSync();
  }
  // Called when a row is removed from the DOM.
  function _posRowRemove($tr) {
    var hasR = !!$tr.data('pos-r'), hasV = !!$tr.data('pos-v');
    if (hasR) { _posCtReceptors--; if (hasV) _posCtFilled--; }
    _posTableSync();
  }
  // Full recount after bulk table rebuilds (snapshot restore, demo fill, clear).
  function _posCountersReset() {
    _posCtReceptors = 0; _posCtFilled = 0;
    $('#mapper-cluster-input-tbody tr').each(function () {
      var hasR = ($(this).find('.mapper-core-receptor-entry').val() || '').trim() !== '';
      var $pos = $(this).find('.mapper-cluster-pos');
      var posVal = $pos.val() || '';
      var hasV = _isValidPos(posVal);
      $pos.toggleClass('pos-invalid', posVal.trim() !== '' && !hasV);
      $pos.toggleClass('pos-empty',   posVal.trim() === '');
      $(this).data('pos-r', hasR).data('pos-v', hasV);
      if (hasR) { _posCtReceptors++; if (hasV) _posCtFilled++; }
    });
    _posTableSync();
  }

  function mapperClusterGetPosFingerprint() {
    var parts = [];
    $('#mapper-cluster-input-tbody tr').each(function () {
      var entry = ($(this).find('.mapper-core-receptor-entry').val() || '').trim();
      if (!entry) return;
      var pos = ($(this).find('.mapper-cluster-pos').val() || '').trim();
      parts.push(entry + ':' + pos);
    });
    return parts.join('|');
  }

  function mapperClusterTryRestorePositionCache() {
    var cache = positionResultCache[plotMode];
    if (!cache || !cache.points || !cache.points.length) return false;
    var fp = mapperClusterGetPosFingerprint();
    if (fp !== cache.fingerprint) {
      positionResultCache[plotMode] = null;
      return false;
    }
    // Shallow-copy the cached points so cosmetic redraw doesn't mutate the cache
    currentClusterData = cache.points.map(function (pt) {
      return Object.assign({}, pt);
    });
    inPositionMode = true;
    positionProcessed = true;
    hidePositionOverlay();
    hidePlaceholder();
    // Apply current mode's fill/userCategory onto the restored geometry
    var inputData = mapperClusterBuildData();
    currentClusterData.forEach(function (point) {
      var key = point.label + '_human';
      var row = inputData[key] || inputData[point.label];
      if (!row) return;
      point.fill = row.Value1 != null ? row.Value1 : 0;
      if (plotMode === 'categories') point.userCategory = row.userCategory || 'Unknown';
    });
    updateMultiClassWarning(inputData);
    return true;
  }

  function mapperClusterSyncFirstRowPlaceholder() {
    var hasAny = false;
    $('#mapper-cluster-input-tbody tr').each(function () {
      if ($(this).find('.mapper-core-receptor-entry').val()) { hasAny = true; return false; }
    });
    $('#mapper-cluster-input-table').toggleClass('mapper-core-receptors-compact', hasAny);
    $('#mapper-cluster-input-tbody tr').each(function (idx) {
      var $inp = $(this).find('.mapper-core-in-receptor');
      if (!$inp.length) return;
      if (idx === 0 && !hasAny) $inp.attr('placeholder', MAPPER_CLUSTER_PLACEHOLDER);
      else                      $inp.removeAttr('placeholder');
    });
  }

  // ── Autocomplete ────────────────────────────────────────────────────────────

  // Returns the label-type-appropriate seed text for a resolved receptor's edit field.
  function mapperClusterEditSeed(entryId) {
    var sid  = entryId != null ? String(entryId).trim() : '';
    var meta = (window.MAPPER_CORE_ENTRY_META && window.MAPPER_CORE_ENTRY_META[sid]) || {};
    var ltype = getActiveLabelType();
    if (ltype === 'Entrez'  && meta.gene)    return meta.gene;
    if (ltype === 'UniProt' && meta.uniprot) return meta.uniprot;
    return meta.name_plain || '';
  }

  function mapperClusterFilterLocal(term) {
    var t     = term.trim().toUpperCase();
    var meta  = window.MAPPER_CORE_ENTRY_META || {};
    var ltype = getActiveLabelType(); // 'IUPHAR' | 'Entrez' | 'UniProt'
    var results = [];
    Object.keys(meta).forEach(function (entry_name) {
      var m = meta[entry_name];
      // Always search across all name forms so the user can find by any schema
      var searchIn = [(m.name_plain || ''), (m.gene || ''), (m.uniprot || ''), entry_name].join(' ').toUpperCase();
      if (searchIn.indexOf(t) !== -1) {
        var dispHtml;
        if (ltype === 'Entrez' && m.gene) {
          dispHtml = $('<span/>').text(m.gene).html();
        } else if (ltype === 'UniProt' && m.uniprot) {
          dispHtml = $('<span/>').text(m.uniprot).html();
        } else {
          dispHtml = m.name_html || $('<span/>').text(m.name_plain || entry_name).html();
        }
        results.push({ label: m.name_plain || entry_name, html: dispHtml, value: entry_name, sortKey: m.name_plain || '' });
      }
    });
    results.sort(function (a, b) { return a.sortKey.localeCompare(b.sortKey); });
    return results.slice(0, 40);
  }

  function mapperClusterResolvedDisplay(entryId) {
    var sid  = entryId != null ? String(entryId).trim() : '';
    var meta = (window.MAPPER_CORE_ENTRY_META && window.MAPPER_CORE_ENTRY_META[sid]) || {};
    var ltype = getActiveLabelType(); // 'IUPHAR' | 'Entrez' | 'UniProt'
    if (ltype === 'Entrez' && meta.gene)    return $('<span/>').text(meta.gene).html();
    if (ltype === 'UniProt' && meta.uniprot) return $('<span/>').text(meta.uniprot).html();
    return meta.name_html ? String(meta.name_html) : $('<span/>').text(sid || '').html();
  }

  function mapperClusterSetResolved($tr, entry_name) {
    var meta = window.MAPPER_CORE_ENTRY_META && window.MAPPER_CORE_ENTRY_META[entry_name];
    if (!meta) return;
    $tr.find('.mapper-core-receptor-entry').val(entry_name);
    $tr.find('.mapper-core-in-receptor').val('').hide();
    $tr.find('.mapper-core-receptor-html-view').html(mapperClusterResolvedDisplay(entry_name)).show();
    $tr.find('.mapper-core-receptor-input-wrap').removeClass('is-empty');
    $tr.addClass('has-receptor').removeClass('receptor-unknown');
    _posRowUpdate($tr, true, _isValidPos($tr.find('.mapper-cluster-pos').val()));
    positionResultCache[plotMode] = null;
    if (inPositionMode) positionProcessed = false;
    if (plotMode === 'categories') mapperClusterInitRowSwatch($tr);
    mapperClusterEnsureTrailingBlankRow();
    mapperClusterCheckPositionMode();
    if (!inPositionMode) mapperClusterScheduleRedraw();
  }

  function mapperClusterBindAc($ta) {
    if (!$.ui || !$.ui.autocomplete) return;
    $ta.autocomplete({
      minLength: 1,
      classes: { 'ui-autocomplete': 'mapper-core-receptor-ac-menu' },
      source: function (req, resp) {
        resp(mapperClusterFilterLocal(req.term).map(function (r) {
          return { label: r.label, html: r.html, value: r.value };
        }));
      },
      select: function (event, ui) {
        event.preventDefault();
        var $tr = $ta.closest('tr');
        mapperClusterSetResolved($tr, ui.item.value);
      }
    });
    $ta.data('ui-autocomplete') && ($ta.data('ui-autocomplete')._renderItem = function (ul, item) {
      var inner = item.html || $('<span/>').text(item.label || '').html();
      return $('<li>').append($('<div class="mapper-core-ac-item-label">').html(inner)).appendTo(ul);
    });
  }

  // ── Demo fill ──────────────────────────────────────────────────────────────
  function mapperClusterFillDemo() {
    suppressRedraw = true;

    // Set default palette colors for the demo categories
    CLUSTER_LABEL_COLORS['Group A'] = '#3c5488';
    CLUSTER_LABEL_COLORS['Group B'] = '#e64b35';
    CLUSTER_LABEL_COLORS['Group C'] = '#00a087';

    // Resolve all demo entries first
    var demoResolved = MAPPER_CLUSTER_DEMO_ROWS.map(function (row) {
      var entry = window.MAPPER_CORE_RESOLVE && window.MAPPER_CORE_RESOLVE[row.receptor.toUpperCase()];
      if (!entry) return null;
      return { entry: entry, gradient: String(row.gradient || ''), category: String(row.category || '') };
    }).filter(Boolean);

    // Populate only the current mode's snapshot — Numbers and Categories stay fully separate
    MAPPER_CLUSTER_MODE_SNAPSHOTS[plotMode] = {
      rows: demoResolved.map(function (r) {
        return plotMode === 'numbers'
          ? { entry: r.entry, receptorText: '', gradient: r.gradient, pos: '' }
          : { entry: r.entry, receptorText: '', category: r.category, pos: '' };
      })
    };

    // Rebuild DOM for the currently active mode
    $('#mapper-cluster-input-tbody tr').each(function () {
      var $ta = $(this).find('.mapper-core-in-receptor');
      try { if ($ta.data('ui-autocomplete')) $ta.autocomplete('destroy'); } catch (e) {}
    });
    $('#mapper-cluster-input-tbody').empty();
    var curSnap = MAPPER_CLUSTER_MODE_SNAPSHOTS[plotMode];
    curSnap.rows.forEach(function (r) {
      var rowOpts = { entry: r.entry };
      if (plotMode === 'numbers') rowOpts.gradient  = r.gradient;
      else                        rowOpts.category = r.category;
      mapperClusterAppendRow(rowOpts);
    });

    mapperClusterEnsureTrailingBlankRow();
    _posCountersReset();
    suppressRedraw = false;
    mapperClusterCheckPositionMode();
    if (!inPositionMode) mapperClusterScheduleRedraw();
  }

  // ── Download ───────────────────────────────────────────────────────────────
  function downloadClusterPlot(format) {
    var plotElement = document.getElementById('plotContainer_cluster');
    if (!plotElement || !currentClusterData.length) return;
    var actualFormat = (format === 'tiff') ? 'png' : (format === 'jpg') ? 'jpeg' : format;
    Plotly.toImage(plotElement, { format: actualFormat, height: 700, width: 1024, scale: 3 })
      .then(function (dataUrl) {
        var link = document.createElement('a');
        link.href = dataUrl;
        link.download = 'Cluster.' + format;
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
      });
  }
  window.downloadClusterPlot = downloadClusterPlot;

  // ── Boot ───────────────────────────────────────────────────────────────────
  function mapperClusterBoot() {
    // Remove loading overlay — page is now interactive
    $('.mapper-cluster-booting').removeClass('mapper-cluster-booting');
    mapperClusterBuildLabelConverter();

    // Init multi-class warning tooltip
    $('#cluster-multiclass-warning').tooltip({
      template: '<div class="tooltip cluster-warning-tooltip" role="tooltip"><div class="tooltip-arrow"></div><div class="tooltip-inner"></div></div>',
      placement: 'bottom',
      container: 'body'
    });

    // Initial blank row
    mapperClusterAppendRow();
    mapperClusterSyncFirstRowPlaceholder();

    showPlaceholder();

    // Numbers / Categories toggle
    $('#mapper-cluster-numbers').on('click',    function () { setPlotMode('numbers'); });
    $('#mapper-cluster-categories').on('click', function () { setPlotMode('categories'); });

    // Demo
    $('#mapper-cluster-demo-rows').on('click', mapperClusterFillDemo);

    // Clear
    function mapperClusterDestroyAcWidgets() {
      $('#mapper-cluster-input-tbody tr').each(function () {
        var $ta = $(this).find('.mapper-core-in-receptor');
        try { if ($ta.data('ui-autocomplete')) $ta.autocomplete('destroy'); } catch (e) {}
      });
    }
    function mapperClusterResetPositionState() {
      clearTimeout(processTimer); processTimer = null;
      if (currentPositionXhr) { currentPositionXhr.abort(); currentPositionXhr = null; }
      inPositionMode = false;
      positionProcessed = false;
      positionResultCache[plotMode] = null;
      hidePositionOverlay();
    }
    function mapperClusterDoWholeTableClear() {
      suppressRedraw = true;
      mapperClusterDestroyAcWidgets();
      currentClusterData = [];
      $('#mapper-cluster-input-tbody').empty();
      MAPPER_CLUSTER_MODE_SNAPSHOTS.numbers    = null;
      MAPPER_CLUSTER_MODE_SNAPSHOTS.categories = null;
      mapperClusterResetPositionState();
      _posCtReceptors = 0; _posCtFilled = 0;
      $('#cluster-multiclass-warning').hide();
      showPlaceholder();
      mapperClusterAppendRow();
      mapperClusterSyncFirstRowPlaceholder();
      suppressRedraw = false;
      $('#mapper-cluster-clear-rows').addClass('mapper-core-clear-clean').blur();
    }
    function mapperClusterDoCurrentModeClear() {
      var curKey = plotMode === 'numbers' ? 'numbers' : 'categories';
      suppressRedraw = true;
      mapperClusterDestroyAcWidgets();
      currentClusterData = [];
      $('#mapper-cluster-input-tbody').empty();
      MAPPER_CLUSTER_MODE_SNAPSHOTS[curKey] = null;
      mapperClusterResetPositionState();
      $('#cluster-multiclass-warning').hide();
      showPlaceholder();
      mapperClusterAppendRow();
      mapperClusterSyncFirstRowPlaceholder();
      suppressRedraw = false;
      $('#mapper-cluster-clear-rows').addClass('mapper-core-clear-clean').blur();
    }
    $('#mapper-cluster-clear-whole').on('click', function() {
      mapperClusterDoWholeTableClear();
    });
    $('#mapper-cluster-clear-mode').on('click', function() {
      mapperClusterDoCurrentModeClear();
    });
    $('#mapper-cluster-clear-gradient').on('click', function() {
      $('#mapper-cluster-input-tbody tr').each(function() {
        $(this).find('.mapper-cluster-gradient').val('').trigger('input');
      });
      mapperClusterScheduleRedraw();
      $('#mapper-cluster-clear-rows').removeClass('mapper-core-clear-clean');
    });
    $('#mapper-cluster-clear-cat').on('click', function() {
      $('#mapper-cluster-input-tbody tr').each(function() {
        $(this).find('.mapper-cluster-cat').val('').trigger('input');
      });
      mapperClusterScheduleRedraw();
      $('#mapper-cluster-clear-rows').removeClass('mapper-core-clear-clean');
    });
    $('#mapper-cluster-clear-pos').on('click', function() {
      suppressRedraw = true;
      $('#mapper-cluster-input-tbody tr').each(function() {
        $(this).find('.mapper-cluster-pos').val('').removeClass('pos-invalid');
      });
      suppressRedraw = false;
      mapperClusterResetPositionState();
      _posCountersReset();
      mapperClusterScheduleRedraw();
    });

    // Cancel position calculation — abort XHR, clear position column, return to live mode
    $('#cluster-pos-cancel-btn').on('click', function() {
      suppressRedraw = true;
      $('#mapper-cluster-input-tbody tr').each(function() {
        $(this).find('.mapper-cluster-pos').val('').removeClass('pos-invalid');
      });
      suppressRedraw = false;
      mapperClusterResetPositionState();
      _posCountersReset();
      mapperClusterScheduleRedraw();
    });

    // Plot type toggle (Dot ↔ Text)
    var plotTypeBtn = document.getElementById('mapper-cluster-plottype-btn');
    if (plotTypeBtn) {
      plotTypeBtn.addEventListener('click', function () {
        labelsVisible = !labelsVisible;
        this.textContent = labelsVisible ? 'Plot type: Text' : 'Plot type: Dot';
        $('.mapper-cluster-text-only').toggle(labelsVisible);
        $('.mapper-cluster-dot-only').toggle(!labelsVisible);
        if (currentClusterData.length) updatePlotWithAnnotations();
      });
    }

    // Text color toggle
    var textColorBtn = document.getElementById('mapper-cluster-textcolor-btn');
    if (textColorBtn) {
      textColorBtn.addEventListener('click', function () {
        cluster_DataStyling.textColorEnabled = !cluster_DataStyling.textColorEnabled;
        this.textContent  = cluster_DataStyling.textColorEnabled ? 'Shown' : 'Hidden';
        this.className    = 'btn btn-xs ' + (cluster_DataStyling.textColorEnabled ? 'btn-success' : 'btn-danger');
        if (currentClusterData.length) updatePlotWithAnnotations();
      });
    }

    // Label font size
    $('#mapper-cluster-label-fontsize').on('input', function () {
      var v = parseInt($(this).val(), 10);
      $('#mapper-cluster-label-fontsize-val').text(v);
      cluster_DataStyling.labelFontSize = v;
      if (currentClusterData.length) updatePlotWithAnnotations();
    });

    // Marker size
    $('#mapper-cluster-marker-size').on('input', function () {
      var v = parseInt($(this).val(), 10);
      $('#mapper-cluster-marker-size-val').text(v);
      cluster_DataStyling.markerSize = v;
      if (currentClusterData.length) updatePlotWithAnnotations();
    });

    // Stroke width
    $('#mapper-cluster-stroke-width').on('input', function () {
      var v = parseFloat($(this).val());
      $('#mapper-cluster-stroke-width-val').text(v);
      cluster_DataStyling.strokeWidth = v;
      if (currentClusterData.length) updatePlotWithAnnotations();
    });

    // Clustering slider — always switches coloring to cluster
    $('#mapper-cluster-slider').on('input', function () {
      var v = parseInt($(this).val(), 10);
      $('#mapper-cluster-slider-val').text(v);
      if (currentClusterData.length) {
        setActiveColorButton('cluster-colorByCluster');
        performClustering(v);
      }
    });

    // Color buttons
    document.querySelectorAll('#mapper-cluster-colors-menu .btn').forEach(function (btn) {
      btn.addEventListener('click', function () {
        setActiveColorButton(this.id);
        if (currentClusterData.length) updatePlotWithAnnotations();
      });
    });

    // Receptor name buttons
    document.querySelectorAll('.cluster-leaf-btn').forEach(function (btn) {
      btn.addEventListener('click', function () {
        document.querySelectorAll('.cluster-leaf-btn').forEach(function (b) {
          b.classList.remove('btn-primary');
          b.classList.add('btn-outline-primary');
        });
        this.classList.add('btn-primary');
        this.classList.remove('btn-outline-primary');
        // Refresh input table chips to show the new label name
        $('#mapper-cluster-input-tbody tr').each(function () {
          var sid = ($(this).find('.mapper-core-receptor-entry').val() || '').trim();
          if (!sid) return;
          var $view = $(this).find('.mapper-core-receptor-html-view');
          if ($view.length && $view.is(':visible')) $view.html(mapperClusterResolvedDisplay(sid));
        });
        if (currentClusterData.length) switchClusterLabels(this.getAttribute('data-value'));
      });
    });

    // Value input → redraw + track empty/invalid state for orange/red border
    $('#mapper-cluster-input-tbody').on('input change', '.mapper-cluster-gradient, .mapper-cluster-cat', function () {
      var val = $(this).val().trim();
      $(this).toggleClass('val-empty', val === '');
      if ($(this).hasClass('mapper-cluster-gradient')) {
        $(this).toggleClass('val-invalid', val !== '' && !_isValidPos(val));
        setActiveColorButton('cluster-colorByGradient');
      } else {
        // Category changed — update this row's swatch to the color for the new label
        var $tr = $(this).closest('tr');
        try { $tr.find('.mapper-cluster-label-swatch').spectrum('set', mapperClusterGetLabelColor(val)); } catch(e) {}
        setActiveColorButton('cluster-colorByUserCategory');
      }
      mapperClusterScheduleRedraw();
    });

    // Tab within each column: jump to the next empty cell in the SAME column.
    // Receptor column: → next empty receptor.
    // Gradient/category column: → next empty gradient/cat.
    // Position column: → next empty position.
    $('#mapper-cluster-input-tbody')
      .off('keydown.clusterTab')
      .on('keydown.clusterTab',
          '.mapper-core-in-receptor, .mapper-cluster-gradient, .mapper-cluster-cat, .mapper-cluster-pos',
          function (e) {
        if (e.which !== 9 || e.shiftKey || e.altKey || e.ctrlKey || e.metaKey) return;
        // Let AC dropdown handle Tab when an item is highlighted
        if ($(this).hasClass('mapper-core-in-receptor') && $(this).data('ui-autocomplete')) {
          var $menu = $(this).autocomplete('widget');
          if ($menu && $menu.is(':visible') && $menu.find('.ui-state-focus, .ui-state-active').length) return;
        }
        var $el = $(this);
        var $all, idx, $target;
        if ($el.hasClass('mapper-core-in-receptor')) {
          $all = $('#mapper-cluster-input-tbody .mapper-core-in-receptor:visible');
          idx  = $all.index(this);
          for (var i = idx + 1; i < $all.length; i++) {
            var $row = $all.eq(i).closest('tr');
            if (!$row.find('.mapper-core-receptor-entry').val().trim() && !$all.eq(i).val().trim()) {
              $target = $all.eq(i); break;
            }
          }
        } else if ($el.hasClass('mapper-cluster-gradient') || $el.hasClass('mapper-cluster-cat')) {
          var sel = $el.hasClass('mapper-cluster-gradient')
            ? '#mapper-cluster-input-tbody tr .mapper-cluster-gradient:visible'
            : '#mapper-cluster-input-tbody tr .mapper-cluster-cat:visible';
          $all = $(sel); idx = $all.index(this);
          for (var j = idx + 1; j < $all.length; j++) {
            if ($all.eq(j).val().trim() === '') { $target = $all.eq(j); break; }
          }
        } else { // position column
          $all = $('#mapper-cluster-input-tbody .mapper-cluster-pos');
          idx  = $all.index(this);
          for (var k = idx + 1; k < $all.length; k++) {
            if ($all.eq(k).val().trim() === '') { $target = $all.eq(k); break; }
          }
        }
        if ($target) { e.preventDefault(); $target.focus(); }
      });

    // Receptor textarea events
    $('#mapper-cluster-input-tbody').on('input', '.mapper-core-in-receptor', function () {
      var $ta  = $(this);
      var val  = $ta.val().trim();
      var $wrap = $ta.closest('.mapper-core-receptor-input-wrap');
      var $row = $ta.closest('tr');
      $row.removeClass('receptor-unknown'); // clear on any edit
      if (!val) {
        $row.find('.mapper-core-receptor-entry').val('');
        $row.removeClass('has-receptor');
        $wrap.addClass('is-empty');
        _posRowUpdate($row, false, _isValidPos($row.find('.mapper-cluster-pos').val()));
        positionResultCache[plotMode] = null;
        if (inPositionMode) positionProcessed = false;
        mapperClusterEnsureTrailingBlankRow();
        mapperClusterCheckPositionMode();
        if (!inPositionMode) mapperClusterScheduleRedraw();
        return;
      }
      $wrap.removeClass('is-empty'); // has text: hide the empty-state styling
      var resolved = window.MAPPER_CORE_RESOLVE && window.MAPPER_CORE_RESOLVE[val.toUpperCase()];
      if (resolved) mapperClusterSetResolved($ta.closest('tr'), resolved);
    });

    // Mark unresolved receptor red on blur; clear on focus.
    // On blur, try a silent re-resolve first — this correctly handles the
    // chip re-edit case where the entry field was cleared but the textarea
    // still holds a valid IUPHAR/gene name that we can look up.
    $('#mapper-cluster-input-tbody')
      .on('blur',  '.mapper-core-in-receptor', function () {
        var $ta   = $(this);
        var val   = $ta.val().trim();
        var $tr   = $ta.closest('tr');
        var entry = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
        if (!val || entry) { $tr.removeClass('receptor-unknown'); return; }
        // Try MAPPER_CORE_RESOLVE (covers entry names, gene names, IUPHAR names)
        var resolved = window.MAPPER_CORE_RESOLVE && window.MAPPER_CORE_RESOLVE[val.toUpperCase()];
        if (resolved) { mapperClusterSetResolved($tr, resolved); return; }
        // UniProt accessions are not in MAPPER_CORE_RESOLVE — scan ENTRY_META directly
        if (window.MAPPER_CORE_ENTRY_META) {
          var vUp = val.toUpperCase();
          var found = Object.keys(window.MAPPER_CORE_ENTRY_META).filter(function (k) {
            var u = (window.MAPPER_CORE_ENTRY_META[k].uniprot || '').toUpperCase();
            return u && u === vUp;
          })[0];
          if (found) { mapperClusterSetResolved($tr, found); return; }
        }
        $tr.addClass('receptor-unknown');
      })
      .on('focus', '.mapper-core-in-receptor', function () {
        $(this).closest('tr').removeClass('receptor-unknown');
      });

    // Clear receptor button
    $('#mapper-cluster-input-tbody').on('click', '.mapper-core-receptor-clear', function () {
      var $tr = $(this).closest('tr');
      var $wrap = $tr.find('.mapper-core-receptor-input-wrap');
      $wrap.find('.mapper-core-receptor-entry').val('');
      $wrap.find('.mapper-core-receptor-html-view').hide().empty();
      $wrap.find('.mapper-core-in-receptor').val('').show();
      $wrap.addClass('is-empty');
      $tr.removeClass('has-receptor receptor-unknown');
      _posRowUpdate($tr, false, _isValidPos($tr.find('.mapper-cluster-pos').val()));
      positionResultCache[plotMode] = null;
      if (inPositionMode) positionProcessed = false;
      mapperClusterEnsureTrailingBlankRow();
      mapperClusterCheckPositionMode();
      if (!inPositionMode) mapperClusterScheduleRedraw();
    });

    // Click resolved chip to re-enter edit mode
    $(document)
      .off('click.mapperClusterHtmlEdit', '#mapper-cluster-input-tbody .mapper-core-receptor-html-view')
      .on('click.mapperClusterHtmlEdit', '#mapper-cluster-input-tbody .mapper-core-receptor-html-view', function () {
        var $tr   = $(this).closest('tr');
        $tr.removeClass('has-receptor');
        var $oldTa = $tr.find('.mapper-core-in-receptor');
        try { if ($oldTa.hasClass('ui-autocomplete-input')) $oldTa.autocomplete('destroy'); } catch (e) {}
        var hid   = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
        $(this).hide().empty();
        var $inp2 = $('<textarea class="form-control input-sm mapper-core-in-receptor" rows="1" autocomplete="off" spellcheck="false"></textarea>');
        var seedText = mapperClusterEditSeed(hid); // label-type-aware seed (Gene/UniProt/IUPHAR)
        $inp2.val(seedText);
        $tr.find('.mapper-core-receptor-input-wrap .mapper-core-in-receptor').remove();
        var $wrapEdit = $tr.find('.mapper-core-receptor-input-wrap').prepend($inp2);
        if (!seedText) $wrapEdit.addClass('is-empty'); else $wrapEdit.removeClass('is-empty');
        $tr.find('.mapper-core-receptor-entry').val('');
        _posRowUpdate($tr, false, _isValidPos($tr.find('.mapper-cluster-pos').val()));
        positionResultCache[plotMode] = null;
        if (inPositionMode) positionProcessed = false;
        mapperClusterBindAc($inp2);
        $inp2.show().focus();
        // Immediately open the autocomplete dropdown so the user sees matching results
        window.setTimeout(function () {
          var t = ($inp2.val() || '').trim();
          if (t.length >= 1 && $inp2.data('ui-autocomplete')) $inp2.autocomplete('search', t);
        }, 0);
        mapperClusterEnsureTrailingBlankRow();
        mapperClusterSyncFirstRowPlaceholder();
        mapperClusterCheckPositionMode();
        if (!inPositionMode) mapperClusterScheduleRedraw();
      });

    // Paste handler (tab-separated receptor + value)
    $('#mapper-cluster-input-tbody').on('paste', '.mapper-core-in-receptor', function (e) {
      var pastedData = (e.originalEvent || e).clipboardData.getData('text');
      if (!pastedData || pastedData.indexOf('\n') === -1) return;
      e.preventDefault();
      suppressRedraw = true;
      var lines = pastedData.replace(/\r\n/g, '\n').split('\n');
      lines.forEach(function (line) {
        if (!line.trim()) return;
        var parts    = line.split('\t');
        var receptor = parts[0].trim();
        var gradient = parts[1] ? parts[1].trim() : '';
        if (!receptor) return;
        var resolved = window.MAPPER_CORE_RESOLVE && window.MAPPER_CORE_RESOLVE[receptor.toUpperCase()];
        if (resolved) {
          mapperClusterAppendRow({ entry: resolved, gradient: gradient });
        } else {
          var $row = mapperClusterFindBlankRow();
          if (!$row.length) { mapperClusterAppendRow(); $row = $('#mapper-cluster-input-tbody tr').last(); }
          $row.find('.mapper-core-in-receptor').val(receptor);
          $row.find('.mapper-cluster-gradient').val(gradient);
        }
      });
      suppressRedraw = false;
      mapperClusterEnsureTrailingBlankRow();
      _posCountersReset();
      mapperClusterCheckPositionMode();
      if (!inPositionMode) mapperClusterScheduleRedraw();
    });

    // Zoom / relayout handler — only act on user-initiated axis changes.
    // Programmatic relayouts (hovermode, margin.r, etc.) must NOT trigger a re-render
    // or they create a feedback loop: relayout → plotly_relayout → render → relayout → …
    $('#plotContainer_cluster').on('plotly_relayout', function (evt, eventdata) {
      if (!currentClusterData.length) return;
      if (!eventdata) return;
      if (eventdata['xaxis.autorange'] === true) {
        _clusterZoomRange = null;
      } else if (eventdata['xaxis.range[0]'] !== undefined) {
        _clusterZoomRange = [
          [eventdata['xaxis.range[0]'], eventdata['xaxis.range[1]']],
          [eventdata['yaxis.range[0]'], eventdata['yaxis.range[1]']]
        ];
      } else {
        return; // programmatic relayout (hovermode, margin.r, etc.) — skip
      }
      switchClusterLabels(getActiveLabelType());
    });

    // Receptor picker modal
    if (typeof window.mapperCoreInitGpcromePickerModal === 'function') {
      window.mapperCoreInitGpcromePickerModal({
        pickerRows: $.isArray(window.MAPPER_CORE_GPCROME_PICKER_ROWS) ? window.MAPPER_CORE_GPCROME_PICKER_ROWS : [],
        maxRows:    MAPPER_CLUSTER_MAX_ROWS,
        onAdd: function (entryIds, meta) {
          suppressRedraw = true;
          var isNumberMode = (plotMode === 'numbers');
          var numberMap = (meta && meta.groupNameById && isNumberMode)
            ? window.mapperCoreBuildSequentialNumberMap(meta.groupNameById)
            : null;
          (entryIds || []).forEach(function (id) {
            var $row = mapperClusterFindBlankRow();
            if (!$row.length) { mapperClusterAppendRow(); $row = $('#mapper-cluster-input-tbody tr').last(); }
            mapperClusterSetResolved($row, id);
            if (meta && meta.groupNameById && meta.groupNameById[id] != null) {
              var $target = isNumberMode ? $row.find('.mapper-cluster-gradient') : $row.find('.mapper-cluster-cat');
              var assignedValue = isNumberMode ? numberMap[id] : meta.groupNameById[id];
              if (assignedValue != null) {
                $target.val(String(assignedValue)).trigger('input');
              }
            }
          });
          suppressRedraw = false;
          mapperClusterEnsureTrailingBlankRow();
          mapperClusterScheduleRedraw();
        }
      });
    }

    mapperClusterSyncClearDropdown();
    mapperClusterCheckPositionMode();

    // ── Colors panel init ─────────────────────────────────────────────────
    if ($.fn.spectrum) {
      // Gradient pickers — initialised once, visible in panel when gradient is active
      var _gradSpec = function(id, key) {
        $('#' + id).spectrum({
          color: CLUSTER_GRADIENT_COLORS[key],
          showPalette: true, showInput: true, showButtons: false, preferredFormat: 'hex',
          appendTo: '#cluster-colors-panel',
          containerClassName: 'mapper-core-lcat-sp-container',
          replacerClassName: 'cluster-grad-replacer',
          palette: CLUSTER_COLOR_PALETTE,
          change: function(c) { CLUSTER_GRADIENT_COLORS[key] = c ? c.toHexString() : CLUSTER_GRADIENT_COLORS[key]; mapperClusterScheduleRedraw(); },
          move:   function(c) { CLUSTER_GRADIENT_COLORS[key] = c ? c.toHexString() : CLUSTER_GRADIENT_COLORS[key]; mapperClusterScheduleRedraw(); }
        });
      };
      _gradSpec('cluster-grad-min-color', 'min');
      _gradSpec('cluster-grad-mid-color', 'mid');
      _gradSpec('cluster-grad-max-color', 'max');
    }

    // Number of colors select — show/hide min/mid pickers and update stop count
    function _syncGradStops(stops) {
      CLUSTER_GRADIENT_STOPS = stops;
      var showMin = stops >= 2;
      var showMid = stops >= 3;
      $('#cluster-grad-min-lbl').css('visibility', showMin ? 'visible' : 'hidden');
      $('#cluster-grad-mid-lbl').css('visibility', showMid ? 'visible' : 'hidden');
      $('#cluster-grad-min-wrap').css('visibility', showMin ? 'visible' : 'hidden');
      $('#cluster-grad-mid-wrap').css('visibility', showMid ? 'visible' : 'hidden');
    }
    $('#cluster-grad-num-colors').on('change', function() {
      _syncGradStops(parseInt($(this).val(), 10));
      mapperClusterScheduleRedraw();
    });
    _syncGradStops(3); // initial state — all three visible

    // Bootstrap dropdown opens Colors panel — rebuild on show
    $('#cluster-colors-btn').closest('.dropdown').on('shown.bs.dropdown', function() {
      mapperClusterRebuildColorsPanel();
    });

    // Hook called by datamapper.js after each plot update
    window.mapperClusterOnPlotUpdated = function() {
      if (!$('#cluster-colors-btn').closest('.dropdown').hasClass('open')) return;
      // Skip rebuild if the user has a label Spectrum picker open (avoids closing it mid-pick)
      if ($('.mapper-core-lcat-sp-container.sp-container:visible').length) return;
      mapperClusterRebuildColorsPanel();
    };
  }

  $(document).ready(function () { mapperClusterBoot(); });

})(jQuery);

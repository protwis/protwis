/**
 * Mapper 2 wheel: Numeric values vs categorical text Labels (Colours per unique label string).
 */
var mapperWheelWheelInputMode = 'numeric';
var MAPPER_CORE_LABEL_COLORS = {};
var MAPPER_CORE_LABEL_ENABLED = {};
/** Captured snapshots for round-trips keyed by wheel input mode. */
var MAPPER_WHEEL_MODE_SNAPSHOTS = {
  numeric: null,
  text: null
};
var MAPPER_WHEEL_SORT_STATE = {
  column: null,
  dir: 'asc'
};
var MAPPER_WHEEL_DEMO_ROWS = {
  '5HT1A': { numeric: 1, text: 'Aminergic' },
  '5HT1B': { numeric: 2, text: 'Aminergic' },
  '5HT1D': { numeric: 3, text: 'Aminergic' },
  '5HT1E': { numeric: 4, text: 'Aminergic' },
  '5HT1F': { numeric: 5, text: 'Aminergic' },
  '5HT2A': { numeric: 6, text: 'Aminergic' },
  '5HT2B': { numeric: 7, text: 'Aminergic' },
  '5HT2C': { numeric: 8, text: 'Aminergic' },
  'ACKR1': { numeric: 9, text: 'Chemokine' },
  'ACKR2': { numeric: 10, text: 'Chemokine' },
  'ACKR3': { numeric: 11, text: 'Chemokine' },
  'ACKR4': { numeric: 12, text: 'Chemokine' },
  'ACM1': { numeric: 13, text: 'Cholinergic' },
  'ACM2': { numeric: 14, text: 'Cholinergic' },
  'ACM3': { numeric: 15, text: 'Cholinergic' },
  'ADA1A': { numeric: 16, text: 'Adrenergic' },
  'ADA1B': { numeric: 17, text: 'Adrenergic' },
  'ADA1D': { numeric: 18, text: 'Adrenergic' },
  'ADRB1': { numeric: 19, text: 'Adrenergic' },
  'ADRB2': { numeric: 20, text: 'Adrenergic' }
};

function mapperWheelIsTextMode() {
  return mapperWheelWheelInputMode === 'text';
}

function mapperWheelMessageReceptorRequiredForLabelCell() {
  return 'Enter a receptor (use search) before typing a label in this row.';
}

function mapperWheelMessageLabelMissingForReceptor() {
  return 'Enter a label for this receptor to show it with a colour on the plot.';
}

// FNV1a32 hash + default label colour moved to MapperPageCore (mapper_page_core.js).
function mapperWheelDefaultHexForLabelKey(lbl) {
  return MapperPageCore.defaultColorForLabel(lbl);
}

function mapperWheelSerializeDomRows() {
  var rows = [];
  $('#mapper-wheel-input-tbody tr').each(function() {
    var $tr = $(this);
    var entry = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
    var unmatched = $tr.data('mapperCoreUnmatchedRaw');
    var unmatchedStr = unmatched != null ? String(unmatched) : '';
    var ta = mapperWheelNormalizeTypedReceptorInput($tr.find('.mapper-core-in-receptor').val() || '');
    var val = ($tr.find('.mapper-wheel-in-value').val() || '');
    rows.push({ entry: entry, unmatched: unmatchedStr, ta: ta, value: val });
  });
  return rows;
}

function mapperWheelNaturalSortCompare(a, b) {
  var collator = mapperWheelNaturalSortCompare._collator;
  if (!collator) {
    collator = new Intl.Collator(undefined, { numeric: true, sensitivity: 'base' });
    mapperWheelNaturalSortCompare._collator = collator;
  }
  return collator.compare(mapperWheelNormalizeSortText(a), mapperWheelNormalizeSortText(b));
}

function mapperWheelNormalizeSortText(value) {
  var s = value == null ? '' : String(value);
  if (/[&<]/.test(s)) {
    var txt = document.createElement('textarea');
    txt.innerHTML = s.replace(/<[^>]*>/g, '');
    s = txt.value || s;
  }
  var greekMap = {
    'α': 'a', 'Α': 'a',
    'β': 'b', 'Β': 'b',
    'γ': 'g', 'Γ': 'g',
    'δ': 'd', 'Δ': 'd',
    'ε': 'e', 'Ε': 'e',
    'ζ': 'z', 'Ζ': 'z',
    'η': 'h', 'Η': 'h',
    'θ': 't', 'Θ': 't',
    'ι': 'i', 'Ι': 'i',
    'κ': 'k', 'Κ': 'k',
    'λ': 'l', 'Λ': 'l',
    'μ': 'm', 'Μ': 'm',
    'ν': 'n', 'Ν': 'n',
    'ξ': 'x', 'Ξ': 'x',
    'ο': 'o', 'Ο': 'o',
    'π': 'p', 'Π': 'p',
    'ρ': 'r', 'Ρ': 'r',
    'σ': 's', 'ς': 's', 'Σ': 's',
    'τ': 't', 'Τ': 't',
    'υ': 'u', 'Υ': 'u',
    'φ': 'f', 'Φ': 'f',
    'χ': 'c', 'Χ': 'c',
    'ψ': 'p', 'Ψ': 'p',
    'ω': 'o', 'Ω': 'o'
  };
  return s.replace(/[αΑβΒγΓδΔεΕζΖηΗθΘιΙκΚλΛμΜνΝξΞοΟπΠρΡσςΣτΤυΥφΦχΧψΨωΩ]/g, function(ch) {
    return greekMap[ch] || ch;
  });
}

function mapperWheelSortHeaderLabel(column) {
  if (column === 'value') {
    return '';
  }
  return 'Receptor';
}

function mapperWheelUpdateSortHeaders() {
  $('#mapper-core-input-table th.mapper-core-sortable-head').each(function() {
    var $th = $(this);
    var col = $th.attr('data-mapper-wheel-sort-col');
    var active = MAPPER_WHEEL_SORT_STATE.column === col;
    var dir = active ? MAPPER_WHEEL_SORT_STATE.dir : null;
    $th.attr('aria-sort', active ? (dir === 'desc' ? 'descending' : 'ascending') : 'none');
    $th.contents().filter(function() {
      return this.nodeType === 3;
    }).remove();
    var lbl = mapperWheelSortHeaderLabel(col);
    if (lbl) { $th.prepend(document.createTextNode(lbl + ' ')); }
    $th.find('.mapper-core-sort-indicator').text(active ? (dir === 'desc' ? '▼' : '▲') : '↕');
  });
}

function mapperWheelRowSortKeys($tr, ix) {
  var entry = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
  var unmatched = $tr.data('mapperCoreUnmatchedRaw');
  var unmatchedStr = unmatched != null ? String(unmatched) : '';
  var ta = mapperWheelNormalizeTypedReceptorInput($tr.find('.mapper-core-in-receptor').val() || '');
  var val = ($tr.find('.mapper-wheel-in-value').val() || '');
  var displayText = ($tr.find('.mapper-core-receptor-html-view:visible').text() || '').replace(/\s+/g, ' ').trim();
  var meta = entry && MAPPER_CORE_ENTRY_META ? (MAPPER_CORE_ENTRY_META[entry] || {}) : {};
  var receptorSort = displayText || meta.gene || meta.name_plain || entry || unmatchedStr || ta;
  return {
    row: { entry: entry, unmatched: unmatchedStr, ta: ta, value: val },
    receptorSort: receptorSort,
    valueSort: mapperWheelTrimmedValueFromRowRaw(val),
    originalIndex: ix
  };
}

function mapperWheelSortSerializedRows(column, explicitDir) {
  column = column === 'value' ? 'value' : 'receptor';
  if (explicitDir === 'asc' || explicitDir === 'desc') {
    MAPPER_WHEEL_SORT_STATE.column = column;
    MAPPER_WHEEL_SORT_STATE.dir = explicitDir;
  } else if (MAPPER_WHEEL_SORT_STATE.column === column) {
    MAPPER_WHEEL_SORT_STATE.dir = MAPPER_WHEEL_SORT_STATE.dir === 'asc' ? 'desc' : 'asc';
  } else {
    MAPPER_WHEEL_SORT_STATE.column = column;
    MAPPER_WHEEL_SORT_STATE.dir = 'asc';
  }

  var rows = [];
  $('#mapper-wheel-input-tbody tr').each(function(ix) {
    var $tr = $(this);
    if (mapperWheelRowIsContentEmpty($tr)) {
      return;
    }
    rows.push(mapperWheelRowSortKeys($tr, ix));
  });
  rows.sort(function(a, b) {
    var primary = column === 'value' ? 'valueSort' : 'receptorSort';
    var secondary = column === 'value' ? 'receptorSort' : 'valueSort';
    var c = mapperWheelNaturalSortCompare(a[primary], b[primary]);
    if (c !== 0) {
      return MAPPER_WHEEL_SORT_STATE.dir === 'desc' ? -c : c;
    }
    c = mapperWheelNaturalSortCompare(a[secondary], b[secondary]);
    if (c !== 0) {
      return c;
    }
    return a.originalIndex - b.originalIndex;
  });
  mapperWheelRestoreRowsFromSerialized(rows.map(function(x) { return x.row; }));
  mapperWheelCaptureSnapshotForActiveMode();
  mapperWheelUpdateSortHeaders();
}

function mapperWheelDestroyLabelColorSpectrumWidgets() {
  $('#mapper-wheel-label-color-pickers .mapper-core-lcat-spectrum').each(function() {
    var $inp = $(this);
    try {
      if ($inp.data('spectrum')) {
        $inp.spectrum('destroy');
      }
    } catch (eIgnore) {}
  });
}

/** Nudge Spectrum popup so it stays inside the viewport (close control remains reachable). */
function mapperWheelClampSpectrumPickerToViewport($input) {
  if (!$input || !$input.length) {
    return;
  }
  function adjust() {
    var $container;
    try {
      $container = $input.spectrum('container');
    } catch (eAdj) {
      return;
    }
    if (!$container || !$container.length || $container.hasClass('sp-hidden')) {
      return;
    }
    var pad = 10;
    var sl = window.pageXOffset || document.documentElement.scrollLeft || 0;
    var st = window.pageYOffset || document.documentElement.scrollTop || 0;
    var vw = window.innerWidth || document.documentElement.clientWidth || 0;
    var vh = window.innerHeight || document.documentElement.clientHeight || 0;
    var maxR = sl + vw - pad;
    var maxB = st + vh - pad;
    var off = $container.offset();
    if (!off) {
      return;
    }
    var left = off.left;
    var top = off.top;
    var cw = $container.outerWidth();
    var ch = $container.outerHeight();
    if (left + cw > maxR) {
      left = Math.max(sl + pad, maxR - cw);
    }
    if (top + ch > maxB) {
      top = Math.max(st + pad, maxB - ch);
    }
    if (left < sl + pad) {
      left = sl + pad;
    }
    if (top < st + pad) {
      top = st + pad;
    }
    $container.offset({ top: top, left: left });
  }
  window.requestAnimationFrame(function() {
    window.requestAnimationFrame(adjust);
  });
}

function mapperWheelDistinctSortedLabelsFromGrid() {
  var seen = {};
  mapperWheelCollectRows().forEach(function(r) {
    var t = mapperWheelTrimmedValueFromRowRaw(r.value_raw);
    if (t) {
      seen[t] = true;
    }
  });
  var arr = Object.keys(seen);
  arr.sort(function(a, b) {
    return a.localeCompare(b, undefined, { numeric: true, sensitivity: 'base' });
  });
  return arr;
}

function mapperWheelEnsureMapsForDistinctLabels(sortedLabels) {
  (sortedLabels || []).forEach(function(lbl) {
    if (!lbl) {
      return;
    }
    if (!Object.prototype.hasOwnProperty.call(MAPPER_CORE_LABEL_COLORS, lbl)) {
      MAPPER_CORE_LABEL_COLORS[lbl] = mapperWheelDefaultHexForLabelKey(lbl);
    }
    if (!Object.prototype.hasOwnProperty.call(MAPPER_CORE_LABEL_ENABLED, lbl)) {
      MAPPER_CORE_LABEL_ENABLED[lbl] = true;
    }
  });
}

function mapperWheelSetRowPickerVisual($rowPicker, hex, skipSpectrumSet) {
  if (!$rowPicker || !$rowPicker.length || !hex) {
    return;
  }
  $rowPicker.val(hex);
  $rowPicker.css('background-color', hex);
  $rowPicker.next('.mapper-wheel-row-swatch-replacer').find('.sp-preview-inner').css('background-color', hex);
  if (!skipSpectrumSet && $rowPicker.data('spectrum')) {
    try {
      $rowPicker.spectrum('set', hex);
      $rowPicker.spectrum('enable').css('opacity', '1');
    } catch (eRowSet) {}
  }
}

function mapperWheelUpdateWheelLabelColorInPlace(lbl, nextHex) {
  if (!lbl || !nextHex || !GPCRome_WheelDict) {
    return;
  }
  function walk(node) {
    if (!node || typeof node !== 'object') {
      return;
    }
    if (!Array.isArray(node) && Object.prototype.hasOwnProperty.call(node, 'Data')) {
      if (mapperWheelTrimmedValueFromRowRaw(node.Data) === lbl) {
        node.Color = nextHex;
      }
    }
    Object.keys(node).forEach(function(k) {
      walk(node[k]);
    });
  }
  walk(GPCRome_WheelDict);
}

function mapperWheelCollectDisplayLabelsForTextValue(lbl) {
  var labels = {};
  var labelType = GPCRomes_styling && GPCRomes_styling.LabelType ? GPCRomes_styling.LabelType : 'Protein';
  function walk(node, keyName) {
    if (!node || typeof node !== 'object') {
      return;
    }
    if (!Array.isArray(node) && Object.prototype.hasOwnProperty.call(node, 'Data')) {
      if (mapperWheelTrimmedValueFromRowRaw(node.Data) === lbl) {
        if (labelType === 'Uniprot') {
          labels[String(node.EntryName || keyName || '')] = true;
        } else if (labelType === 'Entrez') {
          labels[String(node.Entrez || keyName || '')] = true;
        } else {
          labels[String(keyName || '')] = true;
        }
      }
    }
    Object.keys(node).forEach(function(k) {
      walk(node[k], k);
    });
  }
  walk(GPCRome_WheelDict, '');
  return labels;
}

function mapperWheelRecolorExistingWheelTextLabel(lbl, nextHex) {
  var svg = d3v4.select('#' + GPCRome_location + '_svg');
  if (!svg || svg.empty()) {
    return false;
  }
  var labelSet = mapperWheelCollectDisplayLabelsForTextValue(lbl);
  var changed = false;
  svg.selectAll('path[class^="large-hollow-pie-"], path[class*=" large-hollow-pie-"]')
    .filter(function(d) {
      return d && Object.prototype.hasOwnProperty.call(labelSet, String(d.data || ''));
    })
    .style('fill', function() {
      changed = true;
      return nextHex;
    })
    .style('stroke', 'black');
  return changed;
}

function mapperWheelApplyLabelColor(lbl, nextHex, sourceInput) {
  lbl = mapperWheelTrimmedValueFromRowRaw(lbl);
  if (!lbl || !nextHex) {
    return;
  }
  MAPPER_CORE_LABEL_COLORS[lbl] = nextHex;
  MAPPER_CORE_LABEL_ENABLED[lbl] = true;
  $('.mapper-core-lcat-spectrum').each(function() {
    var $inp = $(this);
    if ($inp[0] === sourceInput || $inp.attr('data-mapper-wheel-label') !== lbl || !$inp.data('spectrum')) {
      return;
    }
    try {
      $inp.spectrum('set', nextHex);
      $inp.spectrum('enable').css('opacity', '1');
    } catch (ePanel) {}
  });
  $('#mapper-wheel-input-tbody .mapper-core-row-color-picker').each(function() {
    var $rowPicker = $(this);
    var rowLbl = $rowPicker.attr('data-mapper-wheel-label') || mapperWheelTrimmedValueFromRowRaw($rowPicker.closest('tr').find('.mapper-wheel-in-value').val());
    if (rowLbl !== lbl) {
      return;
    }
    $rowPicker.attr('data-mapper-wheel-label', lbl);
    mapperWheelSetRowPickerVisual($rowPicker, nextHex, $rowPicker[0] === sourceInput);
  });
  mapperWheelUpdateWheelLabelColorInPlace(lbl, nextHex);
  if (!mapperWheelRecolorExistingWheelTextLabel(lbl, nextHex)) {
    updateGPCRome();
  }
  mapperWheelScheduleRebuildPlot();
}

function mapperWheelDestroyRowColorSpectrum($sw) {
  if (!$sw || !$sw.length || !$sw.data('spectrum')) {
    return;
  }
  try {
    $sw.spectrum('destroy');
  } catch (eDestroy) {}
}

function mapperWheelProtectRowColorPickerFromClose($sw) {
  if (!$sw || !$sw.length || !$sw.data('spectrum')) {
    return;
  }
  var events = 'mousedown.mapperWheelRowColor mouseup.mapperWheelRowColor click.mapperWheelRowColor touchstart.mapperWheelRowColor touchend.mapperWheelRowColor pointerdown.mapperWheelRowColor pointerup.mapperWheelRowColor';
  function stopPickerEvent(e) {
    e.stopPropagation();
  }
  try {
    $sw.spectrum('container').off(events).on(events, stopPickerEvent);
  } catch (eContainer) {}
  $('.mapper-core-wheel-left')
    .off(events, '.mapper-wheel-row-swatch-sp-container, .mapper-wheel-row-swatch-replacer')
    .on(events, '.mapper-wheel-row-swatch-sp-container, .mapper-wheel-row-swatch-replacer', stopPickerEvent);
}

function mapperWheelInstallRowColorClosePolicy($sw) {
  if (!$sw || !$sw.length || !$sw.data('spectrum')) {
    return;
  }
  var downEvent = 'mousedown.mapperWheelRowColorPolicy touchstart.mapperWheelRowColorPolicy pointerdown.mapperWheelRowColorPolicy';
  var pickerId = String($sw.data('spectrum') || Math.floor(Math.random() * 1000000));
  var ns = '.mapperWheelRowColorPolicy' + pickerId;

  $sw.off(downEvent).on(downEvent, function() {
    var $replacer = $sw.next('.mapper-wheel-row-swatch-replacer');
    if ($replacer.hasClass('sp-active')) {
      $sw.data('mapperWheelAllowSpectrumClose', true);
    } else {
      $sw.data('mapperWheelAllowSpectrumClose', false);
    }
  });

  $(document).off('mousedown' + ns + ' touchstart' + ns + ' pointerdown' + ns)
    .on('mousedown' + ns + ' touchstart' + ns + ' pointerdown' + ns, function(e) {
      var $target = $(e.target);
      var $container = $();
      try {
        $container = $sw.spectrum('container');
      } catch (eContainer) {}
      if (
        $target.closest('.mapper-wheel-row-swatch-sp-container').length ||
        ($container && $container.length && $target.closest($container).length) ||
        $target.closest($sw.next('.mapper-wheel-row-swatch-replacer')).length
      ) {
        return;
      }
      $sw.data('mapperWheelAllowSpectrumClose', true);
    });
}

function mapperWheelMaybeReopenRowColorPicker($sw) {
  if (!$sw || !$sw.length || !$sw.data('spectrum')) {
    return;
  }
  if ($sw.data('mapperWheelAllowSpectrumClose') === true || !mapperWheelIsTextMode()) {
    $sw.removeData('mapperWheelAllowSpectrumClose');
    return;
  }
  window.setTimeout(function() {
    if (!$sw.length || !$sw.closest('body').length || !$sw.data('spectrum')) {
      return;
    }
    try {
      $sw.spectrum('show');
      mapperWheelProtectRowColorPickerFromClose($sw);
    } catch (eShow) {}
  }, 0);
}

function mapperWheelEnsureRowColorSpectrum($sw, lbl, hex) {
  if (!$sw || !$sw.length || typeof $sw.spectrum !== 'function') {
    return;
  }
  $sw.attr('data-mapper-wheel-label', lbl);
  if (!$sw.data('spectrum')) {
    $sw.spectrum({
      color: hex,
      showPalette: true,
      showInput: true,
      showButtons: false,
      preferredFormat: 'hex',
      appendTo: '.mapper-core-wheel-left',
      containerClassName: 'mapper-core-lcat-sp-container mapper-wheel-row-swatch-sp-container',
      replacerClassName: 'mapper-core-lcat-replacer mapper-wheel-row-swatch-replacer',
      palette: [
        ['#000', '#FF0000', '#00FF00', '#0000FF', '#FFFF00'],
        ['#FF00FF', '#00FFFF', '#FFFFFF', '#C0C0C0', '#808080'],
        ['#800000', '#808000', '#008000', '#800080', '#008080'],
        ['#000080']
      ],
      show: function() {
        $sw.data('mapperWheelAllowSpectrumClose', false);
        mapperWheelRowColorPickerOpenCount += 1;
        mapperWheelProtectRowColorPickerFromClose($sw);
      },
      hide: function() {
        var allowClose = $sw.data('mapperWheelAllowSpectrumClose') === true || !mapperWheelIsTextMode();
        mapperWheelRowColorPickerOpenCount = Math.max(0, mapperWheelRowColorPickerOpenCount - 1);
        mapperWheelMaybeReopenRowColorPicker($sw);
        if (allowClose && mapperWheelRowColorPickerOpenCount === 0 && mapperWheelRowColorPickerPendingRebuild) {
          mapperWheelRowColorPickerPendingRebuild = false;
          mapperWheelRebuildPlotImmediate();
        }
      },
      change: function(c) {
        var activeLabel = mapperWheelTrimmedValueFromRowRaw($sw.closest('tr').find('.mapper-wheel-in-value').val());
        var nextHex = c && c.toHexString ? c.toHexString() : (c && c.toRgbString ? c.toRgbString() : null);
        mapperWheelApplyLabelColor(activeLabel, nextHex, $sw[0]);
      },
      move: function(c) {
        var activeLabel = mapperWheelTrimmedValueFromRowRaw($sw.closest('tr').find('.mapper-wheel-in-value').val());
        var nextHex = c && c.toHexString ? c.toHexString() : (c && c.toRgbString ? c.toRgbString() : null);
        mapperWheelApplyLabelColor(activeLabel, nextHex, $sw[0]);
      }
    });
    mapperWheelProtectRowColorPickerFromClose($sw);
    mapperWheelInstallRowColorClosePolicy($sw);
  } else {
    try {
      $sw.spectrum('set', hex);
      $sw.spectrum('enable').css('opacity', '1');
      mapperWheelProtectRowColorPickerFromClose($sw);
      mapperWheelInstallRowColorClosePolicy($sw);
    } catch (eSet) {}
  }
}

function mapperWheelSyncRowSwatch($tr) {
  if (!$tr || !$tr.length) {
    return;
  }
  var $sw = $tr.find('.mapper-wheel-label-swatch');
  if (!$sw.length) {
    return;
  }
  if (!mapperWheelIsTextMode()) {
    mapperWheelDestroyRowColorSpectrum($sw);
    $sw.removeAttr('data-mapper-wheel-label');
    $sw.css('background-color', 'transparent');
    return;
  }
  var lbl = mapperWheelTrimmedValueFromRowRaw($tr.find('.mapper-wheel-in-value').val());
  if (!lbl) {
    mapperWheelDestroyRowColorSpectrum($sw);
    $sw.removeAttr('data-mapper-wheel-label');
    $sw.css('background-color', '#f5f5f5');
    return;
  }
  mapperWheelEnsureMapsForDistinctLabels([lbl]);
  var hex = MAPPER_CORE_LABEL_COLORS[lbl] || mapperWheelDefaultHexForLabelKey(lbl);
  mapperWheelSetRowPickerVisual($sw, hex, true);
  mapperWheelEnsureRowColorSpectrum($sw, lbl, hex);
}

function mapperWheelRefreshAllLabelSwatches() {
  $('#mapper-wheel-input-tbody tr').each(function() {
    mapperWheelSyncRowSwatch($(this));
  });
}

function mapperWheelSpectrumGetHex(sel) {
  var $el = $(sel);
  if (!$el.length || typeof $el.spectrum !== 'function') {
    return null;
  }
  try {
    var c = $el.spectrum('get');
    return c ? c.toHexString() : null;
  } catch (eCatch) {
    return null;
  }
}

function mapperWheelCloneObj(obj) {
  try {
    return JSON.parse(JSON.stringify(obj));
  } catch (eC) {
    return obj;
  }
}

function mapperWheelCaptureNumericUiSnapshotExtras() {
  return {
    colorSchemeSelect: ($('#GPCRome_color_styling').val() || 'One'),
    stylingNums: {
      colorStart: GPCRomes_styling.colorStart,
      colorAvg: GPCRomes_styling.colorAvg,
      colorEnd: GPCRomes_styling.colorEnd,
      ColorSetup: GPCRomes_styling.ColorSetup,
      showIcon: GPCRomes_styling.showIcon,
      ShowLegend: GPCRomes_styling.ShowLegend,
      LabelType: GPCRomes_styling.LabelType
    },
    specHex: {
      min: mapperWheelSpectrumGetHex('#GPCRome_colorPicker_min'),
      avg: mapperWheelSpectrumGetHex('#GPCRome_colorPicker_avg'),
      max: mapperWheelSpectrumGetHex('#GPCRome_colorPicker_max')
    }
  };
}

function mapperWheelCaptureTextUiSnapshotExtras() {
  return {
    labelColors: mapperWheelCloneObj(MAPPER_CORE_LABEL_COLORS),
    labelEnabled: mapperWheelCloneObj(MAPPER_CORE_LABEL_ENABLED),
    legendLayout: mapperWheelCloneObj(GPCRomes_styling.LegendLayout || { mode: 'row', columns: '1', sorted: 'Vertically' }),
    stylingFrag: {
      showIcon: GPCRomes_styling.showIcon,
      ShowLegend: GPCRomes_styling.ShowLegend,
      LabelType: GPCRomes_styling.LabelType
    }
  };
}

function mapperWheelCaptureSnapshotForActiveMode() {
  var snap = {
    rows: mapperWheelSerializeDomRows(),
    receptorTallPlaceholderIntro: mapperWheelReceptorTallPlaceholderIntro
  };
  if (mapperWheelIsTextMode()) {
    var xt = mapperWheelCaptureTextUiSnapshotExtras();
    Object.assign(snap, xt);
    MAPPER_WHEEL_MODE_SNAPSHOTS.text = snap;
  } else {
    var nx = mapperWheelCaptureNumericUiSnapshotExtras();
    Object.assign(snap, nx);
    MAPPER_WHEEL_MODE_SNAPSHOTS.numeric = snap;
  }
}

function mapperWheelEnsureOppositeModeSeed(targetMode) {
  if (MAPPER_WHEEL_MODE_SNAPSHOTS[targetMode] != null) {
    return;
  }
  if (targetMode === 'text') {
    var baseSnap = MAPPER_WHEEL_MODE_SNAPSHOTS.numeric || {};
    MAPPER_WHEEL_MODE_SNAPSHOTS.text = {
      rows: [],
      receptorTallPlaceholderIntro: baseSnap.receptorTallPlaceholderIntro,
      labelColors: {},
      labelEnabled: {},
      legendLayout: { mode: 'row', columns: '1', sorted: 'Vertically' },
      stylingFrag: {
        showIcon: GPCRomes_styling.showIcon,
        ShowLegend: GPCRomes_styling.ShowLegend,
        LabelType: GPCRomes_styling.LabelType
      }
    };
  } else if (targetMode === 'numeric') {
    var tSnap = MAPPER_WHEEL_MODE_SNAPSHOTS.text || {};
    var n0 = MAPPER_WHEEL_MODE_SNAPSHOTS.numeric || {};
    MAPPER_WHEEL_MODE_SNAPSHOTS.numeric = {
      rows: [],
      receptorTallPlaceholderIntro: tSnap.receptorTallPlaceholderIntro,
      colorSchemeSelect: n0.colorSchemeSelect || ($('#GPCRome_color_styling').val() || 'One'),
      stylingNums: n0.stylingNums ? mapperWheelCloneObj(n0.stylingNums) : null,
      specHex: mapperWheelCloneObj(n0.specHex || {})
    };
    if (!MAPPER_WHEEL_MODE_SNAPSHOTS.numeric.stylingNums) {
      delete MAPPER_WHEEL_MODE_SNAPSHOTS.numeric.stylingNums;
    }
  }
}

function mapperWheelHydrateSpectrumTriplet(specHexObj) {
  if (!specHexObj) {
    return;
  }
  function set(id, hex) {
    if (!hex) {
      return;
    }
    var $e = $('#' + id);
    if ($e.length && $e.data('spectrum')) {
      try {
        $e.spectrum('set', hex);
      } catch (eI) {}
    }
  }
  set('GPCRome_colorPicker_min', specHexObj.min);
  set('GPCRome_colorPicker_avg', specHexObj.avg);
  set('GPCRome_colorPicker_max', specHexObj.max);
}

function mapperWheelRestoreRowsFromSerialized(rows) {
  mapperWheelRestoringDomRows = true;
  try {
  window.clearTimeout(mapperWheelPlotRebuildTimer);
  mapperWheelPlotRebuildTimer = null;
  mapperWheelDestroyRowReceptorWidgets();
  $('#mapper-wheel-input-tbody').empty();
  var list = rows && rows.length ? rows : [{ entry: '', unmatched: '', ta: '', value: '' }];
  list.forEach(function(ro) {
    mapperWheelAppendRow('', ro.value != null ? String(ro.value) : '', true);
    var $tr = $('#mapper-wheel-input-tbody tr').last();
    if ((ro.entry || '').trim()) {
      mapperWheelApplyCanonicalToRowSelect($tr, ro.entry.trim());
    } else if ((ro.unmatched || '').trim()) {
      mapperWheelSetReceptorRawOnRow($tr, String(ro.unmatched).trim());
    } else if ((ro.ta || '').trim()) {
      mapperWheelDestroyReceptorAutocomplete($tr.find('.mapper-core-in-receptor'));
      var $td = $tr.find('td.mapper-core-receptor-cell');
      $td.empty();
      var $inp2 = mapperWheelCreateReceptorCell($td);
      mapperWheelBindReceptorAutocomplete($inp2);
      $tr.find('.mapper-core-receptor-entry').val('');
      $inp2.val(ro.ta).show();
      $tr.removeClass('mapper-core-row-invalid');
      $tr.removeData('mapperCoreUnmatchedRaw');
      mapperWheelSyncReceptorClearBtn($tr);
    }
    mapperWheelSetValueCellRawOnRow($tr, ro.value);
    mapperWheelSyncRowSwatch($tr);
    mapperWheelSyncValueCellHint($tr);
  });
  mapperWheelSyncFirstRowReceptorPlaceholder();
  mapperWheelSyncReceptorTextareaHeights();
  mapperWheelEnsureTrailingBlankRow();
  } finally {
    mapperWheelRestoringDomRows = false;
  }
}

function mapperWheelSyncClearDropdown() {
  var isText = mapperWheelIsTextMode();
  $('#mapper-wheel-clear-mode-label').text(isText ? 'Categories' : 'Numbers');
  $('.mapper-wheel-clear-col-num').toggle(!isText);
  $('.mapper-wheel-clear-col-text').toggle(isText);
}

function mapperWheelApplyChromeForInputMode() {
  var text = mapperWheelIsTextMode();
  $('#mapper-core-input-table').toggleClass('mapper-core-text-mode', text);
  $('#mapper-wheel-datalabels-wrap').toggleClass('mapper-wheel-hide-when-numeric', !text);
  var $cw = $('.mapper-core-colors-panel-wrap');
  $cw.toggleClass('mapper-core-colors-show-numeric', !text);
  $cw.toggleClass('mapper-core-colors-show-text', text);

  $('#mapper-wheel-mode-numeric-btn')
    .toggleClass('btn-primary', !text)
    .toggleClass('btn-default', text)
    .toggleClass('active', !text)
    .attr('aria-pressed', !text ? 'true' : 'false');
  $('#mapper-wheel-mode-labels-btn')
    .toggleClass('btn-primary', text)
    .toggleClass('btn-default', !text)
    .toggleClass('active', text)
    .attr('aria-pressed', text ? 'true' : 'false');
  // Update value column header label
  $('.mapper-wheel-value-col-label').text(text ? 'Category' : 'Number');
  mapperWheelUpdateSortHeaders();
  if (typeof mapperWheelSyncClearDropdown === 'function') { mapperWheelSyncClearDropdown(); }
}

function mapperWheelGraphicalLegendSyncActiveButtons() {
  var graphicalLegendOptions = document.querySelectorAll('.graphical-legend-option');
  var iconShown = !!(GPCRomes_styling.showIcon);
  var legShown = !!(GPCRomes_styling.ShowLegend);
  graphicalLegendOptions.forEach(function(btn) {
    var iconValue = btn.getAttribute('data-icon') === 'true';
    var legendValue = btn.getAttribute('data-legend') === 'true';
    var match = iconValue === iconShown && legendValue === legShown;
    btn.classList.remove('btn-primary', 'btn-outline-primary');
    if (match) {
      btn.classList.add('btn-primary');
    } else {
      btn.classList.add('btn-outline-primary');
    }
  });
}

function mapperWheelSetWheelInputMode(nextMode, skipRedraw) {
  if (nextMode !== 'numeric' && nextMode !== 'text') {
    return;
  }
  if (nextMode === mapperWheelWheelInputMode) {
    return;
  }
  mapperWheelCaptureSnapshotForActiveMode();
  mapperWheelEnsureOppositeModeSeed(nextMode);
  mapperWheelWheelInputMode = nextMode;
  var snap = MAPPER_WHEEL_MODE_SNAPSHOTS[nextMode] || {};
  mapperWheelReceptorTallPlaceholderIntro = snap.receptorTallPlaceholderIntro !== undefined ? !!snap.receptorTallPlaceholderIntro : mapperWheelReceptorTallPlaceholderIntro;
  mapperWheelRestoreRowsFromSerialized(snap.rows || []);

  if (mapperWheelIsTextMode()) {
    MAPPER_CORE_LABEL_COLORS = snap.labelColors && typeof snap.labelColors === 'object' ? mapperWheelCloneObj(snap.labelColors) : {};
    MAPPER_CORE_LABEL_ENABLED = snap.labelEnabled && typeof snap.labelEnabled === 'object' ? mapperWheelCloneObj(snap.labelEnabled) : {};
    if (snap.legendLayout && typeof snap.legendLayout === 'object') {
      GPCRomes_styling.LegendLayout = mapperWheelCloneObj(snap.legendLayout);
    }
    if (snap.stylingFrag && typeof snap.stylingFrag === 'object') {
      if (snap.stylingFrag.LabelType !== undefined) {
        GPCRomes_styling.LabelType = snap.stylingFrag.LabelType;
      }
      if (snap.stylingFrag.showIcon !== undefined) {
        GPCRomes_styling.showIcon = !!snap.stylingFrag.showIcon;
      }
      if (snap.stylingFrag.ShowLegend !== undefined) {
        GPCRomes_styling.ShowLegend = !!snap.stylingFrag.ShowLegend;
      }
    }
    GPCRomes_styling.DataType = 'Text';
    GPCRomes_styling.FontStyle = 'Arial';
    if (!GPCRomes_styling.LegendLayout) {
      GPCRomes_styling.LegendLayout = { mode: 'row', columns: '1', sorted: 'Vertically' };
    }
    mapperWheelSyncLegendLayoutUiFromStyling();
  } else {
    if (snap.stylingNums && typeof snap.stylingNums === 'object') {
      GPCRomes_styling.colorStart = snap.stylingNums.colorStart != null ? snap.stylingNums.colorStart : GPCRomes_styling.colorStart;
      GPCRomes_styling.colorAvg = snap.stylingNums.colorAvg != null ? snap.stylingNums.colorAvg : GPCRomes_styling.colorAvg;
      GPCRomes_styling.colorEnd = snap.stylingNums.colorEnd != null ? snap.stylingNums.colorEnd : GPCRomes_styling.colorEnd;
      GPCRomes_styling.ColorSetup = snap.stylingNums.ColorSetup || GPCRomes_styling.ColorSetup;
      if (snap.stylingNums.LabelType !== undefined) {
        GPCRomes_styling.LabelType = snap.stylingNums.LabelType;
      }
      if (snap.stylingNums.showIcon !== undefined) {
        GPCRomes_styling.showIcon = !!snap.stylingNums.showIcon;
      }
      if (snap.stylingNums.ShowLegend !== undefined) {
        GPCRomes_styling.ShowLegend = !!snap.stylingNums.ShowLegend;
      }
    }
    if (snap.colorSchemeSelect) {
      $('#GPCRome_color_styling').val(snap.colorSchemeSelect).trigger('change.select2');
    }
    GPCRomes_styling.DataType = 'Numeric';
    window.setTimeout(function() {
      mapperWheelHydrateSpectrumTriplet(snap.specHex || {});
      if (typeof window.Mapper20UpdateNumericColorPanels === 'function') {
        window.Mapper20UpdateNumericColorPanels();
      }
    }, 0);
    mapperWheelRecalcNumericStyling();
  }

  mapperWheelApplyChromeForInputMode();
  mapperWheelRefreshResolvedReceptorDisplaysInGrid();
  mapperWheelGraphicalLegendSyncActiveButtons();
  mapperWheelRebuildLabelColorCustomizePanel(false);
  if (!skipRedraw) {
    mapperWheelRebuildPlotImmediate();
  }

  $('.GPCRome-label-btn').each(function() {
    var v = $(this).data('value') || '';
    $(this).toggleClass('btn-primary', v === GPCRomes_styling.LabelType).toggleClass('btn-outline-primary', v !== GPCRomes_styling.LabelType);
  });
}

function mapperWheelInitLegendLayoutSelectUi() {
  var $sel = $('#mapper-wheel-legend-layout-select');
  if (!$sel.length || $sel.children().length > 0) {
    return;
  }
  $sel.append('<option value="row">Continuous text</option>');
  var i;
  for (i = 1; i <= 4; i++) {
    $sel.append($('<option>').attr('value', 'columns-' + String(i)).text(String(i) + ' column(s)'));
  }
}

function mapperWheelSyncLegendLayoutUiFromStyling() {
  var ll = GPCRomes_styling.LegendLayout || { mode: 'row', columns: '1', sorted: 'Vertically' };
  if (!ll.sorted) {
    ll.sorted = 'Vertically';
  }
  var $sel = $('#mapper-wheel-legend-layout-select');
  if (!$sel.length) {
    return;
  }
  if (ll.mode === 'columns' && ll.columns != null) {
    $sel.val('columns-' + String(parseInt(ll.columns, 10) || 1));
  } else {
    $sel.val('row');
  }
  var sortBtn = $('#mapper-wheel-toggle-legend-sorting');
  if (sortBtn.length) {
    sortBtn.text(ll.sorted === 'Horizontally' ? 'Horizontally' : 'Vertically');
  }
}

function mapperWheelApplyLegendLayoutSelectionFromDropdown() {
  var v = ($('#mapper-wheel-legend-layout-select').val() || 'row');
  if (!GPCRomes_styling.LegendLayout) {
    GPCRomes_styling.LegendLayout = {};
  }
  var sortedKeep = GPCRomes_styling.LegendLayout.sorted || 'Vertically';
  if (v === 'row') {
    GPCRomes_styling.LegendLayout.mode = 'row';
    GPCRomes_styling.LegendLayout.columns = '1';
  } else if (typeof v === 'string' && v.indexOf('columns-') === 0) {
    GPCRomes_styling.LegendLayout.mode = 'columns';
    GPCRomes_styling.LegendLayout.columns = String(parseInt(v.split('-')[1], 10) || 1);
  }
  GPCRomes_styling.LegendLayout.sorted = sortedKeep;
}

function mapperWheelRebuildLabelColorCustomizePanel(forceRedrawWheel) {
  mapperWheelDestroyLabelColorSpectrumWidgets();
  var $mount = $('#mapper-wheel-label-color-pickers');
  if (!$mount.length) {
    return;
  }
  $mount.empty();
  if (!mapperWheelIsTextMode()) {
    if (forceRedrawWheel) {
      mapperWheelRebuildPlotImmediate();
    }
    return;
  }

  var labels = mapperWheelDistinctSortedLabelsFromGrid();
  mapperWheelEnsureMapsForDistinctLabels(labels);

  var hdr = $('<div class="mapper-core-label-cat-head">').append(
    $('<input type="checkbox" id="mapper-wheel-label-cat-master" aria-label="Toggle all label colours enabled">'),
    $('<div class="mapper-core-label-cat-head-title">Categories</div>')
  );
  var grid = $('<div class="color-grid" id="mapper-wheel-label-color-grid"></div>');
  $mount.append(hdr, grid);

  function domId(ix, kind) {
    return 'mapperWheel_l_' + kind + '_' + String(ix);
  }

  function updateMasterCheckbox() {
    var allOn = labels.length > 0 && labels.every(function(l) {
      return MAPPER_CORE_LABEL_ENABLED[l] !== false;
    });
    var allOff = labels.length > 0 && labels.every(function(l) {
      return MAPPER_CORE_LABEL_ENABLED[l] === false;
    });
    var mx = $('#mapper-wheel-label-cat-master')[0];
    if (!mx) {
      return;
    }
    mx.checked = !!allOn;
    mx.indeterminate = !allOn && !allOff && labels.length > 0;
  }

  function refreshGridRowUi(ix, lbl) {
    var chkId = domId(ix, 'chk');
    var spId = domId(ix, 'spe');
    var en = MAPPER_CORE_LABEL_ENABLED[lbl] !== false;
    $('#' + chkId).prop('checked', en);
    var $inp = $('#' + spId);
    if (!$inp.length || !$inp.data('spectrum')) {
      return;
    }
    try {
      if (en) {
        $inp.spectrum('enable').css('opacity', '1');
      } else {
        $inp.spectrum('disable').css('opacity', '0.5');
      }
    } catch (eR) {}
  }

  labels.forEach(function(lbl, ix) {
    var chkId = domId(ix, 'chk');
    var pickerId = domId(ix, 'spe');
    var hex = MAPPER_CORE_LABEL_COLORS[lbl] || mapperWheelDefaultHexForLabelKey(lbl);
    MAPPER_CORE_LABEL_COLORS[lbl] = hex;
    var en = MAPPER_CORE_LABEL_ENABLED[lbl] !== false;

    var $spe = $('<input type="text">')
      .addClass('mapper-core-lcat-spectrum form-control input-sm')
      .attr('id', pickerId)
      .attr('data-mapper-wheel-label', lbl);
    grid.append($('<div class="color-item">').append(
      $('<input type="checkbox">').attr('id', chkId).prop('checked', !!en),
      $('<label>').attr('for', chkId).addClass('color-label').text(lbl),
      $spe
    ));

    $('#' + chkId).on('change', function() {
      MAPPER_CORE_LABEL_ENABLED[lbl] = !!$(this).prop('checked');
      refreshGridRowUi(ix, lbl);
      mapperWheelRefreshAllLabelSwatches();
      mapperWheelScheduleRebuildPlot();
    });

    $spe.spectrum({
      color: hex,
      showPalette: true,
      showInput: true,
      showButtons: false,
      preferredFormat: 'hex',
      appendTo: '#mapper-wheel-colors-dropdown-menu',
      containerClassName: 'mapper-core-lcat-sp-container',
      replacerClassName: 'mapper-core-lcat-replacer',
      palette: [
        ['#000', '#FF0000', '#00FF00', '#0000FF', '#FFFF00'],
        ['#FF00FF', '#00FFFF', '#FFFFFF', '#C0C0C0', '#808080']
      ],
      show: function() {
        mapperWheelClampSpectrumPickerToViewport($spe);
      },
      change: function(c) {
        if (MAPPER_CORE_LABEL_ENABLED[lbl] === false) {
          return;
        }
        var nextHex = hex;
        if (c && c.toHexString) {
          nextHex = c.toHexString();
        } else if (c && c.toRgbString) {
          nextHex = c.toRgbString();
        }
        mapperWheelApplyLabelColor(lbl, nextHex || hex, $spe[0]);
      },
      move: function(c) {
        if (MAPPER_CORE_LABEL_ENABLED[lbl] === false) {
          return;
        }
        var nextHex = hex;
        if (c && c.toHexString) {
          nextHex = c.toHexString();
        } else if (c && c.toRgbString) {
          nextHex = c.toRgbString();
        }
        mapperWheelApplyLabelColor(lbl, nextHex || hex, $spe[0]);
      }
    });

    if (!en) {
      $spe.spectrum('disable').css('opacity', '0.5');
    }
  });

  $('#mapper-wheel-label-cat-master').on('change', function() {
    var setOn = !!$(this).prop('checked');
    labels.forEach(function(lbl0, ix0) {
      MAPPER_CORE_LABEL_ENABLED[lbl0] = setOn;
      $('#' + domId(ix0, 'chk')).prop('checked', setOn);
      refreshGridRowUi(ix0, lbl0);
    });
    mapperWheelRefreshAllLabelSwatches();
    mapperWheelScheduleRebuildPlot();
    updateMasterCheckbox();
  });

  updateMasterCheckbox();

  if (forceRedrawWheel) {
    mapperWheelRebuildPlotImmediate();
  }
}

function mapperWheelWheelLabelType() {
  if (typeof GPCRomes_styling !== 'undefined' && GPCRomes_styling && GPCRomes_styling.LabelType) {
    return GPCRomes_styling.LabelType;
  }
  return 'Protein';
}

function mapperWheelFallbackUniprotFromEntry(entryId, meta) {
  var sid = entryId != null ? String(entryId).trim() : '';
  var u = meta && meta.uniprot ? String(meta.uniprot).trim().toUpperCase() : '';
  if (u) {
    return u;
  }
  if (sid.indexOf('_') !== -1) {
    return sid.split('_')[0].trim().toUpperCase();
  }
  return sid.toUpperCase();
}

function mapperWheelResolvedReceptorDisplayHtml(entryId) {
  var sid = entryId != null ? String(entryId).trim() : '';
  if (!sid) {
    return '';
  }
  var meta = MAPPER_CORE_ENTRY_META[sid];
  var lt = mapperWheelWheelLabelType();
  if (lt === 'Entrez') {
    var g = meta && meta.gene ? String(meta.gene).trim() : '';
    return g ? $('<span/>').text(g).html() : $('<span/>').text(sid).html();
  }
  if (lt === 'Uniprot') {
    var uni = mapperWheelFallbackUniprotFromEntry(sid, meta);
    return $('<span/>').text(uni || sid).html();
  }
  if (meta && meta.name_html) {
    return String(meta.name_html);
  }
  return $('<span/>').text(sid).html();
}

/** Plain string for autocomplete `label` (accessibility / matching). */
function mapperWheelReceptorAutocompleteLabelPlain(it) {
  if (!it || it.id == null) {
    return '';
  }
  var meta = MAPPER_CORE_ENTRY_META[it.id];
  var lt = mapperWheelWheelLabelType();
  if (lt === 'Entrez') {
    var g = meta && meta.gene ? String(meta.gene).trim() : '';
    return g || it.name_plain || it.text || it.id;
  }
  if (lt === 'Uniprot') {
    return mapperWheelFallbackUniprotFromEntry(it.id, meta) || it.id;
  }
  return it.name_plain || it.text || it.id;
}

function mapperWheelReceptorEditSeed(id) {
  if (!id) {
    return '';
  }
  var meta = MAPPER_CORE_ENTRY_META[id] || {};
  var lt = mapperWheelWheelLabelType();
  if (lt === 'Entrez') {
    return meta.gene || meta.name_plain || '';
  }
  if (lt === 'Uniprot') {
    return mapperWheelFallbackUniprotFromEntry(id, meta);
  }
  return meta.name_plain || '';
}

function mapperWheelRefreshResolvedReceptorDisplaysInGrid() {
  $('#mapper-wheel-input-tbody tr').each(function() {
    var $tr = $(this);
    var id = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
    if (!id) {
      return;
    }
    var $view = $tr.find('.mapper-core-receptor-html-view');
    if (!$view.length || !$view.is(':visible')) {
      return;
    }
    $view.html(mapperWheelResolvedReceptorDisplayHtml(id));
  });
}

function mapperWheelQuoteForMsg(s) {
  var t = s == null ? '' : String(s);
  if (t.trim() === '') {
    return '(empty)';
  }
  return '"' + t + '"';
}

function mapperWheelMessageInvalidReceptor(x) {
  return mapperWheelQuoteForMsg(x) + ' is an invalid receptor for the GPCR plot — try the search option for correct receptor name input.';
}

function mapperWheelMessageInvalidNumeric(x) {
  return mapperWheelQuoteForMsg(x) + ' is not a valid numeric value.';
}

function mapperWheelMessageReceptorRequiredForValue() {
  return 'Enter a receptor (use search) before applying a numeric value in this row.';
}

function mapperWheelMessageNumericMissingForReceptor() {
  return 'Enter a numeric value for this receptor to show it on the plot.';
}

var mapperWheelPlotDebounceMs = 200;
var mapperWheelPlotRebuildTimer = null;
var mapperWheelRestoringDomRows = false;
var mapperWheelLabelPanelDebouncer = null;
var mapperWheelRowColorPickerOpenCount = 0;
var mapperWheelRowColorPickerPendingRebuild = false;

function mapperWheelScheduleRebuildPlot() {
  if (mapperWheelRowColorPickerOpenCount > 0) {
    mapperWheelRowColorPickerPendingRebuild = true;
    return;
  }
  window.clearTimeout(mapperWheelPlotRebuildTimer);
  mapperWheelPlotRebuildTimer = window.setTimeout(function() {
    mapperWheelPlotRebuildTimer = null;
    mapperWheelRebuildPlotFromTableInternal();
  }, mapperWheelPlotDebounceMs);
}

function mapperWheelRebuildPlotImmediate() {
  if (mapperWheelRowColorPickerOpenCount > 0) {
    mapperWheelRowColorPickerPendingRebuild = true;
    return;
  }
  window.clearTimeout(mapperWheelPlotRebuildTimer);
  mapperWheelPlotRebuildTimer = null;
  mapperWheelRebuildPlotFromTableInternal();
}

function mapperWheelTrimmedValueFromRowRaw(valueRaw) {
  if (valueRaw === true || valueRaw === false) {
    return String(valueRaw).trim();
  }
  if (valueRaw === null || valueRaw === undefined) {
    return '';
  }
  if (typeof valueRaw === 'number') {
    return String(valueRaw);
  }
  return String(valueRaw).trim();
}

function mapperWheelNumericStringAcceptable(trimmed) {
  if (trimmed === '') {
    return true;
  }
  var n = Number(trimmed);
  return Number.isFinite(n);
}

/**
 * Turn locale-style numbers into dot-decimal form for JS Number().
 * 1,222.2 -> 1222.2 · 1.222,2 -> 1222.2 · 1,2 -> 1.2 · 12,345 -> 12345 (US-style thousands)
 */
function mapperWheelNormalizeLocalizedNumberString(raw) {
  var orig = raw == null ? '' : String(raw);
  var t = orig.trim();
  if (!t) {
    return '';
  }
  if (/[eE]/.test(t)) {
    return t;
  }
  var compact = t.replace(/\s+/g, '');
  var signCh = '';
  if (/^[+-]/.test(compact)) {
    signCh = compact.charAt(0);
    compact = compact.slice(1);
  }
  if (!compact.length) {
    return orig.trim();
  }
  if (!/[,.]/.test(compact)) {
    return signCh + compact;
  }
  if (!/^[\d,.]+$/.test(compact)) {
    return orig.trim();
  }

  var s = compact;
  var hasComma = s.indexOf(',') >= 0;
  var hasDot = s.indexOf('.') >= 0;

  if (hasComma && hasDot) {
    if (s.lastIndexOf(',') > s.lastIndexOf('.')) {
      s = s.replace(/\./g, '').replace(',', '.');
    } else {
      s = s.replace(/,/g, '');
    }
    return signCh + s;
  }

  if (hasComma && !hasDot) {
    var cp = s.split(',');
    if (cp.length > 2) {
      s = cp.join('');
    } else if (cp.length === 2) {
      var ca = cp[0];
      var cb = cp[1];
      if (!/^\d+$/.test(ca) || !/^\d+$/.test(cb)) {
        return orig.trim();
      }
      if (cb.length <= 2) {
        s = ca + '.' + cb;
      } else if (cb.length === 3) {
        s = ca + cb;
      } else {
        s = ca + cb;
      }
    }
    return signCh + s;
  }

  if (!hasComma && hasDot) {
    var dp = s.split('.');
    if (dp.length > 2 && dp.every(function(p) {
      return /^\d+$/.test(p);
    })) {
      var decFrag = dp.pop();
      if (decFrag.length <= 2) {
        s = dp.join('') + '.' + decFrag;
      } else {
        s = dp.join('');
      }
    }
    return signCh + s;
  }

  return signCh + s;
}

/** Returns a string acceptable to Number(...), preserving input when already valid or unrecognized. */
function mapperWheelNormalizedValueForJsNumber(trimmedStr) {
  if (trimmedStr === '') {
    return '';
  }
  if (/[eE]/.test(trimmedStr)) {
    return trimmedStr;
  }
  if (!/[,.]/.test(trimmedStr)) {
    return trimmedStr;
  }
  var n = mapperWheelNormalizeLocalizedNumberString(trimmedStr);
  return mapperWheelNumericStringAcceptable(n) ? n : trimmedStr;
}

function mapperWheelApplyOrangeMissingValueHints(rawRows, errReceptorRows, errValRows) {
  $('#mapper-wheel-input-tbody tr').each(function(ix) {
    var idx = ix + 1;
    var row = rawRows[ix];
    if (!row) {
      return;
    }
    if (errReceptorRows[idx] || errValRows[idx]) {
      return;
    }
    var rs = row.receptor_raw != null ? String(row.receptor_raw).trim() : '';
    var vTrim;
    if (mapperWheelIsTextMode()) {
      vTrim = mapperWheelTrimmedValueFromRowRaw(row.value_raw);
    } else {
      vTrim = mapperWheelNormalizedValueForJsNumber(mapperWheelTrimmedValueFromRowRaw(row.value_raw));
    }
    if (vTrim !== '') {
      return;
    }
    if (!mapperWheelNormalizeReceptorForWheel(rs)) {
      return;
    }
    var $tr = $(this);
    if ($tr.hasClass('mapper-core-row-invalid')) {
      return;
    }
    var $tdV = $tr.find('td.mapper-core-value-cell');
    mapperWheelSetCellWarnHint($tdV, mapperWheelIsTextMode() ? mapperWheelMessageLabelMissingForReceptor() : mapperWheelMessageNumericMissingForReceptor());
  });
}

function mapperWheelClearCellHoverHint($td) {
  if (!$td || !$td.length) {
    return;
  }
  $td.removeClass('mapper-wheel-cell-hint-error mapper-wheel-cell-hint-warn');
  $td.removeAttr('title');
}

function mapperWheelSetCellHoverHint($td, message) {
  if (!$td || !$td.length) {
    return;
  }
  $td.removeClass('mapper-wheel-cell-hint-warn');
  $td.addClass('mapper-wheel-cell-hint-error');
  $td.attr('title', message);
}

function mapperWheelSetCellWarnHint($td, message) {
  if (!$td || !$td.length) {
    return;
  }
  $td.removeClass('mapper-wheel-cell-hint-error');
  $td.addClass('mapper-wheel-cell-hint-warn');
  $td.attr('title', message);
}

function mapperWheelClearRowApplyHints($tr) {
  if (!$tr || !$tr.length) {
    return;
  }
  mapperWheelClearCellHoverHint($tr.find('td.mapper-core-receptor-cell'));
  mapperWheelClearCellHoverHint($tr.find('td.mapper-core-value-cell'));
}

function mapperWheelClearAllPlotRowUiHints() {
  $('#mapper-wheel-input-tbody tr').each(function() {
    var $tr = $(this);
    $tr.removeClass('mapper-wheel-row-error');
    mapperWheelClearRowApplyHints($tr);
  });
}

/** Live value-/label-cell hint while typing; full grid hints come from mapperWheelRebuildPlotFromTableInternal. */
function mapperWheelSyncValueCellHint($tr) {
  if (!$tr || !$tr.length) {
    return;
  }
  var $td = $tr.find('td.mapper-core-value-cell');
  var rawTrim = mapperWheelTrimmedValueFromRowRaw($tr.find('.mapper-wheel-in-value').val());
  if (mapperWheelIsTextMode()) {
    mapperWheelClearCellHoverHint($td);
    return;
  }
  var s = mapperWheelNormalizedValueForJsNumber(rawTrim);
  if (s === '') {
    mapperWheelClearCellHoverHint($td);
    return;
  }
  if (!mapperWheelNumericStringAcceptable(s)) {
    mapperWheelSetCellHoverHint($td, mapperWheelMessageInvalidNumeric(s));
  } else {
    mapperWheelClearCellHoverHint($td);
  }
}

/**
 * Validates rows like Apply, merges valid points into wheel data, redraws plot.
 * Applies red hover hints where Apply would; orange warn on value cell when receptor is OK but value empty (no conflicting row-level errors).
 * Canonicalizes receptors only when there are zero validation errors across the grid.
 */
function mapperWheelRebuildPlotFromTableInternal() {
  $('#mapper-core-wheel-errors').empty();
  mapperWheelClearAllPlotRowUiHints();

  var rawRows = mapperWheelCollectRows();
  if (rawRows.length > MAPPER_WHEEL_MAX_ROWS) {
    return;
  }

  var Data = {};
  var errors = [];

  if (!mapperWheelIsTextMode()) {

    rawRows.forEach(function(row, ix) {
      var idx = ix + 1;
      var rs = row.receptor_raw != null ? String(row.receptor_raw).trim() : '';
      var valueRaw = row.value_raw;

      var valueStrOrig = mapperWheelTrimmedValueFromRowRaw(valueRaw);
      var valueStr = mapperWheelNormalizedValueForJsNumber(valueStrOrig);

      var valueNonempty = valueStr !== '';

      if (!rs && !valueNonempty) {
        return;
      }
      if (!rs && valueNonempty) {
        errors.push({ row: idx, field: 'receptor', hint: mapperWheelMessageReceptorRequiredForValue() });
        return;
      }

      var receptor = mapperWheelNormalizeReceptorForWheel(rs);
      if (!receptor) {
        errors.push({ row: idx, field: 'receptor', hint: mapperWheelMessageInvalidReceptor(rs) });
        return;
      }

      if (!valueNonempty) {
        return;
      }
      if (!mapperWheelNumericStringAcceptable(valueStr)) {
        errors.push({ row: idx, field: 'value', hint: mapperWheelMessageInvalidNumeric(valueStrOrig) });
        return;
      }

      var valNum = Number(valueStr);
      Data[receptor] = { Value1: valNum };
    });
  } else {

    rawRows.forEach(function(row, ix) {
      var idx = ix + 1;
      var rs = row.receptor_raw != null ? String(row.receptor_raw).trim() : '';
      var valueRaw = row.value_raw;

      var labelTrim = mapperWheelTrimmedValueFromRowRaw(valueRaw);
      var labelNonempty = labelTrim !== '';

      if (!rs && !labelNonempty) {
        return;
      }
      if (!rs && labelNonempty) {
        errors.push({ row: idx, field: 'receptor', hint: mapperWheelMessageReceptorRequiredForLabelCell() });
        return;
      }

      var receptor = mapperWheelNormalizeReceptorForWheel(rs);
      if (!receptor) {
        errors.push({ row: idx, field: 'receptor', hint: mapperWheelMessageInvalidReceptor(rs) });
        return;
      }

      if (!labelNonempty) {
        return;
      }

      mapperWheelEnsureMapsForDistinctLabels([labelTrim]);
      var enabled = MAPPER_CORE_LABEL_ENABLED[labelTrim] !== false;
      var col = '#ffffff';
      if (enabled) {
        col = MAPPER_CORE_LABEL_COLORS[labelTrim] || mapperWheelDefaultHexForLabelKey(labelTrim);
        MAPPER_CORE_LABEL_COLORS[labelTrim] = col;
      }
      Data[receptor] = { Value1: labelTrim, Value2: col };
    });
  }

  errors.forEach(function(err) {
    if (err.row == null) {
      return;
    }
    var $trErr = $('#mapper-wheel-input-tbody tr').eq(err.row - 1);
    if (!$trErr.length) {
      return;
    }
    $trErr.addClass('mapper-wheel-row-error');
    var $tdR = $trErr.find('td.mapper-core-receptor-cell');
    var $tdV = $trErr.find('td.mapper-core-value-cell');
    if (err.field === 'receptor') {
      mapperWheelClearCellHoverHint($tdV);
      mapperWheelSetCellHoverHint($tdR, err.hint);
      mapperWheelUpgradeReceptorToSelect2($trErr);
    } else if (err.field === 'value') {
      mapperWheelClearCellHoverHint($tdR);
      mapperWheelSetCellHoverHint($tdV, err.hint);
    }
  });

  var errReceptorRows = {};
  var errValRows = {};
  errors.forEach(function(e) {
    if (e.row != null && e.field === 'receptor') {
      errReceptorRows[e.row] = true;
    }
    if (e.row != null && e.field === 'value') {
      errValRows[e.row] = true;
    }
  });

  if (errors.length === 0) {
    $('#mapper-wheel-input-tbody tr').each(function() {
      var $tr = $(this);
      var val = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
      if (!val) {
        val = mapperWheelNormalizeTypedReceptorInput($tr.find('.mapper-core-in-receptor').val() || '');
      }
      if (!val) {
        val = ($tr.data('mapperCoreUnmatchedRaw') || '').trim();
      }
      mapperWheelApplyCanonicalToRowSelect($tr, val);
    });
    rawRows = mapperWheelCollectRows();
  }

  mapperWheelApplyOrangeMissingValueHints(rawRows, errReceptorRows, errValRows);

  var base = JSON.parse(JSON.stringify(MAPPER_WHEEL_GPCROME_BASE));
  GPCRome_WheelDict = mapperWheelUpdateNestedGPCRomeData(base, Data);
  if (!mapperWheelIsTextMode()) {
    mapperWheelRecalcNumericStyling();
  }
  mapperWheelRefreshAllLabelSwatches();
  if (mapperWheelIsTextMode()) {
    mapperWheelScheduleLabelColorPanelRebuildSoon();
  }
  updateGPCRome();
}

/** Collapses line breaks when reading typed receptor paste (textarea allows wrapped placeholder editing). */
function mapperWheelNormalizeTypedReceptorInput(val) {
  if (val == null) {
    return '';
  }
  return String(val).replace(/\r\n/g, '\n').replace(/[\r\n]+/g, ' ').replace(/\s+/g, ' ').trim();
}

function mapperWheelAutoSizeReceptorTextarea($ta) {
  if (!$ta || !$ta.length || !$ta[0] || String($ta[0].tagName || '').toUpperCase() !== 'TEXTAREA') {
    return;
  }
  var el = $ta[0];
  var minPx = 44;
  var maxPx = 200;
  window.requestAnimationFrame(function() {
    el.style.height = 'auto';
    var sh = el.scrollHeight;
    var next = Math.max(minPx, Math.min(sh, maxPx));
    el.style.height = next + 'px';
    el.style.overflowY = sh > maxPx ? 'auto' : 'hidden';
  });
}

// Set once the compact-mode full-table cleanup below has run with nothing left to clean —
// avoids re-touching every row's textarea style on every keystroke once the table has
// settled into compact mode (that pass was a no-op past the initial tall->compact
// transition, but still forced a layout reflow across every row each time it ran).
var mapperWheelReceptorTextareaCompactClean = false;

function mapperWheelSyncReceptorTextareaHeights() {
  $('#mapper-core-input-table').toggleClass('mapper-core-receptors-compact', !mapperWheelReceptorTallPlaceholderIntro);
  if (mapperWheelReceptorTallPlaceholderIntro) {
    mapperWheelReceptorTextareaCompactClean = false;
    $('#mapper-wheel-input-tbody .mapper-core-in-receptor').each(function() {
      this.style.overflowY = '';
      mapperWheelAutoSizeReceptorTextarea($(this));
    });
    return;
  }
  if (mapperWheelReceptorTextareaCompactClean) {
    return;
  }
  $('#mapper-wheel-input-tbody .mapper-core-in-receptor').each(function() {
    this.style.height = '';
    this.style.overflowY = 'hidden';
  });
  mapperWheelReceptorTextareaCompactClean = true;
}

function mapperWheelCompactReceptorRowsAfterInput() {
  if (!mapperWheelReceptorTallPlaceholderIntro) {
    return;
  }
  mapperWheelReceptorTallPlaceholderIntro = false;
  mapperWheelSyncFirstRowReceptorPlaceholder();
  mapperWheelSyncReceptorTextareaHeights();
}

function mapperWheelFilterReceptorsLocal(term, limit) {
  var t = (term || '').trim().toUpperCase();
  if (!t) {
    return [];
  }
  var filtered = receptorSelect2Data.filter(function(item) {
    var st = (item.search_text || item.text || '').toUpperCase();
    if (st.indexOf(t) !== -1) {
      return true;
    }
    if (item.id && String(item.id).toUpperCase().indexOf(t) !== -1) {
      return true;
    }
    return false;
  });
  filtered.sort(function(a, b) {
    function rank(item) {
      var id = String(item.id || '').toUpperCase();
      var st = (item.search_text || item.text || '').toUpperCase();
      if (id.indexOf(t) === 0) {
        return 0;
      }
      if (st.indexOf(t) === 0) {
        return 1;
      }
      if (id.indexOf(t) !== -1) {
        return 2;
      }
      if (st.indexOf(t) !== -1) {
        return 3;
      }
      return 4;
    }
    var ra = rank(a);
    var rb = rank(b);
    if (ra !== rb) {
      return ra - rb;
    }
    return String(a.text || '').localeCompare(String(b.text || ''));
  });
  var lim = limit || 80;
  if (filtered.length > lim) {
    filtered = filtered.slice(0, lim);
  }
  return filtered.map(function(item) {
    return {
      id: item.id,
      text: item.text,
      name_html: item.name_html || '',
      name_plain: item.name_plain || ''
    };
  });
}

function mapperWheelSyncReceptorClearBtn($tr) {
  var $wrap = $tr.find('.mapper-core-receptor-input-wrap');
  if (!$wrap.length) {
    return;
  }
  var entry = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
  var iv = mapperWheelNormalizeTypedReceptorInput($tr.find('.mapper-core-in-receptor').val() || '');
  var um = $tr.data('mapperCoreUnmatchedRaw');
  var umS = um != null ? String(um).trim() : '';
  var $view = $tr.find('.mapper-core-receptor-html-view');
  var viewOn = $view.is(':visible') && ($view.text() || '').replace(/\s+/g, ' ').trim().length > 0;
  var has = !!(entry || iv || umS || viewOn);
  $wrap.toggleClass('is-empty', !has);
}

function mapperWheelCreateReceptorCell($tdR) {
  var $hid = $('<input type="hidden" class="mapper-core-receptor-entry" value="">');
  var $wrap = $('<div class="mapper-core-receptor-input-wrap is-empty">');
  var $inp = $('<textarea class="form-control input-sm mapper-core-in-receptor" rows="1" autocomplete="off" spellcheck="false"></textarea>');
  var $view = $('<div class="form-control input-sm mapper-core-receptor-html-view" tabindex="0"></div>');
  var $clr = $('<button type="button" class="mapper-core-receptor-clear" aria-label="Clear receptor" title="Clear receptor">&times;</button>');
  $wrap.append($inp, $view, $clr);
  $tdR.append($hid, $wrap);
  $view.hide();
  return $inp;
}

function mapperWheelSetReceptorResolvedUI($tr, entryId) {
  var sid = entryId != null ? String(entryId).trim() : '';
  var $inp = $tr.find('.mapper-core-in-receptor');
  var $hid = $tr.find('.mapper-core-receptor-entry');
  var $view = $tr.find('.mapper-core-receptor-html-view');
  mapperWheelDestroyReceptorAutocomplete($inp);
  if (!sid) {
    $hid.val('');
    $view.hide().empty();
    $inp.val('').show();
    mapperWheelBindReceptorAutocomplete($inp);
    mapperWheelSyncReceptorClearBtn($tr);
    return;
  }
  mapperWheelCompactReceptorRowsAfterInput();
  $hid.val(sid);
  $view.html(mapperWheelResolvedReceptorDisplayHtml(sid));
  $inp.val('');
  $inp.hide();
  $view.show();
  mapperWheelBindReceptorAutocomplete($inp);
  mapperWheelSyncReceptorClearBtn($tr);
}

function mapperWheelFocusValueCellForRow($tr) {
  if (!$tr || !$tr.length) {
    return;
  }
  window.setTimeout(function() {
    var $value = $tr.find('.mapper-wheel-in-value');
    if (!$value.length) {
      return;
    }
    $value.focus();
    if ($value[0] && typeof $value[0].select === 'function') {
      $value[0].select();
    }
  }, 0);
}

function mapperWheelRowHasReceptorContent($tr) {
  if (!$tr || !$tr.length) {
    return false;
  }
  var entry = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
  var typed = mapperWheelNormalizeTypedReceptorInput($tr.find('.mapper-core-in-receptor').val() || '');
  var unmatched = $tr.data('mapperCoreUnmatchedRaw');
  var unmatchedStr = unmatched != null ? String(unmatched).trim() : '';
  var viewText = ($tr.find('.mapper-core-receptor-html-view:visible').text() || '').replace(/\s+/g, ' ').trim();
  return !!(entry || typed || unmatchedStr || viewText);
}

function mapperWheelFocusReceptorCellForRow($tr) {
  if (!$tr || !$tr.length) {
    return;
  }
  window.setTimeout(function() {
    var $inp = $tr.find('.mapper-core-in-receptor');
    var $view = $tr.find('.mapper-core-receptor-html-view:visible');
    var $target = $inp.is(':visible') ? $inp : $view;
    if (!$target.length) {
      return;
    }
    $target.focus();
    if ($target[0] && typeof $target[0].select === 'function') {
      $target[0].select();
    }
  }, 0);
}

function mapperWheelFocusNextTableCellFrom($source) {
  var $tr = $source.closest('#mapper-wheel-input-tbody tr');
  if (!$tr.length) {
    return;
  }
  if ($source.closest('.mapper-core-in-receptor, .mapper-core-receptor-html-view').length) {
    mapperWheelFocusValueCellForRow($tr);
    return;
  }
  if (!$source.closest('.mapper-wheel-in-value').length) {
    return;
  }
  mapperWheelEnsureTrailingBlankRow();
  var $next = $tr.next('tr');
  if (!$next.length) {
    mapperWheelAppendRow('', '', true);
    $next = $('#mapper-wheel-input-tbody tr').last();
  }
  if (mapperWheelRowHasReceptorContent($next)) {
    mapperWheelFocusValueCellForRow($next);
  } else {
    mapperWheelFocusReceptorCellForRow($next);
  }
}

function mapperWheelBeginEditReceptor($tr) {
  var $hid = $tr.find('.mapper-core-receptor-entry');
  var id = ($hid.val() || '').trim();
  var $inp = $tr.find('.mapper-core-in-receptor');
  var $view = $tr.find('.mapper-core-receptor-html-view');
  var seed = mapperWheelReceptorEditSeed(id);
  $hid.val('');
  $view.hide().empty();
  $inp.val(seed).show().focus();
  window.setTimeout(function() {
    var t = ($inp.val() || '').trim();
    if (t.length >= 1 && $inp.data('ui-autocomplete')) {
      $inp.autocomplete('search', t);
    }
  }, 0);
  mapperWheelSyncReceptorClearBtn($tr);
  mapperWheelSyncReceptorTextareaHeights();
}

function mapperWheelDestroyReceptorAutocomplete($inp) {
  if ($inp && $inp.length && $inp.data('ui-autocomplete')) {
    try {
      $inp.autocomplete('destroy');
    } catch (e) {}
  }
}

function mapperWheelBindReceptorAutocomplete($inp) {
  if (!$inp || !$inp.length) {
    return;
  }
  mapperWheelDestroyReceptorAutocomplete($inp);
  $inp.autocomplete({
    minLength: 1,
    delay: 0,
    appendTo: '.mapper-core-wheel-left',
    source: function(request, response) {
      var items = mapperWheelFilterReceptorsLocal(request.term, 80);
      response($.map(items, function(it) {
        return {
          label: mapperWheelReceptorAutocompleteLabelPlain(it),
          value: it.id,
          html: mapperWheelResolvedReceptorDisplayHtml(it.id)
        };
      }));
    },
    select: function(event, ui) {
      var $tr = $inp.closest('tr');
      $tr.removeClass('mapper-core-row-invalid mapper-wheel-row-error');
      $tr.removeData('mapperCoreUnmatchedRaw');
      mapperWheelSetReceptorResolvedUI($tr, ui.item.value);
      mapperWheelEnsureTrailingBlankRow();
      mapperWheelFocusValueCellForRow($tr);
      return false;
    },
    open: function() {
      $(this).autocomplete('widget').addClass('mapper-core-receptor-ac-menu');
    },
    close: function() {
      var w = $(this).autocomplete('widget');
      if (w && w.length) {
        w.removeClass('mapper-core-receptor-ac-menu');
      }
    },
    create: function() {
      var widget = $(this).data('ui-autocomplete');
      widget._renderItem = function(ul, item) {
        var inner = item.html || $('<span/>').text(item.label || '').html();
        return $('<li>')
          .append($('<div class="mapper-core-ac-item-label">').html(inner))
          .appendTo(ul);
      };
    }
  });
  $inp.off('input.mapper-wheel-ta-grow').on('input.mapper-wheel-ta-grow', function() {
    if (mapperWheelReceptorTallPlaceholderIntro) {
      mapperWheelAutoSizeReceptorTextarea($inp);
    }
  });
  mapperWheelSyncReceptorTextareaHeights();
  $inp.off('focus.mapper-wheel-ac').on('focus.mapper-wheel-ac', function() {
    var t = ($inp.val() || '').trim();
    if (t.length >= 1) {
      $inp.autocomplete('search', t);
    }
  });
}

function mapperWheelEnsureReceptorInputCanonical($inp, id) {
  mapperWheelSetReceptorResolvedUI($inp.closest('tr'), id);
}

function mapperWheelDestroyRowReceptorWidgets() {
  $('#mapper-wheel-input-tbody .mapper-core-in-receptor').each(function() {
    mapperWheelDestroyReceptorAutocomplete($(this));
  });
  $('#mapper-wheel-input-tbody .mapper-core-row-color-picker').each(function() {
    mapperWheelDestroyRowColorSpectrum($(this));
  });
}

function mapperWheelDestroyRowReceptorWidgetForTr($tr) {
  mapperWheelDestroyReceptorAutocomplete($tr.find('.mapper-core-in-receptor'));
  mapperWheelDestroyRowColorSpectrum($tr.find('.mapper-core-row-color-picker'));
}

function mapperWheelApplyCanonicalToRowSelect($tr, entryId) {
  var id = entryId != null ? String(entryId).trim() : '';
  var $inp = $tr.find('.mapper-core-in-receptor');
  if (!$inp.length) {
    return;
  }
  $tr.removeClass('mapper-core-row-invalid mapper-wheel-row-error');
  $tr.removeData('mapperCoreUnmatchedRaw');
  mapperWheelEnsureReceptorInputCanonical($inp, id);
}

function mapperWheelUpgradeReceptorToSelect2($tr) {
}

function mapperWheelResolveEntry(raw) {
  if (raw == null || !String(raw).trim()) {
    return null;
  }
  var u = String(raw).trim().toUpperCase();
  return MAPPER_CORE_RESOLVE[u] || null;
}

function mapperWheelRowIsContentEmpty($tr) {
  var entry = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
  var v = mapperWheelNormalizeTypedReceptorInput($tr.find('.mapper-core-in-receptor').val() || '');
  var unmatched = $tr.data('mapperCoreUnmatchedRaw');
  var um = unmatched != null ? String(unmatched).trim() : '';
  var valCell = ($tr.find('.mapper-wheel-in-value').val() || '').trim();
  return !entry && !v && !um && !valCell;
}

function mapperWheelSyncFirstRowReceptorPlaceholder() {
  var hint = MAPPER_WHEEL_FIRST_ROW_RECEPTOR_PLACEHOLDER;
  $('#mapper-wheel-input-tbody tr').each(function(idx) {
    var $inp = $(this).find('.mapper-core-in-receptor');
    if (!$inp.length) {
      return;
    }
    if (idx === 0) {
      $inp.attr('placeholder', hint);
    } else {
      $inp.removeAttr('placeholder');
    }
  });
}

// Batches the two full-table cosmetic syncs (placeholder text, remove-button visibility)
// so they run once per debounce window instead of once per keystroke/click when
// mapperWheelEnsureTrailingBlankRow(true) is called from the hottest, highest-frequency
// actions (typing in a value cell, deleting a row) at large row counts.
var _wheelCosmeticSyncDebounced = MapperPageCore.debounce(function () {
  mapperWheelSyncFirstRowReceptorPlaceholder();
  mapperWheelSyncRemoveRowButtons();
}, 200);

function mapperWheelEnsureTrailingBlankRow(deferCosmetics) {
  var $tbody = $('#mapper-wheel-input-tbody');

  function rows() {
    return $tbody.children('tr');
  }

  var $rows = rows();
  if (!$rows.length) {
    mapperWheelAppendRow('', '', true);
  } else {
    while ($rows.length >= 2) {
      var $last = $rows.last();
      var $prev = $last.prev();
      if (mapperWheelRowIsContentEmpty($last) && mapperWheelRowIsContentEmpty($prev)) {
        mapperWheelDestroyRowReceptorWidgetForTr($last);
        $last.remove();
        $rows = rows();
        continue;
      }
      break;
    }
    $rows = rows();
    var $lastOne = $rows.last();
    if (!mapperWheelRowIsContentEmpty($lastOne) && $rows.length < MAPPER_WHEEL_MAX_ROWS) {
      mapperWheelAppendRow('', '', true);
    }
  }

  mapperWheelSyncReceptorTextareaHeights();
  if (deferCosmetics) {
    _wheelCosmeticSyncDebounced.schedule();
  } else {
    mapperWheelSyncFirstRowReceptorPlaceholder();
    mapperWheelSyncRemoveRowButtons();
  }
  if (!mapperWheelRestoringDomRows) {
    mapperWheelScheduleRebuildPlot();
  }
}

function mapperWheelSyncRemoveRowButtons() {
  $('#mapper-wheel-input-tbody tr').each(function() {
    var $tr = $(this);
    var hide = mapperWheelRowIsContentEmpty($tr);
    $tr.find('td.mapper-core-remove-cell').toggleClass('is-remove-hidden', hide);
  });
}

function mapperWheelAppendRow(receptorText, value, skipTrailingEnsure) {
  if (skipTrailingEnsure === undefined) {
    skipTrailingEnsure = false;
  }
  var tr = $('<tr>');
  tr.append(
    $('<td class="mapper-core-remove-cell is-remove-hidden">').append(
      $('<button type="button" class="mapper-core-remove-row" aria-label="Remove row" title="Remove row">&times;</button>')
    )
  );
  var $tdR = $('<td class="mapper-core-receptor-cell">');
  var $inp = mapperWheelCreateReceptorCell($tdR);
  tr.append($tdR);
  tr.append($('<td class="mapper-core-value-cell">').append($('<input type="text" class="form-control input-sm mapper-wheel-in-value" autocomplete="off">').val(value != null ? value : '')));
  tr.append($('<td class="mapper-wheel-swatch-cell">').append('<input type="text" class="mapper-wheel-label-swatch mapper-core-row-color-picker" aria-label="Label colour">'));
  $('#mapper-wheel-input-tbody').append(tr);
  mapperWheelBindReceptorAutocomplete($inp);
  var rt = receptorText != null ? String(receptorText).trim() : '';
  var $tr = tr;
  if (rt) {
    var resolved = mapperWheelResolveEntry(rt);
    var knownId = receptorSelect2Data.some(function(x) {
      return String(x.id) === rt;
    });
    if (resolved || knownId) {
      var rid = resolved || rt;
      mapperWheelEnsureReceptorInputCanonical($inp, rid);
      $tr.removeClass('mapper-core-row-invalid');
      $tr.removeData('mapperCoreUnmatchedRaw');
    } else {
      $tr.find('.mapper-core-receptor-entry').val('');
      $tr.find('.mapper-core-receptor-html-view').hide().empty();
      $inp.val(rt).show();
      $tr.addClass('mapper-core-row-invalid');
      $tr.data('mapperCoreUnmatchedRaw', rt);
    }
  }
  mapperWheelSyncReceptorClearBtn($tr);
  mapperWheelSyncRowSwatch($tr);
  mapperWheelSyncValueCellHint($tr);
  mapperWheelSyncRemoveRowButtons();
  if (!skipTrailingEnsure) {
    mapperWheelEnsureTrailingBlankRow();
  }
}

function mapperWheelAppendRowFromPaste(rawR, valCell) {
  mapperWheelAppendRow('', valCell, true);
  var $tr = $('#mapper-wheel-input-tbody tr').last();
  mapperWheelSetReceptorRawOnRow($tr, rawR);
  mapperWheelEnsureTrailingBlankRow();
}

function mapperWheelFillDemoRows() {
  mapperWheelReceptorTallPlaceholderIntro = false;
  mapperWheelDestroyRowReceptorWidgets();
  $('#mapper-wheel-input-tbody').empty();
  $('#mapper-core-wheel-errors').empty();
  if (mapperWheelIsTextMode()) {
    MAPPER_CORE_LABEL_COLORS = {};
    MAPPER_CORE_LABEL_ENABLED = {};
  }
  Object.keys(MAPPER_WHEEL_DEMO_ROWS).forEach(function(receptor) {
    var row = MAPPER_WHEEL_DEMO_ROWS[receptor] || {};
    var value = mapperWheelIsTextMode() ? row.text : row.numeric;
    mapperWheelAppendRow(receptor, value, true);
  });
  mapperWheelEnsureTrailingBlankRow();
  mapperWheelRefreshAllLabelSwatches();
  mapperWheelRebuildLabelColorCustomizePanel(false);
  mapperWheelCaptureSnapshotForActiveMode();
  mapperWheelRebuildPlotImmediate();
}

function mapperWheelSetValueCellRawOnRow($tr, valCell) {
  if (!$tr || !$tr.length) {
    return;
  }
  $tr.find('.mapper-wheel-in-value').val(valCell != null ? String(valCell) : '');
}

function mapperWheelSetValueOnRow($tr, valCell) {
  mapperWheelSetValueCellRawOnRow($tr, valCell);
  mapperWheelSyncRowSwatch($tr);
  mapperWheelSyncValueCellHint($tr);
}

function mapperWheelSetReceptorRawOnRow($tr, rawR) {
  if (!$tr || !$tr.length) {
    return;
  }
  rawR = rawR != null ? String(rawR).trim() : '';
  var resolved = rawR ? mapperWheelResolveEntry(rawR) : null;
  var knownId = rawR && receptorSelect2Data.some(function(x) {
    return String(x.id) === rawR;
  });
  var $td = $tr.find('td.mapper-core-receptor-cell');
  mapperWheelDestroyReceptorAutocomplete($tr.find('.mapper-core-in-receptor'));
  $td.empty();
  var $inp = mapperWheelCreateReceptorCell($td);
  mapperWheelBindReceptorAutocomplete($inp);
  $tr.removeClass('mapper-wheel-row-error');
  if (!rawR) {
    $tr.removeClass('mapper-core-row-invalid');
    $tr.removeData('mapperCoreUnmatchedRaw');
    $inp.show();
    mapperWheelSyncReceptorClearBtn($tr);
    return;
  }
  mapperWheelCompactReceptorRowsAfterInput();
  if (resolved || knownId) {
    mapperWheelEnsureReceptorInputCanonical($inp, resolved || rawR);
    $tr.removeClass('mapper-core-row-invalid');
    $tr.removeData('mapperCoreUnmatchedRaw');
  } else {
    $td.find('.mapper-core-receptor-entry').val('');
    $td.find('.mapper-core-receptor-html-view').hide().empty();
    $inp.val(rawR).show();
    $tr.addClass('mapper-core-row-invalid');
    $tr.data('mapperCoreUnmatchedRaw', rawR);
  }
  mapperWheelSyncReceptorClearBtn($tr);
}

function mapperWheelPasteAnchorFromEvent(e) {
  var $t = $(e.target);
  var $tr = $t.closest('#mapper-wheel-input-tbody tr');
  if ($t.closest('.mapper-wheel-in-value').length && $tr.length) {
    return { $tr: $tr, col: 'value' };
  }
  if ($tr.length && ($t.closest('.mapper-core-in-receptor').length || $t.closest('.mapper-core-receptor-html-view').length || $t.closest('.mapper-core-receptor-clear').length || $t.closest('.mapper-core-receptor-input-wrap').length || $t.closest('td.mapper-core-receptor-cell').length)) {
    return { $tr: $tr, col: 'receptor' };
  }
  return { $tr: $tr, col: null };
}

function mapperWheelNormalizePastedLines(text) {
  var lines = String(text || '').split(/\r?\n/);
  while (lines.length && lines[lines.length - 1] === '') {
    lines.pop();
  }
  return lines;
}

/** Prefer TAB columns (spreadsheet); else first `;` splits receptor/value (CSV / Excel semicolon export). */
function mapperWheelSplitPasteLineIntoReceptorAndValue(line) {
  var s = String(line || '').trim();
  if (!s) {
    return { rawR: '', valCell: '' };
  }
  var ti = s.indexOf('\t');
  if (ti !== -1) {
    var p = s.split('\t');
    var rawRt = p[0] ? p[0].trim() : '';
    var rest = p.length > 1 ? p.slice(1).join('\t').trim() : '';
    return { rawR: rawRt, valCell: rest };
  }
  var si = s.indexOf(';');
  if (si !== -1) {
    return {
      rawR: s.slice(0, si).trim(),
      valCell: s.slice(si + 1).trim()
    };
  }
  return { rawR: s.trim(), valCell: '' };
}

/** When clipboard has no TAB, treat nonempty lines that all contain ';' as receptor;value rows (incl. one-line CSV). */
function mapperWheelPasteLooksSemicolonCsvGrid(contentLines, hasTab) {
  if (hasTab || !contentLines || !contentLines.length) {
    return false;
  }
  return contentLines.every(function(line) {
    return String(line).trim().indexOf(';') !== -1;
  });
}

/** Align with Classification wheel: rename a few receptor-family keys and refresh leaf metadata. */
function mapperWheelClassicWheelNormGroupKey(fam) {
  var s = String(fam == null ? '' : fam).trim();
  if (!s) return s;
  var lower = s.toLowerCase();
  var familyRename = {
    'opsins': 'Vision receptors',
    'opsins receptors': 'Vision receptors',
    'class a orphans': 'Orphans',
    'class c orphans': 'Orphans',
    'class c orphans ': 'Orphans',
    'peptide p518 receptors': 'QRFP receptors'
  };
  return familyRename[lower] || s;
}

function mapperWheelPreprocessClassicWheelData(data) {
  if (!data || typeof data !== 'object') {
    return data;
  }
  Object.keys(data).forEach(function(circleKey) {
    var circle = data[circleKey];
    if (!circle || typeof circle !== 'object') {
      return;
    }
    Object.keys(circle).forEach(function(classKey) {
      var fams = circle[classKey];
      if (!fams || typeof fams !== 'object') {
        return;
      }
      var newFams = {};
      Object.keys(fams).forEach(function(oldFamKey) {
        var newFamKey = mapperWheelClassicWheelNormGroupKey(oldFamKey);
        var receptors = fams[oldFamKey];
        if (!receptors || typeof receptors !== 'object') {
          return;
        }
        Object.keys(receptors).forEach(function(recName) {
          var r = receptors[recName];
          if (r && typeof r === 'object' && Object.prototype.hasOwnProperty.call(r, 'Receptor family')) {
            r['Receptor family'] = newFamKey;
          }
        });
        if (!newFams[newFamKey]) {
          newFams[newFamKey] = {};
        }
        Object.keys(receptors).forEach(function(rk) {
          newFams[newFamKey][rk] = receptors[rk];
        });
      });
      circle[classKey] = newFams;
    });
  });
  return data;
}

function mapperWheelWheelBadgeNormKey(code) {
  if (code === undefined || code === null) {
    return '';
  }
  var s = String(code).trim();
  if (!s || s.toLowerCase() === 'nan') {
    return '';
  }
  return s;
}

function mapperWheelWheelClassDisplayShortLabel(code) {
  var k = mapperWheelWheelBadgeNormKey(code);
  if (!k) {
    return '';
  }
  if (k === 'Unclassified') {
    return 'U';
  }
  if (/^(A|B1|B2|C|F|T2|V)$/i.test(k)) {
    return k.toUpperCase();
  }
  return k;
}

function mapperWheelWheelUiD3() {
  return typeof d3v4 !== 'undefined' ? d3v4 : d3;
}

function mapperWheelMapperClassBadgeFill(classCode) {
  /* Class badges stay neutral; Numeric/Text data colours should be the only coloured accents. */
  return '#ffffff';
}

function mapperWheelAddClassBadgePillsAfterDraw(locationId) {
  var d3pick = mapperWheelWheelUiD3();
  var svgSel = d3pick.select('#' + locationId + '_svg');
  if (!svgSel || svgSel.empty()) {
    return;
  }
  var UNCLASSIFIED_BADGE_DX = -15;

  svgSel.selectAll('text')
    .filter(function() {
      var el = this;
      var cls = (el.getAttribute && el.getAttribute('class')) ? el.getAttribute('class') : '';
      /* datamapper: class ring labels use `GPCRome-text-{level}-highlight` only */
      return /GPCRome-text-\d+-highlight/.test(cls) && cls.indexOf('GPCRome-family-label') === -1;
    })
    .each(function(d) {
      try {
        var txt = d3pick.select(this);
        var rawClass = mapperWheelWheelBadgeNormKey(d);
        if (!rawClass) {
          rawClass = mapperWheelWheelBadgeNormKey(this.textContent);
        }
        if (!rawClass) {
          return;
        }
        var classKey = mapperWheelWheelBadgeNormKey(rawClass);
        txt.text(mapperWheelWheelClassDisplayShortLabel(rawClass));
        txt.style('fill', '#000');

        var fill = mapperWheelMapperClassBadgeFill(rawClass);
        var fillOpacity = 1;
        var node = txt.node();
        if (!node) {
          return;
        }
        var bb = node.getBBox();
        var padX = 3;
        var padY = 1;
        var g = node.parentNode;
        if (!g || !g.insertBefore) {
          return;
        }
        var rectNode = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
        rectNode.setAttribute('x', String(bb.x - padX));
        rectNode.setAttribute('y', String(bb.y - padY));
        rectNode.setAttribute('width', String(bb.width + padX * 2));
        rectNode.setAttribute('height', String(bb.height + padY * 2));
        rectNode.setAttribute('rx', '9');
        rectNode.setAttribute('ry', '9');
        rectNode.setAttribute('fill', fill);
        rectNode.setAttribute('fill-opacity', String(fillOpacity));
        rectNode.setAttribute('stroke', '#000');
        rectNode.setAttribute('stroke-width', '0.75px');
        g.insertBefore(rectNode, node);

        if (classKey === 'Unclassified') {
          var tPrev = txt.attr('transform') || '';
          txt.attr('transform', (tPrev ? (tPrev + ' ') : '') + 'translate(' + UNCLASSIFIED_BADGE_DX + ',0)');
          rectNode.setAttribute('transform', 'translate(' + UNCLASSIFIED_BADGE_DX + ',0)');
        }
      } catch (e) {
        /* ignore pill layout errors */
      }
    });
}

function mapperWheelDrawGPCRomeWithClassBadges(data, locationId, styling) {
  DrawGPCRomeWheel(data, locationId, styling);
  mapperWheelAddClassBadgePillsAfterDraw(locationId);
}

mapperWheelPreprocessClassicWheelData(GPCRome_WheelDict);
/** Deep copy of empty wheel JSON; Apply clones this and merges user values without re-fetching. */
var MAPPER_WHEEL_GPCROME_BASE = JSON.parse(JSON.stringify(GPCRome_WheelDict));
var GPCRome_location = "GPCRome_plot";
var GPCRomes_styling = {};

function mapperWheelNormalizeGpcromeMergeKey(rawKey) {
  var s = String(rawKey || '');
  var i = s.indexOf('_');
  return (i >= 0 ? s.slice(0, i) : s).toUpperCase();
}

/** Mirrors DataMapperHome.update_nested_GPCRome_data (mapper/views.py). Mutates structureDict. */
function mapperWheelUpdateNestedGPCRomeData(structureDict, rawData) {
  var normalizedRawData = {};
  Object.keys(rawData || {}).forEach(function(key) {
    var val = rawData[key];
    if (val && typeof val === 'object' && Object.prototype.hasOwnProperty.call(val, 'Value1')) {
      normalizedRawData[mapperWheelNormalizeGpcromeMergeKey(key)] = {
        Data: val.Value1,
        Color: val.Value2
      };
    }
  });
  function recursiveUpdate(d) {
    if (d && typeof d === 'object') {
      if (!Array.isArray(d)) {
        if ('EntryName' in d && 'Data' in d) {
          var entry = String(d.EntryName).toUpperCase();
          if (Object.prototype.hasOwnProperty.call(normalizedRawData, entry)) {
            d.Data = normalizedRawData[entry].Data;
            var color = normalizedRawData[entry].Color;
            if (color != null) {
              d.Color = color;
            }
          }
        }
        Object.keys(d).forEach(function(k) {
          recursiveUpdate(d[k]);
        });
      } else {
        d.forEach(function(item) {
          recursiveUpdate(item);
        });
      }
    }
  }
  recursiveUpdate(structureDict);
  return structureDict;
}

function mapperWheelNormalizeReceptorForWheel(rs) {
  var receptor = mapperWheelResolveEntry(rs);
  if (!receptor || !MAPPER_WHEEL_TREE_IDS[receptor]) {
    return null;
  }
  return receptor;
}

function extractDataValues(obj) {
  var values = [];
  function findData(node) {
    if (typeof node === "object" && node !== null) {
      for (var key in node) {
        if (!Object.prototype.hasOwnProperty.call(node, key)) continue;
        var value = node[key];
        if (value && typeof value === "object" && "Data" in value) {
          var dataValue = Number(value.Data);
          if (!isNaN(dataValue)) {
            values.push(dataValue);
          }
        } else {
          findData(value);
        }
      }
    }
  }
  findData(obj);
  return values;
}

function determineDecimalPlaces(minValue, maxValue) {
  var allValues = [minValue, maxValue];
  var maxDecimals = 0;
  allValues.forEach(function(value) {
    if (value == null) return;
    var valueStr = Math.abs(value).toString();
    if (valueStr.indexOf('e') !== -1) {
      var parts = valueStr.split('e-');
      if (parts.length === 2) {
        var exp = parseInt(parts[1], 10);
        maxDecimals = Math.max(maxDecimals, exp);
      }
    } else if (valueStr.indexOf('.') !== -1) {
      var decimals = valueStr.split('.')[1].length;
      maxDecimals = Math.max(maxDecimals, decimals);
    }
  });
  return Math.min(maxDecimals, 6);
}

function mapperWheelRecalcNumericStyling() {
  if (mapperWheelIsTextMode()) {
    return;
  }
  var GPCRomeWheelValues = extractDataValues(GPCRome_WheelDict)
    .filter(function(value) { return value !== null && value !== ''; });
  var GPCRomeMin = GPCRomeWheelValues.length > 0 ? Math.min.apply(null, GPCRomeWheelValues) : null;
  var GPCRomeMax = GPCRomeWheelValues.length > 0 ? Math.max.apply(null, GPCRomeWheelValues) : null;
  var GPCRomeAvg = GPCRomeWheelValues.length > 0 ? (GPCRomeMin + GPCRomeMax) / 2 : null;
  GPCRomes_styling.GPCRomeMin = GPCRomeMin;
  GPCRomes_styling.GPCRomeMax = GPCRomeMax;
  GPCRomes_styling.GPCRomeAvg = GPCRomeAvg;
  GPCRomes_styling.LegendbarDigit = determineDecimalPlaces(GPCRomeMin, GPCRomeMax);
}

function mapperWheelCollectDisplayLabelsForNumericValues() {
  var values = {};
  var labelType = GPCRomes_styling && GPCRomes_styling.LabelType ? GPCRomes_styling.LabelType : 'Protein';
  function walk(node, keyName) {
    if (!node || typeof node !== 'object') {
      return;
    }
    if (!Array.isArray(node) && Object.prototype.hasOwnProperty.call(node, 'Data')) {
      var num = Number(node.Data);
      if (!isNaN(num)) {
        var label;
        if (labelType === 'Uniprot') {
          label = String(node.EntryName || keyName || '');
        } else if (labelType === 'Entrez') {
          label = String(node.Entrez || keyName || '');
        } else {
          label = String(keyName || '');
        }
        values[label] = num;
      }
    }
    Object.keys(node).forEach(function(k) {
      walk(node[k], k);
    });
  }
  walk(GPCRome_WheelDict, '');
  return values;
}

function mapperWheelNumericColorScaleForCurrentStyling() {
  var minValue = GPCRomes_styling.GPCRomeMin || 0;
  var maxValue = GPCRomes_styling.GPCRomeMax || 1;
  var avgValue = GPCRomes_styling.GPCRomeAvg || 0.5;
  var colorMin = GPCRomes_styling.colorStart || '#FFFFFF';
  var colorMax = GPCRomes_styling.colorEnd || '#000000';
  var colorAvg = '#FFFFFF';
  if (minValue === maxValue) {
    var eps = 1e-9 * (Math.abs(minValue) || 1);
    minValue -= eps;
    maxValue += eps;
  }
  if (GPCRomes_styling.ColorSetup === 'Two') {
    return d3v4.scaleLinear().domain([minValue, maxValue]).range([colorMin, colorMax]);
  }
  if (GPCRomes_styling.ColorSetup === 'Three') {
    return d3v4.scaleLinear().domain([minValue, avgValue, maxValue]).range([colorMin, colorAvg, colorMax]);
  }
  return d3v4.scaleLinear().domain([minValue, maxValue]).range([colorAvg, colorMax]);
}

function mapperWheelUpdateNumericLegendGradientInPlace() {
  var gradient = d3v4.select('#gradient-bar-' + GPCRome_location);
  if (!gradient || gradient.empty()) {
    return;
  }
  var colorMin = GPCRomes_styling.colorStart || '#FFFFFF';
  var colorMax = GPCRomes_styling.colorEnd || '#000000';
  gradient.selectAll('stop').remove();
  if (GPCRomes_styling.ColorSetup === 'Three') {
    gradient.append('stop').attr('offset', '0%').attr('stop-color', colorMin);
    gradient.append('stop').attr('offset', '50%').attr('stop-color', '#FFFFFF');
    gradient.append('stop').attr('offset', '100%').attr('stop-color', colorMax);
  } else if (GPCRomes_styling.ColorSetup === 'Two') {
    gradient.append('stop').attr('offset', '0%').attr('stop-color', colorMin);
    gradient.append('stop').attr('offset', '100%').attr('stop-color', colorMax);
  } else {
    gradient.append('stop').attr('offset', '0%').attr('stop-color', '#FFFFFF');
    gradient.append('stop').attr('offset', '100%').attr('stop-color', colorMax);
  }
}

function mapperWheelRecolorExistingWheelNumericValues() {
  if (mapperWheelIsTextMode()) {
    return false;
  }
  var svg = d3v4.select('#' + GPCRome_location + '_svg');
  if (!svg || svg.empty()) {
    return false;
  }
  var valuesByLabel = mapperWheelCollectDisplayLabelsForNumericValues();
  var colorScale = mapperWheelNumericColorScaleForCurrentStyling();
  var changed = false;
  svg.selectAll('path[class^="large-hollow-pie-"], path[class*=" large-hollow-pie-"]')
    .style('fill', function(d) {
      if (!d || !Object.prototype.hasOwnProperty.call(valuesByLabel, String(d.data || ''))) {
        return d3v4.select(this).style('fill') || '#ffffff';
      }
      changed = true;
      return colorScale(valuesByLabel[String(d.data || '')]);
    })
    .style('stroke', function(d) {
      if (!d || !Object.prototype.hasOwnProperty.call(valuesByLabel, String(d.data || ''))) {
        return d3v4.select(this).style('stroke') || 'black';
      }
      return 'black';
    });
  mapperWheelUpdateNumericLegendGradientInPlace();
  return changed;
}

function updateGPCRome() {
  d3.select("#" + GPCRome_location).select("svg").remove();
  mapperWheelDrawGPCRomeWithClassBadges(GPCRome_WheelDict, GPCRome_location, GPCRomes_styling);
}

if (true) {
  GPCRomes_styling = {
    GPCRomeMin: null,
    GPCRomeMax: null,
    GPCRomeAvg: null,
    LegendbarDigit: 2,
    LegendbarLength: 200,
    LegendbarFontsize: "11px",
    FontStyle: "Arial",
    FontsizeGlobal: "11px",
    FontsizeClass: "20px",
    DataType: "Numeric",
    colorStart: "#e64b35",
    colorAvg: "#FFFFFF",
    colorEnd: "#3c5488",
    ColorSetup: 'One',
    showIcon: true,
    LabelType: "Protein",
    ShowLegend: true,
    LegendLayout: { mode: 'row', columns: '1', sorted: 'Vertically' }
  };
  mapperWheelRecalcNumericStyling();
  mapperWheelDrawGPCRomeWithClassBadges(GPCRome_WheelDict, GPCRome_location, GPCRomes_styling);
}

document.addEventListener('DOMContentLoaded', function () {
  var labelButtons = document.querySelectorAll('.GPCRome-label-btn');
  labelButtons.forEach(function(btn) {
    if (btn.getAttribute('data-value') === GPCRomes_styling.LabelType) {
      btn.classList.remove('btn-outline-primary');
      btn.classList.add('btn-primary');
    }
    btn.addEventListener('click', function () {
      var selectedValue = this.getAttribute('data-value');
      GPCRomes_styling.LabelType = selectedValue;
      updateGPCRome();
      $('#mapper-wheel-input-tbody .mapper-core-in-receptor').each(function() {
        var $ta = $(this);
        try {
          if ($ta.autocomplete && $ta.data('ui-autocomplete')) {
            $ta.autocomplete('close');
          }
        } catch (e) {}
      });
      mapperWheelRefreshResolvedReceptorDisplaysInGrid();
      labelButtons.forEach(function(b) {
        b.classList.remove('btn-primary');
        b.classList.add('btn-outline-primary');
      });
      this.classList.remove('btn-outline-primary');
      this.classList.add('btn-primary');
    });
  });
});

document.addEventListener('DOMContentLoaded', function () {
  var graphicalLegendOptions = document.querySelectorAll('.graphical-legend-option');
  graphicalLegendOptions.forEach(function(btn) {
    var iconValue = btn.getAttribute('data-icon') === 'true';
    var legendValue = btn.getAttribute('data-legend') === 'true';
    if (iconValue === GPCRomes_styling.showIcon && legendValue === GPCRomes_styling.ShowLegend) {
      btn.classList.remove('btn-outline-primary');
      btn.classList.add('btn-primary');
    }
    btn.addEventListener('click', function () {
      GPCRomes_styling.showIcon = iconValue;
      GPCRomes_styling.ShowLegend = legendValue;
      updateGPCRome();
      graphicalLegendOptions.forEach(function(b) {
        b.classList.remove('btn-primary');
        b.classList.add('btn-outline-primary');
      });
      this.classList.remove('btn-outline-primary');
      this.classList.add('btn-primary');
    });
  });
});

document.addEventListener('DOMContentLoaded', function () {
  var downloadBtn = document.querySelector('#DownloadDropdownToggle');
  var downloadMenu = downloadBtn ? downloadBtn.nextElementSibling : null;
  if (downloadBtn && downloadMenu) {
    downloadBtn.addEventListener('click', function (e) {
      e.stopPropagation();
      downloadMenu.style.display = (downloadMenu.style.display === 'block') ? 'none' : 'block';
    });
    document.addEventListener('click', function () {
      downloadMenu.style.display = 'none';
    });
    downloadMenu.addEventListener('click', function (e) {
      e.stopPropagation();
    });
  }
});

if (true) {
  function UpdateGPCRomeColorPickers() {
    if (mapperWheelIsTextMode()) {
      return;
    }
    var selectedValue = $('#GPCRome_color_styling').val();
    var colorPresets = {
      "One": { setup: "One", colorStart: "#ffffff", colorAvg: "#ffffff", colorEnd: "#707070" },
      "Two": { setup: "Two", colorStart: "#97a6c4", colorAvg: "#ffffff", colorEnd: "#384860" },
      "Three_RWB": { setup: "Three", colorStart: "#a00000", colorAvg: "#ffffff", colorEnd: "#1a80bb" },
      "Three_TWM": { setup: "Three", colorStart: "#298c8c", colorAvg: "#ffffff", colorEnd: "#800074" }
    };
    var preset = colorPresets[selectedValue];
    if (preset) {
      GPCRomes_styling.ColorSetup = preset.setup;
      GPCRomes_styling.colorStart = preset.colorStart;
      GPCRomes_styling.colorAvg = preset.colorAvg;
      GPCRomes_styling.colorEnd = preset.colorEnd;
      try {
        $("#GPCRome_colorPicker_min").spectrum("set", preset.colorStart);
        $("#GPCRome_colorPicker_avg").spectrum("set", preset.colorAvg);
        $("#GPCRome_colorPicker_max").spectrum("set", preset.colorEnd);
      } catch (eCol) {}
    }
    var showMin = GPCRomes_styling.ColorSetup !== 'One';
    var showAvg = GPCRomes_styling.ColorSetup === 'Three';
    document.getElementById("GPCRome_Color_min_container").style.visibility = showMin ? 'visible' : 'hidden';
    document.getElementById("GPCRome_Color_avg_container").style.visibility = showAvg ? 'visible' : 'hidden';
    document.getElementById("GPCRome_Color_max_container").style.visibility = 'visible';
    document.getElementById("GPCRome_Color_min_label").style.visibility = showMin ? 'visible' : 'hidden';
    document.getElementById("GPCRome_Color_avg_label").style.visibility = showAvg ? 'visible' : 'hidden';
    if (!mapperWheelRecolorExistingWheelNumericValues()) {
      updateGPCRome();
    }
  }

  function updateGPCRomeVisualization() {
    if (mapperWheelIsTextMode()) {
      return;
    }
    var colorMin = $("#GPCRome_colorPicker_min").spectrum("get").toHexString();
    var colorAvg = $("#GPCRome_colorPicker_avg").spectrum("get").toHexString();
    var colorMax = $("#GPCRome_colorPicker_max").spectrum("get").toHexString();
    GPCRomes_styling.colorStart = colorMin;
    GPCRomes_styling.colorAvg = colorAvg;
    GPCRomes_styling.colorEnd = colorMax;
    if (!mapperWheelRecolorExistingWheelNumericValues()) {
      updateGPCRome();
    }
  }

  $(document).ready(function() {
    function initializeGPCRomeColorPicker(elementId, start_color) {
      $("#" + elementId).spectrum({
        color: start_color,
        showPalette: true,
        showInput: true,
        preferredFormat: "hex",
        palette: [
          ["#000", "#FF0000", "#00FF00", "#0000FF", "#FFFF00"],
          ["#FF00FF", "#00FFFF", "#FFFFFF", "#C0C0C0", "#808080"],
          ["#800000", "#808000", "#008000", "#800080", "#008080"],
          ["#000080"]
        ],
        change: function() { updateGPCRomeVisualization(); },
        move: function() { updateGPCRomeVisualization(); }
      });
    }
    initializeGPCRomeColorPicker("GPCRome_colorPicker_min", GPCRomes_styling.colorStart);
    initializeGPCRomeColorPicker("GPCRome_colorPicker_avg", GPCRomes_styling.colorAvg);
    $("#GPCRome_colorPicker_avg").spectrum("disable");
    initializeGPCRomeColorPicker("GPCRome_colorPicker_max", GPCRomes_styling.colorEnd);
  });

  $(document).ready(function () {
    function formatColorScheme(option) {
      if (!option.id) return option.text;
      var colorSlots = {
        "One": [null, null, "#707070"],
        "Two": ["#97a6c4", null, "#384860"],
        "Three_RWB": ["#a00000", "#ffffff", "#1a80bb"],
        "Three_TWM": ["#298c8c", "#ffffff", "#800074"]
      };
      var colors = colorSlots[option.id];
      if (colors) {
        return $('<span style="display:flex; align-items:center;"><span style="min-width: 50px; display:inline-block;text-align: center; padding-right: 5px;">' + option.text + '</span>' +
          colors.map(function(color) {
            var bg = color || 'transparent';
            var opacity = color ? '' : 'opacity:0;';
            return '<span style="display:inline-block;width:14px;height:14px;margin-left:5px;border:1px solid #ccc;border-radius:2px;background:' + bg + ';' + opacity + '"></span>';
          }).join('') + '</span>');
      }
      return option.text;
    }
    $('#GPCRome_color_styling').select2({
      templateResult: formatColorScheme,
      templateSelection: formatColorScheme,
      width: 'resolve'
    });
    $('#GPCRome_color_styling').on('change', function () {
      UpdateGPCRomeColorPickers();
    });
    UpdateGPCRomeColorPickers();
  });
  window.Mapper20UpdateNumericColorPanels = UpdateGPCRomeColorPickers;
}

function mapperWheelScheduleLabelColorPanelRebuildSoon() {
  if (!mapperWheelIsTextMode()) {
    return;
  }
  window.clearTimeout(mapperWheelLabelPanelDebouncer);
  mapperWheelLabelPanelDebouncer = window.setTimeout(function() {
    mapperWheelLabelPanelDebouncer = null;
    mapperWheelRebuildLabelColorCustomizePanel(false);
  }, 260);
}

function mapperWheelCollectRows() {
  var rows = [];
  $('#mapper-wheel-input-tbody tr').each(function() {
    var $tr = $(this);
    var r = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
    if (!r) {
      var $inp = $tr.find('.mapper-core-in-receptor');
      r = mapperWheelNormalizeTypedReceptorInput($inp.val() || '');
    }
    if (!r) {
      r = ($tr.data('mapperCoreUnmatchedRaw') || '').trim();
    }
    var v = $tr.find('.mapper-wheel-in-value').val();
    rows.push({ receptor_raw: r, value_raw: v });
  });
  return rows;
}

$(document).ready(function() {
  mapperWheelInitLegendLayoutSelectUi();
  mapperWheelAppendRow('', '');
  mapperWheelApplyChromeForInputMode();
  mapperWheelSyncLegendLayoutUiFromStyling();
  mapperWheelCaptureSnapshotForActiveMode();

  $('#mapper-wheel-mode-numeric-btn').on('click', function(e) {
    e.preventDefault();
    mapperWheelSetWheelInputMode('numeric');
  });
  $('#mapper-wheel-mode-labels-btn').on('click', function(e) {
    e.preventDefault();
    mapperWheelSetWheelInputMode('text');
  });

  $('#mapper-wheel-legend-layout-select').on('change', function() {
    if (!mapperWheelIsTextMode()) {
      return;
    }
    mapperWheelApplyLegendLayoutSelectionFromDropdown();
    updateGPCRome();
  });

  $('#mapper-wheel-toggle-legend-sorting').on('click', function() {
    if (!mapperWheelIsTextMode()) {
      return;
    }
    if (!GPCRomes_styling.LegendLayout) {
      GPCRomes_styling.LegendLayout = {};
    }
    var next = (GPCRomes_styling.LegendLayout.sorted === 'Horizontally') ? 'Vertically' : 'Horizontally';
    GPCRomes_styling.LegendLayout.sorted = next;
    $(this).text(next === 'Horizontally' ? 'Horizontally' : 'Vertically');
    updateGPCRome();
  });

  $('#mapper-wheel-input-tbody').on('input', '.mapper-core-in-receptor', function() {
    var $tr = $(this).closest('tr');
    if (mapperWheelNormalizeTypedReceptorInput($(this).val() || '')) {
      mapperWheelCompactReceptorRowsAfterInput();
    }
    $tr.removeClass('mapper-core-row-invalid mapper-wheel-row-error');
    $tr.removeData('mapperCoreUnmatchedRaw');
    mapperWheelClearCellHoverHint($tr.find('td.mapper-core-receptor-cell'));
    mapperWheelSyncReceptorClearBtn($tr);
    mapperWheelEnsureTrailingBlankRow();
  });

  $('#mapper-wheel-input-tbody').on('click', '.mapper-core-receptor-clear', function(e) {
    e.preventDefault();
    e.stopPropagation();
    var $tr = $(this).closest('tr');
    var $inp = $tr.find('.mapper-core-in-receptor');
    mapperWheelDestroyReceptorAutocomplete($inp);
    $tr.find('.mapper-core-receptor-entry').val('');
    $tr.find('.mapper-core-receptor-html-view').hide().empty();
    $tr.removeClass('mapper-core-row-invalid mapper-wheel-row-error');
    $tr.removeData('mapperCoreUnmatchedRaw');
    mapperWheelClearCellHoverHint($tr.find('td.mapper-core-receptor-cell'));
    $inp.val('').show().focus();
    mapperWheelBindReceptorAutocomplete($inp);
    mapperWheelSyncReceptorClearBtn($tr);
    mapperWheelEnsureTrailingBlankRow();
  });

  $('#mapper-wheel-input-tbody').on('click', '.mapper-core-remove-row', function(e) {
    e.preventDefault();
    e.stopPropagation();
    var $tr = $(this).closest('tr');
    if (!$tr.length || mapperWheelRowIsContentEmpty($tr)) {
      return;
    }
    mapperWheelDestroyRowReceptorWidgetForTr($tr);
    $tr.remove();
    mapperWheelEnsureTrailingBlankRow(true);
  });

  $('#mapper-wheel-input-tbody').on('input', '.mapper-wheel-in-value', function() {
    var $tr = $(this).closest('tr');
    mapperWheelSyncValueCellHint($tr);
    mapperWheelSyncRowSwatch($tr);
    if (mapperWheelIsTextMode()) {
      mapperWheelScheduleLabelColorPanelRebuildSoon();
    }
    mapperWheelEnsureTrailingBlankRow(true);
  });

  $('#mapper-wheel-input-tbody').on('paste', '.mapper-wheel-in-value', function(e) {
    if (mapperWheelIsTextMode()) {
      return;
    }
    var ev = e.originalEvent || e;
    var text = ev.clipboardData ? ev.clipboardData.getData('text/plain') : '';
    if (text != null && (text.indexOf('\t') !== -1 || /\r|\n/.test(text))) {
      return;
    }
    var take = text == null ? '' : String(text);
    take = mapperWheelTrimmedValueFromRowRaw(take.replace(/\r\n/g, '\n').split('\n')[0] || '');
    if (!take || !/[,.]/.test(take)) {
      return;
    }
    e.preventDefault();
    e.stopPropagation();
    var normalized = mapperWheelNormalizeLocalizedNumberString(take);
    var insertMe = mapperWheelNumericStringAcceptable(normalized) ? normalized : take;
    var el = this;
    var rawVal = $(el).val();
    rawVal = rawVal != null ? String(rawVal) : '';
    var start = typeof el.selectionStart === 'number' ? el.selectionStart : rawVal.length;
    var end = typeof el.selectionEnd === 'number' ? el.selectionEnd : start;
    var newVal = rawVal.slice(0, start) + insertMe + rawVal.slice(end);
    $(el).val(newVal).trigger('input');
    var pos = start + insertMe.length;
    if (typeof el.setSelectionRange === 'function') {
      window.requestAnimationFrame(function() {
        try {
          el.setSelectionRange(pos, pos);
        } catch (err) {}
      });
    }
  });

  $('#mapper-wheel-input-tbody').on('blur', '.mapper-wheel-in-value', function() {
    var $inp = $(this);
    var $tr = $inp.closest('tr');
    var raw = mapperWheelTrimmedValueFromRowRaw($inp.val());
    if (!mapperWheelIsTextMode() && raw && /[,.]/.test(raw)) {
      var norm = mapperWheelNormalizeLocalizedNumberString(raw);
      if (norm !== raw && mapperWheelNumericStringAcceptable(norm)) {
        $inp.val(norm);
      }
    }
    mapperWheelSyncValueCellHint($tr);
    mapperWheelSyncRowSwatch($tr);
    mapperWheelScheduleRebuildPlot();
  });

  $('#mapper-wheel-input-tbody').on('blur', '.mapper-core-in-receptor', function() {
    mapperWheelScheduleRebuildPlot();
  });

  $('#mapper-wheel-input-tbody').on('click', '.mapper-core-receptor-html-view', function(e) {
    e.preventDefault();
    mapperWheelBeginEditReceptor($(this).closest('tr'));
  });

  $('#mapper-wheel-input-tbody').on('keydown', '.mapper-core-receptor-html-view', function(e) {
    if (e.which === 13 || e.which === 32) {
      e.preventDefault();
      mapperWheelBeginEditReceptor($(this).closest('tr'));
    }
  });

  $('#mapper-wheel-input-tbody').on('keydown', '.mapper-core-in-receptor, .mapper-core-receptor-html-view, .mapper-wheel-in-value', function(e) {
    if (e.which !== 9 || e.shiftKey || e.altKey || e.ctrlKey || e.metaKey) {
      return;
    }
    var $target = $(this);
    if ($target.hasClass('mapper-core-in-receptor') && $target.data('ui-autocomplete')) {
      var $menu = $target.autocomplete('widget');
      if ($menu && $menu.is(':visible') && $menu.find('.ui-state-focus, .ui-state-active').length) {
        return;
      }
    }
    e.preventDefault();
    mapperWheelFocusNextTableCellFrom($target);
  });

  function mapperWheelClearWhole() {
    mapperWheelReceptorTallPlaceholderIntro = true;
    mapperWheelDestroyRowReceptorWidgets();
    $('#mapper-wheel-input-tbody').empty();
    $('#mapper-core-wheel-errors').empty();
    MAPPER_WHEEL_MODE_SNAPSHOTS.numeric = null;
    MAPPER_WHEEL_MODE_SNAPSHOTS.text = null;
    if (mapperWheelIsTextMode()) { MAPPER_CORE_LABEL_COLORS = {}; MAPPER_CORE_LABEL_ENABLED = {}; }
    mapperWheelAppendRow('', '');
    $('#mapper-wheel-clear-rows').addClass('mapper-core-clear-clean').blur();
  }
  function mapperWheelClearCurrentMode() {
    var curKey = mapperWheelIsTextMode() ? 'text' : 'numeric';
    mapperWheelReceptorTallPlaceholderIntro = true;
    mapperWheelDestroyRowReceptorWidgets();
    $('#mapper-wheel-input-tbody').empty();
    $('#mapper-core-wheel-errors').empty();
    MAPPER_WHEEL_MODE_SNAPSHOTS[curKey] = null;
    if (mapperWheelIsTextMode()) { MAPPER_CORE_LABEL_COLORS = {}; MAPPER_CORE_LABEL_ENABLED = {}; }
    mapperWheelAppendRow('', '');
    $('#mapper-wheel-clear-rows').addClass('mapper-core-clear-clean').blur();
  }
  $('#mapper-wheel-clear-whole').on('click', function(e) { e.preventDefault(); mapperWheelClearWhole(); });
  $('#mapper-wheel-clear-mode').on('click', function(e) { e.preventDefault(); mapperWheelClearCurrentMode(); });
  $('.mapper-core-clear-menu').on('click', '[data-col]', function(e) {
    e.preventDefault();
    $('#mapper-wheel-input-tbody tr').each(function() {
      $(this).find('.mapper-wheel-in-value').val('').trigger('input');
    });
    $('#mapper-wheel-clear-rows').removeClass('mapper-core-clear-clean');
  });
  mapperWheelSyncClearDropdown();

  $('#mapper-core-input-table').on('input change', '.mapper-wheel-receptor-input, .mapper-wheel-value-input, .mapper-core-row-color-picker', function() {
    $('#mapper-wheel-clear-rows').removeClass('mapper-core-clear-clean');
  });

  $('#mapper-core-input-table').on('click', 'th.mapper-core-sortable-head', function() {
    mapperWheelSortSerializedRows($(this).attr('data-mapper-wheel-sort-col'));
  });

  $('#mapper-core-input-table').on('keydown', 'th.mapper-core-sortable-head', function(e) {
    if (e.which === 13 || e.which === 32) {
      e.preventDefault();
      mapperWheelSortSerializedRows($(this).attr('data-mapper-wheel-sort-col'));
    }
  });

  $('#mapper-wheel-demo-rows').on('click', function(e) {
    e.preventDefault();
    $('#mapper-wheel-clear-rows').removeClass('mapper-core-clear-clean');
    mapperWheelFillDemoRows();
  });

  $('#mapper-core-input-table').on('paste', function(e) {
    $('#mapper-wheel-clear-rows').removeClass('mapper-core-clear-clean');
    var ev = e.originalEvent || e;
    var text = ev.clipboardData.getData('text/plain');
    if (!text) {
      return;
    }

    var lines = mapperWheelNormalizePastedLines(text);
    var contentLines = lines.filter(function(line) {
      return String(line).trim() !== '';
    });
    if (!contentLines.length) {
      return;
    }

    var hasTab = text.indexOf('\t') !== -1;
    var hasNewline = text.indexOf('\n') !== -1 || text.indexOf('\r') !== -1;
    var hasSemiGrid = mapperWheelPasteLooksSemicolonCsvGrid(contentLines, hasTab);
    if (!hasTab && !hasNewline && !hasSemiGrid) {
      return;
    }

    var dualColumn = hasTab || hasSemiGrid;
    if (!dualColumn && !hasNewline) {
      return;
    }

    var anchor = mapperWheelPasteAnchorFromEvent(e);
    var startIdx = anchor.$tr && anchor.$tr.length ? anchor.$tr.index() : -1;
    var maxRows = MAPPER_WHEEL_MAX_ROWS;
    e.preventDefault();

    function finishPasteBatch() {
      mapperWheelReceptorTallPlaceholderIntro = false;
      mapperWheelEnsureTrailingBlankRow();
    }

    if (!dualColumn) {
      var colSingle = anchor.col || 'value';
      if (startIdx < 0) {
        var curCountSingle = $('#mapper-wheel-input-tbody tr').length;
        if (curCountSingle + contentLines.length > maxRows) {
          alert('Cannot paste: would exceed ' + maxRows + ' rows.');
          return;
        }
        contentLines.forEach(function(cellLine) {
          var cell = String(cellLine).trim();
          if (colSingle === 'receptor') {
            mapperWheelAppendRowFromPaste(cell, '');
          } else {
            mapperWheelAppendRow('', cell, true);
          }
        });
        finishPasteBatch();
        return;
      }
      if (startIdx + contentLines.length > maxRows) {
        alert('Cannot paste: would exceed ' + maxRows + ' rows.');
        return;
      }
      while ($('#mapper-wheel-input-tbody tr').length < startIdx + contentLines.length) {
        mapperWheelAppendRow('', '', true);
      }
      for (var si = 0; si < contentLines.length; si++) {
        var cellOne = String(contentLines[si]).trim();
        var $trSi = $('#mapper-wheel-input-tbody tr').eq(startIdx + si);
        if (colSingle === 'receptor') {
          mapperWheelSetReceptorRawOnRow($trSi, cellOne);
        } else {
          mapperWheelSetValueOnRow($trSi, cellOne);
        }
      }
      finishPasteBatch();
      return;
    }

    if (startIdx < 0) {
      var curCountDc = $('#mapper-wheel-input-tbody tr').length;
      if (curCountDc + contentLines.length > maxRows) {
        alert('Cannot paste: would exceed ' + maxRows + ' rows.');
        return;
      }
      contentLines.forEach(function(line0) {
        var sp = mapperWheelSplitPasteLineIntoReceptorAndValue(line0);
        if (!sp.rawR && !sp.valCell) {
          return;
        }
        mapperWheelAppendRowFromPaste(sp.rawR, sp.valCell);
      });
      finishPasteBatch();
      return;
    }
    if (startIdx + contentLines.length > maxRows) {
      alert('Cannot paste: would exceed ' + maxRows + ' rows.');
      return;
    }
    while ($('#mapper-wheel-input-tbody tr').length < startIdx + contentLines.length) {
      mapperWheelAppendRow('', '', true);
    }
    for (var ji = 0; ji < contentLines.length; ji++) {
      var sp2 = mapperWheelSplitPasteLineIntoReceptorAndValue(contentLines[ji]);
      var $trJj = $('#mapper-wheel-input-tbody tr').eq(startIdx + ji);
      mapperWheelSetReceptorRawOnRow($trJj, sp2.rawR);
      mapperWheelSetValueOnRow($trJj, sp2.valCell);
    }
    finishPasteBatch();
  });

  if (typeof window.mapperCoreInitGpcromePickerModal === 'function') {
    window.mapperCoreInitGpcromePickerModal({
      pickerRows: $.isArray(MAPPER_CORE_GPCROME_PICKER_ROWS) ? MAPPER_CORE_GPCROME_PICKER_ROWS : [],
      maxRows: MAPPER_WHEEL_MAX_ROWS,
      onAdd: function(entryIds, meta) {
        if (!entryIds || !entryIds.length) {
          return;
        }

        function findFirstEmptyRow() {
          var $hit = $([]);
          $('#mapper-wheel-input-tbody tr').each(function() {
            if (mapperWheelRowIsContentEmpty($(this))) {
              $hit = $(this);
              return false;
            }
          });
          return $hit;
        }

        function tbodyAtOrOverCap() {
          return $('#mapper-wheel-input-tbody tr').length >= MAPPER_WHEEL_MAX_ROWS;
        }

        var numberMap = (meta && meta.groupNameById && !mapperWheelIsTextMode())
          ? window.mapperCoreBuildSequentialNumberMap(meta.groupNameById)
          : null;

        var added = 0;
        for (var ix = 0; ix < entryIds.length; ix++) {
          var id = entryIds[ix];
          if (tbodyAtOrOverCap() && !findFirstEmptyRow().length) {
            alert('Stopped at ' + MAPPER_WHEEL_MAX_ROWS + ' rows. Added ' + added + ' receptor(s); remaining selection was skipped.');
            break;
          }
          var $tr = findFirstEmptyRow();
          if (!$tr.length) {
            mapperWheelAppendRow('', '', true);
            $tr = $('#mapper-wheel-input-tbody tr').last();
          }
          mapperWheelSetReceptorResolvedUI($tr, id);
          if (meta && meta.groupNameById && meta.groupNameById[id] != null) {
            var assignedValue = mapperWheelIsTextMode() ? meta.groupNameById[id] : numberMap[id];
            if (assignedValue != null) {
              $tr.find('.mapper-wheel-in-value').val(String(assignedValue));
            }
          }
          mapperWheelSyncValueCellHint($tr);
          mapperWheelSyncRowSwatch($tr);
          added++;
        }
        mapperWheelEnsureTrailingBlankRow();
      }
    });
  }
  $('.mapper-core-booting').removeClass('mapper-core-booting');
});
function mapperWheelDownloadTable() {
  var isText = $('#mapper-core-input-table').hasClass('mapper-core-text-mode');
  var headers = ['Receptor', isText ? 'Category' : 'Number'];
  var rows = [headers];
  $('#mapper-wheel-input-tbody tr').each(function() {
    var $tr = $(this);
    var entryId = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
    var unmatched = $tr.data('mapperCoreUnmatchedRaw');
    var ta = ($tr.find('.mapper-core-in-receptor').val() || '').trim();
    var receptor = '';
    if (entryId && MAPPER_CORE_ENTRY_META && MAPPER_CORE_ENTRY_META[entryId]) {
      receptor = MAPPER_CORE_ENTRY_META[entryId].name_plain || MAPPER_CORE_ENTRY_META[entryId].text || entryId;
    } else if (unmatched != null && String(unmatched).trim()) {
      receptor = String(unmatched).trim();
    } else {
      receptor = ta;
    }
    if (!receptor) return;
    var val = ($tr.find('.mapper-wheel-in-value').val() || '').trim();
    rows.push([receptor, val]);
  });
  var wb = XLSX.utils.book_new();
  var ws = XLSX.utils.aoa_to_sheet(rows);
  XLSX.utils.book_append_sheet(wb, ws, 'GPCRome Wheel');
  XLSX.writeFile(wb, 'GPCRome_Wheel_data.xlsx');
}

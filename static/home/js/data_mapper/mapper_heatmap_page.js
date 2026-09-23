/**
 * Mapper 2.0 — Heatmap page.
 * Left panel: receptor + V1-V5 numerical input table.
 * Right panel: Heatmap() from datamapper.js.
 * Data keys are entry_names (e.g. HTR1A_HUMAN); Heatmap() uses label_converter to show names.
 */
(function ($) {
  'use strict';

  var MAPPER_HEATMAP_MAX_ROWS  = 500;
  var MAPPER_HEATMAP_PLACEHOLDER = 'Paste in 1–6 columns of data\nor type';
  var DEBOUNCE_MS = 200;

  var suppressRedraw = false;

  // ── Global data structures ────────────────────────────────────────────────
  /** {entry_name → name_plain} and {entry_name → gene} — accessed as global by Heatmap() */
  var LabelConverter = {
    UniProt_to_IUPHAR_converter: {},
    IUPHAR_to_UniProt_converter: {},
    UniProt_to_Gene_converter:   {}
  };

  /** V1-V5 column display labels */
  var LabelXConverter = { Value1: 'Dataset 1', Value2: 'Dataset 2', Value3: 'Dataset 3', Value4: 'Dataset 4', Value5: 'Dataset 5' };

  /** HeatmapDataStyling instance — initialised at boot from heatmap_DataStyling() */
  var HeatmapDataStyling;

  var SORT_STATE = { col: null, dir: 'asc' };

  var mapperHeatmapSpeciesActive = false;
  var mapperHeatmapSpeciesNameFormat = 'common'; // 'common' | 'latin' | 'both'

  // Demo rows
  var MAPPER_HEATMAP_DEMO_ROWS = [
    { receptor: '5HT1A',  vals: [1,   73,  -10.1,  0.8,  15] },
    { receptor: '5HT1B',  vals: [2,   72,  -10.0,  0.7,  14] },
    { receptor: '5HT2A',  vals: [6,   68,   -9.6,  0.6,  12] },
    { receptor: 'ACKR1',  vals: [13,  61,   -8.9,  0.5,  10] },
    { receptor: 'ACKR2',  vals: [14,  60,   -8.8,  0.4,   9] },
    { receptor: 'ACM1',   vals: [17,  57,   -8.5,  0.3,   8] },
    { receptor: 'ACM2',   vals: [18,  56,   -8.4,  0.2,   7] },
    { receptor: 'ADA1A',  vals: [23,  51,   -7.9,  0.1,   5] },
    { receptor: 'ADRB1',  vals: [29,  45,   -7.3,  0.9,   3] },
    { receptor: 'ADRB2',  vals: [30,  44,   -7.2,  1.0,   2] }
  ];

  // ── Helpers ──────────────────────────────────────────────────────────────
  function parseNumLoose(s) {
    if (s == null || String(s).trim() === '') return null;
    var x = parseFloat(String(s).trim().replace(',', '.'));
    return isFinite(x) ? x : null;
  }

  // ── Label converter ───────────────────────────────────────────────────────
  function mapperHeatmapBuildLabelConverter() {
    LabelConverter = { UniProt_to_IUPHAR_converter: {}, IUPHAR_to_UniProt_converter: {}, UniProt_to_Gene_converter: {} };
    var meta = window.MAPPER_CORE_ENTRY_META || {};
    Object.keys(meta).forEach(function (entry_name) {
      var m = meta[entry_name];
      if (m.name_plain) {
        LabelConverter.UniProt_to_IUPHAR_converter[entry_name] = m.name_plain;
        LabelConverter.IUPHAR_to_UniProt_converter[m.name_plain] = entry_name;
      }
      if (m.gene) LabelConverter.UniProt_to_Gene_converter[entry_name] = m.gene;
    });
  }

  // ── Species helpers ───────────────────────────────────────────────────────
  function mapperHeatmapSpeciesLabel(o) {
    if (mapperHeatmapSpeciesNameFormat === 'latin') return o.latin || o.common || '';
    if (mapperHeatmapSpeciesNameFormat === 'both') {
      var c = (o.common || '').trim(), l = (o.latin || '').trim();
      return (c && l && c !== l) ? c + ' (' + l + ')' : c || l;
    }
    return o.common || o.latin || '';
  }

  function mapperHeatmapFormatSpeciesTag(o) {
    var c = (o.common || '').trim(), l = (o.latin || '').trim();
    if (mapperHeatmapSpeciesNameFormat === 'latin') return l ? '(' + l + ')' : '';
    if (mapperHeatmapSpeciesNameFormat === 'both') {
      if (c && l && c !== l) return '(' + c + ', ' + l + ')';
      return (c || l) ? '(' + (c || l) + ')' : '';
    }
    return (c || l) ? '(' + (c || l) + ')' : '';
  }

  function mapperHeatmapSpeciesLabelForEntry(specEntry) {
    var sd = window.MAPPER_HEATMAP_SPECIES_DATA;
    if (!sd) return '';
    var stem = String(specEntry || '').split('_')[0];
    var pools = (sd.by_stem[stem] || []).concat(sd.nonhuman_only || []);
    for (var i = 0; i < pools.length; i++) {
      if (pools[i].entry === specEntry) return mapperHeatmapSpeciesLabel(pools[i]);
    }
    return '';
  }

  function mapperHeatmapPopulateSpeciesSelect($tr, entryId, preserveSelection) {
    var $sel = $tr.find('.mapper-heatmap-species-select');
    if (!$sel.length || !entryId || !window.MAPPER_HEATMAP_SPECIES_DATA) return;
    var prevVal = preserveSelection ? ($sel.val() || '').trim() : '';
    if ($sel.data('select2')) { try { $sel.select2('destroy'); } catch (e) {} }
    $sel.empty();
    var sd   = window.MAPPER_HEATMAP_SPECIES_DATA;
    var stem = String(entryId).split('_')[0];
    var orthologs = (sd.by_stem[stem] || []).slice();
    if (!orthologs.length) {
      (sd.nonhuman_only || []).forEach(function (nho) { if (nho.stem === stem) orthologs.push(nho); });
    }
    orthologs.forEach(function (o) {
      $('<option>').val(o.entry).text(mapperHeatmapSpeciesLabel(o)).appendTo($sel);
    });
    var validPrev = prevVal && orthologs.some(function (o) { return o.entry === prevVal; });
    if (validPrev) {
      $sel.val(prevVal);
    } else {
      var humanOpt = null;
      orthologs.forEach(function (o) { if (o.is_human) humanOpt = o; });
      $sel.val(humanOpt ? humanOpt.entry : (orthologs[0] ? orthologs[0].entry : ''));
    }
    $sel.select2({ width: 'resolve', dropdownAutoWidth: true, minimumResultsForSearch: 6, dropdownParent: $('body') });
  }

  function mapperHeatmapInitAllSpeciesDropdowns() {
    $('#mapper-heatmap-input-tbody tr').each(function () {
      var entry = ($(this).find('.mapper-core-receptor-entry').val() || '').trim();
      if (entry) mapperHeatmapPopulateSpeciesSelect($(this), entry, true);
    });
  }

  function mapperHeatmapApplySpeciesPanelClass() {
    $('.mapper-core-wheel-wrap.mapper-heatmap-page').toggleClass('mapper-heatmap-species-active', mapperHeatmapSpeciesActive);
  }

  function mapperHeatmapExtendLabelConverterForSpecies() {
    var sd = window.MAPPER_HEATMAP_SPECIES_DATA;
    if (!sd || !sd.by_stem) return;
    function addEntry(o) {
      var stem = o.stem;
      var humanOrthologs = (sd.by_stem[stem] || []).filter(function (x) { return x.is_human; });
      var humanEntry = humanOrthologs.length ? humanOrthologs[0].entry : null;
      var baseMeta = humanEntry && window.MAPPER_CORE_ENTRY_META && window.MAPPER_CORE_ENTRY_META[humanEntry];
      var baseIUPHAR = (baseMeta && baseMeta.name_plain) || stem.toUpperCase();
      var gene = (baseMeta && baseMeta.gene) || '';
      var tag  = mapperHeatmapSpeciesActive ? mapperHeatmapFormatSpeciesTag(o) : (o.is_human ? '' : mapperHeatmapFormatSpeciesTag(o));
      var fullName = tag ? baseIUPHAR + ' ' + tag : baseIUPHAR;
      LabelConverter.UniProt_to_IUPHAR_converter[o.entry] = fullName;
      LabelConverter.IUPHAR_to_UniProt_converter[fullName] = o.entry;
      if (gene) LabelConverter.UniProt_to_Gene_converter[o.entry] = tag ? gene + ' ' + tag : gene;
    }
    Object.keys(sd.by_stem).forEach(function (stem) { (sd.by_stem[stem] || []).forEach(addEntry); });
    (sd.nonhuman_only || []).forEach(addEntry);
  }

  // ── Data construction ─────────────────────────────────────────────────────
  function mapperHeatmapBuildData() {
    var data = {};
    $('#mapper-heatmap-input-tbody tr').each(function () {
      var $tr   = $(this);
      var entry = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
      if (!entry) {
        var ta       = ($tr.find('.mapper-core-in-receptor').val() || '').trim();
        var unmatched = ($tr.data('mapperCoreUnmatchedRaw') || '');
        var up = String(ta || unmatched).trim().toUpperCase();
        if (up && window.MAPPER_CORE_RESOLVE) entry = window.MAPPER_CORE_RESOLVE[up] || '';
      }
      if (!entry) return;

      var v1 = parseNumLoose($tr.find('.mapper-heatmap-val1').val());
      var v2 = parseNumLoose($tr.find('.mapper-heatmap-val2').val());
      var v3 = parseNumLoose($tr.find('.mapper-heatmap-val3').val());
      var v4 = parseNumLoose($tr.find('.mapper-heatmap-val4').val());
      var v5 = parseNumLoose($tr.find('.mapper-heatmap-val5').val());
      if (v1 == null && v2 == null && v3 == null && v4 == null && v5 == null) return;

      var dp = {};
      if (v1 != null) dp.Value1 = v1;
      if (v2 != null) dp.Value2 = v2;
      if (v3 != null) dp.Value3 = v3;
      if (v4 != null) dp.Value4 = v4;
      if (v5 != null) dp.Value5 = v5;
      var dataKey = entry;
      if (mapperHeatmapSpeciesActive) {
        var specVal = ($tr.find('.mapper-heatmap-species-select').val() || '').trim();
        if (specVal) dataKey = specVal;
      }
      data[dataKey] = dp;  // key = entry_name (species-aware when active)
    });
    return data;
  }

  // ── Placeholder ───────────────────────────────────────────────────────────
  function mapperHeatmapShowPlaceholder() {
    $('#mapper-heatmap-plot').empty().append(
      $('<div class="mapper-heatmap-plot-placeholder">').append(
        $('<p class="mapper-heatmap-plot-placeholder-title">').text('No receptors mapped yet'),
        $('<p class="mapper-heatmap-plot-placeholder-hint">').text(
          'Add receptors in the left panel and enter values (V1–V5) — the heatmap renders automatically.'
        )
      )
    );
  }

  function mapperHeatmapSyncFirstRowPlaceholder() {
    $('#mapper-heatmap-input-tbody tr').each(function (idx) {
      var $inp = $(this).find('.mapper-core-in-receptor');
      if (!$inp.length) return;
      if (idx === 0) $inp.attr('placeholder', MAPPER_HEATMAP_PLACEHOLDER);
      else           $inp.removeAttr('placeholder');
    });
  }

  // ── Render ────────────────────────────────────────────────────────────────
  function mapperHeatmapRedrawNow() {
    var data = mapperHeatmapBuildData();
    if (!Object.keys(data).length) { mapperHeatmapShowPlaceholder(); return; }

    $('#mapper-heatmap-plot').empty();

    // Expose label_converter as global — Heatmap() reads it directly
    window.label_converter = LabelConverter;

    try {
      Heatmap(data, 'mapper-heatmap-plot', HeatmapDataStyling, LabelXConverter);
    } catch (e) {
      mapperHeatmapShowPlaceholder();
    }
  }

  var _heatmapRedrawDebounced = MapperPageCore.debounce(function () { mapperHeatmapRedrawNow(); }, DEBOUNCE_MS);
  function mapperHeatmapScheduleRedraw() {
    if (suppressRedraw) return;
    _heatmapRedrawDebounced.schedule();
  }

  // ── Colour pickers ────────────────────────────────────────────────────────
  var HM_PRESETS = {
    One:       { setup: 'One',   start: '#ffffff', mid: '#ffffff', end: '#707070' },
    Two:       { setup: 'Two',   start: '#97a6c4', mid: '#ffffff', end: '#384860' },
    Three_RWB: { setup: 'Three', start: '#1a80bb', mid: '#ffffff', end: '#a00000' },
    Three_TWM: { setup: 'Three', start: '#298c8c', mid: '#ffffff', end: '#800074' }
  };

  function applyHmColorPreset(key, skipRedraw) {
    var p = HM_PRESETS[key] || HM_PRESETS.Three_RWB;
    HeatmapDataStyling.Number_of_colors = p.setup;
    HeatmapDataStyling.min_color    = p.start;
    HeatmapDataStyling.middle_color = p.mid;
    HeatmapDataStyling.max_color    = p.end;

    var isOne   = p.setup === 'One';
    var isThree = p.setup === 'Three';
    $('#hm-picker-min-wrap').css('visibility', isOne   ? 'hidden' : 'visible');
    $('#hm-label-min').css('visibility',       isOne   ? 'hidden' : 'visible');
    $('#hm-picker-mid-wrap').css('visibility', isThree ? 'visible' : 'hidden');
    $('#hm-label-mid').css('visibility',       isThree ? 'visible' : 'hidden');

    if ($.fn.spectrum) {
      try { $('#mapper-heatmap-cpicker-min').spectrum('set', p.start); } catch (e) {}
      try { $('#mapper-heatmap-cpicker-mid').spectrum('set', p.mid);   } catch (e) {}
      try { $('#mapper-heatmap-cpicker-max').spectrum('set', p.end);   } catch (e) {}
    }
    if (!skipRedraw) mapperHeatmapScheduleRedraw();
  }

  function formatHmColorOption(option) {
    if (!option.id) return option.text;
    var cmap = { One: [null,null,'#707070'], Two: ['#97a6c4',null,'#384860'],
                 Three_RWB: ['#1a80bb','#ffffff','#a00000'], Three_TWM: ['#298c8c','#ffffff','#800074'] };
    var colors = cmap[option.id];
    if (!colors) return option.text;
    return $('<span style="display:flex;align-items:center;">').append(
      $('<span style="min-width:50px;text-align:center;margin-right:5px;">').text(option.text),
      colors.map(function (c) {
        return $('<span style="display:inline-block;width:12px;height:12px;margin-left:4px;border:1px solid #ccc;border-radius:2px;background:' +
          (c || 'transparent') + ';' + (c ? '' : 'opacity:0') + ';">');
      })
    );
  }

  function mapperHeatmapInitColorPickers() {
    if (!$.fn.spectrum || !$.fn.select2) return;

    // Style select2
    $('#mapper-heatmap-color-style').select2({
      templateResult: formatHmColorOption, templateSelection: formatHmColorOption, width: 'resolve'
    }).on('change', function () { applyHmColorPreset($(this).val()); });

    var spOpts = {
      showInput: true, showPalette: false, preferredFormat: 'hex',
      palette: [['#1a80bb','#a00000','#298c8c','#800074','#ffffff','#000000']],
      change: function () {
        HeatmapDataStyling.min_color    = $('#mapper-heatmap-cpicker-min').spectrum('get').toHexString();
        HeatmapDataStyling.middle_color = $('#mapper-heatmap-cpicker-mid').spectrum('get').toHexString();
        HeatmapDataStyling.max_color    = $('#mapper-heatmap-cpicker-max').spectrum('get').toHexString();
        mapperHeatmapScheduleRedraw();
      }
    };
    $.extend(true, {}, spOpts);
    $('#mapper-heatmap-cpicker-min').spectrum($.extend({}, spOpts, { color: HM_PRESETS.Three_RWB.start }));
    $('#mapper-heatmap-cpicker-mid').spectrum($.extend({}, spOpts, { color: HM_PRESETS.Three_RWB.mid  }));
    $('#mapper-heatmap-cpicker-max').spectrum($.extend({}, spOpts, { color: HM_PRESETS.Three_RWB.end  }));
    // Add move handler
    $('#mapper-heatmap-cpicker-min, #mapper-heatmap-cpicker-mid, #mapper-heatmap-cpicker-max').spectrum('option', 'move', function () {
      HeatmapDataStyling.min_color    = $('#mapper-heatmap-cpicker-min').spectrum('get').toHexString();
      HeatmapDataStyling.middle_color = $('#mapper-heatmap-cpicker-mid').spectrum('get').toHexString();
      HeatmapDataStyling.max_color    = $('#mapper-heatmap-cpicker-max').spectrum('get').toHexString();
      mapperHeatmapScheduleRedraw();
    });

    // Apply default preset
    $('#mapper-heatmap-color-style').val('Three_RWB').trigger('change');
  }

  // ── Row management ────────────────────────────────────────────────────────
  function mapperHeatmapCreateReceptorTd($td) {
    var $hid  = $('<input type="hidden" class="mapper-core-receptor-entry" value="">');
    var $wrap = $('<div class="mapper-core-receptor-input-wrap is-empty">');
    var $inp  = $('<textarea class="form-control input-sm mapper-core-in-receptor" rows="1" autocomplete="off" spellcheck="false"></textarea>');
    var $view = $('<div class="form-control input-sm mapper-core-receptor-html-view" tabindex="0"></div>');
    var $clr  = $('<button type="button" class="mapper-core-receptor-clear" aria-label="Clear receptor">&times;</button>');
    $wrap.append($inp, $view, $clr);
    $td.append($hid, $wrap);
    $view.hide();
    $clr.on('click', function () { mapperHeatmapSetResolved($clr.closest('tr'), ''); });
    return $inp;
  }

  function mapperHeatmapRowBlank($tr) {
    var entry   = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
    var typed   = ($tr.find('.mapper-core-in-receptor').val() || '').trim();
    var hasVals = [1,2,3,4,5].some(function (i) { return (($tr.find('.mapper-heatmap-val' + i).val() || '').trim() !== ''); });
    return !entry && !typed && !hasVals && !$tr.data('mapperCoreUnmatchedRaw');
  }

  function mapperHeatmapSyncClearBtn($tr) {
    var has = !!($tr.find('.mapper-core-receptor-entry').val() || '').trim()
           || !!($tr.find('.mapper-core-in-receptor').val() || '').trim()
           || !!($tr.data('mapperCoreUnmatchedRaw'));
    $tr.find('.mapper-core-receptor-input-wrap').toggleClass('is-empty', !has);
  }

  function mapperHeatmapSyncRemoveButtons() {
    $('#mapper-heatmap-input-tbody tr').each(function () {
      $(this).find('.mapper-core-remove-cell').toggleClass('is-remove-hidden', mapperHeatmapRowBlank($(this)));
    });
  }

  function mapperHeatmapCompactReceptors() {
    var hasContent = false;
    $('#mapper-heatmap-input-tbody tr').each(function () { if (!mapperHeatmapRowBlank($(this))) { hasContent = true; return false; } });
    $('#mapper-heatmap-input-table').toggleClass('mapper-core-receptors-compact', hasContent);
  }

  function mapperHeatmapAppendRow(skipTrail) {
    var $tr = $('<tr>');
    $tr.append($('<td class="mapper-core-remove-cell is-remove-hidden">').append(
      $('<button type="button" class="mapper-core-remove-row" aria-label="Remove row">&times;</button>')
    ));
    var $tdR = $('<td class="mapper-core-receptor-cell mapper-core-value-cell">');
    var $inp = mapperHeatmapCreateReceptorTd($tdR);
    $tr.append($tdR);
    // Species cell (hidden until species toggle is active)
    var $specSel = $('<select class="form-control input-sm mapper-heatmap-species-select">');
    $tr.append($('<td class="mapper-heatmap-species-cell">').css('display', mapperHeatmapSpeciesActive ? '' : 'none').append($specSel));

    [1,2,3,4,5].forEach(function (i) {
      $tr.append($('<td class="mapper-heatmap-val-cell mapper-core-value-cell">').append(
        $('<input type="text" class="form-control input-sm mapper-heatmap-val mapper-heatmap-val' + i + '" autocomplete="off">')
      ));
    });
    $('#mapper-heatmap-input-tbody').append($tr);
    mapperHeatmapBindAc($inp);
    mapperHeatmapSyncClearBtn($tr);
    mapperHeatmapSyncRemoveButtons();
    if (!skipTrail) mapperHeatmapEnsureTrailingBlankRow();
  }

  function mapperHeatmapEnsureTrailingBlankRow() {
    var $tb = $('#mapper-heatmap-input-tbody');
    while ($tb.children().length >= 2) {
      var $last = $tb.children().last(), $prev = $last.prev();
      if (mapperHeatmapRowBlank($last) && mapperHeatmapRowBlank($prev)) { mapperHeatmapDestroyRowAc($last); $last.remove(); continue; }
      break;
    }
    if (!$tb.children().length) mapperHeatmapAppendRow(true);
    if (!mapperHeatmapRowBlank($tb.children().last()) && $tb.children().length < MAPPER_HEATMAP_MAX_ROWS) mapperHeatmapAppendRow(true);
  }

  function mapperHeatmapDestroyAc($inp) { try { if ($inp.hasClass('ui-autocomplete-input')) $inp.autocomplete('destroy'); } catch (e) {} }
  function mapperHeatmapDestroyRowAc($tr) {
    mapperHeatmapDestroyAc($tr.find('.mapper-core-in-receptor'));
    var $ss = $tr.find('.mapper-heatmap-species-select');
    if ($ss.length && $ss.data('select2')) { try { $ss.select2('destroy'); } catch (e) {} }
  }
  function mapperHeatmapDestroyAllRows() { $('#mapper-heatmap-input-tbody tr').each(function () { mapperHeatmapDestroyRowAc($(this)); }); }

  // ── Autocomplete ──────────────────────────────────────────────────────────
  function mapperHeatmapResolvedDisplay(entryId) {
    var meta = window.MAPPER_CORE_ENTRY_META && window.MAPPER_CORE_ENTRY_META[entryId] || {};
    var lt   = HeatmapDataStyling ? HeatmapDataStyling.LabelType : 'IUPHAR';
    if (lt === 'Gene' && meta.gene) return $('<span/>').text(meta.gene).html();
    if (lt === 'UniProt' && meta.uniprot) return $('<span/>').text(meta.uniprot).html();
    return meta.name_html ? String(meta.name_html) : $('<span/>').text(entryId || '').html();
  }

  function mapperHeatmapEditSeed(entryId) {
    var sid = entryId != null ? String(entryId).trim() : '';
    var meta = (window.MAPPER_CORE_ENTRY_META && window.MAPPER_CORE_ENTRY_META[sid]) || {};
    var lt = HeatmapDataStyling ? HeatmapDataStyling.LabelType : 'IUPHAR';
    if (lt === 'Gene' && meta.gene) { return meta.gene; }
    if (lt === 'UniProt' && meta.uniprot) { return meta.uniprot; }
    return meta.name_plain || '';
  }

  function mapperHeatmapSetResolved($tr, id) {
    var sid  = id != null ? String(id).trim() : '';
    var $inp  = $tr.find('.mapper-core-in-receptor');
    var $hid  = $tr.find('.mapper-core-receptor-entry');
    var $view = $tr.find('.mapper-core-receptor-html-view');
    mapperHeatmapDestroyAc($inp);
    if (!sid) {
      $hid.val(''); $view.hide().empty(); $inp.val('').show();
      mapperHeatmapBindAc($inp); $tr.removeClass('mapper-core-row-invalid').removeData('mapperCoreUnmatchedRaw');
      mapperHeatmapSyncClearBtn($tr); mapperHeatmapSyncRemoveButtons(); mapperHeatmapCompactReceptors();
      mapperHeatmapScheduleRedraw(); return;
    }
    $hid.val(sid);
    if (mapperHeatmapSpeciesActive) mapperHeatmapPopulateSpeciesSelect($tr, sid, false);
    $view.html(mapperHeatmapResolvedDisplay(sid)).show(); $inp.val('').hide();
    mapperHeatmapBindAc($inp);
    $tr.removeClass('mapper-core-row-invalid').removeData('mapperCoreUnmatchedRaw');
    mapperHeatmapSyncClearBtn($tr); mapperHeatmapSyncRemoveButtons(); mapperHeatmapCompactReceptors();
    if (!suppressRedraw) {
      window.setTimeout(function () { var $v = $tr.find('.mapper-heatmap-val1:visible').first(); if ($v.length) $v.focus().select(); }, 0);
    }
    mapperHeatmapScheduleRedraw();
  }

  function mapperHeatmapFilterLocal(term) {
    var t = (term || '').trim().toUpperCase();
    if (!t || !window.receptorSelect2Data) return [];
    var meta = window.MAPPER_CORE_ENTRY_META || {};
    var lt = HeatmapDataStyling ? HeatmapDataStyling.LabelType : 'IUPHAR';
    var results = [];
    window.receptorSelect2Data.forEach(function (item) {
      var m = meta[item.id] || {};
      var searchIn = [(m.name_plain || item.name_plain || item.text || ''), (m.gene || ''), (m.uniprot || ''), item.id].join(' ').toUpperCase();
      if (searchIn.indexOf(t) === -1) return;
      var dispHtml;
      if (lt === 'Gene' && m.gene) {
        dispHtml = $('<span/>').text(m.gene).html();
      } else if (lt === 'UniProt' && m.uniprot) {
        dispHtml = $('<span/>').text(m.uniprot).html();
      } else {
        dispHtml = m.name_html || item.name_html || $('<span/>').text(m.name_plain || item.text || item.id).html();
      }
      results.push({ label: m.name_plain || item.name_plain || item.text || item.id,
        value: item.id, id: item.id, html: dispHtml, name_html: item.name_html || '', name_plain: item.name_plain || '' });
    });
    if (results.length > 80) results = results.slice(0, 80);
    return results;
  }

  function mapperHeatmapBindAc($inp) {
    mapperHeatmapDestroyAc($inp);
    $inp.autocomplete({
      minLength: 1,
      source: function (req, resp) { resp(mapperHeatmapFilterLocal(req.term)); },
      focus: function () { return false; },
      select: function (ev, ui) { mapperHeatmapSetResolved($inp.closest('tr'), ui.item.id); ev.preventDefault(); }
    });
    var w = $inp.data('ui-autocomplete');
    if (w) {
      w._renderItem = function (ul, item) {
        var inner = item.html || item.name_html || $('<span/>').text(item.label || '').html();
        return $('<li>').append($('<div class="mapper-core-ac-item-label">').html(inner)).appendTo(ul);
      };
      if (w.menu && w.menu.element) w.menu.element.addClass('mapper-core-receptor-ac-menu');
    }
    $inp.on('keyup', function () {
      var $tr = $inp.closest('tr');
      if (!$tr.find('.mapper-core-receptor-entry').val()) { $tr.removeData('mapperCoreUnmatchedRaw'); $tr.removeClass('mapper-core-row-invalid'); }
      mapperHeatmapSyncClearBtn($tr);
    });
    $inp.on('blur.mapperHeatmap', function () {
      var $tr = $inp.closest('tr');
      window.setTimeout(function () {
        if (!$inp.is(':visible')) return;
        var raw = ($inp.val() || '').trim();
        if (!raw || $tr.find('.mapper-core-receptor-entry').val()) return;
        var rid = window.MAPPER_CORE_RESOLVE && window.MAPPER_CORE_RESOLVE[raw.toUpperCase()];
        if (rid) { mapperHeatmapSetResolved($tr, rid); return; }
        var hit = (window.receptorSelect2Data||[]).filter(function(x){ return String(x.id)===raw; });
        if (hit.length) { mapperHeatmapSetResolved($tr, hit[0].id); return; }
        if (raw) { $tr.addClass('mapper-core-row-invalid'); $tr.data('mapperCoreUnmatchedRaw', raw); }
        mapperHeatmapSyncClearBtn($tr);
      }, 170);
    });
  }

  // ── Sort ──────────────────────────────────────────────────────────────────
  function mapperHeatmapSerializeRows() {
    var rows = [];
    $('#mapper-heatmap-input-tbody tr').each(function (ix) {
      var $tr = $(this);
      rows.push({
        ix: ix,
        entry: ($tr.find('.mapper-core-receptor-entry').val()||'').trim(),
        typed: ($tr.find('.mapper-core-in-receptor').val()||'').trim(),
        unmatched: ($tr.data('mapperCoreUnmatchedRaw')||''),
        invalid: $tr.hasClass('mapper-core-row-invalid'),
        v1: ($tr.find('.mapper-heatmap-val1').val()||'').trim(),
        v2: ($tr.find('.mapper-heatmap-val2').val()||'').trim(),
        v3: ($tr.find('.mapper-heatmap-val3').val()||'').trim(),
        v4: ($tr.find('.mapper-heatmap-val4').val()||'').trim(),
        v5: ($tr.find('.mapper-heatmap-val5').val()||'').trim()
      });
    });
    return rows;
  }

  function mapperHeatmapApplySerializedRows(rows) {
    suppressRedraw = true;
    mapperHeatmapDestroyAllRows();
    $('#mapper-heatmap-input-tbody').empty();
    rows.forEach(function (r) {
      mapperHeatmapAppendRow(true);
      var $tr = $('#mapper-heatmap-input-tbody tr').last();
      if (r.entry) { mapperHeatmapSetResolved($tr, r.entry); }
      else if (r.typed) {
        $tr.find('.mapper-core-in-receptor').val(r.typed);
        if (r.unmatched) $tr.data('mapperCoreUnmatchedRaw', r.unmatched);
        if (r.invalid)   $tr.addClass('mapper-core-row-invalid');
      }
      $tr.find('.mapper-heatmap-val1').val(r.v1);
      $tr.find('.mapper-heatmap-val2').val(r.v2);
      $tr.find('.mapper-heatmap-val3').val(r.v3);
      $tr.find('.mapper-heatmap-val4').val(r.v4);
      $tr.find('.mapper-heatmap-val5').val(r.v5);
    });
    suppressRedraw = false;
    mapperHeatmapEnsureTrailingBlankRow();
    mapperHeatmapSyncRemoveButtons();
    mapperHeatmapCompactReceptors();
    mapperHeatmapSyncFirstRowPlaceholder();
  }

  function mapperHeatmapSortRows(col) {
    if (SORT_STATE.col === col) SORT_STATE.dir = SORT_STATE.dir === 'asc' ? 'desc' : 'asc';
    else { SORT_STATE.col = col; SORT_STATE.dir = 'asc'; }
    var rows = mapperHeatmapSerializeRows();
    rows.sort(function (a, b) {
      var va, vb, key = col.replace('val','v');
      if (col === 'receptor') { va = a.entry || a.typed || ''; vb = b.entry || b.typed || ''; }
      else { va = parseFloat(a[key]) || 0; vb = parseFloat(b[key]) || 0; return SORT_STATE.dir === 'asc' ? va - vb : vb - va; }
      var c = va.localeCompare(vb, undefined, { numeric: true, sensitivity: 'base' });
      return SORT_STATE.dir === 'asc' ? c : -c;
    });
    mapperHeatmapApplySerializedRows(rows);
    mapperHeatmapUpdateSortHeaders();
    mapperHeatmapScheduleRedraw();
  }

  function mapperHeatmapUpdateSortHeaders() {
    $('#mapper-heatmap-input-table th.mapper-core-sortable-head').each(function () {
      var col    = $(this).attr('data-mapper-heatmap-sort-col');
      var active = SORT_STATE.col === col;
      var dir    = active ? SORT_STATE.dir : null;
      $(this).attr('aria-sort', active ? (dir==='desc'?'descending':'ascending') : 'none');
      $(this).find('.mapper-core-sort-indicator').text(active ? (dir==='desc'?'▼':'▲') : '↕');
    });
  }

  // ── Demo ──────────────────────────────────────────────────────────────────
  function mapperHeatmapFillDemo() {
    suppressRedraw = true;
    mapperHeatmapDestroyAllRows();
    $('#mapper-heatmap-input-tbody').empty();
    MAPPER_HEATMAP_DEMO_ROWS.forEach(function (r) {
      mapperHeatmapAppendRow(true);
      var $tr = $('#mapper-heatmap-input-tbody tr').last();
      var rid = window.MAPPER_CORE_RESOLVE && window.MAPPER_CORE_RESOLVE[r.receptor.toUpperCase()];
      if (rid) mapperHeatmapSetResolved($tr, rid);
      else $tr.find('.mapper-core-in-receptor').val(r.receptor);
      [1,2,3,4,5].forEach(function (i) { $tr.find('.mapper-heatmap-val' + i).val(r.vals[i-1]); });
    });
    suppressRedraw = false;
    $('#mapper-heatmap-clear-rows').removeClass('mapper-core-clear-clean');
    mapperHeatmapEnsureTrailingBlankRow();
    mapperHeatmapSyncRemoveButtons();
    mapperHeatmapCompactReceptors();
    mapperHeatmapSyncFirstRowPlaceholder();
    // Apply Three_RWB preset for demo
    if ($.fn.select2) { try { $('#mapper-heatmap-color-style').val('Three_RWB').trigger('change'); } catch (e) {} }
    mapperHeatmapRedrawNow();
  }

  // ── Row find blank ────────────────────────────────────────────────────────
  function mapperHeatmapFindBlankRow() {
    var $found = $();
    $('#mapper-heatmap-input-tbody tr').each(function () { if (mapperHeatmapRowBlank($(this))) { $found = $(this); return false; } });
    return $found;
  }

  // ── Boot ──────────────────────────────────────────────────────────────────
  function mapperHeatmapBoot() {
    $('.mapper-core-booting').removeClass('mapper-core-booting');
    // Build converters
    mapperHeatmapBuildLabelConverter();
    try { mapperHeatmapExtendLabelConverterForSpecies(); } catch (e) { /* non-critical */ }

    // Extend autocomplete with non-human-only entries (receptors with no human ortholog)
    try {
      if (window.MAPPER_HEATMAP_SPECIES_DATA) {
        var _seenNhoStems = {};
        (window.MAPPER_HEATMAP_SPECIES_DATA.nonhuman_only || []).forEach(function (nho) {
          var stem = nho.stem;

          // Register MAPPER_CORE_ENTRY_META for every individual species entry so that
          // whichever species is later picked in the row's species dropdown, the cell
          // always shows "STEM (no human ortholog)" in all Receptor Names modes.
          var stemU = stem.toUpperCase();
          var TAG = '(no human ortholog)';
          var dispText = stemU + ' ' + TAG;
          var nhoHtml = '<span>' + $('<span>').text(stemU).html() + ' <em>' + TAG + '</em></span>';
          if (window.MAPPER_CORE_ENTRY_META && !window.MAPPER_CORE_ENTRY_META[nho.entry]) {
            window.MAPPER_CORE_ENTRY_META[nho.entry] = {
              name_html:  nhoHtml,
              name_plain: dispText,
              uniprot:    dispText
            };
          }

          if (_seenNhoStems[stem]) return;
          _seenNhoStems[stem] = true;
          if (window.MAPPER_CORE_RESOLVE && (window.MAPPER_CORE_RESOLVE[stem.toUpperCase()] ||
              window.MAPPER_CORE_RESOLVE[(stem + '_human').toUpperCase()])) return;
          var entryText = (nho.common || stem) + ' (no human ortholog)';
          if (window.receptorSelect2Data) {
            window.receptorSelect2Data.push({
              id: nho.entry, text: entryText, name_plain: entryText,
              name_html: $('<span>').text(entryText).html(),
              search_text: (stem + ' ' + (nho.common || '') + ' ' + (nho.latin || '')).toUpperCase()
            });
          }
          if (window.MAPPER_CORE_RESOLVE) {
            window.MAPPER_CORE_RESOLVE[stem.toUpperCase()] = nho.entry;
            window.MAPPER_CORE_RESOLVE[nho.entry.toUpperCase()] = nho.entry;
          }
        });
      }
    } catch (e) { /* non-critical */ }

    // Init HeatmapDataStyling
    if (typeof heatmap_DataStyling === 'function') {
      HeatmapDataStyling = heatmap_DataStyling();
    } else {
      HeatmapDataStyling = { Number_of_colors: 'Three', min_color: '#1a80bb', middle_color: '#ffffff', max_color: '#a00000',
        rotation: 90, label_position: 'Bottom', label_fontsize: 12, receptor_fontsize: 12,
        datalabels: true, data_border: true, data_fontsize: 12, legend_label: 'Value intensity', LabelType: 'IUPHAR' };
    }

    // Set up table
    $('#mapper-heatmap-input-tbody').empty();
    mapperHeatmapAppendRow();
    mapperHeatmapUpdateSortHeaders();
    mapperHeatmapSyncFirstRowPlaceholder();

    // Init colour pickers
    mapperHeatmapInitColorPickers();

    // Legend label inputs
    var lblMap = [['V1','Value1'],['V2','Value2'],['V3','Value3'],['V4','Value4'],['V5','Value5']];
    lblMap.forEach(function (m, idx) {
      $('#mapper-heatmap-lbl-' + m[0]).val(LabelXConverter[m[1]]).on('input', function () {
        LabelXConverter[m[1]] = $(this).val() || ('Dataset ' + (idx + 1));
        mapperHeatmapScheduleRedraw();
      });
    });

    // Receptor names buttons
    $('.mapper-heatmap-leaf-btn').on('click.mapperHeatmap', function () {
      var v = ($(this).attr('data-value') || 'IUPHAR').trim();
      HeatmapDataStyling.LabelType = v;
      $('.mapper-heatmap-leaf-btn').each(function () {
        var ok = ($(this).attr('data-value') || '') === v;
        $(this).toggleClass('btn-primary', ok).toggleClass('btn-outline-primary', !ok);
      });
      // Refresh resolved receptor names in table
      $('#mapper-heatmap-input-tbody tr').each(function () {
        var sid = ($( this).find('.mapper-core-receptor-entry').val() || '').trim();
        if (!sid) return;
        var $view = $(this).find('.mapper-core-receptor-html-view');
        if ($view.length && $view.is(':visible')) $view.html(mapperHeatmapResolvedDisplay(sid));
      });
      mapperHeatmapScheduleRedraw();
    });

    // Layout controls
    $('#mapper-heatmap-toggle-orient').on('click.mapperHeatmap', function () {
      HeatmapDataStyling.rotation = HeatmapDataStyling.rotation === 90 ? 0 : 90;
      $(this).text(HeatmapDataStyling.rotation === 90 ? 'Vertical' : 'Horizontal');
      mapperHeatmapScheduleRedraw();
    });
    $('#mapper-heatmap-toggle-pos').on('click.mapperHeatmap', function () {
      HeatmapDataStyling.label_position = HeatmapDataStyling.label_position === 'Bottom' ? 'Top' : 'Bottom';
      $(this).text(HeatmapDataStyling.label_position);
      mapperHeatmapScheduleRedraw();
    });
    $('#mapper-heatmap-label-fontsize').on('input', function () {
      $('#mapper-heatmap-label-fontsize-val').text($(this).val());
      HeatmapDataStyling.label_fontsize = $(this).val();
      mapperHeatmapScheduleRedraw();
    });
    $('#mapper-heatmap-receptor-fontsize').on('input', function () {
      $('#mapper-heatmap-receptor-fontsize-val').text($(this).val());
      HeatmapDataStyling.receptor_fontsize = $(this).val();
      mapperHeatmapScheduleRedraw();
    });
    $('#mapper-heatmap-toggle-datalabels').on('click.mapperHeatmap', function () {
      HeatmapDataStyling.datalabels = !HeatmapDataStyling.datalabels;
      $(this).text(HeatmapDataStyling.datalabels ? 'Shown' : 'Hidden')
             .toggleClass('btn-success', HeatmapDataStyling.datalabels)
             .toggleClass('btn-danger',  !HeatmapDataStyling.datalabels);
      mapperHeatmapScheduleRedraw();
    });
    $('#mapper-heatmap-toggle-border').on('click.mapperHeatmap', function () {
      HeatmapDataStyling.data_border = !HeatmapDataStyling.data_border;
      $(this).text(HeatmapDataStyling.data_border ? 'Shown' : 'Hidden')
             .toggleClass('btn-success', HeatmapDataStyling.data_border)
             .toggleClass('btn-danger',  !HeatmapDataStyling.data_border);
      mapperHeatmapScheduleRedraw();
    });
    $('#mapper-heatmap-data-fontsize').on('input', function () {
      $('#mapper-heatmap-data-fontsize-val').text($(this).val());
      HeatmapDataStyling.data_fontsize = $(this).val();
      mapperHeatmapScheduleRedraw();
    });

    // Demo
    $('#mapper-heatmap-demo-rows').on('click.mapperHeatmap', function () { mapperHeatmapFillDemo(); });

    // Clear
    $('#mapper-heatmap-clear-whole').on('click.mapperHeatmap', function(e) {
      e.preventDefault();
      mapperHeatmapDestroyAllRows();
      $('#mapper-heatmap-input-tbody').empty();
      mapperHeatmapAppendRow();
      mapperHeatmapCompactReceptors();
      mapperHeatmapSyncFirstRowPlaceholder();
      $('#mapper-heatmap-clear-rows').addClass('mapper-core-clear-clean').blur();
      mapperHeatmapRedrawNow();
    });
    $('.mapper-core-clear-menu').on('click.mapperHeatmap', '[data-heatmap-col]', function(e) {
      e.preventDefault();
      var col = $(this).data('heatmap-col');
      $('#mapper-heatmap-input-tbody tr').each(function() {
        $(this).find('.mapper-heatmap-' + col).val('').trigger('input');
      });
      $('#mapper-heatmap-clear-rows').removeClass('mapper-core-clear-clean');
      mapperHeatmapScheduleRedraw();
    });

    // Remove row
    $('#mapper-heatmap-input-table').on('click.mapperHeatmapRemove', '.mapper-core-remove-row', function () {
      var $tr = $(this).closest('tr');
      if (mapperHeatmapRowBlank($tr) && $tr.is(':last-child')) return;
      mapperHeatmapDestroyRowAc($tr); $tr.remove();
      $('#mapper-heatmap-clear-rows').removeClass('mapper-core-clear-clean');
      mapperHeatmapEnsureTrailingBlankRow(); mapperHeatmapCompactReceptors(); mapperHeatmapScheduleRedraw();
    });

    // Click resolved chip to re-enter edit mode
    $(document)
      .off('click.mapperHeatmapHtmlEdit', '#mapper-heatmap-input-tbody .mapper-core-receptor-html-view')
      .on('click.mapperHeatmapHtmlEdit', '#mapper-heatmap-input-tbody .mapper-core-receptor-html-view', function () {
        var $tr  = $(this).closest('tr');
        mapperHeatmapDestroyAc($tr.find('.mapper-core-in-receptor'));
        var hid  = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
        $(this).hide().empty();
        var $inp2 = $('<textarea class="form-control input-sm mapper-core-in-receptor" rows="1" autocomplete="off" spellcheck="false"></textarea>');
        var seed  = mapperHeatmapEditSeed(hid);
        $inp2.val(seed);
        if (!seed) { $tr.find('.mapper-core-receptor-input-wrap').addClass('is-empty'); }
        $tr.find('.mapper-core-receptor-input-wrap .mapper-core-in-receptor').remove();
        $tr.find('.mapper-core-receptor-input-wrap').prepend($inp2);
        $tr.find('.mapper-core-receptor-entry').val('');
        mapperHeatmapBindAc($inp2);
        $inp2.show().focus();
        window.setTimeout(function () {
          var t = ($inp2.val() || '').trim();
          if (t.length >= 1 && $inp2.data('ui-autocomplete')) { $inp2.autocomplete('search', t); }
        }, 0);
        mapperHeatmapSyncClearBtn($tr);
        mapperHeatmapSyncRemoveButtons();
        mapperHeatmapCompactReceptors();
        mapperHeatmapEnsureTrailingBlankRow();
        mapperHeatmapScheduleRedraw();
      });

    // Sort headers
    $('#mapper-heatmap-input-table').on('click.mapperHeatmapSort', 'th.mapper-core-sortable-head', function () {
      mapperHeatmapSortRows($(this).attr('data-mapper-heatmap-sort-col'));
    });

    // Input changes
    $(document).on('input.mapperHeatmap blur.mapperHeatmap change.mapperHeatmap',
      '#mapper-heatmap-input-tbody input, #mapper-heatmap-input-tbody textarea',
      function () {
        var $tr = $(this).closest('tr');
        mapperHeatmapSyncClearBtn($tr);
        $('#mapper-heatmap-clear-rows').removeClass('mapper-core-clear-clean');
        mapperHeatmapEnsureTrailingBlankRow();
        mapperHeatmapSyncRemoveButtons();
        mapperHeatmapCompactReceptors();
        mapperHeatmapSyncFirstRowPlaceholder();
        mapperHeatmapScheduleRedraw();
      }
    );

    // Paste
    $('#mapper-heatmap-input-table').on('paste.mapperHeatmapPaste', function (ePz) {
      var ev   = ePz.originalEvent || ePz;
      var text = ev.clipboardData ? ev.clipboardData.getData('text/plain') : '';
      if (!text || text.indexOf('\t') === -1) return;
      ePz.preventDefault();
      $('#mapper-heatmap-clear-rows').removeClass('mapper-core-clear-clean');
      suppressRedraw = true;
      text.split(/\r?\n/).forEach(function (ln) {
        if (!ln.trim()) return;
        var parts = ln.split('\t');
        var rawR  = (parts[0] || '').trim();
        var $tr   = mapperHeatmapFindBlankRow();
        if (!$tr.length) { mapperHeatmapAppendRow(true); $tr = $('#mapper-heatmap-input-tbody tr').last(); }
        $tr.find('.mapper-core-in-receptor').val(rawR);
        var rid = window.MAPPER_CORE_RESOLVE && window.MAPPER_CORE_RESOLVE[rawR.toUpperCase()];
        if (rid) mapperHeatmapSetResolved($tr, rid);
        [1,2,3,4,5].forEach(function (i) { if (parts[i]) $tr.find('.mapper-heatmap-val' + i).val((parts[i]||'').trim()); });
      });
      suppressRedraw = false;
      mapperHeatmapEnsureTrailingBlankRow();
      mapperHeatmapSyncRemoveButtons(); mapperHeatmapCompactReceptors();
      mapperHeatmapSyncFirstRowPlaceholder(); mapperHeatmapScheduleRedraw();
    });

    // Receptor lookup modal
    if (typeof window.mapperCoreInitGpcromePickerModal === 'function') {
      window.mapperCoreInitGpcromePickerModal({
        pickerRows: $.isArray(window.MAPPER_CORE_GPCROME_PICKER_ROWS) ? window.MAPPER_CORE_GPCROME_PICKER_ROWS : [],
        maxRows: MAPPER_HEATMAP_MAX_ROWS,
        onAdd: function (entryIds, meta) {
          suppressRedraw = true;
          // Heatmap has no category/text mode — assigned values are always sequential numbers.
          var numberMap = (meta && meta.groupNameById)
            ? window.mapperCoreBuildSequentialNumberMap(meta.groupNameById)
            : null;
          (entryIds || []).forEach(function (id) {
            var $row = mapperHeatmapFindBlankRow();
            if (!$row.length) { mapperHeatmapAppendRow(true); $row = $('#mapper-heatmap-input-tbody tr').last(); }
            mapperHeatmapSetResolved($row, id);
            if (numberMap && numberMap[id] != null) {
              $row.find('.mapper-heatmap-val1').val(String(numberMap[id])).trigger('input');
            }
          });
          suppressRedraw = false;
          $('#mapper-heatmap-clear-rows').removeClass('mapper-core-clear-clean');
          mapperHeatmapEnsureTrailingBlankRow();
          mapperHeatmapCompactReceptors();
          mapperHeatmapScheduleRedraw();
        }
      });
    }

    // Species toggle
    $('#mapper-heatmap-species-toggle').on('click.mapperHeatmap', function () {
      mapperHeatmapSpeciesActive = !mapperHeatmapSpeciesActive;
      $(this).toggleClass('btn-primary', mapperHeatmapSpeciesActive)
             .toggleClass('btn-default', !mapperHeatmapSpeciesActive);
      mapperHeatmapApplySpeciesPanelClass();
      $('#mapper-heatmap-input-tbody .mapper-heatmap-species-cell').toggle(mapperHeatmapSpeciesActive);
      $('.mapper-heatmap-species-h').toggle(mapperHeatmapSpeciesActive);
      $('#mapper-heatmap-species-names-wrap').toggle(mapperHeatmapSpeciesActive);
      mapperHeatmapExtendLabelConverterForSpecies();
      if (mapperHeatmapSpeciesActive) mapperHeatmapInitAllSpeciesDropdowns();
      mapperHeatmapScheduleRedraw();
    });

    // Species dropdown change
    $(document).on('change.mapperHeatmapSpecies',
      '#mapper-heatmap-input-tbody .mapper-heatmap-species-select',
      function () { mapperHeatmapScheduleRedraw(); }
    );

    // Species name format buttons — use direct binding (document delegation is blocked by
    // the stopPropagation on .dropdown-menu in Mapper_Heatmap.html)
    $('.mapper-heatmap-species-name-btn').off('click.mapperHeatmapSpeciesName').on('click.mapperHeatmapSpeciesName', function (eSpec) {
      eSpec.preventDefault();
      mapperHeatmapSpeciesNameFormat = $(this).data('value');
      $('.mapper-heatmap-species-name-btn').each(function () {
        var ok = $(this).data('value') === mapperHeatmapSpeciesNameFormat;
        $(this).toggleClass('btn-primary', ok).toggleClass('btn-outline-primary', !ok);
      });
      mapperHeatmapExtendLabelConverterForSpecies();
      mapperHeatmapInitAllSpeciesDropdowns();
      mapperHeatmapScheduleRedraw();
    });

    mapperHeatmapShowPlaceholder();
  }

  // ── Entry point ───────────────────────────────────────────────────────────
  $(document).ready(function () { mapperHeatmapBoot(); });

})(jQuery);

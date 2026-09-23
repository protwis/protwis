/**
 * Mapper 2.0 — List mapper page.
 * Left panel: receptor input table (numerical: Val1-4, categorical: Category+Swatch).
 * Right panel: list plot rendered via datamapper.js functions.
 */
(function ($) {
  'use strict';

  // ── Constants ────────────────────────────────────────────────────────────────
  var MAPPER_LIST_MAX_ROWS = 500;
  var MAPPER_LIST_PLACEHOLDER_NUMERIC = 'Paste in 1–5 columns of data\nor type';
  var MAPPER_LIST_PLACEHOLDER_CATMODE = 'Paste in 1–2 columns of data\nor type';
  var DEBOUNCE_MS = 200;

  // ── Module state ─────────────────────────────────────────────────────────────
  var suppressRedraw = false;
  var mapperListLastLabelSet = '';

  // ── Mode snapshots (mirrors tree) ────────────────────────────────────────────
  var MAPPER_LIST_MODE_SNAPSHOTS = { numeric: null, categorical: null };

  // Pre-defined demo colours for categorical labels
  var DEMO_CAT_COLORS = {
    'Serotonergic': '#1f77b4',
    'Chemokine':    '#ff7f0e',
    'Cholinergic':  '#2ca02c',
    'Adrenergic':   '#d62728',
    'Adhesion':     '#9467bd'
  };

  window.mapperListInputMode = 'numeric';
  window.MAPPER_CORE_LABEL_COLORS  = window.MAPPER_CORE_LABEL_COLORS  || {};
  window.MAPPER_CORE_LABEL_ENABLED = window.MAPPER_CORE_LABEL_ENABLED || {};

  var mapperListSpeciesActive = false;
  var mapperListSpeciesNameFormat = 'common';

  // ── Persistent config ────────────────────────────────────────────────────────
  /** Receptor info from Django: entry_name → {class, ligandtype, family, name_plain, gene, uniprot} */
  var ReceptorInfo = {};

  /** Label conversion dicts built at boot (same format the list plot JS expects). */
  var LabelConversionDicts = {
    UniProt_to_IUPHAR_converter: {},
    IUPHAR_to_UniProt_converter: {},
    UniProt_to_Gene_converter:   {},
    IUPHAR_to_Gene_converter:    {}
  };

  var LabelNames = 'Protein'; // 'Protein' | 'Gene' | 'UniProt'

  var Layout_dict = { columns: 2, col_max_label: 'Auto', Col_break_number: 20,
                       legend_top_space: 50, legend_position: 'Top',
                       skip_gradient_bars: true }; // use createListLegendBars instead

  var Styling_Option_dict = {
    Class:          { Bold: true,  Italic: false, Underline: true,  Fontsize: '20px', Color: '#000' },
    LigandType:     { Bold: false, Italic: false, Underline: false, Fontsize: '18px', Color: '#000' },
    ReceptorFamily: { Bold: false, Italic: true,  Underline: false, Fontsize: '16px', Color: '#000' },
    Receptor:       { Bold: false, Italic: false, Underline: false, Fontsize: '14px', Color: '#000' }
  };

  // Colour config for each numerical column (persists across re-renders).
  var ColorConfig = {
    Col1: { style: 'One', color1: '#ffffff', colorMid: '#ffffff', color2: '#707070' },
    Col2: { style: 'One', color1: '#ffffff', colorMid: '#ffffff', color2: '#0000ff' },
    Col3: { style: 'One', color1: '#ffffff', colorMid: '#ffffff', color2: '#ff0000' },
    Col4: { style: 'One', color1: '#ffffff', colorMid: '#ffffff', color2: '#008000' }
  };

  var Barlabels = { Col1: 'Dot 1', Col2: 'Dot 2', Col3: 'Dot 3', Col4: 'Dot 4' };

  var ShowLegend = true;
  var Textlegend_styling = {
    layoutMode: 'row', columns: 2, sortDirection: 'Vertically',
    TreeLegendPosition: 'Top', Fontsize: '14px'
  };

  var Data_styling = {}; // rebuilt each render from initializeDataStyling()

  // Demo rows
  var MAPPER_LIST_DEMO_ROWS = [
    { receptor: '5HT1A', vals: [1,  73, -10.1, -25], text: 'Serotonergic' },
    { receptor: '5HT1B', vals: [2,  72, -10.0, -24], text: 'Serotonergic' },
    { receptor: '5HT2A', vals: [6,  68,  -9.6, -20], text: 'Serotonergic' },
    { receptor: 'ACKR1', vals: [13, 61,  -8.9, -13], text: 'Chemokine' },
    { receptor: 'ACKR2', vals: [14, 60,  -8.8, -12], text: 'Chemokine' },
    { receptor: 'ACM1',  vals: [17, 57,  -8.5,  -9], text: 'Cholinergic' },
    { receptor: 'ACM2',  vals: [18, 56,  -8.4,  -8], text: 'Cholinergic' },
    { receptor: 'ADA1A', vals: [23, 51,  -7.9,  -3], text: 'Adrenergic' },
    { receptor: 'ADRB1', vals: [29, 45,  -7.3,   3], text: 'Adrenergic' },
    { receptor: 'ADRB2', vals: [30, 44,  -7.2,   4], text: 'Adrenergic' },
    { receptor: 'ADGRA1', vals: [31, 43, -7.1,   5], text: 'Adhesion' },
    { receptor: 'ADGRA2', vals: [32, 42, -7.0,   6], text: 'Adhesion' }
  ];

  // ── Helpers ──────────────────────────────────────────────────────────────────
  function mapperListIsTextMode() { return window.mapperListInputMode === 'text'; }

  // FNV1a32 hash + default label colour moved to MapperPageCore (mapper_page_core.js).
  function mapperListDefaultColorForLabel(lbl) {
    return MapperPageCore.defaultColorForLabel(lbl);
  }

  function parseNumLoose(s) {
    if (s == null || String(s).trim() === '') return null;
    var x = parseFloat(String(s).trim().replace(',', '.'));
    return isFinite(x) ? x : null;
  }

  // ── Label conversion ─────────────────────────────────────────────────────────
  /**
   * Build from MAPPER_CORE_ENTRY_META (always populated from receptor_select2_json),
   * mirroring exactly what the wheel and tree do for their name-conversion dicts.
   * The list plot's add_text() uses IUPHAR_to_Gene_converter[name_plain] and
   * IUPHAR_to_UniProt_converter[name_plain], so keys must be name_plain values.
   */
  function mapperListBuildLabelConversion() {
    LabelConversionDicts = {
      UniProt_to_IUPHAR_converter: {},
      IUPHAR_to_UniProt_converter: {},
      UniProt_to_Gene_converter:   {},
      IUPHAR_to_Gene_converter:    {}
    };
    var meta = window.MAPPER_CORE_ENTRY_META || {};
    Object.keys(meta).forEach(function (entry_name) {
      var m      = meta[entry_name];
      var iuphar = m.name_plain || '';   // IUPHAR display name (markup-stripped)
      var gene   = m.gene       || '';
      var uni    = m.uniprot    || '';   // short UniProt code, e.g. "HTR1A"
      if (!iuphar) return;
      LabelConversionDicts.UniProt_to_IUPHAR_converter[entry_name] = iuphar;
      LabelConversionDicts.IUPHAR_to_UniProt_converter[iuphar]     = entry_name;
      if (gene) {
        LabelConversionDicts.UniProt_to_Gene_converter[entry_name] = gene;
        LabelConversionDicts.IUPHAR_to_Gene_converter[iuphar]      = gene;
      }
      if (uni) {
        // Also index by short UniProt (already handled via entry_name path above for full form)
        LabelConversionDicts.UniProt_to_Gene_converter[uni] = gene;
      }
    });
  }

  // ── Species helpers ───────────────────────────────────────────────────────────
  function mapperListSpeciesLabel(o) {
    if (mapperListSpeciesNameFormat === 'latin') return o.latin || o.common || '';
    if (mapperListSpeciesNameFormat === 'both') {
      var c = (o.common || '').trim(), l = (o.latin || '').trim();
      return (c && l && c !== l) ? c + ' (' + l + ')' : c || l;
    }
    return o.common || o.latin || '';
  }

  function mapperListFormatSpeciesTag(o) {
    var c = (o.common || '').trim(), l = (o.latin || '').trim();
    if (mapperListSpeciesNameFormat === 'latin') return l ? '(' + l + ')' : '';
    if (mapperListSpeciesNameFormat === 'both') {
      if (c && l && c !== l) return '(' + c + ', ' + l + ')';
      return (c || l) ? '(' + (c || l) + ')' : '';
    }
    return (c || l) ? '(' + (c || l) + ')' : '';
  }

  function mapperListSpeciesLabelForEntry(specEntry) {
    var sd = window.MAPPER_LIST_SPECIES_DATA;
    if (!sd) return '';
    var stem = String(specEntry || '').split('_')[0];
    var pools = (sd.by_stem[stem] || []).concat(sd.nonhuman_only || []);
    for (var i = 0; i < pools.length; i++) {
      if (pools[i].entry === specEntry) return mapperListSpeciesLabel(pools[i]);
    }
    return '';
  }

  function mapperListPopulateSpeciesSelect($tr, entryId, preserveSelection) {
    var $sel = $tr.find('.mapper-list-species-select');
    if (!$sel.length || !entryId || !window.MAPPER_LIST_SPECIES_DATA) return;
    var prevVal = preserveSelection ? ($sel.val() || '').trim() : '';
    if ($sel.data('select2')) { try { $sel.select2('destroy'); } catch (e) {} }
    $sel.empty();
    var sd   = window.MAPPER_LIST_SPECIES_DATA;
    var stem = String(entryId).split('_')[0];
    var orthologs = (sd.by_stem[stem] || []).slice();
    if (!orthologs.length) {
      (sd.nonhuman_only || []).forEach(function (nho) { if (nho.stem === stem) orthologs.push(nho); });
    }
    orthologs.forEach(function (o) {
      $('<option>').val(o.entry).text(mapperListSpeciesLabel(o)).appendTo($sel);
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

  function mapperListInitAllSpeciesDropdowns() {
    $('#mapper-list-input-tbody tr').each(function () {
      var entry = ($(this).find('.mapper-core-receptor-entry').val() || '').trim();
      if (entry) mapperListPopulateSpeciesSelect($(this), entry, true);
    });
  }

  function mapperListApplyLeftPanelWidth() {
    var text = window.mapperListInputMode === 'text';
    var speciesExtra = mapperListSpeciesActive ? 88 : 0;
    $('.mapper-core-wheel-wrap.mapper-list-page .mapper-core-wheel-left').css(
      'width', text ? (400 + speciesExtra) + 'px' : ''
    );
    $('.mapper-core-wheel-wrap.mapper-list-page').toggleClass('mapper-list-species-active', mapperListSpeciesActive);
  }

  function mapperListApplySpeciesPanelClass() {
    mapperListApplyLeftPanelWidth();
  }

  function mapperListExtendLabelConversionForSpecies() {
    var sd = window.MAPPER_LIST_SPECIES_DATA;
    if (!sd || !sd.by_stem) return;
    function addEntry(o) {
      var stem = o.stem;
      var humanOrthologs = (sd.by_stem[stem] || []).filter(function (x) { return x.is_human; });
      var humanEntry = humanOrthologs.length ? humanOrthologs[0].entry : null;
      var baseMeta = humanEntry && window.MAPPER_CORE_ENTRY_META && window.MAPPER_CORE_ENTRY_META[humanEntry];
      var baseIUPHAR = (baseMeta && baseMeta.name_plain) || stem.toUpperCase();
      var gene = (baseMeta && baseMeta.gene) || '';
      var tag  = mapperListSpeciesActive ? mapperListFormatSpeciesTag(o) : (o.is_human ? '' : mapperListFormatSpeciesTag(o));
      var fullName = tag ? baseIUPHAR + ' ' + tag : baseIUPHAR;
      LabelConversionDicts.UniProt_to_IUPHAR_converter[o.entry] = fullName;
      LabelConversionDicts.IUPHAR_to_UniProt_converter[fullName] = o.entry;
      if (gene) {
        var geneLbl = tag ? gene + ' ' + tag : gene;
        LabelConversionDicts.UniProt_to_Gene_converter[o.entry]  = geneLbl;
        LabelConversionDicts.IUPHAR_to_Gene_converter[fullName]  = geneLbl;
      }
    }
    Object.keys(sd.by_stem).forEach(function (stem) { (sd.by_stem[stem] || []).forEach(addEntry); });
    (sd.nonhuman_only || []).forEach(addEntry);
  }

  // ── Data construction ─────────────────────────────────────────────────────────
  /**
   * Builds the two data structures the list plot JS needs:
   *  listdata      – hierarchical: Class → LigandType → Family → iuphar_name → {Value1,…}
   *  list_data_wow – flat:  iuphar_name → {Value1,…}
   * Keys are always name_plain (IUPHAR); Label_conversion_dicts handles display.
   */
  function mapperListBuildData() {
    var listdata     = {};
    var list_data_wow = {};

    $('#mapper-list-input-tbody tr').each(function () {
      var $tr    = $(this);
      var entry  = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
      // Fallback: try to resolve via MAPPER_CORE_RESOLVE if entry hidden-input is empty
      if (!entry) {
        var ta       = ($tr.find('.mapper-core-in-receptor').val() || '').trim();
        var unmatched = ($tr.data('mapperCoreUnmatchedRaw') || '');
        var up = String(ta || unmatched).trim().toUpperCase();
        if (up && window.MAPPER_CORE_RESOLVE) entry = window.MAPPER_CORE_RESOLVE[up] || '';
      }
      if (!entry) return;
      var info = ReceptorInfo[entry];
      if (!info) {
        // Non-human entries (e.g. 5ht5b_mouse) may not be in ReceptorInfo — try human stem fallback
        var stemFallback = entry.split('_')[0] + '_human';
        info = ReceptorInfo[stemFallback] || null;
      }
      if (!info) return;

      var cls     = info['class']   || 'Other';
      var ligtype = info.ligandtype || 'Other';
      var family  = info.family     || 'Other';
      // nameKey: species-aware when species mode is active
      var metaEntry = (window.MAPPER_CORE_ENTRY_META && window.MAPPER_CORE_ENTRY_META[entry]) || {};
      var baseNamePlain = metaEntry.name_plain || info.name_plain || entry;
      var nameKey;
      if (mapperListSpeciesActive) {
        var specVal = ($tr.find('.mapper-list-species-select').val() || '').trim();
        var dataEntryId = specVal || entry;
        nameKey = LabelConversionDicts.UniProt_to_IUPHAR_converter[dataEntryId] || baseNamePlain;
      } else {
        nameKey = baseNamePlain;
      }

      // Validate values BEFORE touching listdata — family level is an ARRAY of names
      // (Data_resorter expects familyObj.receptors.forEach, so must be an array)
      var dataEntry;
      if (mapperListIsTextMode()) {
        var cat = ($tr.find('.mapper-list-cat').val() || '').trim();
        if (!cat) return;
        if (!window.MAPPER_CORE_LABEL_COLORS[cat]) window.MAPPER_CORE_LABEL_COLORS[cat] = mapperListDefaultColorForLabel(cat);
        var en       = window.MAPPER_CORE_LABEL_ENABLED[cat] !== false;
        var colorVal = en ? window.MAPPER_CORE_LABEL_COLORS[cat] : '#ffffff';
        dataEntry = { Value1: cat, ColorValue: colorVal };
      } else {
        var v1 = parseNumLoose($tr.find('.mapper-list-val1').val());
        var v2 = parseNumLoose($tr.find('.mapper-list-val2').val());
        var v3 = parseNumLoose($tr.find('.mapper-list-val3').val());
        var v4 = parseNumLoose($tr.find('.mapper-list-val4').val());
        if (v1 == null && v2 == null && v3 == null && v4 == null) return;
        dataEntry = {};
        if (v1 != null) dataEntry.Value1 = v1;
        if (v2 != null) dataEntry.Value2 = v2;
        if (v3 != null) dataEntry.Value3 = v3;
        if (v4 != null) dataEntry.Value4 = v4;
      }

      // Now commit to listdata (array) and list_data_wow
      if (!listdata[cls]) listdata[cls] = {};
      if (!listdata[cls][ligtype]) listdata[cls][ligtype] = {};
      if (!listdata[cls][ligtype][family]) listdata[cls][ligtype][family] = [];
      if (listdata[cls][ligtype][family].indexOf(nameKey) === -1) {
        listdata[cls][ligtype][family].push(nameKey);
      }
      list_data_wow[nameKey] = dataEntry;
    });
    return { listdata: listdata, list_data_wow: list_data_wow };
  }

  // ── Apply colour config to Data_styling ────────────────────────────────────
  function mapperListApplyColorConfig() {
    ['Col1', 'Col2', 'Col3', 'Col4'].forEach(function (col) {
      if (!Data_styling[col] || Data_styling[col].Data !== 'Yes') return;
      var cfg = ColorConfig[col];
      Data_styling[col].Data_color1           = cfg.color1;
      Data_styling[col].Data_color2           = cfg.color2;
      Data_styling[col].data_color_complexity = (cfg.style === 'One') ? 'One' : (cfg.style === 'Two') ? 'Two' : 'Three';
    });
  }

  // ── Placeholder ───────────────────────────────────────────────────────────────
  function mapperListSyncFirstRowPlaceholder() {
    var hint = mapperListIsTextMode() ? MAPPER_LIST_PLACEHOLDER_CATMODE : MAPPER_LIST_PLACEHOLDER_NUMERIC;
    $('#mapper-list-input-tbody tr').each(function (idx) {
      var $inp = $(this).find('.mapper-core-in-receptor');
      if (!$inp.length) return;
      if (idx === 0) $inp.attr('placeholder', hint);
      else           $inp.removeAttr('placeholder');
    });
  }

  // ── Placeholder plot ──────────────────────────────────────────────────────────
  function mapperListShowPlaceholder() {
    var $host = $('#mapper-list-plot');
    $host.empty();
    $host.append(
      $('<div class="mapper-list-plot-placeholder">').append(
        $('<p class="mapper-list-plot-placeholder-title">').text('No receptors mapped yet'),
        $('<p class="mapper-list-plot-placeholder-hint">').text(
          'Add receptors in the left panel and enter values (Numerical) or category labels (Categorical) — the list renders automatically.'
        )
      )
    );
  }

  // ── Render ───────────────────────────────────────────────────────────────────
  function mapperListRedrawNow() {
    var built = mapperListBuildData();

    if (!Object.keys(built.listdata).length) {
      mapperListShowPlaceholder();
      return;
    }

    // Normalise hierarchy
    var modded;
    try {
      modded = Initialize_Data(built.listdata);
    } catch (e) {
      mapperListShowPlaceholder(); return;
    }

    // Data types & styling
    var dt = mapperListIsTextMode()
      ? { Col1: 'Discrete', Col2: 'Discrete', Col3: 'Discrete', Col4: 'Discrete' }
      : { Col1: 'Continuous', Col2: 'Continuous', Col3: 'Continuous', Col4: 'Continuous' };

    try {
      Data_styling = initializeDataStyling(built.list_data_wow, dt);
    } catch (e) {}
    mapperListApplyColorConfig();

    // Expose as globals — datamapper.js reads Data_styling and list_data_wow directly
    window.Data_styling   = Data_styling;
    window.list_data_wow  = built.list_data_wow;

    // Sort/flatten
    var result;
    try {
      result = Data_resorter(modded);
    } catch (e) {
      mapperListShowPlaceholder(); return;
    }
    var Data_array          = result.final_array;
    var Data_category_array = result.category_array;
    if (!Data_array || !Data_array.length) {
      mapperListShowPlaceholder(); return;
    }

    // Column break
    var cols = Layout_dict.columns;
    var breakN = Layout_dict.col_max_label === 'Auto'
      ? Math.ceil(Data_array.length / cols)
      : (Layout_dict.Col_break_number || Math.ceil(Data_array.length / cols));
    Layout_dict.Col_break_number = breakN;

    // Clear placeholder and any existing SVG before rendering
    $('#mapper-list-plot').empty();

    // Render labels
    try {
      RenderListPlot_Labels(
        Data_array, Data_category_array, 'mapper-list-plot',
        Styling_Option_dict, Layout_dict, LabelConversionDicts, LabelNames
      );
    } catch (e) {
      mapperListShowPlaceholder(); return;
    }

    // Calculate spacing for shapes
    var spacing;
    try {
      spacing = Calculate_dimension(
        Data_array, Data_category_array, breakN, cols,
        LabelConversionDicts, LabelNames, Styling_Option_dict
      );
    } catch (e) {
      spacing = {};
    }

    // Render shapes
    try {
      data_visualization(
        Data_array, Data_category_array, 'mapper-list-plot',
        Layout_dict, Data_styling, spacing, 30, 6, 14, Barlabels
      );
    } catch (e) {}

    // Gradient legend bars (numerical) — uses tree-style createListLegendBars
    if (!mapperListIsTextMode() && typeof createListLegendBars === 'function') {
      try {
        createListLegendBars(
          'mapper-list-plot',
          built.list_data_wow,
          ColorConfig,
          Barlabels,
          dt,
          Layout_dict.legend_position || 'Top',
          Layout_dict.legend_top_space
        );
      } catch (e) {}
    }

    // Categorical legend
    if (mapperListIsTextMode() && ShowLegend) {
      try {
        CreateTextLegend_list('mapper-list-plot', built.list_data_wow, Textlegend_styling);
      } catch (e) {}
    }

    mapperListSyncColorDataRows();
  }

  var _listRedrawDebounced = MapperPageCore.debounce(function () { mapperListRedrawNow(); }, DEBOUNCE_MS);
  function mapperListScheduleRedraw() {
    if (suppressRedraw) return;
    _listRedrawDebounced.schedule();
  }

  // Batches the full-table cosmetic syncs triggered by the two hottest, highest-frequency
  // actions (typing in a value cell, deleting a row) so they run once per debounce window
  // instead of once per keystroke/click at large row counts. Structural row-count upkeep
  // (mapperListEnsureTrailingBlankRow) stays synchronous — it's already bounded to the last
  // two rows, not a full-table scan. Every other call site of these sync functions
  // (restore, demo-fill, paste, clear) is untouched and keeps calling them directly.
  var _listCosmeticSyncDebounced = MapperPageCore.debounce(function () {
    mapperListRefreshSwatches();
    mapperListSyncCatHints();
    mapperListSyncRemoveButtons();
    mapperListCompactReceptors();
    mapperListSyncFirstRowPlaceholder();
  }, DEBOUNCE_MS);

  // Show / hide which color rows have data
  function mapperListSyncColorDataRows() {
    if (mapperListIsTextMode()) return;
    ['Col1','Col2','Col3','Col4'].forEach(function(col) {
      var hasData = Data_styling[col] && Data_styling[col].Data === 'Yes';
      // Always show all rows for the list mapper
    });
    var anyData = ['Col1','Col2','Col3','Col4'].some(function(c){ return Data_styling[c] && Data_styling[c].Data === 'Yes'; });
    $('#mapper-list-colors-no-data').toggle(!anyData);
  }

  // ── Input mode switch ─────────────────────────────────────────────────────────
  // ── Snapshot helpers ──────────────────────────────────────────────────────────
  function mapperListCaptureSnapshot() {
    var mode = mapperListIsTextMode() ? 'categorical' : 'numeric';
    var rows = mapperListSerializeRows();
    if (mode === 'categorical') {
      rows = rows.map(function (r) {
        return { entry: r.entry, typed: r.typed, unmatched: r.unmatched, invalid: r.invalid,
                 v1: '', v2: '', v3: '', v4: '', cat: r.cat, species: r.species };
      });
      MAPPER_LIST_MODE_SNAPSHOTS.categorical = {
        rows: rows,
        labelColors: $.extend({}, window.MAPPER_CORE_LABEL_COLORS || {})
      };
    } else {
      rows = rows.map(function (r) {
        return { entry: r.entry, typed: r.typed, unmatched: r.unmatched, invalid: r.invalid,
                 v1: r.v1, v2: r.v2, v3: r.v3, v4: r.v4, cat: '', species: r.species };
      });
      MAPPER_LIST_MODE_SNAPSHOTS.numeric = { rows: rows };
    }
  }

  function mapperListSeedOppositeMode(targetMode) {
    if (MAPPER_LIST_MODE_SNAPSHOTS[targetMode] != null) return;
    if (targetMode === 'categorical') {
      MAPPER_LIST_MODE_SNAPSHOTS.categorical = { rows: [], labelColors: {} };
    } else {
      MAPPER_LIST_MODE_SNAPSHOTS.numeric = { rows: [] };
    }
  }

  function mapperListSyncClearDropdown() {
    var isText = mapperListIsTextMode();
    $('#mapper-list-clear-mode-label').text(isText ? 'Categories' : 'Numbers');
    $('.mapper-list-clear-col-num').toggle(!isText);
    $('.mapper-list-clear-col-text').toggle(isText);
  }

  function mapperListSetInputMode(mode) {
    var nextMode = mode === 'text' ? 'text' : 'numeric';
    var prevMode = window.mapperListInputMode;
    if (prevMode === nextMode) return;

    // Save current state before switching (skip on initial boot)
    if (prevMode === 'text' || prevMode === 'numeric') {
      mapperListCaptureSnapshot();
    }

    window.mapperListInputMode = nextMode;
    var text = nextMode === 'text';

    // Seed the target snapshot from current if not yet visited
    mapperListSeedOppositeMode(text ? 'categorical' : 'numeric');

    // Restore rows from snapshot — always rebuild, even when the target mode has
    // never been visited (snap.rows is []), so a genuinely empty mode shows as empty
    // rather than leaving the previous mode's rows sitting in the DOM.
    var snap = MAPPER_LIST_MODE_SNAPSHOTS[text ? 'categorical' : 'numeric'] || {};
    if (text) {
      window.MAPPER_CORE_LABEL_COLORS = snap.labelColors ? $.extend({}, snap.labelColors) : {};
    }
    suppressRedraw = true;
    mapperListDestroyAllRows();
    $('#mapper-list-input-tbody').empty();
    mapperListApplySerializedRows(snap.rows || []);
    suppressRedraw = false;

    // Button states
    $('#mapper-list-mode-numeric-btn, #mapper-list-mode-labels-btn').each(function () {
      var isNum = $(this).attr('id') === 'mapper-list-mode-numeric-btn';
      var on = text ? !isNum : isNum;
      $(this).toggleClass('active', on).attr('aria-pressed', on ? 'true' : 'false');
    });

    // Table class
    var $tbl = $('#mapper-list-input-table');
    $tbl.toggleClass('mapper-core-text-mode', text);
    $tbl.find('.mapper-list-val1,.mapper-list-val2,.mapper-list-val3,.mapper-list-val4').prop('disabled', text);

    // Toolbar visibility
    $('#mapper-list-legend-labels-wrap').toggle(!text);
    $('#mapper-list-colors-numeric').toggle(!text);
    $('#mapper-list-colors-text').toggle(text);
    $('#mapper-list-legend-wrap').toggle(text);
    $('#mapper-list-colors-menu').toggleClass('mapper-list-colors-text-mode', text);

    // Shrink/grow left panel
    $('.mapper-core-wheel-wrap.mapper-list-page').toggleClass('mapper-list-catmode', text);
    mapperListApplyLeftPanelWidth();

    // Category hints
    mapperListSyncCatHints();
    mapperListRefreshSwatches();
    mapperListSyncColorsPanel();
    mapperListSyncFirstRowPlaceholder();
    mapperListScheduleRedraw();
    if (typeof mapperListSyncClearDropdown === 'function') { mapperListSyncClearDropdown(); }
  }

  // ── Category cell orange hints ────────────────────────────────────────────────
  function mapperListSyncCatHints() {
    var text = mapperListIsTextMode();
    $('#mapper-list-input-tbody tr').each(function () {
      var $tr = $(this);
      var $td = $tr.find('td.mapper-list-cat-cell');
      if (!$td.length) return;
      if (!text) { $td.removeClass('mapper-list-cat-hint'); return; }
      var hasRec = !!($tr.find('.mapper-core-receptor-entry').val() || '').trim();
      var hasCat = !!($tr.find('.mapper-list-cat').val() || '').trim();
      $td.toggleClass('mapper-list-cat-hint', hasRec && !hasCat);
    });
  }

  // ── Colour swatches (categorical) ─────────────────────────────────────────────
  function mapperListDestroySwatchSpectrum($picker) {
    if (!$picker || !$picker.length || !$.fn.spectrum) return;
    try { if ($picker.data('spectrum.id') != null || $picker.hasClass('sp-replaced')) $picker.spectrum('destroy'); } catch (e) {}
  }

  function mapperListApplyLabelColor(label, color, $activePicker) {
    var key = String(label || '').trim();
    if (!key || !color) return;
    window.MAPPER_CORE_LABEL_COLORS[key] = color;
    $('#mapper-list-input-tbody tr').each(function () {
      var $tr = $(this);
      if (($tr.find('.mapper-list-cat').val() || '').trim() !== key) return;
      var $sw = $tr.find('.mapper-core-row-color-picker');
      $sw.val(color).css('background-color', color);
      if ($activePicker && $activePicker.length && $sw[0] === $activePicker[0]) return;
      if ($.fn.spectrum && ($sw.data('spectrum.id') != null || $sw.hasClass('sp-replaced'))) {
        try { $sw.spectrum('set', color); } catch (e) {}
      }
    });
    // Sync panel picker
    var $pp = $('#mapper-list-label-color-pickers .mapper-core-lcat-spectrum[data-mapper-list-lcat-label="' + key + '"]');
    if ($pp.length && (!$activePicker || !$activePicker.length || $pp[0] !== $activePicker[0])) {
      if ($.fn.spectrum && ($pp.data('spectrum') || $pp.hasClass('sp-replaced'))) {
        try { $pp.spectrum('set', color); } catch (e) {}
      }
    }
    mapperListScheduleRedraw();
  }

  function mapperListEnsureSwatchSpectrum($picker, label, color) {
    if (!$picker.length || !$.fn.spectrum) { $picker.css('background-color', color || '#f5f5f5'); return; }
    if ($picker.data('spectrum.id') != null || $picker.hasClass('sp-replaced')) {
      if ($picker.attr('data-mapper-list-label') === label) {
        $picker.spectrum('set', color || '#f5f5f5'); return;
      }
      mapperListDestroySwatchSpectrum($picker);
    }
    $picker.attr('data-mapper-list-label', label || '').val(color || '#f5f5f5').css('background-color', color || '#f5f5f5');
    $picker.spectrum({
      color: color || '#f5f5f5', preferredFormat: 'hex', showInput: true,
      showPalette: true, showSelectionPalette: true, clickoutFiresChange: true,
      containerClassName: 'mapper-core-row-color-spectrum',
      replacerClassName: 'mapper-list-row-swatch-replacer',
      palette: [
        ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd'],
        ['#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf'],
        ['#000000','#666666','#aaaaaa','#ffffff']
      ],
      move:   function (t) { mapperListApplyLabelColor(label, t ? t.toHexString() : color, $picker); },
      change: function (t) { mapperListApplyLabelColor(label, t ? t.toHexString() : color, $picker); }
    });
  }

  function mapperListRefreshSwatches() {
    if (!mapperListIsTextMode()) {
      $('#mapper-list-input-tbody .mapper-core-row-color-picker').each(function () { mapperListDestroySwatchSpectrum($(this)); });
      $('#mapper-list-input-tbody .mapper-list-label-swatch').css('background-color', 'transparent');
      return;
    }
    $('#mapper-list-input-tbody tr').each(function () {
      var $sw  = $(this).find('.mapper-list-label-swatch');
      var cat  = ($(this).find('.mapper-list-cat').val() || '').trim();
      if (!cat) { mapperListDestroySwatchSpectrum($sw); $sw.css('background-color', '#f5f5f5'); return; }
      var hex = window.MAPPER_CORE_LABEL_COLORS[cat] || mapperListDefaultColorForLabel(cat);
      $sw.css('background-color', hex);
      mapperListEnsureSwatchSpectrum($sw, cat, hex);
    });
    mapperListSyncColorsPanel();
  }

  // ── Categorical colours panel ─────────────────────────────────────────────────
  function mapperListSyncColorsPanel() {
    var $mount = $('#mapper-list-label-color-pickers');
    if (!$mount.length) return;

    if (!mapperListIsTextMode()) {
      $mount.find('.mapper-core-lcat-spectrum').each(function () { try { if ($(this).data('spectrum')) $(this).spectrum('destroy'); } catch (e) {} });
      $mount.empty(); mapperListLastLabelSet = ''; return;
    }

    var seen = {}, labels = [];
    $('#mapper-list-input-tbody tr').each(function () {
      var lbl = ($(this).find('.mapper-list-cat').val() || '').trim();
      if (lbl && !seen[lbl]) { seen[lbl] = true; labels.push(lbl); }
    });
    labels.sort(function (a, b) { return a.localeCompare(b, undefined, { numeric: true, sensitivity: 'base' }); });

    var newKey = labels.join('\n');
    if (newKey === mapperListLastLabelSet) return;

    $mount.find('.mapper-core-lcat-spectrum').each(function () { try { if ($(this).data('spectrum')) $(this).spectrum('destroy'); } catch (e) {} });
    $mount.empty();
    mapperListLastLabelSet = newKey;

    if (!labels.length) {
      $mount.append($('<p style="color:#888; font-size:11px; text-align:center; margin:4px 0;">Enter categorical labels to configure colours.</p>'));
      return;
    }

    // Initialise colors/enabled
    labels.forEach(function (lbl) {
      if (!window.MAPPER_CORE_LABEL_COLORS[lbl]) window.MAPPER_CORE_LABEL_COLORS[lbl] = mapperListDefaultColorForLabel(lbl);
      if (!Object.prototype.hasOwnProperty.call(window.MAPPER_CORE_LABEL_ENABLED, lbl)) window.MAPPER_CORE_LABEL_ENABLED[lbl] = true;
    });

    function domId(ix, kind) { return 'mapper_list_l_' + kind + '_' + ix; }

    function updateMaster() {
      var allOn  = labels.every(function (l) { return window.MAPPER_CORE_LABEL_ENABLED[l] !== false; });
      var allOff = labels.every(function (l) { return window.MAPPER_CORE_LABEL_ENABLED[l] === false; });
      var mx = $('#mapper-list-label-cat-master')[0];
      if (!mx) return;
      mx.checked = !!allOn;
      mx.indeterminate = !allOn && !allOff && labels.length > 0;
    }

    function refreshRow(ix, lbl) {
      var $inp = $('#' + domId(ix, 'spe'));
      if (!$inp.length || !$inp.data('spectrum')) return;
      var en = window.MAPPER_CORE_LABEL_ENABLED[lbl] !== false;
      try { en ? $inp.spectrum('enable').css('opacity','1') : $inp.spectrum('disable').css('opacity','0.5'); } catch (e) {}
    }

    var $hdr = $('<div class="mapper-core-label-cat-head">').append(
      $('<input type="checkbox" id="mapper-list-label-cat-master" aria-label="Toggle all">'),
      $('<div class="mapper-core-label-cat-head-title">Categories</div>')
    );
    var $grid = $('<div class="color-grid">');
    $mount.append($hdr, $grid);
    updateMaster();

    $('#mapper-list-label-cat-master').on('change.mapperListColors', function () {
      var on = !!$(this).prop('checked');
      labels.forEach(function (lbl, ix) { window.MAPPER_CORE_LABEL_ENABLED[lbl] = on; refreshRow(ix, lbl); });
      mapperListRefreshSwatches(); mapperListScheduleRedraw();
    });

    labels.forEach(function (lbl, ix) {
      var chkId = domId(ix, 'chk'), spId = domId(ix, 'spe');
      var hex = window.MAPPER_CORE_LABEL_COLORS[lbl] || mapperListDefaultColorForLabel(lbl);
      window.MAPPER_CORE_LABEL_COLORS[lbl] = hex;
      var en = window.MAPPER_CORE_LABEL_ENABLED[lbl] !== false;

      var $spe = $('<input type="text">').addClass('mapper-core-lcat-spectrum form-control input-sm')
        .attr({ id: spId, 'data-mapper-list-lcat-label': lbl });

      $grid.append($('<div class="color-item">').append(
        $('<input type="checkbox">').attr('id', chkId).prop('checked', !!en),
        $('<label>').attr('for', chkId).addClass('color-label').text(lbl),
        $spe
      ));

      $('#' + chkId).on('change.mapperListColors', function () {
        window.MAPPER_CORE_LABEL_ENABLED[lbl] = !!$(this).prop('checked');
        refreshRow(ix, lbl); updateMaster();
        mapperListRefreshSwatches(); mapperListScheduleRedraw();
      });

      $spe.spectrum({
        color: hex, showPalette: true, showInput: true, showButtons: false, preferredFormat: 'hex',
        appendTo: '#mapper-list-colors-menu',
        containerClassName: 'mapper-core-lcat-sp-container',
        replacerClassName: 'mapper-core-lcat-replacer',
        palette: [
          ['#000','#FF0000','#00FF00','#0000FF','#FFFF00'],
          ['#FF00FF','#00FFFF','#FFFFFF','#C0C0C0','#808080']
        ],
        change: function (c) {
          if (window.MAPPER_CORE_LABEL_ENABLED[lbl] === false) return;
          mapperListApplyLabelColor(lbl, c && c.toHexString ? c.toHexString() : hex, $spe);
        }
      });
      if (!en) { try { $spe.spectrum('disable').css('opacity','0.5'); } catch (e) {} }
    });
  }

  // ── Numerical colour pickers ──────────────────────────────────────────────────
  var COLOR_PRESETS = {
    One:      { setup: 'One',   color1: '#ffffff', colorMid: '#ffffff', color2: '#707070' },
    Two:      { setup: 'Two',   color1: '#97a6c4', colorMid: '#ffffff', color2: '#384860' },
    Three_RWB:{ setup: 'Three', color1: '#a00000', colorMid: '#ffffff', color2: '#1a80bb' },
    Three_TWM:{ setup: 'Three', color1: '#298c8c', colorMid: '#ffffff', color2: '#800074' }
  };

  function formatColorSchemePicker(option) {
    if (!option.id) return option.text;
    var colorMap = {
      One:       [null, null, '#707070'],
      Two:       ['#97a6c4', null, '#384860'],
      Three_RWB: ['#a00000','#ffffff','#1a80bb'],
      Three_TWM: ['#298c8c','#ffffff','#800074']
    };
    var colors = colorMap[option.id];
    if (!colors) return option.text;
    return $('<span style="display:flex; align-items:center;">').append(
      $('<span style="min-width:50px; text-align:center; margin-right:5px;">').text(option.text),
      colors.map(function (c) {
        return $('<span style="display:inline-block; width:12px; height:12px; margin-left:4px; border:1px solid #ccc; border-radius:2px; background:' + (c||'transparent') + '; ' + (c?'':'opacity:0') + ';">');
      })
    );
  }

  function applyColStylePreset(col, styleName, skipRedraw) {
    var preset = COLOR_PRESETS[styleName];
    if (!preset) return;
    ColorConfig[col].style    = styleName;
    ColorConfig[col].color1   = preset.color1;
    ColorConfig[col].colorMid = preset.colorMid;
    ColorConfig[col].color2   = preset.color2;

    var isOne   = preset.setup === 'One';
    var isThree = preset.setup === 'Three';
    $('#mapper-list-cpicker-min-' + col).closest('.mapper-list-cpicker-min-wrap').css('visibility', isOne ? 'hidden' : 'visible');
    $('#mapper-list-cpicker-mid-' + col).closest('.mapper-list-cpicker-mid-wrap').css('visibility', isThree ? 'visible' : 'hidden');

    if ($.fn.spectrum) {
      try { $('#mapper-list-cpicker-min-' + col).spectrum('set', preset.color1); } catch (e) {}
      try { $('#mapper-list-cpicker-max-' + col).spectrum('set', preset.color2); } catch (e) {}
    }
    if (!skipRedraw) mapperListScheduleRedraw();
  }

  function mapperListInitColorPickers() {
    ['Col1','Col2','Col3','Col4'].forEach(function (col) {
      var cfg = ColorConfig[col];

      // Style select2
      $('#mapper-list-cstyle-' + col).select2({
        templateResult: formatColorSchemePicker, templateSelection: formatColorSchemePicker,
        width: 'resolve', dropdownParent: $('#mapper-list-colors-menu')
      }).on('change', function () { applyColStylePreset(col, $(this).val()); });

      // Min picker
      $('#mapper-list-cpicker-min-' + col).spectrum({
        color: cfg.color1, showInput: true, showPalette: false, preferredFormat: 'hex',
        change: function (c) { ColorConfig[col].color1 = c ? c.toHexString() : cfg.color1; mapperListScheduleRedraw(); },
        move:   function (c) { ColorConfig[col].color1 = c ? c.toHexString() : cfg.color1; mapperListScheduleRedraw(); }
      });
      // Mid picker
      $('#mapper-list-cpicker-mid-' + col).spectrum({
        color: cfg.colorMid, showInput: true, showPalette: false, preferredFormat: 'hex',
        change: function (c) { ColorConfig[col].colorMid = c ? c.toHexString() : cfg.colorMid; mapperListScheduleRedraw(); },
        move:   function (c) { ColorConfig[col].colorMid = c ? c.toHexString() : cfg.colorMid; mapperListScheduleRedraw(); }
      });
      // Max picker
      $('#mapper-list-cpicker-max-' + col).spectrum({
        color: cfg.color2, showInput: true, showPalette: false, preferredFormat: 'hex',
        change: function (c) { ColorConfig[col].color2 = c ? c.toHexString() : cfg.color2; mapperListScheduleRedraw(); },
        move:   function (c) { ColorConfig[col].color2 = c ? c.toHexString() : cfg.color2; mapperListScheduleRedraw(); }
      });

      // Apply default preset silently
      applyColStylePreset(col, 'One', true);
    });
  }

  // ── Label styling pickers ─────────────────────────────────────────────────────
  function mapperListInitLabelStylePickers() {
    var layers = ['Class','LigandType','ReceptorFamily','Receptor'];
    layers.forEach(function (layer) {
      $('#colorPicker_' + layer).spectrum({
        color: Styling_Option_dict[layer].Color, showInput: true, showPalette: true,
        preferredFormat: 'hex',
        palette: [['#000','#FF0000','#00FF00','#0000FF'],['#FFFFFF','#808080']],
        change: function (c) { Styling_Option_dict[layer].Color = c ? c.toHexString() : '#000'; mapperListScheduleRedraw(); },
        move:   function (c) { Styling_Option_dict[layer].Color = c ? c.toHexString() : '#000'; mapperListScheduleRedraw(); }
      });
      $('#boldCheckbox_' + layer).prop('checked', Styling_Option_dict[layer].Bold);
      $('#italicCheckbox_' + layer).prop('checked', Styling_Option_dict[layer].Italic);
      $('#underlineCheckbox_' + layer).prop('checked', Styling_Option_dict[layer].Underline);
      $('#fontSizeSlider_' + layer).val(parseInt(Styling_Option_dict[layer].Fontsize));
      $('#fontSizeSliderLabel_' + layer).text(parseInt(Styling_Option_dict[layer].Fontsize));

      $('#boldCheckbox_'+layer+', #italicCheckbox_'+layer+', #underlineCheckbox_'+layer).on('change', function () {
        Styling_Option_dict[layer].Bold      = $('#boldCheckbox_'+layer).prop('checked');
        Styling_Option_dict[layer].Italic    = $('#italicCheckbox_'+layer).prop('checked');
        Styling_Option_dict[layer].Underline = $('#underlineCheckbox_'+layer).prop('checked');
        mapperListScheduleRedraw();
      });
      $('#fontSizeSlider_' + layer).on('input', function () {
        var v = $(this).val();
        $('#fontSizeSliderLabel_' + layer).text(v);
        Styling_Option_dict[layer].Fontsize = v + 'px';
        mapperListScheduleRedraw();
      });
    });
  }

  // ── Row management ────────────────────────────────────────────────────────────
  function mapperListCreateReceptorTd($td) {
    var $hid  = $('<input type="hidden" class="mapper-core-receptor-entry" value="">');
    var $wrap = $('<div class="mapper-core-receptor-input-wrap is-empty">');
    var $inp  = $('<textarea class="form-control input-sm mapper-core-in-receptor" rows="1" autocomplete="off" spellcheck="false"></textarea>');
    var $view = $('<div class="form-control input-sm mapper-core-receptor-html-view" tabindex="0"></div>');
    var $clr  = $('<button type="button" class="mapper-core-receptor-clear" aria-label="Clear receptor">&times;</button>');
    $wrap.append($inp, $view, $clr);
    $td.append($hid, $wrap);
    $view.hide();
    $clr.on('click', function () { mapperListSetResolved($clr.closest('tr'), ''); });
    return $inp;
  }

  function mapperListRowBlank($tr) {
    var entry   = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
    var cat     = ($tr.find('.mapper-list-cat').val() || '').trim();
    var typed   = ($tr.find('.mapper-core-in-receptor').val() || '').trim();
    var hasVals = ['.mapper-list-val1','.mapper-list-val2','.mapper-list-val3','.mapper-list-val4'].some(function (s) {
      return (($tr.find(s).val() || '').trim() !== '');
    });
    return !entry && !cat && !typed && !hasVals && !$tr.data('mapperCoreUnmatchedRaw');
  }

  function mapperListSyncClearBtn($tr) {
    var has = !!($tr.find('.mapper-core-receptor-entry').val() || '').trim()
           || !!($tr.find('.mapper-core-in-receptor').val() || '').trim()
           || !!($tr.data('mapperCoreUnmatchedRaw'));
    $tr.find('.mapper-core-receptor-input-wrap').toggleClass('is-empty', !has);
  }

  function mapperListSyncRemoveButtons() {
    $('#mapper-list-input-tbody tr').each(function () {
      $(this).find('.mapper-core-remove-cell').toggleClass('is-remove-hidden', mapperListRowBlank($(this)));
    });
  }

  function mapperListCompactReceptors() {
    var hasContent = false;
    $('#mapper-list-input-tbody tr').each(function () { if (!mapperListRowBlank($(this))) { hasContent = true; return false; } });
    $('#mapper-list-input-table').toggleClass('mapper-core-receptors-compact', hasContent);
  }

  function mapperListAppendRow(skipTrail) {
    var $tr = $('<tr>');
    $tr.append($('<td class="mapper-core-remove-cell is-remove-hidden">').append(
      $('<button type="button" class="mapper-core-remove-row" aria-label="Remove row">&times;</button>')
    ));
    var $tdR = $('<td class="mapper-core-receptor-cell mapper-core-value-cell">');
    var $inp = mapperListCreateReceptorTd($tdR);
    $tr.append($tdR);
    // Species cell (hidden until species toggle is active)
    var $specSel = $('<select class="form-control input-sm mapper-list-species-select">');
    $tr.append($('<td class="mapper-list-species-cell">').css('display', mapperListSpeciesActive ? '' : 'none').append($specSel));

    // Val 1-4
    ['val1','val2','val3','val4'].forEach(function (suf) {
      var $v = $('<input type="text" class="form-control input-sm mapper-list-' + suf + ' mapper-list-val" autocomplete="off">');
      if (mapperListIsTextMode()) $v.prop('disabled', true);
      $tr.append($('<td class="mapper-list-val-cell mapper-core-value-cell">').append($v));
    });
    // Category
    var $cat = $('<input type="text" class="form-control input-sm mapper-list-cat" autocomplete="off">');
    $tr.append($('<td class="mapper-list-cat-cell mapper-core-value-cell">').append($cat));
    // Swatch
    $tr.append($('<td class="mapper-list-swatch-cell">').append(
      '<input type="text" class="mapper-list-label-swatch mapper-core-row-color-picker" readonly="readonly" aria-label="Label color">'
    ));

    $('#mapper-list-input-tbody').append($tr);
    mapperListBindAc($inp);
    mapperListSyncClearBtn($tr);
    mapperListSyncRemoveButtons();
    if (!skipTrail) mapperListEnsureTrailingBlankRow();
  }

  function mapperListEnsureTrailingBlankRow() {
    var $tb = $('#mapper-list-input-tbody');
    while ($tb.children().length >= 2) {
      var $last = $tb.children().last();
      var $prev = $last.prev();
      if (mapperListRowBlank($last) && mapperListRowBlank($prev)) {
        mapperListDestroyRowAc($last); $last.remove(); continue;
      }
      break;
    }
    if (!$tb.children().length) mapperListAppendRow(true);
    if (!mapperListRowBlank($tb.children().last()) && $tb.children().length < MAPPER_LIST_MAX_ROWS) {
      mapperListAppendRow(true);
    }
  }

  function mapperListDestroyAc($inp) {
    try { if ($inp.hasClass('ui-autocomplete-input')) $inp.autocomplete('destroy'); } catch (e) {}
  }
  function mapperListDestroyRowAc($tr) {
    mapperListDestroyAc($tr.find('.mapper-core-in-receptor'));
    $tr.find('.mapper-core-row-color-picker').each(function () { mapperListDestroySwatchSpectrum($(this)); });
    var $ss = $tr.find('.mapper-list-species-select');
    if ($ss.length && $ss.data('select2')) { try { $ss.select2('destroy'); } catch (e) {} }
  }
  function mapperListDestroyAllRows() {
    $('#mapper-list-input-tbody tr').each(function () { mapperListDestroyRowAc($(this)); });
  }

  // ── Autocomplete ──────────────────────────────────────────────────────────────
  function mapperListResolvedDisplay(entryId) {
    var meta = window.MAPPER_CORE_ENTRY_META && window.MAPPER_CORE_ENTRY_META[entryId] || {};
    if (LabelNames === 'Gene' && meta.gene) return $('<span/>').text(meta.gene).html();
    if (LabelNames === 'UniProt' && meta.uniprot) return $('<span/>').text(meta.uniprot).html();
    return meta.name_html ? String(meta.name_html) : $('<span/>').text(entryId || '').html();
  }

  function mapperListEditSeed(entryId) {
    var sid = entryId != null ? String(entryId).trim() : '';
    var meta = (window.MAPPER_CORE_ENTRY_META && window.MAPPER_CORE_ENTRY_META[sid]) || {};
    if (LabelNames === 'Gene' && meta.gene) { return meta.gene; }
    if (LabelNames === 'UniProt' && meta.uniprot) { return meta.uniprot; }
    return meta.name_plain || '';
  }

  function mapperListSetResolved($tr, id) {
    var sid = id != null ? String(id).trim() : '';
    var $inp  = $tr.find('.mapper-core-in-receptor');
    var $hid  = $tr.find('.mapper-core-receptor-entry');
    var $view = $tr.find('.mapper-core-receptor-html-view');
    mapperListDestroyAc($inp);
    if (!sid) {
      $hid.val(''); $view.hide().empty(); $inp.val('').show();
      mapperListBindAc($inp);
      $tr.removeClass('mapper-core-row-invalid').removeData('mapperCoreUnmatchedRaw');
      mapperListSyncClearBtn($tr); mapperListSyncRemoveButtons(); mapperListCompactReceptors();
      mapperListScheduleRedraw(); return;
    }
    $hid.val(sid);
    if (mapperListSpeciesActive) mapperListPopulateSpeciesSelect($tr, sid, false);
    $view.html(mapperListResolvedDisplay(sid)).show();
    $inp.val('').hide();
    mapperListBindAc($inp);
    $tr.removeClass('mapper-core-row-invalid').removeData('mapperCoreUnmatchedRaw');
    mapperListSyncClearBtn($tr); mapperListSyncRemoveButtons(); mapperListCompactReceptors();
    if (!suppressRedraw) {
      window.setTimeout(function () {
        var $cat = $tr.find('.mapper-list-cat:visible, .mapper-list-val1:visible').first();
        if ($cat.length) $cat.focus().select();
      }, 0);
    }
    mapperListScheduleRedraw();
  }

  function mapperListFilterLocal(term) {
    var t = (term || '').trim().toUpperCase();
    if (!t || !window.receptorSelect2Data) return [];
    var meta = window.MAPPER_CORE_ENTRY_META || {};
    var results = [];
    window.receptorSelect2Data.forEach(function (item) {
      var m = meta[item.id] || {};
      var searchIn = [(m.name_plain || item.name_plain || item.text || ''), (m.gene || ''), (m.uniprot || ''), item.id].join(' ').toUpperCase();
      if (searchIn.indexOf(t) === -1) return;
      var dispHtml;
      if (LabelNames === 'Gene' && m.gene) {
        dispHtml = $('<span/>').text(m.gene).html();
      } else if (LabelNames === 'UniProt' && m.uniprot) {
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

  function mapperListBindAc($inp) {
    mapperListDestroyAc($inp);
    $inp.autocomplete({
      minLength: 1,
      source: function (req, resp) { resp(mapperListFilterLocal(req.term)); },
      focus: function () { return false; },
      select: function (ev, ui) { mapperListSetResolved($inp.closest('tr'), ui.item.id); ev.preventDefault(); }
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
      mapperListSyncClearBtn($tr);
    });
    $inp.on('blur.mapperList', function () {
      var $tr = $inp.closest('tr');
      window.setTimeout(function () {
        if (!$inp.is(':visible')) return;
        var raw = ($inp.val() || '').trim();
        if (!raw || $tr.find('.mapper-core-receptor-entry').val()) return;
        var up = raw.toUpperCase();
        var rid = window.MAPPER_CORE_RESOLVE && window.MAPPER_CORE_RESOLVE[up];
        if (rid) { mapperListSetResolved($tr, rid); return; }
        var hit = (window.receptorSelect2Data||[]).filter(function(x){return String(x.id)===raw;});
        if (hit.length) { mapperListSetResolved($tr, hit[0].id); return; }
        if (raw) { $tr.addClass('mapper-core-row-invalid'); $tr.data('mapperCoreUnmatchedRaw', raw); }
        mapperListSyncClearBtn($tr);
      }, 170);
    });
  }

  // ── Sort ──────────────────────────────────────────────────────────────────────
  var SORT_STATE = { col: null, dir: 'asc' };

  function mapperListSerializeRows() {
    var rows = [];
    $('#mapper-list-input-tbody tr').each(function (ix) {
      var $tr = $(this);
      rows.push({
        ix: ix,
        entry: ($tr.find('.mapper-core-receptor-entry').val() || '').trim(),
        typed: ($tr.find('.mapper-core-in-receptor').val() || '').trim(),
        unmatched: ($tr.data('mapperCoreUnmatchedRaw') || ''),
        invalid: $tr.hasClass('mapper-core-row-invalid'),
        v1: ($tr.find('.mapper-list-val1').val() || '').trim(),
        v2: ($tr.find('.mapper-list-val2').val() || '').trim(),
        v3: ($tr.find('.mapper-list-val3').val() || '').trim(),
        v4: ($tr.find('.mapper-list-val4').val() || '').trim(),
        cat: ($tr.find('.mapper-list-cat').val() || '').trim(),
        species: ($tr.find('.mapper-list-species-select').val() || '').trim()
      });
    });
    return rows;
  }

  function mapperListApplySerializedRows(rows) {
    suppressRedraw = true;
    mapperListDestroyAllRows();
    $('#mapper-list-input-tbody').empty();
    rows.forEach(function (r) {
      mapperListAppendRow(true);
      var $tr = $('#mapper-list-input-tbody tr').last();
      if (r.entry) {
        mapperListSetResolved($tr, r.entry);
        if (r.species && mapperListSpeciesActive) {
          var $specSel = $tr.find('.mapper-list-species-select');
          $specSel.val(r.species);
          if ($specSel.data('select2')) { $specSel.trigger('change.select2'); }
        }
      } else if (r.typed) {
        $tr.find('.mapper-core-in-receptor').val(r.typed);
        if (r.unmatched) $tr.data('mapperCoreUnmatchedRaw', r.unmatched);
        if (r.invalid)   $tr.addClass('mapper-core-row-invalid');
      }
      $tr.find('.mapper-list-val1').val(r.v1);
      $tr.find('.mapper-list-val2').val(r.v2);
      $tr.find('.mapper-list-val3').val(r.v3);
      $tr.find('.mapper-list-val4').val(r.v4);
      $tr.find('.mapper-list-cat').val(r.cat);
    });
    suppressRedraw = false;
    mapperListEnsureTrailingBlankRow();
    mapperListSyncRemoveButtons();
    mapperListCompactReceptors();
    mapperListRefreshSwatches();
    mapperListSyncCatHints();
    mapperListSyncFirstRowPlaceholder();
  }

  function mapperListSortRows(col) {
    if (SORT_STATE.col === col) SORT_STATE.dir = SORT_STATE.dir === 'asc' ? 'desc' : 'asc';
    else { SORT_STATE.col = col; SORT_STATE.dir = 'asc'; }

    var rows = mapperListSerializeRows();
    rows.sort(function (a, b) {
      var va, vb;
      if (col === 'receptor') { va = a.entry || a.typed || ''; vb = b.entry || b.typed || ''; }
      else if (col === 'cat') { va = a.cat; vb = b.cat; }
      else {
        va = parseFloat(a['v' + col.replace('val','')]) || 0;
        vb = parseFloat(b['v' + col.replace('val','')]) || 0;
        return SORT_STATE.dir === 'asc' ? va - vb : vb - va;
      }
      var c = va.localeCompare(vb, undefined, { numeric: true, sensitivity: 'base' });
      return SORT_STATE.dir === 'asc' ? c : -c;
    });
    mapperListApplySerializedRows(rows);
    mapperListUpdateSortHeaders();
    mapperListScheduleRedraw();
  }

  function mapperListUpdateSortHeaders() {
    $('#mapper-list-input-table th.mapper-core-sortable-head').each(function () {
      var col    = $(this).attr('data-mapper-list-sort-col');
      var active = SORT_STATE.col === col;
      var dir    = active ? SORT_STATE.dir : null;
      $(this).attr('aria-sort', active ? (dir==='desc'?'descending':'ascending') : 'none');
      $(this).find('.mapper-core-sort-indicator').text(active ? (dir==='desc'?'▼':'▲') : '↕');
    });
  }

  // ── Demo ──────────────────────────────────────────────────────────────────────
  function mapperListApplyDemoColorPresets() {
    if (!$.fn.spectrum) return;
    var DEMO = [
      { col: 'Col1', style: 'Three_RWB' },
      { col: 'Col2', style: 'Two'       },
      { col: 'Col3', style: 'Three_TWM' },
      { col: 'Col4', style: 'One'       }
    ];
    // Custom colours that override each preset's defaults (mirrors tree approach)
    var CUSTOM = {
      Col1: ['#a00000', '#1a80bb'],
      Col2: ['#97a6c4', '#1a2b3c'],
      Col3: ['#298c8c', '#800074'],
      Col4: ['#ffffff', '#2ca02c']
    };

    // Step 1: trigger change on each select2 — updates display widget AND calls
    //         applyColStylePreset (sets picker visibility + ColorConfig from preset)
    DEMO.forEach(function (d) {
      try { $('#mapper-list-cstyle-' + d.col).val(d.style).trigger('change'); } catch (e) {}
    });

    // Step 2: override ColorConfig with the custom demo palette
    Object.keys(CUSTOM).forEach(function (col) {
      ColorConfig[col].color1 = CUSTOM[col][0];
      ColorConfig[col].color2 = CUSTOM[col][1];
    });

    // Step 3: sync spectrum pickers to the custom colours
    DEMO.forEach(function (d) {
      try { $('#mapper-list-cpicker-min-' + d.col).spectrum('set', CUSTOM[d.col][0]); } catch (e) {}
      try { $('#mapper-list-cpicker-max-' + d.col).spectrum('set', CUSTOM[d.col][1]); } catch (e) {}
    });
  }

  function mapperListFillDemo() {
    suppressRedraw = true;
    mapperListDestroyAllRows();
    $('#mapper-list-input-tbody').empty();
    var text = mapperListIsTextMode();

    // Pre-assign categorical colours so each label gets a distinct colour
    if (text) {
      window.MAPPER_CORE_LABEL_COLORS = {};
      window.MAPPER_CORE_LABEL_ENABLED = {};
      Object.keys(DEMO_CAT_COLORS).forEach(function (cat) {
        window.MAPPER_CORE_LABEL_COLORS[cat] = DEMO_CAT_COLORS[cat];
        window.MAPPER_CORE_LABEL_ENABLED[cat] = true;
      });
    }

    MAPPER_LIST_DEMO_ROWS.forEach(function (r) {
      mapperListAppendRow(true);
      var $tr = $('#mapper-list-input-tbody tr').last();
      var up = r.receptor.toUpperCase();
      var rid = window.MAPPER_CORE_RESOLVE && window.MAPPER_CORE_RESOLVE[up];
      if (rid) {
        mapperListSetResolved($tr, rid);
      } else {
        $tr.find('.mapper-core-in-receptor').val(r.receptor);
      }
      // Fill only the current mode's fields — Numbers and Categories stay fully separate
      if (text) {
        $tr.find('.mapper-list-cat').val(r.text);
      } else {
        $tr.find('.mapper-list-val1').val(r.vals[0]);
        $tr.find('.mapper-list-val2').val(r.vals[1]);
        $tr.find('.mapper-list-val3').val(r.vals[2]);
        $tr.find('.mapper-list-val4').val(r.vals[3]);
      }
    });
    suppressRedraw = false;
    $('#mapper-list-clear-rows').removeClass('mapper-core-clear-clean');
    mapperListEnsureTrailingBlankRow();
    mapperListSyncRemoveButtons();
    mapperListCompactReceptors();
    mapperListRefreshSwatches();
    mapperListSyncCatHints();
    mapperListSyncFirstRowPlaceholder();
    if (!text) mapperListApplyDemoColorPresets();
    mapperListRedrawNow();
  }

  // ── Paste ─────────────────────────────────────────────────────────────────────
  function mapperListFindBlankRow() {
    var $found = $();
    $('#mapper-list-input-tbody tr').each(function () { if (mapperListRowBlank($(this))) { $found = $(this); return false; } });
    return $found;
  }

  // ── Boot ──────────────────────────────────────────────────────────────────────
  function mapperListBoot() {
    $('.mapper-core-booting').removeClass('mapper-core-booting');
    // Load receptor info from Django context
    ReceptorInfo = window.MAPPER_LIST_RECEPTOR_INFO || {};
    mapperListBuildLabelConversion();
    try { mapperListExtendLabelConversionForSpecies(); } catch (e) { /* non-critical */ }

    // Extend autocomplete with non-human-only entries (receptors with no human ortholog)
    try {
      if (window.MAPPER_LIST_SPECIES_DATA) {
        var _seenNhoStems = {};
        (window.MAPPER_LIST_SPECIES_DATA.nonhuman_only || []).forEach(function (nho) {
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

    // Set up initial table
    $('#mapper-list-input-tbody').empty();
    mapperListAppendRow();
    mapperListSetInputMode('numeric'); // sets table state
    mapperListSyncFirstRowPlaceholder(); // explicit call since no-op on first boot
    mapperListUpdateSortHeaders();

    // Init colour pickers
    mapperListInitColorPickers();
    mapperListInitLabelStylePickers();

    // Legend label inputs
    [['Val1','Col1'],['Val2','Col2'],['Val3','Col3'],['Val4','Col4']].forEach(function (m) {
      $('#mapper-list-lbl-' + m[0]).val(Barlabels[m[1]]).on('input', function () {
        Barlabels[m[1]] = $(this).val() || Barlabels[m[1]];
        mapperListScheduleRedraw();
      });
    });

    // Receptor names buttons — update display in table AND schedule plot redraw (mirrors tree)
    $('.mapper-list-leaf-btn').on('click.mapperListLeaf', function () {
      var v = ($(this).attr('data-value') || 'Protein').trim();
      LabelNames = v;
      $('.mapper-list-leaf-btn').each(function () {
        var ok = ($(this).attr('data-value') || '') === v;
        $(this).toggleClass('btn-primary', ok).toggleClass('btn-outline-primary', !ok);
      });
      // Refresh resolved receptor names shown in the input table html-view
      $('#mapper-list-input-tbody tr').each(function () {
        var sid = ($( this).find('.mapper-core-receptor-entry').val() || '').trim();
        if (!sid) return;
        var $view = $(this).find('.mapper-core-receptor-html-view');
        if ($view.length && $view.is(':visible')) $view.html(mapperListResolvedDisplay(sid));
      });
      mapperListScheduleRedraw();
    });

    // Layout controls
    $('#mapper-list-columns').on('change', function () {
      Layout_dict.columns = parseInt($(this).val()) || 2;
      mapperListScheduleRedraw();
    });
    $('#mapper-list-max-label-sel').on('change', function () {
      Layout_dict.col_max_label = $(this).val();
      $('#mapper-list-max-label-num').prop('disabled', $(this).val() !== 'Custom');
      mapperListScheduleRedraw();
    });
    $('#mapper-list-max-label-num').on('input', function () {
      Layout_dict.Col_break_number = parseInt($(this).val()) || 1;
      mapperListScheduleRedraw();
    });
    // Layer checkboxes
    d3.selectAll('#toggle-layer-1, #toggle-layer-2, #toggle-layer-3').on('change', function () { mapperListScheduleRedraw(); });

    // Mode switch
    $('#mapper-list-mode-numeric-btn').on('click.mapperList', function () { mapperListSetInputMode('numeric'); });
    $('#mapper-list-mode-labels-btn').on('click.mapperList', function () { mapperListSetInputMode('text'); });

    // Demo
    $('#mapper-list-demo-rows').on('click.mapperList', function () { mapperListFillDemo(); });

    // Clear
    function mapperListDoWholeTableClear() {
      mapperListDestroyAllRows();
      $('#mapper-list-input-tbody').empty();
      MAPPER_LIST_MODE_SNAPSHOTS.numeric = null;
      MAPPER_LIST_MODE_SNAPSHOTS.categorical = null;
      mapperListAppendRow();
      mapperListCompactReceptors();
      mapperListSyncFirstRowPlaceholder();
      $('#mapper-list-clear-rows').addClass('mapper-core-clear-clean').blur();
      mapperListRedrawNow();
    }
    function mapperListDoCurrentModeClear() {
      var curKey = mapperListIsTextMode() ? 'categorical' : 'numeric';
      mapperListDestroyAllRows();
      $('#mapper-list-input-tbody').empty();
      MAPPER_LIST_MODE_SNAPSHOTS[curKey] = null;
      mapperListAppendRow();
      mapperListCompactReceptors();
      mapperListSyncFirstRowPlaceholder();
      $('#mapper-list-clear-rows').addClass('mapper-core-clear-clean').blur();
      mapperListRedrawNow();
    }
    $('#mapper-list-clear-whole').on('click.mapperList', function(e) {
      e.preventDefault(); mapperListDoWholeTableClear();
    });
    $('#mapper-list-clear-mode').on('click.mapperList', function(e) {
      e.preventDefault(); mapperListDoCurrentModeClear();
    });
    $('.mapper-core-clear-menu').on('click.mapperList', '[data-list-col]', function(e) {
      e.preventDefault();
      var col = $(this).data('list-col');
      $('#mapper-list-input-tbody tr').each(function() {
        $(this).find('.mapper-list-' + col).val('').trigger('input');
      });
      $('#mapper-list-clear-rows').removeClass('mapper-core-clear-clean');
      mapperListScheduleRedraw();
    });

    // Remove row
    $('#mapper-list-input-table').on('click.mapperListRemove', '.mapper-core-remove-row', function () {
      var $tr = $(this).closest('tr');
      if (mapperListRowBlank($tr) && $tr.is(':last-child')) return;
      mapperListDestroyRowAc($tr); $tr.remove();
      $('#mapper-list-clear-rows').removeClass('mapper-core-clear-clean');
      mapperListEnsureTrailingBlankRow();
      _listCosmeticSyncDebounced.schedule();
      mapperListScheduleRedraw();
    });

    // Click resolved chip to re-enter edit mode
    $(document)
      .off('click.mapperListHtmlEdit', '#mapper-list-input-tbody .mapper-core-receptor-html-view')
      .on('click.mapperListHtmlEdit', '#mapper-list-input-tbody .mapper-core-receptor-html-view', function () {
        var $tr  = $(this).closest('tr');
        mapperListDestroyAc($tr.find('.mapper-core-in-receptor'));
        var hid  = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
        $(this).hide().empty();
        var $inp2 = $('<textarea class="form-control input-sm mapper-core-in-receptor" rows="1" autocomplete="off" spellcheck="false"></textarea>');
        var seed  = mapperListEditSeed(hid);
        $inp2.val(seed);
        if (!seed) { $tr.find('.mapper-core-receptor-input-wrap').addClass('is-empty'); }
        $tr.find('.mapper-core-receptor-input-wrap .mapper-core-in-receptor').remove();
        $tr.find('.mapper-core-receptor-input-wrap').prepend($inp2);
        $tr.find('.mapper-core-receptor-entry').val('');
        mapperListBindAc($inp2);
        $inp2.show().focus();
        window.setTimeout(function () {
          var t = ($inp2.val() || '').trim();
          if (t.length >= 1 && $inp2.data('ui-autocomplete')) { $inp2.autocomplete('search', t); }
        }, 0);
        mapperListSyncClearBtn($tr);
        mapperListSyncRemoveButtons();
        mapperListCompactReceptors();
        mapperListEnsureTrailingBlankRow();
        mapperListScheduleRedraw();
      });

    // Sort headers
    $('#mapper-list-input-table').on('click.mapperListSort', 'th.mapper-core-sortable-head', function () {
      mapperListSortRows($(this).attr('data-mapper-list-sort-col'));
    });

    // Input changes
    $(document).on('input.mapperList blur.mapperList change.mapperList',
      '#mapper-list-input-tbody input, #mapper-list-input-tbody textarea',
      function () {
        var $tr = $(this).closest('tr');
        mapperListSyncClearBtn($tr);   // keep is-empty in sync → removes orange border when filled
        $('#mapper-list-clear-rows').removeClass('mapper-core-clear-clean');
        mapperListEnsureTrailingBlankRow();
        _listCosmeticSyncDebounced.schedule();
        mapperListScheduleRedraw();
      }
    );

    // Paste
    $('#mapper-list-input-table').on('paste.mapperListPaste', function (ePz) {
      var ev   = ePz.originalEvent || ePz;
      var text = ev.clipboardData ? ev.clipboardData.getData('text/plain') : '';
      if (!text || text.indexOf('\t') === -1) return;
      ePz.preventDefault();
      $('#mapper-list-clear-rows').removeClass('mapper-core-clear-clean');
      suppressRedraw = true;
      text.split(/\r?\n/).forEach(function (ln) {
        if (!ln.trim()) return;
        var parts = ln.split('\t');
        var rawR  = (parts[0] || '').trim();
        var $tr   = mapperListFindBlankRow();
        if (!$tr.length) { mapperListAppendRow(true); $tr = $('#mapper-list-input-tbody tr').last(); }
        $tr.find('.mapper-core-in-receptor').val(rawR);
        var up = rawR.toUpperCase();
        var rid = window.MAPPER_CORE_RESOLVE && window.MAPPER_CORE_RESOLVE[up];
        if (rid) { mapperListSetResolved($tr, rid); }
        if (!mapperListIsTextMode()) {
          if (parts[1]) $tr.find('.mapper-list-val1').val((parts[1]||'').trim());
          if (parts[2]) $tr.find('.mapper-list-val2').val((parts[2]||'').trim());
          if (parts[3]) $tr.find('.mapper-list-val3').val((parts[3]||'').trim());
          if (parts[4]) $tr.find('.mapper-list-val4').val((parts[4]||'').trim());
        } else {
          if (parts[1]) $tr.find('.mapper-list-cat').val((parts[1]||'').trim());
        }
      });
      suppressRedraw = false;
      mapperListEnsureTrailingBlankRow();
      mapperListSyncRemoveButtons(); mapperListCompactReceptors();
      mapperListRefreshSwatches(); mapperListSyncCatHints();
      mapperListSyncFirstRowPlaceholder(); mapperListScheduleRedraw();
    });

    // Legend (categorical)
    $('#mapper-list-legend-toggle').on('click.mapperList', function () {
      ShowLegend = !ShowLegend;
      $(this).text(ShowLegend ? 'Shown' : 'Hidden')
             .toggleClass('btn-success', ShowLegend)
             .toggleClass('btn-danger', !ShowLegend);
      mapperListRedrawNow();
    });
    $('#mapper-list-legend-layout').on('change', function () {
      var v = $(this).val();
      if (v.startsWith('columns')) { Textlegend_styling.layoutMode = 'columns'; Textlegend_styling.columns = parseInt(v.split('-')[1]); }
      else Textlegend_styling.layoutMode = 'row';
      mapperListRedrawNow();
    });
    $('#mapper-list-legend-sorting').on('click.mapperList', function () {
      Textlegend_styling.sortDirection = Textlegend_styling.sortDirection === 'Vertically' ? 'Horizontally' : 'Vertically';
      $(this).text(Textlegend_styling.sortDirection); mapperListRedrawNow();
    });
    $('#mapper-list-legend-position').on('click.mapperList', function () {
      Textlegend_styling.TreeLegendPosition = Textlegend_styling.TreeLegendPosition === 'Top' ? 'Bottom' : 'Top';
      Layout_dict.legend_position  = Textlegend_styling.TreeLegendPosition;
      Layout_dict.legend_top_space = Textlegend_styling.TreeLegendPosition === 'Bottom' ? 5 : 50;
      $(this).text(Textlegend_styling.TreeLegendPosition); mapperListRedrawNow();
    });
    $('#mapper-list-legend-fontsize').on('input', function () {
      var v = $(this).val();
      $('#mapper-list-legend-fontsize-val').text(v);
      Textlegend_styling.Fontsize = v + 'px';
      mapperListRedrawNow();
    });

    // Receptor lookup modal
    if (typeof window.mapperCoreInitGpcromePickerModal === 'function') {
      window.mapperCoreInitGpcromePickerModal({
        pickerRows: $.isArray(window.MAPPER_CORE_GPCROME_PICKER_ROWS) ? window.MAPPER_CORE_GPCROME_PICKER_ROWS : [],
        maxRows: MAPPER_LIST_MAX_ROWS,
        onAdd: function (entryIds, meta) {
          suppressRedraw = true;
          var isNumberMode = !mapperListIsTextMode();
          var numberMap = (meta && meta.groupNameById && isNumberMode)
            ? window.mapperCoreBuildSequentialNumberMap(meta.groupNameById)
            : null;
          (entryIds || []).forEach(function (id) {
            var $row = mapperListFindBlankRow();
            if (!$row.length) { mapperListAppendRow(true); $row = $('#mapper-list-input-tbody tr').last(); }
            mapperListSetResolved($row, id);
            if (meta && meta.groupNameById && meta.groupNameById[id] != null) {
              var $target = isNumberMode ? $row.find('.mapper-list-val1') : $row.find('.mapper-list-cat');
              var assignedValue = isNumberMode ? numberMap[id] : meta.groupNameById[id];
              if (assignedValue != null) {
                $target.val(String(assignedValue)).trigger('input');
              }
            }
          });
          suppressRedraw = false;
          $('#mapper-list-clear-rows').removeClass('mapper-core-clear-clean');
          mapperListEnsureTrailingBlankRow();
          mapperListCompactReceptors();
          mapperListScheduleRedraw();
        }
      });
    }

    // Species toggle
    $('#mapper-list-species-toggle').on('click.mapperList', function () {
      mapperListSpeciesActive = !mapperListSpeciesActive;
      $(this).toggleClass('btn-primary', mapperListSpeciesActive)
             .toggleClass('btn-default', !mapperListSpeciesActive);
      mapperListApplySpeciesPanelClass();
      $('#mapper-list-input-tbody .mapper-list-species-cell').toggle(mapperListSpeciesActive);
      $('.mapper-list-species-h').toggle(mapperListSpeciesActive);
      $('#mapper-list-species-names-wrap').toggle(mapperListSpeciesActive);
      mapperListExtendLabelConversionForSpecies();
      if (mapperListSpeciesActive) mapperListInitAllSpeciesDropdowns();
      mapperListScheduleRedraw();
    });

    // Species dropdown change
    $(document).on('change.mapperListSpecies',
      '#mapper-list-input-tbody .mapper-list-species-select',
      function () { mapperListScheduleRedraw(); }
    );

    // Species name format buttons — use direct binding (document delegation is blocked by
    // the stopPropagation on .dropdown-menu at line 835 of Mapper_List.html)
    $('.mapper-list-species-name-btn').off('click.mapperListSpeciesName').on('click.mapperListSpeciesName', function (eSpec) {
      eSpec.preventDefault();
      mapperListSpeciesNameFormat = $(this).data('value');
      $('.mapper-list-species-name-btn').each(function () {
        var ok = $(this).data('value') === mapperListSpeciesNameFormat;
        $(this).toggleClass('btn-primary', ok).toggleClass('btn-outline-primary', !ok);
      });
      mapperListExtendLabelConversionForSpecies();
      mapperListInitAllSpeciesDropdowns();
      mapperListScheduleRedraw();
    });

    mapperListShowPlaceholder();
    mapperListSyncClearDropdown();
  }

  // ── Entry point ───────────────────────────────────────────────────────────────
  $(document).ready(function () { mapperListBoot(); });

})(jQuery);

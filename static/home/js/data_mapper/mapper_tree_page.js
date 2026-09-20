/**
 * Mapper 2.0 classification tree — client-side keep_by_names filter, circles payload,
 * debounced redraw (mapper_classification_tree + DrawCircles).
 */
(function ($) {
  'use strict';

  var DEBOUNCE_MS = 140;
  var suppressRedraw = false;

  var ORIG_SKEL = null;
  var ORIG_OPTS = null;
  var skipLigandtypeLevel = false;
  var skipClassLevel = false;
  var classLevelAutoOff = false; // true when class was hidden automatically (single class)

  var ServerReceptorDict = {};
  var ServerGeneDict = {};

  window.mapperTreeInputMode = 'numeric';
  window.mapperTreeRedrawNow = function () {};

  window.MAPPER_CORE_LABEL_COLORS = window.MAPPER_CORE_LABEL_COLORS || {};
  window.MAPPER_CORE_LABEL_ENABLED = window.MAPPER_CORE_LABEL_ENABLED || {};
  /** Stem-upper → label arrays for SVG leaf relabelling (`custom_changeLeavesLabels`). */
  window.mapperTreeStemLabelDicts = { IUPHAR: {}, Gene: {}, UniProt: {} };

  var MAPPER_TREE_MAX_ROWS = 150;
  var MAPPER_TREE_MAX_RECEPTORS = 150;
  var MAPPER_TREE_PLACEHOLDER_NUMERIC  = 'Paste in 1–6 columns of data\nor type';
  var MAPPER_TREE_PLACEHOLDER_CATMODE  = 'Paste in 1–2 columns of data\nor type';
  /** Mapper tree page uses circular layout only (no curved/straight dendrogram UI). */
  var MAPPER_TREE_LAYOUT = 'Tree - Circular';
  var MAPPER_TREE_RADIAL_LEAF_LABEL_GAP_BASE = 10;
  var MAPPER_TREE_DEMO_ROWS = [
    { receptor: '5HT1A', numeric: [1, 73, -10.1, -25, 1], text: 'Serotonergic' },
    { receptor: '5HT1B', numeric: [2, 72, -10, -24, 2], text: 'Serotonergic' },
    { receptor: '5HT1D', numeric: [3, 71, -9.9, -23, 3], text: 'Serotonergic' },
    { receptor: '5HT1E', numeric: [4, 70, -9.8, -22, 4], text: 'Serotonergic' },
    { receptor: '5HT1F', numeric: [5, 69, -9.7, -21, 5], text: 'Serotonergic' },
    { receptor: '5HT2A', numeric: [6, 68, -9.6, -20, 6], text: 'Serotonergic' },
    { receptor: '5HT2B', numeric: [7, 67, -9.5, -19, 7], text: 'Serotonergic' },
    { receptor: '5HT2C', numeric: [8, 66, -9.4, -18, 8], text: 'Serotonergic' },
    { receptor: 'ACKR1', numeric: [13, 61, -8.9, -13, 13], text: 'Chemokine' },
    { receptor: 'ACKR2', numeric: [14, 60, -8.8, -12, 14], text: 'Chemokine' },
    { receptor: 'ACKR3', numeric: [15, 59, -8.7, -11, 15], text: 'Chemokine' },
    { receptor: 'ACKR4', numeric: [16, 58, -8.6, -10, 16], text: 'Chemokine' },
    { receptor: 'ACM1', numeric: [17, 57, -8.5, -9, 17], text: 'Cholinergic' },
    { receptor: 'ACM2', numeric: [18, 56, -8.4, -8, 18], text: 'Cholinergic' },
    { receptor: 'ACM3', numeric: [19, 55, -8.3, -7, 19], text: 'Cholinergic' },
    { receptor: 'ADA1A', numeric: [23, 51, -7.9, -3, 23], text: 'Adrenergic' },
    { receptor: 'ADA1B', numeric: [24, 50, -7.8, -2, 24], text: 'Adrenergic' },
    { receptor: 'ADA1D', numeric: [25, 49, -7.7, -1, 25], text: 'Adrenergic' },
    { receptor: 'ADRB1', numeric: [29, 45, -7.3, 3, 29], text: 'Adrenergic' },
    { receptor: 'ADRB2', numeric: [30, 44, -7.2, 4, 30], text: 'Adrenergic' },
    { receptor: 'ADGRA1', numeric: [31, 43, -7.1, 5, 31], text: 'Adhesion' },
    { receptor: 'ADGRA2', numeric: [32, 42, -7.0, 6, 32], text: 'Adhesion' },
    { receptor: 'ADGRA3', numeric: [33, 41, -6.9, 7, 33], text: 'Adhesion' }
  ];

  var Tree_circles;
  var Tree_colors;
  var Tree_circle_styling_dict;
  var Label_dict;
  var Tree_datatypes_dict;
  var Tree_textlegend_styling;
  var ShowLegend;
  var TreeLegendPosition;
  var styling_circles;
  var maxLeafNodeLength_scaler;
  var mapperTreeLastLabelSet = '';

  var mapperTreeSpeciesActive = false;
  var mapperTreeSpeciesNameFormat = 'common'; // 'common' | 'latin' | 'both'

  function deepClone(obj) {
    return JSON.parse(JSON.stringify(obj));
  }

  function mapperJsKeepByNames(node, namesToKeep) {
    if (Array.isArray(node)) {
      var kept = [];
      for (var i = 0; i < node.length; i++) {
        var xi = mapperJsKeepByNames(node[i], namesToKeep);
        if (xi != null) {
          kept.push(xi);
        }
      }
      return kept.length ? kept : null;
    }
    if (!node || typeof node !== 'object') {
      return node;
    }
    var nm = node.name;
    var hasKids = !!(node.children && node.children.length);
    if (!Object.prototype.hasOwnProperty.call(namesToKeep, nm)) {
      if (hasKids) {
        var ch = mapperJsKeepByNames(node.children, namesToKeep);
        if (!ch || !ch.length) {
          return null;
        }
        var o = deepClone(node);
        o.children = ch;
        return o;
      }
      return null;
    }
    var pay = namesToKeep[nm];
    var out = deepClone(node);
    if (pay && Object.prototype.hasOwnProperty.call(pay, 'Inner')) {
      out.value = pay.Inner;
    }
    if (out.children && out.children.length) {
      var k2 = mapperJsKeepByNames(out.children, namesToKeep);
      if (k2 && k2.length) {
        out.children = k2;
      } else {
        delete out.children;
      }
    }
    return out;
  }

  function mapperTreeMaybePromoteRoot(md, opts) {
    var o = deepClone(opts);
    if (md && md.children && md.children.length === 1) {
      return {
        tree: deepClone(md.children[0]),
        opts: $.extend(o, {
          depth: 3,
          branch_length: { 1: 'Alicarboxylic acid', 2: 'Gonadotrophin-releasing hormone', 3: '' }
        })
      };
    }
    return { tree: md, opts: o };
  }

  function mapperTreeStripLigandtypeLevel(skel) {
    if (!skel || !skel.children) { return skel; }
    var stripped = deepClone(skel);
    stripped.children = (skel.children || []).map(function (classNode) {
      var newClass = deepClone(classNode);
      var hoisted = [];
      (classNode.children || []).forEach(function (ltNode) {
        (ltNode.children || []).forEach(function (familyNode) {
          hoisted.push(deepClone(familyNode));
        });
      });
      newClass.children = hoisted;
      return newClass;
    });
    return stripped;
  }

  function mapperTreeStripClassLevel(skel) {
    // Root → Class → Chemotype → RF → Leaf  becomes  Root → Chemotype → RF → Leaf
    if (!skel || !skel.children) { return skel; }
    var stripped = deepClone(skel);
    var hoisted = [];
    (skel.children || []).forEach(function (classNode) {
      (classNode.children || []).forEach(function (chemotypeNode) {
        var n = deepClone(chemotypeNode);
        n._classNodeName = classNode.name; // preserve for class-color lookup
        hoisted.push(n);
      });
    });
    stripped.children = hoisted;
    return stripped;
  }

  function mapperTreeStripBothLevels(skel) {
    // Root → Class → Chemotype → RF → Leaf  becomes  Root → RF → Leaf
    if (!skel || !skel.children) { return skel; }
    var stripped = deepClone(skel);
    var hoisted = [];
    (skel.children || []).forEach(function (classNode) {
      (classNode.children || []).forEach(function (chemotypeNode) {
        (chemotypeNode.children || []).forEach(function (rfNode) {
          var n = deepClone(rfNode);
          n._classNodeName = classNode.name; // preserve for class-color lookup
          hoisted.push(n);
        });
      });
    });
    stripped.children = hoisted;
    return stripped;
  }

  function mapperStemFromEntry(entryId) {
    return String(entryId || '').trim().replace(/_human$/i, '');
  }

  function mapperTreeStemKeyUpper(entryId) {
    return mapperStemFromEntry(entryId).toUpperCase();
  }

  function mapperFallbackUniprotFromEntry(entryId, meta) {
    var sid = entryId != null ? String(entryId).trim().toUpperCase() : '';
    var u = meta && meta.uniprot ? String(meta.uniprot).trim().toUpperCase() : '';
    if (u) {
      return u;
    }
    if (sid.indexOf('_') !== -1) {
      return sid.split('_')[0].trim().toUpperCase();
    }
    return sid;
  }

  function mapperTreeBuildStemLabelDicts() {
    window.mapperTreeStemLabelDicts = { IUPHAR: {}, Gene: {}, UniProt: {} };
    if (!window.receptorSelect2Data || !window.receptorSelect2Data.length) {
      return;
    }
    window.receptorSelect2Data.forEach(function (item) {
      var id = item.id;
      if (id == null) {
        return;
      }
      var meta = (window.MAPPER_CORE_ENTRY_META && window.MAPPER_CORE_ENTRY_META[id]) || {};
      var stemU = mapperTreeStemKeyUpper(id);
      if (!stemU) {
        return;
      }
      var uni = mapperFallbackUniprotFromEntry(id, meta);

      window.mapperTreeStemLabelDicts.UniProt[stemU] = [uni];

      var iuphar = '';
      if (meta.name_html) {
        iuphar = meta.name_html;
      } else if (meta.name_plain) {
        iuphar = meta.name_plain;
      } else if (item.name_html) {
        iuphar = String(item.name_html);
      } else if (item.text) {
        iuphar = String(item.text);
      }
      if (iuphar) {
        window.mapperTreeStemLabelDicts.IUPHAR[stemU] = [iuphar];
      }

      var g = meta.gene || item.gene;
      if (g) {
        window.mapperTreeStemLabelDicts.Gene[stemU] = [String(g).trim()];
      }

      var rd = ServerReceptorDict[uni];
      if (!meta.name_html && rd && rd.length) {
        window.mapperTreeStemLabelDicts.IUPHAR[stemU] = rd;
      }
      var eg = ServerGeneDict[uni];
      if (eg && eg.length) {
        window.mapperTreeStemLabelDicts.Gene[stemU] = eg;
      }
    });
  }

  function mapperTreeLeafLabelLookupBuild() {
    window.tree_leaf_label_lookup = {};
    if (!window.receptorSelect2Data || !window.receptorSelect2Data.length) {
      return;
    }
    window.receptorSelect2Data.forEach(function (item) {
      var id = item.id;
      if (id == null) {
        return;
      }
      var meta = (window.MAPPER_CORE_ENTRY_META && window.MAPPER_CORE_ENTRY_META[id]) || {};
      var uni = mapperFallbackUniprotFromEntry(id, meta);
      var stemKey = mapperTreeStemKeyUpper(id);
      function row() {
        var proteinHtml = meta.name_html || item.name_html || '';
        var proteinPlain = (meta.name_plain || item.name_plain || item.text || id).trim();
        return {
          Protein: proteinPlain,
          ProteinHtml: proteinHtml || $('<span/>').text(proteinPlain).html(),
          Gene: (meta.gene || '').trim(),
          UniProt: uni || stemKey || ''
        };
      }
      var r = row();
      if (stemKey) {
        window.tree_leaf_label_lookup[stemKey] = r;
      }
      if (uni) {
        window.tree_leaf_label_lookup[uni] = r;
      }
    });
  }

  function mapperTreeBareStemUpper(entryId) {
    return String(entryId || '').trim().split('_')[0].toUpperCase();
  }

  // Returns the dedup key for a row (matches mapperTreeBuildTreeCircles key logic).
  // Returns null for blank / unresolvable rows.
  function mapperTreeRowUniKey($tr) {
    var entry = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
    var ta    = (($tr.find('.mapper-core-in-receptor').val() || '') + '').trim();
    var unmatched = ($tr.data('mapperCoreUnmatchedRaw') || '') + '';
    if (!entry && !ta && !String(unmatched).trim()) { return null; }
    var id = entry || (window.MAPPER_CORE_RESOLVE &&
             window.MAPPER_CORE_RESOLVE[String(ta || unmatched).trim().toUpperCase()]) || '';
    if (!id) { return null; }
    if (mapperTreeSpeciesActive) {
      var specVal = ($tr.find('.mapper-tree-species-select').val() || '').trim();
      return specVal ? specVal.toUpperCase() : mapperTreeBareStemUpper(id);
    }
    return mapperTreeBareStemUpper(id);
  }

  // Mark rows that are duplicates or push the unique-leaf count past MAPPER_TREE_MAX_RECEPTORS.
  // Runs synchronously (fast DOM scan) for immediate visual feedback.
  function mapperTreeUpdateOverLimitMarks() {
    var seen = {};
    var unique = 0;
    $('#mapper-tree-input-tbody tr').each(function () {
      var $tr = $(this);
      var key = mapperTreeRowUniKey($tr);
      var over = false;
      if (key !== null) {
        if (Object.prototype.hasOwnProperty.call(seen, key) || unique >= MAPPER_TREE_MAX_RECEPTORS) {
          over = true;
        } else {
          seen[key] = true;
          unique++;
        }
      }
      $tr.toggleClass('mapper-tree-row-over-limit', over);
    });
  }

  function mapperTreeFormatSpeciesTag(specInfo) {
    if (!specInfo) return '';
    var c = (specInfo.common || '').trim();
    var l = (specInfo.latin || '').trim();
    if (mapperTreeSpeciesNameFormat === 'latin') return l ? '(' + l + ')' : '';
    if (mapperTreeSpeciesNameFormat === 'both') {
      if (c && l && c !== l) return '(' + c + ', ' + l + ')';
      return (c || l) ? '(' + (c || l) + ')' : '';
    }
    return (c || l) ? '(' + (c || l) + ')' : '';
  }

  function mapperTreeSpeciesLabel(speciesObj) {
    if (!speciesObj) return '';
    var c = (speciesObj.common || '').trim();
    var l = (speciesObj.latin || '').trim();
    if (mapperTreeSpeciesNameFormat === 'latin') return l || c;
    if (mapperTreeSpeciesNameFormat === 'both') {
      if (c && l && c !== l) return c + ' (' + l + ')';
      return c || l;
    }
    return c || l;
  }

  function mapperTreeExtendLabelDictsForSpecies() {
    if (!window.MAPPER_TREE_SPECIES_DATA) return;
    var sld = window.mapperTreeStemLabelDicts;
    if (!sld) return;
    var allEntries = window.MAPPER_TREE_SPECIES_DATA.entry_to_species || {};
    Object.keys(allEntries).forEach(function (entryName) {
      var specInfo = allEntries[entryName];
      var stemKey = entryName.toUpperCase();
      var baseStem = entryName.split('_')[0].toUpperCase();
      var tag = mapperTreeFormatSpeciesTag(specInfo);
      var suffix = tag ? ' ' + tag : '';
      var iup = sld.IUPHAR[baseStem];
      sld.IUPHAR[stemKey] = [iup ? iup[0] + suffix : baseStem + suffix];
      var gen = sld.Gene[baseStem];
      sld.Gene[stemKey] = [gen ? gen[0] + suffix : baseStem + suffix];
      var uni = sld.UniProt[baseStem];
      sld.UniProt[stemKey] = [uni ? uni[0] + suffix : baseStem + suffix];
    });
    (window.MAPPER_TREE_SPECIES_DATA.nonhuman_only || []).forEach(function (nho) {
      var stemKey = nho.entry.toUpperCase();
      var baseStem = nho.entry.split('_')[0].toUpperCase();
      var tag = mapperTreeFormatSpeciesTag({ common: nho.common, latin: nho.latin });
      var suffix = tag ? ' ' + tag : '';
      var iup = sld.IUPHAR[baseStem];
      sld.IUPHAR[stemKey] = [iup ? iup[0] + suffix : baseStem + suffix];
      var gen = sld.Gene[baseStem];
      sld.Gene[stemKey] = [gen ? gen[0] + suffix : baseStem + suffix];
      var uni = sld.UniProt[baseStem];
      sld.UniProt[stemKey] = [uni ? uni[0] + suffix : baseStem + suffix];
    });
  }

  function mapperTreePopulateSpeciesSelect($tr, entryId, preserveSelection) {
    var $sel = $tr.find('.mapper-tree-species-select');
    if (!$sel.length || !entryId || !window.MAPPER_TREE_SPECIES_DATA) return;
    var prevVal = preserveSelection ? ($sel.val() || '').trim() : '';
    if ($sel.data('select2')) { $sel.select2('destroy'); }
    $sel.empty();
    var stem = entryId.split('_')[0];
    var orthologs = (window.MAPPER_TREE_SPECIES_DATA.by_stem[stem] || []).slice();
    if (!orthologs.length) {
      // non-human-only: collect ALL species for this stem (not just the stored entry)
      (window.MAPPER_TREE_SPECIES_DATA.nonhuman_only || []).forEach(function (nho) {
        if (nho.stem === stem) orthologs.push(nho);
      });
    }
    orthologs.forEach(function (o) {
      $('<option>').val(o.entry).text(mapperTreeSpeciesLabel(o)).appendTo($sel);
    });
    var validPrev = prevVal && orthologs.some(function (o) { return o.entry === prevVal; });
    if (validPrev) {
      $sel.val(prevVal);
    } else {
      var humanOpt = null;
      orthologs.forEach(function (o) { if (o.is_human) humanOpt = o; });
      $sel.val(humanOpt ? humanOpt.entry : (orthologs[0] ? orthologs[0].entry : ''));
    }
    $sel.select2({
      width: 'resolve',
      dropdownAutoWidth: true,
      minimumResultsForSearch: 6,
      dropdownParent: $('body')
    });
  }

  function mapperTreeInitAllSpeciesDropdowns() {
    $('#mapper-tree-input-tbody tr').each(function () {
      var $tr = $(this);
      var entry = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
      if (entry) mapperTreePopulateSpeciesSelect($tr, entry, true);
    });
  }

  function mapperTreeInjectSpeciesLeavesForStem(node, baseStem, speciesEntries) {
    var ch = node.children;
    if (!ch || !ch.length) return false;
    for (var i = 0; i < ch.length; i++) {
      var child = ch[i];
      if (!child.children || !child.children.length) {
        var leafName = String(child.name || '').trim();
        if (leafName.toLowerCase() === baseStem.toLowerCase()) {
          ch[i] = $.extend(true, {}, child, { name: baseStem + '_human' });
          var insertIdx = i + 1;
          speciesEntries.forEach(function (spEntry) {
            if (spEntry.toLowerCase() !== baseStem.toLowerCase() + '_human') {
              var newLeaf = $.extend(true, {}, child, { name: spEntry });
              ch.splice(insertIdx, 0, newLeaf);
              insertIdx++;
            }
          });
          return true;
        }
      } else {
        if (mapperTreeInjectSpeciesLeavesForStem(child, baseStem, speciesEntries)) return true;
      }
    }
    return false;
  }

  function mapperTreeBuildLiveSkelWithSpecies() {
    var stemSpeciesMap = {};
    $('#mapper-tree-input-tbody tr').each(function () {
      var $tr = $(this);
      var entry = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
      if (!entry) return;
      var specVal = ($tr.find('.mapper-tree-species-select').val() || '').trim() || entry;
      var baseStem = entry.split('_')[0];
      if (!stemSpeciesMap[baseStem]) stemSpeciesMap[baseStem] = [];
      if (stemSpeciesMap[baseStem].indexOf(specVal) === -1) stemSpeciesMap[baseStem].push(specVal);
    });
    var liveSkel = deepClone(ORIG_SKEL);
    Object.keys(stemSpeciesMap).forEach(function (baseStem) {
      var entries = stemSpeciesMap[baseStem];
      // Always rename/inject so the leaf name matches the species-suffixed circles key.
      // Human-only: renames bare leaf 5ht1a -> 5ht1a_human (no siblings added).
      mapperTreeInjectSpeciesLeavesForStem(liveSkel, baseStem, entries);
    });
    return liveSkel;
  }

  function mapperTreeFindLeafNameExact(skel, entryId, taRaw, unmatchedRaw) {
    var prefer = mapperStemFromEntry(entryId);
    if (prefer && skel) {
      var tgt = prefer.toLowerCase();
      var found = null;
      function walk(node) {
        if (!node) {
          return;
        }
        var ch = node.children;
        if (!ch || !ch.length) {
          var nm = node.name != null ? String(node.name).trim() : '';
          if (nm.replace(/_human$/i, '').toLowerCase() === tgt) {
            found = nm;
          }
          return;
        }
        ch.forEach(walk);
      }
      walk(skel);
      if (found) {
        return found;
      }
      // Fallback: try bare stem for non-human entries (e.g. 5ht5b_mouse -> leaf 5ht5b)
      var bareStem = tgt.split('_')[0];
      if (bareStem !== tgt) {
        var walkBare = function (node) {
          if (!node) return;
          var bch = node.children;
          if (!bch || !bch.length) {
            var bnm = node.name != null ? String(node.name).trim() : '';
            if (bnm.replace(/_human$/i, '').toLowerCase() === bareStem) {
              found = bnm;
            }
            return;
          }
          bch.forEach(walkBare);
        };
        walkBare(skel);
        if (found) return found;
      }
    }
    var rawCandidate = unmatchedRaw || taRaw || '';
    if (window.mapperTreeResolveEntry && typeof window.mapperTreeResolveEntry === 'function') {
      var resolved = window.mapperTreeResolveEntry(rawCandidate);
      if (resolved) {
        return mapperTreeFindLeafNameExact(skel, resolved, '', '');
      }
    } else if (window.MAPPER_CORE_RESOLVE) {
      var up = String(rawCandidate || '').trim().toUpperCase();
      var rid = window.MAPPER_CORE_RESOLVE[up];
      if (rid) {
        return mapperTreeFindLeafNameExact(skel, rid, '', '');
      }
    }
    return null;
  }

  function mapperTreeCollectNamesPayload(skel, circlesObj) {
    var map = {};
    $('#mapper-tree-input-tbody tr').each(function () {
      var $tr = $(this);
      var entry = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
      var ta = (($tr.find('.mapper-core-in-receptor').val() || '') + '').trim();
      var unmatched = ($tr.data('mapperCoreUnmatchedRaw') || '') + '';
      if (!entry && !ta && !String(unmatched).trim()) {
        return;
      }
      var effectiveEntry = entry;
      if (mapperTreeSpeciesActive && entry) {
        var specValC = ($tr.find('.mapper-tree-species-select').val() || '').trim();
        if (specValC) effectiveEntry = specValC;
      }
      var ln = mapperTreeFindLeafNameExact(skel, effectiveEntry, ta, unmatched);
      if (!ln) {
        return;
      }
      var keyUpper = mapperTreeSpeciesActive ? ln.toUpperCase() : ln.replace(/_human$/i, '').toUpperCase();
      var cir = circlesObj[keyUpper];
      if (!cir) {
        return;
      }
      map[ln] = cir;
    });
    return map;
  }

  function mapperTreeIsTextMode() {
    return window.mapperTreeInputMode === 'text';
  }

  function mapperTreeSyncFirstRowPlaceholder() {
    var hint = mapperTreeIsTextMode()
      ? MAPPER_TREE_PLACEHOLDER_CATMODE
      : MAPPER_TREE_PLACEHOLDER_NUMERIC;
    $('#mapper-tree-input-tbody tr').each(function (idx) {
      var $inp = $(this).find('.mapper-core-in-receptor');
      if (!$inp.length) { return; }
      if (idx === 0) { $inp.attr('placeholder', hint); }
      else            { $inp.removeAttr('placeholder'); }
    });
  }

  // FNV1a32 hash + default label colour moved to MapperPageCore (mapper_page_core.js).
  function mapperTreeDefaultHexForLabelKey(lbl) {
    return MapperPageCore.defaultColorForLabel(lbl);
  }

  function mapperTreeParseNumLoose(s) {
    if (s == null || String(s).trim() === '') {
      return null;
    }
    var t = String(s).trim().replace(',', '.');
    var x = parseFloat(t);
    return isFinite(x) ? x : null;
  }

  function mapperTreeBuildTreeCircles() {
    var out = {};
    var textMode = mapperTreeIsTextMode();
    $('#mapper-tree-input-tbody tr').each(function () {
      var $tr = $(this);
      if ($tr.hasClass('mapper-tree-row-over-limit')) { return; }
      var entry = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
      var ta = (($tr.find('.mapper-core-in-receptor').val() || '') + '').trim();
      var unmatched = ($tr.data('mapperCoreUnmatchedRaw') || '') + '';
      if (!entry && !ta && !String(unmatched).trim()) {
        return;
      }
      var id = entry;
      if (!id && window.MAPPER_CORE_RESOLVE) {
        id =
          window.MAPPER_CORE_RESOLVE[String(ta || unmatched).trim().toUpperCase()] || '';
      }
      if (!id) {
        return;
      }
      var uniKey;
      if (mapperTreeSpeciesActive) {
        var specValB = ($tr.find('.mapper-tree-species-select').val() || '').trim();
        uniKey = specValB ? specValB.toUpperCase() : mapperTreeBareStemUpper(id);
      } else {
        uniKey = mapperTreeBareStemUpper(id);
      }
      if (!uniKey) {
        return;
      }
      var innerRaw = (($tr.find('.mapper-tree-inner').val() || '') + '').trim();
      var o = {};
      if (!textMode) {
        var innerN = mapperTreeParseNumLoose(innerRaw);
        var o1 = mapperTreeParseNumLoose($tr.find('.mapper-tree-o1').val());
        var o2 = mapperTreeParseNumLoose($tr.find('.mapper-tree-o2').val());
        var o3 = mapperTreeParseNumLoose($tr.find('.mapper-tree-o3').val());
        var o4 = mapperTreeParseNumLoose($tr.find('.mapper-tree-o4').val());
        if (innerN == null && o1 == null && o2 == null && o3 == null && o4 == null) {
          return;
        }
        if (innerN != null) {
          o.Inner = innerN;
        }
        if (o1 != null) {
          o.Outer1 = o1;
        }
        if (o2 != null) {
          o.Outer2 = o2;
        }
        if (o3 != null) {
          o.Outer3 = o3;
        }
        if (o4 != null) {
          o.Outer4 = o4;
        }
      } else {
        if (!innerRaw) {
          return;
        }
        if (!window.MAPPER_CORE_LABEL_COLORS[innerRaw]) {
          window.MAPPER_CORE_LABEL_COLORS[innerRaw] = mapperTreeDefaultHexForLabelKey(innerRaw);
        }
        var enabled = window.MAPPER_CORE_LABEL_ENABLED[innerRaw] !== false;
        o.Inner = innerRaw;
        o.ColorValue = enabled ? window.MAPPER_CORE_LABEL_COLORS[innerRaw] || mapperTreeDefaultHexForLabelKey(innerRaw) : '#ffffff';
      }
      out[uniKey] = o;
    });
    return out;
  }

  function mapperTreeFilterSkeleton(skel, namesPayload) {
    if (!namesPayload || !Object.keys(namesPayload).length) {
      return deepClone(skel);
    }
    var k = mapperJsKeepByNames(deepClone(skel), namesPayload);
    return k;
  }

  function mergedFontSizesFromDom() {
    var def = (ORIG_OPTS && ORIG_OPTS.fontSize) || {};
    function px(id, fb) {
      var $e = $('#' + id);
      return $e.length ? $e.val() + 'px' : fb;
    }
    return {
      class: px('classFontSizeSlider', def.class || '15px'),
      ligandtype: px('ligandTypeFontSizeSlider', def.ligandtype || '14px'),
      receptorfamily: px('receptorFamilyFontSizeSlider', def.receptorfamily || '13px'),
      receptor: px('receptorFontSizeSlider', def.receptor || '12px')
    };
  }

  function parsePxFromFontSize(css) {
    var m = /(\d+)/.exec(String(css || ''));
    return m ? parseInt(m[1], 10) : 12;
  }

  function mapperTreeRadialLeafLabelGap() {
    var circleSize = Number(styling_circles && styling_circles.circle_size);
    return MAPPER_TREE_RADIAL_LEAF_LABEL_GAP_BASE + (isFinite(circleSize) ? circleSize : 3);
  }

  function mapperTreeSyncCircleSizingFromUi() {
    var csEl = $('#mapper-tree-circle-size-slider');
    if (csEl.length) {
      styling_circles.circle_size = 3 + 1 * Number(csEl.val());
    }
    var spEl = $('#mapper-tree-circle-spacer-slider');
    if (spEl.length) {
      styling_circles.circle_spacer =
        styling_circles.circle_size * (2 + 0.5 * Number(spEl.val())) + 1;
    }
  }

  function mapperTreeActiveLeafDropdownVal() {
    var $b = $('.mapper-tree-leaf-btn.btn-primary');
    var dv = ($b.attr('data-value') || '').trim();
    return dv === 'UniProt' || dv === 'Gene' ? dv : 'IUPHAR';
  }

  function mapperTreeSyncLeafLabelUi(btnVal) {
    var mapUi = {
      IUPHAR: 'Protein',
      Gene: 'Gene',
      UniProt: 'UniProt'
    };
    window.TREE_UI = window.TREE_UI || {};
    window.TREE_UI.leafLabelType = mapUi[btnVal] || 'Protein';
    $('.mapper-tree-leaf-btn').each(function () {
      var ok = ($(this).attr('data-value') || '') === btnVal;
      $(this).toggleClass('btn-outline-primary', !ok).toggleClass('btn-primary', ok);
    });
    // Refresh resolved receptor capsules in the input table to match the new label type
    $('#mapper-tree-input-tbody tr').each(function () {
      var sid = ($(this).find('.mapper-core-receptor-entry').val() || '').trim();
      if (!sid) { return; }
      var $view = $(this).find('.mapper-core-receptor-html-view');
      if ($view.length && $view.is(':visible')) {
        $view.html(mapperTreeResolvedDisplay(sid));
      }
    });
  }

  function mapperTreeActiveStemDict(btnVal) {
    if (btnVal === 'UniProt') {
      return window.mapperTreeStemLabelDicts.UniProt;
    }
    if (btnVal === 'Gene') {
      return window.mapperTreeStemLabelDicts.Gene;
    }
    return window.mapperTreeStemLabelDicts.IUPHAR;
  }

  /** No mapped rows yet: empty plot host (no skeleton draw). */
  function mapperTreeShowPlaceholderPlot(kind, msg) {
    var host = $('#tree_plot');
    if (!host.length) {
      return;
    }
    host.empty();
    var wrap = $('<div class="mapper-tree-plot-placeholder" role="region" aria-label="Tree plot placeholder"/>');
    if (kind === 'nomatch') {
      wrap.append(
        $('<p class="mapper-tree-plot-placeholder-title"/>').text(
          'Receptors present but none match this phylogeny'
        )
      );
      wrap.append($('<p class="mapper-tree-plot-placeholder-warn"/>').text(msg || 'Check receptors against the phylogenetic scaffold.'));
    } else {
      wrap.append($('<p class="mapper-tree-plot-placeholder-title"/>').text('No receptors mapped yet'));
      wrap.append(
        $('<p class="mapper-tree-plot-placeholder-hint"/>').text(
          'Pick or paste receptors at left and add Numeric values (Inner + optional O1–O4) or a Text Inner label — the phylogeny renders incrementally row by row.'
        )
      );
    }
    host.append(wrap);
  }

  function mapperTreeFitFinalSvgViewBox() {
    var svgNode = d3.select('#tree_plot svg').node();
    if (!svgNode) {
      return;
    }
    var bbox;
    try {
      bbox = svgNode.getBBox();
    } catch (err) {
      return;
    }
    if (!bbox || !isFinite(bbox.x) || !isFinite(bbox.y) || !isFinite(bbox.width) || !isFinite(bbox.height) || bbox.width <= 0 || bbox.height <= 0) {
      return;
    }
    var pad = 28;
    var minX = bbox.x - pad;
    var minY = bbox.y - pad;
    var vbWidth = bbox.width + pad * 2;
    var vbHeight = bbox.height + pad * 2;
    d3.select(svgNode)
      .attr('viewBox', minX + ' ' + minY + ' ' + vbWidth + ' ' + vbHeight)
      .attr('preserveAspectRatio', 'xMidYMid meet');
  }

  function mapperTreeSetCircleStarterFromLeafLabels() {
    var maxOuterEdge = 0;
    var labelGap = 3;
    var circleSize = Number(styling_circles && styling_circles.circle_size);
    var circleSpacer = Number(styling_circles && styling_circles.circle_spacer);
    if (!(isFinite(circleSize) && circleSize > 0)) {
      circleSize = 3;
    }
    if (!(isFinite(circleSpacer) && circleSpacer > 0)) {
      circleSpacer = 10;
    }
    d3.select('#tree_plot').selectAll('g.node[id]').each(function (d) {
      if (!d || d.depth !== window.tree_options_draw.depth) {
        return;
      }
      var textNode = d3.select(this).select('text').node();
      if (!textNode) {
        return;
      }
      var bbox;
      var matrix;
      try {
        bbox = textNode.getBBox();
        var transform = textNode.transform && textNode.transform.baseVal
          ? textNode.transform.baseVal.consolidate()
          : null;
        matrix = transform ? transform.matrix : null;
      } catch (err) {
        return;
      }
      if (!bbox || !isFinite(bbox.x) || !isFinite(bbox.width)) {
        return;
      }
      var corners = [
        { x: bbox.x, y: bbox.y },
        { x: bbox.x + bbox.width, y: bbox.y },
        { x: bbox.x, y: bbox.y + bbox.height },
        { x: bbox.x + bbox.width, y: bbox.y + bbox.height }
      ];
      corners.forEach(function (corner) {
        var x = corner.x;
        if (matrix) {
          x = matrix.a * corner.x + matrix.c * corner.y + matrix.e;
        }
        if (isFinite(x)) {
          maxOuterEdge = Math.max(maxOuterEdge, x);
        }
      });
    });
    if (maxOuterEdge > 0) {
      styling_circles.starter = Math.max(1, maxOuterEdge + circleSize + labelGap - (2 * circleSpacer));
    }
  }

  function mapperTreeRedrawNow() {
    if (!ORIG_SKEL || !ORIG_OPTS || typeof window.mapperClassificationRedraw !== 'function') {
      mapperTreeShowPlaceholderPlot('empty');
      $('#mapper-tree-messages').empty();
      return;
    }

    mapperTreeUpdateOverLimitMarks();
    $('#mapper-tree-messages').empty();

    Tree_circles = mapperTreeBuildTreeCircles();

    var circlesKeys = Tree_circles && Object.keys(Tree_circles).length;
    if (!circlesKeys) {
      mapperTreeShowPlaceholderPlot('empty');
      return;
    }

    $('#tree_plot').empty();

    var workSkel = mapperTreeSpeciesActive ? mapperTreeBuildLiveSkelWithSpecies() : ORIG_SKEL;
    var namesPayload = mapperTreeCollectNamesPayload(workSkel, Tree_circles);

    /** Values present but no leaves resolved on the scaffold (wrong / unresolved receptors). */
    if (!Object.keys(namesPayload).length) {
      mapperTreeShowPlaceholderPlot(
        'nomatch',
        'Entries have values but no receptors match a leaf on this phylogeny. Use search / picker until each row resolves cleanly.'
      );
      return;
    }

    var filteredSk = mapperTreeFilterSkeleton(workSkel, namesPayload);

    // Auto-manage class level: hide when only 1 class present, restore when 2+ appear
    var rawClassCount = (filteredSk && filteredSk.children) ? filteredSk.children.length : 0;
    if (rawClassCount >= 1) {
      if (rawClassCount <= 1 && !skipClassLevel) {
        skipClassLevel = true;
        classLevelAutoOff = true;
        mapperTreeUpdateSliderStates();
      } else if (rawClassCount > 1 && classLevelAutoOff && skipClassLevel) {
        skipClassLevel = false;
        classLevelAutoOff = false;
        mapperTreeUpdateSliderStates();
      }
    }

    if (skipClassLevel) {
      // Single-class view: keep chemotype only if the sole class is Class A (uniform leaf depth).
      // Non-Class-A single-class: strip both levels so all leaves land at uniform depth.
      var singleClassName = (filteredSk.children && filteredSk.children[0] && filteredSk.children[0].name) || '';
      if (singleClassName.indexOf('Class A') !== -1) {
        skipLigandtypeLevel = false;
        filteredSk = mapperTreeStripClassLevel(filteredSk);
      } else {
        filteredSk = mapperTreeStripBothLevels(filteredSk);
        skipLigandtypeLevel = true;
      }
    } else {
      // Multi-class view: strip all chemotype so all leaves land at a uniform depth.
      filteredSk = mapperTreeStripLigandtypeLevel(filteredSk);
      skipLigandtypeLevel = true;
    }

    var labelDdVal = mapperTreeActiveLeafDropdownVal();
    mapperTreeSyncLeafLabelUi(labelDdVal);

    if (
      filteredSk == null ||
      (filteredSk.children && filteredSk.children.length === 0) ||
      ((filteredSk.name === '' || filteredSk.name === null) &&
        (!filteredSk.children || !filteredSk.children.length))
    ) {
      mapperTreeShowPlaceholderPlot(
        'nomatch',
        'Filtered tree would be empty for this receptor subset. Adjust rows or mappings.'
      );
      return;
    }

    var btnVal = mapperTreeActiveLeafDropdownVal();

    if (mapperTreeSpeciesActive) {
      mapperTreeExtendLabelDictsForSpecies();
    }

    window.TREE_UI = window.TREE_UI || {};
    window.TREE_UI.layout = MAPPER_TREE_LAYOUT;

    var fontMerge = mergedFontSizesFromDom();
    mapperTreeSyncCircleSizingFromUi();
    var baseOptsIn = $.extend(deepClone(ORIG_OPTS), {
      fontSize: fontMerge,
      fontFamily: ORIG_OPTS.fontFamily || 'Palatino',
      radialLeafLabelGap: mapperTreeRadialLeafLabelGap(),
      leafEndDotRadius:
        ORIG_OPTS.leafEndDotRadius != null && isFinite(Number(ORIG_OPTS.leafEndDotRadius))
          ? Number(ORIG_OPTS.leafEndDotRadius)
          : 2
    });

    var wouldPromote = !!(filteredSk && filteredSk.children && filteredSk.children.length === 1);
    var promLocal = mapperTreeMaybePromoteRoot(deepClone(filteredSk), baseOptsIn);

    if (skipClassLevel || skipLigandtypeLevel) {
      var classPresent = !skipClassLevel;
      var chemotypePresent = !skipLigandtypeLevel;
      var newDepth, fontMap, newBL;

      if (!wouldPromote) {
        if (classPresent && !chemotypePresent) {
          // Root → Class → RF → Leaf
          newDepth = 3;
          fontMap = { 1: 'class', 2: 'receptorfamily' };
          newBL = { 1: (ORIG_OPTS.branch_length && ORIG_OPTS.branch_length[1]) || 'Class', 2: 'Receptor family', 3: '' };
        } else if (!classPresent && chemotypePresent) {
          // Root → Chemotype → RF → Leaf
          newDepth = 3;
          fontMap = { 1: 'ligandtype', 2: 'receptorfamily' };
          newBL = { 1: 'Chemotype', 2: 'Receptor family', 3: '' };
        } else {
          // Root → RF → Leaf (both stripped)
          newDepth = 2;
          fontMap = { 1: 'receptorfamily' };
          newBL = { 1: 'Receptor family', 2: '' };
        }
      } else {
        if (!classPresent && !chemotypePresent) {
          // RF → Leaf (both stripped, single RF promoted)
          newDepth = 1;
          fontMap = {};
          newBL = { 1: '' };
        } else {
          // Any other promoted case: RF at depth 1
          newDepth = 2;
          fontMap = { 1: 'receptorfamily' };
          newBL = { 1: 'Receptor family', 2: '' };
        }
      }

      promLocal.opts.depth = newDepth;
      promLocal.opts.fontSize_depth_map = fontMap;
      promLocal.opts.branch_length = newBL;
    }

    var td = deepClone(promLocal.tree);
    window.tree_options_draw = window.mapperClassificationRedraw(td, deepClone(promLocal.opts), window.TREE_UI.layout);

    if (mapperTreeIsTextMode()) {
      styling_circles.mode = 'Text';
      Tree_datatypes_dict = $.extend({}, window.mapperTreeDiscreteDatatypes);
    } else {
      styling_circles.mode = 'Numeric';
      Tree_datatypes_dict = $.extend({}, window.mapperTreeNumericDatatypes);
    }

    styling_circles.starter =
      (maxLeafNodeLength_scaler || 10) * parsePxFromFontSize(fontMerge.receptor || '12px');

    /** Circle size / spacer already synced before tree draw so label gap uses the same effective size. */

    var dictStem = mapperTreeActiveStemDict(btnVal);

    custom_changeLeavesLabels(
      'tree_plot',
      btnVal === 'UniProt' ? 'UniProt' : btnVal === 'Gene' ? 'Gene' : 'IUPHAR',
      dictStem,
      styling_circles
    );
    mapperTreeSetCircleStarterFromLeafLabels();

    DrawCircles('tree_plot', Tree_circles, Tree_colors, styling_circles, Tree_circle_styling_dict);

    d3.select('#tree_plot svg').selectAll('.legend-group').remove();
    if (ShowLegend) {
      if (mapperTreeIsTextMode() && typeof CreateTextLegend === 'function') {
        CreateTextLegend('tree_plot', Tree_circles, Tree_textlegend_styling);
      } else if (typeof createLegendBars === 'function') {
        createLegendBars(
          'tree_plot',
          Tree_circles,
          Tree_colors,
          Tree_circle_styling_dict,
          Tree_datatypes_dict,
          Label_dict,
          TreeLegendPosition || 'Top'
        );
      }
    }
    mapperTreeFitFinalSvgViewBox();
    mapperTreeSyncColorRings();
  }

  window.mapperTreeRedrawNow = mapperTreeRedrawNow;

  var _treeRedrawDebounced = MapperPageCore.debounce(function () { mapperTreeRedrawNow(); }, DEBOUNCE_MS);
  function mapperTreeScheduleRedraw() {
    if (suppressRedraw) {
      return;
    }
    mapperTreeUpdateOverLimitMarks();
    _treeRedrawDebounced.schedule();
  }

  // Batches the full-table cosmetic syncs triggered by the two hottest, highest-frequency
  // actions (typing in a value cell, deleting a row) so they run once per debounce window
  // instead of once per keystroke/click at large row counts. Structural row-count upkeep
  // (mapperTreeEnsureTrailingBlankRow) stays synchronous — it's already bounded to the last
  // two rows, not a full-table scan. Every other call site of these sync functions
  // (restore, demo-fill, paste, clear) is untouched and keeps calling them directly.
  var _treeCosmeticSyncDebounced = MapperPageCore.debounce(function () {
    mapperTreeRefreshInnerSwatches();
    mapperTreeSyncInnerHints();
    mapperTreeSyncRemoveButtons();
    mapperTreeCompactReceptorRowsAfterInput();
    mapperTreeSyncFirstRowPlaceholder();
  }, DEBOUNCE_MS);

  // ── Mode snapshots ──────────────────────────────────────────────────────
  var MAPPER_TREE_MODE_SNAPSHOTS = { numeric: null, categorical: null };

  function mapperTreeCaptureSnapshot() {
    var mode = mapperTreeIsTextMode() ? 'categorical' : 'numeric';
    var rows = mapperTreeSerializeRows();
    // In categorical mode, strip o1-o4 from the saved rows (they are disabled/irrelevant)
    if (mode === 'categorical') {
      rows = rows.map(function (r) {
        return { entry: r.entry, receptorText: r.receptorText, unmatched: r.unmatched,
                 invalid: r.invalid, inner: r.inner, o1: '', o2: '', o3: '', o4: '',
                 speciesEntry: r.speciesEntry };
      });
      MAPPER_TREE_MODE_SNAPSHOTS.categorical = {
        rows: rows,
        labelColors: $.extend({}, window.MAPPER_CORE_LABEL_COLORS || {})
      };
    } else {
      MAPPER_TREE_MODE_SNAPSHOTS.numeric = { rows: rows };
    }
  }

  function mapperTreeSeedOppositeMode(targetMode) {
    if (MAPPER_TREE_MODE_SNAPSHOTS[targetMode] != null) { return; }
    if (targetMode === 'categorical') {
      MAPPER_TREE_MODE_SNAPSHOTS.categorical = { rows: [], labelColors: {} };
    } else {
      MAPPER_TREE_MODE_SNAPSHOTS.numeric = { rows: [] };
    }
  }

  function mapperTreeApplyLeftPanelWidth() {
    var text = mapperTreeIsTextMode();
    var speciesExtra = mapperTreeSpeciesActive ? 88 : 0;
    $('.mapper-core-wheel-wrap.mapper-tree-page .mapper-core-wheel-left').css(
      'width', text ? (400 + speciesExtra) + 'px' : ''
    );
    $('.mapper-core-wheel-wrap.mapper-tree-page').toggleClass('mapper-tree-species-active', mapperTreeSpeciesActive);
  }

  function mapperTreeSyncClearDropdown() {
    var isText = mapperTreeIsTextMode();
    $('#mapper-tree-clear-mode-label').text(isText ? 'Categories' : 'Numbers');
    $('.mapper-tree-clear-col-num').toggle(!isText);
    $('.mapper-tree-clear-col-text').toggle(isText);
  }

  function mapperTreeSetInputMode(mode) {
    var nextMode = mode === 'text' ? 'text' : 'numeric';
    var prevMode = window.mapperTreeInputMode;

    // No-op if mode unchanged, except on initial boot (prevMode undefined)
    if (prevMode === nextMode) { return; }

    // Save current state before switching (skip on initial boot)
    if (prevMode === 'text' || prevMode === 'numeric') {
      mapperTreeCaptureSnapshot();
    }

    window.mapperTreeInputMode = nextMode;
    var text = mapperTreeIsTextMode();

    // Seed the target mode's snapshot if this is the first visit
    mapperTreeSeedOppositeMode(text ? 'categorical' : 'numeric');

    // Restore rows from snapshot — always rebuild, even when the target mode has
    // never been visited (snap.rows is []), so a genuinely empty mode shows as empty
    // rather than leaving the previous mode's rows sitting in the DOM.
    var snap = MAPPER_TREE_MODE_SNAPSHOTS[text ? 'categorical' : 'numeric'] || {};
    if (text) {
      window.MAPPER_CORE_LABEL_COLORS = snap.labelColors ? $.extend({}, snap.labelColors) : {};
    }
    mapperTreeApplySerializedRows(snap.rows || []);

    // Update button active states
    $('#mapper-tree-mode-numeric-btn, #mapper-tree-mode-labels-btn').each(function () {
      var isNum = $(this).attr('id') === 'mapper-tree-mode-numeric-btn';
      var on = text ? !isNum : isNum;
      $(this).toggleClass('active', on).attr('aria-pressed', on ? 'true' : 'false');
    });

    var $tbl = $('#mapper-tree-input-table');
    $tbl.toggleClass('mapper-core-text-mode', text);

    // Disable hidden outer inputs (tab / screen reader skip)
    $tbl.find('.mapper-tree-o1,.mapper-tree-o2,.mapper-tree-o3,.mapper-tree-o4').prop('disabled', text);

    // Swap inner column header label
    $('.mapper-tree-inner-header-label').text(text ? 'Category' : 'Dot 1');

    // Show/hide Colors and Labels dropdown panels
    $('#mapper-tree-colors-numeric, #mapper-tree-labels-numeric').toggle(!text);
    $('#mapper-tree-colors-text,   #mapper-tree-labels-text').toggle(text);

    // Hide Legend labels button and spacer slider in categorical mode; shrink panel
    $('#mapper-tree-labels-dropdown').toggle(!text);
    $('#mapper-tree-spacer-section').toggle(!text);
    $('.mapper-core-wheel-wrap.mapper-tree-page').toggleClass('mapper-tree-catmode', text);
    // Drive the left panel to an explicit pixel width so the flex column actually shrinks
    mapperTreeApplyLeftPanelWidth();

    // Resize colors menu and sync categorical colours panel + placeholder
    $('#mapper-tree-colors-menu').toggleClass('mapper-tree-colors-text-mode', text);
    mapperTreeSyncColorsPanel();
    mapperTreeSyncFirstRowPlaceholder();

    // Orange hint on Category cells that need filling
    mapperTreeSyncInnerHints();

    mapperTreeRefreshInnerSwatches();
    mapperTreeScheduleRedraw();
    if (typeof mapperTreeSyncClearDropdown === 'function') { mapperTreeSyncClearDropdown(); }
  }

  function mapperTreeSyncInnerHints() {
    var text = mapperTreeIsTextMode();
    $('#mapper-tree-input-tbody tr').each(function () {
      var $tr = $(this);
      var $td = $tr.find('td.mapper-tree-cell-inner');
      if (!$td.length) { return; }
      if (!text) {
        $td.removeClass('mapper-tree-inner-hint');
        return;
      }
      var hasReceptor = !!($tr.find('.mapper-core-receptor-entry').val() || '').trim() ||
                        !!($tr.find('.mapper-core-in-receptor').val() || '').trim();
      var hasInner = !!($tr.find('.mapper-tree-inner').val() || '').trim();
      $td.toggleClass('mapper-tree-inner-hint', hasReceptor && !hasInner);
    });
  }

  function mapperTreeDestroyRowColorSpectrum($picker) {
    if (!$picker || !$picker.length || !$.fn.spectrum) {
      return;
    }
    try {
      if ($picker.data('spectrum.id') != null || $picker.hasClass('sp-replaced')) {
        $picker.spectrum('destroy');
      }
    } catch (e2) {}
  }

  function mapperTreeApplyLabelColor(label, color, $activePicker) {
    var key = String(label || '').trim();
    if (!key || !color) {
      return;
    }
    window.MAPPER_CORE_LABEL_COLORS[key] = color;
    $('#mapper-tree-input-tbody tr').each(function () {
      var $tr = $(this);
      if (($tr.find('.mapper-tree-inner').val() || '').trim() !== key) {
        return;
      }
      var $sw = $tr.find('.mapper-core-row-color-picker');
      $sw.val(color).css('background-color', color);
      if ($activePicker && $activePicker.length && $sw[0] === $activePicker[0]) {
        return;
      }
      if ($.fn.spectrum && ($sw.data('spectrum.id') != null || $sw.hasClass('sp-replaced'))) {
        try {
          $sw.spectrum('set', color);
        } catch (e3) {}
      }
    });
    // Sync the categorical Colors panel picker for this label
    var $pp = $('#mapper-tree-label-color-pickers .mapper-core-lcat-spectrum[data-mapper-tree-lcat-label="' + key + '"]');
    if ($pp.length && (!$activePicker || !$activePicker.length || $pp[0] !== $activePicker[0])) {
      if ($.fn.spectrum && ($pp.data('spectrum') || $pp.hasClass('sp-replaced'))) {
        try { $pp.spectrum('set', color); } catch (e4) {}
      }
    }
    mapperTreeScheduleRedraw();
  }

  function mapperTreeEnsureRowColorSpectrum($picker, label, color) {
    if (!$picker.length || !$.fn.spectrum) {
      $picker.css('background-color', color || '#f5f5f5');
      return;
    }
    var currentLabel = $picker.attr('data-mapper-tree-label') || '';
    if ($picker.data('spectrum.id') != null || $picker.hasClass('sp-replaced')) {
      if (currentLabel === label) {
        $picker.spectrum('set', color || '#f5f5f5');
        return;
      }
      mapperTreeDestroyRowColorSpectrum($picker);
    }
    $picker.attr('data-mapper-tree-label', label || '');
    $picker.val(color || '#f5f5f5').css('background-color', color || '#f5f5f5');
    $picker.spectrum({
      color: color || '#f5f5f5',
      preferredFormat: 'hex',
      showInput: true,
      showPalette: true,
      showSelectionPalette: true,
      clickoutFiresChange: true,
      containerClassName: 'mapper-core-row-color-spectrum',
      replacerClassName: 'mapper-tree-row-swatch-replacer',
      palette: [
        ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd'],
        ['#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'],
        ['#000000', '#666666', '#aaaaaa', '#ffffff']
      ],
      move: function (tiny) {
        mapperTreeApplyLabelColor(label, tiny ? tiny.toHexString() : color, $picker);
      },
      change: function (tiny) {
        mapperTreeApplyLabelColor(label, tiny ? tiny.toHexString() : color, $picker);
      }
    });
  }

  function mapperTreeRefreshInnerSwatches() {
    if (!mapperTreeIsTextMode()) {
      $('#mapper-tree-input-tbody .mapper-core-row-color-picker').each(function () {
        mapperTreeDestroyRowColorSpectrum($(this));
      });
      $('#mapper-tree-input-tbody .mapper-tree-label-swatch').css('background-color', 'transparent');
      return;
    }
    $('#mapper-tree-input-tbody tr').each(function () {
      var $sw = $(this).find('.mapper-tree-label-swatch');
      var inner = (($(this).find('.mapper-tree-inner').val() || '') + '').trim();
      if (!inner) {
        mapperTreeDestroyRowColorSpectrum($sw);
        $sw.css('background-color', '#f5f5f5');
        return;
      }
      var hex = window.MAPPER_CORE_LABEL_COLORS[inner] || mapperTreeDefaultHexForLabelKey(inner);
      $sw.css('background-color', hex);
      mapperTreeEnsureRowColorSpectrum($sw, inner, hex);
    });
    mapperTreeSyncColorsPanel();
  }

  function mapperTreeSyncColorsPanel() {
    var $mount = $('#mapper-tree-label-color-pickers');
    if (!$mount.length) { return; }

    if (!mapperTreeIsTextMode()) {
      $mount.find('.mapper-core-lcat-spectrum').each(function () {
        try { if ($(this).data('spectrum')) { $(this).spectrum('destroy'); } } catch (eIgnore) {}
      });
      $mount.empty();
      mapperTreeLastLabelSet = '';
      return;
    }

    // Collect distinct labels from the input table
    var seen = {};
    var labels = [];
    $('#mapper-tree-input-tbody tr').each(function () {
      var lbl = ($(this).find('.mapper-tree-inner').val() || '').trim();
      if (lbl && !seen[lbl]) { seen[lbl] = true; labels.push(lbl); }
    });
    labels.sort(function (a, b) {
      return a.localeCompare(b, undefined, { numeric: true, sensitivity: 'base' });
    });

    var newSetKey = labels.join('\n');
    if (newSetKey === mapperTreeLastLabelSet) {
      // Label set unchanged — colours are kept in sync by mapperTreeApplyLabelColor; nothing to rebuild
      return;
    }

    // Label set changed: full rebuild — destroy existing widgets first
    $mount.find('.mapper-core-lcat-spectrum').each(function () {
      try { if ($(this).data('spectrum')) { $(this).spectrum('destroy'); } } catch (eIgnore) {}
    });
    $mount.empty();

    // Resize the dropdown menu for the current mode
    $('#mapper-tree-colors-menu').toggleClass('mapper-tree-colors-text-mode', true);

    // Initialize colour and enabled state for each label
    labels.forEach(function (lbl) {
      if (!window.MAPPER_CORE_LABEL_COLORS[lbl]) {
        window.MAPPER_CORE_LABEL_COLORS[lbl] = mapperTreeDefaultHexForLabelKey(lbl);
      }
      if (!Object.prototype.hasOwnProperty.call(window.MAPPER_CORE_LABEL_ENABLED, lbl)) {
        window.MAPPER_CORE_LABEL_ENABLED[lbl] = true;
      }
    });

    mapperTreeLastLabelSet = labels.join('\n');

    if (!labels.length) {
      $mount.append($('<p style="color:#888; font-size:11px; text-align:center; margin:4px 0;">Enter categorical labels to configure colours.</p>'));
      return;
    }

    // ── Build panel (matches wheel's mapperWheelRebuildLabelColorCustomizePanel) ──

    function domId(ix, kind) {
      return 'mapper_tree_l_' + kind + '_' + String(ix);
    }

    function updateMasterCheckbox() {
      var allOn  = labels.length > 0 && labels.every(function (l) { return window.MAPPER_CORE_LABEL_ENABLED[l] !== false; });
      var allOff = labels.length > 0 && labels.every(function (l) { return window.MAPPER_CORE_LABEL_ENABLED[l] === false; });
      var mx = $('#mapper-tree-label-cat-master')[0];
      if (!mx) { return; }
      mx.checked = !!allOn;
      mx.indeterminate = !allOn && !allOff && labels.length > 0;
    }

    function refreshGridRowUi(ix, lbl) {
      var $inp = $('#' + domId(ix, 'spe'));
      if (!$inp.length || !$inp.data('spectrum')) { return; }
      var en = window.MAPPER_CORE_LABEL_ENABLED[lbl] !== false;
      try {
        if (en) { $inp.spectrum('enable').css('opacity', '1'); }
        else    { $inp.spectrum('disable').css('opacity', '0.5'); }
      } catch (eR) {}
    }

    var $hdr = $('<div class="mapper-core-label-cat-head">').append(
      $('<input type="checkbox" id="mapper-tree-label-cat-master" aria-label="Toggle all label colours enabled">'),
      $('<div class="mapper-core-label-cat-head-title">Categories</div>')
    );
    var $grid = $('<div class="color-grid" id="mapper-tree-label-color-grid">');
    $mount.append($hdr, $grid);

    updateMasterCheckbox();

    $('#mapper-tree-label-cat-master').on('change.mapperTreeColors', function () {
      var checked = !!$(this).prop('checked');
      labels.forEach(function (lbl, ix) {
        window.MAPPER_CORE_LABEL_ENABLED[lbl] = checked;
        refreshGridRowUi(ix, lbl);
      });
      mapperTreeRefreshInnerSwatches();
      mapperTreeScheduleRedraw();
    });

    labels.forEach(function (lbl, ix) {
      var chkId   = domId(ix, 'chk');
      var spId    = domId(ix, 'spe');
      var hex = window.MAPPER_CORE_LABEL_COLORS[lbl] || mapperTreeDefaultHexForLabelKey(lbl);
      window.MAPPER_CORE_LABEL_COLORS[lbl] = hex;
      var en = window.MAPPER_CORE_LABEL_ENABLED[lbl] !== false;

      var $spe = $('<input type="text">')
        .addClass('mapper-core-lcat-spectrum form-control input-sm')
        .attr('id', spId)
        .attr('data-mapper-tree-lcat-label', lbl);

      $grid.append($('<div class="color-item">').append(
        $('<input type="checkbox">').attr('id', chkId).prop('checked', !!en),
        $('<label>').attr('for', chkId).addClass('color-label').text(lbl),
        $spe
      ));

      $('#' + chkId).on('change.mapperTreeColors', function () {
        window.MAPPER_CORE_LABEL_ENABLED[lbl] = !!$(this).prop('checked');
        refreshGridRowUi(ix, lbl);
        updateMasterCheckbox();
        mapperTreeRefreshInnerSwatches();
        mapperTreeScheduleRedraw();
      });

      $spe.spectrum({
        color: hex,
        showPalette: true,
        showInput: true,
        showButtons: false,
        preferredFormat: 'hex',
        appendTo: '#mapper-tree-colors-menu',
        containerClassName: 'mapper-core-lcat-sp-container',
        replacerClassName: 'mapper-core-lcat-replacer',
        palette: [
          ['#000', '#FF0000', '#00FF00', '#0000FF', '#FFFF00'],
          ['#FF00FF', '#00FFFF', '#FFFFFF', '#C0C0C0', '#808080']
        ],
        change: function (c) {
          if (window.MAPPER_CORE_LABEL_ENABLED[lbl] === false) { return; }
          var nextHex = c && c.toHexString ? c.toHexString() : hex;
          mapperTreeApplyLabelColor(lbl, nextHex, $spe);
        }
      });

      if (!en) {
        try { $spe.spectrum('disable').css('opacity', '0.5'); } catch (e) {}
      }
    });
  }

  function mapperTreeResolvedDisplay(entryId) {
    var sid = entryId != null ? String(entryId).trim() : '';
    var meta = (window.MAPPER_CORE_ENTRY_META && window.MAPPER_CORE_ENTRY_META[sid]) || {};
    var ltype = mapperTreeActiveLeafDropdownVal();
    if (ltype === 'Gene' && meta.gene) {
      return $('<span/>').text(meta.gene).html();
    }
    if (ltype === 'UniProt' && meta.uniprot) {
      return $('<span/>').text(meta.uniprot).html();
    }
    return meta.name_html ? String(meta.name_html) : $('<span/>').text(sid || '').html();
  }

  function mapperTreeDestroyAc($inp) {
    try {
      if ($inp.hasClass('ui-autocomplete-input')) {
        $inp.autocomplete('destroy');
      }
    } catch (e1) {}
  }

  function mapperTreeSyncClearBtn($tr) {
    var entry = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
    var typed = (($tr.find('.mapper-core-in-receptor').val() || '') + '').trim();
    var um = $tr.data('mapperCoreUnmatchedRaw');
    var has = !!(entry || typed || (um != null && String(um).trim()));
    $tr.find('.mapper-core-receptor-input-wrap').toggleClass('is-empty', !has);
  }

  function mapperTreeFocusInnerCell($tr) {
    window.setTimeout(function () {
      var $inner = $tr.find('.mapper-tree-inner:visible').first();
      if ($inner.length) {
        $inner.focus().select();
      }
    }, 0);
  }

  function mapperTreeSyncRemoveButtons() {
    $('#mapper-tree-input-tbody tr').each(function () {
      var $tr = $(this);
      $tr.find('.mapper-core-remove-cell').toggleClass('is-remove-hidden', mapperTreeRowBlank($tr));
    });
  }

  function mapperTreeCompactReceptorRowsAfterInput() {
    var hasContent = false;
    $('#mapper-tree-input-tbody tr').each(function () {
      if (!mapperTreeRowBlank($(this))) {
        hasContent = true;
        return false;
      }
    });
    $('#mapper-tree-input-table').toggleClass('mapper-core-receptors-compact', hasContent);
  }

  function mapperTreeSetResolved($tr, id) {
    var sid = id != null ? String(id).trim() : '';
    var $inp = $tr.find('.mapper-core-in-receptor');
    var $hid = $tr.find('.mapper-core-receptor-entry');
    var $view = $tr.find('.mapper-core-receptor-html-view');
    mapperTreeDestroyAc($inp);
    if (!sid) {
      $hid.val('');
      $view.hide().empty();
      $inp.val('').show();
      mapperTreeBindAc($inp);
      mapperTreeSyncClearBtn($tr);
      mapperTreeSyncRemoveButtons();
      mapperTreeCompactReceptorRowsAfterInput();
      mapperTreeScheduleRedraw();
      return;
    }
    $hid.val(sid);
    $view.html(mapperTreeResolvedDisplay(sid));
    $inp.val('').hide();
    $view.show();
    mapperTreeBindAc($inp);
    $tr.removeClass('mapper-core-row-invalid');
    $tr.removeData('mapperCoreUnmatchedRaw');
    mapperTreeSyncClearBtn($tr);
    mapperTreeSyncRemoveButtons();
    mapperTreeCompactReceptorRowsAfterInput();
    if (mapperTreeSpeciesActive && sid) {
      mapperTreePopulateSpeciesSelect($tr, sid);
    }
    if (!suppressRedraw) {
      mapperTreeFocusInnerCell($tr);
    }
    mapperTreeScheduleRedraw();
  }

  /** Returns the label-type-appropriate seed text for a resolved receptor's edit field. */
  function mapperTreeEditSeed(entryId) {
    var sid = entryId != null ? String(entryId).trim() : '';
    var meta = (window.MAPPER_CORE_ENTRY_META && window.MAPPER_CORE_ENTRY_META[sid]) || {};
    var ltype = mapperTreeActiveLeafDropdownVal();
    if (ltype === 'Gene' && meta.gene) { return meta.gene; }
    if (ltype === 'UniProt' && meta.uniprot) { return meta.uniprot; }
    return meta.name_plain || '';
  }

  function mapperTreeFilterLocal(term, limit) {
    var t = (term || '').trim().toUpperCase();
    if (!t || !window.receptorSelect2Data) {
      return [];
    }
    var ltype = mapperTreeActiveLeafDropdownVal();
    var meta = window.MAPPER_CORE_ENTRY_META || {};
    var results = [];
    window.receptorSelect2Data.forEach(function (item) {
      var m = meta[item.id] || {};
      // Always search across all name forms so any alias finds the receptor
      var searchIn = [(m.name_plain || item.name_plain || item.text || ''), (m.gene || ''), (m.uniprot || ''), item.id].join(' ').toUpperCase();
      if (searchIn.indexOf(t) === -1) { return; }
      var dispHtml;
      if (ltype === 'Gene' && m.gene) {
        dispHtml = $('<span/>').text(m.gene).html();
      } else if (ltype === 'UniProt' && m.uniprot) {
        dispHtml = $('<span/>').text(m.uniprot).html();
      } else {
        dispHtml = m.name_html || item.name_html || $('<span/>').text(m.name_plain || item.text || item.id).html();
      }
      results.push({
        label: m.name_plain || item.name_plain || item.text || item.id,
        value: item.id,
        id: item.id,
        html: dispHtml,
        name_html: item.name_html || '',
        name_plain: item.name_plain || ''
      });
    });
    var lim = limit || 80;
    if (results.length > lim) { results = results.slice(0, lim); }
    return results;
  }

  function mapperTreeBindAc($inp) {
    mapperTreeDestroyAc($inp);
    $inp.autocomplete({
      minLength: 1,
      source: function (request, response) {
        response(mapperTreeFilterLocal(request.term, 80));
      },
      focus: function () {
        return false;
      },
      select: function (event, ui) {
        mapperTreeSetResolved($inp.closest('tr'), ui.item.id);
        event.preventDefault();
      }
    });
    var w = $inp.data('ui-autocomplete');
    if (w) {
      w._renderItem = function (ul, item) {
        var inner = item.html || item.name_html || $('<span/>').text(item.label || '').html();
        return $('<li>')
          .append($('<div class="mapper-core-ac-item-label">').html(inner))
          .appendTo(ul);
      };
      if (w.menu && w.menu.element) {
        w.menu.element.addClass('mapper-tree-ac-menu');
      }
    }

    $inp.on('keyup', function () {
      var $tr = $inp.closest('tr');
      if (!$tr.find('.mapper-core-receptor-entry').val()) {
        $tr.removeData('mapperCoreUnmatchedRaw');
        $tr.removeClass('mapper-core-row-invalid');
      }
      mapperTreeSyncClearBtn($tr);
    });
    $inp.on('blur.mapper-tree', function () {
      var $tr = $inp.closest('tr');
      window.setTimeout(function () {
        if (!$inp.is(':visible')) {
          return;
        }
        var raw = ($inp.val() || '').trim();
        if (!raw || $tr.find('.mapper-core-receptor-entry').val()) {
          return;
        }
        if (window.mapperTreeResolveEntry) {
          var canon = window.mapperTreeResolveEntry(raw);
          if (canon) {
            mapperTreeSetResolved($tr, canon);
          } else if (window.receptorSelect2Data) {
            var hit = window.receptorSelect2Data.filter(function (x) {
              return String(x.id) === raw;
            });
            if (hit.length) {
              mapperTreeSetResolved($tr, hit[0].id);
              return;
            }
            var up = raw.toUpperCase();
            if (window.MAPPER_CORE_RESOLVE && window.MAPPER_CORE_RESOLVE[up]) {
              mapperTreeSetResolved($tr, window.MAPPER_CORE_RESOLVE[up]);
              return;
            }
            var known = window.receptorSelect2Data.some(function (x) {
              return String(x.id) === raw;
            });
            if (!known && raw) {
              $tr.addClass('mapper-core-row-invalid');
              $tr.data('mapperCoreUnmatchedRaw', raw);
            }
          }
        }
        mapperTreeSyncClearBtn($tr);
      }, 170);
    });
  }

  function mapperTreeCreateReceptorTd($td) {
    var $hid = $('<input type="hidden" class="mapper-core-receptor-entry" value="">');
    var $wrap = $('<div class="mapper-core-receptor-input-wrap is-empty">');
    var $inp = $(
      '<textarea class="form-control input-sm mapper-core-in-receptor" rows="1" autocomplete="off" spellcheck="false"></textarea>'
    );
    var $view = $('<div class="form-control input-sm mapper-core-receptor-html-view" tabindex="0"></div>');
    var $clr = $('<button type="button" class="mapper-core-receptor-clear" aria-label="Clear receptor">&times;</button>');
    $wrap.append($inp, $view, $clr);
    $td.append($hid, $wrap);
    $view.hide();
    $clr.on('click', function () {
      mapperTreeSetResolved($clr.closest('tr'), '');
    });
    return $inp;
  }

  function mapperTreeRowBlank($tr) {
    var entry = ($tr.find('.mapper-core-receptor-entry').val() || '').trim();
    var inner = (($tr.find('.mapper-tree-inner').val() || '') + '').trim();
    var typed = (($tr.find('.mapper-core-in-receptor').val() || '') + '').trim();
    var unmatched = ($tr.data('mapperCoreUnmatchedRaw') || '') + '';
    var outerHas = ['.mapper-tree-o1', '.mapper-tree-o2', '.mapper-tree-o3', '.mapper-tree-o4'].some(function (sel) {
      return (($tr.find(sel).val() || '') + '').trim() !== '';
    });
    return !entry && !inner && !typed && !String(unmatched).trim() && !outerHas;
  }

  function mapperTreeDestroyRowAc($tr) {
    mapperTreeDestroyAc($tr.find('.mapper-core-in-receptor'));
    $tr.find('.mapper-core-row-color-picker').each(function () {
      mapperTreeDestroyRowColorSpectrum($(this));
    });
  }

  function mapperTreeAppendRow(skipTrail) {
    var tr = $('<tr>');
    tr.append(
      $('<td class="mapper-core-remove-cell is-remove-hidden">').append(
        $('<button type="button" class="mapper-core-remove-row" aria-label="Remove row">&times;</button>')
      )
    );
    var $tdR = $('<td class="mapper-core-receptor-cell">');
    var $inp = mapperTreeCreateReceptorTd($tdR);
    tr.append($tdR);
    tr.append(
      $('<td class="mapper-tree-species-cell">').append(
        $('<select class="form-control input-sm mapper-tree-species-select"></select>')
      ).css('display', mapperTreeSpeciesActive ? 'table-cell' : 'none')
    );
    tr.append(
      $('<td class="mapper-tree-cell-inner mapper-core-value-cell">').append(
        $('<input type="text" class="form-control input-sm mapper-tree-inner" autocomplete="off">')
      )
    );
    $.each(['o1', 'o2', 'o3', 'o4'], function (_, suf) {
      var $outerInput = $('<input type="text" class="form-control input-sm mapper-tree-' + suf + '" autocomplete="off">');
      if (mapperTreeIsTextMode()) {
        $outerInput.prop('disabled', true);
      }
      tr.append(
        $('<td class="mapper-tree-outer-cell mapper-core-value-cell">').append($outerInput)
      );
    });
    tr.append(
      $('<td class="mapper-tree-swatch-cell">').append(
        '<input type="text" class="mapper-tree-label-swatch mapper-core-row-color-picker" readonly="readonly" aria-label="Label color">'
      )
    );
    $('#mapper-tree-input-tbody').append(tr);
    mapperTreeBindAc($inp);
    mapperTreeSyncClearBtn(tr);
    mapperTreeSyncRemoveButtons();
    if (!skipTrail) {
      mapperTreeEnsureTrailingBlankRow();
    }
  }

  function mapperTreeEnsureTrailingBlankRow() {
    var $tb = $('#mapper-tree-input-tbody');
    while ($tb.children().length >= 2) {
      var $last = $tb.children().last();
      var $prev = $last.prev();
      if (mapperTreeRowBlank($last) && mapperTreeRowBlank($prev)) {
        mapperTreeDestroyRowAc($last);
        $last.remove();
        continue;
      }
      break;
    }
    if (!$tb.children().length) {
      mapperTreeAppendRow(true);
    }
    var $lastOne = $tb.children().last();
    if (!mapperTreeRowBlank($lastOne) && $tb.children().length < MAPPER_TREE_MAX_ROWS) {
      mapperTreeAppendRow(true);
    }
    mapperTreeSyncRemoveButtons();
  }

  function mapperTreePasteSplit(line) {
    var parts = line.split(/\t+/);
    if (parts.length > 2) {
      return { r: (parts[0] || '').trim(), rest: parts.slice(1) };
    }
    if (parts.length === 2) {
      return { r: (parts[0] || '').trim(), rest: [(parts[1] || '').trim()] };
    }
    if (line.indexOf(';') !== -1) {
      var si = line.indexOf(';');
      return { r: line.slice(0, si).trim(), rest: line.slice(si + 1).split(';').map(function (x) { return x.trim(); }) };
    }
    return { r: (parts[0] || '').trim(), rest: [] };
  }

  function mapperTreeNormalizeSortText(value) {
    var text = $('<div/>').html(String(value || '')).text();
    return text
      .replace(/α|Α/g, 'a')
      .replace(/β|Β/g, 'b')
      .replace(/γ|Γ/g, 'g')
      .replace(/δ|Δ/g, 'd')
      .replace(/κ|Κ/g, 'k')
      .replace(/μ|Μ/g, 'm')
      .replace(/&[a-z]+;/gi, '')
      .replace(/<[^>]*>/g, '')
      .replace(/\s+/g, ' ')
      .trim()
      .toLowerCase();
  }

  function mapperTreeNaturalCompare(a, b) {
    if (!mapperTreeNaturalCompare.collator && window.Intl && Intl.Collator) {
      mapperTreeNaturalCompare.collator = new Intl.Collator(undefined, { numeric: true, sensitivity: 'base' });
    }
    if (mapperTreeNaturalCompare.collator) {
      return mapperTreeNaturalCompare.collator.compare(a, b);
    }
    return a < b ? -1 : a > b ? 1 : 0;
  }

  var mapperTreeSortState = { col: null, dir: 'asc' };

  function mapperTreeUpdateSortHeaders() {
    $('#mapper-tree-input-table th.mapper-core-sortable-head').each(function () {
      var $th = $(this);
      var col = $th.attr('data-mapper-tree-sort-col');
      var active = mapperTreeSortState.col === col;
      var dir = active ? mapperTreeSortState.dir : null;
      $th.attr('aria-sort', active ? (dir === 'asc' ? 'ascending' : 'descending') : 'none');
      $th.find('.mapper-core-sort-indicator').text(active ? (dir === 'asc' ? '↑' : '↓') : '↕');
    });
  }

  function mapperTreeSerializeRows() {
    var rows = [];
    $('#mapper-tree-input-tbody tr').each(function () {
      var $tr = $(this);
      rows.push({
        entry: ($tr.find('.mapper-core-receptor-entry').val() || '').trim(),
        receptorText: (($tr.find('.mapper-core-in-receptor').val() || '') + '').trim(),
        unmatched: ($tr.data('mapperCoreUnmatchedRaw') || '') + '',
        invalid: $tr.hasClass('mapper-core-row-invalid'),
        inner: (($tr.find('.mapper-tree-inner').val() || '') + '').trim(),
        o1: (($tr.find('.mapper-tree-o1').val() || '') + '').trim(),
        o2: (($tr.find('.mapper-tree-o2').val() || '') + '').trim(),
        o3: (($tr.find('.mapper-tree-o3').val() || '') + '').trim(),
        o4: (($tr.find('.mapper-tree-o4').val() || '') + '').trim(),
        speciesEntry: ($tr.find('.mapper-tree-species-select').val() || '').trim()
      });
    });
    return rows;
  }

  function mapperTreeRowHasContent(row) {
    return !!(
      row.entry ||
      row.receptorText ||
      String(row.unmatched || '').trim() ||
      row.inner ||
      row.o1 ||
      row.o2 ||
      row.o3 ||
      row.o4
    );
  }

  function mapperTreeSortKey(row, col) {
    if (col === 'inner') {
      var num = mapperTreeParseNumLoose(row.inner);
      return num == null ? mapperTreeNormalizeSortText(row.inner) : num;
    }
    if (row.entry) {
      return mapperTreeNormalizeSortText(mapperTreeResolvedDisplay(row.entry));
    }
    return mapperTreeNormalizeSortText(row.receptorText || row.unmatched || '');
  }

  function mapperTreePopulateRow($tr, row) {
    if (row.entry) {
      mapperTreeSetResolved($tr, row.entry);
      if (row.speciesEntry && mapperTreeSpeciesActive) {
        var $specSel = $tr.find('.mapper-tree-species-select');
        $specSel.val(row.speciesEntry);
        if ($specSel.data('select2')) { $specSel.trigger('change.select2'); }
      }
    } else {
      var $inp = $tr.find('.mapper-core-in-receptor');
      $inp.val(row.receptorText || row.unmatched || '');
      if (row.unmatched) {
        $tr.data('mapperCoreUnmatchedRaw', row.unmatched);
      }
      $tr.toggleClass('mapper-core-row-invalid', !!row.invalid);
      mapperTreeSyncClearBtn($tr);
    }
    $tr.find('.mapper-tree-inner').val(row.inner || '');
    $tr.find('.mapper-tree-o1').val(row.o1 || '');
    $tr.find('.mapper-tree-o2').val(row.o2 || '');
    $tr.find('.mapper-tree-o3').val(row.o3 || '');
    $tr.find('.mapper-tree-o4').val(row.o4 || '');
  }

  function mapperTreeApplySerializedRows(rows) {
    suppressRedraw = true;
    mapperTreeDestroyAllRows();
    $('#mapper-tree-input-tbody').empty();
    rows.forEach(function (row) {
      mapperTreeAppendRow(true);
      mapperTreePopulateRow($('#mapper-tree-input-tbody tr').last(), row);
    });
    mapperTreeEnsureTrailingBlankRow();
    mapperTreeRefreshInnerSwatches();
    mapperTreeCompactReceptorRowsAfterInput();
    mapperTreeUpdateSortHeaders();
    suppressRedraw = false;
  }

  function mapperTreeSortSerializedRows(col) {
    if (mapperTreeSortState.col === col) {
      mapperTreeSortState.dir = mapperTreeSortState.dir === 'asc' ? 'desc' : 'asc';
    } else {
      mapperTreeSortState.col = col;
      mapperTreeSortState.dir = 'asc';
    }
    var rows = mapperTreeSerializeRows();
    var filled = rows.filter(mapperTreeRowHasContent);
    filled.sort(function (a, b) {
      var ak = mapperTreeSortKey(a, col);
      var bk = mapperTreeSortKey(b, col);
      var cmp;
      if (typeof ak === 'number' && typeof bk === 'number') {
        cmp = ak - bk;
      } else {
        cmp = mapperTreeNaturalCompare(String(ak), String(bk));
      }
      return mapperTreeSortState.dir === 'desc' ? -cmp : cmp;
    });
    mapperTreeApplySerializedRows(filled);
  }

  function mapperTreeDestroyAllRows() {
    $('#mapper-tree-input-tbody tr').each(function () {
      mapperTreeDestroyRowAc($(this));
    });
  }

  function mapperTreeFocusReceptorCell($tr) {
    var $inp = $tr.find('.mapper-core-in-receptor:visible').first();
    if ($inp.length) {
      $inp.focus();
      return;
    }
    $tr.find('.mapper-core-receptor-html-view:visible').first().focus();
  }

  function mapperTreeFocusableCellsForRow($tr) {
    var selectors = ['.mapper-core-in-receptor:visible', '.mapper-core-receptor-html-view:visible', '.mapper-tree-inner:visible'];
    if (!mapperTreeIsTextMode()) {
      selectors.push('.mapper-tree-o1:visible', '.mapper-tree-o2:visible', '.mapper-tree-o3:visible', '.mapper-tree-o4:visible');
    }
    return $tr.find(selectors.join(','));
  }

  function mapperTreeFocusNextTableCellFrom($target) {
    var $tr = $target.closest('tr');
    var $cells = mapperTreeFocusableCellsForRow($tr);
    var idx = $cells.index($target);
    if (idx > -1 && idx < $cells.length - 1) {
      $cells.eq(idx + 1).focus().select();
      return;
    }
    var $next = $tr.next('tr');
    if (!$next.length) {
      mapperTreeEnsureTrailingBlankRow();
      $next = $tr.next('tr');
    }
    if (!$next.length) {
      return;
    }
    if (mapperTreeRowBlank($next)) {
      mapperTreeFocusReceptorCell($next);
    } else {
      mapperTreeFocusInnerCell($next);
    }
  }

  function mapperTreeFindBlankRow() {
    var $hit = $();
    $('#mapper-tree-input-tbody tr').each(function () {
      if (mapperTreeRowBlank($(this))) {
        $hit = $(this);
        return false;
      }
    });
    return $hit;
  }

  function mapperTreeSetDemoReceptor($tr, receptor) {
    var raw = String(receptor || '').trim();
    var resolved = null;
    if (raw && window.mapperTreeResolveEntry) {
      resolved = window.mapperTreeResolveEntry(raw);
    }
    if (!resolved && raw && window.MAPPER_CORE_RESOLVE) {
      resolved = window.MAPPER_CORE_RESOLVE[raw.toUpperCase()] || null;
    }
    if (resolved) {
      mapperTreeSetResolved($tr, resolved);
      return;
    }
    $tr.find('.mapper-core-in-receptor').val(raw);
    $tr.find('.mapper-core-receptor-entry').val('');
    $tr.removeData('mapperCoreUnmatchedRaw');
    $tr.removeClass('mapper-core-row-invalid');
    mapperTreeSyncClearBtn($tr);
  }

  function mapperTreeFillDemoRows() {
    var textMode = mapperTreeIsTextMode();
    suppressRedraw = true;
    mapperTreeDestroyAllRows();
    $('#mapper-tree-input-tbody').empty();
    $('#mapper-tree-messages').empty();
    MAPPER_TREE_DEMO_ROWS.forEach(function (row) {
      mapperTreeAppendRow(true);
      var $tr = $('#mapper-tree-input-tbody tr').last();
      mapperTreeSetDemoReceptor($tr, row.receptor);
      if (textMode) {
        $tr.find('.mapper-tree-inner').val(row.text || '');
        $tr.find('.mapper-tree-o1, .mapper-tree-o2, .mapper-tree-o3, .mapper-tree-o4').val('');
        if (row.text && !window.MAPPER_CORE_LABEL_COLORS[row.text]) {
          window.MAPPER_CORE_LABEL_COLORS[row.text] = mapperTreeDefaultHexForLabelKey(row.text);
        }
      } else {
        var values = row.numeric || [];
        $tr.find('.mapper-tree-inner').val(values[0] != null ? values[0] : '');
        $tr.find('.mapper-tree-o1').val(values[1] != null ? values[1] : '');
        $tr.find('.mapper-tree-o2').val(values[2] != null ? values[2] : '');
        $tr.find('.mapper-tree-o3').val(values[3] != null ? values[3] : '');
        $tr.find('.mapper-tree-o4').val(values[4] != null ? values[4] : '');
      }
    });
    suppressRedraw = false;
    mapperTreeEnsureTrailingBlankRow();
    mapperTreeRefreshInnerSwatches();
    mapperTreeSyncRemoveButtons();
    mapperTreeCompactReceptorRowsAfterInput();
    mapperTreeSyncInnerHints();
    $('#mapper-tree-clear-rows').removeClass('mapper-core-clear-clean');
    if (!textMode) {
      mapperTreeApplyDemoColorPresets();
    }
    mapperTreeRedrawNow();
  }

  function mapperTreeApplyDemoColorPresets() {
    if (!$.fn.spectrum) { return; }
    var DEMO = [
      { ring: 'Inner', style: 'Three_RWB' },
      { ring: 'O1',    style: 'Two'       },
      { ring: 'O2',    style: 'Three_TWM' },
      { ring: 'O3',    style: 'Three_RWB' },
      { ring: 'O4',    style: 'Two'       }
    ];
    var CUSTOM = {
      Inner:  ['#a00000', '#1a80bb'],
      Outer1: ['#97a6c4', '#1a2b3c'],
      Outer2: ['#298c8c', '#800074'],
      Outer3: ['#a00000', '#1a80bb'],
      Outer4: ['#ffffff', '#2ca02c']
    };
    // Apply preset (sets styling dict + UI visibility) then override colors
    DEMO.forEach(function (d) {
      try { $('#mapper-tree-cstyle-' + d.ring).val(d.style).trigger('change'); } catch (e) {}
    });
    // Override Tree_colors with the custom demo palette
    if (Tree_colors) {
      Object.keys(CUSTOM).forEach(function (k) { Tree_colors[k] = CUSTOM[k]; });
    }
    // Sync Spectrum pickers to show the custom colors
    DEMO.forEach(function (d) {
      var tKey = { Inner:'Inner', O1:'Outer1', O2:'Outer2', O3:'Outer3', O4:'Outer4' }[d.ring];
      try { $('#mapper-tree-cpicker-min-' + d.ring).spectrum('set', CUSTOM[tKey][0]); } catch (e) {}
      try { $('#mapper-tree-cpicker-max-' + d.ring).spectrum('set', CUSTOM[tKey][1]); } catch (e) {}
    });
  }

  // ── Color pickers (Colors dropdown) ──────────────────────────────────────
  function mapperTreeInitColorPickers() {
    if (!$.fn.spectrum || !$.fn.select2) { return; }

    var RINGS  = ['Inner', 'O1', 'O2', 'O3', 'O4'];
    var KEYS   = { Inner: 'Inner', O1: 'Outer1', O2: 'Outer2', O3: 'Outer3', O4: 'Outer4' };
    var PRESETS = {
      One:       { setup: 'One',   start: '#ffffff', end: '#707070' },
      Two:       { setup: 'Two',   start: '#97a6c4', end: '#384860' },
      Three_RWB: { setup: 'Three', start: '#a00000', end: '#1a80bb' },
      Three_TWM: { setup: 'Three', start: '#298c8c', end: '#800074' }
    };

    var SPECTRUM_PALETTE = [
      ['#000000', '#ffffff', '#a00000', '#1a80bb', '#298c8c', '#800074', '#97a6c4', '#384860'],
      ['#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22']
    ];

    function formatColorScheme(option) {
      if (!option.id) { return option.text; }
      var colorMap = {
        One:       [null, null, '#707070'],
        Two:       ['#97a6c4', null, '#384860'],
        Three_RWB: ['#a00000', '#ffffff', '#1a80bb'],
        Three_TWM: ['#298c8c', '#ffffff', '#800074']
      };
      var colors = colorMap[option.id];
      if (!colors) { return option.text; }
      var chips = colors.map(function (c) {
        return '<span style="display:inline-block;width:12px;height:12px;margin-left:4px;border:1px solid #ccc;border-radius:2px;background:' +
          (c || 'transparent') + ';' + (c ? '' : 'opacity:0') + '"></span>';
      }).join('');
      return $('<span style="display:flex;align-items:center;"><span style="min-width:50px;text-align:center;margin-right:5px;">' +
        option.text + '</span>' + chips + '</span>');
    }

    function applyPreset(ring, presetKey, silent) {
      var p = PRESETS[presetKey] || PRESETS.One;
      var k = KEYS[ring];
      if (Tree_colors)             { Tree_colors[k]             = [p.start, p.end]; }
      if (Tree_circle_styling_dict) { Tree_circle_styling_dict[k] = p.setup; }
      var $minWrap = $('#mapper-tree-cpicker-min-wrap-' + ring);
      var $midWrap = $('#mapper-tree-cpicker-mid-wrap-' + ring);
      if ($minWrap.length) { $minWrap.css('visibility', p.setup !== 'One' ? 'visible' : 'hidden'); }
      if ($midWrap.length) { $midWrap.css('visibility', p.setup === 'Three' ? 'visible' : 'hidden'); }
      try { $('#mapper-tree-cpicker-min-' + ring).spectrum('set', p.start); } catch (e) {}
      try { $('#mapper-tree-cpicker-max-' + ring).spectrum('set', p.end); } catch (e) {}
      if (!silent) { mapperTreeScheduleRedraw(); }
    }

    RINGS.forEach(function (ring) {
      var k = KEYS[ring];
      var initStart = (Tree_colors && Tree_colors[k]) ? Tree_colors[k][0] : '#ffffff';
      var initEnd   = (Tree_colors && Tree_colors[k]) ? Tree_colors[k][1] : '#000000';

      // Wrap min/mid cells in span with ID for easy toggle — the HTML uses class;
      // we add the id attribute here so applyPreset can find them.
      $('#mapper-tree-cpicker-min-' + ring).closest('.mapper-tree-cpicker-min-wrap')
        .attr('id', 'mapper-tree-cpicker-min-wrap-' + ring);
      $('#mapper-tree-cpicker-mid-' + ring).closest('.mapper-tree-cpicker-mid-wrap')
        .attr('id', 'mapper-tree-cpicker-mid-wrap-' + ring);

      // Min picker
      $('#mapper-tree-cpicker-min-' + ring).spectrum({
        color: initStart, preferredFormat: 'hex', showInput: true, showPalette: true,
        palette: SPECTRUM_PALETTE,
        change: function (c) { if (Tree_colors) { Tree_colors[KEYS[ring]][0] = c.toHexString(); } mapperTreeScheduleRedraw(); },
        move:   function (c) { if (Tree_colors) { Tree_colors[KEYS[ring]][0] = c.toHexString(); } mapperTreeScheduleRedraw(); }
      });

      // Mid picker — always disabled (white midpoint for Three-color modes)
      $('#mapper-tree-cpicker-mid-' + ring).spectrum({
        color: '#ffffff', preferredFormat: 'hex', showInput: true, disabled: true
      });

      // Max picker
      $('#mapper-tree-cpicker-max-' + ring).spectrum({
        color: initEnd, preferredFormat: 'hex', showInput: true, showPalette: true,
        palette: SPECTRUM_PALETTE,
        change: function (c) { if (Tree_colors) { Tree_colors[KEYS[ring]][1] = c.toHexString(); } mapperTreeScheduleRedraw(); },
        move:   function (c) { if (Tree_colors) { Tree_colors[KEYS[ring]][1] = c.toHexString(); } mapperTreeScheduleRedraw(); }
      });

      // Style select — Select2 with colour chips
      $('#mapper-tree-cstyle-' + ring).select2({
        templateResult: formatColorScheme,
        templateSelection: formatColorScheme,
        width: 'resolve',
        dropdownParent: $('#mapper-tree-colors-menu')
      }).on('change', function () {
        applyPreset(ring, $(this).val());
      });

      // Apply default "One" preset silently
      applyPreset(ring, 'One', true);
    });
  }

  function mapperTreeSyncColorRings() {
    var KEYS = { Inner: 'Inner', O1: 'Outer1', O2: 'Outer2', O3: 'Outer3', O4: 'Outer4' };
    var anyVisible = false;
    ['Inner', 'O1', 'O2', 'O3', 'O4'].forEach(function (ring) {
      var k = KEYS[ring];
      var hasData = !!(Tree_circles && Object.keys(Tree_circles).some(function (id) {
        return Tree_circles[id] && Object.prototype.hasOwnProperty.call(Tree_circles[id], k);
      }));
      $('#mapper-tree-color-row-' + ring).toggle(hasData);
      if (hasData) { anyVisible = true; }
    });
    $('#mapper-tree-colors-no-data').toggle(!anyVisible);
  }

  function mapperTreeUpdateSliderStates() {
    var $cg = $('#classFontSizeSlider').closest('.form-group');
    $cg.css('opacity', skipClassLevel ? '0.4' : '');
    $('#classFontSizeSlider').prop('disabled', !!skipClassLevel);
  }

  function mapperBoot() {
    $('.mapper-core-booting').removeClass('mapper-core-booting');
    var boot = window.MAPPER_TREE_BOOT || {};
    ORIG_SKEL = boot.tree_skeleton;
    if (typeof ORIG_SKEL === 'string') {
      ORIG_SKEL = JSON.parse(ORIG_SKEL);
    }
    ORIG_OPTS = boot.tree_options;
    if (typeof ORIG_OPTS === 'string') {
      ORIG_OPTS = JSON.parse(ORIG_OPTS);
    }

    ServerReceptorDict =
      typeof boot.receptor_dict === 'string'
        ? JSON.parse(boot.receptor_dict || '{}')
        : boot.receptor_dict || {};

    ServerGeneDict =
      typeof boot.gene_dict === 'string' ? JSON.parse(boot.gene_dict || '{}') : boot.gene_dict || {};

    if (!ORIG_OPTS) {
      ORIG_OPTS = {};
    }
    ORIG_OPTS.anchor = ORIG_OPTS.anchor || 'tree_plot';

    mapperTreeBuildStemLabelDicts();
    mapperTreeLeafLabelLookupBuild();

    // Extend autocomplete and resolve map with non-human-only receptors.
    // One autocomplete entry per stem (first species entry as id); all species
    // entries share the same "(no human ortholog)" display across all label modes.
    if (window.MAPPER_TREE_SPECIES_DATA && window.MAPPER_TREE_SPECIES_DATA.nonhuman_only) {
      var addedNhoStems = {};
      window.MAPPER_TREE_SPECIES_DATA.nonhuman_only.forEach(function (nho) {
        var stemU = nho.stem.toUpperCase();
        var TAG = '(no human ortholog)';
        var dispText = stemU + ' ' + TAG;
        var nhoHtml = '<span>' + $('<span>').text(stemU).html() + ' <em>' + TAG + '</em></span>';

        // Register MAPPER_CORE_ENTRY_META for every individual species entry so that
        // whichever species is later picked in the dropdown, the cell always shows
        // "STEM (no human ortholog)" in all Receptor Names modes.
        if (window.MAPPER_CORE_ENTRY_META && !window.MAPPER_CORE_ENTRY_META[nho.entry]) {
          window.MAPPER_CORE_ENTRY_META[nho.entry] = {
            name_html:  nhoHtml,
            name_plain: dispText,
            uniprot:    dispText   // UniProt mode also shows "(no human ortholog)"
          };
        }
        if (window.MAPPER_CORE_RESOLVE) {
          window.MAPPER_CORE_RESOLVE[nho.entry.toUpperCase()] = nho.entry;
        }

        // Only one autocomplete item per stem (keyed on the first species entry)
        if (addedNhoStems[stemU]) { return; }
        addedNhoStems[stemU] = true;

        var alreadyIn = (window.receptorSelect2Data || []).some(function (x) { return x.id === nho.entry; });
        if (!alreadyIn) {
          window.receptorSelect2Data = window.receptorSelect2Data || [];
          window.receptorSelect2Data.push({
            id: nho.entry,
            text: dispText,
            name_plain: dispText,
            name_html: nhoHtml
          });
          if (window.MAPPER_CORE_RESOLVE) {
            window.MAPPER_CORE_RESOLVE[stemU] = nho.entry;
          }
        }
      });
    }

    Tree_circles = {};
    Tree_colors = {
      Inner: ['#FFFFFF', '#000000'],
      Outer1: ['#FFFFFF', '#0000FF'],
      Outer2: ['#FFFFFF', '#FF0000'],
      Outer3: ['#FFFFFF', '#22c7b1'],
      Outer4: ['#FFFFFF', '#008000'],
      Outer5: ['#FFFFFF', '#FFA500']
    };

    Tree_circle_styling_dict = {
      Inner: 'One',
      Outer1: 'One',
      Outer2: 'One',
      Outer3: 'One',
      Outer4: 'One',
      Outer5: 'One'
    };

    Label_dict = {
      Inner: 'Dot 1',
      Outer1: 'Dot 2',
      Outer2: 'Dot 3',
      Outer3: 'Dot 4',
      Outer4: 'Dot 5',
      Outer5: 'Dot 5'
    };

    window.mapperTreeNumericDatatypes = {
      Inner: 'Continuous',
      Outer1: 'Continuous',
      Outer2: 'Continuous',
      Outer3: 'Continuous',
      Outer4: 'Continuous',
      Outer5: 'Continuous'
    };
    window.mapperTreeDiscreteDatatypes = {
      Inner: 'Discrete',
      Outer1: 'Discrete',
      Outer2: 'Discrete',
      Outer3: 'Discrete',
      Outer4: 'Discrete',
      Outer5: 'Discrete'
    };
    Tree_datatypes_dict = $.extend({}, window.mapperTreeNumericDatatypes);

    Tree_textlegend_styling = {
      layoutMode: 'row',
      columns: 2,
      sortDirection: 'Vertically',
      TreeLegendPosition: 'Top',
      Fontsize: '14px'
    };

    ShowLegend = true;
    TreeLegendPosition = 'Top';

    styling_circles = {
      starter: 1,
      clean: true,
      gradient: true,
      circle_spacer: 22,
      circle_size: 12,
      mode: 'Numeric',
      skipNumericViewBoxAdjust: true
    };

    maxLeafNodeLength_scaler = 10;

    window.TREE_UI = window.TREE_UI || { layout: 'Tree - Circular', leafLabelType: 'Protein' };
    $('#mapper-tree-input-tbody').empty();
    mapperTreeAppendRow();
    mapperTreeSetInputMode('numeric');
    mapperTreeSyncFirstRowPlaceholder(); // mode-switch is a no-op on boot; set placeholder explicitly
    mapperTreeSyncLeafLabelUi('IUPHAR');
    mapperTreeUpdateSortHeaders();
    mapperTreeInitColorPickers();

    // Wire gradient bar label inputs to Label_dict
    var LABEL_INPUT_MAP = [
      { ring: 'Inner', key: 'Inner'  },
      { ring: 'O1',    key: 'Outer1' },
      { ring: 'O2',    key: 'Outer2' },
      { ring: 'O3',    key: 'Outer3' },
      { ring: 'O4',    key: 'Outer4' }
    ];
    LABEL_INPUT_MAP.forEach(function (m) {
      var $inp = $('#mapper-tree-lbl-' + m.ring);
      if (!$inp.length) { return; }
      $inp.val(Label_dict && Label_dict[m.key] ? Label_dict[m.key] : m.key);
      $inp.off('input.mapperTreeLabel').on('input.mapperTreeLabel', function () {
        if (Label_dict) {
          Label_dict[m.key] = $(this).val() || m.key;
        }
        mapperTreeScheduleRedraw();
      });
    });

    /** Legend toggle */
    $('#legendToggleBtn')
      .off('click.mapperTree')
      .on('click.mapperTree', function () {
        ShowLegend = !ShowLegend;
        $(this).text(ShowLegend ? 'Shown' : 'Hidden');
        $(this).toggleClass('btn-success', ShowLegend);
        $(this).toggleClass('btn-danger', !ShowLegend);
        mapperTreeRedrawNow();
      });

    $('#toggleLegendPosition')
      .off('click.mapperTree')
      .on('click.mapperTree', function () {
        var cur = TreeLegendPosition === 'Bottom' ? 'Bottom' : 'Top';
        var next = cur === 'Top' ? 'Bottom' : 'Top';
        TreeLegendPosition = next;
        if (Tree_textlegend_styling) {
          Tree_textlegend_styling.TreeLegendPosition = next;
        }
        $(this).text(next);
        mapperTreeRedrawNow();
      });
    $('#toggleLegendPosition').text(TreeLegendPosition || 'Top');

    skipClassLevel = false;
    skipLigandtypeLevel = false;
    classLevelAutoOff = false;
    mapperTreeUpdateSliderStates();

    $('#mapper-tree-circle-size-slider')
      .off('input.mapperTree')
      .on('input.mapperTree', function () {
        $('#mapper-tree-circle-size-val').text(String($(this).val()));
        mapperTreeScheduleRedraw();
      });
    $('#mapper-tree-circle-spacer-slider')
      .off('input.mapperTree')
      .on('input.mapperTree', function () {
        $('#mapper-tree-circle-spacer-val').text(String($(this).val()));
        mapperTreeScheduleRedraw();
      });

    $('#classFontSizeSlider, #ligandTypeFontSizeSlider, #receptorFamilyFontSizeSlider, #receptorFontSizeSlider')
      .off('input.mapperTree')
      .on('input.mapperTree', function () {
        var idBase = $(this).attr('id').replace('Slider', 'Value');
        $('#' + idBase).text($(this).val());
        mapperTreeScheduleRedraw();
      });

    $(document).off(
      'input.mapperTree blur.mapperTree change.mapperTree',
      '#mapper-tree-input-tbody input, #mapper-tree-input-tbody textarea'
    ).on(
      'input.mapperTree blur.mapperTree change.mapperTree',
      '#mapper-tree-input-tbody input, #mapper-tree-input-tbody textarea',
      function () {
        $('#mapper-tree-clear-rows').removeClass('mapper-core-clear-clean');
        mapperTreeEnsureTrailingBlankRow();
        _treeCosmeticSyncDebounced.schedule();
        mapperTreeScheduleRedraw();
      }
    );

    function mapperTreeDoWholeTableClear() {
      mapperTreeDestroyAllRows();
      $('#mapper-tree-input-tbody').empty();
      $('#mapper-tree-messages').empty();
      MAPPER_TREE_MODE_SNAPSHOTS.numeric = null;
      MAPPER_TREE_MODE_SNAPSHOTS.categorical = null;
      mapperTreeAppendRow();
      mapperTreeCompactReceptorRowsAfterInput();
      mapperTreeSyncFirstRowPlaceholder();
      mapperTreeRedrawNow();
      $('#mapper-tree-clear-rows').addClass('mapper-core-clear-clean').blur();
    }
    function mapperTreeDoCurrentModeClear() {
      var curKey = mapperTreeIsTextMode() ? 'categorical' : 'numeric';
      mapperTreeDestroyAllRows();
      $('#mapper-tree-input-tbody').empty();
      $('#mapper-tree-messages').empty();
      MAPPER_TREE_MODE_SNAPSHOTS[curKey] = null;
      mapperTreeAppendRow();
      mapperTreeCompactReceptorRowsAfterInput();
      mapperTreeSyncFirstRowPlaceholder();
      mapperTreeRedrawNow();
      $('#mapper-tree-clear-rows').addClass('mapper-core-clear-clean').blur();
    }
    $('#mapper-tree-clear-whole').off('click.mapperTree').on('click.mapperTree', function(e) {
      e.preventDefault(); mapperTreeDoWholeTableClear();
    });
    $('#mapper-tree-clear-mode').off('click.mapperTree').on('click.mapperTree', function(e) {
      e.preventDefault(); mapperTreeDoCurrentModeClear();
    });
    $('.mapper-core-clear-menu').on('click.mapperTree', '[data-tree-col]', function(e) {
      e.preventDefault();
      var col = $(this).data('tree-col');
      $('#mapper-tree-input-tbody tr').each(function() {
        $(this).find('.mapper-tree-' + col).val('').trigger('input');
      });
      $('#mapper-tree-clear-rows').removeClass('mapper-core-clear-clean');
      mapperTreeScheduleRedraw();
    });

    $('#mapper-tree-mode-numeric-btn')
      .off('click.mapperTree')
      .on('click.mapperTree', function () {
        mapperTreeSetInputMode('numeric');
      });
    $('#mapper-tree-mode-labels-btn')
      .off('click.mapperTree')
      .on('click.mapperTree', function () {
        mapperTreeSetInputMode('text');
      });

    $('#mapper-tree-demo-rows')
      .off('click.mapperTree')
      .on('click.mapperTree', function (eDemo) {
        eDemo.preventDefault();
        mapperTreeFillDemoRows();
      });

    $('#mapper-tree-species-toggle')
      .off('click.mapperTreeSpecies')
      .on('click.mapperTreeSpecies', function () {
        mapperTreeSpeciesActive = !mapperTreeSpeciesActive;
        $(this).toggleClass('btn-primary', mapperTreeSpeciesActive).toggleClass('btn-default', !mapperTreeSpeciesActive);
        var speciesCellDisp = mapperTreeSpeciesActive ? 'table-cell' : 'none';
        $('.mapper-tree-species-cell').css('display', speciesCellDisp);
        $('.mapper-tree-species-th').css('display', speciesCellDisp);
        $('#mapper-tree-species-names-wrap').toggle(mapperTreeSpeciesActive);
        mapperTreeApplyLeftPanelWidth();
        if (mapperTreeSpeciesActive) {
          mapperTreeInitAllSpeciesDropdowns();
        }
        mapperTreeScheduleRedraw();
      });

    $('.mapper-tree-species-name-btn')
      .off('click.mapperTreeSpeciesName')
      .on('click.mapperTreeSpeciesName', function (eSpec) {
        eSpec.preventDefault();
        mapperTreeSpeciesNameFormat = $(this).attr('data-value') || 'common';
        $('.mapper-tree-species-name-btn').each(function () {
          var isActive = $(this).attr('data-value') === mapperTreeSpeciesNameFormat;
          $(this).toggleClass('btn-primary', isActive).toggleClass('btn-outline-primary', !isActive);
          if (isActive) { $(this).addClass('active'); } else { $(this).removeClass('active'); }
        });
        mapperTreeInitAllSpeciesDropdowns();
        mapperTreeScheduleRedraw();
      });

    $(document)
      .off('change.mapperTreeSpeciesSel', '#mapper-tree-input-tbody .mapper-tree-species-select')
      .on('change.mapperTreeSpeciesSel', '#mapper-tree-input-tbody .mapper-tree-species-select', function () {
        mapperTreeScheduleRedraw();
      });

    $('#mapper-tree-input-table')
      .off('click.mapperTreeRemove', '.mapper-core-remove-row')
      .on('click.mapperTreeRemove', '.mapper-core-remove-row', function () {
        var $trRemove = $(this).closest('tr');
        if (mapperTreeRowBlank($trRemove) && $trRemove.is(':last-child')) {
          return;
        }
        mapperTreeDestroyRowAc($trRemove);
        $trRemove.remove();
        $('#mapper-tree-clear-rows').removeClass('mapper-core-clear-clean');
        mapperTreeEnsureTrailingBlankRow();
        _treeCosmeticSyncDebounced.schedule();
        mapperTreeScheduleRedraw();
      });

    $('#mapper-tree-input-table')
      .off('click.mapperTreeSort', 'th.mapper-core-sortable-head')
      .on('click.mapperTreeSort', 'th.mapper-core-sortable-head', function () {
        mapperTreeSortSerializedRows($(this).attr('data-mapper-tree-sort-col'));
      });

    $('#mapper-tree-input-table')
      .off('keydown.mapperTreeSort', 'th.mapper-core-sortable-head')
      .on('keydown.mapperTreeSort', 'th.mapper-core-sortable-head', function (eSort) {
        if (eSort.which === 13 || eSort.which === 32) {
          eSort.preventDefault();
          mapperTreeSortSerializedRows($(this).attr('data-mapper-tree-sort-col'));
        }
      });

    $('#mapper-tree-input-table')
      .off('keydown.mapperTreeTab', '.mapper-core-in-receptor, .mapper-core-receptor-html-view, .mapper-tree-inner, .mapper-tree-o1, .mapper-tree-o2, .mapper-tree-o3, .mapper-tree-o4')
      .on('keydown.mapperTreeTab', '.mapper-core-in-receptor, .mapper-core-receptor-html-view, .mapper-tree-inner, .mapper-tree-o1, .mapper-tree-o2, .mapper-tree-o3, .mapper-tree-o4', function (eTab) {
        if (eTab.which !== 9 || eTab.shiftKey) {
          return;
        }
        if ($(this).hasClass('mapper-core-in-receptor') && $(this).data('ui-autocomplete')) {
          var $menu = $(this).autocomplete('widget');
          if ($menu && $menu.is(':visible') && $menu.find('.ui-state-focus, .ui-state-active').length) {
            return;
          }
        }
        eTab.preventDefault();
        mapperTreeFocusNextTableCellFrom($(this));
      });

    $('.mapper-tree-leaf-btn').on('click.mapperTreeLeaf', function (eX) {
      eX.preventDefault();
      var v = ($(this).attr('data-value') || 'IUPHAR').trim();
      mapperTreeSyncLeafLabelUi(v);
      mapperTreeScheduleRedraw();
    });

    $('#mapper-tree-input-table').on(
      'paste.mapperTreePaste',
      function (ePz) {
        var evPz = ePz.originalEvent || ePz;
        var textPz = evPz.clipboardData ? evPz.clipboardData.getData('text/plain') : '';
        if (!textPz || textPz.indexOf('\t') === -1) {
          return;
        }
        if (typeof textPz === 'undefined') {
          return;
        }
        ePz.preventDefault();
        $('#mapper-tree-clear-rows').removeClass('mapper-core-clear-clean');
        suppressRedraw = true;
        textPz.split(/\r?\n/).forEach(function (ln) {
          if (!ln.trim()) {
            return;
          }
          var $trPZ = mapperTreeFindBlankRow();
          if (!$trPZ.length) {
            mapperTreeAppendRow(true);
            $trPZ = $('#mapper-tree-input-tbody tr').last();
          }
          var ps = mapperTreePasteSplit(ln);
          mapperTreeDestroyAc($trPZ.find('.mapper-core-in-receptor'));
          $trPZ.find('.mapper-core-receptor-input-wrap .mapper-core-in-receptor').remove();
          var $tx = $(
            '<textarea class="form-control input-sm mapper-core-in-receptor" rows="1" autocomplete="off" spellcheck="false"></textarea>'
          ).val(ps.r || '');
          $trPZ.find('.mapper-core-receptor-input-wrap').prepend($tx);
          mapperTreeBindAc($tx);

          var restPZ = [];
          if (Array.isArray(ps.rest)) {
            restPZ = ps.rest;
          } else if (ps.rest != null && String(ps.rest).trim()) {
            restPZ = String(ps.rest)
              .split(';')
              .map(function (s) {
                return s.trim();
              });
          }
          if (restPZ[0]) {
            $trPZ.find('.mapper-tree-inner').val(restPZ[0]);
          }
          if (restPZ[1]) {
            $trPZ.find('.mapper-tree-o1').val(restPZ[1]);
          }
          if (restPZ[2]) {
            $trPZ.find('.mapper-tree-o2').val(restPZ[2]);
          }
          if (restPZ[3]) {
            $trPZ.find('.mapper-tree-o3').val(restPZ[3]);
          }
          if (restPZ[4]) {
            $trPZ.find('.mapper-tree-o4').val(restPZ[4]);
          }

          var upU = ps.r ? ps.r.trim().toUpperCase() : '';
          if (upU && window.MAPPER_CORE_RESOLVE && window.MAPPER_CORE_RESOLVE[upU]) {
            mapperTreeSetResolved($trPZ, window.MAPPER_CORE_RESOLVE[upU]);
          }
        });
        suppressRedraw = false;
        mapperTreeEnsureTrailingBlankRow();
        mapperTreeSyncRemoveButtons();
        mapperTreeCompactReceptorRowsAfterInput();
        mapperTreeScheduleRedraw();
      }
    );

    if (typeof window.mapperCoreInitGpcromePickerModal === 'function') {
      window.mapperCoreInitGpcromePickerModal({
        pickerRows: $.isArray(window.MAPPER_CORE_GPCROME_PICKER_ROWS) ? window.MAPPER_CORE_GPCROME_PICKER_ROWS : [],
        maxRows: MAPPER_TREE_MAX_ROWS,
        onAdd: function (entryIds, meta) {
          suppressRedraw = true;
          var isNumberMode = !mapperTreeIsTextMode();
          var numberMap = (meta && meta.groupNameById && isNumberMode)
            ? window.mapperCoreBuildSequentialNumberMap(meta.groupNameById)
            : null;
          for (var ix = 0; ix < (entryIds || []).length; ix++) {
            var idPz = entryIds[ix];
            if ($('#mapper-tree-input-tbody tr').length >= MAPPER_TREE_MAX_ROWS && !mapperTreeFindBlankRow().length) {
              break;
            }
            var $rowPz = mapperTreeFindBlankRow();
            if (!$rowPz.length) {
              mapperTreeAppendRow(true);
              $rowPz = $('#mapper-tree-input-tbody tr').last();
            }
            mapperTreeSetResolved($rowPz, idPz);
            if (meta && meta.groupNameById && meta.groupNameById[idPz] != null) {
              var assignedValue = isNumberMode ? numberMap[idPz] : meta.groupNameById[idPz];
              if (assignedValue != null) {
                $rowPz.find('.mapper-tree-inner').val(String(assignedValue)).trigger('input');
              }
            }
          }
          suppressRedraw = false;
          $('#mapper-tree-clear-rows').removeClass('mapper-core-clear-clean');
          mapperTreeEnsureTrailingBlankRow();
          mapperTreeCompactReceptorRowsAfterInput();
          mapperTreeScheduleRedraw();
        }
      });
    }

    mapperTreeRedrawNow();

    /** Resume editing receptor from resolved HTML capsule (delegated once). */
    $(document)
      .off('click.mapperTreeHtmlEdit', '#mapper-tree-input-tbody .mapper-core-receptor-html-view')
      .on(
        'click.mapperTreeHtmlEdit',
        '#mapper-tree-input-tbody .mapper-core-receptor-html-view',
        function () {
          var $trEv = $(this).closest('tr');
          mapperTreeDestroyAc($trEv.find('.mapper-core-in-receptor'));
          var hidEv = ($trEv.find('.mapper-core-receptor-entry').val() || '').trim();
          $(this).hide().empty();
          var $inp2 = $(
            '<textarea class="form-control input-sm mapper-core-in-receptor" rows="1" autocomplete="off" spellcheck="false"></textarea>'
          );
          var seedEv = mapperTreeEditSeed(hidEv);
          $inp2.val(seedEv);
          if (!seedEv) { $trEv.find('.mapper-core-receptor-input-wrap').addClass('is-empty'); }
          $trEv.find('.mapper-core-receptor-input-wrap .mapper-core-in-receptor').remove();
          $trEv.find('.mapper-core-receptor-input-wrap').prepend($inp2);
          $trEv.find('.mapper-core-receptor-entry').val('');
          mapperTreeBindAc($inp2);
          $inp2.show().focus();
          window.setTimeout(function () {
            var t = ($inp2.val() || '').trim();
            if (t.length >= 1 && $inp2.data('ui-autocomplete')) { $inp2.autocomplete('search', t); }
          }, 0);
          mapperTreeSyncClearBtn($trEv);
          mapperTreeSyncRemoveButtons();
          mapperTreeCompactReceptorRowsAfterInput();
        }
      );

    mapperTreeSyncClearDropdown();
  }

  if (typeof window.mapperTreeResolveEntry !== 'function' && window.MAPPER_CORE_RESOLVE) {
    window.mapperTreeResolveEntry = function (raw) {
      if (!raw || !String(raw).trim()) {
        return null;
      }
      return window.MAPPER_CORE_RESOLVE[String(raw).trim().toUpperCase()] || null;
    };
  }

  $(function () {
    if (!window.MAPPER_TREE_BOOT) {
      window.MAPPER_TREE_BOOT = {};
    }
    mapperBoot();
  });
})(jQuery);

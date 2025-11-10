/* ==========================================================================
   Structure Browser – Modern JS
   - All inline JS moved here
   - Grouped by sections with loud headers
   Dependencies: jQuery, DataTables + Buttons, BusyLoad, Select2,
                 NorgesDTFilterBuilder, XLSX, Bootstrap (popover)
   ========================================================================== */

(function ($, window, document) {
  "use strict";

  // ===== 0) GLOBALS & SHARED STATE ==========================================
  // Selection store exposed for other scripts that might rely on it
  window.globalSelectedIds = window.globalSelectedIds || new Set();

  // Track "representatives only" toggle
  let repFilterOn = false;

  // ===== 1) RENDERING HELPERS (cell render/format) ===========================
  function expand1LineRender(d, type) {
    const decode = (s) => {
      const div = document.createElement('div');
      div.innerHTML = s || '';
      return div.textContent || div.innerText || '';
    };

    const toItems = (val) => {
      if (val == null) return [];
      if (Array.isArray(val)) {
        return val
          .map(x => typeof x === 'string' ? x : (x?.name ?? x?.text ?? x?.label ?? ''))
          .map(s => decode(s).trim())
          .filter(Boolean);
      }
      const s = String(val).trim();
      if (!s || s === '-') return [];
      const dec = decode(s.replace(/<br\s*\/?>/gi, '\n'));
      return dec.split(/\n+/).map(t => t.trim()).filter(Boolean);
    };

    const items = toItems(d);

    if (type === 'filter' || type === 'sort') return items.join('|');
    if (type && type !== 'display') return items.join('\n');

    if (!items.length) return '-';
    const first = items[0];
    const more  = items.length - 1;
    const hint  = more > 0 ? ` <span class="muted">(+${more} more)</span>` : '';
    const expanded = items.join('<br>');

    return (
      '<span class="expand_1line">' +
        '<span class="collapsed">' + first + hint + '</span>' +
        '<span class="expanded">'  + expanded + '</span>' +
      '</span>'
    );
  }

  // Show first item, reveal full list of links on hover
  function expandFirstThenHoverLinkList(d, type, hrefBuilder, { showCountHint = true } = {}) {
    const arr = Array.isArray(d) ? d.filter(x => x && (x.name || x.id != null)) : [];
    if (!arr.length) return (type === 'display' ? '-' : '');

    const decode = (s) => {
      const div = document.createElement('div');
      div.innerHTML = s || '';
      return div.textContent || div.innerText || '';
    };

    if (type && type !== 'display') {
      return arr.map(x => decode(x.name || '')).join('\n');
    }

    const first = arr[0];
    const firstHref = (first.id != null && first.id !== '-' && first.id !== '')
      ? (hrefBuilder ? hrefBuilder(first.id, first) : `https://gpcrdb.org/ligand/${first.id}/info`)
      : null;

    const collapsed = firstHref
      ? `<a href="${firstHref}" target="_blank" rel="noopener">${first.name || ''}</a>`
      : (first.name || '');

    const more = arr.length - 1;
    const hint = (showCountHint && more > 0) ? ` <span class="muted">(+${more} more)</span>` : '';

    const expanded = arr.map(x => {
      const nm = x.name || '';
      if (x.id == null || x.id === '' || x.id === '-') return nm;
      const href = hrefBuilder ? hrefBuilder(x.id, x) : `https://gpcrdb.org/ligand/${x.id}/info`;
      return `<a href="${href}" target="_blank" rel="noopener">${nm}</a>`;
    }).join('<br>');

    return (
      '<span class="expand_1line">' +
        '<span class="collapsed">' + collapsed + hint + '</span>' +
        '<span class="expanded">'  + expanded  + '</span>' +
      '</span>'
    );
  }

  // ===== 2) DATATABLES – COLUMNS & DEFS ======================================
  const columnDefs = [
    {
      targets: "_all",
      className: "dt-head-nowrap dt-head-center dt-center dt-body-center",
      createdCell: function (td) {
        td.classList.add("expand_1line");
      }
    },
    { targets: [0], orderable: false }, // checkbox column
    {
      targets: [0, 1, 2, 3, 4],
      createdCell: function (td) { td.classList.add("sticky-col"); },
      className: "sticky-col"
    }
  ];

  const columns = [
    // RECEPTOR
    ({ data: "id", name: "Select", orderable: false, searchable: false, render: (data) => `<input type="checkbox" class="select-row" value="${data}">` }),
    ({ data: "gpcrdb_link", name: "GPCRdb", orderable: false, searchable: false, render: (d,t) => (t === 'display' && d ? `<a href="${d}" target="_blank" rel="noopener"><img class="gpcrdb-link" src="/static/home/logo/gpcr/main.png" width="12" height="12" alt="GPCRdb"></a>` : '') }),
    { data: "Gene", name: "Gene" },
    ({ data: "entry_name", name: "UniProt", render: (d,t,r) => t!=='display' ? d : (r && r.uniprot_link ? `<a href="${r.uniprot_link}" target="_blank">${d ? d.split("_")[0].toUpperCase() : "-"}</a>` : (d ? d.split("_")[0].toUpperCase() : "-")) }),
    ({ data: "iuphar_name", name: "Protein", render: (d,t,r) => { if (t !== 'display') return d; const url = r && r.iuphar_link; return (url && url !== "-") ? `<a href="${url}" target="_blank" rel="noopener">${d}</a>` : d; } }),
    { data: "id", name: "ID", visible: false },

    // CLASSIFICATION
    { data: "family", name: "Receptor family" },
    { data: "class", name: "Class" },
    { data: "species", name: "Species" },

    // STRUCTURE
    { data: "method", name: "Method" },
    ({ data: "pdb", name: "PDB", render: (d,t) => t!=='display' ? d : (d && d!=="-" ? `<a href="https://gpcrdb.org/structure/${d}" target="_blank" rel="noopener">${d}</a>` : d) }),
    ({ data: "refined", name: "Refined structure", render: (d,t) => d === "-" ? "-" : (t === 'display' ? `<a href="${d}" target="_blank">${d.replace("refined/", "").toUpperCase()}_refined</a>` : d.replace("refined/", "").toUpperCase() + "_refined") }),
    ({ data: "resolution", name: "Resolution", render: d => d ? parseFloat(d).toFixed(1) : "-" }),
    { data: "preferred_chain", name: "Preferred chain" },
    { data: "state", name: "State" },
    ({ data: "active_pct", name: "Degree active (%)", render: d => d != null ? Math.round(d) : "-" }),
    { data: "coverage", name: "% of Seq" },

    // SIGNAL PROTEIN
    { data: "arrestin_family", name: "Family" },
    ({ data: "arrestin_name", name: "Subtype", render: (d,t,r) => t!=='display' ? d : (r?.arrestin_entry && r.arrestin_entry!=="-" ? `<a href="https://gpcrdb.org/signprot/${r.arrestin_entry}/" target="_blank" rel="noopener">${d}</a>` : d) }),
    { data: "arrestin_note", name: "Note" },
    { data: "arrestin_coverage", name: "% of Seq" },

    // AUXILIARY PROTEIN
    ({ data: "fusions", name: "Fusion", render: (d,t) => expand1LineRender(d,t) }),
    ({ data: "antibodies", name: "Antibodies", render: (d,t) => expand1LineRender(d,t) }),

    // STRUCTURE LIGAND
    ({ data: "ligands", name: "Name", render: (d,t) => expandFirstThenHoverLinkList(d,t,id => `https://gpcrdb.org/ligand/${id}/info`, { showCountHint: true }) }),
    ({ data: "ligand_type", name: "Type", render: (d,t) => expand1LineRender(d,t) }),
    ({ data: "ligand_role", name: "Modality", render: (d,t) => expand1LineRender(d,t) }),

    // PHYSIOLOGICAL LIGAND
    ({ data: "endo_ligands", name: "Name", render: (d,t) => expandFirstThenHoverLinkList(d,t,id => `https://gpcrdb.org/ligand/${id}/info`, { showCountHint: true }) }),
    ({ data: "endo_type", name: "Type", render: (d,t) => expand1LineRender(d,t) }),

    // SODIUM ION SITE
    { data: "sodium_site", name: "D2x50 S3x39" },
    { data: "sodium", name: "Sodium in structure" },

    // REFERENCE
    { data: "authors", name: "Authors" },
    { data: "reference", name: "Reference" },
    { data: "pub_date", name: "Publication date" }
  ];

  // ===== 3) COLVIS (Table Control) ===========================================
  function controlTableButtonWithGroups(table) {
    const managedIdx = table.columns().indexes().toArray().filter(i => i >= 6); // lock 0..5
    return {
      extend: "colvis",
      text: "Table Control",
      className: "btn-table-control cluster-a",
      columns: managedIdx,
      collectionLayout: "fixed",
      collectionTitle: "Column Visibility"
    };
  }

  function buildColvisMenu($menu, table) {
    $menu.empty();

    const G = [
      { label: 'CLASSIFICATION', items: [ { i:6,t:'Receptor family' }, { i:7,t:'Class' }, { i:8,t:'Species' } ] },
      { label: 'STRUCTURE',      items: [ { i:9,t:'Method' },{ i:10,t:'PDB' },{ i:11,t:'Refined structure' },{ i:12,t:'Resolution' },{ i:13,t:'Preferred chain' },{ i:14,t:'State' },{ i:15,t:'Degree active (%)' },{ i:16,t:'% of Seq' } ] },
      { label: 'SIGNAL PROTEIN', items: [ { i:17,t:'Family' },{ i:18,t:'Subtype' },{ i:19,t:'Note' },{ i:20,t:'% of Seq' } ] },
      { label: 'AUXILIARY PROTEIN', items: [ { i:21,t:'Fusion' },{ i:22,t:'Antibodies' } ] },
      { label: 'STRUCTURE LIGAND', items: [ { i:23,t:'Name' },{ i:24,t:'Type' },{ i:25,t:'Modality' } ] },
      { label: 'PHYSIOLOGICAL LIGAND', items: [ { i:26,t:'Name' },{ i:27,t:'Type' } ] },
      { label: 'SODIUM ION SITE', items: [ { i:28,t:'D2x50 S3x39' },{ i:29,t:'Sodium in structure' } ] },
      { label: 'REFERENCE', items: [ { i:30,t:'Authors' },{ i:31,t:'Reference' },{ i:32,t:'Publication date' } ] }
    ];

    const headings = [];

    const makeBtn = (idx, label) => {
      const $b = $(`<button type="button" class="dt-button buttons-columnVisibility" data-dt-column="${idx}"><span>${label}</span></button>`);
      const sync = () => $b.toggleClass('dt-button-active', table.column(idx).visible());
      $b.on('click', function (e) {
        e.preventDefault(); e.stopPropagation();
        const v = table.column(idx).visible();
        table.column(idx).visible(!v, false);
        table.columns.adjust().draw(false);
        sync();
        syncHeadings();
      });
      sync();
      return $b;
    };

    const syncHeadings = () => {
      headings.forEach(h => {
        const allVisible = h.idxs.every(i => table.column(i).visible());
        h.$el.toggleClass('dt-button-active', allVisible);
      });
    };

    G.forEach(group => {
      const $col = $('<div class="cv-col"></div>');
      const idxs = group.items.map(x => x.i);
      const $heading = $(`<button type="button" class="dt-button dt-colvis-heading"><span>${group.label}</span></button>`).data('idxs', idxs);
      headings.push({ $el: $heading, idxs });
      $col.append($heading);
      group.items.forEach(x => $col.append(makeBtn(x.i, x.t)));
      $menu.append($col);
    });

    $menu.off('click.cvGroupToggle').on('click.cvGroupToggle', '.dt-colvis-heading', function (e) {
      e.preventDefault(); e.stopPropagation();
      const idxs = $(this).data('idxs') || [];
      const show = idxs.some(i => !table.column(i).visible());
      idxs.forEach(i => table.column(i).visible(show, false));
      table.columns.adjust().draw(false);
      syncHeadings();
      $menu.find(`[data-dt-column]`).each(function(){
        const idx = $(this).data('dt-column');
        $(this).toggleClass('dt-button-active', table.column(idx).visible());
      });
    });

    // activate our grid layout
    $menu.addClass('cv-columns');

    // equalize heading heights
    const equalize = () => {
      const $heads = $menu.find('.dt-colvis-heading').css('height', 'auto');
      const maxH = Math.max(0, ...$heads.map((i, el) => $(el).outerHeight()).get());
      $heads.css('height', maxH + 'px');
    };
    equalize();
    $(window).off('resize.cvColvis').on('resize.cvColvis', equalize);

    // Footer with "Show all" + "Hide all"
    // keep these columns hidden regardless of "Show all"
    const ALWAYS_HIDDEN = [5]; // ID column

    // Footer with "Show all" + "Hide all"
    $menu.siblings('.cv-footer').remove();
    const MANAGED = [].concat(...G.map(g => g.items.map(x => x.i)));
    const $footer  = $('<div class="cv-footer"></div>');
    const $showAll = $('<button type="button" class="dt-button"><span>Show all columns</span></button>');
    const $hideAll = $('<button type="button" class="dt-button"><span>Hide all columns</span></button>');

    const refreshStates = () => {
      syncHeadings();
      $menu.find('[data-dt-column]').each(function () {
        const idx = $(this).data('dt-column');
        $(this).toggleClass('dt-button-active', table.column(idx).visible());
      });
      applyGroupBorders(table);
    };

    $showAll.on('click', function (e) {
      e.preventDefault(); e.stopPropagation();
      MANAGED.forEach(i => table.column(i).visible(true,  false));
      // re-hide always-hidden columns
      ALWAYS_HIDDEN.forEach(i => table.column(i).visible(false, false));
      table.columns.adjust().draw(false);
      refreshStates();
    });

    $hideAll.on('click', function (e) {
      e.preventDefault(); e.stopPropagation();
      MANAGED.forEach(i => table.column(i).visible(false, false));
      // also ensure always-hidden stay hidden (no-op but explicit)
      ALWAYS_HIDDEN.forEach(i => table.column(i).visible(false, false));
      table.columns.adjust().draw(false);
      refreshStates();
    });

    $footer.append($hideAll, $showAll);
    $menu.after($footer);


    syncHeadings();
  }

  // ===== 4) GROUP BORDERS (left divider at group starts) ======================
  const BORDER_GROUPS = {
    CLASSIFICATION:         [6,7,8],
    STRUCTURE:              [9,10,11,12,13,14,15,16],
    "SIGNAL PROTEIN":       [17,18,19,20],
    "AUXILIARY PROTEIN":    [21,22],
    "STRUCTURE LIGAND":     [23,24,25],
    "PHYSIOLOGICAL LIGAND": [26,27],
    "SODIUM ION SITE":      [28,29],
    REFERENCE:              [30,31,32]
  };

  function applyGroupBorders(table) {
    const $theadDisp = $(table.table().header());
    const $rowsDisp  = $theadDisp.find('tr');
    const $labelDisp = $rowsDisp.eq(-2); // titles
    const $boundDisp = $rowsDisp.eq(-1); // filters

    const $theadOrig = $('#StructureBrowserTable thead');
    const $rowsOrig  = $theadOrig.find('tr');
    const $labelOrig = $rowsOrig.eq(-2);
    const $boundOrig = $rowsOrig.eq(-1);

    [$labelDisp, $boundDisp, $labelOrig, $boundOrig].forEach($r => $r.children().removeClass('border-left'));
    $(table.cells().nodes()).removeClass('border-left');

    Object.values(BORDER_GROUPS).forEach(idxs => {
      const firstVisible = idxs.find(i => table.column(i).visible());
      if (firstVisible == null) return;

      const $boundTH = $(table.column(firstVisible).header());
      const pos = $boundTH.index();

      $boundTH.addClass('border-left');
      $labelDisp.children().eq(pos).addClass('border-left');

      $boundOrig.children().eq(pos).addClass('border-left');
      $labelOrig.children().eq(pos).addClass('border-left');

      $(table.column(firstVisible).nodes()).addClass('border-left');
    });
  }

  // ===== 5) DATATABLE INIT & LIFECYCLE =======================================
  function initializeDataTable(tableSelector, data) {
    $(tableSelector).DataTable({
      dom:
        "<'row'<'col-sm-9 dt-btn'B><'col-sm-3 text-end'f>>" +
        "<'row'<'col-sm-12'tr>>" +
        "<'row'<'col-sm-12'ip>>",
      order: [2, "asc"],
      pageLength: 50,
      data,
      columns,
      columnDefs,
      autoWidth: true,
      processing: false,
      deferRender: true,
      paging: true,
      scrollX: true,
      scrollY: "55vh",
      scrollCollapse: true,
      fixedHeader: false,
      stateSave: true,
      stateDuration: 0,
      stateSaveParams: function (settings, state) {
        state.columns = state.columns.map(col => ({ visible: col.visible }));
        state.order = settings.aaSorting;
        state.search = { search: "" };
        state.start = 0;
        state.length = 50;
        state.columns.forEach(c => c.search = { search: "" });
      },
      stateLoadParams: function (settings, state) {
        state.search = { search: "" };
        state.columns.forEach(c => { c.search = { search: "" }; });
        state.start = 0;
      },
      buttons: [
        {
          text: 'Reset filters',
          className: 'cluster-a',
          attr: { id: 'btnReset' },
          action: function (e, dt, node) {
            const $btn = $(node).closest('.dt-button');
            showBtnSpinner($btn);
            requestAnimationFrame(() => {
              requestAnimationFrame(() => {
                reset_all().finally(() => hideBtnSpinner($btn));
              });
            });
          }
        },
        {
          text: 'Show Representatives Only',
          attr: { id: "rep_toggle_btn" },
          className: 'cluster-a',
          action: function (e, dt, node) {
            repFilterOn = !repFilterOn;
            dt.draw();
            if (repFilterOn) {
              $(node).addClass("active-rep").text("Show All Structures");
            } else {
              $(node).removeClass("active-rep").text("Show Representatives Only");
            }
          }
        },
        {
          text: 'Align seqs',
          className: 'cluster-sep cluster-b',
          action: function () { handleAlignmentClick("#StructureBrowserTable"); }
        },
        {
          text: 'Download PDBs',
          className: 'cluster-b',
          action: function () { handleDownloadClick("#StructureBrowserTable"); }
        },
        (function () {
          const config = getSuperposeButtonConfig();
          return {
            text: config.label,
            attr: {
              id: "superpose_btn",
              class: (config.class ? config.class + " cluster-b" : "cluster-b"),
              "data-type": config.type
            },
            action: function () { handleSuperpositionClick("#StructureBrowserTable", config.type); }
          };
        })(),
        {
          extend: 'collection',
          text: 'Export Excel',
          className: 'cluster-b',
          autoClose: true,
          name: 'exportExcel',
          buttons: [
            { text: 'Full table (filtered)',  name: 'exportAll',
              action: function (e, dt) { runWithExportSpinner(dt, () => exportToXlsx('all')); } },
            { text: 'Selected rows',          name: 'exportSelected',
              action: function (e, dt) { runWithExportSpinner(dt, () => exportToXlsx('selected')); } }
          ]
        }
      ]
    });
  }

  // ===== 6) FILTERS ===========================================================
  function createFilters() {
    const dt = $("#StructureBrowserTable").DataTable();
    let column_filters = [];
    column_filters = column_filters.concat(CreateColumnFilters(dt,  2,10,"Multi-select-exact"));
    column_filters = column_filters.concat(CreateColumnFilters(dt, 12, 1,"Range-float-vertical"));
    column_filters = column_filters.concat(CreateColumnFilters(dt, 13, 2,"Multi-select-exact"));
    column_filters = column_filters.concat(CreateColumnFilters(dt, 15, 2,"Range-float-vertical"));
    column_filters = column_filters.concat(CreateColumnFilters(dt, 17, 3,"Multi-select-exact"));
    column_filters = column_filters.concat(CreateColumnFilters(dt, 20, 1,"Range-float-vertical"));
    column_filters = column_filters.concat(CreateColumnFilters(dt, 21, 7,"Multi-select-exact-filter"));
    column_filters = column_filters.concat(CreateColumnFilters(dt, 28, 2,"Multi-select-exact"));
    column_filters = column_filters.concat(CreateColumnFilters(dt, 30, 2,"Multi-select-unspecific"));
    column_filters = column_filters.concat(CreateColumnFilters(dt, 32, 1,"Range-select-vertical"));
    createDropdownFilters(dt, column_filters);
  }

  // ===== 7) EXPORT TO EXCEL ===================================================
  const EXCLUDE_EXPORT_COLS = [0, 1, 5]; // checkbox, GPCRdb, hidden ID
  const SUPER_GROUPS = [
    { label: 'RECEPTOR',             idxs: [2,3,4] },
    { label: 'CLASSIFICATION',       idxs: [6,7,8] },
    { label: 'STRUCTURE',            idxs: [9,10,11,12,13,14,15,16] },
    { label: 'SIGNAL PROTEIN',       idxs: [17,18,19,20] },
    { label: 'AUXILIARY PROTEIN',    idxs: [21,22] },
    { label: 'STRUCTURE LIGAND',     idxs: [23,24,25] },
    { label: 'PHYSIOLOGICAL LIGAND', idxs: [26,27] },
    { label: 'SODIUM ION SITE',      idxs: [28,29] },
    { label: 'REFERENCE',            idxs: [30,31,32] }
  ];
  const REFINED_COL   = 11;
  const REFERENCE_COL = 31;

  function buildGroupRowAndMerges(visibleCols) {
    const row = Array(visibleCols.length).fill(null);
    const merges = [];
    SUPER_GROUPS.forEach(g => {
      const vis = visibleCols.filter(c => g.idxs.indexOf(c) !== -1);
      if (!vis.length) return;
      const cStart = visibleCols.indexOf(vis[0]);
      const cEnd   = visibleCols.indexOf(vis[vis.length - 1]);
      row[cStart] = g.label;
      if (cEnd > cStart) merges.push({ s:{ r:0, c:cStart }, e:{ r:0, c:cEnd } });
    });
    return { groupRow: row, merges };
  }

  function exportToXlsx(mode) {
    const table    = $('#StructureBrowserTable').DataTable();
    const settings = table.settings()[0];
    const REFINED_BASE  = 'https://gpcrdb.org/structure/';

    const visibleCols = table.columns(':visible').indexes().toArray()
      .filter(i => EXCLUDE_EXPORT_COLS.indexOf(i) === -1);
    if (!visibleCols.length) { alert('No visible columns to export.'); return; }

    const rowIdxs = (mode === 'selected')
      ? (() => {
          const idxs = [];
          table.rows({ search:'applied' }).every(function(){
            const d = this.data();
            if (d && globalSelectedIds.has(String(d.id))) idxs.push(this.index());
          });
          return idxs;
        })()
      : table.rows({ search:'applied' }).indexes().toArray();
    if (!rowIdxs.length) { alert(mode === 'selected' ? 'No selected rows.' : 'No rows to export.'); return; }

    const htmlToText = (html) => {
      const el = document.createElement('div');
      el.innerHTML = String(html ?? '').replace(/<br\s*\/?>/gi, '\n');
      return (el.textContent || el.innerText || '').trim();
    };
    const hrefFromAnchor = (html) => {
      const m = /<a [^>]*href=["']([^"']+)["'][^>]*>/i.exec(html || '');
      return m ? m[1] : null;
    };
    const toAbsoluteRefined = (href) => {
      if (!href) return '';
      if (/^https?:\/\//i.test(href)) return href;
      const t = href.replace(/^\/+/, '');
      if (/^structure\//i.test(t)) return 'https://gpcrdb.org/' + t;
      if (/^refined\//i.test(t))   return REFINED_BASE + t;
      return REFINED_BASE + t;
    };

    const { groupRow, merges } = buildGroupRowAndMerges(visibleCols);

    const $theadOrig = $('#StructureBrowserTable thead');
    const $titleRow  = $theadOrig.find('tr').eq(-2);
    const getHeaderLabel = (ci) => {
      const colCfg = settings.aoColumns[ci] || {};
      if (colCfg.sName) return String(colCfg.sName);
      const $ths = $titleRow.children('th');
      if (ci < $ths.length) {
        const $th = $($ths[ci]).clone();
        $th.find('button, i, .glyphicon, .icon-button').remove();
        const t = $th.text().replace(/\s+/g, ' ').trim();
        if (t) return t;
      }
      const th = colCfg.nTh || table.column(ci).header();
      return th ? $(th).text().replace(/\s+/g, ' ').trim() : `Column ${ci}`;
    };
    const headerTitles = visibleCols.map(getHeaderLabel);

    const body = rowIdxs.map(ri =>
      visibleCols.map(ci => {
        const disp = table.cell(ri, ci).render('display');
        const exp  = table.cell(ri, ci).render('export');

        if (ci === REFINED_COL) {
          let href = null;
          if (typeof disp === 'string' && /<a\s/i.test(disp)) href = hrefFromAnchor(disp);
          if (!href && typeof exp  === 'string' && /<a\s/i.test(exp))  href = hrefFromAnchor(exp);
          if (href) return toAbsoluteRefined(href);
          return htmlToText(exp);
        }
        if (ci === REFERENCE_COL) return htmlToText(exp);
        return htmlToText(exp);
      })
    );

    const aoa = [groupRow, headerTitles, ...body];
    const ws  = XLSX.utils.aoa_to_sheet(aoa);
    ws['!merges'] = merges;
    ws['!cols'] = headerTitles.map((h, j) => {
      const maxLen = Math.max(String(h||'').length, ...body.map(r => String(r[j]||'').length));
      return { wch: Math.min(Math.max(10, maxLen + 2), 60) };
    });

    const wb = XLSX.utils.book_new();
    XLSX.utils.book_append_sheet(wb, ws, 'Structures');

    const stamp = new Date().toISOString().slice(0,10);
    const fname = mode === 'selected'
      ? `GPCRdb_structures_selected_${stamp}.xlsx`
      : `GPCRdb_structures_full_${stamp}.xlsx`;

    XLSX.writeFile(wb, fname);
  }

  // make export function reachable by DT button actions
  window.exportToXlsx = exportToXlsx;

  // ===== 8) ROW SELECTION & CHECKBOX SYNC ====================================
  function setRowSelected($tr, selected){
    const table = $("#StructureBrowserTable").DataTable();
    const row   = table.row($tr);
    const data  = row.data();
    if (!data) return;

    const id = String(data.id);
    $tr.toggleClass('row-selected', selected);
    $tr.find('input.select-row').prop('checked', selected);

    if (selected) globalSelectedIds.add(id);
    else          globalSelectedIds.delete(id);

    updateExportSelectedState();
  }

  function syncSelectAll(){
    const table = $("#StructureBrowserTable").DataTable();
    const $rows = $(table.rows({ search:'applied', page:'current' }).nodes());
    const allSelected = $rows.length > 0 && $rows.filter('.row-selected').length === $rows.length;
    $('#select-all').prop('checked', allSelected);
  }

  function applySelectionStyles(){
    const table = $("#StructureBrowserTable").DataTable();
    table.rows({ page: 'current' }).every(function(){
      const id   = String(this.data()?.id);
      const sel  = globalSelectedIds.has(id);
      const $tr  = $(this.node());
      $tr.toggleClass('row-selected', sel);
      $tr.find('input.select-row').prop('checked', sel);
    });
  }

  // ===== 9) PAGE ACTIONS (align/download/superpose + clipboard) ==============
  function getSelectedRowData(tableId) {
    const tableSelector = typeof tableId === "string" ? tableId : `#${tableId.id || tableId}`;
    const table = $(tableSelector).DataTable();
    const allData = table.rows().data().toArray();
    return allData.filter(row => globalSelectedIds.has(String(row.id)));
  }

  function handleAlignmentClick(tableId) {
    const selectedData = getSelectedRowData(tableId);
    ClearSelection("targets");
    if (selectedData.length === 0) { showAlert("No entries selected for sequence alignment", "danger"); return; }
    selectedData.forEach(row => AddToSelection("targets", "structure", row.pdb));
    window.location.href = "/structure/selection_convert";
  }

  function handleDownloadClick(tableId) {
    const selectedData = getSelectedRowData(tableId);
    ClearSelection("targets");
    if (selectedData.length === 0) { showAlert("No structures selected for download", "danger"); return; }
    if (selectedData.length > 100){ showAlert("Maximum number of selected entries is 100", "warning"); return; }
    const ids = selectedData.map(row => row.pdb.replace(/\s+/g, ""));
    AddToSelection("targets", "structure_many", ids.join(","));
    window.location.href = "/structure/pdb_download";
  }

  function getSuperposeButtonConfig() {
    const hash = window.location.hash;
    return {
      label: hash === "#keepselectionreference" ? "Add reference to superposition"
           : hash === "#keepselectiontargets"   ? "Add targets to superposition"
                                                : "Superposition",
      type:  hash === "#keepselectionreference" ? "reference"
           : hash === "#keepselectiontargets"   ? "targets"
                                                : "default",
      class: (hash === "#keepselectionreference" || hash === "#keepselectiontargets")
              ? "dt-button active-superpose" : "dt-button"
    };
  }

  function handleSuperpositionClick(tableId, type = "default") {
    const selectedData = getSelectedRowData(tableId);

    if (type === "reference") {
      if (selectedData.length !== 1) { showAlert("Please select exactly one entry as reference", "danger"); return; }
      const refPdb = selectedData[0].pdb.replace(/\s+/g, "");
      AddToSelection("reference", "structure", refPdb);
      window.location.href = "/structure/superposition_workflow_index#keepselection";
      return;
    }

    if (type === "targets") {
      if (selectedData.length === 0) { showAlert("No targets selected", "danger"); return; }
      if (selectedData.length > 50){ showAlert("Maximum number of selected targets is 50", "warning"); return; }
      const ids = selectedData.map(row => row.pdb.replace(/\s+/g, ""));
      AddToSelection("targets", "structure_many", ids.join(","));
      window.location.href = "/structure/superposition_workflow_index#keepselection";
      return;
    }

    // default: open modal
    superpositionModern(tableId,
      ["pdb", "entry_short", "iuphar_name", "family", "class", "species", "state", "pub_date"],
      "structure_browser", "gpcr", "pdb"
    );
  }

  function superpositionModern(tableId, columns, site, source = 'gpcr', structure_column_key = "pdb", hidden_columns = []) {
    ClearSelection('targets'); ClearSelection('reference');

    const selectedData = getSelectedRowData(tableId);
    if (selectedData.length === 0) { showAlert("No entries selected for superposition", "danger"); return; }
    if (selectedData.length > 50){ showAlert("Maximum number of selected entries is 50", "warning"); return; }

    const selected_ids = [];
    let selection_type = '';

    if (site === 'structure_browser') {
      selectedData.forEach(row => { if (row[structure_column_key]) selected_ids.push(row[structure_column_key].replace(/\s+/g, '')); });
      selection_type = 'structure_many';
    }

    // Modal show/hide
    const modal = document.getElementById("superposition-modal");
    const span = document.getElementById("close_superposition_modal");
    modal.style.display = "block";
    span.onclick = () => modal.style.display = "none";
    window.onclick = (e) => { if (e.target === modal) modal.style.display = "none"; };

    // Populate modal table
    const tbody = document.querySelector("#superposition_modal_table tbody");
    tbody.innerHTML = "";
    selectedData.forEach(row => {
      const tr = document.createElement("tr");
      const checkbox = document.createElement("input");
      checkbox.type = "checkbox";
      const checkboxTd = document.createElement("td");
      checkboxTd.appendChild(checkbox);
      tr.appendChild(checkboxTd);
      columns.forEach(key => {
        const td = document.createElement("td");
        td.innerHTML = row[key] || '';
        if (hidden_columns.includes(key)) td.style.display = "none";
        tr.appendChild(td);
      });
      tbody.appendChild(tr);
    });

    // Click row → set reference + targets then redirect
    $("#superposition_modal_table tbody tr").on("click", function () {
      const cells = $(this).children();
      const ref_id = cells.eq(columns.indexOf(structure_column_key) + 1).text().replace(/\s+/g, '');
      AddToSelection('reference', 'structure', ref_id);
      const targets = selected_ids.filter(id => id !== ref_id);
      AddToSelection('targets', selection_type, targets.join(","));
      $(this).children(':first').find("input").prop("checked", true);
      const redirectUrl = (source === 'gpcr') ? '/structure/superposition_workflow_index#keepselection'
                                              : '/structure/superposition_workflow_gprot_index#keepselection';
      window.location.href = redirectUrl;
    });
  }

  // Clipboard helpers (delegated so it works regardless of load order)
  function copySelectedToClipboard(key, label) {
    const table = $("#StructureBrowserTable").DataTable();
    const selectedRows = [];
    table.rows().every(function () {
      const data = this.data();
      if (data && globalSelectedIds.has(String(data.id))) selectedRows.push(data);
    });
    if (!selectedRows.length) { showAlert("No entries selected for copying", "danger"); return; }

    const values = selectedRows.map(r => r[key]).filter(v => v && String(v).trim() !== "");
    if (!values.length) { showAlert(`No ${label} available for selected rows`, "warning"); return; }

    const out = values.join("\n");
    if (navigator.clipboard?.writeText) {
      navigator.clipboard.writeText(out)
        .then(() => showAlert(`${label} copied to clipboard!`, "success"))
        .catch(() => fallbackCopy(out, label));
    } else {
      fallbackCopy(out, label);
    }
  }
  function fallbackCopy(text, label) {
    const textarea = document.createElement("textarea");
    textarea.value = text;
    document.body.appendChild(textarea);
    textarea.select();
    document.execCommand("copy");
    document.body.removeChild(textarea);
    showAlert(`${label} copied to clipboard!`, "success");
  }

  // Toggle Export "Selected rows" enablement based on selection size
  function updateExportSelectedState(){
    const table = $("#StructureBrowserTable").DataTable();
    const hasAny = window.globalSelectedIds.size > 0;
    const btn = table.button('exportSelected:name');
    if (btn && btn.length) btn.enable(hasAny);
  }
  window.updateExportSelectedState = updateExportSelectedState;

  // ===== 10) MISC UI UTILITIES ===============================================
  function reset_all() {
    return new Promise((resolve) => {
      $('.select2').val(null).trigger('change');
      const table = $("#StructureBrowserTable").DataTable();
      table.one('draw.dt', () => resolve());
      table.draw(false);
      setTimeout(resolve, 200);
    });
  }
  function showBtnSpinner($btn){
    if (!$btn.children('.btn-spinner').length){
      $btn.append('<span class="btn-spinner" aria-hidden="true"></span>');
    }
    $btn.addClass('is-busy');
  }
  function hideBtnSpinner($btn){ $btn.removeClass('is-busy'); }
  function runWithExportSpinner(dt, workFn) {
    const api   = dt.button('exportExcel:name');
    const $main = api ? $(api.node()) : $();

    if ($main.length && $main.hasClass('is-busy')) {
      // Already busy: return a resolved promise so caller can still chain .then()
      return Promise.resolve();
    }

    if ($main.length) {
      showBtnSpinner($main);
    }

    return new Promise(res => requestAnimationFrame(() => requestAnimationFrame(res)))
      .then(() => Promise.resolve(workFn()))
      .finally(() => {
        setTimeout(() => {
          if ($main.length) hideBtnSpinner($main);
        }, 150);
      });
  }


  // Expose a few needed funcs (used by DT button configs)
  window.handleAlignmentClick     = handleAlignmentClick;
  window.handleDownloadClick      = handleDownloadClick;
  window.getSuperposeButtonConfig = getSuperposeButtonConfig;
  window.handleSuperpositionClick = handleSuperpositionClick;
  window.runWithExportSpinner     = runWithExportSpinner;
  window.reset_all                = reset_all;

  // ===== 11) INITIALIZATION ===================================================
  $(function () {
    // Busy loader while fetching data
    $("#Init_loader").show().busyLoad("show", {
      spinner: "accordion",
      text: "Loading data...",
      fontSize: "2.5rem",
      textClass: "custom-loader-text",
      textPosition: "top",
      textMargin: "-6rem",
      color: "black",
      background: "white"
    });

    // Custom filter hook (representatives only)
    $.fn.dataTable.ext.search.push(function(settings, data, dataIndex, rowData) {
      if (!repFilterOn) return true;
      return rowData.rep === 1;
    });

    // Fetch & build table
    fetch("/structure/data/")
      .then(resp => { if (!resp.ok) throw new Error(`HTTP ${resp.status}`); return resp.json(); })
      .then(structure_data => {
        initializeDataTable("#StructureBrowserTable", structure_data);
        const table = $("#StructureBrowserTable").DataTable();

        // Add the Table Control (ColVis) and wire up custom menu
        table.button().add(1, controlTableButtonWithGroups(table));
        table.buttons('.btn-table-control').nodes().each(function () {
          $(this).off('click.colvis8').on('click.colvis8', function () {
            setTimeout(function () {
              const $col  = $('div.dt-button-collection:visible');
              const $menu = $col.find('> div[role="menu"]');
              $col.find('.dt-button-collection-title').remove();
              $col.addClass('columns-8');
              buildColvisMenu($menu, table);
            }, 0);
          });
        });

        // Keep dividers in sync with visibility toggles
        table.on('column-visibility.dt', function(){ applyGroupBorders(table); });

        // Build the column filters
        createFilters();

        // Finish UI
        $("#Init_loader").busyLoad("hide").hide();
        $("#StructureBrowserTableContainer").show();
        table.columns.adjust().draw(false);

        // initial export state
        updateExportSelectedState();
      })
      .catch(err => {
        console.error("Failed to fetch structure data:", err);
        $("#StructureBrowserTable").after(
          `<div class="alert alert-danger mt-3" role="alert">
             Could not load structure data. Please try again later or check your network connection.
           </div>`
        );
      });

    // Row click toggles selection (but ignore interactive elements)
    $('#StructureBrowserTable tbody').on('click', 'tr', function(e){
      if ($(e.target).closest('a,button,input,label,.select2-container,.dt-button').length) return;
      setRowSelected($(this), !$(this).hasClass('row-selected'));
      syncSelectAll();
    });

    // Header "select all"
    $('#select-all').on('click', function () {
      const table = $("#StructureBrowserTable").DataTable();
      const checked = this.checked;
      table.rows({ search: 'applied', page: 'current' }).every(function () {
        setRowSelected($(this.node()), checked);
      });
      syncSelectAll();
      updateExportSelectedState();
    });

    // Individual checkbox
    $(document).on('change', '.select-row', function (e) {
      setRowSelected($(this).closest('tr'), this.checked);
      syncSelectAll();
      updateExportSelectedState();
      e.stopPropagation();
    });

    // Clipboard buttons (delegated)
    // Intercept on capture phase so DT never sees the event (stops sorting/reorder)
    (function attachExportCopyInterceptors() {
    // Click → do the copy and block sort
    document.addEventListener('click', function (e) {
        const hit = e.target.closest('#uniprot_copy, #pdb_copy, .uniprot-export, .pdb-export');
        if (!hit) return;

        e.preventDefault();
        e.stopPropagation();
        if (typeof e.stopImmediatePropagation === 'function') e.stopImmediatePropagation();

        if (hit.id === 'uniprot_copy' || hit.classList.contains('uniprot-export')) {
        copySelectedToClipboard('entry_name', 'UniProt IDs');
        } else {
        copySelectedToClipboard('pdb', 'PDB IDs');
        }
    }, /* useCapture: */ true);

    // Mousedown → block column drag/reorder (if ColReorder is present)
    document.addEventListener('mousedown', function (e) {
        const hit = e.target.closest('#uniprot_copy, #pdb_copy, .uniprot-export, .pdb-export');
        if (!hit) return;

        e.preventDefault();
        e.stopPropagation();
        if (typeof e.stopImmediatePropagation === 'function') e.stopImmediatePropagation();
    }, /* useCapture: */ true);
    })();

    // Bootstrap popovers
    $('[data-toggle="popover"]').popover({
      container: 'body',
      html: true,
      sanitize: false,
      template:
        '<div class="popover popover-wide" role="tooltip">' +
          '<div class="arrow"></div>' +
          '<h3 class="popover-title"></h3>' +
          '<div class="popover-content"></div>' +
        '</div>'
    });
    $(document).on('shown.bs.popover', function () {
      $('.popover').last().addClass('popover-wide');
    });
  });

})(jQuery, window, document);

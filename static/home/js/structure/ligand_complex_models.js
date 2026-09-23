/* ==========================================================================
   Ligand Complex Models – Modern JS
   Data comes from /structure/ligand_complex_models/data/ (JSON), rendered
   client-side with DataTables + NorgesDTFilterBuilder -- same building
   blocks as static/home/js/structure/structure_browser_modern.js, trimmed
   down for this page (no ColVis, no sticky columns, no Excel export; only
   Reset filters / Align / Download, wired up as DataTables Buttons so they
   get the same spread-out styling as the main structure browser).
   Dependencies: jQuery, DataTables (NewDrugsBrowser_datatables bundle),
                 NorgesDTFilterBuilder, Select2, BusyLoad, gpcrdb.js
                 (showAlert), browser_functions.js (ClearSelection/AddToSelection).
   ========================================================================== */

(function ($, window, document) {
  "use strict";

  const TABLE_SELECTOR = "#homology_models";
  const selectedIds = new Set();

  // ===== 1) RENDERING HELPERS =================================================
  function escapeAttr(s) {
    return String(s == null ? "" : s)
      .replace(/&/g, "&amp;")
      .replace(/"/g, "&quot;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;");
  }

  function capfirst(s) {
    return s ? s.charAt(0).toUpperCase() + s.slice(1) : s;
  }

  function ligandDisplayName(name) {
    if (!name) return "";
    return name.length >= 6 ? capfirst(name) : name.toUpperCase();
  }

  function receptorShort(name) {
    if (!name) return "";
    let out = name;
    if (!out.startsWith("mGlu")) out = out.charAt(0).toUpperCase() + out.slice(1);
    return out.replace(/ receptor/g, "").replace(/-adrenoceptor/g, "");
  }

  function cutClassname(name) {
    if (!name) return "-";
    return name.startsWith("Class ") ? name.slice(5) : name;
  }

  function scoreDisplay(d) {
    return d === null || d === undefined || d === "" ? "-" : Number(d).toFixed(2);
  }

  // Collapse to one line with ellipsis, full text on hover -- see td.expand_1line CSS
  function wrapExpand(html) {
    if (!html) return html;
    return `<span class="collapsed">${html}</span><span class="expanded">${html}</span>`;
  }

  // ===== 2) DATATABLES – COLUMNS & DEFS =======================================
  const columnDefs = [
    { targets: [0, 1], orderable: false, searchable: false },
    { targets: [2], createdCell: (td) => td.classList.add("expand_1line") },
  ];

  const columns = [
    // checkbox
    {
      data: "id",
      name: "Select",
      width: "34px",
      render: (d, t) => (t === "display" ? `<input type="checkbox" class="select-row" value="${d}">` : d),
    },
    // model icon
    {
      data: null,
      name: "Model",
      width: "34px",
      render: (d, t, row) => {
        if (t !== "display") return "";
        const icon = row.has_signprot ? window.GPCR_STRUCTURE_LOGO_URL : window.MODEL_SURFACE_URL;
        const href = row.pdb_code_index ? `ligand_complex_models/${row.pdb_code_index}` : "#";
        return `<a target="_blank" href="${href}"><img width="14px" height="20px" class="model-link" src="${icon}"></a>`;
      },
    },

    // ----- LIGAND -----
    {
      data: "ligand",
      name: "Name",
      width: "160px",
      render: (d, t) => {
        if (!d) return t === "display" ? "-" : "";
        if (t !== "display") return d.name || "";
        const label = ligandDisplayName(d.name);
        const link = d.type_slug === "peptide"
          ? `<a class="struct" href="/ligand/${d.id}/group" data-ligand-type="peptide" data-sequence="${escapeAttr(d.sequence)}">${label}</a>`
          : `<a class="struct" href="/ligand/${d.id}/group" data-ligand-type="${d.type_slug || "other"}" data-smiles="${escapeAttr(d.smiles)}" rel="${d.picture || "Not_available"}">${label}</a>`;
        return wrapExpand(link);
      },
    },
    { data: "mol_modality", name: "Mol. modality", width: "80px" },
    { data: "pharm_modality", name: "Pharm. modality", width: "90px" },
    { data: "physiological", name: "Physiological", width: "70px" },
    { data: "clinical", name: "Clinical", width: "90px" },

    // ----- RECEPTOR -----
    {
      data: "gene",
      name: "Gene",
      width: "90px",
      render: (d, t) => {
        if (!d) return t === "display" ? "-" : "";
        if (t !== "display") return d.name || "";
        return d.entrez_weblink ? `<a href="${d.entrez_weblink}" target="_blank">${d.name}</a>` : d.name;
      },
    },
    {
      data: "protein",
      name: "Protein",
      width: "140px",
      render: (d, t) => {
        if (!d) return "";
        const label = receptorShort(d.name);
        return t === "display" ? `<a href="/protein/${d.entry_name}">${label}</a>` : label;
      },
    },
    { data: "family", name: "Family", width: "140px" },
    { data: "class", name: "Class", width: "110px", render: (d) => cutClassname(d) },
    { data: "modality", name: "Modality", width: "140px" },
    { data: "chemotype", name: "Chemotype", width: "120px" },

    // ----- STRUCTURE -----
    { data: "state", name: "State", width: "70px", render: (d) => d || "-" },
    { data: "signal_protein_family", name: "Signal protein family", width: "90px", render: (d) => d || "-" },
    { data: "signal_protein_subtype", name: "Signal protein subtype", width: "90px", render: (d) => d || "-" },

    // ----- LIGAND MODEL SCORE -----
    { data: "af2_score", name: "AF2 (PAE mean)", width: "70px", render: (d) => scoreDisplay(d) },
    { data: "boltz2_score", name: "Boltz2 (Ligand pLDDT)", width: "80px", render: (d) => scoreDisplay(d) },
    { data: "rfaa_score", name: "RFAA (pLDDT mean)", width: "70px", render: (d) => scoreDisplay(d) },

    // ----- DATE -----
    { data: "publication_date", name: "Date", width: "90px", render: (d) => d || "-" },
  ];

  // ===== 3) BUTTON SPINNER HELPERS ============================================
  function showBtnSpinner($btn) {
    if (!$btn.children(".btn-spinner").length) {
      $btn.append('<span class="btn-spinner" aria-hidden="true"></span>');
    }
    $btn.addClass("is-busy");
  }
  function hideBtnSpinner($btn) {
    $btn.removeClass("is-busy");
  }

  function reset_all() {
    return new Promise((resolve) => {
      $(".select2").val(null).trigger("change");
      const table = $(TABLE_SELECTOR).DataTable();
      table.one("draw.dt", () => resolve());
      table.draw(false);
      setTimeout(resolve, 200);
    });
  }

  // ===== 4) ROW SELECTION ======================================================
  function setRowSelected($tr, selected) {
    const table = $(TABLE_SELECTOR).DataTable();
    const data = table.row($tr).data();
    if (!data) return;
    const id = String(data.id);
    $tr.toggleClass("row-selected", selected);
    $tr.find("input.select-row").prop("checked", selected);
    if (selected) selectedIds.add(id);
    else selectedIds.delete(id);
  }

  function syncSelectAll() {
    const table = $(TABLE_SELECTOR).DataTable();
    const $rows = $(table.rows({ search: "applied", page: "current" }).nodes());
    const allSelected = $rows.length > 0 && $rows.filter(".row-selected").length === $rows.length;
    $("#select-all").prop("checked", allSelected);
  }

  function getSelectedRowData() {
    const table = $(TABLE_SELECTOR).DataTable();
    return table.rows().data().toArray().filter((row) => selectedIds.has(String(row.id)));
  }

  // ===== 5) PAGE ACTIONS (align/download) =====================================
  function handleAlignClick() {
    const selected = getSelectedRowData();
    if (selected.length === 0) {
      showAlert("No entries selected for alignment", "danger");
      return;
    }
    ClearSelection("targets");
    selected.forEach((row) => {
      if (row.pdb_code_index) AddToSelection("targets", "structure", row.pdb_code_index);
    });
    window.location.href = "/structure/selection_convert";
  }

  function handleDownloadClick() {
    const selected = getSelectedRowData();
    if (selected.length === 0) {
      showAlert("No models selected for download", "danger");
      return;
    }
    const ids = selected.map((row) => row.id);
    window.location.href = "/structure/lig_complexmod_download?ids=" + ids.join(",");
  }

  // ===== 6) DATATABLE INIT =====================================================
  function initializeDataTable(tableSelector, data) {
    $(tableSelector).DataTable({
      dom:
        "<'row'<'col-sm-9 dt-btn'B><'col-sm-3 text-end'f>>" +
        "<'row'<'col-sm-12'tr>>" +
        "<'row'<'col-sm-12'ip>>",
      data,
      columns,
      columnDefs,
      order: [[2, "asc"], [7, "asc"]],
      paging: true,
      pageLength: 50,
      scrollX: true,
      scrollY: "65vh",
      scrollCollapse: true,
      autoWidth: true,
      deferRender: true,
      buttons: [
        {
          text: "Reset filters",
          className: "cluster-a",
          attr: { id: "btnReset" },
          action: function (e, dt, node) {
            const $btn = $(node).closest(".dt-button");
            showBtnSpinner($btn);
            requestAnimationFrame(() => {
              requestAnimationFrame(() => {
                reset_all().finally(() => hideBtnSpinner($btn));
              });
            });
          },
        },
        {
          text: "Align",
          className: "cluster-sep cluster-b",
          action: function () { handleAlignClick(); },
        },
        {
          text: "Download",
          className: "cluster-b",
          action: function () { handleDownloadClick(); },
        },
      ],
    });
  }

  // ===== 7) FILTERS ============================================================
  function createFilters() {
    const dt = $(TABLE_SELECTOR).DataTable();
    let column_filters = [];
    column_filters = column_filters.concat(CreateColumnFilters(dt, 2, 1, "Multi-select-exact-filter"));   // Name (html)
    column_filters = column_filters.concat(CreateColumnFilters(dt, 3, 4, "Multi-select-exact"));          // Mol./Pharm. modality, Physiological, Clinical
    column_filters = column_filters.concat(CreateColumnFilters(dt, 7, 1, "Multi-select-exact"));          // Gene
    column_filters = column_filters.concat(CreateColumnFilters(dt, 8, 1, "Multi-select-exact-filter"));   // Protein (html)
    column_filters = column_filters.concat(CreateColumnFilters(dt, 9, 4, "Multi-select-exact"));          // Family, Class, Modality, Chemotype
    column_filters = column_filters.concat(CreateColumnFilters(dt, 13, 3, "Multi-select-exact"));         // State, Signal protein family/subtype
    column_filters = column_filters.concat(CreateColumnFilters(dt, 16, 3, "Range-float-vertical"));       // Scores
    createDropdownFilters(dt, column_filters);
  }

  // ===== 8) INITIALIZATION =====================================================
  $(function () {
    $("#Init_loader").show().busyLoad("show", {
      spinner: "accordion",
      text: "Loading data...",
      fontSize: "2.5rem",
      textClass: "custom-loader-text",
      textPosition: "top",
      textMargin: "-6rem",
      color: "black",
      background: "white",
    });

    fetch("/structure/ligand_complex_models/data/")
      .then((resp) => {
        if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
        return resp.json();
      })
      .then((data) => {
        initializeDataTable(TABLE_SELECTOR, data);
        createFilters();

        $("#Init_loader").busyLoad("hide").hide();
        $("#homology_models_div").show();
        $(TABLE_SELECTOR).DataTable().columns.adjust();
      })
      .catch((err) => {
        console.error("Failed to fetch ligand complex model data:", err);
        $("#Init_loader").busyLoad("hide").hide();
        $("#browser").append(
          '<div class="alert alert-danger mt-3" role="alert">Could not load ligand complex model data. Please try again later or check your network connection.</div>'
        );
      });

    $(TABLE_SELECTOR).on("click", "tbody tr", function (e) {
      if ($(e.target).closest("a,button,input,label,.select2-container").length) return;
      setRowSelected($(this), !$(this).hasClass("row-selected"));
      syncSelectAll();
    });

    $(document).on("change", ".select-row", function (e) {
      setRowSelected($(this).closest("tr"), this.checked);
      syncSelectAll();
      e.stopPropagation();
    });

    $("#select-all").on("click", function () {
      const table = $(TABLE_SELECTOR).DataTable();
      const checked = this.checked;
      table.rows({ search: "applied", page: "current" }).every(function () {
        setRowSelected($(this.node()), checked);
      });
      syncSelectAll();
    });
  });
})(jQuery, window, document);

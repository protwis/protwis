/* ==========================================================================
   Interaction Landing Page – Structure Picker
   A reduced, single-select subset of the structure browser table, used to
   pick a PDB code before jumping to /interaction/<pdb>.
   Dependencies: jQuery, DataTables, BusyLoad, NorgesDTFilterBuilder
   ========================================================================== */

(function ($, window, document) {
  "use strict";

  let selectedPdb = null;

  // ===== RENDERING HELPERS (trimmed copies of structure_browser_modern.js) ===
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
      ? (hrefBuilder ? hrefBuilder(first.id, first) : `/ligand/${first.id}/info`)
      : null;

    const collapsed = firstHref
      ? `<a href="${firstHref}" target="_blank" rel="noopener">${first.name || ''}</a>`
      : (first.name || '');

    const more = arr.length - 1;
    const hint = (showCountHint && more > 0) ? ` <span class="muted">(+${more} more)</span>` : '';

    const expanded = arr.map(x => {
      const nm = x.name || '';
      if (x.id == null || x.id === '' || x.id === '-') return nm;
      const href = hrefBuilder ? hrefBuilder(x.id, x) : `/ligand/${x.id}/info`;
      return `<a href="${href}" target="_blank" rel="noopener">${nm}</a>`;
    }).join('<br>');

    return (
      '<span class="expand_1line">' +
        '<span class="collapsed">' + collapsed + hint + '</span>' +
        '<span class="expanded">'  + expanded  + '</span>' +
      '</span>'
    );
  }

  // ===== COLUMNS ===============================================================
  const columnDefs = [
    {
      targets: "_all",
      className: "dt-head-nowrap dt-head-center dt-center dt-body-center",
      createdCell: function (td) { td.classList.add("expand_1line"); }
    },
    { targets: [0], orderable: false, searchable: false },
  ];

  const columns = [
    // SELECT
    {
      data: "pdb",
      name: "Select",
      orderable: false,
      searchable: false,
      render: (d) => `<input type="radio" name="interaction-select" class="select-row" value="${d}">`,
    },

    // RECEPTOR
    {
      data: "gpcrdb_link",
      name: "GPCRdb",
      orderable: false,
      searchable: false,
      render: (d, t) =>
        t === "display" && d
          ? `<a href="${d}" target="_blank" rel="noopener"><img class="gpcrdb-link" src="/static/home/logo/gpcr/main.png" width="12" height="12" alt="GPCRdb"></a>`
          : "",
    },
    {
      data: "entry_name",
      name: "UniProt",
      render: (d, t, r) =>
        t === "display"
          ? r?.uniprot_link
            ? `<a href="${r.uniprot_link}" target="_blank" rel="noopener">${d ? d.split("_")[0].toUpperCase() : "-"}</a>`
            : d
              ? d.split("_")[0].toUpperCase()
              : "-"
          : t === "filter" || t === "sort"
            ? d
              ? d.split("_")[0].toUpperCase()
              : ""
            : d || "",
    },
    {
      data: "gene",
      name: "Gene",
      render: (d, t) => {
        if (!d) return "";
        const name = d.name || "";
        const url = d.entrez_url;
        if (t !== "display" || !url || url === "-") return name;
        return `<a href="${url}" target="_blank" rel="noopener">${name}</a>`;
      },
    },
    {
      data: "iuphar_name",
      name: "Protein",
      render: (d, t, r) => {
        if (t !== "display") return d;
        const url = r && r.iuphar_link;
        return url && url !== "-"
          ? `<a href="${url}" target="_blank" rel="noopener">${d}</a>`
          : d;
      },
    },
    { data: "id", name: "ID", visible: false },

    // CLASSIFICATION
    { data: "family", name: "Receptor family" },
    { data: "class", name: "Class" },
    { data: "species", name: "Species" },

    // STRUCTURE LIGANDS
    {
      data: "ligands",
      name: "Name",
      render: (d, t) =>
        expandFirstThenHoverLinkList(
          d, t, (id) => `/ligand/${id}/info`, { showCountHint: true },
        ),
    },
    {
      data: "ligand_type",
      name: "Type",
      render: (d, t) => expand1LineRender(d, t),
    },
    {
      data: "ligand_role",
      name: "Modality",
      render: (d, t) => expand1LineRender(d, t),
    },

    // STRUCTURE
    { data: "method", name: "Method" },
    {
      data: "pdb",
      name: "PDB",
      render: (d, t) =>
        t !== "display"
          ? d
          : d && d !== "-"
            ? `<a href="/structure/${d}" target="_blank" rel="noopener">${d}</a>`
            : d,
    },
    {
      data: "resolution",
      name: "Resolution",
      render: (d) => (d ? parseFloat(d).toFixed(1) : "-"),
    },
    { data: "state", name: "State" },
    {
      data: "active_pct",
      name: "Degree active (%)",
      render: (d) => (d != null ? Math.round(d) : "-"),
    },

    // SIGNAL PROTEIN
    { data: "arrestin_family", name: "Family" },
    {
      data: "arrestin_name",
      name: "Subtype",
      render: (d, t, r) =>
        t !== "display"
          ? d
          : r?.arrestin_entry && r.arrestin_entry !== "-"
            ? `<a href="/signprot/${r.arrestin_entry}/" target="_blank" rel="noopener">${d}</a>`
            : d,
    },

    // REFERENCE
    { data: "reference", name: "Reference" },
  ];

  // ===== DATATABLE INIT =========================================================
  function initializeDataTable(data) {
    $("#InteractionPickerTable").DataTable({
      dom:
        "<'row'<'col-sm-12'tr>>" +
        "<'row'<'col-sm-12'ip>>",
      order: [[2, "asc"]],
      pageLength: 25,
      data,
      columns,
      columnDefs,
      autoWidth: true,
      processing: false,
      deferRender: true,
      paging: true,
      scrollX: true,
      scrollY: "50vh",
      scrollCollapse: true,
    });
  }

  // ===== FILTERS =================================================================
  function createFilters() {
    const dt = $("#InteractionPickerTable").DataTable();
    let column_filters = [];
    column_filters = column_filters.concat(CreateColumnFilters(dt, 2, 7, "Multi-select-exact"));
    column_filters = column_filters.concat(CreateColumnFilters(dt, 9, 3, "Multi-select-exact-filter"));
    column_filters = column_filters.concat(CreateColumnFilters(dt, 12, 2, "Multi-select-exact"));
    column_filters = column_filters.concat(CreateColumnFilters(dt, 14, 1, "Range-float-vertical"));
    column_filters = column_filters.concat(CreateColumnFilters(dt, 15, 1, "Multi-select-exact"));
    column_filters = column_filters.concat(CreateColumnFilters(dt, 16, 1, "Range-float-vertical"));
    column_filters = column_filters.concat(CreateColumnFilters(dt, 17, 2, "Multi-select-exact"));
    column_filters = column_filters.concat(CreateColumnFilters(dt, 19, 1, "Multi-select-unspecific"));
    createDropdownFilters(dt, column_filters);
  }

  // ===== SELECTION ================================================================
  function setRowSelected($tr, pdb) {
    $("#InteractionPickerTable tbody tr.row-selected")
      .removeClass("row-selected")
      .find(".select-row").prop("checked", false);

    if (pdb && pdb !== "-") {
      $tr.addClass("row-selected");
      $tr.find(".select-row").prop("checked", true);
      selectedPdb = pdb;
    } else {
      selectedPdb = null;
    }

    $("#btn-go-interaction, #btn-open-interaction").toggleClass("btn-inactive", !selectedPdb);
  }

  function goToInteraction(newTab) {
    if (!selectedPdb) return;
    const url = "/interaction/" + selectedPdb;
    if (newTab) {
      window.open(url, "_blank");
    } else {
      window.location.href = url;
    }
  }

  // ===== INITIALIZATION ===========================================================
  $(function () {
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

    fetch("/structure/data/")
      .then(resp => { if (!resp.ok) throw new Error(`HTTP ${resp.status}`); return resp.json(); })
      .then(structure_data => {
        const rows = structure_data.filter(d => d.has_ligand_interactions);

        initializeDataTable(rows);
        const table = $("#InteractionPickerTable").DataTable();

        createFilters();

        $("#Init_loader").busyLoad("hide").hide();
        $("#InteractionPickerTableContainer").show();
        table.columns.adjust().draw(false);
      })
      .catch(err => {
        console.error("Failed to fetch structure data:", err);
        $("#Init_loader").busyLoad("hide").hide();
        $("#InteractionPickerTable").after(
          `<div class="alert alert-danger mt-3" role="alert">
             Could not load structure data. Please try again later or check your network connection.
           </div>`
        );
      });

    // Row click selects it (but ignore genuine interactive elements)
    $("#InteractionPickerTable tbody").on("click", "tr", function (e) {
      if ($(e.target).closest("a,button,label,.select2-container,.dt-button").length) return;
      const table = $("#InteractionPickerTable").DataTable();
      const rowData = table.row(this).data();
      if (!rowData) return;
      setRowSelected($(this), rowData.pdb);
    });

    $("#btn-go-interaction").on("click", function () { goToInteraction(false); });
    $("#btn-open-interaction").on("click", function () { goToInteraction(true); });

    // Bootstrap popovers (State / Degree active info icons)
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

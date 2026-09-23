/**
 * GPCRome wheel — receptor picker modal (DataTables + Norges filters).
 * Expects globals: CreateColumnFilters, createDropdownFilters (NorgesDTFilterBuilder.js), jQuery, DataTables.
 */
(function ($) {
  'use strict';

  var pickerState = {
    initialized: false,
    dt: null,
    rowsById: {},
    assignSource: 'none'
  };

  var ASSIGN_SOURCE_LABELS = {
    none: 'None',
    class: 'Class',
    ligandtype: 'Chemotype',
    family: 'Receptor family'
  };

  // Sort the distinct group labels alphabetically and assign each a stable
  // sequential rank (1, 2, 3, ...) — shared by every page's onAdd so the
  // "Numbers" assignment logic isn't duplicated 5 times.
  window.mapperCoreBuildSequentialNumberMap = function (groupNameById) {
    var seen = {};
    var labels = [];
    Object.keys(groupNameById || {}).forEach(function (id) {
      var label = groupNameById[id];
      if (label != null && label !== '' && !seen[label]) {
        seen[label] = true;
        labels.push(label);
      }
    });
    labels.sort(function (a, b) { return String(a).localeCompare(String(b)); });
    var rank = {};
    labels.forEach(function (label, i) { rank[label] = i + 1; });
    var out = {};
    Object.keys(groupNameById || {}).forEach(function (id) {
      var label = groupNameById[id];
      out[id] = (label != null && label !== '') ? rank[label] : null;
    });
    return out;
  };

  function syncSelectAllCheckbox() {
    if (!pickerState.dt) {
      return;
    }
    var dt = pickerState.dt;
    var rows = dt.rows({ search: 'applied' });
    var checked = 0;
    rows.every(function () {
      var $cb = $(this.node()).find('.mapper-core-gpcrome-pick-cb');
      if ($cb.prop('checked')) {
        checked += 1;
      }
    });
    var total = rows.count();
    var $master = $('#mapper-core-gpcrome-picker-select-all');
    if (!total || checked === 0) {
      $master.prop({ checked: false, indeterminate: false });
    } else if (checked === total) {
      $master.prop({ checked: true, indeterminate: false });
    } else {
      $master.prop({ checked: false, indeterminate: true });
    }
  }

  function collectCheckedEntryIds() {
    if (!pickerState.dt) {
      return [];
    }
    var out = [];
    pickerState.dt.rows().every(function () {
      var $cb = $(this.node()).find('.mapper-core-gpcrome-pick-cb');
      if ($cb.prop('checked')) {
        var v = $cb.val();
        if (v != null && v !== '') {
          out.push(String(v));
        }
      }
    });
    return out;
  }

  function escapeAttr(val) {
    return String(val == null ? '' : val).replace(/&/g, '&amp;').replace(/"/g, '&quot;');
  }

  // Column data (e.g. GtoPdb names) can contain markup/entities (<sub>, &alpha;, ...).
  // Filtering still matches on the raw value, but the dropdown label should show the
  // decoded, human-readable text — same as the table cell itself renders it.
  var $htmlDecodeScratch = $('<div>');
  function decodeHtmlToText(str) {
    return $htmlDecodeScratch.html(str == null ? '' : String(str)).text();
  }

  function buildPickerFilterRow(dt) {
    // DataTables 2.x moves the original <thead> into a fixed scrollHead wrapper and
    // puts a visibility-hidden clone in the scrollBody. dt.column(0).header() always
    // returns the TH from the live visible header, so we navigate from there.
    var $filterCells = $(dt.column(0).header()).closest('thead').find('tr').last().find('th');
    var esc = $.fn.dataTable.util.escapeRegex;
    // Columns 2-6 get multi-select filters; 0-1 (checkbox, GPCRdb link) are left empty
    [2, 3, 4, 5, 6].forEach(function (colIdx) {
      var $th = $filterCells.eq(colIdx);
      var filterId = 'mapper-core-gpcrome-picker-table_Filter' + colIdx;
      var $sel = $('<select multiple="multiple" style="width:100%;">')
        .attr('id', filterId);

      var seen = {};
      var options = [];
      dt.column(colIdx).data().unique().each(function (d) {
        var raw = String(d == null ? '' : d).trim();
        if (!raw || seen[raw]) { return; }
        seen[raw] = true;
        options.push({ value: raw, label: decodeHtmlToText(raw) });
      });
      options.sort(function (a, b) { return a.label.localeCompare(b.label); });
      options.forEach(function (o) {
        $sel.append($('<option>').val(o.value).text(o.label));
      });

      $th.append($sel);

      (function (col) {
        $sel.on('change', function () {
          var vals = $(this).val() || [];
          if (!vals.length) {
            dt.column(col).search('').draw();
          } else {
            var regex = vals.map(function (v) { return '^' + esc(v) + '$'; }).join('|');
            dt.column(col).search(regex, { regex: true, smart: false }).draw();
          }
        });
      })(colIdx);

      $sel.select2({
        multiple: true,
        closeOnSelect: true,
        placeholder: { text: 'Filter' },
        dropdownAutoWidth: true,
        width: 'element',
        dropdownParent: $('#mapper-core-gpcrome-picker-modal')
      });
    });
  }

  function buildDataTable(rows) {
    var columns = [
      {
        data: 'id',
        name: '_sel',
        orderable: false,
        searchable: false,
        render: function (entryId, type /*, row*/) {
          if (type !== 'display') {
            return '';
          }
          var vid = escapeAttr(entryId == null ? '' : entryId);
          return (
            '<input type="checkbox" class="mapper-core-gpcrome-pick-cb" value="' +
            vid +
            '" aria-label="Select receptor">'
          );
        }
      },
      {
        data: 'gpcrdb_link',
        name: 'GPCRdb',
        orderable: false,
        searchable: false,
        defaultContent: '',
        render: function (link, type /*, row*/) {
          if (type !== 'display') {
            return '';
          }
          if (!link) {
            return '';
          }
          return (
            '<a href="' +
            escapeAttr(link) +
            '" target="_blank" rel="noopener" title="GPCRdb receptor page">' +
            '<img class="gpcrdb-link" src="/static/home/logo/gpcr/main.png" width="12" height="12" alt="GPCRdb"></a>'
          );
        }
      },
      {
        data: 'name_plain',
        name: 'GtoPdb',
        render: function (plain, type, row) {
          if (type === 'display') {
            if (row && row.name_html) {
              return row.name_html;
            }
            return plain || (row && row.id) || '';
          }
          return plain || '';
        }
      },
      { data: 'gene', name: 'Gene' },
      {
        data: 'uniprot',
        name: 'UniProt',
        render: function (d, type, row) {
          if (type !== 'display') {
            return d || '';
          }
          if (!d) {
            return '-';
          }
          var href =
            (row && row.uniprot_link) || ('https://www.uniprot.org/uniprot/' + d);
          return '<a href="' + escapeAttr(href) + '" target="_blank" rel="noopener">' + d + '</a>';
        }
      },
      { data: 'family', name: 'family' },
      { data: 'class', name: 'class' }
    ];

    var columnDefs = [
      { targets: '_all', className: 'dt-head-center dt-body-center dt-center' },
      { targets: [0], width: '38px', orderable: false, searchable: false },
      {
        targets: [1],
        width: '36px',
        orderable: false,
        searchable: false,
        className: 'dt-head-center dt-body-center dt-center'
      }
    ];

    pickerState.dt = $('#mapper-core-gpcrome-picker-table').DataTable({
      dom: "<'row'<'col-sm-12'tr>>" + "<'row'<'col-sm-12'i>>",
      data: rows || [],
      columns: columns,
      columnDefs: columnDefs,
      autoWidth: true,
      processing: false,
      deferRender: true,
      paging: false,
      scrollX: true,
      scrollY: '52vh',
      scrollCollapse: true,
      orderCellsTop: true,
      order: [[4, 'asc']]
    });
    var dt = pickerState.dt;
    buildPickerFilterRow(dt);

    $('#mapper-core-gpcrome-picker-table').on('draw.dt', syncSelectAllCheckbox);
    $('#mapper-core-gpcrome-picker-table tbody').on('change', '.mapper-core-gpcrome-pick-cb', syncSelectAllCheckbox);
    // Click anywhere on a row to toggle its checkbox
    $('#mapper-core-gpcrome-picker-table tbody').on('click', 'tr', function (e) {
      if ($(e.target).is('input[type="checkbox"], a, img')) return;
      var $cb = $(this).find('.mapper-core-gpcrome-pick-cb');
      $cb.prop('checked', !$cb.prop('checked')).trigger('change');
    });
  }

  window.mapperCoreInitGpcromePickerModal = function (opts) {
    opts = opts || {};
    if (!$('#mapper-core-gpcrome-picker-table').length) {
      return;
    }

    pickerState.pendingRows = $.isArray(opts.pickerRows) ? opts.pickerRows : [];
    pickerState.rowsById = {};
    pickerState.pendingRows.forEach(function (r) {
      if (r && r.id != null) {
        pickerState.rowsById[r.id] = r;
      }
    });

    $('#mapper-core-gpcrome-picker-modal').on('shown.bs.modal', function () {
      // Build lazily on first open (while the modal is actually visible) so DataTables/
      // select2 measure real widths instead of the 0-width hidden container at page load.
      if (!pickerState.initialized) {
        buildDataTable(pickerState.pendingRows);
        pickerState.initialized = true;
      }
      if (pickerState.dt) {
        pickerState.dt.columns.adjust().draw(false);
      }
      syncSelectAllCheckbox();
    });

    $('#mapper-core-gpcrome-picker-modal').on('hidden.bs.modal', function () {
      $('#mapper-core-gpcrome-picker-select-all').prop({ checked: false, indeterminate: false });
      if (pickerState.dt) {
        pickerState.dt.$('.mapper-core-gpcrome-pick-cb').prop('checked', false);
      }
    });

    $('#mapper-core-gpcrome-picker-modal').on(
      'change',
      '#mapper-core-gpcrome-picker-select-all',
      function () {
        if (!pickerState.dt) {
          return;
        }
        var on = $(this).prop('checked');
        pickerState.dt.rows({ search: 'applied' }).every(function () {
          $(this.node()).find('.mapper-core-gpcrome-pick-cb').prop('checked', on);
        });
        syncSelectAllCheckbox();
      }
    );

    // Bound directly on the dropdown-menu (not the modal) because a site-wide
    // convention calls e.stopPropagation() on every .dropdown-menu click (to keep
    // other dropdowns open on interaction), which would otherwise swallow this
    // click before it ever reaches a handler delegated from the modal ancestor.
    $('#mapper-core-gpcrome-picker-assign-menu').on('click', '.mapper-core-picker-assign-src', function () {
      var src = $(this).attr('data-src') || 'none';
      pickerState.assignSource = src;
      $(this).siblings('.mapper-core-picker-assign-src').addBack()
        .toggleClass('btn-outline-primary', true).toggleClass('btn-primary', false);
      $(this).toggleClass('btn-primary', true).toggleClass('btn-outline-primary', false);
      $(this).closest('.dropdown').find('.mapper-core-picker-assign-toggle-label')
        .text('Assign value: ' + (ASSIGN_SOURCE_LABELS[src] || 'None'));
    });

    $('#mapper-core-gpcrome-picker-add').on('click', function () {
      if (!pickerState.dt) {
        return;
      }
      var ids = collectCheckedEntryIds();
      if (!ids.length) {
        alert('Select at least one receptor (filtered rows — use tick boxes).');
        return;
      }
      if (typeof opts.onAdd === 'function') {
        var meta = null;
        if (pickerState.assignSource && pickerState.assignSource !== 'none') {
          var groupNameById = {};
          ids.forEach(function (id) {
            var row = pickerState.rowsById[id];
            groupNameById[id] = row ? (row[pickerState.assignSource] || '') : '';
          });
          meta = { groupNameById: groupNameById };
        }
        opts.onAdd(ids, meta);
      }
      $('#mapper-core-gpcrome-picker-modal').modal('hide');
    });
  };
})(jQuery);

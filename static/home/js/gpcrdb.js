/*eslint complexity: ["error", 15]*/
/*eslint wrap-iife: ["error", "outside"]*/
/*eslint quotes: ["error", "double", { "avoidEscape": true }]*/

function select_all(e) {
  let checkedStatus = $(e).prop("checked");

  $(".select-all  ").each(function() {
    $(this).prop("checked", checkedStatus);
  });

  $(".alt").each(function() {
    $(this).prop("checked", checkedStatus);
  });
}

/**
 * Default Bootstrap alert function that can be used on all GPCRdb pages
 */
function showAlert(message, alerttype) {
  // Alerttype: success, info, warning, danger
  // See https://getbootstrap.com/docs/3.3/components/#alerts
  $("#gpcrdb_alert_placeholder").append('<div id="gpcrdb_alert" class="alert alert-' +
    alerttype + '"><button type="button" class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
    message + "</div>");

  setTimeout(function() {
    $("#gpcrdb_alert").remove();
  }, 4000);
}

/**
 * This Array filter function can be used to deduplicate an Array
 */
function onlyUnique(value, index, self) {
  return self.indexOf(value) === index;
}

/**
 * This function exports a list in a separated fashion to clipboard.
 */
function copyListToClipboard(selected, delimiter = " ") {
  if (selected.length > 0) {
    // Deduplicate and sort list
    selected = selected.filter(onlyUnique).sort();

    // Place values in textbox and copy to clipboard
    var copybox = $('<input type="text">');
    $(document.body).append(copybox);
    copybox.val(selected.join(delimiter));
    copybox.select();
    copybox[0].setSelectionRange(0, 99999);
    document.execCommand("copy");
    showAlert("The selection has been copied to your clipboard.", "info");
    copybox.remove();
  } else {
    showAlert("There is no selection to copy to clipboard.", "warning");
  }
}

/**
 * This helper function can create a range of identical YADCF filters with a single call
 *
 * Arg list: createYADCFfilters(start_column, num_cols, filter_type, select_type*, filter_default_label*, filter_reset_button_text*, filter_match_mode*, column_data_type*, width*, html5_data*, ignore_char*, html_data_type*, html5_data*, other*)
 * Arguments:
 *   other: a javascript object with column parameters of Yet Another DataTables Column Filter - (yadcf) (jquery.dataTables.yadcf.js)
 * The asterisk indicates an optional function argument
 */
function createYADCFfilters(start_column, num_cols, filter_type, select_type = null, filter_default_label = "", filter_reset_button_text = false, filter_match_mode = null, column_data_type = null, width = null, html5_data = null, ignore_char = null, html_data_type = null, other = null) {
  
  let filters = [];
  for (let i = 0; i < num_cols; i++) {
    let filter = {
      "column_number": start_column + i,
      filter_type,
      filter_default_label,
      filter_reset_button_text
    };
    if (select_type !== null) {
      filter["select_type"] = select_type;
    }
    if (filter_match_mode !== null) {
      filter["filter_match_mode"] = filter_match_mode;
    }
    if (column_data_type !== null) {
      filter["column_data_type"] = column_data_type;
    }
    if (width !== null) {
      filter["select_type_options"] = {
        width
      };
    }
    if (html5_data !== null) {
      filter["html5_data"] = html5_data;
    }
    if (ignore_char !== null) {
      filter["ignore_char"] = ignore_char;
    }
    if (html_data_type !== null) {
      filter["html_data_type"] = html_data_type;
    }

    if (other !== null) {
      $.each( other, function( key, value ) {
        
        filter[key] = value;
      });
    }
    filters.push(filter);
  }
  
  return filters;
}

 


/**
 * Exports a specified table to an Excel file.
 * @param {string} tableId - The ID of the table element to be exported.
 * @param {string} sheetName - The name to be given to the Excel sheet.
 * @param {string} fileName - The name of the file to be created.
 */
function GlobalTableToExcel(tableId, sheetName, fileName) {
  // Get the table element by ID.
  const tableElement = document.getElementById(tableId);

  // If the table doesn't exist, log an error and exit the function.
  if (!tableElement) {
    console.error(`Table with ID "${tableId}" not found.`);
    return;
  }

  // Create a new Excel workbook.
  const wb = XLSX.utils.book_new();

  // Clean the table and prepare it for export.
  const clonedTable = cleanTable(tableElement.cloneNode(true));
  const exportTable = prepareExportTable(clonedTable);

  // Convert the prepared table to an Excel worksheet.
  const ws_data = XLSX.utils.table_to_sheet(exportTable);

  // Add the worksheet to the workbook and write it to a file.
  XLSX.utils.book_append_sheet(wb, ws_data, sheetName);
  XLSX.writeFile(wb, fileName, { bookSST: true });
}

/**
 * Cleans the table by removing unwanted elements, cleaning headers, and processing cells.
 * @param {HTMLElement} table - The table to be cleaned.
 * @return {HTMLElement} The cleaned table.
 */
function cleanTable(table) {
  // Remove unwanted elements, clean headers, and process cells.
  removeUnwantedElements(table);
  cleanTableHeaders(table);
  processTableCells(table);

  // Return the cleaned table.
  return table;
}

/**
 * Removes specific unwanted elements from the table.
 * @param {HTMLElement} table - The table to be cleaned.
 */
function removeUnwantedElements(table) {
  // Remove elements with class '.yadcf-filter-wrapper' and all buttons.
  table.querySelectorAll('.yadcf-filter-wrapper').forEach(element => element.closest('tr').remove());
  table.querySelectorAll('button').forEach(button => button.remove());
}

/**
 * Cleans the headers of the table.
 * @param {HTMLElement} table - The table with headers to be cleaned.
 */
function cleanTableHeaders(table) {
  // Process each header, remove 'span' elements and replace 'br' with space.
  table.querySelectorAll('th').forEach(th => {
    const clonedTh = th.cloneNode(true);
    clonedTh.querySelectorAll('span').forEach(span => span.remove());
    clonedTh.querySelectorAll('br').forEach(br => br.replaceWith(' '));
    th.textContent = clonedTh.textContent;
  });
}

/**
 * Processes each cell in the table.
 * @param {HTMLElement} table - The table with cells to be processed.
 */
function processTableCells(table) {
  // Iterate over each cell to remove comment nodes and process links and line breaks.
  table.querySelectorAll('td').forEach(td => {
    Array.from(td.childNodes).forEach(node => {
      if (node.nodeType === Node.COMMENT_NODE) {
        node.remove();
      }
    });
    handleImageLinks(td);
    replaceLineBreaks(td);
    replaceRefinedLinks(td);
  });
}

/**
 * Converts image links in a table cell to text.
 * @param {HTMLElement} td - The table cell to process.
 */
function handleImageLinks(td) {
  // If an image within a link is found, replace its content with the link's URL.
  const img = td.querySelector('img');
  if (img && img.parentNode.tagName === 'A') {
    td.textContent = img.parentNode.href;
  }
}

/**
 * Replaces line breaks in a table cell and wraps content in a paragraph tag.
 * @param {HTMLElement} td - The table cell to process.
 */
function replaceLineBreaks(td) {
  // Get specific HTML content, replace line breaks, and wrap in a paragraph tag.
  const htmlContent = td.querySelector('a i.glyphicon.glyphicon-file.simple-popover.gpcrdb-link')?.getAttribute('data-content');
  if (htmlContent) {
    let newTextContent = htmlContent.replace(/<br>/g, "_X_").replace(/<[^>]*>/g, '').replace(/_X_/g, '<br>' );
    td.innerHTML = `<p>${newTextContent}</p>`;
  }
}

/**
 * Converts refined links in a table cell to plain text.
 * @param {HTMLElement} td - The table cell to process.
 */
function replaceRefinedLinks(td) {
  // Replace the text of certain links with their URLs.
  td.querySelectorAll('a').forEach(a => {
    if (a.textContent.endsWith('_refined')) {
      a.textContent = a.href;
    }
  });
}

/**
 * Prepares the table for export, selecting specific rows or the entire table.
 * @param {HTMLElement} clonedTable - The cloned table to be prepared for export.
 * @return {HTMLElement} The table prepared for export.
 */
function prepareExportTable(clonedTable) {
  // Create a new table with selected rows or return the whole table.
  const selectedRows = clonedTable.querySelectorAll('tr.alt_selected');
  if (selectedRows.length > 0) {
    const exportTable = document.createElement('table');
    exportTable.appendChild(clonedTable.querySelector('thead').cloneNode(true));
    selectedRows.forEach(row => exportTable.appendChild(row.cloneNode(true)));
    return exportTable;
  } else {
    return clonedTable;
  }
}

/**
 * Ajax call for populating reference modal after clicking reference icons in the menu system
 * @param {string} keys multiple keys separated by |
 */
function getReference(keys, callback) {
  $.ajax({
    url: '/common/getreference/',
    data: { keys: keys.toUpperCase() },
    type: 'GET',
    async: true,
    dataType: 'json',
    success: function(data) {
      callback(data.html);  // pass HTML to the callback
    },
    error: function(xhr, status, error) {
      console.error("AJAX Error:", xhr.responseText, error);
      callback("");  // still call back with empty string on error
    }
  });
}

(function (global) {
  'use strict';

  /*
   * Simple generic DataTables → XLSX exporter
   * ------------------------------------------------------------
   * Signature:
   *   ExportDataTablesToXlsx(tableSelector, sheetName?, fileName?)
   *
   * Inputs:
   *   - tableSelector:  CSS selector or DOM node for a DataTable
   *   - sheetName:      optional worksheet name (default: 'Export')
   *   - fileName:       optional file name (default: `${sheetName}_filtered_YYYY-MM-DD.xlsx`)
   *
   * Behavior:
   *   - Detects THEAD structure assumed across your pages:
   *       [0] dummy, [1] super header (merged via colspans), [2] titles, [3] filters
   *   - Uses *currently visible* columns (robust visibility detection)
   *   - Exports filtered + ordered rows exactly as on screen
   *   - Cleans header cells from icons/buttons before reading titles
   *   - Uses title names from column defenition (when generating data)
   *   - Auto column widths; merges super-header ranges
   *   - Depends only on DataTables and SheetJS (xlsx.full.min.js)
   */

  /** Ensure jQuery DataTables and SheetJS are present. */
  function assertDeps() {
    if (!global.jQuery || !jQuery.fn || !jQuery.fn.dataTable) {
      throw new Error('DataTables not found. Load jQuery DataTables first.');
    }
    if (!(global.XLSX && XLSX.utils)) {
      throw new Error('SheetJS not found. Include xlsx.full.min.js first.');
    }
  }

  /**
   * Sanitize a worksheet name to satisfy Excel constraints.
   * - Removes invalid characters and trims to 31 chars.
   */
  function safeSheetName(name) {
    const cleaned = String(name || 'Export')
      .replace(/[:\\\/\?\*\[\]]/g, '–')
      .slice(0, 31)
      .trim();
    return cleaned || 'Export';
  }

  /**
   * Generate a default filename based on the sheet name and today’s date.
   * Example: "Export_filtered_2025-10-03.xlsx"
   */
  function defaultFileName(base) {
    const stamp = new Date().toISOString().slice(0, 10);
    const b = (base || 'Export').replace(/\s+/g, '_');
    return `${b}_filtered_${stamp}.xlsx`;
  }

  /**
   * Convert an HTML snippet to clean, single-line text.
   * - Replaces <br> with newlines
   * - Removes interactive/icons from headers
   * - Collapses whitespace
   */
  function htmlToText(html) {
    const div = document.createElement('div');
    div.innerHTML = String(html ?? '').replace(/<br\s*\/?>(?!\n)/gi, '\n');
    div.querySelectorAll('button, i, .glyphicon, .icon-button, .dt-column-order')
       .forEach(n => n.remove());
    return (div.textContent || div.innerText || '').replace(/\s+/g, ' ').trim();
  }

  /**
   * Best-effort numeric coercion:
   * - Keeps true numbers as numbers
   * - Trims and normalizes thousands/decimal separators (e.g., "1,234.56" or "1.234,56")
   * - Returns original text if it’s not a plain number
   */
  function coerceMaybeNumber(v) {
    if (v == null) return '';
    if (typeof v === 'number') return v;
    const s = String(v).trim();
    if (!s) return '';
    const norm = s
      .replace(/\s+/g, '')
      .replace(/(?<=\d),(?=\d{3}\b)/g, '')   // 1,234 -> 1234
      .replace(/,(\d{1,})$/, '.$1');         // 1,23 -> 1.23
    return (/^[+-]?\d*(?:\.\d+)?$/).test(norm) ? Number(norm) : s;
  }

  /**
   * Locate relevant THEAD rows according to the shared layout contract:
   * [0] dummy, [1] super header, [2] titles, [3] filters.
   */
  function getHeadRows(tableEl) {
    const thead = tableEl.tHead;
    const rows = thead ? Array.from(thead.rows) : [];
    return {
      rows,
      dummyRow: rows[0] || null,               // [0] dummy (often TD colspan)
      superRow: rows[1] || null,               // [1] super header (group labels with colspans)
      titleRow: rows[rows.length - 2] || null, // [-2] titles (leaf text for columns)
      filterRow: rows[rows.length - 1] || null // [-1] filters
    };
  }

  /**
   * Parse the super-header row into an ordered list of groups.
   * Uses the original HTML `colspan` attribute (stable even when columns are hidden)
   * rather than the live `th.colSpan`, which DataTables/DOM might mutate.
   * @returns Array<{label: string, span: number}>
   */
  function getSuperGroups(superRow) {
    const groups = [];
    if (!superRow) return groups;

    Array.from(superRow.children).forEach(th => {
      const label = htmlToText(th.innerHTML);
      const attrSpan = Number(th.getAttribute('colspan')) || 1; // prefer attribute
      const liveSpan = Number(th.colSpan) || attrSpan || 1;     // safe fallback
      const span = Math.max(1, attrSpan || liveSpan);
      groups.push({ label, span });
    });

    return groups;
  }

  /**
   * Determine the DataTables column indexes that are truly visible.
   * Combines DataTables API visibility with DOM checks to ignore hidden/invisible headers.
   */
  function buildVisibleCols(dt, settings) {
    const idxs = dt.columns().indexes().toArray();
    return idxs.filter(i => {
      try { if (dt.column(i).visible() === false) return false; } catch (e) {}
      const th = (settings.aoColumns[i] && settings.aoColumns[i].nTh) || dt.column(i).header();
      if (th && th instanceof HTMLElement) {
        const style = getComputedStyle(th);
        if (style.display === 'none' || style.visibility === 'hidden') return false;
        if (th.offsetParent === null && style.position !== 'sticky') return false;
      }
      return true;
    });
  }

  /**
   * Build the top “group row” and merge ranges by *packing* super groups
   * left-to-right across the current visible columns.
   * This yields contiguous merges (e.g., 5/5/5) regardless of which concrete
   * leaf columns are hidden, matching what users expect visually.
   *
   * @param groups       Array<{label, span}> from the super header
   * @param visibleCount Number of currently visible columns
   * @returns { groupRow: any[], merges: Array<{s:{r,c}, e:{r,c}}>}
   */
  function buildPackedGroupRow(groups, visibleCount) {
    const groupRow = Array(visibleCount).fill(null);
    const merges = [];

    let pos = 0;
    for (const g of groups) {
      if (pos >= visibleCount) break;
      const span = Math.min(g.span, visibleCount - pos);
      if (span <= 0) continue;

      if (g.label) groupRow[pos] = g.label;
      if (span > 1) {
        merges.push({ s: { r: 0, c: pos }, e: { r: 0, c: pos + span - 1 } });
      }
      pos += span;
    }

    return { groupRow, merges };
  }

  /**
   * Export the visible state of a DataTable to an XLSX file.
   * - Honors current search/filter + ordering
   * - Uses only visible columns
   * - Builds a super-header row with merged cells based on the super-row colspans
   */
  function ExportDataTablesToXlsx(tableSelector, sheetName, fileName) {
    assertDeps();

    const $ = jQuery;
    const dt = $(tableSelector).DataTable();
    const settings = dt.settings()[0];

    // Resolve table element safely
    const tableEl = (typeof tableSelector === 'string')
      ? document.querySelector(tableSelector)
      : tableSelector;
    if (!tableEl) throw new Error('Table element not found for selector: ' + tableSelector);

    // THEAD rows (we read group labels from superRow; titles from titleRow)
    const { superRow, titleRow } = getHeadRows(tableEl);
    if (!titleRow) throw new Error('Titles row not found in THEAD.');

    // Parse groups and current visible column order
    const groups = getSuperGroups(superRow);
    const visibleDTIdxs = buildVisibleCols(dt, settings);
    if (!visibleDTIdxs.length) { alert('No columns to export.'); return; }

    // Build the grouping row and merges sized to the number of visible columns
    const { groupRow, merges } = buildPackedGroupRow(groups, visibleDTIdxs.length);

    // Row indexes after current search/order
    const rowIdxs = dt.rows({ search: 'applied', order: 'applied' }).indexes().toArray();

    // Header titles for visible columns (strip icons etc.)
    const headerTitles = visibleDTIdxs.map(ci => {
      const colCfg = (settings.aoColumns && settings.aoColumns[ci]) || {};
      if (colCfg.sName) return String(colCfg.sName);
      const th = dt.column(ci).header() || colCfg.nTh;
      if (th && th.innerHTML != null) return htmlToText(th.innerHTML);
      if (titleRow && titleRow.children && ci < titleRow.children.length) {
        return htmlToText(titleRow.children[ci].innerHTML);
      }
      return `Column ${ci}`;
    });

    // Body values (prefer "export" render; fallback to "display"), with numeric coercion
    const body = rowIdxs.map(ri =>
      visibleDTIdxs.map(ci => {
        const exp = dt.cell(ri, ci).render('export');
        const val = (exp != null && exp !== '') ? exp : dt.cell(ri, ci).render('display');
        return coerceMaybeNumber(htmlToText(val));
      })
    );

    // Build worksheet (group row + header titles + data)
    const aoa = [groupRow, headerTitles, ...body];
    const ws = XLSX.utils.aoa_to_sheet(aoa);
    if (merges.length) ws['!merges'] = merges;

    // Auto column widths based on header + data contents
    ws['!cols'] = headerTitles.map((h, j) => {
      const maxLen = Math.max(String(h || '').length, ...body.map(r => String(r[j] || '').length));
      return { wch: Math.min(Math.max(10, maxLen + 2), 60) };
    });

    // Workbook + file write
    const sheet = safeSheetName(sheetName || 'Export');
    const wb = XLSX.utils.book_new();
    XLSX.utils.book_append_sheet(wb, ws, sheet);

    const fname = (typeof fileName === 'string' && fileName) ? fileName : defaultFileName(sheet);
    XLSX.writeFile(wb, fname);
  }

  // expose
  global.ExportDataTablesToXlsx = ExportDataTablesToXlsx;

})(window);
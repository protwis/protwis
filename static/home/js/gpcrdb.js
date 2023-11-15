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
 * Arg list: createYADCFfilters(start_column, num_cols, filter_type, select_type*, filter_default_label*, filter_reset_button_text*, filter_match_mode*, column_data_type*, width*, html5_data*)
 * The asterisk indicates an optional function argument
 */
function createYADCFfilters(start_column, num_cols, filter_type, select_type = null, filter_default_label = "", filter_reset_button_text = false, filter_match_mode = null, column_data_type = null, width = null, html5_data = null, ignore_char = null) {
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

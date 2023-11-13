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
 * GlobalTableToExcel: Exports an HTML table to an Excel file.
 * 
 * @param {string} table - The ID of the HTML table element to be exported.
 * @param {string} [name='SheetJS'] - The name of the Excel sheet. Defaults to 'SheetJS'.
 * @param {string} [filename='out.xlsx'] - The name of the exported Excel file. Defaults to 'out.xlsx'.
 */

function GlobalTableToExcel(table, name, filename) {
  var wb = XLSX.utils.book_new();
  var ws_name = name || 'SheetJS';
  
  // Retrieve the table element using its ID
  var tableElement = document.getElementById(table);
  
  // Create a deep clone of the table to manipulate without affecting the original
  var clonedTable = tableElement.cloneNode(true);

  // Remove elements related to the 'yadcf' filter and any buttons
  $(clonedTable).find('.yadcf-filter-wrapper').closest('tr').remove();
  $(clonedTable).find('button').remove();
  
  // Iterate through each cell to find images within hyperlinks
  // If found, replace the cell's content with the hyperlink
  $(clonedTable).find('td').each(function() {
    var img = $(this).find('img');
    if (img.length > 0 && img[0].parentNode.tagName === 'A') {
      $(this).text(img[0].parentNode.href);
    }
  });
  
  // Check for rows that are selected (with class 'alt_selected')
  var selectedRows = $(clonedTable).find('tr.alt_selected');
  var exportTable;
  if(selectedRows.length > 0) {
      // If there are selected rows, create a new table for exporting
      // This new table includes the header from the original and the selected rows
      exportTable = document.createElement('table');
      $(exportTable).append($(clonedTable).find('thead').clone());
      $(exportTable).append(selectedRows);
  } else {
      // If no rows are selected, use the entire cloned table for exporting
      exportTable = clonedTable;
  }
  // Convert the table (whether it's the full table or just selected rows) to a worksheet
  var ws_data = XLSX.utils.table_to_sheet(exportTable);
  
  // Add the worksheet to the workbook with the specified sheet name
  XLSX.utils.book_append_sheet(wb, ws_data, ws_name);
  
  // Write the workbook and trigger the download of the Excel file
  XLSX.writeFile(wb, filename || 'out.xlsx');
}

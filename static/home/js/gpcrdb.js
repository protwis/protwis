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
 * This helper function translates a table into an excel sheets using an external called template
 *
 * Arg list: GlobalTableToExcel(table, name, filename)
 */
function GlobalTableToExcel(table, name, filename) {
    var uri = "data:application/vnd.ms-excel;base64,",
    template = "<html xmlns:o='urn:schemas-microsoft-com:office:office' xmlns:x='urn:schemas-microsoft-com:office:excel' xmlns='http://www.w3.org/TR/REC-html40'><head><!--[if gte mso 9]><xml><x:ExcelWorkbook><x:ExcelWorksheets><x:ExcelWorksheet><x:Name>{worksheet}</x:Name><x:WorksheetOptions><x:DisplayGridlines/></x:WorksheetOptions></x:ExcelWorksheet></x:ExcelWorksheets></x:ExcelWorkbook></xml><![endif]--></head><body><table>{table}</table></body></html>",
    base64 = function (s) {
        return window.btoa(unescape(encodeURIComponent(s)));
        }, format = function (s, c) {
        return s.replace(/{(\w+)}/g, function (m, p) {
            return c[p];
            });
        };
    table= $("#"+table).clone();
    $("#excel_table").html(table);
    // Clean up table to remove yadcf stuff
    $("#excel_table thead tr").css("height","");
    $("#excel_table thead th").css("height","");
    $("#excel_table thead div").css("height","");
    $("#excel_table thead .yadcf-filter-wrapper").remove();
    $("#excel_table thead button").remove();
    var tr = $("#excel_table thead tr:eq(1)");
    // reattach th titles
    tr.find("th").each (function( column, th) {
        if ($(th).attr("title")) {
            $(th).html($(th).attr("title"))
        }
    });

    var ctx = {
        worksheet: name || "Worksheet",
        table: $("#excel_table").html()
    };
    $("#excel_table").html("");
    document.getElementById("dlink").href = uri + base64(format(template, ctx));
    document.getElementById("dlink").download = filename;
    document.getElementById("dlink").click();
}

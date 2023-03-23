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
function createYADCFfilters(start_column, num_cols, filter_type, select_type = null, filter_default_label = "", filter_reset_button_text = false, filter_match_mode = null, column_data_type = null, width = null, html5_data=null) {
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

    filters.push(filter);
  }
  return filters;
}

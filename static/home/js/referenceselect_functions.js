/*eslint complexity: ["error", 8]*/
/*eslint quotes: ["error", "double", { "avoidEscape": true }]*/
/*global showAlert */
/*global yadcf*/

var referenceTable;
var selectedTargets = new Set();

/**
 * This function mains the shown and hidden selected targets
 * and updates an information message after each update/filter.
 */
var previousTargetCount = 0;
function updateTargetCount(){
  // Counting the selected targets matching the current filters
  var numTargets = $("table#uniprot_selection tbody input:checked").length;

  var message = selectedTargets.size.toString();
  if (numTargets === 1){
    message += " receptor selected";
  } else {
    message += " receptors selected";
  }

  var filtered = selectedTargets.size - numTargets;
  if (filtered > 0) {
    message += " ("+filtered+" currently filtered)";
  }

  $("#selection_table_info").html(message);
  if (previousTargetCount !== selectedTargets.size) {
    if (!$("#selection_table_info").is(":animated")) {
      $("#selection_table_info").effect("highlight", {color:"#FFAAAA"}, 1000 );
    }
    previousTargetCount = selectedTargets.size;
  }
}

/**
 * This function adds the target of the checkbox to the selected targets
 * @param {object} checkbox Checkbox element of the target to select
 */
function addTarget(checkbox){
    var slug = $(checkbox).attr("entry-value");
    $(checkbox).prop("checked", true);
    $(checkbox).closest("tr").addClass("selected");
    selectedTargets.add(slug);

}

/**
 * This function remove the target of the checkbox from the selected targets
 * @param {object} checkbox Checkbox element of the target to remove
 */
function removeTarget(checkbox){
  var slug = $(checkbox).attr("entry-value");
  $(checkbox).prop("checked", false);
  $(checkbox).closest("tr").removeClass("selected");

  selectedTargets.delete(slug);
}

/**
 * This function clear all filters on a the YADCF target table
 */
function clearFilters(){
  // Redrawing is slow: only reset filters when actually active
  if ($("#uniprot_selection_info").text().includes("filtered")){
    yadcf.exResetAllFilters(referenceTable);

    // Clear status "only_selected" buttons if used
    $("#only_selected").text("Only selected");
  }
}

/**
 * This function clear all currently selected targets and resets the filters
 * @returns {boolean} false to prevent event propagation
 */
function clearTargetSelection(){
  // clear filters, otherwise there will be mismatches
  // clearFilters();

  // uncheck all selected targets
  $("table#uniprot_selection tbody tr [type=checkbox]:checked").each(function() {
    removeTarget(this);
  });

  // update information message
  updateTargetCount();
  return false;
}

/**
 * This is a custom YADCF function that checks if a row is selected or not
 * In this case this is done using the checkbox in the first column
 * @param {object} filterVal Value to filter on (not applicable)
 * @param {object} columnVal Element in the filtered column
 * @param {object} rowValues All elements in this row (not used)
 * @param {object} stateVal Current DOM state of the row (not sufficient in this case)
 * @returns {boolean} true if row contains selected target otherwise false
 */
function selectedTargetFilter(filterVal, columnVal, rowValues, stateVal){
  var checkboxID = columnVal.match(/id="(.*?)"/g);
  if (checkboxID.length > 0){
      checkboxID = checkboxID[0].replace(/id="/g, "").replace('"', "");
      return $("#"+checkboxID).prop("checked");
  } else {
    return false;
  }
}

/**
 * This is a custom YADCF function that checks if a row is selected or not
 * In this case this is done using the checkbox in the first column
 * @param {object} filterVal Value to filter on (not applicable)
 * @param {object} columnVal Element in the filtered column
 * @param {object} rowValues All elements in this row (not used)
 * @param {object} stateVal Current DOM state of the row (not sufficient in this case)
 * @returns {boolean} true if row contains selected target otherwise false
 */
function checkAllTargets(){
  var changedTargetBoxes = 0;
  $("table#uniprot_selection tbody tr").each(function() {
    if (!$(this).hasClass("selected")){
      addTarget($(this).find("[type=checkbox]")[0]);
      changedTargetBoxes++;
    }
  });

  if (changedTargetBoxes === 0){
    $("table#uniprot_selection tbody tr").each(function() {
      removeTarget($(this).find("[type=checkbox]")[0]);
    });
  }

  updateTargetCount();
  return false;
}

/**
 * This function enables protein family selections from the tree browser
 * @param {string} slug Protein family slug of the targets that should be selected
 * @returns {boolean} false to prevent event propagation
 */
function selectInTable(slug){
  // clear filters, otherwise there will be mismatches
  clearFilters();

  $("table#uniprot_selection tbody tr [type=checkbox]").each(function() {
    if ($(this).attr("id") && $(this).attr("id").startsWith(slug)) {
      addTarget(this);
    }
  });

  updateTargetCount();
  return false;
}

/**
 * This function imports a user provided list of targets from an input textbox
 * parses this list and tries to select the matching targets
 */
function importTargets(){
  // clear filters, otherwise there will be mismatches
  clearFilters();

  // process the input table
  var messageMultiple = "";
  var msgTypeMultiple = "info";
  var inputEntries = $("#copyboxTargets").val();
  var splitEntries = inputEntries.split(/[ ,:;]+/);
  if (splitEntries.length > 1){
    messageMultiple = "Only <b>one receptor can be imported</b>. Please provide a single uniprot name to be imported.";
    msgTypeMultiple = "warning";
    showAlert(messageMultiple, msgTypeMultiple);
  }

  // Keep track of matches and misses
  var notFound = [];
  var parsed = 0;
  for (var i = 0; i < splitEntries.length; i++) {
    splitEntries[parseInt(i, 10)] = splitEntries[parseInt(i, 10)].trim().toLowerCase();
    splitEntries[parseInt(i, 10)] = splitEntries[parseInt(i, 10)].split("_")[0];

    // Check minimum protein name length
    if (splitEntries[parseInt(i, 10)].length >= 2){
      // Find checkbox with correct entry
      var items = $("table#uniprot_selection").find(`input[data-entry='${splitEntries[parseInt(i, 10)]}']`);
      if (items.length > 0){
        parsed++;
        addTarget(items[0]);
      } else {
        notFound.push(splitEntries[parseInt(i, 10)]);
      }
    }
  }

  // Add summary on message
  var message = "";
  var msgType = "info";
  if (parsed > 0){
    message = "<b>Successfully</b> imported "+parsed+" targets.<br>";
    if (notFound.length > 0){
      message += "<br>The following name(s) could <i>not</i> be matched:<br>&#8226;&nbsp;&nbsp;" + notFound.join("<br>&#8226;&nbsp;&nbsp;");
    }
  } else {
    message = "The target selection import was <b>not successful</b>. Please make sure you are using the uniprot target names.";
    msgType = "warning";
  }
  showAlert(message, msgType);
  updateTargetCount();
}

/**
 * This function exports the current taret selection to a textbox and copies the
 * selection to clipboard.
 */
function exportTargets(){
  var selected = [];
  var entries = $("table#uniprot_selection tbody tr [type=checkbox]:checked");
  if (entries.length > 0){
    for (var i = 0; i < entries.length; i++) {
      var checkbox = $(entries[parseInt(i, 10)]);
      selected.push(checkbox.attr("data-entry"));
    }
    // Place in textbox and copy to clipboard
    var copybox = $("#copyboxTargets");
    copybox.val(selected.join(" "));
    copybox.select();
    copybox[0].setSelectionRange(0, 99999);
    document.execCommand("copy");
    showAlert("The selected receptor has been copied to your clipboard.", "info");
  } else {
    showAlert("There are currenlty no selected targets.", "warning");
  }
}

/**
 * This function submits the selected targets to the backend and moves to the
 * next page. This function should be executed upon pressing the green button
 * after target selection.
 * @param {string} url The url to go to after synchronizing the target selection
 */
function submitSelection(url, minimum = 1) {
  if (selectedTargets.size >= minimum) {
    // set CSRF csrf_token
    $.ajaxSetup({
        headers:
        { "X-CSRF-TOKEN": $("meta[name=\"csrf-token\"]").attr("content") }
    });

    // Submit proteins to target selection
    var group = Array.from(selectedTargets);
    $.post("/common/referenceformread", { "input-targets": group.join("\r") },  function (data) {
      // On success go to url page
      window.location.href = url;
    })
    .fail(
      function(){
        showAlert("Something went wrong, please try again or contact us.", "danger");
      });
  } else {
      showAlert("Please select a target.", "warning");
    }
}

/**
 * This function initializes the YADCF datatables for a specific element
 * @param {string} elementID The identifier of the container containing the table
 */
function initTargetTable(elementID, pathway) {

    if (pathway){
      if (!$.fn.DataTable.isDataTable(elementID + " table")) {
          referenceTable = $(elementID + " table").DataTable({
              deferRender: true,
              scrollY: "50vh",
              scrollX: true,
              scrollCollapse: true,
              scroller: true,
              paging: false,
              aaSorting: [],
              autoWidth: false,
              bInfo: true,
              columnDefs: [{
                  targets: 0,
                  orderable: false,
                  className: "select-checkbox"
              },{
                  targets: 1,
                  className: "text-center"
              },
            ],
          });

          yadcf.init(referenceTable,
              [
                  {
                      column_number: 0,
                      filter_type: "custom_func",
                      custom_func: selectedTargetFilter,
                      filter_container_id: "hidden_filter_container",
                      html5_data: "data-search", // which does not exist - prevent warning logs
                  },
                  {
                      column_number: 1,
                      filter_type: "multi_select",
                      select_type: "select2",
                      filter_default_label: "Class",
                      filter_reset_button_text: false,
                      select_type_options: {
                          "width": "70px",
                      },
                  },
                  {
                      column_number: 2,
                      filter_type: "multi_select",
                      select_type: "select2",
                      filter_default_label: "Ligand",
                      filter_reset_button_text: false,
                      filter_match_mode : "exact",
                  },
                  {
                      column_number: 3,
                      filter_type: "multi_select",
                      select_type: "select2",
                      filter_default_label: "Family",
                      filter_reset_button_text: false,
                      filter_match_mode : "exact",
                  },
                  {
                      column_number: 4,
                      filter_type: "multi_select",
                      select_type: "select2",
                      filter_default_label: "Species",
                      filter_reset_button_text: false,
                      filter_match_mode : "exact",
                      select_type_options: {
                          "width": "110px",
                      }
                  },
                  {
                      column_number: 5,
                      filter_type: "multi_select",
                      select_type: "select2",
                      column_data_type: "html",
                      filter_default_label: "Uniprot",
                      filter_reset_button_text: false,
                      select_type_options: {
                          "width": "110px",
                      }
                  },
                  {
                      column_number: 6,
                      filter_type: "multi_select",
                      select_type: "select2",
                      column_data_type: "html",
                      html_data_type: "text",
                      filter_default_label: "GtP",
                      filter_match_mode : "exact",
                      filter_reset_button_text: false,
                      select_type_options: {
                          "width": "110px",
                      }
                  },
                  {
                      column_number: 7,
                      filter_type: "range_number",
                      filter_default_label: ["From", "To"],
                      filter_reset_button_text: false,
                      //style_class: "center"
                      //html5_data: "data-search",
                      column_data_type: "html",
                  },
              ], {
                  // cumulative_filtering: true,
                  filters_tr_index: 1
              }
          );
          yadcf.exFilterColumn(referenceTable, [[4, ["Homo sapiens"]]], true);
      }
    }else{
      if (!$.fn.DataTable.isDataTable(elementID + " table")) {
          referenceTable = $(elementID + " table").DataTable({
  //            dom: "ftip",
              deferRender: true,
              scrollY: "50vh",
              scrollX: true,
              scrollCollapse: true,
              scroller: true,
              paging: false,
              aaSorting: [],
              autoWidth: false,
              bInfo: true,
              columnDefs: [{
                  targets: 0,
                  orderable: false,
                  className: "select-checkbox"
              },{
                  targets: 1,
                  className: "text-center"
              },
            ],
          });

          yadcf.init(referenceTable,
              [
                  {
                      column_number: 0,
                      filter_type: "custom_func",
                      custom_func: selectedTargetFilter,
                      filter_container_id: "hidden_filter_container",
                      html5_data: "data-search", // which does not exist - prevent warning logs
                  },
                  {
                      column_number: 1,
                      filter_type: "multi_select",
                      select_type: "select2",
                      filter_default_label: "Class",
                      filter_reset_button_text: false,
                      select_type_options: {
                          "width": "70px",
                      },
                  },
                  {
                      column_number: 2,
                      filter_type: "multi_select",
                      select_type: "select2",
                      filter_default_label: "Ligand",
                      filter_reset_button_text: false,
                      filter_match_mode : "exact",
                  },
                  {
                      column_number: 3,
                      filter_type: "multi_select",
                      select_type: "select2",
                      filter_default_label: "Family",
                      filter_reset_button_text: false,
                      filter_match_mode : "exact",
                  },
                  {
                      column_number: 4,
                      filter_type: "multi_select",
                      select_type: "select2",
                      filter_default_label: "Species",
                      // filter_delay: 10,
                      filter_reset_button_text: false,
                      filter_match_mode : "exact",
                  },
                  {
                      column_number: 5,
                      filter_type: "multi_select",
                      select_type: "select2",
                      column_data_type: "html",
                      filter_default_label: "Uniprot",
                      filter_reset_button_text: false,
                      select_type_options: {
                          "width": "110px",
                      }
                  },
                  {
                      column_number: 6,
                      filter_type: "multi_select",
                      select_type: "select2",
                      column_data_type: "html",
                      html_data_type: "text",
                      filter_default_label: "GtP",
                      filter_match_mode : "exact",
                      filter_reset_button_text: false,
                      select_type_options: {
                          "width": "110px",
                      }
                  },
                  {
                      column_number: 7,
                      filter_type: "range_number",
                      filter_default_label: ["From", "To"],
                      filter_reset_button_text: false,
                      column_data_type: "html",
                  },
                  {
                      column_number: 8,
                      filter_type: "range_number",
                      filter_default_label: ["From", "To"],
                      filter_reset_button_text: false,
                      // style_class: "center",
                      // column_data_type: "html",
                      html5_data: "data-search",
                  },
                  {
                      column_number: 9,
                      filter_type: "range_number",
                      filter_default_label: ["From", "To"],
                      filter_reset_button_text: false,
                      // style_class: "center",
                      // column_data_type: "html",
                      html5_data: "data-search",
                  },
                  {
                      column_number: 10,
                      filter_type: "range_number",
                      filter_default_label: ["From", "To"],
                      filter_reset_button_text: false,
                      // style_class: "center",
                      // column_data_type: "html",
                      html5_data: "data-search",
                  },
                  /*{
                      column_number: 7,
                      filter_type: "multi_select",
                      select_type: "select2",
                      filter_default_label: "Approved",
                  },
                  {
                      column_number: 8,
                      filter_type: "multi_select",
                      select_type: "select2",
                      filter_default_label: "Clinical trial",
                  },*/
                  // {
                  //     column_number: 9,
                  //     filter_type: "multi_select",
                  //     select_type: "select2",
                  //     filter_default_label: "Gs",
                  //     filter_reset_button_text: false,
                  // },
                  // {
                  //     column_number: 10,
                  //     filter_type: "multi_select",
                  //     select_type: "select2",
                  //     filter_default_label: "Gi/o",
                  //     filter_reset_button_text: false,
                  // },
                  // {
                  //     column_number: 11,
                  //     filter_type: "multi_select",
                  //     select_type: "select2",
                  //     filter_default_label: "Gq/11",
                  //     filter_reset_button_text: false,
                  // },
                  // {
                  //     column_number: 12,
                  //     filter_type: "multi_select",
                  //     select_type: "select2",
                  //     filter_default_label: "G12/13",
                  //     filter_reset_button_text: false,
                  // },

              ], {
                  // cumulative_filtering: true,
                  filters_tr_index: 1
              }
          );
          yadcf.exFilterColumn(referenceTable, [[4, ["Homo sapiens"]]], true);
      }
    }
    // When redrawing update the information selection message
    referenceTable.on("draw.dt", function(e, oSettings) {
        updateTargetCount();
    });

    // Put top scroller
    // https://stackoverflow.com/questions/35147038/how-to-place-the-datatables-horizontal-scrollbar-on-top-of-the-table
    $(".dataTables_scrollHead").css({
        "overflow-x": "scroll"
    }).on("scroll", function(e){
        var scrollBody = $(this).parent().find(".dataTables_scrollBody").get(0);
        scrollBody.scrollLeft = this.scrollLeft;
        $(scrollBody).trigger("scroll");
    });

    // Ready to draw the table
    referenceTable.draw();

}


/**
 * This function initializes the YADCF datatables for a specific element
 * @param {string} elementID The identifier of the container containing the table
 */
function initLigandCountTable(elementID) {

    if (!$.fn.DataTable.isDataTable(elementID + " table")) {
        referenceTable = $(elementID + " table").DataTable({
            deferRender: true,
            scrollY: "50vh",
            scrollX: true,
            scrollCollapse: true,
            scroller: true,
            paging: false,
            aaSorting: [],
            bStateSave: true,
            autoWidth: false,
            bInfo: true,
            columnDefs: [{
                targets: 0,
                orderable: false,
                className: "select-checkbox"
            },{
                targets: 1,
                className: "text-center"
            },
          ],
        });

        yadcf.init(referenceTable,
            [
                {
                    column_number: 0,
                    filter_type: "custom_func",
                    custom_func: selectedTargetFilter,
                    filter_container_id: "hidden_filter_container",
                    html5_data: "data-search", // which does not exist - prevent warning logs
                },
                {
                    column_number: 1,
                    filter_type: "multi_select",
                    select_type: "select2",
                    filter_default_label: "Class",
                    filter_reset_button_text: false,
                    select_type_options: {
                        "width": "70px",
                    },
                },
                {
                    column_number: 2,
                    filter_type: "multi_select",
                    select_type: "select2",
                    filter_default_label: "Ligand",
                    filter_reset_button_text: false,
                    filter_match_mode : "exact",
                },
                {
                    column_number: 3,
                    filter_type: "multi_select",
                    select_type: "select2",
                    filter_default_label: "Family",
                    filter_reset_button_text: false,
                    filter_match_mode : "exact",
                },
                {
                    column_number: 4,
                    filter_type: "multi_select",
                    select_type: "select2",
                    filter_default_label: "Species",
                    filter_reset_button_text: false,
                    filter_match_mode : "exact",
                },
                {
                    column_number: 5,
                    filter_type: "multi_select",
                    select_type: "select2",
                    column_data_type: "html",
                    filter_default_label: "Uniprot",
                    filter_reset_button_text: false,
                    select_type_options: {
                        "width": "110px",
                    }
                },
                {
                    column_number: 6,
                    filter_type: "multi_select",
                    select_type: "select2",
                    column_data_type: "html",
                    html_data_type: "text",
                    filter_default_label: "GtP",
                    filter_match_mode : "exact",
                    filter_reset_button_text: false,
                    select_type_options: {
                        "width": "110px",
                    }
                },
                {
                    column_number: 7,
                    filter_type: "range_number",
                    filter_default_label: ["From", "To"],
                    filter_reset_button_text: false,
                    //style_class: "center"
                    //html5_data: "data-search",
                    column_data_type: "html",
                },
                {
                    column_number: 8,
                    filter_type: "range_number",
                    filter_default_label: ["From", "To"],
                    filter_reset_button_text: false,
                    //style_class: "center"
                    //html5_data: "data-search",
                    column_data_type: "html",
                },
                {
                    column_number: 9,
                    filter_type: "range_number",
                    filter_default_label: ["From", "To"],
                    filter_reset_button_text: false,
                    //style_class: "center"
                    //html5_data: "data-search",
                    column_data_type: "html",
                },
            ], {
                // cumulative_filtering: true,
                filters_tr_index: 1
            }
        );
        yadcf.exFilterColumn(referenceTable, [[4, ["Homo sapiens"]]], true);
    }
    // When redrawing update the information selection message
    referenceTable.on("draw.dt", function(e, oSettings) {
        updateTargetCount();
    });

    // Put top scroller
    // https://stackoverflow.com/questions/35147038/how-to-place-the-datatables-horizontal-scrollbar-on-top-of-the-table
    $(".dataTables_scrollHead").css({
        "overflow-x": "scroll"
    }).on("scroll", function(e){
        var scrollBody = $(this).parent().find(".dataTables_scrollBody").get(0);
        scrollBody.scrollLeft = this.scrollLeft;
        $(scrollBody).trigger("scroll");
    });

    // Ready to draw the table
    referenceTable.draw();

}
/**
 * This function submits the selected target slug to the backend and opens
 * the ligands overview page in a new window
* @param {string} slug The slug of the target for which to show the ligands
 */
/*function showLigands(slug) {
  // set CSRF csrf_token
  $.ajaxSetup({
      headers:
      { "X-CSRF-TOKEN": $("meta[name=\"csrf-token\"]").attr("content") }
  });

  // Submit the slug to target selection
  var group = [slug];
  $.post("/common/targetformread", { "input-targets": group.join("\r") },  function (data) {
    // On success go to ligand page
    window.open("/ligand/targets", "_blank");
  })
  .fail(
    function(){
      showAlert("Something went wrong, please try again or contact us.", "danger");
    });
}*/

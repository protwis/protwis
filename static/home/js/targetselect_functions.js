/*eslint complexity: ["error", 8]*/
/*eslint quotes: ["error", "double", { "avoidEscape": true }]*/
/*global showAlert */

var targetTable;
var selected_targets = new Set();
var imported_targets = new Set();

/**
 * This function mains the shown and hidden selected targets
 * and updates an information message after each update/filter.
 */
var previous_target_count = 0;
function updateTargetCount(){
  // Counting the selected targets matching the current filters
  var numTargets = $("table#uniprot_selection tbody input:checked").length;

  var message = (selected_targets.size + imported_targets.size).toString();

  if (numTargets === 1){
    message += " target selected";
  } else {
    message += " targets selected";
  }

  var filtered = selected_targets.size - numTargets;
  if (filtered > 0) {
    message += " ("+filtered+" currently filtered)";
  }

  $("#selection_table_info").html(message);
  if (previous_target_count !== selected_targets.size) {
    if (!$("#selection_table_info").is(":animated")) {
      $("#selection_table_info").effect("highlight", {color:"#FFAAAA"}, 1000 );
    }
    previous_target_count = selected_targets.size;
  }
}

/**
 * This function adds the target of the checkbox to the selected targets
 * @param {object} checkbox Checkbox element of the target to select
 */
function addTarget(checkbox, add_to_targets=true){
  var species = $("div#filters-species a.btn.active")[0].innerText.trim();

  if (add_to_targets) {
    if (species === "Human" && $(checkbox).attr("data-human")==="No"){
      showAlert("Note that <b>"+$(checkbox).attr("data-entry")+"</b> does not have a human ortholog.<br>Adjust the species selection to include it.", "warning");
    } else {
      var slug = $(checkbox).attr("id");
      $(checkbox).prop("checked", true);
      $(checkbox).closest("tr").addClass("selected");
      selected_targets.add(slug);
    }
  } else {
    var slug = $(checkbox).attr("id");
    $(checkbox).prop("checked", true);
    $(checkbox).closest("tr").addClass("selected");
  }

}

/**
 * This function remove the target of the checkbox from the selected targets
 * @param {object} checkbox Checkbox element of the target to remove
 */
function removeTarget(checkbox){
  var slug = $(checkbox).attr("id");
  $(checkbox).prop("checked", false);
  $(checkbox).closest("tr").removeClass("selected");

  selected_targets.delete(slug);
}

/**
 * This function clear all filters on a the YADCF target table
 */
function clearFilters(){
  // Redrawing is slow: only reset filters when actually active
  if ($("#uniprot_selection_info").text().includes("filtered")){
    yadcf.exResetAllFilters(targetTable);

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
  clearFilters();

  // uncheck all selected targets
  $("table#uniprot_selection tbody tr [type=checkbox]:checked").each(function() {
    removeTarget(this);
  });

  // Clean all (handling automated population by browser when using back button)
  selected_targets = new Set();
  imported_targets = new Set();

  // update information message
  updateTargetCount();
  return false;
}

/**
 * This function will filter all unselected targets from the table
 * and toggle the button to Show all targets (+ vice versa)
 * @param {object} target The button element that triggers this function
 */
function onlySelectedTargets(target){
  var msg1 = "Only selected";
  var msg2 = "Show all";

  if (target.textContent === "Only selected"){
    yadcf.exFilterColumn(targetTable, [[0]]);
    target.textContent = msg2;
  } else {
    // clear filters + show allTargets
    yadcf.exResetAllFilters(targetTable);
    target.textContent = msg1;
    clearFilters();
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
function selectedTargetFilter(filterVal, columnVal, rowValues, stateVal){
  var checkbox_id = columnVal.match(/id="(.*?)"/g);
  if (checkbox_id.length > 0){
      checkbox_id = checkbox_id[0].replace(/id="/g, "").replace('"', "");
      return $("#"+checkbox_id).prop("checked");
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
function check_all_targets(){
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
  var input_entries = $("#copyboxTargets").val();
  var split_entries = input_entries.split(/[ ,:;]+/);

  // Keep track of matches and misses
  var not_found = [];
  var full_names = [];
  var parsed = 0;
  var selected_items = [];
  for (var i = 0; i < split_entries.length; i++) {
    split_entries[parseInt(i, 10)] = split_entries[parseInt(i, 10)].trim().toLowerCase();
    full_names.push(split_entries[parseInt(i, 10)]);
    split_entries[parseInt(i, 10)] = split_entries[parseInt(i, 10)].split("_")[0];

    // Check minimum protein name length
    if (split_entries[parseInt(i, 10)].length >= 2){
      // Find checkbox with correct entry
      var items = $("table#uniprot_selection").find(`input[data-entry='${split_entries[parseInt(i, 10)]}']`);
      if (items.length > 0){
        addTarget(items[0], false);
        selected_items.push(items[0]);
      } else {
        not_found.push(split_entries[parseInt(i, 10)]);
      }
    }
  }
  // AJAX to add imported to selection
  var out = importTargetSelection(full_names.join(','))
  var slugs = out['slugs'];
  var species = out['species'];
  var entries = out['found_entries'];
  parsed = out['found_entries'].length;
  for (var i=0; i<entries.length;i++) {
    if (!imported_targets.has(entries[i])) {
      imported_targets.add(entries[i]);
    }
  }
  // Toggle species filter
  if (species.length>1 || species[0]!==1) {
    $("#filters-species .active").removeClass('active');
    $("a[data-target='#SpeciesSelector']").addClass('active');
    for (var i=0; i<species.length; i++) {
      $("a[species-id="+species[i]+"]").addClass('active');
    }
  }
  // Toggle annotation filter
  if (out['source']!=="Swissprot") {
    $("#AnnotationSP").removeClass('active');
    $("#AnnotationAll").addClass('active');
  }
  
  // Update not_found list
  var remove_indeces = [];
  for (i=0;i<slugs.length;i++) {
    var this_item = $('#'+slugs[i])[0];
      if (!selected_items.includes(this_item)) {
        addTarget(this_item, false);
      }
      remove_indeces.push(i);
  }
  var new_not_found = [];
  for (i=0;i<not_found.length;i++) {
    if (!remove_indeces.includes(i)) {
      new_not_found.push(not_found[i])
    }
  }
  for (i=0;i<full_names.length;i++) {
    if (!entries.includes(full_names[i])) {
      new_not_found.push(full_names[i]);
    }
  }
  not_found = new_not_found;

  // Add summary on message
  var message = "";
  var msg_type = "info";
  if (parsed > 0){
    message = "<b>Successfully</b> imported "+parsed+" entries.<br>";
    if (not_found.length > 0){
      message += "<br>The following name(s) could <i>not</i> be matched:<br>&#8226;&nbsp;&nbsp;" + not_found.join("<br>&#8226;&nbsp;&nbsp;");
    }
  } else {
    message = "The target selection import was <b>not successful</b>. Please make sure you are using the uniprot target names.";
    msg_type = "warning";
  }
  showAlert(message, msg_type);
  updateTargetCount();
}

function importTargetSelection(entry_names, response){
  $.ajax({
      'url': '/common/importtargetselection',
      'data': {
          entry_names: entry_names
      },
      'dataType': 'JSON',
      'type': 'GET',
      'async': false,
      'success': function(out) {
        response = JSON.parse(out);
      },
      'error': function(data) {
        console.log('error')
      }
  });
  return response;
}

/**
 * This function exports the current target selection to a textbox and copies the
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
    showAlert("Target selection has been copied to your clipboard.", "info");
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
function submitSelection(url, minimum = 1, maximum = 0) {
  // Check species selection and assume correct # species if not human
  let species = $("div#filters-species a.btn.active")[0].innerText.trim();
  let species_exemption = false;
  if (species !== "Human"){
    species_exemption = true;
  }
  if (species_exemption || (selected_targets.size >= minimum && (maximum === 0 || selected_targets.size <= maximum)) || imported_targets > 0) {
    // set CSRF csrf_token
    $.ajaxSetup({
        headers:
        { "X-CSRF-TOKEN": $("meta[name=\"csrf-token\"]").attr("content") }
    });

    // Submit proteins to target selection
    var group = Array.from(selected_targets);
    $.post("/common/targetformread", { "input-targets": group.join("\r") },  function (data) {
      // On success go to alignment page
      window.location.href = url;
    })
    .fail(
      function(){
        showAlert("Something went wrong, please try again or contact us.", "danger");
      });
  } else {
    if (maximum > 0 && selected_targets.size > maximum){
      if (maximum > 1) {
        showAlert("Please select a maximum of " + maximum + " targets.", "warning");
      } else {
        showAlert("Please select only one target.", "warning");
      }
    } else if (minimum === 1) {
      showAlert("Please select at least one target.", "warning");
    } else {
      showAlert("Please select at least "+minimum+" targets.", "warning");
    }
  }
}

/**
 * This function submits the selected targets to the backend and moves to the
 * next page. This function should be executed upon pressing the green button
 * after target selection.
 * @param {string} url The url to go to after synchronizing the target selection
 */
function onTheFlySelection(url, subtype=false, pathway=false, minimum = 1, maximum = 0) {
  // Check species selection and assume correct # species if not human
  let species = $("div#filters-species a.btn.active")[0].innerText.trim();
  let species_exemption = false;
  if (species !== "Human"){
    species_exemption = true;
  }
  console.log(selected_targets);
  if (species_exemption || (selected_targets.size >= minimum && (maximum === 0 || selected_targets.size <= maximum))) {
    // set CSRF csrf_token
    $.ajaxSetup({
        headers:
        { "X-CSRF-TOKEN": $("meta[name=\"csrf-token\"]").attr("content") }
    });

    // Submit proteins to target selection
    var group = Array.from(selected_targets);
    $.post("/common/targetformread", { "input-targets": group.join("\r") },  function (data) {
      // On success go to alignment page
      window.location.href = url;
    })
    .fail(
      function(){
        showAlert("Something went wrong, please try again or contact us.", "danger");
      });
  } else {
    if (maximum > 0 && selected_targets.size > maximum){
      if (maximum > 1) {
        showAlert("Please select a maximum of " + maximum + " targets.", "warning");
      } else {
        showAlert("Please select only one target.", "warning");
      }
    } else if (minimum === 1) {
      showAlert("Please select at least one target.", "warning");
    } else {
      showAlert("Please select at least "+minimum+" targets.", "warning");
    }
  }
}

/**
 * This function initializes the YADCF datatables for a specific element
 * @param {string} elementID The identifier of the container containing the table
 */
function initTargetTable(elementID) {

    if (!$.fn.DataTable.isDataTable(elementID + " table")) {
        targetTable = $(elementID + " table").DataTable({
//            dom: "ftip",
            deferRender: true,
            scrollY: "50vh",
            scrollX: true,
            scrollCollapse: true,
            scroller: true,
            paging: false,
//            bSortCellsTop: false, //prevent sort arrows going on bottom row
            aaSorting: [],
            autoWidth: false,
            bInfo: true,
//            order: [[ 1, "asc" ], [ 2, "asc" ], [ 3, "asc" ], [ 4, "asc" ]],
            columnDefs: [{
                targets: 0,
                orderable: false,
                className: "select-checkbox"
            },{
                targets: 1,
                className: "text-center"
            },{
                targets: [6,7,9,10,11,12],
                className: "text-center"
            },],
        });

        yadcf.init(targetTable,
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
                    column_data_type: "html",
                    filter_default_label: "Uniprot",
                    filter_reset_button_text: false,
                    select_type_options: {
                        "width": "110px",
                    }
                },
                {
                    column_number: 5,
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
                    column_number: 6,
                    filter_type: "range_number",
                    filter_default_label: ["From", "To"],
                    filter_reset_button_text: false,
                    //style_class: "center"
                    //html5_data: "data-search",
                    column_data_type: "html",
                },
                {
                    column_number: 7,
                    filter_type: "range_number",
                    filter_default_label: ["From", "To"],
                    filter_reset_button_text: false,
                    style_class: "center",
                },
                {
                    column_number: 8,
                    filter_type: "text",
                    select_type: "select2",
                    html5_data: "data-search",
                    filter_default_label: "PDB",
                    filter_reset_button_text: false,
                    select_type_options: {
                        "width": "110px",
                    }
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
                {
                    column_number: 9,
                    filter_type: "multi_select",
                    select_type: "select2",
                    filter_default_label: "Gs",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 10,
                    filter_type: "multi_select",
                    select_type: "select2",
                    filter_default_label: "Gi/o",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 11,
                    filter_type: "multi_select",
                    select_type: "select2",
                    filter_default_label: "Gq/11",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 12,
                    filter_type: "multi_select",
                    select_type: "select2",
                    filter_default_label: "G12/13",
                    filter_reset_button_text: false,
                },

            ], {
//                cumulative_filtering: true,
                filters_tr_index: 1
            }
        );
    }

    // When redrawing update the information selection message
    targetTable.on("draw.dt", function(e, oSettings) {
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
    targetTable.draw();
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

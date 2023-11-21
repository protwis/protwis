/*global yadcf*/
/*eslint complexity: ["error", 20]*/
/*eslint wrap-iife: ["error", "outside"]*/
/*eslint quotes: ["error", "double", { "avoidEscape": true }]*/

function hideColumns(table, columns) {
  table.columns(columns).visible(false, true);
  table.draw();
}

function resetHiddenColumns(table) {
  let col_length = table.columns()[0].length;
  let columns = Array.from(new Array(col_length - 1), (x, i) => i);
  table.columns(columns).visible(true, false);
  table.draw();
}

function reset_tab(table) {
  $("input.yadcf-filter-range-number").val("");
  yadcf.exResetAllFilters(table);
  resetHiddenColumns(table);
}

// Calculation of normalized rank for a given column
function createRank(table_id, column) {
  // Set default values for all cells
  $(table_id + " tbody tr td").filter(":nth-child(" + column + ")").each(function() {
    let cell_span = $(this.firstChild);
    let cell_value = cell_span.text();
    cell_span.attr("data-raw", cell_value);
    cell_span.attr("data-column-nr", column - 1);
    $(this).attr("data-sort", cell_value);
  });

  // Looping over the support values (GPCRdb # 1-3)
  for (let i = 1; i <= 3; i++) {
    // Step 1 - collect all values for a given column
    let min_max = [];
    $(table_id + " tbody tr[data-source='" + i + "'] td").filter(":nth-child(" + column + ")").each(function() {
      var cell_value = $(this).text();
      if (/^-?\d*(\.\d+)?$/.test(cell_value) && cell_value !== "-") {
        min_max.push(parseFloat(cell_value));
      }
    });

    // Step 2 - normalize all values and add them to a data attribute
    let min = Math.min(...min_max);
    let max = Math.max(...min_max);
    $(table_id + " tbody tr[data-source='" + i + "'] td").filter(":nth-child(" + column + ")").each(function() {
      let cell_span = $(this.firstChild);
      let cell_value = cell_span.text();

      // Check if data source is a number, if so confidence value for GPCRdb source
      if (/^-?\d*(\.\d+)?$/.test(cell_value) && cell_value !== "-") {
        cell_span.attr("data-normalized", Math.round((1 - (parseFloat(cell_value) - min) / (max - min)) * 100), 0);
      } else {
        // not value - empty init
        cell_span.attr("data-normalized", "-");
      }
    });
  }
}

/**
 * This is a custom YADCF function that checks min, max, and rank values
 * To enable this functionality, make sure to perform the createRank on the desired columns
 * And add the necessary input filters to the layout + create a hidden filter
 * @param {object} filterVal Value to filter on (not applicable)
 * @param {object} columnVal Element in the filtered column
 * @param {object} rowValues All elements in this row (not used)
 * @param {object} stateVal Current DOM state of the row (not sufficient in this case)
 * @returns {boolean} true if row contains selected target otherwise false
 */
let counter = 0;
let lastRangeRankFilter = "";

function rankedRangeFilter(filterVal, columnVal, rowValues, stateVal, tableNum) {
  let table_nr = tableNum;
  let column_value = $(columnVal).text();
  let column_nr = $(columnVal).attr("data-column-nr");
  let min_filtering = parseFloat($("#ranked_range_min" + table_nr + "_" + column_nr).val());
  let max_filtering = parseFloat($("#ranked_range_max" + table_nr + "_" + column_nr).val());
  let rank_filtering = parseFloat($("#ranked_range_rank" + table_nr + "_" + column_nr).val());

  if (!isNaN(rank_filtering) && lastRangeRankFilter !== "max" && lastRangeRankFilter !== "min") {
    // If filtering on rank - clean range filter
    $("#ranked_range_min" + table_nr + "_" + column_nr).val("");
    $("#ranked_range_max" + table_nr + "_" + column_nr).val("");

    let ranked_value = parseFloat($(columnVal).attr("data-normalized"));
    if (isNaN(ranked_value)) {
      return false;
    } else {
      return ranked_value <= rank_filtering;
    }
  } else if (!isNaN(min_filtering) || !isNaN(max_filtering)) {
    // If filtering on range - clean rank filter
    $("#ranked_range_rank" + table_nr + "_" + column_nr).val("");

    // Filter range on current columnVal
    let range_value = parseFloat(column_value);
    if (isNaN(range_value)) {
      return false;
    } else {
      if (!isNaN(min_filtering) && !isNaN(max_filtering)) {
        return range_value >= min_filtering && range_value <= max_filtering;
      } else if (!isNaN(min_filtering)) {
        return range_value >= min_filtering;
      } else if (!isNaN(max_filtering)) {
        return range_value <= max_filtering;
      } else {
        // Should never happen
      }
    }
  } else {
    return true;
  }
}

/**
 * This is a custom YADCF function that checks ....
 * ....
 * @param {object} filterVal Value to filter on (not applicable)
 * @param {object} columnVal Element in the filtered column
 * @param {object} rowValues All elements in this row (not used)
 * @param {object} stateVal Current DOM state of the row (not sufficient in this case)
 * @returns {boolean} true if row contains selected target otherwise false
 */
function supportFilter(filterVal, columnVal, rowValues, stateVal) {
  //console.log(!/^\d+$/.test(columnVal), columnVal, filterVal);
  //console.log(columnVal === filterVal);
  return (!/^\d+$/.test(columnVal) || columnVal === filterVal);
}

/**
 * Function copied from contactbrowser-tabtables.js
 * When there's a need to repeat the same yadcf filter_type one can use this function to concatenate
 * the range_number filter_type.
 */
function make_rank_col_filters(start_column, repeat_number, cont_prefix, default_filter) {
  let repeated_filter = [];

  for (let i = start_column; i < start_column + repeat_number; i++) {
    let column_info = {
      filter_type: "custom_func",
      column_data_type: "text",
      column_number: i,
      filter_container_id: cont_prefix + i,
      custom_func: default_filter,
    };
    repeated_filter.push(column_info);
  }
  return repeated_filter;
}

// =============================================================================
// START OVERLAY COLUMNS CODE HERE
// =============================================================================

let isScrolling;
let leftOverlayOffset = [];

function createFixedColumnsOverlay(table_id, clone_cols = 6) {
  // When called after filtering: remove current overlay rows
  $("#overlay_" + table_id + " table tbody tr").remove();

  // Clone first X columns for each row
  $("#" + table_id + " tbody tr").each(function() {
    let tds = $(this).children();
    let overlay_row = $("<tr></tr>");

    for (let i = 0; i < clone_cols; i++) {
      overlay_row.append(tds.eq(i).clone());
    }
    $("#overlay_" + table_id + " table tbody").append(overlay_row);
    overlay_row.height($(this).height());
  });
}

function initFixedColumnsOverlay(table_id) {
  // Create overlay container and hide
  $("#" + table_id).closest(".dataTables_scrollBody").append('<div id="overlay_' + table_id + '" class="table_overlay"><table class="row-border text-center compact dataTable no-footer text-nowrap"><tbody></tbody></table></div>');
  $("#overlay_" + table_id).hide();

  // Populate overlay container
  createFixedColumnsOverlay(table_id);

  // Detect filtering events on table and repopulate the overlay
  $("#" + table_id).on("draw.dt", function() {
    createFixedColumnsOverlay(table_id);
  });

  // Detect scrolling events on table and update overlay position and visibility
  leftOverlayOffset[table_id] = 0;
  $("#" + table_id).closest(".dataTables_scrollBody").scroll(function() {
    let left_tmp = $("#" + table_id).closest(".dataTables_scrollBody").scrollLeft();
    if (left_tmp !== leftOverlayOffset[table_id]) {
      leftOverlayOffset[table_id] = left_tmp;

      // Clear our timeout while still scrolling
      window.clearTimeout(isScrolling);

      // Hide overlay if still visible
      if ($("#overlay_" + table_id).is(":visible")) {
        $("#overlay_" + table_id).hide();
      }

      // Run overlay update only when scrolling ended
      isScrolling = setTimeout(function() {
        // If user scrolls and it's > 100px from left, then attach fixed columns overlay
        let left_scroll = $("#" + table_id).closest(".dataTables_scrollBody").scrollLeft();
        if (left_scroll > 50) {
          $("#overlay_" + table_id).css({
            left: left_scroll + "px"
          });
          if ($("#overlay_" + table_id).is(":hidden")) {
            $("#overlay_" + table_id).show();
          }
        }
      }, 70);
    }
  });
}

// =============================================================================
// END OVERLAY COLUMNS CODE HERE
// =============================================================================

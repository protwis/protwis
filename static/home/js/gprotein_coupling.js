/*global yadcf, rankedRangeFilter, createRank, createYADCFfilters, supportFilter, make_rank_col_filters, lastRangeRankFilter, initFixedColumnsOverlay, copyListToClipboard*/
/*eslint complexity: ["error", 20]*/
/*eslint wrap-iife: ["error", "outside"]*/
/*eslint quotes: ["error", "double", { "avoidEscape": true }]*/
let oTable1 = [];
let oTable2 = [];

function rankedRangeFiltert1(filterVal, columnVal, rowValues, stateVal) {
  return rankedRangeFilter(filterVal, columnVal, rowValues, stateVal, "1");
}

function rankedRangeFiltert2(filterVal, columnVal, rowValues, stateVal) {
  return rankedRangeFilter(filterVal, columnVal, rowValues, stateVal, "2");
}

$(document).ready(function() {
  // Activate tooltips and popovers from Bootstrap
  $("[data-toggle='tooltip']").tooltip();
  $("[data-toggle='popover']").popover();

  // Create the ranks for the families table
  for (let i = 14; i <= 24; i++) {
    createRank("#familiestabletab", i); // GS
  }

  // =============================================================================
  // Families Table
  // =============================================================================
  console.time("table1load");
  $("#familiestabletab").show();
  oTable1 = $("#familiestabletab").DataTable({
    deferRender: true,
    scrollY: "50vh",
    scrollX: true,
    scrollCollapse: true,
    scroller: true,
    paging: false,
    bSortCellsTop: false, //prevent sort arrows going on bottom row
    aaSorting: [],
    order: [
      [2, "asc"],
      [4, "asc"],
      [24, "asc"],
    ],
    autoWidth: false,
    bInfo: true,
    columnDefs: [{
      targets: [24],
      visible: false
    }],
  });

  let column_filters = [];
  // Selector column
  column_filters = column_filters.concat(createYADCFfilters(0, 1, "none"));
  // Receptor section
  column_filters = column_filters.concat(createYADCFfilters(1, 1, "multi_select", "select2", "", false, null, "html", "80px"));
  // column_filters = column_filters.concat(createYADCFfilters(2, 1, "multi_select", "select2", "", false, null, null, "80px"));
  column_filters = column_filters.concat(createYADCFfilters(2, 1, "multi_select", "select2", "", false, null, null, "40px"));
  column_filters = column_filters.concat(createYADCFfilters(3, 1, "multi_select", "select2", "", false, null, null, "100px"));
  column_filters = column_filters.concat(createYADCFfilters(4, 1, "multi_select", "select2", "", false, "exact", "html", "80px"));
  column_filters = column_filters.concat(createYADCFfilters(5, 1, "multi_select", "select2", "", false, "exact", "html", "80px"));
  // Ligands section
  column_filters = column_filters.concat(createYADCFfilters(6, 1, "multi_select", "select2", "", false, "exact", "html", "100px"));
  column_filters = column_filters.concat(createYADCFfilters(7, 1, "multi_select", "select2", "", false, null, null, "50px"));
  // Guide to Pharmacology section
  column_filters = column_filters.concat(createYADCFfilters(8, 4, "multi_select", "select2", "", false, null, null, "40px"));
  // Coupling data filters
  column_filters = column_filters.concat(make_rank_col_filters(12, 12, "hide_rankfam", rankedRangeFiltert1));
  // Hidden GPCRdb support type column calls customized function
  column_filters = column_filters.concat([
    {
      column_number: 24,
      filter_type: "custom_func",
      custom_func: supportFilter,
      filter_container_id: "hide_filter1",
    }]);

  yadcf.init(oTable1, column_filters, {
    cumulative_filtering: false
  });

  // Initialize ranked Range Filtering options
  $(".ranked_range_min1, .ranked_range_max1, .ranked_range_rank1").on("click", function(event) {
    event.stopPropagation();
  });

  $(".ranked_range_min1, .ranked_range_max1, .ranked_range_rank1").on("input", function(event) {
    // Get column #
    let column_nr = event.target.id.split("_")[3];

    // Store current type of filtering globally
    lastRangeRankFilter = event.target.id.split("_")[2];

    // SELECT one option in real YADCF filter to trick YADCF into calling the filter function
    let adjust_node1 = $("#yadcf-filter--familiestabletab-" + column_nr + " option").filter(":nth-child(2)").first();
    adjust_node1.prop("selected", true);

    // Invoke filtering
    $("#yadcf-filter--familiestabletab-" + column_nr).change();

    // Clean filter type
    lastRangeRankFilter = "";
  });


  // This prefilters the value 2 in the hidden column 29 which corresponds to being in at least two of the supporting GPCRdb
  // datasets
  yadcf.exFilterColumn(oTable1, [
    [24, 2]
  ]);
  //    yadcf.exResetAllFilters(oTable1);

  //  Select clicked-on boxes for families table
  $("#familiestabletab" + " > tbody > tr").click(function(event) {
    if (event.target.type !== "checkbox") {
      $(":checkbox", this).trigger("click");
      $(this).eq(0).toggleClass("alt_selected");
      $(this).find("td").toggleClass("highlight");
    }
    $(this).eq(0).toggleClass("alt_selected");
    $(this).find("td").toggleClass("highlight");
  });
  console.timeEnd("table1load");

  // Select all boxes in table
  $(".select-all").click(function() {
    $(":checkbox", this).trigger("click");
    if ($(this).prop("checked") === true) {
      $(".alt").prop("checked", true);
      $(".alt").parent().parent().addClass("alt_selected");
      $(".alt").parent().parent().find("td").addClass("highlight");
    }
    if ($(this).prop("checked") === false) {
      $(".alt").prop("checked", false);
      $(".alt").parent().parent().removeClass("alt_selected");
      $(".alt").parent().parent().find("td").removeClass("highlight");
    }
  });

  // Hide column button for table1
  $(".hide_columns1").click(function(evt) {
    var columns = $(this).attr("columns").split(",");
    oTable1.columns(columns).visible(false, false);
    oTable1.draw();
  });

  // Put top scroller
  // https://stackoverflow.com/questions/35147038/how-to-place-the-datatables-horizontal-scrollbar-on-top-of-the-table
  $(".dataTables_scrollHead").css({
    "overflow-x": "scroll"
  }).on("scroll", function(e) {
    var scrollBody = $(this).parent().find(".dataTables_scrollBody").get(0);
    scrollBody.scrollLeft = this.scrollLeft;
    $(scrollBody).trigger("scroll");
  });

  // Enable columns overlay
  initFixedColumnsOverlay("familiestabletab");

  // Enable copy to clipboard option
  $("#uniprot_copy1").click(function() {
    let selected = [];
    $(".alt_selected > .uniprot1 > a").each(function() {
      selected.push($(this).text());
    });
    copyListToClipboard(selected);
  });

  // Uncheck every row when using back button on browser
  $(".alt_selected").prop("checked", false);
  $(".alt").prop("checked", false);
  $(".select-all").prop("checked", false);
});

// =============================================================================
// Subtypes Table
// =============================================================================

let table2_initialized = false;
function initCouplingTable2() {
  if (table2_initialized) {
    return;
  }

  table2_initialized = true;

  // Create the ranks for the subtypes table
  for (let i = 14; i <= 68; i++) {
    createRank("#subtypestabletab", i); // GS
  }

  console.time("table2load");
  $("#subtypestabletab").show();
  oTable2 = $("#subtypestabletab").DataTable({
    deferRender: true,
    scrollY: "50vh",
    scrollX: true,
    scrollCollapse: true,
    scroller: true,
    paging: false,
    bSortCellsTop: false, //prevent sort arrows going on bottom row
    aaSorting: [],
    order: [
      [3, "asc"],
      [5, "asc"],
      [68, "asc"]
    ],
    autoWidth: false,
    bInfo: true,
    columnDefs: [{
      targets: [68],
      visible: false
    }],
  });

  let column_filters = [];
  // Selector column
  column_filters = column_filters.concat(createYADCFfilters(0, 1, "none"));
  // Receptor section
  column_filters = column_filters.concat(createYADCFfilters(1, 1, "multi_select", "select2", "", false, null, "html", "80px"));
  // column_filters = column_filters.concat(createYADCFfilters(2, 1, "multi_select", "select2", "", false, null, null, "80px"));
  column_filters = column_filters.concat(createYADCFfilters(2, 1, "multi_select", "select2", "", false, null, null, "40px"));
  column_filters = column_filters.concat(createYADCFfilters(3, 1, "multi_select", "select2", "", false, null, null, "100px"));
  column_filters = column_filters.concat(createYADCFfilters(4, 1, "multi_select", "select2", "", false, "exact", "html", "80px"));
  column_filters = column_filters.concat(createYADCFfilters(5, 1, "multi_select", "select2", "", false, "exact", "html", "80px"));
  // Ligands section
  column_filters = column_filters.concat(createYADCFfilters(6, 1, "multi_select", "select2", "", false, "exact", "html", "100px"));
  column_filters = column_filters.concat(createYADCFfilters(7, 1, "multi_select", "select2", "", false, null, null, "50px"));
  // Guide to Pharmacology section
  column_filters = column_filters.concat(createYADCFfilters(8, 4, "multi_select", "select2", "", false, null, null, "40px"));
  // Coupling data filters
  column_filters = column_filters.concat(make_rank_col_filters(12, 56, "hide_ranksub", rankedRangeFiltert2));
  // Hidden GPCRdb support type column calls customized function
  column_filters = column_filters.concat([
    {
      column_number: 68,
      filter_type: "custom_func",
      custom_func: supportFilter,
      filter_container_id: "hide_filter2",
    }]);

  yadcf.init(oTable2, column_filters, {
    cumulative_filtering: false
  });

  // Initialize ranked Range Filtering options
  $(".ranked_range_min2, .ranked_range_max2, .ranked_range_rank2").on("click", function(event) {
    event.stopPropagation();
  });

  $(".ranked_range_min2, .ranked_range_max2, .ranked_range_rank2").on("input", function(event) {
    // Get column #
    let column_nr = event.target.id.split("_")[3];

    // Store current type of filtering globally
    lastRangeRankFilter = event.target.id.split("_")[2];

    // SELECT one option in real YADCF filter to trick YADCF into calling the filter function
    let adjust_node1 = $("#yadcf-filter--subtypestabletab-" + column_nr + " option").filter(":nth-child(2)").first();
    adjust_node1.prop("selected", true);

    // Invoke filtering
    $("#yadcf-filter--subtypestabletab-" + column_nr).change();

    // Clean filter type
    lastRangeRankFilter = "";
  });
  yadcf.exFilterColumn(oTable2, [
    [68, 2]
  ]);

  //  Select clicked-on boxes for subtypes table
  $("#subtypestabletab" + " > tbody > tr").click(function(event) {
    if (event.target.type !== "checkbox") {
      $(":checkbox", this).trigger("click");
      $(this).eq(0).toggleClass("alt_selected");
      $(this).find("td").toggleClass("highlight");
    }
    $(this).eq(0).toggleClass("alt_selected");
    $(this).find("td").toggleClass("highlight");
  });

  // Hide column button for table2
  $(".hide_columns2").click(function(evt) {
    var columns = $(this).attr("columns").split(",");
    oTable2.columns(columns).visible(false, false);
    oTable2.draw();
  });

  initFixedColumnsOverlay("subtypestabletab");

  // Enable copy to clipboard option
  $("#uniprot_copy2").click(function() {
    let selected = [];
    $(".alt_selected > .uniprot2 > a").each(function() {
      selected.push($(this).text());
    });
    copyListToClipboard(selected);
  });

  console.timeEnd("table2load");
}

/*global yadcf, rankedRangeFilter, createRank, createYADCFfilters, supportFilter, make_rank_col_filters, lastRangeRankFilter, initFixedColumnsOverlay, copyListToClipboard*/
/*eslint complexity: ["error", 20]*/
/*eslint wrap-iife: ["error", "outside"]*/
/*eslint quotes: ["error", "double", { "avoidEscape": true }]*/
let oTable1 = [];

$(document).ready(function() {
  // Activate tooltips and popovers from Bootstrap
  $("[data-toggle='tooltip']").tooltip();
  $("[data-toggle='popover']").popover();

  // =============================================================================
  // Families Table
  // =============================================================================
  console.time("table1load");
  $("#familiestabletab").show();
  oTable1 = $("#familiestabletab").DataTable({
    deferRender: true,
    scrollY: '65vh',
    scrollX: true,
    scrollCollapse: true,
    scroller: true,
    paging: false,
    bSortCellsTop: false, //prevent sort arrows going on bottom row
    aaSorting: [],
    order: [
      [5, "asc"],
      [1, "asc"]
    ],
    autoWidth: true,
    bInfo: true,
  });

  let column_filters = [];
  // Selector column
  column_filters = column_filters.concat(createYADCFfilters(0, 1, "none"));
  // Receptor section
  
  column_filters = column_filters.concat(createYADCFfilters(1, 1, "multi_select", "select2", "", false, "exact", "html", "80px"));
  column_filters = column_filters.concat(createYADCFfilters(2, 1, "multi_select", "select2", "", false, "exact", "html", "80px"));
  column_filters = column_filters.concat(createYADCFfilters(3, 1, "none"));
  column_filters = column_filters.concat(createYADCFfilters(4, 1, "multi_select", "select2", "", false, "exact", null, "80px"));
  column_filters = column_filters.concat(createYADCFfilters(5, 1, "multi_select", "select2", "", false, null, "html", "45px"));
  column_filters = column_filters.concat(createYADCFfilters(6, 1, "multi_select", "select2", "", false, null, null, "100px"));
  //Source section
  column_filters = column_filters.concat(createYADCFfilters(7, 2, "multi_select", "select2", "", false, null, null, "90px"));
  column_filters = column_filters.concat(createYADCFfilters(9, 1, "multi_select", "select2", "", false, null, null, "35px"));
  column_filters = column_filters.concat(createYADCFfilters(10, 1, "multi_select", "select2", "", false, null, null, "40px"));
  column_filters = column_filters.concat(createYADCFfilters(11, 1, "none"));
  // Ligands section
  column_filters = column_filters.concat(createYADCFfilters(12, 1, "multi_select", "select2", "", false, "exact", "html", "100px"));
  column_filters = column_filters.concat(createYADCFfilters(13, 1, "multi_select", "select2", "", false, null, null, "50px"));
  // Supporting datasets, Guide to Pharmacology and family rank orders section
  column_filters = column_filters.concat(createYADCFfilters(14, 4, "multi_select", "select2", "", false, null, "html", "40px"));
  column_filters = column_filters.concat(createYADCFfilters(18, 11, "multi_select", "select2", "", false, null, null, "40px"));
  column_filters = column_filters.concat(createYADCFfilters(29, 4, "range_number", null, ["Min", "Max"], false, null, null, "40px", null, "-"));
  column_filters = column_filters.concat(createYADCFfilters(33, 1, "multi_select", "select2", "", false, null, "html", "90px"));
  column_filters = column_filters.concat(createYADCFfilters(34, 4, "range_number", null, ["Min", "Max"], false, null, null, "40px", null, "-"));

  column_filters = column_filters.concat(createYADCFfilters(38, 1, "multi_select", "select2", "", false, null, null, "60px"));
  column_filters = column_filters.concat(createYADCFfilters(39, 16, "range_number", null, ["Min", "Max"], false, null, null, "40px", null, "-"));
  column_filters = column_filters.concat(createYADCFfilters(55, 1, "multi_select", "select2", "", false, null, "html", "90px"));
  column_filters = column_filters.concat(createYADCFfilters(56, 16, "range_number", null, ["Min", "Max"], false, null, null, "40px", null, "-"));
  // Coupling data filters

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
    oTable1.columns.adjust();
    $(".no-sort").removeClass("sorting");
  });

  oTable1.columns.adjust();

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

  $(".no-sort").removeClass("sorting");

  $("#familiestabletab").on("order.dt", function() {
    $(".no-sort").removeClass("sorting");
  })

  $(".G1213").addClass("rightborder");
  $(".G13").addClass("rightborder");

});

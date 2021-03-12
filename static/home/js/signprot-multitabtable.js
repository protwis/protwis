/*global yadcf*/
/*eslint complexity: ["error", 15]*/
/*eslint wrap-iife: ["error", "outside"]*/
/*eslint quotes: ["error", "double", { "avoidEscape": true }]*/
let oTable1 = [];
let oTable2 = [];


var tableToExcel = (function () {
//function tableToExcel() {
    var uri = "data:application/vnd.ms-excel;base64,",
        template = "<html xmlns:o='urn:schemas-microsoft-com:office:office' xmlns:x='urn:schemas-microsoft-com:office:excel' xmlns='http://www.w3.org/TR/REC-html40'><head><!--[if gte mso 9]><xml><x:ExcelWorkbook><x:ExcelWorksheets><x:ExcelWorksheet><x:Name>{worksheet}</x:Name><x:WorksheetOptions><x:DisplayGridlines/></x:WorksheetOptions></x:ExcelWorksheet></x:ExcelWorksheets></x:ExcelWorkbook></xml><![endif]--></head><body><table>{table}</table></body></html>",
        base64 = function (s) {
            return window.btoa(unescape(encodeURIComponent(s)));
        }, format = function (s, c) {
            return s.replace(/{(\w+)}/g, function (m, p) {
                return c[p];
            });
        };
    return function (table, name, filename) {
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
                $(th).html($(th).attr("title"));
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
    };
}());

function select_all(e) {
    var checkedStatus = $(e).prop("checked");

    $(".select-all  ").each(function () {
        $(this).prop("checked", checkedStatus);
    });

    $(".alt").each(function () {
        $(this).prop("checked", checkedStatus);
    });
}

function resetHidden1() {
    var columns = Array.from(new Array(17), (x,i) => i + 3);
    columns.forEach(function(column) {
        column = oTable1.column( column );
        try {
            column.visible( true, false );
        }
        catch(err) {
            column.visible( true, false );
        }
    });
    oTable1.draw();
}

function resetHidden2() {
    var columns = Array.from(new Array(28), (x,i) => i + 3);
    columns.forEach(function(column) {
//        console.log("columns variable " + columns);
        column = oTable2.column( column );
        try {
            column.visible( true, false );
        }
        catch(err) {
            column.visible( true, false );
        }
    });
    oTable2.draw();
}

function reset_tab() {
// Just a button to go back to the main page.
// TODO: A better per tab (table object) reset.
    window.location.href = "/signprot/couplings";
}

// Calculation of normalized rank for a given column
function createRank(table_id, column) {
    // Set default values for all cells
    $(table_id+" tbody tr td").filter(":nth-child(" + column + ")").each( function() {
        let cell_span = $(this.firstChild);
        let cell_value = cell_span.text();
        cell_span.attr("data-raw", cell_value);
        cell_span.attr("data-column-nr", column-1);
        $(this).attr("data-sort", cell_value);
    });

    // Looping over the support values (GPCRdb # 1-3)
    for (let i=1; i <= 3; i++){
      // Step 1 - collect all values for a given column
      let min_max = [];
      $(table_id+" tbody tr[data-source='" + i + "'] td").filter(":nth-child(" + column + ")").each( function() {
          var cell_value = $(this).text();
          if (/^-?\d*(\.\d+)?$/.test(cell_value) && cell_value!=="-"){
              min_max.push(parseFloat(cell_value));
          }
      });

      // Step 2 - normalize all values and add them to a data attribute
      let min = Math.min(...min_max);
      let max = Math.max(...min_max);
      $(table_id+" tbody tr[data-source='" + i + "'] td").filter(":nth-child(" + column + ")").each( function() {
          let cell_span = $(this.firstChild);
          let cell_value = cell_span.text();

          // Check if data source is a number, if so confidence value for GPCRdb source
          if (/^-?\d*(\.\d+)?$/.test(cell_value) && cell_value!=="-"){
              cell_span.attr("data-normalized", Math.round((1-(parseFloat(cell_value)-min)/(max-min))*100),0);
          } else {
              // not value - empty init
              cell_span.attr("data-normalized","-");
          }
      });
    }
}

// # custom rankedRangeFilter draft for YADCF
/**
 * This is a custom YADCF function that checks ....
 * ....
 * @param {object} filterVal Value to filter on (not applicable)
 * @param {object} columnVal Element in the filtered column
 * @param {object} rowValues All elements in this row (not used)
 * @param {object} stateVal Current DOM state of the row (not sufficient in this case)
 * @returns {boolean} true if row contains selected target otherwise false
 */
let counter = 0;
function rankedRangeFiltert1(filterVal, columnVal, rowValues, stateVal) {
    // DEBUG
    /*if (counter < 1) {
        counter++;
        console.log("FILTERING FOR", filterVal);
        console.log(columnVal);
        console.log(rowValues);
        console.log(stateVal);
    }*/

    let column_value = $(columnVal).text();
    let column_nr = $(columnVal).attr("data-column-nr");
    let min_filtering = parseFloat($("#ranked_range_min1_" + column_nr).val());
    let max_filtering = parseFloat($("#ranked_range_max1_" + column_nr).val());
    let rank_filtering = parseFloat($("#ranked_range_rank1_" + column_nr).val());

    if (!isNaN(rank_filtering) && lastRangeRankFilter!=="max" && lastRangeRankFilter!=="min") {
      // If filtering on rank - clean range filter
      $("#ranked_range_min1_" + column_nr).val("");
      $("#ranked_range_max1_" + column_nr).val("");

      let ranked_value = parseFloat($(columnVal).attr("data-normalized"));
      if (isNaN(ranked_value)) {
          return false;
      } else {
          return ranked_value <= rank_filtering;
      }
    } else if (!isNaN(min_filtering) || !isNaN(max_filtering)) {
        // If filtering on range - clean rank filter
        $("#ranked_range_rank1_" + column_nr).val("");

        // Filter range on current columnVal
        let range_value = parseFloat(column_value);
         if (isNaN(range_value)) {
              return false;
          } else {
             if (!isNaN(min_filtering) && !isNaN(max_filtering)) {
                return range_value >= min_filtering && range_value <= max_filtering;
             } else if (!isNaN(min_filtering)) {
                  return range_value >= min_filtering;
              } else if(!isNaN(max_filtering)) {
                 return range_value <= max_filtering;
             } else {
                  // Should never happen
              }
          }
    } else {
        return true;
    }
}

function rankedRangeFiltert2(filterVal, columnVal, rowValues, stateVal) {
    // DEBUG
    /*if (counter < 1) {
        counter++;
        console.log("FILTERING FOR", filterVal);
        console.log(columnVal);
        console.log(rowValues);
        console.log(stateVal);
    }*/

    let column_value = $(columnVal).text();
    let column_nr = $(columnVal).attr("data-column-nr");
    let min_filtering = parseFloat($("#ranked_range_min2_" + column_nr).val());
    let max_filtering = parseFloat($("#ranked_range_max2_" + column_nr).val());
    let rank_filtering = parseFloat($("#ranked_range_rank2_" + column_nr).val());

    if (!isNaN(rank_filtering) && lastRangeRankFilter!=="max" && lastRangeRankFilter!=="min") {
      // If filtering on rank - clean range filter
      $("#ranked_range_min2_" + column_nr).val("");
      $("#ranked_range_max2_" + column_nr).val("");

      let ranked_value = parseFloat($(columnVal).attr("data-normalized"));
      if (isNaN(ranked_value)) {
          return false;
      } else {
          return ranked_value <= rank_filtering;
      }
    } else if (!isNaN(min_filtering) || !isNaN(max_filtering)) {
        // If filtering on range - clean rank filter
        $("#ranked_range_rank2_" + column_nr).val("");

        // Filter range on current columnVal
        let range_value = parseFloat(column_value);
         if (isNaN(range_value)) {
              return false;
          } else {
             if (!isNaN(min_filtering) && !isNaN(max_filtering)) {
                return range_value >= min_filtering && range_value <= max_filtering;
             } else if (!isNaN(min_filtering)) {
                  return range_value >= min_filtering;
              } else if(!isNaN(max_filtering)) {
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
function supportFilter(filterVal, columnVal, rowValues, stateVal){
    //console.log(!/^\d+$/.test(columnVal), columnVal, filterVal);
    //console.log(columnVal === filterVal);
    return (!/^\d+$/.test(columnVal) || columnVal === filterVal);
}

/**
 * Function copied from contactbrowser-tabtables.js
 * When there's a need to repeat the same yadcf filter_type one can use this function to concatenate
 * the range_number filter_type.
 */
function make_range_number_cols(start_column, repeat_number, tab) {
    let from_to1 = {
        filter_type: "custom_func",
        custom_func: rankedRangeFiltert1,
    };
    let from_to2 = {
        filter_type: "custom_func",
        custom_func: rankedRangeFiltert2,
    };
    let repeated_from_to1 = [];
    let repeated_from_to2 = [];
    if (tab == "famtab") {
        for (let i = start_column; i < start_column + repeat_number; i++) {
            let column_info = Object.assign({}, from_to1);
            column_info["column_number"] = i;
            column_info["filter_container_id"] = "hide_rankfam" + i;
            repeated_from_to1.push(column_info);
        }
        return repeated_from_to1;
    } else if (tab == "subtab") {
        for (let i = start_column; i < start_column + repeat_number; i++) {
            let column_info = Object.assign({}, from_to2);
            column_info["column_number"] = i;
            column_info["filter_container_id"] = "hide_ranksub" + i;
            repeated_from_to2.push(column_info);
        }
        return repeated_from_to2;
    } else {
        //
    }
}

repfilterfamtab = make_range_number_cols(10, 12, "famtab");
repfiltersubtab = make_range_number_cols(10, 40, "subtab");

let lastRangeRankFilter = "";

$(document).ready(function() {
// Activate tooltips and popovers from Bootstrap   * Bootstrap v3.3.7 (http://getbootstrap.com)
    $("[data-toggle='tooltip']").tooltip();
    $("[data-toggle='popover']").popover();

// Create the ranks for the families table
for (let i=11; i <= 22; i++) {
    createRank("#familiestabletab", i); // GS
}

// ===============
// Families Table
// ===============
    console.time("table1load");
    oTable1 = $("#familiestabletab").DataTable({
        deferRender: true,
        scrollY: "50vh",
        scrollX: true,
        scrollCollapse: true,
        scroller: true,
        paging: false,
        bSortCellsTop: false, //prevent sort arrows going on bottom row
        aaSorting: [],
        order: [[4, "asc"], [22, "desc"]],
        autoWidth: false,
        bInfo: true,
        columnDefs: [
            {
                targets: [22],
                visible: false
            }
        ],
    });


    yadcf.init(oTable1,
        [
            {
                column_number: 0,
                filter_type: "none",
                filter_default_label: "",
                filter_reset_button_text: false,
            },
            {
                column_number: 1,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "80px",
                }
            },

            {
                column_number: 2,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "40px",
                }
            },
            {
                column_number: 3,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "200px",
                }
            },
            {
                column_number: 4,
                filter_type: "multi_select",
                select_type: "select2",
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "UniProt",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "60px",
                }
            },
            {
                column_number: 5,
                filter_type: "multi_select",
                select_type: "select2",
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "80px",
                }
            },

// Guide to Pharmacology
            {
                column_number: 6,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "40px"
                },
            },
            {
                column_number: 7,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "40px"
                },
            },
            {
                column_number: 8,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "40px"
                },
            },
            {
                column_number: 9,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "40px"
                },
            },

// Hidden GPCRdb support type column calls customized function
            {
                column_number: 22,
                filter_type: "custom_func",
                custom_func: supportFilter,
                filter_container_id: "hide_filter1",
            },
        ].concat(repfilterfamtab),

        {filters_tr_index: 2},

        {
            cumulative_filtering: false
        }
    );

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


// This prefilters the value 2 in the hidden column 22 which corresponds to being in at least two of the supporting GPCRdb
// datasets
    yadcf.exFilterColumn(oTable1, [[22, 2]]);
//    yadcf.exResetAllFilters(oTable1);

//  Select clicked-on boxes for families table
    $("#familiestabletab"+" > tbody > tr").click(function(event) {
        if (event.target.type !== "checkbox") {
            $(":checkbox", this).trigger("click");
            $(this).eq(0).toggleClass("alt_selected");
            $(this).find("td").toggleClass("highlight");
        }
        $(this).eq(0).toggleClass("alt_selected");
        $(this).find("td").toggleClass("highlight");
    });

    console.timeEnd("table1load");








// Create the ranks for the subtypes table
for (let i=11; i <= 49; i++) {
    createRank("#subtypestabletab", i); // GS
}
// ===============
// Subtypes Table
// ===============
    console.time("table2load");
    oTable2 = $("#subtypestabletab").DataTable({
        deferRender: true,
        scrollY: "50vh",
        scrollX: true,
        scrollCollapse: true,
        scroller: true,
        paging: false,
        bSortCellsTop: false, //prevent sort arrows going on bottom row
        aaSorting: [],
        order: [4,"asc"],
        autoWidth: false,
        bInfo: true,
        columnDefs: [
            {
                targets: [49],
                visible: false
            }
        ],
    });

    yadcf.init(oTable2,
        [
            {
                column_number: 0,
                filter_type: "none",
                filter_default_label: "",
                filter_reset_button_text: false,
            },
            {
                column_number: 1,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "80px",
                }
            },

            {
                column_number: 2,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "40px",
                }
            },
            {
                column_number: 3,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "200px",
                }
            },
            {
                column_number: 4,
                filter_type: "multi_select",
                select_type: "select2",
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "UniProt",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "60px",
                }
            },
            {
                column_number: 5,
                filter_type: "multi_select",
                select_type: "select2",
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "80px",
                }
            },

// Guide to Pharmacology
            {
                column_number: 6,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "40px"
                },
            },
            {
                column_number: 7,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "40px"
                },
            },
            {
                column_number: 8,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "40px"
                },
            },
            {
                column_number: 9,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "40px"
                },
            },

// Hidden GPCRdb support type column calls customized function
            {
                column_number: 49,
                filter_type: "custom_func",
                custom_func: supportFilter,
                filter_container_id: "hide_filter2",
            },
        ].concat(repfiltersubtab),

        {
            filters_tr_index: 3
        },

        {
            cumulative_filtering: false
        }
    );

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



    yadcf.exFilterColumn(oTable2, [[49, 2]]);

//  Select clicked-on boxes for subtypes table
    $("#subtypestabletab"+" > tbody > tr").click(function(event) {
        if (event.target.type !== "checkbox") {
            $(":checkbox", this).trigger("click");
            $(this).eq(0).toggleClass("alt_selected");
            $(this).find("td").toggleClass("highlight");
        }
        $(this).eq(0).toggleClass("alt_selected");
        $(this).find("td").toggleClass("highlight");
    });

// Select all boxes in table
    $(".select-all").click(function() {
        $(":checkbox", this).trigger("click");
        if ($(this).prop("checked")===true) {
            $(".alt").prop("checked", true);
            $(".alt").parent().parent().addClass("alt_selected");
            $(".alt").parent().parent().find("td").addClass("highlight");
        }
        if ($(this).prop("checked")===false) {
            $(".alt").prop("checked", false);
            $(".alt").parent().parent().removeClass("alt_selected");
            $(".alt").parent().parent().find("td").removeClass("highlight");
        }
    });

//    yadcf.exResetAllFilters(oTable2);
    console.timeEnd("table2load");







// =============================================================================
// GENERAL OPTIONS
// =============================================================================

// By default display the first tab. If this is not ON, one has to click on the tab for display.
    $("#couplingtabs a:first").tab("show");
//    $('#couplingtabs a[href="#table_1"]').tab('show');

// Hide column button for table1
    $(".hide_columns1").click(function(evt) {
        var columns = $(this).attr("columns").split(",");
        columns.forEach(function(column) {
            column = oTable1.column( column );
            try {
                column.visible( false, false );
            }
            catch(err) {
                column.visible( false, false );
            }
        });
        oTable1.draw();
    } );

// Hide column button for table2
    $(".hide_columns2").click(function(evt) {
        var columns = $(this).attr("columns").split(",");
        columns.forEach(function(column) {
            column = oTable2.column( column );
            try {
                column.visible( false, false );
            }
            catch(err) {
                column.visible( false, false );
            }
        });
        oTable2.draw();
    } );

// Put top scroller
// https://stackoverflow.com/questions/35147038/how-to-place-the-datatables-horizontal-scrollbar-on-top-of-the-table
//    console.time("scroll to top");
    $(".dataTables_scrollHead").css({
        "overflow-x":"scroll"
    }).on("scroll", function(e){
        var scrollBody = $(this).parent().find(".dataTables_scrollBody").get(0);
        scrollBody.scrollLeft = this.scrollLeft;
        $(scrollBody).trigger("scroll");
    });
//    console.timeEnd("scroll to top");

// =============================================================================
// START OVERLAY COLUMNS CODE HERE
// =============================================================================
    let toggle_enabled = false;
    $("#toggle_fixed_btn1").click(function() {
        if (toggle_enabled) {
            toggle_enabled = false;
            $("#overlay1").hide();
            $("#toggle_fixed_btn1").attr("value","Enable fixed columns");
            $("#toggle_fixed_btn1").addClass("clicked_button");
        } else {
            toggle_enabled = true;
            $("#familiestabletab").closest(".dataTables_scrollBody").scroll();
            $("#toggle_fixed_btn1").attr("value","Disable fixed columns");
            $("#toggle_fixed_btn1").removeClass("clicked_button");
        }
    });

    $("#toggle_fixed_btn2").click(function() {
        if (toggle_enabled) {
            toggle_enabled = false;
            $("#overlay2").hide();
            $("#toggle_fixed_btn2").attr("value","Enable fixed columns");
            $("#toggle_fixed_btn2").addClass("clicked_button");
        } else {
            toggle_enabled = true;
            $("#subtypestabletab").closest(".dataTables_scrollBody").scroll();
            $("#toggle_fixed_btn2").attr("value","Disable fixed columns");
            $("#toggle_fixed_btn2").removeClass("clicked_button");
        }
    });

// --------- overlay for table 1 ---------
    var left1 = 0;
    var old_left1 = 0;
    $("#familiestabletab").closest(".dataTables_scrollBody").scroll(function(){
        // If user scrolls and it's > 100px from left, then attach fixed columns overlay
        left1 = $("#familiestabletab").closest(".dataTables_scrollBody").scrollLeft();
        if (left1!==old_left1) {
            $("#overlay1").hide();
        }
        old_left1 = left1;

        if (left1 > 50 && toggle_enabled) {
            $("#overlay1").css({ left: left1 + "px" });
            if ($("#overlay1").is(":hidden")) {
                $("#overlay1").show();
            }
        }
    });

    $("#familiestabletab").closest(".dataTables_scrollBody").append('<div id="overlay1"><table id="overlay_table1" class="row-border text-center compact dataTable no-footer text-nowrap"><tbody></tbody></table></div>');

    function create_overlay1() {
        // This function fires upon filtering, to update what rows to show as an overlay
        $("#overlay_table1 tbody tr").remove();
        var $target = $("#overlay_table1 tbody");

        $("#familiestabletab tbody tr").each(function() {
            var $tds = $(this).children(),
                $row = $("<tr></tr>");
            $row.append($tds.eq(0).clone()).append($tds.eq(1).clone()).append($tds.eq(2).clone()).append($tds.eq(3).clone()).append($tds.eq(4).clone()).appendTo($target);
            $row.height($(this).height());
            //$row.font_size("10");
            //$row.height("31px");
            //$row.append($tds.height("2.5")).appendTo($target);
        });
        $("#overlay_table1 .border-right").removeClass("border-right");
    }

    // Function that detects filtering events
    $("#familiestabletab").on( "draw.dt", function (e, oSettings) {
        create_overlay1();
    });

    create_overlay1();
    $("#overlay1").hide();
// --------- End of overlay for table1 ---------


// ---------  start overlay for table2 ---------
    var left2 = 0;
    var old_left2 = 0;
    $("#subtypestabletab").closest(".dataTables_scrollBody").scroll(function(){
        // If user scrolls and it's > 100px from left, then attach fixed columns overlay
        left2 = $("#subtypestabletab").closest(".dataTables_scrollBody").scrollLeft();
        if (left2!==old_left2) {
            $("#overlay2").hide();
        }
        old_left2 = left2;

        if (left2 > 50 && toggle_enabled) {
            $("#overlay2").css({ left: left2 + "px" });
            if ($("#overlay2").is(":hidden")) {
                $("#overlay2").show();
            }
        }
    });

    $("#subtypestabletab").closest(".dataTables_scrollBody").append('<div id="overlay2"><table id="overlay_table2" class="row-border text-center compact dataTable no-footer text-nowrap"><tbody></tbody></table></div>');

    function create_overlay2() {
        // This function fires upon filtering, to update what rows to show as an overlay
        $("#overlay_table2 tbody tr").remove();
        var $target = $("#overlay_table2 tbody");

        $("#subtypestabletab tbody tr").each(function() {
            var $tds = $(this).children(),
                $row = $("<tr></tr>");
            $row.append($tds.eq(0).clone()).append($tds.eq(1).clone()).append($tds.eq(2).clone()).append($tds.eq(3).clone()).append($tds.eq(4).clone()).appendTo($target);
            $row.height($(this).height());
            //$row.font_size("10");
            //$row.height("31px");
            //$row.append($tds.height("2.5")).appendTo($target);
        });
        $("#overlay_table2 .rightborder").removeClass("rightborder");
    }

    // Function that detects filtering events
    $("#subtypestabletab").on( "draw.dt", function (e, oSettings) {
        create_overlay2();
    });

    create_overlay2();
    $("#overlay2").hide();
// End of overlay for table2
// =============================================================================
// END OVERLAY COLUMNS CODE HERE
// =============================================================================

// Gaspar's functions to copy to clipboard selected checkboxes as a newline separated list.
// copied from structure_browser.js and browser_functions.js. Notice that they depend on
// the jquery plugin PowerTip.js

    function copyToClipboard(array, delimiter, data_name, powertip_object=false) {
        var link = array;
//        console.log(link);
        var out = "";
        link.each(function() {
            var ele = $(this).attr("href").split("/");
            out+=ele[ele.length-1]+delimiter;
        });
        if (out.length===0) {
            window.alert("No entries selected for copying");
            return 0;
        }
        var textArea = document.createElement("textarea");
        textArea.value = out;
        document.body.appendChild(textArea);
        textArea.focus();
        textArea.select();
        try {
            var successful = document.execCommand("copy");
            var msg = successful ? "Successful" : "Unsuccessful";
            if (powertip_object!==false) {
                $.powerTip.hide();
                powertip_object.data("powertipjq", $([
                    "<p>"+array.length+" ID's copied to clipboard!</p>"
                ].join("\n")));
                powertip_object.powerTip("show");
                setTimeout(function() {
                    powertip_object.data("powertipjq", $([
                        "<p>Export "+data_name+"</p>"
                    ].join("\n")));
                },1000);
            }
        } catch (err) {
            window.alert("Oops, unable to copy");
        }
        document.body.removeChild(textArea);
    }

    $(".uniprot-export1").data("powertipjq", $([
        "<p>Export UniProt IDs</p>"
    ].join("\n")));

    $(".uniprot-export2").data("powertipjq", $([
        "<p>Export UniProt IDs</p>"
    ].join("\n")));

    $(".glyphicon-export").powerTip({
        placement: "n",
        smartPlacement: true
    });

    $("#uniprot_copy1").click(function () {
        copyToClipboard($(".alt_selected > .uniprot1 > a"), "\n", "UniProt IDs", $(".uniprot-export1"));
    });

    $("#uniprot_copy2").click(function () {
        copyToClipboard($(".alt_selected > .uniprot2 > a"), "\n", "UniProt IDs", $(".uniprot-export2"));
    });

//Uncheck every row when using back button on browser
    $(".alt_selected").prop("checked", false);
    $(".alt").prop("checked", false);
    $(".select-all").prop("checked", false);

});

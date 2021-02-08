/*global yadcf*/
/*eslint complexity: ["error", 8]*/
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
        console.log("columns variable " + columns);
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

// Draft for calculation of normalized rank for a given column
// NOTE: the support is currently not taken into account
function createRank(table_id, column) {
    // Step 1 - collect all values for a given column
    var min_max = [];
    $(table_id+" tbody tr td").filter(":nth-child("+column+")").each( function() {
        var cell_value = $(this).text();
        if (/^-?\d*(\.\d+)?$/.test(cell_value)){
            min_max.push(parseFloat(cell_value));
        }
    });

    // Step 2 - normalize all values and add them to a data attribute
    var min = Math.min(...min_max);
    var max = Math.max(...min_max);
    console.log(min_max, min, max);

    $(table_id+" tbody tr td").filter(":nth-child("+column+")").each( function() {
        var cell_value = $(this).text();
        if (/^-?\d*(\.\d+)?$/.test(cell_value)) {
            $(this).attr("data-normalized", Math.round((parseFloat(cell_value)-min)/(max-min)*100),0);
            $(this).attr("data-raw", cell_value);
        }
    });
}

// # use a place holder for data-normalized

// // # custom rankedRangeFilter draft for YADCF
// /**
//  * This is a custom YADCF function that checks ....
//  * ....
//  * @param {object} filterVal Value to filter on (not applicable)
//  * @param {object} columnVal Element in the filtered column
//  * @param {object} rowValues All elements in this row (not used)
//  * @param {object} stateVal Current DOM state of the row (not sufficient in this case)
//  * @returns {boolean} true if row contains selected target otherwise false
//  */
// function rankedRangeFilter(filterVal, columnVal, rowValues, stateVal){
//   // Check range or rank filter
//   if (filterVal.contains("range")){
//       var range_value = parseFloat($(columnVal).text());
//       var actual_filtering = parseFloat(filterVal.split("_")[2]);
//       if (filterVal.contains("min")){
//         // Example filterVal would be range_min_5
//         return range_value >= actual_filtering;
//       } else {
//         // Example filterVal would be range_max_10
//         return range_value <= actual_filtering;
//       }
//   } else {
//       var ranked_value = parseFloat($(columnVal).attr("data-normalized"));
//       // Example filterVal would be 15 (always a number)
//       return ranked_value<=filterVal;
//   }
// }
// // # YADCF setting for ranked range column
//  {
//                   column_number: X,
//                   filter_type: "custom_func",
//                   custom_func: rankedRangeFilter,
//                 }
// //# Add two input boxes above the automatically added filter
// //# Link boxes to filter by, for example:
// yadcf.exFilterColumn(targetTable, [[X, "range_min_5"]]);



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
function make_range_number_cols(start_column, repeat_number) {
    var from_to = {
        filter_type: "range_number",
        filter_default_label: ["Min", "Max"],
        filter_reset_button_text: false
    };
    var repeated_from_to = [];
    for (var i = start_column; i < start_column + repeat_number; i++) {
        var column_info = Object.assign({}, from_to);
        column_info["column_number"] = i;
        repeated_from_to.push(column_info);
    }
    return repeated_from_to;
}

$(document).ready(function() {
// Activate tooltips and popovers from Bootstrap   * Bootstrap v3.3.7 (http://getbootstrap.com)
    $("[data-toggle='tooltip']").tooltip();
    $("[data-toggle='popover']").popover();

//     $('a[data-toggle="tab"]').on("shown.bs.tab", function (e) {
// //      console.log( 'show tab' );
//         $($.fn.dataTable.tables(true)).DataTable()
//             .columns.adjust().responsive;
//     });

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
        order: [[4,"asc"], [22, "desc"]],
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
                column_number: [0],
                filter_type: "none",
                filter_default_label: "",
                filter_reset_button_text: false,
            },
            {
                column_number: [1],
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "80px",
                }
            },

            {
                column_number: [2],
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "40px",
                }
            },
            {
                column_number: [3],
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "200px",
                }
            },
            {
                column_number: [4],
                filter_type: "multi_select",
                select_type: "select2",
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "60px",
                }
            },
            {
                column_number: [5],
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
                column_number: [6],
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "40px"
                },
            },
            {
                column_number: [7],
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "40px"
                },
            },
            {
                column_number: [8],
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "40px"
                },
            },
            {
                column_number: [9],
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "40px"
                },
            },

// log(Emax/EC50)
            {
                column_number : [10],
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : [11],
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : [12],
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : [13],
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },

// pEC50
            {
                column_number : [14],
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : [15],
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : [16],
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : [17],
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },

// Emax
            {
                column_number : [18],
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : [19],
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : [20],
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : [21],
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },


// Hidden column calling a customized function
            {
                column_number: [22],
                filter_type: "custom_func",
                custom_func: supportFilter,
                filter_container_id: "hide_filter1",
            },

        ],

        {filters_tr_index: 2},

        {
            cumulative_filtering: true
        }
    );

    // Try intializing the rank
//    createRank("#familiestabletab", 11); // GS

    yadcf.exFilterColumn(oTable1, [[22, 2]]);

//    yadcf.exResetAllFilters(oTable1);
//    setTimeout(() => {
//        console.timeEnd("table1load");
//    }, );

    $("#familiestabletab"+" > tbody > tr").click(function(event) {
        if (event.target.type !== "checkbox") {
            $(":checkbox", this).trigger("click");
            $(this).eq(0).toggleClass("alt_selected");
            $(this).find("td").toggleClass("highlight");
        }
        $(this).eq(0).toggleClass("alt_selected");
        $(this).find("td").toggleClass("highlight");
    });

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

    console.timeEnd("table1load");

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

//    repeated_from_to_1 = make_range_number_cols(35, 14);

    yadcf.init(oTable2,
        [
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
                filter_default_label: "",
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

// log(Emax/EC50)
            {
                column_number : 10,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 11,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 12,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 13,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 14,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 15,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 16,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 17,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 18,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 19,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 20,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 21,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 22,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },

// pEC50
            {
                column_number : 23,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 24,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 25,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 26,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 27,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 28,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 29,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 30,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 31,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 32,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 33,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 34,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 35,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },

// Emax
            {
                column_number : 36,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 37,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 38,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 39,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 40,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 41,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 42,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 43,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 44,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 45,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 46,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 47,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 48,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },


            {
                column_number: 49,
                filter_type: "custom_func",
                custom_func: supportFilter,
                filter_container_id: "hide_filter2",
            },
        ],
//        ].concat(repeated_from_to_1),

        {
            filters_tr_index: 2
        },

        {
            cumulative_filtering: true
        }
    );

    yadcf.exFilterColumn(oTable2, [[49, 2]]);

    $("#subtypestabletab"+" > tbody > tr").click(function(event) {
        if (event.target.type !== "checkbox") {
            $(":checkbox", this).trigger("click");
            $(this).eq(0).toggleClass("alt_selected");
            $(this).find("td").toggleClass("highlight");
        }
        $(this).eq(0).toggleClass("alt_selected");
        $(this).find("td").toggleClass("highlight");
    });

//    yadcf.exResetAllFilters(oTable2);
//    setTimeout(() => {
//        console.timeEnd("table2load");
//    }, );
    console.timeEnd("table2load");


// By default display the first tab. If this is not ON, one has to click on the tab for display.
    $("#couplingtabs a:first").tab("show");
//    $('#couplingtabs a[href="#table_1"]').tab('show');

// Just a button to go back to the main page.
    $("#reset_tab1").click(function () {
        window.location.href = "/signprot/couplings2";
    });

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
    let toggle_enabled = true;
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
        console.log(link);
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
                    "<p>Copied to clipboard!</p>"
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

    $(".uniprot-export").data("powertipjq", $([
        "<p>Export UniProt IDs</p>"
    ].join("\n")));

    $(".glyphicon-export").powerTip({
        placement: "n",
        smartPlacement: true
    });

    $("#uniprot_copy").click(function () {
        copyToClipboard($(".alt_selected > .uniprot > a"), "\n", "UniProt IDs", $(".uniprot-export"));
    });

    //Uncheck every row when using back button on browser
    // $(".alt_selected").prop("checked",false);
    // $(".alt").prop("checked",false);
    // $(".select-all").prop("checked",false);


});

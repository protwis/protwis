/*eslint complexity: ["error", 8]*/
/*eslint quotes: ["error", "double", { "avoidEscape": true }]*/

let table1data;

var tableToExcel = (function() {
    var uri = "data:application/vnd.ms-excel;base64,",
        template = "<html xmlns:o=\"urn:schemas-microsoft-com:office:office\" xmlns:x=\"urn:schemas-microsoft-com:office:excel\" xmlns=\"http://www.w3.org/TR/REC-html40\"><head><!--[if gte mso 9]><xml><x:ExcelWorkbook><x:ExcelWorksheets><x:ExcelWorksheet><x:Name>{worksheet}</x:Name><x:WorksheetOptions><x:DisplayGridlines/></x:WorksheetOptions></x:ExcelWorksheet></x:ExcelWorksheets></x:ExcelWorkbook></xml><![endif]--></head><body><table>{table}</table></body></html>",
        base64 = function(s) {
            return window.btoa(unescape(encodeURIComponent(s)));
        },
        format = function(s, c) {
            return s.replace(/{(\w+)}/g, function(m, p) {
                return c[parseInt(p, 10)];
            });
        };
    return function(table, name, filename) {
        var table_obj = $("#" + table).clone();
        $("#excel_table").html(table_obj);
        // Clean up table to remove yadcf stuff
        $("#excel_table thead tr").css("height", "");
        $("#excel_table thead th").css("height", "");
        $("#excel_table thead div").css("height", "");
        $("#excel_table thead .yadcf-filter-wrapper").remove();
        $("#excel_table thead button").remove();
        var tr = $("#excel_table thead tr:eq(1)");
        // reattach th titles
        tr.find("th").each(function(column, th) {
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

function reset_tab1() {
// Just a button to go back to the main page.
    window.location.href = "/signprot/couplings2";
}

function reset_tab2() {
// Just a button to go back to the main page.
    window.location.href = "/signprot/couplings2";
}

//this.element.addEventListener(t, e, { passive: true} )
//$(document.addEventListener('touchstart', null, { passive: true})).ready(function () {
//$(document).ready(function () {
$(function() {
    let oTable1 = [];
    let oTable2 = [];

//     $('a[data-toggle="tab"]').on('shown.bs.tab', function (e) {
// //      console.log( 'show tab' );
//         $($.fn.dataTable.tables(true)).DataTable()
//             .columns.adjust().responsive.recalc();
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
        autoWidth: false,
        bInfo: true,
    });

    let yadcf1 = yadcf;
    yadcf1.init(oTable1,
        [
            {
                column_number: 0,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "Source",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                }
            },
            {
                column_number: 1,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "Cl",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '40px',
                }
            },
            {
                column_number: 2,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "Family",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '200px',
                }
            },
            {
                column_number: 3,
                filter_type: "multi_select",
                select_type: "select2",
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "uniprot",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '60px',
                }
            },
            {
                column_number: 4,
                filter_type: "multi_select",
                select_type: "select2",
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "IUPHAR",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                }
            },

// Guide to Pharmacology
            {
                column_number: 5,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "40px"
                },
            },
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

        ],

        {filters_tr_index: 2},

        {
            cumulative_filtering: true
        }
    );

    yadcf1.exResetAllFilters(oTable1);
//    setTimeout(() => {
//        console.timeEnd("table1load");
//    }, );
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
        autoWidth: false,
        bInfo: true,
    });

    yadcf1.init(oTable2,
        [
            {
                column_number: 0,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "Source",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                }
            },
            {
                column_number: 1,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "Cl",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '40px',
                }
            },
            {
                column_number: 2,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "Family",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '200px',
                }
            },
            {
                column_number: 3,
                filter_type: "multi_select",
                select_type: "select2",
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "uniprot",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '60px',
                }
            },
            {
                column_number: 4,
                filter_type: "multi_select",
                select_type: "select2",
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "IUPHAR",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                }
            },

// log(Emax/EC50)
            {
                column_number : 5,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 6,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 7,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 8,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 9,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
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

// pEC50
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

// Emax
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

        ],
        {filters_tr_index: 2},

        {
            cumulative_filtering: true
        }
    );

    yadcf1.exResetAllFilters(oTable2);
//    setTimeout(() => {
//        console.timeEnd("table2load");
//    }, );
    console.timeEnd("table2load");


// By default display the first tab. If this is not ON, one has to click on the tab for display.
    $("#myTab a:first").tab("show");
//    $('#myTab a[href="#table_1"]').tab('show');

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
    $("#toggle_fixed_btn").click(function() {
        if (toggle_enabled) {
            toggle_enabled = false;
            $("#overlay1").hide();
            $("#overlay2").hide();
            $("#toggle_fixed_btn").attr("value","Enable fixed columns");
            $("#toggle_fixed_btn").addClass("clicked_button");
        } else {
            toggle_enabled = true;
            $("#familiestabletab").closest(".dataTables_scrollBody").scroll();
            $("#subtypestabletab").closest(".dataTables_scrollBody").scroll();

            $("#toggle_fixed_btn").attr("value","Disable fixed columns");
            $("#toggle_fixed_btn").removeClass("clicked_button");
        }
    });

    var left1 = 0;
    var old_left1 = 0;
    $("#subtypestabletab").closest(".dataTables_scrollBody").scroll(function(){
        // If user scrolls and it's > 100px from left, then attach fixed columns overlay
        left1 = $("#subtypestabletab").closest(".dataTables_scrollBody").scrollLeft();
        if (left1!==old_left1) {
            $("#overlay2").hide();
        }
        old_left1 = left1;

        if (left1 > 100 && toggle_enabled) {
            $("#overlay2").css({ left: left1 + "px" });
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


    create_overlay2();
    $("#overlay2").hide();
// =============================================================================
// END OVERLAY COLUMNS CODE HERE
// =============================================================================

});

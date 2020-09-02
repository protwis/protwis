let oTable1 = [];
let oTable2 = [];
let oTable3 = [];
let oTable4 = [];

let table1data;

var tableToExcel = (function () {
    var uri = 'data:application/vnd.ms-excel;base64,',
        template = '<html xmlns:o="urn:schemas-microsoft-com:office:office" xmlns:x="urn:schemas-microsoft-com:office:excel" xmlns="http://www.w3.org/TR/REC-html40"><head><!--[if gte mso 9]><xml><x:ExcelWorkbook><x:ExcelWorksheets><x:ExcelWorksheet><x:Name>{worksheet}</x:Name><x:WorksheetOptions><x:DisplayGridlines/></x:WorksheetOptions></x:ExcelWorksheet></x:ExcelWorksheets></x:ExcelWorkbook></xml><![endif]--></head><body><table>{table}</table></body></html>',
        base64 = function (s) {
            return window.btoa(unescape(encodeURIComponent(s)))
        }, format = function (s, c) {
            return s.replace(/{(\w+)}/g, function (m, p) {
                return c[p];
            })
        }
    return function (table, name, filename) {
        var table= $("#"+table).clone();
        $("#excel_table").html(table);

        // Clean up table to remove yadcf stuff
        $("#excel_table thead tr").css('height','');
        $("#excel_table thead th").css('height','');
        $("#excel_table thead div").css('height','');
        $("#excel_table thead .yadcf-filter-wrapper").remove();
        $("#excel_table thead button").remove();
        var tr = $("#excel_table thead tr:eq(1)");

        // reattach th titles
        tr.find('th').each (function( column, th) {
            if ($(th).attr('title')) $(th).html($(th).attr('title'));
        });

        var ctx = {
            worksheet: name || 'Worksheet',
            table: $("#excel_table").html()
        }
        $("#excel_table").html("");
        document.getElementById("dlink").href = uri + base64(format(template, ctx));
        document.getElementById("dlink").download = filename;
        document.getElementById("dlink").click();
    }
})()

function select_all(e) {
    var checkedStatus = $(e).prop("checked");

    $('.select-all  ').each(function () {
        $(this).prop('checked', checkedStatus);
    });

    $('.alt').each(function () {
        $(this).prop('checked', checkedStatus);
    });
};

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
        console.log('columns variable ' + columns);
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

function resetHidden3() {
    var columns = Array.from(new Array(70), (x,i) => i + 3);
    columns.forEach(function(column) {
        console.log('columns variable ' + columns);
        column = oTable3.column( column );
        try {
            column.visible( true, false );
        }
        catch(err) {
            column.visible( true, false );
        }
    });
    oTable3.draw();
}

function resetHidden4() {
    var columns = Array.from(new Array(70), (x,i) => i + 3);
    columns.forEach(function(column) {
        console.log('columns variable ' + columns);
        column = oTable4.column( column );
        try {
            column.visible( true, false );
        }
        catch(err) {
            column.visible( true, false );
        }
    });
    oTable4.draw();
}

function reset_tab1() {
// Just a button to go back to the main page.
    window.location.href = '/signprot/couplings#table_1';
}

function reset_tab2() {
// Just a button to go back to the main page.
    window.location.href = '/signprot/couplings#table_2';
}

//this.element.addEventListener(t, e, { passive: true} )
//$(document.addEventListener('touchstart', null, { passive: true})).ready(function () {
$(document).ready(function () {

    $('a[data-toggle="tab"]').on('shown.bs.tab', function (e) {
//      console.log( 'show tab' );
        $($.fn.dataTable.tables(true)).DataTable()
            .columns.adjust().responsive.recalc();
    });

    console.time("table1load");
    oTable1 = $("#familiestabletab").DataTable({
//        data: table1data,
//        serverSide: true,
        deferRender: true,
        scrollY: '50vh',
        scrollX: true,
        scrollCollapse: true,
        scroller: true,
        paging: false,
//        lengthMenu: [[10, 25, 50, -1], [10, 25, 50, "All"]],
        bSortCellsTop: false, //prevent sort arrows going on bottom row
        aaSorting: [],
        autoWidth: false,
//        pageLength: -1,
        bInfo: true,
        columnDefs: [
            {
                targets: [6,8,10,12],
                visible: false
            }
        ],
        // columns: [
        //     null,
        //     null,
        //     null,
        //     null,
        //     null, // 4
        //     {width: "20%"}, // 5
        //     null, // 6
        //     null, // 7
        //     null, // 8
        //     null, // 9
        //     null, // 10
        //     null, // 11
        //     null,
        //     null,
        //     null,
        //     null,
        //     null,
        //     null,
        //     null,
        //     null,
        //     null,
        //     null,
        //     null,
        //     null,
        //     null,
        //     null
        // ],
    });

    yadcf.init(oTable1,
        [
            {
                column_number: 1,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Class",
                filter_reset_button_text: false,
            },
            {
                column_number: 2,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Family",
                filter_reset_button_text: false,
            },
            {
                column_number: 3,
                filter_type: "multi_select",
                select_type: 'select2',
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "uniprot",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
            },
            {
                column_number: 4,
                filter_type: "multi_select",
                select_type: 'select2',
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "IUPHAR",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
            },

            {
                column_number : 5,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },

            {
                column_number: 6,
                filter_type: "multi_select",
                filter_container_id: "gs_drop",
                select_type: 'select2',
                filter_default_label: "Gs",
                select_type_options: {
                    width: '55%'
                },
            },

            {
                column_number : 7,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },

            {
                column_number: 8,
                filter_type: "multi_select",
                filter_container_id: "gio_drop",
                select_type: 'select2',
                filter_default_label: "Gi/Go",
                select_type_options: {
                    width: '55%'
                },
            },

            {
                column_number : 9,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },

            {
                column_number: 10,
                filter_type: "multi_select",
                filter_container_id: "gq11_drop",
                select_type: 'select2',
                filter_default_label: "Gq/G11",
                select_type_options: {
                    width: '55%'
                },
            },

            {
                column_number : 11,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },

            {
                column_number: 12,
                filter_type: "multi_select",
                filter_container_id: "g1213_drop",
                select_type: 'select2',
                filter_default_label: "G12/13",
                select_type_options: {
                    width: '55%'
                },
            },

            {
                column_number : 13,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },

            {
                column_number : 14,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },

            {
                column_number : 15,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },

            {
                column_number : 16,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },

            {
                column_number : 17,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },

            {
                column_number : 18,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },

            {
                column_number :19,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },

            {
                column_number : 20,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },

            {
                column_number : 21,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },

            {
                column_number : 22,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },

            {
                column_number : 23,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },

            {
                column_number : 24,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },

            {
                column_number : 25,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },

        ],

        {filters_tr_index: 1},

        {
            cumulative_filtering: true
        }
    );

    yadcf.exResetAllFilters(oTable1);
//    setTimeout(() => {
//        console.timeEnd("table1load");
//    }, );
    console.timeEnd("table1load");

    console.time("table2load");
    oTable2 = $("#subtypestabletab").DataTable({
//        data: table2data,
//        serverSide: true,
        deferRender: true,
        scrollY: '50vh',
        scrollX: true,
        scrollCollapse: true,
        scroller: true,
        paging: false,
//        lengthMenu: [[10, 25, 50, -1], [10, 25, 50, "All"]],
        bSortCellsTop: false, //prevent sort arrows going on bottom row
        aaSorting: [],
        autoWidth: false,
//        fixedColumns:   {
//            heightMatch: 'none'
//        },
//        order: [[4, 'asc'],[26,'desc']],
//        columnDefs: [
//            {
//                targets: 'no-sort',
//                orderable: false
//            }
//            ],
        pageLength: -1,
        bInfo: true
    });
    yadcf.init(oTable2,
        [
            {
                column_number: 1,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Class",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '90px',
                },
            },
            {
                column_number: 2,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Family",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '100px',
                },
            },
            {
                column_number: 3,
                filter_type: "multi_select",
                select_type: 'select2',
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "UniProt",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 4,
                filter_type: "multi_select",
                select_type: 'select2',
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "IUPHAR",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '100px',
                },
            },

            {
                column_number : 5,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 6,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 7,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },


            {
                column_number : 8,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 9,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 10,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 11,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 12,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 13,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 14,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 15,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 16,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },


            {
                column_number : 17,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 18,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 19,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 20,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 21,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 22,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 23,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 24,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },


            {
                column_number : 25,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 26,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 27,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 28,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },


            {
                column_number : 29,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 30,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
            },

        ],
        {filters_tr_index: 2},

        {
            cumulative_filtering: false
        }
    );

    yadcf.exResetAllFilters(oTable2);
//    setTimeout(() => {
//        console.timeEnd("table2load");
//    }, );
    console.timeEnd("table2load");

    console.time("table3load");
    oTable3 = $("#bouviertabletab").DataTable({
//        data: table3data,
        deferRender: true,
        scrollY:  '50vh',
        scrollX: true,
        scrollCollapse: true,
        scroller: true,
        paging: false,
//        lengthMenu: [[10, 25, 50, -1], [10, 25, 50, "All"]],
        bSortCellsTop: false, //prevent sort arrows going on bottom row
        aaSorting: [],
        autoWidth: false,
//        pageLength: -1,
        bInfo: true
    });
    yadcf.init(oTable3,
        [
            {
                column_number: 1,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Class",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '90px',
                },
            },
            {
                column_number: 2,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Family",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '100px',
                },
            },
            {
                column_number: 3,
                filter_type: "multi_select",
                select_type: 'select2',
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "UniProt",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 4,
                filter_type: "multi_select",
                select_type: 'select2',
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "IUPHAR",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '100px',
                },
            },

            {
                column_number: 5,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 6,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 7,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 8,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 9,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '80px',
                },
            },

            {
                column_number: 10,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 11,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 12,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 13,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 14,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 15,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 16,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 17,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 18,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 19,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 20,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },


            {
                column_number: 21,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 22,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 23,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 24,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 25,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 26,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 27,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 28,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 29,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 30,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 31,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },


            {
                column_number: 32,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 33,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 34,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 35,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 36,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 37,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 38,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 39,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 40,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 41,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 42,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },


            {
                column_number: 43,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 44,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 45,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 46,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 47,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 48,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 49,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 50,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 51,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 52,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 53,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },


            {
                column_number: 54,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 55,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 56,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 57,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 58,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 59,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 60,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 61,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 62,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 63,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 64,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },


            {
                column_number: 65,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 66,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 67,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 68,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 69,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 70,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 71,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 72,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 73,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 74,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 75,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },


        ],
        {filters_tr_index: 1},
        {
            cumulative_filtering: true
        }
    );

    yadcf.exResetAllFilters(oTable3);
//    setTimeout(() => {
//        console.timeEnd("table3load");
//    }, );
    console.timeEnd("table3load");

    console.time("table4load");
    oTable4 = $("#inouetabletab").DataTable({
//        data: table4data,
        //data: data1,
        //"ordering": false,
        //"fixedHeader":true,
        //"searching": false,
        //"info": false,
        deferRender: true,
        scrollY:  '50vh',
        scrollX: true,
//        scrollCollapse: true,
        scroller: true,
        paging: false,
        lengthMenu: [[10, 25, 50, -1], [10, 25, 50, "All"]],
        bSortCellsTop: false, //prevent sort arrows going on bottom row
        aaSorting: [],
        autoWidth: false,
//        pageLength: -1,
        bInfo: true
    });

    yadcf.init(oTable4,
        [
            {
                column_number: 1,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Class",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '90px',
                },
            },
            {
                column_number: 2,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Family",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '100px',
                },
            },
            {
                column_number: 3,
                filter_type: "multi_select",
                select_type: 'select2',
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "UniProt",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 4,
                filter_type: "multi_select",
                select_type: 'select2',
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "IUPHAR",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '100px',
                },
            },


            {
                column_number: 5,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 6,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 7,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 8,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '80px',
                },
            },


            {
                column_number: 9,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 10,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 11,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 12,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 13,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 14,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 15,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 16,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 17,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 18,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 19,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },


            {
                column_number: 20,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 21,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 22,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 23,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 24,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 25,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 26,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 27,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 28,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 29,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 30,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },


            {
                column_number: 31,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 32,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 33,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 34,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 35,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 36,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 37,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 38,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 39,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 40,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 41,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },


            {
                column_number: 42,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 43,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 44,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 45,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 46,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 47,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 48,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 49,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 50,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 51,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 52,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },


            {
                column_number: 53,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 54,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 55,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 56,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 57,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 58,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 59,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 60,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 61,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 62,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 63,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },


            {
                column_number: 64,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 65,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 66,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 67,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 68,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 69,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 70,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 71,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 72,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 73,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },
            {
                column_number: 74,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                select_type_options: {
                    width: '70px',
                },
            },



        ],

        {filters_tr_index: 1},

        {
            cumulative_filtering: true
        }
    );

    yadcf.exResetAllFilters(oTable4);

//    setTimeout(() => {
//        console.timeEnd("table4load");
//    }, );
    console.timeEnd("table4load");

// By default display the first tab. If this is not ON, one has to click on the tab for display.
    $('#myTab a:first').tab('show');
//    $('#myTab a[href="#table_1"]').tab('show');

// Just a button to go back to the main page.
    $('#reset_tab1').click(function () {
        window.location.href = '/signprot/couplings';
    });

    $('.hide_columns1').click(function(evt) {
        var columns = $(this).attr('columns').split(",");
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

    $('.hide_columns2').click(function(evt) {
        var columns = $(this).attr('columns').split(",");
        columns.forEach(function(column) {
            var column = oTable2.column( column );
            try {
                column.visible( false, false );
            }
            catch(err) {
                column.visible( false, false );
            }
        });
        oTable2.draw();
    } );

    $('.hide_columns3').click(function(evt) {
        var columns = $(this).attr('columns').split(",");
        columns.forEach(function(column) {
            var column = oTable3.column( column );
            try {
                column.visible( false, false );
            }
            catch(err) {
                column.visible( false, false );
            }
        });
        oTable3.draw();
    } );

    $('.hide_columns4').click(function(evt) {
        var columns = $(this).attr('columns').split(",");
        columns.forEach(function(column) {
            var column = oTable4.column( column );
            try {
                column.visible( false, false );
            }
            catch(err) {
                column.visible( false, false );
            }
        });
        oTable4.draw();
    } );


// Put top scroller
// https://stackoverflow.com/questions/35147038/how-to-place-the-datatables-horizontal-scrollbar-on-top-of-the-table
//    console.time("scroll to top");
    $('.dataTables_scrollHead').css({
        'overflow-x':'scroll'
    }).on('scroll', function(e){
        var scrollBody = $(this).parent().find('.dataTables_scrollBody').get(0);
        scrollBody.scrollLeft = this.scrollLeft;
        $(scrollBody).trigger('scroll');
    });
//    console.timeEnd("scroll to top");



// =============================================================================
// START OVERLAY COLUMNS CODE HERE
// =============================================================================
    toggle_enabled = true;
    $('#toggle_fixed_btn').click(function() {
        if (toggle_enabled) {
            toggle_enabled = false;
            $("#overlay1").hide();
            $("#overlay2").hide();
            $("#overlay3").hide();
            $("#toggle_fixed_btn").attr('value',"Enable fixed columns");
            $("#toggle_fixed_btn").addClass('clicked_button');
        } else {
            toggle_enabled = true;
            $('#subtypestabletab').closest('.dataTables_scrollBody').scroll();
            $('#bouviertabletab').closest('.dataTables_scrollBody').scroll();
            $('#inouetabletab').closest('.dataTables_scrollBody').scroll();
            $("#toggle_fixed_btn").attr('value',"Disable fixed columns");
            $("#toggle_fixed_btn").removeClass('clicked_button');
        }
    });

    var left1 = 0;
    var old_left1 = 0;
    $('#subtypestabletab').closest('.dataTables_scrollBody').scroll(function(){
        // If user scrolls and it's > 100px from left, then attach fixed columns overlay
        left1 = $('#subtypestabletab').closest('.dataTables_scrollBody').scrollLeft();
        if (left1!=old_left1) $("#overlay1").hide();
        old_left1 = left1;

        if (left1 > 100 && toggle_enabled) {
            $("#overlay1").css({ left: left1 + 'px' });
            if ($("#overlay1").is(":hidden")) $("#overlay1").show();
        }
    });

    var left2 = 0;
    var old_left2 = 0;
    $('#bouviertabletab').closest('.dataTables_scrollBody').scroll(function(){
        // If user scrolls and it's > 100px from left, then attach fixed columns overlay
        left2 = $('#bouviertabletab').closest('.dataTables_scrollBody').scrollLeft();
        if (left2!=old_left2) $("#overlay2").hide();
        old_left2 = left2;

        if (left2 > 100 && toggle_enabled) {
            $("#overlay2").css({ left: left2 + 'px' });
            if ($("#overlay2").is(":hidden")) $("#overlay2").show();
        }
    });

    var left3 = 0;
    var old_left3 = 0;
    $('#inouetabletab').closest('.dataTables_scrollBody').scroll(function(){
        // If user scrolls and it's > 100px from left, then attach fixed columns overlay
        left3 = $('#inouetabletab').closest('.dataTables_scrollBody').scrollLeft();
        if (left3!=old_left3) $("#overlay3").hide();
        old_left3 = left3;

        if (left3 > 100 && toggle_enabled) {
            $("#overlay3").css({ left: left3 + 'px' });
            if ($("#overlay3").is(":hidden")) $("#overlay3").show();
        }
    });

    $('#subtypestabletab').closest('.dataTables_scrollBody').append('<div id="overlay1"><table' +
        ' class="row-border text-center compact dataTable no-footer text-nowrap" ' +
        ' id="overlay_table1"><tbody></tbody></table></div>');

    $('#bouviertabletab').closest('.dataTables_scrollBody').append('<div id="overlay2"><table' +
        ' class="row-border text-center compact dataTable no-footer text-nowrap" ' +
        ' id="overlay_table2"><tbody></tbody></table></div>');

    $('#inouetabletab').closest('.dataTables_scrollBody').append('<div id="overlay3"><table' +
        ' class="row-border text-center compact dataTable no-footer text-nowrap" ' +
        ' id="overlay_table3"><tbody></tbody></table></div>');

    function create_overlay1() {
        // This function fires upon filtering, to update what rows to show as an overlay
        $("#overlay_table1 tbody tr").remove();
        var $target = $("#overlay_table1 tbody");

        $("#subtypestabletab tbody tr").each(function() {
            var $tds = $(this).children(),
                $row = $("<tr></tr>");
            $row.append($tds.eq(1).clone()).append($tds.eq(2).clone()).append($tds.eq(3).clone()).append($tds.eq(4).clone()).appendTo($target);
            $row.height($(this).height());
            //$row.font_size("10");
            //$row.height("31px");
            //$row.append($tds.height("2.5")).appendTo($target);
        });
        $("#overlay_table1 .rightborder").removeClass("rightborder");
    }

    function create_overlay2() {
        // This function fires upon filtering, to update what rows to show as an overlay
        $("#overlay_table2 tbody tr").remove();
        var $target = $("#overlay_table2 tbody");

        $("#bouviertabletab tbody tr").each(function() {
            var $tds = $(this).children(),
                $row = $("<tr></tr>");
            $row.append($tds.eq(1).clone()).append($tds.eq(2).clone()).append($tds.eq(3).clone()).append($tds.eq(4).clone()).appendTo($target);
            $row.height($(this).height());
        });
        $("#overlay_table2 .rightborder").removeClass("rightborder");
    }

    function create_overlay3() {
        // This function fires upon filtering, to update what rows to show as an overlay
        $("#overlay_table3 tbody tr").remove();
        var $target = $("#overlay_table3 tbody");

        $("#inouetabletab tbody tr").each(function() {
            var $tds = $(this).children(),
                $row = $("<tr></tr>");
            $row.append($tds.eq(1).clone()).append($tds.eq(2).clone()).append($tds.eq(3).clone()).append($tds.eq(4).clone()).appendTo($target);
            $row.height($(this).height());
        });
        $("#overlay_table3 .rightborder").removeClass("rightborder");
    }

    create_overlay1();
    create_overlay2();
    create_overlay3();
    $("#overlay1").hide();
    $("#overlay2").hide();
    $("#overlay3").hide();
// =============================================================================
// END OVERLAY COLUMNS CODE HERE
// =============================================================================

});

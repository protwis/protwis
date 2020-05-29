var oTable1 = [];
var oTable2 = [];
var oTable3 = [];
var oTable4 = [];

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

$(document).ready(function () {
    $('a[data-toggle="tab"]').on('shown.bs.tab', function (e) {
        console.log( 'show tab' );
        $($.fn.dataTable.tables(true)).DataTable()
            .columns.adjust()
//            .responsive.recalc();
    });

    console.time("table1load");
    oTable1 = $("#familiestabletab").DataTable({
        "scrollY": $(window).height() - 450,
        "scrollX": true,
        "scrollCollapse": true,
        "scroller": true,
        "paging": false,
//        "lengthMenu": [[10, 25, 50, -1], [10, 25, 50, "All"]],
        "bSortCellsTop": false, //prevent sort arrows going on bottom row
        "aaSorting": [],
        "autoWidth": true,
        "pageLength": -1,
        "bInfo": true,
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
                select_type_options: {
                    width: '95px',
                },
            },
            {
                column_number: 5,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Bouvier",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '60px',
                },
            },
            {
                column_number: 6,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Inoue",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '60px',
                },
            },
            {
                column_number: 7,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gtp",
                filter_reset_button_text: false,
            },
            {
                column_number: 8,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Bouvier",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '60px',
                },
            },
            {
                column_number: 9,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Inoue",
                filter_reset_button_text: false,
            },
            {
                column_number: 10,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GtP",
                filter_reset_button_text: false,
            },
            {
                column_number: 11,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Bouvier",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '60px',
                },
            },
            {
                column_number: 12,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Inoue",
                filter_reset_button_text: false,
            },
            {
                column_number: 13,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GtP",
                filter_reset_button_text: false,
            },
            {
                column_number: 14,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Bouvier",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '60px',
                },
            },
            {
                column_number: 15,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Inoue",
                filter_reset_button_text: false,
            },
            {
                column_number: 16,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gtp",
                filter_reset_button_text: false,
            },
            {
                column_number: 17,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Bouvier",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '100px',
                },
            },
        ],
        {filters_tr_index: 1},
        {
            cumulative_filtering: true
        }
    );

//    yadcf.exResetAllFilters(oTable1);
    setTimeout(() => {
        console.timeEnd("table1load");
    }, );

    console.time("table2load");
    oTable2 = $("#subtypestabletab").DataTable({
        "scrollY": $(window).height() - 450,
        "scrollX": true,
        "scrollCollapse": true,
        "scroller": true,
        "paging": false,
        "bSortCellsTop": false, //prevent sort arrows going on bottom row
        "autoWidth": true,
        "pageLength": -1,
        "bInfo": true,
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
                filter_default_label: "UniProt",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 4,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "IUPHAR",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '100px',
                },
            },
            {
                column_number: 5,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GNAS2",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 6,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GNAS2",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 7,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GNAL",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 8,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GNAI1",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 9,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GNAI2",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 10,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GNAO",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 11,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GNAZ",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 12,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GNAI1",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 13,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GNAI3",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 14,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GNAO",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 15,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GNAZ",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 16,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GNAQ",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 17,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GNA11",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 18,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GNA14",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 19,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GNA15",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 20,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GNAQ",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 21,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GNA14",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 22,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GNA15",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 23,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GNA12",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 24,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GNA13",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 25,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GNA12",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 26,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GNA13",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 27,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "BARR1/GRK2",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 28,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "BARR2/GRK2",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
        ],
        {filters_tr_index: 2},
        {
            cumulative_filtering: true
        }
    );

//    yadcf.exResetAllFilters(oTable2);
    setTimeout(() => {
        console.timeEnd("table2load");
    }, );

    console.time("table3load");
    oTable3 = $("#bouviertabletab").DataTable({
        "scrollY": $(window).height() - 450,
        "scrollX": true,
        "scrollCollapse": true,
        "scroller": true,
        "paging": false,
        "bSortCellsTop": false, //prevent sort arrows going on bottom row
        "autoWidth": true,
        "pageLength": -1,
        "bInfo": true,
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
                filter_default_label: "UniProt",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 4,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "IUPHAR",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '100px',
                },
            },
            {
                column_number: 5,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gs",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 6,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi/o",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 7,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "G12/13",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 8,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gq/11",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 9,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Arrestin",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 10,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gs",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 11,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi1",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 12,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi2",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 13,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GoA",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 14,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gs",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 15,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi1",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 16,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi2",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 17,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GoA",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 18,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gs",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 19,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi1",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 20,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi2",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 21,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GoA",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 22,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gs",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 23,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi1",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 24,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi2",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 25,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GoA",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 26,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gs",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 27,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi1",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 28,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi2",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 29,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GoA",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 30,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gs",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 31,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi1",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 32,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi2",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 33,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GoA",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 34,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gs",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 35,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi1",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 36,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi2",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 37,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GoA",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
        ],
        {filters_tr_index: 1},
        {
            cumulative_filtering: true
        }
    );

//    yadcf.exResetAllFilters(oTable3);
    setTimeout(() => {
        console.timeEnd("table3load");
    }, );

    console.time("table4load");
    oTable4 = $("#inouetabletab").DataTable({
        //data: data1,
        //"ordering": false,
        //"fixedHeader":true,
        //"searching": false,
        //"info": false,
        "deferRender": true,
        "scrollY": $(window).height() - 450,
        "scrollX": true,
        "scrollCollapse": true,
        "scroller": true,
        "paging": false,
        "bSortCellsTop": false, //prevent sort arrows going on bottom row
        "autoWidth": true,
        "pageLength": -1,
        "bInfo": true,
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
                filter_default_label: "UniProt",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 4,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "IUPHAR",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '100px',
                },
            },
            {
                column_number: 5,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gs",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 6,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi/o",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 7,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "G12/13",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 8,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gq/11",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                },
            },
            {
                column_number: 9,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gs",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 10,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi1",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 11,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi2",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 12,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GoA",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 13,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gs",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 14,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi1",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 15,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi2",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 16,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GoA",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 17,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gs",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 18,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi1",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 19,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi2",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 20,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GoA",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 21,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gs",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 22,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi1",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 23,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi2",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 24,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GoA",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 25,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gs",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 26,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi1",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 27,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi2",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 28,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GoA",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 29,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gs",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 30,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi1",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 31,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi2",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 32,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GoA",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 33,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gs",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 34,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi1",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 35,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gi2",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
            {
                column_number: 36,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GoA",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '30px',
                },
            },
        ],
        {filters_tr_index: 1},
        {
            cumulative_filtering: true
        }
    );

//    yadcf.exResetAllFilters(oTable4);
    setTimeout(() => {
        console.timeEnd("table4load");
    }, );

    // By default display the first tab. If this is not ON, one has to click on the tab for display.
    $('#myTab a:first').tab('show');

// Just a button to go back to the main page.
    $('#reset_tab1').click(function () {
        window.location.href = '/signprot/couplings/#table_1';
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

});

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
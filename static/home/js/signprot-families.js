$(function () {

    var oTable1 = $('#familiesDT').DataTable({
        "scrollY": $(window).height() - 450,
        "scrollX": true,
        "scrollCollapse": true,
        "scroller": true,
        "paging": false,
        "bSortCellsTop": true, //prevent sort arrows going on bottom row
        "aaSorting": [],
        "autoWidth": true,
        "lengthMenu": [[10, 25, 50, -1], [10, 25, 50, "All"]],
        "pageLength": -1,
        //"dom": 'Blfrtip',
        //"buttons": ['copy', 'csv', 'excel'],
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
                filter_default_label: "uniprot",
                filter_reset_button_text: false,
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
                filter_default_label: "Arrestin",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '60px',
                },
            },
        ],
        {filters_tr_index: 1},
        {
            cumulative_filtering: true
        }
    );

    yadcf.exResetAllFilters(oTable1);


    // By default display the first tab. If this is not on one has to click on the tab for display.
    $('#myTab a:first').tab('show');

    $('#familiesDT' + ' > tbody > tr').click(function (event) {
        if (event.target.type !== 'checkbox') {
            $(':checkbox', this).trigger('click');
            $(this).eq(0).toggleClass('alt_selected');
            $(this).find('td').toggleClass('highlight');
        }
        $(this).eq(0).toggleClass('alt_selected');
        $(this).find('td').toggleClass('highlight');
    });

    $(".select-all").click(function () {
        $(':checkbox', this).trigger('click');
        if ($(this).prop('checked') === true) {
            $('.alt').prop('checked', true);
            $('.alt').parent().parent().addClass('alt_selected');
            $('.alt').parent().parent().find('td').addClass('highlight');
        }
        if ($(this).prop('checked') === false) {
            $('.alt').prop('checked', false);
            $('.alt').parent().parent().removeClass('alt_selected');
            $('.alt').parent().parent().find('td').removeClass('highlight');
        }
    });

    // $('.wrapper').find('div').width($(".yadcf-datatables-table--familiesDT").width());

    $('.hide_columns').click(function (evt) {
        var columns = $(this).attr('columns').split(",");
        columns.forEach(function (column) {
            var column = oTable1.column(column);
            try {
                column.visible(false, false);
            } catch (err) {
                column.visible(false, false);
            }
        });
        oTable1.draw();
    });

    toggle_enabled = true;
    $('#toggle_fixed_btn').click(function () {
        if (toggle_enabled) {
            toggle_enabled = false;
            $("#overlay").hide();
            $("#toggle_fixed_btn").attr('value', "Enable fixed columns");
            $("#toggle_fixed_btn").addClass('clicked_button');
        } else {
            toggle_enabled = true;
            $('.dataTables_scrollBody').scroll();
            $("#toggle_fixed_btn").attr('value', "Disable fixed columns");
            $("#toggle_fixed_btn").removeClass('clicked_button');
        }
    });

    $("#toggle_columns_btn").click(function () {
        var columns = Array.from(new Array(24), (x, i) => i + 6);
        columns.forEach(function (column) {
            var column = oTable1.column(column);
            try {
                column.visible(true, false);
            } catch (err) {
                column.visible(true, false);
            }
        });
        oTable1.draw();
    });

    $('#uniprot_copy').click(function () {
        copyToClipboard('uniprot');
    });

    $('#reset_filters_btn').click(function () {
        window.location.href = '/signprot/couplings'
    });

    $('.dataTables_scrollBody').append('<div id=overlay><table id="overlay_table" ' +
        'class="row-border text-center compact dataTable no-footer text-nowrap"><tbody></tbody></table></div>');

    function create_overlay() {
        // This function fires upon filtering, to update what rows to show as an overlay
        $("#overlay_table tbody tr").remove();
        var $target = $("#overlay_table tbody");
        $("#familiesDT tbody tr").each(function () {
            var $tds = $(this).children(),
                $row = $("<tr></tr>");
            // $row.append($tds.eq(0).clone()).append($tds.eq(1).clone()).appendTo($target);
            $row.append($tds.eq(1).clone()).append($tds.eq(2).clone()).appendTo($target);
            $row.height($(this).height());
        });
        $("#overlay_table .border-right").removeClass("border-right");
    }

    // Function that detects filtering events
    $('#familiesDT').on('draw.dt', function (e, oSettings) {
        create_overlay();
    });

    create_overlay();
    $("#overlay").hide();

    var left = 0;
    var old_left = 0;
    toggle_enabled = true;
    $('.dataTables_scrollBody').scroll(function () {
        // If user scrolls and it's >100px from left, then attach fixed columns overlay
        left = $('.dataTables_scrollBody').scrollLeft();
        if (left != old_left) $("#overlay").hide();
        old_left = left;

        if (left > 100 && toggle_enabled) {
            $("#overlay").css({left: left + 'px'});
            if ($("#overlay").is(":hidden")) $("#overlay").show();
        }
    });
    // console.log($('#yadcf-filter--familiesDT-from-12').width());
    // $('#yadcf-filter--familiesDT-from-12').width(10);
    // console.log($('#yadcf-filter--familiesDT-from-12').width());

});



var tableToExcel = (function () {
    var uri = 'data:application/vnd.ms-excel;base64,',
        template = '<html xmlns:o="urn:schemas-microsoft-com:office:office" ' +
            'xmlns:x="urn:schemas-microsoft-com:office:excel" xmlns="http://www.w3.org/TR/REC-html40"><head><!--[if gte mso 9]><xml><x:ExcelWorkbook><x:ExcelWorksheets><x:ExcelWorksheet><x:Name>{worksheet}</x:Name><x:WorksheetOptions><x:DisplayGridlines/></x:WorksheetOptions></x:ExcelWorksheet></x:ExcelWorksheets></x:ExcelWorkbook></xml><![endif]--></head><body><table>{table}</table></body></html>',
        base64 = function (s) {
            return window.btoa(unescape(encodeURIComponent(s)))
        }, format = function (s, c) {
            return s.replace(/{(\w+)}/g, function (m, p) {
                return c[p];
            })
        }
    return function (table, name, filename) {
        var table = $("#" + table).clone();
        $("#excel_table").html(table);
        // Clean up table to remove yadcf stuff
        $("#excel_table thead tr").css('height', '');
        $("#excel_table thead th").css('height', '');
        $("#excel_table thead div").css('height', '');
        $("#excel_table thead .yadcf-filter-wrapper").remove();
        $("#excel_table thead button").remove();
        var tr = $("#excel_table thead tr:eq(1)");
        // reattach th titles
        tr.find('th').each(function (column, th) {
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

function copyDropdown() {
    document.getElementById("Dropdown").classList.toggle("show");
}

window.onclick = function (event) {
    if (!event.target.matches('.dropbtn')) {
        var dropdowns = document.getElementsByClassName("dropdown-content");
        var i;
        for (i = 0; i < dropdowns.length; i++) {
            var openDropdown = dropdowns[i];
            if (openDropdown.classList.contains('show')) {
                openDropdown.classList.remove('show');
            }
        }
    }
}

function copyToClipboard(data) {
    var link = $('.alt_selected > .' + data + ' > a');
    var out = [];
    link.each(function () {
        var ele = $(this).attr('href').split('/');
        out.push(ele[ele.length - 1])
    });
    console.log(out);
    if (out.length === 0) {
        window.alert('No entries selected for copying')
        return 0
    }
    var textArea = document.createElement("textarea");
    textArea.value = out;
    document.body.appendChild(textArea);
    textArea.focus();
    textArea.select();
    var successful = document.execCommand('copy');
    document.body.removeChild(textArea);
}


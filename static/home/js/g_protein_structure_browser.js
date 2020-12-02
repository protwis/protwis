function gproteinstructurebrowser() {
// $(document).ready(function () {
    // 'use strict';

    var oTable2 = $('#structures_scrollable').DataTable({
        "scrollY":        "65vh",
        "scrollX":        true,
        "scrollCollapse": true,
        "scroller": true,
        "paging":         false,
        // "bSortCellsTop": true,
        "aaSorting": [],
        "autoWidth": false,
        "order": [[29,'desc'],[1,'asc']],
        "columnDefs": [
            { "targets": 'no-sort', "orderable": false }
            ],
        "columns": [
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            {"width": "20%"},
            null,
            null,
            null,
            null,
            null,
            null,
            null //Not displayed, storing protein id
        ],
        "bInfo" : true,
    });

    var prev_ids = Array()
    var current_align_ids = Array()

    //Uncheck every row when using back button on browser
    $('.alt_selected').prop('checked',false)
    $('.alt').prop('checked',false)
    $('.select-all').prop('checked',false)
    //
    ClearSelection('targets');
    ClearSelection('reference');

    $("#loading_div").hide();

    yadcf.init(oTable2,
        [
            {
                column_number : 1,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Fam.",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '50px',
                }
            },
            {
                column_number : 2,
                filter_type: "multi_select",
                select_type: 'select2',
                column_data_type: "html",
                filter_default_label: "&alpha;",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '50px',
                }
            },
            {
                column_number: 3,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Species",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '55px',
                }
            },
            {
                column_number : 4,
                filter_type: "text",
                select_type: 'select2',
                filter_default_label: "Note",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                }
            },
            {
                column_number : 5,
                filter_type: "range_number",
                filter_default_label: ["From", "To"],
                select_type_options: {
                    width: '30px',
                }
            },
            {
                column_number : 6,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "&beta;",
                html_data_type: "text",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '40px',
                }
            },
            {
                column_number: 7,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Species",
                html_data_type: "text",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '55px',
                }
            },
            {
                column_number : 8,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "&gamma;",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '40px',
                }
            },
            {
                column_number: 9,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Species",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '55px',
                }
            },
            {
                column_number : 10,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Method",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '60px',
                }
            },
            {
                column_number : 11,
                filter_type: "multi_select",
                select_type: 'select2',
                column_data_type: "html",
                filter_default_label: "Select",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '50px',
                }
            },
            {
                column_number : 12,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Select",
                column_data_type: "html",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '50px',
                }
            },
            {
                column_number : 13,
                filter_type: "range_number",
                filter_default_label: ["From", "To"],
            },
            {
                column_number : 14,
                filter_type: "multi_select",
                select_type: 'select2',
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "UniProt",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '60px',
                }
            },
            {
                column_number : 15,
                filter_type: "multi_select",
                select_type: 'select2',
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "IUPHAR",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '60px',
                }
            },
            {
                column_number : 16,
                filter_type: "multi_select",
                select_type: 'select2',
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "Receptor family",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '120px',
                }
            },
            {
                column_number: 17,
                filter_type: "multi_select",
                select_type: 'select2',
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "Class",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                }
            },
            {
                column_number: 18,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Species",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '55px',
                }
            },
            {
                column_number : 19,
                filter_type: "text",
                select_type: 'select2',
                filter_default_label: "Receptor fusion",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '100px',
                }
            },
            {
                column_number : 20,
                filter_type: "text",
                select_type: 'select2',
                filter_default_label: "Antibodies",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '100px',
                }
            },
            {
                column_number : 21,
                filter_type: "text",
                select_type: 'select2',
                filter_default_label: "Other",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '100px',
                }
            },
            {
                column_number : 22,
                filter_type: "text",
                select_type: 'select2',
                html_data_type: "text",
                filter_default_label: "Ligand name",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '100px',
                }
            },
            {
                column_number : 23,
                filter_type: "multi_select",
                select_type: 'select2',
                html_data_type: "text",
                filter_default_label: "Ligand type",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '100px',
                },
                data: ['none', 'peptide', 'peptidesmall molecule', 'protein', 'small molecule', 'small moleculesmall molecule']
            },
            {
                column_number : 24,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Modality",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '100px',
                },
                data: ['AgonistPAM', 'Agonist', 'Apo (no ligand)']
            },
            {
                column_number : 25,
                filter_type: "text",
                select_type: 'select2',
                filter_default_label: "Ligand name",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '100px',
                }
            },
            {
                column_number : 26,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Ligand type",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '80px',
                }
            },
            {
                column_number : 27,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Last author",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '100px',
                }
            },
            {
                column_number : 28,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Reference",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '140px',
                }
            },
            {
                column_number : 29,
                filter_type: "range_date",
                date_format: "yyyy-mm-dd",
                select_type: 'select2',
                filter_default_label: ["From", "To"],
                // filter_reset_button_text: false,
            },
        ],
        {
            cumulative_filtering: false
        }
    );

    yadcf.exResetAllFilters(oTable2);

    // $(function(){
    //     $(".wrapper").scroll(function(){
    //         $(".dataTables_scrollBody").eq(0).scrollLeft($(".wrapper").scrollLeft());
    //     });
    //     $(".dataTables_scrollBody").eq(0).scroll(function(){
    //         $(".wrapper").scrollLeft($(".dataTables_scrollBody").eq(0).scrollLeft());
    //     });
    // });

    $('#structures_scrollable'+' > tbody > tr').click(function(event) {
        if (event.target.type !== 'checkbox') {
            $(':checkbox', this).trigger('click');
            $(this).eq(0).toggleClass('alt_selected');
            $(this).find('td').toggleClass('highlight');
        }
        $(this).eq(0).toggleClass('alt_selected');
        $(this).find('td').toggleClass('highlight');
    });

    $(".select-all").click(function() {
        $(':checkbox', this).trigger('click');
        if ($(this).prop('checked')===true) {
            $('.alt').prop('checked', true);
            $('.alt').parent().parent().addClass('alt_selected');
            $('.alt').parent().parent().find('td').addClass('highlight');
        }
        if ($(this).prop('checked')===false) {
            $('.alt').prop('checked', false);
            $('.alt').parent().parent().removeClass('alt_selected');
            $('.alt').parent().parent().find('td').removeClass('highlight');
        }
    });

    // $('.wrapper').find('div').width($(".yadcf-datatables-table--structures_scrollable").width());

    $('.hide_columns').click(function(evt) {
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

    toggle_enabled = true;
    $('#toggle_fixed_btn').click(function() {
        if (toggle_enabled) {
            toggle_enabled = false;
            $("#overlay").hide();
            $("#toggle_fixed_btn").attr('value',"Enable fixed columns");
            $("#toggle_fixed_btn").addClass('clicked_button');
        } else {
            toggle_enabled = true;
            $('.dataTables_scrollBody').scroll();
            $("#toggle_fixed_btn").attr('value',"Disable fixed columns");
            $("#toggle_fixed_btn").removeClass('clicked_button');
        }
    });

    $("#toggle_columns_btn").click(function() {
        var columns = Array.from(new Array(24), (x,i) => i + 6);
        columns.forEach(function(column) {
            var column = oTable2.column( column );
            try {
                column.visible( true, false );
            }
            catch(err) {
                column.visible( true, false );
            }
        });
        oTable2.draw();
    });

    $('#representative_btn').click(function () {
        $(this).toggleClass('toggled');
        if ($(this).hasClass('toggled')) {
            $("#representative_btn").addClass('clicked_button');
            $("#representative_btn").attr('value','All structures');
        }
        else {
            $("#representative_btn").removeClass('clicked_button');
            $("#representative_btn").attr('value','Representative structures (state & receptor)');
        }
        oTable2.draw();
    });

    $.fn.dataTable.ext.search.push(
        function (settings, data, dataIndex) {
            if ($('#representative_btn').hasClass('toggled')) {
                if ($(oTable2.row(dataIndex).node()).hasClass("repr-st") && $(oTable2.row(dataIndex).node()).hasClass("repr-st")) {
                    return true;
                }
                return false;
            }
            else {
                return true;
            }
    });

    $('#align_btn_g_prot').click(function () {
        var checked_data = oTable2.rows('.alt_selected').data();
        ClearSelection('targets');
        console.log(checked_data);
        for (i = 0; i < checked_data.length; i++) {
            var div = document.createElement("div");
            div.innerHTML = checked_data[i][30];
            if (typeof div.innerText !== "undefined") {
                AddToSelection('targets', 'protein', div.innerText.replace(/\s+/g, ''));
            } else {
                AddToSelection('targets', 'protein', div.textContent.replace(/\s+/g, ''));
            }
        }
        window.location.href = '/alignment/segmentselectiongprot';
    });

    $('#superpose_btn').click(function() {
        superposition(oTable2, [7,1,2,3,4,5,11,28], 'structure_browser');
    });

    $('#download_btn').click(function () {
        ClearSelection('targets');
        var checked_data = oTable2.rows('.alt_selected').data();
        for (i = 0; i < checked_data.length; i++) {
            var div = document.createElement("div");
            div.innerHTML = checked_data[i][7];
            if (typeof div.innerText !== "undefined") {
                AddToSelection('targets', 'structure',  div.innerText.replace(/\s+/g, '') );
            } else {
                AddToSelection('targets', 'structure', div.textContent.replace(/\s+/g, ''));
            }
        }
        window.location.href = '/structure/pdb_download_index';
    });

    // $('.glyphicon-export').mouseover(function() {
    //     window.alert($(this));
    // })
    $('.uniprot-export').data('powertipjq', $([
        '<p>Export UniProt IDs</p>'
        ].join('\n')));
    $('.pdb-export').data('powertipjq', $([
        '<p>Export PDB IDs</p>'
        ].join('\n')));
    $('.glyphicon-export').powerTip({
        placement: 'n',
        smartPlacement: true
    });
    $('#uniprot_copy').click(function () {
        copyToClipboard($('.alt_selected > .uniprot > a'), '\n', 'UniProt IDs', $('.uniprot-export'));
    });
    $('#pdb_copy').click(function () {
        copyToClipboard($('.alt_selected > .pdb > a'), '\n', 'PDB IDs', $('.pdb-export'));
    });

    $('#reset_filters_btn').click(function () {
        window.location.href = '/structure/'
    });

    $('.close_modal').click(function () {
        var modal = document.getElementById('myModal');
        modal.style.display = "none";
    });

    $('.dataTables_scrollBody').append('<div id=overlay><table id="overlay_table" class="row-border text-center compact dataTable no-footer text-nowrap"><tbody></tbody></table></div>');

    function create_overlay() {
        // This function fires upon filtering, to update what rows to show as an overlay
        $("#overlay_table tbody tr").remove();
        var $target = $("#overlay_table tbody");
        $("#structures_scrollable tbody tr").each(function() {
            var $tds = $(this).children(),
                $row = $("<tr></tr>");
            // $row.append($tds.eq(0).clone()).append($tds.eq(1).clone()).appendTo($target);
            $row.append($tds.eq(1).clone()).append($tds.eq(2).clone()).appendTo($target);
            $row.height($(this).height());
        });
        $("#overlay_table .border-right").removeClass("border-right");
    }

    // Function that detects filtering events
    $('#structures_scrollable').on( 'draw.dt', function (e,oSettings) {
        create_overlay();
    });

    create_overlay();
    $("#overlay").hide();

    var left = 0;
    var old_left = 0;
    $('.dataTables_scrollBody').scroll(function(){
        // If user scrolls and it's >100px from left, then attach fixed columns overlay
        left = $('.dataTables_scrollBody').scrollLeft();
        if (left!=old_left) $("#overlay").hide();
        old_left = left;

        if (left>100 && toggle_enabled) {
            $("#overlay").css({ left: left+'px' });
            if ($("#overlay").is(":hidden")) $("#overlay").show();
        }
    });
    // console.log($('#yadcf-filter--structures_scrollable-from-12').width());
    // $('#yadcf-filter--structures_scrollable-from-12').width(10);
    // console.log($('#yadcf-filter--structures_scrollable-from-12').width());

};

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

function copyDropdown() {
    document.getElementById("Dropdown").classList.toggle("show");
}

window.onclick = function(event) {
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

function copyToClipboard(array, delimiter, data_name, powertip_object=false) {
    var link = array;
    var out = '';
    link.each(function() {
        var ele = $(this).attr('href').split('/');
        out+=ele[ele.length-1]+delimiter;
    });
    if (out.length===0) {
        window.alert('No entries selected for copying')
        return 0
    }
    var textArea = document.createElement("textarea");
    textArea.value = out;
    document.body.appendChild(textArea);
    textArea.focus();
    textArea.select();
    try {
        var successful = document.execCommand('copy');
        var msg = successful ? 'Successful' : 'Unsuccessful';
        if (powertip_object!==false) {
            $.powerTip.hide();
            powertip_object.data('powertipjq', $([
                '<p>Copied to clipboard!</p>'
                ].join('\n')))
            powertip_object.powerTip('show');
            setTimeout(function() {
            powertip_object.data('powertipjq', $([
                '<p>Export '+data_name+'</p>'
                ].join('\n')))
            },1000);
        }
    } catch (err) {
        window.alert('Oops, unable to copy');
    }
    document.body.removeChild(textArea);
}
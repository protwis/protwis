let receptor_selection = [];
let oTable = [];


/**
 * Connects to an endpoint served by a django view
 */
function initializeTargetChooserTables() {
    $.get('/common/targettabledata', function (data) {
        $('#targetselect-modal-table .tableview').html(data);
        targettabledata = data;
    });

}

function AddMultipleTargets(receptor_selection) { // instead of a form an array which should be called receptor_selection
    // Deals with csrf token
    function getCookie(c_name) {
        if (document.cookie.length > 0) {
            c_start = document.cookie.indexOf(c_name + "=");
            if (c_start != -1) {
                c_start = c_start + c_name.length + 1;
                c_end = document.cookie.indexOf(";", c_start);
                if (c_end == -1) c_end = document.cookie.length;
                return unescape(document.cookie.substring(c_start, c_end));
            }
        }
        return "";
    }
    $.ajaxSetup({
        headers: { "X-CSRFToken": getCookie("csrftoken") }
    });
    rec_sel = Array.from(new Set(receptor_selection))
    //Actual post
    $.ajax({
        type: 'POST',
        url: '/common/targetformread',
        data: {
            "input-targets": rec_sel.join("\r")
        },
        cache: false,
        processData: false,
        contentType: false,
        'success': function (data) {
            $("#selection-targets").html(data);
        },
    }).fail(function (jqXHR, textStatus, error) {
        alert("Request failed: " + textStatus + error);
    });
    return false;
}

function processTableSelection(){

}

function thisTARGET(elem) {
    let mode = $('ul#mode_nav').find('li.active').find('a').text().trim();
   $('.pdb_selected:checked', oTable[mode].cells().nodes()).each(function() {
       receptor_selection.push($(this).attr('id'));
   });
    receptor_selection = Array.from(new Set(receptor_selection))
    console.log(receptor_selection)

}

function resetselection(not_update = false, reset_filters = false) {
    var mode = $('ul#mode_nav').find('li.active').find('a').text().trim();
    group = $('.tableview:visible').attr('group-number');
    if (group) mode = mode + group;

    $('.check_all:visible').prop('checked', false);

    $('input', oTable[mode].cells().nodes()).prop('checked', false);

//    if (!not_update) update_text_in_modal();

    if (reset_filters) yadcf.exResetAllFilters(oTable[mode]);
}

function check_all(elem, button) {
    var mode = $('ul#mode_nav').find('li.active').find('a').text().trim();
    show_all = $('.check_all:visible').prop("checked");
    if (button) {
        if (show_all) {
            $('.check_all:visible').prop("checked", false);
            show_all = false;
        } else {
            $('.check_all:visible').prop("checked", true);
            show_all = true;
        }
    }
}

$.fn.dataTable.ext.order['dom-checkbox'] = function(settings, col) {
    return this.api().column(col, {
        order: 'index'
    }).nodes().map(function(td, i) {
        return $('input', td).prop('checked') ? '1' : '0';
    });
};

function pastePDBs() {
    var mode = $('ul#mode_nav').find('li.active').find('a').text().trim();
    group = $('.tableview:visible').attr('group-number');
    if (group) mode = mode + group;
    pdbs = $('.pastePDBs:visible').val().toUpperCase().trim();
    pdbs = pdbs.split(/[ ,]+/);
    resetselection(1);
    $('.pdb_selected', oTable[mode].cells().nodes()).each(function() {
        pdb = $(this).attr('id')
        check = $.inArray(pdb, pdbs);
        if (check !== -1) {
            $(this).prop("checked", true);
            pdbs.splice(check, 1);
        }
    });
    if (pdbs.length) {
        var popOverSettings = {
            placement: 'bottom',
            container: 'body',
            title: 'One or more PDB(s) not found:', //Sepcify the selector here
            content: pdbs
        }
        $('.pastePDBs').popover(popOverSettings);
        $('.pastePDBs').popover('show');
        $('.pastePDBs').on('shown.bs.popover',function() {
            setTimeout(function() {
                $('.pastePDBs').popover("destroy")}, 3000);
        })
    } /*else {
        $('.pastePDBs').popover('destroy');
    }*/
    oTable[mode].order([
        [1, 'desc']
    ]);
    oTable[mode].draw();
}

function exportPDBs() {
    var mode = $('ul#mode_nav').find('li.active').find('a').text().trim();
    group = $('.tableview:visible').attr('group-number');
    if (group) mode = mode + group;
    var pdbs = [];
    $('.pdb_selected:checked', oTable[mode].cells().nodes()).each(function() {
        pdbs.push($(this).attr('id'));
    });


    var textArea = document.createElement("textarea");
    textArea.value = pdbs;
    document.body.appendChild(textArea);
    textArea.focus();
    textArea.select();
    try {
        var successful = document.execCommand('copy');
        var msg = successful ? 'Successful' : 'Unsuccessful';
        $(".export_pdbs").html(msg);
        setTimeout("$('.export_pdbs').html('Export selected PDB codes');", 4000);
    } catch (err) {
        $(".export_pdbs").html('Oops, unable to copy');
    }

    document.body.removeChild(textArea);
}

function showTARGETtable(element) {
    var mode = $('ul#mode_nav').find('li.active').find('a').text().trim();
    group = $(element + ' .tableview').attr('group-number');
    if (group) mode = mode + group;
    console.log(element, mode, group);

    mode_without_space = mode.replace(/ /g, '_');

    if (!$.fn.DataTable.isDataTable(element + ' .tableview table')) {
        console.log(mode);

        $(element + ' .modal-header .pastePDBs').keypress(function(event) {
            var keycode = (event.keyCode ? event.keyCode : event.which);
            if(keycode == '13'){
                pastePDBs();
            }
        });

        // Exports selected id's
        $(element + ' .modal-header').append(' | <span><button type="button" onclick="exportPDBs();" ' +
            'class="btn btn-xs btn-primary export_pdbs">Export</button></span>');

        console.time('DataTable');
        oTable[mode] = $(element + ' .tableview table').DataTable({
            dom: "ftip",
            deferRender: true,
            scrollY: '50vh',
            scrollX: true,
            scrollCollapse: true,
            scroller: true,
            paging: false,
//        lengthMenu: [[10, 25, 50, -1], [10, 25, 50, "All"]],
            bSortCellsTop: false, //prevent sort arrows going on bottom row
            aaSorting: [],
            autoWidth: true,
//            pageLength: -1,
            bInfo: true,
            columnDefs: [{
                targets: 0,
                orderable: false,
                className: 'select-checkbox'
            },],
        });
        console.timeEnd('DataTable');
        console.time('yadcf');
        yadcf.init(oTable[mode],
            [
               {
                   column_number : 0,
                   column_data_type: "html",
                   html_data_type: "text",
               },
                {
                    column_number: 1,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Class",
//                    filter_reset_button_text: false,
                },
                {
                    column_number: 2,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Family",
//                    filter_reset_button_text: true,
                },
                {
                    column_number: 3,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Uniprot",
//                    filter_reset_button_text: true,
                },
                {
                    column_number: 4,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "IUPHAR",
                    html_data_type: "text",
//                    filter_reset_button_text: true,
                },

            ], {
//                cumulative_filtering: true,
                filters_tr_index: 1
            }
        );
        console.timeEnd('yadcf');
//        yadcf.exResetAllFilters(oTable);

    };

    oTable[mode].draw();
}

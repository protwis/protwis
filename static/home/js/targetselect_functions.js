let receptor_selection = [];
let oTable = [];


/**
 * Connects to an endpoint served by a django view
 */
function initializeTargetChooserTables() {
    $.get('/common/targettabledata', function (data) {
        $('#targetselect-modal-table .tableview').html(data);
        let targettabledata = data;
    });
}

/**
 * POST array of targets (receptor_selection) to receiver endpoint /common/targetformread
 */
function AddMultipleTargets() {
    $.post('/common/targetformread',
        { "input-targets": receptor_selection.join("\r") },
        function (data) {
            $("#selection-targets").html(data);
        }).fail(function (jqXHR, textStatus, error) {
        alert("Request failed: " + textStatus + error);
    });
}

/**
 * Onclick function triggered in input element in html table
 * @param name
 * @returns receptor_selection
 */
function thisTARGET(name) {
    const checkboxes = document.querySelectorAll(`input[name="targets"]:checked`);
    receptor_selection =[];
    checkboxes.forEach((checkbox) => {
        receptor_selection.push(checkbox.id);
    });
    return receptor_selection;
}

function resetselection(not_update = false, reset_filters = false) {
    $('.check_all:visible').prop('checked', false);
    $('input', oTable.cells().nodes()).prop('checked', false);
    if (reset_filters) yadcf.exResetAllFilters(oTable);
}

function check_all(elem, button) {
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
    pdbs = $('.pastePDBs:visible').val().toUpperCase().trim();
    pdbs = pdbs.split(/[ ,]+/);
    resetselection(1);
    $('.pdb_selected', oTable.cells().nodes()).each(function() {
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
            title: 'One or more PDB(s) not found:', //Specify the selector here
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
    oTable.order([
        [1, 'desc']
    ]);
    oTable.draw();
}

function exportPDBs() {
    var pdbs = [];
    $('.pdb_selected:checked', oTable.cells().nodes()).each(function() {
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

    if (!$.fn.DataTable.isDataTable(element + ' .tableview table')) {

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
        oTable = $(element + ' .tableview table').DataTable({
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
        yadcf.init(oTable,
            [
                {
                    column_number: 1,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Class",
                    filter_reset_button_text: true,
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
                    column_data_type: "html",
                    html_data_type: "text",
                    filter_default_label: "IUPHAR",
                    filter_match_mode : "exact",
                },

            ], {
//                cumulative_filtering: true,
                filters_tr_index: 1
            }
        );
        console.timeEnd('yadcf');
    };

    oTable.draw();

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



}

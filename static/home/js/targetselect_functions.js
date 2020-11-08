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

/**
 * Similar to thisTARGET but selects all checkboxes in one go, that is, check_all
 */
function check_all() {
    show_all = $('.check_all:visible').prop("checked");

    if (show_all) {
        $('.pdb_selected:visible').prop("checked", true);
    } else {
        $('.pdb_selected:visible').prop("checked", false);
    }

    const checkboxes = document.querySelectorAll(`input[name="targets"]:checked`);
    receptor_selection =[];
    checkboxes.forEach((checkbox) => {
        receptor_selection.push(checkbox.id);
    });

}

function resetselection(not_update = false, reset_filters = false) {
    $('.check_all:visible').prop('checked', false);
    $('input', oTable.cells().nodes()).prop('checked', false);
    if (reset_filters) yadcf.exResetAllFilters(oTable);
}



$.fn.dataTable.ext.order['dom-checkbox'] = function(settings, col) {
    return this.api().column(col, {
        order: 'index'
    }).nodes().map(function(td, i) {
        return $('input', td).prop('checked') ? '1' : '0';
    });
};

function pasteTargets() {
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

function exportSelectedTargets() {
    var targets = [];
    // TODO Replace with uniprot IDs
    /*$('.pdb_selected:checked', oTable.cells().nodes()).each(function() {
        pdbs.push($(this).attr('id'));
    });*/

    // TODO add export button
    var textArea = document.createElement("textarea");
    textArea.value = pdbs;
    document.body.appendChild(textArea);
    textArea.focus();
    textArea.select();
    try {
        var successful = document.execCommand('copy');
        var msg = successful ? 'Successful' : 'Unsuccessful';
        $(".export_targets").html(msg);
        setTimeout("$('.export_targets').html('Export selected targets');", 4000);
    } catch (err) {
        $(".export_targets").html('Oops, unable to copy');
    }

    // remove area for copying
    document.body.removeChild(textArea);
}

var targetTable;
function initTargetTable(element) {

    if (!$.fn.DataTable.isDataTable(element + ' table')) {
        targetTable = $(element + ' table').DataTable({
            dom: "ftip",
            deferRender: true,
            scrollY: '50vh',
            scrollX: true,
            scrollCollapse: true,
            scroller: true,
            paging: false,
            bSortCellsTop: false, //prevent sort arrows going on bottom row
            aaSorting: [],
            autoWidth: true,
            bInfo: true,
//            order: [[ 1, "asc" ], [ 2, "asc" ], [ 3, "asc" ], [ 4, "asc" ]],
            columnDefs: [{
                targets: 0,
                orderable: false,
                className: 'select-checkbox'
            },{
                targets: 1,
                orderable: false,
                className: 'text-center'
            },],
        });

        yadcf.init(targetTable,
            [
                {
                    column_number: 1,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Class",
                    filter_reset_button_text: false,
                    style_class: "center",
                },
                {
                    column_number: 2,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Ligand",
                    filter_reset_button_text: false,
                    filter_match_mode : "exact",
                },
                {
                    column_number: 3,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Family",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 4,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Uniprot",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 5,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    column_data_type: "html",
                    html_data_type: "text",
                    filter_default_label: "GtP",
                    filter_match_mode : "exact",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 6,
                    filter_type: "text",
                    select_type: 'select2',
                    html5_data: "data-search",
                    filter_default_label: "PDB",
                    filter_reset_button_text: false,
                },
                /*{
                    column_number: 7,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Approved",
                },
                {
                    column_number: 8,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Clinical trial",
                },*/
                {
                    column_number: 7,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Gs",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 8,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Gi/o",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 9,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Gq/11",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 10,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "G12/13",
                    filter_reset_button_text: false,
                },

            ], {
//                cumulative_filtering: true,
                filters_tr_index: 1
            }
        );
    };

    targetTable.draw();

    // Put top scroller
    // https://stackoverflow.com/questions/35147038/how-to-place-the-datatables-horizontal-scrollbar-on-top-of-the-table
    $('.dataTables_scrollHead').css({
        'overflow-x':'scroll'
    }).on('scroll', function(e){
        var scrollBody = $(this).parent().find('.dataTables_scrollBody').get(0);
        scrollBody.scrollLeft = this.scrollLeft;
        $(scrollBody).trigger('scroll');
    });

    // Add selection counter
    updateTargetCount();
}


// Add to buttons
//  yadcf.exResetAllFilters(targetTable);

/*function clearTargetSelection(){
  $('table#uniprot_selection.tbody.checkbox').prop('checked',false)
}*/


var changedTargetBoxes = 0;
function check_all_targets(){
  changedTargetBoxes = 0;
  $("table#uniprot_selection tbody tr").each(function() {
    if (!$(this).hasClass("selected")){
      $(this).addClass("selected");
      var checkbox = $(this).find("[type=checkbox]");
      checkbox.prop("checked", $(this).hasClass("selected"));
      changedTargetBoxes++;
    }
  });

  if (changedTargetBoxes==0){
    $("table#uniprot_selection tbody tr").each(function() {
      $(this).removeClass("selected");
      var checkbox = $(this).find("[type=checkbox]");
      checkbox.prop("checked", false);
    });
  }

  updateTargetCount();
  return false;
}

// TODO: maintain using array or count using the DOM?
function updateTargetCount(){
  // Counting the selected targets matching the current filters
  var numTargets = $('table#uniprot_selection tbody input:checked').length;
  var message = " ("+numTargets+" targets selected )";
  if (numTargets == 1)
    message = " ("+numTargets+" target selected )";

  // Tweak selection message + integrate filtering and selection
  if ($("#uniprot_selection_info_targets").length == 0)
    $("#uniprot_selection_info").append("<span id=\"uniprot_selection_info_targets\"></span>");

  $("#uniprot_selection_info_targets").html(message);

  //var allTargets = $('table#uniprot_selection tbody input').length;
}

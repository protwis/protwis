function update_text_in_modal() {

    var pdbs = [];
    var mode = $('ul#mode_nav').find('li.active').find('a').text().trim();
    group = $('.tableview:visible').attr('group-number');
    if (group) mode = mode + group;
    $('.pdb_selected', oTable[mode].cells().nodes()).each(function() {
        if ($(this).prop("checked")) {
            $(this).closest("tr").addClass("selected");
            $(this).closest(".dataTables_scrollBody").find("#overlay_" + $(this).attr('id')).addClass("selected");
            pdbs.push($(this).attr('id'));
        } else {
            $(this).closest("tr").removeClass("selected");
            $(this).closest(".dataTables_scrollBody").find("#overlay_" + $(this).attr('id')).removeClass("selected");
        }
    });
    var mode = $('ul#mode_nav').find('li.active').find('a').text().trim();

    if (mode == 'Single group of structures' || $("#single-group-tree-tab").length) {
        var total_selected = pdbs.length
        var selected_visible = $('.dataTables_scrollBody:visible .pdb_selected:checked').length
        var ModalpdbsCountSelector = '#single-crystal-group-pdbs-modal-text';

        if (total_selected == selected_visible) {
            $(ModalpdbsCountSelector).html(total_selected + ' structure(s) selected');
        } else {
            $(ModalpdbsCountSelector).html(total_selected + ' structure(s) selected (' + (total_selected - selected_visible) + ' currently filtered)');
        }

        var pdbsInputSelector = '#single-crystal-group-tab .crystal-pdb';
        var pdbsCountSelector = '#single-crystal-group-tab .crystal-count';
        $(pdbsInputSelector).val(JSON.stringify(pdbs));
        $(pdbsCountSelector).html(pdbs.length);

    } else if (mode == 'Two groups of structures') {
        var total_selected = pdbs.length;
        var selected_visible = $('.pdb_selected:checked:visible').length;
        var ModalpdbsCountSelector = '#two-crystal-group-pdbs-modal-' + group + '-text';

        if (total_selected == selected_visible) {
            $(ModalpdbsCountSelector).html(total_selected + ' structure(s) selected');
        } else {
            $(ModalpdbsCountSelector).html(total_selected + ' structure(s) selected (' + (total_selected - selected_visible) + ' currently filtered)');
        }

        var pdbsInputSelector = '#two-crystal-groups-tab .crystal-group-' + group + '-pdbs';
        var pdbsCountSelector = '#two-crystal-groups-tab .crystal-count-' + group;
        $(pdbsInputSelector).val(JSON.stringify(pdbs));
        $(pdbsCountSelector).html(pdbs.length);

    }
}

function thisPDB(elem) {
    var mode = $('ul#mode_nav').find('li.active').find('a').text().trim();
    var ReceptorName = $(elem).attr('long');
    var pdbName = $(elem).attr('id');
    // console.log('thisPDB',pdbName);
    if (mode == 'Single structure') {
        $('input', oTable[mode].cells().nodes()).not(elem).prop('checked', false);
        var pdbs = [];
        if ($(elem).prop("checked")) {
            pdbs.push(pdbName);
            // Update view
            $(".crystal-count:visible").html(ReceptorName + ' - ' + pdbName + ' selected.');
        } else {
            // Update view
            $(".crystal-count:visible").html('No structure selected.');
        }
        $(".crystal-count:visible").parent().parent().find('.crystal-pdb').val(JSON.stringify(pdbs));
    }
    update_text_in_modal();
}

function resetselection(not_update) {
    var mode = $('ul#mode_nav').find('li.active').find('a').text().trim();
    group = $('.tableview:visible').attr('group-number');
    if (group) mode = mode + group;

    $('.check_all:visible').prop('checked', false);

    $('input', oTable[mode].cells().nodes()).prop('checked', false);

    if (!not_update) update_text_in_modal();
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

    if (mode == 'Single group of structures' || $("#single-group-tree-tab").length) {
        var pdbs = [];

        // REMOVE EXISITING? Probably not, more logical that filtering adds more
        // $('input', oTable.cells().nodes()).prop('checked',false);

        if (show_all) {
            $('.pdb_selected:visible').prop("checked", true);
        } else {
            $('.pdb_selected:visible').prop("checked", false);
        }
    } else if (mode == 'Two groups of structures') {
        group = $(elem).closest('.tableview').attr('group-number');
        if (group) mode = mode + group;
        var pdbs = [];

        if (show_all) {
            $('.pdb_selected:visible').prop("checked", true);
        } else {
            $('.pdb_selected:visible').prop("checked", false);
        }
    }
    if (button) {
        update_text_in_modal();
    }
}

function check_all_representatives() {
    var mode = $('ul#mode_nav').find('li.active').find('a').text().trim();
    group = $('.tableview:visible').attr('group-number');
    if (group) mode = mode + group;
    $('input', oTable[mode].cells().nodes()).prop('checked', false);
    $('input[representative="Yes"]:visible').each(function() {
        $(this).prop("checked", true);
    });
    update_text_in_modal();
}

function check_all_distance_representatives() {
    var mode = $('ul#mode_nav').find('li.active').find('a').text().trim();
    group = $('.tableview:visible').attr('group-number');
    if (group) mode = mode + group;
    $('input', oTable[mode].cells().nodes()).prop('checked', false);
    $('input[distance_representative="Yes"]:visible').each(function() {
        // pdbs.push($(this).attr('id'));
        $(this).prop("checked", true);
    });
    update_text_in_modal();
}

function check_all_class_representatives() {
    var mode = $('ul#mode_nav').find('li.active').find('a').text().trim();
    group = $('.tableview:visible').attr('group-number');
    if (group) mode = mode + group;
    $('input', oTable[mode].cells().nodes()).prop('checked', false);
    $('input[class_consensus_based_representative="Yes"]:visible').each(function() {
        // pdbs.push($(this).attr('id'));
        $(this).prop("checked", true);
    });
    update_text_in_modal();
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
    pdbs = $('.pastePDBs:visible').val().toUpperCase();
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
            title: 'PDBs not found', //Sepcify the selector here
            content: pdbs
        }
        $('.pastePDBs').popover(popOverSettings);
        $('.pastePDBs').popover('show');
    } else {
        $('.pastePDBs').popover('destroy');
    }
    oTable[mode].order([
        [20, 'desc']
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

var oTable = [];

function showPDBtable(element) {
    var mode = $('ul#mode_nav').find('li.active').find('a').text().trim();
    group = $(element + ' .tableview').attr('group-number');
    if (group) mode = mode + group;
    // console.log(element,mode,group);

    if (!$.fn.DataTable.isDataTable(element + ' .tableview table')) {
        console.log(mode);

        $(element + ' .tableview').before('<span><button type="button" onclick="check_all(this,1);" class="btn btn-xs btn-primary reset-selection">Select all displayed</button></span>');
        $(element + ' .tableview').before(' | <span><input type=text class="pastePDBs" placeholder="Paste pdbs with comma- or space-separated"><button type="button" onclick="pastePDBs();" class="btn btn-xs btn-primary reset-selection">Load PDB codes</button></span>');
        $(element + ' .tableview').before(' | <span><button type="button" onclick="exportPDBs();" class="btn btn-xs btn-primary export_pdbs">Export selected PDB codes</button></span>');
        if (window.location.href.endsWith("contactnetwork/clustering") || window.location.href.endsWith("contactnetwork/clustering#"))
          $(element + ' .tableview').before(' | <span>Structure shortest distance to all other structures of the same receptor and same state: <button type="button" onclick="check_all_distance_representatives();" class="btn btn-xs btn-primary">Distance Representative</button></span>');
        else {
          $(element + ' .tableview').before(' | <span>Structure with highest % identity to GPCR’s contact consensus: <button type="button" onclick="check_all_representatives();" class="btn btn-xs btn-primary">Contact Representative</button></span>');
          $(element + ' .tableview').before(' | <span>Structure sharing either highest/lowest diff between fraction of active/inactive class consensus contacts, or for intermediate the one closes to a 0 diff: <button type="button" onclick="check_all_class_representatives();" class="btn btn-xs btn-primary">New Representative</button></span>');
        }

        oTable[mode] = $(element + ' .tableview table').DataTable({
            'scrollX': true,
            // 'paging': true,
            // 'autoWidth': true,

            scrollY: '75vh',
            // scrollCollapse: true,
            paging: false,
            columnDefs: [{
                targets: 'no-sort',
                orderable: false
            }],
            "aaSorting": [],
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
                {
                    "orderDataType": "dom-checkbox"
                }
            ]
        });
        console.log('done datatable');
        yadcf.init(oTable[mode],
            [{
                    column_number: 0,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "UniProt",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 1,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    column_data_type: "html",
                    html_data_type: "text",
                    filter_default_label: "Receptor",
                    filter_match_mode: "exact",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 2,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    html_data_type: "text",
                    select_type_options: {
                        width: '150px'
                    },
                    filter_default_label: "Rec Family",
                    filter_match_mode: "exact",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 3,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Class",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 4,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Species",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 5,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Method",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 6,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    select_type_options: {
                        width: '70px'
                    },
                    filter_default_label: "PDB",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 7,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    select_type_options: {
                        width: '70px'
                    },
                    filter_default_label: "Res (Å)",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 8,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "State",
                    filter_match_mode: "exact",
                    filter_reset_button_text: false,

                },
                // {
                //     column_number : 9,
                //     filter_type: "multi_select",
                //     select_type: 'select2',
                //     filter_default_label: "",
                //     filter_match_mode : "exact",
                //     filter_reset_button_text: false,

                // },
                // {
                //     column_number : 10,
                //     filter_type: "multi_select",
                //     select_type: 'select2',
                //     filter_default_label: "",
                //     filter_match_mode : "exact",
                //     filter_reset_button_text: false,

                // },
                {
                    column_number: 12,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Contact rep.",
                    filter_reset_button_text: false,

                },
                {
                    column_number: 13,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "7TM Open IC (Å)",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 14,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "G protein",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 15,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "B arrestin",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 16,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Fusion",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 17,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Antibody",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 18,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Ligand",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 19,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Ligand function",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 20,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Ligand type",
                    filter_reset_button_text: false,
                },
            ], {
                cumulative_filtering: false
            }
        );

        yadcf.exResetAllFilters(oTable[mode]);
        console.log('done yadcf');

        oTable[mode].on('draw.dt', function(e, oSettings) {
            create_overlay(element + ' .structure_selection');
            update_text_in_modal();
            $('.structure_overlay tr').css('cursor', 'pointer');
        });


        $(element + ' .dataTables_scrollBody').append('<div class="structure_overlay"><table id="overlay_table" class="overlay_table row-border text-center compact dataTable no-footer text-nowrap"><tbody></tbody></table></div>');
        $(element + " .structure_overlay").hide();

        $(element + ' .dataTables_scrollBody').before("<div class='top_scroll'><div>&nbsp;</div></div>");

        $('.top_scroll').css({
            'width': '100%',
            'overflow-x': 'scroll',
            'overflow-y': 'auto',
        });
        $('.top_scroll div').css({
            // 'background-color': 'red',
            'font-size': '1px',
            'line-height': '1px',
        });

        $('.structure_overlay').css({
            'top': '0px',
            'position': 'absolute',
            'background': '#f8f8f8',
            '-webkit-box-shadow': '5px 0 2px -2px #888',
            'box-shadow': '5px 0 2px -2px #888',
        });

        $('.structure_overlay tbody tr').css({
            'background-color': '#f8f8f8',
        });

        create_overlay(element + ' .structure_selection');
        track_scrolling(element);
        console.log('done overlays');


        $(function() {
            var tableContainer = $(".dataTables_scrollBody");
            var table = $(".dataTables_scrollBody table");
            var fakeContainer = $(".top_scroll");
            var fakeDiv = $(".top_scroll div");

            var tableWidth = table.width();
            fakeDiv.width(tableWidth);

            fakeContainer.scroll(function() {
                tableContainer.scrollLeft(fakeContainer.scrollLeft());
            });
            tableContainer.scroll(function() {
                fakeContainer.scrollLeft(tableContainer.scrollLeft());
            });
        })
        console.log('done scrolling');

        $('.dataTables_scrollBody tr').css('cursor', 'pointer');
        $('.dataTables_scrollBody tr').click(function(event) {
            if (event.target.type !== 'checkbox') {
                $(':checkbox', this).trigger('click');
                if ($(':checkbox', this).length === 0) {
                    pdb_id = $(this).attr('id').split("_")[1];
                    var checkbox = $(this).closest(".dataTables_scrollBody").find('#' + pdb_id);
                    checkbox.trigger('click');
                }
            }
        });
        console.log('done click event');

    };
    $(element+" .loading_overlay").hide();
}

function create_overlay(table_id) {
    // This function fires upon filtering, to update what rows to show as an overlay
    $(".overlay_table tbody tr").remove();
    var $target = $(".overlay_table tbody");
    $(table_id + " tbody tr").each(function() {
        var $tds = $(this).children(),
            $row = $("<tr id='overlay_" + $tds.eq(6).html() + "'></tr>");
        // $row.append($tds.eq(0).clone()).append($tds.eq(1).clone()).appendTo($target);
        $row.append($tds.eq(0).clone()).append($tds.eq(1).clone()).append($tds.eq(2).clone()).appendTo($target);
    });
    $(".overlay_table .border-right").removeClass("border-right");

    // rebind click event
    $('.structure_overlay2 tr').click(function(event) {
        console.log("clicking overlay tr");
        if (event.target.type !== 'checkbox') {
            $(':checkbox', this).trigger('click');
            pdb_id = $(this).attr('id').split("_")[1];
            console.log(pdb_id);
            if ($(':checkbox', this).length === 0) {
                pdb_id = $(this).attr('id').split("_")[1];
                console.log("2", pdb_id, '#' + pdb_id + ':checkbox:visible', $('#' + pdb_id + ':checkbox:visible'));
                // $('#'+pdb_id+':checkbox:visible').find("tr").trigger('click');
            }
        }
    });
}

function track_scrolling(element) {
    var left = 0;
    var old_left = 0;
    toggle_enabled = true;
    $(element + ' .dataTables_scrollBody').scroll(function() {
        // If user scrolls and it's >100px from left, then attach fixed columns overlay
        left = $(element + ' .dataTables_scrollBody').scrollLeft();
        if (left != old_left) $(".structure_overlay").hide();
        old_left = left;

        if (left > 100 && toggle_enabled) {
            $(".structure_overlay").css({
                left: left + 'px'
            });
            if ($(".structure_overlay").is(":hidden")) $(".structure_overlay").show();
        }
    });
}

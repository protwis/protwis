function update_text_in_modal() {

    var pdbs = [];
    var mode = $('ul#mode_nav').find('li.active').find('a').text().trim();

    group = $('.tableview:visible').attr('group-number');
    if (group) mode = mode + group;
    $('.pdb_selected', oTable[mode].cells().nodes()).each(function() {
        if ($(this).prop("checked")) {
            $(this).closest("tr").addClass("selected");
            $(this).closest(".dataTables_scroll").find("#overlay_" + $(this).attr('id')).addClass("selected");
            selector = $(this).closest(".dataTables_scrollBody").find("#overlay_" + $(this).attr('id')).find("checkbox");
            selector.prop("checked", !selector.prop("checked"));
            pdbs.push($(this).attr('id'));
        } else {
            $(this).closest("tr").removeClass("selected");
            $(this).closest(".dataTables_scroll").find("#overlay_" + $(this).attr('id')).removeClass("selected");
        }
    });
    var mode = $('ul#mode_nav').find('li.active').find('a').text().trim();

    if (mode == 'Single set of structures' || $("#single-group-tree-tab").length) {
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

    } else if (mode == 'Two sets of structures') {
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

    var referenceObject = $(elem).closest(".dataTables_scroll")

    // Fixed/frozen columns check
    if (elem.id.startsWith("overlaycheck_")){
        elem = $(".dataTables_scrollBody:visible").find("#" + elem.id.replace("overlaycheck_", ""))[0]

       $(elem).prop("checked", !elem.checked)
    } else {
      other = referenceObject.find("#overlaycheck_" + elem.id)
      other.prop("checked", elem.checked)
    }

    var pdbName = $(elem).attr('id');
    if (mode == 'Single structure') {
        if ($('input', oTable[mode].cells().nodes()).filter(":checkbox").not(elem).length > 0) {
          $('input', oTable[mode].cells().nodes()).filter(":checkbox").not(elem).each(function(i,e) { if (referenceObject.find("#overlaycheck_" + e.id).length > 0) referenceObject.find("#overlaycheck_" + e.id)[ 0 ].checked = false}) // deselect from overlay
          $('input', oTable[mode].cells().nodes()).filter(":checkbox").not(elem).prop('checked', false); // deselect from original table
        }
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

function resetselection(not_update = false, reset_filters = false) {
    var mode = $('ul#mode_nav').find('li.active').find('a').text().trim();
    group = $('.tableview:visible').attr('group-number');
    if (group) mode = mode + group;

    $('.check_all:visible').prop('checked', false);

    $('input', oTable[mode].cells().nodes()).prop('checked', false);

    if (!not_update) update_text_in_modal();

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

    if (mode == 'Single set of structures' || $("#single-group-tree-tab").length) {
        var pdbs = [];

        // REMOVE EXISITING? Probably not, more logical that filtering adds more
        // $('input', oTable.cells().nodes()).prop('checked',false);

        if (show_all) {
            $('.pdb_selected:visible').prop("checked", true);
        } else {
            $('.pdb_selected:visible').prop("checked", false);
        }
    } else if (mode == 'Two sets of structures') {
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
        $(element + ' .tableview').before(' | <span><input type=text class="pastePDBs" placeholder="Paste pdbs with comma- or space-separated"><button type="button" onclick="pastePDBs();" class="btn btn-xs btn-primary reset-selection">Load PDB codes</button></span>');
        $(element + ' .tableview').before(' | <span><button type="button" onclick="exportPDBs();" class="btn btn-xs btn-primary export_pdbs">Export selected PDB codes</button></span>');
        if (window.location.href.endsWith("contactnetwork/clustering") || window.location.href.endsWith("contactnetwork/clustering#"))
          $(element + ' .tableview').before(' | <span>Structure shortest distance to all other structures of the same receptor and same state: <button type="button" onclick="check_all_distance_representatives();" class="btn btn-xs btn-primary">Distance Representative</button></span>');
        else {
          // a$(element + ' .tableview').before(' | <span>Structure with highest % identity to GPCR’s contact consensus: <button type="button" onclick="check_all_representatives();" class="btn btn-xs btn-primary">Contact Representative</button></span>');
          // $(element + ' .tableview').before(' | <span>Structure sharing either highest/lowest diff between fraction of active/inactive class consensus contacts, or for intermediate the one closes to a 0 diff: <button type="button" onclick="check_all_class_representatives();" class="btn btn-xs btn-primary">New Representative</button></span>');
        }
        $(element + ' .tableview').before(' | <div class="externalfilters" style="display: inline-block;"><span id=external_filter_container_0></span></div>');
        $(element + ' .tableview').before('<div class="externalfilters" style="display: inline-block;"><span id=external_filter_container_1></span></div>');

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
            },
                {"targets": [ -1, -2 ],
                "visible": false}],
            "aaSorting": [],
            "columns": [
                {
                  "orderDataType": "dom-checkbox"
                },
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
                null
            ]
        });
        console.log('done datatable');
        yadcf.init(oTable[mode],
            [{
                    column_number: 1,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "UniProt",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 2,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    column_data_type: "html",
                    html_data_type: "text",
                    filter_default_label: "Receptor",
                    filter_match_mode: "exact",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 3,
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
                    column_number: 4,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Class",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 5,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Species",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 6,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Method",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 7,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    select_type_options: {
                        width: '70px'
                    },
                    filter_default_label: "PDB",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 8,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    select_type_options: {
                        width: '70px'
                    },
                    filter_default_label: "Res (Å)",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 9,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "State",
                    filter_match_mode: "exact",
                    filter_reset_button_text: false,

                },
                // {
                //     column_number : 10,
                //     filter_type: "multi_select",
                //     select_type: 'select2',
                //     filter_default_label: "",
                //     filter_match_mode : "exact",
                //     filter_reset_button_text: false,

                // },
                // {
                //     column_number : 11,
                //     filter_type: "multi_select",
                //     select_type: 'select2',
                //     filter_default_label: "",
                //     filter_match_mode : "exact",
                //     filter_reset_button_text: false,

                // },
                {
                    column_number: 13,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Contact rep.",
                    filter_reset_button_text: false,

                },
                /*{
                    column_number: 13,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "7TM Open IC (Å)",
                    filter_reset_button_text: false,
                },*/
                {
                    column_number: 15,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "G prot",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 16,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "B arr",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 17,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Fusion",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 18,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Antibody",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 19,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Ligand",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 20,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Ligand function",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 21,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Ligand type",
                    filter_reset_button_text: false,
                },
                {
                    column_number: 22, 
                    filter_container_id: 'external_filter_container_0',
                    html_data_type: "text", 
                    select_type: 'select2',
                    // filter_type: "multi_select",
                    filter_default_label: "All species",
                    filter_reset_button_text: false,
                    select_type_options: {
                        width: '250px',
                        minimumResultsForSearch: -1 // remove search box
                    },
                },
                {
                    column_number: 23, 
                    filter_container_id: 'external_filter_container_1',
                    html_data_type: "text", 
                    select_type: 'select2',
                    // filter_type: "multi_select",
                    filter_default_label: "All Structures",
                    filter_reset_button_text: false,
                    select_type_options: {
                        width: '250px',
                        minimumResultsForSearch: -1 // remove search box
                    },
                },
            ], {
                cumulative_filtering: false
            }
        );
        yadcf.exResetAllFilters(oTable[mode]);
        yadcf.exFilterColumn(oTable[mode], [
            [22, "Only show mammalian receptor structures (even if the non-mammalian is the only)"],
            [23, "Only show structures from human or the closest species (for each receptor and state)"]
          ]); 
        console.log('done yadcf');

        oTable[mode].on('draw.dt', function(e, oSettings) {
            create_overlay(element + ' .structure_selection');
            update_text_in_modal();
        });


        // $(element + ' .dataTables_scrollBody').append('<div class="structure_overlay"><table id="overlay_table" class="overlay_table row-border text-center compact dataTable no-footer text-nowrap"><tbody></tbody></table></div>');
        $(element + ' .dataTables_scroll').append('<div class="structure_overlay"><table id="overlay_table" class="overlay_table row-border text-center compact dataTable no-footer text-nowrap"><tbody></tbody></table></div>');
        $(element + " .structure_overlay").hide();

        $(element + ' .dataTables_scrollBody').before("<div class='top_scroll'><div>&nbsp;</div></div>");

        var dataTables_scrollBody_height = $(element + ' .dataTables_scrollBody')[0].offsetHeight;
        var bodyRect = $(element + ' .dataTables_scroll')[0].getBoundingClientRect(),
        elemRect = $(element + ' .dataTables_scrollBody')[0].getBoundingClientRect(),
        offset = elemRect.top - bodyRect.top;
        
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
            // 'top': '0px',
            'top': offset+'px',
            'position': 'absolute',
            'background': '#f8f8f8',
            '-webkit-box-shadow': '5px 0 2px -2px #888',
            'box-shadow': '5px 0 2px -2px #888',
            'height': (dataTables_scrollBody_height-17) + 'px',
            'overflow-y': 'scroll',
            'scrollbar-width': 'none',
        });
        

        $('.structure_overlay tbody tr').css({
            'background-color': '#f8f8f8',
        });

        create_overlay(element + ' .structure_selection');
        // track_scrolling(element);
        new_track_scrolling(element);
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
            $row = $("<tr id='overlay_" + $tds.eq(7).html() + "'></tr>");
        // $row.append($tds.eq(0).clone()).append($tds.eq(1).clone()).appendTo($target);
        $row.append($tds.eq(0).clone()).append($tds.eq(1).clone()).append($tds.eq(2).clone()).append($tds.eq(3).clone()).appendTo($target);
    });
    $(".overlay_table .border-right").removeClass("border-right");

    // rename checkboxes for overlay to link to the original and remove class
    $(".overlay_table tbody tr :checkbox").each(function(i, e){
      e.id = "overlaycheck_" + e.id;
      e.classList.remove("pdb_selected");
    });
    // rebind click event
    /*$('.structure_overlay tr').click(function(event) {
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
    });*/
    $('.structure_overlay tr').click(function(event) {
        if (event.target.type !== 'checkbox') {
            $(':checkbox', this).trigger('click');
            if ($(':checkbox', this).length === 0) {
                pdb_id = $(this).attr('id').split("_")[1];
                console.log(pdb_id);
                var checkbox = $(".dataTables_scrollBody").find('#' + pdb_id);
                checkbox.trigger('click');
            }
        }
    });
    $('.structure_overlay tr').css('cursor', 'pointer');
}

function track_scrolling(element) {
    var left = 0;
    var old_left = 0;
    toggle_enabled = true;
    var d = new Date();
    last_change = d.getTime();
    $(element + ' .dataTables_scrollBody').scroll(function () {

        var d = new Date();
        now = d.getTime();

        if (now>(last_change+1000)) {

            // If user scrolls and it's >100px from left, then attach fixed columns overlay
            left = $(element + ' .dataTables_scrollBody').scrollLeft();
            // if (left != old_left) $(".structure_overlay").hide();
            // old_left = left;


            console.log(now, last_change);

            if (left > 100 && toggle_enabled && left != old_left) {
                console.log('SCROLL ', left, now);
                $(".structure_overlay").css({
                    left: left + 'px'
                });
                if ($(".structure_overlay").is(":hidden")){
                    // link check boxes
                    $(".structure_overlay").show();
                }
            } else if (left < 100) {
                if ($(".structure_overlay").is(":visible")){
                // link check boxes
                $(".structure_overlay").hide();
                }
            }
            old_left = left;
            last_change = now;
        }
    });
}

function new_track_scrolling(element) {
    var left = 0;
    var old_left = 0;
    toggle_enabled = true;
    old_height = 0;
    old_top = 0;
    var d = new Date();
    last_change = d.getTime();
    scroll_visible = false;

    // Init the size to begin with.. It is updated once in a while incase of window resizing.
    var dataTables_scrollBody_height = $(element + ' .dataTables_scrollBody')[0].offsetHeight;
    var bodyRect = $(element + ' .dataTables_scroll')[0].getBoundingClientRect(),
        elemRect = $(element + ' .dataTables_scrollBody')[0].getBoundingClientRect(),
        offset = elemRect.top - bodyRect.top;
    $(".structure_overlay").css({
        height: (dataTables_scrollBody_height - 17) + 'px',
        top: offset + 'px'
    });

    $(element + ' .dataTables_scrollBody').scroll(function () {
        var d = new Date();
        now = d.getTime();

        // This isnt dom heavy..
        scrollTop = $(element + ' .dataTables_scrollBody').scrollTop();
        if (scrollTop != old_top) {
            $(".structure_overlay").scrollTop(scrollTop);
            old_top = scrollTop
        }

        left = $(element + ' .dataTables_scrollBody').scrollLeft();
        if ( left<100 && scroll_visible){
            $(".structure_overlay").hide();
            scroll_visible = false;
        } else if ( left>100 && !scroll_visible){
            $(".structure_overlay").show();
            scroll_visible = true;
        }

        if (now > (last_change + 1000)) {

            // If user scrolls and it's >100px from left, then attach fixed columns overlay
            if (left > 100 && toggle_enabled) {

                var dataTables_scrollBody_height = $(element + ' .dataTables_scrollBody')[0].offsetHeight;
                if (dataTables_scrollBody_height != old_height) {
                    var bodyRect = $(element + ' .dataTables_scroll')[0].getBoundingClientRect(),
                        elemRect = $(element + ' .dataTables_scrollBody')[0].getBoundingClientRect(),
                        offset = elemRect.top - bodyRect.top;
                    $(".structure_overlay").css({
                        height: (dataTables_scrollBody_height - 17) + 'px',
                        top: offset + 'px'
                    });
                    old_height = dataTables_scrollBody_height;
                }
                if (!scroll_visible) {
                    $(".structure_overlay").show();
                    scroll_visible = true;
                }
            }
            last_change = now;
        }
    });
    $('.structure_overlay').scroll(function () {
        scrollTop = $('.structure_overlay').scrollTop();
        if (scrollTop != old_top) {
            $(element + ' .dataTables_scrollBody').scrollTop(scrollTop);
            old_top = scrollTop
        }
    });
}
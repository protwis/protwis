function update_text_in_modal() {
  var mode = $('ul#mode_nav').find('li.active').find('a').text().trim();
    // console.log('update modal',mode);

  if (mode=='Single group of structures') {
    var total_selected = $('.pdb_selected:checked', oTable[mode].cells().nodes()).length
    var selected_visible = $('.pdb_selected:checked').length
    var ModalpdbsCountSelector = '#single-crystal-group-pdbs-modal-text';

    if (total_selected==selected_visible) {
      $(ModalpdbsCountSelector).html(total_selected +' structure(s) selected');
    } else {
      $(ModalpdbsCountSelector).html(total_selected +' structure(s) selected ('+(total_selected-selected_visible)+' currently filtered)');
    }
  } else if (mode=='Two groups of structures') {
    group = $('.tableview:visible').attr('group-number');
    if (group) mode = mode + group;
    //#FIXME
    var ModalpdbsCountSelector = '#two-crystal-group-pdbs-modal-'+group+'-text';
    var total_selected = $('.pdb_selected:checked', oTable[mode].cells().nodes()).length
    var selected_visible = $('.pdb_selected:checked:visible').length
    if (total_selected==selected_visible) {
      $(ModalpdbsCountSelector).html(total_selected +' structure(s) selected');
    } else {
      $(ModalpdbsCountSelector).html(total_selected +' structure(s) selected ('+(total_selected-selected_visible)+' currently filtered)');
    }
  }


}

function thisPDB(elem) {
  var mode = $('ul#mode_nav').find('li.active').find('a').text().trim();
  var ReceptorName = $(elem).attr('long');
  var pdbName = $(elem).attr('id');
  // console.log('thisPDB',pdbName);
  if (mode=='Single structure') {
    $('.pdb_selected').not(elem).prop("checked",false);
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
  } else if (mode=='Single group of structures') {
    var pdbs = [];
    $('.pdb_selected:checked', oTable[mode].cells().nodes()).each(function() {
        pdbs.push($(this).attr('id'));
    });
    var pdbsInputSelector = '#single-crystal-group-tab .crystal-pdb';
    var pdbsCountSelector = '#single-crystal-group-tab .crystal-count';
    var ModalpdbsCountSelector = '#single-crystal-group-pdbs-modal-text';

    $(pdbsInputSelector).val(JSON.stringify(pdbs));
    // Update view
    $(pdbsCountSelector).html(pdbs.length);
  }  else if (mode=='Two groups of structures') {

    group = $(elem).closest('.tableview').attr('group-number');
    if (group) mode = mode + group;

    var pdbs = [];
    $('.pdb_selected:checked', oTable[mode].cells().nodes()).each(function() {
        pdbs.push($(this).attr('id'));
    });

    var pdbsInputSelector = '#two-crystal-groups-tab .crystal-group-'+group+'-pdbs';
    var pdbsCountSelector = '#two-crystal-groups-tab .crystal-count-'+group;
    $(pdbsInputSelector).val(JSON.stringify(pdbs));
    // Update view
    $(pdbsCountSelector).html(pdbs.length);
  }
  update_text_in_modal();
}

function resetselection(elem) {
  var mode = $('ul#mode_nav').find('li.active').find('a').text().trim();

  $('.check_all:visible').prop('checked',false);

  if (mode=='Single group of structures') {
    var pdbs = [];

    $('input', oTable[mode].cells().nodes()).prop('checked',false);

    var pdbsInputSelector = '#single-crystal-group-tab .crystal-pdb';
    var pdbsCountSelector = '#single-crystal-group-tab .crystal-count';
    $(pdbsInputSelector).val(JSON.stringify(pdbs));
    $(pdbsCountSelector).html(pdbs.length);
  }  else if (mode=='Two groups of structures') {

    group = $('.tableview:visible').attr('group-number');
    if (group) mode = mode + group;
    var pdbs = [];

    $('input', oTable[mode].cells().nodes()).prop('checked',false);

    var pdbsInputSelector = '#two-crystal-groups-tab .crystal-group-'+group+'-pdbs';
    var pdbsCountSelector = '#two-crystal-groups-tab .crystal-count-'+group;
    $(pdbsInputSelector).val(JSON.stringify(pdbs));
    $(pdbsCountSelector).html(pdbs.length);
  }

  update_text_in_modal();
}

function check_all(elem) {
  var mode = $('ul#mode_nav').find('li.active').find('a').text().trim();
  show_all = $(elem).prop("checked");


  if (mode=='Single group of structures') {
    var pdbs = [];

    // REMOVE EXISITING? Probably not, more logically that filtering adds more
    // $('input', oTable.cells().nodes()).prop('checked',false);

    if (show_all) {
      $('.pdb_selected:visible').prop("checked",true);
    } else {
      $('.pdb_selected:visible').prop("checked",false);
    }

    $('.pdb_selected:checked', oTable[mode].cells().nodes()).each(function() {
        pdbs.push($(this).attr('id'));
    });

    var pdbsInputSelector = '#single-crystal-group-tab .crystal-pdb';
    var pdbsCountSelector = '#single-crystal-group-tab .crystal-count';
    $(pdbsInputSelector).val(JSON.stringify(pdbs));
    // Update view
    $(pdbsCountSelector).html(pdbs.length);
  }  else if (mode=='Two groups of structures') {

    group = $(elem).closest('.tableview').attr('group-number');
    if (group) mode = mode + group;
    var pdbs = [];

    if (show_all) {
      $('.pdb_selected:visible').prop("checked",true);
    } else {
      $('.pdb_selected:visible').prop("checked",false);
    }

    $('.pdb_selected:checked', oTable[mode].cells().nodes()).each(function() {
        pdbs.push($(this).attr('id'));
    });

    var pdbsInputSelector = '#two-crystal-groups-tab .crystal-group-'+group+'-pdbs';
    var pdbsCountSelector = '#two-crystal-groups-tab .crystal-count-'+group;
    $(pdbsInputSelector).val(JSON.stringify(pdbs));
    // Update view
    $(pdbsCountSelector).html(pdbs.length);
  }

  update_text_in_modal();
}

$.fn.dataTable.ext.order['dom-checkbox'] = function  ( settings, col )
  {
      return this.api().column( col, {order:'index'} ).nodes().map( function ( td, i ) {
          return $('input', td).prop('checked') ? '1' : '0';
      } );
  };

var oTable = [];
function showPDBtable(element) {
  var mode = $('ul#mode_nav').find('li.active').find('a').text().trim();
  group = $(element+' .tableview').attr('group-number');
  if (group) mode = mode + group;
  // console.log(element,mode,group);
  if ( ! $.fn.DataTable.isDataTable( element+' .tableview table' ) ) {
      oTable[mode] = $(element+' .tableview table').DataTable({
        'scrollX': true,
        // 'paging': true,
        // 'autoWidth': true,

        scrollY:        '80vh',
        // scrollCollapse: true,
        paging:         false,
        columnDefs: [
          { targets: 'no-sort', orderable: false }
        ],
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
                        { "orderDataType": "dom-checkbox" }
                    ]
      });

      yadcf.init(oTable[mode],
        [
            {
                column_number : 0,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "UniProt",
                filter_reset_button_text: false,
            },
            {
                column_number : 1,
                filter_type: "multi_select",
                select_type: 'select2',
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "Receptor",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
            },
            {
                column_number : 2,
                filter_type: "multi_select",
                select_type: 'select2',
                html_data_type: "text",
                select_type_options: {
                  width: '150px'
                },
                filter_default_label: "Family",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
            },
            {
                column_number : 3,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Class",
                filter_reset_button_text: false,
            },
            {
                column_number : 4,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Species",
                filter_reset_button_text: false,
            },
            {
                column_number : 5,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Method",
                filter_reset_button_text: false,
            },
            {
                column_number : 6,
                filter_type: "multi_select",
                select_type: 'select2',
                select_type_options: {
                  width: '70px'
                },
                filter_default_label: "PDB",
                filter_reset_button_text: false,
            },
            {
                column_number : 7,
                filter_type: "multi_select",
                select_type: 'select2',
                select_type_options: {
                  width: '70px'
                },
                filter_default_label: "Res (Å)",
                filter_reset_button_text: false,
            },
            {
                column_number : 8,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "State",
                filter_reset_button_text: false,

            },
            {
                column_number : 9,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "7TM Open IC (Å)",
                filter_reset_button_text: false,
            },
            {
                column_number : 10,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "G protein",
                filter_reset_button_text: false,
            },
            {
                column_number : 11,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "B arrestin",
                filter_reset_button_text: false,
            },
            {
                column_number : 12,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Fusion",
                filter_reset_button_text: false,
            },
            {
                column_number : 13,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Antibody",
                filter_reset_button_text: false,
            },
            {
                column_number : 14,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Ligand",
                filter_reset_button_text: false,
            },
            {
                column_number : 15,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Ligand function",
                filter_reset_button_text: false,
            },
            {
                column_number : 16,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Ligand type",
                filter_reset_button_text: false,
            },
        ],
        {
            cumulative_filtering: false
        }
    );

    yadcf.exResetAllFilters(oTable[mode]);

    oTable[mode].on('draw.dt', function (e, oSettings) {
        update_text_in_modal();
        create_overlay(element+' .structure_selection');
    });

    $(element+' .dataTables_scrollBody').append('<div class="structure_overlay"><table id="overlay_table" class="overlay_table row-border text-center compact dataTable no-footer text-nowrap"><tbody></tbody></table></div>');

    $(element+" .structure_overlay").hide();
    
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

    create_overlay(element+' .structure_selection');
    track_scrolling(element);

  };
}

function create_overlay(table_id) {
    // This function fires upon filtering, to update what rows to show as an overlay
    $(".overlay_table tbody tr").remove(); 
    var $target = $(".overlay_table tbody");
    $(table_id+" tbody tr").each(function() {
        var $tds = $(this).children(),
            $row = $("<tr></tr>");
        // $row.append($tds.eq(0).clone()).append($tds.eq(1).clone()).appendTo($target);
        $row.append($tds.eq(0).clone()).append($tds.eq(1).clone()).append($tds.eq(2).clone()).appendTo($target);
    });
    $(".overlay_table .border-right").removeClass("border-right");
}

function track_scrolling(element) {
    var left = 0;
    var old_left = 0;
    toggle_enabled = true;
    $(element+' .dataTables_scrollBody').scroll(function(){
        // If user scrolls and it's >100px from left, then attach fixed columns overlay
        left = $(element+' .dataTables_scrollBody').scrollLeft();
        if (left!=old_left) $(".structure_overlay").hide();
        old_left = left;

        if (left>100 && toggle_enabled) {
            $(".structure_overlay").css({ left: left+'px' });
            if ($(".structure_overlay").is(":hidden")) $(".structure_overlay").show();
        }
    });
}
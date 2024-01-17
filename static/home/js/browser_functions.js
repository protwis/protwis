/*global showAlert*/
/*eslint no-undef: "error"*/

function superposition(oTable, columns, site, source='gpcr', structure_column_index, hidden_columns=[]) {
    // oTable: DataTable object of table of entries
    // columns: Column indeces of oTable to be extracted to build table for reference selection. First column has to be structure/model string used for superposition workflow
    // site: Structure browser or Homology model browser (add new logic when expanding to new sites)
    // source: Through which db is the call launced: gpcr, gprot
    // structure_column_index: column index of the object identifier, int
    // hidden_columns: array of indexes from the columns parameter that should be hidden in the modal 

    ClearSelection('targets');
    ClearSelection('reference');

    var checked_data = oTable.rows('.alt_selected').data();
    if (checked_data.length===0) {
        showAlert("No entries selected for superposition", "danger");
        return 0;
    }
    else if (checked_data.length > 50) {
        showAlert("Maximum number of selected entries is 50", "warning");
        return 0;
    }
    var selected_ids = [];
    var selection_type = '';

    if (site==='structure_browser') {
        for (i = 0; i < checked_data.length; i++) {
            var div = document.createElement("div");
            div.innerHTML = checked_data[i][structure_column_index];
            if (typeof div.innerText !== "undefined") {
                selected_ids.push(div.innerText.replace(/\s+/g, ''));
            } else {
                selected_ids.push(div.textContent.replace(/\s+/g, ''));
            }
        }
        selection_type = 'structure_many';
    }

    else if (site==='homology_model_browser') {
        for (i = 0; i < checked_data.length; i++) {
            var div = document.createElement("div");
            div.innerHTML = checked_data[i][structure_column_index];
            var state = checked_data[i][3];
            if (checked_data[i][4]==='Yes') {
                div.innerHTML = checked_data[i][structure_column_index];
                if (typeof div.innerText !== "undefined") {
                    selected_ids.push(div.innerText.replace(/\s+/g, '')+"_refined");
                } else {
                    selected_ids.push(div.textContent.replace(/\s+/g, '')+"_refined");
                }
            }
            else {
                if (typeof div.innerText !== "undefined") {
                    selected_ids.push(div.innerText.replace(/\s+/g, '')+"_"+state);
                } else {
                    selected_ids.push(div.textContent.replace(/\s+/g, '')+"_"+state);
                }
            }
        }
        selection_type = 'structure_model_many';
    }

    else if (site==='complex_models') {
        for (i = 0; i < checked_data.length; i++) {
            var div = document.createElement("div");
            div.innerHTML = checked_data[i][structure_column_index];
            if (typeof div.innerText !== "undefined") {
                selected_ids.push(div.innerText.replace(/\s+/g, ''));
            } else {
                selected_ids.push(div.textContent.replace(/\s+/g, ''));
            }
        }
        selection_type = 'structure_many';
    }

    else if (site == "g_protein_structure_browser"){
      for (i = 0; i < checked_data.length; i++) {
          var div = document.createElement("div");
          div.innerHTML = checked_data[i][structure_column_index];
          if (typeof div.innerText !== "undefined") {
              selected_ids.push(div.innerText.replace(/\s+/g, ''));
          } else {
              selected_ids.push(div.textContent.replace(/\s+/g, ''));
          }
      }
      selection_type = 'signprot_many';
    }
    // add new logic here for new site

    $('#superposition_modal_table tbody').empty();
    var modal = document.getElementById("superposition-modal");
    var span = document.getElementById("close_superposition_modal");

    modal.style.display = "block";
    span.onclick = function() {
        modal.style.display = "none";
    }
    window.onclick = function(event) {
        if (event.target == modal) {
            modal.style.display = "none";
        }
    }
    var checked_data = oTable.rows('.alt_selected').data();
    // Create modal table HTML
    for (i=0; i<checked_data.length; i++) {
        var div = document.createElement("div");
        row = document.createElement('tr');
        var checkbox = document.createElement('input');
        checkbox.type = "checkbox";
        row.appendChild(checkbox);
        columns.forEach(function(column) {
            div.innerHTML = checked_data[i][column];
            cell = document.createElement('td');
            if (hidden_columns.includes(column)) {
                cell.style.display = "none";
            }
            textnode = document.createTextNode(div.innerText.replace(/\s+/g, ''));
            cell.appendChild(textnode);
            row.appendChild(cell);
        });
        $('#superposition_modal_table tbody').append(row)
    }
    $("#superposition_modal_table tbody tr").click(function() {
        var ref_id = '';

        if (site==='structure_browser') {
            ref_id = $(this).children().eq(columns.indexOf(structure_column_index)+1).text();
            AddToSelection('reference', 'structure', ref_id);
        }

        else if (site==='homology_model_browser') {
            if ($(this).children().eq(8).text()==='Yes') {
                ref_id = $(this).children().eq(11).text()+"_refined";
                AddToSelection('reference', 'structure', ref_id);
            }
            else {
                var state = $(this).children().eq(7).text();
                ref_id = $(this).children().eq(columns.indexOf(structure_column_index)+1).text()+"_"+state;
                AddToSelection('reference', 'structure_model', ref_id);
            }
        }

        else if (site==='complex_models') {
            AddToSelection('reference', 'structure', $(this).children().eq(columns.indexOf(structure_column_index)+1).text());
        }

        else if (site==='g_protein_structure_browser') {
            ref_id = $(this).children().eq(columns.indexOf(structure_column_index)+1).text();
            AddToSelection('reference', 'signprot', ref_id);
        }
        // Add logic here for new site

        // Bulk add targets to selection
        selected_ids.splice(selected_ids.indexOf(ref_id), 1);
        AddToSelection('targets', selection_type, selected_ids.join(","));

        $(this).children(':first').prop("checked",true);
        if (source=='gpcr'){
          window.location.href = '/structure/superposition_workflow_index#keepselection';
        } else {
          window.location.href = '/structure/superposition_workflow_gprot_index#keepselection';
        }

    })
};

function direct_superposition(oTable, source="gpcr", structure_column_index, selection_type) {
    var data_type = document.getElementById("superpose_template_btn").dataset.type;
    var checked_data = oTable.rows('.alt_selected').data();
    var selected_ids = [];
    if (checked_data.length == 0) {
        showAlert("No entries selected", "danger")
        return 0;
    }
    if (data_type == "reference") {
        if (checked_data.length > 1) {
            showAlert("Only select one entry as reference", "danger");
            return 0;
        }
        var div = document.createElement("div");
        div.innerHTML = checked_data[0][structure_column_index];
        if (typeof div.innerText !== "undefined") {
            AddToSelection('reference', selection_type, div.innerText.replace(/\s+/g, ''));
        } else {
            AddToSelection('reference', selection_type, div.textContent.replace(/\s+/g, ''));
        }
    }
    else if (data_type == "targets") {
        if (checked_data.length > 50) {
            showAlert("Maximum number of selected targets is 50");
            return 0;
        }
        for (i = 0; i < checked_data.length; i++) {
            var div = document.createElement("div");
            div.innerHTML = checked_data[i][structure_column_index];
            selected_ids.push(div.textContent.replace(/\s+/g, ''));
        }
        AddToSelection('targets', selection_type+"_many", selected_ids.join(","));
    }
    if (source=="gprot") {
        window.location.href = '/structure/superposition_workflow_gprot_index#keepselection';
    } else {
        window.location.href = '/structure/superposition_workflow_index#keepselection';
    }
}

function updateSuperpositionButtonConfiguration() {
    const button = document.getElementById("superpose_template_btn");

    switch (window.location.hash) {
        case "#keepselectionreference":
            // Update the button for #keepselectionreference
            button.innerText = "Add reference to superposition";
            button.id = "superpose_template_btn";
            button.classList.remove("btn-default");
            button.classList.add("btn-primary");
            button.dataset.type = "reference";
            break;
        case "#keepselectiontargets":
            // Update the button for #keepselectiontargets
            button.innerText = "Add targets to superposition";
            button.id = "superpose_template_btn";
            button.classList.remove("btn-default");
            button.classList.add("btn-primary");
            button.dataset.type = "targets";
            break;
        default:
            // Default configuration (if no recognized hash)
            button.innerText = "Superposition";
            button.id = "superpose_btn";
            button.classList.remove("btn-primary");
            button.classList.add("btn-default");
            button.dataset.type = "default";
          break;
    }
}

function Sort(clicked_th, table1, table2) {
    if (clicked_th.parent().parent().parent().attr('class').includes('scrollable')) {
        var table = table1;
    }
    else if (clicked_th.parent().parent().parent().attr('class').includes('frozen')) {
        var table = table2;
    }
    if (clicked_th.attr('class').includes('asc')) {
        table.order([clicked_th.index(),'asc']);
    }
    else if (clicked_th.attr('class').includes('desc')) {
        table.order([clicked_th.index(),'desc']);
    }
    table.draw();
}

function ClearSelection(selection_type) {
    $.ajax({
        'url': '/common/clearselection',
        'data': {
            selection_type: selection_type
        },
        'type': 'GET',
        'async': false,
        'success': function (data) {
            $("#selection-" + selection_type).html(data);
        }
    });
}

function AddToSelection(selection_type, selection_subtype, selection_id) {
    $.ajax({
        'url': '/common/addtoselection',
        'data': {
            selection_type: selection_type,
            selection_subtype: selection_subtype,
            selection_id: selection_id
        },
        'type': 'GET',
        'async': false,
        'success': function(data) {
            $("#selection-" + selection_type).html(data);
        },
    });
}

function RemoveFromSelection(selection_type, selection_subtype, selection_id) {
    $.ajax({
        'url': '/common/removefromselection',
        'data': {
            selection_type: selection_type,
            selection_subtype: selection_subtype,
            selection_id: selection_id
        },
        'type': 'GET',
        'async': false,
        'success': function(data) {
            $("#selection-" + selection_type).html(data);
        },
    });
}

function CheckSelection(selection_type) {
    var result = null;

    $.ajax({
        'url': '/common/checkselection',
        'data': {
            selection_type: selection_type
        },
        'type': 'GET',
        'dataType': 'json',  // Expecting JSON response from the server
        'async': false,
        'success': function(response) {
            result = response.total;
        },
        'error': function(error) {
            console.error("An error occurred:", error);
        }
    });

    return result;
}


function select_all(e, table1_id, table2_id) {
    $('#loading_div').show();
    var rows = $("#"+table1_id).dataTable().$('tr', {filter:"applied"});
    var checkedStatus = $(e).prop("checked");
    $(rows).each(function() {
        if (this.hidden==false) {
            $(this).prop("checked", checkedStatus);
            $(this).eq(0).toggleClass('alt_selected', checkedStatus);
            // if (checkedStatus==true) {
            $(this).find(":checkbox").prop("checked", checkedStatus)
            // }
            // $(this).find(":checkbox").trigger('click');
        }
    });
    $('#loading_div').hide();
};

function assign_to_row1(table1_id, table2_id){
  $('#'+table1_id+' > tbody > tr').click(function(event) {
    if (event.target.type !== 'checkbox') {
        $(':checkbox', this).trigger('click');
        $(':checkbox', $("#"+table2_id+" > tbody > tr[model_id=" + $(this).attr('model_id') + "]")).trigger('click');
        $(this).eq(0).toggleClass('alt_selected');
        console.log($(this).eq(0));
    }
    else {
        $(this).eq(0).toggleClass('alt_selected');
    }
  });
}

function assign_to_row2(table1_id, table2_id){
  $('#'+table2_id+' > tbody > tr').click(function(event) {
    if (event.target.type !== 'checkbox') {
        $(':checkbox', this).trigger('click');
        $(':checkbox', $("#"+table1_id+" > tbody > tr[model_id=" + $(this).attr('model_id') + "]")).trigger('click');
        $("#"+table1_id+" > tbody > tr[model_id=" + $(this).attr('model_id') + "]").eq(0).toggleClass('alt_selected');
        console.log($("#"+table1_id+" > tbody > tr[model_id=" + $(this).attr('model_id') + "]").eq(0));
    }
    else {
        $("#"+table1_id+" > tbody > tr[model_id=" + $(this).attr('model_id') + "]").eq(0).toggleClass('alt_selected');
    }
  });
}

function match_rowheight(table1_id, table2_id) {
    $("#"+table2_id).find('tr').each(function(i, e) {
        $("#"+table1_id).find('tr').eq(i).height($(this).height());
    });
    $("#"+table1_id).find('tr').each(function(i, e) {
        $("#"+table2_id).find('tr').eq(i).height($(this).height());
    });
}

var trigger_draw = true

function catch_filter(table1_id, table2_id, oTable1, oTable2) {
    if (trigger_draw==true) {
        $('#'+table2_id).on( 'draw.dt', function (e,oSettings) {
            $('#loading_div').show();
            update_tables(table1_id, table2_id, 2);
            var scrollable_height = $('.yadcf-datatables-table--'+table2_id+' > thead > tr').eq(1).eq(0).height();
            $('.yadcf-datatables-table--'+table1_id+' > thead > tr').eq(1).css('height', scrollable_height);
            $('#loading_div').hide();
        });
        $('#'+table1_id).on( 'draw.dt', function (e,oSettings) {
            $('#loading_div').show();
            update_tables(table1_id, table2_id, 1);
            var frozen_height = $('.yadcf-datatables-table--'+table1_id+' > thead > tr').eq(1).eq(0).height();
            $('.yadcf-datatables-table--'+table2_id+' > thead > tr').eq(1).css('height', frozen_height);
            trigger_draw = false;
            oTable2.draw();
            trigger_draw = true;
            $('#loading_div').hide();
        });
    }
}

function update_tables(table1_id, table2_id, source_table) {
    if (source_table===2) {
        $('#'+table1_id+' > tbody > tr').prop('hidden', true);
        ids = Array();
        var rows = $("#"+table2_id).dataTable().$('tr', {filter:"applied"});
        rows.each(function(i,e) {
            $("#"+table1_id+" > tbody > tr[model_id=" + $(e).attr('model_id') + "]").attr('class', $(e).attr('class'));
            $("#"+table1_id+" > tbody > tr[model_id=" + $(e).attr('model_id') + "]").prop('hidden', false);
        });
    }
    else if (source_table===1) {
        $('#'+table2_id+' > tbody > tr').prop('hidden', true);
        ids = Array();
        var rows = $("#"+table1_id).dataTable().$('tr', {"filter":"applied"});
        rows.each(function(i,e) {
            $("#"+table2_id+" > tbody > tr[model_id=" + $(e).attr('model_id') + "]").attr('class', $(e).attr('class'));
            $("#"+table2_id+" > tbody > tr[model_id=" + $(e).attr('model_id') + "]").prop('hidden', false);
        });
    }
}

//Match row highlight in other table
function match_row_highlight(table1_id, table2_id) {
    var prev_index = -1;
    $('#'+table2_id+' > tbody > tr').mouseover(function(){
        if (prev_index!==$(this).index()){
            $(this).css('background-color','rgb(245,245,245)');
            $("#"+table1_id+" > tbody > tr[model_id=" + $(this).attr('model_id') + "]").css('background-color','rgb(245,245,245)');
        }
        prev_index = $(this).index();
    });
    $('#'+table1_id+' > tbody > tr').mouseover(function(){
        if (prev_index!==$(this).index()){
            $(this).css('background-color','rgb(245,245,245)');
            $("#"+table2_id+" > tbody > tr[model_id=" + $(this).attr('model_id') + "]").css('background-color','rgb(245,245,245)');
        }
        prev_index = $(this).index();
    });
    $('#'+table2_id+' > tbody > tr').mouseout(function(){
        if ($(this).hasClass('odd')){
            $(this).css('background-color','rgb(249,249,249)');
            $("#"+table1_id+" > tbody > tr[model_id=" + $(this).attr('model_id') + "]").css('background-color','rgb(249,249,249)');
        }
        else if ($(this).hasClass('even')){
            $(this).css('background-color','rgb(255,255,255)');
            $("#"+table1_id+" > tbody > tr[model_id=" + $(this).attr('model_id') + "]").css('background-color','rgb(255,255,255)');
        }
    });
    $('#'+table1_id+' > tbody > tr').mouseout(function(){
        if ($(this).hasClass('odd')){
            $(this).css('background-color','rgb(249,249,249)');
            $("#"+table2_id+" > tbody > tr[model_id=" + $(this).attr('model_id') + "]").css('background-color','rgb(249,249,249)');
        }
        else if ($(this).hasClass('even')){
            $(this).css('background-color','rgb(255,255,255)');
            $("#"+table2_id+" > tbody > tr[model_id=" + $(this).attr('model_id') + "]").css('background-color','rgb(255,255,255)');
        }
    });
}

//Matching scrolls
function match_scroll_position() {
    var tables = document.querySelectorAll(".dataTables_scrollBody");
    var frozen_table = tables[0];
    var scrollable_table = tables[1];

    var isLeftScrollTopCalled = false;
    $(frozen_table).scroll(function (e) {
        if (isRightScrollTopCalled) {
            return isRightScrollTopCalled = false;
        }
        $(scrollable_table).scrollTop($(this).scrollTop());
        isLeftScrollTopCalled = true;
    });

    var isRightScrollTopCalled = false;
    $(scrollable_table).scroll(function (e) {
        if (isLeftScrollTopCalled) {
            return isLeftScrollTopCalled = false;
        }
        $(frozen_table).scrollTop($(this).scrollTop());
        isRightScrollTopCalled = true;
    });
}

function copyToClipboard(array, delimiter, data_name, powertip_object=false) {
    var link = array;
    var out = "";
    link.each(function() {
        var ele = $(this).attr("href").split("/");
        out+=ele[ele.length-1]+delimiter;
    });
    if (out.length===0) {
        window.alert("No entries selected for copying");
        return 0;
    }
    var textArea = document.createElement("textarea");
    textArea.value = out;
    document.body.appendChild(textArea);
    textArea.focus();
    textArea.select();
    try {
        var successful = document.execCommand("copy");
        var msg = successful ? "Successful" : "Unsuccessful";
        if (powertip_object!==false) {
            $.powerTip.hide();
            powertip_object.data("powertipjq", $([
                "<p>Copied to clipboard!</p>"
                ].join("\n")));
            powertip_object.powerTip("show");
            setTimeout(function() {
            powertip_object.data("powertipjq", $([
                "<p>Export "+data_name+"</p>"
                ].join("\n")));
            },1000);
        }
    } catch (err) {
        window.alert("Oops, unable to copy");
    }
    document.body.removeChild(textArea);
}

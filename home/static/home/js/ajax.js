function AddToSelection(selection_type, selection_subtype, selection_id) {
    $.ajax({
        'url': '/common/addtoselection',
        'data': {
            selection_type: selection_type,
            selection_subtype: selection_subtype,
            selection_id: selection_id
        },
        'type': 'GET',
        'success': function(data) {
            $("#selected-data").html(data);
        }
    });
}

function ClearSelection() {
    $.ajax({
        'url': '/common/clearselection',
        'type': 'GET',
        'success': function(data) {
            $("#selected-data").html(data);
        }
    });
}
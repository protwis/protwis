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
        'success': function(data) {
            $("#selection-" + selection_type).html(data);
        }
    });
}

function ClearSelection(selection_type) {
    $.ajax({
        'url': '/common/clearselection',
        'data': {
            selection_type: selection_type
        },
        'type': 'GET',
        'success': function(data) {
            $("#selection-" + selection_type).html(data);
        }
    });
}

function ToggleFamilyTreeNode(action, type_of_selection, node_id, tree_indent_level) {
    $.ajax({
        'url': '/common/togglefamilytreenode',
        'data': {
            action: action,
            type_of_selection: type_of_selection,
            node_id: node_id,
            tree_indent_level: tree_indent_level
        },
        'type': 'GET',
        'success': function(data) {
            $("#" + node_id).html(data);
        }
    });
}
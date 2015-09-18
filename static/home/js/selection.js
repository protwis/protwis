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

function SelectFullSequence(selection_type) {
    $.ajax({
        'url': '/common/selectfullsequence',
        'data': {
            selection_type: selection_type
        },
        'type': 'GET',
        'success': function(data) {
            $("#selection-" + selection_type).html(data);
        },
    });
}

function SelectAlignableSegments(selection_type) {
    $.ajax({
        'url': '/common/selectalignablesegments',
        'data': {
            selection_type: selection_type
        },
        'type': 'GET',
        'success': function(data) {
            $("#selection-" + selection_type).html(data);
        },
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

function SelectionAnnotation(protein_source) {
    $.ajax({
        'url': '/common/selectionannotation',
        'data': {
            protein_source: protein_source
        },
        'type': 'GET',
        'success': function(data) {
            $("#filters-annotation").html(data);
        }
    });
}

function SelectionSpeciesPredefined(species) {
    $.ajax({
        'url': '/common/selectionspeciespredefined',
        'data': {
            species: species
        },
        'type': 'GET',
        'success': function(data) {
            $("#filters-species").html(data);
        }
    });
}

function SelectionSpeciesToggle(species_id) {
    $.ajax({
        'url': '/common/selectionspeciestoggle',
        'data': {
            species_id: species_id
        },
        'type': 'GET',
        'success': function(data) {
            $("#filters-species-selector").html(data);
        }
    });
}

function ExpandSegment(segment_id, position_type, scheme) {
    $.ajax({
        'url': '/common/expandsegment',
        'data': {
            segment_id: segment_id,
            position_type: position_type,
            numbering_scheme: (typeof scheme === 'undefined') ? false : scheme
        },
        'type': 'GET',
        'success': function(data) {
            $("#segment-generic-numbers").html(data);
        }
    });
}

function SelectionSchemesPredefined(numbering_schemes) {
    $.ajax({
        'url': '/common/selectionschemespredefined',
        'data': {
            numbering_schemes: numbering_schemes
        },
        'type': 'GET',
        'success': function (data) {
            $("#filters-schemes").html(data);
        }
    });
}

function SelectionSchemesToggle(numbering_scheme_id) {
    $.ajax({
        'url': '/common/selectionschemestoggle',
        'data': {
            numbering_scheme_id: numbering_scheme_id
        },
        'type': 'GET',
        'success': function (data) {
            $("#filters-schemes").html(data);
        }
    });
}

function SetTreeSelection(option_no, option_id) {
    console.log(option_no)
    console.log(option_id)
    $.ajax({
        'url': '/common/settreeselection',
        'data': {
            option_no: option_no,
            option_id: option_id
        },
        'type': 'GET',
        'success': function (data) {
            $("#tree_buttons").html(data);
        }
    });
}

function SelectResidueFeature(selection_type, selection_subtype, selection_id, feature) {
    $.ajax({
        'url': '/common/selectresiduefeature',
        'data': {
            selection_type: selection_type,
            selection_subtype: selection_subtype,
            selection_id: selection_id,
            feature: feature
        },
        'type': 'GET',
        'success': function(data) {
            $("#selection-" + selection_type).html(data);
        },
    });
}
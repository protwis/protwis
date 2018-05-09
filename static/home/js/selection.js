$(function () {
    $('#selection-button').click(function () {
        toggleButtonClass('selection-button');
    });
});

function toggleButtonClass(button_id) {
    $('#'+button_id).toggleClass('active')
}

function ToggleSegments() {
    $('#sequence_segments').slideToggle("fast");
}

function ToggleResidueSets() {
    $('#residue_sets').slideToggle("fast");
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
        'async': false,
        'success': function (data) {
            $("#selection-" + selection_type).html(data);
        }
    });
}

function SelectRange(selection_type, selection_subtype, range_start, range_end) {
    $.ajax({
        'url': '/common/selectrange',
        'data': {
            selection_type: selection_type,
            selection_subtype: selection_subtype,
            range_start: range_start,
            range_end: range_end,
        },
        'type': 'GET',
        'async': false,
        'success': function (data) {
            $("#selection-" + selection_type).html(data);
        },
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

function SelectFullSequenceGprotein(selection_type, protein_type) {
    $.ajax({
        'url': '/common/selectfullsequence',
        'data': {
            selection_type: selection_type,
            protein_type: protein_type
        },
        'type': 'GET',
        'success': function(data) {
            $("#selection-" + selection_type).html(data);
        },
    });
}

function SelectStructuredGprotein(selection_type, protein_type) {
    $.ajax({
        'url': '/common/selectalignablesegments',
        'data': {
            selection_type: selection_type,
            protein_type: protein_type
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


function SelectAlignableResidues(selection_type) {
    $.ajax({
        'url': '/common/selectalignableresidues',
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

function SelectionGproteinPredefined(g_protein, pref) {
    $.ajax({
        'url': '/common/selectiongproteinpredefined',
        'data': {
            g_protein: g_protein,
            preferred: pref
        },
        'type': 'GET',
        'success': function (data) {
            if (pref == true) {
                $("#filters-pref-gproteins").html(data);
            }
            else {
                $("#filters-gproteins").html(data);
            }
        }
    });
}

function SelectionGproteinToggle(g_protein_id, pref) {
    $.ajax({
        'url': '/common/selectiongproteintoggle',
        'data': {
            g_protein_id: g_protein_id,
            preferred: pref
        },
        'type': 'GET',
        'success': function (data) {
            if (pref == true) {
                $("#filters-pref-gproteins-selector").html(data);
            }
            else {
                $("#filters-gproteins-selector").html(data);
            }
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

function ReadTargetsForm(form) {
    //Dealing with csrf token
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
    //Actual post
    var fd = new FormData(form);
    $.ajax({
        type: 'POST',
        url: '/common/targetformread',
        data: fd,
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

function SetTreeSelection(option_no, option_id) {
    $.ajax({
        'url': '/common/settreeselection',
        'data': {
            option_no: option_no,
            option_id: option_id
        },
        'type': 'GET',
        'success': function (data) {
            $("#tree-options").html(data);
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

function AddResidueGroup(selection_type) {
    $.ajax({
        'url': '/common/addresiduegroup',
        'data': {
            selection_type: selection_type
        },
        'type': 'GET',
        'success': function(data) {
            $("#selection-" + selection_type).html(data);
        },
    });
}

function SelectResidueGroup(selection_type, group_id) {
    $.ajax({
        'url': '/common/selectresiduegroup',
        'data': {
            selection_type: selection_type,
            group_id: group_id
        },
        'type': 'GET',
        'success': function(data) {
            $("#selection-" + selection_type).html(data);
        },
    });
}

function RemoveResidueGroup(selection_type, group_id) {
    $.ajax({
        'url': '/common/removeresiduegroup',
        'data': {
            selection_type: selection_type,
            group_id: group_id
        },
        'type': 'GET',
        'success': function(data) {
            $("#selection-" + selection_type).html(data);
        },
    });
}

function SetGroupMinMatch(selection_type, group_id, min_match) {
    $.ajax({
        'url': '/common/setgroupminmatch',
        'data': {
            selection_type: selection_type,
            group_id: group_id,
            min_match: min_match
        },
        'type': 'GET',
        'success': function(data) {
            $("#selection-" + selection_type).html(data);
        },
    });
}

function ReadDefinitionFromFile(form, url) {
    //Dealing with csrf token
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
    //Actual post
    var fd = new FormData();
    fd.append('xml_file', form, form.name)
    $.ajax({
        type: 'POST',
        'url': url,
        data: fd,
        cache: false,
        processData: false,
        contentType: false,
        'success': function (data) {
            $("#selection-segments").html(data);
        },
    }).fail(function (jqXHR, textStatus, error) {
        alert("Request failed: " + textStatus + error);
    });
    return false;
}

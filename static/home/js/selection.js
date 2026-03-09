/*global showAlert*/


function resetButtonState() {
    $('.has-spinner.active').removeClass('active').prop('disabled', false);
}

$(document).ready(function() {
    resetButtonState();
    updateSegmentAvailability();  // sync GAIN + D1 state with current targets on load
    syncSegmentButtonSelectionState();
});

window.addEventListener('pageshow', function(event) {
    if (event.persisted) {
        resetButtonState();
    }
});



$(function () {

    $('.has-spinner').on('click', function(e) {
        if ($(this).is('[disabled]')) {
            e.preventDefault();
            return;
        }
        $(this).toggleClass('active');
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
            updateSegmentAvailability();
            syncSegmentButtonSelectionState();
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
            updateSegmentAvailability();
            syncSegmentButtonSelectionState();
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
            updateSegmentAvailability();
            syncSegmentButtonSelectionState();
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
            updateSegmentAvailability();
            syncSegmentButtonSelectionState();
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
            updateSegmentAvailability();
            syncSegmentButtonSelectionState();
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

var verifiedResponse;
function VerifyMinimumSelection(selection_type, minimum) {
    $.ajax({
        'url': '/common/verifyminimumselection',
        'data': {
            selection_type: selection_type,
            minimum: minimum
        },
        'type': 'GET',
        'async': false,
        'success': function(data) {
            if (data == "True"){
              verifiedResponse = true;
            } else {
              showAlert("Please select at least " + minimum + " " + selection_type, "danger");
              verifiedResponse = false;
            }
        },
        'error': function() {
            showAlert("Something went wrong, please try again later", "", "danger");
            return false;
        },
    });

    if (!verifiedResponse)
      toggleButtonClass('selection-button');
    return verifiedResponse;
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

function VerifyMinSegmentSelection() {
    if ($("#selection-segments .target-selection").length === 0){
      showAlert("Please select at least 1 segment item to continue", "danger");
      // Remove active button class => stop spinner after short timeout
      setTimeout(function(){ $(".has-spinner.active").removeClass("active"); }, 1000);
      return false;
    }
    return true;
}

function ToggleSubGroup(id) {
    var el = document.getElementById(id);
    if (!el) return;

    // Toggle visibility
    var isHidden = (el.style.display === 'none' || !el.style.display);
    el.style.display = isHidden ? 'block' : 'none';

    // Find the arrow icon inside the toggle button
    // `event.currentTarget` is the <a> element that was clicked
    var arrow = event.currentTarget.querySelector('.glyphicon');

    if (arrow) {
        if (isHidden) {
            // Now opened → show arrow up
            arrow.classList.remove('glyphicon-arrow-down');
            arrow.classList.add('glyphicon-arrow-up');
        } else {
            // Now closed → show arrow down
            arrow.classList.remove('glyphicon-arrow-up');
            arrow.classList.add('glyphicon-arrow-down');
        }
    }
}

function updateSegmentAvailability() {
    $.ajax({
        url: '/common/checkselectionstatus/',
        type: 'GET',
        success: function(data) {

            let hasB2 = data.has_class_b2;
            let hasD1 = data.has_class_d1;

            // ============================
            //      GAIN DOMAIN (B2)
            // ============================
            let gainHeading = $("#gain-heading");              // make sure you added this id in HTML
            let gainButtons = $("#gain-container .gain-segment");

            if (hasB2) {
                // enable buttons
                gainButtons.removeClass("segment-disabled");

                // restore normal heading
                if (gainHeading.length) {
                    gainHeading.text("Class B2 – GAIN domain");
                }

            } else {
                // disable buttons
                gainButtons.addClass("segment-disabled");

                // show info in heading
                if (gainHeading.length) {
                    gainHeading.text("Class B2 – GAIN domain (No Class B2 in selection)");
                }

                // remove any selected GAIN segments on the right side
                removeGainSegmentsFromSelection();
            }


            // ============================
            //      CLASS D1 DOMAIN
            // ============================
            var d1Buttons = $("[data-segment-domain='D1']");

            if (!hasD1) {
                d1Buttons.addClass("segment-disabled");
                removeD1SegmentsFromSelection();
            } else {
                d1Buttons.removeClass("segment-disabled");
            }


            // ==========================================
            //   WHOLE "CLASS-SPECIFIC SEGMENTS" DROPDOWN
            // ==========================================
            let classBtn   = $("#class-specific-toggle");
            let classPanel = $("#class-specific");

            let hasAnyClassSpecific = hasB2 || hasD1;

            if (!hasAnyClassSpecific) {
                // Disable the toggle and close the panel
                classBtn.addClass("sub-dropdown-disabled disabled");
                classPanel.hide();

                // Replace label with info message
                classBtn.html(
                    'Class-specific segments ' +
                    '<span style="font-size:85%;">(no class-specific targets in selection)</span>' +
                    '<span class="glyphicon glyphicon-arrow-down pull-right"></span>'
                );

            } else {
                // Enable the toggle
                classBtn.removeClass("sub-dropdown-disabled disabled");

                // Restore normal label (don’t auto-open)
                classBtn.html(
                    'Class-specific segments ' +
                    '<span class="glyphicon glyphicon-arrow-down pull-right"></span>'
                );
            }
        }
    });
}



function removeGainSegmentsFromSelection() {
    $("#selection-segments .target-selection[data-segment-domain='GAIN']").each(function() {
        RemoveFromSelection('segments',
            $(this).data("selection-subtype"),
            $(this).data("segment-id"));
    });
}

function removeD1SegmentsFromSelection() {
    $("#selection-segments .target-selection[data-segment-domain='D1']").each(function() {
        RemoveFromSelection('segments',
            $(this).data("selection-subtype"),
            $(this).data("segment-id"));
    });
}

// Toggle a segment when clicking its left-hand button
function ToggleSegmentSelection(btn) {
    var $btn    = $(btn);
    var segId   = $btn.data("segment-id");
    var domain  = $btn.data("segment-domain");
    var subtype = $btn.data("selection-subtype");   // e.g. helix/loop/terminus

    // Check if this segment is already selected (exists on the right side)
    var $existing = $("#selection-segments .target-selection[data-segment-id='" + segId + "'][data-segment-domain='" + domain + "']");

    if ($existing.length) {
        // Already selected -> remove it
        var existingSubtype = $existing.data("selection-subtype") || subtype || 'segment';
        RemoveFromSelection('segments', existingSubtype, segId);
    } else {
        // Not yet selected -> add it
        AddToSelection('segments', subtype, segId);
    }
}

// Sync the "selected" outline on left buttons with the right-hand selection list
function syncSegmentButtonSelectionState() {
    // Reset previous selection state
    $(".segment-btn").removeClass("segment-selected");
    $(".btn-group").removeClass("segment-group-selected");

    // For each selected item on the right side
    $("#selection-segments .target-selection[data-segment-id]").each(function() {
        var segId  = $(this).data("segment-id");
        var domain = $(this).data("segment-domain");

        // Find label button on left
        var selector = ".segment-btn[data-segment-id='" + segId + "'][data-segment-domain='" + domain + "']";

        var $labelBtn = $(selector).not(".arrow-btn");

        if ($labelBtn.length) {
            // add highlight FOR THE WHOLE BUTTON GROUP
            $labelBtn.closest(".btn-group").addClass("segment-group-selected");

            // (optional) still mark the label button as "selected"
            $labelBtn.addClass("segment-selected");
        }
    });
}

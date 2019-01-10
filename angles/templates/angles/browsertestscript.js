window.zoomHeatmap = {};

function toggleFullScreen(fullScreenElement) {
    if (!document.mozFullScreen && !document.webkitFullScreen) {
        if (fullScreenElement.mozRequestFullScreen) {
            fullScreenElement.mozRequestFullScreen();
        } else {
            fullScreenElement.webkitRequestFullScreen(Element.ALLOW_KEYBOARD_INPUT);
        }
    } else {
        if (document.mozCancelFullScreen) {
            document.mozCancelFullScreen();
        } else {
            document.webkitCancelFullScreen();
        }
    }
}

function hidePopovers() {
    $('.popover').each(function(){
        $(this).remove();
    });
}

function HSVtoRGB(h, s, v) {
    var r, g, b, i, f, p, q, t;
    if (arguments.length === 1) {
        s = h.s, v = h.v, h = h.h;
    }
    i = Math.floor(h * 6);
    f = h * 6 - i;
    p = v * (1 - s);
    q = v * (1 - f * s);
    t = v * (1 - (1 - f) * s);
    switch (i % 6) {
        case 0: r = v, g = t, b = p; break;
        case 1: r = q, g = v, b = p; break;
        case 2: r = p, g = v, b = t; break;
        case 3: r = p, g = q, b = v; break;
        case 4: r = t, g = p, b = v; break;
        case 5: r = v, g = p, b = q; break;
    }
    return {
        r: Math.round(r * 255),
        g: Math.round(g * 255),
        b: Math.round(b * 255)
    };
}

function rgb2hex(r,g,b) {
    r = Math.round(r).toString(16);
    g = Math.round(g).toString(16);
    b = Math.round(b).toString(16);

    if (r.length == 1)
        r = '0' + r;

    if (g.length == 1)
        g = '0' + g;

    if (b.length == 1)
        b = '0' + b;

    return '#' + r + g + b;
}

function getInteractionStrength(i) {
    switch (i) {
        case getFriendlyInteractionName('polarsidechainsidechaininteraction'):
        case getFriendlyInteractionName('polarbackbonesidechaininteraction'):
        case'polarsidechainsidechaininteraction':
        case'polarbackbonesidechaininteraction':
            return 4;
        case getFriendlyInteractionName('facetofaceinteraction'):
        case getFriendlyInteractionName('facetoedgeinteraction'):
        case getFriendlyInteractionName('picationinteraction'):
        case'facetofaceinteraction':
        case'facetoedgeinteraction':
        case'picationinteraction':
            return 3;
        case getFriendlyInteractionName('hydrophobicinteraction'):
        case'hydrophobicinteraction':
            return 2;
        case getFriendlyInteractionName('vanderwaalsinteraction'):
        case'vanderwaalsinteraction':
            return 1;
        default:
            return 0;
    }
}

function getColorStrongestInteraction(interactions, rgb = true) {
    var maxStrength = 0;
    for (var i = 0; i < interactions.length; i++)
        maxStrength = Math.max(maxStrength, getInteractionStrength(interactions[i].replace(/-/g, ' ')));

    return getInteractionColor(maxStrength, rgb);
}

function getFrequencyColor(frequency, rgb = true) {
    return getGradientColor(-1*frequency, rgb);
}

function getFlareGradientColor(fDiff, rgb){
    var color;
    var shift = 80;
    var basal = 255 - shift;

    if (fDiff <= 0)
        // If fDiff is close to -1, we want a red color
        color = { r: basal + (fDiff * -1*shift), g: basal-basal*(-fDiff), b: basal-basal*(-fDiff) };
    else
        // If fDiff is close to 1 we want a blue color
        color = { r: basal-basal*fDiff, g: basal-basal*fDiff, b: basal + (fDiff * shift)};

    if (rgb)
        return color;
    else
        return rgb2hex(color.r, color.g, color.b);
}

function getGradientColor(fDiff, rgb){
    var color;
    if (fDiff <= 0)
        // If fDiff is close to -1, we want a red color
        color = { r: 255, g: 255-255*(-fDiff), b: 255-255*(-fDiff) };
    else
        // If fDiff is close to 1 we want a blue color
        color = { r: 255-255*fDiff, g: 255-255*fDiff, b: 255 };

    if (rgb)
        return color;
    else
        return rgb2hex(color.r, color.g, color.b);
}

function getStrongestInteractionType(interactions) {
    if ($.inArray('polarsidechainsidechaininteraction', interactions) > -1)
        return 'polarsidechainsidechaininteraction';

    if ($.inArray('polarbackbonesidechaininteraction', interactions) > -1)
        return 'polarbackbonesidechaininteraction';

    if ($.inArray('facetofaceinteraction', interactions) > -1)
        return 'facetofaceinteraction';

    if ($.inArray('facetoedgeinteraction', interactions) > -1)
        return 'facetoedgeinteraction';

    if ($.inArray('picationinteraction', interactions) > -1)
        return 'picationinteraction';

    if ($.inArray('hydrophobicinteraction', interactions) > -1)
        return 'hydrophobicinteraction';

    if ($.inArray('vanderwaalsinteraction', interactions) > -1)
        return 'vanderwaalsinteraction';

    return 'undefined';
}

function getStrongestInteractionTypeFromPdbObject(obj) {

    var interactions = [];

    for (var key in obj) {
        if (Object.prototype.hasOwnProperty.call(obj, key)) {
            var strongestInteraction = getStrongestInteractionType(obj[key]);
            interactions.push(strongestInteraction);
        }
    }

    return getStrongestInteractionType(interactions);
}

function getInteractionTypesFromPdbObject(obj) {

    var interactions = new Set();
    for (var key in obj) {
        Object.keys(obj[key]).forEach(function(k,index) {
                interactions.add(obj[key][k]);
        });
        // if (Object.prototype.hasOwnProperty.call(obj, key)) {
        //     for (var k in obj[key])
        //         interactions.add(obj[key][k]);
        // }
    }

    // Sort according to strength
    interactions = Array.from(interactions);
    interactions.sort(function (i1, i2) {
        return  getInteractionStrength(i1) - getInteractionStrength(i2);
    });

    return interactions;
}


function getInteractionColor(interaction, rgb = true) {
    var r, g, b;

    switch (interaction) {
        case 'polarsidechainsidechaininteraction':
        case 'polarbackbonesidechaininteraction':
        case getFriendlyInteractionName('polarsidechainsidechaininteraction'):
        case getFriendlyInteractionName('polarbackbonesidechaininteraction'):
        case 4:
            //r = 254; g = 0; b = 16;
            r = 255; g = 98; b = 108;
            break;
        case 'facetofaceinteraction':
        case 'facetoedgeinteraction':
        case 'picationinteraction':
        case getFriendlyInteractionName('facetofaceinteraction'):
        case getFriendlyInteractionName('facetoedgeinteraction'):
        case getFriendlyInteractionName('picationinteraction'):
        case 3:
            //r = 94; g = 241; b = 242;
            r = 255; g = 166; b = 98;
            break;
        case 'hydrophobicinteraction':
        case getFriendlyInteractionName('hydrophobicinteraction'):
        case 2:
            //r = 0; g = 117; b = 220;
            r = 5; g = 200; b = 90;
            break;
        case 'vanderwaalsinteraction':
        case getFriendlyInteractionName('vanderwaalsinteraction'):
        case 1:
            //r = 89; g = 252; b = 197;
            r = 100; g = 100; b = 100;
            break;
        default:
            r = 0; g = 0; b = 0;
    }

    if (rgb)
        return { r: r, g: g, b: b };
    else
        return rgb2hex(r, g, b);
}

function getFriendlyInteractionName(interaction) {
    switch (interaction) {
        case 'polarsidechainsidechaininteraction':
        case 'polarbackbonesidechaininteraction':
            return 'Polar';
        case 'facetofaceinteraction':
        case 'facetoedgeinteraction':
        case 'picationinteraction':
            return 'Aromatic';
        case 'hydrophobicinteraction':
            return 'Hydrophobic';
        case 'vanderwaalsinteraction':
            return 'Van der Waals';
        default:
            return 'Unknown';
    }
}

/*function getFriendlyInteractionName(interaction) {
    switch (interaction) {
        case 'polarsidechainsidechaininteraction':
        return 'Polar (SC-SC)';
        case 'polarbackbonesidechaininteraction':
            return 'Polar (BB-SC)';
        case 'facetofaceinteraction':
            return 'Aromatic (F-F)';
        case 'facetoedgeinteraction':
            return 'Aromatic (F-E)';
        case 'picationinteraction':
            return 'Cation - pi';
        case 'hydrophobicinteraction':
            return 'Hydrophobic';
        case 'vanderwaalsinteraction':
            return 'Van der Waals';
        default:
            return 'Unknown';
    }
}*/

function getSegmentColor(segmentName) {
    var r, g, b;

    switch (segmentName) {
        case 'N-term':
        case 'C-term':
            r = 190; g = 190; b = 190;
            //r = 255; g = 150; b = 150;
            break;
        case 'TM1':
        case 'TM2':
        case 'TM3':
        case 'TM4':
        case 'TM5':
        case 'TM6':
        case 'TM7':
        case 'H8':
            r = 230; g = 230; b = 230;
            //r = 150; g = 255; b = 150;
            break;
        case 'ECL1':
        case 'ECL2':
        case 'ECL3':
            r = 190; g = 190; b = 190;
            //r = 150; g = 150; b = 255;
            break;
        case 'ICL1':
        case 'ICL2':
        case 'ICL3':
            r = 190; g = 190; b = 190;
            //r = 150; g = 150; b = 255;
            break;
        default:
            r = 0; g = 0; b = 0;
    }

    return { r: r, g: g, b: b };
}

function getAminoAcidOneLetterCode(threeLetterCode) {
    switch (threeLetterCode.toUpperCase()) {
        case 'ALA':
            return 'A';
        case 'ARG':
            return 'R';
        case 'ASN':
            return 'N';
        case 'ASP':
            return 'D';
        case 'CYS':
            return 'C';
        case 'GLN':
            return 'Q';
        case 'GLU':
            return 'E';
        case 'GLY':
            return 'G';
        case 'HIS':
            return 'H';
        case 'ILE':
            return 'I';
        case 'LEU':
            return 'L';
        case 'LYS':
            return 'K';
        case 'MET':
            return 'M';
        case 'PHE':
            return 'F';
        case 'PRO':
            return 'P';
        case 'SER':
            return 'S';
        case 'THR':
            return 'T';
        case 'TRP':
            return 'W';
        case 'TYR':
            return 'Y';
        case 'VAL':
            return 'V';
        default:
            return null;
    }
}

function downloadSVG(svgSelector, name) {
    var svgClone = $(svgSelector).clone();
    svgClone.find('.svg-pan-zoom_viewport').attr('transform', 'matrix(2.2,0,0,2.2,295,140)');

    var escapedSVG = new XMLSerializer().serializeToString(svgClone.get(0));

    downloadURI('data:image/svg+xml;base64,' + window.btoa(escapedSVG), name);
}

function downloadSingleCrystalCSV(singleCrystalSvgSelector, name) {
    var data = [];
    var header = ['Residue number 1', 'Residue number 2', 'Segment 1', 'Segment 2',  'Generic number 1', 'Generic number 2', 'Amino acid 1', 'Amino acid 2', 'Interaction type'];
    data.push(header);

    $(singleCrystalSvgSelector + ' rect[data-interaction-type]').each(function(e) {
        var rect = $(this);
        var resNo1 = rect.data('res-no-1');
        var resNo2 = rect.data('res-no-2');
        var seg1 = rect.data('seg-1');
        var seg2 = rect.data('seg-2');
        var genNo1 = rect.data('gen-no-1');
        var genNo2 = rect.data('gen-no-2');
        var aa1 = rect.data('aa-1');
        var aa2 = rect.data('aa-2');
        var iType = rect.data('interaction-type');
        data.push([resNo1, resNo2, seg1, seg2, genNo1, genNo2, aa1, aa2, iType]);
    });

    // Convert to CSV
    var csv = Papa.unparse(data);

    // Download file
    downloadURI('data:text/csv;charset=UTF-8,' + encodeURI(csv), name);
}

function populateTable(selector) {
    var table = [];
        rows = []
    if (selector=='single-table') {
        var header = ['Residue number 1', 'Residue number 2', 'Segment 1', 'Segment 2',  'Generic number 1', 'Generic number 2', 'Amino acid 1', 'Amino acid 2', 'Interaction type'];
        $('#single-crystal-tab .heatmap-container rect[data-interaction-type]').each(function(e) {
            var rect = $(this);
            var resNo1 = rect.data('res-no-1');
            var resNo2 = rect.data('res-no-2');
            var seg1 = rect.data('seg-1');
            var seg2 = rect.data('seg-2');
            var genNo1 = rect.data('gen-no-1');
            var genNo2 = rect.data('gen-no-2');
            var aa1 = rect.data('aa-1');
            var aa2 = rect.data('aa-2');
            var iType = rect.data('interaction-type');
            rows.push('<tr><td>'+resNo1+'</td><td>'+resNo2+'</td><td>'+seg1+'</td><td>'+seg2+'</td><td>'+genNo1+'</td><td>'+genNo2+'</td><td>'+aa1+'</td><td>'+aa2+'</td><td>'+iType+'</td>')
        });
    } else if (selector=='single-group-table') {
    var header = ['Generic number 1', 'Generic number 2', 'Segment 1', 'Segment 2',  'Number of interactions', 'Number of crystals', 'Frequency'];
    $('#single-crystal-group-tab .heatmap-container rect[data-frequency]').each(function(e) {
        var rect = $(this);
        var genNo1 = rect.data('gen-no-1');
        var genNo2 = rect.data('gen-no-2');
        var seg1 = rect.data('seg-1');
        var seg2 = rect.data('seg-2');
        var nInteractions = rect.data('num-interactions');
        var nTotalInteractions = rect.data('total-possible-interactions');
        var frequency = rect.data('frequency');
    //data.push([resNo1, resNo2, seg1, seg2, genNo1, genNo2, aa1, aa2, iType]);
    rows.push('<tr><td>'+genNo1+'</td><td>'+genNo2+'</td><td>'+seg1+'</td><td>'+seg2+'</td><td>'+nInteractions+'</td><td>'+nTotalInteractions+'</td><td>'+frequency+'</td>')
    });

    } else if (selector=='two-groups-table') {
    var header = ['Generic number 1', 'Generic number 2', 'Segment 1', 'Segment 2', 'Interactions group 1', 'Interactions group 2', 'Frequency group 1', 'Frequency group 2', 'Frequency Difference'];
    $('#two-crystal-groups-tab .heatmap-container rect[data-frequency-diff]').each(function(e) {
        var rect = $(this);
        var genNo1 = rect.data('gen-no-1');
        var genNo2 = rect.data('gen-no-2');
        var seg1 = rect.data('seg-1');
        var seg2 = rect.data('seg-2');
        var numIntsGroup1 = rect.data('group-1-num-ints');
        var numIntsGroup2 = rect.data('group-2-num-ints');
        var numPdbsGroup1 = rect.data('group-1-num-pdbs');
        var numPdbsGroup2 = rect.data('group-2-num-pdbs');
        var freqGroup1 = rect.data('group-1-freq');
        var freqGroup2 = rect.data('group-2-freq');
        var fDiff = rect.data('frequency-diff').toFixed(2);
        rows.push('<tr><td>'+genNo1+'</td><td>'+genNo2+'</td><td>'+seg1+'</td><td>'+seg2+'</td><td>'+numIntsGroup1+'</td><td>'+numIntsGroup2+'</td><td>'+freqGroup1+'</td><td>'+freqGroup2+'</td><td>'+fDiff+'</td>')
    });

    }
    table.push('<table id="'+selector+'" class="table display" width="100%"><thead><tr>');
    header.forEach(function(h) {
    table.push('<th></th>'); //'+h+'
    });
    table.push('</tr</thead>');
    table.push('<tbody>');
    table = table.concat(rows);
    table.push('</tbody></table>');
    $("#"+selector+"-tab").html(table.join( "" ));

    $("#"+selector+"-tab-link").click(function(e){
        setTimeout(renderTable,100,"#"+selector);
    });

}


function renderTable(mode) {
    if ( ! $.fn.DataTable.isDataTable(mode ) ) {
    oTable[mode] = $(mode).DataTable({
        'scrollX': true,
        scrollY:        '80vh',
        paging:         false,
        columnDefs: [
            { targets: 'no-sort', orderable: false }
        ],
        "aaSorting": [],
        });
        if (mode=='#single-table') {
            yadcf.init(oTable[mode],
            [
                {
                    column_number : 0,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Residue #1",
                    filter_reset_button_text: false,
                },
                {
                    column_number : 1,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Residue #2",
                    filter_reset_button_text: false,
                },
                {
                    column_number : 2,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Segment #1",
                    filter_reset_button_text: false,
                },
                {
                    column_number : 3,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Segment #2",
                    filter_reset_button_text: false,
                },
                {
                    column_number : 4,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "GN #1",
                    filter_reset_button_text: false,
                },
                {
                    column_number : 5,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "GN #2",
                    filter_reset_button_text: false,

                },
                {
                    column_number : 6,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Amino acid #1",
                    filter_reset_button_text: false,

                },
                {
                    column_number : 7,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Amino acid #2",
                    filter_reset_button_text: false,
                },
                {
                    column_number : 8,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Interaction Type",
                    filter_reset_button_text: false,

                },
            ],
            {
                cumulative_filtering: false
            }

        );
        } else if (mode=='#single-group-table') {
            yadcf.init(oTable[mode],
            [
                {
                    column_number : 0,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "GN #1",
                    filter_reset_button_text: false,
                },
                {
                    column_number : 1,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "GN #2",
                    filter_reset_button_text: false,
                },
                {
                    column_number : 2,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Segment #1",
                    filter_reset_button_text: false,
                },
                {
                    column_number : 3,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Segment #2",
                    filter_reset_button_text: false,
                },
                {
                    column_number : 4,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Number of interactions",
                    filter_reset_button_text: false,
                },
                {
                    column_number : 5,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Total PDBs",
                    filter_reset_button_text: false,

                },
                {
                    column_number : 6,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Frequency",
                    filter_reset_button_text: false,

                },
            ],
            {
                cumulative_filtering: false
            }

        );
        } else if (mode=='#two-groups-table') {
            yadcf.init(oTable[mode],
            [
                {
                    column_number : 0,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "GN #1",
                    filter_reset_button_text: false,
                },
                {
                    column_number : 1,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "GN #2",
                    filter_reset_button_text: false,
                },
                {
                    column_number : 2,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Segment #1",
                    filter_reset_button_text: false,
                },
                {
                    column_number : 3,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Segment #2",
                    filter_reset_button_text: false,
                },
                {
                    column_number : 4,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Interactions group 1",
                    filter_reset_button_text: false,
                },
                {
                    column_number : 5,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Interactions group 2",
                    filter_reset_button_text: false,

                },
                {
                    column_number : 6,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Frequency group 1",
                    filter_reset_button_text: false,

                },
                {
                    column_number : 7,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Frequency group 2",
                    filter_reset_button_text: false,

                },
                {
                    column_number : 8,
                    filter_type: "multi_select",
                    select_type: 'select2',
                    filter_default_label: "Frequency difference",
                    filter_reset_button_text: false,

                },
            ],
            {
                cumulative_filtering: false
            }

        );
        }

    yadcf.exResetAllFilters(oTable[mode]);
    }
}

function downloadSingleCrystalGroupCSV(singleGroupSvgSelector, name) {
    var data = [];
    var header = ['Generic number 1', 'Generic number 2', 'Segment 1', 'Segment 2', 'Frequency',  'Number of interactions', 'Number of crystals'];
    data.push(header);

    $(singleGroupSvgSelector + ' rect[data-frequency]').each(function(e) {
        var rect = $(this);
        var genNo1 = rect.data('gen-no-1');
        var genNo2 = rect.data('gen-no-2');
        var seg1 = rect.data('seg-1');
        var seg2 = rect.data('seg-2');
        var nInteractions = rect.data('num-interactions');
        var nTotalInteractions = rect.data('total-possible-interactions');
        var frequency = rect.data('frequency');
        data.push([genNo1, genNo2, seg1, seg2, nInteractions, nTotalInteractions, frequency]);
    });

    // Convert to CSV
    var csv = Papa.unparse(data);

    // Download file
    downloadURI('data:text/csv;charset=UTF-8,' + encodeURI(csv), name);
}

function downloadTwoCrystalGroupsCSV(twoGroupsSvgSelector, name) {
    var data = [];
    var header = ['Generic number 1', 'Generic number 2', 'Segment 1', 'Segment 2', 'Interactions group 1', 'Interactions group 2', 'Crystals group 1', 'Crystals group 2', 'Frequency group 1', 'Frequency group 2', 'Frequency Difference'];
    data.push(header);

    $(twoGroupsSvgSelector + ' rect[data-frequency-diff]').each(function(e) {
        var rect = $(this);
        var genNo1 = rect.data('gen-no-1');
        var genNo2 = rect.data('gen-no-2');
        var seg1 = rect.data('seg-1');
        var seg2 = rect.data('seg-2');
        var numIntsGroup1 = rect.data('group-1-num-ints');
        var numIntsGroup2 = rect.data('group-2-num-ints');
        var numPdbsGroup1 = rect.data('group-1-num-pdbs');
        var numPdbsGroup2 = rect.data('group-2-num-pdbs');
        var freqGroup1 = rect.data('group-1-freq');
        var freqGroup2 = rect.data('group-2-freq');
        var fDiff = rect.data('frequency-diff').toFixed(2);
        data.push([genNo1, genNo2, seg1, seg2, numIntsGroup1, numIntsGroup2, numPdbsGroup1, numPdbsGroup2, freqGroup1, freqGroup2, fDiff]);
    });

    // Convert to CSV
    var csv = Papa.unparse(data);

    // Download file
    downloadURI('data:text/csv;charset=UTF-8,' + encodeURI(csv), name);
}

function downloadURI(uri, name) {
    var link = document.createElement("a");
    link.download = name;
    link.href = uri;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    delete link;
}

function renderSingleCrystalHeatmap(data, heatMapSelector) {

    // Destroy old zoom on heatmap
    if (window.zoomHeatmap[heatMapSelector] != null) {
        window.zoomHeatmap[heatMapSelector].destroy();
        delete window.zoomHeatmap[heatMapSelector];
    }

    // Destroy old legend content
    $(heatMapSelector + ' .heatmap-legend').empty();

    // Destroy all previous contents
    $(heatMapSelector + ' .heatmap').empty();
    var heatmap = Snap(heatMapSelector + ' .heatmap');

    // Draw heatmap
    var interactions = data.interactions;
    var segment_map = data.segment_map;
    var sequence_numbers = data.sequence_numbers;
    var aa_map = data.aa_map[Object.keys(data.aa_map)[0]];
    var gen_map = data.generic_map;
    var num_seq_numbers = Object.keys(data.sequence_numbers).length;

    x = 0; wi = num_seq_numbers;
    y = 0; hi = num_seq_numbers;

    heatmap.attr({viewBox:[x,y,wi,hi].join(',')});
    heatmap.attr({viewBox:[x,y,wi,hi].join(' ')});

    // Contains all labels
    var labelGroup = heatmap.g();

    // Contains all labels
    var contentGroup = heatmap.g();

    // Compute segment offsets
    var i;

    var segments = [];

    var seg, prevSeg = segment_map[sequence_numbers[0]];
    var seqStart = 0;

    for (i = 0; i < num_seq_numbers; i++) {
        seg = segment_map[sequence_numbers[i]];

        if (seg === prevSeg) {
            continue;
        }

        segments.push({
            seg: prevSeg,
            start: seqStart,
            end: i-1
        });

        seqStart = i;
        prevSeg = seg;
    }

    // Push last segment
    segments.push({
        seg: prevSeg,
        start: seqStart,
        end: i-1
    });

    // Draw segments
    segments.forEach(function(s) {
        var rgb = getSegmentColor(s.seg);

        // Place the segments vertically.
        var cell = heatmap.rect(s.start, 0, s.end - s.start + 1, num_seq_numbers);
        var line = heatmap.line(s.start, 0, s.start, num_seq_numbers);

        cell.attr({
            'fill': "rgb(" + [rgb.r, rgb.g, rgb.b].join(',') + ")",
            'fill-opacity': "0.5"
        });

        line.attr({
            'stroke': "rgb(150,150,150)",
            'strokeWidth': "0.1"
        });

        // Add text
        var label = heatmap.text(s.start + (s.end - s.start + 1)/2 - 6, 9 + cell.getBBox().height, s.seg);

        label.attr({
            'text-anchor': 'start',
            'font-size': 5
        });

        contentGroup.add(cell);
        contentGroup.add(line);
        labelGroup.add(label);

        label.transform("r270s-1,1");
    });

    segments.forEach(function(s) {
        rgb = getSegmentColor(s.seg);

        // Place the segments horizontally.
        cell = heatmap.rect(0, s.start, num_seq_numbers, s.end - s.start + 1);
        line = heatmap.line(0, s.end+1, num_seq_numbers, s.end+1);

        cell.attr({
            'fill': "rgb(" + [rgb.r, rgb.g, rgb.b].join(',') + ")",
            'fill-opacity': "0.5"
        });

        line.attr({
            'stroke': "rgb(150,150,150)",
            'strokeWidth': "0.1"
        });

        var label = heatmap.text(-2, s.start + (s.end - s.start + 1)/2, s.seg);

        label.attr({
            'text-anchor': 'end',
            'font-size': 5,
            'alignment-baseline': 'middle'
        });

        label.transform("r180s-1,1");

        contentGroup.add(cell);
        contentGroup.add(line);
        labelGroup.add(label);
    });


    // Draw cells
    for (i = 0; i < num_seq_numbers; i++) {
        for (var j = 0; j < num_seq_numbers; j++) {
            // Get the sequence numbers
            var seq_i = data.sequence_numbers[i];
            var seq_j = data.sequence_numbers[j];

            // Only draw if an interaction exists
            var num = seq_i + "," + seq_j;

            if (num in interactions) {
                var cells = [];
                // Get all interactions at a given coordinate
                getInteractionTypesFromPdbObject(interactions[num]).forEach(function(interaction) {
                    var rgb = getInteractionColor(interaction);
                    var cell = heatmap.rect(i, j, 1, 1);

                    var interactionsString = getInteractionTypesFromPdbObject(interactions[num]).map(getFriendlyInteractionName).filter(function(item, pos, self) {
                        return self.indexOf(item) == pos;
                    }).join(", ");

                    var content, title = 'Residues ' + aa_map[seq_i] + seq_i + '-' + aa_map[seq_j] + seq_j + '<br />'
                            + 'Interactions: ' + interactionsString + '<br />'
                            + 'Segments: ' + segment_map[seq_i] + ', ' + segment_map[seq_j] + '<br />';

                    // Add generic numbers where applicable
                    if (seq_i in gen_map) {
                        title += 'Res. 1 gen. no: ' + gen_map[seq_i] + '<br />';
                    }

                    if (seq_j in gen_map) {
                        title += 'Res. 2 gen. no: ' + gen_map[seq_j] + '<br />';
                    }

                    var genStrI = gen_map[seq_i];
                    var genStrJ = gen_map[seq_j];

                    if (typeof gen_map[seq_i] == 'undefined') {
                        genStrI = '-';
                    }

                    if (typeof gen_map[seq_j] == 'undefined') {
                        genStrJ = '-';
                    }

                    var popoverTable = '<table class="table">'
                    + '<thead>'
                    + '<tr>'
                    + '<th>Residue</th>'
                    + '<th>' + aa_map[seq_i] + seq_i + '</th>'
                    + '<th>' + aa_map[seq_j] + seq_j + '</th>'
                    + '</tr>'
                    + '</thead>'
                    + '<tbody>'
                    + '<td>Segment</td>'
                    + '<td>' + segment_map[seq_i] + '</td>'
                    + '<td>' + segment_map[seq_j] + '</td>'
                    + '</tr>'
                    + '<tr>'
                    + '<td>Gen. no.</td>'
                    + '<td>' + genStrI + '</td>'
                    + '<td>' + genStrJ + '</td>'
                    + '</tr>'
                    + '</tbody>'
                    + '</table>'
                    + '<table class="table">'
                    + '<thead>'
                    + '<tr>'
                    + '<th>Interactions</th>'
                    + '</tr>'
                    + '<tr>'
                    + '</thead>'
                    + '<tbody>'
                    + '<tr>'
                    + '<td>' + interactionsString + '</td>'
                    + '</tr>'
                    + '</tbody>'
                    + '</table>'

                    var interactionId = 'single-i-' + num.replace(',','-');

                    // Create cell for interaction
                    cell.attr({
                        'fill': "rgb(" + [rgb.r, rgb.g, rgb.b].join(',') + ")",
                        'data-interaction-type': interaction,
                        'data-res-no-1': seq_i,
                        'data-res-no-2': seq_j,
                        'data-aa-1': aa_map[seq_i],
                        'data-aa-2': aa_map[seq_j],
                        'data-seg-1': segment_map[seq_i],
                        'data-seg-2': segment_map[seq_j],
                        'data-gen-no-1': genStrI,
                        'data-gen-no-2': genStrJ,
                        'class': 'heatmap-interaction' + ' ' + getFriendlyInteractionName(interaction).replace(/ /g,"-"),
                        'id': interactionId
                    });

                    contentGroup.add(cell);

                    cell = $(heatMapSelector + ' rect.heatmap-interaction#' + interactionId);

                    // Add tooltip to cell
                    cell.tooltip({
                        'container': heatMapSelector,
                        'placement': 'top',
                        'delay': 0,
                        'html': true,
                        'title': title
                    });

                    // Add popover to cells
                    cell.popover({
                        'container': heatMapSelector,
                        'placement': 'bottom',
                        'animation': true,
                        'html': true,
                        'title': 'Interactions at ' + aa_map[seq_i] + seq_i + ', ' + aa_map[seq_j] + seq_j,
                        'content': popoverTable,
                        'tabindex': '0'
                    });

                    cells.push(cell);
                });
            }
        }
    }

    // Add cover-up triangle
    var bbox = contentGroup.getBBox();
    contentGroup.add(heatmap.polygon([bbox.x, bbox.y, bbox.x2, bbox.y, bbox.x2, bbox.y2]).attr({ fill: "white" }));

    // Add bottom line
    var line = heatmap.line(bbox.x, bbox.y, bbox.x2, bbox.y2);
    line.attr({
        'stroke': "rgb(150,150,150)",
        'strokeWidth': "0.1"
    });
    contentGroup.add(line);

    // Rotate contents
    var g = heatmap.g();
    g.add(contentGroup);
    g.add(labelGroup);
    g.transform("r225s-1,1");


    // Populate heatmap legend
    var interactionTypes = new Set();

    $(heatMapSelector + ' .heatmap-interaction').each(function() {
        var friendlyName = getFriendlyInteractionName($(this).data('interaction-type'));
        interactionTypes.add(friendlyName);
    });

    // Add interactions color legend
    var legendHtml = '<ul>';

    interactionTypes = Array.from(interactionTypes).sort(function (i1, i2) {
        return getInteractionStrength(i2) - getInteractionStrength(i1);
    });

    interactionTypes.forEach(function(i) {
        var rgb = getInteractionColor(i);
        legendHtml = legendHtml
                + '<li>'
                + '<div class="color-box" style="background-color: ' + rgb2hex(rgb.r, rgb.g, rgb.b) + '">' + '<input type="checkbox" data-interaction-type="' + i.replace(/ /g,"-") +'"></input>' + '</div><p>' + i + '</p>'
                + '</li>';
    });
    legendHtml += '</ul>';

    // Add SVG download button
    legendHtml += '<button onclick="downloadSVG(\'' + heatMapSelector + ' .heatmap\', \'interactions.svg\')" type="button" class="btn btn-primary pull-right svg-download-button" aria-label="Left Align">' +
                    '<span class="glyphicon glyphicon-download" aria-hidden="true"></span> Download SVG' +
                    '</button>';

    // Add CSV download button
    legendHtml += '<br /><button onclick="downloadSingleCrystalCSV(\'' + heatMapSelector + ' .heatmap\', \'interactions.csv\')" type="button" class="btn btn-success pull-right csv-download-button" aria-label="Left Align"><span class="glyphicon glyphicon-download" aria-hidden="true"></span> Download CSV' +
        '</button>';

    $(heatMapSelector + ' .heatmap-legend').append(legendHtml);

    $(heatMapSelector + ' .heatmap-legend input[type=checkbox]').each(function() {
        $(this).prop('checked', true);
        $(this).change(function() {
            var interactionType = $(this).data('interaction-type');
            var rects = $(heatMapSelector + ' .heatmap rect.' + interactionType);
            if ($(this).is(':checked')) {
                rects.show();
            } else {
                rects.hide();
            }
        });
    });


    // Make zoomable
    window.zoomHeatmap[heatMapSelector] = svgPanZoom(heatMapSelector + ' .heatmap', {
        zoomEnabled: true,
        // controlIconsEnabled: true,
        fit: true,
        center: true,
        minZoom: 0.75,
        maxZoom: 50,
        zoomScaleSensitivity: 0.25,
        dblClickZoomEnabled: true,
        beforeZoom: hidePopovers,
        beforePan: hidePopovers
    });

    // Set initial zoom level
    window.zoomHeatmap[heatMapSelector].zoom(1);
    window.zoomHeatmap[heatMapSelector].pan({x: 325, y: 140});

    // Close popovers on clicking elsewhere
    $('html').on('mousedown', function(e) {
        if(!$(e.target).closest('.popover').length) {
            if ($(e.target).closest(heatMapSelector).length) {
                hidePopovers();
            }
        }
    });

    // Create table view
    populateTable('single-table');
}

function renderSingleCrystalGroupHeatmap(data, heatMapSelector) {

    // Destroy old zoom on heatmap
    if (window.zoomHeatmap[heatMapSelector] != null) {
        window.zoomHeatmap[heatMapSelector].destroy();
        delete window.zoomHeatmap[heatMapSelector];
    }

    // Destroy old legend content
    $(heatMapSelector + ' .heatmap-legend').empty();

    // Destroy all previous contents
    $(heatMapSelector + ' .heatmap').empty();
    var heatmap = Snap(heatMapSelector + ' .heatmap');

    // Draw heatmap
    var interactions = data.interactions;
    var segment_map = data.segment_map;
    var sequence_numbers = data.sequence_numbers;
    var num_seq_numbers = Object.keys(data.sequence_numbers).length;

    x = 0; wi = num_seq_numbers;
    y = 0; hi = num_seq_numbers;

    heatmap.attr({viewBox:[x,y,wi,hi].join(',')});
    heatmap.attr({viewBox:[x,y,wi,hi].join(' ')});

    // Contains all labels
    var labelGroup = heatmap.g();

    // Contains all content
    var contentGroup = heatmap.g();

    // Compute segment offsets
    var i;

    var segments = [];

    var seg, prevSeg = segment_map[sequence_numbers[0]];
    var seqStart = 0;

    for (i = 0; i < num_seq_numbers; i++) {
        seg = segment_map[sequence_numbers[i]];

        if (seg === prevSeg) {
            continue;
        }

        segments.push({
            seg: prevSeg,
            start: seqStart,
            end: i-1
        });

        seqStart = i;
        prevSeg = seg;
    }

    // Push last segment
    segments.push({
        seg: prevSeg,
        start: seqStart,
        end: i-1
    });

    // Draw segments
    segments.forEach(function(s) {
        var rgb = getSegmentColor(s.seg);

        // Place the segments vertically.
        var cell = heatmap.rect(s.start, 0, s.end - s.start + 1, num_seq_numbers);
        var line = heatmap.line(s.start, 0, s.start, num_seq_numbers);

        cell.attr({
            'fill': "rgb(" + [rgb.r, rgb.g, rgb.b].join(',') + ")",
            'fill-opacity': "0.5"
        });

        line.attr({
            'stroke': "rgb(150,150,150)",
            'strokeWidth': "0.1"
        });


        // Add text
        var label = heatmap.text(s.start + (s.end - s.start + 1)/2 - 6, 9 + cell.getBBox().height, s.seg);

        label.attr({
            'text-anchor': 'start',
            'font-size': 5
        });

        contentGroup.add(cell);
        contentGroup.add(line);
        labelGroup.add(label);

        label.transform("r270s-1,1");
    });

    segments.forEach(function(s) {
        rgb = getSegmentColor(s.seg);

        // Place the segments horizontally.
        cell = heatmap.rect(0, s.start, num_seq_numbers, s.end - s.start + 1);
        line = heatmap.line(0, s.end+1, num_seq_numbers, s.end+1);

        cell.attr({
            'fill': "rgb(" + [rgb.r, rgb.g, rgb.b].join(',') + ")",
            'fill-opacity': "0.5"
        });

        line.attr({
            'stroke': "rgb(150,150,150)",
            'strokeWidth': "0.1"
        });

        var label = heatmap.text(-2, s.start + (s.end - s.start + 1)/2, s.seg);

        label.attr({
            'text-anchor': 'end',
            'font-size': 5,
            'alignment-baseline': 'middle'
        });

        label.transform("r180s-1,1");

        contentGroup.add(cell);
        contentGroup.add(line);
        labelGroup.add(label);
    });

    // Draw cells
    for (i = 0; i < num_seq_numbers; i++) {
        for (var j = 0; j < num_seq_numbers; j++) {
            // Get the sequence numbers
            var seq_i = data.sequence_numbers[i];
            var seq_j = data.sequence_numbers[j];

            // Only draw if an interaction exists
            var num = seq_i + "," + seq_j;

            if (num in interactions) {
                // Get the strongest interaction.
                var nInteractions = Object.keys(interactions[num]).length;
                var frequency = nInteractions / data.pdbs.length;

                var rgb = { r: 255, g: 255-frequency*255, b: 255-frequency*255 };
                var cell = heatmap.rect(i, j, 1, 1);

                var title =  'Residues ' + seq_i + ', ' + seq_j + '<br />'
                        + 'Interaction count: ' + nInteractions + '<br />'
                        + segment_map[seq_i] + ', ' + segment_map[seq_j];

                var popoverTable = '<table class="table">'
                    + '<thead>'
                    + '<tr>'
                    + '<th>Residue #</th>'
                    + '<th>' + seq_i + '</th>'
                    + '<th>' + seq_j + '</th>'
                    + '</tr>'
                    + '</thead>'
                    + '<tbody>'
                    + '<td>Segment</td>'
                    + '<td>' + segment_map[seq_i] + '</td>'
                    + '<td>' + segment_map[seq_j] + '</td>'
                    + '</tr>'
                    + '</tbody>'
                    + '</table>'
                    + 'Interaction count: ' + nInteractions + '<br />'
                    + 'Interaction frequency: ' + frequency.toFixed(2)

                cell.attr({
                    'fill': "rgb(" + [rgb.r, rgb.g, rgb.b].join(',') + ")",
                    'data-num-interactions': nInteractions,
                    'data-total-possible-interactions': data.pdbs.length,
                    'data-frequency': frequency,
                    'data-gen-no-1': seq_i,
                    'data-gen-no-2': seq_j,
                    'data-seg-1': segment_map[seq_i],
                    'data-seg-2': segment_map[seq_j],
                    'class': 'heatmap-interaction'
                });

                $(heatMapSelector + ' rect.heatmap-interaction').tooltip({
                    'container': heatMapSelector,
                    'placement': 'top',
                    'delay': 75,
                    'title': title,
                    'html': true
                });

                // Add popover to cells
                $(heatMapSelector + ' rect.heatmap-interaction').popover({
                    'container': heatMapSelector,
                    'placement': 'bottom',
                    'animation': true,
                    'html': true,
                    'title': 'Interactions at ' + seq_i + ', ' + seq_j,
                    'content': popoverTable,
                    'tabindex': '0'
                });

                contentGroup.add(cell);
            }
        }
    }

    // Add cover-up triangle
    var bbox = contentGroup.getBBox();
    contentGroup.add(heatmap.polygon([bbox.x, bbox.y, bbox.x2, bbox.y, bbox.x2, bbox.y2]).attr({ fill: "white" }));

    // Rotate contents
    var g = heatmap.g();
    g.add(contentGroup);
    g.add(labelGroup);
    g.transform("r225s-1,1");

    // Populate heatmap legend
    var legendHtml = '<h4 class="center">Interaction count</h4>'
        + '<p>Range: <span id="pdbs-range">0 - ' + data.pdbs.length + '</span></p>'
        + '<div class="slider-range" data-text-id="pdbs-range" id="pdbs-range-slider"></div>'
        + '<div class="temperature-scale">'
        + '<span class="white-to-red"></span>'
        + '</div>';

    // Add SVG download button
    legendHtml += '<button onclick="downloadSVG(\'' + heatMapSelector + ' .heatmap\', \'interactions.svg\')" type="button" class="btn btn-primary pull-right svg-download-button" aria-label="Left Align">' +
                    '<span class="glyphicon glyphicon-download" aria-hidden="true"></span> Download SVG' +
                    '</button>';

    // Add CSV download button
    legendHtml += '<br /><button onclick="downloadSingleCrystalGroupCSV(\'' + heatMapSelector + ' .heatmap\', \'interactions.csv\')" type="button" class="btn btn-success pull-right csv-download-button" aria-label="Left Align"><span class="glyphicon glyphicon-download" aria-hidden="true"></span> Download CSV' +
        '</button>';


    $(heatMapSelector + ' .heatmap-legend').append(legendHtml);

    $( function() {
        $( heatMapSelector+" .slider-range" ).slider({
        range: true,
        min: 0,
        max: data.pdbs.length,
        step: 1,
        values: [0,data.pdbs.length],
        slide: function( event, ui ) {
            $( "#"+$(this).attr("data-text-id") ).html( ui.values[ 0 ] + " - " + ui.values[ 1 ] );
            setTimeout(getRangeChangeFunction,100)
        }
        });
    } );

    function getRangeChangeFunction() {

            var [tMin,tMax] = $(heatMapSelector + ' .heatmap-legend #pdbs-range-slider').slider( "option", "values" );

            // Hide all below min treshold
            $(heatMapSelector + ' rect').each(function() {
                var n = $(this).data("num-interactions");
                if (n < tMin || tMax < n) {
                    $(this).hide();
                } else {
                    $(this).show();
                }
            });

    }

    // Make zoomable
    window.zoomHeatmap[heatMapSelector] = svgPanZoom(heatMapSelector + ' .heatmap', {
        zoomEnabled: true,
        // controlIconsEnabled: true,
        fit: true,
        center: true,
        minZoom: 0.75,
        maxZoom: 50,
        zoomScaleSensitivity: 0.25,
        dblClickZoomEnabled: true,
        beforeZoom: hidePopovers,
        beforePan: hidePopovers
    });

    // Set initial zoom level
    window.zoomHeatmap[heatMapSelector].zoom(1);
    window.zoomHeatmap[heatMapSelector].pan({x: 325, y: 140});

    // Close popovers on clicking elsewhere
    $('html').on('mousedown', function(e) {
        if(!$(e.target).closest('.popover').length) {
            if ($(e.target).closest(heatMapSelector).length) {
                hidePopovers();
            }
        }
    });
    // Create table view
    populateTable('single-group-table');
}

function renderTwoCrystalGroupsHeatmap(data1, data2, data3, heatMapSelector) {

    // Destroy old zoom on heatmap
    if (window.zoomHeatmap[heatMapSelector] != null) {
        window.zoomHeatmap[heatMapSelector].destroy();
        delete window.zoomHeatmap[heatMapSelector];
    }

    // Destroy old legend content
    $(heatMapSelector + ' .heatmap-legend').empty();

    // Destroy all previous contents
    $(heatMapSelector + ' .heatmap').empty();
    var heatmap = Snap(heatMapSelector + ' .heatmap');

    // Draw heatmap
    var interactions = data3.interactions;
    var segment_map = data3.segment_map;
    var sequence_numbers = data3.sequence_numbers;
    var num_seq_numbers = Object.keys(data3.sequence_numbers).length;

    x = 0; wi = num_seq_numbers;
    y = 0; hi = num_seq_numbers;

    heatmap.attr({viewBox:[x,y,wi,hi].join(',')});
    heatmap.attr({viewBox:[x,y,wi,hi].join(' ')});

    // Contains all labels
    var labelGroup = heatmap.g();

    // Contains all labels
    var contentGroup = heatmap.g();

    // Compute segment offsets
    var i;

    var segments = [];

    var seg, prevSeg = segment_map[sequence_numbers[0]];
    var seqStart = 0;

    for (i = 0; i < num_seq_numbers; i++) {
        seg = segment_map[sequence_numbers[i]];

        if (seg === prevSeg) {
            continue;
        }

        segments.push({
            seg: prevSeg,
            start: seqStart,
            end: i-1
        });

        seqStart = i;
        prevSeg = seg;
    }

    // Push last segment
    segments.push({
        seg: prevSeg,
        start: seqStart,
        end: i-1
    });

    // Draw segments
    segments.forEach(function(s) {
        var rgb = getSegmentColor(s.seg);

        // Place the segments vertically.
        var cell = heatmap.rect(s.start, 0, s.end - s.start + 1, num_seq_numbers);
        var line = heatmap.line(s.start, 0, s.start, num_seq_numbers);

        cell.attr({
            'fill': "rgb(" + [rgb.r, rgb.g, rgb.b].join(',') + ")",
            'fill-opacity': "0.5"
        });

        line.attr({
            'stroke': "rgb(150,150,150)",
            'strokeWidth': "0.1"
        });

        // Add text
        var label = heatmap.text(s.start + (s.end - s.start + 1)/2 - 6, 9 + cell.getBBox().height, s.seg);

        label.attr({
            'text-anchor': 'start',
            'font-size': 5
        });

        contentGroup.add(cell);
        contentGroup.add(line);
        labelGroup.add(label);

        label.transform("r270s-1,1");
    });

    segments.forEach(function(s) {
        rgb = getSegmentColor(s.seg);

        // Place the segments horizontally.
        cell = heatmap.rect(0, s.start, num_seq_numbers, s.end - s.start + 1);
        line = heatmap.line(0, s.end+1, num_seq_numbers, s.end+1);

        cell.attr({
            'fill': "rgb(" + [rgb.r, rgb.g, rgb.b].join(',') + ")",
            'fill-opacity': "0.5"
        });

        line.attr({
            'stroke': "rgb(150,150,150)",
            'strokeWidth': "0.1"
        });

        var label = heatmap.text(-2, s.start + (s.end - s.start + 1)/2, s.seg);

        label.attr({
            'text-anchor': 'end',
            'font-size': 5,
            'alignment-baseline': 'middle'
        });

        label.transform("r180s-1,1");

        contentGroup.add(cell);
        contentGroup.add(line);
        labelGroup.add(label);
    });

    // Draw cells
    for (i = 0; i < num_seq_numbers; i++) {
        for (var j = 0; j < num_seq_numbers; j++) {
            // Get the sequence numbers
            var seq_i = data3.sequence_numbers[i];
            var seq_j = data3.sequence_numbers[j];

            // Only draw if an interaction exists
            var num = seq_i + "," + seq_j;

            var n1 = 0;
            var n2 = 0;

            if (num in data1.interactions) {
                n1 = Object.keys(data1.interactions[num]).length;
            }

            if (num in data2.interactions) {
                n2 = Object.keys(data2.interactions[num]).length;
            }

            if ((num in data1.interactions) || (num in data2.interactions)) {
                // Difference in frequencies
                var f1 = (n1 / data1.pdbs.length);
                var f2 = (n2 / data2.pdbs.length);
                var fDiff = (n1 / data1.pdbs.length) - (n2 / data2.pdbs.length);

                var rgb = getGradientColor(fDiff, true);

                var cell = heatmap.rect(i, j, 1, 1);

                var title =  'Residues ' + seq_i + ', ' + seq_j + '<br />'
                            + 'Frequency group 1: ' + f1.toFixed(2)+ '<br />'
                            + 'Frequency group 2: ' + f2.toFixed(2)+ '<br />'
                            + 'Frequency difference: ' + fDiff.toFixed(2)+ '<br />' ;

                var popoverTable = '<table class="table">'
                    + '<thead>'
                    + '<tr>'
                    + '<th>Residue #</th>'
                    + '<th>' + seq_i + '</th>'
                    + '<th>' + seq_j + '</th>'
                    + '</tr>'
                    + '</thead>'
                    + '<tbody>'
                    + '<td>Segment</td>'
                    + '<td>' + segment_map[seq_i] + '</td>'
                    + '<td>' + segment_map[seq_j] + '</td>'
                    + '</tr>'
                    + '</tbody>'
                    + '</table>'
                    + 'Group 1 freq: ' + f1.toFixed(2) + '<br />'
                    + 'Group 2 freq: ' + f2.toFixed(2) + '<br />'
                    + 'Frequency difference: ' + fDiff.toFixed(2)

                cell.attr({
                    'fill': "rgb(" + [rgb.r, rgb.g, rgb.b].join(',') + ")",
                    'data-frequency-diff': fDiff,
                    'data-gen-no-1': seq_i,
                    'data-gen-no-2': seq_j,
                    'data-seg-1': segment_map[seq_i],
                    'data-seg-2': segment_map[seq_j],
                    'data-group-1-num-ints': n1,
                    'data-group-2-num-ints': n2,
                    'data-group-1-num-pdbs': data1.pdbs.length,
                    'data-group-2-num-pdbs': data2.pdbs.length,
                    'data-group-1-freq': f1.toFixed(2),
                    'data-group-2-freq': f2.toFixed(2),
                    'class': 'heatmap-interaction'
                });

                $(heatMapSelector + ' rect.heatmap-interaction').tooltip({
                    'container': heatMapSelector,
                    'placement': 'top',
                    'delay': 75,
                    'html': true,
                    'title': title
                });

                // Add popover to cells
                $(heatMapSelector + ' rect.heatmap-interaction').popover({
                    'container': heatMapSelector,
                    'placement': 'bottom',
                    'animation': true,
                    'html': true,
                    'title': 'Interactions at ' + seq_i + ', ' + seq_j,
                    'content': popoverTable,
                    'tabindex': '0'
                });

                contentGroup.add(cell);
            }
        }
    }

    // Add cover-up triangle
    var bbox = contentGroup.getBBox();
    contentGroup.add(heatmap.polygon([bbox.x, bbox.y, bbox.x2, bbox.y, bbox.x2, bbox.y2]).attr({ fill: "white" }));

    // Rotate contents
    var g = heatmap.g();
    g.add(contentGroup);
    g.add(labelGroup);
    g.transform("r225s-1,1");

    // Populate heatmap legend
    var legendHtml = '<h4 class="center">Frequency</h4>'
        + '<p>Group 1 range: <span id="freq-range-1">0 - 1</span></p>'
        + '<div class="slider-range" data-text-id="freq-range-1" id="freq-slider-range-1"></div>'
        + '<p>Group 2 range: <span id="freq-range-2">0 - 1</span></p>'
        + '<div class="slider-range" data-text-id="freq-range-2" id="freq-slider-range-2"></div>'
        + '<p>Freq difference range: <span id="freq-range-3">-1 - 1</span></p>'
        + '<div class="slider-range-diff" data-text-id="freq-range-3" id="freq-slider-range-3"></div>'
        + '<div class="temperature-scale">'
        + '<span class="red-to-white"></span>'
        + '<span class="white-to-blue"></span>'
        + '</div>'
        + '</div>';

    // Add SVG download button
    legendHtml += '<button onclick="downloadSVG(\'' + heatMapSelector + ' .heatmap\', \'interactions.svg\')" type="button" class="btn btn-primary pull-right svg-download-button" aria-label="Left Align">' +
                    '<span class="glyphicon glyphicon-download" aria-hidden="true"></span> Download SVG' +
                    '</button>';

    // Add CSV download button
    legendHtml += '<br /><button onclick="downloadTwoCrystalGroupsCSV(\'' + heatMapSelector + ' .heatmap\', \'interactions.csv\')" type="button" class="btn btn-success pull-right csv-download-button" aria-label="Left Align"><span class="glyphicon glyphicon-download" aria-hidden="true"></span> Download CSV' +
        '</button>';

    $(heatMapSelector + ' .heatmap-legend').append(legendHtml);

    $( function() {
        $( ".slider-range" ).slider({
        range: true,
        min: 0,
        max: 1,
        step: 0.01,
        values: [0,1],
        slide: function( event, ui ) {
            $( "#"+$(this).attr("data-text-id") ).html( ui.values[ 0 ] + " - " + ui.values[ 1 ] );
            setTimeout(getRangeChangeFunction,100)
        }
        });
    } );
    $( function() {
        $( ".slider-range-diff" ).slider({
        range: true,
        min: -1,
        max: 1,
        step: 0.01,
        values: [-1,1],
        slide: function( event, ui ) {
            $( "#"+$(this).attr("data-text-id") ).html( ui.values[ 0 ] + " - " + ui.values[ 1 ] );
            setTimeout(getRangeChangeFunction,100)
        }
        });
    } );


    function getRangeChangeFunction() {

            var r1 = $(heatMapSelector + ' .heatmap-legend #freq-slider-range-1').slider( "option", "values" );
            var r2 = $(heatMapSelector + ' .heatmap-legend #freq-slider-range-2').slider( "option", "values" );
            var r3 = $(heatMapSelector + ' .heatmap-legend #freq-slider-range-3').slider( "option", "values" );

            // Hide all below min treshold
            $(heatMapSelector + ' rect').each(function() {
                var f1 = $(this).data("group-1-freq");
                var f2 = $(this).data("group-2-freq");
                var f3 = $(this).data("frequency-diff");
                if ( (f1 < r1[0] || r1[1] < f1) || (f2 < r2[0] || r2[1] < f2) || (f3 < r3[0] || r3[1] < f3) ) {
                    $(this).hide();
                } else {
                    $(this).show();
                }
            });
    }

    // Make zoomable
    window.zoomHeatmap[heatMapSelector] = svgPanZoom(heatMapSelector + ' .heatmap', {
        zoomEnabled: true,
        // controlIconsEnabled: true,
        fit: true,
        center: true,
        minZoom: 0.75,
        maxZoom: 50,
        zoomScaleSensitivity: 0.25,
        dblClickZoomEnabled: true,
        beforeZoom: hidePopovers,
        beforePan: hidePopovers
    });

    // Set initial zoom level
    window.zoomHeatmap[heatMapSelector].zoom(1);
    window.zoomHeatmap[heatMapSelector].pan({x: 325, y: 140});

    // Close popovers on clicking elsewhere
    $('html').on('mousedown', function(e) {
        if(!$(e.target).closest('.popover').length) {
            if ($(e.target).closest(heatMapSelector).length) {
                hidePopovers();
            }
        }
    });
    populateTable('two-groups-table');
}

function initializeSegmentButtons(selector) {
    // Initialize segment buttons.
    $(selector + ' .segments-panel button').each(function() {
        var s = $(this).attr('data-segment');

        // Return if no segment data
        if (s == null) {
            return;
        }

        $(this).click(function() {
            var segments = [];
            $(this).toggleClass('active');
            $(selector + ' .segments-panel button.active').each(function() {
                segments = segments.concat($(this).data('segment').split(' '));
            });
            $(selector + ' .segments-input').val(JSON.stringify(segments));
        });
    });

    // Initialize 'all' buttons.
    $(selector  + ' .segments-panel .all-button').each(function() {
        $(this).click(function() {
            if ($(this).html() === 'All') {
                $(this).html('None');
                $(this).parent().find('button').each(function() {
                    var s = $(this).attr('data-segment');

                    // Return if no segment data
                    if (s == null) {
                        return;
                    }

                    if (!$(this).hasClass('active')) {
                        $(this).trigger('click');
                    }
                });
            } else {
                $(this).html('All');
                $(this).parent().find('button').each(function() {
                    var s = $(this).attr('data-segment');

                    // Return if no segment data
                    if (s == null) {
                        return;
                    }

                    if ($(this).hasClass('active')) {
                        $(this).trigger('click');
                    }
                });
            }

        });

        // Update data
        var segments = [];
        $(selector + ' .segments-panel button.active').each(function() {
            segments.push($(this).data('segment'));
        });
        $(selector + ' .segments-panel .segments-input').val(JSON.stringify(segments));

        // Trigger click on initialization
        $(this).trigger('click');
    });
}

function initializeInteractionButtons(selector) {
    // Initialize interaction buttons.
    $(selector + ' .interactions-panel button').each(function() {
        var s = $(this).attr('data-interaction-type');

        // Return if no segment data
        if (s == null) {
            return;
        }

        $(this).click(function() {
            var interactions = [];
            $(this).toggleClass('active');
            $(selector + ' .interactions-panel button.active').each(function() {
                interactions = interactions.concat($(this).data('interaction-type').split(' '));
            });
            $(selector + ' .interactions-input').val(JSON.stringify(interactions));
        });
    });

    // Initialize 'all' buttons.
    $(selector  + ' .interactions-panel .all-button').each(function() {
        $(this).click(function() {
            if ($(this).html() === 'All') {
                $(this).html('None');
                $(this).parent().find('button').each(function() {
                    var s = $(this).attr('data-interaction-type');

                    // Return if no segment data
                    if (s == null) {
                        return;
                    }

                    if (!$(this).hasClass('active')) {
                        $(this).trigger('click');
                    }
                });
            } else {
                $(this).html('All');
                $(this).parent().find('button').each(function() {
                    var s = $(this).attr('data-interaction-type');

                    // Return if no segment data
                    if (s == null) {
                        return;
                    }

                    if ($(this).hasClass('active')) {
                        $(this).trigger('click');
                    }
                });
            }
        });

        // Update data
        var interactions = [];
        $(selector + ' .interactions-panel button.active').each(function() {
            interactions.push($(this).data('interaction-type'));
        });
        $(selector + ' .interactions-input').val(JSON.stringify(interactions));

        // Trigger click on initialization
        $(this).trigger('click');
    });
}

function initializeGoButton(selector, heatmapFunction, generic=false) {
    $(selector + ' .go-button').click(function() {
        var pdb = JSON.parse($(selector + ' .crystal-pdb').val());
        //var segments = JSON.parse($(selector + ' .segments-input').val());
        var segments = ['TM1','TM2','TM3','TM4','TM5','TM6','TM7','TM1','ICL1','ECL1','ICL2','ECL2','ICL3','ECL3','N-term','C-term'];
        if (pdb.length > 0 && segments.length > 0) {
            var interactionTypes = JSON.parse($(selector + ' .interactions-input').val());
            $(".heatmap").hide();
            $(".heatmap-legend").hide();
            $(".matrix-tab:visible").click();

            $(selector + ' .heatmap-container').append('<span id=svgloading>Loading...</span>');
            if (!$(selector + ' .interactions-input').val() == null)
                interactionTypes = JSON.parse($(selector + ' .interactions-input').val());

            $.getJSON( '/contactnetwork/interactiondata',
            {
                'segments': segments,
                'generic': generic,
                'pdbs': pdb,
                'interaction_types': interactionTypes
            },
            function( data ) {
                // Re-render heatmap
                $(".heatmap").show();
                $(".heatmap-legend").show();
                heatmapFunction(data, selector + ' .heatmap-container');
                $("#svgloading").remove()
                // Re-render flareplot
                createFlareplotBox(data, selector + " .flareplot-container");

                if (selector == '#single-crystal-group-tab') {
                    createSchematicPlot(data, selector + " .schematic-con-container", {type: 'singleCrystalGroup'}); //.schematic-container
                    createSchematicPlot(data, selector + " .schematic-non-container", {isContiguousPlot: false, type: 'singleCrystalGroup'});
                    createNGLview("single-group",pdb[0], pdb);
                } else {
                    createHiveplotBox(data, selector + " .hiveplot-container");
                    createSchematicPlot(data, selector + " .schematic-con-container"); //.schematic-container
                    createSchematicPlot(data, selector + " .schematic-non-container", {isContiguousPlot: false});
                    createNGLview("single",pdb[0]);
                }
            });
        }
    });
}

function initializeGoButtonTwoCrystalGroups(selector, heatmapFunction, generic=false) {
    $(selector + ' .go-button').click(function() {
        var pdbs1 = JSON.parse($(selector + ' .crystal-group-1-pdbs').val());
        var pdbs2 = JSON.parse($(selector + ' .crystal-group-2-pdbs').val());
        //var segments = JSON.parse($(selector + ' .segments-input').val());
        var segments = ['TM1','TM2','TM3','TM4','TM5','TM6','TM7','TM1','ICL1','ECL1','ICL2','ECL2','ICL3','ECL3','N-term','C-term'];
        if (pdbs1.length > 0 && pdbs2.length > 0 && segments.length > 0) {
            var interactionTypes = JSON.parse($(selector + ' .interactions-input').val());
            $(".heatmap").hide();
            $(".heatmap-legend").hide();
            $(".matrix-tab:visible").click();
            $(selector + ' .heatmap-container').append('<span id=svgloading>Loading... (0%)</span>');

            $.getJSON( '/contactnetwork/interactiondata',
            {
                'segments': segments,
                'generic': generic,
                'pdbs': pdbs1,
                'interaction_types': interactionTypes
            },
            function( data1 ) {
                $("#svgloading").remove();
                $(selector + ' .heatmap-container').append('<span id="svgloading">Loading... (33%)</span>');
                $.getJSON( '/contactnetwork/interactiondata',
                {
                    'segments': segments,
                    'generic': generic,
                    'pdbs': pdbs2,
                    'interaction_types': interactionTypes
                },
                function( data2 ) {
                    $("#svgloading").remove();
                    $(selector + ' .heatmap-container').append('<span id="svgloading">Loading... (66%)</span>');
                    // NOTE: this call seems redundant, shouldn't we already have all the data we need from the previous two calls
                    $.getJSON( '/contactnetwork/interactiondata',
                    {
                        'segments': segments,
                        'generic': generic,
                        'pdbs': pdbs1.concat(pdbs2),
                        'interaction_types': interactionTypes
                    }, function ( data3 ) {
                        // Re-render heatmap
                        $(".heatmap").show();
                        $(".heatmap-legend").show();
                        $("#svgloading").remove()
                        heatmapFunction(data1, data2, data3, selector + ' .heatmap-container');
                        createSchematicPlot(data3, selector + " .schematic-con-container", {type: 'twoCrystalGroups'}, data1, data2); //.schematic-container
                        createSchematicPlot(data3, selector + " .schematic-non-container", {isContiguousPlot: false, type: 'twoCrystalGroups'}, data1, data2)

                        // Re-render flareplot
                        createTwoGroupFlareplotBox(data1, data2, data3, selector + " .flareplot-container");
                        createNGLview("two-groups",pdbs1[0], pdbs1.concat(pdbs2));
                    });
                });
            });
        }
    });
}

function initializeFullscreenButton(selector) {
    var fullScreenElement = $(selector + ' .heatmap-container').get(0);
    $(selector + ' .btn-fullscreen').click(function() {
        toggleFullScreen(fullScreenElement);
    });
}


function update_text_in_modal() {
    var mode = $('ul#mode_nav').find('li.active').find('a').text().trim();

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
                        { "orderDataType": "dom-checkbox" }
                    ]
        });

        yadcf.init(oTable[mode],
        [
            {
                column_number : 0,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Class",
                filter_reset_button_text: false,
            },
            {
                column_number : 1,
                filter_type: "multi_select",
                select_type: 'select2',
                select_type_options: {
                    width: '70px'
                },
                filter_default_label: "PDB",
                filter_reset_button_text: false,
            },
            {
                column_number : 2,
                filter_type: "multi_select",
                select_type: 'select2',
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "Receptor",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
            },
            {
                column_number : 3,
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
                select_type_options: {
                    minimumResultsForSearch: -1 // remove search box
                },
                filter_default_label: "State",
                column_data_type: "html",
                html_data_type: "text",
                filter_match_mode : "exact",
                filter_reset_button_text: false,

            },
            {
                column_number : 6,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Representative",
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
    });

    };
}

function createTwoGroupFlareplotBox(data1, data2, data3, container, toggle = false) {
    // prepare two group data for visualization
    var data = data3;

    // frequency + count holder
    data["frequency"] = {};
    data["count"] = {};

    Object.keys(data.interactions).forEach(function(pair) {
    var f1 = 0, f2 = 0, c1 = 0, c2 = 0;;
    if (pair in data1.interactions) {
        c1 = Object.keys(data1.interactions[pair]).length;
        f1 = c1/data1.pdbs.length;
    }
    if (pair in data2.interactions) {
        c2 = Object.keys(data2.interactions[pair]).length;
        f2 = c2/data2.pdbs.length;
    }
    var f3 = f1 - f2;
    var c3 = c1 + c2;
    data["frequency"][pair] = [f1, f2, f3];
    data["count"][pair] = [c1, c2, c3];
    });

    createFlareplotBox(data, container, toggle = false);
}

var flareplot = {};
var contiguous = false;
var interactionsToggleList = [];
function createFlareplotBox(data, container, toggle = false) {
    // clean
    if (toggle){
        $(container).children().last().remove();
    } else {
        // in case refresh with new parameters => reset
        contiguous = false;
        $(container).html("");
    }

    // add menu
    if (!toggle){
        var newDiv = document.createElement("div");
        newDiv.setAttribute("class", "flareplot-legend");

        var content = '<div class="controls">'
                            +'<h3>Controls</h3>';

        // only possible with more than 4 segments, otherwise it will become a mess
        if (data.segments.length > 4)
            content += '<p>Contacts contiguous segments outwards: <input type=checkbox id="flareplot_contiguous"></p>';

        content += '<p>Line colors: <select id="flareplot_color">'
                +'<option value="none">None (gray)</option>'
                +'<option value="segment">GPCR segment</option>';

        // if single structure - use interaction coloring
        if (container.indexOf("single-crystal-tab") >= 0)
            content += '<option value="interactions" selected>Interaction Type</option>'
        // if group(s) of structures - use frequency coloring (gradient)
        else
            content += '<option value="frequency" selected>Interaction Frequency/Count</option>'

        content += '</select></p>';


        // Populate heatmap legend
        if (container.indexOf("single-crystal-tab") >= 0) {
            // TODO Optimize/generalize interaction set selection - code duplication
            var interactionTypes = new Set(["Polar", "Aromatic", "Hydrophobic", "Van der Waals"]);

            // Add interactions color legend
            content += '<h3>Toggle interactions</h3><ul>';

            interactionTypes = Array.from(interactionTypes).sort(function (i1, i2) {
                return getInteractionStrength(i2) - getInteractionStrength(i1);
            });

            interactionTypes.forEach(function(i) {
                var rgb = getInteractionColor(i, false);
                content += '<li>'
                        + '<div class="color-box" style="background-color: ' + rgb + '">' + '<input type="checkbox" data-interaction-type="' + i.replace(/ /g,"-") +'"></input>' + '</div><p>' + i + '</p>'
                        + '</li>';
            });
            content += '</ul></div>';
        } else if (container.indexOf("single-crystal-group-tab") >= 0) {
            // slider 0 to count
            content += '<h4 class="center">Frequency (#PDBs)</h4>'
                + '<p>Range: <span id="pdbs-range-flare">0 - ' + data.pdbs.length + '</span></p>'
                + '<div class="slider-range" data-text-id="pdbs-range-flare" id="pdbs-range-flare-slider"></div>'
                + '<div class="temperature-scale">'
                + '<span class="gray-to-red"></span>'
                + '</div>';

            $( function() {
                $( container+" .slider-range" ).data({ "referenceContainer" : container })
                $( container+" .slider-range" ).slider({
                range: true,
                min: 0,
                max: data.pdbs.length,
                step: 1,
                values: [0,data.pdbs.length],
                slide: function( event, ui ) {
                    $( "#"+$(this).attr("data-text-id") ).html( ui.values[ 0 ] + " - " + ui.values[ 1 ] );
                    flareplot[$(this).data("referenceContainer")].updateRange(ui.values[ 0 ], ui.values[ 1 ]);
                }
                });
            } );
        } else {
            content += '<h4 class="center">Frequency</h4>'
                + '<p>Group 1 range: <span id="freq-flare-range-1">0 - 1</span></p>'
                + '<div class="slider-range" data-text-id="freq-flare-range-1" id="freq-flare-slider-range-1"></div>'
                + '<p>Group 2 range: <span id="freq-flare-range-2">0 - 1</span></p>'
                + '<div class="slider-range" data-text-id="freq-flare-range-2" id="freq-flare-slider-range-2"></div>'
                + '<p>Freq difference range: <span id="freq-flare-range-3">-1 - 1</span></p>'
                + '<div class="slider-range-diff" data-text-id="freq-flare-range-3" id="freq-flare-slider-range-3"></div>'
                + '<div class="temperature-scale">'
                + '<span class="red-to-gray"></span>'
                + '<span class="gray-to-blue"></span>'
                + '</div>'
                + '</div>';

                $( function() {
                    $( container+" #freq-flare-slider-range-1" ).data({ "referenceContainer" : container });
                    $( container+" #freq-flare-slider-range-2" ).data({ "referenceContainer" : container });
                    $( container+" #freq-flare-slider-range-3" ).data({ "referenceContainer" : container });
                    $( container+" #freq-flare-slider-range-1" ).slider({
                    range: true,
                    min: 0,
                    max: 1,
                    step: 0.01,
                    values: [0,1],
                    slide: function( event, ui ) {
                        $( "#"+$(this).attr("data-text-id") ).html( ui.values[ 0 ] + " - " + ui.values[ 1 ] );
                        updateTwoGroupSliders($(this)[0].id, ui);
                    }
                    });
                    $( container+" #freq-flare-slider-range-2" ).slider({
                    range: true,
                    min: 0,
                    max: 1,
                    step: 0.01,
                    values: [0,1],
                    slide: function( event, ui ) {
                        $( "#"+$(this).attr("data-text-id") ).html( ui.values[ 0 ] + " - " + ui.values[ 1 ] );
                        updateTwoGroupSliders($(this)[0].id, ui);
                    }
                    });
                    $( container+" #freq-flare-slider-range-3" ).slider({
                    range: true,
                    min: -1,
                    max: 1,
                    step: 0.01,
                    values: [-1,1],
                    slide: function( event, ui ) {
                        $( "#"+$(this).attr("data-text-id") ).html( ui.values[ 0 ] + " - " + ui.values[ 1 ] );
                        updateTwoGroupSliders($(this)[0].id, ui);
                    }
                    });
                });
        }

        newDiv.innerHTML = content;

        $(container).append(newDiv);

        $(container+" #flareplot_contiguous").click(function(e){
            $(function() {
                contiguous = !contiguous;
                createFlareplotBox(data, container, true);
            });
        });


        $(container+" #flareplot_color").data({ "referenceContainer" : container })
        $(container+" #flareplot_color").change(function(e){
            flareplot[$(this).data("referenceContainer")].updateColors($(this).val(), interactionsToggleList);
        });

        $(container + ' .flareplot-legend .color-box input[type=checkbox]').each(function() {
            $(this).prop('checked', true);

            // init toggle list
            interactionsToggleList.push($(this).data('interaction-type'));
            $(this).data({ "referenceContainer" : container })

            $(this).change(function() {
                // toggle interactions in flareplot
                interactionsToggleList = [];

                $(container + ' .flareplot-legend input[type=checkbox]').each(function() {
                    if ($(this).prop('checked'))
                        interactionsToggleList.push($(this).data('interaction-type'));
                });

                // toggle interactions in flareplot
                flareplot[$(this).data("referenceContainer")].showInteractions(interactionsToggleList);
            });
        });
    }

    // create flareplot
    flareplot[container] = createFlareplot(800, parseGPCRdb2flare(data), container, contiguous);

    // update coloring and visibility if toggled
    if (toggle) {
        var range = $( container+" .slider-range" ).slider("values");
        flareplot[container].updateRange(range[0], range[1]);

        flareplot[container].updateColors($(container+" #flareplot_color").val(), interactionsToggleList);
        if ($(container + ' .flareplot-legend .color-box input[type=checkbox]').length > 0)
            flareplot[container].showInteractions(interactionsToggleList);
        if (container.indexOf("two-crystal-groups-tab") >= 0)
            updateTwoGroupSliders("skip", []);
    } else {
        flareplot[container].updateColors($(container+" #flareplot_color").val(), interactionsToggleList);
    }
}

function createHiveplotBox(data, container) {
    // clean contents
    $(container).html("");

    createHiveplot(data, container);
}

function updateTwoGroupSliders(origin, ui) {
    // grab data from all sliders
    var g1 = $( "#freq-flare-slider-range-1").slider("values");
    var g2 = $( "#freq-flare-slider-range-2" ).slider("values");
    var diff = $( "#freq-flare-slider-range-3" ).slider("values");

    // update with current slide values
    if (origin == "freq-flare-slider-range-1")
        g1 = ui.values;
    else if (origin == "freq-flare-slider-range-2")
        g2 = ui.values;
    else if (origin == "freq-flare-slider-range-3")
        diff = ui.values;

    // apply cutoffs
    flareplot[$( "#freq-flare-slider-range-1").data("referenceContainer")].updateRangeTwoGroups(g1[0], g1[1], g2[0], g2[1], diff[0], diff[1]);
}

var stage = [];
var color_schemes = [];
var schemeId_grey
function createNGLview(mode,pdb, pdbs = false) {
    var gpcr_rep
    $("#ngl-"+mode).html("");
    stage[mode] = new NGL.Stage( "ngl-"+mode, { backgroundColor: "white" } );

    var pdb_data
    var original_o
    var blue_colors = ['#f7fcf0','#e0f3db','#ccebc5', '#a8ddb5',    '#7bccc4',    '#4eb3d3', '#2b8cbe',    '#0868ac',    '#084081']
    var reps = {} // store ngl representations

    $.getJSON( "pdb/"+pdb,
        function( data ) {
        var highlight = ['TM1', 'TM2', 'TM3', 'TM4', 'TM5', 'TM6', 'TM7', 'H8'];
        var segments_sets = {}
        highlight.forEach( function(e){
            segments_sets[e] = ((e in data['segments']) ? data['segments'][e].join(", ") : "")
        });

        pdb_data = data;
        color_schemes['blue'] = NGL.ColormakerRegistry.addSelectionScheme([
                [blue_colors[1], segments_sets[highlight[0]]],
                [blue_colors[2], segments_sets[highlight[1]]],
                [blue_colors[3], segments_sets[highlight[2]]],
                [blue_colors[4], segments_sets[highlight[3]]],
                [blue_colors[5], segments_sets[highlight[4]]],
                [blue_colors[6], segments_sets[highlight[5]]],
                [blue_colors[7], segments_sets[highlight[6]]],
                [blue_colors[8], segments_sets[highlight[7]]],
                [ "white", "*" ]
                ])

        color_schemes['grey'] = NGL.ColormakerRegistry.addSelectionScheme([
                    ["#ccc", segments_sets[highlight[0]]],
                    ["#bbb", segments_sets[highlight[1]]],
                    ["#aaa", segments_sets[highlight[2]]],
                    ["#888", segments_sets[highlight[3]]],
                    ["#666", segments_sets[highlight[4]]],
                    ["#444", segments_sets[highlight[5]]],
                    ["#333", segments_sets[highlight[6]]],
                    ["#111", segments_sets[highlight[7]]],
                    [ "white", "*" ]
                    ])


        var stringBlob = new Blob( [ pdb_data['pdb'] ], { type: 'text/plain'} );
        stage[mode].loadFile( stringBlob, { ext: "pdb" }  ).then( function( o ){
            original_o = o

            gpcr_rep = o.addRepresentation( "cartoon", {
            sele: ":"+pdb_data['chain'],
            // radiusType: '',
            radiusSize: 1,
            radiusScale: 0.7,
            // color: "atomindex",
            // colorScale: "Accent",
            // color: "residueindex",
            // colorScale: "greys",
            color: color_schemes['grey'],
            metalness: 0,
            colorMode: "hcl",
            roughness: 1,
            opacity: 1,
            depthWrite: true
            });

            var links = []
            var links_gn = []
            var res_int = []

            if (mode=='single') {
                $('#single-crystal-tab rect[data-interaction-type]').each(function(e) {
                    var rect = $(this);
                    var resNo1 = rect.data('res-no-1');
                    var resNo2 = rect.data('res-no-2');
                    var seg1 = rect.data('seg-1');
                    var seg2 = rect.data('seg-2');
                    var genNo1 = rect.data('gen-no-1');
                    var genNo2 = rect.data('gen-no-2');
                    var aa1 = rect.data('aa-1');
                    var aa2 = rect.data('aa-2');
                    var iType = rect.data('interaction-type');
                    if (getInteractionStrength(iType)<3) return
                    res_int.push(resNo1);
                    res_int.push(resNo2);
                    links.push({"atoms": [resNo1+":"+pdb_data['chain']+".CA",resNo2+":"+pdb_data['chain']+".CA"], "data":{"color":getInteractionColor(iType)}, "resID":resNo1+"-"+resNo2})

                    if ((genNo1=='-') || (genNo2=='-')) return
                    links_gn.push({"atoms": [resNo1+":"+pdb_data['chain']+".CA",resNo2+":"+pdb_data['chain']+".CA"], "data":{"color":getInteractionColor(iType)}, "resID":resNo1+"-"+resNo2})
                });
            } else if (mode=='single-group') {
                // single group
                $('#single-crystal-group-tab .heatmap-container rect[data-frequency]').each(function(e) {
                    var rect = $(this);
                    var genNo1 = rect.data('gen-no-1');
                    var genNo2 = rect.data('gen-no-2');
                    var seg1 = rect.data('seg-1');
                    var seg2 = rect.data('seg-2');
                    var nInteractions = rect.data('num-interactions');
                    var nTotalInteractions = rect.data('total-possible-interactions');
                    var frequency = rect.data('frequency');

                    if ((genNo1=='-') || (genNo2=='-')) return

                    // Adjust GN numbering to the shown structure
                    var resNo1 = pdb_data['only_gn'][pdb_data['gn_map'].indexOf(genNo1)];
                    var resNo2 = pdb_data['only_gn'][pdb_data['gn_map'].indexOf(genNo2)];

                    if ((typeof resNo1=='undefined') || (typeof resNo2=='undefined')) return

                    // Push interactions
                    res_int.push(resNo1);
                    res_int.push(resNo2);

                    links_gn.push({"atoms": [resNo1+":"+pdb_data['chain']+".CA",resNo2+":"+pdb_data['chain']+".CA"], "data":{"color":getFrequencyColor(frequency)}, "resID":resNo1+"-"+resNo2})
                });
                links = links_gn;
            } else {
                // two groups
                // single group
                $('#two-crystal-groups-tab .heatmap-container rect[data-frequency-diff]').each(function(e) {
                    var rect = $(this);
                    var genNo1 = rect.data('gen-no-1');
                    var genNo2 = rect.data('gen-no-2');
                    var seg1 = rect.data('seg-1');
                    var seg2 = rect.data('seg-2');
                    var nInteractions = rect.data('num-interactions');
                    var nTotalInteractions = rect.data('total-possible-interactions');
                    var frequency = rect.data('frequencyDiff');

                    if ((genNo1=='-') || (genNo2=='-')) return

                    // Adjust GN numbering to the shown structure
                    var resNo1 = pdb_data['only_gn'][pdb_data['gn_map'].indexOf(genNo1)];
                    var resNo2 = pdb_data['only_gn'][pdb_data['gn_map'].indexOf(genNo2)];

                    if ((typeof resNo1=='undefined') || (typeof resNo2=='undefined')) return

                    // Push interactions
                    res_int.push(resNo1);
                    res_int.push(resNo2);

                    links_gn.push({"atoms": [resNo1+":"+pdb_data['chain']+".CA",resNo2+":"+pdb_data['chain']+".CA"], "data":{"color":getFrequencyColor(-1*frequency)}, "resID":resNo1+"-"+resNo2})
                });
                links = links_gn;
            }

            var linkMap = {}
            links.forEach(function (link) {
                linkMap[link.resID] = link
            })

            var initColourSchemes = function () {
                var linkColourScheme = function () {
                this.bondColor = function (b) {
                    var origLink = linkMap[b.atom1.resno + "-" + b.atom2.resno]
                    if (origLink) {
//                             console.log('FOUND!',origLink,origLink.data.color);
                    r = origLink.data.color
                    return (r["r"] << 16) + (r["g"] << 8) + r["b"]
                    }
                    return (8 << 16) + (8 << 8) + 8 // (128 << 16) + (128 << 8) + 128 // grey default
                }
                }
                reps.linkColourScheme = NGL.ColormakerRegistry.addScheme(linkColourScheme, "xlink")
            }
            initColourSchemes()

            // return unique atom indices as a selection from a set of pairs of atom indices
            // var makeAtomSelection = function (someLinks) {
            //   var atomSet = new Set()
            //   someLinks.forEach(function (link) {
            //     atomSet.add(link.atoms[0].split(".")[0])
            //     atomSet.add(link.atoms[1].split(".")[0])
            //   })
            //   var atomSelection = "(" + Array.from(atomSet).join(", ") + ") and .CA"
            //   return atomSelection
            // }
            // var startAtomSel = makeAtomSelection(links)
            // console.log(startAtomSel,"startAtomSel")

            var baseLinkScale = 20
            reps.links = o.addRepresentation("distance", {
                atomPair: links.map(function (l) {
                return l.atoms
                }),
                // colorValue: "yellow",
                colorScheme: reps.linkColourScheme,
                useCylinder: false,
                labelVisible: false,
                linewidth: 2
            })
            reps.links_gn = o.addRepresentation("distance", {
                atomPair: links_gn.map(function (l) {
                return l.atoms
                }),
                // colorValue: "yellow",
                colorScheme: reps.linkColourScheme,
                useCylinder: false,
                labelVisible: false,
                linewidth: 2,
                visible: false
            })

            reps.int_res = o.addRepresentation( "spacefill", {
                sele: ":"+pdb_data['chain']+" and ("+res_int.join(", ")+") and (.CA)",
                color: "red",
                // colorScale: ["#44f", "#444"],
                radiusScale: .2,
                name: "res",
                visible: false
            });


            res_int_gn = Object.assign([], res_int);
            res_int_gn = intersect(res_int_gn, pdb_data['only_gn']);
            reps.int_res_gn = o.addRepresentation( "spacefill", {
                sele: ":"+pdb_data['chain']+" and ("+res_int_gn.join(", ")+") and (.CA)",
                color: "red",
                // colorScale: ["#44f", "#444"],
                radiusScale: .2,
                name: "res",
                visible: false
            });

            reps.ball_all = o.addRepresentation("ball+stick", {
                sele: ":"+pdb_data['chain']+" and sidechainAttached",
                visible: false
                })

            reps.ball = o.addRepresentation("ball+stick", {
                sele: ":"+pdb_data['chain']+" and ("+pdb_data['only_gn'].join(", ")+") and sidechainAttached",
                visible: false
                })

            reps.ball_int = o.addRepresentation("ball+stick", {
                sele: ":"+pdb_data['chain']+" and ("+res_int.join(", ")+") and sidechainAttached",
                visible: false
                })

            reps.ball_int_gn = o.addRepresentation("ball+stick", {
                sele: ":"+pdb_data['chain']+" and ("+res_int_gn.join(", ")+") and sidechainAttached",
                visible: false
                })


            // reps.residues = o.addRepresentation("spacefill", {
            //   sele: ".CA",
            //   color: "#000",
            //   // colorScale: ["#44f", "#444"],
            //   radiusScale: .2,
            //   name: "res"
            // })


            // reps.residues = o.addRepresentation("spacefill", {
            //   sele: startAtomSel,
            //   color: "#ccc",
            //   // colorScale: ["#44f", "#444"],
            //   radiusScale: .2,
            //   name: "res"
            // })


            // o.addRepresentation( "distance", {
            //   atomPair: links.map(function (l) {
            //     return l.atoms
            //   }),
            //   radiusScale: 1,
            //   labelVisible: false,
            //   useCylinder: true,
            // } );


            reps.ngl_contacts = o.addRepresentation("contact", {
                sele: ":"+pdb_data['chain']+" and ("+pdb_data['only_gn'].join(", ")+")",
                radiusSize: 0.07,
                weakHydrogenBond: false,
                waterHydrogenBond: false,
                backboneHydrogenBond: false,
                visible:false
            })

            o.autoView();


        } );

    });

    var newDiv = document.createElement("div");
    newDiv.setAttribute("style", "position: absolute; top: 50px; left: 20px")
    var controls = '<div class="controls">'
                    + '<h3>Controls</h3>';

    if (pdbs){
            controls += '<p>Structure: <select id="ngl_pdb_'+mode+'_ref">';
            for (var i = 0; i < pdbs.length; i++){
                if (pdbs[i]==pdb)
                    controls += '<option value="'+pdbs[i]+'" SELECTED>'+pdbs[i]+'</option>';
                else
                    controls += '<option value="'+pdbs[i]+'">'+pdbs[i]+'</option>';
            }
            controls += '</select></p>';
    }

    controls += '<p>Colors: <select id="ngl_color"><option value="grey">greys</option><option value="blue">blue</option></select></p>'
                        +'<p>Only GNs: <input type=checkbox id="ngl_only_gns"></p>'
                        +'<p>Highlight interacting res: <input type=checkbox id="highlight_res"></p>'
                        +'<p>Hide interaction lines: <input type=checkbox id="toggle_interactions"></p>'
                        +'<p>Show all side-chains: <input type=checkbox id="toggle_sidechains"></p>'
                        +'<p>Show interacting side-chains: <input type=checkbox id="toggle_sidechains_int"></p>'
//                              +'<p>Show NGL derived contacts: <input type=checkbox id="ngl_contacts"></p>'
                        +'</div>';

    newDiv.innerHTML = controls;

    $("#ngl-"+mode).append(newDiv);

    $("#ngl_pdb_"+mode+"_ref").change(function(e){
        createNGLview(mode, $(this).val(), pdbs);
    });

    $("#"+mode+"-NGL-link").click(function(e){
        $(function() {
            stage[mode].handleResize();
        });
    });


    $("#ngl-"+mode+" #ngl_color").change(function(e){
        gpcr_rep.setParameters({
            color: color_schemes[$(this).val()]
        });
    });

    // TODO: cleanup all NGL toggles and make generic function to handle toggle-dependent showing/hiding
    $("#ngl-"+mode+" #ngl_only_gns").change(function(e){
        if ($(this).prop('checked')) {
            sele = ":"+pdb_data['chain']+" and ("+pdb_data['only_gn'].join(", ")+")";
            // toggle edges
            if (!$("#ngl-"+mode+" #toggle_interactions").prop('checked')){
            reps.links.setVisibility(false);
            reps.links_gn.setVisibility(true);
            }
            // toggle CA spheres
            if ($("#ngl-"+mode+" #highlight_res").prop('checked')){
            reps.int_res_gn.setVisibility(true);
            reps.int_res.setVisibility(false);
            }
            // toggle sidechains
            if ($("#ngl-"+mode+" #toggle_sidechains").prop('checked')){
            reps.ball.setVisibility(true);
            reps.ball_all.setVisibility(false);
            }
            // toggle interacting toggle_sidechains
            if ($("#ngl-"+mode+" #toggle_sidechains_int").prop('checked')){
            reps.ball_int_gn.setVisibility(true);
            reps.ball_int.setVisibility(false);
            }
        } else {
            sele = ":"+pdb_data['chain'];
            // toggle edges
            if (!$("#ngl-"+mode+" #toggle_interactions").prop('checked')){
            reps.links.setVisibility(true);
            reps.links_gn.setVisibility(false);
            }
            // toggle CA spheres
            if ($("#ngl-"+mode+" #highlight_res").prop('checked')){
            reps.int_res_gn.setVisibility(false);
            reps.int_res.setVisibility(true);
            }

            // toggle sidechains
            if ($("#ngl-"+mode+" #toggle_sidechains").prop('checked')){
            reps.ball.setVisibility(false);
            reps.ball_all.setVisibility(true);
            }
            // toggle interacting toggle_sidechains
            if ($("#ngl-"+mode+" #toggle_sidechains_int").prop('checked')){
            reps.ball_int_gn.setVisibility(false);
            reps.ball_int.setVisibility(true);
            }
        }

        gpcr_rep.setSelection(sele);
        original_o.autoView();
    });

    $("#ngl-"+mode+" #highlight_res").change(function(e){
        // TODO check the GN toggle
        if ($(this).prop('checked')) {
            if ($("#ngl-"+mode+" #ngl_only_gns").prop('checked'))
            reps.int_res_gn.setVisibility(true);
            else
            reps.int_res.setVisibility(true);
        } else {
            reps.int_res_gn.setVisibility(false);
            reps.int_res.setVisibility(false);
        }
    });

    $("#ngl-"+mode+" #toggle_interactions").change(function(e){
        if ($(this).prop('checked')) {
            reps.links.setVisibility(false);
            reps.links_gn.setVisibility(false);
        } else {
            if ($("#ngl-"+mode+" #ngl_only_gns").prop('checked'))
            reps.links_gn.setVisibility(true);
            else
            reps.links.setVisibility(true);
        }
    });

    $("#ngl-"+mode+" #toggle_sidechains").change(function(e){
        if ($(this).prop('checked')) {
            if ($("#ngl-"+mode+" #ngl_only_gns").prop('checked'))
            reps.ball.setVisibility(true);
            else
            reps.ball_all.setVisibility(true);
        } else {
            reps.ball_all.setVisibility(false);
            reps.ball.setVisibility(false);
        }
    });

    $("#ngl-"+mode+" #toggle_sidechains_int").change(function(e){
        if ($(this).prop('checked')) {
            if ($("#ngl-"+mode+" #ngl_only_gns").prop('checked'))
            reps.ball_int_gn.setVisibility(true);
            else
            reps.ball_int.setVisibility(true);
        } else {
            reps.ball_int.setVisibility(false);
            reps.ball_int_gn.setVisibility(false);
        }
    });

    $("#ngl-"+mode+" #ngl_contacts").change(function(e){
        if ($(this).prop('checked')) {
            reps.ngl_contacts.setVisibility(true);
        } else {
            reps.ngl_contacts.setVisibility(false);
        }
    });

}

function intersect(a, b) {
    var t;
    if (b.length > a.length) t = b, b = a, a = t; // indexOf to loop over shorter
    return a.filter(function (e) {
        return b.indexOf(e) > -1;
    });
}

$('#single-crystal-pdb-modal-table').on('shown.bs.modal', function (e) {
    showPDBtable('#single-crystal-pdb-modal-table');
})
$('#single-crystal-group-pdbs-modal-table').on('shown.bs.modal', function (e) {
    showPDBtable('#single-crystal-group-pdbs-modal-table');
})
$('#two-crystal-group-pdbs-modal-1-table').on('shown.bs.modal', function (e) {
    showPDBtable('#two-crystal-group-pdbs-modal-1-table');
})
$('#two-crystal-group-pdbs-modal-2-table').on('shown.bs.modal', function (e) {
    showPDBtable('#two-crystal-group-pdbs-modal-2-table');
})


function initializePdbChooserTables() {
    $.get('pdbtabledata', function ( data ) {
    $('#single-crystal-pdb-modal-table .tableview').html(data);
    $('#single-crystal-group-pdbs-modal-table .tableview').html(data);
    $('#two-crystal-group-pdbs-modal-1-table .tableview').html(data);
    $('#two-crystal-group-pdbs-modal-2-table .tableview').html(data);
    pdbtabledata = data;
    });
}


function initalizeSingleCrystalView() {
//            initializeSegmentButtons('#single-crystal-tab');
    initializeGoButton('#single-crystal-tab', renderSingleCrystalHeatmap);
    initializeFullscreenButton('#single-crystal-tab');
}

function initializeSingleGroupCrystalView() {
//            initializeSegmentButtons('#single-crystal-group-tab');
    initializeGoButton('#single-crystal-group-tab', renderSingleCrystalGroupHeatmap, true);
    initializeFullscreenButton('#single-crystal-group-tab');
    initializeInteractionButtons('#single-crystal-group-tab');
}

function initializeTwoCrystalGroupsView() {
//            initializeSegmentButtons('#two-crystal-groups-tab');
    initializeGoButtonTwoCrystalGroups('#two-crystal-groups-tab', renderTwoCrystalGroupsHeatmap, true);
    initializeFullscreenButton('#two-crystal-groups-tab');
    initializeInteractionButtons('#two-crystal-groups-tab');
}
$(document).ready(function() {
    // Get PDBs for table build
    initializePdbChooserTables();

    // Single PDB files
    initalizeSingleCrystalView();

    // Single group of PDB files
    initializeSingleGroupCrystalView();

    // Two groups of PDB files
    initializeTwoCrystalGroupsView();
}); 

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

function populateTable(selector, data) {
    console.log(data)
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
        setTimeout(renderTable,100,"#"+selector,data);
    });

}

function initializeGoButton(selector, generic=false) {
    $(selector + ' .go-button').click(function() {
        var pdb = JSON.parse($(selector + ' .crystal-pdb').val());
        //var segments = JSON.parse($(selector + ' .segments-input').val());
        var segments = ['TM1','TM2','TM3','TM4','TM5','TM6','TM7','TM1','ICL1','ECL1','ICL2','ECL2','ICL3','ECL3','N-term','C-term'];
        if (pdb.length > 0 && segments.length > 0) {
            var interactionTypes = JSON.parse($(selector + ' .interactions-input').val());

            if (!$(selector + ' .interactions-input').val() == null)
                interactionTypes = JSON.parse($(selector + ' .interactions-input').val());
                
            renderTable(pdb);

            if (selector == '#single-crystal-group-tab') {
                createNGLview("single-group",pdb[0], pdb);
            } else {
                createNGLview("single",pdb[0]);
            }
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
    var ReceptorName = $(elem).attr('long');
    var pdbName = $(elem).attr('id');
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

var stage = [];
var color_schemes = [];
var schemeId_grey
function createNGLview(mode,pdb, pdbs = false) {
    var gpcr_rep
    $("#ngl-"+mode).html("");
    stage[mode] = new NGL.Stage( "ngl-"+mode, { backgroundColor: "white" } );

    var pdb_data
    var blue_colors = ['#f7fcf0','#e0f3db','#ccebc5', '#a8ddb5',    '#7bccc4',    '#4eb3d3', '#2b8cbe',    '#0868ac',    '#084081']
    var reps = {} // store ngl representations
    var original_o

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
        
        var angle_color = [];
        var colorstring
        
        // TODO: median angle, difference to median angle
        // groups
        
        residue_data.forEach(function(e){
            
            if(e[2]<60){
                colorstring = "#" + "0000" + Math.floor((180-e[2])*1.41).toString(16) 
            } else if (e[2] < 120) {
                colorstring = "#" + Math.floor(e[2]*1.41).toString(16) + Math.floor((180-Math.abs(e[2]-90))*1.41).toString(16) + Math.floor((180-e[2])*1.41).toString(16)
            } else {
                colorstring = "#" + Math.floor(e[2]*1.41).toString(16) + "0000"
            }
            // +Math.floor(e[2]*1.41).toString(16)+Math.floor(e[2]*1.41).toString(16)+Math.floor(e[2]*1.41).toString(16)
            angle_color.push([colorstring , ""+e[1]])
        });
        
        
        angle_color.push([ "white", "*" ]);
        
        color_schemes['angles'] = NGL.ColormakerRegistry.addSelectionScheme(angle_color)


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


            reps.ngl_contacts = o.addRepresentation("contact", {
                sele: ":"+pdb_data['chain']+" and ("+pdb_data['only_gn'].join(", ")+")",
                radiusSize: 0.07,
                weakHydrogenBond: false,
                waterHydrogenBond: false,
                backboneHydrogenBond: false,
                visible:false
            })

            o.autoView();
            
            // mousover and click on datatable row to highlight residue in NGL viewer
            var temprepr
            $("#single-table-tab-table tbody").on("mouseover", "tr", function(event){
                temprepr = o.addRepresentation("ball+stick", {sele: ""+residuetable.row(this).data()[1]});
            }).mouseout(function(event){
                o.removeRepresentation(temprepr)
            });
            
            var repr_dict = {}
            $("#single-table-tab-table tbody").on("click", "tr", function(event){
                if(residuetable.row(this).data()[1] in repr_dict){
                    o.removeRepresentation(repr_dict[residuetable.row(this).data()[1]])
                    delete repr_dict[residuetable.row(this).data()[1]]
                    console.log(this)
                    $(this).removeClass("table-selected")
                }else{
                    repr_dict[residuetable.row(this).data()[1]] = o.addRepresentation("ball+stick", {sele: ""+residuetable.row(this).data()[1]});
                    $(this).addClass("table-selected")
                }
            });


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

    controls += '<p>Colors: <select id="ngl_color"><option value="grey">greys</option><option value="blue">blue</option><option value="angles">angles</option></select></p>'
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
        } else {!$("#ngl-"+mode+" #toggle_interactions").prop('checked')
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


var residue_data

function renderTable(pdb) {
    console.log(pdb)
    $.get('angledat?pdbs[]='+pdb[0], function(newDataArray) {
    console.log()
    residue_data = newDataArray["data"]
    residuetable.clear();
    residuetable.rows.add(residue_data);
    residuetable.draw();
    });

}


function initializeResidueTable() {
    residuetable = $("#single-table-tab-table").DataTable();
    yadcf.init(residuetable,
    [
        {
            column_number : 0,
            filter_type: "multi_select",
            select_type: 'select2',
            filter_default_label: "Generic number",
            filter_reset_button_text: false,
        },
        {
            column_number : 1,
            filter_type: "multi_select",
            select_type: 'select2',
            filter_default_label: "angle",
            filter_reset_button_text: false,
        },
        {
            column_number : 2,
            filter_type: "multi_select",
            select_type: 'select2',
            filter_default_label: "angle",
            filter_reset_button_text: false,
        }
    ],
    {
        cumulative_filtering: false
    });

    yadcf.exResetAllFilters(residuetable);
    
//     $("#single-table-tab-table tbody").on("mouseover", "tr", function(event){
//         original_o.addRepresentation("ball+stick", ""+{sele: residuetable.row(this).data()[1]});
//     });
    
//     original_o.addRepresentation("ball+stick", {sele: residuetable.row(this).data()[1]});
}

$('#single-crystal-pdb-modal-table').on('shown.bs.modal', function (e) {
    showPDBtable('#single-crystal-pdb-modal-table');
})

function initializePdbChooserTables() {
    $.get('/contactnetwork/pdbtabledata', function ( data ) {
    $('#single-crystal-pdb-modal-table .tableview').html(data);
    pdbtabledata = data;
    });
}


function initalizeSingleCrystalView() {
    initializeResidueTable()
    initializeGoButton('#single-crystal-tab');
    initializeFullscreenButton('#single-crystal-tab');
}

var selectortable
var residuetable


$(document).ready(function() {
    // Get PDBs for table build
    initializePdbChooserTables();

    // Single PDB files
    initalizeSingleCrystalView();
}); 

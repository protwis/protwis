window.zoomHeatmap = {};
window.zoomHiveplot = {};

// create NGL tooltip
var ngl_tooltip = document.createElement("div")
Object.assign(ngl_tooltip.style, {
    display: "none",
    position: "fixed",
    zIndex: 10,
    pointerEvents: "none",
    backgroundColor: "rgba( 0, 0, 0, 0.6 )",
    color: "lightgrey",
    padding: "8px",
    fontFamily: "sans-serif"
})
document.body.appendChild(ngl_tooltip)

var is_fullscreen = false;

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
    $('.popover').each(function () {
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
        case 0:
            r = v, g = t, b = p;
            break;
        case 1:
            r = q, g = v, b = p;
            break;
        case 2:
            r = p, g = v, b = t;
            break;
        case 3:
            r = p, g = q, b = v;
            break;
        case 4:
            r = t, g = p, b = v;
            break;
        case 5:
            r = v, g = p, b = q;
            break;
    }
    return {
        r: Math.round(r * 255),
        g: Math.round(g * 255),
        b: Math.round(b * 255)
    };
}

function rgb2hex(r, g, b) {
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
    switch (i.toLowerCase()) {
        case "ionic":
            return 5;
        case "polar":
            return 4;
        case "aromatic":
            return 3;
        case "hydrophobic":
            return 2;
        case "van-der-waals":
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
    return getGradientColor(-1 * frequency, rgb);
}

function getFlareGradientColor(fDiff, rgb = true, hide_zero = true) {
    var color;
    var shift = 80;
    var basal = 255 - shift;
    if (hide_zero && fDiff == 0) {
        // to hide zero, make paths completely white.
        return rgb2hex(255, 255, 255);
    }

    if (fDiff <= 0)
        // If fDiff is close to -1, we want a red color
        color = { r: basal + (fDiff * -1 * shift), g: basal - basal * (-fDiff), b: basal - basal * (-fDiff) };
    else
        // If fDiff is close to 1 we want a blue color
        color = { r: basal - basal * fDiff, g: basal - basal * fDiff, b: basal + (fDiff * shift) };

    if (rgb)
        return color;
    else
        return rgb2hex(color.r, color.g, color.b);
}

function getGradientColor(fDiff, rgb = true) {
    var color;
    if (fDiff <= 0)
        // If fDiff is close to -1, we want a red color
        color = { r: 255, g: 255 - 255 * (-fDiff), b: 255 - 255 * (-fDiff) };
    else
        // If fDiff is close to 1 we want a blue color
        color = { r: 255 - 255 * fDiff, g: 255 - 255 * fDiff, b: 255 };

    if (rgb)
        return color;
    else
        return rgb2hex(color.r, color.g, color.b);
}

function getStrongestInteractionType(interactions) {
    if ($.inArray('ionic', interactions) > -1)
        return 'ionic';
    else if ($.inArray('polar', interactions) > -1)
        return 'polar';
    else if ($.inArray('aromatic', interactions) > -1)
        return 'aromatic';
    else if ($.inArray('hydrophobic', interactions) > -1)
        return 'hydrophobic';
    else if ($.inArray('van-der-waals', interactions) > -1)
        return 'van-der-waals';

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
        Object.keys(obj[key]).forEach(function (k, index) {
            interactions.add(obj[key][k]);
        });
        // if (Object.prototype.hasOwnProperty.call(obj, key)) {
        //     for (var k in obj[key])
        //         interactions.add(obj[key][k]);
        // }
    }

    // Sort according to strength - now done by backend
    interactions = Array.from(interactions);
    interactions.sort(function (i1, i2) {
        return getInteractionStrength(i1) - getInteractionStrength(i2);
    });

    return interactions;
}


function getInteractionColor(interaction, rgb = true) {
    var r, g, b;

    value = interaction
    if (isNaN(value))
        value = value.toLowerCase()

    switch (value) {
        case 'ionic':
        case 5:
            r = 197;
            g = 66;
            b = 244;
            break;
        case 'polar':
        case 4:
            //r = 254; g = 0; b = 16;
            r = 255;
            g = 98;
            b = 108;
            break;
        case 'aromatic':
        case 3:
            //r = 94; g = 241; b = 242;
            r = 255;
            g = 166;
            b = 98;
            break;
        case 'hydrophobic':
        case 2:
            //r = 0; g = 117; b = 220;
            r = 5;
            g = 200;
            b = 90;
            break;
        case 'van-der-waals':
        case 1:
            //r = 89; g = 252; b = 197;
            r = 100;
            g = 100;
            b = 100;
            break;
        default:
            r = 0;
            g = 0;
            b = 0;
    }

    if (rgb)
        return { r: r, g: g, b: b };
    else
        return rgb2hex(r, g, b);
}

function getFriendlyInteractionName(interaction) {
    /*switch (interaction) {
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
    }*/
    return interaction;
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
            r = 190;
            g = 190;
            b = 190;
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
            r = 230;
            g = 230;
            b = 230;
            //r = 150; g = 255; b = 150;
            break;
        case 'ECL1':
        case 'ECL2':
        case 'ECL3':
            r = 190;
            g = 190;
            b = 190;
            //r = 150; g = 150; b = 255;
            break;
        case 'ICL1':
        case 'ICL2':
        case 'ICL3':
            r = 190;
            g = 190;
            b = 190;
            //r = 150; g = 150; b = 255;
            break;
        default:
            r = 0;
            g = 0;
            b = 0;
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

// function downloadSVG2(svgSelector, name) {
//   var svgClone = $(svgSelector).clone();
//   svgClone.find('.svg-pan-zoom_viewport').attr('transform', 'matrix(2.2,0,0,2.2,295,140)');

//   var escapedSVG = new XMLSerializer().serializeToString(svgClone.get(0));

//   downloadURI('data:image/svg+xml;base64,' + window.btoa(escapedSVG), name);
// }

// function downloadSingleCrystalCSV(singleCrystalSvgSelector, name) {
//     var data = [];
//     var header = ['Residue number 1', 'Residue number 2', 'Segment 1', 'Segment 2',  'Generic number 1', 'Generic number 2', 'Amino acid 1', 'Amino acid 2', 'Interaction type'];
//     data.push(header);

//     $(singleCrystalSvgSelector + ' rect[data-interaction-type]').each(function(e) {
//       var rect = $(this);
//       var resNo1 = rect.data('res-no-1');
//       var resNo2 = rect.data('res-no-2');
//       var seg1 = rect.data('seg-1');
//       var seg2 = rect.data('seg-2');
//       var genNo1 = rect.data('gen-no-1');
//       var genNo2 = rect.data('gen-no-2');
//       var aa1 = rect.data('aa-1');
//       var aa2 = rect.data('aa-2');
//       var iType = rect.data('interaction-type');
//       data.push([resNo1, resNo2, seg1, seg2, genNo1, genNo2, aa1, aa2, iType]);
//     });

//     // Convert to CSV
//     var csv = Papa.unparse(data);

//     // Download file
//     downloadURI('data:text/csv;charset=UTF-8,' + encodeURI(csv), name);
// }

// function downloadSingleCrystalGroupCSV(singleGroupSvgSelector, name) {
//     var data = [];
//     var header = ['Generic number 1', 'Generic number 2', 'Segment 1', 'Segment 2', 'Frequency',  'Number of interactions', 'Number of crystals'];
//     data.push(header);

//     $(singleGroupSvgSelector + ' rect[data-frequency]').each(function(e) {
//       var rect = $(this);
//       var genNo1 = rect.data('gen-no-1');
//       var genNo2 = rect.data('gen-no-2');
//       var seg1 = rect.data('seg-1');
//       var seg2 = rect.data('seg-2');
//       var nInteractions = rect.data('num-interactions');
//       var nTotalInteractions = rect.data('total-possible-interactions');
//       var frequency = rect.data('frequency');
//       data.push([genNo1, genNo2, seg1, seg2, nInteractions, nTotalInteractions, frequency]);
//     });

//     // Convert to CSV
//     var csv = Papa.unparse(data);

//     // Download file
//     downloadURI('data:text/csv;charset=UTF-8,' + encodeURI(csv), name);
// }

// function downloadTwoCrystalGroupsCSV(twoGroupsSvgSelector, name) {
//     var data = [];
//     var header = ['Generic number 1', 'Generic number 2', 'Segment 1', 'Segment 2', 'Interactions group 1', 'Interactions group 2', 'Crystals group 1', 'Crystals group 2', 'Frequency group 1', 'Frequency group 2', 'Frequency Difference'];
//     data.push(header);

//     $(twoGroupsSvgSelector + ' rect[data-frequency-diff]').each(function(e) {
//       var rect = $(this);
//       var genNo1 = rect.data('gen-no-1');
//       var genNo2 = rect.data('gen-no-2');
//       var seg1 = rect.data('seg-1');
//       var seg2 = rect.data('seg-2');
//       var numIntsGroup1 = rect.data('group-1-num-ints');
//       var numIntsGroup2 = rect.data('group-2-num-ints');
//       var numPdbsGroup1 = rect.data('group-1-num-pdbs');
//       var numPdbsGroup2 = rect.data('group-2-num-pdbs');
//       var freqGroup1 = rect.data('group-1-freq');
//       var freqGroup2 = rect.data('group-2-freq');
//       var fDiff = rect.data('frequency-diff').toFixed(2);
//       data.push([genNo1, genNo2, seg1, seg2, numIntsGroup1, numIntsGroup2, numPdbsGroup1, numPdbsGroup2, freqGroup1, freqGroup2, fDiff]);
//     });

//     // Convert to CSV
//     var csv = Papa.unparse(data);

//     // Download file
//     downloadURI('data:text/csv;charset=UTF-8,' + encodeURI(csv), name);
// }

// function downloadURI(uri, name) {
//     var link = document.createElement("a");
//     link.download = name;
//     link.href = uri;
//     document.body.appendChild(link);
//     link.click();
//     document.body.removeChild(link);
//     delete link;
// }

function get_current_mode() {
    return $(".main_option:visible").attr("id").replace("-tab", "");
}


function redraw_renders() {

    // Makes sure diagrams fit sizes
    console.time('redraw renders');
    var visible_svg = $('svg:visible');
    var svg_class = visible_svg.attr("class");
    if (window.innerHeight == screen.height || is_fullscreen) {
        // browser is fullscreen
        console.log('dont redraw in full screen', is_fullscreen);

        if (svg_class == 'flareplot') {
            $("svg.flareplot").css('height', screen.height);
        } else {
            visible_svg.css('height', screen.height);
        }

        console.timeEnd('redraw renders');
        return
    }


    $('div.dataTables_scrollBody:visible').height('50vh');

    // Make sure browser-tables are not too wide.
    browser_table_div_width = $('.contact-browser:visible').width();
    if (browser_table_div_width > 2060) {
        browser_table_width = 2030;
    } else {
        browser_table_width = browser_table_div_width - 30;
    }
    $('.contact-browser .dataTables_wrapper').width(browser_table_width);
    // console.log(browser_table_div_width,browser_table_width);

    var width_svg = visible_svg.width();

    // Don't go too high
    if (width_svg > screen.height) width_svg = screen.height;
    width_svg = 500;
    visible_svg.height(width_svg);

    // Resize NGL
    var ngl = $('.ngl-container:visible');
    var width_ngl = ngl.width();

    // Don't go too high
    if (width_ngl > screen.height) width_ngl = screen.height;
    ngl.height(width_svg);

    // console.log('redraw',svg_class,width_svg,screen.height);

    if (svg_class == 'heatmap') {
        // If heatmap being resized, reset zoom

        // Destroy old zoom on heatmap
        var heatMapSelector = "#" + $('.main_option:visible').attr('id');
        heatMapSelector = heatMapSelector + ' .heatmap-container';
        if (window.zoomHeatmap[heatMapSelector] != null) {
            window.zoomHeatmap[heatMapSelector].destroy();
            delete window.zoomHeatmap[heatMapSelector];

            window.zoomHeatmap[heatMapSelector] = svgPanZoom(heatMapSelector + ' .heatmap', {
                zoomEnabled: true,
                // controlIconsEnabled: true,
                fit: true,
                center: true,
                minZoom: 0.40,
                maxZoom: 50,
                zoomScaleSensitivity: 0.25,
                dblClickZoomEnabled: true,
                beforeZoom: hidePopovers,
                beforePan: hidePopovers
            });

            window.zoomHeatmap[heatMapSelector].zoom(0.85);
        }
    }

    if (svg_class == 'hiveplot') {
        // If heatmap being resized, reset zoom

        // Destroy old zoom on heatmap
        var container = "#" + $('.main_option:visible').attr('id');
        container = container + ' .hiveplot-container';
        if (window.zoomHiveplot[container] != null) {
            window.zoomHiveplot[container].destroy();
            delete window.zoomHiveplot[container];

            window.zoomHiveplot[container] = svgPanZoom(container + ' .hiveplot', {
                zoomEnabled: true,
                // controlIconsEnabled: true,
                fit: true,
                center: true,
                minZoom: 0.75,
                maxZoom: 50,
                zoomScaleSensitivity: 0.25,
                dblClickZoomEnabled: true
            });

            window.zoomHiveplot[container].zoom(0.85);
        }
    }
    console.timeEnd('redraw renders');

}

function initializeSegmentButtons(selector) {
    // Initialize segment buttons.
    $(selector + ' .segments-panel button').each(function () {
        var s = $(this).attr('data-segment');

        // Return if no segment data
        if (s == null) {
            return;
        }

        $(this).click(function () {
            var segments = [];
            $(this).toggleClass('active');
            $(selector + ' .segments-panel button.active').each(function () {
                segments = segments.concat($(this).data('segment').split(' '));
            });
            $(selector + ' .segments-input').val(JSON.stringify(segments));
        });
    });

    // Initialize 'all' buttons.
    $(selector + ' .segments-panel .all-button').each(function () {
        $(this).click(function () {
            if ($(this).html() === 'All') {
                $(this).html('None');
                $(this).parent().find('button').each(function () {
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
                $(this).parent().find('button').each(function () {
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
        $(selector + ' .segments-panel button.active').each(function () {
            segments.push($(this).data('segment'));
        });
        $(selector + ' .segments-panel .segments-input').val(JSON.stringify(segments));

        // Trigger click on initialization
        $(this).trigger('click');
    });
}

function initializeInteractionButtons(selector) {
    // Initialize interaction buttons.
    $(selector + ' .interactions-panel button').each(function () {
        var s = $(this).attr('data-interaction-type');

        // Return if no segment data
        if (s == null) {
            return;
        }

        $(this).click(function () {
            var interactions = [];
            $(this).toggleClass('active');
            $(selector + ' .interactions-panel button.active').each(function () {
                interactions = interactions.concat($(this).data('interaction-type').split(' '));
            });
            $(selector + ' .interactions-input').val(JSON.stringify(interactions));
        });
    });

    // Initialize 'all' buttons.
    $(selector + ' .interactions-panel .all-button').each(function () {
        $(this).click(function () {
            if ($(this).html() === 'All') {
                $(this).html('None');
                $(this).parent().find('button').each(function () {
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
                $(this).parent().find('button').each(function () {
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
        $(selector + ' .interactions-panel button.active').each(function () {
            interactions.push($(this).data('interaction-type'));
        });
        $(selector + ' .interactions-input').val(JSON.stringify(interactions));

        // Trigger click on initialization
        $(this).trigger('click');
    });
}


function initializeGoButton(selector, generic = false) {
    $(selector + ' .go-button').click(function () {
        var pdb = JSON.parse($(selector + ' .crystal-pdb').val());
        loadPDBsView(pdb, selector, generic)
    });
}

function drawPlotPanel(plot_type, plot_div) {
    plot_id = plot_div.attr('id');
    // console.log(plot_type);
    // console.log(plot_div);
    // console.log(plot_id);

    // Delete whatever is already there
    plot_div.find('.plot-container').html('');
    plot_div.find('.plot-container').attr('class', 'plot-container');
    var mode = get_current_mode();

    plot_div.find('.plot-title').html("&nbsp;" + display_plot_names[plot_type]);

    ngl_color_mode = false;

    if (plot_type.startsWith("ngl") && plot_type.includes("_")) {
        ngl_color_mode = plot_type.split("_")[1];
        plot_type = "ngl";
    }

    console.log("SET UP PLOT", plot_type, plot_id, mode,'ngl_mode',ngl_color_mode);
    switch (mode) {
        case "two-crystal-groups":
            raw_data = two_sets_data;
            break;
        case "single-crystal":
            raw_data = single_crystal_data;
            break;
        case "single-crystal-group":
            raw_data = single_set_data;
            break;
    }
    switch (plot_type) {
        case "ngl":
            plot_div.find('.plot-container').removeClass('none');
            plot_div.find('.plot-container').addClass('ngl-container');
            plot_div.find('.plot-container').attr('id', 'ngl-' + plot_id);
            plot_div.find('.plot-container').attr('style', 'margin: auto; width: 100%; overflow: hidden;');

            // TODO add method to select PDB with highest number of GNs (in set)

            switch (mode) {
                case "two-crystal-groups":
                    createNGLview(plot_id, two_sets_pdbs1[0], two_sets_pdbs1, two_sets_pdbs2, two_sets_pdbs2[0]);
                    break;
                case "single-crystal":
                    createNGLview(plot_id, single_crystal_pdb);
                    break;
                case "single-crystal-group":
                    createNGLview(plot_id, single_set_pdbs[0], single_set_pdbs);
                    break;
            }
            break;
        case "heatmap":
            plot_div.find('.plot-container').removeClass('none');
            plot_div.find('.plot-container').addClass('heatmap-container');
            plot_div.find('.plot-container').attr('id', "heatmapcontainer-" + plot_id);
            plot_div.find('.plot-container').html('<svg class="heatmap" xmlns="http://www.w3.org/2000/svg" preserveAspectRatio="xMinYMin" id="heatmap-' + plot_id + '" style="height: 500px;"></svg>');

            renderHeatmap(raw_data, '#heatmapcontainer-' + plot_id);
            break;
        case "heatmap_distances":
            plot_div.find('.plot-container').removeClass('none');
            plot_div.find('.plot-container').addClass('heatmap_distances-container');
            plot_div.find('.plot-container').attr('id', "heatmap_distancescontainer-" + plot_id);
            plot_div.find('.plot-container').html('<svg class="heatmap_distances" xmlns="http://www.w3.org/2000/svg" preserveAspectRatio="xMinYMin" id="heatmap_distances-' + plot_id + '" style="height: 500px;"></svg>');

            renderHeatmap_distances(raw_data, '#heatmap_distancescontainer-' + plot_id);
            break;
        case "flareplot":
            plot_div.find('.plot-container').removeClass('none');
            plot_div.find('.plot-container').addClass('flareplot-container');
            plot_div.find('.plot-container').attr('id', 'flareplot-' + plot_id);

            createFlareplotBox(raw_data, '#flareplot-' + plot_id, toggle = false);
            break;
        case "flareplot_subset":
            plot_div.find('.plot-container').removeClass('none');
            plot_div.find('.plot-container').addClass('flareplot-container-subset');
            plot_div.find('.plot-container').attr('id', 'flareplot-' + plot_id);

            createFlareplotBox(raw_data, '#flareplot-' + plot_id, toggle = false, subset = true);
            break;
        case "flareplot_segments":
            plot_div.find('.plot-container').removeClass('none');
            plot_div.find('.plot-container').addClass('flareplot-container-segments');
            plot_div.find('.plot-container').attr('id', 'flareplot-' + plot_id);

            createFlareplot_segment(raw_data, 1000, '', '#flareplot-' + plot_id);
            break;
        case "force_network":
            plot_div.find('.plot-container').removeClass('none');
            plot_div.find('.plot-container').addClass('networkPlot');
            plot_div.find('.plot-container').attr('id', 'network-' + plot_id);

            createNetworkPlot(raw_data, 400, '', '#network-' + plot_id, segment_view = false);
            break;
        case "force_network_segment":
            plot_div.find('.plot-container').removeClass('none');
            plot_div.find('.plot-container').addClass('networkPlot');
            plot_div.find('.plot-container').attr('id', 'network-' + plot_id);

            createNetworkPlot(raw_data, 400, '', '#network-' + plot_id, segment_view = true);
            break;
        case "force_network_3d":
            plot_div.find('.plot-container').removeClass('none');
            plot_div.find('.plot-container').addClass('networkPlot');
            plot_div.find('.plot-container').attr('id', 'network-' + plot_id);

            createNetworkPlot3D(raw_data, 400, '', '#network-' + plot_id, segment_view = false);
            break;
        case "force_network_3d_segment":
            plot_div.find('.plot-container').removeClass('none');
            plot_div.find('.plot-container').addClass('networkPlot');
            plot_div.find('.plot-container').attr('id', 'network-' + plot_id);

            createNetworkPlot3D(raw_data, 400, '', '#network-' + plot_id, segment_view = true);
            break;
        case "boxplot":
            plot_div.find('.plot-container').removeClass('none');
            plot_div.find('.plot-container').addClass('boxplot-container');
            plot_div.find('.plot-container').attr('id', 'boxplot-' + plot_id);

            createBoxPlot(raw_data, 'boxplot-' + plot_id);
            break;
        case "boxplot_angles":
            plot_div.find('.plot-container').removeClass('none');
            plot_div.find('.plot-container').addClass('boxplot-container');
            plot_div.find('.plot-container').attr('id', 'boxplot-' + plot_id);

            createBoxPlot(raw_data, 'boxplot-' + plot_id, 'angles');
            break;
        case "snakeplot":
            plot_div.find('.plot-container').removeClass('none');
            plot_div.find('.plot-container').addClass('snakeplot-container');
            plot_div.find('.plot-container').attr('id', 'snakeplot-' + plot_id);

            createSnakeplot(raw_data, '#snakeplot-' + plot_id);
            break;
        case "scatterplot":
            plot_div.find('.plot-container').removeClass('none');
            plot_div.find('.plot-container').addClass('scatterplot-container');
            plot_div.find('.plot-container').attr('id', 'scatterplot-' + plot_id);

            createScatterplot(raw_data, 'scatterplot-' + plot_id);
            break;
        case "tm7_plot_major":
            plot_div.find('.plot-container').removeClass('none');
            plot_div.find('.plot-container').addClass('tm_movment-container');
            plot_div.find('.plot-container').attr('id', 'tm_movment-' + plot_id);
            tm7_plot('#tm_movment-' + plot_id, raw_data["tm_movement_2D"]["classA_ligands"],raw_data["tm_movement_2D"]["viewbox_size"]);
            break;
        case "tm7_plot_middle":
            plot_div.find('.plot-container').removeClass('none');
            plot_div.find('.plot-container').addClass('tm_movment-container');
            plot_div.find('.plot-container').attr('id', 'tm_movment-' + plot_id);
            tm7_plot('#tm_movment-' + plot_id, raw_data["tm_movement_2D"]["membrane_mid"],raw_data["tm_movement_2D"]["viewbox_size"]);
            break;
        case "tm7_plot_intra":
            plot_div.find('.plot-container').removeClass('none');
            plot_div.find('.plot-container').addClass('tm_movment-container');
            plot_div.find('.plot-container').attr('id', 'tm_movment-' + plot_id);
            tm7_plot('#tm_movment-' + plot_id, raw_data["tm_movement_2D"]["intracellular"],raw_data["tm_movement_2D"]["viewbox_size"]);
            break;
        case "tm7_plot_extra":
            plot_div.find('.plot-container').removeClass('none');
            plot_div.find('.plot-container').addClass('tm_movment-container');
            plot_div.find('.plot-container').attr('id', 'tm_movment-' + plot_id);
            tm7_plot('#tm_movment-' + plot_id, raw_data["tm_movement_2D"]["extracellular"],raw_data["tm_movement_2D"]["viewbox_size"]);
            break;
        case "tm7_plot_3d_major":
            plot_div.find('.plot-container').removeClass('none');
            plot_div.find('.plot-container').addClass('tm_movment-container');
            plot_div.find('.plot-container').attr('id', 'tm_movment-' + plot_id);
            tm7_plot_3d('#tm_movment-' + plot_id, raw_data["tm_movement_2D"]["classA_ligands"]);
            break;
        case "tm7_plot_3d_middle":
            plot_div.find('.plot-container').removeClass('none');
            plot_div.find('.plot-container').addClass('tm_movment-container');
            plot_div.find('.plot-container').attr('id', 'tm_movment-' + plot_id);
            tm7_plot_3d('#tm_movment-' + plot_id, raw_data["tm_movement_2D"]["membrane_mid"]);
            break;
        case "tm7_plot_3d_intra":
            plot_div.find('.plot-container').removeClass('none');
            plot_div.find('.plot-container').addClass('tm_movment-container');
            plot_div.find('.plot-container').attr('id', 'tm_movment-' + plot_id);
            tm7_plot_3d('#tm_movment-' + plot_id, raw_data["tm_movement_2D"]["intracellular"]);
            break;
        case "tm7_plot_3d_extra":
            plot_div.find('.plot-container').removeClass('none');
            plot_div.find('.plot-container').addClass('tm_movment-container');
            plot_div.find('.plot-container').attr('id', 'tm_movment-' + plot_id);
            tm7_plot_3d('#tm_movment-' + plot_id, raw_data["tm_movement_2D"]["extracellular"]);
            break;
    }

}

var plotting_options = {
    'TM1-7 segment movement' : {
        'cytosolic': [
            //['tm7_heatmap_intra','Heatmap'],
            ['tm7_plot_intra', 'Segment plot (2D)'],
            ['tm7_plot_3d_intra', 'Segment plot (3D)']
        ],
        'extracellular': [
            //['tm7_heatmap_extra','Heatmap'],
            ['tm7_plot_extra', 'Segment plot (2D)'],
            ['tm7_plot_3d_extra','Segment plot (3D)']
        ],
        // 'class A major pocket': [
        //     ['tm7_plot_major', '2D plot'],
        //     ['tm7_plot_3d_major','3D plot'],
        //     ['tm7_heatmap_major','Heatmap']
        // ],
        'Middle of membrane': [
            //['tm7_heatmap_middle','Heatmap'],
            ['tm7_plot_middle', 'Segment plot (2D)'],
            ['tm7_plot_3d_middle','Segment plot (3D)']
        ],
    },
    'Contacts between generic residue positions': [
        ['flareplot', 'Flare Plot (all contacts)'],
        ['flareplot_subset', 'Flare Plot (shown contacts)'],
        ['heatmap', 'Heatmap'],
        ['force_network', 'Network (2D)'],
        ['force_network_3d', 'Network (3D)'],
        ['ngl', 'Structure (3D)'],
        // ['schematic_non', 'Schematic (Non-consecutive)'],
        // ['schematic_con', 'Schematic (Consecutive)'],
    ],
    'Contacts between segments (TM1-7, H8 & loops)': [
        ['flareplot_segments', 'Flare Plot'],
        ['force_network_segment', 'Network (2D)'],
        ['force_network_3d_segment', 'Network (3D)'],
    ],
    'Contact frequencies': [
        ['boxplot', 'Box plot'],
    ],
    'Residue Properties': [
        ['boxplot_angles', 'Box plot (distribution)'],
        ['heatmap_distances', 'Heatmap (distance)'],
        ['scatterplot', 'Scatter plot (correlation)'],
        ['snakeplot', 'Snakeplot (2D, topology)'],
        ['ngl_distances', 'Structure (3D, movement)']],
    // '3D structure': { '3D structures': [['ngl_distances', 'Distances'], ['ngl_angles', 'Angles']] },
};

display_plot_names = {}



function generate_display_options() {
    dropdown_html = '<div class="dropdown" style="display: inline;"> \
                          <button class="btn btn-xs btn-primary dropdown-toggle" type="button" data-toggle="dropdown"> \
                          Select plot \
                          <span class="caret"></span></button><ul class="dropdown-menu">';
    var after_first = false;
    for (let key in plotting_options) {

        if (after_first) dropdown_html += '<li class="divider"></li>';
        if ($.isArray(plotting_options[key])) {
            dropdown_html += '<li class="dropdown-header text-uppercase"><strong>' + key + '</strong></li>'
            plotting_options[key].forEach(function (opt) {
                display_plot_names[opt[0]] = '<span class="text-uppercase"><strong>'+key + '</strong></span> - ' + opt[1];
                dropdown_html += '<li><a class="plot_selection" href="#" plot_type="' + opt[0] + '">' + opt[1] + '</a></li>'
            });
        } else {
            dropdown_html += '<li class="dropdown-submenu dropleft"><a tabindex="0" href="#">' + key + '</a>'
            dropdown_html += '<ul class="dropdown-menu">'
            for (let key2 in plotting_options[key]) {
                dropdown_html += '<li class="dropdown-header text-uppercase"><strong>' + key2 + '</strong></li>'
                plotting_options[key][key2].forEach(function (opt) {
                    display_plot_names[opt[0]] = '<span class="text-uppercase"><strong>'+key + '</strong></span> - <span class="text-uppercase"><strong>'+ key2 + '</strong></span> - ' + opt[1];
                    dropdown_html += '<li><a class="plot_selection" href="#" plot_type="' + opt[0] + '">' + opt[1] + '</a></li>'
                });
            }
            dropdown_html += '</ul></li>'
        }
        after_first = true;
    }
    dropdown_html += '</ul></div>';
    $('.plot-select').each(function (e) {
        $(this).html(dropdown_html);
    });
    $('.plot_selection').click(function () {
        plot_type = $(this).attr('plot_type');
        plot_div = $(this).closest('.panel');

        drawPlotPanel(plot_type, plot_div);
    });
}

var single_set_pdbs = '';
var single_crystal_pdb = '';
var single_set_data = '';
var single_crystal_data = '';
function loadPDBsView(pdb, selector, generic) {
    console.time('Get loadPDBsView Data');
    $(".main_loading_overlay").show();
    // $(".main_loading_overlay").show();
    //var segments = JSON.parse($(selector + ' .segments-input').val());
    //var segments = ['TM1', 'TM2', 'TM3', 'TM4', 'TM5', 'TM6', 'TM7', 'H8', 'ICL1', 'ECL1', 'ICL2', 'ECL2', 'ICL3', 'ECL3', 'N-term', 'C-term'];
    if (pdb.length > 0 /*&& segments.length > 0*/) {
        var interactionTypes = JSON.parse($(selector + ' .interactions-input').val());
        $(".heatmap").hide();
        // $(".heatmap-legend").hide();
        $(".matrix-tab:visible").click();

        $(selector + ' .heatmap-container').append('<span id=svgloading>Loading...</span>');
        //                if (!$(selector + ' .interactions-input').val() == null)
        //                    interactionTypes = JSON.parse($(selector + ' .interactions-input').val());

        //normalized = $(".normalized:visible input").prop("checked");
        $.ajax({
            url: '/contactnetwork/browserdata',
            dataType: 'json',
            data: {
                // 'segments': segments,
                'csrfmiddlewaretoken': csrf_token,
                'generic': generic,
                'pdbs': pdb,
                //'normalized': normalized,
                'interaction_types': currentSettings[currentTab]["types"],
                'strict_interactions': currentSettings[currentTab]["strict"],
                'options': currentSettings[currentTab]["options"]
            },
            async: true,
            method: "POST",
            success: function (data) {
                console.timeEnd('Get loadPDBsView Data');
                if (data['error']) {
                    alert('ERROR, multiple classes selected (' + data['error'] + ') - Please change your selection and retry');
                    $(".main_loading_overlay").hide();
                    return
                }
                // Re-render heatmap
                data_browser = data;

                var mode = get_current_mode();
                switch (mode) {
                    case "single-crystal-group":
                        single_set_data = data;
                        single_set_pdbs = pdb;
                        break;
                    case "single-crystal":
                        single_crystal_data = data;
                        single_crystal_pdb = pdb;
                        break;
                }
                renderBrowser(data);
                renderBrowser_2(data);
                //renderBrowser_3(data);
                renderBrowser_4(data);
                renderBrowser_5(data);
                generate_display_options();
                browser_visible = $(".nav-browsers:visible li.active a").attr('id');
                renderDataTablesYadcf(browser_visible);
                $(".main_loading_overlay").hide();

                // Set up default visualisation
                initilizeInitialPlots();
            }
        });

    }
}


function initializeGoButtonTwoCrystalGroups(selector, generic = false) {
    $(selector + ' .go-button').click(function () {
        var pdbs1 = JSON.parse($(selector + ' .crystal-group-1-pdbs').val());
        var pdbs2 = JSON.parse($(selector + ' .crystal-group-2-pdbs').val());
        loadTwoPDBsView(pdbs1, pdbs2, selector, generic)


    });
}

var two_sets_pdbs1 = '';
var two_sets_pdbs2 = '';
var two_sets_data = '';
function loadTwoPDBsView(pdbs1, pdbs2, selector, generic) {
    console.time('Get loadTwoPDBsView Data');
    $(".main_loading_overlay").show();
    //var segments = JSON.parse($(selector + ' .segments-input').val());
    //var segments = ['TM1', 'TM2', 'TM3', 'TM4', 'TM5', 'TM6', 'TM7', 'H8', 'ICL1', 'ECL1', 'ICL2', 'ECL2', 'ICL3', 'ECL3', 'N-term', 'C-term'];
    if (pdbs1.length > 0 && pdbs2.length > 0 /*&& segments.length > 0*/) {
        var interactionTypes = JSON.parse($(selector + ' .interactions-input').val());
        $(".heatmap").hide();
        // $(".heatmap-legend").hide();
        $(".matrix-tab:visible").click();
        $(selector + ' .heatmap-container').append('<span id=svgloading>Loading... (0%)</span>');

        two_sets_pdbs1 = pdbs1;
        two_sets_pdbs2 = pdbs2;

        //normalized = $(".normalized:visible input").prop("checked");

        $.ajax({
            url: '/contactnetwork/browserdata',
            dataType: 'json',
            data: {
                // 'segments': segments,
                'csrfmiddlewaretoken': csrf_token,
                'generic': generic,
                'pdbs1': pdbs1,
                'pdbs2': pdbs2,
                //'normalized': normalized,
                'interaction_types': currentSettings[currentTab]["types"],
                'strict_interactions': currentSettings[currentTab]["strict"],
                'options': currentSettings[currentTab]["options"]
            },
            async: true,
            method: "POST",
            success: function (data) {
                console.timeEnd('Get loadTwoPDBsView Data');
                if (data['error']) {
                    alert('ERROR, multiple classes selected (' + data['error'] + ') - Please change your selection and retry');
                    $(".main_loading_overlay").hide();
                    return
                }
                // Re-render heatmap
                two_sets_data = data;
                renderBrowser(data);
                renderBrowser_2(data);
                // renderBrowser_3(data); //Disabled for now..
                renderBrowser_4(data);
                renderBrowser_5(data);
                browser_visible = $(".nav-browsers:visible li.active a").attr('id');
                renderDataTablesYadcf(browser_visible);
                generate_display_options();
                $(".main_loading_overlay").hide();

                // Set up default visualisation
                initilizeInitialPlots();
            }
        });
    }
    // $(".main_loading_overlay").hide();
}

function initilizeInitialPlots() {
    default_plot_types = ['force_network', 'flareplot', 'ngl'];
    var mode = get_current_mode();
    // if single structure - use interaction coloring
    if (mode == "two-crystal-groups") {
        default_plot_types = ['tm7_plot_extra', 'tm7_plot_middle', 'tm7_plot_intra'];
        // default_plot_types = ['scatterplot', 'snakeplot', ''];
    }

    $(".plot_row:visible").find(".panel").each(function (i) {
        plot_type = default_plot_types[i];
        plot_div = $(this);
        drawPlotPanel(plot_type, plot_div);
    })
}


function initializeFullscreenButton(selector) {
    // var fullScreenElement = $(selector + ' .heatmap-container').get(0);
    $(selector + ' .btn-fullscreen').click(function () {
        is_fullscreen = true;
        // console.log($(this).parent().parent().next().children().first());
        //                console.log($(this).attr('id'));
        fullScreenElement = $(this).parent().parent().next().children().first().find("canvas");
        if (fullScreenElement.length) {
            fullScreenElement.css('background-color', 'white');
            console.log('who to fullscreen?', fullScreenElement.closest('.plot-container').attr('id'));
            plot_id = fullScreenElement.closest('.plot-container').attr('id');
            ngl_stage = plot_id.replace("ngl-", "");
            if (ngl_stage in stage) {
                $(function () {
                    stage[ngl_stage].toggleFullscreen()();
                });
                return
            }
        } else {
            fullScreenElement = $(this).closest(".panel-default").find(".plot-container");
            fullScreenElement.css('background-color', 'white');
            var cp = fullScreenElement.find(".controls-panel");
            cp.toggleClass("fullscreen");
        }

        toggleFullScreen(fullScreenElement.get(0));

        if (fullScreenElement.attr('id')) {
            if (fullScreenElement.attr('id').startsWith('DataTable')) {
                top_height = $('div.dataTables_scrollHead:visible').outerHeight();
                bottom_height = $('div.dataTables_info:visible').outerHeight();
                scrollbody_height = screen.height - top_height - bottom_height;
                $('div.dataTables_scrollBody:visible').height(scrollbody_height + 'px');
            }
        }
    });
}


function createTwoGroupFlareplotBox(data1, data2, data3, container, toggle = false) {
    // prepare two group data for visualization
    var data = data3;

    // frequency + count holder
    data["frequency"] = {};
    data["count"] = {};

    Object.keys(data.interactions).forEach(function (pair) {
        var f1 = 0,
            f2 = 0,
            c1 = 0,
            c2 = 0;;
        if (pair in data1.interactions) {
            c1 = Object.keys(data1.interactions[pair]).length;
            f1 = c1 / data1.pdbs.length;
        }
        if (pair in data2.interactions) {
            c2 = Object.keys(data2.interactions[pair]).length;
            f2 = c2 / data2.pdbs.length;
        }
        var f3 = f1 - f2;
        var c3 = c1 + c2;
        data["frequency"][pair] = [f1, f2, f3];
        data["count"][pair] = [c1, c2, c3];
    });

    createFlareplotBox(data, container, toggle = false);
}

var flareplot = {};
var contiguous = true;
var interactionsToggleList = [];

function createFlareplotBox(data, container, toggle = false, subset = false) {
    // clean
    if (toggle) {
        $(container).children().last().remove();
    } else {
        // in case refresh with new parameters => reset
        contiguous = true;
        $(container).html("");
    }

    // add menu
    if (!toggle) {
        var newDiv = document.createElement("div");
        newDiv.setAttribute("class", "flareplot-legend");

        var content = '<div class="controls">'
        //                                  +'<h4>Controls</h4>';

        // only possible with more than 4 segments, otherwise it will become a mess
        if (data.segments.length > 4 && !subset)
            content += '<p>Consecutive segment contacts on outside: <input type=checkbox id="flareplot_contiguous" checked></p>';

        content += '<p>Line colors: <select id="flareplot_color">' +
            '<option value="none">None (gray)</option>' +
            '<option value="rainbow">GPCR rainbow</option>' +
            '<option value="segment">GPCR segment</option>';

        var mode = get_current_mode();
        // if single structure - use interaction coloring
        if (mode == "single-crystal") {
            content += '<option value="interactions" selected>Interaction Type</option>';
            // if single group of structures - use frequency coloring (gradient)
        } else if (mode == "single-crystal-group") {
            content += '<option value="frequency" selected>Interaction Frequency/Count</option>';
            // if group(s) of structures - use frequency coloring (gradient)
        } else {
            content += '<option value="frequency" selected>Frequency difference Gr1 - Gr2</option>';
            content += '<option value="frequency_1">Frequency group 1</option>';
            content += '<option value="frequency_2">Frequency group 2</option>';
        }
        content += '</select></p>';


        // Populate heatmap legend
        if (container.indexOf("single-crystal-tab") >= 0) {
            /*
            // TODO Optimize/generalize interaction set selection - code duplication
            var interactionTypes = new Set(["Ionic", "Polar", "Aromatic", "Hydrophobic", "Van-der-Waals"]);

            // Add interactions color legend
            //content += '<h4>Toggle interactions</h4><ul>';
            content += '<ul>';

            interactionTypes = Array.from(interactionTypes).sort(function (i1, i2) {
                return getInteractionStrength(i2) - getInteractionStrength(i1);
            });

            interactionTypes.forEach(function(i) {
                var rgb = getInteractionColor(i, false);
                content += '<li>'
                        + '<div class="color-box" style="background-color: ' + rgb + '">' + '<input type="checkbox" data-interaction-type="' + i.toLowerCase() +'"></input>' + '</div><p>' + i + '</p>'
                        + '</li>';
            });
            content += '</ul></div>';*/
        } else if (container.indexOf("single-crystal-group-tab") >= 0) {
            // slider 1 to count
            /*content += '<h4 class="center">Frequency (#PDBs)</h4>'
                + '<p>Range: <span id="pdbs-range-flare">1 - ' + data.pdbs.length + '</span></p>'
                + '<div class="slider-range" data-text-id="pdbs-range-flare" id="pdbs-range-flare-slider"></div>'
                + '<div class="temperature-scale">'
                + '<span class="gray-to-red"></span>'
                + '</div>';

            $( function() {
              $( container+" .slider-range" ).data({ "referenceContainer" : container })
              $( container+" .slider-range" ).slider({
                range: true,
                min: 1,
                max: data.pdbs.length,
                step: 1,
                values: [0,data.pdbs.length],
                slide: function( event, ui ) {
                  $( "#"+$(this).attr("data-text-id") ).html( ui.values[ 0 ] + " - " + ui.values[ 1 ] );
                  flareplot[$(this).data("referenceContainer")].updateRange(ui.values[ 0 ], ui.values[ 1 ]);
                }
              });
            } );*/
        } else {
            /*content += '<h4 class="center">Frequency</h4>'
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
                });*/
        }

        newDiv.innerHTML = content;

        $(container).append(newDiv);

        $(container + " #flareplot_contiguous").click(function (e) {
            $(function () {
                contiguous = !contiguous;
                createFlareplotBox(data, container, toggle = true);
                redraw_renders();
            });
        });


        $(container + " #flareplot_color").data({ "referenceContainer": container })
        $(container + " #flareplot_color").change(function (e) {
            flareplot[$(this).data("referenceContainer")].updateColors($(this).val(), interactionsToggleList);
        });

        $(container + ' .flareplot-legend .color-box input[type=checkbox]').each(function () {
            $(this).prop('checked', true);

            // init toggle list
            interactionsToggleList.push($(this).data('interaction-type'));
            $(this).data({ "referenceContainer": container })

            $(this).change(function () {
                // toggle interactions in flareplot
                interactionsToggleList = [];

                $(container + ' .flareplot-legend input[type=checkbox]').each(function () {
                    if ($(this).prop('checked'))
                        interactionsToggleList.push($(this).data('interaction-type'));
                });

                // toggle interactions in flareplot
                flareplot[$(this).data("referenceContainer")].showInteractions(interactionsToggleList);
            });
        });
    }

    // create flareplot
    if (subset) {
        flareplot[container] = createFlareplot_subset(300, parseGPCRdb2flare(data, subset = true), container, false);
        $(container + " #flareplot_color").val("rainbow");
    } else {
        flareplot[container] = createFlareplot(1000, parseGPCRdb2flare(data), container, contiguous);
    }

    // update coloring and visibility if toggled
    if (toggle) {
        // var range = $( container+" .slider-range" ).slider("values");
        // flareplot[container].updateRange(range[0], range[1]);

        flareplot[container].updateColors($(container + " #flareplot_color").val(), interactionsToggleList);
        //if ($(container + ' .flareplot-legend .color-box input[type=checkbox]').length > 0)
        //  flareplot[container].showInteractions(interactionsToggleList);
        // if (container.indexOf("two-crystal-groups-tab") >= 0)
        //   updateTwoGroupSliders("skip", []);

        updateGeneralControls()
    } else {
        flareplot[container].updateColors($(container + " #flareplot_color").val(), interactionsToggleList);
    }
}

function createHiveplotBox(data, container) {
    // clean contents
    $(container).html("");

    createHiveplot(data, container);

    // display in the back to enable SVGPan
    $("#single-hiveplot-tab").show();

    // Make zoomable
    window.zoomHiveplot[container] = svgPanZoom(container + ' .hiveplot', {
        zoomEnabled: true,
        // controlIconsEnabled: true,
        fit: true,
        center: true,
        minZoom: 0.75,
        maxZoom: 50,
        zoomScaleSensitivity: 0.25,
        dblClickZoomEnabled: true
    });

    // remove display value on element
    $("#single-hiveplot-tab").css("display", "")
}


// TODO: make this function obsolete and merge remaining code with *createNGLRepresentations*
function updateStructureRepresentations(mode) {
    console.log('updateStructureRepresentations', mode);
    var structures = 1;
    if (mode.startsWith("two_sets"))
        structures = 2;

    for (var key = 0; key < structures; key++) {

        hide_structure = $("#ngl-" + mode + " #hide_pdb" + (key + 1)).prop('checked');
        var o = reps[mode][key].structureComponent;
        if (hide_structure) {
            o.setVisibility(false);
            break;
        } else {
            o.setVisibility(true);
        }
        // toggle edges
        reps[mode][key].links.setVisibility($("#ngl-" + mode + " #toggle_interactions").prop('checked'));

        // toggle CA spheres
        reps[mode][key].int_res.setVisibility($("#ngl-" + mode + " #highlight_res").prop('checked'));

        // toggle interacting toggle_sidechains
        reps[mode][key].ball_int.setVisibility($("#ngl-" + mode + " #toggle_sidechains_int").prop('checked'));

        // toggle segment movement spheres
        if (mode_short == 'two-groups')
          reps[mode][key].movement_spheres.setVisibility($("#ngl-" + mode + " #toggle_movement").prop('checked'));

        // Update cartoon using selection
        checked = $("#ngl-" + mode + " #ngl_only_gns").prop('checked');
        sele = ":" + pdb_data[mode][key]['chain'];
        int_sele = sele
        if (checked)
            sele = ":" + pdb_data[mode][key]['chain'] + " and (" + pdb_data[mode][key]['only_gn'].join(", ") + ")";
        int_sele =

            gpcr_rep[mode][key].setSelection(sele);
    }
}


function redraw_ngl() {
    var mode = $('ul#mode_nav').find('li.active').find('a').text().trim();
    var mode = $('.ngl-container:visible').attr('id').replace('ngl-', '');
    if (mode in stage) {
        $(function () {
            stage[mode].handleResize();
        });
    }
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
    $.get('/contactnetwork/pdbtabledata', function (data) {
        $('#single-crystal-pdb-modal-table .tableview').html(data);
        $('#single-crystal-group-pdbs-modal-table .tableview').html(data);
        $('#two-crystal-group-pdbs-modal-1-table .tableview').html(data);
        $('#two-crystal-group-pdbs-modal-2-table .tableview').html(data);
        pdbtabledata = data;
    });
    $(".main_loading_overlay").hide();
}


function initalizeSingleCrystalView() {
    //            initializeSegmentButtons('#single-crystal-tab');
    initializeGoButton('#single-crystal-tab');
    initializeFullscreenButton('#single-crystal-tab');
}

function initializeSingleGroupCrystalView() {
    //            initializeSegmentButtons('#single-crystal-group-tab');
    initializeGoButton('#single-crystal-group-tab', true);
    initializeFullscreenButton('#single-crystal-group-tab');
    initializeInteractionButtons('#single-crystal-group-tab');
}

function initializeTwoCrystalGroupsView() {
    //            initializeSegmentButtons('#two-crystal-groups-tab');
    initializeGoButtonTwoCrystalGroups('#two-crystal-groups-tab', true);
    initializeFullscreenButton('#two-crystal-groups-tab');
    initializeInteractionButtons('#two-crystal-groups-tab');
}

function updateMatrix() {
    switch (currentTab) {
        case "single-crystal-tab":
            // $('#' + currentTab+ ' .controls-panel input[type=checkbox]').each(function() {
            //     var interactionType = $(this).data('interaction-type');
            //     var rects = $('#' + currentTab + ' .heatmap rect.' + interactionType);
            //     if ($(this).is(':checked')) {
            //         rects.show();
            //     } else {
            //         rects.hide();
            //     }
            // });

            if (!filtered_gn_pairs.length) break;

            // Hide all below min treshold
            $('.tab-pane.main_option.active .heatmap .heatmap-interaction').each(function () {
                var pair = $(this).data("res-no-1") + "," + $(this).data("res-no-2");
                var pair_reverse = $(this).data("res-no-2") + "," + $(this).data("res-no-1");
                if (filtered_gn_pairs.includes(pair) || filtered_gn_pairs.includes(pair_reverse)) {
                    $(this).show();
                } else {
                    $(this).hide();
                }
            });

            break;
        case "single-crystal-group-tab":
            // var [tMin,tMax] = $('#' + currentTab + '-crystal-tab #pdbs-range-slider').slider( "option", "values" );
            var [tMin, tMax] = $('#pdbs-range-slider').slider("option", "values");

            if (!filtered_gn_pairs.length) break;

            // Hide all below min treshold
            $('.tab-pane.main_option.active .heatmap .heatmap-interaction').each(function () {
                var pair = $(this).data("gen-no-1") + "," + $(this).data("gen-no-2");
                var pair_reverse = $(this).data("gen-no-2") + "," + $(this).data("gen-no-1");
                if (filtered_gn_pairs.includes(pair) || filtered_gn_pairs.includes(pair_reverse)) {
                    $(this).show();
                } else {
                    $(this).hide();
                }
            });
            break;
        case "two-crystal-groups-tab":

            if (!filtered_gn_pairs.length) break;

            $('.tab-pane.main_option.active .heatmap .heatmap-interaction').each(function () {
                var pair = $(this).data("gen-no-1") + "," + $(this).data("gen-no-2");
                var pair_reverse = $(this).data("gen-no-2") + "," + $(this).data("gen-no-1");
                if (filtered_gn_pairs.includes(pair) || filtered_gn_pairs.includes(pair_reverse)) {
                    $(this).show();
                } else {
                    $(this).hide();
                }
            });
            break;
    }
}

// TODO update for other tabs
function updateFlareplot() {

    var paths = $('.main_option:visible .flareplot-container path').each(function () {
        // var f1 = $(this).data("group-1-freq");
        // var f2 = $(this).data("group-2-freq");
        // var f3 = $(this).data("frequency-diff");
        // if ( (f1 < r1[0] || r1[1] < f1) || (f2 < r2[0] || r2[1] < f2) || (f3 < r3[0] || r3[1] < f3) ) {
        //     $(this).hide();
        // } else {
        //     $(this).show();
        // }
        if ($(this).attr("class")) {
            var path_class = $(this).attr("class").split(' '); //[1].replace("edge-","")
            var gn1 = path_class[1].replace("source-", "");
            var gn2 = path_class[2].replace("target-", "");
            var pair = gn1 + "," + gn2;
            // console.log(pair);
            if (filtered_gn_pairs.includes(pair)) {
                $(this).show();
            } else {
                $(this).hide();
            }
        }
    });

    // var container = '#' + currentTab + ' .flareplot-container'


    // switch (currentTab) {
    //     case "single-crystal-tab":
    //         // interactionsToggleList = [];
    //         // $('#' + currentTab+ ' .controls-panel input[type=checkbox]').each(function() {
    //         //     // init toggle list
    //         //     if ($(this).is(':checked'))
    //         //         interactionsToggleList.push($(this).data('interaction-type'));
    //         // });

    //         // // toggle interactions in flareplot
    //         // if (container in flareplot)
    //         //     flareplot[container].showInteractions(interactionsToggleList);
    //         if (!filtered_gn_pairs.length) break;

    //         var paths = $(container + ' path').each(function() {
    //             // var f1 = $(this).data("group-1-freq");
    //             // var f2 = $(this).data("group-2-freq");
    //             // var f3 = $(this).data("frequency-diff");
    //             // if ( (f1 < r1[0] || r1[1] < f1) || (f2 < r2[0] || r2[1] < f2) || (f3 < r3[0] || r3[1] < f3) ) {
    //             //     $(this).hide();
    //             // } else {
    //             //     $(this).show();
    //             // }
    //             if ($(this).attr("class")) {
    //                 var path_class = $(this).attr("class").split(' '); //[1].replace("edge-","")
    //                 var gn1 = path_class[1].replace("source-", "");
    //                 var gn2 = path_class[2].replace("target-", "");
    //                 var pair = gn1 + "," + gn2;
    //                 // console.log(pair);
    //                 if (filtered_gn_pairs.includes(pair)) {
    //                     $(this).show();
    //                 } else {
    //                     $(this).hide();
    //                 }
    //             }
    //         });
    //         break;
    //     case "single-crystal-group-tab":
    //         // var [tMin,tMax] = $('#' + currentTab + ' #pdbs-range-slider').slider( "option", "values" );
    //         // if (container in flareplot)
    //         //     flareplot[container].updateRange(tMin, tMax);
    //         break;
    //     case "two-crystal-groups-tab":

    //         if (!filtered_gn_pairs.length) break;

    //         var paths = $(container + ' path').each(function() {
    //             // var f1 = $(this).data("group-1-freq");
    //             // var f2 = $(this).data("group-2-freq");
    //             // var f3 = $(this).data("frequency-diff");
    //             // if ( (f1 < r1[0] || r1[1] < f1) || (f2 < r2[0] || r2[1] < f2) || (f3 < r3[0] || r3[1] < f3) ) {
    //             //     $(this).hide();
    //             // } else {
    //             //     $(this).show();
    //             // }
    //             if ($(this).attr("class")) {
    //                 var path_class = $(this).attr("class").split(' '); //[1].replace("edge-","")
    //                 var gn1 = path_class[1].replace("source-", "");
    //                 var gn2 = path_class[2].replace("target-", "");
    //                 var pair = gn1 + "," + gn2;
    //                 // console.log(pair);
    //                 if (filtered_gn_pairs.includes(pair)) {
    //                     $(this).show();
    //                 } else {
    //                     $(this).hide();
    //                 }
    //             }
    //         });
    //         break;
    // }
}

// TODO create interaction toggle for Hiveplot
function updateHiveplot() {
    // PLACEHOLDER
}

// TODO update for other tabs
function updateSchematic() {
    switch (currentTab) {
        case "single-crystal-tab":
            $('#' + currentTab + ' .controls-panel input[type=checkbox]').each(function () {
                var interactionType = $(this).data('interaction-type');
                var paths = $('#' + currentTab + '-crystal-tab  .' + currentViz + '-container path.' + interactionType);
                if ($(this).is(':checked')) {
                    paths.show();
                } else {
                    paths.hide();
                }
            });
            break;
        case "single-crystal-group-tab":
            var [tMin, tMax] = $('#' + currentTab + ' #pdbs-range-slider').slider("option", "values");

            // Hide all below min treshold
            var paths = $('#' + currentTab + ' .' + currentViz + '-container path').each(function () {
                var n = $(this).data("num-interactions");
                if (n < tMin || tMax < n) {
                    $(this).hide();
                } else {
                    $(this).show();
                }
            });
            break;
        case "two-crystal-groups-tab":
            // var r1 = $('#' + currentTab + '#freq-slider-range-1').slider( "option", "values" );
            // var r2 = $('#' + currentTab + '#freq-slider-range-2').slider( "option", "values" );
            // var r3 = $('#' + currentTab + '#freq-slider-range-3').slider( "option", "values" );

            if (!filtered_gn_pairs.length) break;

            // // Hide all below min or above treshold
            var paths = $('#' + currentTab + ' .' + currentViz + '-container path').each(function () {
                // var f1 = $(this).data("group-1-freq");
                // var f2 = $(this).data("group-2-freq");
                // var f3 = $(this).data("frequency-diff");
                // if ( (f1 < r1[0] || r1[1] < f1) || (f2 < r2[0] || r2[1] < f2) || (f3 < r3[0] || r3[1] < f3) ) {
                //     $(this).hide();
                // } else {
                //     $(this).show();
                // }
                if ($(this).attr("class")) {
                    path_class = $(this).attr("class").split(' ')[1].replace("edge-", "");
                    // console.log(path_class);
                    if (filtered_gn_pairs.includes(path_class)) {
                        $(this).show();
                    } else {
                        $(this).hide();
                    }

                }
            });
            // $('#' + currentTab+ ' .' + currentViz + '-container path').hide();
            // $.each(filtered_gn_pairs,function(i,v) {
            //   console.log(v,'#' + currentTab+ ' .' + currentViz + '-container .edge-'+v,$('#' + currentTab+ ' .' + currentViz + '-container .edge-'+v));
            //   $('#' + currentTab+ ' .' + currentViz + '-container .edge-'+v).show();
            // });
            break;
    }
}

function updateGeneralControls(ignore_ngl = false) {
    // update current vizualization
    console.log('ignore_ngl', ignore_ngl);

    updateFlareplot();
    switch (currentViz.toLowerCase()) {
        case "matrix":
            updateMatrix()
            break;
        case "flareplot":
            updateFlareplot()
            break;
        case "hiveplot":
            updateHiveplot()
            break;
        case "schematic_con":
        case "schematic_non":
            updateSchematic()
            break;
        case "table":
            // do nothing
            break;
        case "ngl":
            // DEPRECATED: do nothing
            break;
        default:
            console.log("Error: missing representation update for " + currentViz)
    }

    // Always invoke NGL update
    if (!ignore_ngl) {
        ngl_plots = [];

        $(".ngl-container:visible").each(function () {
            // statements
            ngl_plot = $(this).attr('id').replace("ngl-", "");
            createNGLRepresentations(ngl_plot, 0, ngl_plot)
            if (ngl_plot.startsWith('two_sets')) createNGLRepresentations(ngl_plot, 1, ngl_plot)
        });
        // Do not update when simply changing viz tabs.
    }
}

var currentViz = "matrix";

function updateCurrentTab(id) {
    alt = $('.main_option:visible').attr('id');
    var now = id.replace("-link", "")
    now = now.replace("-tab", "")

    // update settings
    currentTab = now.substr(0, now.lastIndexOf('-'))
    currentViz = now.replace(currentTab + "-", "")
    currentTab = alt;
    redraw_renders();
}

var settingsSubmit = false;
function updateInteractionSettings() {
    settingsSubmit = false;
    // TODO modulate size of modal: display in center and resize to content?
    $("#resModal").find(".modal-dialog").removeClass("modal-wide")
    $("#resModal").find(".modal-title").html("Settings")
    $("#resModal").find(".modal-body").html("<div id='interaction_settings' style='height:100%;width:100%;display: inline-block;'></div>");
    //$("#resModal").find(".modal-footer .btn-default").addClass("hidden")

    // Add OTHER options
    $("#interaction_settings").append('<h5 class="border-bottom">Options</h5>')
    var option_content = '<ul class="list-group">'

    // Normalize the data for the group analysis
    if (currentTab.includes("group")) {
        checked = currentSettings[currentTab]["options"].indexOf("normalize") >= 0 ? "checked" : "";
        var tooltip = '<span class="glyphicon glyphicon-info-sign" data-html="true" data-toggle="popover" data-trigger="hover" data-placement="below" data-content="This option merges the data of structures from same receptor and utilizes the average as a representative for that receptor. These averages of the receptors are subsequently utilized in this tool."></span>';
        option_content += '<li class="list-group-item">' + tooltip + ' Normalize data<div class="material-switch pull-right"><input id="option-normalize" name="option-toggles" ' + checked + ' type="checkbox"/><label for="option-normalize" class="label-primary"></label></div></li>';
    }

    // Only use class A
    checked = currentSettings[currentTab]["options"].indexOf("classa") >= 0 ? "checked" : "";
    var classtooltip = '<span class="glyphicon glyphicon-info-sign" data-html="true" data-toggle="popover" data-trigger="hover" data-placement="below" data-content="When enabled, generic numbers will be in Class A regardless of class of selection"></span>';
    option_content += '<li class="list-group-item">' + classtooltip + ' Class A numbering<div class="material-switch pull-right"><input id="option-classa" name="option-toggles" ' + checked + ' type="checkbox"/><label for="option-classa" class="label-primary"></label></div></li>';

    // Only between helices
    //checked = currentSettings[currentTab]["options"].indexOf("intrahelical") >= 0 ? "checked" : "";
    //var heltooltip = '<span class="glyphicon glyphicon-info-sign" data-html="true" data-toggle="popover" data-trigger="hover" data-placement="below" data-content="When enabled, interactions between residues within the same segment (i.e. TM1-7, H8 or a loop) are also included in the analysis."></span>';
    //option_content += '<li class="list-group-item">' + heltooltip + ' Intrasegment contacts<div class="material-switch pull-right"><input id="option-intrahelical" name="option-toggles" ' + checked + ' type="checkbox"/><label for="option-intrahelical" class="label-primary"></label></div></li>';

    // Toggle backbone interactions
    //checked = currentSettings[currentTab]["options"].indexOf("backbone") >= 0 ? "checked" : "";
    //var bbtooltip = '<span class="glyphicon glyphicon-info-sign" data-html="true" data-toggle="popover" data-trigger="hover" data-placement="below" data-content="When enabled, interactions with and between backbone atoms are also included in the analysis."></span>';
    //option_content += '<li class="list-group-item">' + bbtooltip + ' Backbone contacts<div class="material-switch pull-right"><input id="option-backbone" name="option-toggles" ' + checked + ' type="checkbox"/><label for="option-backbone" class="label-primary"></label></div></li>';

    option_content += "</ul>"
    $("#interaction_settings").append(option_content);

    // Setup for splitting all interactions
    $("#interaction_settings").append('<h5 class="border-bottom">Intersegment contacts</h5>')
    var option_content = '<ul class="list-group">'
    checked = currentSettings[currentTab]["options"].indexOf("inter_scsc") >= 0 ? "checked" : "";
    option_content += '<li class="list-group-item">Sidechain-sidechain<div class="material-switch pull-right"><input id="option-inter_scsc" name="option-toggles" ' + checked + ' type="checkbox"/><label for="option-inter_scsc" class="label-primary"></label></div></li>';
    checked = currentSettings[currentTab]["options"].indexOf("inter_scbb") >= 0 ? "checked" : "";
    option_content += '<li class="list-group-item">Sidechain-backbone<div class="material-switch pull-right"><input id="option-inter_scbb" name="option-toggles" ' + checked + ' type="checkbox"/><label for="option-inter_scbb" class="label-primary"></label></div></li>';
    checked = currentSettings[currentTab]["options"].indexOf("inter_bbbb") >= 0 ? "checked" : "";
    option_content += '<li class="list-group-item">Backbone-backbone<div class="material-switch pull-right"><input id="option-inter_bbbb" name="option-toggles" ' + checked + ' type="checkbox"/><label for="option-inter_bbbb" class="label-primary"></label></div></li>';
    option_content += "</ul>"
    $("#interaction_settings").append(option_content);

    // Setup for splitting all interactions
    $("#interaction_settings").append('<h5 class="border-bottom">Intrasegment contacts</h5>')
    var option_content = '<ul class="list-group">'
    checked = currentSettings[currentTab]["options"].indexOf("intra_scsc") >= 0 ? "checked" : "";
    option_content += '<li class="list-group-item">Sidechain-sidechain<div class="material-switch pull-right"><input id="option-intra_scsc" name="option-toggles" ' + checked + ' type="checkbox"/><label for="option-intra_scsc" class="label-primary"></label></div></li>';
    checked = currentSettings[currentTab]["options"].indexOf("intra_scbb") >= 0 ? "checked" : "";
    option_content += '<li class="list-group-item">Sidechain-backbone<div class="material-switch pull-right"><input id="option-intra_scbb" name="option-toggles" ' + checked + ' type="checkbox"/><label for="option-intra_scbb" class="label-primary"></label></div></li>';
    checked = currentSettings[currentTab]["options"].indexOf("intra_bbbb") >= 0 ? "checked" : "";
    option_content += '<li class="list-group-item">Backbone-backbone<div class="material-switch pull-right"><input id="option-intra_bbbb" name="option-toggles" ' + checked + ' type="checkbox"/><label for="option-intra_bbbb" class="label-primary"></label></div></li>';
    option_content += "</ul>"
    $("#interaction_settings").append(option_content);

    // Add interaction type selection
    $("#interaction_settings").append('<h5 class="border-bottom">Interaction types</h5>')
    var types_content = '<ul class="list-group">'
    var interaction_options = ["Ionic", "Polar", "Aromatic", "Hydrophobic", "van der Waals"]
    for (var i = 0; i < interaction_options.length; i++) {
        var checked = currentSettings[currentTab]["types"].indexOf(interaction_options[i].replace(/\s+/g, '-')) >= 0 ? "checked" : "";
        types_content += '<li class="list-group-item">' + interaction_options[i] + '<div class="material-switch pull-right"><input id="types-' + interaction_options[i].replace(/\s+/g, '-') + '" name="types-toggles" ' + checked + ' type="checkbox"/><label for="types-' + interaction_options[i].replace(/\s+/g, '-') + '" class="label-primary"></label></div></li>';
    }
    types_content += "</ul>"
    $("#interaction_settings").append(types_content);

    // Add strict toggles for H-bond, aromatic, Hydrophobic and VdW
    $("#interaction_settings").append('<h5 class="border-bottom">Strict interaction cut-offs</h5>')
    $(function () { $('div#interaction_settings span.glyphicon.glyphicon-info-sign').popover() })

    var strict_toggles = ["Polar", "Aromatic", "Hydrophobic", "van der Waals"];
    var strict_tooltips = ["<b>Enabled:</b><br/> donor-acceptor dist. ≤3.5Å + donor angle ≥120°<br/><b>Disabled:</b><br/> donor-acceptor dist. ≤4Å", "<b>Enabled:</b><br/>Face-to-face:<br/>dist. ≤4.4Å + angle ≤30°<br/>Face-to-edge:<br/> dist. ≤5.5Å + angle >30°<br/>Cation-π:<br/>dist. to cat. ≤6.6Å + angle ≤30°<br/><b>Disabled:</b><br/> dist. ≤5.5Å OR Cation-π<br/></br>Calculated from ring center", "<b>Enabled:</b><br/> ≥3 atom pairs<br/><b>Disabled:</b><br/> ≥1 atom pair", "<b>Enabled:</b><br/> ≥3 atom pairs<br/><b>Disabled:</b><br/> ≥1 atom pair"]
    var strict_content = '<ul class="list-group">'
    for (var i = 0; i < strict_toggles.length; i++) {
        var checked = currentSettings[currentTab]["strict"].indexOf(strict_toggles[i].replace(/\s+/g, '-')) >= 0 ? "checked" : "";
        var tooltip = '<span class="glyphicon glyphicon-info-sign" data-html="true" data-toggle="popover" data-trigger="hover" data-placement="below" data-content="' + strict_tooltips[i] + '"></span>';
        strict_content += '<li class="list-group-item">' + tooltip + strict_toggles[i] + '<div class="material-switch pull-right"><input id="strict-' + strict_toggles[i].replace(/\s+/g, '-') + '" name="strict-toggles" ' + checked + ' type="checkbox"/><label for="strict-' + strict_toggles[i].replace(/\s+/g, '-') + '" class="label-primary"></label></div></li>';
    }
    strict_content += "</ul>"
    $("#interaction_settings").append(strict_content);

    $("#resModal").modal();

    // Link to save settings when closing
    $('#resModal').on('hidden.bs.modal', closeInteractionSettings)

    // Enable toggling by click on LI
    $("#interaction_settings li").click(function (e) {
        if (e.target == this) {
            var cb = $(this).find(":checkbox")[0];
            cb.checked = !cb.checked;
        }
    });

    // Add submit button
    $("#resModal").find(".modal-footer").prepend('<button type="button" class="btn btn-success btn-process" data-dismiss="modal">Close & Go</button>')
    $("#resModal").find(".modal-footer .btn-process").click(function (e) {
        settingsSubmit = true;
    });
}

function closeInteractionSettings(e) {
    // Reset modal
    $(e.currentTarget).off('hidden'); //DEPRECATED: $(e.currentTarget).unbind();
    $("#resModal").find(".modal-dialog").addClass("modal-wide")
    $("#resModal").find(".modal-footer .btn-process").remove()

    // Save settings
    var strict = []
    $("#resModal input[name=strict-toggles]:checked").each(function (i, node) { strict.push(node.id.replace("strict-", "")) })

    var types = []
    $("#resModal input[name=types-toggles]:checked").each(function (i, node) { types.push(node.id.replace("types-", "")) })

    var options = []
    $("#resModal input[name=option-toggles]:checked").each(function (i, node) { options.push(node.id.replace("option-", "")) })

    currentSettings[currentTab]["strict"] = strict;
    currentSettings[currentTab]["types"] = types;
    currentSettings[currentTab]["options"] = options;

    if (settingsSubmit) {
        settingsSubmit = false
        $("#" + currentTab + ' .go-button').click()
    }
}

// Show plotting options panel
var selecttab = { "single-crystal-tab": "structure", "single-crystal-group-tab": "single", "two-crystal-groups-tab": "double" }
function showVisualizationPanel(plot_destination, table_number, datatype, column_number) {
    // TODO modulate size of modal: display in center and resize to content?
    $("#resModal").addClass("plot_settings")
    $("#resModal").find(".modal-dialog").removeClass("modal-wide")
    $("#resModal").find(".modal-title").html("Plotting options")
    $("#resModal").find(".modal-body").html("<div id='interaction_settings' style='height:100%;width:100%;display: inline-block;'></div>");
    $("#resModal").find(".modal-footer").addClass("hidden");

    // Link to save settings when closing
    $('#resModal').on('hidden.bs.modal', closeVisualizationPanel)

    // collect options for this column
    var options = plot_options['tab' + table_number][selecttab[currentTab]][datatype];

    // Data options : show them when multiple columns or underlying data
    if (options[0].length > 1 || options[1][0].endsWith("_original")) {

        // Collect options and show (button)
        if (options[1][0].endsWith("_original")) {
            // select underlying data options
            $("#interaction_settings").append('<h5 class="border-bottom">Select which data to plot:</h5>')
            $("#interaction_settings").append('<button type="button" class="btn btn-primary btn-block">Set 1</button><br/>')
            $("#interaction_settings").append('<button type="button" class="btn btn-primary btn-block">Set 2</button><br/>')
            $("#interaction_settings").append('<button type="button" class="btn btn-primary btn-block">Difference</button><br/>')
        } else {
            // show column options
            $("#interaction_settings").append('<h5 class="border-bottom">Select which column to plot:</h5>')
            var header = $(".main_option:visible .browser-table-" + table_number)[0].childNodes[0].children[3];
            var num_columns = options[0].reduce((a, b) => a + b, 0);
            for (var i = 0; i < num_columns; i++) {
                //console.log(header.children[parseInt(columnNumber) + i].innerText);
                $("#interaction_settings").append('<button type="button" class="btn btn-primary btn-block" column_selector="' + i + '">' + header.children[parseInt(column_number) + i].innerText + '</button><br/>')
            }
        }
        // add click event
        $("#resModal").find(".btn-block").on('click', (function (a, b, c, d) { return function (e) { showPlottingPanel(a, b, c, d, e); } })(plot_destination, table_number, datatype, column_number))
        $("#resModal").modal();

    } else {
        showPlottingPanel(plot_destination, table_number, datatype, 0);
    }
}

var plottingData = []
function showPlottingPanel(plot_destination, table_number, datatype, column_number, event = "") {
    // collect options for this column
    var options = plot_options['tab' + table_number][selecttab[currentTab]][datatype];

    // TODO Grab plotting data
    plottingData = []

    // Option 1: all data from all columns
    console.log(event)
    if (event == "") {
        // Just plot all data in the current column
        // FOLLOW same idea as colorByData function - maybe merge function?
        console.log("Option 1")
        // pair or position?
        if (options[1][column_number].startsWith("residuepair")) {

        } else {

        }
    } else if (event.target.hasAttribute("column_selector")) {
        // Option 2: all data from 1 of the columns (find column number in column_selector attribute)
        console.log("Option 2 " + selected)
        var selected = event.target.getAttribute("column_selector");
        // FOLLOW same idea as colorByData function - maybe merge options?

    } else {
        // Option 3: grab data from raw input (TODO: filter based on residues in list table)
        console.log("Option 3")
        var selected = event.target.innerText;

    }

    // Show plotting options
    $("#interaction_settings").html('<h5 class="border-bottom">Select the plot type:</h5>')

    // Check current plot - use for coloring current plot type
    var plots = $('.main_option:visible').find(".plot-container").not(".plotly");
    var plotType = plots[plot_destination].id

    // TODO: color current plot option green (success)
    if (options[1][column_number].startsWith("residuepair")) {
        // Visualizing data for a residue pair
        $("#interaction_settings").append('<button type="button" class="btn btn-primary btn-block" table_number="' + plot_destination + '" plot_type="heatmap">Matrix of interactions</button><br/>')
        $("#interaction_settings").append('<button type="button" class="btn btn-primary btn-block" table_number="' + plot_destination + '" plot_type="flareplot">Flare Plot</button><br/>')
        $("#interaction_settings").append('<button type="button" class="btn btn-primary btn-block" table_number="' + plot_destination + '" plot_type="ngl">3D View</button><br/>')
        //$("#interaction_settings").append('<button type="button" class="btn btn-primary btn-block">Schematic (Non-consecutive)</button><br/>')
        //$("#interaction_settings").append('<button type="button" class="btn btn-primary btn-block">Schematic (Consecutive)</button><br/>')
    } else {
        // Visualizing data for a residue position
        $("#interaction_settings").append('<button type="button" class="btn btn-primary btn-block" table_number="' + plot_destination + '" plot_type="ngl">3D View</button><br/>')
        $("#interaction_settings").append('<button type="button" class="btn btn-primary btn-block" table_number="' + plot_destination + '" plot_type="snakeplot">Snake Plot</button><br/>')
        $("#interaction_settings").append('<button type="button" class="btn btn-primary btn-block" table_number="' + plot_destination + '" plot_type="boxplot">Box Plot</button><br/>')
    }
    // add click event
    $("#resModal").find(".btn-block").on('click', initiatePlotWithData)

    $("#resModal").modal();
}

var prefix_list = { "single-crystal-tab": "single_", "single-crystal-group-tab": "single_group_", "two-crystal-groups-tab": "two_sets_" }
function initiatePlotWithData(e) {
    // close modal panel
    $("#resModal").modal('hide');

    // Initialize plot in panel
    var prefix = prefix_list[currentTab];
    var origin = $(e.target);
    var currentPanel = prefix + (parseInt(origin.attr('table_number')) + 1);
    drawPlotPanel(origin.attr('plot_type'), $("#" + currentPanel))

    // TODO: Wait until ready and send data to visualize
    //colorByData(mode, tableNumber, columnNumber, type)
}

function closeVisualizationPanel(e) {
    // Show Close button again
    $("#resModal").find(".modal-footer").removeClass("hidden");
    $("#resModal").removeClass("plot_settings")

    // Reset modal
    $(e.currentTarget).off('hidden');
    $("#resModal").find(".modal-dialog").addClass("modal-wide")
    $("#resModal").find(".modal-footer .btn-process").remove()
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

function separatePair(stringPair) {
    var regex = /([0-9x]+),([0-9x]+)/g;
    var m;
    matches = regex.exec(stringPair);

    return [matches[1], matches[2]];
}

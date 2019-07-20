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

function initializeFullscreenButton(selector) {
    $(selector + ' .btn-fullscreen').click(function() {
        fullScreenElement = $(this).parent().next().children().first();
        fullScreenElement.css('background-color', 'white');
        toggleFullScreen(fullScreenElement.get(0));
    });
}

var stage = {};
var pdb_data = {};
var color_schemes = {};
var reps = {} // store ngl representations
var gpcr_rep = {}
var selectColoring = {} // store coloring settings
var int_labels = []
var basecolor = { "ngl-1": 'red', "ngl-2": 'blue', "ngl-3": 'grey' }

function createNGLview(mode, pdb, pdbs = false) {
    var mode = "ngl-" + mode
    $("#" + mode).html("");
    stage[mode] = new NGL.Stage(mode, { backgroundColor: "white" });
    color_schemes[mode] = [
        [],
        []
    ];
    reps[mode] = [{}, {}]
    gpcr_rep[mode] = [];
    pdb_data[mode] = [];
    int_labels[mode] = []

    var first_structure;
    var num_set1;
    var gn_num_set1;
    var two_structures = false;

    var blue_colors = ['#f0fcfa', '#D3E4EA', '#B6CDDB', '#99B5CC', '#7C9EBD', '#5F86AE', '#426F9F', '#255790', '#084081']
    var red_colors = ['#fbf0fc', '#f3dbec', '#ebc5df', '#dda8bc', '#cc7b7f', '#d3574e', '#be372b', '#ac1808', '#811808']

    $.getJSON("/contactnetwork/pdb/" + pdb,
        function(data) {
            var highlight = ['TM1', 'TM2', 'TM3', 'TM4', 'TM5', 'TM6', 'TM7', 'H8'];
            var segments_sets = {}
            int_labels[mode][0] = {}

            highlight.forEach(function(e) {
                segments_sets[e] = ((e in data['segments']) ? data['segments'][e].join(", ") : "")
            });

            pdb_data[mode][0] = data;
            color_schemes[mode][0]['blue'] = NGL.ColormakerRegistry.addSelectionScheme([
                [blue_colors[1], segments_sets[highlight[0]]],
                [blue_colors[2], segments_sets[highlight[1]]],
                [blue_colors[3], segments_sets[highlight[2]]],
                [blue_colors[4], segments_sets[highlight[3]]],
                [blue_colors[5], segments_sets[highlight[4]]],
                [blue_colors[6], segments_sets[highlight[5]]],
                [blue_colors[7], segments_sets[highlight[6]]],
                [blue_colors[8], segments_sets[highlight[7]]],
                ["white", "*"]
            ])

            color_schemes[mode][0]['grey'] = NGL.ColormakerRegistry.addSelectionScheme([
                ["#ccc", segments_sets[highlight[0]]],
                ["#bbb", segments_sets[highlight[1]]],
                ["#aaa", segments_sets[highlight[2]]],
                ["#888", segments_sets[highlight[3]]],
                ["#666", segments_sets[highlight[4]]],
                ["#444", segments_sets[highlight[5]]],
                ["#333", segments_sets[highlight[6]]],
                ["#111", segments_sets[highlight[7]]],
                ["white", "*"]
            ])

            color_schemes[mode][0]['red'] = NGL.ColormakerRegistry.addSelectionScheme([
                [red_colors[1], segments_sets[highlight[0]]],
                [red_colors[2], segments_sets[highlight[1]]],
                [red_colors[3], segments_sets[highlight[2]]],
                [red_colors[4], segments_sets[highlight[3]]],
                [red_colors[5], segments_sets[highlight[4]]],
                [red_colors[6], segments_sets[highlight[5]]],
                [red_colors[7], segments_sets[highlight[6]]],
                [red_colors[8], segments_sets[highlight[7]]],
                ["white", "*"]
            ])


            pdb_data[mode][0] = data;



            chain_set1 = pdb_data[mode][0]['chain'];
            num_set1 = pdb_data[mode][0]['only_gn'];
            gn_num_set1 = pdb_data[mode][0]['gn_map'];

            var stringBlob = new Blob([pdb_data[mode][0]['pdb']], { type: 'text/plain' });
            stage[mode].loadFile(stringBlob, { ext: "pdb" }).then(function(o) {
                first_structure = o

                // Cleanup
                pdb_data[mode][0]['pdb'] = null;

                // TODO: cartoon for ECL2 is currently not shown (only 3 res, but 4 needed for cartoon) - show ribbon for ICL/ECLs?
                gpcr_rep[mode][0] = o.addRepresentation("cartoon", {
                    sele: ":" + pdb_data[mode][0]['chain'] + " and (" + pdb_data[mode][0]['only_gn'].join(", ") + ")",
                    radiusSize: 1,
                    radiusScale: 0.6,
                    color: color_schemes[mode][0][basecolor[mode]],
                    metalness: 0,
                    colorMode: "hcl",
                    roughness: 1,
                    opacity: 1,
                    side: "front",
                    depthWrite: true
                });

                // REMOVE default hoverPick mouse action
                stage[mode].mouseControls.remove("hoverPick")

                // Add residue labels for GN residues
                pdb_data[mode][0]['only_gn'].forEach(function(resNo, index) {
                    var genNo = pdb_data[mode][0]['gn_map'][index]
                    int_labels[mode][0][resNo] = genNo
                })

                // listen to `hovered` signal to move tooltip around and change its text
                // TODO: implement this to handle two structures
                stage[mode].signals.hovered.add(function(pickingProxy) {
                    var label_index = false

                    if (pickingProxy && pickingProxy.atom) {
                        label_index = pickingProxy.atom.resno;
                        ngl_tooltip.innerText = "RESIDUE:"
                    }

                    if (label_index && (label_index in int_labels[mode][0])) {
                        var mp = pickingProxy.mouse.position
                        ngl_tooltip.innerText += " " + int_labels[mode][0][label_index]
                        ngl_tooltip.style.bottom = window.innerHeight - mp.y + 3 + "px"
                        ngl_tooltip.style.left = mp.x + 3 + "px"
                        ngl_tooltip.style.display = "block"
                    } else {
                        ngl_tooltip.style.display = "none"
                    }
                })

                reps[mode][0].structureComponent = o
                createNGLRepresentations(mode, 0, false)

                // Automatic GPCR positioning
                if ("translation" in pdb_data[mode][0]) {
                    var translation = JSON.parse(pdb_data[mode][0]["translation"])
                    var center_axis = JSON.parse(pdb_data[mode][0]["center_axis"])

                    // calculate rotation and apply
                    v1 = new NGL.Vector3(0, 1, 0)
                    v2 = new NGL.Vector3(center_axis[0], center_axis[1], center_axis[2])
                    var quaternion = new NGL.Quaternion(); // create one and reuse it
                    quaternion.setFromUnitVectors(v2, v1)
                    o.setRotation(quaternion)

                    // calculate translation and apply
                    var v = new NGL.Vector3(-1 * translation[0], -1 * translation[1], -1 * translation[2])
                    v.applyMatrix4(o.matrix)
                    o.setPosition([-1 * v.x, -1 * v.y, -1 * v.z])

                    // calculate H8 position (based on TM1)
                    var tm1_vector
                    var ref_tm1 = pdb_data[mode][0]["only_gn"][pdb_data[mode][0]["gn_map"].indexOf("1x46")]
                    o.structure.eachAtom(function(ap) {
                        tm1_vector = new NGL.Vector3(ap.x, ap.y, ap.z)
                        tm1_vector.applyMatrix4(o.matrix)
                    }, new NGL.Selection(":" + pdb_data[mode][0]['chain'] + " and " + ref_tm1 + " and .CA"))

                    tm1_vector.y = 0 // height position doesn't matter
                    tm1_vector.normalize()

                    // calculate rotation angle around Y-axis (as the GPCR is now upright)
                    v3 = new NGL.Vector3(-1, 0, 0)
                    var m = new NGL.Matrix4()
                    if (tm1_vector.z < 0)
                        m.makeRotationY(v3.angleTo(tm1_vector))
                    else if (tm1_vector.z > 0)
                        m.makeRotationY(-1 * v3.angleTo(tm1_vector))

                    o.setTransform(m)
                }

                var view = new NGL.Matrix4()
                view.set(144.5, 12.5, 8.8, 0, -13.1, 144.3, 10.0, 0, -7.9, -10.8, 144.7, 0, 4, -0.4, -1.1, 1)
                stage[mode].viewerControls.orient(view)
            });

        })

    var newDiv = document.createElement("div");
    newDiv.setAttribute("style", "position: absolute; top: 45px; left: 20px; background-color: #DDD; opacity: .8; padding: 10px;")
    var controls = '<div class="controls">' +
        '<h4>Controls</h4>';

    if (pdbs) {
        controls += '<p>Structure: <select id="ngl_pdb_' + mode + '_ref">';
        for (var i = 0; i < pdbs.length; i++) {
            if (pdbs[i] == pdb)
                controls += '<option value="' + pdbs[i] + '" SELECTED>' + pdbs[i] + '</option>';
            else
                controls += '<option value="' + pdbs[i] + '">' + pdbs[i] + '</option>';
        }
        controls += '</select></p>';
    }

    controls += '<p>Only GNs: <input type=checkbox id="ngl_only_gns" checked></p>'
    controls += '</select></p>'
    newDiv.innerHTML = controls;

    $("#" + mode).append(newDiv);

    selectColoring[mode] = "diff"
    $('#ngl_coloring_' + mode).change(function(e) {
        selectColoring[mode] = $(this).val()
        updateStructureRepresentations(mode)
    });

    // structure selection
    $("#ngl_pdb_" + mode + "_ref").change(function(e) {
        createNGLview(mode, $(this).val(), pdbs, pdbs_set2);
    });

    $("#" + mode + " #ngl_color").change(function(e) {
        gpcr_rep[mode][0].setParameters({
            color: color_schemes[mode][0][$(this).val()]
        });
    });


    $("#" + mode + "-NGL-tab-link").click(function(e) {
        $(function() {
            stage[mode].handleResize();
        });
    });

    $("#" + mode + " #ngl_only_gns").change(function(e) {
        updateStructureRepresentations(mode);
    });

    // Link the 3D viewers together
    stage[mode].mouseObserver.signals.dragged.add(function (){linkNGLMouseControls(mode)});
    stage[mode].mouseObserver.signals.scrolled.add(function (){linkNGLMouseControls(mode)});
    // Click signals not fully functional because of animation -> disabled for now
    //stage[mode].mouseObserver.signals.clicked.add(function (){linkNGLMouseControls(mode)});
    //stage[mode].mouseObserver.signals.doubleClicked.add(function (){linkNGLMouseControls(mode)});

    // Prevent scrolling of document
    document.getElementById( mode ).onwheel = function(event){event.preventDefault();};
    document.getElementById( mode ).onmousewheel = function(event){event.preventDefault();};
}

// Linking the viewers together
function linkNGLMouseControls(origin){
  var mode = origin.substring(0, origin.length - 1)

  for (var graph in stage){
    if (graph!=origin && graph.startsWith(mode)){
      stage[graph].viewerControls.orient(stage[origin].viewerControls.getOrientation());
    }
  }
}

var linkColourScheme = {}

function createNGLRepresentations(mode, structureNumber, update = false) {
    var res_int = []
    if (mode in reps && structureNumber in reps[mode] && reps[mode][structureNumber].structureComponent)
        var o = reps[mode][structureNumber].structureComponent;
    else
        return // NGL not initialized

    // initialize linkMap + colorScheme
    if (!(mode in linkColourScheme)) linkColourScheme[mode] = {}

    var gnOnly = !update || $("#" + mode + " #ngl_only_gns").prop('checked');

    // Empty? Update selection with a fake residue -> hide everything
    if (res_int.length == 0) res_int.push("9999999")

    if (update) {
        reps[mode][structureNumber].int_res.setSelection(":" + pdb_data[mode][structureNumber]['chain'] + " and (" + res_int.join(", ") + ") and (.CA)")
        reps[mode][structureNumber].ball_int.setSelection(":" + pdb_data[mode][structureNumber]['chain'] + " and (" + res_int.join(", ") + ") and sidechainAttached")
    } else {
        // selection representation
        reps[mode][structureNumber].hover = o.addRepresentation("ball+stick", { sele: ":" + pdb_data[mode][structureNumber]['chain'] + " and (" + res_int.join(", ") + ") and (.CA)" });

        reps[mode][structureNumber].int_res = o.addRepresentation("spacefill", {
            sele: ":" + pdb_data[mode][structureNumber]['chain'] + " and (" + res_int.join(", ") + ") and (.CA)",
            color: (structureNumber == 0 ? "#084081" : "#811808"),
            // colorScale: ["#44f", "#444"],
            radiusScale: .2,
            name: "res",
            visible: true
        });

        reps[mode][structureNumber].ball_int = o.addRepresentation("ball+stick", {
            sele: ":" + pdb_data[mode][structureNumber]['chain'] + " and (" + res_int.join(", ") + ") and sidechainAttached",
            color: "element",
            colorValue: (structureNumber == 0 ? "#99B5CC" : "#DDA8BC"),
            visible: false
        })
    }

    // update show/hide when updating representations
    updateStructureRepresentations(mode)
}

// TODO: make this function obsolete and merge remaining code with *createNGLRepresentations*
function updateStructureRepresentations(mode) {
    var structures = 1;

    if (gpcr_rep[mode] != undefined && gpcr_rep[mode][0] != undefined) {
        for (var key = 0; key < structures; key++) {
            var color_scheme = []
            var coloringData = []

            // TODO: currently not working - example code for gradient coloring using external data
            if (coloringData[mode] != undefined) {
                // find selected coloring option
                var colorIndex = 0

                var colorType = $(tabMap[mode] + " .variable_selector").val()
                for (gn_key in coloringData[mode]) {
                    // select residue
                    gn = pdb_data[mode][key]["gn_map"].indexOf(gn_key)
                    ngl_selection = ":" + pdb_data[mode][key]['chain'] + " and " + pdb_data[mode][key]["only_gn"][gn]

                    // add color
                    // 0 - diff | 1- AvgG1/value | 2 - rangeG1 | 3 - avgG2 - 4 - rangeG2
                    switch (colorType) {
                        case "outer_angle":
                        case "bb_angle":
                        case "sc_angle":
                            if (colorIndex == 1 || colorIndex == 3)
                                color_scheme.push([numberToColor2(140, coloringData[mode][gn_key][colorIndex] - 40, false), ngl_selection])
                            else if (colorIndex == 2 || colorIndex == 4)
                                color_scheme.push([numberToColor2(90, coloringData[mode][gn_key][colorIndex], true), ngl_selection])
                            else
                                color_scheme.push([numberToColor2(45, coloringData[mode][gn_key][colorIndex], true), ngl_selection])
                            break;
                        default:
                            break;
                    }
                }

                color_scheme.push(["white", "*"]);
                diff_coloring = NGL.ColormakerRegistry.addSelectionScheme(color_scheme)
                gpcr_rep[mode][key].setParameters({
                    color: diff_coloring
                });
            }

            // Update cartoon using selection
            checked = $("#" + mode + " #ngl_only_gns").prop('checked');
            sele = ":" + pdb_data[mode][key]['chain'];
            int_sele = sele
            if (checked)
                sele = ":" + pdb_data[mode][key]['chain'] + " and (" + pdb_data[mode][key]['only_gn'].join(", ") + ")";

            gpcr_rep[mode][key].setSelection(sele);
        }
    }
}

function numberToColor2(max, value, neg_and_pos = false) {
    if (neg_and_pos) {
        value = value + max
        max = max * 2
    }

    if (value > max)
        value = max
    if (value < 0)
        value = 0

    return colorGradient(value / max, { red: 255, green: 255, blue: 255 }, { red: 0, green: 0, blue: 255 })
}

function colorGradient(fadeFraction, rgbColor1, rgbColor2, rgbColor3) {
    var color1 = rgbColor1;
    var color2 = rgbColor2;
    var fade = fadeFraction;

    // Do we have 3 colors for the gradient? Need to adjust the params.
    if (rgbColor3) {
        fade = fade * 2;

        // Find which interval to use and adjust the fade percentage
        if (fade >= 1) {
            fade -= 1;
            color1 = rgbColor2;
            color2 = rgbColor3;
        }
    }

    var diffRed = color2.red - color1.red;
    var diffGreen = color2.green - color1.green;
    var diffBlue = color2.blue - color1.blue;

    var gradient = {
        red: parseInt(Math.floor(color1.red + (diffRed * fade)), 10),
        green: parseInt(Math.floor(color1.green + (diffGreen * fade)), 10),
        blue: parseInt(Math.floor(color1.blue + (diffBlue * fade)), 10),
    };

    //return 'rgb(' + gradient.red + ',' + gradient.green + ',' + gradient.blue + ')';
    return rgb2hex(gradient.red.toString(16), gradient.green.toString(16), gradient.blue.toString(16));
}

function redraw_ngl() {
    var mode = $('ul#mode_nav').find('li.active').find('a').text().trim();
    var mode = $('.ngl-container:visible').attr('id').replace('ngl-', '');
    if (mode in stage) {
        $(function() {
            stage[mode].handleResize();
        });
    }
}

function intersect(a, b) {
    var t;
    if (b.length > a.length) t = b, b = a, a = t; // indexOf to loop over shorter
    return a.filter(function(e) {
        return b.indexOf(e) > -1;
    });
}

function swap(dict) {
    var ret = {};
    for (var key in dict) {
        ret[dict[key]] = key;
    }
    return ret;
}

var hotspotsdata = {}

function getBackendData() {
    $(".main_loading_overlay").show();

    console.time('Grab data backend')
    $.ajax({
        url: '/hotspots/hotspotsdata',
        dataType: 'json',
        data: {},
        async: true,
        success: function(data) {
            // Re-render heatmap
            hotspotsdata = data;
            console.table(data['residue_matrix']);
            console.timeEnd('Grab data backend');

            // FULL SET
            // renderHeatMapTable();

            // DEMO - Get n number of positions
            n = 40
            var gns = hotspotsdata['sorted_gns'].slice(0);;
            const shuffled = gns.sort(() => 0.5 - Math.random());
            // Get sub-array of first n elements after shuffled
            let selected = shuffled.slice(0, n);
            ordered_selected = [];
            $.each(hotspotsdata['sorted_gns'], function(i, gn) {
                if (selected.includes(gn)) ordered_selected.push(gn);
            })
            renderHeatMapTable(ordered_selected);

            $(".main_loading_overlay").hide();
        }
    });
}

function renderHeatMapTable(gns = false) {

    console.time('Populate table');
    div = $("#heatmap-container");

    if ($.fn.DataTable.isDataTable("#heatmap_table")) {
        $("#heatmap_table").DataTable().destroy();
    }
    div.html('<table class="row-border text-center dataTable no-footer text-nowrap" width="100%" id="heatmap_table"><thead></thead><tbody></tbody></table>');
    table = div.find("table");
    thead = div.find("thead");
    th_gns_1 = '';
    th_gns_2 = '';

    if (!gns) {
        gns = hotspotsdata['sorted_gns'];
    }
    $.each(gns, function(i, gn) {
        if (gn.includes("x")) {
            th_gns_1 += `<th colspan=3 class="text-center">${gn}</th>`;
            th_gns_2 += `<th></th><th></th><th></th>`;
        }
    });

    tr = `<tr>
                    <th rowspan="2">entry_name</th>
                    <th rowspan="2">receptor_family</th>
                    <th rowspan="2">ligand_type</th>
                    ${th_gns_1}
                  </tr>
                  <tr>
                  ${th_gns_2}
                  </tr>`;
    thead.append(tr);

    tbody = div.find("tbody");
    tbody.hide();
    trs = ''
    $.each(hotspotsdata['residue_matrix'], function(entry_name, row) {
        td_gns = '';
        $.each(gns, function(i, gn) {
            v = row[gn];
            if (gn.includes("x")) {
                td_gns += `<td class='border-left num_cell' bgcolor='${v[2][0]}'>${v[2][1]}</td>
                          <td class='num_cell' bgcolor='${v[3][0]}'>${v[3][1]}</td>
                          <td class='num_cell' bgcolor='${v[4][0]}'>${v[4][1]}</td>`;
            }
        });

        tr = `<tr id="${entry_name}">
                      <td nowrap>${entry_name}</td>
                      <td nowrap>${row['receptor_family']}</td>
                      <td nowrap>${row['ligand_type']}</td>
                      ${td_gns}
                    </tr>`;
        trs += tr;
    });
    tbody.append(trs);
    console.timeEnd('Populate table');
    console.time('Apply DataTable');
    dtable = table.DataTable({
        'scrollX': true,
        scrollY: '50vh',
        // "sDom": 't', // To disable the pages on the button..
        paging: false,
        pageLength: 200,
        "bLengthChange": false,
        "bPaginate": false,
        "bInfo": false,
        "order": [],
    });
    tbody.show();
    dtable.columns.adjust().draw();
    console.timeEnd('Apply DataTable');

    console.time('Apply overlay');
    init_overlay();
    console.timeEnd('Apply overlay');
}

function init_overlay() {
    div = $("#heatmap-container");
    element = "#heatmap-container";
    if (!$.fn.DataTable.isDataTable("#heatmap_table")) {
        return
    }

    table = div.find("table").DataTable();

    table.on('draw.dt', function(e, oSettings) {
        create_overlay("#heatmap_table");
    });


    $(element + ' .dataTables_scrollBody').append('<div class="fixed_cols_overlay"><table id="overlay_table" class="overlay_table row-border text-center dataTable no-footer text-nowrap"><tbody></tbody></table></div>');
    $(element + " .fixed_cols_overlay").hide();

    $(element + ' .dataTables_scrollBody').before("<div class='top_scroll'><div>&nbsp;</div></div>");

    $('.top_scroll').css({
        'width': '100%',
        'overflow-x': 'scroll',
        'overflow-y': 'auto',
    });
    $('.top_scroll div').css({
        // 'background-color': 'red',
        'font-size': '1px',
        'line-height': '1px',
    });

    $('.fixed_cols_overlay').css({
        'top': '0px',
        'position': 'absolute',
        'background': '#f8f8f8',
        '-webkit-box-shadow': '5px 0 2px -2px #888',
        'box-shadow': '5px 0 2px -2px #888',
    });

    $('.fixed_cols_overlay tbody tr').css({
        'background-color': '#f8f8f8',
    });

    create_overlay(element + ' #heatmap_table');
    track_scrolling(element);

    $(function() {
        var tableContainer = $(".dataTables_scrollBody");
        var table = $(".dataTables_scrollBody table");
        var fakeContainer = $(".top_scroll");
        var fakeDiv = $(".top_scroll div");

        var tableWidth = table.width();
        fakeDiv.width(tableWidth);

        fakeContainer.scroll(function() {
            tableContainer.scrollLeft(fakeContainer.scrollLeft());
        });
        tableContainer.scroll(function() {
            fakeContainer.scrollLeft(tableContainer.scrollLeft());
        });
    })
}

function create_overlay(table_id) {
    // This function fires upon filtering, to update what rows to show as an overlay
    $(".overlay_table tbody tr").remove();
    var $target = $(".overlay_table tbody");
    $(table_id + " tbody tr").each(function() {
        var $tds = $(this).children(),
            $row = $("<tr id='overlay_" + $tds.eq(6).html() + "'></tr>");
        $row.append($tds.eq(0).clone()).append($tds.eq(1).clone()).append($tds.eq(2).clone()).appendTo($target);
    });
    $(".overlay_table .border-right").removeClass("border-right");

}

function track_scrolling(element) {
    var left = 0;
    var old_left = 0;
    toggle_enabled = true;
    $(element + ' .dataTables_scrollBody').scroll(function() {
        // If user scrolls and it's >100px from left, then attach fixed columns overlay
        left = $(element + ' .dataTables_scrollBody').scrollLeft();
        if (left != old_left) $(".fixed_cols_overlay").hide();
        old_left = left;

        if (left > 100 && toggle_enabled) {
            $(".fixed_cols_overlay").css({
                left: left + 'px'
            });
            if ($(".fixed_cols_overlay").is(":hidden")) $(".fixed_cols_overlay").show();
        }
    });
}

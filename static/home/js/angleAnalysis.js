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

function initializeGoButton(selector, generic=false) {
    $(selector + ' .go-button').click(function() {
        var pdb = JSON.parse($('#pdb-input').val());
        //var segments = JSON.parse($(selector + ' .segments-input').val());
        var segments = ['TM1','TM2','TM3','TM4','TM5','TM6','TM7','TM1','ICL1','ECL1','ICL2','ECL2','ICL3','ECL3','N-term','C-term'];
        if (pdb.length > 0 && segments.length > 0) {

            renderTable(pdb);

            var second = JSON.parse($('#second-input').val());

            createNGLview("single",pdb[0], second);
        }
    });
}

function initializeFullscreenButton(selector) {
    var fullScreenElement = $(selector + ' .heatmap-container').get(0);
    $(selector + ' .btn-fullscreen').click(function() {
        toggleFullScreen(fullScreenElement);
    });
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
    $('#pdb-input').val(JSON.stringify(pdbs));
}

function thisPDB2(elem) {

    var ReceptorName = $(elem).attr('long');
    var pdbName = $(elem).attr('id');
    $('.pdb_selected').not(elem).prop("checked",false);
    var pdbs = [];
    if ($(elem).prop("checked")) {
        pdbs.push(pdbName);
        // Update view
        $("#second-count").html(ReceptorName + ' - ' + pdbName + ' selected.');
    } else {
        // Update view
        $("#second-count").html('No structure selected.');
    }
    $('#second-input').val(JSON.stringify(pdbs));
}


$.fn.dataTable.ext.order['dom-checkbox'] = function  ( settings, col )
    {
        return this.api().column( col, {order:'index'} ).nodes().map( function ( td, i ) {
            return $('input', td).prop('checked') ? '1' : '0';
        } );
    };

var oTable = [];
function showPDBtable(element, table) {
    if ( ! $.fn.DataTable.isDataTable( element+' .tableview table' ) ) {
        oTable[table] = $(element+' .tableview table').DataTable({
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

        yadcf.init(oTable[table],
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

    yadcf.exResetAllFilters(oTable[table]);

    };
}


function rgb2hex(r,g,b) {

    if (r.length == 1)
        r = '0' + r;

    if (g.length == 1)
        g = '0' + g;

    if (b.length == 1)
        b = '0' + b;

    return '#' + r + g + b;
}

function numberToColor(max,value) {

    var hexc = 255/max;

    if(value<(max/3)){
        return rgb2hex("00","00",Math.floor((max-value)*hexc).toString(16));
    } else if (value < ((2*max)/3)) {
        return rgb2hex(Math.floor(value*hexc).toString(16),Math.floor((max-Math.abs(value-max/2))*hexc).toString(16),Math.floor((max-value)*hexc).toString(16));
    } else {
        return rgb2hex(Math.floor(value*hexc).toString(16),"00","00");
    }
}



var stage = [];
var color_schemes = [];
var schemeId_grey

function createNGLview(mode,pdb, pdb2, pdbs = false) {
    console.log(pdb)
    console.log(pdb2)

    var gpcr_rep
    $("#ngl-"+mode).html("");
    stage[mode] = new NGL.Stage( "ngl-"+mode, { backgroundColor: "white" } );

    var pdb_data
    var blue_colors = ['#f7fcf0','#e0f3db','#ccebc5', '#a8ddb5',    '#7bccc4',    '#4eb3d3', '#2b8cbe',    '#0868ac',    '#084081']
    var reps = {} // store ngl representations
    var original_o

    $.getJSON( "/contactnetwork/pdb/"+pdb,
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
        var mediancolor = [];
        var signifcolor = [];
        var medianwindowcolor = [];
        var hsecolor = [];
        var sasacolor = [];
        var colorstring;
        var ctemp;
        var medmax;
        var score = [0,0,0,0,0];
        var axis = 100; //something larger than possible
        var j = 0;


//         make variable
        var wlength = 5;
        var cutoff = 100;

        wlength = $('#window-input').val();
        cutoff  = $('#score-input').val();

        console.log(wlength)
        console.log(cutoff)
//         var max = 0;
//         residue_data.forEach(function(e){
//             if (max > e[3]){
//                 max = e[3]
//             }
//         }

        // groups
        residue_data.forEach(function(e){
            angle_color.push([numberToColor(180,e[2]) , ""+e[1]])

            mediancolor.push([numberToColor(90,e[3]) , ""+e[1]])

            signifcolor.push([numberToColor(0.5,e[4]-0.5) , ""+e[1]])

            hsecolor.push([numberToColor(40,e[5]) , ""+e[1]])

            sasacolor.push([numberToColor(100,e[6]) , ""+e[1]])

//             if loop data added e[0][0] gives wrong coloring
            if(axis != Number(e[0][0])){
                axis = Number(e[0][0])
                j = 0
                score = [0,0,0,0,0]
            }

            score[j%wlength] = e[3]
            ctemp = score.reduce((a,b)=>a+b);
            if(ctemp > cutoff){
                medianwindowcolor.push([numberToColor(260,ctemp) , ""+e[1]])
            } else {
                medianwindowcolor.push(["#00F" , ""+e[1]])
            }
            j +=1;
        });



        angle_color.push([ "white", "*" ]);
        mediancolor.push([ "white", "*" ]);
        signifcolor.push([ "white", "*" ]);
        medianwindowcolor.push([ "white", "*" ]);
        hsecolor.push([ "white", "*" ]);
        sasacolor.push([ "white", "*" ]);

        color_schemes['angles'] = NGL.ColormakerRegistry.addSelectionScheme(angle_color)
        color_schemes['mediancolor'] = NGL.ColormakerRegistry.addSelectionScheme(mediancolor)
        color_schemes['medianwindowcolor'] = NGL.ColormakerRegistry.addSelectionScheme(medianwindowcolor)
        color_schemes['significance'] = NGL.ColormakerRegistry.addSelectionScheme(signifcolor)
        color_schemes['hse'] = NGL.ColormakerRegistry.addSelectionScheme(hsecolor)
        color_schemes['sasa'] = NGL.ColormakerRegistry.addSelectionScheme(sasacolor)



        if (pdb2.length > 0) {

            $.get('angledata?pdbs[]='+pdb2[0], function(secondArray) {
                var second_residues = secondArray["data"];

                scnd_angle = []
                scnd_hse   = []
                scnd_sasa  = []

                second_residues.forEach(function(scnd){
                    residue_data.forEach(function(e){
                        if (e[0] == scnd[0]){
                            scnd_angle.push([numberToColor(100, Math.abs(e[2] - scnd[2])) , ""+e[1]]);
                            scnd_hse.push([numberToColor(40, Math.abs(e[5] - scnd[5])) , ""+e[1]]);
                            scnd_sasa.push([numberToColor(100, Math.abs(e[6] - scnd[6])) , ""+e[1]]);
                        }
                    });
                });

                scnd_angle.push([ "white", "*" ]);
                scnd_hse.push([ "white", "*" ]);
                scnd_sasa.push([ "white", "*" ]);

                color_schemes['scnd_angle'] = NGL.ColormakerRegistry.addSelectionScheme(scnd_angle)
                color_schemes['scnd_hse'] = NGL.ColormakerRegistry.addSelectionScheme(scnd_hse)
                color_schemes['scnd_sasa'] = NGL.ColormakerRegistry.addSelectionScheme(scnd_sasa)

            });

        }


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

            reps.ball_all = o.addRepresentation("ball+stick", {
                sele: ":"+pdb_data['chain']+" and sidechainAttached",
                visible: false
                })

            reps.ball = o.addRepresentation("ball+stick", {
                sele: ":"+pdb_data['chain']+" and ("+pdb_data['only_gn'].join(", ")+") and sidechainAttached",
                visible: false
                })

            o.autoView();

            // mousover and click on datatable row to highlight residue in NGL viewer
            var temprepr;
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

    var optional = ''
    if (pdb2.length >0){
        optional = '<option value="scnd_angle">Difference angle from'+ pdb2[0] +'</option><option value="scnd_hse">Difference hse from'+ pdb2[0] +'</option><option value="scnd_sasa">Difference sasa from'+ pdb2[0] +'</option>'
    }

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

    controls += '<p>Colors: <select id="ngl_color"><option value="grey">greys</option><option value="blue">blue</option><option value="angles">angles</option><option value="mediancolor">mediancolor</option><option value="medianwindowcolor">medianwindowcolor</option><option value="significance">significance</option><option value="hse">hse</option><option value="sasa">sasa</option>'+ optional +'</select></p>'
                        +'<p>Only GNs: <input type=checkbox id="ngl_only_gns"></p>'
                        +'<p>Show all side-chains: <input type=checkbox id="toggle_sidechains"></p>'
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
        } else {
            sele = ":"+pdb_data['chain'];
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

}

var residue_data

function renderTable(pdb) {
    $.get('angledata?pdbs[]='+pdb[0], function(newDataArray) {
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
        },
        {
            column_number : 3,
            filter_type: "multi_select",
            select_type: 'select2',
            filter_default_label: "diff med",
            filter_reset_button_text: false,
        },
        {
            column_number : 4,
            filter_type: "multi_select",
            select_type: 'select2',
            filter_default_label: "sig med",
            filter_reset_button_text: false,
        },
        {
            column_number : 5,
            filter_type: "multi_select",
            select_type: 'select2',
            filter_default_label: "hse",
            filter_reset_button_text: false,
        },
        {
            column_number : 6,
            filter_type: "multi_select",
            select_type: 'select2',
            filter_default_label: "sasa",
            filter_reset_button_text: false,
        }
    ],
    {
        cumulative_filtering: false
    });

    yadcf.exResetAllFilters(residuetable);
}

$('#single-crystal-pdb-modal-table').on('shown.bs.modal', function (e) {
    showPDBtable('#single-crystal-pdb-modal-table', "firsttable");
})

$('#second-structure-pdb-modal-table').on('shown.bs.modal', function (e) {
    showPDBtable('#second-structure-pdb-modal-table',"secondtable");
})


function initializePdbChooserTables() {
    $.get('/contactnetwork/pdbtabledata', function ( data ) {
    $('#single-crystal-pdb-modal-table .tableview').html(data);
    pdbtabledata = data;
    });
}

function initializeSecondTable() {
    $.get('/contactnetwork/pdbtabledata', function ( data ) {
    $('#second-structure-pdb-modal-table .tableview').html(data);
    secondtable = data;
    });
}

function initalizeSingleCrystalView() {

    initializeGoButton('#single-crystal-tab');
    initializeFullscreenButton('#single-crystal-tab');
}

var selectortable
var residuetable


$(document).ready(function() {
    // residue Table
    initializeResidueTable();
    // Get PDBs for table build
    initializePdbChooserTables();
    initializeSecondTable();

    // Single PDB files
    initalizeSingleCrystalView();

});

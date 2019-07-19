var selectedPDBs = []
function initializeGoButton(selector, generic=false) {
    $(selector + ' .go-button').click(function() {
        var pdb = JSON.parse($('#pdb-input').val());
        if (pdb.length == 0)
          pdb = ['3VW7', '5IUB', '6AKX', '6IGL', '4RWD']

        selectedPDBs = pdb
        console.log(selectedPDBs)
        renderTable();
    });
}

function initializeFullscreenButton(selector) {
    var fullScreenElement = $(selector + ' .heatmap-container').get(0);
    $(selector + ' .btn-fullscreen').click(function() {
        toggleFullScreen(fullScreenElement);
    });
}

/*function thisPDB(elem) {
    var group = $(elem).closest('.tableview').attr('group-number');
    var ReceptorName = $(elem).attr('long');
    var pdbName = $(elem).attr('id');
    //$('.pdb_selected').not(elem).prop("checked",false);

    var pdbs = [];
    if (group==0){
      $('.pdb_selected:checked', oTable["firsttable"].cells().nodes()).each(function() {
          pdbs.push($(this).attr('id'));
      });
      if (pdbs.length==1 && $(elem).prop("checked")) {
          // Update view
          $(".crystal-count:visible").html(ReceptorName + ' - ' + pdbName + ' selected.');
      } else if (pdbs.length>=1){
          // Update view
          $(".crystal-count:visible").html(pdbs.length + ' structures selected.');
      } else {
          // Update view
          $(".crystal-count:visible").html('No structure selected.');
      }
      $('#pdb-input').val(JSON.stringify(pdbs));
    }
}*/
/*
$.fn.dataTable.ext.order['dom-checkbox'] = function  ( settings, col )
    {
        return this.api().column( col, {order:'index'} ).nodes().map( function ( td, i ) {
            return $('input', td).prop('checked') ? '1' : '0';
        } );
    };

function check_all(elem) {
  show_all = $(elem).prop("checked");

  if (show_all) {
    $('.pdb_selected:visible').prop("checked",true);
  } else {
    $('.pdb_selected:visible').prop("checked",false);
  }

  thisPDB(elem)
}*/

var stage = [];
var color_schemes = [];
var schemeId_grey
var chain_selection = ""
function createNGLview(mode, pdb, pdb2, pdbs = false) {
    var gpcr_rep
    $("#ngl-"+mode).html("");
    console.log(pdb)
    console.log(pdb2)
    console.log(mode)
    stage[mode] = new NGL.Stage( "ngl-"+mode, { backgroundColor: "white" } );

    var pdb_data
    var blue_colors = ['#f7fcf0','#e0f3db','#ccebc5', '#a8ddb5',    '#7bccc4',    '#4eb3d3', '#2b8cbe',    '#0868ac',    '#084081']
    var rb_colors = ['#736DA7','#5EB7B7','#CE9AC6', '#DD7D7E', '#E6AF7C', '#DEDB75', '#80B96F', '#000000']  // #C897B8
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

        color_schemes['rainbow'] = NGL.ColormakerRegistry.addSelectionScheme([
                [rb_colors[0], segments_sets[highlight[0]]],
                [rb_colors[1], segments_sets[highlight[1]]],
                [rb_colors[2], segments_sets[highlight[2]]],
                [rb_colors[3], segments_sets[highlight[3]]],
                [rb_colors[4], segments_sets[highlight[4]]],
                [rb_colors[5], segments_sets[highlight[5]]],
                [rb_colors[6], segments_sets[highlight[6]]],
                [rb_colors[7], segments_sets[highlight[7]]],
                [ "white", "*" ]
                ])

        var stringBlob = new Blob( [ pdb_data['pdb'] ], { type: 'text/plain'} );
        stage[mode].loadFile( stringBlob, { ext: "pdb" }  ).then( function( o ){
            original_o = o
            // workaround
            pdb_data['chain'] = pdb_data['chain'].substring(0,1)
            chain_selection = ":" + pdb_data['chain'] + " and "

            gpcr_rep = o.addRepresentation( "cartoon", {
              sele: ":"+pdb_data['chain']+" and ("+pdb_data['only_gn'].join(", ")+") and (.CA)",
              // radiusType: '',
              radiusSize: 1,
              radiusScale: 0.7,
              // color: "atomindex",
              // colorScale: "Accent",
              // color: "residueindex",
              // colorScale: "greys",
              color: color_schemes['rainbow'],
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

            // alignment of GPCR structure
            if ("translation" in pdb_data){
              var translation = JSON.parse(pdb_data["translation"])
              var center_axis = JSON.parse(pdb_data["center_axis"])

              // calculate rotation and apply
              v1 = new NGL.Vector3(0,1,0)
              v2 = new NGL.Vector3(center_axis[0], center_axis[1], center_axis[2])
              var quaternion = new NGL.Quaternion(); // create one and reuse it
              quaternion.setFromUnitVectors( v2, v1 )
              o.setRotation(quaternion)

              // calculate translation and apply
              var v = new NGL.Vector3( -1*translation[0], -1*translation[1], -1*translation[2])
              v.applyMatrix4(o.matrix)
              o.setPosition([-1*v.x, -1*v.y, -1*v.z])

              // calculate H8 position (based on TM1)
              var tm1_vector
              var ref_tm1 = pdb_data["only_gn"][pdb_data["gn_map"].indexOf("1x46")]
              o.structure.eachAtom(function (ap) {
                tm1_vector = new NGL.Vector3(ap.x, ap.y, ap.z)
                tm1_vector.applyMatrix4(o.matrix)
              }, new NGL.Selection(":"+pdb_data['chain']+" and "+ ref_tm1 +" and .CA"))

              tm1_vector.y = 0 // height position doesn't matter
              tm1_vector.normalize()

              // calculate rotation angle around Y-axis (as the GPCR is now upright)
              v3 = new NGL.Vector3(-1, 0, 0)
              var m = new NGL.Matrix4()
              if (tm1_vector.z < 0)
                m.makeRotationY(v3.angleTo(tm1_vector))
              else if (tm1_vector.z > 0)
                m.makeRotationY(-1*v3.angleTo(tm1_vector))

              o.setTransform(m)
            }

            o.autoView(":"+pdb_data['chain']+" and ("+pdb_data['only_gn'].join(", ")+") and (.CA)")

            // mousover and click on datatable row to highlight residue in NGL viewer
            var temprepr;
            $("#single-table-tab-table tbody").on("mouseover", "tr", function(event){
                gn = pdb_data["gn_map"].indexOf(residuetable.row(this).data()[0])
                if (gn > -1) {
                  ngl_selection = ":" + pdb_data['chain'] + " and " + pdb_data["only_gn"][gn]
                  temprepr = o.addRepresentation("ball+stick", {sele: ngl_selection});
                }
            }).mouseout(function(event){
                o.removeRepresentation(temprepr)
            });


        } );

    });

    var newDiv = document.createElement("div");
    newDiv.setAttribute("style", "position: absolute; top: 50px; left: 20px")
    var controls = '<div class="controls">'
                    + '<h3>Controls</h3>';

    controls += '<p>Colors: <select id="ngl_color"><option value="rainbow">rainbow</option><option value="grey">greys</option><option value="blue">blue</option></select></p>'
                        +'<p>Only GNs: <input type=checkbox id="ngl_only_gns" checked></p>'
                        +'<p>Show all side-chains: <input type=checkbox id="toggle_sidechains"></p>'
                        +'</div>';

    newDiv.innerHTML = controls;

    $("#ngl-"+mode).append(newDiv);

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

var current = 0
function renderTable() {
  var pdb = selectedPDBs[current]
  $('#ngl-title').html(pdb.toUpperCase())
  createNGLview("single", pdb, undefined);
}

$('#single-crystal-group-pdbs-modal-table').on('shown.bs.modal', function (e) {
    showPDBtable('#single-crystal-group-pdbs-modal-table');
})

var selectortable
var residuetable
$(document).ready(function() {
    // Get PDBs for table build
    $.get('/contactnetwork/pdbtabledata', function ( data ) {
      $('#single-crystal-group-pdbs-modal-table .tableview').html(data);

      pdbtabledata = data;
    });

    // Single PDB files
    initializeGoButton('#single-crystal-group-tab');
    initializeFullscreenButton('#single-crystal-group-tab');

     // toggle
     $('#forward').click(function() {
       current++;
       if (current >= selectedPDBs.length)
         current = 0;
       renderTable();
     });
     $('#backward').click(function() {
        current--;
        if (current < 0)
          current = selectedPDBs.length - 1;

        renderTable();
     });
});

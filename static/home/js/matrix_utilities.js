let pdb_sel = [];

// Global Variable
let data;
let svg;
let xScale;
let yScale;
let fscale;
let xAxis;
let yAxis;
let xAxisGrid;
let yAxisGrid;
let pdbScale;
let sigScale;
let colScale;
let tooltip;

let sigmatch_table;
let sigmatch_data;

let old_sets = [];
let pos_set = [];
let neg_set = [];

const get_gn = function(){
  const segments = [];
  const orig = [];
  const patt = /(^\d{1,})|(x\d{2,})/g

  const selection = data.receptor
  const value = 'rec_gn'

  for (let index = 0; index < selection.length; index++) {
    let comb_gn = '';
    curr = selection[index][value];
    match = patt.exec(curr);
    while (match != null) {
      comb_gn += match[0];
      match = patt.exec(selection[index][value]);
    }
    if (comb_gn.length > 1) {
      segments.push(comb_gn);
      orig.push(curr);
    }
  };
  return {
    orig: _.uniq(orig),
    label: _.uniq(segments)
  };
};

const get_ignore = function(){
  // get all residues for all valid receptor proteins
  let data_non = non_interactions;

  // filter the residues based on the valid generic numbers
  data_non = _.filter(data_non, function(d) {
    return xScale(d.rec_gn);
  });

  // filter the residues based on teh valid pdbs
  data_non = _.filter(data_non, function(d) {
    return pdbScale(d.pdb_id);
  });

  // reformat the data such that rec_gn are keys and
  // values are an array of entry_names with an
  // interaction at that position
  const data_int = data.receptor.reduce( function (acc, current) {
    acc[current.rec_gn] = [current.entry_name].concat(acc[current.rec_gn])
    return acc;
  }, {});

  // reformat the non interacting data in the same way
  // but also exclude all proteins that already
  // have an interaction registered in the data_int object
  const ignore_markers = data_non.reduce(function (acc, current) {
    let protein_array = data_int[current.rec_gn]
    if ( !protein_array.includes(current.entry_name) ) {
      acc[current.rec_gn] = [current.entry_name].concat(acc[current.rec_gn])
    }
    return acc;
  }, {});

  // example return
  // ignore_markers = {
  //    '2.39x39': [ "oprm_mouse", "oprm_mouse" ],
  //}

  return ignore_markers;
}

const run_seq_sig = function(){
  let segments = get_gn();
  let ignore_markers = get_ignore();
  //let pos_set = ["5ht2c_human", "acm4_human", "drd1_human"];
  //let neg_set = ["agtr1_human", "ednrb_human", "gnrhr_human"];
  //console.log(segments);

  if (pos_set.length < 1){
    alert('No valid set selected.');
    return;
  };

  if (_.isEqual(old_sets.sort(), pos_set.sort())){
    alert('The selected set is identical to the previously selected one.');
    return;
  };

  // api request
  let req = $.ajax({
    type: 'POST',
    url: '/signprot/matrix/seqsig/',
    data: {
      csrfmiddlewaretoken: csrf_token,
      pos: pos_set,
      seg: segments.label,
      ignore: JSON.stringify(ignore_markers),
    },
    beforeSend: function(){
      old_sets = pos_set;
      document.querySelector("#calc_spin").style.display = "inline-block";
    },
    success: function(data){
      $('.svg-content.seqsig').remove();
      $('.svg-content.conseq').remove();

      // sorting data correctly, first by frequency; then by length
      let tmp_data = {}
      for (let rec_gn in data.feat) {
        let position = _.orderBy(
          data.feat[rec_gn],
          ['freq', 'sort_score'],
          ['desc', 'asc']
        )
        tmp_data[rec_gn] = position
      }
      data.feat = tmp_data
      console.log(data)
      document.querySelector('#seqsig-container').style.display = "block";
      document.querySelector('#conseq-container').style.display = "block";
      //document.querySelector('#seqsig').scrollIntoView({behavior: 'smooth'});
      // d3 draw
      svg = signprotmat.d3.setup("div#seqsig-svg", 'seqsig');
      signprotmat.d3.draw_seq_sig(
        data,
        svg,
        xScale,
      );
      svg = signprotmat.d3.setup("div#conseq-svg", 'conseq');
      signprotmat.d3.draw_seq_cons(
        data,
        svg,
        xScale,
        xAxis
      );
      initialize_consensus(data.feat);
    },
    error: function(error){
      console.log(error)
      alert(error);
    },
    complete: function(){
      document.querySelector("#calc_spin").style.display = "none";
    }
  });
};

const initialize_consensus = function(data){
  const row_height = 30;
  for (let key of Object.keys(data)) {
    con_seq[key] = [data[key][0]];
  };
  signprotmat.d3.conSeqUpdate(row_height);
};

const run_sig_match = function(){
  let cutoff = $('#cutoff-val').val();
  if (cutoff == '') {
    cutoff = 0;
  }

  cutoff = 0;
  let segments = get_gn();

  let req = $.ajax({
    type: 'POST',
    url: '/signprot/matrix/sigmat/',
    data: {
      csrfmiddlewaretoken: csrf_token,
      pos: pos_set,
      seg: segments.label,
      cutoff: cutoff,
    },
    beforeSend: function(){
      document.querySelector("#sigm_spin").style.display = "inline-block";
    },
    success: function(data){
      console.log(data)
      document.querySelector('#sigmatch-container').style.display = "inline-block";
      sigmatch_data = Object.keys(data).map(key => data[key])
      sigmatch_table = $('#sigmatch_table').DataTable({
        dom: 'Bfrtip',
        data: sigmatch_data,
        scrollY: '60vh',
        destroy: true,
        columnDefs: [
          {
            data: null,
            targets: 0,
            defaultContent: '',
            orderable: false,
            className: 'select-checkbox',
          }, {
            data: 'class',
            title: 'Class',
            targets: 1,
          }, {
            data: 'prot',
            title: 'Protein',
            targets: 2,
          }, {
            data: 'nscore',
            title: 'Score',
            targets: 4,
          // }, {
          //   data: 'score',
          //   title: 'Score',
          //   targets: 4,
          }, {
            data: 'entry',
            title: 'Entry Name',
            targets: 3,
          }, {
            data: 'GuideToPharma.Gs',
            title: '   Gs   ',
            targets: 5,
            className: 'gtop dt-center',
            visible: true,
          }, {
            data: 'GuideToPharma.Gi/Go',
            title: 'Gi / Go ',
            targets: 6,
            className: 'gtop dt-center',
            visible: true,
          }, {
            data: 'GuideToPharma.Gq/G11',
            title: 'Gq / G11',
            targets: 7,
            className: 'gtop dt-center',
            visible: true,
          }, {
            data: 'GuideToPharma.G12/G13',
            title: 'G12 / G13',
            targets: 8,
            className: 'gtop dt-center',
            visible: true,
          }, {
            data: 'Aska.Gs',
            title: '   Gs   ',
            targets: 9,
            className: 'aska dt-center',
            visible: false,
          }, {
            data: 'Aska.Gi/Go',
            title: 'Gi / Go ',
            targets: 10,
            className: 'aska dt-center',
            visible: false,
          }, {
            data: 'Aska.Gq/G11',
            title: 'Gq / G11',
            targets: 11,
            className: 'aska dt-center',
            visible: false,
          }, {
            data: 'Aska.G12/G13',
            title: 'G12 / G13',
            targets: 12,
            className: 'aska dt-center',
            visible: false,
          },
        ],
        order: [[ 4, "desc" ]],
        select: {
          style: 'single',
        },
        paging: false,
        buttons: [
          {
              text: 'Toggle <b>GuideToPharma</b> / Asuka Inoue',
              className: 'toggle-source',
              action: function () {
                let gtopText = "<span>Toggle <b>GuideToPharma</b> / Asuka Inoue</span>"
                let askaText = "<span>Toggle GuideToPharma / <b>Asuka Inoue</b></span>"
                let currText = $('.dt-button.toggle-source').html()

                if (currText == gtopText) {
                  $('.dt-button.toggle-source').html(askaText)
                } else if (currText == askaText) {
                  $('.dt-button.toggle-source').html(gtopText)
                }

                var columns = sigmatch_table.columns('.gtop');
                columns.visible(!columns.visible()[0])

                columns = sigmatch_table.columns('.aska');
                columns.visible(!columns.visible()[0])
              }
          },
          {
            text: "Export to Excel",
            action: function() {
              tableToExcel('sigmatch_table', 'Signature Match data', 'SignatureMatch_coupling.xls')
            }
          },
          {
            text: 'Deselect',
            action: function () {
              sigmatch_table.rows().deselect();
            }
          },
        ]
      })

      sigmatch_table.on('select', function (e, dt, type, indexes) {
        const entry = sigmatch_table.rows(indexes).data().toArray()[0];
        console.log(entry);
        svg = d3.select('svg.svg-content.conseq');
        $("g#sigmatch_frame").remove();
        signprotmat.d3.draw_seq_cons(
          entry,
          svg,
          xScale,
          xAxis,
          true
        );
      })

    },
    error: function(error){
      console.log(error);
      document.querySelector("#sigm_spin").style.display = "none";
    },
    complete: function(){
      document.querySelector("#sigm_spin").style.display = "none";
    }
  })
}


const replace_filter_value = function(d) {
  const num = parseInt($('#currentpairs').text())
  if (d === 'dec') {
    num > 1 ? $('#currentpairs').text(num - 1) : null
  } else if (d === 'inc') {
    $('#currentpairs').text(num + 1)
  }
};

const filter_pairs = function(floor, ceiling) {
  // const num = parseInt($('#currentpairs').text())
  d3.select('g#interact')
    .selectAll("rect")
    .style('display', function(d){
      return (floor <= d.pairs.length && ceiling >= d.pairs.length ? 'block' : 'none')
    })
}


const initialize_filter_slider = function() {
  // initializing range slider
  $( "#slider-range" ).slider({
    range: true,
    min: 1,
    max: 20,
    values: [ 1, 20 ],
    slide: function( event, ui ) {
      const floor = parseInt(ui.values[0])
      const ceil = parseInt(ui.values[1])
      $( "#amount" ).val( 'From: ' + floor + " To: " + ceil );
      filter_pairs(floor, ceil)
    }
  });
};

const reset_slider = function() {
  // reset the slider to the possible min and max values
  const min = $( "#slider-range" ).slider( "option", "min" );
  const max = $( "#slider-range" ).slider( "option", "max" );
  $( "#slider-range" ).slider( "values", [ min, max ] );
}

const set_slider_max_value = function(){
  // setting the max value of the interaction filter slider according to selected data
  const max_val = parseInt($("#interface-count").text().match(/\d+/))
  $( "#slider-range" ).slider( "option", "max", max_val );
}

const update_slider_label = function() {
    // update the text span above the slider
  $( "#amount" ).val( 'From: ' + $( "#slider-range" ).slider( "values", 0 ) +
    " To: " + $( "#slider-range" ).slider( "values", 1 ) );
}

var tableToExcel = (function () {
  var uri = 'data:application/vnd.ms-excel;base64,',
    template = '<html xmlns:o="urn:schemas-microsoft-com:office:office" xmlns:x="urn:schemas-microsoft-com:office:excel" xmlns="http://www.w3.org/TR/REC-html40"><head><!--[if gte mso 9]><xml><x:ExcelWorkbook><x:ExcelWorksheets><x:ExcelWorksheet><x:Name>{worksheet}</x:Name><x:WorksheetOptions><x:DisplayGridlines/></x:WorksheetOptions></x:ExcelWorksheet></x:ExcelWorksheets></x:ExcelWorkbook></xml><![endif]--></head><body><table>{table}</table></body></html>',
    base64 = function (s) {
      return window.btoa(unescape(encodeURIComponent(s)))
    }, format = function (s, c) {
      return s.replace(/{(\w+)}/g, function (m, p) {
        return c[p];
      })
    }
  return function (table, name, filename) {
    var table= $("#"+table).clone();
    $("#excel_table").html(table);
    // Clean up table to remove yadcf stuff
    $("#excel_table thead tr").css('height','');
    $("#excel_table thead th").css('height','');
    $("#excel_table thead div").css('height','');
    $("#excel_table thead .yadcf-filter-wrapper").remove();
    $("#excel_table thead button").remove();
    var tr = $("#excel_table thead tr:eq(1)");
    // reattach th titles
    tr.find('th').each (function( column, th) {
      if ($(th).attr('title')) $(th).html($(th).attr('title'));
    });     

    var ctx = {
      worksheet: name || 'Worksheet',
      table: $("#excel_table").html()
    }
    $("#excel_table").html("");
    document.getElementById("dlink").href = uri + base64(format(template, ctx));
    document.getElementById("dlink").download = filename;
    document.getElementById("dlink").click();
  }
})()


$(document).ready(function () {
  non_interactions = signprotmat.data.annotateNonInteractionData(interactions_metadata, non_interactions);

  $('[data-toggle="tooltip"]').tooltip();

  const table = $('#table-interface').DataTable({
    dom: 'Bfrtip',
    data: interactions_metadata,
    columnDefs: [
      {
        data: null,
        targets: 0,
        defaultContent: '',
        orderable: false,
        className: 'select-checkbox',
      }, {
        data: 'name',
        title: 'Name',
        targets: 1,
      }, {
        data: 'family',
        title: 'Family',
        targets: 2,
      }, {
        data: 'class',
        title: 'Class',
        targets: 3,
      }, {
        data: 'organism',
        title: 'Organism',
        targets: 4,
      }, {
        data: 'pdb_id',
        title: 'PDB',
        targets: 5,
      },
      {
        data: 'gprot',
        title: 'G-Protein',
        targets: 6,
      },
      {
        data: 'gprot_class',
        title: 'G-Protein Class',
        targets: 7,
      }],
    order: [[ 7, "desc" ],[ 6, "asc" ],[ 3, "asc" ]],
    select: {
      style: 'os',
    },
    paging: false,
    scrollY: '60vh',
    scrollCollapse: true,
    buttons: [
      {
        text: 'Select all',
        action: function () {
          table.rows().select();
        }
      },
      {
        text: 'Select none',
        action: function () {
          table.rows().deselect();
        }
      },
    ]
  });

  table
    .on('select', function (e, dt, type, indexes) {
      const row_count = table.rows({ selected: true }).count()
      if (row_count >= 2 || row_count === 0) {
        $('#interface-count').text(row_count + ' interfaces selected.');
      } else {
        $('#interface-count').text(row_count + ' interface selected.');
      }
    })
    .on('deselect', function (e, dt, type, indexes) {
      const row_count = table.rows({ selected: true }).count()
      if (row_count >= 2 || row_count === 0) {
        $('#interface-count').text(row_count + ' interfaces selected.');
      } else {
        $('#interface-count').text(row_count + ' interface selected.');
      }
    });

  // Default: Select all rows on page load
  table.rows().select();
  let selection = table.rows({ selected: true }).data();
  pdb_sel = signprotmat.data.select_by_value(selection, 'pdb_id');
  pos_set = signprotmat.data.select_by_value(selection, 'entry_name')

  initialize_filter_slider()
  set_slider_max_value()
  update_slider_label()

  $('#interface-modal-table').on('hidden.bs.modal', function (e) {
    selection = table.rows({ selected: true }).data();
    let old_pdb_sel = pdb_sel;
    pdb_sel = signprotmat.data.select_by_value(selection, 'pdb_id')
    pos_set = signprotmat.data.select_by_value(selection, 'entry_name')

    if (!_.isEqual(old_pdb_sel.sort(), pdb_sel.sort())){
      $('.svg-content').remove();
      con_seq = {};
      document.querySelector('#seqsig-container').style.display = "none";
      document.querySelector('#conseq-container').style.display = "none";
      document.getElementById("interface-svg").className = "collapse in";

      data = signprotmat.data.dataTransformationWrapper(interactions, keys, pdb_sel);
      svg = signprotmat.d3.setup("div#interface-svg");
      xScale = signprotmat.d3.xScale(data.transformed);
      yScale = signprotmat.d3.yScale(data.transformed, gprot);
      xAxis = signprotmat.d3.xAxis(xScale);
      yAxis = signprotmat.d3.yAxis(yScale);
      xAxisGrid = signprotmat.d3.xAxisGrid(xScale, yScale);
      yAxisGrid = signprotmat.d3.yAxisGrid(xScale, yScale);
      pdbScale = signprotmat.d3.pdbScale(data.transformed, interactions_metadata);
      sigScale = signprotmat.d3.sigScale(data.transformed, interactions_metadata);
      colScale = signprotmat.d3.colScale(data.inttypes);
      tooltip = signprotmat.d3.tooltip(svg);
      signprotmat.d3.renderData(
        svg,
        data,
        non_interactions,
        interactions_metadata,
        xScale,
        yScale,
        xAxis,
        yAxis,
        xAxisGrid,
        yAxisGrid,
        colScale,
        pdbScale,
        sigScale,
        tooltip
      );
      document.querySelector("#intbut").classList.add("active");
      document.querySelector("#resbut").classList.remove("active");

      set_slider_max_value()   
      reset_slider()
      update_slider_label()

      // interface_data_table.clear();
      // interface_data_table.rows.add(data.transformed);
      // interface_data_table.draw();

      // const row_selection = pdb_table.rows({ selected: true }).data();
      // const xvals = xScale.domain();
      // const prids = signprotmat.data.select_by_value(row_selection, 'entry_name');

      // receptor_data = signprotmat.data.get_additional_receptors(rs, xvals, prids)
      // console.log(receptor_data)
      // signprotmat.d3.addReceptor(receptor_data, data, svg);
    };
  });

  let keys = [
    "rec_chain",
    "rec_aa",
    "rec_pos",
    "rec_gn",
    "sig_chain",
    "sig_aa",
    "sig_pos",
    "sig_gn",
    "int_ty",
    "gprot",
    "entry_name",
    "pdb_id"
  ];

  data = signprotmat.data.dataTransformationWrapper(interactions, keys, pdb_sel);
  svg = signprotmat.d3.setup("div#interface-svg");
  xScale = signprotmat.d3.xScale(data.transformed);
  yScale = signprotmat.d3.yScale(data.transformed, gprot);
  xAxis = signprotmat.d3.xAxis(xScale);
  yAxis = signprotmat.d3.yAxis(yScale);
  xAxisGrid = signprotmat.d3.xAxisGrid(xScale, yScale);
  yAxisGrid = signprotmat.d3.yAxisGrid(xScale, yScale);
  pdbScale = signprotmat.d3.pdbScale(data.transformed, interactions_metadata);
  sigScale = signprotmat.d3.sigScale(data.transformed, interactions_metadata);
  colScale = signprotmat.d3.colScale(data.inttypes);
  tooltip = signprotmat.d3.tooltip(svg);
  signprotmat.d3.renderData(
    svg,
    data,
    non_interactions,
    interactions_metadata,
    xScale,
    yScale,
    xAxis,
    yAxis,
    xAxisGrid,
    yAxisGrid,
    colScale,
    pdbScale,
    sigScale,
    tooltip
  );

  // PDB TABLE MODAL
  // const pdb_table = $('#table-pdb').DataTable({
  //     dom: 'Bfrtip',
  //     data: ps,
  //     scrollX: true,
  //     columnDefs: [
  //         {
  //             data: null,
  //             targets: 0,
  //             defaultContent: '',
  //             orderable: false,
  //             className: 'select-checkbox',
  //         }, {
  //             data: 'name',
  //             title: 'Name',
  //             targets: 1,
  //         }, {
  //             data: 'protein_class',
  //             title: 'Class',
  //             targets: 2,
  //         }, {
  //             data: 'protein_family',
  //             title: 'Family',
  //             targets: 3,
  //         }, {
  //             data: 'ligand',
  //             title: 'Ligand',
  //             targets: 4,
  //         }],
  //     select: {
  //         style: 'os',
  //     },
  //     paging: true,
  //     buttons: [
  //         {
  //             text: 'Select all',
  //             action: function () {
  //                 pdb_table.rows().select();
  //             }
  //         },
  //         {
  //             text: 'Select none',
  //             action: function () {
  //                 pdb_table.rows().deselect();
  //             }
  //         },
  //     ]
  // });


  // pdb_table
  //     .on('select', function (e, dt, type, indexes) {
  //         const row_count = pdb_table.rows({ selected: true }).count()
  //         if (row_count >= 2 || row_count === 0) {
  //             $('#pdb-count').text(row_count + ' PDBs selected.')
  //         } else {
  //             $('#pdb-count').text(row_count + ' PDB selected.')
  //         }
  //     })
  //     .on('deselect', function (e, dt, type, indexes) {
  //         const row_count = pdb_table.rows({ selected: true }).count()
  //         if (row_count >= 2 || row_count === 0) {
  //             $('#pdb-count').text(row_count + ' PDBs selected.')
  //         } else {
  //             $('#pdb-count').text(row_count + ' PDB selected.')
  //         }
  //     });

  // $('#pdb-modal-table').on('hidden.bs.modal', function (e) {
  //     const row_selection = pdb_table.rows({ selected: true }).data();
  //     const xvals = xScale.domain();
  //     const prids = signprotmat.data.select_by_value(row_selection, 'entry_name');

  //     neg_set = signprotmat.data.select_by_value(row_selection, 'entry_name');

  //     console.log(xvals)
  //     console.log(prids)
  //     receptor_data = signprotmat.data.get_additional_receptors(rs, xvals, prids)
  //     console.log(receptor_data)
  //     signprotmat.d3.addReceptor(receptor_data, data, svg);
  // });

  // let interface_data_table = $('#table-data').DataTable({
  //     data: data.transformed,
  //     columns: [
  //         { data: 'int_ty' },
  //         { data: 'pdb_id' },
  //         { data: 'rec_aa' },
  //         { data: 'rec_gn' },
  //         { data: 'rec_pos' },
  //         { data: 'sig_aa' },
  //         { data: 'sig_gn' }
  //     ]
  // });

  // https://www.datatables.net/examples/api/tabs_and_scrolling.html
  $(document).on('shown.bs.modal', function (e) {
    $.fn.dataTable.tables({ visible: true, api: true }).columns.adjust();
  });

});

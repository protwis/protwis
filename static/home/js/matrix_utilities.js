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
    acc[current.rec_gn] = [current.pdb_id].concat(acc[current.rec_gn])
    return acc;
  }, {});

  // reformat the non interacting data in the same way
  // but also exclude all proteins that already
  // have an interaction registered in the data_int object
  const ignore_markers = data_non.reduce(function (acc, current) {
    let protein_array = data_int[current.rec_gn]
    if ( !protein_array.includes(current.pdb_id) ) {
      acc[current.rec_gn] = [current.pdb_id.toLowerCase()].concat(acc[current.rec_gn])
    }
    return acc;
  }, {});

  // example return
  // ignore_markers = {
  //    '2.39x39': [ "oprm_mouse", "oprm_mouse" ],
  //}

  return ignore_markers;
}

const get_receptor_classes = function(receptor_metadata, pdb_id_array){
  return receptor_metadata.filter(x => pdb_id_array.includes(x.pdb_id)).map(x => x.class)
}

const run_seq_sig = function(){
  let segments = get_gn();
  let ignore_markers = get_ignore();
  let selected_receptor_classes = get_receptor_classes(interactions_metadata, pdb_sel);
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
      selectedreceptorclasses: selected_receptor_classes,
      ignore: JSON.stringify(ignore_markers),
    },
    beforeSend: function(){
      old_sets = pos_set;
      // document.querySelector("#calc_spin").style.display = "inline-block";
      // $("#calc_spin").addClass("fa-spin");
      // $("#calc_spin").addClass("fa-spinner");
      // $("#calc_spin").removeClass("fa-times");
      // $("#sigm_spin").addClass("fa-spin");
      // $("#sigm_spin").addClass("fa-spinner");
      // $("#sigm_spin").removeClass("fa-times");
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
      // svg = signprotmat.d3.setup("div#conseq-svg", 'conseq');
      // signprotmat.d3.draw_seq_cons(
      //   data,
      //   svg,
      //   xScale,
      //   xAxis
      // );
      initialize_consensus(data.feat);

      // Once done run the signature match
      run_sig_match();
    },
    error: function(error){
      // $("#calc_spin").addClass("fa-times");
      // $("#calc_spin").removeClass("fa-spinner");
      // $("#calc_spin").removeClass("fa-spin");
      // $("#sigm_spin").addClass("fa-times");
      // $("#sigm_spin").removeClass("fa-spinner");
      // $("#sigm_spin").removeClass("fa-spin");
      console.log(error)
      alert(error);
    },
    complete: function(){
      // $("#calc_spin").addClass("fa-check");
      // $("#calc_spin").removeClass("fa-spinner");
      // $("#calc_spin").removeClass("fa-spin");
    }
  });
};

const initialize_consensus = function(cons_data){
  const row_height = 30;
  const AMINO_ACID_GROUPS = {"HY_any": ["A", "C", "F", "I", "L", "M", "P", "V", "W", "Y"], "HY_4-5": ["F", "M", "Y"], "HA_any": ["A", "I", "L", "M", "V"], "HA_1-2": ["A", "V"], "HA_2-3": ["I", "L", "V"], "HA_3": ["I", "L"], "HA_3-4": ["I", "L", "M"], "M_4": ["M"], "A_1": ["A"], "I_": ["I"], "L_": ["L"], "V_": ["V"], "HR_any": ["F", "H", "W", "Y"], "HR_4-5": ["F", "H", "Y"], "HR_5": ["F", "Y"], "HR_5-6": ["F", "W", "Y"], "W_6": ["W"], "Y_?": ["Y"], "F_": ["F"], "Hb_any": ["D", "E", "H", "K", "N", "Q", "R", "S", "T", "Y"], "Hb_2": ["S", "T"], "Hb_3": ["H", "N"], "N_": ["N"], "Q_": ["Q"], "S_": ["S"], "T_": ["T"], "Hu_2-3": ["N", "S", "T"], "Hu_3-4": ["H", "N", "Q"], "Ha_2-3": ["D", "N", "S", "T"], "Ha_3": ["D", "N"], "Ha_3-4": ["D", "E", "H", "N", "Q"], "Ha_4": ["E", "H", "Q"], "Hd_4": ["H", "Q"], "Hd_4-5": ["H", "K", "Q"], "Hd_5-6": ["K", "R", "Y"], "Hd_6": ["R", "Y"], "+-_any": ["D", "E", "H", "K", "R"], "+-_3-4": ["D", "E", "H"], "+-_4-5": ["E", "H", "K"], "+_any": ["H", "K", "R"], "H_4": ["H"], "+_4-5": ["H", "K"], "K_5": ["K"], "+_5-6": ["K", "R"], "R_6": ["R"], "-_any": ["D", "E"], "D_3": ["D"], "E_4": ["E"], "Sm_0-2": ["A", "C", "G", "S"], "Sm_0-1": ["A", "G"], "Sm_1-2": ["A", "C", "S"], "aH_Hig": ["A", "K", "L", "M", "R"], "aH_Low": ["G", "P"], "G_0": ["G"], "P_2": ["P"], "C_2": ["C"], "-": "-"}

  for (let key of Object.keys(cons_data)) {
    con_seq[key] = [cons_data[key][0]];
  };

  non_interactions.push(...data.receptor)
  var seq_cons_data = []
  for (gn of xScale.domain()){
    if (Object.keys(con_seq).includes(gn)){
      var aa = con_seq[gn][0]['aa']
      var seq_cons = parseInt(
        (_.uniqBy(non_interactions.filter(x => x.rec_gn === gn), "pdb_id")
          .filter(x => pdbScale.domain().includes(x.pdb_id))
          .filter(x => x.rec_aa === aa).length /
          pdbScale.domain().length) *
          100
      );

      var sort_code = con_seq[gn][0]['sort_code'].replace('α', 'a')
      var gn_aa = _.uniqBy(non_interactions.filter(x => x.rec_gn === gn), "pdb_id").filter(x => pdbScale.domain().includes(x.pdb_id)).map(x => x.rec_aa)
      prop_seq_cons = 0
      ignore_sort_code = ['-_', 'Sm_any']
      for (var aa of gn_aa){
        if (!ignore_sort_code.includes(sort_code) && AMINO_ACID_GROUPS[sort_code].includes(aa)){
          prop_seq_cons += 1
        }
      }
      prop_seq_cons = parseInt(prop_seq_cons / pdbScale.domain().length * 100)

      con_seq[gn][0]['seq_cons'] = seq_cons
      con_seq[gn][0]['prop_seq_cons'] = prop_seq_cons
    }
  }
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
      // document.querySelector("#sigm_spin").style.display = "inline-block";
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
            visible: false,
          }, {
            data: 'class',
            title: 'Class',
            targets: 1,
          }, {
            data: 'family',
            title: 'Family',
            targets: 2,
          }, {
            data: 'subfamily',
            title: 'Sub-Family',
            targets: 3,
            visible: false,
          }, {
            data: 'prot',
            title: 'IUPHAR',
            targets: 4,
          }, {
            data: 'nscore',
            title: 'Interface Conservation (%)',
            targets: 6,
          }, {
            data: 'entry',
            title: 'UniProt',
            targets: 5,
          }, {
            data: 'GuideToPharma.Gs.html',
            title: '   Gs   ',
            targets: 7,
            className: 'gtop dt-center',
            visible: true,
          }, {
            data: 'GuideToPharma.Gi/Go.html',
            title: 'Gi / Go ',
            targets: 8,
            className: 'gtop dt-center',
            visible: true,
          }, {
            data: 'GuideToPharma.Gq/G11.html',
            title: 'Gq / G11',
            targets: 9,
            className: 'gtop dt-center',
            visible: true,
          }, {
            data: 'GuideToPharma.G12/G13.html',
            title: 'G12 / G13',
            targets: 10,
            className: 'gtop dt-center',
            visible: true,
          }, {
            data: 'Aska.Gs.html',
            title: '   Gs   ',
            targets: 11,
            className: 'aska dt-center',
            visible: false,
          }, {
            data: 'Aska.Gi/Go.html',
            title: 'Gi / Go ',
            targets: 12,
            className: 'aska dt-center',
            visible: false,
          }, {
            data: 'Aska.Gq/G11.html',
            title: 'Gq / G11',
            targets: 13,
            className: 'aska dt-center',
            visible: false,
          }, {
            data: 'Aska.G12/G13.html',
            title: 'G12 / G13',
            targets: 14,
            className: 'aska dt-center',
            visible: false,
          }, {
            data: 'Merged.Gs.html',
            title: '   Gs   ',
            targets: 15,
            className: 'merg dt-center',
            visible: false,
          }, {
            data: 'Merged.Gi/Go.html',
            title: 'Gi / Go ',
            targets: 16,
            className: 'merg dt-center',
            visible: false,
          }, {
            data: 'Merged.Gq/G11.html',
            title: 'Gq / G11',
            targets: 17,
            className: 'merg dt-center',
            visible: false,
          }, {
            data: 'Merged.G12/G13.html',
            title: 'G12 / G13',
            targets: 18,
            className: 'merg dt-center',
            visible: false,
          },
        ],
        order: [[ 6, "desc" ]],
        select: {
          style: 'single',
        },
        paging: false,
        buttons: [
          {
              text: 'Toggle <b>GuideToPharma</b> / Asuka Inoue / Merged Data',
              className: 'toggle-source',
              action: function () {
                let gtopText = "<span>Toggle <b>GuideToPharma</b> / Asuka Inoue / Merged Data</span>"
                let askaText = "<span>Toggle GuideToPharma / <b>Asuka Inoue</b> / Merged Data</span>"
                let mergText = "<span>Toggle GuideToPharma / Asuka Inoue / <b>Merged Data</b></span>"
                let currText = $('.dt-button.toggle-source').html()

                if (currText == gtopText) {
                  $('.dt-button.toggle-source').html(askaText)
                    var columns = sigmatch_table.columns('.gtop');
                    columns.visible(!columns.visible()[0])

                    columns = sigmatch_table.columns('.aska');
                    columns.visible(!columns.visible()[0])
                } else if (currText == askaText) {
                  $('.dt-button.toggle-source').html(mergText)
                    var columns = sigmatch_table.columns('.aska');
                    columns.visible(!columns.visible()[0])

                    columns = sigmatch_table.columns('.merg');
                    columns.visible(!columns.visible()[0])
                } else if (currText == mergText) {
                  $('.dt-button.toggle-source').html(gtopText)
                    var columns = sigmatch_table.columns('.merg');
                    columns.visible(!columns.visible()[0])

                    columns = sigmatch_table.columns('.gtop');
                    columns.visible(!columns.visible()[0])
                }

              }
          },
          {
            text: "Export to Excel",
            action: function() {
              tableToExcel('sigmatch_table', 'Signature Match data', 'SignatureMatch_coupling.xls')
            }
          },
          {
            text: "Export to CSV",
            action: function ( e, dt, button, config ) {
              var table_data = sigmatch_table.data().toArray();

              var export_data = []
              for (let i of Object.values(table_data)){
                let r = {}
                r['name'] = i['entry']
                r['family'] = i['family']
                r['subfamily'] = i['subfamily']
                r['score'] = i['nscore']
                r['aska_gs'] = i['Aska']['Gs']['text']
                r['aska_gio'] = i['Aska']['Gi/Go']['text']
                r['aska_gq11'] = i['Aska']['Gq/G11']['text']
                r['aska_g1213'] = i['Aska']['G12/G13']['text']
                r['gtop_gs'] = i['GuideToPharma']['Gs']['text']
                r['gtop_gio'] = i['GuideToPharma']['Gi/Go']['text']
                r['gtop_gq11'] = i['GuideToPharma']['Gq/G11']['text']
                r['gtop_g1213'] = i['GuideToPharma']['G12/G13']['text']
                r['merg_gs'] = i['Merged']['Gs']['text']
                r['merg_gio'] = i['Merged']['Gi/Go']['text']
                r['merg_gq11'] = i['Merged']['Gq/G11']['text']
                r['merg_g1213'] = i['Merged']['G12/G13']['text']
                export_data.push(r)
              }

              export_data = Papa.unparse(export_data)
              
              $('<a></a>')
                .attr('id','downloadFile')
                .attr('href','data:text/csv;charset=utf8,' + encodeURIComponent(export_data))
                .attr('download', 'export.csv')
                .appendTo('body');
              
              $('#downloadFile').ready(function() {
                $('#downloadFile').get(0).click();
                $('#downloadFile').remove();
              });
            }
          },
          {
            text: 'Reset All Filters',
            action: function () {
              yadcf.exResetAllFilters(sigmatch_table)
            }
          },
          {
            text: 'Show Alignment',
            className: 'score-button',
            action: function () {
              sigmatch_table.rows().deselect();
            }
          },
        ]
      })

      text_col_filter = {
        filter_type: "multi_select",
        select_type: 'select2',
        filter_reset_button_text: false,
      }

      range_col_filter = {
        filter_type: "range_number_slider",
        filter_delay: 70,
        filter_reset_button_text: false,
      }

      coupl_col_filter = {
        filter_type: "multi_select",
        select_type: 'select2',
        filter_reset_button_text: false,
        column_data_type: "html",
        html_data_type: "text",
      }

      yadcf.init(sigmatch_table, [
        {column_number : 2, ...text_col_filter},
        {column_number : 4, ...text_col_filter},
        {column_number : 6, ...range_col_filter},
        {column_number : 7, ...coupl_col_filter},
        {column_number : 8, ...coupl_col_filter},
        {column_number : 9, ...coupl_col_filter},
        {column_number : 10, ...coupl_col_filter},
        {column_number : 11, ...coupl_col_filter},
        {column_number : 12, ...coupl_col_filter},
        {column_number : 13, ...coupl_col_filter},
        {column_number : 14, ...coupl_col_filter},
        {column_number : 15, ...coupl_col_filter},
        {column_number : 16, ...coupl_col_filter},
        {column_number : 17, ...coupl_col_filter},
        {column_number : 18, ...coupl_col_filter},
      ]);

      $('.score-button').click( function () {
          const render_url = window.location.origin + '/signprot/matrix/render_sigmat/';
          window.open(render_url,'_blank');
      });

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
      // $("#sigm_spin").addClass("fa-times");
      // $("#sigm_spin").removeClass("fa-spinner");
      // $("#sigm_spin").removeClass("fa-spin");
      console.log(error)
      alert(error);
    },
    complete: function(){
      // $("#sigm_spin").addClass("fa-check");
      // $("#sigm_spin").removeClass("fa-spinner");
      // $("#sigm_spin").removeClass("fa-spin");
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
  const max_val = get_max_interface_count()
  d3.select('g#interact')
    .selectAll("g")
    .style('display', function(d){
      var ratio = d.pairs.length / max_val * 100
      return (floor <= ratio && ceiling >= ratio ? 'block' : 'none')
    })
}

const set_min = function () {
  $('#amount_min').html( $('#slider-range').slider('values', 0)).position({
      my: 'center bottom',
      at: 'center top',
      of: $('#slider-range span:eq(0)'),
  });
}

const set_max = function () {
    $('#amount_max').html( $('#slider-range').slider('values', 1)).position({
      my: 'center bottom',
      at: 'center top',
      of: $('#slider-range span:eq(1)'),
  });
}

const initialize_filter_slider = function() {
  // initializing range slider
  $( "#slider-range" ).slider({
    range: true,
    min: 0,
    max: 100,
    step: 1,
    values: [ 0, 100 ],
    slide: function( event, ui ) {
      const floor = parseInt(ui.values[0])
      const ceil = parseInt(ui.values[1])

      var label = ui.handleIndex == 0 ? '#amount_min' : '#amount_max';
      $(label).html(ui.value).position({
          my: 'center bottom',
          at: 'center top',
          of: ui.handle,
      });

      filter_pairs(floor, ceil)
    }
  });

  setTimeout(set_min, setTimeout(set_max, 1), 1)
};

const reset_slider = function() {
  // reset the slider to the possible min and max values
  const min = $( "#slider-range" ).slider( "option", "min" );
  const max = $( "#slider-range" ).slider( "option", "max" );
  $( "#slider-range" ).slider( "values", [ min, max ] );

  setTimeout(set_min, setTimeout(set_max, 1), 1)
}

const get_max_interface_count = function() {
  // return parseInt($("#single-crystal-group-pdbs-modal-text").text().match(/\d+/))
  return pos_set.length
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
  $.get('/signprot/pdbtabledata', { exclude_non_interacting: true}, function(data) {
    $('#interface-modal-table .tableview').html(data);
  })

  signprotmat.d3.print_resScaleColor_legend()
  $('[data-toggle="tooltip"]').tooltip();
  
  let table = $($.fn.dataTable.tables()[0]).DataTable();

  initialize_filter_slider();

  $('#interface-modal-table').on('shown.bs.modal', function(e) {
    showPDBtable('#interface-modal-table');
  })

  $('#interface-modal-table').on('hidden.bs.modal', function (e) {
    table = $($.fn.dataTable.tables()[0]).DataTable();
    selection = table.rows('.selected').data();
   
    let old_pdb_sel = pdb_sel;
    pdb_sel = [];
    pos_set = [];

    // get selected pdb ids
    $('.pdb_selected:checked').each(function( index ) {pdb_sel.push(($( this ).attr('id')))});
  
    // get corresponding protein entry_name values
    for (var int_meta of interactions_metadata){
      if (pdb_sel.indexOf(int_meta['pdb_id']) != -1){
        pos_set = [int_meta['conf_id'], ...pos_set]
      }
    }

    if (!_.isEqual(old_pdb_sel.sort(), pdb_sel.sort())){
      $('.svg-content').remove();
      con_seq = {};
      document.querySelector('#seqsig-container').style.display = "none";
      document.querySelector('#conseq-container').style.display = "none";
      document.getElementById("interface-svg").className = "collapse in";

      data = signprotmat.data.dataTransformationWrapper(interactions, pdb_sel);
      svg = signprotmat.d3.setup("div#interface-svg");
      xScale = signprotmat.d3.xScale(data.transformed, receptor);
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

      reset_slider();
      run_seq_sig();
    };
  });
});

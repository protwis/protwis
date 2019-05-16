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
            data: 'Gs',
            title: '   Gs   ',
            targets: 5,
            className: 'dt-center',
          }, {
            data: 'Gi/Go',
            title: 'Gi / Go ',
            targets: 6,
            className: 'dt-center',
          }, {
            data: 'Gq/G11',
            title: 'Gq / G11',
            targets: 7,
            className: 'dt-center',
          }, {
            data: 'G12/G13',
            title: 'G12 / G13',
            targets: 8,
            className: 'dt-center',
          },
        ],
        order: [[ 4, "desc" ]],
        select: {
          style: 'single',
        },
        paging: false,
        buttons: [
          // {
          //   text: 'Select all',
          //   action: function () {
          //     sigmatch_table.rows().select();
          //   }
          // },
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

const filter_pairs = function() {
  const num = parseInt($('#currentpairs').text())
  d3.select('g#interact')
    .selectAll("rect")
    .style('display', function(d){return (num <= d.pairs.length ? 'block' : 'none') })
}


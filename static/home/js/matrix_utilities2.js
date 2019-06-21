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

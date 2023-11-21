/*global Papa,d3,yadcf,createYADCFfilters,GlobalTableToExcel,showPDBtable,showAlert, _, signprotmat, con_seq, non_interactions,interactions_metadata,interactions,gprot,receptor,csrf_token*/
/*eslint complexity: ["error", 20]*/

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

var selected_interactions;
var non_interactions;

let old_sets = [];
let pos_set = [];
let neg_set = [];

var filtering_particle = false;
try {
  filtering_particle = document.getElementById("htmlVariable").getAttribute( "data-url" );
} catch (e) {
  filtering_particle = false;
}

const get_gn = function() {
  const segments = [];
  const orig = [];
  const patt = /(^\d{1,})|(x\d{2,})/g;

  const selection = data.receptor;
  const value = "rec_gn";

  for (let index = 0; index < selection.length; index++) {
    let comb_gn = "";
    let curr = selection[parseInt(index)][String(value)];
    let match = patt.exec(curr);
    while (match !== null) {
      comb_gn += match[0];
      match = patt.exec(selection[index][value]);
    }
    if (comb_gn.length > 1) {
      segments.push(comb_gn);
      orig.push(curr);
    }
  }
  return {
    orig: _.uniq(orig),
    label: _.uniq(segments)
  };
};

const get_ignore = function() {
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
  const data_int = data.receptor.reduce(function(acc, current) {
    acc[current.rec_gn] = [current.pdb_id].concat(acc[current.rec_gn]);
    return acc;
  }, {});

  // reformat the non interacting data in the same way
  // but also exclude all proteins that already
  // have an interaction registered in the data_int object
  const ignore_markers = data_non.reduce(function(acc, current) {
    let protein_array = data_int[current.rec_gn];
    if (!protein_array.includes(current.pdb_id)) {
      acc[current.rec_gn] = [current.pdb_id.toLowerCase()].concat(acc[current.rec_gn]);
    }
    return acc;
  }, {});

  // example return
  // ignore_markers = {
  //    "2.39x39": [ "oprm_mouse", "oprm_mouse" ],
  //}

  return ignore_markers;
};

const get_receptor_classes = function(receptor_metadata, pdb_id_array) {
  return receptor_metadata.filter((x) => pdb_id_array.includes(x.pdb_id)).map((x) => x.class);
}

const initialize_consensus = function(cons_data) {
  const row_height = 30;
  const AMINO_ACID_GROUPS = {
    "HY_any": ["A", "C", "F", "I", "L", "M", "P", "V", "W", "Y"],
    "HY_4-5": ["F", "M", "Y"],
    "HA_any": ["A", "I", "L", "M", "V"],
    "HA_1-2": ["A", "V"],
    "HA_2-3": ["I", "L", "V"],
    "HA_3": ["I", "L"],
    "HA_3-4": ["I", "L", "M"],
    "M_4": ["M"],
    "A_1": ["A"],
    "I_": ["I"],
    "L_": ["L"],
    "V_": ["V"],
    "HR_any": ["F", "H", "W", "Y"],
    "HR_4-5": ["F", "H", "Y"],
    "HR_5": ["F", "Y"],
    "HR_5-6": ["F", "W", "Y"],
    "W_6": ["W"],
    "Y_?": ["Y"],
    "F_": ["F"],
    "Hb_any": ["D", "E", "H", "K", "N", "Q", "R", "S", "T", "Y"],
    "Hb_2": ["S", "T"],
    "Hb_3": ["H", "N"],
    "N_": ["N"],
    "Q_": ["Q"],
    "S_": ["S"],
    "T_": ["T"],
    "Hu_2-3": ["N", "S", "T"],
    "Hu_3-4": ["H", "N", "Q"],
    "Ha_2-3": ["D", "N", "S", "T"],
    "Ha_3": ["D", "N"],
    "Ha_3-4": ["D", "E", "H", "N", "Q"],
    "Ha_4": ["E", "H", "Q"],
    "Hd_4": ["H", "Q"],
    "Hd_4-5": ["H", "K", "Q"],
    "Hd_5-6": ["K", "R", "Y"],
    "Hd_6": ["R", "Y"],
    "+-_any": ["D", "E", "H", "K", "R"],
    "+-_3-4": ["D", "E", "H"],
    "+-_4-5": ["E", "H", "K"],
    "+_any": ["H", "K", "R"],
    "H_4": ["H"],
    "+_4-5": ["H", "K"],
    "K_5": ["K"],
    "+_5-6": ["K", "R"],
    "R_6": ["R"],
    "-_any": ["D", "E"],
    "D_3": ["D"],
    "E_4": ["E"],
    "Sm_0-2": ["A", "C", "G", "S"],
    "Sm_0-1": ["A", "G"],
    "Sm_1-2": ["A", "C", "S"],
    "aH_Hig": ["A", "K", "L", "M", "R"],
    "aH_Low": ["G", "P"],
    "G_0": ["G"],
    "P_2": ["P"],
    "C_2": ["C"],
    "-": "-"
  }

  for (let key of Object.keys(cons_data)) {
    con_seq[key] = [cons_data[key][0]];
  };

  non_interactions.push(...data.receptor);
  var gn;
  var seq_cons_data = [];
  for (gn of xScale.domain()) {
    if (Object.keys(con_seq).includes(gn)) {
      var aa = con_seq[String(gn)][0]["aa"];
      var seq_cons = parseInt(
        (_.uniqBy(non_interactions.filter((x) => {
            x.rec_gn === String(gn);
          }), "pdb_id")
          .filter((x) => {
            pdbScale.domain().includes(x.pdb_id);
          })
          .filter((x) => {
            x.rec_aa === aa;
          }).length /
          pdbScale.domain().length) *
        100
      );

      var sort_code = con_seq[String(gn)][0]["sort_code"].replace("α", "a");
      var gn_aa = _.uniqBy(non_interactions.filter((x) => {
        x.rec_gn === String(gn);
      }), "pdb_id").filter((x) => {
        pdbScale.domain().includes(x.pdb_id);
      }).map((x) => {
        x.rec_aa;
      });
      var prop_seq_cons = 0;
      var ignore_sort_code = ["-_", "Sm_any"];
      for (var aa_item of gn_aa) {
        if (!ignore_sort_code.includes(sort_code) && AMINO_ACID_GROUPS[String(sort_code)].includes(aa_item)) {
          prop_seq_cons += 1;
        }
      }
      prop_seq_cons = parseInt(prop_seq_cons / pdbScale.domain().length * 100);

      con_seq[String(gn)][0]["seq_cons"] = seq_cons;
      con_seq[String(gn)][0]["prop_seq_cons"] = prop_seq_cons;
    }
  }
  signprotmat.d3.conSeqUpdate(row_height);
};

const run_seq_sig = function(interface_data) {
  let segments = get_gn();
  let ignore_markers = get_ignore();
  var selected_receptor_classes = get_receptor_classes(interactions_metadata, pdb_sel);

  if (pos_set.length < 1) {
    showAlert("No valid set selected.", "warning");
    return;
  };

  if (_.isEqual(old_sets.sort(), pos_set.sort())) {
    showAlert("The selected set is identical to the previously selected one.", "warning");
    return;
  };

  // api request
  let req = $.ajax({
    type: "POST",
    url: "/signprot/matrix/seqsig/",
    data: {
      csrfmiddlewaretoken: csrf_token,
      pos: pos_set,
      seg: segments.label,
      selectedreceptorclasses: selected_receptor_classes,
      ignore: JSON.stringify(ignore_markers),
    },
    beforeSend() {
      old_sets = pos_set;
    },
    success(data) {
      $(".svg-content.seqsig").remove();
      $(".svg-content.conseq").remove();

      // sorting data correctly, first by frequency; then by length
      let tmp_data = {};
      for (let rec_gn in data.feat) {
        let position = _.orderBy(
          data.feat[rec_gn],
          ["freq", "sort_score"],
          ["desc", "asc"]
        );
        tmp_data[rec_gn] = position;
      }

      data.feat = tmp_data;
      document.querySelector("#seqsig-container").style.display = "block";
      document.querySelector("#conseq-container").style.display = "block";

      // Calculate interaction conservation
      let int_cons_list = {};
      for (const interaction of interface_data["transformed"]) {
        let struct = interaction.pdb_id;
        let gn = interaction.rec_gn;
        if (typeof(int_cons_list[gn]) == "undefined") {
          int_cons_list[gn] = [struct];
        } else {
          if (int_cons_list[gn].indexOf(struct) === -1) {
            int_cons_list[gn].push(struct);
          }
        }
      }

      // Add interaction frequencies to residue property data
      let struct_count = interface_data["pdbids"].length;
      for (let i in data.feat){
        if (Object.prototype.hasOwnProperty.call(data.feat, i)) {
          for (const obj of data.feat[i]){
            if (typeof(int_cons_list[obj.gn]) == "undefined") {
              obj.int_freq = 0;
            } else {
              obj.int_freq = Math.round(int_cons_list[obj.gn].length/struct_count*100,0);
            }
          }
        }
      }

      // d3 draw
      svg = signprotmat.d3.setup("div#seqsig-svg", "seqsig");
      signprotmat.d3.draw_seq_sig(
        data,
        svg,
        xScale,
      );
      var startTime = performance.now();
      // Once done run the signature match
      run_sig_match();
      var endTime = performance.now();
      var elapsed = endTime - startTime;
      console.log("run_sig_match execution Time: " + elapsed);

      startTime = performance.now();
      initialize_consensus(data.feat);
      endTime = performance.now();
      elapsed = endTime - startTime;
      console.log("initialize_consensus execution Time: " + elapsed);
    },
    error(jqXHR, exception) {
      console.log(jqXHR);
      console.log(exception);
    },
    complete() {}
  });
};

const run_sig_match = function() {
  let cutoff = $("#cutoff-val").val();
  if (cutoff === "") {
    cutoff = 0;
  }

  cutoff = 0;
  let segments = get_gn();

  let req = $.ajax({
    type: "POST",
    url: "/signprot/matrix/sigmat/",
    data: {
      csrfmiddlewaretoken: csrf_token,
      pos: pos_set,
      seg: segments.label,
      filtering_particle,
      cutoff,
    },
    beforeSend() {
      // document.querySelector("#sigm_spin").style.display = "inline-block";
    },
    success(data) {
      // console.log(data);
      document.querySelector("#sigmatch-container").style.display = "inline-block";
      //sigmatch_data = Object.keys(data).map((key) => {data[String(key)];});
      let column_filters = [];
      let columns_to_add = [];
      // Family column
      column_filters = column_filters.concat(createYADCFfilters(2, 1, "multi_select", "select2", "Select", false, null, null, "120px"));
      // IUPHAR column
      column_filters = column_filters.concat(createYADCFfilters(3, 1, "multi_select", "select2", "Select", false, null, null, "80px"));
      // Interface Conservation % column
      column_filters = column_filters.concat(createYADCFfilters(5, 1, "range_number_slider", null, null, false, null, null, "40px"));
      let columns_definition = [{
                data: null,
                targets: 0,
                defaultContent: "",
                orderable: false,
                className: "select-checkbox",
                visible: false,
              }, {
                data: "class",
                title: "Class",
                targets: 1,
              }, {
                data: "family",
                title: "Family",
                targets: 2,
              }, {
                data: "prot",
                title: "IUPHAR",
                targets: 3,
              }, {
                data: "entry",
                title: "UniProt",
                targets: 4,
              }, {
                data: "nscore",
                title: "Interface Conservation (%)",
                targets: 5,
              }];
      sigmatch_data = Object.keys(data).map(key => data[key]);
      if (filtering_particle === "G alpha") {
        columns_to_add = [{
              data: "Gs.html",
              title: "Gs",
              targets: 6,
            }, {
              data: "Gi/o.html",
              title: "Gi/o",
              targets: 7,
            }, {
              data: "Gq/11.html",
              title: "Gq/11",
              targets: 8,
            }, {
              data: "G12/13.html",
              title: "G12/13",
              targets: 9,
            }, {
              data: "Gs_emax.html",
              title: "Gs",
              targets: 10,
            }, {
              data: "Gi/o_emax.html",
              title: "Gi/o",
              targets: 11,
            }, {
              data: "Gq/11_emax.html",
              title: "Gq/11",
              targets: 12,
            }, {
              data: "G12/13_emax.html",
              title: "G12/13",
              targets: 13,
            }, ];
        columns_definition = columns_definition.concat(columns_to_add);
        // Gprots columns (4 GtoP + 4 GPCRdb Mean)
        column_filters = column_filters.concat(createYADCFfilters(6, 8, "range_number", null, ["Min", "Max"], false, null, "html", "40px"));
      } else {
        columns_to_add = [{
              data: "arrb1.html",
              title: "ARRB1",
              targets: 6,
            }, {
              data: "arrb2.html",
              title: "ARRB2",
              targets: 7,
            },];
            // {
            //   data: "arrb1_emax.html",
            //   title: "ARRB1 (GPCRdb Mean)",
            //   targets: 8,
            //   visible: false,
            // }, {
            //   data: "arrb2.html",
            //   title: "ARRB2 (GPCRdb Mean)",
            //   targets: 9,
            //   visible: false,
            // },];
        // Arrestins columns (2 GtoP + 2 GPCRdb Mean)
        column_filters = column_filters.concat(createYADCFfilters(6, 2,"range_number", null, ["Min", "Max"], false, null, "html", "40px"));
      }
      columns_definition = columns_definition.concat(columns_to_add);
      if (filtering_particle === "G alpha") {
        sigmatch_table = $("#sigmatch_table").DataTable({
          dom: "Bfrtip",
          data: sigmatch_data,
          // scrollY: "60vh",
          scrollY: true,
          scrollX: true,
          scrollCollapse: true,
          scroller: true,
          destroy: true,
          columnDefs: columns_definition,
          order: [
            [5, "desc"]
          ],
          select: {
            style: "single",
          },
          paging: true,
          buttons: [
            {
              text: "Export to Excel",
              action() {
                GlobalTableToExcel("sigmatch_table", "Signature Match data", "SignatureMatch_coupling.xls");
              }
            },
            {
              text: "Export to CSV",
              action(e, dt, button, config) {
                var table_data = sigmatch_table.data().toArray();

                var export_data = [];
                for (let item of Object.values(table_data)) {
                  let record = {};
                  record["name"] = item["entry"];
                  record["family"] = item["family"];
                  record["subfamily"] = item["subfamily"];
                  record["score"] = item["nscore"];
                  record["gtop_gs"] = item["Gs"]["text"];
                  record["gtop_gio"] = item["Gi/o"]["text"];
                  record["gtop_gq11"] = item["Gq/11"]["text"];
                  record["gtop_g1213"] = item["G12/13"]["text"];
                  record["GPCRdb_mean_gs"] = item["Gs_emax"]["text"];
                  record["GPCRdb_mean_gio"] = item["Gi/o_emax"]["text"];
                  record["GPCRdb_mean_gq11"] = item["Gq/11_emax"]["text"];
                  record["GPCRdb_mean_g1213"] = item["G12/13_emax"]["text"];
                  export_data.push(record);
                }

                export_data = Papa.unparse(export_data);

                $("<a></a>")
                  .attr("id", "downloadFile")
                  .attr("href", "data:text/csv;charset=utf8," + encodeURIComponent(export_data))
                  .attr("download", "export.csv")
                  .appendTo("body");

                $("#downloadFile").ready(function() {
                  $("#downloadFile").get(0).click();
                  $("#downloadFile").remove();
                });
              }
            },
            {
              text: "Reset All Filters",
              action() {
                yadcf.exResetAllFilters(sigmatch_table);
              }
            },
            {
              text: "Show Alignment",
              className: "score-button",
              action() {
                sigmatch_table.rows().deselect();
              }
            },
          ]
        });
      } else {
        sigmatch_table = $("#sigmatch_table").DataTable({
          dom: "Bfrtip",
          data: sigmatch_data,
          // scrollY: "60vh",
          scrollY: true,
          scrollX: true,
          scrollCollapse: true,
          scroller: true,
          destroy: true,
          columnDefs: columns_definition,
          order: [
            [5, "desc"]
          ],
          select: {
            style: "single",
          },
          paging: true,
          buttons: [
            {
              text: "Export to Excel",
              action() {
                GlobalTableToExcel("sigmatch_table", "Signature Match data", "SignatureMatch_coupling.xls");
              }
            },
            {
              text: "Export to CSV",
              action(e, dt, button, config) {
                var table_data = sigmatch_table.data().toArray();

                var export_data = [];
                for (let item of Object.values(table_data)) {
                  let record = {};
                  record["name"] = item["entry"];
                  record["family"] = item["family"];
                  record["subfamily"] = item["subfamily"];
                  record["score"] = item["nscore"];
                  record["GPCRdb_mean_arrb1"] = item["arrb1"]["text"];
                  record["GPCRdb_mean_arrb2"] = item["arrb2"]["text"];
                  export_data.push(record);
                }

                export_data = Papa.unparse(export_data);

                $("<a></a>")
                  .attr("id", "downloadFile")
                  .attr("href", "data:text/csv;charset=utf8," + encodeURIComponent(export_data))
                  .attr("download", "export.csv")
                  .appendTo("body");

                $("#downloadFile").ready(function() {
                  $("#downloadFile").get(0).click();
                  $("#downloadFile").remove();
                });
              }
            },
            {
              text: "Reset All Filters",
              action() {
                yadcf.exResetAllFilters(sigmatch_table);
              }
            },
            {
              text: "Show Alignment",
              className: "score-button",
              action() {
                sigmatch_table.rows().deselect();
              }
            },
          ]
        });
      }
      yadcf.init(sigmatch_table.draw(), column_filters, {
        cumulative_filtering: false
      });

      $(".score-button").click(function() {
        //const render_url = window.location.origin + "/signprot/matrix/render_sigmat/";
        window.open("/signprot/matrix/render_sigmat/", "_blank");
      });

      sigmatch_table.on("select", function(e, dt, type, indexes) {
        const entry = sigmatch_table.rows(indexes).data().toArray()[0];
        //console.log(entry);
        svg = d3.select("svg.svg-content.conseq");
        $("g#sigmatch_frame").remove();
        signprotmat.d3.draw_seq_cons(
          entry,
          svg,
          xScale,
          xAxis,
          true
        );
      });
      // Recalculate table layout incl. column widths
      yadcf.exResetAllFilters(sigmatch_table);
    },
    error(jqXHR, exception) {
      console.log(jqXHR);
      console.log(exception);
    },
    complete() {
      // $("#sigm_spin").addClass("fa-check");
      // $("#sigm_spin").removeClass("fa-spinner");
      // $("#sigm_spin").removeClass("fa-spin");
    }
  })
}

const get_max_interface_count = function() {
  // return parseInt($("#single-crystal-group-pdbs-modal-text").text().match(/\d+/))
  return pos_set.length;
}

const filter_pairs = function(floor, ceiling) {
  const max_val = get_max_interface_count();
  d3.select("g#interact")
    .selectAll("g")
    .style("display", function(d) {
      var ratio = d.pairs.length / max_val * 100;
      return (floor <= ratio && ceiling >= ratio ? "block" : "none");
    })
}

const set_min = function() {
  $("#amount_min").html($("#slider-range").slider("values", 0)).position({
    my: "center bottom",
    at: "center top",
    of: $("#slider-range span:eq(0)"),
  });
}

const set_max = function() {
  $("#amount_max").html($("#slider-range").slider("values", 1)).position({
    my: "center bottom",
    at: "center top",
    of: $("#slider-range span:eq(1)"),
  });
}

const initialize_filter_slider = function() {
  // initializing range slider
  $("#slider-range").slider({
    range: true,
    min: 0,
    max: 100,
    step: 1,
    values: [0, 100],
    slide(event, ui) {
      const floor = parseInt(ui.values[0]);
      const ceil = parseInt(ui.values[1]);

      var label = ui.handleIndex === 0 ? "#amount_min" : "#amount_max";
      $(label).html(ui.value).position({
        my: "center bottom",
        at: "center top",
        of: ui.handle,
      });

      filter_pairs(floor, ceil);
    }
  });

  setTimeout(set_min, setTimeout(set_max, 1), 1);
};

const reset_slider = function() {
  // reset the slider to the possible min and max values
  const min = $("#slider-range").slider("option", "min");
  const max = $("#slider-range").slider("option", "max");
  $("#slider-range").slider("values", [min, max]);

  setTimeout(set_min, setTimeout(set_max, 1), 1);
}


$(document).ready(function() {
  // Calling the PDBTableData
  $.get("/signprot/pdbtabledata", {
    exclude_non_interacting: true,
    effector: filtering_particle,
  }, function(data) {
    $("#interface-modal-table .tableview").html(data);
  })

  signprotmat.d3.print_resScaleColor_legend();
  $("[data-toggle=\"tooltip\"]").tooltip();

  let table = $($.fn.dataTable.tables()[0]).DataTable();

  initialize_filter_slider();

  $("#interface-modal-table").on("shown.bs.modal", function(e) {
    showPDBtable("#interface-modal-table");
  })

  $("#interface-modal-table").on("hidden.bs.modal", function(e) {
    table = $($.fn.dataTable.tables()[0]).DataTable();
    let selection = table.rows(".selected").data();

    let old_pdb_sel = pdb_sel;
    pdb_sel = [];
    pos_set = [];

    // get selected pdb ids
    $(".pdb_selected:checked").each(function(index) {
      pdb_sel.push(($(this).attr("id")));
    });

    // api request
    $.ajax({
      type: "POST",
      url: "/signprot/matrix/AJAX_Interactions/",
      async: false,
      data: {
        selected_pdbs: pdb_sel,
        effector: filtering_particle,
        csrfmiddlewaretoken: csrf_token,
      },
      success(data) {
        non_interactions = data[0];
        selected_interactions = data[1];
      }
    });
    // get corresponding protein entry_name values
    for (var int_meta of interactions_metadata) {
      if (pdb_sel.indexOf(int_meta["pdb_id"]) !== -1) {
        pos_set = [int_meta["conf_id"], ...pos_set];
      }
    }

    if (!_.isEqual(old_pdb_sel.sort(), pdb_sel.sort())) {
      $(".svg-content").remove();
      con_seq = {};
      document.querySelector("#seqsig-container").style.display = "none";
      document.querySelector("#conseq-container").style.display = "none";
      document.getElementById("interface-svg").className = "collapse in";
      // data = signprotmat.data.dataTransformationWrapper(interactions, pdb_sel);
      data = signprotmat.data.dataTransformationWrapper(selected_interactions, pdb_sel);
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

      var startTime = performance.now();
      run_seq_sig(data);
      var endTime = performance.now();
      var elapsed = endTime - startTime;
      console.log("run_seq_sig execution Time: " + elapsed);
    };
  });
});

// * CONSTANTS
const margin = { top: 40, right: 200, bottom: 180, left: 200 };
let w = 1200 - margin.left - margin.right,
  h = 1000- margin.top - margin.bottom;

// array for data in infobox
let info_data = [];

const signprotmat = {
  // * DATA TRANSFORMATION FUNCTIONS
  data: {
    extractPdbIDs: function(dataset: object) {
      let ret = [];
      Object.keys(dataset).forEach(e => ret.push(e.toUpperCase()));

      return ret;
    },

    objectToArray: function(dataset: object) {
      return Object.keys(dataset).map(key => dataset[key]);
    },

    moveKeyToArray: function(dataset: Array, pdb_ids: Array) {
      for (let i = 0; i < pdb_ids.length; i++) {
        const pdb = pdb_ids[i];
        for (let j = 0; j < dataset[i].length; j++) {
          dataset[i][j].push(pdb);
        }
      }
      return dataset;
    },

    // https://stackoverflow.com/questions/10865025/merge-flatten-an-array-of-arrays-in-javascript/25804569#comment50580016_10865042
    flattenOnce: array => {
      return [].concat(...array);
    },

    labelData: function(data, keys) {
      const data_labeled = data.map(function(e) {
        let obj = {};
        keys.forEach(function(key, i) {
          // comment this out later
          if (key === "sig_gn") {
            obj[key] = Math.floor(e[i] / 10);
            return;
          }
          obj[key] = e[i];
        });
        return obj;
      });

      return data_labeled;
    },

    getInteractionTypes: function(dataset: Array) {
      let int_ty = [];
      for (let i = 0; i < dataset.length; i++) {
        const int_arr = dataset[i].int_ty;
        int_ty.push(int_arr);
      }
      int_ty = signprotmat.data
        .flattenOnce(int_ty)
        .filter((v, i, a) => a.indexOf(v) === i);
      let rm_index = int_ty.indexOf("undefined");
      if (rm_index > -1) {
        int_ty.splice(rm_index, 1);
      }

      return int_ty;
    },

    get_additional_receptors: function(data, xvals, prids) {
      let new_receptor_data = [];

      for (let index = 0; index < data.length; index++) {
        const e1 = data[index]["rec_gn"];
        const e2 = data[index]["rec_id"];
        if (xvals.includes(e1) && prids.includes(e2)) {
          new_receptor_data.push(data[index]);
        }
      }
      return new_receptor_data;
    },

    extractRecSigData: function(data, which_component: string) {
      if (which_component === "rec") {
        return _.uniqBy(data, t => [t.rec_gn, t.pdb_id].join());
      } else if (which_component === "sig") {
        return _.uniqBy(data, t => [t.sig_gn, t.pdb_id].join());
      } else {
        console.log("No component specified...");
      }
    },

    select_by_value: function (selection, value) {
      let ret_sel = []
      for (let index = 0; index < selection.length; index++) {
          ret_sel.push(selection[index][value]);
      };
      return ret_sel;
    },

    dataTransformationWrapper: function(dataset, keys, pdb_sel) {
      dataset = _.pick(dataset, pdb_sel);
      let pdb_ids = signprotmat.data.extractPdbIDs(dataset);
      let data_t = signprotmat.data.objectToArray(dataset);
      data_t = signprotmat.data.moveKeyToArray(data_t, pdb_ids);
      data_t = signprotmat.data.flattenOnce(data_t);
      data_t = signprotmat.data.labelData(data_t, keys);

      let data_t_rec = signprotmat.data.extractRecSigData(data_t, "rec");
      let data_t_sig = signprotmat.data.extractRecSigData(data_t, "sig");
      let int_ty = signprotmat.data.getInteractionTypes(data_t);

      let return_data = {
        transformed: data_t,
        receptor: data_t_rec,
        signprot: data_t_sig,
        inttypes: int_ty,
        pdbids: pdb_ids
      };

      return return_data;
    }
  },

  // * D3 DRAW FUNCTIONS
  d3: {
    // * SETTING UP SVG FOR OUTPUT
    setup: function(div, loc) {
      if (loc === 'seqsig'){
        h = 1500
      } else if (loc === 'conseq') {
        h = 300
      };
      let svg = d3
        .select("body")
        .select("div#content")
        .select(div)
        .append("svg")
        .attr("preserveAspectRatio", "xMinYMin meet")
        .attr(
          "viewBox",
          "0 0 " +
            (w + margin.left + margin.right) +
            " " +
            (h + 200 + margin.top + margin.bottom)
        )
        // .classed("svg-content", true) //class to make it responsive
        .attr("class", (typeof loc !== 'undefined' ? "svg-content " + loc :"svg-content"))
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

      return svg;
    },

    // * SETTING THE X/Y SCALE
    xScale: function(data) {
      let xScale = d3
        .scaleBand()
        .domain(
          d3
            .map(data, (d: any) => d.rec_gn)
            .keys()
            .sort(d3.ascending)
        )
        .range([0, w])
        // .round(true)
        .padding(1);

      return xScale;
    },

    yScale: function(data) {
      let yScale = d3
        .scaleBand()
        .domain(
          d3
            .map(data, (d: any) => d.sig_gn)
            .keys()
            .sort(d3.descending)
        )
        .range([h, 0])
        // .round(true)
        .padding(1);

      return yScale;
    },

    // * SETTING THE PDB/SIG-PROT SCALE
    pdbScale: function(data) {
      let pdbScale = d3
        .scaleBand()
        .domain(
          d3
            .map(data, (d: any) => d.pdb_id)
            .keys()
            .sort(d3.descending)
        )
        .range([300, 0])
        .padding(1);

      return pdbScale;
    },

    sigScale: function(data) {
      let sigScale = d3
        .scaleBand()
        .domain(
          d3
            .map(data, (d: any) => d.pdb_id)
            .keys()
            .sort(d3.descending)
        )
        .range([120, 0])
        .padding(1);

      return sigScale;
    },

    // * SETTING THE COLOR SCALE
    colScale: function(data) {
      let colScale = d3
        .scaleOrdinal()
        .domain(data)
        .range(d3.schemeSet3);

        return colScale;
      },

      // * seqsig
      // * SETTING THE FEATURE SCALE
      fScale: function(data) {

        const feats = [
          "Hydrophobic",
          "Aliphatic",
          "Aliphatic small",
          "Aliphatic medium",
          "Aliphatic large",
          "Aliphatic extra large",
          "Aromatic",
          "Aromatic medium",
          "Aromatic large",
          "Charge",
          "Charge small",
          "Charge negative",
          "Charge positive",
          "Charge positive small",
          "Charge positive large",
          "Helix propencity extra large",
          "Helix propencity large",
          "Helix propencity medium",
          "Helix propencity small",
          "Helix propencity extra small",
          "Polar or charged",
          "Pro (helix kink)",
          "Polar uncharged extra small",
          "Polar uncharged small",
          "Polar uncharged large",
          "Polar or charged extralarge",
          "Polar or negative small",
          "Polar or negative medium",
          "Polar or negative large",
          "Polar or negative extra large",
          "Polar or positive extra large",
          "Any tiny",
          "Any extra small",
          "Any small",
          "Any medium",
          "Any large",
          "Any extra large",
          "Gly",
          "Ala",
          "Cys",
          "Ser",
          "Thr",
          "Val",
          "Asp",
          "Ile",
          "Leu",
          "Asn",
          "Met",
          "Glu",
          "Gln",
          "His",
          "Lys",
          "Phe",
          "Arg",
          "Tyr",
          "Trp",
          "Gap"]

      let fScale = d3
        .scaleBand()
        .domain(feats)
        .range([0, h])
        // .round(true)
        .padding(1);

      return fScale;
    },

    cScale: function(data){
      // const values = d3
      //   .map(data, (d: any) => d.cons)
      //   .keys()
      //   .map(Number)
      //   .sort(d3.ascending)
      // const min = values.pop()
      // const max = values[0]

      // conservation is calculated to be between -1 and 10 by python
      let cScale = d3
      .scaleSequential(d3.interpolateRdBu)
      .domain([0,10]);

    return cScale;
    },

    // * DEFINING AXIS FOR X/Y AND GRID
    xAxis: function(xScale) {
      let xAxis = d3
        .axisBottom(xScale)
        .tickSize(0)
        .tickPadding(8);

      return xAxis;
    },

    yAxis: function(yScale) {
      let yAxis = d3
        .axisRight(yScale)
        .tickSize(0)
        .tickPadding(8);

      return yAxis;
    },

    xAxisGrid: function(xScale, yScale) {
      let xAxisGrid = d3
        .axisTop(xScale)
        .tickSize(h - yScale.step())
        .tickFormat(d => "");

      return xAxisGrid;
    },

    yAxisGrid: function(xScale, yScale) {
      let yAxisGrid = d3
        .axisRight(yScale)
        .tickSize(w - xScale.step())
        .tickFormat(d => "");

      return yAxisGrid;
    },

    // * ADD TOOLTIP FUNCTIONALITY
    tooltip: function(svg) {
      let tip = d3
        .tip()
        .attr("class", "d3-tip")
        .html(function(d) {
          let pair_string = '';
          d.pairs.forEach(element => {
            pair_string += element.pdb_id + '<br>';
          });
          return "Receptor: " + d.rec_gn +
          "<br>" + "Signaling Protein: " + d.sig_gn +
          "<br>" + 'PDBs:' + '<br>' + pair_string;
        });
      svg.call(tip);

      return tip;
    },

    // * RENDER DATA
    renderData: function(
      svg,
      data,
      xScale,
      yScale,
      xAxis,
      yAxis,
      xAxisGrid,
      yAxisGrid,
      colScale,
      pdbScale,
      sigScale,
      tip
    ) {
      let shift_left: number = 7 / 8;
      let shift_top: number = 1 / 8;
      let scale_size: number = shift_left - shift_top;
      let offset: number = 1;
      let each_res;

      let helper = {};
      let result = data.transformed.reduce(function (r, o) {
        let key = o.rec_gn + '-' + o.sig_gn;

        if (!helper[key]) {
          let tmp = {
            "rec_gn": o.rec_gn,
            "sig_gn": o.sig_gn,
            "pairs": [{
              'pdb_id': o.pdb_id,
              'rec_aa': o.rec_aa,
              'sig_aa': o.sig_aa,
              'int_ty': o.int_ty,
            }]
          };
          helper[key] = tmp;
          r.push(helper[key]);
        } else {
          helper[key].pairs.push({
            'pdb_id': o.pdb_id,
            'rec_aa': o.rec_aa,
            'sig_aa': o.sig_aa,
            'int_ty': o.int_ty,
          })
        }

        return r;
      }, []);

      let bwScale = d3.scaleSequential(d3.interpolateGreys).domain([0,pdbScale.domain().length/2])

      svg
        .append("g")
        .attr("id", "interact")
        .selectAll("rects")
        .data(result)
        .enter()
        .append("rect")
        .attr("x", function(d: any) {
          return xScale(d.rec_gn) - shift_left * xScale.step() + offset;
        })
        .attr("y", function(d: any) {
          return yScale(d.sig_gn) + shift_top * yScale.step() + offset;
        })
        .attr("rx", function() {
          if (data.transformed.length < 15) {
            return 5;
          } else {
            return 3;
          }
        })
        .attr("ry", function() {
          if (data.transformed.length < 15) {
            return 5;
          } else {
            return 3;
          }
        })
        .attr("width", xScale.step() * scale_size)
        .attr("height", yScale.step() * scale_size)
        .attr("fill", (d: any) => bwScale(d.pairs.length))
        .attr("class", (d: any) => "p" + d.pairs.length)
        .on("mouseover", function(d) {
          tip.show(d);
        })
        .on("mouseout", function(d) {
          tip.hide();
        })
        .on("click", function(d) {
          let index;
          // let rect_x = d3.event.target.getAttribute('x')
          // let rect_y = d3.event.target.getAttribute('y')
          // console.log(rect_x, rect_y)

          // https://stackoverflow.com/a/20251369/8160230
          // select the rect under cursor
          let curr = d3.select(this);

          // Determine if current rect was clicked before
          let active = d.active ? false : true;

          // Update whether or not the elements are active
          d.active = active;

          // set style in regards to active
          if (d.active) {
            curr.style("stroke", "black").style("stroke-width", 2);
            info_data.push(d);
          } else {
            curr.style("stroke", "none").style("stroke-width", 2);
            index = info_data.indexOf(d);
            info_data.splice(index, 1);
          }
          signprotmat.d3.infoBoxUpdate();
          signprotmat.d3.colorRecResidues(d);
          signprotmat.d3.colorSigResidues(d);
        });

      // * DRAWING AXES
      svg
        .append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(" + -xScale.step() / 2 + "," + h + ")")
        .call(xAxis)
        .selectAll("text")
        .attr("text-anchor", "end")
        .attr("font-size", "12px")
        .attr("dx", "-5px")
        .attr("dy", "-5px")
        .attr("transform", "rotate(-90)");

      svg
        .append("g")
        .attr("class", "y axis")
        .attr(
          "transform",
          "translate(" + (w - xScale.step()) + "," + yScale.step() / 2 + ")"
        )
        .call(yAxis)
        .selectAll("text")
        .attr("font-size", "12px");

      // * DRAWING GRIDLINES
      svg
        .append("g")
        .attr("class", "x grid")
        .attr("transform", "translate(" + 0 + "," + h + ")")
        .call(xAxisGrid);

      svg
        .append("g")
        .attr("class", "y grid")
        .attr("transform", "translate(" + 0 + "," + yScale.step() + ")")
        .call(yAxisGrid);

      // * ADDITIONAL FIGURE LINES
      // top x line
      svg
        .append("line")
        .style("stroke", "black")
        .attr("x1", 0)
        .attr("y1", yScale.step())
        .attr("x2", w - xScale.step())
        .attr("y2", yScale.step());

      // left y line
      svg
        .append("line")
        .style("stroke", "black")
        .attr("x1", 0)
        .attr("y1", yScale.step())
        .attr("x2", 0)
        .attr("y2", h);

      // * ADD AXIS LABELS
      svg
        .append("text")
        .attr("class", "x axis_label")
        .attr("text-anchor", "end")
        .attr("x", 0)
        .attr("y", h + 15)
        .text("GPCR");

      svg
        .append("text")
        .attr("class", "y axis_label")
        .attr("text-anchor", "begin")
        .attr("x", w - 0.8 * xScale.step())
        .attr("y", 0.8 * yScale.step())
        .text("G-Protein");

      // * ADD INFOBOX ELEMENT
      svg
        .append("g")
        .attr("id", "infobox")
        .attr(
          "transform",
          "translate(-15," + (data.inttypes.length + 2) * 20 + ")"
        );

      // * ADDING COLOR LEGEND
      svg
        .append("g")
        .attr("class", "legendOrdinal")
        .attr("transform", "translate(-30," + yScale.step() + ")");

      let legendOrdinal = d3
        .legendColor()
        .cells(data.inttypes.length)
        .scale(colScale)
        // .cellFilter(function (d) { return d.label !== "undefined" })
        .orient("vertical")
        .labelOffset(-20);

      svg
        .select(".legendOrdinal")
        .call(legendOrdinal)
        .selectAll("rect")
        .attr("rx", 3)
        .attr("ry", 3);
      svg
        .select(".legendOrdinal")
        .selectAll("text")
        .attr("class", "legend")
        .attr("text-anchor", "end");

      // * APPENDING COL TICK ANNOTATION FOR RECEPTOR GNs
      svg
        .append("g")
        .attr("id", "recPDB")
        .attr("transform", "translate(" + 0 + "," + h + ")")
        .selectAll("text")
        .data(data.pdbids)
        .enter()
        .append("text")
        .attr("class", "y seq_label")
        .attr("x", -10)
        .attr("y", function(d: any) {
          return pdbScale(d) - pdbScale.step() / 2;
        })
        .attr("text-anchor", "end")
        .attr("dy", 75)
        .text(function(d: any) {
          return d;
        });

      // * APPENDING ROW TICK ANNOTATION FOR SIGPROT GNs
      svg
        .append("g")
        .attr("id", "sigPDB")
        .attr(
          "transform",
          "translate(" + w + "," + yScale.step() + ")rotate(-90)"
        )
        .selectAll("text")
        .data(data.pdbids)
        .enter()
        .append("text")
        .attr("class", "x seq_label")
        .attr("x", function(d: any, i) {
          return 10;
        })
        .attr("y", function(d: any, i) {
          return sigScale.step() * (i + 1);
        })
        .attr("text-anchor", "begin")
        .attr("dy", 65)
        .text(function(d: any) {
          return d;
        });

      // * APPENDING AMINOACID SEQUENCE [RECEPTOR]
      svg
        .append("g")
        .attr("id", "recAA")
        .attr("transform", "translate(" + -xScale.step() / 2 + "," + h + ")")
        .append("rect")
        .attr("class", "border-bg")
        .style("fill", "#ffffff")
        .attr("x", xScale.step() / 2)
        .attr("y", 75)
        .attr("width", xScale.range()[1] - xScale.step())
        .attr("height", pdbScale.range()[0] - pdbScale.step());

      each_res = svg
        .select("g#recAA")
        .selectAll("text")
        .data(data.receptor)
        .enter()
        .append("g")
        .attr("class", (d: any) => 'R_' + _.replace(d.rec_gn, '.', 'p') + '_P_' + d.pdb_id );

      each_res
        .append("rect")
        .attr("class", "res_rect")
        .style("fill", (d: any) => colScale(d.int_ty[0]))
        .attr("x", (d: any) => xScale(d.rec_gn) - xScale.step() / 2)
        .attr("y", (d: any) => 75 + pdbScale(d.pdb_id) - pdbScale.step())
        .attr("width", xScale.step())
        .attr("height", pdbScale.step());

      each_res
        .append("text")
        .attr("class", "res_label")
        .attr("x", (d: any) => xScale(d.rec_gn))
        .attr("y", (d: any) => pdbScale(d.pdb_id) - pdbScale.step() / 2)
        .attr("text-anchor", "middle")
        .attr("dy", 75)
        .text((d: any) => d.rec_aa);

      d3.select("g#recAA")
        .append("rect")
        .attr("class", "border")
        .style("stroke", "black")
        .style("fill", "none")
        .attr("x", xScale.step() / 2)
        .attr("y", 75)
        .attr("width", xScale.range()[1] - xScale.step())
        .attr("height", pdbScale.range()[0] - pdbScale.step());

      // * APPENDING AMINOACID SEQUENCE [SIGPROT]
      svg
        .append("g")
        .attr("id", "sigAA")
        .attr(
          "transform",
          "translate(" +
            (w + (1 / 3) * margin.right) +
            "," +
            yScale.step() / 2 +
            ")"
        )
        .append("rect")
        .style("fill", "#ffffff")
        .attr("x", 0 + sigScale.step() / 2)
        .attr("y", yScale.step() / 2)
        .attr("width", sigScale.range()[0] - sigScale.step())
        .attr("height", yScale.range()[0] - yScale.step());

      each_res = svg
        .select("g#sigAA")
        .selectAll("text")
        .data(data.signprot)
        .enter()
        .append("g")
        .attr("class", (d: any) => 'S_' + _.replace(d.sig_gn, '.', 'p') + '_P_' + d.pdb_id );

      each_res
        .append("rect")
        .style("fill", (d: any) => colScale(d.int_ty[0]))
        .attr("x", (d: any) => sigScale(d.pdb_id) - sigScale.step() / 2)
        .attr("y", (d: any) => yScale(d.sig_gn) - yScale.step() / 2)
        .attr("width", sigScale.step())
        .attr("height", yScale.step());

      each_res
        .append("text")
        .attr("class", "res_label")
        .attr("x", (d: any) => sigScale(d.pdb_id))
        .attr("y", (d: any) => yScale(d.sig_gn))
        .attr("text-anchor", "middle")
        .text((d: any) => d.sig_aa);

      d3.select("g#sigAA")
        .append("rect")
        .attr("class", "border")
        .style("stroke", "black")
        .style("fill", "none")
        .attr("x", 0 + sigScale.step() / 2)
        .attr("y", yScale.step() / 2)
        .attr("width", sigScale.range()[0] - sigScale.step())
        .attr("height", yScale.range()[0] - yScale.step());

      return svg;
    },

    addReceptor: function(new_data, data, svg) {
      data = _.union(data.transformed, new_data);
      data = signprotmat.data.extractRecSigData(data, "rec");
      let pdb_ids = [];
      _.forEach(_.uniqBy(data, "pdb_id"), value => {
        pdb_ids.push(value["pdb_id"]);
      });

      let pdbScale = signprotmat.d3.pdbScale(data);
      let xScale = signprotmat.d3.xScale(data);

      let selection = svg
        .select("g#recAA")
        .selectAll("text.res_label")
        .data(data);

      let selection_rect = svg
        .select("g#recAA")
        .selectAll("rect.res_rect")
        .data(data);

      let selection_enter = selection.enter().append("g");

      selection_enter
        .append("rect")
        .attr("class", "res_rect")
        .style("fill", "slategrey")
        .attr("x", (d: any) => xScale(d.rec_gn) - xScale.step() / 2)
        .attr("y", (d: any) => 75 + pdbScale(d.pdb_id) - pdbScale.step())
        .attr("width", xScale.step())
        .attr("height", pdbScale.step())
        .merge(selection_rect)
        .transition()
        .duration(500)
        .attr("x", (d: any) => xScale(d.rec_gn) - xScale.step() / 2)
        .attr("y", (d: any) => 75 + pdbScale(d.pdb_id) - pdbScale.step())
        .attr("width", xScale.step())
        .attr("height", pdbScale.step());

      selection_enter
        .append("text")
        .attr("class", "res_label")
        .style("fill", "white")
        .attr("x", (d: any) => xScale(d.rec_gn))
        .attr("y", (d: any) => pdbScale(d.pdb_id) - pdbScale.step() / 2)
        .attr("text-anchor", "middle")
        .attr("dy", 75)
        .text((d: any) => d.rec_aa)
        .merge(selection)
        .transition()
        .duration(500)
        .attr("x", (d: any) => xScale(d.rec_gn))
        .attr("y", (d: any) => pdbScale(d.pdb_id) - pdbScale.step() / 2);

      selection
        .exit()
        .transition()
        .duration(500)
        .remove();

      selection_rect
        .exit()
        .transition()
        .duration(500)
        .remove();

      selection = svg
        .select("g#recPDB")
        .selectAll("text")
        .data(pdb_ids);

      selection_enter = selection.enter();

      selection_enter
        .append("text")
        .attr("class", "y seq_label")
        .attr("x", -10)
        .attr("y", function(d: any) {
          return pdbScale(d) - pdbScale.step() / 2;
        })
        .attr("text-anchor", "end")
        .attr("dy", 75)
        .text(function(d: any) {
          return d;
        })
        .merge(selection)
        .transition()
        .duration(500)
        .attr("x", -10)
        .attr("y", function(d: any) {
          return pdbScale(d) - pdbScale.step() / 2;
        })
        .attr("dy", 75);

      selection
        .exit()
        .transition()
        .duration(500)
        .remove();

      d3.select("g#recAA")
        .selectAll("rect.border, rect.border-bg")
        .transition()
        .duration(500)
        .attr("x", xScale.step() / 2)
        .attr("y", 75)
        .attr("width", xScale.range()[1] - xScale.step())
        .attr("height", pdbScale.range()[0] - pdbScale.step());
    },

    colorRecResidues: function(d) {
      const rec_gn = _.replace(d.rec_gn, '.', 'p')
      const pdb_list = d.pairs.map((x) => x['pdb_id']);

      // select the rect in the g that corresponds to this rec_gn and pdb_id
      pdb_list.forEach(pdb => {
        d3.select('g.' + 'R_' + rec_gn + '_P_' + pdb)
        .select('rect')
        .classed('activeRes', d.active ? true : false);
      });
    },

    colorSigResidues: function(d) {
      const sig_gn = _.replace(d.sig_gn, '.', 'p')
      const pdb_list = d.pairs.map((x) => x['pdb_id']);

      // select the rect in the g that corresponds to this rec_gn and pdb_id
      pdb_list.forEach(pdb => {
        d3.select('g.' + 'S_' + sig_gn + '_P_' + pdb)
        .select('rect')
        .classed('activeRes', d.active ? true : false);
      });
    },

    infoBoxUpdate: function() {
      // create selection and bind data
      let info_box = d3
        .select("g#infobox")
        .selectAll("text")
        .data(info_data);

      // update existing nodes
      info_box
        .attr("y", function(d, i) {
          return i * 15;
        })
        .attr("text-anchor", "end")
        .attr("class", "legend");

      // create nodes for new data
      info_box
        .enter()
        .append("text")
        .attr("y", function(d, i) {
          return i * 15;
        })
        .attr("text-anchor", "end")
        .attr("class", "legend")
        .text(function(d) {
          return d.rec_gn + " : " + d.sig_gn;
        });

      // discard removed nodes
      info_box.exit().remove();

      // print the data again in case it changed
      info_box.text(function(d) {
        return d.rec_gn + " : " + d.sig_gn;
      });
    },

    draw_seq_sig: function(data_in, svg, xScale){
      console.log('running draw seq sig')
      let data = data_in.feat;
      let fScale = signprotmat.d3.fScale(data);
      let cScale = signprotmat.d3.cScale(data);
      let uniq_feats = _.uniq(_.map(data, 'feature'));



      svg
        .append("g")
        .attr("id", "seqsig_feature")
        .attr("transform", "translate(" + 0 + "," + -30 + ")")
        .selectAll("text")
        .data(uniq_feats)
        .enter()
        .append("text")
        .attr("class", "y seq_label")
        .attr("x", -10)
        .attr("y", function(d: any) {
          return fScale(d) - fScale.step() / 2;
        })
        .attr("text-anchor", "end")
        .attr("dy", 75)
        .text(function(d: any) {
          return d;
        });

      svg
        .append("g")
        .attr("id", "seqsig_mat")
        .attr("transform", "translate(" + -xScale.step() / 2 + "," + -30 + ")")
        .append("rect")
        .attr("class", "border-bg")
        .style("fill", "#ffffff")
        .attr("x", xScale.step() / 2)
        .attr("y", 75)
        .attr("width", xScale.range()[1] - xScale.step())
        .attr("height", fScale.range()[1] - fScale.step());

      let each_res = svg
        .select("g#seqsig_mat")
        .selectAll("text")
        .data(data)
        .enter()
        .append("g");

      // the rectangles, colored by conservation
      each_res
        .append("rect")
        .attr("class", "res_rect")
        .style("fill",function(d: any) {
          if(d.cons === -1){
            return '#ffffff';
          } else {return cScale(d.cons);}
        })
        .attr("x", (d: any) => xScale(d.gn) - xScale.step() / 2)
        .attr("y", (d: any) => 75 + fScale(d.feature) - fScale.step())
        .attr("width", xScale.step())
        .attr("height", fScale.step());

      // adding the frequency text to each rectangle
      each_res
        .append("text")
        .attr("class", "res_label")
        .attr("x", (d: any) => xScale(d.gn))
        .attr("y", (d: any) => fScale(d.feature) - fScale.step() / 2)
        .style("fill", (d: any) => {
          if(Math.abs(d.freq) >= 50) {
            return '#eaeaea';
          } else if (Math.abs(d.freq) < 50) {
            return '#000000';
          }
        })
        .attr("text-anchor", "middle")
        .attr("dy", 75)
        .text((d: any) => d.freq);
        // .text((d: any) => _.round(d.freq/100, 1));

      // adding the explanation tooltip to each rectangle



      // putting a black border around the signature
      d3.select("g#seqsig_mat")
        .append("rect")
        .attr("class", "border")
        .style("stroke", "black")
        .style("fill", "none")
        .attr("x", xScale.step() / 2)
        .attr("y", 75)
        .attr("width", xScale.range()[1] - xScale.step())
        .attr("height", fScale.range()[1] - fScale.step());

    },


    draw_seq_cons: function(data_in, svg, xScale){
      console.log('running draw seq cons')
      let data = data_in.cons;
      let fScale = signprotmat.d3.fScale(data);
      let cScale = signprotmat.d3.cScale(data);
      let uniq_feats = _.uniq(_.map(data, 'feature'));


      svg
        .append("g")
        .attr("id", "conseq_mat")
        .attr("transform", "translate(" + -xScale.step() / 2 + "," + -30 + ")")
        .append("rect")
        .attr("class", "border-bg")
        .style("fill", "#ffffff")
        .attr("x", xScale.step() / 2)
        .attr("y", 75)
        .attr("width", xScale.range()[1] - xScale.step())
        .attr("height", 75);

      let each_res = svg
        .select("g#conseq_mat")
        .selectAll("text")
        .data(data)
        .enter()
        .append("g");

      // the rectangles, colored by conservation
      each_res
        .append("rect")
        .attr("class", "res_rect")
        .style("fill",function(d: any) {
          if(d.cons === -1){
            return '#ffffff';
          } else {return cScale(d.cons);}
        })
        .attr("x", (d: any) => xScale(d.gn) - xScale.step() / 2)
        .attr("y", (d: any) => 75)
        .attr("width", xScale.step())
        .attr("height", 75);

      // adding the frequency text to each rectangle
      each_res
        .append("text")
        .attr("class", "res_label")
        // .attr("x", (d: any) => xScale(d.gn))
        // .attr("y", (d: any) => 50)
        .attr("transform", (d: any) => "translate(" + xScale(d.gn) + ",112.5)rotate(270)")
        .style("fill", (d: any) => {
          if(Math.abs(d.score) >= 50) {
            return '#eaeaea';
          } else if (Math.abs(d.score) < 50) {
            return '#000000';
          }
        })
        .attr("text-anchor", "middle")
        .text((d: any) => d.code);
        // .text((d: any) => _.round(d.freq/100, 1));

      // adding the explanation tooltip to each rectangle



      // putting a black border around the signature
      d3.select("g#conseq_mat")
        .append("rect")
        .attr("class", "border")
        .style("stroke", "black")
        .style("fill", "none")
        .attr("x", xScale.step() / 2)
        .attr("y", 75)
        .attr("width", xScale.range()[1] - xScale.step())
        .attr("height", 75);
    }
  }
};

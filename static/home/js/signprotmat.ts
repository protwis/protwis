// * CONSTANTS
const margin = { top: 60, right: 200, bottom: 180, left: 200 };
let w = 1200 - margin.left - margin.right,
  h = 1000 - margin.top - margin.bottom;
// change to 600 for more compact view
const non_int_col = "#fff";

// array for data in infobox
let info_data = [];
let con_seq = {};

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
      const data_labeled = data.map(e => {
        let obj = {};
        keys.forEach(function(key, i) {
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
        const e2 = data[index]["entry_name"];
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

    select_by_value: function(selection, value) {
      let ret_sel = [];
      for (let index = 0; index < selection.length; index++) {
        ret_sel.push(selection[index][value]);
      }
      return ret_sel;
    },

    removeUndefinedGN: function(dataset) {
      return _.filter(dataset, function(o) {
        return o.rec_gn !== "-";
      });
    },

    dataTransformationWrapper: function(dataset, keys, pdb_sel) {
      // dataset = _.pick(dataset, pdb_sel);
      // let pdb_ids = signprotmat.data.extractPdbIDs(dataset);
      // let data_t = signprotmat.data.objectToArray(dataset);
      // data_t = signprotmat.data.moveKeyToArray(data_t, pdb_ids);
      // data_t = signprotmat.data.flattenOnce(data_t);
      let data_t = signprotmat.data.labelData(dataset, keys);
      data_t = signprotmat.data.removeUndefinedGN(data_t);
      data_t = _.filter(data_t, d => pdb_sel.includes(d.pdb_id));

      let data_t_rec = signprotmat.data.extractRecSigData(data_t, "rec");
      let data_t_sig = signprotmat.data.extractRecSigData(data_t, "sig");
      let int_ty = signprotmat.data.getInteractionTypes(data_t);
      let pdb_ids = _.uniqBy(data_t, "pdb_id");
      pdb_ids = _.map(pdb_ids, d => d.pdb_id);

      let return_data = {
        transformed: data_t,
        receptor: data_t_rec,
        signprot: data_t_sig,
        inttypes: int_ty,
        pdbids: pdb_ids
      };

      return return_data;
    },

    annotateNonInteractionData: function(meta, data) {
      data.forEach(element => {
        const tmp = _.find(meta, d => d.entry_name === element.entry_name);
        element["pdb_id"] = tmp.pdb_id;
      });
      return data;
    }
  },

  // * D3 DRAW FUNCTIONS
  d3: {
    // * SETTING UP SVG FOR OUTPUT
    setup: function(div, loc) {
      if (loc === "seqsig") {
        h = 1200;
        let mt = 0;
      } else if (loc === "conseq") {
        h = 0;
        let mt = 30;
      } else {
        h = 1000 - margin.top - margin.bottom;
        let mt = 60;
      }
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
        .attr(
          "class",
          typeof loc !== "undefined" ? "svg-content " + loc : "svg-content"
        )
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + mt + ")");

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

    yScale: function(data, gprot) {
      const domain = d3
        .map(data, (d: any) => d.sig_gn)
        .keys()
        // .sort(d3.descending);
        .sort(function(a, b) {
          const a_patt = /\.(\S*)\./g;
          const b_patt = /\.(\S*)\./g;
          const a_match = a_patt.exec(a);
          const b_match = b_patt.exec(b);
          const a_obj = _.find(gprot, d => d.slug === a_match[1]);
          const b_obj = _.find(gprot, d => d.slug === b_match[1]);
          // console.log(a_obj.id);
          // console.log(a_obj.id + (+a.slice(-2)/100));
          return d3.descending(
            a_obj.id + +a.slice(-2) / 100,
            b_obj.id + +b.slice(-2) / 100
          );
          // return d3.descending(a_obj.id, b_obj.id);
        });

      let yScale = d3
        .scaleBand()
        .domain(domain)
        .range([h, 0])
        .padding(1);

      return yScale;
    },

    // * SETTING THE PDB/SIG-PROT SCALE
    pdbScale: function(data, meta) {
      let pdbScale = d3
        .scaleBand()
        .domain(
          d3
            .map(data, (d: any) => d.pdb_id)
            .keys()
            .sort(function(a, b) {
              const a_obj = _.find(meta, d => d.pdb_id === a);
              const b_obj = _.find(meta, d => d.pdb_id === b);
              return d3.descending(a_obj.entry_name, b_obj.entry_name);
            })
        )
        .range([300, 0])
        .padding(1);

      return pdbScale;
    },

    sigScale: function(data, meta) {
      let sigScale = d3
        .scaleBand()
        .domain(
          d3
            .map(data, (d: any) => d.pdb_id)
            .keys()
            .sort(function(a, b) {
              const a_obj = _.find(meta, d => d.pdb_id === a);
              const b_obj = _.find(meta, d => d.pdb_id === b);
              return d3.descending(a_obj.gprot, b_obj.gprot);
            })
        )
        .range([120, 0])
        .padding(1);

      return sigScale;
    },

    // * SETTING THE COLOR SCALE
    colScale: function(f) {
      const scale = {
        "van-der-waals": "#d9d9d9",
        "edge-to-face": "#969696",
        "water-mediated": "#7DB144",
        hydrophobic: "#93d050",
        "polar-sidechain-sidechain": "#EAA91F",
        "polar-sidechain-backbone": "#C38E1A",
        "polar-backbone-sidechain": "#C3A563",
        "h-bond donor-acceptor": "#7030a0",
        "h-bond acceptor-donor": "#B24DFF",
        "cation-pi": "#0070c0",
        "pi-cation": "#005693",
        ionic: "#00B9BF"
      };
      var colScale = d3
        .scaleOrdinal()
        .domain(Object.keys(scale))
        .range(Object.values(scale));

      return colScale;
    },

    // * seqsig
    // * SETTING THE FEATURE SCALE
    fScale: function(data) {
      const f = _.map(data, function(d) {
        const length_text = d.length != "" ? " (" + d.length + ")" : "";
        return d.feature + length_text;
      });

      let fScale = d3
        .scaleBand()
        .domain(f)
        .range([0, h])
        // .round(true)
        .padding(1);

      return fScale;
    },

    resScaleColor: function(f) {
      const scale = {
        A: { bg_color: "#E6E600", font_color: "#000000" },
        C: { bg_color: "#B2B548", font_color: "#000000" },
        D: { bg_color: "#E60A0A", font_color: "#FDFF7B" },
        E: { bg_color: "#E60A0A", font_color: "#FDFF7B" },
        F: { bg_color: "#18FF0B", font_color: "#000000" },
        G: { bg_color: "#FF00F2", font_color: "#000000" },
        H: { bg_color: "#0093DD", font_color: "#000000" },
        I: { bg_color: "#E6E600", font_color: "#000000" },
        K: { bg_color: "#145AFF", font_color: "#FDFF7B" },
        L: { bg_color: "#E6E600", font_color: "#000000" },
        M: { bg_color: "#E6E600", font_color: "#000000" },
        N: { bg_color: "#A70CC6", font_color: "#FDFF7B" },
        P: { bg_color: "#CC0099", font_color: "#FDFF7B" },
        Q: { bg_color: "#A70CC6", font_color: "#FDFF7B" },
        R: { bg_color: "#145AFF", font_color: "#FDFF7B" },
        S: { bg_color: "#A70CC6", font_color: "#FDFF7B" },
        T: { bg_color: "#A70CC6", font_color: "#FDFF7B" },
        V: { bg_color: "#E6E600", font_color: "#000000" },
        W: { bg_color: "#0BCF00", font_color: "#000000" },
        Y: { bg_color: "#18FF0B", font_color: "#000000" },
        "-": { bg_color: "#FFFFFF", font_color: "#000000" },
        _: { bg_color: "#EDEDED", font_color: "#000000" },
        "+": { bg_color: "#FFFFFF", font_color: "#000000" }
      };
      return scale[f];
    },

    fScaleColor: function(f) {
      if (f === "αH") {
        f = "aH";
      }
      const scale = {
        HY: { bg_color: "#93d050" },
        HA: { bg_color: "#ffff00" },
        M: { bg_color: "#ffff00" },
        A: { bg_color: "#ffff00" },
        I: { bg_color: "#ffff00" },
        L: { bg_color: "#ffff00" },
        V: { bg_color: "#ffff00" },
        HR: { bg_color: "#07b050" },
        W: { bg_color: "#07b050" },
        Y: { bg_color: "#07b050" },
        F: { bg_color: "#07b050" },
        Hb: { bg_color: "#7030a0", font_color: "#ffffff" },
        N: { bg_color: "#7030a0", font_color: "#ffffff" },
        Q: { bg_color: "#7030a0", font_color: "#ffffff" },
        S: { bg_color: "#7030a0", font_color: "#ffffff" },
        T: { bg_color: "#7030a0", font_color: "#ffffff" },
        Hu: { bg_color: "#7030a0", font_color: "#ffffff" },
        Ha: { bg_color: "#7030a0", font_color: "#ff0000" },
        Hd: { bg_color: "#7030a0", font_color: "#02b0f0" },
        "+-": { bg_color: "#0070c0", font_color: "#ff0000" },
        "+": { bg_color: "#0070c0", font_color: "#000000" },
        H: { bg_color: "#0070c0", font_color: "#000000" },
        K: { bg_color: "#0070c0", font_color: "#000000" },
        R: { bg_color: "#0070c0", font_color: "#000000" },
        "-": { bg_color: "#ff0000" },
        D: { bg_color: "#ff0000" },
        E: { bg_color: "#ff0000" },
        Sm: { bg_color: "#ffffff" },
        aH: { bg_color: "#d9d9d9" },
        G: { bg_color: "#ff02ff" },
        P: { bg_color: "#d603ff", font_color: "#ffffff" },
        C: { bg_color: "#bf8f00" }
      };
      return scale[f];
    },

    colorBySwitch: function(which, colScale) {
      if (which === "res") {
        const other = "int";
        let res = d3.select("g#recAA").selectAll("g");

        res.selectAll("rect").style("fill", function(d) {
          const col = signprotmat.d3.resScaleColor(d.rec_aa);
          if (typeof col != "undefined") {
            return col.bg_color;
          } else {
            return null;
          }
        });

        res.selectAll("text").style("fill", function(d) {
          const col = signprotmat.d3.resScaleColor(d.rec_aa);
          if (typeof col != "undefined") {
            if (typeof col.font_color != "undefined") {
              return col.font_color;
            } else {
              return "#000000";
            }
          } else {
            return "#000000";
          }
        });

        res = d3.select("g#sigAA").selectAll("g");
        res.selectAll("rect").style("fill", function(d) {
          const col = signprotmat.d3.resScaleColor(d.sig_aa);
          if (typeof col != "undefined") {
            return col.bg_color;
          } else {
            return null;
          }
        });

        res.selectAll("text").style("fill", function(d) {
          const col = signprotmat.d3.resScaleColor(d.sig_aa);
          if (typeof col != "undefined") {
            if (typeof col.font_color != "undefined") {
              return col.font_color;
            } else {
              return "#000000";
            }
          } else {
            return "#000000";
          }
        });
      } else if (which === "int") {
        const other = "res";
        let res = d3.select("g#recAA").selectAll("g");

        res.selectAll("rect").style("fill", function(d) {
          if (typeof d.int_ty != "undefined") {
            return colScale(d.int_ty[0]);
          } else {
            return non_int_col;
          }
        });
        res.selectAll("text").style("fill", null);
        res = d3.select("g#sigAA").selectAll("g");
        res.selectAll("rect").style("fill", function(d) {
          if (typeof d.int_ty != "undefined") {
            return colScale(d.int_ty[0]);
          } else {
            return "#ccc";
          }
        });
        res.selectAll("text").style("fill", null);
      }
      document.querySelector("#" + which + "but").classList.add("active");
      document.querySelector("#" + other + "but").classList.remove("active");
    },

    cScale: function(data) {
      // const values = d3
      //   .map(data, (d: any) => d.cons)
      //   .keys()
      //   .map(Number)
      //   .sort(d3.ascending)
      // const min = values.pop()
      // const max = values[0]

      // conservation is calculated to be between -1 and 10 by python
      let cScale = d3.scaleSequential(d3.interpolateGreys).domain([0, 100]);

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
          let pair_string = "";
          d.pairs.forEach(e => {
            pair_string +=
              e.pdb_id +
              ": " +
              e.rec_aa +
              " vs. " +
              e.sig_aa +
              "(" +
              e.int_ty +
              ")" +
              "<br>";
          });
          return (
            "Receptor: " +
            d.rec_gn +
            "<br>" +
            "Signaling Protein: " +
            d.sig_gn +
            "<br>" +
            "PDBs: " +
            "<br>" +
            pair_string
          );
        });
      svg.call(tip);

      return tip;
    },

    // * RENDER DATA
    renderData: function(
      svg,
      data,
      data_non,
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
      tip
    ) {
      let shift_left: number = 7 / 8;
      let shift_top: number = 1 / 8;
      let scale_size: number = shift_left - shift_top;
      let offset: number = 1;
      let each_res;

      let helper = {};
      let result = data.transformed.reduce(function(r, o) {
        let key = o.rec_gn + "-" + o.sig_gn;

        if (!helper[key]) {
          let tmp = {
            rec_gn: o.rec_gn,
            sig_gn: o.sig_gn,
            pairs: [
              {
                pdb_id: o.pdb_id,
                rec_aa: o.rec_aa,
                sig_aa: o.sig_aa,
                int_ty: o.int_ty
              }
            ]
          };
          helper[key] = tmp;
          r.push(helper[key]);
        } else {
          helper[key].pairs.push({
            pdb_id: o.pdb_id,
            rec_aa: o.rec_aa,
            sig_aa: o.sig_aa,
            int_ty: o.int_ty
          });
        }

        return r;
      }, []);

      let bwScale = d3
        .scaleSequential(d3.interpolateGreys)
        .domain([0, pdbScale.domain().length]);

      svg
        .append("g")
        .attr("id", "interact")
        .selectAll("rects")
        .data(result)
        .enter()
        .append("rect")
        .attr("x", function(d: any) {
          // return xScale(d.rec_gn) - shift_left * xScale.step() + offset;
          return xScale(d.rec_gn) - xScale.step();
        })
        .attr("y", function(d: any) {
          // return yScale(d.sig_gn) + shift_top * yScale.step() + offset;
          return yScale(d.sig_gn);
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
        // .attr("width", xScale.step() * scale_size)
        // .attr("height", yScale.step() * scale_size)
        .attr("width", xScale.step())
        .attr("height", yScale.step())
        .attr("fill", (d: any) => {
          if (pdbScale.domain().length > 1) {
            return bwScale(d.pairs.length);
          } else {
            return colScale(d.pairs[0].int_ty[0]);
          }
        })
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
            curr.style("stroke", "yellow").style("stroke-width", 2);
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

      // * ADDING Interaction Type LEGEND
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
          const i_obj = _.find(interactions_metadata, e => e.pdb_id === d);
          let text = i_obj.name.replace("&beta;", "\u03B2"); // beta
          text = text.replace("&mu;", "\u03BC"); // mu
          return text.replace(/<[^>]*>/g, "") + " (" + d.toUpperCase() + ")";
          // return d;
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
          // return sigScale(d) - sigScale.step() / 2;
          return sigScale(d);
          // return sigScale.step() * (i + 1);
        })
        .attr("text-anchor", "begin")
        .attr("dy", 68)
        .text(function(d: any) {
          const i_obj = _.find(interactions_metadata, e => e.pdb_id === d);
          // let text = i_obj.gprot.replace('Engineered', 'Eng.')
          let text = i_obj.gprot.replace("Engineered", "E.");
          // text = text.replace('protein', 'prot.')
          text = text.replace("protein", "p.");
          return text.replace(/<[^>]*>/g, "") + " (" + d.toUpperCase() + ")";
        });

      // * APPENDING AMINOACID SEQUENCE [RECEPTOR]
      let recTip = d3
        .tip()
        .attr("class", "d3-tip")
        .html(function(d) {
          return (
            "Signal Prot. AA: " +
            d.sig_aa +
            "<br>" +
            "Interaction type: " +
            d.int_ty
          );
        });

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

      data_non = _.filter(data_non, function(d) {
        return xScale(d.rec_gn);
      });

      data_non = _.filter(data_non, function(d) {
        return pdbScale(d.pdb_id);
      });

      // data.receptor.push(...data_non)
      data_non.push(...data.receptor);

      each_res = svg
        .select("g#recAA")
        .selectAll("text")
        .data(data_non)
        .enter()
        .append("g")
        .attr(
          "class",
          (d: any) => "R_" + _.replace(d.rec_gn, ".", "p") + "_P_" + d.pdb_id
        )
        .call(recTip)
        .on("mouseover", function(d) {
          recTip.show(d);
        })
        .on("mouseout", function(d) {
          recTip.hide();
        });

      each_res
        .append("rect")
        .attr("class", "res_rect")
        .style("fill", (d: any) =>
          typeof d.int_ty !== "undefined" ? colScale(d.int_ty[0]) : non_int_col
        )
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
      let sigTip = d3
        .tip()
        .attr("class", "d3-tip")
        .html(function(d) {
          return (
            "Receptor AA: " +
            d.rec_aa +
            "<br>" +
            "Interaction type: " +
            d.int_ty
          );
        });

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
        .attr(
          "class",
          (d: any) => "S_" + d.sig_gn.replace(/\./gi, "p") + "_P_" + d.pdb_id
        )
        .call(sigTip)
        .on("mouseover", function(d) {
          sigTip.show(d);
        })
        .on("mouseout", function(d) {
          sigTip.hide();
        });

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
      const rec_gn = _.replace(d.rec_gn, ".", "p");
      const pdb_list = d.pairs.map(x => x["pdb_id"]);

      // select the rect in the g that corresponds to this rec_gn and pdb_id
      pdb_list.forEach(pdb => {
        d3.select("g." + "R_" + rec_gn + "_P_" + pdb)
          .select("rect")
          .classed("activeRes", d.active ? true : false);
      });
    },

    colorSigResidues: function(d) {
      const sig_gn = d.sig_gn.replace(/\./gi, "p");
      const pdb_list = d.pairs.map(x => x["pdb_id"]);

      // select the rect in the g that corresponds to this rec_gn and pdb_id
      pdb_list.forEach(pdb => {
        d3.select("g." + "S_" + sig_gn + "_P_" + pdb)
          .select("rect")
          .classed("activeRes", d.active ? true : false);
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

    conSeqUpdate: function(row_height) {
      let cScale = signprotmat.d3.cScale();

      let seqsigTip = d3
        .tip()
        .attr("class", "d3-tip")
        .html(function(d) {
          return (
            "Generic Residue No.: " +
            d.gn +
            "<br>" +
            "Feature: " +
            d.feature +
            "<br>" +
            "Length: " +
            d.length +
            "<br>" +
            // "Score: " +
            // d.expl +
            // "<br>" +
            "Frequency: " +
            d.freq +
            "<br>"
          );
        });

      // create selection and bind data
      let con_seq_mat = d3
        .select("g#con_seq_mat")
        .selectAll("g")
        .data(Object.values(con_seq))
        .join("g");

      const svg = d3.select("svg.svg-content.seqsig");
      const t = svg.transition().duration(750);

      // PROPERTY CODES
      con_seq_mat
        .selectAll("rect.res_rect")
        .data(d => d)
        .join(
          enter =>
            enter
              .append("rect")
              .attr("class", "res_rect")
              .call(seqsigTip)
              .on("mouseover", function(d) {
                seqsigTip.show(d);
              })
              .on("mouseout", function(d) {
                seqsigTip.hide();
              })
              .style("fill", function(d: any) {
                const gcol = signprotmat.d3.fScaleColor(d.feature_code);
                if (typeof gcol != "undefined") {
                  return gcol.bg_color;
                } else {
                  return null;
                }
              })
              //.style("fill", '#ffffff')
              //.call(enter => enter.transition(t)
              //.style("fill", function(d: any) {
              //const gcol = signprotmat.d3.fScaleColor(d.feature_code);
              //if (typeof gcol != "undefined") {
              //return gcol.bg_color;
              //} else {
              //return null;
              //}
              //}))
              .style("stroke", "black")
              .attr("x", (d: any) => xScale(d.gn) - xScale.step() / 2)
              .attr("y", 75)
              .attr("width", xScale.step())
              .attr("height", row_height),
          update =>
            update
              .attr("x", (d: any) => xScale(d.gn) - xScale.step() / 2)
              .attr("y", 75)
              //.style("fill", '#ffffff')
              .call(update =>
                update.transition(t).style("fill", function(d: any) {
                  const gcol = signprotmat.d3.fScaleColor(d.feature_code);
                  if (typeof gcol != "undefined") {
                    return gcol.bg_color;
                  } else {
                    return null;
                  }
                })
              ),
          exit => exit.remove()
        );

      con_seq_mat
        .selectAll("text.res_label")
        .data(d => d)
        .join(
          enter =>
            enter
              .append("text")
              .attr("class", "res_label")
              .attr("text-anchor", "middle")
              .attr("x", (d: any) => xScale(d.gn))
              .attr("y", 75)
              .attr("dy", row_height / 2)
              .style("fill", (d: any) => {
                const gcol = signprotmat.d3.fScaleColor(d.feature_code);
                if (typeof gcol != "undefined") {
                  if (typeof gcol.font_color != "undefined") {
                    return gcol.font_color;
                  } else {
                    return "#000000";
                  }
                } else {
                  return "#000000";
                }
              })
              .text(d => d.feature_code),
          update =>
            update
              .attr("x", (d: any) => xScale(d.gn))
              .attr("y", 75)
              .style("fill", (d: any) => {
                const gcol = signprotmat.d3.fScaleColor(d.feature_code);
                if (typeof gcol != "undefined") {
                  if (typeof gcol.font_color != "undefined") {
                    return gcol.font_color;
                  } else {
                    return "#000000";
                  }
                } else {
                  return "#000000";
                }
              })
              .text(d => d.feature_code),
          exit => exit.remove()
        );

      // PROPERTY LENGTH
      con_seq_mat
        .selectAll("rect.len_rect")
        .data(d => d)
        .join(
          enter =>
            enter
              .append("rect")
              .attr("class", "len_rect")
              .style("fill", "#ffffff")
              .style("stroke", "black")
              .attr("x", (d: any) => xScale(d.gn) - xScale.step() / 2)
              .attr("y", (d: any) => 75 + row_height)
              .attr("width", xScale.step())
              .attr("height", row_height / 2),
          update =>
            update
              .attr("x", (d: any) => xScale(d.gn) - xScale.step() / 2)
              .attr("y", (d: any) => 75 + row_height),
          exit => exit.remove()
        );

      con_seq_mat
        .selectAll("text.len_label")
        .data(d => d)
        .join(
          enter =>
            enter
              .append("text")
              .attr("class", "len_label")
              .attr("text-anchor", "middle")
              .attr("x", (d: any) => xScale(d.gn))
              .attr("y", 75 + row_height)
              .attr("dy", row_height / 4)
              .style("fill", "#000000")
              .text(d => d.length),
          update =>
            update
              .attr("x", (d: any) => xScale(d.gn))
              .attr("y", 75 + row_height)
              .style("fill", "#000000")
              .text(d => d.length),
          exit => exit.remove()
        );

      // CONSERVATION
      con_seq_mat
        .selectAll("rect.cons_rect")
        .data(d => d)
        .join(
          enter =>
            enter
              .append("rect")
              .attr("class", "cons_rect")
              .style("fill", function(d: any) {
                if (d.cons === -1) {
                  return "#ffffff";
                } else {
                  return cScale(d.freq);
                }
              })
              .style("stroke", "black")
              .attr("x", (d: any) => xScale(d.gn) - xScale.step() / 2)
              .attr("y", (d: any) => 75 + 1.5 * row_height)
              .attr("width", xScale.step())
              .attr("height", row_height / 2),
          update =>
            update
              .style("fill", function(d: any) {
                if (d.cons === -1) {
                  return "#ffffff";
                } else {
                  return cScale(d.freq);
                }
              })
              .attr("x", (d: any) => xScale(d.gn) - xScale.step() / 2)
              .attr("y", (d: any) => 75 + 1.5 * row_height),
          exit => exit.remove()
        );

      con_seq_mat
        .selectAll("text.cons_label")
        .data(d => d)
        .join(
          enter =>
            enter
              .append("text")
              .attr("class", "cons_label")
              .attr("text-anchor", "middle")
              .attr("x", (d: any) => xScale(d.gn))
              .attr("y", 75 + 1.5 * row_height)
              .attr("dy", row_height / 4)
              .style("fill", (d: any) => {
                if (Math.abs(d.freq) >= 50) {
                  return "#eaeaea";
                } else if (Math.abs(d.freq) < 50) {
                  return "#000000";
                }
              })
              .text(d => d.freq),
          update =>
            update
              .attr("x", (d: any) => xScale(d.gn))
              .attr("y", 75 + 1.5 * row_height)
              .style("fill", (d: any) => {
                if (Math.abs(d.freq) >= 50) {
                  return "#eaeaea";
                } else if (Math.abs(d.freq) < 50) {
                  return "#000000";
                }
              })
              .text(d => d.freq),
          exit => exit.remove()
        );
    },
    draw_seq_sig: function(data_in, svg, xScale) {
      let data = data_in.feat;
      let cScale = signprotmat.d3.cScale();
      let feats = [];

      for (let key of Object.keys(data)) {
        if (xScale(key) == null) {
          data = _.omit(data, key);
        }
      }

      let col_lengths = [];
      for (let elem of Object.keys(data)) {
        col_lengths.push(data[elem].length);
      }
      const row_height = 30;
      const area_height = _.max(col_lengths) * row_height;

      let seqsigTip = d3
        .tip()
        .attr("class", "d3-tip")
        .html(function(d) {
          return (
            "Generic Residue No.: " +
            d.gn +
            "<br>" +
            "Feature: " +
            d.feature +
            "<br>" +
            "Length: " +
            d.length +
            "<br>" +
            // "Score: " +
            // d.expl +
            // "<br>" +
            "Frequency: " +
            d.freq +
            "<br>"
          );
        });
     
      let viewbox_svg = d3.select('.svg-content.seqsig')
     
      const viewbox = viewbox_svg.attr('viewBox')
      const viewbox_1 = viewbox.slice(0, 9)
      
      const new_vb = viewbox_1 + (area_height + Math.round(area_height / 4) + 150)
      viewbox_svg.attr('viewBox', new_vb)

      // generating the white backdrop for all the properties
      svg
        .append("g")
        .attr("id", "seqsig_mat")
        .attr("transform", "translate(" + -xScale.step() / 2 + "," + 90 + ")")
        .append("rect")
        .attr("class", "border-bg")
        .style("fill", "#ffffff")
        .attr("x", xScale.step() / 2)
        .attr("y", 75)
        .attr("width", xScale.range()[1] - xScale.step())
        .attr("height", area_height);

      let each_res = svg
        .select("g#seqsig_mat")
        .selectAll("g")
        .data(Object.values(data))
        .enter()
        .append("g")
        .selectAll("rect")
        .data(function(d) {
          return d;
        })
        .enter();

      /*
       *            each_res
       *                .append('rect')
       *                .call(seqsigTip)
       *                .on("mouseover", function(d) {
       *                    if (d.freq !== 0) {
       *                        seqsigTip.show(d);
       *                    }
       *                })
       *                .on("mouseout", function(d) {
       *                    seqsigTip.hide();
       *                })
       *                .attr("class", "res_rect")
       *                .style("fill", d => cScale(d.freq))
       *                .attr("x", (d: any) => xScale(d.gn) - xScale.step() / 2)
       *                .attr("y", (d: any, i: number) => 75 + (i * row_height))
       *                .attr("width", xScale.step())
       *                .attr("height", row_height);
       *
       *
       *            each_res
       *                .append("text")
       *                .attr("class", "res_label")
       *                .attr("x", (d: any) => xScale(d.gn))
       *                .attr("y", (d: any, i: number) => 75 + (i * row_height))
       *                .style("fill", (d: any) => {
       *                    if(Math.abs(d.freq) >= 50) {
       *                        return '#eaeaea';
       *                    } else if (Math.abs(d.freq) < 50) {
       *                        return '#000000';
       *                    }
       *                })
       *                .attr("text-anchor", "middle")
       *                .attr("dy", row_height / 2)
       *                .text((d: any) => d.freq);
       */

      each_res
        .append("rect")
        .call(seqsigTip)
        .on("mouseover", function(d) {
          if (d.freq !== 0) {
            seqsigTip.show(d);
          }
        })
        .on("mouseout", function(d) {
          seqsigTip.hide();
        })
        // .on("click", function(d) {
        //   let index;
        //   // let rect_x = d3.event.target.getAttribute('x')
        //   // let rect_y = d3.event.target.getAttribute('y')
        //   // console.log(rect_x, rect_y)
        //
        //   // https://stackoverflow.com/a/20251369/8160230
        //   // select the rect under cursor
        //   let curr = d3.select(this);
        //   d.active = true;
        //
        //   // set style in regards to active
        //   if (
        //     d.active &&
        //     con_seq[d.gn].length > 0 &&
        //     d.key == con_seq[d.gn][0].key
        //   ) {
        //     // remove prev. sele. rectangle
        //     con_seq[d.gn].active = false;
        //     //con_seq = _.omit(con_seq, d.gn)
        //     con_seq[d.gn] = [];
        //     curr.style("stroke", "black").style("stroke-width", 1);
        //   } else if (d.active && con_seq[d.gn].length > 0) {
        //     // change the selected property to the clicked one
        //     // remove the active value from the prev. selec. property
        //     con_seq[d.gn].active = false;
        //     // remove the previously selected property
        //     //con_seq = _.omit(con_seq, d.gn)
        //     con_seq[d.gn] = [];
        //     // add the new prop to the object
        //     con_seq[d.gn] = [d];
        //     // apply appropriate syling
        //     curr.style("stroke", "yellow").style("stroke-width", 2);
        //     d3.select("g#seqsig_mat")
        //       .selectAll("rect.res_rect")
        //       .filter(e => e.active)
        //       .filter(e => e.gn == d.gn)
        //       .filter(
        //         e => e.feature_code + e.length != d.feature_code + d.length
        //       )
        //       .style("stroke", "black")
        //       .style("stroke-width", 1);
        //   } else if (d.active) {
        //     // add a newly clicked rectangle
        //     con_seq[d.gn] = [d];
        //     curr.style("stroke", "yellow").style("stroke-width", 2);
        //   } else {
        //     // edge case rect removal
        //     con_seq[d.gn].active = false;
        //     con_seq = _.omit(con_seq, d.gn);
        //     curr.style("stroke", "black").style("stroke-width", 1);
        //   }
        //   //console.log(con_seq)
        //
        //   signprotmat.d3.conSeqUpdate(row_height);
        // })
        .attr("class", "res_rect")
        .style("fill", function(d: any) {
          const gcol = signprotmat.d3.fScaleColor(d.feature_code);
          if (typeof gcol != "undefined") {
            return gcol.bg_color;
          } else {
            return null;
          }
        })
        .style("stroke", "black")
        .attr("x", (d: any) => xScale(d.gn) - xScale.step() / 2)
        .attr("y", (d: any, i: number) => 75 + i * row_height)
        .attr("width", xScale.step())
        .attr("height", row_height);

      each_res
        .append("text")
        .attr("class", "res_label")
        .attr("text-anchor", "middle")
        .attr("x", (d: any) => xScale(d.gn))
        .attr("y", (d: any, i: number) => 75 + i * row_height)
        .attr("dy", row_height / 2)
        .style("fill", (d: any) => {
          const gcol = signprotmat.d3.fScaleColor(d.feature_code);
          if (typeof gcol != "undefined") {
            if (typeof gcol.font_color != "undefined") {
              return gcol.font_color;
            } else {
              return "#000000";
            }
          } else {
            return "#000000";
          }
        })
        .text(function(d: any) {
          return d.feature_code;
        });

      d3.select("svg.svg-content.seqsig")
        .select("g")
        .append("g")
        .attr("id", "con_seq_mat")
        .attr("transform", "translate(" + -xScale.step() / 2 + "," + -10 + ")");
      //.append("rect")
      //.attr("class", "border-bg")
      //.style("fill", "#ffffff")
      //.attr("x", xScale.step() / 2)
      //.attr("y", 75)
      //.attr("width", xScale.range()[1] - xScale.step())
      //.attr("height", row_height * 1.5)

      //d3.select("g#con_seq_mat")

      //d3.select("g#con_seq_mat")
      //.append("rect")
      //.attr("class", "border")
      //.style("stroke", "black")
      //.style("fill", "none")
      //.attr("x", xScale.step() / 2)
      //.attr("y", 75)
      //.attr("width", xScale.range()[1] - xScale.step())
      //.attr("height", row_height * 1.5);

      d3.select("svg.svg-content.seqsig")
        .select("g")
        .append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(" + -xScale.step() / 2 + "," + 50 + ")")
        .call(xAxis)
        .selectAll("text")
        .attr("text-anchor", "start")
        .attr("font-size", "12px")
        .attr("dx", "-5px")
        .attr("dy", "-5px")
        .attr("transform", "rotate(-90)");

      svg
        .select(".x.axis")
        .selectAll("path")
        .remove();

      // putting a black border around the signature
      /*
       *d3.select("g#seqsig_mat")
       *    .append("rect")
       *    .attr("class", "border")
       *    .style("stroke", "black")
       *    .style("fill", "none")
       *    .attr("x", xScale.step() / 2)
       *    .attr("y", 75)
       *    .attr("width", xScale.range()[1] - xScale.step())
       *    .attr("height", area_height);
       */
    },

    draw_seq_cons: function(data_in, svg, xScale, xAxis, sigmatch) {
      let data = data_in.cons;
      let fScale = signprotmat.d3.fScale(data);
      let cScale = signprotmat.d3.cScale(data);
      let uniq_feats = _.uniq(_.map(data, "feature"));

      // filter out NA generic numbers based on xScale
      data = _.filter(data, function(d) {
        return xScale(d.gn);
      });

      let conseqTip = d3
        .tip()
        .attr("class", "d3-tip")
        .html(function(d) {
          return (
            "Generic Residue No.: " +
            d.gn +
            "<br>" +
            "Feature: " +
            d.feature +
            "<br>" +
            "Length: " +
            d.length +
            "<br>" +
            "Score: " +
            d.score +
            "<br>"
          );
        });

      if (sigmatch) {
        svg = svg
          .append("g")
          .attr("id", "sigmatch_frame")
          .attr(
            "transform",
            "translate(" + margin.left + "," + (margin.top + 140) + ")"
          );

        svg
          .append("g")
          .attr("id", "sigmatch_mat")
          .attr(
            "transform",
            "translate(" + -xScale.step() / 2 + "," + -30 + ")"
          )
          .append("rect")
          .attr("class", "border-bg")
          .style("fill", "#ffffff")
          .attr("x", xScale.step() / 2)
          .attr("y", 75)
          .attr("width", xScale.range()[1] - xScale.step())
          .attr("height", 75);

        svg
          .append("text")
          .attr("class", "y seq_label")
          .attr("text-anchor", "end")
          .attr("x", -10)
          .attr("y", 65)
          .text("Property");

        svg
          .append("text")
          .attr("class", "y seq_label")
          .attr("text-anchor", "end")
          .attr("x", -10)
          .attr("y", 102)
          .text("Conservation");

        let group = "g#sigmatch_mat";
      } else {
        svg
          .append("g")
          .attr("id", "conseq_mat")
          .attr(
            "transform",
            "translate(" + -xScale.step() / 2 + "," + -30 + ")"
          )
          .append("rect")
          .attr("class", "border-bg")
          .style("fill", "#ffffff")
          .attr("x", xScale.step() / 2)
          .attr("y", 75)
          .attr("width", xScale.range()[1] - xScale.step())
          .attr("height", 75);

        svg
          .append("g")
          .attr("class", "x axis")
          .attr("transform", "translate(" + -xScale.step() / 2 + "," + 35 + ")")
          .call(xAxis)
          .selectAll("text")
          .attr("text-anchor", "start")
          .attr("font-size", "12px")
          .attr("dx", "-5px")
          .attr("dy", "-5px")
          .attr("transform", "rotate(-90)");

        svg
          .select(".x.axis")
          .selectAll("path")
          .remove();

        svg
          .append("text")
          .attr("class", "y seq_label")
          .attr("text-anchor", "end")
          .attr("x", -10)
          .attr("y", 65)
          .text("Property");

        svg
          .append("text")
          .attr("class", "y seq_label")
          .attr("text-anchor", "end")
          .attr("x", -10)
          .attr("y", 102)
          .text("Conservation");

        let group = "g#conseq_mat";
      }

      let each_res = svg
        .select(group)
        .selectAll("text")
        .data(data)
        .enter()
        .append("g")
        .call(conseqTip)
        .on("mouseover", function(d) {
          conseqTip.show(d);
        })
        .on("mouseout", function(d) {
          conseqTip.hide();
        });

      // the rectangles, colored by feature
      each_res
        .append("rect")
        .attr("class", "res_rect")
        .style("fill", function(d: any) {
          const gcol = signprotmat.d3.fScaleColor(d.code);
          if (typeof gcol != "undefined") {
            return gcol.bg_color;
          } else {
            return null;
          }
        })
        .attr("x", (d: any) => xScale(d.gn) - xScale.step() / 2)
        .attr("y", (d: any) => 75)
        .attr("width", xScale.step())
        .attr("height", 37.5);

      // the rectangles, colored by conservation
      each_res
        .append("rect")
        .attr("class", "res_rect")
        .style("fill", function(d: any) {
          if (d.cons === -1) {
            return "#ffffff";
          } else {
            return cScale(d.score);
          }
        })
        .attr("x", (d: any) => xScale(d.gn) - xScale.step() / 2)
        .attr("y", (d: any) => 75 + 37.5)
        .attr("width", xScale.step())
        .attr("height", 37.5);

      // adding the feature text to each rectangle
      each_res
        .append("text")
        .attr("class", "res_label")
        // .attr("x", (d: any) => xScale(d.gn))
        // .attr("y", (d: any) => 50)
        .attr(
          "transform",
          (d: any) => "translate(" + xScale(d.gn) + ",93.75)" // + "rotate(270)"
        )
        .style("fill", (d: any) => {
          const gcol = signprotmat.d3.fScaleColor(d.code);
          if (typeof gcol != "undefined") {
            if (typeof gcol.font_color != "undefined") {
              return gcol.font_color;
            } else {
              return "#000000";
            }
          } else {
            return "#000000";
          }
        })
        .attr("text-anchor", "middle")
        .text((d: any) => d.code);

      // adding the conservation value to each rectangle
      each_res
        .append("text")
        .attr("class", "res_label")
        // .attr("x", (d: any) => xScale(d.gn))
        // .attr("y", (d: any) => 50)
        .attr(
          "transform",
          (d: any) =>
            "translate(" + xScale(d.gn) + "," + (75 + 37.5 + 37.5 / 2) + ")" // + "rotate(270)"
        )
        .style("fill", (d: any) => {
          if (Math.abs(d.score) >= 50) {
            return "#eaeaea";
          } else if (Math.abs(d.score) < 50) {
            return "#000000";
          }
        })
        .attr("text-anchor", "middle")
        .text((d: any) => d.score);

      // putting a black border around the signature
      d3.select(group)
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

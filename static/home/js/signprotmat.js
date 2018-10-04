// * CONSTANTS
var margin = { top: 40, right: 200, bottom: 180, left: 130 };
var w = 1200 - margin.left - margin.right, h = 900 - margin.top - margin.bottom;
var signprotmat = {
    data: {
        // * DATA TRANSFORMATION FUNCTIONS
        extractPdbIDs: function (dataset) {
            return Object.keys(dataset);
        },
        objectToArray: function (dataset) {
            return Object.keys(dataset).map(function (key) { return dataset[key]; });
        },
        moveKeyToArray: function (dataset, pdb_ids) {
            for (var i = 0; i < pdb_ids.length; i++) {
                var pdb = pdb_ids[i];
                for (var j = 0; j < dataset[i].length; j++) {
                    dataset[i][j].push(pdb);
                }
            }
            ;
            return dataset;
        },
        // https://stackoverflow.com/questions/10865025/merge-flatten-an-array-of-arrays-in-javascript/25804569#comment50580016_10865042
        flattenOnce: function (array) { return [].concat.apply([], array); },
        labelData: function (data, keys) {
            var data_labeled = data.map(function (e) {
                var obj = {};
                keys.forEach(function (key, i) {
                    // comment this out later
                    if (key === 'sig_gn') {
                        obj[key] = Math.floor(e[i] / 10);
                        return;
                    }
                    obj[key] = e[i];
                });
                return obj;
            });
            return data_labeled;
        },
        getInteractionTypes: function (dataset) {
            var int_ty = [];
            for (var i = 0; i < dataset.length; i++) {
                var int_arr = dataset[i].int_ty;
                int_ty.push(int_arr);
            }
            int_ty = signprotmat.data.flattenOnce(int_ty).filter(function (v, i, a) { return a.indexOf(v) === i; });
            var rm_index = int_ty.indexOf("undefined");
            if (rm_index > -1) {
                int_ty.splice(rm_index, 1);
            }
            return int_ty;
        },
        dataTransformationWrapper: function (dataset, keys, pdb_sel) {
            dataset = _.pick(dataset, pdb_sel);
            var pdb_ids = signprotmat.data.extractPdbIDs(dataset);
            var data_t = signprotmat.data.objectToArray(dataset);
            data_t = signprotmat.data.moveKeyToArray(data_t, pdb_ids);
            data_t = signprotmat.data.flattenOnce(data_t);
            data_t = signprotmat.data.labelData(data_t, keys);
            var data_t_rec = _.uniqBy(data_t, function (t) { return [t.rec_gn, t.pdb_id].join(); });
            var data_t_sig = _.uniqBy(data_t, function (t) { return [t.sig_gn, t.pdb_id].join(); });
            var int_ty = signprotmat.data.getInteractionTypes(data_t);
            var return_data = {
                transformed: data_t,
                receptor: data_t_rec,
                signprot: data_t_sig,
                inttypes: int_ty
            };
            return return_data;
        }
    },
    // * D3 DRAW FUNCTIONS
    d3: {
        // * SETTING UP SVG FOR OUTPUT
        setup: function () {
            var svg = d3
                .select("body")
                .select("div#content")
                .append("div")
                .classed("svg-container", true) //container class to make it responsive
                .append("svg")
                .attr("preserveAspectRatio", "xMinYMin meet")
                .attr("viewBox", "0 0 " +
                (w + margin.left + margin.right) +
                " " +
                (h + 200 + margin.top + margin.bottom))
                .classed("svg-content", true) //class to make it responsive
                .append("g")
                .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
            return svg;
        },
        // * SETTING THE X/Y SCALE
        xScale: function (data) {
            var xScale = d3
                .scaleBand()
                .domain(d3
                .map(data.transformed, function (d) { return d.rec_gn; })
                .keys()
                .sort(d3.ascending))
                .range([0, w])
                // .round(true)
                .padding(1);
            return xScale;
        },
        yScale: function (data) {
            var yScale = d3
                .scaleBand()
                .domain(d3
                .map(data.transformed, function (d) { return d.sig_gn; })
                .keys()
                .sort(d3.descending))
                .range([h, 0])
                // .round(true)
                .padding(1);
            return yScale;
        },
        // * SETTING THE PDB/SIG-PROT SCALE
        pdbScale: function (data) {
            var pdbScale = d3
                .scaleBand()
                .domain(d3
                .map(data.transformed, function (d) { return d.pdb_id; })
                .keys()
                .sort(d3.descending))
                .range([180, 0])
                .padding(1);
            return pdbScale;
        },
        sigScale: function (data) {
            var sigScale = d3
                .scaleBand()
                .domain(d3
                .map(data.transformed, function (d) { return d.pdb_id; })
                .keys()
                .sort(d3.descending))
                .range([120, 0])
                .padding(1);
            return sigScale;
        },
        // * SETTING THE COLOR SCALE
        colScale: function (data) {
            var colScale = d3
                .scaleOrdinal()
                .domain(data.inttypes)
                .range(d3.schemeSet3);
            return colScale;
        },
        // * DEFINING AXIS FOR X/Y AND GRID
        xAxis: function (xScale) {
            var xAxis = d3
                .axisBottom(xScale)
                .tickSize(0)
                .tickPadding(8);
            return xAxis;
        },
        yAxis: function (yScale) {
            var yAxis = d3
                .axisRight(yScale)
                .tickSize(0)
                .tickPadding(8);
            return yAxis;
        },
        xAxisGrid: function (xScale, yScale) {
            var xAxisGrid = d3
                .axisTop(xScale)
                .tickSize(h - yScale.step())
                .tickFormat(function (d) { return ""; });
            return xAxisGrid;
        },
        yAxisGrid: function (xScale, yScale) {
            var yAxisGrid = d3
                .axisRight(yScale)
                .tickSize(w - xScale.step())
                .tickFormat(function (d) { return ""; });
            return yAxisGrid;
        }
    }
};
// * ADD TOOLTIP FUNCTIONALITY
var tip = d3
    .tip()
    .attr("class", "d3-tip")
    .html(function (d) {
    return d.rec_gn + "<br>" + d.sig_gn + "<br>" + d.int_ty;
});
svg.call(tip);
// * RENDER DATA
var shift_left = 7 / 8;
var shift_top = 1 / 8;
var scale_size = shift_left - shift_top;
var offset = 1;
// array for data in infobox
var info_data = [];
svg
    .append("g")
    .attr("id", "interact")
    .selectAll("rects")
    .data(data.transformed)
    .enter()
    .append("rect")
    .attr("x", function (d) {
    return xScale(d.rec_gn) - shift_left * xScale.step() + offset;
})
    .attr("y", function (d) {
    return yScale(d.sig_gn) + shift_top * yScale.step() + offset;
})
    .attr("rx", function () {
    if (data.transformed.length < 15) {
        return 5;
    }
    else {
        return 3;
    }
})
    .attr("ry", function () {
    if (data.transformed.length < 15) {
        return 5;
    }
    else {
        return 3;
    }
})
    .attr("width", xScale.step() * scale_size)
    .attr("height", yScale.step() * scale_size)
    .attr("fill", function (d) {
    if (d.int_ty === undefined) {
        return "none";
    }
    else {
        return colScale(d.int_ty[0]);
    }
})
    .on("mouseover", function (d) {
    tip.show(d);
})
    .on("mouseout", function (d) {
    tip.hide();
})
    .on("click", function (d) {
    var index;
    // let rect_x = d3.event.target.getAttribute('x')
    // let rect_y = d3.event.target.getAttribute('y')
    // console.log(rect_x, rect_y)
    // https://stackoverflow.com/a/20251369/8160230
    // select the rect under cursor
    var curr = d3.select(this);
    // Determine if current rect was clicked before
    var active = d.active ? false : true;
    // Update whether or not the elements are active
    d.active = active;
    // set style in regards to active
    if (d.active) {
        curr.style("stroke", "black").style("stroke-width", 2);
        info_data.push(d);
    }
    else {
        curr.style("stroke", "none").style("stroke-width", 2);
        index = info_data.indexOf(d);
        info_data.splice(index, 1);
    }
    infoBoxUpdate();
});
// * ADD INFOBOX ELEMENT
svg
    .append("g")
    .attr("id", "infobox")
    .attr("transform", "translate(-15," + (data.inttypes.length + 2) * 20 + ")");
function infoBoxUpdate() {
    // create selection and bind data
    var info_box = d3
        .select("g#infobox")
        .selectAll("text")
        .data(info_data);
    // update existing nodes
    info_box
        .attr("y", function (d, i) {
        return i * 15;
    })
        .attr("text-anchor", "end")
        .attr("class", "legend");
    // create nodes for new data
    info_box
        .enter()
        .append("text")
        .attr("y", function (d, i) {
        return i * 15;
    })
        .attr("text-anchor", "end")
        .attr("class", "legend")
        .text(function (d) {
        return d.rec_gn + " : " + d.sig_gn;
    });
    // discard removed nodes
    info_box.exit().remove();
    // print the data again in case it changed
    info_box.text(function (d) {
        return d.rec_gn + " : " + d.sig_gn;
    });
}
// * ADDING COLOR LEGEND
svg
    .append("g")
    .attr("class", "legendOrdinal")
    .attr("transform", "translate(-30," + yScale.step() + ")");
var legendOrdinal = d3
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
    .data(Object.keys(dataset))
    .enter()
    .append("text")
    .attr("class", "x axis_label")
    .attr("x", -10)
    .attr("y", function (d) {
    return pdbScale(d);
})
    .attr("text-anchor", "end")
    .attr("dy", 75)
    .text(function (d) {
    return d;
});
// * APPENDING ROW TICK ANNOTATION FOR SIGPROT GNs
svg
    .append("g")
    .attr("id", "sigPDB")
    .attr("transform", "translate(" + w + "," + yScale.step() + ")rotate(-45)")
    .selectAll("text")
    .data(Object.keys(dataset))
    .enter()
    .append("text")
    .attr("class", "x axis_label")
    .attr("x", function (d) {
    return sigScale(d);
})
    .attr("y", function (d) {
    return sigScale(d);
})
    .attr("text-anchor", "begin")
    .attr("dx", 45)
    .attr("dy", 45)
    .text(function (d) {
    return d;
});
// * APPENDING AMINOACID SEQUENCE [RECEPTOR]
svg
    .append("g")
    .attr("id", "recAA")
    .attr("transform", "translate(" + -xScale.step() / 2 + "," + h + ")")
    .append("rect")
    .style("fill", "#eaeaea")
    .attr("x", xScale.step() / 2)
    .attr("y", 75)
    .attr("width", xScale.range()[1] - xScale.step())
    .attr("height", pdbScale.range()[0] - pdbScale.step());
var each_res = svg
    .select('g#recAA')
    .selectAll("text")
    .data(data.receptor)
    .enter()
    .append('g');
each_res
    .append('rect')
    .style("fill", function (d) { return colScale(d.int_ty); })
    .attr("x", function (d) { return xScale(d.rec_gn) - xScale.step() / 2; })
    .attr("y", function (d) { return 75 + pdbScale(d.pdb_id) - pdbScale.step(); })
    .attr('width', xScale.step())
    .attr('height', pdbScale.step());
each_res
    .append("text")
    .attr("class", "res_label")
    .attr("x", function (d) { return xScale(d.rec_gn); })
    .attr("y", function (d) { return pdbScale(d.pdb_id); })
    .attr("text-anchor", "middle")
    .attr("dy", 75)
    .text(function (d) { return d.rec_aa; });
d3.select('g#recAA').append('rect')
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
    .attr("transform", "translate(" + (w + (1 / 3) * margin.right) + "," + yScale.step() / 2 + ")")
    .append("rect")
    .style("fill", "#eaeaea")
    .attr("x", 0 + sigScale.step() / 2)
    .attr("y", yScale.step() / 2)
    .attr("width", sigScale.range()[0] - sigScale.step())
    .attr("height", yScale.range()[0] - yScale.step());
var each_res = svg
    .select('g#sigAA')
    .selectAll("text")
    .data(data.signprot)
    .enter()
    .append('g');
each_res
    .append('rect')
    .style("fill", function (d) { return colScale(d.int_ty); })
    .attr("x", function (d) { return sigScale(d.pdb_id) - sigScale.step() / 2; })
    .attr("y", function (d) { return yScale(d.sig_gn) - yScale.step() / 2; })
    .attr('width', sigScale.step())
    .attr('height', yScale.step());
each_res
    .append("text")
    .attr("class", "res_label")
    .attr("x", function (d) { return sigScale(d.pdb_id); })
    .attr("y", function (d) { return yScale(d.sig_gn); })
    .attr("text-anchor", "middle")
    .attr("dy", 5)
    .text(function (d) { return d.sig_aa; });
d3.select('g#sigAA').append("rect")
    .style("stroke", "black")
    .style("fill", "none")
    .attr("x", 0 + sigScale.step() / 2)
    .attr("y", yScale.step() / 2)
    .attr("width", sigScale.range()[0] - sigScale.step())
    .attr("height", yScale.range()[0] - yScale.step());
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
    .attr("transform", "translate(" + (w - xScale.step()) + "," + yScale.step() / 2 + ")")
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
//# sourceMappingURL=signprotmat.js.map
// Like d3.svg.diagonal.radial, but with square corners.
function step(startAngle, startRadius, endAngle, endRadius) {
    var c0 = Math.cos(startAngle = (startAngle - 90) / 180 * Math.PI),
        s0 = Math.sin(startAngle),
        c1 = Math.cos(endAngle = (endAngle - 90) / 180 * Math.PI),
        s1 = Math.sin(endAngle);
    return "M" + startRadius * c0 + "," + startRadius * s0
        + (endAngle === startAngle ? "" : "A" + startRadius + "," + startRadius + " 0 0 " + (endAngle > startAngle ? 1 : 0) + " " + startRadius * c1 + "," + startRadius * s1)
        + "L" + endRadius * c1 + "," + endRadius * s1;
}

function string_pixlen(text, depth) {
    var canvas = document.createElement('canvas');
    var ctx = canvas.getContext("2d");
    if (depth < 2) {
        ctx.font = "20px Palatino"
    } else if (depth == 2) {
        ctx.font = "14px Palatino"
    } else {
        ctx.font = "12px Palatino"
    }
    return parseInt(ctx.measureText(text).width) + 40;
}

function draw_tree(data, options) {

    var branches = {};
    var branch_offset = 0;
    for (var key in options.branch_length) {
        branch_offset = branch_offset + string_pixlen(options.branch_length[key], key);
        branches[key] = branch_offset;
    }
    branches[options.depth] = branch_offset + 10;

    //var color = d3.scale.category20();

    var diameter = 2 * branches[options.depth] + 150;

    var tree = d3.layout.tree()
        .size([360, diameter / 2])
        .separation(function (a, b) { return (a.parent == b.parent ? 1 : 2) / a.depth; });

    var diagonal = d3.svg.diagonal.radial()
        .projection(function (d) { return [d.y, d.x / 180 * Math.PI]; });

    var svg = d3.select('#'+options.anchor).append("svg")
        .attr("width", diameter)
        .attr("height", diameter)
        .attr("id", options.anchor+"_svg")
        .attr("xmlns", "http://www.w3.org/2000/svg");

    svg.append("rect")
        .attr("width", "100%")
        .attr("height", "100%")
        .attr("fill", "white");

    var svg_g = svg.append("g")
        .attr("transform", "translate(" + diameter / 2 + "," + diameter / 2 + ")");

    var nodes = tree.nodes(data);

    nodes.forEach(function (d) {
        if (d.depth == 0) {
            d.y = 0
        } else {
            d.y = branches[d.depth]
        }
    });

    var links = tree.links(nodes);

    var link = svg_g.append("g")
        .attr("class", "links")
        .selectAll("path")
        .data(links)
        .enter().append("path")
        .each(function (d) { d.target.linkNode = this; })
        .attr("d", diagonal) //function (d) { return step(d.source.x, d.source.y, d.target.x, d.target.y) })
        .style("stroke", function (d) { return d.target.color; })
        .style("stroke-width", function (d) { if (d.target.depth > 0) { return 4 - d.target.depth; } else { return 0; } })
        .style("opacity", function (d) {
            if ((d.target.interactions > 0 && d.target.mutations_an > 0) || 1 == 1) { return 0.8 } //|| 1==1
            else if (d.target.interactions > 0) { return 0.5 }
            else if (d.target.mutations_an > 0) { return 0.5 }
            else { return 0.1 };
        });

    var node = svg_g.selectAll(".node")
        .data(nodes)
        .enter().append("g")
        .attr("class", "node")
        .attr("transform", function (d) { if (d.name == '') { return "rotate(" + (d.x) + ")translate(" + d.y + ")"; } else { return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")"; } })

    node.filter(function (d) { return (d.depth == options.depth) }).append("circle")
        .attr("r", function (d) { if (d.name == '') { return "0" } else { return "4.0" } })
        .style("fill", function (d) {
            if (d.color && d.depth < options.depth) { return d.color }
            else if (d.value > 0) {
                return "FireBrick";
            }
            //Here should go code for ligands and mutations
            else { return "#eee" };
        })
        .style("opacity", .99);

    node.append("text")
        .attr("dy", ".31em")
        .attr("text-anchor", function (d) { return d.x < 180 ? "end" : "start"; })
        .attr("transform", function (d) {
            if (d.depth == 3) {
                return d.x < 180 ? "translate(50)" : "rotate(180)translate(-50)";
            } else {
                return d.x < 180 ? "translate(-12)" : "rotate(180)translate(12)";
            }
        })
        .text(function (d) {
            if (d.depth == options.depth) {
                return d.name.toUpperCase();
            } else if (options.label_free.includes(d.depth)) {
                return "";
            } else if (d.depth > 0) {
                return d.name;
            } else {
                return "";
            }
        })

        .style("font-size", function (d) { if (d.depth < 2) { return "14px" } else if (d.depth == 2) { return "12px" } else { return "10px" } })
        .style("font-family", "Palatino")
        .style("fill", function (d) {
            if (d.color) { return "#111" }
            //else if (d.interactions > 0 && d.mutations_an > 0 && 1 == 2) { return "green" }
            //else if (d.interactions > 0 && 1 == 2) { return "Olive" }
            //else if (d.mutations_an > 0 && 1 == 2) { return "palegreen" }
            else { return "#222" };
        }).call(getBB);
    node.filter(function (d) { return (d.depth != options.depth) }).insert("rect", "text")
        .attr("x", function (d) { return d.x < 180 ? d.bbox.x - 12 : d.bbox.x - d.bbox.width - 12; })
        .attr("y", function (d) { return d.bbox.y })
        .attr("width", function (d) { return d.bbox.width })
        .attr("height", function (d) { return d.bbox.height })
        .style("fill", "#FFF");

    function getBB(selection) {
        selection.each(function (d) { d.bbox = this.getBBox(); })
    }
}
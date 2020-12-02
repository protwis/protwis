function make_coverage_tree() {

    // Like d3.svg.diagonal.radial, but with square corners.
    function step(startAngle, startRadius, endAngle, endRadius) {
        var c0 = Math.cos(startAngle = (startAngle - 90) / 180 * Math.PI),
            s0 = Math.sin(startAngle),
            c1 = Math.cos(endAngle = (endAngle - 90) / 180 * Math.PI),
            s1 = Math.sin(endAngle);
        return "M" + startRadius * c0 + "," + startRadius * s0 +
            (endAngle === startAngle ? "" : "A" + startRadius + "," + startRadius + " 0 0 " + (endAngle > startAngle ? 1 : 0) + " " + startRadius * c1 + "," + startRadius * s1) +
            "L" + endRadius * c1 + "," + endRadius * s1;
    }

    var color = d3.scale.category20();

    var diameter = 1350;

    function getBB(selection) {
        selection.each(function(d) { d.bbox = this.getBBox(); })
    }

    function getAngle(d) {
        // Offset the angle by 90 deg since the '0' degree axis for arc is Y axis, while
        // for text it is the X axis.
        var thetaDeg = (180 / Math.PI * (arc.startAngle()(d) + arc.endAngle()(d)) / 2 - 90);
        // If we are rotating the text by more than 90 deg, then "flip" it.
        // This is why "text-anchor", "middle" is important, otherwise, this "flip" would
        // a little harder.
        return (thetaDeg > 90) ? thetaDeg - 180 : thetaDeg;
    }

    // Stash the old values for transition.
    function stash(d) {
        d.x0 = d.x;
        d.dx0 = d.dx;
    }

    //////////////////////////////

    var color = d3.scale.category20c();

    var only_data = 0;

    var width = 1350,
        height = 1350,
        radius = Math.min(width, height) / 2 - 50;

    var x = d3.scale.linear()
        .range([0, 2 * Math.PI]);

    var y = d3.scale.linear()
        .range([0, radius]);

    var svg3 = d3.select("#svg")
        .append("svg")
        // .attr("width", width)
        // .attr("height", height)
        .attr("xmlns", "http://www.w3.org/2000/svg")
        .attr("version", "1.1")
        .attr("width", "100%")
        .attr("viewBox", "0 0 " + width + " " + height)
        .append("g")
        .attr("id", "scan")
        .attr("transform", "translate(" + width / 2 + "," + (height / 2) + ")");

    var partition = d3.layout.partition()
        .sort(function(a, b) { return b.sort - a.sort; })
        .value(function(d) {
            if (only_data) {
                if (d.depth != 4) { return d.receptor_i + d.receptor_m_an; } else { return Math.min(1, d.receptor_i + d.receptor_m_an); }
            } else {
                return d.receptor_t;
            }
        });
    // .value(function(d) { return d.receptor_t });

    var arc = d3.svg.arc()
        .startAngle(function(d) { return Math.max(0, Math.min(2 * Math.PI, x(d.x))); })
        .endAngle(function(d) { return Math.max(0, Math.min(2 * Math.PI, x(d.x + d.dx))); })
        .innerRadius(function(d) {

            // console.log("inner" + d.y+" "+Math.max(0, y(d.y))); return Math.max(0, y(d.y));
            if (d.y < 0.3) {
                return 2;
            } else if (d.y < 0.5) {
                return 180;
            } else if (d.y < 0.7) {
                return 150 + (radius - 45 - 150) * 0.4;
            } else if (d.y < 0.9) {
                return radius - 45;
            }

        })
        .outerRadius(function(d) {

            // console.log("outer" + d.y+" "+Math.max(0, y(d.y + d.dy))); return Math.max(0, y(d.y + d.dy));
            if (d.y < 0.3) {
                return 180;
            } else if (d.y < 0.5) {
                return 150 + (radius - 45 - 150) * 0.4;
            } else if (d.y < 0.7) {
                return radius - 45;
            } else if (d.y < 0.9) {
                return radius + 10; // determines box size
            }
        });

    var lineFunction = d3.svg.line()
        .x(function(d) { return d.x; })
        .y(function(d) { return d.y; })
        .interpolate("linear");

    var count;

    function getDescendants(node) {
        if (!node.children) {
            return 0;
        }
        var total = 0;
        node.children.forEach(function(d) {
            total += 1 + getDescendants(d);
        })
        return total;
    }

    var tip = d3.tip()
        .attr('class', 'd3-tip')
        .style("visibility", "visible")
        .offset([0, 0])
        .html(function(d) {
            if (d.depth === 4 && d.number_of_variants > 0) { return d.name.toUpperCase() + "<br>" + "<span style='color:#fff'> Number of isoforms: " + d.number_of_variants} else if (d.depth < 4 && d.number_of_variants > 0) {
                return d.name.toUpperCase() + "<br>" + "<span style='color:#fff'> Number of isoforms: " + d.number_of_variants;
            } else { return d.name + "<br>" + "<span style='color:#fff'> NO DATA AVAILABLE"; }
        });
    tip(svg3.append("g"));

    var path = svg3.datum(interactions).selectAll("path")
        .data(partition.nodes)
        .enter().append("g")

    path.append("path")
        .attr("display", function(d) { return d.depth ? null : "none"; }) // hide inner ring
        .attr("d", arc)
        .style("fill", function(d) {
            if (d.depth < 4) { return "#c40100" } //&& d.depth<4 && 1==2
            else { return "#c40100" }
        })
        .style("stroke-width", function(d) { return (d.depth < 4 ? "2px" : "1px") })
        .style("stroke", function(d) {
            if (d.depth < 4) { return "#fff" } //&& d.depth<4 && 1==2
            else { return "#fff" }
        })

        .style("opacity", function(d) {
            if (d.depth === 4) { return d.density_of_variants } // || 1==1
            else { return (d.density_of_variants) };
        })
        // .style("fill-rule", "evenodd")
        .on('mouseover', function(d) {
            tip.show(d)
            if (d.depth === 4) { d3.select(this).style("cursor", "pointer") }
        })
        .on('mouseout', function(d) {
            tip.hide(d)
            d3.select(this).style("cursor", "cursor")
        })
        .on("click", function(d) {
            if (d.depth === 4) { filter_this_receptor(d.name); };
        })
        .each(stash);

    path.filter(function(d) { return (d.depth == 1 || d.name == 'Other GPCRs') }).each(function(d, i) {
        svg3.append("line")
            .attr("x1", function() { return 2 * Math.cos(x(d.x) - Math.PI / 2); })
            .attr("y1", function() { return 2 * Math.sin(x(d.x) - Math.PI / 2); })
            .attr("x2", function() { return (radius - 36) * Math.cos(x(d.x) - Math.PI / 2); })
            .attr("y2", function() { return (radius - 36) * Math.sin(x(d.x) - Math.PI / 2); })
            .attr("stroke-width", 2)
            .attr("stroke", "#fff");
    });


    function a(t) {
        return .299 * t.r + .587 * t.g + .114 * t.b
    }

    function n(t) {
        if (t.children) {
            var e = t.children.map(n),
                r = d3.hsl(e[0]),
                a = d3.hsl(e[1]);
            return d3.hsl((r.h + a.h) / 2, 1.2 * r.s, r.l / 1.2)
        }
        return t.colour || "#fff"
    }

    var i = radius * 2,
        l = i,
        o = i / 2,
        d = d3.scale.linear().range([0, 2 * Math.PI]),
        u = d3.scale.pow().exponent(1.3).domain([0, 1]).range([0, o]),
        c = 5,
        s = 1e3

    path.each(function(d, i) {

        svg3.append("text")
            .text(function() {
                return (((d.depth == 3 && d.value == 1 && only_data != 1) || d.value == 0 || d.depth == 0) ? "" : (d.depth == 4 ? d.name.toUpperCase() : d.name));
            })
            .classed("label2", true)
            .attr("x", function() { if (d.depth == 4) { return d.x; } else { return d.x; } })
            // .attr("text-anchor", "end")
            .attr("text-anchor", function() {
                // console.log(d.name+ " "+arc.centroid(d));
                if (d.depth == 4) {
                    return arc.centroid(d)[0] < 0 ? "start" : "end";
                } else if (d.depth > 0) {
                    return arc.centroid(d)[0] < 0 ? "start" : "end";
                } else {
                    return "middle";
                }
            })
            // translate to the desired point and set the rotation
            .attr("transform", function() {
                if (d.depth > 0 && d.depth < 4) {
                    return "translate(" + arc.centroid(d) + ")" +
                        "rotate(" + getAngle(d) + ")";
                } else if (d.depth == 4) {
                    // return d.x < 180 ? "translate(8)" : "rotate(180)translate(-8)";
                    return "translate(" + arc.centroid(d) + ")" +
                        "rotate(" + getAngle(d) + ")";
                } else {
                    return null;
                }
            })
            // .attr("dx", "1") // margin
            .attr("dx", function() {
                if (d.depth == 4) {
                    return arc.centroid(d)[0] < 0 ? "-25px" : "25px";
                } else if (d.depth == 1) {
                    return arc.centroid(d)[0] < 0 ? "-65px" : "65px";
                } else if (d.depth == 2) {
                    return arc.centroid(d)[0] < 0 ? "-65px" : "65px";
                } else if (d.depth == 3) {
                    return arc.centroid(d)[0] < 0 ? "-100px" : "100px";
                } else { return "0px"; }
            }) // vertical-align
            .attr("dy", ".35em") // vertical-align

            // .style("font", function () { return (d.depth==1 ? "14px" : d.depth!=4 ? "12px" : "10px") + " Palatino" ;})
            //
            // .style("font", "Palatino")
            .style("font", function() {
                if (d.depth < 2) { return "21px" + " Palatino"; } else if (d.depth == 2) { return "17px" + " Palatino"; } else if (d.depth == 3) { return "16px" + " Palatino"; } else { return "14px" + " Palatino"; }
            })

            .style("cursor", function() {
                if (d.depth == 4) { return "pointer"; } else { return "cursor"; }
            })

            .on('mouseover', function() {
                tip.show(d)
                if (d.depth === 4) { d3.select(this).style("cursor", "pointer") }
            })
            .on('mouseout', function() {
                tip.hide(d)
                d3.select(this).style("cursor", "cursor")
            })

            .on("click", function() {
                if (d.depth === 4) { filter_this_receptor(d.name); };
            })

            .style("fill", function() { return "#000" })
    })
}
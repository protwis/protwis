// highlight link and connected nodes on mouseover
var hivesvg;
function linkMouseoverHP(d) {
  hivesvg.selectAll(".link")
    .classed("turnedOn", function(dl) {
      return dl === d;
    })
    .classed("turnedOff", function(dl) {
      return !(dl === d);
    })
  hivesvg.selectAll(".node")
    .classed("turnedOn", function(dl) {
      return dl === d.source || dl === d.target;
    })
}

// highlight node and connected links on mouseover
function nodeMouseoverHP(d) {
  hivesvg.selectAll(".link")
    .classed("turnedOn", function(dl) {
      return dl.source === d || dl.target === d;
    })
    .classed("turnedOff", function(dl) {
      return !(dl.source === d || dl.target === d)
    });
  d3.select(this)
    .classed("turnedOn", true);
}

// clear highlighted nodes or links
function mouseoutHP() {
  hivesvg.selectAll(".turnedOn").classed("turnedOn", false);
  hivesvg.selectAll(".turnedOff").classed("turnedOff", false);
}

function degreesHP(radians) {
  return radians / Math.PI * 180 - 90;
}

function createHiveplot(data, container) {
    var width = 800,
        height = 800,
        innerRadius = 0.15*width,
        outerRadius = 0.45*width;

    var angle = d3.scale.ordinal()
                  .domain(d3.range(8))
                  .rangePoints([0, 2 * Math.PI]),
        radius = d3.scale.linear()
                  .range([innerRadius, outerRadius]),
        color = d3.scale.category10()
                  .domain(d3.range(20));

    var nodes = [
      {x: 0, y: .1},
      {x: 0, y: .9},
      {x: 1, y: .2},
      {x: 1, y: .3},
      {x: 2, y: .1},
      {x: 2, y: .8},
      {x: 3, y: .1},
      {x: 3, y: .8},
      {x: 4, y: .1},
      {x: 4, y: .8},
      {x: 5, y: .1},
      {x: 5, y: .8},
      {x: 6, y: .1},
      {x: 6, y: .8}
    ];

    var links = [
      {source: nodes[0], target: nodes[2]},
      {source: nodes[1], target: nodes[3]},
      {source: nodes[2], target: nodes[4]},
      {source: nodes[2], target: nodes[5]},
      {source: nodes[5], target: nodes[10]},
      {source: nodes[5], target: nodes[8]}
    ];

    hivesvg = d3.select(container).append("svg")
        .attr("width", width)
        .attr("height", height)
      .append("g")
        .attr("transform", "translate(" + width/2 + "," + height/2 + ")");

    hivesvg.selectAll(".axis")
        .data(d3.range(7))
      .enter().append("line")
        .attr("class", "axis")
        .attr("transform", function(d) { return "rotate(" + degreesHP(angle(d)) + ")" })
        .attr("x1", radius.range()[0])
        .attr("x2", radius.range()[1]);

    // draw links
    hivesvg.selectAll(".link")
        .data(links)
      .enter().append("path")
        .attr("class", "link")
        .attr("d", d3.hive.link()
          .angle(function(d) { return angle(d.x); })
          .radius(function(d) { return radius(d.y); }))
        .style("stroke", function(d) { return color(d.source.x); })
        .on("mouseover", linkMouseoverHP)
        .on("mouseout", mouseoutHP);

    // draw nodes
    hivesvg.selectAll(".node")
        .data(nodes)
      .enter().append("circle")
        .attr("class", "node")
        .attr("transform", function(d) { return "rotate(" + degreesHP(angle(d.x)) + ")"; })
        .attr("cx", function(d) { return radius(d.y); })
        .attr("r", 5)
        .style("fill", function(d) { return color(d.x); })
        .on("mouseover", nodeMouseoverHP)
        .on("mouseout", mouseoutHP);
}


// D3 hive extension
d3.hive = {}, d3.hive.link = function() {
    function t(t, s) {
        var u, h = a(r, this, t, s),
            i = a(n, this, t, s);
        h.a > i.a && (u = i, i = h, h = u), i.a - h.a > Math.PI && (h.a += 2 * Math.PI);
        var e = h.a + (i.a - h.a) / 3,
            c = i.a - (i.a - h.a) / 3;
        return h.r0 - h.r1 || i.r0 - i.r1 ? "M" + Math.cos(h.a) * h.r0 + "," + Math.sin(h.a) * h.r0 + "L" + Math.cos(h.a) * h.r1 + "," + Math.sin(h.a) * h.r1 + "C" + Math.cos(e) * h.r1 + "," + Math.sin(e) * h.r1 + " " + Math.cos(c) * i.r1 + "," + Math.sin(c) * i.r1 + " " + Math.cos(i.a) * i.r1 + "," + Math.sin(i.a) * i.r1 + "L" + Math.cos(i.a) * i.r0 + "," + Math.sin(i.a) * i.r0 + "C" + Math.cos(c) * i.r0 + "," + Math.sin(c) * i.r0 + " " + Math.cos(e) * h.r0 + "," + Math.sin(e) * h.r0 + " " + Math.cos(h.a) * h.r0 + "," + Math.sin(h.a) * h.r0 : "M" + Math.cos(h.a) * h.r0 + "," + Math.sin(h.a) * h.r0 + "C" + Math.cos(e) * h.r1 + "," + Math.sin(e) * h.r1 + " " + Math.cos(c) * i.r1 + "," + Math.sin(c) * i.r1 + " " + Math.cos(i.a) * i.r1 + "," + Math.sin(i.a) * i.r1
    }

    function a(t, a, r, n) {
        var e = t.call(a, r, n),
            c = +("function" == typeof s ? s.call(a, e, n) : s) + i,
            o = +("function" == typeof u ? u.call(a, e, n) : u),
            M = u === h ? o : +("function" == typeof h ? h.call(a, e, n) : h);
        return {
            r0: o,
            r1: M,
            a: c
        }
    }
    var r = function(t) {
            return t.source
        },
        n = function(t) {
            return t.target
        },
        s = function(t) {
            return t.angle
        },
        u = function(t) {
            return t.radius
        },
        h = u,
        i = -Math.PI / 2;
    return t.source = function(a) {
        return arguments.length ? (r = a, t) : r
    }, t.target = function(a) {
        return arguments.length ? (n = a, t) : n
    }, t.angle = function(a) {
        return arguments.length ? (s = a, t) : s
    }, t.radius = function(a) {
        return arguments.length ? (u = h = a, t) : u
    }, t.startRadius = function(a) {
        return arguments.length ? (u = a, t) : u
    }, t.endRadius = function(a) {
        return arguments.length ? (h = a, t) : h
    }, t
};

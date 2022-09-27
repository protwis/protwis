var diameter = 1200;

var tree = d3.layout.tree()
    .size([360, diameter / 2 - 50])
    .separation(function(a, b) { return (a.parent == b.parent ? 1 : 2) / a.depth; });

var diagonal = d3.svg.diagonal.radial()
    .projection(function(d) { return [d.y, d.x / 180 * Math.PI]; });

var svg = d3.select("#variantdata").append("svg")
    .attr("width", diameter)
    .attr("height", diameter)
    .attr("id", "drugmapping_svg")
    .attr("xmlns", "http://www.w3.org/2000/svg")

    svg.append("rect")
    .attr("width", "100%")
    .attr("height", "100%")
    .attr("fill", "white")

var svg_g = svg.append("g")
    .attr("transform", "translate(" + diameter / 2 + "," + diameter / 2 + ")");

var color = d3.scale.category10();

var data = variantdata;

var nodes = tree.nodes(data),
  links = tree.links(nodes);

var link = svg_g.selectAll(".link")
  .data(links)
.enter().append("path")

  // no color:
  .attr("style", function(d) {
    if (d.source.depth === 0) {return "fill: none; opacity: 0.8; stroke-width: " + 10 + "; stroke: #C0C0C0"}
    else if (d.source.depth === 1) {return "fill: none; opacity: 0.8; stroke-width: " + 4 + "; stroke: #C0C0C0"}
    else if (d.source.depth === 2) {return "fill: none; opacity: 0.8; stroke-width: " + 2 + "; stroke: #C0C0C0"}
    else if (d.source.depth === 3) {return "fill: none; opacity: 0.8; stroke-width: " + 1 + "; stroke: #C0C0C0"}})

  .attr("d", diagonal);

function Get(yourUrl){
  var Httpreq = new XMLHttpRequest(); // a new request
  Httpreq.open("GET",yourUrl,false);
  Httpreq.send(null);
  return Httpreq.responseText;}

var tip = d3.tip()
        .attr('class', 'd3-tip')
        .style("visibility","visible")
        .offset([0, 0])
        .html(function(d) {
          if (d.depth === 4) {return d.name.toUpperCase() + "<br>" + "<span style='color:#d7ded7'> Number of variants: " + d.number_of_variants + "<br>" + "<span style='color:#d7ded7'>Density: " + d.density_of_variants;}
          else {return d.name + "<br>" + "<span style='color:#d7ded7'> Number of variants: " + d.family_sum_trials + "<br>" + "<span style='color:#d7ded7'> Density: " + d.family_sum_approved;}
        });

tip(svg_g.append("g"));

var node = svg_g.selectAll(".node")
  .data(nodes)
.enter().append("g")
  .attr("class", "node")
  .attr("transform", function(d) { return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")"; })
  .on('mouseover', function(d) {
    tip.show(d)
    if (d.depth === 4) {d3.select(this).style("cursor", "pointer")}
  })
  .on('mouseout', function(d) {
    tip.hide(d)
    d3.select(this).style("cursor", "default")
  })

node.append("circle")
  .attr("r", function(d) {
    if (d.number_of_variants >= 1) {return 6.0}
    else {return 0.0} })
  .style("fill", function(d) {return "#43A047"})
  .style("fill-opacity", function(d) {return d.density_of_variants * 2});

node.append("circle")
  .attr("r", 1.0)
  .style("fill-opacity", 1.0)
  .style("fill", function(d) { return "black" })

node.append("text")
  .attr("dy", ".31em")
  .attr("text-anchor", function(d) { return d.x < 180 ? "start" : "end"; })
  .attr("transform", function(d) { return d.x < 180 ? "translate(8)" : "rotate(180)translate(-8)"; })
  .text(function(d) { if (d.depth==4) { return d.name.toUpperCase() ; } else if (d.depth>0) { return d.name;} else { return ""; } })
  // .style("font-size", function(d) { if (d.depth<3) {return "1px"} else if (d.depth==3) {return "1px"} else { return "1px" } })
  .style("font-family", "Palatino")
  .style("font", function(d) {
    if (d.depth === 1) {return "20px sans-serif"}
    else if (d.depth === 2) {return "16px sans-serif"}
    else if (d.depth === 3) {return "12px sans-serif"}
    else if (d.depth === 4) {return "9px sans-serif"}
    })

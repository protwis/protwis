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

var color = d3.scale.category20();

var diameter = 1350;

var tree = d3.layout.tree()
    .size([360, diameter / 2 - 120])
    .sort(function(a,b) { return b.sort - a.sort; })
    .separation(function(a, b) { return (a.parent == b.parent ? 1 : 2) / a.depth; });
    // .separation(function(a, b) { return ((a.parent == root) && (b.parent == root)) ? 3 : 1; });

var diagonal = d3.svg.diagonal.radial()
    .projection(function(d) { return [d.y, d.x / 180 * Math.PI]; });

var svg = d3.select("#coverage").append("svg")
    .attr("width", diameter)
    .attr("height", diameter)
    .attr("id", "coverage_svg")
    .attr("xmlns", "http://www.w3.org/2000/svg");
    
    svg.append("rect")
    .attr("width", "100%")
    .attr("height", "100%")
    .attr("fill", "white");

var svg_g = svg.append("g")
    .attr("transform", "translate(" + diameter / 2 + "," + diameter / 2 + ")");

    // svg.append("rect")
    // .attr("width", "100%")
    // .attr("height", "100%")
    // .attr("transform", "translate(" + diameter / 2 + "," + diameter / 2 + ")");


var nodes = tree.nodes(coverage)
      

    nodes.forEach(function(d) { 
    if (d.depth == 0) {
        d.y =  0
      } else if (d.depth == 1) {
        d.y =  130
      } else if (d.depth == 2  ) {
        d.y =  300
      } else if (d.depth == 3  ) {
        d.y =  550
      } else if (d.depth == 4  ) {
        d.y =  600
      }  else {
        d.y =  d.depth*150
      }
    });

var links = tree.links(nodes);

var link = svg_g.append("g")
      .attr("class", "links")
    .selectAll("path")
      .data(links)
    .enter().append("path")
      .each(function(d) { d.target.linkNode = this; })
      .attr("d", function(d) { return step(d.source.x, d.source.y, d.target.x, d.target.y) })
      .style("stroke", function(d) { return d.target.color; })
      .style("stroke-width", function(d) { if (d.target.depth>0) {return 5-d.target.depth;} else { return 0;} })
      .style("opacity", function(d) {
          if ((d.target.interactions > 0 && d.target.mutations_an > 0) || 1==1) {return 0.8 } //|| 1==1
          else if (d.target.interactions > 0) {return 0.5 }
          else if (d.target.mutations_an > 0) {return 0.5 }
          else {return 0.1}; });


 /// The more fluid version
  // var link = svg.selectAll(".link")
  //     .data(links)
  //     .enter().append("path")
  //     .attr("class", "link")
  //     .attr("d", diagonal)
  //     .style("stroke", function(d) { return d.target.color; })
  //     .style("stroke-width", function(d) { if (d.target.depth>0) {return 17-d.target.depth*d.target.depth;} else { return 0;} })
  //     .style("opacity", function(d) {
  //         if ((d.target.interactions > 0 && d.target.mutations_an > 0) || 1==1) {return 0.8 } //|| 1==1
  //         else if (d.target.interactions > 0) {return 0.5 }
  //         else if (d.target.mutations_an > 0) {return 0.5 }
  //         else {return 0.1}; })

  var node = svg_g.selectAll(".node")
      .data(nodes)
    .enter().append("g")
      .attr("class", "node")
      .attr("transform", function(d) { if (d.name=='GPCRs' ) {return "rotate(" + (d.x) + ")translate(" + d.y + ")";} else { return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")"; } })

  node.filter(function(d) {return (d.depth==4)}).append("circle")
      .attr("r", function(d) { if (d.name=='GPCRs' ) {return "0"} else { return "4.5" } })
      .style("fill", function(d) {
          if (d.color && d.depth<4) {return d.color } 
          // else if (d.interactions > 0 && d.mutations_an > 0) {return "Olive" }
          else if (d.interactions > 0) {return "FireBrick" }
          // else if (d.mutations_an > 0) {return "Chocolate" }
          else {return "#eee"}; })
      .style("opacity",.99);

  node.append("text")
      .attr("dy", ".31em")
      .attr("text-anchor", function(d) { return d.x < 180 ? "end" : "start"; })
      .attr("transform", function(d) { if (d.depth==4) {
            return d.x < 180 ? "translate(50)" : "rotate(180)translate(-50)";
        } else {
            return d.x < 180 ? "translate(-12)" : "rotate(180)translate(12)";
        }
         })
      .text(function(d) { if (d.depth==4) { return d.name.toUpperCase() ; } else if (d.depth>0) { return d.name;} else { return ""; } })
      // .style("font-size", function(d) { if (d.depth<3) {return "12px"} else { return "9px" } })

      .style("font-size", function(d) { if (d.depth<3) {return "16px"} else if (d.depth==3) {return "14px"} else { return "12px" } })
      .style("font-family", "Palatino")
      .style("fill", function(d) {
          if (d.color) {return "#111" } 
          else if (d.interactions > 0 && d.mutations_an > 0 && 1==2) {return "green" }
          else if (d.interactions > 0 && 1==2) {return "Olive" }
          else if (d.mutations_an > 0 && 1==2) {return "palegreen" }
          else {return "#222"}; }).call(getBB);
  node.filter(function(d) {return (d.depth!=4)}).insert("rect","text")
    .attr("x", function(d){return d.x < 180 ? d.bbox.x-12 : d.bbox.x-d.bbox.width-12; })
    .attr("y", function(d){return d.bbox.y})
    .attr("width", function(d){return d.bbox.width})
    .attr("height", function(d){return d.bbox.height})
    .style("fill", "#FFF");

function getBB(selection) {
    selection.each(function(d){d.bbox = this.getBBox();})
}
 
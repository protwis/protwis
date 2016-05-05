function parseNewick(a){for(var e=[],r={},s=a.split(/\s*(;|\(|\)|,|:)\s*/),t=0;t<s.length;t++){var n=s[t];switch(n){case"(":var c={};r.branchset=[c],e.push(r),r=c;break;case",":var c={};e[e.length-1].branchset.push(c),r=c;break;case")":r=e.pop();break;case":":break;default:var h=s[t-1];")"==h||"("==h||","==h?r.name=n:":"==h&&(r.length=parseFloat(n))}}return r}

phylo_tm = "((((((CaS:0.43546,GPRC6:0.43546)blah2:0.36894,((mGlu1:0.13746,mGlu5:0.13746):0.26302,((mGlu2:0.16096,mGlu3:0.16096):0.17358,((mGlu4:0.09402,(mGlu7:0.06962,mGlu8:0.06962):0.02440):0.05175,mGlu6:0.14577):0.18877):0.06594)mGlu:0.40392):0.31168,((TAS1R1:0.77719,TAS1R2:0.77719):0.22828,TAS1R3:1.00548):0.11060)Nexttest:0.14325,((GPRC5A:0.34164,GPRC5D:0.34164):0.17194,(GPCR5B:0.40816,GPRC5C:0.40816):0.10542)blah:0.74575):0.10652,((GABAb1:0.50615,GABAb2:0.50615):0.32881,GPR156:0.83496)Eukaryota:0.53089):0.39150,(GPR158:0.33108,GPR179:0.33108)Archaea:1.42627);";


var outerRadius = 900 / 2,
    innerRadius = outerRadius - 100;

var color = d3.scale.category10()
    .domain(["Bacteria", "Eukaryota", "Archaea","Nexttest","mGlu","blah","blah2"]);

var cluster = d3.layout.cluster()
    .size([360, innerRadius])
    .children(function(d) { return d.branchset; })
    .value(function(d) { return 1; })
    .sort(function(a, b) { return (a.value - b.value) || d3.ascending(a.length, b.length); })
    .separation(function(a, b) { return 1; });


var diagonal2 = d3.svg.diagonal.radial()
    .projection(function(d) { return [d.y, d.x / 180 * Math.PI]; });

var diagonal3 = d3.svg.diagonal()
    .source(function(d) { return {"x":d.target.x, "y":d.target.y}; })            
    .target(function(d) { return {"x":d.target.x, "y":d.source.y}; })
    .projection(function(d) { return [d.y, d.x / 180 * Math.PI]; });

var diagonal4 = d3.svg.diagonal.radial()
    .projection(function(d) { var r = d.y, a = (d.x - 90) / 180 * Math.PI;
                              return [r * Math.cos(a), r * Math.sin(a)]; });    

var svg = d3.select("#svg_tree").append("svg")
    .attr("width", outerRadius * 2)
    .attr("height", outerRadius * 2);

// var legend = svg.append("g")
//     .attr("class", "legend")
//   .selectAll("g")
//     .data(color.domain())
//   .enter().append("g")
//     .attr("transform", function(d, i) { return "translate(" + (outerRadius * 2 - 10) + "," + (i * 20 + 10) + ")"; });

// legend.append("rect")
//     .attr("x", -18)
//     .attr("width", 18)
//     .attr("height", 18)
//     .style("fill", color);

// legend.append("text")
//     .attr("x", -24)
//     .attr("y", 9)
//     .attr("dy", ".35em")
//     .style("text-anchor", "end")
//     .text(function(d) { return d; });

var chart = svg.append("g")
    .attr("transform", "translate(" + outerRadius + "," + outerRadius + ")");


  var root = parseNewick(phylo_tm),
      nodes = cluster.nodes(root),
      links = cluster.links(nodes),
      input = d3.select("#show-length input").on("change", changed),
      timeout = setTimeout(function() { input.property("checked", true).each(changed); }, 2000);

  setRadius(root, root.length = 0, innerRadius / maxLength(root));
  var total = totalFinalChildren(root);
  var arc_size = 360/total;
  setColor(root);

  var linkExtension = chart.append("g")
      .attr("class", "link-extensions")
    .selectAll("path")
      .data(links.filter(function(d) { return !d.target.children; }))
    .enter().append("path")
      .each(function(d) { d.target.linkExtensionNode = this; })
      .attr("d", function(d) { return step(d.target.x, d.target.y, d.target.x, innerRadius); });
      // .attr("d", diagonal4); 

  var link = chart.append("g")
      .attr("class", "links")
    .selectAll("path")
      .data(links)
    .enter().append("path")
      .each(function(d) { d.target.linkNode = this; })
      .attr("d", function(d) { return step(d.source.x, d.source.y, d.target.x, d.target.y) })
      // .attr("d", diagonal4)
      .style("stroke", function(d) { return d.target.color; })
      .style("stroke-width", "2px");

  chart.append("g")
      .attr("class", "labels")
    .selectAll("text")
      .data(nodes.filter(function(d) { return !d.children; }))
    .enter().append("text")
      .attr("dy", ".31em")
      .attr("transform", function(d) { return "rotate(" + (d.x - 90) + ")translate(" + (innerRadius + 4) + ",0)" + (d.x < 180 ? "" : "rotate(180)"); })
      .style("text-anchor", function(d) { return d.x < 180 ? "start" : "end"; })
      // .text(function(d) { return d.name.replace(/_/g, " "); })
      .on("mouseover", mouseovered(true))
      .on("mouseout", mouseovered(false));

  var arc = d3.svg.arc()
    .innerRadius(function(d) {return innerRadius; })
    .outerRadius(innerRadius+20)
    .startAngle(function(d) { return Math.round(d.x-arc_size/2+1) * (Math.PI/180);}) //converting from degs to radians
    .endAngle(function(d) { return  Math.round(d.x+arc_size/2-1)* (Math.PI/180);}) //just radians
    //   .startAngle(function(d) { return Math.max(0, Math.min(2 * Math.PI, x(d.x))); })
    //   .endAngle(function(d) { return Math.max(0, Math.min(2 * Math.PI, x(d.x + d.dx))); })

  chart.append("g")
      .attr("class", "labels")
    .selectAll("text")
      .data(nodes.filter(function(d) { return !d.children; }))
    .enter().append("path")
    .attr("class", "donutArcs")
    .attr("d", arc)
    .attr("fill",function(d) { return d.color; })
    .attr("opacity",function(d){  return (d.name=="mGlu1" || d.name=="mGlu5" ? "0.9" : "0.9") ;}).each(function(d,i) {
      //Search pattern for everything between the start and the first capital L
      var firstArcSection = /(^.+?)L/;  

      //Grab everything up to the first Line statement
      var newArc = firstArcSection.exec( d3.select(this).attr("d") )[1];
      //Replace all the comma's so that IE can handle it
      newArc = newArc.replace(/,/g , " ");
      
      //If the end angle lies beyond a quarter of a circle (90 degrees or pi/2) 
      //flip the end and start position
      if (d.x > 90 && d.x<270) {
        var startLoc  = /M(.*?)A/,    //Everything between the first capital M and first capital A
          middleLoc   = /A(.*?)0 0 1/,  //Everything between the first capital A and 0 0 1
          endLoc    = /0 0 1 (.*?)$/; //Everything between the first 0 0 1 and the end of the string (denoted by $)
        //Flip the direction of the arc by switching the start en end point (and sweep flag)
        //of those elements that are below the horizontal line
        var newStart = endLoc.exec( newArc )[1];
        var newEnd = startLoc.exec( newArc )[1];
        var middleSec = middleLoc.exec( newArc )[1];
        
        //Build up the new arc notation, set the sweep-flag to 0
        newArc = "M" + newStart + "A" + middleSec + "0 0 0 " + newEnd;
      }//if
      
      //Create a new invisible arc that the text can flow along
      chart.append("path")
        .attr("class", "hiddenDonutArcs")
        .attr("id", "donutArc"+i)
        .attr("d", newArc)
        .style("fill", "none");
    });

    chart.append("g")
    .attr("class", "labels")
    .selectAll("text")
    .data(nodes.filter(function(d) { return !d.children; }))
    .enter().append("text")
      // .attr("class", "donutText")
      .attr("class", "labels")
      // .attr("x", Math.round(arc_size) * (Math.PI/180)*innerRadius*0.5) //Move the text from the start angle of the arc
      // .attr("dy", 15) //Move the text down
      .attr("dy", function(d) { return (d.x > 90 && d.x<270 ? -5 : +15); })
      // .attr("transform", function(d) { console.log(d.x +" " + d.name); return (d.x < 270 && d.x>90 ? "translate(180,180)" : ""); })

      .append("textPath")
      .attr("fill", function(d){  return (d.name=="mGlu1123" ? "#000" : "#FFF") ;})

      .attr("text-decoration", function(d){  return (d.name=="mGlu1" ? "" : "") ;})
      .attr("startOffset","50%")  
      .style("text-anchor", "middle")
      .attr("xlink:href",function(d,i){return "#donutArc"+i;})
      .text(function(d){  return d.name;})
      .on("mouseover", mouseovered(true))
      .on("mouseout", mouseovered(false));

    //   .append("text")
    // .attr("class", "donutText")
    // //Move the labels below the arcs for those slices with an end angle greater than 90 degrees
    // .attr("dy", function(d,i) { return (d.endAngle > 90 * Math.PI/180 ? 18 : -11); })
    //  .append("textPath")
    // .attr("startOffset","50%")
    // .style("text-anchor","middle")
    // .attr("xlink:href",function(d,i){return "#"+d.name;})
    // .text(function(d){return d.name;});


  function changed() {
    clearTimeout(timeout);
    var checked = this.checked;
    d3.transition().duration(750).each(function() {
      linkExtension.transition().attr("d", function(d) { return step(d.target.x, checked ? d.target.radius : d.target.y, d.target.x, innerRadius); });
      link.transition().attr("d", function(d) { return step(d.source.x, checked ? d.source.radius : d.source.y, d.target.x, checked ? d.target.radius : d.target.y) });
      // link.transition().attr("d", diagonal4)
    });
  }

  function mouseovered(active) {
    return function(d) {
      d3.select(this).classed("label--active", active);
      d3.select(d.linkExtensionNode).classed("link-extension--active", active).each(moveToFront);
      do d3.select(d.linkNode).classed("link--active", active).each(moveToFront); while (d = d.parent);
    };
  }

  function moveToFront() {
    this.parentNode.appendChild(this);
  }

// Compute the maximum cumulative length of any node in the tree.
function maxLength(d) {
  return d.length + (d.children ? d3.max(d.children, maxLength) : 0);
}

// Compute the maximum cumulative length of any node in the tree.
function totalFinalChildren(d) {
  if (d.children){
    var count = 0
    d.children.forEach(function(d) { count += totalFinalChildren(d); }); 
    return count;
  } else {
    return 1;
  }
}

// Set the radius of each node by recursively summing and scaling the distance from the root.
function setRadius(d, y0, k) {
  d.radius = (y0 += d.length) * k;
  if (d.children) d.children.forEach(function(d) { setRadius(d, y0, k); });
}

// Set the color of each node by recursively inheriting.
function setColor(d) {
  d.color = color.domain().indexOf(d.name) >= 0 ? color(d.name) : d.parent ? d.parent.color : null;
  if (d.children) d.children.forEach(setColor);
}

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
var hivesvg, hiveTooltip;
function createHiveplot(data, container) {
    var width = 800,
        height = 1000,
        innerRadius = 0.15*(height),
        outerRadius = 0.36*(height);

    var angleLineCore = d3.scale.ordinal()
                      .domain(d3.range(71))
                      .rangePoints([2 * Math.PI, 0]),
        radius = d3.scale.linear()
                  .range([innerRadius, outerRadius]),
        color = d3.scale.category10()
                  .domain(d3.range(10));

      var angle = function(e){
        // TODO optimize placement of TMs based on angle, but possibly also starting position
        if (e >= 6.5)
          e = 6.5 + (e-7.0);
        else if (e < 0){
          e += 7;
        }

        return angleLineCore( (e*10+18) % 70 );
      };

    // create TM residue count
    // replace later with data.segments
    var TM_segments = ["TM1", "TM2", "TM3", "TM4", "TM5", "TM6", "TM7", "H8"];
    var resCount = [];
    for (var i = 0; i < 8; i++)
        resCount[i] = 0;
    var currentCount = resCount.slice(0);

    // count residues in each segment
    for (key in data.segment_map){
        var index = TM_segments.indexOf(data.segment_map[key]);
        if (index >= 0)
            resCount[index] += 1;
    }

    // x is TM number - 1 and y is ordered position in TM
    var nodes = [];
    var indexList = Object.assign({}, data.generic_map);
          var last = -1;
    for (residue in data.generic_map){
      // residue contains original number
      // data.generic_map[residue] contains GN
      var gn = data.generic_map[residue].split('x');
      var tm = gn[0];
      var position = gn[1];
      if (tm <= 8){
          currentCount[tm-1] += 1;
          // Reverse upward helices so EC in inside
          if (tm == 8 || (tm % 2) == 1)
              nodes.push({x: tm-1, y: currentCount[tm-1], gn: data.generic_map[residue], number: residue});
          else
              nodes.push({x: tm-1, y: resCount[tm-1]-currentCount[tm-1]+1, gn: data.generic_map[residue], number: residue});
          indexList[residue] = nodes.length - 1;
          //console.log("{name: \"" + data.generic_map[residue] + "\", group:" + tm + "},");
//          if (last == tm)
//              console.log("{\"source\": " + (nodes.length - 2) + ", \"target\":" + (nodes.length - 1) + ", \"value\": 0},");
          last =  tm;
      } else
          indexList[residue] = null;
    }

    // create links between nodes
    var links = [];
    for (contact in data.interactions){
        var pair = contact.split(',');
        var res1 = pair[0];
        var res2 = pair[1];
        // Create link if nodes will be drawn
        if (indexList[res1] != null && indexList[res2] != null) {
            // only is consecutive
            var tmDistance = Math.abs(nodes[indexList[res1]].x - nodes[indexList[res2]].x);
            if (tmDistance > 0) {
                if (tmDistance == 1 || (tmDistance == 6 && nodes[indexList[res2]].x != 7) || tmDistance == 7)
                    links.push({source: nodes[indexList[res1]], target: nodes[indexList[res2]], consecutive: true, tmDistance: tmDistance});
                else
                    links.push({source: nodes[indexList[res1]], target: nodes[indexList[res2]], consecutive: false, tmDistance: tmDistance});
            }
//            console.log("{\"source\": " + (indexList[res1]) + ", \"target\":" + (indexList[res2]) + ", \"value\": 5},");
        }
    }

    hiveTooltip = d3.select("body").append("div")
                  .style("display", "none")
                  .style("position", "absolute")
//                  .attr("class", "rounded")
                  .attr("id", "hiveTooltip");

    hivesvg = d3.select(container).append("svg")
        .attr("class", "hiveplot")
        // .attr("width", width)
        .attr("width", "100%")
        // .attr("height", height)
        // .attr('viewBox', "0 0 "+width+" "+height)
      .append("g")
        .attr("transform", "translate(" + width/1.6 + "," + height/2.2 + ")");

    // draw lines based on length TM helix
    hivesvg.selectAll(".hive-axis")
        .data(d3.range(8))
      .enter().append("line")
        .attr("class", "hive-axis")
        .attr("transform", function(d) { return "rotate(" + degreesHP(angle(d)) + ")" })
        .attr("x1", function(d) { if (resCount[d] == 0){ return 0; } else if (d == 7) { return (resCount[0] * 1.2) * (outerRadius-innerRadius)/20+innerRadius; } else { return innerRadius; }})
        .attr("x2", function(d) { var extra = 0; if (resCount[d] == 0){ return 0; } else if (d == 7  && resCount[d]>0){ extra = resCount[0] * 1.2}; return (outerRadius-innerRadius)/20*(extra+resCount[d]+1)+innerRadius; });
        //.attr("x2", function(d) { var extra = 0; if (d.x == 7){ extra = resCount[0] * 1.2}; return (outerRadius-innerRadius)/20*(resCount[d]+extra+1)+innerRadius; });

    // draw links between nodes
    /*
    hivesvg.selectAll(".link")
        .data(links)
      .enter().append("path")
        .attr("class", "link")
        .attr("d", d3.hive.link()
          .angle(function(d) { return angle(d.x); })
          .radius(function(d) { return d.y*(outerRadius-innerRadius)/20+innerRadius; }))
        .style("stroke", function(d) { if (d.tmDistance == 6) return color(d.target.x); else return color(d.source.x); })
        .on("mouseover", linkMouseoverHP)
        .on("mouseout", mouseoutHP);*/
    var lineCore = d3.svg.line.radial()
        .interpolate( "bundle" )
        .tension( 0.85 )
//        .radius(function(d) { return d.y*(outerRadius-innerRadius)/20+innerRadius; })
        .radius(function(d) { var extra = 0; if (d.x >= 6.5){ extra = resCount[0] * 1.2}; return (d.y+extra)*(outerRadius-innerRadius)/20+innerRadius; })
        .angle(function(d) { /*if (d.x < 0) d.x = 7 + d.x;*/ return angle(d.x); });

    hivesvg.selectAll(".hive-link")
        .data(links)
      .enter().append("path")
        .attr("class", "hive-link")
        .attr("d", function(d, i) {
            var placer = 0.3;
//            if (d.source.x == 7 || d.target.x == 7)
//                placer = placer/3;

            var margin = 10;
            if (d.tmDistance > 3)
                placer *= -1;
            var path = [d.source];
            path.push({x: d.source.x+placer, y: d.source.y});

            if (!d.consecutive) {
                var middle;
                var inside = true;
                if (d.tmDistance == 2 || d.tmDistance >= 5)
                    inside = false;

                var via = 0;
                if (d.tmDistance > 3)
                    //via = d.target.x + (7-d.tmDistance)/2;
                    via = d.source.x - (7-d.tmDistance)/2;
                else
                    via = d.source.x + d.tmDistance/2;
//                via = via % 7;

                // draw nice path outside or inside
                if (inside){
                    middle = {x: via, y: -1*margin/2};
                    path.push({x: middle.x - placer, y: middle.y + margin/3});
                    path.push(middle);
                    path.push({x: middle.x + placer, y: middle.y + margin/3});
                } else {
                    if (via == 0.5)
                      via = 0;
                    middle = {x: via, y: resCount[via] + margin };
                    path.push({x: middle.x - placer, y: middle.y - margin/2});
                    path.push(middle);
                    path.push({x: middle.x + placer, y: middle.y - margin/2});
                }
            }
            path.push({x: d.target.x-placer, y: d.target.y});
            path.push(d.target);

            return lineCore(path);
          })
        //.attr("d", d3.hive.link()
        //  .angle(function(d) { return angle(d.x); })
        //  .radius(function(d) { return d.y*(outerRadius-innerRadius)/20+innerRadius; }))
        .style("stroke", function(d) { if (d.tmDistance == 6) return color(d.target.x); else return color(d.source.x); })
        .on("mouseover", linkMouseoverHP)
        .on("mouseout", mouseoutHP);



    // draw nodes
    hivesvg.selectAll(".hive-node")
        .data(nodes)
      .enter().append("circle")
        .attr("class", "hive-node")
        .attr("transform", function(d) { return "rotate(" + degreesHP(angle(d.x)) + ")"; })
        .attr("cx", function(d) { var extra = 0; if (d.x == 7){ extra = resCount[0] * 1.2}; return (d.y+extra)*(outerRadius-innerRadius)/20+innerRadius; })
        .attr("r", 5)
        .style("fill", function(d) { return color(d.x); })
        .on("mouseover", nodeMouseoverHP)
        .on("mouseout", mouseoutHP)
// fix nice bootstrap tooltips
//        .each(function(d,i){ console.log("enable for "+d.gn);$(this).tooltip({title: d.gn + " (" + d.number + ")", html: false, placement: 'bottom', container: 'body'});});
}


// Based on Mike Bostock's example: https://bl.ocks.org/mbostock/2066415
// highlight link and connected nodes on mouseover
function linkMouseoverHP(d) {
  hivesvg.selectAll(".hive-link")
    .classed("turnedOn", function(dl) {
      return dl === d;
    })
    .classed("turnedOff", function(dl) {
      return !(dl === d);
    })
  hivesvg.selectAll(".hive-node")
    .classed("turnedOn", function(dl) {
      return dl === d.source || dl === d.target;
    })
}

// highlight node and connected links on mouseover
function nodeMouseoverHP(d) {
  hivesvg.selectAll(".hive-link")
    .classed("turnedOn", function(dl) {
      return dl.source === d || dl.target === d;
    })
    .classed("turnedOff", function(dl) {
      return !(dl.source === d || dl.target === d)
    });
  d3.select(this)
    .classed("turnedOn", true);

  // show hiveTooltip
  hiveTooltip.style("display", "block")
               .style("left",(d3.event.pageX+20)+"px")
               .style("top",(d3.event.pageY+10)+"px")
               .html(d.gn + " (" + d.number + ")");
}

// clear highlighted nodes or links
function mouseoutHP() {
  hivesvg.selectAll(".turnedOn").classed("turnedOn", false);
  hivesvg.selectAll(".turnedOff").classed("turnedOff", false);
  // hide hiveTooltip
  hiveTooltip.style("display", "none");
}

function degreesHP(radians) {
  return radians / Math.PI * 180 - 90;
}

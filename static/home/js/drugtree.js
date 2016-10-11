var diameter = 1200;

var tree = d3.layout.tree()
    .size([360, diameter / 2 - 50])
    .separation(function(a, b) { return (a.parent == b.parent ? 1 : 2) / a.depth; });

var diagonal = d3.svg.diagonal.radial()
    .projection(function(d) { return [d.y, d.x / 180 * Math.PI]; });


var svg = d3.select("#drugdata").append("svg")
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

var data = drugdata;

var nodes = tree.nodes(data),
  links = tree.links(nodes);

var link = svg_g.selectAll(".link")
  .data(links)
.enter().append("path")

  // .attr("style", function(d) { 
  //   if (d.source.depth === 0) {return "fill: none; opacity: 0.6; stroke-width: 10.0px; stroke: " + color(d.target.group)}
  //   else if (d.source.depth === 1) {return "fill: none; opacity: 0.6; stroke-width: 5.0px; stroke: " + color(d.source.group)}
  //   else if (d.source.depth === 2) {return "fill: none; opacity: 0.6; stroke-width: 2.5px; stroke: " + color(d.source.parent.group)}
  //   else if (d.source.depth === 3) {return "fill: none; opacity:" + (0.6 - (1/((d.target.trials + d.target.approved) + 2))) + "; stroke-width: 1.0px; stroke: " + color(d.source.parent.parent.group)};})

  // .attr("style", function(d) { 
  //   if (d.source.depth === 0) {return "fill: none; opacity: 0.25; stroke-width: " + 10 + "; stroke: " + color(d.target.group)}
  //   else if (d.source.depth === 1) {return "fill: none; opacity: 0.25; stroke-width: " + 4 + "; stroke: " + color(d.source.group)}
  //   else if (d.source.depth === 2) {return "fill: none; opacity: 0.25; stroke-width: " + 2 + "; stroke: " + color(d.source.parent.group)}
  //   else if (d.source.depth === 3) {return "fill: none; opacity: 0.25; stroke-width: " + 1 + "; stroke: " + color(d.source.parent.parent.group)};})

  // no color:
  .attr("style", function(d) { 
    if (d.source.depth === 0) {return "fill: none; opacity: 0.8; stroke-width: " + 10 + "; stroke: #C0C0C0"}
    else if (d.source.depth === 1) {return "fill: none; opacity: 0.8; stroke-width: " + 4 + "; stroke: #C0C0C0"}
    else if (d.source.depth === 2) {return "fill: none; opacity: 0.8; stroke-width: " + 2 + "; stroke: #C0C0C0"}
    else if (d.source.depth === 3) {return "fill: none; opacity: 0.8; stroke-width: " + 1 + "; stroke: #C0C0C0"};})

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
          if (d.depth === 4) {return d.name.toUpperCase() + "<br>" + "<span style='color:#43A047'> Trials: " + d.trials + "<br>" + "<span style='color:#d62728'> Approved: " + d.approved;}
          else {return d.name + "<br>" + "<span style='color:#43A047'> Trials: " + d.family_sum_trials + "<br>" + "<span style='color:#d62728'> Approved: " + d.family_sum_approved;}                  
        });
                        
tip(svg_g.append("g"));

// LEGEND


var t = $('#clickdata').DataTable({
            "scrollX": true,
            // 'scrollY': $(window).height()/5,
            'paging': false,
            'autoWidth': true,
            'bScrollCollapse': true,
            'orderCellsTop': true,
            'dom': 'T<"clear">lfrtip',
            "aoColumns": [ 
                        {"sClass": "center"},
                        {"sClass": "center"},
                        {"sClass": "center"},
                        {"sClass": "center"},
                        {"sClass": "center"},
                        {"sClass": "center"}],
            // 'aoColumnDefs': [],
            "order": [[ 2, "asc" ], [ 4, "desc" ]],
            'tableTools': {
                "sRowSelect": "multi",
                "aButtons": []
            },
            "language": {
            "zeroRecords": "No data available in table - please click on a given receptor name to load drug data.",
            "infoEmpty": "No records available"
            },
                initComplete: function () {
                    $('#clickdata').dataTable().columnFilter()
                }
            // ,
            // initComplete: function () {
            //     $('#clickdata').dataTable().columnFilter({
            //         sPlaceHolder: "head:after",
            //         aoColumns: [
            //             { type: "text" },
            //             { type: "text" },
            //             { type: "select" },
            //             { type: "select" },
            //             { type: "select" },
            //             { type: "select" },
            //         ]
            //     });
            // }

        });

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
  .on("click",  function(d) {
    url = "/services/drugs/"+d.name+"_human"
    var json_obj = JSON.parse(Get(url))

    var tableCap = document.getElementById("caption");
    tableCap.innerHTML = "<h5>Drug details for: <a href='/protein/"+d.name+"_human'>"+d.name.toUpperCase()+"</h5>" ;
    t.clear();

    for (object in json_obj){
      t.row.add( [
              json_obj[object]['name'],
              json_obj[object]['indication'],
              json_obj[object]['status'],
              json_obj[object]['approval'],
              json_obj[object]['drugtype'],
              json_obj[object]['novelty']
          ] ).draw( false );
    }

    t.draw();

    });

var xlegend = 370;
var ylegend = -520;

var legend1 = svg_g.selectAll("circle")
.data([ylegend, ylegend-40]);

var legendEnter = legend1.enter();

legendEnter.append("circle")
            .attr("cy", function(d) { return d; })
           .attr("cx", xlegend)
           .attr("r", 10)
           .style("fill", function(d) { 
      if (d == ylegend) {return "#d62728"}
      else {return "#43A047"};});

legendEnter.append("circle")
           .attr("cy", function(d) { return d; })
           .attr("cx", xlegend)
           .attr("r", 1.5)
           .style("fill", "black");

legendEnter.append("text")
            .attr("x", xlegend + 30)
            .attr("y", ylegend + 5)         
            .style({"font": "15lspx sans-serif"})
            .text("Target established");

legendEnter.append("text")
            .attr("x", xlegend + 30)
            .attr("y", ylegend - 40 + 5)         
            .style({"font": "15lspx sans-serif"})
            .text("Target in trials");

node.append("circle")
  .attr("r", function(d) { 
    if (d.trials >= 1 || d.approved >= 1 ) {return 6.0}
    else {return 0.0} ;})
  .style("fill-opacity", 1.0)
  .style("fill", function(d) { 
    if (d.establishment == 4) {return "#d62728"}
    else {return "#43A047"};});

// node.append("circle")
//     .attr("r", function(d) { return d.approved})
//     // .attr("r", 3.0)
//     .style("fill-opacity", 0.50)
//     // .style("fill", function(d) { return color(d.establishment); })
//     .style("fill", function(d) { return "green" })

// node.append("circle")
//     .attr("r", function(d) { return d.trials})
//     // .attr("r", 3.0)
//     .style("fill-opacity", 0.50)
//     .style("fill", function(d) { return "red" })

node.append("circle")
  // .attr("r", function(d) {if (d.depth === 3) return  d.family_sum_approved})
  .attr("r", function(d) {
  if (d.depth === 3 && d.family_sum_approved >= 1 && d.family_sum_approved <= 5) {return 5.0}
  else if (d.depth === 3 && d.family_sum_approved > 5 && d.family_sum_approved <= 10) {return 8.0}
  else if (d.depth === 3 && d.family_sum_approved > 10 && d.family_sum_approved <= 20) {return 15.0}
  else if (d.depth === 3 && d.family_sum_approved > 20 && d.family_sum_approved <= 30) {return 20.0}
  else if (d.depth === 3 && d.family_sum_approved > 30 && d.family_sum_approved <= 40) {return 25.0}
  else if (d.depth === 3 && d.family_sum_approved > 40 && d.family_sum_approved <= 50) {return 30.0}
  else if (d.depth === 3 && d.family_sum_approved > 50 && d.family_sum_approved <= 75) {return 40.0}
  else if (d.depth === 3 && d.family_sum_approved > 75 ) {return 50.0}
  else {return 0.0} ;})
  // .attr("r", 3.0)

  .style("fill-opacity", 0.70)
  // .style("fill", function(d) { return color(d.establishment); })
  .style("fill", function(d) { return "#d62728" })

node.append("circle")
  // .attr("r", function(d) {if (d.depth === 3) return d.family_sum_trials})
  .attr("r", function(d) {
  if (d.depth === 3 && d.family_sum_trials >= 1 && d.family_sum_trials <= 5) {return 5.0}
  else if (d.depth === 3 && d.family_sum_trials > 5 && d.family_sum_trials <= 10) {return 8.0}
  else if (d.depth === 3 && d.family_sum_trials > 10 && d.family_sum_trials <= 20) {return 15.0}
  else if (d.depth === 3 && d.family_sum_trials > 20 && d.family_sum_trials <= 30) {return 20.0}
  else if (d.depth === 3 && d.family_sum_trials > 30 && d.family_sum_trials <= 40) {return 25.0}
  else if (d.depth === 3 && d.family_sum_trials > 40 && d.family_sum_trials <= 50) {return 30.0}
  else if (d.depth === 3 && d.family_sum_trials > 50 && d.family_sum_trials <= 75) {return 40.0}
  else if (d.depth === 3 && d.family_sum_trials > 75 ) {return 50.0}
  else {return 0.0} ;})
  // .attr("r", 3.0)
  .style("fill-opacity", 0.60)
  // .style("fill", function(d) { return color(d.establishment); })
  .style("fill", function(d) { return "#43A047" })

node.append("circle")
  .attr("r", 1.0)
  .style("fill-opacity", 1.0)
  .style("fill", function(d) { return "black" })

node.append("text")
  .attr("dy", ".31em")
  .attr("text-anchor", function(d) { return d.x < 180 ? "start" : "end"; })
  .attr("transform", function(d) { return d.x < 180 ? "translate(8)" : "rotate(180)translate(-8)"; })
  .text(function(d) { if (d.depth==4) { return d.name.toUpperCase() ; } else if (d.depth>0) { return d.name;} else { return ""; } })
  .style("font-size", function(d) { if (d.depth<3) {return "16px"} else if (d.depth==3) {return "14px"} else { return "12px" } })
  .style("font-family", "Palatino")
  .style("font", function(d) {
    if (d.depth === 1) {return "16px sans-serif"} 
    else if (d.depth === 2) {return "15px sans-serif"}
    else if (d.depth === 3) {return "12px sans-serif"}
    else if (d.depth === 4) {return "8px sans-serif"}
    })
  // .style("fill", function(d) {
  //     if (d.crystal_structure > 0) {return "#00B8D4" }
  //     else {return "black"}; })
  .style("font-weight", function(d) {
      if (d.crystal_structure > 0) {return "bold" }
      else {return "padding"}; })



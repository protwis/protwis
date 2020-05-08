function renderHeatmap_distances(data, heatMapSelector) {

    var margin = {top: 80, right: 0, bottom: 30, left: 80},
        width = 500,
        height = 500;
    var margin = { top: 20, right: 30, bottom: 40, left: 30 },
    width = 520 - margin.left - margin.right,
    height = 520 - margin.top - margin.bottom;

    var x = d3.scale.ordinal().rangeBands([0, width]),
        z = d3.scale.linear().domain([0, 4]).clamp(true);

    var svg = d3v4.select(heatMapSelector).select(".heatmap_distances")
        .attr("viewBox", 0 + " " + 0 + " " + (width+ margin.left + margin.right) + " " + (height+ margin.top + margin.bottom))
        .attr("width", "100%")
        .attr("style", "height: 500px")
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    var matrix = [], nodes = [],
        distances = data['distances'];
    var dis_min = 0, dis_max = 0;
      
    var gn_index = data.segm_lookup_ordered;

    // Remove empty distances
    filtered = gn_index.filter(function (value, index, arr) { return distances[value];});
    gn_index = filtered;

    var segments = [];
    var segment_start = 0;
    var previous_segment = data['segment_map'][gn_index[0]];
    $.each(gn_index, function (i, gn) {
        seg = data['segment_map'][gn];
        if (seg != previous_segment || i == gn_index.length-1) {
            segments.push({ segment: previous_segment, start: segment_start, end: i })
            segment_start = i;
            previous_segment = seg;
        }
        nodes.push({ name: gn, segment: seg });
    })
    var n = nodes.length;
    label_font_size = Math.round(1*n/20);
    
    
  // Compute index per node.
  nodes.forEach(function(node, i) { 
      node.index = i;
      node.count = 0;
      matrix[i] = d3.range(n).map(function (j) { return { x: j, y: i, z: 0 }; });
      
      counts = distances[node.name]
      $.each(counts, function (target, value) { 
          if (target !='avg') {
              var j = gn_index.indexOf(target);
              if (j==-1) console.log('ERRORR',node.name,target)
              matrix[i][j].z = value;
              dis_min = value < dis_min ? value : dis_min;
              dis_max = value > dis_max ? value : dis_max;
            }
      })

  });

  // Precompute the orders.
  var orders = {
    init: d3.range(n),
    name: d3.range(n).sort(function(a, b) { return d3.ascending(nodes[a].name, nodes[b].name); }),
    count: d3.range(n).sort(function(a, b) { return nodes[b].count - nodes[a].count; }),
    group: d3.range(n).sort(function(a, b) { return nodes[b].group - nodes[a].group; })
  };
    
  // Color scale.
  c = d3.scale.linear().domain([dis_min, 0, dis_max]).range(['red', 'white', 'blue']);
  // The default sort order.
  x.domain(orders.init);

  svg.append("rect")
      .attr("class", "background")
      .attr("fill", "white")
      .attr("width", width)
      .attr("height", height);

  var row = svg.selectAll(".row")
      .data(matrix)
    .enter().append("g")
      .attr("class", "row")
      .attr("transform", function(d, i) { return "translate(0," + x(i) + ")"; })
      .each(row);

  row.append("line")
      .attr("x2", width);
    
  var row_labels = svg.selectAll(".row_label")
        .data(segments)
    .enter().append("text")
    .attr("x", -6)
      .attr("y", function (d) { return x(Math.round((d.start+d.end)/2)); })
    .attr("dy", ".32em")
    .attr("text-anchor", "end")
    .attr("font-size",label_font_size)
    .text(function(d, i) { return d.segment; })
    .attr("transform", function(d, i) { return "rotate(-90,-6,"+x(Math.round((d.start+d.end)/2))+")"; });

  var column = svg.selectAll(".column")
      .data(matrix)
    .enter().append("g")
      .attr("class", "column")
      .attr("transform", function(d, i) { return "translate(" + x(i) + ")rotate(-90)"; });

  column.append("line")
      .attr("x1", -width);

  var cell_labels = svg.selectAll(".cell_label")
        .data(segments)
        .enter().append("text")
        .attr("y", -6)
        .attr("x", function (d) { return x(Math.round((d.start + d.end) / 2)); })
        .attr("dy", ".32em")
        .attr("text-anchor", "middle")
        .attr("font-size", label_font_size)
        .text(function (d, i) { return d.segment; });

  function row(row) {
    var cell = d3.select(this).selectAll(".cell")
      .data(row.filter(function (d) { return d.z; }))
      .enter().append("rect")
      .attr("class", "cell")
      .attr("x", function (d) { return x(d.x); })
      .attr("width", x.rangeBand())
      .attr("height", x.rangeBand())
      .style("fill", function (d) { return c(d.z) })
      .on("click", click);
  }

  
  if ('max_dispersion' in data) {

      // Populate heatmap legend
      var legendHtml = '<div class="heatmap-legend">'
          + '<span>Range: 0 to ' + Math.round(data.max_dispersion*100)/100 + '</span>'
          + '<div class="temperature-scale">'
          + '<span class="gray-to-red"></span>'
          + '</div></div>';

  } else {

      // Populate heatmap legend
    
      var legend_def = svg.append("defs")
        .append("svg:linearGradient")
        .attr("id", "gradient")
        .attr("x1", "0%")
        .attr("y1", "100%")
        .attr("x2", "100%")
        .attr("y2", "100%")
        .attr("spreadMethod", "pad");

      legend_def.append("stop")
        .attr("offset", "0%")
        .attr("stop-color", "red")
        .attr("stop-opacity", 1);

      legend_def.append("stop")
        .attr("offset", "50%")
        .attr("stop-color", "white")
        .attr("stop-opacity", 1);

      legend_def.append("stop")
        .attr("offset", "100%")
        .attr("stop-color", "blue")
        .attr("stop-opacity", 1);
    
      label = svg.append("g")
      label.append("text")
        .text('Range: ' + Math.round(dis_min * 10) / 10 + 'Å to ' + Math.round(dis_max * 10) / 10 + 'Å')
        .attr("y", height+20)
        .attr("x", 20)
    
      label.append("rect")
        .attr("y", height+25)
        .attr("x", 20)
        .attr("width", 200)
        .attr("height", 10)
        .style("fill", "url(#gradient)");

  }
  
  function click(p) {

    $this = $(this);
    // Create the data for popover
    $(this).attr("title","" + gn_index[p.x] + " to "  + gn_index[p.y]);
    $(this).attr("data-content","Mean distance: " + p.z + "Å");
    $this.popover({
      'container': heatMapSelector,
      'placement': 'bottom',
      'animation': false,
      'html': true,
      'tabindex': '0'
    }).popover('show');
  }

  // be sure to remove popovers when user clicks elsewhere.
  $('html').on('mousedown', function(e) {
      if(!$(e.target).closest('.popover').length) {
          if ($(e.target).closest(heatMapSelector).length) {
            $('.popover').each(function(){
              $(this).remove();
            });
          }
      }
  });

  function mouseout() {
    d3.selectAll("text").classed("active", false);
  }

}
// HEATMAP REPRESENTATION

function simple_heatmap(data, location, element_id, legend_label) {
  // Set the dimensions and margins of the graph
  var margin = {top: 30, right: 30, bottom: 30, left: 30};

  // Extract rows and columns from the new data format
  var rows = Object.keys(data);
  var cols = Object.keys(data[rows[0]]);

  // Prepare chart data
  var chartData = [];
  var highest_value = 0;

  for (var row in data) {
    for (var col in data[row]) {
      chartData.push({ row: row, col: col, value: data[row][col] });
      if (data[row][col] > highest_value) {
        highest_value = data[row][col];
      }
    }
  }

  var width = (35 * cols.length);
  var height = (20 * rows.length);

  // Append the SVG object to the body of the page
  var svg_home = d3.select("#" + location)
    .append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + (margin.top * 5))
    .attr("transform", "translate(0,-" + margin.bottom + ")")
    .attr("id", element_id);

  var redscale = ['#ffcccc', '#ff0000']; // Changed from greyscale

  var legend = svg_home.append('defs')
    .append('linearGradient')
    .attr('id', 'grad_' + element_id)
    .attr('x1', '0%')
    .attr('x2', '100%')
    .attr('y1', '0%')
    .attr('y2', '0%');

  legend.selectAll('stop')
    .data(redscale) // Changed from greyscale
    .enter()
    .append('stop')
    .style('stop-color', function(d) { return d; })
    .attr('offset', function(d, i) {
      return 100 * (i / (redscale.length - 1)) + '%'; // Changed from greyscale
    });


  var legend_svg = svg_home.append("g")
    .attr("transform", "translate(0," + (height + 100) + ")");

  var color_svg = svg_home.append("g")
    .attr("transform", "translate(" + (margin.left * 1.5) + ",0)");

  var svg = svg_home.append("g")
    .attr("transform", "translate(" + (margin.left * 2) + ",0)");

  legend_svg.append("text")
    .attr('y', 10)
    .style("font", "14px sans-serif")
    .style("font-weight", "bold")
    .text(legend_label);

  legend_svg.append("text")
    .attr('x', 30)
    .attr('y', 30)
    .style("font", "14px sans-serif")
    .style("font-weight", "bold")
    .text('0');

  legend_svg.append('rect')
    .attr('x', 40)
    .attr('y', 15)
    .attr('width', (width * 0.9) + margin.left + margin.right)
    .attr('height', 20)
    .style('fill', 'url(#grad_' + element_id + ')');

  legend_svg.append("text")
    .attr('x', width + margin.left + (margin.right * 1.75))
    .attr('y', 30)
    .style("font", "14px sans-serif")
    .style("font-weight", "bold")
    .text(highest_value);

  // Using each to ensure the text element is fully created
  var legendLabel;
  legend_svg.selectAll("text").each(function() {
    legendLabel = this.getBBox().width * 1.05 + 0.5 * 10;
  });

  legend_svg.select("text")
    .attr("x", (width + margin.left + margin.right - legendLabel) / 2 + margin.left);

  // Build X scales and axis
  var x = d3.scale.ordinal()
    .rangeBands([0, width], 0.01)
    .domain(cols);

  svg.append("g")
    .attr("transform", "translate(0," + height + ")")
    .attr('id', 'Xaxis')
    .call(d3.svg.axis().scale(x).orient("bottom").tickSize(0))
    .selectAll("text")
    .style("text-anchor", "end")
    .attr("transform", "rotate(-45)")
    .attr("dx", "-0.8em")
    .attr("dy", "0.15em");

  // Build Y scales and axis
  var y = d3.scale.ordinal()
    .rangeBands([height, 0], 0.01)
    .domain(rows);

  svg.append("g")
    .attr('id', 'Yaxis')
    .call(d3.svg.axis().scale(y).orient("left").tickSize(0));


  // Build color scale
  var myColor = d3.scale.linear()
    .range(["#ffcccc", "#ff0000"]) // Changed from ["white", "black"]
    .domain([1, highest_value]);

  // Read the data
  svg.selectAll("rect")
    .data(chartData, function(d) { return d.row + ':' + d.col; })
    .enter()
    .append("rect")
    .attr("x", function(d) { return x(d.col); })
    .attr("y", function(d) { return y(d.row); })
    .attr("width", x.rangeBand())
    .attr("height", y.rangeBand())
    .style("fill", function(d) { return myColor(d.value); });

  svg.select('#Xaxis')
    .attr('text-anchor', 'start');

  d3.selectAll("#Yaxis>.tick>text")
    .each(function(d, i) {
      d3.select(this).style("font-size", "1.3em");
    });

  d3.selectAll("#Xaxis>.tick>text")
    .each(function(d, i) {
      d3.select(this).style("font-size", "1.3em");
    });

  var count = 1;
  var ticks = svg.select('#Xaxis').selectAll('.tick');
  ticks.each(function(d) {
    var value = (count & 1) ? "odd" : "even";
    var text = d3.select(this).select('text');
    var textSize = Math.floor(text.node().getBBox().width * 1.05 + 0.5 * 10);
    text.attr("transform", null);
    text.attr("x", "-" + textSize / 2 + "px");
    if (value === "even") {
      text.attr("y", "30");
    }
    count = count + 1;
  });
}

// LIST REPRESENTATION

function renderDataVisualization(data, svgSelector='svg') {
    const svg = d3.select(svgSelector);
    let yOffset = 30;

    function drawItems(items, xOffset) {
        items.forEach(item => {
            if (typeof item === 'object' && item.url && item.text) {
                // Draw linkable item
                const link = svg.append('a')
                                .attr('xlink:href', item.url)
                                .attr('target', '_blank');
                link.append('text')
                    .attr('x', xOffset)
                    .attr('y', yOffset)
                    .text(`- ${item.text}`);
            } else {
                // Draw regular item
                svg.append('text')
                   .attr('x', xOffset)
                   .attr('y', yOffset)
                   .text(`- ${item}`);
            }
            yOffset += 20; // Increase Y offset for next item
        });
    }

    function processNode(node, xOffset, depth = 0) {
        Object.entries(node).forEach(([key, value]) => {
            // Append text for the node
            svg.append('text')
               .attr('x', xOffset + 20)
               .attr('y', yOffset)
               .text(key);

            // Append visual markers based on depth
            if (depth === 0) {
                // Red circle for root nodes
                svg.append('circle')
                   .attr('cx', xOffset + 10)
                   .attr('cy', yOffset - 4)
                   .attr('r', 4)
                   .attr('fill', 'red');
            } else if (depth === 1) {
                // Blue square for first nested level
                svg.append('circle')
                   .attr('cx', xOffset + 10)
                   .attr('cy', yOffset - 4)
                   .attr('r', 4)
                   .attr('fill', 'blue');
            } else if (depth === 2) {
                // Black square for second nested level
                svg.append('circle')
                   .attr('cx', xOffset + 10)
                   .attr('cy', yOffset - 4)
                   .attr('r', 4)
                   .attr('fill', 'black');
            }

            // Increment Y offset after title
            yOffset += 30;

            if (Array.isArray(value)) {
                // Leaf node, draw items
                drawItems(value, xOffset + 40);
            } else {
                // Non-leaf node, recursive process
                processNode(value, xOffset + 20, depth + 1);
            }
        });
    }
    // Start processing from the root
    processNode(data, 10);
}

function renderClusterPlot(data, selector, method) {

    const width = 800;
    const height = 600;
    const margin = { top: 50, right: 50, bottom: 50, left: 50 };

    // Clear the previous plot
    d3.select(selector).html("");

    const svg = d3.select(selector)
        .append("svg")
        .attr("width", width)
        .attr("height", height)
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    const x = d3.scale.linear()
        .domain(d3.extent(data, d => +d.x))
        .range([0, width - margin.left - margin.right]);

    const y = d3.scale.linear()
        .domain(d3.extent(data, d => +d.y))
        .range([height - margin.top - margin.bottom, 0]);

    svg.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + (height - margin.top - margin.bottom) + ")")
        .call(d3.svg.axis().scale(x).orient("bottom"));

    svg.append("g")
        .attr("class", "y axis")
        .call(d3.svg.axis().scale(y).orient("left"));

    svg.selectAll(".dot")
        .data(data)
        .enter().append("circle")
        .attr("class", "dot")
        .attr("cx", d => x(+d.x))
        .attr("cy", d => y(+d.y))
        .attr("r", 5)
        .style("fill", "steelblue");

    svg.selectAll(".text")
        .data(data)
        .enter().append("text")
        .attr("x", d => x(+d.x))
        .attr("y", d => y(+d.y))
        .attr("dy", -10)
        .text(d => d.label)
        .style("font-size", "10px")
        .style("text-anchor", "middle");

    // Add title based on the method
    svg.append("text")
        .attr("x", (width - margin.left - margin.right) / 2)
        .attr("y", -20)
        .attr("text-anchor", "middle")
        .style("font-size", "20px")
        .text("Cluster Plot (" + method.toUpperCase() + ")");
}

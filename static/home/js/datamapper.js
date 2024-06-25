// HEATMAP REPRESENTATION

function simple_heatmap(data, location, element_id, legend_label) {
  // Set the dimensions and margins of the graph
  var margin = { top: 30, right: 100, bottom: 30, left: 60 }; // Increased left margin to accommodate row labels

  // Extract rows and columns from the new data format
  var rows = Object.keys(data);
  var cols = Object.keys(data[rows[0]]);

  // Calculate the width required for row labels
  var rowLabelWidth = d3.max(rows, function(d) {
    return d.length * 25; // Adjust the multiplier as needed
  });

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

  var width = (45 * cols.length) + rowLabelWidth; // Adjust the multiplier as needed
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
    .attr("transform", "translate(30," + (height + 100) + ")");

  var color_svg = svg_home.append("g")
    .attr("transform", "translate(" + (margin.left * 1.5) + ",0)");

  var svg = svg_home.append("g")
    .attr("transform", "translate(" + (margin.left * 2) + ",0)");

  legend_svg.append("text")
    .attr('y', 10)
    .attr('x', -50)
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
    .attr('width', (width * 0.7) + margin.left + margin.right)
    .attr('height', 20)
    .style('fill', 'url(#grad_' + element_id + ')');

  legend_svg.append("text")
    .attr('x', width + margin.left)
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
    .attr("x", (width + margin.left + margin.right - legendLabel) / 2 + margin.left - 120);

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
    .attr("dy", "0.15em")
    .attr("class", "column-label"); // Add a class for styling

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
    // text.attr("x", "-" + textSize / 2 + "px");
    text.attr("x","42")
    text.attr("y","15")
    // if (value === "even") {
    //   text.attr("y", "0");
    // }
    count = count + 1;
  });
}

// LIST REPRESENTATION

function renderDataVisualization(data, location) {
    // Define width, height, and margins for the SVG
    const width = 800; // Adjust as needed
    const height = 600; // Adjust as needed
    const margin = { top: 20, right: 20, bottom: 20, left: 20 }; // Adjust as needed

    // Append SVG to the specified location
    const svg_home = d3.select("#" + location)
        .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .attr("id", "visualization");

    // Initialize yOffset
    let yOffset = margin.top;

    // Define drawItems function
    function drawItems(items, xOffset) {
        items.forEach(item => {
            const key = Object.keys(item)[0];
            const values = item[key];

            // Draw item key
            const textElement = svg_home.append('text')
                .attr('x', xOffset)
                .attr('y', yOffset)
                .text(`- ${key}`);

            // Get the width of the text element to adjust the xOffset for values
            const keyWidth = textElement.node().getBBox().width;
            let valueXOffset = xOffset + keyWidth + 10;

            values.forEach(value => {
                if (typeof value === 'string') {
                    // Try to convert string to number
                    const numericValue = parseFloat(value);
                    if (!isNaN(numericValue)) {
                        value = numericValue;
                    }
                }

                if (typeof value === 'number') {
                    // Numeric value
                    svg_home.append('text')
                        .attr('x', valueXOffset * 2.75)
                        .attr('y', yOffset)
                        .style('font-weight', 'bold')
                        .text(`${value}`);
                    valueXOffset += 5;
                } else if (value === 'yes') {
                    // Boolean yes
                    svg_home.append('text')
                        .attr('x', valueXOffset * 2.75)
                        .attr('y', yOffset)
                        .style('fill', 'green')
                        .text(' ✓');
                    valueXOffset += 5;
                } else if (value === 'no') {
                    // Boolean no
                    svg_home.append('text')
                        .attr('x', valueXOffset * 2.75)
                        .attr('y', yOffset)
                        .style('fill', 'red')
                        .text(' ✗');
                    valueXOffset += 5;
                }
            });

            yOffset += 20; // Increase Y offset for next item
        });
    }

    // Define processNode function
    function processNode(node, xOffset, depth = 0) {
        Object.entries(node).forEach(([key, value]) => {
            // Append text for the node
            svg_home.append('text')
                .attr('x', xOffset + 20)
                .attr('y', yOffset)
                .text(key);

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
    processNode(data, margin.left);

    // Return the SVG element
    return svg_home;
}

// RENDER THE CLUSTER

function renderClusterPlot(data, selector, method) {
    const width = 800;
    const height = 600;
    const margin = { top: 80, right: 150, bottom: 50, left: 50 };

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

    const color = d3.scale.category10();

    let selectedCluster = null;

    const circles = svg.selectAll(".dot")
        .data(data)
        .enter().append("circle")
        .attr("class", "dot")
        .attr("cx", d => x(+d.x))
        .attr("cy", d => y(+d.y))
        .attr("r", 5)
        .attr("cluster", d => d.cluster)
        .style("fill", d => color(d.cluster))
        .attr("id", d => "circle" + d.label.replace(/\[|\]|\(|\)|\s|\,|\'/g,""))
        .on("click", handleClick);

    const labels = svg.selectAll(".label")
        .data(data)
        .enter().append("text")
        .attr("class", "label")
        .attr("x", d => x(+d.x))
        .attr("y", d => y(+d.y))
        .attr("dy", -10)
        .attr("cluster", d => d.cluster)
        .text(d => d.label)
        .style("font-size", "11px")
        .style("text-anchor", "middle")
        .on("click", function(d, i) {
            d3.event.stopPropagation();
            handleClick.call(circles[0][i], d);
        });

    // Add leader lines from circles to labels
    const lines = svg.selectAll(".line")
        .data(data)
        .enter().append("line")
        .attr("class", "line")
        .style("stroke", "black")
        .style("stroke-width", 1);

    // Add title based on the method with more spacing
    svg.append("text")
        .attr("x", (width - margin.left - margin.right) / 2)
        .attr("y", -50)
        .attr("text-anchor", "middle")
        .style("font-size", "20px")
        .text("Cluster Plot (" + method.toUpperCase() + ")");

    // Add force simulation to spread out circles
    const simulation = d3v4.forceSimulation(data)
        .force("x", d3v4.forceX(d => x(d.x)).strength(1))
        .force("y", d3v4.forceY(d => y(d.y)).strength(1))
        .force("collide", d3.forceCollide(10)) // Set to twice the circle radius for spacing
        .on("tick", ticked);

    function ticked() {
        circles
            .attr("cx", d => d.x = Math.max(5, Math.min(width - 5, d.x))) // Prevent circles from going outside SVG
            .attr("cy", d => d.y = Math.max(5, Math.min(height - 5, d.y))); // Prevent circles from going outside SVG

        labels
            .attr("x", d => d.x)
            .attr("y", d => d.y - 10);

        lines
            .attr("x1", d => x(d.x))
            .attr("y1", d => y(d.y))
            .attr("x2", d => d.x)
            .attr("y2", d => d.y - 10);
    }

    // Function to handle click event
    function handleClick(d) {
        d3.event.stopPropagation();
        const cluster = d.cluster;

        if (selectedCluster === cluster) {
            selectedCluster = null;
            resetOpacity();
        } else {
            resetOpacity();
            selectedCluster = cluster;
            circles.style("opacity", 0.2);
            labels.style("opacity", 0.2);
            lines.style("opacity", 0.2);

            circles.filter(c => c.cluster === cluster)
                .style("opacity", 1)
                .style("stroke", "black")
                .style("stroke-width", 2);

            labels.filter(c => c.cluster === cluster)
                .style("opacity", 1);

            lines.filter(c => c.cluster === cluster)
                .style("opacity", 1);
        }
    }

    // Function to reset the opacity when clicking outside circles
    function resetOpacity() {
        circles.style("opacity", 1).style("stroke", "none").style("stroke-width", 0);
        labels.style("opacity", 1).style("font-weight", "normal");
        lines.style("opacity", 1);
    }

    // Add an overlay to capture clicks outside circles
    svg.append("rect")
        .attr("width", width - margin.left - margin.right)
        .attr("height", height - margin.top - margin.bottom)
        .style("fill", "none")
        .style("pointer-events", "all")
        .on("click", resetOpacity);

    // Create a legend based on clusters
    const uniqueClusters = [...new Set(data.map(d => d.cluster))].sort((a, b) => a - b);

    const legend = svg.selectAll(".legend")
        .data(uniqueClusters)
        .enter().append("g")
        .attr("class", "legend")
        .attr("transform", (d, i) => `translate(${width - margin.right + 10},${i * 20})`)
        .on("click", function(d) {
            selectedCluster = d;
            circles.style("opacity", 0.2);

            circles.filter(c => c.cluster === d)
                .style("opacity", 1)
                .style("stroke", "black")
                .style("stroke-width", 2);

            d3.event.stopPropagation();
        });

    legend.append("circle")
        .attr("cx", 10)
        .attr("cy", 9)
        .attr("r", 5)
        .style("fill", d => color(d));

    legend.append("text")
        .attr("x", 20)
        .attr("y", 14)
        .text(d => `Cluster ${d}`)
        .style("font-size", "12px");

    // Adjust label positions to prevent overlap
    adjustLabelPositions(labels);
}

// Function to adjust label positions to prevent overlap
function adjustLabelPositions(labels) {
    const labelPadding = 2; // Padding between labels to avoid overlap
    let iterations = 10; // Maximum iterations to adjust labels

    while (iterations-- > 0) {
        let overlaps = false;

        labels.each(function() {
            const label = d3.select(this);
            const bbox = label.node().getBBox();

            labels.each(function() {
                const otherLabel = d3.select(this);
                if (label.node() !== otherLabel.node()) {
                    const otherBBox = otherLabel.node().getBBox();
                    if (intersect(bbox, otherBBox)) {
                        overlaps = true;
                        const dx = (bbox.x + bbox.width / 2) - (otherBBox.x + otherBBox.width / 2);
                        const dy = (bbox.y + bbox.height / 2) - (otherBBox.y + otherBBox.height / 2);

                        const offset = labelPadding / Math.sqrt(dx * dx + dy * dy);
                        const moveX = dx * offset;
                        const moveY = dy * offset;

                        label.attr("x", parseFloat(label.attr("x")) + moveX)
                            .attr("y", parseFloat(label.attr("y")) + moveY);
                        otherLabel.attr("x", parseFloat(otherLabel.attr("x")) - moveX)
                            .attr("y", parseFloat(otherLabel.attr("y")) - moveY);

                        bbox.x += moveX;
                        bbox.y += moveY;
                    }
                }
            });
        });

        if (!overlaps) break; // Stop if no overlaps found
    }
}

// Function to check if two bounding boxes intersect
function intersect(bbox1, bbox2) {
    return !(bbox2.x > bbox1.x + bbox1.width ||
             bbox2.x + bbox2.width < bbox1.x ||
             bbox2.y > bbox1.y + bbox1.height ||
             bbox2.y + bbox2.height < bbox1.y);
}


// function renderClusterPlot(data, selector, method) {
//     const width = 800;
//     const height = 600;
//     const margin = { top: 80, right: 50, bottom: 50, left: 50 };
//
//     // Clear the previous plot
//     d3.select(selector).html("");
//
//     const svg = d3.select(selector)
//         .append("svg")
//         .attr("width", width)
//         .attr("height", height)
//         .append("g")
//         .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
//
//     const x = d3.scale.linear()
//         .domain(d3.extent(data, d => +d.x))
//         .range([0, width - margin.left - margin.right]);
//
//     const y = d3.scale.linear()
//         .domain(d3.extent(data, d => +d.y))
//         .range([height - margin.top - margin.bottom, 0]);
//
//     const color = d3.scale.category10();
//
//     const circles = svg.selectAll(".dot")
//         .data(data)
//         .enter().append("circle")
//         .attr("class", "dot")
//         .attr("cx", d => x(+d.x))
//         .attr("cy", d => y(+d.y))
//         .attr("r", 5)
//         .attr("cluster", d => d.cluster)
//         .style("fill", d => color(d.cluster))
//         .attr("id", d => "circle" + d.label.replace(/\[|\]|\(|\)|\s|\,|\'/g,""))
//         .on("click", handleClick);
//
//     const labels = svg.selectAll(".text")
//         .data(data)
//         .enter().append("text")
//         .attr("x", d => x(+d.x))
//         .attr("y", d => y(+d.y))
//         .attr("dy", -10)
//         .attr("cluster", d => d.cluster)
//         .text(d => d.label)
//         .style("font-size", "11px")
//         .style("text-anchor", "middle");
//         // .on("click", function(d, i) {
//         //     d3.event.stopPropagation();
//         //     handleClick.call(circles[0][i], d);
//         // });
//
//     // Add title based on the method with more spacing
//     svg.append("text")
//         .attr("x", (width - margin.left - margin.right) / 2)
//         .attr("y", -50)
//         .attr("text-anchor", "middle")
//         .style("font-size", "20px")
//         .text("Cluster Plot (" + method.toUpperCase() + ")");
//
//     // Function to handle click event
//     function handleClick(d) {
//         const cluster = d.cluster;
//         circles.style("opacity", 0.2);
//         labels.style("opacity", 0.2);
//
//         circles.filter(c => c.cluster === cluster)
//             .style("opacity", 1)
//             .style("stroke", "black")
//             .style("stroke-width", 2);
//
//         labels.filter(c => c.cluster === cluster)
//             .style("opacity", 1)
//             .style("font-weight", "bold");
//     }
//
//     // Function to reset the opacity when clicking outside circles
//     function resetOpacity() {
//         circles.style("opacity", 1).style("stroke", "none").style("stroke-width", 0);
//         labels.style("opacity", 1).style("font-weight", "normal");
//     }
//
//     // Add an overlay to capture clicks outside circles
//     // svg.append("rect")
//     //     .attr("width", width - margin.left - margin.right)
//     //     .attr("height", height - margin.top - margin.bottom)
//     //     .style("fill", "none")
//     //     .style("pointer-events", "all")
//     //     .on("click", resetOpacity);
// }

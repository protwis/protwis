// #################
// ###   TREE    ###
// #################

function draw_tree(data, options,circle_size) {

    // Remove existing SVG if present
    d3.select('#' + options.anchor + "_svg").remove();

    var branches = {};
    var branch_offset = 0;
    var thickness = options.depth + 1;
    for (var key in options.branch_length) {
        if (key == options.depth) { continue };
        if (options.label_free.includes(parseInt(key))) {
            branch_offset = branch_offset + 10;
        } else {
            if (options.branch_trunc != 0) {
                branch_offset = branch_offset + 2 * options.branch_trunc + 10;
            } else {
                branch_offset = branch_offset + string_pixlen(options.branch_length[key], key);
            }
        }
        branches[key] = branch_offset;
    }
    branches[options.depth] = branch_offset + options.leaf_offset;

    var diameter = 2 * branches[options.depth] + 140;

    var tree = d3.layout.tree()
        .size([360, diameter / 2])
        .separation(function (a, b) { return (a.parent == b.parent ? 1 : 2) / a.depth; });

    var diagonal = d3.svg.diagonal.radial()
        .projection(function (d) { return [d.y, d.x / 180 * Math.PI]; });

    var svg = d3.select('#' + options.anchor).append("svg")
        .attr("width", diameter)
        .attr("height", diameter)
        .attr("id", options.anchor + "_svg")
        .attr("xmlns", "http://www.w3.org/2000/svg");

    var svg_g = svg.append("g")
        .attr("transform", "translate(" + diameter / 2 + "," + diameter / 2 + ")");

    var nodes = tree.nodes(data);

    nodes.forEach(function (d) {
        if (d.depth == 0) {
            d.y = 0
        } else {
            d.y = branches[d.depth]
        }
    });

    var links = tree.links(nodes);

    var link = svg_g.append("g")
        .attr("class", "links")
        .selectAll("path")
        .data(links)
        .enter().append("path")
        .each(function (d) { d.target.linkNode = this; })
        .attr("d", diagonal) //function (d) { return step(d.source.x, d.source.y, d.target.x, d.target.y) })
        .style("stroke", function (d) { return d.target.color; })
        .style("stroke-width", function (d) { if (d.target.depth > 0) { return thickness - d.target.depth; } else { return 0; } })
        .style("fill-opacity", 0)
        .style("opacity", function (d) {
            if ((d.target.interactions > 0 && d.target.mutations_an > 0) || 1 == 1) { return 0.8 } //|| 1==1
            else if (d.target.interactions > 0) { return 0.5 }
            else if (d.target.mutations_an > 0) { return 0.5 }
            else { return 0.1 };
        });

    var node = svg_g.selectAll(".node")
        .data(nodes)
        .enter().append("g")
        .attr("class", "node")
        .attr("transform", function (d) { if (d.name == '') { return "rotate(" + (d.x) + ")translate(" + d.y + ")"; } else { return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")"; } })
//TODO: add a check to remove circles when nothing is passed (?)
    node.filter(function (d) { return (d.depth == options.depth) })
        .filter(function (d) { return (d.value !== 3000) })
        .append("circle")
        .attr("r", function (d) { if (d.name == '') { return "0" } else { return circle_size } })
        .style("stroke", "black")
        .style("stroke-width", ".3px")
        .style("fill", function (d) {
            if (d.color && d.depth < options.depth) { return d.color }
            else if (d.value === 1) {
                return "FireBrick";
            }
            else if (d.value === 10) {
                return "LightGray";
            }
            else if (d.value === 20) {
                return "DarkGray";
            }
            else if (d.value === 30) {
                return "Gray";
            }
            else if (d.value === 40) {
                return "Black";
            }
            else if (d.value === 100) {
                return 'LightGray';
            }
            else if (d.value === 500) {
                return 'DarkGray';
            }
            else if (d.value === 1000) {
                return 'Gray';
            }
            else if (d.value === 2000) {
                return 'Black';
            }
            else { return "White" };
        })
        .style("opacity", .99);

    node.filter(function (d) { return (d.depth == options.depth) })
        .attr("id", function (d) { if (d.name == '') { return "innerNode" } else { return 'X' + d.name.toUpperCase() } });

    node.append("text")
        .attr("dy", ".31em")
        .attr("name", function (d) { if (d.name == '') { return "branch" } else { return d.name } })
        .attr("text-anchor", function (d) {
            if (d.depth == options.depth) {
                return d.x < 180 ? "start" : "end";
            } else {
                return d.x < 180 ? "end" : "start";
            }
        })
        .attr("transform", function (d) {
            var labelOffset = parseFloat(circle_size) + 4;  // Adjust the offset as needed
            if (d.depth == options.depth) {
                return d.x < 180 ? `translate(${labelOffset})` : `rotate(180)translate(-${labelOffset})`;
            } else {
                return d.x < 180 ? "translate(-12)" : "rotate(180)translate(12)";
            }
        })
        .text(function (d) {
            if (d.depth == options.depth) {
                return d.name.toUpperCase();
            } else if (options.label_free.includes(d.depth)) {
                return "";
            } else if (d.depth > 0) {
                return d.name;
            } else {
                return "";
            }
        })
        .call(wrap, options.branch_trunc)
        .style("font-size", function (d) { 
            // Use the custom font size from options
            if (d.depth == 1) { return options.fontSize.class; } 
            else if (d.depth == 2) { return options.fontSize.ligandtype; } 
            else if (d.depth == 3) { return options.fontSize.receptorfamily; } 
            else { return options.fontSize.receptor; } 
        })
        .style("font-family", "Palatino")
        .style("fill", function (d) {
            if (d.color) { return "#111"; }
            else { return "#222"; };
        }).call(getBB);

    node.filter(function (d) { return (d.depth != options.depth) }).insert("rect", "text")
        .attr("x", function (d) { return d.x < 180 ? d.bbox.x - 12 : d.bbox.x - d.bbox.width - 12; })
        .attr("y", function (d) { return d.bbox.y; })
        .attr("width", function (d) { return d.bbox.width; })
        .attr("height", function (d) { return d.bbox.height; })
        .style("fill", "#FFF");

    function step(startAngle, startRadius, endAngle, endRadius) {
        var c0 = Math.cos(startAngle = (startAngle - 90) / 180 * Math.PI),
            s0 = Math.sin(startAngle),
            c1 = Math.cos(endAngle = (endAngle - 90) / 180 * Math.PI),
            s1 = Math.sin(endAngle);
        return "M" + startRadius * c0 + "," + startRadius * s0
            + (endAngle === startAngle ? "" : "A" + startRadius + "," + startRadius + " 0 0 " + (endAngle > startAngle ? 1 : 0) + " " + startRadius * c1 + "," + startRadius * s1)
            + "L" + endRadius * c1 + "," + endRadius * s1;
    }

    function string_pixlen(text, depth) {
        var canvas = document.createElement('canvas');
        var ctx = canvas.getContext("2d");
        // Use the custom font size from options
        if (depth == 1) {
            ctx.font = options.fontSize.class + " Palatino";
        } else if (depth == 2) {
            ctx.font = options.fontSize.ligandtype + " Palatino";
        } else if (depth == 3) {
            ctx.font = options.fontSize.receptorfamily + " Palatino";
        } else {
            ctx.font = options.fontSize.receptor + " Palatino";
        }
        return parseInt(ctx.measureText(text).width) + 40;
    }

    function getBB(selection) {
        selection.each(function (d) { d.bbox = this.getBBox(); })
    }

    function wrap(text, width) {
        if (width == 0) {
            return;
        }
        text.each(function () {
            var text = d3.select(this),
                words = text.text().split(/\s+/).reverse(),
                word,
                line = [],
                lineNumber = 0,
                lineHeight = 1.1, // ems
                y = text.attr("y"),
                dy = parseFloat(text.attr("dy")),
                tspan = text.text(null).append("tspan").attr("x", 0).attr("y", y).attr("dy", dy + "em");
            while (word = words.pop()) {
                line.push(word);
                tspan.text(line.join(" "));
                if (tspan.node().getComputedTextLength() > width) {
                    line.pop();
                    tspan.text(line.join(" "));
                    line = [word];
                    tspan = text.append("tspan").attr("x", 0).attr("y", y).attr("dy", ++lineNumber * lineHeight + dy + "em").text(word);
                }
            }
        });
    }
    // === Centering Logic ===
    var scaleFactor = 0.8;  // Adjust this as needed
    var svgElement = d3.select('#' + options.anchor + ' svg');
    var svgWidth = +svgElement.attr('width');
    var svgHeight = +svgElement.attr('height');

    // Calculate the center point
    var cx = svgWidth / 2;
    var cy = svgHeight / 2;

    // Adjust the translation to keep the tree centered after scaling
    var translateX = cx;
    var translateY = cy;

    // Apply the transform to the 'g' element to scale and center it
    svgElement.select('g')
        .attr('transform', `translate(${translateX},${translateY}) scale(${scaleFactor},${scaleFactor})`);
}


function changeLeavesLabels(location, value, dict){
    // Initialize leaf node length
    maxLeafNodeLenght = 0;
    // Find longest label
    gNodes = d3.select('#'+location).selectAll('g');
    gNodes.each(function(d) {
      if (d3.select(this).attr("id") !== null) {
        name = d3.select(this).attr("id").substring(1);
        labelName = dict[name][0];
        // replaces labels derived from view
        labelName = labelName.replace("-adrenoceptor", '');
        labelName = labelName.replace(" receptor-", '-');
        labelName = labelName.replace("<sub>", '</tspan><tspan baseline-shift = "sub">');
        labelName = labelName.replace("</sub>", '</tspan><tspan>');
        labelName = labelName.replace("<i>", '</tspan><tspan font-style = "italic">');
        labelName = labelName.replace("</i>", '</tspan><tspan>');
        labelName = labelName.replace("Long-wave-sensitive",'LWS');
        labelName = labelName.replace("Medium-wave-sensitive",'MWS');
        labelName = labelName.replace("Short-wave-sensitive",'SWS');
        labelName = labelName.replace("Olfactory", 'OLF');
        labelName = labelName.replace("Calcitonin -like", 'CLR');
        node = d3.select('#'+location).select('#X'+name);
        if (node.size() !== 0){
          if (value === "IUPHAR"){
            node.selectAll("text")[0].forEach(
              function(node_label){
                node_label.innerHTML = labelName;
                labelSize = node_label.getBBox().width*1.05 + 0.5 * 10
                if (labelSize > maxLeafNodeLenght){
                  // change initialization label length, needed for outer circles
                  maxLeafNodeLenght = labelSize
                }
              });
          } else if (value === "UniProt"){
            node.selectAll("text")[0].forEach(
              function(node_label){
                node_label.innerHTML = name;
                labelSize = node_label.getBBox().width*1.05 + 0.5 * 10
                if (labelSize > maxLeafNodeLenght){
                  maxLeafNodeLenght = labelSize
                }
              });
          }
        }
      }
    });
  }

  function DrawCircles(location, data, starter, dict, clean = true, gradient = true, circle_styling_dict,circle_spacer,circle_size) {
    var svg = d3.select('#' + location);
    var node = svg.selectAll(".node");

    if (clean === true) {
        node.selectAll("circle.outerCircle").remove();
    }

    var spacer = circle_spacer;

    // Initialize a dictionary to store min and max values for each unit
    var minMaxValues = {};

    // First pass: Determine min and max values for each unit
    for (var x in data) {
        var keys = Object.keys(data[x]);
        for (var unit of keys) {
            if (!minMaxValues[unit]) {
                minMaxValues[unit] = { min: Infinity, max: -Infinity };
            }
            var value = data[x][unit];
            if (value < minMaxValues[unit].min) {
                minMaxValues[unit].min = value;
            }
            if (value > minMaxValues[unit].max) {
                minMaxValues[unit].max = value;
            }
        }
    }

    // Second pass: Draw the circles
    for (var x in data) {
        var keys = Object.keys(data[x]);

        for (var unit of keys) {
            if (dict[unit]) {
                var value = data[x][unit];
                var minValue = minMaxValues[unit].min;
                var maxValue = minMaxValues[unit].max;

                // Determine the styling

                var styling = circle_styling_dict[unit] || "Two";

                // Create color scale based on min and max values
                var colorScale;
                if (styling === "Three") {
                    // Three-color gradient with white in the middle
                    colorScale = d3.scale.linear()
                        .domain([minValue, (minValue + maxValue) / 2, maxValue])
                        .range([dict[unit][0], "#FFFFFF", dict[unit][1]]);
                } else {
                    // Two-color gradient
                    colorScale = d3.scale.linear()
                        .domain([minValue, maxValue])
                        .range(dict[unit]);
                }

                // Calculate color using the color scale
                var color = gradient ? colorScale(value) : dict[unit][0]; // Use the first color if no gradient

                var multiply = 1 + Object.keys(dict).indexOf(unit);
                var leaf = svg.selectAll('g[id=X' + x + ']');

                leaf.append("circle")
                    .attr("r", circle_size)
                    .attr("class", "outerCircle") // Add class to distinguish outer circles
                    .style("stroke", "black") // Use the first color for the stroke
                    .style("stroke-width", 0.8)
                    .style("fill", color)
                    .attr("transform", "translate(" + (Math.ceil(starter) + multiply * spacer) + ",0)");
            }
        }
    }
}

function createLegendBars(location, data, conversion, circle_styling_dict) {
    var svg = d3.select('#' + location + ' svg');

    // Clear existing content
    svg.selectAll("*").remove();

    var margin = { top: 20, right: 20, bottom: 30, left: 60 };
    var width = +svg.attr("width") - margin.left - margin.right;
    var height = +svg.attr("height") - margin.top - margin.bottom;

    // Flatten data to get categories and their max values
    var categoryMax = {};

    Object.values(data).forEach(item => {
        Object.entries(item).forEach(([category, value]) => {
            if (category in conversion) {
                if (!categoryMax[category]) {
                    categoryMax[category] = { min: value, max: value };
                } else {
                    if (value > categoryMax[category].max) {
                        categoryMax[category].max = value;
                    }
                    if (value < categoryMax[category].min) {
                        categoryMax[category].min = value;
                    }
                }
            }
        });
    });

    var existingCategories = Object.keys(categoryMax).filter(cat => categoryMax[cat].max > 0);

    // Create a color scale function for each category
    var colorScales = {};
    existingCategories.forEach(category => {
        var colors = conversion[category];
        var gradientId = `gradient-${category}`;

        var gradient = svg.append("defs")
            .append("linearGradient")
            .attr("id", gradientId)
            .attr("x1", "0%")
            .attr("x2", "100%")
            .attr("y1", "0%")
            .attr("y2", "0%");

        var styling = circle_styling_dict[category] || "Two";
        if (styling === "Three") {
            gradient.append("stop")
                .attr("offset", "0%")
                .attr("stop-color", colors[0]);
            gradient.append("stop")
                .attr("offset", "50%")
                .attr("stop-color", "#FFFFFF");
            gradient.append("stop")
                .attr("offset", "100%")
                .attr("stop-color", colors[1]);
        } else {
            gradient.append("stop")
                .attr("offset", "0%")
                .attr("stop-color", colors[0]);
            gradient.append("stop")
                .attr("offset", "100%")
                .attr("stop-color", colors[1]);
        }

        colorScales[category] = gradientId;
    });

    // var barWidth = width / existingCategories.length - 20; // Adjust the width of each bar
    var barWidth = 120;
    var barHeight = 15;
    var spacing = 50; // Horizontal spacing between bars
    // console.log(barWidth);
    // Adjust SVG width if needed
    // var totalWidth = margin.left + (existingCategories.length * (barWidth + spacing)) + margin.right;
    var totalWidth = 1200;
    svg.attr("width", totalWidth);

    existingCategories.forEach((category, index) => {
        var gradientId = colorScales[category];
        var minValue = categoryMax[category].min;
        var maxValue = categoryMax[category].max;
        var midValue = (minValue + maxValue) / 2;

        var xPosition = margin.left + index * (barWidth + spacing); // Calculate the x position

        // Add a rectangle to represent the gradient
        svg.append("rect")
            .attr("x", xPosition)
            .attr("y", margin.top)
            .attr("width", barWidth)
            .attr("height", barHeight)
            .style("fill", `url(#${gradientId})`)
            .style("stroke", "black")
            .style("stroke-width", "1px");

        // Add a label for the gradient bar
        svg.append("text")
            .attr("x", xPosition + barWidth / 2)
            .attr("y", margin.top - 10)
            .attr("text-anchor", "middle")
            .text(category);

        // Add min and max value labels
        svg.append("text")
            .attr("x", xPosition)
            .attr("y", margin.top + barHeight + 15)
            .attr("text-anchor", "start")
            .text(`${minValue}`);

        svg.append("text")
            .attr("x", xPosition + barWidth)
            .attr("y", margin.top + barHeight + 15)
            .attr("text-anchor", "end")
            .text(`${maxValue}`);

        // Add mid value label
        svg.append("text")
            .attr("x", xPosition + barWidth / 2)
            .attr("y", margin.top + barHeight + 15)
            .attr("text-anchor", "middle")
            .text(`${midValue}`);
    });
}


// #################
// ###  CLUSTER  ###
// #################

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

// #################
// ###   LIST    ###
// #################

// ## Initialization of global data for list ##
function initializeDataStyling(list_data_wow, data_types_list) {

    // Function to check if any of the specified keys are present in a given object
    function checkKeysPresent(obj, keys) {
        return keys.some(k => k in obj);
    }

    // Define the keys to check for each column
    var col1Keys = ["Value1", "Value2"];
    var col2Keys = ["Value3", "Value4"];
    var col3Keys = ["Value5", "Value6"];
    var col4Keys = ["Value7", "Value8"];


    var Data_styling = {
        Col1: {Data: "No",Datatype: 'None',Data_min: Infinity,Data_max: -Infinity, Data_coloring: 'All', Data_color1: '#000000', Data_color2: '#000000', data_color_complexity: 'One'},
        Col2: {Data: "No",Datatype: 'None',Data_min: Infinity,Data_max: -Infinity, Data_coloring: 'All', Data_color1: '#000000', Data_color2: '#000000', data_color_complexity: 'One'},
        Col3: {Data: "No",Datatype: 'None',Data_min: Infinity,Data_max: -Infinity, Data_coloring: 'All', Data_color1: '#000000', Data_color2: '#000000', data_color_complexity: 'One'},
        Col4: {Data: "No",Datatype: 'None',Data_min: Infinity,Data_max: -Infinity, Data_coloring: 'All', Data_color1: '#000000', Data_color2: '#000000', data_color_complexity: 'One'}
        };

    // Iterate through initialData and update Data_styling accordingly
    Object.keys(list_data_wow).forEach(function(key) {
        var currentObj = list_data_wow[key];

        // Check for Col1 if not already set to "Yes"
        if (Data_styling.Col1.Data !== "Yes" && checkKeysPresent(currentObj, col1Keys)) {
            Data_styling.Col1.Data = "Yes";
        }

        // Check for Col2 if not already set to "Yes"
        if (Data_styling.Col2.Data !== "Yes" && checkKeysPresent(currentObj, col2Keys)) {
            Data_styling.Col2.Data = "Yes";
        }

        // Check for Col3 if not already set to "Yes"
        if (Data_styling.Col3.Data !== "Yes" && checkKeysPresent(currentObj, col3Keys)) {
            Data_styling.Col3.Data = "Yes";
        }

        // Check for Col4 if not already set to "Yes"
        if (Data_styling.Col4.Data !== "Yes" && checkKeysPresent(currentObj, col4Keys)) {
            Data_styling.Col4.Data = "Yes";
        }

        // If both Col1 and Col2 are "Yes", no need to continue iterating
        if (Data_styling.Col1.Data === "Yes" && Data_styling.Col2.Data === "Yes" && Data_styling.Col3.Data === "Yes" && Data_styling.Col4.Data === "Yes") {
            return; // Exit forEach loop early
        }
    });


    // Check for Col1 has data
    // Update Datatype for each column based on data_types
    ['Col1', 'Col2', 'Col3', 'Col4'].forEach(col => {
        if (Data_styling[col].Data === "Yes") {
            if (data_types_list[col] === 'Discrete') {
                Data_styling[col].Datatype = 'Discrete';
            } else if (data_types_list[col] === 'Continuous') {
                Data_styling[col].Datatype = 'Continuous';
            }
        }
    });


    // Run through data and set min and max values

    function updateDataMinMax(column, key, currentObj) {
        if (Data_styling[column].Data === "Yes" && Data_styling[column].Datatype === 'Continuous' && currentObj.hasOwnProperty(key)) {
            if (currentObj[key] > Data_styling[column].Data_max) {
                Data_styling[column].Data_max = currentObj[key];
            }
            if (currentObj[key] < Data_styling[column].Data_min) {
                Data_styling[column].Data_min = currentObj[key];
            }
        }
    }

    Object.keys(list_data_wow).forEach(function(key) {
        var currentObj = list_data_wow[key];

        updateDataMinMax('Col1', 'Value2', currentObj);
        updateDataMinMax('Col2', 'Value4', currentObj);
        updateDataMinMax('Col3', 'Value6', currentObj);
        updateDataMinMax('Col4', 'Value8', currentObj);
    });

    return Data_styling;
}
// ### Visualization Function ###
function renderDataVisualization(data, location,styling_option,Layout_dict,data_styling,label_conversion_dict,label_names) {

    // ######################
    // ## Initializzation  ##
    // ######################

    // Check the visibility of each layer
    const Layer1_isChecked = d3.select(`#toggle-layer-1`).property('checked');
    const Layer2_isChecked = d3.select(`#toggle-layer-2`).property('checked');
    const Layer3_isChecked = d3.select(`#toggle-layer-3`).property('checked');

    let Checklist = [Layer1_isChecked, Layer2_isChecked, Layer3_isChecked];

    // Indentation

    // const Indentation_toggle = d3.select(`#toggle-indentation`).property('checked');


    // ## Global Variables ##
    let columns = Layout_dict.columns;
    let col_breaker = Layout_dict.col_breaker;
    let col_max_label = Layout_dict.col_max_label;
    let indent = Layout_dict.indentation;


    // Calculate total count
    const total_count = calculateTotalCount(data,Checklist);

    // set col breaker number
    var Col_break_number;

    if (col_max_label == "Auto") {
        Col_break_number = Math.ceil(total_count / columns);
    } else if (col_max_label == "Custom") {
        Col_break_number = Layout_dict.Col_break_number
    }

    // Dimensions
    let width = 0;
    const height = 200;
    const margin = { top: 40, right: 20, bottom: 20, left: 20 };

    // ## Calculate all spacing and dimensions ##
    let spacing_dict = Calculate_dimension(data,Checklist,Col_break_number,columns,label_conversion_dict,label_names,styling_option,margin);
    
    // update spacing_dict
    spacing_dict.Col1 = spacing_dict.Col1 === -Infinity ? spacing_dict.Col1 : spacing_dict.Col1 + (indent == "No" ? 80 : 0); // 90 is the xoffset of the labels due to datapoints (this might change)
    spacing_dict.Col2 = spacing_dict.Col2 === -Infinity ? spacing_dict.Col2 : spacing_dict.Col2 + (indent == "No" ? 80 : 0);
    spacing_dict.Col3 = spacing_dict.Col3 === -Infinity ? spacing_dict.Col3 : spacing_dict.Col3 + (indent == "No" ? 80 : 0);
    spacing_dict.Col4 = spacing_dict.Col4 === -Infinity ? spacing_dict.Col4 : spacing_dict.Col4 + (indent == "No" ? 80 : 0);
    // console.log(spacing_dict);
    // Update width

    col_list = ['Col1','Col2','Col3','Col4'];
    for (let i = 0; i < columns; i++) {
        width = width+spacing_dict[col_list[i]];
    }
    width = width+45;
    
    // Recalculate width if its less than legend bars
    let number_of_bar_legends = 0;
    Object.keys(Data_styling).forEach(function(column) {
        if (Data_styling[column].Data === "Yes" && Data_styling[column].Datatype === 'Continuous') {
            number_of_bar_legends++;
        }
    });
    let minimum_width = 0;
    if (number_of_bar_legends > 0) {
        minimum_width = 230*number_of_bar_legends+margin.left+margin.right;
        if (width < minimum_width) {
            width = minimum_width;
        }
    }

    // X & Y coordinates
    let ylist_start = margin.top+5+80;
    let yOffset = ylist_start;
    let yOffset_max = 0
    let xOffset = 0;

    // Create svg element
    const svg = d3.select("#" + location)
        .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .attr("id", "ListPlot_visualization");


    // ###############
    // ## Functions ##
    // ###############


    // ## Function to calculate the total count of keys and values ##
    function calculateTotalCount(obj,checklist) {
        let totalCount = 0;
        // Recursive function to traverse the nested object
        for (const key in obj) {
            // Increment count for each key
            checklist[0] ? totalCount++ : totalCount;
            for (const subkey1 in obj[key]) {
                checklist[1] ? totalCount++ : totalCount;
                for (const subkey2 in obj[key][subkey1]) {
                    checklist[2] ? totalCount++ : totalCount;
                    for (const item in obj[key][subkey1][subkey2]) {
                        totalCount++;
                    }
                }
            }
        } return totalCount;
        }


    // ## Shape Funciton ##
    // Function to add different shapes
    function addShape(shapeType, x, y, size, fillColor) {
        switch (shapeType) {
            case 'circle':
                svg.append('circle')
                    .attr('cx', x)
                    .attr('cy', y)
                    .attr('r', size)
                    .style('stroke', 'black')
                    .style('stroke-width', 1)
                    .style('fill', fillColor);
                break;
            case 'rect':
                svg.append('rect')
                    .attr('x', x - size)
                    .attr('y', y - size)
                    .attr('width', size * 2)
                    .attr('height', size * 2)
                    .style('stroke', 'black')
                    .style('stroke-width', 1)
                    .style('fill', fillColor);
                break;
            case 'triangle':
                svg.append('path')
                    .attr('d', `M ${x} ${y - size} L ${x - size} ${y + size} L ${x + size} ${y + size} Z`)
                    .style('stroke', 'black')
                    .style('stroke-width', 1)
                    .style('fill', fillColor);
                break;
            case 'diamond':
                svg.append('path')
                    .attr('d', `M ${x} ${y - size} L ${x - size} ${y} L ${x} ${y + size} L ${x + size} ${y} Z`)
                    .style('stroke', 'black')
                    .style('stroke-width', 1)
                    .style('fill', fillColor);
                break;
            case 'star':
                const points = 5;
                const outerRadius = size;
                const innerRadius = size / 2.5;
                let starPath = '';
                for (let i = 0; i < points * 2; i++) {
                    const angle = (Math.PI / points) * i;
                    const radius = i % 2 === 0 ? outerRadius : innerRadius;
                    const xPoint = x + Math.cos(angle) * radius;
                    const yPoint = y + Math.sin(angle) * radius;
                    starPath += i === 0 ? `M ${xPoint} ${yPoint}` : `L ${xPoint} ${yPoint}`;
                }
                starPath += 'Z';
                svg.append('path')
                    .attr('d', starPath)
                    .style('stroke', 'black')
                    .style('stroke-width', 1)
                    .style('fill', fillColor);
                break;
            default:
                console.log('Unknown shape type');
        }
    }

    function Calculate_dimension(data,Checklist,Col_break_number,columns,label_conversion_dict,label_names,styling_option,margin) {

        temp_col_state = 1;

        label_max_dict = {col1_label_max: {'level1':-Infinity,'level2':-Infinity,'level3':-Infinity,'level4':-Infinity},
                            col2_label_max: {'level1':-Infinity,'level2':-Infinity,'level3':-Infinity,'level4':-Infinity},
                            col3_label_max: {'level1':-Infinity,'level2':-Infinity,'level3':-Infinity,'level4':-Infinity},
                            col4_label_max: {'level1':-Infinity,'level2':-Infinity,'level3':-Infinity,'level4':-Infinity}};

        Col_spacing_dict = {'Col1': -Infinity,'Col2': -Infinity,'Col3': -Infinity,'Col4': -Infinity};

        const label_dict_keys = Object.keys(label_max_dict);

        const htmlEntities = [
            '&alpha;', '&beta;', '&gamma;', '&delta;', '&epsilon;',
            '&zeta;', '&eta;', '&theta;', '&iota;', '&kappa;',
            '&lambda;', '&mu;', '&nu;', '&xi;', '&omicron;',
            '&pi;', '&rho;', '&sigma;', '&tau;', '&upsilon;',
            '&phi;', '&chi;', '&psi;', '&omega;'
        ];
        function replaceHtmlEntities(str) {
            // Create a regular expression to match any of the HTML entities as keys
            const entityRegex = new RegExp(htmlEntities.join('|'), 'gi');

            // Replace occurrences of HTML entities with a placeholder character (e.g., 'a')
            return str.replace(entityRegex, 'a');
        }


        // initialize variables for length calculation
        let label_dim_counter = 0;
        let label_length = 0;

        // Functions for length updates
        // Function to update column state and label max length
        function updateLabelMaxLength(key,level) {
            label_dim_counter++;

            // Update current column state
            if (label_dim_counter > Col_break_number && temp_col_state < columns) {
                label_dim_counter = 1;
                temp_col_state++;
            }

            
            if (level == 'level1') {
                // Trimming Class
                key = key.split(" (")[0];
            } else if (level == 'level2') {
                key = key;
            } else if (level == 'level3') {
                // Trimming receptor family
                key = key.replace(/( receptors|neuropeptide )/g, '');
                key = key.split(" (")[0];
            }

            // Dynamic key for the current column
            const colKey = `col${temp_col_state}_label_max`;
            label_length = key.length;

            // Update the max length for the current column
            if (label_max_dict[colKey][level] < label_length) {
                label_max_dict[colKey][level] = label_length;
            }
        }

        // Function to update label max length for array items with specific conversions
        function updateLabelMaxLengthForItems(item) {
            label_dim_counter++;

            // Update current column state
            if (label_dim_counter > Col_break_number && temp_col_state < columns) {
                label_dim_counter = 1;
                temp_col_state++;
            }

            // Dynamic key for the current column
            const colKey = `col${temp_col_state}_label_max`;
            let label = item;

            if (label_names == 'UniProt') {
                label = label_conversion_dict[item];
                // label = label.replace(/_human/g, '');
            } else if (label_names == 'IUPHAR') {
                label = item;
                label = replaceHtmlEntities(label)
                label = label.replace(/<\/?i>|(-adrenoceptor| receptor)|<\/?sub>/g, '');
            }

            label_length = label.length;
            // console.log(label);
            // Update the max length for the current column
            if (label_max_dict[colKey]['level4'] < label_length) {
                label_max_dict[colKey]['level4'] = label_length;
            }
        }
        Object.entries(data).forEach(([key, value]) => {
            // If Class is there
            if (Checklist[0]) {
                updateLabelMaxLength(key,'level1');
                // transverse next level 2
                if (Checklist[1] && typeof value === 'object') {
                    Object.entries(value).forEach(([subKey1, subValue1]) => {
                        updateLabelMaxLength(subKey1,'level2');
                        // transverse next level 3
                        if (Checklist[2] && typeof subValue1 === 'object') {
                            Object.entries(subValue1).forEach(([subKey2, subValue2]) => {
                                updateLabelMaxLength(subKey2,'level3');
                                // transverse final level 4
                                if (Array.isArray(subValue2)) {
                                    subValue2.forEach(item => {
                                        updateLabelMaxLengthForItems(item)
                                    });
                                }
                            });
                        } else if (!Checklist[2]) {
                            Object.entries(subValue1).forEach(([subKey2, subValue2]) => {
                                if (Array.isArray(subValue2)) {
                                    subValue2.forEach(item => {
                                        updateLabelMaxLengthForItems(item)
                                    });
                                }
                            });
                        }
                    });
                } else if (!Checklist[1] && Checklist[2]) {
                    // transverse final level 2
                    Object.entries(value).forEach(([subKey1, subValue1]) => {
                        if (Checklist[2] && typeof subValue1 === 'object') {
                            // transverse final level 3
                            Object.entries(subValue1).forEach(([subKey2, subValue2]) => {
                                updateLabelMaxLength(subKey2,'level3');
                                // transverse final level 4
                                if (Array.isArray(subValue2)) {
                                    subValue2.forEach(item => {
                                        updateLabelMaxLengthForItems(item)
                                    });
                                }
                            });
                        }
                    });
                } else if (!Checklist[1] && !Checklist[2]) {
                    // transverse final level 2
                    Object.entries(value).forEach(([subKey1, subValue1]) => {
                        // transverse final level 3
                        Object.entries(subValue1).forEach(([subKey2, subValue2]) => {
                            // transverse final level 4
                            if (Array.isArray(subValue2)) {
                                subValue2.forEach(item => {
                                    updateLabelMaxLengthForItems(item)
                                });
                            }
                        });
                    });
                }
            } else if (!Checklist[0] && Checklist[1]) {
                // transverse final level 2
                Object.entries(value).forEach(([subKey1, subValue1]) => {
                    updateLabelMaxLength(subKey1,'level2');
                    if (Checklist[2] && typeof subValue1 === 'object') {
                        // transverse final level 3
                        Object.entries(subValue1).forEach(([subKey2, subValue2]) => {
                            updateLabelMaxLength(subKey2,'level3');
                            // transverse final level 4
                            if (Array.isArray(subValue2)) {
                                subValue2.forEach(item => {
                                    updateLabelMaxLengthForItems(item)
                                });
                            }
                        });
                    } else if (!Checklist[2]) {
                        // transverse final level 3
                        Object.entries(subValue1).forEach(([subKey2, subValue2]) => {
                            // transverse final level 4
                            if (Array.isArray(subValue2)) {
                                subValue2.forEach(item => {
                                    updateLabelMaxLengthForItems(item)
                                });
                            }
                        });
                    }
                });
            } else if (!Checklist[0] && !Checklist[1] && Checklist[2]) {
                // transverse final level 2
                Object.entries(value).forEach(([subKey1, subValue1]) => {
                    // transverse final level 3
                    Object.entries(subValue1).forEach(([subKey2, subValue2]) => {
                        updateLabelMaxLength(subKey2,'level3');
                        // transverse final level 4
                        if (Array.isArray(subValue2)) {
                            subValue2.forEach(item => {
                                updateLabelMaxLengthForItems(item)
                            });
                        }
                    });
                });
            } else if (!Checklist[0] && !Checklist[1] && !Checklist[2]){
                // transverse final level 2
                Object.entries(value).forEach(([subKey1, subValue1]) => {
                    // transverse final level 3
                    Object.entries(subValue1).forEach(([subKey2, subValue2]) => {
                        // transverse final level 4
                        if (Array.isArray(subValue2)) {
                            subValue2.forEach(item => {
                                updateLabelMaxLengthForItems(item)
                            });
                        }
                    });
                });
            }
        });

        // calculate the longest entry for each column
        const levels = ['level1', 'level2', 'level3', 'level4'];
        const laters = ['Layer1', 'Layer2', 'Layer3', 'Layer4'];
        const cols = ['Col1', 'Col2', 'Col3', 'Col4'];


        for (let i = 0; i < columns; i++) {
            const key = label_dict_keys[i];
            const max_label_values = label_max_dict[key];
            // console.log(label_dict_keys);
            for (let j = 0; j < levels.length; j++) {
                const levelKey = levels[j];
                const layerKey = laters[j];
                const columnKey = cols[i];

                if (max_label_values[levelKey] !== -Infinity) {
                    // Create a dummy text element to measure its size
                    const dummyText = d3.select("body")
                                    .append("svg")
                                    .attr("class", "dummy-text")
                                    .append("text")
                                    .attr("font-size", styling_option[layerKey].Fontsize)
                                    .attr("font-weight", styling_option[layerKey].Bold ? "bold" : "normal")
                                    .text("X".repeat(max_label_values[levelKey]));

                    // Measure the bounding box of the dummy text
                    const bbox = dummyText.node().getBBox();
                    let estimatedLength = 0;
                    if (layerKey == 'Layer4') {
                        estimatedLength = bbox.width*0.8+margin.left+80;
                    } else {
                        estimatedLength = bbox.width*0.8+margin.left;
                    }

                    // Remove the dummy text element
                    d3.select(".dummy-text").remove();

                    // Update Col_spacing_dict if the estimated length is greater
                    if (estimatedLength > Col_spacing_dict[columnKey]) {
                        Col_spacing_dict[columnKey] = estimatedLength;
                    }
                }
            }
        }
        
        return Col_spacing_dict
    }
    

    // ## Column breaker function ##
    function ColumnsBreakerFunc (columns,label_counter,Col_break_number,state) {
        if (columns == 2 && label_counter > Col_break_number && state.current_col == 1){
            state.current_col = 2;
            xOffset = spacing_dict.Col1;
            if (yOffset > yOffset_max) {yOffset_max = yOffset;}
            yOffset = ylist_start;
        }
        if (columns == 3 && label_counter > Col_break_number && state.current_col == 1){
            state.current_col = 2;
            xOffset = spacing_dict.Col1;
            if (yOffset > yOffset_max) {yOffset_max = yOffset;}
            yOffset = ylist_start;
        } else if (columns == 3 && label_counter > Col_break_number*2 && state.current_col == 2){
            state.current_col = 3;
            xOffset =  spacing_dict.Col1+spacing_dict.Col2;
            if (yOffset > yOffset_max) {yOffset_max = yOffset;}
            yOffset = ylist_start;
        }
        if (columns == 4 && label_counter > Col_break_number && state.current_col == 1){
            state.current_col = 2;
            xOffset = spacing_dict.Col1;
            if (yOffset > yOffset_max) {yOffset_max = yOffset;}
            yOffset = ylist_start;
        } else if (columns == 4 && label_counter > Col_break_number*2 && state.current_col == 2){
            state.current_col = 3;
            xOffset = spacing_dict.Col1+spacing_dict.Col2;
            if (yOffset > yOffset_max) {yOffset_max = yOffset;}
            yOffset = ylist_start;
        } else if (columns == 4 && label_counter > Col_break_number*3 && state.current_col == 3){
            state.current_col = 4;
            xOffset = spacing_dict.Col1+spacing_dict.Col2+spacing_dict.Col3;
            if (yOffset > yOffset_max) {yOffset_max = yOffset;}
            yOffset = ylist_start;
        }
    }

    // #####################
    // ## List processing ##
    // #####################

    function processList(list) {
        // Check the visibility of each layer
        const Layer1_isChecked = d3.select(`#toggle-layer-1`).property('checked');
        const Layer2_isChecked = d3.select(`#toggle-layer-2`).property('checked');
        const Layer3_isChecked = d3.select(`#toggle-layer-3`).property('checked');

        let Checklist = [Layer1_isChecked, Layer2_isChecked, Layer3_isChecked];


        let col_max_label = Layout_dict.col_max_label;
        var Col_break_number = 20;

        if (col_max_label == "Auto") {
            Col_break_number = Math.ceil(total_count / columns);
        } else if (col_max_label == "Custom") {
            Col_break_number = Layout_dict.Col_break_number
        }

        var state = { current_col: 1 };
        // Setup amount of cols
        let label_counter = 1;

        // Set color gradient based on input

        // Define Color_dict to store color scales for each column and shape
        var Color_dict = {};

        // # Setup color styling if "All" is chosen #
        Object.keys(data_styling).forEach(function(col) {
            Color_dict[col] = {};
            if (data_styling[col].data_color_complexity == 'One') {
                Color_dict[col]['All'] = d3.scale.linear().domain([data_styling[col].Data_min,data_styling[col].Data_max]).range(["#FFFFFF", data_styling[col].Data_color1]);
            } else if (data_styling[col].data_color_complexity == 'Two') {
                Color_dict[col]['All'] = d3.scale.linear()
                .domain([data_styling[col].Data_min, (data_styling[col].Data_min + data_styling[col].Data_max) / 2, data_styling[col].Data_max])
                .range([data_styling[col].Data_color1, "#FFFFFF", data_styling[col].Data_color2]);
            }
        });

        // ##############################
        // ## Handle all list elements ##
        // ##############################

        Object.entries(list).forEach(([key, value]) => {
            // Text input function
            function add_text(label, yOffset, layer, bold = false, italic = false, underline = false, fontSize = '16px', color = 'black') {

                // Function to break label into multiple lines based on length
                function breakLabel(label) {
                    const maxLength = 15;
                    let result = '';
                    let start = 0;

                    while (start < label.length) {
                        if (label.length - start <= maxLength) {
                            result += label.slice(start);
                            break;
                        }

                        let breakPosition = label.lastIndexOf(' ', start + maxLength);
                        if (breakPosition === -1 || breakPosition < start) {
                            breakPosition = label.indexOf('/', start);
                            if (breakPosition === -1 || breakPosition > start + maxLength) {
                                breakPosition = start + maxLength; // Default break position if no space or '/' found
                            }
                        }

                        result += label.slice(start, breakPosition) + '\n';
                        start = breakPosition + 1;
                    }

                    return result;
                }

                label_key = label;

                // # Different types of label handling #
                if (label_names == 'IUPHAR') {

                    const htmlEntities = {
                        '&alpha;': 'α',
                        '&beta;': 'β',
                        '&gamma;': 'γ',
                        '&delta;': 'δ',
                        '&epsilon;': 'ε',
                        '&zeta;': 'ζ',
                        '&eta;': 'η',
                        '&theta;': 'θ',
                        '&iota;': 'ι',
                        '&kappa;': 'κ',
                        '&lambda;': 'λ',
                        '&mu;': 'μ',
                        '&nu;': 'ν',
                        '&xi;': 'ξ',
                        '&omicron;': 'ο',
                        '&pi;': 'π',
                        '&rho;': 'ρ',
                        '&sigma;': 'σ',
                        '&tau;': 'τ',
                        '&upsilon;': 'υ',
                        '&phi;': 'φ',
                        '&chi;': 'χ',
                        '&psi;': 'ψ',
                        '&omega;': 'ω'
                    };

                    // Replace HTML entities with corresponding Unicode characters
                    for (const [entity, char] of Object.entries(htmlEntities)) {
                        label = label.replace(new RegExp(entity, 'g'), char);
                    }

                    const textElement = svg.append('text')
                        .attr('x', margin.left + xOffset + 80 + ((indent === 'Yes' && layer === 'Layer-4') ? 85 : 0))
                        .attr('y', yOffset)
                        .attr('class', layer)
                        .style('dominant-baseline', 'middle') // Set vertical alignment
                        .style('font-weight', bold ? 'bold' : 'normal')
                        .style('font-style', italic ? 'italic' : 'normal')
                        .style('text-decoration', underline ? 'underline' : 'none')
                        .style('font-size', fontSize)
                        .style('fill', color); // Set text color

                    // Calculate the subscript font size (e.g., 75% of the main font size)
                    const mainFontSize = parseFloat(fontSize);
                    const subFontSize = mainFontSize * 0.75 + 'px';
                    
                    // Trimming Class
                    if (layer == 'Layer-1') {
                        label = label.split(" (")[0];
                    }

                    // Trimming receptor family
                    if (layer == 'Layer-3') {
                        label = label.replace(/( receptors|neuropeptide )/g, '');
                        label = label.split(" (")[0];
                    }

                    // Remove <i> and </i> tags from the label
                    // label = label.replace(/<\/?i>/g, '');
                    if (layer == 'Layer-4') {
                        label = label.replace(/<\/?i>|(-adrenoceptor| receptor)/g, '');
                    }
                    const parts = label.split(/(<sub>|<\/sub>)/); // Split label into parts including <sub> tags
                    let inSub = false; // Flag to track if we are inside a subscript

                    // Handle subscripts
                    parts.forEach(part => {
                        if (part === '<sub>') {
                            inSub = true;
                        } else if (part === '</sub>') {
                            inSub = false;
                        } else {
                            const tspan = textElement.append('tspan')
                                .text(part);

                            if (inSub) {
                                tspan.attr('dy', '0.3em') // Subscript positioning
                                    .attr('font-size', subFontSize);
                            } else {
                                tspan.attr('dy', '-0.3em') // Reset positioning
                                    .attr('font-size', fontSize);
                            }
                        }
                    });

                } else if (label_names == 'UniProt') {

                    // Trimming Class
                    if (layer == 'Layer-1') {
                        label = label_key.split(" (")[0];
                    } else if (layer == 'Layer-3') {
                        label = label_key.replace(/( receptors|neuropeptide )/g, '');
                        label = label.split(" (")[0];
                    // } If UniProt, handle receptor names but keep classes etc. #
                    } else if (layer == 'Layer-4') {
                    label = label_conversion_dict[label_key];
                    label = label.replace(/_human/g, '').toUpperCase();
                    } else {
                        label = label_key
                    }

                    const textElement = svg.append('text')
                        .attr('x', margin.left + xOffset + 80 + ((indent === 'Yes' && layer === 'Layer-4') ? 85 : 0))
                        .attr('y', yOffset)
                        .attr('class', layer)
                        .attr('dy', '-0.3em') // Adjust this value to move the text higher
                        .style('dominant-baseline', 'middle') // Set vertical alignment
                        .style('font-weight', bold ? 'bold' : 'normal')
                        .style('font-style', italic ? 'italic' : 'normal')
                        .style('text-decoration', underline ? 'underline' : 'none')
                        .style('font-size', fontSize)
                        .style('fill', color)
                        .text(label);


                }

                // #############################
                // ### Implement data points ###
                // #############################

                if (layer == 'Layer-4') {
                    if (list_data_wow.hasOwnProperty(label_key)) {

                        // Initialize all variables for better overview

                        // Col data checker
                        var Col1_data_checker = data_styling.Col1.Data == "Yes";
                        var Col2_data_checker = data_styling.Col2.Data == "Yes";
                        var Col3_data_checker = data_styling.Col3.Data == "Yes";
                        var Col4_data_checker = data_styling.Col4.Data == "Yes";

                        // Col lebels
                        var Col1_shape = list_data_wow[label_key].hasOwnProperty('Value1') ? list_data_wow[label_key].Value1.toLowerCase() : false;
                        var Col2_shape = list_data_wow[label_key].hasOwnProperty('Value3') ? list_data_wow[label_key].Value3.toLowerCase() : false;
                        var Col3_shape = list_data_wow[label_key].hasOwnProperty('Value5') ? list_data_wow[label_key].Value5.toLowerCase() : false;
                        var Col4_shape = list_data_wow[label_key].hasOwnProperty('Value7') ? list_data_wow[label_key].Value7.toLowerCase() : false;

                        // Col data shapes checker
                        var Col1_data = list_data_wow[label_key].hasOwnProperty('Value2') ? list_data_wow[label_key].Value2 : false;
                        var Col2_data = list_data_wow[label_key].hasOwnProperty('Value4') ? list_data_wow[label_key].Value4 : false;
                        var Col3_data = list_data_wow[label_key].hasOwnProperty('Value6') ? list_data_wow[label_key].Value6 : false;
                        var Col4_data = list_data_wow[label_key].hasOwnProperty('Value8') ? list_data_wow[label_key].Value8 : false;

                        // Col Datatype
                        var Col1_datatype = data_styling.Col1.Datatype;
                        var Col2_datatype = data_styling.Col2.Datatype;
                        var Col3_datatype = data_styling.Col3.Datatype;
                        var Col4_datatype = data_styling.Col4.Datatype;

                        // Col data coloring scheme

                        const Col1_coloring = data_styling.Col1.Data_coloring;
                        const Col2_coloring = data_styling.Col2.Data_coloring;
                        const Col3_coloring = data_styling.Col3.Data_coloring;
                        const Col4_coloring = data_styling.Col4.Data_coloring;

                        const Color_list = ['black', 'red', 'blue', 'green'];
                        const Shape_list = ['circle','rect','triangle','star','diamond']

                        const col1_XoffSet = 0 + (indent === 'Yes' ? 85 : 0)
                        const col2_XoffSet = 20 + (indent === 'Yes' ? 85 : 0)
                        const col3_XoffSet = 40 + (indent === 'Yes' ? 85 : 0)
                        const col4_XoffSet = 60 + (indent === 'Yes' ? 85 : 0)

                        // ## Col 1 ##
                        if (Col1_data_checker && (Col1_shape || Col1_data)) {
                            if (Col1_shape) {
                                if (Col1_data && Col1_datatype == 'Discrete') {
                                    shape_color = Color_list.includes(Col1_data.toLowerCase()) ? Col1_data : 'black';
                                } else if (Col1_data && Col1_datatype == 'Continuous') {
                                    if (Col1_coloring == 'All') {
                                        shape_color = Color_dict.Col1.All(Col1_data);
                                    } else if (Col1_coloring == 'Shapes') {
                                        shape_color = Color_dict['Col1'][Col1_shape](Col1_data);
                                    }
                                } else {
                                    shape_color = 'black';
                                }
                                shape_value = Col1_shape === "square" ? 'rect' : Col1_shape;
                                addShape(Shape_list.includes(shape_value) ? shape_value : 'circle', margin.left + xOffset  + col1_XoffSet, yOffset - (parseFloat(fontSize) * 0.50), 6, shape_color);
                            } else {
                                if (Col1_data && Col1_datatype == 'Discrete') {
                                    shape_color = Color_list.includes(Col1_data.toLowerCase()) ? Col1_data : 'black';
                                    addShape('circle', margin.left + xOffset + col1_XoffSet, yOffset - (parseFloat(fontSize) * 0.50), 6, shape_color);
                                } else if (Col1_data && Col1_datatype == 'Continuous') {
                                    shape_color = Color_dict.Col1.All(Col1_data);
                                    addShape('circle', margin.left + xOffset + col1_XoffSet, yOffset - (parseFloat(fontSize) * 0.50), 6, shape_color);
                                }
                            }
                        }

                        // ## Col 2 ##
                        if (Col2_data_checker && (Col2_shape || Col2_data)) {
                            if (Col2_shape) {


                                if (Col2_data && Col2_datatype == 'Discrete') {
                                    shape_color = Color_list.includes(Col2_data.toLowerCase()) ? Col2_data : 'black';
                                } else if (Col2_data && Col2_datatype == 'Continuous') {

                                    if (Col2_coloring == 'All') {
                                        shape_color = Color_dict.Col2.All(Col2_data);
                                    } else if (Col2_coloring == 'Shapes') {
                                        shape_color = Color_dict['Col2'][Col2_shape](Col2_data);
                                    }
                                } else {
                                    shape_color = 'black';
                                }
                                shape_value = Col2_shape === "square" ? 'rect' : Col2_shape;
                                addShape(Shape_list.includes(shape_value) ? shape_value : 'circle', margin.left + xOffset + col2_XoffSet, yOffset - (parseFloat(fontSize) * 0.50), 6, shape_color);
                            } else {
                                if (Col2_data && Col2_datatype == 'Discrete') {
                                    shape_color = Color_list.includes(Col2_data.toLowerCase()) ? Col2_data : 'black';
                                    addShape('circle', margin.left + xOffset + col2_XoffSet, yOffset - (parseFloat(fontSize) * 0.50), 6, shape_color);
                                } else if (Col2_data && Col2_datatype == 'Continuous') {
                                    shape_color = Color_dict.Col2.All(Col2_data);
                                    addShape('circle', margin.left + xOffset + col2_XoffSet, yOffset - (parseFloat(fontSize) * 0.50), 6, shape_color);
                                }
                            }
                        }

                        // ## Col 3 ##
                        if (Col3_data_checker && (Col3_shape || Col3_data)) {
                            if (Col3_shape) {
                                if (Col3_data && Col3_datatype == 'Discrete') {
                                    shape_color = Color_list.includes(Col3_data.toLowerCase()) ? Col3_data : 'black';
                                } else if (Col3_data && Col3_datatype == 'Continuous') {
                                    if (Col3_coloring == 'All') {
                                        shape_color = Color_dict.Col3.All(Col3_data);
                                    } else if (Col3_coloring == 'Shapes') {
                                        shape_color = Color_dict['Col3'][Col3_shape](Col3_data);
                                    }
                                } else {
                                    shape_color = 'black';
                                }
                                shape_value = Col3_shape === "square" ? 'rect' : Col3_shape;
                                addShape(Shape_list.includes(shape_value) ? shape_value : 'circle', margin.left + xOffset + col3_XoffSet, yOffset - (parseFloat(fontSize) * 0.50), 6, shape_color);
                            } else {
                                if (Col3_data && Col3_datatype == 'Discrete') {
                                    shape_color = Color_list.includes(Col3_data.toLowerCase()) ? Col3_data : 'black';
                                    addShape('circle', margin.left + xOffset + col3_XoffSet, yOffset - (parseFloat(fontSize) * 0.50), 6, shape_color);
                                } else if (Col3_data && Col3_datatype == 'Continuous') {
                                    shape_color = Color_dict.Col3.All(Col3_data);
                                    addShape('circle', margin.left + xOffset + col3_XoffSet, yOffset - (parseFloat(fontSize) * 0.50), 6, shape_color);
                                }
                            }
                        }

                        // ## Col 4 ##
                        if (Col4_data_checker && (Col4_shape || Col4_data)) {
                            if (Col4_shape) {
                                if (Col4_data && Col4_datatype == 'Discrete') {
                                    shape_color = Color_list.includes(Col4_data.toLowerCase()) ? Col4_data : 'black';
                                } else if (Col4_data && Col4_datatype == 'Continuous') {
                                    if (Col4_coloring == 'All') {
                                        shape_color = Color_dict.Col4.All(Col4_data);
                                    } else if (Col4_coloring == 'Shapes') {
                                        shape_color = Color_dict['Col4'][Col4_shape](Col4_data);
                                    }
                                } else {
                                    shape_color = 'black';
                                }
                                shape_value = Col4_shape === "square" ? 'rect' : Col4_shape;
                                addShape(Shape_list.includes(shape_value) ? shape_value : 'circle', margin.left + xOffset + col4_XoffSet, yOffset - (parseFloat(fontSize) * 0.50), 6, shape_color);
                            } else {
                                if (Col4_data && Col4_datatype == 'Discrete') {
                                    shape_color = Color_list.includes(Col4_data.toLowerCase()) ? Col4_data : 'black';
                                    addShape('circle', margin.left + xOffset + col4_XoffSet, yOffset - (parseFloat(fontSize) * 0.50), 6, shape_color);
                                } else if (Col4_data && Col4_datatype == 'Continuous') {
                                    shape_color = Color_dict.Col4.All(Col4_data);
                                    addShape('circle', margin.left + xOffset + col4_XoffSet, yOffset - (parseFloat(fontSize) * 0.50), 6, shape_color);
                                }
                            }
                        }
                    } else {
                        console.log(label_key);
                    }
                }
            }

            // # Check and append text for Layer 1 breaker if it's visible #

            if (Checklist[0]) {
                if (col_breaker == "Layer1" || col_breaker == "Custom") {
                    ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                }
                add_text(key, yOffset, "Layer-1",styling_option.Layer1.Bold,styling_option.Layer1.Italic,styling_option.Layer1.Underline,styling_option.Layer1.Fontsize,styling_option.Layer1.Color);
                yOffset += 30;
                label_counter++;

                // Check and append text for Layer 2 if it's visible and has sublayers
                if (Checklist[1] && typeof value === 'object') {
                    Object.entries(value).forEach(([subKey, subValue]) => {
                        if (col_breaker == "Layer2" || col_breaker == "Custom") {
                            ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                        }
                        add_text(subKey, yOffset, "Layer-2",styling_option.Layer2.Bold,styling_option.Layer2.Italic,styling_option.Layer2.Underline,styling_option.Layer2.Fontsize,styling_option.Layer2.Color);
                        yOffset += 30;
                        label_counter++;

                        // Check and append text for Layer 3 if it's visible and has sublayers
                        if (Checklist[2] && typeof subValue === 'object') {
                            Object.entries(subValue).forEach(([subSubKey, subSubValue]) => {
                                if (col_breaker == "Layer3" || col_breaker == "Custom") {
                                    ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                                }
                                add_text(subSubKey, yOffset, "Layer-3",styling_option.Layer3.Bold,styling_option.Layer3.Italic,styling_option.Layer3.Underline,styling_option.Layer3.Fontsize,styling_option.Layer3.Color);
                                yOffset += 30;
                                label_counter++;

                                // Check and append text for Layer 4 if it's visible and is an array
                                if (Array.isArray(subSubValue)) {
                                    subSubValue.forEach(item => {
                                        if (col_breaker == "Layer4" || col_breaker == "Custom") {
                                            ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                                        }
                                        add_text(item, yOffset, "Layer-4",styling_option.Layer4.Bold,styling_option.Layer4.Italic,styling_option.Layer4.Underline,styling_option.Layer4.Fontsize,styling_option.Layer4.Color);

                                        yOffset += 30;
                                        label_counter++;

                                    });
                                }
                            });
                        } else {
                            Object.entries(subValue).forEach(([subSubKey, subSubValue]) => {
                                // Check and append text for Layer 4 if it's visible and is an array
                                if (Array.isArray(subSubValue)) {
                                    subSubValue.forEach(item => {
                                        if (col_breaker == "Layer4" || col_breaker == "Custom") {
                                            ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                                        }
                                        add_text(item, yOffset, "Layer-4",styling_option.Layer4.Bold,styling_option.Layer4.Italic,styling_option.Layer4.Underline,styling_option.Layer4.Fontsize,styling_option.Layer4.Color);
                                        yOffset += 30;
                                        label_counter++;

                                    });
                                }
                            });
                        }
                    });
                } else {
                    Object.entries(value).forEach(([subKey, subValue]) => {
                        // Check and append text for Layer 3 if it's visible and has sublayers
                        if (Checklist[2] && typeof subValue === 'object') {
                            Object.entries(subValue).forEach(([subSubKey, subSubValue]) => {
                                if (col_breaker == "Layer3" || col_breaker == "Custom") {
                                    ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                                }
                                add_text(subSubKey, yOffset, "Layer-3",styling_option.Layer3.Bold,styling_option.Layer3.Italic,styling_option.Layer3.Underline,styling_option.Layer3.Fontsize,styling_option.Layer3.Color);
                                yOffset += 30;
                                label_counter++;

                                // Check and append text for Layer 4 if it's visible and is an array
                                if (Array.isArray(subSubValue)) {
                                    subSubValue.forEach(item => {
                                        if (col_breaker == "Layer4" || col_breaker == "Custom") {
                                            ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                                        }
                                        add_text(item, yOffset, "Layer-4",styling_option.Layer4.Bold,styling_option.Layer4.Italic,styling_option.Layer4.Underline,styling_option.Layer4.Fontsize,styling_option.Layer4.Color);
                                        yOffset += 30;
                                        label_counter++;

                                    });
                                }
                            });
                        } else {
                            Object.entries(subValue).forEach(([subSubKey, subSubValue]) => {
                                // Check and append text for Layer 4 if it's visible and is an array
                                if (Array.isArray(subSubValue)) {
                                    subSubValue.forEach(item => {
                                        if (col_breaker == "Layer4" || col_breaker == "Custom") {
                                            ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                                        }
                                        add_text(item, yOffset, "Layer-4",styling_option.Layer4.Bold,styling_option.Layer4.Italic,styling_option.Layer4.Underline,styling_option.Layer4.Fontsize,styling_option.Layer4.Color);
                                        yOffset += 30;
                                        label_counter++;

                                    });
                                }
                            });
                        }
                    });
                }
            } else if (!Checklist[0] && Checklist[1]) {
                // Check and append text for Layer 2 if it's visible and has sublayers
                if (Checklist[1] && typeof value === 'object') {
                    Object.entries(value).forEach(([subKey, subValue]) => {
                        if (col_breaker == "Layer2" || col_breaker == "Custom") {
                            ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                        }
                        add_text(subKey, yOffset, "Layer-2",styling_option.Layer2.Bold,styling_option.Layer2.Italic,styling_option.Layer2.Underline,styling_option.Layer2.Fontsize,styling_option.Layer2.Color);
                        yOffset += 30;
                        label_counter++;

                        // Check and append text for Layer 3 if it's visible and has sublayers
                        if (Checklist[2] && typeof subValue === 'object') {
                            Object.entries(subValue).forEach(([subSubKey, subSubValue]) => {
                                if (col_breaker == "Layer3" || col_breaker == "Custom") {
                                    ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                                }
                                add_text(subSubKey, yOffset, "Layer-3",styling_option.Layer3.Bold,styling_option.Layer3.Italic,styling_option.Layer3.Underline,styling_option.Layer3.Fontsize,styling_option.Layer3.Color);
                                yOffset += 30;
                                label_counter++;

                                // Check and append text for Layer 4 if it's visible and is an array
                                if (Array.isArray(subSubValue)) {
                                    subSubValue.forEach(item => {
                                        if (col_breaker == "Layer4" || col_breaker == "Custom") {
                                            ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                                        }
                                        add_text(item, yOffset, "Layer-4",styling_option.Layer4.Bold,styling_option.Layer4.Italic,styling_option.Layer4.Underline,styling_option.Layer4.Fontsize,styling_option.Layer4.Color);
                                        yOffset += 30;
                                        label_counter++;

                                    });
                                }
                            });
                        } else {
                            Object.entries(subValue).forEach(([subSubKey, subSubValue]) => {
                                // Check and append text for Layer 4 if it's visible and is an array
                                if (Array.isArray(subSubValue)) {
                                    subSubValue.forEach(item => {
                                        if (col_breaker == "Layer4" || col_breaker == "Custom") {
                                            ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                                        }
                                        add_text(item, yOffset, "Layer-4",styling_option.Layer4.Bold,styling_option.Layer4.Italic,styling_option.Layer4.Underline,styling_option.Layer4.Fontsize,styling_option.Layer4.Color);
                                        yOffset += 30;
                                        label_counter++;

                                    });
                                }
                            });
                        }
                    });
                }
            } else if (!Checklist[0] && !Checklist[1] && Checklist[2]) {
                Object.entries(value).forEach(([subKey, subValue]) => {
                    if (Checklist[2] && typeof subValue === 'object') {
                        Object.entries(subValue).forEach(([subSubKey, subSubValue]) => {
                            if (col_breaker == "Layer3" || col_breaker == "Custom") {
                                ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                            }
                            add_text(subSubKey, yOffset, "Layer-3",styling_option.Layer3.Bold,styling_option.Layer3.Italic,styling_option.Layer3.Underline,styling_option.Layer3.Fontsize,styling_option.Layer3.Color);
                            yOffset += 30;
                            label_counter++;

                            // Check and append text for Layer 4 if it's visible and is an array
                            if (Array.isArray(subSubValue)) {
                                subSubValue.forEach(item => {
                                    if (col_breaker == "Layer4" || col_breaker == "Custom") {
                                        ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                                    }
                                    add_text(item, yOffset, "Layer-4",styling_option.Layer4.Bold,styling_option.Layer4.Italic,styling_option.Layer4.Underline,styling_option.Layer4.Fontsize,styling_option.Layer4.Color);
                                    yOffset += 30;
                                    label_counter++;

                                });
                            }
                        });
                    }
                });
            } else if (!Checklist[0] && !Checklist[1] && !Checklist[2]) {
                Object.entries(value).forEach(([subKey, subValue]) => {
                    Object.entries(subValue).forEach(([subSubKey, subSubValue]) => {
                        if (Array.isArray(subSubValue)) {
                            subSubValue.forEach(item => {
                                if (col_breaker == "Layer4" || col_breaker == "Custom") {
                                    ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                                }
                                add_text(item, yOffset, "Layer-4",styling_option.Layer4.Bold,styling_option.Layer4.Italic,styling_option.Layer4.Underline,styling_option.Layer4.Fontsize,styling_option.Layer4.Color);
                                yOffset += 30;
                                label_counter++;

                            });
                        }
                    });
                });
            }
        });

    // # check if y-off-set is more than max, if so expand canvas #
    if (yOffset_max < yOffset) {
        yOffset_max = yOffset;
    }

    // Rerender height of plot as the last thing
    svg.attr("height",yOffset_max + margin.top + margin.bottom)

    // ## Add legend bars on top ##
    let bar_index = 0
    Object.keys(Data_styling).forEach(function(column) {
        if (Data_styling[column].Data === "Yes" && Data_styling[column].Datatype === 'Continuous') {
            const legendWidth = 200; // Width of the legend bar
            const data_fontsize = 14; // Adjust as needed
            const lowest_value = Data_styling[column].Data_min;
            const highest_value = Data_styling[column].Data_max;
            const spacing_bar = 30;
            const bar_height = 20;
            const text_off_set = 35;
            const x_off_set = 100;
            const y_off_set = 20;
            // Calculate the x position for the current bar
            const x_position = bar_index * (legendWidth + spacing_bar)+x_off_set;

            // Increment bar index for the next bar
            bar_index++;

            const legend_svg = svg.append("g"); // Append group for the legend bar

            const gradientId = `Gradient_${column}`;
            const defs = svg.append('defs');

            // Add centered text on top of the bar specifying the data column
            legend_svg.append("text")
            .attr('x', x_position + legendWidth / 2) // Center the text horizontally over the bar
            .attr('y', bar_height - 10 + y_off_set) // Adjust Y position as needed, e.g., 10px above the bar
            .style("font-size", `${data_fontsize}px`)
            .style("font-family", "sans-serif")
            .style("text-anchor", "middle") // Center the text horizontally
            .text(`Data column ${column.substr(3)}`); // Use the column name as the text

            // Gradient definition
            const gradient = defs.append('linearGradient')
            .attr('id', gradientId)
            .attr('x1', '0%')
            .attr('x2', '100%')
            .attr('y1', '0%')
            .attr('y2', '0%');

            // Starting color stop
            gradient.append('stop')
            .attr('offset', '0%');

            // Adjust gradient stops based on color complexity
            if (Data_styling[column].data_color_complexity !== 'One') {
                gradient.select('stop')
                    .attr('stop-color', Data_styling[column].Data_color1); // Start with color1 at 0%

                gradient.append('stop')
                    .attr('offset', '100%')
                    .attr('stop-color', Data_styling[column].Data_color2); // End with color2 at 100%
            } else {
                gradient.select('stop')
                    .attr('stop-color', '#FFFFFF'); // Start with white at 0%

                gradient.append('stop')
                    .attr('offset', '100%')
                    .attr('stop-color', Data_styling[column].Data_color1); // End with color1 at 100%
            }

            // Border around the rectangle
            legend_svg.append('rect')
            .attr('x', x_position) // Align with the gradient rectangle
            .attr('y', bar_height + y_off_set) // Align with the gradient rectangle
            .attr('width', legendWidth) // Same width as the gradient rectangle
            .attr('height', bar_height) // Same height as the gradient rectangle
            .style('fill', 'none')
            .style('stroke', 'black')
            .style('stroke-width', 2);

            // Rectangle for the gradient
            legend_svg.append('rect')
            .attr('x', x_position) // Horizontal spacing between bars
            .attr('y', bar_height + y_off_set) // Adjust Y position as needed
            .attr('width', legendWidth) // Legend bar width
            .attr('height', bar_height) // Legend bar height
            .style('fill', `url(#${gradientId})`);

            // Minimum value text
            legend_svg.append("text")
            .attr('x', x_position) // Align with the bar
            .attr('y', bar_height+text_off_set + y_off_set)
            .style("font-size", `${data_fontsize}px`)
            .style("font-family", "sans-serif")
            .text(lowest_value);

            // Middle value text (if applicable)
            if (Data_styling[column].data_color_complexity === 'Three') {
                legend_svg.append("text")
                .attr('x', x_position + legendWidth / 2)
                .attr('y', bar_height+text_off_set + y_off_set)
                .style("font-size", `${data_fontsize}px`)
                .style("font-family", "sans-serif")
                .style("text-anchor", "middle")
                .text((highest_value + lowest_value) / 2); // Middle value
            }

            // Maximum value text
            legend_svg.append("text")
            .attr('x', x_position + legendWidth)
            .attr('y', bar_height+text_off_set + y_off_set)
            .style("font-size", `${data_fontsize}px`)
            .style("font-family", "sans-serif")
            .style("text-anchor", "end")
            .text(highest_value);
        }
    });


}
    // # Run list process and return svg #
    processList(data,styling_option);
    return svg;
}

// #################
// ###  HEATMAP  ###
// #################


function heatmap_DataStyling() {
    var heatmap_DataStyling = {Number_of_colors: "Three",
        min_color: "#1a80bb",
        middle_color: "#FFFFFF",
        max_color: "#a00000",
        rotation: 90,
        label_position: 'Bottom',
        label_fontsize: 12,
        receptor_fontsize: 12,
        datalabels: true,
        data_border: true,
        data_fontsize: 12,
        legend_label: "Value intensity",
        LabelType: 'UniProt'

    }
    return heatmap_DataStyling
}

// Function to handle row labels based on label type
function handleRowLabels(textElement, label, labelType, fontSize) {
    const htmlEntities = {
        '&alpha;': 'α',
        '&beta;': 'β',
        '&gamma;': 'γ',
        '&delta;': 'δ',
        '&epsilon;': 'ε',
        '&zeta;': 'ζ',
        '&eta;': 'η',
        '&theta;': 'θ',
        '&iota;': 'ι',
        '&kappa;': 'κ',
        '&lambda;': 'λ',
        '&mu;': 'μ',
        '&nu;': 'ν',
        '&xi;': 'ξ',
        '&omicron;': 'ο',
        '&pi;': 'π',
        '&rho;': 'ρ',
        '&sigma;': 'σ',
        '&tau;': 'τ',
        '&upsilon;': 'υ',
        '&phi;': 'φ',
        '&chi;': 'χ',
        '&psi;': 'ψ',
        '&omega;': 'ω'
    };

    // Handle label transformations based on labelType
    let transformedLabel = label; // Initialize transformedLabel with the original label

    if (labelType === 'UniProt') {
        transformedLabel = label.replace(/_human/g, '').toUpperCase();
        textElement.text(transformedLabel);
    } else if (labelType === 'IUPHAR') {

        // Retrieve the IUPHAR label from your converter object
        transformedLabel = label_converter.UniProt_to_IUPHAR_converter[label];
        // Replace HTML entities with corresponding Unicode characters
        for (const [entity, char] of Object.entries(htmlEntities)) {
            transformedLabel = transformedLabel.replace(new RegExp(entity, 'g'), char);
        }
        // Remove <i> and </i> tags from the label
        transformedLabel = transformedLabel.replace(/<\/?i>/g, '');
        // remove receptor and adrenoceptor
        transformedLabel = transformedLabel.replace("-adrenoceptor", '').replace(" receptor", '');
        // Clear existing content in textElement if needed
        textElement.text('');
        // textElement.text(transformedLabel);

        // Additional transformations specific to IUPHAR labels
        const parts = transformedLabel.split(/(<sub>|<\/sub>)/); // Split label into parts including <sub> tags
        let inSub = false; // Flag to track if we are inside a subscript

        // Calculate the subscript font size (e.g., 75% of the main font size)
        const mainFontSize = parseFloat(fontSize);
        const subFontSize = parseInt(mainFontSize * 0.75);

        // Handle each part of the label
        parts.forEach(part => {
            if (part === '<sub>') {
                inSub = true;
            } else if (part === '</sub>') {
                inSub = false;
            } else {
                const tspan = textElement.append('tspan').text(part).style("font-family", "sans-serif");

                // Adjust font size and positioning for subscripts
                // tspan.attr('font-size', inSub ? subFontSize : fontSize);
                // Adjust dy based on subscript flag
                if (inSub) {
                    tspan.attr('dy', '0.5em').attr('font-size', subFontSize); // Adjust dy for subscript
                } else {
                    tspan.attr('dy', '0.1em').attr('font-size', mainFontSize); // Default dy for regular text
                }
            }
        });
    }
}

function Heatmap(data, location, heatmap_DataStyling,label_x_converter) {

    const margin = { top: 30, right: 100, bottom: 30, left: 60 }; // Adjusted margin for row labels
    const rows = Object.keys(data);
    const cols = Object.keys(data[rows[0]]);
    const col_labels = cols.map(col => label_x_converter[col]);

    const rotation = heatmap_DataStyling.rotation;
    const data_labels = heatmap_DataStyling.datalabels;
    const labelType = heatmap_DataStyling.LabelType;
    const label_position = heatmap_DataStyling.label_position;
    const label_fontsize = heatmap_DataStyling.label_fontsize;
    const receptor_fontsize = heatmap_DataStyling.receptor_fontsize;
    const data_fontsize = heatmap_DataStyling.data_fontsize;
    let legend_y_position = 30;

    let rowLabelWidth;
    if (rotation === 90 || rotation === 45) {
        // Adjust rowLabelWidth to be slightly wider than the font size height of x-axis labels
        rowLabelWidth = 20; // Adjust this value based on your font size and label requirements
    } else {
        rowLabelWidth = d3.max(col_labels, d => d.length* 40 * (label_fontsize/14)); // Default width calculation
    }

    // Calculate the longest column label length
    const longestLabel = d3.max(col_labels, d => d.length);
    if (label_position === 'Top') {
        // Calculate the longest column label length from the data
        const fontSize = label_fontsize; // Assuming a font size of 14px
        // Adjust margin calculation with a scaling factor
        if (rotation === 90) {
            margin.top = longestLabel * (fontSize * 0.6); // Adding extra space for padding
        }
        legend_y_position = 0;
    } else if (label_position === 'Bottom') {
        if (rotation === 90) {
             // Calculate the longest column label length from the data
            const fontSize = label_fontsize; // Assuming a font size of 14px
            legend_y_position = longestLabel * (fontSize * 0.55)
        }
    }

    const chartData = [];
    let highest_value = -Infinity;
    let lowest_value = Infinity;

    rows.forEach(row => {
      cols.forEach(col => {
        const value = data[row][col];
        chartData.push({ row, col, value });
        if (value > highest_value) {
          highest_value = value;
        }
        if (value < lowest_value) {
            lowest_value = value;
        }
      });
    });

    const width = (20 * cols.length) + rowLabelWidth;
    const height = (20 * rows.length);
    const svg_home = d3.select("#" + location)
      .append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + (rotation === 90 ? (margin.top * longestLabel/3) : (margin.top*2))) // Needs to account for label length or something like it.
      .attr("id", "Heatmap_plot");

    //  x / y scale transformers
    const x = d3.scale.ordinal()
        .rangeBands([0, width])
        .domain(cols);

    const y = d3.scale.ordinal()
      .rangeBands([height, 0])
      .domain(rows.reverse());

    if (heatmap_DataStyling.Number_of_colors === 'Three') {
        var myColor = d3.scale.linear()
        .range([heatmap_DataStyling.min_color, heatmap_DataStyling.middle_color, heatmap_DataStyling.max_color]) // Adjust colors here
        .domain([lowest_value, (highest_value + lowest_value) / 2, highest_value]); // Adjust domain based on your data
    } else if (heatmap_DataStyling.Number_of_colors === 'One' || heatmap_DataStyling.Number_of_colors === 'Two') {
        var myColor = d3.scale.linear()
        .range([heatmap_DataStyling.min_color, heatmap_DataStyling.max_color]) // Adjust colors here
        .domain([lowest_value, highest_value]); // Adjust domain based on your data
    }

    const svg = svg_home.append("g")
        .attr("transform", `translate(${margin.left * 2}, ${margin.top})`);

    let xAxis;
    if (label_position === 'Bottom') {
    xAxis = svg.append("g")
        .attr("transform", `translate(0, ${height})`)
        .attr('id', 'Xaxis')
        .call(d3.svg.axis().scale(x).orient("bottom").tickSize(0));
    } else if (label_position === 'Top') {
    xAxis = svg.append("g")
        .attr("transform", `translate(0, -30)`)
        .attr('id', 'Xaxis')
        .call(d3.svg.axis().scale(x).orient("top").tickSize(0));
    }

    xAxis.selectAll("text")
    .attr("class", "column-label")
    .style("font-size", `${label_fontsize}px`)
    .style("font-family", "sans-serif")
    .each(function(d, i) {
        const col = cols[i]; // Assuming cols is defined somewhere in your code
        const text = d3.select(this);

        // Set text content based on label_x_converter
        text.text(label_x_converter[col]);

        if (rotation === 0) {
            text.style("text-anchor", "middle")
                .attr("transform", "rotate(0)")
                .attr("dx", `0`)
                .attr("dy", `${label_fontsize * 1.5}px`); // Adjust dy based on font size
        } else if (rotation === 45) {
            text.style("text-anchor", label_position === 'Bottom' ? "end" : "start")
                .attr("transform", "rotate(-45)")
                .attr("dx", label_position === 'Bottom' ? `${-0.5 * label_fontsize}px` : `${-1.1 * label_fontsize}px`)
                .attr("dy", label_position === 'Bottom' ? `${0.5 * label_fontsize}px` : `${1.1 * label_fontsize}px`);
        } else if (rotation === 90) {
            text.style("text-anchor", label_position === 'Bottom' ? "end" : "start")
                .attr("transform", "rotate(-90)")
                .attr("dx", label_position === 'Bottom' ? `${-0.3 * label_fontsize}px` : `${-1.3 * label_fontsize}px`)
                .attr("dy", label_position === 'Bottom' ? `0px` : `${0.6 * label_fontsize}px`);
        }
    });

    // Y labels
    const yAxis = svg.append("g")
        .attr('id', 'Yaxis')
        .call(d3.svg.axis().scale(y).orient("left").tickSize(0))
        .selectAll("text")
        .style("font-size", `${receptor_fontsize}px`)
        .style("font-family", "sans-serif");

    yAxis.each(function(d) {
        const textElement = d3.select(this);

        // Call function to handle row labels based on label type
        handleRowLabels(textElement, d, labelType,receptor_fontsize);
    });

    // Heatmap data boxes

    const rects = svg.selectAll("rect")
      .data(chartData, d => `${d.row}:${d.col}`)
      .enter()
      .append("rect")
      .attr("x", d => x(d.col))
      .attr("y", d => y(d.row))
      .attr("width", x.rangeBand())
      .attr("height", y.rangeBand())
      .style("fill", d => myColor(d.value))
      .each(function(d) {
        if (heatmap_DataStyling.data_border) {
            d3.select(this)
                .style('stroke', 'black')
                .style('stroke-width', 0.5);
        }
    });


    if (data_labels) {
      rects.each(function(d) {
        const rect = d3.select(this);
        const textColor = getContrastColor(myColor(d.value));

        svg.append("text")
          .attr("x", +rect.attr("x") + rect.attr("width") / 2)
          .attr("y", +rect.attr("y") + rect.attr("height") / 2)
          .attr("dy", ".35em")
          .attr("text-anchor", "middle")
          .style("font-size", `${data_fontsize}px`)
          .style("font-family", "sans-serif")
          .style("fill", textColor)
          .text(d.value);
      });
    }

    function getContrastColor(hexColor) {
      const rgb = d3.rgb(hexColor);
      const brightness = (rgb.r * 0.299 + rgb.g * 0.587 + rgb.b * 0.114);
      return brightness > 128 ? "black" : "white";
    }

    // ##########
    // # Legend #
    // ##########

    // Adjust for the width of the legend to center it
    let legendX;
    if (cols.length % 2 === 0) {
        // Even number of columns
        const leftMiddleColIndex = (cols.length / 2) - 1;
        const rightMiddleColIndex = cols.length / 2;
        const leftMiddleX = x(cols[leftMiddleColIndex]) + x.rangeBand() / 2;
        const rightMiddleX = x(cols[rightMiddleColIndex]) + x.rangeBand() / 2;
        legendX = (leftMiddleX + rightMiddleX) / 2;
    } else {
        // Odd number of columns
        const middleColIndex = Math.floor(cols.length / 2);
        legendX = x(cols[middleColIndex]) + x.rangeBand() / 2; // Center of the middle column
    }

    // Adjust for the width of the legend to center it
    let legendWidth;
    if (width <= 300) {
        legendWidth = width;
    } else {
        legendWidth = 300;
    }
    const adjustedLegendX = legendX - (legendWidth / 2);

    // Position the legend using the adjustedLegendX
    const legend_svg = svg.append("g")
        .attr("transform", `translate(${adjustedLegendX}, ${height + legend_y_position})`);

    const gradientId = "Heatmap_gradient";
    const defs = svg_home.append('defs');

    // Gradient definition
    const gradient = defs.append('linearGradient')
      .attr('id', gradientId)
      .attr('x1', '0%')
      .attr('x2', '100%')
      .attr('y1', '0%')
      .attr('y2', '0%');

    gradient.append('stop')
      .attr('offset', '0%')
      .attr('stop-color', heatmap_DataStyling.min_color); // Start with red at 0%

    if (heatmap_DataStyling.Number_of_colors === 'Three') {
        gradient.append('stop')
        .attr('offset', '50%')
        .attr('stop-color', heatmap_DataStyling.middle_color); // Middle point is white
    }
    gradient.append('stop')
      .attr('offset', '100%')
      .attr('stop-color', heatmap_DataStyling.max_color); // End with green at 100%

    const gradientRect = legend_svg.append('rect')
      .attr('x', 1) // Adjust for border
      .attr('y', 16) // Adjust for border
      .attr('width', legendWidth - 2) // Adjust for border
      .attr('height', 14) // Adjust for border
      .style('fill', `url(#${gradientId})`);

    if (heatmap_DataStyling.data_border) {
      const borderRect = legend_svg.append('rect')
      .attr('x', 0)
      .attr('y', 15)
      .attr('width', legendWidth)
      .attr('height', 15)
      .style('fill', 'none')
      .style('stroke', 'black')
      .style('stroke-width', 1);
    }
    legend_svg.append("text")
      .attr('x', 0)
      .attr('y', 50)
      .style("font-size", `${data_fontsize}px`)
      .style("font-family", "sans-serif")
      .text(lowest_value);

    if (heatmap_DataStyling.Number_of_colors === 'Three') {
        legend_svg.append("text")
        .attr('x', legendWidth / 2)
        .attr('y', 50)
        .style("font-size", `${data_fontsize}px`)
        .style("font-family", "sans-serif")
        .style("text-anchor", "middle")
        .text((highest_value + lowest_value) / 2); // Middle value
    }

    legend_svg.append("text")
      .attr('x', legendWidth)
      .attr('y', 50)
      .style("font-size", `${data_fontsize}px`)
      .style("font-family", "sans-serif")
      .style("text-anchor", "end")
      .text(highest_value);

    // Rerender height of plot as the last thing
    let label_length_final;
    if (rotation == 90) {
        // label_length_final = Math.ceil(longestLabel*0.8 + 13.25*label_fontsize)
        label_length_final = Math.ceil(-128.69+6.61*longestLabel+11.19*label_fontsize)
        svg_home.attr("height",height + margin.bottom + label_length_final+55);
    } else {
        svg_home.attr("height",height + margin.bottom + 55 + 55);
    }
    
  }

// #################
// ###  Target  ###
// #################

function Target_initializeData(data) {
    // Initialize the Targets
    let Targets = {
        Target_1: {},
        Target_2: {},
        Target_3: {},
        Target_4: {},
        Target_5: {},
        Target_6: {},
        Target_7: {}
    };

    // Helper function to add items to the Target using receptor family as key
    function addItemsToCircle(Target, items) {
        Object.keys(items).forEach(ligandType => {
            const receptorFamilies = items[ligandType];
            
            // Ensure receptorFamilies is an object
            if (typeof receptorFamilies === 'object' && receptorFamilies !== null) {
                Object.keys(receptorFamilies).forEach(family => {
                    if (!Target[family]) {
                        Target[family] = [];
                    }

                    const receptors = receptorFamilies[family];
                    if (Array.isArray(receptors)) {
                        Target[family].push(...receptors);
                    } else {
                        console.warn(`Unexpected data format for family ${family}`);
                    }
                });
            } else {
                console.warn(`Unexpected data format for ligand type ${ligandType}`);
            }
        });
    }

    // Process the data
    Object.keys(data).forEach(className => {
        const classData = data[className];

        // Handle cases where classData is not an object
        if (typeof classData !== 'object' || classData === null) {
            console.warn(`Expected object but got ${typeof classData} for class ${className}`);
            return;
        }

        // Circle 1: All of "Class A (Rhodopsin)" excluding "Orphan receptors"
        if (className === "Class A (Rhodopsin)") {
            Object.keys(classData).forEach(ligandType => {
                if (ligandType !== "Orphan receptors" && ligandType !== "Olfactory receptors") {
                    addItemsToCircle(Targets.Target_1, { [ligandType]: classData[ligandType] });
                }
            });
        }

        // Circle 2: "Orphan receptors" from "Class A (Rhodopsin)" 
        if (className === "Class A (Rhodopsin)" && classData["Orphan receptors"]) {
            addItemsToCircle(Targets.Target_2, { "Orphan receptors": classData["Orphan receptors"] });
        }

        // Circle 3: "Class B1 (Secretin)" or "Class B2 (Adhesion)"
        if (className === "Class B1 (Secretin)" || className === "Class B2 (Adhesion)") {
            addItemsToCircle(Targets.Target_3, classData);
        }

        // Circle 4: "Class C (Glutamate)"
        if (className === "Class C (Glutamate)") {
            addItemsToCircle(Targets.Target_4, classData);
        }

        // Circle 5: "Class F (Frizzled)"
        if (className === "Class F (Frizzled)") {
            addItemsToCircle(Targets.Target_5, classData);
        }

        // Circle 6: "Class T (Taste 2)"
        if (className === "Class T (Taste 2)") {
            addItemsToCircle(Targets.Target_6, classData);
        }

        // Circle 7: "Other GPCRs"
        if (className === "Other GPCRs") {
            addItemsToCircle(Targets.Target_7, classData);
        }
    });

    // Convert the arrays to unique values
    Object.keys(Targets).forEach(TargetKey => {
        Object.keys(Targets[TargetKey]).forEach(familyKey => {
            Targets[TargetKey][familyKey] = Array.from(new Set(Targets[TargetKey][familyKey]));
        });
    });

    // Print out Target name and total number of values within each Target
    // Object.keys(Targets).forEach(TargetKey => {
    //     let totalValues = 0;
    //     Object.keys(Targets[TargetKey]).forEach(familyKey => {
    //         totalValues += Targets[TargetKey][familyKey].length;
    //     });
    //     console.log(`${TargetKey}: Total number of values = ${totalValues}`);
    // });

    return Targets;
}

function formatTextWithHTML(text) {
    // Apply all the replacements step by step
    return text
        .replace(" receptor", '')
        .replace("-adrenoceptor", '')
        .replace(" receptor-", '-')
        .replace("<sub>", '</tspan><tspan baseline-shift="sub">')
        .replace("</sub>", '</tspan><tspan>')
        .replace("<i>", '</tspan><tspan font-style="italic">')
        .replace("</i>", '</tspan><tspan>')
        .replace("Long-wave-sensitive", 'LWS')
        .replace("Medium-wave-sensitive", 'MWS')
        .replace("Short-wave-sensitive", 'SWS')
        .replace("Olfactory", 'OLF')
        .replace("calcitonin-like receptor", 'CLR');
}


function Draw_Targets(layout_data, fill_data, location,Target_styling) {

    dimensions = { height: 1000, width: 1000 };

    const Spacing = Target_styling.Spacing;

    // Append SVG to the div with the provided id
    const svg = d3.select("#" + location)
        .append("svg")
        .attr("width", dimensions.width)
        .attr("height", dimensions.height);

    // Draw concentric Target for layout_data
    Draw_a_Target(layout_data.Target_1, fill_data, 0, dimensions,Spacing);
    Draw_a_Target(layout_data.Target_2, fill_data, 1, dimensions,Spacing);
    Draw_a_Target(layout_data.Target_3, fill_data, 2, dimensions,Spacing);
    Draw_a_Target(layout_data.Target_4, fill_data, 3, dimensions,Spacing);
    Draw_a_Target(layout_data.Target_5, fill_data, 4, dimensions,Spacing);
    Draw_a_Target(layout_data.Target_6, fill_data, 5, dimensions,Spacing);
    Draw_a_Target(layout_data.Target_7, fill_data, 6, dimensions,Spacing);

    function Draw_a_Target(label_data, fill_data, level, dimensions,Spacing) {

        // Define SVG dimensions
        const width = dimensions.width;
        const height = dimensions.height;
        const label_offset = 7; // Increased offset to push labels outward
        

        const Target_radius = Math.min(width, height) / 2 - 60 - (70 * level); // Radius for each Target

        let values = [];
        if (Spacing && Object.keys(label_data).length > 1) {
            // If Spacing is true and there are multiple keys
            for (const key in label_data) {
                if (label_data.hasOwnProperty(key)) {
                    values = values.concat(label_data[key], ""); // Add an empty string as a placeholder
                }
            }
        } else {
            // If Spacing is false or there is only one key
            values = Object.values(label_data).flat(); // Flatten the data without placeholders
        }

        // Add in Class

        if (level == 0) {
            values.unshift("CLASS A");
        } else if (level == 1) {
            values.unshift("A ORPHANS");
        } else if (level == 2) {
            if (Spacing) {
                values.splice(19, 0, "CLASS B2"); 
            } else {
                values.splice(15, 0, "CLASS B2");
            }
            values.unshift("CLASS B1");
        } else if (level == 3) {
            values.unshift("CLASS C");
        } else if (level == 4) {
            values.unshift("CLASS F");
        } else if (level == 5) {
            values.unshift("CLASS T");
        } else if (level == 6) {
            values.unshift("CLASSLESS");
        }

        function calculatePositionAndAngle(index, total) {
            // Offset the angle calculation by -90 degrees (or -π/2 radians) to start at 12 o'clock
            const angle = -((index / total) * 2 * Math.PI) + (Math.PI / 2);
        
            // Position on the Target's border with or without label offset
            const x = width / 2 + Math.cos(angle) * (index === 0 ? (Target_radius-3) : (Target_radius + label_offset));
            const y = height / 2 - Math.sin(angle) * (index === 0 ? (Target_radius-3) : (Target_radius + label_offset));
        
            // Rotation angle so that the text faces outward
            const rotation = -(angle * 180 / Math.PI);
        
            return { x, y, rotation };
        }

        function calculateDataPositionAndAngle(index, total) {
            // Offset the angle calculation by -90 degrees (or -π/2 radians) to start at 12 o'clock ("-" its clockwise)
            const angle = -((index / total) * 2 * Math.PI) + (Math.PI / 2);


            // Position on the Target's border with label offset
            const x = width / 2 + Math.cos(angle) * (Target_radius);
            const y = height / 2 - Math.sin(angle) * (Target_radius);

            // Rotation angle so that the text faces outward
            const rotation = -(angle * 180 / Math.PI);

            return { x, y, rotation };
        }
        if (level % 2 === 0) {  // Check if it's every other circle (i.e., even index)
            // Append the outer circle (larger radius)
            svg.append("circle")
            .attr("cx", width / 2)  // X position (center of the Target)
            .attr("cy", height / 2) // Y position (center of the Target)
            .attr("r", Target_radius+60)  // Outer radius of the band
            .style("fill", "lightgray")  // Light gray fill color
            .style("opacity", 0.5)       // 50% transparency
            .style("stroke", "none");    // No stroke for the outer circle

            // Append the inner circle (smaller radius) with the same fill to cover the inner part
            svg.append("circle")
            .attr("cx", width / 2)  // X position (center of the Target)
            .attr("cy", height / 2) // Y position (center of the Target)
            .attr("r", Target_radius-10)  // Inner radius of the band
            .style("fill", "white")      // White fill to cover the inner area
            .style("opacity", 1)         // No transparency for the inner circle
            .style("stroke", "none");    // No stroke for the inner circle
        }
        // console.log(calculatePositionAndAngle(61, values.length),values.length);

        // Bind data and append text elements for the specific Target
        svg.selectAll(`.Target-text-${level}`)  // Unique selection class for each level
        .data(values)
        .enter()
        .append("text")
        .attr("class", (d, i) => 
            i === 0 ? `Target-text-${level} Target-text-${level}-first` : `Target-text-${level}`
        )  // Assign an additional class to the first element
        .attr("x", (d, i) => {
            const pos = calculatePositionAndAngle(i, values.length);
            return pos.x;
        })
        .attr("y", (d, i) => {
            const pos = calculatePositionAndAngle(i, values.length);
            return pos.y;
        })
        .attr("text-anchor", (d, i) => {
            const angle = (i / values.length) * 360 - 90;
            return (angle >= -90 && angle < 90) ? "start" : "end";
        })
        .attr("dominant-baseline", "middle")
        .attr("transform", (d, i) => {
            const pos = calculatePositionAndAngle(i, values.length);
            const angle = (i / values.length) * 360 - 90;
            const rotation = angle >= -90 && angle < 90  ? 0 : 180;
            return `rotate(${pos.rotation + rotation}, ${pos.x}, ${pos.y})`;
        })
        .html(d => formatTextWithHTML(d))
        .style("font-size", (d, i) => i === 0 ? "10px" : "9px")  // Larger font size for the first element
        .style("font-family", "Palatino")
        .attr("font-weight", (d, i) => i === 0 ? "bold" : "normal")
        .style("fill", "black");

         // Add data Targets
        svg.selectAll(`.data-Target-${level}`) // Unique selection class for each level
        .data(values) // Use values to ensure correct positions
        .enter()
        .append("circle")
        .attr("class", `data-Target-${level}`)  // Unique class to each Target element
        .attr("cx", (d, i) => {
            const pos = calculateDataPositionAndAngle(i, values.length);
            return pos.x;
        })
        .attr("cy", (d, i) => {
            const pos = calculateDataPositionAndAngle(i, values.length);
            return pos.y;
        })
        .attr("r", 3) // Adjust the radius of the Targets as needed
        .style("fill", d => {
            // Use the label to get the corresponding value in fill_data
            const value = fill_data[d]?.Value1;
            if (value === "Yes") return "black";
            if (value === "No") return "none";
            return "none"; // For empty values
        })
        .style("stroke", "black")
        .style("stroke-width", d => {
            const value = fill_data[d]?.Value1;
            return value ? 1 : 0; // Set stroke-width to 0 if there's no value
        });
    }
}







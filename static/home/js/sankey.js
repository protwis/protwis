d3v4.sankey = function() {
  var allLinks = [];  // Store all links for click event tracing
  var sankey = {},
      nodeWidth = 24,
      nodePadding = 8,
      size = [1, 1],
      nodes = [],
      links = [];

  sankey.nodeWidth = function(_) {
    if (!arguments.length) return nodeWidth;
    nodeWidth = +_;
    return sankey;
  };

  sankey.nodePadding = function(_) {
    if (!arguments.length) return nodePadding;
    nodePadding = +_;
    return sankey;
  };

  sankey.nodes = function(_) {
    if (!arguments.length) return nodes;
    nodes = _;
    return sankey;
  };

  sankey.links = function(_) {
    if (!arguments.length) return links;
      allLinks = _;  // Store all links (primary + interconnected)
      links = _.filter(function(link) {
          return link.linkage_key === "primary";  // Only use primary for layout
      });
      return sankey;
  };

  sankey.size = function(_) {
    if (!arguments.length) return size;
    size = _;
    return sankey;
  };

  sankey.layout = function(iterations) {
    computeNodeLinks();
    computeNodeValues();
    adjustNodeHeightForLabels(); // Measure text height and store as minHeight
    computeNodeBreadths();
    computeNodeDepths(iterations);
    computeLinkDepths();
    calculateLinkHeights();
    return sankey;
  };

  sankey.relayout = function() {
    computeLinkDepths();
    return sankey;
  };
  

  function calculateLinkHeights() {
    links.forEach(function(link) {
      // Calculate start and end heights for the link
      link.startHeight = (link.value / link.source.value) * link.source.dy;
      link.endHeight = (link.value / link.target.value) * link.target.dy;
    });
  }

  // Measure Text Height and Store as minHeight
  function adjustNodeHeightForLabels() {
    nodes.forEach(function(node) {
      var requiredHeight = measureTextHeight(node) + 5;
      node.minHeight = requiredHeight;  // Store it as minHeight
    });
  }

  // function to map out height of text
  function measureTextHeight(node) {
    var tempDiv = d3v4.select("body").append("div")
        .style("position", "absolute")
        .style("visibility", "hidden")
        .style("width", "140px")
        .style("font-weight", "bold")
        .style("font-size", "12px")
        .style("line-height", "1.1em")
        .style("overflow-wrap", "break-word")
        .style("white-space", "pre-line")
        .style("hyphens", "none")
        .text(node.name);

    var height = tempDiv.node().getBoundingClientRect().height;
    tempDiv.remove();
    return height;
  }

  // Populate the sourceLinks and targetLinks for each node.
  // Also, if the source and target are not objects, assume they are indices.
  function computeNodeLinks() {
    nodes.forEach(function(node) {
      node.sourceLinks = [];
      node.targetLinks = [];
    });
    links.forEach(function(link) {
      var source = link.source,
          target = link.target;
      if (typeof source === "number") source = link.source = nodes[link.source];
      if (typeof target === "number") target = link.target = nodes[link.target];
      source.sourceLinks.push(link);
      target.targetLinks.push(link);
    });
  }

  // Compute the value (size) of each node by summing the associated links.
  function computeNodeValues() {
    nodes.forEach(function(node) {
      node.value = Math.max(
        d3v4.sum(node.sourceLinks, value),
        d3v4.sum(node.targetLinks, value)
      );
    });
  }

  // Iteratively assign the breadth (x-position) for each node.
  // Nodes are assigned the maximum breadth of incoming neighbors plus one;
  // nodes with no incoming links are assigned breadth zero, while
  // nodes with no outgoing links are assigned the maximum breadth.
  function computeNodeBreadths() {
    var remainingNodes = nodes,
        nextNodes,
        x = 0;

    while (remainingNodes.length) {
      nextNodes = [];
      remainingNodes.forEach(function(node) {
        node.x = x;
        node.dx = nodeWidth;
        node.sourceLinks.forEach(function(link) {
          nextNodes.push(link.target);
        });
      });
      remainingNodes = nextNodes;
      ++x;
    }

    //
    moveSinksRight(x);
    scaleNodeBreadths((size[0] - nodeWidth) / (x - 1));
  }

  // function moveSourcesRight() {
  //   nodes.forEach(function(node) {
  //     if (!node.targetLinks.length) {
  //       node.x = d3v4.min(node.sourceLinks, function(d) { return d.target.x; }) - 1;
  //     }
  //   });
  // }

  function moveSinksRight(x) {
    nodes.forEach(function(node) {
      if (!node.sourceLinks.length) {
        node.x = x - 1;
      }
    });
  }

  function scaleNodeBreadths(kx) {
    nodes.forEach(function(node) {
      node.x *= kx;
    });
  }

  function computeNodeDepths(iterations) {
    var nodesByBreadth = d3v4.nest()
        .key(function(d) { return d.x; })
        .sortKeys(d3v4.ascending)
        .entries(nodes)
        .map(function(d) { return d.values; });

    //
    initializeNodeDepth();
    resolveCollisions();
    for (var alpha = 1; iterations > 0; --iterations) {
      relaxRightToLeft(alpha *= 0.99);
      resolveCollisions();
      relaxLeftToRight(alpha);
      resolveCollisions();
    }

    function initializeNodeDepth() {
      var ky = d3v4.min(nodesByBreadth, function(nodes) {
        return (size[1] - (nodes.length - 1) * nodePadding) / d3v4.sum(nodes, value);
      });

      nodesByBreadth.forEach(function(nodes) {
        nodes.forEach(function(node, i) {
          node.y = i;
          node.dy = node.value * ky;

          // Adjust node.dy if its less than minHeight
          if (node.minHeight) {
            node.dy = Math.max(node.dy, node.minHeight);
          }
        });
      });

      links.forEach(function(link) {
        link.dy = link.value * ky;
      });
    }

    function relaxLeftToRight(alpha) {
      nodesByBreadth.forEach(function(nodes, breadth) {
        nodes.forEach(function(node) {
          if (node.targetLinks.length) {
            var y = d3v4.sum(node.targetLinks, weightedSource) / d3v4.sum(node.targetLinks, value);
            node.y += (y - center(node)) * alpha;
          }
        });
      });

      function weightedSource(link) {
        return center(link.source) * link.value;
      }
    }

    function relaxRightToLeft(alpha) {
      nodesByBreadth.slice().reverse().forEach(function(nodes) {
        nodes.forEach(function(node) {
          if (node.sourceLinks.length) {
            var y = d3v4.sum(node.sourceLinks, weightedTarget) / d3v4.sum(node.sourceLinks, value);
            node.y += (y - center(node)) * alpha;
          }
        });
      });

      function weightedTarget(link) {
        return center(link.target) * link.value;
      }
    }

    function resolveCollisions() {
      nodesByBreadth.forEach(function(nodes) {
        var node,
            dy,
            y0 = 0,
            n = nodes.length,
            i;

        // Push any overlapping nodes down.
        nodes.sort(ascendingDepth);
        for (i = 0; i < n; ++i) {
          node = nodes[i];
          dy = y0 - node.y;
          if (dy > 0) node.y += dy;
          y0 = node.y + node.dy + nodePadding;
        }

        // If the bottommost node goes outside the bounds, push it back up.
        dy = y0 - nodePadding - size[1];
        if (dy > 0) {
          y0 = node.y -= dy;

          // Push any overlapping nodes back up.
          for (i = n - 2; i >= 0; --i) {
            node = nodes[i];
            dy = node.y + node.dy + nodePadding - y0;
            if (dy > 0) node.y -= dy;
            y0 = node.y;
          }
        }
      });
    }

    function ascendingDepth(a, b) {
      return a.y - b.y;
    }
  }

  function computeLinkDepths() {
    nodes.forEach(function(node) {
      node.sourceLinks.sort(ascendingTargetDepth);
      node.targetLinks.sort(ascendingSourceDepth);
    });
    nodes.forEach(function(node) {
      var sy = 0, ty = 0;
      node.sourceLinks.forEach(function(link) {
        link.sy = sy;
        sy += link.dy;
      });
      node.targetLinks.forEach(function(link) {
        link.ty = ty;
        ty += link.dy;
      });
    });

    function ascendingSourceDepth(a, b) {
      return a.source.y - b.source.y;
    }

    function ascendingTargetDepth(a, b) {
      return a.target.y - b.target.y;
    }
  }

  function center(node) {
    return node.y + node.dy / 2;
  }

  function value(link) {
    return link.value;
  }

  return sankey;
};


function adjustSVGHeight(sankey_data) {
  // 1) Find minY and maxY from node data
  let minY = Infinity;
  let maxY = -Infinity;

  sankey_data.nodes.forEach(node => {
    const top = node.y;
    const bottom = node.y + node.dy;

    minY = Math.min(minY, top);
    maxY = Math.max(maxY, bottom);
  });

  

  // 2) Calculate how tall the diagram needs to be (with some padding)
  // e.g., 50 px padding
  let neededHeight = maxY - minY + 50;

  // 3) Grab the original svg height from your #sankey element
  let originalHeight = +d3v4.select("#sankey").attr("height");

  // 4) If the neededHeight is bigger, we’ll expand
  if (neededHeight > originalHeight) {
    d3v4.select("#sankey")
      .attr("height", neededHeight)
      .style("height", neededHeight + "px");
  }

  // 5) Return the values so we can use them to shift if needed
  return { minY, maxY, neededHeight, originalHeight };
}

function applyYOffset(sankey_data) {
  // We call adjustSVGHeight and get the calculations
  let { minY, neededHeight, originalHeight } = adjustSVGHeight(sankey_data);

  // If minY < 20, let’s shift everything down so the topmost node is at y=20
  const desiredTopPadding = 20; 
  const offset = desiredTopPadding - minY; // e.g. if minY= -5, offset=25

  // If offset <= 0, we don’t need to move anything down
  if (offset <= 0) return;

  // Shift the entire <g> group 
  // (Assuming your top-level <g> is the only child of #sankey)
  d3v4.select("#sankey g")
    .attr("transform", function() {
      // We also might keep the original margin. If your code had something like:
      // .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
      // then we can preserve that like:
      const marginLeft = 0;  // or whatever your margin was
      const marginTop  = 0;  // or whatever your margin was
      return `translate(${marginLeft}, ${marginTop + offset})`;
    });
}



/**
 * Creates a Sankey plot based on data, adjust heights based on highest number of nodes
 * @param {jsondata} sankey_data - the data for generating the sankey plot
 * @param {integer} top_nodes    - the number of nodes in the most populated section
 * @param {integer} totalPoints  - the number of total nodes in the plot
 * @param {string} location      - location where to draw the plot
 * @param {integer} [fix_width]  - (Optional) Fixed width for the Sankey plot in pixels. Defaults to 1000 if not provided.
 */

function SankeyPlot(sankey_data, location, top_nodes, totalPoints, fix_width, pathMatrix){
  var margin = {top: 10, right: 10, bottom: 10, left: 10},
      height = (top_nodes * 60) - margin.top - margin.bottom;
  
  var containerWidth = document.getElementById(location).offsetWidth;
  var width = (typeof fix_width !== 'undefined' && fix_width !== null) 
      ? fix_width - margin.left - margin.right 
      : containerWidth - margin.left - margin.right;

  // append the svg object to the specified location
  var svg = d3v4.select('#' + location).append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
      .attr("id", "sankey")
      .append("g")
      .attr("transform",
            "translate(" + margin.left + "," + margin.top + ")");

  // Color scale used
  var color = d3v4.scaleSequential(d3v4.interpolateBlues)
      .domain([0, totalPoints - 1]); // Adjust the domain based on your data points

  // Set the sankey diagram properties
  var sankey = d3v4.sankey()
      .nodeWidth(25)
      .nodePadding(10)
      .size([width, height]);
  
  // Constructs a new Sankey generator with the default settings.
  sankey
      .nodes(sankey_data.nodes)
      .links(sankey_data.links)
      .layout(100);
  
  // Add in the links
  // Add in the links using custom polygon path
  var link = svg.append("g")
  .selectAll(".link")
  .data(sankey_data.links.filter(function(d) {
    return d.linkage_key === "primary";
    }))  // Filter primary links only for drawing
  .enter()
  .append("path")
  .attr("class", "link")
  .attr("d", function(d) {
    // Calculate points for curved polygon path
    var x0 = d.source.x + d.source.dx,
        x1 = d.target.x,
        xi = d3v4.interpolateNumber(x0, x1),
        x2 = xi(0.5),           // Control point 1 for curvature
        x3 = xi(1 - 0.5),       // Control point 2 for curvature
        y0 = d.source.y + d.sy + d.startHeight / 2,  // Use startHeight
        y1 = d.target.y + d.ty + d.endHeight / 2;    // Use endHeight

    // Calculate the four corner points for the path
    var topLeft     = [x0, y0 - d.startHeight / 2];
    var bottomLeft  = [x0, y0 + d.startHeight / 2];
    var bottomRight = [x1, y1 + d.endHeight / 2];
    var topRight    = [x1, y1 - d.endHeight / 2];

    // Create path with smooth curves and variable width
    var path = "M" + topLeft[0] + "," + topLeft[1] +
               "C" + x2 + "," + (topLeft[1]) + 
               " " + x3 + "," + (topRight[1]) + 
               " " + topRight[0] + "," + topRight[1] +

               "L" + bottomRight[0] + "," + bottomRight[1] +

               "C" + x3 + "," + (bottomRight[1]) + 
               " " + x2 + "," + (bottomLeft[1]) + 
               " " + bottomLeft[0] + "," + bottomLeft[1] +

               "Z";

    return path;
  })
  .attr("link-id", function(d) { return d.link_identifier; })
  .attr("lig", function(d) { return d.ligtrace; })
  .attr("prot", function(d) { return d.prottrace; })
  .attr("Source_name", function(d) { return d.source.name; })
  .attr("target_name", function(d) { return d.target.name; })
  .style("fill", "#969696")   // Changed to #969696
  .style("opacity", 0.1)      // Set opacity to 0.1
  .sort(function(a, b) { return b.dy - a.dy; });


  // Add in the nodes
  var node = svg.append("g")
    .selectAll(".node")
    .data(sankey_data.nodes)
    .enter().append("g")
      .attr("class", "node")
      .attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; });

  var paths = svg.selectAll("path");

  // Add the rectangles for the nodes
  node
    .append("rect")
      .attr("height", function(d) { return d.dy; })
      .attr("width", sankey.nodeWidth())
      .style("fill", function(d, i) {
         return color(i); // Apply color based on index
       })
       .on("click", function(d) {
        var isActive = d3v4.select(this).classed("active");
        
        // Reset all paths to default state
        paths.style("opacity", 0.1).style("fill", '#969696');
        var linkIdsToHighlight = new Set();
    
        // If the node is active and clicked again, reset all highlights
        if (isActive) {
            paths.style("opacity", 0.1).style("fill", function(d) {
                return '#969696'; // Reset to default color
            });
        } else {
            // Identify Node Type and Filter Matrix
            if (d.column === "x1") {
                // x1 Node
                var filteredRows = filterMatrixByNodeName(d.name, 'x1', pathMatrix);
                linkIdsToHighlight = collectLinkIdentifiersFromMatrix(filteredRows);
            } else if (d.column === "x2") {
                // x2 Node
                var filteredRows = filterMatrixByNodeName(d.name, 'x2', pathMatrix);
                linkIdsToHighlight = collectLinkIdentifiersFromMatrix(filteredRows);
            } else if (d.column === "x3") {
                // x3 Node
                var filteredRows = filterMatrixByNodeName(d.name, 'x3', pathMatrix);
                linkIdsToHighlight = collectLinkIdentifiersFromMatrix(filteredRows);
            } else if (d.column === "x4") {
                // x4 Node
                var filteredRows = filterMatrixByNodeName(d.name, 'x4', pathMatrix);
                linkIdsToHighlight = collectLinkIdentifiersFromMatrix(filteredRows);
            }
    
            console.log("Link IDs to Highlight:", Array.from(linkIdsToHighlight));  // Debugging
    
            // Highlight Only the Relevant Paths
            paths
                .style("opacity", 0) // Make all paths invisible first
                .filter(function(d) {
                    return linkIdsToHighlight.has(d.link_identifier);
                })
                .style("opacity", 0.8)
                .style("fill", "rgb(214, 230, 244)");
        }
    
        // Toggle the active state of the node
        d3v4.select(this).classed("active", !isActive);
    })

      function filterMatrixByNodeName(nodeName, column, pathMatrix) {
        return pathMatrix.filter(function(row) {
            return row[column + '_name'] === nodeName;
        });
    }
    
    
    function collectLinkIdentifiersFromMatrix(pathMatrix) {
      var linkIdsToHighlight = new Set();
  
      // Loop through each row in the pathMatrix
      pathMatrix.forEach(function(row) {
          // Go through each stage of the path
          var stages = [
              [row.x1_name, row.x2_name],  // x1 -> x2
              [row.x2_name, row.x3_name],  // x2 -> x3
              [row.x3_name, row.x4_name]   // x3 -> x4
          ];
  
          // For each stage, find matching links
          stages.forEach(function(stage) {
              var sourceName = stage[0];
              var targetName = stage[1];
  
              // Find matching links by comparing Source_name and target_name
              d3v4.selectAll(".link")
                  .filter(function(d) {
                      return d.source.name === sourceName && d.target.name === targetName;
                  })
                  .each(function(d) {
                      linkIdsToHighlight.add(d.link_identifier);
                  });
          });
      });
  
      return linkIdsToHighlight;
  }
    
  
  node.append("foreignObject")
  .attr("x", function(d) {
    return d.x < width / 2 ? sankey.nodeWidth() + 6 : -150;
  })
  .attr("y", 0)
  .attr("width", 140)
  .attr("height", function(d) { return d.dy; })
  .style("overflow", "visible")
  .append("xhtml:div")
  .style("display", "flex")
  .style("flex-direction", "column")
  .style("justify-content", "center")
  .style("align-items", function(d) {
    return d.x < width / 2 ? "flex-start" : "flex-end";
  })
  .style("text-align", function(d) {
    return d.x < width / 2 ? "left" : "right";
  })
  .style("width", "100%")
  .style("height", "100%")
  .style("line-height", "normal")
  .style("font-weight", "bold")
  .style("font-size", "12px")
  .style("white-space", "pre-wrap")
  .style("hyphens", "none")
  .html(function(d) {
    return `<span style="cursor: pointer; color: black; text-decoration: none;" onmouseover="this.style.color='#007bff'; this.style.textDecoration='underline';" onmouseout="this.style.color='black'; this.style.textDecoration='none';" onclick="window.open('${d.url}', '_blank')">${d.name}</span>`;
  })

  applyYOffset(sankey_data); // This will call adjustSVGHeight internally, expand the SVG if needed, and then shift if needed.
}

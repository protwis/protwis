d3v4.sankey = function() {
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
    links = _;
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
    computeNodeBreadths();
    computeNodeDepths(iterations);
    computeLinkDepths();
    return sankey;
  };

  sankey.relayout = function() {
    computeLinkDepths();
    return sankey;
  };

  sankey.link = function() {
    var curvature = .5;

    function link(d) {
      var x0 = d.source.x + d.source.dx,
          x1 = d.target.x,
          xi = d3v4.interpolateNumber(x0, x1),
          x2 = xi(curvature),
          x3 = xi(1 - curvature),
          y0 = d.source.y + d.sy + d.dy / 2,
          y1 = d.target.y + d.ty + d.dy / 2;
      return "M" + x0 + "," + y0
           + "C" + x2 + "," + y0
           + " " + x3 + "," + y1
           + " " + x1 + "," + y1;
    }

    link.curvature = function(_) {
      if (!arguments.length) return curvature;
      curvature = +_;
      return link;
    };

    return link;
  };

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

  function moveSourcesRight() {
    nodes.forEach(function(node) {
      if (!node.targetLinks.length) {
        node.x = d3v4.min(node.sourceLinks, function(d) { return d.target.x; }) - 1;
      }
    });
  }

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
      relaxRightToLeft(alpha *= .99);
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




/**
 * Creates a Sankey plot based on data, adjust heights based on highest number of nodes
 * @param {jsondata} sankey_data - the data for generating the sankey plot
 * @param {integer} top_nodes    - the number of nodes in the most populated section
 * @param {integer} totalPoints  - the number of total nodes in the plot
 * @param {string} location      - location where to draw the plot
 */

function SankeyPlot(sankey_data, location, top_nodes, totalPoints){
  var margin = {top: 10, right: 10, bottom: 10, left: 10},
      width = 1000 - margin.left - margin.right,
      height = (top_nodes*30) - margin.top - margin.bottom;

  // append the svg object to the body of the page
  var svg = d3v4.select('#'+location).append("svg")
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
      .nodePadding(5)
      .size([width, height]);

  // Constructs a new Sankey generator with the default settings.
  sankey
      .nodes(sankey_data.nodes)
      .links(sankey_data.links)
      .layout(100);

  // add in the links
  var link = svg.append("g")
    .selectAll(".link")
    .data(sankey_data.links)
    .enter()
    .append("path")
      .attr("class", "link")
      .attr("d", sankey.link())
      .attr("lig", function(d) { return d.ligtrace; })
      .attr("prot", function(d) { return d.prottrace; })
      .style("stroke-width", function(d) { return Math.max(1, d.dy); })
      .sort(function(a, b) { return b.dy - a.dy; });

  // add in the nodes
  var node = svg.append("g")
    .selectAll(".node")
    .data(sankey_data.nodes)
    .enter().append("g")
      .attr("class", "node")
      .attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; })
      // .call(d3v4.drag()
      //   .subject(function(d) { return d; })
      //   .on("start", function() { this.parentNode.appendChild(this); })
      //   .on("drag", dragmove));

  var paths = svg.selectAll("path");
  // add the rectangles for the nodes
  node
    .append("rect")
      .attr("height", function(d) { return d.dy; })
      .attr("width", sankey.nodeWidth())
      .style("fill", function(d, i) {
         return color(i); // Apply color based on index
       })
      .on("click", function(d,i) {
         var isActive = d3v4.select(this).classed("active");
         var currentFillColor = d3v4.select(this).style("fill");
         // Reset all paths to default opacity
         paths.style("stroke-opacity", 0.1);
         paths.style("stroke", '#969696');
         if (!isActive) {
             paths.style("stroke-opacity", 0);
             // Highlight the paths with the same name
             paths.filter(function() {
                 var ligPathName = d3v4.select(this).attr("lig");
                 var protPathName = d3v4.select(this).attr("prot");
                 return ligPathName === d.name || protPathName === d.name;
             })
             .style("stroke-opacity", 0.5)
             // .style("stroke", currentFillColor);
             .style("stroke", "rgb(214, 230, 244)");
         }
         // Toggle the 'active' class on the clicked node
         d3v4.select(this).classed("active", !isActive);
       })
    // Add hover text
    .append("title")
      .text(function(d) { return d.name + "\n" + "There is " + d.value + " stuff in this node"; });


  // add in the title for the nodes
  node.append("text")
      .attr("x", -6)
      .attr("y", function(d) { return d.dy / 2; })
      .attr("dy", ".35em")
      .attr("text-anchor", "end")
      .attr("transform", null)
      .style("font-weight", "bold")
      .on("click", function(d) {window.open(d.url, '_blank')})
      .each(function(d) {
          // Replace HTML entities with Unicode characters
          var formattedName = d.name.replace(/&alpha;/g, "α").replace(/&beta;/g, "β")
                                    .replace(/&kappa;/g, "κ").replace(/&delta;/g, "δ")
                                    .replace(/&mu;/g, "µ").replace(/&gamma;/g, "γ");

          var parts = formattedName.split("<sub>"),
              mainText = parts[0],
              subText = parts[1] ? parts[1].split("</sub>")[0] : '';

          var supParts = mainText.split("<sup>"),
              supText = supParts[1] ? supParts[1].split("</sup>")[0] : '';
              mainText = supParts[0]; // Update mainText to exclude supText

          var italicParts = mainText.split("<i>"),
              italicText = italicParts[1] ? italicParts[1].split("</i>")[0] : '';
              mainText = italicParts[0]; // Update mainText to exclude italicText

          var text = d3v4.select(this);
          text.text(mainText); // Set the main text

          if(italicText) {
              text.append("tspan")
                  .style("font-style", "italic")
                  .text(italicText);

              // Check for remaining text after italic tag
              if(italicParts[1] && italicParts[1].split("</i>")[1]) {
                  text.append("tspan")
                      .style("font-style", "normal") // Reset style to normal for following text
                      .text(italicParts[1].split("</i>")[1]);
              }
          }

          if(supText) {
              text.append("tspan") // Add the superscript part
                  .attr("dy", "-0.5em") // Move the superscript up
                  .attr("dx", "-0.1em") // Nudge back to align after main text
                  .style("font-size", "smaller") // Make superscript smaller
                  .text(supText);
          }

          if(subText) {
              text.append("tspan") // Add the subscript part
                  .attr("dy", "0.7em") // Move the subscript down
                  .attr("dx", "-0.1em") // Move the subscript back to align correctly after the main text
                  .style("font-size", "smaller") // Make subscript smaller
                  .text(subText);
          }

          if(parts[1] && parts[1].split("</sub>")[1]) {
              // Add remaining text after the subscript, if any
              text.append("tspan")
                  .attr("dy", "-0.7em") // Move back up to align with the main text baseline
                  .attr("dx", "0.0em") // Move forward past the subscript
                  .text(parts[1].split("</sub>")[1]);
          }
      })
      .filter(function(d) { return d.x < width / 2; })
          .attr("x", 6 + sankey.nodeWidth())
          .attr("text-anchor", "start");


  // the function for moving the nodes
  function dragmove(d) {
    d3v4.select(this)
      .attr("transform",
            "translate("
               + d.x + ","
               + (d.y = Math.max(
                  0, Math.min(height - d.dy, d3v4.event.y))
                 ) + ")");
    sankey.relayout();
    link.attr("d", sankey.link() );
  }
}

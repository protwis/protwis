function RadarChart(id, data, options, name) {
  var cfg = {
    w: 200, //Width of the circle
    h: 200, //Height of the circle
    margin: {
      top: 50,
      right: 50,
      bottom: 50,
      left: 50
    }, //The margins of the SVG
    levels: 4, //How many levels or inner circles should there be drawn
    maxValue: 0, //What is the value that the biggest circle will represent
    minValue: 0, //What is the value that the first circle will represent
    labelFactor: 1.05, //How much farther than the radius of the outer circle should the labels be placed
    wrapWidth: 60, //The number of pixels after which a label needs to be given a new line
    opacityArea: 0.35, //The opacity of the area of the blob
    dotRadius: 4, //The size of the colored circles of each blog
    opacityCircles: 0.1, //The opacity of the circles of each blob
    strokeWidth: 2, //The width of the stroke around each blob
    roundStrokes: false, //If true the area and stroke will follow a round path (cardinal-closed)
    title: "Title", //Set the title name of the plot
    filling: "red" //Color function
  };

  //Put all of the options into a variable called cfg
  if ('undefined' !== typeof options) {
    for (var i in options) {
      if ('undefined' !== typeof options[i]) {
        cfg[i] = options[i];
      }
    } //for i
  } //if

  //create the div to host the plot
  var parentDiv = document.getElementById(id);
  var nestedDiv = document.createElement("div");
  nestedDiv.setAttribute("id", name);
  nestedDiv.setAttribute("style", "display: inline-block");
  parentDiv.appendChild(nestedDiv);
  //If the supplied maxValue is smaller than the actual one, replace by the max in the data
  var maxValue = Math.max(cfg.maxValue, d3.max(data, function(i) {
    return d3.max(i.map(function(o) {
      return o.value;
    }))
  }));
  var minValue = Math.floor(Math.min(cfg.minValue, d3.min(data, function(i) {
    return d3.min(i.map(function(o) {
      return o.value;
    }))
  })));

  cfg.levels = Math.abs(maxValue) + Math.abs(minValue);

  var check = Math.abs(minValue);

  var allAxis = (data[0].map(function(i, j) {
      return i.axis
    })), //Names of each axis
    total = allAxis.length, //The number of different axes
    radius = Math.min(cfg.w / 2, cfg.h / 2), //Radius of the outermost circle
    Format = d3.format('%'), //Percentage formatting
    angleSlice = Math.PI * 2 / total; //The width in radians of each "slice"

  //Scale for the radius
  var rScale = d3.scale.linear()
    .range([0, radius])
    .domain([minValue, maxValue]);

  /////////////////////////////////////////////////////////
  //////////// Create the container SVG and g /////////////
  /////////////////////////////////////////////////////////

  //Remove whatever chart with the same id/class was present before
  d3.select(id).html("");
  //Initiate the radar chart SVG
  var svg = d3.select('#' + name).append("svg")
    .attr("width", cfg.w + cfg.margin.left + cfg.margin.right)
    .attr("height", cfg.h + cfg.margin.top + cfg.margin.bottom)
    .attr("class", "radar_" + name);
  //Append a g element
  var g = svg.append("g")
    .attr("transform", "translate(" + (cfg.w / 2 + cfg.margin.left) + "," + (cfg.h / 2 + cfg.margin.top) + ")");

  /////////////////////////////////////////////////////////
  ////////// Glow filter for some extra pizzazz ///////////
  /////////////////////////////////////////////////////////

  //Filter for the outside glow
  var filter = g.append('defs').append('filter').attr('id', 'glow'),
    feGaussianBlur = filter.append('feGaussianBlur').attr('stdDeviation', '2.5').attr('result', 'coloredBlur'),
    feMerge = filter.append('feMerge'),
    feMergeNode_1 = feMerge.append('feMergeNode').attr('in', 'coloredBlur'),
    feMergeNode_2 = feMerge.append('feMergeNode').attr('in', 'SourceGraphic');

  /////////////////////////////////////////////////////////
  /////////////// Draw the Circular grid //////////////////
  /////////////////////////////////////////////////////////

  //Wrapper for the grid & axes
  var axisGrid = g.append("g").attr("class", "axisWrapper");

  //Draw the background circles
  axisGrid.selectAll(".levels")
    .data(d3.range(1, (cfg.levels)).reverse())
    .enter()
    .append("circle")
    .attr("class", "gridCircle")
    .attr("r", function(d, i) {
      return radius / cfg.levels * d;
    })
    .style("fill", "white")
    .style("stroke", function(d) {
      if (d == check) {
        return "gray"
      } else {
        return "lightgray"
      }
    })
    .style("fill-opacity", cfg.opacityCircles)
    .style("filter", "url(#glow)");

  //Text indicating at what % each level is
  axisGrid.selectAll(".axisLabel")
    .data(d3.range(minValue, maxValue))
    .enter().append("text")
    .attr("class", "axisLabel")
    .attr("x", 4)
    .attr("y", function(d, i) {
      return -i * radius / cfg.levels;
    })
    .attr("dy", "0.4em")
    .style("font-size", "10px")
    .attr("fill", "black")
    .text(function(d, i) {
      return (d);
    }); //ToDo: need to adjust this function

  //Creating the title label
  axisGrid
    .append("text")
    .attr("class", "title")
    .attr("x", 4)
    .attr("y", -140)
    .attr("dy", "0.35em")
    .attr("text-anchor", "middle")
    .style("font-size", "14px")
    .style("font-weight", "bold")
    .attr("fill", "black")
    .text(cfg.title.split(".,")[0]);

  /////////////////////////////////////////////////////////
  //////////////////// Draw the axes //////////////////////
  /////////////////////////////////////////////////////////

  //Create the straight lines radiating outward from the center
  var axis = axisGrid.selectAll(".axis")
    .data(allAxis)
    .enter()
    .append("g")
    .attr("class", "axis");
  //Append the lines
  axis.append("line")
    .attr("x1", 0)
    .attr("y1", 0)
    .attr("x2", function(d, i) {
      return rScale(maxValue * 1) * Math.cos(angleSlice * i - Math.PI / 2);
    })
    .attr("y2", function(d, i) {
      return rScale(maxValue * 0.9) * Math.sin(angleSlice * i - Math.PI / 2);
    })
    .attr("class", "line")
    .style("stroke", "lightgrey")
    .style("stroke-width", "1px");

  //Append the labels at each axis
  axis.append("text")
    .attr("class", "legend")
    .style("font-size", "11px")
    .attr("text-anchor", "middle")
    .attr("dy", "0.35em")
    .attr("x", function(d, i) {
      return rScale(maxValue * cfg.labelFactor) * Math.cos(angleSlice * i - Math.PI / 2);
    })
    .attr("y", function(d, i) {
      return rScale(maxValue * cfg.labelFactor) * Math.sin(angleSlice * i - Math.PI / 2);
    })
    .text(function(d) {
      return d
    })
    .call(wrap, cfg.wrapWidth);

  /////////////////////////////////////////////////////////
  ///////////// Draw the radar chart blobs ////////////////
  /////////////////////////////////////////////////////////

  //The radial line function
  var radarLine = d3.svg.line.radial()
    .interpolate("linear-closed")
    .radius(function(d) {
      return rScale(d.value);
    })
    .angle(function(d, i) {
      return i * angleSlice;
    });

  if (cfg.roundStrokes) {
    radarLine.interpolate("cardinal-closed");
  }

  //Create a wrapper for the blobs
  var blobWrapper = g.selectAll(".radarWrapper")
    .data(data)
    .enter().append("g")
    .attr("class", "radarWrapper");

  //Append the backgrounds
  blobWrapper
    .append("path")
    .attr("class", "radarArea")
    .attr("d", function(d, i) {
      return radarLine(d);
    })
    .style("fill", function(d) {
      return cfg.filling;
    })
    .style("fill-opacity", cfg.opacityArea)
    .on('mouseover', function(d, i) {
      //Dim all blobs
      d3.selectAll(".radarArea")
        .transition().duration(200)
        .style("fill-opacity", 0.1);
      //Bring back the hovered over blob
      d3.select(this)
        .transition().duration(200)
        .style("fill-opacity", 0.7);
    })
    .on('mouseout', function() {
      //Bring back all blobs
      d3.selectAll(".radarArea")
        .transition().duration(200)
        .style("fill-opacity", cfg.opacityArea);
    });

  //Create the outlines
  blobWrapper.append("path")
    .attr("class", "radarStroke")
    .attr("d", function(d, i) {
      return radarLine(d);
    })
    .style("stroke-width", cfg.strokeWidth + "px")
    .style("stroke", function(d) {
      return cfg.filling;
    })
    .style("fill", "none")
    .style("filter", "url(#glow)");

  //Append the circles
  blobWrapper.selectAll(".radarCircle")
    .data(function(d, i) {
      return d;
    })
    .enter().append("circle")
    .attr("class", "radarCircle")
    .attr("r", cfg.dotRadius)
    .attr("cx", function(d, i) {
      return rScale(d.value) * Math.cos(angleSlice * i - Math.PI / 2);
    })
    .attr("cy", function(d, i) {
      return rScale(d.value) * Math.sin(angleSlice * i - Math.PI / 2);
    })
    .style("fill", function(d, i) {
      return cfg.filling;
    })
    .style("fill-opacity", 0.8);

  /////////////////////////////////////////////////////////
  //////// Append invisible circles for tooltip ///////////
  /////////////////////////////////////////////////////////

  //Wrapper for the invisible circles on top
  var blobCircleWrapper = g.selectAll(".radarCircleWrapper")
    .data(data)
    .enter().append("g")
    .attr("class", "radarCircleWrapper");

  //Append a set of invisible circles on top for the mouseover pop-up
  blobCircleWrapper.selectAll(".radarInvisibleCircle")
    .data(function(d, i) {
      return d;
    })
    .enter().append("circle")
    .attr("class", "radarInvisibleCircle")
    .attr("r", cfg.dotRadius * 1.5)
    .attr("cx", function(d, i) {
      return rScale(d.value) * Math.cos(angleSlice * i - Math.PI / 2);
    })
    .attr("cy", function(d, i) {
      return rScale(d.value) * Math.sin(angleSlice * i - Math.PI / 2);
    })
    .style("fill", "none")
    .style("pointer-events", "all")
    .on("mouseover", function(d, i) {
      newX = parseFloat(d3.select(this).attr('cx')) - 10;
      newY = parseFloat(d3.select(this).attr('cy')) - 10;

      tooltip
        .attr('x', newX)
        .attr('y', newY)
        .text(d.value)
        .transition().duration(200)
        .style('opacity', 1);
    })
    .on("mouseout", function() {
      tooltip.transition().duration(200)
        .style("opacity", 0);
    });

  //Set up the small tooltip for when you hover over a circle
  var tooltip = g.append("text")
    .attr("class", "tooltip")
    .style("opacity", 0);

  /////////////////////////////////////////////////////////
  /////////////////// Helper Function /////////////////////
  /////////////////////////////////////////////////////////

  //Taken from http://bl.ocks.org/mbostock/7555321
  //Wraps SVG text
  function wrap(text, width) {
    text.each(function() {
      var text = d3.select(this),
        words = text.text().split(/\s+/).reverse(),
        word,
        line = [],
        lineNumber = 0,
        lineHeight = 1.4, // ems
        y = text.attr("y"),
        x = text.attr("x"),
        dy = parseFloat(text.attr("dy")),
        tspan = text.text(null).append("tspan").attr("x", x).attr("y", y).attr("dy", dy + "em");

      while (word = words.pop()) {
        line.push(word);
        tspan.text(line.join(" "));
        if (tspan.node().getComputedTextLength() > width) {
          line.pop();
          tspan.text(line.join(" "));
          line = [word];
          tspan = text.append("tspan").attr("x", x).attr("y", y).attr("dy", ++lineNumber * lineHeight + dy + "em").text(word);
        }
      }
    });
  } //wrap

} //RadarChart

function ShowRadarPlot(type) {
  var svgNS = "http://www.w3.org/2000/svg";
  var parentDiv = document.getElementById("radarcontainer_" + type);
  var mergedDiv = document.getElementById("mergedcontainer_" + type);
  parentDiv.querySelectorAll('*').forEach(n => n.remove());
  mergedDiv.querySelectorAll('*').forEach(n => n.remove());
  var title = document.createElement("div");
  title.setAttribute("class", "RadarTitle");
  parentDiv.appendChild(title);
  var ids = new Array();
  if (type === 'ligand') {
    // $('#RadarPublicationsList').val('');
    var lig = document.getElementById('RadarLigandList').value;
    var pub_hash = document.querySelector("#select-ligand option[value='" + lig + "']").dataset.value;
    var header = document.createTextNode(lig);
    title.appendChild(header);
    for (var pub in spiders) {
      name = pub.replace(/\[|\]|\(|\)|\s|\,|\'|\./g, "")
      for (var web in spiders[pub]) {
        if (web === pub_hash) {
          spiders[pub][web]["Options"]["title"] = pub;
          spiders[pub][web]["Options"]["filling"] = lig_color;
          RadarChart('radarcontainer_' + type, spiders[pub][web]["Data"], spiders[pub][web]["Options"], "radar" + name);
        }
      }
    }
  } else {
    // $('#RadarLigandList').val('');
    var pub = document.getElementById('RadarPublicationsList').value;
    var header = document.createTextNode(pub);
    title.appendChild(header);
    for (var web in spiders[pub]) {
      spiders[pub][web]["Options"]["filling"] = pub_color;
      RadarChart('radarcontainer_' + type, spiders[pub][web]["Data"], spiders[pub][web]["Options"], "radar" + web);
    }
  }

  // create a merged-div where we are going to merge the svgs
  var merged = document.createElement('div');
  merged.setAttribute('id', 'merged-div_' + type);
  mergedDiv.appendChild(merged);

  $('#radarcontainer_' + type + ' > div').map(function() {
    ids.push(this.id);
  });
  var mergedHeight = 300 * Math.round(ids.length / 3);
  if (ids.length > 3) {
    var mergedWidth = 900;
  } else {
    var mergedWidth = 300 * (ids.length - 1);
  }

  // createElementNS for svg
  var mergedSvg = document.createElementNS(svgNS, 'svg');
  mergedSvg.setAttribute('id', 'merged_' + type);
  // keep the viewBox of the chart
  mergedSvg.setAttribute('width', mergedWidth);
  mergedSvg.setAttribute('height', mergedHeight);
  merged.appendChild(mergedSvg);

  var xMultiply = -1;
  var yMultiply = 1;
  for (let i = 1; i < ids.length; i++) {
    if (xMultiply < 2) {
      xMultiply = xMultiply + 1;
    } else {
      xMultiply = 0;
      yMultiply = yMultiply + 1;
    }
    var plot = document.getElementById(ids[i]);
    var plotSvg = plot.getElementsByTagName('svg')[0];
    var plotContent = Array.from(plotSvg.childNodes)[0];
    var xTranslate = (150 * (xMultiply + 1)) + (130 * xMultiply);
    var yTranslate = (150 * yMultiply) + (150 * (yMultiply - 1));

    plotContent.setAttribute('transform', 'translate(' + xTranslate + ',' + yTranslate + ')');

    mergedSvg.appendChild(plotContent);
  }
  gNodes = d3.selectAll('.axisWrapper').selectAll('.axis');
  gNodes.each(function(d) {
    old = d3.select(this).select('text')[0][0].innerHTML;
    old = old.replace("G12/13", 'G</tspan><tspan baseline-shift = "sub">12/13</tspan><tspan>');
    old = old.replace("Gi/o", 'G</tspan><tspan baseline-shift = "sub">i/o</tspan><tspan>');
    old = old.replace("Gq/11", 'G</tspan><tspan baseline-shift = "sub">q/11</tspan><tspan>');
    old = old.replace("Gs", 'G</tspan><tspan baseline-shift = "sub">s</tspan><tspan>');
    node = d3.select(this).select('text')[0][0];
    node.innerHTML = old;
  });
}


function DrawMultiLineChart(Data, BaseDiv, Keys, ID, linkTitle, reference, linkPath, ylabel) {

    var parentDiv = document.getElementById(BaseDiv)

    var title = document.createElement("div");
        title.setAttribute("id", "title_"+ID);
        title.setAttribute("style", "padding-top: 15px;");

    var a = document.createElement("a");
        a.setAttribute("target", "_blank");

    var linkText = document.createTextNode(linkTitle);
        a.appendChild(linkText);
        a.innerHTML = linkTitle;
        a.href = linkPath;
        title.appendChild(a);

    var nestedDiv = document.createElement("div");
      nestedDiv.setAttribute("id", ID);
      parentDiv.appendChild(nestedDiv);

    var downloadDiv = document.createElement("div");
      downloadDiv.setAttribute("id", "Download_"+ID);
      downloadDiv.setAttribute("class", "btn-group");
      parentDiv.appendChild(downloadDiv);

    var button = document.createElement("BUTTON");
        button.setAttribute("class", "btn btn-primary dropdown-toggle");
        button.setAttribute("aria-haspopup", "true");
        button.setAttribute("aria-expanded", "false");
        button.setAttribute("data-toggle", "dropdown");
        button.innerHTML = "Download";
        downloadDiv.appendChild(button);

    var optionList = document.createElement("ul");
        optionList.setAttribute("class", "dropdown-menu");

    var li_png = document.createElement("li");
    var text_png = document.createElement("a");
        text_png.textContent = "PNG";
        text_png.setAttribute("href", "javascript:saveSvgAsPng(document.getElementById('"+ID+"').getElementsByTagName('svg')[0], '"+ reference + "_" + ID +".png');");
        li_png.appendChild(text_png);
        optionList.appendChild(li_png);

    var li_jpg = document.createElement("li");
    var text_jpg = document.createElement("a");
        text_jpg.textContent = "JPG";
        text_jpg.setAttribute("href", "javascript:saveSvgAsJpg(document.getElementById('"+ID+"').getElementsByTagName('svg')[0], '"+ reference + "_" + ID +".jpg');");
        li_jpg.appendChild(text_jpg);
        optionList.appendChild(li_jpg);

    var li_tiff = document.createElement("li");
    var text_tiff = document.createElement("a");
        text_tiff.textContent = "TIFF";
        text_tiff.setAttribute("href", "javascript:saveSvgAsTiff(document.getElementById('"+ID+"').getElementsByTagName('svg')[0], '"+ reference + "_" + ID +".tiff');");
        li_tiff.appendChild(text_tiff);
        optionList.appendChild(li_tiff);

    var li_svg = document.createElement("li");
    var text_svg = document.createElement("a");
        text_svg.textContent = "SVG";
        text_svg.setAttribute("href", "javascript:saveSvg(document.getElementById('"+ID+"'), '"+ reference + "_" + ID +".svg');");
        li_svg.appendChild(text_svg);
        optionList.appendChild(li_svg);

        downloadDiv.appendChild(optionList);

    var margin = { top: 20, right: 80, bottom: 30, left: 50 },
     width = 1120 - margin.left - margin.right,
     height = 350 - margin.top - margin.bottom;

    function ResetOpacity(){
       d3.selectAll("circle")
         .style("opacity", 1);
       d3.selectAll("rect")
         .style("opacity", 1);
       d3.selectAll("g.segment path")
         .style("opacity", 1);
     };

    document.getElementById(ID).onclick = ResetOpacity;

    var x = d3.scale.ordinal()
            .rangeRoundBands([0, width]);

    var y = d3.scale.linear()
        .range([height, 0])

    var color = d3.scale.category20();
    var legendColor = d3.scale.category20();

    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("bottom");

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left")
        .ticks(10);

    // xData gives an array of distinct "Pathways" for which trends chart is going to be made.
    var xData = Data[0].PathwaysData.map(function (d) { return d.Pathway; });

    var real_points = Data[0].PathwaysData.filter(function(d) {
      return d.value[1] == "REAL"; });

    var artificial_points = Data[0].PathwaysData.filter(function(d) {
      return d.value[1] == "ARTIFICIAL"; });

    var line = d3.svg.line()
        //.interpolate("basis")
        .x(function (d) { return (x(d.Pathway) + x.rangeBand() / 2)-100; })
        .y(function (d) { return y(d.value[0]); });

    // document.BaseDiv.appendChild(div);

    var svg = d3.select("#" + ID)
        .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom + 30)
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + (margin.top + 30) + ")");

    color.domain(Data.map(function (d) { return d.name; }));

    x.domain(xData);

    var valueMax = d3.max(Data, function (r) { return d3.max(r.PathwaysData, function (d) { return d.value[0]; }) });
    var valueMin = d3.min(Data, function (r) { return d3.min(r.PathwaysData, function (d) { return d.value[0]; }) });
    y.domain([valueMin, valueMax]);

    //Drawing title label
    svg.append("g")
            // .attr("transform", "translate(- " + margin.top + ",- " + margin.top + ")")
            // .append("foreignObject")
            //   .attr("width", width)
            //   .attr("height", 40)
            //   .attr("class", "title")
            .append("text")
              .attr("y", -20)
              .attr("x", 30)
              .style("fill", "#357db5")
              .style("font-size", "15px")
              .attr("font-family", "Arial")
              .style("font-style", "italic")
              .style("padding-bottom", "3px")
              .style("padding-top", "15px")
              .text(linkTitle)
              .on("click", function(){
                  window.open(a, '_blank').focus();
              });


    //Drawing X Axis
    svg.append("g")
            .attr("class", "x axis")
            .style("font", "12px sans-serif")
            .attr("transform", "translate(-100," + height + ")")
            .call(xAxis);

    // Drawing Horizontal grid lines.
    svg.append("g")
        .attr("class", "GridX")
      .selectAll("line.grid").data(y.ticks()).enter()
        .append("line")
        .attr(
        {
            "class": "grid",
            "x1": x(xData[0]),
            "x2": x(xData[xData.length - 1]) + x.rangeBand() / 2,
            "y1": function (d) { return y(d); },
            "y2": function (d) { return y(d); }
        });
    // Drawing Y Axis
    svg.append("g")
        .attr("class", "y axis")
        .style("font", "10px sans-serif")
        .call(yAxis)
        .append("text")
            .attr("transform", "rotate(-90)")
            .attr("y", -35)
            .attr("x", -50)
            .attr("dy", ".10em")
            .style("font", "12px sans-serif")
            .style("text-anchor", "end")
            .text(ylabel + " relative to " + reference);

    // Drawing Lines for each segments
    var segment = svg.selectAll(".segment")
                    .data(Data)
                    .enter().append("g")
                    .attr("class", function(d) { return "segment "+ d.name;});
                    // .attr("id", function (d) { return d.name; });

    segment.append("path")
            .attr("class", "line")
            .attr("visible",1)
            .attr("d", function (d) { return line(d.PathwaysData); })
            .style("stroke", function (d) { return color(d.name); });

    // Creating Dots on line
    segment.selectAll("dot")
            .data(function (d) { return d.PathwaysData; } )
            .enter()
            .append("circle")
            .filter(function(d) { return d.value[1] == "REAL"})
            .attr("class",function(d) { return d.value[1] })
            .attr("r", 5)
            .attr("cx", function (d) { return (x(d.Pathway) + x.rangeBand() / 2)-100; })
            .attr("cy", function (d) { return y(d.value[0]); })
            .style("stroke", "black")
            .style("fill", function (d) { return color(this.parentNode.__data__.name); })
            .on("mouseover", mouseover)
            .on("mousemove", function (d) {
                divToolTip
                // .text(d.value[0].toFixed(2))
                .html(d.tooltip)
                .style("left", (d3.event.pageX + 15) + "px")
                .style("top", (d3.event.pageY - 10) + "px");
            })
            .on("mouseout", mouseout);

    // Creating Squares on line
    segment.selectAll("dot")
            .data(function (d) { return d.PathwaysData; } )
            .enter()
            .append("rect")
            .filter(function(d) { return d.value[1] == "ARTIFICIAL"})
            .attr("class",function(d) { return d.value[1] })
            .attr("x", function (d) { return (x(d.Pathway) + x.rangeBand()/2)-105; })
            .attr("y", function (d) { return y(d.value[0]) -5; })
            .attr("width", 10)
            .attr("height", 10)
            .style("stroke", "black")
            .style("fill", "black")
            .on("mouseover", mouseover)
            .on("mousemove", function (d) {
                divToolTip
                // .text(d.value[0].toFixed(2))
                .html(d.tooltip)
                .style("left", (d3.event.pageX + 15) + "px")
                .style("top", (d3.event.pageY - 10) + "px");
            })
            .on("mouseout", mouseout);

     // Adding Tooltip
    var divToolTip = d3.select("body").append("div")
                .attr("class", "tooltip")
                .style("opacity", 1e-6);

    function mouseover() {
        divToolTip.transition()
            .duration(500)
            .style("opacity", 1);
    }
    function mouseout() {
        divToolTip.transition()
            .duration(500)
            .style("opacity", 1e-6);
    }

  if(Keys.length < 20){
    function position(d,i) {
      var c = 1;   // number of columns
      var h = Keys.length;  // height of each entry
      var w = 90; // width of each entry (so you can position the next column)
      var x = width + i % c * w - 130;
      var y = i*13;
      return "translate(" + x + "," + y + ")";
    }
   }else if(Keys.length > 40){
     function position(d,i) {
       var c = 3;   // number of columns
       var h = Keys.length;  // height of each entry
       var w = 110; // width of each entry (so you can position the next column)
       var x = width + i % c * w - 150;
       var y = Math.floor(i / c)*13;
       return "translate(" + x + "," + y + ")";
     }
   }else{
      function position(d,i) {
        var c = 2;   // number of columns
        var h = Keys.length;  // height of each entry
        var w = 110; // width of each entry (so you can position the next column)
        var x = width + i % c * w - 150;
        var y = Math.floor(i / c)*13;
        return "translate(" + x + "," + y + ")";
      }
  }

// setting the legend positions depending on paths numbering
  var path_nr = xData.length;

  if(path_nr >= 3){
    var xSeed = 6;
  }else{
    var xSeed = -76;
  }

    svg.append("g")
       .attr("class", "ytitle")
       .attr("transform", position)
          .append("text")
            .attr("x", xSeed)
            .attr("y", -10)
            .text("Ligands tested for bias by")
            .attr("text-anchor", "left")
            .attr("font-weight", "bold")
            .style("font", "12px sans-serif")
            .style("alignment-baseline", "middle");

    svg.append("g")
       .attr("class", "ytitle")
       .attr("transform", position)
          .append("text")
            .attr("x", xSeed)
            .attr("y", 4)
            .text("decreasing " + ylabel)
            .attr("text-anchor", "left")
            .attr("font-weight", "bold")
            .style("font", "12px sans-serif")
            .style("alignment-baseline", "middle");

    svg.append("g")
       .attr("class", "ytitle")
       .attr("transform", position)
          .append("text")
            .attr("x", xSeed)
            .attr("y", 21)
            .text("Endogenous ligand: " + reference + " (0.0)")
            .attr("text-anchor", "left")
            .attr("font-weight", "bold")
            .style("font", "12px sans-serif")
            .style("alignment-baseline", "middle");


    svg.append("g")
       .attr("class", "ytitle")
       .attr("transform", position)
          .append("text")
            .attr("x", xSeed+12)
            .attr("y", 38)
            .text("Qualitative data point")
            .attr("text-anchor", "left")
            .style("font", "12px sans-serif")
            .style("alignment-baseline", "middle");

    svg.append("g")
       .attr("class", "ytitle")
       .attr("transform", position)
          .append("rect")
            .attr("x", xSeed)
            .attr("y", 33)
            .attr("width", 9)
            .attr("height", 9)
            .style("stroke", "black")
            .style("fill", "black");

var legend = svg.selectAll("mylabels")
      .data(Keys)
      .enter()
      .append("g")
      .attr("transform", position);

      legend.append("text")
            .attr("x", xSeed+12)
            .attr("y", 53)
            .style("fill", function(a){ return legendColor(a)})
            .text(function(d) {
                    if (d[1].length > 16) {
                        return d[1].substring(0,16)+"..."
                    }else {
                        return d[1]
                    }
            })
            .attr("id", function(d) { return d[0]})
            .style("font", "12px sans-serif")
            .attr("text-anchor", "left")
            .style("alignment-baseline", "middle")
            .on("click", function (d) {
                var tempId = d3.select(this).attr("id");
                d3.selectAll("g.segment circle")
                   .style("opacity", 0.2);
                d3.selectAll("g.segment rect")
                   .style("opacity", 0.2);
                d3.selectAll("g.segment path")
                   .style("opacity", 0.2);
                d3.selectAll("g.segment." + tempId + " path")
                  .style("opacity", 1);
                d3.selectAll("g.segment." + tempId + " circle")
                  .style("opacity", 1);
                d3.event.stopPropagation();
            });

legend.append("circle")
  .attr("cx",xSeed+4)
  .attr("cy",53)
  .attr("r", 5)
  .style("stroke", "black")
  .attr("fill",function(a) { return legendColor(a); })
  .on("click", function (d) {
      var tempId = d3.select(this).attr("id");
      d3.selectAll("g.segment circle")
         .style("opacity", 0.2);
      d3.selectAll("g.segment rect")
         .style("opacity", 0.2);
      d3.selectAll("g.segment path")
         .style("opacity", 0.2);
      d3.selectAll("g.segment." + tempId + " path")
        .style("opacity", 1);
      d3.selectAll("g.segment." + tempId + " circle")
        .style("opacity", 1);
      d3.event.stopPropagation();
  });
}

function DotScatter(data, BaseDiv, ID, colors, legendData, header, ylabel, qualitative){

  // check if there are qualitative values in the dataset
  // for specify the y labels later
  var hb = null;
  var fb = null;

  function step_rounder(value, step) {
      step || (step = 1.0);
      var inv = 1.0 / step;
      return Math.round(value * inv) / inv;
  }

  function highlight(object) {
    d3.selectAll("circle")
       .style("opacity", 0.2);
    d3.selectAll("path")
       .style("opacity", 0.2);
    d3.selectAll("rect")
       .style("opacity", 0.2);
    d3.selectAll("circle#" + object)
      .style("opacity", 1);
    d3.selectAll("path#" + object)
      .style("opacity", 1);
    d3.selectAll("rect#" + object)
      .style("opacity", 1);
    d3.selectAll(".Legend")
      .style("opacity", 1);
    d3.event.stopPropagation();
  }

  function ResetOpacity(){
    d3.selectAll("circle")
      .style("opacity", 1.0);
    d3.selectAll("rect")
      .style("opacity", 1.0);
    d3.selectAll("path")
      .style("opacity", 1.0);
  }

  function mouseover() {
      divToolTipTest.transition()
          .duration(500)
          .style("opacity", 1);
  }

  function mouseout() {
      divToolTipTest.transition()
          .duration(500)
          .style("opacity", 0);
  }
  //change this because no longer relevant
  if(legendData.length < 24){
    var margin = {top: 20, right: 200, bottom: 100, left: 60};
    width = 1000 - margin.left - margin.right;
    xSeed = margin.right;
    var jitterWidth = 30;
    function position(d,i) {
      var c = 1;   // number of columns
      var h = colors.length;  // height of each entry
      var w = 90; // width of each entry (so you can position the next column)
      var x = width + i % c * w - 130;
      var y = i*13;
      return "translate(" + x + "," + y + ")";
    }
  }else if(legendData.length > 40){
    var margin = {top: 20, right: 365, bottom: 100, left: 60};
    width = 1000 - margin.left - margin.right;
    xSeed = 260;
    var jitterWidth = 5;
     function position(d,i) {
       var c = 4;   // number of columns
       var h = colors.length;  // height of each entry
       var w = 80; // width of each entry (so you can position the next column)
       var x = width + i % c * w - 150;
       var y = Math.floor(i / c)*13;
       return "translate(" + x + "," + y + ")";
     }
   }else{
     var margin = {top: 20, right: 200, bottom: 100, left: 60};
     width = 1000 - margin.left - margin.right;
     xSeed = margin.right;
     var jitterWidth = 15;
      function position(d,i) {
        var c = 2;   // number of columns
        var h = colors.length;  // height of each entry
        var w = 90; // width of each entry (so you can position the next column)
        var x = width + i % c * w - 150;
        var y = Math.floor(i / c)*13;
        return "translate(" + x + "," + y + ")";
      }
  }

    height = 500 - margin.top - margin.bottom;

    var parentDiv = document.getElementById(BaseDiv)
    var title = document.createElement("div");
        title.setAttribute("id", "title_"+ID);
        parentDiv.appendChild(title);
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
        text_png.setAttribute("href", "javascript:saveSvgAsPng(document.getElementById('"+ID+"').getElementsByTagName('svg')[0], '"+ header + "_" + ID +".png');");
        li_png.appendChild(text_png);
        optionList.appendChild(li_png);

    var li_jpg = document.createElement("li");
    var text_jpg = document.createElement("a");
        text_jpg.textContent = "JPG";
        text_jpg.setAttribute("href", "javascript:saveSvgAsJpg(document.getElementById('"+ID+"').getElementsByTagName('svg')[0], '"+ header + "_" + ID +".jpg');");
        li_jpg.appendChild(text_jpg);
        optionList.appendChild(li_jpg);

    var li_tiff = document.createElement("li");
    var text_tiff = document.createElement("a");
        text_tiff.textContent = "TIFF";
        text_tiff.setAttribute("href", "javascript:saveSvgAsTiff(document.getElementById('"+ID+"').getElementsByTagName('svg')[0], '"+ header + "_" + ID +".tiff');");
        li_tiff.appendChild(text_tiff);
        optionList.appendChild(li_tiff);

    var li_svg = document.createElement("li");
    var text_svg = document.createElement("a");
        text_svg.textContent = "SVG";
        text_svg.setAttribute("href", "javascript:saveSvg(document.getElementById('"+ID+"'), '"+ header + "_" + ID +".svg');");
        li_svg.appendChild(text_svg);
        optionList.appendChild(li_svg);

        downloadDiv.appendChild(optionList);

    document.getElementById(ID).onclick = ResetOpacity;

    d3.select("#title_"+ID)
          .append("foreignObject")
            .attr("width", width)
            .attr("height", height)
            .attr("class", "title")
          .append("xhtml:body")
            .style("font", "15px 'Arial'")
            .style("padding-bottom", "3px")
            .style("padding-top", "15px")
            .html(header);


  var x = d3.scale.ordinal()
    .domain(data.map(function(d) {
      return d[0];
    }))
    .rangePoints([0, width], 0.5);


    var y = d3.scale.linear()
      .domain([
        d3.min(data, function(d) {
        return d[1];
      }),
        d3.max(data, function(d) {
        return d[1];
      })])
      .range([height, 0]);

  var chart = d3.select("#" + ID)
    .append("svg:svg")
    .attr("width", width + margin.right + margin.left)
    .attr("height", height + margin.top + margin.bottom)
    .attr("class", "chart")

  var main = chart.append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
    .attr("width", width)
    .attr("height", height)
    .attr("class", "main");

  // draw the x axis
  var xAxis = d3.svg.axis()
    .scale(x)
    .orient("bottom");

  main.append("g")
    .attr("transform", "translate(0," + height + ")")
    .attr("class", "main axis date")
    .style("font", "10px sans-serif")
    .call(xAxis)
    .selectAll("text")
      .attr("y", 0)
      .attr("x", 9)
      .attr("dy", ".35em")
      .attr("transform", "rotate(45)")
      .style("text-anchor", "start");

  // draw the y axis
  var yAxis = d3.svg.axis()
    .scale(y)
    .orient("left");

  var step = yAxis.scale().ticks()[1];

  for(i = 0; i < data.length; i++) {
    if(data[i][4] == 'High Bias'){
      tmp = data[i][1];
      hb = step_rounder(tmp, step);
    } else if(data[i][4] == 'Full Bias'){
      tmp = data[i][1];
      fb = step_rounder(tmp, step);
    }
  }

  yAxis.tickFormat(function(d) {
          if(d === fb){
            return "Full Bias";
          } else if(d === hb){
            return "High Bias";
          } else if(d === hb+step){
            return "";
          } else {
            return d;
          }
        });

  main.append("g")
      .attr("class", "y axis")
      .style("font", "10px sans-serif")
      .call(yAxis)
      .append("text")
          .attr("transform", "rotate(-90)")
          .attr("y", -40)
          .attr('x', -150)
          .attr("dy", ".71em")
          .style("font", "12px sans-serif")
          .style("text-anchor", "end")
          .text(ylabel);

  var divToolTipTest = d3.select("body")
              .append("div")
              .attr("class", "tooltip")
              .style("opacity", 0);

  var g = main.append("svg:g");

  // Creating Full Bias Squares
  g.selectAll("scatter-dots")
          .data(data)
          .enter()
          // .append("rect")
          .append("path")
          .filter(function(d) { return d[4] == 'Full Bias'})
          .attr("class","FullBias")
          .attr("d", d3.svg.symbol()
            .size( function(d) { return 6*6 })
            .type( function(d) { return "diamond" }))
            .style("fill", function(d) {return d[2];})
            .style("opacity", 1.0)
            // .style("stroke", "black")
            .attr("transform", function (d, i) {
                var xTrans = x(d[0]) - jitterWidth/2 + Math.random()*jitterWidth;
                return "translate(" + xTrans + ",0)"
            })
          .attr("id", function(d) {return "LC" + d[3].replace(/\[|\]|\(|\)|\s|\,|\'/g,"");})
          // .attr("x", function(d) {return x(d[0]) - jitterWidth/2 + Math.random()*jitterWidth ;})
          // .attr("y", function(d) {return y(d[1]);})
          // .attr("width", 7)
          // .attr("height", 7)
          // .style("stroke", "black")
          // // .style("fill", "red")
          // .style("fill", function(d) {return d[2];})
          // .style("opacity", 1.0)
          .on("mouseover", mouseover)
          .on("mousemove", function (d) {
              divToolTipTest
                .style("left", (d3.event.pageX) + "px")
                .style("top", (d3.event.pageY) + "px")
                .html("<b>Compound Name:</b> " + d[3]
                    + "<br\><b>Compared Pathway:</b> " + d[5]
                    + "<br\><b>Plotted Value:</b> Full Bias");
          })
          .on("mouseout", mouseout)
          .on("click", function (d) {
              var tempId = d3.select(this).attr("id");
              highlight(tempId);
          });

    // Creating Full Bias Squares
    g.selectAll("scatter-dots")
            .data(data)
            .enter()
            .append("rect")
            .filter(function(d) { return d[4] == "High Bias"})
            .attr("class", "HighBias")
            .attr("id", function(d) {return "LC" + d[3].replace(/\[|\]|\(|\)|\s|\,/g,"");})
            .attr("x", function(d) {return x(d[0]) - jitterWidth/2 + Math.random()*jitterWidth ;})
            .attr("y", function(d) {return y(d[1]);})
            .attr("width", 7)
            .attr("height", 7)
            // .style("stroke", "black")
            // .style("fill", "orange")
            .style("fill", function(d) {return d[2];})
            .style("opacity", 1.0)
            .on("mouseover", mouseover)
            .on("mousemove", function (d) {
                divToolTipTest
                  .style("left", (d3.event.pageX) + "px")
                  .style("top", (d3.event.pageY) + "px")
                  .html("<b>Compound Name:</b> " + d[3]
                      + "<br\><b>Compared Pathway:</b> " + d[5]
                      + "<br\><b>Plotted Value:</b> High Bias");
            })
            .on("mouseout", mouseout)
            .on("click", function (d) {
                var tempId = d3.select(this).attr("id");
                highlight(tempId);
            });

  g.selectAll("scatter-dots")
    .data(data)
    .enter().append("svg:circle")
    .filter(function(d) { return d[4] == 'REAL'})
    .attr("cx", function(d) {return x(d[0]) - jitterWidth/2 + Math.random()*jitterWidth ;})
    .attr("cy", function(d) {return y(d[1]);})
    .attr("r", 4)
    .attr("id", function(d) {return "LC" + d[3].replace(/\[|\]|\(|\)|\s|\,/g,"");})
    .style("fill", function(d) {return d[2];})
    .style("opacity", 1.0)
    .on("mouseover", mouseover)
    .on("mousemove", function (d) {
        divToolTipTest
          .style("left", (d3.event.pageX) + "px")
          .style("top", (d3.event.pageY) + "px")
          .html("<b>Compound Name:</b> " + d[3]
              + "<br\><b>Compared Pathway:</b> " + d[5]
              + "<br\><b>Plotted Value:</b> " + d[1]);
    })
    .on("mouseout", mouseout)
    .on("click", function (d) {
        var tempId = d3.select(this).attr("id");
        highlight(tempId);
    });
    // .on("click", function (d) {
    //     console.log('clicked!');
    //     divToolTipTest
    //     .text(d[3])
    //     .style("left", (d3.event.pageX) + "px")
    //     .style("top", (d3.event.pageY) + "px");
    // });

    chart.append('g')
       .attr('class', 'ytitle')
       .attr("transform", position)
          .append("text")
            .attr("x", xSeed - 5)
            .attr("y", margin.top)
            .text("Top 20 Ligands") //ΔΔLog(Emax/EC50)
            .attr("text-anchor", "left")
            .attr("font-weight", "bold")
            .style("font", "10px sans-serif")
            .style("alignment-baseline", "middle");

    chart.append('g')
       .attr('class', 'ytitle')
       .attr("transform", position)
          .append("text")
            .attr("x", xSeed - 5)
            .attr("y", margin.top + 10)
            .text("sorted by decreasing value")
            .attr("text-anchor", "left")
            .attr("font-weight", "bold")
            .style("font", "10px sans-serif")
            .style("alignment-baseline", "middle");

// set the gradient for the quantitative circle
// in the legend fixed
  var gradient = chart.append("svg:defs")
      .append("svg:linearGradient")
      .attr("id", "gradient")
      .attr("x1", "0%")
      .attr("y1", "0%")
      .attr("x2", "100%")
      .attr("y2", "100%")
      .attr("spreadMethod", "pad");

  gradient.append("svg:stop")
     .attr("offset", "0%")
     .attr("stop-color", "#66ffff")
     .attr("stop-opacity", 1);

  gradient.append("svg:stop")
      .attr("offset", "100%")
      .attr("stop-color", "#ff00ff")
      .attr("stop-opacity", 1);

  var legend = chart.selectAll("mylabels")
        .data(legendData)
        .enter()
        .append("g")
        .attr("transform", position);

  if(qualitative === true){

    chart.append('g')
       .attr("transform", position)
          // .append("rect")
          .append("path")
            .attr("d", d3.svg.symbol()
              .size( function(d) { return 6*6 })
              .type( function(d) { return "diamond" }))
              .style("stroke", "black")
              .style("fill", "#FAFAFA") //'url(#gradient)'
              .attr("class", "Legend")
              .attr("transform", "translate(" + (xSeed + 10) +"," + (margin.top + 28) +")")
              // .style("stroke", "black")
              // .attr("x", xSeed + 6)
              // .attr("y", margin.top + 25)
            // .attr("x", xSeed + 6)
            // .attr("y", margin.top + 25)
            // .attr('width', 8)
            // .attr('height', 8)
            // .style("stroke", "black")
            // .style("fill", "red")
            // .attr("class", "FullBias")
            // .style("stroke", "black")
            .on("click", function (d) {
                var tempId = d3.select(this).attr("id");
                highlight(tempId);
            });

    chart.append('g')
       .attr('class', 'ytitle')
       .attr("transform", position)
          .append("text")
            .attr("x", xSeed + 20)
            .attr("y", margin.top + 30)
            .text("Full Bias qualitative point")
            .attr("text-anchor", "left")
            .attr("font-weight", "bold")
            .style("font", "10px sans-serif")
            .style("alignment-baseline", "middle")
            .on("click", function (d) {
                var tempId = d3.select(this).attr("id");
                highlight(tempId);
            });

    chart.append('g')
       .attr("transform", position)
          .append("rect")
            .attr("x", xSeed + 6)
            .attr("y", margin.top + 38)
            .attr('width', 8)
            .attr('height', 8)
            .style("fill", "#FAFAFA") //'url(#gradient)'
            .attr("class", "Legend")
            .style("stroke", "black")
            .on("click", function (d) {
                var tempId = d3.select(this).attr("id");
                highlight(tempId);
            });

    chart.append('g')
       .attr('class', 'ytitle')
       .attr("transform", position)
          .append("text")
            .attr("x", xSeed + 20)
            .attr("y", margin.top + 43)
            .text("High Bias qualitative point")
            .attr("text-anchor", "left")
            .attr("font-weight", "bold")
            .style("font", "10px sans-serif")
            .style("alignment-baseline", "middle")
            .on("click", function (d) {
                var tempId = d3.select(this).attr("id");
                highlight(tempId);
            });

    chart.append('g')
       .attr("transform", position)
          .append("circle")
            .attr("cx", xSeed + 10)
            .attr("cy", margin.top + 57)
            .attr("r", 4)
            .style('stroke', 'black')
            .attr("class", "Legend")
            .attr("fill", "#FAFAFA"); //'url(#gradient)'

    chart.append('g')
       .attr('class', 'ytitle')
       .attr("transform", position)
          .append("text")
            .attr("x", xSeed + 20)
            .attr("y", margin.top + 57)
            .attr("class", "Legend")
            .text("Quantitative point")
            .attr("text-anchor", "left")
            .attr("font-weight", "bold")
            .style("font", "10px sans-serif")
            .style("alignment-baseline", "middle");

    legend.append("text")
      .attr("x", xSeed + 20)
      .attr("y", margin.top + 75)
      .style("fill", function(d){ return colors[d]})
      .text(function(d) {
              if (d.length > 25) {
                  return d.substring(0,25)+'...'
              }else {
                  return d
              }
      })
      .attr("text-anchor", "left")
      .style("font", "10px sans-serif")
      .style("alignment-baseline", "middle")
      .on("click", function (d) {
          d3.selectAll("circle")
             .style("opacity", 0.2);
          d3.selectAll("rect")
             .style("opacity", 0.2);
          d3.selectAll("path")
            .style("opacity", 0.2);
          d3.selectAll("circle#LC" + d.replace(/\[|\]|\(|\)|\s|\,|\'/g,""))
            .style("opacity", 1);
          d3.selectAll("rect#LC" + d.replace(/\[|\]|\(|\)|\s|\,|\'/g,""))
            .style("opacity", 1);
          d3.selectAll("path#LC" + d.replace(/\[|\]|\(|\)|\s|\,|\'/g,""))
            .style("opacity", 1);
          d3.selectAll(".Legend")
            .style("opacity", 1);
          d3.event.stopPropagation();
      });

  } else {

    chart.append('g')
       .attr("transform", position)
          .append("circle")
            .attr("class", "Legend")
            .attr("cx", xSeed + 10)
            .attr("cy", margin.top + 23)
            .attr("r", 4)
            .style("stroke", "black")
            .attr("fill","#FAFAFA"); //'url(#gradient)'

    chart.append('g')
       .attr('class', 'ytitle')
       .attr("transform", position)
          .append("text")
            .attr("x", xSeed + 20)
            .attr("y", margin.top + 24)
            .text("Quantitative point")
            .attr("text-anchor", "left")
            .attr("font-weight", "bold")
            .style("font", "10px sans-serif")
            .style("alignment-baseline", "middle");

    legend.append("text")
      .attr("x", xSeed + 20)
      .attr("y", margin.top + 40)
      .style("fill", function(d){ return colors[d]})
      .text(function(d) {
              if (d.length > 25) {
                  return d.substring(0,25)+'...'
              }else {
                  return d
              }
      })
      .attr("text-anchor", "left")
      .style("font", "10px sans-serif")
      .style("alignment-baseline", "middle")
      .on("click", function (d) {
          d3.selectAll("circle")
             .style("opacity", 0.2);
          d3.selectAll("rect")
             .style("opacity", 0.2);
          d3.selectAll("path")
            .style("opacity", 0.2);
          d3.selectAll("circle#LC" + d.replace(/\[|\]|\(|\)|\s|\,|\'/g,""))
            .style("opacity", 1);
          d3.selectAll("rect#LC" + d.replace(/\[|\]|\(|\)|\s|\,|\'/g,""))
            .style("opacity", 1);
          d3.selectAll("path#LC" + d.replace(/\[|\]|\(|\)|\s|\,|\'/g,""))
            .style("opacity", 1);
          d3.selectAll(".Legend")
            .style("opacity", 1);
          d3.event.stopPropagation();
      });
  }


// This is the dot for the ligands which is commented right now
  // legend.append("circle")
  //   .attr("cx", xSeed + 10)
  //   .attr("cy", margin.top + 55)
  //   .attr("r", 5)
  //   .attr("id", function(d) {return "LC" + d.replace(/\[|\]|\(|\)|\s|\,/g,"");})
  //   .style('stroke', 'black')
  //   .attr("fill",function(d) { return colors[d] })
  //   .on("click", function (d) {
  //       d3.selectAll("circle")
  //          .style("opacity", 0.2);
  //       d3.selectAll("rect")
  //          .style("opacity", 0.2);
  //       d3.selectAll("circle#LC" + d.replace(/\[|\]|\(|\)|\s|\,/g,""))
  //         .style("opacity", 1);
  //       d3.event.stopPropagation();
  //   });

}

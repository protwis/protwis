/*global d3*/
/*eslint complexity: ["error", 20]*/

function DotScatter(data, BaseDiv, ID, colors, legendData, header, ylabel, qualitative){

  // check if there are qualitative values in the dataset
  // for specify the y labels later
  var additional_gap;
  var y;
  var hb;
  var fb;
  var spacers;
  var idRemove;
  var pubs = new Array();
  var first = ylabel.replace(/\Δ|\(|\)|[Log]/g,"").split("/")[0];
  var second = ylabel.replace(/\Δ|\(|\)|[Log]/g,"").split("/")[1];
  for (var i = 0; i < data.length; i++){
    pubs.indexOf(data[i][0]) === -1 ? pubs.push(data[i][0]) : console.log();
  }

  function step_rounder(value, step) {
      step || (step = 1.0);
      var inv = 1.0 / step;
      return Math.round(value * inv) / inv;
  }

  function highlight(object) {
    d3.selectAll("circle")
       .style("opacity", 0.2);
    d3.selectAll("rect")
       .style("opacity", 0.2);
    d3.selectAll("circle#" + object)
      .style("opacity", 1);
    d3.selectAll("rect#" + object)
      .style("opacity", 1);
    d3.selectAll(".Legend")
      .style("opacity", 1);
    d3.selectAll(".divider")
      .style("opacity", 1);
    d3.selectAll(".domain")
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

  var margin = {top: 20, right: 200, bottom: 100, left: 150};
  // width = 1000 - margin.left - margin.right;
  var width = pubs.length * 62;
  var xSeed = 0;
  var jitterWidth = 30;
  function position(d,i) {
    var c = 1;   // number of columns
    var h = colors.length;  // height of each entry
    var w = 90; // width of each entry (so you can position the next column)
    var x = 56; // adjusted for being on the left
    var y = i*13;
    return "translate(" + x + "," + y + ")";
  }

    var height = 600 - margin.top - margin.bottom;

    var parentDiv = document.getElementById(BaseDiv);

    var nestedDiv = document.createElement("div");
      nestedDiv.setAttribute("id", ID);
      nestedDiv.setAttribute("style", "margin-top: 50px;");
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

  var x = d3.scale.ordinal()
    .domain(data.map(function(d) {
      return d[0];
    }))
    .rangePoints([0, width], 0.5);

  if(qualitative === true){
    y = d3.scale.linear()
      .domain([0,
        d3.max(data, function(d) {
        return d[1];
      })])
      .range([height, 0]);

  }else{
    y = d3.scale.linear()
      .domain([
        d3.min(data, function(d) {
        return d[1];
      }),
        d3.max(data, function(d) {
        return d[1];
      })])
      .range([height, 0]);

  }
  var chart = d3.select("#" + ID)
    .append("svg:svg")
    .attr("width", width + margin.right + margin.left)
    .attr("height", height + margin.top + margin.bottom)
    .attr("class", "chart")

  var main = chart.append("g")
    .attr("transform", "translate(" + (margin.left + 120) + "," + (margin.top/2) + ")")
    .attr("width", width)
    .attr("height", height  + margin.top + margin.bottom)
    .attr("class", "main");

  // draw the x axis
  var xAxis = d3.svg.axis()
    .scale(x)
    .orient("bottom");

  main.append("g")
    .attr("transform", "translate(0," + height + ")")
    .attr("class", "main axis date")
    .style("font", "8px sans-serif")
    .call(xAxis)
    .selectAll("text")
      .attr("y", 5)
      .attr("x", 9)
      .attr("class", "column_label")
      .attr("dy", ".35em")
      .attr("transform", "rotate(45)")
      .style("text-anchor", "start");

  // draw the y axis
  var yAxis = d3.svg.axis()
    .scale(y)
    .orient("left");

  var tickArray = yAxis.scale().ticks();
  // console.log(tickArray);
  // console.log(tickArray.slice(-1)[0]);
  var step = yAxis.scale().ticks()[1];
  var count_ticks = yAxis.scale().ticks().length;

  for(let i = 0; i < data.length; i++) {
    if(data[i][4] == "High Bias"){
      let tmp = data[i][1];
      if(tickArray.includes(tmp)){
        hb = step_rounder(tmp, step);
      }else{
        hb = tickArray.slice(-1)[0];
      }
    } else if(data[i][4] == "Full Bias"){
      let tmp = data[i][1];
      if(tickArray.includes(tmp)){
        fb = step_rounder(tmp, step);
      }else{
        fb = tickArray.slice(-1)[0];
      }
      // hb = step_rounder((tmp - step), step);
    }
  }

  yAxis.tickFormat(function(d) {
          if(d === fb){
            return "Full Bias";
          } else if(d === hb){
            return "High Bias";
          } else {
            return d;
          }
        });

  var gap = Math.round(380/count_ticks);
  var dividers = Math.round(380 - gap);

  main.append("g")
      .attr("class", "y axis")
      .style("font", "10px sans-serif")
      .call(yAxis)
      .append("text")
          .attr("transform", "rotate(-90)")
          .attr("y", -40)
          .attr("x", -150)
          .attr("dy", ".71em")
          .style("font", "12px sans-serif")
          .style("text-anchor", "end")
          .text(ylabel);

  // additional gap calculated for ticks not on the top of the axis
  var node_checks = d3.selectAll(".y.axis").selectAll(".tick").each(function(d){
      additional_gap = d3.select(this).attr("transform");
      additional_gap = Math.round(additional_gap.replace(/\(|\)/g,"").split(",")[1]);
  });

  if(hb === null && fb !== null){
    spacers = gap * (tickArray.indexOf(fb) - 1);
    main.append("rect")
        .attr("class", "divider")
        .attr("x", -5)
        .attr("y", dividers - spacers + additional_gap)
        .attr("height", 6)
        .attr("width", 10);

    main.append("rect")
        .attr("class", "divider")
        .attr("x", -5)
        .attr("y", (dividers - spacers + additional_gap) + 1.5)
        .attr("height", 3)
        .attr("width", 10)
        .style("fill", "white");
  }else if(hb !== null && fb === null){
    spacers = gap * (tickArray.indexOf(hb) - 1);
    main.append("rect")
        .attr("class", "divider")
        .attr("x", -5)
        .attr("y", dividers - spacers + additional_gap)
        .attr("height", 6)
        .attr("width", 10);

    main.append("rect")
        .attr("class", "divider")
        .attr("x", -5)
        .attr("y", (dividers - spacers + additional_gap) + 1.5)
        .attr("height", 3)
        .attr("width", 10)
        .style("fill", "white");
  }else if(hb !== null && fb !== null){
    spacers = gap * (tickArray.indexOf(hb) - 1);
    main.append("rect")
        .attr("class", "divider")
        .attr("x", -5)
        .attr("y", dividers - spacers + additional_gap)
        .attr("height", 6)
        .attr("width", 10);

    main.append("rect")
        .attr("class", "divider")
        .attr("x", -5)
        .attr("y", (dividers - spacers + additional_gap) + 1.5)
        .attr("height", 3)
        .attr("width", 10)
        .style("fill", "white");
  }

  // remove tick in between full bias and high bias IF it exists
  if((tickArray.indexOf(fb) - tickArray.indexOf(hb)) === 2){
    idRemove = tickArray[tickArray.indexOf(fb) - 1];
    nodes = d3.selectAll(".y.axis").selectAll(".tick").each(function(d){
      if(d3.select(this).text() == idRemove){
        d3.select(this)[0][0].innerHTML = "";
      };
    });
  };

  var divToolTipTest = d3.select("body")
              .append("div")
              .attr("class", "tooltip")
              .style("opacity", 0);

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

  var g = main.append("svg:g");

    // Creating Full Bias Squares
    g.selectAll("scatter-dots")
            .data(data)
            .enter()
            .append("rect")
            .filter(function(d) { return d[4] === "Full Bias";})
            .attr("class", "FullBias")
            .attr("id", function(d) {return "LC" + d[3].replace(/\[|\]|\(|\)|\s|\,/g,"");})
            .attr("x", function(d) {return x(d[0]) - jitterWidth/2 + Math.random()*jitterWidth ;})
            .attr("y", function(d) {return y(d[1]);})
            .attr("width", 7)
            .attr("height", 7)
            .style("fill", function(d) {return d[2];})
            .style("opacity", 1.0)
            .on("mouseover", mouseover)
            .on("mousemove", function (d) {
                divToolTipTest
                  .style("left", (d3.event.pageX) + "px")
                  .style("top", (d3.event.pageY) + "px")
                  .html(d[5]);
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
            .filter(function(d) { return d[4] === "High Bias";})
            .attr("class", "HighBias")
            .attr("id", function(d) {return "LC" + d[3].replace(/\[|\]|\(|\)|\s|\,/g,"");})
            .attr("x", function(d) {return x(d[0]) - jitterWidth/2 + Math.random()*jitterWidth ;})
            .attr("y", function(d) {return y(d[1]);})
            .attr("width", 7)
            .attr("height", 7)
            .style("fill", function(d) {return d[2];})
            .style("opacity", 1.0)
            .on("mouseover", mouseover)
            .on("mousemove", function (d) {
                divToolTipTest
                  .style("left", (d3.event.pageX) + "px")
                  .style("top", (d3.event.pageY) + "px")
                  .html(d[5]);
            })
            .on("mouseout", mouseout)
            .on("click", function (d) {
                var tempId = d3.select(this).attr("id");
                highlight(tempId);
            });

  g.selectAll("scatter-dots")
    .data(data)
    .enter().append("svg:circle")
    .filter(function(d) { return (d[4] !== "Full Bias" && d[4] !== "High Bias");})
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
          .html(d[5]);
    })
    .on("mouseout", mouseout)
    .on("click", function (d) {
        var tempId = d3.select(this).attr("id");
        highlight(tempId);
    });

  chart.append("g")
     .attr("class", "header")
     .attr("transform", position)
        .append("text")
          .attr("x", xSeed)
          .attr("y", margin.top)
          .text(header) //ΔΔLog(Emax/EC50)
          .attr("text-anchor", "left")
          .style("font", "14px sans-serif")
          .style("alignment-baseline", "middle");

    chart.append("g")
       .attr("class", "ytitle")
       .attr("transform", position)
          .append("text")
            .attr("x", xSeed)
            .attr("y", margin.top + 25)
            .text("Top 20 Ligands") //ΔΔLog(Emax/EC50)
            .attr("text-anchor", "left")
            .attr("font-weight", "bold")
            .style("font", "10px sans-serif")
            .style("alignment-baseline", "middle");

    chart.append("g")
       .attr("class", "ytitle")
       .attr("transform", position)
          .append("text")
            .attr("x", xSeed)
            .attr("y", margin.top + 35)
            .text("sorted by decreasing value")
            .attr("text-anchor", "left")
            .attr("font-weight", "bold")
            .style("font", "10px sans-serif")
            .style("alignment-baseline", "middle");

  var legend = chart.selectAll("mylabels")
        .data(legendData)
        .enter()
        .append("g")
        .attr("transform", position);

  if(qualitative === true){

    chart.append("g")
       .attr("transform", position)
          .append("rect")
            .attr("x", xSeed + 6)
            .attr("y", margin.top + 48)
            .attr("width", 8)
            .attr("height", 8)
            .style("fill", "#FAFAFA") //'url(#gradient)'
            .attr("class", "Legend")
            .style("stroke", "black")
            .on("click", function (d) {
                var tempId = d3.select(this).attr("id");
                highlight(tempId);
            });

    chart.append("g")
       .attr("class", "ytitle")
       .attr("transform", position)
          .append("text")
            .attr("x", xSeed + 20)
            .attr("y", margin.top + 53)
            .text("Qualitative data point")
            .attr("text-anchor", "left")
            .attr("font-weight", "bold")
            .style("font", "10px sans-serif")
            .style("alignment-baseline", "middle")
            .on("click", function (d) {
                var tempId = d3.select(this).attr("id");
                highlight(tempId);
            });

    chart.append("g")
       .attr("transform", position)
          .append("circle")
            .attr("cx", xSeed + 10)
            .attr("cy", margin.top + 67)
            .attr("r", 4)
            .style("stroke", "black")
            .attr("class", "Legend")
            .attr("fill", "#FAFAFA"); //'url(#gradient)'

    chart.append("g")
       .attr("class", "ytitle")
       .attr("transform", position)
          .append("text")
            .attr("x", xSeed + 20)
            .attr("y", margin.top + 67)
            .attr("class", "Legend")
            .text("Quantitative data point")
            .attr("text-anchor", "left")
            .attr("font-weight", "bold")
            .style("font", "10px sans-serif")
            .style("alignment-baseline", "middle");

    legend.append("rect")
      .attr("x", xSeed + 6)
      .attr("y", margin.top + 83.5)
      .attr("width", 8)
      .attr("height", 2)
      .attr("id", function(d) {return "LC" + d.replace(/\[|\]|\(|\)|\s|\,/g,"");})
      .style("stroke", function(d) { return colors[d];})
      .attr("fill", function(d) { return colors[d];})
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
          d3.selectAll(".domain")
            .style("opacity", 1);
          d3.selectAll(".divider")
            .style("opacity", 1);
          d3.event.stopPropagation();
      });

    legend.append("text")
      .attr("x", xSeed + 20)
      .attr("y", margin.top + 85)
      // .style("fill", function(d){ return colors[d]})
      .text(function(d) {
              if (d.length > 25) {
                  return d.substring(0,25)+"...";
              }else {
                  return d;
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
          d3.selectAll(".domain")
            .style("opacity", 1);
          d3.selectAll(".divider")
            .style("opacity", 1);
          d3.event.stopPropagation();
      });

  } else {

    chart.append("g")
       .attr("transform", position)
          .append("circle")
            .attr("class", "Legend")
            .attr("cx", xSeed + 10)
            .attr("cy", margin.top + 48)
            .attr("r", 4)
            .style("stroke", "black")
            .attr("fill","#FAFAFA"); //'url(#gradient)'

    chart.append("g")
       .attr("class", "ytitle")
       .attr("transform", position)
          .append("text")
            .attr("x", xSeed + 20)
            .attr("y", margin.top + 49)
            .text("Quantitative data point")
            .attr("text-anchor", "left")
            .attr("font-weight", "bold")
            .style("font", "10px sans-serif")
            .style("alignment-baseline", "middle");

    legend.append("rect")
      .attr("x", xSeed + 6)
      .attr("y", margin.top + 63)
      .attr("width", 8)
      .attr("height", 2)
      .attr("id", function(d) {return "LC" + d.replace(/\[|\]|\(|\)|\s|\,/g,"");})
      .style("stroke", function(d) { return colors[d] })
      .attr("fill", function(d) { return colors[d] })
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
          d3.selectAll(".domain")
            .style("opacity", 1);
          d3.selectAll(".divider")
            .style("opacity", 1);
          d3.event.stopPropagation();
      });

    legend.append("text")
      .attr("x", xSeed + 20)
      .attr("y", margin.top + 65)
      // .style("fill", function(d){ return colors[d]})
      .text(function(d) {
              if (d.length > 25) {
                  return d.substring(0,25)+"..."
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
          d3.selectAll(".domain")
            .style("opacity", 1);
          d3.selectAll(".divider")
            .style("opacity", 1);
          d3.event.stopPropagation();
      });
  }

}

function DotScatter(data, BaseDiv, ID, colors, legend, header){

  var margin = {
      top: 20, right: 100, bottom: 100, left: 60
    },
    width = 1000 - margin.left - margin.right,
    height = 500 - margin.top - margin.bottom;

    var parentDiv = document.getElementById(BaseDiv)
    var title = document.createElement("div");
        title.setAttribute("id", "title_"+ID);
        parentDiv.appendChild(title);
    var nestedDiv = document.createElement("div");
      nestedDiv.setAttribute("id", ID);
      parentDiv.appendChild(nestedDiv);

    function ResetOpacity(){
      d3.selectAll("circle")
        .style("opacity", 1)
    };

    document.getElementById(ID).onclick = ResetOpacity;

    d3.select("#title_"+ID)
          .append("foreignObject")
            .attr("width", width)
            .attr("height", height)
            .attr("class", "title")
          .append("xhtml:body")
            .style("font", "15px "Arial"")
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

  var jitterWidth = 10;

  main.append("g")
      .attr("class", "y axis")
      .call(yAxis)
      .append("text")
          .attr("transform", "rotate(-90)")
          .attr("y", -40)
          .attr("dy", ".71em")
          .style("text-anchor", "end")
          .text("ΔΔLog(Emax/EC50)");

  main.append("g")
    .attr("transform", "translate(0,0)")
    .attr("class", "main axis date")
    .call(yAxis);

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
          .style("opacity", 1e-6);
  }

  var g = main.append("svg:g");

  g.selectAll("scatter-dots")
    .data(data)
    .enter().append("svg:circle")
    .attr("cx", function(d) {return x(d[0]) - jitterWidth/2 + Math.random()*jitterWidth ;})
    .attr("cy", function(d) {return y(d[1]);})
    .attr("r", 8)
    .attr("id", function(d) {return "LC" + d[3].replace(/\[|\]|\(|\)|\s/g,"");})
    .style("stroke", "black")
    .style("fill", function(d) {return d[2];})
    .on("mouseover", mouseover)
    .on("mousemove", function (d) {
        divToolTipTest
        .text(d[3])
        .style("left", (d3.event.pageX) + "px")
        .style("top", (d3.event.pageY) + "px");
    })
    .on("mouseout", mouseout)
    .on("click", function (d) {
        var tempId = d3.select(this).attr("id");
        d3.selectAll("circle")
           .style("opacity", 0.2);
        d3.selectAll("circle#" + tempId)
          .style("opacity", 1);
        d3.event.stopPropagation();
    });

    if(legend.length < 24){

      function position(d,i) {
        var c = 1;   // number of columns
        var h = colors.length;  // height of each entry
        var w = 90; // width of each entry (so you can position the next column)
        var x = width + i % c * w - 130;
        var y = i*13;
        return "translate(" + x + "," + y + ")";
      }

    }else if(legend.length > 40){

       function position(d,i) {
         var c = 4;   // number of columns
         var h = colors.length;  // height of each entry
         var w = 90; // width of each entry (so you can position the next column)
         var x = width + i % c * w - 150;
         var y = Math.floor(i / c)*13;
         return "translate(" + x + "," + y + ")";
       }
     }else{

        function position(d,i) {
          var c = 2;   // number of columns
          var h = colors.length;  // height of each entry
          var w = 90; // width of each entry (so you can position the next column)
          var x = width + i % c * w - 150;
          var y = Math.floor(i / c)*13;
          return "translate(" + x + "," + y + ")";
        }
    }

}

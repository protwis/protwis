
function DrawMultiLineChart(Data, BaseDiv, RevenueName, Keys, ID, header) {

    var parentDiv = document.getElementById(BaseDiv)
    var title = document.createElement("div");
        title.setAttribute("id", "title_"+ID);
        parentDiv.appendChild(title);
    var nestedDiv = document.createElement("div");
      nestedDiv.setAttribute("id", ID);
      parentDiv.appendChild(nestedDiv);

    var margin = { top: 20, right: 80, bottom: 30, left: 50 },
     width = 1000 - margin.left - margin.right,
     height = 350 - margin.top - margin.bottom;

    d3.select('#title_'+ID)
          .append("foreignObject")
            .attr("width", width)
            .attr("height", height)
            .attr('class', 'title')
          .append("xhtml:body")
            .style("font", "15px 'Arial'")
            .style("padding-bottom", "3px")
            .style("padding-top", "15px")
            .html(header);

    var x = d3.scale.ordinal()
            .rangeRoundBands([0, width]);

    var y = d3.scale.linear()
        .range([height, 0])

    var color = d3.scale.category10();

    var xAxis = d3.svg.axis()
        .scale(x)
        .orient('bottom');

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient('left')
        .ticks(10);

    // xData gives an array of distinct 'Pathways' for which trends chart is going to be made.
    var xData = Data[0].PathwaysData.map(function (d) { return d.Pathway; });
    //console.log(xData);

    var line = d3.svg.line()
        //.interpolate('basis')
        .x(function (d) { return x(d.Pathway) + x.rangeBand() / 2; })
        .y(function (d) { return y(d.value); });

    // document.BaseDiv.appendChild(div);

    var svg = d3.select('#' + ID)
        .append('svg')
        .attr('width', width + margin.left + margin.right)
        .attr('height', height + margin.top + margin.bottom)
        .append('g')
        .attr('transform', 'translate(' + margin.left + ',' + margin.top + ')');

    color.domain(Data.map(function (d) { return d.name; }));

    x.domain(xData);

    var valueMax = d3.max(Data, function (r) { return d3.max(r.PathwaysData, function (d) { return d.value; }) });
    var valueMin = d3.min(Data, function (r) { return d3.min(r.PathwaysData, function (d) { return d.value; }) });
    y.domain([valueMin, valueMax]);

    //Drawing X Axis
    svg.append('g')
            .attr('class', 'x axis')
            .attr('transform', 'translate(0,' + height + ')')
            .call(xAxis);

    // Drawing Horizontal grid lines.
    svg.append('g')
        .attr('class', 'GridX')
      .selectAll('line.grid').data(y.ticks()).enter()
        .append('line')
        .attr(
        {
            'class': 'grid',
            'x1': x(xData[0]),
            'x2': x(xData[xData.length - 1]) + x.rangeBand() / 2,
            'y1': function (d) { return y(d); },
            'y2': function (d) { return y(d); }
        });
    // Drawing Y Axis
    svg.append('g')
        .attr('class', 'y axis')
        .call(yAxis)
        .append('text')
            .attr('transform', 'rotate(-90)')
            .attr('y', 6)
            .attr('dy', '.71em')
            .style('text-anchor', 'end')
            .text(RevenueName);

    // Drawing Lines for each segments
    var segment = svg.selectAll('.segment')
                    .data(Data)
                    .enter().append('g')
                    .attr('class', function(d) { return "segment "+ d.name;});
                    // .attr('id', function (d) { return d.name; });

    segment.append('path')
            .attr('class', 'line')
            .attr('visible',1)
            .attr('d', function (d) { return line(d.PathwaysData); })
            .style('stroke', function (d) { return color(d.name); });

    // Creating Dots on line
    segment.selectAll('dot')
            .data(function (d) { return d.PathwaysData; })
            .enter().append('circle')
            .attr('r', 5)
            .attr('cx', function (d) { return x(d.Pathway) + x.rangeBand() / 2; })
            .attr('cy', function (d) { return y(d.value); })
            .style('stroke', 'black')
            .style('fill', function (d) { return color(this.parentNode.__data__.name); })
            .on('mouseover', mouseover)
            .on('mousemove', function (d) {
                divToolTip
                .text(this.parentNode.__data__.name +' : '+ d.value)
                .style('left', (d3.event.pageX + 15) + 'px')
                .style('top', (d3.event.pageY - 10) + 'px');
            })
            .on('mouseout', mouseout);

     // Adding Tooltip
    var divToolTip = d3.select('body').append('div')
                .attr('class', 'tooltip')
                .style('opacity', 1e-6);

    function mouseover() {
        divToolTip.transition()
            .duration(500)
            .style('opacity', 1);
    }
    function mouseout() {
        divToolTip.transition()
            .duration(500)
            .style('opacity', 1e-6);
    }

   // Add one dot in the legend for each name.
   svg.selectAll("mydots")
     .append('g')
     .data(Keys)
     .enter()
     .append("circle")
       .attr("cx", width-50)
       .attr("cy", function(d,i){ return i*13})
       .attr("r", 5)
       .style('stroke', 'black')
       .style("fill", function(d){ return color(d)});

   // Add one dot in the legend for each name.
   svg.selectAll("mylabels")
     .data(Keys)
     .enter()
     .append("text")
       .attr("x", width-30)
       .attr("y", function(d,i){ return i*13})
       .style("fill", function(a){ return color(a)})
       .text(function(d){ return d[1]})
       .attr("id", function(d) { return d[0]})
       .attr("text-anchor", "left")
       .style("alignment-baseline", "middle")
       .on('click', function (d) {
           var tempId = d3.select(this).attr('id');
           d3.selectAll('g.segment circle')
              .style('opacity', 0.2);
           d3.selectAll('g.segment path')
              .style('opacity', 0.2);
           // Hide or show the elements
           d3.selectAll('g.segment.' + tempId + ' path')
             .style('opacity', 1)
           d3.selectAll('g.segment.' + tempId + ' circle')
             .style('opacity', 1)
       });
}

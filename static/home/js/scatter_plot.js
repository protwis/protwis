function draw_scatter_plot(data, options) {
    // Draw a scatter plot of a difference between two protein profiles.
    // Radius of the circles corresponds to the conservation of the given feature.
    // The data contains a list of generic residue positions
    // and a list of pairs [index, color, value] where index is an index in
    // a list of generic residue positions, feature_color defines a color of the circle
    // and value is a conservation defference affecting the size of a circle.
    var w = 800,
        h = 300,
        pad = 20,
        left_pad = 100;

    var svg = d3.select('#'+options.anchor)
                .append("svg")
                .attr("width", w)
                .attr("height", h);

    var x = d3.scaleLinear().domain([0, options.xticks.length]).range([left_pad, w-pad]),
        y = d3.scaleLinear().domain([0, 2]).range([pad, h-pad*2]);
    console.info(options.xticks)
    var xAxis = d3.axisTop(x)
                .ticks(options.xticks.length)
                .tickValues(options.xticks);
        // Not really needing that atm
        yAxis = d3.axisLeft(y)
                //.ticks(1)
                .ticks(['Profile',]);

    svg.append("g")
        .attr("class", "axis")
        .call(xAxis);

    svg.append("g")
        .attr("class", "axis")
        .call(yAxis);

    svg.selectAll('circle')
        .data(data)
        .enter()
        .append('circle')
        .attr('class', 'circle')
        .attr('cx', function (d) { return x(d[0]); })
        .attr('cy', (h-pad*2)/2)
        .attr('r', function (d) {
            if (d[2] < 20) {
                return 1;
            }
            else if (d[2] < 50) {
                return 3;
            }
            else {
                return 5;
            }
        })
        .attr('fill', function (d) { return d[1]; })

}
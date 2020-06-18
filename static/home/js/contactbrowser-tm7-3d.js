function tm7_plot_3d(containerSelector, plot_data) {


    var set1_data = JSON.parse(JSON.stringify(plot_data["coordinates_set1"]));
    var set2_data = JSON.parse(JSON.stringify(plot_data["coordinates_set2"]));
    set1_data.forEach(function (entry) {
        entry.id = entry.label
        entry.color = "set1"
    });

    set2_data.forEach(function (entry) {
        entry.id = entry.label
        entry.color = "set2"
      });

    set_colors = [];
    set_colors['set1'] = "#99ccff";
    set_colors['set2'] = "#ff9999";
    // Set colors
    set_1_color = "#99ccff";
    set_2_color = "#ff9999";

    circle_r = 23;
    line_widths = 3;
    path_r = circle_r + 3 + line_widths;
    line_distance_from_center = path_r;
    values_font_size = 16;
    values_font_size_hiding = 6;
    tm_font_size = 8;

    minimum_angle_to_show = 10;
    minimum_distance_to_show = 1;

    padding = 53; //33
    min_y = Math.min.apply(Math, [...set1_data, ...set2_data].map(a => a.y)) - padding;
    max_y = Math.max.apply(Math, [...set1_data, ...set2_data].map(a => a.y)) + padding;

    min_x = Math.min.apply(Math, [...set1_data, ...set2_data].map(a => a.x)) - padding;
    max_x = Math.max.apply(Math, [...set1_data, ...set2_data].map(a => a.x)) + padding;

    $(containerSelector).html('')
    $(containerSelector).addClass("tm7Plot");
    $(containerSelector).css("position","relative");

    var svgContainer = d3v4.select(containerSelector).append("svg")
        // .attr("viewBox", min_x + " " + min_y + " " + (max_x - min_x) + " " + (max_y - min_y))
        .attr("viewBox", 0 + " " + 0 + " " + 600 + " " + 600)
        .attr("width", "100%")
        .attr("style", "height: 500px");


    var defs = svgContainer.append("defs");

    defs.append("marker")
        .attr("id", "arrowhead")
        .attr("refX", 2) /*must be smarter way to calculate shift*/
        .attr("refY", 2)
        .attr("markerWidth", 6)
        .attr("markerHeight", 4)
        .attr("viewBox", "0 0 10 10")
        .attr("orient", "auto")
        .append("path")
        .attr("d", "M 2,0 V 4 L6,2 Z")
        .attr("fill", "grey");

    defs.append("marker")
        .attr("id", "arrowhead-rev")
        .attr("refX", 5)
        .attr("refY", 5)
        .attr("markerWidth", 2)
        .attr("markerHeight", 2)
        .attr("viewBox", "0 0 10 10")
        .attr("orient", "auto-start-reverse")
        .append("path")
        .attr("d", "M 0 0 L 10 5 L 0 10 z")
        .attr("fill", "grey");


    //Filter for the outside glow
    var filter = defs.append("filter")
        .attr("id","glow");
    filter.append("feGaussianBlur")
        .attr("stdDeviation","1.5")
        .attr("result","coloredBlur");
    var feMerge = filter.append("feMerge");
    feMerge.append("feMergeNode")
        .attr("in","coloredBlur");
    feMerge.append("feMergeNode")
        .attr("in", "SourceGraphic");

    var filter2 = defs.append("filter")
        .attr("id", "shadow")
        .append("feDropShadow")
        .attr("dx", 0.2)
        .attr("dy", 1)
        .attr("stdDeviation",0.2);

    var filter3 = defs.append('filter').attr('id', 'lightMe');

    filter3.append('feDiffuseLighting')
        .attr('in', 'SourceGraphic')
        .attr('result', 'light')
        .attr('lighting-color', 'white')
        .append('fePointLight')
            .attr('x', 0.85 * 35)
            .attr('y', 0.85 * 35)
            .attr('z', 50);

    filter3.append('feComposite')
        .attr('in', 'SourceGraphic')
        .attr('in2', 'light')
        .attr('operator', 'arithmetic')
        .attr('k1', '1')
        .attr('k2', '0')
        .attr('k3', '0')
        .attr('k4', '0');

    var origin = [250, 250],
     j = 10,
     scale = 10,
     scatter = [],
     yLine = [],
     xGrid = [],
     beta = 0,
     alpha = 0,
     key = function(d){ return d.id; },
     startAngle = 0;
    var svgContainer    = svgContainer.call(d3v4.drag().on('drag', dragged).on('start', dragStart).on('end', dragEnd)).append('g');
    var color  = d3v4.scaleOrdinal(d3v4.schemeCategory10);
    var mx, my, mouseX, mouseY;

    var grid3d = d3v4._3d()
        .shape('GRID', 20)
        .origin(origin)
        .rotateY( startAngle)
        .rotateX(-startAngle)
        .scale(scale);

    var point3d = d3v4._3d()
        .x(function(d){ return d.x; })
        .y(function(d){ return d.y; })
        .z(function(d){ return d.z; })
        .origin(origin)
        .rotateY( startAngle)
        .rotateX(-startAngle)
        .scale(scale);

    var yScale3d = d3v4._3d()
        .shape('LINE_STRIP')
        .origin(origin)
        .rotateY( startAngle)
        .rotateX(-startAngle)
        .scale(scale);

    var lines;

    var latest_data;
    var sets;
    var g_sets, circles;
    function processData(data, tt, initial) {
        latest_data = data;
        for (var i of Object.keys(set1_data)) {
            ii = parseInt(i) + 7;
            distance = get_distance(data[1][i], data[1][ii]);
            direction = get_vector(data[1][i], data[1][ii]);
            cos_direction = direction[0] / direction[1];
            rad = Math.atan2(direction[1], direction[0]);
            degrees = -rad * (180 / Math.PI);
            degrees = degrees ? -degrees : 0;

            // scaled_rotation = Math.abs(set2_data[i].rotation) < 20 ? 20 : Math.abs(set2_data[i].rotation);
            // from = - scaled_rotation / 2;
            // to = scaled_rotation / 2;
            // cos1 = parseInt(Math.cos(toRadians(from)) * path_r);
            // sin1 = parseInt(Math.sin(toRadians(from)) * path_r);

            // cos2 = parseInt(Math.cos(toRadians(to)) * path_r);
            // sin2 = parseInt(Math.sin(toRadians(to)) * path_r);

            // arc_flag = Math.sign(set2_data[i].rotation) > 0 ? 1 : 0;
            // arc_flag = 1;

            // arc_path = "M" + cos1 + " " + sin1 + " A "+path_r+" "+path_r+" 1 0 " + arc_flag + " " + cos2 + " " + sin2; //Q-32,-32

            set2_data[i].rotate_text = degrees;
            set2_data[i].movement = distance;
            // set2_data[i].arc_path = arc_path;

            distance_sign = -1;
            if (degrees < -90 && degrees > -180) distance_sign = 1;
            if (degrees > 90 && degrees < 180) distance_sign = 1;
            set2_data[i].distance_sign = distance_sign;

            set1_data[i].dis_x1 = data[1][i].projected.x + Math.cos(toRadians(degrees + 90*distance_sign)) * line_distance_from_center;
            set1_data[i].dis_y1 = data[1][i].projected.y + Math.sin(toRadians(degrees + 90*distance_sign)) * line_distance_from_center;
            set1_data[i].dis_x2 = data[1][ii].projected.x + Math.cos(toRadians(degrees + 90*distance_sign)) * line_distance_from_center;
            set1_data[i].dis_y2 = data[1][ii].projected.y + Math.sin(toRadians(degrees + 90 * distance_sign)) * line_distance_from_center;



            set1_data[i].i = parseInt(i);
            set2_data[i].i = parseInt(i)+7;


        }


        /* ----------- GRID ----------- */

        var xGrid = svgContainer.selectAll('path.grid').data(data[0], key);

        xGrid
            .enter()
            .append('path')
            .attr('class', '_3d grid')
            .merge(xGrid)
            .attr('stroke', 'black')
            .attr('stroke-width', 0.3)
            .attr('fill', function(d){ return d.ccw ? 'lightgrey' : '#eeeeee'; })
            .attr('fill-opacity', 0.1)
            .attr('d', grid3d.draw);

        xGrid.exit().remove();

    /* ----------- POINTS ----------- */

        if (initial) {
            sets = svgContainer.selectAll('test').data(data[1], key);
            g_sets = sets
                .enter()
                .append('g')
                .classed("_3d", true)
                .classed("tm7", true)
                .classed("set1", function (d) { return d.color == "set1"; })
                .classed("set2", function (d) { return d.color == "set2"; })
                .attr('id', function (d) { return d.label+"_"+d.color })
                .attr('opacity', 1)
                .attr("transform", function (d) { return "translate(" + d.projected.x + "," + d.projected.y + ")" })


            circles = g_sets.append("circle")
                .attr("r", 32)
                .attr('id', function (d) { return d.label+"_"+d.color })
                .attr('opacity', function (d) { return d.color == "set1" ? 1 : 1; })
                .attr('fill', function (d) { return set_colors[d.color]; })
                .style("stroke-width", "1px")
                .attr("stroke", "#555")
                .style("stroke-dasharray", function (d,i) { return d.color=="set1" ? 4 : 0; })
                .attr('filter', 'url(#lightMe)')

            texts = g_sets.append("text")
                .attr("class","sets_labels")
                .attr("dy", ".35em")
                .attr("text-anchor", 'middle')
                .text(function (d) { return d.label })
        } else {
            g_sets.attr("transform", function (d, i) { dd = data[1][i]; return "translate(" + dd.projected.x + "," + dd.projected.y + ")";})
        }

    /* ----------- distance-lines ----------- */
        line_distance_from_center = 35;
        if (initial) {
            alllines = svgContainer.append("g").selectAll(".distance_line").data(set1_data).enter().append("g");

            lines = alllines.append("line")
                .attr("class", "distance_line _3d")
                .attr("x1", function (d, i) { return d.dis_x1 })
                .attr("y1", function (d, i) { return d.dis_y1 })
                .attr("x2", function (d, i) { return d.dis_x2 })
                .attr("y2", function (d, i) { return d.dis_y2 })
                .style("stroke", "grey")
                .attr("stroke-width", 3)
                ;
        } else {
            lines
                .attr("class", "distance_line")
                .attr("x1", function (d, i) { return d.dis_x1 })
                .attr("y1", function (d, i) { return d.dis_y1 })
                .attr("x2", function (d, i) { return d.dis_x2 })
                .attr("y2", function (d, i) { return d.dis_y2 })

        }


        /* ----------- y-Scale ----------- */

        // var yScale = svg.selectAll('path.yScale').data(data[2]);

        // yScale
        //     .enter()
        //     .append('path')
        //     .attr('class', '_3d yScale')
        //     .merge(yScale)
        //     .attr('stroke', 'black')
        //     .attr('stroke-width', .5)
        //     .attr('d', yScale3d.draw);

        // yScale.exit().remove();

        // /* ----------- y-Scale Text ----------- */

        // var yText = svg.selectAll('text.yText').data(data[2][0]);

        // yText
        //     .enter()
        //     .append('text')
        //     .attr('class', '_3d yText')
        //     .attr('dx', '.3em')
        //     .merge(yText)
        //     .each(function(d){
        //         d.centroid = {x: d.rotated.x, y: d.rotated.y, z: d.rotated.z};
        //     })
        //     .attr('x', function(d){ return d.projected.x; })
        //     .attr('y', function(d){ return d.projected.y; })
        //     .text(function (d) { return d[2] <= 0 ? d[2] : ''; });

        // yText.exit().remove();

        svgContainer.selectAll('._3d').sort(d3v4._3d().sort);
    }

    function posPointX(d){
        return d.projected.x;
    }

    function posPointY(d){
        return d.projected.y;
    }

	function init(){
        var cnt = 0;
        xGrid = [], scatter = [], yLine = [];

        for (var i of Object.keys(set1_data)) {
            n = set1_data[i];
            // xGrid.push([n.x, 1, n.z]);
            scatter.push({x: n.x, y: n.y, z: n.z, id: n.label, color: 'set1'});
        }
        for(var z = -j; z < j; z++){
            for(var x = -j; x < j; x++){
                xGrid.push([x, z, -10]);
            }
        }
        d3v4.range(-1, 11, 1).forEach(function(d){ yLine.push([-j, -j, -d]); });

        var data = [
            grid3d(xGrid),
            point3d([...set1_data, ...set2_data]),
            yScale3d([yLine])
        ];
        processData(data, 1000, true);
    }

    function dragStart(){
        mx = d3v4.event.x;
        my = d3v4.event.y;
    }

    function dragged(){
        mouseX = mouseX || 0;
        mouseY = mouseY || 0;
        beta   = (d3v4.event.x - mx + mouseX) * Math.PI / 230 ;
        alpha  = (d3v4.event.y - my + mouseY) * Math.PI / 230  * (-1);
        var data = [
            grid3d.rotateY(beta + startAngle).rotateX(alpha - startAngle)(xGrid),
            point3d.rotateY(beta + startAngle).rotateX(alpha - startAngle)([...set1_data, ...set2_data]),
            yScale3d.rotateY(beta + startAngle).rotateX(alpha - startAngle)([yLine]),
        ];
        processData(data, 0);
    }

    function dragEnd(){
        mouseX = d3v4.event.x - mx + mouseX;
        mouseY = d3v4.event.y - my + mouseY;
    }

    init();

    // animate_movement();
    create_overlay();
    return true;


    // Scaling factor from Å distance to pixel distance..
    var scaling_factor = 8;

    // Set colors
    set_1_color = "#99ccff";
    set_2_color = "#ff9999";

    // Grab plotting data from call
    var set1_data = plot_data["coordinates_set1"]
    var set2_data = plot_data["coordinates_set2"]

    // Correct first ref2 compared to ref1
    /*desired_distance = get_desired_distance(ref[0], ref[1]); //matrix_set1[ref[0]][ref[1]] * scaling_factor;
    current_distance = get_distance(set2_data[ref[0]], set2_data[ref[1]]);
    move_to_fit_distance(set2_data[ref[0]], set2_data[ref[1]], desired_distance)

    // Start moving remaining set2_data compared to ref1 and ref2
    for (var key of Object.keys(set2_data)) {
        if (ref.includes(parseInt(key))) continue;

        distance_1 = get_desired_distance(ref[0], key);
        distance_2 = get_desired_distance(ref[1], key);

        intersections = intersection(set2_data[ref[0]].x, set2_data[ref[0]].y, distance_1, set2_data[ref[1]].x, set2_data[ref[1]].y, distance_2);

        delta_1 = get_distance(set2_data[key], { "x": intersections[0], "y": intersections[2] })
        delta_2 = get_distance(set2_data[key], { "x": intersections[1], "y": intersections[3] })

        if (delta_1 < delta_2) {
            set2_data[key].x = Math.round(intersections[0]);
            set2_data[key].y = Math.round(intersections[2]);
        } else {
            set2_data[key].x = Math.round(intersections[1]);
            set2_data[key].y = Math.round(intersections[3]);
        }
    }
    var myJSON = JSON.stringify(set2_data);*/

    // Scaling of the coordinates
    set1_data.forEach(function (entry) {
      entry.x = entry.x*scaling_factor
      entry.y = entry.y*scaling_factor
    });

    set2_data.forEach(function (entry) {
      entry.x = entry.x*scaling_factor
      entry.y = entry.y*scaling_factor
    });

    padding = 53; //33
    min_y = Math.min.apply(Math, [...set1_data, ...set2_data].map(a => a.y)) - padding;
    max_y = Math.max.apply(Math, [...set1_data, ...set2_data].map(a => a.y)) + padding;

    min_x = Math.min.apply(Math, [...set1_data, ...set2_data].map(a => a.x)) - padding;
    max_x = Math.max.apply(Math, [...set1_data, ...set2_data].map(a => a.x)) + padding;

    circle_r = 33;
    line_widths = 3;
    path_r = circle_r + 3 + line_widths;
    line_distance_from_center = path_r;
    values_font_size = 16;
    values_font_size_hiding = 6;
    tm_font_size = 8;

    minimum_angle_to_show = 10;
    minimum_distance_to_show = 1;

    for (var i of Object.keys(set1_data)) {
        distance = get_distance(set1_data[i], set2_data[i])/scaling_factor;
        direction = get_vector(set1_data[i], set2_data[i]);
        cos_direction = direction[0] / direction[1];
        rad = Math.atan2(direction[1], direction[0]);
        degrees = -rad * (180 / Math.PI);
        degrees = degrees ? -degrees : 0;

        scaled_rotation = Math.abs(set2_data[i].rotation) < 20 ? 20 : Math.abs(set2_data[i].rotation);
        from = - scaled_rotation / 2;
        to = scaled_rotation / 2;
        cos1 = parseInt(Math.cos(toRadians(from)) * path_r);
        sin1 = parseInt(Math.sin(toRadians(from)) * path_r);

        cos2 = parseInt(Math.cos(toRadians(to)) * path_r);
        sin2 = parseInt(Math.sin(toRadians(to)) * path_r);

        arc_flag = Math.sign(set2_data[i].rotation) > 0 ? 1 : 0;
        arc_flag = 1;

        arc_path = "M" + cos1 + " " + sin1 + " A "+path_r+" "+path_r+" 1 0 " + arc_flag + " " + cos2 + " " + sin2; //Q-32,-32

        set2_data[i].rotate_text = degrees;
        set2_data[i].movement = distance;
        set2_data[i].arc_path = arc_path;

        distance_sign = -1;
        if (degrees < -90 && degrees > -180) distance_sign = 1;
        if (degrees > 90 && degrees < 180) distance_sign = 1;
        set2_data[i].distance_sign = distance_sign;

        console.log(i,degrees,set2_data[i].rotation)

    }

    $(containerSelector).html('')
    $(containerSelector).addClass("tm7Plot");
    $(containerSelector).css("position","relative");

    var svgContainer = d3v4.select(containerSelector).append("svg")
        .attr("viewBox", min_x + " " + min_y + " " + (max_x - min_x) + " " + (max_y - min_y))
        .attr("width", "100%")
        .attr("style", "height: 500px");




    var animate_run = 0;
    var repeat_animate = false;
    function animate_movement() {
        var delay = 500;
        var duration = 2000;

        // initial
        svgContainer.selectAll(".set2.tm7")
            .transition()
            .attr("transform", function (d, i) { return "translate(" + latest_data[1][d.i-7].projected.x + "," + latest_data[1][d.i-7].projected.y + ")" })
        svgContainer.selectAll(".set2.tm7").attr("opacity", 0.9);
        svgContainer.selectAll(".set2.circle").attr("fill", set_1_color);

        svgContainer.selectAll(".angles").attr("d", "M 35 0 A 35 35 1 0 1 35 0").attr("opacity", 0)
        svgContainer.selectAll(".labelText").attr("font-size",values_font_size_hiding)
        svgContainer.selectAll(".angles_text")
            .attr("transform", function (d, i) {
                var x = set1_data[i].x + Math.cos(toRadians(d.rotate_text)) * (line_distance_from_center + 2);
                var y = set1_data[i].y + Math.sin(toRadians(d.rotate_text)) * (line_distance_from_center + 2);
                return "translate(" + x + "," + y + ") rotate(" + (d.rotate_text + 90) + ")";
            })
            .attr("opacity", 0)
            .attr("font-size",values_font_size_hiding)

        svgContainer.selectAll(".distance_line")
            .attr("opacity", 0)
            .attr("x2", function (d, i) { return d.dis_x1 })
            .attr("y2", function (d, i) { return d.dis_y1 })

        svgContainer.selectAll(".distance_text")
            .attr("opacity", 0)
            .attr("transform", function (d, i) {
                var x = (set1_data[i].x) + Math.cos(toRadians(d.rotate_text + 90*d.distance_sign)) * (line_distance_from_center + line_widths*1.5);
                var y = (set1_data[i].y) + Math.sin(toRadians(d.rotate_text + 90*d.distance_sign)) * (line_distance_from_center + line_widths*1.5);
                text_rotation = d.distance_sign == -1 ? 0 : 180;
                return "translate(" + x + "," + y + ") rotate(" + (d.rotate_text + text_rotation) + ")";
            })
            .attr("font-size",values_font_size_hiding)

        svgContainer.selectAll(".set2_labels")
            .attr("transform", function (d) { return "rotate(" + (-d.rotate_text) + ")" })

        // move out to fully extended
        // https://bl.ocks.org/d3noob/1ea51d03775b9650e8dfd03474e202fe
        var n = 10002;
        var t0 = svgContainer.transition().delay(delay).ease(d3v4.easeExp).duration(duration);

        // move circles and rotate
        t0.selectAll(".set2.tm7")
            .attr("opacity", "0.9")
            .attr("transform", function (d, i) { return "translate(" + latest_data[1][d.i].projected.x + "," + latest_data[1][d.i].projected.y + ")" })
        t0.selectAll(".set2.circle").attr("fill", set_2_color);
        t0.selectAll(".set1.tm7").attr("opacity", "0.2");
         // make angle arc
        t0.selectAll(".angles")
            .attr("opacity", 1)
            .attr("d", function (d, i) { return d.arc_path })
        t0.selectAll(".labelText")
            .tween("text", function (d, ii) {
                var i = d3.interpolate(0, Math.abs(d.rotation));
                return function (t) {
                    d3.select(this).text(i(t).toFixed(0) + "°");
                };
            })
            .attr("font-size",values_font_size)

        t0.selectAll(".angles_text")
            .attr("transform", function (d, i) {
                var x = d.x + Math.cos(toRadians(d.rotate_text)) * (line_distance_from_center + 2);
                var y = d.y + Math.sin(toRadians(d.rotate_text)) * (line_distance_from_center + 2);
                return "translate(" + x + "," + y + ") rotate(" + (d.rotate_text + 90) + ")";
            })
            .tween("text", function (d, ii) {
                var i = d3.interpolate(0, Math.abs(d.rotation));
                return function (t) {
                    d3.select(this).text(i(t).toFixed(0) + "°");
                };
            })
            .attr("opacity", 1)
            .attr("font-size",values_font_size)
        // make distance lines
        t0.selectAll(".distance_line")
            .attr("marker-end", "url(#arrowhead)")
            .attr("x2", function (d, i) { return d.dis_x2 })
            .attr("y2", function (d, i) { return d.dis_y2 })
            .attr("opacity", 1)
        t0.selectAll(".distance_text")
            .attr("opacity", 1)
            .attr("transform", function (d, i) {
                var x = (d.x + set1_data[i].x) / 2 + Math.cos(toRadians(d.rotate_text + 90*d.distance_sign)) * (line_distance_from_center + line_widths*1.5);
                var y = (d.y + set1_data[i].y) / 2 + Math.sin(toRadians(d.rotate_text + 90*d.distance_sign)) * (line_distance_from_center + line_widths*1.5);
                text_rotation = d.distance_sign == -1 ? 0 : 180;
                return "translate(" + x + "," + y + ") rotate(" + (d.rotate_text + text_rotation) + ")";
            })
            .tween("text", function (d, ii) {
                var i = d3.interpolate(0, d.movement);
                return function (t) {
                    d3.select(this).text(i(t).toFixed(1) + "Å");
                };
            })
            .attr("font-size",values_font_size)
        t0.selectAll(".set2_labels")
            .attr("transform", function (d) { return "rotate(" + (-d.rotate_text) + ")" })

        if (!repeat_animate) return

        // move back to initial
        var t1 = t0.transition().delay(delay).duration(duration);

        // move circles and rotate
        t1.selectAll(".set2.tm7")
            .attr("transform", function (d, i) { return "translate(" + latest_data[1][d.i-7].projected.x + "," + latest_data[1][d.i-7].projected.y + ")" })
        t1.selectAll(".set2.circle").attr("fill", set_1_color);
        // make angle arc small again and disapear
        t1.selectAll(".angles")
            .attr("opacity", 0)
            .attr("d", "M 35 0 A 35 35 1 0 1 35 0")
        t1.selectAll(".labelText")
            .tween("text", function (d, ii) {
                var i = d3.interpolate(Math.abs(d.rotation), 0);
                return function (t) {
                    d3.select(this).text(i(t).toFixed(0) + "°");
                };
            })
            .attr("font-size",values_font_size_hiding)
        t1.selectAll(".angles_text")
            .attr("transform", function (d, i) {
                var x = set1_data[i].x + Math.cos(toRadians(d.rotate_text)) * (line_distance_from_center + 2);
                var y = set1_data[i].y + Math.sin(toRadians(d.rotate_text)) * (line_distance_from_center + 2);
                return "translate(" + x + "," + y + ") rotate(" + (d.rotate_text + 90) + ")";
            })
            .tween("text", function (d, ii) {
                var i = d3.interpolate(Math.abs(d.rotation), 0);
                return function (t) {
                    d3.select(this).text(i(t).toFixed(0) + "°");
                };
            })
            .attr("opacity", 0)
            .attr("font-size",values_font_size_hiding)
        // make distance lines
        t1.selectAll(".distance_line")
        .attr("x2", function (d, i) { return d.dis_x1 })
        .attr("y2", function (d, i) { return d.dis_y1 })
            .attr("opacity", 0)
        t1.selectAll(".distance_text")
            .attr("transform", function (d, i) {
                var x = (set1_data[i].x) + Math.cos(toRadians(d.rotate_text + 90*d.distance_sign)) * (line_distance_from_center + line_widths*1.5);
                var y = (set1_data[i].y) + Math.sin(toRadians(d.rotate_text + 90*d.distance_sign)) * (line_distance_from_center + line_widths*1.5);
                text_rotation = d.distance_sign == -1 ? 0 : 180;
                return "translate(" + x + "," + y + ") rotate(" + (d.rotate_text + text_rotation) + ")";
            })
            .tween("text", function (d, ii) {
                var i = d3.interpolate(d.movement, 0);
                return function (t) {
                    d3.select(this).text(i(t).toFixed(1) + "Å");
                };
            })
            .attr("opacity", 0)
            .attr("font-size",values_font_size_hiding)
        t1.selectAll(".set2_labels")
            .attr("transform", function (d) { return "rotate(" + (-d.rotate_text) + ")" })

        // Repeat!
        if (repeat_animate) t1.on("end", function (d, i) { return i == 0 ? animate_movement() : '' });

        animate_run += 1;
        // console.log('animate..',animate_run)
    }


    function create_overlay() {
        var newDiv = document.createElement("div");

        $(containerSelector).find(".controls-panel").remove();

        newDiv.setAttribute("class", "controls-panel");
        content = '<span class="pull-right network_controls_toggle" style="cursor: pointer;"><span class="glyphicon glyphicon-option-horizontal btn-download png"></span></span><span class="options" style="display: block; min-width: 120px;">' +
            'Animate <input id="animate" type="checkbox" '+ (repeat_animate ? 'checked' : '') +'><br>' +
            'Hide low numbers <input id="hide_low" type="checkbox" checked><br>' +
        '</span>';
        newDiv.innerHTML = content;

        $(containerSelector).prepend(newDiv);
        $(containerSelector).find(".options").toggle();

        $(containerSelector).find(".network_controls_toggle").click(function() {
            $(containerSelector).find(".options").slideToggle();
        });

        d3v4.select(containerSelector).select("#animate").on("change", function () {
            repeat_animate = d3.select(containerSelector).select("#animate").property("checked");
            animate_movement();
        });
        d3.select(containerSelector).select("#hide_low").on("change", function () {
            hide_low = d3.select(containerSelector).select("#hide_low").property("checked");
            $(containerSelector).find(".low_number").fadeToggle();
        });

    }

    function toRadians(angle) {
        return angle * (Math.PI / 180);
    }

    function get_distance(node1, node2) {
        let dx = parseInt(node2.projected.x) - parseInt(node1.projected.x);
        let dy = parseInt(node2.projected.y) - parseInt(node1.projected.y);
        return Math.sqrt(dx * dx + dy * dy);
    }

    function get_desired_distance(from, to) {
        distance = matrix_set2[from][to];
        if (parseInt(distance) == 0) distance = matrix_set2[to][from];
        return distance * scaling_factor;
    }

    function move_to_fit_distance(fixed_node, moving_node, desired_distance) {
        current_distance = get_distance(fixed_node, moving_node);
        console.log('current_distance', current_distance);
        vector = get_vector(fixed_node, moving_node);
        console.log('vector', vector);
        scaling = (desired_distance / current_distance);
        new_vector = vector.map(function (x) { return x * scaling; });
        console.log('new_vector', new_vector);
        moving_node.x = fixed_node.x + new_vector[0];
        moving_node.y = fixed_node.y + new_vector[1];

    }

    function get_vector(node1, node2) {
        return [node2.projected.x - node1.projected.x, node2.projected.y - node1.projected.y];
    }

    function move_to_fit_refs(ref1, ref2, node) {
        intersections = intersection(ref1.x, ref1.y)
    }

    function intersection(x0, y0, r0, x1, y1, r1) {
        var a, dx, dy, d, h, rx, ry;
        var x2, y2;

        /* dx and dy are the vertical and horizontal distances between
         * the circle centers.
         */
        dx = x1 - x0;
        dy = y1 - y0;

        /* Determine the straight-line distance between the centers. */
        d = Math.sqrt((dy * dy) + (dx * dx));

        /* Check for solvability. */
        if (d > (r0 + r1)) {
            /* no solution. circles do not intersect. */
            return false;
        }
        if (d < Math.abs(r0 - r1)) {
            /* no solution. one circle is contained in the other */
            return false;
        }

        /* 'point 2' is the point where the line through the circle
         * intersection points crosses the line between the circle
         * centers.
         */

        /* Determine the distance from point 0 to point 2. */
        a = ((r0 * r0) - (r1 * r1) + (d * d)) / (2.0 * d);

        /* Determine the coordinates of point 2. */
        x2 = x0 + (dx * a / d);
        y2 = y0 + (dy * a / d);

        /* Determine the distance from point 2 to either of the
         * intersection points.
         */
        h = Math.sqrt((r0 * r0) - (a * a));

        /* Now determine the offsets of the intersection points from
         * point 2.
         */
        rx = -dy * (h / d);
        ry = dx * (h / d);

        /* Determine the absolute intersection points. */
        var xi = x2 + rx;
        var xi_prime = x2 - rx;
        var yi = y2 + ry;
        var yi_prime = y2 - ry;

        return [xi, xi_prime, yi, yi_prime];
    }
}

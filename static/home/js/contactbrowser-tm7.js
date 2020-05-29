function tm7_plot(containerSelector, plot_data, viewBox) {


    // Ideas
    // https://en.wikipedia.org/wiki/Kabsch_algorithm <- get rotation matrix by using 3-4 least moving tms
    // Hint to record https://stackoverflow.com/questions/20864874/creating-animated-gif-files-out-of-d3-js-animations

    /*matrix_set1 = JSON.parse(matrix_set1);
    matrix_set2 = JSON.parse(matrix_set2);

    // Scaling factor from Å distance to pixel distance..
    var scaling_factor = 8;

    var TMS = ["TM1", "TM2", "TM3", "TM4", "TM5", "TM6", "TM7"];

    var set2_data = [
        { "label": "TM1", "x": 300, "y": 250, "rotation": 5 },
        { "label": "TM2", "x": 190, "y": 200, "rotation": 10 },
        { "label": "TM3", "x": 120, "y": 190, "rotation": 30 },
        { "label": "TM4", "x": 70, "y": 100, "rotation": -10 },
        { "label": "TM5", "x": 30, "y": 230, "rotation": 20 },
        { "label": "TM6", "x": 80, "y": 300, "rotation": -90 },
        { "label": "TM7", "x": 200, "y": 350, "rotation": 40 }];

        var set1_data = [{ "label": "TM1", "x": 300, "y": 250 },
            { "label": "TM2", "x": 215, "y": 271 },
            { "label": "TM3", "x": 106, "y": 248 },
            { "label": "TM4", "x": 200, "y": 185 },
            { "label": "TM5", "x": 101, "y": 347 },
            { "label": "TM6", "x": 190, "y": 377 },
            { "label": "TM7", "x": 276, "y": 361 }]

    var nodes_initial = [
        { "label": "TM1", "x": 300, "y": 250 },
        { "label": "TM2", "x": 230, "y": 180 },
        { "label": "TM3", "x": 150, "y": 140 },
        { "label": "TM4", "x": 70, "y": 100 },
        { "label": "TM5", "x": 30, "y": 230 },
        { "label": "TM6", "x": 80, "y": 300 },
        { "label": "TM7", "x": 200, "y": 350 }];
    */

    /*if (containerSelector=='#single_2'){
        for (var key of Object.keys(set1_data)) {
            set1_data[key].x -= 5;
        }
    }
    if (containerSelector=='#single_3'){
        for (var key of Object.keys(set1_data)) {
            set1_data[key].x -= 10;
        }
    }*/

    // Scaling factor from Å distance to pixel distance..
    var scaling_factor = 8;

    // Set colors
    set_1_color = "#99ccff";
    set_2_color = "#ff9999";

    // Grab plotting data from call
    // var set1_data = Object.assign({}, plot_data["coordinates_set1"]);
    // var set2_data = Object.assign({}, plot_data["coordinates_set2"]);
    var set1_data = JSON.parse(JSON.stringify(plot_data["coordinates_set1"]));
    var set2_data = JSON.parse(JSON.stringify(plot_data["coordinates_set2"]));

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

    console.log('viewBox', viewBox)
    viewBox_x = viewBox['diff_x']*scaling_factor+2*padding
    viewBox_y = viewBox['diff_y'] * scaling_factor + 2 * padding

    // fix min_x,min_y to reflect the viewBox from input
    min_y = min_y - (viewBox_y-(max_y-min_y))/2
    min_x = min_x - (viewBox_x-(max_x-min_x))/2

    circle_r = 23;
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
        .attr("viewBox", min_x + " " + min_y + " " + viewBox_x + " " + viewBox_y)
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


    var set1 = svgContainer.selectAll("set1")
        .data(set1_data)

    /*Create and place the "blocks" containing the circle and the text */
    var elemEnter = set1.enter()
        .append("g")
        .attr("class", "set1 tm7")
        .attr("transform", function (d) { return "translate(" + d.x + "," + d.y + ")" })
        .attr("opacity", 0)

    /*Create the circle for each block */
    var circle1 = elemEnter.append("circle")
        .attr("r", circle_r)
        .attr("stroke", "black")
        .style("stroke-dasharray","4")
        .attr("class", "set1 circle")
        .attr("fill", set_1_color);

    /* Create the text for each block */
    elemEnter.append("text")
        .attr("dy", ".35em")
        .attr("text-anchor", 'middle')
        // .text(function (d) { return d.label })


    var set2 = svgContainer.selectAll("set2")
        .data(set2_data)

    /*Create and place the "blocks" containing the circle and the text */
    var elemEnter2 = set2.enter()
        .append("g")
        .attr("class", "set2 tm7")
        .attr("transform", function (d, i) {
            return "translate(" + d.x + "," + d.y + ") rotate(" + d.rotate_text + ")"
        })
        .attr("opacity", 1)

    /*Create the circle for each block */
    var circle2 = elemEnter2.append("circle")
        .attr("r", circle_r)
        .attr("class", "set2 circle")
        .style("stroke-width", "1px")
        .attr("stroke", "#555")
        .attr("fill", set_2_color);

        // Append images
    // var images = elemEnter2.append("svg:image")
    // .attr("xlink:href",  "http://localhost:8010/static/home/images/helix.png")
    // .attr("x", function(d) { return -30;})
    // .attr("y", function(d) { return -30;})
    // .attr("height", 60)
    // .attr("width", 60);

    /* Create the text for each block */
    text2 = elemEnter2.append("text")
        .attr("class","set2_labels")
        .attr("dy", ".35em")
        .attr("text-anchor", 'middle')
        .attr("transform", function (d, i) {
            return "rotate(" + (-d.rotate_text) + ")"
        })
        .text(function (d) { return d.label })

    angles = elemEnter2.append("path")
        .attr("class", "angles")
        .classed("low_number", function (d, i) { return (Math.abs(d.rotation)<minimum_angle_to_show);})
        .attr("marker-end", function (d, i) { return Math.sign(d.rotation) > 0 ? "url(#arrowhead)" : "" })
        .attr("marker-start", function (d, i) { return Math.sign(d.rotation) > 0 ? "" : "url(#arrowhead-rev)" })
        .attr("d", function (d, i) {
            return d.arc_path
        })
        .attr("stroke", "grey")
        .attr("stroke-width", line_widths)
        .attr("fill", "transparent")
        .attr("id", function (d, i) { return containerSelector + "linkId_" + i; })
        .attr("display", function (d, i) {
            return Math.abs(d.rotation)>=minimum_angle_to_show ? "" : "none";
        });

    var labelParent = elemEnter2.append("text")
        .attr("class", "labelParent")
        .attr("dx", 0)
        .attr("dy", -4)
        .style("fill", "black")
        .style("opacity", 0.5)
        .attr("id", function (d, i) { return "labelText_" + i; });

    var labelText = labelParent.append("textPath")
        .attr("class", "labelText")
        .attr("xlink:href", function (d, i) { return "#" + containerSelector + "linkId_" + i; })
        .attr("startOffset", "50%")
        .attr("text-anchor", "Middle")
        .attr("font-size", values_font_size)
        .text(function (d, i) { return d.rotation ? Math.abs(d.rotation) + "°" : "" })
        .attr("display", function (d, i) {
            return Math.abs(d.rotation) >= 40 && Math.abs(d.rotation) > 0 ? "" : "none";
        });

    alllines = svgContainer.append("g").selectAll("lines").data(set2_data).enter().append("g");
    var angles_text_alt = alllines.append("text")
        .attr("class","angles_text")
        .classed("low_number", function (d, i) { return (Math.abs(d.rotation) < minimum_angle_to_show);})
        .attr("text-anchor", "middle")
        .attr("dx", function (d) { return Math.sign(d.rotation)==1 ? 2 : 5; })
        .attr("transform", function (d, i) {
            var x = d.x + Math.cos(toRadians(d.rotate_text)) * (line_distance_from_center + 2);
            var y = d.y + Math.sin(toRadians(d.rotate_text)) * (line_distance_from_center + 2);
            return "translate(" + x + "," + y + ") rotate(" + (d.rotate_text + 90) + ")";
        })
        .attr("font-size", values_font_size)
        .attr("fill", "grey")
        .text(function (d, i) { return d.rotation ? Math.abs(d.rotation) + "°" : "" })
        .attr("display", function (d, i) {
            return Math.abs(d.rotation) < 40 && Math.abs(d.rotation) >= minimum_angle_to_show  ? "" : "none";
        });


    lines = alllines.append("line")
        .attr("class","distance_line")
        .classed("low_number", function (d, i) { return d.movement <= minimum_distance_to_show;})
        .attr("x1", function (d, i) { return set1_data[i].x + Math.cos(toRadians(d.rotate_text + 90*d.distance_sign)) * line_distance_from_center })
        .attr("y1", function (d, i) { return set1_data[i].y + Math.sin(toRadians(d.rotate_text + 90*d.distance_sign)) * line_distance_from_center })
        .attr("x2", function (d) { return d.x + Math.cos(toRadians(d.rotate_text + 90*d.distance_sign)) * line_distance_from_center })
        .attr("y2", function (d) { return d.y + Math.sin(toRadians(d.rotate_text + 90*d.distance_sign)) * line_distance_from_center })
        .style("stroke", "grey")
        .attr("stroke-width", line_widths)
        .attr("display", function (d, i) {
            return d.movement > minimum_distance_to_show ? "" : "none";
        });

    lines_text = alllines.append("text")
        .attr("class","distance_text")
        .classed("low_number", function (d, i) { return d.movement <= minimum_distance_to_show;})
        .attr("text-anchor", "middle")
        .attr("transform", function (d, i) {
            var x = (d.x + set1_data[i].x) / 2 + Math.cos(toRadians(d.rotate_text + 90*d.distance_sign)) * (line_distance_from_center + line_widths*1.5);
            var y = (d.y + set1_data[i].y) / 2 + Math.sin(toRadians(d.rotate_text + 90*d.distance_sign)) * (line_distance_from_center + line_widths*1.5);
            return "translate(" + x + "," + y + ") rotate(" + (d.rotate_text + 180) + ")";
        })
        .attr("font-size", values_font_size)
        .attr("fill", "grey")
        .text(function (d, i) {
            return d.movement > 0 ? d.movement.toFixed(1) + "Å" : "";
        })
        .attr("display", function (d, i) {
            return d.movement > minimum_distance_to_show ? "" : "none";
        });

    d3.selectAll(".set2.circle")
        .style("filter", "url(#glow)");
    // d3.selectAll(".set1.circle")
    //     .style("filter", "url(#shadow)");

    var animate_run = 0;
    var repeat_animate = false;
    function animate_movement() {
        var delay = 500;
        var duration = 2000;

        // initial
        svgContainer.selectAll(".set2.tm7")
            .transition()
            .attr("transform", function (d, i) { return "translate(" + set1_data[i].x + "," + set1_data[i].y + ") rotate(" + (d.rotate_text - d.rotation) + ")" })
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
            .attr("x2", function (d, i) { return set1_data[i].x + Math.cos(toRadians(d.rotate_text + 90*d.distance_sign)) * line_distance_from_center })
            .attr("y2", function (d, i) { return set1_data[i].y + Math.sin(toRadians(d.rotate_text + 90*d.distance_sign)) * line_distance_from_center })

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
            .attr("transform", function (d) { return "translate(" + d.x + "," + d.y + ") rotate(" + (d.rotate_text) + ")" })
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
            .attr("x2", function (d) { return d.x + Math.cos(toRadians(d.rotate_text + 90*d.distance_sign)) * line_distance_from_center })
            .attr("y2", function (d) { return d.y + Math.sin(toRadians(d.rotate_text + 90*d.distance_sign)) * line_distance_from_center })
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
            .attr("transform", function (d, i) { return "translate(" + set1_data[i].x + "," + set1_data[i].y + ") rotate(" + (d.rotate_text - d.rotation) + ")" });
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
            .attr("x2", function (d, i) { return set1_data[i].x + Math.cos(toRadians(d.rotate_text + 90*d.distance_sign)) * line_distance_from_center })
            .attr("y2", function (d, i) { return set1_data[i].y + Math.sin(toRadians(d.rotate_text + 90*d.distance_sign)) * line_distance_from_center })
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

    animate_movement();
    create_overlay();

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
        let dx = parseInt(node2.x) - parseInt(node1.x);
        let dy = parseInt(node2.y) - parseInt(node1.y);
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
        return [node2.x - node1.x, node2.y - node1.y];
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

function createScatterplot(data,containerSelector) {

    containerSelector_hash = "#" + containerSelector;

    $(containerSelector_hash).html('')
    $(containerSelector_hash).addClass("scatterplot");
    $(containerSelector_hash).css("position","relative");

    // var svgContainer = d3v4.select(containerSelector).append("svg")
    //     .attr("viewBox", min_x + " " + min_y + " " + (max_x - min_x) + " " + (max_y - min_y))
    //     .attr("width", "100%")
    //     .attr("style", "height: 500px");
    
    // console.log('making snakeplot!')
    // console.log(containerSelector)
    // console.log(data['snakeplot'])
    var colors = {}
    colors['distance'] = {}
    colors['distance_abs'] = {}


    index_names = { 0: 'core_distance', 1: 'a_angle', 2: 'outer_angle', 3: 'tau', 4: 'phi', 5: 'psi', 6: 'sasa', 7: 'rsa', 8: 'theta', 9: 'hse', 10: 'tau_angle' }
    neg_and_positives = ['core_distance','sasa','rsa', 'hse']
    nice_index_names = {
        'a_angle' : 'Angle to helix&7TM axes',
        'outer_angle': 'Rotamer',
        'tau': 'Tau',
        'phi': 'Phi',
        'psi': 'Psi',
        'sasa': 'SASA',
        'sasa_abs': 'SASA (abs)',
        'rsa': 'RSA',
        'rsa_abs': 'RSA (abs)',
        'theta': 'Theta',
        'hse': 'HSE',
        'hse_abs': 'HSE (abs)',
        'tau_angle': 'Tau dihedral',
        'distance': 'Distance',
        'distance_abs': 'Distance (abs)',
        'core_distance': 'Distance to 7TM axis',
        'core_distance_abs': 'Distance to 7TM axis (abs)',  
    }

    $.each(data['distances'], function (gn, dis) {
        seq_pos = gn;
        // console.log(gn, dis['avg'], seq_pos);
        value = dis['avg'];
        scale = Math.abs(value) / data['ngl_max_diff_distance'];
        if (value < 0) {
            // if the header is a set two, then make it red
            scale_abs = (1-scale)*0.5;
        } else if (value >= 0) {
            // Positive numbers are blue either cos they are set 1 or cos "set 1 has most"
            // This is also used for single set/structure
            scale_abs = (scale)*0.5+0.5;
        }
        colors['distance'][seq_pos] = [value,scale_abs,data['ngl_max_diff_distance']];
        colors['distance_abs'][seq_pos] = [Math.abs(value), scale, data['ngl_max_diff_distance']];
        
    });
    // console.log(colors)
    // get maximum values
    var max_values = {}
    $.each(data['tab4'], function (gn, v) {
        $.each(v["angles"], function (i, a) {
            if (!(i in max_values)) {
                max_values[i] = 0;
                colors[index_names[i]] = {};
                if (neg_and_positives.includes(index_names[i])) {
                    colors[index_names[i]+"_abs"] = {};
                }
            }
            if (Math.abs(a[0])>max_values[i]) max_values[i] = Math.abs(a[0])
        });
    });
    $.each(data['tab4'], function (gn, v) {
        $.each(v["angles"], function (i, a) {
            seq_pos = gn
            value = a[0];
            scale = Math.abs(value) / max_values[i];
            neg_and_positive = false;
            if (neg_and_positives.includes(index_names[i])) {
                neg_and_positive = true;
            }
            if (value < 0) {
                // if the header is a set two, then make it red
                scale_abs = (1-scale)*0.5;
            } else if (value >= 0) {
                // Positive numbers are blue either cos they are set 1 or cos "set 1 has most"
                // This is also used for single set/structure
                scale_abs = (scale)*0.5+0.5;
            }
            if (neg_and_positive) {
                colors[index_names[i]][seq_pos] = [value,scale_abs,max_values[i]];
                colors[index_names[i] + '_abs'][seq_pos] = [Math.abs(value), scale, max_values[i]];              
            } else {
                colors[index_names[i]][seq_pos] = [value, scale, max_values[i]];
            }
            if (!(i in max_values)) max_values[i] = 0;
        });
    });

    // var colors = {}
    // colors['distance'] = {}

    // $.each(data['distances'], function (gn, dis) {
    //     seq_pos = data['snakeplot_lookup'][gn];
    //     // console.log(gn, dis['avg'], seq_pos);
    //     value = dis['avg'];
    //     scale = Math.abs(value) / data['ngl_max_diff_distance'];
    //     var color = { r: 255, g: 255, b: 255 };
    //     var color2 = { r: 255, g: 255, b: 255 };
    //     if (value < 0) {
    //         // if the header is a set two, then make it red
    //         color = { r: 255, g: 255-(255-153)*scale, b: 255-(255-153)*scale }; //red
    //         color2 = { r: 255-(255-128)*scale, g: 255-(255)*scale, b: 255-(255-128)*scale }; //purple
    //     } else if (value > 0) {
    //         // Positive numbers are blue either cos they are set 1 or cos "set 1 has most"
    //         // This is also used for single set/structure
    //         color = { r: 255-(255-153)*scale, g: 255-(255-204)*scale, b: 255 }; //blue
    //         color2 = { r: 255-(255)*scale, b: 255-(255)*scale, g: 128 }; //green
    //     }
    //     var hex = rgb2hex(color.r, color.g, color.b);
    //     var hex2 = rgb2hex(color2.r, color2.g, color2.b);
    //     // grey
    //     var color_grey = { r: 255*(1-scale), g: 255*(1-scale), b: 255*(1-scale) };
    //     var hex_grey = rgb2hex(color_grey.r, color_grey.g, color_grey.b);
    //     colors['distance'][seq_pos] = [hex,value,scale,hex_grey,hex2];
    // });

    // index_names = { 0: 'core_distance', 1: 'a_angle', 2: 'outer_angle', 3: 'tau', 4: 'phi', 5: 'psi', 6: 'sasa', 7: 'rsa', 8: 'theta', 9: 'hse', 10: 'tau_angle' }
    // // get maximum values
    // var max_values = {}
    // $.each(data['tab4'], function (gn, v) {
    //     $.each(v["angles"], function (i, a) {
    //         if (!(i in max_values)) {
    //             max_values[i] = 0;
    //             colors[index_names[i]] = {};
    //         }
    //         if (Math.abs(a[0])>max_values[i]) max_values[i] = Math.abs(a[0])
    //         // console.log(gn, i, a);
    //     });
    // });
    // $.each(data['tab4'], function (gn, v) {
    //     $.each(v["angles"], function (i, a) {
    //         seq_pos = data['snakeplot_lookup'][gn];
    //         value = a[0];
    //         scale = Math.abs(value) / max_values[i];
    //         var color = { r: 255, g: 255, b: 255 };
    //         var color2 = { r: 255, g: 255, b: 255 };
    //         if (value < 0) {
    //             // if the header is a set two, then make it red
    //             color = { r: 255, g: 255-(255-153)*scale, b: 255-(255-153)*scale }; //red
    //             color2 = { r: 255-(255-128)*scale, g: 255-(255)*scale, b: 255-(255-128)*scale }; //purple
    //         } else if (value > 0) {
    //             // Positive numbers are blue either cos they are set 1 or cos "set 1 has most"
    //             // This is also used for single set/structure
    //             color = { r: 255-(255-153)*scale, g: 255-(255-204)*scale, b: 255 }; //blue
    //             color2 = { r: 255-(255-153)*scale, b: 255-(255-204)*scale, g: 255 }; //green
    //         }
    //         var hex = rgb2hex(color.r, color.g, color.b);
    //         var hex2 = rgb2hex(color2.r, color2.g, color2.b);
    //         // grey
    //         var color_grey = { r: 255*(1-scale), g: 255*(1-scale), b: 255*(1-scale) }; 
    //         var hex_grey = rgb2hex(color_grey.r, color_grey.g, color_grey.b);
    //         colors[index_names[i]][seq_pos] = [hex,value,scale,hex_grey,hex2];
    //         if (!(i in max_values)) max_values[i] = 0;
    //         // console.log(gn, i, a);
    //     });
    // });
    // console.log(colors);

    
    function get_values(x, y) {
        index_names = { 0: 'core_distance', 1: 'a_angle', 2: 'outer_angle', 3: 'tau', 4: 'phi', 5: 'psi', 6: 'sasa', 7: 'rsa', 8: 'theta', 9: 'hse', 10: 'tau_angle' }
        const index_names_rev = {};
        Object.keys(index_names).forEach(key => {
            index_names_rev[index_names[key]] = key;
        });
        console.log(index_names_rev,index_names)
        
        var X = [];
        var Y = [];
        var names = [];
        $.each(filtered_gns, function (i, gn) {
            // console.log('filtered', gn,x,index_names_rev[x]);
            // console.log('filtered', gn,y,index_names_rev[y]);
            if (x == 'distance') {
                if (!(gn in data['distances'])) return true;
                X.push(data['distances'][gn]['avg']);
            } else {
                X.push(data['tab4'][gn]['angles'][index_names_rev[x]][0]);
            }
            if (y == 'distance') {
                if (!(gn in data['distances'])) return true;
                Y.push(data['distances'][gn]['avg']);
            } else {
                Y.push(data['tab4'][gn]['angles'][index_names_rev[y]][0]);
            }
            names.push(gn);
        });
        return [X, Y, names];
    }

    index_names = { 0: 'core_distance', 1: 'a_angle', 2: 'outer_angle', 3: 'tau', 4: 'phi', 5: 'psi', 6: 'sasa', 7: 'rsa', 8: 'theta', 9: 'hse', 10: 'tau_angle' }
    neg_and_positives = ['core_distance','sasa','rsa', 'hse']
    nice_index_names = {
        'a_angle' : 'Angle to helix/7TM axes',
        'outer_angle': 'Rotamer',
        'tau': 'Tau',
        'phi': 'Phi',
        'psi': 'Psi',
        'sasa': 'SASA',
        'sasa_abs': 'SASA (abs)',
        'rsa': 'RSA',
        'rsa_abs': 'RSA (abs)',
        'theta': 'Theta',
        'hse': 'HSE',
        'hse_abs': 'HSE (abs)',
        'tau_angle': 'Tau dihedral',
        'distance': 'Distance',
        'distance_abs': 'Distance (abs)',
        'core_distance': 'Distance to 7TM axis',
        'core_distance_abs': 'Distance to 7TM axis (abs)',  
    }

    var margin = { top: 10, right: 30, bottom: 40, left: 50 },
    width = 520 - margin.left - margin.right,
    height = 520 - margin.top - margin.bottom;
    var scatter_svg = d3v4.select(containerSelector_hash)
        .append("svg")
        .attr("viewBox", 0 + " " + 0 + " " + (width+ margin.left + margin.right) + " " + (height+ margin.top + margin.bottom))
        .attr("width", "100%")
        .attr("style", "height: 500px")
        // .attr("width", width + margin.left + margin.right)
        // .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("id","scatter_id")
        .attr("transform",
            "translate(" + margin.left + "," + margin.top + ")")
    
    function draw_plot() {
        x_axis_type = $(containerSelector_hash + " #change_x").val();
        y_axis_type = $(containerSelector_hash + " #change_y").val();
        sp_color = $(containerSelector_hash + " #sp_color").val();
        sp_size = $(containerSelector_hash + " #sp_size").val();
        console.log('draw scatter', x_axis_type, y_axis_type, sp_color, sp_size);
        plot_data = get_values(x_axis_type, y_axis_type)
        // set the dimensions and margins of the graph
        scatter_svg.html("");
        // append the svg object to the body of the page
        // .attr("width", "100%")
        // .attr("style", "height: 500px");
    


        // Add the grey background that makes ggplot2 famous
        scatter_svg
            .append("rect")
            .attr("x", 0)
            .attr("y", 0)
            .attr("height", height)
            .attr("width", height)
            // .style("fill", "EBEBEB")
            .style("fill", "white")
    
        console.log(plot_data);

        // Add X axis
        var x = d3v4.scaleLinear()
            .domain([Math.min(...plot_data[0]), Math.max(...plot_data[0])])
            .range([0, width])
        scatter_svg.append("g")
            .attr("transform", "translate(0," + height + ")")
            .call(d3v4.axisBottom(x).tickSize(-height * 1.3).ticks(10))
            .select(".domain").remove()

        // Add Y axis
        var y = d3v4.scaleLinear()
            .domain([Math.min(...plot_data[1]), Math.max(...plot_data[1])])
            .range([height, 0])
            .nice()
        scatter_svg.append("g")
            .call(d3v4.axisLeft(y).tickSize(-width * 1.3).ticks(7))
            .select(".domain").remove()

        // Customization
        scatter_svg.selectAll(".tick line").attr("stroke", "#EBEBEB")

        // Add X axis label:
        scatter_svg.append("text")
            .attr("text-anchor", "end")
            .attr("x", width / 2 + margin.left)
            .attr("y", height + margin.top + 20)
            .text(nice_index_names[x_axis_type]);

        // Y axis label:
        scatter_svg.append("text")
            .attr("text-anchor", "end")
            .attr("transform", "rotate(-90)")
            .attr("y", -margin.left + 20)
            .attr("x", -margin.top - height / 2 + 20)
            .text(nice_index_names[y_axis_type])

        // Color scale: give me a specie name, I return a color
        var color = d3v4.scaleOrdinal()
            .domain(["setosa", "versicolor", "virginica"])
            .range(["#F8766D", "#00BA38", "#619CFF"])

        // generate data
        var data = []
        plot_data[0].forEach(function (d_x, i) {
            d_y = plot_data[1][i]
            name = plot_data[2][i]
            if (d_y)
                data.push({ x: d_x, y: d_y, name: name });
        });
        console.log(data)
        // Add dots
        scatter_svg.append('g')
            .selectAll("dot")
            .data(data)
            .enter()
            .append("circle")
            .attr("cx", function (d) { return x(d.x); })
            .attr("cy", function (d) { return y(d.y); })
            .attr("opacity",0.7)
            .attr("r", function (d) { return d.name in colors[sp_size] ? 2+colors[sp_size][d.name][1] * 8 : 5;})
            .style("fill", function (d) { return d.name in colors[sp_size] ? color_by_scale(colors[sp_color][d.name][1], "red", "blue") : "black";} )
    
    scatter_svg.append("g")
        .attr("font-family", "sans-serif")
        .attr("font-size", 10)
        .selectAll("text")
        .data(data)
        .join("text")
        .attr("dy", "0.35em")
        .attr("x", d => x(d.x) + 7)
        .attr("y", d => y(d.y))
        .text(d => d.name);
    

}


   
    nice_index_names = {
        'a_angle' : 'Angle to helix&7TM axes',
        'outer_angle': 'Rotamer',
        'tau': 'Tau',
        'phi': 'Phi',
        'psi': 'Psi',
        'sasa': 'SASA',
        'sasa_abs': 'SASA (abs)',
        'rsa': 'RSA',
        'rsa_abs': 'RSA (abs)',
        'theta': 'Theta',
        'hse': 'HSE',
        'hse_abs': 'HSE (abs)',
        'tau_angle': 'Tau dihedral',
        'distance': 'Distance',
        'distance_abs': 'Distance (abs)',
        'core_distance': 'Distance to 7TM axis',
        'core_distance_abs': 'Distance to 7TM axis (abs)',  
    }
    var select_data_options = ''
    $.each(nice_index_names, function (key, description) {
        select_data_options += '<option value="' + key + '">' + description + '</option>';
    });

    function create_overlay() {
        var newDiv = document.createElement("div");

        $(containerSelector_hash).find(".controls-panel").remove();
        newDiv.setAttribute("class", "controls-panel");
        content = '<span class="pull-right snakeplot_controls_toggle" style="cursor: pointer;"><span class="glyphicon glyphicon-option-horizontal btn-download png"></span></span><span class="options" style="display: block; min-width: 120px;">';
            // 'Generic number with AA<input id="generic" type="checkbox"><br>' +
            // 'Only plot kept positions <input id="color_filtered" type="checkbox" checked><br>';
        
        index_names = { 0: 'core_distance', 1: 'a_angle', 2: 'outer_angle', 3: 'tau', 4: 'phi', 5: 'psi', 6: 'sasa', 7: 'rsa', 8: 'theta', 9: 'hse', 10: 'tau_angle' }
    
        content += '<table><tr><th>X-axis</th><th>Y-axis</th></tr><tr><td><select id="change_x" class="change_axis">' +
            select_data_options +
            '</select></td>'
            ;
        content += '<td><select id="change_y" class="change_axis">' +
            select_data_options +
            '</select></td></tr>';
            
        content += '<tr><th>Color</th><th>Size</th></tr><tr>';
        content += '<td><select id="sp_color" class="change_axis">' +
            select_data_options +
            '</select></td>';
        content += '<td><select id="sp_size" class="change_axis">' +
            select_data_options +
            '</select></td>';
        content += '</table > '
                ;
        content += '</span>';
        newDiv.innerHTML = content;

        $(containerSelector_hash).prepend(newDiv);
        $(containerSelector_hash).find(".options").toggle();

        $(containerSelector_hash).find(".snakeplot_controls_toggle").click(function() {
            $(containerSelector_hash).find(".options").slideToggle();
        });

        d3v4.select(containerSelector_hash).select("#generic").on("change", function () {
            generic = d3.select(containerSelector_hash).select("#generic").property("checked");
            $(containerSelector_hash).find('.rtext').each(function () {
                original_title = $(this).attr('original_title');
                if (original_title.split(" ")[1].length) {
                    // this has GN
                    gn = original_title.split(" ")[1].split("x")[1];
                    AA = original_title[0];
                } else {
                    gn = "-";
                    AA = original_title[0];
                }
                $(this).attr("font-size", generic ? 12 : 16);
                $(this).html(generic ? "<tspan dy=-6>" +AA + "</tspan><tspan dy=9 dx=-10>" + gn +"</tspan>" : AA);
            });
        });



        $(containerSelector_hash + " .change_axis").on("change", function () {
            console.log('change axis!');
            draw_plot();
        });

        

    }


    create_overlay();
    draw_plot();
    
}

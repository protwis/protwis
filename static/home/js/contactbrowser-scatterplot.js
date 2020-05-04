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
    colors['segment'] = {}
    colors['network'] = {}
    colors['set_presence'] = {}


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
        'set_presence' : 'Set specific presense'
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
    console.log(colors)

    $.each(filtered_gns_presence, function (gn, value) {
        seq_pos = gn;
        scale = value / 1
        colors['set_presence'][seq_pos] = [value,scale,1];
    })

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
    var highlight = ['TM1', 'TM2', 'TM3', 'TM4', 'TM5', 'TM6', 'TM7', 'H8'];
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
        var seg = data['segm_lookup'][gn];
        colors["segment"][seq_pos] = rb_colors[highlight.indexOf(seg)];
        
    });
    //console.log('colors',colors)
    
    function get_values(x, y) {
        index_names = { 0: 'core_distance', 1: 'a_angle', 2: 'outer_angle', 3: 'tau', 4: 'phi', 5: 'psi', 6: 'sasa', 7: 'rsa', 8: 'theta', 9: 'hse', 10: 'tau_angle' }
        const index_names_rev = {};
        Object.keys(index_names).forEach(key => {
            index_names_rev[index_names[key]] = key;
        });
        //console.log(index_names_rev, index_names)
        var x_abs = false;
        var y_abs = false;
        if (x.includes("_abs")) {
            x = x.slice(0, -4);
            x_abs = true;
        }
        if (y.includes("_abs")) {
            y = y.slice(0, -4);
            y_abs = true;
        }
        
        var X = [];
        var Y = [];
        var names = [];
        $.each(filtered_gns, function (i, gn) {
            // console.log('filtered', gn,x,index_names_rev[x]);
            // console.log('filtered', gn,y,index_names_rev[y]);
            if (x == 'distance') {
                if (!(gn in data['distances'])) return true;
                x_abs ? X.push(Math.abs(data['distances'][gn]['avg'])) : X.push(data['distances'][gn]['avg']);
            } else if (x == 'set_presence') {
                if (!(gn in filtered_gns_presence)) return true;
                x_abs ? X.push(Math.abs(filtered_gns_presence[gn])) : X.push(filtered_gns_presence[gn]);
            } else {
                x_abs ? X.push(Math.abs(data['tab4'][gn]['angles'][index_names_rev[x]][0])) : X.push(data['tab4'][gn]['angles'][index_names_rev[x]][0]);
            }
            if (y == 'distance') {
                if (!(gn in data['distances'])) return true;
                //Y.push(data['distances'][gn]['avg']);
                y_abs ? Y.push(Math.abs(data['distances'][gn]['avg'])) : Y.push(data['distances'][gn]['avg']);
            } else if (y == 'set_presence') {
                if (!(gn in filtered_gns_presence)) return true;
                y_abs ? Y.push(Math.abs(filtered_gns_presence[gn]['avg'])) : Y.push(filtered_gns_presence[gn]['avg']);
            } else {
                //Y.push(data['tab4'][gn]['angles'][index_names_rev[y]][0]);
                y_abs ? Y.push(Math.abs(data['tab4'][gn]['angles'][index_names_rev[y]][0])) : Y.push(data['tab4'][gn]['angles'][index_names_rev[y]][0]);
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
        'set_presence' : 'Set specific presense'
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
        label_only = $(containerSelector_hash + " #label_only").prop("checked");
        label = $(containerSelector_hash + " #change_label").val();


        color_id1  = $(containerSelector_hash+" #color1").val();
        color_id2  = $(containerSelector_hash+" #color2").val();
        color_id3  = $(containerSelector_hash+" #color3").val();

        
        console.log('draw scatter', x_axis_type, y_axis_type, sp_color, sp_size,'label_only',label_only);
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
            .attr("width", width)
            // .style("fill", "EBEBEB")
            .style("fill", "white")
    
        //console.log(plot_data);

        // Add X axis
        var x = d3v4.scaleLinear()
            .domain([Math.min(...plot_data[0]), Math.max(...plot_data[0])])
            .range([0, width])
        scatter_svg.append("g")
            .attr("transform", "translate(0," + height + ")")
            .call(d3v4.axisBottom(x).tickSize(-height).ticks(10))
            .select(".domain").remove()

        // Add Y axis
        var y = d3v4.scaleLinear()
            .domain([Math.min(...plot_data[1]), Math.max(...plot_data[1])])
            .range([height, 0])
            .nice()
        scatter_svg.append("g")
            .call(d3v4.axisLeft(y).tickSize(-width).ticks(10))
            .select(".domain").remove()

        // Customization
        scatter_svg.selectAll(".tick line").attr("stroke", "#EBEBEB")

        // Add X axis label:
        scatter_svg.append("text")
            .attr("text-anchor", "middle")
            .attr("x", width / 2)
            .attr("y", height + margin.top + 20)
            .text(nice_index_names[x_axis_type]);

        // Y axis label:
        scatter_svg.append("text")
            .attr("text-anchor", "middle")
            .attr("transform", "rotate(-90)")
            .attr("y", -margin.left + 20)
            .attr("x", -margin.top - height / 2 + 20)
            .text(nice_index_names[y_axis_type])

        // Color scale: give me a specie name, I return a color
        // var color = d3v4.scaleOrdinal()
        //     .domain(["setosa", "versicolor", "virginica"])
        //     .range(["#F8766D", "#00BA38", "#619CFF"])

        // generate data
        var scatter_data = []
        plot_data[0].forEach(function (d_x, i) {
            d_y = plot_data[1][i]
            name = plot_data[2][i]
            // console.log(name, data['tab4']);
            // console.log(name,data['tab4'][name]['set1_seq_cons'],data['tab4'][name]['set2_seq_cons'],data['tab4'][name]['all_seq_cons'],data['snakeplot_lookup_aa'][name])
            if (d_y)
            scatter_data.push({ x: d_x, y: d_y, name: name, set1_cons : data['tab4'][name]['set1_seq_cons'][0]+name, set2_cons : data['tab4'][name]['set2_seq_cons'][0]+name, comb_cons: data['tab4'][name]['all_seq_cons'][0]+name});
        });
        //console.log(scatter_data)
        // Add dots
        label_offset = 0;
        text_anchor = 'middle';
        font_weight = 'bold';
        if (!label_only) {
            label_offset = 7;
            text_anchor = 'left';
            font_weight = 'normal';
            scatter_svg.append('g')
                .selectAll("dot")
                .data(scatter_data)
                .enter()
                .append("circle")
                .attr("cx", function (d) { return x(d.x); })
                .attr("cy", function (d) { return y(d.y); })
                .attr("opacity",0.7)
                .attr("r", function (d) { return sp_size!='none' && d.name in colors[sp_size] ? 2+colors[sp_size][d.name][1] * 8 : 5;})
                .style("fill", function (d) {
                    if (sp_color == "segment") {
                        return colors[sp_color][d.name];
                    } else if (sp_color == "network") {
                        network_group = filtered_cluster_groups.filter(l => l.includes(d.name));
                        network_group_id = filtered_cluster_groups.indexOf(network_group[0]);
                        if (filtered_cluster_groups.length > 10 || 1==1) {
                            return d3v4.interpolateTurbo(network_group_id/filtered_cluster_groups.length);
                        } else {
                            return d3v4.schemePaired[network_group_id]
                        }
                    } else {
                        console.log(sp_color);
                        console.log(colors[sp_color]);
                        return d.name in colors[sp_color] ? color_by_scale(colors[sp_color][d.name][1], color_id1, color_id2, color_id3) : "black";
                    }
                })
        }
    
        scatter_svg.append("g")
            .attr("font-family", "sans-serif")
            .selectAll("text")
            .data(scatter_data)
            .join("text")
            .attr("class","marker_label")
            .attr("dy", "0.35em")
            .attr("x", function (d) { return sp_size!='none' && d.name in colors[sp_size] ? x(d.x) + 2 + (2+colors[sp_size][d.name][1] * 8) : x(d.x) + label_offset;} )
            .attr("y", d => y(d.y))
            .text(function (d) {
                if (label == 'set1') {
                    return d.set1_cons;
                } else if (label == 'set2') {
                    return d.set2_cons;
                } else if (label == 'combined') {
                    return d.comb_cons;
                } else {
                    return d.name;
                }
            })
            // .attr("font-size", function (d) { return sp_size != 'none' && label_only && d.name in colors[sp_size] ? 5 + colors[sp_size][d.name][1] * 10 : 8; })
            .attr("font-size", 8)
            .attr("font-weight", font_weight)
            .attr("text-anchor",text_anchor)
            .style("fill", function (d) {
                if (sp_color == "segment") {
                    return colors[sp_color][d.name];
                } else if (sp_color == "network") {
                    network_group = filtered_cluster_groups.filter(l => l.includes(d.name));
                    network_group_id = filtered_cluster_groups.indexOf(network_group[0]);
                    if (filtered_cluster_groups.length >  10 || 1==1) {
                        return d3v4.interpolateTurbo(network_group_id/filtered_cluster_groups.length);
                    } else {
                        return d3v4.schemePaired[network_group_id]
                    }
                } else {
                    return label_only && d.name in colors[sp_color] ? color_by_scale(colors[sp_color][d.name][1], color_id1, color_id2, color_id3) : "black";
                }
            });
    

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
        'set_presence' : 'Set specific presense'
    }
    var select_data_options = ''
    $.each(nice_index_names, function (key, description) {
        select_data_options += '<option value="' + key + '">' + description + '</option>';
    });
    var color_options = {
        'w': 'white',
        'r': 'red',
        'b': 'blue',
        'black': 'black',
        'g': 'green',
        'y': 'yellow',
        'o': 'orange',
        'p': 'purple',
    }
    var select_color_options = ''
    $.each(color_options, function (key, description) {
        checked = '';
        select_color_options += '<option value="' + description + '">' + description + '</option>';
    });

    function create_overlay() {
        var newDiv = document.createElement("div");

        $(containerSelector_hash).find(".controls-panel").remove();
        newDiv.setAttribute("class", "controls-panel");
        content = '<span class="pull-right scatter_controls_toggle" style="cursor: pointer;"><span class="glyphicon glyphicon-option-horizontal btn-download png"></span></span>' +
            '<div class="scatter_options" style="display: block; min-width: 120px;">' +
            '<div><strong>Label Format</strong></div><div></div>' +
            '<div>Text</div><div><select id="change_label" class="change_axis">' +
            '<option value="gn">Generic Number</option>' +
            '<option value="combined">Consensus AA (combined sets) + Generic Number</option>' +
            '<option value="set1">Consensus AA (set1) + Generic Number</option>' +
            '<option value="set2">Consensus AA (set2) + Generic Number</option>' +
            '</select></div> ' +
            '<div>Font Size</div><div><input id="change_text_size" style="width:80px;" type="range" min="0" max="50" step="any" value="12"></div>' +
            '<div>Replace marker with label</div><div><input id="label_only" class="change_axis" type="checkbox"></div>';
            // 'Only plot kept positions <input id="color_filtered" type="checkbox" checked><br>';
        
        index_names = { 0: 'core_distance', 1: 'a_angle', 2: 'outer_angle', 3: 'tau', 4: 'phi', 5: 'psi', 6: 'sasa', 7: 'rsa', 8: 'theta', 9: 'hse', 10: 'tau_angle' }
    
        content += '<div><strong>Visualisation</strong></div><div></div>' +
                   '<div>X-axis</div><div><select id="change_x" class="change_axis">' +
            select_data_options +
            '</select></div>'
            ;
        content += '<div>Y-axis</div><div><select id="change_y" class="change_axis">' +
            select_data_options +
            '</select></div>';
            
        content += '<div>Marker size</div><div><select id="sp_size" class="change_axis">' +
        '<option value="none">Fixed</option>' +
        select_data_options +
        '</select></div>';
        content += '<div>Label color</div>';
        content += '<div><select id="sp_color" class="change_axis">' +
            '<option value="segment">Segment</option>' +
            '<option value="network">Network group</option>' +
            select_data_options +
            '</select></div>' +
            '<div>Color scale</div><div>' +
                '<select id=color1 class=color>' +
                select_color_options +
                '</select>' +
                '<select id=color2 class=color>' +
                select_color_options +
                '</select>' +
                '<select id=color3 class=color>' +
                '<option value="none">None</option>' +
                select_color_options +
                '</select>'+
            '</div > ';
        content += '</div>';
        newDiv.innerHTML = content;

        $(containerSelector_hash).prepend(newDiv);
        $(containerSelector_hash).find(".options").toggle();

        $(containerSelector_hash).find(".scatter_controls_toggle").click(function() {
            $(containerSelector_hash).find(".scatter_options").slideToggle();
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

        d3v4.select(containerSelector_hash).select("#change_text_size")
            .on("input", change_text_size);
        
        function change_text_size() {
            console.log('change text size!',this.value)
            original_size = 8;
            scaled_size = original_size * this.value / 12;
            console.log('change text size', this.value,scaled_size);
            $(containerSelector_hash).find('.marker_label').each(function () {
                $(this).attr("font-size", scaled_size);
            });
        }

        $(containerSelector_hash + " .change_axis").on("change", function () {
            console.log('change axis!');
            draw_plot();
        });

        $(containerSelector_hash).find(".scatter_options").css("display", "grid");
        // $(containerSelector_hash).find(".scatter_options").css("grid-template-columns", "50px 50px 50px");
        // $(containerSelector_hash).find(".scatter_options").css("grid-template-rows", "auto");
        $(containerSelector_hash).find(".scatter_options").css("grid-template-areas", "'. .'");

        ;

    }


    create_overlay();

    // Init
    x_axis_type = $(containerSelector_hash + " #change_x").val("distance");
    y_axis_type = $(containerSelector_hash + " #change_y").val("outer_angle");
    // sp_color = $(containerSelector_hash + " #sp_color").val();
    // sp_size = $(containerSelector_hash + " #sp_size").val();
    // label_only = $(containerSelector_hash + " #label_only").prop("checked");
    // label = $(containerSelector_hash + " #change_label").val();


    color_id1  = $(containerSelector_hash+" #color1").val("blue");
    color_id2  = $(containerSelector_hash+" #color2").val("red");
    // color_id3  = $(containerSelector_hash+" #color3").val();


    draw_plot();
    
}

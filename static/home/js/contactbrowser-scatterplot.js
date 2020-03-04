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

    $.each(data['distances'], function (gn, dis) {
        seq_pos = data['snakeplot_lookup'][gn];
        // console.log(gn, dis['avg'], seq_pos);
        value = dis['avg'];
        scale = Math.abs(value) / data['ngl_max_diff_distance'];
        var color = { r: 255, g: 255, b: 255 };
        var color2 = { r: 255, g: 255, b: 255 };
        if (value < 0) {
            // if the header is a set two, then make it red
            color = { r: 255, g: 255-(255-153)*scale, b: 255-(255-153)*scale }; //red
            color2 = { r: 255-(255-128)*scale, g: 255-(255)*scale, b: 255-(255-128)*scale }; //purple
        } else if (value > 0) {
            // Positive numbers are blue either cos they are set 1 or cos "set 1 has most"
            // This is also used for single set/structure
            color = { r: 255-(255-153)*scale, g: 255-(255-204)*scale, b: 255 }; //blue
            color2 = { r: 255-(255)*scale, b: 255-(255)*scale, g: 128 }; //green
        }
        var hex = rgb2hex(color.r, color.g, color.b);
        var hex2 = rgb2hex(color2.r, color2.g, color2.b);
        // grey
        var color_grey = { r: 255*(1-scale), g: 255*(1-scale), b: 255*(1-scale) };
        var hex_grey = rgb2hex(color_grey.r, color_grey.g, color_grey.b);
        colors['distance'][seq_pos] = [hex,value,scale,hex_grey,hex2];
    });

    index_names = { 0: 'core_distance', 1: 'a_angle', 2: 'outer_angle', 3: 'tau', 4: 'phi', 5: 'psi', 6: 'sasa', 7: 'rsa', 8: 'theta', 9: 'hse', 10: 'tau_angle' }
    // get maximum values
    var max_values = {}
    $.each(data['tab4'], function (gn, v) {
        $.each(v["angles"], function (i, a) {
            if (!(i in max_values)) {
                max_values[i] = 0;
                colors[index_names[i]] = {};
            }
            if (Math.abs(a[0])>max_values[i]) max_values[i] = Math.abs(a[0])
            // console.log(gn, i, a);
        });
    });
    $.each(data['tab4'], function (gn, v) {
        $.each(v["angles"], function (i, a) {
            seq_pos = data['snakeplot_lookup'][gn];
            value = a[0];
            scale = Math.abs(value) / max_values[i];
            var color = { r: 255, g: 255, b: 255 };
            var color2 = { r: 255, g: 255, b: 255 };
            if (value < 0) {
                // if the header is a set two, then make it red
                color = { r: 255, g: 255-(255-153)*scale, b: 255-(255-153)*scale }; //red
                color2 = { r: 255-(255-128)*scale, g: 255-(255)*scale, b: 255-(255-128)*scale }; //purple
            } else if (value > 0) {
                // Positive numbers are blue either cos they are set 1 or cos "set 1 has most"
                // This is also used for single set/structure
                color = { r: 255-(255-153)*scale, g: 255-(255-204)*scale, b: 255 }; //blue
                color2 = { r: 255-(255-153)*scale, b: 255-(255-204)*scale, g: 255 }; //green
            }
            var hex = rgb2hex(color.r, color.g, color.b);
            var hex2 = rgb2hex(color2.r, color2.g, color2.b);
            // grey
            var color_grey = { r: 255*(1-scale), g: 255*(1-scale), b: 255*(1-scale) }; 
            var hex_grey = rgb2hex(color_grey.r, color_grey.g, color_grey.b);
            colors[index_names[i]][seq_pos] = [hex,value,scale,hex_grey,hex2];
            if (!(i in max_values)) max_values[i] = 0;
            // console.log(gn, i, a);
        });
    });
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
            console.log('filtered', gn,x,index_names_rev[x]);
            console.log('filtered', gn,y,index_names_rev[y]);
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

    x_axis_type = 'distance';
    y_axis_type = 'outer_angle';
    draw_plot(x_axis_type, y_axis_type);
    function draw_plot(x_axis_type,y_axis_type) {
        plot_data = get_values(x_axis_type,y_axis_type)
        console.log('plot_data', plot_data);

        var trace1 = {
            x: plot_data[0],
            y: plot_data[1],
            mode: 'markers+text',
            type: 'scatter',
            name: 'Team A',
            text: plot_data[2],
            textposition: 'top center',
            textfont: {
            family:  'Raleway, sans-serif'
            },
            marker: { size: 12 }
        };
        
        var data = [ trace1 ];
        
        var layout = {
            xaxis: {
                title: {
                    text: nice_index_names[x_axis_type],
                    font: {
                        family: 'Courier New, monospace',
                        size: 18,
                        color: '#7f7f7f'
                    }
                }
            },
            yaxis: {
                title: {
                    text: nice_index_names[y_axis_type],
                    font: {
                        family: 'Courier New, monospace',
                        size: 18,
                        color: '#7f7f7f'
                    }
                }
            },
            legend: {
            y: 0.5,
            yref: 'paper',
            font: {
                family: 'Arial, sans-serif',
                size: 20,
                color: 'grey',
            }
            },
            title:'Scatterplot of residue properities'
        };
        
        Plotly.newPlot(containerSelector, data, layout);
    }

    create_overlay();

    function create_overlay() {
        var newDiv = document.createElement("div");

        $(containerSelector_hash).find(".controls-panel").remove();

        newDiv.setAttribute("class", "controls-panel");
        content = '<span class="pull-right snakeplot_controls_toggle" style="cursor: pointer;"><span class="glyphicon glyphicon-option-horizontal btn-download png"></span></span><span class="options" style="display: block; min-width: 120px;">' +
            // 'Generic number with AA<input id="generic" type="checkbox"><br>' +
            'Only plot kept positions <input id="color_filtered" type="checkbox" checked><br>';
        
        index_names = { 0: 'core_distance', 1: 'a_angle', 2: 'outer_angle', 3: 'tau', 4: 'phi', 5: 'psi', 6: 'sasa', 7: 'rsa', 8: 'theta', 9: 'hse', 10: 'tau_angle' }
    
        content += '<table><tr><th>X-axis</th><th>Y-axis</th></tr><tr><td><select id="change_x" class="change_axis">' +
            '<option value="distance">Distance to all</option>' +
            '<option value="core_distance">Distance to 7TM axis</option>' +
            '<option value="a_angle">Angle to helix</option>' + 
            '<option value="outer_angle">Rotamer</option>' + 
            '<option value="tau_angle">Tau angle</option>' + 
            '<option value="theta">Theta angle</option>' + 
            '<option value="phi">Phi dihedral</option>' + 
            '<option value="psi">Psi dihedral</option>' + 
            '<option value="tau">Theta dihedral</option>' + 
            '<option value="hse">HSE</option>' + 
            '<option value="sasa">SASA</option>' + 
            '<option value="rsa">RSA</option>' + 
            '<option value="presense">Position Presence</option>' +
            '</select></td>'
            ;
        content += '<td><select id="change_y" class="change_axis">' +
            '<option value="distance">Distance to all</option>' +
            '<option value="core_distance">Distance to 7TM axis</option>' +
            '<option value="a_angle">Angle to helix</option>' + 
            '<option value="outer_angle" selected>Rotamer</option>' + 
            '<option value="tau_angle">Tau angle</option>' + 
            '<option value="theta">Theta angle</option>' + 
            '<option value="phi">Phi dihedral</option>' + 
            '<option value="psi">Psi dihedral</option>' + 
            '<option value="tau">Theta dihedral</option>' + 
            '<option value="hse">HSE</option>' + 
            '<option value="sasa">SASA</option>' + 
            '<option value="rsa">RSA</option>' + 
            '<option value="presense">Position Presence</option>' +
            '</select></td></tr></table>'
                ;
        content += '</span>';
        newDiv.innerHTML = content;

        $(containerSelector_hash).prepend(newDiv);
        $(containerSelector_hash).find(".options").toggle();

        $(containerSelector_hash).find(".snakeplot_controls_toggle").click(function() {
            $(containerSelector_hash).find(".options").slideToggle();
        });

        d3.select(containerSelector_hash).select("#generic").on("change", function () {
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
            x_axis_type = $(containerSelector_hash + " #change_x").val();
            y_axis_type = $(containerSelector_hash + " #change_y").val();
            console.log('change axis!',x_axis_type,y_axis_type);
            draw_plot(x_axis_type, y_axis_type);
        });

        

    }
    
}

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

    x_axis_type = 'distance';
    y_axis_type = 'outer_angle';

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
                  text: x_axis_type,
                  font: {
                      family: 'Courier New, monospace',
                      size: 18,
                      color: '#7f7f7f'
                  }
              }
        },
        yaxis: {
            title: {
                text: y_axis_type,
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


    create_overlay();

    function create_overlay() {
        var newDiv = document.createElement("div");

        $(containerSelector_hash).find(".controls-panel").remove();

        newDiv.setAttribute("class", "controls-panel");
        content = '<span class="pull-right snakeplot_controls_toggle" style="cursor: pointer;"><span class="glyphicon glyphicon-option-horizontal btn-download png"></span></span><span class="options" style="display: block; min-width: 120px;">' +
            'Generic number with AA<input id="generic" type="checkbox"><br>' +
            'Only color filtered set <input id="color_filtered" type="checkbox" checked><br>';
        
        index_names = { 0: 'core_distance', 1: 'a_angle', 2: 'outer_angle', 3: 'tau', 4: 'phi', 5: 'psi', 6: 'sasa', 7: 'rsa', 8: 'theta', 9: 'hse', 10: 'tau_angle' }
    
        content += '<table><tr><td>Fill color:</td><td><select id="snakeplot_color">' +
            '<option value="none">None</option>' +
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
            '</select></td><td>' +
            '<select id=fill_color>' +
            '<option value=0>Red-Blue</option>'+
            '<option value=4>Purple-Green</option>'+
            '<option value=3>Grey (Abs values)</option>'+
            '</select></td></tr>'
            ;
        content += '<tr><td>Border color:</td><td><select id="snakeplot_color_border">' +
            '<option value="none">None</option>' +
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
            '</select></td><td>' +
            '<select id=border_color>' +
            '<option value=0>Red-Blue</option>'+
            '<option value=4 >Purple-Green</option>'+
            '<option value=3 selected>Grey (Abs values)</option>'+
            '</select></td></tr>'
            ;
        content += '<tr><td>Text color:</td><td><select id="snakeplot_color_text">' +
                '<option value="none">None</option>' +
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
                '</select></td><td>' +
                '<select id=text_color>' +
                '<option value=0>Red-Blue</option>'+
                '<option value=4 selected>Purple-Green</option>'+
                '<option value=3>Grey (Abs values)</option>'+
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


        d3.select(containerSelector_hash).select("#color_filtered").on("change", function () {
            change_fill();
            change_stroke();
            change_text();
        });



        d3.select(containerSelector_hash).select("#snakeplot_color").on("change", function () {
            change_fill();
        });

        d3.select(containerSelector_hash).select("#fill_color").on("change", function () {
            change_fill();
        });

        d3.select(containerSelector_hash).select("#border_color").on("change", function () {
            change_stroke();
        });

        d3.select(containerSelector_hash).select("#snakeplot_color_border").on("change", function () {
            change_stroke();
        });

        d3.select(containerSelector_hash).select("#text_color").on("change", function () {
            change_text();
        });

        d3.select(containerSelector_hash).select("#snakeplot_color_text").on("change", function () {
            change_text();
        });

        function change_fill() {
            fill_color = $(containerSelector_hash + " #snakeplot_color").val();
            console.log('change fill color to!', fill_color);
            
            color_filtered = d3.select(containerSelector_hash).select("#color_filtered").property("checked");
            color_id  = $(containerSelector_hash+" #fill_color").val();
            fill_grey = color_id == '3' ? true : false;

            $(containerSelector_hash).find('.rcircle').each(function () {

                pos_id = $(this).attr('id');
                original_title = $(this).attr('original_title');
                gn = false;
                if (original_title.split(" ")[1].length) {
                    // this has GN
                    seg = original_title.split(" ")[1].split(".")[0];
                    gn = original_title.split(" ")[1].split("x")[1];
                }
                if (color_filtered && (!gn || !filtered_gns.includes(seg+"x"+gn))) {
                    $(this).attr("fill", "#fff");
                    $(this).attr("stroke", "#ccc");  
                    $(this).attr("stroke-width", 1); 
                    $(containerSelector_hash).find('#' + pos_id + 't').attr("fill", "#ddd");
                    return true;
                }

                if (fill_color!='none') {
                    color = pos_id in colors[fill_color] ? colors[fill_color][pos_id] : ["#fff",0,0,"#fff","#fff"];
                } else {
                    color = ["#fff", 0, 0,"#fff","#fff"];     
                }
                
                $(this).attr("fill", color[color_id]);
                $(this).attr("fill_value", color[2]);  

                if ($(containerSelector_hash + " #snakeplot_color_text").val() == 'none') {
                    $(containerSelector_hash).find('#'+pos_id+'t').removeAttr("fill");
                    $(containerSelector_hash).find('#' + pos_id + 't').attr("fill", "#000");
                    if (fill_grey && color[2]>0.5) {
                        $(containerSelector_hash).find('#'+pos_id+'t').attr("fill", "#fff");
                    }
                }

                $(this).css("opacity", 1);
                $(this).removeAttr("fill-opacity");
                $(containerSelector_hash).find('#'+pos_id+'t').css("opacity", 1);
                // if (color[2] < 0.05 && $(this).attr("stroke_value")<0.05) {
                //     $(this).css("opacity", 0.3);
                //     $(containerSelector_hash).find('#'+pos_id+'t').css("opacity", 0.1);
                // }
                             
            });
        }

        function change_stroke() {
            fill_color = $(containerSelector_hash + " #snakeplot_color_border").val();
            color_id  = $(containerSelector_hash+" #border_color").val();
            fill_grey = color_id == '3' ? true : false;
            console.log('change stroke color to!');

            $(containerSelector_hash).find('.rcircle').each(function () {
                original_title = $(this).attr('original_title');
                gn = false;
                pos_id = $(this).attr('id');
                if (original_title.split(" ")[1].length) {
                    // this has GN
                    seg = original_title.split(" ")[1].split(".")[0];
                    gn = original_title.split(" ")[1].split("x")[1];
                }
                if (color_filtered && (!gn || !filtered_gns.includes(seg+"x"+gn))) {
                    $(this).attr("fill", "#fff");
                    $(this).attr("stroke", "#ccc");  
                    $(this).attr("stroke-width", 1); 
                    $(containerSelector_hash).find('#' + pos_id + 't').attr("fill", "#ddd");
                    return true;
                }
                if (fill_color!='none') {
                    color = pos_id in colors[fill_color] ? colors[fill_color][pos_id] : ["#fff",0,0,"#000","#000"];
                } else {
                    color = ["#000", 0, 0,"#000","#000"];     
                }
                
                
                $(this).attr("stroke", color[color_id]);
                $(this).attr("stroke_value", color[2]);          
                $(this).attr("stroke-width", 1); 
                // $(this).css("opacity", 0.3);
                if (color[2] > 0.1) {    
                    $(this).attr("stroke-width", 3); 
                    $(this).css("opacity", 1);
                } else {
                    $(this).attr("stroke", "#ccc");
                }
            });
        }


        function change_text() {
            fill_color = $(containerSelector_hash + " #snakeplot_color_text").val();
            color_id  = $(containerSelector_hash+" #text_color").val();
            fill_grey = color_id == '3' ? true : false;
            console.log('change text color to!');

            $(containerSelector_hash).find('.rcircle').each(function () {
                original_title = $(this).attr('original_title');
                gn = false;
                pos_id = $(this).attr('id');
                if (original_title.split(" ")[1].length) {
                    // this has GN
                    seg = original_title.split(" ")[1].split(".")[0];
                    gn = original_title.split(" ")[1].split("x")[1];
                }
                if (color_filtered && (!gn || !filtered_gns.includes(seg+"x"+gn))) {
                    $(containerSelector_hash).find('#' + pos_id + 't').attr("fill", "#ddd");
                    return true;
                }
                if (fill_color!='none') {
                    color = pos_id in colors[fill_color] ? colors[fill_color][pos_id] : ["#000",0,0,"#000","#000"];
                } else {
                    color = ["#000", 0, 0,"#000","#000"];     
                }
                
                $(containerSelector_hash).find('#' + pos_id + 't').attr("fill", color[color_id]);
                // $(this).css("opacity", 0.3);
                // if (color[2] > 0.1) {    
                    
                // } else {
                //     $(containerSelector_hash).find('#' + pos_id + 't').attr("fill", "#ddd");
                // }
            });
        }

    }
    
}

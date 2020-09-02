
var color_scale_colors = {}

color_scale_colors.red = {red:255, green:0, blue: 0};
color_scale_colors.red = {red:195, green:74, blue: 54};
color_scale_colors.blue = {red:0, green:0, blue: 255};
color_scale_colors.blue = {red:0, green:140, blue: 204};
color_scale_colors.green = {red:0, green:255, blue: 0};
color_scale_colors.green = {red:0, green:201, blue: 167};
color_scale_colors.white = {red:255, green:255, blue: 255};
color_scale_colors.yellow = {red:255, green:255, blue: 0};
color_scale_colors.yellow = {red:255, green:255, blue: 0};
color_scale_colors.black = { red: 0, green: 0, blue: 0 };
color_scale_colors.grey = { red: 211, green: 211, blue: 211 };
color_scale_colors.orange = { red: 255, green: 150, blue: 113 };
color_scale_colors.purple = { red: 128, green: 0, blue: 128 };
color_scale_colors.brown = { red: 165, green: 42, blue: 42 };
color_scale_colors.olive = { red: 128, green: 128, blue: 0 };
color_scale_colors.magenta = { red: 255, green: 0, blue: 255 };
color_scale_colors.pink = { red: 255, green: 20, blue: 147 };
color_scale_colors.false = false;

function color_by_category(value, possibilities) {
    rgb = d3.rgb(d3.interpolateSpectral(value / possibilities));
    return rgb.hex();
}

function prepare_residue_colors(data) {
    var colors_gn = {}
    // var colors = {}
    colors_gn['distance'] = {}
    colors_gn['distance_abs'] = {}
    colors_gn['network'] = {}
    colors_gn['set_presence'] = {}
    colors_gn['ligand'] = {}
    colors_gn['complex'] = {}
    colors_gn['ligandcomplex'] = {}
    colors_gn['mutations'] = {}
    colors_gn['conservation'] = {}

    index_names = { 0: 'core_distance', 1: 'a_angle', 2: 'outer_angle', 3: 'tau', 4: 'phi', 5: 'psi', 6: 'sasa', 7: 'rsa', 8: 'theta', 9: 'hse', 10: 'tau_angle', 11:'rotation_angle' }
    neg_and_positives = ['core_distance','sasa','rsa', 'hse']
    nice_index_names = {
        'a_angle' : 'Angle to helix&7TM axes',
        'outer_angle': 'Rotamer',
        'rotation_angle': 'Rotation angle',
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
        'network' : 'Network group no.',
        'set_presence' : 'Set specific presense',
        'ligand' : 'Ligand interactions freq',
        'complex' : 'G protein interactions',
        'ligandcomplex' : 'Ligand and G protein interactions',
        'mutations' : 'Mutations with >5 fold effect',
        'conservation' : 'Conservation of set(s) consensus AA in class'
    }

    $.each(data['snakeplot_lookup_aa_cons'], function (gn, cons) {
        seq_pos = data['snakeplot_lookup'][gn];
        scale = cons / 100
        colors_gn['conservation'][gn] = [cons,scale,100];
    } )


    $.each(filtered_cluster_groups, function (id, group) {
        $.each(group, function (i, gn) {
            scale = id / filtered_cluster_groups.length
            colors_gn['network'][gn] = [id,scale,filtered_cluster_groups.length];
        })
    })

    $.each(filtered_gns_presence, function (gn, value) {
        scale = value / 1
        colors_gn['set_presence'][gn] = [value,scale,1];
    })

    max_ligand_interactions = Math.max(...Object.values(data['class_ligand_interactions'])) 
    $.each(data['class_ligand_interactions'], function (gn, value) {
        scale = value / max_ligand_interactions;
        colors_gn['ligand'][gn] = [value,scale,max_ligand_interactions];
    })


    max_complex_interactions = Math.max(...Object.values(data['class_complex_interactions'])) 
    $.each(data['class_complex_interactions'], function (gn, value) {
        scale = value / max_complex_interactions;
        colors_gn['complex'][gn] = [value,scale,max_complex_interactions];
    })

    $.each(data['snakeplot_lookup'], function (gn, seq_pos) {
        complex_value = gn in data['class_complex_interactions'] ? data['class_complex_interactions'][gn] : 0;
        ligand_value = gn in data['class_ligand_interactions'] ? data['class_ligand_interactions'][gn] : 0;

        if (complex_value == 0 && ligand_value == 0) {
            // either do nothing or put in 0 value
            // colors_gn['ligandcomplex'][seq_pos] = [value,scale,max_mutations];
        } else {
            if (complex_value/max_complex_interactions >= ligand_value/max_ligand_interactions) {
                // treat complex values 0.5-1 range
                value = complex_value
                scale = 0.5+(0.5)*(value / max_complex_interactions);
                colors_gn['ligandcomplex'][gn] = [value,scale,max_complex_interactions];
            } else {
                // treat ligand values 0-0.5 range
                value = ligand_value
                scale = 0.5-(0.5)*(value / max_ligand_interactions);
                colors_gn['ligandcomplex'][gn] = [value,scale,max_ligand_interactions];
            }
        }
    })

    max_mutations = Math.max(...Object.values(data['class_mutations'])) 
    $.each(data['class_mutations'], function (gn, value) {
        scale = value / max_mutations;
        colors_gn['mutations'][gn] = [value,scale,max_mutations];
    })

    $.each(data['distances'], function (gn, dis) {
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
        colors_gn['distance'][gn] = [value,scale_abs,data['ngl_max_diff_distance']];
        colors_gn['distance_abs'][gn] = [Math.abs(value), scale, data['ngl_max_diff_distance']];

    });

    // get maximum values
    var max_values = {}
    $.each(data['tab4'], function (gn, v) {
        $.each(v["angles"], function (i, a) {
            value = $.isNumeric(a) ? a : a[0]; // deal with single structure and two sets
            if (!(i in max_values)) {
                max_values[i] = 0;
                colors_gn[index_names[i]] = {};
                if (neg_and_positives.includes(index_names[i])) {
                    colors_gn[index_names[i]+"_abs"] = {};
                }
            }
            if (Math.abs(value)>max_values[i]) max_values[i] = Math.abs(value)
        });
    });
    $.each(data['tab4'], function (gn, v) {
        $.each(v["angles"], function (i, a) {
            value = $.isNumeric(a) ? a : a[0]; // deal with single structure and two sets
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
                colors_gn[index_names[i]][gn] = [value,scale_abs,max_values[i]];
                colors_gn[index_names[i] + '_abs'][gn] = [Math.abs(value), scale, max_values[i]];
            } else {
                colors_gn[index_names[i]][gn] = [value, scale, max_values[i]];
            }
            if (!(i in max_values)) max_values[i] = 0;
        });
    });

    console.log('colors_gn', colors_gn);
    return colors_gn;
}

function color_by_scale(scale, color1, color2, color3) {

    // case "rwb": // red-white-blue
    // case "bwr": // blue-white-red
    // case "ryg": // red-yellow-green
    // case "gyr": // green-yellow-red
    // case "rgb":
    // case "wr": // white-red
    // case "wg": // white-green
    // case "wb": // white-blue
    // case "rb": // red-blue
    // case "br": // blue-red
    // case "wp": // white-purple
    // case "grey": // grey

    // return numberToColorGradient(scale, 1, color)

    // var colors = {}

    // colors.red = {red:255, green:0, blue: 0};
    // colors.red = {red:195, green:74, blue: 54};
    // colors.blue = {red:0, green:0, blue: 255};
    // colors.blue = {red:0, green:140, blue: 204};
    // colors.green = {red:0, green:255, blue: 0};
    // colors.green = {red:0, green:201, blue: 167};
    // colors.white = {red:255, green:255, blue: 255};
    // colors.yellow = {red:255, green:255, blue: 0};
    // colors.yellow = {red:255, green:255, blue: 0};
    // colors.black = { red: 0, green: 0, blue: 0 };
    // colors.orange = { red: 255, green: 150, blue: 113 };
    // colors.purple = { red: 128, green: 0, blue: 128 };
    // colors.false = false;
    
    return colorGradient(scale, color_scale_colors[color1], color_scale_colors[color2], color_scale_colors[color3])
}

function createSnakeplot(data, containerSelector) {


    $(containerSelector).html('')
    $(containerSelector).addClass("snakeplot");
    $(containerSelector).css("position","relative");

    $(containerSelector).html(data['snakeplot']);

    max_x = $(containerSelector).find("svg").attr("width");
    max_y = $(containerSelector).find("svg").attr("height");
    
    var svgContainer = $(containerSelector).find("svg").closest("svg")
        .attr("viewBox", 0 + " " + 0 + " " + (max_x) + " " + (max_y))
        .attr("width", "100%")
        .attr("style", "height: 500px");
    
    $(containerSelector).find(".long").hide();
    maxmin();
    create_legend();


    var colors = {}
    colors['distance'] = {}
    colors['distance_abs'] = {}
    colors['network'] = {}
    colors['set_presence'] = {}
    colors['ligand'] = {}
    colors['complex'] = {}
    colors['ligandcomplex'] = {}
    colors['mutations'] = {}
    colors['conservation'] = {}
    colors['gn_aa'] = {}

    index_names = { 0: 'core_distance', 1: 'a_angle', 2: 'outer_angle', 3: 'tau', 4: 'phi', 5: 'psi', 6: 'sasa', 7: 'rsa', 8: 'theta', 9: 'hse', 10: 'tau_angle', 11:'rotation_angle' }
    neg_and_positives = ['core_distance','sasa','rsa', 'hse']
    nice_index_names = {
        'a_angle' : 'Angle to helix&7TM axes',
        'outer_angle': 'Rotamer',
        'rotation_angle': 'Rotation angle',
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
        'network' : 'Network group no.',
        'set_presence' : 'Set specific presense',
        'ligand' : 'Ligand interactions freq',
        'complex' : 'G protein interactions',
        'ligandcomplex' : 'Ligand and G protein interactions',
        'mutations' : 'Mutations with >5 fold effect',
        'conservation' : 'Conservation of set(s) consensus AA in class'
    }
    path_groups = {}
    path_groups_lookup = {}
    $(containerSelector).find('.helix_path').each(function () {
        path_id = $(this).attr('id');
        res_ids = $(this).attr('previous_res');
        path_groups[path_id] = { 'distance': [], 'distance_abs': [] };
        
        $.each(index_names, function (i, name) {
            path_groups[path_id][name] = [];
            if (neg_and_positives.includes(name)) {
                path_groups[path_id][name+"_abs"] = [];
            }
        })

        $.each(res_ids.split(","), function (i, res_id) {
            path_groups_lookup[res_id] = path_id;
        });
    });

    $.each(data['snakeplot_lookup_aa_cons'], function (gn, cons) {
        seq_pos = data['snakeplot_lookup'][gn];
        scale = cons / 100
        colors['conservation'][seq_pos] = [cons,scale,100];
    } )


    $.each(filtered_cluster_groups, function (id, group) {
        $.each(group, function (i, gn) {
            seq_pos = data['snakeplot_lookup'][gn];
            scale = id / filtered_cluster_groups.length
            colors['network'][seq_pos] = [id,scale,filtered_cluster_groups.length];
        })
    })

    $.each(filtered_gns_presence, function (gn, value) {
        seq_pos = data['snakeplot_lookup'][gn];
        scale = value / 1
        colors['set_presence'][seq_pos] = [value,scale,1];
    })

    max_ligand_interactions = Math.max(...Object.values(data['class_ligand_interactions'])) 
    $.each(data['class_ligand_interactions'], function (gn, value) {
        seq_pos = data['snakeplot_lookup'][gn];
        scale = value / max_ligand_interactions;
        colors['ligand'][seq_pos] = [value,scale,max_ligand_interactions];
    })


    max_complex_interactions = Math.max(...Object.values(data['class_complex_interactions'])) 
    $.each(data['class_complex_interactions'], function (gn, value) {
        seq_pos = data['snakeplot_lookup'][gn];
        scale = value / max_complex_interactions;
        colors['complex'][seq_pos] = [value,scale,max_complex_interactions];
    })

    $.each(data['snakeplot_lookup'], function (gn, seq_pos) {
        complex_value = gn in data['class_complex_interactions'] ? data['class_complex_interactions'][gn] : 0;
        ligand_value = gn in data['class_ligand_interactions'] ? data['class_ligand_interactions'][gn] : 0;

        if (complex_value == 0 && ligand_value == 0) {
            // either do nothing or put in 0 value
            // colors['ligandcomplex'][seq_pos] = [value,scale,max_mutations];
        } else {
            if (complex_value/max_complex_interactions >= ligand_value/max_ligand_interactions) {
                // treat complex values 0.5-1 range
                value = complex_value
                scale = 0.5+(0.5)*(value / max_complex_interactions);
                colors['ligandcomplex'][seq_pos] = [value,scale,max_complex_interactions];
            } else {
                // treat ligand values 0-0.5 range
                value = ligand_value
                scale = 0.5-(0.5)*(value / max_ligand_interactions);
                colors['ligandcomplex'][seq_pos] = [value,scale,max_ligand_interactions];
            }
            // console.log(gn,complex_value/max_complex_interactions,ligand_value/max_ligand_interactions, colors['ligandcomplex'][seq_pos])
        }
    })

    max_mutations = Math.max(...Object.values(data['class_mutations'])) 
    $.each(data['snakeplot_lookup'], function (gn, seq_pos) {
        if (gn in data['class_mutations']) {
            value = data['class_mutations'][gn];
        } else {
            value = 0;
        }
        scale = value / max_mutations;
        colors['mutations'][seq_pos] = [value,scale,max_mutations];
    })

    $.each(data['distances'], function (gn, dis) {
        seq_pos = data['snakeplot_lookup'][gn];
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
        
        if (seq_pos in path_groups_lookup) {
            path_groups[path_groups_lookup[seq_pos]]['distance'].push(scale_abs);
            path_groups[path_groups_lookup[seq_pos]]['distance_abs'].push(scale);
        }
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
        colors["gn_aa"][gn] = v["all_seq_cons"][0]
        $.each(v["angles"], function (i, a) {
            seq_pos = data['snakeplot_lookup'][gn];
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
                
                if (seq_pos in path_groups_lookup) {
                    path_groups[path_groups_lookup[seq_pos]][index_names[i]].push(scale_abs);
                    path_groups[path_groups_lookup[seq_pos]][index_names[i] + '_abs'].push(scale);
                }
            } else {
                colors[index_names[i]][seq_pos] = [value, scale, max_values[i]];
                if (seq_pos in path_groups_lookup) {
                    path_groups[path_groups_lookup[seq_pos]][index_names[i]].push(scale);
                }
            }
            if (!(i in max_values)) max_values[i] = 0;
        });
    });

    $(containerSelector).find('.rtext').each(function () {
        original_title = $(this).attr('original_title');
        if (original_title.split(" ")[1].length) {
            // this has GN
            pos = original_title.split(" ")[0].substring(1);
            gn = original_title.split(" ")[1];
            AA = colors["gn_aa"][gn];
        } else {
            gn = "-";
            seg = "";
            AA = original_title[0];
        }
        c_id = $(this).attr('id').slice(0, -1);
        $(this).attr('original_title', AA+pos+" "+gn)
        $(this).attr('title', AA+pos+" "+gn)
        $(this).text(AA);
        $(containerSelector).find('#'+c_id).attr('original_title', AA+pos+" "+gn)
        $(containerSelector).find('#'+c_id).attr('title', AA+pos+" "+gn)
    });

    $(containerSelector).find("text").tooltip({
        'container': 'body',
        'placement': 'top',
        'animation': false,
        'html' : true
    });

    $(containerSelector).find("circle").tooltip({
        'container': 'body',
        'placement': 'top',
        'animation': false,
        'html' : true
    });

    $(containerSelector).find("circle").hover(function(){
        $('.tooltip').css('top',parseInt($('.tooltip').css('top')) + 2.8 + 'px')
    });
    console.log('colors', colors);


    var color_options = {
        'rwb': 'red-white-blue',
        'bwr': 'blue-white-red',
        'ryg': 'red-yellow-green',
        'gyr': 'green-yellow-red',
        'wr': 'white-red',
        'wg': 'white-green',
        'wb': 'white-blue',
        'wy': 'white-yellow',
        'wo': 'white-orange',
        'wp': 'white-purple',
        'rb': 'red-blue',
        'br': 'blue-red',
        'grey': 'grey'
    }


    var color_options = {
        'black': 'black',
        'b': 'blue',
        'brown': 'brown',
        'g': 'green',
        'grey' : 'grey',
        'magenta': 'magenta',
        'olive': 'olive',
        'o': 'orange',
        'pink' : 'pink',
        'p': 'purple',
        'r': 'red',
        'w': 'white',
        'y': 'yellow',
    }

    var select_color_options_white = ''
    $.each(color_options, function (key, description) {
        checked = '';
        if (description == 'white') checked = 'selected';
        select_color_options_white += '<option value="' + description + '" '+checked+'>' + description + '</option>';
    });

    var select_color_options_red = ''
    $.each(color_options, function (key, description) {
        checked = '';
        if (description == 'red') checked = 'selected';
        select_color_options_red += '<option value="' + description + '" '+checked+'>' + description + '</option>';
    });

    var select_color_options_black = ''
    $.each(color_options, function (key, description) {
        checked = '';
        if (description == 'black') checked = 'selected';
        select_color_options_black += '<option value="' + description + '" '+checked+'>' + description + '</option>';
    });

    var select_color_options = ''
    $.each(color_options, function (key, description) {
        checked = '';
        select_color_options += '<option value="' + description + '">' + description + '</option>';
    });

    var select_data_options = ''
    $.each(nice_index_names, function (key, description) {
        select_data_options += '<option value="' + key + '">' + description + '</option>';
    });

    var select_data_options_backbone = ''
    $.each(nice_index_names, function (key, description) {
        if (description.includes('Distance') || key == 'a_angle')
            select_data_options_backbone += '<option value="' + key + '">' + description + '</option>';
    });

    var select_data_options_backbone_shift = ''
    $.each(nice_index_names, function (key, description) {
        if (description.includes('Distance'))
        select_data_options_backbone_shift += '<option value="' + key + '">' + description + '</option>';
    });
    // console.log(colors);

    text_content = {'aa': 'Amino Acid', 'gn': 'Generic Residue Number', 'both':'Both', 'none':'Empty'}
    var text_content_options = ''
    $.each(text_content, function (key, description) {
        text_content_options += '<option value="' + key + '">' + description + '</option>';
    });
    create_overlay();

    function create_overlay() {
        var newDiv = document.createElement("div");

        $(containerSelector).find(".controls-panel").remove();

        newDiv.setAttribute("class", "controls-panel");
        content = '<span class="pull-right snakeplot_controls_toggle" style="cursor: pointer;"><span class="glyphicon glyphicon-option-horizontal btn-download png"></span></span><span class="options" style="display: block; min-width: 120px;">';
        
        index_names = { 0: 'core_distance', 1: 'a_angle', 2: 'outer_angle', 3: 'tau', 4: 'phi', 5: 'psi', 6: 'sasa', 7: 'rsa', 8: 'theta', 9: 'hse', 10: 'tau_angle' }
    
        content += '<table><tr><td>Hide loops</td><td><input id=hide_loops type=checkbox></td></tr><tr><td>Residue text (kept positions):</td><td><select id="text_included" class="snakeplot_property_select">' +
        text_content_options +
        '</select>' + 
        '<div class="btn-group btn-toggle" id="included_color">' +
        '    <button class="btn btn-xs btn-primary active" value="black">Black</button>' +
        '    <button class="btn btn-xs btn-default" value="grey">Grey</button>' +
            '</div>' +
        '</td ></tr > <tr><td>Residue text (filtered out positions):</td><td><select id="text_excluded" class="snakeplot_property_select">' +
        text_content_options +
        '</select>' + 
        '<div class="btn-group btn-toggle" id="excluded_color">' +
        '    <button class="btn btn-xs btn-primary active" value="black">Black</button>' +
        '    <button class="btn btn-xs btn-default" value="grey">Grey</button>' +
            '</div>' +
        '</td ></tr ></table > '
        ;

        content += '<table><tr><th>Area</th><th>Property</th><th>Color1</th><th>Color2</th><th>Color3</th><th>Positions</th></tr>' +
            '<tr><td>Residue fill:</td><td><select id="snakeplot_color" class="residue_fill snakeplot_property_select">' +
            '<option value="none">None</option>' +
            select_data_options +
            '</select></td>' +
            '<td><input type="text" id="fill_color1" class="togglePaletteOnly_red residue_fill snakeplot_color_select" value="red" /></td>' +
            '<td><input type="text" id="fill_color2" class="togglePaletteOnly_blue residue_fill snakeplot_color_select" value="blue" /></td>' +
            '<td><input type="text" id="fill_color3" class="togglePaletteOnly_empty residue_fill snakeplot_color_select" value="" /></td>' +
            // '<td>' +
            // '<select id=fill_color1 class="fill_color residue_fill snakeplot_color_select">' +
            // select_color_options_white +
            // '</select></td><td>' +
            // '<select id=fill_color2 class="fill_color residue_fill snakeplot_color_select">' +
            // select_color_options_red +
            // '</select></td><td>' +
            // '<select id=fill_color3 class="fill_color residue_fill snakeplot_color_select">' +
            // '<option value="none">None</option>' +
            // select_color_options +
            // '</select></td>' +
            '<td>' + 
            '<div class="btn-group btn-toggle residue_fill" id="fill_filtered">' +
            '    <button class="btn btn-xs btn-primary active" value="true">Kept</button>' +
            '    <button class="btn btn-xs btn-default" value="false">All</button>' +
            '</div>' +
            '</td></tr> '
            ;
        content += '<tr><td>Residue border:</td><td><select id="snakeplot_color_border" class="residue_border snakeplot_property_select">' +
            '<option value="none">None</option>' +
            select_data_options +
            '</select></td>' +

            '<td><input type="text" id="border_color1" class="togglePaletteOnly_red border_color residue_border snakeplot_color_select" value="red" /></td>' +
            '<td><input type="text" id="border_color2" class="togglePaletteOnly_blue border_color residue_border snakeplot_color_select" value="blue" /></td>' +
            '<td><input type="text" id="border_color3" class="togglePaletteOnly_empty border_color residue_border snakeplot_color_select" value="" /></td>' +
            // '<td>' +
            // '<select id=border_color1 class="border_color residue_border snakeplot_color_select">' +
            // select_color_options_white +
            // '</select></td><td>' +
            // '<select id=border_color2 class="border_color residue_border snakeplot_color_select">' +
            // select_color_options_red +
            // '</select></td><td>' +
            // '<select id=border_color3 class="border_color residue_border snakeplot_color_select">' +
            // '<option value="none">None</option>' +
            // select_color_options +
            // '</select></td>' +
            '<td>' + 
            '<div class="btn-group btn-toggle residue_border" id="border_filtered">' +
            '    <button class="btn btn-xs btn-primary active" value="true">Kept</button>' +
            '    <button class="btn btn-xs btn-default" value="false">All</button>' +
            '</div>' +
            '</td></tr> '
            ;
        content += '<tr><td>Border thickness:</td><td><select id="snakeplot_border_stroke" class="residue_border">' +
            '<option>1</option>' +
            '<option>2</option>' +
            '<option SELECTED>3</option>' +
            '<option>4</option>' +
            '<option>5</option>' +
            '</td></tr>';
        content += '<tr><td>Residue text:</td><td><select id="snakeplot_color_text" class="residue_text snakeplot_property_select">' +
                '<option value="none">None</option>' +
                select_data_options +
                '</select></td>' +

                '<td><input type="text" id="text_color1" class="togglePaletteOnly_red text_color residue_border snakeplot_color_select" value="red" /></td>' +
                '<td><input type="text" id="text_color2" class="togglePaletteOnly_blue text_color residue_border snakeplot_color_select" value="blue" /></td>' +
                '<td><input type="text" id="text_color3" class="togglePaletteOnly_empty text_color residue_border snakeplot_color_select" value="" /></td>' +
                // '<td>' +
                // '<select id=text_color1 class="text_color residue_text snakeplot_color_select">' +
                // select_color_options_white +
                // '</select></td><td>' +
                // '<select id=text_color2 class="text_color residue_text snakeplot_color_select">' +
                // select_color_options_red +
                // '</select></td><td>' +
                // '<select id=text_color3 class="text_color residue_text snakeplot_color_select">' +
                // '<option value="none">None</option>' +
                // select_color_options +
                // '</select></td>' +
                '<td>' + 
                '<div class="btn-group btn-toggle residue_text" id="text_filtered">' +
                '    <button class="btn btn-xs btn-primary active" value="true">Kept</button>' +
                '    <button class="btn btn-xs btn-default" value="false">All</button>' +
                '</div>' +
                '</td></tr> '
                ;
        content += '<tr><td>Residue rotation:</td><td colspan=4><select id="snakeplot_color_rotation" class="residue_rotation snakeplot_property_select">' +
            '<option value="none">None</option>' +
            select_data_options +
            '</select></td>' +
            '<td>' + 
            '<div class="btn-group btn-toggle residue_rotation" id="rotation_filtered">' +
            '    <button class="btn btn-xs btn-primary active" value="true">Kept</button>' +
            '    <button class="btn btn-xs btn-default" value="false">All</button>' +
            '</div>' +
            '</td></tr> '
        content += '<tr><td colspan=6><hr></td></tr><tr><td>Backbone line:</td><td><select id="snakeplot_color_backbone" class="snakeplot_property_select">' +
                '<option value="none">None</option>' +
                select_data_options_backbone +
                '</select></td>' +

                '<td><input type="text" id="backbone_color1" class="togglePaletteOnly_red backbone_color snakeplot_color_select" value="red" /></td>' +
                '<td><input type="text" id="backbone_color2" class="togglePaletteOnly_blue backbone_color snakeplot_color_select" value="blue" /></td>' +
                '<td><input type="text" id="backbone_color3" class="togglePaletteOnly_empty backbone_color snakeplot_color_select" value="" /></td>' +
                // '<td>' +
                // '<select id=backbone_color1 class="backbone_color snakeplot_color_select">' +
                // select_color_options +
                // '</select></td><td>' +
                // '<select id=backbone_color2 class="backbone_color snakeplot_color_select">' +
                // select_color_options_black +
                // '</select></td><td>' +
                // '<select id=backbone_color3 class="backbone_color snakeplot_color_select">' +
                // '<option value="none">None</option>' +
                // select_color_options +
                // '</select></td>' +
                '</tr>'
                ;
        content += '<tr><td>Backbone shift:</td><td><select id="snakeplot_move_circle" class="snakeplot_property_select">' +
                '<option value="none">None</option>' +
                select_data_options_backbone_shift +
                '</select></td><td>' +
                '</td></tr></table>'
                ;
        content += '</span>';
        newDiv.innerHTML = content;

        $(containerSelector).prepend(newDiv);

        
        color_palette = [
            ["#000", "#444", "#666", "#999", "#ccc", "#eee", "#f3f3f3", "#fff"],
            ["#f00", "#f90", "#ff0", "#0f0", "#0ff", "#00f", "#90f", "#f0f"],
            ["#f4cccc", "#fce5cd", "#fff2cc", "#d9ead3", "#d0e0e3", "#cfe2f3", "#d9d2e9", "#ead1dc"],
            ["#ea9999", "#f9cb9c", "#ffe599", "#b6d7a8", "#a2c4c9", "#9fc5e8", "#b4a7d6", "#d5a6bd"],
            ["#e06666", "#f6b26b", "#ffd966", "#93c47d", "#76a5af", "#6fa8dc", "#8e7cc3", "#c27ba0"],
            ["#c00", "#e69138", "#f1c232", "#6aa84f", "#45818e", "#3d85c6", "#674ea7", "#a64d79"],
            ["#900", "#b45f06", "#bf9000", "#38761d", "#134f5c", "#0b5394", "#351c75", "#741b47"],
            ["#600", "#783f04", "#7f6000", "#274e13", "#0c343d", "#073763", "#20124d", "#4c1130"] //
        ];



        $(containerSelector+" .togglePaletteOnly_red").spectrum({
            showPaletteOnly: true,
            togglePaletteOnly: true,
            hideAfterPaletteSelect:true,
            togglePaletteMoreText: 'more',
            togglePaletteLessText: 'less',
            color: 'red',
            palette: color_palette
        });
        $(containerSelector+" .togglePaletteOnly_blue").spectrum({
            showPaletteOnly: true,
            togglePaletteOnly: true,
            hideAfterPaletteSelect:true,
            togglePaletteMoreText: 'more',
            togglePaletteLessText: 'less',
            color: 'blue',
            palette: color_palette
        });
        $(containerSelector+" .togglePaletteOnly_empty").spectrum({
            showPalette: true,
            togglePaletteOnly: true,
            hideAfterPaletteSelect:true,
            togglePaletteMoreText: 'more',
            togglePaletteLessText: 'less',
            allowEmpty:true,
            palette: color_palette
        });

        $(containerSelector).find(".options").toggle();

        $(containerSelector + ' .btn-toggle').click(function() {
            $(this).find('.btn').toggleClass('active');
            $(this).find('.btn').toggleClass('btn-primary');
            $(this).find('.btn').toggleClass('btn-default');
            id = $(this).attr("id");
            if (id=='fill_filtered') change_fill();
            if (id=='border_filtered') change_stroke();
            if (id == 'text_filtered') change_text();
            if (id == 'rotation_filtered') change_rotation();
            
            if (id == 'excluded_color') $(containerSelector + " #text_excluded")[0].dispatchEvent(new Event("change"));
            if (id == 'included_color') $(containerSelector + " #text_included")[0].dispatchEvent(new Event("change"));
            
            // $(this).find('.active').html($(this).find('.active').attr("value"));
            // $(this).find('.btn-default').html('&nbsp;');
            // toggle_best($(this).attr("mode"), $(this).attr("column"),$(this).find('.active').attr("value"))
        })


        $(containerSelector).find(".snakeplot_controls_toggle").click(function() {
            $(containerSelector).find(".options").slideToggle();
        });

        $(containerSelector).find("#hide_loops").click(function() {
            hide_loops = $(containerSelector + " #hide_loops").prop("checked");
            console.log("#hide_loops", hide_loops);
            loops = ['ICL1', 'ICL2', 'ICL3', 'ECL1', 'ECL2','ECL3'];
            $.each(loops, function (i, l) {
                if (hide_loops) {
                    $(containerSelector).find("." + l).hide();
                } else {
                    $(containerSelector).find("." + l+".short").show();
                }
            })
            maxmin();
        });

        

        d3.select(containerSelector).select("#text_included").on("change", function () {
            text_format = $(containerSelector + " #text_included").val();
            color = $(containerSelector + " #included_color").find(".active").attr('value');
            console.log('change to', text_format);
            $(containerSelector).find('.rtext').each(function () {
                original_title = $(this).attr('original_title');
                if (original_title.split(" ")[1].length) {
                    // this has GN
                    seg = original_title.split(" ")[1].split("x")[0];
                    gn = original_title.split(" ")[1].split("x")[1];
                    AA = original_title[0];
                    // AA = colors["gn_aa"][gn];
                } else {
                    gn = "-";
                    seg = "";
                    AA = original_title[0];
                }

                if ((!gn || !filtered_gns.includes(seg+"x"+gn))) {
                    return true;
                } else {
                }
                var fontsize = 16;
                var text = '';
                switch (text_format) {
                    case "gn":
                        text = gn;
                    break;
                    case "both":
                        fontsize = 12;

                        var isSafari = /^((?!chrome|android).)*safari/i.test(navigator.userAgent);
                        y_shift = isSafari ? 2 : 5;
                        
                        text = "<tspan x='"+$(this).attr('x')+"' y='"+(parseInt($(this).attr('y'))-y_shift)+"'>" + AA + "</tspan><tspan dy=9  x='"+$(this).attr('x')+"'>" + gn + "</tspan>";
                    break;
                    case "aa":
                            text = AA;
                    break;
                    case "none":
                            text = '';
                    break;
                }

                $(this).attr("font-size", fontsize);
                $(this).attr("fill", color);
                $(this).html(text);
            });
            create_legend();
        });

        d3.select(containerSelector).select("#text_excluded").on("change", function () {
            
            color = $(containerSelector + " #excluded_color").find(".active").attr('value');
            text_format = $(containerSelector + " #text_excluded").val();
            console.log('change to', text_format);
            $(containerSelector).find('.rtext').each(function () {
                original_title = $(this).attr('original_title');
                if (original_title.split(" ")[1].length) {
                    seg = original_title.split(" ")[1].split("x")[0];
                    // this has GN
                    gn = original_title.split(" ")[1].split("x")[1];
                    AA = original_title[0];
                } else {
                    gn = "-";
                    seg = "";
                    AA = original_title[0];
                }

                if ((!gn || !filtered_gns.includes(seg+"x"+gn))) {
                    
                } else {
                    return true;
                }
                var fontsize = 16;
                var text = '';
                switch (text_format) {
                    case "gn":
                        text = gn;
                    break;
                    case "both":
                        fontsize = 12;
                        text = "<tspan dy=-6>" + AA + "</tspan><tspan dy=9 dx=-10>" + gn + "</tspan>";
                    break;
                    case "aa":
                            text = AA;
                    break;
                    case "none":
                            text = '';
                    break;
                }

                $(this).attr("font-size", fontsize);
                $(this).attr("fill", color);
                $(this).html(text);
            });
        });

        d3.select(containerSelector).select("#color_filtered").on("change", function () {
            change_fill();
            change_stroke();
            change_text();
            create_legend();
        });


        $(containerSelector+ " .residue_fill").on("change", function () {
            change_fill();
            create_legend();
        });

        $(containerSelector+ " .residue_border").on("change", function () {
            change_stroke();
            create_legend();
        });


        $(containerSelector+ " .residue_text").on("change", function () {
            change_text();
            create_legend();
        });


        $(containerSelector+ " .backbone_color").on("change", function () {
            change_backbone();
            create_legend();
        });

        $(containerSelector+ " .residue_rotation").on("change", function () {
            change_rotation();
            create_legend();
        });

        d3.select(containerSelector).select("#snakeplot_color_backbone").on("change", function () {
            change_backbone();
            create_legend();
        });

        d3.select(containerSelector).select("#snakeplot_move_circle").on("change", function () {
            change_movement();
            create_legend();
        });

        function change_fill() {
            fill_color = $(containerSelector + " #snakeplot_color").val();
            
            //color_filtered = d3.select(containerSelector).select("#fill_filtered").property("checked");
            color_filtered = ($(containerSelector + " #fill_filtered").find(".active").attr('value') == 'true');
            console.log('change fill color to!', fill_color,'color_filtered',color_filtered);

            // color_id1  = $(containerSelector+" #fill_color1").val();
            // color_id2  = $(containerSelector+" #fill_color2").val();
            // color_id3  = $(containerSelector+" #fill_color3").val();
            color_id1  = $(containerSelector+" #fill_color1").spectrum("get").toHexString();
            if ($(containerSelector + " #fill_color2").spectrum("get") && $(containerSelector + " #fill_color2").length) {
                color_id2 = $(containerSelector + " #fill_color2").spectrum("get").toHexString();
            } else {
                color_id2 = false;
            }
            if ($(containerSelector + " #fill_color3").spectrum("get") && $(containerSelector + " #fill_color3").length) {
                color_id3 = $(containerSelector + " #fill_color3").spectrum("get").toHexString();
            } else {
                color_id3 = false;
            }
            console.log('change fill color to!', fill_color, color_id1, color_id2, color_id3);
            // if no color 3 make only linear between two colors.
            if (color_id3) {
                var color_range = d3v4.scaleLinear()
                    .domain([0, 0.5, 1])
                    .range([color_id1, color_id2, color_id3]);
            } else {
                var color_range = d3v4.scaleLinear()
                    .domain([0, 1])
                    .range([color_id1, color_id2]);
            }
            $(containerSelector).find('.rcircle').each(function () {

                pos_id = $(this).attr('id');
                original_title = $(this).attr('original_title');
                gn = false;
                if (original_title.split(" ")[1].length) {
                    // this has GN
                    seg = original_title.split(" ")[1].split("x")[0];
                    gn = original_title.split(" ")[1].split("x")[1];
                }

                if (fill_color == 'none') {
                    $(this).attr("fill", "#fff");
                    return true;
                }

                if (color_filtered && (!gn || !filtered_gns.includes(seg+"x"+gn))) {
                    $(this).attr("fill", "#fff");
                    // $(this).attr("stroke", "#ccc");  
                    // $(this).attr("stroke-width", 1); 
                    // $(containerSelector).find('#' + pos_id + 't').attr("fill", "#ddd");
                    return true;
                }

                if (pos_id in colors[fill_color]) {
                    var scale = colors[fill_color][pos_id][1];
                    if (fill_color == 'network') {
                        // Network is more of a 'categorical' color scale, so needs different code
                        $(this).attr("fill", color_by_category(colors[fill_color][pos_id][0],colors[fill_color][pos_id][2]));
                    } else {
                        $(this).attr("fill", color_range(scale)); //color_by_scale(scale,color_id1,color_id2,color_id3));
                    }
                } else {


                    if (fill_color.includes("distance")) {
                        closest = find_closest_value(seg, gn, fill_color);
                        //console.log('found closest', closest, seg, gn);
                        var scale = closest[1];
                    } else {
                        var closest = 'not_found';
                    }
                    if (closest == 'not_found') {
                        $(this).attr("fill", "#fff"); 
                    } else {
                        $(this).attr("fill", color_range(scale)); //color_by_scale(scale,color_id1,color_id2,color_id3));
                    }
                }
                // console.log(pos_id,color[1],color_by_scale(color_id,color[1]))
                // $(this).attr("fill_value", color[2]);  

                if ($(containerSelector + " #snakeplot_color_text").val() == 'none') {
                    $(containerSelector).find('#'+pos_id+'t').removeAttr("fill");
                    $(containerSelector).find('#' + pos_id + 't').attr("fill", "#000");
                }

                $(this).css("opacity", 1);
                $(this).removeAttr("fill-opacity");
                $(containerSelector).find('#'+pos_id+'t').css("opacity", 1);
                // if (color[2] < 0.05 && $(this).attr("stroke_value")<0.05) {
                //     $(this).css("opacity", 0.3);
                //     $(containerSelector).find('#'+pos_id+'t').css("opacity", 0.1);
                // }
                             
            });
        }

        function change_stroke() {
            fill_color = $(containerSelector + " #snakeplot_color_border").val();
            // color_id1  = $(containerSelector+" #border_color1").val();
            // color_id2  = $(containerSelector+" #border_color2").val();
            // color_id3  = $(containerSelector+" #border_color3").val();
            color_filtered = d3.select(containerSelector).select("#border_filtered").property("checked");
            color_filtered = ($(containerSelector + " #border_filtered").find(".active").attr('value') == 'true');
            var stroke_width = $(containerSelector + " #snakeplot_border_stroke").val();
            console.log('change stroke color to!');
            color_id1  = $(containerSelector+" #border_color1").spectrum("get").toHexString();
            if ($(containerSelector + " #border_color2").spectrum("get") && $(containerSelector + " #border_color2").length) {
                color_id2 = $(containerSelector + " #border_color2").spectrum("get").toHexString();
            } else {
                color_id2 = false;
            }
            if ($(containerSelector + " #border_color3").spectrum("get") && $(containerSelector + " #border_color3").length) {
                color_id3 = $(containerSelector + " #border_color3").spectrum("get").toHexString();
            } else {
                color_id3 = false;
            }
            // if no color 3 make only linear between two colors.
            if (color_id3) {
                var color_range = d3v4.scaleLinear()
                    .domain([0, 0.5, 1])
                    .range([color_id1, color_id2, color_id3]);
            } else {
                var color_range = d3v4.scaleLinear()
                    .domain([0, 1])
                    .range([color_id1, color_id2]);
            }

            $(containerSelector).find('.rcircle').each(function () {
                original_title = $(this).attr('original_title');
                gn = false;
                pos_id = $(this).attr('id');
                if (original_title.split(" ")[1].length) {
                    // this has GN
                    seg = original_title.split(" ")[1].split("x")[0];
                    gn = original_title.split(" ")[1].split("x")[1];
                }


                if (fill_color == 'none') {
                    $(this).attr("stroke", "#ccc");
                    return true;
                }

                if (color_filtered && (!gn || !filtered_gns.includes(seg+"x"+gn))) {
                    // $(this).attr("fill", "#fff");
                    $(this).attr("stroke", "#ccc");  
                    $(this).attr("stroke-width", 1); 
                    // $(containerSelector).find('#' + pos_id + 't').attr("fill", "#ddd");
                    return true;
                }

                if (pos_id in colors[fill_color]) {
                    var scale = colors[fill_color][pos_id][1];
                    if (fill_color == 'network') {
                        // Network is more of a 'categorical' color scale, so needs different code
                        $(this).attr("stroke", color_by_category(colors[fill_color][pos_id][0],colors[fill_color][pos_id][2]));
                    } else {
                        $(this).attr("stroke", color_range(scale));
                    }
                    $(this).attr("stroke-width", stroke_width); 
                } else {
                    if (fill_color.includes("distance")) {
                        closest = find_closest_value(seg, gn, fill_color);
                        //console.log('found closest', closest, seg, gn);
                        var scale = closest[1];
                    } else {
                        var closest = 'not_found';
                    }
                    if (closest == 'not_found') {
                        $(this).attr("stroke", "#ccc");
                        $(this).attr("stroke-width", 1); 
                    } else {
                        $(this).attr("stroke", color_range(scale));
                        $(this).attr("stroke-width", stroke_width); 
                    }
                }
                
                
                // $(this).attr("stroke", color[color_id]);
                // $(this).attr("stroke_value", color[2]);          
                // $(this).attr("stroke-width", 1); 
                // // $(this).css("opacity", 0.3);
                // if (color[2] > 0.1) {    
                //     $(this).attr("stroke-width", 3); 
                //     $(this).css("opacity", 1);
                // } else {
                //     $(this).attr("stroke", "#ccc");
                // }
            });
        }


        function change_text() {
            fill_color = $(containerSelector + " #snakeplot_color_text").val();
            // color_id1  = $(containerSelector+" #text_color1").val();
            // color_id2  = $(containerSelector+" #text_color2").val();
            // color_id3  = $(containerSelector+" #text_color3").val();
            color_filtered = d3.select(containerSelector).select("#text_filtered").property("checked");
            color_filtered = ($(containerSelector + " #text_filtered").find(".active").attr('value') == 'true');
            console.log('change text color to!', fill_color, 'color_filtered', color_filtered);
            color_id1  = $(containerSelector+" #text_color1").spectrum("get").toHexString();
            if ($(containerSelector + " #text_color2").spectrum("get") && $(containerSelector + " #text_color2").length) {
                color_id2 = $(containerSelector + " #text_color2").spectrum("get").toHexString();
            } else {
                color_id2 = false;
            }
            if ($(containerSelector + " #text_color3").spectrum("get") && $(containerSelector + " #text_color3").length) {
                color_id3 = $(containerSelector + " #text_color3").spectrum("get").toHexString();
            } else {
                color_id3 = false;
            }
            console.log('change fill color to!', fill_color, color_id1, color_id2, color_id3);
            // if no color 3 make only linear between two colors.
            if (color_id3) {
                var color_range = d3v4.scaleLinear()
                    .domain([0, 0.5, 1])
                    .range([color_id1, color_id2, color_id3]);
            } else {
                var color_range = d3v4.scaleLinear()
                    .domain([0, 1])
                    .range([color_id1, color_id2]);
            }
            console.log(color_id1, color_id2, color_id3);

            $(containerSelector).find('.rcircle').each(function () {
                original_title = $(this).attr('original_title');
                gn = false;
                pos_id = $(this).attr('id');
                if (original_title.split(" ")[1].length) {
                    // this has GN
                    seg = original_title.split(" ")[1].split("x")[0];
                    gn = original_title.split(" ")[1].split("x")[1];
                }

                if (fill_color == 'none') {
                    $(containerSelector).find('#' + pos_id + 't').attr("fill", "#000");
                    $(containerSelector).find('#' + pos_id + 't').attr("font-weight", 0);
                    return true;
                }

                if (color_filtered && (!gn || !filtered_gns.includes(seg+"x"+gn))) {
                    $(containerSelector).find('#' + pos_id + 't').attr("fill", "#ddd");
                    $(containerSelector).find('#' + pos_id + 't').attr("font-weight", 0);
                    return true;
                }


                if (pos_id in colors[fill_color]) {
                    var scale = colors[fill_color][pos_id][1];
                    if (fill_color == 'network') {
                        // Network is more of a 'categorical' color scale, so needs different code
                        $(containerSelector).find('#' + pos_id + 't').attr("fill", color_by_category(colors[fill_color][pos_id][0],colors[fill_color][pos_id][2]));
                    } else {
                        $(containerSelector).find('#' + pos_id + 't').attr("fill", color_range(scale));
                    }
                    $(containerSelector).find('#' + pos_id + 't').attr("font-weight", 1000);
                    // $(containerSelector).find('#' + pos_id + 't').attr("stroke", "#000");
                    // $(containerSelector).find('#' + pos_id + 't').attr("stroke-width", 1);
                } else {

                    if (fill_color.includes("distance")) {
                        closest = find_closest_value(seg, gn, fill_color);
                        //console.log('found closest', closest, seg, gn);
                        var scale = closest[1];
                    } else {
                        var closest = 'not_found';
                    }
                    if (closest == 'not_found') {
                        $(containerSelector).find('#' + pos_id + 't').attr("fill", "#000");
                        $(containerSelector).find('#' + pos_id + 't').attr("font-weight", 0);
                    } else {
                        $(containerSelector).find('#' + pos_id + 't').attr("fill", color_range(scale));
                        $(containerSelector).find('#' + pos_id + 't').attr("font-weight", 1000);
                    }
                    
                }
                
                // $(containerSelector).find('#' + pos_id + 't').attr("fill", color[color_id]);
                // $(this).css("opacity", 0.3);
                // if (color[2] > 0.1) {    
                    
                // } else {
                //     $(containerSelector).find('#' + pos_id + 't').attr("fill", "#ddd");
                // }
            });
        }

        function change_rotation() {
            fill_color = $(containerSelector + " #snakeplot_color_rotation").val();
            color_filtered = ($(containerSelector + " #rotation_filtered").find(".active").attr('value') == 'true');
            console.log('change rotation to!', fill_color, 'color_filtered', color_filtered);

            $(containerSelector).find('.rcircle').each(function () {
                original_title = $(this).attr('original_title');
                gn = false;
                pos_id = $(this).attr('id');
                if (original_title.split(" ")[1].length) {
                    // this has GN
                    seg = original_title.split(" ")[1].split("x")[0];
                    gn = original_title.split(" ")[1].split("x")[1];
                }

                if (fill_color == 'none') {
                    $(containerSelector).find('#' + pos_id + 't').attr("transform", "");
                    return true;
                }

                if (color_filtered && (!gn || !filtered_gns.includes(seg+"x"+gn))) {
                    $(containerSelector).find('#' + pos_id + 't').attr("transform", "");
                    return true;
                }


                if (pos_id in colors[fill_color]) {
                    var scale = colors[fill_color][pos_id][1];
                    var rotation_value = Math.round(scale * 180 - 90);
                    if (fill_color == 'outer_angle' || fill_color == 'rotation_angle') {
                        // For these 'actual' rotation values, use the 'true value' as the rotation.
                        var rotation_value = Math.round(colors[fill_color][pos_id][0]);
                    }
                    var x = $(containerSelector).find('#' + pos_id + 't').attr("x");
                    var y = $(containerSelector).find('#' + pos_id + 't').attr("y");
                    // d3.select("#t37t").transition().duration(10000).attr("transform","rotate(45)")
                    $(containerSelector).find('#' + pos_id + 't').attr("transform", "rotate("+rotation_value+" " + x + " " + y + ")");
                } else {
                    $(containerSelector).find('#' + pos_id + 't').attr("transform", "");
                }
                
                // $(containerSelector).find('#' + pos_id + 't').attr("fill", color[color_id]);
                // $(this).css("opacity", 0.3);
                // if (color[2] > 0.1) {    
                    
                // } else {
                //     $(containerSelector).find('#' + pos_id + 't').attr("fill", "#ddd");
                // }
            });
        }

        function change_backbone() {
            fill_color = $(containerSelector + " #snakeplot_color_backbone").val();
            // color_id1  = $(containerSelector+" #backbone_color1").val();
            // color_id2  = $(containerSelector+" #backbone_color2").val();
            // color_id3  = $(containerSelector+" #backbone_color3").val();
            color_id1  = $(containerSelector+" #backbone_color1").spectrum("get").toHexString();
            if ($(containerSelector + " #backbone_color2").spectrum("get") && $(containerSelector + " #backbone_color2").length) {
                color_id2 = $(containerSelector + " #backbone_color2").spectrum("get").toHexString();
            } else {
                color_id2 = false;
            }
            if ($(containerSelector + " #backbone_color3").spectrum("get") && $(containerSelector + " #backbone_color3").length) {
                color_id3 = $(containerSelector + " #backbone_color3").spectrum("get").toHexString();
            } else {
                color_id3 = false;
            }
            // if no color 3 make only linear between two colors.
            if (color_id3) {
                var color_range = d3v4.scaleLinear()
                    .domain([0, 0.5, 1])
                    .range([color_id1, color_id2, color_id3]);
            } else {
                var color_range = d3v4.scaleLinear()
                    .domain([0, 1])
                    .range([color_id1, color_id2]);
            }
            console.log('change backbone color to!', fill_color);
            var path_max = 0;
            $(containerSelector).find('.helix_path').each(function () {
                path_id = $(this).attr('id');
                fill_scale = path_groups[path_id][fill_color];
                const fill_sum = fill_scale.reduce((a, b) => a + b, 0);
                const fill_avg = (fill_sum / fill_scale.length) || 0;
                path_max = fill_avg > path_max ? fill_avg : path_max;
            });

            $(containerSelector).find('.helix_path').each(function () {
                path_id = $(this).attr('id');

                // if (color_filtered && (!gn || !filtered_gns.includes(seg+"x"+gn))) {
                //     $(containerSelector).find('#' + pos_id + 't').attr("fill", "#ddd");
                //     return true;
                // }

                if (fill_color == 'none') {
                    $(this).attr("stroke", "grey");
                    $(this).attr("stroke-width", 2);
                    return
                }

                fill_scale = path_groups[path_id][fill_color];

                const fill_sum = fill_scale.reduce((a, b) => a + b, 0);
                const fill_avg = (fill_sum / fill_scale.length) || 0;
                if (fill_color != 'none' && fill_scale.length!=0) {
                        $(this).attr("stroke", color_range(fill_avg));
                        $(this).attr("stroke-width", 6);
                } 
                
            });
        }

        function find_closest_value(seg, gn, fill_color) {
            color = "not_found";   
            if (gn && parseInt(gn.substring(0,2)) >= 50) {
                for (i = parseInt(gn.substring(0,2)); i > 0; i--) {
                    seq_pos = data['snakeplot_lookup'][seg+"x"+i]
                    if (seq_pos in colors[fill_color]) {
                        color = colors[fill_color][seq_pos];
                        break;
                    }
                }
            } else if (gn) {
                for (i = parseInt(gn.substring(0,2)); i < 100; i++) {
                    seq_pos = data['snakeplot_lookup'][seg+"x"+i]
                    if (seq_pos in colors[fill_color]) {
                        color = colors[fill_color][seq_pos];
                        break;
                    }
                }
            }
            return color;
        }

        function change_movement() {
            fill_color = $(containerSelector + " #snakeplot_move_circle").val();
            console.log('change movement to!', fill_color);
            
            // color_filtered = d3.select(containerSelector).select("#color_filtered").property("checked");

            $(containerSelector).find('.rcircle').each(function () {

                pos_id = $(this).attr('id');
                original_title = $(this).attr('original_title');
                gn = false;
                if (original_title.split(" ")[1].length) {
                    // this has GN
                    seg = original_title.split(" ")[1].split("x")[0];
                    gn = original_title.split(" ")[1].split("x")[1];
                }
                //if (color_filtered && (!gn || !filtered_gns.includes(seg+"x"+gn))) {
                    // return true;
                //}

                if (fill_color != 'none') {
                    if (pos_id in colors[fill_color]) {
                        color = colors[fill_color][pos_id];
                    } else {
                        console.log("no movement for ", seg,gn, pos_id);
                        // find closest value
                        
                        color = find_closest_value(seg, gn, fill_color);
                        if (color == 'not_found')
                            color = ["#fff", 0, 0, "#fff", "#fff"];   
                    }
                    //color = pos_id in colors[fill_color] ? colors[fill_color][pos_id] : ["#fff",0,0,"#fff","#fff"];
                } else {
                    color = ["#fff", 0, 0, "#fff", "#fff"];   
                }

                max_momement = 25;

                current_x = parseInt($(this).attr("original_cx"));
                current_x_text = parseInt($(containerSelector).find('#' + pos_id + 't').attr("original_x"));
                $(this).attr("cx", current_x + color[1] * max_momement);
                $(containerSelector).find('#' + pos_id + 't').attr("x",current_x_text + color[1] * max_momement)
                $(containerSelector).find('#' + pos_id + 't').find("tspan").attr("x",current_x_text + color[1] * max_momement)
                             
            });
        }

    }

    function create_legend() {
        console.log('create legend!');
        // var newDiv = document.createElement("div");

        // $(containerSelector).find(".snakeplot-legend").remove();

        // newDiv.setAttribute("class", "snakeplot-legend");
        // newDiv.innerHTML = 'test legend';
        // $(containerSelector).append(newDiv);
        var fill_color;
        legends = [];
        // Deduce which legends to make
        fill_color = $(containerSelector + " #snakeplot_color").val();
        if (fill_color != "none" && fill_color) {

            var color1  = $(containerSelector+" #fill_color1").spectrum("get").toHexString();
            var color2  = $(containerSelector+" #fill_color2").spectrum("get").toHexString();
            var color3  = $(containerSelector+" #fill_color3").val() == "" ? 'none' : $(containerSelector+" #fill_color3").spectrum("get").toHexString();
            var colors = [color1,color2,color3].filter(item => !(item == 'none'));
            legends.push({ icon: 'fill', value: nice_index_names[fill_color],colors:colors })
        }
        fill_color = $(containerSelector + " #snakeplot_color_border").val();
        if (fill_color!="none" && fill_color) {

            var color1  = $(containerSelector+" #border_color1").val();
            var color2  = $(containerSelector+" #border_color2").val();
            var color3  = $(containerSelector+" #border_color3").val();
            var colors = [color1,color2,color3].filter(item => !(item == ''));
            legends.push({ icon: 'border', value: nice_index_names[fill_color],colors:colors })
        }
        fill_color = $(containerSelector + " #snakeplot_color_text").val();
        if (fill_color!="none" && fill_color) {

            var color1  = $(containerSelector+" #text_color1").val();
            var color2  = $(containerSelector+" #text_color2").val();
            var color3  = $(containerSelector+" #text_color3").val();
            var colors = [color1, color2, color3].filter(item => !(item == ''));
            legends.push({ icon: 'text', value: nice_index_names[fill_color],colors:colors })
        }
        fill_color = $(containerSelector + " #snakeplot_color_rotation").val();
        if (fill_color!="none" && fill_color) legends.push({ icon: 'rotation', value: nice_index_names[fill_color]})
        fill_color = $(containerSelector + " #snakeplot_color_backbone").val();
            if (fill_color != "none" && fill_color) {

                var color1 = $(containerSelector + " #backbone_color1").val();
                var color2 = $(containerSelector + " #backbone_color2").val();
                var color3 = $(containerSelector + " #backbone_color3").val();
                var colors = [color1, color2, color3].filter(item => !(item == ''));
                legends.push({ icon: 'backbone_line', value: nice_index_names[fill_color],colors:colors })
            }
        fill_color = $(containerSelector + " #snakeplot_move_circle").val();
        if (fill_color!="none" && fill_color) legends.push({ icon: 'backbone', value: nice_index_names[fill_color]})

        console.log(legends);
        // legends = [{ icon: 'fill', value: 'rotamer' },
        //            { icon: 'border', value: 'distance' },
        //            { icon: 'text', value: 'distance' },
        //            { icon: 'rotation', value: 'distance' },
        //            { icon: 'backbone_line', value: 'distance' },
        //            { icon: 'backbone', value: 'distance' }];

        var dataL = 0;
        var offset = 45;
        d3.select(containerSelector).select(".legend").remove();
        var legend = d3.select(containerSelector).select("svg").append("g").attr("class","legend");
        legend.attr("transform", "translate(15,"+(newheight-70)+")");
        var legend4 = legend.selectAll('.legends')
            .data(legends)
            .enter().append('g')
            .attr("transform", function (d, i) {
                if (i === 0) {
                    dataL = Math.max(50,d.value.length*6) + offset
                    return "translate(0,0)"
                } else {
                    var newdataL = dataL
                    console.log(d.value,"length",d.value.length)
                    dataL += Math.max(50,d.value.length*6)  + offset
                    return "translate(" + (newdataL) + ",0)"
                }
            })

        legend4.append('text')
            .attr("x", 25)
            .attr("y", function (d) { return 'colors' in d ? 1 : 6 })
            .text(function (d, i) {return d.value})
            .style("text-anchor", "start")
            .style("font-size", 15)

        legend4.append("svg:image")
            .attr('x', 0)
            .attr('y', -12)
            .attr('width', 20)
            .attr('height', 24)
            .attr("xlink:href", function (d, i) { return "/static/home/images/legends/" + d.icon + ".png" })
        
        legendcolorscales = legend4.append('g').attr('class','colorscale')
        legendcolorscales.append('rect')
            .attr('x', 25)
            .attr('y', 3)
            .attr('width', 50)
            .attr('height', 8)
            .style("stroke", "black")
            .style("stroke-width", 1)
            .style("fill", "white")
        
        // console.log(legendcolorscales)

        // for (element in legendcolorscales[0]) {
        //     console.log("test",element)
        // }
        legendcolorscales.each(function(d) {
            // your update code here as it was in your example
            test = d3.select(this) // Transform to d3 Object
            if (!('colors' in d)) {
                test.remove();
                return;
            }
            var data_legend = d3.range(48);
            range_colors = []
            d.colors.forEach(c => {
                // var rgb = color_scale_colors[c]
                // var hex = rgb2hex(rgb.red, rgb.green, rgb.blue);
                range_colors.push(c);
            });
            
            if (d.value == nice_index_names['network']) {
                var colors = d3v4.scaleSequential(d3v4.interpolateSpectral)
                    .domain([0, 48]);
            } else {
                if (range_colors.length == 2) {
                    var colors = d3v4.scaleLinear()
                        .domain([0,48])
                        .range(range_colors);
                } else { // three colors
                    var colors = d3v4.scaleLinear()
                        .domain([0,24,48])
                        .range(range_colors);
                }
            }
            console.log(colors, range_colors);
            var rects = test.selectAll(".colorinterval")
                .data(data_legend)
                .enter()
                .append("rect")
                .attr("y", 4)
                .attr("height", 6)
                .attr("x", (d,i)=>26 + i)
                .attr("width", 1)
                .attr("fill", d => colors(d))
                .attr("id", d => d)
                .attr("stroke-width", 0);
                            
          });


        legends_width = legend.node().getBBox().width
        legend.attr("transform", "translate("+((max_x-legends_width)/2)+","+(newheight-70)+")");
        
        // legend4.append('rect')
        //     .attr("x", 0)
        //     .attr("y", 0)
        //     .attr("width", 10)
        //     .attr("height", 10)
        //     .style("fill", function (d, i) {
        //         return color(i)
        //     })
        
        // legend4.append('text')
        //     .attr("x", 20)
        //     .attr("y", 10)
        //     //.attr("dy", ".35em")
        //     .text(function (d, i) {
        //         return d
        //     })
        //     .attr("class", "textselected")
        //     .style("text-anchor", "start")
        //     .style("font-size", 15)
            
        }
    
}

// Override "normal max min"
function maxmin() {

    if (!$('#snake').length) return

    $('.snakeplot').each(function(i, obj) {
        // console.log('fix maxmin of obj!',obj,$(obj).find('svg').children('.rtext'));
        margin = 70;
        svgmax = 0;
        svgmin = 0;
        x_svgmax = 0;
        count = 0;
        classmax = '';
        classmin = '';
        counter = 0;
    
        // console.log("temp",y_max,y_min);
        $(obj).find('svg').find('g').children('.rtext').each(function () {
            counter += 1;
            y = parseInt($(this).attr( "y" ));
            x = parseInt($(this).attr( "x" ));
            classtext = $(this).attr("class");
            // test = $(this).attr("original_title");
            // test2 = $(this).css("display");
            if ($(this).css("display")!='none') {
                count = count +1;
                if (y<svgmin) {
                    svgmin = y;
                    classmin = classtext;
                    }
                if (y>svgmax) {
                    classmax = classtext;
                    svgmax= y;
                }
                if (x>x_svgmax) x_svgmax = x;

            }
        });
        // console.log("temp",y_max,y_min);
        $(obj).find('svg').find('g').children('.segment').each(function () {
            counter += 1;
            y = parseInt($(this).attr( "y" ));
            x = parseInt($(this).attr( "x" ));
            classtext = $(this).attr("class");
            // test = $(this).attr("original_title");
            // test2 = $(this).css("display");
            if ($(this).css("display")!='none') {
                count = count +1;
                if (y<svgmin) {
                    svgmin = y;
                    classmin = classtext;
                    }
                if (y>svgmax) {
                    classmax = classtext;
                    svgmax= y;
                }
                if (x>x_svgmax) x_svgmax = x;

            }
        });
        // if (svgmin>y_min) svgmin = y_min;
        // if (svgmax<y_max) svgmax = y_max;

        // console.log('max '+svgmax+' '+classmax+' min'+svgmin+' '+classmin+' count'+count);

        var svg = $(obj).find('svg');

        check = svg.attr('viewBox');
        // console.log('check',check)
        if (typeof check !== typeof undefined && check !== false && check !== null ) {
        oldheight = check.split(" ")[3];
        width = check.split(" ")[2];
        } else {
        // console.log('not found it');
        oldheight = $(svg).attr('height');
        width = $(svg).attr('width');
        }

        newheight = (svgmax-svgmin+margin*2);

        // console.log('New height:'+ newheight +' old height:'+oldheight);
        // console.log("Prev attr"+$('#snake').attr("transform"));
        if (newheight!=oldheight) {
            svg.attr('height', (svgmax-svgmin+margin*2));
            svg.find('g#snake').attr("transform", "translate(0," + (-svgmin+margin/2) + ")");
            svg.find('g.legend').attr("transform", "translate(15," + (newheight-70) + ")");

            // $('#snakeplot')[0].attr("viewBox", "0 0 " + width + " " + newheight);
            svg.attr("viewBox", "0 0 " + width + " " + newheight);

            svg.attr('height', "100%");
            svg.attr('width', "100%");
        }
        
    });
}
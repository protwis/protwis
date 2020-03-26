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

    var colors = {}

    colors.red = {red:255, green:0, blue: 0};
    colors.red = {red:195, green:74, blue: 54};
    colors.blue = {red:0, green:0, blue: 255};
    colors.blue = {red:0, green:140, blue: 204};
    colors.green = {red:0, green:255, blue: 0};
    colors.green = {red:0, green:201, blue: 167};
    colors.white = {red:255, green:255, blue: 255};
    colors.yellow = {red:255, green:255, blue: 0};
    colors.yellow = {red:255, green:255, blue: 0};
    colors.black = { red: 0, green: 0, blue: 0 };
    colors.orange = { red: 255, green: 150, blue: 113 };
    colors.purple = { red: 128, green: 0, blue: 128 };
    colors.false = false;
    
    return colorGradient(scale, colors[color1], colors[color2], colors[color3])
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
    // console.log(path_groups, path_groups_lookup);
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
        'w': 'white',
        'r': 'red',
        'b': 'blue',
        'black': 'black',
        'g': 'green',
        'y': 'yellow',
        'o': 'orange',
        'p': 'purple',
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
    
        content += '<table><tr><td>Residue text (kept positions):</td><td><select id="text_included">' +
        text_content_options +
        '</select>' + 
        '<div class="btn-group btn-toggle" id="included_color">' +
        '    <button class="btn btn-xs btn-primary active" value="black">Black</button>' +
        '    <button class="btn btn-xs btn-default" value="grey">Grey</button>' +
            '</div>' +
        '</td ></tr > <tr><td>Residue text (filtered out positions):</td><td><select id="text_excluded">' +
        text_content_options +
        '</select>' + 
        '<div class="btn-group btn-toggle" id="excluded_color">' +
        '    <button class="btn btn-xs btn-primary active" value="black">Black</button>' +
        '    <button class="btn btn-xs btn-default" value="grey">Grey</button>' +
            '</div>' +
        '</td ></tr ></table > '
        ;

        content += '<table><tr><th>Area</th><th>Property</th><th>Color1</th><th>Color2</th><th>Color3</th><th>Positions</th></tr>' +
            '<tr><td>Residue fill:</td><td><select id="snakeplot_color" class="residue_fill">' +
            '<option value="none">None</option>' +
            select_data_options +
            '</select></td><td>' +
            '<select id=fill_color1 class="fill_color residue_fill">' +
            select_color_options_white +
            '</select></td><td>' +
            '<select id=fill_color2 class="fill_color residue_fill">' +
            select_color_options_red +
            '</select></td><td>' +
            '<select id=fill_color3 class="fill_color residue_fill">' +
            '<option value="none">None</option>' +
            select_color_options +
            '</select></td>' +
            '<td>' + 
            '<div class="btn-group btn-toggle residue_fill" id="fill_filtered">' +
            '    <button class="btn btn-xs btn-primary active" value="true">Kept</button>' +
            '    <button class="btn btn-xs btn-default" value="false">All</button>' +
            '</div>' +
            '</td></tr > '
            ;
        content += '<tr><td>Residue border:</td><td><select id="snakeplot_color_border" class="residue_border">' +
            '<option value="none">None</option>' +
            select_data_options +
            '</select></td><td>' +
            '<select id=border_color1 class="border_color residue_border">' +
            select_color_options_white +
            '</select></td><td>' +
            '<select id=border_color2 class="border_color residue_border">' +
            select_color_options_red +
            '</select></td><td>' +
            '<select id=border_color3 class="border_color residue_border">' +
            '<option value="none">None</option>' +
            select_color_options +
            '</select></td>' +
            '<td>' + 
            '<div class="btn-group btn-toggle residue_border" id="border_filtered">' +
            '    <button class="btn btn-xs btn-primary active" value="true">Kept</button>' +
            '    <button class="btn btn-xs btn-default" value="false">All</button>' +
            '</div>' +
            '</td></tr> '
            ;
        content += '<tr><td>Residue text:</td><td><select id="snakeplot_color_text" class="residue_text">' +
                '<option value="none">None</option>' +
                select_data_options +
                '</select></td><td>' +
                '<select id=text_color1 class="text_color residue_text">' +
                select_color_options_white +
                '</select></td><td>' +
                '<select id=text_color2 class="text_color residue_text">' +
                select_color_options_red +
                '</select></td><td>' +
                '<select id=text_color3 class="text_color residue_text">' +
                '<option value="none">None</option>' +
                select_color_options +
                '</select></td>' +
                '<td>' + 
                '<div class="btn-group btn-toggle residue_text" id="text_filtered">' +
                '    <button class="btn btn-xs btn-primary active" value="true">Kept</button>' +
                '    <button class="btn btn-xs btn-default" value="false">All</button>' +
                '</div>' +
                '</td></tr> '
                ;
        content += '<tr><td>Residue rotation:</td><td colspan=4><select id="snakeplot_color_rotation" class="residue_rotation">' +
            '<option value="none">None</option>' +
            select_data_options +
            '</select></td>' +
            '<td>' + 
            '<div class="btn-group btn-toggle residue_rotation" id="rotation_filtered">' +
            '    <button class="btn btn-xs btn-primary active" value="true">Kept</button>' +
            '    <button class="btn btn-xs btn-default" value="false">All</button>' +
            '</div>' +
            '</td></tr> '
        content += '<tr><td colspan=6><hr></td></tr><tr><td>Backbone line:</td><td><select id="snakeplot_color_backbone">' +
                '<option value="none">None</option>' +
                select_data_options_backbone +
                '</select></td><td>' +
                '<select id=backbone_color1 class=backbone_color>' +
                select_color_options +
                '</select></td><td>' +
                '<select id=backbone_color2 class=backbone_color>' +
                select_color_options +
                '</select></td><td>' +
                '<select id=backbone_color3 class=backbone_color>' +
                '<option value="none">None</option>' +
                select_color_options +
                '</select></td></tr>'
                ;
        content += '<tr><td>Backbone shift:</td><td><select id="snakeplot_move_circle">' +
                '<option value="none">None</option>' +
                select_data_options_backbone_shift +
                '</select></td><td>' +
                '</td></tr></table>'
                ;
        content += '</span>';
        newDiv.innerHTML = content;

        $(containerSelector).prepend(newDiv);
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

        d3.select(containerSelector).select("#text_included").on("change", function () {
            text_format = $(containerSelector + " #text_included").val();
            color = $(containerSelector + " #included_color").find(".active").attr('value');
            console.log('change to', text_format);
            $(containerSelector).find('.rtext').each(function () {
                original_title = $(this).attr('original_title');
                if (original_title.split(" ")[1].length) {
                    // this has GN
                    seg = original_title.split(" ")[1].split(".")[0];
                    gn = original_title.split(" ")[1].split("x")[1];
                    AA = original_title[0];
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
        });

        d3.select(containerSelector).select("#text_excluded").on("change", function () {
            
            color = $(containerSelector + " #excluded_color").find(".active").attr('value');
            text_format = $(containerSelector + " #text_excluded").val();
            console.log('change to', text_format);
            $(containerSelector).find('.rtext').each(function () {
                original_title = $(this).attr('original_title');
                if (original_title.split(" ")[1].length) {
                    seg = original_title.split(" ")[1].split(".")[0];
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
        });


        $(containerSelector+ " .residue_fill").on("change", function () {
            change_fill();
        });

        $(containerSelector+ " .residue_border").on("change", function () {
            change_stroke();
            console.log('residue_border');
        });


        $(containerSelector+ " .residue_text").on("change", function () {
            change_text();
        });


        $(containerSelector+ " .backbone_color").on("change", function () {
            change_backbone();
        });

        $(containerSelector+ " .residue_rotation").on("change", function () {
            change_rotation();
        });

        d3.select(containerSelector).select("#snakeplot_color_backbone").on("change", function () {
            change_backbone();
        });

        d3.select(containerSelector).select("#snakeplot_move_circle").on("change", function () {
            change_movement();
        });

        function change_fill() {
            fill_color = $(containerSelector + " #snakeplot_color").val();
            
            //color_filtered = d3.select(containerSelector).select("#fill_filtered").property("checked");
            color_filtered = ($(containerSelector + " #fill_filtered").find(".active").attr('value') == 'true');
            console.log('change fill color to!', fill_color,'color_filtered',color_filtered);

            color_id1  = $(containerSelector+" #fill_color1").val();
            color_id2  = $(containerSelector+" #fill_color2").val();
            color_id3  = $(containerSelector+" #fill_color3").val();

            $(containerSelector).find('.rcircle').each(function () {

                pos_id = $(this).attr('id');
                original_title = $(this).attr('original_title');
                gn = false;
                if (original_title.split(" ")[1].length) {
                    // this has GN
                    seg = original_title.split(" ")[1].split(".")[0];
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
                    $(this).attr("fill", color_by_scale(scale,color_id1,color_id2,color_id3));
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
                        $(this).attr("fill", color_by_scale(scale,color_id1,color_id2,color_id3));
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
            color_id1  = $(containerSelector+" #border_color1").val();
            color_id2  = $(containerSelector+" #border_color2").val();
            color_id3  = $(containerSelector+" #border_color3").val();
            color_filtered = d3.select(containerSelector).select("#border_filtered").property("checked");
            color_filtered = ($(containerSelector + " #border_filtered").find(".active").attr('value') == 'true');
            console.log('change stroke color to!');

            $(containerSelector).find('.rcircle').each(function () {
                original_title = $(this).attr('original_title');
                gn = false;
                pos_id = $(this).attr('id');
                if (original_title.split(" ")[1].length) {
                    // this has GN
                    seg = original_title.split(" ")[1].split(".")[0];
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
                    $(this).attr("stroke", color_by_scale(scale,color_id1,color_id2,color_id3));
                    $(this).attr("stroke-width", 5); 
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
                        $(this).attr("stroke", color_by_scale(scale,color_id1,color_id2,color_id3));
                        $(this).attr("stroke-width", 5); 

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
            color_id1  = $(containerSelector+" #text_color1").val();
            color_id2  = $(containerSelector+" #text_color2").val();
            color_id3  = $(containerSelector+" #text_color3").val();
            color_filtered = d3.select(containerSelector).select("#text_filtered").property("checked");
            color_filtered = ($(containerSelector + " #text_filtered").find(".active").attr('value') == 'true');
            console.log('change text color to!', fill_color, 'color_filtered', color_filtered);
            console.log(color_id1, color_id2, color_id3);

            $(containerSelector).find('.rcircle').each(function () {
                original_title = $(this).attr('original_title');
                gn = false;
                pos_id = $(this).attr('id');
                if (original_title.split(" ")[1].length) {
                    // this has GN
                    seg = original_title.split(" ")[1].split(".")[0];
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
                    $(containerSelector).find('#' + pos_id + 't').attr("fill", color_by_scale(scale,color_id1,color_id2,color_id3));
                    $(containerSelector).find('#' + pos_id + 't').attr("font-weight", 1000);
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
                        $(containerSelector).find('#' + pos_id + 't').attr("fill", color_by_scale(scale,color_id1,color_id2,color_id3));
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
                    seg = original_title.split(" ")[1].split(".")[0];
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
            color_id1  = $(containerSelector+" #backbone_color1").val();
            color_id2  = $(containerSelector+" #backbone_color2").val();
            color_id3  = $(containerSelector+" #backbone_color3").val();
            console.log('change backbone color to!',fill_color);
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
                        $(this).attr("stroke", color_by_scale(fill_avg,color_id1,color_id2,color_id3));
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
                    seg = original_title.split(" ")[1].split(".")[0];
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
            svg.find('g').attr("transform", "translate(0," + (-svgmin+margin) + ")");

            // $('#snakeplot')[0].attr("viewBox", "0 0 " + width + " " + newheight);
            svg.attr("viewBox", "0 0 " + width + " " + newheight);

            svg.attr('height', "100%");
            svg.attr('width', "100%");
        }
        
    });
}
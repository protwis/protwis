function createSnakeplot(data,containerSelector) {


    $(containerSelector).html('')
    $(containerSelector).addClass("snakeplot");
    $(containerSelector).css("position","relative");

    // var svgContainer = d3v4.select(containerSelector).append("svg")
    //     .attr("viewBox", min_x + " " + min_y + " " + (max_x - min_x) + " " + (max_y - min_y))
    //     .attr("width", "100%")
    //     .attr("style", "height: 500px");
    
    // console.log('making snakeplot!')
    // console.log(containerSelector)
    // console.log(data['snakeplot'])

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

    $.each(data['distances'], function (gn, dis) {
        seq_pos = data['snakeplot_lookup'][gn];
        // console.log(gn, dis['avg'], seq_pos);
        value = dis['avg'];
        scale = Math.abs(value) / data['ngl_max_diff_distance'];
        var color = { r: 255, g: 255, b: 255 };
        if (value < 0) {
            // if the header is a set two, then make it red
            color = { r: 255, g: 255-(255-153)*scale, b: 255-(255-153)*scale }; //red
        } else if (value > 0) {
            // Positive numbers are blue either cos they are set 1 or cos "set 1 has most"
            // This is also used for single set/structure
            color = { r: 255-(255-153)*scale, g: 255-(255-204)*scale, b: 255 }; //blue
        }
        var hex = rgb2hex(color.r, color.g, color.b);
        // grey
        var color_grey = { r: 255*(1-scale), g: 255*(1-scale), b: 255*(1-scale) };
        var hex_grey = rgb2hex(color_grey.r, color_grey.g, color_grey.b);
        colors['distance'][seq_pos] = [hex,value,scale,hex_grey];
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
            if (value < 0) {
                // if the header is a set two, then make it red
                color = { r: 255, g: 255-(255-153)*scale, b: 255-(255-153)*scale }; //red
            } else if (value > 0) {
                // Positive numbers are blue either cos they are set 1 or cos "set 1 has most"
                // This is also used for single set/structure
                color = { r: 255-(255-153)*scale, g: 255-(255-204)*scale, b: 255 }; //blue
            }
            var hex = rgb2hex(color.r, color.g, color.b);
            // grey
            var color_grey = { r: 255*(1-scale), g: 255*(1-scale), b: 255*(1-scale) }; 
            var hex_grey = rgb2hex(color_grey.r, color_grey.g, color_grey.b);
            colors[index_names[i]][seq_pos] = [hex,value,scale,hex_grey];
            if (!(i in max_values)) max_values[i] = 0;
            // console.log(gn, i, a);
        });
    });
    // console.log(colors);


    
    create_overlay();

    function create_overlay() {
        var newDiv = document.createElement("div");

        $(containerSelector).find(".controls-panel").remove();

        newDiv.setAttribute("class", "controls-panel");
        content = '<span class="pull-right snakeplot_controls_toggle" style="cursor: pointer;"><span class="glyphicon glyphicon-option-horizontal btn-download png"></span></span><span class="options" style="display: block; min-width: 120px;">' +
            'Generic number instead of AA<input id="generic" type="checkbox"><br>' +
            'Hide low numbers <input id="hide_low" type="checkbox" checked><br>';
        
        index_names = { 0: 'core_distance', 1: 'a_angle', 2: 'outer_angle', 3: 'tau', 4: 'phi', 5: 'psi', 6: 'sasa', 7: 'rsa', 8: 'theta', 9: 'hse', 10: 'tau_angle' }
    
        content += 'Fill color: <select id="snakeplot_color">' +
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
            '</select><br>'
            ;
            content += 'Grey fill scale <input id="fill_grey" type="checkbox"><br>';
        content += 'Border color: <select id="snakeplot_color_border">' +
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
            '</select><br>'
            ;
        content += 'Grey border scale <input id="border_grey" type="checkbox" checked><br>';
        content += '</span>';
        newDiv.innerHTML = content;

        $(containerSelector).prepend(newDiv);
        $(containerSelector).find(".options").toggle();

        $(containerSelector).find(".snakeplot_controls_toggle").click(function() {
            $(containerSelector).find(".options").slideToggle();
        });

        d3.select(containerSelector).select("#generic").on("change", function () {
            generic = d3.select(containerSelector).select("#generic").property("checked");
            $(containerSelector).find('.rtext').each(function () {
                original_title = $(this).attr('original_title');
                if (original_title.split(" ")[1].length) {
                    // this has GN
                    gn = original_title.split(" ")[1].split("x")[1];
                    AA = original_title[0];
                } else {
                    gn = "-";
                    AA = original_title[0];
                }
                $(this).text(generic ? gn : AA);
            });
        });


        d3.select(containerSelector).select("#snakeplot_color").on("change", function () {
            change_fill();
        });

        d3.select(containerSelector).select("#fill_grey").on("change", function () {
            change_fill();
        });

        d3.select(containerSelector).select("#border_grey").on("change", function () {
            change_stroke();
        });

        d3.select(containerSelector).select("#snakeplot_color_border").on("change", function () {
            change_stroke();
        });

        function change_fill() {
            fill_color = $(containerSelector + " #snakeplot_color").val();
            console.log('change fill color to!', fill_color);

            $(containerSelector).find('.rcircle').each(function () {
                pos_id = $(this).attr('id');
                if (fill_color!='none') {
                    color = pos_id in colors[fill_color] ? colors[fill_color][pos_id] : ["#fff",0,0];
                } else {
                    color = ["#fff", 0, 0];     
                }
                fill_grey  = d3.select(containerSelector).select("#fill_grey").property("checked");
                
                $(this).attr("fill", fill_grey ? color[3] : color[0]);
                $(this).attr("fill_value", color[2]);  
                $(containerSelector).find('#'+pos_id+'t').removeAttr("fill");

                if (fill_grey && color[2]>0.5) {
                    $(containerSelector).find('#'+pos_id+'t').attr("fill", "#fff");
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
            console.log('change stroke color to!', fill_color);

            $(containerSelector).find('.rcircle').each(function () {
                pos_id = $(this).attr('id');
                if (fill_color!='none') {
                    color = pos_id in colors[fill_color] ? colors[fill_color][pos_id] : ["#fff",0,0];
                } else {
                    color = ["#ccc", 0, 0];     
                }
                fill_grey = d3.select(containerSelector).select("#border_grey").property("checked");
                
                
                $(this).attr("stroke", fill_grey ? color[3] : color[0]);
                $(this).attr("stroke_value", color[2]);          
                $(this).attr("stroke-width", 1); 
                // $(this).css("opacity", 0.3);
                if (color[2] > 0.1) {    
                    $(this).attr("stroke-width", 5); 
                    $(this).css("opacity", 1);
                } else {
                    $(this).attr("stroke", "#ccc");
                }
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
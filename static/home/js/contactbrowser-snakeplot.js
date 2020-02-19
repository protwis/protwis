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
        $(containerSelector).find('#'+seq_pos).css("fill", hex);

    });

    
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
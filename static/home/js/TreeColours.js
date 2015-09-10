
 function colours() {

    presetColors = {
        'D': ['#E60A0A', '#FDFF7B'], 'E': ['#E60A0A', '#FDFF7B'],
        'K': ['#145AFF', '#FDFF7B'], 'R': ['#145AFF', '#FDFF7B'],
        'S': ['#A70CC6', '#FDFF7B'], 'T': ['#A70CC6', '#FDFF7B'],
        'N': ['#A70CC6', '#FDFF7B'], 'Q': ['#A70CC6', '#FDFF7B'],
        'V': ['#E6E600', '#000000'], 'L': ['#E6E600', '#000000'],
        'I': ['#E6E600', '#000000'], 'A': ['#E6E600', '#000000'],
        'M': ['#E6E600', '#000000'], 'F': ['#18FF0B', '#000000'],
        'Y': ['#18FF0B', '#000000'], 'W': ['#0BCF00', '#000000'],
        'H': ['#0093DD', '#000000'], 'P': ['#CC0099', '#FDFF7B'],
        'C': ['#B2B548', '#000000'], 'G': ['#FF00F2', '#000000'],
        '-': ['#FFFFFF', '#000000'],
        '+': ['#FFFFFF', '#000000']
    };

    $('#svgCanvas').find('svg').attr('id', 'svgtree');
    $('#svgCanvas').find('svg').attr('class', 'svgtree');
    $('#svgCanvas').find('svg').removeAttr('xmlns:xlink');
    $('#svgCanvas').find('svg').removeAttr('xmlns');


    $(".bgfield").click(function () {
        parentid = $(this).closest('svg').attr('id');
        newcolor = $(".pick-color.selected").attr('id');
        newcolor = newcolor.split('-');
        $(this).css("fill", newcolor[1]);
    });
    $("[class^=chart]").click(function () {
        if ($(this).attr('id')) {
            id = $(this).attr('id');
            newcolor = $(".pick-color.selected").attr('id');
            newcolor = newcolor.split('-');
            $("[class^=chart]").each(function (index) {
                if ($(this).attr('id')) {
                    if ($(this).attr("id") == id && id.length != 1) {
                        $(this).css("fill", newcolor[1]);
                    };
                };
            });

        }
    });
 }

 function colour_lines() {
     $(".path").click(function () {
         newcolor = $(".pick-color.selected").attr('id');
         newcolor = newcolor.split('-');

         $(this).css("stroke", newcolor[1]);
         $(this).css("stroke-width", '4');
     });
     $(".path2").click(function () {
         newcolor = $(".pick-color.selected").attr('id');
         newcolor = newcolor.split('-');

         $(this).css("stroke", newcolor[1]);
         $(this).css("stroke-width", '4');
     });
 };
 function predefined_colours(defs, colours) {
     resetColors(defs);
     $('#svgCanvas').find(".bgfield").each(function (index) {
         if (colours['proteins'].indexOf($(this).attr("id")) > -1) {
             $(this).css("fill", colours['colours'][0]);
         };
     });
 };
 function resetColors(defs) {
     $('#svgCanvas').find(".path").each(function (index) {
         $(this).css("stroke", '');
         $(this).css("stroke-width", '');
     });
     $('#svgCanvas').find(".path2").each(function (index) {
         $(this).css("stroke", '');
         $(this).css("stroke-width", '');
     });
     $('#svgCanvas').find(".bgfield").each(function (index) {
         $(this).css("fill", '');
         $(this).css("stroke", '');
     });
     $("[class^=chart]").each(function (index) {
         if ($(this).attr('id')) {
             $(this).css("fill", defs[$(this).attr('id')][0]);
         };
     });


 };

 function toggleLegend() {
     var chart0 = $('#svgCanvas').find(".chart0")
     var chart1 = $('#svgCanvas').find(".chart1")
     var chart2 = $('#svgCanvas').find(".chart2")
     if ($(chart0).css("visibility") == 'hidden' && $(chart1).css("visibility") == 'hidden' && $(chart2).css("visibility") == 'hidden') {
         $(chart0).css("visibility", 'visible');
         $(chart1).css("visibility", 'visible');
         $(chart2).css("visibility", 'visible');
     } else {
         $(chart0).css("visibility", 'hidden');
         $(chart1).css("visibility", 'hidden');
         $(chart2).css("visibility", 'hidden');
     };
 };


 function toggleRings(ring) {
     $('#svgCanvas').find("." + ring).each(function (index) {
         if ($(this).css("visibility") == 'hidden') {
             $(this).css("visibility", 'visible');
         } else {
             $(this).css("visibility", 'hidden');
         };
     });
 };

 function mergeSVG() {
     var SVG = $('#svgCanvas').find('svg')[0];
     h = parseInt($(SVG).attr('height'));
     w = parseInt($(SVG).attr('width'));
     var legend = $('#legend').find('svg')[0];
     h2 = parseInt($(legend).attr('height'));
     w2 = parseInt($(legend).attr('width'));
     new_w = (w-w2)/2
     SVG.setAttribute('height', (h + h2));
     for (i = 0; i < legend.children.length; i++) {
         legend.children[i].setAttribute('transform', 'translate ('+new_w.toString()+' ' + h.toString()+')');
         $(SVG).append(legend.children[i]);
     }

 };
$( document ).ready(function (){

    $(".pick-color").click(function () {
        $(".pick-color").css('borderWidth', '2px');
        $(".pick-color").css('height', '20px');
        $(".pick-color").removeClass('selected');
        $(this).css('borderWidth', '3px');
        $(this).css('height', '22px');
        $(this).addClass('selected');

    });



    $("#tree_svg_link").click(function () {
        svgAsDataUri(document.getElementById("svgtree"), {}, function (uri) {
            $("#tree_svg_link").attr('href', uri);
        });
    });


});



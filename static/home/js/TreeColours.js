
 function colours() {

    console.log($(".bgfield").length);
    console.log('were done loading!');
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
        console.log('klik');
        parentid = $(this).closest('svg').attr('id');
        newcolor = $(".pick-color.selected").attr('id');
        newcolor = newcolor.split('-');

        $(this).css("fill", newcolor[1]);


        //$(this).next().css("fill", newcolor[2]);
    });
    $(".chart").click(function () {
        console.log('klik');
        id = $(this).attr('id');
         
        newcolor = $(".pick-color.selected").attr('id');
        newcolor = newcolor.split('-');
        $('#svgCanvas').find(".chart").each(function (index) {

            if ($(this).attr("id") == id && id.length != 1) {
                $(this).css("fill", newcolor[1]);
            };
            // $(this).next().css("fill", 'black');
        });
        $('#legend').find(".chart").each(function (index) {

            if ($(this).attr("id") == id && id.length != 1) {
                $(this).css("fill", newcolor[1]);
            };
            // $(this).next().css("fill", 'black');
        });
    });

 
 }
 function resetColors(defs) {
     console.log($(".bgfield").length);

     $('#svgCanvas').find(".bgfield").each(function (index) {
         $(this).css("fill", '#EEE');
     });
     $('#svgCanvas').find(".chart").each(function (index) {
         $(this).css("fill", defs[$(this).attr('id')][0]);
     });
     $('#legend').find(".chart").each(function (index) {
         $(this).css("fill", defs[$(this).attr('id')][0]);
     });
 };

 function toggleLegend() {
     $('#svgCanvas').find(".chart").each(function (index) {
         console.log($(this).css('visibility'));
         if ($(this).css("visibility") == 'hidden') {
             $(this).css("visibility", 'visible');
         } else {
             $(this).css("visibility", 'hidden');
         };
     });
    $('#svgCanvas').find(".legend").each(function (index) {
        console.log($(this).css('visibility'));
        if ($(this).css("visibility")== 'hidden') {
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
     // new_h = h+
     var legend = $('#legend').find('svg')[0];
    // legend.setAttribute('transform', 'translate (' + h + 20 + ' 0)');
     h2 = parseInt($(legend).attr('height'));
     w2 = parseInt($(legend).attr('width'));
     new_w = (w-w2)/2
     SVG.setAttribute('height', (h + h2));
     for (i = 0; i < legend.children.length; i++) {
         legend.children[i].setAttribute('transform', 'translate ('+new_w.toString()+' ' + h.toString()+')');
         console.log(legend.children[i]);
         $(SVG).append(legend.children[i]);
     }
     
     console.log(legend.children);
     console.log($(SVG).attr('height'));

     //console.log(svg1.css('height'));
     //console.log(svg1);
     //console.log(newSVG);
 };
$( document ).ready(function (){

    $(".pick-color").click(function () {
        console.log($(this).attr('id'));
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



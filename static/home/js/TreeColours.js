
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
        $(this).css("stroke", newcolor[1]);
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
     resetColors();
     $('#svgCanvas').find(".bgfield").each(function (index) {
         if (defs.indexOf($(this).attr("id")) > -1) {
             $(this).css("fill", colours);
         };
     });
 };
 
 function resetColors() {
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
             $(this).css("fill",'');
             $(this).css("stroke", '');
         };
     });


 };

function toggleAll() {
     var style = ''
     $('#svgCanvas').find("[class^=chart]").each(function (index) {
         if ($(this).css("visibility") == 'visible'){
             style = 'hidden'};
     });
     $('#svgCanvas').find("[class^=chart]").each(function (index) {
         $(this).css("visibility", style);
    });
 };
function toggleLegend() {
    var style = 'visible'
    $('#svgCanvas').find(".legend").each(function (index) {
        if ($(this).css("visibility") == 'visible' || $(this).css("visibility") == '') {
            style = 'hidden'
        };
    });
    $('#svgCanvas').find(".legend").each(function (index) {
        $(this).css("visibility", style);
    });
};
 //function toggleRings(ring) {
 //    $('#svgCanvas').find("." + ring).each(function (index) {
 //        if ($(this).css("visibility") == 'hidden') {
 //            $(this).css("visibility", 'visible');
 //        } else {
 //            $(this).css("visibility", 'hidden');
 //        };
 //    });
 //};
 
  function SelectSubmenu(name) {
     $('.button_container').find('.btn-group').each(function (index) {
         if ($(this).attr('id') == name || $(this).attr('id') == 'types' || $(this).attr('id') == 'choice') {
             $(this).css("display", '');
         } else {
             $(this).css("display", 'none');
         };
     });
 };
 function SelectButtons(name) {
     $('.button_container').find('.btn-default').each(function (index) {
         if ($(this).attr('id') == name) {
             $(this).css("display", '');
         };
     });
 };
  function DeselectButtons(name) {
     $('.button_container').find('.btn-default').each(function (index) {
         if ($(this).attr('id') == name) {
             $(this).css("display", 'none');
         };
     });
 };
 function GetRings() {
     var args = [];
     var values = [];
     $('.button_container').find('.btn-default').each(function (index) {
         if ($(this).css("display") != 'none'){
             args.push($(this).attr('id'));
             values.push('True');
         } else if ($(this).css("display") == 'none'){
             args.push($(this).attr('id'));
             values.push('False');
         };
     });
    $("#svgCanvas").empty(); 
    $.ajax({
    'url': '/phylogenetic_trees/showrings',
    'data': {
        arg: args,
        value: values
    },
    'type': 'GET',
    'success': function(data) {
           $("#main").html(data);
       }
     });
    $.ajax({
    'url': '/phylogenetic_trees/get_buttons',
    'type': 'GET',
    'success': function(data) {
           $("#ring_buttons").html(data);
       }
     });
       };



 function mergeSVG() {
     var SVG = $('#svgCanvas').find('svg')[0];
     h = parseInt($(SVG).attr('height'));
     w = parseInt($(SVG).attr('width'));
     var legend = $('#legend').find('svg')[0];
     h2 = parseInt($(legend).attr('height'));
     w2 = parseInt($(legend).attr('width'));
     leg_w = (w-w2)/2
     
     SVG.setAttribute('height', (h + h2));
     if (w2 > w) {
         SVG.setAttribute('width', (w2));
         leg_w = 0 
         svg_w = Math.abs(w-w2)/2
     } else {
         leg_w = Math.abs(w-w2)/2 
         svg_w = 0
     };
     for (i = 0; i < legend.children.length; i++) {
         legend.children[i].setAttribute('transform', 'translate ('+leg_w.toString()+' ' + h.toString()+')');
         $(SVG).append(legend.children[i]);
     };


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



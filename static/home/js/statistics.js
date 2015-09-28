 function mergeSVG(div) {
     var SVG = $('#'+div).find('svg')[0];
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
window.onload = function(){
    $(".chart_type").click(function () {
        Clear_all()
        $(this).css("fill", '#000000');
        point = $('#' + $(this).attr('id')).find('svg')
        $(point).css("visibility", 'hidden');
        $('#'+$(this).attr('id') + '.chart_container').css("display", '');
    });
    $(".Class_phylo").click(function () {
        Clear_phylo()
        var id = $(this).attr('id');
        if ($('#' + id + '.container').find("svg").length == 0) {
            console.log('Drawing...');
            draw(id.substr(id.length - 1));
        };

        $(this).css("fill", '#000000');
        //point = $('#' + $(this).attr('id')).find('svg')
        //$(point).css("visibility", 'hidden');
        $('#' + $(this).attr('id') + '.container').css("display", '');
    });

    function Clear_all() {
        $('#charts').find(".chart_type").each(function (index) {
            $(this).css("fill", '');
        });
        $('#charts').find(".chart_container").each(function (index) {
            $(this).css("display", 'none');
        });
    };
    function Clear_phylo(){
        $('#phylos').find(".Class_phylo").each(function (index) {
            $(this).css("fill", '');
            $('#'+$(this).attr('id')+".container").css("display", 'none');
        });
    };
    function draw(number) {
        $.ajax({
            'url': '/structure/showtrees',
            'data': {
                number: number
            },
            'type': 'get',
            'success': function (data) {
                $("#phylos_div").html(data);
            }
        });
    };
    
 
    $(document).ready(function () {
        $('#phylo_1.Class_phylo').css("fill", '#000000');
        $('#phylo_1.container').css("display", '');
        $('#unique.chart_type').css("fill", '#000000');
        $('#unique.chart_container').css("display", '');
        draw('1');

    });
        
};

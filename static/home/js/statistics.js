window.onload = function(){
    $(".Class").click(function () {
        Clear_phylos()
        $(this).css("fill", '#000000');
        point = $('#' + $(this).attr('id')).find('svg')
        $(point).css("visibility", 'hidden');
        $('#'+$(this).attr('id') + '.container').css("display", '');
    });
    $(".chart_type").click(function () {
        Clear_charts()
        $(this).css("fill", '#000000');
        point = $('#' + $(this).attr('id')).find('svg')
        $(point).css("visibility", 'hidden');
        $('#' + $(this).attr('id') + '.chart_container').css("display", '');
    });

    function Clear_phylos(){
        $('#phylos').find(".Class").each(function (index) {
            $(this).css("fill", '');
        });
        $('#phylos').find(".container").each(function (index) {
            $(this).css("display", 'none');
        });
    };
    function Clear_charts(){
        $('#charts').find(".chart_type").each(function (index) {
            $(this).css("fill", '');
        });
        $('#charts').find(".chart_container").each(function (index) {
            $(this).css("display", 'none');
        });
    };
    $(document).ready(function () {
        $('#A.Class').css("fill", '#000000');
        $('#A.container').css("display", '');
        $('#unique.chart_type').css("fill", '#000000');
        $('#unique.chart_container').css("display", '');

    });
        
};

window.onload = function(){
    $(".Class").click(function () {
        Clear_all()
        $(this).css("fill", '#000000');
        point = $('#' + $(this).attr('id')).find('svg')
        $(point).css("visibility", 'hidden');
        $('#'+$(this).attr('id') + '.container').css("display", '');
    });
    $(".Class_phylo").click(function () {
        Clear_phylo()
        $(this).css("fill", '#000000');
        point = $('#' + $(this).attr('id')).find('svg')
        $(point).css("visibility", 'hidden');
        $('#' + $(this).attr('id') + '.container').css("display", '');
    });

    function Clear_all() {
        $('#charts').find(".Class").each(function (index) {
            $(this).css("fill", '');
        });
        $('#charts').find(".container").each(function (index) {
            $(this).css("display", 'none');
        });
    };
    function Clear_phylo(){
        $('#phylos').find(".Class_phylo").each(function (index) {
            $(this).css("fill", '');
        });
        $('#phylos').find(".container").each(function (index) {
            $(this).css("display", 'none');
        });
    };
    $(document).ready(function () {
        $('#phylo_A.Class_phylo').css("fill", '#000000');
        $('#phylo_A.container').css("display", '');
        $('#unique.Class').css("fill", '#000000');
        $('#unique.container').css("display", '');

    });
        
};

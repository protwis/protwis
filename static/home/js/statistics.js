window.onload = function(){
    $(".Class").click(function () {
        Clear_all()
        $(this).css("fill", '#000000');
        point = $('#' + $(this).attr('id')).find('svg')
        $(point).css("visibility", 'hidden');
        $('#'+$(this).attr('id') + '.container').css("display", '');
    });

    function Clear_all(){
        $('#phylos').find(".Class").each(function (index) {
            $(this).css("fill", '');
        });
        $('#phylos').find(".container").each(function (index) {
            $(this).css("display", 'none');
        });
    };
    $(document).ready(function () {
        $('#A.Class').css("fill", '#000000');
        $('#A.container').css("display", '');
    });
        
};

$(function () {
    $('[data-toggle="tooltip"]').tooltip()
})

$(function(){
    $('.ali-scroll-div').scroll(function(){
        $('.ali-main-div')
            .scrollLeft($('.ali-scroll-div').scrollLeft());
    });
    $('.ali-main-div').scroll(function(){
        $('.ali-scroll-div')
            .scrollLeft($('.ali-main-div').scrollLeft());
    });
});

$(window).load(function () {
    $('.internal-scroll-div').css('width', $('.dynamic-div').outerWidth() );
});
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

function ScoreBreakdown (protein_conformation, cutoff) {
    $.ajax({
        'url': '/alignment/score_breakdown',
        'data': {
            protein_conformation: protein_conformation,
            cutoff: cutoff,
        },
        'type': 'GET',
        'async': false,
        'success': function(data) {
            $("#pconf-" + protein_conformation).html(data);
        },
    });
}

$('#cutoff-apply').click( function() {
    var cutoff = parseInt($('#cutoff-val').val());
    var row = [];
    $('#first-row td').each(function(){
        row.push($(this).attr('id'));
    });
    //console.info(row.length);
    for (i=0; i < row.length; i++) {
        if (row[i] == 'anchor'){
            continue;
        }
        var cell = Math.abs(row[i].replace('cutoff-', ''));
        j = i + 1;
        if (cell <= cutoff){
            $('td:nth-child('+j+')').hide();
        }
        else {
            $('td:nth-child('+j+')').show();
        }
    };
});
$(window).on("load", function () {
    $('.internal-scroll-div').css('width', $('.dynamic-div').outerWidth() );

    var row = [];
    $('#first-row td').each(function(){
        row.push($(this).attr('id'));
    });
    // console.info(row.length);
    for (i=0; i < row.length; i++) {
        if (row[i] == 'anchor'){
            continue;
        }
        var qq = Math.abs(row[i].replace('cutoff-', ''));
        // console.info(qq);
        if (qq < 70){
            j = i + 1;
            $('td:nth-child('+j+')').hide();
        }
    };


    // console.info(row);

    // var w = 0;
    // $('#test tr td:first').each(function(){
    //     if($(this).width() > w) {
    //         w = $(this).width();
    //     }
    // });
    // console.info(w);
});
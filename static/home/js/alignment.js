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
};

function ApplyCutoff (cutoff) {
    var row = [];
    $('#first-row td').each(function(){
        row.push($(this).attr('id'));
    });
    for (i=0; i < row.length; i++) {
        if (row[i] == 'anchor'){
            continue;
        }
        var cell = Math.abs(row[i].replace('cutoff-', ''));
        j = i + 1;
        if (cell < cutoff){
            console.info('Got it!');
            $('#signature-table td:not(".ali-td-segment-title, .ali-td-header-row, .ali-td-anchor"):nth-child('+j+')').hide();//.css('background-color', 'gray');
        }
        else {
            $('#signature-table td:not(".ali-td-segment-title, .ali-td-header-row, .ali-td-anchor"):nth-child('+j+')').show();
        }
    };
    $('#signature-table #segments td:not("#anchor")').each(function(){
        var segment_name = $(this).attr('id').replace('segment-', '');
        var colspan = $('#gns td#gn-'+segment_name+':visible').length;
        $(this).attr('colspan', colspan);
    });
};

$('#cutoff-apply').click( function() {
    var cutoff = parseInt($('#cutoff-val').val());
    ApplyCutoff(cutoff);
});

$(window).on("load", function () {
    $('.internal-scroll-div').css('width', $('.dynamic-div').outerWidth() );
    ApplyCutoff(40);
});
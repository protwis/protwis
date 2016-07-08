$(function () {
    $('[data-toggle="tooltip"]').tooltip({ container: 'body', html: true})
})

$(document).ready(function() {

    $('#seq_c').hide();

    $('input[type=radio][name=sequence]').change(function() {
        if (this.value == 'wt') {
           $('#seq_wt').slideDown();
           $('#seq_c').slideUp();
        }
        else if (this.value == 'construct') {
           $('#seq_c').slideDown();
           $('#seq_wt').slideUp();
        }
    });

    $('#schematic_seq_c').hide();

    $('input[type=radio][name=schematic]').change(function() {
        if (this.value == 'wt') {
           $('#schematic_seq_wt').slideDown();
           $('#schematic_seq_c').slideUp();
        }
        else if (this.value == 'construct') {
           $('#schematic_seq_c').slideDown();
           $('#schematic_seq_wt').slideUp();
        }
    });
});
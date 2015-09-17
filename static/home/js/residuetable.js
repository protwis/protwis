presetColors = {'D': ['#E60A0A', '#FDFF7B'],'E': ['#E60A0A', '#FDFF7B'],
                'K': ['#145AFF', '#FDFF7B'],'R': ['#145AFF', '#FDFF7B'],
                'S': ['#A70CC6', '#FDFF7B'],'T': ['#A70CC6', '#FDFF7B'],
                'N': ['#A70CC6', '#FDFF7B'],'Q': ['#A70CC6', '#FDFF7B'],
                'V': ['#E6E600', '#000000'],'L': ['#E6E600', '#000000'],
                'I': ['#E6E600', '#000000'],'A': ['#E6E600', '#000000'],
                'M': ['#E6E600', '#000000'],'F': ['#18FF0B', '#000000'],
                'Y': ['#18FF0B', '#000000'],'W': ['#0BCF00', '#000000'],
                'H': ['#0093DD', '#000000'],'P': ['#CC0099', '#FDFF7B'],
                'C': ['#B2B548', '#000000'],'G': ['#FF00F2', '#000000'],
                '-': ['#FFFFFF', '#000000'],
                '+': ['#FFFFFF', '#000000']        
                };

function table_applyPresentColors() {

    //console.log( $('#'+target));

    $('.residue').each(function( index ){
          console.log($(this).text());

          aa =  $(this).text().trim().charAt(0);
          console.log(aa);
          $(this).css("background-color", presetColors[aa][0]);
          $(this).css("color", presetColors[aa][1]);
        });

};

function table_resetColors() {


    $('.residue').each(function( index ){
          $(this).css("background-color", 'white');
          $(this).css("color", 'black');
        });

}


function table_ajaxMutants() {

    $('.protein').each(function( index ){
        console.log($(this).find('span').text());
        console.log($(this).index());
        var protein = $(this).find('span').text();
        var protein_index = $(this).index()+1;
        $.getJSON( '/mutations/ajax/'+protein+'/', function( data ) {
                count = 0;
              $.each( data, function( key, val ) {
                count = count + 1;
                 //console.log('#'+protein_index+"#"+key);

                 max = String(Math.max.apply(null, val));
                 min = String(Math.min.apply(null, val));
                 extra = "\n" + String(val.length) + " mutations | "+ max +" maxFold | "+ min +" minFold";

                 //original_title = $('#'+plotid).find("#"+key).attr('original_title')
                 //console.log($('#P'+protein_index+"R"+key).text());
                 $('#P'+protein_index+"R"+key).css("background-color", "#E60A0A");
                 $('#P'+protein_index+"R"+key).css("color", "#FDFF7B");
                 //$('#'+plotid).find("#"+key).attr('title',original_title+extra);
                 //$('#'+plotid).find("#"+key+"t").attr('title',original_title+extra);


              });
              console.log("Mutants for "+protein+ " : "+count);

        });


    });


}

function table_ajaxInteractions() {

    $.getJSON( '/interaction/ajax/'+protein+'/', function( data ) {
      $.each( data, function( key, val ) {

        var flags = [], falgsAA = [], output = [], outputAA = [], l = val.length, i;
        for( i=0; i<l; i++) {
            if( flags[val[i][1]]) continue;
            flags[val[i][1]] = true;
            output.push(val[i][1]);
        }
        for( i=0; i<l; i++) {
            if( flags[val[i][0]]) continue;
            flags[val[i][0]] = true;
            outputAA.push(val[i][0]);
        }
         
         extra = "\n" + String(val.length) + " interactions | Type: "+ output +" | Residue in crystal:"+ outputAA;


         $('#'+plotid).find("#"+key).css("fill", "#E60A0A");
         $('#'+plotid).find("#"+key).next().css("fill", "#FDFF7B");

         original_title = $('#'+plotid).find("#"+key).attr('original_title')


         $('#'+plotid).find("#"+key).attr('title',original_title+extra);
         $('#'+plotid).find("#"+key+"t").attr('title',original_title+extra);


      });
    $("circle").tooltip('fixTitle');
    $("text").tooltip('fixTitle');

    });
}



$(".pick-color").click(function() {
    plottype = $(this).attr('class').split(' ')[1];
    
    console.log($(this).attr('id'));
    $(".pick-color."+plottype).css('borderWidth','2px');
    $(".pick-color."+plottype).css('height','20px');
    $(".pick-color."+plottype).removeClass('selected');
    $(this).css('borderWidth','3px');
    $(this).css('height','22px');
    $(this).addClass('selected');
    
});
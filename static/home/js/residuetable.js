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

function ajaxMutants(plotid,protein) {

                        $.getJSON( '/mutations/ajax/'+protein+'/', function( data ) {
                          $.each( data, function( key, val ) {

                             var ligands = [], bigincreases=0, increases = 0, bigdecreases=0, decreases = 0, unchanged=0, unspecified = 0;
                             
                             
                             $.each( val, function( key, v ) {
                               if( !(ligands[v[1]]) ) ligands[v[1]] = [];
                              ligands[v[1]].push(v[0])
                              if (v[0]>10) {
                                bigincreases ++; //mix-up increase is decrease.
                              } else if (v[0]>5) {
                                increases ++;
                              } else if (v[0]>0) {
                                unchanged ++;
                              }  else if (v[0]<-10) {
                                bigdecreases ++;
                              } else if (v[0]<-5) {
                                decreases ++;
                              } else if (v[0]<0) {
                                unchanged ++;
                              } else if (v[2]=='No effect on') {
                                unchanged ++;
                              } else if (v[2]=='No effect') {
                                unchanged ++;
                              } else if (v[2]=='Abolish') {
                                bigincreases ++;
                              } else if (v[2]=='Abolished effect') {
                                bigincreases ++;
                              } else if (v[2]=='Gain of') {
                                bigdecreases ++;
                              } else if (v[2]=='Increase') {
                                increases ++;
                              } else if (v[2]=='Decrease') {
                                decreases ++;
                              } else {
                                unspecified ++;
                              }
                             });
                             
                             extra = "\n" + String(val[0].length) + " mutations: " +
                              (decreases+bigdecreases) +" increases | " +
                             (increases+bigincreases) +" decreases  |  " +
                              (unchanged) +" Unchanged | " +
                              unspecified + " Unspecified";

                            counts = [(increases+bigincreases),(decreases+bigdecreases),(unchanged)];
                            winner = counts.indexOf(Math.max.apply(window,counts));
                            winner2 = Math.max.apply(window,counts);
                            color = "#D9D7CE";
                            color_letter = "#000";
                            if (winner==0 && winner2) {
                              if (increases>bigincreases) {
                                color = "#F05960";
                                color_letter = "#FDFF7B";
                              } else {
                                color = "#CC434A";
                                color_letter = "#FDFF7B";
                              }
                            } else if (winner==1) {
                              if (decreases>bigdecreases) {
                                color = "#87E88F";
                              } else {
                                color = "#66B36C";
                              }
                            } else if (winner==2) {
                              color = "#F7DA00";
                              color_letter = "#000";
                            }

                            console.log(counts + " " + counts.indexOf(Math.max.apply(window,counts)));


                             original_title = $('#'+plotid).find("#"+key).attr('original_title')

                             $('#'+plotid).find("#"+key).css("fill", color);
                             $('#'+plotid).find("#"+key).next().css("fill",color_letter );
                             $('#'+plotid).find("#"+key).attr('title',original_title+extra);
                             $('#'+plotid).find("#"+key+"t").attr('title',original_title+extra);


                          });
                        $("circle").tooltip('fixTitle');
                        $("text").tooltip('fixTitle');
    
                        });
                    }

function table_ajaxMutants() {


    $('.protein').each(function( index ){
        var protein = $(this).find('span').attr('id');
        var protein_index = $(this).index()+1;
        console.log(protein);
        $.getJSON( '/mutations/ajax/'+protein+'/', function( data ) {
                count = 0;
              $.each( data, function( key, val ) {
                count = count + 1;
                 // console.log('#'+protein_index+"#"+key+" val:"+val);

                            var ligands = [], bigincreases=0, increases = 0, bigdecreases=0, decreases = 0, unchanged=0, unspecified = 0;
                             
                             
                             $.each( val, function( key, v ) {
                               if( !(ligands[v[1]]) ) ligands[v[1]] = [];
                              ligands[v[1]].push(v[0])
                              if (v[0]>10) {
                                bigincreases ++; //mix-up increase is decrease.
                              } else if (v[0]>5) {
                                increases ++;
                              } else if (v[0]>0) {
                                unchanged ++;
                              }  else if (v[0]<-10) {
                                bigdecreases ++;
                              } else if (v[0]<-5) {
                                decreases ++;
                              } else if (v[0]<0) {
                                unchanged ++;
                              } else if (v[2]=='No effect on') {
                                unchanged ++;
                              } else if (v[2]=='No effect') {
                                unchanged ++;
                              } else if (v[2]=='Abolish') {
                                bigincreases ++;
                              } else if (v[2]=='Abolished effect') {
                                bigincreases ++;
                              } else if (v[2]=='Gain of') {
                                bigdecreases ++;
                              } else if (v[2]=='Increase') {
                                increases ++;
                              } else if (v[2]=='Decrease') {
                                decreases ++;
                              } else {
                                unspecified ++;
                              }
                             });
                             
                             extra = "\n" + String(val[0].length) + " mutations: " +
                              (decreases+bigdecreases) +" increases | " +
                             (increases+bigincreases) +" decreases  |  " +
                              (unchanged) +" Unchanged | " +
                              unspecified + " Unspecified";

                            counts = [(increases+bigincreases),(decreases+bigdecreases),(unchanged)];
                            winner = counts.indexOf(Math.max.apply(window,counts));
                            winner2 = Math.max.apply(window,counts);
                            color = "#D9D7CE";
                            color_letter = "#000";
                            if (winner==0 && winner2) {
                              if (increases>bigincreases) {
                                color = "#FF7373";
                                color_letter = "#FFF";
                              } else {
                                color = "#FA1111";
                                color_letter = "#FDFF7B";
                              }
                            } else if (winner==1) {
                              if (decreases>bigdecreases) {
                                color = "#87E88F";
                              } else {
                                color = "#66B36C";
                              }
                            } else if (winner==2) {
                              color = "#F7DA00";
                              color_letter = "#000";
                            }

                 //original_title = $('#'+plotid).find("#"+key).attr('original_title')
                 //console.log($('#P'+protein_index+"R"+key).text());
                 $('#P'+protein_index+"R"+key).css("background-color", color);
                 $('#P'+protein_index+"R"+key).css("color", color_letter);
                 $('#P'+protein_index+"R"+key).attr('title', extra);
                 //$('#'+plotid).find("#"+key).attr('title',original_title+extra);
                 //$('#'+plotid).find("#"+key+"t").attr('title',original_title+extra);


              });
              // console.log("Mutants for "+protein+ " : "+count);

        });


    });


}

function table_ajaxInteractions() {
  $('.protein').each(function( index ){
      var protein = $(this).find('span').attr('id');
      var protein_index = $(this).index()+1;
      console.log(protein);
    $.getJSON( '/interaction/ajax/'+protein+'/', function( data ) {
      $.each( data, function( key, val ) {
        console.log('#P'+protein_index+"R"+key);

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


          $('#P'+protein_index+"R"+key).css("background-color", "#E60A0A");
          $('#P'+protein_index+"R"+key).css("color", "#FDFF7B");

         original_title =  $('#P'+protein_index+"R"+key).attr('original_title')


          $('#P'+protein_index+"R"+key).attr('title',extra);


      });
    $("circle").tooltip('fixTitle');
    $("text").tooltip('fixTitle');

    });
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
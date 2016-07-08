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

                    var translateOffset = 0;
                    function showToolTip(x, y, str,rid,plotid) {
                        var tipElement = document.getElementById('tool-tip-'+plotid);

                        var originalCircle = $('#'+plotid).find("#"+rid);

                        //console.log(originalCircle.attr('extra'));
                        //originalCircle.attr('extra');

                        var rect = tipElement.childNodes[1];
                        var text = tipElement.childNodes[3];

                        while (text.lastChild) {
                           text.removeChild(text.lastChild);
                        }
                        
                        // var NS = "http://www.w3.org/2000/svg";


                        // //text.textContent =  str;
                        //     var text_tspan = document.createElementNS(NS, "tspan");


                        //     rect.setAttribute('height', 25);
                        //     rect.setAttribute('y', -40);
                        //     text_tspan.textContent = String(str);
                        //     text.appendChild(text_tspan);

                        if (originalCircle.attr('extra')) { //Only display if mutateddata flag is on and there is info

                            text.innerHTML =  "<tspan  x=\"0\" y=\"-33\">" + str + "</tspan>";
                            //text.innerHTML =  text.innerHTML + "<tspan  x=\"0\" y=\"-20\">" + residueColor.mutatedCalc[rid][0] + " mutations" + " | " + residueColor.mutatedCalc[rid][1] + " maxFold"+ " | " + residueColor.mutatedCalc[rid][2] + " minFold </tspan>";
                            text.innerHTML =  text.innerHTML + "<tspan  x=\"0\" y=\"-20\">" + originalCircle.attr('extra') + " </tspan>";
                            rect.setAttribute('height', 35);
                            rect.setAttribute('y', -50);
                        } else {
                            text.innerHTML =  "<tspan  x=\"0\" y=\"-23\">" + str + "</tspan>";
                            rect.setAttribute('height', 25);
                            rect.setAttribute('y', -40);
                        }
                        
                        var bbox = text.getBBox();
                        rect.setAttribute('width', bbox.width + 8);
                        rect.setAttribute('x', -bbox.width/2 - 4);

                        var transX = (x <= (bbox.width + 8) / 2) ? (bbox.width + 8) / 2 : x;
                        tipElement.setAttribute('transform', 'translate(' + transX + ',' + (y + translateOffset) + ')');
                        tipElement.setAttribute('visibility', 'visible');
                    }

                    function hideToolTip(plotid) {
                        var tipElement = document.getElementById('tool-tip-'+plotid);
                        tipElement.setAttribute('visibility', 'hidden');
                    }

                    function toggleLoop(id,type, skipmaxmin) {
                        $(id+".long").toggle();
                        $(id+".short").toggle();

                        // $(id+".long").each(function () {
                        //     curr = $(this).css("display");
                        //     if (curr == 'none') $(this).removeAttr("display");
                        //     if (!curr) $(this).css("display", "none");
                        // });

                        // $(id+".short").each(function () {
                        //     curr = $(this).css("display");
                        //     if (curr == 'none') $(this).removeAttr("display");
                        //     if (!curr) $(this).css("display", "none");
                        // });

                        if (skipmaxmin!=1) maxmin();

                        // $(id+".long").toggle();
                        // $(id+".short").toggle();

                        // $(id+".long").fadeToggle();
                        // $(id+".short").fadeToggle();
                    }

                    function applyPresentColors(target) {

                        //console.log( $('#'+target));

                        $('#'+target).find("circle").each(function( index ){
                              //console.log( index + ": " + $( this ).text() );
                              aa =  $(this).next().text().trim();
                              //console.log( index + ": " + aa );
                              $(this).css("fill", presetColors[aa][0]);
                              $(this).next().css("fill", presetColors[aa][1]);
                            });

                    };

                    function resetColors(target) {


                        $('#'+target).find("circle").each(function( index ){
                              //console.log( index + ": " + $( this ).text() );
                              aa =  $(this).next().text();
                              //console.log( index + ": " + aa );
                              $(this).css("fill", 'white');
                              $(this).next().css("fill", 'black');
                            });

                    }

                    function maxmin() {
                        margin = 50;
                        svgmax = 0;
                        svgmin = 0;
                        count = 0;
                        classmax = '';
                        classmin = '';
                        counter = 0;
                        $('#snake').children('.rtext').each(function () {
                            counter += 1;
                            y = parseInt($(this).attr( "y" ));
                            classtext = $(this).attr( "class" );
                            test = $(this).attr("original_title");
                            test2 = $(this).css("display");
                            if ($(this).css("display")!='none') {
                                count = count +1;
                                if (y<svgmin) {
                                    svgmin = y; 
                                    classmin = classtext;
                                    }
                                if (y>svgmax) {
                                    classmax = classtext;
                                    svgmax= y; 
                                 }

                            }   
                        });
                        console.log('max '+svgmax+' '+classmax+' min'+svgmin+' '+classmin+' count'+count);
                        
                        var svg = $('#snake').closest('svg')[0];

                        check = document.getElementById("snakeplot").getAttribute('viewBox');
                        console.log('check',check)
                        if (typeof check !== typeof undefined && check !== false && check !== null ) {
                          oldheight = check.split(" ")[3];
                          width = check.split(" ")[2];
                        } else {
                          console.log('not found it');
                          oldheight = $(svg).attr('height');
                          width = $("#snakeplot").attr('width');
                        }

                        newheight = (svgmax-svgmin+margin*2);
                        svg.setAttribute('height', (svgmax-svgmin+margin*2));

                        console.log('New height:'+ (svgmax-svgmin) +' old height:'+oldheight);
                        console.log("Prev attr"+$('#snake').attr("transform"));

                        $('#snake').attr("transform", "translate(0," + (-svgmin+margin) + ")");
                        
                        // $('#snakeplot')[0].attr("viewBox", "0 0 " + width + " " + newheight);
                        document.getElementById("snakeplot").setAttribute("viewBox", "0 0 " + width + " " + newheight);


                        svg.setAttribute('height', "100%");
                        svg.setAttribute('width', "100%");

                        console.log("New attr"+$('#snake').attr("transform"));
                    }

                    $( document ).ready(function() {    
                        // var elements = document.getElementsByClassName('long')

                        // for (var i = 0; i < elements.length; i++){
                        //     elements[i].style.display = 'none';
                        // }
                        $(".long").hide();


                        $('rect').each(function(){
                            
                            rectclass = $(this).attr('class');
                            if (rectclass) {
                                if (rectclass.indexOf("CL") >= 0 && rectclass.indexOf("long") >= 0) {

                                    numResidues = ($('.'+rectclass.replace(/ /g,".")).length-3)/2

                                    console.log('class:'+rectclass+' count'+numResidues);

                                    if (numResidues<10) {
                                        toggleLoop('.'+rectclass.split(' ')[0],'',1);
                                    }
                                }
                            }
                        });

                        maxmin();

                        $("text").tooltip({
                            'container': 'body',
                            'placement': 'top',
                            'animation': false,
                            'html' : true
                        });


                        $("circle").tooltip({
                            'container': 'body',
                            'placement': 'top',
                            'animation': false,
                            'html' : true
                        });

                        $("circle").hover(function(){
                            $('.tooltip').css('top',parseInt($('.tooltip').css('top')) + 2.8 + 'px')
                        });

                    });

                    $(".rtext").click(function() {
                        parentid = $(this).closest('svg').attr('id');
                        newcolor = $(".pick-color."+parentid+".selected").attr('id');
                        if (newcolor) {
                          newcolor = newcolor.split('-');
                        } else {
                          custom =  $("#custom_color_"+parentid).val();
                          custom_text = getContrast50(custom);
                          newcolor = ["pick",custom,custom_text];
                        }
                      console.log(newcolor);

                      $(this).css("fill", newcolor[2]);
                      $(this).prev().css("fill", newcolor[1]);
                    });
                    
                    $(".rcircle").click(function() {
                        parentid = $(this).closest('svg').attr('id');
                        newcolor = $(".pick-color."+parentid+".selected").attr('id');

                        if (newcolor) {
                          newcolor = newcolor.split('-');
                        } else {
                          custom =  $("#custom_color_"+parentid).val();
                          custom_text = getContrast50(custom);
                          newcolor = ["pick",custom,custom_text];
                        }
                      console.log(newcolor);
                      $(this).css("fill", newcolor[1]);
                      $(this).next().css("fill", newcolor[2]);
                   });



                    $("#snake_svg_link").click(function() {
                        svgAsDataUri(document.getElementById("snakeplot"),{}, function(uri) {
                            $("#snake_svg_link").attr('href',uri);
                        });
                    });
                    $("#helix_svg_link").click(function() {
                        svgAsDataUri(document.getElementById("helixbox"),{}, function(uri) {
                            $("#helix_svg_link").attr('href',uri);
                        });
                    });

                    function ajaxMutants(plotid,protein) {

                      resetColors(plotid);

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

                    function ajaxMutantsPos(plotid) {

                      resetColors(plotid);

                        var pos = jQuery.parseJSON(mutant_json);

                        $.each( pos, function( key, val ) {
                             var ligands = [], bigincreases=0, increases = 0, bigdecreases=0, decreases = 0, unchanged=0, unspecified = 0;
                             
                             
                             $.each( val[0], function( key, v ) {
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


                             original_title = $('#'+plotid).find("#"+key).attr('original_title')

                             $('#'+plotid).find("#"+key).css("fill", color);
                             $('#'+plotid).find("#"+key).next().css("fill",color_letter );
                             $('#'+plotid).find("#"+key).attr('title',original_title+extra);
                             $('#'+plotid).find("#"+key+"t").attr('title',original_title+extra);


                          });
                        $("circle").tooltip('fixTitle');
                        $("text").tooltip('fixTitle');
    
                        
                    }

                    function ajaxInteractions(plotid,protein) {

                      resetColors(plotid);

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

                    function ajaxInteractionsPos(plotid) {

                      resetColors(plotid);


                      var pos = jQuery.parseJSON(interaction_json);

                          $.each( pos, function( key, val ) {

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
  
                    }

                    function construct_annotations(plotid) {

                      resetColors(plotid);

                      var pos = jQuery.parseJSON(annotations_json);

                          $.each( pos, function( key, val ) {

                            // var flags = [], falgsAA = [], output = [], outputAA = [], l = val.length, i;
                            // for( i=0; i<l; i++) {
                            //     if( flags[val[i][1]]) continue;
                            //     flags[val[i][1]] = true;
                            //     output.push(val[i][1]);
                            // }
                            // for( i=0; i<l; i++) {
                            //     if( flags[val[i][0]]) continue;
                            //     flags[val[i][0]] = true;
                            //     outputAA.push(val[i][0]);
                            // }
                             
                            //  extra = "\n" + String(val.length) + " interactions | Type: "+ output +" | Residue in crystal:"+ outputAA;
                             extra = "<br>"+val[1]; //.replace(/<br>/g, '&#013;');


                             $('#'+plotid).find("#"+key).css("fill", val[2]);
                             $('#'+plotid).find("#"+key).next().css("fill", val[3]);

                             original_title = $('#'+plotid).find("#"+key).attr('original_title')

                             $('#'+plotid).find("#"+key).attr('title',original_title+extra);
                             $('#'+plotid).find("#"+key+"t").attr('title',original_title+extra);


                          });
                        $("circle").tooltip('fixTitle');
                        $("text").tooltip('fixTitle');
  
                    }

                    function ajaxInteractionsLigand(protein,ligand) {

                        resetColors('snakeplot');
                        resetColors('helixbox');

                        $.getJSON( '/interaction/ajaxLigand/'+protein+'/'+ligand, function( data ) {
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


                             $('[id='+key+']').css("fill", "#E60A0A");
                             $('[id='+key+']').next().css("fill", "#FDFF7B");

                             original_title = $("#"+key).attr('original_title')

                             $('[id='+key+']').attr('title',original_title+extra);
                             $('[id='+key+'t]').attr('title',original_title+extra);


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

                    $( document ).ready(function() {
                      if ( $( "#cp2_helixbox" ).length ) {
                        $('#cp2_helixbox').colorpicker();
                        $('#cp2_snakeplot').colorpicker();
                      }
                    });

                    function getContrast50(hexcolor){
                         return (parseInt(hexcolor.replace('#', ''), 16) > 0xffffff/2) ? 'black':'white';
                    }
                        
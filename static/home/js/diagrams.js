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
                                    '-': ['#FFFFFF', '#000000']    
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

                    function toggleLoop(id,type) {
                        $(id+".long").toggle();
                        $(id+".short").toggle();
                        maxmin();

                        $(id+".long").toggle();
                        $(id+".short").toggle();

                        $(id+".long").fadeToggle();
                        $(id+".short").fadeToggle();
                    }

                    function applyPresentColors(target) {

                        //console.log( $('#'+target));

                        $('#'+target).find("circle").each(function( index ){
                              //console.log( index + ": " + $( this ).text() );
                              aa =  $(this).next().text();
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
                        classmax = ''
                        classmin = ''
                        $('#snake').children('text').each(function () {
                            if ($(this).is(":visible")) {
                                count = count +1;
                                y = parseInt($(this).attr( "y" ));
                                classtext = $(this).attr( "class" );
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
                        oldheight = $(svg).attr('height');
                        svg.setAttribute('height', (svgmax-svgmin+margin*2));
                        // $(svg).animate(
                          //  {"min-height": (svgmax-svgmin+margin*2)},
                          // {duration: 500,
                          //  step: function( top ){
                          //      this.setAttribute("height", "translate(0,"+Math.round(top)+")");
                          //    }
                         //   });

                        console.log('New height:'+ (svgmax-svgmin) +' old height:'+oldheight);
                        console.log("Prev attr"+$('#snake').attr("transform"));

                        $('#snake').attr("transform", "translate(0," + (-svgmin+margin) + ")");
                        
                        // $('#snake')
                        // .animate(
                         //  {"min-height": (-svgmin+margin)},
                         //  {duration: 500,
                         //   step: function( top ){
                         //       this.setAttribute("transform", "translate(0,"+top+")");
                          //    }
                          //  });

                        console.log("New attr"+$('#snake').attr("transform"));
                        // svgAsDataUri(document.getElementById("snakeplot"),{}, function(uri) {
                        //     console.log(uri);
                        //   $("body").append('<img src="' + uri + '" /><a href-lang="image/svg+xml" href="data:image/svg+xml;base64,\n'+uri+'" title="file.svg">download</a>');
                        //    $('body').append(
                        //     $('<a>')
                        //       .attr('href-lang', 'image/svg+xml')
                        //       .attr('href',  uri)
                        //       .attr('download',  'svg.svg')
                        //       .text('Download new')
                        //   );
                        // });
                    }

                    $( document ).ready(function() {    
                        var elements = document.getElementsByClassName('long')

                        for (var i = 0; i < elements.length; i++){
                            elements[i].style.display = 'none';
                        }
                        maxmin();


                        $('rect').each(function(){
                            
                            rectclass = $(this).attr('class');
                            if (rectclass) {
                                if (rectclass.indexOf("CL") >= 0 && rectclass.indexOf("long") >= 0) {

                                    numResidues = ($('.'+rectclass.replace(/ /g,".")).length-3)/2

                                    console.log('class:'+rectclass+' count'+numResidues);

                                    if (numResidues<10) {
                                        toggleLoop('.'+rectclass.split(' ')[0],'');
                                    }
                                }
                            }
                        });

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
                        newcolor = newcolor.split('-');

                      $(this).css("fill", newcolor[2]);
                      $(this).prev().css("fill", newcolor[1]);
                    });
                    
                    $(".rcircle").click(function() {
                        parentid = $(this).closest('svg').attr('id');
                        newcolor = $(".pick-color."+parentid+".selected").attr('id');
                        newcolor = newcolor.split('-');

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

                        $.getJSON( '/mutations/ajax/'+protein+'/', function( data ) {
                          $.each( data, function( key, val ) {

                             max = String(Math.max.apply(null, val));
                             min = String(Math.min.apply(null, val));
                             extra = "<br>" + String(val.length) + " mutations | "+ max +" maxFold | "+ min +" minFold";


                             original_title = $('#'+plotid).find("#"+key).attr('original_title')

                             $('#'+plotid).find("#"+key).css("fill", "#E60A0A");
                             $('#'+plotid).find("#"+key).next().css("fill", "#FDFF7B");
                             $('#'+plotid).find("#"+key).attr('title',original_title+extra);
                             $('#'+plotid).find("#"+key+"t").attr('title',original_title+extra);


                          });
                        $("circle").tooltip('fixTitle');
                        $("text").tooltip('fixTitle');
    
                        });
                    }

                    function ajaxInteractions(plotid,protein) {

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
                             
                             extra = "<br>" + String(val.length) + " interactions | Type: "+ output +" | Residue in crystal:"+ outputAA;


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



                    

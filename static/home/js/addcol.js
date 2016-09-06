

function addColumn(tableid, myindex) {
   var r_i = 0 ;
    $('#'+tableid+' tr').each(function(){       
          td = $(this).find("td").eq(myindex);
          $(this).append(td.clone(true));
          newtd = $(this).find("td").last();
          cur_index = newtd.index();
          newtd.children().each(function(){
            if ($(this).is('input:text')){
            $(this).val('');
          }
            cur_name = $(this).attr('name');
            new_name = cur_name+"_"+cur_index;
            $(this).attr('name');
            $(this).attr('name',new_name)
            $(this).attr('name'); 
          })
    });
}

function deleteColumn(tblId)
{
	var allRows = document.getElementById(tblId).rows;
	for (var i=0; i<allRows.length; i++) {
		if (allRows[i].cells.length > 2) {
			allRows[i].deleteCell(-1);   
		}
	}
}

function addaux(tableid){
   var r_i= 0 ;
      $("#"+tableid+" tr").each(function(){      
            td = $(this).find("td").eq(0);
            $(this).append(td.clone(true));      
            newtd = $(this).find("td").last();   
            newtd.show();
            cur_index = newtd.index();        
            if (newtd.hasClass("klon1")){
              newtd.attr("class", newtd.attr("class").replace("klon1", "klon"+cur_index));
            }
            newtd.children().each(function(){
             if( $(this).hasClass('hidetext') ){
                    $(this).hide();
                    if ($(this).is('input, textarea')){
                     $(this).val('');
                    }  
              }
            if( $(this).hasClass('signal') ){
              $(this).attr("class", $(this).attr("class").replace("signal", "signal duplicate"));
            }
            if ($(this).is('span,input,textarea')){ 
                if ($(this).hasClass('delcol')){
                    $(this).show();
                    $(this).val('-'+ $(this).parent().index() );
                }
                else{
                $(this).val(''); 
                }
            }
              cur_class = $(this).attr('class');
              new_class = cur_class+"_"+cur_index;
              cur_id = $(this).attr('id');
              new_id=cur_id+"_"+cur_index;
              cur_name=$(this).attr('name');
              new_name= cur_name+"_"+cur_index;
              $(this).attr('class',new_class);
              $(this).attr('name', new_name);
              $(this).attr('id', new_id);
            })
      });
}

function deleteRow(tblId){
 var allRows = document.getElementById(tblId).rows;
    if (allRows.length > 2) {
    $('#'+tblId+' tr:last').remove();
    }
}
  
function addRow(tblid){
var r_i= 0 ;
  $("#"+tblid+" tr").eq(0).show();
  $("#"+tblid+" tr").eq(1).clone(true).show().appendTo("#"+tblid);
  $("#"+tblid+" tr:last").each(function(){
   cur_index = $(this).index();
   newtd=$(this).find('td');   //find the cells for each row
   newtd.children().each(function(){
              if ($(this).is('span,input,textarea')){ 
                  if ($(this).hasClass('checkrow')){
                    button_class=$(this).attr('class');
                    new_button_class=button_class+" new_delrow row";
                    $(this).attr('class',new_button_class); 
                  } 
                  else if ($(this).hasClass('delrow')){ 
                    $(this).show();
                  $(this).val('-');
                  }
                 else{
                  $(this).val('');
                 }
              }
              if ($(this).hasClass('new_delrow')){
                $(this).show();
              }
              cur_class = $(this).attr('class');
              new_class = cur_class+"_"+cur_index;
              cur_id = $(this).attr('id');
              new_id=cur_id+"_"+cur_index;
              cur_name=$(this).attr('name');
              new_name= cur_name+"_"+cur_index; 
              $(this).attr('class',new_class)
              $(this).attr('name', new_name);
              $(this).attr('id', new_id);
            }) 
      });
}

function showOther(dropdown_id, other_class, drop_value){
    $("#"+dropdown_id).on('change', function () {
          if(this.value === drop_value){
            $('.'+other_class).show();  
          } 
          else {
            $("."+other_class).val('')
            $("."+other_class).hide(); 
          }
      });
}

function showOtherAux(dropdown_class, other_class, drop_value){
    $("."+dropdown_class).on('change', function () {
          this_index=$(this).parent().index();
          console.log(this_index);
      if (this_index>1){
          if(this.value === drop_value){
            $('.'+other_class+'.col_id_'+this_index).show();  
          } 
          else {
            $("."+other_class+'.col_id_'+this_index).val('')
            $("."+other_class+'.col_id_'+this_index).hide(); 
          }
      } 
      else{
         if(this.value === drop_value){
            $('.'+other_class+'.col_id').show();  
          } 
          else {
            $("."+other_class+'.col_id').val('')
            $("."+other_class+'.col_id').hide(); 
          }
      }
      });
}

function showOtherDynamic(dropdown_class, other_class, drop_value){
    $("."+dropdown_class).on('change', function () {
          this_index=$(this).parent().parent().index();
          console.log(this_index);
      if (this_index>1){
          if(this.value === drop_value){
            $('.'+other_class+'.row_id_'+this_index).show();  
          } 
          else {
            $("."+other_class+'.row_id_'+this_index).val('')
            $("."+other_class+'.row_id_'+this_index).hide(); 
          }
      } 
      else{
         if(this.value === drop_value){
            $('.'+other_class+'.row_id').show();  
          } 
          else {
            $("."+other_class+'.row_id').val('')
            $("."+other_class+'.row_id').hide(); 
          }
      }
      });
}

function delRow(button_class, table_id){
 $(button_class).on('click', function (){
      if (confirm("Are you sure you want to delete this row?")){
      $(this).closest('tr').remove ();   //this parent().parent().remove()
      $(table_id+' tr').each(function(){
        var my_index=$(this).index();
          $(this).children().each(function(){    
            $(this).children().each(function(){       
              if (my_index>1){
                var out_class=$(this).attr("class").split(' ').pop();
                //console.log("current_class"+out_class);
                $(this).removeClass(out_class);
                $(this).addClass("row_id_"+my_index); 
               if ($(this).hasClass("with_rec_val")){
                }
               else{
                $(this).attr("id", $(this).attr("id").replace(/\d+$/, my_index));
                $(this).attr("name", $(this).attr("name").replace(/\d+$/, my_index));
               }
            }
        });
     });
   });
  }
    else{
    }
  }); 
}

function delLast(button_class, table_id, ins_index){
   $(button_class).on('click', function (){
      if (confirm("Are you sure you want to delete this row?")){
          $(this).closest('tr').remove ();   //this parent().parent().remove()
          $(table_id+' tr').each(function(){
              var my_index=$(this).index();
              $(this).children().each(function(){    //td
                  $(this).children().each(function(){   //elements inside td
                         if (my_index>ins_index){
                      var out_class=$(this).attr('class').split(' ').pop();
                      $(this).removeClass(out_class);
                      $(this).addClass('row_id_'+(my_index-ins_index));
                     //update ids and names excluding those from the insertions table
                     if ($(this).hasClass('with_rec_val')){
                      }
                     else{
                      $(this).attr("id", $(this).attr("id").replace(/\d+$/, (my_index-ins_index)));
                      $(this).attr("name", $(this).attr("name").replace(/\d+$/, (my_index-ins_index) ));
                     }
                      }
                  });
              });
          });

        }
        else{
        }
   });
}

function addLast(button_id, table_id, row_index){
  $(button_id).on('click', function() {
    $(table_id+" tr").eq(row_index).clone(true).appendTo(table_id).addClass('newClass').find('th').html("");
      $(table_id+" tr:last").each(function(){
       cur_index = $(this).index()-row_index;       //so the counting starts from 1
       newtd=$(this).find('td');   //find the cells for each row
       newtd.children().each(function(){
                   if ($(this).is("label")){
                    $(this).hide();
                  }
                  if ($(this).hasClass("chem_enz_delrow")){
                    $(this).show();
                  }
                  if ($(this).hasClass("deterg_delrow")){
                    $(this).show();
                  }
                  if ($(this).hasClass("delrow_db")){
                    $(this).show();
                  }
                  if ($(this).hasClass("chem_enz_remark")){
                    $(this).hide();
                  }
                  if ($(this).hasClass("other_type_deterg")){
                    $(this).hide();
                  }
                  cur_class=$(this).attr('class');
                  cur_id = $(this).attr('id');
                  new_id=cur_id+"_"+cur_index;
                  cur_name=$(this).attr('name');
                  new_name= cur_name+"_"+cur_index;
                  new_class=cur_class+"_"+cur_index;
                  $(this).attr('name', new_name);
                  $(this).attr('id', new_id);
                  $(this).attr('class', new_class);
                })
     });
  });
}

function showOtherChem(dropdown_class, other_class, drop_value, last_row){
    $("."+dropdown_class).on('change', function () {
          this_index=$(this).parent().parent().index();
          console.log(this_index);
        if (this_index>last_row){
            if(this.value === drop_value){
              $('.'+other_class+'.row_id_'+(this_index-last_row)).show();  
            } 
            else {
              $("."+other_class+'.row_id_'+(this_index-last_row)).val('')
              $("."+other_class+'.row_id_'+(this_index-last_row)).hide(); 
            }
        }
        else{
           if(this.value === drop_value){
              $('.'+other_class+'.row_id').show();  
            } 
            else {
              $("."+other_class+'.row_id').val('');
              $("."+other_class+'.row_id').hide(); 
            }
        }
   });
}

//---udpate ids for auxiliary proteins
function UpdateIds(table_original, table_destination){
    //udpate the id's  (of all table and originals)
        $(table_original+' td').each(function(){
          var my_index=$(this).index();
          $(this).children().each(function(){    
              var out_class=$(this).attr('class').split(' ').pop();
              if (out_class==='col_id'){
              $(this).removeClass(out_class);
              $(this).addClass('col_id');
              }
             else {
              $(this).removeClass(out_class);
              $(this).addClass('col_id_'+my_index);
              }
              if (my_index>1){
              $(this).attr("id", $(this).attr("id").replace(/\d+$/, my_index));
              $(this).attr("name", $(this).attr("name").replace(/\d+$/, my_index));
              }
          });
       });  
       ///update also the cloned td ids classes etc
         $(table_original+' td').each(function(){
           var td_class=$(this).attr('class').split(' ').pop();
            if (td_class.match("^klon")){
              console.log(td_class);
              //$("."+td_class).hide();
              klon_td=$(this).index();
              $(table_destination+" td").each(function(){
                  if ( $(this).hasClass(td_class) ){
                    del_row_index=$(this).parent().index();
                      var del_class=$(this).attr('class').split(' ').pop();
                      $(this).children().each(function(){
                        //now my this is children of del td
                        //tr index  
                        $(this).attr("class", $(this).attr("class").replace(/\bins_del_single.*?\b/g, 'ins_del_single'+klon_td));
                        $(this).attr("class", $(this).attr("class").replace(/\bins_del_start.*?\b/g, 'ins_del_start'+klon_td));
                        $(this).attr("class", $(this).attr("class").replace(/\bins_del_end.*?\b/g, 'ins_del_end'+klon_td));
                        $(this).attr("class", $(this).attr("class").replace(/\bcol_id.*?\b/g, 'col_id_'+klon_td));
                        $(this).attr("class", $(this).attr("class").replace(/\brow_id.*?\b/g, 'row_id_'+del_row_index));
                        $(this).attr("id", $(this).attr("id").replace(/\d+$/, klon_td));

                        }); //children
                    $(this).attr("class", $(this).attr("class").replace(/\bklon.*?\b/g, 'klon'+klon_td));
                    $(this).attr("class", $(this).attr("class").replace(/\bklon.*?\b/g, 'klon'+klon_td));
       
                    //update the tr cloned insertion class
                    $(this).parent().each(function(){
                      $(this).attr("class", $(this).attr("class").replace(/\bcloned.*?\b/g, 'cloned_insertion'+klon_td));   
                      $(this).children().each(function(){
                        if ($(this).is('th')){
                          klon_td_minus=klon_td-1;
                          $(this).text("Insertion"+klon_td_minus); //update of th
                        }
                    });               
                    }); 
                  }
            }); ///deletions td
              $(this).attr("class", $(this).attr("class").replace(/\bklon.*?\b/g, 'klon'+klon_td));
           }
        });

}


// //----rearrange of auxiliary proteins based on C-term N-term
// function Rearrange(table) {
//    $(table' tr:first').each(function() {
//     var first_row = $(this);
//     var td1 = tr.find('td:eq(1)'); // indices are zero-based here
//     var td2 = tr.find('td:eq(3)');
//     td1.detach().insertAfter(td2);
// });
//     UpdateIds("#aux_proteins", "#deletions");
// }

//----rearrange of auxiliary proteins based on C-term N-term
function moveColumn(table, from, to) {
    var rows = jQuery('tr', table);
    var cols;
    rows.each(function() {
        cols = jQuery(this).children('th, td');
        cols.eq(from).detach().insertBefore(cols.eq(to));
    });
    UpdateIds("#aux_proteins", "#deletions");
}

//----rearrange of auxiliary proteins based on C-term N-term
function Nterm(){
   $('#aux_proteins tr:first').each(function() {
     // $(this).find('td').each(function() {
          var rows = jQuery('tr', "#aux_proteins");
          var rows_last=$(this).find("td:last").index();
          var cols;
          var n_term=$(this).find("select :selected[value='N-term']");
          
          var n_term_last=$(this).find("select :selected[value='N-term']").last();
          var n_term_first=$(this).find("select :selected[value='N-term']").first();
          var last_index=n_term_last.parent().parent().index();
          var first_index=n_term_first.parent().parent().index();

          var receptor=$(this).find("select :selected[value='Within Receptor']");
          var last_receptor=$(this).find("select :selected[value='Within Receptor']").last();
          var last_receptor_index=last_receptor.parent().parent().index();

          var ind_array = [];
          var c_array=[];
          i = 0;
           n_term.each(function(){
            var nterm_index=$(this).parent().parent().index();
            ind_array[i++] = nterm_index;
            //console.log($(this).parent().parent().index()) ;
           });
           
           var arr_i, len;
           for (arr_i=0, len=ind_array.length; arr_i < len; ++arr_i) {
              rows.each(function() {
                  cols = jQuery(this).children('th, td');
                  //cols.eq(first_index).detach().insertAfter(cols.eq(1));
                  cols.eq(ind_array[arr_i]).detach().insertAfter(cols.eq(1));
                   });

            // console.log(n_term.val());
            // console.log("index"+n_term.parent().parent().index());
          }
         });
      
      UpdateIds("#aux_proteins", "#deletions");
    }
          
//         function Cterm(){
//           var cterm=$(this).find("select :selected[value='C-term']");
//           var cterm_first=$(this).find("select :selected[value='C-term']").first();
//           var cterm_last=$(this).find("select :selected[value='C-term']").last();
//           var cterm_last_index=cterm_last.parent().parent().index();
//           var cterm_first_index=cterm_first.parent().parent().index();

//            cterm.each(function(){
//             var cterm_index=$(this).parent().parent().index();
//             c_array[i++] = cterm_index;
//            });

//           var carr_i, len_c;
//            for (carr_i=0, len_c=c_array.length; carr_i < len_c; ++carr_i) {
//                 rows.each(function() {
//                     cols = jQuery(this).children('th, td');
//                     //cols.eq(cterm_first_index).detach().insertAfter(cols.last());
//                     //cols.eq(cterm_first_index).detach().insertAfter(cols.eq());
//                     console.log("Edw"+last_index);
//                     cols.eq(c_array[carr_i]).detach().insertAfter(cols.eq(last_index));
//                 });
              
//               // console.log("index"+n_term.parent().parent().index());
//            //}); 
//           }
//           //console.log(rows_last);

//           // console.log(cterm_last_index);
        
          
// //     //  $(this).find("td[value=""]").attr('disabled', true);
// //     // var td1 = tr.find('td:eq(1)'); // indices are zero-based here
// //     // var td2 = tr.find('td:eq(3)');
// //     // td1.detach().insertAfter(td2);
// // });
// UpdateIds("#aux_proteins", "#deletions");
// }

 function WithinReceptor(){
    Nterm();
   $('#aux_proteins tr:first').each(function() {
     // $(this).find('td').each(function() {
          var rows = jQuery('tr', "#aux_proteins");
          var rows_last=$(this).find("td:last").index();
          var cols;

          var n_term_last=$(this).find("select :selected[value='N-term']").last();
          var last_index=n_term_last.parent().parent().index();
          
          var receptor=$(this).find("select :selected[value='Within Receptor']");
          var last_receptor=$(this).find("select :selected[value='Within Receptor']").last();
          var last_receptor_index=last_receptor.parent().parent().index();

          var rec_array = [];
          i = 0;
           receptor.each(function(){
            var rec_index=$(this).parent().parent().index();
            rec_array[i++] = rec_index;
            //console.log($(this).parent().parent().index()) ;
           });

          if (n_term_last.length>0){
             var rec_arr_i, len_rec;
             for (rec_arr_i=0, len_rec=rec_array.length; rec_arr_i < len_rec; ++rec_arr_i) {
                rows.each(function() {
                    cols = jQuery(this).children('th, td');
                    //cols.eq(first_index).detach().insertAfter(cols.eq(1));
                    cols.eq(rec_array[rec_arr_i]).detach().insertAfter(cols.eq(last_index));
                     });

              // console.log(n_term.val());
              // console.log("index"+n_term.parent().parent().index());
            }
         }
         else{}
         });
      UpdateIds("#aux_proteins", "#deletions");
    }

function ValidateForm(){
  $("input:visible:not(.optional)").each(function(){
    blank_id=$(this).attr("id");
      if ($(this).val() === "" && !$("#error-"+blank_id).length){
        console.log("missing in "+blank_id);
        blank_id=$(this).attr("id");
        blank_name=$(this).attr("name");
        console.log("nothing in "+blank_id);
        var outer_class=$(this).attr("class").split(' ').pop();
        $(this).after("<label name=error-"+blank_name+" id=error-"+blank_id+" class='error "+outer_class+"'>*</label>");
      } else if ($(this).val() != "" && $("#error-"+blank_id).length) {
        $("#error-"+blank_id).remove();
      }
  });
// after("<tr><th class='cloned_th"+insertion.index()+"'>Insertion"+actual_index+"</th></tr>");
  $("select:visible:not(.optional)").each(function(){
    blank_select_id=$(this).attr("id");
    blank_select_name=$(this).attr("name");
    var outer_sel_class=$(this).attr("class").split(' ').pop();
      if ($(this).val() === "" && $("#error-"+blank_select_id).length == 0){
        $(this).after("<label name=error-"+blank_select_name+" id=error-"+blank_select_id+" class='error "+outer_sel_class+"'>*</label>");
      }
  });
}

function FieldsRequired(){
    if ($("label.error:visible").length>0){
      $("span.required_ind").show(); 
      $("html, body").animate({ scrollTop: $(".span").first() });
    }
    else{
      $("span.required_ind").hide(); 
    }
}

//handle validation error labels 
$('input:not(.optional):not(.unit)').on('keyup', function () {
    this.value === '' ? $(this).nextAll('label.error').first().show() : $(this).nextAll('label.error').first().hide();
});

$('input.unit:not(.optional)').on("change",function () {
    this.value === '' ? $(this).nextAll('label.error').first().show() : $(this).nextAll('label.error').first().hide();
});

$('select:not(.optional)').on('change', function () {
    var select_id=$(this).attr("id");
    this.value === '' ? $("#error-"+select_id).show() : $("#error-"+select_id).hide();
});

function LabelErrors(){
//makes sure that "label" messages are shown only when fields are visible
    $("label.error").each( function(){
    // if ($(this).is("label")){
        var original_label= $(this).attr("id");
        var input_id= original_label.split("-").pop();
        console.log(input_id);
       //  var error_id = original.substr(0, original_label.indexOf('-')); 
       if ($("#"+input_id).is(":visible")){  
       } 
        else{
         $("#"+original_label).hide();
        }
      // }
    });
}

////detach prevention 
// function Submission(){
//   ValidateForm();
//   var tade=ValidateForm();
//   console.log(tade);
// if (tade!=false){
//     $("#xtals_form").unbind("submit");
//       // 
      
//     // $("#id_name_cont").on("click", function(){    
//     }
// }
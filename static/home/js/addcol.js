
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
            //newtd.show();
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
              $(this).attr('class',new_class)
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
  $("#"+tblid+" tr").eq(1).clone(true).appendTo("#"+tblid);
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

function delRow(button_class, table_id){
 $(button_class).on('click', function (){
      if (confirm("Are you sure you want to delete this row?")){
      $(this).closest('tr').remove ();   //this parent().parent().remove()
      $(table_id+' tr').each(function(){
        var my_index=$(this).index();
          $(this).children().each(function(){    
            $(this).children().each(function(){       
              if (my_index>1){
                var out_class=$(this).attr('class').split(' ').pop();
                $(this).removeClass(out_class);
                $(this).addClass('row_id_'+my_index); 
               if ($(this).hasClass('with_rec_val')){
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

function delLast(button_class, table_id){
   $(button_class).on('click', function (){
      if (confirm("Are you sure you want to delete this row?")){
          $(this).closest('tr').remove ();   //this parent().parent().remove()
          $(table_id+' tr').each(function(){
              var my_index=$(this).index();
              $(this).children().each(function(){    //td
                  $(this).children().each(function(){   //elements inside td
                         if (my_index>4){
                      var out_class=$(this).attr('class').split(' ').pop();
                      $(this).removeClass(out_class);
                      $(this).addClass('row_id_'+(my_index-4));
                     //update ids and names excluding those from the insertions table
                     if ($(this).hasClass('with_rec_val')){
                      }
                     else{
                      $(this).attr("id", $(this).attr("id").replace(/\d+$/, (my_index-4)));
                      $(this).attr("name", $(this).attr("name").replace(/\d+$/, (my_index-4) ));
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
                  if ($(this).hasClass("chem_enz_delrow")){
                    $(this).show();
                  }
                  if ($(this).hasClass("chem_enz_remark")){
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

function showOtherChem(dropdown_class, other_class, drop_value){
    $("."+dropdown_class).on('change', function () {
          this_index=$(this).parent().parent().index();
          console.log(this_index);
        if (this_index>4){
            if(this.value === drop_value){
              $('.'+other_class+'.row_id_'+(this_index-4)).show();  
            } 
            else {
              $("."+other_class+'.row_id_'+(this_index-4)).val('')
              $("."+other_class+'.row_id_'+(this_index-4)).hide(); 
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


$(document).ready(function () {

   $('.other_fusion').mousemove(function(e){ 
       var hovertext = $(this).attr('hinttext');  
         $('#hintbox').text(hovertext).show();
    $('#hintbox').css('top',e.clientY+290).css('left',e.clientX+15);})
     .mouseout(function(){
    $('#hintbox').hide();
 });

addaux('aux_proteins');
addRow('chem_comp');

console.log("bind to submit");
console.log($("#xtals_form"));

$("#xtals_form").submit(function (e) {
  // $("html, body").animate({ scrollTop: $(".error").first() });
  // ValidateForm(); 
 ValidateForm(); 
 FieldsRequired();
  if ($("label.error:visible").length>0){
  e.preventDefault();
  }
  else{
    $("#xtals_form").unbind("submit");
  }
// $("#id_name_cont").on("click", function(){    
});

$("#btn_s").on("click", function(){
//   Submission();
 ValidateForm(); 
 FieldsRequired();
});

$(document).on("change", function(){
  LabelErrors();
  // FieldsRequired();
});

 $('.datepicker').datepicker({
     inline: true,
    showOtherMonths: true
 });



$('.numeric').on('keydown',  function(e){
    -1!==$.inArray(e.keyCode,[46,8,9,27,13,110,190])||/65|67|86|88/.test(e.keyCode)&&(!0===e.ctrlKey||!0===e.metaKey)||35<=e.keyCode&&40>=e.keyCode||(e.shiftKey||48>e.keyCode||57<e.keyCode)&&(96>e.keyCode||105<e.keyCode)&&e.preventDefault()});

$("#id_crystal_type").on("change", function() {
  if (this.value === "lipidic cubic phase (LCP)" ){
      $(".lcp").addClass("optional opt").show();     
    }
    else{
      $(".lcp").hide();
    }
});

//--WT,MUT AA cannot have the same value 
$(".wild").on('change', function() {
  row_id=$(this).attr('class').split(' ').pop();
  my_id=$(this).attr("id")
  selected=$("#"+my_id+" :selected").text();
  $(".mutant."+row_id).each(function (){
    console.log("h aksia "+selected);
    console.log($(this).attr("id"));
    $(this).find("option").attr("disabled", false); //reset by enabling all 
   $(this).find("option[value="+ selected +"]").attr('disabled', true);

  }); 
});

$(".mutant").on('change', function() {
  row_id=$(this).attr('class').split(' ').pop();
  my_id=$(this).attr("id")
  selected=$("#"+my_id+" :selected").text();
  $(".wild."+row_id).each(function (){
    $(this).find("option").attr("disabled", false); //reset by enabling all 
   $(this).find("option[value="+ selected +"]").attr('disabled', true);

  }); 
});

//--make optional visible to the user
$('input').on('keyup', function () {
    $("#xtals_form input.opt").filter(function () {
        this.value === '' ? $(this).addClass('optional') : $(this).removeClass('optional')
    });
});

//--if ratio on value then the unit is optional
$("input.half").on("keyup", function() {
    str=$(this).val();
  if (str.indexOf(":") >= 0){
    pl=$(this).nextAll("input.unit").attr('placeholder');
    $(this).nextAll("input.unit").removeAttr('placeholder');
    //$(this).next("input.unit").attr("placeholder", "ratio");
    $(this).nextAll("input.unit").val('ratio').trigger("change");
    // $(this).nextAll("input.unit").removeClass('error');
     // $(this).nextAll("label.error").first().hide();
    $(this).nextAll("input.unit").keydown(function(){
      return false;
    });
    //make sure that it is not necessary
    console.log(pl);

  }
else{
  $(this).nextAll("input.unit").val('').trigger("change");
   $(this).nextAll("input.unit").attr("placeholder", "unit i.e.:%w/v");  //this doesn't apply to id protein conc_unit, adjust to that
   $(this).nextAll("input.unit").off("keydown");
   if ($(this).hasClass("optional")){
   $(this).nextAll("input.unit").addClass("optional");
   // $(this).nextAll("label.error").first().show();
   }
}
});

//---------------------------------------------------------------------------------------------------
      $('.aamod_pos_type').on('change', function () {
        temp_index = $(this).parent().parent().index();     //parent of td is tr
        var aamod_type= ["","single", "pair", "range"];
      for (var k=0; k<aamod_type.length; k+=1){
        var aa=aamod_type[k] ;
          if(this.value === aa){
          if (temp_index>1) {
                    $('.aa_type.row_id_'+temp_index).val(''); 
                    $('.aa_type.row_id_'+temp_index).hide(); 
                    $('.'+aa+'.row_id_'+temp_index).show(); 
                  } else {
                    $('.aa_type.row_id_'+temp_index).val(''); 
                    $('.aa_type.row_id').hide(); 
                    $('.'+aa+'.row_id').show(); 
                  }
          }
     }
  });

$(".insert_pos_type").on('change', function () {
        temp_index = $(this).parent().index();     //parent of td is tr
        var insert_type= ["","ins_single", "ins_range"];
      for (var k=0; k<insert_type.length; k+=1){
        var ins=insert_type[k] ;
          if(this.value === ins){
          if (temp_index>1) {
                    $(".ins_pos_type.col_id_"+temp_index).val('');
                    $(".ins_pos_type.col_id_"+temp_index).hide(); 
                    $('.'+ins+'.col_id_'+temp_index).show(); 
                  } else {
                    $(".ins_pos_type.col_id").val(''); 
                    $(".ins_pos_type.col_id").hide(); 
                    $('.'+ins+'.col_id').show(); 
                  }
          } 
      }
  });

$(".position").on('change', function () {
        temp_index = $(this).parent().index();     
          if(this.value === 'Within Receptor'){
            if (temp_index>1) {
                      $(".with_rec.col_id_"+temp_index).show(); 
                      $(".within_receptor_th").show();
                    } else {
                      $(".with_rec.col_id").show(); 
                      $(".within_receptor_th").show();
                    }
          }   
         else{
             if (temp_index>1) {
                      $(".with_rec.col_id_"+temp_index).prop('selectedIndex',0); 
                      $(".with_rec.col_id_"+temp_index).hide(); 
                      $(".with_rec_val.col_id_"+temp_index).val('');
                      $(".with_rec_val.col_id_"+temp_index).hide();
                      $(".stable.col_id_"+temp_index).show();
                    } 
                    else { 
                      $(".with_rec.col_id").prop('selectedIndex',0); 
                      $(".with_rec.col_id").hide(); 
                      $(".with_rec_val.col_id").val(''); 
                      $(".with_rec_val.col_id").hide();
                      $(".stable.col_id").show();
                    }
         }
  });

var proteins=[ "", "signal", "tag", "fusion", "linker", "prot_cleavage"];
            $(".protein_type").on('change', function () {
             temp_index = $(this).parent().index();
              for (pos=0; pos<proteins.length; pos+=1){
                var i=proteins[pos] ;
                  if(this.value === i){
                    $(".sub_type").show();       
                    if (temp_index>1) {
                      $('.prot_type.col_id_'+temp_index).hide(); 
                      $('.prot_type.col_id_'+temp_index).val(""); 
                      $('.others.col_id_'+temp_index).val("");
                      $('.linker.col_id_'+temp_index).val("");
                      $('.'+i+'.col_id_'+temp_index).show();   
                    } else {
                      $('.prot_type.col_id').hide(); 
                      $('.'+i+'.col_id').show();   
                      $(".prot_type.col_id").val('Please Select');   
                      $('.others.col_id').val("");
                      $('.linker.col_id').val("");
                    }
                  }
             }
        });

addLast("#addurl", "#contact_information", "5");
addLast("#add_treatment", "#solubil_purif", "3");
addLast("#add_deterg", "#solubil_deterg", "0");
// addLast("#add_treatment", "#solubil_purif");

$("#deleteurl").on('click', function() {
  $("#contact_information tr:last").each(function(){
    if ($(this).hasClass("newClass")){
      $(this).remove();
    }
 });
 });

$(".ph").on('change', function () {
  var ph_type= ["", "single_ph", "range_ph"];
  for (var n=0; n<ph_type.length; n+=1){
      var ph_val=ph_type[n] ;
      if(this.value === ph_val){
        $(".ph_type").val('');
        $(".ph_type").hide();
        $("."+ph_val).show();
      }
  }
});

$(".deletion_type").on('change', function () {
        temp_index = $(this).parent().parent().index();     //parent of td is tr
        var delet_type= ["","del_single", "del_range"];
      for (var l=0; l<delet_type.length; l+=1){
        var del=delet_type[l] ;
          if(this.value === del){    
          if (temp_index>1) {
                    $(".del_type.row_id_"+temp_index).val('');
                    $(".del_type.row_id_"+temp_index).hide(); 
                    $('.'+del+'.row_id_'+temp_index).show();
                  } else {
                    $(".del_type.row_id").val(''); 
                    $(".del_type.row_id").hide(); 
                    $('.'+del+'.row_id').show(); 
                  }
          }
     }
  });

showOtherDynamic('mod_other', 'other_mod', 'other')
showOther('id_ligand_activity','ligand_act_oth','other');
showOther('id_ligand_id_type','other_ligand','Other');
showOther('id_detergent','other_det','other [See next field]');
showOther('id_lcp_lipid','other_lcp','other [See next field]');
showOther('id_crystal_type','other_cryst_type','other [See next field]');
showOther('id_lipid','other_lipid','other [See next field]');
showOther('id_crystal_method','other_method','other [See next field]');
showOther('id_expr_method','other_expr','Other [In case of E.Coli or Yeast recombinant expression]');
showOther('id_host_cell_type','other_host_cell','other [See next field]');
showOther('id_host_cell','host_cell_other','other [See next field]');

showOtherAux('signal','other_signal','Other');
showOtherAux('prot_cleavage','other_prot','Other');
showOtherAux('tag','other_tag','Other');
showOtherAux('fusion','other_fusion','Other');

showOtherChem('det_type', 'other_type_deterg', 'other [See next field]','0');
showOtherChem('chem', 'chem_enz_remark', 'Other [See remark]','3');
//call delrow for the tables
delRow(".del_delrow", "#deletions");
delRow(".mut_delrow", "#mutations");
delRow(".chem_delrow", "#chem_comp");
delRow(".mod_delrow", "#modifications");

delLast(".chem_enz_delrow", "#solubil_purif", "3");
delLast(".delrow_db", "#contact_information", "5");
delLast(".deterg_delrow", "#solubil_deterg", "0");


$('.delcol').on('click', function (){
    col_index=$(this).parent().index();      

    if (confirm("Are you sure you want to delete this column?")){
    $("#aux_proteins tr").each(function() {
      $(this).children().eq(col_index).remove();

          console.log("col_index:"+col_index);
          $("#deletions tr").each(function() {
            if ($(this).hasClass("cloned_insertion"+col_index)){
                $(this).remove();
            }
          });
  
    });  

    UpdateIds("#aux_proteins", "#deletions");

  }
  else{}

 });

$(".position").on('click',function(){ 
   var my_index=$(this).parent().index();
   console.log("my_index="+my_index);
if (my_index>1){
    $(".insert_pos_type.col_id_"+my_index).one('click',function(){ 
        cur_index=$(this).parent().index();
        actual_index=cur_index-1;
        console.log("my_index="+cur_index);
        console.log("cur_index="+cur_index);
              var insertion= $(this).closest('td');
              insertion.addClass("klon"+insertion.index());
              var insertion_clone=$(insertion).clone(true);
              insertion_clone.find(".insert_pos_type").remove();
              $("#deletions tr:last").after("<tr><th class='cloned_th"+insertion.index()+"'>Insertion"+actual_index+"</th></tr>");
              $("#deletions tr:last").append(insertion_clone);
              var last_del_index=$("#deletions tr:last").index();
              insertion_clone.parent().addClass("cloned_insertion"+cur_index);
              //avoid duplicates
              $('.cloned_insertion'+cur_index).not(':first').remove();
              $(insertion_clone).each(function(){
                $(this).children().each(function(){
                $(this).addClass("row_id_"+last_del_index); 
                $(this).attr("id", $(this).attr("id").replace(/\d+$/, cur_index));   
                $(this).attr("name", $(this).attr("name").replace(/\d+$/, cur_index));       
                if ($(this).hasClass("ins_start")){
                  $(this).attr("class", $(this).attr("class").replace("ins_start", "ins_del_start"+cur_index));
                }
                else if ($(this).hasClass("ins_end")){
                  $(this).attr("class", $(this).attr("class").replace("ins_end", "ins_del_end"+cur_index));
                }
                else if ($(this).hasClass("ins_single")){
                  $(this).attr("class", $(this).attr("class").replace("to_change", "ins_del_single"+cur_index));
                }
                  });
              });
              $(".ins_start.col_id_"+cur_index).on("keyup paste", function() {
              $(".ins_del_start"+$(this).parent().index()).val($(this).val());
              });
             $('.ins_del_start'+cur_index).keydown(function() {
             return false;
              });
              $(".ins_end.col_id_"+cur_index).on("keyup paste", function() {
              $(".ins_del_end"+$(this).parent().index()).val($(this).val());
              });
             $('.ins_del_end'+cur_index).keydown(function() {
             return false;
              });
              $(".ins_single.col_id_"+cur_index).on("keyup paste", function() {
              $(".ins_del_single"+$(this).parent().index()).val($(this).val());
              });
              $('.ins_del_single'+cur_index).keydown(function() {
              return false;
              });
        });
}
else{
   $(".insert_pos_type.col_id").one('click',function(){ 
      first_index=$(this).parent().index();
      console.log("first index="+first_index);
            var insertion= $(this).closest('td');
            var insertion_clone=$(insertion).clone();
            insertion_clone.find(".insert_pos_type").remove();
            $("#deletions tr:last").after("<tr><th>Insertion"+first_index+"</th></tr>")
            $("#deletions tr:last").append(insertion_clone);
            var del_index=$("#deletions tr:last").index();
            insertion_clone.parent().addClass("cloned_insertion");
             $('.cloned_insertion').not(':first').remove();
            $(insertion_clone).each(function(){
            $(this).children().each(function(){
            $(this).addClass("row_id_"+del_index);      
            if ($(this).hasClass("ins_start")){
              $(this).attr("class", $(this).attr("class").replace("ins_start", "ins_del_start"+first_index));
            }
           else if ($(this).hasClass("ins_end")){
              $(this).attr("class", $(this).attr("class").replace("ins_end", "ins_del_end"+first_index));
            }
           else if ($(this).hasClass("ins_single")){
              $(this).attr("class", $(this).attr("class").replace("to_change", "ins_del_single"+first_index));
            }
                }); 
        });
            $(".ins_start.col_id").on("keyup paste", function() {
            $(".ins_del_start"+$(this).parent().index()).val($(this).val());
            });
           $('.ins_del_start'+first_index).keydown(function() {
           return false;
            });
            $(".ins_end.col_id").on("keyup paste", function() {
            $(".ins_del_end"+$(this).parent().index()).val($(this).val());
            });
           $('.ins_del_end'+first_index).keydown(function() {
           return false;
            });
            $(".ins_single.col_id").on("keyup paste", function() {
            $(".ins_del_single"+$(this).parent().index()).val($(this).val());
            });
            $('.ins_del_single'+first_index).keydown(function() {
            return false;
            });
      });
  }
});

$(".position").on('change', function(){
  var grab_index=$(this).parent().index();
  console.log("grabbed index="+grab_index);
 if ( $(".insert_pos_type.col_id_"+grab_index).is(":visible") ){
            $(".cloned_insertion"+grab_index).show();
          }
 else {
  $(".cloned_insertion"+grab_index).hide();
  //$("th").hasClass("cloned_th"+grab_index).hide();
 } 
});

$('.checked').on('click', function (){
      if (confirm("Are you sure you want to delete these rows?")){
           $('input:checked').each(function(){
            if ($(this).hasClass('checkrow')){
             $(this).closest('tr').remove();
          $("#mutations"+' tr').each(function(){
            var my_index=$(this).index();
            $(this).children().each(function(){    
              $(this).children().each(function(){   
                  if (my_index>1){
                  var out_class=$(this).attr('class').split(' ').pop();
                  $(this).removeClass(out_class);
                  $(this).addClass('row_id_'+my_index);
                 //update ids and names 
                  $(this).attr("id", $(this).attr("id").replace(/\d+$/, my_index));
                  $(this).attr("name", $(this).attr("name").replace(/\d+$/, my_index));
                      }
                  });
                });
              });
            }
        });
     }
 });

//!!!!! closing of document ready function
});










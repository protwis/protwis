function toggleButtonClass(button_id) {
  $('#'+button_id).toggleClass('active')
}
$("html").css('pointer-events','none');
$("html").css('-webkit-user-select','none'); /* Safari */
$("html").css('-webkit-user-select','none'); /* IE 10 and IE 11 */
$("html").css('user-select','none'); /* IE 10 and IE 11 */

function queryByLigandUnlockClick() {
  $("html").css('pointer-events','');
  $("html").css('-webkit-user-select',''); /* Safari */
  $("html").css('-webkit-user-select',''); /* IE 10 and IE 11 */
  $("html").css('user-select',''); /* IE 10 and IE 11 */
}

$(document).ready(function() {
  var error = false;
  $(".error-msg-p").each(function (i,e) {
    
    if ($(e).text() != '') {
      showAlert("No results obtained. Scroll down to the corresponding section for details.","info");
      error = true;
      var hash = '#';
      switch ($(e).attr('id')) {
        case 'ligand_bulk_search_by_name_error_msg':
          hash = "#by-name";
          break;
        case 'ligand_bulk_search_by_names_error_msg':
          hash = "#by-names";
          break;
        case 'ligand_bulk_search_by_id_error_msg':
          hash = "#by-id";
          break;
        case 'ligand_bulk_search_error_msg':
          hash = "#bulk-search";
          break;
        case 'ligand_structural_search_error_msg':
          hash = "#simsub-search";
          break;      
        default:
          break;
      }

      setTimeout(() => {
        location.hash = hash;
        $(document).ready(function() {

          queryByLigandUnlockClick();
        });
      }, 500);

    } 

  });

  if (!error) {
    queryByLigandUnlockClick();
  }

  $("#selection-autocomplete-ligand-by-name")
  .on('change',function (e) {$("#ligand_bulk_search_by_name_error_msg").text('');});

  $("#selection_ligand_bulk_search_by_names_textarea")
  .on('change',function (e) {$("#ligand_bulk_search_by_names_error_msg").text('');});

  $("#selection_ligand_bulk_search_by_id_textarea, #selection_ligand_by_id_fields")
  .on('change',function (e) {$("#ligand_bulk_search_by_id_error_msg").text('');});

  $("#selection_ligand_bulk_search_textarea, #selection_ligand_bulk_search_search_type, #selection_ligand_bulk_search_stereochemistry")
  .on('change',function (e) {$("#ligand_bulk_search_error_msg").text('');});

  $("#ligand_structural_search_smiles,#ligand_structural_search_search_type,#ligand_structural_search_input_type,#ligand_structural_search_similarity_threshold,#ligand_structural_search_stereochemistry")
  .on('change',function (e) {$("#ligand_structural_search_error_msg").text('');});

});

function submitSelectionLigandBulkSearchByName(url) {


  toggleButtonClass('selection-button-ligand-by-name');

  let search_params_data = {
    search_text : $("#selection-autocomplete-ligand-by-name").val().trim(),
    search_type : 'name',
  }

  // set CSRF csrf_token
  $.ajaxSetup({
    headers:
    { "X-CSRFToken": csrf_token }
  });

  // Submit search parameters
  $.post("/ligand/read_input_ligand_bulk_search", search_params_data,  function (data) {
    setTimeout(function(){
      toggleButtonClass('selection-button-ligand-by-name');
      window.location.href = url;
    }, 500);
    
  })
  .fail(
    function(jqXHR, textStatus, errorThrown){
      if (jqXHR.hasOwnProperty('responseJSON')) {
        if (jqXHR.responseJSON.hasOwnProperty('status_code')) {
          if (jqXHR.responseJSON.status_code === 400) {
            if (jqXHR.responseJSON.hasOwnProperty('msg')) {
              location.hash = '#';
              location.reload();
              return;
            }
          }
        }
          
      }
      toggleButtonClass('selection-button-ligand-by-name');
      showAlert("Something went wrong, please try again or contact us.", "danger");
    });


}

function submitSelectionLigandBulkSearchByNames(url) {


  toggleButtonClass('selection-button-ligand-by-names');

  let search_params_data = {
    search_text : $("#selection_ligand_bulk_search_by_names_textarea").val().trim(),
    search_type : 'names',
  }

  // set CSRF csrf_token
  $.ajaxSetup({
    headers:
    { "X-CSRFToken": csrf_token }
  });

  // Submit search parameters
  $.post("/ligand/read_input_ligand_bulk_search", search_params_data,  function (data) {
    setTimeout(function(){
      toggleButtonClass('selection-button-ligand-by-names');
      window.location.href = url;
    }, 500);
    
  })
  .fail(
    function(jqXHR, textStatus, errorThrown){
      if (jqXHR.hasOwnProperty('responseJSON')) {
        if (jqXHR.responseJSON.hasOwnProperty('status_code')) {
          if (jqXHR.responseJSON.status_code === 400) {
            if (jqXHR.responseJSON.hasOwnProperty('msg')) {
              location.hash = '#';
              location.reload();
              return;
            }
          }
        }
          
      }
      toggleButtonClass('selection-button-ligand-by-names');
      showAlert("Something went wrong, please try again or contact us.", "danger");
    });


}

function submitSelectionLigandBulkSearchById(url) {
  toggleButtonClass('selection-button-ligand-by-id');
  let search_params_data = {
    search_text : $("#selection_ligand_bulk_search_by_id_textarea").val().trim(),
    search_type : 'id',
    field: $("#selection_ligand_by_id_fields").val(),
  }


  // set CSRF csrf_token
  $.ajaxSetup({
    headers:
    { "X-CSRFToken": csrf_token }
  });

  // Submit search parameters
  $.post("/ligand/read_input_ligand_bulk_search", search_params_data,  function (data) {
    setTimeout(function(){
      toggleButtonClass('selection-button-ligand-by-id');
      location.hash = '#';
      window.location.href = url;
    }, 500);

  })
  .fail(
    function(jqXHR, textStatus, errorThrown){
      if (jqXHR.hasOwnProperty('responseJSON')) {
        if (jqXHR.responseJSON.hasOwnProperty('status_code')) {
          if (jqXHR.responseJSON.status_code === 400) {
            if (jqXHR.responseJSON.hasOwnProperty('msg')) {
              location.hash = '#';
              location.reload();
              return;
            }
          }
        }
          
      }
      toggleButtonClass('selection-button-ligand-by-id');
      showAlert("Something went wrong, please try again or contact us.", "danger");
    });


}






function submitSelectionLigandBulkSearch(url) {
  toggleButtonClass('selection-button');
  let search_params_data = {
    selection_ligand_chemical_structure_search_search_type_selection : $("#selection_ligand_chemical_search_type").val(),
    search_text : $("#selection_ligand_bulk_search_textarea").val().trim(),
    search_type : $("#selection_ligand_bulk_search_search_type").val(),
    stereochemistry : $("#selection_ligand_bulk_search_stereochemistry").is(":checked"),
  }
  

  // set CSRF csrf_token
  $.ajaxSetup({
    headers:
    { "X-CSRFToken": csrf_token }
  });

  // Submit search parameters
  $.post("/ligand/read_input_ligand_bulk_search", search_params_data,  function (data) {
    setTimeout(function(){
      toggleButtonClass('selection-button');
      window.location.href = url;
    }, 500);
  })
  .fail(
    function(jqXHR, textStatus, errorThrown){
      if (jqXHR.hasOwnProperty('responseJSON')) {
        if (jqXHR.responseJSON.hasOwnProperty('status_code')) {
          if (jqXHR.responseJSON.status_code === 400) {
            if (jqXHR.responseJSON.hasOwnProperty('msg')) {
              location.hash = '#';
              location.reload();
              return;
            }
          }
        }
          
      }
      toggleButtonClass('selection-button');
      showAlert("Something went wrong, please try again or contact us.", "danger");
    });


}

function submitLigandStructuralSearch(url) {
    let search_params_data = {
      selection_ligand_chemical_structure_search_search_type_selection : $("#selection_ligand_chemical_search_type").val(),
      smiles : $("#ligand_structural_search_smiles").val().trim(),
      search_type : $("#selection_ligand_chemical_search_type").val(),
      input_type : $("#ligand_structural_search_input_type").val(),
      similarity_threshold : $("#ligand_structural_search_similarity_threshold").val().trim(),
      stereochemistry : $("#ligand_structural_search_stereochemistry").is(":checked"),
    }
    toggleButtonClass('selection-button-smiles-search');

    // set CSRF csrf_token
    $.ajaxSetup({
      headers:
      { "X-CSRFToken": csrf_token }
    });

    // Submit search parameters
    $.post("/ligand/read_input_ligand_structural_search", search_params_data,  function (data) {
      setTimeout(function(){
        toggleButtonClass('selection-button-smiles-search');
        window.location.href = url;
      }, 500);
    })
    .fail(
      function(jqXHR, textStatus, errorThrown){
        if (jqXHR.hasOwnProperty('responseJSON')) {
          if (jqXHR.responseJSON.hasOwnProperty('status_code')) {
            if (jqXHR.responseJSON.status_code === 400) {
              if (jqXHR.responseJSON.hasOwnProperty('msg')) {
                location.hash = '#';
                location.reload();
                return;
              }
            }
          }
            
        }
        toggleButtonClass('selection-button-smiles-search');
        showAlert("Something went wrong, please try again or contact us.", "danger");
      });


}
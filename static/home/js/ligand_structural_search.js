function toggleButtonClass(button_id) {
  $('#'+button_id).toggleClass('active')
}

function submitSelectionLigandBulkSearchByName(url) {


  toggleButtonClass('selection-button-ligand-by-name');

  let search_params_data = {
    search_text : $("#selection-autocomplete-ligand-by-name").val().trim(),
    search_type : 'name',
    // limit : $("#ligand_structural_search_limit").val().trim()
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

function submitSelectionLigandBulkSearchById(url) {
  toggleButtonClass('selection-button-ligand-by-id');
  let search_params_data = {
    search_text : $("#selection_ligand_bulk_search_by_id_textarea").val().trim(),
    search_type : 'id',
    field: $("#selection_ligand_by_id_fields").val(),
    // limit : $("#ligand_structural_search_limit").val().trim()
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
      window.location.href = url;
    }, 500);

  })
  .fail(
    function(jqXHR, textStatus, errorThrown){
      if (jqXHR.hasOwnProperty('responseJSON')) {
        if (jqXHR.responseJSON.hasOwnProperty('status_code')) {
          if (jqXHR.responseJSON.status_code === 400) {
            if (jqXHR.responseJSON.hasOwnProperty('msg')) {
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
    search_text : $("#selection_ligand_bulk_search_textarea").val().trim(),
    search_type : $("#selection_ligand_bulk_search_search_type").val(),
    stereochemistry : $("#selection_ligand_bulk_search_stereochemistry").is(":checked"),
    // limit : $("#ligand_structural_search_limit").val().trim()
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
      smiles : $("#ligand_structural_search_smiles").val().trim(),
      search_type : $("#ligand_structural_search_search_type").val(),
      input_type : $("#ligand_structural_search_input_type").val(),
      similarity_threshold : $("#ligand_structural_search_similarity_threshold").val().trim(),
      stereochemistry : $("#ligand_structural_search_stereochemistry").is(":checked"),
      // limit : $("#ligand_structural_search_limit").val().trim()
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
function submitSelectionLigandBulkSearchByName(url) {
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
    // On success go to alignment page
    window.location.href = url;
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
      showAlert("Something went wrong, please try again or contact us.", "danger");
    });


}

function submitSelectionLigandBulkSearchById(url) {
  let search_params_data = {
    search_text : $("#selection_ligand_bulk_search_by_id_textarea").val().trim(),
    search_type : 'id',
    field: $("#selection_ligand_by_id_fields").val(),
    // limit : $("#ligand_structural_search_limit").val().trim()
  }
  console.log($("#selection_ligand_by_id_fields").val());

  // set CSRF csrf_token
  $.ajaxSetup({
    headers:
    { "X-CSRFToken": csrf_token }
  });

  // Submit search parameters
  $.post("/ligand/read_input_ligand_bulk_search", search_params_data,  function (data) {
    // On success go to alignment page
    window.location.href = url;
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
      showAlert("Something went wrong, please try again or contact us.", "danger");
    });


}






function submitSelectionLigandBulkSearch(url) {
  let search_params_data = {
    search_text : $("#selection_ligand_bulk_search_textarea").val().trim(),
    search_type : $("#selection_ligand_bulk_search_search_type").val(),
    stereochemistry : $("#selection_ligand_bulk_search_stereochemistry").is(":checked"),
    // limit : $("#ligand_structural_search_limit").val().trim()
  }
  console.log($("#selection_ligand_bulk_search_search_type").val());

  // set CSRF csrf_token
  $.ajaxSetup({
    headers:
    { "X-CSRFToken": csrf_token }
  });

  // Submit search parameters
  $.post("/ligand/read_input_ligand_bulk_search", search_params_data,  function (data) {
    // On success go to alignment page
    window.location.href = url;
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
      showAlert("Something went wrong, please try again or contact us.", "danger");
    });


}

function submitLigandStructuralSearch(url) {
    let search_params_data = {
      smiles : $("#ligand_structural_search_smiles").val().trim(),
      search_type : $("#ligand_structural_search_search_type").val(),
      similarity_threshold : $("#ligand_structural_search_similarity_threshold").val().trim(),
      stereochemistry : $("#ligand_structural_search_stereochemistry").is(":checked"),
      // limit : $("#ligand_structural_search_limit").val().trim()
    }


    // set CSRF csrf_token
    $.ajaxSetup({
      headers:
      { "X-CSRFToken": csrf_token }
    });

    // Submit search parameters
    $.post("/ligand/read_input_ligand_structural_search", search_params_data,  function (data) {
      // On success go to alignment page
      window.location.href = url;
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
        showAlert("Something went wrong, please try again or contact us.", "danger");
      });


}
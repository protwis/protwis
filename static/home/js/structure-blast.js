$(document).ready(function() {
    console.log('Document ready');

    const $form = $('#structure-blast-form');
    const $spinner = $('#loading-spinner');
    const $results = $('#results');
    const $fileInput = $('#input_pdb');

    if (!$form.length) {
        console.error('Form not found');
        return;
    }

    $form.on('submit', function(e) {
        console.log('Form submitted');
        e.preventDefault();

        $spinner.show();
        $results.hide();

        if ($.fn.DataTable.isDataTable('#results_table')) {
            $('#results_table').DataTable().destroy();
        }

        $.ajax({
            timeout: 120000,  // 2 minutes
            url: $form.attr('action'),
            type: $form.attr('method'),
            data: new FormData(this),
            processData: false,
            contentType: false,
            success: function(response) {
                console.log('AJAX request successful');
                $spinner.hide();
                
                var $parsedHtml = $(response);
                var $newResultsTable = $parsedHtml.find('#results_table');
                
                if ($newResultsTable.length) {
                    $results.html($newResultsTable).show();
                    initializeDataTable();
                    $fileInput.val('');
                } else {
                    $results.html('<p>No results found.</p>').show();
                }
            },
            error: function(xhr, status, error) {
                console.error('AJAX request failed');
                console.error('Status:', status);
                console.error('Error:', error);
                console.error('Status code:', xhr.status);
                console.error('Response text:', xhr.responseText);
                
                $spinner.hide();
                $results.html('<p>An error occurred. Please try again. Error details: ' + status + ' - ' + error + '</p>').show();
            }
        });
    });

    function initializeDataTable() {
        console.log('Initializing DataTable');
        if ($('#results_table').length) {
            try {
                var table = $('#results_table').DataTable({
                    "order": [[1, "desc"], [2, "desc"], [3, "asc"]],
                    "pageLength": 20,
                    'lengthChange': false,
                    "columnDefs": [{
                        "className": "text-left",
                        "targets": [9, 10]
                    },
                    {
                        "targets": [1,2,3],  
                        "render": function(data, type, row) {
                            if (type === 'display') {
                                return parseFloat(data).toExponential(1);  // Display with one decimal in exponential notation
                            }
                            return data;  // Return the original data for sorting/filtering
                        }
                    },
                ]
                });


                yadcf.init(table, [
                    { column_number: 1, filter_type: "range_number", filter_reset_button_text:false }, //TM score
                    { column_number: 2, filter_type: "range_number", filter_reset_button_text:false }, // LDDT
                    { column_number: 3, filter_type: "range_number", filter_reset_button_text:false  }, // E value
                    { column_number: 4, filter_type: "multi_select", filter_reset_button_text:false ,select_type: "select2" }, // PDB ID
                    { column_number: 5, filter_type: "multi_select", filter_reset_button_text:false ,select_type: "select2" }, // Structure/model
                    { column_number: 6, filter_type: "multi_select", filter_reset_button_text:false ,select_type: "select2" }, 
                    { column_number: 7, filter_type: "multi_select", filter_reset_button_text:false ,select_type: "select2" }, 
                    { column_number: 8, filter_type: "multi_select", filter_reset_button_text:false ,select_type: "select2", column_data_type: "html", html_data_type: "text", filter_match_mode : "exact" },
                    { column_number: 9, filter_type: "multi_select", filter_reset_button_text:false ,select_type: "select2" }, 
                    { column_number: 10, filter_type: "multi_select", filter_reset_button_text:false ,select_type: "select2" }, 
                    { column_number: 11, filter_type: "multi_select", filter_reset_button_text:false ,select_type: "select2" }, 
                ]);

                console.log('DataTable initialized successfully');
            } catch (error) {
                console.error('Error initializing DataTable:', error);
            }
        } else {
            console.log('Results table not found');
        }
    }

    if ($.fn.dataTable.isDataTable('#results_table')) {
        $('#results_table').DataTable().destroy();
    }
    
    $('#results_table').DataTable({
        "scrollX": true,
    });

    $('.tooltip-wrapper').each(function() {
        var $trigger = $(this).find('.tooltip-trigger');
        var titleText = $trigger.attr('title');
        $trigger.removeAttr('title');
        $(this).find('.tooltiptext').text(titleText);
    });

});
/*! Norges Databtables Filter Builder
 * Verison 1.0.0 - 2024
 * Architect: Søren Norge Andreassen 
*/


// #######################################################
// #### Norges custom filtering scheme for datatables ####
// #######################################################

// ##############################################################
// ## Plugin API method to determine is a column is searchable ##
// ##############################################################

$.fn.dataTable.Api.register('column().searchable()', function() {
    var ctx = this.context[0];
    return ctx.aoColumns[this[0]].bSearchable;
});

// ##############################################################
// ##                 Create Column filters                    ## 
// ##############################################################

function CreateColumnFilters(datatable_selector,column_number, column_range, header_row, filter_type) {
    // ################################################
    // # Array to track the checks for each condition #
    // ################################################
    
    let array_check = [false,false, false, false, false];  
    
    
    // ######################################################################
    // #                   Check if the datatable is valid                  #
    // ######################################################################

    if ($.fn.dataTable.isDataTable(datatable_selector)) {
        array_check[0] = true;
    } else { 
        console.log("Why did you not pass a valid DataTable? You might have mispelled, no worries, or go read up on $.fn.dataTable.isDataTable() and try again.")
    }


    // ########################
    // # Initialize variables #
    // ########################

    if (array_check[0] == true) {
        // # Number of columns in the DataTable (minus 1 to account for index starting at 0 #
        const column_number_max = datatable_selector.columns().nodes().length;
        
        // # Numer of header rows (minus 1 as index start at 0 #
        const Header_count = datatable_selector.tables().header()[0].children.length-1
        
        // # Filter type list of valid inputs # 

        const Filter_type_list = ['Multi-select','Range_filter_float','Range_filter_year']

        
        // ######################################################################
        // # Check if column_number is an integer in the range of the DataTable #
        // ######################################################################
        
        if (Number.isInteger(column_number) && column_number >= 0 && column_number <= column_number_max) {
        
            // # If a correct value for the column number is set, return true # 
            array_check[1] = true;
        
        } else {

            // # Write to console that an incorrect column number is set (array check for column number remains false) #
            console.log(`Column number ${column_number} is not a valid integer within the range of 0 to ${column_number_max}.`);

        }

        // ######################################################################
        // #  Check if header_row is a valid integer and within the datatable   #
        // ######################################################################

        if (Number.isInteger(header_row) && header_row >= 0 && header_row <= Header_count) {
            array_check[2] = true;

        } else {

            console.log(`Header row ${header_row} is not a valid integer. Or selected header is not between 0 and max header count ${Header_count}`);

        }
        // ######################################################################
        // #                Check if filter type is a valid inpunt              #
        // ######################################################################

        if (Filter_type_list.includes(filter_type)) {

            array_check[3] = true;

        } else {

            console.log(`Filter type (${filter_type}) is not included in current valid filter types: ${Filter_type_list}`)
        }

        // ######################################################################
        // #                Check if column_range is an integer                 # 
        // #            that does not exceed the DataTable's max number         #
        // #          of columns from the starting point of column_number       #
        // ######################################################################

        if (Number.isInteger(column_range) && column_range >= 1 && column_range <= column_number_max && (column_range + column_number) <= column_number_max) {
            array_check[4] = true;
            
        } else {
            array_check[1] = true;
            console.log(`Column range ${column_range} is not a valid integer within the range of 0 to ${column_number_max} from the starting point of column number ${column_number}. As it will end at index ${column_number+column_range}`);
        }
        
        
        
        // Check if all conditions are met
        const allConditionsMet = array_check.every(Boolean);
        
        // Construct and return the filter array
        if (allConditionsMet == true) {
            if (column_range > 1) {
                let Multiple_FilterArray = []
                for (i = 0; i < column_range; i++) {
                    if (i == 0) {
                    const FilterArray = [[column_number,1,header_row,filter_type]]
                    Multiple_FilterArray = Multiple_FilterArray.concat(FilterArray)
                    } else {
                        FilterArray = [[column_number+i,1,header_row,filter_type]]
                        Multiple_FilterArray = Multiple_FilterArray.concat(FilterArray)
                    }
                }
                return Multiple_FilterArray
            } else {
            const filterArray = [[column_number, column_range, header_row, filter_type]];
            return filterArray
            }
        } else {
            console.log(`will not pass filter array of: [DataTable_assigned_variable, ${column_number},${column_range},${header_row},${filter_type}]`)
        }
        
    }
}

var prefixSorter = function(results) {
    if (!results || results.length == 0)
      return results
  
    // Find the open select2 search field and get its value
    var term = document.querySelector('.select2-search__field').value.toLowerCase()
    if (term.length == 0)
      return results
  
    return results.sort(function(a, b) {
        const key = term.toLowerCase();
        const isGoodMatchA = a.toLowerCase().startsWith(key);
        const isGoodMatchB = b.toLowerCase().startsWith(key);
    
        if (isGoodMatchA ^ isGoodMatchB) { // XOR
            return isGoodMatchA ? -1 : 1;
        }
    
        return a.localeCompare(b);
    });
  }


// ##############################################################
// ####    Norges Datatables Filter Builder - Lets go!       ####
// ##############################################################

function createDropdownFilters(api,column_filters) {
    // ###################################
    // ## check if the input is correct ##
    // ###################################
    
    // ################################################
    // # Check if the filter elements are all correct #
    // ################################################
    filter_pass = false;
    if (column_filters.length > 0 && column_filters.every(subArray => Array.isArray(subArray) && subArray.length > 0)) {
        if (!column_filters.map(subArray => subArray[0]).some((value, index, array) => array.indexOf(value) !== index)) {
            filter_pass = true;
        } else {
            console.log("Your filters have overlapping columns! Write your code better.. plz..")
        }
    } else {
        console.log("oh lordy lord, your filters are not valid for the filter initialization, please go through the filters and see if you can correct your mistakes.")
    }

    // ######################################################################
    // #                   Check if the datatable is valid                  #
    // ######################################################################

    if ($.fn.dataTable.isDataTable(api)) {
        DT_pass = true;
    } else { 
        console.log("Why did you not pass a valid DataTable? You might have mispelled, no worries, or go read up on $.fn.dataTable.isDataTable() and try again.")
        DT_pass = false;
    }

    // #########################################################
    // # All parameters are passed and the filters are applied #
    // #########################################################

    if (filter_pass == true && DT_pass == true) {
        
        // ####################################################
        // ###     Iterate over all filters and apply       ###
        // ####################################################

        for (const column_filter of column_filters) {
            
            // ###################################
            // ###     Initialize values       ###
            // ###################################
            
            const column_number = column_filter[0];
            const header_row = column_filter[2];
            const filter_type = column_filter[3];

            // #####################################################
            // # Dependent on the filter type -> apply appropriate #
            // # currently avaliable: Multi-select, 
            // # Range_filter_float, and Range_filter_year.
            // #####################################################

            // ############################################
            // ##          Multi select filter           ##
            // ############################################

            if (filter_type == 'Multi-select') {
                api.columns(column_number).every(function() {
                if (this.searchable()) {
                    var that = this;
                    var col = this.index();
                    // Only create if not there or blank
                    var selected = $('thead tr:eq('+header_row+') td:eq(' + col + ') select').val();
                    if (selected === undefined || selected === '') {
                        // // Create the `select` element
                        // $('thead tr:eq('+header_row+') td')
                        //     .eq(col)
                        //     .empty();
                        var select = $('<select id="bob'+col+'" class="select2" style="width: 80%;"></select>')
                            .appendTo($('thead tr:eq('+header_row+') td').eq(col))
                            .on('change', function() {
                            //Get the "text" property from each selected data 
                            //regex escape the value and store in array
                            var data = $.map( $(this).select2('data'), function( value, key ) {
                                // return value.text ? '^' + $.fn.dataTable.util.escapeRegex(value.text) + '$' : null; // exact match
                                return value.text ? $.fn.dataTable.util.escapeRegex(value.text): null; // not exact match --> string in string
                            });
                
                            //if no data selected use ""
                            if (data.length === 0) {
                                data = [""];
                            }
                            
                            //join array into string with regex or (|)
                            var val = data.join('|');
                    
                            //search for the option(s) selected
                            that.search(val ? val : '', true, false ).draw();
                        } );
                        
                
                        select.append('<option>'+''+'</option>'); // <-- this is added to be selected as the first option, so it doesnt select a value.      
                        api.cells(null, col, {
                            search: 'applied'
                        }).data().unique().sort().each(function(d) {
                            select.append($('<option>' + d + '</option>'));
                        });

                        // Remove the empty string value from the selection options -- This is the workaround to get fast loading without selection of the first item.
                        $('#bob'+col).attr('multiple', 'multiple');
                        $('#bob'+col + " option")[0].remove();
                        

                        var select2 = $('#bob'+col).select2({
                            // theme: 'bootstrap-4',
                            multiple: true,
                            closeOnSelect: true,
                            placeholder: {text: "Filter"},
                            dropdownAutoWidth : true,
                            tags: true, // allows for selection of undefined values --> Needs the not exact match to function well.
                            // sorter: data => data.sort((a, b) => a.text.length < b.text.length ? -1 : a.text.length > b.text.length ? 1 : 0)
                    });
                    //initially clear select otherwise first option is selected
                    // $('.select2').val(null).trigger('change');
                        }
                }
            }); // End of muli-selct filter

            // ############################################
            // ##          Range filter float            ##
            // ############################################

            } else if ( filter_type == "Range_filter_float" ){
                api.columns([column_number]).every(function() {
                if (this.searchable()) {
                    var that = this;
                    var col = column_number;
                    // Only create if not there or blank
                    var selected = $('thead tr:eq('+header_row+') td:eq(' + col + ') select').val();
                    if (selected === undefined || selected === '') {
                                // Create the `select` element
                                $('thead tr:eq('+header_row+') td')
                                    .eq(col)
                                    .empty();
                        var select_min = $('<input id="bob'+col+'min" class="pull-left select2" placeholder="min" style="width: 45%;text-align:center;"></input>')
                                    .appendTo($('thead tr:eq('+header_row+') td').eq(col)).on('keyup', function() {
                            // Custom filtering function which will search data in column four between two values
                            $.fn.dataTable.ext.search.push(
                                function( settings, data, dataIndex ) {
                                    var min = parseInt( $('#bob'+col+'min').val(), 10 );
                                    var max = parseInt( $('#bob'+col+'max').val(), 10 );
                                    var age = parseFloat( data[col] ) || 0; // use data float data
                            
                                    if ( ( isNaN( min ) && isNaN( max ) ) ||
                                        ( isNaN( min ) && age <= max ) ||
                                        ( min <= age   && isNaN( max ) ) ||
                                        ( min <= age   && age <= max ) )
                                    {
                                        return true;
                                    }
                                    return false;
                                }
                            
                            );
                        that.draw();
                        });
                        var select_max = $('<input id="bob'+col+'max" class="pull-right select2" placeholder="max" style="width: 45%;text-align:center;"></input>')
                        .appendTo($('thead tr:eq('+header_row+') td').eq(col)).on('keyup', function() {
                        $.fn.dataTable.ext.search.push(
                                function( settings, data, dataIndex ) {
                                    var min = parseInt( $('#bob'+col+'min').val(), 10 );
                                    var max = parseInt( $('#bob'+col+'max').val(), 10 );
                                    var age = parseFloat( data[col] ) || 0;
                            
                                    if ( ( isNaN( min ) && isNaN( max ) ) ||
                                        ( isNaN( min ) && age <= max ) ||
                                        ( min <= age   && isNaN( max ) ) ||
                                        ( min <= age   && age <= max ) )
                                    {
                                        return true;
                                    }
                                    return false;
                                }
                            );
                        that.draw();
                        });
                    }
                }
            }); // End of Range filter float 
            
            // ############################################
            // ##          Range filter Year             ##
            // ############################################
            
            } else if ( filter_type == "Range_filter_year" ) {
                api.columns([column_number]).every(function() {
                if (this.searchable()) {
                    var that = this;
                    var col = column_number;
                    // Only create if not there or blank
                    var selected = $('thead tr:eq('+header_row+') td:eq(' + col + ') select').val();
                    if (selected === undefined || selected === '') {
                        // Create the `select` element
                        $('thead tr:eq('+header_row+') td').eq(col).empty();  
                        var select_min = $('<select id="bob'+col+'min" class="pull-left select2" style="width: 45%"></select>')
                            .appendTo($('thead tr:eq('+header_row+') td').eq(col)).on('change', function() {
                                // Custom filtering function which will search data in column four between two values
                                $.fn.dataTable.ext.search.push(
                                    function( settings, data, dataIndex ) {
                                        var min = parseInt( $('#bob'+col+'min').val(), 10 );
                                        var max = parseInt( $('#bob'+col+'max').val(), 10 );
                                        var age = parseFloat( data[col] ) || 0; // use data for the age column
                                
                                        if ( ( isNaN( min ) && isNaN( max ) ) ||
                                            ( isNaN( min ) && age <= max ) ||
                                            ( min <= age   && isNaN( max ) ) ||
                                            ( min <= age   && age <= max ) )
                                        {
                                            return true;
                                        }
                                        return false;
                                    }
                                
                                );
                            that.draw();
                        }); // select_min variable end
                        var select_max = $('<select id="bob'+col+'max" class="pull-right select2" style="width: 45%;"></select>')
                            .appendTo($('thead tr:eq('+header_row+') td').eq(col)).on('change', function() {
                                $.fn.dataTable.ext.search.push(
                                        function( settings, data, dataIndex ) {
                                            var min = parseInt( $('#bob'+col+'min').val(), 10 );
                                            var max = parseInt( $('#bob'+col+'max').val(), 10 );
                                            var age = parseFloat( data[col] ) || 0; // use data for the age column
                                    
                                            if ( ( isNaN( min ) && isNaN( max ) ) ||
                                                ( isNaN( min ) && age <= max ) ||
                                                ( min <= age   && isNaN( max ) ) ||
                                                ( min <= age   && age <= max ) )
                                            {
                                                return true;
                                            }
                                            return false;
                                        }
                                );
                            that.draw();
                        }); // select_max variable end
                    select_min.append('<option>'+''+'</option>'); // <-- this is added to be selected as the first option, so it doesnt select a value.
                    select_max.append('<option>'+''+'</option>'); // <-- this is added to be selected as the first option, so it doesnt select a value.
                    api.cells(null, col, {
                        search: 'applied'
                    }).data().unique().sort().each(function(d) {
                        select_min.append($('<option>' + d + '</option>'));
                        select_max.append($('<option>' + d + '</option>'));
                    });

                    // // Remove the empty string value from the selection options -- This is the workaround to get fast loading without selection of the first item.
                    // $('#bob'+col+'max').attr('multiple', 'multiple');
                    // $('#bob'+col+'max' + " option")[0].remove();
                    // $('#bob'+col+'min').attr('multiple', 'multiple');
                    // $('#bob'+col+'min' + " option")[0].remove();

                    var select2_min = $('#bob'+col+'min').select2({
                            multiple: false,
                            closeOnSelect: true,
                            dropdownAutoWidth : true,
                            placeholder: {text: "min"}
                            // width: 'resolve'
                    });
                    var select2_max = $('#bob'+col+'max').select2({
                        multiple: false,
                        closeOnSelect: true,
                        dropdownAutoWidth : true,
                        placeholder: {text: "max"}
                        // width: 'resolve'
                    });
                    //initially clear select otherwise first option is selected
                    // $('.select2').val(null).trigger('change');
                }
                }
            }); // End of Range filter year
        }
        if (column_number == 0) {
            const time = Date.now();
            const time_now = time-timerStart;
            var time_s = (time_now*0.001).toFixed(2);
            console.log("col "+column_number+": "+time_s+"s")
        } else {
            const current_time = time_s
            const time = Date.now();
            const time_now = time-timerStart;
            var time_s = time_now*0.001;
            interval_time = (time_s-current_time).toFixed(2);
            console.log("col "+column_number+": "+interval_time+"s");
        }
    }
}
// $('.select2').val(null).trigger('change');
}


// return data.sort(function (a, b) {
//     if (a.text.toLowerCase() === b.text.toLowerCase()) {
//         return -1; // exact match, prioritize a
//     } else if (a.text > b.text) {
//         return 1; // b should come before a alphabetically
//     } else {
//         return -1; // a should come before b alphabetically
//     }
// });
// }
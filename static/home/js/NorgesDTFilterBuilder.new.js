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

function CreateColumnFilters(datatable_selector,column_number, column_range, filter_type, width) {
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

    if (array_check[0] === true) {

        // ####################################################################################
        // ## Number of columns in the DataTable (minus 1 to account for index starting at 0 ##
        // ####################################################################################
        
        const column_number_max = datatable_selector.columns().nodes().length;
        
        // ######################################
        // ## Filter type list of valid inputs ##
        // ######################################

        const Filter_type_list = ['Multi-select-exact','Multi-select-unspecific','Range-float-vertical','Range-float-horizontal','Range-select-vertical','Range-select-horizontal']
        
        // ######################################################################
        // # Check if column_number is an integer in the range of the DataTable #
        // ######################################################################

        if (Number.isInteger(column_number) && column_number >= 0 && column_number <= column_number_max) {
            
            // #################################################################
            // ## If a correct value for the column number is set, return true # 
            // #################################################################

            array_check[1] = true;

        } else {

            // # Write to console that an incorrect column number is set (array check for column number remains false) #
            console.log(`Column number ${column_number} is not a valid integer within the range of 0 to ${column_number_max}.`);

        }

        // ######################################################################
        // #                Check if filter type is a valid inpunt              #
        // ######################################################################

        if (Filter_type_list.includes(filter_type)) {

            array_check[2] = true;

        } else {

            console.log(`Filter type (${filter_type}) is not included in current valid filter types: ${Filter_type_list}`)
        }

        // ######################################################################
        // #                Check if column_range is an integer                 #
        // #            that does not exceed the DataTable's max number         #
        // #          of columns from the starting point of column_number       #
        // ######################################################################

        if (Number.isInteger(column_range) && column_range >= 1 && column_range <= column_number_max && (column_range + column_number) <= column_number_max) {
            array_check[3] = true;
            
        } else {
            console.log(`Column range ${column_range} is not a valid integer within the range of 0 to ${column_number_max} from the starting point of column number ${column_number}. As it will end at index ${column_number+column_range}`);
        }
        
        
         // #############################################
        // #  Check if width is set and correct value   #
        // ##############################################

        if (width === undefined) {

            array_check[4] = true;

        } else {
            // width
            widthsplit = width.match(/-?\d+|[a-zA-Z]+|[^a-zA-Z\d]+/g);
            if (widthsplit.length === 2) {
                if (/^-?\d+$/.test(widthsplit[0]) && parseInt(widthsplit[0], 10) > 0) {
                    if (widthsplit[1] === "px" || widthsplit[1] === "%") {
                        array_check[4] = true;
                    } else {
                        console.log(`${width} (1) is not correctly formatted: needs to be [integer][px/%], e.g 100px or 50%`);
                    }
                } else {
                    console.log(`${width} (1) is not correctly formatted: needs to be [integer][px/%], e.g 100px or 50%`);
                }
            } else {
            console.log(`${width} (2) is not correctly formatted: needs to be [integer][px/%], e.g 100px or 50%`);
        }
        }
        // Check if all conditions are met
        const allConditionsMet = array_check.every(Boolean);

        // Construct and return the filter array
        if (allConditionsMet === true) {
            if (column_range > 1) {
                let Multiple_FilterArray = []
                for (i = 0; i < column_range; i++) {
                    if (i === 0) {
                    const FilterArray = [[column_number,1,filter_type]]
                    Multiple_FilterArray = Multiple_FilterArray.concat(FilterArray)
                    } else {
                        FilterArray = [[column_number+i,1,filter_type]]
                        Multiple_FilterArray = Multiple_FilterArray.concat(FilterArray)
                    }
                }
                return Multiple_FilterArray
            } else {
            const filterArray = [[column_number, column_range, filter_type]];
            return filterArray
            }
        } else {
            console.log(`will not pass filter array of: [${"dt assigned value"},${column_number},${column_range},${filter_type}]`)
        }

    }
}

// ##############################################################
// ####    Norges Datatables Filter Builder - Lets go!       ####
// ##############################################################

function createDropdownFilters(api,column_filters) {

    // ###################################
    // ## check if the input is correct ##
    // ###################################

    filter_pass = false;
    if (column_filters.length > 0 && column_filters.every(subArray => Array.isArray(subArray) && subArray.length > 0)) {
        if (!column_filters.map(subArray => subArray[0]).some((value, index, array) => array.indexOf(value) !== index)) {
            filter_pass = true;
        } else {
            console.log("Your filters have overlapping columns! Write your code better.. plz..");
        }
    } else {
        console.log("oh lordy lord, your filters are not valid for the filter initialization, please go through the filters and see if you can correct your mistakes.");
    }

    DT_pass = false;
    if ($.fn.dataTable.isDataTable(api)) {
        DT_pass = true;
    } else { 
        console.log("Why did you not pass a valid DataTable? You might have mispelled, no worries, or go read up on $.fn.dataTable.isDataTable() and try again.");
    }

    // #########################################################
    // # All parameters are passed and the filters are applied #
    // #########################################################

    if (filter_pass === true && DT_pass === true) {
        
        const Table_id = api.table().node().id
        const visibility_state = api.columns().visible().toArray();
        api.columns().visible(true); // Temporarily make all columns visible for header access

        var $tableHeader = $(api.table().header());
        var $headerContainer = $(api.table().node()).closest('.dataTables_wrapper').find('.dt-scroll-headInner');
        if ($headerContainer.length) {
            $tableHeader = $headerContainer.find('table thead');
        }
        const headerRows = $tableHeader.find('tr');

        for (const column_filter of column_filters) {
            const column_number_from_config = column_filter[0]; // This is the config column number
            const filter_type_from_config = column_filter[2];

            const dtColumn = api.column(column_number_from_config); // Get DT API for this column

            if (filter_type_from_config === 'Multi-select-exact' || filter_type_from_config === 'Multi-select-unspecific') {
                if (dtColumn.searchable()) {
                    
                    var selected_cell_jq;
                    // The `column_number_from_config` is the index for DataTables columns.
                    // The `headerRows[2]` is the third TR element.
                    // `eq(column_number_from_config)` targets the TD at that column index within the third TR.
                    if (headerRows.length >= 3) { 
                        selected_cell_jq = $(headerRows[2]).find('td').eq(column_number_from_config);
                    } else {
                        console.error(`Table ${Table_id} does not have the expected 3 header rows for filter placement in column ${column_number_from_config}.`);
                        continue; 
                    }

                    if (!selected_cell_jq || selected_cell_jq.length === 0) {
                        console.error(`Could not find filter cell for column ${column_number_from_config} in table ${Table_id}.`);
                        continue;
                    }
                    
                    // Use these captured values in the event handler and for element IDs
                    var currentColumnAPI_for_multi = dtColumn; 
                    var currentColIndex_for_multi = dtColumn.index(); 

                    var select_html = '<select id="' + Table_id + '_Filter' + currentColIndex_for_multi + '" class="select2" style="width: 80%;"></select>';
                    selected_cell_jq.html(select_html);

                    $('#' + Table_id + '_Filter' + currentColIndex_for_multi).on('change', (function(capturedColumnInstance, capturedColumnIndex, capturedFilterType) { 
                        return function() { 
                            var data = $.map($(this).select2('data'), function(value, key) {
                                if (capturedFilterType === 'Multi-select-exact') { 
                                    return value.text ? '^' + $.fn.dataTable.util.escapeRegex(value.text) + '$' : null;
                                } else if (capturedFilterType === 'Multi-select-unspecific') {
                                    return value.text ? $.fn.dataTable.util.escapeRegex(value.text) : null;
                                }
                            });

                            if (data.length === 0) data = [""];
                            var val = data.join('|');
                            console.log(`Col ${capturedColumnIndex} - Search Value (val):`, val); 
                            capturedColumnInstance.search(val ? val : '', true, false).draw(); 
                        };
                    })(currentColumnAPI_for_multi, currentColIndex_for_multi, filter_type_from_config)); // Pass filter_type too if it's needed inside


                    $('#' + Table_id + '_Filter' + currentColIndex_for_multi).append('<option>' + '' + '</option>');
                    
                    var renderFunction_multi = currentColumnAPI_for_multi.settings()[0].aoColumns[currentColIndex_for_multi].mRender; 

                    currentColumnAPI_for_multi.data().unique().sort().each(function(d_multi) { 
                         var renderedValue_multi = renderFunction_multi ? renderFunction_multi(d_multi, 'display', undefined, {col: currentColIndex_for_multi, row: -1, settings: currentColumnAPI_for_multi.settings()[0]}) : d_multi;
                         renderedValue_multi = (typeof renderedValue_multi === 'string') ? renderedValue_multi : String(renderedValue_multi);
                        var tempDiv_multi = document.createElement("div");
                        tempDiv_multi.innerHTML = renderedValue_multi;
                        var textContent_multi = tempDiv_multi.textContent || tempDiv_multi.innerText || "";
                        // console.log(`Col ${currentColIndex_for_multi} - Original Data (d):`, d_multi, `Rendered for Display:`, renderedValue_multi, `Text Content for Option:`, textContent_multi);
                        $('#' + Table_id + '_Filter' + currentColIndex_for_multi).append($('<option>' + textContent_multi + '</option>'));
                    });

                    $('#' + Table_id + '_Filter' + currentColIndex_for_multi).attr('multiple', 'multiple');
                    $('#' + Table_id + '_Filter' + currentColIndex_for_multi + " option").first().remove();

                    if (filter_type_from_config === 'Multi-select-unspecific') {
                         $('#' + Table_id + '_Filter' + currentColIndex_for_multi).select2({
                            multiple: true,
                            closeOnSelect: true, 
                            placeholder: "Filter", 
                            dropdownAutoWidth: true,
                            tags: true,
                            allowClear: true 
                        });
                    } else if (filter_type_from_config === 'Multi-select-exact') {
                         $('#' + Table_id + '_Filter' + currentColIndex_for_multi).select2({
                            multiple: true,
                            closeOnSelect: true, 
                            placeholder: "Filter", 
                            dropdownAutoWidth: true,
                            allowClear: true 
                        });
                    }
                }
            }
            else if (filter_type_from_config === "Range-float-vertical" || filter_type_from_config === "Range-float-horizontal") {
                if (dtColumn.searchable()) {
                    var currentColIndex_Range = dtColumn.index(); 

                    var selected_cell_jq_range;
                     if (headerRows.length >= 3) {
                        selected_cell_jq_range = $(headerRows[2]).find('td').eq(currentColIndex_Range);
                    } else {
                        console.error(`Table ${Table_id} does not have the expected 3 header rows for range filter in column ${currentColIndex_Range}.`);
                        continue;
                    }
                    if (!selected_cell_jq_range || selected_cell_jq_range.length === 0) {
                        console.error(`Could not find filter cell for column ${currentColIndex_Range} in table ${Table_id} for range filter.`);
                        continue;
                    }

                    var html_input1 = '<input id="' + Table_id + '_Filter' + currentColIndex_Range + 'min" class="select2" placeholder="min" style="width: 35px;text-align:center;margin-top: 5px;"></input>';
                    var html_input2 = '<input id="' + Table_id + '_Filter' + currentColIndex_Range + 'max" class="select2" placeholder="max" style="width: 35px;text-align:center;"></input>';
                    
                    if (filter_type_from_config === "Range-float-vertical") {
                        selected_cell_jq_range.html(html_input2 + '<br>' + html_input1);
                    } else if (filter_type_from_config === "Range-float-horizontal") {
                        selected_cell_jq_range.html(html_input1 + html_input2);
                    }

                    $('#' + Table_id + '_Filter' + currentColIndex_Range + 'min' + ',' + '#' + Table_id + '_Filter' + currentColIndex_Range + 'max').keyup( (function(capturedColumnIndexForRange) {
                        return function() {
                            for (var i = $.fn.dataTable.ext.search.length - 1; i >= 0; i--) {
                                if ($.fn.dataTable.ext.search[i].name === Table_id + '_range_filter_' + capturedColumnIndexForRange) {
                                    $.fn.dataTable.ext.search.splice(i, 1);
                                }
                            }

                            var range_filter_function = function(settings, data, dataIndex) {
                                if (settings.nTable.id !== Table_id) return true;

                                var min_val_str = $('#' + Table_id + '_Filter' + capturedColumnIndexForRange + 'min').val();
                                var max_val_str = $('#' + Table_id + '_Filter' + capturedColumnIndexForRange + 'max').val();
                                var min = min_val_str !== "" ? parseFloat(min_val_str) : NaN;
                                var max = max_val_str !== "" ? parseFloat(max_val_str) : NaN;
                                var data_col_val = data[capturedColumnIndexForRange];
                                
                                var data_range;
                                if (typeof data_col_val === 'string') {
                                    var clean_data = data_col_val.replace(/<[^>]*>/g, "").trim();
                                    if (clean_data === "-" || clean_data === "None" || clean_data === "" || clean_data.toLowerCase() === "nan" ) {
                                        data_range = NaN; 
                                    } else {
                                        data_range = parseFloat(clean_data);
                                    }
                                } else if (data_col_val === null || data_col_val === undefined) {
                                    data_range = NaN;
                                } else {
                                    data_range = parseFloat(data_col_val);
                                }
                                
                                // If both min and max are not set (NaN), then all rows should pass, regardless of data_range
                                if (isNaN(min) && isNaN(max)) {
                                    return true;
                                }

                                // If data_range could not be parsed as a number, it cannot match a numeric range
                                // (unless both min and max are also NaN, which is handled above)
                                if (isNaN(data_range)) {
                                    return false; 
                                }

                                // Now, data_range is a valid number, proceed with standard range checks
                                if ((isNaN(min) && data_range <= max) ||   // Only max is set
                                    (min <= data_range && isNaN(max)) ||   // Only min is set
                                    (min <= data_range && data_range <= max)) { // Both min and max are set
                                    return true;
                                }
                                return false;
                            };
                            range_filter_function.name = Table_id + '_range_filter_' + capturedColumnIndexForRange;
                            $.fn.dataTable.ext.search.push(range_filter_function);
                            api.draw(); 
                        };
                    })(currentColIndex_Range) ); 
                }
            }
            // NOTE: Range-select-vertical and Range-select-horizontal filter types are not yet updated with explicit closure for event handlers.
            // If these are used and show similar issues, they will need the same IIFE pattern for their event handlers.
            else if ( filter_type_from_config === "Range-select-vertical" || filter_type_from_config === "Range-select-horizontal") {
                api.columns([column_number_from_config]).every(function() { // Still using 'this' context here, might need adjustment if problematic
                if (this.searchable()) {
                    var that_range_select = this; // Potential closure issue if 'this' is not what's expected in handler
                    var col_range_select = column_number_from_config;
                    
                    var selected_cell_jq_range_select;
                    if (headerRows.length >= 3) {
                        selected_cell_jq_range_select = $(headerRows[2]).find('td').eq(col_range_select);
                    } else {
                         console.error(`Table ${Table_id} does not have the expected 3 header rows for range-select filter in column ${col_range_select}.`);
                        return; // continue to next column_filter in the outer loop
                    }
                     if (!selected_cell_jq_range_select || selected_cell_jq_range_select.length === 0) {
                        console.error(`Could not find filter cell for column ${col_range_select} in table ${Table_id} for range-select filter.`);
                        return;
                    }

                    var html_input1_rs, html_input2_rs;
                    if (filter_type_from_config === "Range-select-vertical") {
                        html_input1_rs = '<select id="'+Table_id+'_Filter'+col_range_select+'min" class="select2" data-placeholder="Min" style="width: 75%"></select>'
                        html_input2_rs = '<select id="'+Table_id+'_Filter'+col_range_select+'max" class="select2" data-placeholder="Max" style="width: 75%;"></select>'
                        selected_cell_jq_range_select.html(html_input1_rs+'<br>'+html_input2_rs);
                    } else if (filter_type_from_config === "Range-select-horizontal") {
                        html_input1_rs = '<select id="'+Table_id+'_Filter'+col_range_select+'min" class="pull-left select2" style="width: 45%"></select>'
                        html_input2_rs = '<select id="'+Table_id+'_Filter'+col_range_select+'max" class="pull-right select2" style="width: 45%;"></select>'
                        selected_cell_jq_range_select.html(html_input1_rs+html_input2_rs);
                    }

                    // Apply IIFE for closure on col_range_select and that_range_select (DataTables column API instance)
                    $('#'+Table_id+'_Filter'+col_range_select+'min'+','+'#'+Table_id+'_Filter'+col_range_select+'max').on('select2:select', 
                        (function(capturedCol_rs, capturedThat_rs) {
                            return function() {
                                for (var i = $.fn.dataTable.ext.search.length - 1; i >= 0; i--) {
                                    if ($.fn.dataTable.ext.search[i].name === Table_id + '_range_select_filter_' + capturedCol_rs) {
                                        $.fn.dataTable.ext.search.splice(i, 1);
                                    }
                                }
                                var range_select_filter_func = function( settings, data, dataIndex ) {
                                    if ( settings.nTable.id !== Table_id ) return true;
                                    
                                    var min = parseFloat( $('#'+Table_id+'_Filter'+capturedCol_rs+'min').val(), 10 );
                                    var max = parseFloat( $('#'+Table_id+'_Filter'+capturedCol_rs+'max').val(), 10 );
                                    var data_col_val_rs = data[capturedCol_rs];
                                    var data_range_rs;

                                    if (typeof data_col_val_rs === 'string') {
                                        var clean_data_rs = data_col_val_rs.replace(/<[^>]*>/g, "").trim();
                                        if (clean_data_rs === "-" || clean_data_rs === "None" || clean_data_rs === "" || clean_data_rs.toLowerCase() === "nan") {
                                            data_range_rs = NaN;
                                        } else {
                                            data_range_rs = parseFloat(clean_data_rs);
                                        }
                                    } else if (data_col_val_rs === null || data_col_val_rs === undefined) {
                                        data_range_rs = NaN;
                                    } else {
                                        data_range_rs = parseFloat(data_col_val_rs);
                                    }

                                    if (isNaN(data_range_rs)) return false;

                                    if ( ( isNaN( min ) && isNaN( max ) ) ||
                                         ( isNaN( min ) && data_range_rs <= max ) ||
                                         ( min <= data_range_rs   && isNaN( max ) ) ||
                                         ( min <= data_range_rs   && data_range_rs <= max ) )
                                    {
                                        return true;
                                    }
                                    return false;
                                };
                                range_select_filter_func.name = Table_id + '_range_select_filter_' + capturedCol_rs;
                                $.fn.dataTable.ext.search.push(range_select_filter_func);
                                api.draw(); // Draw main table
                            };
                        })(col_range_select, that_range_select) // Pass current col_range_select and that_range_select
                    );

                    $('#'+Table_id+'_Filter'+col_range_select+'min').append('<option></option>'); 
                    $('#'+Table_id+'_Filter'+col_range_select+'max').append('<option></option>'); 

                    // Use the captured 'that_range_select' (which is dtColumn for this iteration)
                    that_range_select.data().unique().sort().each(function(d_rs) {
                        // It's better to get render function directly from settings if 'that_range_select' is the column API
                        var renderFunction_rs = that_range_select.settings()[0].aoColumns[col_range_select].mRender;
                        var renderedValue_rs = renderFunction_rs ? renderFunction_rs(d_rs, 'display', undefined, {col: col_range_select, row: -1, settings: that_range_select.settings()[0]}) : d_rs;
                        renderedValue_rs = (typeof renderedValue_rs === 'string') ? renderedValue_rs : String(renderedValue_rs);
                        var tempDiv_rs = document.createElement("div");
                        tempDiv_rs.innerHTML = renderedValue_rs;
                        var textContent_rs = tempDiv_rs.textContent || tempDiv_rs.innerText || "";

                        $('#'+Table_id+'_Filter'+col_range_select+'min').append($('<option>' + textContent_rs + '</option>'));
                        $('#'+Table_id+'_Filter'+col_range_select+'max').append($('<option>' + textContent_rs + '</option>'));
                    });

                    $('#'+Table_id+'_Filter'+col_range_select+'min').select2({
                            minimumResultsForSearch: -1, multiple: false, closeOnSelect: true, dropdownAutoWidth : true,
                            placeholder: "Min" // Directly set placeholder
                    });
                    $('#'+Table_id+'_Filter'+col_range_select+'max').select2({
                        minimumResultsForSearch: -1, multiple: false, closeOnSelect: true, dropdownAutoWidth : true,
                        placeholder: "Max" // Directly set placeholder
                    });
                }
            }); 
        }
    }

        // Restore original column visibility
        api.columns().visible(false); 
        for (var j = 0; j < visibility_state.length; ++j) {
             if (visibility_state[j]) { 
                api.column(j).visible(true, false); 
            }
        }
        api.columns.adjust().draw(false); 
    }
}
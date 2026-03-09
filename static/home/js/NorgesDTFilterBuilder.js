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

        const Filter_type_list = ['Multi-select-exact','Multi-select-unspecific','Range-float-vertical','Range-float-horizontal','Range-select-vertical','Range-select-horizontal','Multi-select-exact-filter']
        
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

    // ################################################
    // # Check if the filter elements are all correct #
    // ################################################

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

    // ######################################################################
    // #                   Check if the datatable is valid                  #
    // ######################################################################

    if ($.fn.dataTable.isDataTable(api)) {

        DT_pass = true;

    } else { 

        console.log("Why did you not pass a valid DataTable? You might have mispelled, no worries, or go read up on $.fn.dataTable.isDataTable() and try again.");
        
        DT_pass = false;

    }

    // #########################################################
    // # All parameters are passed and the filters are applied #
    // #########################################################

    if (filter_pass === true && DT_pass === true) {
        
        // #####################
        // ## Global variable ##
        // #####################

        const Table_id = api.table().node().id

        // ###########################
        // ## Handle hidden columns ##
        // ###########################

        const visibility_state = api.columns().visible().toArray();

        // ###############################################################################################
        // ##  Set all columns visibility to true and change them back at the end of the Filter builder ##
        // ## ############################################################################################
        
        api.columns().visible(true);

        // ####################################################
        // ###     Iterate over all filters and apply       ###
        // ####################################################

        for (const column_filter of column_filters) {

            // ###################################
            // ###     Initialize values       ###
            // ###################################

            const column_number = column_filter[0];

            const filter_type = column_filter[2];
            
            
            // #####################################################
            // # Dependent on the filter type -> apply appropriate #
            // # currently avaliable:                              #
            // # Multi-select-exact / Multi-select-unspecific      #
            // # Range-float-vertical / Range-float-horizontal     #
            // # Range_filter_select                               #
            // #####################################################

            // ############################################
            // ##          Multi select filter           ##
            // ############################################

            if (filter_type === 'Multi-select-exact' || filter_type === 'Multi-select-unspecific') {
                api.columns([column_number]).every(function() {
                if (this.searchable()) {
                    var that = this;
                    var col = this.index();
                    // Only create if not there or blank
                    var selected_cell = that.column(col).header();
                    var select = '<select id="'+Table_id+'_Filter'+col+'" class="select2" style="width: 80%;"></select>';
                    selected_cell.innerHTML = select;
                    $('#'+Table_id+'_Filter'+col).on('change', function() {

                        // #####################################################
                        // ## Get the "text" property from each selected data ##
                        // ##    regex escape the value and store in array    ##
                        // #####################################################

                        var data = $.map( $(this).select2('data'), function( value, key ) {
                            if (filter_type === 'Multi-select-exact') {
                                return value.text ? '^' + $.fn.dataTable.util.escapeRegex(value.text) + '$' : null; // exact match
                            } else if (filter_type === 'Multi-select-unspecific') {
                                return value.text ? $.fn.dataTable.util.escapeRegex(value.text): null; // not exact match --> string in string
                            }
                        });

                        // ################################
                        // ## if no data selected use "" ##
                        // ################################
                        if (data.length === 0) {
                            data = [""];
                        }

                        // ##############################################
                        // ## join array into string with regex or (|) ##
                        // ##############################################
                        var val = data.join('|');

                        // #######################################
                        // ## search for the option(s) selected ##
                        // #######################################
                        that.search(val ? val : '', true, false ).draw();
                    }); // End of search function (on change)
                    
                    // ##  add empty to be selected as the first option ##
                    // ##  so it doesnt select a value. (process speed optimazation) ##
                    $('#'+Table_id+'_Filter'+col).append('<option>'+''+'</option>'); // <-- this is added to be selected as the first option, so it doesnt select a value.

                    // Retrieve the render function for the specific column
                    var renderFunction = that.settings()[0].aoColumns[col].mRender;

                    api.cells(null, col, {
                        search: 'applied'
                    }).data().unique().sort().each(function(d) {
                        // Use the render function and convert to string if necessary
                        var renderedValue = renderFunction ? renderFunction(d, 'display') : d; // Fallback to raw data if no render function

                        // Append the rendered value to the select2 options
                        $('#'+Table_id+'_Filter'+col).append($('<option>' + renderedValue + '</option>'));
                    });
                    // ###############################################################
                    // ## Remove the empty string value from the selection options  ##
                    // ##    This is the workaround to get fast loading without     ##
                    // ## selection of the first item. (process speed optimazation) ##
                    // ###############################################################
                    $('#'+Table_id+'_Filter'+col).attr('multiple', 'multiple');
                    $('#'+Table_id+'_Filter'+col+ " option")[0].remove();
                    
                    // ####################################
                    // ## Setup of select2 filter scheme ##
                    // ####################################
                    if (filter_type === 'Multi-select-unspecific') {
                        var select2 = $('#'+Table_id+'_Filter'+col).select2({
                            multiple: true,
                            closeOnSelect: true,
                            placeholder: {text: "Filter"},
                            dropdownAutoWidth : true,
                            tags: true, // allows for selection of undefined values --> Needs the not exact match to function well.
                        });
                    } else if (filter_type === 'Multi-select-exact') {
                        var select2 = $('#'+Table_id+'_Filter'+col).select2({
                            multiple: true,
                            closeOnSelect: true,
                            placeholder: {text: "Filter"},
                            dropdownAutoWidth : true,
                        });
                    }
                } // End of is searchable
            }); // End of multi-selct filter

            } else if (filter_type === 'Multi-select-exact-filter') {
                api.columns([column_number]).every(function () {
                    if (!this.searchable()) return;

                    var that = this;
                    var col = this.index();
                    var selected_cell = that.column(col).header();
                    var selectId = Table_id + '_Filter' + col;

                    // render header control
                    selected_cell.innerHTML = '<select id="'+selectId+'" class="select2" style="width: 80%;"></select>';

                    // helpers
                    function escapeRegex(s){ return String(s||'').replace(/[.*+?^${}()|[\]\\]/g, '\\$&'); }
                    function decodeHTML(s){
                    const div = document.createElement('div');
                    div.innerHTML = String(s || '');
                    return div.textContent || div.innerText || '';
                    }
                    function normalizeText(s){
                    // turn <br> to \n, strip tags, decode entities → plain text
                    const hasTags = /<[^>]+>/.test(s) || /&[a-z#0-9]+;/.test(s);
                    if (!hasTags) return String(s || '');
                    const div = document.createElement('div');
                    div.innerHTML = String(s || '').replace(/<br\s*\/?>/gi, '\n');
                    return div.textContent || div.innerText || '';
                    }

                    // ✅ collect tokens from the raw data when it's an array of names (best),
                    //    otherwise from the orthogonal "filter" text (fallback)
                    var tokenSet = new Set();

                    var cells = api.cells(null, col, { search: 'applied' });
                    cells.every(function () {
                    var raw = this.data();

                    if (Array.isArray(raw)) {
                        // expected for ligands / endo_ligands: [{id, name}, ...]
                        raw.forEach(function(x){
                        var nm = x && (x.name || x.text || x.label);
                        if (!nm) return;
                        nm = decodeHTML(nm).trim();
                        if (nm) tokenSet.add(nm);
                        });
                    } else {
                        // fallback: use rendered "filter" text, then split
                        var v = this.render('filter');
                        v = normalizeText(v);
                        v.split(/\s*(?:\n|,|\|)\s*/g).forEach(function(tok){
                        tok = tok && tok.trim();
                        if (tok) tokenSet.add(tok);
                        });
                    }
                    });

                    var tokens = Array.from(tokenSet).sort(function(a,b){
                    return a.localeCompare(b, undefined, { numeric: true, sensitivity: 'base' });
                    });

                    // build select (as text, not HTML)
                    var $sel = $('#'+selectId);

                    // Make it a real multi-select so the browser doesn't auto-select the first option
                    $sel.prop('multiple', true);

                    // Populate options (no empty option)
                    tokens.forEach(function(tok){
                    $('<option>', { value: tok, text: tok }).appendTo($sel);
                    });

                    // Init Select2
                    $sel.select2({
                    multiple: true,
                    closeOnSelect: true,
                    placeholder: { text: 'Filter' },
                    dropdownAutoWidth: true
                    });

                    // Ensure no selection on init (helps when stateSave or the browser cached a value)
                    $sel.val(null).trigger('change');
                    
                    // selection → exact token match within newline/comma/pipe-delimited list
                    $sel.on('change', function(){
                    var vals = ($(this).val() || []).map(escapeRegex);
                    var rx = '';
                    if (vals.length) {
                        // boundaries: start/end or whitespace/comma/pipe
                        rx = '(?:^|[\\s,|])(?:' + vals.join('|') + ')(?=$|[\\s,|])';
                    }
                    that.search(rx, true, false).draw();
                    });
                });

            // #########################################################
            // ## Range filter float (vertical and horizontal layout) ##
            // #########################################################
            
            } else if (filter_type === "Range-float-vertical" ||  filter_type === "Range-float-horizontal") {
                api.columns([column_number]).every(function() {
                if (this.searchable()) {
                    var that = this;
                    var col = column_number;
                    var selected_cell = that.column(col).header();
                    var html_input1 = '<input id="'+Table_id+'_Filter'+col+'min" class="select2" placeholder="min" style="width: 35px;text-align:center;margin-top: 5px;"></input>';
                    var html_input2 = '<input id="'+Table_id+'_Filter'+col+'max" class="select2" placeholder="max" style="width: 35px;text-align:center;"></input>';
                    // ##############################################
                    // ## Different setup (vertical or horizontal) ##
                    // ##############################################
                    if (filter_type === "Range-float-vertical") {
                    selected_cell.innerHTML = html_input2+'<br>'+html_input1;
                    } else if (filter_type === "Range-float-horizontal"){
                        selected_cell.innerHTML = html_input1+html_input2;
                    }
                    $('#'+Table_id+'_Filter'+col+'min'+','+'#'+Table_id+'_Filter'+col+'max').keyup(function() {
                        $.fn.dataTable.ext.search.push(
                            function( settings, data, dataIndex ) {
                                // ########################################################
                                // ## Don't filter on anything other than specific table ##
                                // ########################################################
                                if ( settings.nTable.id !== Table_id ) {

                                    return true;
                                
                                } else {
                                    
                                    // ##############################
                                    // ## Range filtering function ##
                                    // ##############################

                                    var min = parseFloat( $('#'+Table_id+'_Filter'+col+'min').val(), 10 );
                                    
                                    var max = parseFloat( $('#'+Table_id+'_Filter'+col+'max').val(), 10 );
                                    
                                    var data_range = parseFloat( data[col] ) || 0;
                            
                                    if ( ( isNaN( min ) && isNaN( max ) ) ||
                                        
                                        ( isNaN( min ) && data_range <= max ) ||
                                        
                                        ( min <= data_range   && isNaN( max ) ) ||
                                        
                                        ( min <= data_range   && data_range <= max ) )
                                    {
                                        return true;
                                    }
                                    return false;
                            }
                        });

                        // ########################################################
                        // ## After correct entries have been found redraw table ##
                        // ######################################################## 
                        
                        that.draw();       
                    }); // End of keyup function
                } // End of "is searchable" 
            }); // End of range filter
            
            // ############################################
            // ##          Range filter select           ##
            // ############################################
            
            } else if ( filter_type === "Range-select-vertical" || filter_type === "Range-select-horizontal") {
                api.columns([column_number]).every(function() {
                if (this.searchable()) {
                    var that = this;
                    var col = column_number;
                    var selected_cell = that.column(col).header();
                    
                    if (filter_type === "Range-select-vertical") {
                    var html_input1 = '<select id="'+Table_id+'_Filter'+col+'min" class="select2" data-placeholder="Min" style="width: 75%"></select>'
                    var html_input2 = '<select id="'+Table_id+'_Filter'+col+'max" class="select2" data-placeholder="Max" style="width: 75%;"></select>'
                    selected_cell.innerHTML = html_input1+'<br>'+html_input2;
                    } else if (filter_type === "Range-select-horizontal") {
                        var html_input1 = '<select id="'+Table_id+'_Filter'+col+'min" class="pull-left select2" style="width: 45%"></select>'
                        var html_input2 = '<select id="'+Table_id+'_Filter'+col+'max" class="pull-right select2" style="width: 45%;"></select>'
                        selected_cell.innerHTML = html_input1+html_input2;
                    }
                    $('#'+Table_id+'_Filter'+col+'min'+','+'#'+Table_id+'_Filter'+col+'max').on('select2:select', function() {
                            // Custom filtering function which will search data in column four between two values
                            $.fn.dataTable.ext.search.push(
                                function( settings, data, dataIndex ) {
                                    // ########################################################
                                    // ## Don't filter on anything other than specific table ##
                                    // ########################################################
                                    if ( settings.nTable.id !== Table_id ) {
    
                                        return true;
                                    
                                    } else {
                                        
                                        // ##############################
                                        // ## Range filtering function ##
                                        // ##############################
    
                                        var min = parseFloat( $('#'+Table_id+'_Filter'+col+'min').val(), 10 );
                                        
                                        var max = parseFloat( $('#'+Table_id+'_Filter'+col+'max').val(), 10 );
                                        
                                        var data_range = parseFloat( data[col] ) || 0;
                                
                                        if ( ( isNaN( min ) && isNaN( max ) ) ||
                                            
                                            ( isNaN( min ) && data_range <= max ) ||
                                            
                                            ( min <= data_range   && isNaN( max ) ) ||
                                            
                                            ( min <= data_range   && data_range <= max ) )
                                        {
                                            return true;
                                        }
                                        return false;
                                }
                            });

                            // ########################################################
                            // ## After correct entries have been found redraw table ##
                            // ######################################################## 
                            
                            that.draw();
                    });

                    $('#'+Table_id+'_Filter'+col+'min').append('<option></option>'); // <-- this is added to be selected as the first option, so it doesnt select a value.
                    $('#'+Table_id+'_Filter'+col+'max').append('<option></option>'); // <-- this is added to be selected as the first option, so it doesnt select a value.


                    api.cells(null, col, {search: 'applied'}).data().unique().sort().each(function(d) {
                        $('#'+Table_id+'_Filter'+col+'min').append($('<option>' + d + '</option>'));
                        $('#'+Table_id+'_Filter'+col+'max').append($('<option>' + d + '</option>'));
                    });

                    $('#'+Table_id+'_Filter'+col+'min').select2({
                            minimumResultsForSearch: -1,
                            multiple: false,
                            closeOnSelect: true,
                            dropdownAutoWidth : true,
                            placeholder: function(){
                                $(this).data('placeholder');
                            }
                    });
                    $('#'+Table_id+'_Filter'+col+'max').select2({
                        minimumResultsForSearch: -1,
                        multiple: false,
                        closeOnSelect: true,
                        dropdownAutoWidth : true,
                        placeholder: function(){
                            $(this).data('placeholder');
                        }
                    });
                }
            }); // End of Range filter year
        }

    }

    // ########################################
    // ## reset visability to original state ##
    // ########################################

    for (var j = 0; j < visibility_state.length; ++j) {
        api.column(j).visible(visibility_state[j]);
    }

    }
}


/**
 * @summary     HideEmptyColumns
 * @description Hide any (or specified) columns if no cells in the column(s)
 *              are populated with any values
 * @version     1.2.1
 * @file        dataTables.hideEmptyColumns.js
 * @author      Justin Hyland (http://www.justinhyland.com)
 * @contact     j@linux.com
 * @copyright   Copyright 2015 Justin Hyland
 * @url         https://github.com/jhyland87/DataTables-Hide-Empty-Columns
 *
 * License      MIT - http://datatables.net/license/mit
 *
 * Set the column visibility to hidden for any targeted columns that contain nothing
 * but null or empty values.
 *
 *
 * Parameters:
 *
 * -------------
 * hideEmptyCols
 *      Required:			true
 *      Type:				boolean|array|object
 *      Aliases:            hideEmptyColumns
 *      Description:		Primary setting, either target all columns, or specify an array for a list of cols, or an
 *                          object for advanced settings
 *      Examples:           hideEmptyCols: true
 *                          hideEmptyCols: [ 0, 2, 'age' ]
 *                          hideEmptyCols: { columns: [ 0, 2, 'age' ] }
 *                          hideEmptyCols: { columns: true }
 *
 * hideEmptyCols.columns
 *      Required:           false
 *      Type:               boolean|array
 *      Description:        Either true for all columns, or an array of indexes or dataSources
 *      Examples:           [ 0, 2, 'age' ]  // Column indexes 0 and 2, and the column name 'age'
 *                          true  // All columns
 *
 * hideEmptyCols.whiteList
 *      Required:           false
 *      Type:               boolean
 *      Default:            true
 *      Description:        Specify if the column list is to be targeted, or excluded
 *
 * hideEmptyCols.trim
 *      Required:           false
 *      Type:               boolean
 *      Default:            true
 *      Description:        Determines if column data values should be trimmed before checked for empty values
 *
 * hideEmptyCols.emptyVals
 *      Required:           false
 *      Type:               string|number|array|regex
 *      Description:        Define one or more values that will be interpreted as "empty"
 *      Examples:           [ '<br>', '</br>', '<BR>', '</BR>', '&nbsp;' ]  // HTML Line breaks, and the HTML NBSP char
 *                          /\<\/?br\>/i    // Any possible HTML line break character that matches this pattern
 *                          ['false', '0', /\<\/?br\>/i]   // The string values 'false' and '0', and all HTML breaks
 *
 * hideEmptyCols.onStateLoad
 *      Required:           false
 *      Type:               boolean
 *      Default:            true
 *      Description:        Determines if the main _checkColumns function should execute after the DT state is loaded
 *                          (when the DT stateSave option is enabled). This function will override the column visibility
 *                          state in stateSave
 *
 * hideEmptyCols.perPage
 *      Required:           false
 *      Type:               boolean
 *      Description:        Determine if columns should only be hidden if it has no values on the current page
 *
 *
 * @example
 *    // Target all columns - Hide any columns that contain all null/empty values
 *    $('#example').DataTable({
 *        hideEmptyCols: true
 *    })
 *
 * @example
 *    // Target the column indexes 0 & 2
 *    $('#example').DataTable({
 *        hideEmptyCols: [0,2]
 *    })
 *
 * @example
 *    // Target the column with 'age' data source
 *    $('#example').DataTable({
 *        ajax: 'something.js',
 *        hideEmptyCols: ['age'],
 *        buttons: [ 'columnsToggle' ],
 *        columns: [
 *            { name: 'name',     data: 'name' },
 *            { name: 'position', data: 'position' },
 *            { name: 'age',      data: 'age' }
 *        ]
 *    })
 *
 * @example
 *    // Target everything *except* the columns 1, 2 & 3
 *    $('#example').DataTable({
 *        hideEmptyCols: {
 *              columns: [ 1, 2, 3 ],
 *              whiteList: false
 *        }
 *    })
 *
 * @example
 *    // Target column indexes 1 and 4, adding custom empty values, and only hide the column if empty on current page
 *    $('#example').DataTable({
 *        hideEmptyCols: {
 *              columns: [ 1, 4 ],
 *              perPage: true,
 *              emptyVals: [ '0', /(no|false|disabled)/i ]
 *        }
 *    })
 */
"use strict"
;(function(window, document, $) {
    // On DT Initialization
    $(document).on('init.dt', function(e, dtSettings) {
        if ( e.namespace !== 'dt' )
            return

        // Check for either hideEmptyCols or hideEmptyColumns
        var options = dtSettings.oInit.hideEmptyCols || dtSettings.oInit.hideEmptyColumns

        // If neither of the above settings are found, then call it quits
        if( ! options )
            return

        // Helper function to get the value of a config item
        var _cfgItem = function( item, def ){
            if( $.isPlainObject( options ) && typeof options[ item ] !== 'undefined' )
                return options[ item ]

            return def
        }

        // Gather all the setting values which will be used
        var api         = new $.fn.dataTable.Api( dtSettings ),
            emptyCount  = 0,
            colList     = [],
            isWhiteList = ! _cfgItem( 'whiteList', false ),
            perPage     = _cfgItem( 'perPage' ),
            trimData    = _cfgItem( 'trim', true ),
            onStateLoad = _cfgItem( 'onStateLoad', true )

        // Helper function to determine if a cell is empty (including processing custom empty values)
        var _isEmpty = function( colData ){
            // Trim the data (unless its set to false)
            if( trimData )
                colData = $.trim( colData )

            // Basic check
            if( colData === null || colData.length === 0 )
                return true

            // Default to false, any empty matches will reset to true
            var retVal = false

            var emptyVals = _cfgItem( 'emptyVals' )

            // Internal helper function to check the value against a custom defined empty value (which can be a
            // regex pattern or a simple string)
            var _checkEmpty = function( val, emptyVal ){
                var objType = Object.prototype.toString.call( emptyVal )

                var match = objType.match( /^\[object\s(.*)\]$/ )

                // If its a regex pattern, then handle it differently
                if( match[1] === 'RegExp' )
                    return val.match( emptyVal )

                // Note: Should this comparison maybe use a lenient/loose comparison operator? hmm..
                return val === emptyVal
            }

            // If multiple custom empty values are defined in an array, then check each
            if( $.isArray( emptyVals ) ){
                $.each( emptyVals, function( i, ev ){
                    if( _checkEmpty( colData, ev ) )
                        retVal = true
                })
            }

            // Otherwise, just check the one, if set
            else if( typeof emptyVals !== 'undefined' ){
                if( _checkEmpty( colData, emptyVals ) )
                    retVal = true
            }

            return retVal
        }

        // If the hideEmptyCols setting is an Array (of column indexes to target)
        if( $.isArray( options ) ) {
            // And its populated..
            if( options.length !== 0 ) {
                $.each( options, function( k, i ){
                    // Try to get the real column index from whatever was configured
                    var indx = api.column( i ).index()

                    colList.push( typeof indx !== 'undefined' ? indx : i )
                })
            }
            else {
                // Otherwise, quit! since its just an empty array
                return
            }
        }

        // If hideEmptyCols setting is an Object (of plugin settings)
        else if( $.isPlainObject( options ) ) {
            // If options.columns isnt specifically
            if( typeof options.columns === 'undefined' || options.columns === true ) {
                // Set colList to true, enabling every column as a target
                colList = api.columns().indexes().toArray()
            }

            // If its an array, then it should contain the column indexs, so use that
            else if( $.isArray( options.columns ) ) {
                // Otherwise, set the colList
                colList = options.columns
            }

            // If `options.columns` isn't an array (of indexes) or a boolean (disable/enable all columns),
            // then throw a hissy fit
            else if( typeof options.columns !== 'boolean' ) {
                console.error( '[Hide Empty Columns]: Expected typeof `columns` setting value to be an array, boolean or undefined, but received value type "%s"', typeof options.columns )
                return
            }

            // The only thing left could be if its false, so just stop all together
            else {
                return
            }
        }

        // If its just a basic 'true' targeting all columns..
        else if( options === true ){
            // .. Then get the list of all column indexes
            colList = api.columns().indexes().toArray()
        }

        // Anything else should just go away
        else {
            return
        }

        // Function to check the column rows
        var _checkColumns = function( ) {
            var info = api.page.info(),
                colFilter = ( perPage ? { search: 'applied' } : undefined )

            // Iterate through the table, column by column
            //api.columns({ search: 'applied' }).every(function () {
            api.columns( colFilter ).every(function () {
                emptyCount = 0

                // If the current column is *not* found in the list..
                if( $.inArray( this.index(), colList ) === - 1 // Check column index #
                    && $.inArray( api.column( this.index() ).dataSrc(), colList ) === - 1 )// Check column name (dataSrc)
                {
                    // .. And the list type is whitelist, then skip this loop
                    if( isWhiteList === true ) return
                }
                // If the current column *is* found in the list..
                else {
                    // .. And the list type is blacklist, then skip this loop
                    if( isWhiteList === false ) return
                }

                // This gets ALL data in current column.. Need just the visible rows
                var data     = this.data().toArray(),
                    isVis    = false,
                    intStart = ( perPage === true && info.serverSide === false ? info.start : 0 ),
                    intStop  = ( perPage === true && info.serverSide === false ? info.end   : data.length ),
                    dtState  = api.state.loaded()

                //for( var i = 0; i < data.length; i ++ ) {
                for( var i = intStart; i < intStop; i++ ){
                    if( ! _isEmpty( data[ i ] ) ) {
                        isVis = true
                        break
                    }
                }

                // If the # of empty is the same as the length, then no values in col were found
                api.column( this.index() ).visible( isVis )

            } )
        }

        // If stateSave is enabled in this DT instance, then toggle the column visibility afterwords
        if( onStateLoad === true )
            api.on( 'stateLoadParams.dt', _checkColumns )

        // If were checking for each page, then attach functions to any events that may introduce or remove new
        // columns/rows from the table (page, order, search and length)
        if( perPage === true )
            api
                .on( 'page.dt',   _checkColumns )
                .on( 'search.dt', _checkColumns )
                .on( 'order.dt',  _checkColumns )
                .on( 'length.dt', _checkColumns )
                .on( 'draw.dt',   _checkColumns ) // triggers after data loaded with AJAX

        // Run check for the initial page load
        _checkColumns()
    })
})(window, document, jQuery)
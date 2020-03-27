$(document).ready(function () {
    'use strict';

    var oTable1,
        oTable2;

    oTable1 = $('#example1').DataTable({
        "paging": false,
        stateSave: true,
    });

    //----------------------------------------------
    //notice the new yadcf API for init the filters:
    //----------------------------------------------

    yadcf.init(oTable1, [{
        column_number: 0
    }, {
        column_number: 1,
        filter_type: "range_number_slider"
    }, {
        column_number: 2,
        filter_type: "date"
    }, {
        column_number: 3,
        filter_type: "auto_complete",
        text_data_delimiter: ","
    }, {
        column_number: 4,
        column_data_type: "html",
        html_data_type: "text",
        filter_default_label: "Select tag"
    }]);

    oTable2 = $('#example2').DataTable({
        "responsive": true,
        "processing": true,
        "ajax": "http://0.0.0.0:8000/static/sources/deep.txt",
        "columns": [{
            "data": "engine"
        }, {
            "data": "browser"
        }, {
            "data": "platform.inner"
        }, {
            "data": "platform.details.0"
        }, {
            "data": "platform.details.1"
        }]
    });

    //----------------------------------------------
    //notice the new yadcf API for init the filters:
    //----------------------------------------------

    yadcf.init(oTable2, [{
        column_number: 0
    }, {
        column_number: 1,
        filter_type: "text",
        exclude: true,
        exclude_label: '!(not)'
    }, {
        column_number: 2,
        filter_type: "auto_complete"
    }, {
        column_number: 3,
        filter_type: "range_number_slider",
        ignore_char: "-"
    }, {
        column_number: 4
    }]);

    // By default display the first tab. If this is not on one has to click on the tab for display.
    $('#myTab a:first').tab('show');

//    yadcf.exFilterColumn(oTable2, [[0, "Misc"]]);


    yadcf.initMultipleTables([oTable1, oTable2], [{
        filter_container_id: 'multi-table-filter',
        filter_default_label: 'Filter all tables!'
    }]);

});



$(function () {

    var oTable2 = $('#subtypesDT').DataTable({
        "scrollY": $(window).height() - 450,
        "scrollX": true,
        "scrollCollapse": true,
        "scroller": true,
        "paging": false,
        "bSortCellsTop": true, //prevent sort arrows going on bottom row
        "aaSorting": [],
        "autoWidth": true,
        "lengthMenu": [[10, 25, 50, -1], [10, 25, 50, "All"]],
        "pageLength": -1,
        "bInfo": true,
    });

    yadcf.init(oTable2,
        [
            {
                column_number: 1,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Class",
                filter_reset_button_text: false,
            },
            {
                column_number: 2,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Family",
                filter_reset_button_text: false,
            },
            {
                column_number: 3,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "uniprot",
                filter_reset_button_text: false,
            },
            {
                column_number: 4,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "IUPHAR",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '100px',
                },
            },
            {
                column_number: 5,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Bouvier",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '60px',
                },
            },
            {
                column_number: 6,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Inoue",
                filter_reset_button_text: false,
            },
            {
                column_number: 7,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gtp",
                filter_reset_button_text: false,
            },
            {
                column_number: 8,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Bouvier",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '60px',
                },
            },
            {
                column_number: 9,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Inoue",
                filter_reset_button_text: false,
            },
            {
                column_number: 10,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GtP",
                filter_reset_button_text: false,
            },
            {
                column_number: 11,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Bouvier",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '60px',
                },
            },
            {
                column_number: 12,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Inoue",
                filter_reset_button_text: false,
            },
            {
                column_number: 13,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "GtP",
                filter_reset_button_text: false,
            },
            {
                column_number: 14,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Bouvier",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '60px',
                },
            },
            {
                column_number: 15,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Inoue",
                filter_reset_button_text: false,
            },
            {
                column_number: 16,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Gtp",
                filter_reset_button_text: false,
            },
            {
                column_number: 17,
                filter_type: "multi_select",
                select_type: 'select2',
                filter_default_label: "Arrestin",
                filter_reset_button_text: false,
                select_type_options: {
                    width: '60px',
                },
            },
        ],
        {filters_tr_index: 1},
        {
            cumulative_filtering: true
        }
    );

    yadcf.exResetAllFilters(oTable2);



    yadcf.initMultipleTables([oTable1, oTable2], [{
        filter_container_id: 'multi-table-filter',
        filter_default_label: 'Filter all tables!'
    }]);


});


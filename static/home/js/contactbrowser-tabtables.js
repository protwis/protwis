var filtered_gn_pairs = [];

function filter_browser() {
    old_filtered_gn_pairs = filtered_gn_pairs;
    filtered_gn_pairs = [];

    if ($.fn.DataTable.isDataTable(".browser-table-1:visible")) {
        var table = $(".browser-table-1:visible").DataTable();
        table.rows({
            filter: 'applied'
        }).data().each(function(i) {
            filtered_gn_pairs.push(i['DT_RowId'])
        })
    }
    console.log('filtered positions! ', filtered_gn_pairs.length);

    if (old_filtered_gn_pairs.sort().join(',') === filtered_gn_pairs.sort().join(',')) {
        console.log('no change in filtering');
    } else {
        updateGeneralControls();
    }
}

function make_list_narrow_cols(list, start_column, end_column) {
    list_narrow_cols = Array(end_column - start_column + 1).fill().map((_, idx) => start_column + idx)
    list_narrow_cols = list.concat(list_narrow_cols);
    return list_narrow_cols;
}

function make_range_number_cols(start_column, repeat_number) {
    from_to = {
        filter_type: "range_number",
        filter_reset_button_text: false
    };
    repeated_from_to = []
    for (i = start_column; i < start_column + repeat_number; i++) {
        // var column_info = from_to.slice(0);
        var column_info = Object.assign({}, from_to);
        column_info['column_number'] = i;
        repeated_from_to.push(column_info);
    }
    return repeated_from_to;
}

function renderDataTablesYadcf(element) {

    console.time("renderDataTablesYadcf");

    const tab_number = element.split("-").reverse()[2];
    const selector = "#" + $('.main_option:visible').attr('id');
    const analys_mode = selector.replace('-tab', '');
    var table = $(selector + " .browser-table-" + tab_number);
    // If table is without tbody, then do not init further.
    if (!(table.find("thead").length)) {
        console.timeEnd("renderDataTablesYadcf");
        $(".main_loading_overlay").hide();
        return
    }
    // Do not re init the table.
    if ($.fn.DataTable.isDataTable(selector + " .browser-table-" + tab_number)) {
        console.timeEnd("renderDataTablesYadcf");
        $(".main_loading_overlay").hide();
        return
    }

    // $(".main_loading_overlay").show(0);
    switch (tab_number) {
        case "1":
            // statements_1
            // Specify which columns are to be fixed to 40px
            list_narrow_cols = [];
            if (analys_mode == "#two-crystal-groups") {
                list_narrow_cols = make_list_narrow_cols([2, 3, 4], 6, 6 + 19);
            }

            console.time("Render DataTable " + tab_number)
            btable = table.DataTable({
                'scrollX': true,
                scrollY: '50vh',
                // "sDom": 't', // To disable the pages on the button..
                paging: true,
                pageLength: 200,
                "bLengthChange": false,
                "bPaginate": false,
                "bInfo": false,
                "order": [],
                columnDefs: [{
                        type: "string",
                        targets: 1
                    },
                    {
                        "width": "40px",
                        "targets": list_narrow_cols
                    }
                ]
            });

            console.timeEnd("Render DataTable " + tab_number);

            if (analys_mode == "#two-crystal-groups") {

                repeated_from_to = make_range_number_cols(6, 20);

                yadcf.init(btable,
                    [{
                            column_number: 0,
                            filter_type: "text",
                            // exclude: true,
                            filter_delay: 500,
                            filter_reset_button_text: false,
                        },
                        {
                            column_number: 1,
                            filter_type: "text",
                            filter_default_label: "Pos",
                            filter_reset_button_text: false,
                        },
                        {
                            column_number: 2,
                            filter_type: "range_number",
                            filter_reset_button_text: false,
                        },
                        {
                            column_number: 3,
                            filter_type: "range_number",
                            filter_reset_button_text: false,

                        },
                        {
                            column_number: 4,
                            filter_type: "range_number",
                            filter_reset_button_text: false,

                        },
                        {
                            column_number: 5,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            filter_default_label: "Type",
                            text_data_delimiter: "|",
                            filter_reset_button_text: false,
                        }
                    ].concat(repeated_from_to), {
                        cumulative_filtering: false
                    }

                );
            } else if (analys_mode == "#single-crystal-group") {
                repeated_from_to = make_range_number_cols(4, 18);
                yadcf.init(btable,
                    [{
                            column_number: 0,
                            filter_type: "text",
                            // exclude: true,
                            filter_delay: 500,
                            filter_reset_button_text: false,
                        },
                        {
                            column_number: 1,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            filter_default_label: "Res.",
                            text_data_delimiter: "-",
                            filter_reset_button_text: false,
                        },
                        {
                            column_number: 2,
                            filter_type: "range_number",
                            filter_reset_button_text: false,
                        },
                        {
                            column_number: 3,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            filter_default_label: "Type",
                            text_data_delimiter: "|",
                            filter_reset_button_text: false,
                        }
                    ].concat(repeated_from_to), {
                        cumulative_filtering: false
                    }

                );
            } else if (analys_mode == "#single-crystal") {
                repeated_from_to = make_range_number_cols(4, 14);

                yadcf.init(btable,
                    [{
                            column_number: 0,
                            filter_type: "text",
                            // exclude: true,
                            filter_delay: 500,
                            filter_reset_button_text: false,
                        },
                        {
                            column_number: 1,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            // column_data_type: "html",
                            select_type_options: {
                                width: '60px'
                            },
                            html_data_type: "text",
                            text_data_delimiter: "-",
                            filter_default_label: "Res. No.",
                            filter_reset_button_text: false // hide yadcf reset button
                        },
                        {
                            column_number: 2,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '60px'
                            },
                            filter_default_label: "Res. Gn.",
                            text_data_delimiter: "-",
                            filter_reset_button_text: false,
                        },
                        {
                            column_number: 3,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '60px'
                            },
                            filter_default_label: "Type",
                            text_data_delimiter: "|",
                            filter_reset_button_text: false,
                        }
                    ].concat(repeated_from_to), {
                        cumulative_filtering: false
                    }

                );
            }
            btable.on('draw.dt', function(e, oSettings) {
                filter_browser();
            });
            btable.columns.adjust().draw();
            break;
        case "2":
            // statements_1
            console.log("do tab", tab_number)
            // Specify which columns are to be fixed to 40px
            list_narrow_cols = [];
            if (analys_mode == "#two-crystal-groups") {
                list_narrow_cols = make_list_narrow_cols([3, 4], 6, 6 + 17);
            }

            gray_scale_table(table);

            btable = table.DataTable({
                'scrollX': true,
                scrollY: '50vh',
                // "sDom": 't', // To disable the pages on the button..
                "bLengthChange": false,
                "bPaginate": false,
                "bInfo": false,
                paging: true,
                pageLength: 200,
                "order": [],
                columnDefs: [{
                        type: "string",
                        targets: 1
                    },
                    {
                        "width": "40px",
                        "targets": list_narrow_cols
                    }
                ]
            });


            if (analys_mode == "#two-crystal-groups") {


                repeated_from_to_1 = make_range_number_cols(2, 6);
                repeated_from_to_2 = make_range_number_cols(10, 12);
                repeated_from_to_3 = make_range_number_cols(23, 17);

                yadcf.init(btable,
                    [{
                            column_number: 0,
                            filter_type: "text",
                            // exclude: true,
                            filter_delay: 500,
                            filter_reset_button_text: false,
                        },
                        {
                            column_number: 1,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '80px'
                            },
                            filter_default_label: "Pos",
                            text_data_delimiter: "-",
                            filter_reset_button_text: false,
                        }
                    ].concat(repeated_from_to_1).concat([{
                        column_number: 8,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "AA",
                        filter_reset_button_text: false,
                    }, {
                        column_number: 9,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "AA",
                        filter_reset_button_text: false,
                    }]).concat(repeated_from_to_2).concat([{
                        column_number: 22,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "Type",
                        text_data_delimiter: "|",
                        filter_reset_button_text: false,
                    }]).concat(repeated_from_to_3), {
                        cumulative_filtering: false
                    }

                );
            } else if (analys_mode == "#single-crystal-group") {

                repeated_from_to_1 = make_range_number_cols(2, 2);
                repeated_from_to_2 = make_range_number_cols(6, 6);
                repeated_from_to_3 = make_range_number_cols(13, 17);

                yadcf.init(btable,
                    [{
                            column_number: 0,
                            filter_type: "text",
                            // exclude: true,
                            filter_delay: 500,
                            filter_reset_button_text: false,
                        },
                        {
                            column_number: 1,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '80px'
                            },
                            filter_default_label: "Pos",
                            text_data_delimiter: "-",
                            filter_reset_button_text: false,
                        }
                    ].concat(repeated_from_to_1).concat([{
                        column_number: 4,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "AA",
                        filter_reset_button_text: false,
                    }, {
                        column_number: 5,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "AA",
                        filter_reset_button_text: false,
                    }]).concat(repeated_from_to_2).concat([{
                        column_number: 12,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "Type",
                        text_data_delimiter: "|",
                        filter_reset_button_text: false,
                    }]).concat(repeated_from_to_3), {
                        cumulative_filtering: false
                    }

                );
            } else if (analys_mode == "#single-crystal") {

                repeated_from_to_1 = make_range_number_cols(4, 3);
                repeated_from_to_2 = make_range_number_cols(8, 13);

                yadcf.init(btable,
                    [{
                            column_number: 0,
                            filter_type: "text",
                            // exclude: true,
                            filter_delay: 500,
                            filter_reset_button_text: false,
                        },
                        {
                            column_number: 1,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '80px'
                            },
                            filter_default_label: "Pos",
                            text_data_delimiter: "-",
                            filter_reset_button_text: false,
                        },
                        {
                            column_number: 2,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '60px'
                            },
                            filter_default_label: "AA",
                            filter_reset_button_text: false,
                        }, {
                            column_number: 3,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '60px'
                            },
                            filter_default_label: "AA",
                            filter_reset_button_text: false,
                        }
                    ].concat(repeated_from_to_1).concat([{
                        column_number: 7,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "Type",
                        text_data_delimiter: "|",
                        filter_reset_button_text: false,
                    }]).concat(repeated_from_to_2), {
                        cumulative_filtering: false
                    }

                );
            }
            // btable.on('draw.dt', function(e, oSettings) {
            //     filter_browser();
            // });
            btable.columns.adjust().draw();
            break;

        case "3":
            // statements_1
            console.log("do tab", tab_number)
            // Specify which columns are to be fixed to 40px
            list_narrow_cols = [];
            if (analys_mode == "#two-crystal-groups") {
                list_narrow_cols = make_list_narrow_cols([2], 5, 7);
            }

            gray_scale_table(table);

            btable = table.DataTable({
                'scrollX': true,
                scrollY: '50vh',
                // "sDom": 't', // To disable the pages on the button..
                "bLengthChange": false,
                "bPaginate": false,
                "bInfo": false,
                paging: true,
                pageLength: 200,
                "order": [],
                columnDefs: [{
                        type: "string",
                        targets: 1
                    },
                    {
                        "width": "40px",
                        "targets": list_narrow_cols
                    }
                ]
            });


            if (analys_mode == "#two-crystal-groups") {


                repeated_from_to_1 = make_range_number_cols(2, 7);
                repeated_from_to_2 = make_range_number_cols(11, 3);
                repeated_from_to_3 = make_range_number_cols(15, 5);

                yadcf.init(btable,
                    [{
                            column_number: 0,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '80px'
                            },
                            filter_default_label: "Seg",
                            filter_reset_button_text: false,
                        },
                        {
                            column_number: 1,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '80px'
                            },
                            filter_default_label: "Pos",
                            filter_reset_button_text: false,
                        }
                    ].concat(repeated_from_to_1).concat([{
                        column_number: 9,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "AA",
                        filter_reset_button_text: false,
                    }, {
                        column_number: 10,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "AA",
                        filter_reset_button_text: false,
                    }]).concat(repeated_from_to_2).concat([{
                        column_number: 14,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "AA",
                        filter_reset_button_text: false,
                        filter_reset_button_text: false,
                    }]).concat(repeated_from_to_3), {
                        cumulative_filtering: false
                    }

                );
            } else if (analys_mode == "#single-crystal-group") {

            } else if (analys_mode == "#single-crystal") {

            }
            // btable.on('draw.dt', function(e, oSettings) {
            //     filter_browser();
            // });
            btable.columns.adjust().draw();
            break;
        case "4":
            // statements_1
            console.log("do tab", tab_number)
            // Specify which columns are to be fixed to 40px
            list_narrow_cols = [];
            // if (analys_mode == "#two-crystal-groups") {
            //     list_narrow_cols = make_list_narrow_cols([2], 5, 7);
            // }

            gray_scale_table(table);

            btable = table.DataTable({
                'scrollX': true,
                scrollY: '50vh',
                // "sDom": 't', // To disable the pages on the button..
                "bLengthChange": false,
                "bPaginate": false,
                "bInfo": false,
                paging: true,
                pageLength: 200,
                "order": [],
                columnDefs: [{
                        type: "string",
                        targets: 1
                    },
                    {
                        "width": "40px",
                        "targets": list_narrow_cols
                    }
                ]
            });

            if (analys_mode == "#two-crystal-groups") {
                repeated_from_to_1 = make_range_number_cols(2, 20);
                repeated_from_to_2 = make_range_number_cols(24, 3);
                repeated_from_to_3 = make_range_number_cols(28, 1);

                yadcf.init(btable,
                    [{
                            column_number: 0,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '80px'
                            },
                            filter_default_label: "Seg",
                            filter_reset_button_text: false,
                        },
                        {
                            column_number: 1,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '80px'
                            },
                            filter_default_label: "Pos",
                            filter_reset_button_text: false,
                        }
                    ].concat(repeated_from_to_1).concat([{
                        column_number: 22,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "AA",
                        filter_reset_button_text: false,
                    }, {
                        column_number: 23,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "AA",
                        filter_reset_button_text: false,
                    }]).concat(repeated_from_to_2).concat([{
                        column_number: 27,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "AA",
                        filter_reset_button_text: false,
                    }]).concat(repeated_from_to_3), {
                        cumulative_filtering: false
                    }

                );
            } else if (analys_mode == "#single-crystal-group") {
                repeated_from_to_1 = make_range_number_cols(3, 6);
                repeated_from_to_2 = make_range_number_cols(10, 1);
                repeated_from_to_3 = make_range_number_cols(12, 1);

                yadcf.init(btable,
                    [{
                            column_number: 0,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '80px'
                            },
                            filter_default_label: "Seg",
                            filter_reset_button_text: false,
                        },
                        {
                            column_number: 1,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '80px'
                            },
                            filter_default_label: "Pos",
                            filter_reset_button_text: false,
                        }
                    ].concat(repeated_from_to_1).concat([{
                        column_number: 9,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "AA",
                        filter_reset_button_text: false,
                    }]).concat(repeated_from_to_2).concat([{
                        column_number: 11,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "AA",
                        filter_reset_button_text: false,
                    }]).concat(repeated_from_to_3), {
                        cumulative_filtering: false
                    }

                );
            } else if (analys_mode == "#single-crystal") {
                repeated_from_to_1 = make_range_number_cols(3, 5);
                repeated_from_to_2 = make_range_number_cols(10, 1);

                yadcf.init(btable,
                    [{
                            column_number: 0,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '80px'
                            },
                            filter_default_label: "Seg",
                            filter_reset_button_text: false,
                        },
                        {
                            column_number: 1,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '80px'
                            },
                            filter_default_label: "Pos",
                            filter_reset_button_text: false,
                        }
                    ].concat(repeated_from_to_1).concat([{
                        column_number: 8,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "AA",
                        filter_reset_button_text: false,
                    }, {
                        column_number: 9,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "AA",
                        filter_reset_button_text: false,
                    }]).concat(repeated_from_to_2), {
                        cumulative_filtering: false
                    }

                );
            }
            // btable.on('draw.dt', function(e, oSettings) {
            //     filter_browser();
            // });
            btable.columns.adjust().draw();
            break;
        case "5":
            // statements_1
            console.log("do tab", tab_number)
            // Specify which columns are to be fixed to 40px
            list_narrow_cols = [];
            // if (analys_mode == "#two-crystal-groups") {
            //     list_narrow_cols = make_list_narrow_cols([2], 5, 7);
            // }

            gray_scale_table(table);

            btable = table.DataTable({
                'scrollX': true,
                scrollY: '50vh',
                // "sDom": 't', // To disable the pages on the button..
                "bLengthChange": false,
                "bPaginate": false,
                "bInfo": false,
                paging: true,
                pageLength: 200,
                "order": [],
                columnDefs: [{
                        type: "string",
                        targets: 1
                    },
                    {
                        "width": "40px",
                        "targets": list_narrow_cols
                    }
                ]
            });
            repeated_from_to = make_range_number_cols(2, 3);

            yadcf.init(btable,
                [{
                        column_number: 0,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '80px'
                        },
                        filter_default_label: "Seg",
                        filter_reset_button_text: false,
                    },
                    {
                        column_number: 1,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '80px'
                        },
                        filter_default_label: "Pos",
                        filter_reset_button_text: false,
                    }
                ].concat(repeated_from_to), {
                    cumulative_filtering: false
                }

            );
            if (analys_mode == "#two-crystal-groups") {

            } else if (analys_mode == "#single-crystal-group") {

            } else if (analys_mode == "#single-crystal") {

            }
            // btable.on('draw.dt', function(e, oSettings) {
            //     filter_browser();
            // });
            btable.columns.adjust().draw();
            break;
        default:
            // statements_def
            break;
    }
    // Show hidden tr now that table is rendered. (Faster rending, since less in DOM)
    table.find(".hidden").removeClass("hidden");
    $(".main_loading_overlay").hide();
    console.timeEnd("renderDataTablesYadcf");
}

const types_to_short = { 'ionic': 'Ion', 'aromatic': 'Aro', 'polar': 'Pol', 'hydrophobic': 'Hyd', 'van-der-waals': 'vDw' }

function renderBrowser(data) {
    console.time("RenderBrowser");
    var selector = $('ul#mode_nav').find('li.active').find('a').attr("href");
    var table = $(selector + " .browser-table-1");
    if ($.fn.DataTable.isDataTable(selector + " .browser-table-1")) {
        table.DataTable().destroy();
    }
    table.parent().html('<table class="browser-table-1 compact cell-border" style3="max-width: 2000px !important; width: 2000px !important" width2="2500px" style2="margin-left:0px" style="margin:0 auto"><thead></thead><tbody class="hidden"></tbody></table>');
    var table = $(selector + " .browser-table-1");
    // table.parent().before('<span><button type="button" onclick="filter_browser(this);" class="btn btn-xs btn-primary reset-selection">Filter</button></span>');
    var tbody = table.find('tbody');
    if (data['proteins2']) {

        // <th colspan="2">Ca distance from<br> 7TM axis</th> \
        // <th colspan="2">Backbone Rotation</th> \
        // <th colspan="2">Residue Rotamer</th> \
        // <th colspan="2">Tau angle</th> \
        // <th colspan="2">Phi dihedral</th> \
        // <th colspan="2">Psi dihedral</th> \
        thead = '<tr> \
                      <th colspan="2" class="skip"></th> \
                      <th colspan="3" class="pairselector" datatype="frequency"></th> \
                      <th colspan="1" class="skip"></th> \
                      <th colspan="1" class="selector" datatype="distance_diff"></th> \
                      <th colspan="2" class="selector" datatype="core_distance_diff"></th> \
                      <th colspan="2" class="selector" datatype="rotation_diff"></th> \
                      <th colspan="2" class="selector" datatype="rotamer_diff"></th> \
                      <th colspan="2" class="selector"datatype="SASA_diff"></th> \
                      <th colspan="2" class="selector"datatype="RSA_diff"></th> \
                      <th colspan="2" class="selector"datatype="presence_diff"></th> \
                      <th colspan="2" class="selector"datatype="consensus_SS"></th> \
                      <th colspan="2" class="selector"datatype="consensus_freq"></th> \
                      <th colspan="3" class="skip"></th> \
                  </tr> \
                  <tr> \
                          <th colspan="1" rowspan="2">Segment</th> \
                          <th colspan="1" rowspan="2">Positions</th> \
                          <th colspan="3" rowspan="2"> Frequency (%)</th> \
                          <th rowspan="2">Interactions</th> \
                          <th rowspan="2">Distance (Ca atoms)*</th> \
                          <th colspan="2">Backbone movement</th> \
                          <th colspan="2">(Ca-7TM axis)</th> \
                          <th colspan="2">Sidechain differences</th> \
                          <th colspan="2"></th> \
                          <th colspan="2"></th> \
                          <th colspan="2" rowspan="2">Position presence %</th> \
                          <th colspan="2">Secondary structure</th> \
                          <th colspan="2"></th> \
                          <th rowspan="2" colspan="3">Sum of conservation of contact AA pairs in class (%)</th> \
                        </tr> \
                        <tr> \
                          <th colspan="2">Distance</th> \
                          <th colspan="2">Rotation (Ca angle)</th> \
                          <th colspan="2">Rotamer</th> \
                          <th colspan="2">SASA</th> \
                          <th colspan="2">RSA</th> \
                          <th colspan="2">Consensus SS</th> \
                          <th colspan="2">Frequency %</th> \
                        </tr> \
                        <tr> \
                          <th class="dt-center"></th> \
                          <th class="dt-center">Pos1-Pos2</th> \
                          <th class="narrow_col">Set 1<br></th> \
                          <th class="narrow_col">Set 2<br></th> \
                          <th class="narrow_col">Diff<br></th> \
                          <th></th> \
                          <th class="narrow_col">Pos1-Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Set 1<br></th> \
                          <th class="narrow_col">Set 2<br></th> \
                          <th class="narrow_col">Diff<br></th> \
                        </tr>';
        table.find('thead').html(thead);
        // two groups
        var proteins_1 = data['proteins1'].length
        var proteins_2 = data['proteins2'].length
        var pdbs_1 = data['pdbs1'].length
        var pdbs_2 = data['pdbs2'].length
        $.each(data['interactions'], function(i, v) {
            var gn1 = i.split(",")[0]
            var gn2 = i.split(",")[1]
            var pfreq1 = Math.round(100 * v['proteins1'].length / proteins_1);
            var pfreq2 = Math.round(100 * v['proteins2'].length / proteins_2);
            var diff_pfreq = pfreq1 - pfreq2;
            var sfreq1 = Math.round(100 * v['pdbs1'].length / pdbs_1);
            var sfreq2 = Math.round(100 * v['pdbs2'].length / pdbs_2);
            var diff_sfreq = sfreq1 - sfreq2;
            var class_seq_cons = v['class_seq_cons'];

            var class_seq_cons_diff = class_seq_cons[0] - class_seq_cons[1];

            // var types = v['types'].join(",<br>");
            const types = v['types'].map((t) => types_to_short[t]).join('|');
            var seg1 = data['segm_lookup'][gn1];
            var seg2 = data['segm_lookup'][gn2];
            var distance = v['distance'];
            var angles_1 = v['angles'][0];
            var angles_2 = v['angles'][1];

            var pos1_presence = v['pos1_presence'];
            var pos2_presence = v['pos2_presence'];



            // console.log(gn1,gn2);
            all_angles_1 = data['all_angles'][gn1];
            all_angles_2 = data['all_angles'][gn2];
            ss_pos1_set1 = [];
            ss_pos1_set2 = [];
            ss_pos2_set1 = [];
            ss_pos2_set2 = [];
            pdbs = data['pdbs1'].concat(data['pdbs2']);
            pdbs.forEach(function(pdb){
                pdb_upper = pdb.toUpperCase();
                if (all_angles_1) {
                    let d1 = all_angles_1[pdb_upper];
                    if (d1.length) {
                        if (data['pdbs1'].includes(pdb)) {
                            ss_pos1_set1.push(d1[12]);
                        } else if (data['pdbs2'].includes(pdb)) {
                            ss_pos1_set2.push(d1[12]);
                        }
                    }
                }
                if (all_angles_2) {
                    let d2 = all_angles_2[pdb_upper];
                    if (d2.length) {
                        if (data['pdbs1'].includes(pdb)) {
                            ss_pos2_set1.push(d2[12])
                        } else if (data['pdbs2'].includes(pdb)) {
                            ss_pos2_set2.push(d2[12])
                        }
                    }
                }
            });

            dssp = [];
            [ss_pos1_set1,ss_pos1_set2,ss_pos2_set1,ss_pos2_set2].forEach(function(list){
                if (list.length) {
                    // make count map
                    c = new Map([...new Set(list)].map(
                            x => [x, list.filter(y => y === x).length]
                        ));
                    // make count with highest number first element
                    const mapSort = new Map([...c.entries()].sort((a, b) => b[1] - a[1]));
                    // get first element as [key,value]
                    most = mapSort.entries().next().value;
                    // calculate frequency
                    freq = most[1]/list.length;
                    most = most[0];
                } else {
                    freq = 0;
                    most = 'N/A';
                }
                dssp.push([most,freq]);
            })
            // console.table(dssp);
            dssp_pos1 = '';
            dssp_pos1_freq = '';
            if (dssp[0][0]==dssp[1][0]){
                // if pos1 category same
                dssp_pos1=dssp[0][0];
                dssp_pos1_freq = Math.round(100*(dssp[0][1]-dssp[1][1]));
            }


            dssp_pos2 = '';
            dssp_pos2_freq = '';
            if (dssp[2][0]==dssp[3][0]){
                // if pos2 category same
                dssp_pos2=dssp[2][0];
                dssp_pos2_freq = Math.round(100*(dssp[2][1]-dssp[3][1]));
            }
            // 0 'core_distance',
            // 1 'a_angle',
            // 2 'outer_angle',
            // 3 'tau',
            // 4 'phi',
            // 5 'psi',
            // 6 'sasa',
            // 7 'rsa',
            // 8 'theta',
            // 9 'hse'
            tr = `
                    <tr class="clickable-row filter_rows" id="${i}">
                      <td class="dt-center">${seg1}-${seg2}</td>
                      <td class="dt-center">${gn1}-${gn2}</td>
                      <td class="narrow_col">${sfreq1}</td>
                      <td class="narrow_col">${sfreq2}</td>
                      <td class="narrow_col">${diff_sfreq}</td>
                      <td>${types}</td>
                      <td class="narrow_col">${distance}</td>
                      <td class="narrow_col angles_modal">${angles_1[0]}</td>
                      <td class="narrow_col angles_modal">${angles_2[0]}</td>
                      <td class="narrow_col angles_modal">${angles_1[1]}</td>
                      <td class="narrow_col angles_modal">${angles_2[1]}</td>
                      <td class="narrow_col angles_modal">${angles_1[2]}</td>
                      <td class="narrow_col angles_modal">${angles_2[2]}</td>
                      <td class="narrow_col angles_modal">${angles_1[6]}</td>
                      <td class="narrow_col angles_modal">${angles_2[6]}</td>
                      <td class="narrow_col angles_modal">${angles_1[7]}</td>
                      <td class="narrow_col angles_modal">${angles_2[7]}</td>
                      <td class="narrow_col">${pos1_presence}</td>
                      <td class="narrow_col">${pos2_presence}</td>
                      <td class="narrow_col">${dssp_pos1}</td>
                      <td class="narrow_col">${dssp_pos2}</td>
                      <td class="narrow_col">${dssp_pos1_freq}</td>
                      <td class="narrow_col">${dssp_pos2_freq}</td>
                      <td class="narrow_col">${class_seq_cons[0]}</td>
                      <td class="narrow_col">${class_seq_cons[1]}</td>
                      <td class="narrow_col">${class_seq_cons_diff}</td>
                    </tr>`;
            tbody.append(tr);
        });
    } else if (data['proteins'].length > 1) {
        thead = '<tr> \
                      <th colspan="2" class="skip"></th> \
                      <th colspan="1" class="pairselector" datatype="frequency"></th> \
                      <th colspan="1" class="skip"></th> \
                      <th colspan="1" class="pairselector" datatype="distance"></th> \
                      <th colspan="2" class="selector" datatype="core_distance"></th> \
                      <th colspan="2" class="selector" datatype="rotation"></th> \
                      <th colspan="2" class="selector" datatype="rotamer"></th> \
                      <th colspan="2" class="selector" datatype="SASA"></th> \
                      <th colspan="2" class="selector" datatype="RSA"></th> \
                      <th colspan="2" class="selector" datatype="presence"></th> \
                      <th colspan="2" class="selector" datatype="consensus_SS"></th> \
                      <th colspan="2" class="selector" datatype="consensus_freq"></th> \
                      <th colspan="1" class="skip"></th> \
                  </tr> \
                  <tr> \
                          <th colspan="1" rowspan="2">Segment</th> \
                          <th colspan="1" rowspan="2">Positions</th> \
                          <th colspan="1" rowspan="2"> Frequency (%)</th> \
                          <th rowspan="2">Interactions</th> \
                          <th rowspan="2">Distance (Ca atoms)*</th> \
                          <th colspan="2">Backbone movement</th> \
                          <th colspan="2">(Ca-7TM axis)</th> \
                          <th colspan="2">Sidechain differences</th> \
                          <th colspan="2"></th> \
                          <th colspan="2"></th> \
                          <th colspan="2" rowspan="2">Position presence %</th> \
                          <th colspan="2">Secondary structure</th> \
                          <th colspan="2"></th> \
                          <th rowspan="2">Class Seq Cons(%)</th> \
                        </tr> \
                        <tr> \
                          <th colspan="2">Distance</th> \
                          <th colspan="2">Rotation (Ca angle)</th> \
                          <th colspan="2">Rotamer</th> \
                          <th colspan="2">SASA</th> \
                          <th colspan="2">RSA</th> \
                          <th colspan="2">Consensus SS</th> \
                          <th colspan="2">Frequency %</th> \
                        </tr> \
                        <tr> \
                          <th class="dt-center"></th> \
                          <th class="dt-center">Pos1-Pos2</th> \
                          <th class="narrow_col">Set<br></th> \
                          <th></th> \
                          <th class="narrow_col">Pos1-Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">AA pairs</th> \
                        </tr>';
        table.find('thead').html(thead);
        var proteins = data['proteins'].length
        var pdbs_counts = data['pdbs'].length
        $.each(data['interactions'], function(i, v) {
            var gn1 = i.split(",")[0]
            var gn2 = i.split(",")[1]
            var sfreq1 = Math.round(100 * v['pdbs'].length / pdbs_counts);
            var class_seq_cons = v['class_seq_cons'];
            // var types = v['types'].join(",<br>");
            const types = v['types'].map((t) => types_to_short[t]).join('|');
            var seg1 = data['segm_lookup'][gn1];
            var seg2 = data['segm_lookup'][gn2];
            var distance = v['distance'];
            var angles_1 = v['angles'][0];
            var angles_2 = v['angles'][1];
            var pos1_presence = v['pos1_presence'];
            var pos2_presence = v['pos2_presence'];

            all_angles_1 = data['all_angles'][gn1];
            all_angles_2 = data['all_angles'][gn2];
            ss_pos1_set1 = [];
            ss_pos2_set1 = [];
            pdbs = v['pdbs'];
            pdbs.forEach(function(pdb){
                pdb_upper = pdb.toUpperCase();
                if (all_angles_1) {
                    let d1 = all_angles_1[pdb_upper];
                    if (d1.length) {
                        ss_pos1_set1.push(d1[12]);
                    }
                }
                if (all_angles_2) {
                    let d2 = all_angles_2[pdb_upper];
                    if (d2.length) {
                        ss_pos2_set1.push(d2[12])
                    }
                }
            });

            dssp = [];
            [ss_pos1_set1,ss_pos2_set1].forEach(function(list){
                if (list.length) {
                    // make count map
                    c = new Map([...new Set(list)].map(
                            x => [x, list.filter(y => y === x).length]
                        ));
                    // make count with highest number first element
                    const mapSort = new Map([...c.entries()].sort((a, b) => b[1] - a[1]));
                    // get first element as [key,value]
                    most = mapSort.entries().next().value;
                    // calculate frequency
                    freq = most[1]/list.length;
                    most = most[0];
                } else {
                    freq = 0;
                    most = 'N/A';
                }
                dssp.push([most,freq]);
            })
            // console.table(dssp);
            dssp_pos1 = '';
            dssp_pos1_freq = '';
            if (dssp[0][0]){
                // if pos1 category same
                dssp_pos1=dssp[0][0];
                dssp_pos1_freq = Math.round(100*(dssp[0][1]));
            }


            dssp_pos2 = '';
            dssp_pos2_freq = '';
            if (dssp[1][0]){
                // if pos2 category same
                dssp_pos2=dssp[1][0];
                dssp_pos2_freq = Math.round(100*(dssp[1][1]));
            }

            // 0 'core_distance',
            // 1 'a_angle',
            // 2 'outer_angle',
            // 3 'tau',
            // 4 'phi',
            // 5 'psi',
            // 6 'sasa',
            // 7 'rsa',
            // 8 'theta',
            // 9 'hse'
            tr = `
                    <tr class="clickable-row filter_rows" id="${i}">
                      <td class="dt-center">${seg1}-${seg2}</td>
                      <td class="dt-center">${gn1}-${gn2}</td>
                      <td class="narrow_col">${sfreq1}</td>
                      <td>${types}</td>
                      <td class="narrow_col">${distance}</td>
                      <td class="narrow_col core_distance">${angles_1[0]}</td>
                      <td class="narrow_col core_distance">${angles_2[0]}</td>
                      <td class="narrow_col a_angle">${angles_1[1]}</td>
                      <td class="narrow_col a_angle">${angles_2[1]}</td>
                      <td class="narrow_col outer_angle">${angles_1[2]}</td>
                      <td class="narrow_col outer_angle">${angles_2[2]}</td>
                      <td class="narrow_col sasa">${angles_1[6]}</td>
                      <td class="narrow_col sasa">${angles_2[6]}</td>
                      <td class="narrow_col rsa">${angles_1[7]}</td>
                      <td class="narrow_col rsa">${angles_2[7]}</td>
                      <td class="narrow_col">${pos1_presence}</td>
                      <td class="narrow_col">${pos2_presence}</td>
                      <td class="narrow_col">${dssp_pos1}</td>
                      <td class="narrow_col">${dssp_pos2}</td>
                      <td class="narrow_col">${dssp_pos1_freq}</td>
                      <td class="narrow_col">${dssp_pos2_freq}</td>
                      <td class="narrow_col">${class_seq_cons}</td>
                    </tr>`;
            tbody.append(tr);
        });
    } else {
        thead = '<tr> \
                          <th colspan="1">Segment</th> \
                          <th colspan="1">Res No</th> \
                          <th colspan="1">Generic No</th> \
                          <th>Interaction</th> \
                          <th>Class Seq Cons(%)</th> \
                          <th>Ca distance</th> \
                          <th colspan="2">Ca distance from 7TM axis</th> \
                          <th colspan="2">Backbone<br>Rotation</th> \
                          <th colspan="2">Residue<br>Rotamer</th> \
                          <th colspan="2">Tau angle</th> \
                          <th colspan="2">Phi dihedral</th> \
                          <th colspan="2">Psi dihedral</th> \
                        </tr> \
                        <tr> \
                          <th class="dt-center"></th> \
                          <th class="dt-center"></th> \
                          <th class="dt-center"></th> \
                          <th></th> \
                          <th class="narrow_col">AA pairs</th> \
                          <th class="narrow_col">Res1-Res2</th> \
                          <th class="narrow_col">Res1</th> \
                          <th class="narrow_col">Res2</th> \
                          <th class="narrow_col">Res1</th> \
                          <th class="narrow_col">Res2</th> \
                          <th class="narrow_col">Res1</th> \
                          <th class="narrow_col">Res2</th> \
                          <th class="narrow_col">Res1</th> \
                          <th class="narrow_col">Res2</th> \
                          <th class="narrow_col">Res1</th> \
                          <th class="narrow_col">Res2</th> \
                          <th class="narrow_col">Res1</th> \
                          <th class="narrow_col">Res2</th> \
                        </tr>';

        thead = '<tr> \
                      <th colspan="4" class="skip"></th> \
                      <th colspan="1" class="pairselector" datatype="distance"></th> \
                      <th colspan="2" class="selector" datatype="core_distance"></th> \
                      <th colspan="2" class="selector" datatype="rotation"></th> \
                      <th colspan="2" class="selector" datatype="rotamer"></th> \
                      <th colspan="2" class="selector" datatype="SASA"></th> \
                      <th colspan="2" class="selector" datatype="RSA"></th> \
                      <th colspan="2" class="selector" datatype="consensus_SS"></th> \
                      <th colspan="1" class="skip"></th> \
                  </tr> \
                  <tr> \
                          <th colspan="1" rowspan="2">Segment</th> \
                          <th colspan="1" rowspan="2">Positions</th> \
                          <th colspan="1" rowspan="2">Positions GN</th> \
                          <th rowspan="2">Interaction</th> \
                          <th rowspan="2">Distance (Ca atoms)*</th> \
                          <th colspan="2">Backbone movement</th> \
                          <th colspan="2">(Ca-7TM axis)</th> \
                          <th colspan="2">Sidechain differences</th> \
                          <th colspan="2"></th> \
                          <th colspan="2"></th> \
                          <th colspan="2">Secondary structure</th> \
                          <th rowspan="2">Class Seq Cons(%)</th> \
                        </tr> \
                        <tr> \
                          <th colspan="2">Distance</th> \
                          <th colspan="2">Rotation (Ca angle)</th> \
                          <th colspan="2">Rotamer</th> \
                          <th colspan="2">SASA</th> \
                          <th colspan="2">RSA</th> \
                          <th colspan="2">Consensus SS</th> \
                        </tr> \
                        <tr> \
                          <th class="dt-center"></th> \
                          <th class="dt-center">Pos1-Pos2</th> \
                          <th class="narrow_col">Pos1-Pos2</th> \
                          <th></th> \
                          <th class="narrow_col">Pos1-Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">AA pairs</th> \
                        </tr>';
        table.find('thead').html(thead);
        $.each(data['interactions'], function(i, v) {
            var gn1 = i.split(",")[0]
            var gn2 = i.split(",")[1]
            // var types = v['types'].join(",<br>");
            const types = v['types'].map((t) => types_to_short[t]).join('|');
            var pos1 = v['seq_pos'][0];
            var pos2 = v['seq_pos'][1];
            var class_seq_cons = v['class_seq_cons'];
            var seg1 = data['segm_lookup'][gn1];
            var seg2 = data['segm_lookup'][gn2];
            var distance = v['distance'];
            var angles_1 = v['angles'][0];
            var angles_2 = v['angles'][1];

            all_angles_1 = data['all_angles'][gn1];
            all_angles_2 = data['all_angles'][gn2];

            if (all_angles_1) {
                all_angles_1 = all_angles_1[data['pdbs'][0].toUpperCase()][12];
            } else {
                all_angles_1 = '';
            }

            if (all_angles_2) {
                all_angles_2 = all_angles_2[data['pdbs'][0].toUpperCase()][12];
            } else {
                all_angles_2 = '';
            }
            //id="${pos1},${pos2}"
            tr = `
                    <tr class="clickable-row filter_rows" id="${i}">
                      <td class="dt-center">${seg1}-${seg2}</td>
                      <td class="dt-center"><span>${pos1}</span>-<span>${pos2}</span></td>
                      <td class="dt-center">${gn1}-${gn2}</td>
                      <td>${types}</td>
                      <td class="narrow_col">${distance}</td>
                      <td class="narrow_col angles_modal">${angles_1[0]}</td>
                      <td class="narrow_col angles_modal">${angles_2[0]}</td>
                      <td class="narrow_col angles_modal">${angles_1[1]}</td>
                      <td class="narrow_col angles_modal">${angles_2[1]}</td>
                      <td class="narrow_col angles_modal">${angles_1[2]}</td>
                      <td class="narrow_col angles_modal">${angles_2[2]}</td>
                      <td class="narrow_col angles_modal">${angles_1[6]}</td>
                      <td class="narrow_col angles_modal">${angles_2[6]}</td>
                      <td class="narrow_col angles_modal">${angles_1[7]}</td>
                      <td class="narrow_col angles_modal">${angles_2[7]}</td>
                      <td class="narrow_col">${all_angles_1}</td>
                      <td class="narrow_col">${all_angles_2}</td>
                      <td class="narrow_col">${class_seq_cons}</td>
                    </tr>`;
            tbody.append(tr);
        });
    }

    table.on('click', '.angles_modal', function(event) {
        // $(this)

        gn_pair = $(this).closest("tr").attr('id').split(",");
        gn1 = gn_pair[0];
        gn2 = gn_pair[1];

        all_angles_1 = two_sets_data['all_angles'][gn1];
        all_angles_2 = two_sets_data['all_angles'][gn2];

        $("#resModal").find(".modal-body").html("<div id='modal_plotly_1' style='height:100%;width:50%;display: inline-block;'></div><div id='modal_plotly_2' style='height:100%;width:50%;display: inline-block;'></div>");
        $("#resModal").modal();

        //Slight wait, to be sure modal is open.
        if (typeof all_angles_1 !== 'undefined')
          setTimeout(function(){ createBoxPlotResidue(all_angles_1,'modal_plotly_1','angles') }, 500);
        if (typeof all_angles_2 !== 'undefined')
          setTimeout(function(){ createBoxPlotResidue(all_angles_2,'modal_plotly_2','angles') }, 500);

    });

    console.timeEnd("RenderBrowser");
    gray_scale_table(table);
    enable_hover(table);
    //enable_3Dclick(table)
}

function renderBrowser_2(data) {
    console.time("RenderBrowser2");
    var selector = $('ul#mode_nav').find('li.active').find('a').attr("href");
    var table = $(selector + " .browser-table-2");
    if ($.fn.DataTable.isDataTable(selector + " .browser-table-2")) {
        table.DataTable().destroy();
    }
    table.parent().html('<table class="browser-table-2 compact cell-border" style3="max-width: 2000px !important; width: 2000px !important" width2="2500px" style2="margin-left:0px" style="margin:0 auto"><thead></thead><tbody class="hidden"></tbody></table>');
    var table = $(selector + " .browser-table-2");
    // table.parent().before('<span><button type="button" onclick="filter_browser(this);" class="btn btn-xs btn-primary reset-selection">Filter</button></span>');
    var tbody = table.find('tbody');
    if (data['proteins2']) {
        var proteins_1 = data['proteins1'].length
        var proteins_2 = data['proteins2'].length
        var pdbs_1 = data['pdbs1'].length
        var pdbs_2 = data['pdbs2'].length
        thead = '<tr> \
                          <th colspan="1" rowspan="2">Segment</th> \
                          <th colspan="1" rowspan="2">Positions</th> \
                          <th colspan="3" rowspan="2">Pos pair contact frequency (%)</th> \
                          <th colspan="3" rowspan="2">AA pair contact frequency (%)</th> \
                          <th colspan="2" rowspan="2">Amino acids</th> \
                          <th colspan="9" rowspan="1">AA occurrence in structure sets (%)</th> \
                          <th colspan="3" rowspan="2">Conservation in class (%)</th> \
                          <th rowspan="2">Interactions</th> \
                          <th rowspan="2">Distance (Ca atoms)*</th> \
                          <th colspan="4">Backbone movement (Ca-7TM axis)</th> \
                          <th colspan="6">Sidechain differences</th> \
                          <th colspan="2" rowspan="2">Position presence %</th> \
                          <th colspan="4">Secondary structure</th> \
                        </tr> \
                        <tr> \
                          <th colspan="3">AA1</th> \
                          <th colspan="3">AA2</th> \
                          <th colspan="3">AA Pair</th> \
                          <th colspan="2">Distance</th> \
                          <th colspan="2">Rotation (Ca angle)</th> \
                          <th colspan="2">Rotamer</th> \
                          <th colspan="2">SASA</th> \
                          <th colspan="2">RSA</th> \
                          <th colspan="2">Consensus SS</th> \
                          <th colspan="2">Frequency %</th> \
                        </tr> \
                        <tr> \
                          <th class="dt-center"></th> \
                          <th class="dt-center">Pos1-Pos2</th> \
                          <th class="narrow_col">Set 1<br></th> \
                          <th class="narrow_col">Set 2<br></th> \
                          <th class="narrow_col">Diff<br></th> \
                          <th class="narrow_col">Set 1<br></th> \
                          <th class="narrow_col">Set 2<br></th> \
                          <th class="narrow_col">Diff<br></th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Set 1<br></th> \
                          <th class="narrow_col">Set 2<br></th> \
                          <th class="narrow_col">Diff<br></th> \
                          <th class="narrow_col">Set 1<br></th> \
                          <th class="narrow_col">Set 2<br></th> \
                          <th class="narrow_col">Diff<br></th> \
                          <th class="narrow_col">Set 1<br></th> \
                          <th class="narrow_col">Set 2<br></th> \
                          <th class="narrow_col">Diff<br></th> \
                          <th class="narrow_col">AA1<br></th> \
                          <th class="narrow_col">AA2<br></th> \
                          <th class="narrow_col">Pair<br></th> \
                          <th></th> \
                          <th class="narrow_col">Pos1-Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                        </tr>';
        table.find('thead').html(thead);
        // two groups
        var proteins_1 = data['proteins1'].length
        var proteins_2 = data['proteins2'].length
        var pdbs_1 = data['pdbs1'].length
        var pdbs_2 = data['pdbs2'].length
        tr_list = ''
        $.each(data['tab2'], function(i, v2) {


            // GENERAL THINGS
            var v = data['interactions'][v2['pos_key']]
            var gn1 = v2['pos_key'].split(",")[0]
            var gn2 = v2['pos_key'].split(",")[1]
            var seg1 = data['segm_lookup'][gn1];
            var seg2 = data['segm_lookup'][gn2];
            // var pfreq1 = Math.round(100 * v['proteins1'].length / proteins_1);
            // var pfreq2 = Math.round(100 * v['proteins2'].length / proteins_2);
            // var diff_pfreq = pfreq1 - pfreq2;
            var sfreq1 = Math.round(100 * v['pdbs1'].length / pdbs_1);
            var sfreq2 = Math.round(100 * v['pdbs2'].length / pdbs_2);
            var diff_sfreq = sfreq1 - sfreq2;
            var class_seq_cons = v['class_seq_cons'];
            const types = v2['types'].map((t) => types_to_short[t]).join('|');
            var distance = v['distance'];
            var distance_2 = v2['distance'];
            var angles_1 = v2['angles'][0];
            var angles_2 = v2['angles'][1];

            // TAB-2 THINGS
            var aafreq1 = v2['set1']['interaction_freq'];
            var aafreq2 = v2['set2']['interaction_freq'];
            var diff_aafreq = (aafreq1 - aafreq2).toFixed(1);
            var aa1 = v2['aa1'];
            var aa2 = v2['aa2'];

            var set1_occurance_aa1 = Math.round(100 * v2['set1']['occurance']['aa1'].length / pdbs_1);
            var set1_occurance_aa2 = Math.round(100 * v2['set1']['occurance']['aa2'].length / pdbs_1);
            var set1_occurance_pair = Math.round(100 * v2['set1']['occurance']['pair'].length / pdbs_1);

            var set2_occurance_aa1 = Math.round(100 * v2['set2']['occurance']['aa1'].length / pdbs_2);
            var set2_occurance_aa2 = Math.round(100 * v2['set2']['occurance']['aa2'].length / pdbs_2);
            var set2_occurance_pair = Math.round(100 * v2['set2']['occurance']['pair'].length / pdbs_2);
            if (set1_occurance_aa1 > 100) {
                console.log(set1_occurance_aa1, v2['set1']['occurance']['aa1'].length, v2['set1']['occurance']['aa1'], pdbs_1, data['pdbs1']);
                console.log(set2_occurance_aa1, v2['set2']['occurance']['aa1'].length, v2['set2']['occurance']['aa1'], pdbs_2, data['pdbs2']);
            }
            var occurance_aa1_diff = set1_occurance_aa1 - set2_occurance_aa1;
            var occurance_aa2_diff = set1_occurance_aa2 - set2_occurance_aa2;
            var occurance_pair_diff = set1_occurance_pair - set2_occurance_pair;


            var pos1_presence = v['pos1_presence'];
            var pos2_presence = v['pos2_presence'];

            all_angles_1 = data['all_angles'][gn1];
            all_angles_2 = data['all_angles'][gn2];
            ss_pos1_set1 = [];
            ss_pos1_set2 = [];
            ss_pos2_set1 = [];
            ss_pos2_set2 = [];
            pdbs = data['pdbs1'].concat(data['pdbs2']);
            pdbs.forEach(function(pdb){
                pdb_upper = pdb.toUpperCase();
                if (all_angles_1) {
                    let d1 = all_angles_1[pdb_upper];
                    if (d1.length) {
                        if (v2['set1']['occurance']['aa1'].includes(pdb)) {
                            ss_pos1_set1.push(d1[12]);
                        } else if (v2['set2']['occurance']['aa1'].includes(pdb)) {
                            ss_pos1_set2.push(d1[12]);
                        }
                    }
                }
                if (all_angles_2) {
                    let d2 = all_angles_2[pdb_upper];
                    if (d2.length) {
                        if (v2['set1']['occurance']['aa2'].includes(pdb)) {
                            ss_pos2_set1.push(d2[12])
                        } else if (v2['set2']['occurance']['aa2'].includes(pdb)) {
                            ss_pos2_set2.push(d2[12])
                        }
                    }
                }
            });

            dssp = [];
            [ss_pos1_set1,ss_pos1_set2,ss_pos2_set1,ss_pos2_set2].forEach(function(list){
                if (list.length) {
                    // make count map
                    c = new Map([...new Set(list)].map(
                            x => [x, list.filter(y => y === x).length]
                        ));
                    // make count with highest number first element
                    const mapSort = new Map([...c.entries()].sort((a, b) => b[1] - a[1]));
                    // get first element as [key,value]
                    most = mapSort.entries().next().value;
                    // calculate frequency
                    freq = most[1]/list.length;
                    most = most[0];
                } else {
                    freq = 0;
                    most = 'N/A';
                }
                dssp.push([most,freq]);
            })
            // console.table(dssp);
            dssp_pos1 = '';
            dssp_pos1_freq = '';
            if (dssp[0][0]==dssp[1][0]){
                // if pos1 category same
                dssp_pos1=dssp[0][0];
                dssp_pos1_freq = Math.round(100*(dssp[0][1]-dssp[1][1]));
            }


            dssp_pos2 = '';
            dssp_pos2_freq = '';
            if (dssp[2][0]==dssp[3][0]){
                // if pos2 category same
                dssp_pos2=dssp[2][0];
                dssp_pos2_freq = Math.round(100*(dssp[2][1]-dssp[3][1]));
            }

            tr = ''
            tr_list += `
                    <tr class="clickable-row filter_rows" id="${i}">
                      <td class="dt-center">${seg1}-${seg2}</td>
                      <td class="dt-center">${gn1}-${gn2}</td>
                      <td class="narrow_col">${sfreq1}</td>
                      <td class="narrow_col">${sfreq2}</td>
                      <td class="narrow_col">${diff_sfreq}</td>
                      <td class="narrow_col">${aafreq1}</td>
                      <td class="narrow_col">${aafreq2}</td>
                      <td class="narrow_col">${diff_aafreq}</td>

                      <td class="narrow_col">${aa1}</td>
                      <td class="narrow_col">${aa2}</td>

                      <td class="narrow_col">${set1_occurance_aa1}</td>
                      <td class="narrow_col">${set2_occurance_aa1}</td>
                      <td class="narrow_col">${occurance_aa1_diff}</td>

                      <td class="narrow_col">${set1_occurance_aa2}</td>
                      <td class="narrow_col">${set2_occurance_aa2}</td>
                      <td class="narrow_col">${occurance_aa2_diff}</td>

                      <td class="narrow_col">${set1_occurance_pair}</td>
                      <td class="narrow_col">${set2_occurance_pair}</td>
                      <td class="narrow_col">${occurance_pair_diff}</td>

                      <td class="narrow_col">${v2['class_aa1']}</td>
                      <td class="narrow_col">${v2['class_aa2']}</td>
                      <td class="narrow_col">${v2['class']}</td>

                      <td>${types}</td>
                      <td class="narrow_col">${distance_2}</td>
                      <td class="narrow_col angles_modal">${angles_1[0]}</td>
                      <td class="narrow_col angles_modal">${angles_2[0]}</td>
                      <td class="narrow_col angles_modal">${angles_1[1]}</td>
                      <td class="narrow_col angles_modal">${angles_2[1]}</td>
                      <td class="narrow_col angles_modal">${angles_1[2]}</td>
                      <td class="narrow_col angles_modal">${angles_2[2]}</td>
                      <td class="narrow_col angles_modal">${angles_1[6]}</td>
                      <td class="narrow_col angles_modal">${angles_2[6]}</td>
                      <td class="narrow_col angles_modal">${angles_1[7]}</td>
                      <td class="narrow_col angles_modal">${angles_2[7]}</td>
                      <td class="narrow_col">${pos1_presence}</td>
                      <td class="narrow_col">${pos2_presence}</td>
                      <td class="narrow_col">${dssp_pos1}</td>
                      <td class="narrow_col">${dssp_pos2}</td>
                      <td class="narrow_col">${dssp_pos1_freq}</td>
                      <td class="narrow_col">${dssp_pos2_freq}</td>
                    </tr>`;
            // tbody.append(tr);
        });
        // insert natively for speed increase on Chrome
        tbody[0].innerHTML = tr_list;
    } else if (data['proteins'].length > 1) {
        thead = '<tr> \
                          <th colspan="1" rowspan="2">Segment</th> \
                          <th colspan="1" rowspan="2">Positions</th> \
                          <th colspan="1" rowspan="2">Pos pair contact frequency (%)</th> \
                          <th colspan="1" rowspan="2">AA pair contact frequency (%)</th> \
                          <th colspan="2" rowspan="2">Amino acids</th> \
                          <th colspan="3" rowspan="1">AA occurrence in set (%)</th> \
                          <th colspan="3" rowspan="2">Conservation in class (%)</th> \
                          <th rowspan="2">Interactions</th> \
                          <th rowspan="2">Distance (Ca atoms)*</th> \
                          <th colspan="4">Backbone movement (Ca-7TM axis)</th> \
                          <th colspan="6">Sidechain differences</th> \
                          <th colspan="2" rowspan="2">Position presence %</th> \
                          <th colspan="4">Secondary structure</th> \
                        </tr> \
                        <tr> \
                          <th colspan="1">AA1</th> \
                          <th colspan="1">AA2</th> \
                          <th colspan="1">AA Pair</th> \
                          <th colspan="2">Distance</th> \
                          <th colspan="2">Rotation (Ca angle)</th> \
                          <th colspan="2">Rotamer</th> \
                          <th colspan="2">SASA</th> \
                          <th colspan="2">RSA</th> \
                          <th colspan="2">Consensus SS</th> \
                          <th colspan="2">Frequency %</th> \
                        </tr> \
                        <tr> \
                          <th class="dt-center"></th> \
                          <th class="dt-center">Pos1-Pos2</th> \
                          <th class="narrow_col"><br></th> \
                          <th class="narrow_col"><br></th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col"><br></th> \
                          <th class="narrow_col"><br></th> \
                          <th class="narrow_col"><br></th> \
                          <th class="narrow_col">AA1<br></th> \
                          <th class="narrow_col">AA2<br></th> \
                          <th class="narrow_col">Pair<br></th> \
                          <th></th> \
                          <th class="narrow_col">Pos1-Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                        </tr>';
        table.find('thead').html(thead);
        // two groups
        var proteins = data['proteins'].length
        var pdbs_count = data['pdbs'].length
        tr_list = ''
        $.each(data['tab2'], function(i, v2) {


            // GENERAL THINGS
            var v = data['interactions'][v2['pos_key']]
            var gn1 = v2['pos_key'].split(",")[0]
            var gn2 = v2['pos_key'].split(",")[1]
            var seg1 = data['segm_lookup'][gn1];
            var seg2 = data['segm_lookup'][gn2];
            // var pfreq1 = Math.round(100 * v['proteins1'].length / proteins_1);
            // var pfreq2 = Math.round(100 * v['proteins2'].length / proteins_2);
            // var diff_pfreq = pfreq1 - pfreq2;
            var sfreq = Math.round(100 * v['pdbs'].length / pdbs_count);
            var class_seq_cons = v['class_seq_cons'];
            const types = v2['types'].map((t) => types_to_short[t]).join('|');
            var distance = v['distance'];
            var distance_2 = v2['distance'];
            var angles_1 = v2['angles'][0];
            var angles_2 = v2['angles'][1];

            // TAB-2 THINGS
            var aafreq = v2['set']['interaction_freq'];
            var aa1 = v2['aa1'];
            var aa2 = v2['aa2'];

            var set1_occurance_aa1 = Math.round(100 * v2['set']['occurance']['aa1'].length / pdbs_count);
            var set1_occurance_aa2 = Math.round(100 * v2['set']['occurance']['aa2'].length / pdbs_count);
            var set1_occurance_pair = Math.round(100 * v2['set']['occurance']['pair'].length / pdbs_count);

            var pos1_presence = v['pos1_presence'];
            var pos2_presence = v['pos2_presence'];

            all_angles_1 = data['all_angles'][gn1];
            all_angles_2 = data['all_angles'][gn2];
            ss_pos1_set = [];
            ss_pos2_set = [];
            pdbs = data['pdbs'];
            pdbs.forEach(function(pdb){
                pdb_upper = pdb.toUpperCase();
                if (all_angles_1) {
                    let d1 = all_angles_1[pdb_upper];
                    if (d1.length) {
                        if (v2['set']['occurance']['aa1'].includes(pdb)) {
                            ss_pos1_set.push(d1[12]);
                        }
                    }
                }
                if (all_angles_2) {
                    let d2 = all_angles_2[pdb_upper];
                    if (d2.length) {
                        if (v2['set']['occurance']['aa2'].includes(pdb)) {
                            ss_pos2_set.push(d2[12])
                        }
                    }
                }
            });

            dssp = [];
            [ss_pos1_set,ss_pos2_set].forEach(function(list){
                if (list.length) {
                    // make count map
                    c = new Map([...new Set(list)].map(
                            x => [x, list.filter(y => y === x).length]
                        ));
                    // make count with highest number first element
                    const mapSort = new Map([...c.entries()].sort((a, b) => b[1] - a[1]));
                    // get first element as [key,value]
                    most = mapSort.entries().next().value;
                    // calculate frequency
                    freq = most[1]/list.length;
                    most = most[0];
                } else {
                    freq = 0;
                    most = 'N/A';
                }
                dssp.push([most,freq]);
            })
            // console.table(dssp);
            dssp_pos1 = '';
            dssp_pos1_freq = '';
            if (dssp[0][0]){
                // if pos1 category same
                dssp_pos1=dssp[0][0];
                dssp_pos1_freq = Math.round(100*dssp[0][1]);
            }


            dssp_pos2 = '';
            dssp_pos2_freq = '';
            if (dssp[1][0]){
                // if pos2 category same
                dssp_pos2=dssp[1][0];
                dssp_pos2_freq = Math.round(100*dssp[1][1]);
            }

            tr = ''
            tr_list += `
                    <tr class="clickable-row filter_rows" id="${i}">
                      <td class="dt-center">${seg1}-${seg2}</td>
                      <td class="dt-center">${gn1}-${gn2}</td>
                      <td class="narrow_col">${sfreq}</td>
                      <td class="narrow_col">${aafreq}</td>

                      <td class="narrow_col">${aa1}</td>
                      <td class="narrow_col">${aa2}</td>

                      <td class="narrow_col">${set1_occurance_aa1}</td>
                      <td class="narrow_col">${set1_occurance_aa2}</td>
                      <td class="narrow_col">${set1_occurance_pair}</td>

                      <td class="narrow_col">${v2['class_aa1']}</td>
                      <td class="narrow_col">${v2['class_aa2']}</td>
                      <td class="narrow_col">${v2['class']}</td>

                      <td>${types}</td>
                      <td class="narrow_col">${distance_2}</td>
                      <td class="narrow_col core_distance">${angles_1[0]}</td>
                      <td class="narrow_col core_distance">${angles_2[0]}</td>
                      <td class="narrow_col a_angle">${angles_1[1]}</td>
                      <td class="narrow_col a_angle">${angles_2[1]}</td>
                      <td class="narrow_col outer_angle">${angles_1[2]}</td>
                      <td class="narrow_col outer_angle">${angles_2[2]}</td>
                      <td class="narrow_col sasa">${angles_1[6]}</td>
                      <td class="narrow_col sasa">${angles_2[6]}</td>
                      <td class="narrow_col rsa">${angles_1[7]}</td>
                      <td class="narrow_col rsa">${angles_2[7]}</td>
                      <td class="narrow_col">${pos1_presence}</td>
                      <td class="narrow_col">${pos2_presence}</td>
                      <td class="narrow_col">${dssp_pos1}</td>
                      <td class="narrow_col">${dssp_pos2}</td>
                      <td class="narrow_col">${dssp_pos1_freq}</td>
                      <td class="narrow_col">${dssp_pos2_freq}</td>
                    </tr>`;
            tbody.append(tr);
        });
        // insert natively for speed increase
        tbody[0].innerHTML = tr_list;
    } else {
        thead = '<tr> \
                          <th colspan="1" rowspan="2">Segment</th> \
                          <th colspan="1" rowspan="2">Positions</th> \
                          <th colspan="2" rowspan="2">Amino acids</th> \
                          <th colspan="3" rowspan="2">Conservation in class (%)</th> \
                          <th rowspan="2">Interactions</th> \
                          <th rowspan="2">Distance (Ca atoms)*</th> \
                          <th colspan="4">Backbone movement (Ca-7TM axis)</th> \
                          <th colspan="6">Sidechain differences</th> \
                          <th colspan="2">Secondary structure</th> \
                        </tr> \
                        <tr> \
                          <th colspan="2">Distance</th> \
                          <th colspan="2">Rotation (Ca angle)</th> \
                          <th colspan="2">Rotamer</th> \
                          <th colspan="2">SASA</th> \
                          <th colspan="2">RSA</th> \
                          <th colspan="2">Consensus SS</th> \
                        </tr> \
                        <tr> \
                          <th class="dt-center"></th> \
                          <th class="dt-center">Pos1-Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">AA1<br></th> \
                          <th class="narrow_col">AA2<br></th> \
                          <th class="narrow_col">Pair<br></th> \
                          <th></th> \
                          <th class="narrow_col">Pos1-Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                        </tr>';
        table.find('thead').html(thead);
        // two groups
        //var proteins = data['proteins'].length
        //var pdbs = data['pdbs'].length
        tr_list = ''
        $.each(data['tab2'], function(i, v2) {


            // GENERAL THINGS
            var v = data['interactions'][v2['pos_key']]
            var gn1 = v2['pos_key'].split(",")[0]
            var gn2 = v2['pos_key'].split(",")[1]
            var seg1 = data['segm_lookup'][gn1];
            var seg2 = data['segm_lookup'][gn2];
            // var pfreq1 = Math.round(100 * v['proteins1'].length / proteins_1);
            // var pfreq2 = Math.round(100 * v['proteins2'].length / proteins_2);
            // var diff_pfreq = pfreq1 - pfreq2;
            var class_seq_cons = v['class_seq_cons'];
            const types = v2['types'].map((t) => types_to_short[t]).join('|');
            var distance = v['distance'];
            var distance_2 = v2['distance'];
            var angles_1 = v2['angles'][0];
            var angles_2 = v2['angles'][1];

            // TAB-2 THINGS
            var aa1 = v2['aa1'];
            var aa2 = v2['aa2'];

            all_angles_1 = data['all_angles'][gn1];
            all_angles_2 = data['all_angles'][gn2];

            if (all_angles_1) {
                all_angles_1 = all_angles_1[data['pdbs'][0].toUpperCase()][12];
            } else {
                all_angles_1 = '';
            }

            if (all_angles_2) {
                all_angles_2 = all_angles_2[data['pdbs'][0].toUpperCase()][12];
            } else {
                all_angles_2 = '';
            }

            tr = ''
            tr_list += `
                    <tr class="clickable-row filter_rows" id="${i}">
                      <td class="dt-center">${seg1}-${seg2}</td>
                      <td class="dt-center">${gn1}-${gn2}</td>

                      <td class="narrow_col">${aa1}</td>
                      <td class="narrow_col">${aa2}</td>

                      <td class="narrow_col">${v2['class_aa1']}</td>
                      <td class="narrow_col">${v2['class_aa2']}</td>
                      <td class="narrow_col">${v2['class']}</td>

                      <td>${types}</td>
                      <td class="narrow_col">${distance_2}</td>
                      <td class="narrow_col core_distance">${angles_1[0]}</td>
                      <td class="narrow_col core_distance">${angles_2[0]}</td>
                      <td class="narrow_col a_angle">${angles_1[1]}</td>
                      <td class="narrow_col a_angle">${angles_2[1]}</td>
                      <td class="narrow_col outer_angle">${angles_1[2]}</td>
                      <td class="narrow_col outer_angle">${angles_2[2]}</td>
                      <td class="narrow_col sasa">${angles_1[6]}</td>
                      <td class="narrow_col sasa">${angles_2[6]}</td>
                      <td class="narrow_col rsa">${angles_1[7]}</td>
                      <td class="narrow_col rsa">${angles_2[7]}</td>
                      <td class="narrow_col">${all_angles_1}</td>
                      <td class="narrow_col">${all_angles_2}</td>
                    </tr>`;
            tbody.append(tr);
        });
        // insert natively for speed increase
        tbody[0].innerHTML = tr_list;
    }
    table.on('click', '.angles_modal', function(event) {
        // $(this)
        i = $(this).closest("tr").attr('id');
        v2 = data['tab2'][i];
        var gn1 = v2['pos_key'].split(",")[0]
        var gn2 = v2['pos_key'].split(",")[1]
        var pdbs_aa1 = v2['set1']['occurance']['aa1'].concat(v2['set2']['occurance']['aa1']);
        var pdbs_aa2 = v2['set1']['occurance']['aa2'].concat(v2['set2']['occurance']['aa2']);
        var aa1 = v2['aa1'];
        var aa2 = v2['aa2'];

        all_angles_1 = two_sets_data['all_angles'][gn1];
        all_angles_2 = two_sets_data['all_angles'][gn2];

        $("#resModal").find(".modal-body").html("<div id='modal_plotly_1' style='height:100%;width:50%;display: inline-block;'></div><div id='modal_plotly_2' style='height:100%;width:50%;display: inline-block;'></div>");
        $("#resModal").modal();

        //Slight wait, to be sure modal is open.
        console.log(pdbs_aa1,pdbs_aa2);
        setTimeout(function(){ createBoxPlotResidue(all_angles_1,'modal_plotly_1','angles',pdbs_aa1,aa1) }, 500);
        setTimeout(function(){ createBoxPlotResidue(all_angles_2,'modal_plotly_2','angles',pdbs_aa2,aa2) }, 500);

    });

    console.timeEnd("RenderBrowser2");
}

function renderBrowser_3(data) {
    console.time("RenderBrowser3");
    var selector = $('ul#mode_nav').find('li.active').find('a').attr("href");
    var table = $(selector + " .browser-table-3");
    if ($.fn.DataTable.isDataTable(selector + " .browser-table-3")) {
        table.DataTable().destroy();
    }
    table.parent().html('<table class="browser-table-3 compact cell-border" style3="max-width: 2000px !important; width: 2000px !important" width2="2500px" style2="margin-left:0px" style="margin:0 auto"><thead></thead><tbody class="hidden"></tbody></table>');
    var table = $(selector + " .browser-table-3");
    // table.parent().before('<span><button type="button" onclick="filter_browser(this);" class="btn btn-xs btn-primary reset-selection">Filter</button></span>');
    var tbody = table.find('tbody');
    var proteins_1 = data['proteins1'].length
    var proteins_2 = data['proteins2'].length
    var pdbs_1 = data['pdbs1'].length
    var pdbs_2 = data['pdbs2'].length
    if (data['proteins2']) {
        thead = '<tr> \
                      <th colspan="2" class="skip"></th> \
                      <th colspan="3" class="selector" datatype="contacts"></th> \
                      <th colspan="3" class="selector" datatype="contacts"></th> \
                      <th colspan="1" class="selector" datatype="contacts"></th> \
                      <th colspan="2" class="skip"></th> \
                      <th colspan="3" class="selector" datatype="conservation"></th> \
                      <th colspan="1" class="skip"></th> \
                      <th colspan="1" class="selector" datatype="conservation"></th> \
                      <th colspan="1" class="selector" datatype="core_distance_diff"></th> \
                      <th colspan="1" class="selector" datatype="rotation_diff"></th> \
                      <th colspan="1" class="selector" datatype="rotamer_diff"></th> \
                      <th colspan="1" class="selector"datatype="SASA_diff"></th> \
                  </tr> \
                  <tr> \
                          <th colspan="1" rowspan="2">Segment</th> \
                          <th colspan="1" rowspan="2">Positions</th> \
                          <th colspan="3" rowspan="2">Avg no contact pairs</th> \
                          <th colspan="3" rowspan="1">No distinct or shared contacts</th> \
                          <th colspan="1" rowspan="2">Avg freq difference of all set1-2 contacts</th> \
                          <th colspan="5" rowspan="1">Seq consensus</th> \
                          <th colspan="2" rowspan="1">Class seq consensus</th> \
                          <th colspan="3">Backbone movement (Ca-7TM axis)</th> \
                          <th colspan="2">Sidechain differences</th> \
                        </tr> \
                        <tr> \
                          <th colspan="2">Distinct</th> \
                          <th colspan="1">Shared</th> \
                          <th colspan="2">AA</th> \
                          <th colspan="3">Conservation (%)</th> \
                          <th colspan="1">AA</th> \
                          <th colspan="1">Cons (%)</th> \
                          <th colspan="1">Distance*</th> \
                          <th colspan="1">Rotation</th> \
                          <th colspan="1">Rotamer</th> \
                          <th colspan="1">SASA</th> \
                        </tr> \
                        <tr> \
                          <th class="dt-center"></th> \
                          <th class="dt-center">Pos</th> \
                          <th class="narrow_col">Set 1<br></th> \
                          <th class="narrow_col">Set 2<br></th> \
                          <th class="narrow_col">Diff<br></th> \
                          <th class="narrow_col">Set 1<br></th> \
                          <th class="narrow_col">Set 2<br></th> \
                          <th class="narrow_col">Set 1&2<br></th> \
                          <th class="narrow_col">Set1-2</th> \
                          <th class="narrow_col">Set 1<br></th> \
                          <th class="narrow_col">Set 2<br></th> \
                          <th class="narrow_col">Set 1<br></th> \
                          <th class="narrow_col">Set 2<br></th> \
                          <th class="narrow_col">Diff<br></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                        </tr>';
        table.find('thead').html(thead);
        // two groups
        var proteins_1 = data['proteins1'].length
        var proteins_2 = data['proteins2'].length
        var pdbs_1 = data['pdbs1'].length
        var pdbs_2 = data['pdbs2'].length
        tr_list = ''
        $.each(data['tab3'], function(i, v) {

            // console.log(i,v);
            var seg = data['segm_lookup'][i];

            var avg_set1 = Math.round(10 * v['set1_count'].length / pdbs_1) / 10;
            var avg_set2 = Math.round(10 * v['set2_count'].length / pdbs_2) / 10;
            var avg_diff = Math.round(10 * (avg_set1 - avg_set2)) / 10;

            var unique_set1 = v['set1'].filter(function(obj) { return v['set2'].indexOf(obj) == -1; });
            var unique_set2 = v['set2'].filter(function(obj) { return v['set1'].indexOf(obj) == -1; });
            var common = v['set1'].filter(value => -1 !== v['set2'].indexOf(value))

            var avg_freq_diff_sets = Math.round(10 * v['avg_freq_diff_sets']) / 10;

            var set1_seq_cons_aa = v['set1_seq_cons'][0];
            var set2_seq_cons_aa = v['set2_seq_cons'][0];
            var set1_seq_cons_freq = Math.round(100 * v['set1_seq_cons'][1] / pdbs_1);
            var set2_seq_cons_freq = Math.round(100 * v['set2_seq_cons'][2] / pdbs_2);
            var diff_seq_cons_freq = Math.round((set1_seq_cons_freq - set2_seq_cons_freq));

            var class_cons_aa = v['class_cons'][0];
            var class_cons_freq = Math.round(100 * v['class_cons'][1]);


            var angles = v['angles'];
            // 0 'core_distance',
            // 1 'a_angle',
            // 2 'outer_angle',
            // 3 'tau',
            // 4 'phi',
            // 5 'psi',
            // 6 'sasa',
            // 7 'rsa',
            // 8 'theta',
            // 9 'hse'

            tr = ''
            tr_list += `
                    <tr class="clickable-row filter_rows" id="${i}" data_gns='${JSON.stringify({'set1':unique_set1, 'set2':unique_set2, 'common':common})}'>
                      <td class="dt-center">${seg}</td>
                      <td class="dt-center">${i}</td>

                      <td class="narrow_col">${avg_set1}</td>
                      <td class="narrow_col">${avg_set2}</td>
                      <td class="narrow_col">${avg_diff}</td>

                      <td class="narrow_col">${unique_set1.length}</td>
                      <td class="narrow_col">${unique_set2.length}</td>
                      <td class="narrow_col">${common.length}</td>

                      <td class="narrow_col">${avg_freq_diff_sets}</td>

                      <td class="narrow_col">${set1_seq_cons_aa}</td>
                      <td class="narrow_col">${set2_seq_cons_aa}</td>
                      <td class="narrow_col">${set1_seq_cons_freq}</td>
                      <td class="narrow_col">${set2_seq_cons_freq}</td>
                      <td class="narrow_col">${diff_seq_cons_freq}</td>

                      <td class="narrow_col">${class_cons_aa}</td>
                      <td class="narrow_col">${class_cons_freq}</td>

                      <td class="narrow_col">${angles[0]}</td>
                      <td class="narrow_col">${angles[1]}</td>

                      <td class="narrow_col">${angles[2]}</td>
                      <td class="narrow_col">${angles[6]}</td>
                    </tr>`;
            // tbody.append(tr);
        });
        // insert natively for speed increase on Chrome
        tbody[0].innerHTML = tr_list;
    } else if (data['proteins'].length > 1) {

    } else {

    }

    $(selector + " .browser-table-3" + " .clickable-row").click(function() {
        var $this = $(this);
        var set_values = JSON.parse($this.attr('data_gns'));
        var column = $(this).children("td").eq(6);

        $('.popover').each(function() {
            $(this).hide();
        });

        //if not already initialized
        if (!column.data('bs.popover')) {
            set1 = set_values['set1'].join("<br>");
            set2 = set_values['set2'].join("<br>");
            common = set_values['common'].join("<br>");
            table = `<table>
                        <thead><tr><th style='width:60px'>Set 1</th><th style='width:60px'>Set 2</th><th style='width:60px'>Common</th></tr></thead>
                        <tbody><tr><td>${set1}</td><td>${set2}</td><td>${common}</td></tr></tbody></table>`;
            column.attr("data-content", table);
            column.popover({
                'container': selector + ' .dataTables_scrollBody',
                'placement': 'bottom',
                'animation': true,
                'html': true,
                'tabindex': '0'
            }).popover('show');
        } else {
            column.popover('show');
        }
    });

    enable_hover(table);
    console.timeEnd("RenderBrowser3");
}

function renderBrowser_4(data) {
    console.time("RenderBrowser4");
    var selector = $('ul#mode_nav').find('li.active').find('a').attr("href");
    var table = $(selector + " .browser-table-4");
    if ($.fn.DataTable.isDataTable(selector + " .browser-table-4")) {
        table.DataTable().destroy();
    }
    table.parent().html('<table class="browser-table-4 compact cell-border" style3="max-width: 2000px !important; width: 2000px !important" width2="2500px" style2="margin-left:0px" style="margin:0 auto"><thead></thead><tbody class="hidden"></tbody></table>');
    var table = $(selector + " .browser-table-4");
    // table.parent().before('<span><button type="button" onclick="filter_browser(this);" class="btn btn-xs btn-primary reset-selection">Filter</button></span>');
    var tbody = table.find('tbody');
    if (data['proteins2']) {
        var proteins_1 = data['proteins1'].length
        var proteins_2 = data['proteins2'].length
        var pdbs_1 = data['pdbs1'].length
        var pdbs_2 = data['pdbs2'].length
        thead = '<tr> \
                      <th colspan="2" class="skip"></th> \
                      <th colspan="2" class="selector" datatype="consensus_SS"></th> \
                      <th colspan="3" class="selector" datatype="consensus_freq"></th> \
                      <th colspan="3" class="selector" datatype="phi"></th> \
                      <th colspan="3" class="selector" datatype="psi"></th> \
                      <th colspan="3" class="selector" datatype="tau"></th> \
                      <th colspan="3" class="selector" datatype="theta"></th> \
                      <th colspan="3" class="selector" datatype="theta"></th> \
                      <th colspan="2" class="skip"></th> \
                      <th colspan="3" class="selector" datatype="conservation"></th> \
                      <th colspan="1" class="skip"></th> \
                      <th colspan="1" class="selector" datatype="conservation"></th> \
                  </tr> \
                  <tr> \
                          <th colspan="1" rowspan="2">Segment</th> \
                          <th colspan="1" rowspan="2">Positions</th> \
                          <th colspan="5" rowspan="1">Secondary structure</th> \
                          <th colspan="6" rowspan="1">Residue angles</th> \
                          <th colspan="9" rowspan="1">Helix turn angle</th> \
                          <th colspan="5" rowspan="1">Seq consensus</th> \
                          <th colspan="2" rowspan="1">Class seq consensus</th> \
                        </tr> \
                        <tr> \
                          <th colspan="2">Consensus SS</th> \
                          <th colspan="3">Frequency (%)</th> \
                          <th colspan="3">Phi (N(+1)-C-Ca-N)</th> \
                          <th colspan="3">Psi (C-Ca-N-C(-1))</th> \
                          <th colspan="3">Tau (Ca(+1)-Ca-Ca(-1)-Ca(-2)-)</th> \
                          <th colspan="3">Theta (Ca(+1)-Ca-Ca(-1))</th> \
                          <th colspan="3">Next Theta (Ca(+2)-Ca(+1)-Ca)</th> \
                          <th colspan="2">AA</th> \
                          <th colspan="3">Conservation (%)</th> \
                          <th colspan="1">AA</th> \
                          <th colspan="1">Cons (%)</th> \
                        </tr> \
                        <tr> \
                          <th class="dt-center"></th> \
                          <th class="dt-center"></th> \
                          <th class="narrow_col">Set 1<br></th> \
                          <th class="narrow_col">Set 2<br></th> \
                          <th class="narrow_col">Set 1<br></th> \
                          <th class="narrow_col">Set 2<br></th> \
                          <th class="narrow_col">Diff<br></th> \
                          <th class="narrow_col">Set 1<br></th> \
                          <th class="narrow_col">Set 2<br></th> \
                          <th class="narrow_col">Diff<br></th> \
                          <th class="narrow_col">Set 1<br></th> \
                          <th class="narrow_col">Set 2<br></th> \
                          <th class="narrow_col">Diff<br></th> \
                          <th class="narrow_col">Set 1<br></th> \
                          <th class="narrow_col">Set 2<br></th> \
                          <th class="narrow_col">Diff<br></th> \
                          <th class="narrow_col">Set 1<br></th> \
                          <th class="narrow_col">Set 2<br></th> \
                          <th class="narrow_col">Diff<br></th> \
                          <th class="narrow_col">Set 1<br></th> \
                          <th class="narrow_col">Set 2<br></th> \
                          <th class="narrow_col">Diff<br></th> \
                          <th class="narrow_col">Set 1<br></th> \
                          <th class="narrow_col">Set 2<br></th> \
                          <th class="narrow_col">Set 1<br></th> \
                          <th class="narrow_col">Set 2<br></th> \
                          <th class="narrow_col">Diff<br></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                        </tr>';
        table.find('thead').html(thead);
        tr_list = ''
        $.each(data['tab3'], function(i, v) {

            var seg = data['segm_lookup'][i];
            if (seg == 'ECL1' || seg == 'ECL2') return true;

            var angles1 = v['angles_set1'];
            var angles2 = v['angles_set2'];
            var angles_diff = v['angles'];
            // 0 'core_distance',
            // 1 'a_angle',
            // 2 'outer_angle',
            // 3 'tau',
            // 4 'phi',
            // 5 'psi',
            // 6 'sasa',
            // 7 'rsa',
            // 8 'theta',
            // 9 'hse'
            var set1_seq_cons_aa = v['set1_seq_cons'][0];
            var set2_seq_cons_aa = v['set2_seq_cons'][0];
            var set1_seq_cons_freq = Math.round(100 * v['set1_seq_cons'][1] / pdbs_1);
            var set2_seq_cons_freq = Math.round(100 * v['set2_seq_cons'][2] / pdbs_2);
            var diff_seq_cons_freq = Math.round((set1_seq_cons_freq - set2_seq_cons_freq));

            var class_cons_aa = v['class_cons'][0];
            var class_cons_freq = Math.round(100 * v['class_cons'][1]);

            all_angles_1 = data['all_angles'][i];
            ss_pos1_set1 = [];
            ss_pos1_set2 = [];
            pdbs = data['pdbs1'].concat(data['pdbs2']);
            pdbs.forEach(function(pdb){
                pdb_upper = pdb.toUpperCase();
                if (all_angles_1) {
                    let d1 = all_angles_1[pdb_upper];
                    if (d1.length) {
                        if (data['pdbs1'].includes(pdb)) {
                            ss_pos1_set1.push(d1[12]);
                        } else if (data['pdbs2'].includes(pdb)) {
                            ss_pos1_set2.push(d1[12]);
                        }
                    }
                }
            });

            dssp = [];
            [ss_pos1_set1,ss_pos1_set2].forEach(function(list){
                if (list.length) {
                    // make count map
                    c = new Map([...new Set(list)].map(
                            x => [x, list.filter(y => y === x).length]
                        ));
                    // make count with highest number first element
                    const mapSort = new Map([...c.entries()].sort((a, b) => b[1] - a[1]));
                    // get first element as [key,value]
                    most = mapSort.entries().next().value;
                    // calculate frequency
                    freq = most[1]/list.length;
                    most = most[0];
                } else {
                    freq = 0;
                    most = 'N/A';
                }
                dssp.push([most,freq]);
            })
            // console.table(dssp);
            // dssp_pos1 = '';
            // dssp_pos1_freq = '';
            // if (dssp[0][0]==dssp[1][0]){
            //     // if pos1 category same
            //     dssp_pos1=dssp[0][0];
            //     dssp_pos1_freq = Math.round(100*(dssp[0][1]-dssp[1][1]));
            // }


            // dssp_pos2 = '';
            // dssp_pos2_freq = '';
            // if (dssp[2][0]==dssp[3][0]){
            //     // if pos2 category same
            //     dssp_pos2=dssp[2][0];
            //     dssp_pos2_freq = Math.round(100*(dssp[2][1]-dssp[3][1]));
            // }

            tr = ''
            tr_list += `
                    <tr class="clickable-row filter_rows" id="${i}">
                      <td class="dt-center">${seg}</td>
                      <td class="dt-center">${i}</td>

                      <td class="narrow_col">${dssp[0][0]}</td>
                      <td class="narrow_col">${dssp[1][0]}</td>
                      <td class="narrow_col">${Math.round(100*dssp[0][1])}</td>
                      <td class="narrow_col">${Math.round(100*dssp[1][1])}</td>
                      <td class="narrow_col">${Math.round(100*(dssp[0][1]-dssp[1][1]))}</td>

                      <td class="narrow_col">${angles1[4]}</td>
                      <td class="narrow_col">${angles2[4]}</td>
                      <td class="narrow_col">${angles_diff[4]}</td>

                      <td class="narrow_col">${angles1[5]}</td>
                      <td class="narrow_col">${angles2[5]}</td>
                      <td class="narrow_col">${angles_diff[5]}</td>

                      <td class="narrow_col">${angles1[3]}</td>
                      <td class="narrow_col">${angles2[3]}</td>
                      <td class="narrow_col">${angles_diff[3]}</td>

                      <td class="narrow_col">${angles1[8]}</td>
                      <td class="narrow_col">${angles2[8]}</td>
                      <td class="narrow_col">${angles_diff[8]}</td>

                      <td class="narrow_col"></td>
                      <td class="narrow_col"></td>
                      <td class="narrow_col"></td>

                      <td class="narrow_col">${set1_seq_cons_aa}</td>
                      <td class="narrow_col">${set2_seq_cons_aa}</td>
                      <td class="narrow_col">${set1_seq_cons_freq}</td>
                      <td class="narrow_col">${set2_seq_cons_freq}</td>
                      <td class="narrow_col">${diff_seq_cons_freq}</td>

                      <td class="narrow_col">${class_cons_aa}</td>
                      <td class="narrow_col">${class_cons_freq}</td>

                    </tr>`;
            // tbody.append(tr);
        });
        // insert natively for speed increase on Chrome
        tbody[0].innerHTML = tr_list;


    } else if (data['proteins'].length > 1) {
        var proteins = data['proteins'].length
        var pdbs_count = data['pdbs'].length
        thead = '<tr> \
                      <th colspan="2" class="skip"></th> \
                      <th colspan="1" class="selector" datatype="consensus_SS"></th> \
                      <th colspan="1" class="selector" datatype="consensus_freq"></th> \
                      <th colspan="1" class="selector" datatype="phi"></th> \
                      <th colspan="1" class="selector" datatype="psi"></th> \
                      <th colspan="1" class="selector" datatype="tau"></th> \
                      <th colspan="1" class="selector" datatype="theta"></th> \
                      <th colspan="1" class="selector" datatype="theta"></th> \
                      <th colspan="1" class="skip"></th> \
                      <th colspan="1" class="selector" datatype="conservation"></th> \
                      <th colspan="1" class="skip"></th> \
                      <th colspan="1" class="selector" datatype="conservation"></th> \
                  </tr> \
                  <tr> \
                          <th colspan="1" rowspan="2">Segment</th> \
                          <th colspan="1" rowspan="2">Positions</th> \
                          <th colspan="2" rowspan="1">Secondary structure</th> \
                          <th colspan="2" rowspan="1">Residue angles</th> \
                          <th colspan="3" rowspan="1">Helix turn angle</th> \
                          <th colspan="2" rowspan="1">Seq consensus</th> \
                          <th colspan="2" rowspan="1">Class seq consensus</th> \
                        </tr> \
                        <tr> \
                          <th colspan="1">Consensus SS</th> \
                          <th colspan="1">Frequency (%)</th> \
                          <th colspan="1">Phi (N(+1)-C-Ca-N)</th> \
                          <th colspan="1">Psi (C-Ca-N-C(-1))</th> \
                          <th colspan="1">Tau (Ca(+1)-Ca-Ca(-1)-Ca(-2)-)</th> \
                          <th colspan="1">Theta (Ca(+1)-Ca-Ca(-1))</th> \
                          <th colspan="1">Next Theta (Ca(+2)-Ca(+1)-Ca)</th> \
                          <th colspan="1">AA</th> \
                          <th colspan="1">Conservation (%)</th> \
                          <th colspan="1">AA</th> \
                          <th colspan="1">Cons (%)</th> \
                        </tr> \
                        <tr> \
                          <th class="dt-center"></th> \
                          <th class="dt-center"></th> \
                          <th class="narrow_col">Set<br></th> \
                          <th class="narrow_col">Set<br></th> \
                          <th class="narrow_col">Set<br></th> \
                          <th class="narrow_col">Set<br></th> \
                          <th class="narrow_col">Set<br></th> \
                          <th class="narrow_col">Set<br></th> \
                          <th class="narrow_col">Set<br></th> \
                          <th class="narrow_col">Set<br></th> \
                          <th class="narrow_col">Set<br></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                        </tr>';
        table.find('thead').html(thead);
        tr_list = ''
        $.each(data['tab3'], function(i, v) {

            var seg = data['segm_lookup'][i];
            if (seg == 'ECL1' || seg == 'ECL2') return true;

            var angles = v['angles_set'];
            // 0 'core_distance',
            // 1 'a_angle',
            // 2 'outer_angle',
            // 3 'tau',
            // 4 'phi',
            // 5 'psi',
            // 6 'sasa',
            // 7 'rsa',
            // 8 'theta',
            // 9 'hse'
            var set_seq_cons_aa = v['set_seq_cons'][0];
            var set_seq_cons_freq = Math.round(100 * v['set_seq_cons'][1] / pdbs_count);

            var class_cons_aa = v['class_cons'][0];
            var class_cons_freq = Math.round(100 * v['class_cons'][1]);

            all_angles_1 = data['all_angles'][i];
            ss_pos1_set = [];
            pdbs = data['pdbs'];
            pdbs.forEach(function(pdb){
                pdb_upper = pdb.toUpperCase();
                if (all_angles_1) {
                    let d1 = all_angles_1[pdb_upper];
                    if (d1.length) {
                        ss_pos1_set.push(d1[12]);
                    }
                }
            });

            dssp = [];
            [ss_pos1_set].forEach(function(list){
                if (list.length) {
                    // make count map
                    c = new Map([...new Set(list)].map(
                            x => [x, list.filter(y => y === x).length]
                        ));
                    // make count with highest number first element
                    const mapSort = new Map([...c.entries()].sort((a, b) => b[1] - a[1]));
                    // get first element as [key,value]
                    most = mapSort.entries().next().value;
                    // calculate frequency
                    freq = most[1]/list.length;
                    most = most[0];
                } else {
                    freq = 0;
                    most = 'N/A';
                }
                dssp.push([most,freq]);
            })

            tr = ''
            tr_list += `
                    <tr class="clickable-row filter_rows" id="${i}">
                      <td class="dt-center">${seg}</td>
                      <td class="dt-center">${i}</td>

                      <td class="narrow_col">${dssp[0][0]}</td>
                      <td class="narrow_col">${Math.round(100*dssp[0][1])}</td>

                      <td class="narrow_col">${angles[4]}</td>

                      <td class="narrow_col">${angles[5]}</td>

                      <td class="narrow_col">${angles[3]}</td>

                      <td class="narrow_col">${angles[8]}</td>

                      <td class="narrow_col"></td>

                      <td class="narrow_col">${set_seq_cons_aa}</td>
                      <td class="narrow_col">${set_seq_cons_freq}</td>

                      <td class="narrow_col">${class_cons_aa}</td>
                      <td class="narrow_col">${class_cons_freq}</td>

                    </tr>`;
        });
        // insert natively for speed increase on Chrome
        tbody[0].innerHTML = tr_list;

    } else {
        //var proteins = data['proteins'].length
        //var pdbs = data['pdbs'].length
        thead = '<tr> \
                      <th colspan="2" class="skip"></th> \
                      <th colspan="1" class="selector" datatype="consensus_SS"></th> \
                      <th colspan="1" class="selector" datatype="phi"></th> \
                      <th colspan="1" class="selector" datatype="psi"></th> \
                      <th colspan="1" class="selector" datatype="tau"></th> \
                      <th colspan="1" class="selector" datatype="theta"></th> \
                      <th colspan="1" class="selector" datatype="theta"></th> \
                      <th colspan="1" class="skip"></th> \
                      <th colspan="1" class="skip"></th> \
                      <th colspan="1" class="selector" datatype="conservation"></th> \
                  </tr> \
                  <tr> \
                          <th colspan="1" rowspan="2">Segment</th> \
                          <th colspan="1" rowspan="2">Positions</th> \
                          <th colspan="1" rowspan="1">Secondary structure</th> \
                          <th colspan="2" rowspan="1">Residue angles</th> \
                          <th colspan="3" rowspan="1">Helix turn angle</th> \
                          <th colspan="1" rowspan="1">Seq</th> \
                          <th colspan="2" rowspan="1">Class seq consensus</th> \
                        </tr> \
                        <tr> \
                          <th colspan="1">SS</th> \
                          <th colspan="1">Phi (N(+1)-C-Ca-N)</th> \
                          <th colspan="1">Psi (C-Ca-N-C(-1))</th> \
                          <th colspan="1">Tau (Ca(+1)-Ca-Ca(-1)-Ca(-2)-)</th> \
                          <th colspan="1">Theta (Ca(+1)-Ca-Ca(-1))</th> \
                          <th colspan="1">Next Theta (Ca(+2)-Ca(+1)-Ca)</th> \
                          <th colspan="1">AA</th> \
                          <th colspan="1">AA</th> \
                          <th colspan="1">Cons (%)</th> \
                        </tr> \
                        <tr> \
                          <th class="dt-center"></th> \
                          <th class="dt-center"></th> \
                          <th class="narrow_col">Set<br></th> \
                          <th class="narrow_col">Set<br></th> \
                          <th class="narrow_col">Set<br></th> \
                          <th class="narrow_col">Set<br></th> \
                          <th class="narrow_col">Set<br></th> \
                          <th class="narrow_col">Set<br></th> \
                          <th class="narrow_col">Set<br></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                        </tr>';
        table.find('thead').html(thead);
        tr_list = ''
        $.each(data['tab3'], function(i, v) {

            var seg = data['segm_lookup'][i];
            if (seg == 'ECL1' || seg == 'ECL2') return true;

            var angles = v['angles_set'];
            // 0 'core_distance',
            // 1 'a_angle',
            // 2 'outer_angle',
            // 3 'tau',
            // 4 'phi',
            // 5 'psi',
            // 6 'sasa',
            // 7 'rsa',
            // 8 'theta',
            // 9 'hse'
            var set_seq_cons_aa = v['set_seq_cons'][0];
            var class_cons_aa = v['class_cons'][0];
            var class_cons_freq = Math.round(100 * v['class_cons'][1]);

            all_angles_1 = data['all_angles'][i];

            if (all_angles_1) {
                all_angles_1 = all_angles_1[data['pdbs'][0].toUpperCase()][12];
            } else {
                all_angles_1 = '';
            }

            tr = ''
            tr_list += `
                    <tr class="clickable-row filter_rows" id="${i}">
                      <td class="dt-center">${seg}</td>
                      <td class="dt-center">${i}</td>

                      <td class="narrow_col">${all_angles_1}</td>

                      <td class="narrow_col">${angles[4]}</td>

                      <td class="narrow_col">${angles[5]}</td>

                      <td class="narrow_col">${angles[3]}</td>

                      <td class="narrow_col">${angles[8]}</td>

                      <td class="narrow_col"></td>

                      <td class="narrow_col">${set_seq_cons_aa}</td>

                      <td class="narrow_col">${class_cons_aa}</td>
                      <td class="narrow_col">${class_cons_freq}</td>

                    </tr>`;
        });
        // insert natively for speed increase on Chrome
        tbody[0].innerHTML = tr_list;

    }

    enable_hover(table);
    console.timeEnd("RenderBrowser4");
}

function renderBrowser_5(data) {
    console.time("RenderBrowser5");
    var selector = $('ul#mode_nav').find('li.active').find('a').attr("href");
    var table = $(selector + " .browser-table-5");
    if ($.fn.DataTable.isDataTable(selector + " .browser-table-5")) {
        table.DataTable().destroy();
    }
    table.parent().html('<table class="browser-table-5 compact cell-border"><thead></thead><tbody class="hidden"></tbody></table>');
    var table = $(selector + " .browser-table-5");
    // table.parent().before('<span><button type="button" onclick="filter_browser(this);" class="btn btn-xs btn-primary reset-selection">Filter</button></span>');
    var tbody = table.find('tbody');

    var thead;
    if (data['proteins2']) {
      thead = '<tr> \
                    <th colspan="2" class="skip"></th> \
                    <th colspan="1" class="selector" datatype="core_distance_diff"></th> \
                    <th colspan="1" class="selector" datatype="rotation_diff"></th> \
                    <th colspan="1" class="selector" datatype="HSE_diff"></th> \
                </tr>';
    } else {
      thead = '<tr> \
                    <th colspan="2" class="skip"></th> \
                    <th colspan="1" class="selector" datatype="core_distance"></th> \
                    <th colspan="1" class="selector" datatype="rotation"></th> \
                    <th colspan="1" class="selector" datatype="HSE"></th> \
                </tr>';
    }

    thead += '<tr> \
                      <th colspan="1" rowspan="2">Segment</th> \
                      <th colspan="1" rowspan="2">Positions</th> \
                      <th colspan="2">Backbone movement (Ca-7TM axis)</th> \
                      <th colspan="1" rowspan="2">Ca half-sphere exposure</th> \
                    </tr> \
                    <tr> \
                      <th colspan="1">Distance*</th> \
                      <th colspan="1">Rotation</th> \
                    </tr> \
                    <tr> \
                      <th class="dt-center"></th> \
                      <th class="dt-center"></th> \
                      <th class="narrow_col"></th> \
                      <th class="narrow_col"></th> \
                      <th class="narrow_col"></th> \
                    </tr>';
    table.find('thead').html(thead);
    tr_list = ''
    $.each(data['tab3'], function(i, v) {

        // console.log(i,v);
        var seg = data['segm_lookup'][i];
        var angles = v['angles'];
        // 0 'core_distance',
        // 1 'a_angle',
        // 2 'outer_angle',
        // 3 'tau',
        // 4 'phi',
        // 5 'psi',
        // 6 'sasa',
        // 7 'rsa',
        // 8 'theta',
        // 9 'hse'

        tr = ''
        tr_list += `
                <tr class="clickable-row filter_rows" id="${i}">
                  <td class="dt-center">${seg}</td>
                  <td class="dt-center">${i}</td>
                  <td class="narrow_col">${angles[0]}</td>
                  <td class="narrow_col">${angles[1]}</td>
                  <td class="narrow_col">${angles[9]}</td>
                </tr>`;
        // tbody.append(tr);
    });
    // insert natively for speed increase on Chrome
    tbody[0].innerHTML = tr_list;

    enable_hover(table)
    console.timeEnd("RenderBrowser5");
}

function gray_scale_table(table) {
    console.time('Greyscale');
    // Create grey scale of values.
    var cols = []
    for (let [i, row] of [...table.find("tbody")[0].rows].entries()) {
        for (let [j, cell] of [...row.cells].entries()) {
            cols[j] = cols[j] || [];
            cols[j].push(cell.innerText)
        }
    }
    maxmin = [];
    cols.forEach(function(col, index) {
        var max = Math.max.apply(null, col);
        var min = Math.min.apply(null, col);
        maxmin.push([max, min]);
    });
    // console.time('Greyscale cells');
    var cell_count = 0;
    for (let [i, row] of [...table.find("tbody")[0].rows].entries()) {
        for (let [j, cell] of [...row.cells].entries()) {
            c_maxmin = maxmin[j];
            value = parseFloat(cell.innerText);
            if (!(isNaN(value) || isNaN(c_maxmin[0]) || isNaN(c_maxmin[1]))) {
                // console.log(`[${i},${j}] = ${cell.innerText} ${c_maxmin}`);
                scale = 1 - (value - c_maxmin[1]) / (c_maxmin[0] - c_maxmin[1]);
                frequency = 0.5 - scale * .5;
                color_255 = Math.round(255 - frequency * 255);
                var rgb = {
                    r: color_255,
                    g: color_255,
                    b: color_255
                };
                var hex = rgb2hex(rgb.r, rgb.g, rgb.b);
                cell.setAttribute("bgcolor", hex);
                cell_count++;
            }
        }
    }
    // console.timeEnd('Greyscale cells');
    console.log(cell_count, 'cells greyscaled');
    console.timeEnd('Greyscale');
}

var currentHover = -1;
function enable_hover(table){
    table[0].children[0].addEventListener("mouseover", function(e){
      var th = e.target
      while (th.nodeName != "TH") {
        th = th.parentNode
      }
      var columnNumber = $(th).cellPos().left;

      // Get correct selector cell
      var selectorHeader = th.parentNode.parentNode.children[0]
      var selector = selectorHeader.children[0]
      var columnSelector = 0
      for (var i = 0; i < selectorHeader.children.length; i++) {
        if ($(selectorHeader.children[i]).cellPos().left > columnNumber)
          break
        selector = selectorHeader.children[i]
        columnSelector = $(selectorHeader.children[i]).cellPos().left
      }

      if (currentHover != columnSelector && selector.className!="skip" && !selector.className.includes("keep")) {
        // other variables
        var tableNumber = th.parentNode.parentNode.parentNode.className.split(" ")[0]
        var tableNumber = tableNumber.substr(-1)

        // grab graph options
        var plots = $('.main_option:visible').find(".plot-container");
        for (var i = 0; i < plots.length; i++){
          var plotType = plots[i].id

          var button = document.createElement("span")
          button.className = "glyphicon glyphicon-stats toggle"+i
          selector.appendChild(button)

          var found = true;
          if (selector.className=="pairselector") {
              if (plotType.startsWith("heatmapcontainer") || plotType.startsWith("flareplot") || plotType.startsWith("boxplot")) {
                button.addEventListener("click", (function(a, b, c, d){ return function(){colorByData(a, b, c, d);}})(plotType, tableNumber, columnSelector, selector.getAttribute("datatype")))
              } else {
                found = false;
              }
          } else if (selector.className=="selector") {
            if (plotType.startsWith("ngl") || plotType.startsWith("snakeplot")) {
              button.addEventListener("click", (function(a, b, c, d){ return function(){colorByData(a, b, c, d);}})(plotType, tableNumber, columnSelector, selector.getAttribute("datatype")))
            } else {
              found = false;
            }
          }

          if (found){
            button.addEventListener("click", function(e){
              var targetClasses = e.target.className.split(" ")
              var object = $(e.target)
              if (!object.hasClass("red")){
                // Remove toggle and keep from other header if present
                $(".glyphicon-stats.red."+targetClasses[targetClasses.length -1]).each( function(i, other){
                    $(other).removeClass("red")
                    // clean header
                    if ($(other).parent().find(".red").length == 0) {
                      $(other).parent().removeClass("keep")
                    }
                })

                // Keep header enabled
                if (!object.parent().hasClass("keep"))
                  object.parent().addClass("keep");

                // Toggle icon color
                object.addClass("red")
              }
            });
          } else {
            // Grayout button if not available
            button.className = button.className + " gray"
          }
        }

        currentHover = columnSelector;
      }
    });

    table[0].children[0].addEventListener("mouseout", function(e){
      classes = e.target.className
      if (!(classes.includes("glyphicon") || classes.includes("selector") || classes.includes("pairselector"))) {
        clearGraphHeader(e)
      }
    });

    header = table[0].children[0].children[0];
    for (var i = 0; i < header.children.length; i++){
      $(header.children[i]).mouseleave( clearGraphHeader );
    }
}

function clearGraphHeader(e){
  // clear selector header on mouse out
  var header = e.target
  while (header.nodeName != "THEAD") {
    header = header.parentNode
  }

  // cleanup with smarter class selector
  header = header.children[0];
  for (var i = 0; i < header.children.length; i++){
      if (header.children[i].innerHTML.length > 0 && !header.children[i].className.includes("keep")){
          header.children[i].innerHTML = ""
      }
  }
  currentHover = -1;
}

function colorByData(mode, tableNumber, columnNumber, type) {
    var defaultColor = "#666";

    // TODO: make a switch for the different data tables (or handle in the click handler function?)
    /*     switch (mode) {
             case "single-crystal-group":
               break;
             case "single-crystal":
               break;
             case "two-crystal-groups":
               break;
         }*/

    // Grab data from table
    var rows = []
    const selector = "#" + $('.main_option:visible').attr('id');
    const analys_mode = selector.replace('-tab', '');
    if (analys_mode == "#single-crystal" && tableNumber==1) {
      rows = getDateFromTable(tableNumber, [2, columnNumber].flat());
    } else {
      rows = getDateFromTable(tableNumber, [1, columnNumber].flat());
    }

    // Residue positions
    var residue_positions = getColumn(rows, 0);

    // If residue pairs -> split and concatenate
    if (residue_positions.length == 0)
      return;

    if (residue_positions[0].indexOf("-")>=1){
      residue_positions1 = residue_positions.map(function(e){return e.split("-")[0]});
      residue_positions2 = residue_positions.map(function(e){return e.split("-")[1]});
      residue_positions = residue_positions1.concat(residue_positions2)
    }

    // Filter NaNs
    var residue_values = getColumn(rows, 1)
    if (type!="consensus_SS"){
      residue_values = residue_values.map(function(e){ return parseInt(e);})
      if (Array.isArray(columnNumber) && columnNumber.length > 1)
         residue_values = residue_values.concat(getColumn(rows, 2).map(function(e){ return parseInt(e);}));

      residue_positions = residue_positions.filter( function(value, index){ return !isNaN(residue_values[index]); })
      residue_values = residue_values.filter( function(value, index){ return !isNaN(residue_values[index]); })
    } else {
      // remove positions with no data
      residue_positions = residue_positions.filter( function(value, index){ return !(residue_values==""); })
      residue_values = residue_values.filter( function(value, index){ return !(residue_values==""); })
    }

    // Identify range
    var valMax = Math.max(...residue_values)
    var valMin = Math.min(...residue_values)
    var palette = "rwb"
    if (type in dataType){
      valMax = dataType[type][2]
      valMin = dataType[type][1]
      palette = dataType[type][3]
    } else {
      console.log("TYPE not found: "+type)
    }

    var data_colors = []
    for (var i = 0; i < residue_positions.length; i++) {
      // get color
      var newColor = defaultColor;
      if (palette=="ss_coloring") { // coloring for SS
          if (residue_values[i] in SS_COLORING)
            newColor = SS_COLORING[residue_values[i]]
      } else {
        if (valMin < 0)
          newColor = numberToColorGradient(residue_values[i], valMax, palette, neg_and_pos = true)
        else
          newColor = numberToColorGradient(residue_values[i], valMax, palette, neg_and_pos = false)
      }

      data_colors.push(newColor)
    }

    // Toggle for different plots
    if (mode.startsWith("ngl")) { // 3D
      colorNGLByData(mode.replace("ngl-",""), residue_positions, data_colors, defaultColor);
    } else if (mode.startsWith("heatmapcontainer")) { // heatmap
      //somefunction(mode, residue_positions, data_colors, defaultColor);
    } else if (mode.startsWith("flareplot")) {  // flareplot
      //somefunction(mode, residue_positions, data_colors, defaultColor);
    } else if (mode.startsWith("boxplot")) { // Box-plot (not available right now)
      //somefunction(mode, residue_positions, data_colors, defaultColor);
    } else if (mode.startsWith("snakeplot")) { // Snakeplot
      //somefunction(mode, residue_positions, data_colors, defaultColor);
    } else {
      // missing plot type?
      console.log("Missing plot type")
    }
}

// Reference values for table values
// entry: display (boolean), min, max, coloring
// Values for single residue (skipping the pairs)
dataType = {}
dataType["no_viewing"]     = [false, 0, 0, ""]
dataType["core_distance"]  = [true, 0, 15, "wb"]
dataType["rotation"]       = [true, 0, 180, "wb"]
dataType["rotamer"]        = [true, 0, 180, "wb"]
dataType["consensus_freq"] = [true, 0, 100, "wb"]
dataType["consensus_SS"]   = [true, 0, 100, "ss_coloring"]
dataType["SASA"]           = [true, 0, 200, "wb"]
dataType["RSA"]            = [true, 0, 100, "wb"]
dataType["HSE"]            = [true, 0, 20, "wb"]
dataType["conservation"]   = [true, 0, 100, "wb"]
dataType["presence"]       = [true, 0, 100, "wb"]
dataType["phi"]            = [true, 0, 180, "wb"]
dataType["psi"]            = [true, 0, 180, "wb"]
dataType["tau"]            = [true, 0, 180, "wb"]
dataType["theta"]          = [true, 0, 180, "wb"]

// Same for differences between groups
// TODO: finetune to relative variation per data type
dataType["core_distance_diff"] = [true, -5, 5, "rwb"]
dataType["rotation_diff"]      = [true, -45, 45, "rwb"]
dataType["rotamer_diff"]       = [true, -45, 45, "rwb"]
dataType["ss_freq_diff"]       = [true, -50, 50, "rwb"]
dataType["SASA_diff"]          = [true, -50, 50, "rwb"]
dataType["RSA_diff"]           = [true, -25, 25, "rwb"]
dataType["HSE_diff"]           = [true, -10, 10, "rwb"]
//dataType["conservation"]  = [true, 0, 100, "rwb"]
dataType["phi_diff"]           = [true, -45, 45, "rwb"]
dataType["psi_diff"]           = [true, -45, 45, "rwb"]
dataType["tau_diff"]           = [true, -45, 45, "rwb"]
dataType["theta_diff"]         = [true, -45, 45, "rwb"]

// Secondary structure coloringData
var SS_COLORING = {
  "H": "#FF0000", // alpha helix
  "G": "#FF00FF", // 3-10 helix
  "I": "#00FF00", // pi helix
  "B": "#ffd700", // isolated bridge
  "E": "#00a8ff", // strand/extended conformation
  "T": "#ff5700", // turn
  "S": "#00ffff", // bend (unique to DSSP)
  "C": "#5d8aa8", // coil/other (STRIDE)
  "-": "#5d8aa8", // coil/other (DSSP)
}

// TODO: add smart handling of minimum values (now based on max)
function numberToColorGradient(value, max, palette, neg_and_pos = false) {
    if (neg_and_pos) {
      value = value + max
      max = max*2
    }

    if (value > max)
      value = max
    if (value < 0)
      value = 0

    switch(palette){
        case "rwb": // red-white-blue
          return colorGradient(value/max, {red:255, green:0, blue: 0}, {red:255, green:255, blue: 255}, {red:0, green:0, blue: 255})
          break;
        case "bwr": // blue-white-red
          return colorGradient(value/max, {red:0, green:0, blue: 255}, {red:255, green:255, blue: 255}, {red:255, green:0, blue: 0})
          break;
        case "ryg": // red-yellow-green
          return colorGradient(value/max, {red:255, green:0, blue: 0}, {red:0, green:255, blue: 0}, {red:0, green:255, blue: 0})
          break;
        case "gyr": // green-yellow-red
          return colorGradient(value/max, {red:255, green:0, blue: 0}, {red:255, green:255, blue: 0}, {red:0, green:255, blue: 0})
          break;
        case "rgb":
          return colorGradient(value/max, {red:255, green:0, blue: 0}, {red:255, green:255, blue: 255}, {red:0, green:0, blue: 255})
          break;
        case "wr": // white-red
          return colorGradient(value/max, {red:255, green:255, blue: 255}, {red:255, green:0, blue: 0})
          break;
        case "wg": // white-green
          return colorGradient(value/max, {red:255, green:255, blue: 255}, {red:0, green:255, blue: 0})
          break;
        case "wb": // white-blue
          return colorGradient(value/max, {red:255, green:255, blue: 255}, {red:0, green:0, blue: 255})
          break;
        case "rb": // red-blue
          return colorGradient(value/max, {red:255, green:0, blue: 0}, {red:0, green:0, blue: 255})
          break;
        // ADDON if you're missing gradient values
        case "br": // blue-red
        default:
          return colorGradient(value/max, {red:0, green:0, blue: 255}, {red:255, green:0, blue: 0})
          break;
    }
}

function colorGradient(fadeFraction, rgbColor1, rgbColor2, rgbColor3) {
    var color1 = rgbColor1;
    var color2 = rgbColor2;
    var fade = fadeFraction;

    // Do we have 3 colors for the gradient? Need to adjust the params.
    if (rgbColor3) {
      fade = fade * 2;

      // Find which interval to use and adjust the fade percentage
      if (fade >= 1) {
        fade -= 1;
        color1 = rgbColor2;
        color2 = rgbColor3;
      }
    }

    var diffRed = color2.red - color1.red;
    var diffGreen = color2.green - color1.green;
    var diffBlue = color2.blue - color1.blue;

    var gradient = {
      red: parseInt(Math.floor(color1.red + (diffRed * fade)), 10),
      green: parseInt(Math.floor(color1.green + (diffGreen * fade)), 10),
      blue: parseInt(Math.floor(color1.blue + (diffBlue * fade)), 10),
    };

    return rgb2hexCG(gradient.red.toString(16),gradient.green.toString(16),gradient.blue.toString(16));
}

function rgb2hexCG(r,g,b) {
    if (r.length == 1)
        r = '0' + r;
    if (g.length == 1)
        g = '0' + g;
    if (b.length == 1)
        b = '0' + b;

    return '#' + r + g + b;
}

/*function enable_3Dclick(table){
  for (header in table[0].children[0].children[1].children){
    var th = table[0].children[0].children[1].children[header]
    if (typeof th === 'object')
      th.addEventListener("click", function(e){

        // filter keys for current mode (single/single_group/two_sets)
        const analys_mode = $('.main_option:visible').attr('id').replace('-tab', '');
        var cmode = "single_"
        if (analys_mode=="two-crystal-groups")
          cmode = "two_sets_"
        else if (analys_mode=="single-crystal-group")
          cmode = "single_group_"

        // BUG: single_ and single_group both match the single_ string
        var viewers = Object.keys(stage).filter(function(x){ return x.startsWith(cmode)})
        if (viewers.length > 0) {
          var mode = viewers[0];
          // TODO: select which 3D view if more than one

          // Table data
          var th = e.target
          var tableNumber = th.parentNode.parentNode.parentNode.className.split(" ")[0];
          var tableNumber = tableNumber.substr(-1)
          var columnNumber = $(th).cellPos().left;

          // Color 3D viewer
          if ( th.colSpan == 3 ){
            // Toggle between group 1/2 values and group differences

          } else if (th.colSpan == 2 ){
            colorByData(mode, tableNumber, [columnNumber, columnNumber+1])
          } else {
            colorByData(mode, tableNumber, columnNumber)
          }
        }
      })
  }
}*/

/*  cellPos jQuery plugin
    ---------------------
    Get visual position of cell in HTML table (or its block like thead).
    Return value is object with "top" and "left" properties set to row and column index of top-left cell corner.
    Example of use:
        $("#myTable tbody td").each(function(){
            $(this).text( $(this).cellPos().top +", "+ $(this).cellPos().left );
        });
*/
(function($){
    /* scan individual table and set "cellPos" data in the form { left: x-coord, top: y-coord } */
    function scanTable( $table ) {
        var m = [];
        $table.children( "tr" ).each( function( y, row ) {
            $( row ).children( "td, th" ).each( function( x, cell ) {
                var $cell = $( cell ),
                    cspan = $cell.attr( "colspan" ) | 0,
                    rspan = $cell.attr( "rowspan" ) | 0,
                    tx, ty;
                cspan = cspan ? cspan : 1;
                rspan = rspan ? rspan : 1;
                for( ; m[y] && m[y][x]; ++x );  //skip already occupied cells in current row
                for( tx = x; tx < x + cspan; ++tx ) {  //mark matrix elements occupied by current cell with true
                    for( ty = y; ty < y + rspan; ++ty ) {
                        if( !m[ty] ) {  //fill missing rows
                            m[ty] = [];
                        }
                        m[ty][tx] = true;
                    }
                }
                var pos = { top: y, left: x };
                $cell.data( "cellPos", pos );
            } );
        } );
    };

    /* plugin */
    $.fn.cellPos = function( rescan ) {
        var $cell = this.first(),
            pos = $cell.data( "cellPos" );
        if( !pos || rescan ) {
            var $table = $cell.closest( "table, thead, tbody, tfoot" );
            scanTable( $table );
        }
        pos = $cell.data( "cellPos" );
        return pos;
    }
})(jQuery);

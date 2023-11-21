var filtered_gn_pairs = [];
var filtered_cluster_groups = [];
var filtered_gns = [];
var filtered_gns_abs_diff_values = {};
var filtered_cluster_groups_set = {};
function filter_browser() {
    old_filtered_gn_pairs = filtered_gn_pairs;
    filtered_gn_pairs = [];
    filtered_gns = [];
    filtered_gns_abs_diff_values = {};
    filtered_gns_presence = {}; // Whether a filtered gn takes part in contacts only in set1,set2 or in both
    pos_contacts_count = {};
    filtered_cluster_groups = [];
    const selector = "#" + $('.main_option:visible').attr('id');
    const analys_mode = selector.replace('-tab', '');

    if ($.fn.DataTable.isDataTable(".browser-table-1:visible")) {
        var table = $(".browser-table-1:visible").DataTable();
        table.rows({
            filter: 'applied'
        }).data().each(function(i) {
            filtered_gn_pairs.push(i['DT_RowId'])
            gns = separatePair(i['DT_RowId']);
            filtered_gns.push(gns[0]);
            filtered_gns.push(gns[1]);

            if (analys_mode == "#two-crystal-groups") {
                if (!(gns[0] in filtered_gns_abs_diff_values)) filtered_gns_abs_diff_values[gns[0]] = [];
                if (!(gns[1] in filtered_gns_abs_diff_values)) filtered_gns_abs_diff_values[gns[1]] = [];
                // BEWARE this 4th index can change if the column changes.. only on relevant in 2 group
                diff_value = i['2'] - i['3'];
                filtered_gns_abs_diff_values[gns[0]].push(diff_value);
                filtered_gns_abs_diff_values[gns[1]].push(diff_value);
            }

            // see if there is a key for gns1
            if (!(gns[0] in pos_contacts_count)) pos_contacts_count[gns[0]] = 0;
            pos_contacts_count[gns[0]] += 1;
            // see if there is a key for gns2
            if (!(gns[1] in pos_contacts_count)) pos_contacts_count[gns[1]] = 0;
            pos_contacts_count[gns[1]] += 1;

            // Track network groups

            test1 = filtered_cluster_groups.filter(l => l.includes(gns[0]));
            test2 = filtered_cluster_groups.filter(l => l.includes(gns[1]));
            if (!test1.length && !test2.length) {
                filtered_cluster_groups.push([gns[0], gns[1]]);
            } else if (test1.length && !test2.length) {
                i1 = filtered_cluster_groups.indexOf(test1[0])
                filtered_cluster_groups[i1].push(gns[1]);
            } else if (!test1.length && test2.length) {
                i2 = filtered_cluster_groups.indexOf(test2[0])
                filtered_cluster_groups[i2].push(gns[0]);
            } else if (test1.length && test2.length) {
                i1 = filtered_cluster_groups.indexOf(test1[0])
                i2 = filtered_cluster_groups.indexOf(test2[0])
                //i1 = filtered_cluster_groups.indexOfForArrays(test1[0]);
                if (i1!=i2) {
                    filtered_cluster_groups[i1] = test1[0].concat(test2[0])
                    filtered_cluster_groups.splice(i2, 1);
                }
            }

        })
        filtered_gns_presence = {};
        filtered_cluster_groups_set = {};
        if (analys_mode == "#two-crystal-groups") {
            $.each(filtered_gns_abs_diff_values, function (i, v) {
                // Go through all the diff values. If both negative and positive diff numbers exist
                // then label position as "both". Otherwise the correct set. This gives information whether
                // the position is only participating in interactions in one set or the other..
                let max = Math.max.apply(null, v);
                let min = Math.min.apply(null, v);
                var in_set_1 = true ? max >= 0 : false;
                var in_set_2 = true ? min <= 0 : false;
                if (in_set_1 && in_set_2) {
                    filtered_gns_presence[i] = 0.5; //both, middle
                } else if (in_set_1) {
                    filtered_gns_presence[i] = 0; //set1
                } else if (in_set_2) {
                    filtered_gns_presence[i] = 1; //set2
                }

            })

            $.each(filtered_cluster_groups, function (i, gns) {

                // console.log('filter id', i);
                var sum = 0;
                for (var ii = 0; ii < gns.length; ii++){
                    // console.log(ii, gns[ii],filtered_gns_presence[gns[ii]]);
                    sum += filtered_gns_presence[gns[ii]];
                }
                var avg = sum / gns.length;
                group_set = "both";
                if (avg == 0) {
                    group_set = "set1";
                } else if (avg == 1) {
                    group_set = "set2";
                }
                // console.log('filter id', i, avg, group_set);
                filtered_cluster_groups_set[i] = group_set;
            })
        }

        console.time('Update network');
        if (old_filtered_gn_pairs.sort().join(',') !== filtered_gn_pairs.sort().join(',')) {
            // only update this if there are new filtered things..
            rowIndexes = table.rows({ filter: 'applied' }).indexes();
            table.rows({
                filter: 'applied'
            }).data().each(function (i, index) {
                rowindex = rowIndexes[index];
                gns = separatePair(i['DT_RowId']);

                network_group = filtered_cluster_groups.filter(l => l.includes(gns[0]));
                network_group_id = filtered_cluster_groups.indexOf(network_group[0])

                if (analys_mode == "#two-crystal-groups") {
                    column_ids = [5, 6, 7];
                } else if (analys_mode == "#single-crystal-group") {
                    column_ids = [3, 4, 5];

                } else {
                    column_ids = [3, 4, 5];

                }
                table.cell({ row: rowindex, column: column_ids[0] }).data(pos_contacts_count[gns[0]]);
                table.cell({ row: rowindex, column: column_ids[1] }).data(pos_contacts_count[gns[1]]);
                table.cell({ row: rowindex, column: column_ids[2] }).data("#"+(network_group_id+1));

            })
        }

        console.timeEnd('Update network');
    } else {
        console.log('filter_browser requested, but tab-1 not visible.');
        console.log('reset filtered.')
    }

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

jQuery.extend( jQuery.fn.dataTableExt.oSort, {
    "non-empty-string-asc": function (str1, str2) {
        if(str1 == "")
            return 1;
        if(str2 == "")
            return -1;
        if ($.isNumeric(str1))
            str1 = parseFloat(str1);
        if ($.isNumeric(str2))
            str2 = parseFloat(str2);
        return ((str1 < str2) ? -1 : ((str1 > str2) ? 1 : 0));
    },

    "non-empty-string-desc": function (str1, str2) {
        if(str1 == "")
            return 1;
        if(str2 == "")
            return -1;
        if ($.isNumeric(str1))
            str1 = parseFloat(str1);
        if ($.isNumeric(str2))
            str2 = parseFloat(str2);
        return ((str1 < str2) ? 1 : ((str1 > str2) ? -1 : 0));
    }
} );

function renderDataTablesYadcf(element) {

    console.time("renderDataTablesYadcf");

    const tab_number = element.split("-").reverse()[2];
    const selector = "#" + $('.main_option:visible').attr('id');
    const analys_mode = selector.replace('-tab', '');
    var table = $(selector + " .browser-table-" + tab_number);
    var heading = $(selector + " .tab-content .panel-title:visible");
    if (analys_mode == "#two-crystal-groups") {
        heading.find(".abs_button").remove();
        heading.append(' <button type="button"  onclick="make_abs_values(this,\''+selector + " .browser-table-" + tab_number+'\');" class="btn btn-primary btn-xs abs_button" changed=0>Change negative to absolute values</button>');
        //heading.addClass("button_added");
    }
    // If table is without tbody, then do not init further.
    if (!(table.find("thead").length)) {
        console.timeEnd("renderDataTablesYadcf");
        $(".main_loading_overlay").hide();
        return
    }
    // Do not re init the table.
    if ($.fn.DataTable.isDataTable(selector + " .browser-table-" + tab_number)) {
        console.timeEnd("renderDataTablesYadcf");
        var btable = $(selector + " .browser-table-" + tab_number).DataTable();
        btable.columns.adjust().draw();
        $(".main_loading_overlay").hide();
        return
    }

    // $(".main_loading_overlay").show(0);
    var buttonCommon = {}
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
                "bInfo": true,
                "fnInfoCallback": function (oSettings, iStart, iEnd, iMax, iTotal, sPre) {
                    filtered = iMax - iTotal;
                    filtered_text = filtered ? " (" + filtered + " contact-pairs filtered out)" : "";
                    var cols = []
                    var table = $(selector + ' .dataTables_scrollBody');
                    cols_of_interest = [0, 1];
                    for (let [i, row] of [...table.find("tbody")[0].rows].entries()) {
                        for (let [j, cell] of [...row.cells].entries()) {
                            if (cols_of_interest.includes(j)) {
                                cols[j] = cols[j] || [];
                                cols[j].push(cell.innerText)
                            }
                        }
                    }
                    distinctPositions = [...new Set(cols[1].map((val, i) => val.split("-")).flat())]
                    //console.log(cols);
                    // distinctReceptors = [...new Set(cols[1])];
                    // distinctReceptorState = [...new Set(cols[1].map((val, i) => [cols[11]].reduce((a, arr) => [...a, arr[i]], [val])))];
                    // distinctReceptorState = [...new Set(distinctReceptorState.map(x => x[0] + "_" + x[1]))]
                    //console.log(iStart, iEnd, iMax, iTotal, sPre)
                    return "Showing " + iTotal + " contact-pairs covering "+distinctPositions.length+" positions"+filtered_text;
                  },
                "order": [],
                columnDefs: [{
                        type: "string",
                        targets: 1
                    },
                        {type: 'non-empty-string', targets: "_all"},
                    {
                        "width": "40px",
                        "targets": list_narrow_cols
                    }
                ]
            });

            console.timeEnd("Render DataTable " + tab_number);

            if (analys_mode == "#two-crystal-groups") {

                repeated_from_to_1 = make_range_number_cols(8, 18);
                repeated_from_to_2 = make_range_number_cols(28, 10);

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
                        // {
                        //     column_number: 5,
                        //     filter_type: "multi_select",
                        //     select_type: 'select2',
                        //     filter_default_label: "Type",
                        //     text_data_delimiter: "|",
                        //     filter_reset_button_text: false,
                        // }
                    ].concat(repeated_from_to_1).concat([{
                        column_number: 26,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "AA",
                        filter_reset_button_text: false,
                    }, {
                        column_number: 27,
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
            } else if (analys_mode == "#single-crystal-group") {
                repeated_from_to_1 = make_range_number_cols(6, 18);
                repeated_from_to_2 = make_range_number_cols(26, 7);
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
                        // {
                        //     column_number: 3,
                        //     filter_type: "multi_select",
                        //     select_type: 'select2',
                        //     filter_default_label: "Type",
                        //     text_data_delimiter: "|",
                        //     filter_reset_button_text: false,
                        // },
                        {
                            column_number: 24,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '60px'
                            },
                            filter_default_label: "AA",
                            filter_reset_button_text: false,
                        }, {
                            column_number: 25,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '60px'
                            },
                            filter_default_label: "AA",
                            filter_reset_button_text: false,
                        }
                    ].concat(repeated_from_to_1).concat(repeated_from_to_2), {
                        cumulative_filtering: false
                    }

                );
            } else if (analys_mode == "#single-crystal") {
                repeated_from_to_1 = make_range_number_cols(7, 11);
                repeated_from_to_2 = make_range_number_cols(20, 1);

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
                            column_number: 6,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '60px'
                            },
                            filter_default_label: "Type",
                            text_data_delimiter: "|",
                            filter_reset_button_text: false,
                        },
                        {
                            column_number: 18,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '60px'
                            },
                            filter_default_label: "AA",
                            filter_reset_button_text: false,
                        }, {
                            column_number: 19,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '60px'
                            },
                            filter_default_label: "AA",
                            filter_reset_button_text: false,
                        }
                    ].concat(repeated_from_to_1).concat(repeated_from_to_2), {
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
                "bInfo": true,
                "fnInfoCallback": function (oSettings, iStart, iEnd, iMax, iTotal, sPre) {
                    filtered = iMax - iTotal;
                    filtered_text = filtered ? " (" + filtered + " contact-pairs filtered out)" : "";
                    var cols = []
                    var table = $(selector + ' .dataTables_scrollBody');
                    cols_of_interest = [0, 1];
                    for (let [i, row] of [...table.find("tbody")[0].rows].entries()) {
                        for (let [j, cell] of [...row.cells].entries()) {
                            if (cols_of_interest.includes(j)) {
                                cols[j] = cols[j] || [];
                                cols[j].push(cell.innerText)
                            }
                        }
                    }
                    distinctPositions = [...new Set(cols[1].map((val, i) => val.split("-")).flat())]
                    //console.log(cols);
                    // distinctReceptors = [...new Set(cols[1])];
                    // distinctReceptorState = [...new Set(cols[1].map((val, i) => [cols[11]].reduce((a, arr) => [...a, arr[i]], [val])))];
                    // distinctReceptorState = [...new Set(distinctReceptorState.map(x => x[0] + "_" + x[1]))]
                    //console.log(iStart, iEnd, iMax, iTotal, sPre)
                    return "Showing " + iTotal + " contact-pairs covering "+distinctPositions.length+" positions"+filtered_text;
                  },
                paging: true,
                pageLength: 200,
                "order": [],
                columnDefs: [{
                        type: "string",
                        targets: 1
                    },
                        {type: 'non-empty-string', targets: "_all"},
                    {
                        "width": "40px",
                        "targets": list_narrow_cols
                    }
                ]
            });


            if (analys_mode == "#two-crystal-groups") {


                repeated_from_to_1 = make_range_number_cols(2, 6);
                repeated_from_to_2 = make_range_number_cols(10, 12);
                repeated_from_to_3 = make_range_number_cols(27, 13);
                repeated_from_to_4 = make_range_number_cols(42, 6);

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
                    }]).concat(repeated_from_to_2).concat(repeated_from_to_3).concat([{
                        column_number: 36,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "AA",
                        filter_reset_button_text: false,
                    }, {
                        column_number: 37,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "AA",
                        filter_reset_button_text: false,
                    }]).concat(repeated_from_to_4), {
                        cumulative_filtering: false
                    }

                );
            } else if (analys_mode == "#single-crystal-group") {

                repeated_from_to_1 = make_range_number_cols(2, 2);
                repeated_from_to_2 = make_range_number_cols(6, 6);
                repeated_from_to_3 = make_range_number_cols(17, 13);
                repeated_from_to_4 = make_range_number_cols(32, 6);

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
                    },{
                        column_number: 30,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "AA",
                        filter_reset_button_text: false,
                    }, {
                        column_number: 31,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "AA",
                        filter_reset_button_text: false,
                    }]).concat(repeated_from_to_2).concat([
                    //     {
                    //     column_number: 12,
                    //     filter_type: "multi_select",
                    //     select_type: 'select2',
                    //     select_type_options: {
                    //         width: '60px'
                    //     },
                    //     filter_default_label: "Type",
                    //     text_data_delimiter: "|",
                    //     filter_reset_button_text: false,
                    // }
                    ]).concat(repeated_from_to_3).concat(repeated_from_to_4), {
                            cumulative_filtering: false
                     }

                );
            } else if (analys_mode == "#single-crystal") {

                repeated_from_to_1 = make_range_number_cols(4, 3);
                repeated_from_to_2 = make_range_number_cols(8, 11);
                repeated_from_to_3 = make_range_number_cols(21, 4);

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
                        },
                        {
                            column_number: 19,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '60px'
                            },
                            filter_default_label: "AA",
                            filter_reset_button_text: false,
                        }, {
                            column_number: 20,
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
                    }]).concat(repeated_from_to_2).concat(repeated_from_to_3), {
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
                "bInfo": true,
                "fnInfoCallback": function (oSettings, iStart, iEnd, iMax, iTotal, sPre) {
                    filtered = iMax - iTotal;
                    filtered_text = filtered ? " (" + filtered + " positions filtered out)" : "";
                    var cols = []
                    var table = $(selector + ' .dataTables_scrollBody');
                    cols_of_interest = [0, 1];
                    for (let [i, row] of [...table.find("tbody")[0].rows].entries()) {
                        for (let [j, cell] of [...row.cells].entries()) {
                            if (cols_of_interest.includes(j)) {
                                cols[j] = cols[j] || [];
                                cols[j].push(cell.innerText)
                            }
                        }
                    }
                    distinctPositions = [...new Set(cols[1])]
                    //console.log(cols);
                    // distinctReceptors = [...new Set(cols[1])];
                    // distinctReceptorState = [...new Set(cols[1].map((val, i) => [cols[11]].reduce((a, arr) => [...a, arr[i]], [val])))];
                    // distinctReceptorState = [...new Set(distinctReceptorState.map(x => x[0] + "_" + x[1]))]
                    //console.log(iStart, iEnd, iMax, iTotal, sPre)
                    return "Showing " + iTotal + " positions"+filtered_text;
                  },
                paging: true,
                pageLength: 200,
                "order": [],
                columnDefs: [{
                        type: "string",
                        targets: 1
                    },
                        {type: 'non-empty-string', targets: "_all"},
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
                                width: '40px'
                            },
                            filter_default_label: "Seg",
                            filter_reset_button_text: false,
                        },
                        {
                            column_number: 1,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '40px'
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
                "bInfo": true,
                "fnInfoCallback": function (oSettings, iStart, iEnd, iMax, iTotal, sPre) {
                    filtered = iMax - iTotal;
                    filtered_text = filtered ? " (" + filtered + " positions filtered out)" : "";
                    var cols = []
                    var table = $(selector + ' .dataTables_scrollBody');
                    cols_of_interest = [0, 1];
                    for (let [i, row] of [...table.find("tbody")[0].rows].entries()) {
                        for (let [j, cell] of [...row.cells].entries()) {
                            if (cols_of_interest.includes(j)) {
                                cols[j] = cols[j] || [];
                                cols[j].push(cell.innerText)
                            }
                        }
                    }
                    distinctPositions = [...new Set(cols[1])]
                    //console.log(cols);
                    // distinctReceptors = [...new Set(cols[1])];
                    // distinctReceptorState = [...new Set(cols[1].map((val, i) => [cols[11]].reduce((a, arr) => [...a, arr[i]], [val])))];
                    // distinctReceptorState = [...new Set(distinctReceptorState.map(x => x[0] + "_" + x[1]))]
                    //console.log(iStart, iEnd, iMax, iTotal, sPre)
                    return "Showing " + iTotal + " positions"+filtered_text;
                  },
                paging: true,
                pageLength: 200,
                "order": [],
                columnDefs: [{
                        type: "string",
                        targets: 1
                    },
                        {type: 'non-empty-string', targets: "_all"},
                    {
                        "width": "40px",
                        "targets": list_narrow_cols
                    }
                ]
            });

            if (analys_mode == "#two-crystal-groups") {
                repeated_from_to_1 = make_range_number_cols(4, 24);
                repeated_from_to_2 = make_range_number_cols(28, 3);
                repeated_from_to_3 = make_range_number_cols(31, 2);
                repeated_from_to_4 = make_range_number_cols(35, 1);

                yadcf.init(btable,
                    [{
                            column_number: 0,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '40px'
                            },
                            filter_default_label: "Seg",
                            filter_reset_button_text: false,
                        },
                        {
                            column_number: 1,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '40px'
                            },
                            filter_default_label: "Pos",
                            filter_reset_button_text: false,
                        }
                    ].concat([{
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
                    }]).concat(repeated_from_to_1).concat(repeated_from_to_2).concat([{
                        column_number: 29,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "AA",
                        filter_reset_button_text: false,
                    }, {
                        column_number: 30,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "AA",
                        filter_reset_button_text: false,
                    },
                    {
                        column_number: 34,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '60px'
                        },
                        filter_default_label: "AA",
                        filter_reset_button_text: false,
                    }]).concat(repeated_from_to_3).concat(repeated_from_to_4), {
                        cumulative_filtering: false
                    }

                );
            } else if (analys_mode == "#single-crystal-group") {
                repeated_from_to_1 = make_range_number_cols(3, 7);
                repeated_from_to_2 = make_range_number_cols(11, 1);
                repeated_from_to_3 = make_range_number_cols(13, 1);

                yadcf.init(btable,
                    [{
                            column_number: 0,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '40px'
                            },
                            filter_default_label: "Seg",
                            filter_reset_button_text: false,
                        },
                        {
                            column_number: 1,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '40px'
                            },
                            filter_default_label: "Pos",
                            filter_reset_button_text: false,
                        },
                        {
                            column_number: 2,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '60px'
                            },
                            filter_default_label: "SS",
                            filter_reset_button_text: false,
                        }
                    ].concat(repeated_from_to_1).concat([{
                        column_number: 10,
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
                        filter_default_label: "AA",
                        filter_reset_button_text: false,
                    }]).concat(repeated_from_to_3), {
                        cumulative_filtering: false
                    }

                );
            } else if (analys_mode == "#single-crystal") {
                repeated_from_to_1 = make_range_number_cols(3, 6);
                repeated_from_to_2 = make_range_number_cols(11, 1);

                yadcf.init(btable,
                    [{
                            column_number: 0,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '40px'
                            },
                            filter_default_label: "Seg",
                            filter_reset_button_text: false,
                        },
                        {
                            column_number: 1,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '40px'
                            },
                            filter_default_label: "Pos",
                            filter_reset_button_text: false,
                        },
                        {
                            column_number: 2,
                            filter_type: "multi_select",
                            select_type: 'select2',
                            select_type_options: {
                                width: '60px'
                            },
                            filter_default_label: "SS",
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
                "bInfo": true,
                "fnInfoCallback": function (oSettings, iStart, iEnd, iMax, iTotal, sPre) {
                    filtered = iMax - iTotal;
                    filtered_text = filtered ? " (" + filtered + " positions filtered out)" : "";
                    var cols = []
                    var table = $(selector + ' .dataTables_scrollBody');
                    cols_of_interest = [0, 1];
                    for (let [i, row] of [...table.find("tbody")[0].rows].entries()) {
                        for (let [j, cell] of [...row.cells].entries()) {
                            if (cols_of_interest.includes(j)) {
                                cols[j] = cols[j] || [];
                                cols[j].push(cell.innerText)
                            }
                        }
                    }
                    distinctPositions = [...new Set(cols[1])]
                    //console.log(cols);
                    // distinctReceptors = [...new Set(cols[1])];
                    // distinctReceptorState = [...new Set(cols[1].map((val, i) => [cols[11]].reduce((a, arr) => [...a, arr[i]], [val])))];
                    // distinctReceptorState = [...new Set(distinctReceptorState.map(x => x[0] + "_" + x[1]))]
                    //console.log(iStart, iEnd, iMax, iTotal, sPre)
                    return "Showing " + iTotal + " positions"+filtered_text;
                  },
                paging: true,
                pageLength: 200,
                "order": [],
                columnDefs: [{
                        type: "string",
                        targets: 1
                    },
                        {type: 'non-empty-string', targets: "_all"},
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
                            width: '40px'
                        },
                        filter_default_label: "Seg",
                        filter_reset_button_text: false,
                    },
                    {
                        column_number: 1,
                        filter_type: "multi_select",
                        select_type: 'select2',
                        select_type_options: {
                            width: '40px'
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
    btable.columns.adjust().draw();

    $(".main_loading_overlay").hide();
    $('div.dataTables_scrollBody:visible').height('50vh');

    // Make sure browser-tables are not too wide.
    browser_table_div_width = $('.contact-browser:visible').width();
    if (browser_table_div_width > 2060) {
        browser_table_width = 2030;
    } else {
        browser_table_width = browser_table_div_width - 30;
    }
    $('.contact-browser .dataTables_wrapper').width(browser_table_width);
    console.timeEnd("renderDataTablesYadcf");
}

const types_to_short = { 'ionic': 'Ion', 'aromatic': 'Aro', 'polar': 'Pol', 'hydrophobic': 'Hyd', 'van-der-waals': 'vdW' }

var plot_options = {'tab1' : {}, 'tab2' : {}, 'tab3' : {}, 'tab4' : {}, 'tab5' : {}}
// First array contains number of columns per property that will be visualized
// 1,1,1 indicates three columns with individual coloring Options
// 2,1 indicates three columns of which the first two are combined and the last is individual
// type options: residuepair or residue from datatable or original data (for more options)
// TAB1 plot options - double
plot_options['tab1']['double'] = {}
plot_options['tab1']['double']['frequency'] = [[1,1,1], ['residuepair_datatable','residuepair_datatable','residuepair_datatable']]
plot_options['tab1']['double']['distance_diff'] = [[1], ['residuepair_original']]
plot_options['tab1']['double']['core_distance_diff'] = [[2], ['residue_original']]
plot_options['tab1']['double']['rotation_diff'] = [[2], ['residue_original']]
plot_options['tab1']['double']['rotamer_diff'] = [[2], ['residue_original']]
plot_options['tab1']['double']['SASA_diff'] = [[2], ['residue_original']]
plot_options['tab1']['double']['RSA_diff'] = [[2], ['residue_original']]
plot_options['tab1']['double']['presence_diff'] = [[2], ['residue_original']]
plot_options['tab1']['double']['consensus_SS'] = [[2], ['residue_datatable']]
plot_options['tab1']['double']['consensus_freq'] = [[2], ['residue_original']]
plot_options['tab1']['double']['no_gn'] = [[2], ['residue_datatable']]
plot_options['tab1']['double']['no_3d'] = [[2], ['residue_datatable']]
plot_options['tab1']['double']['class_conservation'] = [[1,1,1], ['residuepair_datatable','residuepair_datatable','residuepair_datatable']]

// TAB1 plot options - single
plot_options['tab1']['single'] = {}
plot_options['tab1']['single']['frequency'] = [[1], ['residuepair_datatable']]
plot_options['tab1']['single']['distance'] = [[1], ['residuepair_datatable']]
plot_options['tab1']['single']['core_distance'] = [[2], ['residue_datatable']]
plot_options['tab1']['single']['rotation'] = [[2], ['residue_datatable']]
plot_options['tab1']['single']['rotamer'] = [[2], ['residue_datatable']]
plot_options['tab1']['single']['SASA'] = [[2], ['residue_datatable']]
plot_options['tab1']['single']['RSA'] = [[2], ['residue_datatable']]
plot_options['tab1']['single']['presence'] = [[2], ['residue_datatable']]
plot_options['tab1']['single']['consensus_SS'] = [[2], ['residue_datatable']]
plot_options['tab1']['single']['consensus_freq'] = [[2], ['residue_datatable']]
plot_options['tab1']['single']['no_gn'] = [[2], ['residue_datatable']]
plot_options['tab1']['single']['no_3d'] = [[2], ['residue_datatable']]
plot_options['tab1']['single']['class_conservation'] = [[1], ['residuepair_datatable']]

// TAB1 plot options - single structure
plot_options['tab1']['structure'] = plot_options['tab1']['single']

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

        thead = '<tr> \
                          <th colspan="1" rowspan="2">Segment</th> \
                          <th colspan="1" rowspan="2">Positions</th> \
                          <th colspan="3" rowspan="2">Contact Frequency (%)</th> \
                          <th colspan="2" rowspan="2">Position no. contacts (in filtered rows)</th> \
                          <th colspan="1" rowspan="2">Net-<br>work no.</th> \
                          <th colspan="5" rowspan="2">Interaction types (%)</th> \
                          <th rowspan="2">Contact Ca distance (Å)</th> \
                          <th colspan="4">Backbone Ca movement</th> \
                          <th colspan="6">Sidechain differences</th> \
                          <th colspan="2" rowspan="2">Position presence %</th> \
                          <th colspan="2">Secondary structure</th> \
                          <th colspan="2"></th> \
                          <th colspan="4" rowspan="1">Absence in receptor or structure (%)</th> \
                          <th rowspan="2" colspan="4">Contact AA pair sequence conservation in class (%)</th> \
                        </tr> \
                        <tr> \
                          <th colspan="2">Distance to all other pos.</th> \
                          <th colspan="2">Angle to helix<br/>and 7TM axes</th> \
                          <th colspan="2">Rotamer</th> \
                          <th colspan="2">SASA</th> \
                          <th colspan="2">RSA</th> \
                          <th colspan="2">Consensus SS</th> \
                          <th colspan="2">Frequency %</th> \
                          <th colspan="2">No generic number (gap pos)</th> \
                          <th colspan="2">No 3D coordinates</th> \
                        </tr> \
                        <tr> \
                          <th class="dt-center"></th> \
                          <th class="dt-center">Pos1-Pos2</th> \
                          <th class="narrow_col">Set 1<br></th> \
                          <th class="narrow_col">Set 2<br></th> \
                          <th class="narrow_col">Diff<br></th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">No.</th> \
                          <th style="narrow_col">Ion</th> \
                          <th style="narrow_col">Pol</th> \
                          <th style="narrow_col">Aro</th> \
                          <th style="narrow_col">Hyd</th> \
                          <th style="narrow_col">vdW</th> \
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
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Set 1<br></th> \
                          <th class="narrow_col">Set 2<br></th> \
                          <th class="narrow_col">Diff<br></th> \
                          <th class="narrow_col">Max<br></th> \
                        </tr>';
        table.find('thead').html(thead);
        // two groups
        var proteins_1 = data['proteins1'].length
        var proteins_2 = data['proteins2'].length
        var pdbs_1 = data['pdbs1'].length
        var pdbs_2 = data['pdbs2'].length
        var pfs_1 = data['pfs1'].length
        var pfs_2 = data['pfs2'].length
        var normalized = data['normalized'];
        console.log('normalized!?',normalized);
        $.each(data['interactions'], function(i, v) {
            var gn1 = i.split(",")[0]
            var gn2 = i.split(",")[1]
            var sfreq1 = Math.round(100 * (normalized ? v['pf_freq_1'] : v['pdbs_freq_1']));
            var sfreq2 = Math.round(100* (normalized ? v['pf_freq_2'] : v['pdbs_freq_2']));
            var class_seq_cons = v['class_seq_cons'];
            if (normalized) {
                pos1_missing_1 = data['pfs1'].filter(x => data['missing'][gn1]['present'].includes(x)).length / data['pfs1'].length;
                pos1_missing_2 = data['pfs2'].filter(x => data['missing'][gn1]['present'].includes(x)).length / data['pfs2'].length;
                pos1_missing = Math.round(100*(pos1_missing_2-pos1_missing_1));
                pos2_missing_1 = data['pfs1'].filter(x => data['missing'][gn2]['present'].includes(x)).length / data['pfs1'].length;
                pos2_missing_2 = data['pfs2'].filter(x => data['missing'][gn2]['present'].includes(x)).length / data['pfs2'].length;
                pos2_missing = Math.round(100*(pos2_missing_2-pos2_missing_1));
            } else {
                pos1_missing_1 = data['pdbs1'].filter(x => data['missing'][gn1]['present'].includes(x)).length / data['pdbs1'].length;
                pos1_missing_2 = data['pdbs2'].filter(x => data['missing'][gn1]['present'].includes(x)).length / data['pdbs2'].length;
                pos1_missing = Math.round(100*(pos1_missing_2-pos1_missing_1));
                pos2_missing_1 = data['pdbs1'].filter(x => data['missing'][gn2]['present'].includes(x)).length / data['pdbs1'].length;
                pos2_missing_2 = data['pdbs2'].filter(x => data['missing'][gn2]['present'].includes(x)).length / data['pdbs2'].length;
                pos2_missing = Math.round(100*(pos2_missing_2-pos2_missing_1));

            }
            var diff_sfreq = sfreq1 - sfreq2;



            var class_seq_cons_diff = class_seq_cons[0] - class_seq_cons[1];
            var class_seq_cons_max = Math.max(class_seq_cons[0], class_seq_cons[1]);

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

            all_angles_1_set1 = data['all_angles_set1'][gn1];
            all_angles_1_set2 = data['all_angles_set2'][gn1];
            all_angles_2_set1 = data['all_angles_set1'][gn2];
            all_angles_2_set2 = data['all_angles_set2'][gn2];

            all_angles_1 = data['all_angles'][gn1];
            ss_pos1_set1 = [];
            ss_pos1_set2 = [];
            ss_pos2_set1 = [];
            ss_pos2_set2 = [];

            if (normalized) {
                pdbs = data['pfs1'].concat(data['pfs2']);
                set_1 = data['pfs1'];
                set_2 = data['pfs2'];
            } else {
                pdbs = data['pdbs1'].concat(data['pdbs2']);
                set_1 = data['pdbs1'];
                set_2 = data['pdbs2'];
            }

            types_count = {};
            Object.entries(v['types_count']).forEach(([key,value])=>{
                types_count_set1 = Math.round(100* (normalized ? value[0]['pf_freq'] : value[0]['pdb_freq'])); //set1
                types_count_set2 = Math.round(100* (normalized ? value[1]['pf_freq'] : value[1]['pdb_freq'])); //set2
                types_count[key] = [types_count_set1,types_count_set2,types_count_set1-types_count_set2];
            })

            // console.log(gn1, all_angles_1_set1, all_angles_1_set2)
            if (all_angles_1_set1) ss_pos1_set1 = Object.entries(all_angles_1_set1).filter(x => x[1].length > 6).map(x => x[1][12]);
            if (all_angles_1_set2) ss_pos1_set2 = Object.entries(all_angles_1_set2).filter(x => x[1].length > 6).map(x => x[1][12]);
            if (all_angles_2_set1) ss_pos2_set1 = Object.entries(all_angles_2_set1).filter(x => x[1].length > 6).map(x => x[1][12]);
            if (all_angles_2_set2) ss_pos2_set2 = Object.entries(all_angles_2_set2).filter(x => x[1].length > 6).map(x => x[1][12]);

            // pdbs.forEach(function(pdb){
            //     pdb_upper = pdb.toUpperCase();
            //     if (normalized) pdb_upper = pdb; //using pfs.. do not uppercase
            //     console.log(gn1,gn2,pdb_upper)
            //     if (all_angles_1_set1) {

            //         if (all_angles_1_set1.includes())

            //         let d1 = all_angles_1_set1[pdb_upper];
            //         if (d1.length) {
            //             if (set_1.includes(pdb)) {
            //                 ss_pos1_set1.push(d1[12]);
            //             } else if (set_2.includes(pdb)) {
            //                 ss_pos1_set2.push(d1[12]);
            //             }
            //         }
            //     }
            //     if (all_angles_2) {
            //         let d2 = all_angles_2[pdb_upper];
            //         if (d2.length) {
            //             if (set_1.includes(pdb)) {
            //                 ss_pos2_set1.push(d2[12])
            //             } else if (set_2.includes(pdb)) {
            //                 ss_pos2_set2.push(d2[12])
            //             }
            //         }
            //     }
            // });
            // console.log(gn1,gn2,ss_pos1_set1,ss_pos1_set2,ss_pos2_set1,ss_pos2_set2)
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
            // console.table(gn1,gn2,dssp);
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

            // avg distance ''
            distance_all_gn1 = '';
            if (gn1 in data['distances']) {
                distance_all_gn1 = data['distances'][gn1]['avg'];
            }
            distance_all_gn2 = '';
            if (gn2 in data['distances']) {
                distance_all_gn2 = data['distances'][gn2]['avg'];
            }

            tr = `
                    <tr class="clickable-row filter_rows" id="${i}">
                      <td class="dt-center">${seg1}-${seg2}</td>
                      <td class="dt-center">${gn1}-${gn2}</td>
                      <td class="narrow_col">${sfreq1}</td>
                      <td class="narrow_col">${sfreq2}</td>
                      <td class="narrow_col">${diff_sfreq}</td>
                      <td class="narrow_col pos1 count"></td>
                      <td class="narrow_col pos2 count"></td>
                      <td class="narrow_col group"></td>
                      <td class="dt-center angles_tooltip" data-set1="${types_count['ionic'][0]}" data-set2="${types_count['ionic'][1]}">${types_count['ionic'][2]}</td>
                      <td class="dt-center angles_tooltip" data-set1="${types_count['polar'][0]}" data-set2="${types_count['polar'][1]}">${types_count['polar'][2]}</td>
                      <td class="dt-center angles_tooltip" data-set1="${types_count['aromatic'][0]}" data-set2="${types_count['aromatic'][1]}">${types_count['aromatic'][2]}</td>
                      <td class="dt-center angles_tooltip" data-set1="${types_count['hydrophobic'][0]}" data-set2="${types_count['hydrophobic'][1]}">${types_count['hydrophobic'][2]}</td>
                      <td class="dt-center angles_tooltip" data-set1="${types_count['van-der-waals'][0]}" data-set2="${types_count['van-der-waals'][1]}">${types_count['van-der-waals'][2]}</td>
                      <td class="narrow_col">${distance}</td>
                      <td class="narrow_col" data-type="distance_all_avg">${distance_all_gn1}</td>
                      <td class="narrow_col" data-type="distance_all_avg">${distance_all_gn2}</td>
                      <td class="narrow_col angles_modal angles_tooltip" data-type="a_angle" data-pos="0" data-set1="${angles_1[1][1]}" data-set2="${angles_1[1][2]}">${angles_1[1][0]}</td>
                      <td class="narrow_col angles_modal angles_tooltip" data-type="a_angle" data-pos="1" data-set1="${angles_2[1][1]}" data-set2="${angles_2[1][2]}">${angles_2[1][0]}</td>
                      <td class="narrow_col angles_modal angles_tooltip" data-type="outer_angle" data-pos="0" data-set1="${angles_1[2][1]}" data-set2="${angles_1[2][2]}">${angles_1[2][0]}</td>
                      <td class="narrow_col angles_modal angles_tooltip" data-type="outer_angle" data-pos="1" data-set1="${angles_2[2][1]}" data-set2="${angles_2[2][2]}">${angles_2[2][0]}</td>
                      <td class="narrow_col angles_modal angles_tooltip" data-type="sasa" data-pos="0" data-set1="${angles_1[6][1]}" data-set2="${angles_1[6][2]}">${angles_1[6][0]}</td>
                      <td class="narrow_col angles_modal angles_tooltip" data-type="sasa" data-pos="1" data-set1="${angles_2[6][1]}" data-set2="${angles_2[6][2]}">${angles_2[6][0]}</td>
                      <td class="narrow_col angles_modal angles_tooltip" data-type="rsa" data-pos="0" data-set1="${angles_1[7][1]}" data-set2="${angles_1[7][2]}">${angles_1[7][0]}</td>
                      <td class="narrow_col angles_modal angles_tooltip" data-type="rsa" data-pos="1" data-set1="${angles_2[7][1]}" data-set2="${angles_2[7][2]}">${angles_2[7][0]}</td>
                      <td class="narrow_col">${pos1_presence}</td>
                      <td class="narrow_col">${pos2_presence}</td>
                      <td class="narrow_col">${dssp_pos1}</td>
                      <td class="narrow_col">${dssp_pos2}</td>
                      <td class="narrow_col">${dssp_pos1_freq}</td>
                      <td class="narrow_col">${dssp_pos2_freq}</td>
                      <td class="narrow_col">${pos1_missing}</td>
                      <td class="narrow_col">${pos2_missing}</td>
                      <td class="narrow_col">-</td>
                      <td class="narrow_col">-</td>
                      <td class="narrow_col">${class_seq_cons[0]}</td>
                      <td class="narrow_col">${class_seq_cons[1]}</td>
                      <td class="narrow_col">${class_seq_cons_diff}</td>
                      <td class="narrow_col" data_comparison="${class_seq_cons_diff}">${class_seq_cons_max}</td>
                    </tr>`;
            tbody.append(tr);
        });
    } else if ((data['proteins'].length > 1 && normalized) || (data['pdbs'].length > 1 && !normalized)) {
        thead = '<tr> \
                          <th colspan="1" rowspan="2">Segment</th> \
                          <th colspan="1" rowspan="2">Positions</th> \
                          <th colspan="1" rowspan="2">Contact Frequency (%)</th> \
                          <th colspan="2" rowspan="2">Position no. contacts (in filtered rows)</th> \
                          <th colspan="1" rowspan="2">Net-<br>work no.</th> \
                          <th rowspan="2" colspan="5">Interaction types (%)</th> \
                          <th rowspan="2">Contact Ca distance (Å)</th> \
                          <th colspan="4">Backbone Ca movement</th> \
                          <th colspan="6">Sidechain differences</th> \
                          <th colspan="2" rowspan="2">Position presence %</th> \
                          <th colspan="2">Secondary structure</th> \
                          <th colspan="2"></th> \
                          <th colspan="4" rowspan="1">Absence in receptor or structure (%)</th> \
                          <th rowspan="2">Class Seq Cons(%)</th> \
                        </tr> \
                        <tr> \
                          <th colspan="2">Distance to all other pos.</th> \
                          <th colspan="2">Angle to helix<br/>and 7TM axes</th> \
                          <th colspan="2">Rotamer</th> \
                          <th colspan="2">SASA</th> \
                          <th colspan="2">RSA</th> \
                          <th colspan="2">Consensus SS</th> \
                          <th colspan="2">Frequency %</th> \
                          <th colspan="2">No generic number (gap pos)</th> \
                          <th colspan="2">No 3D coordinates</th> \
                        </tr> \
                        <tr> \
                          <th class="dt-center"></th> \
                          <th class="dt-center">Pos1-Pos2</th> \
                          <th class="narrow_col">Set<br></th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">No.</th> \
                          <th style="narrow_col">Ion</th> \
                          <th style="narrow_col">Pol</th> \
                          <th style="narrow_col">Aro</th> \
                          <th style="narrow_col">Hyd</th> \
                          <th style="narrow_col">vdW</th> \
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
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">AA pairs</th> \
                        </tr>';
        table.find('thead').html(thead);
        var proteins = data['proteins'].length
        var pdbs_counts = data['pdbs'].length
        var normalized = data['normalized'];
        var pfs = data['pfs'].length
        $.each(data['interactions'], function(i, v) {
            var gn1 = i.split(",")[0]
            var gn2 = i.split(",")[1]
            var sfreq1 = Math.round(100* (normalized ? v['pf_freq'] : v['pdbs_freq']));

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


            types_count = {};
            Object.entries(v['types_count']).forEach(([key,value])=>{
                types_count[key] = Math.round(100* (normalized ? value['pf_freq'] : value['pdb_freq']));
            })

            all_angles_1 = data['all_angles'][gn1];
            all_angles_2 = data['all_angles'][gn2];
            ss_pos1_set1 = [];
            ss_pos2_set1 = [];
            if (normalized) {
                pdbs = data['pfs'];
            } else {
                pdbs = data['pdbs'];
            }
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
                      <td class="narrow_col" title="${v['pdbs']}">${sfreq1}</td>
                      <td class="narrow_col pos1 count"></td>
                      <td class="narrow_col pos2 count"></td>
                      <td class="narrow_col group"></td>
                      <td class="dt-center">${types_count['ionic']}</td>
                      <td class="dt-center">${types_count['polar']}</td>
                      <td class="dt-center">${types_count['aromatic']}</td>
                      <td class="dt-center">${types_count['hydrophobic']}</td>
                      <td class="dt-center">${types_count['van-der-waals']}</td>
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
                      <td class="narrow_col">-</td>
                      <td class="narrow_col">-</td>
                      <td class="narrow_col">-</td>
                      <td class="narrow_col">-</td>
                      <td class="narrow_col">${class_seq_cons}</td>
                    </tr>`;
            tbody.append(tr);
        });
    } else {
        thead = '<tr> \
                          <th colspan="1" rowspan="2">Segment</th> \
                          <th colspan="1" rowspan="2">Positions</th> \
                          <th colspan="1" rowspan="2">Positions GN</th> \
                          <th colspan="2" rowspan="2">Position no. contacts (in filtered rows)</th> \
                          <th colspan="1" rowspan="2">Net-<br>work no.</th> \
                          <th rowspan="2">Interaction types (%)</th> \
                          <th rowspan="2">Contact Ca distance (Å)</th> \
                          <th colspan="4">Backbone Ca movement</th> \
                          <th colspan="6">Sidechain differences</th> \
                          <th colspan="2">Secondary structure</th> \
                          <th rowspan="2">Class Seq Cons(%)</th> \
                        </tr> \
                        <tr> \
                          <th colspan="2">Distance to all other pos.</th> \
                          <th colspan="2">Angle to helix<br/>and 7TM axes</th> \
                          <th colspan="2">Rotamer</th> \
                          <th colspan="2">SASA</th> \
                          <th colspan="2">RSA</th> \
                          <th colspan="2">Consensus SS</th> \
                        </tr> \
                        <tr> \
                          <th class="dt-center"></th> \
                          <th class="dt-center">Pos1-Pos2</th> \
                          <th class="narrow_col">Pos1-Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">No.</th> \
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

            if (all_angles_1 && data['pdbs'][0].toUpperCase() in all_angles_1) {
                all_angles_1 = all_angles_1[data['pdbs'][0].toUpperCase()][12];
            } else {
                all_angles_1 = '';
            }

            if (all_angles_2 && data['pdbs'][0].toUpperCase() in all_angles_2) {
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
                      <td class="narrow_col pos1 count"></td>
                      <td class="narrow_col pos2 count"></td>
                      <td class="narrow_col group"></td>
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

    $(".angles_tooltip").hover(function() {
                                    hover_text = "Set1: "+Math.round($(this).data('set1')*10)/10 + " Set2: "+Math.round($(this).data('set2')*10)/10;
                                    $(this).css('cursor','pointer').attr('title', hover_text);
                                }, function() {
                                    $(this).css('cursor','auto');
                                }
                              );

    table.on('click', '.angles_modal', function(event) {
        // $(this)

        // figure out which cell is selected
        cell_index = $(this).index();
        data_type = $(this).data("type");
        data_pos = $(this).data("pos");

        gn_pair = $(this).closest("tr").attr('id').split(",");
        gn1 = gn_pair[0];
        gn2 = gn_pair[1];

        all_angles_1 = two_sets_data['all_angles'][gn1];
        all_angles_2 = two_sets_data['all_angles'][gn2];

        // Commented out is if both pairs tobe shown at the same time. Probably not relevant.
        // $("#resModal").find(".modal-body").html("<div id='modal_plotly_1' style='height:100%;width:100%;display: inline-block;'></div><div id='modal_plotly_2' style='height:100%;width:100%;display: inline-block;'></div>");
        $("#resModal").find(".modal-body").html("<div id='modal_plotly_1' style='height:100%;width:100%;display: inline-block;'></div>");
        $("#resModal").modal();

        //Slight wait, to be sure modal is open.

        if (data_pos == 0) {
            // odd cell number is pos1
            if (typeof all_angles_1 !== 'undefined')
              setTimeout(function(){ createBoxPlotResidue(gn1,'modal_plotly_1','angles',data_type) }, 500);

        } else {
            if (typeof all_angles_2 !== 'undefined')
              setTimeout(function(){ createBoxPlotResidue(gn2,'modal_plotly_1','angles',data_type) }, 500);
          }

    });

    console.timeEnd("RenderBrowser");
    gray_scale_table(table);
    //enable_hover(table);
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
        var pfs_1 = data['pfs1'].length
        var pfs_2 = data['pfs2'].length
        var normalized = data['normalized'];
        thead = '<tr> \
                          <th colspan="1" rowspan="2">Segment</th> \
                          <th colspan="1" rowspan="2">Positions</th> \
                          <th colspan="3" rowspan="2">Pos pair contact frequency (%)</th> \
                          <th colspan="3" rowspan="2">AA pair contact frequency (%)</th> \
                          <th colspan="2" rowspan="2">Amino acids</th> \
                          <th colspan="9" rowspan="1">AA occurrence in structure sets (%)</th> \
                          <th colspan="3" rowspan="2">Sequence conservation in class (%)</th> \
                          <th rowspan="2" colspan="5">Interaction types (%)</th> \
                          <th rowspan="2">Distance (Ca, Å)</th> \
                          <th colspan="4">Backbone Ca movement</th> \
                          <th colspan="6">Sidechain differences</th> \
                          <th colspan="2" rowspan="2">Position presence %</th> \
                          <th colspan="4">Secondary structure</th> \
                          <th colspan="4" rowspan="1">Absence in receptor or structure (%)</th> \
                        </tr> \
                        <tr> \
                          <th colspan="3">AA1</th> \
                          <th colspan="3">AA2</th> \
                          <th colspan="3">AA Pair</th> \
                          <th colspan="2">Distance to<br/>7TM axis (Å)</th> \
                          <th colspan="2">Angle to helix<br/>and 7TM axes</th> \
                          <th colspan="2">Rotamer</th> \
                          <th colspan="2">SASA</th> \
                          <th colspan="2">RSA</th> \
                          <th colspan="2">Consensus SS</th> \
                          <th colspan="2">Frequency %</th> \
                          <th colspan="2">No generic number (gap pos)</th> \
                          <th colspan="2">No 3D coordinates</th> \
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
                          <th style="narrow_col">Ion</th> \
                          <th style="narrow_col">Pol</th> \
                          <th style="narrow_col">Aro</th> \
                          <th style="narrow_col">Hyd</th> \
                          <th style="narrow_col">vdW</th> \
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
        var pfs_1 = data['pfs1'].length
        var pfs_2 = data['pfs2'].length
        var normalized = data['normalized'];
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
            if (normalized) {
                var aafreq1 = v2['set1']['interaction_freq_pf'];
                var aafreq2 = v2['set2']['interaction_freq_pf'];
            } else {
                var aafreq1 = v2['set1']['interaction_freq'];
                var aafreq2 = v2['set2']['interaction_freq'];
            }
            var diff_aafreq = (aafreq1 - aafreq2).toFixed(0);
            var aa1 = v2['aa1'];
            var aa2 = v2['aa2'];

            if (normalized) {
                denominator1 =  pfs_1;
                denominator2 =  pfs_2;
            } else {
                denominator1 =  pdbs_1;
                denominator2 =  pdbs_2;
            }
            var set1_occurance_aa1 = Math.round(100 * v2['set1']['occurance']['aa1'].length / denominator1);
            var set1_occurance_aa2 = Math.round(100 * v2['set1']['occurance']['aa2'].length / denominator1);
            var set1_occurance_pair = Math.round(100 * v2['set1']['occurance']['pair'].length / denominator1);

            var set2_occurance_aa1 = Math.round(100 * v2['set2']['occurance']['aa1'].length / denominator2);
            var set2_occurance_aa2 = Math.round(100 * v2['set2']['occurance']['aa2'].length / denominator2);
            var set2_occurance_pair = Math.round(100 * v2['set2']['occurance']['pair'].length / denominator2);
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
            if (normalized) {
                pdbs = data['pfs1'].concat(data['pfs2']);
            } else {
                pdbs = data['pdbs1'].concat(data['pdbs2']);
            }
            pdbs.forEach(function(pdb){
                pdb_upper = pdb.toUpperCase();
                if (normalized) pdb_upper = pdb; //using pfs.. do not uppercase
                if (all_angles_1) {
                    let d1 = all_angles_1[pdb_upper];
                    // console.log(d1[12],d1);
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
            // console.log([ss_pos1_set1,ss_pos1_set2,ss_pos2_set1,ss_pos2_set2]);
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

                      <td>${v2['types_freq']['ionic'][2]}</td>
                      <td>${v2['types_freq']['polar'][2]}</td>
                      <td>${v2['types_freq']['aromatic'][2]}</td>
                      <td>${v2['types_freq']['hydrophobic'][2]}</td>
                      <td>${v2['types_freq']['van-der-waals'][2]}</td>
                      <td class="narrow_col">${distance_2}</td>
                      <td class="narrow_col angles_modal angles_tooltip" data-type="core_distance" data-pos="0" data-set1="${angles_1[0][1]}" data-set2="${angles_1[0][2]}">${angles_1[0][0]}</td>
                      <td class="narrow_col angles_modal angles_tooltip" data-type="core_distance" data-pos="1" data-set1="${angles_2[0][1]}" data-set2="${angles_2[0][2]}">${angles_2[0][0]}</td>
                      <td class="narrow_col angles_modal angles_tooltip" data-type="a_angle" data-pos="0" data-set1="${angles_1[1][1]}" data-set2="${angles_1[1][2]}">${angles_1[1][0]}</td>
                      <td class="narrow_col angles_modal angles_tooltip" data-type="a_angle" data-pos="1" data-set1="${angles_2[1][1]}" data-set2="${angles_2[1][2]}">${angles_2[1][0]}</td>
                      <td class="narrow_col angles_modal angles_tooltip" data-type="outer_angle" data-pos="0" data-set1="${angles_1[2][1]}" data-set2="${angles_1[2][2]}">${angles_1[2][0]}</td>
                      <td class="narrow_col angles_modal angles_tooltip" data-type="outer_angle" data-pos="1" data-set1="${angles_2[2][1]}" data-set2="${angles_2[2][2]}">${angles_2[2][0]}</td>
                      <td class="narrow_col angles_modal angles_tooltip" data-type="sasa" data-pos="0" data-set1="${angles_1[6][1]}" data-set2="${angles_1[6][2]}">${angles_1[6][0]}</td>
                      <td class="narrow_col angles_modal angles_tooltip" data-type="sasa" data-pos="1" data-set1="${angles_2[6][1]}" data-set2="${angles_2[6][2]}">${angles_2[6][0]}</td>
                      <td class="narrow_col angles_modal angles_tooltip" data-type="rsa" data-pos="0" data-set1="${angles_1[7][1]}" data-set2="${angles_1[7][2]}">${angles_1[7][0]}</td>
                      <td class="narrow_col angles_modal angles_tooltip" data-type="rsa" data-pos="1" data-set1="${angles_2[7][1]}" data-set2="${angles_2[7][2]}">${angles_2[7][0]}</td>
                      <td class="narrow_col">${pos1_presence}</td>
                      <td class="narrow_col">${pos2_presence}</td>
                      <td class="narrow_col">${dssp_pos1}</td>
                      <td class="narrow_col">${dssp_pos2}</td>
                      <td class="narrow_col">${dssp_pos1_freq}</td>
                      <td class="narrow_col">${dssp_pos2_freq}</td>

                      <td class="narrow_col"></td>
                      <td class="narrow_col"></td>
                      <td class="narrow_col"></td>
                      <td class="narrow_col"></td>
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
                          <th colspan="3" rowspan="2">Sequence conservation in class (%)</th> \
                          <th rowspan="2" colspan="5">Interaction types (%)</th> \
                          <th rowspan="2">Distance (Ca, Å)</th> \
                          <th colspan="4">Backbone Ca movement</th> \
                          <th colspan="6">Sidechain differences</th> \
                          <th colspan="2" rowspan="2">Position presence %</th> \
                          <th colspan="4">Secondary structure</th> \
                          <th colspan="4" rowspan="1">Absence in receptor or structure (%)</th> \
                        </tr> \
                        <tr> \
                          <th colspan="1">AA1</th> \
                          <th colspan="1">AA2</th> \
                          <th colspan="1">AA Pair</th> \
                          <th colspan="2">Distance to<br/>7TM axis (Å)</th> \
                          <th colspan="2">Angle to helix<br/>and 7TM axes</th> \
                          <th colspan="2">Rotamer</th> \
                          <th colspan="2">SASA</th> \
                          <th colspan="2">RSA</th> \
                          <th colspan="2">Consensus SS</th> \
                          <th colspan="2">Frequency %</th> \
                          <th colspan="2">No generic number (gap pos)</th> \
                          <th colspan="2">No 3D coordinates</th> \
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
                          <th style="narrow_col">Ion</th> \
                          <th style="narrow_col">Pol</th> \
                          <th style="narrow_col">Aro</th> \
                          <th style="narrow_col">Hyd</th> \
                          <th style="narrow_col">vdW</th> \
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
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                          <th class="narrow_col">Pos1</th> \
                          <th class="narrow_col">Pos2</th> \
                        </tr>';
        table.find('thead').html(thead);
        // two groups
        var proteins = data['proteins'].length
        var pdbs_count = data['pdbs'].length
        var pfs = data['pfs'].length
        var normalized = data['normalized'];
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

            if (normalized) {
                pdbs = data['pfs'];
            } else {
                pdbs = data['pdbs'];
            }

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

                      <td>${v2['types_freq']['ionic']}</td>
                      <td>${v2['types_freq']['polar']}</td>
                      <td>${v2['types_freq']['aromatic']}</td>
                      <td>${v2['types_freq']['hydrophobic']}</td>
                      <td>${v2['types_freq']['van-der-waals']}</td>
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

                      <td class="narrow_col"></td>
                      <td class="narrow_col"></td>
                      <td class="narrow_col"></td>
                      <td class="narrow_col"></td>
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
                          <th rowspan="2">Interaction types (%)</th> \
                          <th rowspan="2">Distance (Ca, Å)</th> \
                          <th colspan="4">Backbone Ca movement</th> \
                          <th colspan="6">Sidechain differences</th> \
                          <th colspan="2">Secondary structure</th> \
                          <th colspan="4" rowspan="1">Absence in receptor or structure (%)</th> \
                        </tr> \
                        <tr> \
                          <th colspan="2">Distance to<br/>7TM axis (Å)</th> \
                          <th colspan="2">Angle to helix<br/>and 7TM axes</th> \
                          <th colspan="2">Rotamer</th> \
                          <th colspan="2">SASA</th> \
                          <th colspan="2">RSA</th> \
                          <th colspan="2">Consensus SS</th> \
                          <th colspan="2">No generic number (gap pos)</th> \
                          <th colspan="2">No 3D coordinates</th> \
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

            if (all_angles_1 && data['pdbs'][0].toUpperCase() in all_angles_1) {
                all_angles_1 = all_angles_1[data['pdbs'][0].toUpperCase()][12];
            } else {
                all_angles_1 = '';
            }

            if (all_angles_2 && data['pdbs'][0].toUpperCase() in all_angles_2) {
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

                      <td class="narrow_col"></td>
                      <td class="narrow_col"></td>
                      <td class="narrow_col"></td>
                      <td class="narrow_col"></td>
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

        data_type = $(this).data("type");
        data_pos = $(this).data("pos");

        all_angles_1 = two_sets_data['all_angles'][gn1];
        all_angles_2 = two_sets_data['all_angles'][gn2];

        $("#resModal").find(".modal-body").html("<div id='modal_plotly_1' style='height:100%;width:50%;display: inline-block;'></div><div id='modal_plotly_2' style='height:100%;width:50%;display: inline-block;'></div>");
        $("#resModal").modal();

        //Slight wait, to be sure modal is open.
        // console.log(pdbs_aa1,pdbs_aa2);
        if (data_pos == 0) {
            // odd cell number is pos1
            if (typeof all_angles_1 !== 'undefined')
              setTimeout(function(){ createBoxPlotResidue(gn1,'modal_plotly_1','angles',data_type,pdbs_aa1,aa1) }, 500);

        } else {
            if (typeof all_angles_2 !== 'undefined')
              setTimeout(function(){ createBoxPlotResidue(gn2,'modal_plotly_1','angles',data_type,pdbs_aa2,aa2) }, 500);
          }
        // setTimeout(function(){ createBoxPlotResidue(all_angles_1,'modal_plotly_1','angles',pdbs_aa1,aa1) }, 500);
        // setTimeout(function(){ createBoxPlotResidue(all_angles_2,'modal_plotly_2','angles',pdbs_aa2,aa2) }, 500);

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
                          <th colspan="1" rowspan="2">Segment</th> \
                          <th colspan="1" rowspan="2">Position</th> \
                          <th colspan="3" rowspan="2">Avg no. contact pairs</th> \
                          <th colspan="3" rowspan="1">No distinct or shared contacts</th> \
                          <th colspan="1" rowspan="2">Avg freq difference of all set1-2 contacts</th> \
                          <th colspan="5" rowspan="1">Seq consensus</th> \
                          <th colspan="2" rowspan="1">Class seq consensus</th> \
                          <th colspan="3">Backbone Ca movement</th> \
                          <th colspan="2">Sidechain differences</th> \
                        </tr> \
                        <tr> \
                          <th colspan="2">Distinct</th> \
                          <th colspan="1">Shared</th> \
                          <th colspan="2">AA</th> \
                          <th colspan="3">Conservation (%)</th> \
                          <th colspan="1">AA</th> \
                          <th colspan="1">Cons (%)</th> \
                          <th colspan="1">Distance to<br/>7TM axis (Å)</th> \
                          <th colspan="1">Angle to helix<br/>and 7TM axes</th> \
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
      thead = '<tr> \
                        <th colspan="1" rowspan="2">Segment</th> \
                        <th colspan="1" rowspan="2">Position</th> \
                        <th colspan="3" rowspan="2">Avg no. contact pairs</th> \
                        <th colspan="3" rowspan="1">No distinct or shared contacts</th> \
                        <th colspan="1" rowspan="2">Avg freq difference of all set1-2 contacts</th> \
                        <th colspan="5" rowspan="1">Seq consensus</th> \
                        <th colspan="2" rowspan="1">Class seq consensus</th> \
                        <th colspan="3">Backbone Ca movement</th> \
                        <th colspan="2">Sidechain differences</th> \
                      </tr> \
                      <tr> \
                        <th colspan="2">Distinct</th> \
                        <th colspan="1">Shared</th> \
                        <th colspan="2">AA</th> \
                        <th colspan="3">Conservation (%)</th> \
                        <th colspan="1">AA</th> \
                        <th colspan="1">Cons (%)</th> \
                        <th colspan="1">Distance to<br/>7TM axis (Å)</th> \
                        <th colspan="1">Angle to helix<br/>and 7TM axes</th> \
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

    //enable_hover(table);
    console.timeEnd("RenderBrowser3");
}

// TAB4 plot options - double
plot_options['tab4']['double'] = {}
plot_options['tab4']['double']['consensus_SS'] =  [[1,1], ['residue_datatable', 'residue_datatable']]
plot_options['tab4']['double']['consensus_freq'] = [[1,1,1], ['residue_datatable','residue_datatable','residue_datatable']]
plot_options['tab4']['double']['no_gn'] =  [[1,1], ['residue_datatable', 'residue_datatable']]
plot_options['tab4']['double']['no_3d'] =  [[1,1], ['residue_datatable', 'residue_datatable']]
plot_options['tab4']['double']['phi'] = [[1,1,1], ['residue_datatable', 'residue_datatable', 'residue_datatable']]
plot_options['tab4']['double']['psi'] = [[1,1,1], ['residue_datatable', 'residue_datatable', 'residue_datatable']]
plot_options['tab4']['double']['tau_angle'] = [[1,1,1], ['residue_datatable', 'residue_datatable', 'residue_datatable']]
plot_options['tab4']['double']['tau'] = [[1,1,1], ['residue_datatable', 'residue_datatable', 'residue_datatable']]
plot_options['tab4']['double']['theta'] = [[1,1,1], ['residue_datatable', 'residue_datatable', 'residue_datatable']]
plot_options['tab4']['double']['conservation'] = [[1,1,1], ['residue_datatable', 'residue_datatable', 'residue_datatable']]
plot_options['tab4']['double']['class_conservation'] = [[1], ['residue_datatable']]

// TAB4 plot options - single
plot_options['tab4']['single'] = {}
plot_options['tab4']['single']['consensus_SS'] = [[1], ['residue_datatable']]
plot_options['tab4']['single']['consensus_freq'] = [[1], ['residue_datatable']]
plot_options['tab4']['single']['no_gn'] = [[1], ['residue_datatable']]
plot_options['tab4']['single']['no_3d'] = [[1], ['residue_datatable']]
plot_options['tab4']['single']['phi'] = [[1], ['residue_datatable']]
plot_options['tab4']['single']['psi'] = [[1], ['residue_datatable']]
plot_options['tab4']['single']['tau_angle'] = [[1], ['residue_datatable']]
plot_options['tab4']['single']['tau'] = [[1], ['residue_datatable']]
plot_options['tab4']['single']['theta'] = [[1], ['residue_datatable']]
plot_options['tab4']['single']['conservation'] = [[1], ['residue_datatable']]
plot_options['tab4']['single']['class_conservation'] = [[1], ['residue_datatable']]

// TAB4 plot options - single structure
plot_options['tab4']['structure'] = plot_options['tab4']['single']


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
        var normalized = data['normalized'];
        thead = '<tr> \
                          <th colspan="1" rowspan="2">Segment</th> \
                          <th colspan="1" rowspan="2">Pos</th> \
                          <th colspan="5" rowspan="1">Secondary structure</th> \
                          <th colspan="4" rowspan="1">Absence in receptor or structure (%)</th> \
                          <th colspan="9" rowspan="1">Residue angles and dihedrals</th> \
                          <th colspan="9" rowspan="1">Helix turn angles and dihedrals</th> \
                          <th colspan="5" rowspan="1">Seq consensus</th> \
                          <th colspan="2" rowspan="1">Class seq consensus</th> \
                        </tr> \
                        <tr> \
                          <th colspan="2">Consensus SS</th> \
                          <th colspan="3">Frequency (%)</th> \
                          <th colspan="2">No generic number (gap pos)</th> \
                          <th colspan="2">No 3D coordinates</th> \
                          <th colspan="3">Phi dihedral<br/><span class="small">(N(+1)-C-Ca-N)</span></th> \
                          <th colspan="3">Psi dihedral<br/><span class="small">(C-Ca-N-C(-1))</span></th> \
                          <th colspan="3">Tau angle<br/><span class="small">(N-Ca-C)</span></th> \
                          <th colspan="3">Tau dihedral<br/><span class="small">(Ca(+1)-Ca-Ca(-1)-Ca(-2))</span></th> \
                          <th colspan="3">Next tau dihedral<br/><span class="small">(Ca(+2)-Ca(+1)-Ca-Ca(-1))</span></th> \
                          <th colspan="3">Theta angle<br/><span class="small">(Ca(+1)-Ca-Ca(-1))</span></th> \
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
        $.each(data['tab4'], function(i, v) {

            var seg = data['segm_lookup'][i];
            if (seg == 'ECL1' || seg == 'ECL2') return true;

            var angles1 = v['angles_set1'];
            var angles2 = v['angles_set2'];
            var angles_diff = v['angles'];
            // index_names = {0:'core_distance',1:'a_angle',2:'outer_angle',3:'tau',4:'phi',5:'psi',6: 'sasa',7: 'rsa',8:'theta',9:'hse',10:'tau_angle'}
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
            if (normalized) {
                pdbs = data['pfs1'].concat(data['pfs2']);
                dssp_set1 = data['pfs1'];
                dssp_set2 = data['pfs2'];
            } else {
                pdbs = data['pdbs1'].concat(data['pdbs2']);
                dssp_set1 = data['pdbs1'];
                dssp_set2 = data['pdbs2'];
            }

            // missing_1 = [...new Set([...data['missing'][i]['present'], ...dssp_set1])].length / dssp_set1.length;
            // missing_2 = [...new Set([...data['missing'][i]['present'], ...dssp_set2])].length / dssp_set2.length;

            missing_1 = Math.round(100*dssp_set1.filter(x => data['missing'][i]['present'].includes(x)).length / dssp_set1.length);
            missing_2 = Math.round(100 * dssp_set2.filter(x => data['missing'][i]['present'].includes(x)).length / dssp_set2.length);


            all_angles_1_set1 = data['all_angles_set1'][i];
            all_angles_1_set2 = data['all_angles_set2'][i];
            if (all_angles_1_set1) ss_pos1_set1 = Object.entries(all_angles_1_set1).filter(x => x[1].length > 6).map(x => x[1][12]);
            if (all_angles_1_set2) ss_pos1_set2 = Object.entries(all_angles_1_set2).filter(x => x[1].length > 6).map(x => x[1][12]);


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

                      <td class="narrow_col">${missing_1}</td>
                      <td class="narrow_col">${missing_2}</td>
                      <td class="narrow_col"></td>
                      <td class="narrow_col"></td>

                      <td class="narrow_col">${angles1[4]}</td>
                      <td class="narrow_col">${angles2[4]}</td>
                      <td class="narrow_col">${angles_diff[4][0]}</td>

                      <td class="narrow_col">${angles1[5]}</td>
                      <td class="narrow_col">${angles2[5]}</td>
                      <td class="narrow_col">${angles_diff[5][0]}</td>

                      <td class="narrow_col">${angles1[10]}</td>
                      <td class="narrow_col">${angles2[10]}</td>
                      <td class="narrow_col">${Math.abs(Math.round(angles1[10]-angles2[10]))}</td>

                      <td class="narrow_col">${angles1[3]}</td>
                      <td class="narrow_col">${angles2[3]}</td>
                      <td class="narrow_col">${angles_diff[3][0]}</td>

                      <td class="narrow_col"></td>
                      <td class="narrow_col"></td>
                      <td class="narrow_col"></td>

                      <td class="narrow_col">${angles1[8]}</td>
                      <td class="narrow_col">${angles2[8]}</td>
                      <td class="narrow_col">${angles_diff[8][0]}</td>

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
        var proteins = data['proteins'].length;
        var pdbs_count = data['pdbs'].length;
        var normalized = data['normalized'];
        thead = '<tr> \
                          <th colspan="1" rowspan="2">Segment</th> \
                          <th colspan="1" rowspan="2">Pos</th> \
                          <th colspan="2" rowspan="1">Secondary structure</th> \
                          <th colspan="3" rowspan="1">Residue angles</th> \
                          <th colspan="3" rowspan="1">Helix turn angle</th> \
                          <th colspan="2" rowspan="1">Seq consensus</th> \
                          <th colspan="2" rowspan="1">Class seq consensus</th> \
                        </tr> \
                        <tr> \
                          <th colspan="1">Consensus SS</th> \
                          <th colspan="1">Frequency (%)</th> \
                          <th colspan="1">Phi dihedral<br/><span class="small">(N(+1)-C-Ca-N)</span></th> \
                          <th colspan="1">Psi dihedral<br/><span class="small">(C-Ca-N-C(-1))</span></th> \
                          <th colspan="1">Tau angle<br/><span class="small">(N-Ca-C)</span></th> \
                          <th colspan="1">Tau dihedral<br/><span class="small">(Ca(+1)-Ca-Ca(-1)-Ca(-2))</span></th> \
                          <th colspan="1">Next tau dihedral<br/><span class="small">(Ca(+2)-Ca(+1)-Ca-Ca(-1))</span></th> \
                          <th colspan="1">Theta angle<br/><span class="small">(Ca(+1)-Ca-Ca(-1))</span></th> \
                          <th colspan="1">AA</th> \
                          <th colspan="1">Conservation (%)</th> \
                          <th colspan="1">AA</th> \
                          <th colspan="1">Cons (%)</th> \
                        </tr> \
                        <tr> \
                          <th class="dt-center"></th> \
                          <th class="dt-center"></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                        </tr>';
        table.find('thead').html(thead);
        tr_list = ''
        $.each(data['tab4'], function(i, v) {

            var seg = data['segm_lookup'][i];
            if (seg == 'ECL1' || seg == 'ECL2') return true;

            var angles = v['angles_set'];
            // index_names = {0:'core_distance',1:'a_angle',2:'outer_angle',3:'tau',4:'phi',5:'psi',6: 'sasa',7: 'rsa',8:'theta',9:'hse',10:'tau_angle'}
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
            if (normalized) {
                pdbs = data['pfs'];
            } else {
                pdbs = data['pdbs'];
            }
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
                    freq = Math.round(100*most[1]/list.length);
                    most = most[0];
                } else {
                    freq = '';
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
                      <td class="narrow_col">${dssp[0][1]}</td>

                      <td class="narrow_col">${angles[4]}</td>

                      <td class="narrow_col">${angles[5]}</td>

                      <td class="narrow_col">${angles[10]}</td>

                      <td class="narrow_col">${angles[3]}</td>

                      <td class="narrow_col"></td>

                      <td class="narrow_col">${angles[8]}</td>

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
                          <th colspan="1" rowspan="2">Segment</th> \
                          <th colspan="1" rowspan="2">Pos</th> \
                          <th colspan="1" rowspan="1">Secondary structure</th> \
                          <th colspan="3" rowspan="1">Residue angles</th> \
                          <th colspan="3" rowspan="1">Helix turn angle</th> \
                          <th colspan="1" rowspan="1">Seq</th> \
                          <th colspan="2" rowspan="1">Class seq consensus</th> \
                        </tr> \
                        <tr> \
                          <th colspan="1">SS</th> \
                          <th colspan="1">Phi dihedral<br/><span class="small">(N(+1)-C-Ca-N)</span></th> \
                          <th colspan="1">Psi dihedral<br/><span class="small">(C-Ca-N-C(-1))</span></th> \
                          <th colspan="1">Tau angle<br/><span class="small">(N-Ca-C)</span></th> \
                          <th colspan="1">Tau dihedral<br/><span class="small">(Ca(+1)-Ca-Ca(-1)-Ca(-2))</span></th> \
                          <th colspan="1">Next tau dihedral<br/><span class="small">(Ca(+2)-Ca(+1)-Ca-Ca(-1))</span></th> \
                          <th colspan="1">Theta angle<br/><span class="small">(Ca(+1)-Ca-Ca(-1))</span></th> \
                          <th colspan="1">AA</th> \
                          <th colspan="1">AA</th> \
                          <th colspan="1">Cons (%)</th> \
                        </tr> \
                        <tr> \
                          <th class="dt-center"></th> \
                          <th class="dt-center"></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                          <th class="narrow_col"></th> \
                        </tr>';
        table.find('thead').html(thead);
        tr_list = ''
        $.each(data['tab4'], function(i, v) {

            var seg = data['segm_lookup'][i];
            if (seg == 'ECL1' || seg == 'ECL2') return true;

            var angles = v['angles_set'];
            // index_names = {0:'core_distance',1:'a_angle',2:'outer_angle',3:'tau',4:'phi',5:'psi',6: 'sasa',7: 'rsa',8:'theta',9:'hse',10:'tau_angle'}
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

            if (all_angles_1 && data['pdbs'][0].toUpperCase() in all_angles_1) {
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

                      <td class="narrow_col">${angles[10]}</td>

                      <td class="narrow_col">${angles[3]}</td>

                      <td class="narrow_col"></td>

                      <td class="narrow_col">${angles[8]}</td>

                      <td class="narrow_col">${set_seq_cons_aa}</td>

                      <td class="narrow_col">${class_cons_aa}</td>
                      <td class="narrow_col">${class_cons_freq}</td>

                    </tr>`;
        });
        // insert natively for speed increase on Chrome
        tbody[0].innerHTML = tr_list;

    }

    //enable_hover(table);
    console.timeEnd("RenderBrowser4");
}

// TAB5 plot options - double
plot_options['tab5']['double'] = {}
plot_options['tab5']['double']['core_distance_diff'] = [[1], ['residue_original']]
plot_options['tab5']['double']['rotation_diff'] = [[1], ['residue_original']]
plot_options['tab5']['double']['HSE_diff'] = [[1], ['residue_original']]

// TAB5 plot options - single
plot_options['tab5']['single'] = {}
plot_options['tab5']['single']['core_distance'] = [[1], ['residue_original']]
plot_options['tab5']['single']['rotation'] = [[1], ['residue_original']]
plot_options['tab5']['single']['HSE_diff'] = [[1], ['residue_original']]

// TAB5 plot options - single structure
plot_options['tab5']['structure'] = plot_options['tab5']['single']

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

    var thead = "";
    if (data['proteins2']) {
        thead += '<tr> \
                        <th colspan="1" rowspan="2">Seg-<br>ment</th> \
                        <th colspan="1" rowspan="2">Pos</th> \
                        <th colspan="3">Backbone Ca movement</th> \
                        <th colspan="1" rowspan="2">Ca half-sphere exposure (Å&sup2;)</th> \
                        <th colspan="3">Sidechain differences</th> \
                        <th colspan="5" rowspan="1">Seq consensus</th> \
                        <th colspan="2" rowspan="1">Class seq consensus</th> \
                        </tr> \
                        <tr> \
                        <th colspan="1">Avg distance to<br/>residues</th> \
                        <th colspan="1">Distance to<br/>7TM axis (Å)</th> \
                        <th colspan="1">Angle to helix<br/>and 7TM axes</th> \
                        <th colspan="1">Rotamer</th> \
                        <th colspan="1">SASA (Å&sup2;)</th> \
                        <th colspan="1">RSA (Å&sup2;)</th> \
                        <th colspan="2">AA</th> \
                        <th colspan="3">Conservation (%)</th> \
                        <th colspan="1">AA</th> \
                        <th colspan="1">Cons (%)</th> \
                        </tr> \
                        <tr> \
                        <th class="dt-center"></th> \
                        <th class="dt-center"></th> \
                        <th class="narrow_col"></th> \
                        <th class="narrow_col"></th> \
                        <th class="narrow_col"></th> \
                        <th class="narrow_col"></th> \
                        <th class="narrow_col"></th> \
                        <th class="narrow_col"></th> \
                        <th class="narrow_col"></th> \
                        <th class="narrow_col">Set 1<br></th> \
                        <th class="narrow_col">Set 2<br></th> \
                        <th class="narrow_col">Set 1<br></th> \
                        <th class="narrow_col">Set 2<br></th> \
                        <th class="narrow_col">Diff<br></th> \
                        <th class="narrow_col"></th> \
                        <th class="narrow_col"></th> \
                        </tr>';
    } else {
        thead += '<tr> \
                <th colspan="1" rowspan="2">Seg-<br>ment</th> \
                <th colspan="1" rowspan="2">Pos</th> \
                <th colspan="3">Backbone Ca movement</th> \
                <th colspan="1" rowspan="2">Ca half-sphere exposure (Å&sup2;)</th> \
                <th colspan="3">Sidechain differences</th> \
                <th colspan="2" rowspan="1">Seq consensus</th> \
                <th colspan="2" rowspan="1">Class seq consensus</th> \
                </tr> \
                <tr> \
                <th colspan="1">Avg distance to<br/>residues</th> \
                <th colspan="1">Distance to<br/>7TM axis (Å)</th> \
                <th colspan="1">Angle to helix<br/>and 7TM axes</th> \
                <th colspan="1">Rotamer</th> \
                <th colspan="1">SASA (Å&sup2;)</th> \
                <th colspan="1">RSA (Å&sup2;)</th> \
                <th colspan="1">AA</th> \
                <th colspan="1">Conservation (%)</th> \
                <th colspan="1">AA</th> \
                <th colspan="1">Cons (%)</th> \
                </tr> \
                <tr> \
                <th class="dt-center"></th> \
                <th class="dt-center"></th> \
                <th class="narrow_col"></th> \
                <th class="narrow_col"></th> \
                <th class="narrow_col"></th> \
                <th class="narrow_col"></th> \
                <th class="narrow_col"></th> \
                <th class="narrow_col"></th> \
                <th class="narrow_col"></th> \
                <th class="narrow_col"><br></th> \
                <th class="narrow_col"><br></th> \
                <th class="narrow_col"></th> \
                <th class="narrow_col"></th> \
                </tr>';
    }
    table.find('thead').html(thead);
    tr_list = ''
    if (data['proteins2']) {
        $.each(data['tab4'], function(i, v) {

            // console.log(i,v);
            var seg = v['ps'];
            var angles = v['angles'];
            // index_names = {0:'core_distance',1:'a_angle',2:'outer_angle',3:'tau',4:'phi',5:'psi',6: 'sasa',7: 'rsa',8:'theta',9:'hse',10:'tau_angle'}
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

            var pdbs_1 = data['pdbs1'].length
            var pdbs_2 = data['pdbs2'].length

            var set1_seq_cons_aa = v['set1_seq_cons'][0];
            var set2_seq_cons_aa = v['set2_seq_cons'][0];
            var set1_seq_cons_freq = Math.round(100 * v['set1_seq_cons'][1] / pdbs_1);
            var set2_seq_cons_freq = Math.round(100 * v['set2_seq_cons'][2] / pdbs_2);
            var diff_seq_cons_freq = Math.round((set1_seq_cons_freq - set2_seq_cons_freq));

            var class_cons_aa = v['class_cons'][0];
            var class_cons_freq = Math.round(100 * v['class_cons'][1])

            if (i in data['distances']) {
                distance = data['distances'][i]['avg'];
            } else {
                // console.log('no ', i, 'in distances');
                distance = '';
            }

            tr = ''
            tr_list += `
                    <tr class="clickable-row filter_rows" id="${i}">
                    <td class="dt-center">${seg}</td>
                    <td class="dt-center">${i}</td>
                    <td class="narrow_col">${distance}</td>
                    <td class="narrow_col">${angles[0][0]}</td>
                    <td class="narrow_col">${angles[1][0]}</td>
                    <td class="narrow_col">${angles[9][0]}</td>
                    <td class="narrow_col">${angles[2][0]}</td>
                    <td class="narrow_col">${angles[6][0]}</td>
                    <td class="narrow_col">${angles[7][0]}</td>

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
    } else {
        $.each(data['tab4'], function(i, v) {

            // console.log(i,v);
            var seg = v['ps'];
            var angles = v['angles'];
            // index_names = {0:'core_distance',1:'a_angle',2:'outer_angle',3:'tau',4:'phi',5:'psi',6: 'sasa',7: 'rsa',8:'theta',9:'hse',10:'tau_angle'}
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

            var pdbs_count = data['pdbs'].length

            var set_seq_cons_aa = v['set_seq_cons'][0];
            var set_seq_cons_freq = Math.round(100 * v['set_seq_cons'][1] / pdbs_count);

            var class_cons_aa = v['class_cons'][0];
            var class_cons_freq = Math.round(100 * v['class_cons'][1]);

            tr = ''
            tr_list += `
                    <tr class="clickable-row filter_rows" id="${i}">
                    <td class="dt-center">${seg}</td>
                    <td class="dt-center">${i}</td>
                    <td class="narrow_col">pair</td>
                    <td class="narrow_col">${angles[0]}</td>
                    <td class="narrow_col">${angles[1]}</td>
                    <td class="narrow_col">${angles[9]}</td>
                    <td class="narrow_col">${angles[2]}</td>
                    <td class="narrow_col">${angles[6]}</td>
                    <td class="narrow_col">${angles[7]}</td>

                    <td class="narrow_col">${set_seq_cons_aa}</td>
                    <td class="narrow_col">${set_seq_cons_freq}</td>

                    <td class="narrow_col">${class_cons_aa}</td>
                    <td class="narrow_col">${class_cons_freq}</td>

                    </tr>`;
            // tbody.append(tr);
        });
    }
    // insert natively for speed increase on Chrome
    tbody[0].innerHTML = tr_list;

    //enable_hover(table)
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
            if (cell.innerText.charAt(0) == '-' && cell.innerText.length > 1) {
                $(cell).addClass("minus");
            }
        }
    }
    var maxmin = [];
    cols.forEach(function(col, index) {
        var max = Math.max.apply(null, col);
        var min = Math.min.apply(null, col);
        var abs_max = Math.max.apply(null, [max, min].map(Math.abs));
        maxmin.push([max, min,abs_max]);
    });
    // console.time('Greyscale cells');

    // Get the header texts to find out which are "set specific"
    var h_cols = []
    for (let [i, row] of [...table.find("thead")[0].rows].entries()) {
        for (let [j, cell] of [...row.cells].entries()) {
            // h_cols[j] = h_cols[j] || [];
            h_cols[j] = cell.innerText;
        }
    }
    var cell_count = 0;
    for (let [i, row] of [...table.find("tbody")[0].rows].entries()) {
        for (let [j, cell] of [...row.cells].entries()) {
            c_maxmin = maxmin[j];
            c_header = h_cols[j];
            value = parseFloat(cell.innerText);
            comparison_toggle = cell.hasAttribute("data_comparison")
            if (!(isNaN(value) || isNaN(c_maxmin[0]) || isNaN(c_maxmin[1]))) {
                scale = Math.abs(value) / c_maxmin[2];
                var color = { r: 255, g: 255, b: 255 };
                if ((c_header.includes('Set 2') || value < 0 || (comparison_toggle && parseFloat(cell.getAttribute("data_comparison")) < 0)) && !(c_header.includes('Set 1'))) {
                    // if the header is a set two, then make it red
                    color = { r: 255, g: 255-(255-153)*scale, b: 255-(255-153)*scale }; //red
                } else if (value > 0) {
                    // Positive numbers are blue either cos they are set 1 or cos "set 1 has most"
                    // This is also used for single set/structure
                    color = { r: 255-(255-153)*scale, g: 255-(255-204)*scale, b: 255 }; //blue
                }
                var hex = rgb2hex(color.r, color.g, color.b);
                cell.setAttribute("bgcolor", hex);
                cell_count++;
            }
        }
    }
    // console.timeEnd('Greyscale cells');
    console.log(cell_count, 'cells greyscaled');
    console.timeEnd('Greyscale');
}

function make_abs_values(e,table) {
    $(".main_loading_overlay").show();
    console.time('Abs values')

    if ($(e).attr('changed') == '0') {
        $(e).html('Change back to original values');
        $(e).attr('changed','1')
    } else {
        $(e).html('Change negative to absolute values');
        $(e).attr('changed','0')

    }
    console.log(table);

    myVar = setTimeout(function () {
            var dt_table = $(table).DataTable();

            c = 0;
            dt_table.cells('.minus').every(function () {
                d = String(this.data());
                c += 1;
                if (d.charAt(0) == '-' && d.length > 1) {
                    this.data(d.substr(1));
                    $(this.node()).addClass("minus_removed");
                } else if ($(this.node()).hasClass("minus_removed")) {
                    this.data("-" + d);
                    $(this.node()).removeClass("minus_removed");
                }
            })

            dt_table.draw(false);
            console.timeEnd('Abs values')
            console.log(c + ' cells changed');
            $(".main_loading_overlay").hide();
        }
        , 100);
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

    // get residue_values
    var residue_values = getColumn(rows, 1)

    if (residue_positions[0].indexOf("-")>=1){
      residue_positions1 = residue_positions.map(function(e){return e.split("-")[0]});
      residue_positions2 = residue_positions.map(function(e){return e.split("-")[1]});
      residue_positions = residue_positions1.concat(residue_positions2)

      // copy values from the same array
      if (!Array.isArray(columnNumber) || columnNumber.length == 1)
        residue_values = residue_values.concat(residue_values)
    }

    // Filter NaNs
    if (type!="consensus_SS"){
      residue_values = residue_values.map(function(e){ return parseInt(e);})
      if (Array.isArray(columnNumber) && columnNumber.length > 1)
         residue_values = residue_values.concat(getColumn(rows, 2).map(function(e){ return parseInt(e);}));

      residue_positions = residue_positions.filter( function(value, index){ return !isNaN(residue_values[index]); })
      residue_values = residue_values.filter( function(value, index){ return !isNaN(residue_values[index]); })
    } else {
      // remove positions with no data
      residue_positions = residue_positions.filter( function(value, index){ return !(residue_values[index]==""); })
      residue_values = residue_values.filter( function(value, index){ return !(residue_values[index]==""); })
    }

    // Remove duplicates
    var tmp_rp = []
    var tmp_v = []
    for (var i = 0; i < residue_positions.length; i++) {
      if (tmp_rp.indexOf(residue_positions[i]) == -1){
          tmp_rp.push(residue_positions[i])
          tmp_v.push(residue_values[i])
      }
    }
    residue_positions = tmp_rp
    residue_values = tmp_v


    // Identify range
    var valMax = Math.max(...residue_values)
    var valMin = Math.min(...residue_values)
    var palette = "rwb"
    if (type in dataType){
      valMax = dataType[type][2]
      valMin = dataType[type][1]
      palette = dataType[type][3]
    } else {
      console.log("TYPE not found: " + type)
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
dataType["frequency"]      = [true, 0, 100, "wb"]
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


    var red = {red:255, green:0, blue: 0}
    var red = {red:195, green:74, blue: 54}
    var blue = {red:0, green:0, blue: 255}
    var blue = {red:0, green:140, blue: 204}
    var green = {red:0, green:255, blue: 0}
    var green = {red:0, green:201, blue: 167}
    var white = {red:255, green:255, blue: 255}
    var yellow = {red:255, green:255, blue: 0}
    var yellow = {red:255, green:255, blue: 0}
    var black = {red:0, green:0, blue: 0}

    switch(palette){
        case "rwb": // red-white-blue
          return colorGradient(value/max, red, white, blue)
          break;
        case "bwr": // blue-white-red
          return colorGradient(value/max, blue, white, red)
          break;
        case "ryg": // red-yellow-green
          return colorGradient(value/max, red, yellow, green)
          break;
        case "gyr": // green-yellow-red
          return colorGradient(value/max, green, yellow, red)
          break;
        case "rgb":
          return colorGradient(value/max, red, green, blue)
          break;
        case "wr": // white-red
          return colorGradient(value/max, white, red)
          break;
        case "wg": // white-green
          return colorGradient(value/max, white, green)
          break;
        case "wb": // white-blue
          return colorGradient(value/max, white, blue)
          break;
        case "wy": // white-yellow
            return colorGradient(value/max, white, yellow)
            break;
        case "wo": // white-orange
            return colorGradient(value/max, white, {red:255, green:150, blue: 113})
            break;
        case "rb": // red-blue
          return colorGradient(value/max, red, blue)
          break;
        case "wp": // white-purple
            return colorGradient(value / max, white, { red: 128, green: 0, blue: 128 })
            break;
        case "grey": // grey
            return colorGradient(value / max, white, black)
            break;
        // ADDON if you're missing gradient values
        case "br": // blue-red
        default:
          return colorGradient(value/max, blue, red)
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

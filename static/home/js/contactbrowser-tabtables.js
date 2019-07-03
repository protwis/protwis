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
    console.log('filtered positions!', filtered_gn_pairs);

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

            console.log("do tab", tab_number)
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
                // function myCustomFilterFunction(filterVal, columnVal) {
                //     var found;
                //     if (columnVal === '') {
                //         return true;
                //     }
                //     switch (filterVal) {
                //     case 'happy':
                //         found = columnVal.search(/:-\]|:\)|Happy|JOY|:D/g);
                //         break;
                //     case 'sad':
                //         found = columnVal.search(/:\(|Sad|:'\(/g);
                //         break;
                //     case 'angry':
                //         found = columnVal.search(/!!!|Arr\.\.\./g);
                //         break;
                //     case 'lucky':
                //         found = columnVal.search(/777|Bingo/g);
                //         break;
                //     case 'january':
                //         found = columnVal.search(/01|Jan/g);
                //         break;
                //     default:
                //         found = 1;
                //         break;
                //     }

                //     if (found !== -1) {
                //         return true;
                //     }
                //     return false;
                // }


                repeated_from_to = make_range_number_cols(4, 16);

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

            } else if (analys_mode == "#single-crystal") {

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

            } else if (analys_mode == "#single-crystal") {

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

            if (analys_mode == "#two-crystal-groups") {
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
                          <th colspan="1" rowspan="2">Segment</th> \
                          <th colspan="1" rowspan="2">Positions</th> \
                          <th colspan="3" rowspan="2"> Frequency (%)</th> \
                          <th rowspan="2">Interactions</th> \
                          <th rowspan="2">Distance (Ca atoms)*</th> \
                          <th colspan="4">Backbone movement (Ca-7TM axis)</th> \
                          <th colspan="6">Sidechain differences</th> \
                          <th colspan="2" rowspan="2">Position presence %</th> \
                          <th colspan="4">Secondary structure</th> \
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

            var class_seq_cons_diff = class_seq_cons[0]-class_seq_cons[1];

            // var types = v['types'].join(",<br>");
            const types = v['types'].map((t) => types_to_short[t]).join('|');
            var seg1 = data['segm_lookup'][gn1];
            var seg2 = data['segm_lookup'][gn2];
            var distance = v['distance'];
            var angles_1 = v['angles'][0];
            var angles_2 = v['angles'][1];

            var pos1_presence = v['pos1_presence'];
            var pos2_presence = v['pos2_presence'];
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
                      <td class="narrow_col"> </td>
                      <td class="narrow_col"> </td>
                      <td class="narrow_col"> </td>
                      <td class="narrow_col"> </td>
                      <td class="narrow_col">${class_seq_cons[0]}</td>
                      <td class="narrow_col">${class_seq_cons[1]}</td>
                      <td class="narrow_col">${class_seq_cons_diff}</td>
                    </tr>`;
            tbody.append(tr);
        });
    } else if (data['proteins'].length > 1) {
        thead = '<tr> \
                          <th colspan="1" rowspan="2">Segment</th> \
                          <th colspan="1" rowspan="2">Positions</th> \
                          <th colspan="1" rowspan="2"> Frequency (%)</th> \
                          <th rowspan="2">Interactions</th> \
                          <th rowspan="2">Distance (Ca atoms)*</th> \
                          <th colspan="4">Backbone movement (Ca-7TM axis)</th> \
                          <th colspan="6">Sidechain differences</th> \
                          <th colspan="2" rowspan="2">Position presence %</th> \
                          <th colspan="4">Secondary structure</th> \
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
        var pdbs = data['pdbs'].length
        $.each(data['interactions'], function(i, v) {
            var gn1 = i.split(",")[0]
            var gn2 = i.split(",")[1]
            var sfreq1 = Math.round(100 * v['pdbs'].length / pdbs);
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
                      <td class="narrow_col"> </td>
                      <td class="narrow_col"> </td>
                      <td class="narrow_col"> </td>
                      <td class="narrow_col"> </td>
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
                          <th colspan="1" rowspan="2">Segment</th> \
                          <th colspan="1" rowspan="2">Positions</th> \
                          <th colspan="1" rowspan="2">Positions GN</th> \
                          <th rowspan="2">Interaction</th> \
                          <th rowspan="2">Distance (Ca atoms)*</th> \
                          <th colspan="4">Backbone movement (Ca-7TM axis)</th> \
                          <th colspan="6">Sidechain differences</th> \
                          <th colspan="4">Secondary structure</th> \
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
            tr = `
                    <tr class="clickable-row filter_rows" id="${pos1},${pos2}">
                      <td class="dt-center">${seg1}-${seg2}</td>
                      <td class="dt-center"><span>${pos1}</span>-<span>${pos2}</span></td>
                      <td class="dt-center">${gn1}-${gn2}</td>
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
                      <td class="narrow_col"> </td>
                      <td class="narrow_col"> </td>
                      <td class="narrow_col"> </td>
                      <td class="narrow_col"> </td>
                      <td class="narrow_col">${class_seq_cons}</td>
                    </tr>`;
            tbody.append(tr);
        });
    }

    // table.on('click', '.clickable-row', function(event) {
    //   if($(this).hasClass('active')){
    //     $(this).removeClass('active'); 
    //     $(selector + " .secondary-table").find('tbody').empty();
    //   } else {
    //     $(this).addClass('active').siblings().removeClass('active');
    //     renderSecondary($(this).attr('id'))
    //   }
    // });

    console.timeEnd("RenderBrowser");
    gray_scale_table(table);
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
    var proteins_1 = data['proteins1'].length
    var proteins_2 = data['proteins2'].length
    var pdbs_1 = data['pdbs1'].length
    var pdbs_2 = data['pdbs2'].length
    if (data['proteins2']) {
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
                      <td class="narrow_col"> </td>
                      <td class="narrow_col"> </td>
                      <td class="narrow_col"> </td>
                      <td class="narrow_col"> </td>
                    </tr>`;
            // tbody.append(tr);
        });
        // insert natively for speed increase on Chrome
        tbody[0].innerHTML = tr_list;
    } else if (data['proteins'].length > 1) {

    } else {

    }

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
        console.log('render 4');
        var proteins_1 = data['proteins1'].length
        var proteins_2 = data['proteins2'].length
        var pdbs_1 = data['pdbs1'].length
        var pdbs_2 = data['pdbs2'].length
        thead = '<tr> \
                          <th colspan="1" rowspan="2">Segment</th> \
                          <th colspan="1" rowspan="2">Positions</th> \
                          <th colspan="5" rowspan="1">Secondary structure</th> \
                          <th colspan="9" rowspan="1">Residue angles</th> \
                          <th colspan="6" rowspan="1">Helix turn angle</th> \
                          <th colspan="5" rowspan="1">Seq consensus</th> \
                          <th colspan="2" rowspan="1">Class seq consensus</th> \
                        </tr> \
                        <tr> \
                          <th colspan="2">Consensus SS</th> \
                          <th colspan="3">Frequency (%)</th> \
                          <th colspan="3">Tau (N-Ca-C)</th> \
                          <th colspan="3">Phi (N(+1)-C-Ca-N)</th> \
                          <th colspan="3">Psi (C-Ca-N-C(-1))</th> \
                          <th colspan="3">Theta (Ca(+2)-Ca(+1)-Ca-Ca(-1))</th> \
                          <th colspan="3">Next Theta (Ca(+1)-CA-CA(-1)-CA(-2))</th> \
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
            if (seg=='ECL1' || seg=='ECL2') return true;

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

            tr = ''
            tr_list += `
                    <tr class="clickable-row filter_rows" id="${i}">
                      <td class="dt-center">${seg}</td>
                      <td class="dt-center">${i}</td>

                      <td class="narrow_col"></td>
                      <td class="narrow_col"></td>
                      <td class="narrow_col"></td>
                      <td class="narrow_col"></td>
                      <td class="narrow_col"></td>

                      <td class="narrow_col">${angles1[3]}</td>
                      <td class="narrow_col">${angles2[3]}</td>
                      <td class="narrow_col">${angles_diff[3]}</td>

                      <td class="narrow_col">${angles1[4]}</td>
                      <td class="narrow_col">${angles2[4]}</td>
                      <td class="narrow_col">${angles_diff[4]}</td>

                      <td class="narrow_col">${angles1[5]}</td>
                      <td class="narrow_col">${angles2[5]}</td>
                      <td class="narrow_col">${angles_diff[5]}</td>

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

    } else {

    }

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
    var proteins_1 = data['proteins1'].length
    var proteins_2 = data['proteins2'].length
    var pdbs_1 = data['pdbs1'].length
    var pdbs_2 = data['pdbs2'].length
    if (data['proteins2']) {
       thead = '<tr> \
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
    } else if (data['proteins'].length > 1) {

    } else {

    }

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
    console.time('Greyscale cells');
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
    console.timeEnd('Greyscale cells');
    console.log(cell_count, 'cells greyscaled');
    console.timeEnd('Greyscale', table.attr('id'));
}
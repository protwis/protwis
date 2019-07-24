function createBoxPlot(data, element, plottype) {
    var mode = get_current_mode();
    var layout = {};
    switch (mode) {
        case "single-crystal-group":

            var data = [{
                y: [0, 1, 1, 2, 3, 5, 8, 13, 21],
                boxpoints: 'all',
                jitter: 0.3,
                pointpos: -1.8,
                type: 'box',
                text: ['Text A', 'Text B', 'Text C', 'Text D', 'Text E'],
            }];

            break;
        case "single-crystal":

            var data = [{
                y: [0, 1, 1, 2, 3, 5, 8, 13, 21],
                boxpoints: 'all',
                jitter: 0.3,
                pointpos: -1.8,
                type: 'box',
                text: ['Text A', 'Text B', 'Text C', 'Text D', 'Text E'],
            }];

            break;
        case "two-crystal-groups":

            switch (plottype) {

                case "angles":
                    var rows = getDateFromTable(1, [1, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]);

                    var text_array = getColumn(rows, 0).concat(getColumn(rows, 0));
                    var y = Array(getColumn(rows, 1).length).fill('Pos 1').concat(Array(getColumn(rows, 1).length).fill('Pos 2'))

                    var trace1 = {
                        y: getColumn(rows, 1).concat(getColumn(rows, 2)),
                        x: y,
                        name: 'Distance',
                        type: 'box',
                        boxmean: false,
                        text: text_array
                    };

                    var trace2 = {
                        y: getColumn(rows, 3).concat(getColumn(rows, 4)),
                        x: y,
                        name: 'Rotation (Ca angle)',
                        type: 'box',
                        boxmean: false,
                        text: text_array
                    };


                    var trace3 = {
                        y: getColumn(rows, 5).concat(getColumn(rows, 6)),
                        x: y,
                        name: 'Rotamer',
                        type: 'box',
                        boxmean: false,
                        text: text_array
                    };

                    var trace4 = {
                        y: getColumn(rows, 7).concat(getColumn(rows, 8)),
                        x: y,
                        name: 'SASA',
                        type: 'box',
                        boxmean: false,
                        text: text_array
                    };

                    var trace5 = {
                        y: getColumn(rows, 9).concat(getColumn(rows, 10)),
                        x: y,
                        name: 'RSA',
                        type: 'box',
                        boxmean: false,
                        text: text_array
                    };

                    var trace6 = {
                        y: getColumn(rows, 11).concat(getColumn(rows, 12)),
                        x: y,
                        name: 'Presence',
                        type: 'box',
                        boxmean: false,
                        text: text_array
                    };


                    var data = [trace1, trace2, trace3, trace4, trace5, trace6];

                    var layout = {
                        title: 'Grouped Horizontal Box Plot',
                        xaxis: {
                            title: 'Angles',
                            zeroline: false
                        },
                        boxmode: 'group'
                    };

                    break;
                default:
                    var rows = getDateFromTable(1, [1, 2, 3, 4]);

                    var text_array = getColumn(rows, 0);
                    var y0 = getColumn(rows, 1);
                    var y1 = getColumn(rows, 2);
                    var y2 = getColumn(rows, 3);


                    var trace1 = {
                        y: y0,
                        type: 'box',
                        text: text_array,
                        name: 'Set 1',
                        boxpoints: 'all',
                        jitter: 0.3,
                        pointpos: -1.8
                    };

                    var trace2 = {
                        y: y1,
                        type: 'box',
                        text: text_array,
                        name: 'Set 2',
                        boxpoints: 'all',
                        jitter: 0.3,
                        pointpos: -1.8
                    };

                    var trace3 = {
                        y: y2,
                        type: 'box',
                        text: text_array,
                        name: 'Diff',
                        boxpoints: 'all',
                        jitter: 0.3,
                        pointpos: -1.8
                    };

                    var data = [trace1, trace2, trace3];
            }

            break;
    }


    Plotly.newPlot(element, data, layout);

}

function createBoxPlotResidue(data, element, plottype, limit_pdbs = false, aa = false) {
    var mode = get_current_mode();
    var layout = {};
    switch (mode) {
        case "single-crystal-group":

            var data = [{
                y: [0, 1, 1, 2, 3, 5, 8, 13, 21],
                boxpoints: 'all',
                jitter: 0.3,
                pointpos: -1.8,
                type: 'box',
                text: ['Text A', 'Text B', 'Text C', 'Text D', 'Text E'],
            }];

            break;
        case "single-crystal":

            var data = [{
                y: [0, 1, 1, 2, 3, 5, 8, 13, 21],
                boxpoints: 'all',
                jitter: 0.3,
                pointpos: -1.8,
                type: 'box',
                text: ['Text A', 'Text B', 'Text C', 'Text D', 'Text E'],
            }];

            break;
        case "two-crystal-groups":

            switch (plottype) {

                case "angles":

                    x = [];
                    ys = {};
                    pdbs = [];
                    pos = '';
                    pdbs = two_sets_pdbs1.concat(two_sets_pdbs2);
                    // If only use a subset of pdbs.
                    if (limit_pdbs) pdbs = limit_pdbs;
                    pdbs_shown = []
                    pdbs.forEach(function(pdb){
                        pdb = pdb.toUpperCase();
                        let d = data[pdb];
                        if (d.length > 0) {
                            pos = d[0];
                            pdbs_shown.push(pdb)
                            if (two_sets_pdbs1.includes(pdb)) {
                                x.push('Set 1');
                            } else if (two_sets_pdbs2.includes(pdb)) {
                                x.push('Set 2');
                            }
                            for (var i = 2; i < d.length; i++) {
                                if (!ys[i - 2]) ys[i - 2] = [];
                                ys[i - 2].push(d[i]);
                            }
                        }
                    });

                    if (aa) pos = pos+" "+aa;

                    names = ['core_distance','a_angle','outer_angle','tau','phi','psi', 'sasa', 'rsa','theta','hse']

                    var traces = [];
                    for (const key in ys) {
                        trace = {
                            y: ys[key],
                            x: x,
                            name: names[key],
                            type: 'box',
                            boxmean: true,
                            boxpoints: 'all',
                            text: pdbs_shown
                        };
                        traces.push(trace);
                    }

                    var data = traces;

                    var layout = {
                        title: 'Angles data for position '+pos,
                        xaxis: {
                            zeroline: false
                        },
                        yaxis: {
                            zeroline: false
                        },
                        boxmode: 'group'
                    };

                    break;
                default:
            }

            break;
    }

    Plotly.newPlot(element, data, layout);
}


function getDateFromTable(browser_tab, columns) {
    var mode = get_current_mode();
    var mode_tab = mode + "-tab";
    var tab_id = "#" + mode_tab + " .browser-table-" + browser_tab;
    console.log(tab_id);
    var rows = [];
    if ($.fn.DataTable.isDataTable(tab_id)) {
        var table = $(tab_id).DataTable();
        table.rows({
            filter: 'applied'
        }).data().each(function(i) {
            var temp = [];
            columns.forEach(function(c, ii) {
                temp.push(i[c]);
            });
            rows.push(temp);
        })
    }
    return rows;

}

function getColumn(anArray, columnNumber) {
    return anArray.map(function(row) {
        return row[columnNumber];
    });
}

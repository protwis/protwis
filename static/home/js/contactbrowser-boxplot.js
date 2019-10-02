function sum(a) {
    var s = 0;
    for (var i = 0; i < a.length; i++) s += a[i];
    return s;
} 
 
function degToRad(a) {
    return Math.PI / 180 * a;
}
 
function meanAngleDeg(a) {
    return 180 / Math.PI * Math.atan2(
        sum(a.map(degToRad).map(Math.sin)) / a.length,
        sum(a.map(degToRad).map(Math.cos)) / a.length
    );
}

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
                    var rows = getDateFromTable(1, [1, 7, 9, 11, 13, 15, 17, 8, 10, 12, 14, 16, 18]);

                    var pos_titles = [];
                    var values = [];
                    for (var i = 0; i < rows.length; i++) {
                        title = rows[i][0].split("-");
                        if (!pos_titles.includes(title[0])) {
                            values.push(rows[i].slice(1, 7));
                            pos_titles.push(title[0]);
                        }
                        if (!pos_titles.includes(title[1])) {
                            values.push(rows[i].slice(7));
                            pos_titles.push(title[1]);
                        }
                    }

                    var y = Array(values.length).fill('Angles')

                    var trace1 = {
                        y: getColumn(values, 0),
                        name: 'Distance to 7TM axis',
                        type: 'box',
                        boxmean: false,
                        text: pos_titles
                    };

                    var trace2 = {
                        y: getColumn(values, 1),
                        name: 'Rotation (Ca angle)',
                        type: 'box',
                        boxmean: false,
                        text: pos_titles
                    };

                    var trace3 = {
                        y: getColumn(values, 2),
                        name: 'Rotamer',
                        type: 'box',
                        boxmean: false,
                        text: pos_titles
                    };

                    var trace4 = {
                        y: getColumn(values, 3),
                        name: 'SASA',
                        type: 'box',
                        boxmean: false,
                        text: pos_titles
                    };

                    var trace5 = {
                        y: getColumn(values, 4),
                        name: 'RSA',
                        type: 'box',
                        boxmean: false,
                        text: pos_titles
                    };

                    var trace6 = {
                        y: getColumn(values, 5),
                        name: 'Presence',
                        type: 'box',
                        boxmean: false,
                        text: pos_titles
                    };

                    var data = [trace1, trace2, trace3, trace4, trace5, trace6];

                    var layout = {
                        title: 'Grouped Horizontal Box Plot',
                        xaxis: {
                            // title: 'Angles',
                            zeroline: false
                        },
                        // grid: {rows: 1, columns: 6, pattern: 'independent'},
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

function createBoxPlotResidue(data, element, plottype, cell_index, limit_pdbs = false, aa = false) {
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
                    ys1 = {};
                    ys2 = {};
                    pdbs = [];
                    pos = '';
                    pdbs = two_sets_pdbs1.concat(two_sets_pdbs2);
                    // If only use a subset of pdbs.
                    if (limit_pdbs) pdbs = limit_pdbs;
                    pdbs_shown = []
                    pdbs.forEach(function(pdb) {
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
                                if (!ys[i - 2]) {
                                    ys[i - 2] = [];
                                    ys1[i - 2] = [];
                                    ys2[i - 2] = [];
                                }
                                ys[i - 2].push(d[i]);

                                if (two_sets_pdbs1.includes(pdb))
                                    ys1[i - 2].push(d[i]);

                                if (two_sets_pdbs2.includes(pdb))
                                    ys2[i - 2].push(d[i]);
                            }
                        }
                    });

                    if (aa) pos = pos + " " + aa;

                    names = [['core_distance',false], ['a_angle',true], ['outer_angle',true], ['tau',true], ['phi',true], ['psi',true], ['sasa',false], ['rsa',false], ['theta',true], ['hse',false], ['dssp',false]]

                    name_index = Math.floor((cell_index - 6) / 2);
                    name_index = {
                        7: 'core_distance',
                        8: 'core_distance',
                        9: 'a_angle',
                        10: 'a_angle',
                        11: 'outer_angle',
                        12: 'outer_angle',
                        13: 'sasa',
                        14: 'sasa',
                        15: 'rsa',
                        16: 'rsa'
                    };

                    // console.log(cell_index,name_index,name_index[cell_index]);

                    var traces = [];
                    is_angle = true;
                    for (const key in ys) {
                        if (name_index[cell_index] != names[key][0]) continue;
                        is_angle = names[key][1]
                        new_x = []
                        for (var i = 0; i < x.length; i++) {
                            new_x.push(names[key][0] + "<br>" + x[i]);
                        }
                        visible_trace = 'legendonly';
                        if (name_index[cell_index] == names[key][0]) visible_trace = true;
                        if (names[key][1] == true) {
                            mean1 = meanAngleDeg(ys1[key]);
                            mean2 = meanAngleDeg(ys2[key]);
                            traces = [{
                                    type: "scatterpolar",
                                    mode: "markers",
                                    r: Array(ys1[key].length).fill(1),
                                    theta: ys1[key],
                                    name: "Set 1",
                                    marker: {
                                        color: "#8090c7",
                                        symbol: "circle",
                                        size: 15
                                    },
                                    opacity: 0.8,
                                },
                                {
                                    type: "scatterpolar",
                                    mode: "lines",
                                    r: [0, 0.9],
                                    theta: [0, mean1],
                                    name: "Set 1 mean",
                                    line: {
                                        color: "#8090c7",
                                        width: 5,
                                    },
                                },
                                {
                                    type: "scatterpolar",
                                    mode: "markers",
                                    r: Array(ys2[key].length).fill(1),
                                    theta: ys2[key],
                                    name: "Set 2",
                                    marker: {
                                        color: "red",
                                        symbol: "circle",
                                        size: 15
                                    },
                                    opacity: 0.8,
                                    subplot: "polar2"
                                },
                                {
                                    type: "scatterpolar",
                                    mode: "lines",
                                    r: [0, 0.9],
                                    theta: [0, mean2],
                                    name: "Set 2 mean",
                                    line: {
                                        color: "red",
                                        width: 5,
                                    },
                                    subplot: "polar2"
                                }];



                        } else {
                            trace = {
                                y: ys[key],
                                x: new_x,
                                name: names[key][0],
                                type: 'box',
                                boxmean: true,
                                // boxpoints: 'all',
                                text: pdbs_shown,
                                jitter: 0.5,
                                whiskerwidth: 0.2,
                                fillcolor: 'cls',
                                marker: {
                                    size: 2
                                },
                                line: {
                                    width: 1
                                },
                                visible: visible_trace
                            };
                            // Only show one data type on a chart.
                            traces.push(trace);

                        }
                    }

                    var data = traces;
                    if (is_angle == true) {

                        var layout = {
                            title: 'Angles data for position ' + pos,
                            polar: {
                                radialaxis: {
                                    //autorange: true,
                                    showgrid: true,
                                    zeroline: false,
                                    showline: false,
                                    autotick: true,
                                    ticks: '',
                                    showticklabels: false,
                                    range: [0, 1.1]

                                },
                                angularaxis: {
                                    tickfont: {
                                        size: 8
                                    }
                                }
                                },
                            polar2: {
                                radialaxis: {
                                    //autorange: true,
                                    showgrid: true,
                                    zeroline: false,
                                    showline: false,
                                    autotick: true,
                                    ticks: '',
                                    showticklabels: false,
                                    range: [0, 1.1]

                                },
                                angularaxis: {
                                    tickfont: {
                                        size: 8
                                    }
                                }
                                }
                        }
                    } else {
                        var layout = {
                            title: 'Angles data for position ' + pos,
                            // xaxis: {
                            //     zeroline: false
                            // },
                            xaxis: {
                                showgrid: false,
                                zeroline: false,
                            },
                            yaxis: {
                                zeroline: false
                            },
                            // boxmode: 'group'
                        };
                    }
                    break;
                default:
            }

            break;
    }
    Plotly.newPlot(element, data, layout, { showLink: true }); //editable: true,
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
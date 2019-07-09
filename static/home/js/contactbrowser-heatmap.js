function renderHeatmap(data, heatMapSelector) {

    var start = new Date().getTime();
    // Destroy old zoom on heatmap
    if (window.zoomHeatmap[heatMapSelector] != null) {
        window.zoomHeatmap[heatMapSelector].destroy();
        delete window.zoomHeatmap[heatMapSelector];
    }

    var end = new Date().getTime();
    // console.log('Remove old zoom',end-start);

    // Destroy old legend content
    // $(heatMapSelector + ' .heatmap-legend').empty();

    // Destroy all previous contents
    var heatmap_id = $(heatMapSelector + ' .heatmap').attr("id");
    $("#" + heatmap_id + "_content").remove();
    $(heatMapSelector + ' .heatmap-interaction').remove();
    $(heatMapSelector + ' .heatmap').empty();
    var heatmap = Snap(heatMapSelector + ' .heatmap');
    var heatmap_id = $(heatMapSelector + ' .heatmap').attr("id");
    console.log('ID is', heatmap_id, 'from', heatMapSelector);

    if (data.length == 3) {
        data1 = data[0];
        data2 = data[1];
        data3 = data[2];
        var interactions = data3.interactions;
        var segment_map = data3.segment_map;
        var sequence_numbers = data3.sequence_numbers;
        var num_seq_numbers = Object.keys(data3.sequence_numbers).length;
    } else {
        // Draw heatmap
        var interactions = data.interactions;
        var segment_map = data.segment_map;
        var sequence_numbers = data.sequence_numbers;
        var aa_map = data.aa_map[Object.keys(data.aa_map)[0]];
        var gen_map = data.generic_map;
        var num_seq_numbers = Object.keys(data.sequence_numbers).length;

        var pdbs = data.pdbs;
        var pdbs1 = data.pdbs1;
        var pdbs2 = data.pdbs2;
    }
    x = 0;
    wi = num_seq_numbers;
    y = 0;
    hi = num_seq_numbers;

    heatmap.attr({ viewBox: [x, y, wi, hi].join(',') });
    heatmap.attr({ viewBox: [x, y, wi, hi].join(' ') });

    // Contains all labels
    var labelGroup = heatmap.g();

    // Contains all content
    var contentGroup = heatmap.g();
    contentGroup.attr('id', heatmap_id + '_content');

    // Compute segment offsets
    var i;

    var segments = [];

    var seg, prevSeg = segment_map[sequence_numbers[0]];
    var seqStart = 0;

    for (i = 0; i < num_seq_numbers; i++) {
        seg = segment_map[sequence_numbers[i]];

        if (seg === prevSeg) {
            continue;
        }

        segments.push({
            seg: prevSeg,
            start: seqStart,
            end: i - 1
        });

        seqStart = i;
        prevSeg = seg;
    }

    // Push last segment
    segments.push({
        seg: prevSeg,
        start: seqStart,
        end: i - 1
    });

    // Determine font-size
    label_font_size = Math.round(1 * num_seq_numbers / 20);
    label_font_size = "5";
    console.log(num_seq_numbers, 'num_seq_numbers', label_font_size);

    lines = Array();
    // Draw segments
    last_end = 0
    segments.forEach(function(s) {
        var rgb = getSegmentColor(s.seg);
        // Place the segments vertically.
        var cell = heatmap.rect(s.start, 0, s.end - s.start + 1, 0);
        var line = heatmap.line(s.start, 0, s.start, num_seq_numbers);
        last_end = s.end;
        cell.attr({
            'fill': "rgb(" + [rgb.r, rgb.g, rgb.b].join(',') + ")",
            // 'fill': "#FFFFFF",
            'fill-opacity': "0.5"
        });

        line.attr({
            'stroke': "rgb(150,150,150)",
            'strokeWidth': "0.1"
        });


        // Add text
        var label = heatmap.text(s.start + (s.end - s.start + 1) / 2 - label_font_size, -1 + cell.getBBox().height, s.seg);
        label.attr({
            'text-anchor': 'start',
            'font-size': label_font_size
        });

        // contentGroup.add(cell);
        lines.push(line);
        labelGroup.add(label);

        // label.transform("r270s-1,1");
    });

    var line = heatmap.line(last_end + 1, 0, last_end + 1, num_seq_numbers);
    line.attr({
        'stroke': "rgb(150,150,150)",
        'strokeWidth': "0.1"
    });
    lines.push(line);
    var line = heatmap.line(0, 0, num_seq_numbers, 0);
    line.attr({
        'stroke': "rgb(150,150,150)",
        'strokeWidth': "0.1"
    });
    lines.push(line);

    segments.forEach(function(s) {
        rgb = getSegmentColor(s.seg);

        // Place the segments horizontally.
        cell = heatmap.rect(0, s.start, num_seq_numbers, s.end - s.start + 1);
        line = heatmap.line(0, s.end + 1, num_seq_numbers, s.end + 1);

        cell.attr({
            'fill': "rgb(" + [rgb.r, rgb.g, rgb.b].join(',') + ")",
            // 'fill': "#FFFFFF",
            'fill-opacity': "0.5"
        });

        line.attr({
            'stroke': "rgb(150,150,150)",
            'strokeWidth': "0.1"
        });

        var label = heatmap.text(-1, 1 + s.start + (s.end - s.start) / 2, s.seg);

        label.attr({
            'text-anchor': 'end',
            'font-size': label_font_size,
            'alignment-baseline': 'middle'
        });

        // label.transform("r270");

        contentGroup.add(cell);
        lines.push(line);
        labelGroup.add(label);
    });
    // console.log('segments drawn')
    var start = new Date().getTime();
    cells = Array();
    rects = Array();
    // Draw cells
    var svgns = "http://www.w3.org/2000/svg";
    // console.log(relevant_interactions);
    if (pdbs.length == 1) {
        heatmap_mode = 'single';
    } else if (pdbs.length > 1 && !(pdbs2)) {
        heatmap_mode = 'single_group';
    } else if (pdbs2.length > 1) {
        heatmap_mode = 'two_groups';
    } 
    console.log(heatmap_mode, 'pdbs', pdbs, 'pdbs1', pdbs1, 'pdbs2', pdbs2);
    for (i = 0; i < num_seq_numbers; i++) {
        for (var j = 0; j < num_seq_numbers; j++) {

            // Get the sequence numbers
            var seq_i = data.sequence_numbers[i];
            var seq_j = data.sequence_numbers[j];

            // Only draw if an interaction exists
            var num = seq_i + "," + seq_j;
            var num2 = seq_j + "," + seq_i;

            if (num2 in interactions) num = num2;

            // console.log(num,'num of freq',interactions.length);
            if (num in interactions) {

                if (heatmap_mode == 'single') {
                    // getInteractionTypesFromPdbObject(interactions[num]).forEach(function(interaction) {
                    interaction = interactions[num];
                    var rgb = getInteractionColor(interaction.types[0]);
                    // var cell = heatmap.rect(i, j, 1, 1);

                    var rect = document.createElementNS(svgns, 'rect');
                    rect.setAttributeNS(null, 'x', i);
                    rect.setAttributeNS(null, 'y', j);
                    rect.setAttributeNS(null, 'height', '1');
                    rect.setAttributeNS(null, 'width', '1');

                    var interactionsString = interaction.types.join(", ");

                    var aa_i = aa_map[seq_i];
                    var aa_j = aa_map[seq_j];

                    var seq_pos_i = interaction.seq_pos[0];
                    var seq_pos_j = interaction.seq_pos[1];

                    var content, title = 'Residues ' + aa_map[seq_i] + seq_pos_i + '-' + aa_map[seq_j] + seq_pos_j + '<br />' +
                        'Interactions: ' + interactionsString + '<br />' +
                        'Segments: ' + segment_map[seq_i] + ', ' + segment_map[seq_j] + '<br />';

                    // // Add generic numbers where applicable
                    // if (seq_i in gen_map) {
                    //     title += 'Res. 1 gen. no: ' + gen_map[seq_i] + '<br />';
                    // }

                    // if (seq_j in gen_map) {
                    //     title += 'Res. 2 gen. no: ' + gen_map[seq_j] + '<br />';
                    // }

                    // var genStrI = gen_map[seq_i];
                    // var genStrJ = gen_map[seq_j];

                    // if (typeof gen_map[seq_i] == 'undefined') {
                    //     genStrI = '-';
                    // }

                    // if (typeof gen_map[seq_j] == 'undefined') {
                    //     genStrJ = '-';
                    // }

                    var popoverTable = '<table class="table">' +
                        '<thead>' +
                        '<tr>' +
                        '<th>Residue</th>' +
                        '<th>' + aa_map[seq_i] + seq_pos_i + '</th>' +
                        '<th>' + aa_map[seq_j] + seq_pos_j + '</th>' +
                        '</tr>' +
                        '</thead>' +
                        '<tbody>' +
                        '<td>Segment</td>' +
                        '<td>' + segment_map[seq_i] + '</td>' +
                        '<td>' + segment_map[seq_j] + '</td>' +
                        '</tr>' +
                        '<tr>' +
                        '<td>Gen. no.</td>' +
                        '<td>' + seq_i + '</td>' +
                        '<td>' + seq_j + '</td>' +
                        '</tr>' +
                        '</tbody>' +
                        '</table>' +
                        '<table class="table">' +
                        '<thead>' +
                        '<tr>' +
                        '<th>Interactions</th>' +
                        '</tr>' +
                        '<tr>' +
                        '</thead>' +
                        '<tbody>' +
                        '<tr>' +
                        '<td>' + interactionsString + '</td>' +
                        '</tr>' +
                        '</tbody>' +
                        '</table>';

                    rect.setAttributeNS(null, 'fill', "rgb(" + [rgb.r, rgb.g, rgb.b].join(',') + ")");
                    rect.setAttributeNS(null, 'class', "heatmap-interaction");
                    rect.setAttributeNS(null, 'title', "Interaction between " + seq_i + ", " + seq_j);
                    rect.setAttributeNS(null, 'data-content', popoverTable);
                    rect.setAttributeNS(null, 'data-res-no-1', seq_i);
                    rect.setAttributeNS(null, 'data-res-no-2', seq_j);
                    // rect.setAttributeNS(null, 'data-gen-no-1', genStrI);
                    // rect.setAttributeNS(null, 'data-gen-no-2', genStrJ);
                    rect.setAttributeNS(null, 'data-interaction-type', interaction);
                    document.getElementById(heatmap_id).appendChild(rect);
                    // });

                } else if (heatmap_mode == 'single_group') {
                    var rect = document.createElementNS(svgns, 'rect');
                    rect.setAttributeNS(null, 'x', i);
                    rect.setAttributeNS(null, 'y', j);
                    rect.setAttributeNS(null, 'height', '1');
                    rect.setAttributeNS(null, 'width', '1');

                    var nInteractions = Object.keys(interactions[num].pdbs).length;
                    var frequency = nInteractions / data.pdbs.length;

                    var rgb = { r: 255, g: Math.round(255 - frequency * 255), b: Math.round(255 - frequency * 255) };

                    var popoverTable = '<table class="table">' +
                        '<thead>' +
                        '<tr>' +
                        '<th>Residue #</th>' +
                        '<th>' + seq_i + '</th>' +
                        '<th>' + seq_j + '</th>' +
                        '</tr>' +
                        '</thead>' +
                        '<tbody>' +
                        '<td>Segment</td>' +
                        '<td>' + segment_map[seq_i] + '</td>' +
                        '<td>' + segment_map[seq_j] + '</td>' +
                        '</tr>' +
                        '</tbody>' +
                        '</table>' +
                        'Interaction count: ' + nInteractions + '<br />' +
                        'Interaction frequency: ' + frequency.toFixed(2)


                    rect.setAttributeNS(null, 'fill', "rgb(" + [rgb.r, rgb.g, rgb.b].join(',') + ")");
                    rect.setAttributeNS(null, 'class', "heatmap-interaction");
                    rect.setAttributeNS(null, 'title', "Interaction between " + seq_i + ", " + seq_j);
                    rect.setAttributeNS(null, 'data-content', popoverTable);
                    rect.setAttributeNS(null, 'data-gen-no-1', seq_i);
                    rect.setAttributeNS(null, 'data-gen-no-2', seq_j);
                    rect.setAttributeNS(null, 'data-frequency', frequency);
                    // rect.setAttributeNS(null, 'data-extra', JSON.stringify(interactions[num][2]));
                    document.getElementById(heatmap_id).appendChild(rect);
                } else if (heatmap_mode == 'two_groups') {

                    // Only draw if an interaction exists
                    var num = seq_i + "," + seq_j;
                    var num2 = seq_j + "," + seq_i;

                    if (num2 in interactions) num = num2;
                    // console.log(interactions[num]);

                    // // Difference in frequencies
                    var f1 = (interactions[num].pdbs1.length / pdbs1.length);
                    var f2 = (interactions[num].pdbs2.length / pdbs2.length);
                    var fDiff = f1 - f2;

                    var rgb = getGradientColor(fDiff, true);

                    var rect = document.createElementNS(svgns, 'rect');
                    rect.setAttributeNS(null, 'x', i);
                    rect.setAttributeNS(null, 'y', j);
                    rect.setAttributeNS(null, 'height', '1');
                    rect.setAttributeNS(null, 'width', '1');

                    var title = 'Residues ' + seq_i + ', ' + seq_j + '<br />' +
                        'Frequency group 1: ' + f1.toFixed(2) + '<br />' +
                        'Frequency group 2: ' + f2.toFixed(2) + '<br />' +
                        'Frequency difference: ' + fDiff.toFixed(2) + '<br />';

                    var popoverTable = '<table class="table">' +
                        '<thead>' +
                        '<tr>' +
                        '<th>Residue #</th>' +
                        '<th>' + seq_i + '</th>' +
                        '<th>' + seq_j + '</th>' +
                        '</tr>' +
                        '</thead>' +
                        '<tbody>' +
                        '<td>Segment</td>' +
                        '<td>' + segment_map[seq_i] + '</td>' +
                        '<td>' + segment_map[seq_j] + '</td>' +
                        '</tr>' +
                        '</tbody>' +
                        '</table>' +
                        'Group 1 freq: ' + f1.toFixed(2) + '<br />' +
                        'Group 2 freq: ' + f2.toFixed(2) + '<br />' +
                        'Frequency difference: ' + fDiff.toFixed(2)

                    rect.setAttributeNS(null, 'fill', "rgb(" + [rgb.r, rgb.g, rgb.b].join(',') + ")");
                    rect.setAttributeNS(null, 'class', "heatmap-interaction");
                    rect.setAttributeNS(null, 'title', "Interaction between " + seq_i + ", " + seq_j);
                    rect.setAttributeNS(null, 'data-content', popoverTable);
                    rect.setAttributeNS(null, 'data-frequency-diff', fDiff.toFixed(2));
                    rect.setAttributeNS(null, 'data-gen-no-1', seq_i);
                    rect.setAttributeNS(null, 'data-gen-no-2', seq_j);
                    document.getElementById(heatmap_id).appendChild(rect);



                }
            }
        }
    }

    // Start new g() here so it comes on top of cells
    var contentLines = heatmap.g();
    contentLines.add(lines);
    contentLines.add(labelGroup);

    var end = new Date().getTime();

    $(heatMapSelector + ' rect.heatmap-interaction').click(function() {
        // alert( "Handler for .click() called." );
        var $this = $(this);
        //if not already initialized
        if (!$this.data('bs.popover')) {
            $this.popover({
                'container': heatMapSelector,
                'placement': 'bottom',
                'animation': true,
                'html': true,
                'tabindex': '0'
            }).popover('show');
        }
    });


    var end2 = new Date().getTime();


    // // Make zoomable
    window.zoomHeatmap[heatMapSelector] = svgPanZoom(heatMapSelector + ' .heatmap', {
        zoomEnabled: true,
        // controlIconsEnabled: true,
        fit: true,
        center: true,
        minZoom: 0.75,
        maxZoom: 50,
        zoomScaleSensitivity: 0.25,
        dblClickZoomEnabled: true,
        beforeZoom: hidePopovers,
        beforePan: hidePopovers
    });

    // // Set initial zoom level

    // Close popovers on clicking elsewhere
    $('html').on('mousedown', function(e) {
        if (!$(e.target).closest('.popover').length) {
            if ($(e.target).closest(heatMapSelector).length) {
                hidePopovers();
            }
        }
    });
    old_heatmap_width = 0;

}
function createFlareplot_segment(raw_data,width, inputGraph, containerSelector) {

    var w = width;
    var h = w;

    var rx = w * 1 * 0.5;
    var ry = rx;
    innerRadius = rx*0.9;
    outerRadius = rx;


    var cx = w * 0.5;
    var cy = cx;


    var opacityValueBase = 0.7;
    var opacityValue = 1;
    
    
    /*//////////////////////////////////////////////////////////
    ////////////////// Set up the Data /////////////////////////
    //////////////////////////////////////////////////////////*/


    var protein_segments = ["N-term", "TM1", "ICL1", "TM2", "ECL1", "TM3", "ICL2", "TM4", "ECL2", "TM5", "ICL3", "TM6", "ECL3", "TM7", "H8", "C-term"];
    
    var matrix = [];
    var colors = [];
    // prepare colors and matrix arrays
    protein_segments.forEach(function (key, index) {
        colors.push(assignRainbowColor(key))
        matrix.push(Array(protein_segments.length).fill(0))
    });

    // Populate matrix for interactions between segments
    $.each(raw_data['interactions'], function (i, v) {

        if (!filtered_gn_pairs.includes(i)) return;
        gns = separatePair(i);
        seg1 = raw_data['segment_map'][gns[0]];
        seg2 = raw_data['segment_map'][gns[1]];
        if (seg1 != seg2) {
            seg1_i = protein_segments.indexOf(seg1);
            seg2_i = protein_segments.indexOf(seg2);
            matrix[seg1_i][seg2_i] += 1;
            matrix[seg2_i][seg1_i] += 1;
        }
    });

    // Purge segments with no counts..
    protein_segments.slice().reverse().forEach(function (key, index, object) {
        index = protein_segments.indexOf(key);
        summed = matrix[index].reduce((a, b) => a + b, 0);
        if (summed == 0) {
            matrix.splice(index, 1);
            colors.splice(index, 1);
            // object.splice(index, 1);
            matrix.forEach(a => a.splice(index, 1));
        }
    });
    

    /*Initiate the color scale*/
    var fill = d3.scale.ordinal()
        .domain(d3.range(protein_segments.length))
        .range(colors);
        
    /*//////////////////////////////////////////////////////////
    /////////////// Initiate Chord Diagram /////////////////////
    //////////////////////////////////////////////////////////*/

    d3.select(containerSelector).style("position","relative");

    div = d3.select(containerSelector).insert("div")
        .attr("class", "flareplot")
        .style("width", "100%")
        // .style("margin-top","100px")
        // .style("height", h + "px")
        .style("-webkit-backface-visibility", "hidden");
    

    svg = div.append("svg:svg")
        .attr("viewBox", "0 0 " + w + " " + h )
        .attr("width", "100%")
        .attr("style", "height: 500px")
        .attr("class", "flareplot")
        .append("svg:g")
        .attr("transform", "translate(" + cx + "," + cy + ")");
    
    var chord = d3.layout.chord()
        .padding(.05)
        // .sortSubgroups(d3.descending) /*sort the chords inside an arc from high to low*/
        .sortChords(d3.ascending) /*which chord should be shown on top when chords cross. Now the biggest chord is at the bottom*/
        .matrix(matrix);
        

    /*//////////////////////////////////////////////////////////
    ////////////////// Draw outer Arcs /////////////////////////
    //////////////////////////////////////////////////////////*/

    var arc = d3.svg.arc()
        .innerRadius(innerRadius)
        .outerRadius(outerRadius);
        
    var g = svg.selectAll("g.group")
        .data(chord.groups)
        .enter().append("svg:g")
        .attr("class", function(d) {return "group " + protein_segments[d.index];})
        .on("mouseover", fade(.2))
        .on("mouseout", fade(opacityValue));
        
    g.append("svg:path")
        .attr("class", "arc")
        .style("stroke", function(d) { return fill(d.index); })
        .style("fill", function(d) { return fill(d.index); })
        .attr("d", arc)
        .style("opacity", 1);

        
    /*//////////////////////////////////////////////////////////
    ////////////////// Initiate Names //////////////////////////
    //////////////////////////////////////////////////////////*/

    g.append("svg:text")
    .each(function(d) { d.angle = (d.startAngle + d.endAngle) / 2; })
    .attr("dy", ".35em")
    .attr("class", "segmentElement")
    .attr("transform", function(d) {
            return "rotate(" + (d.angle * 180 / Math.PI - 90) + ")"
            + "translate(" + (innerRadius + 20) + ")"
            + "rotate(90)";
    })
    .attr('opacity', 1)
    .style('fill','white')
    .style("text-anchor", "middle")
    .style("font-size", "30px")
    .text(function(d,i) { return protein_segments[i]; });  

    /*//////////////////////////////////////////////////////////
    //////////////// Initiate inner chords /////////////////////
    //////////////////////////////////////////////////////////*/

    var chords = svg.selectAll("path.chord")
        .data(chord.chords)
        .enter().append("svg:path")
        .attr("class", "chord")
        .style("stroke", function(d) { return d3.rgb(fill(d.source.index)).darker(); })
        .style("fill", function(d) { return fill(d.source.index); })
        .attr("d", d3.svg.chord().radius(innerRadius))
        .attr('opacity', opacityValue);


    /*Make mouse over and out possible*/
    // d3.selectAll(".group")
    // .on("mouseover", fade(.02))
    // .on("mouseout", fade(.80));
        
    /*Show all the text*/
    d3.selectAll("g.group").selectAll("line")
        .style("stroke", "#000");
    
    // /*And the Names of each Arc*/	
    // svg.selectAll("g.group")
    //     .transition().duration(100)
    //     .selectAll(".segmentElement").style("opacity",1);		



    /*//////////////////////////////////////////////////////////
    ////////////////// Extra Functions /////////////////////////
    //////////////////////////////////////////////////////////*/

    /*Returns an event handler for fading a given chord group*/
    function fade(opacity) {
    return function(d, i) {
        svg.selectAll("path.chord")
            .filter(function(d) { return d.source.index != i && d.target.index != i; })
            .transition()
            .style("stroke-opacity", opacity)
            .style("fill-opacity", opacity);
    };
    };/*fade*/

    function assignRainbowColor(segment) {
        var segmentRainbowColors2 = {
            "1": "#736DA7",
            "2": "#5EB7B7",
            "3": "#CE9AC6",
            "4": "#DD7D7E",
            "5": "#E6AF7C",
            "6": "#DEDB75",
            "7": "#80B96F",
            "8": "#000000",
            "0": "#EEE"
        };
        var color = "";
        switch (segment) {
        case "TM1":
            color = segmentRainbowColors2["1"];
            break;
        case "TM2":
            color = segmentRainbowColors2["2"];
            break;
        case "TM3":
            color = segmentRainbowColors2["3"];
            break;
        case "TM4":
            color = segmentRainbowColors2["4"];
            break;
        case "TM5":
            color = segmentRainbowColors2["5"];
            break;
        case "TM6":
            color = segmentRainbowColors2["6"];
            break;
        case "TM7":
            color = segmentRainbowColors2["7"];
            break;
        case "H8":
            color = segmentRainbowColors2["8"];
            break;
        default:
            color = segmentRainbowColors2["0"];
        }
        return color;
    }
    
    function separatePair(stringPair) {
        var regex = /([0-9x]+),([0-9x]+)/g;
        var m;

        matches = regex.exec(stringPair);

        return [matches[1], matches[2]];
    }
}
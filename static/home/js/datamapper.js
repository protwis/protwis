// HEATMAP REPRESENTATION

function simple_heatmap(data, location, element_id, legend_label) {
  // Set the dimensions and margins of the graph
  var margin = { top: 30, right: 100, bottom: 30, left: 60 }; // Increased left margin to accommodate row labels

  // Extract rows and columns from the new data format
  var rows = Object.keys(data);
  var cols = Object.keys(data[rows[0]]);

  // Calculate the width required for row labels
  var rowLabelWidth = d3.max(rows, function(d) {
    return d.length * 25; // Adjust the multiplier as needed
  });

  // Prepare chart data
  var chartData = [];
  var highest_value = 0;

  for (var row in data) {
    for (var col in data[row]) {
      chartData.push({ row: row, col: col, value: data[row][col] });
      if (data[row][col] > highest_value) {
        highest_value = data[row][col];
      }
    }
  }

  var width = (45 * cols.length) + rowLabelWidth; // Adjust the multiplier as needed
  var height = (20 * rows.length);

  // Append the SVG object to the body of the page
  var svg_home = d3.select("#" + location)
    .append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + (margin.top * 5))
    .attr("transform", "translate(0,-" + margin.bottom + ")")
    .attr("id", element_id);

  var redscale = ['#ffcccc', '#ff0000']; // Changed from greyscale

  var legend = svg_home.append('defs')
    .append('linearGradient')
    .attr('id', 'grad_' + element_id)
    .attr('x1', '0%')
    .attr('x2', '100%')
    .attr('y1', '0%')
    .attr('y2', '0%');

  legend.selectAll('stop')
    .data(redscale) // Changed from greyscale
    .enter()
    .append('stop')
    .style('stop-color', function(d) { return d; })
    .attr('offset', function(d, i) {
      return 100 * (i / (redscale.length - 1)) + '%'; // Changed from greyscale
    });


  var legend_svg = svg_home.append("g")
    .attr("transform", "translate(30," + (height + 100) + ")");

  var color_svg = svg_home.append("g")
    .attr("transform", "translate(" + (margin.left * 1.5) + ",0)");

  var svg = svg_home.append("g")
    .attr("transform", "translate(" + (margin.left * 2) + ",0)");

  legend_svg.append("text")
    .attr('y', 10)
    .attr('x', -50)
    .style("font", "14px sans-serif")
    .style("font-weight", "bold")
    .text(legend_label);

  legend_svg.append("text")
    .attr('x', 30)
    .attr('y', 30)
    .style("font", "14px sans-serif")
    .style("font-weight", "bold")
    .text('0');

  legend_svg.append('rect')
    .attr('x', 40)
    .attr('y', 15)
    .attr('width', (width * 0.7) + margin.left + margin.right)
    .attr('height', 20)
    .style('fill', 'url(#grad_' + element_id + ')');

  legend_svg.append("text")
    .attr('x', width + margin.left)
    .attr('y', 30)
    .style("font", "14px sans-serif")
    .style("font-weight", "bold")
    .text(highest_value);

  // Using each to ensure the text element is fully created
  var legendLabel;
  legend_svg.selectAll("text").each(function() {
    legendLabel = this.getBBox().width * 1.05 + 0.5 * 10;
  });

  legend_svg.select("text")
    .attr("x", (width + margin.left + margin.right - legendLabel) / 2 + margin.left - 120);

  // Build X scales and axis
  var x = d3.scale.ordinal()
    .rangeBands([0, width], 0.01)
    .domain(cols);

  svg.append("g")
    .attr("transform", "translate(0," + height + ")")
    .attr('id', 'Xaxis')
    .call(d3.svg.axis().scale(x).orient("bottom").tickSize(0))
    .selectAll("text")
    .style("text-anchor", "end")
    .attr("transform", "rotate(-45)")
    .attr("dx", "-0.8em")
    .attr("dy", "0.15em")
    .attr("class", "column-label"); // Add a class for styling

  // Build Y scales and axis
  var y = d3.scale.ordinal()
    .rangeBands([height, 0], 0.01)
    .domain(rows);

  svg.append("g")
    .attr('id', 'Yaxis')
    .call(d3.svg.axis().scale(y).orient("left").tickSize(0));


  // Build color scale
  var myColor = d3.scale.linear()
    .range(["#ffcccc", "#ff0000"]) // Changed from ["white", "black"]
    .domain([1, highest_value]);

  // Read the data
  svg.selectAll("rect")
    .data(chartData, function(d) { return d.row + ':' + d.col; })
    .enter()
    .append("rect")
    .attr("x", function(d) { return x(d.col); })
    .attr("y", function(d) { return y(d.row); })
    .attr("width", x.rangeBand())
    .attr("height", y.rangeBand())
    .style("fill", function(d) { return myColor(d.value); });

  svg.select('#Xaxis')
    .attr('text-anchor', 'start');

  d3.selectAll("#Yaxis>.tick>text")
    .each(function(d, i) {
      d3.select(this).style("font-size", "1.3em");
    });

  d3.selectAll("#Xaxis>.tick>text")
    .each(function(d, i) {
      d3.select(this).style("font-size", "1.3em");
    });

  var count = 1;
  var ticks = svg.select('#Xaxis').selectAll('.tick');
  ticks.each(function(d) {
    var value = (count & 1) ? "odd" : "even";
    var text = d3.select(this).select('text');
    var textSize = Math.floor(text.node().getBBox().width * 1.05 + 0.5 * 10);
    text.attr("transform", null);
    // text.attr("x", "-" + textSize / 2 + "px");
    text.attr("x","42")
    text.attr("y","15")
    // if (value === "even") {
    //   text.attr("y", "0");
    // }
    count = count + 1;
  });
}


// ####################
// ### LIST PLOT ######
// ####################

// Initialization of global functions

function initializeDataStyling(list_data_wow, data_types) {

    // Function to check if any of the specified keys are present in a given object
    function checkKeysPresent(obj, keys) {
        return keys.some(k => k in obj);
    }

    // Define the keys to check for each column
    var col1Keys = ["Value1", "Value2"];
    var col2Keys = ["Value3", "Value4"];
    var col3Keys = ["Value5", "Value6"];
    var col4Keys = ["Value7", "Value8"];


    var Data_styling = {
        Col1: {Data: "No",Datatype: 'None',Data_min: Infinity,Data_max: -Infinity, Data_coloring: 'All', Data_color1: '#000000', Data_color2: '#000000', data_color_complexity: 'One'},
        Col2: {Data: "No",Datatype: 'None',Data_min: Infinity,Data_max: -Infinity, Data_coloring: 'All', Data_color1: '#000000', Data_color2: '#000000', data_color_complexity: 'One'},
        Col3: {Data: "No",Datatype: 'None',Data_min: Infinity,Data_max: -Infinity, Data_coloring: 'All', Data_color1: '#000000', Data_color2: '#000000', data_color_complexity: 'One'},
        Col4: {Data: "No",Datatype: 'None',Data_min: Infinity,Data_max: -Infinity, Data_coloring: 'All', Data_color1: '#000000', Data_color2: '#000000', data_color_complexity: 'One'}
        };

    // Iterate through initialData and update Data_styling accordingly
    Object.keys(list_data_wow).forEach(function(key) {
        var currentObj = list_data_wow[key];

        // Check for Col1 if not already set to "Yes"
        if (Data_styling.Col1.Data !== "Yes" && checkKeysPresent(currentObj, col1Keys)) {
            Data_styling.Col1.Data = "Yes";
        }

        // Check for Col2 if not already set to "Yes"
        if (Data_styling.Col2.Data !== "Yes" && checkKeysPresent(currentObj, col2Keys)) {
            Data_styling.Col2.Data = "Yes";
        }

        // Check for Col3 if not already set to "Yes"
        if (Data_styling.Col3.Data !== "Yes" && checkKeysPresent(currentObj, col3Keys)) {
            Data_styling.Col3.Data = "Yes";
        }

        // Check for Col4 if not already set to "Yes"
        if (Data_styling.Col4.Data !== "Yes" && checkKeysPresent(currentObj, col4Keys)) {
            Data_styling.Col4.Data = "Yes";
        }

        // If both Col1 and Col2 are "Yes", no need to continue iterating
        if (Data_styling.Col1.Data === "Yes" && Data_styling.Col2.Data === "Yes" && Data_styling.Col3.Data === "Yes" && Data_styling.Col4.Data === "Yes") {
            return; // Exit forEach loop early
        }
    });


    // Check for Col1 has data
    // Update Datatype for each column based on data_types
    ['Col1', 'Col2', 'Col3', 'Col4'].forEach(col => {
        if (Data_styling[col].Data === "Yes") {
            if (data_types['Listplot'][col] === 'Discrete') {
                Data_styling[col].Datatype = 'Discrete';
            } else if (data_types['Listplot'][col] === 'Continuous') {
                Data_styling[col].Datatype = 'Continuous';
            }
        }
    });


    // Run through data and set min and max values

    function updateDataMinMax(column, key, currentObj) {
        if (Data_styling[column].Data === "Yes" && Data_styling[column].Datatype === 'Continuous' && currentObj.hasOwnProperty(key)) {
            if (currentObj[key] > Data_styling[column].Data_max) {
                Data_styling[column].Data_max = currentObj[key];
            }
            if (currentObj[key] < Data_styling[column].Data_min) {
                Data_styling[column].Data_min = currentObj[key];
            }
        }
    }

    Object.keys(list_data_wow).forEach(function(key) {
        var currentObj = list_data_wow[key];

        updateDataMinMax('Col1', 'Value2', currentObj);
        updateDataMinMax('Col2', 'Value4', currentObj);
        updateDataMinMax('Col3', 'Value6', currentObj);
        updateDataMinMax('Col4', 'Value8', currentObj);
    });

    return Data_styling;
}


// ##############################
// ### Visualization Function ###
// ##############################
function renderDataVisualization(data, location,styling_option,Layout_dict,data_styling,label_conversion_dict,label_names) {

    // ######################
    // ## Initializzation  ##
    // ######################

    // ## Global Variables ##

    let columns = Layout_dict.columns;
    let col_breaker = Layout_dict.col_breaker;
    let col_spacing = Layout_dict.col_spacing;

    // Calculate total count
    const total_count = calculateTotalCount(data);
    const width = 420*columns;
    const height = 200;
    const margin = { top: 40, right: 20, bottom: 20, left: 20 };

    // X & Y coordinates
    let yOffset = margin.top+5;
    let yOffset_max = 0
    let xOffset = 0;

    // Create svg element
    const svg = d3.select("#" + location)
        .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .attr("id", "ListPlot_visualization");


    // ###############
    // ## Functions ##
    // ###############


    // ## Function to calculate the total count of keys and values ##
    function calculateTotalCount(obj) {
        let totalCount = 0;

        // Recursive function to traverse the nested object
        function traverse(obj) {
            for (const key in obj) {
                // Increment count for each key
                totalCount++;
                if (typeof obj[key] === 'object') {
                    // Recursively traverse nested objects
                    traverse(obj[key]);
                } else if (Array.isArray(obj[key])) {
                    // Increment count for each array item
                    totalCount += obj[key].length;
                }
            }
        }

        traverse(obj);
        return totalCount;
    }

    // ## Shape Funciton ##
    // Function to add different shapes
    function addShape(shapeType, x, y, size, fillColor) {
        switch (shapeType) {
            case 'circle':
                svg.append('circle')
                    .attr('cx', x)
                    .attr('cy', y)
                    .attr('r', size)
                    .style('stroke', 'black')
                    .style('stroke-width', 1)
                    .style('fill', fillColor);
                break;
            case 'rect':
                svg.append('rect')
                    .attr('x', x - size)
                    .attr('y', y - size)
                    .attr('width', size * 2)
                    .attr('height', size * 2)
                    .style('stroke', 'black')
                    .style('stroke-width', 1)
                    .style('fill', fillColor);
                break;
            case 'triangle':
                svg.append('path')
                    .attr('d', `M ${x} ${y - size} L ${x - size} ${y + size} L ${x + size} ${y + size} Z`)
                    .style('stroke', 'black')
                    .style('stroke-width', 1)
                    .style('fill', fillColor);
                break;
            case 'diamond':
                svg.append('path')
                    .attr('d', `M ${x} ${y - size} L ${x - size} ${y} L ${x} ${y + size} L ${x + size} ${y} Z`)
                    .style('stroke', 'black')
                    .style('stroke-width', 1)
                    .style('fill', fillColor);
                break;
            case 'star':
                const points = 5;
                const outerRadius = size;
                const innerRadius = size / 2.5;
                let starPath = '';
                for (let i = 0; i < points * 2; i++) {
                    const angle = (Math.PI / points) * i;
                    const radius = i % 2 === 0 ? outerRadius : innerRadius;
                    const xPoint = x + Math.cos(angle) * radius;
                    const yPoint = y + Math.sin(angle) * radius;
                    starPath += i === 0 ? `M ${xPoint} ${yPoint}` : `L ${xPoint} ${yPoint}`;
                }
                starPath += 'Z';
                svg.append('path')
                    .attr('d', starPath)
                    .style('stroke', 'black')
                    .style('stroke-width', 1)
                    .style('fill', fillColor);
                break;
            default:
                console.log('Unknown shape type');
        }
    }

    // ## Column breaker function ##
    function ColumnsBreakerFunc (columns,label_counter,Col_break_number,state) {
        if (columns == 2 && label_counter > Col_break_number && state.current_col == 1){
            state.current_col = 2;
            xOffset = col_spacing;
            if (yOffset > yOffset_max) {yOffset_max = yOffset;}
            yOffset = margin.top+5;
        }
        if (columns == 3 && label_counter > Col_break_number && state.current_col == 1){
            state.current_col = 2;
            xOffset = col_spacing;
            if (yOffset > yOffset_max) {yOffset_max = yOffset;}
            yOffset = margin.top+5;
        } else if (columns == 3 && label_counter > Col_break_number*2 && state.current_col == 2){
            state.current_col = 3;
            xOffset = col_spacing*2;
            if (yOffset > yOffset_max) {yOffset_max = yOffset;}
            yOffset = margin.top+5;
        }
        if (columns == 4 && label_counter > Col_break_number && state.current_col == 1){
            state.current_col = 2;
            xOffset = col_spacing;
            if (yOffset > yOffset_max) {yOffset_max = yOffset;}
            yOffset = margin.top+5;
        } else if (columns == 4 && label_counter > Col_break_number*2 && state.current_col == 2){
            state.current_col = 3;
            xOffset = col_spacing*2;
            if (yOffset > yOffset_max) {yOffset_max = yOffset;}
            yOffset = margin.top+5;
        } else if (columns == 4 && label_counter > Col_break_number*3 && state.current_col == 3){
            state.current_col = 4;
            xOffset = col_spacing*3;
            if (yOffset > yOffset_max) {yOffset_max = yOffset;}
            yOffset = margin.top+5;
        }
    }

    // #####################
    // ## List processing ##
    // #####################

    function processList(list) {
        // Check the visibility of each layer
        const Layer1_isChecked = d3.select(`#toggle-layer-1`).property('checked');
        const Layer2_isChecked = d3.select(`#toggle-layer-2`).property('checked');
        const Layer3_isChecked = d3.select(`#toggle-layer-3`).property('checked');
        const Layer4_isChecked = d3.select(`#toggle-layer-4`).property('checked');

        let Checklist = [Layer1_isChecked, Layer2_isChecked, Layer3_isChecked, Layer4_isChecked];


        let col_max_label = Layout_dict.col_max_label;
        var Col_break_number = 20;

        if (col_max_label == "Auto") {
            Col_break_number = Math.ceil(total_count / columns);
        } else if (col_max_label == "Custom") {
            Col_break_number = Layout_dict.Col_break_number
        }

        var state = { current_col: 1 };
        // Setup amount of cols
        let label_counter = 1;

        // Set color gradient based on input

        // Define Color_dict to store color scales for each column and shape
        var Color_dict = {};

        // # Setup color styling if "All" is chosen #
        Object.keys(data_styling).forEach(function(col) {
            Color_dict[col] = {};
            if (data_styling[col].data_color_complexity == 'One') {
                Color_dict[col]['All'] = d3.scale.linear().domain([data_styling[col].Data_min,data_styling[col].Data_max]).range(["#FFFFFF", data_styling[col].Data_color1]);
            } else if (data_styling[col].data_color_complexity == 'Two') {
                Color_dict[col]['All'] = d3.scale.linear()
                .domain([data_styling[col].Data_min, (data_styling[col].Data_min + data_styling[col].Data_max) / 2, data_styling[col].Data_max])
                .range([data_styling[col].Data_color1, "#FFFFFF", data_styling[col].Data_color2]);
            }
        });

        // ##############################
        // ## Handle all list elements ##
        // ##############################

        Object.entries(list).forEach(([key, value]) => {
            // Text input function
            function add_text(label, yOffset, layer, bold = false, italic = false, underline = false, fontSize = '16px', color = 'black') {

                label_key = label;

                // # Different types of label handling #
                if (label_names == 'IUPHAR') {

                    const htmlEntities = {
                        '&alpha;': 'α',
                        '&beta;': 'β',
                        '&gamma;': 'γ',
                        '&delta;': 'δ',
                        '&epsilon;': 'ε',
                        '&zeta;': 'ζ',
                        '&eta;': 'η',
                        '&theta;': 'θ',
                        '&iota;': 'ι',
                        '&kappa;': 'κ',
                        '&lambda;': 'λ',
                        '&mu;': 'μ',
                        '&nu;': 'ν',
                        '&xi;': 'ξ',
                        '&omicron;': 'ο',
                        '&pi;': 'π',
                        '&rho;': 'ρ',
                        '&sigma;': 'σ',
                        '&tau;': 'τ',
                        '&upsilon;': 'υ',
                        '&phi;': 'φ',
                        '&chi;': 'χ',
                        '&psi;': 'ψ',
                        '&omega;': 'ω'
                    };

                    // Replace HTML entities with corresponding Unicode characters
                    for (const [entity, char] of Object.entries(htmlEntities)) {
                        label = label.replace(new RegExp(entity, 'g'), char);
                    }

                    const textElement = svg.append('text')
                        .attr('x', margin.left + xOffset + 80)
                        .attr('y', yOffset)
                        .attr('class', layer)
                        .style('dominant-baseline', 'middle') // Set vertical alignment
                        .style('font-weight', bold ? 'bold' : 'normal')
                        .style('font-style', italic ? 'italic' : 'normal')
                        .style('text-decoration', underline ? 'underline' : 'none')
                        .style('font-size', fontSize)
                        .style('fill', color); // Set text color

                    // Calculate the subscript font size (e.g., 75% of the main font size)
                    const mainFontSize = parseFloat(fontSize);
                    const subFontSize = mainFontSize * 0.75 + 'px';

                    // Remove <i> and </i> tags from the label
                    label = label.replace(/<\/?i>/g, '');
                    if (layer == 'Layer-4') {
                        label = label.replace("-adrenoceptor", '').replace(" receptor", '');
                    }
                    const parts = label.split(/(<sub>|<\/sub>)/); // Split label into parts including <sub> tags
                    let inSub = false; // Flag to track if we are inside a subscript

                    // Handle subscripts
                    parts.forEach(part => {
                        if (part === '<sub>') {
                            inSub = true;
                        } else if (part === '</sub>') {
                            inSub = false;
                        } else {
                            const tspan = textElement.append('tspan')
                                .text(part);

                            if (inSub) {
                                tspan.attr('dy', '0.3em') // Subscript positioning
                                    .attr('font-size', subFontSize);
                            } else {
                                tspan.attr('dy', '-0.3em') // Reset positioning
                                    .attr('font-size', fontSize);
                            }
                        }
                    });

                } else if (label_names == 'UniProt') {

                    // # If UniProt, handle receptor names but keep classes etc. #
                    if (layer == 'Layer-4') {
                    label = label_conversion_dict[label_key];
                    label = label.replace(/_human/g, '').toUpperCase();
                    } else {
                        label = label_key
                    }

                    // # Add label element #

                    const textElement = svg.append('text')
                        .attr('x', margin.left + xOffset + 80)
                        .attr('y', yOffset)
                        .attr('class', layer)
                        .attr('dy', '-0.3em') // Adjust this value to move the text higher
                        .style('dominant-baseline', 'middle') // Set vertical alignment
                        .style('font-weight', bold ? 'bold' : 'normal')
                        .style('font-style', italic ? 'italic' : 'normal')
                        .style('text-decoration', underline ? 'underline' : 'none')
                        .style('font-size', fontSize)
                        .style('fill', color)
                        .text(label);


                }

                // #############################
                // ### Implement data points ###
                // #############################

                if (layer == 'Layer-4') {
                    if (list_data_wow.hasOwnProperty(label_key)) {

                        // Initialize all variables for better overview

                        // Col data checker
                        var Col1_data_checker = data_styling.Col1.Data == "Yes";
                        var Col2_data_checker = data_styling.Col2.Data == "Yes";
                        var Col3_data_checker = data_styling.Col3.Data == "Yes";
                        var Col4_data_checker = data_styling.Col4.Data == "Yes";

                        // Col lebels
                        var Col1_shape = list_data_wow[label_key].hasOwnProperty('Value1') ? list_data_wow[label_key].Value1.toLowerCase() : false;
                        var Col2_shape = list_data_wow[label_key].hasOwnProperty('Value3') ? list_data_wow[label_key].Value3.toLowerCase() : false;
                        var Col3_shape = list_data_wow[label_key].hasOwnProperty('Value5') ? list_data_wow[label_key].Value5.toLowerCase() : false;
                        var Col4_shape = list_data_wow[label_key].hasOwnProperty('Value7') ? list_data_wow[label_key].Value7.toLowerCase() : false;

                        // Col data shapes checker
                        var Col1_data = list_data_wow[label_key].hasOwnProperty('Value2') ? list_data_wow[label_key].Value2 : false;
                        var Col2_data = list_data_wow[label_key].hasOwnProperty('Value4') ? list_data_wow[label_key].Value4 : false;
                        var Col3_data = list_data_wow[label_key].hasOwnProperty('Value6') ? list_data_wow[label_key].Value6 : false;
                        var Col4_data = list_data_wow[label_key].hasOwnProperty('Value8') ? list_data_wow[label_key].Value8 : false;

                        // Col Datatype
                        var Col1_datatype = data_styling.Col1.Datatype;
                        var Col2_datatype = data_styling.Col2.Datatype;
                        var Col3_datatype = data_styling.Col3.Datatype;
                        var Col4_datatype = data_styling.Col4.Datatype;

                        // Col data coloring scheme

                        const Col1_coloring = data_styling.Col1.Data_coloring;
                        const Col2_coloring = data_styling.Col2.Data_coloring;
                        const Col3_coloring = data_styling.Col3.Data_coloring;
                        const Col4_coloring = data_styling.Col4.Data_coloring;

                        const Color_list = ['black', 'red', 'blue', 'green'];
                        const Shape_list = ['circle','rect','triangle','star','diamond']

                        const col1_XoffSet = 0
                        const col2_XoffSet = 20
                        const col3_XoffSet = 40
                        const col4_XoffSet = 60

                        // ## Col 1 ##
                        if (Col1_data_checker && (Col1_shape || Col1_data)) {
                            if (Col1_shape) {
                                if (Col1_data && Col1_datatype == 'Discrete') {
                                    shape_color = Color_list.includes(Col1_data.toLowerCase()) ? Col1_data : 'black';
                                } else if (Col1_data && Col1_datatype == 'Continuous') {
                                    if (Col1_coloring == 'All') {
                                        shape_color = Color_dict.Col1.All(Col1_data);
                                    } else if (Col1_coloring == 'Shapes') {
                                        shape_color = Color_dict['Col1'][Col1_shape](Col1_data);
                                    }
                                } else {
                                    shape_color = 'black';
                                }
                                shape_value = Col1_shape === "square" ? 'rect' : Col1_shape;
                                addShape(Shape_list.includes(shape_value) ? shape_value : 'circle', margin.left + xOffset + col1_XoffSet, yOffset - (parseFloat(fontSize) * 0.50), 6, shape_color);
                            } else {
                                if (Col1_data && Col1_datatype == 'Discrete') {
                                    shape_color = Color_list.includes(Col1_data.toLowerCase()) ? Col1_data : 'black';
                                    addShape('circle', margin.left + xOffset + col1_XoffSet, yOffset - (parseFloat(fontSize) * 0.50), 6, shape_color);
                                } else if (Col1_data && Col1_datatype == 'Continuous') {
                                    shape_color = Color_dict.Col1.All(Col1_data);
                                    addShape('circle', margin.left + xOffset + col1_XoffSet, yOffset - (parseFloat(fontSize) * 0.50), 6, shape_color);
                                }
                            }
                        }

                        // ## Col 2 ##
                        if (Col2_data_checker && (Col2_shape || Col2_data)) {
                            if (Col2_shape) {


                                if (Col2_data && Col2_datatype == 'Discrete') {
                                    shape_color = Color_list.includes(Col2_data.toLowerCase()) ? Col2_data : 'black';
                                } else if (Col2_data && Col2_datatype == 'Continuous') {

                                    if (Col2_coloring == 'All') {
                                        shape_color = Color_dict.Col2.All(Col2_data);
                                    } else if (Col2_coloring == 'Shapes') {
                                        shape_color = Color_dict['Col2'][Col2_shape](Col2_data);
                                    }
                                } else {
                                    shape_color = 'black';
                                }
                                shape_value = Col2_shape === "square" ? 'rect' : Col2_shape;
                                addShape(Shape_list.includes(shape_value) ? shape_value : 'circle', margin.left + xOffset + col2_XoffSet, yOffset - (parseFloat(fontSize) * 0.50), 6, shape_color);
                            } else {
                                if (Col2_data && Col2_datatype == 'Discrete') {
                                    shape_color = Color_list.includes(Col2_data.toLowerCase()) ? Col2_data : 'black';
                                    addShape('circle', margin.left + xOffset + col2_XoffSet, yOffset - (parseFloat(fontSize) * 0.50), 6, shape_color);
                                } else if (Col2_data && Col2_datatype == 'Continuous') {
                                    shape_color = Color_dict.Col2.All(Col2_data);
                                    addShape('circle', margin.left + xOffset + col2_XoffSet, yOffset - (parseFloat(fontSize) * 0.50), 6, shape_color);
                                }
                            }
                        }

                        // ## Col 3 ##
                        if (Col3_data_checker && (Col3_shape || Col3_data)) {
                            if (Col3_shape) {
                                if (Col3_data && Col3_datatype == 'Discrete') {
                                    shape_color = Color_list.includes(Col3_data.toLowerCase()) ? Col3_data : 'black';
                                } else if (Col3_data && Col3_datatype == 'Continuous') {
                                    if (Col3_coloring == 'All') {
                                        shape_color = Color_dict.Col3.All(Col3_data);
                                    } else if (Col3_coloring == 'Shapes') {
                                        shape_color = Color_dict['Col3'][Col3_shape](Col3_data);
                                    }
                                } else {
                                    shape_color = 'black';
                                }
                                shape_value = Col3_shape === "square" ? 'rect' : Col3_shape;
                                addShape(Shape_list.includes(shape_value) ? shape_value : 'circle', margin.left + xOffset + col3_XoffSet, yOffset - (parseFloat(fontSize) * 0.50), 6, shape_color);
                            } else {
                                if (Col3_data && Col3_datatype == 'Discrete') {
                                    shape_color = Color_list.includes(Col3_data.toLowerCase()) ? Col3_data : 'black';
                                    addShape('circle', margin.left + xOffset + col3_XoffSet, yOffset - (parseFloat(fontSize) * 0.50), 6, shape_color);
                                } else if (Col3_data && Col3_datatype == 'Continuous') {
                                    shape_color = Color_dict.Col3.All(Col3_data);
                                    addShape('circle', margin.left + xOffset + col3_XoffSet, yOffset - (parseFloat(fontSize) * 0.50), 6, shape_color);
                                }
                            }
                        }

                        // ## Col 4 ##
                        if (Col4_data_checker && (Col4_shape || Col4_data)) {
                            if (Col4_shape) {
                                if (Col4_data && Col4_datatype == 'Discrete') {
                                    shape_color = Color_list.includes(Col4_data.toLowerCase()) ? Col4_data : 'black';
                                } else if (Col4_data && Col4_datatype == 'Continuous') {
                                    if (Col4_coloring == 'All') {
                                        shape_color = Color_dict.Col4.All(Col4_data);
                                    } else if (Col4_coloring == 'Shapes') {
                                        shape_color = Color_dict['Col4'][Col4_shape](Col4_data);
                                    }
                                } else {
                                    shape_color = 'black';
                                }
                                shape_value = Col4_shape === "square" ? 'rect' : Col4_shape;
                                addShape(Shape_list.includes(shape_value) ? shape_value : 'circle', margin.left + xOffset + col4_XoffSet, yOffset - (parseFloat(fontSize) * 0.50), 6, shape_color);
                            } else {
                                if (Col4_data && Col4_datatype == 'Discrete') {
                                    shape_color = Color_list.includes(Col4_data.toLowerCase()) ? Col4_data : 'black';
                                    addShape('circle', margin.left + xOffset + col4_XoffSet, yOffset - (parseFloat(fontSize) * 0.50), 6, shape_color);
                                } else if (Col4_data && Col4_datatype == 'Continuous') {
                                    shape_color = Color_dict.Col4.All(Col4_data);
                                    addShape('circle', margin.left + xOffset + col4_XoffSet, yOffset - (parseFloat(fontSize) * 0.50), 6, shape_color);
                                }
                            }
                        }
                    } else {
                        console.log(label_key);
                    }
                }
            }

            // # Check and append text for Layer 1 breaker if it's visible #

            if (Checklist[0]) {
                if (col_breaker == "Layer1" || col_breaker == "Custom") {
                    ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                }
                add_text(key, yOffset, "Layer-1",styling_option.Layer1.Bold,styling_option.Layer1.Italic,styling_option.Layer1.Underline,styling_option.Layer1.Fontsize,styling_option.Layer1.Color);
                yOffset += 30;
                label_counter++;

                // Check and append text for Layer 2 if it's visible and has sublayers
                if (Checklist[1] && typeof value === 'object') {
                    Object.entries(value).forEach(([subKey, subValue]) => {
                        if (col_breaker == "Layer2" || col_breaker == "Custom") {
                            ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                        }
                        add_text(subKey, yOffset, "Layer-2",styling_option.Layer2.Bold,styling_option.Layer2.Italic,styling_option.Layer2.Underline,styling_option.Layer2.Fontsize,styling_option.Layer2.Color);
                        yOffset += 30;
                        label_counter++;

                        // Check and append text for Layer 3 if it's visible and has sublayers
                        if (Checklist[2] && typeof subValue === 'object') {
                            Object.entries(subValue).forEach(([subSubKey, subSubValue]) => {
                                if (col_breaker == "Layer3" || col_breaker == "Custom") {
                                    ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                                }
                                add_text(subSubKey, yOffset, "Layer-3",styling_option.Layer3.Bold,styling_option.Layer3.Italic,styling_option.Layer3.Underline,styling_option.Layer3.Fontsize,styling_option.Layer3.Color);
                                yOffset += 30;
                                label_counter++;

                                // Check and append text for Layer 4 if it's visible and is an array
                                if (Checklist[3] && Array.isArray(subSubValue)) {
                                    subSubValue.forEach(item => {
                                        if (col_breaker == "Layer4" || col_breaker == "Custom") {
                                            ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                                        }
                                        add_text(item, yOffset, "Layer-4",styling_option.Layer4.Bold,styling_option.Layer4.Italic,styling_option.Layer4.Underline,styling_option.Layer4.Fontsize,styling_option.Layer4.Color);

                                        yOffset += 30;
                                        label_counter++;

                                    });
                                }
                            });
                        } else {
                            Object.entries(subValue).forEach(([subSubKey, subSubValue]) => {
                                // Check and append text for Layer 4 if it's visible and is an array
                                if (Checklist[3] && Array.isArray(subSubValue)) {
                                    subSubValue.forEach(item => {
                                        if (col_breaker == "Layer4" || col_breaker == "Custom") {
                                            ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                                        }
                                        add_text(item, yOffset, "Layer-4",styling_option.Layer4.Bold,styling_option.Layer4.Italic,styling_option.Layer4.Underline,styling_option.Layer4.Fontsize,styling_option.Layer4.Color);
                                        yOffset += 30;
                                        label_counter++;

                                    });
                                }
                            });
                        }
                    });
                } else {
                    Object.entries(value).forEach(([subKey, subValue]) => {
                        // Check and append text for Layer 3 if it's visible and has sublayers
                        if (Checklist[2] && typeof subValue === 'object') {
                            Object.entries(subValue).forEach(([subSubKey, subSubValue]) => {
                                if (col_breaker == "Layer3" || col_breaker == "Custom") {
                                    ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                                }
                                add_text(subSubKey, yOffset, "Layer-3",styling_option.Layer3.Bold,styling_option.Layer3.Italic,styling_option.Layer3.Underline,styling_option.Layer3.Fontsize,styling_option.Layer3.Color);
                                yOffset += 30;
                                label_counter++;

                                // Check and append text for Layer 4 if it's visible and is an array
                                if (Checklist[3] && Array.isArray(subSubValue)) {
                                    subSubValue.forEach(item => {
                                        if (col_breaker == "Layer4" || col_breaker == "Custom") {
                                            ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                                        }
                                        add_text(item, yOffset, "Layer-4",styling_option.Layer4.Bold,styling_option.Layer4.Italic,styling_option.Layer4.Underline,styling_option.Layer4.Fontsize,styling_option.Layer4.Color);
                                        yOffset += 30;
                                        label_counter++;

                                    });
                                }
                            });
                        } else {
                            Object.entries(subValue).forEach(([subSubKey, subSubValue]) => {
                                // Check and append text for Layer 4 if it's visible and is an array
                                if (Checklist[3] && Array.isArray(subSubValue)) {
                                    subSubValue.forEach(item => {
                                        if (col_breaker == "Layer4" || col_breaker == "Custom") {
                                            ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                                        }
                                        add_text(item, yOffset, "Layer-4",styling_option.Layer4.Bold,styling_option.Layer4.Italic,styling_option.Layer4.Underline,styling_option.Layer4.Fontsize,styling_option.Layer4.Color);
                                        yOffset += 30;
                                        label_counter++;

                                    });
                                }
                            });
                        }
                    });
                }
            } else if (!Checklist[0] && Checklist[1]) {
                // Check and append text for Layer 2 if it's visible and has sublayers
                if (Checklist[1] && typeof value === 'object') {
                    Object.entries(value).forEach(([subKey, subValue]) => {
                        if (col_breaker == "Layer2" || col_breaker == "Custom") {
                            ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                        }
                        add_text(subKey, yOffset, "Layer-2",styling_option.Layer2.Bold,styling_option.Layer2.Italic,styling_option.Layer2.Underline,styling_option.Layer2.Fontsize,styling_option.Layer2.Color);
                        yOffset += 30;
                        label_counter++;

                        // Check and append text for Layer 3 if it's visible and has sublayers
                        if (Checklist[2] && typeof subValue === 'object') {
                            Object.entries(subValue).forEach(([subSubKey, subSubValue]) => {
                                if (col_breaker == "Layer3" || col_breaker == "Custom") {
                                    ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                                }
                                add_text(subSubKey, yOffset, "Layer-3",styling_option.Layer3.Bold,styling_option.Layer3.Italic,styling_option.Layer3.Underline,styling_option.Layer3.Fontsize,styling_option.Layer3.Color);
                                yOffset += 30;
                                label_counter++;

                                // Check and append text for Layer 4 if it's visible and is an array
                                if (Checklist[3] && Array.isArray(subSubValue)) {
                                    subSubValue.forEach(item => {
                                        if (col_breaker == "Layer4" || col_breaker == "Custom") {
                                            ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                                        }
                                        add_text(item, yOffset, "Layer-4",styling_option.Layer4.Bold,styling_option.Layer4.Italic,styling_option.Layer4.Underline,styling_option.Layer4.Fontsize,styling_option.Layer4.Color);
                                        yOffset += 30;
                                        label_counter++;

                                    });
                                }
                            });
                        } else {
                            Object.entries(subValue).forEach(([subSubKey, subSubValue]) => {
                                // Check and append text for Layer 4 if it's visible and is an array
                                if (Checklist[3] && Array.isArray(subSubValue)) {
                                    subSubValue.forEach(item => {
                                        if (col_breaker == "Layer4" || col_breaker == "Custom") {
                                            ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                                        }
                                        add_text(item, yOffset, "Layer-4",styling_option.Layer4.Bold,styling_option.Layer4.Italic,styling_option.Layer4.Underline,styling_option.Layer4.Fontsize,styling_option.Layer4.Color);
                                        yOffset += 30;
                                        label_counter++;

                                    });
                                }
                            });
                        }
                    });
                }
            } else if (!Checklist[0] && !Checklist[1] && Checklist[2]) {
                Object.entries(value).forEach(([subKey, subValue]) => {
                    if (Checklist[2] && typeof subValue === 'object') {
                        Object.entries(subValue).forEach(([subSubKey, subSubValue]) => {
                            if (col_breaker == "Layer3" || col_breaker == "Custom") {
                                ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                            }
                            add_text(subSubKey, yOffset, "Layer-3",styling_option.Layer3.Bold,styling_option.Layer3.Italic,styling_option.Layer3.Underline,styling_option.Layer3.Fontsize,styling_option.Layer3.Color);
                            yOffset += 30;
                            label_counter++;

                            // Check and append text for Layer 4 if it's visible and is an array
                            if (Checklist[3] && Array.isArray(subSubValue)) {
                                subSubValue.forEach(item => {
                                    if (col_breaker == "Layer4" || col_breaker == "Custom") {
                                        ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                                    }
                                    add_text(item, yOffset, "Layer-4",styling_option.Layer4.Bold,styling_option.Layer4.Italic,styling_option.Layer4.Underline,styling_option.Layer4.Fontsize,styling_option.Layer4.Color);
                                    yOffset += 30;
                                    label_counter++;

                                });
                            }
                        });
                    }
                });
            } else if (!Checklist[0] && !Checklist[1] && !Checklist[2] && Checklist[3]) {
                Object.entries(value).forEach(([subKey, subValue]) => {
                    Object.entries(subValue).forEach(([subSubKey, subSubValue]) => {
                        if (Checklist[3] && Array.isArray(subSubValue)) {
                            subSubValue.forEach(item => {
                                if (col_breaker == "Layer4" || col_breaker == "Custom") {
                                    ColumnsBreakerFunc(columns,label_counter,Col_break_number,state);
                                }
                                add_text(item, yOffset, "Layer-4",styling_option.Layer4.Bold,styling_option.Layer4.Italic,styling_option.Layer4.Underline,styling_option.Layer4.Fontsize,styling_option.Layer4.Color);
                                yOffset += 30;
                                label_counter++;

                            });
                        }
                    });
                });
            }
        });

    // # check if y-off-set is more than max, if so expand canvas #
    if (yOffset_max < yOffset) {
        yOffset_max = yOffset;
    }

    // Rerender height of plot as the last thing
    svg.attr("height",yOffset_max + margin.top + margin.bottom)
}
    // # Run list process and return svg #
    processList(data,styling_option);
    return svg;
}

// RENDER THE CLUSTER

function renderClusterPlot(data, selector, method) {
    const width = 800;
    const height = 600;
    const margin = { top: 80, right: 150, bottom: 50, left: 50 };

    // Clear the previous plot
    d3.select(selector).html("");

    const svg = d3.select(selector)
        .append("svg")
        .attr("width", width)
        .attr("height", height)
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    const x = d3.scale.linear()
        .domain(d3.extent(data, d => +d.x))
        .range([0, width - margin.left - margin.right]);

    const y = d3.scale.linear()
        .domain(d3.extent(data, d => +d.y))
        .range([height - margin.top - margin.bottom, 0]);

    const color = d3.scale.category10();

    let selectedCluster = null;

    const circles = svg.selectAll(".dot")
        .data(data)
        .enter().append("circle")
        .attr("class", "dot")
        .attr("cx", d => x(+d.x))
        .attr("cy", d => y(+d.y))
        .attr("r", 5)
        .attr("cluster", d => d.cluster)
        .style("fill", d => color(d.cluster))
        .attr("id", d => "circle" + d.label.replace(/\[|\]|\(|\)|\s|\,|\'/g,""))
        .on("click", handleClick);

    const labels = svg.selectAll(".label")
        .data(data)
        .enter().append("text")
        .attr("class", "label")
        .attr("x", d => x(+d.x))
        .attr("y", d => y(+d.y))
        .attr("dy", -10)
        .attr("cluster", d => d.cluster)
        .text(d => d.label)
        .style("font-size", "11px")
        .style("text-anchor", "middle")
        .on("click", function(d, i) {
            d3.event.stopPropagation();
            handleClick.call(circles[0][i], d);
        });

    // Add leader lines from circles to labels
    const lines = svg.selectAll(".line")
        .data(data)
        .enter().append("line")
        .attr("class", "line")
        .style("stroke", "black")
        .style("stroke-width", 1);

    // Add title based on the method with more spacing
    svg.append("text")
        .attr("x", (width - margin.left - margin.right) / 2)
        .attr("y", -50)
        .attr("text-anchor", "middle")
        .style("font-size", "20px")
        .text("Cluster Plot (" + method.toUpperCase() + ")");

    // Add force simulation to spread out circles
    const simulation = d3v4.forceSimulation(data)
        .force("x", d3v4.forceX(d => x(d.x)).strength(1))
        .force("y", d3v4.forceY(d => y(d.y)).strength(1))
        .force("collide", d3.forceCollide(10)) // Set to twice the circle radius for spacing
        .on("tick", ticked);

    function ticked() {
        circles
            .attr("cx", d => d.x = Math.max(5, Math.min(width - 5, d.x))) // Prevent circles from going outside SVG
            .attr("cy", d => d.y = Math.max(5, Math.min(height - 5, d.y))); // Prevent circles from going outside SVG

        labels
            .attr("x", d => d.x)
            .attr("y", d => d.y - 10);

        lines
            .attr("x1", d => x(d.x))
            .attr("y1", d => y(d.y))
            .attr("x2", d => d.x)
            .attr("y2", d => d.y - 10);
    }

    // Function to handle click event
    function handleClick(d) {
        d3.event.stopPropagation();
        const cluster = d.cluster;

        if (selectedCluster === cluster) {
            selectedCluster = null;
            resetOpacity();
        } else {
            resetOpacity();
            selectedCluster = cluster;
            circles.style("opacity", 0.2);
            labels.style("opacity", 0.2);
            lines.style("opacity", 0.2);

            circles.filter(c => c.cluster === cluster)
                .style("opacity", 1)
                .style("stroke", "black")
                .style("stroke-width", 2);

            labels.filter(c => c.cluster === cluster)
                .style("opacity", 1);

            lines.filter(c => c.cluster === cluster)
                .style("opacity", 1);
        }
    }

    // Function to reset the opacity when clicking outside circles
    function resetOpacity() {
        circles.style("opacity", 1).style("stroke", "none").style("stroke-width", 0);
        labels.style("opacity", 1).style("font-weight", "normal");
        lines.style("opacity", 1);
    }

    // Add an overlay to capture clicks outside circles
    svg.append("rect")
        .attr("width", width - margin.left - margin.right)
        .attr("height", height - margin.top - margin.bottom)
        .style("fill", "none")
        .style("pointer-events", "all")
        .on("click", resetOpacity);

    // Create a legend based on clusters
    const uniqueClusters = [...new Set(data.map(d => d.cluster))].sort((a, b) => a - b);

    const legend = svg.selectAll(".legend")
        .data(uniqueClusters)
        .enter().append("g")
        .attr("class", "legend")
        .attr("transform", (d, i) => `translate(${width - margin.right + 10},${i * 20})`)
        .on("click", function(d) {
            selectedCluster = d;
            circles.style("opacity", 0.2);

            circles.filter(c => c.cluster === d)
                .style("opacity", 1)
                .style("stroke", "black")
                .style("stroke-width", 2);

            d3.event.stopPropagation();
        });

    legend.append("circle")
        .attr("cx", 10)
        .attr("cy", 9)
        .attr("r", 5)
        .style("fill", d => color(d));

    legend.append("text")
        .attr("x", 20)
        .attr("y", 14)
        .text(d => `Cluster ${d}`)
        .style("font-size", "12px");

    // Adjust label positions to prevent overlap
    adjustLabelPositions(labels);
}

// Function to adjust label positions to prevent overlap
function adjustLabelPositions(labels) {
    const labelPadding = 2; // Padding between labels to avoid overlap
    let iterations = 10; // Maximum iterations to adjust labels

    while (iterations-- > 0) {
        let overlaps = false;

        labels.each(function() {
            const label = d3.select(this);
            const bbox = label.node().getBBox();

            labels.each(function() {
                const otherLabel = d3.select(this);
                if (label.node() !== otherLabel.node()) {
                    const otherBBox = otherLabel.node().getBBox();
                    if (intersect(bbox, otherBBox)) {
                        overlaps = true;
                        const dx = (bbox.x + bbox.width / 2) - (otherBBox.x + otherBBox.width / 2);
                        const dy = (bbox.y + bbox.height / 2) - (otherBBox.y + otherBBox.height / 2);

                        const offset = labelPadding / Math.sqrt(dx * dx + dy * dy);
                        const moveX = dx * offset;
                        const moveY = dy * offset;

                        label.attr("x", parseFloat(label.attr("x")) + moveX)
                            .attr("y", parseFloat(label.attr("y")) + moveY);
                        otherLabel.attr("x", parseFloat(otherLabel.attr("x")) - moveX)
                            .attr("y", parseFloat(otherLabel.attr("y")) - moveY);

                        bbox.x += moveX;
                        bbox.y += moveY;
                    }
                }
            });
        });

        if (!overlaps) break; // Stop if no overlaps found
    }
}

// Function to check if two bounding boxes intersect
function intersect(bbox1, bbox2) {
    return !(bbox2.x > bbox1.x + bbox1.width ||
             bbox2.x + bbox2.width < bbox1.x ||
             bbox2.y > bbox1.y + bbox1.height ||
             bbox2.y + bbox2.height < bbox1.y);
}


// function renderClusterPlot(data, selector, method) {
//     const width = 800;
//     const height = 600;
//     const margin = { top: 80, right: 50, bottom: 50, left: 50 };
//
//     // Clear the previous plot
//     d3.select(selector).html("");
//
//     const svg = d3.select(selector)
//         .append("svg")
//         .attr("width", width)
//         .attr("height", height)
//         .append("g")
//         .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
//
//     const x = d3.scale.linear()
//         .domain(d3.extent(data, d => +d.x))
//         .range([0, width - margin.left - margin.right]);
//
//     const y = d3.scale.linear()
//         .domain(d3.extent(data, d => +d.y))
//         .range([height - margin.top - margin.bottom, 0]);
//
//     const color = d3.scale.category10();
//
//     const circles = svg.selectAll(".dot")
//         .data(data)
//         .enter().append("circle")
//         .attr("class", "dot")
//         .attr("cx", d => x(+d.x))
//         .attr("cy", d => y(+d.y))
//         .attr("r", 5)
//         .attr("cluster", d => d.cluster)
//         .style("fill", d => color(d.cluster))
//         .attr("id", d => "circle" + d.label.replace(/\[|\]|\(|\)|\s|\,|\'/g,""))
//         .on("click", handleClick);
//
//     const labels = svg.selectAll(".text")
//         .data(data)
//         .enter().append("text")
//         .attr("x", d => x(+d.x))
//         .attr("y", d => y(+d.y))
//         .attr("dy", -10)
//         .attr("cluster", d => d.cluster)
//         .text(d => d.label)
//         .style("font-size", "11px")
//         .style("text-anchor", "middle");
//         // .on("click", function(d, i) {
//         //     d3.event.stopPropagation();
//         //     handleClick.call(circles[0][i], d);
//         // });
//
//     // Add title based on the method with more spacing
//     svg.append("text")
//         .attr("x", (width - margin.left - margin.right) / 2)
//         .attr("y", -50)
//         .attr("text-anchor", "middle")
//         .style("font-size", "20px")
//         .text("Cluster Plot (" + method.toUpperCase() + ")");
//
//     // Function to handle click event
//     function handleClick(d) {
//         const cluster = d.cluster;
//         circles.style("opacity", 0.2);
//         labels.style("opacity", 0.2);
//
//         circles.filter(c => c.cluster === cluster)
//             .style("opacity", 1)
//             .style("stroke", "black")
//             .style("stroke-width", 2);
//
//         labels.filter(c => c.cluster === cluster)
//             .style("opacity", 1)
//             .style("font-weight", "bold");
//     }
//
//     // Function to reset the opacity when clicking outside circles
//     function resetOpacity() {
//         circles.style("opacity", 1).style("stroke", "none").style("stroke-width", 0);
//         labels.style("opacity", 1).style("font-weight", "normal");
//     }
//
//     // Add an overlay to capture clicks outside circles
//     // svg.append("rect")
//     //     .attr("width", width - margin.left - margin.right)
//     //     .attr("height", height - margin.top - margin.bottom)
//     //     .style("fill", "none")
//     //     .style("pointer-events", "all")
//     //     .on("click", resetOpacity);
// }

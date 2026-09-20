// #################
// ###   TREE    ###
// #################

// Backend trees may ship blank per-node `color` (the palette now lives client-side, not in
// Python). draw_tree()'s branch lines read d.target.color directly for stroke -- an empty string
// clears the style, so branch lines would otherwise vanish. These are only a fallback: real colors
// coming from the backend are never overwritten (see the `if (!n.color)` guards below).
var TREE_CLASS_COLORS = {
    'A': 'Red', 'B1': 'Green', 'B2': 'Blue', 'C': 'Purple', 'F': 'Grey',
    'T2': 'Orange', 'O1': '#66CDAA', 'O2': '#3CB371', 'V': '#B8860B', 'U': 'Gold'
};
// Ligand-type ("chemotype") colors; mirrors CSS_COLORS in common/phylogenetic_tree.py, including
// its Class B2 blanket override and its 'Black' default for anything unrecognized.
var TREE_CHEMOTYPE_COLORS = {
    'Adhesion': 'Crimson', 'Alicarboxylic acid': 'Red', 'Aminergic': 'OrangeRed',
    'Amino acid': 'Orange', 'Ion': 'GoldenRod', 'Lipid': 'Gold', 'Melatonin': 'Yellow',
    'Nucleotide': 'YellowGreen', 'Orphan': 'Gold', 'Other': 'Green', 'Peptide': 'SkyBlue',
    'Protein': 'SteelBlue', 'Sensory': 'Indigo', 'Steroid': 'Purple'
};

function tree_class_key(label) {
    // "Unclassified" is the class's real DB name; its canonical short symbol is "U".
    var raw = String(label || '').trim();
    return raw.startsWith('Unclassified') ? 'U' : raw.split(" (")[0].replace(/^Class\s+/i, '').trim();
}

function tree_hash_color(key) {
    // Stable fallback for a class/chemotype this palette doesn't recognize yet (e.g. a brand new
    // class), so it renders a consistent color instead of silently going colorless.
    var hash = 0;
    var str = String(key || '');
    for (var i = 0; i < str.length; i++) {
        hash = ((hash << 5) - hash) + str.charCodeAt(i);
        hash |= 0;
    }
    return 'hsl(' + (Math.abs(hash) % 360) + ', 68%, 48%)';
}

function tree_class_color(label) {
    var key = tree_class_key(label);
    return TREE_CLASS_COLORS[key] || tree_hash_color(key);
}

function ensure_tree_colors(data, depth) {
    var classNodes = (depth === 4) ? (data.children || []) : [data];
    classNodes.forEach(function (cls) {
        var classKey = tree_class_key(cls.name);
        if (!cls.color) { cls.color = tree_class_color(cls.name); }
        (cls.children || []).forEach(function (lt) {
            var chemKey = tree_class_key(lt.name);
            var chemColor = (classKey === 'B2') ? 'LimeGreen' : (TREE_CHEMOTYPE_COLORS[chemKey] || 'Black');
            (function paint(n) {
                if (!n.color) { n.color = chemColor; }
                (n.children || []).forEach(paint);
            })(lt);
        });
    });
    return data;
}

// Restructure tree data if only one Class is present
function update_tree_data(data,depth) {
    data = ensure_tree_colors(data, depth);
    if (depth === 4) {
        // Iterate over each class and update name
        data.children.forEach(Class_child => {
            Class_child.name = Class_child.name.split(" (")[0];
            // iterate over ligand type
            Class_child.children.forEach(LigandType => {
                // iterate over receptor family
                LigandType.children.forEach(ReceptorFamily => {
                    // Trimming receptor family
                    ReceptorFamily.name = ReceptorFamily.name.replace(/( receptors|neuropeptide )/g, '');
                    ReceptorFamily.name = ReceptorFamily.name.split(" (")[0];
                });
            });

        });
    } else if (depth === 3) {
         // iterate over ligand type
         data.children.forEach(LigandType => {
            // iterate over receptor family
            LigandType.children.forEach(ReceptorFamily => {
                // Trimming receptor family
                ReceptorFamily.name = ReceptorFamily.name.replace(/( receptors|neuropeptide )/g, '');
                ReceptorFamily.name = ReceptorFamily.name.split(" (")[0];
            });
        });
    }
    return data
}

// Draw the Tree
function draw_tree(data, options) {

    // Remove existing SVG if present
    d3.select('#' + options.anchor + "_svg").remove();

    var branches = {};
    var branch_offset = 0;
    var thickness = options.depth + 1;
    for (var key in options.branch_length) {
        if (key === options.depth) { continue };
        if (options.label_free.includes(parseInt(key,10))) {
            branch_offset = branch_offset + 10;
        } else {
            if (options.branch_trunc !== 0) {
                branch_offset = branch_offset + 2 * options.branch_trunc + 10;
            } else {
                branch_offset = branch_offset + string_pixlen(options.branch_length[key], key);
            }
        }
        branches[key] = branch_offset;
    }
    branches[options.depth] = branch_offset + options.leaf_offset;

    var diameter = 2 * branches[options.depth] + 140;

    var tree = d3.layout.tree()
        .size([360, diameter / 2])
        .separation(function (a, b) { return (a.parent === b.parent ? 1 : 2) / Math.max(a.depth, 1); });

    var diagonal = d3.svg.diagonal.radial()
        .projection(function (d) { return [d.y, d.x / 180 * Math.PI]; });

    var svg = d3.select('#' + options.anchor).append("svg")
        .attr("width", diameter)
        .attr("height", diameter)
        .attr("id", options.anchor + "_svg")
        .attr("xmlns", "http://www.w3.org/2000/svg");

    var svg_g = svg.append("g")
        .attr("transform", "translate(" + diameter / 2 + "," + diameter / 2 + ")");

    var nodes = tree.nodes(data);

    nodes.forEach(function (d) {
        if (d.depth === 0) {
            d.y = 0
        } else {
            d.y = branches[d.depth]
        }
    });

    var links = tree.links(nodes);

    var link = svg_g.append("g")
        .attr("class", "links")
        .selectAll("path")
        .data(links)
        .enter().append("path")
        .each(function (d) { d.target.linkNode = this; })
        .attr("d", diagonal) //function (d) { return step(d.source.x, d.source.y, d.target.x, d.target.y) })
        .style("stroke", function (d) { return d.target.color; })
        .style("stroke-width", function (d) { if (d.target.depth > 0) { return thickness - d.target.depth; } else { return 0; } })
        .style("fill-opacity", 0)
        .style("opacity", function (d) {
            if ((d.target.interactions > 0 && d.target.mutations_an > 0) || 1 === 1) { return 0.8 } //|| 1==1
            else if (d.target.interactions > 0) { return 0.5 }
            else if (d.target.mutations_an > 0) { return 0.5 }
            else { return 0.1 };
        });

    var node = svg_g.selectAll(".node")
        .data(nodes)
        .enter().append("g")
        .attr("class", "node")
        .attr("transform", function (d) { if (d.name === '') { return "rotate(" + (d.x) + ")translate(" + d.y + ")"; } else { return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")"; } })

    node.filter(function (d) { return (d.depth === options.depth) })
        .attr("id", function (d) { if (d.name === '') { return "innerNode" } else { return 'X' + d.name.toUpperCase() } });

    node.append("text")
        .attr("dy", ".31em")
        .attr("name", function (d) { if (d.name === '') { return "branch" } else { return d.name } })
        .attr("text-anchor", function (d) {
            if (d.depth === options.depth) {
                return d.x < 180 ? "start" : "end";
            } else {
                return d.x < 180 ? "end" : "start";
            }
        })
        .attr("transform", function (d) {
            var labelOffset = parseFloat(circle_size) + 4;  // Adjust the offset as needed
            if (d.depth === options.depth) {
                return d.x < 180 ? `translate(${labelOffset})` : `rotate(180)translate(-${labelOffset})`;
            } else {
                return d.x < 180 ? "translate(-12)" : "rotate(180)translate(12)";
            }
        })
        .text(function (d) {
            if (d.depth === options.depth) {
                return d.name.toUpperCase();
            } else if (options.label_free.includes(d.depth)) {
                return "";
            } else if (d.depth > 0) {
                return d.name;
            } else {
                return "";
            }
        })
        .call(wrap, options.branch_trunc)
        .style("font-size", function (d) {
            // Use the custom font size from options
            if (options.depth === 4) {
            if (d.depth === 1) { return options.fontSize.class; }
            else if (d.depth === 2) { return options.fontSize.ligandtype; }
            else if (d.depth === 3) { return options.fontSize.receptorfamily; }
            else { return options.fontSize.receptor; }
        } else {
            if (d.depth === 1) { return options.fontSize.ligandtype; }
            else if (d.depth === 2) { return options.fontSize.receptorfamily; }
            else { return options.fontSize.receptor; }
        }
        })
        .style("font-family", "Palatino")
        .style("fill", function (d) {
            if (d.color) { return "#111"; }
            else { return "#222"; };
        }).call(getBB);

    node.filter(function (d) { return (d.depth !== options.depth) }).insert("rect", "text")
        .attr("x", function (d) { return d.x < 180 ? d.bbox.x - 12 : d.bbox.x - d.bbox.width - 12; })
        .attr("y", function (d) { return d.bbox.y; })
        .attr("width", function (d) { return d.bbox.width; })
        .attr("height", function (d) { return d.bbox.height; })
        .style("fill", "#FFF");

    function step(startAngle, startRadius, endAngle, endRadius) {
        var c0 = Math.cos(startAngle = (startAngle - 90) / 180 * Math.PI),
            s0 = Math.sin(startAngle),
            c1 = Math.cos(endAngle = (endAngle - 90) / 180 * Math.PI),
            s1 = Math.sin(endAngle);
        return "M" + startRadius * c0 + "," + startRadius * s0
            + (endAngle === startAngle ? "" : "A" + startRadius + "," + startRadius + " 0 0 " + (endAngle > startAngle ? 1 : 0) + " " + startRadius * c1 + "," + startRadius * s1)
            + "L" + endRadius * c1 + "," + endRadius * s1;
    }

    function string_pixlen(text, depth) {
        var canvas = document.createElement('canvas');
        var ctx = canvas.getContext("2d");

        // Check the depth condition
        if (options.depth === 4) {
            // Use the custom font size from options for all levels
            if (depth === 1) {
                ctx.font = options.fontSize.class + " Palatino";
            } else if (depth === 2) {
                ctx.font = options.fontSize.ligandtype + " Palatino";
            } else if (depth === 3) {
                ctx.font = options.fontSize.receptorfamily + " Palatino";
            } else {
                ctx.font = options.fontSize.receptor + " Palatino";
            }
        } else if (options.depth === 3) {
            // Omit class and start from ligandtype
            if (depth === 1) {
                ctx.font = options.fontSize.ligandtype + " Palatino";
            } else if (depth === 2) {
                ctx.font = options.fontSize.receptorfamily + " Palatino";
            } else {
                ctx.font = options.fontSize.receptor + " Palatino";
            }
        }

        // Measure and return the text width plus padding
        return parseInt(ctx.measureText(text).width,10);
    }

    function getBB(selection) {
        selection.each(function (d) { d.bbox = this.getBBox(); });;
    }

    function wrap(text, width) {
        if (width === 0) {
            return;
        }
        text.each(function () {
            var text = d3.select(this),
                words = text.text().split(/\s+/).reverse(),
                word,
                line = [],
                lineNumber = 0,
                lineHeight = 1.1, // ems
                y = text.attr("y"),
                dy = parseFloat(text.attr("dy")),
                tspan = text.text(null).append("tspan").attr("x", 0).attr("y", y).attr("dy", dy + "em");

            word = words.pop(); // Initialize the first word
            while (word !== undefined) {
                line.push(word);
                tspan.text(line.join(" "));

                if (tspan.node().getComputedTextLength() > width) {
                    line.pop(); // Remove the last word as it exceeds the width
                    tspan.text(line.join(" ")); // Update the tspan with the valid words
                    line = [word]; // Start a new line with the word that didn't fit
                    tspan = text.append("tspan")
                                .attr("x", 0)
                                .attr("y", y)
                                .attr("dy", ++lineNumber * lineHeight + dy + "em")
                                .text(word); // Add the word to the new tspan
                }

                word = words.pop(); // Reassign the next word at the end of the loop
            }
        });
    }
    // === Centering Logic ===
    // Calculate the actual bounding box of the content
    var extraPadding = 40

    var contentWidth = diameter + extraPadding;
    var contentHeight = diameter + extraPadding;

    // Set up SVG with viewBox instead of fixed width/height
    svg.attr('width', '100%')
    .attr('height', '100%')
    .attr('viewBox', `0 0 ${contentWidth} ${contentHeight}`);

    // Center the tree in the available space
    var cx = contentWidth / 2;
    var cy = contentHeight / 2;

    svg.select('g')
    .attr('transform', `translate(${cx},${cy})`);
}

// replaces labels derived from view
function formatTextWithHTML(text) {
    // Apply all the replacements step by step
    return text
        .replace(" receptor", '')
        .replace("-adrenoceptor", '')
        .replace(" receptor-", '-')
        .replace("<sub>", '</tspan><tspan baseline-shift="sub">')
        .replace("</sub>", '</tspan><tspan>')
        .replace("<i>", '</tspan><tspan font-style="italic">')
        .replace("</i>", '</tspan><tspan>')
        .replace("Long-wave-sensitive", 'LWS')
        .replace("Medium-wave-sensitive", 'MWS')
        .replace("Short-wave-sensitive", 'SWS')
        .replace("Olfactory", 'OLF')
        .replace("calcitonin-like receptor", 'CLR');
}

// Change the labels
function changeLeavesLabels(location, value, dict, styling_dict) {
    // Initialize leaf node length
    styling_dict.starter = 0;


    // Select all leaf nodes
    let gNodes = d3.select('#' + location).selectAll('g');

    gNodes.each(function(d) {
        let g = d3.select(this);
        if (g.attr("id") !== null) {
            let name = g.attr("id").substring(1);  // skip leading "X"

            let labelName;
            if (value === "UniProt") {
                labelName = name;
            } else if (dict[name]) {
                labelName = dict[name][0];
                labelName = formatTextWithHTML(labelName);
            } else {
                labelName = name; // fallback if label not found
            }

            let node = d3.select('#' + location).select('#X' + name);
            if (node.size() !== 0) {
                node.selectAll("text")[0].forEach(function(node_label) {
                    node_label.innerHTML = labelName;

                    let labelSize = node_label.getBBox().width;
                    if (labelSize > styling_dict.starter) {
                        styling_dict.starter = labelSize ;
                    }
                });
            }
        }
    });
}

// Draw the circles (data) of the tree plot

function DrawCircles(location, data, Tree_colors, circles_styling, circle_styling_dict = {}) {

    const starter = circles_styling.starter || 1;
    const clean = circles_styling.clean || true;
    const gradient = circles_styling.gradient || true;
    const circle_spacer = circles_styling.circle_spacer || 10;
    const circle_size = circles_styling.circle_size || 5;
    const mode = circles_styling.mode ||"Numeric";

    var svg = d3.select('#' + location);
    var node = svg.selectAll(".node");

    if (clean === true) {
        node.selectAll("circle.outerCircle, circle.innerCircle").remove();
    }

    var minMaxValues = {};

    // Compute min/max values only for Numeric mode
    if (mode === "Numeric") {
        for (var x in data) {
            var keys = Object.keys(data[x]).filter(key =>
                ["Inner", "Outer1", "Outer2", "Outer3", "Outer4", "Outer5"].includes(key)
            );
            for (var unit of keys) {
                if (!minMaxValues[unit]) {
                    minMaxValues[unit] = { min: Infinity, max: -Infinity };
                }
                var value = data[x][unit];
                if (value < minMaxValues[unit].min) minMaxValues[unit].min = value;
                if (value > minMaxValues[unit].max) minMaxValues[unit].max = value;
            }
        }
    }

    for (var x in data) {
        // Filter keys based on mode
        var keys = Object.keys(data[x]).filter(key => {
            if (mode === "Numeric") {
                return ["Inner", "Outer1", "Outer2", "Outer3", "Outer4", "Outer5"].includes(key);
            } else if (mode === "Text") {
                return key === "Inner";
            }
            return false;
        });

        for (var unit of keys) {
            var leaf = svg.selectAll('g[id=X' + x + ']');
            var isInner = unit === 'Inner';
            var className = isInner ? "innerCircle" : "outerCircle";

            var transform = isInner
                ? "translate(0,0)"
                : "translate(" + (Math.ceil(starter) + (Object.keys(Tree_colors).indexOf(unit) + 1) * circle_spacer) + ",0)";

            var fillColor = "#FFFFFF"; // Default color

            if (mode === "Numeric") {
                if (Tree_colors[unit]) {
                    var value = data[x][unit];
                    var minValue = minMaxValues[unit].min;
                    var maxValue = minMaxValues[unit].max;
                    if (!(isFinite(value) && isFinite(minValue) && isFinite(maxValue))) {
                        fillColor = Tree_colors[unit][0];
                    } else {
                        var styling = circle_styling_dict[unit] || "Two";
                        var colorScale;

                        if (minValue === maxValue) {
                            var eps = 1e-9 * (Math.abs(minValue) || 1);
                            minValue = minValue - eps;
                            maxValue = maxValue + eps;
                        }

                        if (styling === "One") {
                            colorScale = d3.scale.linear()
                                .domain([minValue, maxValue])
                                .range(["#FFFFFF", Tree_colors[unit][1]]);
                        } else if (styling === "Three") {
                            colorScale = d3.scale.linear()
                                .domain([minValue, (minValue + maxValue) / 2, maxValue])
                                .range([Tree_colors[unit][0], "#FFFFFF", Tree_colors[unit][1]]);
                        } else {
                            colorScale = d3.scale.linear()
                                .domain([minValue, maxValue])
                                .range(Tree_colors[unit]);
                        }

                        fillColor = gradient ? colorScale(value) : Tree_colors[unit][0];
                        if (!gradient && fillColor === undefined) {
                            fillColor = Tree_colors[unit][0];
                        }
                    }
                }
            } else if (mode === "Text") {
                if (data[x] && "ColorValue" in data[x]) {
                    fillColor = data[x].ColorValue;
                }
            }

            leaf.append("circle")
                .attr("r", circle_size)
                .attr("class", className)
                .style("stroke", "black")
                .style("stroke-width", 0.8)
                .style("fill", fillColor)
                .attr("transform", transform);
        }
    }
    // === Optional viewBox tweak after stacking data circles ===
    // Legacy math mixed scaled group bbox with pixel translate → invalid viewBoxes / SVG NaNs on some trees.
    if (mode === "Numeric" && circles_styling && !circles_styling.skipNumericViewBoxAdjust) {
        const elSVG = svg.select('svg');
        const g = elSVG.select('g');

        if (!g.empty()) {
            const bbox = g.node().getBBox();
            const match = g.attr("transform")?.match(/translate\(([^,]+),([^)]+)\)/);

            if (match && bbox && isFinite(bbox.width) && isFinite(bbox.height)) {
                const tx = parseFloat(match[1]);
                const ty = parseFloat(match[2]);
                const padding = 20;

                const viewBoxX = tx - bbox.width / 2 - padding;
                const viewBoxY = ty - bbox.height / 2 - padding;
                const viewBoxWidth = bbox.width + padding * 2;
                const viewBoxHeight = bbox.height + padding * 2;

                if (isFinite(viewBoxX) && isFinite(viewBoxY) && isFinite(viewBoxWidth) && viewBoxWidth > 0 && isFinite(viewBoxHeight) && viewBoxHeight > 0) {
                    elSVG
                        .attr("viewBox", `${viewBoxX} ${viewBoxY} ${viewBoxWidth} ${viewBoxHeight}`)
                        .attr("preserveAspectRatio", "xMidYMid meet");
                }
            }
        }
    }
}

// Create the bar legends
function createLegendBars(location, data, conversion, circle_styling_dict, datatype_dict, label_dict = {}, TreeLegendPosition = "bottom") {
    const svg = d3.select('#' + location + ' svg');
    svg.selectAll(".legend-group").remove();

    const legendGroup = svg.append("g").attr("class", "legend-group");

    const margin = { top: 50, right: 20, bottom: 30, left: 60 };
    const barWidth = 100;
    const barHeight = 15;
    const spacing = 50;

    const categoryMax = {};

    // Collect min/max values
    Object.values(data).forEach(item => {
        Object.entries(item).forEach(([category, value]) => {
            if (category in conversion) {
                if (!categoryMax[category]) {
                    categoryMax[category] = { min: value, max: value };
                } else {
                    categoryMax[category].min = Math.min(categoryMax[category].min, value);
                    categoryMax[category].max = Math.max(categoryMax[category].max, value);
                }
            }
        });
    });

    const existingCategories = Object.keys(categoryMax);

    const colorScales = {};
    let skippedDiscreteCount = 0;

    // Create color gradients
    existingCategories.forEach(category => {
        if (datatype_dict[category] === "Discrete") return;

        const colors = conversion[category];
        const gradientId = `gradient-${category}`;
        const gradient = legendGroup.append("defs")
            .append("linearGradient")
            .attr("id", gradientId)
            .attr("x1", "0%").attr("x2", "100%")
            .attr("y1", "0%").attr("y2", "0%");

        const styling = circle_styling_dict[category] || "Two";
        if (styling === "Three") {
            gradient.append("stop").attr("offset", "0%").attr("stop-color", colors[0]);
            gradient.append("stop").attr("offset", "50%").attr("stop-color", "#FFFFFF");
            gradient.append("stop").attr("offset", "100%").attr("stop-color", colors[1]);
        } else if (styling === "One") {
            gradient.append("stop").attr("offset", "0%").attr("stop-color", "#FFFFFF");
            gradient.append("stop").attr("offset", "100%").attr("stop-color", colors[1]);
        } else {
            gradient.append("stop").attr("offset", "0%").attr("stop-color", colors[0]);
            gradient.append("stop").attr("offset", "100%").attr("stop-color", colors[1]);
        }

        colorScales[category] = gradientId;
    });

    // Draw bars
    existingCategories.forEach((category, index) => {
        if (datatype_dict[category] !== "Continuous") {
            skippedDiscreteCount++;
            return;
        }

        const BarText = label_dict[category] || category.replace(/([a-zA-Z]+)(\d+)/, '$1 $2');
        const gradientId = colorScales[category];
        const min = parseFloat(categoryMax[category].min.toFixed(2));
        const max = parseFloat(categoryMax[category].max.toFixed(2));
        const x = (index - skippedDiscreteCount) * (barWidth + spacing);

        legendGroup.append("rect")
            .attr("x", x)
            .attr("y", margin.top)
            .attr("width", barWidth)
            .attr("height", barHeight)
            .style("fill", `url(#${gradientId})`)
            .style("stroke", "black")
            .style("stroke-width", "1px");

        legendGroup.append("text")
            .attr("x", x + barWidth / 2)
            .attr("y", margin.top - 10)
            .attr("text-anchor", "middle")
            .text(BarText);

        legendGroup.append("text")
            .attr("x", x)
            .attr("y", margin.top + barHeight + 15)
            .attr("text-anchor", "start")
            .text(`${min}`);

        legendGroup.append("text")
            .attr("x", x + barWidth)
            .attr("y", margin.top + barHeight + 15)
            .attr("text-anchor", "end")
            .text(`${max}`);
    });

    const legendBBox = legendGroup.node()?.getBBox();
    if (legendBBox) {
    const svgNode = svg.node();
    const vb = svg.attr("viewBox").split(" ").map(Number);
    let [viewBoxX, viewBoxY, viewBoxWidth, viewBoxHeight] = vb;

    const treeGroup = svg.select('g');
    const treeBBox = treeGroup.node().getBBox();

    // Get tree's transform
    const transform = treeGroup.attr("transform");
    let translateX = 0, translateY = 0;
    const match = transform?.match(/translate\(([^,]+),\s*([^)]+)\)/);
    if (match) {
        translateX = parseFloat(match[1]);
        translateY = parseFloat(match[2]);
    }

    const padding = 40;
    const legendTotalHeight = legendBBox.height + padding;
    const centerX = viewBoxX + viewBoxWidth / 2;

    let yOffset;
    if (TreeLegendPosition.toLowerCase() === "top") {
        viewBoxY -= legendTotalHeight;
        viewBoxHeight += legendTotalHeight;

        // Move tree down
        treeGroup.attr("transform", `translate(${translateX}, ${translateY})`);

        // Place legend above original tree position
        const originalTreeTop = treeBBox.y + translateY;
        yOffset = originalTreeTop - legendBBox.height - padding*1.5;
    } else {
        // Just increase viewBox height
        viewBoxHeight += legendTotalHeight;

        const treeBottom = treeBBox.y + treeBBox.height + translateY;
        yOffset = treeBottom;
    }

    // Set adjusted viewBox
    svg.attr("viewBox", `${viewBoxX} ${viewBoxY} ${viewBoxWidth} ${viewBoxHeight}`);

    // Center the legend
    const centerOffsetX = centerX - (legendBBox.x + legendBBox.width / 2);
    legendGroup.attr("transform", `translate(${centerOffsetX}, ${yOffset})`);

    }
}

/**
 * Gradient legend bars for the List mapper — mirrors createLegendBars visually.
 * @param {string}   location      - DOM id of the plot container (e.g. 'mapper-list-plot')
 * @param {object}   list_data_wow - {receptorName: {Value1:n, Value2:n, ...}, ...}
 * @param {object}   colorConfig   - {Col1:{style,color1,color2}, ...}  from mapper_list_page.js
 * @param {object}   barlabels     - {Col1:'Value 1', ...}
 * @param {object}   datatype_map  - {Col1:'Continuous'|'Discrete', ...}
 * @param {string}   position      - 'Top' or 'Bottom'
 * @param {number}   topReserved   - pixels reserved at top of SVG for bars (Top mode only)
 */
function createListLegendBars(location, list_data_wow, colorConfig, barlabels, datatype_map, position, topReserved) {
    const svg = d3.select('#' + location + ' svg');
    if (!svg.node()) return;
    svg.selectAll('.list-legend-group').remove();

    // Column definitions: internal col key → data key
    const colDef = [
        { col: 'Col1', valKey: 'Value1' },
        { col: 'Col2', valKey: 'Value2' },
        { col: 'Col3', valKey: 'Value3' },
        { col: 'Col4', valKey: 'Value4' }
    ];

    // Compute min/max per column from list_data_wow
    const stats = {};
    Object.values(list_data_wow).forEach(function (item) {
        colDef.forEach(function (d) {
            var v = item[d.valKey];
            if (typeof v === 'number' && isFinite(v)) {
                if (!stats[d.col]) { stats[d.col] = { min: v, max: v }; }
                else {
                    stats[d.col].min = Math.min(stats[d.col].min, v);
                    stats[d.col].max = Math.max(stats[d.col].max, v);
                }
            }
        });
    });

    var activeCols = colDef.filter(function (d) {
        return stats[d.col] && (datatype_map[d.col] || 'Continuous') !== 'Discrete';
    });
    if (!activeCols.length) return;

    const barWidth = 100, barHeight = 15, hSpacing = 50;
    const legendGroup = svg.append('g').attr('class', 'list-legend-group');

    activeCols.forEach(function (d, idx) {
        var cfg    = colorConfig[d.col] || {};
        var raw    = cfg.style || 'One';
        var style  = raw === 'One' ? 'One' : (raw.indexOf('Three') === 0 ? 'Three' : 'Two');
        var min    = parseFloat((stats[d.col].min).toFixed(2));
        var max    = parseFloat((stats[d.col].max).toFixed(2));
        var label  = (barlabels && barlabels[d.col]) || ('Value ' + d.col.substr(3));
        var x      = idx * (barWidth + hSpacing);
        var gradId = 'list-bar-grad-' + d.col;

        var defs = legendGroup.append('defs');
        var grad = defs.append('linearGradient')
            .attr('id', gradId).attr('x1', '0%').attr('x2', '100%').attr('y1', '0%').attr('y2', '0%');

        if (style === 'One') {
            grad.append('stop').attr('offset', '0%').attr('stop-color', '#FFFFFF');
            grad.append('stop').attr('offset', '100%').attr('stop-color', cfg.color2 || '#707070');
        } else if (style === 'Three') {
            grad.append('stop').attr('offset', '0%').attr('stop-color', cfg.color1 || '#a00000');
            grad.append('stop').attr('offset', '50%').attr('stop-color', '#FFFFFF');
            grad.append('stop').attr('offset', '100%').attr('stop-color', cfg.color2 || '#1a80bb');
        } else {
            grad.append('stop').attr('offset', '0%').attr('stop-color', cfg.color1 || '#97a6c4');
            grad.append('stop').attr('offset', '100%').attr('stop-color', cfg.color2 || '#384860');
        }

        // Gradient bar
        legendGroup.append('rect')
            .attr('x', x).attr('y', 0).attr('width', barWidth).attr('height', barHeight)
            .style('fill', 'url(#' + gradId + ')')
            .style('stroke', '#444').style('stroke-width', '0.8px');

        // Label above
        legendGroup.append('text')
            .attr('x', x + barWidth / 2).attr('y', -8)
            .attr('text-anchor', 'middle')
            .style('font-size', '11px').style('font-family', 'sans-serif')
            .text(label);

        // Min value
        legendGroup.append('text')
            .attr('x', x).attr('y', barHeight + 13)
            .attr('text-anchor', 'start')
            .style('font-size', '10px').style('font-family', 'sans-serif')
            .text(String(min));

        // Max value
        legendGroup.append('text')
            .attr('x', x + barWidth).attr('y', barHeight + 13)
            .attr('text-anchor', 'end')
            .style('font-size', '10px').style('font-family', 'sans-serif')
            .text(String(max));
    });

    // Position the legend group — center horizontally, top or bottom
    var svgW    = +svg.attr('width')  || 600;
    var svgH    = +svg.attr('height') || 400;
    var lBBox   = legendGroup.node().getBBox();
    var centerX = svgW / 2;
    var tx      = centerX - (lBBox.x + lBBox.width / 2);

    if (position === 'Bottom') {
        var plotGroup    = svg.select('.main-plot-group').node();
        var contentBottom = svgH - 20;
        if (plotGroup) {
            var pb = plotGroup.getBBox();
            contentBottom = pb.y + pb.height + 25;
        }
        var ty = contentBottom - lBBox.y;
        legendGroup.attr('transform', 'translate(' + tx + ',' + ty + ')');
        var needed = ty + lBBox.height + 20;
        if (needed > svgH) svg.attr('height', needed);
    } else {
        // Top: place within the reserved band at the top of the SVG
        var reserved = topReserved || 50;
        var ty2 = Math.max(8, (reserved - (lBBox.height)) / 2) - lBBox.y;
        legendGroup.attr('transform', 'translate(' + tx + ',' + ty2 + ')');
    }
}

// Updated CreateTextLegend for Tree Plot with layout/sorting logic and centering
function CreateTextLegend(location, circle_data, Layout) {
    const layoutMode = Layout.layoutMode || "row";
    const columns = Layout.columns || 2;
    const sortDirection = Layout.sortDirection || "Vertically";
    const TreeLegendPosition = Layout.TreeLegendPosition || "bottom";
    const LegendFontSize = Layout.Fontsize || "11px";

    const svg = d3.select(`#${location} svg`);
    svg.selectAll(".legend-text-categories").remove();

    const legendGroup = svg.append("g").attr("class", "legend-text-categories");

    const LegendFontStyle = "Arial";
    const spacingY = 25;
    const padding = 10;

    const legendItems = new Set();

    Object.entries(circle_data).forEach(([key, value]) => {
        if (value.Inner && value.ColorValue) {
            const pairKey = `${value.ColorValue}|||${value.Inner}`;
            legendItems.add(pairKey);
        }
    });

    const sortedItems = Array.from(legendItems)
        .filter(pair => {
            const [color, label] = pair.split("|||");
            return !(color === "#FFFFFF" && label === "Empty");
        })
        .sort((a, b) => {
            const labelA = (a.split("|||")[1] || "");
            const labelB = (b.split("|||")[1] || "");
            return labelA.localeCompare(labelB, undefined, { numeric: true, sensitivity: 'base' });
        });

    const tempText = svg.append("text")
        .attr("x", -9999).attr("y", -9999)
        .style("font-size", LegendFontSize).style("font-family", LegendFontStyle);

    const svgNode = svg.node();
    let svgWidth = 800;
    let vb = svg.attr("viewBox");
    if (vb) {
        const [, , w] = vb.split(" ").map(Number);
        svgWidth = w;
    }

    let x = 50, y = 40;

    if (layoutMode === "row") {
        const maxItemWidth = svgWidth - 100;

        sortedItems.forEach(pairKey => {
            let [color, label] = pairKey.split("|||");
            label = label || "";

            tempText.text(label);
            const labelWidth = tempText.node()?.getComputedTextLength() || 0;
            const totalWidth = 6 * 2 + padding + labelWidth + 20;

            if (totalWidth > maxItemWidth) {
                label = "⚠ Too long label";
                tempText.text(label);
            }

            const fixedWidth = 6 * 2 + padding + (tempText.node()?.getComputedTextLength() || 0) + 20;

            if (x + fixedWidth > svgWidth - 50) {
                x = 50;
                y += spacingY;
            }

            const fontSizeNumber = parseFloat(LegendFontSize) || 11;
            const circleRadius = Math.round(fontSizeNumber * 0.4); // e.g. 5px for 12px font

            const centerY = y + fontSizeNumber * 0.35; // 0.35 gives a good vertical centering empirically

            legendGroup.append("circle")
                .attr("cx", x)
                .attr("cy", centerY)
                .attr("r", circleRadius)
                .style("fill", color)
                .style("stroke", "black");

            legendGroup.append("text")
                .attr("x", x + circleRadius + 6)
                .attr("y", centerY)
                .attr("dy", "0.35em") // shift baseline to middle consistently
                .style("font-size", LegendFontSize)
                .style("font-family", LegendFontStyle)
                .text(label);

            x += fixedWidth;
        });

    } else if (layoutMode === "columns") {
        const numCols = columns;
        let colData = [];

        if (sortDirection === "Horizontally") {
            colData = Array.from({ length: numCols }, () => []);
            sortedItems.forEach((item, index) => {
                const colIndex = index % numCols;
                colData[colIndex].push(item);
            });
        } else {
            const perCol = Math.floor(sortedItems.length / numCols);
            const remainder = sortedItems.length % numCols;
            let index = 0;
            for (let i = 0; i < numCols; i++) {
                const count = perCol + (i < remainder ? 1 : 0);
                colData.push(sortedItems.slice(index, index + count));
                index += count;
            }
        }

        const colWidths = colData.map(col => {
            let maxWidth = 0;
            col.forEach(pair => {
                const label = pair.split("|||")[1] || "";
                tempText.text(label);
                const width = tempText.node()?.getComputedTextLength() || 0;
                maxWidth = Math.max(maxWidth, width);
            });
            return maxWidth + 40;
        });

        let startX = 50;
        let startY = 40;

        colData.forEach((col, colIndex) => {
            let x = startX;
            let y = startY;

            col.forEach(pair => {
                const [color, labelRaw] = pair.split("|||");
                const label = labelRaw || "";

                const fontSizeNumber = parseFloat(LegendFontSize) || 11;
                const circleRadius = Math.round(fontSizeNumber * 0.4); // e.g. 5px for 12px font

                const centerY = y + fontSizeNumber * 0.35; // 0.35 gives a good vertical centering empirically

                legendGroup.append("circle")
                    .attr("cx", x)
                    .attr("cy", centerY)
                    .attr("r", circleRadius)
                    .style("fill", color)
                    .style("stroke", "black");

                legendGroup.append("text")
                    .attr("x", x + circleRadius + 6)
                    .attr("y", centerY)
                    .attr("dy", "0.35em") // shift baseline to middle consistently
                    .style("font-size", LegendFontSize)
                    .style("font-family", LegendFontStyle)
                    .text(label);

                y += spacingY;
            });

            startX += colWidths[colIndex];
        });
    }

    tempText.remove();

    const legendBBox = legendGroup.node()?.getBBox();
    if (legendBBox) {
        const centerOffsetX = (svgWidth - legendBBox.width) / 2 - legendBBox.x;
        let yOffset = 0;

        const paddingBelowLegend = 40;
        const extraHeight = legendBBox.height + paddingBelowLegend;

        if (vb) {
            const parts = vb.split(" ").map(Number);

            if (TreeLegendPosition === "Bottom") {
                // Expand viewBox downward, place legend below tree
                parts[3] += extraHeight;
                svg.attr("viewBox", parts.join(" "));
                yOffset = parts[3] - extraHeight;
            } else {
                // Top: place legend above tree's actual bounding box.
                // Use svg.select('g') — same approach as createLegendBars which is known to work.
                const treeGroup = svg.select('g');
                const treeBBoxTop = treeGroup.node() ? treeGroup.node().getBBox() : null;
                let translateX = 0, translateY = 0;
                const treeTransform = treeGroup.attr('transform');
                const treeMatch = treeTransform && treeTransform.match(/translate\(([^,]+),([^)]+)\)/);
                if (treeMatch) {
                    translateX = parseFloat(treeMatch[1]);
                    translateY = parseFloat(treeMatch[2]);
                }
                // Expand viewBox upward so legend is visible above the tree
                parts[1] -= extraHeight;
                parts[3] += extraHeight;
                svg.attr("viewBox", parts.join(" "));
                // treeTop is the topmost Y of all tree content in SVG coordinates.
                // legendBBox.y is the local offset of legend content within the group (~28-40px),
                // so we must subtract it to correctly align the legend bottom with treeTop - gap.
                const treeTop = (treeBBoxTop ? treeBBoxTop.y : 0) + translateY;
                yOffset = treeTop - legendBBox.y - legendBBox.height - 10;
            }
        }

        legendGroup.attr("transform", `translate(${centerOffsetX}, ${yOffset})`);
    }
}



// #################
// ###  CLUSTER  ###
// #################

// Function to naturally sort an array
function naturalSort(a, b) {
    return a.localeCompare(b, undefined, { numeric: true, sensitivity: 'base' });
}

function decodeHtmlEntities(text) {
    const textarea = document.createElement("textarea");
    textarea.innerHTML = text;
    return textarea.value;
}

// Function to create traces for the plot
function createTraces(colorOption, showLabels, colorMapping, textColorEnabled) {
    let traces = [];
    const stroke_width = cluster_DataStyling.strokeWidth;  // Get stroke width from data styling
    const marker_size = cluster_DataStyling.markerSize;    // Get marker size from data styling
    const border_on = cluster_DataStyling.borderOn;        // Get border toggle state from data styling

    const clusterSet = Array.from(new Set(currentClusterData.map(d => d.cluster))).sort((a, b) => a - b);

    // Handle clustering color option
    if (colorOption === 'cluster') {
        clusterSet.forEach(cluster => {
            const clusterData = currentClusterData.filter(d => d.cluster === cluster);

            // Create marker traces, only display markers when labels are not shown
            const markerTrace = {
                x: clusterData.map(d => d.x),
                y: clusterData.map(d => d.y),
                mode: showLabels ? 'none' : 'markers',  // Hide markers if labels are shown
                type: 'scatter',
                hoverinfo: 'text',
                text: clusterData.map(d => `${decodeHtmlEntities(d.label)}`),  // Combine label and fill for hover text
                marker: {
                    size: marker_size,
                    symbol: 'circle',
                    line: {
                        width: border_on ? stroke_width : 0,
                        color: 'black'
                    },
                    color: clusterData.map(d => {
                        const key = 'Cluster ' + (d.cluster + 1);
                        if (CLUSTER_LABEL_COLORS && !CLUSTER_LABEL_COLORS[key]) {
                            CLUSTER_LABEL_COLORS[key] = colorPalette[d.cluster % colorPalette.length];
                        }
                        return (CLUSTER_LABEL_COLORS && CLUSTER_LABEL_COLORS[key]) || colorPalette[d.cluster % colorPalette.length];
                    }),
                },
                // Customize legend text color only when showLabels and textColorEnabled are both true
                name: (showLabels && textColorEnabled)
                ? `<span style="color:${(CLUSTER_LABEL_COLORS && CLUSTER_LABEL_COLORS['Cluster ' + (cluster + 1)]) || colorPalette[cluster % colorPalette.length]}">Cluster ${cluster + 1}</span>`
                : `Cluster ${cluster + 1}`
            };

            traces.push(markerTrace);
        });
    }
    // Handle gradient color option
    else if (colorOption === 'gradient') {
        const fillValues = currentClusterData.map(d => d.fill);
        const minFill = Math.min(...fillValues);
        const maxFill = Math.max(...fillValues);
        const gradientTrace = {
            x: currentClusterData.map(d => d.x),
            y: currentClusterData.map(d => d.y),
            mode: showLabels ? 'none' : 'markers',  // Hide markers if labels are shown
            type: 'scatter',
            hoverinfo: 'text',
            text: currentClusterData.map(d => `${decodeHtmlEntities(d.label)}: ${d.fill}`),  // Combine label and fill for hover text
            marker: {
                size: marker_size,
                symbol: 'circle',
                line: {
                    width: border_on ? stroke_width : 0,
                    color: 'black'
                },
                color: currentClusterData.map(d => d.fill),
                cmin: minFill,
                cmax: maxFill,
                colorscale: (function() {
                    if (!CLUSTER_GRADIENT_COLORS) return 'RdBu';
                    var stops = window.CLUSTER_GRADIENT_STOPS || 3;
                    if (stops === 1) return [[0, '#ffffff'], [1, CLUSTER_GRADIENT_COLORS.max]];
                    if (stops === 2) return [[0, CLUSTER_GRADIENT_COLORS.min], [1, CLUSTER_GRADIENT_COLORS.max]];
                    return [[0, CLUSTER_GRADIENT_COLORS.min], [0.5, CLUSTER_GRADIENT_COLORS.mid], [1, CLUSTER_GRADIENT_COLORS.max]];
                })(),
                colorbar: {
                    thickness: 20,
                    len: 0.5,
                }
            },
            showlegend: false
        };

        traces.push(gradientTrace);
    }
    // Handle Class, Ligand type, or Receptor family color option
    else if (['Class', 'Ligand type', 'Receptor family', 'userCategory'].includes(colorOption)) {

        const uniqueEntries = new Set();

        currentClusterData.forEach(point => {
            const entry = point[colorOption];
            if (!uniqueEntries.has(entry)) {
                uniqueEntries.add(entry);

                // Create marker traces, only display markers when labels are not shown
                const markerTrace = {
                    x: currentClusterData.filter(d => d[colorOption] === point[colorOption]).map(d => d.x),
                    y: currentClusterData.filter(d => d[colorOption] === point[colorOption]).map(d => d.y),
                    mode: showLabels ? 'none' : 'markers',
                    type: 'scatter',
                    hoverinfo: 'text',
                    text: currentClusterData.filter(d => d[colorOption] === point[colorOption]).map(d => d.label.toUpperCase()),
                    marker: {
                        size: marker_size,
                        symbol: 'circle',
                        line: {
                            width: border_on ? stroke_width : 0,
                            color: 'black'
                        },
                        color: colorMapping[entry],
                    },
                    // Customize legend text color only when showLabels and textColorEnabled are both true
                    name: (showLabels && textColorEnabled)
                        ? `<span style="color:${colorMapping[entry]}">${entry}</span>`
                        : entry  // Regular label if conditions are false
                };

                traces.push(markerTrace);
            }
        });

        traces.sort((a, b) => naturalSort(a.name, b.name));
    }

    return traces;
}


// Function to create annotations for labels with optional text coloring based on `textColorEnabled`
function createAnnotations(filteredData, colorOption, textColorEnabled, colorMapping) {
    const annotations = [];

    // Get global min and max values for `fill` from all data (currentClusterData), not just the filtered data
    const fillValues = currentClusterData.map(d => d.fill);  // Use all data, not just filtered
    const minFill = Math.min(...fillValues);
    const maxFill = Math.max(...fillValues);

    // Build the gradient color scale (custom or default RdBu)
    let gradColorScale;
    if (!CLUSTER_GRADIENT_COLORS) {
        gradColorScale = d3v4.scaleSequential(d3v4.interpolateRdBu).domain([maxFill, minFill]);
    } else {
        const gradStops = window.CLUSTER_GRADIENT_STOPS || 3;
        if (gradStops === 1) {
            gradColorScale = d3v4.scaleLinear().domain([minFill, maxFill])
                .range(['#ffffff', CLUSTER_GRADIENT_COLORS.max]);
        } else if (gradStops === 2) {
            gradColorScale = d3v4.scaleLinear().domain([minFill, maxFill])
                .range([CLUSTER_GRADIENT_COLORS.min, CLUSTER_GRADIENT_COLORS.max]);
        } else {
            gradColorScale = d3v4.scaleLinear().domain([minFill, (minFill + maxFill) / 2, maxFill])
                .range([CLUSTER_GRADIENT_COLORS.min, CLUSTER_GRADIENT_COLORS.mid, CLUSTER_GRADIENT_COLORS.max]);
        }
    }

    filteredData.forEach((d) => {
        let textColor;

        if (textColorEnabled) {
            if (colorOption === 'cluster') {
                const clKey = 'Cluster ' + (d.cluster + 1);
                textColor = (CLUSTER_LABEL_COLORS && CLUSTER_LABEL_COLORS[clKey]) || colorPalette[d.cluster % colorPalette.length];
            } else if (colorOption === 'gradient') {
                textColor = gradColorScale(d.fill);
            } else if (['Class', 'Ligand type', 'Receptor family', 'userCategory'].includes(colorOption)) {
                const entry = d[colorOption];
                textColor = colorMapping[entry];
            } else {
                textColor = 'black';
            }
        } else {
            textColor = 'black';
        }

        annotations.push({
            x: d.x,
            y: d.y,
            xref: 'x',
            yref: 'y',
            text: `${decodeHtmlEntities(d.label)}`,
            showarrow: false,
            xanchor: 'center',
            yanchor: 'middle',
            xshift: 0,
            font: {
                family: 'Arial',
                size: cluster_DataStyling.labelFontSize,
                color: textColor
            },
            align: 'center',
            bgcolor: 'rgba(0,0,0,0)',
            borderwidth: 0
        });
    });

    return annotations;
}

// Function to update the plot with markers or labels
function updatePlotWithAnnotations() {
    const colorOption = getActiveColorOption();
    const showLabels = labelsVisible;
    const textColorEnabled = cluster_DataStyling.textColorEnabled;

    // Data bounds
    const xDataMin = Math.min(...currentClusterData.map(d => d.x));
    const xDataMax = Math.max(...currentClusterData.map(d => d.x));
    const yDataMin = Math.min(...currentClusterData.map(d => d.y));
    const yDataMax = Math.max(...currentClusterData.map(d => d.y));
    const xSpan = (xDataMax - xDataMin) || 1;
    const ySpan = (yDataMax - yDataMin) || 1;

    // Compute padding: in text mode measure real label pixel widths via canvas,
    // then use closed-form to guarantee labels at extreme positions stay inside the viewport.
    // Formula derivation: if plot area = P px and data span = S, then after adding xPad on each side
    // the scale is P/(S+2*xPad) px/unit. We need halfLabelPx <= xPad * P/(S+2*xPad), which gives:
    //   xPad = halfLabelPx * S / (P - fullLabelPx)
    let xPad, yPad;
    if (showLabels && !_clusterZoomRange) {
        const _pGd = document.getElementById('plotContainer_cluster');
        const _pW = (_pGd && _pGd.offsetWidth  > 50) ? _pGd.offsetWidth  : 600;
        const _pH = (_pGd && _pGd.offsetHeight > 50) ? _pGd.offsetHeight : _pW * 0.75;
        // Effective Plotly plot area: ~68% of container width (accounts for margins + legend),
        // ~82% of container height (top/bottom margins)
        const plotAreaW = _pW * 0.68;
        const plotAreaH = _pH * 0.82;
        let maxLabelPx = cluster_DataStyling.labelFontSize * 5;
        try {
            const _cv = document.createElement('canvas');
            const _cx = _cv.getContext('2d');
            _cx.font = `${cluster_DataStyling.labelFontSize}px Arial`;
            maxLabelPx = Math.max(...currentClusterData.map(d => _cx.measureText(d.label || '').width));
        } catch(e) {}
        const labelH = cluster_DataStyling.labelFontSize * 1.2;
        // Clamp denominator so it never goes below 10% of plot area (prevents explosion on very long labels)
        const denomX = Math.max(plotAreaW - maxLabelPx, plotAreaW * 0.10);
        const denomY = Math.max(plotAreaH - labelH,     plotAreaH * 0.10);
        xPad = Math.max(xSpan * 0.08, (maxLabelPx / 2) * xSpan / denomX);
        yPad = Math.max(ySpan * 0.08, (labelH     / 2) * ySpan / denomY);
    } else {
        xPad = xSpan * 0.08;
        yPad = ySpan * 0.08;
    }

    // Axis viewport: preserve user zoom if active, else full data range with padding
    const zoom = typeof _clusterZoomRange !== 'undefined' ? _clusterZoomRange : null;
    const xRange = zoom ? zoom[0] : [xDataMin - xPad, xDataMax + xPad];
    const yRange = zoom ? zoom[1] : [yDataMin - yPad, yDataMax + yPad];

    const filterXRange = xRange;
    const filterYRange = yRange;

    let colorMapping = {};
    if (['Class', 'Ligand type', 'Receptor family', 'userCategory'].includes(colorOption)) {
        let uniqueValues = Array.from(new Set(currentClusterData.map(d => d[colorOption])));
        uniqueValues.sort(naturalSort);
        uniqueValues.forEach((value, index) => {
            if (CLUSTER_LABEL_COLORS && !CLUSTER_LABEL_COLORS[value]) {
                CLUSTER_LABEL_COLORS[value] = colorPalette[index % colorPalette.length];
            }
            colorMapping[value] = (CLUSTER_LABEL_COLORS && CLUSTER_LABEL_COLORS[value])
                || colorPalette[index % colorPalette.length];
        });
    }

    const traces = createTraces(colorOption, showLabels, colorMapping, textColorEnabled);

    // Color bar only shown in Text mode with gradient coloring AND text color enabled
    if (colorOption === 'gradient' && showLabels && textColorEnabled) {
        const fillValues = currentClusterData.map(d => d.fill);
        const minFill = Math.min(...fillValues);
        const maxFill = Math.max(...fillValues);
        traces.push({
            z: [[minFill, maxFill], [minFill, maxFill]],
            x: [0, 1], y: [0, 1],
            type: 'heatmap',
            colorscale: (function() {
                if (!CLUSTER_GRADIENT_COLORS) return 'RdBu';
                var stops = window.CLUSTER_GRADIENT_STOPS || 3;
                if (stops === 1) return [[0, '#ffffff'], [1, CLUSTER_GRADIENT_COLORS.max]];
                if (stops === 2) return [[0, CLUSTER_GRADIENT_COLORS.min], [1, CLUSTER_GRADIENT_COLORS.max]];
                return [[0, CLUSTER_GRADIENT_COLORS.min], [0.5, CLUSTER_GRADIENT_COLORS.mid], [1, CLUSTER_GRADIENT_COLORS.max]];
            })(),
            showscale: true,
            colorbar: { thickness: 20, len: 0.5 },
            opacity: 0
        });
    }

    // Annotations only in Text mode; filter to visible viewport
    const filteredData = currentClusterData.filter(d =>
        d.x >= filterXRange[0] && d.x <= filterXRange[1] &&
        d.y >= filterYRange[0] && d.y <= filterYRange[1]
    );
    const annotations = showLabels ? createAnnotations(filteredData, colorOption, textColorEnabled, colorMapping) : [];

    const xAxisCfg = { visible: false, showgrid: false, range: xRange };
    const yAxisCfg = { visible: false, showgrid: false, range: yRange };

    // Always read the container's current pixel width so Plotly.react
    // doesn't reuse a stale stored width from a previous render.
    const _plotGd = document.getElementById('plotContainer_cluster');
    const _containerW = (_plotGd && _plotGd.offsetWidth > 50) ? _plotGd.offsetWidth : undefined;

    const layout = {
        xaxis: xAxisCfg,
        yaxis: yAxisCfg,
        hovermode: 'closest',
        showlegend: !(showLabels && !textColorEnabled),
        annotations: annotations,
        legend: {
            x: 1.01, xanchor: 'left',
            y: 0.5, yanchor: 'middle',
            orientation: 'v', font: { size: 11 }, tracegroupgap: 2,
            bgcolor: 'rgba(0,0,0,0)', borderwidth: 0
        },
        plot_bgcolor: '#FFFFFF',
        autosize: true,
        width: _containerW,
        height: 650,
        margin: { l: 20, r: 80, t: 30, b: 20 },
        shapes: [{
            type: 'rect', xref: 'paper', yref: 'paper',
            x0: 0, y0: 0, x1: 1, y1: 1,
            line: { color: 'black', width: 1 },
            fillcolor: 'rgba(0,0,0,0)'
        }]
    };

    Plotly.react('plotContainer_cluster', traces, layout, { responsive: true });

    // Clip annotation text at the plot-area boundary — but NOT the legend.
    // Both annotations and the legend live in .infolayer, so we clip the individual
    // .annotation groups rather than the whole layer.
    requestAnimationFrame(function () {
        var gd = document.getElementById('plotContainer_cluster');
        if (!gd || !gd._fullLayout) return;

        // Measure the actual right-side element (legend for discrete modes,
        // colorbar for gradient mode) and fit the right margin to it with a
        // 12 px gap.  Only relayout when the difference is > 4 px to prevent
        // an infinite loop.
        var rightEl = gd.querySelector('g.legend') || gd.querySelector('g.colorbar');
        if (rightEl) {
            try {
                var rBB = rightEl.getBBox();
                if (rBB.width > 0) {
                    var neededR = Math.ceil(rBB.width) + 12;
                    if (Math.abs(neededR - gd._fullLayout.margin.r) > 4) {
                        Plotly.relayout(gd, { 'margin.r': neededR });
                        // Annotation clipping will run on the next call once the
                        // margin has settled — skip it this frame.
                        return;
                    }
                }
            } catch (e) { /* getBBox can fail on hidden elements */ }
        }

        var fl = gd._fullLayout;
        var l = fl.margin.l, t = fl.margin.t;
        var w = fl.width - fl.margin.l - fl.margin.r;
        var h = fl.height - fl.margin.t - fl.margin.b;
        if (w <= 0 || h <= 0) return;
        var svg = gd.querySelector('svg.main-svg');
        if (!svg) return;
        var clipId = 'cluster-annotation-clip';
        var defs = svg.querySelector('defs') || svg.insertBefore(document.createElementNS('http://www.w3.org/2000/svg', 'defs'), svg.firstChild);
        var cp = defs.querySelector('#' + clipId);
        if (!cp) {
            cp = document.createElementNS('http://www.w3.org/2000/svg', 'clipPath');
            cp.setAttribute('id', clipId);
            defs.appendChild(cp);
        }
        cp.innerHTML = '<rect x="' + l + '" y="' + t + '" width="' + w + '" height="' + h + '"/>';
        // Clip each annotation group individually — legend (.legend) is left untouched
        gd.querySelectorAll('.infolayer .annotation').forEach(function (ann) {
            ann.setAttribute('clip-path', 'url(#' + clipId + ')');
        });
        // Ensure the layer itself is never clipped (would cut the legend)
        var infolayer = gd.querySelector('.infolayer');
        if (infolayer) infolayer.removeAttribute('clip-path');
    });

    // Notify cluster page (Colors panel refresh, etc.)
    if (typeof window.mapperClusterOnPlotUpdated === 'function') window.mapperClusterOnPlotUpdated();
}


// #################
// ###   LIST    ###
// #################

// ## Initialization of global data for list ##
function initializeDataStyling(list_data_wow, data_types_list) {

    // Function to check if any of the specified keys are present in a given object
    function checkKeysPresent(obj, keys) {
        return keys.some(k => k in obj);
    }

    // Define the keys to check for each column
    var col1Keys = ["Value1"];
    var col2Keys = ["Value2"];
    var col3Keys = ["Value3"];
    var col4Keys = ["Value4"];


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
            if (data_types_list[col] === 'Discrete') {
                Data_styling[col].Datatype = 'Discrete';
            } else if (data_types_list[col] === 'Continuous') {
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

        updateDataMinMax('Col1', 'Value1', currentObj);
        updateDataMinMax('Col2', 'Value2', currentObj);
        updateDataMinMax('Col3', 'Value3', currentObj);
        updateDataMinMax('Col4', 'Value4', currentObj);
    });

    return Data_styling;
}

// Function to modify the data
function Initialize_Data(data) {
    // // 3. Remove "Olfactory receptors" from Layer2
    // Object.keys(data).forEach(layer1Key => {
    //   const layer2Data = data[layer1Key];
    //   if (layer2Data["Olfactory receptors"]) {
    //     delete layer2Data["Olfactory receptors"];
    //   }
    // });

    return data;  // Return the modified data
  }

//   Create Classification dict
function Create_classification_dict(data) {
    // Initialize the arrays for present items
    let present_class = [];
    let present_LigandType = [];
    let present_ReceptorFamilies = [];
    let present_receptors = [];

    // Iterate over the nested dictionary to split into the respective arrays
    for (const classKey in data) {
        present_class.push(classKey);  // Add the class key to the present_class array

        for (const ligandKey in data[classKey]) {
            present_LigandType.push(ligandKey);  // Add the ligand type key to the present_LigandType array

            for (const receptorFamilyKey in data[classKey][ligandKey]) {
                present_ReceptorFamilies.push(receptorFamilyKey);  // Add receptor family key to present_ReceptorFamilies

                // Add the receptors to the present_receptors array
                present_receptors.push(...data[classKey][ligandKey][receptorFamilyKey]);
            }
        }
    }

}

// Create array of sorted entries and classification array for printing
function Data_resorter(data) {

    // Check the visibility of each layer
    const Layer1_isChecked = d3.select(`#toggle-layer-1`).property('checked');
    const Layer2_isChecked = d3.select(`#toggle-layer-2`).property('checked');
    const Layer3_isChecked = d3.select(`#toggle-layer-3`).property('checked');

    // Data sorter depending on layers being shown:

    // Initialize the final array for results and the category array for categories
    let final_array = [];
    let category_array = [];

    // Helper function to sort arrays naturally
    function naturalSort(a, b) {
        return a.localeCompare(b, undefined, { numeric: true, sensitivity: 'base' });
    }

    // Build the sorted Class list (only include it if Layer1 is checked)
    let sortedClasses = Layer1_isChecked ? Object.keys(data).sort(naturalSort) : [];

    // Iterate over sorted Classes and handle layers accordingly
    for (const classKey of sortedClasses) {
        final_array.push(classKey); // Add the Class if Layer1 is checked
        category_array.push("Class"); // Add "Class" category

        let ligandTypeDict = {};

        // Collect all LigandTypes and their families/receptors for the current Class
        for (const ligandKey in data[classKey]) {
            if (!ligandTypeDict[ligandKey]) {
                ligandTypeDict[ligandKey] = [];
            }

            for (const receptorFamilyKey in data[classKey][ligandKey]) {
                ligandTypeDict[ligandKey].push({
                    family: receptorFamilyKey,
                    receptors: data[classKey][ligandKey][receptorFamilyKey]
                });
            }
        }

        // Sort LigandTypes if Layer2 is checked, otherwise skip LigandTypes
        let sortedLigandTypes = Layer2_isChecked ? Object.keys(ligandTypeDict).sort(naturalSort) : [];

        // For each LigandType, process its ReceptorFamilies and Receptors
        sortedLigandTypes.forEach(ligandKey => {
            final_array.push(ligandKey); // Add LigandType if Layer2 is checked
            category_array.push("LigandType"); // Add "LigandType" category

            let receptorFamilies = ligandTypeDict[ligandKey].sort((a, b) => a.family.localeCompare(b.family)); // Sort ReceptorFamilies

            receptorFamilies.forEach(familyObj => {
                if (Layer3_isChecked) {
                    final_array.push(familyObj.family); // Add ReceptorFamily if Layer3 is checked
                    category_array.push("ReceptorFamily"); // Add "ReceptorFamily" category
                }

                // Always sort and add receptors
                familyObj.receptors.forEach(receptor => {
                    final_array.push(receptor);
                    category_array.push("Receptor"); // Add "Receptor" category for each receptor
                });
            });
        });

        // If Layer2 is not checked but Layer3 is checked (Case [True, False, True])
        if (!Layer2_isChecked && Layer3_isChecked) {
            let receptorFamilies = [];
            for (const ligandKey in ligandTypeDict) {
                receptorFamilies.push(...ligandTypeDict[ligandKey]);
            }

            receptorFamilies.sort((a, b) => a.family.localeCompare(b.family));

            receptorFamilies.forEach(familyObj => {
                final_array.push(familyObj.family); // Add ReceptorFamily if Layer3 is checked
                category_array.push("ReceptorFamily"); // Add "ReceptorFamily" category
                familyObj.receptors.forEach(receptor => {
                    final_array.push(receptor);
                    category_array.push("Receptor"); // Add "Receptor" category
                });
            });
        }

        if (!Layer2_isChecked && !Layer3_isChecked) {
            for (const ligandKey in ligandTypeDict) {
                ligandTypeDict[ligandKey].forEach(familyObj => {
                    familyObj.receptors.forEach(receptor => {
                        final_array.push(receptor); // Add sorted Receptors
                        category_array.push("Receptor"); // Add "Receptor" category
                    });
                });
            }
        }
    }

    // Case [False, True, True]: Include LigandType, ReceptorFamily, and Receptors, but no Class
    if (!Layer1_isChecked && Layer2_isChecked && Layer3_isChecked) {
        let allLigandTypes = [];

        // Collect LigandTypes and ReceptorFamilies across all classes
        for (const classKey in data) {
            for (const ligandKey in data[classKey]) {
                if (!allLigandTypes[ligandKey]) {
                    allLigandTypes[ligandKey] = [];
                }
                for (const receptorFamilyKey in data[classKey][ligandKey]) {
                    allLigandTypes[ligandKey].push({
                        family: receptorFamilyKey,
                        receptors: data[classKey][ligandKey][receptorFamilyKey]
                    });
                }
            }
        }

        let sortedLigandTypes = Object.keys(allLigandTypes).sort(naturalSort);

        sortedLigandTypes.forEach(ligandKey => {
            final_array.push(ligandKey); // Add LigandType
            category_array.push("LigandType"); // Add "LigandType" category

            let receptorFamilies = allLigandTypes[ligandKey].sort((a, b) => a.family.localeCompare(b.family));

            receptorFamilies.forEach(familyObj => {
                final_array.push(familyObj.family); // Add ReceptorFamily
                category_array.push("ReceptorFamily"); // Add "ReceptorFamily" category

                familyObj.receptors.forEach(receptor => {
                    final_array.push(receptor);
                    category_array.push("Receptor"); // Add "Receptor" category
                });
            });
        });
    }

    // Case [False, True, False]: Include LigandType and Receptors, but no Class or ReceptorFamily
    if (!Layer1_isChecked && Layer2_isChecked && !Layer3_isChecked) {
        let allLigandTypes = [];

        // Collect LigandTypes and Receptors across all classes
        for (const classKey in data) {
            for (const ligandKey in data[classKey]) {
                if (!allLigandTypes[ligandKey]) {
                    allLigandTypes[ligandKey] = [];
                }
                for (const receptorFamilyKey in data[classKey][ligandKey]) {
                    allLigandTypes[ligandKey].push(...data[classKey][ligandKey][receptorFamilyKey]);
                }
            }
        }

        let sortedLigandTypes = Object.keys(allLigandTypes).sort(naturalSort);

        sortedLigandTypes.forEach(ligandKey => {
            final_array.push(ligandKey); // Add LigandType
            category_array.push("LigandType"); // Add "LigandType" category

            allLigandTypes[ligandKey]
            .slice()                      // optional: clone to avoid mutating original
            .sort(naturalSort)
            .forEach(receptor => {
                final_array.push(receptor);
                category_array.push("Receptor");
            });
        });
    }

    // If no Class or LigandType is checked but Receptors need to be included
    if (!Layer1_isChecked && !Layer2_isChecked && Layer3_isChecked) {
        let allReceptorFamilies = [];

        for (const classKey in data) {
            for (const ligandKey in data[classKey]) {
                for (const receptorFamilyKey in data[classKey][ligandKey]) {
                    // Collect receptor families and their receptors directly from data
                    allReceptorFamilies.push({
                        family: receptorFamilyKey,
                        receptors: data[classKey][ligandKey][receptorFamilyKey]
                    });
                }
            }
        }

        // Sort the receptor families
        allReceptorFamilies.sort((a, b) => a.family.localeCompare(b.family));

        allReceptorFamilies.forEach(familyObj => {
            final_array.push(familyObj.family); // Add ReceptorFamily
            category_array.push("ReceptorFamily"); // Add "ReceptorFamily" category
            familyObj.receptors.forEach(receptor => {
                final_array.push(receptor);
                category_array.push("Receptor"); // Add "Receptor" category
            });
        });
    }


    // If all layers are unchecked, add all sorted receptors
    if (!Layer1_isChecked && !Layer2_isChecked && !Layer3_isChecked) {
        for (const classKey in data) {
            for (const ligandKey in data[classKey]) {
                for (const receptorFamilyKey in data[classKey][ligandKey]) {
                    const receptors = data[classKey][ligandKey][receptorFamilyKey];
                    receptors
                        .slice()               // optional defensive copy
                        .sort(naturalSort)
                        .forEach(receptor => {
                            final_array.push(receptor);
                            category_array.push("Receptor");
                        });
                }
            }
        }
    }
    return { final_array, category_array };
}

// === Calculate position of data and receptor labels ====
function computeDynamicOffsets(data_styling, shape_spacing = 20) {
    const offsetMap = {};
    let currentOffset = 0;

    // Get only active columns, in original order
    const activeCols = ['Col1', 'Col2', 'Col3', 'Col4'].filter(col => data_styling[col]?.Data === 'Yes');

    activeCols.forEach(col => {
        offsetMap[col] = currentOffset;
        currentOffset += shape_spacing;
    });

    const labelOffset = currentOffset;

    return {
        offsetMap,      // e.g. { Col3: 0 }
        labelOffset,    // e.g. 20
        activeCols
    };
}


// Calculate the dimensions of the plot
function Calculate_dimension(data, Category_data, Col_break_number, columns, label_conversion_dicts, label_names, styling_option) {

    // Define the column tracking variables
    let temp_col_state = 1; // Current column
    let label_dim_counter = 0; // Counter for how many labels processed in the current column
    const { labelOffset } = computeDynamicOffsets(Data_styling);

    // Initialize label max width tracking — store actual longest string per column/category
    let label_max_dict = {
        col1_label_max: { 'Class': '', 'LigandType': '', 'ReceptorFamily': '', 'Receptor': '' },
        col2_label_max: { 'Class': '', 'LigandType': '', 'ReceptorFamily': '', 'Receptor': '' },
        col3_label_max: { 'Class': '', 'LigandType': '', 'ReceptorFamily': '', 'Receptor': '' },
        col4_label_max: { 'Class': '', 'LigandType': '', 'ReceptorFamily': '', 'Receptor': '' }
    };

    let Col_spacing_dict = { 'Col1': -Infinity, 'Col2': -Infinity, 'Col3': -Infinity, 'Col4': -Infinity };

    // List of known HTML entities to replace for correct label width calculation (for IUPHAR)
    const htmlEntities = ['&alpha;', '&beta;', '&gamma;', '&delta;', '&epsilon;', '&zeta;', '&eta;', '&theta;', '&iota;', '&kappa;', '&lambda;', '&mu;', '&nu;', '&xi;', '&omicron;', '&pi;', '&rho;', '&sigma;', '&tau;', '&upsilon;', '&phi;', '&chi;', '&psi;', '&omega;'];

    // Helper function to replace HTML entities in label text (only for IUPHAR)
    function replaceHtmlEntities(str) {
        const entityRegex = new RegExp(htmlEntities.join('|'), 'gi');
        return str.replace(entityRegex, 'a'); // Replace with 'a' or some neutral character of approximate size
    }

    // Process each label in the data and check its corresponding category from Category_data
    data.forEach((label, index) => {
        const category = Category_data[index]; // Get the category of the current label ("Class", "LigandType", "ReceptorFamily", or "Receptor")

        // Increment the label dimension counter
        label_dim_counter++;

        // Move to the next column if the Col_break_number is exceeded
        if (label_dim_counter > Col_break_number && temp_col_state < columns) {
            label_dim_counter = 1; // Reset the counter
            temp_col_state++; // Move to the next column
        }

        // Trimming and processing labels based on category
        if (category === 'Class') {
            // Trimming Class (e.g., "Class (some text)" becomes "Class")
            label = label.split(" (")[0];

        } else if (category === 'LigandType') {
            // No additional processing for LigandType
            label = label;

        } else if (category === 'ReceptorFamily') {
            // No trimming — matches RenderListPlot_Labels which shows full family names
            label = label;

        } else if (category === 'Receptor') {
            // For Receptors, apply the label conversion based on 'label_names'
            if (label_names === 'UniProt') {
                // Convert based on UniProt data
                label = label_conversion_dicts.IUPHAR_to_UniProt_converter[label];
                label = label ? label.toUpperCase() : label;
            } else if (label_names === 'Protein') {
                // Apply IUPHAR-specific replacements and clean up
                label = replaceHtmlEntities(label)
                    .replace(/<\/?i>|(-adrenoceptor| receptor)|<\/?sub>/g, '');

                // Subscript handling for accurate label length measurement
                label = label.replace(/<sub>.*?<\/sub>/g, ''); // Simplified subscript removal for length estimation
            } else if (label_names === 'Gene') {
                // Convert based on UniProt data
                label = label_conversion_dicts.IUPHAR_to_Gene_converter[label];;
            }
        }

        // Determine current column
        const colKey = `col${temp_col_state}_label_max`;

        // Store the actual longest label string (strip HTML tags for length comparison)
        const cleanLabel = label ? label.replace(/<[^>]*>/g, '') : '';
        if (cleanLabel.length > label_max_dict[colKey][category].length) {
            label_max_dict[colKey][category] = cleanLabel;
        }
    });

    // Calculate the actual pixel widths for each column and category
    const categories = ['Class', 'LigandType', 'ReceptorFamily', 'Receptor'];
    const cols = ['Col1', 'Col2', 'Col3', 'Col4'];

    for (let i = 0; i < columns; i++) {
        const key = `col${i + 1}_label_max`;
        const max_label_values = label_max_dict[key];

        categories.forEach(category => {
            const columnKey = cols[i];
            const actualLabel = max_label_values[category];

            if (actualLabel !== '') {
                // Measure the actual longest label string with real font settings
                const dummyText = d3.select("body")
                    .append("svg")
                    .attr("class", "dummy-text-measure")
                    .append("text")
                    .attr("font-size", styling_option[category].Fontsize)
                    .attr("font-weight", styling_option[category].Bold ? "bold" : "normal")
                    .attr("font-style", styling_option[category].Italic ? "italic" : "normal")
                    .text(actualLabel);

                const bbox = dummyText.node().getBBox();
                // Use actual measured width + generous padding to prevent overflow
                let estimatedLength = bbox.width + 30;

                if (category === 'Receptor') {
                    estimatedLength += labelOffset; // Additional margin for data shapes
                }

                d3.select(".dummy-text-measure").remove();

                if (estimatedLength > Col_spacing_dict[columnKey]) {
                    Col_spacing_dict[columnKey] = estimatedLength;
                }
            }
        });
    }
    return Col_spacing_dict;
}

// Genreate the printing of the labels
function RenderListPlot_Labels(data, category_data, location, styling_option, Layout_dict, label_conversion_dicts, label_names,y_off_set_variable=30) {
    // ######################
    // ## Initialization   ##
    // ######################

    // Extract necessary values from the layout
    let columns = Layout_dict.columns;
    let col_max_label = Layout_dict.col_max_label;

    // Calculate total count
    const total_count = data.length;

    // Determine Col_break_number (automatic or custom)
    let Col_break_number = col_max_label === "Auto" ? Math.ceil(total_count / columns) : Layout_dict.Col_break_number;


    // ##################
    // ## Dimensions   ##
    // ##################
    let width = 0;
    const initial_height = 200;  // Starting height of the SVG
    const margin = { top: 40, right: 20, bottom: 20, left: 20 };

    // ## Calculate all spacing and dimensions ##
    let spacing_dict = Calculate_dimension(data, category_data, Col_break_number, columns, label_conversion_dicts, label_names, styling_option);

    // Calculate total width based on spacing_dict and columns
    const col_list = ['Col1', 'Col2', 'Col3', 'Col4'];
    for (let i = 0; i < columns; i++) {
         if (spacing_dict[col_list[i]] === -Infinity) {
            spacing_dict[col_list[i]] = 0; // or use 20 if you want a fixed minimum column width
        }
        width += spacing_dict[col_list[i]];
    }
    width += 45; // Add some padding

    // // Recalculate width if its less than legend bars
    let number_of_bar_legends = 0;
    Object.keys(Data_styling).forEach(function(column) {
        if (Data_styling[column].Data === "Yes" && Data_styling[column].Datatype === 'Continuous') {
            number_of_bar_legends++;
        }
    });
    let minimum_width = 0;
    if (number_of_bar_legends > 0) {
        minimum_width = 230*number_of_bar_legends+margin.left+margin.right;
        if (width < minimum_width) {
            width = minimum_width;
        }
    }

    // Initialize SVG container
    const svg = d3.select(`#${location}`)
        .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", initial_height + margin.top + margin.bottom)  // Start with initial height
        .attr("id", "ListPlot_visualization");

    // Adding correct group for appending
    const plotGroup = svg.append("g").attr("class", "main-plot-group");

    // Set initial X & Y positions — reserve top space for gradient bars unless they go to the bottom
    const _legendTopSpace = (Layout_dict && Layout_dict.legend_top_space != null)
        ? Layout_dict.legend_top_space : 50;
    let yOffset = margin.top + 5 + _legendTopSpace;
    let yOffset_max = yOffset; // Track the maximum yOffset
    let xOffset = 0;
    const { labelOffset } = computeDynamicOffsets(Data_styling);

    // ######################
    // ## Label Handling   ##
    // ######################

    // Function to handle text rendering with special cases for IUPHAR and UniProt (only for Receptor category)
    function add_text(label, yOffset, category, label_key, bold = false, italic = false, underline = false, fontSize = '16px', color = 'black', label_names) {
        // Handle special cases for Receptors
        if (category === 'Receptor') {
            if (label_names === 'Protein') {

                const htmlEntities = {
                    '&alpha;': 'α', '&beta;': 'β', '&gamma;': 'γ', '&delta;': 'δ',
                    '&epsilon;': 'ε', '&zeta;': 'ζ', '&eta;': 'η', '&theta;': 'θ',
                    '&iota;': 'ι', '&kappa;': 'κ', '&lambda;': 'λ', '&mu;': 'μ',
                    '&nu;': 'ν', '&xi;': 'ξ', '&omicron;': 'ο', '&pi;': 'π',
                    '&rho;': 'ρ', '&sigma;': 'σ', '&tau;': 'τ', '&upsilon;': 'υ',
                    '&phi;': 'φ', '&chi;': 'χ', '&psi;': 'ψ', '&omega;': 'ω'
                };

                // Replace HTML entities with corresponding Unicode characters
                for (const [entity, char] of Object.entries(htmlEntities)) {
                    label = label.replace(new RegExp(entity, 'g'), char);
                }

                // Remove <i> tags and (-adrenoceptor| receptor) from the label
                label = label.replace(/<\/?i>/g, '')  // Remove <i> tags
                .replace(/(-adrenoceptor| receptor)/g, '');  // Remove adrenoceptor/receptor text


                // Create text element for IUPHAR receptor
                const textElement = plotGroup.append('text')
                    .attr('x', margin.left + xOffset + labelOffset)
                    .attr('y', yOffset)
                    .attr('class', category)
                    .style('dominant-baseline', 'middle') // Set vertical alignment
                    .style('font-weight', bold ? 'bold' : 'normal')
                    .style('font-style', italic ? 'italic' : 'normal')
                    .style('text-decoration', underline ? 'underline' : 'none')
                    .style('font-size', fontSize)
                    .style('fill', color);

                // Subscript handling for IUPHAR receptors
                const mainFontSize = parseFloat(fontSize);
                const subFontSize = mainFontSize * 0.75 + 'px';
                const parts = label.split(/(<sub>|<\/sub>)/);
                let inSub = false;

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

            } else if (label_names === 'UniProt') {
                // Handle UniProt receptor labels
                label = (label_conversion_dicts.IUPHAR_to_UniProt_converter[label_key] || label_key).toUpperCase();
                plotGroup.append('text')
                    .attr('x', margin.left + xOffset + labelOffset)
                    .attr('y', yOffset)
                    .attr('class', category)
                    .attr('dy', '-0.3em') // Adjust this value to move the text higher
                    .style('dominant-baseline', 'middle') // Set vertical alignment
                    .style('font-weight', bold ? 'bold' : 'normal')
                    .style('font-style', italic ? 'italic' : 'normal')
                    .style('text-decoration', underline ? 'underline' : 'none')
                    .style('font-size', fontSize)
                    .style('fill', color)
                    .text(label);
            } else if (label_names === 'Gene') {
                // Handle UniProt receptor labels
                label = label_conversion_dicts.IUPHAR_to_Gene_converter[label_key] || label_key;
                plotGroup.append('text')
                    .attr('x', margin.left + xOffset + labelOffset)
                    .attr('y', yOffset)
                    .attr('class', category)
                    .attr('dy', '-0.3em') // Adjust this value to move the text higher
                    .style('dominant-baseline', 'middle') // Set vertical alignment
                    .style('font-weight', bold ? 'bold' : 'normal')
                    .style('font-style', italic ? 'italic' : 'normal')
                    .style('text-decoration', underline ? 'underline' : 'none')
                    .style('font-size', fontSize)
                    .style('fill', color)
                    .text(label);
            }

        } else if (category === 'ReceptorFamily') {
            // Handle ReceptorFamily (subscript handling but no label_names logic)
            label = label_key;

            // Create text element for ReceptorFamily
            const textElement = plotGroup.append('text')
                .attr('x', margin.left + xOffset)
                .attr('y', yOffset)
                .attr('class', category)
                .style('dominant-baseline', 'middle') // Set vertical alignment
                .style('font-weight', bold ? 'bold' : 'normal')
                .style('font-style', italic ? 'italic' : 'normal')
                .style('text-decoration', underline ? 'underline' : 'none')
                .style('font-size', fontSize)
                .style('fill', color);

            // Subscript handling for ReceptorFamily
            const mainFontSize = parseFloat(fontSize);
            const subFontSize = mainFontSize * 0.75 + 'px';
            const parts = label.split(/(<sub>|<\/sub>)/);
            let inSub = false;

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
        } else {
            // General label handling for other categories
            if (category === 'Class') {
                label = label_key.split(" (")[0]; // Trim Class
            }

            // Render text for non-receptor categories
            plotGroup.append('text')
                .attr('x', margin.left + xOffset)
                .attr('y', yOffset)
                .attr('class', category)
                .style('dominant-baseline', 'middle') // Set vertical alignment
                .style('font-weight', bold ? 'bold' : 'normal')
                .style('font-style', italic ? 'italic' : 'normal')
                .style('text-decoration', underline ? 'underline' : 'none')
                .style('font-size', fontSize)
                .style('fill', color)
                .text(label);
        }
    }


    // ###########################
    // ## Render List of Labels ##
    // ###########################

    let label_counter = 1;
    let current_col = 1; // Track the current column being rendered

    for (let i = 0; i < data.length; i++) {
        const label = data[i];
        const category = category_data[i];
        const label_key = label; // Use original label before conversions
        const style = styling_option[category]; // Get the style options for the category

        // Add text for the label
        add_text(label, yOffset, category, label_key, style.Bold, style.Italic, style.Underline, style.Fontsize, style.Color,label_names);

        yOffset += y_off_set_variable; // Adjust Y-offset for the next label

        // Check if this is the new maximum yOffset encountered
        if (yOffset > yOffset_max) {
            yOffset_max = yOffset; // Update yOffset_max globally
        }

        label_counter++;

        // Handle column breaks when Col_break_number is reached
        if (label_counter > Col_break_number && current_col < columns) {
            current_col++;
            xOffset += spacing_dict[`Col${current_col - 1}`]; // Move to the next column
            label_counter = 1; // Reset label counter for the new column
            yOffset = margin.top + 5 + _legendTopSpace; // Reset to same start as column 1
        }

        // Ensure the global maximum yOffset is tracked
        if (yOffset > yOffset_max) {
            yOffset_max = yOffset; // Update the maximum yOffset encountered
        }
    }

    // Rerender height of plot based on the max yOffset encountered
    svg.attr("height", yOffset_max + margin.top + margin.bottom);

    return svg;
}

// Handles the data visualization
function data_visualization(data, category_data, location, Layout_dict, data_styling, spacing_dict,y_off_set_variable=30,data_size=6,data_fontsize_variable=14,BarLabels={}) {

    // ###########################
    // ## Initialize Variables  ##
    // ###########################

    // Set default margins, xOffset, yOffset, and columns based on Layout_dict
    const margin = { top: 40, right: 20, bottom: 20, left: 20 };
    const _vizLegendTopSpace = (Layout_dict && Layout_dict.legend_top_space != null)
        ? Layout_dict.legend_top_space : 50;
    let yOffset = margin.top + 5 + _vizLegendTopSpace + 5;
    let xOffset = 5;
    let columns = Layout_dict.columns;

    // Get the SVG container
    const svg = d3.select(`#${location} svg`);

    // select correct group for appending
    const plotGroup = svg.select(".main-plot-group");

    // Track the current column and Y-offsets
    let current_col = 1;
    let label_counter = 1;
    const Col_break_number = Layout_dict.Col_break_number || Math.ceil(data.length / columns);

    // Function to add different shapes
    function addShape(shapeType, x, y, size, fillColor) {
        switch (shapeType) {
            case 'circle':
                plotGroup.append('circle')
                    .attr('cx', x)
                    .attr('cy', y)
                    .attr('r', size)
                    .style('dominant-baseline', 'middle') // Set vertical alignment
                    .style('stroke', 'black')
                    .style('stroke-width', 1)
                    .style('fill', fillColor);
                break;
            case 'rectangle':
                plotGroup.append('rect')
                    .attr('x', x - size)
                    .attr('y', y - size)
                    .attr('width', size * 2)
                    .attr('height', size * 2)
                    .style('dominant-baseline', 'middle') // Set vertical alignment
                    .style('stroke', 'black')
                    .style('stroke-width', 1)
                    .style('fill', fillColor);
                break;
            case 'triangle':
                plotGroup.append('path')
                    .attr('d', `M ${x} ${y - size} L ${x - size} ${y + size} L ${x + size} ${y + size} Z`)
                    .style('stroke', 'black')
                    .style('dominant-baseline', 'middle') // Set vertical alignment
                    .style('stroke-width', 1)
                    .style('fill', fillColor);
                break;
            case 'diamond':
                plotGroup.append('path')
                    .attr('d', `M ${x} ${y - size} L ${x - size} ${y} L ${x} ${y + size} L ${x + size} ${y} Z`)
                    .style('stroke', 'black')
                    .style('dominant-baseline', 'middle') // Set vertical alignment
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
                plotGroup.append('path')
                    .attr('d', starPath)
                    .style('stroke', 'black')
                    .style('dominant-baseline', 'middle') // Set vertical alignment
                    .style('stroke-width', 1)
                    .style('fill', fillColor);
                break;
            default:
                break;
        }
    }

    // Function to handle column breaks
    function handleColumnBreak() {
        if (label_counter > Col_break_number && current_col < columns) {
            current_col++;
            xOffset += spacing_dict[`Col${current_col - 1}`]; // Move to the next column
            label_counter = 1; // Reset label counter for the new column
            yOffset = margin.top + 5 + _vizLegendTopSpace + 5; // Reset to same start as column 1
        }
    }

    // ###############################
    // ## Color Logic Function      ##
    // ###############################

    function getShapeColor(column, data_value, receptorData = {}) {
        const column_styling = data_styling[column];
        let color = 'black';

        if (column_styling.Datatype === 'Discrete') {
            // Use ColorValue from receptorData
            if (receptorData && receptorData.ColorValue) {
                color = receptorData.ColorValue;
            } else {
                color = 'black';
            }
        } else if (column_styling.Datatype === 'Continuous') {
            const gradientScale = d3.scale.linear()
                .domain([column_styling.Data_min, column_styling.Data_max]);

            if (column_styling.data_color_complexity === 'One') {
                gradientScale.range(['#FFFFFF', column_styling.Data_color2]);
            } else if (column_styling.data_color_complexity === 'Two') {
                gradientScale.range([column_styling.Data_color1, column_styling.Data_color2]);
            } else if (column_styling.data_color_complexity === 'Three') {
                gradientScale.range([column_styling.Data_color1, '#FFFFFF', column_styling.Data_color2])
                    .domain([column_styling.Data_min, (column_styling.Data_min + column_styling.Data_max) / 2, column_styling.Data_max]);
            }

            color = gradientScale(data_value);
        }

        return color;
    }

    // ###############################
    // ## Iterate through data rows ##
    // ###############################

    for (let i = 0; i < data.length; i++) {
        const label = data[i];
        const category = category_data[i];

        // Check if category is 'Receptor' (only for receptors do we render shapes)
        if (category === 'Receptor') {
            // Check if data exists for the receptor (in the shape data)
            if (list_data_wow.hasOwnProperty(label)) {
                const receptorData = list_data_wow[label];

                // Col data checker
                const Col1_data_checker = data_styling.Col1.Data == "Yes";
                const Col2_data_checker = data_styling.Col2.Data == "Yes";
                const Col3_data_checker = data_styling.Col3.Data == "Yes";
                const Col4_data_checker = data_styling.Col4.Data == "Yes";

                // Col labels and shapes
                const Col1_shape = receptorData.hasOwnProperty('Value5') ? receptorData.Value5.toLowerCase() : false;
                const Col2_shape = receptorData.hasOwnProperty('Value6') ? receptorData.Value6.toLowerCase() : false;
                const Col3_shape = receptorData.hasOwnProperty('Value7') ? receptorData.Value7.toLowerCase() : false;
                const Col4_shape = receptorData.hasOwnProperty('Value8') ? receptorData.Value8.toLowerCase() : false;

                // Col data
                const Col1_data = receptorData.hasOwnProperty('Value1') ? receptorData.Value1 : false;
                const Col2_data = receptorData.hasOwnProperty('Value2') ? receptorData.Value2 : false;
                const Col3_data = receptorData.hasOwnProperty('Value3') ? receptorData.Value3 : false;
                const Col4_data = receptorData.hasOwnProperty('Value4') ? receptorData.Value4 : false;

                // Shapes and data rendering for each column
                const { offsetMap } = computeDynamicOffsets(data_styling);
                // const col1_XoffSet = 0;
                // const col2_XoffSet = 20;
                // const col3_XoffSet = 40;
                // const col4_XoffSet = 60;

                const Shape_list = ['circle', 'rectangle', 'triangle', 'star', 'diamond'];

                // ### Column 1 ###
                if (Col1_data_checker && (Col1_shape || Col1_data) && offsetMap['Col1'] !== undefined) {
                    const shape_color = getShapeColor('Col1', Col1_data, receptorData);
                    addShape(Shape_list.includes(Col1_shape) ? Col1_shape : 'circle', margin.left + xOffset + offsetMap['Col1'], yOffset - 10, data_size, shape_color);
                }

                // ### Column 2 ###
                if (Col2_data_checker && (Col2_shape || Col2_data)  && offsetMap['Col2'] !== undefined) {
                    const shape_color = getShapeColor('Col2', Col2_data, receptorData);
                    addShape(Shape_list.includes(Col2_shape) ? Col2_shape : 'circle', margin.left + xOffset + offsetMap['Col2'], yOffset - 10, data_size, shape_color);
                }

                // ### Column 3 ###
                if (Col3_data_checker && (Col3_shape || Col3_data)  && offsetMap['Col3'] !== undefined) {
                    const shape_color = getShapeColor('Col3', Col3_data, receptorData);
                    addShape(Shape_list.includes(Col3_shape) ? Col3_shape : 'circle', margin.left + xOffset + offsetMap['Col3'], yOffset - 10, data_size, shape_color);
                }

                // ### Column 4 ###
                if (Col4_data_checker && (Col4_shape || Col4_data)  && offsetMap['Col4'] !== undefined) {
                    const shape_color = getShapeColor('Col4', Col4_data, receptorData);
                    addShape(Shape_list.includes(Col4_shape) ? Col4_shape : 'circle', margin.left + xOffset + offsetMap['Col4'], yOffset - 10, data_size, shape_color);
                }
            }

            // Increment yOffset for the next label and shape
            yOffset += y_off_set_variable;
            label_counter++;

            // Handle column break
            handleColumnBreak();

        } else {
            // If the category is not 'Receptor', simply move the Y offset without rendering shapes
            yOffset += y_off_set_variable;
            label_counter++;

            // Handle column break
            handleColumnBreak();
        }
    }

    // #############################
    // ## Render Legend Bars on Top ##
    // #############################

    // === Helper functions ====

    function getDecimals(num) {
        const parts = num.toString().split('.');
        return parts[1]?.length || 0;
    }

    function formatMidpoint(low, high) {
        const decimals = Math.max(getDecimals(low), getDecimals(high));
        return ((low + high) / 2).toFixed(decimals);
    }

    // Allow the list mapper (or any caller) to skip these built-in bars and use createListLegendBars instead
    if (Layout_dict && Layout_dict.skip_gradient_bars) return;

    const _barsAtBottom = Layout_dict && Layout_dict.legend_position === 'Bottom';
    // Compact bar dimensions (similar visual weight to the tree's legend bars)
    const _LW = 120; const _BH = 12; const _SB = 15;
    const _XO = 25;  const _TO = 18; const _DF = Math.min(data_fontsize_variable, 11);
    // When bars go at bottom, place them after the last drawn row
    const _barBaseY = _barsAtBottom ? (yOffset_max + 22) : 10;

    let bar_index = 0;
    Object.keys(data_styling).forEach(function(column) {
        if (data_styling[column].Data === "Yes" && data_styling[column].Datatype === 'Continuous') {
            const legendWidth = _LW;
            const data_fontsize = _DF;
            const lowest_value = data_styling[column].Data_min;
            const highest_value = data_styling[column].Data_max;
            const midpointText = formatMidpoint(lowest_value, highest_value);
            const spacing_bar = _SB;
            const bar_height = _BH;
            const text_off_set = _TO;
            const x_off_set = _XO;
            const y_off_set = _barBaseY;

            // Calculate the x position for the current bar
            const x_position = bar_index * (legendWidth + spacing_bar) + x_off_set;

            // Increment bar index for the next bar
            bar_index++;

            const legend_svg = svg.append("g"); // Append group for the legend bar
            const gradientId = `Gradient_${column}`;
            const defs = svg.append('defs');
            var BarText = (BarLabels && BarLabels[column])
                ? BarLabels[column]
                : `Dataset ${column.substr(3)}`;
            // Add centered text on top of the bar specifying the data column
            legend_svg.append("text")
                .attr('x', x_position + legendWidth / 2) // Center the text horizontally over the bar
                .attr('y', bar_height - 10 + y_off_set) // Adjust Y position as needed, e.g., 10px above the bar
                .style("font-size", `${data_fontsize}px`)
                .style("font-family", "sans-serif")
                .style("text-anchor", "middle") // Center the text horizontally
                .style("font-weight", "bold")   // Make the text bold
                .attr('id', `legend_${column.substr(3)}`)
                .text(BarText); // Use the column name as the text

            // Gradient definition
            const gradient = defs.append('linearGradient')
                .attr('id', gradientId)
                .attr('x1', '0%')
                .attr('x2', '100%')
                .attr('y1', '0%')
                .attr('y2', '0%');

            // Adjust gradient stops based on color complexity
            if (data_styling[column].data_color_complexity === 'One') {
                // One color: white to Data_color1
                gradient.append('stop')
                    .attr('offset', '0%')
                    .attr('stop-color', '#FFFFFF'); // Start with white

                gradient.append('stop')
                    .attr('offset', '100%')
                    .attr('stop-color', data_styling[column].Data_color2); // End with color1
            } else if (data_styling[column].data_color_complexity === 'Two') {
                // Two colors: Data_color1 to Data_color2
                gradient.append('stop')
                    .attr('offset', '0%')
                    .attr('stop-color', data_styling[column].Data_color1); // Start with color1

                gradient.append('stop')
                    .attr('offset', '100%')
                    .attr('stop-color', data_styling[column].Data_color2); // End with color2
            } else if (data_styling[column].data_color_complexity === 'Three') {
                // Three colors: Data_color1 to white in the middle, then to Data_color2
                gradient.append('stop')
                    .attr('offset', '0%')
                    .attr('stop-color', data_styling[column].Data_color1); // Start with color1

                gradient.append('stop')
                    .attr('offset', '50%')
                    .attr('stop-color', '#FFFFFF'); // Middle with white

                gradient.append('stop')
                    .attr('offset', '100%')
                    .attr('stop-color', data_styling[column].Data_color2); // End with color2
            }

            // Border around the rectangle
            legend_svg.append('rect')
                .attr('x', x_position) // Align with the gradient rectangle
                .attr('y', bar_height + y_off_set) // Align with the gradient rectangle
                .attr('width', legendWidth) // Same width as the gradient rectangle
                .attr('height', bar_height) // Same height as the gradient rectangle
                .style('fill', 'none')
                .style('stroke', 'black')
                .style('stroke-width', 2);

            // Rectangle for the gradient
            legend_svg.append('rect')
                .attr('x', x_position) // Horizontal spacing between bars
                .attr('y', bar_height + y_off_set) // Adjust Y position as needed
                .attr('width', legendWidth) // Legend bar width
                .attr('height', bar_height) // Legend bar height
                .style('fill', `url(#${gradientId})`);

            // Minimum value text
            legend_svg.append("text")
                .attr('x', x_position) // Align with the bar
                .attr('y', bar_height + text_off_set + y_off_set)
                .style("font-size", `${data_fontsize}px`)
                .style("font-family", "sans-serif")
                .text(lowest_value);

            // Middle value text (if applicable)
            if (data_styling[column].data_color_complexity === 'Three') {
                legend_svg.append("text")
                    .attr('x', x_position + legendWidth / 2)
                    .attr('y', bar_height + text_off_set + y_off_set)
                    .style("font-size", `${data_fontsize}px`)
                    .style("font-family", "sans-serif")
                    .style("text-anchor", "middle")
                    .text(midpointText); // Middle value
            }

            // Maximum value text
            legend_svg.append("text")
                .attr('x', x_position + legendWidth)
                .attr('y', bar_height + text_off_set + y_off_set)
                .style("font-size", `${data_fontsize}px`)
                .style("font-family", "sans-serif")
                .style("text-anchor", "end")
                .text(highest_value);
        }
    });

    // When bars are at the bottom, extend SVG height to fit them
    if (_barsAtBottom && bar_index > 0) {
        const barsBottom = _barBaseY + _BH + _TO + _DF + 10;
        const curH = +svg.attr('height') || 0;
        if (barsBottom > curH) svg.attr('height', barsBottom);
    }
}

function CreateTextLegend_list(location, data, Layout) {
    const layoutMode = Layout?.layoutMode || "row";
    const columns = Layout?.columns || 2;
    const sortDirection = Layout?.sortDirection || "Vertically";
    const TreeLegendPosition = Layout?.TreeLegendPosition || "Top";
    const LegendFontSize = Layout?.Fontsize || "11px";
    const LegendFontStyle = "Arial";
    const spacingY = 25;
    const padding = 10;

    const svg = d3.select('#' + location + ' svg');

    // Remove existing legend
    svg.selectAll(".legend-group").remove();

    const legendGroup = svg.append("g")
        .attr("class", "legend-group");

    const tempText = svg.append("text")
        .attr("x", -9999)
        .attr("y", -9999)
        .style("font-size", LegendFontSize)
        .style("font-family", LegendFontStyle);

    const svgNode = svg.node();
    let svgWidth = 800;
    let svgHeight = 800;

    if (svgNode) {
        const vb = svg.attr("viewBox");
        if (vb) {
            const [, , w, h] = vb.split(" ").map(Number);
            svgWidth = w;
            svgHeight = h;
        } else {
            const box = svgNode.getBoundingClientRect();
            svgWidth = box.width;
            svgHeight = box.height;
        }
    }

    const legendItems = new Set();

    Object.entries(data).forEach(([key, value]) => {
        if (value.Value1 && value.ColorValue) {
            const pairKey = `${value.ColorValue}|||${value.Value1}`;
            legendItems.add(pairKey);
        }
    });

    const sortedItems = Array.from(legendItems)
        .filter(pair => {
            const [color, label] = pair.split("|||");
            return !(color === "#FFFFFF" && label === "Empty");
        })
        .sort((a, b) => {
            const labelA = a.split("|||")[1];
            const labelB = b.split("|||")[1];
            return labelA.localeCompare(labelB, undefined, { numeric: true, sensitivity: 'base' });
        });

    const fontSizeNum = parseFloat(LegendFontSize) || 11;
    const circleRadius = Math.round(fontSizeNum * 0.4);

    let legendElements = [];

    if (layoutMode === "row") {
        const startX = 50;
        let x = startX;
        let y = 40;
        const maxItemWidth = svgWidth - 100;

        sortedItems.forEach(pairKey => {
            let [color, label] = pairKey.split("|||");

            tempText.text(label);
            const labelWidth = tempText.node().getComputedTextLength();
            const totalWidth = 6 * 2 + padding + labelWidth + 20;

            if (totalWidth > maxItemWidth) {
                label = "⚠ Too long label";
                tempText.text(label);
            }

            const fixedWidth = 6 * 2 + padding + tempText.node().getComputedTextLength() + 20;

            if (x + fixedWidth > svgWidth - 50) {
                x = startX;
                y += spacingY;
            }

            const centerY = y + fontSizeNum * 0.35;

            legendElements.push({
                x, y: centerY, color, label,
                labelX: x + circleRadius + 6,
                labelY: centerY + 0.35 * fontSizeNum
            });

            x += fixedWidth;
        });
    } else if (layoutMode === "columns") {
        const colData = Array.from({ length: columns }, () => []);

        if (sortDirection === "Horizontally") {
            sortedItems.forEach((item, i) => colData[i % columns].push(item));
        } else {
            const perCol = Math.floor(sortedItems.length / columns);
            const remainder = sortedItems.length % columns;
            let index = 0;
            for (let i = 0; i < columns; i++) {
                const count = perCol + (i < remainder ? 1 : 0);
                colData[i] = sortedItems.slice(index, index + count);
                index += count;
            }
        }

        const colWidths = colData.map(col => {
            let maxWidth = 0;
            col.forEach(pair => {
                const label = pair.split("|||")[1];
                tempText.text(label);
                maxWidth = Math.max(maxWidth, tempText.node().getComputedTextLength());
            });
            return maxWidth + 40;
        });

        let x = 50;
        colData.forEach((col, colIndex) => {
            let y = 40;
            col.forEach(pairKey => {
                const [color, label] = pairKey.split("|||");
                const centerY = y + fontSizeNum * 0.35;

                legendElements.push({
                    x, y: centerY, color, label,
                    labelX: x + circleRadius + 6,
                    labelY: centerY + 0.35 * fontSizeNum
                });

                y += spacingY;
            });
            x += colWidths[colIndex];
        });
    }

    tempText.remove();

    legendElements.forEach(d => {
        legendGroup.append("circle")
            .attr("cx", d.x)
            .attr("cy", d.y)
            .attr("r", circleRadius)
            .style("fill", d.color)
            .style("stroke", "black");

        legendGroup.append("text")
            .attr("x", d.labelX)
            .attr("y", d.labelY)
            .attr("text-anchor", "start")
            .style("font-size", LegendFontSize)
            .style("font-family", LegendFontStyle)
            .text(d.label);
    });

    // === Positioning + Offsetting ===
    const legendBBox = legendGroup.node()?.getBBox();

    if (legendBBox) {
        // Center horizontally
        const centerOffsetX = (svgWidth - legendBBox.width) / 2 - legendBBox.x;

        let yOffset = 0;
        const paddingBelowLegend = 40;
        const visualHeight = d3.max(legendElements.map(d => d.labelY)) - d3.min(legendElements.map(d => d.y));
        const extraHeight = visualHeight + paddingBelowLegend;

        // Adjust height or move plot based on legend position
        if (TreeLegendPosition === "Bottom") {
            const currentHeight = +svg.attr("height") || 0;
            const newHeight = Math.max(currentHeight, svgHeight + extraHeight);
            svg.attr("height", newHeight);
            yOffset = newHeight - legendBBox.height - paddingBelowLegend*3;
        } else {
            // Shift main plot group down
            const plotGroup = svg.select(".main-plot-group");
            const currentTransform = plotGroup.attr("transform");
            const match = currentTransform?.match(/translate\(([^,]+),([^)]+)\)/);
            if (match) {
                const currentX = parseFloat(match[1]);
                const currentY = parseFloat(match[2]);
                plotGroup.attr("transform", `translate(${currentX}, ${currentY + extraHeight})`);
            } else {
                plotGroup.attr("transform", `translate(0, ${extraHeight-paddingBelowLegend})`);
            }
        }

        // Position legend
        legendGroup.attr("transform", `translate(${centerOffsetX-20}, ${yOffset})`);
    }
}


// #################
// ###  HEATMAP  ###
// #################

// Initialize the datastyling dict
function heatmap_DataStyling() {
    var heatmap_DataStyling = {Number_of_colors: "Three",
        min_color: "#1a80bb",
        middle_color: "#FFFFFF",
        max_color: "#a00000",
        rotation: 90,
        label_position: 'Bottom',
        label_fontsize: 12,
        receptor_fontsize: 12,
        datalabels: true,
        data_border: true,
        data_fontsize: 12,
        legend_label: "Value intensity",
        LabelType: 'IUPHAR'

    }
    return heatmap_DataStyling
}

// Function to handle row labels based on label type
function handleRowLabels(textElement, label, labelType, fontSize) {
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

    // Handle label transformations based on labelType
    let transformedLabel = label; // Initialize transformedLabel with the original label

    if (labelType === 'UniProt') {
        transformedLabel = label.toUpperCase();
        textElement.text(transformedLabel);
    } else if (labelType === 'Gene') {
        transformedLabel = label_converter.UniProt_to_Gene_converter[label];
        textElement.text(transformedLabel);
    } else if (labelType === 'IUPHAR') {

        // Retrieve the IUPHAR label from your converter object
        transformedLabel = label_converter.UniProt_to_IUPHAR_converter[label];
        // Replace HTML entities with corresponding Unicode characters
        for (const [entity, char] of Object.entries(htmlEntities)) {
            transformedLabel = transformedLabel.replace(new RegExp(entity, 'g'), char);
        }
        // Remove <i> and </i> tags from the label
        transformedLabel = transformedLabel.replace(/<\/?i>/g, '');
        // remove receptor and adrenoceptor
        transformedLabel = transformedLabel.replace("-adrenoceptor", '').replace(" receptor", '');
        // Clear existing content in textElement if needed
        textElement.text('')

        // Additional transformations specific to IUPHAR labels
        const parts = transformedLabel.split(/(<sub>|<\/sub>)/); // Split label into parts including <sub> tags
        let inSub = false; // Flag to track if we are inside a subscript

        // Calculate the subscript font size (e.g., 75% of the main font size)
        const mainFontSize = parseFloat(fontSize);
        const subFontSize = parseInt(mainFontSize * 0.75,10);

        // Handle each part of the label
        parts.forEach(part => {
            if (part === '<sub>') {
                inSub = true;
            } else if (part === '</sub>') {
                inSub = false;
            } else {
                const tspan = textElement.append('tspan').text(part).style("font-family", "sans-serif");

                // Adjust dy based on subscript flag
                if (inSub) {
                    tspan.attr('dy', '0.5em').attr('font-size', subFontSize); // Adjust dy for subscript
                } else {
                    tspan.attr('dy', '0.1em').attr('font-size', mainFontSize); // Default dy for regular text
                }
            }
        });
    }
}

// Create the heatmap
function Heatmap(data, location, heatmap_DataStyling,label_x_converter) {

    // Format a numeric value for display in cells and the colour-bar legend.
    // Shows full precision for numbers whose digit characters (0-9, not counting '-', '.', ',')
    // number ≤ 10; switches to scientific notation (3 d.p.) for longer numbers.
    function fmtHeatmapVal(val) {
        const n = parseFloat(val);
        if (!isFinite(n)) return '';
        const s = n.toString();
        if ((s.match(/\d/g) || []).length <= 12) return s;
        return n.toExponential(3); // e.g. 1.235e+10
    }

    const margin = { top: 30, right: 100, bottom: 30, left: 60 }; // Adjusted margin for row labels
    const rows = Object.keys(data);

    // Collect all column keys that appear in at least one row
    const colSet = new Set();
    rows.forEach(row => {
        Object.keys(data[row] || {}).forEach(col => colSet.add(col));
    });

    // Use a fixed canonical order for Value1-Value5 so that skipped columns
    // don't reorder the remaining ones. Any extra keys are appended after.
    const FIXED_COL_ORDER = ['Value1', 'Value2', 'Value3', 'Value4', 'Value5'];
    const cols = [
        ...FIXED_COL_ORDER.filter(c => colSet.has(c)),
        ...Array.from(colSet).filter(c => !FIXED_COL_ORDER.includes(c))
    ];
    // Blank cells (rows that have no value for a column) are handled below:
    // data[row][col] === undefined → isNaN → fill:'None' (transparent cell).
    const col_labels = cols.map((col, i) => {
        const label = label_x_converter[col];
        return (typeof label === 'string' && label.trim() !== '') ? label : `Dataset ${i + 1}`;
      });

    const rotation = heatmap_DataStyling.rotation;
    const data_labels = heatmap_DataStyling.datalabels;
    const labelType = heatmap_DataStyling.LabelType;
    const label_position = heatmap_DataStyling.label_position;
    const label_fontsize = heatmap_DataStyling.label_fontsize;
    const receptor_fontsize = heatmap_DataStyling.receptor_fontsize;
    const data_fontsize = heatmap_DataStyling.data_fontsize;
    let legend_y_position = 30;

    let rowLabelWidth;
    const baseWidth = 20;  // Minimum column width
    const chartData = [];
    let highest_value = -Infinity;
    let lowest_value = Infinity;

    rows.forEach(row => {
      cols.forEach(col => {
        const value = data[row][col];
        chartData.push({ row, col, value });
        if (value > highest_value) {
          highest_value = value;
        }
        if (value < lowest_value) {
            lowest_value = value;
        }
      });
    });

    // Calculate the length of the longest data value

    function measureTextSize(text, fontSize, fontFamily = 'sans-serif') {
        // Create temporary SVG text element
        const tempSvg = d3.select('body').append('svg')
        .attr('class', 'temp-text-measurer')
        .style('position', 'absolute')
        .style('visibility', 'hidden');

        const tempText = tempSvg.append('text')
        .style('font-size', `${fontSize}px`)
        .style('font-family', fontFamily)
        .text(text);

        const bbox = tempText.node().getBBox();
        tempSvg.remove();  // Clean up

        return { width: bbox.width, height: bbox.height };
    }

    // Measure the longest column label
    let longestLabelSize = { width: 0, height: 0 };

    col_labels.forEach(label => {
    const size = measureTextSize(label, label_fontsize);
    if (size.width > longestLabelSize.width) longestLabelSize.width = size.width;
    if (size.height > longestLabelSize.height) longestLabelSize.height = size.height;
    });

    // Measure every unique formatted data value's width using measureTextSize
    // (same mechanism that successfully measures longestLabelSize above).
    // Character count ≠ pixel width, so we measure each unique string individually.
    const cellPad = 0; // px padding each side inside the cell
    let maxDataValueWidth = 0;
    const seenValStrs = new Set();
    chartData.forEach(function (d) {
        if (!isNaN(d.value)) {
            seenValStrs.add(fmtHeatmapVal(d.value));
        }
    });
    seenValStrs.forEach(function (valStr) {
        const sz = measureTextSize(valStr, data_fontsize);
        if (sz.width > maxDataValueWidth) maxDataValueWidth = sz.width;
    });

    // Column width: all columns equal, driven by the widest measurement.
    // Both rotation modes use maxDataValueWidth as the primary cell-width driver.
    if (rotation === 90 || rotation === 45) {
        // Vertical labels: cell width must fit the data value + padding;
        // longestLabelSize.height is the font line-height (≈ font size), the label's
        // footprint in the column direction when rotated 90°.
        rowLabelWidth = Math.max(baseWidth,
                                 maxDataValueWidth + cellPad * 2,
                                 longestLabelSize.height + 4);
    } else {
        // Horizontal labels: wider of the column-label text or the data value + padding.
        rowLabelWidth = Math.max(baseWidth,
                                 longestLabelSize.width + 8,
                                 maxDataValueWidth + cellPad * 2);
    }

    // Adjust margin and legend position for top/bottom label placement
    if (label_position === 'Top' && rotation === 90) {
        margin.top = longestLabelSize.width + 20;  // Add padding
        legend_y_position = 0;
    } else if (label_position === 'Bottom' && rotation === 90) {
        legend_y_position = longestLabelSize.width;
    } else if (label_position === 'Top' && rotation === 0) {
        legend_y_position = 0;
    }

    // Measure the longest Y-axis (receptor) label to set left margin dynamically
    let longestRowLabelPx = 0;
    rows.forEach(function(row) {
        var lbl;
        if (labelType === 'UniProt') {
            lbl = row.toUpperCase();
        } else if (labelType === 'Gene') {
            lbl = (window.label_converter && window.label_converter.UniProt_to_Gene_converter &&
                   window.label_converter.UniProt_to_Gene_converter[row]) || row;
        } else {
            lbl = (window.label_converter && window.label_converter.UniProt_to_IUPHAR_converter &&
                   window.label_converter.UniProt_to_IUPHAR_converter[row]) || row;
            lbl = String(lbl).replace(/<\/?(?:sub|i)>/g, '');
        }
        var w = measureTextSize(String(lbl || row), receptor_fontsize).width;
        if (w > longestRowLabelPx) longestRowLabelPx = w;
    });
    const yAxisAreaWidth = Math.max(60, Math.ceil(longestRowLabelPx) + 12);

    const width = (rowLabelWidth * cols.length) + rowLabelWidth;
    const height = (20 * rows.length);
    const adjustedTopMargin = (rotation === 90 && label_position === 'Top') ? longestLabelSize.width + 20 : margin.top;

    const svg_home = d3.select("#" + location)
        .append("svg")
        .attr("width", yAxisAreaWidth + width + margin.right)
        .attr("height", height + adjustedTopMargin * 2)
        .attr("id", "Heatmap_plot_svg");

    //  x / y scale transformers
    const x = d3.scale.ordinal()
        .rangeBands([0, width])
        .domain(cols);

    const y = d3.scale.ordinal()
      .rangeBands([height, 0])
      .domain(rows.reverse());

    if (heatmap_DataStyling.Number_of_colors === 'Three') {
        var myColor = d3.scale.linear()
        .range([heatmap_DataStyling.min_color, heatmap_DataStyling.middle_color, heatmap_DataStyling.max_color]) // Adjust colors here
        .domain([lowest_value, (highest_value + lowest_value) / 2, highest_value]); // Adjust domain based on your data
    } else if (heatmap_DataStyling.Number_of_colors === 'One' || heatmap_DataStyling.Number_of_colors === 'Two') {
        var myColor = d3.scale.linear()
        .range([heatmap_DataStyling.min_color, heatmap_DataStyling.max_color]) // Adjust colors here
        .domain([lowest_value, highest_value]); // Adjust domain based on your data
    }

    const svg = svg_home.append("g")
        .attr("transform", `translate(${yAxisAreaWidth}, ${margin.top})`);

    let xAxis;
    if (label_position === 'Bottom') {
    xAxis = svg.append("g")
        .attr("transform", `translate(0, ${height})`)
        .attr('id', 'Xaxis')
        .call(d3.svg.axis().scale(x).orient("bottom").tickSize(0));
    } else if (label_position === 'Top') {
    xAxis = svg.append("g")
        .attr("transform", `translate(0, -30)`)
        .attr('id', 'Xaxis')
        .call(d3.svg.axis().scale(x).orient("top").tickSize(0));
    }

    xAxis.selectAll("text")
    .attr("class", "column-label")
    .style("font-size", `${label_fontsize}px`)
    .style("font-family", "sans-serif")
    .each(function(d, i) {
        const col = cols[i]; // Assuming cols is defined somewhere in your code
        const text = d3.select(this);

        // Set text content based on label_x_converter
        const label = label_x_converter[col];
        text.text((typeof label === 'string' && label.trim() !== '') ? label : `Dataset ${i + 1}`);

        if (rotation === 0) {
            text.style("text-anchor", "middle")
                .attr("transform", "rotate(0)")
                .attr("dx", `0`)
                .attr("dy", `${label_fontsize * 1.5}px`); // Adjust dy based on font size
        } else if (rotation === 45) {
            text.style("text-anchor", label_position === 'Bottom' ? "end" : "start")
                .attr("transform", "rotate(-45)")
                .attr("dx", label_position === 'Bottom' ? `${-0.5 * label_fontsize}px` : `${-1.1 * label_fontsize}px`)
                .attr("dy", label_position === 'Bottom' ? `${0.5 * label_fontsize}px` : `${1.1 * label_fontsize}px`);
        } else if (rotation === 90) {
            text.style("text-anchor", label_position === 'Bottom' ? "end" : "start")
                .attr("transform", "rotate(-90)")
                .attr("dx", label_position === 'Bottom' ? `${-0.3 * label_fontsize}px` : `${-1.3 * label_fontsize}px`)
                .attr("dy", label_position === 'Bottom' ? `0px` : `${0.6 * label_fontsize}px`);
        }
    });

    // Y labels
    const yAxis = svg.append("g")
        .attr('id', 'Yaxis')
        .call(d3.svg.axis().scale(y).orient("left").tickSize(0))
        .selectAll("text")
        .style("font-size", `${receptor_fontsize}px`)
        .style("font-family", "sans-serif");

    yAxis.each(function(d) {
        const textElement = d3.select(this);

        // Call function to handle row labels based on label type
        handleRowLabels(textElement, d, labelType,receptor_fontsize);
    });

    // Heatmap data boxes

    const rects = svg.selectAll("rect")
      .data(chartData, d => `${d.row}:${d.col}`)
      .enter()
      .append("rect")
      .attr("x", d => x(d.col))
      .attr("y", d => y(d.row))
      .attr("width", x.rangeBand())
      .attr("height", y.rangeBand())
      .style("fill", d => isNaN(d.value) ? 'None' : myColor(d.value))
      .each(function(d) {
        if (heatmap_DataStyling.data_border) {
            d3.select(this)
                .style('stroke', 'black')
                .style('stroke-width', 0.5);
        }
    });


    if (data_labels) {
      rects.each(function(d) {

        if (isNaN(d.value)) return;  // Skip label if value is NaN

        const rect = d3.select(this);
        const textColor = getContrastColor(myColor(d.value));

        svg.append("text")
          .attr("x", +rect.attr("x") + rect.attr("width") / 2)
          .attr("y", +rect.attr("y") + rect.attr("height") / 2)
          .attr("dy", ".35em")
          .attr("text-anchor", "middle")
          .style("font-size", `${data_fontsize}px`)
          .style("font-family", "sans-serif")
          .style("fill", textColor)
          .text(fmtHeatmapVal(d.value));
      });
    }

    function getContrastColor(hexColor) {
      const rgb = d3.rgb(hexColor);
      const brightness = (rgb.r * 0.299 + rgb.g * 0.587 + rgb.b * 0.114);
      return brightness > 128 ? "black" : "white";
    }

    // ##########
    // # Legend #
    // ##########

    // Adjust for the width of the legend to center it
    let legendX;
    if (cols.length % 2 === 0) {
        // Even number of columns
        const leftMiddleColIndex = (cols.length / 2) - 1;
        const rightMiddleColIndex = cols.length / 2;
        const leftMiddleX = x(cols[leftMiddleColIndex]) + x.rangeBand() / 2;
        const rightMiddleX = x(cols[rightMiddleColIndex]) + x.rangeBand() / 2;
        legendX = (leftMiddleX + rightMiddleX) / 2;
    } else {
        // Odd number of columns
        const middleColIndex = Math.floor(cols.length / 2);
        legendX = x(cols[middleColIndex]) + x.rangeBand() / 2; // Center of the middle column
    }

    // Adjust for the width of the legend to center it
    let legendWidth;
    if (width <= 300) {
        legendWidth = width;
    } else {
        legendWidth = 300;
    }
    const adjustedLegendX = legendX - (legendWidth / 2);

    // Position the legend using the adjustedLegendX
    const legend_svg = svg.append("g")
        .attr("transform", `translate(${adjustedLegendX}, ${height + legend_y_position})`);

    const gradientId = "Heatmap_gradient";
    const defs = svg_home.append('defs');

    // Gradient definition
    const gradient = defs.append('linearGradient')
      .attr('id', gradientId)
      .attr('x1', '0%')
      .attr('x2', '100%')
      .attr('y1', '0%')
      .attr('y2', '0%');

    gradient.append('stop')
      .attr('offset', '0%')
      .attr('stop-color', heatmap_DataStyling.min_color); // Start with red at 0%

    if (heatmap_DataStyling.Number_of_colors === 'Three') {
        gradient.append('stop')
        .attr('offset', '50%')
        .attr('stop-color', heatmap_DataStyling.middle_color); // Middle point is white
    }
    gradient.append('stop')
      .attr('offset', '100%')
      .attr('stop-color', heatmap_DataStyling.max_color); // End with green at 100%

    const gradientRect = legend_svg.append('rect')
      .attr('x', 1) // Adjust for border
      .attr('y', 16) // Adjust for border
      .attr('width', legendWidth - 2) // Adjust for border
      .attr('height', 14) // Adjust for border
      .style('fill', `url(#${gradientId})`);

    if (heatmap_DataStyling.data_border) {
      const borderRect = legend_svg.append('rect')
      .attr('x', 0)
      .attr('y', 15)
      .attr('width', legendWidth)
      .attr('height', 15)
      .style('fill', 'none')
      .style('stroke', 'black')
      .style('stroke-width', 1);
    }
    legend_svg.append("text")
      .attr('x', 0)
      .attr('y', 50)
      .style("font-size", `${data_fontsize}px`)
      .style("font-family", "sans-serif")
      .text(fmtHeatmapVal(lowest_value));

    if (heatmap_DataStyling.Number_of_colors === 'Three') {
        legend_svg.append("text")
        .attr('x', legendWidth / 2)
        .attr('y', 50)
        .style("font-size", `${data_fontsize}px`)
        .style("font-family", "sans-serif")
        .style("text-anchor", "middle")
        .text(fmtHeatmapVal((highest_value + lowest_value) / 2));
    }

    legend_svg.append("text")
      .attr('x', legendWidth)
      .attr('y', 50)
      .style("font-size", `${data_fontsize}px`)
      .style("font-family", "sans-serif")
      .style("text-anchor", "end")
      .text(fmtHeatmapVal(highest_value));

    // Rerender height of plot as the last thing using real label height
    const extraPadding = 55;  // base padding
    const labelHeightAdjustment = (rotation === 90) ? longestLabelSize.width + 20 : 30;

    svg_home.attr("height", height + margin.bottom + labelHeightAdjustment + extraPadding);

  }

// =========================================
// ===  Draw / generate the GPCRome plot ===
// =========================================

function DrawGPCRomeWheel(Data, location, GPCRome_styling) {

    dimensions = { height: 1000, width: 1000 };

    const showIcon = GPCRome_styling.showIcon;  // Get the icon visibility state
    const FontsizeGlobal = GPCRome_styling.FontsizeGlobal || "11px";
    const FontsizeClass = GPCRome_styling.FontsizeClass || "20px";
    const FontStyle = GPCRome_styling.Fontstyle || "Arial";
    const dataTypeRaw = GPCRome_styling.DataType != null ? String(GPCRome_styling.DataType) : "";
    const DataType = dataTypeRaw.trim().toLowerCase() === "text" ? "Text" : "Numeric";
    const ColorSetup = GPCRome_styling.ColorSetup || "One";
    const MinValue = GPCRome_styling.GPCRomeMin || 0;
    const MaxValue = GPCRome_styling.GPCRomeMax || 1;
    const AvgValue =  GPCRome_styling.GPCRomeAvg || 0.5;
    const ColorMin =  GPCRome_styling.colorStart || "#FFFFFF";
    const ColorMax =  GPCRome_styling.colorEnd || "#000000";
    const ColorAvg = "#FFFFFF";
    const ShowLegend = GPCRome_styling.ShowLegend || false;
    const LegendLabel = GPCRome_styling.LegendLabel || "";
    const LegendMode = GPCRome_styling.LegendLayout?.mode || "row";
    const LegendCols = GPCRome_styling.LegendLayout?.columns || "1";
    const Legendsorting = GPCRome_styling.LegendLayout?.sorted || "Vertically";

    const svg = d3v4.select("#" + location)
    .append("svg")
    .attr("id", location+'_svg')  // Add the id here for download reference
    .attr("width", dimensions.width)
    .attr("height", dimensions.height)
    .attr("xmlns", "http://www.w3.org/2000/svg")  // Add the SVG namespace
    .attr("xmlns:xlink", "http://www.w3.org/1999/xlink");  // Add the xlink namespace for images

    // Get the image URL from the data-attribute
    const imageUrl = document.getElementById("image-container").getAttribute("data-image-url");

    if (showIcon) {
        var img = new Image();
        img.crossOrigin = "Anonymous";  // Ensure cross-origin handling
        img.src = imageUrl;

        img.onload = function() {
            var canvas = document.createElement('canvas');
            canvas.width = img.width;
            canvas.height = img.height;
            var context = canvas.getContext('2d');
            context.drawImage(img, 0, 0);

            // Convert the image to a base64-encoded string
            var dataUrl = canvas.toDataURL('image/png');

            // Append the image to the SVG using 'xlink:href' (D3 v4 compatible)
            svg.append("image")
                .attr("xlink:href", dataUrl)  // Use 'xlink:href' for D3 v4 compatibility
                .attr("x", 0)  // Top-left corner
                .attr("y", -30)  // Top-left corner
                .attr("width", 210)  // Set width for the image
                .attr("height", 210)  // Set height for the image
                .attr("class", "toggle-image");  // Add a class to control visibility
        };
    }
    // Trailing per-circle arguments to Draw_a_GPCRome, tuned by eye (no formula derives any of
    // these from radius): SeamGapStartDegrees/SeamGapStopDegrees are the plain-degree widths of
    // the gap reserved before/after the 12-o'clock seam around whichever class badge sits there
    // (the only badge on a single-class circle; the first badge in CircleHeaders order on a
    // multi-class one) — independent so, e.g., a side next to a two-line family label can get
    // more room than one next to a plain receptor (start=4, stop=6 spans -4deg to +6deg).
    // NonSeamBadgeUnits sizes every *other* class badge on that same circle (e.g. B2 alongside
    // seam-sitting B1) as that many receptor-tick-widths; unused on single-class circles (1/2,
    // all 4 Odorant), since there's no "other" badge there. See computeClassLayout.
    if (Object.keys(Data).length === 5) {
        Draw_a_GPCRome(Data.Circle_1, 0, 440, dimensions, 4.4, 3.9)
        Draw_a_GPCRome(Data.Circle_2, 1, 355, dimensions, 5.3, 4)
        Draw_a_GPCRome(Data.Circle_3, 2, 265, dimensions, 6.2, 3, 3)
        Draw_a_GPCRome(Data.Circle_4, 3, 185, dimensions, 7.2, 1.7, 3)
        Draw_a_GPCRome(Data.Circle_5, 4, 90, dimensions, 14, 14, 4)
    } else if (Object.keys(Data).length === 4) {
        Draw_a_GPCRome(Data.Circle_1, 0, 440, dimensions, 4.4, 3)
        Draw_a_GPCRome(Data.Circle_2, 1, 345, dimensions, 5.2, 3.5)
        Draw_a_GPCRome(Data.Circle_3, 2, 250, dimensions, 7, 3.5)
        Draw_a_GPCRome(Data.Circle_4, 3, 140, dimensions, 10, 6)
    }

    function Draw_a_GPCRome(Data, level, Radius, dimensions, SeamGapStartDegrees, SeamGapStopDegrees, NonSeamBadgeUnits) {

        // Define SVG dimensions
        const width = dimensions.width;
        const height = dimensions.height;
        const label_offset = 7; // Increased offset to push labels outward
        let GPCRome_radius = Radius || Math.min(width, height) / 2 - 60 - ((level === 4) ? (90 * level) : (85 * level));
        // Width, in radians, of the manually-set gap reserved before/after the 12-o'clock seam
        // around this circle's seam-sitting class badge — independently sized so, e.g., a side
        // that lands next to a two-line family label can get more room than one next to a plain
        // receptor (see the Draw_a_GPCRome call sites for the plain-degrees input and rationale).
        const seamGapStartRad = (SeamGapStartDegrees != null ? SeamGapStartDegrees : 4) * Math.PI / 180;
        const seamGapStopRad = (SeamGapStopDegrees != null ? SeamGapStopDegrees : 4) * Math.PI / 180;
        // How many "arc spaces" (receptor-tick widths) every *other* class badge on this circle
        // (e.g. B2 alongside seam-sitting B1) gets — unused on single-class circles, so the
        // default here never matters for those.
        const nonSeamBadgeUnits = NonSeamBadgeUnits != null ? NonSeamBadgeUnits : 4;

        function extractHeaders(circleData) {
            let CircleHeaders = Object.keys(circleData); // Top-level keys (Classes)
            let CircleSubHeaders = [];

            // Collect all Receptor Families (keys inside each Class), but only if there is more than one
            CircleHeaders.forEach(classKey => {
                let receptorFamilies = Object.keys(circleData[classKey]);

                // Only add receptor families if there are more than one
                if (!["T2", "Unclassified"].some(f => CircleHeaders.includes(f))) {
                    CircleSubHeaders.push(...receptorFamilies);
                }
            });

            // Remove duplicates from CircleSubHeaders
            CircleSubHeaders = [...new Set(CircleSubHeaders)];

            return { CircleHeaders, CircleSubHeaders };
        }

        // Family labels forced onto two lines regardless of the usual length threshold — purely a
        // visual call (e.g. "Calcium-sensing" reads better split even though it's under 18 chars).
        const FORCE_TWO_LINE_FAMILIES = new Set(["Calcium-sensing"]);

        // Builds the ordered list of real content (receptor-family headers + receptors) for this
        // circle, with NO blank filler entries and NO class-header entries — class badges are
        // drawn separately by drawClassBadges() once each class's reserved gap angle is known.
        function buildContentItems(circleData, CircleHeaders, CircleSubHeaders) {
            let contentItems = [];
            let DataFill = {};
            const collator = new Intl.Collator(undefined, { numeric: true, sensitivity: 'base' });

            for (const classKey of CircleHeaders) {

                // Sort receptor families naturally before iterating
                let receptorFamilies = Object.keys(circleData[classKey]).sort(collator.compare);

                // Move "Class A orphans" to the end if it exists
                const orphanIndex = receptorFamilies.indexOf("Class A orphans");

                if (orphanIndex !== -1) {
                    receptorFamilies.splice(orphanIndex, 1);         // remove it
                    receptorFamilies.push("Class A orphans");        // add to end
                }

                // Loop through receptor families
                for (const receptorFamily of receptorFamilies) {
                    if (CircleSubHeaders.includes(receptorFamily)) {
                        // Precomputed with the exact same formatting the family-label rendering
                        // pass uses, so the two checks can never drift apart.
                        const formattedFamilyLabel = GPCRome_formatTextWithHTML(receptorFamily, CircleSubHeaders);
                        contentItems.push({
                            type: 'family',
                            label: receptorFamily,
                            classKey,
                            familyKey: receptorFamily,
                            needsSplit: formattedFamilyLabel.length > 18 || FORCE_TWO_LINE_FAMILIES.has(formattedFamilyLabel)
                        });
                    }

                    // Loop through receptors inside each receptor family
                    for (const receptor of Object.keys(circleData[classKey][receptorFamily])) {

                        let label;

                        if (GPCRome_styling.LabelType === "Uniprot") {
                            label = circleData[classKey][receptorFamily][receptor]["EntryName"] || receptor;
                        } else if (GPCRome_styling.LabelType === "Entrez") {
                            label = circleData[classKey][receptorFamily][receptor]["Entrez"] || receptor;
                        } else {
                            label = receptor; // fallback
                        }

                        contentItems.push({
                            type: 'receptor',
                            label,
                            classKey,
                            familyKey: receptorFamily
                        });

                        const receptorObj = circleData[classKey][receptorFamily][receptor];
                        if (DataType === "Text") {
                            DataFill[label] = receptorObj.Color || "White";
                        } else {
                            DataFill[label] = receptorObj.Data || 0;
                        }
                    }
                }
            }

            return { contentItems, DataFill };
        }

        // Walks CircleHeaders in order, reserving one gap per class boundary (the circle is a
        // closed loop, so that's exactly one gap per class) and laying content out contiguously
        // through whatever's left. Mutates each contentItems entry with startAngle/endAngle/midAngle
        // (all plain "cursor" angles: 0 at 12 o'clock, increasing clockwise — matching d3's arc()
        // convention directly, so these feed arcGenerator() for the pie wedges with no extra
        // sign/offset juggling), plus splitFirstAngle/splitSecondAngle (only meaningful for
        // needsSplit family items) giving two evenly-spaced positions one unit in from each edge of
        // their 3-unit slot. Every item is centered within its own weight*anglePerUnit slot — a
        // family header's bigger weight (2 single-line, 3 needing a split) than a plain receptor's
        // (1) is what gives it breathing room on both sides "for free", with no separate gap concept
        // needed for families at all.
        //
        // The class (badge) gap works differently for exactly one class per circle: whichever class
        // sits at the 12-o'clock seam (the only class, on a single-class circle; the first class in
        // CircleHeaders order, on a multi-class one) gets a manually-set gap — seamGapStartRad before
        // the seam, seamGapStopRad after it, independently, since the two sides often want different
        // amounts of room (e.g. one landing next to a two-line family label, the other next to a
        // plain receptor). Every *other* class on the same circle (e.g. B2 alongside seam-sitting B1)
        // is instead inserted as a plain weighted item into the very same pool as families/receptors,
        // sized as classBadgeUnits "arc spaces" — deliberately: giving *every* class on a shared
        // circle its own hand-picked number would reopen the old "B1 gets special seam treatment, B2
        // doesn't" inconsistency, so only the one seam class gets manual numbers.
        //
        // The cursor starts at -seamGapStartRad so the seam class's gap straddles 12 o'clock (not
        // necessarily symmetrically, if start and stop differ) instead of sitting entirely after it.
        function computeClassLayout(CircleHeaders, contentItems, seamGapStartRad, seamGapStopRad, classBadgeUnits) {
            let contentWeightSum = 0;
            for (const item of contentItems) {
                item.weight = (item.type === 'family') ? (item.needsSplit ? 3 : 2) : 1;
                contentWeightSum += item.weight;
            }

            // Every class except the seam one consumes classBadgeUnits from this same pool — there
            // are (CircleHeaders.length - 1) of those (0 for a single-class circle, so this term
            // vanishes and the only gap is the manual seam one).
            const totalUnits = contentWeightSum + classBadgeUnits * (CircleHeaders.length - 1);
            const seamGapRad = seamGapStartRad + seamGapStopRad;
            const anglePerUnit = totalUnits > 0 ? (2 * Math.PI - seamGapRad) / totalUnits : 0;

            const classGapMidAngle = {};

            let cursor = -seamGapStartRad;
            let itemIndex = 0;
            let isSeamClass = true;

            for (const classKey of CircleHeaders) {
                const gapAngleRad = isSeamClass ? seamGapRad : (classBadgeUnits * anglePerUnit);
                const gapEnd = cursor + gapAngleRad;
                // The seam badge always sits at the standard 12-o'clock position (angle 0),
                // regardless of how asymmetric seamGapStart/StopRad are — start/stop only control
                // how much breathing room content gets on each side, not where the badge itself
                // is drawn. Every other badge keeps the usual "centered in its own gap" position.
                classGapMidAngle[classKey] = isSeamClass ? 0 : (cursor + gapEnd) / 2;
                cursor = gapEnd;
                isSeamClass = false;

                while (itemIndex < contentItems.length && contentItems[itemIndex].classKey === classKey) {
                    const item = contentItems[itemIndex];

                    item.startAngle = cursor;
                    item.endAngle = cursor + item.weight * anglePerUnit;
                    item.midAngle = (item.startAngle + item.endAngle) / 2;
                    item.splitFirstAngle = item.startAngle + anglePerUnit;
                    item.splitSecondAngle = item.endAngle - anglePerUnit;

                    cursor = item.endAngle;
                    itemIndex++;
                }
            }

            return { classGapMidAngle };
        }

        let { CircleHeaders, CircleSubHeaders } = extractHeaders(Data);
        let { contentItems, DataFill } = buildContentItems(Data, CircleHeaders, CircleSubHeaders);
        let { classGapMidAngle } = computeClassLayout(
            CircleHeaders, contentItems, seamGapStartRad, seamGapStopRad, nonSeamBadgeUnits
        );

        // Converts a cursor angle (0 = 12 o'clock, increasing clockwise) into an (x, y, rotation)
        // for a text label or split-label anchor point. Headers/badges are positioned separately
        // by drawClassBadges() and no longer flow through this function.
        function positionFromAngle(cursorAngle, kind, isSplit) {
            const theta = (Math.PI / 2) - cursorAngle;
            const adjustedRadius = isSplit
                ? (GPCRome_radius - 8)
                : (kind === 'family' ? (GPCRome_radius - 20) : (GPCRome_radius + label_offset));

            const x = width / 2 + Math.cos(theta) * adjustedRadius;
            const y = height / 2 - Math.sin(theta) * adjustedRadius;
            const rotation = -(theta * 180 / Math.PI);

            return { x, y, rotation };
        }

        // Which angle a label should be anchored at for pass-1 rendering. Split family labels are
        // re-rendered and removed entirely by pass 2, so their pass-1 position never shows. For a
        // non-split item, inset-by-one-anglePerUnit from either edge of its 2-unit slot (the same
        // convention split lines use on their 3-unit slot's splitFirstAngle/splitSecondAngle) lands
        // on the exact same point both ways — the midpoint — which is why plain d.midAngle is
        // already the correct anchor here.
        function labelAnchorAngle(d) {
            return d.midAngle;
        }

        function getLabelText(d) {
            if (typeof d === 'object' && d !== null) {
                return d.name || d.label || '';
            } else {
                return d;
            }
        }

        function GPCRome_formatTextWithHTML(text, Family_list) {
            // Define a dictionary of HTML entity replacements
            const htmlEntities = {
                "&alpha;": "α",
                "&beta;": "β",
                "&gamma;": "γ",
                "&delta;": "δ",
                "&epsilon;": "ε",
                "&zeta;": "ζ",
                "&eta;": "η",
                "&theta;": "θ",
                "&iota;": "ι",
                "&kappa;": "κ",
                "&lambda;": "λ",
                "&mu;": "μ",
                "&nu;": "ν",
                "&xi;": "ξ",
                "&omicron;": "ο",
                "&pi;": "π",
                "&rho;": "ρ",
                "&sigma;": "σ",
                "&tau;": "τ",
                "&upsilon;": "υ",
                "&phi;": "φ",
                "&chi;": "χ",
                "&psi;": "ψ",
                "&omega;": "ω",
                "&ndash;": "-",  // En dash to hyphen
                "&mdash;": "--", // Em dash to double hyphen
                "&nbsp;": " ",   // Non-breaking space to regular space
                "&lt;": "<",
                "&gt;": ">",
                "&amp;": "&",
                "&quot;": '"',
                "&apos;": "'"
            };

            // Apply all the replacements step by step
            let formattedText = text
                .replace(/ receptors/g, '')
                .replace(/ receptor/g, '')
                .replace(/-adrenoceptor/g, '')
                .replace(/ receptor-/g, '-')
                .replace(/<sub>/g, '</tspan><tspan baseline-shift="-20%">')
                .replace(/<\/sub>/g, '</tspan><tspan>')
                .replace(/<i>/g, '</tspan><tspan font-style="italic">')
                .replace(/<\/i>/g, '</tspan><tspan>')
                .replace(/Long-wave-sensitive/g, 'LWS')
                .replace(/Medium-wave-sensitive/g, 'MWS')
                .replace(/Short-wave-sensitive/g, 'SWS')
                .replace(/Olfactory/g, 'OLF')
                .replace(/calcitonin-like receptor/g, 'CLR')
                .replace(/5-Hydroxytryptamine/g, '5-HT');

            // Replace HTML entities
            formattedText = formattedText.replace(/&[a-z]+;/g, match => htmlEntities[match] || match);

            // Capitalize the first letter only if it's not a Greek letter, special entity, or already has uppercase letters
            function capitalizeFirstLetter(str) {
                // If there's any uppercase letter, skip capitalization
                if (/[A-Z]/.test(str)) return str;

                // If the string contains any known HTML entity, skip capitalization
                for (let key in htmlEntities) {
                    if (str.includes(htmlEntities[key])) return str;
                }

                // Match the first actual letter (skip non-letters at the beginning)
                let match = str.match(/^[^a-zA-Z]*([a-zA-Z])/);
                if (match) {
                    let firstLetter = match[1];
                    let index = match.index + match[0].length - 1;
                    return str.slice(0, index) + firstLetter.toUpperCase() + str.slice(index + 1);
                }

                return str;
            }

            // Apply capitalization only if needed
            formattedText = capitalizeFirstLetter(formattedText);

            // Check if the text is in the Family_list for additional formatting
            const isInFamilyList = Family_list.includes(text);

            // Apply additional formatting if the text is in the Family_list
            if (isInFamilyList) {
                formattedText = formattedText
                    .replace(/( receptors|neuropeptide )/g, '') // Remove specific substrings
                    .replace(/(-releasing)/g, '-rel.') // Abbreviate specific substrings
                    .replace(/(-concentrating)/g, '-conc.') // Abbreviate specific substrings
                    .replace(/( and )/g, ' & ') // Replace "and" with "&"
                    .replace(/(GPR18, GPR55 & GPR119)/g, 'GPR18, 55 & 119') // Special case formatting
                    .replace(/(Orphan receptors)/g, 'Orphans') // Replace "Class C Orphans"
                    .split("</tspan>")[0] // Keep only part before the first closing tspan tag
                    .split(" (")[0]; // Keep only the part before the first " (" parenthesis
            }

            return formattedText;
        }

       // Bind data and append text elements for the receptor-family headers and receptors of
       // this circle. Class badges are NOT part of this join anymore — see drawClassBadges().
       svg.selectAll(`.GPCRome-text-${level}`)
           .data(contentItems)
           .enter()
           .append("text")
           .attr("class", (d) => {
               let baseClass = `GPCRome-text GPCRome-text-${level}`;  // Add 'GPCRome-text' as a common class
               // Add a family-specific class if this is a receptor-family header
               if (d.type === 'family') {
                   baseClass += " GPCRome-family-label";  // Add this class for family labels
               }
               return baseClass;
           })
           // positionFromAngle's radius depends on isSplit, not on how many lines a label
           // actually renders as: passing isSplit=true for every family label (matching what the
           // split-label lines below use) keeps short, unsplit labels at the exact same radius as
           // split ones — flush with the arc/tick marks — instead of the much-further-inward
           // radius (`kind==='family' && !isSplit`) that used to make short labels look detached
           // from their own receptors. Receptor ticks are unaffected (isSplit doesn't change
           // their radius at all).
           .attr("x", (d) => positionFromAngle(labelAnchorAngle(d), d.type, d.type === 'family').x)
           .attr("y", (d) => positionFromAngle(labelAnchorAngle(d), d.type, d.type === 'family').y)
           .attr("text-anchor", (d) => {
               const angle = (d.midAngle * 180 / Math.PI) - 90;
               return (angle >= -90 && angle < 90) ? "start" : "end";
           })
           // PowerPoint's SVG importer doesn't honor dominant-baseline (browsers do), which made
           // downloaded-SVG labels shift vertically once pasted there. A plain dy offset is a
           // universally-supported attribute both renderers apply identically.
           .attr("dy", "0.35em")
           .attr("transform", (d) => {
               const pos = positionFromAngle(labelAnchorAngle(d), d.type, d.type === 'family');

               // Calculate the angle and determine the text's side (right or left)
               const angle = (d.midAngle * 180 / Math.PI) - 90;

               // For non-headers, flip the text on the left-hand side by 180 degrees
               const rotation = angle >= -90 && angle < 90 ? 0 : 180;

               // Apply the rotation and positioning
               return `rotate(${pos.rotation + rotation}, ${pos.x}, ${pos.y})`;
           })
           .html(d => {
               const labelText = getLabelText(d);
               return GPCRome_formatTextWithHTML(labelText, CircleSubHeaders);
           })
           .style("font-size", FontsizeGlobal)
           .style("font-family", FontStyle)
           .style("fill", "black")

        // Only split (two-line) family labels need a second pass: pass 1 above already rendered
        // every item — family and receptor alike — centered at its own midAngle, which is exactly
        // the desired final position for a short (single-line) family label now that its slot's
        // extra weight (2 units, vs. a receptor's 1) provides breathing room on both sides just by
        // being centered in a wider-than-needed span. So short labels are left untouched here.
        svg.selectAll(".GPCRome-family-label")
            .each(function(d) {
                // svg.selectAll here is unscoped across the whole SVG, so on every subsequent
                // circle's draw this also re-matches already-finalized labels re-appended (via a
                // plain svg.append, with no bound datum) by a *previous* circle's own pass here —
                // skip those rather than reprocessing/removing already-placed labels.
                if (!d || !d.needsSplit) {
                    return;
                }

                const textElement = d3v4.select(this);

                // Determine if the text anchor should be "start" or "end"
                const angle = (d.midAngle * 180 / Math.PI) - 90;  // Calculate the angle based on the item's position
                const additionalRotation = angle >= -90 && angle < 90 ? 0 : 180;  // Conditional rotation adjustment

                // Format the text before checking the length
                const formattedText = GPCRome_formatTextWithHTML(d.label, CircleSubHeaders);

                // Remove the existing (single-line) text element before appending the split elements
                textElement.remove();

                // fontsize
                let family_fontsize = FontsizeGlobal;

                let splitIndex;
                if (formattedText.includes("-")) {
                    // If the text contains a "-", split after the "-"
                    splitIndex = formattedText.indexOf("-",3) + 1;
                } else {
                    // Otherwise, split at the nearest space
                    splitIndex = formattedText.lastIndexOf(" ", formattedText.length-1);
                }
                const firstPart = formattedText.substring(0, splitIndex);  // First part
                const secondPart = formattedText.substring(splitIndex);  // Second part

                // Two evenly-spaced points, one anglePerUnit in from each edge of this item's
                // 3-unit slot (computeClassLayout) — splitFirstAngle is always < splitSecondAngle.
                const firstPos = positionFromAngle(d.splitFirstAngle, 'family', true);
                const secondPos = positionFromAngle(d.splitSecondAngle, 'family', true);

                const off_set = level+1

                if (angle >= -90 && angle < 90) {
                    // Right-hand side: firstPos gets the first part, secondPos gets the second part

                    svg.append("text")
                        .attr("x", firstPos.x)
                        .attr("y", firstPos.y+off_set)
                        .attr("text-anchor", "start")
                        .attr("dy", "0.35em")
                        .attr("transform", `rotate(${firstPos.rotation + additionalRotation}, ${firstPos.x}, ${firstPos.y})`)
                        .attr("class", "GPCRome-family-label-split")
                        .text(firstPart)
                        .style("font-family", FontStyle)
                        .style("font-size",family_fontsize);

                    svg.append("text")
                        .attr("x", secondPos.x)
                        .attr("y", secondPos.y-off_set)
                        .attr("text-anchor", "start")
                        .attr("dy", "0.35em")
                        .attr("transform", `rotate(${secondPos.rotation + additionalRotation}, ${secondPos.x}, ${secondPos.y})`)
                        .attr("class", "GPCRome-family-label-split")
                        .text(secondPart)
                        .style("font-family", FontStyle)
                        .style("font-size", family_fontsize);

                } else {
                    // Left-hand side: secondPos gets the first part, firstPos gets the second part

                    svg.append("text")
                        .attr("x", secondPos.x)
                        .attr("y", secondPos.y+off_set)
                        .attr("text-anchor", "end")
                        .attr("dy", "0.35em")
                        .attr("transform", `rotate(${secondPos.rotation + additionalRotation}, ${secondPos.x}, ${secondPos.y})`)
                        .attr("class", "GPCRome-family-label-split")
                        .text(firstPart)
                        .style("font-family", FontStyle)
                        .style("font-size",family_fontsize);

                    svg.append("text")
                        .attr("x", firstPos.x)
                        .attr("y", firstPos.y-off_set)
                        .attr("text-anchor", "end")
                        .attr("dy", "0.35em")
                        .attr("transform", `rotate(${firstPos.rotation + additionalRotation}, ${firstPos.x}, ${firstPos.y})`)
                        .attr("class", "GPCRome-family-label-split")
                        .text(secondPart)
                        .style("font-family", FontStyle)
                        .style("font-size",family_fontsize);
                }
            });

        // ################
        // ### Coloring ###
        // ################
        // Define color scale for continuous data
        let colorScale;

        if (DataType === "Numeric") {
            if (ColorSetup === 'One') {
            // White to Max (One color)
            colorScale = d3v4.scaleLinear()
                .domain([MinValue, MaxValue])  // Only two points in the domain
                .range([ColorAvg, ColorMax]);  // White to Max color

            } else if (ColorSetup === 'Two') {
            // Min to Max (Two colors)
            colorScale = d3v4.scaleLinear()
                .domain([MinValue, MaxValue])  // Min to Max in the domain
                .range([ColorMin, ColorMax]);  // Min to Max color in range

            } else if (ColorSetup === 'Three') {
            // Min to White to Max (Three colors)
            colorScale = d3v4.scaleLinear()
                .domain([MinValue, AvgValue, MaxValue])  // Min, Avg, Max in the domain
                .range([ColorMin, ColorAvg, ColorMax]);  // Min to White to Max in range
            }
        }

        // Add large hollow pie chart for the entire level. Wedge angles come directly from the
        // startAngle/endAngle computed for each content item in computeClassLayout() — the same
        // angles used to position that item's label — so wedges and labels can never drift out
        // of sync. No wedge is drawn for the reserved class-boundary gaps (previously an
        // invisible fill:none wedge was drawn there anyway, so no visual change).
        const arcGenerator = d3v4.arc()
            .innerRadius(GPCRome_radius - 7)  // Adjust to control the hollow center size
            .outerRadius(GPCRome_radius);  // Adjust to control the thickness of the pie

        svg.selectAll(`.large-hollow-pie-${level}`)
            .data(contentItems)
            .enter()
            .append("path")
            .attr("class", `large-hollow-pie-${level}`)
            .attr("d", arcGenerator)
            .attr("transform", `translate(${width / 2}, ${height / 2})`)
            .style("fill", (d) => {
                const value = DataFill[d.label];

                if (DataType === "Text") {
                    return value || "none";  // Use Color directly
                }

                const numericValue = parseFloat(value);
                return !isNaN(numericValue) ? colorScale(numericValue) : "none";
            })
            .style("stroke", (d) => {
                const value = DataFill[d.label];
                return value != null ? "black" : "none";
            })
            .style("stroke-width", (d) => {
                const value = DataFill[d.label];
                return value === "" ? 0.5 : 0.5;  // Set stroke-width to 0 if the value is an empty string
            });

        // Class badges ("A", "B1", "T2", ...) are drawn last, as their own annotation pass,
        // centered on each class's reserved gap rather than occupying a slot in contentItems.
        // Emits the same GPCRome-text-{level}-highlight class the header text used to carry, so
        // the pill-drawing post-processors in classification/wheel.js and
        // data_mapper/mapper_gpcrome_page.js keep working unmodified against these nodes.
        function drawClassBadges() {
            function headerRadius(classKey) {
                return classKey === "Unclassified" ? (GPCRome_radius - 10) : (GPCRome_radius + 18);
            }

            function headerX(classKey) {
                const theta = (Math.PI / 2) - classGapMidAngle[classKey];
                return width / 2 + Math.cos(theta) * headerRadius(classKey);
            }

            function headerY(classKey) {
                const theta = (Math.PI / 2) - classGapMidAngle[classKey];
                return height / 2 - Math.sin(theta) * headerRadius(classKey);
            }

            svg.selectAll(null)
                .data(CircleHeaders)
                .enter()
                .append("text")
                .attr("class", (d) => `GPCRome-text GPCRome-text-${level} GPCRome-text-${level}-highlight`)
                .attr("x", (d) => headerX(d))
                .attr("y", (d) => headerY(d))
                .attr("text-anchor", (d) => {
                    if (d !== "Unclassified") {
                        return "middle";
                    }
                    const angle = (classGapMidAngle[d] * 180 / Math.PI) - 90;
                    return (angle >= -90 && angle < 90) ? "start" : "end";
                })
                // See the receptor-label pass above for why dominant-baseline was dropped.
                .attr("dy", "0.35em")
                .attr("transform", (d) => `rotate(0, ${headerX(d)}, ${headerY(d)})`)
                .html((d) => GPCRome_formatTextWithHTML(d, CircleSubHeaders))
                .style("font-size", FontsizeClass)
                .style("font-family", FontStyle)
                .style("font-weight", "950")
                .style("fill", "black");
        }

        drawClassBadges();
    }
    // === Legends ===
    let AddBottomHeight = 0;
    if (ShowLegend === true) {
        if (DataType === "Numeric") {
            // Add gradient bar legend for numeric data
            const legendGroup = svg.append("g").attr("class", "legend-gradient-bar");

            const barWidth = GPCRome_styling.LegendbarLength || 300;
            const barHeight = 15;
            const legendPaddingRight = 20;
            const legendPaddingTop = 40;
            const barX = dimensions.width - barWidth - legendPaddingRight;
            const barY = legendPaddingTop;
            const BarFixedDigit = GPCRome_styling.LegendbarDigit || 2;
            const BarFontSize = GPCRome_styling.LegendbarFontsize || "11px";
            const uniqueGradientId = `gradient-bar-${location}`;

            // Create defs and linearGradient
            const defs = svg.append("defs");
            const gradient = defs.append("linearGradient")
                .attr("id", uniqueGradientId)
                .attr("x1", "0%")
                .attr("x2", "100%")
                .attr("y1", "0%")
                .attr("y2", "0%");

            if (ColorSetup === "Three") {
                gradient.append("stop")
                    .attr("offset", "0%")
                    .attr("stop-color", ColorMin);

                gradient.append("stop")
                    .attr("offset", "50%")
                    .attr("stop-color", ColorAvg);

                gradient.append("stop")
                    .attr("offset", "100%")
                    .attr("stop-color", ColorMax);
            } else if (ColorSetup === "Two") {
                gradient.append("stop")
                    .attr("offset", "0%")
                    .attr("stop-color", ColorMin);

                gradient.append("stop")
                    .attr("offset", "100%")
                    .attr("stop-color", ColorMax);
            } else {
                gradient.append("stop")
                    .attr("offset", "0%")
                    .attr("stop-color", "#FFFFFF");

                gradient.append("stop")
                    .attr("offset", "100%")
                    .attr("stop-color", ColorMax);
            }

            // Draw the gradient bar
            legendGroup.append("rect")
                .attr("x", barX)
                .attr("y", barY)
                .attr("width", barWidth)
                .attr("height", barHeight)
                .style("fill", `url(#${uniqueGradientId})`)
                .style("stroke", "black");

            // Add min, avg (if needed), and max labels
            legendGroup.append("text")
                .attr("x", barX)
                .attr("y", barY + barHeight + 15)
                .attr("text-anchor", "start")
                .style("font-size", BarFontSize)
                .text(() => {
                    return Number.isInteger(MinValue)
                        ? parseInt(MinValue)
                        : parseFloat(MinValue).toFixed(BarFixedDigit);
                });

            legendGroup.append("text")
                .attr("x", barX + barWidth)
                .attr("y", barY + barHeight + 15)
                .attr("text-anchor", "end")
                .style("font-size", BarFontSize)
                .text(() => {
                    return Number.isInteger(MaxValue)
                        ? parseInt(MaxValue)
                        : parseFloat(MaxValue).toFixed(BarFixedDigit);
                });

            // Add label above gradient bar if set
            if (LegendLabel) {
                legendGroup.append("text")
                    .attr("x", barX + barWidth / 2)
                    .attr("y", barY - 10)  // 10px above the bar
                    .attr("text-anchor", "middle")
                    .style("font-size", FontsizeGlobal)
                    .text(LegendLabel);
            }

        } else {
            if (LegendMode === "row") {
                const legendGroup = svg.append("g").attr("class", "legend-text-categories");

                const spacingY = 25;
                const padding = 10;
                const startX = 50;
                const startY = dimensions.height + 40;
                let x = startX;
                let y = startY;

                const LegendFontStyle = GPCRome_styling.FontStyle || "Arial";
                const LegendFontSize = GPCRome_styling.LegendbarFontsize || "11px";

                const legendItems = new Set();

                function extractColorData(obj) {
                    if (typeof obj !== "object" || obj === null) return;
                    if ("Color" in obj && "Data" in obj) {
                        const pairKey = `${obj.Color}|||${obj.Data}`;
                        legendItems.add(pairKey);
                    }
                    for (const key in obj) {
                        extractColorData(obj[key]);
                    }
                }

                Object.values(Data).forEach(circleData => {
                    extractColorData(circleData);
                });

                const collator = new Intl.Collator(undefined, { numeric: true, sensitivity: 'base' });
                const sortedItems = Array.from(legendItems)
                    .filter(pair => {
                        const [color, label] = pair.split("|||");
                        return !(color === "#FFFFFF" && label === "Empty");
                    })
                    .sort((a, b) => {
                        const labelA = a.split("|||")[1];
                        const labelB = b.split("|||")[1];
                        return collator.compare(labelA, labelB);
                    });

                const tempText = svg.append("text")
                    .attr("x", -9999)
                    .attr("y", -9999)
                    .style("font-size", LegendFontSize)
                    .style("font-family", LegendFontStyle);

                const maxItemWidth = dimensions.width - 100;

                sortedItems.forEach(pairKey => {
                    let [color, label] = pairKey.split("|||");

                    tempText.text(label);
                    const labelWidth = tempText.node().getComputedTextLength();
                    const totalWidth = 6 * 2 + padding + labelWidth + 20;

                    // Check if item itself is too wide even on a new line
                    if (totalWidth > maxItemWidth) {
                        label = "⚠ Too long label";
                        tempText.text(label);
                    }

                    const fixedWidth = 6 * 2 + padding + tempText.node().getComputedTextLength() + 20;

                    if (x + fixedWidth > dimensions.width - 50) {
                        x = startX;
                        y += spacingY;
                    }

                    legendGroup.append("circle")
                        .attr("cx", x)
                        .attr("cy", y)
                        .attr("r", 6)
                        .style("fill", color)
                        .style("stroke", "black");

                    legendGroup.append("text")
                        .attr("x", x + 10)
                        .attr("y", y + 4)
                        .attr("text-anchor", "start")
                        .style("font-size", LegendFontSize)
                        .style("font-family", LegendFontStyle)
                        .text(label);

                    x += fixedWidth;
                });
                // Clean up measuring element
                tempText.remove();
                //
                const legendBBox =
                    sortedItems.length > 0 ? svg.select(".legend-text-categories").node()?.getBBox() : null;
                if (legendBBox && sortedItems.length > 0) {
                    const centerOffsetX = (dimensions.width - legendBBox.width) / 2 - legendBBox.x;
                    svg.select(".legend-text-categories")
                        .attr("transform", `translate(${centerOffsetX}, 0)`);

                    AddBottomHeight = legendBBox.y + legendBBox.height - dimensions.height + 40;
                }
           } else if (LegendMode === "columns") {

                const legendGroup = svg.append("g").attr("class", "legend-text-categories");

                const legendItems = new Set();
                const numCols = LegendCols;
                const LegendFontStyle = GPCRome_styling.FontStyle || "Arial";
                const LegendFontSize = GPCRome_styling.LegendbarFontsize || "11px";

                function extractColorData(obj) {
                    if (typeof obj !== "object" || obj === null) return;
                    if ("Color" in obj && "Data" in obj) {
                        const pairKey = `${obj.Color}|||${obj.Data}`;
                        legendItems.add(pairKey);
                    }
                    for (const key in obj) {
                        extractColorData(obj[key]);
                    }
                }

                Object.values(Data).forEach(circleData => {
                    extractColorData(circleData);
                });

                const collator = new Intl.Collator(undefined, { numeric: true, sensitivity: 'base' });
                const sortedItems = Array.from(legendItems)
                    .filter(pair => {
                        const [color, label] = pair.split("|||");
                        return !(color === "#FFFFFF" && label === "Empty");
                    })
                    .sort((a, b) => {
                        const labelA = a.split("|||")[1];
                        const labelB = b.split("|||")[1];
                        return collator.compare(labelA, labelB);
                    });

                // Distribute into columns
                let columns = [];

                if (Legendsorting === "Horizontally") {
                    // Row-wise distribution (across columns first)
                    columns = Array.from({ length: numCols }, () => []);
                    sortedItems.forEach((item, index) => {
                        const colIndex = index % numCols;
                        columns[colIndex].push(item);
                    });
                } else {
                    // Default: Column-wise distribution (down each column)
                    const perColumn = Math.floor(sortedItems.length / numCols);
                    const remainder = sortedItems.length % numCols;

                    let index = 0;
                    for (let i = 0; i < numCols; i++) {
                        let count = perColumn + (i < remainder ? 1 : 0);
                        columns.push(sortedItems.slice(index, index + count));
                        index += count;
                    }
                }

                // Compute column widths
                const tempText = legendGroup.append("text")
                    .attr("x", -9999).attr("y", -9999)
                    .style("font-size", LegendFontSize)
                    .style("font-family", LegendFontStyle);

                const colWidths = columns.map(col => {
                    let maxWidth = 0;
                    col.forEach(pair => {
                        const label = pair.split("|||")[1];
                        tempText.text(label);
                        const width = tempText.node().getComputedTextLength();
                        maxWidth = Math.max(maxWidth, width);
                    });
                    return maxWidth + 40; // Add buffer for circle and spacing
                });

                tempText.remove();

                let startX = 50;
                let startY = dimensions.height + 40;

                // Render each column
                columns.forEach((column, colIndex) => {
                    let x = startX;
                    let y = startY;

                    column.forEach(pair => {
                        const [color, label] = pair.split("|||");

                        legendGroup.append("circle")
                            .attr("cx", x).attr("cy", y).attr("r", 6)
                            .style("fill", color).style("stroke", "black");

                        legendGroup.append("text")
                            .attr("x", x + 10)
                            .attr("y", y + 4)
                            .style("font-size", LegendFontSize)
                            .style("font-family", LegendFontStyle)
                            .text(label);

                        y += 25; // Line height
                    });

                    startX += colWidths[colIndex];
                });

                // Center the legend
                const legendBBox =
                    sortedItems.length > 0 ? svg.select(".legend-text-categories").node()?.getBBox() : null;
                if (legendBBox && sortedItems.length > 0) {
                    const centerOffsetX = (dimensions.width + 50 - legendBBox.width) / 2 - legendBBox.x;
                    svg.select(".legend-text-categories")
                        .attr("transform", `translate(${centerOffsetX}, 0)`);

                    AddBottomHeight = legendBBox.y + legendBBox.height - dimensions.height + 40;
                }
            }
        }
    }

    const stylingSaysText =
        String(GPCRome_styling && GPCRome_styling.DataType != null ? GPCRome_styling.DataType : "")
            .trim()
            .toLowerCase() === "text";
    const textWheelRender = stylingSaysText || DataType === "Text";

    // Add padding and update height to match content
    const padding = 10;
    const vbOuterW = dimensions.width + 2 * padding;
    const artboardH =
        dimensions && typeof dimensions.height === "number" && dimensions.height > 0 ? dimensions.height : 1000;

    let addBottom = Number.isFinite(AddBottomHeight) ? AddBottomHeight : 0;
    addBottom = Math.max(addBottom, 0);

    /*
     * Categorical / Text wheels: degenerate bbox + CSS `height:auto` on `#GPCRome_plot svg`
     * (Mapper template) amplify a short viewBox into a ~40px-tall viewport. Clamp both the SVG
     * height attribute and viewBox outer height against the styling object (rerenders can race const DataType).
     */
    if (textWheelRender) {
        addBottom = Math.max(addBottom, 60);
    }

    let newHeight = artboardH + addBottom;
    if (textWheelRender) {
        newHeight = Math.max(newHeight, artboardH + 60, 1040);
    }

    let vbOuterH = newHeight + 2 * padding;
    if (textWheelRender) {
        vbOuterH = Math.max(vbOuterH, Math.round(vbOuterW * 0.96));
        if (vbOuterH > newHeight + 2 * padding) {
            newHeight = vbOuterH - 2 * padding;
        }
    }

    // width/height attributes must equal the viewBox's own dimensions, not the pre-padding
    // "artboard" size — a viewer that falls back to the width/height attributes for scaling
    // (rather than the viewBox) would otherwise crop content the viewBox's negative origin
    // shifts into the padding border.
    svg
        .attr("data-gpcrome-datatype", textWheelRender ? "text" : "numeric")
        .attr("width", vbOuterW)
        .attr("height", vbOuterH)
        .attr("viewBox", `-${padding} -${padding} ${vbOuterW} ${vbOuterH}`);
}

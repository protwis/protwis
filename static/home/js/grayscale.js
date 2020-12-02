function componentToHex(c) {
    var hex = c.toString(16);
    return hex.length === 1 ? "0" + hex : hex;
}

function rgbToHex(r, g, b) {
    return "#" + componentToHex(r) + componentToHex(g) + componentToHex(b);
}

/* Calculate color based on value and scale.
*  @value: numeric value of cell
*  @scale: values relation to set, preferably normalized
*/
function getColor(value, scale, reverse = false) {
    var color;
    var offset;
    if (reverse) {
        scale = 1-scale;
    }
    color = { r: 255, g: 255, b: 255 };
    if (value > 0) {
        color = { r: Math.round(255-(255-153)*scale), g: Math.round(255-(255-153)*scale), b: Math.round(255-(255-153)*scale) }; //gray
    }
    return color;
}

/* Create gray scale for numberic values in table. Assign either the "color-column" class to a cell (td object) to use min-max values for coloring
*  only in that column, or assign the "color-set[int]" class to use min-max values spanning through multiple columns sharing the class set name. 
*  By default high values are colored dark. You can also add the "color-reverse" class to the cell to reverse this coloring.
*
*  @table: The table object
*  @colorSetIds: An array of color-set[int] class name strings
*/
function gray_scale_table(table, colorSetIds = []) {
    // Collect values for columns and sets
    var cols = [];
    var sets = {};
    var setsmaxmin = {};
    colorSetIds.forEach(function(id) {
        sets[String(id)] = [];
        setsmaxmin[String(id)] = [];
    });

    var colIdColorSet = {};
    for (let [i, row] of [...table.find("tbody")[0].rows].entries()) {
        for (let [j, cell] of [...row.cells].entries()) {
            cols[parseInt(j)] = cols[j] || [];
            if (cell.innerText!=='-' && cell.classList.contains("color-column")) {
                cols[j].push(cell.innerText);
            }
            else {
                for (var k = 0; k < colorSetIds.length; k++) {
                    if (cell.classList.contains(colorSetIds[k])) {
                        if (cell.innerText!=='-') {
                            sets[String(colorSetIds[k])].push(cell.innerText);
                        }
                        if (i===0) {
                            colIdColorSet[parseInt(j)] = colorSetIds[k];
                        }
                        break;
                    }
                }    
            }
        }
    }

    // Calculate max, min and abs max values for columns and sets
    var maxmin = [];
    cols.forEach(function(col, index) {
        var max = Math.max.apply(null, col);
        var min = Math.min.apply(null, col);
        var abs_max = Math.max.apply(null, [max, min].map(Math.abs));
        maxmin.push([max, min, abs_max]);
    });
    for (i = 0; i < colorSetIds.length; i++) {
        setsmaxmin[String(colorSetIds[i])] = [Math.max.apply(null, sets[colorSetIds[i]]), Math.min.apply(null, sets[colorSetIds[i]]), Math.max.apply(null, [Math.max.apply(null, sets[colorSetIds[i]]), Math.min.apply(null, sets[colorSetIds[i]])].map(Math.abs))];
    }
    
    var c_maxmin;
    var value;
    var scale;
    var hex;
    var color;
    for (let [i, row] of [...table.find("tbody")[0].rows].entries()) {
        for (let [j, cell] of [...row.cells].entries()) {
            var calculate_color = false;
            var reverse = false;
            scale = NaN;
            // Find corresponding min-max-abs values
            if (cell.classList.contains("color-column")) {
                c_maxmin = maxmin[j];
                calculate_color = true;
            }
            else if (cell.classList.contains(colIdColorSet[j])) {
                c_maxmin = setsmaxmin[colIdColorSet[j]];
                calculate_color = true;
            }
            // Assign color to cell
            if (calculate_color) {
                value = parseFloat(cell.innerText);
                if (cell.classList.contains('color-reverse')) {
                    reverse = true;
                }
                if (!(isNaN(value) || isNaN(c_maxmin[0]) || isNaN(c_maxmin[1]))) {
                    // Normalize data to get in-set color extremes
                    scale = (Math.abs(value)-Math.abs(c_maxmin[1])) / (c_maxmin[2]-Math.abs(c_maxmin[1]));
                    // Calculate color
                    color = getColor(value, scale, reverse);
                    hex = rgbToHex(color.r, color.g, color.b);
                    cell.setAttribute("bgcolor", hex);
                }
                // Non-numeric cells get colored white
                else {
                    color = color = { r: 255, g: 255, b: 255 };
                    hex = rgbToHex(color.r, color.g, color.b);
                    cell.setAttribute("bgcolor", hex);
                }
            }
        }
    }
}
function gray_scale_table(table) {
    // Create gray scale for numberic values in table. Assign either the 'color-column' class to a cell (td object) to use min-max values for coloring
    // only in that column, or assign the 'color-set[int]' class to use min-max values spanning through multiple columns sharing the class set name.
    
    // Find all color-sets
    var colorSetIds = []
    var i = 1;
    while ($('.color-set'+i).length!==0) {
        if (!colorSetIds.includes('color-set'+i)) {
            colorSetIds.push('color-set'+i);
        }
        i = i+1;
    }

    // Collect values for columns and sets
    var cols = []
    var sets = {}
    var setsmaxmin = {};
    colorSetIds.forEach(function(id) {
        sets[id] = [];
        setsmaxmin[id] = [];
    })
    for (let [i, row] of [...table.find("tbody")[0].rows].entries()) {
        for (let [j, cell] of [...row.cells].entries()) {
            cols[j] = cols[j] || [];
            if (cell.innerText!=='-' && cell.classList.contains('color-column')) {
                cols[j].push(cell.innerText);
            }
            else if (cell.innerText!=='-') {
                for (k = 0; k < colorSetIds.length; k++) {
                    if (cell.classList.contains(colorSetIds[k])) {
                        sets[colorSetIds[k]].push(cell.innerText);
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
        maxmin.push([max, min,abs_max]);
    });
    
    for (i = 0; i < colorSetIds.length; i++) {
        setsmaxmin[colorSetIds[i]] = [Math.max.apply(null, sets[colorSetIds[i]]), Math.min.apply(null, sets[colorSetIds[i]]), Math.max.apply(null, [Math.max.apply(null, sets[colorSetIds[i]]), Math.min.apply(null, sets[colorSetIds[i]])].map(Math.abs))];
    }
     
    var cell_count = 0;
    for (let [i, row] of [...table.find("tbody")[0].rows].entries()) {
        for (let [j, cell] of [...row.cells].entries()) {
            var calculate_color = false;
            // Find corresponding min-max-abs values
            if (cell.classList.contains('color-column')) {
                c_maxmin = maxmin[j];
                calculate_color = true;
            }
            else {
                for (k = 0; k < colorSetIds.length; k++) {
                    if (cell.classList.contains(colorSetIds[k])) {
                        c_maxmin = setsmaxmin[colorSetIds[k]];
                        calculate_color = true;
                        break;
                    }
                }
            }
            // Calculate and assign color to cell
            if (calculate_color) {
                value = parseFloat(cell.innerText);
                if (!(isNaN(value) || isNaN(c_maxmin[0]) || isNaN(c_maxmin[1]))) {
                    scale = Math.abs(value) / c_maxmin[2];
                    var color = { r: 255, g: 255, b: 255 };
                    if (value > 0) {
                        color = { r: Math.round(255-(255-153)*scale), g: Math.round(255-(255-153)*scale), b: Math.round(255-(255-153)*scale) }; //gray
                    }
                    var hex = rgbToHex(color.r, color.g, color.b);
                    cell.setAttribute("bgcolor", hex);
                    cell_count++;
                }
                // Non-numeric cells get colored white
                else {
                    var color = { r: 255, g: 255, b: 255 };
                    var hex = rgbToHex(color.r, color.g, color.b);
                    cell.setAttribute("bgcolor", hex);
                }
            }
        }
    }
}

function componentToHex(c) {
    var hex = c.toString(16);
    return hex.length == 1 ? "0" + hex : hex;
}

function rgbToHex(r, g, b) {
    return "#" + componentToHex(r) + componentToHex(g) + componentToHex(b);
}
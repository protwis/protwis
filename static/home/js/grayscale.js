/*eslint complexity: ["error", 20]*/

// Color list based on https://www.w3schools.com/colors/colors_names.asp
const color_set = {
  "aliceblue": "#f0f8ff",
  "antiquewhite": "#faebd7",
  "aqua": "#00ffff",
  "aquamarine": "#7fffd4",
  "azure": "#f0ffff",
  "beige": "#f5f5dc",
  "bisque": "#ffe4c4",
  "black": "#000000",
  "blanchedalmond": "#ffebcd",
  "blue": "#0000ff",
  "blueviolet": "#8a2be2",
  "brown": "#a52a2a",
  "burlywood": "#deb887",
  "cadetblue": "#5f9ea0",
  "chartreuse": "#7fff00",
  "chocolate": "#d2691e",
  "coral": "#ff7f50",
  "cornflowerblue": "#6495ed",
  "cornsilk": "#fff8dc",
  "crimson": "#dc143c",
  "cyan": "#00ffff",
  "darkblue": "#00008b",
  "darkcyan": "#008b8b",
  "darkgoldenrod": "#b8860b",
  "darkgray": "#a9a9a9",
  "darkgreen": "#006400",
  "darkkhaki": "#bdb76b",
  "darkmagenta": "#8b008b",
  "darkolivegreen": "#556b2f",
  "darkorange": "#ff8c00",
  "darkorchid": "#9932cc",
  "darkred": "#8b0000",
  "darksalmon": "#e9967a",
  "darkseagreen": "#8fbc8f",
  "darkslateblue": "#483d8b",
  "darkslategray": "#2f4f4f",
  "darkturquoise": "#00ced1",
  "darkviolet": "#9400d3",
  "deeppink": "#ff1493",
  "deepskyblue": "#00bfff",
  "dimgray": "#696969",
  "dodgerblue": "#1e90ff",
  "firebrick": "#b22222",
  "floralwhite": "#fffaf0",
  "forestgreen": "#228b22",
  "fuchsia": "#ff00ff",
  "gainsboro": "#dcdcdc",
  "ghostwhite": "#f8f8ff",
  "gold": "#ffd700",
  "goldenrod": "#daa520",
  "gray": "#808080",
  "green": "#008000",
  "greenyellow": "#adff2f",
  "honeydew": "#f0fff0",
  "hotpink": "#ff69b4",
  "indianred": "#cd5c5c",
  "indigo": "#4b0082",
  "ivory": "#fffff0",
  "khaki": "#f0e68c",
  "lavender": "#e6e6fa",
  "lavenderblush": "#fff0f5",
  "lawngreen": "#7cfc00",
  "lemonchiffon": "#fffacd",
  "lightblue": "#add8e6",
  "lightcoral": "#f08080",
  "lightcyan": "#e0ffff",
  "lightgoldenrodyellow": "#fafad2",
  "lightgrey": "#d3d3d3",
  "lightgreen": "#90ee90",
  "lightpink": "#ffb6c1",
  "lightsalmon": "#ffa07a",
  "lightseagreen": "#20b2aa",
  "lightskyblue": "#87cefa",
  "lightslategray": "#778899",
  "lightsteelblue": "#b0c4de",
  "lightyellow": "#ffffe0",
  "lime": "#00ff00",
  "limegreen": "#32cd32",
  "linen": "#faf0e6",
  "magenta": "#ff00ff",
  "maroon": "#800000",
  "mediumaquamarine": "#66cdaa",
  "mediumblue": "#0000cd",
  "mediumorchid": "#ba55d3",
  "mediumpurple": "#9370d8",
  "mediumseagreen": "#3cb371",
  "mediumslateblue": "#7b68ee",
  "mediumspringgreen": "#00fa9a",
  "mediumturquoise": "#48d1cc",
  "mediumvioletred": "#c71585",
  "midnightblue": "#191970",
  "mintcream": "#f5fffa",
  "mistyrose": "#ffe4e1",
  "moccasin": "#ffe4b5",
  "navajowhite": "#ffdead",
  "navy": "#000080",
  "oldlace": "#fdf5e6",
  "olive": "#808000",
  "olivedrab": "#6b8e23",
  "orange": "#ffa500",
  "orangered": "#ff4500",
  "orchid": "#da70d6",
  "palegoldenrod": "#eee8aa",
  "palegreen": "#98fb98",
  "paleturquoise": "#afeeee",
  "palevioletred": "#d87093",
  "papayawhip": "#ffefd5",
  "peachpuff": "#ffdab9",
  "peru": "#cd853f",
  "pink": "#ffc0cb",
  "plum": "#dda0dd",
  "powderblue": "#b0e0e6",
  "purple": "#800080",
  "rebeccapurple": "#663399",
  "red": "#ff0000",
  "rosybrown": "#bc8f8f",
  "royalblue": "#4169e1",
  "saddlebrown": "#8b4513",
  "salmon": "#fa8072",
  "sandybrown": "#f4a460",
  "seagreen": "#2e8b57",
  "seashell": "#fff5ee",
  "sienna": "#a0522d",
  "silver": "#c0c0c0",
  "skyblue": "#87ceeb",
  "slateblue": "#6a5acd",
  "slategray": "#708090",
  "snow": "#fffafa",
  "springgreen": "#00ff7f",
  "steelblue": "#4682b4",
  "tan": "#d2b48c",
  "teal": "#008080",
  "thistle": "#d8bfd8",
  "tomato": "#ff6347",
  "turquoise": "#40e0d0",
  "violet": "#ee82ee",
  "wheat": "#f5deb3",
  "white": "#ffffff",
  "whitesmoke": "#f5f5f5",
  "yellow": "#ffff00",
  "yellowgreen": "#9acd32"
};

// Based on https://css-tricks.com/converting-color-spaces-in-javascript/
function hexToRGB(h) {
  let r = 0,
    g = 0,
    b = 0;

  // 3 digits
  if (h.length === 4) {
    r = "0x" + h[1] + h[1];
    g = "0x" + h[2] + h[2];
    b = "0x" + h[3] + h[3];

    // 6 digits
  } else if (h.length === 7) {
    r = "0x" + h[1] + h[2];
    g = "0x" + h[3] + h[4];
    b = "0x" + h[5] + h[6];
  }

  return {
    red: parseInt(r, 16),
    green: parseInt(g, 16),
    blue: parseInt(b, 16)
  };
}

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
    scale = 1 - scale;
  }
  color = {
    r: 255,
    g: 255,
    b: 255
  };

  color = {
    r: Math.round(255 - (255 - 153) * scale),
    g: Math.round(255 - (255 - 153) * scale),
    b: Math.round(255 - (255 - 153) * scale)
  }; //gray

  return color;
}

/*
 * Calculate linear gradient color based on a value in a specific scale.
 *
 * The function takes two or three reference colors in RGB format as an argument
 * and will calculate the gradient color for the provided fraction value
 * @param {number} fraction - Value between 0 and 1 indicating position in the gradient scale
 * @param {bool} reverse - Indicate if the gradient should be inverted or not
 * @param {dict} color1 - Reference color for the lowest value (fraction is 0)
 * @param {dict} color2 - Reference color for the highest value (fraction is 1) or
 *                          halfway value (fraction is 0.5) when three colors are given.
 * @param {dict} [color3] Optional third gradient color for the highest value
 * @return
 */
function colorGradient(fraction, reverse, color1, color2, color3) {
  var fade = fraction;
  if (reverse){
    fade = 1 - fraction;
  }

  // Do we have 3 colors for the gradient?
  if (color3) {
    fade = fade * 2;

    // Find which interval to use and adjust the fade percentage
    if (fade >= 1) {
      fade -= 1;
      color1 = color2;
      color2 = color3;
    }
  }

  var diffRed = color2.red - color1.red;
  var diffGreen = color2.green - color1.green;
  var diffBlue = color2.blue - color1.blue;

  var gradient = {
    r: Math.round(color1.red + (diffRed * fade)),
    g: Math.round(color1.green + (diffGreen * fade)),
    b: Math.round(color1.blue + (diffBlue * fade)),
  };

  return gradient;
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
  var colorGradientAssigned = {};
  var colorGradients = {};
  for (let [i, row] of [...table.find("tbody")[0].rows].entries()) {
    for (let [j, cell] of [...row.cells].entries()) {
      cols[parseInt(j,10)] = cols[j] || [];
      var colored = false;
      if (cell.innerText !== "-" && cell.classList.contains("color-column")) {
        cols[j].push(cell.innerText);
        colored = true;
      } else {
        for (let k = 0; k < colorSetIds.length; k++) {
          if (cell.classList.contains(colorSetIds[k])) {
            if (cell.innerText !== "-") {
              sets[String(colorSetIds[k])].push(cell.innerText);
            }
            if (i === 0) {
              colIdColorSet[parseInt(j,10)] = colorSetIds[k];
            }
            colored = true;
            break;
          }
        }
      }

      // Check for gradients and prepare them
      if (colored) {
        for (let k = 0; k < cell.classList.length; k++) {
          if (cell.classList[k].startsWith("color-gradient_")) {
            var gradScheme = cell.classList[k];
            var gradientName = gradScheme.split("_");
            gradientName.shift(); // remove first part
            if (! (gradScheme in colorGradients)){
              colorGradients[gradScheme] = [];
              for (var l = 0; l < gradientName.length; l++) {
                if (gradientName[l] in color_set){
                  colorGradients[gradScheme].push(hexToRGB(color_set[gradientName[l]]));
                }
              }
            }

            // Also assign gradient to column for increased performance
            colorGradientAssigned[parseInt(j, 10)] = gradScheme;

            // Found the gradient class so we can stop
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
    maxmin.push([max, min]);
    //var abs_max = Math.max.apply(null, [max, min].map(Math.abs));
    //maxmin.push([max, min, abs_max]);
  });
  for (var i = 0; i < colorSetIds.length; i++) {
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
      var abs = false;
      var gradient = false;
      scale = NaN;
      // Find corresponding min-max-abs values
      if (cell.classList.contains("color-column")) {
        c_maxmin = maxmin[j];
        calculate_color = true;
      } else if (cell.classList.contains(colIdColorSet[j])) {
        c_maxmin = setsmaxmin[colIdColorSet[j]];
        calculate_color = true;
      }
      // Assign color to cell
      if (calculate_color) {
        value = parseFloat(cell.innerText);
        if (cell.classList.contains("color-reverse")) {
          reverse = true;
        }
        if (cell.classList.contains("color-abs")) {
          abs = true;
          c_maxmin = [Math.max(Math.abs(c_maxmin[0]), Math.abs(c_maxmin[1])), 0];
          value = Math.abs(value);
        }
        if (j in colorGradientAssigned && cell.classList.contains(colorGradientAssigned[j])) {
          gradient = colorGradients[colorGradientAssigned[j]];
        }

        if (!(isNaN(value) || isNaN(c_maxmin[0]) || isNaN(c_maxmin[1]))) {
          // Normalize data to get in-set color extremes
          scale = (value - c_maxmin[1]) / (c_maxmin[0] - c_maxmin[1]);
          // Calculate color
          if (Array.isArray(gradient)) {
            if (gradient.length === 2){
              color = colorGradient(scale, reverse, gradient[0], gradient[1]);
            } else if (gradient.length === 3){
              color = colorGradient(scale, reverse, gradient[0], gradient[1], gradient[2]);
            }
          } else {
            color = getColor(value, scale, reverse);
          }
          hex = rgbToHex(color.r, color.g, color.b);
          cell.setAttribute("bgcolor", hex);
          cell.style.backgroundColor = hex;
        }
        // Non-numeric cells get colored white
        else {
          color = color = {
            r: 255,
            g: 255,
            b: 255
          };
          hex = rgbToHex(color.r, color.g, color.b);
          cell.setAttribute("bgcolor", hex);
          cell.style.backgroundColor = hex;
        }
      }
    }
  }
}

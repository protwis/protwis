/*global nv,d3,data_cryst_year_container,data_unique_class_cryst_container,data_unique_class_cryst_year_container,data_cryst_class_year_container,data_unique_cryst_container,data_unique_cryst_year_container*/



 function mergeSVG(div) {
     let SVG = $("#"+div).find("svg")[0];
     let h = parseInt($(SVG).attr("height"), 10);
     let w = parseInt($(SVG).attr("width"), 10);
     let legend = $("#legend").find("svg")[0];
     let h2 = parseInt($(legend).attr("height"), 10);
     let w2 = parseInt($(legend).attr("width"), 10);
     let leg_w = (w-w2)/2;
     SVG.setAttribute("height", (h + h2));
     if (w2 > w) {
         SVG.setAttribute("width", (w2));
         leg_w = 0
         //svg_w = Math.abs(w-w2)/2
     } else {
         leg_w = Math.abs(w-w2)/2
         //svg_w = 0
     }
     for (let i = 0; i < legend.children.length; i++) {
         legend.children[i].setAttribute("transform", "translate ("+leg_w.toString()+" " + h.toString()+")");
         $(SVG).append(legend.children[i]);
     }
 }

 function destroyExistingChart(config) {
    Plotly.purge(config.chart_container_id);
}

function flattenDataForPlotly(data) {
    // This function transforms the data stucture from 
    // data = [
    //     {key: "category1", values: [{x: x1, y: y1}, {x: x2, y: y2}, ...]},
    //     {key: "category2", values: [{x: x1, y: y1}, {x: x2, y: y2}, ...]},
    //     ...
    // ]
    // to 
    // data = {
    //     category1: {x: [x1, x2, ...], y: [y1, y2, ...]},
    //     category2: {x: [x1, x2, ...], y: [y1, y2, ...]},
    //     ...
    // }
    return data.map(function (category) {
        values = category.values.reduce(
            function (accumulator, item) {
            return {
                x: [...accumulator.x, item.x],
                y: [...accumulator.y, item.y],
            };
            },
            { x: [], y: [] },
        )
        return {[category.key] : values};
    }).reduce(function (acc, item) {
        return Object.assign(acc, item);
    }, {});
}

function initialiseChart(config)
{
    plot_object = preparePlotlyConfiguration(config);
   
    layout = preparePlotlyLayout(config);

    plotlyConfig = {
        responsive: true,
        displayModeBar: false,
        displaylogo: false,
    };

    return Plotly.newPlot(config.chart_container_id, plot_object, layout, plotlyConfig);

}

function preparePlotlyConfiguration(config) {
    let i=0;
    data = {...config.chart_data}; //clone data to avoid mutating original config   
    for (let [category, data_entry] of Object.entries(data))
    {
        data_entry.name = category;
        data_entry.type = "bar";
        data_entry.marker = {color: config.color_palette[i] };
        i++;
    }    
    return Object.values(data);
}

function generatePlotTitle(config) {

    let count_type_label = "";
    switch (config.count_type) {
        case "distinct_receptors":
        count_type_label = "Count of distinct receptors";
        break;
        case "distinct_gpcrligandcomplexes":
        count_type_label = "Count of distinct GPCR-ligand complexes";
        break;
        case "all_structures":
        count_type_label = "Total structure count (including redundant receptors)";
        break;
    }
return `${count_type_label} by ${config.category}`;
}

function preparePlotlyLayout(config) {  
  plot_title = generatePlotTitle(config);

  return {
    title: { text: plot_title, font: { size: 16 } },
    autosize: false,
    width: 1200,
    height: 600,
    margin: {l:10, t:30, r:0, b:0},
    xaxis: {
      tickangle: config.x_label_angle * -1,
      showgrid: config.show_x_gridlines || false,
      automargin: true,
      dtick: 1
    },
    yaxis: {
      showgrid: config.show_y_gridlines || false,
      automargin: true,
      side: 'right',
        ticks: 'outside',
        tick0: 0,
        ticklen: 4,
        tickwidth: 1,
        tickcolor: '#000',
        showline: true,
        linecolor: '#000',
        linewidth: 1
    },
    annotations: generateBarTotalAnnotations(config.chart_data),
    barmode: config.stack ? "stack" : "group",
    legend: { x: 1.05 },
    responsive: true
  };
}

//Generates annotations for total values on top of stacked bars
generateBarTotalAnnotations = function (data) {
  //Get unique x values/labels in order of appearance
  const xLabels = Object.values(data).reduce((acc, trace) => {
    trace.x.forEach((v) => {
      if (!acc.includes(v)) acc.push(v);
    });
    return acc;
  }, []);

  //Calculate total y value for each x value across all traces
  const yTotal = xLabels.map((xLab) =>
    Object.values(data).reduce((sum, trace) => {
      const i = trace.x.indexOf(xLab);
      return sum + (i !== -1 ? trace.y[i] : 0);
    }, 0),
  );

  //Calculate max Y value and determine annotation offset
  maxY = Math.max(...yTotal);
  annotationOffset = maxY * 0.05; // 5% of max Y for annotation offset

  stack_value_annotations = xLabels.map((xLab, i) => ({
    x: xLab,
    y: yTotal[i] + annotationOffset,
    text: String(yTotal[i]),
    showarrow: false,
    font: { color: "black", size: 12 },
  }));

  return stack_value_annotations;
};

 $(window).on("load", function () {

      //Add event listener for class selection modal
      var confirmBtn = document.getElementById('confirmClassSelection');
      if (confirmBtn) {
          confirmBtn.addEventListener('click', function() {
              chart_configuration.category = "class";
              var selected = [];
              document.querySelectorAll('#classCheckboxList input[type="checkbox"]:checked').forEach(function(cb) {
                  selected.push(cb.value);
              });
              chart_configuration.selection_list = selected;
              $('#classSelectionModal').modal('hide');
              
              generateStatisticsChart(chart_configuration);
          });
      }

      //Add event listener for chemotype selection modal
      var confirmChemotypeBtn = document.getElementById('confirmChemotypeSelection');
      if (confirmChemotypeBtn) {
          confirmChemotypeBtn.addEventListener('click', function() {
              chart_configuration.category = "chemotype";
              var selected = [];
              document.querySelectorAll('#chemotypeCheckboxList input[type="checkbox"]:checked').forEach(function(cb) {
                  selected.push(cb.value);
              });
              chart_configuration.selection_list = selected;
              $('#chemotypeSelectionModal').modal('hide');

              generateStatisticsChart(chart_configuration);
          });
      }
          //Add event listener for family selection modal
      var confirmFamilyBtn = document.getElementById('confirmFamilySelection');
      if (confirmFamilyBtn) {
          confirmFamilyBtn.addEventListener('click', function() {
              chart_configuration.category = "family";
              var selected = [];
              // Only leaf checkboxes carry a family value; parent (class/super-family) rows are skipped.
              document.querySelectorAll('#familyCheckboxList input[type="checkbox"][value]:checked').forEach(function(cb) {
                  selected.push(cb.value);
              });
              chart_configuration.selection_list = selected;
              $('#familySelectionModal').modal('hide');

              generateStatisticsChart(chart_configuration);
          });
      }
          //Unique crystallized receptors graph
    if (typeof chart_configuration != "undefined" && typeof chart_configuration.chart_data != "undefined"){
      generateStatisticsChart(chart_configuration)
    }
});

function generateStatisticsChart(config) {    
    destroyExistingChart(config);
    config = updateChartData(config);
    return initialiseChart(config);
}

function setDisplayCategory(sender) {
    chart_configuration.category = sender.dataset.category;
    chart_configuration.selection_list = [];
    
    generateStatisticsChart(chart_configuration);
}

function displayClassSelection() {
    var class_options = class__all_structures__data.map(x => x.key).sort();
    var checkboxList = document.getElementById('classCheckboxList');
    checkboxList.innerHTML = '';

    class_options.forEach(function(cls) {
        var isChecked = Array.isArray(chart_configuration.selection_list) &&
                        chart_configuration.selection_list.includes(cls);
        var div = document.createElement('div');
        div.className = 'checkbox';
        div.innerHTML = '<label><input type="checkbox" value="' + cls + '"' +
                        (isChecked ? ' checked' : '') + '> ' + cls + '</label>';
        checkboxList.appendChild(div);
    });

    $('#classSelectionModal').modal('show');
}

function displayChemotypeSelection() {
    var chemotype_options = chemotype__all_structures__data.map(x => x.key).sort();
    var checkboxList = document.getElementById('chemotypeCheckboxList');
    checkboxList.innerHTML = '';

    chemotype_options.forEach(function(chemotype) {
        var isChecked = Array.isArray(chart_configuration.selection_list) &&
                        chart_configuration.selection_list.includes(chemotype);
        var div = document.createElement('div');
        div.className = 'checkbox';
        div.innerHTML = '<label><input type="checkbox" value="' + chemotype + '"' +
                        (isChecked ? ' checked' : '') + '> ' + chemotype + '</label>';
        checkboxList.appendChild(div);
    });

    $('#chemotypeSelectionModal').modal('show');
}

function displayFamilySelection() {
    var checkboxList = document.getElementById('familyCheckboxList');
    checkboxList.innerHTML = '';

    var selection = Array.isArray(chart_configuration.selection_list) ? chart_configuration.selection_list : [];

    // family_hierarchy is structured as { gpcr_class: { chemotype: [family, ...] } }.
    // Leaf checkboxes are families (the level the chart filters on); parent rows
    // (gpcr_class, chemotype) are tri-state and select/deselect everything underneath.
    Object.keys(family_hierarchy).sort().forEach(function(gpcr_class) {
        var classNode = buildFamilyTreeNode(gpcr_class, 0);
        var classChildren = document.createElement('div');

        Object.keys(family_hierarchy[gpcr_class]).sort().forEach(function(chemotype) {
            var chemotypeNode = buildFamilyTreeNode(chemotype, 1);
            var chemotypeChildren = document.createElement('div');

            family_hierarchy[gpcr_class][chemotype].slice().sort().forEach(function(family) {
                var leaf = buildFamilyTreeNode(family, 2, family, selection.includes(family));
                chemotypeChildren.appendChild(leaf);
            });

            chemotypeNode.appendChild(chemotypeChildren);
            classChildren.appendChild(chemotypeNode);
        });

        classNode.appendChild(classChildren);
        checkboxList.appendChild(classNode);
    });

    // Wire up parent <-> child synchronisation and initialise parent states.
    checkboxList.querySelectorAll('input[type="checkbox"]').forEach(function(cb) {
        cb.addEventListener('change', function() {
            var container = cb.closest('.family-tree-node');
            // Apply this checkbox's state to all descendant leaf checkboxes.
            container.querySelectorAll('input[type="checkbox"]').forEach(function(child) {
                if (child !== cb) {
                    child.checked = cb.checked;
                    child.indeterminate = false;
                }
            });
            // Recompute the state of all ancestors.
            updateFamilyParentStates(checkboxList);
        });
    });
    updateFamilyParentStates(checkboxList);

    $('#familySelectionModal').modal('show');
}

// Build one row of the family tree. `value` (leaf only) becomes the checkbox value
// collected into selection_list; parent rows carry no value.
function buildFamilyTreeNode(label, depth, value, isChecked) {
    var node = document.createElement('div');
    node.className = 'family-tree-node';

    var div = document.createElement('div');
    div.className = 'checkbox';
    div.style.marginLeft = (depth * 20) + 'px';

    var valueAttr = (typeof value !== 'undefined') ? ' value="' + value + '"' : '';
    div.innerHTML = '<label><input type="checkbox"' + valueAttr +
                    (isChecked ? ' checked' : '') + '> ' + label + '</label>';

    node.appendChild(div);
    return node;
}

// Set each parent checkbox to checked / unchecked / indeterminate based on its leaf descendants.
function updateFamilyParentStates(root) {
    var parents = Array.prototype.slice.call(root.querySelectorAll('.family-tree-node'));
    // Process deepest nodes first so child states are settled before parents read them.
    parents.reverse().forEach(function(node) {
        var cb = node.querySelector(':scope > .checkbox input[type="checkbox"]');
        if (!cb || cb.hasAttribute('value')) return; // leaf nodes have a value
        var descendants = node.querySelectorAll('input[type="checkbox"][value]');
        if (descendants.length === 0) return;
        var checkedCount = 0;
        descendants.forEach(function(d) { if (d.checked) checkedCount++; });
        if (checkedCount === 0) {
            cb.checked = false;
            cb.indeterminate = false;
        } else if (checkedCount === descendants.length) {
            cb.checked = true;
            cb.indeterminate = false;
        } else {
            cb.checked = false;
            cb.indeterminate = true;
        }
    });
}

function setDisplayDataType(sender) {
    chart_configuration.count_type = sender.dataset.counttype;
    generateStatisticsChart(chart_configuration);
}

function updateChartData(chart_config) {
    switch (chart_config.category) {
        case "class":
            switch(chart_config.count_type) {
                case "distinct_receptors":
                    chart_config.chart_data = flattenDataForPlotly(class__distinct_structure__data);
                    break;
                case "distinct_gpcrligandcomplexes":
                    chart_config.chart_data = flattenDataForPlotly(class__distinct_gpcrligand__data);
                    break;
                case "all_structures":
                    chart_config.chart_data = flattenDataForPlotly(class__all_structures__data);
                    break;
            }
            break;
        case "chemotype":
            switch(chart_config.count_type) {
                case "distinct_receptors":
                    chart_config.chart_data = flattenDataForPlotly(chemotype__distinct_structure__data);
                    break;
                case "distinct_gpcrligandcomplexes":
                    chart_config.chart_data = flattenDataForPlotly(chemotype__distinct_gpcrligand__data);
                    break;
                case "all_structures":
                    chart_config.chart_data = flattenDataForPlotly(chemotype__all_structures__data);
                    break;
            }
            break;
        case "family":
            switch(chart_config.count_type) {
                case "distinct_receptors":
                    chart_config.chart_data = flattenDataForPlotly(family__distinct_structure__data);
                    break;
                case "distinct_gpcrligandcomplexes":
                    chart_config.chart_data = flattenDataForPlotly(family__distinct_gpcrligand__data);
                    break;
                case "all_structures":
                    chart_config.chart_data = flattenDataForPlotly(family__all_structures__data);
                    break;
            }
            break;
    }

    if (Array.isArray(chart_config.selection_list) && chart_config.selection_list.length > 0) {
        chart_config.chart_data = Object.fromEntries(
          Object.entries(chart_config.chart_data)
            .map((entry) =>
              chart_config.selection_list.includes(entry[0]) ? entry : null,
            )
            .filter(Boolean), //removes null entries
        );
    }

    chart_config.chart_data = removeLeadingZeroValueYears(chart_config.chart_data);

    return chart_config;
}

function removeLeadingZeroValueYears(chart_data) {
    //This function removes leading years with zero values for all categories to improve readability of the chart
    let years = Object.values(chart_data)[0].x; //Assumes all categories have the same x values (years)
    let numCategories = Object.keys(chart_data).length;
    let leadingZeroCount = 0;

    for (let i = 0; i < years.length; i++) {
        let allZero = true;
        for (let category in chart_data) {
            if (chart_data[category].y[i] !== 0) {
                allZero = false;
                break;
            }
        }
        if (allZero) {
            leadingZeroCount++;
        } else {
            break;
        }
    }

    if (leadingZeroCount > 0) {
        for (let category in chart_data) {
            chart_data[category].x = chart_data[category].x.slice(leadingZeroCount-1); //keep one leading zero year for context
            chart_data[category].y = chart_data[category].y.slice(leadingZeroCount-1);
        }
    }

    return chart_data;
}

function downloadPlotlyChart(sender) {
    let fileFormat = sender.dataset.format;
    let filename = chart_configuration.category + '_' + chart_configuration.count_type + '.' + fileFormat;

    function triggerDownload(dataUrl) {
        let link = document.createElement('a');
        link.href = dataUrl;
        link.download = filename;
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
    }

    switch(fileFormat) {
        case "png":
        generateStatisticsChart(chart_configuration).then(function(chart_reference) {
            Plotly.toImage(chart_reference, {format: 'png', height: 600, width: 1200}).then(function (dataUrl) {
                triggerDownload(dataUrl);
            });    
        });
            break;
        case "jpg":
            generateStatisticsChart(chart_configuration).then(function(chart_reference) {
                Plotly.toImage(chart_reference, {format: 'jpeg', height: 600, width: 1200}).then(function (dataUrl) {
                    triggerDownload(dataUrl);
                });
            });
            break;
        case "svg":
            generateStatisticsChart(chart_configuration).then(function(chart_reference) {
                Plotly.toImage(chart_reference, {format: 'svg', height: 600, width: 1200}).then(function (dataUrl) {
                    triggerDownload(dataUrl);
                });
            });
            break;
        default:
            console.error("Unsupported file format: " + fileFormat);
    }
}


window.zoomCluster = {};

var phylotree;
var treedata;
var radialTree = true;
var doBranchColoring = true;
var hidePDBs = false;
var labelReorder = false;
var treeAnnotations = [];
var couplingAnnotations = [];
var typeClasses = ["Receptor activation state", "", "", "Receptor family", "Ligand type", "GPCR class", "Structure determination method", "Ligand function", "G-protein coupling", "Primary G proteins", "Secondary G proteins", "TM6 opening", "Degree active"]
var dataClasses, colorClasses;
function renderTree(data) {
    dataClasses = [["active", "inactive", "intermediate", "other"],0,0];
    colorClasses = [["#0F0", "#F00", "#F80", "#888"], 0, 0];
    treedata = data // store globally
    var tree = data["tree"]; // contains tree in Newick format

    // data annotations
    treeAnnotations = data["annotations"];

    // Annotations: 0 state, 1 name, 2 fullname, 3 family, 4 ligand type, 5 class, 6 method, 7 ligand function, 8 g-protein coupling
    for (var i = 3; i < 7; i++){
        dataClasses[i] = new Set()
        for (key in treeAnnotations)
          dataClasses[i].add(treeAnnotations[key][i])
        dataClasses[i] = Array.from(dataClasses[i]).sort();
    }
    dataClasses[7] = ["agonist", "partial-agonist", "pam", "antagonist", "inverse-agonist", "nam"]
    dataClasses[8] = ['Gi/Go family', 'Gq/G11 family', 'Gs family', 'G12/G13 family']
    dataClasses[9] = dataClasses[8]
    dataClasses[10] = dataClasses[8]
    dataClasses[11] = ['0%', '50%', '100%'] // percentage
    dataClasses[12] = ['0%', '50%', '100%'] // percentage

    // Receptor family coloring
    if (dataClasses[3].length > 20 && dataClasses[3].length < 32)
      // Hybrid 32 palette: https://lospec.com/palette-list/hybrid32
      colorClasses[3] = ["#593339", "#903d62", "#ae6253", "#dd9c68", "#edce9f", "#c2c237", "#6aba3b", "#3b8f5b", "#335a5c", "#376129", "#979ea8", "#c0c7b7", "#e4edf5", "#38c2d6", "#296291", "#353456", "#613755", "#955b8d", "#d467a2", "#e5df52", "#ec6b24", "#a83135", "#565299", "#645964", "#2e2a35", "#d96f67", "#9d5a33", "#5095e6", "#526626", "#101b21", "#f21e44"]
    else if (dataClasses[3].length >= 32 && dataClasses[3].length < 56)
      // Juice 56 palette: https://lospec.com/palette-list/juice-56
      colorClasses[3] = ["#000005", "#c8e1eb", "#a5becd", "#7891a5", "#55647d", "#37415a", "#191e3c", "#14465a", "#0f7373", "#0fa569", "#41cd73", "#73ff73", "#dc9b78", "#b26247", "#8c3c32", "#5a1423", "#370a14", "#ffd2a5", "#f5a56e", "#e66e46", "#c3412d", "#8c2323", "#410041", "#7d0041", "#aa143c", "#d72d2d", "#f06923", "#ffaa32", "#ffe65a", "#bed72d", "#64a51e", "#237d14", "#0f5519", "#0f3223", "#82ffe1", "#41d7d7", "#14a0cd", "#1469c3", "#0f379b", "#0f0f69", "#3c1e8c", "#642db4", "#a041d7", "#e65ae6", "#ff8cc8", "#820a64", "#b4236e", "#e65078", "#ff8c8c", "#ffcdb4", "#e69b96", "#be6973", "#96465f", "#6e2850"]
    else if (dataClasses[3].length >= 56)
      // SPF 80 palette: https://lospec.com/palette-list/spf-80
      colorClasses[3] = ["#d2ccf3", "#a392d4", "#615476", "#332f3a", "#3f0d76", "#611894", "#8f4bec", "#d291ff", "#edcaff", "#ffaffc", "#f276ff", "#d63be9", "#951cbc", "#680b76", "#30201a", "#473513", "#67541f", "#a79a5f", "#ffe22c", "#fda414", "#ff8d3e", "#f16e03", "#c3680a", "#e0a186", "#db8060", "#c37053", "#a65133", "#88512b", "#6d4734", "#452d25", "#600119", "#900c47", "#974c7a", "#c02214", "#dd3939", "#ff7693", "#ffb7b7", "#fffcdb", "#ffd887", "#f7a357", "#d7863d", "#cc7037", "#b24e2c", "#823314", "#5a260b", "#3a1603", "#0a2563", "#0d396f", "#1c8393", "#42c39c", "#4cd494", "#aaffd8", "#dafffe", "#d1ffcc", "#b6ff8c", "#5bbe61", "#4dae53", "#118448", "#1b744a", "#10594c", "#084339", "#033017", "#221478", "#2c17a5", "#321cbd", "#343af1", "#2274ff", "#39aeff", "#96daff", "#acdecd", "#90d5bd", "#4e9884", "#26795f", "#a2bcc5", "#69849c", "#435655", "#2c3233", "#101010"]

    for (var i = colorClasses.length; i <= 8; i++) {
        colorClasses[i] = []
        if (i != 7) {
          var class_colors;
          if (dataClasses[i].length > 10) {
            class_colors = d3.scale.category20()
                .domain(dataClasses[i])
          } else {
            class_colors = d3.scale.category10()
                .domain(dataClasses[i])
          }
          for (var j = 0; j < dataClasses[i].length; j++) {
            colorClasses[i].push(class_colors(dataClasses[i][j]))
          }
        }
    }
    colorClasses[7] = ["#0F0", "#0F0", "#0F0", "#F00", "#F00", "#F00"]
    colorClasses[9] = colorClasses[8]
    colorClasses[10] = colorClasses[8]
    colorClasses[11] = ["#F00","#F80","#0B0"]
    colorClasses[12] = colorClasses[11]

    // G-protein coupling
    couplingAnnotations = data["Gprot_coupling"];

    // Adjust default coloring selection based on data
    var defaultInner = "Receptor family" // By default use receptor family
    if (dataClasses[5].length >= 2) // Use GPCR class if 2 or more
      defaultInner = "GPCR class"
    else if (dataClasses[4].length >= 4) // Use Ligand type if 4 or more
      defaultInner = "Ligand type"

    // select item
    var innerMenu = d3.select(".datainner").selectAll("li")[0]
    for (var i = 0; i < innerMenu.length; i++) {
      if (innerMenu[i].innerText==defaultInner){
        toggleDataInner({target: innerMenu[i].firstChild})
        break;
      }
    }

    // Annotate and order coupling data
    for (name in treeAnnotations){
      var slug = treeAnnotations[name][8]
      var opening = treeAnnotations[name][9]
      var gprot_likeness = treeAnnotations[name][10]
      treeAnnotations[name][8] = []
      treeAnnotations[name][9] = []
      treeAnnotations[name][10] = []
      treeAnnotations[name][11] = opening
      treeAnnotations[name][12] = gprot_likeness
      var gproteins = dataClasses[8]
      for (g = 0; g < gproteins.length; g++) {
        if (slug in couplingAnnotations){
          if ("primary" in couplingAnnotations[slug] && couplingAnnotations[slug]["primary"].includes(gproteins[g])) {
            treeAnnotations[name][8].push(gproteins[g])
            treeAnnotations[name][9].push(gproteins[g])
          } else if ("secondary" in couplingAnnotations[slug] && couplingAnnotations[slug]["secondary"].includes(gproteins[g])){
            treeAnnotations[name][8].push(gproteins[g])
            treeAnnotations[name][10].push(gproteins[g])
          }
        }
      }
    }

    /*var r = 1200 / 2;
    var spacing = 240;
    var innerRadius = r - spacing // change inner radius of tree with this argument
    var names = 0; // indexing for all nodes*/

    //var plotsize = [document.getElementById('tree-container').offsetWidth*0.75, document.getElementById('tree-container').offsetWidth*0.9]
    var plotsize = window.innerHeight*0.9;
    // maximum is window height - resize if available width is less
    if (document.getElementById('tree-container').offsetWidth*0.9 < plotsize)
      plotsize = document.getElementById('tree-container').offsetWidth*0.9

    plotsize = [plotsize, plotsize]



    // TODO: dynamic setting of the inner_spacing option
    phylotree = d3.layout.phylotree()
    .svg(d3.select("#clustering-tree"))
    .options({
      'left-right-spacing': 'fit-to-size', // fit to given size left-to-right
      'top-bottom-spacing': 'fit-to-size', // fit to given size top-to-bottom
      //'left-right-spacing': 'fixed-step',
      //'top-bottom-spacing': 'fixed-step',
      //'left-offset': 200,
      'restricted-selectable': 'none', // disable normal selection
      'selectable': false,  // make nodes and branches not selectable
      'collapsible': false, // turn off the menu on internal nodes
      'transitions': false, // d3 animations
      'show-scale': false, // show/hide the scale
      'align-tips': true, // align node tips
      'brush': false, // brush
      'reroot': true, // rerooting
      'hide': false, // hiding a subtree
      'zoom': false, // zooming
      'inner_spacing': (clusterMethod==1 ? 1 : 2) // extra spacing between splits in tree (own custom option)
    })
    .radial(radialTree)
    .node_span ('equal')
    .size(plotsize)
    .separation (function (a,b) {return 0.1;});


    // DRAW
    phylotree(tree)
      .style_nodes(nodeStyler)
      .style_edges(branchStyler)
      .layout(false);

    // REDRAW with the correct leaves
    maximumLeafSize()
    phylotree.layout(true)

    // set consensus values for each node
    phylotree.get_nodes().forEach(function(node){
      pdbs = connectedPDBs(node);

      // loop over dataclasses and assign conserved if present
      node['colorClasses'] = JSON.parse(JSON.stringify(treeAnnotations[pdbs[0]]));
      for (i = 1; i < pdbs.length; i++) {
        var pdb = pdbs[i]
        for (j = 0; j < node['colorClasses'].length; j++) {
          if (Array.isArray(node['colorClasses'][j])) {
              node['colorClasses'][j] = node['colorClasses'][j].filter(x => treeAnnotations[pdbs[i]][j].includes(x));
          } else {
              if (node['colorClasses'][j] != treeAnnotations[pdbs[i]][j]) {
                node['colorClasses'][j] = "-"
              }
          }
        }
      }
    });

    if (doBranchColoring){
      // Refresh twice to ensure correct leaf coloring
      d3.layout.phylotree.trigger_refresh(phylotree);
      d3.layout.phylotree.trigger_refresh(phylotree);
    }

    // add legends
    refreshOuterLegend()
    refreshInnerLegend()

    // add custom menu items
    phylotree.get_nodes()
//              .filter(d3.layout.phylotree.is_leafnode)
      .forEach(function(tree_node) {
        d3.layout.phylotree.add_custom_menu(tree_node, // add to this node
          menuGroup1, // display this text for the menu
          function() { addGroup(0, tree_node); }, // on-click callback include a reference to tree_node via transitive closure
          notInGroup1 // condition on when to display the menu
          );
        d3.layout.phylotree.add_custom_menu(tree_node, // add to this node
          menuGroup2, // display this text for the menu
          function() { addGroup(1, tree_node); }, // on-click callback include a reference to tree_node via transitive closure
          notInGroup2 // condition on when to display the menu
          );
        d3.layout.phylotree.add_custom_menu(tree_node, // add to this node
          menuRemoveGroup1, // display this text for the menu
          function() { removeGroup(0, tree_node); }, // on-click callback include a reference to tree_node via transitive closure
          inGroup1 // condition on when to display the menu
          );
        d3.layout.phylotree.add_custom_menu(tree_node, // add to this node
          menuRemoveGroup2, // display this text for the menu
          function() { removeGroup(1, tree_node); }, // on-click callback include a reference to tree_node via transitive closure
          inGroup2 // condition on when to display the menu
          );
    });
}

var maxLeafNodeLenght = -1
function maximumLeafSize(refresh = true) {
  // set to 0
  maxLeafNodeLenght = 0;

  // Find longest label
  d3.select("#clustering-tree").selectAll("text")[0].forEach(
    function(node_label){
      labelSize = node_label.getBBox().width*1.05
      if (labelSize > maxLeafNodeLenght){
        maxLeafNodeLenght = labelSize
      }
    });

  // redraw labels with new max
  if (refresh){
    d3.layout.phylotree.trigger_refresh(phylotree);
    resizeTree()
  }
}

function menuGroup1(node) {
    return "Add to group 1";
}

function menuGroup2(node) {
  return "Add to group 2";
}

function menuRemoveGroup1(node) {
    return "Remove from group 1";
}

function menuRemoveGroup2(node) {
  return "Remove from group 2";
}


function colorGradient(fadeFraction, rgbColor1, rgbColor2, rgbColor3) {
    var color1 = rgbColor1;
    var color2 = rgbColor2;
    var fade = fadeFraction;

    // Do we have 3 colors for the gradient? Need to adjust the params.
    if (rgbColor3) {
      fade = fade * 2;

      // Find which interval to use and adjust the fade percentage
      if (fade >= 1) {
        fade -= 1;
        color1 = rgbColor2;
        color2 = rgbColor3;
      }
    }

    var diffRed = color2.red - color1.red;
    var diffGreen = color2.green - color1.green;
    var diffBlue = color2.blue - color1.blue;

    var gradient = {
      red: parseInt(Math.floor(color1.red + (diffRed * fade)), 10),
      green: parseInt(Math.floor(color1.green + (diffGreen * fade)), 10),
      blue: parseInt(Math.floor(color1.blue + (diffBlue * fade)), 10),
    };

    return rgb2hexCG(gradient.red.toString(16),gradient.green.toString(16),gradient.blue.toString(16));
}

function rgb2hexCG(r,g,b) {
    if (r.length == 1)
        r = '0' + r;
    if (g.length == 1)
        g = '0' + g;
    if (b.length == 1)
        b = '0' + b;

    return '#' + r + g + b;
}

var group1 = []
var group2 = []
function addGroup(selector, node) {
  // Go through nodes and update selection
  function sel(d) {
      if (d3.layout.phylotree.is_leafnode(d)) {
        if (selector > 0)
          group2.push(d.name)
        else
          group1.push(d.name)
      }

      d['group' + selector ] = true
      if (d.children && d.children.length > 0){
        d.children.forEach(sel);
      }
  }
  sel(node);

  // unique PDB codes
  group1 = Array.from(new Set(group1))
  group2 = Array.from(new Set(group2))
  group1.sort()
  group2.sort()

  // Toggle buttons
  if (group1.length > 0 || group2.length > 0){
    $("#CN-button").removeClass("disabled");
//            $("#DN-button").removeClass("disabled");
  }

  // update selection with PDB codes
  var textarea = document.getElementById('input-targets-'+selector);
  if (selector > 0)
    textarea.value = group2.join("\n");
  else
    textarea.value = group1.join("\n");

  // refresh layout with updated selections
  d3.layout.phylotree.trigger_refresh(phylotree);
}

function toggleTreeType(event){
  switch(event.target.innerText){
      case "Circular dendrogram":
        radialTree = true
        break;
      case "Horizontal dendrogram":
        radialTree = false
        break;
    }

    // draw new tree
    phylotree.radial(radialTree).placenodes().update();

    // resizeTree
    resizeTree()

    // update active label on menu items
    event.target.parentNode.parentElement.querySelectorAll( ".active" ).forEach( e =>
      e.classList.remove( "active" ) );
    event.target.parentNode.classList.add("active")
}

function windowResize(){
  var plotsize = window.innerHeight*0.9;
  // maximum is window height - resize if available width is less
  if (document.getElementById('tree-container').offsetWidth*0.9 < plotsize)
    plotsize = document.getElementById('tree-container').offsetWidth*0.9

  var plot = document.getElementById("clustering-tree")
  plot.style.height = plotsize + "px"
  plot.style.width = plotsize + "px"

  // resize plot area and fit zoom level
  if (window.zoomCluster["#clustering-tree"] != null) {
    window.zoomCluster["#clustering-tree"].resize()
    resizeTree()
  }
}

async function resizeTree( eval_string = ""){
  // Partial workaround for horizontal view
  if (!radialTree){
    window.zoomCluster["#clustering-tree"].updateBBox()
    window.zoomCluster["#clustering-tree"].fit()
    window.zoomCluster["#clustering-tree"].zoom(window.zoomCluster["#clustering-tree"].getZoom()*0.9)
    window.zoomCluster["#clustering-tree"].center()
  } else {
    await sleep(100); // slight wait for fully updated SVG
    window.zoomCluster["#clustering-tree"].updateBBox()
    window.zoomCluster["#clustering-tree"].fit()
    window.zoomCluster["#clustering-tree"].zoom(window.zoomCluster["#clustering-tree"].getZoom()*0.99)
    window.zoomCluster["#clustering-tree"].center()
  }
  if (eval_string.length > 0){
    await sleep(100); // slight wait for redraw (e.g. in case of download)
    eval(eval_string)
  }
}

function zoomInTree(){
  window.zoomCluster["#clustering-tree"].zoom(window.zoomCluster["#clustering-tree"].getZoom()/0.9)
}

function zoomOutTree(){
  window.zoomCluster["#clustering-tree"].zoom(window.zoomCluster["#clustering-tree"].getZoom()*0.9)
}

function sleep(ms) {
  return new Promise(resolve => setTimeout(resolve, ms));
}

var displayName = 0
function toggleNames(event){
  switch(event.target.innerText){
      case "IUPHAR":
        displayName = 1
        hidePDBs = true
        labelReorder = true
        break;
      case "IUPHAR (PDB)":
        displayName = 1
        hidePDBs = false
        labelReorder = true
        break;
      case "IUPHAR [inner] (PDB) [outer]":
        displayName = 1
        hidePDBs = false
        labelReorder = false
        break;
      case "UniProt":
        displayName = 0
        hidePDBs = true
        labelReorder = true
        break;
      case "UniProt (PDB)":
        displayName = 0
        hidePDBs = false
        labelReorder = true
        break;
      default:
      case "UniProt [inner] (PDB) [outer]":
        displayName = 0
        hidePDBs = false
        labelReorder = false
        break;
      case "PDB":
        displayName = 2
        hidePDBs = false
        labelReorder = true
      break;
  }

  // update active label on menu items
  event.target.parentNode.parentElement.querySelectorAll( ".active" ).forEach( e =>
    e.classList.remove( "active" ) );
  event.target.parentNode.classList.add("active")

  // refresh layout with updated selections
  d3.layout.phylotree.trigger_refresh(phylotree);
  maximumLeafSize();
}

var displayData = 0
function toggleDataOuter(event){
  var dataName = event.target.innerText

  displayData = menuItem(dataName)

  // Refresh legend
  refreshOuterLegend()

  // update active label on menu items
  event.target.parentNode.parentElement.querySelectorAll( ".active" ).forEach( e =>
    e.classList.remove( "active" ) );
  event.target.parentNode.classList.add("active")

  // Refresh layout with updated selections
  d3.layout.phylotree.trigger_refresh(phylotree);

  // Refresh twice to ensure correct leaf coloring
  if (doBranchColoring)
    d3.layout.phylotree.trigger_refresh(phylotree);

  resizeTree()
}

function refreshOuterLegend(){
  refreshLegend(".outer_legend", displayData)
}

function refreshInnerLegend(){
  refreshLegend(".inner_legend", displayDataInner)
}

function refreshLegend(div_class, selectData){
  d3.select(div_class).html("");
  if (selectData >= 0) {
    var svg_legend_receptor_class = d3.select(div_class).append("svg")
        .attr("height", (dataClasses[selectData].length+1)*20+"px")
        .attr("id", "class_legend_box");

    // add specific title to legend
    svg_legend_receptor_class.append("text")
        .attr("x", 0)
        .attr("y", 7)
        .attr("dy", ".35em")
        .text(typeClasses[selectData])
        .attr("fill", "black")
        .style("font-size", 14)
        .style("font-weight", "bold");

    //// Vertical Legend ////
    var class_legend = svg_legend_receptor_class.selectAll('.class_legend')
        .data(dataClasses[selectData])
        .enter().append('g')
        .attr("class", "class_legend_group")
        .attr("transform", function (d, i) {
            return "translate(0," + i * 20 + ")"
        });


    // Hide the icon if repetition of the same symbol and color
    var opacityArray = Array(dataClasses[selectData].length).fill(1)
    if (selectData==7){
      var colComparison = colorClasses[selectData][0]
      for (var i = 1; i < colorClasses[selectData].length; i++) {
        if (colComparison==colorClasses[selectData][i]){
          opacityArray[i] = 0
        } else {
          colComparison = colorClasses[selectData][i]
        }
      }
    }
    if (div_class==".inner_legend"){
      class_legend.append("circle")
      .attr("cx", 7)
      .attr("cy", 25)
      .attr("r", 5)
      .style("stroke", "#000")
      .style("stroke-width", function (d) {
          return opacityArray[dataClasses[selectData].indexOf(d)]+"px"
      })
      .style("fill-opacity", function (d) {
          return opacityArray[dataClasses[selectData].indexOf(d)]
      })
      .style("fill", function (d) {
          return colorClasses[selectData][dataClasses[selectData].indexOf(d)]
      });
    } else {
      class_legend.append("rect")
      .attr("x", 2)
      .attr("y", 20)
      .attr("width", 10)
      .attr("height", 10)
      .style("stroke", "#000")
      .style("stroke-width", function (d) {
          return opacityArray[dataClasses[selectData].indexOf(d)]+"px"
      })
      .style("fill-opacity", function (d) {
          return opacityArray[dataClasses[selectData].indexOf(d)]
      })
      .style("fill", function (d) {
          return colorClasses[selectData][dataClasses[selectData].indexOf(d)]
      });
    }

    class_legend.append('text')
        .attr("x", 20)
        .attr("y", 30)
        .text(function (d) {
            return (d.replace(new RegExp(" receptors$"), "").replace("/neuropeptide "," / "))
        })
        .style("text-anchor", "start")
        .style("font-size", 13);
  }
}

var displayDataInner = 3
function toggleDataInner(event){
  var dataName = event.target.innerText
  displayDataInner = menuItem(dataName)

  // Refresh legend
  refreshInnerLegend()

  // update active label on menu items
  event.target.parentNode.parentElement.querySelectorAll( ".active" ).forEach( e =>
    e.classList.remove( "active" ) );
  event.target.parentNode.classList.add("active")

  // Refresh layout with updated selections
  d3.layout.phylotree.trigger_refresh(phylotree);

  // Refresh twice to ensure correct leaf coloring
  if (doBranchColoring)
    d3.layout.phylotree.trigger_refresh(phylotree);

  resizeTree()
}

function menuItem(dataName){
  // Annotations: 0 state, 1 name, 2 fullname, 3 family, 4 ligand type, 5 class, 6 method, 7 ligand function
  switch(dataName){
      case "No data":
        return -1;
      case "Receptor activation state":
        return 0;
      case "GPCR class":
        return 5
      case "Ligand type":
        return 4
      case "Receptor family":
        return 3
      case "Structure determination method":
        return 6
      case "All G proteins":
        return 8
      case "Primary G proteins":
        return 9
      case "Secondary G proteins":
        return 10
      case "Receptor cytosolic opening (%)":
      case "TM6 opening (%)":
        return 11
      case "Degree active (%)":
      case "G-protein bound likeness (%)":
        return 12
      case "Ligand function":
        return 7
      default:
        return -1
  }
}

function connectedPDBs(node){
  var descendants = []
  phylotree.descendants(node).forEach(function(leafnode) {
    descendants.push(leafnode.name);
  });
  return descendants;
}

function inGroup1(node) {
  var selected = connectedPDBs(node);
  let intersection = selected.filter(x => group1.includes(x));

  return intersection.length > 0;
}

function notInGroup1(node) {
  var selected = connectedPDBs(node);
  let distinct = selected.filter(x => !group1.includes(x));

  return distinct.length > 0;
}

function inGroup2(node) {
  var selected = connectedPDBs(node);
  let intersection = selected.filter(x => group2.includes(x));

  return intersection.length > 0;
}

function notInGroup2(node) {
  var selected = connectedPDBs(node);
  let distinct = selected.filter(x => !group2.includes(x));

  return distinct.length > 0;
}

function removeGroup(selector, node) {
  // Go through nodes and update selection
  function sel(d) {
      if (d3.layout.phylotree.is_leafnode(d)) {
        if (selector > 0)
          group2 = remove(group2, d.name)
        else
          group1 = remove(group1, d.name)
      }

      d['group' + selector ] = false
      if (d.children && d.children.length > 0){
        d.children.forEach(sel);
      }
  }
  sel(node);

  // unique PDB codes
  group1 = Array.from(new Set(group1))
  group2 = Array.from(new Set(group2))
  group1.sort()
  group2.sort()

  // Toggle buttons
  if (group1.length == 0 && group2.length == 0){
    $("#CN-button").addClass("disabled");
//            $("#DN-button").addClass("disabled");
  }

  // update selection with PDB codes
  var textarea = document.getElementById('input-targets-'+selector);
  if (selector > 0)
    textarea.value = group2.join("\n");
  else
    textarea.value = group1.join("\n");

  // refresh layout with updated selections
  d3.layout.phylotree.trigger_refresh(phylotree);
}

function remove(array, element) {
  return array.filter(el => el !== element);
}

var splitNodesize = 0
var referenceFontSize = 0
function nodeStyler(element, node){
    var scale = 0.8;
    if (d3.layout.phylotree.is_leafnode(node)) {
        if (node.name in treeAnnotations) {
            var node_label = element.select("text");
            var font_size  = parseFloat(node_label.style("font-size"));
            referenceFontSize = font_size;

            // Add inner markers
            element.selectAll("circle").remove()
            var firstColor = "#888";
            if (displayDataInner >= 0) {
              if (!Array.isArray(treeAnnotations[node.name][displayDataInner]))
                treeAnnotations[node.name][displayDataInner] = [treeAnnotations[node.name][displayDataInner]]

              for (var i = 0; i < treeAnnotations[node.name][displayDataInner].length; i++) {
                var currentType = typeClasses[displayDataInner]
                var currentData = treeAnnotations[node.name][displayDataInner][i]
                var classIndex = dataClasses[displayDataInner].indexOf(currentData)
                var currentColor = colorClasses[displayDataInner][classIndex]
                if (i==0){
                  firstColor = currentColor
                }

                splitNodesize = (font_size*scale)/2
                var annotation = element.append("circle")
                annotation.attr("r", (font_size*scale)/2)
                  .style("stroke", "#000")
                  .style("stroke-width", font_size/20+"px")
                  .style("fill", currentColor)
                  .on("mouseover", function(d) { // add tooltip
                      class_tooltip.transition()
                        .style("opacity", .9);
                      class_tooltip.html("<b>" + typeClasses[displayDataInner] + ":</b> " + $(d3.event.target).data("datalabel"))
                        .style("left", (d3.event.pageX) + "px")
                        .style("top", (d3.event.pageY - 28) + "px");
                    })
                  .on("mouseout", function() {
                      class_tooltip.transition().duration(100)
                          .style("opacity", 0);
                  });

                  annotation.attr("data-datalabel", currentData)

                // Shift circles according to placement label
                if (phylotree.radial ()) {
                    var shifter_rect = phylotree.shift_tip(node)[0];
                    var spacer_rect  = (font_size*scale)/2;

                    var setX = shifter_rect > 0 ? shifter_rect - spacer_rect - i*(font_size*scale): shifter_rect + spacer_rect + i*(font_size*scale);
                    annotation.attr("transform", "rotate (" + node.text_angle + ")")
                        .attr ("cx", setX)
                }
              }
            }

            // Color branch-tracer
            var tracer = element.select("line")[0][0];
            tracer.setAttribute("stroke-dasharray", "3,4");
            tracer.setAttribute("stroke-width", referenceFontSize/4);
            tracer.className.baseVal = ""

            // check if in group selection or has assigned coloring
            if ( 'group0' in node && node['group0'] && 'group1' in node && node['group1']){
              tracer.setAttribute("stroke", "purple");
              tracer.setAttribute("stroke-width", referenceFontSize/3);
            } else if ('group0' in node && node['group0']) {
              tracer.setAttribute("stroke", "red");
              tracer.setAttribute("stroke-width", referenceFontSize/3);
            } else if ('group1' in node && node['group1']) {
              tracer.setAttribute("stroke", "blue");
              tracer.setAttribute("stroke-width", referenceFontSize/3);
            } else if (doBranchColoring) {
              // check length of current and if array
              if (displayDataInner >= 0 && treeAnnotations[node.name][displayDataInner].length > 1 && 'colorClasses' in node && Array.isArray(node['colorClasses'][displayDataInner]) && node['colorClasses'][displayDataInner].length >= 1) {
                classIndex = dataClasses[displayDataInner].indexOf(node['colorClasses'][displayDataInner][0])
                firstColor = colorClasses[displayDataInner][classIndex]
              }

              tracer.setAttribute("stroke", firstColor);
            } else {
              tracer.setAttribute("stroke", "#888");
            }

            var labelName
            if (displayName == 1) {
              // fix italics and subscript annotations
              labelName = treeAnnotations[node.name][2].replace(new RegExp(" receptor$"), "")
              labelName = labelName.replace("-adrenoceptor", '')
              labelName = labelName.replace(" receptor-", '-')

              labelName = labelName.replace("<sub>", '</tspan><tspan baseline-shift = "sub">')
              labelName = labelName.replace("</sub>", '</tspan><tspan>')
              labelName = labelName.replace("<i>", '</tspan><tspan font-style = "italic">')
              labelName = labelName.replace("</i>", '</tspan><tspan>')
            } else
              labelName = treeAnnotations[node.name][1].split("_")[0].toUpperCase()

            label = node_label[0][0]

            // Space left and right of the label (instead of &nbsp/&#160 which breaks SVG)
            /*
            label.setAttribute("margin-left", font_size+"px")
            label.setAttribute("margin-right", font_size+"px")*/

            // Extra spacing of labels
            /*var dx_label = parseFloat(label.getAttribute("dx"))
            console.log("DX")
            console.log(dx_label)
            if (dx_label != null && dx_label < 0)
              label.setAttribute("dx", -1 * font_size/2 + "px")
            else if (dx_label != null && dx_label > 0)
              label.setAttribute("dx", font_size/2 + "px")*/

            if (hidePDBs){
              label.innerHTML = "<tspan>" + labelName + "</tspan>"
            } else if (displayName == 2) {
              label.innerHTML = "<tspan>" + node.name + "</tspan>"
            } else {
              label.innerHTML = "<tspan>" + labelName + " (" + node.name + ")" + "</tspan>"
              // Extra: rotate labels when on the left half
              if (!labelReorder && label.getAttribute("dx") != null && label.getAttribute("dx") < 0)
                label.innerHTML = "<tspan>" + "(" + node.name + ") " + labelName + "</tspan>"
            }




            // color labels
            if ( 'group0' in node && node['group0'] && 'group1' in node && node['group1'])
              element.style("fill", "purple");
            else if ('group0' in node && node['group0'])
              element.style("fill", "red");
            else if ('group1' in node && node['group1'])
              element.style("fill", "blue");
            else {
              if (doBranchColoring)
                element.style("fill", firstColor);
              else
                element.style("fill", "black");
            }

            // Add outer markers
            element.selectAll("rect").remove() // remove old data
            if (displayData >= 11) { // Receptor cytosolic opening using growing bars
              if (!Array.isArray(treeAnnotations[node.name][displayData]))
                treeAnnotations[node.name][displayData] = [treeAnnotations[node.name][displayData]]

                var currentType = typeClasses[displayData]
                var percentage = treeAnnotations[node.name][displayData]
                if (percentage > 100)
                  percentage = 100
                else
                  percentage = parseInt(percentage)

                //var currentColor = colorClasses[displayData][0]
                var currentColor = colorGradient(percentage/100,  {red:255, green:0, blue: 0}, {red:255, green:136, blue: 0}, {red:0, green:187, blue: 0})

                // create block
                var annotation = element.append("rect")
                annotation.attr("width", font_size*scale*0.9*5*(percentage+10)/110) // use percentage to resize bar
                  .style("stroke", "#000")
                  .style("stroke-width", font_size/20+"px")
                  .attr("height", font_size*scale*0.9)
                  .attr ("y", -(font_size*scale*0.9)/2)
                  .attr("fill", currentColor)
                  .on("mouseover", function(d) { // add tooltip
                      class_tooltip.transition()
                        .style("opacity", .9);
                      class_tooltip.html("<b>" + currentType + ":</b> " + $(d3.event.target).data("datalabel"))
                        .style("left", (d3.event.pageX) + "px")
                        .style("top", (d3.event.pageY - 28) + "px");
                    })
                  .on("mouseout", function() {
                      class_tooltip.transition().duration(100)
                          .style("opacity", 0);
                  });

                  // label
                  var toollabel = ""
                  annotation.attr("data-datalabel", percentage+"%")

                if (phylotree.radial ()) {
                    var shifter_rect = phylotree.shift_tip(node)[0];
                    var setX = shifter_rect > 0 ? shifter_rect + maxLeafNodeLenght : shifter_rect - maxLeafNodeLenght - font_size*scale*0.9*5*(percentage+10)/110;
                    annotation.attr("transform", "rotate (" + node.text_angle + ")")
                        .attr ("x", setX)
                } else {
                    var x_shift_rect = phylotree.shift_tip(node)[0] + maxLeafNodeLenght + i*(font_size*scale+2);
                    annotation.attr("transform", null).attr("x", function (d, i) { return  x_shift_rect + font_size * i;})
                }
            } else if (displayData >= 0){

              if (!Array.isArray(treeAnnotations[node.name][displayData]))
                treeAnnotations[node.name][displayData] = [treeAnnotations[node.name][displayData]]

              for (var i = 0; i < treeAnnotations[node.name][displayData].length; i++) {
                var currentType = typeClasses[displayData]
                var currentData = treeAnnotations[node.name][displayData][i]
                var classIndex = dataClasses[displayData].indexOf(currentData)
                var currentColor = colorClasses[displayData][classIndex]

                // create block
                var annotation = element.append("rect")
                annotation.attr("width", font_size*scale*0.9)
                  .style("stroke", "#000")
                  .style("stroke-width", font_size/20+"px")
                  .attr("height", font_size*scale*0.9)
                  .attr ("y", -(font_size*scale*0.9)/2)
                  .attr("fill", currentColor)
                  .on("mouseover", function(d) { // add tooltip
                      class_tooltip.transition()
                        .style("opacity", .9);
                      class_tooltip.html("<b>" + currentType + ":</b> " + $(d3.event.target).data("datalabel"))
                        .style("left", (d3.event.pageX) + "px")
                        .style("top", (d3.event.pageY - 28) + "px");
                    })
                  .on("mouseout", function() {
                      class_tooltip.transition().duration(100)
                          .style("opacity", 0);
                  });

                  // label
                  var toollabel = ""
                  annotation.attr("data-datalabel", currentData)

                if (phylotree.radial ()) {
                    var shifter_rect = phylotree.shift_tip(node)[0];
                    var setX = shifter_rect > 0 ? shifter_rect + maxLeafNodeLenght + i*(font_size*scale+2) : shifter_rect - maxLeafNodeLenght - font_size*scale*0.9 - i*(font_size*scale+2);
                    annotation.attr("transform", "rotate (" + node.text_angle + ")")
                        .attr ("x", setX)
                } else {
                    var x_shift_rect = phylotree.shift_tip(node)[0] + maxLeafNodeLenght + i*(font_size*scale+2);
                    annotation.attr("transform", null).attr("x", function (d, i) { return  x_shift_rect + font_size * i;})
                }
              }
            }
        } else {
          console.log("Structure without annotation: " + node.name)
        }
    } else {
        // restyle innner node by coloring according to silhouette score
        silhouette_color = "#FFFFFF";
        if (!isNaN(node.name)) {
          var score = node.name;
          if (score>1) score=1;

          if (score>0){
            silhouette_color = shadeColor2("#AAAAAA", 100-(score*80-20));
          } else {
            silhouette_color = "#FFAAAA";
          }

        }
        element.selectAll("circle").style("fill", silhouette_color)
                                   .attr("r", splitNodesize*0.9)
                                   .attr("stroke-width", (splitNodesize*2)/scale/20*0.9 + "px")

        // Add tooltip
        element.selectAll("circle").on("mouseover", function(d) { // add tooltip
            class_tooltip.transition()
              .style("opacity", .9);
            class_tooltip.html("<b>Silhouette score:</b> " + score)
              .style("left", (d3.event.pageX) + "px")
              .style("top", (d3.event.pageY - 28) + "px");
          })
        .on("mouseout", function() {
            class_tooltip.transition().duration(100)
                .style("opacity", 0);
        });
    }
}

function branchStyler(element, node){
    if ( 'group0' in node.target && node.target['group0'] && 'group1' in node.target && node.target['group1']) {
      element.style("stroke", "purple");
      element.style("stroke-width", referenceFontSize/3)
    } else if ('group0' in node.target && node.target['group0']) {
      element.style("stroke", "red");
      element.style("stroke-width", referenceFontSize/3)
    } else if ('group1' in node.target && node.target['group1']) {
      element.style("stroke", "blue");
      element.style("stroke-width", referenceFontSize/3)
    } else {
      element.style("stroke-width", referenceFontSize/4)
      element.style("stroke", "#888");
      if (doBranchColoring && displayDataInner >= 0 && 'colorClasses' in node.target){
        // color according to shared feature
        if (node.target['colorClasses'][displayDataInner] != "-" && node.target['colorClasses'][displayDataInner].length > 0 ) {
            if (Array.isArray(node.target['colorClasses'][displayDataInner])) {
              // if array with multiple options - first based on previous connection
              if (node.target['colorClasses'][displayDataInner].length > 1){
                  if (node.source['colorClasses'][displayDataInner] != "-" && node.source['colorClasses'][displayDataInner].length > 0){
                    node.target['colorClasses'][displayDataInner] = node.target['colorClasses'][displayDataInner].filter(x => node.source['colorClasses'][displayDataInner].includes(x))
                  }
              }
              var select = dataClasses[displayDataInner].indexOf(node.target['colorClasses'][displayDataInner][0])
              element.style("stroke", colorClasses[displayDataInner][select]);
            } else {
              var select = dataClasses[displayDataInner].indexOf(node.target['colorClasses'][displayDataInner])
              element.style("stroke", colorClasses[displayDataInner][select]);
            }
        }
      }
    }
}

// Based on https://stackoverflow.com/questions/5560248
function shadeColor2(color, percent) {
    var R = parseInt(color.substring(1,3),16);
    var G = parseInt(color.substring(3,5),16);
    var B = parseInt(color.substring(5,7),16);

    R = parseInt(R * percent/100);
    G = parseInt(G * percent/100);
    B = parseInt(B * percent/100);
    R = (R<255)?R:255;
    G = (G<255)?G:255;
    B = (B<255)?B:255;

    var RR = ((R.toString(16).length==1)?"0"+R.toString(16):R.toString(16));
    var GG = ((G.toString(16).length==1)?"0"+G.toString(16):G.toString(16));
    var BB = ((B.toString(16).length==1)?"0"+B.toString(16):B.toString(16));

    return "#"+RR+GG+BB;
}

function downloadNewick(name){
  // Obtain from download data
  if ("tree" in treedata){
    download(name, treedata["tree"])
  }
}

async function download(filename, text) {
  var element = document.createElement('a');
  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
  element.setAttribute('download', filename);

  element.style.display = 'none';
  document.body.appendChild(element);
  // await addition before click - in some cases otherwise there is no download
  await new Promise(r => setTimeout(r, 500));
  element.click();

  document.body.removeChild(element);
}

function downloadSVGWait(svgSelector, name) {
  resizeTree("downloadSVG('"+svgSelector+"', '"+name+"')")
}

function downloadPDFWait(svgSelector, name) {
  resizeTree("downloadPDF('"+svgSelector+"', '"+name+"')")
}

function downloadSVG(svgSelector, name) {
  if ($("#" + svgSelector).length > 0) {
    var ContainerElements = ["svg","g"];
    var RelevantStyles = {"rect":["fill","stroke","stroke-width"],"path":["fill","stroke","stroke-width"],"circle":["fill","stroke","stroke-width"],"line":["stroke","stroke-width"],"text":["fill","font-size","text-anchor"],"polygon":["stroke","fill"]};
    function read_Element(ParentNode, OrigData){
        var Children = ParentNode.childNodes;
        var OrigChildDat = OrigData.childNodes;

        for (var cd = 0; cd < Children.length; cd++){
            var Child = Children[cd];

            var TagName = Child.tagName;
            if (ContainerElements.indexOf(TagName) != -1){
                read_Element(Child, OrigChildDat[cd])
            } else if (TagName in RelevantStyles){
                var StyleDef = window.getComputedStyle(OrigChildDat[cd]);

                var StyleString = "";
                for (var st = 0; st < RelevantStyles[TagName].length; st++){
                    StyleString += RelevantStyles[TagName][st] + ":" + StyleDef.getPropertyValue(RelevantStyles[TagName][st]) + "; ";
                }

                Child.setAttribute("style",StyleString);
            }
        }
    }

    var SVGElem = document.getElementById(svgSelector);
    var oDOM = SVGElem.cloneNode(true);
    read_Element(oDOM, SVGElem)

    var escapedSVG = new XMLSerializer().serializeToString(oDOM);
    var svg = new Blob([escapedSVG], { type: "image/svg+xml;charset=utf-8" });
    var url = URL.createObjectURL(svg);

    downloadURI(url, name);
  }
}

function downloadPDF(svgSelector, name) {
  if ($("#" + svgSelector).length > 0) {
    var ContainerElements = ["svg","g"];
    var RelevantStyles = {"rect":["fill","stroke","stroke-width"],"path":["fill","stroke","stroke-width"],"circle":["fill","stroke","stroke-width"],"line":["stroke","stroke-width"],"text":["fill","font-size","text-anchor"],"polygon":["stroke","fill"]};
    function read_Element(ParentNode, OrigData){
        var Children = ParentNode.childNodes;
        var OrigChildDat = OrigData.childNodes;

        for (var cd = 0; cd < Children.length; cd++){
            var Child = Children[cd];

            var TagName = Child.tagName;
            if (ContainerElements.indexOf(TagName) != -1){
                read_Element(Child, OrigChildDat[cd])
            } else if (TagName in RelevantStyles){
                var StyleDef = window.getComputedStyle(OrigChildDat[cd]);

                var StyleString = "";
                for (var st = 0; st < RelevantStyles[TagName].length; st++){
                    StyleString += RelevantStyles[TagName][st] + ":" + StyleDef.getPropertyValue(RelevantStyles[TagName][st]) + "; ";
                }

                Child.setAttribute("style",StyleString);
            }
        }
    }

    var SVGElem = document.getElementById(svgSelector);
    var oDOM = SVGElem.cloneNode(true);
    read_Element(oDOM, SVGElem)

    var escapedSVG = new XMLSerializer().serializeToString(oDOM);
    $.post('/common/convertsvg', {dataUrl: escapedSVG, filename: name},  function (data) {
            var blob = new Blob([data], { type: "application/pdf" });
            var url = URL.createObjectURL(blob);
            downloadURI(url, name);
        });

//    var svg = new Blob([escapedSVG], { type: "image/svg+xml;charset=utf-8" });
//    var url = URL.createObjectURL(svg);
//    downloadURI(url, name);
  }
}

function downloadPNGWait(svgSelector, name) {
  resizeTree("downloadPNG('"+svgSelector+"', '"+name+"')")
}

function downloadPNG(pngSelector, name) {
  resizeTree()
  var SVGElem = document.getElementById(pngSelector);
  var oDOM = SVGElem.cloneNode(true);

  // Key: strip xmlns tags otherwise this will result in a double definition
  oDOM.removeAttribute("xmlns");

  // Note: keep scale low enough for Chrome to allow download (max size 1.3MB)
  saveSvgAsPng(oDOM, name, {scale: 2});
}

var lastData;
var newCluster = false;
function initializeButtons(selector, treeFunction) {
    $(selector + ' .go-button').click(function() {
        var pdbs = JSON.parse($(selector + ' .crystal-pdb').val());

        // DEBUG input:
        if (pdbs.length == 0)
          pdbs = ["6A94", "6IIV", "6FJ3", "3V2Y", "5T04", "6BQH", "5DSG", "6E3Y", "4Z35", "6B73", "3RZE", "4ZUD", "4XT1", "5ZTY", "5WIV", "3VW7", "5YWY", "4Z9G", "4XNW", "5ZBH", "5TUD", "5VEW", "5UEN", "6AK3", "4BUO", "5X33", "5DHG", "5CXV", "6D32", "5NX2", "4ZJ8", "4MQT", "4DJH", "5TGZ", "6BD4", "5XSZ", "5ZKP", "4OR2", "5ZKQ", "5ZKC", "6HLL", "3OE9", "6MXT", "6C1R", "5WQC", "6D26", "5XRA", "4PY0", "5OM1", "5NJ6", "5XF1", "5UNF", "6N51", "5CGC", "6IGK", "5TZR", "2YCZ", "6DDF", "2YDV", "4XEE", "4IAR", "5X7D", "5EN0", "5XPR", "4U15", "6DO1", "6AKY", "6GPS", "5UZ7", "4RWD", "5G53", "6BQG", "2G87", "6D9H", "6G79", "6CM4", "6B3J", "4DKL", "6DS0", "5VBL", "6H7O", "6M9T", "3PBL"]

        if (pdbs.length > 1) {
            $("#clustering-tree").html("")
            $("#svgloading").remove();
            // empty the SVG field and reload
            $(selector + ' .tree-container').html('<svg id="clustering-tree"></svg><span id="svgloading">Loading...</span>');
            clusterMethod = 0
            if (this.innerHTML == "Go")
              clusterMethod = 0
            else if (this.innerHTML == "Go (distance pairs - normalized)")
                clusterMethod = 0
            else if (this.innerHTML == "Go (distance to 7TM axis)")
              clusterMethod = 1
            else if (this.innerHTML.startsWith("Go (distance to most stable residue)"))
              clusterMethod = 2
            else if (this.innerHTML.startsWith("Go (distance pairs)"))
              clusterMethod = 3
            else if (this.innerHTML.startsWith("Go (distance to membrane mid)"))
              clusterMethod = 4
            else if (this.innerHTML.startsWith("Go (distance to 7TM axis and membrane mid)"))
              clusterMethod = 5
            else if (this.innerHTML.startsWith("Go (distance to \"origin\")"))
              clusterMethod = 6

            $.getJSON( '/contactnetwork/clusteringdata',
            {
                'pdbs': pdbs.join(","),
                'cluster-method': clusterMethod
            },
            /*$.post('/contactnetwork/clusteringdata',
            {
                'pdbs': pdbs,
            },*/
            function( data ) {
                lastData = data

                // create Tee
                treeFunction(data);

                // reset all elements
                $("#svgloading").remove();
                $("#output-group0").removeClass("hidden");
                $("#input-targets-0").val("");
                $("#output-group1").removeClass("hidden");
                $("#input-targets-1").val("");
                $("#submit-group").removeClass("hidden");
                $("#CN-button").addClass("disabled");
                $(".zoombutton-container").removeClass("hidden");
                $(".tree-toggles").removeClass("hidden");

//                        $("#DN-button").addClass("disabled");
                firstTreeClick = true;

                // zoom + pan - internal zoom of phylotree is buggy
                var container = "#clustering-tree";

                // Destroy old zoom on update
                if (window.zoomCluster[container] != null) {
                  if(typeof window.zoomCluster[container].destroy === 'function') {
                    window.zoomCluster[container].destroy();
                  }
                  delete window.zoomCluster[container];
                }

                // Create svg-pan-zoom container
                window.zoomCluster[container] = svgPanZoom(container, {
                    zoomEnabled: false,
                    panEnabled: true,
                    controlIconsEnabled: false,
                    fit: true,
                    center: true,
                    minZoom: 0.1,
                    maxZoom: 10,
                    zoomScaleSensitivity: 0.25,
                    dblClickZoomEnabled: false
                });
            });
        } else {
            toggleAlert()
        }
    });

    $(selector + ' #CN-button').click(function(){ if (!$(this).hasClass("disabled")) { submitToPage("CN"); }});
//            $(selector + ' #DN-button').click(function(){ if (!$(this).hasClass("disabled")) { submitToPage("DN") }});
}

function submitToPage(destination){
  var url = "/contactnetwork/interactions";
  if (destination=="DN")
      url = "/contactnetwork/distances";

  var form = $('<form action="' + url + '" method="post">' +
      '<textarea name="pdbs1" id="submit-pdbs1" />' +
      '<textarea name="pdbs2" id="submit-pdbs2" />' +
      '<input type="hidden" name="csrfmiddlewaretoken" value="' + csrf_token + '" />' +
      '</form>');

  $('body').append(form);

  // set values
  $("#submit-pdbs1").val($("#input-targets-0").val());
  $("#submit-pdbs2").val($("#input-targets-1").val());
  form.submit();
}

function toggleAlert(){
    $(".alert_pdb").fadeTo(2000, 500).slideUp(500, function(){
        $("#success-alert").slideUp(500);
    });
    return false;
}

function initializeTopButtons(selector) {
    // Fullscreen SVG
    $(selector + ' .btn-fullscreen').click(function() {
        fullScreenElement = $(this).parent().next();
        fullScreenElement.css('background-color','white');
        toggleFullScreen(fullScreenElement.get(0));
        windowResize();
    });
}

function downloadDistanceMatrix(filename){
    // create download table
    if (lastData != undefined){
        var header = lastData['dm_labels']
        header.unshift("PDB-code")

        var data = [];
        data.push(header);

        for (var i = 0; i < lastData["distance_matrix"].length ; i++){
          var row = lastData["distance_matrix"][i]
          row.unshift(header[i+1])
          data.push(row);
        }

        // Convert to CSV
        var csv = Papa.unparse(data);

        // Download file
        downloadURI('data:text/csv;charset=UTF-8,' + encodeURI(csv), filename);
    }
}

function downloadURI(uri, name) {
    var link = document.createElement("a");
    link.download = name;
    link.href = uri;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    delete link;
}

function toggleFullScreen(fullScreenElement) {
    if (!document.mozFullScreen && !document.webkitFullScreen) {
        if (fullScreenElement.mozRequestFullScreen) {
            fullScreenElement.mozRequestFullScreen();
        } else {
            fullScreenElement.webkitRequestFullScreen(Element.ALLOW_KEYBOARD_INPUT);
        }
    } else {
        if (document.mozCancelFullScreen) {
          document.mozCancelFullScreen();
        } else {
          document.webkitCancelFullScreen();
        }
    }
}

$('#single-crystal-group-pdbs-modal-table').on('shown.bs.modal', function (e) {
  showPDBtable('#single-crystal-group-pdbs-modal-table');
})

$(document).ready(function() {
    // Get PDBs for table build
    $.get('/contactnetwork/pdbtabledata', function ( data ) {
      $('#single-crystal-group-pdbs-modal-table .tableview').html(data);
      pdbtabledata = data;
    });

    // Single group of PDB files
    initializeButtons('#single-crystal-group-tab', renderTree);
    initializeTopButtons('#single-crystal-group-tab');

    // init tree toggles
    // Button Toggles
    /*$("#radial-layout").on("click", function(e) {
      var radial = $(this).prop("checked");
      phylotree.radial(radial).placenodes().update();
      if (!radial){
        $("#clustering-tree g.svg-pan-zoom_viewport").css("transform", "matrix(.8, 0, 0, .8, 0, 0)")
      } else {
        window.zoomCluster["#clustering-tree"].zoom(1.1)
        window.zoomCluster["#clustering-tree"].center()
      }
    });*/

    $("#colored-edges").on("click", function(e) {
      doBranchColoring = $("#colored-edges").prop("checked")
      // Refresh twice to ensure correct leaf coloring
      d3.layout.phylotree.trigger_refresh(phylotree);
      d3.layout.phylotree.trigger_refresh(phylotree);
    });

/*            $("#hide-pdbs").on("click", function(e) {
      // refresh layout
      hidePDBs = $("#hide-pdbs").prop("checked")
      d3.layout.phylotree.trigger_refresh(phylotree);
      maximumLeafSize();
    });

    $("#label-order").on("click", function(e) {
      // refresh layout
      labelReorder = $("#label-order").prop("checked")
      d3.layout.phylotree.trigger_refresh(phylotree);
    });*/

    // init dropdowns
    $(".dropdown-menu.names").on("click", "li", toggleNames)
    $(".dropdown-menu.treetypes").on("click", "li", toggleTreeType)
    $(".dropdown-menu.dataouter").on("click", "li", toggleDataOuter)
    $(".dropdown-menu.datainner").on("click", "li", toggleDataInner)

    // resize window events
    window.onresize = windowResize
});

var class_tooltip = d3.select("body").append("div")
                      .attr("class", "class_tooltip")
                      .style("opacity", 0);

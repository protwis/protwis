function draw_tree(data, options) {

    var branches = {};
    var branch_offset = 0;
    for (var key in options.branch_length) {
        if (key == options.depth) { continue };
        if (options.label_free.includes(parseInt(key))) {
            branch_offset = branch_offset + 10;
        } else {
            if (options.branch_trunc != 0) {
                branch_offset = branch_offset + 2*options.branch_trunc + 10;
            } else {
                branch_offset = branch_offset + string_pixlen(options.branch_length[key], key);
            }
        }
        branches[key] = branch_offset;
    }
    branches[options.depth] = branch_offset + options.leaf_offset;

    var diameter = 2 * branches[options.depth] + 100;

    var tree = d3.layout.tree()
        .size([360, diameter / 2])
        .separation(function (a, b) { return (a.parent == b.parent ? 1 : 2) / a.depth; });

    var diagonal = d3.svg.diagonal.radial()
        .projection(function (d) { return [d.y, d.x / 180 * Math.PI]; });

    var svg = d3.select('#'+options.anchor).append("svg")
        .attr("width", diameter)
        .attr("height", diameter)
        .attr("id", options.anchor+"_svg")
        .attr("xmlns", "http://www.w3.org/2000/svg");

    var svg_g = svg.append("g")
        .attr("transform", "translate(" + diameter / 2 + "," + diameter / 2 + ")");

    var nodes = tree.nodes(data);

    nodes.forEach(function (d) {
        if (d.depth == 0) {
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
        .style("stroke-width", function (d) { if (d.target.depth > 0) { return 4 - d.target.depth; } else { return 0; } })
        .style("fill-opacity", 0)
        .style("opacity", function (d) {
            if ((d.target.interactions > 0 && d.target.mutations_an > 0) || 1 == 1) { return 0.8 } //|| 1==1
            else if (d.target.interactions > 0) { return 0.5 }
            else if (d.target.mutations_an > 0) { return 0.5 }
            else { return 0.1 };
        });

    var node = svg_g.selectAll(".node")
        .data(nodes)
        .enter().append("g")
        .attr("class", "node")
        .attr("transform", function (d) { if (d.name == '') { return "rotate(" + (d.x) + ")translate(" + d.y + ")"; } else { return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")"; } })
//TODO: add a check to remove circles when nothing is passed (?)
    node.filter(function (d) { return (d.depth == options.depth) })
        .filter(function (d) { return (d.value !== 3000) })
        .append("circle")
        .attr("r", function (d) { if (d.name == '') { return "0" } else { return "4.0" } })
        .style("stroke", "black")
        .style("stroke-width", ".3px")
        .style("fill", function (d) {
            if (d.color && d.depth < options.depth) { return d.color }
            else if ( d.value === 1) {
                return "FireBrick";
            }
            else if ( d.value === 10) {
                return "LightGray";
            }
            else if ( d.value === 20) {
                return "DarkGray";
            }
            else if ( d.value === 30) {
                return "Gray";
            }
            else if ( d.value === 40) {
                return "Black";
            }
            else if (d.value === 100) {
                return 'LightGray';
            }
            else if (d.value === 500) {
                return 'DarkGray';
            }
            else if (d.value === 1000) {
                return 'Gray';
            }
            else if (d.value === 2000) {
                return 'Black';
            }
            else { return "White" };
        })
        .style("opacity", .99);

    node.filter(function (d) { return (d.depth == options.depth) })
        .attr("id", function (d) { if (d.name == '') { return "innerNode" } else { return 'X'+d.name.toUpperCase() } });

    node.append("text")
        .attr("dy", ".31em")
        .attr("name", function (d) { if (d.name == '') { return "branch" } else { return d.name } })
        .attr("text-anchor", function (d) {
            if (d.depth == 3 ) {
              return d.x < 180 ? "start" : "end";
            } else {
              return d.x < 180 ? "end" : "start";
            }
        })
        .attr("transform", function (d) {
            if (d.depth == 3) {
                return d.x < 180 ? "translate(7)" : "rotate(180)translate(-7)";
            } else {
                return d.x < 180 ? "translate(-12)" : "rotate(180)translate(12)";
            }
        })
        .text(function (d) {
            if (d.depth == options.depth) {
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

        .style("font-size", function (d) { if (d.depth < 2) { return "14px" } else if (d.depth == 2) { return "12px" } else { return "10px" } })
        .style("font-family", "Palatino")
        .style("fill", function (d) {
            if (d.color) { return "#111" }
            else { return "#222" };
        }).call(getBB);
    node.filter(function (d) { return (d.depth != options.depth) }).insert("rect", "text")
        .attr("x", function (d) { return d.x < 180 ? d.bbox.x - 12 : d.bbox.x - d.bbox.width - 12; })
        .attr("y", function (d) { return d.bbox.y })
        .attr("width", function (d) { return d.bbox.width })
        .attr("height", function (d) { return d.bbox.height })
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
        if (depth < 2) {
            ctx.font = "20px Palatino"
        } else if (depth == 2) {
            ctx.font = "14px Palatino"
        } else {
            ctx.font = "12px Palatino"
        }
        return parseInt(ctx.measureText(text).width) + 40;
    }

    function getBB(selection) {
        selection.each(function (d) { d.bbox = this.getBBox(); })
    }

    function wrap(text, width) {
        if (width == 0) {
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
            while (word = words.pop()) {
                line.push(word);
                tspan.text(line.join(" "));
                if (tspan.node().getComputedTextLength() > width) {
                    line.pop();
                    tspan.text(line.join(" "));
                    line = [word];
                    tspan = text.append("tspan").attr("x", 0).attr("y", y).attr("dy", ++lineNumber * lineHeight + dy + "em").text(word);
                }
            }
        });
    }
}


/**
* changeLeavesLabels
*
* Function designed to change label names of phylo tree
*
* @location {string} location - svg in which to draw outer circles
* @value {string} value - either IUPHAR or UniProt for label names
* @dict {dictionary} dict - the translation dictionary for label names (from IUPHAR to UniProt and vice versa)
*/

function changeLeavesLabels(location, value, dict){
  // Initialize leaf node length
  maxLeafNodeLenght = 0;
  // Find longest label
  gNodes = d3.select('#'+location).selectAll('g');
  gNodes.each(function(d) {
    if (d3.select(this).attr("id") !== null) {
      name = d3.select(this).attr("id").substring(1);
      labelName = dict[name][0];
      // replaces labels derived from view
      labelName = labelName.replace("-adrenoceptor", '');
      labelName = labelName.replace(" receptor-", '-');
      labelName = labelName.replace("<sub>", '</tspan><tspan baseline-shift = "sub">');
      labelName = labelName.replace("</sub>", '</tspan><tspan>');
      labelName = labelName.replace("<i>", '</tspan><tspan font-style = "italic">');
      labelName = labelName.replace("</i>", '</tspan><tspan>');
      labelName = labelName.replace("Long-wave-sensitive",'LWS');
      labelName = labelName.replace("Medium-wave-sensitive",'MWS');
      labelName = labelName.replace("Short-wave-sensitive",'SWS');
      labelName = labelName.replace("Olfactory", 'OLF');
      labelName = labelName.replace("Calcitonin -like", 'CLR');
      node = d3.select('#'+location).select('#X'+name);
      if (node.size() !== 0){
        if (value === "IUPHAR"){
          node.selectAll("text")[0].forEach(
            function(node_label){
              node_label.innerHTML = labelName;
              labelSize = node_label.getBBox().width*1.05 + 0.5 * 10
              if (labelSize > maxLeafNodeLenght){
                // change initialization label length, needed for outer circles
                maxLeafNodeLenght = labelSize
              }
            });
        } else if (value === "UniProt"){
          node.selectAll("text")[0].forEach(
            function(node_label){
              node_label.innerHTML = name;
              labelSize = node_label.getBBox().width*1.05 + 0.5 * 10
              if (labelSize > maxLeafNodeLenght){
                maxLeafNodeLenght = labelSize
              }
            });
        }
      }
    }
  });
}

/**
* DrawCircles
*
* Function designed to append data circles on the external part of phylo tree
*
* @location {string} location - svg in which to draw outer circles
* @data {Object} data - data provided by the view (json dict usually)
* @starter {integer} starter - the max length of the leaves, to start drawing the circles (calculated by changeLeavesLabels)
* @dict {dictionary} dict - the translation dictionary for color codes
* @fancy {boolean} fancy - set the option for fancy circles
* @clean {boolean} clean - remove the inner circles in the plot
*/

function DrawCircles(location, data, starter, dict, fancy=false, clean=true){

    const pSBC=(p,c0,c1,l)=>{
        let r,g,b,P,f,t,h,i=parseInt,m=Math.round,a=typeof(c1)=="string";
        if(typeof(p)!="number"||p<-1||p>1||typeof(c0)!="string"||(c0[0]!='r'&&c0[0]!='#')||(c1&&!a))return null;
        if(!this.pSBCr)this.pSBCr=(d)=>{
            let n=d.length,x={};
            if(n>9){
                [r,g,b,a]=d=d.split(","),n=d.length;
                if(n<3||n>4)return null;
                x.r=i(r[3]=="a"?r.slice(5):r.slice(4)),x.g=i(g),x.b=i(b),x.a=a?parseFloat(a):-1
            }else{
                if(n==8||n==6||n<4)return null;
                if(n<6)d="#"+d[1]+d[1]+d[2]+d[2]+d[3]+d[3]+(n>4?d[4]+d[4]:"");
                d=i(d.slice(1),16);
                if(n==9||n==5)x.r=d>>24&255,x.g=d>>16&255,x.b=d>>8&255,x.a=m((d&255)/0.255)/1000;
                else x.r=d>>16,x.g=d>>8&255,x.b=d&255,x.a=-1
            }return x};
        h=c0.length>9,h=a?c1.length>9?true:c1=="c"?!h:false:h,f=this.pSBCr(c0),P=p<0,t=c1&&c1!="c"?this.pSBCr(c1):P?{r:0,g:0,b:0,a:-1}:{r:255,g:255,b:255,a:-1},p=P?p*-1:p,P=1-p;
        if(!f||!t)return null;
        if(l)r=m(P*f.r+p*t.r),g=m(P*f.g+p*t.g),b=m(P*f.b+p*t.b);
        else r=m((P*f.r**2+p*t.r**2)**0.5),g=m((P*f.g**2+p*t.g**2)**0.5),b=m((P*f.b**2+p*t.b**2)**0.5);
        a=f.a,t=t.a,f=a>=0||t>=0,a=f?a<0?t:t<0?a:a*P+t*p:0;
        if(h)return"rgb"+(f?"a(":"(")+r+","+g+","+b+(f?","+m(a*1000)/1000:"")+")";
        else return"#"+(4294967296+r*16777216+g*65536+b*256+(f?m(a*255):0)).toString(16).slice(1,f?undefined:-2)
    }

    function componentToHex(c)
    {
        var hex = c.toString(16);
        return hex.length == 1 ? "0" + hex : hex;
    }

    var svg = d3.select('#'+location);
    var node = svg.selectAll(".node");
    if (clean === true) {
      node.selectAll("circle").remove();
    }
    if (fancy === false) {
    var spacer = 8;
      for (var x in data){
        for (var unit in dict){
          if (data[x].indexOf(unit)>= 0) {
            // variable to set the location of the different circle drawing
            multiply = 1+Object.keys(dict).indexOf(unit);
            var leafwithname = svg.selectAll('g[id=X'+x+']')
                .append("circle")
                .attr("r", 3.25)
                .style("fill", dict[unit])
                .attr("transform", "translate(" + (Math.ceil(starter) + multiply*spacer) + ",0)");
              }
        }
      }
    } else {
      var spacer = 10;
      node.selectAll("ellipse").remove();
      for (var x in data){
        keys = Object.keys(data[x]);
        tot = Object.keys(data[x]).reduce((a, b) => a + b, 0);
        var sum = 0;
        for( var el in data[x] ) {
          if( data[x].hasOwnProperty( el ) ) {
            sum += parseFloat( data[x][el] );
          }
        }
        for (var unit in dict){
          if (keys.indexOf(unit)>= 0) {
            percentage = 1-(data[x][unit]/sum);
            color=pSBC (percentage, dict[unit]);
            multiply = 1+Object.keys(dict).indexOf(unit);
            leaf = svg.selectAll('g[id=X'+x+']');
            var leafwithname = svg.selectAll('g[id=X'+x+']')
                .append("ellipse")
                .attr("rx", 4.25)
                .attr("ry", 4.25)
                .style("stroke", dict[unit])
                .style("stroke-width", 1.2)
                .style("fill", color)
                .attr("transform", "translate(" + (Math.ceil(starter) + multiply*spacer) + ",0)");
            }
        }
      }
    }
  }

/**
* FancyCircles
*
* Function designed to append fancy data circles on the external part of phylo tree with
*
* @location {string} location - svg in which to draw outer circles
* @data {Object} data - data provided by the view (json dict usually)
* @starter {integer} starter - the max length of the leaves, to start drawing the circles (calculated by changeLeavesLabels)
* @dict {dictionary} dict - the translation dictionary for color codes
*/

function FancyCircles(location, data, starter, dict){
    var spacer = 8;
    var svg = d3.select('#'+location);
    var node = svg.selectAll(".node");
    for (var x in data){
      keys = Object.keys(data[x]);
      max = Object.keys(data[x]).reduce(function(a, b){ return data[x][a] > data[x][b] ? a : b });
      for (var unit in dict){
        if (keys.indexOf(unit)>= 0) {
          // variable to set the location of the different circle drawing
          multiply = 1+Object.keys(dict).indexOf(unit);
          if (unit === max) {
            var leafwithname = svg.selectAll('g[id=X'+x+']')
                .append("circle")
                .attr("r", 3.25)
                .style("stroke", dict[unit])
                .style("stroke-width", 1)
                .style("fill", dict[unit])
                .attr("transform", "translate(" + (Math.ceil(starter) + multiply*spacer) + ",0)");
          }
          else {
            var leafwithname = svg.selectAll('g[id=X'+x+']')
                .append("circle")
                .attr("r", 3.25)
                .style("stroke", dict[unit])
                .style("stroke-width", 1)
                .style("fill", "white")
                .attr("transform", "translate(" + (Math.ceil(starter) + multiply*spacer) + ",0)");
              }
          }
      }
    }
  }

/**
* draw_cluster
*
* Function designed to draw horizontal cluster
*
*
* @data {Object} data - data provided by the view (json dict usually)
* @options {dictionary} options - options provided for the sake of depth, data info and text anchor
* @trim {boolean} trim - boolean value to select trimming empty leaves
*/

function draw_cluster(data, options, trim=true) {

    if (trim==true) {
      for (var family in data['children']) {
        for (var child in data['children'][family]['children']){
          var i = data['children'][family]['children'][child]['children'].length;
          while (i-- ){
            if (data['children'][family]['children'][child]['children'][i]['value'] == 0) {
              data['children'][family]['children'][child]['children'].splice(i, 1);
            }
          }
        }
        var k = data['children'][family]['children'].length;
        while (k-- ){
          if (data['children'][family]['children'][k]['children'].length == 0) {
            data['children'][family]['children'].splice(k, 1);
          }
        }
      }
      var j = data['children'].length;
      while (j-- ){
        if (data['children'][j]['children'].length == 0) {
          data['children'].splice(j, 1);
        }
      }
    }

    var branches = {};
    var branch_offset = 0;
    for (var key in options.branch_length) {
        if (key == options.depth) { continue };
        if (options.label_free.includes(parseInt(key))) {
            branch_offset = branch_offset + 10;
        } else {
            if (options.branch_trunc != 0) {
                branch_offset = branch_offset + 2*options.branch_trunc + 10;
            } else {
                branch_offset = branch_offset + string_pixlen(options.branch_length[key], key);
            }
        }
        branches[key] = branch_offset;
    }
    branches[options.depth] = branch_offset + options.leaf_offset;

    var dimension = 2 * branches[options.depth] + 100;

    var cluster = d3.layout.cluster()
        .size([dimension/3, (dimension/3) - 100])
        .separation(function (a, b) { return (a.parent == b.parent ? 1 : 2) / a.depth; });

    var diagonal = d3.svg.diagonal()
        .projection(function (d) {
          return [d.y, d.x];
        });

    var svg = d3.select('#'+options.anchor).append("svg")
        .attr("width", dimension)
        .attr("height", dimension/2)
        .attr("id", options.anchor+"_svg")
        .attr("transform", "translate(0," + (dimension/12) + ")");



    var nodes = cluster.nodes(data);

    nodes.forEach(function (d) {
        if (d.depth == 0) {
            d.y = 0
        } else {
            d.y = branches[d.depth]
        }
    });

    var links = cluster.links(nodes);

    var link = svg.selectAll(".link")
        .data(links)
        .enter()
        .append("path")
        .each(function (d) { d.target.linkNode = this; })
        .attr("class", "link_continuous")
        .attr("d", diagonal)
        .style("stroke", function (d) { return d.target.color; })
        .style("stroke-width", function (d) { if (d.target.depth > 0) { return 4 - d.target.depth; } else { return 0; } })
        .style("fill-opacity", 0);

    var node = svg.selectAll(".node")
        .data(nodes)
        .enter()
        .append("g")
        .attr("class", "node")
        .attr("transform", function (d) {
      return "translate(" + d.y + "," + d.x + ")";
        })

    node.filter(function (d) { return (d.depth == options.depth) })
        .append("circle")
        .attr("r", "4.0")
        .style("stroke", "black")
        .style("stroke-width", ".3px")
        .style("fill", function (d) {
            if (d.color && d.depth < options.depth) { return d.color }
            else if ( d.value === 10) {
                return "LightGray";
            }
            else if ( d.value === 20) {
                return "DarkGray";
            }
            else if ( d.value === 30) {
                return "Gray";
            }
            else if ( d.value === 40) {
                return "Black";
            }
            else { return "White" };
        })
        .style("opacity", .99);

    node.filter(function (d) { return (d.depth == options.depth) })
        .attr("id", function (d) { if (d.name == '') { return "innerNode" } else { return 'X'+d.name.toUpperCase() } });

    node.append("text")
        .attr("dy", "3")
        .attr("dx", function(d) { return d.children ? -8 : 8; })
        .attr("name", function (d) { if (d.name == '') { return "branch" } else { return d.name } })
        .attr("text-anchor", function (d) {
              return d.children ? "end" : "start";
        })
        .text(function (d) {
            if (d.depth == options.depth) {
                return d.name.toUpperCase();
            } else if (options.label_free.includes(d.depth)) {
                return "";
            } else if (d.depth > 0) {
                return d.name;
            } else {
                return "";
            }
        })
        // .call(wrap, options.branch_trunc)
        .style("font-size", "12px")
        .style("font-family", "Palatino")
        .style("fill", function (d) {
            if (d.color) { return "#111" }
            else { return "#222" };
        }).call(getBB);
    node.filter(function (d) { return (d.depth != options.depth) }).insert("rect", "text")
        .attr("x", function (d) { return d.bbox.x })
        .attr("y", function (d) { return d.bbox.y })
        .attr("width", function (d) { return d.bbox.width })
        .attr("height", function (d) { return d.bbox.height })
        .style("fill", "#FFF");

    function string_pixlen(text, depth) {
        var canvas = document.createElement('canvas');
        var ctx = canvas.getContext("2d");
        if (depth < 2) {
            ctx.font = "22px Palatino"
        } else if (depth == 2) {
            ctx.font = "16px Palatino"
        } else {
            ctx.font = "14px Palatino"
        }
        return parseInt(ctx.measureText(text).width) - 100;
    }

    function getBB(selection) {
        selection.each(function (d) { d.bbox = this.getBBox(); })
    }

    function wrap(text, width) {
        if (width == 0) {
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
                dy = parseFloat(text.attr("dy"));
                // tspan = text.text(null).append("tspan").attr("x", 0).attr("y", y).attr("dy", dy + "em");
            while (word = words.pop()) {
                line.push(word);
                tspan.text(line.join(" "));
                if (tspan.node().getComputedTextLength() > width) {
                    line.pop();
                    tspan.text(line.join(" "));
                    line = [word];
                    tspan = text.append("tspan").attr("x", 0).attr("y", y).attr("dy", ++lineNumber * lineHeight + dy + "em").text(word);
                }
            }
        });
    }

}



// Heatmap (circles_bal, heatmap_data, master_dict, buttonObjects);
function draw_heatmap(square_data, data, bible, options, location, element_id, legend_label, trim=true) {
  // set the dimensions and margins of the graph

  var support = (function () {
  	if (!window.DOMParser) return false;
  	var parser = new DOMParser();
  	try {
  		parser.parseFromString('x', 'text/html');
  	} catch(err) {
  		return false;
  	}
  	return true;
  })();

  function oddOrEven(x) {
    return ( x & 1 ) ? "odd" : "even";
  }

  var textToHTML = function (str) {

  	// check for DOMParser support
  	if (support) {
  		var parser = new DOMParser();
  		var doc = parser.parseFromString(str, 'text/html');
  		return doc.body.innerHTML;
  	}

  	// Otherwise, create div and append HTML
  	var dom = document.createElement('div');
  	dom.innerHTML = str;
  	return dom;

  };

  var counting = function (var1, var2, var3) {
    for (var i=0; i < var1.length; i++){
      id = var1[i].replace(' (','').replace(')','');
      var2[id] = (i + 1)*var3;
    }
    return var2;
  };

  var margin = {top: 30, right: 30, bottom: 30, left: 30};

  // Labels of row and columns
  var rows = [];
  var html_converter = [];
  var cols = [];
  var colorCols = ['color'];
  var chartData = [];
  var colorData = [];
  var clone = [];
  if (typeof structuredClone === "function") {
    clone = structuredClone(data);
  } else {
    clone = JSON.parse(JSON.stringify(data));
  }
  var gpcrClasses = [];
  var gpcrLigandType = [];
  var gpcrReceptorFamily = [];
  var gpcrName = [];
  var sorterDict = {};
  var filtered_square_data = {};
  var highest_value = 0;
  var rows_to_uniprot = {};

  if (trim==true) {
    for (var prot_class in clone){
      for (var entry in clone[prot_class]) {
        for (var family in clone[prot_class][entry]['children']) {
          var i = clone[prot_class][entry]['children'][family]['children'].length;
          while (i-- ){
            if (clone[prot_class][entry]['children'][family]['children'][i]['value'] == 0) {
              clone[prot_class][entry]['children'][family]['children'].splice(i, 1);
            }
          }
        }
        var k = clone[prot_class][entry]['children'].length;
        while (k-- ){
          if (clone[prot_class][entry]['children'][k]['children'].length == 0) {
            clone[prot_class][entry]['children'].splice(k, 1);
          }
        }
      }
      var l = clone[prot_class].length;
      while (l-- ){
        if (clone[prot_class][l]['children'].length == 0) {
          clone[prot_class].splice(l, 1);
        }
      }
    }
    for (var prot_class in clone) {
      gpcrClasses.push(prot_class);
      if (options['classClick'].includes(prot_class)){
        for (var entry in clone[prot_class]){
          gpcrLigandType.push(clone[prot_class][entry]['name']);
          for (var family in clone[prot_class][entry]['children']) {
            gpcrReceptorFamily.push(clone[prot_class][entry]['children'][family]['name']);
            for (var i in clone[prot_class][entry]['children'][family]['children']){
              tmp = clone[prot_class][entry]['children'][family]['children'][i]['name'].toUpperCase();
              row_val = bible[clone[prot_class][entry]['children'][family]['children'][i]['name'].toUpperCase()][options['namesClick']];
              gpcrName.push(tmp);
              filtered_square_data[tmp] = square_data[tmp];
              rows_to_uniprot[row_val] = tmp;
              rows.push(row_val);
              html_converter.push(
                {
                  row: bible[clone[prot_class][entry]['children'][family]['children'][i]['name'].toUpperCase()][options['namesClick']],
                  html: textToHTML(bible[clone[prot_class][entry]['children'][family]['children'][i]['name'].toUpperCase()][options['namesClick']])
                }
              );
              // console.log(textToHTML(bible[clone[prot_class][entry]['children'][family]['children'][i]['name'].toUpperCase()][options['namesClick']]));
              colorData.push(
                { row: bible[clone[prot_class][entry]['children'][family]['children'][i]['name'].toUpperCase()][options['namesClick']],
                  col: 'color',
                  value: bible[clone[prot_class][entry]['children'][family]['children'][i]['name'].toUpperCase()][options['colorClick']]
                }
              );
            }
          }
        }
      }
    }
  } else {
    for (var prot_class in data) {
      gpcrClasses.push(prot_class);
      if (options['classClick'].includes(prot_class)){
        for (var entry in data[prot_class]){
          gpcrLigandType.push(data[prot_class][entry]['name']);
          for (var family in data[prot_class][entry]['children']) {
            gpcrReceptorFamily.push(data[prot_class][entry]['children'][family]['name']);
            for (var i in data[prot_class][entry]['children'][family]['children']){
              gpcrName.push(data[prot_class][entry]['children'][family]['children'][i]['name'].toUpperCase());
              rows.push(
                textToHTML(bible[data[prot_class][entry]['children'][family]['children'][i]['name'].toUpperCase()][options['namesClick']])
              );
              colorData.push(
                { row: bible[data[prot_class][entry]['children'][family]['children'][i]['name'].toUpperCase()][options['namesClick']],
                  col: 'color',
                  value: bible[data[prot_class][entry]['children'][family]['children'][i]['name'].toUpperCase()][options['colorClick']]
                }
              );
            }
          }
        }
      }
    }
  }

  sorterDict = counting(gpcrClasses.sort().reverse(), sorterDict, 1000000);
  sorterDict = counting(gpcrLigandType.sort().reverse(), sorterDict, 10000);
  sorterDict = counting(gpcrReceptorFamily.sort().reverse(), sorterDict, 100);
  sorterDict = counting(gpcrName.sort().reverse(), sorterDict, 1);

  var sorted_rows = [];
  var tmpdict = {};
  for (var i=0; i < rows.length; i++){
    k = 0;
    for (var id in options['sortClick']){
      name = options['sortClick'][id];
      fam = bible[rows_to_uniprot[rows[i]]][name].replace(' receptors','').replace('(','').replace(')','');
      k = k + sorterDict[fam];
    }
    tmpdict[rows[i]] = k;

  }
  // from dictionary to array
  var tmpItems = Object.keys(tmpdict).map(function(key) {
    return [key, tmpdict[key]];
  });
  // sorting the generated array based on the calculated values
  tmpItems.sort(function(first, second) {
    return second[1] - first[1];
  });
  // pushing the sorted names to the rows
  for (var i=0; i < tmpItems.length; i++){
    sorted_rows.push(tmpItems[i][0]);
  }
  // reversing the array for plotting purposes
  sorted_rows.reverse();
  // pushing data to the chart to be plotted
  for (var receptor in filtered_square_data){
    for (var subtype in filtered_square_data[receptor]){
      if (subtype != 'null'){
        cols.push(subtype);
        chartData.push({ row: bible[receptor][options['namesClick']], col: subtype, value: filtered_square_data[receptor][subtype]});
        if (filtered_square_data[receptor][subtype] > highest_value){
          highest_value = filtered_square_data[receptor][subtype];
        }
      }
    }
  }

  var cols = cols.filter(function(item, pos){
                          return cols.indexOf(item)== pos;
                          });
  cols.sort();

  var width = (35 * cols.length);
  var height = (20 * rows.length);

  // append the svg object to the body of the page
  var svg_home = d3v4.select("#" + location)
              .append("svg")
              .attr("width", width + margin.left + margin.right)
              .attr("height", height + (margin.top*5))
              .attr("transform", "translate(0,-" + margin.bottom + ")")
              .attr("id", element_id)

  var greyscale = [ 'rgb(255,255,255)', 'rgb(0,0,0)' ];

  var legend = svg_home.append('defs')
                       .append('linearGradient')
                       .attr('id', 'grad_' + element_id )
                       .attr('x1', '0%')
                       .attr('x2', '100%')
                       .attr('y1', '0%')
                       .attr('y2', '0%');

   legend.selectAll('stop')
         .data(greyscale)
         .enter()
         .append('stop')
         .style('stop-color', function(d){ return d; })
         .attr('offset', function(d,i){
                return 100 * (i / (greyscale.length - 1)) + '%';
              })


  var legend_svg = svg_home.append("g")
                    .attr("transform", "translate(0," + (height+100)  + ")");

  var color_svg = svg_home.append("g")
                    .attr("transform", "translate(" + (margin.left*1.5) + ",0)");

  var svg = svg_home.append("g")
                    .attr("transform", "translate(" + (margin.left*2.5) + ",0)");

  legend_svg.append("text")
            // .attr('x', 75)
            .attr('y', 10)
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
            .attr('width', (width*0.9) + margin.left + margin.right)
            .attr('height', 20)
            .style('fill', 'url(#grad_' + element_id + ')');

  legend_svg.append("text")
            .attr('x', width + margin.left + (margin.right*1.75))
            .attr('y', 30)
            .style("font", "14px sans-serif")
            .style("font-weight", "bold")
            .text(highest_value);

  legendLabel = legend_svg.select('text')['_groups'][0][0].getBBox().width*1.05 + 0.5 * 10;

  legend_svg.select("text")
            .attr("x", (width + margin.left + margin.right - legendLabel)/2 + margin.left)

  // Build X scales and axis:
  var x = d3v4.scaleBand()
    .range([ 0, width ])
    .domain(cols)
    .padding(0.01);

  var xColor = d3v4.scaleBand()
    .range([ 0, width ])
    .domain(colorCols)
    .padding(0.01);

  svg.append("g")
    .attr("transform", "translate(0," + height + ")")
    .attr('id', 'Xaxis')
    .call(d3v4.axisBottom(x).tickSizeOuter(0));

  // Build X scales and axis:
  var y = d3v4.scaleBand()
    .range([ height, 0 ])
    .domain(sorted_rows)
    .padding(0.01);

  color_svg.append("g")
    .attr('id', 'Yaxis')
    .call(d3v4.axisLeft(y).tickSize(0));

  // Build color scale
  var myColor = d3v4.scaleLinear()
    .range(["white", "black"])
    .domain([1,highest_value])

  //Read the data
    color_svg.selectAll()
        .data(colorData, function(d) {return d.row+':'+d.col;})
        .enter()
        .append("rect")
        .attr("x", function(d) { return x(d.col) })
        .attr("y", function(d) { return y(d.row) })
        .attr("width", "15px" )
        .attr("height", y.bandwidth() )
        .style("fill", function(d) { return d.value} )

  //Read the data
    svg.selectAll()
        .data(chartData, function(d) {return d.row+':'+d.col;})
        .enter()
        .append("rect")
        .attr("x", function(d) { return x(d.col) })
        .attr("y", function(d) { return y(d.row) })
        .attr("width", x.bandwidth() )
        .attr("height", y.bandwidth() )
        .style("fill", function(d) { return myColor(d.value)} )

    svg.select('#Xaxis')
       .attr('text-anchor', 'start')

    d3v4.selectAll("#Yaxis>.tick>text")
        .each(function(d, i){
          d3.select(this).style("font-size","1.3em");
        });

    d3v4.selectAll("#Xaxis>.tick>text")
        .each(function(d, i){
          d3.select(this).style("font-size","1.3em");
        });

    var count = 1;
    ticks = svg.select('#Xaxis').selectAll('.tick');
    ticks.each(function(d) {
      value = oddOrEven(count);
      text = d3.select(this).select('text');
      textSize = Math.floor(text[0][0].getBBox().width*1.05 + 0.5 * 10);
      text.attr("transform", null);
      text.attr("x", "-" + textSize/2 + "px");
      if(value === "even"){
        // text.attr("x", "-15");
        text.attr("y", "30");
      }
      count = count + 1;
    });

    yTicks = color_svg.select('#Yaxis').selectAll('.tick');
    yTicks.each(function(d) {
      yText = d3.select(this).select('text')[0];
      for (var i=0; i < html_converter.length; i++){
        if (html_converter[i]['row'] === yText[0].innerHTML.replaceAll("&lt;", "<").replaceAll("&gt;", ">").replaceAll("amp;","")){
          yText[0].innerHTML = html_converter[i]['html'];
          labelSize = yText[0].getBBox().width*1.05 + 0.5 * 10
          svg_home.attr("width", width + margin.left + margin.right + labelSize);
          svg.attr("transform", "translate(" + (labelSize+40) + "," + (margin.top*1.5) + ")");
          color_svg.attr("transform", "translate(" + (labelSize+20) + "," + (margin.top*1.5) + ")");
        }
      }
    });

}


function draw_model_scores(location, element_id, score, startValue, endValue, definedValue, gradientChange) {

    var margin = {top: 30, right: 30, bottom: 30, left: 30};
    var width = 200;
    var height = 50;

    // append the svg object to the body of the page
    var svg_home = d3v4.select("#" + location)
                .append("svg")
                .attr("width", width + margin.left + margin.left + margin.right)
                .attr("height", (margin.top))
                .attr("transform", "translate(-30,0)")
                .attr("id", element_id);

    var greyscale = [ 'rgb(255,0,0)', 'rgb(255, 255, 255)', 'rgb(0,255,0)' ];

    var first_legend = svg_home.append('defs')
                         .append('linearGradient')
                         .attr('id', 'grad_' + score)
                         .attr('x1', '0%')
                         .attr('x2', '100%')
                         .attr('y1', '0%')
                         .attr('y2', '0%');

     var stopOffset = (gradientChange - startValue) / (endValue - startValue) * 100;

    first_legend.append('stop')
           .style('stop-color', greyscale[0])
           .attr('offset', '0%');
    first_legend.append('stop')
           .style('stop-color', greyscale[1])
           .attr('offset', stopOffset + '%');
    first_legend.append('stop')
           .style('stop-color', greyscale[1])
           .attr('offset', stopOffset + '%');
    first_legend.append('stop')
           .style('stop-color', greyscale[2])
           .attr('offset', '100%');


    var first_legend_svg = svg_home.append("g")
                      .attr("transform", "translate(0,-15)");

    first_legend_svg.append("text")
              .attr('y', 5)
              .attr('x', 40)
              .style("font", "14px sans-serif")
              .style("font-weight", "bold");

    first_legend_svg.append("text")
              .attr('x', 30)
              .attr('y', 33)
              .style("font", "14px sans-serif")
              .style("font-weight", "bold")
              .text(startValue);

    first_legend_svg.append('rect')
              .attr('x', 40)
              .attr('y', 20) // Adjust the y-coordinate to add space between gradient and text box
              .attr('width', (width*0.9))
              .attr('height', 20)
              .style('fill', 'url(#grad_' + score + ')');

    first_legend_svg.append("text")
              .attr('x', width + 20)
              .attr('y', 33)
              .style("font", "14px sans-serif")
              .style("font-weight", "bold")
              .text(endValue);

    legendLabel = first_legend_svg.select('text')['_groups'][0][0].getBBox().width*1.05 + 0.5 * 10;

    first_legend_svg.select("text")
              .attr("x", (width + margin.left + margin.right - legendLabel)/2);

    // Draw the black vertical line at the defined value
    var lineX = (definedValue - startValue) / (endValue - startValue) * (width * 0.9) + 40;
    svg_home.append("line")
            .attr("x1", lineX)
            .attr("y1", height - 45) // Adjust the y-coordinate of the top arrowhead 90
            .attr("x2", lineX)
            .attr("y2", height - 25) // Adjust the y-coordinate of the bottom arrowhead 125
            .attr("stroke", "black")
            .attr("stroke-width", 2);

    // Add arrowheads
    svg_home.append("path")
            .attr("d", "M " + lineX + " " + (height - 45) + " L " + (lineX - 5) + " " + (height - 50) + " L " + (lineX + 5) + " " + (height - 50) + " Z")
            .attr("fill", "black");

    svg_home.append("path")
            .attr("d", "M " + lineX + " " + (height - 25) + " L " + (lineX - 5) + " " + (height - 20) + " L " + (lineX + 5) + " " + (height - 20) + " Z")
            .attr("fill", "black");

}

// This function generates two concentric circles made of aminoAcids
// then creates lines connecting inner and outer aminoacids, given interactions between them
function draw_interactions_in_circles(location, interactions, inner_data, outer_data) {

  // D3 select the SVG
  const svg2 = d3.select("#" + location);

  // Types of interactions and their colors
  const interactionTypes = {
    'Aromatic': 'red',
    'Hydrophobic': 'blue',
    'Ionic': 'green',
    'Polar': 'purple',
    'Van der waals': 'orange'
  };

  const segmentColors = {
    'TM1': 'red',
    'TM2': 'blue',
    'TM3': 'green',
    'TM4': 'purple',
    'TM5': 'orange',
    'TM6': 'brown',
    'TM7': 'pink',
    'IL1': 'grey',
    'IL2': 'cyan',
    'IL3': 'yellow'
  };

  const aminoAcids = [
    { aminoAcid: "Alanine", singleLetter: "A" },
    { aminoAcid: "Arginine", singleLetter: "R" },
    { aminoAcid: "Asparagine", singleLetter: "N" },
    { aminoAcid: "Aspartic acid", singleLetter: "D" },
    { aminoAcid: "Cysteine", singleLetter: "C" },
    { aminoAcid: "Glutamine", singleLetter: "Q" },
    { aminoAcid: "Glutamic acid", singleLetter: "E" },
    { aminoAcid: "Glycine", singleLetter: "G" },
    { aminoAcid: "Histidine", singleLetter: "H" },
    { aminoAcid: "Isoleucine", singleLetter: "I" },
    { aminoAcid: "Leucine", singleLetter: "L" },
    { aminoAcid: "Lysine", singleLetter: "K" },
    { aminoAcid: "Methionine", singleLetter: "M" },
    { aminoAcid: "Phenylalanine", singleLetter: "F" },
    { aminoAcid: "Proline", singleLetter: "P" },
    { aminoAcid: "Serine", singleLetter: "S" },
    { aminoAcid: "Threonine", singleLetter: "T" },
    { aminoAcid: "Tryptophan", singleLetter: "W" },
    { aminoAcid: "Tyrosine", singleLetter: "Y" },
    { aminoAcid: "Valine", singleLetter: "V" }
  ];

  const presetColors = {'D': ['#E60A0A', '#FDFF7B'],'E': ['#E60A0A', '#FDFF7B'],
                        'K': ['#145AFF', '#FDFF7B'],'R': ['#145AFF', '#FDFF7B'],
                        'S': ['#A70CC6', '#FDFF7B'],'T': ['#A70CC6', '#FDFF7B'],
                        'N': ['#A70CC6', '#FDFF7B'],'Q': ['#A70CC6', '#FDFF7B'],
                        'V': ['#E6E600', '#000000'],'L': ['#E6E600', '#000000'],
                        'I': ['#E6E600', '#000000'],'A': ['#E6E600', '#000000'],
                        'M': ['#E6E600', '#000000'],'F': ['#18FF0B', '#000000'],
                        'Y': ['#18FF0B', '#000000'],'W': ['#0BCF00', '#000000'],
                        'H': ['#0093DD', '#000000'],'P': ['#CC0099', '#FDFF7B'],
                        'C': ['#B2B548', '#000000'],'G': ['#FF00F2', '#000000'],
                        '-': ['#FFFFFF', '#000000'],'+': ['#FFFFFF', '#000000']};


  // // This is an example of the structure of interactions data
  // const mockupInteractions = [
  //   { innerIndex: 0, outerIndex: 30, type: 'Aromatic' },
  //   { innerIndex: 1, outerIndex: 31, type: 'Hydrophobic' },
  //   { innerIndex: 5, outerIndex: 35, type: 'Ionic' },
  //   { innerIndex: 7, outerIndex: 37, type: 'Polar' },
  //   // [...]
  // ];
  //
  // // This is an example of the structure of outer beads data
  // const mockupOuterBeads = [
  //   { aminoAcid: 'Y', segment: "TM3", generic_number: '1.33x55', interaction: 'Yes'},
  //   //[...]
  // ];
  // // This is an example of the structure of inner beads data
  // const mockupOuterBeads = [
  //   { aminoAcid: 'Y', generic_number: '1.33x55', interaction: 'Yes'},
  //   //[...]
  // ];

  // Generate random arrays for inner and outer circles
  const beadInfo = {
    innerCircle: inner_data,
    outerCircle: outer_data
  };

  function sanitizeClassName(name) {
    return name.toLowerCase().replace(/\s+/g, '-');
  }

  const tooltip = d3.select("body").append("div")
    .attr("class", "tooltip")
    .style("opacity", 0);

  // Function to create a bead
  function createBead(cx, cy, aminoAcid, segment, index, circleType) {

    const bead = svg2.append("circle")
      .attr("cx", cx)
      .attr("cy", cy)
      .attr("r", 10)
      .attr("fill", "white")
      .attr("stroke", "black")
      .attr("stroke-width", 1)
      .attr("data-index", index)
      .attr("data-circle", circleType)
      .attr("data-aa", aminoAcid)
      .attr("data-segment", segment)
      .classed(segment, true)
      .on("mouseover", function(d) {
        const event = window.event ? window.event : d3.event;
        tooltip.transition()
          .duration(200)
          .style("opacity", .9);
        tooltip.html("Amino Acid: " + aminoAcid + "<br/>" +
                    "Segment: " + segment + "<br/>" +
                    "Generic Number: " + index + "<br/>" +
                    "Interaction: TBD") // Assuming you will provide interaction later
          .style("left", (event.pageX + 5) + "px")
          .style("top", (event.pageY - 28) + "px");
      })
      .on("mouseout", function(d) {
        const event = d3.event; // get the event object from d3.event
        tooltip.transition()
          .duration(500)
          .style("opacity", 0);
      });

    // Add text inside the circle to represent the amino acid
    svg2.append("text")
      .attr("x", cx)
      .attr("y", cy)
      .attr("dy", 5)
      .attr("text-anchor", "middle")
      .text(aminoAcid);
  }

  // Function to create a circle of beads
  function createCircle(cx, cy, radius, beads, circleType) {
    const numBeads = beads.length;
    beads.forEach((bead, index) => {
      const angle = (index / numBeads) * 2 * Math.PI;
      const x = cx + radius * Math.cos(angle);
      const y = cy + radius * Math.sin(angle);
      createBead(x, y, bead.aminoAcid, bead.segment, index, circleType);
    });
  }

  // Create connection
  function createConnection(innerIndex, outerIndex, type) {
    const innerBead = d3.select(svg2.selectAll("circle").nodes()[innerIndex]);
    const outerBead = d3.select(svg2.selectAll("circle").nodes()[outerIndex]);

    const x1 = +innerBead.attr("cx");
    const y1 = +innerBead.attr("cy");
    const x2 = +outerBead.attr("cx");
    const y2 = +outerBead.attr("cy");

    // Calculate angle between two beads
    const angle = Math.atan2(y2 - y1, x2 - x1);

    // Calculate new coordinates for the line to end at the circle border instead of center
    const new_x1 = x1 + 10 * Math.cos(angle);
    const new_y1 = y1 + 10 * Math.sin(angle);
    const new_x2 = x2 - 10 * Math.cos(angle);
    const new_y2 = y2 - 10 * Math.sin(angle);

    const sanitizedType = sanitizeClassName(type);

    svg2.append("line")
      .attr("x1", new_x1)
      .attr("y1", new_y1)
      .attr("x2", new_x2)
      .attr("y2", new_y2)
      .attr("stroke", interactionTypes[type])
      .classed(sanitizedType, true);
  }

  // Create circles of beads
  createCircle(200, 200, 70, beadInfo.innerCircle, "inner");
  createCircle(200, 200, 150, beadInfo.outerCircle, "outer");

  // Generate legend
  const interactionLegend = svg2.append("g")
    .attr("transform", "translate(390, 30)");

    Object.entries(interactionTypes).forEach(([type, color], i) => {
      const legendItem = interactionLegend.append("g")
        .attr("transform", `translate(0, ${i * 25})`);  // Increased spacing between items

      legendItem.append("rect")
        .attr("width", 20)
        .attr("height", 20)
        .attr("fill", color)
        .attr("stroke", "black")  // Black border
        .attr("stroke-width", 1)  // Border width
        .on("click", function() {
          if(d3.event) {
            d3.event.stopPropagation(); // Prevent the SVG click event when clicking on the legend
          }
          const sanitizedType = sanitizeClassName(type);
          d3.selectAll("line").attr("stroke-opacity", 0.1); // Make other lines faded
          d3.selectAll(`.${sanitizedType}`).attr("stroke-opacity", 1); // Highlight the lines of the clicked type
        });

      legendItem.append("text")
        .attr("x", 30)  // Adjusted the x position so that text doesn't overlap with rectangle
        .attr("y", 15)  // Center text in the y-direction
        .text(type);

      // Header for the interaction legend
      interactionLegend.append("text")
        .attr("x", 0)
        .attr("y", -10)
        .text("Interactions")
        .attr("font-family", "Arial")
        .attr("font-size", "14px")
        .attr("font-weight", "bold")
        .attr("fill", "black");

    });

    // Generate legend for segments
    const segmentLegend = svg2.append("g")
      .attr("transform", "translate(390, 200)");  // Adjust this to fit your needs

    // Header for the segment legend
    segmentLegend.append("text")
      .attr("x", 0)
      .attr("y", -10)
      .text("Segments")
      .attr("font-family", "Arial")
      .attr("font-size", "14px")
      .attr("font-weight", "bold")
      .attr("fill", "black");

  // Filter segments into TMs and ILs
  const tmSegments = Object.entries(segmentColors).filter(([type]) => type.startsWith('TM'));
  const ilSegments = Object.entries(segmentColors).filter(([type]) => type.startsWith('IL'));

  const createLegendColumn = (legendData, xOffset) => {
    legendData.forEach(([type, color], i) => {
      const legendItem = segmentLegend.append("g")
        .attr("transform", `translate(${xOffset}, ${i * 25})`);

      legendItem.append("rect")
        .attr("width", 20)
        .attr("height", 20)
        .attr("fill", color)
        .attr("stroke", "black")
        .attr("stroke-width", 1)
        .on("click", function() {
          if (d3.event) {
            d3.event.stopPropagation(); // Prevent the SVG click event when clicking on the legend
          }

          // Reset all beads to white
          d3.selectAll("circle").attr("fill", "white");

          // Color beads of the clicked segment with the color parameter
          d3.selectAll(`.${type}`).attr("fill", function(d) {
            return d3.select(this).attr("fill") === "white" ? color : "white";
          });
        });

      legendItem.append("text")
        .attr("x", 30)
        .attr("y", 15)
        .text(type);
    });
  }

  // Create TM and IL legend columns
  createLegendColumn(tmSegments, 0);  // TM segments column starts at x=0
  createLegendColumn(ilSegments, 75);  // IL segments column starts at x=150


  interactions.forEach(({ innerIndex, outerIndex, type }) => createConnection(innerIndex, outerIndex, type));

  function applyPresetColors(svg) {
    svg.selectAll("circle")
      .each(function(d, i) {
        const aa = d3.select(this).attr("data-aa"); // Assuming 'data-aa' stores the amino acid type
        d3.select(this).attr("fill", presetColors[aa][0]);

        // Assuming the text element immediately follows the circle
        d3.select(this.nextElementSibling).attr("fill", presetColors[aa][1]);
      });
  }

  document.getElementById("propertiesButton").addEventListener("click", function() {
    applyPresetColors(svg2); // Assuming your SVG is stored in the svg2 variable
  });

  // Function to reset colors
  function resetColors(svg) {
    svg.selectAll("circle")
      .attr("fill", "white")
      .attr("stroke", "black");

    // Assuming the text elements are immediate siblings of the circle elements
    svg.selectAll("text")
      .attr("fill", "black");
  }

  // Add event listener for reset button
  document.getElementById("resetButton").addEventListener("click", function() {
    resetColors(svg2);  // Reset colors of circles
    d3.selectAll("line").attr("stroke-opacity", 1);  // Reset opacity of lines
    d3.selectAll("circle").attr("fill-opacity", 1);  // Reset opacity of lines
  });

  function applySegmentColors(svg) {
    svg.selectAll("circle[data-circle='outer']")
      .each(function(d, i) {
        const segment = d3.select(this).attr("data-segment");
        if (segment) {
          d3.select(this).attr("fill", segmentColors[segment]);
        }
      });
  }

  document.getElementById("segmentButton").addEventListener("click", function() {
    applySegmentColors(svg2);
  });

}

function renderTree(data) {
    var tree = data["tree"]; // contains tree in Newick format
    // Annotations: state, name, family, ligand type, class
    var annotations = data["annotations"];
    var r = 1200 / 2;
    var spacing = 240;
    var innerRadius = r - spacing // change inner radius of tree with this argument
    var names = 0; // indexing for all nodes


    var cluster = d3.layout.cluster()
        .size([360, innerRadius])
        .sort(function(a, b) {
            return (a.value - b.value) || d3.ascending(a.length, b.length);
        })
        .value(function(d) {
            return d.length;
        })
        .children(function(d) {
            return d.branchset;
        })
        .separation(function(a, b) {
            return 1;
        });

    function project(d) {
        var r = d.y,
            a = (d.x - 90) / 180 * Math.PI;
        return [r * Math.cos(a), r * Math.sin(a)];
    }

    function cross(a, b) {
        return a[0] * b[1] - a[1] * b[0];
    }

    function dot(a, b) {
        return a[0] * b[0] + a[1] * b[1];
    }

    function step(d) {
        var s = project(d.source),
            m = project({
                x: d.target.x,
                y: d.source.y
            }),
            t = project(d.target),
            r = d.source.y,
            sweep = d.target.x > d.source.x ? 1 : 0;
        return (
            "M" + s[0] + "," + s[1] +
            "A" + r + "," + r + " 0 0," + sweep + " " + m[0] + "," + m[1] +
            "L" + t[0] + "," + t[1]);
    }

    var rect_dim = (r - 90) * 2;
    var svg_dim = (r - 90) * 2;

    var wrap = d3.select("#clustering-tree").append("svg")
        .attr("width", svg_dim)
        .attr("height", svg_dim)
        .style("-webkit-backface-visibility", "hidden");

    // Catch mouse events in Safari.
    wrap.append("rect")
        .attr("width", rect_dim)
        .attr("height", rect_dim)
        .attr("fill", "none");

    var vis = wrap.append("g")
        .attr("transform", "translate(" + rect_dim / 2 + "," + rect_dim / 2 + ")");

    var start = null,
        rotate = 0,
        div = document.getElementById("clustering-tree");

    //Branch length function
    function phylo(n, offset) {
        if (n.length != null) offset += n.length/10;
        n.y = offset;
        if (n.children)
            n.children.forEach(function(n) {
                phylo(n, offset);
            });
    }

    var nodes = cluster.nodes(newick.parse(tree));
    var reg = new RegExp('^[0-9]+[\.]?[0-9]*$');
    nodes.forEach(function(n) {
        // HACK: the internal nodes names are now cluster scoring values
        if (n.name != "" && reg.test(n.name)) {
            n.score = n.name
            n.name = names.toString();
            names++;
        }

        // ORIGINAL
        /*if (n.name == "") {
            n.name = names.toString();
            names++;
        }*/
    });

    // Utilized to calculate actual branch lengths
    phylo(nodes[0], 0);

    var link = vis.selectAll("path.link")
        .data(cluster.links(nodes))
        .enter().append("path")
        .attr("class", "link")
        .attr("d", step);


    //D3 selection to distinguish between inner and leaf nodes
    var node = vis.selectAll("g.node")
        .data(nodes.filter(function(n) {
            return n.x !== undefined;
        }))
        .enter().append("g")
        .attr("class", function(n) {
            if (n.children) {
                return "inner node";
            } else {
                //return "leaf node";
                return ('X' + n.name + ' terminal-node');
            }
        })
        .attr("transform", function(d) {
            return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")";
        })

    // Add terminal node with coloring
    vis.selectAll("g.terminal-node").append("circle")
        .attr("r", 5)
        .style("fill", function(n) {
          // color based on activity
          if (annotations[n.name]) {
              switch(annotations[n.name][0]){
                  case "active":
                      return "#F00";
                  case "inactive":
                      return "#00F";
                  case "intermediate":
                      return "#F80";
                  default:
                      return "#888";
              }
          } else {
            return "#000";
          }
        });

    // Add terminal node with coloring
    vis.selectAll("g.terminal-node").append("circle")
        .attr("r", 5)
        .attr("transform", "translate(10, 0)")
        .style("fill", function(n) {
          // color based on activity
          if (annotations[n.name]) {
              if (annotations[n.name][7].length > 0){
                  switch(annotations[n.name][7][0]){
                      case "agonist":
                      case "agonist-partial":
                      case "pam":
                          return "#F00";
                      case "antagonist":
                      case "inverse-agonist":
                      case "nam":
                          return "#00F";
                      default:
                          return "#888";
                  }
              } else {
                  return "#888";
              }
          }
        });


    vis.selectAll('g.inner.node')
        .append("circle")
        .attr("r", 5)
        .style("fill", function(n){ if (isNaN(n.score)) return "#FFFFFF"; var score = n.score; if (score>1) score=1; if (score>0){ return shadeColor2("#AAAAAA", 100-(score*80-20))} else { return "#FFAAAA"} })
        .attr('data-length', function(n){ return Math.round(n.y*10) })
        .attr('data-score', function(n){ return n.score })
        .on("mouseover", function(d,i) {
            var distance = d3.select(this).attr("data-length")
            var score = d3.select(this).attr("data-score")
            var label = "Root node"
            if (distance > 0)
                label = "Distance from root: " + distance + "<br/>Silhouette coefficient: " + score;

            tooltip
                .style("background-color", shadeColor("#999999", 50))
                .style("border-color", "#999999")
                .style("display", "block")
                .style("opacity", .9);
            tooltip.html(label)
                .style("left", (d3.event.pageX + 5) + "px")
                .style("top", (d3.event.pageY - 28) + "px");
            })
        .on("mouseout", function(d) {
            tooltip.style("display", "none")
        });

    // Adding annotations
    // Annotations: state, name, family, ligand type, class
    var categoryName = ["State", "Name", "Family", "Ligand family", "GPCR Class", "Method"]
    var spacer = 10;
    var colorscheme = []
    colorscheme[0] = ['#008000','#797f98','#7a97b2','#75afc9','#68c9dc','#50e4ee','#00ffff']
    colorscheme[1] = ["#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf","#999999"]
    colorscheme[2] = ['#ffd700','#dfda5d','#d3d772','#cad381','#c2ce8e','#bcc998','#b7c4a0','#b4c0a6','#b0b9ad','#adb4b2','#abaeb7','#a9a9bb','#a8a2bf','#a89cc1','#a895c4','#a98fc6','#aa87c7','#ad7fc7','#b077c7','#f44191']
    colorscheme[3] = ['#8b0000','#960110','#9e051b','#a80c25','#b1142d','#b91c35','#c1253d','#c92e43','#d03649','#d7404e','#dd4852','#e35256','#e85b59','#ed655d','#f26f60','#f67863','#fa8266','#fd8c69','#ff956d','#ffa072','#ffac77','#ffb57e','#ffbf86','#ffc98f','#ffd399','#ffdca5','#ffe5b2','#ffedbf','#fff6cf','#ffffe0']

    var categories = [];
    var tooltip = d3.select("body").append("div")
                  .attr("class", "tooltip")
                  .style("opacity", 0);
    for  (i=6; i>2; i--) {
        // store categories for legend
        categories[6-i] = []
        var ref = 6-i;
        for (var pdb in annotations){
            // assign color
            if (!categories[ref].includes(annotations[pdb][i]))
              categories[ref].push(annotations[pdb][i])
            colorIndex = categories[ref].indexOf(annotations[pdb][i]);
            color = colorscheme[ref][colorIndex];

            // create tooltip
            var popoverTable = "<h3>" + annotations[pdb][1] + " ("+ pdb + ")" + "</h3>"  + categoryName[i] + ": " + annotations[pdb][i];

            // create annotation
            var rect = vis.selectAll('g.X'+pdb)
              //.append("circle")
              //.attr("r", 3.25)
              .append("rect")
              .attr("width", spacer-1)
              .attr("height", spacer-1)
              .style("fill", color)
              .attr("transform", "translate(" + (110 + Math.abs(i-6)*spacer) + ", " + -1 * (spacer-1)/2 + ")")
              .attr('data-content', popoverTable)
              .attr('data-color', color)
              .on("mouseover", function(d,i) {
                  tooltip
                      .style("background-color", shadeColor(d3.select(this).attr("data-color"), 50))
                      .style("border-color", d3.select(this).attr("data-color"))
                      .style("display", "block")
                      .style("opacity", .9);
                  tooltip.html(d3.select(this).attr("data-content"))
                      .style("left", (d3.event.pageX + 5) + "px")
                      .style("top", (d3.event.pageY - 28) + "px");
                  })
              .on("mouseout", function(d) {
                  tooltip.style("display", "none")
              });
        }
    }

    var label = vis.selectAll("text")
        .data(nodes.filter(function(d) {
            return d.x !== undefined && !d.children;
        }))
        .enter().append("text")
        .attr("dy", ".31em")
        .attr("text-anchor", function(d) {
            return d.x < 180 ? "start" : "end";
        })
        .attr("transform", function(d) {
            return "rotate(" + (d.x - 90) + ")translate(" + (d.y + 18) + ")rotate(" + (d.x < 180 ? 0 : 180) + ")";
        })
        .text(function(d) {
            // add receptor name
            if (annotations[d.name]) {
                var name = annotations[d.name][1].split("_")[0].toUpperCase()
                if (d.x < 180) {
                  return d.name.replace(/_/g, ' ') + ' (' + name + ')';
                } else {
                  return '(' + name + ') ' + d.name.replace(/_/g, ' ');
                }
            }
        });


    //click count to change clade color selection and unhighlight subtree nodes
    var click_count = 0;

    //On click event to get subtree of clicked node
    var firstTreeClick = true;
    node.on("click", function(d) {
        if (firstTreeClick){
          firstTreeClick = false;
          $("#CN-button").removeClass("disabled");
          $("#DN-button").removeClass("disabled");
        }
        doubleclade(d);
    });

    //Highlight clade when mouse is over node
    node.on("mouseover", function(h) {
        traverse(h);
    });

    //Unhighlight clade when mouse is out of node
    node.on("mouseout", function(u) {
        vis.selectAll("path.selected")
            .classed("selected", false);
        vis.selectAll("text").style("font-weight", "normal");
    });

    //iterate over links for (source,target) info to highlight subtree path
    function highlightlink(src, tgt) {
        var link = d3.selectAll("path.link")[0].filter(function(d) {
            var flag = (d3.select(d).data()[0].source.name == src && d3.select(d).data()[0].target.name == tgt);
            return flag;
        });

        d3.selectAll(link)
            .classed("selected", true);
    }

    //iterate over labels for leaf nodes info to highlight labels
    function highlightlabel(tgt) {
        var highlight_labels = vis.selectAll("text")[0].filter(function(t) {
            var th = (d3.select(t).data()[0].name == tgt);
            return th;
        });

        d3.selectAll(highlight_labels)
            .style("font-weight", "bold");
    }

    //Recursive function to highlight all links of a subtree
    function traverse(node) {
        if (node.children) {
            node.children.forEach(function(d) {
                highlightlabel(d.name);
                highlightlink(node.name, d.name);
                traverse(d);
            });
        }
    }

    //Recursive function to traverse a subtree and deselect all nodes
    function deselect(node) {
        node.isSelected = false;
        if (node.children) {
            node.children.forEach(function(d) {
                //Select all children
                d.isSelected = false;
                deselect(d);
            });
        }
    }


    //Recursive function to traverse a node and its subtree
    function selectSubtree(node, subtree) {
        node.isSelected = true;
        if (node.children) {
            node.children.forEach(function(d) {
                d.isSelected = true;
                highlightlabel(d.name)

                if (d.name.length > 3) {
                    subtree.push(d.name);
                }

                selectSubtree(d, subtree);
            });
        }

        return subtree;
    }

    // Based on https://stackoverflow.com/questions/5560248
    function shadeColor(color, percent) {
        var R = parseInt(color.substring(1,3),16);
        var G = parseInt(color.substring(3,5),16);
        var B = parseInt(color.substring(5,7),16);

        //R = parseInt(R * (100 + percent)/100);
        //G = parseInt(G * (100 + percent)/100);
        //B = parseInt(B * (100 + percent)/100);
        //R = (R<255)?R:255;
        //G = (G<255)?G:255;
        //B = (B<255)?B:255;
        R = parseInt(R + (255-R)*percent/100);
        G = parseInt(G + (255-R)*percent/100);
        B = parseInt(B + (255-R)*percent/100);

        var RR = ((R.toString(16).length==1)?"0"+R.toString(16):R.toString(16));
        var GG = ((G.toString(16).length==1)?"0"+G.toString(16):G.toString(16));
        var BB = ((B.toString(16).length==1)?"0"+B.toString(16):B.toString(16));

        return "#"+RR+GG+BB;
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

    //Double clade selection
    function doubleclade(node) {
        var groupId = click_count % 2;
        var subtreeArray = selectSubtree(node, []);

        // Reset previous selection - remove class
        vis.selectAll('circle.group' + groupId)
            .classed('group' + groupId, false)
        vis.selectAll('path.group' + groupId)
            .classed('group' + groupId, false)

        var subnode = vis.selectAll('g.inner.node')
            .attr("id", function(d) {
                if (d.isSelected == true) {
                    return 'sub' + groupId;
                }
            });

        // Color new selection - add class
        vis.selectAll('#sub' + groupId).selectAll("circle")
            .attr("class", "group" + groupId);
        vis.selectAll('path.selected')
            .attr("class", "link group" + groupId);

        // Deselect
        vis.selectAll('g#sub' + groupId +'.inner.node')
            .attr("id", "deselect");
        vis.selectAll("text").style("font-weight", "normal");

        makeUL(subtreeArray, click_count % 2);
        deselect(node);
        click_count++;
    }

    //Fill text area with subtree data
    function makeUL(descendents, group) {
        var textarea = document.getElementById('input-targets-'+group);
        textarea.value = descendents.join("\n");
    }
}

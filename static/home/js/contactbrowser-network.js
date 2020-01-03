function createNetworkPlot(raw_data,original_width, inputGraph, containerSelector, segment_view = true) {
// https://github.com/d3/d3-3.x-api-reference/blob/master/Force-Layout.md
 //https://archive.nytimes.com/www.nytimes.com/interactive/2013/02/20/movies/among-the-oscar-contenders-a-host-of-connections.html

    circle_size = 20;
    max_link_size = 15; 

    // https://bl.ocks.org/mbostock/4600693
    // https://stackoverflow.com/questions/13455510/curved-line-on-d3-force-directed-tree
    // http://bl.ocks.org/mbostock/1153292
    var curved_links = true;

    // var dr = 4,      // default point radius
    // off = 15,    // cluster hull offset
    // expand = {}, // expanded clusters
    // data, net, force, hullg, hull, linkg, link, nodeg, node;

    var w = original_width;
    var h = w;
    var height = w;
    var width;

    var selected_single_cluster = false

    var new_data;
    var plot_specified_filtered = filtered_gn_pairs;
    var cluster_groups = filtered_cluster_groups;
    function prepare_data(single_cluster = false) {
        selected_single_cluster = single_cluster;

        new_data = { "nodes": [], "links": [] };
        var track_gns = [];
        // cluster_groups = [];
        // console.log('single_cluster',single_cluster);
        // // // Populate matrix for interactions between segments
        // 
        // $.each(raw_data['interactions'], function (i, v) {
        //     if (!plot_specified_filtered.includes(i)) return;
        //     gns = separatePair(i);

        //     test1 = cluster_groups.filter(l => l.includes(gns[0]));
        //     test2 = cluster_groups.filter(l => l.includes(gns[1]));
        //     if (!test1.length && !test2.length) {
        //         cluster_groups.push([gns[0], gns[1]]);
        //     } else if (test1.length && !test2.length) {
        //         i1 = cluster_groups.indexOf(test1[0])
        //         cluster_groups[i1].push(gns[1]);
        //     } else if (!test1.length && test2.length) {
        //         i2 = cluster_groups.indexOf(test2[0])
        //         cluster_groups[i2].push(gns[0]);
        //     } else if (test1.length && test2.length) {
        //         i1 = cluster_groups.indexOf(test1[0])
        //         i2 = cluster_groups.indexOf(test2[0])
        //         //i1 = cluster_groups.indexOfForArrays(test1[0]);
        //         if (i1!=i2) {
        //             cluster_groups[i1] = test1[0].concat(test2[0])
        //             cluster_groups.splice(i2, 1);
        //         }
        //     }

        //     // if (seg1 != seg2) {
        //     //     new_data["links"].push({ "source": gns[0], "target": gns[1], "value": 1 })
        //     // }
        // });

        $.each(raw_data['interactions'], function (i, v) {
            if (!plot_specified_filtered.includes(i)) return;
            gns = separatePair(i);
            seg1 = raw_data['segment_map'][gns[0]];
            seg2 = raw_data['segment_map'][gns[1]];

            cg = cluster_groups.filter(l => l.includes(gns[0]));
            cg_index1 = cluster_groups.indexOf(cg[0]);
            cg = cluster_groups.filter(l => l.includes(gns[1]));
            cg_index2 = cluster_groups.indexOf(cg[0]);

            if (single_cluster!==false && (cg_index1!=single_cluster || cg_index2!=single_cluster)) return;

            if (!track_gns.includes(gns[0])) {
                new_data["nodes"].push({ "name": gns[0], "group": seg1, "group2":cg_index1 })
                track_gns.push(gns[0]);
            }
            if (!track_gns.includes(gns[1])) {
                new_data["nodes"].push({ "name": gns[1], "group": seg2, "group2":cg_index2 })
                track_gns.push(gns[1]);
            }

            source = new_data["nodes"].filter(obj => { return obj.name === gns[0] })[0];
            target = new_data["nodes"].filter(obj => { return obj.name === gns[1] })[0];
            new_data["links"].push({ "source": source, "target": target, "value": 1 })
        })


        console.log(cluster_groups);
            // console.log(graph);
        // new_data['nodes'].forEach(function (n) {
        //     cg = cluster_groups.filter(l => l.includes(n.name));
        //     cg_index = cluster_groups.indexOf(cg[0]);
        //     console.log(n, n.size, cg_index, cg);
        //     n.group2 = cg_index;
        // });
        // console.log(new_data['nodes']);

        if (!segment_view) {
            // resize by number of nodes if not segments
            console.log(new_data["nodes"].length, "nodes");
            width = original_width*Math.sqrt(new_data["nodes"].length/8);
        } else {
            width = original_width*0.8;
        }
        w = width;
        h = w;
        height = w;

    }
    prepare_data();




    var expand = {},net, force;

    d3.select(containerSelector).style("position","relative");

    div = d3.select(containerSelector).insert("div")
        .attr("class", "network")
        .style("width", "100%")
        .style("-webkit-backface-visibility", "hidden");

    init();
    
    var path, link, node;
    function init() {
        console.log('init',div, containerSelector);
        div = d3.select(containerSelector).select("div");
        div.select("svg").remove();
        // if (force) force.stop();
        svg = div.append("svg:svg")
            .attr("viewBox", "0 0 " + w + " " + h)
            .attr("width", "100%")
            .attr("style", "height: 500px");
        
        // slightly hide the initial load
        svg.attr("opacity", 1e-6)
            .transition()
            .duration(3000)
            .attr("opacity", 1);

        var graph = new_data;

        if (segment_view) {
            // Group nodes into their "groups" (segments)
            graph = network(new_data, net, getGroup, expand);

            max_size = Math.max.apply(Math, graph['nodes'].map(function(v) {
                return v.size;
            }));

            max_link = Math.max.apply(Math, graph['links'].map(function(v) {
                return v.size;
            }));

            // normalize
            graph['links'].forEach(function (n) {
                n.size = Math.round(max_link_size*n.size/max_link);
            });
        }
        if (segment_view) {
            force = d3.layout.force()
                .size([w, h])
                .gravity(0.05)
                // .charge(-1000)
                .charge(function (d, i) {
                    return ((d.weight-1) * -800) - 50;
                })
                .linkStrength(1)
                .linkDistance(100)
                .friction(0.5) 
                // .friction(0.5)
                .on("tick", tick);
        } else {
            // var force = d3.layout.force()
            //     .size([w, h])
            //     .gravity(0.15)
            //     .charge(-900)
            //     .linkStrength(1)
            //     .linkDistance(50)
            //     // .friction(0.5)
            //     .on("tick", tick);
            // https://unpkg.com/force-in-a-box/dist/forceInABox.js
            force = d3.layout.forceInABox()
                .charge(-100)
                // .charge(function (d, i) {
                //     return ((d.weight-1) * -400) - 100;
                // })
                .linkDistance(20)
                .linkStrength(0.1)
                // .linkStrengthInterCluster(0.2)
                .gravityToFoci(0.01)
                .gravityOverall(0.01)
                .size([w, h])
                .enableGrouping(true)
                .groupBy("group2")
                .on("tick", tick);

        }
        
        rects = svg.append("g").attr("class","treemap");

        var drag = force.drag()
            .on("dragstart", dragstart);

        link = svg.selectAll(".link"),
            node = svg.selectAll(".node");
        
        force
            .nodes(graph.nodes)
            .links(graph.links)
            .start();
        path = svg.selectAll("path")
            .data(force.links());
        
        path.enter().insert("svg:path")
            .attr("class", "link")
            .style("fill", "none")
            // .style("stroke-width", "8")
            .style("stroke-width", function(d) { return d.size || 5; })
            .style("stroke", "#000");
        
        link = link.data(graph.links)
            .enter().append("line")
            .attr("class", "link")
            .style("stroke-width", function(d) { return d.size || 5; });

        link.style("visibility", "hidden");

        var color = d3.scale.category20();
        
        var n = graph.nodes.length;
        node = svg.selectAll(".node")
            .data(graph.nodes)
            .enter()
            .append("g")
            .attr("transform", function(d, i) {
                var x = width / n * i;
                return "translate(" + x + "," + x + ")";
            })
            .on("dblclick", dblclick)
            .call(drag);
        
        // var n = graph.nodes.length;
        // node.forEach(function(d, i) {
        //     d.x = d.y = width / n * i;
            
        // });
        // force.stop();

                    
        node.append("circle")
            // .attr("class", function(d) { return "node" + (d.size?"":" leaf"); })
            .attr("class", "node")
            // .attr("r", function (d) { return d.size ? assignSize(d.group) :  d.weight * 5 + 20; })
            .attr("r", function (d) { return d.size ? assignSize(d.group) :  20; })
            // .attr("r", 20)
            // .style("fill", function (d) { return assignRainbowColor(d.group); })
            .style("fill", function (d) { return color(d.weight); })
            .on('contextmenu', function(d){ 
                // http://bl.ocks.org/jakosz/ce1e63d5149f64ac7ee9
                d3.event.preventDefault();
                if (selected_single_cluster) {
                    console.log('already zoomed in')
                    prepare_data(false);
                } else { 
                    prepare_data(d.group2);
                }
                init();
            })
            
        node.append("text")
            .attr("dy", ".35em")
            .attr("text-anchor", 'middle')
            .attr("fill", "#000")
            .attr("cursor","pointer")
            .attr("font-size", function(d) { return d.size ? assignSize(d.group)-6 : 12; })
            .text(function(d){return d.size?d.group:d.name})
            .on('contextmenu', function(d){ 
                d3.event.preventDefault();
                if (selected_single_cluster) {
                    console.log('already zoomed in')
                    prepare_data(false);
                } else { 
                    prepare_data(d.group2);
                }
                init();
            });

        create_overlay();

        // setTimeout(function () {
        //     console.log('timer!');

        //     force.start();
        //     force.friction(0.8); 
        // }, 2000)
        
        // setTimeout(function () {
        //     console.log('timer!');
        //     force.start();
        //     force.friction(0.5); 
        // },4000)
        // if (!segment_view) {
        //     force.drawTreemap(svg);
        //     bindTreeMap(div);    
        // }

        function tick(e) {
            if (!segment_view) force.onTick(e);
            node.attr("transform", function (d) {
                d.x = Math.max(circle_size, Math.min(width - circle_size, d.x));
                d.y = Math.max(circle_size, Math.min(height - circle_size, d.y));
                return "translate(" + d.x + "," + d.y + ")";
            });

            path.attr("d", function(d) {
                var dx = d.target.x - d.source.x,
                    dy = d.target.y - d.source.y,
                    dr = Math.sqrt(dx * dx + dy * dy);
                return "M" + d.source.x + "," + d.source.y + "A" + dr + "," + dr + " 0 0,1 " + d.target.x + "," + d.target.y;
            });
            link.attr("x1", function(d) { return d.source.x; })
                .attr("y1", function(d) { return d.source.y; })
                .attr("x2", function(d) { return d.target.x; })
                .attr("y2", function (d) { return d.target.y; });
            
            
        }

    }
    // https://bl.ocks.org/mbostock/3750558
    function dblclick(d) {
    d3.select(this).classed("fixed", d.fixed = false);
    }

    function dragstart(d) {
    d3.select(this).classed("fixed", d.fixed = true);
    }

    function bindTreeMap() {
        // d3.selectAll(".cell").each(function (d, i) {
        //     console.log("i", i, "d", d);
        //     d.on("click", function(d){ console.log('clicked1')}) 
        // })
        div.selectAll(".cell").attr("fill", "#fff").on("click", function(d){ 
            if (selected_single_cluster!==false) {
                console.log('already zoomed in')
                prepare_data(false);
            } else {
                console.log('clicked tree, redraw only',d.id) 
                prepare_data(d.id);
            }
            init();
            })
        // d3.selectAll(".node").attr("fill", "#111").on("click", function (d) { console.log('clicked1 node') })
        
        $(containerSelector).find(".cell").appendTo(containerSelector+" .treemap");
        console.log('bind tree map!',containerSelector);
        // console.log($(".cell"));
        // $(".cell").click(function () {
        //     console.log('clicked!');
        // })
    }

    function network(data, prev, index, expand) {
        expand = expand || {};
        var gm = {},    // group map
            nm = {},    // node map
            lm = {},    // link map
            gn = {},    // previous group nodes
            gc = {},    // previous group centroids
            nodes = [], // output nodes
            links = []; // output links
    
        // process previous nodes for reuse or centroid calculation
        if (prev) {
            prev.nodes.forEach(function(n) {
            var i = index(n), o;
            if (n.size > 0) {
                gn[i] = n;
                n.size = 0;
            } else {
                o = gc[i] || (gc[i] = {x:0,y:0,count:0});
                o.x += n.x;
                o.y += n.y;
                o.count += 1;
            }
            });
        }
    
        // determine nodes
        for (var k=0; k<data.nodes.length; ++k) {
            var n = data.nodes[k],
                i = index(n),
                l = gm[i] || (gm[i]=gn[i]) || (gm[i]={group:i, size:0, nodes:[]});
    
            if (expand[i]) {
            // the node should be directly visible
            nm[n.name] = nodes.length;
            nodes.push(n);
            if (gn[i]) {
                // place new nodes at cluster location (plus jitter)
                n.x = gn[i].x + Math.random();
                n.y = gn[i].y + Math.random();
            }
            } else {
            // the node is part of a collapsed cluster
            if (l.size == 0) {
                // if new cluster, add to set and position at centroid of leaf nodes
                nm[i] = nodes.length;
                nodes.push(l);
                if (gc[i]) {
                l.x = gc[i].x / gc[i].count;
                l.y = gc[i].y / gc[i].count;
                }
            }
            l.nodes.push(n);
            }
        // always count group size as we also use it to tweak the force graph strengths/distances
            l.size += 1;
        n.group_data = l;
        }
    
        for (i in gm) { gm[i].link_count = 0; }
    
        // determine links
        for (k=0; k<data.links.length; ++k) {
            var e = data.links[k],
                u = index(e.source),
                v = index(e.target);
        if (u != v) {
            gm[u].link_count++;
            gm[v].link_count++;
        }
            u = expand[u] ? nm[e.source.name] : nm[u];
            v = expand[v] ? nm[e.target.name] : nm[v];
            var i = (u<v ? u+"|"+v : v+"|"+u),
                l = lm[i] || (lm[i] = {source:u, target:v, size:0});
            l.size += 1;
        }
        for (i in lm) { links.push(lm[i]); }
    
        return {nodes: nodes, links: links};
    }
        
    function getGroup(n) { return n.group; }



    function assignSize(segment) {

        switch (segment) {
        case "TM1":
            size = 20;
            break;
        case "TM2":
            size = 20;
            break;
        case "TM3":
            size = 20;
            break;
        case "TM4":
            size = 20;
            break;
        case "TM5":
            size = 20;
            break;
        case "TM6":
            size = 20;
            break;
        case "TM7":
            size = 20;
            break;
        case "H8":
            size = 18;
            break;
        default:
            size = 15;
        }
        return size;
    }
    
    function assignRainbowColor(segment) {
        var segmentRainbowColors2 = {
            "1": "#736DA7",
            "2": "#5EB7B7",
            "3": "#CE9AC6",
            "4": "#DD7D7E",
            "5": "#E6AF7C",
            "6": "#DEDB75",
            "7": "#80B96F",
            "8": "#EEE",
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

    function create_overlay() {
        var newDiv = document.createElement("div");

        $(containerSelector).find(".flareplot-legend").remove();

        newDiv.setAttribute("class", "flareplot-legend");

        var content = '<div class="controls">'
        //                                  +'<h4>Controls</h4>';


        content += '<p>Line colors: <select id="flareplot_color">' +
            '<option value="none">None (gray)</option>' +
            '<option value="rainbow">GPCR rainbow</option>' +
            '<option value="segment">GPCR segment</option>';

        var mode = get_current_mode();
        // if single structure - use interaction coloring
        if (mode == "single-crystal") {
            content += '<option value="interactions" selected>Interaction Type</option>';
            // if single group of structures - use frequency coloring (gradient)
        } else if (mode == "single-crystal-group") {
            content += '<option value="frequency" selected>Interaction Frequency/Count</option>';
            // if group(s) of structures - use frequency coloring (gradient)
        } else {
            content += '<option value="frequency" selected>Frequency difference Gr1 - Gr2</option>';
            content += '<option value="frequency_1">Frequency group 1</option>';
            content += '<option value="frequency_2">Frequency group 2</option>';
        }
        content += '</select></p>';
        content = '<span class="options">' +
        '<input id="checkCurved" type="checkbox" checked>' +
        '<span class="checkboxtext"> Curved links' +
        '</span>' +
        '</input>' +
        '<br><button class="btn btn-primary btn-xs" id="resetfixed">Release fixed</button>' +
        '<br><button class="btn btn-primary btn-xs" id="freeze">Freeze all</button>' +
        '</span>';
        // content = '';
        newDiv.innerHTML = content;

        $(containerSelector).append(newDiv);

        d3.select(containerSelector).select("#checkCurved").on("change", function () {
            curved_links = d3.select(containerSelector).select("#checkCurved").property("checked");
            console.log('they changed!', curved_links, containerSelector)
            if (curved_links) {
                link.style("visibility", "hidden");
                path.style("visibility", "");

            } else {
                path.style("visibility", "hidden");
                link.style("visibility", "");
            }

            // init();
        });
                
        d3.select(containerSelector).select("#resetfixed").on("click", function () {
            node.classed("fixed", function (d) {
                d.fixed = false;
            });
            force.start();
        });

        d3.select(containerSelector).select("#freeze").on("click", function () {
            node.classed('fixed', function(d, i) {
                d.fixed = true;
            });
        });

        
        
    }
    
}

Array.prototype.indexOfForArrays = function(search)
{
  var searchJson = JSON.stringify(search); // "[3,566,23,79]"
  var arrJson = this.map(JSON.stringify); // ["[2,6,89,45]", "[3,566,23,79]", "[434,677,9,23]"]
    console.log("hi!",arrJson, searchJson,arrJson.indexOf(searchJson));
  return arrJson.indexOf(searchJson);
};
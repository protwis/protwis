function createNetworkPlot3D(raw_data,original_width, inputGraph, containerSelector, segment_view = true) {
    // https://github.com/d3/d3-3.x-api-reference/blob/master/Force-Layout.md
     //https://archive.nytimes.com/www.nytimes.com/interactive/2013/02/20/movies/among-the-oscar-contenders-a-host-of-connections.html
    
        // Other ideas 3D: https://bl.ocks.org/vasturiano/f59675656258d3f490e9faa40828c0e7
        
        
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
        var stickyDrag = true;
        var highlightNode = true;
        function prepare_data(single_cluster = false) {
    
            console.log('PREPARE DATA',containerSelector,'cluster_groups',cluster_groups.length, plot_specified_filtered.length )
    
            selected_single_cluster = single_cluster;
    
            new_data = { "nodes": [], "links": [] };
            var track_gns = [];
    
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
                    new_data["nodes"].push({ "name": gns[0], "group": seg1, "group2":cg_index1, "links":0 })
                    track_gns.push(gns[0]);
                }
                if (!track_gns.includes(gns[1])) {
                    new_data["nodes"].push({ "name": gns[1], "group": seg2, "group2":cg_index2, "links":0 })
                    track_gns.push(gns[1]);
                }
    
                source = new_data["nodes"].filter(obj => { return obj.name === gns[0] })[0];
                source.links += 1;
                target = new_data["nodes"].filter(obj => { return obj.name === gns[1] })[0];
                target.links += 1;
                new_data["links"].push({ "source": source, "target": target, "value": 1 })
            })
    
            if (!segment_view) {
                // resize by number of nodes if not segments
                console.log(new_data["nodes"].length, "nodes");
                width = original_width*Math.sqrt(new_data["nodes"].length/2);
            } else {
                width = original_width*0.8;
            }
            w = width;
            h = w;
            height = w;
    
        }
        prepare_data();
    
    
        var expand = {}, net, force;
    
        d3v4.select(containerSelector).style("position","relative");
    
        div = d3v4.select(containerSelector).insert("div")
            .attr("class", "network")
            .style("width", "100%")
            .style("-webkit-backface-visibility", "hidden");
    
        init();
        
        var path, link, node, svg;
        function init() {
            console.log('init',div, containerSelector);
            div = d3v4.select(containerSelector);
            div.classed("3dnetwork", true); 
            // div.select("svg").remove();
            // // if (force) force.stop();
            // svg = div.append("svg:svg")
            //     .attr("id","3d-graph")
            //     .attr("viewBox", "0 0 " + w + " " + h)
            //     .attr("width", "100%")
            //     .attr("style", "height: 500px");
            
            // // slightly hide the initial load
            // svg.attr("opacity", 1e-6)
            //     .transition()
            //     .duration(3000)
            //     .attr("opacity", 1);
    
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
                graph['nodes'].forEach(function (n) {
                    n.name = n.group;
                });

                // normalize
                graph['links'].forEach(function (n) {
                    n.size = Math.round(max_link_size * n.size / max_link);
                    n.source = graph["nodes"][n.source];
                    n.target = graph["nodes"][n.target];
                });
            }
            console.log("3d", graph);
            console.log(div.attr('id'));

            width = $(containerSelector).width();

            const elem = document.getElementById(div.attr('id'));
            var myGraph = ForceGraph3D()
                (document.getElementById(div.attr('id')))
                .width(width)
                .height("500")
                // .backgroundColor("#fff")
                .linkColor("#fff")
                .linkOpacity(1)
                .graphData(graph)
                .nodeLabel('name')
                .nodeAutoColorBy('group')
                .onNodeHover(node => elem.style.cursor = node ? 'pointer' : null)
                .onNodeClick(node => {
                    console.log('node click!');
                    // Aim at node from outside it
                    const distance = 40;
                    const distRatio = 1 + distance / Math.hypot(node.x, node.y, node.z);
                    myGraph.cameraPosition(
                        { x: node.x * distRatio, y: node.y * distRatio, z: node.z * distRatio }, // new position
                        node, // lookAt ({ x, y, z })
                        3000  // ms transition duration
                    );
                });
            
            $(window).resize(function () {
                console.log('3d resize', is_fullscreen);
                setTimeout(
                    function() {
                    if (is_fullscreen) {
                        myGraph.refresh();
                        width = $(containerSelector).width();
                        width = window.innerWidth;
                        height = window.innerHeight;
                        myGraph.width(width);
                        myGraph.height(height);
                    } else {
                        console.log('resize 3d to div');
                        width = $(containerSelector).width();
                        myGraph.width(width);
                        myGraph.height(500);
    
                    }
                },500
                )
                            // Resize 3D network
            // console.log('"check if resize 3d back down');
            // if ($(".3dnetwork").length && !is_fullscreen) {
            //     var myGraph = ForceGraph3D()
            //     (document.getElementById($(".3dnetwork").attr('id')))
            //     if (myGraph) {
            //         width = $(".3dnetwork").width();
            //         myGraph.height("500");
            //         myGraph.width(width);
            //     }
            // }

            });
               
    
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
    }

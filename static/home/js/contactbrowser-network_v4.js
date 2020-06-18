function createNetworkPlot(raw_data,original_width, inputGraph, containerSelector, segment_view = true) {
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
    var highlightNode = segment_view ? false : true;
    var max_freq = 0;
    var colorLinks = true;

    var filter_sets = 'all';

    var nice_index_names = {
        'solid': 'Color',
        'conservation' : 'Conservation of set(s) consensus AA in class',
        'outer_angle': 'Rotamer',
        // 'rotation_angle': 'Rotation angle',
        'phi': 'Phi',
        'psi': 'Psi',
        'tau': 'Tau',
        'distance': 'Distance',
        'distance_abs': 'Distance (abs)',
        'core_distance': 'Distance to 7TM axis',
        'core_distance_abs': 'Distance to 7TM axis (abs)',  
        'a_angle' : 'Angle to helix&7TM axes',
        'theta': 'Theta',
        'tau_angle': 'Tau dihedral',
        // 'sasa': 'SASA',
        // 'sasa_abs': 'SASA (abs)',
        // 'rsa': 'RSA',
        // 'rsa_abs': 'RSA (abs)',
        // 'hse': 'HSE',
        // 'hse_abs': 'HSE (abs)',
        // 'network' : 'Network group no.',
        // 'set_presence' : 'Set specific presense',
        // 'ligand' : 'Ligand interactions freq',
        // 'complex' : 'G protein interactions',
        // 'ligandcomplex' : 'Ligand and G protein interactions',
        'mutations' : 'Mutations with >5 fold effect',
        'cons_prop' : 'Consensus AA property ',
    }



    function prepare_data(single_cluster = false) {

        console.log('PREPARE DATA',containerSelector,'cluster_groups',cluster_groups.length, plot_specified_filtered.length )

        selected_single_cluster = single_cluster;

        new_data = { "nodes": [], "links": [] };
        var track_gns = [];
        // cluster_groups = [];
        // console.log('single_cluster',single_cluster);
        // // // Populate matrix for interactions between segments
        // track_gns = []
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
        if (raw_data['proteins2']) {
            var pdbs_1 = raw_data['pdbs1'].length
            var pdbs_2 = raw_data['pdbs2'].length
            var pfs_1 = raw_data['pfs1'].length
            var pfs_2 = raw_data['pfs2'].length
            var normalized = raw_data['normalized'];
        } else if (raw_data['proteins'].length > 1) {
            var pdbs_counts = raw_data['pdbs'].length
            var normalized = raw_data['normalized'];
            var pfs = raw_data['pfs'].length
        }
        $.each(raw_data['interactions'], function (i, v) {
            if (!plot_specified_filtered.includes(i)) return;
            gns = separatePair(i);
            var freq = 0;
            var sfreq1 = 0;
            var sfreq2 = 0;
            if (raw_data['proteins2']) {
                if (normalized) {
                    sfreq1 = Math.round(100 * v['pfs1'].length / pfs_1);
                    sfreq2 = Math.round(100 * v['pfs2'].length / pfs_2);
                } else {
                    sfreq1 = Math.round(100 * v['pdbs1'].length / pdbs_1);
                    sfreq2 = Math.round(100 * v['pdbs2'].length / pdbs_2);
                }
                freq = sfreq1 - sfreq2;
            } else if (raw_data['proteins'].length > 1) {
                if (normalized) {
                    freq = Math.round(100 * v['pfs'].length / pfs);
                } else {
                    freq = Math.round(100 * v['pdbs'].length / pdbs_counts);
                }
                sfreq1 = Math.round(freq / 2);
                sfreq2 = freq - sfreq1;
            }
            seg1 = raw_data['segment_map'][gns[0]];
            seg2 = raw_data['segment_map'][gns[1]];

            cg = cluster_groups.filter(l => l.includes(gns[0]));
            cg_index1 = cluster_groups.indexOf(cg[0]);
            cg = cluster_groups.filter(l => l.includes(gns[1]));
            cg_index2 = cluster_groups.indexOf(cg[0]);

            set_type = filtered_cluster_groups_set[cg_index1];
            if (filter_sets=='all' || filter_sets==set_type) {
                if (single_cluster!==false && (cg_index1!=single_cluster || cg_index2!=single_cluster)) return;

                if (!track_gns.includes(gns[0])) {
                    new_data["nodes"].push({ "name": gns[0], "group": seg1, "group2":cg_index1, "group3":set_type, "links":0 })
                    track_gns.push(gns[0]);
                }
                if (!track_gns.includes(gns[1])) {
                    new_data["nodes"].push({ "name": gns[1], "group": seg2, "group2":cg_index2, "group3":set_type, "links":0 })
                    track_gns.push(gns[1]);
                }

                source = new_data["nodes"].filter(obj => { return obj.name === gns[0] })[0];
                source.links += 1;
                target = new_data["nodes"].filter(obj => { return obj.name === gns[1] })[0];
                target.links += 1;
                new_data["links"].push({ "source": source, "target": target, "value": 1, "freq": freq, "sfreq1": sfreq1, "sfreq2": sfreq2 })
                if (Math.abs(freq) > max_freq) max_freq = Math.abs(freq);
            }
        })

        // console.log(new_data["nodes"]);
        // console.log(cluster_groups);
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
            width = original_width*Math.sqrt(new_data["nodes"].length/2);
        } else {
            width = original_width*1;
        }
        w = width;
        h = w;
        height = w;

    }
    prepare_data();

    // See snakeplot for this code. Placed there since it reused alot of logic from the snakeplot coloring.
    var colors_gn = prepare_residue_colors(raw_data);




    var expand = {},net, force;

    d3v4.select(containerSelector).style("position","relative");

    div = d3v4.select(containerSelector).insert("div")
        .attr("class", "network")
        .style("width", "100%")
        .style("-webkit-backface-visibility", "hidden");

    init();
    
    var path, link, node, svg,labelText;
    function init() {
        div = d3v4.select(containerSelector).select("div");
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

            // Ensure all segments are there..
            $.each(['TM1', 'TM2', 'TM3', 'TM4', 'TM5', 'TM6', 'TM7', 'H8','ICL1','ICL2','ECL1','ECL2'], function (i, s) {
                if (!(graph['nodes'].find(({ group }) => group === s))) {
                    graph['nodes'].push({'group':s, nodes: [], link_count: 0, size:1})
                }
            })

            max_size = Math.max.apply(Math, graph['nodes'].map(function(v) {
                return v.size;
            }));

            max_link = Math.max.apply(Math, graph['links'].map(function(v) {
                return v.size;
            }));

            // normalize
            total_links = 0;
            graph['links'].forEach(function (n) {
                n.links = n.size;
                total_links += n.links;
                n.size = Math.max(1,Math.round(max_link_size*n.size/max_link));
            });
        }
    
        // var link = svg.append("g")
        //     .attr("class", "links")
        //     .selectAll("line")
        //     .data(graph.links)
        //     .enter()
        //     .append("line")
        //     .attr("class", "link")
        //     .style("stroke-width", function(d) { return d.size || 5; })
        //     .attr("stroke", "black")
    
        var link = svg.append("g")
            .attr("class", "links")
            .selectAll("line")
            .data(graph.links)
            .enter()
            .append("path")
            .attr("class", "link")
            .style("fill", "none")
            .style("stroke-width", function(d) { return d.size || 5; })
            .style("stroke", "#000")
            .attr("id", function (d, i) { return containerSelector + "linkId_" + i; });
        
        if (segment_view) {
            var link1 = svg.append("g")
                .attr("class", "links")
                .selectAll("line")
                .data(graph.links)
                .enter()
                .append("path")
                .attr("class", "link")
                .style("fill", "none")
                .style("stroke-width", function (d) { return 5 })
                .style("stroke", "#F00");
        }
            // .style("fill", "none")
            // // .style("stroke-width", "8")
            // .style("stroke-width", function(d) { return d.size || 5; })
            // .style("stroke", "#000");
        
        var labelParent = svg.selectAll(".labelText")
            .data(graph.links)
            .enter().append("text")
            .attr("class", "labelText")
            .attr("dx", 0)
            .attr("dy", function (d, i) { return d.size ? -d.size / 2 : -5 / 2; })
            .style("fill", "black")
            .style("opacity", 0.5)
            .text(function(d,i) { return  d.links;})
            .attr("text-anchor", "Middle")
            .attr("id", function (d, i) { return "labelText_" + i; });
        
        var labelText = labelParent.append("textPath")
            .attr("xlink:href", function (d, i) { return "#" + containerSelector + "linkId_" + i; })
            .attr("startOffset", "50%").attr("text-anchor", "Middle");
            // .text(function(d,i) { return  d.links;});
        
        var n = graph.nodes.length;
        var node = svg.selectAll(".node")
            .data(graph.nodes)
            .enter()
            .append("g")
            .attr("transform", function(d, i) {
                var x = width / n * i;
                return "translate(" + x + "," + x + ")";
            })
            .call(d3v4.drag()
                .on("start", dragstarted)
                .on("drag", dragged)
                .on("end", dragended))
            .on("dblclick", dblclick)
            .on('contextmenu', function (d) { 
                    if (!segment_view) {
                        d3v4.event.preventDefault();
                        if (selected_single_cluster) {
                            console.log('already zoomed in')
                            prepare_data(false);
                        } else { 
                            prepare_data(d.group2);
                        }
                        init();
                    }
            })
            .on("mouseover", (d) => {
                if (!dragged && highlightNode) {
                    highlight_network(d)
                }
            })
            .on("mouseout", (d) => {
                if (!dragged && highlightNode) {
                    link.style("stroke-opacity", 0.5)
                    node.style("opacity", 1)
                }
            });
        
        
        node.append("circle")
            .attr("class", "node")
            .attr("r", function (d) { return d.size ? assignSize(d.group) : 20; })
            // .attr("r", function (d) { return d.links**2+20; })
            .style("opacity", 1)
            .style("fill", function (d) { return assignRainbowColor(d.group); });
        
        if (segment_view) {
            node.select("circle").attr("visibility", function (d) { return d.group.startsWith("TM") ? "visible" : "hidden"; })
            node.select("circle").style("stroke", "#000");

            node.append("rect")
                .attr("class", "node")
                .attr("visibility", function (d) { return d.group.startsWith("TM") ? "hidden" : "visible" ; }) // H8 and loops are <20..
                .attr("x", function (d) { return d.group=='H8' ? -15 : -13; })
                .attr("y", function (d) { return d.group == 'H8' ? -13 : -6; })
                .attr("rx", function (d) { return d.group=='H8' ? 0 : 6; })
                .attr("ry", function (d) { return d.group=='H8' ? 0 : 6; })
                .attr("width", function (d) { return d.group=='H8' ? 30 : 25; })
                .attr("height", function (d) { return d.group=='H8' ? 26 : 12; })
                .attr("stroke", "#000")
                .attr("stroke-width",function (d) { return d.group=='H8' ? 2 : 1; })
                .style("opacity", 1)
                .style("fill", function (d) { return assignRainbowColor(d.group); });
        }
        
        node.append("text")
            .attr("dy", ".35em")
            .attr("text-anchor", 'middle')
            .attr("fill", "#000")
            .attr("cursor", "pointer")
            .attr("font-size", function (d) { return d.size ? assignSize(d.group) - 6 : 14; })
            .text(function (d) { return d.size ? d.group : d.name });
        
        if (segment_view) {
            // Fix inital layout for segments
            node.attr("transform", function (d) {
                if (assignPosition(d.group)[0] != 0 ) {
                    [d.x, d.y] = assignPosition(d.group);
                }
                d.fixed = true;
                d.fx = d.x;
                d.fy = d.y;
                return "translate(" + d.x + "," + d.y + ")";
            });
        }
        
        var ticked = function () {
            colorLinks = false;
            if ($(containerSelector +" #colorLinks").length) {
                colorLinks = d3.select(containerSelector).select("#colorLinks").property("checked");
            }
            
            node.attr("transform", function (d) {
                size = d.size ? assignSize(d.group) : 20
                size += 2;
                d.x = Math.max(size, Math.min(width - size, d.x));
                d.y = Math.max(size, Math.min(height - size, d.y));
                return "translate(" + d.x + "," + d.y + ")";
            });
            

            link
                .attr("x1", function(d) { return d.source.x; })
                .attr("y1", function(d) { return d.source.y; })
                .attr("x2", function(d) { return d.target.x; })
                .attr("y2", function (d) { return d.target.y; });
            
            labelParent
                .attr("x", function(d) { return (d.source.x+d.target.x)/2; })
                .attr("y", function(d) { return (d.source.y+d.target.y)/2; })
            
            if (segment_view) {
                link1
                    .attr("x1", function (d) { return d.source.x; })
                    .attr("y1", function (d) { return d.source.y; })
                    .attr("x2", function (d) { return d.target.x; })
                    .attr("y2", function (d) { return d.target.y; });

                link1.attr("d", function(d, i) {
                    var dx = d.target.x - d.source.x,
                        dy = d.target.y - d.source.y,
                        dr = Math.sqrt(dx * dx + dy * dy);
                    
                        offset_x = 0;
                        offset_y = 0;
                        
                        if (segment_view && colorLinks) {
                            offset = 10*d.sfreq2 / max_freq;
                            var vector = dx / dy;
                            //var atan = Math.atan(vector);
                            var rad = Math.atan2(dy, dx);
                            var degrees = -rad * (180 / Math.PI);
                            degrees = degrees ? -degrees : 0;
                            offset_x = offset * Math.cos(Math.PI/2 + rad);
                            offset_y = offset * Math.sin(Math.PI/2 + rad);
                            //console.log(d.source.group,d.target.group, vector,degrees,rad);
                        }
                        if (dx > 0) {
                            return "M" + (d.source.x+offset_x) + " " + (d.source.y+offset_y) + " L " + (d.target.x+offset_x) + " " + (d.target.y+offset_y);
                        } else {
                            return "M" + (d.target.x-offset_x) + " " + (d.target.y-offset_y) + " L " + (d.source.x-offset_x) + " " + (d.source.y-offset_y);
                        }
                });
                link1.attr("visibility", function(d,i) { return colorLinks ? "visible" : "hidden"})
            }

            link.attr("d", function(d, i) {
                var dx = d.target.x - d.source.x,
                    dy = d.target.y - d.source.y,
                    dr = Math.sqrt(dx * dx + dy * dy);
                
                offset_x = 0;
                offset_y = 0;
                
                if (segment_view && colorLinks) {
                    offset = 10*d.sfreq1 / max_freq;
                    var vector = dx / dy;
                    //var atan = Math.atan(vector);
                    var rad = Math.atan2(dy, dx);
                    var degrees = -rad * (180 / Math.PI);
                    degrees = degrees ? -degrees : 0;
                    offset_x = offset * Math.cos(Math.PI/2 + rad);
                    offset_y = offset * Math.sin(Math.PI/2 + rad);
                    //console.log(d.source.group,d.target.group, vector,degrees,rad);
                }
                if (dx > 0) {
                    return "M" + (d.source.x-offset_x) + " " + (d.source.y-offset_y) + " L " + (d.target.x-offset_x) + " " + (d.target.y-offset_y);
                } else {
                    return "M" + (d.target.x+offset_x) + " " + (d.target.y+offset_y) + " L " + (d.source.x+offset_x) + " " + (d.source.y+offset_y);
                }
            });
        }  
        // tooltip https://observablehq.com/@skofgar/force-directed-graph-integrated-html

        var simulation = d3v4.forceSimulation()
            // .force("collide",d3v4.forceCollide( function(d){return 30 }).iterations(1) )
            // .force("collide",d3v4.forceCollide( function (d) { return d.links**2+10; }).strength(0.5).iterations(5))
            // .alphaDecay(0)
            .force("charge", d3v4.forceManyBody().strength(-700))
            .force("x", d3v4.forceX(w/2).strength(0.0))
            .force("y", d3v4.forceY(h/2).strength(0.0));
        
        link_distance = 70;
        link_strength = 0.7;
        gravity = 0;
        charge = -700;
        var distanceMax = 0;
        collide = 0;
        var useGroupingForce = false;
        // simulation.alphaDecay(0.001);
        if (selected_single_cluster === false && cluster_groups.length > 1 && !segment_view) {

            max_cluster_size = Math.max.apply(Math, cluster_groups.map(function(v) {
                return v.length;
            }));
            
            distanceMax = max_cluster_size*20;
            simulation.force("charge", d3v4.forceManyBody()
                .strength(charge)
                .distanceMax(distanceMax)
            )
            gravity = 0.1;
            useGroupingForce = true;
            // Instantiate the forceInABox force
            var groupingForce = forceInABox()
                .strength(gravity) // Strength to foci
                .template("treemap") // Either treemap or force
                .groupBy("group2") // Node attribute to group
                // .links(graph.links) // The graph links. Must be called after setting the grouping attribute
                .size([width, height]) // Size of the chart
            simulation
                .force("group", groupingForce)
            simulation.force("collide", d3v4.forceCollide(function (d) { return d.links ** 2 + 20; }).strength(1).iterations(1))
        } else {
            if (segment_view) {
                console.log('segment! collide');
                collide = 50;
                simulation.force("collide", d3v4.forceCollide(function (d) { return 50; }).strength(1).iterations(1))
                link_distance = 90;
                link_strength = 0.7;
                // link_strength = function (l) { return l.size / max_link_size };
                gravity = 0.04;
                simulation
                .force("x", d3v4.forceX(w/2).strength(gravity))
                .force("y", d3v4.forceY(h/2).strength(gravity));

                // simulation.force("charge", d3v4.forceManyBody().strength(function (d) { return -d.size*50 }))
            } else {
                // charge = -1200;
                gravity = 0.05;
                simulation.force("charge", d3v4.forceManyBody().strength(charge))
                    .force("x", d3v4.forceX(w/2).strength(gravity))
                    .force("y", d3v4.forceY(h / 2).strength(gravity));
                
                // simulation.force("charge", d3v4.forceManyBody().strength(function (d) { return -(d.links**3) }))
                simulation.force("collide", d3v4.forceCollide(function (d) { return d.links ** 2 + 20; }).strength(1).iterations(1))
            }
        }
        // console.log('distance_max', distance_max, cluster_groups.length, selected_single_cluster);
        
      
        // simulation.force("collide", d3v4.forceCollide(function (d) { return 0; }).strength(1).iterations(1))
        

        
        // simulation
        //     .nodes(graph.nodes)
        //     .on("tick", ticked);
    
        // simulation.force("link")
        //     .links(graph.links);  

        simulation
            .nodes(graph.nodes)
            .force("link", d3v4.forceLink(graph.links)
                .distance(link_distance)
                .strength(link_strength).iterations(1)
                // .strength(0.8).iterations(1)
            //   .strength(groupingForce.getLinkStrength) // default link force will try to join nodes in the same group stronger than if they are in different groups
            ).on("tick", ticked);
        
        // if (segment_view) simulation.stop();
          
        // https://bl.ocks.org/mbostock/3750558
        function dblclick(d) {
            if (!d3v4.event.active) simulation.alphaTarget(0);
            // d3v4.select(this).classed("fixed", d.fixed = false);
            d3v4.select(this).select("circle").classed("fixed", false);
            d.fx = null;
            d.fy = null;
        }

        function highlight_network(d) {
            ds = [d.name];
            already_done = [];
            link.transition();
            node.transition();
            link.style("stroke-opacity", 0.05)
            node.style("opacity", 0.05);
            delay_step = 0;
            steps = 4;
            for (let step = 0; step < steps; step++) {
                ds2 = []
                link.filter(l => !already_done.includes(l) && ( ds.includes(l.source.name) || ds.includes(l.target.name))).style('stroke-opacity', function (l) {
                        ds2.push(l.source.name, l.target.name);
                        already_done.push(l)
                        return 1-0.8*(step/steps)
                })
                ds = ds2;
                node.filter(d => !already_done.includes(d) && ds.includes(d.name)).style('opacity', function (d) {
                    already_done.push(d)    
                    return 1 - 0.8*(step/steps);
                })
                // console.log(step, ds);
            }
        }
        
        dragged = false;
        function dragstarted(d) {
            dragged = true;
            d3v4.select(this).select("circle").classed("fixed", true);
            if (!d3v4.event.active) simulation.alphaTarget(0.3).restart();
            d.fx = d.x;
            d.fy = d.y;
        }
        function dragged(d) {
            dragged = true;
            d.fx = d3v4.event.x;
            d.fy = d3v4.event.y;
        }
        
        function dragended(d) {
            dragged = false;
            if (!d3v4.event.active) simulation.alphaTarget(0);
            if (!stickyDrag) {
                d.fx = null;
                d.fy = null;
                d3v4.select(this).select("circle").classed("fixed", false);
            }
        } 

        create_overlay();

        d3v4.select(containerSelector).select("#link_strength_change")
            .on("input", link_strength_change);
        
        function link_strength_change() {
            console.log('link strength', this.value);
            simulation.force("link").strength(+this.value);
            if (!d3v4.event.active) simulation.alpha(1).restart();
        }
        
        d3v4.select(containerSelector).select("#link_distance_change")
            .on("input", link_distance_change);
        
        function link_distance_change() {
            console.log('link distance', this.value);
            simulation.force("link").distance(+this.value);
            if (!d3v4.event.active) simulation.alpha(1).restart();
        }
        
        d3v4.select(containerSelector).select("#gravity_change")
            .on("input", gravity_change);
        
        function gravity_change() {
            gravity = this.value;
            console.log('gravity', gravity);
            if (useGroupingForce) {
                var groupingForce = forceInABox()
                    .strength(+this.value) // Strength to foci
                    .template("treemap") // Either treemap or force
                    .groupBy("group2") // Node attribute to group
                    // .links(graph.links) // The graph links. Must be called after setting the grouping attribute
                    .size([width, height]) // Size of the chart
            
                simulation.force("group", groupingForce);
            } else {
                simulation.force("x", d3v4.forceX(w/2).strength(gravity))
                .force("y", d3v4.forceY(h/2).strength(gravity));    
            }
            if (!d3v4.event.active) simulation.alpha(1).restart();
        }
        
        d3v4.select(containerSelector).select("#charge_change")
            .on("input", charge_change);
        
        function charge_change() {
            console.log('charge', this.value);
            if (distanceMax) {
                simulation.force("charge", d3v4.forceManyBody()
                    .strength(+(-1 * this.value))
                    .distanceMax(distanceMax)
                    )
            } else {
                simulation.force("charge", d3v4.forceManyBody()
                        .strength(+(-1*this.value))
                    )
            }
            if (!d3v4.event.active) simulation.alpha(1).restart();
        }

        d3v4.select(containerSelector).select("#collide_change")
            .on("input", collide_change);
        
        function collide_change() {
            console.log('collide', this.value);
            simulation.force("collide", d3v4.forceCollide(this.value).strength(1).iterations(1))
            if (!d3v4.event.active) simulation.alpha(1).restart();
        }

        d3v4.select(containerSelector).select("#size_change")
            .on("input", size_change);
        
        function size_change() {
            width = height = w = h = this.value;
            svg.attr("viewBox", "0 0 " + w + " " + h);
            // simulation.force("center", d3v4.forceCenter(w / 2, h / 2));
            if (useGroupingForce) {
                var groupingForce = forceInABox()
                .strength(gravity) // Strength to foci
                .template("treemap") // Either treemap or force
                .groupBy("group2") // Node attribute to group
                // .links(graph.links) // The graph links. Must be called after setting the grouping attribute
                .size([width, height]) // Size of the chart
            
                simulation.force("group", groupingForce);
            }

            simulation.force("x", d3v4.forceX(w/2))
                      .force("y", d3v4.forceY(h/2));
            simulation.alpha(1).restart();
        }
        
        d3.select(containerSelector).select("#stickyDrag").on("change", function () {
            stickyDrag = d3.select(containerSelector).select("#stickyDrag").property("checked");
        });
        d3.select(containerSelector).select("#highlightNode").on("change", function () {
            highlightNode = d3.select(containerSelector).select("#highlightNode").property("checked");
        });

        d3.select(containerSelector).select("#linkLabel").on("change", function () {
            linkLabel = d3.select(containerSelector).select("#linkLabel").property("checked");
            labelText.attr("visibility", linkLabel ? "visible" : "hidden");
            
        });


        d3.select(containerSelector).select("#addAA").on("change", function () {
            addAA = d3.select(containerSelector).select("#addAA").property("checked");
            node.select("circle").attr("r", function (d) { return addAA ? 20 : 20; })
            node.select("text").html(
                function (d) { return addAA ? "<tspan x=0 dy=-3>"+raw_data["tab4"][d.name]["all_seq_cons"][0]+"</tspan><tspan x=0 dy=12>"+d.name+"</tspan>" : d.name }
            );
        });

        d3.select(containerSelector).select("#change_to_freq").on("change", function () {
            changeFreq = d3.select(containerSelector).select("#change_to_freq").property("checked");
            // labelText.text(function (d) { return changeFreq ? (100*d.links / total_links).toFixed(0)+"%" :  d.links  });
            labelParent.text(function (d) { return changeFreq ? (100*d.links / total_links).toFixed(0)+"%" :  d.links  });
        });

        $(containerSelector+ " .residue_rotation").on("change", function () {
            change_rotation();
        });
        function change_rotation() {
            fill_color = $(containerSelector + " #residue_rotation").val();
            console.log('change rotation to!', fill_color);

            if (fill_color == 'none') {
                node.select("text").attr("transform", "");
            } else {
                // loop over nodes, find name
                node.select("text").attr("transform", function (d) {
                    var name = d.name;
                    var scale = colors_gn[fill_color][name][1];
                    var rotation_value = Math.round(scale * 180 - 90);
                    if (fill_color == 'outer_angle' || fill_color == 'rotation_angle') {
                        // For these 'actual' rotation values, use the 'true value' as the rotation.
                        var rotation_value = Math.round(colors_gn[fill_color][name][0]);
                    }

                    return "rotate("+rotation_value+")";
                })
            }
        }

        d3.select(containerSelector).select(".residue_fill").on("change", function () {
            change_fill();
        });
        function change_fill() {
            fill_color = $(containerSelector + " #node_color").val();

            color_id1  = $(containerSelector+" #fill_color1").spectrum("get").toHexString();
            if ($(containerSelector + " #fill_color2").spectrum("get") && $(containerSelector + " #fill_color2").length) {
                color_id2 = $(containerSelector + " #fill_color2").spectrum("get").toHexString();
            } else {
                color_id2 = false;
            }
            if ($(containerSelector + " #fill_color3").spectrum("get") && $(containerSelector + " #fill_color3").length) {
                color_id3 = $(containerSelector + " #fill_color3").spectrum("get").toHexString();
            } else {
                color_id3 = false;
            }
            console.log('change fill color to!', fill_color, color_id1, color_id2, color_id3);
            
            // if no color 3 make only linear between two colors.
            if (color_id3) {
                var color_range = d3v4.scaleLinear()
                    .domain([0, 0.5, 1])
                    .range([color_id1, color_id2, color_id3]);
            } else {
                var color_range = d3v4.scaleLinear()
                    .domain([0, 1])
                    .range([color_id1, color_id2]);
            }
            if (fill_color == 'segment') {
                node.select("circle").style("fill", function (d) { return assignRainbowColor(d.group); });
                node.select("rect").style("fill", function (d) { return assignRainbowColor(d.group); });
            } else if (fill_color == 'solid') {
                node.select("rect").style("fill", color_id1);
                node.select("circle").style("fill", color_id1);
            } else {
                // loop over nodes, find name
                node.select("circle").style("fill", function (d) {
                    name = d.name;
                    if (colors_gn[fill_color][name]) {
                        return color_range(colors_gn[fill_color][name][1])
                    } else {
                        return "#fff";
                    }
                })
            }
        }
        $(containerSelector +" input.residue_fill").on('move.spectrum', function () { change_fill(); });

        
        d3.select(containerSelector).select(".residue_border").on("change", function () {
            change_border();
        });
        function change_border() {
            fill_color = $(containerSelector + " #residue_border").val();

            color_id1  = $(containerSelector+" #border_color1").spectrum("get").toHexString();
            color_id2 = $(containerSelector + " #border_color2").spectrum("get").toHexString();
            if ($(containerSelector + " #border_color3").spectrum("get")) {
                color_id3 = $(containerSelector + " #border_color3").spectrum("get").toHexString();
            } else {
                color_id3 = false;
            }
            console.log('change border color to!', fill_color, color_id1, color_id2, color_id3);
            
            // if no color 3 make only linear between two colors.
            if (color_id3) {
                var color_range = d3v4.scaleLinear()
                    .domain([0, 0.5, 1])
                    .range([color_id1, color_id2, color_id3]);
            } else {
                var color_range = d3v4.scaleLinear()
                    .domain([0, 1])
                    .range([color_id1, color_id2]);
            }
            if (fill_color == 'none') {
                node.select("rect").style("stroke", "#000");
                node.select("circle").style("stroke", "#000");
            } else if (fill_color == 'solid') {
                node.select("rect").style("stroke", color_id1);
                node.select("circle").style("stroke", color_id1);
            } else {
                // loop over nodes, find name
                node.select("circle").style("stroke", function (d) {
                    name = d.name;
                    if (colors_gn[fill_color][name]) {
                        return color_range(colors_gn[fill_color][name][1]);
                    } else {
                        return color_range(0);
                        // return "#000";
                    }
                })
            }
        }
        $(containerSelector +" input.residue_border").on('move.spectrum', function () { change_border(); });

        d3.select(containerSelector).select(".residue_text").on("change", function () {
            change_text();
        });
        function change_text() {
            fill_color = $(containerSelector + " #residue_text").val();

            color_id1  = $(containerSelector+" #text_color1").spectrum("get").toHexString();
            color_id2 = $(containerSelector + " #text_color2").spectrum("get").toHexString();
            if ($(containerSelector + " #text_color3").spectrum("get")) {
                color_id3 = $(containerSelector + " #text_color3").spectrum("get").toHexString();
            } else {
                color_id3 = false;
            }
            console.log('change text color to!', fill_color, color_id1, color_id2, color_id3);
            
            // if no color 3 make only linear between two colors.
            if (color_id3) {
                var color_range = d3v4.scaleLinear()
                    .domain([0, 0.5, 1])
                    .range([color_id1, color_id2, color_id3]);
            } else {
                var color_range = d3v4.scaleLinear()
                    .domain([0, 1])
                    .range([color_id1, color_id2]);
            }
            if (fill_color == 'none') {
                node.select("text").style("fill", "#000");
            } else if (fill_color == 'solid') {
                node.select("text").style("fill", color_id1);
            } else {
                // loop over nodes, find name
                node.select("text").style("fill", function (d) {
                    name = d.name;
                    if (colors_gn[fill_color][name]) {
                        return color_range(colors_gn[fill_color][name][1]);
                    } else {
                        return color_range(0);
                        // return "#000";
                    }
                })
            }
        }
        $(containerSelector +" input.residue_text").on('move.spectrum', function () { change_text(); });


        d3.select(containerSelector).select("#colorLinks").on("change", function () {
            colorLinks = d3.select(containerSelector).select("#colorLinks").property("checked");
            
            if (!segment_view) {
                link.style("stroke", function (d) {
                    scale = 0.5 + Math.abs(d.freq)*0.5 / 100;
                    // scale = Math.abs(d.freq) / 100;
                    scale = 0.2 + Math.abs(d.freq) * 0.8 / max_freq;
                    //console.log(max_freq, d.freq, d.links);
                    color = { r: 200, g: 200, b: 200 }; //grey
                    if (d.freq < 0) {
                        // if the header is a set two, then make it red
                        // color = { r: 255*scale, g: (153)*scale, b: (153)*scale }; //red
                        color = { r: 255, g: 255-(255-153)*scale, b: 255-(255-153)*scale }; //red
                    } else if (d.freq > 0) {
                        // Positive numbers are blue either cos they are set 1 or cos "set 1 has most"
                        // This is also used for single set/structure
                        // color = { r: (153)*scale, g: (204)*scale, b: 255*scale }; //blue
                        color = { r: 255-(255-153)*scale, g: 255-(255-204)*scale, b: 255 }; //blue
                    }
                    var hex = rgb2hex(color.r, color.g, color.b);
                    
                    return colorLinks ? hex : "#000";
                });
            } else {
                link.style("stroke", function (d) {
                    scale = 0.5 + Math.abs(d.freq)*0.5 / 100;
                    // scale = Math.abs(d.freq) / 100;
                    scale = Math.abs(d.sfreq1) / max_freq;
                    color = { r: 255-(255-153)*scale, g: 255-(255-204)*scale, b: 255 }; //blue
                    var hex = rgb2hex(color.r, color.g, color.b);
                    return colorLinks ? hex : "#000";
                });

            }
            link.style("stroke-width", function (d) {
                scale = 0.5 + Math.abs(d.freq)*0.5 / 100;
                // scale = Math.abs(d.freq) / 100;
                scale = segment_view ? Math.abs(d.sfreq1) / max_freq : Math.abs(d.freq) / 100;
                
                return colorLinks ? Math.round(20*scale) : d.size || 5;
            });
            link.style("stroke-opacity", function (d) { return colorLinks ? 1 : 0.3; })

            if (segment_view) {
                link1.style("stroke", function (d) {
                    scale = 0.5 + Math.abs(d.freq)*0.5 / 100;
                    // scale = Math.abs(d.freq) / 100;
                    scale = Math.abs(d.sfreq2) / max_freq;
                    color = { r: 255, g: 255-(255-153)*scale, b: 255-(255-153)*scale }; //red
                    var hex = rgb2hex(color.r, color.g, color.b);
                    return colorLinks ? hex : "#000";
                });
                link1.style("stroke-width", function (d) {
                    scale = 0.5 + Math.abs(d.freq)*0.5 / 100;
                    // scale = Math.abs(d.freq) / 100;
                    scale = Math.abs(d.sfreq2) / max_freq;
                    
                    return colorLinks ? Math.round(20*scale) : d.size || 5; ;
                });
                link1.style("stroke-opacity", function (d) { return colorLinks ? 1 : 0; })
                if (colorLinks) {
                    svg.style("background-color", "#fff");
                    labelParent.attr("dy",function(d,i) { stroke_width = Math.round(20*Math.abs(d.sfreq1) / max_freq); return  stroke_width ? -stroke_width/2 : -5/2;})
                    labelParent.attr("dy",5)
                } else {
                    svg.style("background-color", "#fff");
                    labelParent.attr("dy", function (d, i) { return d.size ? -d.size / 2 : -5 / 2; })
                    labelParent.attr("dy",5)
                }
            }
            simulation.alpha(1).restart();
        });

        d3.select(containerSelector).select("#set_filter").on("change", function () {
            // filter_sets = d3v4.select(containerSelector).select("#set_filter").val();
            filter_sets = $(containerSelector+" #set_filter").val();
            console.log('changing filtering!', filter_sets);
            prepare_data(false);
            init();
        });
        // init option values

        // console.log("set link_strength_change to", link_strength);
        // d3v4.select(containerSelector).select("#link_strength_change").property("value", link_strength);
        console.log("set link_distance_change to", link_distance);
        d3v4.select(containerSelector).select("#link_distance_change").property("value", link_distance);
        console.log("set charge_change to", charge);
        d3v4.select(containerSelector).select("#charge_change").property("value", -charge);
        console.log("set gravity to", gravity);
        d3v4.select(containerSelector).select("#gravity_change").property("value", gravity);
        console.log("set collide to", collide);
        d3v4.select(containerSelector).select("#collide_change").property("value", gravity);
        console.log("set set_filter to", filter_sets);
        d3v4.select(containerSelector).select("#set_filter").property("value", filter_sets);
        console.log("set colorLinks to", colorLinks);
        d3.select(containerSelector).select("#colorLinks").property("checked",colorLinks);

        // init colors from begining
        d3v4.select(containerSelector).select("#colorLinks").dispatch("change");

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

    function network(data, prev, index, expand) {        expand = expand || {};
        var gm = {},    // group map
            nm = {},    // node map
            lm = {},    // link map
            gn = {},    // previous group nodes
            gc = {},    // previous group centroids
            nodes = [], // output nodes
            links = []; // output links
        max_freq = 0; //reset max_freq
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
                l = lm[i] || (lm[i] = {source:u, target:v, size:0, freq:0, sfreq1:0, sfreq2:0});
            l.size += 1;
            l.freq += e.freq;
            l.sfreq1 += e.sfreq1>e.sfreq2 ? 1 : 0; // Change 1 to e.freq to 'nuance' the number
            l.sfreq2 += e.sfreq2>e.sfreq1 ? 1 : 0;
        }
        for (i in lm) {
            // lm[i].freq /= lm[i].size; // Don't do this as the sum is more "meaningful"
            // if (Math.abs(lm[i].freq) > max_freq) max_freq = Math.abs(lm[i].freq);
            if (lm[i].sfreq1 > max_freq) max_freq = lm[i].sfreq1;
            if (lm[i].sfreq2 > max_freq) max_freq = lm[i].sfreq2;
            links.push(lm[i]);
        }
        console.log('max_freq', max_freq);
        return {nodes: nodes, links: links};
    }
        
    function getGroup(n) { return n.group; }

    function assignPosition(segment) {

        switch (segment) {
        case "TM1":
                x = 323;
                y = 170;
            break;
        case "TM2":
                x = 248;
                y = 63;
            break;
        case "TM3":
                x = 165;
                y = 177;
            break;
        case "TM4":
                x = 58;
                y = 82;
            break;
        case "TM5":
                x = 40;
                y = 186;
            break;
        case "TM6":
                x = 83;
                y = 293;
            break;
        case "TM7":
                x = 250;
                y = 285;
            break;
        case "H8":
                x = 358;
                y = 292;
            break;
        case "ECL1":
                x = 164;
                y = 26;
            break;
        case "ECL2":
                x = 17;
                y = 125;
            break;
        case "ICL1":
                x = 383;
                y = 76;
            break;
        case "ICL2":
                x = 95;
                y = 35;
            break;
        default:
            x = 0;
            y = 0;
        }
        return [x,y];
    }

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
            size = 20;
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

        $(containerSelector).find(".controls-panel").remove();

        newDiv.setAttribute("class", "controls-panel");

        var content = '<div class="controls">'
        //                                  +'<h4>Controls</h4>';

        var select_data_options = ''
        $.each(nice_index_names, function (key, description) {
            select_data_options += '<option value="' + key + '">' + description + '</option>';
        });

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
        content = '<span class="pull-right network_controls_toggle" style="cursor: pointer;"><span class="glyphicon glyphicon-option-horizontal btn-download png"></span></span><span class="options" style="display: block; min-width: 120px;">';
            // '<input id="checkCurved" type="checkbox" checked>' +
            // '<span class="checkboxtext"> Curved links' +
            // '</span>' +
            // '</input>' +
            // '<br><button class="btn btn-primary btn-xs" id="resetfixed">Release fixed</button>' +
            // '<br><button class="btn btn-primary btn-xs" id="freeze">Freeze all</button>' +
            
            // 'Stick drag <input id="stickyDrag" type="checkbox" checked><br>' +
            
            
            // 'Node color <select id="node_color"><option value="segment">Segment</option><option value="grey">Grey</option><option value="white">White</option></select><br>' +

        content += '<Strong>Options</strong><br>';

        // common options

        content += 'Highlight res <input id="highlightNode" type="checkbox" ><br>';
            
        if (mode != "single-crystal") {
            content += 'Color links by frequency <input id="colorLinks" type="checkbox" checked><br>' +
            '% of kept contacts <input id="change_to_freq" type="checkbox"><br>';
        }
        
        if (mode == "two-crystal-groups") {
            content += 'Filter <select id="set_filter"><option value="all">All</option><option value="set1">Set1</option><option value="set2">Set2</option><option value="both">Both</option></select><br>';
        }
                   
        if (segment_view) {
            content += 'Link Label <input id="linkLabel" type="checkbox" checked><br>' +
                'Segment fill:<select id="node_color" class="residue_fill snakeplot_property_select">' +
                '<option value="segment">Segment</option><option value="solid">Color</option></select><input type="text" id="fill_color1" class="togglePaletteOnly_red residue_fill" value="red" /><br>';
                
        } else {
            content += 'Add consensus AA<input id="addAA" type="checkbox"><br>';
        
            content += '<table><tr><th>Area</th><th>Property</th><th>Color1</th><th>Color2</th><th>Color3</th></tr>';

            content += '<tr><td>Residue fill:</td><td><select id="node_color" class="residue_fill snakeplot_property_select">' +
                '<option value="segment">Segment</option>' +
                select_data_options +
                '</select></td>' +
                '<td><input type="text" id="fill_color1" class="togglePaletteOnly_red residue_fill" value="red" /></td>' +
                '<td><input type="text" id="fill_color2" class="togglePaletteOnly_blue residue_fill" value="blue" /></td>' +
                '<td><input type="text" id="fill_color3" class="togglePaletteOnly_empty residue_fill" value="" /></td>' +
                '</tr>'
                ;
            
            content += '<tr><td>Residue border:</td><td><select id="residue_border" class="residue_border snakeplot_property_select">' +
                '<option value="none">None</option>' +
                select_data_options +
                '</select></td>' +
                '<td><input type="text" id="border_color1" class="togglePaletteOnly_red residue_border" value="red" /></td>' +
                '<td><input type="text" id="border_color2" class="togglePaletteOnly_blue residue_border" value="blue" /></td>' +
                '<td><input type="text" id="border_color3" class="togglePaletteOnly_empty residue_border" value="" /></td>' +
                '</tr>'
                ;
            
            content += '<tr><td>Residue text:</td><td><select id="residue_text" class="residue_text snakeplot_property_select">' +
                '<option value="none">None</option>' +
                select_data_options +
                '</select></td>' +
                '<td><input type="text" id="text_color1" class="togglePaletteOnly_red residue_text" value="red" /></td>' +
                '<td><input type="text" id="text_color2" class="togglePaletteOnly_blue residue_text" value="blue" /></td>' +
                '<td><input type="text" id="text_color3" class="togglePaletteOnly_empty residue_text" value="" /></td>' +
                '</tr>'
                ;
            content += '<tr><td>Residue rotation:</td><td colspan=4><select id="residue_rotation" class="residue_rotation snakeplot_property_select">' +
                '<option value="none">None</option>' +
                select_data_options +
                '</select></td>' 
                '</tr> '
            // content += '<tr><td>Residue border:</td><td><select id="snakeplot_color_border" class="residue_border snakeplot_property_select">' +
            //     '<option value="none">None</option>' +
            //     select_data_options +
            //     '</select></td><td>' +
            //     '<select id=border_color1 class="border_color residue_border snakeplot_color_select">' +
            //     select_color_options_white +
            //     '</select></td><td>' +
            //     '<select id=border_color2 class="border_color residue_border snakeplot_color_select">' +
            //     select_color_options_red +
            //     '</select></td><td>' +
            //     '<select id=border_color3 class="border_color residue_border snakeplot_color_select">' +
            //     '<option value="none">None</option>' +
            //     select_color_options +
            //     '</select></td>' +
            //     '<td>' + 
            //     '<div class="btn-group btn-toggle residue_border" id="border_filtered">' +
            //     '    <button class="btn btn-xs btn-primary active" value="true">Kept</button>' +
            //     '    <button class="btn btn-xs btn-default" value="false">All</button>' +
            //     '</div>' +
            //     '</td></tr> '
            //     ;
            // content += '<tr><td>Border thickness:</td><td><select id="snakeplot_border_stroke" class="residue_border">' +
            //     '<option>1</option>' +
            //     '<option>2</option>' +
            //     '<option SELECTED>3</option>' +
            //     '<option>4</option>' +
            //     '<option>5</option>' +
            //     '</td></tr>';
            // content += '<tr><td>Residue text:</td><td><select id="snakeplot_color_text" class="residue_text snakeplot_property_select">' +
            //         '<option value="none">None</option>' +
            //         select_data_options +
            //         '</select></td><td>' +
            //         '<select id=text_color1 class="text_color residue_text snakeplot_color_select">' +
            //         select_color_options_white +
            //         '</select></td><td>' +
            //         '<select id=text_color2 class="text_color residue_text snakeplot_color_select">' +
            //         select_color_options_red +
            //         '</select></td><td>' +
            //         '<select id=text_color3 class="text_color residue_text snakeplot_color_select">' +
            //         '<option value="none">None</option>' +
            //         select_color_options +
            //         '</select></td>' +
            //         '<td>' + 
            //         '<div class="btn-group btn-toggle residue_text" id="text_filtered">' +
            //         '    <button class="btn btn-xs btn-primary active" value="true">Kept</button>' +
            //         '    <button class="btn btn-xs btn-default" value="false">All</button>' +
            //         '</div>' +
            //         '</td></tr> '
            //         ;
            // content += '<tr><td>Residue rotation:</td><td colspan=4><select id="snakeplot_color_rotation" class="residue_rotation snakeplot_property_select">' +
            //     '<option value="none">None</option>' +
            //     select_data_options +
            //     '</select></td>' +
            //     '<td>' + 
            //     '<div class="btn-group btn-toggle residue_rotation" id="rotation_filtered">' +
            //     '    <button class="btn btn-xs btn-primary active" value="true">Kept</button>' +
            //     '    <button class="btn btn-xs btn-default" value="false">All</button>' +
            //     '</div>' +
            //     '</td></tr> '

            content += "</table>";     
        
        }    
        content += '<strong>Network settings</strong> [TODO: Show/hide button]<br><table><tr><td>Link Strength</td><td><input id="link_strength_change" style="width:80px;" type="range" min="0" max="1" step="any" value="0.5"></td>' +
            '<td>Link Distance</td><td><input id="link_distance_change" style="width:80px;" type="range" min="0" max="200" step="any" value="40"></td></tr>' +
            '<tr><td>Charge</td><td><input id="charge_change" style="width:80px;" type="range" min="0" max="1400" step="any" value="700"></td>' +
            '<td>Gravity</td><td><input id="gravity_change" style="width:80px;" type="range" min="0" max="1" step="any" value="0.1"></td></tr>' +
            '<tr><td>Collide</td><td><input id="collide_change" style="width:80px;" type="range" min="0" max="200" step="any" value="30"></td> ' +
            '<td>Space</td><td><input id="size_change" style="width:80px;" type="range" min="200" max="4000" step="any" value="'+w+'"></td></tr></table>' +
            '</span>';
        // content = '';
        newDiv.innerHTML = content;

        $(containerSelector).append(newDiv);

        color_palette = [
            ["#000", "#444", "#666", "#999", "#ccc", "#eee", "#f3f3f3", "#fff"],
            ["#f00", "#f90", "#ff0", "#0f0", "#0ff", "#00f", "#90f", "#f0f"],
            ["#f4cccc", "#fce5cd", "#fff2cc", "#d9ead3", "#d0e0e3", "#cfe2f3", "#d9d2e9", "#ead1dc"],
            ["#ea9999", "#f9cb9c", "#ffe599", "#b6d7a8", "#a2c4c9", "#9fc5e8", "#b4a7d6", "#d5a6bd"],
            ["#e06666", "#f6b26b", "#ffd966", "#93c47d", "#76a5af", "#6fa8dc", "#8e7cc3", "#c27ba0"],
            ["#c00", "#e69138", "#f1c232", "#6aa84f", "#45818e", "#3d85c6", "#674ea7", "#a64d79"],
            ["#900", "#b45f06", "#bf9000", "#38761d", "#134f5c", "#0b5394", "#351c75", "#741b47"],
            ["#600", "#783f04", "#7f6000", "#274e13", "#0c343d", "#073763", "#20124d", "#4c1130"] //
        ];



        $(containerSelector+" .togglePaletteOnly_red").spectrum({
            showPaletteOnly: true,
            togglePaletteOnly: true,
            hideAfterPaletteSelect:true,
            togglePaletteMoreText: 'more',
            togglePaletteLessText: 'less',
            color: 'red',
            palette: color_palette
        });
        $(containerSelector+" .togglePaletteOnly_blue").spectrum({
            showPaletteOnly: true,
            togglePaletteOnly: true,
            hideAfterPaletteSelect:true,
            togglePaletteMoreText: 'more',
            togglePaletteLessText: 'less',
            color: 'blue',
            palette: color_palette
        });
        $(containerSelector+" .togglePaletteOnly_empty").spectrum({
            showPalette: true,
            togglePaletteOnly: true,
            hideAfterPaletteSelect:true,
            togglePaletteMoreText: 'more',
            togglePaletteLessText: 'less',
            allowEmpty:true,
            palette: color_palette
        });
        $(containerSelector).find(".options").toggle();

        $(containerSelector).find(".network_controls_toggle").click(function() {
            $(containerSelector).find(".options").slideToggle();
        });

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
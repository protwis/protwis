function createNetworkPlot(raw_data,width, inputGraph, containerSelector, segment_view = true) {
// https://github.com/d3/d3-3.x-api-reference/blob/master/Force-Layout.md


    circle_size = 20;
    max_link_size = 15;

    // var dr = 4,      // default point radius
    // off = 15,    // cluster hull offset
    // expand = {}, // expanded clusters
    // data, net, force, hullg, hull, linkg, link, nodeg, node;


    // // Populate matrix for interactions between segments
    track_gns = []
    new_data = { "nodes": [], "links": [] }
    $.each(raw_data['interactions'], function (i, v) {
        if (!filtered_gn_pairs.includes(i)) return;
        gns = separatePair(i);
        seg1 = raw_data['segment_map'][gns[0]];
        seg2 = raw_data['segment_map'][gns[1]];

        if (!track_gns.includes(gns[0])) {
            new_data["nodes"].push({ "name": gns[0], "group": seg1 })
            track_gns.push(gns[0]);
        }
        if (!track_gns.includes(gns[1])) {
            new_data["nodes"].push({"name":gns[1], "group":seg2})
            track_gns.push(gns[1]);
        }

        source = new_data["nodes"].filter(obj => {return obj.name === gns[0]})[0];
        target = new_data["nodes"].filter(obj => {return obj.name === gns[1]})[0];

        new_data["links"].push({ "source": source, "target": target, "value": 1 })
        // if (seg1 != seg2) {
        //     new_data["links"].push({ "source": gns[0], "target": gns[1], "value": 1 })
        // }
    });

    if (!segment_view) {
        // resize by number of nodes if not segments
        console.log(new_data["nodes"].length, "nodes");
        width = width*Math.sqrt(new_data["nodes"].length/8);
    } else {
        width = width*0.8;
    }


    var w = width;
    var h = w;
    var height = w;

    var expand = {},net;

    d3.select(containerSelector).style("position","relative");

    div = d3.select(containerSelector).insert("div")
        .attr("class", "network")
        .style("width", "100%")
        .style("-webkit-backface-visibility", "hidden");
    
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

        // console.log(graph);
        // graph['nodes'].forEach(function (n) {
        //     console.log(n, n.size);
        // });

        max_size = Math.max.apply(Math, graph['nodes'].map(function(v) {
            return v.size;
        }));

        max_link = Math.max.apply(Math, graph['links'].map(function(v) {
            return v.size;
        }));
        
        console.log(max_link);

        // var graph = new_data;

        // normalize
        graph['links'].forEach(function (n) {
            console.log(n, n.size);
            n.size = Math.round(max_link_size*n.size/max_link)
            console.log(n, n.size);
        });
    }

    var force = d3.layout.force()
    .size([w, h])
    .gravity(0.2)
    .charge(-1000)
    .linkDistance(100)
    // .friction(0.5)
    .on("tick", tick);

    var drag = force.drag()
        .on("dragstart", dragstart);

    var link = svg.selectAll(".link"),
        node = svg.selectAll(".node");
    
    force
        .nodes(graph.nodes)
        .links(graph.links)
        .start();

    link = link.data(graph.links)
        .enter().append("line")
        .attr("class", "link")
        .style("stroke-width", function(d) { return d.size || 1; });
    
    
    node = svg.selectAll(".node")
        .data(graph.nodes)
        .enter()
        .append("g")
        .on("dblclick", dblclick)
        .call(drag);

                
    node.append("circle")
        // .attr("class", function(d) { return "node" + (d.size?"":" leaf"); })
        .attr("class", "node")
        .attr("r", function(d) { return d.size ? assignSize(d.group) : 20; })
        // .attr("r", 20)
        .style("fill", function (d) { return assignRainbowColor(d.group); })
        
    node.append("text")
        .attr("dy", ".35em")
        .attr("text-anchor", 'middle')
        .attr("fill", "#000")
        .attr("font-size", function(d) { return d.size ? assignSize(d.group)-6 : 12; })
        .text(function(d){return d.size?d.group:d.name});

    function tick() {
        link.attr("x1", function(d) { return d.source.x; })
            .attr("y1", function(d) { return d.source.y; })
            .attr("x2", function(d) { return d.target.x; })
            .attr("y2", function(d) { return d.target.y; });
        
        node.attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; });
    }

    function dblclick(d) {
    d3.select(this).classed("fixed", d.fixed = false);
    }

    function dragstart(d) {
    d3.select(this).classed("fixed", d.fixed = true);
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
    
    function separatePair(stringPair) {
        var regex = /([0-9x]+),([0-9x]+)/g;
        var m;
    
        matches = regex.exec(stringPair);
    
        return [matches[1], matches[2]];
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
    
}
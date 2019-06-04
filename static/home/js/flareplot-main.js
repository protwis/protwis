
/**
 * Creates a flareplot svg in the container matched by the `containerSelector`, sets its width to `width` pixels,
 * and associates it with the contents of `json`
 * @param width
 * @param inputGraph
 * @param containerSelector
 * @returns {{getNumFrames, setFrame, framesIntersect, framesSum, setTrack, setTree, getTreeNames, getTrackNames, addNodeToggleListener, addNodeHoverListener, addEdgeToggleListener, addEdgeHoverListener, graph}}
 */
function createFlareplot(width, inputGraph, containerSelector){
    var w = width;
    var h = w;
    var rx = w * 0.5;
    var ry = w * 0.5;
    var rotate = 0;
    var discRad = 55;

    if( typeof inputGraph == "string" ){
        inputGraph = JSON.parse(inputGraph);
    }

    var svg;
    var div;
    var bundle;
    var line;
    var nodes;
    var splines;
    var links;
    var graph;

    var selectedTree = 0;
    var selectedTrack = 0;
    var toggledNodes = {};
    var visibleEdges = [];
    var splineDico;

    return (function() {

        function create_bundle() {
            cluster = d3.layout.cluster()
                .size([360, ry - discRad])
                .sort(function(a, b) {

                    var structures = ["N-term", "TM1", "ICL1", "TM2", "ECL1", "TM3", "ICL2", "TM4", "ECL2", "TM5", "ICL3", "TM6", "ECL3", "TM7", "H8", "C-term"];                    
                    var startWithNonDigit = /^\D/;

                    if (startWithNonDigit.test(a.key) || startWithNonDigit.test(b.key)){
                        var aRes = structures.indexOf(a.key);
                        var bRes = structures.indexOf(b.key);
                    } else {
                        var aRes = a.key.match(/[0-9]*$/);
                        var bRes = b.key.match(/[0-9]*$/);
                        if (aRes.length==0 || bRes.length==0) {
                            aRes = a.key;
                            bRes = b.key;
                        } else {
                            aRes = parseInt(aRes[0]);
                            bRes = parseInt(bRes[0]);
                        }
                    }    
                    return d3.ascending(aRes, bRes);
                });

            graph = preprocessGraph(inputGraph);
            nodes = cluster.nodes(graph.trees[selectedTree].tree[""]);
            bundle = d3.layout.bundle();

            //links = graph.trees[selectedTree].frames;
            //splines = bundle(links[0]);
            //splineDico = buildSplineIndex(splines);

            links = graph.trees[selectedTree].allEdges;
            splines = bundle(links);
            splineDico = buildSplineIndex(splines);

            line = d3.svg.line.radial()
                .interpolate("bundle")
                .tension(0.85)
                .radius(function(d) { return d.y; })
                .angle(function(d) { return d.x / 180 * Math.PI; });

            d3.select(containerSelector).style("position","relative");

            div = d3.select(containerSelector).insert("div")
                .style("width", w + "px")
                .style("height", h + "px")
                .style("-webkit-backface-visibility", "hidden");

            svg = div.append("svg:svg")
                .attr("width", w)
                .attr("height", h)
                .append("svg:g")
                .attr("transform", "translate(" + rx + "," + ry + ")");

            //// Find the width of the node-name track. Temporarily add all text, go through them and get max-width
            //var tmpTexts = svg.selectAll("g.node")
            //    .data(nodes.filter(function(n) { return !n.children; }), function(d) { return d.key; })
            //    .enter().append("svg:g")
            //    .attr("class", "node")
            //    .attr("id", function(d) { return "node-" + d.key; })
            //    .append("text")
            //    .text(function(d) { return d.key; });
            //var maxTextWidth = d3.max(svg.selectAll("text")[0], function(t){ return t.getBBox().width; });
            //svg.selectAll("g.node").remove();


            var path = svg.selectAll("path.link")
                .data(links, function(d,i){
                    var key = "source-" + d.source.key + "target-" + d.target.key;
                    return key;
                })
                .enter().append("svg:path")
                .attr("class", function(d) {
                    var ret = "link source-" + d.source.key + " target-" + d.target.key;
                    if( d.source.key in toggledNodes || d.target.key in toggledNodes) {
                        ret += " toggled";
                    }
                    return ret;
                })
                .style("stroke-width",function(d){
                    return 0;
                })
                .style("stroke",function(d){ return d.color; })
                .style("fill","none")
                .attr("d", function(d, i) { return line(splines[i]); })
                .on("mouseover", function(d){ fireEdgeHoverListeners(d); })
                .on("click", function(d){ fireEdgeToggleListeners(d); });






            svg.selectAll("g.node")
                .data(nodes.filter(function(n) { return !n.children; }), function(d) { return d.key; })
                .enter().append("svg:g")
                .attr("class", "node")
                .attr("id", function(d) { return "node-" + d.key; })
                .attr("transform", function(d) { return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")"; })
                .append("text")
                .attr("dx", function(d) { return d.x < 180 ? 8 : -8; })
                .attr("dy", ".31em")
                .attr("text-anchor", function(d) { return d.x < 180 ? "start" : "end"; })
                .attr("transform", function(d) { return d.x < 180 ? null : "rotate(180)"; })
                .text(function(d) { return d.key; })
                .on("mouseover", mouseoverNode)
                .on("mouseout", mouseoutNode)
                .on("click", function(d){ toggleNode(d.name); });


            var arcW = 250.0/(graph.nodeNames.length)*Math.PI/360;
            var arc = d3.svg.arc()
                .innerRadius(ry-15)
                .outerRadius(function(d,i){
                  var sz = d.size;
                  if(!sz) { sz = 0.0; }
                  var or = ry-15+sz*15;
                  return or;
                })
                .startAngle(-arcW)
                .endAngle(arcW);

            svg.selectAll("g.trackElement")
                .data(graph.tracks[selectedTrack].trackProperties, function(d){ return d.nodeName; })
                .enter().append("svg:g")
                .attr("class", "trackElement")
                .attr("id", function(d) { return "trackElement-" + d.nodeName; })
                .append("path")
                .attr("transform", function(d) {
                    var x = graph.trees[selectedTree].tree[d.nodeName].x;
                    return "rotate("+x+")" ;
                })
                .style("fill", function(d){ return d.color; })
                .attr("d", arc)
                .on("click", function(d){
                    //Locate corresponding node
                    nodes.filter(function(n){ return n.name==d.nodeName; })
                        .forEach(function(n){ toggleNode(n.name); });

                });

            setFrame(0);
        }

        /**
         * Preprocess the graph and return its reference but with additional fields:
         * - nodeMap: a map associating each node name (string) with a node reference
         * - nodeNames: list of string-representation of nodes in no particular order
         * - edges: list of (not necessarily distinct) interactions with interaction timepoints (frames)
         * - tracks: ...
         * - trees: ...
         */
        function preprocessGraph(graph) {
            "use strict";
            //Create defaults, trees, and tracks if they don't exist
            if (!graph.defaults) {
                graph.defaults = {};
            }
            if (!graph.trees) {
                graph.trees = [{"treeLabel": "default", "treePaths": []}];
            }
            if (!graph.tracks) {
                graph.tracks = [{"trackLabel": "default", "trackProperties": []}];
            }

            // =========== Construct `nodeNames` list =========== \\
            graph.nodeNames = [];

            //Fill nodeNames from edges
            graph.edges.forEach(function (e) {
                if (graph.nodeNames.indexOf(e.name1) == -1) {
                    graph.nodeNames.push(e.name1);
                }
                if (graph.nodeNames.indexOf(e.name2) == -1) {
                    graph.nodeNames.push(e.name2);
                }
            });

            //Fill nodeNames from trees
            graph.trees.forEach(function (t) {
                t.treePaths.forEach(function (p) {
                    var name = p.substring(p.lastIndexOf(".") + 1);
                    if (!(graph.nodeNames.indexOf(name) > -1)) {
                        graph.nodeNames.push(name);
                    }
                });
            });

            //Fill nodeNames from tracks
            graph.tracks.forEach(function (t) {
                t.trackProperties.forEach(function (c) {
                    var name = c.nodeName;
                    if (!(graph.nodeNames.indexOf(name) > -1)) {
                        graph.nodeNames.push(name);
                    }
                });
            });

            // =========== Parse `edges` section ========== \\

            //Go through edges and ensure that default widths have been assigned
            graph.edges.forEach(function (e) {
                e.width = e.width || graph.defaults.edgeWidth || 1;
            });


            // =========== Parse `trees` section ========== \\
            function addToMap(nodeMap, fullName) {
                var i = fullName.lastIndexOf(".");
                var name = fullName.substring(i + 1);
                var node = nodeMap[name];
                if (!node) {
                    node = {name: name, children: []};
                    nodeMap[name] = node;
                    if (name.length) {
                        node.parent = addToMap(nodeMap, fullName.substring(0, i));
                        node.parent.children.push(node);
                        node.key = name;
                    }
                }
                return node;
            }

            graph.trees.forEach(function (t) {
                var addedNames = [];
                //Ensure that each tree-object has a `tree` attribute with the hierarchy
                t.tree = {};
                t.treePaths.forEach(function (p) {
                    addToMap(t.tree, p);
                    addedNames.push(p.substring(p.lastIndexOf(".") + 1));
                });

                //Ensure that even nodes not mentioned in the treePaths are added to the tree
                graph.nodeNames.forEach(function (p) {
                    if (addedNames.indexOf(p) === -1) {
                        addToMap(t.tree, p);
                    }
                });
            });


            //Go through graph.edges and convert name1, name2, and frames to target and source object arrays.
            graph.trees.forEach(function (t) {
                t.frames = [];
                var summaryEdges = {};
                t.allEdges = [];
                graph.edges.forEach(function (e) {
                    //Set source and target of edge
                    var edge = {
                        source: t.tree[e.name1],
                        target: t.tree[e.name2],
                        key: "" + t.tree[e.name1].key + "-" + t.tree[e.name2].key,
                        color: e.color || graph.defaults.edgeColor || "rgba(100,100,100)",
                        opacity: e.opacity || graph.defaults.edgeOpacity || 1,
                        width: e.width || graph.defaults.edgeWidth || 1
                    };

                    var edgeKey = edge.key;
                    if (!summaryEdges[edgeKey]) {
                        summaryEdges[edgeKey] = {
                            source: edge.source,
                            target: edge.target,
                            key: edge.key,
                            color: edge.color,
                            opacity: edge.opacity,
                            width: edge.width
                        };
                        t.allEdges.push({
                            source: edge.source,
                            target: edge.target,
                            key: edge.key,
                            color: edge.color,
                            opacity: edge.opacity,
                            width: edge.width
                        });
                    } else {
                        summaryEdges[edgeKey].width += edge.width;
                    }

                    //edge.source = t.tree[edge.name1]; console.assert(edge.source);
                    //edge.target = t.tree[edge.name2]; console.assert(edge.target);
                    //edge.key = ""+i;


                    //Add interaction frames
                    e.frames.forEach(function (f) {
                        while (t.frames.length <= f) {
                            t.frames.push([]);
                        }
                        t.frames[f].push(edge);
                    });
                });
                t.summaryEdges = summaryEdges;
            });

            // =========== Parse `tracks` section ========== \\

            graph.tracks.forEach(function (track) {
                //Ensure that nodes not mentioned in the trackProperties are created with default values
                var remainingNodeNames = [];
                graph.nodeNames.forEach(function (n) {
                    remainingNodeNames.push(n);
                });

                track.trackProperties.forEach(function (p) {
                    //Remove p.name from remainingNodeNames
                    var idx = remainingNodeNames.indexOf(p.nodeName);
                    if (idx > -1) {
                        remainingNodeNames.splice(idx, 1);
                    }
                });

                remainingNodeNames.forEach(function (n) {
                    var color = graph.defaults.trackColor || "white";
                    var size = graph.defaults.trackSize || 0;
                    track.trackProperties.push({"nodeName": n, "color": color, "size": size});
                });
            });

            // =========== Parse `defaults` section ========== \\

            //From https://stackoverflow.com/questions/566203/changing-css-values-with-javascript
            function insertCSSrule(selector, property, value) {
                for (var i=0; i<document.styleSheets.length;i++) {//Loop through all styles
                    try {
                        document.styleSheets[i].insertRule(selector+ ' {'+property+':'+value+'}', document.styleSheets[i].cssRules.length);
                    } catch(err) {//IE
                        try {
                            document.styleSheets[i].addRule(selector, property+':'+value);
                        } catch(err) {}
                    }
                }
            }

            if (graph.defaults.edgeOpacity) {
                insertCSSrule(".link", "stroke-opacity", graph.defaults.edgeOpacity);
                insertCSSrule(".link", "opacity", graph.defaults.edgeOpacity);
            }

            return graph;
        }


        /**
         *
         * @param clusterDefinition {}
         *
         * keys are the key of the cluster
         * values are array that correspond to node keys
         *
         */
        function assignCluster(clusterDefinition, oldCluster, graph) {

            var nodesMap = clusterDefinition.tree;
            var root = nodesMap[""];
            var rootNodes = root.children;

            // recursively copy x and y propery from the old cluster
            rootNodes.forEach(copyAndGoThruChildren);
            function copyAndGoThruChildren(node) {

                var newNode;
                var nodeKey = node.key;
                if   (nodesMap[nodeKey]) {
                    newNode = nodesMap[nodeKey];
                    var oldNode =  oldCluster.tree[nodeKey];
                    if (oldNode) {
                        newNode.oldX = oldNode.x;
                        newNode.oldY = oldNode.y;
                    }
                } else
                {
                    // it could happen that an new node come (in case of intermediate level)
                    newNode = nodesMap[nodeKey];
                    newNode.clusterName = nodeKey;

                }
                if (newNode.children && newNode.children.length > 0) {
                    newNode.children.forEach(copyAndGoThruChildren);
                }
            }
        }


        /**
         * Find the total number of frames in the graph.
         * @returns {number}
         */
        function getNumFrames(){
            var maxFrame = -1;
            graph.edges.forEach(function(e){
                maxFrame = Math.max(maxFrame,  Math.max.apply(Math, e.frames));
            });
            return maxFrame+1;
        }

        /**
         * Change the state of the flareplot so it reflects the interactions in the indicated frame.
         * @param frameNum a number indicating the frame to set.
         */
        function setFrame(frameNum){
            frameNum = parseInt(frameNum);
            rangeSum(frameNum,frameNum+1);
        }

        /**
         * Update the state of the flareplot so it reflects the intersection over a range.
         * @param rangeStart first frame to include (must be less than `rangeEnd`)
         * @param rangeEnd first frame after `rangeStart` that should not be included
         */
        function rangeIntersect(rangeStart, rangeEnd){
            //splines = bundle(links);
            //splineDico = buildSplineIndex(splines);
            var path = svg.selectAll("path.link");

            path.style("stroke-width",
                function(d,i){
                    var count = graph.edges[i].frames.rangeCount(rangeStart, rangeEnd-1);
                    return count==(rangeEnd-rangeStart)?(2 * graph.edges[i].width):0;
                })
                .attr("class", function(d) {
                    var ret = "link source-" + d.source.key + " target-" + d.target.key;
                    if( d.source.key in toggledNodes || d.target.key in toggledNodes) {
                        ret += " toggled";
                    }
                    return ret;
                })
                //.attr("class", function(d) { return "link source-" + d.source.key + " target-" + d.target.key; })
                .style("stroke",function(d){ return d.color; })
                .attr("d", function(d, i) { return line(splines[i]); });

        }

        /**
         * Update the state of the flareplot so it reflects the sum over a range.
         * @param rangeStart first frame to include (must be less than `rangeEnd`)
         * @param rangeEnd first frame after `rangeStart` that should not be included
         */
        function rangeSum(rangeStart, rangeEnd){
            //splines = bundle(links);
            //splineDico = buildSplineIndex(splines);
            var path = svg.selectAll("path.link");

            var widthScale = d3.scale.linear()
                .domain([1,Math.max(1,rangeEnd-rangeStart)])
                .range([2,10]);

            visibleEdges = [];

            path.style("stroke-width",
                function(d,i){
                    var count = graph.edges[i].frames.rangeCount(rangeStart, rangeEnd-1);
                    if (count>0){
                        var e = {edge:graph.edges[i], weight:count/(rangeEnd-rangeStart)};
                        e.toggled = e.edge.name1 in toggledNodes || e.edge.name2 in toggledNodes;
                        visibleEdges.push(e);

                        return widthScale(count) * graph.edges[i].width;
                    } else {
                        return 0;
                    }
                })
                .attr("class", function(d) {
                    var ret = "link source-" + d.source.key + " target-" + d.target.key;
                    if( d.source.key in toggledNodes || d.target.key in toggledNodes) {
                        ret += " toggled";
                    }
                    return ret;
                })
                //.attr("class", function(d) { return "link source-" + d.source.key + " target-" + d.target.key; })
                .style("stroke",function(d){ return d.color; })
                .attr("d", function(d, i) { return line(splines[i]); });

            // fireFrameListeners(visibleEdges);
            fireFrameListeners({type:"sum", start:rangeStart, end:rangeEnd});
        }

        /**
         * Update the state of the flareplot so it reflects the intersection over a subset.
         * @param subset a list of numbers indicating which frames to include
         */
        function subsetIntersect(subset, subtract){
            //splines = bundle(links);
            //splineDico = buildSplineIndex(splines);
            var path = svg.selectAll("path.link");
            visibleEdges = [];

            path.style("stroke-width",
                function(d,i) {
                    for (var c = 0; c < subset.length; c++) {
                        var frame = subset[c];
                        var iud = graph.edges[i].frames.indexUpDown(frame);
                        if (iud[0] != iud[1]) return 0;
                    }

                    if (subtract) {
                        for (var c = 0; c < subtract.length; c++) {
                            var frame = subtract[c];
                            var iud = graph.edges[i].frames.indexUpDown(frame);
                            if (iud[0] == iud[1]) return 0;
                        }
                    }

                    var e = {edge:graph.edges[i], weight:1};
                    e.toggled = e.edge.name1 in toggledNodes || e.edge.name2 in toggledNodes;
                    visibleEdges.push(e);
                    return 2 * e.width;
                })
                .attr("class", function(d) {
                    var ret = "link source-" + d.source.key + " target-" + d.target.key;
                    if( d.source.key in toggledNodes || d.target.key in toggledNodes)
                        ret+=" toggled";
                    return ret;
                })
                .style("stroke",function(d){ return d.color; })
                .attr("d", function(d, i) { return line(splines[i]); });

            fireFrameListeners({type:"intersect", intersected:subset, excluded:subtract});
        }

        /**
         * Update the state of the flareplot so it reflects the sum over a subset.
         * @param subset a list of numbers indicating which frames to include
         */
        function subsetSum(subset){
            //splines = bundle(links);
            //splineDico = buildSplineIndex(splines);
            var path = svg.selectAll("path.link");
            visibleEdges = [];

            var widthScale = d3.scale.linear()
                .domain([1,subset.length])
                .range([2,10]);

            path.style("stroke-width",
                function(d,i){
                    var count = 0;
                    subset.forEach(function(f){
                        var iud = graph.edges[i].frames.indexUpDown(f);
                        if( iud[0] == iud[1] ){ count++; }
                    });
                    var e = {edge:graph.edges[i], weight:count/subset.length};
                    e.toggled = e.edge.name1 in toggledNodes || e.edge.name2 in toggledNodes;
                    visibleEdges.push(e);

                    return count==0?0:(widthScale(count) * e.width);
                })
                .attr("class", function(d) {
                    var ret = "link source-" + d.source.key + " target-" + d.target.key;
                    if( d.source.key in toggledNodes || d.target.key in toggledNodes) {
                        ret += " toggled";
                    }
                    return ret;
                })
                //.attr("class", function(d) { return "link source-" + d.source.key + " target-" + d.target.key; })
                .style("stroke",function(d){ return d.color; })
                .attr("d", function(d, i) { return line(splines[i]); });

            fireFrameListeners({type:"setsum", set:subset});
        }

        /**
         * Update the state of the flareplot so it reflects the intersection over the specified selection. If
         * `selection` and `optionalSelection` are both numbers then the intersection will be taken over the range of
         * frames spanned by the two. If `selection` is an array of numbers the intersection will be taken over the
         * frames in the array. If the arguments don't satisfy these requirements an error is thrown.
         * @param selection
         * @param optionalSelection
         * @returns {*}
         */
        function framesIntersect(selection, optionalSelection) {
            if (typeof selection === "number" && typeof optionalSelection === "number") {
                return rangeIntersect(selection, optionalSelection);
            }
            if (Object.prototype.toString.call(selection) === "[object Array]" && optionalSelection === undefined) {
                return subsetIntersect(selection);
            }
            throw "framesIntersect must take either two integers (range), or an array (subset) as argument";
        }

        function framesIntersectSubtract(intersectSelection, subtractSelection) {
            return subsetIntersect(intersectSelection, subtractSelection);
        }

        /**
         * Update the state of the flareplot so it reflects the sum over the specified selection. If `selection` and
         * `optionalSelection` are both numbers then the intersection will be taken over the range of frames spanned by
         * the two. If `selection` is an array of numbers the intersection will be taken over the frames in the array.
         * If the arguments don't satisfy these requirements an error is thrown.
         * @param selection
         * @param optionalSelection
         * @returns {*}
         */
        function framesSum(selection, optionalSelection) {
            if (typeof selection === "number" && typeof optionalSelection === "number") {
                return rangeSum(selection, optionalSelection);
            }
            if (Object.prototype.toString.call(selection) === "[object Array]" && optionalSelection === undefined) {
                return subsetSum(selection);
            }
            throw "framesSum must take either two integers (range), or an array (subset) as argument";
        }

        function toggleNode(nodeName){
            var svgNodeElement = svg.selectAll("g.node#node-"+nodeName).node();
            var toggled = !d3.select(svgNodeElement).classed("toggledNode");
            d3.select(svgNodeElement)
                .classed("toggledNode", function(){return toggled; });

            var name = nodeName.substring(nodeName.lastIndexOf(".")+1);
            if(!toggled)
                delete toggledNodes[name];
            else
                toggledNodes[name] = "";

            path = svg.selectAll("path.link")
                .classed("toggled", function(d) {
                    return ( d.source.key in toggledNodes || d.target.key in toggledNodes)
                });

            visibleEdges.forEach(function(e){
                if(e.edge.name1==nodeName || e.edge.name2==nodeName){
                    e.toggled = toggled;
                }
            });

            fireNodeToggleListeners(nodeName);
        }


        function mouseoverNode(d) {
            svg.selectAll("path.link.target-" + d.key)
                .classed("target", true)
                .each(updateNodes("source", true));

            svg.selectAll("path.link.source-" + d.key)
                .classed("source", true)
                .each(updateNodes("target", true));
        }

        function mouseoutNode(d) {
            svg.selectAll("path.link.source-" + d.key)
                .classed("source", false)
                .each(updateNodes("target", false));

            svg.selectAll("path.link.target-" + d.key)
                .classed("target", false)
                .each(updateNodes("source", false));
        }

        function updateNodes(name, value) {
            return function(d) {
                //if (value) this.parentNode.appendChild(this);
                svg.select("#node-" + d[name].key).classed(name, value);
            };
        }

        function getTreeNames(){
            var ret = [];
            for (var t=0; t<graph.trees.length; t++ ){
                ret.push(graph.trees[t].treeLabel);
            }
            return ret;
        }

        function setTree(treeIdx){
            var oldTreeIdx = selectedTree;
            selectedTree = treeIdx;
            assignCluster(graph.trees[selectedTree], graph.trees[oldTreeIdx], graph);
            var recomposedSplines = [];
            nodes = cluster.nodes(graph.trees[selectedTree].tree[""]);

            links = Object.values(graph.trees[selectedTree].summaryEdges);

            svg.selectAll("g.node")
                .data(nodes.filter(function(n) { return !n.children; }), function(d) { return d.key})
                .transition().duration(500)
            //.attr("transform", function(d) { return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")"; })
                .attrTween("transform", function(d) {
                    var oldMatrix = "rotate(" + (d.oldX - 90) + ")translate(" + d.y + ")";
                    var newMatrix = "rotate(" + (d.x - 90) + ")translate(" + d.y + ")";
                    return d3.interpolateString(oldMatrix, newMatrix);
                })
                .select("text")
                .attr("dx", function(d) { return d.x < 180 ? 8 : -8; })
                .attr("text-anchor", function(d) { return d.x < 180 ? "start" : "end"; })
                .attr("transform", function(d) { return d.x < 180 ? null : "rotate(180)"; });

            var arcW = 250.0/(graph.nodeNames.length)*Math.PI/360;
            var arc = d3.svg.arc()
                .innerRadius(ry-15)
                .outerRadius(function(d){
                    var sz = d.size;
                    if(!sz) sz = 0.0;
                    return ry-15+sz*15;
                })
                .startAngle(-arcW)
                .endAngle(arcW);

            svg.selectAll("g.trackElement")
                .select("path")
                .transition().duration(500)
                .attrTween("transform", function(d) {
                    var node = graph.trees[selectedTree].tree[d.nodeName];
                    var oldMatrix = "rotate(" + (node.oldX) + ")";
                    var newMatrix = "rotate(" + (node.x) + ")";
                    return d3.interpolateString(oldMatrix, newMatrix);
                })
                .style("fill", function(d){ return d.color; })
                .attr("d", arc);


            // transition the splines
            var newSplines = bundle(Object.values((graph.trees[selectedTree].summaryEdges)));
            var newSplinesDico = buildSplineIndex(newSplines);


            var done = false;
            var path = svg.selectAll("path.link").data(links, function(d){
                return  d.key;
            });

            // i dont understand how d3 orders the spline array, so we need
            path.transition().attrTween("d",
                function(d, i, a) {

                    //if (i != 2) return;
                    // make a copy of the targeted Spline, and put all x to the value of OldX..
                    var oldSpline = [];
                    var key = d.key;

                    var oldSplineIdx =  splineDico[key];
                    var newSplineIdx = newSplinesDico[key];
                    if (oldSplineIdx === void 0 || newSplineIdx === void 0) {
                        console.log("Not found Spline with key", key);
                        return;
                    }

                    for (var j = 0; j < splines[oldSplineIdx].length; j++) {
                        var s = Object.assign({}, splines[oldSplineIdx][j]);
                        oldSpline.push(s);
                    }
                    oldSpline = oldSpline.map(function(s) {
                        return {x: s.x, y: s.y};
                    });

                    var simpleSpline = newSplines[newSplineIdx].map(function(s) {
                        return {x: s.x, y:s.y, key:s.key}
                    });
                    // now if oldspine is missing controlpoints
                    var delta = simpleSpline.length - oldSpline.length;
                    if (oldSpline.length < simpleSpline.length) {
                        //positive delta
                        var recomposedOldSpline = [];
                        // we make the assumption that we start with 3 control points
                        // but they may be more complicated situations
                        // if delta =  2   0 - 0, 1-0, 2-1, 3-2, 4-2  (3 to 5 )
                        // if delta =  4   0-0 1-0 2-0, 3-1, 4-2, 5-2, 6-2  ( 3 to 7 )
                        // if delta = 2 ( 5 to 7) what happens ?
                        // if delta = 4 ( 5 to 9) what happens ?
                        for (i = 0, currentIndex = 0; i < simpleSpline.length; i++) {
                            recomposedOldSpline[i] = oldSpline[currentIndex];
                            if (i <= delta/2 || currentIndex >= oldSpline.length - 1) { } else {
                                currentIndex++;
                            }
                        }
                    } else if (delta < 0) { // (5 < 3)
                        // newer spline has less target point than older spline
                        var recomposedNewSpline = [];
                        // -2, 5 to 3   => 0 -0, 1-0, 2-1, 3-2,4-2  (simplespline 3, oldSpine = 5)
                        // -4 ,7 to 3   => 0-0, 1-0, 2-0, 3-1, 4-2 5-2 6-2
                        delta = Math.abs(delta);
                        for (i = 0, currentIndex = 0; i < oldSpline.length; i++) {
                            recomposedNewSpline[i] = simpleSpline[currentIndex];
                            if (i <= Math.floor(delta / 2) || currentIndex >= simpleSpline.length - 1) {} else {
                                currentIndex++;
                            }
                        }
                        simpleSpline = recomposedNewSpline;
                        recomposedOldSpline = oldSpline;

                    } else
                    {
                        recomposedOldSpline = oldSpline;
                    }
                    recomposedSplines.push(simpleSpline);
                    var interpolate = d3.interpolate(recomposedOldSpline, simpleSpline);
                    // we can update the splines at the next loop, or it will mess D3
                    setTimeout(function(){
                        if (!done){
                            done = true;
                            splines = recomposedSplines;
                            splineDico = buildSplineIndex(recomposedSplines);
                            // we do not want to rebind data here
                        }

                    }, 500);

                    return function(t) {
                        return line(interpolate(t))
                    };
                })
                .duration(500);

        }

        /**
         * For all splines in the array, create an index that match the key of
         * the link and the index in the spline array.
         */
        function buildSplineIndex(splines) {
            var linkKeyToSplineIdx = {};
            splines.forEach(function(spline, idx){
                var source = spline[0].key;
                var target = spline[spline.length-1].key;
                var key = source + "-" + target;
                linkKeyToSplineIdx[key] = idx;
            });
            return linkKeyToSplineIdx;
        }


        function getTrackNames(){
            var ret = [];
            for(var t=0;t<graph.tracks.length;t++){
                ret.push(graph.tracks[t].trackLabel);
            }
            return ret;
        }

        function setTrack(trackIdx){
            selectedTrack = trackIdx;

            var arcW = 250.0/(graph.nodeNames.length)*Math.PI/360;
            var arc = d3.svg.arc()
                .innerRadius(ry-15)
                .outerRadius(function(d){
                    var sz = d.size;
                    if(!sz) sz = 0.0;
                    return ry-15+sz*15;
                })
                .startAngle(-arcW)
                .endAngle(arcW);

            svg.selectAll("g.trackElement")
                .data(graph.tracks[selectedTrack].trackProperties, function(d){ return d.nodeName; })
                .select("path")
                .transition()
                .style("fill", function(d){ return d.color; })
                .attr("d", arc);

        }

        function setTension(tension){
            line.tension(tension);
            var path = svg.selectAll("path.link")
                .attr("d", function(d, i) { return line(splines[i]); })
        }

        function getEdges(){
            return visibleEdges;
        }

        var nodeToggleListeners = [];
        var nodeHoverListeners  = [];
        var edgeToggleListeners = [];
        var edgeHoverListeners  = [];
        var frameListeners  = [];

        function addNodeToggleListener(l){ nodeToggleListeners.push(l); }
        function addNodeHoverListener(l){  nodeHoverListeners.push(l);  }
        function addEdgeToggleListener(l){ edgeToggleListeners.push(l); }
        function addEdgeHoverListener(l){  edgeHoverListeners.push(l);  }
        function addFrameListener(l){  frameListeners.push(l);  }

        function fireNodeToggleListeners(n){ nodeToggleListeners.forEach(function(l){l(n);}); }
        function fireNodeHoverListeners(n){ nodeHoverListeners.forEach(function(l){l(n);}); }
        function fireEdgeToggleListeners(n){ edgeToggleListeners.forEach(function(l){l(n);}); }
        function fireEdgeHoverListeners(n){ edgeHoverListeners.forEach(function(l){l(n);}); }
        function fireFrameListeners(f){ frameListeners.forEach(function(l){l(f);}); }

        create_bundle();

        return {
            getNumFrames: getNumFrames,
            setFrame: setFrame,
            getEdges: getEdges,
            framesIntersect: framesIntersect,
            framesSum: framesSum,
            framesIntersectSubtract: framesIntersectSubtract,
            setTrack: setTrack,
            setTree: setTree,
            getTreeNames: getTreeNames,
            getTrackNames: getTrackNames,
            setTension: setTension,
            addNodeToggleListener: addNodeToggleListener,
            addNodeHoverListener: addNodeHoverListener,
            addEdgeToggleListener: addEdgeToggleListener,
            addEdgeHoverListener: addEdgeHoverListener,
            addFrameListener: addFrameListener,
            graph: graph//, for debugging purposes
        }
    }) ();
}

function upload_button(el, callback) {
    var uploader = document.getElementById(el);
    var reader = new FileReader();

    reader.onload = function(e) {
        var contents = e.target.result;
        callback(contents);
    };

    uploader.addEventListener("change", handleFiles, false);

    function handleFiles() {
        d3.select("#table").text("loading...");
        var file = this.files[0];
        reader.readAsText(file);
    }
}



/**
 * Gets the index of the value just above and just below `key` in a sorted array.
 * If the exact element was found, the two indices are identical.
 */
function indexUpDown(key) {
  "use strict";

  var minIdx = 0;
  var maxIdx = this.length - 1;
  var curIdx, curElm, resIdx;

  while (minIdx <= maxIdx) {
    resIdx = curIdx = (minIdx + maxIdx) / 2 | 0;
    curElm = this[curIdx];

    if (curElm < key)      minIdx = curIdx + 1;
    else if (curElm > key) maxIdx = curIdx - 1;
    else return [curIdx,curIdx];
  }

  return [minIdx,maxIdx];
}

/** Get the number of entries whose value are greater than or equal to `start`
 * and lower than or equal to `end` in a sorted array*/
function rangeCount(start, end){
  var startIdx = this.indexUpDown(start)[0];
  var endIdx   = this.indexUpDown(end)[1];
  return endIdx-startIdx+1;
}

Array.prototype.indexUpDown = indexUpDown;
Array.prototype.rangeCount = rangeCount;

// var list = [1,2,5,10,15,16];
// function testRange(l,s,e, expected){
//   var res = l.rangeCount(s,e);
//     console.log("["+l+"].count("+s+","+e+") -> "+res+" expects "+expected+(res==expected?" PASS":" FAILED"));
// }
//
// testRange(list,   0,  0, 0);
// testRange(list,   0,  1, 1);
// testRange(list, -10, -1, 0);
// testRange(list,   1,  1, 1);
// testRange(list,   1,  2, 2);
// testRange(list,   2,  2, 1);
// testRange(list,   2,  4, 1);
// testRange(list,   2,  5, 2);
// testRange(list,  16, 16, 1);
// testRange(list,  16, 20, 1);
// testRange(list,  17, 17, 0);

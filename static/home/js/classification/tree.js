// Classification-tree D3 dendrogram — used by the Detail pages that inline the tree content
// (a GPCR class/modality/chemotype selection, not the Newick-based phylogenetic tree used
// elsewhere in the app). Reads its data from window.CLASSIFICATION_TREE_DATA and
// window.ClassificationCore for shared helpers.

// #################
// ### CUSTOM TREE MAPPER FOR CLASSIFICATION TREE ###
// #################
// This is a custom copy that can be modified without affecting datamapper.js

// Generic cleanup / normalization for display.
// - Class nodes: drop the "(...)" suffix
// - Receptor family nodes (second-to-last internal level): strip " receptors" etc.
function custom_update_tree_data(root) {
    function isLeaf(n) { return !n.children || n.children.length === 0; }

    function walk(node, depth, parent) {
        if (!node) return;

        // Default stroke color if none supplied by backend
        if (!node.color) node.color = "#333";

        // If this is the class node (child of root), shorten the label
        if (parent && parent.name === "" && node.name) {
            node.name = node.name.split(" (")[0];
            // For plots where depth-1 nodes are class symbols (A, B1, ...), show "Class X"
            // to match the desired display naming.
            const key = String(node.name).trim();
            if (CLASS_COLORS[key] && !/^Class\s+/i.test(key)) {
                node.name = "Class " + key;
            }
        }

        // If this is a receptor-family node (children are all leaves), shorten label a bit
        if (node.children && node.children.length > 0 && node.children.every(isLeaf) && node.name) {
            node.name = node.name.replace(/( receptors|neuropeptide )/g, '');
            node.name = node.name.split(" (")[0];
        }

        if (node.children) {
            node.children.forEach(child => walk(child, depth + 1, node));
        }
    }

    walk(root, 0, null);
    return root;
}
function custom_get_forced_chemotype_from_meta(stacked_meta) {
    if (!stacked_meta || !stacked_meta.length) return null;
    const hit = stacked_meta.find(x => (x.label || '').toLowerCase() === 'chemotype');
    return hit ? hit.value : null;
}

function custom_has_meta(stacked_meta, labelLower) {
    if (!stacked_meta || !stacked_meta.length) return false;
    return stacked_meta.some(x => (x.label || '').toLowerCase() === labelLower);
}

function custom_apply_tree_colors(root, stacked_meta) {
    // If Chemotype got collapsed into stacked meta, apply a single chemotype color to everything.
    const forcedChem = custom_get_forced_chemotype_from_meta(stacked_meta);
    if (forcedChem) {
        const c = custom_get_chemotype_color(forcedChem);
        (function walk(n) {
            if (!n) return;
            n.color = c;
            if (n.children) n.children.forEach(walk);
        })(root);
        return root;
    }

    // Otherwise, depth 1 nodes are chemotypes after lifting class layer.
    if (root && root.children && root.children.length) {
        root.children.forEach(function (chemNode) {
            const c = custom_get_chemotype_color(chemNode.name);
            (function walk(n) {
                if (!n) return;
                n.color = c;
                if (n.children) n.children.forEach(walk);
            })(chemNode);
        });
    }
    return root;
}

function custom_get_max_depth(node, depth) {
    depth = depth || 0;
    if (!node || !node.children || node.children.length === 0) return depth;
    return Math.max.apply(null, node.children.map(ch => custom_get_max_depth(ch, depth + 1)));
}

// Collapse single-child intermediate levels and return labels to show as stacked meta.
// IMPORTANT: This runs *after* we lift the Class layer into the center badge,
// so the expected structure here is:
// root('') → Chemotype → Family → UniProt (leaves)
function custom_collapse_singletons_after_lift(root) {
    const stacked = [];
    if (!root || !root.children || root.children.length === 0) return { data: root, stacked };

    // Collapse Chemotype if only one (common for e.g. T2 / Retinal)
    if (root.children.length === 1 && root.children[0].children) {
        const chemNode = root.children[0];
        stacked.push({ label: "Chemotype", value: chemNode.name });
        root.children = chemNode.children || [];
    }

    // Collapse Family if only one (after possible chemotype collapse)
    if (root.children.length === 1 && root.children[0].children) {
        const famNode = root.children[0];
        stacked.push({ label: "Family", value: famNode.name });
        root.children = famNode.children || [];
    }

    return { data: root, stacked };
}

// Compute branch_length strings for spacing (longest label per depth).
function custom_compute_branch_lengths(root, maxDepth) {
    const longest = {};
    function walk(node, depth) {
        if (!node) return;
        if (depth > 0 && depth < maxDepth && node.name) {
            const cur = longest[depth] || "";
            if (String(node.name).length > cur.length) longest[depth] = String(node.name);
        }
        if (node.children) node.children.forEach(ch => walk(ch, depth + 1));
    }
    walk(root, 0);
    const out = {};
    for (let d = 1; d < maxDepth; d++) out[d] = longest[d] || "";
    out[maxDepth] = "";
    return out;
}

// Lift the single Class layer into a center label and return a root whose children
// are the former class children (Chemotypes).
function custom_lift_class_layer(root) {
    let centerLabel = "";
    if (root && root.name === "" && root.children && root.children.length === 1) {
        const classNode = root.children[0];
        centerLabel = (classNode.name || "").split(" (")[0];
        root.children = classNode.children || [];
    }
    return { data: root, centerLabel };
}

// Custom draw_tree function - can be modified here
function custom_draw_tree(data, options, stacked_meta, centerLabel) {
    // Remove existing SVG if present
    d3.select('#' + options.anchor + "_svg").remove();

    var branches = {};
    var branch_offset = 0;
    var thickness = options.depth + 1;

    // Scale factor to push intermediate branches closer to outer edge
    // Higher values push branches further out (closer to leaves)
    var intermediate_branch_scale = 1.7; // Scale intermediate branches outward

    // Calculate branch offsets (sorted numerically)
    const depthKeys = Object.keys(options.branch_length || {}).map(k => parseInt(k, 10)).filter(n => !isNaN(n)).sort((a,b) => a-b);
    for (var i = 0; i < depthKeys.length; i++) {
        var k = depthKeys[i];
        if (k === options.depth) { continue; }

        // First ring should start right after the center badge
        if (k === 1) {
            branches[k] = (options.centerBadgeR || 32) + (options.centerBadgePadding || 18) + (options.firstRingExtra || 0);
            branch_offset = branches[k];
            continue;
        }

        var base_offset = 0;
        if (options.label_free.includes(k)) {
            base_offset = 10;
        } else {
            if (options.branch_trunc !== 0) {
                base_offset = 2 * options.branch_trunc + 10;
            } else {
                base_offset = custom_string_pixlen(options.branch_length[k] || "", k, options);
            }
        }
        // Apply scaling to intermediate branches to push them outward
        branch_offset = branch_offset + (base_offset * intermediate_branch_scale);
        branches[k] = branch_offset;
    }
    // Increase leaf_offset to give more space for labels (prevent overlap)
    var adjusted_leaf_offset = options.leaf_offset; // Increase spacing for leaves
    branches[options.depth] = branch_offset + adjusted_leaf_offset;

    // Depth-aware fixups:
    // When options.depth === 1, the loop above skips k===depth so we never set branches[1],
    // which previously put leaves too close to the center badge. Ensure a sane leaf radius.
    if (options.depth === 1) {
        const base = (options.centerBadgeR || 32) + (options.centerBadgePadding || 18) + (options.firstRingExtra || 0);
        branches[1] = base + adjusted_leaf_offset;
    }

    // Also ensure we have a non-zero branch_offset for diameter calculations in low-depth trees.
    if (!branch_offset || branch_offset < 1) {
        branch_offset = branches[options.depth] || adjusted_leaf_offset || 30;
    }

    // Normalize plot radius across tabs/classes.
    // Problem: if SVG intrinsic size varies, CSS scales it differently per tab,
    // making fonts look bigger/smaller. We fix that by scaling branch radii so the
    // leaf radius is consistent (and thus diameter is consistent) across plots.
    var diameterPadForTarget = (options.diameterPad !== undefined) ? options.diameterPad : 40;
    var targetLeafRadius = options.targetLeafRadius;
    if (options.targetSvgSize !== undefined) {
        // leafRadius = (svgSize - diameterPad)/2
        targetLeafRadius = (options.targetSvgSize - diameterPadForTarget) / 2;
    }
    if (targetLeafRadius !== undefined && branches[options.depth] > 0) {
        var m = targetLeafRadius / branches[options.depth];
        if (isFinite(m) && m > 0) {
            Object.keys(branches).forEach(function (k) {
                branches[k] = branches[k] * m;
            });
            branch_offset = branch_offset * m;
        }
    }

    // Per-plot radius scaling (shorten/lengthen all branches without changing angles).
    // This affects both geometry (d.y radii) and the SVG diameter/height.
    var radiusScale = (options.radiusScale !== undefined && options.radiusScale !== null) ? Number(options.radiusScale) : 1.0;
    if (isFinite(radiusScale) && radiusScale > 0) {
        Object.keys(branches).forEach(function (k) {
            branches[k] = branches[k] * radiusScale;
        });
        branch_offset = branch_offset * radiusScale;
    }

    // Calculate diameter based on content.
    // We keep this configurable because different label sizes/classes benefit from different padding.
    var diameterPad = (options.diameterPad !== undefined) ? options.diameterPad : 40; // px total padding (split across both sides)
    var extraPadding = (options.extraPadding !== undefined) ? options.extraPadding : 0; // px extra padding on top of diameterPad
    var diameter = 2 * branches[options.depth] + diameterPad;

    // Content dimensions / viewBox size
    var contentWidth = diameter + extraPadding;
    var contentHeight = diameter + extraPadding;
    var centerX = contentWidth / 2;
    var topSafetyPad = TREE_UI.leafLabelType === "Protein" ? 10 : 0;
    var centerY = (contentHeight / 2) + topSafetyPad;
    contentHeight += topSafetyPad;

    // Calculate scale factor to make tree fill the viewBox
    var scaleFactor = contentWidth / diameter;

    var tree = d3.layout.tree()
        .size([360, diameter / 2])
        .separation(function (a, b) {
            // Controls angular spacing between nodes.
            // Bootstrap: reduce the extra spacing between different receptor families.
            // Defaults: sibling=1, cousin=1 => no extra gap between families.
            var sib = (options.separationSibling !== undefined) ? options.separationSibling : 1;
            var cous = (options.separationCousin !== undefined) ? options.separationCousin : 1;
            var base = (a.parent === b.parent) ? sib : cous;
            return base / a.depth;
        });

    var diagonal = d3.svg.diagonal.radial()
        .projection(function (d) { return [d.y, d.x / 180 * Math.PI]; });

    var svg = d3.select('#' + options.anchor).append("svg")
        .attr("width", contentWidth)
        .attr("height", contentHeight)
        .attr("id", options.anchor + "_svg")
        .attr("xmlns", "http://www.w3.org/2000/svg");

    var svg_g = svg.append("g")
        .attr("transform", "translate(" + centerX + "," + centerY + ") scale(" + scaleFactor + ")");

    var nodes = tree.nodes(data);

    nodes.forEach(function (d) {
        if (d.depth === 0) {
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
        .attr("d", diagonal)
        .style("stroke", function (d) { return d.target.color; })
        .style("stroke-width", function (d) { if (d.target.depth > 0) { return thickness - d.target.depth; } else { return 0; } })
        .style("fill-opacity", 0)
        .style("opacity", 1);

    var node = svg_g.selectAll(".node")
        .data(nodes)
        .enter().append("g")
        .attr("class", "node")
        .attr("transform", function (d) { if (d.name === '') { return "rotate(" + (d.x) + ")translate(" + d.y + ")"; } else { return "rotate(" + (d.x - 90) + ")translate(" + d.y + ")"; } })

    node.filter(function (d) { return (d.depth === options.depth) })
        .attr("id", function (d) { if (d.name === '') { return "innerNode" } else { return 'X' + String(d._leafKey || d.name).toUpperCase() } });

    node.append("text")
        .attr("dy", ".31em")
        .attr("name", function (d) { if (d.name === '') { return "branch" } else { return d.name } })
        .attr("text-anchor", function (d) {
            if (d.depth === options.depth) {
                return d.x < 181 ? "start" : "end";
            } else {
                return d.x < 180 ? "end" : "start";
            }
        })
        .attr("transform", function (d) {
            // Label offsets (no dots/circles on this page)
            var labelOffset = d.depth === options.depth ? 10 : 6;
            if (d.depth === options.depth) {
                return d.x < 181 ? `translate(${labelOffset})` : `rotate(180)translate(-${labelOffset})`;
            } else {
                // Internal label offset. Push selected internal labels outward (text + pill together).
                var innerOffset = 12;
                // Depth 1 is Chemotype only for class-lifted plots (layer1Type === "chemotype").
                // If there is no Family level (depth=2 tree), Chemotype is also the last internal ring,
                // so we use an if/else to avoid applying both adjustments.
                if (d.depth === 1 && options.layer1Type === "chemotype" && !options.chemotypeCollapsed) {
                    innerOffset = innerOffset - (options.chemotypeLabelExtra || 0);
                } else if (d.depth === (options.depth - 1)) {
                    innerOffset = innerOffset - (options.familyLabelExtra || 0);
                }
                // innerOffset can be negative; avoid invalid "translate(--8)" by
                // computing the final signed translate value as a number.
                var tx = (d.x < 180) ? -innerOffset : innerOffset;
                tx = Number(tx);
                return d.x < 180 ? `translate(${tx})` : `rotate(180)translate(${tx})`;
            }
        })
        .text(function (d) {
            if (d.depth === options.depth) {
                return TREE_UI.leafLabelType === "UniProt" ? d.name.toUpperCase() : d.name;
            } else if (options.label_free.includes(d.depth)) {
                return "";
            } else if (d.depth > 0) {
                return d.name;
            } else {
                return "";
            }
        })
        .each(function (d) {
            // Leaves and some internal labels may contain entities/tags (e.g. &kappa;, GABA<sub>B</sub>).
            if (d.depth > 0 && d.name && (d.depth === options.depth || /[<&]/.test(String(d.name)))) {
                this.innerHTML = custom_formatTextWithHTML(d.name);
            }
            custom_adjust_leaf_label_baseline(this, d, options);
        })
        .call(custom_wrap, options.branch_trunc)
        .style("font-size", function (d) {
            // Custom font size logic - can be modified here
            if (options.depth === 4) {
                if (d.depth === 1) { return options.fontSize.class; }
                else if (d.depth === 2) { return options.fontSize.ligandtype; }
                else if (d.depth === 3) { return options.fontSize.receptorfamily; }
                else { return options.fontSize.receptor; }
            } else {
                if (d.depth === 1) { return options.fontSize.ligandtype; }
                else if (d.depth === 2) { return options.fontSize.receptorfamily; }
                else { return options.fontSize.receptor; }
            }
        })
        .style("font-family", function(d) {
            // Custom font family - can be modified here
            return options.fontFamily || "'Palatino Linotype', Georgia, 'Times New Roman', serif";
        })
        .style("fill", function (d) {
            if (d.color) { return "#111"; }
            else { return "#222"; };
        }).call(custom_getBB);

    // Background shapes for internal labels: tight rounded boxes around the text.
    // IMPORTANT: internal labels use a text-specific transform (sometimes rotate(180)).
    // To keep the box behind the text on both sides, copy the text's transform.
    function labelBoxPadding(depth, maxDepth) {
        const padX = (options.labelBoxPadX !== undefined) ? options.labelBoxPadX : 7;
        const padY = (options.labelBoxPadY !== undefined) ? options.labelBoxPadY : 2;
        return { x: padX, y: padY };
    }

    node.filter(function (d) { return (d.depth !== options.depth && d.depth > 0); })
        .each(function (d) {
            const g = d3.select(this);
            const t = g.select("text");
            if (t.empty()) return;

            const bb = t.node().getBBox(); // local bbox (ignores transform)
            const pad = labelBoxPadding(d.depth, options.depth);
            const tr = t.attr("transform"); // copy transform so box aligns with text
            const extraTop = /[A-Z]/.test(String(d.name || "")) ? 1 : 0;

            g.insert("rect", "text")
                .attr("transform", tr || null)
                .attr("x", bb.x - pad.x)
                // Some fonts render slightly "high" relative to the bbox; allow tiny manual nudge.
                .attr("y", bb.y - pad.y + (options.labelBoxDy || 0) + (options.labelBoxTopTrim || 0) - extraTop)
                .attr("width", bb.width + pad.x * 2)
                .attr("height", bb.height + pad.y * 2 - (options.labelBoxTopTrim || 0) - (options.labelBoxBottomTrim || 0) + extraTop)
                .attr("rx", 10)
                .attr("ry", 10)
                .style("fill", "#FFF")
                .style("stroke", ((options.labelBoxStrokeColorMode || "chemotype") === "black") ? "#000" : (d.color || "#999"))
                .style("stroke-width", "1px");
        });

    // Optional: leaf pills (outermost labels)
    // These are separate from internal label pills and can be tuned independently.
    if (options.leafPills) {
        const leafPadX = (options.leafBoxPadX !== undefined) ? options.leafBoxPadX : 4;
        const leafPadY = (options.leafBoxPadY !== undefined) ? options.leafBoxPadY : 2;
        const leafDy = (options.leafBoxDy !== undefined) ? options.leafBoxDy : 0;
        // getBBox() height is driven by the font's full ascent+descent metric, not the actual
        // ink height of these (mostly ascender-free, uppercase/digit) labels, so most of a leaf
        // pill's "extra" height sits above the glyphs. Trim it off the top only, keeping the
        // bottom edge (and its padding) exactly where it was.
        const leafTopTrim = (options.leafBoxTopTrim !== undefined) ? options.leafBoxTopTrim : 0;
        const leafBottomTrim = (options.leafBoxBottomTrim !== undefined) ? options.leafBoxBottomTrim : 0;
        const leafRx = (options.leafBoxRx !== undefined) ? options.leafBoxRx : 4;
        const leafStrokeW = (options.leafBoxStrokeWidth !== undefined) ? options.leafBoxStrokeWidth : 1;
        const leafStrokeMode = (options.leafBoxStrokeColorMode !== undefined) ? options.leafBoxStrokeColorMode : "chemotype"; // "chemotype" | "black"

        node.filter(function (d) { return (d.depth === options.depth); })
            .each(function (d) {
                const g = d3.select(this);
                const t = g.select("text");
                if (t.empty()) return;

                const bb = t.node().getBBox();
                const tr = t.attr("transform"); // copy transform so pill aligns with leaf text
                const extraTop = /[A-Z]/.test(String(d.name || "")) ? 1 : 0;

                g.insert("rect", "text")
                    .attr("transform", tr || null)
                    .attr("x", bb.x - leafPadX)
                    .attr("y", bb.y - leafPadY + leafDy + leafTopTrim - extraTop)
                    .attr("width", bb.width + leafPadX * 2)
                    .attr("height", bb.height + leafPadY * 2 - leafTopTrim - leafBottomTrim + extraTop)
                    .attr("rx", leafRx)
                    .attr("ry", leafRx)
                    .style("fill", "#fff")
                    .style("stroke", (leafStrokeMode === "black") ? "#000" : (d.color || "#999"))
                    .style("stroke-width", leafStrokeW + "px");
            });
    }

    // Center badge removed (no middle label/pill).

    function applyCircularBoundingBoxPadding() {
        var pad = (options.circularBBoxPadding !== undefined) ? Number(options.circularBBoxPadding) : 14;
        var bb;
        var renderedLeft;
        var renderedRight;
        var renderedTop;
        var renderedBottom;
        var addLeft;
        var addRight;
        var addTop;
        var addBottom;
        if (!isFinite(pad) || pad < 0 || !svg_g.node()) {
            return;
        }
        try {
            bb = svg_g.node().getBBox();
        } catch (err) {
            return;
        }
        if (!bb || !isFinite(bb.x) || !isFinite(bb.y) || !isFinite(bb.width) || !isFinite(bb.height)) {
            return;
        }

        renderedLeft = centerX + (bb.x * scaleFactor);
        renderedRight = centerX + ((bb.x + bb.width) * scaleFactor);
        renderedTop = centerY + (bb.y * scaleFactor);
        renderedBottom = centerY + ((bb.y + bb.height) * scaleFactor);

        addLeft = Math.max(0, pad - renderedLeft);
        addRight = Math.max(0, renderedRight + pad - contentWidth);
        addTop = Math.max(0, pad - renderedTop);
        addBottom = Math.max(0, renderedBottom + pad - contentHeight);

        if (!(addLeft || addRight || addTop || addBottom)) {
            return;
        }

        contentWidth += addLeft + addRight;
        contentHeight += addTop + addBottom;
        centerX += addLeft;
        centerY += addTop;
        svg.attr("width", contentWidth)
            .attr("height", contentHeight);
        svg_g.attr("transform", "translate(" + centerX + "," + centerY + ") scale(" + scaleFactor + ")");
    }

    applyCircularBoundingBoxPadding();

    function custom_string_pixlen(text, depth, options) {
        var canvas = document.createElement('canvas');
        var ctx = canvas.getContext("2d");
        var fontFamily = options.fontFamily || "'Palatino Linotype', Georgia, 'Times New Roman', serif";

        if (options.depth === 4) {
            if (depth === 1) {
                ctx.font = options.fontSize.class + " " + fontFamily;
            } else if (depth === 2) {
                ctx.font = options.fontSize.ligandtype + " " + fontFamily;
            } else if (depth === 3) {
                ctx.font = options.fontSize.receptorfamily + " " + fontFamily;
            } else {
                ctx.font = options.fontSize.receptor + " " + fontFamily;
            }
        } else if (options.depth === 3) {
            if (depth === 1) {
                ctx.font = options.fontSize.ligandtype + " " + fontFamily;
            } else if (depth === 2) {
                ctx.font = options.fontSize.receptorfamily + " " + fontFamily;
            } else {
                ctx.font = options.fontSize.receptor + " " + fontFamily;
            }
        }

        return parseInt(ctx.measureText(text).width,10);
    }

    // Set viewBox to scale to container width, align top, allow bottom to extend
    // contentWidth and contentHeight already calculated above
    svg.attr('viewBox', '0 0 ' + contentWidth + ' ' + contentHeight)
        .attr('preserveAspectRatio', 'xMidYMin meet');

    // Don't change the transform - svg_g is already correctly centered at (diameter/2, diameter/2)
    // The tree is calculated to extend from center to radius = diameter/2, so it fills ~95% of the SVG
    // The viewBox will allow it to scale down to fit the container while maintaining aspect ratio
}

// Global helper functions for all drawing functions
function custom_getBB(selection) {
    selection.each(function (d) { d.bbox = this.getBBox(); });
}

function custom_adjust_leaf_label_baseline(textNode, datum, options) {
    if (!textNode || !datum || datum.depth !== options.depth) return;
    const label = String(datum.name || "");
    if (!/[a-z]/.test(label) || /<\s*sub\b/i.test(label)) return;
    d3.select(textNode).attr("dy", ".16em");
}

function custom_wrap(text, width) {
    if (width === 0) {
        return;
    }
    text.each(function () {
        var text = d3.select(this),
            words = text.text().split(/\s+/).reverse(),
            word,
            line = [],
            lineNumber = 0,
            lineHeight = 1.1,
            y = text.attr("y"),
            dy = parseFloat(text.attr("dy")),
            tspan = text.text(null).append("tspan").attr("x", 0).attr("y", y).attr("dy", dy + "em");

        word = words.pop();
        while (word !== undefined) {
            line.push(word);
            tspan.text(line.join(" "));

            if (tspan.node().getComputedTextLength() > width) {
                line.pop();
                tspan.text(line.join(" "));
                line = [word];
                tspan = text.append("tspan")
                            .attr("x", 0)
                            .attr("y", y)
                            .attr("dy", ++lineNumber * lineHeight + dy + "em")
                            .text(word);
            }

            word = words.pop();
        }
    });
}

// Custom draw_dendrogram_curved function - horizontal curved dendrogram
function custom_draw_dendrogram_curved(data, options, stacked_meta, centerLabel) {
    // Remove existing SVG if present
    d3.select('#' + options.anchor + "_svg").remove();

    // Set defaults for options
    if (!options.label_free) options.label_free = [];
    if (options.branch_trunc === undefined) options.branch_trunc = 0;

    // Calculate dimensions - use similar width to tree (targetSvgSize)
    var targetSvgSize = options.targetSvgSize || 800;
    var margin = { top: 20, right: 10, bottom: 40, left: 10 };
    var baseWidth = targetSvgSize;
    var baseHeight = targetSvgSize;

    // Calculate number of leaves for height scaling
    function countLeaves(node) {
        if (!node.children || node.children.length === 0) return 1;
        return node.children.reduce(function(sum, child) { return sum + countLeaves(child); }, 0);
    }
    var leafCount = countLeaves(data);
    var nodeHeight = Math.max(18, Math.min(25, (baseHeight - margin.top - margin.bottom) / Math.max(leafCount, 1)));
    var height = Math.max(leafCount * nodeHeight + margin.top + margin.bottom, 400);

    // Calculate horizontal spacing based on depth
    // Add compression factor to bring root closer to leaves
    var horizontalCompression = (options.dendrogramHorizontalCompression !== undefined) ? options.dendrogramHorizontalCompression : 1.0;
    var depthKeys = Object.keys(options.branch_length || {}).map(k => parseInt(k, 10)).filter(n => !isNaN(n)).sort((a,b) => a-b);
    var maxDepth = options.depth || 3;
    var levelWidth = ((baseWidth - margin.left - margin.right) / (maxDepth + 1)) * horizontalCompression;
    var width = (maxDepth + 1) * levelWidth + margin.left + margin.right;

    var thickness = options.depth + 1;

    var tree = d3.layout.tree()
        .size([height - margin.top - margin.bottom, width - margin.left - margin.right])
        .separation(function (a, b) {
            var sib = (options.separationSibling !== undefined) ? options.separationSibling : 1;
            var cous = (options.separationCousin !== undefined) ? options.separationCousin : 1;
            var base = (a.parent === b.parent) ? sib : cous;
            return base;
        });

    var diagonal = d3.svg.diagonal()
        .projection(function (d) { return [d.y, d.x]; });

    // Initial width - will be recalculated after node positioning
    var initialWidth = width;
    var svg = d3.select('#' + options.anchor).append("svg")
        .attr("width", initialWidth)
        .attr("height", height)
        .attr("id", options.anchor + "_svg")
        .attr("xmlns", "http://www.w3.org/2000/svg");

    var svg_g = svg.append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    var nodes = tree.nodes(data);
    var links = tree.links(nodes);

    // Adjust node positions for horizontal layout
    // Use nested dict-based spacing: levelSpacing[type][maxDepth][level] = spacing multiplier
    // Example: levelSpacing["Class"][3][1] = 0.8 means for Class type, maxDepth=3 trees, level 1 uses 80% spacing
    var levelSpacing = options.dendrogramLevelSpacing || {};
    var plotType = options.plotType || "Class"; // Get the type (Class, Modality, Chemotype)

    // Helper function to get spacing for a specific level at a given maxDepth
    function getLevelSpacing(level) {
        // Check if we have type-specific spacing defined
        if (levelSpacing[plotType] && levelSpacing[plotType][maxDepth] && levelSpacing[plotType][maxDepth][level] !== undefined) {
            return levelSpacing[plotType][maxDepth][level];
        }
        // Fallback: check for type-agnostic spacing (backwards compatibility)
        if (levelSpacing[maxDepth] && levelSpacing[maxDepth][level] !== undefined) {
            return levelSpacing[maxDepth][level];
        }
        // Default: normal spacing
        return 1.0;
    }

    // Calculate X positions by depth level (all nodes at same depth get same X)
    // First, calculate cumulative X positions for each depth level
    var depthPositions = {};
    var cumulativeX = 0;
    for (var depth = 1; depth <= maxDepth; depth++) {
        var spacing = getLevelSpacing(depth);
        cumulativeX += levelWidth * spacing;
        depthPositions[depth] = cumulativeX;
    }

    // Now assign X positions to nodes based on their depth
    nodes.forEach(function (d) {
        if (d.depth === 0) {
            d.y = 0;
        } else {
            d.y = depthPositions[d.depth];
        }
    });

    // Recalculate width based on actual leaf positions
    var maxY = 0;
    nodes.forEach(function(d) {
        if (d.y > maxY) maxY = d.y;
    });
    var newWidth = maxY + margin.left + margin.right + 50; // Add some padding for labels
    if (newWidth !== width) {
        width = newWidth;
        svg.attr("width", width);
    }

    var link = svg_g.append("g")
        .attr("class", "links")
        .selectAll("path")
        .data(links)
        .enter().append("path")
        .each(function (d) { d.target.linkNode = this; })
        .attr("d", diagonal)
        .style("stroke", function (d) { return d.target.color; })
        .style("stroke-width", function (d) { if (d.target.depth > 0) { return thickness - d.target.depth; } else { return 0; } })
        .style("fill-opacity", 0)
        .style("opacity", 1);

    var node = svg_g.selectAll(".node")
        .data(nodes)
        .enter().append("g")
        .attr("class", "node")
        .attr("transform", function (d) { return "translate(" + d.y + "," + d.x + ")"; });

    node.filter(function (d) { return (d.depth === options.depth) })
        .attr("id", function (d) { if (d.name === '') { return "innerNode" } else { return 'X' + String(d._leafKey || d.name).toUpperCase() } });

    node.append("text")
        .attr("dy", ".31em")
        .attr("name", function (d) { if (d.name === '') { return "branch" } else { return d.name } })
        .attr("text-anchor", function (d) {
            if (d.depth === options.depth) {
                return "start";
            } else {
                return "end";
            }
        })
        .attr("transform", function (d) {
            var labelOffset = d.depth === options.depth ? 10 : -6;
            // Add horizontal offsets for chemotype and receptor family labels
            if (d.depth === 1 && options.dendrogramChemotypeOffset !== undefined) {
                labelOffset += options.dendrogramChemotypeOffset;
            } else if (d.depth === 2 && options.dendrogramFamilyOffset !== undefined) {
                labelOffset += options.dendrogramFamilyOffset;
            }
            return "translate(" + labelOffset + ",0)";
        })
        .text(function (d) {
            if (d.depth === options.depth) {
                return TREE_UI.leafLabelType === "UniProt" ? d.name.toUpperCase() : d.name;
            } else if (options.label_free && options.label_free.includes(d.depth)) {
                return "";
            } else if (d.depth > 0) {
                return d.name;
            } else {
                return "";
            }
        })
        .each(function (d) {
            if (d.depth > 0 && d.name && (d.depth === options.depth || /[<&]/.test(String(d.name)))) {
                this.innerHTML = custom_formatTextWithHTML(d.name);
            }
            custom_adjust_leaf_label_baseline(this, d, options);
        })
        .call(custom_wrap, options.branch_trunc)
        .style("font-size", function (d) {
            // Use same font sizes as tree for consistency
            if (options.depth === 4) {
                if (d.depth === 1) { return options.fontSize.class; }
                else if (d.depth === 2) { return options.fontSize.ligandtype; }
                else if (d.depth === 3) { return options.fontSize.receptorfamily; }
                else { return options.fontSize.receptor; }
            } else {
                if (d.depth === 1) { return options.fontSize.ligandtype; }
                else if (d.depth === 2) { return options.fontSize.receptorfamily; }
                else { return options.fontSize.receptor; }
            }
        })
        .style("font-family", function(d) {
            return options.fontFamily || "'Palatino Linotype', Georgia, 'Times New Roman', serif";
        })
        .style("fill", function (d) {
            if (d.color) { return "#111"; }
            else { return "#222"; };
        }).call(custom_getBB);

    // Background shapes for internal labels
    function labelBoxPadding(depth, maxDepth) {
        const padX = (options.labelBoxPadX !== undefined) ? options.labelBoxPadX : 7;
        const padY = (options.labelBoxPadY !== undefined) ? options.labelBoxPadY : 2;
        return { x: padX, y: padY };
    }

    node.filter(function (d) { return (d.depth !== options.depth && d.depth > 0); })
        .each(function (d) {
            const g = d3.select(this);
            const t = g.select("text");
            if (t.empty()) return;

            // Get bbox in local coordinates (before transform)
            const bb = t.node().getBBox();
            const pad = labelBoxPadding(d.depth, options.depth);
            const tr = t.attr("transform");
            const extraTop = /[A-Z]/.test(String(d.name || "")) ? 1 : 0;

            g.insert("rect", "text")
                .attr("transform", tr || null)
                .attr("x", bb.x - pad.x)
                .attr("y", bb.y - pad.y + (options.labelBoxDy || 0) + (options.labelBoxTopTrim || 0) - extraTop)
                .attr("width", bb.width + pad.x * 2)
                .attr("height", bb.height + pad.y * 2 - (options.labelBoxTopTrim || 0) - (options.labelBoxBottomTrim || 0) + extraTop)
                .attr("rx", 10)
                .attr("ry", 10)
                .style("fill", "#FFF")
                .style("stroke", ((options.labelBoxStrokeColorMode || "chemotype") === "black") ? "#000" : (d.color || "#999"))
                .style("stroke-width", "1px");
        });

    // Optional: leaf pills
    if (options.leafPills) {
        const leafPadX = (options.leafBoxPadX !== undefined) ? options.leafBoxPadX : 4;
        const leafPadY = (options.leafBoxPadY !== undefined) ? options.leafBoxPadY : 2;
        const leafDy = (options.leafBoxDy !== undefined) ? options.leafBoxDy : 0;
        // getBBox() height is driven by the font's full ascent+descent metric, not the actual
        // ink height of these (mostly ascender-free, uppercase/digit) labels, so most of a leaf
        // pill's "extra" height sits above the glyphs. Trim it off the top only, keeping the
        // bottom edge (and its padding) exactly where it was.
        const leafTopTrim = (options.leafBoxTopTrim !== undefined) ? options.leafBoxTopTrim : 0;
        const leafBottomTrim = (options.leafBoxBottomTrim !== undefined) ? options.leafBoxBottomTrim : 0;
        const leafRx = (options.leafBoxRx !== undefined) ? options.leafBoxRx : 4;
        const leafStrokeW = (options.leafBoxStrokeWidth !== undefined) ? options.leafBoxStrokeWidth : 1;
        const leafStrokeMode = (options.leafBoxStrokeColorMode !== undefined) ? options.leafBoxStrokeColorMode : "chemotype";

        node.filter(function (d) { return (d.depth === options.depth); })
            .each(function (d) {
                const g = d3.select(this);
                const t = g.select("text");
                if (t.empty()) return;

                const bb = t.node().getBBox();
                const tr = t.attr("transform");
                const extraTop = /[A-Z]/.test(String(d.name || "")) ? 1 : 0;

                g.insert("rect", "text")
                    .attr("transform", tr || null)
                    .attr("x", bb.x - leafPadX)
                    .attr("y", bb.y - leafPadY + leafDy + leafTopTrim - extraTop)
                    .attr("width", bb.width + leafPadX * 2)
                    .attr("height", bb.height + leafPadY * 2 - leafTopTrim - leafBottomTrim + extraTop)
                    .attr("rx", leafRx)
                    .attr("ry", leafRx)
                    .style("fill", "#fff")
                    .style("stroke", (leafStrokeMode === "black") ? "#000" : (d.color || "#999"))
                    .style("stroke-width", leafStrokeW + "px");
            });
    }

    svg.attr('viewBox', '0 0 ' + width + ' ' + height)
        .attr('preserveAspectRatio', 'xMidYMin meet');
}

// Custom draw_dendrogram_straight function - horizontal straight dendrogram with 90-degree angles
function custom_draw_dendrogram_straight(data, options, stacked_meta, centerLabel) {
    // Remove existing SVG if present
    d3.select('#' + options.anchor + "_svg").remove();

    // Set defaults for options
    if (!options.label_free) options.label_free = [];
    if (options.branch_trunc === undefined) options.branch_trunc = 0;

    // Calculate dimensions - use similar width to tree (targetSvgSize)
    var targetSvgSize = options.targetSvgSize || 800;
    var margin = { top: 20, right: 250, bottom: 40, left: 100 };
    var baseWidth = targetSvgSize;
    var baseHeight = targetSvgSize;

    // Calculate number of leaves for height scaling
    function countLeaves(node) {
        if (!node.children || node.children.length === 0) return 1;
        return node.children.reduce(function(sum, child) { return sum + countLeaves(child); }, 0);
    }
    var leafCount = countLeaves(data);
    var nodeHeight = Math.max(18, Math.min(25, (baseHeight - margin.top - margin.bottom) / Math.max(leafCount, 1)));
    var height = Math.max(leafCount * nodeHeight + margin.top + margin.bottom, 400);

    // Calculate horizontal spacing based on depth
    // Add compression factor to bring root closer to leaves
    var horizontalCompression = (options.dendrogramHorizontalCompression !== undefined) ? options.dendrogramHorizontalCompression : 1.0;
    var depthKeys = Object.keys(options.branch_length || {}).map(k => parseInt(k, 10)).filter(n => !isNaN(n)).sort((a,b) => a-b);
    var maxDepth = options.depth || 3;
    var levelWidth = ((baseWidth - margin.left - margin.right) / (maxDepth + 1)) * horizontalCompression;
    var width = (maxDepth + 1) * levelWidth + margin.left + margin.right;

    var thickness = options.depth + 1;

    var tree = d3.layout.tree()
        .size([height - margin.top - margin.bottom, width - margin.left - margin.right])
        .separation(function (a, b) {
            var sib = (options.separationSibling !== undefined) ? options.separationSibling : 1;
            var cous = (options.separationCousin !== undefined) ? options.separationCousin : 1;
            var base = (a.parent === b.parent) ? sib : cous;
            return base;
        });

    // Orthogonal path generator for 90-degree angles
    function orthogonal(d) {
        return "M" + d.source.y + "," + d.source.x
             + "H" + d.target.y
             + "V" + d.target.x;
    }

    var svg = d3.select('#' + options.anchor).append("svg")
        .attr("width", width)
        .attr("height", height)
        .attr("id", options.anchor + "_svg")
        .attr("xmlns", "http://www.w3.org/2000/svg");

    var svg_g = svg.append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    var nodes = tree.nodes(data);
    var links = tree.links(nodes);

    // Use level spacing configuration for straight dendrogram if available
    var levelSpacing = options.dendrogramStraightLevelSpacing || {};
    var plotType = options.plotType || "Class";

    // Helper function to get spacing for a specific level at a given maxDepth
    function getLevelSpacing(level) {
        if (levelSpacing[plotType] && levelSpacing[plotType][maxDepth] && levelSpacing[plotType][maxDepth][level] !== undefined) {
            return levelSpacing[plotType][maxDepth][level];
        }
        if (levelSpacing[maxDepth] && levelSpacing[maxDepth][level] !== undefined) {
            return levelSpacing[maxDepth][level];
        }
        return 1.0;
    }

    // Calculate X positions by depth level with spacing configuration
    var depthPositions = {};
    var cumulativeX = 0;
    for (var depth = 1; depth <= maxDepth; depth++) {
        var spacing = getLevelSpacing(depth);
        cumulativeX += levelWidth * spacing;
        depthPositions[depth] = cumulativeX;
    }

    // Adjust node positions for horizontal layout
    nodes.forEach(function (d) {
        if (d.depth === 0) {
            d.y = 0;
        } else {
            d.y = depthPositions[d.depth];
        }
    });

    // Build a map of siblings to create stacked appearance (prevent overlapping)
    // ONLY for root (depth 0) to level 1 connections (Class/Modality/Chemotype)
    // Do NOT apply stacking for receptor families or leaves
    var rootNode = nodes.find(function(n) { return n.depth === 0; });
    var rootX = rootNode ? rootNode.x : 0;

    var rootChildren = [];
    links.forEach(function(link) {
        var parent = link.source;
        var child = link.target;
        // Only apply stacking from root (depth 0) to level 1
        if (parent.depth === 0 && child.depth === 1) {
            rootChildren.push({parent: parent, child: child});
        }
    });

    // Calculate offsets for siblings to stack them (only for root->level1)
    var branchStackOffset = (options.branchStackOffset !== undefined) ? options.branchStackOffset : 4.0;

    if (rootChildren.length > 0) {
        // Sort siblings by their vertical position (x coordinate)
        rootChildren.sort(function(a, b) { return a.child.x - b.child.x; });

        // Split into above (x <= root.x) and below (x > root.x)
        // Use <= so ties (rare) still get a deterministic side and receive split coordinates.
        var aboveChildren = rootChildren.filter(function(item) { return item.child.x <= rootX; });
        var belowChildren = rootChildren.filter(function(item) { return item.child.x > rootX; });

        // Root split tuning:
        // - branchStackOffset: step size for offsets (e.g. 4)
        // - rootSplitBase: how far from root (in X direction) the split happens
        // - rootStackReverseAbove: reverse ordering for "above" group (to align with bottom)
        var rootSplitBase = (options.rootSplitBase !== undefined) ? options.rootSplitBase : 12.0;
        var rootStackReverseAbove = (options.rootStackReverseAbove !== undefined) ? !!options.rootStackReverseAbove : true;

        // Calculate maxOffset so both halves share the same ladder (0..maxOffset)
        var totalAbove = aboveChildren.length;
        var totalBelow = belowChildren.length;
        var maxOffset = Math.max(
            (totalAbove > 0 ? (totalAbove - 1) : 0),
            (totalBelow > 0 ? (totalBelow - 1) : 0)
        ) * branchStackOffset;

        // For above: assign offsets either normal (0,4,8,...) or reversed (max, max-4, ...)
        aboveChildren.forEach(function(item, idx) {
            item.child._stackOffset = rootStackReverseAbove ? (maxOffset - (idx * branchStackOffset)) : (idx * branchStackOffset);
            item.child._rootSplitX = (rootNode ? rootNode.y : 0) + rootSplitBase + item.child._stackOffset;
        });

        // For below: always reversed (max, max-4, ...)
        belowChildren.forEach(function(item, idx) {
            item.child._stackOffset = maxOffset - (idx * branchStackOffset);
            item.child._rootSplitX = (rootNode ? rootNode.y : 0) + rootSplitBase + item.child._stackOffset;
        });

        // Root "horizontal split" (visual de-overlap):
        // Give each root->level1 branch its own horizontal lane by offsetting it vertically
        // near the root, then merging back into the tree at the splitX.
        var rootStemYOffsetStep = (options.rootStemYOffsetStep !== undefined) ? options.rootStemYOffsetStep : 3.0;
        var rootStemReverseAbove = (options.rootStemReverseAbove !== undefined) ? !!options.rootStemReverseAbove : true;
        var aboveLen = aboveChildren.length;
        aboveChildren.forEach(function(item, idx) {
            // Above: move up (negative y in screen coords)
            // Optionally reverse the ordering so the top half aligns visually with the bottom half.
            var k = rootStemReverseAbove ? (aboveLen - 1 - idx) : idx;
            // Remove the visual "gap" between top and bottom by allowing the top group
            // to occupy the 0-lane (k=0 => offset 0). Bottom still starts at +1*step.
            item.child._rootStemYOffset = -(k) * rootStemYOffsetStep;
        });
        belowChildren.forEach(function(item, idx) {
            // Below: move down (positive y in screen coords)
            item.child._rootStemYOffset = (idx + 1) * rootStemYOffsetStep;
        });

        // Store root split information for path generation
        // The root split should extend horizontally to cover all branches
        rootNode._hasSplit = true;
        rootNode._maxOffset = maxOffset;
    }

    // Ensure all other nodes have no stacking offset
    nodes.forEach(function(node) {
        if (node.depth !== 1 && !node._stackOffset) {
            node._stackOffset = 0;
        }
    });

    // Enhanced orthogonal path generator that stacks branches and adds leaf splits
    function orthogonalStacked(d) {
        var sourceX = d.source.x;
        var sourceY = d.source.y;
        var targetX = d.target.x;
        var targetY = d.target.y;

        // For leaf nodes, add a small horizontal split before the leaf
        var leafSplitLength = (options.leafSplitLength !== undefined) ? options.leafSplitLength : 8;

        if (d.target.depth === options.depth) {
            // Leaf node: create a split branch (horizontal line before vertical to leaf)
            // Path: horizontal from source -> vertical to leaf level -> horizontal split -> vertical to leaf
            // No stacking offset for leaves
            var splitStartY = targetY - leafSplitLength;
            return "M" + sourceY + "," + sourceX
                 + "H" + splitStartY
                 + "V" + targetX
                 + "H" + targetY;
        } else if (d.source.depth === 0 && d.target.depth === 1) {
            // Root → level 1: split close to the root, using a consistent X position:
            // splitX = rootSplitBase + stackOffset (computed earlier and stored on the node).
            // This avoids "back-and-forth" artifacts from mixing coordinate bases.
            var splitX = (d.target && d.target._rootSplitX !== undefined) ? d.target._rootSplitX : targetY;
            // Never let the split go past the child's column, otherwise the path would "backtrack".
            if (splitX > targetY) splitX = targetY;
            var stemDY = (d.target && d.target._rootStemYOffset !== undefined) ? d.target._rootStemYOffset : 0;
            // Small cosmetic stub to the left of the root.
            var rootLeftStub = (options.rootLeftStub !== undefined) ? Number(options.rootLeftStub) : 6.0;
            if (!isFinite(rootLeftStub)) rootLeftStub = 0;
            rootLeftStub = Math.max(0, Math.min(rootLeftStub, 40));
            // Visual horizontal split at root:
            // (optional) vertical spread -> horizontal to splitX -> vertical to child -> horizontal to child level
            // Note: start directly on the stem lane to avoid drawing the initial connector from the root point.
            return "M" + (sourceY - rootLeftStub) + "," + (sourceX + stemDY)
                 + "H" + splitX
                 + "V" + targetX
                 + "H" + targetY;
        } else {
            // Internal node: standard orthogonal dendrogram path (no stacking).
            // Stacking is ONLY for root → level 1.
            var offsetY = targetY;
            return "M" + sourceY + "," + sourceX
                 + "H" + offsetY
                 + "V" + targetX;
        }
    }

    // Uniform stroke width for Tree - Straight (same as level 1: thickness - 1)
    var strokeWidth = Math.max(1, thickness - 2.5);
    var link = svg_g.append("g")
        .attr("class", "links")
        .selectAll("path")
        .data(links)
        .enter().append("path")
        .each(function (d) { d.target.linkNode = this; })
        .attr("d", orthogonalStacked)
        .style("stroke", function (d) { return d.target.color; })
        .style("stroke-width", function (d) { if (d.target.depth > 0) { return strokeWidth; } else { return 0; } })
        .style("fill-opacity", 0)
        .style("opacity", 1);

    var node = svg_g.selectAll(".node")
        .data(nodes)
        .enter().append("g")
        .attr("class", "node")
        .attr("transform", function (d) { return "translate(" + d.y + "," + d.x + ")"; });

    node.filter(function (d) { return (d.depth === options.depth) })
        .attr("id", function (d) { if (d.name === '') { return "innerNode" } else { return 'X' + String(d._leafKey || d.name).toUpperCase() } });

    node.append("text")
        .attr("dy", ".31em")
        .attr("name", function (d) { if (d.name === '') { return "branch" } else { return d.name } })
        .attr("text-anchor", function (d) {
            if (d.depth === options.depth) {
                return "start";
            } else {
                return (options.dendrogramCenterInternalLabels ? "middle" : "end");
            }
        })
        .attr("transform", function (d) {
            var labelOffset = d.depth === options.depth ? 10 : -6;
            // Add horizontal offsets for internal labels only (chemotype, receptor family) — not for leaves.
            // Otherwise in single-level trees leaves get the chemotype offset and sit too far right / get clipped.
            if (d.depth !== options.depth) {
                if (d.depth === 1 && options.dendrogramChemotypeOffset !== undefined) {
                    labelOffset += options.dendrogramChemotypeOffset;
                } else if (d.depth === 2 && options.dendrogramFamilyOffset !== undefined) {
                    labelOffset += options.dendrogramFamilyOffset;
                }
            }
            // Center internal labels/pills on the middle of their outgoing horizontal segment (to the right),
            // instead of sitting on the elbow.
            if (options.dendrogramCenterInternalLabels && d.depth > 0 && d.depth !== options.depth && d.children && d.children.length) {
                // All children are at the next depth, so y is constant for that segment.
                var childY = d.children[0].y;
                var midShift = (childY - d.y) / 2;
                // When centering, don't include the -6 elbow offset; keep only a small optional nudge.
                var nudge = 0;
                if (d.depth === 1 && options.dendrogramChemotypeOffset !== undefined) {
                    nudge += options.dendrogramChemotypeOffset;
                } else if (d.depth === 2 && options.dendrogramFamilyOffset !== undefined) {
                    nudge += options.dendrogramFamilyOffset;
                }
                return "translate(" + (midShift + nudge) + ",0)";
            }
            return "translate(" + labelOffset + ",0)";
        })
        .text(function (d) {
            if (d.depth === options.depth) {
                return TREE_UI.leafLabelType === "UniProt" ? d.name.toUpperCase() : d.name;
            } else if (options.label_free && options.label_free.includes(d.depth)) {
                return "";
            } else if (d.depth > 0) {
                return d.name;
            } else {
                return "";
            }
        })
        .each(function (d) {
            if (d.depth > 0 && d.name && (d.depth === options.depth || /[<&]/.test(String(d.name)))) {
                this.innerHTML = custom_formatTextWithHTML(d.name);
            }
            custom_adjust_leaf_label_baseline(this, d, options);
        })
        .call(custom_wrap, options.branch_trunc)
        .style("font-size", function (d) {
            // Use same font sizes as tree for consistency
            if (options.depth === 4) {
                if (d.depth === 1) { return options.fontSize.class; }
                else if (d.depth === 2) { return options.fontSize.ligandtype; }
                else if (d.depth === 3) { return options.fontSize.receptorfamily; }
                else { return options.fontSize.receptor; }
            } else {
                if (d.depth === 1) { return options.fontSize.ligandtype; }
                else if (d.depth === 2) { return options.fontSize.receptorfamily; }
                else { return options.fontSize.receptor; }
            }
        })
        .style("font-family", function(d) {
            return options.fontFamily || "'Palatino Linotype', Georgia, 'Times New Roman', serif";
        })
        .style("fill", function (d) {
            if (d.color) { return "#111"; }
            else { return "#222"; };
        }).call(custom_getBB);

    // Background shapes for internal labels
    function labelBoxPadding(depth, maxDepth) {
        const padX = (options.labelBoxPadX !== undefined) ? options.labelBoxPadX : 7;
        const padY = (options.labelBoxPadY !== undefined) ? options.labelBoxPadY : 2;
        return { x: padX, y: padY };
    }

    node.filter(function (d) { return (d.depth !== options.depth && d.depth > 0); })
        .each(function (d) {
            const g = d3.select(this);
            const t = g.select("text");
            if (t.empty()) return;

            // Get bbox in local coordinates (before transform)
            const bb = t.node().getBBox();
            const pad = labelBoxPadding(d.depth, options.depth);
            const tr = t.attr("transform");
            const extraTop = /[A-Z]/.test(String(d.name || "")) ? 1 : 0;

            g.insert("rect", "text")
                .attr("transform", tr || null)
                .attr("x", bb.x - pad.x)
                .attr("y", bb.y - pad.y + (options.labelBoxDy || 0) + (options.labelBoxTopTrim || 0) - extraTop)
                .attr("width", bb.width + pad.x * 2)
                .attr("height", bb.height + pad.y * 2 - (options.labelBoxTopTrim || 0) - (options.labelBoxBottomTrim || 0) + extraTop)
                .attr("rx", 10)
                .attr("ry", 10)
                .style("fill", "#FFF")
                .style("stroke", ((options.labelBoxStrokeColorMode || "chemotype") === "black") ? "#000" : (d.color || "#999"))
                .style("stroke-width", "1px");
        });

    // Optional: leaf pills
    if (options.leafPills) {
        const leafPadX = (options.leafBoxPadX !== undefined) ? options.leafBoxPadX : 4;
        const leafPadY = (options.leafBoxPadY !== undefined) ? options.leafBoxPadY : 2;
        const leafDy = (options.leafBoxDy !== undefined) ? options.leafBoxDy : 0;
        // getBBox() height is driven by the font's full ascent+descent metric, not the actual
        // ink height of these (mostly ascender-free, uppercase/digit) labels, so most of a leaf
        // pill's "extra" height sits above the glyphs. Trim it off the top only, keeping the
        // bottom edge (and its padding) exactly where it was.
        const leafTopTrim = (options.leafBoxTopTrim !== undefined) ? options.leafBoxTopTrim : 0;
        const leafBottomTrim = (options.leafBoxBottomTrim !== undefined) ? options.leafBoxBottomTrim : 0;
        const leafRx = (options.leafBoxRx !== undefined) ? options.leafBoxRx : 4;
        const leafStrokeW = (options.leafBoxStrokeWidth !== undefined) ? options.leafBoxStrokeWidth : 1;
        const leafStrokeMode = (options.leafBoxStrokeColorMode !== undefined) ? options.leafBoxStrokeColorMode : "chemotype";

        node.filter(function (d) { return (d.depth === options.depth); })
            .each(function (d) {
                const g = d3.select(this);
                const t = g.select("text");
                if (t.empty()) return;

                const bb = t.node().getBBox();
                const tr = t.attr("transform");
                const extraTop = /[A-Z]/.test(String(d.name || "")) ? 1 : 0;

                g.insert("rect", "text")
                    .attr("transform", tr || null)
                    .attr("x", bb.x - leafPadX)
                    .attr("y", bb.y - leafPadY + leafDy + leafTopTrim - extraTop)
                    .attr("width", bb.width + leafPadX * 2)
                    .attr("height", bb.height + leafPadY * 2 - leafTopTrim - leafBottomTrim + extraTop)
                    .attr("rx", leafRx)
                    .attr("ry", leafRx)
                    .style("fill", "#fff")
                    .style("stroke", (leafStrokeMode === "black") ? "#000" : (d.color || "#999"))
                    .style("stroke-width", leafStrokeW + "px");
            });
    }

    // Calculate actual content bounds after rendering
    // Find the rightmost leaf position to reduce white space
    var maxLeafX = 0;
    node.filter(function(d) { return d.depth === options.depth; }).each(function(d) {
        // Calculate actual position including label width
        var textNode = d3.select(this).select("text").node();
        var labelWidth = textNode ? textNode.getBBox().width : 0;
        var nodeX = d.y + 10 + labelWidth; // leaf offset + label width
        if (nodeX > maxLeafX) maxLeafX = nodeX;
    });

    var bounds = svg_g.node().getBBox();
    // Use actual leaf positions to minimize right margin (minimum 20px padding)
    var rightMargin = 20;
    var actualWidth = maxLeafX + margin.left + rightMargin;
    var actualHeight = Math.max(height, bounds.height + margin.top + margin.bottom + 50);

    svg.attr("width", actualWidth)
        .attr("height", actualHeight)
        .attr('viewBox', '0 0 ' + actualWidth + ' ' + actualHeight)
        .attr('preserveAspectRatio', 'xMidYMin meet');
}

function custom_decodeHtmlEntities(text) {
    return String(text || '').replace(/&[A-Za-z0-9#]+;/g, function(entity) {
        const el = document.createElement('textarea');
        el.innerHTML = entity;
        return el.value || entity;
    });
}

// Format text with HTML (handles subscripts, italics, entities like &kappa;, etc.)
function custom_formatTextWithHTML(text) {
    // Apply all the replacements step by step
    return custom_decodeHtmlEntities(text)
        .replace(" receptor", '')
        .replace("-adrenoceptor", '')
        .replace(" receptor-", '-')
        .replace("<sub>", '</tspan><tspan baseline-shift="-0.18em" font-size="70%">')
        .replace("</sub>", '</tspan><tspan>')
        .replace("<i>", '</tspan><tspan font-style="italic">')
        .replace("</i>", '</tspan><tspan>')
        .replace("Long-wave-sensitive", 'LWS')
        .replace("Medium-wave-sensitive", 'MWS')
        .replace("Short-wave-sensitive", 'SWS')
        .replace("Olfactory", 'OLF')
        .replace("calcitonin-like receptor", 'CLR');
}

// Custom changeLeavesLabels function
function custom_changeLeavesLabels(location, value, dict, styling_dict) {
    styling_dict.starter = 0;

    let gNodes = d3.select('#' + location).selectAll('g');

    gNodes.each(function(d) {
        let g = d3.select(this);
        if (g.attr("id") !== null) {
            let name = g.attr("id").substring(1);
            let labelName;
            if (value === "UniProt") {
                labelName = name;
            } else if (dict[name]) {
                labelName = dict[name][0];
                labelName = custom_formatTextWithHTML(labelName);  // Apply HTML formatting
            } else {
                labelName = name;
            }

            let node = d3.select('#' + location).select('#X' + name);
            if (node.size() !== 0) {
                node.selectAll("text")[0].forEach(function(node_label) {
                    node_label.innerHTML = labelName;  // Use innerHTML to preserve HTML formatting
                    let labelSize = node_label.getBBox().width;
                    if (labelSize > styling_dict.starter) {
                        styling_dict.starter = labelSize;
                    }
                });
            }
        }
    });
}


// ========================
// === Initialize data  ===
// ========================
const { normKey, fnv1a32, stableColorForKey, CLASS_COLORS, CHEMOTYPE_COLORS, CHEMOTYPE_FALLBACK_PALETTE } = window.ClassificationCore;
const custom_norm_key = normKey;
const custom_hash32 = fnv1a32;
function custom_get_chemotype_color(name) {
    const key = custom_norm_key(name);
    if (CHEMOTYPE_COLORS[key]) return CHEMOTYPE_COLORS[key];
    return stableColorForKey(key, CHEMOTYPE_FALLBACK_PALETTE);
}

const TREE_DATA = window.CLASSIFICATION_TREE_DATA || {};
const tree_sets = TREE_DATA.treeSets || {};
const tree_leaf_label_lookup = TREE_DATA.leafLabelLookup || {};
const TREE_INITIAL_TYPE = TREE_DATA.initialType || "Class";
const TREE_INITIAL_SELECTION = TREE_DATA.initialSelection || "";
const TREE_LOCKED_TYPE = TREE_DATA.lockedType || "";
const TREE_LOCKED_SELECTION = TREE_DATA.lockedSelection || "";

// ========================
// === Manual scaling   ===
// ========================
// Manually set per-selection radius scaling.
// - Key: `${Type}:${SelectionKey}` (e.g. "Class:A", "Modality:Peptide receptors")
// - Value: radiusScale multiplier (>0). 1.0 means "no scaling".
//
// If a plot is not listed here, we default to 1.0.
var MANUAL_RADIUS_SCALE = {
    // Class:
    "Class:A": 1.1,
    "Class:B1": 0.6,
    "Class:B2": 0.8,
    "Class:C": 0.55,
    "Class:F": 0.24,
    "Class:T2": 0.35,
    "Class:V": 0.5,

    // Modality:
    "Modality:Orphan receptors": 0.8,
    "Modality:Peptide receptors": 0.9,
    "Modality:Protein receptors": 0.8,
    "Modality:Small molecule receptors": 1.0,

    // Chemotype (same set as CHEMOTYPE_COLORS, excluding "Ion receptors"):
    "Chemotype:Adhesion receptors": 0.5,
    "Chemotype:Alicarboxylic acid receptors": 0.4,
    "Chemotype:Aminergic receptors": 0.4,
    "Chemotype:Amino acid receptors": 0.4,
    "Chemotype:Lipid receptors": 0.5,
    "Chemotype:Melatonin receptors": 0.1,
    "Chemotype:Nucleotide receptors": 0.25,
    "Chemotype:Orphan receptors": 0.8,
    "Chemotype:Peptide receptors": 0.8,
    "Chemotype:Protein receptors": 0.5,
    "Chemotype:Retinal receptors": 0.1,
    "Chemotype:Steroid receptors": 0.2,
    "Chemotype:Tastant receptors": 0.45,
};

function getManualRadiusScale(plotKey) {
    if (!plotKey) return 1.0;
    var v = (MANUAL_RADIUS_SCALE && MANUAL_RADIUS_SCALE[plotKey] !== undefined) ? Number(MANUAL_RADIUS_SCALE[plotKey]) : 1.0;
    if (!isFinite(v) || v <= 0) return 1.0;
    return v;
}

// Function to clear all SVGs in a container
function clearContainer(containerId) {
    const container = document.getElementById(containerId);
    if (container) {
        // Remove all SVG elements using d3 to ensure proper cleanup
        d3.select('#' + containerId + '_svg').remove();
        // Also remove any other SVGs that might exist
        const svgs = container.querySelectorAll('svg');
        svgs.forEach(svg => {
            // Remove any associated d3 selections
            d3.select(svg).remove();
            svg.remove();
        });
        // Clear container content
        container.innerHTML = '';
    }
}

// Only for the Class O2 dual-plot layout (see isDualPlot in renderCurrent): stitches the two
// side-by-side SVGs into one exportable SVG, so "Download" gets both halves instead of just
// the left one. Every other class/tab has no #tree_plot_secondary_svg, so getActiveExportSvg()
// below falls straight back to the normal single-plot export for them.
function combineTwoSvgsSideBySide(svgA, svgB, gap) {
    gap = gap || 16;
    const svgNS = "http://www.w3.org/2000/svg";
    const wA = parseFloat(svgA.getAttribute("width")) || 800;
    const hA = parseFloat(svgA.getAttribute("height")) || 800;
    const wB = parseFloat(svgB.getAttribute("width")) || 800;
    const hB = parseFloat(svgB.getAttribute("height")) || 800;
    const totalW = wA + gap + wB;
    const totalH = Math.max(hA, hB);

    const combined = document.createElementNS(svgNS, "svg");
    combined.setAttribute("xmlns", svgNS);
    combined.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink");
    combined.setAttribute("width", totalW);
    combined.setAttribute("height", totalH);
    combined.setAttribute("viewBox", "0 0 " + totalW + " " + totalH);

    const bg = document.createElementNS(svgNS, "rect");
    bg.setAttribute("width", totalW);
    bg.setAttribute("height", totalH);
    bg.setAttribute("fill", "#ffffff");
    combined.appendChild(bg);

    const cloneA = svgA.cloneNode(true);
    cloneA.setAttribute("x", 0);
    cloneA.setAttribute("y", 0);
    combined.appendChild(cloneA);

    const cloneB = svgB.cloneNode(true);
    cloneB.setAttribute("x", wA + gap);
    cloneB.setAttribute("y", 0);
    combined.appendChild(cloneB);

    return combined;
}

function getActiveExportSvg() {
    const primary = document.getElementById("tree_plot_main_svg");
    const secondary = document.getElementById("tree_plot_secondary_svg");
    const wrap = document.getElementById("tree_plot_scroll_wrap");
    const isDualPlotActive = !!(secondary && wrap && wrap.classList.contains("dual-plot"));
    if (isDualPlotActive && primary) {
        return combineTwoSvgsSideBySide(primary, secondary, 16);
    }
    return primary;
}

function wireDownloadButton() {
    const pngBtn = document.getElementById("download_png_main");
    const svgBtn = document.getElementById("download_svg_main");
    function safeTitle() {
        const plotHost = document.getElementById("tree_plot_main");
        const title = (plotHost && plotHost.getAttribute("data-export-title")) || (document.getElementById("tree_plot_title") || {}).textContent || "Classification_tree";
        return String(title).trim().replace(/\s+/g, "_").replace(/[^A-Za-z0-9_\-]+/g, "") || "Classification_tree";
    }
    if (pngBtn) {
        pngBtn.addEventListener("click", function () {
            const svg = getActiveExportSvg();
            GPCRomeSvgExport.downloadSvgAsPng(svg, safeTitle() + ".png", 2);
        });
    }
    if (svgBtn) {
        svgBtn.addEventListener("click", function () {
            const svg = getActiveExportSvg();
            GPCRomeSvgExport.downloadSvg(svg, safeTitle() + ".svg");
        });
    }
}

// UI controller (single plot)
let TREE_UI = {
    type: "Class",
    selectionKey: "A",
    layout: "Tree - Circular",
    leafPills: false,
    leafLabelType: "Protein",
};

function getTypeOptions(type) {
    const set = tree_sets[type];
    if (!set || !set.options) return [];
    return set.options;
}

function getPlot(type, selectionKey) {
    const set = tree_sets[type];
    if (!set || !set.plots) return null;
    return set.plots[selectionKey] || null;
}

function setActiveType(type) {
    TREE_UI.type = type;
    const group = document.getElementById("tree_type_group");
    if (group) {
        Array.prototype.slice.call(group.querySelectorAll("button[data-type]")).forEach(function (btn) {
            const t = btn.getAttribute("data-type");
            if (t === type) {
                btn.classList.add("active");
                btn.classList.add("btn-primary");
                btn.classList.remove("btn-default");
            } else {
                btn.classList.remove("active");
                btn.classList.remove("btn-primary");
                btn.classList.add("btn-default");
            }
        });
    }
    populateSelectionDropdown(type);
}

function setActiveLayout(layout) {
    TREE_UI.layout = layout;
    const group = document.getElementById("tree_layout_group");
    const dropdownButton = document.getElementById("treeLayoutDropdownBtn");
    let activeLabel = "Circular";
    if (group) {
        Array.prototype.slice.call(group.querySelectorAll("button[data-layout]")).forEach(function (btn) {
            const v = btn.getAttribute("data-layout");
            if (v === layout) {
                btn.classList.add("active");
                btn.classList.add("btn-primary");
                btn.classList.remove("btn-default");
                btn.classList.remove("btn-outline-primary");
                activeLabel = btn.getAttribute("data-layout-label") || btn.textContent.trim() || activeLabel;
            } else {
                btn.classList.remove("active");
                btn.classList.remove("btn-primary");
                btn.classList.remove("btn-default");
                btn.classList.add("btn-outline-primary");
            }
        });
    }
    if (dropdownButton) {
        dropdownButton.innerHTML = "Layout: " + activeLabel + ' <span class="caret"></span>';
    }
    scheduleRender();
}

function syncLeafColorToggle() {
    const btn = document.getElementById("treeLeafColorToggle");
    if (!btn) return;
    const enabled = !!TREE_UI.leafPills;
    btn.classList.toggle("btn-primary", enabled);
    btn.classList.toggle("btn-danger", !enabled);
    btn.setAttribute("aria-pressed", enabled ? "true" : "false");
}

function setLeafPills(enabled) {
    TREE_UI.leafPills = !!enabled;
    syncLeafColorToggle();
    scheduleRender();
}

function syncLeafLabelDropdown() {
    const group = document.getElementById("tree_leaf_label_group");
    const dropdownButton = document.getElementById("treeLeafLabelDropdownBtn");
    const activeType = TREE_UI.leafLabelType || "UniProt";
    if (group) {
        Array.prototype.slice.call(group.querySelectorAll("button[data-label-type]")).forEach(function (btn) {
            const isActive = btn.getAttribute("data-label-type") === activeType;
            btn.classList.toggle("active", isActive);
            btn.classList.toggle("btn-primary", isActive);
            btn.classList.toggle("btn-outline-primary", !isActive);
        });
    }
    if (dropdownButton) {
        dropdownButton.innerHTML = "Receptor names: " + activeType + ' <span class="caret"></span>';
    }
}

function setLeafLabelType(labelType) {
    TREE_UI.leafLabelType = labelType || "UniProt";
    syncLeafLabelDropdown();
    scheduleRender();
}

function applyLeafLabels(root, labelType) {
    const selectedType = labelType || "UniProt";
    function walk(node) {
        if (!node) return;
        const children = node.children || [];
        if (!children.length) {
            const originalName = String(node.name || "").trim();
            const key = originalName.toUpperCase();
            const labels = tree_leaf_label_lookup[key] || {};
            node._leafKey = key || originalName;
            node.name = labels[selectedType] || originalName;
            return;
        }
        children.forEach(walk);
    }
    walk(root);
    return root;
}

function populateSelectionDropdown(type) {
    const sel = document.getElementById("tree_selection");
    if (!sel) return;
    const opts = getTypeOptions(type);
    sel.innerHTML = "";
    let visibleOpts = opts.slice();
    if (TREE_LOCKED_TYPE && TREE_LOCKED_SELECTION && type === TREE_LOCKED_TYPE) {
        visibleOpts = opts.filter(function (o) { return o.key === TREE_LOCKED_SELECTION; });
    }
    visibleOpts.forEach(function (o) {
        const opt = document.createElement("option");
        opt.value = o.key;
        opt.textContent = o.label;
        sel.appendChild(opt);
    });
    // Default to "A" for Class type; otherwise use first option
    var defaultKey = visibleOpts.length ? visibleOpts[0].key : null;
    if (TREE_LOCKED_TYPE && TREE_LOCKED_SELECTION && type === TREE_LOCKED_TYPE) {
        defaultKey = TREE_LOCKED_SELECTION;
    } else if (type === TREE_INITIAL_TYPE && TREE_INITIAL_SELECTION && visibleOpts.some(function (o) { return o.key === TREE_INITIAL_SELECTION; })) {
        defaultKey = TREE_INITIAL_SELECTION;
    } else if (type === "Class" && visibleOpts.some(function (o) { return o.key === "A"; })) {
        defaultKey = "A";
    }
    TREE_UI.selectionKey = defaultKey;
    if (defaultKey !== null) sel.value = defaultKey;
    scheduleRender();
}

function scheduleRender() {
    requestAnimationFrame(function () {
        requestAnimationFrame(function () {
            renderCurrent();
        });
    });
}

function collapse_singletons(root, labels) {
    const stacked = [];
    if (!root || !root.children || root.children.length === 0) return { data: root, stacked };
    if (!labels || !labels.length) return { data: root, stacked };
    for (let i = 0; i < labels.length; i++) {
        if (root.children.length === 1 && root.children[0] && root.children[0].children) {
            const node = root.children[0];
            stacked.push({ label: labels[i], value: node.name });
            root.children = node.children || [];
        } else {
            break;
        }
    }
    return { data: root, stacked };
}

function applyTreeColors(root, stacked_meta, options) {
    const mode = (options && options.colorMode) ? String(options.colorMode) : "chemotype";
    // Special case: Chemotype plots should be uniformly colored by the selected chemotype,
    // even if the tree contains Class/Family levels (which may collapse away).
    if (mode === "chemotype" && options && options.forceChemotype) {
        const c = custom_get_chemotype_color(String(options.forceChemotype));
        (function walk(n) {
            if (!n) return;
            n.color = c;
            if (n.children) n.children.forEach(walk);
        })(root);
        return root;
    }
    if (mode === "fixed") {
        const c = (options && options.fixedColor) ? String(options.fixedColor) : "#333";
        (function walk(n) {
            if (!n) return;
            n.color = c;
            if (n.children) n.children.forEach(walk);
        })(root);
        return root;
    }
    if (mode === "class") {
        // If Class got collapsed into stacked meta, apply that single class color to everything.
        if (stacked_meta && stacked_meta.length) {
            const hit = stacked_meta.find(x => (x.label || '').toLowerCase() === 'class');
            if (hit && hit.value) {
                let key = String(hit.value).trim();
                key = key.replace(/^Class\s+/i, '');
                key = (key.match(/^[A-Za-z0-9]+/) || [""])[0];
                const c = CLASS_COLORS[key] || "#333";
                (function walk(n) {
                    if (!n) return;
                    n.color = c;
                    if (n.children) n.children.forEach(walk);
                })(root);
                return root;
            }
        }
        // Color by class (depth 1 nodes).
        if (root && root.children && root.children.length) {
            root.children.forEach(function (clsNode) {
                let key = String(clsNode.name || "").trim();
                key = key.replace(/^Class\s+/i, '');
                key = (key.match(/^[A-Za-z0-9]+/) || [""])[0];
                const c = CLASS_COLORS[key] || "#333";
                (function walk(n) {
                    if (!n) return;
                    n.color = c;
                    if (n.children) n.children.forEach(walk);
                })(clsNode);
            });
        }
        return root;
    }
    // Default: chemotype coloring (existing logic).
    return custom_apply_tree_colors(root, stacked_meta);
}

function renderOrphanLegendPills(svgId, labels, options) {
    if (!labels || !labels.length) return;
    const svg = d3.select('#' + svgId);
    if (svg.empty()) return;

    const baseW = parseFloat(svg.attr("width")) || 800;
    const baseH = parseFloat(svg.attr("height")) || 800;

    // Layout: 20 columns, as requested
    const cols = 20;
    const cellW = baseW / cols;
    const rowH = 18;
    const headerH = 14; // room for the "Orphans" header above the pill rows
    const headerPadLeft = 4;
    const padTop = headerH + 10; // header sits right after the plot; small gap before the first pill row
    const padBottom = 10;
    const nRows = Math.ceil(labels.length / cols);
    const extraH = padTop + (nRows * rowH) + padBottom;

    const newH = baseH + extraH;
    svg.attr("height", newH);
    svg.attr("viewBox", "0 0 " + baseW + " " + newH);

    const stroke = custom_get_chemotype_color("Orphan receptors");
    const g = svg.append("g")
        .attr("class", "orphan-legend")
        .attr("transform", "translate(0," + (baseH + padTop) + ")");

    g.append("text")
        .attr("class", "orphan-legend-header")
        .attr("x", headerPadLeft)
        .attr("y", -padTop)
        .attr("dy", "0.8em")
        .attr("text-anchor", "start")
        .style("font-family", (options && options.fontFamily) ? options.fontFamily : "'Palatino Linotype', Georgia, 'Times New Roman', serif")
        .style("font-size", "11px")
        .style("font-weight", "bold")
        .style("fill", "#333")
        .text("Orphans");

    // Measure max text bbox so all pills share the same size.
    const padX = 4;
    const padY = 2;
    const fontFamily = (options && options.fontFamily) ? options.fontFamily : "'Palatino Linotype', Georgia, 'Times New Roman', serif";
    const fontSize = (options && options.fontSize && options.fontSize.receptor) ? options.fontSize.receptor : "9px";

    let maxW = 0;
    let maxH = 0;
    const meas = g.append("g").attr("transform", "translate(-99999,-99999)");
    for (let i = 0; i < labels.length; i++) {
        const t = meas.append("text")
            .attr("text-anchor", "middle")
            .attr("dy", "0.95em")
            .style("font-family", fontFamily)
            .style("font-size", fontSize)
            .text(String(labels[i]).toUpperCase());
        const bb = t.node().getBBox();
        if (bb.width > maxW) maxW = bb.width;
        if (bb.height > maxH) maxH = bb.height;
        t.remove();
    }
    meas.remove();

    // Make pills slightly narrower than the cell to leave a small gap between neighbors.
    const gapX = 4; // px total horizontal gap per cell
    const pillW = Math.max(10, Math.min(maxW + padX * 2, cellW - gapX));
    const pillH = maxH + padY * 2;

    for (let i = 0; i < labels.length; i++) {
        const col = i % cols;
        const row = Math.floor(i / cols);
        const x = (col * cellW) + (cellW / 2);
        const y = row * rowH;

        const item = g.append("g").attr("transform", "translate(" + x + "," + y + ")");
        const t = item.append("text")
            .attr("text-anchor", "middle")
            .attr("dy", "0.95em")
            .style("font-family", fontFamily)
            .style("font-size", fontSize)
            .style("fill", "#111")
            .text(String(labels[i]).toUpperCase());

        const bb = t.node().getBBox();
        item.insert("rect", "text")
            .attr("x", -pillW / 2)
            .attr("y", bb.y - padY)
            .attr("width", pillW)
            .attr("height", pillH)
            .attr("rx", 4)
            .attr("ry", 4)
            .style("fill", "#fff")
            .style("stroke", stroke)
            .style("stroke-width", "1px");
    }
}

function renderCurrent() {
    const type = TREE_UI.type;
    const selectionKey = TREE_UI.selectionKey;
    const plot = getPlot(type, selectionKey);
    const containerId = "tree_plot_main";

    clearContainer(containerId);
    const titleEl = document.getElementById("tree_plot_title");
    if (titleEl) titleEl.textContent = "";

    if (!plot || plot.error) {
        const container = document.getElementById(containerId);
        if (container) container.innerHTML = plot && plot.error
            ? `<div class="alert alert-warning">${plot.error}</div>`
            : '<div class="alert alert-info">No data available</div>';
        return;
    }

    const meta = plot.meta || {};
    const plotTitle = meta.title || selectionKey || "";
    if (titleEl) titleEl.textContent = "";
    const plotHost = document.getElementById(containerId);
    if (plotHost) plotHost.setAttribute("data-export-title", plotTitle);

    // Clone inputs so we never mutate the backend-provided JSON.
    let tree_data = JSON.parse(JSON.stringify(plot.tree));
    let tree_options = JSON.parse(JSON.stringify(plot.tree_options || {}));

    let centerLabel = "";
    if (meta.liftClassLayer) {
        const lifted = custom_lift_class_layer(tree_data);
        tree_data = lifted.data;
        centerLabel = lifted.centerLabel || "";
        // The Class layer (and its color-bearing node name) was just lifted out of tree_data
        // above, so applyTreeColors's stacked_meta/depth-1 lookups never see it -- resolve the
        // single class color here instead, from the label custom_lift_class_layer handed back.
        if (tree_options.colorMode === "class") {
            const classKey = centerLabel.replace(/^Class\s+/i, '').trim();
            tree_options = Object.assign({}, tree_options, {
                colorMode: "fixed",
                fixedColor: CLASS_COLORS[classKey] || "#333",
            });
        }
    }

    const collapsed = collapse_singletons(tree_data, meta.collapseLabels || []);
    tree_data = custom_update_tree_data(collapsed.data);
    tree_data = applyTreeColors(tree_data, collapsed.stacked, tree_options);
    tree_data = applyLeafLabels(tree_data, TREE_UI.leafLabelType);

    const maxDepth = custom_get_max_depth(tree_data, 0);
    tree_options.depth = maxDepth;
    tree_options.branch_length = custom_compute_branch_lengths(tree_data, maxDepth);

    // Common visual options
    tree_options.chemotypeLabelExtra = 65;
    tree_options.familyLabelExtra = 6;

    // Tell the renderer what the first ring represents so label offsets behave correctly.
    // - Class-lifted plots: ring1 is Chemotype (or Family for the Class C special-case).
    // - Non-lifted plots (Modality/Chemotype): ring1 is Class.
    if (meta.liftClassLayer) {
        const c0 = (meta.collapseLabels && meta.collapseLabels.length) ? String(meta.collapseLabels[0]) : "";
        tree_options.layer1Type = (c0 === "Family") ? "family" : "chemotype";
    } else {
        tree_options.layer1Type = "class";
    }

    tree_options.labelBoxPadX = 7;
    tree_options.labelBoxPadY = 1;
    tree_options.labelBoxDy = 0.33;
    // getBBox() height includes the font's full ascent (room for tall ascenders these labels
    // mostly don't have), which stacks up as empty space above the text. Trim it off the top
    // only -- purely a visual constant, expect to retune by eye.
    tree_options.labelBoxTopTrim = 2;
    tree_options.labelBoxBottomTrim = 0.5;
    // Modality and class-colored plots should have colored internal pills too.
    tree_options.labelBoxStrokeColorMode = (tree_options.colorMode === "chemotype") ? "chemotype" : "class";

    tree_options.leafPills = !!TREE_UI.leafPills;
    tree_options.leafBoxPadX = 4;
    tree_options.leafBoxPadY = 0;
    tree_options.leafBoxDy = 0;
    tree_options.leafBoxTopTrim = 2;
    tree_options.leafBoxBottomTrim = 0.5;
    tree_options.leafBoxRx = 4;
    tree_options.leafBoxStrokeWidth = 1;
    // Leaf pills should be colored for class-mode plots too.
    tree_options.leafBoxStrokeColorMode = (tree_options.colorMode === "chemotype") ? "chemotype" : "class";

    tree_options.diameterPad = 100;
    tree_options.extraPadding = 0;
    tree_options.targetSvgSize = 800;
    tree_options.stackedMetaLineGap = 6;
    tree_options.separationSibling = 1;
    tree_options.separationCousin = 1;

    tree_options.chemotypeCollapsed = custom_has_meta(collapsed.stacked, "chemotype");
    tree_options.familyCollapsed = custom_has_meta(collapsed.stacked, "family");

    tree_options.circleRadii = { chemotype: 40, family: 34, other: 30 };
    tree_options.fontSize = { 'class': "25px", 'ligandtype': "9px", 'receptorfamily': "9px", 'receptor': "9px" };
    // "Palatino" alone isn't a real font name on Windows (it ships "Palatino Linotype"
    // instead), so the browser and PowerPoint each independently guess a different
    // substitute serif font for it -- and getBBox()-sized label pills baked against one
    // guess don't match the other guess's glyph metrics once reopened elsewhere. Naming a
    // real, cross-platform-resolvable stack keeps both renderers in agreement.
    tree_options.fontFamily = "'Palatino Linotype', Georgia, 'Times New Roman', serif";

    // Dendrogram-specific options
    // Horizontal compression: < 1.0 brings root closer to leaves (e.g., 0.7 = 70% of original spacing)
    tree_options.dendrogramHorizontalCompression = 1.7;

    // Level spacing: nested dict-based approach for flexible spacing per type/depth/level combination
    // Structure: {type: {maxDepth: {level: spacing_multiplier}}}
    // Example: {"Class": {3: {1: 0.8, 2: 1.0, 3: 0.3}}, "Modality": {2: {1: 1.0, 2: 1.0}}}
    // This means:
    //   - For Class type, maxDepth=3 trees: level 1 uses 80% spacing, level 2 uses 100% spacing, level 3 uses 30% spacing
    //   - For Modality type, maxDepth=2 trees: level 1 uses 100% spacing, level 2 uses 100% spacing
    // Values < 1.0 compress spacing, > 1.0 expand spacing
    tree_options.dendrogramLevelSpacing = {
        "Class": {
            // Depth 1: Class -> Receptor (1 level)
            1: {
                1: 0.3,  // Class to Receptor
            },
            // Depth 2: Class -> Receptor family -> Receptor (2 levels)
            2: {
                1: 0.6,  // Class to Receptor family
                2: 0.3   // Receptor family to Receptor
            },
            // Depth 3: Class -> Chemotype -> Receptor family -> Receptor (3 levels)
            3: {
                1: 1.1,  // Class to Chemotype
                2: 1.0,  // Chemotype to Receptor family
                3: 0.3   // Receptor family to Receptor
            },
            // Depth 4: Class -> Chemotype -> Receptor family -> Receptor (4 levels, if exists)
            4: {
                1: 0.8,  // Root to Class
                2: 1.0,  // Class to Chemotype
                3: 0.3,  // Chemotype to Receptor family
                4: 0.3   // Receptor family to Receptor
            }
        },
        "Modality": {
            // Depth 2: Modality -> Receptor family -> Receptor (2 levels)
            2: {
                1: 0.4,  // Modality to Receptor family
                2: 0.6   // Receptor family to Receptor
            },
            // Depth 3: Modality -> Class -> Receptor family -> Receptor (3 levels)
            3: {
                1: 0.6,  // Modality to Class
                2: 1.0,  // Class to Receptor family
                3: 0.5   // Receptor family to Receptor
            }
        },
        "Chemotype": {
            // Depth 1: Chemotype -> Receptor (1 levels)
            1: {
                1: 0.3,  // Chemotype to Receptor
            },
            // Depth 2: Chemotype -> Receptor family -> Receptor (2 levels)
            2: {
                1: 0.6,  // Chemotype to Receptor family
                2: 0.4   // Receptor family to Receptor
            },
            // Depth 3: Chemotype -> Class -> Receptor family -> Receptor (3 levels)
            3: {
                1: 0.6,  // Modality to Class
                2: 1.0,  // Class to Receptor family
                3: 0.5   // Receptor family to Receptor
            }
        }
    };

    // Horizontal offsets for labels (positive = right, negative = left)
    tree_options.dendrogramChemotypeOffset = 20;  // Offset for chemotype labels (depth 1)
    tree_options.dendrogramFamilyOffset = 5;     // Offset for receptor family labels (depth 2)

    // Straight Dendrogram-specific level spacing configuration
    // Structure: {type: {maxDepth: {level: spacing_multiplier}}}
    // Similar to dendrogramLevelSpacing but specifically for straight dendrogram layout
    tree_options.dendrogramStraightLevelSpacing = {
        "Class": {
            // Depth 1: Class -> Receptor (1 level)
            1: {
                1: 0.2,  // Class to Receptor
            },
            // Depth 2: Class -> Receptor family -> Receptor (2 levels)
            2: {
                1: 0.1,  // Class to Receptor family
                2: 1.0   // Receptor family to Receptor
            },
            // Depth 3: Class -> Chemotype -> Receptor family -> Receptor (3 levels)
            3: {
                1: 0.1,  // Class to Chemotype
                2: 1.0,  // Chemotype to Receptor family
                3: 1.5   // Receptor family to Receptor
            },
            // Depth 4: Class -> Chemotype -> Receptor family -> Receptor (4 levels, if exists)
            4: {
                1: 1.0,  // Root to Class
                2: 1.0,  // Class to Chemotype
                3: 1.0,  // Chemotype to Receptor family
                4: 1.0   // Receptor family to Receptor
            }
        },
        "Modality": {
            // Depth 2: Modality -> Receptor family -> Receptor (2 levels)
            2: {
                1: 0.1,  // Modality to Receptor family
                2: 1.0   // Receptor family to Receptor
            },
            // Depth 3: Modality -> Class -> Receptor family -> Receptor (3 levels)
            3: {
                1: 0.1,  // Modality to Class
                2: 1.0,  // Class to Receptor family
                3: 1.5   // Receptor family to Receptor
            }
        },
        "Chemotype": {
            // Depth 1: Chemotype -> Receptor (1 levels)
            1: {
                1: 0.2,  // Chemotype to Receptor
            },
            // Depth 2: Chemotype -> Receptor family -> Receptor (2 levels)
            2: {
                1: 0.1,  // Chemotype to Receptor family
                2: 1.0   // Receptor family to Receptor
            },
            // Depth 3: Chemotype -> Class -> Receptor family -> Receptor (3 levels)
            3: {
                1: 0.1,  // Chemotype to Class
                2: 1.0,  // Class to Receptor family
                3: 1.5   // Receptor family to Receptor
            }
        }
    };

    // Leaf split length for straight dendrogram (horizontal segment before leaf)
    // Twice as long for Tree - Straight (was 8)
    tree_options.leafSplitLength = 16;

    // Branch stacking offset for straight dendrogram (horizontal offset to prevent overlapping)
    // Only applies from root (depth 0) to level 1 (Class/Modality/Chemotype)
    // Increased value for more pronounced stacking effect
    tree_options.branchStackOffset = 4.0;

    // Root split controls (Tree - Straight only)
    // - rootSplitBase: how far from the root (in X direction) the split should happen.
    //   Increase if you want the horizontal split “further right” from the origin.
    // - rootStackReverseAbove: reverse the ordering for the "above" half (top side).
    //   This helps align the top/bottom ladders and prevents the “back-and-forth” look.
    tree_options.rootSplitBase = 4.0;
    tree_options.rootStackReverseAbove = false;
    // - rootStemYOffsetStep: vertical separation (in px) between root branch "lanes"
    //   so the root→level1 horizontal stems don't sit on top of each other.
    tree_options.rootStemYOffsetStep = 4.0;
    // - rootStemReverseAbove: reverse the TOP stem-lane ordering (independent of stackOffset ordering).
    //   Set true if you want the top horizontal stems to be assigned in reverse order.
    tree_options.rootStemReverseAbove = true;
    // - rootLeftStub: small leftward stub (px) before the root split (cosmetic).
    tree_options.rootLeftStub = 6.0;
    // - dendrogramCenterInternalLabels: center internal pills/labels on the outgoing horizontal segment.
    tree_options.dendrogramCenterInternalLabels = true;

    // Manual scaling per selection
    const plotKey = type + ":" + String(selectionKey);
    tree_options.radiusScale = getManualRadiusScale(plotKey);

    tree_options.plotType = type; // Pass type to drawing functions for type-specific spacing

    // Route to appropriate drawing function based on layout
    const layout = TREE_UI.layout || "Tree - Circular";
    const secondaryContainerId = "tree_plot_secondary";
    const scrollWrap = document.getElementById("tree_plot_scroll_wrap");
    // Class O2 only: it's too crowded as one plot, so split it into two side-by-side
    // circular plots -- families 1-4 on the left, the rest on the right.
    const isDualPlot = (layout === "Tree - Circular" && type === "Class" && selectionKey === "O2");

    if (isDualPlot) {
        if (scrollWrap) scrollWrap.classList.add("dual-plot");
        const secondaryEl = document.getElementById(secondaryContainerId);
        if (secondaryEl) secondaryEl.style.display = "";

        const halves = splitFamilyChildrenForDualPlot(tree_data);
        custom_draw_tree(halves.left, Object.assign({}, tree_options, { anchor: containerId }), collapsed.stacked, centerLabel);
        custom_draw_tree(halves.right, Object.assign({}, tree_options, { anchor: secondaryContainerId }), collapsed.stacked, centerLabel);
    } else {
        if (scrollWrap) scrollWrap.classList.remove("dual-plot");
        clearContainer(secondaryContainerId);
        const secondaryEl = document.getElementById(secondaryContainerId);
        if (secondaryEl) secondaryEl.style.display = "none";

        tree_options.anchor = containerId;
        if (layout === "Tree - Organic") {
            custom_draw_dendrogram_curved(tree_data, tree_options, collapsed.stacked, centerLabel);
        } else if (layout === "Tree - Straight") {
            custom_draw_dendrogram_straight(tree_data, tree_options, collapsed.stacked, centerLabel);
        } else {
            custom_draw_tree(tree_data, tree_options, collapsed.stacked, centerLabel);
        }
    }

    // Class A only: append orphan pill legend under the SVG (alphabetical, unrotated)
    // Note: Orphan legend only works with Tree - Circular layout currently
    if (layout === "Tree - Circular" && type === "Class" && selectionKey === "A" && meta.orphanLeafLabels && meta.orphanLeafLabels.length) {
        renderOrphanLegendPills(containerId + "_svg", meta.orphanLeafLabels, tree_options);
    }

}

// Class O2 dual-plot split: family nodes numbered 1-4 go left, everything else goes right.
// Works off the trailing number in the family label (e.g. "Olfactory family 3"), so it's
// resilient to another family-name rewording the way "Odorant family N" -> "Olfactory family N"
// already happened once.
function splitFamilyChildrenForDualPlot(tree_data) {
    const children = (tree_data && tree_data.children) || [];
    const left = [];
    const right = [];
    children.forEach(function (child) {
        const m = /(\d+)\s*$/.exec(String((child && child.name) || "").trim());
        const num = m ? parseInt(m[1], 10) : null;
        (num !== null && num <= 4 ? left : right).push(child);
    });
    return {
        left: Object.assign({}, tree_data, { children: left }),
        right: Object.assign({}, tree_data, { children: right }),
    };
}

$(function () {
    wireDownloadButton();

    const group = document.getElementById("tree_type_group");
    if (group) {
        group.addEventListener("click", function (e) {
            const btn = e.target && e.target.closest ? e.target.closest("button[data-type]") : null;
            if (!btn) return;
            const type = btn.getAttribute("data-type");
            if (!type) return;
            setActiveType(type);
        });
    }

    const sel = document.getElementById("tree_selection");
    if (sel) {
        sel.addEventListener("change", function () {
            TREE_UI.selectionKey = sel.value;
            scheduleRender();
        });
    }

    const layoutGroup = document.getElementById("tree_layout_group");
    if (layoutGroup) {
        layoutGroup.addEventListener("click", function (e) {
            const btn = e.target && e.target.closest ? e.target.closest("button[data-layout]") : null;
            if (!btn) return;
            const layout = btn.getAttribute("data-layout");
            if (!layout || layout === TREE_UI.layout) return;
            setActiveLayout(layout);
        });
    }

    const leafColorToggle = document.getElementById("treeLeafColorToggle");
    if (leafColorToggle) {
        syncLeafColorToggle();
        leafColorToggle.addEventListener("click", function () {
            setLeafPills(!TREE_UI.leafPills);
        });
    }

    const leafLabelGroup = document.getElementById("tree_leaf_label_group");
    if (leafLabelGroup) {
        syncLeafLabelDropdown();
        leafLabelGroup.addEventListener("click", function (e) {
            const btn = e.target && e.target.closest ? e.target.closest("button[data-label-type]") : null;
            if (!btn) return;
            const labelType = btn.getAttribute("data-label-type");
            if (!labelType || labelType === TREE_UI.leafLabelType) return;
            setLeafLabelType(labelType);
        });
    }

    // Initialize default
    setActiveType(TREE_LOCKED_TYPE || TREE_INITIAL_TYPE || "Class");
    setActiveLayout(TREE_UI.layout || "Tree - Circular");
    if (TREE_INITIAL_SELECTION && TREE_UI.selectionKey !== TREE_INITIAL_SELECTION && !TREE_LOCKED_SELECTION) {
        const sel = document.getElementById("tree_selection");
        if (sel) {
            sel.value = TREE_INITIAL_SELECTION;
            TREE_UI.selectionKey = sel.value;
            scheduleRender();
        }
    }
});

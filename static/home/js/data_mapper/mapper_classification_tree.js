/* Mapper 2 classification tree skin — excerpt from Classification_tree prototype (lines ~218–1546). Requires d3 v3 (d3.svg.diagonal.radial). Anchor: typically `tree_plot`. */
var TREE_UI = window.TREE_UI || { layout: 'Tree - Circular', leafLabelType: 'Protein' };
window.TREE_UI = TREE_UI;

/** Map toolbar `fontSize.*` onto tree depths (Circular + dendrogram layouts share this). */
function custom_mapper_tree_font_px(d, options) {
    var fz = options && options.fontSize ? options.fontSize : {};
    var LD = options && options.depth != null ? options.depth : 0;
    function pick(val, fb) {
        var s = val != null && String(val).trim() !== '' ? String(val).trim() : '';
        return s || fb;
    }
    if (!d || d.depth === LD) {
        return pick(fz.receptor, '12px');
    }
    if (d.depth <= 0) {
        return '10px';
    }
    // Explicit depth→fontSizeKey map (set when tree levels are stripped)
    var dm = options && options.fontSize_depth_map;
    if (dm && dm[d.depth] != null) {
        return pick(fz[dm[d.depth]], '12px');
    }
    var md = LD;
    if (md >= 4) {
        if (d.depth === 1) { return pick(fz.class, '14px'); }
        if (d.depth === 2) { return pick(fz.ligandtype, '13px'); }
        if (d.depth === 3) { return pick(fz.receptorfamily, '11px'); }
    } else {
        if (d.depth === 1) { return pick(fz.ligandtype, '13px'); }
        if (d.depth === 2) { return pick(fz.receptorfamily, '11px'); }
    }
    return pick(fz.receptor, '12px');
}

function custom_update_tree_data(root) {
    function isLeaf(n) { return !n.children || n.children.length === 0; }

    function walk(node, depth, parent) {
        if (!node) return;

        // Default stroke color if none supplied by backend
        if (!node.color) node.color = "#333";

        // If this is the class node (child of root), shorten the label
        if (parent && parent.name === "" && node.name) {
            node.name = node.name.split(" (")[0];
            node.name = String(node.name).replace(/^Class\s+/i, '').trim();
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

function custom_tree_plain_label(label) {
    var wrapper = document.createElement("span");
    wrapper.innerHTML = custom_decodeHtmlEntities(String(label || ""));
    return (wrapper.textContent || wrapper.innerText || "").replace(/\s+/g, " ").trim();
}

function custom_tree_branch_spacing_factor(depth, maxDepth, options) {
    var map = options && options.depthSpacingFactors ? options.depthSpacingFactors : null;
    if (map && map[depth] != null && isFinite(Number(map[depth]))) {
        return Number(map[depth]);
    }
    if (maxDepth >= 4) {
        if (depth === 1) { return 0.3; }
        if (depth === 2) { return 1.5; }
        if (depth === 3) { return 0.9; }
    }
    if (maxDepth === 3) {
        if (depth === 1) { return 0.9; }
        if (depth === 2) { return 0.85; }
    }
    if (maxDepth === 2) {
        return 1.05;
    }
    return 1.15;
}

function custom_tree_branch_min_gap(depth, maxDepth, options) {
    var map = options && options.depthMinGaps ? options.depthMinGaps : null;
    if (map && map[depth] != null && isFinite(Number(map[depth]))) {
        return Number(map[depth]);
    }
    if (maxDepth >= 4) {
        if (depth === 1) { return 28; }
        if (depth === 2) { return 30; }
        if (depth === 3) { return 32; }
    }
    if (maxDepth === 3) {
        if (depth === 1) { return 34; }
        if (depth === 2) { return 30; }
    }
    if (maxDepth === 2) {
        return 46;
    }
    return 36;
}

function custom_tree_class_key(label) {
    // "Unclassified" is the class's real DB name (analogous to "Class A (Rhodopsin)"), but its
    // canonical short symbol is "U" -- not "Cl"/"Classless"/the word "Unclassified" itself.
    var raw = String(label || "").trim();
    var key = raw.startsWith("Unclassified") ? "U" : raw.split(" (")[0].replace(/^Class\s+/i, '').trim();
    return CLASS_COLORS[key] ? key : "";
}

function custom_tree_first_ring_radius(maxDepth, options) {
    if (options && options.firstRingRadius != null && isFinite(Number(options.firstRingRadius))) {
        return Number(options.firstRingRadius);
    }
    if (maxDepth >= 4) {
        return 78;
    }
    if (maxDepth === 3) {
        return 74;
    }
    if (maxDepth === 2) {
        return 62;
    }
    return 50;
}

// ------------------------------
// Coloring (Class + Chemotype)
// ------------------------------
// Class colors (used for the center badge stroke)
const CLASS_COLORS = {
    'A':  '#1f78b4',
    'B1': '#33a02c',
    'B2': '#6A3D9A',
    'C':  '#d62728',
    'F':  '#FF7F0E',
    'T2': '#F7B6D2',
    'U':  '#9e9e9e',
    'O1': '#66CDAA',
    'O2': '#3CB371',
    'V':  '#B8860B',
};

// 14-category chemotype colors (based on the legacy set in common/phylogenetic_tree.py)
// Palette: ColorBrewer Paired (12) + 2 extras.
const CHEMOTYPE_COLORS = {
  "Adhesion receptors": "#1F77B4",
  "Alicarboxylic acid receptors": "#AEC7E8",
  "Aminergic receptors": "#FF6B6B",
  "Amino acid receptors": "#98DF8A",
  "Ion receptors": "#D62728",
  "Lipid receptors": "#9467BD",
  "Melatonin receptors": "#DBDB8D",
  "Nucleotide receptors": "#9EDAE5",
  "Orphan receptors": "#a0b8ba",
  "Peptide receptors": "#2CA02C",
  "Protein receptors": "#FF7F0E",
  "Retinal receptors": "#FFBB78",
  "Steroid receptors": "#E377C2",
  "Tastant receptors": "#F7B6D2",
};

function custom_norm_key(s) {
    return (s || '').toString().trim();
}

function custom_hash32(str) {
    // FNV-1a 32-bit
    let h = 0x811c9dc5;
    for (let i = 0; i < str.length; i++) {
        h ^= str.charCodeAt(i);
        h = (h + ((h << 1) + (h << 4) + (h << 7) + (h << 8) + (h << 24))) >>> 0;
    }
    return h >>> 0;
}

const CHEMOTYPE_FALLBACK_PALETTE = [
    "#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c",
    "#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928",
    "#66c2a5","#fc8d62"
];

function custom_get_chemotype_color(name) {
    const key = custom_norm_key(name);
    if (CHEMOTYPE_COLORS[key]) return CHEMOTYPE_COLORS[key];
    // Stable fallback for unexpected/new chemotypes
    const idx = custom_hash32(key.toLowerCase()) % CHEMOTYPE_FALLBACK_PALETTE.length;
    return CHEMOTYPE_FALLBACK_PALETTE[idx];
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
            const label = custom_tree_plain_label(node.name);
            if (label.length > cur.length) longest[depth] = label;
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

    // Leaf label radial gap (pixels): room for stem dot + first data circle stacked along the tangent.
    var leaf_label_offset = isFinite(Number(options.radialLeafLabelGap)) ? Number(options.radialLeafLabelGap) : 18;

    // Calculate branch offsets (sorted numerically)
    const depthKeys = Object.keys(options.branch_length || {}).map(k => parseInt(k, 10)).filter(n => !isNaN(n)).sort((a,b) => a-b);
    for (var i = 0; i < depthKeys.length; i++) {
        var k = depthKeys[i];
        if (k === options.depth) { continue; }

        // First ring should start right after the center badge
        if (k === 1) {
            branches[k] = custom_tree_first_ring_radius(options.depth, options) + (options.firstRingExtra || 0);
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
        var scaledOffset = base_offset * custom_tree_branch_spacing_factor(k, options.depth, options);
        var minGap = custom_tree_branch_min_gap(k, options.depth, options);
        if (!(isFinite(scaledOffset) && scaledOffset > 0)) {
            scaledOffset = 0;
        }
        branch_offset = branch_offset + Math.max(minGap, scaledOffset);
        branches[k] = branch_offset;
    }
    // Increase leaf_offset to give more space for labels (prevent overlap)
    var adjusted_leaf_offset = options.leaf_offset; // Increase spacing for leaves
    branches[options.depth] = branch_offset + adjusted_leaf_offset;

    // Depth-aware fixups:
    // When options.depth === 1, the loop above skips k===depth so we never set branches[1],
    // which previously put leaves too close to the center badge. Ensure a sane leaf radius.
    if (options.depth === 1) {
        const base = custom_tree_first_ring_radius(options.depth, options) + (options.firstRingExtra || 0);
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
    var outerBranch = branches[options.depth];
    if (!(typeof outerBranch === 'number' && isFinite(outerBranch) && outerBranch > 0)) {
        outerBranch =
            Math.max(branch_offset || 0, adjusted_leaf_offset || 55, (options.centerBadgeR || 32) + (options.centerBadgePadding || 18));
        branches[options.depth] = outerBranch;
    }
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
    if (!(isFinite(scaleFactor) && scaleFactor > 0)) {
        scaleFactor = 1;
    }

    var tree = d3.layout.tree()
        .size([360, diameter / 2])
        .separation(function (a, b) {
            // Controls angular spacing between nodes.
            // Bootstrap: reduce the extra spacing between different receptor families.
            // Defaults: sibling=1, cousin=1 => no extra gap between families.
            var sib = (options.separationSibling !== undefined) ? options.separationSibling : 1;
            var cous = (options.separationCousin !== undefined) ? options.separationCousin : 1;
            var base = (a.parent === b.parent) ? sib : cous;
            // Root depth is 0; dividing by a.depth corrupts layouts (Infinity / NaN x positions).
            return base / Math.max(a.depth, 1);
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
            d.y = 0;
        } else {
            var y = branches[d.depth];
            if (typeof y !== 'number' || !isFinite(y)) {
                var outerRing = branches[options.depth];
                if (typeof outerRing !== 'number' || !isFinite(outerRing) || outerRing <= 0) {
                    outerRing = Math.max(branch_offset || 60, adjusted_leaf_offset || 40);
                    branches[options.depth] = outerRing;
                }
                var depthNorm = Math.max(options.depth, 1);
                y = outerRing * (d.depth / depthNorm);
                branches[d.depth] = y;
            }
            d.y = y;
        }
    });

    var links = tree.links(nodes);

    var link = svg_g.append("g")
        .attr("class", "links")
        .selectAll("path")
        .data(links)
        .enter().append("path")
        .each(function (d) { d.target.linkNode = this; })
        .attr("d", function (d) {
            if (!d.source || !d.target) {
                return null;
            }
            if (!(isFinite(d.source.x) && isFinite(d.source.y) && isFinite(d.target.x) && isFinite(d.target.y))) {
                return null;
            }
            return diagonal(d);
        })
        .style("stroke", function (d) { return d.target.color; })
        .style("stroke-width", function (d) { if (d.target.depth > 0) { return thickness - d.target.depth; } else { return 0; } })
        .style("fill-opacity", 0)
        .style("opacity", 1);

    var node = svg_g.selectAll(".node")
        .data(nodes)
        .enter().append("g")
        .attr("class", "node")
        .attr("transform", function (d) {
            var x = isFinite(d.x) ? d.x : 0;
            var y = isFinite(d.y) ? d.y : 0;
            if (d.name === '') {
                return "rotate(" + x + ")translate(" + y + ")";
            }
            return "rotate(" + (x - 90) + ")translate(" + y + ")";
        })

    node.filter(function (d) { return (d.depth === options.depth) })
        .attr("id", function (d) { if (d.name === '') { return "innerNode" } else { return 'X' + String(d._leafKey || d.name).toUpperCase() } });

    custom_add_leaf_end_dots(node, options);

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
            var labelOffset = d.depth === options.depth ? leaf_label_offset : 6;
            if (d.depth === options.depth) {
                return d.x < 181 ? `translate(${labelOffset})` : `rotate(180)translate(-${labelOffset})`;
            } else {
                // Internal label offset only moves the label boxes/text, not the branch geometry.
                // Tweak these options to push specific hierarchy labels outward along their branch.
                var innerOffset = 12;
                var classLabelOut = Number(options.classLabelOut || 0);
                var chemotypeLabelOut = Number(options.chemotypeLabelOut || 0);
                var familyLabelOut = Number(options.familyLabelOut || 0);
                if (d.depth === 1) {
                    innerOffset = innerOffset - (isFinite(classLabelOut) ? classLabelOut : 0);
                    innerOffset = innerOffset - (isFinite(chemotypeLabelOut) ? chemotypeLabelOut : 0);
                } else if (d.depth === 2 && options.depth >= 3) {
                    innerOffset = innerOffset - (isFinite(chemotypeLabelOut) ? chemotypeLabelOut : 0);
                } else if (d.depth === (options.depth - 1)) {
                    innerOffset = innerOffset - (isFinite(familyLabelOut) ? familyLabelOut : 0);
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
            if (d.depth === options.depth && d._labelHtml) {
                this.innerHTML = custom_formatTextWithHTML(d._labelHtml);
            } else if (d.depth > 0 && d.name && (d.depth === options.depth || /[<&]/.test(String(d.name)))) {
                this.innerHTML = custom_formatTextWithHTML(d.name);
            }
            custom_adjust_leaf_label_baseline(this, d, options);
        })
        .call(custom_wrap, options.branch_trunc)
        .style("font-size", function (d) {
            return custom_mapper_tree_font_px(d, options);
        })
        .style("font-family", function(d) {
            // Custom font family - can be modified here
            return options.fontFamily || "Palatino";
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

            g.insert("rect", "text")
                .attr("transform", tr || null)
                .attr("x", bb.x - pad.x)
                // Some fonts render slightly "high" relative to the bbox; allow tiny manual nudge.
                .attr("y", bb.y - pad.y + (options.labelBoxDy || 0))
                .attr("width", bb.width + pad.x * 2)
                .attr("height", bb.height + pad.y * 2)
                .attr("rx", 10)
                .attr("ry", 10)
                .style("fill", "#FFF")
                .style("stroke", ((options.labelBoxStrokeColorMode || "chemotype") === "black") ? "#000" : (d.color || "#999"))
                .style("stroke-width", "1px");
        });

    // Center badge removed (no middle label/pill).

    function applyCircularBoundingBoxPadding() {
        var pad = (options.circularBBoxPadding !== undefined) ? Number(options.circularBBoxPadding) : 22;
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
        var fontFamily = options.fontFamily || "Palatino";

        if (options.depth === 4) {
            if (depth === 1) {
                ctx.font = ((options.fontSize && options.fontSize.class) || '14px') + " " + fontFamily;
            } else if (depth === 2) {
                ctx.font = ((options.fontSize && options.fontSize.ligandtype) || '13px') + " " + fontFamily;
            } else if (depth === 3) {
                ctx.font = ((options.fontSize && options.fontSize.receptorfamily) || '11px') + " " + fontFamily;
            } else {
                ctx.font = ((options.fontSize && options.fontSize.receptor) || '12px') + " " + fontFamily;
            }
        } else if (options.depth === 3) {
            if (depth === 1) {
                ctx.font = ((options.fontSize && options.fontSize.ligandtype) || '13px') + " " + fontFamily;
            } else if (depth === 2) {
                ctx.font = ((options.fontSize && options.fontSize.receptorfamily) || '11px') + " " + fontFamily;
            } else {
                ctx.font = ((options.fontSize && options.fontSize.receptor) || '12px') + " " + fontFamily;
            }
        } else {
            var fb = ((options.fontSize && options.fontSize.ligandtype) || '13px');
            ctx.font = fb + " " + fontFamily;
        }

        var w = ctx.measureText(text || '').width;
        return parseInt(isFinite(w) ? w : 0, 10);
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

function custom_add_leaf_end_dots(node, options) {
    const dotRadius = (options.leafEndDotRadius !== undefined) ? Number(options.leafEndDotRadius) : 3;
    const radius = (isFinite(dotRadius) && dotRadius > 0) ? dotRadius : 3;

    node.filter(function (d) { return d.depth === options.depth; })
        .append("circle")
        .attr("class", "leaf-end-dot")
        .attr("r", radius)
        .style("fill", function (d) { return d.color || "#333"; })
        .style("stroke", "none");
}

function custom_adjust_leaf_label_baseline(textNode, datum, options) {
    if (!textNode || !datum || datum.depth !== options.depth) return;
    const label = String(datum.name || "");
    const text = d3.select(textNode);
    if (/[a-z]/.test(label) && !/<\s*sub\b/i.test(label)) {
        text.attr("dy", ".16em");
    }
    if (custom_leaf_label_needs_lowercase_nudge(datum, options)) {
        text.attr("y", custom_leaf_label_text_dy(datum, options));
    }
}

function custom_plain_leaf_label(text) {
    const wrapper = document.createElement("span");
    wrapper.innerHTML = custom_decodeHtmlEntities(String(text || ""));
    return (wrapper.textContent || wrapper.innerText || "").trim();
}

function custom_leaf_label_text_dy(datum, options) {
    if (!datum || datum.depth !== options.depth) return 0;
    if (!custom_leaf_label_needs_lowercase_nudge(datum, options)) return 0;
    const configuredTextDy = (options.leafLowercaseTextDy !== undefined) ? Number(options.leafLowercaseTextDy) : 0.25;
    return isFinite(configuredTextDy) ? configuredTextDy : 0.25;
}

function custom_leaf_label_needs_lowercase_nudge(datum, options) {
    if (!datum || datum.depth !== options.depth) return false;
    const label = custom_plain_leaf_label(datum.name);
    const compactLabel = label.replace(/\s+/g, "");
    const greekOnly = /^[\u0370-\u03FF]+$/.test(compactLabel);
    const allLowercase = /[a-z]/.test(label) && !/[A-Z]/.test(label) && label === label.toLowerCase();
    return greekOnly || allLowercase;
}

function custom_wrap(text, width) {
    if (width === 0) {
        return;
    }
    text.each(function () {
        if (this.__data__ && this.__data__._labelHtml) {
            return;
        }
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

    custom_add_leaf_end_dots(node, options);

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
            return custom_mapper_tree_font_px(d, options);
        })
        .style("font-family", function(d) {
            return options.fontFamily || "Palatino";
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

            g.insert("rect", "text")
                .attr("transform", tr || null)
                .attr("x", bb.x - pad.x)
                .attr("y", bb.y - pad.y + (options.labelBoxDy || 0))
                .attr("width", bb.width + pad.x * 2)
                .attr("height", bb.height + pad.y * 2)
                .attr("rx", 10)
                .attr("ry", 10)
                .style("fill", "#FFF")
                .style("stroke", ((options.labelBoxStrokeColorMode || "chemotype") === "black") ? "#000" : (d.color || "#999"))
                .style("stroke-width", "1px");
        });

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

    custom_add_leaf_end_dots(node, options);

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
            return custom_mapper_tree_font_px(d, options);
        })
        .style("font-family", function(d) {
            return options.fontFamily || "Palatino";
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

            g.insert("rect", "text")
                .attr("transform", tr || null)
                .attr("x", bb.x - pad.x)
                .attr("y", bb.y - pad.y + (options.labelBoxDy || 0))
                .attr("width", bb.width + pad.x * 2)
                .attr("height", bb.height + pad.y * 2)
                .attr("rx", 10)
                .attr("ry", 10)
                .style("fill", "#FFF")
                .style("stroke", ((options.labelBoxStrokeColorMode || "chemotype") === "black") ? "#000" : (d.color || "#999"))
                .style("stroke-width", "1px");
        });

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
        .replace("<sub>", '</tspan><tspan baseline-shift="-0.25em" font-size="70%">')
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

/** Collapse single-child intermediates using named labels (Mapper uses []). */
function collapse_singletons(root, labels) {
    var stacked = [];
    if (!root || !root.children || root.children.length === 0) return { data: root, stacked: stacked };
    if (!labels || !labels.length) return { data: root, stacked: stacked };
    var i;
    for (i = 0; i < labels.length; i++) {
        if (root.children.length === 1 && root.children[0] && root.children[0].children) {
            var node0 = root.children[0];
            stacked.push({ label: labels[i], value: node0.name });
            root.children = node0.children || [];
        } else {
            break;
        }
    }
    return { data: root, stacked: stacked };
}

/** Stroke colours for branch lines (classification palette). */
function applyTreeColors(root, stacked_meta, options) {
    var mode = (options && options.colorMode) ? String(options.colorMode) : "chemotype";
    if (mode === "chemotype" && options && options.forceChemotype) {
        var fc = custom_get_chemotype_color(String(options.forceChemotype));
        (function walk(n) {
            if (!n) return;
            n.color = fc;
            if (n.children) n.children.forEach(walk);
        })(root);
        return root;
    }
    if (mode === "fixed") {
        var fix = (options && options.fixedColor) ? String(options.fixedColor) : "#333";
        (function walk2(n) {
            if (!n) return;
            n.color = fix;
            if (n.children) n.children.forEach(walk2);
        })(root);
        return root;
    }
    if (mode === "class") {
        if (options && options.forceClassColor) {
            var forcedKey = custom_tree_class_key(options.forceClassColor);
            var forcedClassColor = CLASS_COLORS[forcedKey] || "#333";
            (function walkForcedClass(n) {
                if (!n) return;
                n.color = forcedClassColor;
                if (n.children) n.children.forEach(walkForcedClass);
            })(root);
            return root;
        }
        if (stacked_meta && stacked_meta.length) {
            var hit = stacked_meta.find(function(x) { return (x.label || '').toLowerCase() === 'class'; });
            if (hit && hit.value) {
                var keyHit = String(hit.value).trim();
                keyHit = keyHit.replace(/^Class\s+/i, '');
                var ch = CLASS_COLORS[keyHit] || "#333";
                (function walk3(n) {
                    if (!n) return;
                    n.color = ch;
                    if (n.children) n.children.forEach(walk3);
                })(root);
                return root;
            }
        }
        if (root && root.children && root.children.length) {
            root.children.forEach(function(clsNode) {
                // _classNodeName is set by the tree page when the class level is stripped,
                // so chemotype/RF nodes that were hoisted up retain their original class color.
                var k = custom_tree_class_key(clsNode._classNodeName || clsNode.name);
                var cc = CLASS_COLORS[k] || "#333";
                (function walk4(n) {
                    if (!n) return;
                    n.color = cc;
                    if (n.children) n.children.forEach(walk4);
                })(clsNode);
            });
        }
        return root;
    }
    return custom_apply_tree_colors(root, stacked_meta);
}

function mapperTreeManualRadiusScale() {
    return 1.0;
}

/** Expect `tree_leaf_label_lookup`: { UNIPROT: { Protein: \"…\", Gene: \"…\", UniProt: \"…\" } }. */
function applyLeafLabels(root, labelType) {
    var selectedType = labelType || "Protein";
    function walk(node) {
        if (!node) return;
        var children = node.children || [];
        if (!children.length) {
            var originalName = String(node.name || "").trim();
            var key = originalName.toUpperCase();
            var labels = (window.tree_leaf_label_lookup && window.tree_leaf_label_lookup[key]) || {};
            node._leafKey = key || originalName;
            if (selectedType === "Protein") {
                node._labelHtml = labels.ProteinHtml || "";
                node.name = labels.Protein || custom_tree_plain_label(node._labelHtml) || originalName;
            } else {
                node._labelHtml = "";
                node.name = labels[selectedType] || originalName;
            }
            return;
        }
        children.forEach(walk);
    }
    walk(root);
    return root;
}

/**
 * Draw classification-style tree into `#tree_plot` (+ `_svg`).
 * Returns merged options clone for downstream circle drawing offsets.
 */
window.mapperClassificationRedraw = function(td, optsBase, layoutLabel) {
    var tree_data = JSON.parse(JSON.stringify(td));
    var tree_options = JSON.parse(JSON.stringify(optsBase || {}));
    tree_options.anchor = tree_options.anchor || 'tree_plot';
    if (tree_data) {
        var rootClassKey = custom_tree_class_key(tree_data.name);
        if (rootClassKey) {
            tree_options.forceClassColor = rootClassKey;
        } else if (tree_data.name === "" && tree_data.children && tree_data.children.length === 1) {
            var childClassKey = custom_tree_class_key(tree_data.children[0].name);
            if (childClassKey) {
                tree_options.forceClassColor = childClassKey;
            }
        }
    }
    var collapsed = collapse_singletons(tree_data, []);
    tree_data = custom_update_tree_data(collapsed.data);
    tree_options.colorMode = tree_options.colorMode || 'class';
    tree_options.plotType = tree_options.plotType || 'Class';
    tree_data = applyTreeColors(tree_data, collapsed.stacked, tree_options);
    tree_data = applyLeafLabels(tree_data, (window.TREE_UI && TREE_UI.leafLabelType) || 'Protein');
    var maxDepth = custom_get_max_depth(tree_data, 0);
    tree_options.depth = maxDepth;
    tree_options.branch_length = custom_compute_branch_lengths(tree_data, maxDepth);
    if (tree_options.classLabelOut == null) tree_options.classLabelOut = 0;
    if (tree_options.chemotypeLabelOut == null) tree_options.chemotypeLabelOut = maxDepth >= 4 ? 14 : 10;
    if (tree_options.familyLabelOut == null) tree_options.familyLabelOut = 0;
    if (!tree_options.layer1Type) tree_options.layer1Type = 'class';
    if (tree_options.labelBoxPadX == null) tree_options.labelBoxPadX = 7;
    if (tree_options.labelBoxPadY == null) tree_options.labelBoxPadY = 1;
    if (tree_options.labelBoxDy == null) tree_options.labelBoxDy = 0.33;
    tree_options.labelBoxStrokeColorMode = (tree_options.colorMode === "chemotype") ? "chemotype" : "class";
    if (tree_options.diameterPad == null) tree_options.diameterPad = 100;
    if (tree_options.extraPadding == null) tree_options.extraPadding = 0;
    if (tree_options.targetSvgSize == null) {
        tree_options.targetSvgSize = maxDepth <= 2 ? 650 : maxDepth === 3 ? 730 : 800;
    }
    tree_options.radiusScale = mapperTreeManualRadiusScale();
    if (!tree_options.fontSize) {
        tree_options.fontSize = { 'class': "25px", 'ligandtype': "10px", 'receptorfamily': "10px", 'receptor': "10px" };
    }
    if (!tree_options.fontFamily) tree_options.fontFamily = 'Palatino';
    if (tree_options.leaf_offset == null || !isFinite(Number(tree_options.leaf_offset))) {
        tree_options.leaf_offset = maxDepth <= 2 ? 42 : maxDepth === 3 ? 50 : 58;
    }
    if (tree_options.radialLeafLabelGap == null || !isFinite(Number(tree_options.radialLeafLabelGap))) {
        tree_options.radialLeafLabelGap = 12;
    }
    if (tree_options.leafEndDotRadius == null || !isFinite(Number(tree_options.leafEndDotRadius))) {
        tree_options.leafEndDotRadius = 2;
    }
    var layout = layoutLabel || (window.TREE_UI && TREE_UI.layout) || 'Tree - Circular';
    if (layout === "Tree - Organic") {
        custom_draw_dendrogram_curved(tree_data, tree_options, [], "");
    } else if (layout === "Tree - Straight") {
        custom_draw_dendrogram_straight(tree_data, tree_options, [], "");
    } else {
        custom_draw_tree(tree_data, tree_options, [], "");
    }
    return tree_options;
};

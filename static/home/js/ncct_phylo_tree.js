/**
 * Minimal phylo tree renderer for NewClassClusterTree.
 * Uses the same phylo_library as structure clustering, with simplified styling.
 */
(function(window) {
  'use strict';

  var ncctPhylotree = null;

  function ncctConnectedLeaves(node) {
    var leaves = [];
    if (ncctPhylotree && ncctPhylotree.descendants) {
      try {
        ncctPhylotree.descendants(node).forEach(function(n) {
          if (d3.layout.phylotree.is_leafnode(n)) leaves.push(n.name);
        });
      } catch (e) {}
    }
    return leaves;
  }

  function ncctNodeStyler(element, node, colorBySymbol) {
    if (d3.layout.phylotree.is_leafnode(node)) {
      var color = (colorBySymbol && colorBySymbol[node.name]) ? colorBySymbol[node.name] : '#888';
      element.style('fill', color);
      element.selectAll('circle').style('fill', color).style('stroke', '#444');
    } else {
      element.selectAll('circle').style('fill', '#ccc').style('stroke', 'black');
    }
  }

  function ncctBranchStyler(element, node, colorBySymbol) {
    var color = '#999';
    if (colorBySymbol && node.target && d3.layout.phylotree.is_leafnode(node.target)) {
      color = colorBySymbol[node.target.name] || '#999';
    } else if (colorBySymbol && node.target) {
      var leaves = ncctConnectedLeaves(node.target);
      if (leaves.length > 0) color = colorBySymbol[leaves[0]] || '#999';
    }
    element.style('stroke', color).style('stroke-width', 2);
  }

  window.ncctRenderPhyloTree = function(newick, colorBySymbol, containerSelector, size) {
    if (typeof d3 === 'undefined' || !d3.layout.phylotree) return;

    var container = document.querySelector(containerSelector);
    if (!container) return;

    container.innerHTML = '<svg id="ncct-clustering-tree"></svg>';
    var svgEl = document.getElementById('ncct-clustering-tree');
    if (!svgEl) return;

    var plotsize = size || 580;
    if (container.offsetWidth > 0) {
      plotsize = Math.min(container.offsetWidth * 0.95, Math.max(window.innerHeight * 0.7, 580));
    }
    plotsize = Math.max(plotsize, 500);

    colorBySymbol = colorBySymbol || {};

    ncctPhylotree = d3.layout.phylotree()
      .svg(d3.select('#ncct-clustering-tree'))
      .options({
        'left-right-spacing': 'fit-to-size',
        'top-bottom-spacing': 'fit-to-size',
        'restricted-selectable': 'none',
        'selectable': false,
        'collapsible': false,
        'transitions': false,
        'show-scale': false,
        'align-tips': true,
        'brush': false,
        'reroot': true,
        'hide': false,
        'zoom': false,
        'inner_spacing': 2
      })
      .radial(true)
      .node_span('equal')
      .size([plotsize, plotsize])
      .separation(function() { return 0.1; });

    ncctPhylotree(newick)
      .style_nodes(function(el, node) { ncctNodeStyler(el, node, colorBySymbol); })
      .style_edges(function(el, edge) { ncctBranchStyler(el, edge, colorBySymbol); })
      .layout(false);

    d3.layout.phylotree.trigger_refresh(ncctPhylotree);
    d3.layout.phylotree.trigger_refresh(ncctPhylotree);
  };
})(window);

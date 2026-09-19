// Tree (phylogenetic) tab bootstrap/glue for the receptor-family Detail page. The tree-drawing
// engine itself lives in the vendor phylo_library.js/phylotree.js/classification_phylotree.js
// files (shared with the phylogenetic_trees app, not moved here) — this file only wires the
// leaf-label dropdown, the empty-state message, and the tab-shown resize hook around it.
//
// Note: window.classificationPhyloTree* config flags and window.data stay inline in the
// template, not here — classification_phylotree.js reads them synchronously the moment its
// <script> tag runs, and it loads before this file.
(function () {
  const DATA = window.CLASSIFICATION_FAMILY_TREE_DATA || {};

  var treeEl = document.getElementById("tree-container");
  var treeLeafLabelGroupEl = document.getElementById("cvFamilyTreeLeafLabelGroup");
  var treeLeafLabelBtnEl = document.getElementById("cvFamilyTreeLeafLabelBtn");

  if ((!DATA || !DATA.tree) && treeEl) {
    treeEl.innerHTML = '<div class="cv-family-empty">No persisted family tree payload is available for this receptor family yet.</div>';
  }

  function syncTreeLeafLabelDropdown(labelType) {
    var activeType = labelType || window.classificationPhyloTreeLeafLabelType || "Protein";
    if (treeLeafLabelGroupEl) {
      Array.prototype.forEach.call(treeLeafLabelGroupEl.querySelectorAll("button[data-label-type]"), function(button) {
        var isActive = button.getAttribute("data-label-type") === activeType;
        button.classList.toggle("active", isActive);
        button.classList.toggle("btn-primary", isActive);
        button.classList.toggle("btn-outline-primary", !isActive);
      });
    }
    if (treeLeafLabelBtnEl) {
      treeLeafLabelBtnEl.innerHTML = "Receptor names: " + activeType + ' <span class="caret"></span>';
    }
  }

  syncTreeLeafLabelDropdown(window.classificationPhyloTreeLeafLabelType || "Protein");
  if (treeLeafLabelGroupEl) {
    treeLeafLabelGroupEl.addEventListener("click", function(event) {
      var button = event.target && event.target.closest ? event.target.closest("button[data-label-type]") : null;
      var labelType;
      if (!button) {
        return;
      }
      labelType = button.getAttribute("data-label-type") || "Protein";
      window.classificationPhyloTreeLeafLabelType = labelType;
      syncTreeLeafLabelDropdown(labelType);
      if (typeof window.classificationPhyloTreeSetLeafLabelType === "function") {
        window.classificationPhyloTreeSetLeafLabelType(labelType);
      }
    });
  }

  if (window.jQuery) {
    window.jQuery('a[data-toggle="tab"]').on("shown.bs.tab", function(event) {
      var href = event && event.target ? event.target.getAttribute("href") : "";
      if (href === "#cv-family-tree-tab" && typeof window.resizeTree === "function") {
        window.setTimeout(window.resizeTree, 0);
        window.setTimeout(window.resizeTree, 250);
      }
    });
  }
})();

// Matrix tab (receptor x receptor identity/similarity heatmap) for the receptor-family Detail
// page. Reuses the entities/matrix arrays from the same payload that feeds the Tree tab
// (window.CLASSIFICATION_FAMILY_MATRIX_DATA), since both tabs describe the same receptor set.
(function () {
  const { escapeHtml } = window.ClassificationCore;
  const DATA = window.CLASSIFICATION_FAMILY_MATRIX_DATA || {};

  var matrixScrollEl = document.getElementById("cvFamilyMatrixScroll");
  var matrixTableEl = document.getElementById("cvFamilyMatrixTable");
  var matrixEmptyEl = document.getElementById("cvFamilyMatrixEmpty");
  var matrixHoverEl = document.getElementById("cvFamilyMatrixHover");
  var matrixLeafLabelGroupEl = document.getElementById("cvFamilyMatrixLeafLabelGroup");
  var matrixLeafLabelBtnEl = document.getElementById("cvFamilyMatrixLeafLabelBtn");
  var matrixLeafLabelType = "Protein";

  function decodeHtmlEntities(value) {
    return String(value == null ? "" : value).replace(/&[A-Za-z0-9#]+;/g, function(entity) {
      var el = document.createElement("textarea");
      el.innerHTML = entity;
      return el.value || entity;
    });
  }

  function formattedLabelHtml(value) {
    return escapeHtml(decodeHtmlEntities(value))
      .replace(/&lt;(\/?)sub&gt;/gi, "<$1sub>")
      .replace(/&lt;(\/?)i&gt;/gi, "<$1i>");
  }

  function syncMatrixLeafLabelDropdown(labelType) {
    var activeType = labelType || matrixLeafLabelType || "Protein";
    if (matrixLeafLabelGroupEl) {
      Array.prototype.forEach.call(matrixLeafLabelGroupEl.querySelectorAll("button[data-label-type]"), function(button) {
        var isActive = button.getAttribute("data-label-type") === activeType;
        button.classList.toggle("active", isActive);
        button.classList.toggle("btn-primary", isActive);
        button.classList.toggle("btn-outline-primary", !isActive);
      });
    }
    if (matrixLeafLabelBtnEl) {
      matrixLeafLabelBtnEl.innerHTML = "Receptor names: " + activeType + ' <span class="caret"></span>';
    }
  }

  function matrixEntityLabel(entity) {
    if (!entity) {
      return "";
    }
    if (matrixLeafLabelType === "Gene") {
      return String(entity.gene_label || entity.short_label || entity.symbol || "");
    }
    if (matrixLeafLabelType === "UniProt") {
      return String(entity.short_label || entity.symbol || "");
    }
    return String(entity.name || entity.short_label || entity.symbol || "");
  }

  function matrixMetricValue(cell, metric) {
    var raw = cell && cell[metric];
    if (raw === null || raw === undefined || raw === "") {
      return null;
    }
    var value = Number(raw);
    return isFinite(value) ? value : null;
  }

  function matrixMetricDisplay(value) {
    if (value === null || value === undefined || !isFinite(value)) {
      return "";
    }
    var rounded = Math.round(value * 10) / 10;
    return Math.round(rounded) === rounded ? String(Math.round(rounded)) : rounded.toFixed(1);
  }

  function matrixHeatColor(value, range, type) {
    var t;
    var hue;
    var lightness;
    if (value === null || !range || range.min === null || range.max === null) {
      return "";
    }
    t = range.max === range.min ? 0.6 : (value - range.min) / (range.max - range.min);
    t = Math.max(0, Math.min(1, t));
    hue = type === "identity" ? 130 : 210;
    lightness = 96 - (t * 26);
    return "hsl(" + hue + ",60%," + lightness + "%)";
  }

  function buildMatrixRanges(matrix) {
    var ranges = {
      identity: { min: Infinity, max: -Infinity },
      similarity: { min: Infinity, max: -Infinity }
    };
    matrix.forEach(function(row, i) {
      (row || []).forEach(function(cell, j) {
        var type;
        var value;
        if (i === j) {
          return;
        }
        type = i < j ? "identity" : "similarity";
        value = matrixMetricValue(cell, type);
        if (value === null) {
          return;
        }
        ranges[type].min = Math.min(ranges[type].min, value);
        ranges[type].max = Math.max(ranges[type].max, value);
      });
    });
    ["identity", "similarity"].forEach(function(type) {
      if (ranges[type].min === Infinity) {
        ranges[type].min = null;
        ranges[type].max = null;
      }
    });
    return ranges;
  }

  function matrixMetricLabel(metric) {
    return metric === "identity" ? "Identity" : "Similarity";
  }

  function matrixMetricClass(metric) {
    return metric === "identity" ? "is-identity" : "is-similarity";
  }

  function hideFamilyMatrixHover() {
    if (!matrixHoverEl) {
      return;
    }
    matrixHoverEl.style.display = "none";
    matrixHoverEl.setAttribute("aria-hidden", "true");
  }

  function positionFamilyMatrixHover(event) {
    var width;
    var height;
    var left;
    var top;
    if (!matrixHoverEl || matrixHoverEl.style.display === "none") {
      return;
    }
    width = matrixHoverEl.offsetWidth || 260;
    height = matrixHoverEl.offsetHeight || 120;
    left = (event.clientX || 0) + 16;
    top = (event.clientY || 0) + 16;
    if (left + width + 12 > window.innerWidth) {
      left = (event.clientX || 0) - width - 16;
    }
    if (top + height + 12 > window.innerHeight) {
      top = (event.clientY || 0) - height - 16;
    }
    matrixHoverEl.style.left = Math.max(8, left) + "px";
    matrixHoverEl.style.top = Math.max(8, top) + "px";
  }

  function showFamilyMatrixHover(event, sourceEntity, targetEntity, cell, primaryMetric) {
    var metricOrder = primaryMetric === "identity" ? ["identity", "similarity"] : ["similarity", "identity"];
    var sourceLabel = matrixEntityLabel(sourceEntity);
    var targetLabel = matrixEntityLabel(targetEntity);
    if (!matrixHoverEl) {
      return;
    }
    matrixHoverEl.innerHTML = (
      '<div class="cv-family-matrix-hover-title">' +
        '<span class="cv-family-matrix-hover-name">' + formattedLabelHtml(sourceLabel) + '</span>' +
        '<span>|</span>' +
        '<span class="cv-family-matrix-hover-name">' + formattedLabelHtml(targetLabel) + '</span>' +
      '</div>' +
      '<div class="cv-family-matrix-hover-metrics">' +
        metricOrder.map(function(metric) {
          var value = matrixMetricValue(cell, metric);
          return (
            '<div class="cv-family-matrix-hover-row ' + matrixMetricClass(metric) + '">' +
              '<span class="cv-family-matrix-hover-badge">' + matrixMetricLabel(metric) + '</span>' +
              '<strong>' + (value === null ? "n/a" : matrixMetricDisplay(value) + "%") + '</strong>' +
            '</div>'
          );
        }).join("") +
      '</div>'
    );
    matrixHoverEl.style.display = "block";
    matrixHoverEl.setAttribute("aria-hidden", "false");
    positionFamilyMatrixHover(event);
  }

  function renderFamilyMatrix() {
    var entities = Array.isArray(DATA && DATA.entities) ? DATA.entities : [];
    var matrix = Array.isArray(DATA && DATA.matrix) ? DATA.matrix : [];
    var ranges = buildMatrixRanges(matrix);

    if (!matrixTableEl || !matrixScrollEl || !matrixEmptyEl) {
      return;
    }
    matrixTableEl.querySelector("thead").innerHTML = "";
    matrixTableEl.querySelector("tbody").innerHTML = "";

    if (!entities.length || !matrix.length) {
      matrixScrollEl.style.display = "none";
      matrixEmptyEl.style.display = "";
      return;
    }

    matrixScrollEl.style.display = "";
    matrixEmptyEl.style.display = "none";

    var headerRow = document.createElement("tr");
    var corner = document.createElement("th");
    corner.className = "cv-family-matrix-corner";
    headerRow.appendChild(corner);

    entities.forEach(function(entity) {
      var th = document.createElement("th");
      var span = document.createElement("span");
      span.className = "cv-family-matrix-col-label";
      span.innerHTML = formattedLabelHtml(matrixEntityLabel(entity));
      th.title = matrixEntityLabel(entity);
      th.appendChild(span);
      headerRow.appendChild(th);
    });
    matrixTableEl.querySelector("thead").appendChild(headerRow);

    entities.forEach(function(entity, i) {
      var tr = document.createElement("tr");
      var rowHeader = document.createElement("th");
      rowHeader.className = "cv-family-matrix-row-label";
      rowHeader.innerHTML = formattedLabelHtml(matrixEntityLabel(entity));
      rowHeader.title = matrixEntityLabel(entity);
      tr.appendChild(rowHeader);

      entities.forEach(function(otherEntity, j) {
        var td = document.createElement("td");
        var cell = matrix[i] && matrix[i][j] ? matrix[i][j] : {};
        var metric = i < j ? "identity" : "similarity";
        var value = i === j ? null : matrixMetricValue(cell, metric);
        var displayValue = matrixMetricDisplay(value);

        if (i === j) {
          td.className = "cv-family-matrix-cell cv-family-matrix-diagonal";
        } else if (value === null) {
          td.className = "cv-family-matrix-cell cv-family-matrix-empty";
          td.textContent = "-";
        } else {
          td.className = "cv-family-matrix-cell " + (metric === "identity" ? "is-identity" : "is-similarity");
          td.textContent = displayValue;
          td.style.backgroundColor = matrixHeatColor(value, ranges[metric], metric);
          td.addEventListener("mouseenter", function(event) {
            showFamilyMatrixHover(event, entity, otherEntity, cell, metric);
          });
          td.addEventListener("mousemove", positionFamilyMatrixHover);
          td.addEventListener("mouseleave", hideFamilyMatrixHover);
        }
        tr.appendChild(td);
      });

      matrixTableEl.querySelector("tbody").appendChild(tr);
    });
  }

  syncMatrixLeafLabelDropdown(matrixLeafLabelType);
  if (matrixLeafLabelGroupEl) {
    matrixLeafLabelGroupEl.addEventListener("click", function(event) {
      var button = event.target && event.target.closest ? event.target.closest("button[data-label-type]") : null;
      var labelType;
      if (!button) {
        return;
      }
      labelType = button.getAttribute("data-label-type") || "Protein";
      matrixLeafLabelType = labelType;
      hideFamilyMatrixHover();
      syncMatrixLeafLabelDropdown(labelType);
      renderFamilyMatrix();
    });
  }

  renderFamilyMatrix();

  if (matrixScrollEl) {
    matrixScrollEl.addEventListener("scroll", hideFamilyMatrixHover, { passive: true });
  }
})();

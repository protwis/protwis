// Cluster tab (sequence-similarity t-SNE scatter + similarity side panel) for the receptor-family
// Detail page. Reads window.CLASSIFICATION_FAMILY_CLUSTER_DATA, which carries the cluster payload
// (points/meta) plus the entities/matrix arrays shared with the Tree tab (used for the hover-card
// similarity lookups) and the page title (used to name downloaded images).
(function () {
  const { escapeHtml } = window.ClassificationCore;
  const DATA = window.CLASSIFICATION_FAMILY_CLUSTER_DATA || {};
  const familyClusterPayload = DATA;
  const familyEntities = Array.isArray(DATA.entities) ? DATA.entities : [];
  const familyMatrix = Array.isArray(DATA.matrix) ? DATA.matrix : [];

  var plotEl = document.getElementById("cvFamilyPlot");
  var clusterLabelsEl = document.getElementById("cvFamilyClusterLabels");
  var clusterHoverEl = document.getElementById("cvFamilyClusterHover");
  var clusterMetaEl = document.getElementById("cvFamilyClusterMeta");
  var clusterSimilarityTitleEl = document.getElementById("cvFamilyClusterSimilarityTitle");

  var familyEntityIndexBySymbol = {};
  var activeClusterSymbol = null;
  var hoveredClusterSymbol = null;
  var frozenClusterSymbol = null;
  var suppressPlotBackgroundClick = false;

  familyEntities.forEach(function(entityRow, index) {
    var symbol = entityRow && entityRow.symbol ? String(entityRow.symbol) : "";
    if (symbol) {
      familyEntityIndexBySymbol[symbol] = index;
    }
  });

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

  function downloadClusterImage(format) {
    if (!plotEl || !plotEl.data || typeof Plotly === "undefined") return;
    Plotly.toImage(plotEl, {
      format: format,
      width: 1400,
      height: 900,
      scale: format === "png" ? 2 : 1
    }).then(function(url) {
      var link = document.createElement("a");
      link.href = url;
      link.download = String(DATA.pageTitle || "").trim().replace(/\s+/g, "_").replace(/[^A-Za-z0-9_\-]+/g, "") + "_cluster." + format;
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
    });
  }

  function familyEntityBySymbol(symbol) {
    var entityIndex = familyEntityIndexBySymbol[String(symbol || "")];
    if (entityIndex === undefined) {
      return null;
    }
    return familyEntities[entityIndex] || null;
  }

  function familyMatrixEntry(sourceSymbol, targetSymbol) {
    var sourceIndex = familyEntityIndexBySymbol[String(sourceSymbol || "")];
    var targetIndex = familyEntityIndexBySymbol[String(targetSymbol || "")];
    if (sourceIndex === undefined || targetIndex === undefined) {
      return null;
    }
    if (!familyMatrix[sourceIndex] || !familyMatrix[sourceIndex][targetIndex]) {
      return null;
    }
    return familyMatrix[sourceIndex][targetIndex];
  }

  function familyClusterPointBySymbol(symbol) {
    var match = null;
    (familyClusterPayload.points || []).forEach(function(point, index) {
      if (!match && point.entry_name === String(symbol || "")) {
        match = {
          point: point,
          index: index
        };
      }
    });
    return match;
  }

  function openFamilyClusterProtein(symbol) {
    var pointMatch = familyClusterPointBySymbol(symbol);
    if (!pointMatch || !pointMatch.point || !pointMatch.point.protein_url) {
      return;
    }
    window.open(pointMatch.point.protein_url, "_blank", "noopener");
  }

  function buildFamilyMarkerSizes(points, activeSymbol) {
    return points.map(function(point) {
      return point.entry_name === activeSymbol ? 28 : 18;
    });
  }

  function buildFamilyMarkerOpacities(points, activeSymbol) {
    return points.map(function(point) {
      if (!activeSymbol) {
        return 1;
      }
      return point.entry_name === activeSymbol ? 1 : 0.72;
    });
  }

  function updateFamilySimilarityTitle(symbol) {
    var entity = familyEntityBySymbol(symbol);
    if (!clusterSimilarityTitleEl || !entity) {
      if (clusterSimilarityTitleEl) {
        clusterSimilarityTitleEl.textContent = "Similarity ranking";
      }
      return;
    }
    clusterSimilarityTitleEl.innerHTML = formattedLabelHtml(entity.short_label || entity.symbol || "")
      + ' similarity ranking'
      + (entity.subtitle ? '<span class="cv-family-cluster-side-subtitle">' + formattedLabelHtml(entity.subtitle) + '</span>' : '');
  }

  function renderFamilyClusterHoverCard(symbol) {
    var entity = familyEntityBySymbol(symbol);
    if (!clusterHoverEl) {
      return;
    }
    if (!entity) {
      clusterHoverEl.innerHTML = '<p class="cv-family-cluster-hover-empty">Hover a receptor point to inspect its sequence similarity against the rest of the family.</p>';
      return;
    }

    var peers = [];
    familyEntities.forEach(function(otherEntity) {
      var otherSymbol = otherEntity && otherEntity.symbol ? String(otherEntity.symbol) : "";
      if (!otherSymbol || otherSymbol === entity.symbol) {
        return;
      }
      var cell = familyMatrixEntry(entity.symbol, otherSymbol) || {};
      peers.push({
        symbol: otherSymbol,
        shortLabel: otherEntity.short_label || otherSymbol,
        name: otherEntity.name || otherSymbol,
        color: otherEntity.color || "#808080",
        similarityDisplay: cell.similarity_display || "n/a",
        similarity: cell.similarity
      });
    });

    peers.sort(function(a, b) {
      var similarityA = a.similarity === null || a.similarity === undefined ? -Infinity : Number(a.similarity);
      var similarityB = b.similarity === null || b.similarity === undefined ? -Infinity : Number(b.similarity);
      return similarityB - similarityA;
    });

    clusterHoverEl.innerHTML = '<div class="cv-family-cluster-hover-list">' + peers.map(function(peer) {
      return (
        '<div class="cv-family-cluster-hover-row">' +
          '<div class="cv-family-cluster-hover-label">' +
            '<span class="cv-family-cluster-label-swatch" style="background:' + escapeHtml(peer.color) + ';"></span>' +
            '<div class="cv-family-cluster-hover-text">' +
              '<strong>' + formattedLabelHtml(peer.name) + '</strong>' +
            '</div>' +
          '</div>' +
          '<div class="cv-family-cluster-hover-metric">' +
            '<div class="cv-family-cluster-hover-metric-label">Similarity</div>' +
            '<div class="cv-family-cluster-hover-metric-value">' + escapeHtml(peer.similarityDisplay) + '</div>' +
          '</div>' +
        '</div>'
      );
    }).join("") + '</div>';
  }

  function setActiveFamilyClusterLabel(symbol) {
    var activeButton = null;
    if (!clusterLabelsEl) {
      return;
    }
    Array.prototype.forEach.call(clusterLabelsEl.querySelectorAll(".cv-family-cluster-label"), function(button) {
      var isActive = button.getAttribute("data-symbol") === String(symbol || "");
      button.classList.toggle("is-active", isActive);
      if (isActive) {
        activeButton = button;
      }
    });
    return activeButton;
  }

  function updateFamilyClusterPointHighlight(symbol) {
    if (!plotEl || !plotEl.data || typeof Plotly === "undefined" || !familyClusterPayload.points || !familyClusterPayload.points.length) {
      return;
    }
    Plotly.restyle(plotEl, {
      "marker.size": [buildFamilyMarkerSizes(familyClusterPayload.points, symbol)],
      "marker.opacity": [buildFamilyMarkerOpacities(familyClusterPayload.points, symbol)],
      "marker.line.width": [familyClusterPayload.points.map(function(point) { return point.entry_name === symbol ? 2.5 : 1; })]
    }, [0]);
  }

  function clearFrozenClusterHover() {
    if (!plotEl || typeof Plotly === "undefined" || !Plotly.Fx || typeof Plotly.Fx.unhover !== "function") {
      return;
    }
    Plotly.Fx.unhover(plotEl);
  }

  function syncFrozenClusterHover() {
    var pointMatch;
    if (!frozenClusterSymbol || !plotEl || typeof Plotly === "undefined" || !Plotly.Fx || typeof Plotly.Fx.hover !== "function") {
      return;
    }
    pointMatch = familyClusterPointBySymbol(frozenClusterSymbol);
    if (!pointMatch) {
      return;
    }
    Plotly.Fx.hover(plotEl, [{ curveNumber: 0, pointNumber: pointMatch.index }]);
  }

  function clearFrozenFamilyClusterSelection() {
    frozenClusterSymbol = null;
    clearFrozenClusterHover();
    if (hoveredClusterSymbol && familyEntityBySymbol(hoveredClusterSymbol)) {
      setActiveFamilyClusterSymbol(hoveredClusterSymbol, { scrollLabelIntoView: true });
      return;
    }
    if (activeClusterSymbol && familyEntityBySymbol(activeClusterSymbol)) {
      setActiveFamilyClusterSymbol(activeClusterSymbol);
      return;
    }
    if (familyClusterPayload.points && familyClusterPayload.points.length) {
      setActiveFamilyClusterSymbol(familyClusterPayload.points[0].entry_name);
    }
  }

  function toggleFrozenFamilyClusterSelection(symbol) {
    if (!symbol) {
      return;
    }
    if (frozenClusterSymbol === symbol) {
      clearFrozenFamilyClusterSelection();
      return;
    }
    frozenClusterSymbol = String(symbol);
    setActiveFamilyClusterSymbol(frozenClusterSymbol, { scrollLabelIntoView: true });
    window.setTimeout(syncFrozenClusterHover, 0);
  }

  function setActiveFamilyClusterSymbol(symbol, options) {
    var settings = options || {};
    var activeButton;
    if (!familyEntityBySymbol(symbol)) {
      return;
    }
    activeClusterSymbol = String(symbol);
    activeButton = setActiveFamilyClusterLabel(activeClusterSymbol);
    if (settings.scrollLabelIntoView && activeButton && typeof activeButton.scrollIntoView === "function") {
      activeButton.scrollIntoView({ block: "nearest", inline: "nearest" });
    }
    updateFamilySimilarityTitle(activeClusterSymbol);
    renderFamilyClusterHoverCard(activeClusterSymbol);
    updateFamilyClusterPointHighlight(activeClusterSymbol);
  }

  function renderFamilyClusterLabels() {
    if (!clusterLabelsEl) {
      return;
    }
    clusterLabelsEl.innerHTML = "";
    familyClusterPayload.points.forEach(function(point) {
      var button = document.createElement("button");
      button.type = "button";
      button.className = "cv-family-cluster-label";
      button.setAttribute("data-symbol", point.entry_name);
      button.innerHTML =
        '<span class="cv-family-cluster-label-marker">' +
          '<span class="cv-family-cluster-label-swatch" style="background:' + escapeHtml(point.color || "#808080") + ';"></span>' +
          '<span class="cv-family-cluster-label-symbol">' + formattedLabelHtml(point.label || "") + '</span>' +
        '</span>' +
        '<span class="cv-family-cluster-label-text">' +
          '<strong>' + formattedLabelHtml(point.display_name || point.label || "") + '</strong>' +
          '<span class="cv-family-cluster-secondary">' + formattedLabelHtml(point.class_label || "") + '</span>' +
        '</span>';
      button.addEventListener("mouseenter", function() {
        if (frozenClusterSymbol) {
          return;
        }
        setActiveFamilyClusterSymbol(point.entry_name);
      });
      button.addEventListener("focus", function() {
        if (frozenClusterSymbol) {
          return;
        }
        setActiveFamilyClusterSymbol(point.entry_name);
      });
      button.addEventListener("click", function(event) {
        event.preventDefault();
        event.stopPropagation();
        toggleFrozenFamilyClusterSelection(point.entry_name);
      });
      button.addEventListener("contextmenu", function(event) {
        event.preventDefault();
        event.stopPropagation();
        openFamilyClusterProtein(point.entry_name);
      });
      clusterLabelsEl.appendChild(button);
    });
  }

  function renderCluster() {
    if (!plotEl || typeof Plotly === "undefined") return;
    if (!familyClusterPayload.points || !familyClusterPayload.points.length) {
      plotEl.innerHTML = '<div class="cv-family-empty">No receptor-family cluster data is available for this family.</div>';
      if (clusterLabelsEl) clusterLabelsEl.innerHTML = "";
      if (clusterHoverEl) {
        clusterHoverEl.innerHTML = '<p class="cv-family-cluster-hover-empty">No receptor-family similarity data is available for this family.</p>';
      }
      if (clusterMetaEl) clusterMetaEl.textContent = "";
      return;
    }

    var trace = {
      type: "scatter",
      mode: "markers",
      x: familyClusterPayload.points.map(function(point) { return point.x; }),
      y: familyClusterPayload.points.map(function(point) { return point.y; }),
      customdata: familyClusterPayload.points.map(function(point) {
        return [
          point.entry_name,
          point.display_name,
          point.class_label,
          (point.chemotypes || []).join(", "),
          (point.modality_groups || []).join(", "),
          point.protein_url
        ];
      }),
      hovertemplate: [
        "<b>%{text}</b>",
        "%{customdata[1]}",
        "Class: %{customdata[2]}",
        "Chemotype: %{customdata[3]}",
        "Modality: %{customdata[4]}",
        "<extra></extra>"
      ].join("<br>"),
      text: familyClusterPayload.points.map(function(point) { return point.label; }),
      marker: {
        size: buildFamilyMarkerSizes(familyClusterPayload.points, activeClusterSymbol),
        color: familyClusterPayload.points.map(function(point) { return point.color; }),
        opacity: buildFamilyMarkerOpacities(familyClusterPayload.points, activeClusterSymbol),
        line: { width: 1, color: "rgba(0,0,0,0.25)" }
      }
    };

    var layout = {
      margin: { l: 24, r: 24, t: 12, b: 36 },
      paper_bgcolor: "#ffffff",
      plot_bgcolor: "#ffffff",
      hovermode: "closest",
      showlegend: false,
      xaxis: { zeroline: false, showgrid: false, showticklabels: false, title: "", fixedrange: false },
      yaxis: { zeroline: false, showgrid: false, showticklabels: false, title: "", fixedrange: false },
      dragmode: "zoom",
      font: { family: "Source Code Pro, monospace", color: "#111" }
    };

    Plotly.newPlot(plotEl, [trace], layout, {
      displayModeBar: true,
      responsive: true,
      scrollZoom: true,
      doubleClick: "reset"
    }).then(function() {
      renderFamilyClusterLabels();
      if (typeof plotEl.removeAllListeners === "function") {
        plotEl.removeAllListeners("plotly_hover");
        plotEl.removeAllListeners("plotly_unhover");
        plotEl.removeAllListeners("plotly_click");
      }
      plotEl.on("plotly_hover", function(event) {
        var point = event && event.points && event.points.length ? event.points[0] : null;
        if (!point || !point.customdata || !point.customdata.length) return;
        hoveredClusterSymbol = String(point.customdata[0]);
        if (frozenClusterSymbol) {
          window.setTimeout(syncFrozenClusterHover, 0);
          return;
        }
        setActiveFamilyClusterSymbol(hoveredClusterSymbol, { scrollLabelIntoView: true });
      });
      plotEl.on("plotly_unhover", function() {
        hoveredClusterSymbol = null;
        if (frozenClusterSymbol) {
          window.setTimeout(syncFrozenClusterHover, 0);
          return;
        }
        if (!activeClusterSymbol && familyClusterPayload.points.length) {
          setActiveFamilyClusterSymbol(familyClusterPayload.points[0].entry_name);
          return;
        }
        if (activeClusterSymbol) {
          setActiveFamilyClusterSymbol(activeClusterSymbol);
        }
      });
      plotEl.on("plotly_click", function(event) {
        var point = event && event.points && event.points.length ? event.points[0] : null;
        var symbol = point && point.customdata && point.customdata.length ? String(point.customdata[0]) : "";
        if (!symbol) {
          return;
        }
        suppressPlotBackgroundClick = true;
        window.setTimeout(function() {
          suppressPlotBackgroundClick = false;
        }, 0);
        toggleFrozenFamilyClusterSelection(symbol);
      });
      plotEl.onclick = function(event) {
        var target = event && event.target ? event.target : null;
        if (suppressPlotBackgroundClick || !frozenClusterSymbol || !target) {
          return;
        }
        if (target.closest && target.closest(".modebar")) {
          return;
        }
        clearFrozenFamilyClusterSelection();
      };
      plotEl.onmousedown = function(event) {
        if (!event || event.button !== 2) {
          return;
        }
        event.preventDefault();
        event.stopPropagation();
      };
      plotEl.oncontextmenu = function(event) {
        var targetSymbol = frozenClusterSymbol || hoveredClusterSymbol || activeClusterSymbol;
        if (!targetSymbol) {
          return;
        }
        event.preventDefault();
        event.stopPropagation();
        openFamilyClusterProtein(targetSymbol);
        if (frozenClusterSymbol) {
          window.setTimeout(syncFrozenClusterHover, 0);
        }
      };
      if (!activeClusterSymbol && familyClusterPayload.points.length) {
        setActiveFamilyClusterSymbol(familyClusterPayload.points[0].entry_name);
      } else if (activeClusterSymbol) {
        setActiveFamilyClusterSymbol(activeClusterSymbol);
      }
      if (clusterMetaEl) {
        clusterMetaEl.textContent = "Left click to freeze selection and right click to go to receptor info page.";
      }
    });
  }

  var svgBtn = document.getElementById("cvFamilyDownloadSvgBtn");
  if (svgBtn) {
    svgBtn.addEventListener("click", function() { downloadClusterImage("svg"); });
  }
  var pngBtn = document.getElementById("cvFamilyDownloadPngBtn");
  if (pngBtn) {
    pngBtn.addEventListener("click", function() { downloadClusterImage("png"); });
  }

  renderCluster();

  document.addEventListener("click", function(event) {
    var target = event && event.target ? event.target : null;
    if (!frozenClusterSymbol || !target) {
      return;
    }
    if ((clusterLabelsEl && clusterLabelsEl.contains(target)) || (plotEl && plotEl.contains(target))) {
      return;
    }
    clearFrozenFamilyClusterSelection();
  });

  if (window.jQuery) {
    window.jQuery('a[data-toggle="tab"]').on("shown.bs.tab", function(event) {
      var href = event && event.target ? event.target.getAttribute("href") : "";
      if (href === "#cv-family-cluster-tab" && plotEl && typeof Plotly !== "undefined" && plotEl.data) {
        window.setTimeout(function() { Plotly.Plots.resize(plotEl); }, 0);
      }
    });
  }
})();

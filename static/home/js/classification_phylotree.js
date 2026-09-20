(function() {
  var baseScriptSrc = window.classificationPhyloTreeBaseSrc || "/static/home/js/phylotree.js";
  var radialFitScale = typeof window.classificationPhyloTreeRadialFitScale === "number" ? window.classificationPhyloTreeRadialFitScale : 0.84;
  var horizontalFitScale = typeof window.classificationPhyloTreeHorizontalFitScale === "number" ? window.classificationPhyloTreeHorizontalFitScale : 0.9;
  var showInternalNodes = window.classificationPhyloTreeShowInternalNodes === true;
  var showDistanceRings = window.classificationPhyloTreeShowDistanceRings !== false;
  var showRootCap = window.classificationPhyloTreeShowRootCap !== false;
  var defaultFillBySymbol = {
    "A": "#1f78b4",
    "B1": "#33a02c",
    "B2": "#6A3D9A",
    "C": "#d62728",
    "F": "#FF7F0E",
    "O1": "#17becf",
    "O2": "#bc80bd",
    "T2": "#F7B6D2"
  };
  var entityKind = String(window.classificationPhyloTreeEntityKind || "class").toLowerCase();
  var panelTitle = String(window.classificationPhyloTreePanelTitle || "Distances");
  var outerLabelToggleEnabled = window.classificationPhyloTreeEnableOuterLabelToggle === true;
  var outerLabelToggleId = String(window.classificationPhyloTreeOuterToggleId || "classificationPhyloTreeOuterToggle");
  var outerLabelsActive = window.classificationPhyloTreeOuterLabelsDefault === true;
  var outerToggleLabel = String(window.classificationPhyloTreeOuterToggleLabel || "Layout: Outer");
  var innerToggleLabel = String(window.classificationPhyloTreeInnerToggleLabel || "Layout: Edge");
  var outerLabelViewportPadding = typeof window.classificationPhyloTreeOuterViewportPadding === "number"
    ? window.classificationPhyloTreeOuterViewportPadding
    : 18;
  var autoSelectLayoutEnabled = window.classificationPhyloTreeAutoSelectLayout === true;
  var overlapStatusId = String(window.classificationPhyloTreeOverlapStatusId || "");
  var leafLabelType = String(window.classificationPhyloTreeLeafLabelType || "UniProt");
  var emptyMessage = String(
    window.classificationPhyloTreeEmptyMessage
    || (entityKind === "receptor"
      ? "Hover a receptor pill to inspect its distances and similarities to the other receptors in this family."
      : "Hover a class pill to inspect its distances and similarities to the other GPCR classes.")
  );
  var statusText = String(
    window.classificationPhyloTreeStatusText
    || (entityKind === "receptor"
      ? "Rendered with the phylogenetic tree renderer using receptor-family distances."
      : "Rendered with the phylogenetic tree renderer using classification superfamily distances.")
  );
  var treeEntities = Array.isArray(window.data && window.data.entities)
    ? window.data.entities
    : (Array.isArray(window.data && window.data.classes) ? window.data.classes : []);
  var treeMatrix = Array.isArray(window.data && window.data.matrix) ? window.data.matrix : [];
  var treeEntityIndexBySymbol = {};
  var activeTreeSymbol = null;
  var readyCallbacks = [];
  var originalJQueryReady = null;
  var baseScriptLoaded = false;
  var baseScriptStarted = false;
  var readyExecuted = false;
  var renderSettled = false;
  var pillHeight = 20;
  var pillPaddingX = 8;
  var pillMinWidth = 24;
  var pillTextFont = "600 10px Arial";
  var pillMeasurementCanvas = null;
  var lastEdgeOverlapResult = null;
  var autoLayoutApplied = false;

  function notifyFrame(eventName) {
    try {
      if (window.frameElement) {
        window.frameElement.dispatchEvent(new Event(eventName));
      }
    } catch (err) {}
  }

  function updateStatus() {
    var statusEl = document.getElementById("superfamilyTreeStatus");
    if (statusEl) {
      statusEl.textContent = statusText;
    }
  }

  function markTreeReady() {
    var pageEl = document.getElementById("sftree-page");
    if (pageEl) {
      pageEl.classList.remove("sftree-pending");
    }
  }

  var sideTitleEl = document.querySelector(".sftree-side-title");
  if (sideTitleEl) {
    sideTitleEl.textContent = panelTitle;
  }

  treeEntities.forEach(function(entityRow, index) {
    var symbol = entityRow && entityRow.symbol ? String(entityRow.symbol) : "";
    if (symbol) {
      treeEntityIndexBySymbol[symbol] = index;
    }
  });

  function ensureOverrideStyles() {
    if (document.getElementById("classification-phylotree-overrides")) {
      return;
    }

    var styleEl = document.createElement("style");
    styleEl.id = "classification-phylotree-overrides";
    styleEl.textContent = [
      "#clustering-tree path.branch { fill: none !important; }",
      "#clustering-tree g.internal-node circle { " + (showInternalNodes ? "" : "display: none !important;") + " }",
      "#clustering-tree .classification-distance-rings circle { fill: none; stroke: #6b7280; stroke-width: 3px; stroke-dasharray: 1.5 5; stroke-linecap: round; opacity: 0.2; shape-rendering: geometricPrecision; vector-effect: non-scaling-stroke; }",
      "#clustering-tree .classification-distance-rings text { fill: #8aa0b7; font-size: 9px; font-family: Arial, sans-serif; }",
      "#clustering-tree path.branch { stroke-linecap: round; stroke-linejoin: round; shape-rendering: geometricPrecision; vector-effect: non-scaling-stroke; }",
      "#clustering-tree .classification-root-cap { stroke: none; shape-rendering: geometricPrecision; vector-effect: non-scaling-stroke; }",
      "#clustering-tree .classification-pill-stub { fill: none; stroke: #94a3b8; stroke-width: 1.5px; stroke-dasharray: 2 4; stroke-linecap: round; shape-rendering: geometricPrecision; vector-effect: non-scaling-stroke; opacity: 0.95; }",
      "#clustering-tree .classification-pill rect { stroke: rgba(17, 24, 39, 0.9); stroke-width: 3px; shape-rendering: geometricPrecision; vector-effect: non-scaling-stroke; paint-order: stroke fill; }",
      "#clustering-tree .classification-pill text { fill: #ffffff; font-size: 10px; font-family: Arial, sans-serif; font-weight: 600; dominant-baseline: middle; }"
    ].join("\n");
    document.head.appendChild(styleEl);
  }

  function treeContainerVisible() {
    var treeContainer = document.getElementById("tree-container");
    return !!(treeContainer && treeContainer.offsetWidth > 0 && treeContainer.offsetHeight > 0);
  }

  function isFinitePositive(value) {
    return typeof value === "number" && isFinite(value) && value > 0;
  }

  function elementHasUsableBox(element) {
    var rect;
    if (!element) {
      return false;
    }
    rect = element.getBoundingClientRect();
    return isFinitePositive(rect.width) && isFinitePositive(rect.height);
  }

  function getZoomController() {
    return window.zoomCluster && window.zoomCluster["#clustering-tree"]
      ? window.zoomCluster["#clustering-tree"]
      : null;
  }

  function treeZoomGeometryReady() {
    return treeContainerVisible()
      && elementHasUsableBox(document.getElementById("tree-container"))
      && elementHasUsableBox(document.getElementById("clustering-tree"));
  }

  function recoverZoomController() {
    var container = "#clustering-tree";
    var currentZoom = getZoomController();

    if (!treeZoomGeometryReady() || typeof window.svgPanZoom !== "function") {
      return null;
    }

    try {
      if (currentZoom && typeof currentZoom.destroy === "function") {
        currentZoom.destroy();
      }
    } catch (err) {}

    try {
      window.zoomCluster = window.zoomCluster || {};
      window.zoomCluster[container] = window.svgPanZoom(container, {
        zoomEnabled: false,
        panEnabled: true,
        controlIconsEnabled: false,
        fit: true,
        center: true,
        minZoom: 0.1,
        maxZoom: 10,
        zoomScaleSensitivity: 0.25,
        dblClickZoomEnabled: false
      });
    } catch (err) {
      return null;
    }

    return getZoomController();
  }

  function applySafeTreeZoom(scaleFactor, allowRecovery) {
    var zoomController = getZoomController();
    var currentZoom;
    var nextZoom;

    if (!treeZoomGeometryReady() || !zoomController) {
      return false;
    }

    try {
      zoomController.updateBBox();
      zoomController.fit();
      currentZoom = Number(zoomController.getZoom());
      nextZoom = currentZoom * scaleFactor;
      if (isFinitePositive(nextZoom)) {
        zoomController.zoom(nextZoom);
      }
      zoomController.center();
      return true;
    } catch (err) {
      if (allowRecovery !== false && recoverZoomController()) {
        return applySafeTreeZoom(scaleFactor, false);
      }
    }

    return false;
  }

  function installSafeWindowResize() {
    window.onresize = function() {
      var treeContainer = document.getElementById("tree-container");
      var plot = document.getElementById("clustering-tree");
      var plotsize;
      var zoomController;

      if (!treeContainerVisible() || !treeContainer || !plot) {
        return;
      }

      plotsize = window.innerHeight * 0.9;
      if ((treeContainer.offsetWidth * 0.9) < plotsize) {
        plotsize = treeContainer.offsetWidth * 0.9;
      }

      if (!isFinitePositive(plotsize)) {
        return;
      }

      plot.style.height = plotsize + "px";
      plot.style.width = plotsize + "px";

      if (!treeZoomGeometryReady()) {
        return;
      }

      zoomController = getZoomController();
      if (zoomController && typeof zoomController.resize === "function") {
        try {
          zoomController.resize();
        } catch (err) {
          recoverZoomController();
        }
      }

      if (typeof window.resizeTree === "function") {
        window.resizeTree();
      }
    };
  }

  function parseTranslate(transformValue) {
    if (!transformValue) {
      return null;
    }

    var match = /translate\s*\(\s*([-\d.eE]+)[,\s]+([-\d.eE]+)\s*\)/.exec(transformValue);
    if (!match) {
      return null;
    }

    return [parseFloat(match[1]), parseFloat(match[2])];
  }

  function getPillMeasurementContext() {
    if (!pillMeasurementCanvas) {
      pillMeasurementCanvas = document.createElement("canvas");
    }
    var context = pillMeasurementCanvas.getContext("2d");
    if (!context) {
      return null;
    }
    context.font = pillTextFont;
    return context;
  }

  function estimatePillWidth(label) {
    var pillLabel = plainLabelText(label).trim();
    if (!pillLabel) {
      return pillMinWidth;
    }
    var context = getPillMeasurementContext();
    if (!context) {
      return Math.max(pillMinWidth, (pillLabel.length * 7) + (pillPaddingX * 2));
    }
    return Math.max(pillMinWidth, Math.ceil(context.measureText(pillLabel).width + (pillPaddingX * 2)));
  }

  function entityLeafLabel(entity, fallbackLabel) {
    if (!entity) {
      return String(fallbackLabel || "");
    }
    if (leafLabelType === "Protein") {
      return String(entity.name || entity.short_label || entity.symbol || fallbackLabel || "");
    }
    if (leafLabelType === "Gene") {
      return String(entity.gene_label || entity.short_label || entity.symbol || fallbackLabel || "");
    }
    return String(entity.short_label || entity.symbol || fallbackLabel || "");
  }

  function resolveLeafPillData(leafNode, textNode) {
    var pillLabel = String(textNode && textNode.textContent || "").trim();
    var leafDatum = leafNode && leafNode.__data__ ? leafNode.__data__ : {};
    var leafName = leafDatum.name || pillLabel;
    var pillMarker = null;
    var entity = null;
    var pillDisplay = pillLabel;

    Array.prototype.forEach.call(leafNode.querySelectorAll("circle"), function(circleNode) {
      if (!pillMarker && circleNode.getAttribute("data-datalabel") !== null) {
        pillMarker = circleNode;
      }
    });

    entity = treeEntityBySymbol(leafName) || treeEntityBySymbol(pillLabel);
    if (entity) {
      pillDisplay = entityLeafLabel(entity, pillLabel);
    }

    return {
      entity: entity,
      leafName: leafName,
      pillLabel: pillLabel,
      pillMarker: pillMarker,
      pillDisplay: pillDisplay,
      pillWidth: estimatePillWidth(pillDisplay)
    };
  }

  function collectLeafPillGeometry() {
    var svg = document.getElementById("clustering-tree");
    var rootNode = svg ? svg.querySelector("g.internal-node") : null;
    var leafNodes = svg ? svg.querySelectorAll("g.node") : [];
    var geometry = [];
    var rootCoords;

    if (!rootNode || !leafNodes.length) {
      return null;
    }

    rootCoords = parseTranslate(rootNode.getAttribute("transform"));
    if (!rootCoords) {
      return null;
    }

    Array.prototype.forEach.call(leafNodes, function(leafNode) {
      var leafCoords = parseTranslate(leafNode.getAttribute("transform"));
      var textNode = leafNode.querySelector("text");
      var dx;
      var dy;
      var leafRadius;
      var pillData;
      var radialX;
      var radialY;

      if (!leafCoords || !textNode) {
        return;
      }

      dx = leafCoords[0] - rootCoords[0];
      dy = leafCoords[1] - rootCoords[1];
      leafRadius = Math.sqrt(dx * dx + dy * dy);
      if (!isFinite(leafRadius) || leafRadius <= 0) {
        return;
      }

      pillData = resolveLeafPillData(leafNode, textNode);
      if (!pillData.pillLabel) {
        return;
      }

      radialX = dx / leafRadius;
      radialY = dy / leafRadius;

      geometry.push({
        angle: Math.atan2(dy, dx),
        centerX: leafCoords[0] + (radialX * (pillData.pillWidth / 2)),
        centerY: leafCoords[1] + (radialY * (pillData.pillWidth / 2)),
        halfHeight: pillHeight / 2,
        halfWidth: pillData.pillWidth / 2,
        label: pillData.pillDisplay,
        leafRadius: leafRadius,
        radialX: radialX,
        radialY: radialY,
        symbol: pillData.leafName || pillData.pillLabel
      });
    });

    return geometry;
  }

  function buildPillCorners(geometryItem) {
    var tangentX = -geometryItem.radialY;
    var tangentY = geometryItem.radialX;
    var halfWidth = geometryItem.halfWidth;
    var halfHeight = geometryItem.halfHeight;
    var centerX = geometryItem.centerX;
    var centerY = geometryItem.centerY;

    return [
      {
        x: centerX - (geometryItem.radialX * halfWidth) - (tangentX * halfHeight),
        y: centerY - (geometryItem.radialY * halfWidth) - (tangentY * halfHeight)
      },
      {
        x: centerX + (geometryItem.radialX * halfWidth) - (tangentX * halfHeight),
        y: centerY + (geometryItem.radialY * halfWidth) - (tangentY * halfHeight)
      },
      {
        x: centerX + (geometryItem.radialX * halfWidth) + (tangentX * halfHeight),
        y: centerY + (geometryItem.radialY * halfWidth) + (tangentY * halfHeight)
      },
      {
        x: centerX - (geometryItem.radialX * halfWidth) + (tangentX * halfHeight),
        y: centerY - (geometryItem.radialY * halfWidth) + (tangentY * halfHeight)
      }
    ];
  }

  function projectCorners(corners, axisX, axisY) {
    var minProjection = Infinity;
    var maxProjection = -Infinity;

    corners.forEach(function(corner) {
      var projection = (corner.x * axisX) + (corner.y * axisY);
      minProjection = Math.min(minProjection, projection);
      maxProjection = Math.max(maxProjection, projection);
    });

    return {
      min: minProjection,
      max: maxProjection
    };
  }

  function rectanglesIntersect(cornersA, cornersB) {
    var polygons = [cornersA, cornersB];
    var polygonIndex;
    var cornerIndex;

    for (polygonIndex = 0; polygonIndex < polygons.length; polygonIndex += 1) {
      var polygon = polygons[polygonIndex];
      for (cornerIndex = 0; cornerIndex < polygon.length; cornerIndex += 1) {
        var nextIndex = (cornerIndex + 1) % polygon.length;
        var edgeX = polygon[nextIndex].x - polygon[cornerIndex].x;
        var edgeY = polygon[nextIndex].y - polygon[cornerIndex].y;
        var axisX = -edgeY;
        var axisY = edgeX;
        var axisLength = Math.sqrt((axisX * axisX) + (axisY * axisY));
        var projectionA;
        var projectionB;

        if (!axisLength) {
          continue;
        }

        axisX /= axisLength;
        axisY /= axisLength;
        projectionA = projectCorners(cornersA, axisX, axisY);
        projectionB = projectCorners(cornersB, axisX, axisY);

        if (projectionA.max < projectionB.min || projectionB.max < projectionA.min) {
          return false;
        }
      }
    }

    return true;
  }

  function calculateEdgeOverlap() {
    var geometry = collectLeafPillGeometry();
    var index;

    if (!geometry || geometry.length < 2) {
      return null;
    }

    geometry.sort(function(a, b) {
      return a.angle - b.angle;
    });

    for (index = 0; index < geometry.length; index += 1) {
      var current = geometry[index];
      var next = geometry[(index + 1) % geometry.length];
      var dx = current.centerX - next.centerX;
      var dy = current.centerY - next.centerY;
      var halfDiagonalCurrent = Math.sqrt((current.halfWidth * current.halfWidth) + (current.halfHeight * current.halfHeight));
      var halfDiagonalNext = Math.sqrt((next.halfWidth * next.halfWidth) + (next.halfHeight * next.halfHeight));
      var centerDistance = Math.sqrt((dx * dx) + (dy * dy));

      if (centerDistance <= (halfDiagonalCurrent + halfDiagonalNext)
        && rectanglesIntersect(buildPillCorners(current), buildPillCorners(next))) {
        return {
          hasOverlap: true,
          pairCount: 1,
          samplePairs: [current.label + "/" + next.label]
        };
      }
    }

    return {
      hasOverlap: false,
      pairCount: 0,
      samplePairs: []
    };
  }

  function updateOverlapStatus(mode) {
    var statusEl = overlapStatusId ? document.getElementById(overlapStatusId) : null;
    var currentLayout = outerLabelsActive ? "Default" : "Edge";
    var overlapText = "pending";

    if (!statusEl) {
      return;
    }

    if (!lastEdgeOverlapResult) {
      lastEdgeOverlapResult = calculateEdgeOverlap();
    }

    if (lastEdgeOverlapResult) {
      overlapText = lastEdgeOverlapResult.hasOverlap
        ? "yes (" + String(lastEdgeOverlapResult.pairCount) + " pairs)"
        : "no";
    }

    statusEl.classList.remove("text-danger", "text-success");
    statusEl.classList.add("text-muted");
    if (lastEdgeOverlapResult) {
      statusEl.classList.toggle("text-danger", lastEdgeOverlapResult.hasOverlap);
      statusEl.classList.toggle("text-success", !lastEdgeOverlapResult.hasOverlap);
      statusEl.classList.remove("text-muted");
    }

    statusEl.textContent = (mode === "auto" ? "Layout auto-selected: " : "Current layout: ")
      + currentLayout
      + " | Edge overlap detected: "
      + overlapText;
  }

  function primeOverlapStatus() {
    var statusEl = overlapStatusId ? document.getElementById(overlapStatusId) : null;
    if (!statusEl || String(statusEl.textContent || "").trim()) {
      return;
    }
    statusEl.textContent = "Checking edge overlap...";
  }

  function applyAutomaticLayoutChoice() {
    if (!autoSelectLayoutEnabled || autoLayoutApplied) {
      return;
    }

    lastEdgeOverlapResult = calculateEdgeOverlap();
    if (!lastEdgeOverlapResult) {
      return;
    }

    outerLabelsActive = lastEdgeOverlapResult.hasOverlap;
    autoLayoutApplied = true;
  }

  function formatDistanceLabel(distanceValue) {
    if (!isFinite(distanceValue)) {
      return "";
    }
    return String(Math.round(distanceValue));
  }

  function formatMetricValue(value) {
    if (value === null || value === undefined || value === "") {
      return "n/a";
    }
    var rounded = Math.round(Number(value) * 10) / 10;
    if (!isFinite(rounded)) {
      return "n/a";
    }
    return Math.round(rounded) === rounded ? String(Math.round(rounded)) : rounded.toFixed(1);
  }

  function escapeHtml(value) {
    return String(value == null ? "" : value)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;")
      .replace(/'/g, "&#39;");
  }

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

  function plainLabelText(value) {
    var wrapper = document.createElement("span");
    wrapper.innerHTML = formattedLabelHtml(value);
    return wrapper.textContent || wrapper.innerText || "";
  }

  function formattedSvgLabelHtml(value) {
    return escapeHtml(decodeHtmlEntities(value))
      .replace(/&lt;sub&gt;/gi, '<tspan baseline-shift="-0.25em" font-size="70%">')
      .replace(/&lt;\/sub&gt;/gi, "</tspan>")
      .replace(/&lt;i&gt;/gi, '<tspan font-style="italic">')
      .replace(/&lt;\/i&gt;/gi, "</tspan>");
  }

  function normalizeTreeClassName(name) {
    return String(name || "")
      .replace("Class O1 (fish-like)", "Class O1 (fish-like olfactory receptors)")
      .replace("Class O2 (tetrapod specific)", "Class O2 (tetrapod-specific olfactory receptors)");
  }

  function sentenceCaseFragment(value) {
    var raw = String(value || "").trim().toLowerCase();
    return raw.replace(/([A-Za-z])/, function(match) {
      return match.toUpperCase();
    });
  }

  function formattedClassNameHtml(name) {
    var raw = normalizeTreeClassName(name).trim();
    if (!raw) {
      return "";
    }
    var match = raw.match(/^(.*?)(\s*\(([^)]+)\))$/);
    if (!match) {
      return escapeHtml(raw);
    }
    var mainText = match[1].trim();
    if (/^Class\s+/i.test(mainText)) {
      mainText = mainText.replace(/^Class\s+/i, "");
    }
    var main = "<strong>" + escapeHtml("Class " + mainText) + "</strong>";
    var secondary = escapeHtml("(" + sentenceCaseFragment(match[3].trim()) + ")");
    return main + '<span class="sftree-secondary">' + secondary + "</span>";
  }

  function formattedEntityNameHtml(entity) {
    if (!entity) {
      return "";
    }
    var name = String(entity.name || entity.symbol || "").trim();
    if (!name) {
      return "";
    }
    if (entityKind === "class") {
      return formattedClassNameHtml(name);
    }
    var main = "<strong>" + formattedLabelHtml(name) + "</strong>";
    var secondaryParts = [];
    if (entity.short_label) {
      secondaryParts.push(String(entity.short_label));
    }
    if (entity.subtitle) {
      secondaryParts.push(String(entity.subtitle));
    }
    if (!secondaryParts.length) {
      return main;
    }
    return main + '<span class="sftree-secondary">' + formattedLabelHtml(secondaryParts.join(" | ")) + "</span>";
  }

  function treeEntityBySymbol(symbol) {
    var classIndex = treeEntityIndexBySymbol[String(symbol || "")];
    if (classIndex === undefined) {
      return null;
    }
    return treeEntities[classIndex] || null;
  }

  function treeMatrixEntry(sourceSymbol, targetSymbol) {
    var sourceIndex = treeEntityIndexBySymbol[String(sourceSymbol || "")];
    var targetIndex = treeEntityIndexBySymbol[String(targetSymbol || "")];
    if (sourceIndex === undefined || targetIndex === undefined) {
      return null;
    }
    if (!treeMatrix[sourceIndex] || !treeMatrix[sourceIndex][targetIndex]) {
      return null;
    }
    return treeMatrix[sourceIndex][targetIndex];
  }

  function renderTreeHoverCard(symbol) {
    var hoverEl = document.getElementById("sftree-hover");
    var activeEntity = treeEntityBySymbol(symbol);
    if (!hoverEl) {
      return;
    }
    if (!activeEntity) {
      hoverEl.innerHTML = '<p class="sftree-hover-empty">' + escapeHtml(emptyMessage) + '</p>';
      return;
    }

    var peers = [];
    treeEntities.forEach(function(otherEntity) {
      var otherSymbol = otherEntity && otherEntity.symbol ? String(otherEntity.symbol) : "";
      if (!otherSymbol || otherSymbol === activeEntity.symbol) {
        return;
      }
      var cell = treeMatrixEntry(activeEntity.symbol, otherSymbol) || {};
      peers.push({
        symbol: otherSymbol,
        label: otherEntity.short_label || otherSymbol,
        color: otherEntity.color || defaultFillBySymbol[otherSymbol] || "#808080",
        distance: cell.distance,
        similarityDisplay: cell.similarity_display || "n/a"
      });
    });

    peers.sort(function(a, b) {
      var distanceA = a.distance === null || a.distance === undefined ? Infinity : Number(a.distance);
      var distanceB = b.distance === null || b.distance === undefined ? Infinity : Number(b.distance);
      return distanceA - distanceB;
    });

    var headerHtml = (
      '<div class="sftree-hover-header">' +
        '<div class="sftree-hover-header-name">' + formattedEntityNameHtml(activeEntity) + '</div>' +
      '</div>'
    );
    var rowsHtml = peers.map(function(peer) {
      return (
        '<div class="sftree-hover-row">' +
          '<div class="sftree-hover-label">' +
            '<span class="sftree-hover-swatch" style="background:' + escapeHtml(peer.color) + ';"></span>' +
            '<div class="sftree-hover-text">' + escapeHtml(peer.label) + '</div>' +
          '</div>' +
          '<div class="sftree-hover-metrics">' +
            '<div class="sftree-hover-metric-label">Distance</div>' +
            '<div class="sftree-hover-distance">' + escapeHtml(formatMetricValue(peer.distance)) + '</div>' +
          '</div>' +
          '<div class="sftree-hover-metrics is-similarity">' +
            '<div class="sftree-hover-metric-label">Similarity</div>' +
            '<div class="sftree-hover-similarity">' + escapeHtml(peer.similarityDisplay) + '</div>' +
          '</div>' +
        '</div>'
      );
    }).join("");
    hoverEl.innerHTML = headerHtml + '<div class="sftree-hover-list">' + rowsHtml + '</div>';
  }

  function updateActivePillState() {
    var svg = document.getElementById("clustering-tree");
    var pillGroups = svg ? svg.querySelectorAll(".classification-pill") : [];
    Array.prototype.forEach.call(pillGroups, function(pillGroup) {
      pillGroup.classList.toggle("is-active", pillGroup.getAttribute("data-class-symbol") === activeTreeSymbol);
    });
  }

  function setActiveTreeSymbol(symbol) {
    if (!treeEntityBySymbol(symbol)) {
      return;
    }
    activeTreeSymbol = String(symbol);
    updateActivePillState();
    renderTreeHoverCard(activeTreeSymbol);
    notifyFrame("cv:content-resize");
  }

  window.classificationPhyloTreeSetLeafLabelType = function(labelType) {
    leafLabelType = String(labelType || "UniProt");
    window.classificationPhyloTreeLeafLabelType = leafLabelType;
    renderLeafPills();
    updateActivePillState();
    updateOverlapStatus(autoLayoutApplied ? "auto" : "manual");
    if (typeof window.resizeTree === "function") {
      window.resizeTree();
      window.setTimeout(window.resizeTree, 200);
    }
    notifyFrame("cv:content-resize");
  };

  function attachTreeHoverInteractions() {
    var svg = document.getElementById("clustering-tree");
    var pillGroups = svg ? svg.querySelectorAll(".classification-pill") : [];
    var firstSymbol = null;
    Array.prototype.forEach.call(pillGroups, function(pillGroup) {
      var symbol = pillGroup.getAttribute("data-class-symbol") || "";
      if (!symbol) {
        return;
      }
      if (!firstSymbol) {
        firstSymbol = symbol;
      }
      pillGroup.setAttribute("tabindex", "0");
      pillGroup.setAttribute("role", "button");
      pillGroup.setAttribute("aria-label", symbol + " " + entityKind + " distances");
      pillGroup.onmouseenter = function() {
        setActiveTreeSymbol(symbol);
      };
      pillGroup.onfocus = function() {
        setActiveTreeSymbol(symbol);
      };
      pillGroup.onclick = function() {
        setActiveTreeSymbol(symbol);
      };
    });

    if (activeTreeSymbol && treeEntityBySymbol(activeTreeSymbol)) {
      setActiveTreeSymbol(activeTreeSymbol);
    } else if (firstSymbol) {
      setActiveTreeSymbol(firstSymbol);
    } else if (treeEntities.length) {
      setActiveTreeSymbol(treeEntities[0].symbol);
    }
  }

  function parseNewickTree(newickText) {
    if (!newickText || typeof newickText !== "string") {
      return null;
    }

    var index = 0;

    function skipWhitespace() {
      while (index < newickText.length && /\s/.test(newickText.charAt(index))) {
        index += 1;
      }
    }

    function readToken() {
      skipWhitespace();
      var start = index;
      while (index < newickText.length && !/[,:();]/.test(newickText.charAt(index))) {
        index += 1;
      }
      return newickText.slice(start, index).trim();
    }

    function readLength() {
      skipWhitespace();
      if (newickText.charAt(index) !== ":") {
        return 0;
      }

      index += 1;
      skipWhitespace();
      var start = index;
      while (index < newickText.length && !/[,);]/.test(newickText.charAt(index))) {
        index += 1;
      }

      var parsed = parseFloat(newickText.slice(start, index).trim());
      return isFinite(parsed) ? parsed : 0;
    }

    function parseNode() {
      skipWhitespace();
      var node = { name: "", length: 0, children: [] };

      if (newickText.charAt(index) === "(") {
        index += 1;
        while (index < newickText.length) {
          node.children.push(parseNode());
          skipWhitespace();
          if (newickText.charAt(index) === ",") {
            index += 1;
            continue;
          }
          if (newickText.charAt(index) === ")") {
            index += 1;
            break;
          }
          break;
        }
        node.name = readToken();
        node.length = readLength();
        return node;
      }

      node.name = readToken();
      node.length = readLength();
      return node;
    }

    try {
      return parseNode();
    } catch (err) {
      return null;
    }
  }

  function getMaxTreeDistance(treeNode, runningDistance) {
    if (!treeNode) {
      return 0;
    }

    var accumulated = (runningDistance || 0) + (isFinite(treeNode.length) ? treeNode.length : 0);
    if (!treeNode.children || !treeNode.children.length) {
      return accumulated;
    }

    var childDistances = treeNode.children.map(function(childNode) {
      return getMaxTreeDistance(childNode, accumulated);
    });
    return Math.max.apply(Math, childDistances);
  }

  function renderLeafPills() {
    var svg = document.getElementById("clustering-tree");
    var rootNode = svg ? svg.querySelector("g.internal-node") : null;
    var leafNodes = svg ? svg.querySelectorAll("g.node") : [];
    if (!rootNode || !leafNodes.length) {
      return;
    }

    var rootCoords = parseTranslate(rootNode.getAttribute("transform"));
    if (!rootCoords) {
      return;
    }

    var maxLeafRadius = 0;
    Array.prototype.forEach.call(leafNodes, function(leafNode) {
      var leafCoords = parseTranslate(leafNode.getAttribute("transform"));
      if (!leafCoords) {
        return;
      }
      var leafDx = leafCoords[0] - rootCoords[0];
      var leafDy = leafCoords[1] - rootCoords[1];
      maxLeafRadius = Math.max(maxLeafRadius, Math.sqrt(leafDx * leafDx + leafDy * leafDy));
    });

    Array.prototype.forEach.call(leafNodes, function(leafNode) {
      var leafCoords = parseTranslate(leafNode.getAttribute("transform"));
      var textNode = leafNode.querySelector("text");
      if (!leafCoords || !textNode) {
        return;
      }

      var dx = leafCoords[0] - rootCoords[0];
      var dy = leafCoords[1] - rootCoords[1];
      var leafRadius = Math.sqrt(dx * dx + dy * dy);
      var computedStyle = window.getComputedStyle(textNode);
      var connectorLength = outerLabelsActive ? Math.max(0, maxLeafRadius - leafRadius) : 0;
      var pillGap = outerLabelsActive ? connectorLength : 0;
      var pillGroup = leafNode.querySelector("g.classification-pill");
      var pillContent;
      var pillText;
      var pillRect;
      var pillData = resolveLeafPillData(leafNode, textNode);
      var pillLabel = pillData.pillLabel;
      var leafName = pillData.leafName;
      var pillMarker = pillData.pillMarker;
      var entity = pillData.entity;
      var pillFill = (entity && entity.color) || defaultFillBySymbol[pillLabel] || defaultFillBySymbol[leafName] || "";
      var pillDisplay = pillData.pillDisplay;
      if (!pillFill) {
        pillFill = pillMarker ? (window.getComputedStyle(pillMarker).fill || pillMarker.getAttribute("fill") || pillMarker.style.fill) : "";
      }
      if (!pillFill) {
        pillFill = leafNode.style.fill || computedStyle.fill || "#4b7bbb";
      }
      var angleDeg = Math.atan2(dy, dx) * 180 / Math.PI;
      var flipText = dx < 0;

      if (!pillLabel) {
        return;
      }

      if (!pillGroup) {
        pillGroup = document.createElementNS("http://www.w3.org/2000/svg", "g");
        pillGroup.setAttribute("class", "classification-pill");
        pillContent = document.createElementNS("http://www.w3.org/2000/svg", "g");
        pillRect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
        pillText = document.createElementNS("http://www.w3.org/2000/svg", "text");
        pillContent.appendChild(pillRect);
        pillContent.appendChild(pillText);
        pillGroup.appendChild(pillContent);
        leafNode.appendChild(pillGroup);
      } else {
        pillContent = pillGroup.querySelector("g");
        pillRect = pillGroup.querySelector("rect");
        pillText = pillGroup.querySelector("text");
      }

      Array.prototype.forEach.call(leafNode.querySelectorAll("line.classification-pill-stub"), function(stubNode) {
        if (stubNode.parentNode) {
          stubNode.parentNode.removeChild(stubNode);
        }
      });

      var tracer = leafNode.querySelector("line");
      if (tracer) {
        tracer.style.display = "none";
      }
      Array.prototype.forEach.call(leafNode.querySelectorAll("circle"), function(circleNode) {
        circleNode.style.display = "none";
      });
      textNode.style.display = "none";

      pillText.innerHTML = formattedSvgLabelHtml(pillDisplay);
      pillText.setAttribute("text-anchor", "middle");
      pillText.setAttribute("x", 0);
      pillText.setAttribute("y", 0.5);

      var pillWidth = pillData.pillWidth;
      var pillCenter = pillGap + pillWidth / 2;

      if (outerLabelsActive && connectorLength > 0.5) {
        var stubNode = document.createElementNS("http://www.w3.org/2000/svg", "line");
        stubNode.setAttribute("class", "classification-pill-stub");
        stubNode.setAttribute("x1", 0);
        stubNode.setAttribute("y1", 0);
        stubNode.setAttribute("x2", connectorLength);
        stubNode.setAttribute("y2", 0);
        stubNode.setAttribute("transform", "rotate(" + angleDeg + ")");
        leafNode.appendChild(stubNode);
      }

      pillRect.setAttribute("width", pillWidth);
      pillRect.setAttribute("height", pillHeight);
      pillRect.setAttribute("rx", pillHeight / 2);
      pillRect.setAttribute("ry", pillHeight / 2);
      pillRect.setAttribute("fill", pillFill);
      pillRect.setAttribute("x", -pillWidth / 2);
      pillRect.setAttribute("y", -pillHeight / 2);
      pillText.setAttribute("x", 0);
      pillText.setAttribute("y", 1);

      leafNode.setAttribute("data-class-symbol", leafName || pillLabel);
      pillGroup.setAttribute("data-class-symbol", leafName || pillLabel);
      pillGroup.setAttribute("transform", "rotate(" + angleDeg + ") translate(" + pillCenter + ",0)");
      pillContent.setAttribute("transform", flipText ? "rotate(180)" : "");

      leafNode.appendChild(pillGroup);
    });
  }

  function outerLabelFitScale() {
    if (!outerLabelsActive || !window.radialTree) {
      return radialFitScale;
    }

    var svg = document.getElementById("clustering-tree");
    var rootNode = svg ? svg.querySelector("g.internal-node") : null;
    var pillRects = svg ? svg.querySelectorAll("g.classification-pill rect") : [];
    if (!svg || !rootNode || !pillRects.length || !svg.clientWidth || !svg.clientHeight) {
      return radialFitScale;
    }

    var rootCoords = parseTranslate(rootNode.getAttribute("transform"));
    if (!rootCoords) {
      return radialFitScale;
    }

    var baseRadius = 0;
    var maxExtraX = 0;
    var maxExtraY = 0;

    Array.prototype.forEach.call(pillRects, function(pillRect) {
      var pillGroup = pillRect.parentNode && pillRect.parentNode.parentNode;
      var pillWidth = parseFloat(pillRect.getAttribute("width")) || 0;
      var pillHeight = parseFloat(pillRect.getAttribute("height")) || 0;
      var pillCoords = pillGroup ? parseTranslate(pillGroup.getAttribute("transform")) : null;
      if (!pillCoords) {
        return;
      }

      var dx = pillCoords[0] - rootCoords[0];
      var dy = pillCoords[1] - rootCoords[1];
      var centerRadius = Math.sqrt(dx * dx + dy * dy);
      var angle = Math.atan2(dy, dx);
      var projectedHalfX = Math.abs(Math.cos(angle)) * (pillWidth / 2) + Math.abs(Math.sin(angle)) * (pillHeight / 2);
      var projectedHalfY = Math.abs(Math.sin(angle)) * (pillWidth / 2) + Math.abs(Math.cos(angle)) * (pillHeight / 2);

      baseRadius = Math.max(baseRadius, Math.max(0, centerRadius - (pillWidth / 2)));
      maxExtraX = Math.max(maxExtraX, projectedHalfX);
      maxExtraY = Math.max(maxExtraY, projectedHalfY);
    });

    if (!isFinite(baseRadius) || baseRadius <= 0) {
      return radialFitScale;
    }

    var availableHalfWidth = Math.max(1, (svg.clientWidth / 2) - outerLabelViewportPadding);
    var availableHalfHeight = Math.max(1, (svg.clientHeight / 2) - outerLabelViewportPadding);
    var widthScale = availableHalfWidth / (baseRadius + maxExtraX);
    var heightScale = availableHalfHeight / (baseRadius + maxExtraY);
    var fittedScale = radialFitScale * Math.min(1, widthScale, heightScale);

    return Math.max(0.55, fittedScale);
  }

  function syncOuterLabelToggle() {
    if (!outerLabelToggleEnabled) {
      return;
    }
    var toggleButton = document.getElementById(outerLabelToggleId);
    if (!toggleButton) {
      return;
    }
    toggleButton.classList.toggle("active", outerLabelsActive);
    toggleButton.setAttribute("aria-pressed", outerLabelsActive ? "true" : "false");
    toggleButton.textContent = outerLabelsActive ? outerToggleLabel : innerToggleLabel;
  }

  function rerenderLeafPills() {
    lastEdgeOverlapResult = calculateEdgeOverlap();
    renderLeafPills();
    attachTreeHoverInteractions();
    updateActivePillState();
    updateOverlapStatus("manual");
    if (typeof window.resizeTree === "function") {
      window.resizeTree();
      window.setTimeout(window.resizeTree, 200);
    }
    notifyFrame("cv:content-resize");
  }

  function mountOuterLabelToggle() {
    if (!outerLabelToggleEnabled) {
      return;
    }
    var toggleButton = document.getElementById(outerLabelToggleId);
    if (!toggleButton || toggleButton.getAttribute("data-tree-toggle-bound") === "1") {
      syncOuterLabelToggle();
      return;
    }
    toggleButton.setAttribute("data-tree-toggle-bound", "1");
    syncOuterLabelToggle();
    toggleButton.addEventListener("click", function() {
      outerLabelsActive = !outerLabelsActive;
      autoLayoutApplied = false;
      syncOuterLabelToggle();
      rerenderLeafPills();
    });
  }

  function renderDistanceRings() {
    if (!showDistanceRings || !window.radialTree) {
      return;
    }

    var svg = document.getElementById("clustering-tree");
    var container = svg ? svg.querySelector(".phylotree-container") : null;
    var rootNode = svg ? svg.querySelector("g.internal-node") : null;
    var leafNodes = svg ? svg.querySelectorAll("g.node") : [];
    if (!container || !rootNode || !leafNodes.length) {
      return;
    }

    var existingLayer = container.querySelector(".classification-distance-rings");
    if (existingLayer) {
      existingLayer.parentNode.removeChild(existingLayer);
    }

    var rootCoords = parseTranslate(rootNode.getAttribute("transform"));
    if (!rootCoords) {
      return;
    }

    var maxRadius = 0;
    Array.prototype.forEach.call(leafNodes, function(leafNode) {
      var leafCoords = parseTranslate(leafNode.getAttribute("transform"));
      if (!leafCoords) {
        return;
      }

      var dx = leafCoords[0] - rootCoords[0];
      var dy = leafCoords[1] - rootCoords[1];
      maxRadius = Math.max(maxRadius, Math.sqrt(dx * dx + dy * dy));
    });

    if (!isFinite(maxRadius) || maxRadius <= 0) {
      return;
    }

    var parsedTree = parseNewickTree(window.data && window.data.tree);
    var maxDistance = getMaxTreeDistance(parsedTree, 0);
    var roundedMaxDistance = Math.max(1, Math.round(maxDistance));
    var svgNs = "http://www.w3.org/2000/svg";
    var ringLayer = document.createElementNS(svgNs, "g");
    ringLayer.setAttribute("class", "classification-distance-rings");
    container.insertBefore(ringLayer, container.firstChild);

    var ringCount = 4;
    var renderedRingValues = {};
    for (var ringIndex = 1; ringIndex <= ringCount; ringIndex += 1) {
      var ratio = ringIndex / ringCount;
      var distanceValue = Math.max(1, Math.round(roundedMaxDistance * ratio));
      if (renderedRingValues[distanceValue]) {
        continue;
      }
      renderedRingValues[distanceValue] = true;

      var radius = maxRadius * (distanceValue / roundedMaxDistance);

      var circle = document.createElementNS(svgNs, "circle");
      circle.setAttribute("cx", rootCoords[0]);
      circle.setAttribute("cy", rootCoords[1]);
      circle.setAttribute("r", radius);
      ringLayer.appendChild(circle);

      var label = document.createElementNS(svgNs, "text");
      label.setAttribute("x", rootCoords[0]);
      label.setAttribute("y", rootCoords[1] - radius - 6);
      label.setAttribute("text-anchor", "middle");
      label.textContent = formatDistanceLabel(distanceValue);
      ringLayer.appendChild(label);
    }
  }

  function renderBranchDecorations() {
    var svg = document.getElementById("clustering-tree");
    var container = svg ? svg.querySelector(".phylotree-container") : null;
    var rootNode = svg ? svg.querySelector("g.internal-node") : null;
    var branchNodes = container ? container.querySelectorAll("path.branch") : [];
    if (!container || !rootNode || !branchNodes.length) {
      return;
    }

    var rootCoords = parseTranslate(rootNode.getAttribute("transform"));
    if (!rootCoords) {
      return;
    }

    var existingRootCap = container.querySelector(".classification-root-cap");
    if (existingRootCap) {
      existingRootCap.parentNode.removeChild(existingRootCap);
    }

    var svgNs = "http://www.w3.org/2000/svg";
    var branchColor = "#111827";

    Array.prototype.forEach.call(branchNodes, function(branchNode) {
      branchNode.setAttribute("shape-rendering", "geometricPrecision");
      branchNode.setAttribute("vector-effect", "non-scaling-stroke");
      branchNode.style.shapeRendering = "geometricPrecision";
      branchNode.style.stroke = branchColor;
      branchNode.setAttribute("stroke", branchColor);
    });

    if (!showRootCap) {
      return;
    }

    var rootCap = document.createElementNS(svgNs, "circle");
    rootCap.setAttribute("class", "classification-root-cap");
    rootCap.setAttribute("cx", rootCoords[0]);
    rootCap.setAttribute("cy", rootCoords[1]);
    rootCap.setAttribute("r", 2.125);
    rootCap.setAttribute("fill", branchColor);
    container.appendChild(rootCap);
  }

  function installClassificationPatches() {
    if (window.__classificationPhyloTreePatched) {
      return;
    }

    var originalResizeTree = window.resizeTree;
    var originalNodeStyler = window.nodeStyler;
    var originalBranchStyler = window.branchStyler;

    window.resizeTree = async function(eval_string) {
      var expression = typeof eval_string === "string" ? eval_string : "";
      var scaleFactor;

      if (radialTree) {
        await sleep(100);
      }

      scaleFactor = radialTree ? outerLabelFitScale() : horizontalFitScale;
      if (!isFinitePositive(scaleFactor)) {
        scaleFactor = radialTree ? radialFitScale : horizontalFitScale;
      }
      applySafeTreeZoom(scaleFactor, true);

      if (expression.length > 0) {
        await sleep(100);
        eval(expression);
      }
    };

    window.nodeStyler = function(element, node) {
      if (!showInternalNodes && !d3.layout.phylotree.is_leafnode(node)) {
        element.selectAll("circle")
          .style("fill", "none")
          .style("stroke-width", "0px")
          .attr("r", 0);
        return;
      }

      return originalNodeStyler(element, node);
    };

    window.branchStyler = function(element, node) {
      element.style("fill", "none");
      return originalBranchStyler(element, node);
    };

    window.__classificationPhyloTreePatched = {
      originalResizeTree: originalResizeTree
    };
  }

  function interceptJQueryReady() {
    if (!window.jQuery || !window.jQuery.fn || typeof window.jQuery.fn.ready !== "function") {
      return false;
    }

    originalJQueryReady = window.jQuery.fn.ready;
    window.jQuery.fn.ready = function(callback) {
      readyCallbacks.push(callback);
      return this;
    };
    return true;
  }

  function restoreJQueryReady() {
    if (originalJQueryReady && window.jQuery && window.jQuery.fn) {
      window.jQuery.fn.ready = originalJQueryReady;
    }
  }

  function finalizeRender() {
    if (renderSettled) {
      return;
    }

    if (!window.zoomCluster || !window.zoomCluster["#clustering-tree"]) {
      window.setTimeout(finalizeRender, 100);
      return;
    }

    renderSettled = true;
    ensureOverrideStyles();
    updateStatus();
    renderBranchDecorations();
    renderDistanceRings();
    applyAutomaticLayoutChoice();
    renderLeafPills();
    attachTreeHoverInteractions();
    mountOuterLabelToggle();
    updateOverlapStatus(autoLayoutApplied ? "auto" : "manual");

    if (typeof window.resizeTree === "function") {
      window.resizeTree();
      window.setTimeout(window.resizeTree, 200);
      window.setTimeout(window.resizeTree, 600);
    }

    window.setTimeout(markTreeReady, 700);

    notifyFrame("cv:content-resize");
    window.setTimeout(function() { notifyFrame("cv:content-resize"); }, 200);
    window.setTimeout(function() { notifyFrame("cv:content-resize"); }, 600);
    notifyFrame("cv:content-ready");
  }

  function runDeferredReadyCallbacks() {
    if (readyExecuted) {
      return;
    }

    if (!baseScriptLoaded || !treeContainerVisible()) {
      window.setTimeout(runDeferredReadyCallbacks, 150);
      return;
    }

    readyExecuted = true;
    ensureOverrideStyles();
    installClassificationPatches();

    readyCallbacks.forEach(function(callback) {
      callback.call(document, window.jQuery);
    });

    installSafeWindowResize();
    window.setTimeout(finalizeRender, 0);
  }

  function loadBaseScript() {
    if (baseScriptStarted) {
      return;
    }
    baseScriptStarted = true;

    ensureOverrideStyles();
    primeOverlapStatus();

    if (!interceptJQueryReady()) {
      notifyFrame("cv:content-error");
      return;
    }

    var scriptEl = document.createElement("script");
    scriptEl.src = baseScriptSrc;
    scriptEl.async = false;
    scriptEl.onload = function() {
      baseScriptLoaded = true;
      restoreJQueryReady();
      runDeferredReadyCallbacks();
    };
    scriptEl.onerror = function() {
      restoreJQueryReady();
      notifyFrame("cv:content-error");
    };

    document.head.appendChild(scriptEl);
  }

  loadBaseScript();
})();

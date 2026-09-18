// GPCR superfamily class-cluster scatter (Plotly) + similarity side panel — used by the
// Cluster tab inlined into ClassificationSuperfamilyDetail.html. Fetches its own data via
// window.CLASSIFICATION_CLUSTER_CONFIG.embedUrl (a JSON API served by NewClassClusterTree).
//
// Only the superfamily 3-column layout is implemented — the plain/non-superfamily layout
// and the max/top3/top5/top10-mean distance-variant switcher were never reachable in
// practice (the one caller always requests the superfamily layout with the max variant),
// so they were dropped rather than carried forward as dead code.
(function () {
  const { escapeHtml } = window.ClassificationCore;
  const CONFIG = window.CLASSIFICATION_CLUSTER_CONFIG || {};
  const EMBED_URL = CONFIG.embedUrl || "";

  const clusterPlotEl = document.getElementById('ncct-cluster-plot');
  const clusterMetaEl = document.getElementById('ncct-cluster-meta');
  const statusEl = document.getElementById('ncct-status');
  const superfamilyLabelsEl = document.getElementById('ncct-superfamily-labels');
  const superfamilySimilarityTitleEl = document.getElementById('ncct-superfamily-similarity-title');
  const superfamilyHoverEl = document.getElementById('ncct-superfamily-hover');

  let currentPayload = null;
  let hoveredClassId = null;

  function titleCaseFragment(value) {
    const normalizedValue = String(value || '').trim().toLowerCase();
    if (normalizedValue === 'tetrapod-specific olfactory receptors') {
      return 'Tetrapod-specific olfactory receptors';
    }
    if (normalizedValue === 'fish-like olfactory receptors') {
      return 'Fish-like olfactory receptors';
    }
    return String(value || '').replace(/\b([A-Za-z][A-Za-z'-]*)/g, function(match) {
      return match.charAt(0).toUpperCase() + match.slice(1).toLowerCase();
    });
  }

  function normalizeClassName(name) {
    return String(name || '')
      .replace('Class O1 (fish-like)', 'Class O1 (fish-like olfactory receptors)')
      .replace('Class O2 (tetrapod specific)', 'Class O2 (tetrapod-specific olfactory receptors)');
  }

  function formattedClassNameHtml(name) {
    const raw = normalizeClassName(name).trim();
    if (!raw) return '';
    const match = raw.match(/^(.*?)(\s*\(([^)]+)\))$/);
    if (!match) return escapeHtml(raw);
    const main = escapeHtml(match[1].trim());
    const secondary = escapeHtml('(' + titleCaseFragment(match[3].trim()) + ')');
    return main + '<span class="ncct-superfamily-secondary">' + secondary + '</span>';
  }

  function formatPercent(value) {
    if (value === null || value === undefined || value === '') return 'n/a';
    const rounded = Math.round(Number(value) * 10) / 10;
    if (!isFinite(rounded)) return 'n/a';
    return (Math.round(rounded) === rounded ? String(Math.round(rounded)) : rounded.toFixed(1)) + '%';
  }

  function matrixEntryFor(sourceId, targetId) {
    const matrix = (currentPayload || {}).matrix || [];
    if (!matrix[sourceId] || !matrix[sourceId][targetId]) return null;
    return matrix[sourceId][targetId];
  }

  function scatterPoints() {
    const methods = ((currentPayload || {}).scatter || {}).methods || {};
    return (methods.tsne || {}).points || [];
  }

  function buildMarkerSizes(points, activeClassId) {
    return points.map(function(point) {
      return point.id === activeClassId ? 48 : 36;
    });
  }

  function buildMarkerOpacities(points, activeClassId) {
    return points.map(function(point) {
      if (activeClassId === null || activeClassId === undefined) return 1;
      return point.id === activeClassId ? 1 : 0.72;
    });
  }

  function buildBasePlotLabels(points) {
    return points.map(function(point) {
      return '<b>' + point.symbol + '</b>';
    });
  }

  function buildLabelSizes(points, activeClassId) {
    return points.map(function(point) {
      return point.id === activeClassId ? 28 : 14;
    });
  }

  function setActiveSuperfamilyLabel(classId) {
    if (!superfamilyLabelsEl) return;
    const buttons = superfamilyLabelsEl.querySelectorAll('.ncct-superfamily-label');
    Array.prototype.forEach.call(buttons, function(button) {
      const isActive = Number(button.getAttribute('data-class-id')) === classId;
      button.classList.toggle('is-active', isActive);
    });
  }

  function updateSuperfamilySimilarityTitle(classId) {
    if (!superfamilySimilarityTitleEl) return;
    const classes = (currentPayload || {}).classes || [];
    const classRow = classes[classId];
    if (!classRow) {
      superfamilySimilarityTitleEl.textContent = 'Similarity ranking';
      return;
    }
    const classSymbol = escapeHtml(classRow.symbol || '');
    const formattedName = formattedClassNameHtml(classRow.name || classRow.symbol || '');
    const secondaryMatch = formattedName.match(/<span class="ncct-superfamily-secondary">([\s\S]+?)<\/span>/);
    const secondaryHtml = secondaryMatch ? secondaryMatch[0] : '';
    if (classSymbol) {
      superfamilySimilarityTitleEl.innerHTML = 'Class ' + classSymbol + ' similarity ranking' + secondaryHtml;
      return;
    }
    superfamilySimilarityTitleEl.textContent = 'Class similarity ranking';
  }

  function renderSuperfamilyHoverCard(classId) {
    if (!superfamilyHoverEl) return;
    const classes = (currentPayload || {}).classes || [];
    const classRow = classes[classId];
    if (!classRow) {
      superfamilyHoverEl.innerHTML = '<p class="ncct-superfamily-hover-empty">Hover a class point to inspect its max similarity against the rest of the GPCR superfamily.</p>';
      return;
    }

    const peers = [];
    classes.forEach(function(otherClass, otherId) {
      if (otherId === classId) return;
      const cell = matrixEntryFor(classId, otherId) || {};
      peers.push({
        id: otherId,
        symbol: otherClass.symbol,
        name: otherClass.name,
        color: otherClass.color || '#808080',
        similarity: cell.similarity,
        similarityDisplay: cell.similarity_display,
        distance: cell.distance,
      });
    });
    peers.sort(function(a, b) {
      const simA = (a.similarity === null || a.similarity === undefined) ? -Infinity : Number(a.similarity);
      const simB = (b.similarity === null || b.similarity === undefined) ? -Infinity : Number(b.similarity);
      return simB - simA;
    });

    const rowsHtml = peers.map(function(peer) {
      return (
        '<div class="ncct-superfamily-hover-row">' +
          '<div class="ncct-superfamily-hover-label">' +
            '<span class="ncct-superfamily-label-swatch" style="background:' + peer.color + ';"></span>' +
            '<div class="ncct-superfamily-hover-text"><strong>' + peer.symbol + '</strong></div>' +
          '</div>' +
          '<div class="ncct-superfamily-hover-metric ncct-superfamily-hover-value-box">' +
            '<div class="ncct-superfamily-hover-metric-label">Similarity</div>' +
            '<div class="ncct-superfamily-hover-metric-value ncct-superfamily-hover-value">' + escapeHtml(peer.similarityDisplay || formatPercent(peer.similarity)) + '</div>' +
          '</div>' +
        '</div>'
      );
    }).join('');
    superfamilyHoverEl.innerHTML = '<div class="ncct-superfamily-hover-list">' + rowsHtml + '</div>';
  }

  function updateSuperfamilyPointHighlight(classId) {
    if (!clusterPlotEl || !currentPayload) return;
    const points = scatterPoints();
    if (!points.length) return;
    const markerSizes = buildMarkerSizes(points, classId);
    const markerLineWidths = points.map(function(point) { return point.id === classId ? 2.5 : 1; });
    const markerOpacities = buildMarkerOpacities(points, classId);
    Plotly.restyle(clusterPlotEl, {
      'marker.size': [markerSizes],
      'marker.line.width': [markerLineWidths],
      'marker.opacity': [markerOpacities],
    }, [0]);
    Plotly.restyle(clusterPlotEl, {
      text: [buildBasePlotLabels(points)],
      'textfont.size': [buildLabelSizes(points, classId)],
    }, [1]);
  }

  function setHoveredClass(classId) {
    hoveredClassId = classId;
    setActiveSuperfamilyLabel(classId);
    updateSuperfamilySimilarityTitle(classId);
    renderSuperfamilyHoverCard(classId);
    updateSuperfamilyPointHighlight(classId);
  }

  function renderSuperfamilyLabels(classes) {
    if (!superfamilyLabelsEl) return;
    superfamilyLabelsEl.innerHTML = '';
    (classes || []).forEach(function(cls, idx) {
      const button = document.createElement('button');
      button.type = 'button';
      button.className = 'ncct-superfamily-label';
      button.setAttribute('data-class-id', idx);
      button.innerHTML =
        '<span class="ncct-superfamily-label-marker">' +
          '<span class="ncct-superfamily-label-swatch" style="background:' + (cls.color || '#808080') + ';"></span>' +
          '<span class="ncct-superfamily-label-symbol">' + (cls.symbol || '') + '</span>' +
        '</span>' +
        '<span class="ncct-superfamily-label-text">' + formattedClassNameHtml(cls.name) + '</span>';
      button.addEventListener('mouseenter', function() {
        setHoveredClass(idx);
      });
      button.addEventListener('focus', function() {
        setHoveredClass(idx);
      });
      superfamilyLabelsEl.appendChild(button);
    });
  }

  function renderCluster(payload) {
    const points = scatterPoints();
    const activeClassId = hoveredClassId;

    const trace = {
      type: 'scatter',
      mode: 'markers',
      x: points.map(p => p.x),
      y: points.map(p => p.y),
      customdata: points.map(p => p.id),
      hoverinfo: 'none',
      marker: {
        size: buildMarkerSizes(points, activeClassId),
        color: points.map(p => p.color),
        opacity: buildMarkerOpacities(points, activeClassId),
        line: {
          width: 1,
          color: '#222',
        },
      },
    };

    const labelTrace = {
      type: 'scatter',
      mode: 'text',
      x: points.map(p => p.x),
      y: points.map(p => p.y),
      text: buildBasePlotLabels(points),
      textposition: 'middle center',
      textfont: {
        size: buildLabelSizes(points, activeClassId),
        color: '#111',
        family: 'Source Code Pro, monospace',
      },
      hoverinfo: 'skip',
    };

    const layout = {
      margin: { l: 24, r: 24, t: 12, b: 36 },
      paper_bgcolor: '#fff',
      plot_bgcolor: '#fff',
      showlegend: false,
      hovermode: 'closest',
      xaxis: { showgrid: false, zeroline: false, showticklabels: false, title: '', fixedrange: true },
      yaxis: { showgrid: false, zeroline: false, showticklabels: false, title: '', fixedrange: true },
      dragmode: false,
      font: {
        family: 'Source Code Pro, monospace',
        color: '#111',
      },
    };

    const config = {
      responsive: true,
      displayModeBar: false,
      scrollZoom: false,
      doubleClick: false,
    };

    Plotly.newPlot(clusterPlotEl, [trace, labelTrace], layout, config);
    if (typeof clusterPlotEl.removeAllListeners === 'function') {
      clusterPlotEl.removeAllListeners('plotly_hover');
      clusterPlotEl.removeAllListeners('plotly_unhover');
    }
    if (typeof clusterPlotEl.on === 'function') {
      clusterPlotEl.on('plotly_hover', function(eventData) {
        const pointData = eventData && eventData.points && eventData.points[0];
        if (!pointData) return;
        setHoveredClass(Number(pointData.customdata));
      });
      clusterPlotEl.on('plotly_unhover', function() {
        if (hoveredClassId === null || hoveredClassId === undefined) return;
        setHoveredClass(hoveredClassId);
      });
    }
    clusterMetaEl.textContent = '';
    if (hoveredClassId === null || hoveredClassId === undefined) {
      hoveredClassId = points.length ? Number(points[0].id) : null;
    }
    setHoveredClass(hoveredClassId);
  }

  function renderStatus(payload) {
    const meta = payload.meta || {};
    const notes = [];
    notes.push((meta.n_classes || 0) + ' classes');
    if (meta.missing_pairs) {
      notes.push(String(meta.missing_pairs) + ' missing pairs filled with max observed distance');
    }
    statusEl.textContent = notes.join(' | ');
  }

  async function loadPage() {
    try {
      const res = await fetch(EMBED_URL, { headers: { 'Accept': 'application/json' } });
      if (!res.ok) {
        let detail = '';
        try {
          const err = await res.json();
          if (err && err.error) detail = String(err.error);
        } catch (e) {}
        throw new Error('HTTP ' + res.status + (detail ? ': ' + detail : ''));
      }

      const payload = await res.json();
      currentPayload = payload;
      renderSuperfamilyLabels(payload.classes || []);
      renderSuperfamilyHoverCard(null);
      renderCluster(payload);
      renderStatus(payload);
    } catch (err) {
      statusEl.textContent = 'Failed to load class cluster data.';
      clusterMetaEl.textContent = String(err && err.message ? err.message : err);
    }
  }

  $('a[data-toggle="tab"]').on('shown.bs.tab', function () {
    try { Plotly.Plots.resize(clusterPlotEl); } catch (e) {}
  });

  window.addEventListener('resize', function() {
    try { Plotly.Plots.resize(clusterPlotEl); } catch (e) {}
  });

  loadPage();
})();

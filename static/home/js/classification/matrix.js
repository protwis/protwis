// Cross-class similarity matrix (heatmap) — used by the Detail page that inlines
// the matrix content. Reads its data from window.CLASSIFICATION_MATRIX_DATA and
// window.ClassificationCore for shared helpers.

(function () {
  const { escapeHtml } = window.ClassificationCore;
  const MATRIX_DATA = window.CLASSIFICATION_MATRIX_DATA || {};
  var classes = MATRIX_DATA.classes || [];
  var matrix  = MATRIX_DATA.matrix || [];
  const GPCRDB_ICON_URL = MATRIX_DATA.gpcrdbIconUrl || "";
  /* ========= DOM cache ========= */
  var axisWrap    = document.querySelector('.ccs-axis-wrap');
  var empty       = document.getElementById('ccsEmpty');
  var table       = document.getElementById('ccsTable');
  var thead       = table.querySelector('thead');
  var tbody       = table.querySelector('tbody');
  var scroller    = document.querySelector('.ccs-scroller');
  var scopeToggle = document.getElementById('ccsScopeToggle');

  var ccsFly   = document.getElementById('ccsFly');
  var ccsDT    = null;
  var overlay  = document.getElementById('ccsOverlay');

  var flyAnchor = null;
  var isLocked  = false;

  var currentRanges = null;  // { id:{min,max}, sim:{min,max} }
  var lastFlyIdFirst = null;

  var currentScope = 'human';
  var showDashInBlanks = true;

  const FLY_MAX_WIDTH_PX = 780;

  /* ========= Utilities ========= */
  function clearTable(){
    thead.innerHTML=''; tbody.innerHTML='';
    var cg=table.querySelector('colgroup'); if(cg) cg.remove();
  }
  function displayClassLabel(raw){
    return String(raw || '');
  }
  function symbolClassLabel(raw){
    var s = displayClassLabel(raw).replace(/^\s*class\s+/i, '').trim();
    var m = s.match(/^([A-Za-z0-9-]+)/);
    return m ? m[1] : s;
  }
  function twoLineLabel(raw){
    return symbolClassLabel(raw);
  }
  function isUnclassified(name){ return /\(Unclassified\)\s*$/i.test(name); }
  function flatClassLabel(raw){ return symbolClassLabel(raw); }

  var baseOrder = Array.from({length: classes.length}, (_, i) => i);

  function computeHeatRanges(indices){
    var idMin=Infinity,idMax=-Infinity, simMin=Infinity,simMax=-Infinity;
    function readVal(i,j,isUpper){
      var a=isUpper?Math.min(i,j):Math.max(i,j);
      var b=isUpper?Math.max(i,j):Math.min(i,j);
      var cell=(matrix[a]&&matrix[a][b])||{}; return cell.value;
    }
    for (var r=0;r<indices.length;r++){
      for (var c=0;c<indices.length;c++){
        if (r===c) continue;
        var i=indices[r], j=indices[c], isU=r<c, v=readVal(i,j,isU);
        if (v==null) continue;
        if (isU){ idMin=Math.min(idMin,v); idMax=Math.max(idMax,v); }
        else    { simMin=Math.min(simMin,v); simMax=Math.max(simMax,v); }
      }
    }
    if (idMin===Infinity){ idMin=null; idMax=null; }
    if (simMin===Infinity){ simMin=null; simMax=null; }
    return { id:{min:idMin,max:idMax}, sim:{min:simMin,max:simMax} };
  }
  function shade(val, range, kind){
    if (val==null || !range || range.min==null || range.max==null) return null;
    var t = (range.max===range.min)?0.6:(val-range.min)/(range.max-range.min);
    t = Math.max(0, Math.min(1,t));
    var hue = (kind==='id')?130:210, sat=60, L_hi=96, L_lo=70;
    var L = L_hi - t*(L_hi-L_lo);
    return 'hsl('+hue+','+sat+'%,'+L+'%)';
  }

  /* ========= Main toolbar dropdown (table/badges) ========= */
  (function(){
    var toggle = document.getElementById('DownloadDropdownToggle');
    var menu   = document.getElementById('DownloadDropdownMenu');

    function closeMenu(){ menu.style.display = 'none'; }
    function openMenu(){ menu.style.display = 'block'; }

    toggle.addEventListener('click', function(e){
      e.stopPropagation();
      if (menu.style.display === 'block') closeMenu(); else openMenu();
    });
    document.addEventListener('click', function(e){
      if (!document.getElementById('ccsDownloadDropdown').contains(e.target)) closeMenu();
    });
    document.addEventListener('keydown', function(e){
      if (e.key === 'Escape') closeMenu();
    });
  })();

  /* ========= Flyout helpers ========= */
  function esc(s){ return escapeHtml(s); }
  function allowSubsOnly(htmlish){ return esc(htmlish||'').replace(/&lt;(\/?)sub&gt;/gi,'<$1sub>'); }
  function linkOrDash(url,label){ return url?('<a href="'+esc(url)+'" target="_blank" rel="noopener">'+esc(label)+'</a>'):'<span class="muted">—</span>'; }
  function linkGtoPdbHtml(url,labelHtml){ return url?('<a href="'+esc(url)+'" target="_blank" rel="noopener">'+allowSubsOnly(labelHtml)+'</a>'):'<span class="muted">—</span>'; }
  function valOrDash(s){ return s?esc(s):'<span class="muted">—</span>'; }
  function iconLinkOrDash(url){ return url?('<a class="ccs-ico" href="'+esc(url)+'" target="_blank" rel="noopener" title="Open in GPCRdb"><img src="'+esc(GPCRDB_ICON_URL)+'" alt="GPCRdb"></a>'):'<span class="muted">—</span>'; }

  function initFlyDTIfPossible(){
    if (!window.jQuery || !jQuery.fn || !jQuery.fn.dataTable) return false;
    if (jQuery.fn.dataTable.isDataTable('#ccsFlyTable')) {
      ccsDT = jQuery('#ccsFlyTable').DataTable();
      return true;
    }
    ccsDT = jQuery('#ccsFlyTable').DataTable({
      dom:'t', paging:false, ordering:true, searching:false, info:false,
      autoWidth:false, deferRender:true, orderClasses:false,
      order: [[0,'desc'], [1,'desc'], [3,'asc']],
      columnDefs: [
        { targets:[0,1],
          render:function(data, type){
            if (type === 'sort' || type === 'type') {
              var n = (data === '' || data == null) ? NaN : parseFloat(data);
              return isNaN(n) ? -Infinity : n;
            }
            return data;
          }
        },
        { targets: 1, className: 'divider' },
        { targets: 5, className: 'divider' }
      ]
    });
    jQuery('#ccsFlyTable').on('draw.dt', function(){
      if (lastFlyIdFirst != null) colorFlyMetrics(lastFlyIdFirst);
    });
    return true;
  }

  function positionFlyout(cellEl, rowsCount){
    var wrap = axisWrap.getBoundingClientRect();
    var cell = cellEl.getBoundingClientRect();
    var viewportW = window.innerWidth || document.documentElement.clientWidth || 0;
    var viewportH = window.innerHeight || document.documentElement.clientHeight || 0;

    var W = Math.min(FLY_MAX_WIDTH_PX, Math.max(280, viewportW - 24));
    ccsFly.style.width = W + 'px';

    var rightSpace = viewportW - cell.right;
    var leftSpace  = cell.left;

    var place;
    if (rightSpace >= W + 12) place = 'right';
    else if (leftSpace >= W + 12) place = 'left';
    else {
      var tr  = cellEl.closest('tr');
      var tbd = tr && tr.parentNode;
      var idx = tbd ? Array.prototype.indexOf.call(tbd.children, tr) : 0;
      var mid = Math.floor(((rowsCount||1)-1)/2);
      place = (idx <= mid) ? 'bottom' : 'top';
    }

    var top, left;
    if (place==='right'){
      top  = cell.top - 6;
      left = cell.right + 8;
    } else if (place==='left'){
      top  = cell.top - 6;
      left = cell.left - 8 - W;
    } else if (place==='top'){
      top  = cell.top - 8 - ccsFly.offsetHeight;
      left = (cell.left + cell.right)/2 - W/2;
    } else { // bottom
      top  = cell.bottom + 8;
      left = (cell.left + cell.right)/2 - W/2;
    }

    var minLeft = 8, maxLeft = viewportW - 8 - W;
    left = Math.max(minLeft, Math.min(maxLeft, left));

    var h = ccsFly.offsetHeight;
    var minTop = 8, maxTop = viewportH - 8 - h;
    top = Math.max(minTop, Math.min(maxTop, top));

    ccsFly.style.left = (left - wrap.left) + 'px';
    ccsFly.style.top = (top - wrap.top) + 'px';
  }

  function onEscClose(e){ if (e.key === 'Escape') unlockFlyout(); }
  function showOverlay(){ if (overlay) overlay.classList.add('active'); }
  function hideOverlay(){ if (overlay) overlay.classList.remove('active'); }

  // ===== Flyout Download button + menu =====
  var flyDlBtn = document.getElementById('ccsFlyDownloadBtn');
  var flyMenu  = null;

  function ensureFlyMenu() {
    if (flyMenu) return flyMenu;
    flyMenu = document.createElement('div');
    flyMenu.id = 'ccsFlyDownloadMenu';
    flyMenu.className = 'dropdown-menu';
    flyMenu.innerHTML = `
      <div class="dropdown-item-text text-muted" style="font-weight:600; font-size:12px; text-align:center; pointer-events:none; padding:2px 0 6px;">
        — Flyout —
      </div>
      <button type="button" class="btn btn-primary btn-block" id="flyDLjpg">JPG</button>
      <button type="button" class="btn btn-primary btn-block" id="flyDLpng">PNG</button>
      <button type="button" class="btn btn-primary btn-block" id="flyDLsvg">SVG</button>
    `;
    axisWrap.appendChild(flyMenu);

    // Export: JPG
    flyMenu.querySelector('#flyDLjpg').addEventListener('click', () => {
      DataTablesImage.export('#ccsFly', 'ccs-fly', {
        format: 'jpg', scale: 3, quality: 0.95, background: '#ffffff',
        filename: 'GPCRdb_CrossClass_Flyout.jpg'
      });
      hideFlyMenu();
    });
    // Export: PNG
    flyMenu.querySelector('#flyDLpng').addEventListener('click', () => {
      DataTablesImage.export('#ccsFly', 'ccs-fly', {
        format: 'png', scale: 3, background: '#ffffff',
        filename: 'GPCRdb_CrossClass_Flyout.png'
      });
      hideFlyMenu();
    });
    // SVG (vector) — flyout
    flyMenu.querySelector('#flyDLsvg').addEventListener('click', async () => {
      hideFlyMenu();
      try {
        await DataTablesSVG.exportFlyout('#ccsFly', '#ccsFlyTable', 'GPCRdb_CrossClass_Flyout', {
          filename: 'GPCRdb_CrossClass_Flyout.svg'
        });
      } catch (e) {
        console.error('Flyout SVG export failed (vector). Falling back to raster-in-SVG.', e);
        // Fallback: rasterize then wrap in <svg>
        try {
          const canvas = await DataTablesImage.render('#ccsFly', { scale: 3, background: '#ffffff' });
          const w = canvas.width, h = canvas.height;
          const url = canvas.toDataURL('image/png');
          const svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}">
            <image href="${url}" x="0" y="0" width="${w}" height="${h}"/>
          </svg>`;
          const blob = new Blob([svg], { type: 'image/svg+xml' });
          const dl  = URL.createObjectURL(blob);
          const a = document.createElement('a'); a.href = dl; a.download = 'GPCRdb_CrossClass_Flyout.svg';
          document.body.appendChild(a); a.click(); a.remove(); URL.revokeObjectURL(dl);
        } catch (e2) {
          console.error('Flyout SVG raster fallback failed:', e2);
        }
      }
    });

    return flyMenu;
  }

  function positionFlyDownloadBtn(){
    if (!flyDlBtn || !ccsFly || !axisWrap) return;

    const flyVisible = ccsFly.style.display !== 'none';
    if (!flyVisible) { flyDlBtn.style.display = 'none'; return; }

    // measure button
    flyDlBtn.style.visibility = 'hidden';
    flyDlBtn.style.display = 'block';

    const wrapRect = axisWrap.getBoundingClientRect();
    const flyRect  = ccsFly.getBoundingClientRect();
    const gap = 8;

    let left = (flyRect.right - wrapRect.left) - flyDlBtn.offsetWidth;
    let top  = (flyRect.top   - wrapRect.top)  - flyDlBtn.offsetHeight - gap;

    left = Math.max(8, left);
    top  = Math.max(8, top);

    flyDlBtn.style.left = left + 'px';
    flyDlBtn.style.top  = top  + 'px';

    flyDlBtn.style.visibility = isLocked ? 'visible' : 'hidden';
    if (!isLocked) flyDlBtn.style.display = 'none';
  }

  function showFlyMenu() {
    const btn  = flyDlBtn;
    const menu = ensureFlyMenu();
    if (!btn || !menu) return;

    menu.style.display = 'block';
    menu.style.visibility = 'hidden';

    const wrapRect = axisWrap.getBoundingClientRect();
    const btnRect  = btn.getBoundingClientRect();

    const left = (btnRect.left - wrapRect.left) - menu.offsetWidth - 8; // open to the LEFT
    const top  = (btnRect.top  - wrapRect.top) + (btn.offsetHeight - menu.offsetHeight) / 2;

    menu.style.left = Math.max(8, left) + 'px';
    menu.style.top  = Math.max(8, top)  + 'px';
    menu.style.visibility = 'visible';
  }
  function hideFlyMenu(){ if (flyMenu) flyMenu.style.display = 'none'; }

  if (flyDlBtn && !flyDlBtn._wired){
    flyDlBtn.addEventListener('click', function(e){
      e.stopPropagation();
      const isOpen = flyMenu && flyMenu.style.display === 'block';
      hideFlyMenu();
      if (!isOpen) showFlyMenu();
    });
    flyDlBtn._wired = true;
  }
  // close the mini menu on outside click / Esc
  document.addEventListener('click', (e) => {
    if (!flyMenu || flyMenu.style.display !== 'block') return;
    if (flyMenu.contains(e.target) || (flyDlBtn && flyDlBtn.contains(e.target))) return;
    hideFlyMenu();
  });
  document.addEventListener('keydown', (e) => { if (e.key === 'Escape') hideFlyMenu(); });

  /* ===== Lock / unlock flyout (incl. button+menu lifecycle) ===== */
  function lockFlyout(){
    if (!ccsFly) return;
    isLocked = true;
    ccsFly.classList.add('locked');
    showOverlay();
    flyAnchor = null;
    document.addEventListener('keydown', onEscClose);
    if (flyDlBtn){
      flyDlBtn.style.display = 'block';
      positionFlyDownloadBtn();
    }
  }
  function unlockFlyout(){
    isLocked = false;
    if (ccsFly) ccsFly.classList.remove('locked');
    hideOverlay();
    document.removeEventListener('keydown', onEscClose);
    hideFlyout(true);
    if (flyDlBtn) flyDlBtn.style.display = 'none';
    hideFlyMenu();
  }
  function hideFlyout(force){
    if (!ccsFly) return;
    if (isLocked && !force) return;
    ccsFly.style.display = 'none';
    flyAnchor = null;
    if (flyDlBtn) flyDlBtn.style.display = 'none';
    hideFlyMenu();
  }

  // grey overlay click closes again
  if (overlay && !overlay._wired){
    overlay.addEventListener('click', unlockFlyout);
    overlay._wired = true;
  }

  function hideFlyoutIfAnchor(el){ if (flyAnchor === el) hideFlyout(); }

  /* ========= Populate flyout ========= */
  function rowsFromItems(items, opts){
    var idFirst = !!(opts && opts.idFirst);
    return items.map(function(it){
      var L = it.ref || {}, R = it.target || {};
      var sim = (typeof it.similarity === 'number') ? it.similarity : '';
      var idn = (typeof it.identity   === 'number') ? it.identity   : '';
      var m0 = idFirst ? idn : sim;
      var m1 = idFirst ? sim : idn;
      return [
        m0, m1,
        iconLinkOrDash(L.gpcrdb_link),
        linkGtoPdbHtml(L.gtopdb_link, L.display_name || 'Open'),
        linkOrDash(L.uniprot_link, L.uniprot || 'Open'),
        valOrDash(L.gene),
        iconLinkOrDash(R.gpcrdb_link),
        linkGtoPdbHtml(R.gtopdb_link, R.display_name || 'Open'),
        linkOrDash(R.uniprot_link, R.uniprot || 'Open'),
        valOrDash(R.gene)
      ];
    });
  }

  function setFlyHeaders(leftLabel, rightLabel, idFirst){
    var thead = document.querySelector('#ccsFlyTable thead');
    if (!thead) return;
    var L = esc(flatClassLabel(leftLabel));
    var R = esc(flatClassLabel(rightLabel));
    var m0Title = idFirst ? 'ID'  : 'SIM';
    var m1Title = idFirst ? 'SIM' : 'ID';
    var superRow =
      '<tr class="ccs-super">' +
        '<th class="metric-cap">' + m0Title + '</th>' +
        '<th class="metric-cap">' + m1Title + '</th>' +
        '<th class="group-cap" colspan="4">' + L + '</th>' +
        '<th class="group-cap" colspan="4">' + R + '</th>' +
      '</tr>';
    var baseRow =
      '<tr class="metric-row">' +
        '<th class="metric-unit">(%)</th>' +
        '<th class="metric-unit">(%)</th>' +
        '<th>GPCRdb</th><th>GtoPdb</th><th>UniProt</th><th>Gene</th>' +
        '<th>GPCRdb</th><th>GtoPdb</th><th>UniProt</th><th>Gene</th>' +
      '</tr>';
    thead.innerHTML = superRow + baseRow;
  }

  function openFlyoutForCell(cellEl, items, rowsCount, opts){
    if (!ccsFly) return;

    var idFirst    = !!(opts && opts.idFirst);
    var leftTitle  = (opts && opts.leftLabel)  || '';
    var rightTitle = (opts && opts.rightLabel) || '';

    setFlyHeaders(leftTitle, rightTitle, idFirst);

    if (!ccsDT) initFlyDTIfPossible();

    var defaultOrder = [[0,'desc'], [1,'desc'], [3,'asc']];

    if (ccsDT){
      ccsDT.clear().rows.add(rowsFromItems(items, { idFirst:idFirst }));
      ccsDT.order(defaultOrder).draw(false);
      lastFlyIdFirst = idFirst;
      colorFlyMetrics(idFirst);
    } else {
      var sorted = items.slice().sort(function(a,b){
        var a0 = idFirst ? Number(a.identity)   : Number(a.similarity);
        var b0 = idFirst ? Number(b.identity)   : Number(b.similarity);
        if ((b0||-Infinity) !== (a0||-Infinity)) return (b0||-Infinity) - (a0||-Infinity);
        var a1 = idFirst ? Number(a.similarity) : Number(a.identity);
        var b1 = idFirst ? Number(b.similarity) : Number(b.identity);
        if ((b1||-Infinity) !== (a1||-Infinity)) return (b1||-Infinity) - (a1||-Infinity);
        var la = (a && a.ref && a.ref.display_name) ? a.ref.display_name : '';
        var lb = (b && b.ref && b.ref.display_name) ? b.ref.display_name : '';
        return la.localeCompare(lb, undefined, { numeric:true, sensitivity:'base' });
      });
      var tb = document.querySelector('#ccsFlyTable tbody');
      tb.innerHTML = sorted.map(function(it){
        var row = rowsFromItems([it], { idFirst:idFirst })[0];
        return '<tr>' + row.map(function(c,i){
          var cls = (i===1 || i===5) ? ' class="divider"' : '';
          return '<td'+cls+'>' + c + '</td>';
        }).join('') + '</tr>';
      }).join('');
      lastFlyIdFirst = idFirst;
      colorFlyMetrics(idFirst);
    }

    flyAnchor = cellEl;
    ccsFly.style.display = 'block';
    positionFlyout(cellEl, rowsCount);
    positionFlyDownloadBtn();
  }

  function colorFlyMetrics(idFirst){
    if (!currentRanges) return;
    var tb = document.querySelector('#ccsFlyTable tbody');
    if (!tb) return;

    var kind0 = idFirst ? 'id'  : 'sim';
    var kind1 = idFirst ? 'sim' : 'id';

    Array.prototype.forEach.call(tb.rows, function(tr){
      var c0 = tr.cells[0], c1 = tr.cells[1];
      if (!c0 || !c1) return;

      var v0 = parseFloat((c0.textContent || '').replace(/[^\d.+-]/g,''));
      var v1 = parseFloat((c1.textContent || '').replace(/[^\d.+-]/g,''));

      var bg0 = isNaN(v0) ? null : shade(v0, currentRanges[kind0], kind0);
      var bg1 = isNaN(v1) ? null : shade(v1, currentRanges[kind1], kind1);

      c0.style.background = bg0 || '';
      c1.style.background = bg1 || '';
    });
  }

  /* ========= Scope + badges ========= */
  function applyScopeWidth(scope){
    axisWrap.classList.remove('scope-human', 'scope-human-unclassified');
    if (scope === 'human') axisWrap.classList.add('scope-human');
    else if (scope === 'human-unclassified') axisWrap.classList.add('scope-human-unclassified');
  }
  function indicesForScope(scope){
    var allowed=[]; baseOrder.forEach(function(i){
      var name=classes[i], isUnclassifiedClass=isUnclassified(name);
      if (scope==='human'){ if (!isUnclassifiedClass) allowed.push(i); }
      else if (scope==='human-unclassified'){ allowed.push(i); }
    });
    return allowed;
  }
  function positionBadges(){
    var container = document.querySelector('.ccs-axis-wrap');
    var scroller  = document.querySelector('.ccs-scroller');
    var xWrap = document.querySelector('.axis-label-x');
    var yWrap = document.querySelector('.axis-label-y');
    if (!container || !scroller || !xWrap || !yWrap) return;

    var xBadge = xWrap.querySelector('.ccs-axis-badge');
    var yBadge = yWrap.querySelector('.ccs-axis-badge');
    if (!xBadge || !yBadge) return;

    var cRect = container.getBoundingClientRect();
    var sRect = scroller.getBoundingClientRect();

    var sLeft = sRect.left - cRect.left;
    var sTop  = sRect.top  - cRect.top;

    var gapX = 8;
    xWrap.style.left = sLeft + 'px';
    xWrap.style.top  = (sTop - xBadge.offsetHeight - gapX) + 'px';
    xWrap.style.width = scroller.offsetWidth + 'px';
    xWrap.style.textAlign = 'center';

    var gapY = 10;
    var yH = yBadge.offsetHeight;
    yWrap.style.top  = (sTop + (scroller.offsetHeight - yH)/2) + 'px';
    yWrap.style.left = (sLeft - yBadge.offsetWidth - gapY) + 'px';
  }

  scopeToggle.addEventListener('click', function(){
    currentScope = (currentScope === 'human') ? 'human-unclassified' : 'human';
    // Label shows the action you can take next (not the current state)
    scopeToggle.textContent = (currentScope === 'human') ? 'Include unclassified receptors' : 'Exclude unclassified receptors';
    applyScopeWidth(currentScope);
    renderTable(indicesForScope(currentScope));
    positionBadges();
  });

  (scroller || window).addEventListener('scroll', function(){
    hideFlyout();
    positionFlyDownloadBtn();
  }, { passive:true });
  window.addEventListener('resize', function(){
    hideFlyout();
    positionBadges();
    positionFlyDownloadBtn();
  }, { passive:true });

  /* ========= Render heatmap ========= */
  if (!Array.isArray(classes) || !Array.isArray(matrix) || !classes.length) {
    if (axisWrap) axisWrap.style.display='none';
    if (empty) empty.style.display='';
    window.gpcrdbCrossClassSimilarityRefreshLayout = function(){};
    return;
  }

  function renderTable(indices){
    clearTable();
    hideFlyout();

    if (!indices.length) { axisWrap.style.display='none'; empty.style.display=''; return; }
    axisWrap.style.display=''; empty.style.display='none';

    var ranges = computeHeatRanges(indices);
    currentRanges = ranges;

    // colgroup: corner + N
    var totalCols = 1 + indices.length;
    var cg = document.createElement('colgroup');
    var col0 = document.createElement('col');
    col0.style.width = 'var(--ccs-row-header-w)'; cg.appendChild(col0);
    for (var i = 1; i < totalCols; i++){
      var col = document.createElement('col');
      col.style.width = 'var(--ccs-cell)'; cg.appendChild(col);
    }
    table.insertBefore(cg, table.firstChild);

    // THEAD
    var trHead = document.createElement('tr');
    var thCorner = document.createElement('th'); trHead.appendChild(thCorner);
    indices.forEach(function(j){
      var th = document.createElement('th');
      th.innerHTML = '<span class="ccs-col-label">' + esc(twoLineLabel(classes[j])) + '</span>';
      trHead.appendChild(th);
    });
    thead.appendChild(trHead);

    // TBODY
    indices.forEach(function(origI, rIdx){
      var tr = document.createElement('tr');

      var thRow = document.createElement('th');
      thRow.innerHTML = esc(twoLineLabel(classes[origI]));
      tr.appendChild(thRow);

      indices.forEach(function(origJ, cIdx){
        var td = document.createElement('td');

        if (rIdx === cIdx){
          td.className = 'cell-empty cell-diagonal';
          td.textContent = '';
        } else {
          var isUpper = rIdx < cIdx;
          var lookI   = isUpper ? Math.min(origI, origJ) : Math.max(origI, origJ);
          var lookJ   = isUpper ? Math.max(origI, origJ) : Math.min(origI, origJ);

          var cell = (matrix[lookI] && matrix[lookI][lookJ]) || {};
          var val  = cell.value;

          if (val == null){
            td.className = 'cell-empty';
            td.textContent = showDashInBlanks ? '—' : '';
          } else {
            td.textContent = val;
            td.className   = isUpper ? 'cell-id' : 'cell-sim';
            var color = shade(val, isUpper ? ranges.id : ranges.sim, isUpper ? 'id' : 'sim');
            if (color) td.style.setProperty('--ccs-bg', color);

            if (Array.isArray(cell.items) && cell.items.length){
              (function(tdRef, itemsRef, isUpperRef, labI, labJ){
                var leftHeader  = isUpperRef ? labI : labJ;
                var rightHeader = isUpperRef ? labJ : labI;
                var opts = {
                  idFirst:    !!isUpperRef,
                  leftLabel:  flatClassLabel(leftHeader),
                  rightLabel: flatClassLabel(rightHeader)
                };
                tdRef.addEventListener('mouseenter', function(){
                  openFlyoutForCell(tdRef, itemsRef, indices.length, opts);
                });
                tdRef.addEventListener('mouseleave', function(){ hideFlyoutIfAnchor(tdRef); });
                tdRef.addEventListener('click', function(e){
                  openFlyoutForCell(tdRef, itemsRef, indices.length, opts);
                  lockFlyout();
                  e.stopPropagation();
                });
              })(td, cell.items, isUpper, classes[origI], classes[origJ]);
            }
          }
        }

        tr.appendChild(td);
      });

      tbody.appendChild(tr);
    });

    positionBadges();
    positionFlyDownloadBtn();
  }

    /* ========= How-to modal ========= */
    var howToBtn   = document.getElementById('ccsHowToBtn');
    var howToModal = document.getElementById('ccsHowToModal');
    var howToBack  = document.getElementById('ccsHowToBackdrop');
    var howToClose = howToModal && howToModal.querySelector('.ccs-modal-close');
    var howToOk    = document.getElementById('ccsHowToOk');

    function openHowTo(){
      if (!howToModal || !howToBack) return;
      howToModal.style.display = 'block';
      howToBack.style.display  = 'block';
      howToModal.setAttribute('aria-hidden','false');
      howToBack.setAttribute('aria-hidden','false');
      // optional: prevent background scroll
      document.documentElement.style.overflow = 'hidden';
    }
    function closeHowTo(){
      if (!howToModal || !howToBack) return;
      howToModal.style.display = 'none';
      howToBack.style.display  = 'none';
      howToModal.setAttribute('aria-hidden','true');
      howToBack.setAttribute('aria-hidden','true');
      document.documentElement.style.overflow = '';
    }

    if (howToBtn)   howToBtn.addEventListener('click', openHowTo);
    if (howToBack)  howToBack.addEventListener('click', closeHowTo);
    if (howToClose) howToClose.addEventListener('click', closeHowTo);
    if (howToOk)    howToOk.addEventListener('click', closeHowTo);
    document.addEventListener('keydown', function(e){
      if (e.key === 'Escape' && howToModal && howToModal.style.display === 'block') closeHowTo();
    });


  /* ========= Initial render ========= */
  applyScopeWidth(currentScope);
  renderTable(indicesForScope(currentScope));
  positionBadges();

  // Exposed so a parent page can re-run layout once this content becomes
  // visible (e.g. after switching into its tab from display:none).
  window.gpcrdbCrossClassSimilarityRefreshLayout = function(){
    positionBadges();
    positionFlyDownloadBtn();
  };
})();

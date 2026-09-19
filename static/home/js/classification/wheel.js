// GPCRome wheel rendering — used by the Detail page that inlines the wheel content.
// Reads its data from window.CLASSIFICATION_WHEEL_DATA (set by a small inline bootstrap
// script in the Detail page template) and window.ClassificationCore for shared helpers.

// Prevent dropdown close on click inside
$(document).on('click', '.dropdown-menu', function (e) { e.stopPropagation(); });

const { normKey, fnv1a32, stableColorForKey, CLASS_COLORS: CLASS_COLORS_TREE, CHEMOTYPE_COLORS: CHEMOTYPE_COLORS_TREE, CHEMOTYPE_FALLBACK_PALETTE, getChemotypeColor } = window.ClassificationCore;

const WHEEL_DATA = window.CLASSIFICATION_WHEEL_DATA || {};
const wheelData = {
  classic: WHEEL_DATA.classicData || {},
  odorant: WHEEL_DATA.odorantData || {}
};

const ICON_URL = WHEEL_DATA.iconUrl || "";

const wheelStyling = {
  classic: { DataType:"Text", FontStyle:"Arial", FontsizeGlobal:"11px", FontsizeClass:"20px", showIcon:true, LabelType:"Protein", ShowLegend:true, LegendLayout:{mode:"row", columns:1, sorted:"Vertically"}, IconUrl: ICON_URL, ColorBy: "Class" },
  odorant: { DataType:"Text", FontStyle:"Arial", FontsizeGlobal:"11px", FontsizeClass:"20px", showIcon:true, LabelType:"Protein", ShowLegend:true, LegendLayout:{mode:"row", columns:1, sorted:"Vertically"}, IconUrl: ICON_URL, ColorBy: "Class" }
};

function buildSchemeFromCategories(categories, getColorFn) {
  const out = {};
  (categories || []).forEach(cat => {
    const k = normKey(cat);
    if (!k) return;
    out[k] = getColorFn(k);
  });
  // Ensure a stable fallback bucket exists for missing values.
  if (!out["Unknown"]) out["Unknown"] = "#cccccc";
  return out;
}

// =============================
// === Classic family renames ===
// =============================
// Renames are applied to the *wheel structure* (grouping) and also written back into
// each leaf receptor's `Receptor family` field (for coloring/legend consistency).
function preprocessClassicWheelData(data) {
  if (!data || typeof data !== "object") return data;

  const familyRename = {
    "opsins": "Vision receptors",
    "opsins receptors": "Vision receptors",
    "class a orphans": "Orphans",
    "class c orphans": "Orphans",
    "class c orphans ": "Orphans",
    "peptide p518 receptors": "QRFP receptors",
  };

  function renameFamilyKey(fam) {
    const s = normKey(fam);
    if (!s) return s;
    const mapped = familyRename[s.toLowerCase()];
    return mapped || s;
  }

  // Data shape: Circle_X -> ClassKey -> FamilyKey -> ReceptorName -> receptorObj
  Object.keys(data).forEach(circleKey => {
    const circle = data[circleKey];
    if (!circle || typeof circle !== "object") return;

    Object.keys(circle).forEach(classKey => {
      const fams = circle[classKey];
      if (!fams || typeof fams !== "object") return;

      const newFams = {};
      Object.keys(fams).forEach(oldFamKey => {
        const newFamKey = renameFamilyKey(oldFamKey);
        const receptors = fams[oldFamKey];
        if (!receptors || typeof receptors !== "object") return;

        // Update leaf metadata `Receptor family` to the renamed family.
        Object.keys(receptors).forEach(recName => {
          const r = receptors[recName];
          if (r && typeof r === "object") {
            r["Receptor family"] = newFamKey;
          }
        });

        // Merge if collision after renaming
        if (!newFams[newFamKey]) newFams[newFamKey] = {};
        Object.assign(newFams[newFamKey], receptors);
      });

      circle[classKey] = newFams;
    });
  });

  return data;
}

// Apply Classic preprocessing once, early (before extracting categories / building schemes)
wheelData.classic = preprocessClassicWheelData(wheelData.classic);

function generateCategoryColors(categories) {
  if (!categories || categories.length === 0) return {};
  const colors = {};
  const hueStep = 360 / categories.length;
  categories.forEach((cat, i) => {
    const hue = (i * hueStep) % 360;
    colors[cat] = `hsl(${hue}, 70%, 50%)`; // nice saturated rainbow
  });
  return colors;
}

function extractUniqueCategories(type, fieldName) {
  const vals = new Set();
  const data = wheelData[type] || {};
  for (const circleKey in data) {
    const circle = data[circleKey];
    for (const classKey in circle) {
      const families = circle[classKey];
      for (const familyKey in families) {
        const receptors = families[familyKey];
        for (const recKey in receptors) {
          const receptor = receptors[recKey];
          const v = receptor ? receptor[fieldName] : null;
          if (v !== undefined && v !== null && String(v).trim() !== "") {
            vals.add(String(v).trim());
          }
        }
      }
    }
  }
  return Array.from(vals).sort((a, b) => a.localeCompare(b));
}

const CLASSIC_MODALITIES = extractUniqueCategories("classic", "Modality");
const ODORANT_MODALITIES = extractUniqueCategories("odorant", "Modality");
const CLASSIC_CHEMOTYPES = extractUniqueCategories("classic", "Chemotype");
const ODORANT_CHEMOTYPES = extractUniqueCategories("odorant", "Chemotype");
const CLASSIC_SENSES = extractUniqueCategories("classic", "Sense");
const ODORANT_SENSES = extractUniqueCategories("odorant", "Sense");
const CLASSIC_RECEPTOR_FAMILIES = extractUniqueCategories("classic", "Receptor family");
const ODORANT_RECEPTOR_FAMILIES_DYNAMIC = extractUniqueCategories("odorant", "Receptor family");

// Fixed palettes (curated)
const MODALITY_COLORS_FIXED = {
  "Small molecule receptors": "#1f78b4",
  "Peptide receptors": "#2CA02C",
  "Protein receptors": "#FF7F0E",
  "Orphan receptors": "#a0b8ba",
};

const SENSE_COLORS_FIXED = {
  "Vision": "#ff7f00",
  "Taste": "#377eb8",
  "Odorant": "#4daf4a",
  "Non-sensory": "#33a02c",
  "Unknown": "#cccccc",
};

function buildSenseScheme(categories) {
  const out = {};
  (categories || []).forEach(s => {
    const k = normKey(s);
    if (!k) return;
    out[k] = SENSE_COLORS_FIXED[k] || stableColorForKey(k, CHEMOTYPE_FALLBACK_PALETTE);
  });
  // Ensure common fixed keys exist even if absent in data
  Object.keys(SENSE_COLORS_FIXED).forEach(k => {
    if (!out[k]) out[k] = SENSE_COLORS_FIXED[k];
  });
  return out;
}

const coloringSchemes = {
  "Class": {
    "A": CLASS_COLORS_TREE["A"],
    "B1": CLASS_COLORS_TREE["B1"],
    "B2": CLASS_COLORS_TREE["B2"],
    "C": CLASS_COLORS_TREE["C"],
    "F": CLASS_COLORS_TREE["F"],
    "T2": CLASS_COLORS_TREE["T2"],
    "V": CLASS_COLORS_TREE["V"],
    "U": CLASS_COLORS_TREE["U"],
  },
  "Chemotype": buildSchemeFromCategories(CLASSIC_CHEMOTYPES, getChemotypeColor),
  "Modality": buildSchemeFromCategories(CLASSIC_MODALITIES, (m) => MODALITY_COLORS_FIXED[m] || stableColorForKey(m, CHEMOTYPE_FALLBACK_PALETTE)),
  "Receptor family": generateCategoryColors(CLASSIC_RECEPTOR_FAMILIES),
  "Sense": buildSenseScheme(CLASSIC_SENSES),
};

// Add odorant-specific mappings
const ODORANT_RECEPTOR_FAMILIES = [
  "Odorant family 1", "Odorant family 2", "Odorant family 3", "Odorant family 4",
  "Odorant family 5", "Odorant family 6", "Odorant family 7", "Odorant family 8",
  "Odorant family 9", "Odorant family 10", "Odorant family 11", "Odorant family 12",
  "Odorant family 13", "Odorant family 14", "Odorant family 51", "Odorant family 52", "Odorant family 56"
];

const coloringSchemesOdorant = {
  "Class": {
    "O1": CLASS_COLORS_TREE["O1"],
    "O2": CLASS_COLORS_TREE["O2"]
  },
  "Chemotype": buildSchemeFromCategories(ODORANT_CHEMOTYPES, getChemotypeColor),
  "Modality": buildSchemeFromCategories(ODORANT_MODALITIES, (m) => MODALITY_COLORS_FIXED[m] || stableColorForKey(m, CHEMOTYPE_FALLBACK_PALETTE)),
  // Keep odorant receptor-family coloring the same behavior as before (generated).
  "Receptor family": generateCategoryColors(ODORANT_RECEPTOR_FAMILIES),
  "Sense": buildSenseScheme(ODORANT_SENSES),
};

// ----------------------------
// Legend/category ordering
// ----------------------------
function getLastLegendLabelForColorBy(colorBy) {
  const cb = normKey(colorBy);
  if (cb === "Class") return "U";
  if (cb === "Modality" || cb === "Chemotype") return "Orphan receptors";
  if (cb === "Sense") return "Unknown";
  return null;
}

function sortLegendLabels(labels, colorBy) {
  const last = getLastLegendLabelForColorBy(colorBy);
  const collator = new Intl.Collator(undefined, { numeric: true, sensitivity: "base" });
  const arr = (labels || []).slice();
  arr.sort((a, b) => {
    const A = normKey(a);
    const B = normKey(b);
    const aLast = last && A.toLowerCase() === String(last).toLowerCase();
    const bLast = last && B.toLowerCase() === String(last).toLowerCase();
    if (aLast !== bLast) return aLast ? 1 : -1; // last goes to the end
    return collator.compare(A, B);
  });
  return arr;
}

function displayLegendLabel(label, colorBy) {
  const raw = normKey(label);
  if (!raw) return raw;
  return raw;
}

function applyColoring(type) {
  const schemeSet = (type === "odorant") ? coloringSchemesOdorant : coloringSchemes;
  const colorBy = wheelStyling[type].ColorBy;
  const scheme = schemeSet[colorBy];
  if (!scheme) return;

  const enabledMap = (activeCategories[type] && activeCategories[type][colorBy]) || {};

  for (const circleKey in wheelData[type]) {
    const circle = wheelData[type][circleKey];
    for (const classKey in circle) {
      const families = circle[classKey];
      for (const familyKey in families) {
        const receptors = families[familyKey];
        for (const recKey in receptors) {
          const receptor = receptors[recKey];
          let categoryValue = normKey(receptor ? receptor[colorBy] : "");
          // Normalize missing values into a consistent bucket so they are colorable/togglable.
          if (!categoryValue && (colorBy === "Chemotype" || colorBy === "Modality" || colorBy === "Sense")) {
            categoryValue = "Unknown";
          }

          // Always keep data so the wheel can render even if everything is "off"
          receptor.Data = categoryValue || "";

          // Only change the color when disabled
          const enabled = enabledMap[categoryValue] !== false; // default true
          receptor.Color = enabled ? (scheme[categoryValue] || "#ffffff") : "#ffffff";
        }
      }
    }
  }
}




function drawWheel(type) {
  const locationId = (type === "classic") ? "GPCRome_classic_plot" : "GPCRome_odorant_plot";
  d3v4.select(`#${locationId}`).select("svg").remove();
  CustomDrawGPCRomeWheel(wheelData[type], locationId, wheelStyling[type]);
}

function shouldShowBottomLegend(wheelType) {
  const st = wheelStyling[wheelType] || {};
  return st.ShowLegend !== false && st.ColorBy !== "Class";
}

// Rebuild the text-category legend with custom ordering rules.
// (The datamapper renderer builds a legend too, but it sorts alphabetically; we replace it.)
function rebuildWheelLegend(locationId, wheelType) {
  const st = wheelStyling[wheelType] || {};
  if (!shouldShowBottomLegend(wheelType)) return;

  const colorBy = st.ColorBy;
  const schemeSet = (wheelType === "odorant") ? coloringSchemesOdorant : coloringSchemes;
  const scheme = schemeSet[colorBy];
  if (!scheme) return;

  const svg = d3v4.select("#" + locationId + "_svg");
  if (svg.empty()) return;

  // Remove any existing legend (from datamapper).
  svg.selectAll(".legend-text-categories").remove();

  const legendGroup = svg.append("g").attr("class", "legend-text-categories");

  // Match the original wheel sizing logic from datamapper.js:
  // - Base plot is 1000x1000
  // - Legend is drawn just below the wheel with a compact gap
  const baseWidth = parseFloat(svg.attr("width")) || 1000;
  const baseHeight = 1000;
  const startX = 50;
  const legendTopGap = 22;
  const legendBottomGap = 18;
  const startY = baseHeight + legendTopGap;

  // Collect the actually-used categories from the (already colored) wheelData.
  const used = new Set();
  const data = wheelData[wheelType] || {};
  for (const circleKey in data) {
    const circle = data[circleKey];
    for (const classKey in circle) {
      const families = circle[classKey];
      for (const famKey in families) {
        const receptors = families[famKey];
        for (const recKey in receptors) {
          const r = receptors[recKey];
          if (!r) continue;
          // Prefer the plotted label value (Data) since applyColoring normalizes blanks -> "Unknown"
          const v = normKey(r.Data || r[colorBy] || "");
          if (v) used.add(v);
        }
      }
    }
  }

  const labels = sortLegendLabels(Array.from(used), colorBy);
  if (!labels.length) return;

  const fontSize = st.FontsizeGlobal || "11px";
  const fontFamily = st.FontStyle || "Arial";
  const padding = 10;
  const spacingY = 25;
  const mode = st.LegendLayout?.mode || "row";
  const cols = parseInt(st.LegendLayout?.columns || "1", 10) || 1;
  const sortDirection = st.LegendLayout?.sorted || "Vertically";

  const tempText = svg.append("text")
    .attr("x", -9999).attr("y", -9999)
    .style("font-size", fontSize)
    .style("font-family", fontFamily);

  function colorForLabel(lbl) {
    const k = normKey(lbl);
    return scheme[k] || "#ffffff";
  }

  const fontSizeNumber = parseFloat(fontSize) || 11;
  const circleRadius = Math.round(fontSizeNumber * 0.4);

  if (mode === "row") {
    let x = startX, y = startY;
    const maxItemWidth = baseWidth - 100;
    labels.forEach(label => {
      const displayLabel = displayLegendLabel(label, colorBy);
      tempText.text(displayLabel);
      const labelWidth = tempText.node()?.getComputedTextLength() || 0;
      const fixedWidth = circleRadius * 2 + padding + labelWidth + 20;
      if (x + fixedWidth > baseWidth - 50) {
        x = startX;
        y += spacingY;
      }
      // Follow wheel legend convention: place circle at y, text at y+4.
      const cy = y;
      legendGroup.append("circle")
        .attr("cx", x)
        .attr("cy", cy)
        .attr("r", circleRadius)
        .style("fill", colorForLabel(label))
        .style("stroke", "black");

      legendGroup.append("text")
        .attr("x", x + circleRadius + 6)
        .attr("y", cy + 4)
        .style("font-size", fontSize)
        .style("font-family", fontFamily)
        .text(displayLabel);

      x += fixedWidth;
    });
  } else {
    // columns
    const numCols = Math.max(1, cols);
    let colData = [];
    if (sortDirection === "Horizontally") {
      colData = Array.from({ length: numCols }, () => []);
      labels.forEach((label, index) => {
        colData[index % numCols].push(label);
      });
    } else {
      const perCol = Math.floor(labels.length / numCols);
      const remainder = labels.length % numCols;
      let idx = 0;
      for (let i = 0; i < numCols; i++) {
        const count = perCol + (i < remainder ? 1 : 0);
        colData.push(labels.slice(idx, idx + count));
        idx += count;
      }
    }

    const colWidths = colData.map(col => {
      let maxW = 0;
      col.forEach(label => {
        tempText.text(displayLegendLabel(label, colorBy));
        maxW = Math.max(maxW, tempText.node()?.getComputedTextLength() || 0);
      });
      return maxW + 40;
    });

    let colStartX = startX;
    colData.forEach((col, colIndex) => {
      let x = colStartX;
      let y = startY;
      col.forEach(label => {
        const displayLabel = displayLegendLabel(label, colorBy);
        const cy = y;
        legendGroup.append("circle")
          .attr("cx", x)
          .attr("cy", cy)
          .attr("r", circleRadius)
          .style("fill", colorForLabel(label))
          .style("stroke", "black");

        legendGroup.append("text")
          .attr("x", x + circleRadius + 6)
          .attr("y", cy + 4)
          .style("font-size", fontSize)
          .style("font-family", fontFamily)
          .text(displayLabel);

        y += spacingY;
      });
      colStartX += colWidths[colIndex];
    });
  }

  tempText.remove();

  // Center the legend group horizontally (like datamapper does).
  const legendBBox = legendGroup.node()?.getBBox();
  if (legendBBox) {
    const centerOffsetX = (baseWidth - legendBBox.width) / 2 - legendBBox.x;
    legendGroup.attr("transform", `translate(${centerOffsetX}, 0)`);

    // Expand SVG width/height/viewBox so the legend sits below the plot (same idea as
    // datamapper). width/height attributes must equal the viewBox's own dimensions — a
    // viewer that scales off the attributes rather than the viewBox would otherwise crop
    // content the viewBox's negative origin shifts into the padding border.
    const addBottomHeight = (legendBBox.y + legendBBox.height) - baseHeight + legendBottomGap;
    const newHeight = baseHeight + Math.max(0, addBottomHeight);
    const vbPad = 10;
    const vbOuterW = baseWidth + 2 * vbPad;
    const vbOuterH = newHeight + 2 * vbPad;
    svg
      .attr("width", vbOuterW)
      .attr("height", vbOuterH)
      .attr("viewBox", `-${vbPad} -${vbPad} ${vbOuterW} ${vbOuterH}`);
  }
}

// =========================================================
// === Custom wheel wrapper (keeps datamapper.js untouched) ===
// =========================================================
function classDisplayLabel(code) {
  const k = normKey(code);
  if (!k) return "";
  if (k === "Unclassified") return "U";
  if (/^(A|B1|B2|C|F|O1|O2|T2|V)$/i.test(k)) return k.toUpperCase();
  return k;
}

function classStrokeColor(type, classCode) {
  const k = normKey(classCode);
  if (!k) return "#999";
  if (type === "odorant") {
    const cmap = coloringSchemesOdorant["Class"] || {};
    return cmap[k] || "#999";
  }
  return (coloringSchemes["Class"] || {})[k] || "#999";
}

function classPillFillColor(type, classCode) {
  return wheelStyling[type].ColorBy === "Class" ? classStrokeColor(type, classCode) : "#ffffff";
}

function classPillFillOpacity(type) {
  return wheelStyling[type].ColorBy === "Class" ? 0.75 : 1;
}

function addClassPills(locationId, wheelType) {
  const svg = d3v4.select("#" + locationId + "_svg");
  if (svg.empty()) return;

  // Per-class cosmetic nudges for the badge pill+text (px, SVG space: +x = right, +y = down).
  // Purely visual fine-tuning for a few badges that otherwise sit a little awkwardly on their
  // circle -- expect these to keep changing by eye as the page gets polished.
  const CLASS_BADGE_NUDGE = {
    B2: { dx: -10, dy: 5 },
    F:  { dx: -5,  dy: -10 },
  };

  // Class headers are tagged with the highlight class in datamapper's wheel renderer.
  svg.selectAll("text")
    .filter(function(d) { return d !== undefined && d !== null; })
    .filter(function(d) {
      const el = this;
      const cls = (el.getAttribute && el.getAttribute("class")) ? el.getAttribute("class") : "";
      // Only the class headers get the highlight marker
      return cls.indexOf("highlight") !== -1;
    })
    .each(function(d) {
      try {
        const txt = d3v4.select(this);
        // Replace header text (and special-case Unclassified).
        txt.text(classDisplayLabel(d));

        const fill = classPillFillColor(wheelType, d);
        const fillOpacity = classPillFillOpacity(wheelType);

        // Text styling (keep black text; no outline)
        txt.style("fill", "#000");

        // Compute bbox after text/style updates for better centering.
        const bb = this.getBBox();
        // Tight pill padding (as small as possible without clipping).
        const padX = 3;
        const padY = 1;

        // Insert behind text
        const rect = d3v4.select(this.parentNode).insert("rect", () => this)
          .attr("x", bb.x - padX)
          .attr("y", bb.y - padY)
          .attr("width", bb.width + padX * 2)
          .attr("height", bb.height + padY * 2)
          .attr("rx", 9)
          .attr("ry", 9)
          .style("fill", fill)
          .style("fill-opacity", fillOpacity)
          .style("stroke", "#000")
          .style("stroke-width", "0.75px");

        const nudge = CLASS_BADGE_NUDGE[normKey(d)];
        if (nudge) {
          const tTr = txt.attr("transform") || "";
          txt.attr("transform", (tTr ? (tTr + " ") : "") + `translate(${nudge.dx},${nudge.dy})`);
          rect.attr("transform", `translate(${nudge.dx},${nudge.dy})`);
        }
      } catch (e) {
        // ignore
      }
    });
}

function CustomDrawGPCRomeWheel(Data, locationId, GPCRome_styling) {
  const wheelType = (locationId.indexOf("odorant") !== -1) ? "odorant" : "classic";
  const effectiveStyling = Object.assign({}, GPCRome_styling, {
    ShowLegend: shouldShowBottomLegend(wheelType)
  });

  // Use the existing renderer, then extend/patch the output.
  DrawGPCRomeWheel(Data, locationId, effectiveStyling);

  addClassPills(locationId, wheelType);
  rebuildWheelLegend(locationId, wheelType);
}
// Tracks enable/disable per wheel and per scheme (Class, Chemotype, etc.)
const activeCategories = {
  classic: {},
  odorant: {}
};

// Initialize enable-state for every category in every scheme (both wheels)
function initActiveCategories() {
  ["classic", "odorant"].forEach(type => {
    const schemeSet = (type === "odorant") ? coloringSchemesOdorant : coloringSchemes;
    Object.keys(schemeSet).forEach(schemeName => {
      if (!activeCategories[type][schemeName]) activeCategories[type][schemeName] = {};
      Object.keys(schemeSet[schemeName]).forEach(category => {
        if (activeCategories[type][schemeName][category] === undefined) {
          activeCategories[type][schemeName][category] = true;
        }
      });
    });
  });
}

const customizeMenuWidths = {
  "Class": 200,
  "Chemotype": 350,
  "Modality": 260,
  "Receptor family": 350,
  "Sense": 230
};

// Build the customize panel for the current ColorBy on a given wheel & container
function initCustomizeColorPickers(type, containerId) {
  const pickerContainer = document.getElementById(containerId);
  if (!pickerContainer) return;

  const colorBy = wheelStyling[type].ColorBy;

  // Per-scheme width
  const menuEl = pickerContainer.closest(".customize-menu");
  if (menuEl) {
    const w = customizeMenuWidths[colorBy] || 350;
    menuEl.style.width = w + "px";
    menuEl.style.minWidth = w + "px";
  }

  const schemeSet = (type === "odorant") ? coloringSchemesOdorant : coloringSchemes;
  const scheme = schemeSet[colorBy];
  if (!scheme) { pickerContainer.innerHTML = "<em>No categories</em>"; return; }

  // Ensure enable-state exists
  if (!activeCategories[type][colorBy]) activeCategories[type][colorBy] = {};
  const orderedCats = sortLegendLabels(Object.keys(scheme), colorBy);
  orderedCats.forEach(cat => {
    if (activeCategories[type][colorBy][cat] === undefined) {
      activeCategories[type][colorBy][cat] = true;
    }
  });

  const safeColorBy = colorBy.replace(/[^a-zA-Z0-9_-]/g, "_");

  // Header (checkbox aligns with rows; title spans name+picker columns)
  pickerContainer.innerHTML = `
    <div class="panel-header panel-grid">
      <input type="checkbox" id="${containerId}-selectall" aria-label="Toggle all">
      <div class="header-title">${colorBy}</div>
    </div>
    <div class="color-grid" id="${containerId}-grid"></div>
  `;
  const grid = document.getElementById(`${containerId}-grid`);
  const selectAllEl = document.getElementById(`${containerId}-selectall`);

  // Master checkbox state updater
  function updateMasterCheckbox() {
    const vals = Object.values(activeCategories[type][colorBy] || {});
    const all  = vals.length > 0 && vals.every(Boolean);
    const none = vals.length > 0 && vals.every(v => !v);
    selectAllEl.checked = all;
    selectAllEl.indeterminate = !all && !none;
  }

  // Build rows
  orderedCats.forEach(category => {
    const safeCategory = category.replace(/[^a-zA-Z0-9_-]/g, "_");
    const safeId = `${type}_${safeColorBy}_picker_${safeCategory}`;
    const checkboxId = `check_${safeId}`;

    const item = document.createElement("div");
    item.className = "color-item";
    item.innerHTML = `
      <input type="checkbox" id="${checkboxId}" ${activeCategories[type][colorBy][category] ? "checked" : ""}>
      <label for="${checkboxId}" class="color-label">${category}</label>
      <input type="text" id="${safeId}">
    `;
    grid.appendChild(item);

    initializeGenericColorPicker({
      elementId: safeId,
      startColor: scheme[category],
      wheelType: type,
      schemeName: colorBy,
      categoryKey: category
    });

    document.getElementById(checkboxId).addEventListener("change", function () {
      activeCategories[type][colorBy][category] = this.checked;

      const $inp = $("#" + safeId);
      if (this.checked) {
        $inp.spectrum("enable").css("opacity", "1");
      } else {
        $inp.spectrum("disable").css("opacity", "0.5");
      }

      applyColoring(type);
      drawWheel(type);
      updateMasterCheckbox();
    });

    if (!activeCategories[type][colorBy][category]) {
      $("#" + safeId).spectrum("disable").css("opacity", "0.5");
    }
  });

  // Master toggle
  selectAllEl.addEventListener("change", function () {
    const makeChecked = this.checked; // clicking clears indeterminate
    Object.keys(activeCategories[type][colorBy]).forEach(category => {
      activeCategories[type][colorBy][category] = makeChecked;

      const safeCategory = category.replace(/[^a-zA-Z0-9_-]/g, "_");
      const safeId = `${type}_${safeColorBy}_picker_${safeCategory}`;

      const rowBox = document.getElementById(`check_${safeId}`);
      if (rowBox) rowBox.checked = makeChecked;

      const $inp = $("#" + safeId);
      if (makeChecked) { $inp.spectrum("enable").css("opacity","1"); }
      else { $inp.spectrum("disable").css("opacity","0.5"); }
    });

    applyColoring(type);
    drawWheel(type);
    updateMasterCheckbox();
  });

  updateMasterCheckbox();
}


// Spectrum init (shared)
function initializeGenericColorPicker({ elementId, startColor, wheelType, schemeName, categoryKey }) {
  const $input = $("#" + elementId);
  const appendTarget = $input.closest(".customize-menu");   // <— NEW

  $input.spectrum({
    color: startColor,
    showPalette: true,
    showInput: true,
    showButtons: false,                 // <— trims popup height
    preferredFormat: "hex",
    appendTo: 'body',
    palette: [
      ["#000", "#FF0000", "#00FF00", "#0000FF", "#FFFF00"],
      ["#FF00FF", "#00FFFF", "#FFFFFF", "#C0C0C0", "#808080"],
      ["#800000", "#808000", "#008000", "#800080", "#008080"],
      ["#000080"]
    ],
    change: function (color) {
      if (activeCategories[wheelType][schemeName][categoryKey]) {
        const schemeSet = (wheelType === "odorant") ? coloringSchemesOdorant : coloringSchemes;
        schemeSet[schemeName][categoryKey] = color.toHexString();
        if (wheelStyling[wheelType].ColorBy === schemeName) {
          applyColoring(wheelType);
          drawWheel(wheelType);
        }
      }
    },
    move: function (color) {
      if (activeCategories[wheelType][schemeName][categoryKey]) {
        const schemeSet = (wheelType === "odorant") ? coloringSchemesOdorant : coloringSchemes;
        schemeSet[schemeName][categoryKey] = color.toHexString();
        if (wheelStyling[wheelType].ColorBy === schemeName) {
          applyColoring(wheelType);
          drawWheel(wheelType);
        }
      }
    }
  });
}

function initControls(type) {
  // Label button
  const buttons = document.querySelectorAll(`.GPCRome-label-btn[data-wheel="${type}"]`);
  buttons.forEach(btn => {
    btn.addEventListener('click', function () {
      buttons.forEach(b => { b.classList.remove('btn-primary'); b.classList.add('btn-outline-primary'); });
      this.classList.remove('btn-outline-primary');
      this.classList.add('btn-primary');
      wheelStyling[type].LabelType = this.dataset.value;
      drawWheel(type);
    });
  });
  const defaultBtn = document.querySelector(`.GPCRome-label-btn[data-wheel="${type}"][data-value="${wheelStyling[type].LabelType}"]`);
  if (defaultBtn) { defaultBtn.classList.remove('btn-outline-primary'); defaultBtn.classList.add('btn-primary'); }

  // Color-by buttons
  const colorButtons = document.querySelectorAll(`.GPCRome-colorby-btn[data-wheel="${type}"]`);
  colorButtons.forEach(btn => {
    btn.addEventListener('click', function () {
      colorButtons.forEach(b => { b.classList.remove('btn-primary'); b.classList.add('btn-outline-primary'); });
      this.classList.remove('btn-outline-primary');
      this.classList.add('btn-primary');
      wheelStyling[type].ColorBy = this.dataset.value;

      applyColoring(type);
      drawWheel(type);

      // Rebuild the corresponding customize panel
      const containerId = (type === "classic") ? "classic-color-pickers" : "odorant-color-pickers";
      initCustomizeColorPickers(type, containerId);
    });
  });

  const defaultColorBtn = document.querySelector(`.GPCRome-colorby-btn[data-wheel="${type}"][data-value="${wheelStyling[type].ColorBy}"]`);
  if (defaultColorBtn) { defaultColorBtn.classList.remove('btn-outline-primary'); defaultColorBtn.classList.add('btn-primary'); }

  // Icon button
  document.getElementById(`${type}-toggleIcon`).addEventListener('click', function() {
    const st = wheelStyling[type]; st.showIcon = !st.showIcon;
    this.textContent = 'Graphical legend';
    this.classList.toggle('btn-primary', st.showIcon);
    this.classList.toggle('btn-danger', !st.showIcon);
    drawWheel(type);
  });

  $(`#${type}-legend-layout-select`).select2({ minimumResultsForSearch: Infinity, dropdownAutoWidth: true })
    .on('change', function() {
      const v = this.value;
      if (v === 'row') { wheelStyling[type].LegendLayout = {mode:'row', columns:1, sorted:'Vertically'}; }
      else if (v.startsWith('columns')) { wheelStyling[type].LegendLayout = {mode:'columns', columns:parseInt(v.split('-')[1], 10), sorted:'Vertically'}; }
      drawWheel(type);
    });

  document.getElementById(`${type}-toggleLegendSorting`).addEventListener('click', function() {
    const st = wheelStyling[type];
    st.LegendLayout.sorted = (st.LegendLayout.sorted === 'Vertically') ? 'Horizontally' : 'Vertically';
    this.textContent = st.LegendLayout.sorted;
    drawWheel(type);
  });
}

function initWheelPage() {
  initActiveCategories();

  applyColoring("classic");
  drawWheel("classic");
  initControls("classic");
  initCustomizeColorPickers("classic", "classic-color-pickers");

  applyColoring("odorant");
  initControls("odorant");
  initCustomizeColorPickers("odorant", "odorant-color-pickers");

  $('a[href="#GPCRomewheel_odorant_tab"]').on('shown.bs.tab', function () {
    applyColoring("odorant");
    drawWheel("odorant");
    initCustomizeColorPickers("odorant", "odorant-color-pickers");
  });
}

$(function () {
  initWheelPage();
});

// Shared helpers used by more than one classification visualization (wheel/tree/matrix/cluster).
// Keep this file small: only things proven identical (or safely unifiable) across pages belong
// here — anything with real per-visualization behavior stays in its own file.
(function (window) {
  "use strict";

  // Normalize a value to a trimmed string, treating missing/NaN-ish values as "".
  function normKey(v) {
    if (v === undefined || v === null) return "";
    const s = String(v).trim();
    if (!s || s.toLowerCase() === "nan") return "";
    return s;
  }

  // FNV-1a 32-bit hash, used to deterministically pick a fallback color for a key.
  function fnv1a32(str) {
    let h = 0x811c9dc5;
    for (let i = 0; i < str.length; i++) {
      h ^= str.charCodeAt(i);
      h = (h + ((h << 1) + (h << 4) + (h << 7) + (h << 8) + (h << 24))) >>> 0;
    }
    return h >>> 0;
  }

  function stableColorForKey(key, palette) {
    const k = normKey(key);
    if (!k) return "#cccccc";
    const idx = fnv1a32(k.toLowerCase()) % palette.length;
    return palette[idx];
  }

  function escapeHtml(value) {
    return String(value == null ? "" : value)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;")
      .replace(/'/g, "&#39;");
  }

  // Class-symbol -> color, same symbol used everywhere else (CLASS_VISUALIZATION_CONFIG, etc.).
  const CLASS_COLORS = {
    A: "#1f78b4",
    B1: "#33a02c",
    B2: "#6A3D9A",
    C: "#d62728",
    F: "#FF7F0E",
    O1: "#9D4EDD",
    O2: "#2A9D8F",
    T2: "#F7B6D2",
    V: "#B8860B",
    U: "#9e9e9e",
  };

  // 14-category chemotype colors (ColorBrewer Paired 12 + 2 extras). Shared verbatim
  // between the wheel and the classification tree — keep in sync if changed.
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

  const CHEMOTYPE_FALLBACK_PALETTE = [
    "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
    "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928",
    "#66c2a5", "#fc8d62",
  ];

  function getChemotypeColor(chemotype) {
    const k = normKey(chemotype);
    if (!k) return "#cccccc";
    if (CHEMOTYPE_COLORS[k]) return CHEMOTYPE_COLORS[k];
    return stableColorForKey(k, CHEMOTYPE_FALLBACK_PALETTE);
  }

  window.ClassificationCore = {
    normKey,
    fnv1a32,
    stableColorForKey,
    escapeHtml,
    CLASS_COLORS,
    CHEMOTYPE_COLORS,
    CHEMOTYPE_FALLBACK_PALETTE,
    getChemotypeColor,
  };
})(window);

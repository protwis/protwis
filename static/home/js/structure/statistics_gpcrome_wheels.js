// Drives the 3 GPCRome wheels on /structure/statistics (structure coverage, receptor-ligand
// complex counts, olfactory complex counts). All 3 share the same underlying D3 renderer
// (DrawGPCRomeWheel, datamapper.js) the classification wheel page also uses. The generic
// dropdown/legend-toggle/colors-panel/download wiring itself lives in the shared
// gpcrome_wheel_controls.js (also used by /mutations/statistics and /ligand/coverage) -- this
// file is just this page's own data/config plus the bits that really are page-specific: the
// coverage tab's categorical status swatches and its own drawn legend, and which tab lazily
// redraws when shown.
(function () {
    var DATA = window.STRUCTURE_GPCROME_WHEEL_DATA || {};

    // Starting colors for the customizable categories. "Empty" (no structure) is deliberately
    // not offered here -- it isn't a meaningful thing to recolor.
    var COVERAGE_CATEGORIES = [
        { key: "Active", label: "Active", color: "#0066ff" },
        { key: "Inactive", label: "Inactive", color: "#ff3300" },
        { key: "Both", label: "Both", color: "#ad4aad" },
    ];

    var SPECTRUM_PALETTE = [
        ["#000", "#FF0000", "#00FF00", "#0000FF", "#FFFF00"],
        ["#FF00FF", "#00FFFF", "#FFFFFF", "#C0C0C0", "#808080"],
        ["#800000", "#808000", "#008000", "#800080", "#008080"],
        ["#000080"],
    ];

    // Numeric-wheel gradient presets -- only "One"/"Two" (no "Three"/mid-color scheme; this page
    // doesn't need it). "Two" is blue (low) -> red (high), matching this page's Active/Inactive
    // categorical colors -- also the default preset shared with /mutations/statistics and
    // /ligand/coverage's wheels.
    var NUMERIC_COLOR_PRESETS = {
        One: { setup: "One", colorStart: "#ffffff", colorEnd: "#707070" },
        Two: { setup: "Two", colorStart: "#0066ff", colorEnd: "#ff3300" },
    };
    var DEFAULT_NUMERIC_PRESET_KEY = "Two";

    var complexesStats = DATA.complexesStats || { min: 0, max: 0, avg: 0 };
    var olfactoryStats = DATA.olfactoryStats || { min: 0, max: 0, avg: 0 };

    var WHEELS = {
        coverage: {
            data: DATA.coverageData,
            locationId: "GPCRome_plot",
            filenameBase: "GPCRome_plot",
            categorical: true,
            categories: COVERAGE_CATEGORIES,
            categoryColors: {},
            styling: {
                DataType: "Text", FontStyle: "Arial", FontsizeGlobal: "11px", FontsizeClass: "20px",
                showIcon: true, LabelType: "Protein", ShowLegend: false,
                LegendLayout: { mode: "row", columns: "1", sorted: "Vertically" },
            },
            drawn: false,
        },
        complexes: {
            data: DATA.complexesData,
            locationId: "GPCRome_complexes_plot",
            filenameBase: "GPCRome_complexes_plot",
            categorical: false,
            styling: {
                FontStyle: "Arial", FontsizeGlobal: "11px", FontsizeClass: "20px", DataType: "Numeric",
                showIcon: true, LabelType: "Protein",
                GPCRomeMin: complexesStats.min, GPCRomeMax: complexesStats.max, GPCRomeAvg: complexesStats.avg,
                LegendbarDigit: 0, LegendbarLength: 200, LegendbarFontsize: "11px",
                ShowLegend: true, LegendLabel: "Number of receptor-ligand complexes",
            },
            drawn: false,
        },
        olfactory: {
            data: DATA.olfactoryData,
            locationId: "GPCRome_olfactory_plot",
            filenameBase: "GPCRome_olfactory_plot",
            categorical: false,
            styling: {
                FontStyle: "Arial", FontsizeGlobal: "11px", FontsizeClass: "20px", DataType: "Numeric",
                showIcon: true, LabelType: "Protein",
                GPCRomeMin: olfactoryStats.min, GPCRomeMax: olfactoryStats.max, GPCRomeAvg: olfactoryStats.avg,
                LegendbarDigit: 0, LegendbarLength: 200, LegendbarFontsize: "11px",
                ShowLegend: true, LegendLabel: "Number of structures",
            },
            drawn: false,
        },
    };

    COVERAGE_CATEGORIES.forEach(function (cat) {
        WHEELS.coverage.categoryColors[cat.key] = cat.color;
    });

    // Seed the actual wheel data with these colors up front -- the server bakes in its own
    // status_color_map colors (structure/views.py), so without this the first render (and
    // anything before a user touches a picker) would still show the server's colors, not ours.
    function seedCoverageColors() {
        var categoryColorByKey = {};
        COVERAGE_CATEGORIES.forEach(function (cat) { categoryColorByKey[cat.key] = cat.color; });
        walkLeaves(WHEELS.coverage.data, function (leaf) {
            if (Object.prototype.hasOwnProperty.call(categoryColorByKey, leaf.Data)) {
                leaf.Color = categoryColorByKey[leaf.Data];
            }
        });
    }

    // Recursively find leaf receptor nodes (the skeleton is class -> ... -> receptor family ->
    // receptor, so a plain nested-object walk works regardless of how many levels deep that is).
    function walkLeaves(node, callback) {
        if (!node || typeof node !== "object") return;
        if (Object.prototype.hasOwnProperty.call(node, "EntryName") && Object.prototype.hasOwnProperty.call(node, "Data")) {
            callback(node);
            return;
        }
        Object.keys(node).forEach(function (key) { walkLeaves(node[key], callback); });
    }

    // ---- Coverage wheel's own legend: top-right corner (blank space above the circular plot),
    // with the "Structure coverage" title restored as part of the drawn legend instead of a
    // separate HTML block.
    function drawCoverageLegend(locationId, title, categories, styling) {
        var svg = d3v4.select("#" + locationId + "_svg");
        if (svg.empty()) return;
        svg.selectAll(".coverage-legend").remove();

        var baseWidth = parseFloat(svg.attr("width")) || 1000;
        var fontSize = (styling && styling.FontsizeGlobal) || "11px";
        var fontFamily = (styling && styling.FontStyle) || "Arial";

        // Top-right corner, in from the plot's own edge.
        var blockRight = baseWidth - 60;
        var blockY = 45; // a little clearance above the wheel's own content
        var rowHeight = 20;
        var circleRadius = 5;
        var circleLabelGap = 8;

        var legend = svg.append("g").attr("class", "coverage-legend");

        // Measure every row so the block's left edge -- where every row starts -- can be found,
        // then center the title over that same left-aligned block.
        var measure = legend.append("text").attr("x", -9999).attr("y", -9999)
            .style("font-size", fontSize).style("font-family", fontFamily);

        measure.style("font-weight", "bold").text(title);
        var titleWidth = measure.node() ? measure.node().getComputedTextLength() : 0;
        measure.style("font-weight", null);

        var rowWidths = categories.map(function (cat) {
            measure.text(cat.label);
            var labelWidth = measure.node() ? measure.node().getComputedTextLength() : 0;
            return circleRadius * 2 + circleLabelGap + labelWidth;
        });
        measure.remove();

        var maxRowWidth = rowWidths.reduce(function (m, w) { return Math.max(m, w); }, 0);
        var blockWidth = Math.max(maxRowWidth, titleWidth);
        var blockX = blockRight - blockWidth;

        legend.append("text")
            .attr("x", blockX + (blockWidth - titleWidth) / 2).attr("y", blockY)
            .style("font-size", fontSize).style("font-family", fontFamily)
            .style("font-weight", "bold")
            .text(title);

        categories.forEach(function (cat, i) {
            var y = blockY + 20 + i * rowHeight;
            legend.append("circle")
                .attr("cx", blockX + circleRadius).attr("cy", y).attr("r", circleRadius)
                .style("fill", cat.color).style("stroke", "black");
            legend.append("text")
                .attr("x", blockX + circleRadius * 2 + circleLabelGap).attr("y", y + 4)
                .style("font-size", fontSize).style("font-family", fontFamily)
                .text(cat.label);
        });
    }

    $(function () {
        seedCoverageColors();

        var controls = GPCRomeWheelControls.init({
            wheels: WHEELS,
            numericPresets: NUMERIC_COLOR_PRESETS,
            defaultPresetKey: DEFAULT_NUMERIC_PRESET_KEY,
            spectrumPalette: SPECTRUM_PALETTE,
            onAfterRedraw: function (key, wheel) {
                if (key === "coverage") {
                    var legendCategories = COVERAGE_CATEGORIES.map(function (cat) {
                        return { key: cat.key, label: cat.label, color: wheel.categoryColors[cat.key] };
                    });
                    drawCoverageLegend(wheel.locationId, "Structure coverage", legendCategories, wheel.styling);
                }
            },
        });

        controls.redraw("coverage"); // active tab on load; the other two draw lazily below

        $('a[href="#receptor_ligand_wheel_pane"]').on("shown.bs.tab", function () {
            if (!WHEELS.complexes.drawn) controls.redraw("complexes");
        });
        $('a[href="#olfactory_wheel_pane"]').on("shown.bs.tab", function () {
            if (!WHEELS.olfactory.drawn) controls.redraw("olfactory");
        });
    });
})();

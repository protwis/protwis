// Drives the 3 GPCRome wheels on /structure/statistics (structure coverage, receptor-ligand
// complex counts, olfactory complex counts). All 3 share the same underlying D3 renderer
// (DrawGPCRomeWheel, datamapper.js) the classification wheel page also uses -- this file is
// this page's own thin control layer around it, not a classification concern.
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

    var complexesStats = DATA.complexesStats || { min: 0, max: 0, avg: 0 };
    var olfactoryStats = DATA.olfactoryStats || { min: 0, max: 0, avg: 0 };

    var WHEELS = {
        coverage: {
            data: DATA.coverageData,
            locationId: "GPCRome_plot",
            filenameBase: "GPCRome_plot",
            categorical: true,
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

    // Numeric-wheel gradient presets -- only "One"/"Two" (no "Three"/mid-color scheme; this page
    // doesn't need it). "Two" is swapped from data_mapper/mapper_gpcrome_page.js's own colors to
    // blue (low) -> red (high), matching this page's Active/Inactive categorical colors.
    var NUMERIC_COLOR_PRESETS = {
        One: { setup: "One", colorStart: "#ffffff", colorEnd: "#707070" },
        Two: { setup: "Two", colorStart: "#0066ff", colorEnd: "#ff3300" },
    };
    var DEFAULT_NUMERIC_PRESET_KEY = "Two";

    // ---- Class-ring badge pills, ported from data_mapper/mapper_gpcrome_page.js's
    // mapperWheelAddClassBadgePillsAfterDraw (that page has the exact same need: neutral white
    // pills behind each class-ring label, since neither page colors "by class"). Kept inline
    // here rather than a shared file -- it's ~30 lines, not worth a whole module for.
    function badgeNormKey(code) {
        if (code === undefined || code === null) return "";
        var s = String(code).trim();
        if (!s || s.toLowerCase() === "nan") return "";
        return s;
    }

    function classDisplayShortLabel(code) {
        var k = badgeNormKey(code);
        if (!k) return "";
        if (k === "Unclassified") return "U";
        if (/^(A|B1|B2|C|F|T2|V)$/i.test(k)) return k.toUpperCase();
        return k;
    }

    // Per-class cosmetic nudges for the badge pill+text (px, SVG space: +x = right, +y = down).
    // Purely visual fine-tuning, expected to keep changing by eye.
    var CLASS_BADGE_NUDGE = {
        A: { dx: 0, dy: 1 },
        Unclassified: { dx: -15, dy: 0 },
    };

    function addClassBadgePills(locationId) {
        var svg = d3v4.select("#" + locationId + "_svg");
        if (svg.empty()) return;

        svg.selectAll("text")
            .filter(function () {
                var cls = (this.getAttribute && this.getAttribute("class")) ? this.getAttribute("class") : "";
                // datamapper: class ring labels use `GPCRome-text-{level}-highlight` only
                return /GPCRome-text-\d+-highlight/.test(cls) && cls.indexOf("GPCRome-family-label") === -1;
            })
            .each(function (d) {
                try {
                    var txt = d3v4.select(this);
                    var rawClass = badgeNormKey(d) || badgeNormKey(this.textContent);
                    if (!rawClass) return;

                    txt.text(classDisplayShortLabel(rawClass));
                    txt.style("fill", "#000");

                    var node = txt.node();
                    if (!node) return;
                    var bb = node.getBBox();
                    var padX = 3, padY = 1;
                    var g = node.parentNode;
                    if (!g || !g.insertBefore) return;

                    var rectNode = document.createElementNS("http://www.w3.org/2000/svg", "rect");
                    rectNode.setAttribute("x", String(bb.x - padX));
                    rectNode.setAttribute("y", String(bb.y - padY));
                    rectNode.setAttribute("width", String(bb.width + padX * 2));
                    rectNode.setAttribute("height", String(bb.height + padY * 2));
                    rectNode.setAttribute("rx", "9");
                    rectNode.setAttribute("ry", "9");
                    rectNode.setAttribute("fill", "#ffffff");
                    rectNode.setAttribute("fill-opacity", "1");
                    rectNode.setAttribute("stroke", "#000");
                    rectNode.setAttribute("stroke-width", "0.75px");
                    g.insertBefore(rectNode, node);

                    var nudge = CLASS_BADGE_NUDGE[rawClass];
                    if (nudge) {
                        var tPrev = txt.attr("transform") || "";
                        txt.attr("transform", (tPrev ? (tPrev + " ") : "") + "translate(" + nudge.dx + "," + nudge.dy + ")");
                        rectNode.setAttribute("transform", "translate(" + nudge.dx + "," + nudge.dy + ")");
                    }
                } catch (e) {
                    // ignore -- purely cosmetic
                }
            });
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

    function redraw(key) {
        var wheel = WHEELS[key];
        if (!wheel || !wheel.data) return;
        d3.select("#" + wheel.locationId).select("svg").remove();
        DrawGPCRomeWheel(wheel.data, wheel.locationId, wheel.styling);

        addClassBadgePills(wheel.locationId);

        if (key === "coverage") {
            var legendCategories = COVERAGE_CATEGORIES.map(function (cat) {
                return { key: cat.key, label: cat.label, color: wheel.categoryColors[cat.key] };
            });
            drawCoverageLegend(wheel.locationId, "Structure coverage", legendCategories, wheel.styling);
        }

        wheel.drawn = true;
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

    function wireLabelButtons(key) {
        var buttons = document.querySelectorAll('.GPCRomeStats-label-btn[data-wheel="' + key + '"]');
        buttons.forEach(function (btn) {
            btn.addEventListener("click", function () {
                WHEELS[key].styling.LabelType = btn.getAttribute("data-value");
                buttons.forEach(function (b) {
                    b.classList.remove("btn-primary");
                    b.classList.add("btn-outline-primary");
                });
                btn.classList.remove("btn-outline-primary");
                btn.classList.add("btn-primary");
                redraw(key);
            });
        });
        // Same "sync the matching button to the styling's current value" step
        // classification/wheel.js's initControls does -- HTML alone can't express "active"
        // for a state that only exists in JS.
        var defaultBtn = document.querySelector('.GPCRomeStats-label-btn[data-wheel="' + key + '"][data-value="' + WHEELS[key].styling.LabelType + '"]');
        if (defaultBtn) {
            defaultBtn.classList.remove("btn-outline-primary");
            defaultBtn.classList.add("btn-primary");
        }
    }

    function wireLegendToggle(key) {
        var btn = document.getElementById(key + "-toggleIcon");
        if (!btn) return;
        btn.addEventListener("click", function () {
            var wheel = WHEELS[key];
            wheel.styling.showIcon = !wheel.styling.showIcon;
            btn.classList.toggle("btn-primary", wheel.styling.showIcon);
            btn.classList.toggle("btn-danger", !wheel.styling.showIcon);
            redraw(key);
        });
    }

    function applyCategoryColor(key, categoryKey, colorValue) {
        var wheel = WHEELS[key];
        wheel.categoryColors[categoryKey] = colorValue;
        walkLeaves(wheel.data, function (leaf) {
            if (leaf.Data === categoryKey) leaf.Color = colorValue;
        });
        redraw(key);
    }

    // No enable/disable checkboxes here -- unlike classification's page (many categories, worth
    // hiding some), this is a fixed 3-status key, always all shown, so a per-row toggle wouldn't
    // mean anything.
    function initCategoricalColorPicker(key) {
        var wheel = WHEELS[key];
        var container = document.getElementById(key + "-color-pickers");
        if (!container) return;

        container.innerHTML =
            '<div class="panel-header"><div class="header-title">Colors</div></div>' +
            '<div class="color-grid" id="' + key + '-grid"></div>';

        // Half classification/wheel.css's .customize-menu default (350px) -- only a label and
        // one swatch per row here, no room needed for anything wider.
        var menuEl = container.closest(".customize-menu");
        if (menuEl) {
            menuEl.style.width = "175px";
            menuEl.style.minWidth = "175px";
        }

        var grid = document.getElementById(key + "-grid");

        COVERAGE_CATEGORIES.forEach(function (cat) {
            var safeId = key + "_picker_" + cat.key;

            var item = document.createElement("div");
            item.className = "color-item";
            item.innerHTML =
                '<span></span>' +
                '<label class="color-label">' + cat.label + '</label>' +
                '<input type="text" id="' + safeId + '">';
            grid.appendChild(item);

            $("#" + safeId).spectrum({
                color: wheel.categoryColors[cat.key],
                showPalette: true, showInput: true, showButtons: false, preferredFormat: "hex",
                appendTo: "body", palette: SPECTRUM_PALETTE,
                change: function (color) { applyCategoryColor(key, cat.key, color.toHexString()); },
                move: function (color) { applyCategoryColor(key, cat.key, color.toHexString()); },
            });
        });
    }

    // Same "Number of colors" preset dropdown + Min/Max pickers as
    // data_mapper/mapper_gpcrome_page.js's UpdateGPCRomeColorPickers/updateGPCRomeVisualization,
    // ported wheel-for-wheel (ids namespaced by `key` instead of the single "GPCRome_*" ids
    // that page hardcodes), minus the "Three"/mid-color scheme this page doesn't need.
    // select2 only treats the return value as markup when it's a jQuery/DOM object -- a plain
    // string gets escaped and inserted as text (which is exactly what was happening: the
    // selected-value box was printing the literal "<span...>" source instead of rendering it).
    // Mapper's own formatColorScheme sidesteps this by using ONE function, always $(...)-wrapped,
    // for both templateResult and templateSelection -- same fix here.
    function formatColorPresetOption(opt) {
        var preset = NUMERIC_COLOR_PRESETS[opt.id];
        if (!preset) return opt.text;
        var swatchColors = [preset.setup !== "One" ? preset.colorStart : null, preset.colorEnd];
        var swatches = swatchColors.map(function (c) {
            var bg = c || "transparent";
            var opacity = c ? "" : "opacity:0;";
            return '<span style="display:inline-block;width:14px;height:14px;margin-left:5px;border:1px solid #ccc;border-radius:2px;background:' + bg + ';' + opacity + '"></span>';
        }).join("");
        return $('<span style="display:flex; align-items:center;"><span style="min-width:50px; display:inline-block; text-align:center; padding-right:5px;">' + opt.text + '</span>' + swatches + '</span>');
    }

    function updateNumericColorPanelVisibility(key) {
        var wheel = WHEELS[key];
        var showMin = wheel.styling.ColorSetup !== "One";
        document.getElementById(key + "_Color_min_container").style.visibility = showMin ? "visible" : "hidden";
        document.getElementById(key + "_Color_min_label").style.visibility = showMin ? "visible" : "hidden";
    }

    // Only redraws if this wheel has already been drawn at least once -- called during initial
    // setup too (to seed wheel.styling before the tab is ever shown), where a redraw would just
    // be wasted work on a still-hidden plot.
    function applyNumericColorPreset(key, presetKey) {
        var wheel = WHEELS[key];
        var preset = NUMERIC_COLOR_PRESETS[presetKey];
        if (!preset) return;
        wheel.styling.ColorSetup = preset.setup;
        wheel.styling.colorStart = preset.colorStart;
        wheel.styling.colorEnd = preset.colorEnd;
        try {
            $("#" + key + "_colorPicker_min").spectrum("set", preset.colorStart);
            $("#" + key + "_colorPicker_max").spectrum("set", preset.colorEnd);
        } catch (e) { /* pickers not initialized yet */ }
        updateNumericColorPanelVisibility(key);
        if (wheel.drawn) redraw(key);
    }

    function syncNumericPickersToStyling(key) {
        var wheel = WHEELS[key];
        wheel.styling.colorStart = $("#" + key + "_colorPicker_min").spectrum("get").toHexString();
        wheel.styling.colorEnd = $("#" + key + "_colorPicker_max").spectrum("get").toHexString();
        if (wheel.drawn) redraw(key);
    }

    function initContinuousColorPicker(key, label) {
        var container = document.getElementById(key + "-color-pickers");
        if (!container) return;

        // One label row (Number of colors / Min / Max, all the same size/weight), then one
        // controls row underneath with each control vertically aligned under its own label --
        // Min and Max sit right next to each other (no reserved space for a mid/avg color,
        // since "Three" doesn't exist here).
        var LABEL_STYLE = "font-weight:bold; font-size:12px; display:block; margin-bottom:5px;";
        container.innerHTML =
            '<div class="panel-header"><div class="header-title">' + label + '</div></div>' +
            '<div style="margin:10px 0;">' +
            '<div style="display:flex; gap:15px;">' +
            '<div style="width:150px; text-align:center;"><label style="' + LABEL_STYLE + '">Number of colors</label></div>' +
            '<div style="display:flex; gap:4px;">' +
            '<div id="' + key + '_Color_min_label" style="width:80px; text-align:center;"><label style="' + LABEL_STYLE + '">Min</label></div>' +
            '<div style="width:80px; text-align:center;"><label style="' + LABEL_STYLE + '">Max</label></div>' +
            '</div>' +
            '</div>' +
            '<div style="display:flex; gap:15px;">' +
            '<div style="width:150px; display:flex; justify-content:center;">' +
            '<select id="' + key + '_color_styling" class="form-control" style="width:150px;">' +
            '<option value="One">One</option>' +
            '<option value="Two" selected>Two</option>' +
            '</select>' +
            '</div>' +
            '<div style="display:flex; gap:4px;">' +
            '<div id="' + key + '_Color_min_container" style="width:80px; display:flex; justify-content:center;"><input type="text" id="' + key + '_colorPicker_min"></div>' +
            '<div style="width:80px; display:flex; justify-content:center;"><input type="text" id="' + key + '_colorPicker_max"></div>' +
            '</div>' +
            '</div>' +
            '</div>';

        // Mapper's own numeric colors panel (#mapper-wheel-colors-dropdown-menu.mapper-core-
        // colors-show-numeric in data_mapper/wheel.css) sets min-width:420px on the dropdown-menu
        // itself for its 3-swatch (Min/Mid/Max) layout -- ours only ever needs 2 swatches now, so
        // a narrower fixed width fits without the shared classification/wheel.css .customize-menu
        // default (350px) cramping it.
        var menuEl = container.closest(".customize-menu");
        if (menuEl) {
            menuEl.style.width = "380px";
            menuEl.style.minWidth = "380px";
        }

        var defaultPreset = NUMERIC_COLOR_PRESETS[DEFAULT_NUMERIC_PRESET_KEY];
        function initPicker(elId, startColor) {
            $("#" + elId).spectrum({
                color: startColor, showPalette: true, showInput: true, showButtons: false,
                preferredFormat: "hex", appendTo: "body", palette: SPECTRUM_PALETTE,
                change: function () { syncNumericPickersToStyling(key); },
                move: function () { syncNumericPickersToStyling(key); },
            });
        }
        initPicker(key + "_colorPicker_min", defaultPreset.colorStart);
        initPicker(key + "_colorPicker_max", defaultPreset.colorEnd);

        $("#" + key + "_color_styling").select2({
            templateResult: formatColorPresetOption,
            templateSelection: formatColorPresetOption,
            // "resolve" measures the live-rendered <select> to size itself -- but the
            // complexes/olfactory tabs are still display:none at this point (lazy-drawn), so it
            // measures 0 and the widget effectively never renders. "style" reads the inline
            // width instead, which works regardless of visibility.
            width: "style",
        });
        $("#" + key + "_color_styling").on("change", function () {
            applyNumericColorPreset(key, $(this).val());
        });

        applyNumericColorPreset(key, DEFAULT_NUMERIC_PRESET_KEY);
    }

    function wireDownload(key) {
        var wheel = WHEELS[key];
        var svgBtn = document.getElementById(key + "-download-svg");
        var pngBtn = document.getElementById(key + "-download-png");
        if (svgBtn) {
            svgBtn.addEventListener("click", function () {
                GPCRomeSvgExport.downloadSvg(document.getElementById(wheel.locationId + "_svg"), wheel.filenameBase + ".svg");
            });
        }
        if (pngBtn) {
            pngBtn.addEventListener("click", function () {
                GPCRomeSvgExport.downloadSvgAsPng(document.getElementById(wheel.locationId + "_svg"), wheel.filenameBase + ".png", 2);
            });
        }
    }

    function initWheelControls(key) {
        wireLabelButtons(key);
        wireLegendToggle(key);
        wireDownload(key);
        if (WHEELS[key].categorical) {
            initCategoricalColorPicker(key);
        } else {
            initContinuousColorPicker(key, WHEELS[key].styling.LegendLabel);
        }
    }

    $(function () {
        // Same global rule Mapper_GPCRomeWheel.html uses: without it, clicking a swatch/checkbox
        // inside a dropdown-menu bubbles to Bootstrap's document-level "click outside closes the
        // dropdown" listener, closing the whole Customize colors panel mid-click. Spectrum's own
        // popup and select2's own dropdown are both appended straight to <body> (`appendTo:
        // "body"`), so they're outside .dropdown-menu's DOM subtree and need the same treatment.
        $(document).on("mousedown click", ".dropdown-menu, .sp-container, .select2-container, .select2-dropdown", function (event) {
            event.stopPropagation();
        });

        seedCoverageColors();
        initWheelControls("coverage");
        initWheelControls("complexes");
        initWheelControls("olfactory");

        redraw("coverage"); // active tab on load; the other two draw lazily below

        $('a[href="#receptor_ligand_wheel_pane"]').on("shown.bs.tab", function () {
            if (!WHEELS.complexes.drawn) redraw("complexes");
        });
        $('a[href="#olfactory_wheel_pane"]').on("shown.bs.tab", function () {
            if (!WHEELS.olfactory.drawn) redraw("olfactory");
        });
    });
})();

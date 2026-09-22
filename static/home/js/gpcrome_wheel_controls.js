// Shared control-panel wiring for a page's GPCRome wheel(s): the "Receptor names" dropdown, the
// "Graphical legend" toggle button, the "Colors" customize panel (gradient preset + min/max
// pickers, or a fixed category swatch list), and the SVG/PNG download buttons.
//
// Originally this lived duplicated per page (structure/statistics, and copy-pasted variants
// would have followed for mutations/ligand coverage) -- pulled out once here instead, following
// the same "don't re-implement a shared wheel concern per page" lesson as the class-badge pills
// in datamapper.js. Each page passes its own small `wheels` config (data/locationId/styling per
// wheel) and gets back `{ redraw, wheels }` to drive when each wheel actually draws (immediately,
// lazily on tab-show, etc -- that timing is page-specific, so it's not decided in here).
window.GPCRomeWheelControls = (function () {
    function walkLeaves(node, callback) {
        if (!node || typeof node !== "object") return;
        if (Object.prototype.hasOwnProperty.call(node, "EntryName") && Object.prototype.hasOwnProperty.call(node, "Data")) {
            callback(node);
            return;
        }
        Object.keys(node).forEach(function (key) { walkLeaves(node[key], callback); });
    }

    function init(opts) {
        var WHEELS = opts.wheels || {};
        var NUMERIC_COLOR_PRESETS = opts.numericPresets || {};
        var DEFAULT_NUMERIC_PRESET_KEY = opts.defaultPresetKey;
        var SPECTRUM_PALETTE = opts.spectrumPalette || [];
        var onAfterRedraw = opts.onAfterRedraw || function () {};

        function redraw(key) {
            var wheel = WHEELS[key];
            if (!wheel || !wheel.data) return;
            d3.select("#" + wheel.locationId).select("svg").remove();
            DrawGPCRomeWheel(wheel.data, wheel.locationId, wheel.styling);
            onAfterRedraw(key, wheel);
            wheel.drawn = true;
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
            // HTML alone can't express "active" for a state that only exists in JS -- sync the
            // matching button to the styling's current value on load.
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

        function applyCategoryColor(key, categoryKey, colorValue) {
            var wheel = WHEELS[key];
            wheel.categoryColors[categoryKey] = colorValue;
            walkLeaves(wheel.data, function (leaf) {
                if (leaf.Data === categoryKey) leaf.Color = colorValue;
            });
            redraw(key);
        }

        // No enable/disable checkboxes here -- unlike classification's page (many categories,
        // worth hiding some), a categorical wheel here has a small fixed status key, always all
        // shown, so a per-row toggle wouldn't mean anything.
        function initCategoricalColorPicker(key) {
            var wheel = WHEELS[key];
            var categories = wheel.categories || [];
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

            categories.forEach(function (cat) {
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

        // select2 only treats the return value as markup when it's a jQuery/DOM object -- a plain
        // string gets escaped and inserted as text. Always $(...)-wrap, for both templateResult
        // and templateSelection.
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

        // Only redraws if this wheel has already been drawn at least once -- called during
        // initial setup too (to seed wheel.styling before a lazily-drawn tab is ever shown),
        // where a redraw would just be wasted work on a still-hidden plot.
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
            // controls row underneath with each control vertically aligned under its own label.
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

            // A narrower fixed width fits without classification/wheel.css's .customize-menu
            // default (350px) cramping it -- only 2 swatches (Min/Max) here, no third/mid one.
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
                // "resolve" measures the live-rendered <select> to size itself -- but a lazily
                // drawn wheel's tab may still be display:none at this point, so it measures 0 and
                // the widget effectively never renders. "style" reads the inline width instead,
                // which works regardless of visibility.
                width: "style",
            });
            $("#" + key + "_color_styling").on("change", function () {
                applyNumericColorPreset(key, $(this).val());
            });

            applyNumericColorPreset(key, DEFAULT_NUMERIC_PRESET_KEY);
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

        // Runs synchronously (not deferred to a jQuery ready callback): every caller loads this
        // module via a bottom-of-body <script>, by which point the DOM is already fully parsed,
        // and callers rely on being able to call redraw() themselves right after init() returns
        // -- deferring this setup to another ready callback would race that (jQuery's ready
        // resolves via a microtask once already-ready, so it would run *after* the caller's own
        // next line, not before).

        // Without this, clicking a swatch/checkbox inside a dropdown-menu bubbles to Bootstrap's
        // document-level "click outside closes the dropdown" listener, closing the whole
        // Customize colors panel mid-click. Spectrum's own popup and select2's own dropdown are
        // both appended straight to <body> (`appendTo: "body"`), so they're outside
        // .dropdown-menu's DOM subtree and need the same treatment.
        $(document).on("mousedown click", ".dropdown-menu, .sp-container, .select2-container, .select2-dropdown", function (event) {
            event.stopPropagation();
        });

        Object.keys(WHEELS).forEach(initWheelControls);

        return { redraw: redraw, wheels: WHEELS };
    }

    return { init: init };
})();

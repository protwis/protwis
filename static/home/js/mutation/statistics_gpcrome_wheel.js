// Drives the single GPCRome wheel on /mutations/statistics. Same generic control panel as
// /structure/statistics's "Receptor-Ligand Complexes" wheel and /ligand/coverage's wheel, wired
// through the shared gpcrome_wheel_controls.js -- this file only holds what's specific to this
// page: its one wheel's data/styling.
(function () {
    var DATA = window.MUTATION_GPCROME_WHEEL_DATA || {};
    var stats = DATA.stats || { min: 0, max: 0, avg: 0 };

    var SPECTRUM_PALETTE = [
        ["#000", "#FF0000", "#00FF00", "#0000FF", "#FFFF00"],
        ["#FF00FF", "#00FFFF", "#FFFFFF", "#C0C0C0", "#808080"],
        ["#800000", "#808000", "#008000", "#800080", "#008080"],
        ["#000080"],
    ];

    // Same presets as /structure/statistics's wheels (blue/red "Two" is the default here, per
    // the "start at 1, blue-red like the structures setup" request), and shared with
    // /ligand/coverage's wheel.
    var NUMERIC_COLOR_PRESETS = {
        One: { setup: "One", colorStart: "#ffffff", colorEnd: "#707070" },
        Two: { setup: "Two", colorStart: "#0066ff", colorEnd: "#ff3300" },
    };
    var DEFAULT_NUMERIC_PRESET_KEY = "Two";

    var WHEELS = {
        mutations: {
            data: DATA.data,
            locationId: "GPCRome_plot",
            filenameBase: "GPCRome_plot",
            categorical: false,
            styling: {
                FontStyle: "Arial", FontsizeGlobal: "11px", FontsizeClass: "20px", DataType: "Numeric",
                showIcon: true, LabelType: "Protein",
                GPCRomeMin: stats.min, GPCRomeMax: stats.max, GPCRomeAvg: stats.avg,
                LegendbarDigit: 0, LegendbarLength: 200, LegendbarFontsize: "11px",
                ShowLegend: true, LegendLabel: "Mutations count",
                badgeNudge: { A: { dx: 0, dy: 1 }, Unclassified: { dx: -15, dy: 0 } },
            },
            drawn: false,
        },
    };

    var controls = GPCRomeWheelControls.init({
        wheels: WHEELS,
        numericPresets: NUMERIC_COLOR_PRESETS,
        defaultPresetKey: DEFAULT_NUMERIC_PRESET_KEY,
        spectrumPalette: SPECTRUM_PALETTE,
    });

    controls.redraw("mutations");
})();

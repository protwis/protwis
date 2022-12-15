$(document).on( "hidden.bs.modal", "#SpeciesSelector", function (e) {
    SelectionSpeciesPredefined("")
});
$(document).on("hidden.bs.modal", "#PrefGProteinsSelector", function (e) {
    SelectionGproteinPredefined("", true)
});
$(document).on("hidden.bs.modal", "#GProteinsSelector", function (e) {
    SelectionGproteinPredefined("", true)
});
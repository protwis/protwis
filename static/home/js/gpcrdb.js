/*eslint complexity: ["error", 15]*/
/*eslint wrap-iife: ["error", "outside"]*/
/*eslint quotes: ["error", "double", { "avoidEscape": true }]*/

/**
 * This is a general function to convert html tables to excel
 * Actually it converts them to xml, but excel understands it.
 * @return {tableidtag, name, filename}      <table id=tableidtag>
 */
var tableToExcel = function () {
    var uri = "data:application/vnd.ms-excel;base64,",
        template = "<html xmlns:o='urn:schemas-microsoft-com:office:office' xmlns:x='urn:schemas-microsoft-com:office:excel' xmlns='http://www.w3.org/TR/REC-html40'><head><!--[if gte mso 9]><xml><x:ExcelWorkbook><x:ExcelWorksheets><x:ExcelWorksheet><x:Name>{worksheet}</x:Name><x:WorksheetOptions><x:DisplayGridlines/></x:WorksheetOptions></x:ExcelWorksheet></x:ExcelWorksheets></x:ExcelWorkbook></xml><![endif]--></head><body><table>{table}</table></body></html>",
        base64 = function (s) {
            return window.btoa(unescape(encodeURIComponent(s)));
        }, format = function (s, c) {
            return s.replace(/{(\w+)}/g, function (m, p) {
                return c[p];
            });
        };
    return function (table, name, filename) {
            table= $("#"+table).clone();
            $("#excel_table").html(table);
            // Clean up table to remove yadcf stuff
            $("#excel_table thead tr").css("height","");
            $("#excel_table thead th").css("height","");
            $("#excel_table thead div").css("height","");
            $("#excel_table thead .yadcf-filter-wrapper").remove();
            $("#excel_table thead button").remove();
            var tr = $("#excel_table thead tr:eq(1)");
            // reattach th titles
            tr.find("th").each (function( column, th) {
                if ($(th).attr("title")) {
                    $(th).html($(th).attr("title"))
                }
            });

        var ctx = {
            worksheet: name || "Worksheet",
            table: $("#excel_table").html()
        };
        $("#excel_table").html("");
        document.getElementById("dlink").href = uri + base64(format(template, ctx));
        document.getElementById("dlink").download = filename;
        document.getElementById("dlink").click();
    }
}();


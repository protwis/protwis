/*eslint complexity: ["error", 15]*/
/*eslint wrap-iife: ["error", "outside"]*/
/*eslint quotes: ["error", "double", { "avoidEscape": true }]*/
/**
 * Gaspar's function to convert html tables to excel.
 * @return {tableidtag, name, filename}      <table id=tableidtag>
 */
let tableToExcel = function () {
    let uri = "data:application/vnd.ms-excel;base64,",
        template = "<html xmlns:o='urn:schemas-microsoft-com:office:office' xmlns:x='urn:schemas-microsoft-com:office:excel' xmlns='http://www.w3.org/TR/REC-html40'><head><!--[if gte mso 9]><xml><x:ExcelWorkbook><x:ExcelWorksheets><x:ExcelWorksheet><x:Name>{worksheet}</x:Name><x:WorksheetOptions><x:DisplayGridlines/></x:WorksheetOptions></x:ExcelWorksheet></x:ExcelWorksheets></x:ExcelWorkbook></xml><![endif]--></head><body><table>{table}</table></body></html>",
        base64 = function (s) {
            return window.btoa(unescape(encodeURIComponent(s)));
        }, format = function (s, c) {
            return s.replace(/{(\w+)}/g, function (m, p) {
                return c[p];
            });
        };
    return function (table, name, filename) {
        table = $("#" + table).clone();
        $("#excel_table").html(table);
        // Clean up table to remove yadcf stuff
        $("#excel_table thead tr").css("height", "");
        $("#excel_table thead th").css("height", "");
        $("#excel_table thead div").css("height", "");
        $("#excel_table thead .yadcf-filter-wrapper").remove();
        $("#excel_table thead button").remove();
        let tr = $("#excel_table thead tr:eq(1)");
        // reattach th titles
        tr.find("th").each(function (column, th) {
            if ($(th).attr("title")) {
                $(th).html($(th).attr("title"))
            }
        });

        let ctx = {
            worksheet: name || "Worksheet",
            table: $("#excel_table").html()
        };
        $("#excel_table").html("");
        document.getElementById("dlink").href = uri + base64(format(template, ctx));
        document.getElementById("dlink").download = filename;
        document.getElementById("dlink").click();
    }
}();

function select_all(e) {
    let checkedStatus = $(e).prop("checked");

    $(".select-all  ").each(function () {
        $(this).prop("checked", checkedStatus);
    });

    $(".alt").each(function () {
        $(this).prop("checked", checkedStatus);
    });
}

/**
 * Gaspar's function to copy to clipboard selected checkboxes as a newline separated list.
 * copied from structure_browser.js and browser_functions.js. Notice the dependency on
 * the jquery plugin PowerTip.js
 */
function copyToClipboard(array, delimiter, data_name, powertip_object = false) {
    let link = array;
//        console.log(link);
    let out = "";
    link.each(function () {
        let ele = $(this).attr("href").split("/");
        out += ele[ele.length - 1] + delimiter;
    });
    if (out.length === 0) {
        window.alert("No entries selected for copying");
        return 0;
    }
    let textArea = document.createElement("textarea");
    textArea.value = out;
    document.body.appendChild(textArea);
    textArea.focus();
    textArea.select();
    try {
        let successful = document.execCommand("copy");
        let msg = successful ? "Successful" : "Unsuccessful";
        if (powertip_object !== false) {
            $.powerTip.hide();
            powertip_object.data("powertipjq", $([
                "<p>" + array.length + " ID's copied to clipboard!</p>"
            ].join("\n")));
            powertip_object.powerTip("show");
            setTimeout(function () {
                powertip_object.data("powertipjq", $([
                    "<p>Export " + data_name + "</p>"
                ].join("\n")));
            }, 1000);
        }
    } catch (err) {
        window.alert("Oops, unable to copy");
    }
    document.body.removeChild(textArea);
}

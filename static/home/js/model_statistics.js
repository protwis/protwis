function model_statistics() {
// $(document).ready(function () {
    // "use strict";

    var oTable2 = $("#structures_scrollable").DataTable({
        "scrollY":        "65vh",
        "scrollX":        true,
        "scrollCollapse": true,
        "scroller": true,
        "paging":         false,
        "aaSorting": [],
        "autoWidth": false,
        "order": [[11,"desc"],[13,"asc"]],
        "columnDefs": [
            { "targets": "no-sort", "orderable": false }
            ],
        "columns": [
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
            null,
        ],
        "bInfo" : true,
    });

    //Uncheck every row when using back button on browser
    $(".alt_selected").prop("checked",false)
    $(".alt").prop("checked",false)
    $(".select-all").prop("checked",false)
    
    ClearSelection("targets");
    ClearSelection("reference");

    $("#loading_div").hide();

    yadcf.init(oTable2,
        [
            {
                column_number : 1,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["From", "To"],
            },
            {
                column_number : 2,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["From", "To"],
            },
            {
                column_number : 3,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["From", "To"],
            },
            {
                column_number : 4,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["From", "To"],
            },
            {
                column_number : 5,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["From", "To"],
            },
            {
                column_number : 6,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["From", "To"],
            },
            {
                column_number : 7,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["From", "To"],
            },
            {
                column_number : 8,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["From", "To"],
            },
            {
                column_number : 9,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["From", "To"],
            },
            {
                column_number : 10,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["From", "To"],
            },
            {
                column_number : 11,
                filter_type: "range_date",
                filter_reset_button_text: false,
                date_format: "yyyy-mm-dd",
                select_type: "select2",
                filter_default_label: ["From", "To"],
                // filter_reset_button_text: false,
            },
            {
                column_number : 12,
                filter_type: "text",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                // filter_match_mode : "exact",
                select_type_options: {
                    width: "70px",
                }
            },
            {
                column_number : 13,
                filter_type: "multi_select",
                select_type: "select2",
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "UniProt",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "60px",
                }
            },
            {
                column_number : 14,
                filter_type: "multi_select",
                select_type: "select2",
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "IUPHAR",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "60px",
                }
            },
            {
                column_number : 15,
                filter_type: "multi_select",
                select_type: "select2",
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "Receptor family",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "120px",
                }
            },
            {
                column_number: 16,
                filter_type: "multi_select",
                select_type: "select2",
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "Class",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "80px",
                }
            },
            {
                column_number : 17,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                column_data_type: "html",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "50px",
                }
            },
            {
                column_number : 18,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "Select",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "50px",
                }
            },
            {
                column_number : 19,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["From", "To"],
            },
            {
                column_number : 20,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["From", "To"],
            },
            {
                column_number : 21,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["From", "To"],
            },
            {
                column_number : 22,
                filter_type: "multi_select",
                select_type: "select2",
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "UniProt",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "60px",
                }
            },
            {
                column_number : 23,
                filter_type: "multi_select",
                select_type: "select2",
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "IUPHAR",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "60px",
                }
            },
            {
                column_number : 24,
                filter_type: "multi_select",
                select_type: "select2",
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "Receptor family",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "120px",
                }
            },
            {
                column_number: 25,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "Species",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "80px",
                }
            },
            {
                column_number : 26,
                filter_type: "multi_select",
                select_type: "select2",
                column_data_type: "html",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "50px",
                }
            },
            {
                column_number : 27,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "Select",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "50px",
                }
            },
            {
                column_number : 28,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["From", "To"],
            },
            {
                column_number : 29,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["From", "To"],
            },
        ],
        {
            cumulative_filtering: false
        }
    );

    yadcf.exResetAllFilters(oTable2);

    // Init gray scale
    gray_scale_table($("#structures_scrollable"));

    $("#structures_scrollable"+" > tbody > tr").click(function(event) {
        if (event.target.type !== "checkbox") {
            $(":checkbox", this).trigger("click");
            $(this).eq(0).toggleClass("alt_selected");
            $(this).find("td").toggleClass("highlight");
        }
        $(this).eq(0).toggleClass("alt_selected");
        $(this).find("td").toggleClass("highlight");
    });

    $(".select-all").click(function() {
        $(":checkbox", this).trigger("click");
        if ($(this).prop("checked")===true) {
            $(".alt").prop("checked", true);
            $(".alt").parent().parent().addClass("alt_selected");
            $(".alt").parent().parent().find("td").addClass("highlight");
        }
        if ($(this).prop("checked")===false) {
            $(".alt").prop("checked", false);
            $(".alt").parent().parent().removeClass("alt_selected");
            $(".alt").parent().parent().find("td").removeClass("highlight");
        }
    });

    // $(".wrapper").find("div").width($(".yadcf-datatables-table--structures_scrollable").width());

    $(".hide_columns").click(function(evt) {
        var columns = $(this).attr("columns").split(",");
        columns.forEach(function(column) {
            var column = oTable2.column( column );
            try {
                column.visible( false, false );
            }
            catch(err) {
                column.visible( false, false );
            }
            });
            oTable2.draw();
    });

    // Remove odd and even classes from rows to avoid alternating background. Has to be called on every table draw
    $(".odd").removeClass("odd");
    $(".even").removeClass("even");
    $("#structures_scrollable").on("draw.dt", function(e, oSettings) {
        $(".odd").removeClass("odd");
        $(".even").removeClass("even");
    });
    
    $(".uniprot-export").data("powertipjq", $([
        "<p>Export UniProt IDs</p>"
        ].join("\n")));
    $(".pdb-export").data("powertipjq", $([
        "<p>Export PDB IDs</p>"
        ].join("\n")));
    $(".glyphicon-export").powerTip({
        placement: "n",
        smartPlacement: true
    });
    $("#uniprot_copy").click(function () {
        copyToClipboard($(".alt_selected > .uniprot > a"), "\n", "UniProt IDs", $(".uniprot-export"));
    });
    $("#pdb_copy").click(function () {
        copyToClipboard($(".alt_selected > .pdb > a"), "\n", "PDB IDs", $(".pdb-export"));
    });
};

var tableToExcel = (function () {
    var uri = "data:application/vnd.ms-excel;base64,",
        template = "<html xmlns:o='urn:schemas-microsoft-com:office:office' xmlns:x='urn:schemas-microsoft-com:office:excel' xmlns='http://www.w3.org/TR/REC-html40'><head><!--[if gte mso 9]><xml><x:ExcelWorkbook><x:ExcelWorksheets><x:ExcelWorksheet><x:Name>{worksheet}</x:Name><x:WorksheetOptions><x:DisplayGridlines/></x:WorksheetOptions></x:ExcelWorksheet></x:ExcelWorksheets></x:ExcelWorkbook></xml><![endif]--></head><body><table>{table}</table></body></html>",
        base64 = function (s) {
            return window.btoa(unescape(encodeURIComponent(s)))
        }, format = function (s, c) {
            return s.replace(/{(\w+)}/g, function (m, p) {
                return c[p];
            })
        }
    return function (table, name, filename) {
            var table= $("#"+table).clone();
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
              if ($(th).attr("title")) $(th).html($(th).attr("title"));
            });

        var ctx = {
            worksheet: name || "Worksheet",
            table: $("#excel_table").html()
        }
        $("#excel_table").html("");
        document.getElementById("dlink").href = uri + base64(format(template, ctx));
        document.getElementById("dlink").download = filename;
        document.getElementById("dlink").click();
    }
})()

function copyToClipboard(array, delimiter, data_name, powertip_object=false) {
    var link = array;
    var out = "";
    link.each(function() {
        var ele = $(this).attr("href").split("/");
        out+=ele[ele.length-1]+delimiter;
    });
    if (out.length===0) {
        window.alert("No entries selected for copying")
        return 0
    }
    var textArea = document.createElement("textarea");
    textArea.value = out;
    document.body.appendChild(textArea);
    textArea.focus();
    textArea.select();
    try {
        var successful = document.execCommand("copy");
        var msg = successful ? "Successful" : "Unsuccessful";
        if (powertip_object!==false) {
            $.powerTip.hide();
            powertip_object.data("powertipjq", $([
                "<p>Copied to clipboard!</p>"
                ].join("\n")))
            powertip_object.powerTip("show");
            setTimeout(function() {
            powertip_object.data("powertipjq", $([
                "<p>Export "+data_name+"</p>"
                ].join("\n")))
            },1000);
        }
    } catch (err) {
        window.alert("Oops, unable to copy");
    }
    document.body.removeChild(textArea);
}


/*global yadcf, ClearSelection, gray_scale_table, copyToClipboard*/
/*eslint no-undef: "error"*/

function model_statistics() {
    var oTable2 = $("#structures_scrollable").DataTable({
        "scrollY":        "65vh",
        "scrollX":        true,
        "scrollCollapse": true,
        "scroller": true,
        "paging":         false,
        "aaSorting": [],
        "autoWidth": false,
        "order": [[18,"desc"],[12,"asc"]],
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
            // null,
            // null,
            // null,
            // null,
            // null,
            // null,
            // null,
            // null,
            // null,
            // null,
        ],
        "bInfo" : true,
    });

    //Uncheck every row when using back button on browser
    $(".alt_selected").prop("checked",false);
    $(".alt").prop("checked",false);
    $(".select-all").prop("checked",false);
    
    ClearSelection("targets");
    ClearSelection("reference");

    $("#loading_div").hide();

    yadcf.init(oTable2,
        [
            {
                column_number : 1,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 2,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 3,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 4,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 5,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 6,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 7,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 8,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 9,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 10,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["Min", "Max"],
            },
            
            {
                column_number : 11,
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
                column_number : 12,
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
                column_number : 13,
                filter_type: "multi_select",
                select_type: "select2",
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "Select",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "80px",
                }
            },
            {
                column_number: 14,
                filter_type: "multi_select",
                select_type: "select2",
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "35px",
                }
            },
            {
                column_number : 15,
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
                column_number : 16,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "Select",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "50px",
                }
            },
            {
                column_number : 17,
                filter_type: "range_number",
                filter_reset_button_text: false,
                filter_default_label: ["Min", "Max"],
            },
            {
                column_number : 18,
                filter_type: "range_date",
                filter_reset_button_text: false,
                date_format: "yyyy-mm-dd",
                select_type: "select2",
                filter_default_label: ["Min", "Max"],
                // filter_reset_button_text: false,
            },
            {
                column_number : 19,
                filter_type: "text",
                select_type: "select2",
                filter_default_label: "Note",
                filter_reset_button_text: false,
                // filter_match_mode : "exact",
                select_type_options: {
                    width: "70px",
                }
            },
            // {
            //     column_number : 20,
            //     filter_type: "range_number",
            //     filter_reset_button_text: false,
            //     filter_default_label: ["Min", "Max"],
            // },
            // {
            //     column_number : 21,
            //     filter_type: "range_number",
            //     filter_reset_button_text: false,
            //     filter_default_label: ["Min", "Max"],
            // },
            
            // {
            //     column_number : 22,
            //     filter_type: "multi_select",
            //     select_type: "select2",
            //     column_data_type: "html",
            //     html_data_type: "text",
            //     filter_default_label: "UniProt",
            //     filter_match_mode : "exact",
            //     filter_reset_button_text: false,
            //     select_type_options: {
            //         width: "60px",
            //     }
            // },
            // {
            //     column_number : 23,
            //     filter_type: "multi_select",
            //     select_type: "select2",
            //     column_data_type: "html",
            //     html_data_type: "text",
            //     filter_default_label: "IUPHAR",
            //     filter_match_mode : "exact",
            //     filter_reset_button_text: false,
            //     select_type_options: {
            //         width: "60px",
            //     }
            // },
            // {
            //     column_number : 24,
            //     filter_type: "multi_select",
            //     select_type: "select2",
            //     column_data_type: "html",
            //     html_data_type: "text",
            //     filter_default_label: "Select",
            //     filter_match_mode : "exact",
            //     filter_reset_button_text: false,
            //     select_type_options: {
            //         width: "80px",
            //     }
            // },
            // {
            //     column_number: 25,
            //     filter_type: "multi_select",
            //     select_type: "select2",
            //     filter_default_label: "Species",
            //     filter_reset_button_text: false,
            //     select_type_options: {
            //         width: "60px",
            //     }
            // },
            // {
            //     column_number : 26,
            //     filter_type: "multi_select",
            //     select_type: "select2",
            //     column_data_type: "html",
            //     filter_default_label: "",
            //     filter_reset_button_text: false,
            //     select_type_options: {
            //         width: "50px",
            //     }
            // },
            // {
            //     column_number : 27,
            //     filter_type: "multi_select",
            //     select_type: "select2",
            //     filter_default_label: "Select",
            //     filter_reset_button_text: false,
            //     select_type_options: {
            //         width: "50px",
            //     }
            // },
            // {
            //     column_number : 28,
            //     filter_type: "range_number",
            //     filter_reset_button_text: false,
            //     filter_default_label: ["Min", "Max"],
            // },
            // {
            //     column_number : 29,
            //     filter_type: "range_number",
            //     filter_reset_button_text: false,
            //     filter_default_label: ["Min", "Max"],
            // },
        ],
        {
            cumulative_filtering: false
        }
    );

    yadcf.exResetAllFilters(oTable2);

    // Init gray scale
    gray_scale_table($("#structures_scrollable"), ["color-set1"]);

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
        var column;
        columns.forEach(function(column) {
            column = oTable2.column( column );
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
    
    // Tooltip popups
    $(".uniprot-export1").data("powertipjq", $([
        "<p>Export UniProt IDs</p>"
        ].join("\n")));
    $(".pdb-export1").data("powertipjq", $([
        "<p>Export PDB IDs</p>"
        ].join("\n")));
    $(".uniprot-export2").data("powertipjq", $([
        "<p>Export UniProt IDs</p>"
        ].join("\n")));
    $(".pdb-export2").data("powertipjq", $([
        "<p>Export PDB IDs</p>"
        ].join("\n")));
    $(".gn-seqid").data("powertipjq", $([
        "<p>Sequence identity in % between structure and main template for GPCRdb aligned residues that have a generic number</p>"
        ].join("\n")));
    $(".gn-seqsim").data("powertipjq", $([
        "<p>Sequence similarity in % between structure and main template for GPCRdb aligned residues that have a generic number</p>"
        ].join("\n")));
    $(".glyphicon").powerTip({
        placement: "n",
        smartPlacement: true
    });
    $("#uniprot_copy1").click(function () {
        copyToClipboard($(".alt_selected > .uniprot > a:even"), "\n", "UniProt IDs", $(".uniprot-export1"));
    });
    $("#uniprot_copy2").click(function () {
        copyToClipboard($(".alt_selected > .uniprot > a:odd"), "\n", "UniProt IDs", $(".uniprot-export2"));
    });
    $("#pdb_copy1").click(function () {
        copyToClipboard($(".alt_selected > .pdb > a:even"), "\n", "PDB IDs", $(".pdb-export1"));
    });
    $("#pdb_copy2").click(function () {
        copyToClipboard($(".alt_selected > .pdb > a:odd"), "\n", "PDB IDs", $(".pdb-export2"));
    });
}

/*global yadcf*/
/*eslint complexity: ["error", 8]*/
/*eslint wrap-iife: ["error", "outside"]*/
/*eslint quotes: ["error", "double", { "avoidEscape": true }]*/
let oTable1 = [];

function select_all(e) {
    var checkedStatus = $(e).prop("checked");

    $(".select-all  ").each(function () {
        $(this).prop("checked", checkedStatus);
    });

    $(".alt").each(function () {
        $(this).prop("checked", checkedStatus);
    });
}

function resetHidden1() {
    var columns = Array.from(new Array(17), (x,i) => i + 3);
    columns.forEach(function(column) {
        column = oTable1.column( column );
        try {
            column.visible( true, false );
        }
        catch(err) {
            column.visible( true, false );
        }
    });
    oTable1.draw();
}

function reset_tab() {
// Just a button to go back to the main page.
    window.location.href = "/signprot/arrestincouplings";
}

$(document).ready(function() {
// Activate tooltips and popovers from Bootstrap   * Bootstrap v3.3.7 (http://getbootstrap.com)
    $("[data-toggle='tooltip']").tooltip();
    $("[data-toggle='popover']").popover();

//     $('a[data-toggle="tab"]').on("shown.bs.tab", function (e) {
// //      console.log( 'show tab' );
//         $($.fn.dataTable.tables(true)).DataTable()
//             .columns.adjust().responsive;
//     });

    console.time("table1load");
    oTable1 = $("#arrestintable").DataTable({
        deferRender: true,
        scrollY: "50vh",
        scrollX: true,
        scrollCollapse: true,
        scroller: true,
        paging: false,
        bSortCellsTop: false, //prevent sort arrows going on bottom row
        aaSorting: [],
        order: [[4,"asc"]],
        autoWidth: false,
        bInfo: true,
        columnDefs: [
            {
                targets: [0],
                visible: true,
                orderable: false
            }
        ],
    });

    yadcf.init(oTable1,
        [
            {
                column_number: 0,
                filter_type: "none",
                filter_default_label: "",
                filter_reset_button_text: false,
            },
            {
                column_number: 1,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "80px",
                }
            },

            {
                column_number: 2,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "40px",
                }
            },
            {
                column_number: 3,
                filter_type: "multi_select",
                select_type: "select2",
                filter_default_label: "",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "200px",
                }
            },
            {
                column_number: 4,
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
                column_number: 5,
                filter_type: "multi_select",
                select_type: "select2",
                column_data_type: "html",
                html_data_type: "text",
                filter_default_label: "",
                filter_match_mode : "exact",
                filter_reset_button_text: false,
                select_type_options: {
                    width: "80px",
                }
            },

// log(Emax/EC50)
            {
                column_number : 6,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 7,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },

// pEC50
            {
                column_number : 8,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 9,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },

// Emax
            {
                column_number : 10,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },
            {
                column_number : 11,
                filter_type: "range_number",
                filter_default_label: ["Min", "Max"],
                filter_reset_button_text: false,
            },

        ],

        {filters_tr_index: 2},

        {
            cumulative_filtering: true
        }
    );


    $("#arrestintable"+" > tbody > tr").click(function(event) {
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

    console.timeEnd("table1load");


// Hide column button for table1
    $(".hide_columns1").click(function(evt) {
        var columns = $(this).attr("columns").split(",");
        columns.forEach(function(column) {
            column = oTable1.column( column );
            try {
                column.visible( false, false );
            }
            catch(err) {
                column.visible( false, false );
            }
        });
        oTable1.draw();
    } );


// Put top scroller
// https://stackoverflow.com/questions/35147038/how-to-place-the-datatables-horizontal-scrollbar-on-top-of-the-table
//    console.time("scroll to top");
    $(".dataTables_scrollHead").css({
        "overflow-x":"scroll"
    }).on("scroll", function(e){
        var scrollBody = $(this).parent().find(".dataTables_scrollBody").get(0);
        scrollBody.scrollLeft = this.scrollLeft;
        $(scrollBody).trigger("scroll");
    });
//    console.timeEnd("scroll to top");



// =============================================================================
// START OVERLAY COLUMNS CODE HERE
// =============================================================================
    let toggle_enabled = true;
    $("#toggle_fixed_btn1").click(function() {
        if (toggle_enabled) {
            toggle_enabled = false;
            $("#overlay1").hide();
            $("#toggle_fixed_btn1").attr("value","Enable fixed columns");
            $("#toggle_fixed_btn1").addClass("clicked_button");
        } else {
            toggle_enabled = true;
            $("#arrestintable").closest(".dataTables_scrollBody").scroll();
            $("#toggle_fixed_btn1").attr("value","Disable fixed columns");
            $("#toggle_fixed_btn1").removeClass("clicked_button");
        }
    });


// --------- overlay for table 1 ---------
    var left1 = 0;
    var old_left1 = 0;
    $("#arrestintable").closest(".dataTables_scrollBody").scroll(function(){
        // If user scrolls and it's > 100px from left, then attach fixed columns overlay
        left1 = $("#arrestintable").closest(".dataTables_scrollBody").scrollLeft();
        if (left1!==old_left1) {
            $("#overlay1").hide();
        }
        old_left1 = left1;

        if (left1 > 50 && toggle_enabled) {
            $("#overlay1").css({ left: left1 + "px" });
            if ($("#overlay1").is(":hidden")) {
                $("#overlay1").show();
            }
        }
    });

    $("#arrestintable").closest(".dataTables_scrollBody").append('<div id="overlay1"><table id="overlay_table1" class="row-border text-center compact dataTable no-footer text-nowrap"><tbody></tbody></table></div>');

    function create_overlay1() {
        // This function fires upon filtering, to update what rows to show as an overlay
        $("#overlay_table1 tbody tr").remove();
        var $target = $("#overlay_table1 tbody");

        $("#arrestintable tbody tr").each(function() {
            var $tds = $(this).children(),
                $row = $("<tr></tr>");
            $row.append($tds.eq(0).clone()).append($tds.eq(1).clone()).append($tds.eq(2).clone()).append($tds.eq(3).clone()).append($tds.eq(4).clone()).appendTo($target);
            $row.height($(this).height());
            //$row.font_size("10");
            //$row.height("31px");
            //$row.append($tds.height("2.5")).appendTo($target);
        });
        $("#overlay_table1 .border-right").removeClass("border-right");
    }

    // Function that detects filtering events
    $("#arrestintable").on( "draw.dt", function (e, oSettings) {
        create_overlay1();
    });

    create_overlay1();
    $("#overlay1").hide();
// --------- End of overlay for table1 ---------


// =============================================================================
// END OVERLAY COLUMNS CODE HERE
// =============================================================================

// Gaspar's functions to copy to clipboard selected checkboxes as a newline separated list.
// copied from structure_browser.js and browser_functions.js. Notice that they depend on
// the jquery plugin PowerTip.js

    function copyToClipboard(array, delimiter, data_name, powertip_object=false) {
        var link = array;
//        console.log(link);
        var out = "";
        link.each(function() {
            var ele = $(this).attr("href").split("/");
            out+=ele[ele.length-1]+delimiter;
        });
        if (out.length===0) {
            window.alert("No entries selected for copying");
            return 0;
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
                    "<p>"+array.length+" ID's copied to clipboard!</p>"
                ].join("\n")));
                powertip_object.powerTip("show");
                setTimeout(function() {
                    powertip_object.data("powertipjq", $([
                        "<p>Export "+data_name+"</p>"
                    ].join("\n")));
                },1000);
            }
        } catch (err) {
            window.alert("Oops, unable to copy");
        }
        document.body.removeChild(textArea);
    }

    $(".uniprot-export1").data("powertipjq", $([
        "<p>Export UniProt IDs</p>"
    ].join("\n")));

    $(".glyphicon-export").powerTip({
        placement: "n",
        smartPlacement: true
    });

    $("#uniprot_copy1").click(function () {
        copyToClipboard($(".alt_selected > .uniprot1 > a"), "\n", "UniProt IDs", $(".uniprot-export1"));
    });

    //Uncheck every row when using back button on browser
    $(".alt_selected").prop("checked",false);
    $(".alt").prop("checked",false);
    $(".select-all").prop("checked",false);


});

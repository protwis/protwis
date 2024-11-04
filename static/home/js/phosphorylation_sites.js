function showColumn(colClass, event) {
    // Hide all columns
    $('.icl2, .icl3, .c-term, .icl2_icl3, .icl2_c-term, .icl3_c-term, .icl2_icl3_c-term').hide();

    const classList = ['icl2', 'icl3', 'c-term']

    if (classList.includes(colClass)) {
        $('.main').attr('colspan', '10')
    }
    else {
        $('.main').attr('colspan', '9')
    }

    // Show the selected column
    $('.' + colClass).show();
}

let oTable = []

// function hideColumns(table, columns) {
//     table.columns(columns).visible(false, true);
//     table.draw();
//   }
  
// function resetHiddenColumns(table) {
// let col_length = table.columns()[0].length;
// let columns = Array.from(new Array(col_length - 1), (x, i) => i);
// table.columns(columns).visible(true, false);
// table.draw();
// }
  
// function reset_tab(table) {
// $("input.yadcf-filter-range-number").val("");
// yadcf.exResetAllFilters(table);
// resetHiddenColumns(table);
// }

// Truncate text content in elements with class 'truncate-js'

$(document).ready(function() {


    $('.table').show() 

    $('.truncate').each(function() {
        let fullText = $(this).text();
        $(this).attr('title', fullText);

        $(this).css({
            'max-width': '4rem', 
            'white-space': 'nowrap',
            'overflow': 'hidden',
            'text-overflow': 'ellipsis'
        });
    });

    oTable = $('#table').DataTable({
        // deferRender: true,
        // scrollY: '50vh',
        // scrollX: true,
        // scrollCollapse: true,
        // scroller: true,
        // paging: false,
        // autoWidth: true,
        // bInfo: true,

        deferRender: true,
        scrollY: '65vh',
        scrollX: true,
        scrollCollapse: true,
        scroller: true,
        paging: false,
        bSortCellsTop: false, //prevent sort arrows going on bottom row
        aaSorting: [],
    });
    let column_filters = [];
    // Selector column
    column_filters = column_filters.concat(createYADCFfilters(0, 1, "none"));
    // Receptor section
    column_filters = column_filters.concat(createYADCFfilters(1, 2, "multi_select", "select2", "", false, "exact", "html", "80px"));
    column_filters = column_filters.concat(createYADCFfilters(3, 1, "multi_select", "select2", "", false, "exact"));
    column_filters = column_filters.concat(createYADCFfilters(4, 1, "multi_select", "select2", "", false, "exact"));
    column_filters = column_filters.concat(createYADCFfilters(71, 3, "range_number", null, ['Min', 'Max'], false))

    yadcf.init(oTable, 
        column_filters,
    {
        cumulative_filtering: false
    });

    // First tab initiation  
    showColumn('icl2');
    // Changing span of phosphorylation cell to match the tabs
    $('.main').attr('colspan', '10')
    oTable.draw()

    // Refresh the table with the current tab
    oTable.on('draw.dt', function() {
        let activeColumn = $('.nav-item.active').attr('data-column');
        console.log(activeColumn , 'Current collumn')
        showColumn(activeColumn)
    });

    $("#table" + " > tbody > tr").click(function(event) {
        if (event.target.type !== "checkbox") {
          $(":checkbox", this).trigger("click");
          $(this).eq(0).toggleClass("alt_selected");
          $(this).find("td").toggleClass("highlight");
        }
        $(this).eq(0).toggleClass("alt_selected");
        $(this).find("td").toggleClass("highlight");
      });
});


$(document).ready(function () {

   $('a[data-toggle="tab"]').on('shown.bs.tab', function (e) {
      console.log( 'show tab' );
      $($.fn.dataTable.tables(true)).DataTable()
          .columns.adjust()
          .responsive.recalc();
   });

   $("#familiestabletab").dataTable({
        "scrollY": $(window).height() - 450,
        "scrollX": true,
        "scrollCollapse": true,
        "scroller": true,
        "paging": false,
        "bSortCellsTop": true, //prevent sort arrows going on bottom row
        "aaSorting": [],
        "autoWidth": true,
        "lengthMenu": [[10, 25, 50, -1], [10, 25, 50, "All"]],
        "pageLength": -1,
        "dom": 'Blfrtip',
        "buttons": ['copy', 'csv', 'excel'],
        "bInfo": true,
   });

   $("#subtypestabletab").dataTable({
      //data: data1,
      "ordering": false,
      "paging": false,
      "fixedHeader":true,
      "searching": false,
      "info": false,
      //"deferRender": true,
   });

   $("#bouviertabletab").dataTable({
      //data: data1,
      "ordering": false,
      "paging": false,
      "fixedHeader":true,
      "searching": false,
      "info": false,
      //"deferRender": true,
   });

   $("#inouetabletab").dataTable({
      //data: data1,
      "ordering": false,
      "paging": false,
      "fixedHeader":true,
      "searching": false,
      "info": false,
      //"deferRender": true,
   });

   // By default display the first tab. If this is not ON, one has to click on the tab for display.
   $('#myTab a:first').tab('show');

});


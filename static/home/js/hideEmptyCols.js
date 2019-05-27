function hideEmptyCols(table) {
  //count # of columns
  var numCols = $("th", table).length;
  for ( var i=1; i<=numCols; i++ ) {
      var empty = true;
      //grab all the <td>'s of the column at i
      $("td:nth-child(" + i + ")", table).each(function(index, el) {
          //check if the <span> of this <td> is empty
          if ( $("td", el).text() != "No data" ) {
              empty = false;
              return false; //break out of each() early
          }
      });
      if ( empty ) {
          $("td:nth-child(" + i + ")", table).hide(); //hide <td>'s
          $("th:nth-child(" + i + ")", table).hide(); //hide header <th>
      }
  }
}

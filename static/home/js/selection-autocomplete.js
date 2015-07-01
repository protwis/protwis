$.widget( "custom.catcomplete", $.ui.autocomplete, {
    _create: function() {
      this._super();
      this.widget().menu( "option", "items", "> :not(.ui-autocomplete-category)" );
    },
    _renderMenu: function( ul, items ) {
      var that = this,
        currentCategory = "";
      $.each( items, function( index, item ) {
        var li;
        if ( item.category != currentCategory ) {
          ul.append( "<li class='ui-autocomplete-category'>" + item.category + "</li>" );
          currentCategory = item.category;
        }
        li = that._renderItemData( ul, item );
        if ( item.category ) {
          li.attr( "aria-label", item.category + " : " + item.label );
        }
      });
    }
  });

$(function() {
    $("#selection-autocomplete").catcomplete({
        source: "/protein/autocomplete?type_of_selection=" + type_of_selection,
        minLength: 2,
        autoFocus: true,
        delay: 500,
        create: function(event, ui) { this.focus();return false; },
        focus: function(event, ui) { return false; },
        select: function(event, ui) {
            AddToSelection(type_of_selection, ui.item['type'], ui.item['id']);
            $( '#selection-autocomplete' ).val('');
            return false;
        }
    }).data("custom-catcomplete")._renderItem = function (ul, item) {
        return $("<li></li>")
        .data("item.autocomplete", item)
        .append("<a>" + item.label + "</a>")
        .appendTo(ul);
    };
});
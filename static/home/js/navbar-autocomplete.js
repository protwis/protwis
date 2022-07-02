$.widget( "custom.catcomplete", $.ui.autocomplete, {
    _create() {
      this._super();
      this.widget().menu( "option", "items", "> :not(.ui-autocomplete-category)" );
    },
    _renderMenu( ul, items ) {
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
    redirect_on_select =''
    $("#nav-selection-autocomplete").catcomplete({
        source: "/protein/autocomplete?type_of_selection=navbar",
        minLength: 2,
        autoFocus: true,
        delay: 500,
        create(event, ui) { this.focus();return false; },
        focus(event, ui) { return false; },
        select(event, ui) {
            $( '#selection-autocomplete' ).val('');
            console.log(ui.item['id'],ui.item['type'],ui.item['label']);
            redirect_url = '/protein/'+ui.item['slug'];
            // redirect if select a target/family to browse
            setTimeout(function(){window.location = redirect_url;}, 1);
            return false;
        }
    }).data("custom-catcomplete")._renderItem = function (ul, item) {
        return $("<li></li>")
        .data("item.autocomplete", item)
        .append("<a>" + item.label + "</a>")
        .appendTo(ul);
    };
});
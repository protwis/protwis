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
var redirect_on_select_new_scope = redirect_on_select;
$(function() {
    $("#selection-autocomplete").catcomplete({
        source: "/protein/autocomplete?selection_only_receptors="+selection_only_receptors+"&type_of_selection=" + type_of_selection,
        minLength: 2,
        autoFocus: true,
        delay: 500,
        create: function(event, ui) { this.focus();return false; },
        focus: function(event, ui) { return false; },
        select: function(event, ui) {
            $( "#selection-autocomplete" ).val("");

            // redirect if select a target/family to browse
            if (type_of_selection === "browse") {
                AddToSelection("targets", ui.item["type"], ui.item["id"]);
                toggleButtonClass("selection-button"); // loading effect on button
                setTimeout(function(){window.location = "/" + ui.item["type"] + "/" + ui.item["slug"];}, 200);

            } else if (type_of_selection === "ginterface") {
                //custom for ginterface
                AddToSelection("targets", ui.item["type"], ui.item["id"]);
                toggleButtonClass("selection-button"); // loading effect on button
                setTimeout(function(){window.location = "/signprot/ginterface/" + ui.item["slug"];}, 200);
            } else if (type_of_selection === "browse_gprot") {
                //custom for ginterface
                AddToSelection("targets", ui.item["type"], ui.item["id"]);
                toggleButtonClass("selection-button"); // loading effect on button
                setTimeout(function(){window.location = "/signprot/" + ui.item["slug"];}, 200);
            } else if (type_of_selection === "ligands") {
                //custom for ligands
                //AddToSelection("targets", ui.item["type"], ui.item["id"]);
                toggleButtonClass("selection-button"); // loading effect on button
                setTimeout(function(){window.location = "/ligand/" + ui.item["id"] + "/info";}, 200);
            } else {
                // add to selection
                AddToSelection(type_of_selection, ui.item["type"], ui.item["id"]);
                // redirect the user if only one target can be selected
                if (type_of_selection === "reference" && redirect_on_select === "True") {
                    toggleButtonClass("selection-button"); // loading effect on button
                    setTimeout(function(){window.location = redirect_url;}, 200);
                }
                if (type_of_selection === "targets" && redirect_on_select_new_scope === "True") {
                    console.log("try to go!");
                    toggleButtonClass("selection-button"); // loading effect on button
                    $("#selection-button").click();
                    if (redirect_url) {
                        setTimeout(function(){window.location = redirect_url;}, 200);
                    }
                }
            }

            return false;
        }
    }).data("custom-catcomplete")._renderItem = function (ul, item) {
        return $("<li></li>")
        .data("item.autocomplete", item)
        .append("<a>" + item.label + "</a>")
        .appendTo(ul);
    };
});

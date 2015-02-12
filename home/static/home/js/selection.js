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

function AddToSelection(selection_type, selection_subtype, selection_id) {
    $.ajax({
        'url': '/common/addtoselection',
        'data': {
            selection_type: selection_type,
            selection_subtype: selection_subtype,
            selection_id: selection_id
        },
        'type': 'GET',
        'success': function(data) {
            $("#selection-" + selection_type).html(data);
        },
    });
}

function RemoveFromSelection(selection_type, selection_subtype, selection_id) {
    $.ajax({
        'url': '/common/removefromselection',
        'data': {
            selection_type: selection_type,
            selection_subtype: selection_subtype,
            selection_id: selection_id
        },
        'type': 'GET',
        'success': function(data) {
            $("#selection-" + selection_type).html(data);
        }
    });
}

function ClearSelection(selection_type) {
    $.ajax({
        'url': '/common/clearselection',
        'data': {
            selection_type: selection_type
        },
        'type': 'GET',
        'success': function(data) {
            $("#selection-" + selection_type).html(data);
        }
    });
}

function ToggleFamilyTreeNode(action, type_of_selection, node_id, tree_indent_level) {
    $.ajax({
        'url': '/common/togglefamilytreenode',
        'data': {
            action: action,
            type_of_selection: type_of_selection,
            node_id: node_id,
            tree_indent_level: tree_indent_level
        },
        'type': 'GET',
        'success': function(data) {
            $("#" + node_id).html(data);
        }
    });
}
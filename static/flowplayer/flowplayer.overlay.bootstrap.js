(function() {
  /* global flowplayer */

  var $ = window.jQuery;

  flowplayer.overlay.bootstrap = function(api, root) {
    var conf = api.conf.overlay
      , trigger = conf.trigger
      , modalSize = conf.size
      , modalTitle = conf.title
      , modal = $(
      '<div class="modal">'
        + '<div class="modal-dialog ' + (modalSize ? 'modal-' + modalSize : '') +  '">'
        + '<div class="modal-content">'
        + '<div class="modal-header">'
        + '<button type="button" class="close" data-dismiss="modal"><span>&times;</span></button>'
        + '<h4 class="modal-title">' + (modalTitle ? modalTitle : 'Video') + '</h4>'
        + '</div>'
        + '<div class="modal-body"></div>'
        + '</div>'
        + '</div>'
        + '</div>'
      ).appendTo('body');

    modal.find('.modal-body').append(root);

    if (conf.keyboard !== false) {
      $(document).on('keydown', function(ev) {
        if (ev.keyCode === 27) api.unload();
      });
    }

    modal.on('shown.bs.modal', function() {
      api.load();
      $(root).addClass('is-open');
    }).on('hidden.bs.modal', function() {
      api.unload();
      $(root).removeClass('is-open');
    });

    api.on('load ready', function(e, api, video) {
      if (!modalTitle) {
        if (/l/.test(e.type)) {
          modal.find('.modal-title').text(video.title);
        } else {
          $('.fp-title', root).hide();
        }
      }
    }).on('unload', function() {
      modal.modal('hide');
    });

    $(trigger).on('click', function() {
      modal.modal();
    });

  };
})();
(function () {
  var doctype = '<?xml version="1.0" standalone="no"?>\n' +
    '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" ' +
    '"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">';

  function isExternal(url) {
    return url && url.startsWith('http') && !url.includes(window.location.host);
  }

  function inlineImages(el, callback) {
    var images = el.querySelectorAll('image');
    var left = images.length;
    if (left === 0) return callback();

    for (var i = 0; i < images.length; i++) {
      (function (image) {
        var href = image.getAttributeNS("http://www.w3.org/1999/xlink", "href") || image.getAttribute("href");
        if (isExternal(href)) {
          console.warn("External image skipped:", href);
          left--;
          if (left === 0) callback();
          return;
        }

        var img = new Image();
        img.src = href;
        img.onload = function () {
          var canvas = document.createElement("canvas");
          canvas.width = img.width;
          canvas.height = img.height;
          canvas.getContext("2d").drawImage(img, 0, 0);
          image.setAttributeNS("http://www.w3.org/1999/xlink", "href", canvas.toDataURL("image/png"));
          left--;
          if (left === 0) callback();
        };
        img.onerror = function () {
          console.warn("Image failed to load:", href);
          left--;
          if (left === 0) callback();
        };
      })(images[i]);
    }
  }

  function styles(el) {
    var css = "", sheets = document.styleSheets;
    for (var i = 0; i < sheets.length; i++) {
      var rules;
      try {
        rules = sheets[i].cssRules;
      } catch (e) {
        continue;
      }
      if (!rules) continue;

      for (var j = 0; j < rules.length; j++) {
        var rule = rules[j];
        if (rule instanceof CSSStyleRule && el.querySelector(rule.selectorText)) {
          css += rule.selectorText + " { " + rule.style.cssText + " }\n";
        } else if (rule.cssText && rule.cssText.startsWith("@font-face")) {
          css += rule.cssText + "\n";
        }
      }
    }
    return css;
  }

  window.saveSvgAsPngTree = function (el, name, options) {
    options = options || {};
    var xmlns = "http://www.w3.org/2000/xmlns/";
    var scale = options.scale || 2;

    inlineImages(el, function () {
      var clone = el.cloneNode(true);
      var viewBox = clone.getAttribute("viewBox");
      if (!viewBox) {
        console.error("SVG must have a viewBox attribute.");
        return;
      }

      // Parse viewBox
      var [x, y, width, height] = viewBox.split(" ").map(parseFloat);

      // Add padding to avoid cropping
      var pad = 50;
      var newX = x - pad;
      var newY = y - pad;
      var newWidth = width + pad * 2;
      var newHeight = height + pad * 2;

      clone.setAttribute("version", "1.1");
      clone.setAttributeNS(xmlns, "xmlns", "http://www.w3.org/2000/svg");
      clone.setAttributeNS(xmlns, "xmlns:xlink", "http://www.w3.org/1999/xlink");
      clone.setAttribute("width", newWidth * scale);
      clone.setAttribute("height", newHeight * scale);
      clone.setAttribute("viewBox", `${newX} ${newY} ${newWidth} ${newHeight}`);

      var outer = document.createElement("div");
      outer.appendChild(clone);

      var css = styles(el);
      var styleEl = document.createElement("style");
      styleEl.setAttribute("type", "text/css");
      styleEl.innerHTML = "<![CDATA[\n" + css + "\n]]>";
      var defs = document.createElement("defs");
      defs.appendChild(styleEl);
      clone.insertBefore(defs, clone.firstChild);

      var svgStr = doctype + outer.innerHTML;
      var uri = "data:image/svg+xml;base64," + window.btoa(unescape(encodeURIComponent(svgStr)));

      var img = new Image();
      img.onload = function () {
        var canvas = document.createElement("canvas");
        canvas.width = img.width;
        canvas.height = img.height;
        canvas.getContext("2d").drawImage(img, 0, 0);

        var a = document.createElement("a");
        a.download = name;
        a.href = canvas.toDataURL("image/png");
        document.body.appendChild(a);
        a.click();
        a.remove();
      };
      img.src = uri;
    });
  };
})();

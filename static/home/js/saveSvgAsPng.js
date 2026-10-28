(function() {
  var out$ = typeof exports != 'undefined' && exports || this;

  var doctype = '<?xml version="1.0" standalone="no"?><!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">';

  function isExternal(url) {
    return url && url.lastIndexOf('http',0) == 0 && url.lastIndexOf(window.location.host) == -1;
  }

  function inlineImages(el, callback) {
    var images = el.querySelectorAll('image');
    var left = images.length;
    if (left == 0) {
      callback();
    }
    for (var i = 0; i < images.length; i++) {
      (function(image) {
        var href = image.getAttributeNS("http://www.w3.org/1999/xlink", "href");
        if (href) {
          if (isExternal(href.value)) {
            console.warn("Cannot render embedded images linking to external hosts: "+href.value);
            return;
          }
        }
        var canvas = document.createElement('canvas');
        var ctx = canvas.getContext('2d');
        var img = new Image();
        href = href || image.getAttribute('href');
        img.src = href;
        img.onload = function() {
          canvas.width = img.width;
          canvas.height = img.height;
          ctx.drawImage(img, 0, 0);
          image.setAttributeNS("http://www.w3.org/1999/xlink", "href", canvas.toDataURL('image/png'));
          left--;
          if (left == 0) {
            callback();
          }
        }
        img.onerror = function() {
          console.log("Could not load "+href);
          left--;
          if (left == 0) {
            callback();
          }
        }
      })(images[i]);
    }
  }

  function styles(el, selectorRemap) {
    var css = "";
    var sheets = document.styleSheets;
    for (var i = 0; i < sheets.length; i++) {
      if (isExternal(sheets[i].href)) {
        console.warn("Cannot include styles from other hosts: "+sheets[i].href);
        continue;
      }
      var rules = sheets[i].cssRules;
      if (rules != null) {
        for (var j = 0; j < rules.length; j++) {
          var rule = rules[j];
          if (typeof(rule.style) != "undefined") {
            var match = null;
            try {
              match = el.querySelector(rule.selectorText);
            } catch(err) {
              console.warn('Invalid CSS selector "' + rule.selectorText + '"', err);
            }
            if (match) {
              var selector = selectorRemap ? selectorRemap(rule.selectorText) : rule.selectorText;
              css += selector + " { " + rule.style.cssText + " }\n";
            } else if(rule.cssText.match(/^@font-face/)) {
              css += rule.cssText + '\n';
            }
          }
        }
      }
    }
    return css;
  }

  // Illustrator, Inkscape, PowerPoint etc. ignore dominant-/alignment-baseline and often the
  // embedded CSS too, so bake what the browser actually renders into the clone: elements hidden
  // on screen (e.g. the snake plot's ".long" loops) are dropped, and baseline offsets become
  // plain translates. The PNG/JPG/TIFF exports render this same clone, so they stay identical.
  var BASELINE_PROPS = ['dominant-baseline', 'alignment-baseline'];
  var GRAPHIC_TAG = /^(g|text|circle|rect|path|line|polyline|polygon|ellipse|image|use)$/i;
  var NON_RENDERED_ANCESTOR = 'defs, clipPath, mask, marker, pattern, symbol';

  function bboxY(node) {
    try { return node.getBBox().y; } catch (err) { return NaN; }  // Firefox throws if not rendered
  }

  function bakeForExport(el, clone) {
    var liveTexts = el.querySelectorAll('text');
    var cloneTexts = clone.querySelectorAll('text');
    var shifted = [];
    for (var i = 0; i < liveTexts.length; i++) {
      var cs = out$.getComputedStyle(liveTexts[i]);
      if (cs.display == 'none') continue;
      var custom = BASELINE_PROPS.some(function(p) {
        var v = cs.getPropertyValue(p);
        return v && v != 'auto' && v != 'baseline';
      });
      if (custom) shifted.push({live: liveTexts[i], copy: cloneTexts[i], y: bboxY(liveTexts[i])});
    }

    // Measure everything with the baseline reset in one pass (a single relayout), then restore.
    shifted.forEach(function(s) {
      s.saved = BASELINE_PROPS.map(function(p) {
        return [p, s.live.style.getPropertyValue(p), s.live.style.getPropertyPriority(p)];
      });
      BASELINE_PROPS.forEach(function(p) { s.live.style.setProperty(p, 'auto', 'important'); });
    });
    shifted.forEach(function(s) { s.shift = s.y - bboxY(s.live); });
    shifted.forEach(function(s) {
      s.saved.forEach(function(v) {
        if (v[1]) s.live.style.setProperty(v[0], v[1], v[2]);
        else s.live.style.removeProperty(v[0]);
      });
    });

    // No measurable shift (e.g. the svg isn't on screen) -- leave that text untouched.
    shifted.forEach(function(s) {
      if (!s.copy || !(Math.abs(s.shift) > 0.01)) return;
      BASELINE_PROPS.forEach(function(p) {
        s.copy.removeAttribute(p);
        s.copy.style.setProperty(p, 'auto');  // beats any exported CSS rule in browsers
      });
      var t = s.copy.getAttribute('transform');
      s.copy.setAttribute('transform', (t ? t + ' ' : '') + 'translate(0,' + Math.round(s.shift * 100) / 100 + ')');
    });

    var liveAll = el.querySelectorAll('*');
    var cloneAll = clone.querySelectorAll('*');
    var hidden = [];
    for (var j = 0; j < liveAll.length; j++) {
      var node = liveAll[j];
      if (GRAPHIC_TAG.test(node.tagName) && !node.closest(NON_RENDERED_ANCESTOR) &&
          out$.getComputedStyle(node).display == 'none') {
        hidden.push(cloneAll[j]);
      }
    }
    hidden.forEach(function(n) { if (n.parentNode) n.parentNode.removeChild(n); });
  }

  out$.svgAsDataUri = function(el, options, cb) {
    options = options || {};
    options.scale = options.scale || 1.6;
    var xmlns = "http://www.w3.org/2000/xmlns/";

    inlineImages(el, function() {
      var outer = document.createElement("div");
      var clone = el.cloneNode(true);
      bakeForExport(el, clone);
      var width, height;
      if(el.tagName == 'svg') {
        var box = el.getBoundingClientRect();
        width = parseInt(clone.getAttribute('width') ||
          box.width ||
          clone.style.width ||
          out$.getComputedStyle(el).getPropertyValue('width'));
        height = parseInt(clone.getAttribute('height') ||
          box.height ||
          clone.style.height ||
          out$.getComputedStyle(el).getPropertyValue('height'));

        if (clone.getAttribute('viewBox')) {
          width = clone.getAttribute('viewBox').split(" ")[2];
          height = clone.getAttribute('viewBox').split(" ")[3];
        }


        if (width === undefined ||
            width === null ||
            isNaN(parseFloat(width))) {
          width = 0;
        }
        if (height === undefined ||
            height === null ||
            isNaN(parseFloat(height))) {
          height = 0;
        }
      } else {
        var box = el.getBBox();
        width = box.x + box.width;
        height = box.y + box.height;
        clone.setAttribute('transform', clone.getAttribute('transform').replace(/translate\(.*?\)/, ''));

        var svg = document.createElementNS('http://www.w3.org/2000/svg','svg')
        svg.appendChild(clone)
        clone = svg;
      }

      clone.setAttribute("version", "1.1");
      clone.setAttributeNS(xmlns, "xmlns", "http://www.w3.org/2000/svg");
      clone.setAttributeNS(xmlns, "xmlns:xlink", "http://www.w3.org/1999/xlink");
      clone.setAttribute("width", width * options.scale);
      clone.setAttribute("height", height * options.scale);
      clone.setAttribute("viewBox", "0 0 " + width + " " + height);

      outer.appendChild(clone);

      var css = styles(el, options.selectorRemap);
      var s = document.createElement('style');
      s.setAttribute('type', 'text/css');
      s.innerHTML = "<![CDATA[\n" + css + "\n]]>";
      var defs = document.createElement('defs');
      defs.appendChild(s);
      clone.insertBefore(defs, clone.firstChild);

      var svg = doctype + outer.innerHTML;
      var uri = 'data:image/svg+xml;base64,' + window.btoa(unescape(encodeURIComponent(svg)));
      if (cb) {
        cb(uri);
      }
    });
  }

  out$.saveSvgAsPng = function(el, name, options) {
    options = options || {};
    out$.svgAsDataUri(el, options, function(uri) {
      var image = new Image();
      image.onload = function() {
        var canvas = document.createElement('canvas');
        canvas.width = image.width;
        canvas.height = image.height;
        var context = canvas.getContext('2d');
        context.drawImage(image, 0, 0);

        var a = document.createElement('a');
        a.download = name;
        a.href = canvas.toDataURL('image/png');
        document.body.appendChild(a);
        a.addEventListener("click", function(e) {
          a.parentNode.removeChild(a);
        });
        a.click();
      }
      image.src = uri;
    });
  }

  out$.saveSvgAsJpg = function(el, name, options) {
    options = options || {};
    out$.svgAsDataUri(el, options, function(uri) {
      var image = new Image();
      image.onload = function() {
        var canvas = document.createElement('canvas');
        canvas.width = image.width;
        canvas.height = image.height;
        var context = canvas.getContext('2d');
        context.drawImage(image, 0, 0);

        var a = document.createElement('a');
        a.download = name;
        a.href = canvas.toDataURL('image/jpg');
        document.body.appendChild(a);
        a.addEventListener("click", function(e) {
          a.parentNode.removeChild(a);
        });
        a.click();
      }
      image.src = uri;
    });
  }

  out$.saveSvgAsTiff = function(el, name, options) {
    options = options || {};
    out$.svgAsDataUri(el, options, function(uri) {
      var image = new Image();
      image.onload = function() {
        var canvas = document.createElement('canvas');
        canvas.width = image.width;
        canvas.height = image.height;
        var context = canvas.getContext('2d');
        context.drawImage(image, 0, 0);

        var a = document.createElement('a');
        a.download = name;
        a.href = canvas.toDataURL('image/tiff');
        document.body.appendChild(a);
        a.addEventListener("click", function(e) {
          a.parentNode.removeChild(a);
        });
        a.click();
      }
      image.src = uri;
    });
  }
})();

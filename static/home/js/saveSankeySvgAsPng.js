(function() {
  var out$ = typeof exports !== 'undefined' && exports || this;

  // Minimal SVG doctype
  var doctype = '<?xml version="1.0" standalone="no"?>'
    + '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" '
    + '"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">';

  function isExternal(url) {
    return url && url.indexOf('http') === 0 && url.indexOf(window.location.host) === -1;
  }

  function inlineImages(el, callback) {
    const images = el.querySelectorAll('image');
    let remaining = images.length;
    if (remaining === 0) {
      callback();
      return;
    }
    for (let i = 0; i < images.length; i++) {
      (function(image) {
        let href = image.getAttributeNS("http://www.w3.org/1999/xlink", "href") ||
                   image.getAttribute('href');
        if (href && isExternal(href)) {
          console.warn("Skipping external image:", href);
          remaining--; // Decrement separately
          if (remaining === 0) {
              callback();
          }
          return;
        }
        const canvas = document.createElement('canvas');
        const ctx = canvas.getContext('2d');
        const img = new Image();
        img.src = href;
        img.crossOrigin = 'anonymous';
        img.onload = function() {
          canvas.width = img.width;
          canvas.height = img.height;
          ctx.drawImage(img, 0, 0);
          const png = canvas.toDataURL('image/png');
          image.setAttributeNS("http://www.w3.org/1999/xlink", "href", png);
          remaining--; // Decrement separately
          if (remaining === 0) {
              callback();
          }
        };
        img.onerror = function() {
          console.warn("Could not load image at " + href);
          remaining--; // Decrement separately
          if (remaining === 0) {
              callback();
          }
        };
      })(images[i]);
    }
  }

  function inlineStyles(el) {
    let css = "";
    const sheets = Array.from(document.styleSheets);

    sheets.forEach((sheet) => {
        // Try-catch block to avoid SecurityErrors
        try {
            if (sheet.cssRules) {
                const rules = Array.from(sheet.cssRules);

                rules.forEach((rule) => {
                    // Only include styles that are relevant to the SVG
                    if (rule instanceof CSSStyleRule && el.querySelector(rule.selectorText)) {
                        css += rule.cssText + "\n";
                    }
                    // Include @font-face and other at-rules
                    else if (rule instanceof CSSFontFaceRule || rule.cssText.match(/^@font-face/)) {
                        css += rule.cssText + "\n";
                    }
                });
            }
        } catch (e) {
            console.warn("Skipping cross-origin stylesheet:", sheet.href, e);
        }
    });

    return css;
}


  /**
   * 1. Clone the original <svg>
   * 2. Inline images & CSS
   * 3. Return a data URI (base64) for the resulting SVG
   */
  out$.svgAsDataUriExact = function(originalSvg, callback) {
    // 1) Clone
    const clone = originalSvg.cloneNode(true);

    // 2) Convert local <image> references to data URIs
    inlineImages(clone, function() {
      // gather matching CSS
      let css = inlineStyles(clone);

      // Example: forcibly add link styling if you want

      if (css.trim()) {
        const styleElem = document.createElement('style');
        styleElem.setAttribute('type', 'text/css');
        styleElem.textContent = css;
        const defs = document.createElement('defs');
        defs.appendChild(styleElem);
        clone.insertBefore(defs, clone.firstChild);
      }

      // 3) Serialize
      const serializer = new XMLSerializer();
      let rawSvg = serializer.serializeToString(clone);
      rawSvg = doctype + rawSvg;

      // 4) Base64
      let base64;
      try {
        base64 = btoa(unescape(encodeURIComponent(rawSvg)));
      } catch (e) {
        console.warn("Falling back to direct btoa. Some chars may break. Error:", e);
        base64 = btoa(rawSvg);
      }
      const uri = 'data:image/svg+xml;base64,' + base64;
      callback(uri);
    });
  };

  /**
   * Helper: given an <svg>, produce a <canvas> image, then download it as PNG
   */
  out$.saveSankeySvgAsPng = function(svgEl, fileName, scaleFactor = 4) {
    out$.svgAsDataUriExact(svgEl, function(svgUri) {
      const img = new Image();
      img.crossOrigin = 'anonymous';
      img.onload = function() {
        // Use naturalWidth / naturalHeight rather than a second canvas:
        const originalWidth = img.naturalWidth;
        const originalHeight = img.naturalHeight;
        
        // Scale up
        const canvas = document.createElement('canvas');
        canvas.width = originalWidth * scaleFactor;
        canvas.height = originalHeight * scaleFactor;
  
        // High-quality rendering
        const ctx = canvas.getContext('2d');
        ctx.imageSmoothingEnabled = true;
        ctx.imageSmoothingQuality = 'high';
  
        // Draw the image scaled up
        ctx.drawImage(
          img,
          0, 0, originalWidth, originalHeight,   // source
          0, 0, canvas.width, canvas.height     // destination
        );
  
        // Convert canvas to PNG and download
        const link = document.createElement('a');
        link.download = fileName || 'sankey.png';
        link.href = canvas.toDataURL('image/png');
        document.body.appendChild(link);
        link.addEventListener('click', () => {
          link.parentNode.removeChild(link);
        });
        link.click();
      };
      img.onerror = err => {
        console.error('Could not load SVG into an image', err);
      };
      img.src = svgUri;
    });
  };
  

  /* ====================== ADDITIONAL EXPORT METHODS BELOW ======================== */

  /**
   * Save as raw .SVG (text). This will simply download the inlined SVG itself.
   * It uses the same svgAsDataUriExact, which returns a "data:image/svg+xml;base64,..."
   */
  out$.saveSankeySvgAsSVG = function(svgEl, fileName) {
    out$.svgAsDataUriExact(svgEl, function(uri) {
      const link = document.createElement('a');
      link.download = fileName || 'sankey.svg';
      link.href = uri;
      document.body.appendChild(link);
      link.addEventListener('click', () => {
        link.parentNode.removeChild(link);
      });
      link.click();
    });
  };

  /**
   * Save as .JPG using the Canvas "image/jpeg" format.
   * Note you can pass a second parameter for quality, e.g. toDataURL('image/jpeg', 0.9).
   */
  out$.saveSankeySvgAsJpg = function(svgEl, fileName) {
    out$.svgAsDataUriExact(svgEl, function(uri) {
      const image = new Image();
      image.crossOrigin = 'anonymous';
      image.onload = function() {
        const canvas = document.createElement('canvas');
        canvas.width = image.width;
        canvas.height = image.height;
        const ctx = canvas.getContext('2d');
        ctx.drawImage(image, 0, 0);

        const link = document.createElement('a');
        link.download = fileName || 'sankey.jpg';
        // The second parameter (0.92) is the quality for JPEG
        link.href = canvas.toDataURL('image/jpeg', 0.92);
        document.body.appendChild(link);
        link.addEventListener('click', () => {
          link.parentNode.removeChild(link);
        });
        link.click();
      };
      image.onerror = function(err) {
        console.error("Could not load SVG data into an image. Error:", err);
      };
      image.src = uri;
    });
  };

  /**
   * Save as .TIFF - Requires UTIF.js (https://github.com/photopea/UTIF.js)
   * If UTIF is not loaded, this will fail. 
   * Make sure to include UTIF.js before this script:
   * <script src="path/to/UTIF.js"></script>
   */
  out$.saveSankeySvgAsTiff = function(svgEl, fileName) {
    if (typeof UTIF === 'undefined') {
      console.error("UTIF.js not found. Cannot save as TIFF.");
      return;
    }
    out$.svgAsDataUriExact(svgEl, function(uri) {
      const image = new Image();
      image.crossOrigin = 'anonymous';
      image.onload = function() {
        const canvas = document.createElement('canvas');
        canvas.width = image.width;
        canvas.height = image.height;
        const ctx = canvas.getContext('2d');
        ctx.drawImage(image, 0, 0);

        // Get the RGBA data from canvas
        const imageData = ctx.getImageData(0, 0, canvas.width, canvas.height);
        // UTIF.encode takes raw RGBA, width, height
        const tiffBuffer = UTIF.encode(imageData.data, canvas.width, canvas.height);

        // Turn the ArrayBuffer returned by UTIF into a Blob
        const blob = new Blob([tiffBuffer], { type: 'image/tiff' });

        // Create a download link
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.download = fileName || 'sankey.tiff';
        link.href = url;
        document.body.appendChild(link);
        link.addEventListener('click', () => {
          link.parentNode.removeChild(link);
          URL.revokeObjectURL(url);
        });
        link.click();
      };
      image.onerror = function(err) {
        console.error("Could not load SVG data into an image. Error:", err);
      };
      image.src = uri;
    });
  };

})();

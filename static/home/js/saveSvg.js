/**
 * Collect page CSS rules that match any element inside the given SVG element.
 * Returns a string of CSS rules safe to embed in an SVG <style> tag.
 */
function collectSvgCss(svgEl) {
    var css = '';
    var sheets = document.styleSheets;
    for (var i = 0; i < sheets.length; i++) {
        var rules;
        try { rules = sheets[i].cssRules; } catch (e) { continue; }
        if (!rules) { continue; }
        for (var j = 0; j < rules.length; j++) {
            var rule = rules[j];
            try {
                if (rule instanceof CSSStyleRule && svgEl.querySelector(rule.selectorText)) {
                    css += rule.selectorText + ' { ' + rule.style.cssText + ' }\n';
                } else if (rule.cssText && rule.cssText.startsWith('@font-face')) {
                    css += rule.cssText + '\n';
                }
            } catch (e) { /* skip selectors that querySelector can't parse */ }
        }
    }
    return css;
}

function saveSvg(svgEl, name) {
    // cloneNode(true) copies all D3-applied inline styles (stroke, fill, etc.) already.
    // Do NOT call inlineStyles() on a detached clone — getComputedStyle returns empty
    // values for detached elements and would overwrite those styles with blank strings.
    var clone = svgEl.cloneNode(true);

    // Embed any CSS-class-based rules as a <style> block so the file is self-contained.
    var svgNS = 'http://www.w3.org/2000/svg';
    var css = collectSvgCss(svgEl);
    if (css) {
        var defs = clone.querySelector('defs');
        if (!defs) {
            defs = document.createElementNS(svgNS, 'defs');
            clone.insertBefore(defs, clone.firstChild);
        }
        var styleEl = document.createElementNS(svgNS, 'style');
        styleEl.setAttribute('type', 'text/css');
        styleEl.textContent = css;
        defs.insertBefore(styleEl, defs.firstChild);
    }

    clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    clone.setAttribute('xmlns:xlink', 'http://www.w3.org/1999/xlink');

    var svgData = clone.outerHTML;
    var preface = '<?xml version="1.0" standalone="no"?>\r\n';
    var svgBlob = new Blob([preface, svgData], { type: 'image/svg+xml;charset=utf-8' });
    var svgUrl = URL.createObjectURL(svgBlob);

    var a = document.createElement('a');
    a.href = svgUrl;
    a.download = name;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(svgUrl);
}

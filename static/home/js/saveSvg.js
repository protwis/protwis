function inlineStyles(svgEl) {
    // Get all elements within the SVG
    const elements = svgEl.querySelectorAll('*');
    
    // Loop through each element and inline its computed styles
    elements.forEach(el => {
        const computedStyle = window.getComputedStyle(el);
        let styleString = '';
        
        // Copy each computed style to a string
        for (let i = 0; i < computedStyle.length; i++) {
            const key = computedStyle[i];
            const value = computedStyle.getPropertyValue(key);
            styleString += `${key}:${value};`;
        }
        
        // Set the inline style for the element
        el.setAttribute('style', styleString);
    });
}

function saveSvg(svgEl, name) {
    // Ensure styles are inlined before saving
    inlineStyles(svgEl);
    
    // Add necessary namespaces for SVG
    svgEl.setAttribute("xmlns", "http://www.w3.org/2000/svg");
    svgEl.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink");

    // Serialize the SVG and create a Blob for downloading
    var svgData = svgEl.outerHTML;
    var preface = '<?xml version="1.0" standalone="no"?>\r\n';
    var svgBlob = new Blob([preface, svgData], {type: "image/svg+xml;charset=utf-8"});
    var svgUrl = URL.createObjectURL(svgBlob);

    // Create a download link and trigger the download
    var downloadLink = document.createElement("a");
    downloadLink.href = svgUrl;
    downloadLink.download = name;
    document.body.appendChild(downloadLink);
    downloadLink.click();
    document.body.removeChild(downloadLink);
}

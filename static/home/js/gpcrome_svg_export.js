// Shared SVG/PNG export helpers for GPCRome-style D3 plots (classification trees, the
// classification wheel, and the structure/statistics wheels). Deliberately not nested under
// classification/ -- non-classification pages depend on it too.
//
// Kept minimal on purpose: a plain white-background PNG rasterize and a plain SVG-file
// download. Anything with real per-page behavior (combining multiple SVGs, picking which SVG
// is "active", building the download filename, etc.) stays in the page's own script.
var GPCRomeSvgExport = (function () {
    function downloadSvg(svgElement, filename) {
        if (!svgElement) return;
        const clone = svgElement.cloneNode(true);
        clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
        clone.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink");
        const serializer = new XMLSerializer();
        const blob = new Blob([serializer.serializeToString(clone)], { type: "image/svg+xml;charset=utf-8" });
        const url = URL.createObjectURL(blob);
        const a = document.createElement("a");
        a.href = url;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        setTimeout(function () { URL.revokeObjectURL(url); }, 1000);
    }

    function downloadSvgAsPng(svgElement, filename, scale) {
        if (!svgElement) return;
        const s = scale || 2;

        // Clone and ensure namespaces
        const svgClone = svgElement.cloneNode(true);
        svgClone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
        svgClone.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink");

        const width = parseFloat(svgClone.getAttribute("width")) || 800;
        const height = parseFloat(svgClone.getAttribute("height")) || 800;

        // Serialize
        const serializer = new XMLSerializer();
        const svgString = serializer.serializeToString(svgClone);
        const svgBlob = new Blob([svgString], { type: "image/svg+xml;charset=utf-8" });
        const url = URL.createObjectURL(svgBlob);

        const img = new Image();
        img.onload = function () {
            const canvas = document.createElement("canvas");
            canvas.width = Math.round(width * s);
            canvas.height = Math.round(height * s);
            const ctx = canvas.getContext("2d");

            // White background
            ctx.fillStyle = "#ffffff";
            ctx.fillRect(0, 0, canvas.width, canvas.height);

            ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
            URL.revokeObjectURL(url);

            const pngUrl = canvas.toDataURL("image/png");
            const a = document.createElement("a");
            a.href = pngUrl;
            a.download = filename;
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
        };
        img.onerror = function () {
            URL.revokeObjectURL(url);
        };
        img.src = url;
    }

    return {
        downloadSvg: downloadSvg,
        downloadSvgAsPng: downloadSvgAsPng,
    };
})();

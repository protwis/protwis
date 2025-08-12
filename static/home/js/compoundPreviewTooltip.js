/**
 * compoundPreviewTooltip.js (v2.4)
 *
 * This module provides a tooltip for chemical structures.
 * - For Peptides: Uses PeptideLayoutEngine if sequence is standard and engine is available.
 * - For Small Molecules (and others):
 *   - If `data-picture` or `rel` is "Generate_from_SMILES" AND `data-smiles` and OCL are available, generates SVG.
 *   - Else if `data-picture` or `rel` is a valid image URL, displays that image.
 *   - Else, falls back to a default "No image available" picture.
 *
 * Usage:
 *   compoundPreview("a.struct");
 *
 * Dependencies:
 *   - jQuery
 *   - OpenChemLib (OCL) - must be loaded if SMILES rendering is needed.
 *   - PeptideLayoutEngine - must be loaded if peptide rendering is needed.
 *
 * Note:
 *   - `window.NO_IMAGE_URL` should be defined.
 *   - Data attributes for peptides: `data-ligand-type="peptide" data-sequence="SEQUENCE_STRING"`
 *   - Data attributes for small molecules:
 *     - SMILES rendering: `data-smiles="SMILES_STRING" [data-picture="Generate_from_SMILES" or rel="Generate_from_SMILES"]`
 *     - Pre-rendered image: `[data-picture="IMAGE_URL" or rel="IMAGE_URL"]`
 */
function compoundPreview(selector = "a.struct", noImageUrl = window.NO_IMAGE_URL) {
    if (!noImageUrl) {
        noImageUrl = "/static/home/images/No_image_available.svg";
    }
    const xOffset = 30, yOffset = -10; // Original yOffset makes it appear slightly above cursor line.

    function isStandardPeptideSequence(seqStr) {
        if (!seqStr) return false;
        return /^[ABCDEFGHIKLMNOPQRSTUVWXYZ\s]+$/i.test(seqStr);
    }

    function calculateAndSetPosition(e, $tooltip) {
        if (!$tooltip || !$tooltip.length) return;

        const tooltipWidth = $tooltip.outerWidth();
        const tooltipHeight = $tooltip.outerHeight();
        const windowWidth = $(window).width();
        const windowHeight = $(window).height();
        // e.pageX/pageY are document coordinates.
        // e.clientX/clientY are viewport coordinates.

        let newLeft = e.pageX + xOffset;
        // Original logic for yOffset: e.pageY - yOffset means e.pageY - (-10) = e.pageY + 10
        // So the top of the tooltip is 10px BELOW the mouse y-coordinate.
        let newTop = e.pageY - yOffset;

        // Check right boundary (using clientX for viewport-relative check)
        if (e.clientX + xOffset + tooltipWidth > windowWidth) {
            newLeft = e.pageX - tooltipWidth - xOffset; // Place to the left of cursor
        }

        // Check bottom boundary (using clientY for viewport-relative check)
        // e.clientY - yOffset + tooltipHeight: e.clientY - (-10) + tooltipHeight = e.clientY + 10 + tooltipHeight
        if (e.clientY - yOffset + tooltipHeight > windowHeight) {
            // Place above the cursor.
            // newTop becomes: e.pageY (cursor) + yOffset (push up by 10) - tooltipHeight
            newTop = e.pageY + yOffset - tooltipHeight;
        }

        // Prevent tooltip from going off the top or left of the document after adjustments
        // (though less likely with this logic, good for robustness)
        if (newLeft < $(window).scrollLeft()) {
            newLeft = $(window).scrollLeft();
        }
        if (newTop < $(window).scrollTop()) {
            newTop = $(window).scrollTop();
        }

        $tooltip.css({ top: newTop + "px", left: newLeft + "px" });
    }

    $('body').on('mouseenter', selector, function(e) {
        let $link = $(this);
        let tooltipHtml = "";

        const ligandType = $link.data('ligand-type');
        const sequence = $link.data('sequence');

        if (ligandType === 'peptide') {
            if (isStandardPeptideSequence(sequence) && typeof PeptideLayoutEngine !== 'undefined') {
                try {
                    const peptideDrawer = new PeptideLayoutEngine({
                        residueRadius: 16, overlapAmount: 6, padding: 10,
                        fontSize: 12, residuesPerRow: 10, maxLength: 50,
                        rowSpacingFactor: 1.5
                    });
                    const cleanedSequence = sequence.replace(/\s+/g, '').toUpperCase();
                    if (cleanedSequence) {
                        const svgElement = peptideDrawer.generateSvgElement(cleanedSequence);
                        if (svgElement) {
                            tooltipHtml = svgElement.outerHTML;
                        } else {
                            console.warn("PeptideLayoutEngine.generateSvgElement returned null for sequence:", cleanedSequence);
                            tooltipHtml = `<img style="max-width: 200px; height: auto" src="${noImageUrl}" alt="Preview not available" />`;
                        }
                    } else {
                        console.warn("Peptide sequence became empty after cleaning:", sequence);
                        tooltipHtml = `<img style="max-width: 200px; height: auto" src="${noImageUrl}" alt="Preview not available" />`;
                    }
                } catch(err) {
                    console.error("PeptideLayoutEngine error in tooltip:", err, "Original sequence:", sequence);
                    tooltipHtml = `<img style="max-width: 200px; height: auto" src="${noImageUrl}" alt="Preview not available" />`;
                }
            } else {
                // Fallback for non-standard peptide or missing engine
                if (sequence && !isStandardPeptideSequence(sequence)) {
                     console.log("Tooltip: Peptide sequence non-standard, falling back. Sequence:", sequence);
                } else if (!sequence) {
                     console.log("Tooltip: Peptide sequence data attribute is missing or empty.");
                } else if (typeof PeptideLayoutEngine === 'undefined') {
                     console.warn("Tooltip: PeptideLayoutEngine library not found for peptide rendering.");
                }
                tooltipHtml = `<img style="max-width: 200px; height: auto" src="${noImageUrl}" alt="Preview not available" />`;
            }
        } else { // Not 'peptide', so treat as small molecule or other
            let smiles = $link.data("smiles") || "";
            // Prioritize data-picture, then rel, for pictureSource info
            let pictureSource = $link.data("picture") || $link.attr("rel");

            if (pictureSource === "Generate_from_SMILES" && smiles && typeof OCL !== 'undefined') {
                // Explicit request to generate from SMILES (similar to v2.0's primary SMILES condition)
                try {
                    let mol = OCL.Molecule.fromSmiles(smiles);
                    mol.inventCoordinates();
                    const svgSize = 200;
                    const svgOptions = {
                        autoCrop: true, noStereoProblem: true,
                        suppressChiralText: true, suppressCIPParity: true, suppressESR: true,
                    };
                    tooltipHtml = mol.toSVG(svgSize, svgSize, null, svgOptions);
                } catch (err) {
                    console.error("OCL error in tooltip (explicit SMILES gen):", err, "SMILES:", smiles);
                    tooltipHtml = `<img style="max-width: 200px; height: auto" src="${noImageUrl}" alt="Preview not available" />`; // Changed from "Invalid SMILES" div to image
                }
            } else if (pictureSource && pictureSource !== "Not_available" && pictureSource !== "Generate_from_SMILES") {
                // If pictureSource is a valid image URL (not "Generate_from_SMILES" or "Not_available")
                tooltipHtml = `<img style="max-width: 200px; height: auto" src="${pictureSource}" alt="Ligand preview" />`;
            } else {
                // Final fallback:
                // This will be hit if:
                // 1. pictureSource was "Not_available"
                // 2. pictureSource was undefined (neither data-picture nor rel was set)
                // 3. pictureSource was "Generate_from_SMILES" but SMILES was missing or OCL was undefined.
                if (smiles && pictureSource === "Generate_from_SMILES" && typeof OCL === 'undefined') {
                    console.warn("Tooltip: SMILES present and Generate_from_SMILES requested, but OCL library not found.");
                } else if (smiles && pictureSource === "Generate_from_SMILES" && !smiles) {
                     console.warn("Tooltip: Generate_from_SMILES requested, but no SMILES data found.");
                }
                tooltipHtml = `<img style="max-width: 200px; height: auto" src="${noImageUrl}" alt="Preview not available" />`;
            }
        }

        // Create the tooltip, but don't show it yet
        $("body").append(`<p class="tooltip-struct">${tooltipHtml}</p>`);
        const $tooltipInstance = $(".tooltip-struct");

        // Calculate position (this function now handles viewport collision)
        // Note: For accurate dimensions, the tooltip must be in the DOM.
        // If its content is dynamic and affects size significantly, you might need to
        // make it visible but off-screen (e.g., left: -9999px) to measure, then position.
        // However, with max-width:200px on images, current approach should be mostly fine.
        calculateAndSetPosition(e, $tooltipInstance);

        $tooltipInstance.fadeIn("fast");

    }).on('mouseleave', selector, function() {
        $(".tooltip-struct").remove();
    }).on('mousemove', selector, function(e) {
        // Update position on mouse move, also handling viewport collision
        calculateAndSetPosition(e, $(".tooltip-struct"));
    });
}
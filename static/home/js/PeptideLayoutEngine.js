// PeptideLayoutEngine.js
// Author: Binghan
// A robust engine for generating SVG representations of peptide sequences.

const PEPTIDE_LAYOUT_EPSILON = Number('1e-9');

/**
 * A robust engine for generating SVG representations of peptide sequences
 * in a smooth, snake-like layout.
 * @class PeptideLayoutEngine
 */
class PeptideLayoutEngine {
    /**
     * @param {object} options - Configuration options for the layout engine.
     * @param {number} [options.residueRadius=15] - The radius of each residue circle.
     * @param {number} [options.overlapAmount=8] - The pixel amount of overlap between adjacent residues.
     * @param {number} [options.residuesPerRow=10] - The preferred number of residues per row. Auto-adjusted for fixed widths.
     * @param {number} [options.rowSpacingFactor=1.2] - Multiplier for vertical space between rows (1.0 = touching, >1.0 = more space).
     * @param {number} [options.turnFactor=0.55] - Controls the "roundness" of turns. A value around 0.5 is recommended.
     * @param {object} [options.colorMap={}] - A map for custom residue colors, e.g., {'A': '#C8C8C8', 'K': '#0000FF'}.
     * @param {object} [options.tooltipProfile={}] - A map for custom tooltips, keyed by 1-based residue index, e.g., {1: 'N-terminus', 10: 'Active site'}.
     * @param {number|null} [options.fixedWidth=null] - If set, enforces a fixed width for the SVG canvas.
     * @param {number|null} [options.fixedHeight=null] - If set, enforces a fixed height for the SVG canvas.
     * @param {number|null} [options.maxLength=null] - An absolute maximum number of residues to display.
     * @param {number} [options.padding=20] - Padding in pixels around the peptide drawing.
     * @param {string} [options.fontFamily="Arial, sans-serif"] - Font family for the residue letters.
     * @param {number} [options.fontSize=12] - Font size for the residue letters.
     * @param {string} [options.defaultFillColor="white"] - Default fill color for residues not in the colorMap.
     * @param {string} [options.defaultStrokeColor="black"] - Default stroke color for the residue circles.
     * @param {string} [options.textColor="black"] - Color of the residue letters.
     */
    constructor({
        residueRadius = 15, overlapAmount = 8, residuesPerRow = 8,
        fontFamily = "Arial, sans-serif", fontSize = 12,
        defaultFillColor = "white",
        defaultStrokeColor = "black", textColor = "black",
        padding = 20,
        turnFactor = 0.55,
        rowSpacingFactor = 1.3,
        colorMap = {},
        tooltipProfile = {}, 
        fixedWidth = null, fixedHeight = null, maxLength = null
    } = {}) {
        this.residueRadius = residueRadius;
        this.overlapAmount = overlapAmount;
        this.residuesPerRow = Math.max(2, residuesPerRow);
        this.fontFamily = fontFamily;
        this.fontSize = fontSize;
        this.defaultFillColor = defaultFillColor;
        this.defaultStrokeColor = defaultStrokeColor;
        this.textColor = textColor;
        this.padding = padding;
        this.turnFactor = turnFactor;
        this.rowSpacingFactor = rowSpacingFactor;
        this.colorMap = colorMap;
        this.tooltipProfile = tooltipProfile;
        this.residueDiameter = this.residueRadius * 2;
        this.stepDistance = this.residueDiameter - this.overlapAmount;
        this.svgNS = "http://www.w3.org/2000/svg";
        this.fixedWidth = fixedWidth;
        this.fixedHeight = fixedHeight;
        this.maxLength = maxLength;
    }
    
    _sanitizeForId(inputString) {
        let s = String(inputString).replace(/[^\w.-]/g, '_');
        if (!s || !/^[a-zA-Z]/.test(s[0])) { s = "id_" + s; }
        return s;
    }

    _drawSingleResidue(dwgGroup, xCenter, yCenter, aaChar, seqNum, ligandIdPrefix) {
        const sanitizedPrefix = this._sanitizeForId(ligandIdPrefix);
        const resUniqueId = `${sanitizedPrefix}-res-${seqNum}`;
        const txtUniqueId = `${sanitizedPrefix}-txt-${seqNum}`;
        const tooltipText = this.tooltipProfile[seqNum] || `${aaChar}${seqNum}`;
        
        const circle = document.createElementNS(this.svgNS, "circle");
        circle.setAttribute("cx", xCenter);
        circle.setAttribute("cy", yCenter);
        circle.setAttribute("r", this.residueRadius);

        const fillColor = this.colorMap[aaChar.toUpperCase()] || this.defaultFillColor;
        circle.setAttribute("fill", fillColor);
        
        circle.setAttribute("stroke", this.defaultStrokeColor);
        circle.setAttribute("stroke-width", "1.5px");
        circle.setAttribute("class", "peptide-residue-circle");
        circle.setAttribute("id", resUniqueId);
        const circleTitle = document.createElementNS(this.svgNS, "title");
        circleTitle.textContent = tooltipText;
        circle.appendChild(circleTitle);
        dwgGroup.appendChild(circle);

        const text = document.createElementNS(this.svgNS, "text");
        text.setAttribute("x", xCenter);
        text.setAttribute("y", yCenter);
        text.setAttribute("text-anchor", "middle");
        text.setAttribute("dominant-baseline", "middle");
        text.setAttribute("font-family", this.fontFamily);
        text.setAttribute("font-size", `${this.fontSize}px`);
        text.setAttribute("fill", this.textColor);
        text.setAttribute("class", "peptide-residue-text");
        text.setAttribute("id", txtUniqueId);
        text.textContent = aaChar;
        const textTitle = document.createElementNS(this.svgNS, "title");
        textTitle.textContent = tooltipText;
        text.appendChild(textTitle);
        dwgGroup.appendChild(text);
    }
    
    _drawEllipsisResidue(dwgGroup, xCenter, yCenter, numTruncated) {
        const tooltipText = `... and ${numTruncated} more residue(s)`;
        const circle = document.createElementNS(this.svgNS, "circle");
        circle.setAttribute("cx", xCenter);
        circle.setAttribute("cy", yCenter);
        circle.setAttribute("r", this.residueRadius);
        circle.setAttribute("fill", this.defaultFillColor);
        circle.setAttribute("stroke", "grey");
        circle.setAttribute("stroke-dasharray", "2 2");
        dwgGroup.appendChild(circle);
        const circleTitle = document.createElementNS(this.svgNS, "title");
        circleTitle.textContent = tooltipText;
        circle.appendChild(circleTitle);
        const text = document.createElementNS(this.svgNS, "text");
        text.setAttribute("x", xCenter);
        text.setAttribute("y", yCenter);
        text.setAttribute("text-anchor", "middle");
        text.setAttribute("dominant-baseline", "middle");
        text.setAttribute("font-family", this.fontFamily);
        text.setAttribute("font-size", `${this.fontSize * 1.2}px`);
        text.setAttribute("fill", "grey");
        text.style.letterSpacing = "1px";
        text.textContent = "...";
        dwgGroup.appendChild(text);
        const textTitle = document.createElementNS(this.svgNS, "title");
        textTitle.textContent = tooltipText;
        text.appendChild(textTitle);
    }
    
    _calculateLayout(numToLayout, residuesPerRow) {
        if (numToLayout <= 0) return [];
        
        let d = "M 0 0";
        let x = 0, y = 0, direction = 1;
        const verticalStep = this.residueDiameter * this.rowSpacingFactor;
        const handleLength = verticalStep * this.turnFactor;
        
        let residuesLeft = numToLayout;
        while (residuesLeft > 0) {
            const residuesInThisRow = Math.min(residuesLeft, residuesPerRow);
            if (residuesInThisRow > 1) {
                const thisLineLength = (residuesInThisRow - 1) * this.stepDistance;
                d += ` l ${direction * thisLineLength} 0`;
                x += direction * thisLineLength;
            }
            residuesLeft -= residuesInThisRow;
            if (residuesLeft > 0) {
                const p1x = x + direction * handleLength;
                const p1y = y;
                const p2x = x + direction * handleLength;
                const p2y = y + verticalStep;
                const p3x = x;
                const p3y = y + verticalStep;
                d += ` C ${p1x} ${p1y}, ${p2x} ${p2y}, ${p3x} ${p3y}`;
                y += verticalStep;
                direction *= -1;
            }
        }

        const path = document.createElementNS(this.svgNS, 'path');
        path.setAttribute('d', d);
        const pathLength = path.getTotalLength();
        
        const centers = [];
        const numPointsToSample = pathLength > 0 ? Math.floor(pathLength / this.stepDistance + PEPTIDE_LAYOUT_EPSILON) + 1 : numToLayout;
        const finalNumPoints = Math.min(numToLayout, numPointsToSample);

        for (let i = 0; i < finalNumPoints; i++) {
            const point = path.getPointAtLength(i * this.stepDistance);
            centers.push({ x: point.x, y: point.y });
        }
        
        while (centers.length < numToLayout) {
             if (centers.length === 0) { centers.push({x:0, y:0}); continue; }
             const lastPoint = path.getPointAtLength(pathLength);
             if (centers.length > 0) {
                 const lastCenter = centers[centers.length - 1];
                 const dist = Math.sqrt(Math.pow(lastPoint.x - lastCenter.x, 2) + Math.pow(lastPoint.y - lastCenter.y, 2));
                 if (dist < this.stepDistance / 2) break;
             }
             centers.push({x: lastPoint.x, y: lastPoint.y});
             break; 
        }
        while (centers.length > numToLayout) { centers.pop(); }

        return centers;
    }

    generateSvgElement(peptideSequenceString, ligandId = "peptide") {
        if (!peptideSequenceString) { return document.createElementNS(this.svgNS, "svg");}
        
        let residuesPerRow = this.residuesPerRow;
        const fullLength = peptideSequenceString.length;
        let effectiveMaxLength = fullLength;

        if (this.maxLength !== null) {
            effectiveMaxLength = Math.min(fullLength, this.maxLength);
        }

        if (this.fixedWidth && this.fixedHeight) {
            const availableWidth = this.fixedWidth - (2 * this.padding);
            const turnWidthGuess = this.residueDiameter * this.turnFactor * 0.5;
            let calculatedResPerRow = Math.floor((availableWidth - turnWidthGuess) / this.stepDistance) + 1;
            calculatedResPerRow = Math.max(2, calculatedResPerRow);
            residuesPerRow = Math.min(this.residuesPerRow, calculatedResPerRow);
            
            const availableHeight = this.fixedHeight - (2 * this.padding);
            const rowHeight = this.residueDiameter * this.rowSpacingFactor;
            const maxRows = rowHeight > 0 ? Math.floor((availableHeight - this.residueDiameter) / rowHeight) + 1 : 1;
            
            const maxByFit = maxRows * residuesPerRow;
            effectiveMaxLength = Math.min(effectiveMaxLength, maxByFit);
        }
        
        const isTruncated = fullLength > effectiveMaxLength;
        const numToDraw = isTruncated ? effectiveMaxLength - 1 : effectiveMaxLength;
        const numTruncated = fullLength - numToDraw;

        if (numToDraw <= 0) {
             const svg = document.createElementNS(this.svgNS, "svg");
             if (this.fixedWidth && this.fixedHeight) { svg.setAttribute("width", this.fixedWidth); svg.setAttribute("height", this.fixedHeight); }
             return svg;
        }

        const layoutToDraw = this._calculateLayout(numToDraw, residuesPerRow);
        const layoutForEllipsis = this._calculateLayout(numToDraw + 1, residuesPerRow);
        const residueCenters = layoutToDraw;
        let ellipsisCenter = isTruncated && layoutForEllipsis.length > numToDraw ? layoutForEllipsis[layoutForEllipsis.length - 1] : null;

        const allPoints = [...residueCenters];
        if (ellipsisCenter) allPoints.push(ellipsisCenter);
        if (allPoints.length === 0) return document.createElementNS(this.svgNS, "svg");

        const allX = allPoints.map(p => p.x);
        const allY = allPoints.map(p => p.y);
        const minContentX = Math.min(...allX) - this.residueRadius;
        const maxContentX = Math.max(...allX) + this.residueRadius;
        const minContentY = Math.min(...allY) - this.residueRadius;
        const maxContentY = Math.max(...allY) + this.residueRadius;
        const contentWidth = maxContentX - minContentX;
        const contentHeight = maxContentY - minContentY;
        
        const svgRoot = document.createElementNS(this.svgNS, "svg");
        const mainGroup = document.createElementNS(this.svgNS, "g");
        mainGroup.setAttribute("class", "content-group");

        let transformOffsetX, transformOffsetY;
        if (this.fixedWidth && this.fixedHeight) {
            svgRoot.setAttribute("width", this.fixedWidth);
            svgRoot.setAttribute("height", this.fixedHeight);
            transformOffsetX = (this.fixedWidth - contentWidth) / 2 - minContentX;
            transformOffsetY = (this.fixedHeight - contentHeight) / 2 - minContentY;
        } else {
            svgRoot.setAttribute("width", contentWidth + 2 * this.padding);
            svgRoot.setAttribute("height", contentHeight + 2 * this.padding);
            transformOffsetX = -minContentX + this.padding;
            transformOffsetY = -minContentY + this.padding;
        }
        mainGroup.setAttribute("transform", `translate(${transformOffsetX}, ${transformOffsetY})`);
        svgRoot.appendChild(mainGroup);

        const sanitizedLigandId = this._sanitizeForId(ligandId);
        residueCenters.forEach((pos, i) => {
            this._drawSingleResidue(mainGroup, pos.x, pos.y, peptideSequenceString[i], i + 1, sanitizedLigandId);
        });
        if (isTruncated && ellipsisCenter) {
            this._drawEllipsisResidue(mainGroup, ellipsisCenter.x, ellipsisCenter.y, numTruncated);
        }
        return svgRoot;
    }
}
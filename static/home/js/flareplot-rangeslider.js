/*jslint browser: true */
/*jslint this */

/**
 * Create a time slider specifically for a flareplot.
 * Calls the `createD3RangeSlider` and returns the result, but also sets up bindings.
 * @param flareplot
 * @param containerSelector
 * @returns {{range: (function(number, number)), onChange: (function(Function))}}
 */
function createFlareplotSlider (flareplot, containerSelector) {
    var slider = createD3RangeSlider(0, flareplot.getNumFrames()-1, containerSelector, true);
    slider.onChange(function (newRange) {
        flareplot.framesSum(newRange.begin, newRange.end+1);
    });

    return slider;
}


/**
 * Create a d3 range slider that selects ranges between `rangeMin` and `rangeMax`, and add it to the
 * `containerSelector`. The contents of the container is laid out as follows
 * <code>
 * <div class="drag">
 *     <div class="handle WW"></div>
 *     <div class="handle EE"></div>
 * </div>
 * </code>
 * The appearance can be changed with CSS, but the `position` must be `relative`, and the width of `.slider` should be
 * left unaltered.
 *
 * @param rangeMin Minimum value of the range
 * @param rangeMax Maximum value of the range
 * @param containerSelector A CSS selection indicating exactly one element in the document
 * @returns {{range: function(number, number), onChange: function(function)}}
 */
function createD3RangeSlider (rangeMin, rangeMax, containerSelector, playButton) {
    "use strict";

    var minWidth = 10;

    var sliderRange = {begin: rangeMin, end: rangeMin};
    var changeListeners = [];
    var container = d3.select(containerSelector);
    var playing = false;
    var resumePlaying = false; // Used by drag-events to resume playing on release
    var playingRate = 30;
    var containerHeight = container.node().offsetHeight;

    // Set up play button if requested
    if (playButton) {
        // Wrap an additional container inside the main one, and set up a box-layout, see also
        // http://stackoverflow.com/questions/14319097/css-auto-resize-div-to-fit-container-width
        var box = container.append("div")
            .style("display", "box")
            .style("display", "-moz-box")
            .style("display", "-webkit-box")
            .style("box-orient", "horizontal")
            .style("-moz-box-orient", "horizontal")
            .style("-webkit-box-orient", "horizontal");

        var playBox = box.append("div")
            .style("width", containerHeight + "px")
            .style("height", containerHeight + "px")
            .style("margin-right", "10px")
            .style("box-flex", "0")
            .style("-moz-box-flex", "0")
            .style("-webkit-box-flex", "0")
            .classed("play-container", true);

        var sliderBox = box.append("div")
            .style("position", "relative")
            .style("min-width", (minWidth*2) + "px")
            .style("height", containerHeight + "px")
            .style("box-flex", "1")
            .style("-moz-box-flex", "1")
            .style("-webkit-box-flex", "1")
            .classed("slider-container", true);

        var playSVG = playBox.append("svg")
            .attr("width", containerHeight + "px")
            .attr("height", containerHeight + "px")
            .style("overflow", "visible");

        var circleSymbol = playSVG.append("circle")
            .attr("cx", containerHeight / 2)
            .attr("cy", containerHeight / 2)
            .attr("r", containerHeight / 2)
            .classed("button", true);

        var h = containerHeight;
        var stopSymbol = playSVG.append("rect")
            .attr("x", 0.3*h)
            .attr("y", 0.3*h)
            .attr("width", 0.4*h)
            .attr("height", 0.4*h)
            .style("visibility", "hidden")
            .classed("stop", true);

        var playSymbol = playSVG.append("polygon")
            .attr("points", (0.37*h) + "," + (0.2*h) + " " + (0.37*h) + "," + (0.8*h) + " " + (0.75*h) + "," + (0.5*h))
            .classed("play", true);

        //Circle that captures mouse interactions
        playSVG.append("circle")
            .attr("cx", containerHeight / 2)
            .attr("cy", containerHeight / 2)
            .attr("r", containerHeight / 2)
            .style("fill-opacity", "0.0")
            .style("cursor", "pointer")
            .on("click", togglePlayButton)
            .on("mouseenter", function(){
                circleSymbol
                    .transition()
                    .attr("r", 1.2 * containerHeight / 2)
                    .transition()
                    .attr("r", containerHeight / 2);
            });


    } else {
        sliderBox = container.append("div")
            .style("position", "relative")
            .style("height", containerHeight + "px")
            .style("min-width", (minWidth*2) + "px")
            .classed("slider-container", true);
    }

    //Create elements in container
    var slider = sliderBox
        .append("div")
        .attr("class", "slider");
    var handleW = slider.append("div").attr("class", "handle WW");
    var handleE = slider.append("div").attr("class", "handle EE");

    /** Update the `left` and `width` attributes of `slider` based on `sliderRange` */
    function updateUIFromRange () {
        var conW = sliderBox.node().clientWidth;
        var rangeW = sliderRange.end - sliderRange.begin;
        var slope = (conW - minWidth) / (rangeMax - rangeMin);
        var uirangeW = minWidth + rangeW * slope;
        var ratio = sliderRange.begin / (rangeMax - rangeMin - rangeW);
        if (isNaN(ratio)) {
            ratio = 0;
        }
        var uirangeL = ratio * (conW - uirangeW);

        slider
            .style("left", uirangeL + "px")
            .style("width", uirangeW + "px");
    }

    /** Update the `sliderRange` based on the `left` and `width` attributes of `slider` */
    function updateRangeFromUI () {
        var uirangeL = parseFloat(slider.style("left"));
        var uirangeW = parseFloat(slider.style("width"));
        var conW = sliderBox.node().clientWidth; //parseFloat(container.style("width"));
        var slope = (conW - minWidth) / (rangeMax - rangeMin);
        var rangeW = (uirangeW - minWidth) / slope;
        if (conW == uirangeW) {
            var uislope = 0;
        } else {
            var uislope = (rangeMax - rangeMin - rangeW) / (conW - uirangeW);
        }
        var rangeL = rangeMin + uislope * uirangeL;
        sliderRange.begin = Math.round(rangeL);
        sliderRange.end = Math.round(rangeL + rangeW);

        //Fire change listeners
        changeListeners.forEach(function (callback) {
            callback({begin: sliderRange.begin, end: sliderRange.end});
        });
    }

    // configure drag behavior for handles and slider
    var dragResizeE = d3.behavior.drag()
        .on("dragstart", function () {
            d3.event.sourceEvent.stopPropagation();
            resumePlaying = playing;
            playing = false;
        })
        .on("dragend", function () {
            if (resumePlaying) {
                startPlaying();
            }
        })
        .on("drag", function () {
            var dx = d3.event.dx;
            if (dx == 0) return;
            var conWidth = sliderBox.node().clientWidth; //parseFloat(container.style("width"));
            var newLeft = parseInt(slider.style("left"));
            var newWidth = parseFloat(slider.style("width")) + dx;
            newWidth = Math.max(newWidth, minWidth);
            newWidth = Math.min(newWidth, conWidth - newLeft);
            slider.style("width", newWidth + "px");
            updateRangeFromUI();
        });

    var dragResizeW = d3.behavior.drag()
        .on("dragstart", function () {
            this.startX = d3.mouse(this)[0];
            d3.event.sourceEvent.stopPropagation();
            resumePlaying = playing;
            playing = false;
        })
        .on("dragend", function () {
            if (resumePlaying) {
                startPlaying();
            }
        })
        .on("drag", function () {
            var dx = d3.mouse(this)[0] - this.startX;
            if (dx==0) return;
            var newLeft = parseFloat(slider.style("left")) + dx;
            var newWidth = parseFloat(slider.style("width")) - dx;

            if (newLeft < 0) {
                newWidth += newLeft;
                newLeft = 0;
            }
            if (newWidth < minWidth) {
                newLeft -= minWidth - newWidth;
                newWidth = minWidth;
            }

            slider.style("left", newLeft + "px");
            slider.style("width", newWidth + "px");

            updateRangeFromUI();
        });

    var dragMove = d3.behavior.drag()
        .on("dragstart", function () {
            d3.event.sourceEvent.stopPropagation();
            resumePlaying = playing;
            playing = false;
        })
        .on("dragend", function () {
            if (resumePlaying) {
                startPlaying();
            }
        })
        .on("drag", function () {
            var dx = d3.event.dx;
            var conWidth = sliderBox.node().clientWidth; //parseInt(container.style("width"));
            var newLeft = parseInt(slider.style("left")) + dx;
            var newWidth = parseInt(slider.style("width"));

            newLeft = Math.max(newLeft, 0);
            newLeft = Math.min(newLeft, conWidth - newWidth);
            slider.style("left", newLeft + "px");

            updateRangeFromUI();
        });

    handleE.call(dragResizeE);
    handleW.call(dragResizeW);
    slider.call(dragMove);

    //Click on bar
    sliderBox.on("mousedown", function (ev) {
        var x = d3.mouse(sliderBox.node())[0];
        var props = {};
        var sliderWidth = parseFloat(slider.style("width"));
        var conWidth = sliderBox.node().clientWidth; //parseFloat(container.style("width"));
        props.left = Math.min(conWidth - sliderWidth, Math.max(x - sliderWidth / 2, 0));
        props.left = Math.round(props.left);
        props.width = Math.round(props.width);
        slider.style("left", props.left + "px")
            .style("width", props.width + "px");
        updateRangeFromUI();
    });

    //Reposition slider on window resize
    window.addEventListener("resize", function () {
        updateUIFromRange();
    });

    function onChange(callback){
        changeListeners.push(callback);
        return this;
    }

    function setRange (b, e) {
        sliderRange.begin = b;
        sliderRange.end = e;

        updateUIFromRange();

        //Fire change listeners
        changeListeners.forEach(function (callback) {
            callback({begin: sliderRange.begin, end: sliderRange.end});
        });
    }


    /**
     * Returns or sets the range depending on arguments.
     * If `b` and `e` are both numbers then the range is set to span from `b` to `e`.
     * If `b` is a number and `e` is undefined the beginning of the slider is moved to `b`.
     * If both `b` and `e` are undefined the currently set range is returned as an object with `begin` and `end`
     * attributes.
     * If any arguments cause the range to be outside of the `rangeMin` and `rangeMax` specified on slider creation
     * then a warning is printed and the range correspondingly clamped.
     * @param b beginning of range
     * @param e end of range
     * @returns {{begin: number, end: number}}
     */
    function range(b, e) {
        var rLower;
        var rUpper;

        if (typeof b === "number" && typeof e === "number") {

            rLower = Math.min(b, e);
            rUpper = Math.max(b, e);

            //Check that lower and upper range are within their bounds
            if (rLower < rangeMin || rUpper > rangeMax) {
                console.log("Warning: trying to set range (" + rLower + "," + rUpper + ") which is outside of bounds (" + rangeMin + "," + rangeMax + "). ");
                rLower = Math.max(rLower, rangeMin);
                rUpper = Math.min(rUpper, rangeMax);
            }

            //Set the range
            setRange(rLower, rUpper);
        } else if (typeof b === "number") {

            rLower = b;
            var dif = sliderRange.end - sliderRange.begin;
            rUpper = rLower + dif;

            if (rLower < rangeMin) {
                console.log("Warning: trying to set range (" + rLower + "," + rUpper + ") which is outside of bounds (" + rangeMin + "," + rangeMax + "). ");
                rLower = rangeMin;
            }
            if(rUpper > rangeMax){
                console.log("Warning: trying to set range (" + rLower + "," + rUpper + ") which is outside of bounds (" + rangeMin + "," + rangeMax + "). ");
                rLower = rangeMax - dif;
                rUpper = rangeMax;
            }

            setRange(rLower, rUpper);
        }

        return {begin: sliderRange.begin, end: sliderRange.end};
    }

    function togglePlayButton () {
        if (playing) {
            stopPlaying();
        } else {
            startPlaying();
        }
    }

    function frameTick() {
        if (!playing) {
            return;
        }

        var limitWidth = rangeMax - rangeMin + 1;
        var rangeWidth = sliderRange.end - sliderRange.begin + 1;
        var delta = Math.min(Math.ceil(rangeWidth / 10), Math.ceil(limitWidth / 100));

        // Check if playback has reached the end
        if (sliderRange.end + delta > rangeMax) {
            delta = rangeMax - sliderRange.end;
            stopPlaying();
        }

        setRange(sliderRange.begin + delta, sliderRange.end + delta);

        setTimeout(frameTick, playingRate);
    }

    function startPlaying(rate) {
        if (rate !== undefined) {
            playingRate = rate;
        }

        if (playing) {
            return;
        }

        playing = true;
        if (playButton) {
            playSymbol.style("visibility", "hidden");
            stopSymbol.style("visibility", "visible");
        }
        frameTick();
    }

    function stopPlaying() {
        playing = false;
        if (playButton) {
            playSymbol.style("visibility", "visible");
            stopSymbol.style("visibility", "hidden");
        }
    }

    setRange(sliderRange.begin, sliderRange.end);

    return {
        range: range,
        startPlaying: startPlaying,
        stopPlaying: stopPlaying,
        onChange: onChange
    };
}


/**
 * Creates a tree and track selector for `flareplot` within `container`.
 *
 * If there are more than one tree in the graph a div will be added containing
 * radio buttons to switch between them. If the graph has more than one track
 * they'll be added in a div after. If there is only one tree and one track the
 * container will be left untouched. The resulting layout will look as follows:
 * <code>
 *   <div id="container" class="selector-container">
 *       <div class="treeselector">
 *           <label for="tree-radio-0">
 *               <input type="radio" name="tree-radios" id="tree-radio-0" checked />
 *               <span>Helical</span>
 *           </label>
 *           <label for="tree-radio-1">
 *               <input type="radio" name="tree-radios" id="tree-radio-1" />
 *               <span>Binding pocket</span>
 *           </label>
 *
 *       </div>
 *       <div class="trackselector">
 *           <label for="track-radio-0">
 *               <input type="radio" name="track-radios" id="track-radio-0" />
 *               <span>Helix colors</span>
 *           </label>
 *           <label for="track-radio-1">
 *               <input type="radio" name="track-radios" id="track-radio-1" />
 *               <span>Node centrality</span>
 *           </label>
 *       </div>
 *   </div>
 * </code>
 * No styling will be applied, but the `flareplot-selectors.css` file can be
 * used as is or as a template. Note that the container will have the class
 * "selector-container" added.
 *
 * @param flareplot A flareplot object created by calling `createFlareplot`.
 * @param container A div or DOM container assumed to be empty.
 * @returns {Function}
 */
function createSelectors(flareplot, container){

    return (function() {
        function initialize() {
            d3.select(container)
                .classed("selector-container",true);

            var treeNames = flareplot.getTreeNames();
            if(treeNames.length > 1){
                var div = d3.select(container)
                    .append("div")
                    .classed("treeselector", true);

                for (var i=0; i<treeNames.length; i++){
                    var label = div.append("label")
                        .attr("type","radio")
                        .attr("for","tree-radio-"+i)
                        .style("display","block");

                    label.append("input")
                        .attr("type","radio")
                        .attr("name","tree-radios")
                        .attr("id","tree-radio-"+i)
                        .on("click",function(){
                            var treeId = parseInt(this.id.split("-")[2]);
                            flareplot.setTree(treeId);
                        });

                    label.append("span")
                        .text(treeNames[i]);
                }
                div.select(":first-child input")
                    .attr("checked", "true");
            }

            var trackNames = flareplot.getTrackNames();
            if(trackNames.length > 1){
                var div = d3.select(container)
                    .append("div")
                    .classed("trackselector", true);

                for (var i=0; i<trackNames.length; i++){
                    var label = div.append("label")
                        .attr("type","radio")
                        .attr("for","track-radio-"+i);

                    label.append("input")
                        .attr("type","radio")
                        .attr("name","track-radios")
                        .attr("id","track-radio-"+i)
                        .on("click",function(){
                            var trackId = parseInt(this.id.split("-")[2]);
                            flareplot.setTrack(trackId);
                        });

                    label.append("span")
                        .text(trackNames[i]);
                }
                div.select(":first-child input")
                    .attr("checked", "true");
            }
        }

        initialize();

        return flareplot.getTrackNames().length > 1 || flareplot.getTreeNames().length > 1;
    }) ();
}

        function generateStyleDefs(svgDomElement) {
          // console.log(svgDomElement.className, svgDomElement.className=='hiveplot');
          var styleDefs = "";
          var sheets = document.styleSheets;
          for (var i = 0; i < sheets.length; i++) {
            // console.log(sheets[i]);
            var rules = sheets[i].cssRules;
            for (var j = 0; j < rules.length; j++) {
              var rule = rules[j];
              if (rule.style) {
                var selectorText = rule.selectorText;
                var elems = svgDomElement.querySelectorAll(selectorText);
                // console.log(selectorText,elems.length);

                if (elems.length) {
                  styleDefs += selectorText + " { " + rule.style.cssText + " }\n";
                }
              }
            }
          }

          var s = document.createElement('style');
          s.setAttribute('type', 'text/css');
          s.innerHTML = styleDefs;

          var defs = document.createElement('defs');
          defs.appendChild(s);
          svgDomElement.insertBefore(defs, svgDomElement.firstChild);
        }



        function downloadSVG(svgSelector, filetype) {
          var options = {};
          svg_class = svgSelector.closest(".panel-body").find("div").attr("class");
          svg_title = svgSelector.closest(".panel").find(".panel-heading").find("h3").text();
          var svgClone = svgSelector.clone();
          // remove data stuff that makes it slow
          console.log('download diagram',svg_class,svg_title);
          if (svg_class=='heatmap-container' && svg_title=='Interactions') {
            svgClone.find('rect').removeAttr('data-content');
            svgClone.find('rect').removeAttr('title');
            svgClone.find('rect').removeAttr('data-extra');

            // lots of logic here to figure out size of diagram which is made hard due to rotation and scaling.
            test = svgSelector.find('.svg-pan-zoom_viewport').find('g').get(0).getBBox();
            max_x = test.width;
            // calculate the width of the rotated diagram
            c = Math.sqrt(2*max_x**2);
            viewBox_x = parseInt(c);
            viewBox_y = parseInt(c*0.6);
            // svgClone.removeAttr('style');
            svgClone.removeAttr('preserveAspectRatio');
            svgClone.attr('height', '1200');
            svgClone.attr('width', '2000');
            svgClone.attr('viewBox', '0 0 '+(viewBox_x)+' '+(viewBox_y));

            // remove former transformations to ensure control
            svgClone.find('.svg-pan-zoom_viewport').removeAttr('transform');
            svgClone.find('.svg-pan-zoom_viewport').find('g').removeAttr('transform');
            svgClone.find('.svg-pan-zoom_viewport').attr('transform','translate(10,'+parseInt(max_x*0.7)+') scale(-1,1) rotate(135)');

          } else if (svg_class=='heatmap-container' ) {
            svgClone.find('rect').removeAttr('data-content');
            svgClone.find('rect').removeAttr('title');
            svgClone.find('rect').removeAttr('data-extra');

            // Figure out how much to translate svg to offset for the headers
            textNode1 = svgSelector.find('text').first().get(0);
            bbox1 = textNode1.getBBox();
            extra_x = parseInt(bbox1.height)-2;
            max_x = parseInt(svgClone.find('rect').last().attr('x'))+extra_x;
            svgClone.removeAttr('style');
            svgClone.removeAttr('preserveAspectRatio');
            svgClone.attr('height', '2000');
            svgClone.attr('width', '2000');
            svgClone.attr('viewBox', '0 0 '+(max_x+2)+' '+(max_x+2));
            svgClone.find('.svg-pan-zoom_viewport').removeAttr('transform');
            svgClone.find('.svg-pan-zoom_viewport').attr('transform','translate('+extra_x+','+extra_x+')');
          } else if (svg_class=='flareplot-container') {

            svgClone.attr("version", "1.1");
            svgClone.attr("xmlns", "http://www.w3.org/2000/svg");
            svgClone.attr("xmlns:xlink", "http://www.w3.org/1999/xlink");
            svgClone.attr('height', '3000');
            svgClone.attr('width', '3000');
            svgClone.attr('viewBox', '0 0 800 800');
          } else if (svg_class=='ngl-container') {
            // Requires canvas-toBlob.js
            svgSelector.get(0).toBlob(function(blob) {
                saveAs(blob, "ngl_view_distances.png");
            });
            return
          } else if (svg_class=='hiveplot-container') {
            svgClone.attr("version", "1.1");
            svgClone.attr("xmlns", "http://www.w3.org/2000/svg");
            svgClone.attr("xmlns:xlink", "http://www.w3.org/1999/xlink");
            svgClone.attr('height', '3000');
            svgClone.attr('width', '3000');
            svgClone.attr('viewBox', '0 0 900 1000');
            svgClone.find('.svg-pan-zoom_viewport').removeAttr('transform');

          }

          // svgClone.css('background-color','white');

          // Get styles
          generateStyleDefs(svgClone.get(0));
          var escapedSVG = new XMLSerializer().serializeToString(svgClone.get(0));

          // Escape special characters
          escapedSVG = escapedSVG.replace(/°/g,"&#176;");
          escapedSVG = escapedSVG.replace(/Å/g,"&#197;");

          downloadURI('data:image/svg+xml;base64,' + window.btoa(escapedSVG), 'distances_'+svg_class+'.svg');

          // If SVG, end here.
          if (filetype=='svg') {
            downloadURI('data:image/svg+xml;base64,' + window.btoa(escapedSVG), 'distances_'+svg_class+'.svg');
            return
          }

          var imgsrc = 'data:image/svg+xml;base64,'+ btoa( unescape( encodeURIComponent( escapedSVG ) ) ); // Convert SVG string to data URL
          var image = new Image();
          image.onload = function() {
            var canvas = document.createElement("canvas");
            canvas.width = image.width;
            canvas.height = image.height;
            var context = canvas.getContext('2d');
            context.drawImage(image, 0, 0);
            // // Requires canvas-toBlob.js
            canvas.toBlob( function(blob) {
              saveAs( blob, 'distances_'+svg_class+'.png' ); // FileSaver.js function
            });
          };
          image.src = imgsrc;

        }

      var filename_prefix = {"single-crystal-tab" : "single", "single-crystal-group-tab" : "group", "two-crystal-groups-tab" : "comparison"}
      function downloadCurrentTableCSV(){
          // grab all data from current table
          var table = $("#" + currentTab + " .contact-browser.active .dataTable").DataTable();

          // TODO change this function to download nice Excel/CSV including headers
          // There is an issue in our implementation with the Datatables buttons option

          // Convert current dataTable to CSV
          csv = Papa.unparse(table.rows().data().toArray());

          // remove HTML tags
          csv = csv.replace(/(<([^>]+)>)/ig, "")

          // Download file
          downloadURI('data:text/csv;charset=UTF-8,' + encodeURI(csv), "Structure_analyzer-" + filename_prefix[currentTab]+".csv");
      }

      function downloadCurrentTableExcel(){
          GlobalTableToExcel(currentTab + " .contact-browser.active .dataTable", 'GPCRdb structure analyzer', "Structure_analyzer-" + filename_prefix[currentTab]+".xls")
      }


      $(document).ready(function() {
        $('.btn-download.png').click(function() {
            DownloadElement = $(this).closest(".panel-heading").next().find("svg,canvas");
            downloadSVG(DownloadElement,'png');
        });
      });

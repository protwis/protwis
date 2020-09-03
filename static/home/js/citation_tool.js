function citation_tool(url) {
	// Modal part
	$('#modal_table tbody').empty();
    var modal = document.getElementById('citation-tool');
    var span = document.getElementsByClassName("close")[0];
    modal.classList.add("modal");
    modal.style.display = "block";
    span.onclick = function() {
        modal.style.display = "none";
        $(".article_option").remove();
        $(".article").remove();
    }
    window.onclick = function(event) {
        if (event.target == modal) {
            modal.style.display = "none";
            $(".article_option").remove();
            $(".article").remove();
        }
    }

    // Process location
    var this_site = parse_url(url);
    var highlight_main = false;
    if (window.location.pathname==='/') {
    	highlight_main = true
    }

    // AJAX to create citation list
    var dropdown_articles = document.getElementById("dropdown_articles");
    var article_list = document.getElementById("article_list");
    var env = url.split('//')[1].split('/')[0];
    var main_ref_id = false;
    var cit_request = new XMLHttpRequest();
    cit_request.open('GET', url.split('/')[0] + '/citations');
    cit_request.onload = function() {
		var data = JSON.parse(cit_request.responseText)
		// Create HTML once per call
		if (dropdown_articles.length===0) {
			var articles = {};
			for (i = 0; i < data.length; i++) {
				var site = parse_url(data[i][0]);

				var option = document.createElement("option");
				option.text = data[i][4];
				option.setAttribute("id", site.replace('#', '_').replace('.', '_'));
				option.classList.add("article_option");
				if (site===this_site) {
					option.setAttribute("selected", "selected");
				}
				dropdown_articles.add(option);

				// Reformat data
				if (!(data[i][5] in articles)) {
					articles[data[i][5]] = {'tools':{}};
					articles[data[i][5]]['tools'][site] = option.text;
					articles[data[i][5]]['main'] = data[i][3];
					articles[data[i][5]]['authors'] = data[i][6];
					articles[data[i][5]]['year'] = data[i][7];
					articles[data[i][5]]['reference'] = data[i][8];
					articles[data[i][5]]['journal'] = data[i][9];
					articles[data[i][5]]['doi'] = data[i][10];
					if (data[i][3]) {
						main_ref_id = site;
					}
				}
				else {
					articles[data[i][5]]['tools'][site] = option.text;
					if (data[i][3]) {
						main_ref_id = site;
					}
				}
			}
			// Create HTML
			for (var key in articles) {
				var entry = document.createElement("ul");
				entry.setAttribute("id", articles[key]['doi']);
				var entry_holder = document.createElement("div");
				var tools_list = document.createElement("div");
				// Loop through tools list to add them with ID
				var selected = false;
				for (var tool_key in articles[key]['tools']) {
					if (tool_key===this_site) {
						selected = true;
					}
					var tool = document.createElement("div");
					tool.innerHTML = articles[key]['tools'][tool_key].bold();
					tool.setAttribute("id", tool_key.replace('#', '_').replace('.', '_')+'_tag');
					tool.style.display = 'inline-block';
					tools_list.appendChild(tool);
					var separator = document.createElement("div");
					separator.innerHTML = '\u00A0&bull;\u00A0';
					separator.style.display = 'inline-block';
					tools_list.appendChild(separator);
				}
				tools_list.removeChild(tools_list.lastChild);
				tools_list.style.display = 'block';

				var ref = document.createElement("div");
				var p = document.createElement("p");
				p.innerHTML = articles[key]['authors'] + ", " + '"' + key + '", '
				// Link
				var a = document.createElement("a");
				a.innerHTML = articles[key]['journal'].italics() + ", " + articles[key]['year'] + ", " + articles[key]['reference'];
      			a.title = key;
      			a.href = "https://doi.org/"+articles[key]['doi'];
      			a.target = "_blank";
      			p.appendChild(a);
				ref.appendChild(p);
				ref.style.display = 'block';

				entry_holder.appendChild(tools_list);
				entry_holder.appendChild(ref);
				entry.appendChild(entry_holder);
				entry.classList.add('article');
				if (selected) {
					entry.classList.add('highlight_reference');
				}
				article_list.appendChild(entry);
			}
			if (highlight_main || $('.highlight_reference').length==0) {
				highlight_article(main_ref_id);
				$('#dropdown_articles option#'+main_ref_id).attr('selected', 'selected');
			}
		}
    };
    cit_request.send();
}

function highlight_article(key) {
	$('.highlight_reference').removeClass('highlight_reference');
	$('#'+key+'_tag').parent().parent().parent().addClass('highlight_reference');
}

function parse_url(url) {
	var url_split = url.split('/');
	if (url_split[url_split.length-1]==="") {
    	var this_site = url_split[url_split.length-2];
    	if (this_site==='targetselection') {
    		this_site = url_split[url_split.length-3]+'-'+this_site;
    	}
    }
    else {
    	var this_site = url_split[url_split.length-1];
    	if (this_site==='targetselection') {
    		this_site = url_split[url_split.length-2]+'-'+this_site;
    	}
    }
    return this_site
}

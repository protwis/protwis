function citation_tool(url, cite_id) {
	// Modal part
	$('#citation_modal_table tbody').empty();
    var modal = document.getElementById('citation-tool');
    var span = document.getElementById("close-citation");
    modal.classList.add("modal");
    setTimeout(function(){
    	modal.style.display = "block";
    }, 200)

    span.onclick = function() {
        modal.style.display = "none";
        $(".article_option").remove();
        $(".article").remove();
        $(".cite-submenu").remove();
    }
    window.onclick = function(event) {
        if (event.target == modal) {
            modal.style.display = "none";
            $(".article_option").remove();
            $(".article").remove();
            $(".cite-submenu").remove();
        }
    }

    // Process location
    var this_site = parse_url_long(url);
    var highlight_main = false;
    if (window.location.pathname==='/') {
    	highlight_main = true
    }

    // AJAX to create citation list
    var dropdown_articles = document.getElementById("dropdown_articles");
    var article_list = document.getElementById("article_list");
    var env = url.split('//')[1].split('/')[0];
    var main_ref_id = [];
    var cit_request = new XMLHttpRequest();
    var domains = ["gpcrdb.org", "gproteindb.org", "arrestindb.org", "biasedsignalingatlas.org"];
    var filter_for = false;
		let cite_id_missing = typeof cite_id === "undefined";
    if ((cite_id_missing && env==="gpcrdb.org") || cite_id==="cite_gpcrdb") {
    	filter_for = "gpcrdb";
    }
    if ((cite_id_missing && env==="gproteindb.org") || cite_id==="cite_gproteindb") {
    	filter_for = "gproteindb";
    }
    if ((cite_id_missing && env==="arrestindb.org") || cite_id==="cite_arrestindb") {
    	filter_for = "arrestindb";
    }
    if ((cite_id_missing && env==="biasedsignalingatlas.org") || cite_id==="cite_biasedsignalingatlas") {
    	filter_for = "biasedsignalingatlas";
    }
    cit_request.open('GET', url.split('/')[0] + '/citations');
    cit_request.onload = function() {
		var data = JSON.parse(cit_request.responseText)
		for (let i = 0; i < data.length; i++) {
			var link = "";
			if (!domains.includes(env)) {
				link = "/" + data[i][0].split("//")[1].split("/").slice(1).join("/");
			}
			else {
				link = data[i][0];
			}
			if (link[link.length-1]==='/') {
				link = link.substring(0, link.length-1);
			}
			if (link==='/GPCR-Gprotein_Mutations.xlsx') {
				link = 'https://files.gpcrdb.org/GPCR-Gprotein_Mutations.xlsx'
			}
			else if (link==='/mutations.html#mutation-data-submission') {
				link = 'https://docs.gpcrdb.org/mutations.html#mutation-data-submission'
			}
			if ($('a[href="'+link+'"]').first().parent().parent().parent().first().hasClass('dropdown-submenu')) {
				data[i].push($('a[href="'+link+'"]').first().parent().parent().parent().parent().parent().first().find(">:first-child").text().trim());
			}
			else {
				if (link==="/protein/ste2_yeast") {
					data[i].push("GPCRdb");
				}
				else {
					data[i].push($('a[href="'+link+'"]').first().parent().parent().parent().first().find(">:first-child").text().trim());
				}
			}
		}
		// Create HTML once per call
		var articles = {};
		var tags = [];
		for (let i = 0; i < data.length; i++) {
			if (data[i][11]==='') {
				continue;
			}
			if (filter_for==="gpcrdb" && data[i][11]!=="GPCRdb") {
				continue;
			}
			else if (filter_for==="gproteindb" && data[i][11]!=="GproteinDb" && data[i][7]===2022) {
				continue;
			}
			else if (filter_for==="arrestindb" && data[i][11]!=="ArrestinDb") {
				continue;
			}
			else if (filter_for==="biasedsignalingatlas" && data[i][11]!=="Biased Signaling Atlas") {
				continue;
			}
			else if (filter_for==="arrestindb" && data[i][11]==="ArrestinDb") {
				data[i][5] = "The arrestin database, ArrestinDb";
				data[i][6] = "Jimmy Caroli, Gáspár Pándy-Szekeres, Alexander S. Hauser, György M. Keserű, Albert J. Kooistra and David E. Gloriam";
				data[i][7] = 2022;
				data[i][8] = "";
				data[i][9] = "Manuscript";
				data[i][10] = "";
			}
			else if (filter_for==="biasedsignalingatlas" && data[i][11]==="Biased Signaling Atlas") {
				data[i][5] = "A community Biased Signaling Atlas";
				data[i][6] = "Jimmy Caroli, Alibek Mamyrbekov, Kasper Harpsøe, Sahar Gardizi, Linda Dörries, Eshan Ghosh, Alexander S. Hauser, Albert J. Kooistra, and David E. Gloriam";
				data[i][7] = 2023;
				data[i][8] = "19, 531–535";
				data[i][9] = "Nature Chemical Biology";
				data[i][10] = "https://doi.org/10.1038/s41589-023-01292-8";
			}
			var site = parse_url_long(data[i][0]);
			tags.push(site);
			// Dropdown menus
			var a_sub, submenu, submenu_ul;
			if (document.getElementById(data[i][11])===null) {
				submenu = document.createElement("li");
				submenu.setAttribute("id", data[i][11]);
				submenu.classList.add("dropdown-submenu");
				submenu.classList.add("cite-submenu");
				a_sub = document.createElement("a");
				a_sub.innerHTML = data[i][11];
				a_sub.setAttribute('tabindex', '0');
				submenu.appendChild(a_sub);
				submenu_ul = document.createElement("ul");
				submenu_ul.classList.add('dropdown-menu');
				submenu_ul.classList.add('dropdown-auto-overflow')
				submenu.appendChild(submenu_ul);
			}
			else {
				submenu = document.getElementById(data[i][11])
				submenu_ul = submenu.getElementsByTagName('ul')[0];
			}

			// Dropdown options
			var option = document.createElement("li");
			var option_link = document.createElement("a");
			option_link.innerHTML = data[i][4];
			option_link.classList.add("article_option");
			option.appendChild(option_link);
			option.setAttribute("id", site.replace('#', '_').replace('.', '_'));
			option.setAttribute("role", "button");

			if (site===this_site) {
				option.setAttribute("selected", "selected");
				$('#page_select_button').html(data[i][4]);
			}

			submenu_ul.appendChild(option);
			dropdown_articles.appendChild(submenu);

			// Reformat data
			if (!(data[i][5] in articles)) {
				articles[data[i][5]] = {'tools':{}};
				articles[data[i][5]]['tools'][site] = option.innerHTML;
				articles[data[i][5]]['main'] = data[i][3];
				articles[data[i][5]]['authors'] = data[i][6];
				articles[data[i][5]]['year'] = data[i][7];
				articles[data[i][5]]['reference'] = data[i][8];
				articles[data[i][5]]['journal'] = data[i][9];
				articles[data[i][5]]['doi'] = data[i][10];
				articles[data[i][5]]['menu'] = data[i][11];
				if (data[i][3] === "GPCRdb" || data[i][3] === "GproteinDb") {
					main_ref_id.push(site);
				}
			}
			else {
				articles[data[i][5]]['tools'][site] = option.innerHTML;
				if (data[i][3]=="GPCRdb" || data[i][3]=="GproteinDb") {
					main_ref_id.push(site);
				}
			}
		}

		// Create HTML
		for (var key in articles) {
			var main_ref = false
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
				tool.style.display = 'none';
				tools_list.appendChild(tool);
				// var separator = document.createElement("div");
				// separator.innerHTML = '\u00A0&bull;\u00A0';
				// separator.style.display = 'inline-block';
				// tools_list.appendChild(separator);
				if (main_ref_id.includes(tool_key)) {
					main_ref = true;
				}
			}
			// tools_list.removeChild(tools_list.lastChild);
			tools_list.style.display = 'block';

			var ref = document.createElement("div");
			var p = document.createElement("p");
			var d1 = document.createElement("div");
			var d2 = document.createElement("div");
			var d3 = document.createElement("div");
			d1.innerHTML = key;
			d2.innerHTML = articles[key]['authors'];

			// Link
			var a = document.createElement("a");
			a.innerHTML = articles[key]['journal'].italics() + ", " + articles[key]['year'] + ", " + articles[key]['reference'];
  			a.title = key;
  			a.href = "https://doi.org/"+articles[key]['doi'];
  			a.target = "_blank";
  			d3.appendChild(a);
			p.appendChild(d1);
			p.appendChild(d2);
			p.appendChild(d3);
			ref.appendChild(p);
			ref.style.display = 'block';

			entry_holder.appendChild(tools_list);
			entry_holder.appendChild(ref);
			entry.appendChild(entry_holder);
			entry.classList.add('article');
			if (selected) {
				entry.classList.add('highlight_reference');
			}
			if (main_ref) {
				var main_ref_list = document.getElementById("main_reference");

				main_ref_list.appendChild(entry);
			}
			else {
				article_list.appendChild(entry);
			}
		}
		if (env=="gproteindb.org") {
			main_ref_id = main_ref_id[0];
		}
		else {
			main_ref_id = main_ref_id[1];
		}
		/*var color = "#BE00BE";
		if (filter_for==="gproteindb") {
			color = "#F46615";
		}
		else if (filter_for==="arrestindb") {
			color = "#55AE35";
		}
		$(".article").first().css("border-top","4px solid "+color);*/
		// $(".article").eq(1).css("border-top","4px solid #BE00BE");
		if (highlight_main || $('.highlight_reference').length==0) {
			highlight_article(main_ref_id);
			$('#dropdown_articles option#'+main_ref_id).attr('selected', 'selected');
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

function parse_url_long(url) {
	var url_split = url.split('/');
	if (url_split[url_split.length-1]==="") {
		if (url_split.length>5) {
			return url_split[3]+"-"+url_split[4];
		}
		else {
			return url_split[3];
		}
	}
	else if (url_split.length===4) {
		return url_split[3];
	}
	else {
		return url_split[3]+"-"+url_split[4];
	}
}

function toggle_widget() {
	$("#ref_widget").animate({width: "toggle"});
    if ($("#widget_chevron").hasClass("glyphicon glyphicon-chevron-left")) {
        $("#widget_chevron").removeClass("glyphicon glyphicon-chevron-left");
        $("#widget_chevron").addClass("glyphicon glyphicon-chevron-right");
        $("#widget_infosign").show();
    }
    else {
        $("#widget_chevron").removeClass("glyphicon glyphicon-chevron-right");
        $("#widget_chevron").addClass("glyphicon glyphicon-chevron-left");
        $("#widget_infosign").hide();
        setTimeout(function() {
        	$("#ref_widget").animate({width: "toggle"});
        	$("#widget_chevron").removeClass("glyphicon glyphicon-chevron-left");
        	$("#widget_chevron").addClass("glyphicon glyphicon-chevron-right");
        	$("#widget_infosign").show();
        }, 8000)
    }
}

async function check_for_video(url) {
	var this_site = parse_url_long(url);
	var cit_request = new XMLHttpRequest();
    cit_request.open('GET', url.split('/')[0] + '/citations');
    cit_request.onload = function() {
		var data = JSON.parse(cit_request.responseText)
		var video = false;
		for (i = 0; i < data.length; i++) {
			var site = parse_url_long(data[i][0]);
			if (site===this_site && data[i][1]!=null) {
				video = data[i][1];
				$('#icon_video').attr('href', video);
				$('#icon_video').css('display','inline');
			}
		}
	}
	cit_request.send();
}

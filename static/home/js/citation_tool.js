function citation_tool(url) {
	// Modal part
	$('#modal_table tbody').empty();
    var modal = document.getElementById('citation-tool');
    var span = document.getElementsByClassName("close-citation")[0];
    modal.classList.add("modal");
    setTimeout(function(){
    	modal.style.display = "block";
    }, 200)
    
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
		for (i = 0; i < data.length; i++) {
			var link = '/'+data[i][0].split('//')[1].split('/').slice(1).join('/');
			if (link[link.length-1]==='/') {
				link = link.substring(0, link.length-1);
			}
			if (link==='/GPCR-Gprotein_Mutations.xlsx') {
				link = 'https://files.gpcrdb.org/GPCR-Gprotein_Mutations.xlsx'
			}
			else if (link==='/mutations.html#mutation-data-submission') {
				link = 'https://docs.gpcrdb.org/mutations.html#mutation-data-submission'
			}
			if ($('a[href="'+link+'"]').parent().parent().parent().first().hasClass('dropdown-submenu')) {
				data[i].push($('a[href="'+link+'"]').parent().parent().parent().parent().parent().first().find(">:first-child").text().trim());
			}
			else {
				data[i].push($('a[href="'+link+'"]').parent().parent().parent().first().find(">:first-child").text().trim());
			}
		}

		// Create HTML once per call
		var articles = {};
		var tags = [];
		for (i = 0; i < data.length; i++) {
			if (data[i][11]==='') {
				continue;
			}
			var site = parse_url(data[i][0]);
			if (tags.includes(site)) {
				site = site+i.toString();
				tags.push(site);
			}
			else {
				tags.push(site);
			}

			// Dropdown menus
			if (document.getElementById(data[i][11])===null) {
				var submenu = document.createElement("li");
				submenu.setAttribute("id", data[i][11]);
				submenu.classList.add("dropdown-submenu");
				var a_sub = document.createElement("a");
				a_sub.innerHTML = data[i][11];
				a_sub.setAttribute('tabindex', '0');
				submenu.appendChild(a_sub);
				var submenu_ul = document.createElement("ul");
				submenu_ul.classList.add('dropdown-menu');
				submenu_ul.classList.add('dropdown-auto-overflow')
				submenu.appendChild(submenu_ul);
			}
			else {
				var submenu = document.getElementById(data[i][11])
				var submenu_ul = submenu.getElementsByTagName('ul')[0];
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
				if (data[i][3]) {
					main_ref_id = site;
				}
			}
			else {
				articles[data[i][5]]['tools'][site] = option.innerHTML;
				if (data[i][3]) {
					main_ref_id = site;
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
				if (tool_key===main_ref_id) {
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

function toggle_widget() {
	$('#ref_widget').animate({width: 'toggle'});
    if ($('#ref_widget_openclose').hasClass('glyphicon glyphicon-chevron-left')) {
        $('#ref_widget_openclose').removeClass('glyphicon glyphicon-chevron-left')
        $('#ref_widget_openclose').addClass('glyphicon glyphicon-chevron-right')
    }
    else {
        $('#ref_widget_openclose').removeClass('glyphicon glyphicon-chevron-right')
        $('#ref_widget_openclose').addClass('glyphicon glyphicon-chevron-left')
        setTimeout(function() {
        	$('#ref_widget').animate({width: 'toggle'});
        	$('#ref_widget_openclose').removeClass('glyphicon glyphicon-chevron-left')
        	$('#ref_widget_openclose').addClass('glyphicon glyphicon-chevron-right')
        }, 8000)
    }
}

function check_for_video(url) {
	var this_site = parse_url(url);
	var cit_request = new XMLHttpRequest();
    cit_request.open('GET', url.split('/')[0] + '/citations');
    cit_request.onload = function() {
		var data = JSON.parse(cit_request.responseText)
		var video = false;
		for (i = 0; i < data.length; i++) {
			var site = parse_url(data[i][0]);
			if (site===this_site && data[i][1]!=null) {
				video = data[i][1];
				$('#icon_video').attr('href', video);
				$('#icon_video').css('display','inline');
			}
		}
	}
	cit_request.send();
}

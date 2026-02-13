function get_publication_html(publications) {
    var publication_html = '<ul>';
    var publications_found = false;
    $.each(publications, function(index, value) {
        publication_html += '<li>';
        if (value['title']) {
            publication_html += value['title'];
            if (value['title'].slice(-1) != '.') {
                publication_html += '.'
            }
        }
        if (value['authors']) {
            publication_html += '<br/>' + value['authors'];
            if (value['authors'].slice(-1) != '.') {
                publication_html += '.'
            }
        }
        if (value['doi'] && value['year'] && !value['journal_name']) {
            publication_html +='<br/><u><a href="https://doi.org/'+value['doi']+'" target="_blank">'+
                'https://doi.org/'+value['doi']+'</a></u>'+', '+value['year'];
        } else if (value['doi'] && !value['year'] && !value['journal_name']) {
            publication_html +='<br/><u><a href="https://doi.org/'+value['doi']+'" target="_blank">'+
                'https://doi.org/'+value['doi']+'</a></u>';
        } else if (value['doi'] && value['year'] && value['journal_name']) {
            publication_html +='<br/><u><a href="https://doi.org/'+value['doi']+'" target="_blank">'+'<i>'+value['journal_name']+'</i>, '+value['year']+'</a></u>';
        } else if (!value['doi']) {
            let journal_name_available = false;
            if (value['journal_name']) {
                publication_html += '<br><i>'+value['journal_name']+'</i>';
                journal_name_available = true;
            } 
            if (value['year']) {
                if (journal_name_available) {
                    publication_html += ', '
                } else {
                    publication_html += '<br>'
                }
                publication_html += value['year'];
            }
            publication_html += '.';
        }
        publication_html += '</li>';
        publications_found = true;
    });
    publication_html += '</ul>';
    return [publication_html,publications_found]
}

function numberToWords(num) {
    const ones = ["", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine"];
    const teens = ["ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen", 
                    "sixteen", "seventeen", "eighteen", "nineteen"];
    const tens = ["", "", "twenty", "thirty", "forty", "fifty", 
                    "sixty", "seventy", "eighty", "ninety"];

    if (num === 0) return "zero";

    if (num < 10) return ones[num];
    if (num < 20) return teens[num - 10];
    if (num < 100) 
        return tens[Math.floor(num / 10)] + (num % 10 ? "-" + ones[num % 10] : "");
    if (num < 1000) 
        return ones[Math.floor(num / 100)] + " hundred" + (num % 100 ? " " + numberToWords(num % 100) : "");
    if (num < 10000) 
        return ones[Math.floor(num / 1000)] + " thousand" + (num % 1000 ? " " + numberToWords(num % 1000) : "");

    return "number too big";
}



function getCitationByUrl(url, success, settings=null) {
    var new_settings = $.extend(true,{},settings,{
        url: '/citation_by_url/object/',
        data: { url: removeQuery(url) },
        success: function() {
            var args = Array.prototype.slice.call(arguments); // ES5 Convert array-like
                                                    // object arguments to a real array

            var data = args[0]; // data[0] = citations for the specific tool page
            var publication_html = '';
            var default_publication_html = '';
            var main = data[1]; // DB name alias
            var only_default_citation = data[2]; // Are citations the same for
                                                // all the tools of the DB?

            var no_citation = data[3]; // Must citations be hidden for the url?

            var default_citation_data = data[4]; // Common citations for
                                                // all the tools in the DB.
            var publications_found = false;
            var default_publications_found = false;
            var msg = '';
            var gpcrdb_is_main = false;
            var html = '';

            if (!no_citation) {

                if (data[0].length > 0) {
                    let r_pub = get_publication_html(data[0][0]['publication']);
                    publication_html = r_pub[0];
                    publications_found = r_pub[1];
                }
                if (default_citation_data.length > 0) {
                    let r_pub = get_publication_html(default_citation_data[0]['publication']);
                    default_publication_html = r_pub[0];
                    default_publications_found = r_pub[1];
                }
                var lower_main = main.toLowerCase();

                if (default_publications_found || publications_found) {

                    if (main) {
                        
                        if (lower_main.indexOf("gpcr") !== -1) {
                            msg = 'Please support GPCRdb development by citing the both the latest GPCRdb paper as well as the paper(s) describing the given page.';
                            gpcrdb_is_main = true;
                        } else if (lower_main.indexOf("arrestin") !== -1) {
                            msg = 'When using ArrestinDb for publications, please support the development by citing the paper.';
                            if (default_citation_data[0]['publication'].length > 1) {
                                msg = 'When using ArrestinDb for publications, please support the development by citing its '+numberToWords(default_citation_data[0]['publication'].length)+' papers, as each paper describes different data and tools.'
                            }
                        } else if (lower_main.indexOf("gprotein") !== -1) {
                            msg =  'When using GproteinDb for publications, please support the development by citing its '+numberToWords(default_citation_data[0]['publication'].length)+' papers, as each paper describes different data and tools.'
                        } else if (lower_main.indexOf("bias") !== -1) {
                            msg = 'When using BSA for publications, please support the development by citing the paper.';
                        }
                    }

                    if (gpcrdb_is_main) {

                        if (!publications_found) {
                            msg = 'Please support GPCRdb development by citing the latest GPCRdb paper as well as the paper(s) describing the given page once they are available.';
                            html = msg+'<br><br>' +
                            '<b>Latest GPCRdb paper:</b><br>'+
                            default_publication_html
                        } else {
                            html = msg+'<br><br>' +
                            '<b>Latest GPCRdb paper:</b><br>'+
                            default_publication_html +
                            '<b>Current page paper(s):</b><br>'+
                            publication_html;
                        }

                    } else {
                        html = msg+'<br><br>' +
                            default_publication_html;
                    }
                }
            }
            args.unshift(no_citation); // arguments in reverse order
                                        // as they are added at the beginning of args list
            args.unshift(html);
            success.apply(this, args);
        }, 
    });

    return $.get(new_settings);
}
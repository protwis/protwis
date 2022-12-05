function joinus(type) {
	var title = "Join us - ";
	var message = "Testing";
	switch (type){
		case "annotation":
			title += "Data annotation";
			message = "The GPCRdb team is looking for students and researchers that can help with annotation of literature or structure data for joint coming publications. Please contact our lead developer Gáspár Pándy-Szekeres  (<a href='mailto:gaspar@sund.ku.dk'>gaspar@sund.ku.dk</a>) if you would like to know more.";
			break;
		case "development":
			title += "Database development";
			message = "The GPCRdb team is constantly looking for talented programmers either for remote collaboration, to fill an open position or to apply together for a scientific project.  Please contact our lead developer lead developer Gáspár Pándy-Szekeres  (<a href='mailto:gaspar@sund.ku.dk'>gaspar@sund.ku.dk</a>) if you would like to know more.";
			break;
		case "collaboration":
			title += "Scientific collaboration";
			message = "If you have an experimental dataset that you would like to make available in GPCRdb or would like assistance with computational analyses, please send a short description of the suggested collaboration to the head of GPCRdb David Gloriam (<a href='mailto:david.gloriam@sund.ku.dk'>david.gloriam@sund.ku.dk</a>).";
			break;
		case "media":
			title += "Social media";
			message = "To follow the latest updates and news from GPCRdb join our <a href='https://www.linkedin.com/groups/8104679/' target='_blank'>LinkedIn group</a> and follow us on Twitter (<a href='https://twitter.com/gpcrdb' target='_blank'>@gpcrdb</a>). The latest releases are also listed on the main page of GPCRdb.";
			break;
		case "feedback":
		default:
			title += "Give feedback";
			message = "We appreciate suggestions of new and feedback on existing resources to <a href='mailto:info@gpcrdb.org'>info@gpcrdb.org</a> and bug reports to GitHub repository (<a href='https://github.com/protwis/protwis/issues' target='_blank'>https://github.com/protwis/protwis/issues</a>). While we do not have the resources to implement all, the suggestions are prioritized in regular meetings of the development team.";
			break;
	}
	// create join us modal and fill with text + show
	$("#join-us .modal-title").html(title);
	$("#join-us .modal-body p").html(message);

	// trigger modal
	$('#join-us').modal('show');
	// Move as a temporary workaround
	$("#join-us").appendTo("body");
}

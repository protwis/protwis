function superposition(oTable, columns, site, hide_first_column) {
    // oTable: DataTable object of table of entries
    // columns: Column indeces of oTable to be extracted to build table for reference selection. First column has to be structure/model string used for superposition workflow
    // site: Structure browser or Homology model browser (add new logic when expanding to new sites)
    if(window.location.hash === "#keepselection") {}
    else {
        ClearSelection('targets');
        ClearSelection('reference');
    }

    var checked_data = oTable.rows('.alt_selected').data();
    if (checked_data.length===0) {
        window.alert('No entries selected for superposition')
        return 0;
    }
    var selected_ids = []
    if (site==='structure_browser') {
        for (i = 0; i < checked_data.length; i++) {
            var div = document.createElement("div");
            div.innerHTML = checked_data[i][6];
            if (typeof div.innerText !== "undefined") {
                selected_ids.push(div.innerText.replace(/\s+/g, ''));
            } else {
                selected_ids.push(div.textContent.replace(/\s+/g, ''));
            }
        }
        AddToSelection('targets', 'structure_many', selected_ids.join(","));
    } else if (site==='homology_model_browser') {
        for (i = 0; i < checked_data.length; i++) {
            var div = document.createElement("div");
            div.innerHTML = checked_data[i][12];
            var state = checked_data[i][3];
            if (checked_data[i][4]==='Yes') {
                div.innerHTML = checked_data[i][11];
                if (typeof div.innerText !== "undefined") {
                    selected_ids.push(div.innerText.replace(/\s+/g, '')+"_refined");
                } else {
                    selected_ids.push(div.textContent.replace(/\s+/g, '')+"_refined");
                }
            }
            else {
                if (typeof div.innerText !== "undefined") {
                    selected_ids.push(div.innerText.replace(/\s+/g, '')+"_"+state);
                    console.log(div.textContent.replace(/\s+/g, '')+"_"+state);
                } else {
                    selected_ids.push(div.textContent.replace(/\s+/g, '')+"_"+state);
                    console.log(div.textContent.replace(/\s+/g, '')+"_"+state);
                }
            }
        }
        AddToSelection('targets', 'structure_models_many', selected_ids.join(","));
    } // add new logic here for new site

    $('#modal_table tbody').empty();
    var modal = document.getElementById('myModal');
    var span = document.getElementsByClassName("close")[0];
    modal.classList.add("modal");
    // console.log('hasclass',$('#myModal').hasClass('modal'));
    // console.log(modal);
    modal.style.display = "block";
    span.onclick = function() {
        modal.style.display = "none";
    }
    window.onclick = function(event) {
        if (event.target == modal) {
            modal.style.display = "none";
        }
    }
    var checked_data = oTable.rows('.alt_selected').data();
    // var needed_columns = [6,1,2,3,4,5,10,26]
    for (i=0; i<checked_data.length; i++) {
        var div = document.createElement("div");
        row = document.createElement('tr');
        var checkbox = document.createElement('input');
        checkbox.type = "checkbox";
        row.appendChild(checkbox);
        var column_count = 0;
        columns.forEach(function(column) {
            div.innerHTML = checked_data[i][column];
            cell = document.createElement('td');
            textnode = document.createTextNode(div.innerText.replace(/\s+/g, ''));
            if(column_count===0 && typeof hide_first_column==="undefined") {}
            else if (column_count===0 && hide_first_column===true) {
                cell.style.display = "none";
            }
            cell.appendChild(textnode);
            row.appendChild(cell);
            column_count++;
        });
        $('#modal_table tbody').append(row)
    }
    $('#modal_table tbody tr').click(function() {
        if (site==='structure_browser') {
            AddToSelection('reference', 'structure', $(this).children().eq(1).text());
        } else if (site==='homology_model_browser') {
            console.log($(this).children().eq(9));
            console.log($(this).children().eq(8));
            if ($(this).children().eq(9).text()==='Yes') {
                AddToSelection('reference', 'structure',  $(this).children().eq(11).text()+"_refined");
                console.log($(this).children().eq(11).text()+"_refined");
            }
            else {
                var state = $(this).children().eq(8).text();
                AddToSelection('reference', 'structure_model', $(this).children().eq(1).text()+"_"+state);
                console.log($(this).children().eq(1).text()+"_"+state);
            }
        } // Add logic here for new site
        $(this).children(':first').prop("checked",true);
        window.location.href = '/structure/superposition_workflow_index';
    })
};
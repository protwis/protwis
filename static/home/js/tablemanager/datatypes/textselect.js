// This class represents the "textselect" data type used in the table manager.
// It provides a select input filtering interface, allowing users to filter table data based on options selected from a dropdown list.
// The select2 library is used to enhance the select input with features such as multi-select, search, and tagging.
// The options for the select input can be populated either from the client-side data or fetched from the server-side endpoint, depending on the configuration of the DataTable.

import { DataTypeBase } from "./datatypebase.js";

class TextSelect extends DataTypeBase {

    select2options = {
                            multiple: true,
                            closeOnSelect: true, 
                            dropdownAutoWidth: true,
                            placeholder: "Select...",
                            tags: true,
                            allowClear: true 
                        }    
    
    createFilterInterface(tableManagerReference) {  
      if (!this.allow_filter) {
        return null;
      }

      //-----------------//
      //  HTML Elements  //
      //-----------------//

      let filterInput = this._createInputForSelectFilter();
      
      //-----------------------------------------------------------------//
      /// Callbacks for post processing and event handler registration ////
      //-----------------------------------------------------------------//

      // Add function to populate select options as a callback to be run once data is loaded into the table,
      // as we need access to the dataset to determine the unique values for the options.
      if(tableManagerReference.dataTableReference.settings()[0].oFeatures.bServerSide === true) {
        tableManagerReference.dataLoadedCallBacks.push( () => this._populateTextSelectFilterOptionsServerSide(this.json_id, this.col_idx, filterInput, tableManagerReference.dataTableReference));
      } else {
        tableManagerReference.dataLoadedCallBacks.push( () => this._populateTextSelectFilterOptionsClientSide(this.json_id, this.col_idx, filterInput, tableManagerReference.dataTableReference));
      }

      // Add select2 initialization as a callback to be run when document is ready to ensure the relevant DOM element is available  
      tableManagerReference.documentReadyCallBacks.push( 
        () => $(`#${filterInput.attributes.id.value}`).select2(this.select2options)
      );

      //Add event handler registration to callback for when document is ready
      tableManagerReference.documentReadyCallBacks.push( () => $(`#${filterInput.attributes.id.value}`).on('change', (e) => this.onSelectChangeHandler(e, this.col_idx, tableManagerReference)));

      return filterInput;
  }

    //------------------//
    //  Event Handlers  //
    //------------------//
    onSelectChangeHandler = (e, col_idx, tableManagerReference) => {
        let currentSearchValue = tableManagerReference.dataTableReference.columns(col_idx).search()[0] || "%";

        let newSearchValues = $('#' + e.target.id).select2('data').map( (o) => o.text)

        const allOptions = $('#' + e.target.id).find('option').toArray()

        let useRegex = newSearchValues.length > 1; //DataTables regex search is required for multiple values, as we will join them with the OR operator "|"
        let tagState;
        if (this.select2options.tags) {
          tagState = allOptions.map((o) => o.dataset['select2Tag'] === 'true' || false)
          if (tagState.includes(true)) {
            useRegex = true;
          }
        }

        if (useRegex) {
          newSearchValues = newSearchValues.map((v) => v.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')); // escape regex special characters
          for (let i = 0; i < newSearchValues.length; i++) {
            const optIdx = allOptions.map( (o) => o.text).indexOf(newSearchValues[i]);
            const isTag = tagState[optIdx];
            if (isTag) {
              newSearchValues[i] = "(.*" + newSearchValues[i] + ".*)"; // wrap tag values with wildcards for "contains" search
            } else {
              newSearchValues[i] = "(^" + newSearchValues[i] + "$)"; // wrap non-tag values with anchors for "exact match" search
            }
          }
        }

        let newSearchValue = newSearchValues.join("|"); // join multiple values with OR operator for regex search

        if (newSearchValue !== currentSearchValue) {
          tableManagerReference.dataTableReference.columns(col_idx).search(newSearchValue, useRegex).draw();
        }
    }

    _createInputForSelectFilter() {
        let filterInput = document.createElement("select");
        filterInput.setAttribute("data-field", this.json_id);
        filterInput.setAttribute("id", this.tableId + "_" + this.json_id + "_search");
        filterInput.setAttribute("name", this.tableId + "_" + this.json_id + "_select");

        filterInput.classList.add("col-filter");
        filterInput.classList.add("select-filter");

        if (this.cssClassFilterInput) {
          filterInput.classList.add(this.cssClassFilterInput);
        }

        return filterInput;
    }

    _populateTextSelectFilterOptionsFromData(data, ColumnJsonId, filterInput, dataTableReference) {
      
      if (!data || data.length === 0) {
        console.warn("No data available to populate select filter options for column " + ColumnJsonId);
        return;
      }

      const uniqueValues = new Set([...data]);
      
      uniqueValues.forEach(value => {
        let option = document.createElement("option");
        option.textContent = value;
        // set option value to empty string for placeholder option to allow clearing the filter when selected
        option.value = value; 
        filterInput.appendChild(option);
      });
    }

    _populateTextSelectFilterOptionsClientSide(ColumnJsonId, columnIdx, filterInput, dataTableReference) {
      const data = dataTableReference.column(columnIdx).data().toArray()
      if (!data || data.length === 0) {
        console.warn("No data available to populate select filter options for column " + ColumnJsonId);
        return;
      }

      const uniqueValues = new Set([...data]);
      
      uniqueValues.forEach(value => {
        let option = document.createElement("option");
        option.textContent = value;
        // set option value to empty string for placeholder option to allow clearing the filter when selected
        option.value = value; 
        filterInput.appendChild(option);
      });
    }

    _populateTextSelectFilterOptionsServerSide(ColumnJsonId, columnIdx, filterInput, dataTableReference) {
      const dataEndpoint = dataTableReference.ajax.url()
      const optionsEndpoint = dataEndpoint + "/options/" + ColumnJsonId + "/"

      fetch(optionsEndpoint)
        .then(response => {
          if (!response.ok) {
            throw new Error(`HTTP error when fetching options for column ${ColumnJsonId}! Status: ${response.status}`);
          }
          return response.json();
        })
        .then(data => {
          this._populateTextSelectFilterOptionsFromData(data, ColumnJsonId, filterInput, dataTableReference);
        });

    }
}

export { TextSelect };
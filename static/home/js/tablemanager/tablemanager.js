import { Column } from "./column.js";
import { DataTypeFactory } from "./datatypefactory.js";

// This class represents the TableManager, which is responsible for managing the data table and its associated columns, filters, and rendering logic.
// The TableManager interacts with the back end TableProvider to fetch configuration and content data, then initializes a DataTable with the appropriate settings and column definitions.
// TableManager supports both server-side and client-side processing modes, allowing for flexible handling of large datasets and dynamic filtering capabilities.

class TableManager {

  /**
   * Constructs a new TableManager instance.
   * @param {Array} columnSpecification - An array of column specification objects defining the properties and behaviors of each column in the data table.
   * @param {string} tableId - The ID of the HTML table element to be managed by the TableManager.
   * @param {DataTypeFactory} DataTypeFactoryReference - An optional reference to a DataTypeFactory instance for initializing data types. Only required if you want to override the default factory with a custom one.
  **/
  constructor(
    columnSpecification,
    tableId, 
    DataTypeFactoryReference
  ) {
    this.dataTypeFactory = DataTypeFactoryReference || new DataTypeFactory();
    this.columns = this.initialiseColumns(columnSpecification);
    this.tableId = tableId;
    //Array to store functions that need to be run after data is loaded into the table, 
    // e.g. to populate filter interfaces that require access to the dataset to determine their options (e.g. select dropdowns populated with unique values from the dataset)
    this.dataLoadedCallBacks = []; 
    this.documentReadyCallBacks = []; 
  }

  /**
   * Initialize the columns based on the provided column specifications. 
   * @param {Array} columnSpecification - An array of column specification objects defining the properties and behaviors of each column in the data table.
   * @returns {Array} An array of initialized Column instances.
  **/   
  initialiseColumns(columnSpecification) {
    return columnSpecification.map(
      (colSpec) =>
        new Column(colSpec, this.dataTypeFactory),
    );
  }

  /**
   * This method initializes the DataTable with the provided options, setting up the table structure, headers, and filters based on the column specifications.
   * It also handles the execution of callbacks for data loading and document readiness, ensuring that filter interfaces are properly populated and initialized.
   * @param {Object} dataTableOptions -  Object containing the configuration options for initializing the DataTable, including column definitions, server-side processing settings, and any custom callbacks.
   */
  initialiseDataTable(dataTableOptions) {   
    
    try {
      if (this._columnIndexesAbsent()) {
        // If no columns have col_idx, we assume they are in the correct order and assign col_idx based on their position in the array
        this._createColumnIndexes();
      } else {
        this.validateColumnSpecifications();
      }
    } catch (error) {
      alert(
        "Error in column specifications during DataTable initialisation:" + error.message        
      );
      return;
    }


    //Returns [table, tableHead,  primaryHeaderRow, tableBody] HTML object references after adding them to DOM
    let tableComponents = this.createPrimaryTableStructure();

    this.populateColumnTitleHeaders(tableComponents.primaryHeaderRow);

    if (!dataTableOptions.columns) {
      dataTableOptions.columns = this.columns.map((col) =>
        this.columnToDataTableColDef(col),
      );
    } else {
      console.log(
        "DataTable options already contains column definitions, skipping automatic column definition generation from specifications.",
      );
    }

    // We need to wait until the initial data has been loaded into the table before populating the filter headers, 
    // as the filter interfaces (e.g. select dropdowns) may need to be populated with unique values from the dataset which are not available until the data is loaded. 
    const userInitComplete = dataTableOptions.initComplete; // store any initComplete function that was supplied in config, may be undefined
    dataTableOptions.initComplete = (settings, json) => {
      this.dataLoadedCallBacks.forEach((callBack) => callBack());

      // invoke original callback, if one was provided
      if (typeof userInitComplete === "function") {
        userInitComplete(settings, json);
      }
    }    

    this.dataTableReference = $(`#${this.tableId}`).DataTable(dataTableOptions);

    tableComponents = this.createSecondaryTableStructure(tableComponents);

    this.populateSuperHeaders(tableComponents.superHeaderRow);

    this.populateFilterHeaders(tableComponents.filterHeaderRow);

    //We need to wait until the document is fully loaded before initializing any filter interface plugins (e.g. select2) to ensure the relevant DOM elements are available
    // so we store the initialization of such plugins as callbacks to be run on document ready, and push them into the documentReadyCallBacks array
    $(document).ready(() => {
      this.documentReadyCallBacks.forEach((callBack) => callBack());
    });

  }

/**
 * Populates the primary table structure by creating and appending the table header, primary header row, and table body elements to the specified HTML table element.
 * @returns - A dictionary-like object containing references to the table, table header, primary header row, and table body HTML elements
 */
  createPrimaryTableStructure() {
    let table = document.getElementById(this.tableId);

    if(!table) {
      throw new Error(`No table element found with id "${this.tableId}". Please ensure an HTML <table> element with this id exists in the DOM before initializing the TableManager.`);
    }

    let tableHead = this._createTableHeaderSection();

    let primaryHeaderRow = document.createElement("tr");
    this.primaryHeaderRowId = this.tableId + "_primaryHeaders";
    primaryHeaderRow.id = this.primaryHeaderRowId;
    tableHead.append(primaryHeaderRow);

    let tableBody = document.createElement("tbody");
    tableBody.id = this.tableId + "_body";
    table.prepend(tableBody);

    table.prepend(tableHead);
    table.append(tableBody);

    return {
      table: table,
      tableHead: tableHead,
      primaryHeaderRow: primaryHeaderRow,
      tableBody: tableBody,
    };
  }

  /**
   * Creates the secondary table structure by adding super-headers and filter headers to the table.
   * @param {Object} TableComponents - A dictionary-like object containing references to the table, table header, primary header row, and table body HTML elements.
   * @returns {Object} - The updated table components with secondary structure added.
   */
  createSecondaryTableStructure(TableComponents) {
    let superHeaderRow = null;
    if (this._hasSuperHeaders()) {
      this.superHeaderRowId = this.tableId + "_superHeaders";
      superHeaderRow = this._createSuperHeaderRow(this.superHeaderRowId);
      TableComponents.tableHead.prepend(superHeaderRow);
    }

    let filterHeaderRow = null;
    if (this._hasFilters()) {
      this.filterHeaderRowId = this.tableId + "_filterHeaders";
      filterHeaderRow = this._createFilterHeaderRow(this.filterHeaderRowId);
      TableComponents.tableHead.append(filterHeaderRow);
    }

    return {
      ...TableComponents,
      superHeaderRow: superHeaderRow,
      filterHeaderRow: filterHeaderRow,
    };
  }

  //-------------------------/
  //                         /
  //  Table header creation  /
  //                         /
  //-------------------------/

  // Creates the <thead> section of the table and assigns it an ID based on the table's ID.
  _createTableHeaderSection() {
    const tableHead = document.createElement("thead");
    this.headerSetId = this.tableId + "_headerSet";
    tableHead.id = this.headerSetId;
    return tableHead;
  }

  //---------------------------/
  //  Super(/grouped) headers  /
  //---------------------------/

  // Checks if all columns have a column_group property defined, indicating that super-headers should be created for the table.
  _hasSuperHeaders() {
    return this.columns.every((col) => col.column_group);
  }

  // Checks if any columns have allow_filter set to true, indicating that a filter row should be created for the table.
  _hasFilters() {
    return this.columns.some((col) => col.allow_filter);
  }
  
  /**
   * Uses column_group property of column specifications to create and append a super-header row to the table,
   * grouping columns under shared headers as specified. 
   * @Example 
   * If columns 0-4 have column_group 'Classification' and columns 5-7 have column_group 'State',
   * a super-header row is created with 'Classification' spanning columns 0-4 and 'State' spanning columns 5-7.
   * @param {HTMLElement} superHeaderRow - The <tr> element representing the super-header row to which the grouped headers will be appended.
  **/
  populateSuperHeaders(superHeaderRow) {
    //Add super-header row if all columns have a column groups defined in specifications, otherwise remove super-header row (which would just be empty)
    if (!this._hasSuperHeaders()) {
      console.log(
        "Column groups incomplete or absent in specifications, skipping super-header creation.",
      );
      return;
    }

    const superHeaderStructure = this._getSuperHeaderStructure(
      this.columns,
    );

    superHeaderStructure.forEach(({ name, colspan }) => {
      let th = document.createElement("th");
      th.textContent = name;
      th.colSpan = colspan;
      th.classList.add("super_header_cell");
      superHeaderRow.append(th);
    });
  }

  // Creates a super-header row element with the specified ID and returns it.
  _createSuperHeaderRow(superHeaderRowId) {
    const superHeaderRow = document.createElement("tr");
    superHeaderRow.id = superHeaderRowId;
    return superHeaderRow;
  }

  /**
   * Computes the required colspan values for super-headers based on the column_group properties of the column specifications
   * @param {Array} columnsSpecifications 
   * @returns - An array of { name, colspan } objects in col_idx order, e.g. [{ name: 'Classification', colspan: 5 }, { name: 'State', colspan: 3 }, ...] 
   */
  _getSuperHeaderStructure(columnsSpecifications) {
    try {
      this.validateColumnSpecifications();
    } catch (error) {
      console.error(
        "Error in column specification when getting super header structure:",
        error.message,
      );
      return [];
    }

    const groups = [];

    for (const col of columnsSpecifications) {
      const last = groups[groups.length - 1];
      if (last && last.name === col.column_group) {
        last.colspan++;
      } else {
        groups.push({ name: col.column_group, colspan: 1 });
      }
    }

    return groups;
  }

  //-------------------/
  //  Primary headers  /
  //-------------------/

  /**
   * Populates the primary header row with column title cells.
   * @param {HTMLElement} primaryHeaderRow - The <tr> element representing the primary header row.
   */
  populateColumnTitleHeaders(primaryHeaderRow) {
    this.columns.forEach((col) => {
      // add header
      const th = document.createElement("th");
      th.textContent = col.title;
      
      if (col.cssClassHeaderCell) {
        th.classList.add(col.cssClassHeaderCell);
      } else {
        th.classList.add("primary_header_cell");
      }
      
      primaryHeaderRow.append(th);
    });
  }

  //--------------/
  //  Filter row  /
  //--------------/

  // Filter Header row

  // Creates a filter header row element with the specified ID and returns it.
  _createFilterHeaderRow(filterHeaderRowId) {
    const filterHeaderRow = document.createElement("tr");
    filterHeaderRow.id = filterHeaderRowId;
    return filterHeaderRow;
  }

  /**
   * Populates the filter header row with filter interface elements for each column that allows filtering.
   * @param {HTMLElement} filterHeaderRow - The <tr> element representing the filter header row.
   */
  populateFilterHeaders(filterHeaderRow) {
    this.columns.forEach((col) => {
      const filterTh = document.createElement("th");

      if (col.cssClassFilterCell) {
        filterTh.classList.add(col.cssClassFilterCell);
      } else {
        filterTh.classList.add("filter_header_cell");
      }

      if (col.allow_filter) {
        const filterInterface = col.createFilterInterface(this);
        if (filterInterface) {
          //Add filter <input> elements  to filter header cell
          if (!Array.isArray(filterInterface)) {
            filterTh.appendChild(filterInterface);
          } else {
            filterInterface.forEach((el) => filterTh.appendChild(el));
          }
        }
      }
      filterHeaderRow.append(filterTh);
    });
  }

  //------------------------------------/
  //  Column specification validation   /
  //------------------------------------/

  validateColumnSpecifications() {
    if (!this._checkColumnsIdxContinuous(this.columns)) {
      throw new Error(
        "Column specifications must have continuous col_idx values starting from 0 with no gaps",
      );
    }

    if (!this._checkColumnsIdxSorted(this.columns)) {
      throw new Error(
        "Column specifications must be sorted in ascending order of col_idx",
      );
    }

    if (!this._checkConsecutiveColumnGroups(this.columns)) {
      throw new Error(
        "Column specifications must have column groups that form contiguous blocks of columns",
      );
    }
  }

  _createColumnIndexes() {
    this.columns.forEach((col, index) => col.col_idx = index);
  }

  _columnIndexesAbsent() {
    return this.columns.every((col) => !col.hasOwnProperty("col_idx"));
  }

  // Checks if the col_idx values of the column specifications are continuous, meaning that they form a sequence without any gaps.
  _checkColumnsIdxContinuous(columnsSpecifications) {
    let someColumnsMissingColIdx = columnsSpecifications.some((col) => !col.hasOwnProperty("col_idx"));

    if (someColumnsMissingColIdx) {      
        throw new Error("Either all or none of the column specifications should have a col_idx property");
    }    

    for (let i = 1; i < columnsSpecifications.length; i++) {
      if (
        columnsSpecifications[i].col_idx - columnsSpecifications[i - 1].col_idx !== 1
      ) {
        return false;
      }
    }
    return true;
  }

  // Checks if the col_idx values of the column specifications are sorted in ascending order.
  _checkColumnsIdxSorted(columnsSpecifications) {
    for (let i = 1; i < columnsSpecifications.length; i++) {
      if (
        columnsSpecifications[i].col_idx < columnsSpecifications[i - 1].col_idx
      ) {
        return false;
      }
    }
    return true;
  }

  /**
   * Checks if the column groups are consecutive, meaning that they form contiguous blocks of columns.
   * @param {Array} columnSpecification - An array of column specification objects defining the properties and behaviors of each column in the data table.
   * @returns Boolean indicating whether the column groups are consecutive/contiguous
  */
  _checkConsecutiveColumnGroups(columnsSpecifications) {
    const sorted = [...columnsSpecifications].sort(
      (a, b) => a.col_idx - b.col_idx,
    );
    const closedGroups = new Set();
    let currentGroup = null;

    for (const col of sorted) {
      if (col.column_group !== currentGroup) {
        // true = end of group
        if (closedGroups.has(col.column_group)) {
          //group previously closed
          throw new Error(
            `column_group "${col.column_group}" is non-contiguous: it appears at col_idx ${col.col_idx} after being interrupted by other groups`,
          );
        }
        if (currentGroup !== null) closedGroups.add(currentGroup);
        currentGroup = col.column_group;
      }
    }
    return true;
  }

  //------------------------------/
  //  Data Table Initialisation   /
  //------------------------------/

  /**
    * Converts a column specification into a DataTables column definition object.
    * @param {Column} column - The Column instance representing the column specification.
    * @returns {Object} - The DataTables column definition object with the appropriate data and render properties.
  */
  columnToDataTableColDef(column) {
    const colDef = { data: column.json_id };
    colDef.render = (data, type, row, meta) => column.dataTableRenderer(data, type, row, meta);
    return colDef;
  }


  //-------------------------/
  //  Configuration request  /
  //-------------------------/

  /**
   * Fetches the configuration for a data table and its columns.
   * @param {string} dataTableConfigName - The name of the data table configuration.
   * @param {string} dataTableConfigVariant - The variant of the data table configuration.
   * @param {string} columnConfigName - The name of the column configuration.
   * @param {string} columnConfigVariant - The variant of the column configuration.
   * @returns {Promise<Object>} - A promise resolving to the fetched configuration objects.
   */
  static async fetchConfiguration(dataTableConfigName, dataTableConfigVariant, columnConfigName, columnConfigVariant,) {
    
     const [cols, datatable] = await Promise.all([
      TableManager._fetchConfigurationAsync(
        "column",
        columnConfigName,
        columnConfigVariant),      
      TableManager._fetchConfigurationAsync(
        "datatable",
        dataTableConfigName,
        dataTableConfigVariant,
      )]);      

      return {
        columnsSpecifications: cols,
        datatableConfiguration: datatable,
      };

  }

  // Fetch table configuration JSON from API endpoint
  // 
  /**
   * Fetches a configuration asynchronously from the API endpoint. Endpoint format: /table_provider/configuration/{configurationType}/{tableName}/{configurationVariant}
   * @param {string} configurationType - The type of the configuration.
   * @param {string} tableName - The name of the table.
   * @param {string} configurationVariant - The variant of the configuration.
   * @returns {Promise<Object>} - A promise resolving to the fetched configuration object.
   */
  static async _fetchConfigurationAsync(
    configurationType,
    tableName,
    configurationVariant,
  ) {
    if (!configurationType || !tableName || !configurationVariant) {
      throw new Error(
        "configurationType, tableName and configurationVariant are required",
      );
    }

    const url = `/table_provider/configuration/${encodeURIComponent(
      configurationType,
    )}/${encodeURIComponent(tableName)}/${encodeURIComponent(
      configurationVariant,
    )}`;

    try {
      const resp = await fetch(url, {
        method: "GET",
        headers: { Accept: "application/json" },
        credentials: "same-origin",
      });

      if (!resp.ok) {
        const text = await resp.text();
        throw new Error(
          `Failed to fetch configuration: ${resp.status} ${resp.statusText} - ${text}`,
        );
      }

      const data = await resp.json();
      return data;
    } catch (err) {
      console.error("Error fetching table configuration:", err);
      throw err;
    }
  }
}

export { TableManager };
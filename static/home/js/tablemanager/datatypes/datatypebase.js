// This is the base class for all data types used in the table manager.
// It provides a common interface and default implementations for methods 
// that can be overridden by subclasses to provide specific functionality for each data type.

class DataTypeBase {
  // class methods
  constructor() { }
  
  createFilterInterface(tableManagerReference) { 
    // This method should be overridden by subclasses to create the appropriate filter interface for the data type.
    // Default implementation returns null, indicating no filter interface is provided, resulting in no filter controls being added/built.
    return false;
  }
  
  dataTableRenderer(data, type, row, meta) { 
    // This method should be overridden by subclasses to provide the appropriate rendering logic for the data type (if required).
    // Default implementation returns the data as-is.
    return data;
  } 
  
}

export { DataTypeBase };
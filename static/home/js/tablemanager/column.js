// This column class represents a data column definition for the table manager. 
// It is responsible for managing the properties and behaviors of a specific column in the data table, 
// including its data type, filtering capabilities, and rendering logic.



class Column{
    
    // The constructor takes a column specification object (colSpec) and a data type factory (dataTypeFactory) as parameters.
    // The datafactory is used to initialize the appropriate data type class for the column based on its data_type property.
    // The datatype returned from the factory is then set as the prototype of the Column instance, 
    // allowing it to inherit the methods and properties of the specific data type class.
    constructor(colSpec, dataTypeFactory) {
        //Assign all properties from colSpec to the Column instance
        Object.assign(this, { ... colSpec } )
        // Set the prototype of the Column instance to the appropriate data type class based on the data_type property
        Object.setPrototypeOf(this, dataTypeFactory.initialiseDataType(this));
    }
    
}

export { Column }
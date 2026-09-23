import { Text } from "./datatypes/text.js";
import { Numeric } from "./datatypes/numeric.js";
import { TextSelect } from "./datatypes/textselect.js";
import { WebLinkText } from "./datatypes/weblinktext.js";

// This class is responsible for managing the different data types used in the table manager. 
// It facilitate dynamic initialization of data types by calling the constructor for the 
// data type class associated with the data_type string supplied in the column specification.
// New data types can be registered with the factory, allowing for easy extension of the table manager's 
// capabilities without modifying the core codebase.

class DataTypeFactory {
    dataTypes = {};

    constructor() {
        this.registerDataType("text", Text);
        this.registerDataType("numeric", Numeric);
        this.registerDataType("textselect", TextSelect);
        this.registerDataType("weblinktext", WebLinkText);
    }

    registerDataType(dataTypeName, dataTypeClass) {
        this.dataTypes[dataTypeName] = dataTypeClass;
    }

    initialiseDataType(Column) {
        if (!this.dataTypes[Column.data_type]) {
            throw new Error(`Data type "${Column.data_type}" is not registered.`);
        }
        return new this.dataTypes[Column.data_type](Column);
    }
}

export { DataTypeFactory }
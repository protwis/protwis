// This class represents the "weblink" data type used in the table manager.
// It extends the "text" data type, providing additional functionality for rendering and filtering web links in the table after deserialization. 
// The data for this data type is expected to be an object containing two properties: "anchor_href" and "anchor_content".


import { Text } from "./text.js";

class WebLinkText extends Text {  

    dataTableRenderer(data, type, row, meta) {
        if (type === "display" && data) {
            if (!data.anchor_href) {
                return data.anchor_content || ""; // fallback to displaying anchor_content as plain text if href is missing
            }
            return `<a href="${data.anchor_href}" target="_blank">${data.anchor_content}</a>`;
        }
        // Search, order and type return
        return data ? data.anchor_content : ""
    }

    blurHandler = (e, col_idx, tableManagerReference) => {
        function stripHtmlTags(str) {
            return str.replace(/<[^>]*>/g, '');
        }

        const currentSearchValue = tableManagerReference.dataTableReference.columns(col_idx).search()[0];
        
        let newSearchValue = e.target.value;
        
        const dataOptions = tableManagerReference.dataTableReference.column(col_idx).data(0).toArray().map((item) => stripHtmlTags(item.anchor_content));

        const useRegex = !dataOptions.includes(newSearchValue)

        if (useRegex) {
              newSearchValue = newSearchValue.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'); // escape regex special characters
              newSearchValue = ".*" + newSearchValue + ".*"; // wrap tag values with wildcards for "contains" search              
        } 

        if (newSearchValue !== currentSearchValue) {
            tableManagerReference.dataTableReference.columns(col_idx).search(newSearchValue, useRegex).draw();
        }
    }
    
}

export { WebLinkText };
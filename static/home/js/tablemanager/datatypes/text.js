// This class represents the "text" data type used in the table manager.
// It provides simple text filtering capabilities, allowing users to filter table data based on text input.

import { DataTypeBase } from "./datatypebase.js";

class Text extends DataTypeBase {

    createFilterInterface(tableManagerReference) {
        if (!this.allow_filter) {
            return null;
        }

        let filterInput = document.createElement("input");
        filterInput.setAttribute("data-field", this.json_id);
        filterInput.setAttribute("id", this.json_id + "_search");

        filterInput.classList.add("col-filter");
        filterInput.classList.add("text-filter");

        if (this.cssClassFilterInput) {
            filterInput.classList.add(this.cssClassFilterInput);
        }

        filterInput.addEventListener("blur", (e) => this.blurHandler(e, this.col_idx, tableManagerReference));
        filterInput.addEventListener("keydown", (e) => {
            if (e.key === "Enter") {
                this.blurHandler(e, this.col_idx, tableManagerReference);
            }
        });

        const wrapper = document.createElement("div");
        wrapper.style.cssText = "position:relative;display:inline-flex;align-items:center;gap:2px;";
        wrapper.appendChild(filterInput);

        if (this.show_symbol_popup) {
            const { symbolBtn, popup } = this.createSymbolPopup(filterInput);
            wrapper.appendChild(symbolBtn);
            wrapper.appendChild(popup);
        }

        return wrapper;
    }

    createSymbolPopup(filterInput) {
        const greekSymbols = [
            { label: "α", name: "alpha" },
            { label: "β", name: "beta" },
            { label: "δ", name: "delta" },
            { label: "κ", name: "kappa" },
            { label: "μ", name: "mu" },
        ];

        const popup = document.createElement("div");
        let popupCss = "display:none;"
        popupCss += "position:absolute;"
        popupCss += "z-index:9999;"
        popupCss += "background:#fff;"
        popupCss += "border:1px solid #ccc;"
        popupCss += "border-radius:4px;"
        popupCss += "padding:4px;"
        popupCss += "box-shadow:0 2px 6px rgba(0,0,0,0.2);"
        popupCss += "display:none;"
        popupCss += "grid-template-columns:repeat(3,1fr);"
        popupCss += "gap:3px;"
        popup.style.cssText = popupCss;

        greekSymbols.forEach(({ label, name }) => {
            const cell = document.createElement("button");
            cell.type = "button";
            cell.title = name;
            cell.textContent = label;
            let cellCss = "cursor:pointer;"
            cellCss += "border:1px solid #ddd;"
            cellCss += "border-radius:3px;"
            cellCss += "background:#f8f8f8;"
            cellCss += "padding:2px 5px;"
            cellCss += "font-size:1em;"
            cellCss += "line-height:1.4;"
            cell.style.cssText = cellCss;
            cell.addEventListener("mousedown", (e) => {
                e.preventDefault();
                filterInput.value += label;
                filterInput.focus();
            });
            popup.appendChild(cell);
        });

        const symbolBtn = document.createElement("button");
        symbolBtn.type = "button";
        symbolBtn.textContent = "Ω";
        let symbolBtnCss = "cursor:pointer;"
        symbolBtnCss += "border:1px solid #ccc;"
        symbolBtnCss += "border-radius:3px;"
        symbolBtnCss += "background:#f0f0f0;"
        symbolBtnCss += "padding:1px 4px;"
        symbolBtnCss += "font-size:0.8em;"
        symbolBtnCss += "white-space:nowrap;"
        symbolBtn.style.cssText = symbolBtnCss;

        let popupVisible = false;

        symbolBtn.addEventListener("click", (e) => {
            e.stopPropagation();
            popupVisible = !popupVisible;
            popup.style.display = popupVisible ? "grid" : "none";
        });

        document.addEventListener("click", () => {
            if (popupVisible) {
                popupVisible = false;
                popup.style.display = "none";
            }
        }, true);

        return { symbolBtn, popup };
    }

    blurHandler = (e, col_idx, tableManagerReference) => {
        function stripHtmlTags(str) {
            return str.replace(/<[^>]*>/g, '');
        }

        const currentSearchValue = tableManagerReference.dataTableReference.columns(col_idx).search()[0];

        let newSearchValue = e.target.value;

        const dataOptions = tableManagerReference.dataTableReference.column(col_idx).data(0).toArray().map((item) => stripHtmlTags(item));

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

export { Text };
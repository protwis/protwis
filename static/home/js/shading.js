class MyShadingTable {

	constructor(table, columns, darkGrey, lightGrey) {
        this.darkGrey = darkGrey;
		this.lightGrey = lightGrey;
        this.table = table;
        this.minValues = [];
        this.maxValues = [];
        this.setColumns(columns);
    }

    setColumns(columns) {
        this.columnsIndexes = new Set(columns.filter(word => typeof word === "number"));
        let columnsStrings = columns.filter(word => this.isValidSelector(word));
        this.columnsClasses = columnsStrings.filter(word => word[0] == ".");
        this.columnsIds = columnsStrings.filter(word => word[0] == "#");
    }

    getIndexes() {
        let indexes = new Set(this.columnsIndexes);

        for (let i = 0; i < this.columnsIds.length; i++) {
            this.table.column(this.columnsIds[i]).every(function() {
                indexes.add(this.index() + 1);
            });
        }
        for (let i = 0; i < this.columnsClasses.length; i++) {
            this.table.columns(this.columnsClasses[i]).every(function() {
                indexes.add(this.index() + 1);
            });
        }
        return indexes;
    }

    calculateLimitValues(data) {
        let columnsLength = this.table.columns().nodes().length;

        // minValues and maxValues will have the minimum and maximum
        // values for the selected columns, and null for the rest.
        this.minValues = new Array(columnsLength + 1).fill(null);
        this.maxValues = new Array(columnsLength + 1).fill(null);

        let indexes = this.getIndexes();
        for (let n of indexes) {
            this.minValues[n] = this.getMin(n, data);
            this.maxValues[n] = this.getMax(n, data);
        }
    }

    shade(rows) {
        let that = this;
        this.calculateLimitValues(rows.data());
        rows.every(function(){
            that.shadeMyTableRow(this.node(), this.data());
        });
	}

    isValidSelector(selector) {
        if (typeof(selector) !== "string")
            return false;

        try {
            var $element = $(selector);
        } catch(error) {
            return false;
        }

        return true;
    }

    getValidNumber(text, currentIndex) {
        if (!this.isValidSelector(text))
            return null;

        // This handles nested font elements properly.
        if ($(text).is("font"))
            text = $(text).text();

        if (text == null || text == "" || isNaN(text))
            return null;
        else
            return text;
    }

    getMax(numberColumn, data) {
        let that = this;
        return data.reduce(function(total, currentValue, currentIndex, arr) {
            let value = that.getValidNumber(currentValue[numberColumn-1], currentIndex);
            if (value)
                return Math.max(total, value);
            else
                return total;
		}, Number.NEGATIVE_INFINITY);
	}

	getMin(numberColumn, data) {
        let that = this;
        return data.reduce(function(total, currentValue, currentIndex, arr) {
            let value = that.getValidNumber(currentValue[numberColumn-1], currentIndex);
            if (value)
                return Math.min(total, value);
            else
                return total;
		}, Number.POSITIVE_INFINITY);
	}

	shadeMyTableRow(row, data, index) {
        for (let i = 0; i < this.minValues.length; i++)
            if (this.minValues[i] != null)
                this.paintCellOfColumn(i, this.minValues[i], this.maxValues[i], row, data);
            else
                this.paintCellOfColumn(i, 0, 0, row, data);
	}

	paintCellOfColumn(numberColumn, min, max, row, data) {
        if (max <= min) {
            $(row).find("td:eq(" + (numberColumn-1) + ")").css("background",
                        "none");
            return;
        }

        let value = this.getValidNumber(data[numberColumn-1], 0);
        if (!value)
            return;

		let factor = (value - min) / (max - min);
		let color = this.lightGrey - factor * (this.lightGrey - this.darkGrey);
        let n = Math.round(color);

        $(row).find("td:eq(" + (numberColumn-1) + ")").css("background",
                    "rgb(" + n + "," + n + "," + n + ")");
	}
}

// - The first argument is the table instance.
// - The second are the columns, e.g. [25, 26, ".niceColumn", "#firstColumn"].
// - The fourth and the fifth are the darkest and lightest grey colors:
//   valid values are from 0 to 255.
function shadeTable(table, columnIndexes, darkGrey, lightGrey)
{   let myShadingTable = null;
     myShadingTable = new MyShadingTable(table, columnIndexes, darkGrey, lightGrey);
     myShadingTable.shade(table.rows({page: "current"}));
     table.on("draw.dt", function(e, settings) {
        let api = new $.fn.dataTable.Api(settings);
        myShadingTable.shade(api.rows({page: "current"}));
    });


}

class MyShadingTable {

	constructor(table, columns, darkGrey, lightGrey) {
        this.darkGrey = darkGrey
		this.lightGrey = lightGrey
        this.table = table
        this.minValues = []
        this.maxValues = []
        this.columnsIndexes = new Set(columns.filter(word => typeof word === 'number'))

        let columnsStrings = columns.filter(word => this.isValidSelector(word))
        let columnsClasses = columnsStrings.filter(word => word[0] == '.')
        let columnsIds = columnsStrings.filter(word => word[0] == '#')
        let that = this

        for (let i = 0; i < columnsIds.length; i++) {
            this.table.column(columnsIds[i]).every(function() {
                that.columnsIndexes.add(this.index() + 1)
            });
        }

        for (let i = 0; i < columnsClasses.length; i++) {
            this.table.columns(columnsClasses[i]).every(function() {
                that.columnsIndexes.add(this.index() + 1)
            });
        }
    }

    calculateLimitValues(data) {
        if (this.columnsIndexes.size <= 0)
            return

        this.maximumColumn = Math.max.apply(null, [...this.columnsIndexes])

        // minValues and maxValues will have the minimum and maximum
        // values for the selected columns, and null for the rest.
        this.minValues = new Array(this.maximumColumn + 1).fill(null);
        this.maxValues = new Array(this.maximumColumn + 1).fill(null);

        for (let n of this.columnsIndexes) {
            this.minValues[n] = this.getMin(n, data)
            this.maxValues[n] = this.getMax(n, data)
        }
    }

    shade(rows) {
        let that = this
        this.calculateLimitValues(rows.data())
        rows.every(function(){
            that.shadeMyTableRow(this.node(), this.data())
        });
	}

    isValidSelector(selector) {
        if (typeof(selector) !== 'string')
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
        if ($(text).is('font'))
            text = $(text).text()

        if (text == null || text == '' || isNaN(text))
            return null
        else
            return text
    }

    getMax(numberColumn, data) {
        let that = this
        return data.reduce(function(total, currentValue, currentIndex, arr) {
            let value = that.getValidNumber(currentValue[numberColumn-1], currentIndex)
            if (value)
                return Math.max(total, value);
            else
                return total;
		}, Number.NEGATIVE_INFINITY);
	}

	getMin(numberColumn, data) {
        let that = this
        return data.reduce(function(total, currentValue, currentIndex, arr) {
            let value = that.getValidNumber(currentValue[numberColumn-1], currentIndex)
            if (value)
                return Math.min(total, value);
            else
                return total;
		}, Number.POSITIVE_INFINITY);
	}

	shadeMyTableRow(row, data, index){
        for (let i = 0; i < this.minValues.length; i++)
            if (this.minValues[i] != null)
                this.paintCellOfColumn(i, this.minValues[i], this.maxValues[i], row, data)
	}

	paintCellOfColumn(numberColumn, min, max, row, data){
        if (max <= min) {
            $(row).find('td:eq(' + (numberColumn-1) + ')').css('background',
                        'none');
            return;
        }

        let value = this.getValidNumber(data[numberColumn-1], 0)
        if (!value)
            return;

		let factor = (value - min) / (max - min)
		let color = this.lightGrey - factor * (this.lightGrey - this.darkGrey)
        let n = Math.round(color)

        $(row).find('td:eq(' + (numberColumn-1) + ')').css('background',
                    'rgb(' + n + ',' + n + ',' + n + ')');
	}
}

function shadeTable(table, columnIndexes, darkGrey, lightGrey)
{
    let myShadingTable = new MyShadingTable(table, columnIndexes, darkGrey, lightGrey);

    table.on('draw.dt', function(e, settings) {
        let api = new $.fn.dataTable.Api(settings)
        myShadingTable.shade(api.rows({page: 'current'}))
    });

    myShadingTable.shade(table.rows({page: 'current'}))
}

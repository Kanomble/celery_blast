var _table_ = document.createElement('table'),
    _tr_ = document.createElement('tr'),
    _th_ = document.createElement('th'),
    _td_ = document.createElement('td');


// Builds the HTML Table out of element list
function buildHtmlTable(arr) {
    var table = _table_.cloneNode(false),
    columns = addAllColumnHeaders(arr, table);
    for (var i = 0, maxi = arr.length; i < maxi; ++i) {
        var tr = _tr_.cloneNode(false);
        for (var j = 0, maxj = columns.length; j < maxj; ++j) {
            var td = _td_.cloneNode(false);
            cellValue = arr[i][columns[j]];

            if(arr[i][columns[j]].includes('href')){
                console.log(arr[i][columns[j]]);
                var link = document.createElement('a');
                link.innerHTML = arr[i][columns[j]];
                td.appendChild(link)
            } else {
                td.appendChild(document.createTextNode(arr[i][columns[j]] || ''));
            }

            tr.appendChild(td);
            }
        table.appendChild(tr);
        }
    return table;
}

// Adds a header row to the table and returns the set of columns.
// arr needs to be set up as dictionary
function addAllColumnHeaders(arr, table) {
    var columnSet = [],
    tr = _tr_.cloneNode(false);
    for (var i = 0, l = arr.length; i < l; i++) {
        for (var key in arr[i]) {
            if (arr[i].hasOwnProperty(key) && columnSet.indexOf(key) === -1) {
                columnSet.push(key);
                var th = _th_.cloneNode(false);
                th.appendChild(document.createTextNode(key));
                tr.appendChild(th);
                }
            }
        }

    table.appendChild(tr);
    return columnSet;
}

function createElementListForTable(jsonResponse){
    var elementList = [];
    for(var i = 0;i < jsonResponse["QUERIES"].length;i++) {
        //var pElem = document.createElement('p');
        //pElem.innerHTML = jsonResponse["QUERIES"][i];
        var links = {'TIGR': 'no entries', 'PFAM': 'no entries',
            'REFSEQ' : `<a href="${jsonResponse['REFSEQ'][jsonResponse["QUERIES"][i]]}">${jsonResponse['REFSEQ'][jsonResponse["QUERIES"][i]]}</a>`,
            'WP_NUMBER': jsonResponse["QUERIES"][i]}

        if (jsonResponse['TIGR'][jsonResponse["QUERIES"][i]]) {
            var link = `<a href="${jsonResponse['TIGR'][jsonResponse["QUERIES"][i]]}">${jsonResponse['TIGR'][jsonResponse["QUERIES"][i]]}</a>`
            links['TIGR'] = link;
        }

        if (jsonResponse['PFAM'][jsonResponse["QUERIES"][i]]) {
            var link = `<a href="${jsonResponse['PFAM'][jsonResponse["QUERIES"][i]]}">${jsonResponse['PFAM'][jsonResponse["QUERIES"][i]]}</a>`
            links['PFAM'] = link;
        }

        elementList.push(links)
    }

    return elementList;
}
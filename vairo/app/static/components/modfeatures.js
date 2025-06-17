const objectStoreName = 'data';
const nameDatabase = 'vairoFeaturesDatabase';
const objectKey = 'dataInfo';

async function runFeaturesInfo(configurationDict){
    let resultDict = {};
    try{
        resultDict = await postData('/read-features-info', configurationDict);
    } catch (error) {
        alert('It has not been possible to read features file');
    }
    return resultDict;
}

async function runDeleteFeatures(configurationDict){
    let resultDict = {};
    try{
        resultFile = await postData('/modify-pkl', configurationDict, false);
        if (resultFile) {
            const url = window.URL.createObjectURL(resultFile);
            const a = document.createElement('a');
            a.href = url;
            a.download = 'modified_features.pkl';
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            window.URL.revokeObjectURL(url);
            alert('File downloaded successfully');
        } else {
            alert('File data not found in the response');
        }
    } catch (error) {
        alert(error.message);
        alert('It has not been possible to modify features file');
    }
    return resultDict;
}

function createCoveragePlot(infoDict, regionsList){
    const generalPlotsDiv = document.getElementById("modfeatures-coverage-plot");
    const dataShape = [];
    const maxIndex = Math.max(infoDict.coverage.msa_coverage.length, infoDict.coverage.templates_coverage.length);

    function processCoverage(coverage, yOffset) {
        coverage.forEach((value, index) => {
            const color = gradient(value);
            dataShape.push({
                type: 'rect',
                x0: index + 1,
                x1: index + 2,
                y0: yOffset + 0.3,
                y1: yOffset + 0.9,
                fillcolor: color,
                line: {
                    width: 1,
                    color: color,
                },
            });
        });
    }
    processCoverage(infoDict.coverage.msa_coverage, 0);
    processCoverage(infoDict.coverage.templates_coverage, 0.9);
    dataShape.push({
        type: 'rect',
        x0: 1-maxIndex*0.03,
        x1: maxIndex+maxIndex*0.03,
        y0: 0,
        y1: 2.2,
        fillcolor: "rgb(251, 247, 230, 1)",
        layer: 'below',
        line: {
            width: 0
        }
    });
    regionsList.forEach(value => {
        if(value < maxIndex){
            dataShape.push({
                type: 'rect',
                x0: value,
                x1: value+1,
                y0: 0,
                y1: 2.2,
                layer: 'below',
                fillcolor: 'rgb(255, 181, 181)',
                line: {
                    width: 1,
                    color: 'rgb(255,181,181)',
                },
            });
        }
    });


    const layout = {
        title: 'Coverage of templates and MSA sequences',
        xaxis: {
            title: 'Residue number',
            range: [1-maxIndex*0.03, maxIndex+maxIndex*0.03],
            tickangle: 45,
            tick0: 1,
            dtick: (Math.ceil(maxIndex / 100) * 10 - 1),
            gridcolor: 'black',
            zeroline: false,
            showline: false,
        },
        yaxis: {
            title: 'Type',
            tickvals: [0.6, 1.5],
            ticktext: ['MSA', 'Templates'],
            range: [0, 2.2],
            showgrid: false,
            zeroline: false,
            showline: false,
            tickmode: 'array',
            line: {
                width: 0
            }
        },
        showlegend: false,
        shapes: dataShape
    };

    const config = {
        displayModeBar: false
    };
    Plotly.newPlot(generalPlotsDiv, [], layout, config);
}

function createDescriptionPlot(infoDict){
    const generalPlotsDiv = document.getElementById("modfeatures-description-plot");
    var keys = Object.keys(infoDict.uniprot_description_statistics);
    var values = Object.values(infoDict.uniprot_description_statistics);
    var descriptions = values.map(item => item.description);
    var identities = values.map(item => item.identity);
    var maxVal = Math.max(...identities);
    var colorScale = 'Viridis';
    var hoverText = values.map((val, index) => `Identity: ${identities[index]}, Description: ${descriptions[index]}`);

    var trace = {
        x: keys,
        y: descriptions,
        type: 'bar',
        text: descriptions.map(String),
        hovertext: hoverText,
        textposition: 'inside',
        marker: {
            color: identities,
            colorscale: colorScale,
            cmin: 0,
            cmax: maxVal,
            colorbar: {
                title: 'Identity',
            },
        },
    };
    var plotData = [trace];
    var layout = {
        autosize: true,
        title: 'Sequence descriptions',
        xaxis: {
            title: 'Description',
            tickangle: -45,
        },
        yaxis: {
              title: 'Total',
        },
    };
    Plotly.newPlot(generalPlotsDiv, plotData, layout);
}

function createOrganismPlot(infoDict){
    const generalPlotsDiv = document.getElementById("modfeatures-organism-plot");
    var keys = Object.keys(infoDict.uniprot_organism_statistics);
    var values = Object.values(infoDict.uniprot_organism_statistics);
    var organisms = values.map(item => item.organism);
    var identities = values.map(item => item.identity);
    var maxVal = Math.max(...identities);
    var colorScale = 'Viridis';
    var hoverText = values.map((val, index) => `Identity: ${identities[index]}, Description: ${organisms[index]}`);

    var trace = {
        x: keys,
        y: organisms,
        type: 'bar',
        text: organisms.map(String),
        hovertext: hoverText,
        textposition: 'inside',
        marker: {
            color: identities,
            colorscale: colorScale,
            cmin: 0,
            cmax: maxVal,
            colorbar: {
                title: 'Identity',
            },
        },
    };
    var plotData = [trace];
    var layout = {
        autosize: true,
        title: 'Sequence organism',
        xaxis: {
            title: 'Organism',
            tickangle: -45,
        },
        yaxis: {
              title: 'Total',
        },
    };
    Plotly.newPlot(generalPlotsDiv, plotData, layout);
}

function createTable(configurationDict, typeDiv){
    const sortedKeys = Object.keys(configurationDict).sort((a, b) => configurationDict[b].identity - configurationDict[a].identity);
    const table = document.createElement('table');
    table.classList.add('table', 'table-bordered', 'table-hover')
    const thead = document.createElement('thead');
    const headerRow = document.createElement('tr');
    const headers = ['ID', 'Global Identity (%)', 'Identity (%)', 'Coverage (%)'];
    headers.forEach(headerText => {
        const th = document.createElement('th');
        th.textContent = headerText;
        headerRow.appendChild(th);
    });
    thead.appendChild(headerRow);
    table.appendChild(thead);
    const tbody = document.createElement('tbody');
    sortedKeys.forEach(key => {
        const row = document.createElement('tr');
        row.classList.add('table-hover');
        row.style.cursor = 'pointer';
        const cell1 = document.createElement('td');
        cell1.appendChild(document.createTextNode(key));
        row.appendChild(cell1);
        const cell2 = document.createElement('td');
        cell2.appendChild(document.createTextNode(configurationDict[key].global_identity));
        row.appendChild(cell2);
        const cell3 = document.createElement('td');
        cell3.appendChild(document.createTextNode(configurationDict[key].identity));
        row.appendChild(cell3);
        const cell4 = document.createElement('td');
        cell4.appendChild(document.createTextNode(configurationDict[key].coverage));
        row.appendChild(cell4);
        row.classList.add('table-row-hover');
        row.setAttribute('data-key', key);
        row.addEventListener('click', handleRowClick);
        tbody.appendChild(row);
    });
    table.appendChild(tbody);
    typeDiv.appendChild(table);
    function handleRowClick(event) {
        const key = event.currentTarget.getAttribute('data-key');
        const url = `modfeaturesinfo?key=${key}`;
        window.open(url, '_blank');
    }
}

async function updateModFeatures(){
    document.getElementById("modfeatures-run-button").setAttribute("disabled", "true");

    const fileFeatures = document.getElementById('modfeatures-pkl')?.files?.[0];
    const rangeFeatures = document.getElementById('modfeatures-range').value;
    const iniIdentity = document.getElementById('modfeatures-indentity-ini').value;
    const endIdentity = document.getElementById('modfeatures-indentity-end').value;
    const runUniprot = document.getElementById('modfeatures-uniprot').checked;
    inputDict = {
        'fileFeatures': fileFeatures,
        'rangeFeatures': rangeFeatures,
        'iniIdentity': iniIdentity,
        'endIdentity': endIdentity,
        'runUniprot': runUniprot
    }
    resultsFeatures = await runFeaturesInfo(inputDict);
    const regionsList = extendedNumbers(rangeFeatures);
    document.getElementById('modfeatures-templates-table').innerHTML = '';
    document.getElementById('modfeatures-msa-table').innerHTML = '';
    document.getElementById('modfeatures-download-button').innerHTML = '';
    document.getElementById('modfeatures-general-info').innerHTML = '';
    document.getElementById('modfeatures-description-plot').innerHTML = '';
    document.getElementById('modfeatures-organism-plot').innerHTML = '';
    document.getElementById('modfeatures-coverage-plot').innerHTML = '';

    if(Object.keys(resultsFeatures).length){
        indexedDB.deleteDatabase(nameDatabase);
        const idb = indexedDB.open(nameDatabase, 1);
        idb.onupgradeneeded = function (evt) {
            const dbobject = evt.target.result;
            dbobject.createObjectStore(objectStoreName);
        }
        idb.onsuccess = function(evt) {
            const dbobject = evt.target.result;
            const transaction = dbobject.transaction([objectStoreName], 'readwrite');
            const objectStore = transaction.objectStore(objectStoreName);
            const new_object = Object.assign({}, resultsFeatures['templates'], resultsFeatures['msa']);
            const request = objectStore.add(new_object, objectKey);
            createCoveragePlot(resultsFeatures, regionsList);
            createDescriptionPlot(resultsFeatures);
            createOrganismPlot(resultsFeatures);
            document.getElementById('accordion-table').removeAttribute('hidden');
            document.getElementById('modfeatures-templates-text').textContent = `TEMPLATES (${resultsFeatures.coverage.num_templates} models)`;
            document.getElementById('modfeatures-msa-text').textContent = `MSA (${resultsFeatures.coverage.num_msa} sequences)`;
            createTable(resultsFeatures['templates'], document.getElementById('modfeatures-templates-table'))
            createTable(resultsFeatures['msa'], document.getElementById('modfeatures-msa-table'))

            const generalMessage = `
                <p class="text-break">The full query sequence is:<br> ${resultsFeatures.general_information.query_sequence}</p>
                <p class="text-break">The search region of the query sequence is:<br> ${resultsFeatures.general_information.query_search}<p>
            `;
            createDownloadButton(resultsFeatures);
            document.getElementById("modfeatures-general-info").innerHTML = generalMessage;
            document.getElementById("modfeatures-run-button").removeAttribute("disabled");
        }
    }
}

// Function to dynamically create a button and append it to the DOM
function createDownloadButton(resultsFeatures) {
    const button = document.createElement('button');
    button.id = 'modfeatures-download-json';
    button.textContent = 'Download JSON results';
    document.getElementById('modfeatures-download-button').appendChild(button);
    button.addEventListener('click', function() {
        const jsonString = JSON.stringify(resultsFeatures);
        const blob = new Blob([jsonString], { type: 'application/json' });
        const link = document.createElement('a');
        link.href = URL.createObjectURL(blob);
        link.download = 'results.json';
        document.body.appendChild(link);
        link.click();
    });
}


function updateModFeaturesInfo(key, modFeaturesParams) {
    const modFeaturesInfoDiv = document.getElementById('modfeaturesinfo-data');
    modFeaturesInfoDiv.innerHTML = `
        <p><strong>ID:</strong> ${key}</p>
        <p><strong>Global Identity (%):</strong> ${modFeaturesParams.global_identity}</p>
        <p class="text-break"><strong>Identity (%):</strong> ${modFeaturesParams.identity}</p>
        <p class="text-break"><strong>Coverage (%):</strong> ${modFeaturesParams.coverage}</p>
        <p class="text-break"><strong>Full MSA sequence:</strong><br> ${modFeaturesParams.seq.replace(/-/g, "\u2011")}</p>
        <p class="text-break"><strong>Region MSA sequence:</strong><br> ${modFeaturesParams.seq_msa.replace(/-/g, "\u2011")}</p>
        <p class="text-break"><strong>Query sequence search:</strong><br> ${modFeaturesParams.seq_query.replace(/-/g, "\u2011")}</p>`

    if(modFeaturesParams.hasOwnProperty("uniprot")){
        modFeaturesInfoDiv.innerHTML += "<p><strong>UniProt information </strong></p>";
        for (let i = 0; i < modFeaturesParams.uniprot.length; i++) {
            const entry = modFeaturesParams.uniprot[i];
            modFeaturesInfoDiv.innerHTML += `<ul>
                        <li><strong>Accession ID:</strong> <a href="https://www.uniprot.org/uniprotkb/${entry.uniprot_accession_id}/entry">${entry.uniprot_accession_id}</a></li>                    <ul>
                        <li><strong>Annotation Score:</strong> ${entry.uniprot_annotation_score}</li>
                        <li><strong>E-value:</strong> ${entry.uniprot_evalue}</li>
                        <li><strong>Identity:</strong> ${entry.uniprot_identity}</li>
                        <li><strong>Organism:</strong> ${entry.uniprot_organism}</li>
                        <li><strong>Protein Description:</strong> ${entry.uniprot_protein_description}</li>
                    </ul>
                </li>
            </ul>`;
        }
    }
    if(modFeaturesParams.hasOwnProperty("pdb")){
        const stage = new NGL.Stage('viewport');
        const file = new Blob([modFeaturesParams["pdb"]], { type: 'text/plain' });
        NGL.autoLoad(file, { ext: "pdb" }).then(function(structure) {
            const o = stage.addComponentFromObject(structure)
            stage.defaultFileRepresentation(o)
        });
    }
}

async function deleteModFeatures(){
    document.getElementById("modfeatures-delete-button").setAttribute("disabled", "true");
    const fileFeatures = document.getElementById('modfeatures-pkl')?.files?.[0];
    const deleteList = document.getElementById('modfeatures-remove-list').value;
    const maxIdentity = document.getElementById('modfeatures-identity-max').value;
    const minIdentity = document.getElementById('modfeatures-identity-min').value;
    inputDict = {
        'fileFeatures': fileFeatures,
        'deleteList': deleteList,
        'iniIdentity': maxIdentity,
        'endIdentity': minIdentity
    }
    resultsFeatures = await runDeleteFeatures(inputDict);
    document.getElementById("modfeatures-delete-button").removeAttribute("disabled");
}

document.addEventListener("DOMContentLoaded", function() {
    if (window.location.pathname.endsWith("modfeaturesinfo")) {
        const urlParams = new URLSearchParams(window.location.search);
        const key = urlParams.get("key");
        const idb = indexedDB.open(nameDatabase, 1);
        idb.onsuccess = function (evt) {
            const dbobject = evt.target.result;
            const transaction = dbobject.transaction([objectStoreName], 'readonly');
            const objectStore = transaction.objectStore(objectStoreName);
            const getRequest = objectStore.get(objectKey);
            getRequest.onsuccess = function(event) {
                const dataInfo = event.target.result;
                console.log(dataInfo);
                updateModFeaturesInfo(key, dataInfo[key]);
            };
        }
    }
});
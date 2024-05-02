const featuresDict = {}

async function readFeatures(id){
    const featuresFile = document.getElementById(`feature-pkl-${id}`)?.files?.[0];
    try{
        const resultDict = await postData('/read-pkl', {'featuresFile': featuresFile});
        featuresDict[id] = resultDict;
    } catch (error) {
        alert('It has not been possible to read features');
        console.error('Error:', error);
        featuresDict[id] = {};
    }
    updatePlot();
}

function checkFeaturesMSA(id){
    const check = document.getElementById(`feature-addmsa-${id}`).checked;
    const divHide = document.getElementById(`feature-mask-div-${id}`);
    divHide.classList.toggle('hidden', !check);
    updatePlot();
}

function checkFeaturesTemplates(id){
    const check = document.getElementById(`feature-addtemplates-${id}`).checked;
    const divHide = document.getElementById(`feature-fasta-div-${id}`);
    divHide.classList.toggle('hidden', !check);
    updatePlot();
}

function deleteFeatures(id){
    delete featuresDict[id];
    updatePlot();
}

class featureTable extends HTMLElement {
    static formAssociated = true;
    static observedAttributes = ['value'];

    constructor() {
        super();
        this.attachInternals();
        this.featureID = this.getAttribute('featureID');
    }

    connectedCallback() {
        this.render();
    }

    render() {
        const optionsHTML = selectPositionsArray.map(option => `<option value="${option}">${option}</option>`).join('');
        this.innerHTML =  `            
            <fieldset name="feature-field" class="row g-3"> 
                <div class="row row-margin">
                    <div class="form-group">
                        <label class="form-label" for="feature-pkl-${this.featureID}">Select an AlphaFold2 input file</label>
                        <input type="file" accept=".pkl" class="form-control" name="feature-pkl-${this.featureID}" id="feature-pkl-${this.featureID}" title="Choose an AlphaFold2 input file (pkl format)" onchange="readFeatures('${this.featureID}')">
                    </div>
                </div>
                <div class="row row-margin" id="feature-add-templates-div-${this.featureID}">
                    <div class="col-md-auto">
                        <input type="checkbox" id="feature-addtemplates-${this.featureID}" name="feature-addtemplates-${this.featureID}" value="true" onchange="checkFeaturesTemplates('${this.featureID}')" checked>
                        <label class="form-label" for="feature-addtemplates-${this.featureID}"> Add to templates</label>
                    </div>
                </div>
                <div class="row">
                    <div class="col-md-auto">
                        <input type="checkbox" id="feature-addmsa-${this.featureID}" name="feature-addmsa-${this.featureID}" value="true" onchange="checkFeaturesMSA('${this.featureID}')" checked>
                        <label class="form-label" for="feature-addmsa-${this.featureID}"> Add to MSA</label>
                    </div>
                </div>
                <div class="row row-margin">
                    <div class="col-md-auto">
                        <label class="form-label" for="feature-pos-${this.featureID}">Insert into query sequence position</label>
                        <select class="form-select" id="feature-pos-${this.featureID}" name="feature-pos-${this.featureID}" title="Select the position in the query sequence" onchange="updatePlot()">
                            ${optionsHTML}
                        </select>
                    </div>
                </div>
                <div class="row row-margin">
                    <div class="col-md-auto">
                        <label class="form-label" for="feature-regionfeat-${this.featureID}">Select regions from features</label>
                        <input type="text" class="form-control" id="feature-regionfeat-${this.featureID}" name="feature-regionfeat-${this.featureID}" placeholder="1, 3-10" title="Select regions from features" onchange="updatePlot()">
                    </div>
                    <div class="col-md-auto">
                        <label class="form-label" for="feature-regionquery-${this.featureID}">Insert those regions in the following positions of the query sequence</label>
                        <input type="text" class="form-control" id="feature-regionquery-${this.featureID}" name="feature-regionquery-${this.featureID}" placeholder="Integer numbers separated by comma e.g. 1,3" title="Insert those regions in the following positions of the query sequence" onchange="updatePlot()">
                    </div>
                </div>
                <div class="row row-margin" id="feature-fasta-div-${this.featureID}">
                    <div class="col-md-auto">
                        <label class="form-label" for="feature-fasta-${this.featureID}">Replace templates sequence with the following sequence</label>
                        <input class="form-control" type="file" accept=".fasta, .seq, .sequence" name="feature-fasta-${this.featureID}" id="feature-fasta-${this.featureID}" title="Choose sequence file (fasta format) to replace the templates sequences" onchange="updatePlot()">
                    </div>
                </div>
                <div class="row row-margin" id="feature-mask-div-${this.featureID}">
                    <div class="col-md-auto">
                        <label class="form-label" for="feature-mask-${this.featureID}">Mask residues in the MSA</label>
                        <input type="text" class="form-control" id="feature-mask-${this.featureID}" name="feature-mask-${this.featureID}" placeholder="1, 3-10" title="Select the residues to mask of the features">
                    </div>
                </div>
            </fieldset>
        `;
    }
}
customElements.define('feature-component', featureTable);
const featuresDict = {};

class FeatureTable extends HTMLElement {
    static formAssociated = true;
    static observedAttributes = ['value'];

    constructor() {
        super();
        this.internals = this.attachInternals();
        this.featureID = this.getAttribute('feature-id');
    }

    connectedCallback() {
        this.render();
        this.attachEventListeners();
    }

    async readFeatures(event) {
        const id = this.featureID;
        const featuresFile = event.target.files?.[0];
        try {
            const resultDict = await postData('/read-pkl', { featuresFile });
            featuresDict[id] = resultDict;
        } catch (error) {
            alert('It has not been possible to read features');
            console.error('Error:', error);
            featuresDict[id] = {};
        }
        updatePlot();
    }

    checkFeaturesMSA(event) {
        const id = this.featureID;
        const check = event.target.checked;
        const divHide = this.querySelector(`#feature-mask-div-${id}`);
        divHide.classList.toggle('hidden', !check);
        updatePlot();
    }

    checkFeaturesTemplates(event) {
        const id = this.featureID;
        const check = event.target.checked;
        const divHide = this.querySelector(`#feature-fasta-div-${id}`);
        divHide.classList.toggle('hidden', !check);
        updatePlot();
    }

    deleteFeatures() {
        delete featuresDict[this.featureID];
        updatePlot();
    }

    togglePositionFeatures(event) {
        const featureID = this.featureID;
        const posSection = this.querySelector(`#pos-section-${featureID}`);
        const regionSection = this.querySelector(`#region-section-${featureID}`);
        const posRadio = this.querySelector(`#mode-pos-${featureID}`);
        if (posRadio.checked) {
            posSection.style.display = '';
            regionSection.style.display = 'none';
        } else {
            posSection.style.display = 'none';
            regionSection.style.display = '';
        }
    }

    attachEventListeners() {
        const id = this.featureID;

        this.querySelector(`#feature-pkl-${id}`).addEventListener('change', this.readFeatures.bind(this));
        this.querySelector(`#feature-addmsa-${id}`).addEventListener('change', this.checkFeaturesMSA.bind(this));
        this.querySelector(`#feature-addtemplates-${id}`).addEventListener('change', this.checkFeaturesTemplates.bind(this));
        this.querySelector(`#mode-pos-${id}`).addEventListener('change', this.togglePositionFeatures.bind(this));
        this.querySelector(`#mode-region-${id}`).addEventListener('change', this.togglePositionFeatures.bind(this));
    }

    render() {
        const id = this.featureID;
        const optionsHTML = (window.selectPositionsArray || []).map(option =>`<option value="${option}">${option}</option>`).join('');
        this.innerHTML =  `
            <fieldset name="feature-field" class="row g-3">
                <div class="row row-margin">
                    <div class="form-group">
                        <label class="form-label" for="feature-pkl-${id}">Select an AlphaFold2 input file</label>
                        <input type="file" accept=".pkl" class="form-control" name="feature-pkl-${id}" id="feature-pkl-${id}" title="Choose an AlphaFold2 input file (pkl format)" aria-describedby="feature-pkl-desc-${id}">
                        <small id="feature-pkl-desc-${id}" class="form-text text-muted">
                            Choose an AlphaFold2 input file in PKL format.
                        </small>
                    </div>
                </div>
                <div class="row row-margin" id="feature-add-templates-div-${id}">
                    <div class="col-md-auto">
                        <input type="checkbox" id="feature-addtemplates-${id}" name="feature-addtemplates-${id}" value="true" checked>
                        <label class="form-label" for="feature-addtemplates-${id}"> Add to templates</label>
                    </div>
                </div>
                <div class="row">
                    <div class="col-md-auto">
                        <input type="checkbox" id="feature-addmsa-${id}" name="feature-addmsa-${id}" value="true" checked>
                        <label class="form-label" for="feature-addmsa-${id}"> Add to MSA</label>
                    </div>
                </div>
                <div class="row row-margin">
                  <div class="col-md-auto">
                    <div class="form-check form-check-inline">
                      <input class="form-check-input" type="radio" name="inputMode-${id}" id="mode-pos-${id}" value="pos" checked>
                      <label class="form-check-label" for="mode-pos-${id}">Insert by position</label>
                    </div>
                    <div class="form-check form-check-inline">
                      <input class="form-check-input" type="radio" name="inputMode-${id}" id="mode-region-${id}" value="region">
                      <label class="form-check-label" for="mode-region-${id}">Insert by region</label>
                    </div>
                  </div>
                </div>
                <div class="row row-margin" id="pos-section-${id}">
                  <div class="col-md-auto">
                    <label class="form-label" for="feature-pos-${id}">Insert into query sequence position</label>
                    <select class="form-select" id="feature-pos-${id}" name="feature-pos-${id}" title="Select the position in the query sequence" onchange="updatePlot()">
                      ${optionsHTML}
                    </select>
                  </div>
                </div>
                <div class="row row-margin" id="region-section-${id}" style="display: none;">
                  <div class="col-md-auto">
                    <label class="form-label" for="feature-regionfeat-${id}">Select regions from features</label>
                    <input type="text" class="form-control" id="feature-regionfeat-${id}" name="feature-regionfeat-${id}" placeholder="1, 3-10" title="Select regions from features" onchange="updatePlot()">
                  </div>
                  <div class="col-md-auto">
                    <label class="form-label" for="feature-regionquery-${id}">Insert those regions in the following query sequence numbering</label>
                    <input type="text" class="form-control" id="feature-regionquery-${id}" name="feature-regionquery-${id}" placeholder="Integer numbers separated by comma e.g. 1,3" title="Insert those regions in the following query sequence numbering" onchange="updatePlot()">
                  </div>
                </div>
                <div class="row row-margin" id="feature-fasta-div-${id}">
                    <div class="col-md-auto">
                        <label class="form-label" for="feature-fasta-${id}">Replace templates sequence with the following sequence</label>
                        <input class="form-control" type="file" accept=".fasta, .seq, .sequence" name="feature-fasta-${id}" id="feature-fasta-${id}" title="Choose sequence file (fasta format) to replace the templates sequences" onchange="updatePlot()" aria-describedby="feature-fasta-desc-${id}">
                        <small id="feature-fasta-desc-${id}" class="form-text text-muted">
                            Choose a fasta file.
                        </small>
                    </div>
                </div>
                <div class="row row-margin" id="feature-mask-div-${id}">
                    <div class="col-md-auto">
                        <label class="form-label" for="feature-mask-${id}">Mask residues in the MSA</label>
                        <input type="text" class="form-control" id="feature-mask-${id}" name="feature-mask-${id}" placeholder="1, 3-10" title="Select the residues to mask of the features">
                    </div>
                </div>
            </fieldset>
        `;
    }
}

customElements.define('feature-component', FeatureTable);

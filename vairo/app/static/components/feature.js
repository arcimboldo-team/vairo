class FeatureTable extends HTMLElement {
  static formAssociated = true;
  static observedAttributes = ['value'];

  constructor() {
    super();
    this.internals = this.attachInternals();
    this.featureID = this.getAttribute('feature-id');
    this.features = {};

    this.boundHandlers = {
      readFeatures: this.readFeatures.bind(this),
      checkFeaturesMSA: this.checkFeaturesMSA.bind(this),
      checkFeaturesTemplates: this.checkFeaturesTemplates.bind(this),
      togglePositionFeatures: this.togglePositionFeatures.bind(this),
      handleInputChange: this.handleInputChange.bind(this),
    };

    this.handlers = [
      ['fileInput', 'change', 'readFeatures'],
      ['addmsa', 'change', 'checkFeaturesMSA'],
      ['addtemplates', 'change', 'checkFeaturesTemplates'],
      ['modePos', 'change', 'togglePositionFeatures'],
      ['modeRegion', 'change', 'togglePositionFeatures'],
      ['positionSelect', 'change', 'handleInputChange'],
      ['featureRegionfeat', 'change', 'handleInputChange'],
      ['featureRegionquery', 'change', 'handleInputChange'],
      ['featureFasta', 'change', 'handleInputChange'],
    ];
  }

  cacheElements() {
    const id = this.featureID;
    this._elements = {
      fileInput: this.querySelector(`#feature-pkl-${id}`),
      addmsa: this.querySelector(`#feature-addmsa-${id}`),
      addtemplates: this.querySelector(`#feature-addtemplates-${id}`),
      modePos: this.querySelector(`#mode-pos-${id}`),
      modeRegion: this.querySelector(`#mode-region-${id}`),
      positionSelect: this.querySelector(`#feature-pos-${id}`),
      featureRegionfeat: this.querySelector(`#feature-regionfeat-${id}`),
      featureRegionquery: this.querySelector(`#feature-regionquery-${id}`),
      featureFasta: this.querySelector(`#feature-fasta-${id}`),
      maskDiv: this.querySelector(`#feature-mask-div-${id}`),
      fastaDiv: this.querySelector(`#feature-fasta-div-${id}`),
      posSection: this.querySelector(`#pos-section-${id}`),
      regionSection: this.querySelector(`#region-section-${id}`),
      posRadio: this.querySelector(`#mode-pos-${id}`),
    };
  }

  connectedCallback() {
    this.render();
    this.cacheElements();
    this.attachEventListeners();
  }

  disconnectedCallback() {
    if (!this._elements) return;
    this.handlers.forEach(([key, evt, handler]) => {
      const el = this._elements[key];
      if (el && typeof this.boundHandlers[handler] === 'function') {
        el.removeEventListener(evt, this.boundHandlers[handler]);
      }
    });
  }

  attachEventListeners() {
    if (!this._elements) return;
    this.handlers.forEach(([key, evt, handler]) => {
      const el = this._elements[key];
      if (el && typeof this.boundHandlers[handler] === 'function') {
        el.addEventListener(evt, this.boundHandlers[handler]);
      }
    });
  }

  async readFeatures(event) {
    const id = this.featureID;
    const featuresFile = event.target.files?.[0];
    if (!featuresFile) return;
    try {
      const resultDict = await postData('/read-pkl', { featuresFile });
      featuresDict[id] = resultDict;
      this.features = resultDict;
    } catch (error) {
      this.showError('It has not been possible to read features');
      featuresDict[id] = {};
      this.features = {};
    }
    this.triggerUpdatePlot();
  }

  checkFeaturesMSA(event) {
    const check = event.target.checked;
    this._elements.maskDiv?.classList.toggle('hidden', !check);
    this.triggerUpdatePlot();
  }

  checkFeaturesTemplates(event) {
    const check = event.target.checked;
    this._elements.fastaDiv?.classList.toggle('hidden', !check);
    this.triggerUpdatePlot();
  }

  togglePositionFeatures() {
    const { posRadio, posSection, regionSection } = this._elements;
    if (posRadio?.checked) {
      posSection.style.display = '';
      regionSection.style.display = 'none';
    } else {
      posSection.style.display = 'none';
      regionSection.style.display = '';
    }
    this.triggerUpdatePlot();
  }

  handleInputChange() {
    this.triggerUpdatePlot();
  }

  triggerUpdatePlot() {
    updatePlot()
  }

  render() {
    const id = this.featureID;
    const optionsHTML = selectPositionsArray.map(option => `<option value="${option}">${option}</option>`).join('');
    this.innerHTML = `
      <fieldset name="feature-field">
          <div class="form-group mb-2">
            <label for="feature-pkl-${id}">PKL file</label>
            <input type="file" accept=".pkl" class="form-control" name="feature-pkl-${id}" id="feature-pkl-${id}" title="Choose an AlphaFold2 input file (pkl format)" aria-describedby="feature-pkl-desc-${id}">
            <small id="feature-pkl-desc-${id}" class="form-text text-muted">
              Choose an AlphaFold2 input file in PKL format.
            </small>
          </div>
          <div class="form-check mb-2" id="feature-add-templates-div-${id}">
            <input class="form-check-input" type="checkbox" id="feature-addtemplates-${id}" name="feature-addtemplates-${id}" value="true" checked>
            <label class="form-check-label" for="feature-addtemplates-${id}"> Add to templates</label>
          </div>
          <div class="form-check mb-2">
            <input class="form-check-input" type="checkbox" id="feature-addmsa-${id}" name="feature-addmsa-${id}" value="true" checked>
            <label class="form-check-label" for="feature-addmsa-${id}"> Add to MSA</label>
          </div>
          <div class="form-group mb-2">
            <div class="form-check form-check-inline">
              <input class="form-check-input" type="radio" name="inputMode-${id}" id="mode-pos-${id}" value="pos" checked>
              <label class="form-check-label" for="mode-pos-${id}">Insert by position</label>
            </div>
            <div class="form-check form-check-inline">
              <input class="form-check-input" type="radio" name="inputMode-${id}" id="mode-region-${id}" value="region">
              <label class="form-check-label" for="mode-region-${id}">Insert by region</label>
            </div>
          </div>
          <div class="form-group mb-2" id="pos-section-${id}">
            <label for="feature-pos-${id}">Insert into query sequence position</label>
                <select class="form-select"
                        style="width: 80px;"
                        id="feature-pos-${id}"
                        name="feature-pos-${id}"
                        title="Select the position in the query sequence">
                  ${optionsHTML}
                </select>
          </div>
        <div class="form-group mb-2 row" id="region-section-${id}" style="display: none;">
          <div class="col-md-auto">
            <label for="feature-regionfeat-${id}">Select regions from features</label>
            <input type="text" class="form-control" id="feature-regionfeat-${id}" aria-describedby="feature-regionfeat-desc-${id}" name="feature-regionfeat-${id}" placeholder="e.g. 1, 3-10" title="Select specific regions from the features file">
            <small id="feature-regionfeat-desc-${id}" class="form-text text-muted">
              Select specific regions from the features file.
            </small>
          </div>
          <div class="col-md-auto">
            <label for="feature-regionquery-${id}">Insert those regions in the following query sequence numbering</label>
            <input type="text" class="form-control" id="feature-regionquery-${id}" name="feature-regionquery-${id}" placeholder="e.g.  1,3" aria-describedby="feature-regionquery-desc-${id}" title="Insert those regions in the corresponding query sequence numbering">
            <small id="feature-regionquery-desc-${id}" class="form-text text-muted">
              Insert those regions in the corresponding query sequence numbering.
            </small>
          </div>
        </div>

        <div class="form-group mb-2" id="feature-fasta-div-${id}">
            <label for="feature-fasta-${id}">Replace templates sequence with the following sequence</label>
            <input class="form-control" type="file" accept=".fasta, .seq, .sequence" name="feature-fasta-${id}" id="feature-fasta-${id}" aria-describedby="feature-fasta-desc-${id}" title="Choose sequence file (fasta format) to replace the templates sequences">
            <small id="feature-fasta-desc-${id}" class="form-text text-muted">
              Choose a fasta file.
            </small>
        </div>

          <div class="form-group mb-2" id="feature-mask-div-${id}">
            <label for="feature-mask-${id}">Mask residues in the MSA</label>
            <input type="text" class="form-control" id="feature-mask-${id}" name="feature-mask-${id}" placeholder="e.g. 1, 3-10" aria-describedby="feature-mask-desc-${id}" title="Select the numbering of residues to mask">
            <small id="feature-mask-desc-${id}" class="form-text text-muted">
              Select the numbering of residues to mask.
            </small>
          </div>
      </fieldset>
    `;
  }
}

customElements.define('feature-component', FeatureTable);
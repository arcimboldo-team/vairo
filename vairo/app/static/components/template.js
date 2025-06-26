class TemplateComponent extends HTMLElement {
    static formAssociated = true;
    static observedAttributes = ['value'];

    constructor() {
        super();
        this.attachInternals();
        this.templateID = this.getAttribute('template-id');
        this.modificationID = 1;
        this.templateName = "";
        this.templateData = "";
        this.chainSeq = {};
    }

    connectedCallback() {
        this.render();
        this.cacheElements();
        this.attachEventListeners();
    }

      cacheElements() {
        const id = this.templateID;
        this._elements = {
            radioFile: this.querySelector(`input[type="radio"][value="file"]`),
            fileSection: this.querySelector(`#template-file-section-${id}`),
            codeSection: this.querySelector(`#template-code-section-${id}`),
            templateCode: this.querySelector(`#template-code-${id}`),
            templateFile: this.querySelector(`#template-file-${id}`),
            templateMultimer: this.querySelector(`#template-multimer-${id}`),

        };
      }

    attachEventListeners() {
        if (!this._elements) return;
        this.querySelector(`#add-modification-${this.templateID}`).addEventListener('click', () => {
            this.addModification();
        });

        document.querySelectorAll(`[name="template-radio-${this.templateID}"]`).forEach(radio => {
            radio.addEventListener("change", async () => {
                this.handleRadioChange();
                await this.readTemplate();
            });
        });

        this._elements.templateCode.addEventListener("change", async () => {
            await this.readTemplate();
        });

        this._elements.templateFile.addEventListener("change", async () => {
            await this.readTemplate();
        });

        this._elements.templateMultimer.addEventListener("change", async () => {
            await this.generateMultimer();
        });
    }
      handleRadioChange() {
        const id = this.templateID;
        const { radioFile, fileSection, codeSection } = this._elements;
        if (radioFile?.checked) {
          fileSection.style.display = '';
          codeSection.style.display = 'none';
        } else {
          fileSection.style.display = 'none';
          codeSection.style.display = '';
        }
      }

    async readTemplate() {
        const id = this.templateID;
        if (this._elements.radioFile?.checked) {
            const pdbFileInput = this._elements.templateFile?.files?.[0];
            if (pdbFileInput === null || pdbFileInput === undefined) return;
            this.templateName  = pdbFileInput.name.split('.')[0];
            if (id in templatesDict && templatesDict[id].templateName === this.templateName) return;
            this.templateData = await readFile(pdbFileInput);
        } else {
            this.templateName = this._elements.templateCode.value;
            if ((id in templatesDict && templatesDict[id].templateName === this.templateName) || (this.templateName === "")) return;
            const fetchPromise = fetchPDB(this.templateName);
            const timeoutPromise = new Promise((resolve, reject) => {
                setTimeout(() => reject(new Error('Timeout')), 60000);
            });
            try {
                this.templateData = await Promise.race([fetchPromise, timeoutPromise]);
            } catch (error) {
                console.error('Error fetching PDB:', error);
                alert('The request has timed out. Please try again later.');
                this._elements.templateCode.value = "";
                return;
            }
        }
        this.generateMultimer();
    }

    async generateMultimer() {
        const id = this.templateID;
        if(!templatesDict.hasOwnProperty(id) && this.templateData === "undefined") return;
        if(this.templateData === "undefined"){
            this.templateData = templatesDict[id].templateData;
            this.templateName = templatesDict[id].templateName;
        }
        this.chainSeq = getSequences(this.templateData);
        for (let key in this.chainSeq) {
            this.chainSeq[key] = [this.chainSeq[key]];
        }
        if(this._elements.templateMultimer.checked){
            try{
                this.chainSeq = await postData('/generate-multimer', {'templateData': this.templateData});
            } catch (error) {
                alert('It has not been possible to generate the multimer. Unchecking it');
                check.checked = false;
                console.error('Error:', error);
            }
        }
        templatesDict[id] = {"templateName": this.templateName, "templateData": this.templateData, "chainSequences": this.chainSeq};
        this.triggerUpdatePlot();
        this.populateChainSelect();
    }

    populateChainSelect() {
        const id = this.templateID;
        let options = '';
        if(templatesDict.hasOwnProperty(id)){
            const chains = templatesDict[id].chainSequences;
            const uniqueChains = Object.keys(chains)
            options = `<option value="all">All chains</option>` + uniqueChains.map(chain => `<option value="${chain}">${chain} (${chains[chain].length} copies)</option>`).join('');
        } else {
            options = `<option value="all">All chains</option>`;
        }
        const modifications = document.querySelectorAll(`select[id^=template-modify-where-${id}-]`);
        modifications.forEach(modification => {
            const oldValue = modification.value;
            modification.innerHTML = options;
            const optionExists = Array.from(modification.options).some(option => option.value === oldValue);
            if(optionExists){
                modification.value = oldValue;
            }
        });
    }

    addModification() {
        const modificationsContainer = this.querySelector(`#modifications-container-${this.templateID}`);
        const modificationComponent = document.createElement('modification-component');
        modificationComponent.setAttribute('template-id', this.templateID);
        modificationComponent.setAttribute('modification-id', this.modificationID);
        modificationsContainer.appendChild(modificationComponent);
        this.populateChainSelect()
        ++this.modificationID;
    }

    triggerUpdatePlot() {
        updatePlot();
    }

    render() {
        const id = this.templateID;
        this.innerHTML = `
            <fieldset name="template-field">
                <div class="form-group mb-2">
                  <label class="me-3">Insert template PDB from</label>
                <div class="form-check form-check-inline">
                  <input
                    class="form-check-input"
                    type="radio"
                    name="template-radio-${id}"
                    id="template-radio-file-${id}"
                    value="file"
                    checked
                  >
                  <label class="form-check-label" for="template-radio-file-${id}">File</label>
                </div>
                <div class="form-check form-check-inline">
                  <input
                    class="form-check-input"
                    type="radio"
                    name="template-radio-${id}"
                    id="template-radio-code-${id}"
                    value="code"
                  >
                  <label class="form-check-label" for="template-radio-code-${id}">Code</label>
                </div>

         <div id="template-code-section-${id}" class="form-group mb-2" style="display: none;">
              <label class="form-check-label" for="template-code-${id}">PDB code</label>
              <input
                class="form-control"
                name="template-code-${id}"
                id="template-code-${id}"
                type="text"
                style="width: 100px;"
                placeholder="e.g. 1ixc"
                title="Choose an existing PDB code"
                aria-describedby="template-code-desc-${id}"
              >
              <small id="template-code-desc-${id}" class="form-text text-muted">
                Choose an existing PDB code.
              </small>
          </div>
          <div id="template-file-section-${id}" class="form-group mb-2">
                <label class="form-check-label" for="template-file-${id}">PDB file</label>
              <input
                type="file"
                accept=".pdb"
                class="form-control"
                name="template-file-${id}"
                id="template-file-${id}"
                title="Choose a PDB file"
                aria-describedby="template-file-desc-${id}"
              >
              <small id="template-file-desc-${id}" class="form-text text-muted">
                Choose a PDB file.
              </small>
          </div>
                <div class="form-check mb-2">
                        <input class="form-check-input" type="checkbox" id="template-addtemplates-${id}" name="template-addtemplates-${id}" value="true" onchange="updatePlot()" checked>
                        <label class="form-check-label" for="template-addtemplates-${id}"> Add to templates</label>
                </div>
                <div class="form-check mb-2">
                        <input class="form-check-input" type="checkbox" id="template-addmsa-${id}" name="template-addmsa-${id}" value="true" onchange="updatePlot()">
                        <label class="form-check-label" for="template-addmsa-${id}"> Add to MSA</label>
                </div>
                <div class="form-check mb-2">
                        <input class="form-check-input" type="checkbox" id="template-multimer-${id}" name="template-multimer-${id}" value="true"">
                        <label class="form-check-label" for="template-multimer-${id}"> Generate multimer</label>
                </div>
                <div class="form-check mb-2">
                        <input class="form-check-input" type="checkbox" id="template-aligned-${id}" name="template-aligned-${id}" value="true">
                        <label class="form-check-label" for="template-aligned-${id}"> Align template to query sequence</label>
                </div>
                <div class="mb-2">
                  <a class="link-opacity-100" style="display: inline-block; margin-top: 10px;" id="add-modification-${id}" href="javascript:void(0)">Add template modification</a>
                  <div id="modifications-container-${id}" role="region" aria-live="polite"></div>
                </div>
            </fieldset>
        `;
    }
}
customElements.define('template-component', TemplateComponent);


class ModificationComponent extends HTMLElement {
    static formAssociated = true;
    static observedAttributes = ['template-id', 'modification-id', 'value'];

    constructor() {
        super();
        this.attachInternals();
        this.modAminoacidsID = 1;
    }

    connectedCallback() {
        this.templateID = this.getAttribute('template-id');
        this.modificationID = this.getAttribute('modification-id');

        this.render();
        this.attachEventListeners();
    }

    attributeChangedCallback(name, oldValue, newValue) {
        if (oldValue !== newValue) {
            this.render();
            this.attachEventListeners();
        }
    }

    attachEventListeners() {
        this.querySelector(`#template-modify-where-${this.templateID}-${this.modificationID}`).addEventListener('change', () => {
            this.toggleChainSelection();
        });


        this.querySelector(`#add-modaminoacids-${this.templateID}-${this.modificationID}`).addEventListener('click', () => {
            this.addModAminoacid();
        });

        this.querySelector(`.delete-modification`).addEventListener('click', () => {
            this.remove();
            this.triggerUpdatePlot();
        });
    }

    toggleChainSelection() {
        const selection = this.querySelector(`#template-modify-where-${this.templateID}-${this.modificationID}`).value;
        const specificChar = this.querySelectorAll(`div[name=chain-div-${this.templateID}-${this.modificationID}]`);
        specificChar.forEach(div => {
            div.classList.toggle('hidden', selection === 'all');
        });
        this.triggerUpdatePlot();
    }

    toggleModify() {
        const selection = this.querySelector(`#template-modify-residues-${this.templateID}-${this.modificationID}`).checked;
        const divHide = this.querySelector(`#modaminoacids-div-${this.templateID}-${this.modificationID}`);
        divHide.classList.toggle('hidden', !selection);
        this.triggerUpdatePlot();
    }

    addModAminoacid() {
        const aminoacidsContainer = this.querySelector(`#aminoacids-container-${this.templateID}-${this.modificationID}`);
        const modAminoacidComponent = document.createElement('mod-aminoacids-component');

        modAminoacidComponent.setAttribute('template-id', this.templateID);
        modAminoacidComponent.setAttribute('modification-id', this.modificationID);
        modAminoacidComponent.setAttribute('aminoacid-id', this.modAminoacidsID);

        aminoacidsContainer.appendChild(modAminoacidComponent);
        this.triggerUpdatePlot();
        ++this.modAminoacidsID;
    }

    triggerUpdatePlot() {
        updatePlot();
    }

    render() {
        this.innerHTML = `
            <li class="row p-3" id="modify-${this.templateID}-${this.modificationID}">
                <div class="row mb-2">
                    <div class="col-md-auto">
                        <label for="template-modify-where-${this.templateID}-${this.modificationID}">Modify chain</label>
                        <select class="form-select" id="template-modify-where-${this.templateID}-${this.modificationID}" name="template-modify-where-${this.templateID}-${this.modificationID}">
                            <option value="all">All chains</option>
                        </select>
                    </div>
                    <div class="col-md-auto">
                        <label for="template-modify-delete-${this.templateID}-${this.modificationID}"> Delete residues</label>
                        <input type="text" class="form-control" id="template-modify-delete-${this.templateID}-${this.modificationID}" name="template-modify-delete-${this.templateID}-${this.modificationID}" placeholder="e.g. 1, 3-10" title="Select residue numbers to delete in the chain, the rest will be kept" onchange="updatePlot()">
                    </div>
                    <div class="hidden col-md-auto" name="chain-div-${this.templateID}-${this.modificationID}">
                        <label for="template-modify-pos-${this.templateID}-${this.modificationID}">Position</label>
                        <select class="form-select" id="template-modify-pos-${this.templateID}-${this.modificationID}" name="template-modify-pos-${this.templateID}-${this.modificationID}" title="Choose position of the query sequence to insert the chain" onchange="updatePlot()">
                        </select>
                    </div>
                    <div class="col-md-auto delete-mutations">
                        <span class="fa fa-trash-alt delete-icon-format delete-icon delete-modification"></span>
                    </div>
                </div>
                <div class="mb-2">
                    <a class="link-opacity-100" style="display: inline-block; margin-top: 10px;" id="add-modaminoacids-${this.templateID}-${this.modificationID}" href="javascript:void(0)">Replace amino acids</a>
                </div>
                <div class="mb-2" id="modaminoacids-div-${this.templateID}-${this.modificationID}">
                  <div id="aminoacids-container-${this.templateID}-${this.modificationID}" role="region" aria-live="polite"></div>
                </div>

                <hr class="solid" style="margin-bottom: 0px; margin-top: 0px">
            </li>
        `;
    }
}
customElements.define('modification-component', ModificationComponent);



class ModAminoacidsComponent extends HTMLElement {
    static formAssociated = true;
    static observedAttributes = ['template-id', 'modification-id', 'aminoacid-id', 'value'];

    constructor() {
        super();
        this.attachInternals();

    }

    connectedCallback() {
        this.templateID = this.getAttribute('template-id');
        this.modificationID = this.getAttribute('modification-id');
        this.aminoacidID = this.getAttribute('aminoacid-id');
        this.id = `${this.templateID}-${this.modificationID}-${this.aminoacidID}`;
        this.render();
        this.attachEventListeners();
    }

    attachEventListeners() {
        this.querySelector(`#template-modify-amino-select-${this.id}`).addEventListener('change', () => {
            this.selectModify();
        });
        this.querySelector(`#template-modify-amino-pos-${this.id}`).addEventListener('change', () => {
            this.triggerUpdatePlot();
        });

        this.querySelector('.delete-aminoacid').addEventListener('click', () => {
            this.remove();
            this.triggerUpdatePlot();
        });
    }

    selectModify() {
        const selection = this.querySelector(`#template-modify-amino-select-${this.id}`).value;
        const fastaHide = this.querySelector(`#modify-div-fasta-${this.id}`);
        const resnameHide = this.querySelector(`#modify-div-resname-${this.id}`);
        fastaHide.classList.toggle('hidden', selection === 'residue');
        resnameHide.classList.toggle('hidden', selection !== 'residue');
        this.triggerUpdatePlot();
    }

    triggerUpdatePlot() {
        updatePlot();
    }

    render() {
        this.innerHTML = `
            <li class="row mb-2" id="modaminoacids-${this.id}">
                <div class="col-md-auto">
                    <label for="template-modify-amino-pos-${this.id}">Residue numbers</label>
                    <input type="text" class="form-control" id="template-modify-amino-pos-${this.id}" name="template-modify-amino-pos-${this.id}" placeholder="e.g. 1, 3-10" title="Select residue numbers to be replaced with">
                </div>
                <div class="col-md-auto">
                    <label for="template-modify-amino-select-${this.id}">Replace by</label>
                    <select class="form-select" id="template-modify-amino-select-${this.id}" name="template-modify-amino-select-${this.id}">
                        <option selected value="residue">Residue</option>
                        <option value="fasta">Sequence</option>
                    </select>
                </div>
                <div id="modify-div-fasta-${this.id}" class="col-md-auto hidden">
                    <label for="template-modify-amino-fasta-${this.id}">Sequence</label>
                    <textarea class="form-control" id="template-modify-amino-fasta-${this.id}" name="template-modify-amino-fasta-${this.id}" placeholder="e.g. MQQLRFLHYIDA" title="Sequence to replace" rows="1"></textarea>
                </div>
                <div id="modify-div-resname-${this.id}" class="col-md-auto">
                    <label for="template-modify-amino-resname-${this.id}">Type</label>
                    <select class="form-select" id="template-modify-amino-resname-${this.id}" name="template-modify-amino-resname-${this.id}" class="form-control" title="Three letter amino acid to replace">
                        ${aminoacidSelect}
                    </select>
                </div>
                <div class="col-md-auto delete-mutations">
                    <span class="fa fa-trash-alt delete-icon-format delete-icon delete-aminoacid"></span>
                </div>
            </li>`;
    }
}
customElements.define('mod-aminoacids-component', ModAminoacidsComponent);
class sequenceTable extends HTMLElement {
    static formAssociated = true;
    static observedAttributes = ['value'];

    constructor() {
        super();
        this.attachInternals();
        this.sequenceID = this.getAttribute('sequence-id');
        this.mutationID = 1;
    }

    connectedCallback() {
        this.render();
        this.attachEventListeners();
    }

    attachEventListeners() {
        const id =  this.sequenceID;

        this.querySelector(`#add-mutation-${id}`).addEventListener('click', () => {
            this.addMutationLine();
        });
        document.querySelectorAll(`[name="sequence-input-${id}"]`).forEach(radio => {
            radio.addEventListener("change", async () => {
                this.handleRadioChange();
                this.triggerUpdatePlot();
            });
        });
        this.querySelector(`#sequence-fasta-${id}`).addEventListener("change", async () => {
            this.triggerUpdatePlot();
        });
        this.querySelector(`#sequence-text-${id}`).addEventListener("change", async () => {
            this.triggerUpdatePlot();
        });
    }

    triggerUpdatePlot() {
        updatePlot();
    }


      handleRadioChange() {
        const id = this.sequenceID;
        const radioFile = this.querySelector(`input[type="radio"][value="file"]`);
        const fastaSection = this.querySelector(`#sequence-fasta-section-${id}`);
        const textSection = this.querySelector(`#sequence-text-section-${id}`);
        if (radioFile?.checked) {
          fastaSection.style.display = '';
          textSection.style.display = 'none';
        } else {
          fastaSection.style.display = 'none';
          textSection.style.display = '';
        }
      }

    render() {
        const id = this.sequenceID;
        this.innerHTML =  `            
            <fieldset name="sequence-field">
                <div class="form-group mb-2">
                  <label class="me-3">Read sequence from</label>
                <div class="form-check form-check-inline">
                  <input
                    class="form-check-input"
                    type="radio"
                    name="sequence-input-${id}"
                    id="sequence-radio-fasta-${id}"
                    value="file"
                    checked
                  >
                  <label class="form-check-label" for="sequence-radio-fasta-${id}">File</label>
                </div>
                <div class="form-check form-check-inline">
                  <input
                    class="form-check-input"
                    type="radio"
                    name="sequence-input-${id}"
                    id="sequence-radio-text-${id}"
                    value="text"
                  >
                  <label class="form-check-label" for="sequence-radio-text-${id}">Text</label>
                </div>

         <div id="sequence-fasta-section-${id}" class="form-group mb-2">
              <label class="form-check-label" for="sequence-fasta-${id}">Sequence file</label>
              <input
                class="form-control"
                name="sequence-fasta-${id}"
                id="sequence-fasta-${id}"
                type="file"
                title=" Choose a sequence file"
                accept=".seq, .fasta, .sequence"
                aria-describedby="sequence-fasta-desc-${id}"
              >
              <small id="sequence-fasta-desc-${id}" class="form-text text-muted">
                Choose a sequence file.
              </small>
          </div>
          <div id="sequence-text-section-${id}" class="form-group mb-2" style="display: none;">
                <label class="form-check-label" for="sequence-text-${id}">Sequence</label>
              <input
                type="text"
                class="form-control"
                name="sequence-text-${id}"
                id="sequence-text-${id}"
                title="Insert a sequence"
                aria-describedby="sequence-text-desc-${id}"
                placeholder="e.g. MQHLRFLHYIDAVARCGSIRAA"
              >
              <small id="sequence-text-desc-${id}" class="form-text text-muted">
                Insert a sequence.
              </small>
          </div>
                <div class="row mb-2">
                    <div class="col-md-auto">
                        <label for="sequence-copies-${id}">Number of copies</label>
                        <input type="text" class="form-control" id="sequence-copies-${id}" style="width:100px;" name="sequence-copies-${id}" min="1" data-bind="value:replyNumber" value="1" maxlength="1" onchange="updatePlot()" title="Number of copies of the query sequence" required>
                        <div class="invalid-feedback">
                            Mandatory field
                        </div>
                    </div>
                </div>
                <div class="mb-2">
                  <a class="link-opacity-100" style="display: inline-block; margin-top: 10px;" id="add-mutation-${id}" name="add-mutation-${id}" href="javascript:void(0)">Add sequence mutation</a>
                  <div id="mutations-container-${id}" role="region" aria-live="polite"></div>
                </div>
            </fieldset>
        `;
    }

    addMutationLine() {
        const addMutationButton = this.querySelector(`#mutations-container-${this.sequenceID}`);
        const mutationLine = document.createElement('li');
        mutationLine.classList.add('row', 'p-2');
        mutationLine.id = `li-mutation-${this.sequenceID}-${this.mutationID}`;
        mutationLine.innerHTML = `
            <div class="col-md-auto">
                <label for="sequence-mutations-res-${this.sequenceID}-${this.mutationID}"> Type</label>
                <select class="form-select" name="sequence-mutations-res-${this.sequenceID}-${this.mutationID}" id="sequence-mutations-res-${this.sequenceID}-${this.mutationID}" class="form-control" title="Aminoacid (e.g. G)">
                    ${aminoacidSelect}
                </select>
            </div>
            <div class="col-md-auto">
                <label for="sequence-mutations-pos-${this.sequenceID}-${this.mutationID}"> Residue number</label>
                <input type="text" name="sequence-mutations-pos-${this.sequenceID}-${this.mutationID}" id="sequence-mutations-pos-${this.sequenceID}-${this.mutationID}" class="form-control" placeholder="E.g. 1, 2-100" title="Residues to replace (e.g. 1, 2-100)" onchange="updatePlot()" required>
            </div>
            <div class="col-md-auto delete-mutations">
                <span onclick="this.parentNode.parentNode.remove(); updatePlot()" class="fa fa-trash-alt delete-icon-format delete-icon" ></span>
            </div>
        `;
        addMutationButton.appendChild(mutationLine);
        ++this.mutationID;
    }
}

customElements.define('sequence-component', sequenceTable);

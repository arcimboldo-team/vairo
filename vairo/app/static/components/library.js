class LibraryTable extends HTMLElement {
  static formAssociated = true;
  static observedAttributes = ['value'];

  constructor() {
    super();
    this.attachInternals();
    this.libraryID = this.getAttribute('library-id');

    this.boundHandlers = {
      handleRadioChange: this.handleRadioChange.bind(this),
      handleFolderChange: this.handleFolderChange.bind(this),
      handleFastaChange: this.handleFastaChange.bind(this),
      handleAddMsaChange: this.handleAddMsaChange.bind(this),
      handleAddTemplatesChange: this.handleAddTemplatesChange.bind(this),
      handleQueryChange: this.handleQueryChange.bind(this),
    };

    this.handlers = [
      ['radioFolder', 'change', 'handleRadioChange'],
      ['radioFasta', 'change', 'handleRadioChange'],
      ['folderInput', 'change', 'handleFolderChange'],
      ['fastaInput', 'change', 'handleFastaChange'],
      ['addMsa', 'change', 'handleAddMsaChange'],
      ['addTemplates', 'change', 'handleAddTemplatesChange'],
      ['queryPositions', 'change', 'handleQueryChange'],
      ['libraryPositions', 'change', 'handleQueryChange'],
    ];
  }

  connectedCallback() {
    this.render();
    this.cacheElements();
    this.attachEventListeners();
  }

  disconnectedCallback() {
    this.detachEventListeners();
  }

  cacheElements() {
    const id = this.libraryID;
    this._elements = {
      radioFolder: this.querySelector(`input[type="radio"][value="folder"]`),
      radioFasta: this.querySelector(`input[type="radio"][value="fasta"]`),
      folderInput: this.querySelector(`#library-folder-${id}`),
      fastaInput: this.querySelector(`#library-fasta-${id}`),
      addMsa: this.querySelector(`#library-addmsa-${id}`),
      addTemplates: this.querySelector(`#library-addtemplates-${id}`),
      templatesDiv: this.querySelector(`#library-addtemplates-div-${id}`),
      libPositions: this.querySelector(`#library-lib-${id}`),
      queryPositions: this.querySelector(`#library-query-${id}`),
      libraryPositions: this.querySelector(`#library-lib-${id}`),
      folderSection: this.querySelector(`#library-folder-section-${id}`),
      fastaSection: this.querySelector(`#library-fasta-section-${id}`),
    };
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

  disconnectedCallback() {
    if (!this._elements) return;
    this.handlers.forEach(([key, evt, handler]) => {
      const el = this._elements[key];
      if (el && typeof this.boundHandlers[handler] === 'function') {
        el.removeEventListener(evt, this.boundHandlers[handler]);
      }
    });
  }

  async readLibrary() {
    const id = this.libraryID;
    const selectedRadio = this.querySelector(`input[type="radio"][name="library-input-${id}"]:checked`);
    let counterMsa = 0;
    let counterTemplates = 0;
    const addMsa = this._elements.addMsa.checked;

    if (selectedRadio.value === 'folder') {
      const folder = this._elements.folderInput;
      const addTemplates = this._elements.addTemplates.checked;
      const files = folder.files;
      if (files) {
        for (let i = 0; i < files.length; i++) {
          const file = files[i];
          if (file.name.endsWith('.pdb')) {
            if (addMsa) counterMsa++;
            if (addTemplates) counterTemplates++;
          } else if ((file.name.endsWith('.fasta') || file.name.endsWith('.seq')) && addMsa) {
            counterMsa++;
          }
        }
      }
    } else {
      const fileInput = this._elements.fastaInput?.files?.[0];
      if (fileInput) {
        const fileData = await readFile(fileInput);
        const lines = fileData.split('\n');
        if (addMsa) {
          for (let i = 0; i < lines.length; i++) {
            if (lines[i][0] !== '>') {
              counterMsa++;
            }
          }
        }
      }
    }
    librariesDict[id] = { "msa": counterMsa, "templates": counterTemplates };
    this.triggerUpdatePlot();
  }

  handleRadioChange() {
    const { radioFolder, folderSection, fastaSection } = this._elements;
    if (radioFolder?.checked) {
      folderSection.style.display = '';
      fastaSection.style.display = 'none';
    } else {
      folderSection.style.display = 'none';
      fastaSection.style.display = '';
    }
    this.insertModify();
  }

  insertModify() {
    const id = this.libraryID;
    const selectedRadio = this.querySelector(`input[type="radio"][name="library-input-${id}"]:checked`);
    this._elements.templatesDiv.classList.toggle('hidden', selectedRadio.value === 'fasta');

    if (!this._elements.addMsa.checked) {
      this._elements.addMsa.checked = true;
    }
    this.readLibrary();
  }

  handleFolderChange() {
    this.readLibrary();
  }

  handleFastaChange() {
    this.readLibrary();
  }

  handleAddMsaChange() {
    this.readLibrary();
  }

  handleAddTemplatesChange() {
    this.readLibrary();
  }

  handleQueryChange() {
    this.triggerUpdatePlot();
  }

  triggerUpdatePlot() {
    updatePlot()
  }

  render() {
    const id = this.libraryID;
    this.innerHTML = `
      <fieldset name="library-field">
        <div class="form-group mb-2">
          <label class="me-3">Insert library from</label>
        <div class="form-check form-check-inline">
          <input
            class="form-check-input"
            type="radio"
            name="library-input-${id}"
            id="library-radio-folder-${id}"
            value="folder"
            checked
          >
          <label class="form-check-label" for="library-radio-folder-${id}">Folder</label>
        </div>
        <div class="form-check form-check-inline">
          <input
            class="form-check-input"
            type="radio"
            name="library-input-${id}"
            id="library-radio-fasta-${id}"
            value="fasta"
          >
          <label class="form-check-label" for="library-radio-fasta-${id}">FASTA</label>
        </div>
          <div id="library-folder-section-${id}" class="form-group mb-2">
              <input
                class="form-control"
                name="library-folder-${id}"
                id="library-folder-${id}"
                type="file"
                webkitdirectory directory
                multiple
                title="Choose a folder containing FASTA files or PDBs"
                aria-describedby="library-folder-desc-${id}"
              >
              <small id="library-folder-desc-${id}" class="form-text text-muted">
                Choose a folder containing FASTA files or PDBs.
              </small>
          </div>
          <div id="library-fasta-section-${id}" class="form-group mb-2" style="display: none;">
              <input
                type="file"
                accept=".fasta"
                class="form-control"
                name="library-fasta-${id}"
                id="library-fasta-${id}"
                title="Choose a FASTA file containing sequences"
                aria-describedby="library-fasta-desc-${id}"
              >
              <small id="library-fasta-desc-${id}" class="form-text text-muted">
                Choose a FASTA file containing sequences.
              </small>
          </div>

        <div class="form-check mb-2">
            <input class="form-check-input" type="checkbox" id="library-addmsa-${id}" name="library-addmsa-${id}" value="true">
            <label class="form-check-label" for="library-addmsa-${id}"> Add to MSA</label>
        </div>
        <div class="form-check mb-2" id="library-addtemplates-div-${id}">
            <input class="form-check-input" type="checkbox" id="library-addtemplates-${id}" name="library-addtemplates-${id}" value="true" checked>
            <label class="form-check-label" for="library-addtemplates-${id}"> Add to templates</label>
        </div>

        <div class="form-group mb-2 row">
          <div class="col-md-auto">
            <label for="library-query-${id}">Query sequence positions</label>
            <input type="text" class="form-control" id="library-query-${id}" name="library-query-${id}" placeholder="e.g. 1, 30" title="Select the positions to insert it in the query sequence">
          </div>
          <div class="col-md-auto">
            <label for="library-lib-${id}">Library positions</label>
            <input type="text" class="form-control" id="library-lib-${id}" name="library-lib-${id}" placeholder="e.g. 1-10, 40-50" title="Select the residue numbers from the library">
          </div>
        </div>
      </fieldset>
    `;
  }
}

customElements.define('library-component', LibraryTable);

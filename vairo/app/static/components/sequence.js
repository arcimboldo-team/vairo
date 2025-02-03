class sequenceTable extends HTMLElement {
    static formAssociated = true;
    static observedAttributes = ['value'];

    constructor() {
        super();
        this.attachInternals();
        this.sequenceID = this.getAttribute('sequenceID');
        this.mutationID = 1;
    }

    connectedCallback() {
        this.render();
        const addButton = this.querySelector(`#add-mutation-${this.sequenceID}`);
        addButton.addEventListener('click', this.addMutationLine.bind(this));
    }

    render() {
        this.innerHTML =  `            
            <fieldset name="sequence-field" class="row g-3"> 
                <div class="row row-margin">
                    <div class="form-group">
                        <label class="form-label" for="sequence-fasta-${this.sequenceID}">Insert sequence</label>
                        <div class="form-check radio-container">
                            <input class="form-check-input" type="radio" name="sequence-input-${this.sequenceID}" value="file" onchange="updatePlot()" checked>
                            <input class="form-control" type="file" name="sequence-fasta-${this.sequenceID}" accept=".fasta, .seq, .sequence" id="sequence-fasta-${this.sequenceID}" title="Choose sequence file (fasta format)" onchange="updatePlot()">
                        </div>
                        <div class="form-check radio-container">
                            <input class="form-check-input" type="radio" name="sequence-input-${this.sequenceID}" value="text" onchange="updatePlot()">
                            <input class="form-control" type="text" name="sequence-text-${this.sequenceID}" id="sequence-text-${this.sequenceID}"  title="Write amino acid sequence" placeholder="Sequence: MQHLRFLHYIDAVARCGSIRAA" onchange="updatePlot()">
                        </div>
                    </div>
                </div>
                <div class="row row-margin">
                    <div class="col-md-auto">
                        <label class="form-label" for="sequence-copies-${this.sequenceID}">Number of copies</label>
                        <input type="text" class="form-control" id="sequence-copies-${this.sequenceID}" name="sequence-copies-${this.sequenceID}" min="1" data-bind="value:replyNumber" value="1" maxlength="1" onchange="updatePlot()" title="Number of copies of the query sequence" required>
                        <div class="invalid-feedback">
                            Mandatory field
                        </div>
                    </div>
                    <div class="col-md-auto">
                        <label class="form-label" for="sequence-positions-${this.sequenceID}">Positions</label>
                        <input type="text" class="form-control" name="sequence-positions-${this.sequenceID}" id="sequence-positions-${this.sequenceID}" name="sequence-positions-${this.sequenceID}" onchange="updatePlot()" placeholder="1, 3, 5" title="(Optional) Specify any position inside the query sequence.">
                    </div>
                </div>
                <div class="row row-margin">
                    <label>Mutations:</label>
                    <div class="col-md-auto">
                        <ul id="ul-mutation-${this.sequenceID}"></ul>
                        <a class="link-opacity-100 link-line-padding" name="add-mutation-${this.sequenceID}" id="add-mutation-${this.sequenceID}" href="javascript:void(0)">Add mutation</a>
                    </div>
                </div>
            </fieldset>
        `;
    }

    addMutationLine() {
        const addMutationButton = this.querySelector(`#ul-mutation-${this.sequenceID}`);
        const mutationLine = document.createElement('li');
        mutationLine.classList.add('row', 'g-3', 'element-line-padding');
        mutationLine.id = `li-mutation-${this.sequenceID}-${this.mutationID}`;
        mutationLine.innerHTML = `
            <div class="col-md-auto">
                <label class="form-label" for="sequence-mutations-res-${this.sequenceID}-${this.mutationID}"> Type</label>
                <select class="form-select" name="sequence-mutations-res-${this.sequenceID}-${this.mutationID}" id="sequence-mutations-res-${this.sequenceID}-${this.mutationID}" class="form-control" title="Aminoacid (e.g. G)">
                    ${aminoacidSelect}
                </select>
            </div>
            <div class="col-md-auto">
                <label class="form-label" for="sequence-mutations-pos-${this.sequenceID}-${this.mutationID}"> Residue number</label>
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

const librariesDict = {};


function deleteLibrary(id){
    delete librariesDict[id];
    updatePlot();
}

async function readLibrary(id){
    const selectedRadio = document.querySelector(`input[type=radio][name="library-input-${id}"]:checked`);
    let counterMsa = 0;
    let counterTemplates = 0;
    if(selectedRadio.value === 'folder'){
        const folder = document.getElementById(`library-folder-${id}`);
        const addMsa = document.getElementById(`library-addmsa-${id}`).checked;
        const addTemplates = document.getElementById(`library-addtemplates-${id}`).checked;
        var files = folder.files;
        for (let i = 0; i < files.length; i++) {
            const file = files[i];            
            if (file.name.endsWith('.pdb')) {
                if(addMsa)
                    counterMsa += 1;
                if(addTemplates)
                    counterTemplates += 1;
            }
            else if ((file.name.endsWith('.fasta') ||  file.name.endsWith('.seq')) && addMsa){
                counterMsa += 1;
            }
        }
    } else {
        const fileInput = document.getElementById(`library-folder-${id}`)?.files?.[0];
        const fileData = await readFile(fileInput);
        const lines = fileData.split('\n');
        if(addMsa){
            for (let i = 0; i < lines.length; i++) {
                if (lines[i][0] !== '>') {
                    counterMsa++;
                }
            }
        }
    }
    librariesDict[id] = {"msa": counterMsa, "templates": counterTemplates};
    updatePlot();
}

function insertModify(id){
    const selectedRadio = document.querySelector(`input[type=radio][name="library-input-${id}"]:checked`);
    const divHide = document.getElementById(`library-addtemplates-div-${id}`);
    divHide.classList.toggle('hidden', selectedRadio.value === 'fasta');
    
    const addMsa = document.getElementById(`library-addmsa-${id}`);
    if(!addMsa.checked){
        addMsa.checked = true;
    }
    readLibrary(id);
}


class libraryTable extends HTMLElement {
    static formAssociated = true;
    static observedAttributes = ['value'];

    constructor() {
        super();
        this.attachInternals();
        this.libraryID = this.getAttribute('libraryID');
    }

    connectedCallback() {
        this.render();
    }

    render() {
        this.innerHTML =  `            
            <fieldset name="library-field" class="row g-3"> 
                <div class="form-group">
                    <label class="form-label" for="library-folder-${this.libraryID}">Insert library from</label>
                    <div class="form-check radio-container">
                        <div class="col-md-auto">
                            <input class="form-check-input" type="radio" name="library-input-${this.libraryID}" value="folder" onchange="insertModify('${this.libraryID}')" checked>
                            <label class="form-label" for="library-folder-${this.libraryID}" style="margin-right: 10px;">Folder with PDB or FASTA Files</label>
                        </div>
                        <div class="col-md-9">
                            <input class="form-control" name="library-folder-${this.libraryID}" id="library-folder-${this.libraryID}" type="file" webkitdirectory directory multiple title="Choose a folder containing PDB or FASTA files" onchange="readLibrary('${this.libraryID}')">
                        </div>
                    </div>
                    <div class="form-check radio-container">
                        <div class="col-md-auto">
                            <input class="form-check-input" type="radio" name="library-input-${this.libraryID}" value="fasta" onchange="insertModify('${this.libraryID}')">
                            <label class="form-label" for="library-fasta-${this.libraryID}" style="margin-right: 10px;">FASTA file with sequences</label>
                        </div>
                        <div class="col-md-9">
                            <input type="file" accept=".fasta" class="form-control" name="library-fasta-${this.libraryID}" id="library-fasta-${this.libraryID}" title="Choose fasta file with sequences" onchange="readLibrary('${this.libraryID}')">
                        </div>  
                    </div>
                </div>    
                <div class="row row-margin" id="library-add-templates-div-${this.libraryID}">
                    <div class="col-md-auto">
                        <input type="checkbox" id="library-addtemplates-${this.libraryID}" name="library-addtemplates-${this.libraryID}" value="true" onchange="readLibrary('${this.libraryID}')" checked>
                        <label class="form-label" for="library-addtemplates-${this.libraryID}"> Add to templates</label>
                    </div>
                </div>
                <div class="row">
                    <div class="col-md-auto">
                        <input type="checkbox" id="library-addmsa-${this.libraryID}" name="library-addmsa-${this.libraryID}" value="true" onchange="readLibrary('${this.libraryID}')">
                        <label class="form-label" for="library-addmsa-${this.libraryID}"> Add to MSA</label>
                    </div>
                </div>
                <div class="row row-margin">
                    <div class="col-md-auto">
                        <label class="form-label" for="library-lib-${this.libraryID}">Library positions</label>
                        <input type="text" class="form-control" id="library-lib-${this.libraryID}" name="library-lib-${this.libraryID}" placeholder="1, 3-10" title="Select the residues numbers from the library">
                    </div>
                    <div class="col-md-auto">
                        <label class="form-label" for="library-query-${this.libraryID}">Query positions</label>
                        <input type="text" class="form-control" id="library-query-${this.libraryID}" name="library-query-${this.libraryID}" onchange="updatePlot()" placeholder="1, 3-10" title="Select the positions to insert it in the query sequence">
                    </div>
                </div>
            </fieldset>
        `;
    }
}

customElements.define('library-component', libraryTable);
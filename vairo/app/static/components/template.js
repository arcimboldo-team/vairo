const templatesDict = {};


function toggleModify(templateID, modificationID) {
    const selection = document.getElementById(`template-modify-residues-${templateID}-${modificationID}`).checked;
    const divHide = document.getElementById(`modaminoacids-div-${templateID}-${modificationID}`);
    divHide.classList.toggle('hidden', !selection);
    updatePlot();
}

function toggleChainSelection(templateID, modificationID) {
    const selection = document.getElementById(`template-modify-where-${templateID}-${modificationID}`).value;
    const specificChar = document.querySelectorAll(`div[name=chain-div-${templateID}-${modificationID}]`);
    specificChar.forEach(div => {
        div.classList.toggle('hidden', selection === 'all');
    });
    updatePlot();
}

function selectModify(id) {
    const selection = document.getElementById(`template-modify-amino-select-${id}`).value;
    const fastaHide = document.getElementById(`modify-div-fasta-${id}`);
    const resnameHide = document.getElementById(`modify-div-resname-${id}`);
    fastaHide.classList.toggle('hidden', selection === 'residue');
    resnameHide.classList.toggle('hidden', selection !== 'residue');
    updatePlot();
}

async function generateMultimer(key, templateName="undefined", templateData="undefined") {
    if(!templatesDict.hasOwnProperty(key) && templateData === "undefined") return;
    if(templateData === "undefined"){
        templateData = templatesDict[key].templateData;
        templateName = templatesDict[key].templateName;
    }
    let chainSeq = getSequences(templateData);
    for (let key in chainSeq) {
        chainSeq[key] = [chainSeq[key]];
    }
    check = document.getElementById(`template-multimer-${key}`);
    if(check.checked){
        try{
            chainSeq = await postData('/generate-multimer', {'templateData': templateData});
        } catch (error) {
            alert('It has not been possible to generate the multimer. Unchecking it');
            check.checked = false;
            console.error('Error:', error);
        }
    }
    templatesDict[key] = {"templateName": templateName, "templateData": templateData, "chainSequences": chainSeq};
    updatePlot();
    populateChainSelect(key);
}

function deleteTemplate(id){
    delete templatesDict[id];
    updatePlot();
}

async function readTemplate(id){
    var templateData;
    var templateName;
    const selected = document.querySelector(`input[name="template-input-${id}"]:checked`);
    if (selected.value === "file") {
        const pdbFileInput = document.getElementById(`template-file-${id}`)?.files?.[0];
        if (pdbFileInput === null || pdbFileInput === undefined) return;
        templateName  = pdbFileInput.name.split('.')[0];
        if (id in templatesDict && templatesDict[id].templateName === templateName) return;
        templateData = await readFile(pdbFileInput);
    } else {
        templateName = document.getElementById(`template-code-${id}`).value;
        if ((id in templatesDict && templatesDict[id].templateName === templateName) || (templateName === "")) return;
        const fetchPromise = fetchPDB(templateName);
        const timeoutPromise = new Promise((resolve, reject) => {
            setTimeout(() => reject(new Error('Timeout')), 60000);
        });
        try {
            templateData = await Promise.race([fetchPromise, timeoutPromise]);
        } catch (error) {
            console.error('Error fetching PDB:', error);
            alert('The request has timed out. Please try again later.');
            document.getElementById(`template-code-${id}`).value = "";
            return;
        }
    }
    generateMultimer(id, templateName, templateData);
}


function getRestrictions(templateID){
    const divs = document.querySelectorAll(`li[id^=modify-${templateID}-]`);
    const modificationsDict = {};
    let positionsAccepted = false;
    divs.forEach(div => {
        const selection = div.querySelector(`select[id^=template-modify-where-${templateID}-]`).value;
        if(selection === "") return;
        modificationsDict[selection] = modificationsDict[selection] || {};
        let positionKey = ''
        if(selection !== "all"){
            positionKey = div.querySelector(`select[id^=template-modify-pos-${templateID}-]`).value;
            if(positionKey !== 'ANY'){
                positionsAccepted = true;
            }
        }
        else{
            positionKey = 'ANY'
        }
        modificationsDict[selection][positionKey] = modificationsDict[selection][positionKey] || {};

        const residues = div.querySelector(`input[id^=template-modify-delete-${templateID}-]`).value.trim();
        if(residues){
            if (!modificationsDict[selection][positionKey].hasOwnProperty('delete')) {
                modificationsDict[selection][positionKey]["delete"] = new Set();
            }
            const oldValues = modificationsDict[selection][positionKey]["delete"];
            const newValues = extendedNumbers(residues).sort(function(a, b) { return a - b });
            console.log(oldValues);
            console.log(newValues);
            modificationsDict[selection][positionKey]["delete"] = new Set(...oldValues, [...newValues]);
        }
        const changeAminoacids = div.querySelectorAll(`li[id^=modaminoacids-${templateID}-]`);
        changeAminoacids.forEach(change => {
            const positions = change.querySelector(`input[id^=template-modify-amino-pos-${templateID}]`).value.trim();
            const sequence = change.querySelector(`input[id^=template-modify-amino-fasta-${templateID}-]`);
            const choose = change.querySelector(`select[id^=template-modify-amino-select-${templateID}-]`).value;
            if(positions && ((choose === "residue") || (choose === "fasta" && sequence.files.length > 0))){
                const positionsArray = extendedNumbers(positions);
                modificationsDict[selection][positionKey][choose] = [...new Set([...(modificationsDict[selection][positionKey][choose] || []), ...positionsArray])];
            }
        });
    });
    return [modificationsDict, positionsAccepted];
}

function populateChainSelect(templateID) {
    let options = '';
    if(templatesDict.hasOwnProperty(templateID)){
        const chains = templatesDict[templateID].chainSequences;
        const uniqueChains = Object.keys(chains)
        options = `<option value="all">All chains</option>` + uniqueChains.map(chain => `<option value="${chain}">${chain} (${chains[chain].length} copies)</option>`).join('');
    } else {
        options = `<option value="all">All chains</option>`;
    }
    const modifications = document.querySelectorAll(`select[id^=template-modify-where-${templateID}-]`);
    modifications.forEach(modification => {
        const oldValue = modification.value;
        modification.innerHTML = options;
        const optionExists = Array.from(modification.options).some(option => option.value === oldValue);
        if(optionExists){
            modification.value = oldValue;
        }
    });
}

class templateTable extends HTMLElement {
    static formAssociated = true;
    static observedAttributes = ['value'];

    constructor() {
        super();
        this.attachInternals();
        this.templateID = this.getAttribute('templateID');
        this.modificationID = 1;
        this.modAminoacidsID = 1;
    }

    connectedCallback() {
        this.render();
        this.querySelector(`#add-modification-${this.templateID}`).addEventListener('click', () => this.addModificationLine());
        document.getElementById(`template-code-${this.templateID}`).addEventListener("change", async () => {
            await readTemplate(this.templateID); 
        });
        document.getElementById(`template-file-${this.templateID}`).addEventListener("change", async () => {
            await readTemplate(this.templateID); 
        });      
        const radioInputTemplate = this.querySelectorAll(`input[name=template-input-${this.templateID}]`);
        radioInputTemplate.forEach(radio => {
            radio.addEventListener("change", async () => {
                await readTemplate(this.templateID);
            });
        });
    }

    render() {
        this.innerHTML =  `
            <fieldset name="template-field" class="row g-3"> 
                <div class="row row-margin">
                    <div class="form-group">
                        <label class="form-label" for="template-code-${this.templateID}">Insert template PDB</label>
                        <div class="form-check radio-container">
                            <input class="form-check-input" type="radio" name="template-input-${this.templateID}" value="code" checked>
                            <input type="text" class="form-control" name="template-code-${this.templateID}" id="template-code-${this.templateID}" placeholder="Insert PDB code (e.g 1ixc)" title="Insert PDB code (e.g 1ixc)">
                        </div>
                        <div class="form-check radio-container">
                            <input class="form-check-input" type="radio" name="template-input-${this.templateID}" value="file">
                            <input type="file" name="template-file-${this.templateID}" accept=".pdb" class="form-control" id="template-file-${this.templateID}" title="Choose pdb file (pdb format)">
                        </div>
                    </div>
                </div>
                <div class="row row-margin">
                    <div class="col-md-auto">
                        <input type="checkbox" id="template-addtemplates-${this.templateID}" name="template-addtemplates-${this.templateID}" value="true" onchange="updatePlot()" checked>
                        <label class="form-label" for="template-addtemplates-${this.templateID}"> Add to templates</label>
                    </div>
                </div>
                <div class="row">
                    <div class="col-md-auto">
                        <input type="checkbox" id="template-addmsa-${this.templateID}" name="template-addmsa-${this.templateID}" value="true" onchange="updatePlot()">
                        <label class="form-label" for="template-addmsa-${this.templateID}"> Add to MSA</label>
                    </div>
                </div>
                <div class="row">
                    <div class="col-md-auto">
                        <input type="checkbox" id="template-multimer-${this.templateID}" name="template-multimer-${this.templateID}" value="true" onchange="generateMultimer(${this.templateID})">
                        <label class="form-label" for="template-multimer-${this.templateID}"> Generate multimer</label>
                    </div>
                </div>
                <div class="row">
                    <div class="col-md-auto">
                        <input type="checkbox" id="template-aligned-${this.templateID}" name="template-aligned-${this.templateID}" value="true">
                        <label class="form-label" for="template-aligned-${this.templateID}"> Align template to query sequence</label>
                    </div>
                </div>
                <div class="row row-margin-after-check">
                    <label>Modify template:</label>
                    <div id="modification-res-${this.templateID}">
                        <ul id="ul-modification-${this.templateID}"></ul>
                        <a class="link-opacity-100 link-line-padding" id="add-modification-${this.templateID}" href="javascript:void(0)">Add modification</a>
                    </div>
                </div>
            </fieldset>
        `;
    }

    addModAminoacidsLine(modificationID, firstParamater = false) {
        const addModAminoacidsButton = this.querySelector(`#ul-modaminoacids-${this.templateID}-${modificationID}`);
        const modAminoacidsLine = document.createElement('li');
        const id = `${this.templateID}-${modificationID}-${this.modAminoacidsID}`;
        modAminoacidsLine.classList.add('row', 'g-3', 'element-line-padding');
        modAminoacidsLine.id = `modaminoacids-${id}`;
        modAminoacidsLine.innerHTML = `
                <div class="col-md-auto">
                    <label class="form-label" for="template-modify-amino-pos-${id}">Residue numbers</label>
                    <input type="text" class="form-control" id="template-modify-amino-pos-${id}" name="template-modify-amino-pos-${id}" placeholder="1, 3-10" title="Specify any position inside the template" onchange="updatePlot()">
                </div>
                <div class="col-md-auto">
                    <label class="form-label" for="template-modify-amino-select-${id}">Replace by</label>
                    <select class="form-select" id="template-modify-amino-select-${id}" name="template-modify-amino-select-${id}" onchange="selectModify('${id}')">
                        <option selected value="residue">Residue</option>
                        <option value="fasta">Sequence</option>
                    </select>
                </div>
                <div id="modify-div-fasta-${id}" class="col-md-auto hidden">
                    <label class="form-label" for="template-modify-amino-fasta-${id}">Insert sequence</label>
                    <input class="form-control" type="file" accept=".fasta, .seq, .sequence" name="template-modify-amino-fasta-${id}" id="template-modify-amino-fasta-${id}" title="Choose sequence file (fasta format)" onchange="updatePlot()">
                </div>
                <div id="modify-div-resname-${id}" class="col-md-auto">
                    <label class="form-label" for="template-modify-amino-resname-${id}">Type</label>
                    <select class="form-select" id="template-modify-amino-resname-${id}" name="template-modify-amino-resname-${id}" class="form-control" title="Three letter aminoa acid to replace">
                        ${aminoacidSelect}
                    </select>
                </div>`;
        if (!firstParamater) {
            modAminoacidsLine.innerHTML += `
                <div class="col-md-auto delete-mutations">
                    <span onclick="this.parentNode.parentNode.remove(); updatePlot()" class="fa fa-trash-alt delete-icon-format delete-icon" ></span>
                </div>`;
        }
        addModAminoacidsButton.appendChild(modAminoacidsLine);
        ++this.modAminoacidsID;
    }

    addModificationLine() {
        const optionsHTML = selectPositionsArray.map(option => `<option value="${option}">${option}</option>`).join('');
        const addModificationButton = this.querySelector(`#ul-modification-${this.templateID}`);
        const modificationLine = document.createElement('li');
        const id = `${this.templateID}-${this.modificationID}`;
        const currentModificationValue = this.modificationID;
        modificationLine.classList.add('row', 'g-3', 'element-line-padding');
        modificationLine.id = `modify-${id}`;
        modificationLine.innerHTML = `
            <div class="row row-margin" style="margin-top: 16px">
                <div class="col-md-auto">
                    <label class="form-label" for="template-modify-where-${id}">Apply modifications</label>
                    <select class="form-select" id="template-modify-where-${id}" name="template-modify-where-${id}" onchange="toggleChainSelection(${this.templateID}, ${currentModificationValue})">
                        <option value="all">All chains</option>
                    </select>
                </div>
                <div class="col-md-auto">
                    <label class="form-label" for="template-modify-delete-${id}"> Delete residues</label>
                    <input type="text" class="form-control" id="template-modify-delete-${id}" name="template-modify-delete-${id}" placeholder="Default: ALL. 1, 3-10" title="Select residue numbers to delete in the chain, the rest will be deleted" onchange="updatePlot()">
                </div>
                <div class="hidden col-md-auto" name="chain-div-${id}">
                    <label class="form-label" for="template-modify-pos-${id}">Position</label>
                    <select class="form-select" id="template-modify-pos-${id}" name="template-modify-pos-${id}" title="Choose position of the query sequence to insert the chain" onchange="updatePlot()">
                        ${optionsHTML}
                    </select>
                </div>
            </div>
            <div class="row row-margin">
                <div class="col-md-auto">
                    <input type="checkbox" name="template-modify-residues-${id}" id="template-modify-residues-${id}" onchange="toggleModify(${this.templateID}, ${currentModificationValue})" value="true">
                    <label class="form-label" for="template-modify-residues-${id}"> Modify amino acids</label>
                </div>
            </div>
            <div class="aminoacid-line-padding hidden" style="margin-top: 0px" id="modaminoacids-div-${id}">
                <ul id="ul-modaminoacids-${id}" style="margin-bottom: 10px"></ul>
                <a class="link-opacity-100 link-line-padding" id="add-modaminoacids-${id}" href="javascript:void(0)">Add amino acid change</a>
            </div>
            <div class="row">
                <div class="col-md-1 offset-md-6 delete-mutations">
                    <span onclick="this.parentNode.parentNode.parentNode.remove(); updatePlot()" class="fa fa-trash-alt delete-icon-format delete-icon" ></span>
                </div>
            </div>
            <hr class="solid" style="margin-bottom: 0px; margin-top: 0px">

        `;
        addModificationButton.appendChild(modificationLine);
        this.querySelector(`#add-modaminoacids-${this.templateID}-${this.modificationID}`).addEventListener('click', () => this.addModAminoacidsLine(currentModificationValue));
        this.addModAminoacidsLine(this.modificationID, true)
        ++this.modificationID;
        populateChainSelect(this.templateID, true);
    }
}
customElements.define('template-component', templateTable);

function createLi(text) {
    const li = document.createElement('li');
    li.textContent = text;
    return li;
}

const featuresDict = {};
const librariesDict = {};
const templatesDict = {};

function deleteInDict(type, id) {
    switch(type) {
        case 'template':
            delete templatesDict[id];
            break;
        case 'library':
            delete librariesDict[id];
            break;
        case 'feature':
            delete featuresDict[id];
            break;
    }
    updatePlot();
}

class Summary {
    constructor() {
        this.wrapper = "summary";

        $(`#${this.wrapper}`).append(`
            <div class="card-header">
                <h2 class="card-title">Summary</h2>
            </div>
            <div class="card-body">
                <div id="summary-info">

                </div>
                <div class="plot" id="gantt-plot">
                </div>
            </div>
        `);
    }


    updateSummary(infoDict) {
        const summaryDiv = document.getElementById("summary-info");
        summaryDiv.innerHTML = ''; // Clear the div
    
        const summaryUl = document.createElement('ul');
        summaryUl.appendChild(createLi(infoDict['totalNumberCopies'] === 0 ? "No query sequence information added." : `The query sequence will be composed by ${infoDict['totalNumberCopies']} sequences concatenated with glycine as linkers.`));

        const templatesLi = document.createElement('li');
        templatesLi.textContent = infoDict['msaTemplatesNum'] === 0 && infoDict['templatesTemplatesNum'] === 0 ? "No information added in templates section." : "Templates information:";
    
        if (infoDict['msaTemplatesNum'] !== 0 || infoDict['templatesTemplatesNum'] !== 0) {
            const innerUl = document.createElement('ul');
            innerUl.appendChild(createLi(infoDict['msaTemplatesNum'] === 0 ? "No MSA sequences added using templates section." : `${infoDict['msaTemplatesNum']} MSA sequences added.`));
            innerUl.appendChild(createLi(infoDict['templatesTemplatesNum'] === 0 ? "No templates added using templates section." : `${infoDict['templatesTemplatesNum']} templates added.`));
            templatesLi.appendChild(innerUl);
        }
        summaryUl.appendChild(templatesLi);
        const librariesLi = document.createElement('li');
        librariesLi.textContent = infoDict['msaLibrariesNum'] === 0 && infoDict['templatesLibrariesNum'] === 0 ? "No information added in libraries section." : "Libraries information:";
        if (infoDict['msaLibrariesNum'] !== 0 || infoDict['templatesLibrariesNum'] !== 0) {
            const innerUl = document.createElement('ul');
            innerUl.appendChild(createLi(infoDict['msaLibrariesNum'] === 0 ? "No MSA sequences added using libraries section." : `${infoDict['msaLibrariesNum']} MSA sequences added.`));
            innerUl.appendChild(createLi(infoDict['templatesLibrariesNum'] === 0 ? "No templates added using libraries section." : `${infoDict['templatesLibrariesNum']} templates added.`));
            librariesLi.appendChild(innerUl);
        }        
        summaryUl.appendChild(librariesLi);
        const featuresLi = document.createElement('li');
        featuresLi.textContent = infoDict['msaFeaturesNum'] === 0 && infoDict['templatesFeaturesNum'] === 0 ? "No information added in features section." : "Features information:";
        if (infoDict['msaFeaturesNum'] !== 0 || infoDict['templatesFeaturesNum'] !== 0) {
            const innerUl = document.createElement('ul');
            innerUl.appendChild(createLi(infoDict['msaFeaturesNum'] === 0 ? "No MSA sequences added using features section." : `${infoDict['msaFeaturesNum']} MSA sequences added.`));
            innerUl.appendChild(createLi(infoDict['templatesFeaturesNum'] === 0 ? "No templates added using features section." : `${infoDict['templatesFeaturesNum']} templates added.`));
            featuresLi.appendChild(innerUl);
        }        
        summaryUl.appendChild(featuresLi);
        summaryDiv.appendChild(summaryUl);
    }
}

function getTemplateRestrictions(templateID){
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
            modificationsDict[selection][positionKey]["delete"] = new Set(...oldValues, [...newValues]);
        }
        const changeAminoacids = div.querySelectorAll(`li[id^=modaminoacids-${templateID}-]`);
        changeAminoacids.forEach(change => {
            const positions = change.querySelector(`input[id^=template-modify-amino-pos-${templateID}]`).value.trim();
            const sequence = change.querySelector(`textarea[id^=template-modify-amino-fasta-${templateID}-]`).value;
            const choose = change.querySelector(`select[id^=template-modify-amino-select-${templateID}-]`).value;
            if(positions && ((choose === "residue") || (choose === "fasta" && sequence.length > 0))){
                const positionsArray = extendedNumbers(positions);
                modificationsDict[selection][positionKey][choose] = [...new Set([...(modificationsDict[selection][positionKey][choose] || []), ...positionsArray])];
            }
        });
    });
    return [modificationsDict, positionsAccepted];
}

async function updatePlot() {
    const gatherInfo = {};
    let msaTemplatesNum = 0;
    let templatesTemplatesNum = 0;
    let msaFeaturesValues = 0;
    let templatesFeaturesValues = 0;
    let maxY = 10;
    const sequenceDiv = document.getElementById('sequence');
    const sequenceComponents = Array.from(sequenceDiv.getElementsByTagName('sequence-component'));
    let querySequence = '';
    let totalNumberCopies = 0;
    let queryArray = [];
    const linkerLength = 50;
    const seqColor = "rgba(46,117,182,1)";
    const linkerColor = "rgba(127,127,127,1)";
    const backgroundQueryColor = "rgb(251, 247, 230, 1)";
    const backgroundTemplateColor = "rgb(220, 237, 233)";
    const backgroudMSAColor = "rgb(215, 233, 203)";
    const xTicks = [];
    const valTicks = [0.5];
    const textTicks = ["Sequence"];
    const guidedMode = document.getElementById('guided-radio').checked;
    let numberResidues = 1;

    for (const component of sequenceComponents) {
        const sequenceID = component.getAttribute('sequence-id');
        const selected = component.querySelector(`input[name=sequence-input-${sequenceID}]:checked`);
        let sequenceInput = '';
        let nameQuery = `seq-${sequenceID}`;
        if (selected.value === 'file') {
            const fastaFileInput = document.getElementById(`sequence-fasta-${sequenceID}`)?.files?.[0];;
            if (fastaFileInput != null && fastaFileInput != 'undefined') {
                try {
                    nameQuery = fastaFileInput.name.split('.')[0].substring(0, 5);
                    const readSeq = await readFile(fastaFileInput);
                    const lines = readSeq.split('\n');
                    if (lines.length > 0) {
                        sequenceInput = lines[0].startsWith('>') ? lines[1] : lines[0];
                    }
                } catch (e) {
                    console.error(e);
                }
            }
        } else {
            sequenceInput = document.getElementById(`sequence-text-${sequenceID}`).value;
        }
        let numberCopiesInput = document.getElementById(`sequence-copies-${sequenceID}`).value;
        totalNumberCopies += parseInt(numberCopiesInput);

        if(sequenceInput === "") continue;

        let mutationsArray = [];
        const ulMutation = document.getElementById(`mutations-container-${sequenceID}`);
        const posmutationInputs = ulMutation.querySelectorAll(`input[id^=sequence-mutations-pos-${sequenceID}-]`);
        posmutationInputs.forEach(input => {
            inputValue = input.value;
            if(input.value !== ""){
                const resultArray = extendedNumbers(input.value);
                mutationsArray = [...new Set([...mutationsArray, ...resultArray])];
            }
        });
        const sequenceDict = {"name": nameQuery, "seq": sequenceInput, "mutations": mutationsArray.sort(function(a, b){return a-b})};
        const numUndefined = queryArray.filter(queryArray => queryArray === undefined).length;
        if (numUndefined < numberCopiesInput) {
            queryArray = queryArray.concat(Array(numberCopiesInput - numUndefined).fill(undefined));
        }
        queryArray = queryArray.reduce((acc, curr, i) => {
            if (curr === undefined && numberCopiesInput > 0) {
                acc[i] = sequenceDict;
                numberCopiesInput--;
            } else {
                acc[i] = curr;
            }
            return acc;
        }, []);
    }

    updatePositionsModifications(totalNumberCopies);

    queryArray = queryArray.filter(item => item);
    const dataShape = queryArray.reduce((acc, queryDict, index) => {
        acc.push({
            type: 'rect',
            x0: numberResidues,
            x1: numberResidues + queryDict.seq.length,
            y0: 0,
            y1: 1,
            fillcolor: seqColor,
            label: {
                text: queryDict.name,
                font: { size: 12, color: 'black' },
                textposition: 'center',
            },
            line: {
                width: 1,
                color: seqColor,
            },
        });
        queryDict.mutations.forEach(value => {
            if(value < numberResidues + queryDict.seq.length){
                acc.push({
                    type: 'rect',
                    x0: numberResidues+value-1,
                    x1: numberResidues+value,
                    y0: 0.8,
                    y1: 1.05,
                    fillcolor: "yellow",
                    line: {
                        width: 1,
                        color: "yellow",
                    },
                });
            }
        });
        xTicks.push(numberResidues)
        numberResidues += queryDict.seq.length;
        querySequence += queryDict.seq;
        xTicks.push(numberResidues)
        if (index !== queryArray.length - 1) {
            acc.push({ 
                type: 'rect',
                x0: numberResidues,
                x1: numberResidues + linkerLength,
                y0: 0.5,
                y1: 0.5,
                fillcolor: linkerColor,
                line: {
                    width: 5,
                    color: linkerColor,
                    dash: 'dot'
                },
            });
            numberResidues += linkerLength;
            querySequence += '-'.repeat(linkerLength);
        }
        return acc;
    }, []);

    if(dataShape.length > 0){
        if(guidedMode){
            let msaSequenceTemplates = new Array(numberResidues).fill(0);
            const arrayTemplates = [];
            const msaSequences = [];
            for (const [key, values] of Object.entries(templatesDict)) {
                const addMSA = document.getElementById(`template-addmsa-${key}`).checked;
                const addTemplate = document.getElementById(`template-addtemplates-${key}`).checked;
                const chainsDict = Object.assign({}, values.chainSequences);
                const resultsArray = Array(totalNumberCopies).fill("");
                const resultsChainArray = Array(totalNumberCopies).fill("");
                let changesResidue = Array(totalNumberCopies).fill([]);
                let changesFasta = Array(totalNumberCopies).fill([]);
                let deleteResidues = Array(totalNumberCopies).fill(new Set());

                const resValues = getTemplateRestrictions(key);
                const positionsAccepted = resValues[1]
                const restrictionsDict = resValues[0]

                if(!positionsAccepted){
                    let currentIndex = 0;
                    Object.keys(chainsDict).forEach((key) => {
                        const valueList = chainsDict[key];
                        for (let i = 0; i < valueList.length && currentIndex < totalNumberCopies; i++) {
                            resultsArray[currentIndex] = valueList[i];
                            resultsChainArray[currentIndex] = key;
                            currentIndex++;
                        }
                    });

                }

                for (const [restrictChain, restrictions] of Object.entries(restrictionsDict)) {
                    if (restrictChain === 'all') {
                        if (restrictions['ANY'].hasOwnProperty('residue')){
                            changesResidue = changesResidue.map(change =>
                                new Set([...change, ...restrictions['ANY']['residue']])
                            );
                        }
                        if (restrictions['ANY'].hasOwnProperty('fasta')){
                            changesFasta = changesFasta.map(change =>
                                new Set([...change, ...restrictions['ANY']['fasta']])
                            );
                        }
                        if (restrictions['ANY'].hasOwnProperty('delete')){
                            deleteResidues = deleteResidues.map(deleteres =>
                                new Set([...deleteres, ...restrictions['ANY']['delete']])
                            );
                        }
                        if (restrictions['ANY'].hasOwnProperty('residue')){
                            changesResidue = changesResidue.map(changes =>
                                new Set([...changes, ...restrictions['ANY']['residue']])
                            );
                        }
                        if (restrictions['ANY'].hasOwnProperty('fasta')){
                            changesFasta = changesFasta.map(fasta =>
                                new Set([...fasta, ...restrictions['ANY']['fasta']])
                            );
                        }
                    }
                    else{
                        for (let [positionStr, values] of Object.entries(restrictions)) {
                            if (!positionsAccepted) {
                                for (let i = 0; i < resultsChainArray.length; i++) {
                                    if (resultsChainArray[i] === restrictChain) {
                                        if (values.hasOwnProperty('residue')) {
                                            changesResidue[i] = new Set([...changesResidue[i], ...values['residue']]);
                                        }
                                        if (values.hasOwnProperty('fasta')) {
                                            changesFasta[i] = new Set([...changesFasta[i], ...values['fasta']]);
                                        }
                                        if (values.hasOwnProperty('delete')) {
                                            deleteResidues[i] = new Set([...deleteResidues[i], ...values['delete']]);
                                        }
                                    }
                                }
                            } else {
                                if(positionStr !== 'ANY'){
                                    let position = parseInt(positionStr) - 1;
                                    const sequence = chainsDict[restrictChain].shift();
                                    chainsDict[restrictChain].push(sequence);
                                    resultsArray[position] = sequence;
                                    resultsChainArray[position] = restrictChain;
                                    if (values.hasOwnProperty('residue')) {
                                        changesResidue[position] = new Set([...changesResidue[position], ...values['residue']]);
                                    }
                                    if (values.hasOwnProperty('fasta')) {
                                        changesFasta[position] = new Set([...changesFasta[position], ...values['fasta']]);
                                    }
                                    if (values.hasOwnProperty('delete')) {
                                        deleteResidues[position] = new Set([...deleteResidues[position], ...values['delete']]);
                                    }
                                }
                            }
                        }
                    }
                }

                let trimmedSequences = '';
                for (let i = 0; i < xTicks.length; i += 2) {
                    changesResidue[i/2] = Array.from(changesResidue[i/2]).map(value => value + (xTicks[i]-1));
                    changesFasta[i/2] = Array.from(changesFasta[i/2]).map(value => value + (xTicks[i]-1));
                    const chainName = resultsChainArray[i / 2];
                    let sequence = resultsArray[i / 2];
                    sequence = Array.from(sequence).map((char, index) => deleteResidues[i / 2].has(index + 1) ? '-' : char).join('');
                    const seqLength = xTicks[i+1] - xTicks[i];
                    let trimmedSeq = chainName !== '' && sequence
                    ? seqLength >= sequence.length ? sequence + '-'.repeat(seqLength - sequence.length) : sequence.substring(0, seqLength)
                    : '-'.repeat(seqLength);
                    trimmedSequences += trimmedSeq + '-'.repeat(linkerLength);
                }

                if(addTemplate){
                    templatesTemplatesNum += 1;
                    let sequence = compareSequences(querySequence, trimmedSequences);
                    sequence = sequence.map(value => value / 6);
                    arrayTemplates.push({
                        'name':values.templateName,
                        'seq': sequence,
                        'changesRes': changesResidue.flat().sort(function(a, b){return a-b}),
                        'changesFasta': changesFasta.flat().sort(function(a, b){return a-b})
                    });
                }
                if(addMSA){
                    msaTemplatesNum += 1;
                    const msaSeq = compareSequences(querySequence, trimmedSequences);
                    msaSeq.forEach((value, index) => msaSequenceTemplates[index] += value);
                }
            }

            if(msaTemplatesNum > 0){
                msaSequenceTemplates = scaleValues(msaSequenceTemplates);
                msaSequences.push([msaTemplatesNum,  msaSequenceTemplates]);
            }


            for (const [key, values] of Object.entries(librariesDict)) {
                const positionLibrary =  document.getElementById(`library-lib-${key}`).value;
                const positionQuery =  document.getElementById(`library-query-${key}`).value;
                let seqArray = Array(numberResidues).fill(0);
                let positionLibraryArray = [];
                if(positionLibrary.trim() !== "" && positionQuery.trim() !== ""){
                    positionLibrarySplit = positionLibrary.split(',');
                    const queryPositions = extendedNumbers(positionQuery);
                    queryPositions.forEach((start, i) => {
                        const rangeLength = extendedNumbers(positionLibrarySplit[i]).length;
                        for (let j = 0; j < rangeLength; j++) {
                            const pos = start + j;
                            positionLibraryArray.push(pos);
                        }
                    });
                } else {
                    positionLibraryArray = extendedNumbers(`1-${numberResidues-1}`);
                }
                positionLibraryArray.forEach(index => seqArray[index-1] = 1);
                arrayTemplates.push({
                    'name': `lib-${key}`,
                    'seq': seqArray,
                    'changesRes': [],
                    'changesFasta': []
                });
            }
            for (const [key, values] of Object.entries(featuresDict)) {
                const addMSA =  document.getElementById(`feature-addmsa-${key}`).checked;
                const addTemplates =  document.getElementById(`feature-addtemplates-${key}`).checked;
                const featPosition = document.getElementById(`feature-pos-${key}`).value;
                const shiftingArray = new Array(xTicks[(featPosition-1)*2]-1).fill(0);
                const rangeFeatures = document.getElementById(`feature-regionfeat-${key}`).value;
                const rangeQuery =  document.getElementById(`feature-regionquery-${key}`).value;
                let arrayTemp = new Array(numberResidues).fill(0);
                let arrayMsa = new Array(numberResidues).fill(0);

                if(rangeFeatures && rangeQuery){
                    posQuery = rangeQuery.split(',').map(Number);
                    posFeat = rangeFeatures.split(',');
                    for (let i = 0; i < posQuery.length; i++) {
                        let rangeFeat = posFeat[i].split('-');
                        let startFeat = parseInt(rangeFeat[0]) - 1;
                        let endFeat = parseInt(rangeFeat[1])-1 || startFeat;
                        let posQ = posQuery[i] - 1;
                        if (startFeat < 0) startFeat = 0;
                        if (endFeat > numberResidues) endFeat = numberResidues;
                        for (let j = startFeat; j <= endFeat; j++) {
                            let newIndex = posQ + j - startFeat;
                            if (j >= 0 && j < values['templates_coverage'].length) {
                                arrayTemp[newIndex] = values['templates_coverage'][j];
                                arrayMsa[newIndex] = values['msa_coverage'][j];
                            }
                        }
                    }
                }
                else{
                    arrayTemp = [...shiftingArray, ...values['templates_coverage'].map(value => value)];
                    arrayMsa = [...shiftingArray, ...values['msa_coverage'].map(value => value)];
                }
                if(addTemplates){
                    arrayTemplates.push({
                        'name': `feat-${key}`,
                        'seq': arrayTemp,
                        'changesRes': [],
                        'changesFasta': []
                    });
                    templatesFeaturesValues += values['num_templates'];
                }
                if(addMSA){
                    const numMsa = values['num_msa'];
                    msaFeaturesValues += values['num_msa'];
                    msaSequences.push([numMsa, arrayMsa]);
                }
            }

            const hasMSA = msaSequences.length > 0;
            if(hasMSA){
                const totalSum = msaSequences.reduce((sum, item) => sum + item[0], 0);
                const aggregatedMSA = msaSequences.reduce((result, [scale, array]) => {
                    const scaledArray = array.map(num => num * (scale / totalSum));
                    return result.map((sum, idx) => sum + (scaledArray[idx] || 0));
                }, new Array(numberResidues).fill(0));
                const posVertical = valTicks[valTicks.length - 1] + 2;
                aggregatedMSA.forEach((value, index) => {
                    const color = gradient(value);
                    if(value !== 0 && index+1 < numberResidues){
                        dataShape.push({
                            type: 'rect',
                            x0: index+1,
                            x1: index+2,
                            y0: posVertical-0.5,
                            y1: posVertical+0.5,
                            fillcolor: color,
                            line: {
                                width: 1,
                                color: color,
                            },
                        });
                    }
                });
                valTicks.push(posVertical);
                textTicks.push('MSA');
            }

            const hasTemplates = arrayTemplates.length > 0;
            if(hasTemplates){
                arrayTemplates.forEach(templateDict => {
                    const posVertical = valTicks[valTicks.length - 1] + 2;
                    templateDict.seq.forEach((seqValue, index) => {
                        if(seqValue !== 0 && index+1 < numberResidues){
                            dataShape.push({
                                type: 'rect',
                                x0: index+1,
                                x1: index+2,
                                y0: posVertical-0.5,
                                y1: posVertical+0.5,
                                fillcolor: gradient(seqValue),
                                line: {
                                    width: 1,
                                    color: gradient(seqValue),
                                },
                            });
                            if(templateDict.changesRes.includes(index+1)){
                                dataShape.push({
                                    type: 'rect',
                                    x0: index+1,
                                    x1: index+2,
                                    y0: posVertical+0.3,
                                    y1: posVertical+0.55,
                                    fillcolor: "yellow",
                                    line: {
                                        width: 1,
                                        color: "yellow",
                                    },
                                });
                            }
                            if(templateDict.changesFasta.includes(index+1)){
                                dataShape.push({
                                    type: 'rect',
                                    x0: index+1,
                                    x1: index+2,
                                    y0: posVertical+0.3,
                                    y1: posVertical+0.55,
                                    fillcolor: "red",
                                    line: {
                                        width: 1,
                                        color: "red",
                                    },
                                });
                            }
                        }
                    });
                    valTicks.push(posVertical);
                    textTicks.push(templateDict.name.substring(0,6));
                });
            }

            maxY = valTicks[valTicks.length - 1] > 8 ? valTicks[valTicks.length - 1] + 2 : maxY;
            dataShape.push({
                type: 'rect',
                x0: 1-numberResidues*0.03,
                x1: numberResidues+numberResidues*0.03,
                y0: -0.5,
                y1: hasMSA || hasTemplates ? 1.5 : maxY,
                fillcolor: backgroundQueryColor,
                layer: 'below',
                line: {
                    width: 0
                }
            });

            if(hasMSA){
                dataShape.push({
                    type: 'rect',
                    x0: 1-numberResidues*0.03,
                    x1: numberResidues+numberResidues*0.03,
                    y0: 1.5,
                    y1: hasTemplates ? 3.5 : maxY,
                    fillcolor: backgroudMSAColor,
                    layer: 'below',
                    line: {
                        width: 0
                    }
                });
            }
            if(hasTemplates){
                dataShape.push({
                    type: 'rect',
                    x0: 1-numberResidues*0.03,
                    x1: numberResidues+numberResidues*0.03,
                    y0: hasMSA ? 3.5 : 1.5,
                    y1: maxY,
                    fillcolor: backgroundTemplateColor,
                    layer: 'below',
                    line: {
                        width: 0
                    }
                });
            }
        } else {
            dataShape.push({
                type: 'rect',
                x0: 1-numberResidues*0.03,
                x1: numberResidues+numberResidues*0.03,
                y0: -0.5,
                y1: maxY,
                fillcolor: backgroundQueryColor,
                layer: 'below',
                line: {
                    width: 0
                }
            });
        }
    }
    const layout = {
        title: 'TEMPLATES and ALIGNED SEQUENCES (MSA)',
        xaxis: {
            title: 'Residue number',
            range: [1-numberResidues*0.03, numberResidues+numberResidues*0.03],
            tickvals: xTicks,
            tickmode: 'array',
            tickangle: 45,
            gridcolor: 'black',
            zeroline: false,
            showline: false,
        },
        yaxis: {
            title: 'Information',
            tickvals: valTicks,
            ticktext: textTicks,
            range: [-0.5, maxY],
            showgrid: false,
            zeroline: false,
            showline: false,
            autotick: false,
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
    
    Plotly.newPlot('gantt-plot', [], layout, config);

    let {msaLibraryValues, templatesLibraryValues} = Object.values(librariesDict).reduce((acc, val) => {
        if (val.msa) acc.msaLibraryValues += val.msa;
        if (val.templates) acc.templatesLibraryValues += val.templates;
        return acc;
    }, {msaLibraryValues: 0, templatesLibraryValues: 0});

    gatherInfo["msaLibrariesNum"] = msaLibraryValues;
    gatherInfo["templatesLibrariesNum"] = templatesLibraryValues;
    gatherInfo["msaTemplatesNum"] = msaTemplatesNum;
    gatherInfo["templatesTemplatesNum"] = templatesTemplatesNum;
    gatherInfo["msaFeaturesNum"] = msaFeaturesValues;
    gatherInfo["templatesFeaturesNum"] = templatesFeaturesValues;
    gatherInfo["totalNumberCopies"] = totalNumberCopies;

    mainModule.getSummaryClass().updateSummary(gatherInfo);
}

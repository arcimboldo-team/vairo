const aminoacidMap = {
  'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
  'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
  'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
  'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
};
const aminoacidKeys = Object.keys(aminoacidMap);
const aminoacidSelect = aminoacidKeys.map((aminoacid) => `<option value="${aminoacid}">${aminoacid}</option>`).join('');

function convertToOneLetter(threeLetterCode) {
  return aminoacidMap[threeLetterCode] || '-';
}

let selectPositionsArray = ['ANY'];
function updatePositionsModifications(totalPositions){
    selectPositionsArray = ['ANY', ...Array.from({length: totalPositions}, (_, i) => i + 1)];
    const optionsHTML = selectPositionsArray.map(option => `<option value="${option}">${option}</option>`).join('');
    const selectionsTemplate = document.querySelectorAll(`select[id^=template-modify-pos-]`);
    selectionsTemplate.forEach(select => {
        const oldValue = select.value;
        select.innerHTML = optionsHTML;
        if (selectPositionsArray.map(String).includes(String(oldValue))) {
            select.value = oldValue;
        } else{
            select.value = 'ANY';
        }
    });

    const selectionsFeature = document.querySelectorAll(`select[id^=feature-pos-]`);
    const selectFeatureArray = selectPositionsArray = [...Array.from({length: totalPositions}, (_, i) => i + 1)];
    const optionsFeatureHTML = selectFeatureArray.map(option => `<option value="${option}">${option}</option>`).join('');
    selectionsFeature.forEach(select => {
        const oldValue = select.value;
        select.innerHTML = optionsFeatureHTML;
        if (selectFeatureArray.map(String).includes(String(oldValue))) {
            select.value = oldValue;
        } else{
            select.value = '1';
        }
    });
}


function extendedNumbers(input){
  inputArray = input.split(',');
  return inputArray.flatMap((str) => {
    if (str.includes("-")) {
      const [start, end] = str.split("-");
      return Array.from({ length: end - start + 1 }, (_, i) => i + parseInt(start));
    } else {
      return parseInt(str);
    }
  });
}

function getGroup(res) {
  let groups = ['GAVLI', 'FYW', 'CM', 'ST', 'KRH', 'DENQ', 'P'];
  let group = groups.filter(s => s.includes(res));
  if (group.length > 0) {
      return group[0];
  }
  return res;
}

function compareSequences(sequence1, sequence2) {
  let return_list = [];
  for (let i = 0; i < sequence1.length; i++) {
      if (i < sequence2.length) {
          let res1 = sequence1.charAt(i);
          let res2 = sequence2.charAt(i);
          if (res2 === '-') {
              return_list.push(0);
          }
          else if (res1 === res2) {
              return_list.push(6);
          } else if (getGroup(res1) === getGroup(res2)) {
              return_list.push(4);
          } else {
              return_list.push(2);
          }
      } else {
          return_list.push(0);
      }
  }
  return return_list;
}

function gradient(value) {
  let color = Math.floor((1-value) * 255);
  return "rgb(" + color + "," + color + "," + color + ")";
}


function scaleValues(inputList) {
  let max_value = Math.max(...inputList);
  let new_list = [];
  for (let value of inputList) {
    let new_value;
    if (value <= 0) {
        new_value = 0;
    } else if (value >= max_value) {
        new_value = 1;
    } else {
        new_value = Math.round(Math.log(value + 1) / Math.log(max_value + 1) * 100) / 100;
    }
    new_list.push(new_value);
  }
  return new_list;
}

function getSequences(pdbContent) {
  const lines=pdbContent.split('\n');
  const sequences={};
  let currentChain='';
  for(let i=0;i<lines.length;i++){
      if(lines[i].startsWith('ATOM')||lines[i].startsWith('HETATM')){
          const chainID=lines[i][21].trim();
          const residueNumber=+lines[i].substring(22,26).trim();
          let residueName=lines[i].substring(17,20).trim();
          if(lines[i].startsWith('HETATM')){
              if (residueName !== 'MSE') {
                  continue;
              }
              residueName = 'MET';            
          }
          const oneLetterCode=convertToOneLetter(residueName);
          if(chainID!==currentChain){
              currentChain=chainID;
              sequences[currentChain]='';
          }
          if(sequences[currentChain].length<residueNumber){
              sequences[currentChain]+='-'.repeat(residueNumber-sequences[currentChain].length-1)+oneLetterCode;
          }
      }
  }
  return sequences;
}

async function readFile(file, binary=false) {
  return new Promise((resolve, reject) => {
      let fileReader = new FileReader();
      if(binary){
          fileReader.readAsArrayBuffer(file);
      } else {
          fileReader.readAsText(file);
      }
      fileReader.onload = (e) => resolve(fileReader.result);
      fileReader.onerror = (e) => reject(e);
  });
}

async function fetchPDB(pdbCode) {    
    try {
        const response = await fetch(`https://files.rcsb.org/download/${pdbCode}.pdb`);
        const pdbData = await response.text();
        return pdbData;
    } catch (error) {
        console.error('Error fetching PDB:', error);
        alert('Error fetching PDB. Please try again later.');
        document.getElementById(`template-code-${id}`).value = "";
    }
}

async function postData(url = '', dataDict, jsonData = true, alreadyForm = false) {
    if (!url) throw new Error('URL is required');
    let formData;
    if (alreadyForm) {
        formData = dataDict;
    } else if (typeof dataDict === 'object' && dataDict !== null) {
        formData = new FormData();
        for (const [key, value] of Object.entries(dataDict)) {
            formData.append(key, value);
        }
    } else {
        throw new Error('dataDict must be an object or FormData');
    }
    let response;
    try {
        response = await fetch(url, {
            method: 'POST',
            body: formData,
        });
    } catch (error) {
        throw new Error(`Network error: ${error.message}`);
    }
    if (!response.ok) {
        const status = response.status;
        let errorMessage = `Server error: ${status} ${response.statusText}`;
        throw new Error(errorMessage);
    }
    try {
        return jsonData ? await response.json() : await response.blob();
    } catch (error) {
        throw new Error(`Error parsing response: ${error.message}`);
    }
}

function enableNavButtons(ymlExists, htmlExists) {
    localStorage.setItem('inputButtonEnabled', ymlExists);
    localStorage.setItem('outputButtonEnabled', htmlExists);
    const inputButton = document.getElementById('input-button');
    const outputButton = document.getElementById('output-button');
    if (ymlExists) {
        inputButton.classList.remove('disabled');
    } else {
        inputButton.classList.add('disabled');
    }
    if (htmlExists) {
        outputButton.classList.remove('disabled');
    } else {
        outputButton.classList.add('disabled');
    }
}

function restoreButtonStates() {
    const inputButton = document.getElementById('input-button');
    const outputButton = document.getElementById('output-button');

    const inputEnabled = localStorage.getItem('inputButtonEnabled') === 'true';
    const outputEnabled = localStorage.getItem('outputButtonEnabled') === 'true';

    if (inputButton && inputEnabled) {
        inputButton.classList.remove('disabled');
    }

    if (outputButton && outputEnabled) {
        outputButton.classList.remove('disabled');
    }
}

document.addEventListener('DOMContentLoaded', restoreButtonStates);




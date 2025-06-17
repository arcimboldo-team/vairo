window.addEventListener('DOMContentLoaded', event => {
    const vairoForm = document.getElementById('vairo-form');
    if (vairoForm) {
        vairoForm.addEventListener('submit', handleFormSubmission);
    }
    loadSessionData();
});

window.addEventListener('beforeunload', async event => {
    const form = document.getElementById('vairo-form');
    const formData = new FormData(form);
    const results = await postData('/save-form-data', formData, true, true);
});


function checkMode(){
    var guidedRadio = document.getElementById('guided-radio');
    var naiveRadio = document.getElementById('naive-radio');
    var element1 = document.getElementById('template');
    var element2 = document.getElementById('library');
    var element3 = document.getElementById('feature');

    if (naiveRadio.checked) {
        element1.classList.add('hidden');
        element2.classList.add('hidden');
        element3.classList.add('hidden');
    } else {
        element1.classList.remove('hidden');
        element2.classList.remove('hidden');
        element3.classList.remove('hidden');
    }
    updatePlot();
}

async function handleFormSubmission(event) {
    event.preventDefault();
    clearErrors();
    const errors = await validateForm(event.target);
    if (errors.length > 0) {
        displayErrors(errors);
        document.getElementById("error-container").scrollIntoView({ behavior: "smooth", block: "end" });
        return;
    }
    const formData = new FormData(event.target);
    const results = await postData('/form-vairo', formData, true, true);
    if (results) {
        window.location.href = "/input";
    } else {
        displayErrors(["A server error occurred while submitting the form. Please try again later."]);
        document.getElementById("error-container").scrollIntoView({ behavior: "smooth", block: "end" });
        return;
    }
}

async function validateForm(form) {
  const errors = [];
  const inputDir = document.getElementById('general-output').value;
  if(inputDir === ""){
    errors.push("Output folder not found.");
  }
  const folderPath = document.getElementById('general-databases').value;
  const resultDict = await postData('/check-databases', { 'folder': folderPath });
  if(!resultDict.exists){
    errors.push("AlphaFold2 database not found.");
  }

  const sequenceComponents = Array.from(form.querySelectorAll('sequence-component'));
  if (sequenceComponents.length === 0) {
    errors.push("No query sequences added. There has to be at least one sequence.");
  }
  return errors;
}

function displayErrors(errors) {
  const errorContainer = document.getElementById('error-container');
  const errorMessages = document.getElementById('error-messages');
  if (errorContainer && errorMessages) {
    let html = '<ul>';
    errors.forEach(error => {
      html += `<li>${error}</li>`;
    });
    html += '</ul>';
    errorMessages.innerHTML = html;
    errorContainer.style.display = 'block';
  }
}

function clearErrors() {
  const errorContainer = document.getElementById('error-container');
  const errorMessages = document.getElementById('error-messages');
  if (errorContainer && errorMessages) {
    errorMessages.innerHTML = '';
    errorContainer.style.display = 'none';
  }
}

function restoreNestedFields(section, data) {
    for (const [id, fields] of Object.entries(data)) {
        for (const [fieldName, value] of Object.entries(fields)) {
            const element = document.querySelector(`[name="${section}[${id}][${fieldName}]"]`);
            if (element) {
                element.value = value;
                if (element.type === 'checkbox' || element.type === 'radio') {
                    element.checked = true;
                }
            }
        }
    }
}

async function loadSessionData() {
    try {
        const response = await fetch('/load-form-data');
        const result = await response.json();
        for (const [key, value] of Object.entries(result.data)) {
                if (typeof value === 'object' && value !== null) {
                    restoreNestedFields(key, value);
                } else {
                    const element = document.getElementById(key) || document.querySelector(`[name="${key}"]`);
                    if (element) {
                        element.value = value;
                        if (element.type === 'checkbox' || element.type === 'radio') {
                            element.checked = true;
                        }
                    }
                }
            }
        parameterizationFolderChanged();
        databaseFolderChanged();
    } catch (error) {
        console.error('Error loading session data:', error);
    }
}

const mainModule = (function () {
    let sequenceClass, templateClass, libraryClass, featureClass, summaryClass;

    $(document).ready(function(){
        sequenceClass = new Accordion("sequence", "Query sequences");
        templateClass = new Accordion("template", "Templates");
        libraryClass = new Accordion("library", "Libraries");
        featureClass = new Accordion("feature", "Features");
        summaryClass = new Summary();
        updatePlot();
  });

    return {
        getSequenceClass: function() { return sequenceClass; },
        getTemplateClass: function() { return templateClass; },
        getLibraryClass: function() { return libraryClass; },
        getFeatureClass: function() { return featureClass; },
        getSummaryClass: function() { return summaryClass; }
    };
})();

async function parameterizationFolderChanged() {
    const inputElement = document.getElementById('general-output');
    const folderPath = inputElement.value.trim();
    clearValidationMessage('general-output');
    const resultDict = await postData('/check-output', { folder: folderPath });
    if (resultDict.yml_exists || resultDict.html_exists) {
        updateValidationMessage('general-output', 'warning', 'Folder already exists, check Input and Output tab to see the input file and the results. The results will be overwritten.');
    }
    enableNavButtons(resultDict.yml_exists, resultDict.html_exists);
}

async function databaseFolderChanged() {
    const inputElement = document.getElementById('general-databases');
    const folderPath = inputElement.value.trim();
    clearValidationMessage('general-databases');
    const resultDict = await postData('/check-databases', { 'folder': folderPath });
    if(resultDict.exists){
        updateValidationMessage('general-databases', 'ok', 'Databases found');
    } else {
        updateValidationMessage('general-databases', 'nok', 'Error validating databases, check folder.');
    }
}

function updateValidationMessage(fieldId, type, message) {
    const inputElement = document.getElementById(fieldId);
    let messageElement = document.getElementById(fieldId + '-validation');

    if (!messageElement) {
        messageElement = document.createElement('small');
        messageElement.id = fieldId + '-validation';
        messageElement.className = 'form-text';
        inputElement.parentNode.appendChild(messageElement);
    }
    if (type == 'ok') {
        messageElement.className = 'form-text text-success';
        messageElement.innerHTML = '✓ ' + message;
    } else if (type == 'warning') {
        messageElement.className = 'form-text text-warning';
        messageElement.innerHTML = '⚠ ' + message;
    } else {
        messageElement.className = 'form-text text-danger';
        messageElement.innerHTML = '✗ ' + message;
    }
}

function clearValidationMessage(fieldId) {
    const messageElement = document.getElementById(fieldId + '-validation');
    if (messageElement) {
        messageElement.remove();
    }
}
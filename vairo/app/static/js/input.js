
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
}

async function handleFormSubmission(event) {
  event.preventDefault();
  clearErrors();
  const errors = await validateForm(event.target);
  console.log(errors);
  console.log(errors.length > 0);
  if (errors.length > 0) {
    console.log('entra');
    displayErrors(errors);
    document.getElementById("error-container").scrollIntoView({ behavior: "smooth", block: "end" });
    return;
  }
  //const formData = new FormData(event.target);
  //const resultFile = await postData('/form-vairo', formData, true, true);
  //if (resultFile) {
  //  const vairoPath = document.getElementById('general-output').value;
  //  const encodedPath = encodeURIComponent(vairoPath);
  //  window.location.href = `/output?data=${encodedPath}`;
  //}
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
  console.log(errorContainer);
  console.log(errorMessages);
  if (errorContainer && errorMessages) {
    console.log('entra')
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

window.addEventListener('DOMContentLoaded', event => {
    const vairoForm = document.getElementById('vairo-form');
    if (vairoForm) {
        vairoForm.addEventListener('submit', handleFormSubmission);
    }    
    const sidebarToggle = document.body.querySelector('#sidebarToggle');
    if (sidebarToggle) {
        sidebarToggle.addEventListener('click', event => {
            event.preventDefault();
            document.body.classList.toggle('sb-sidenav-toggled');
            localStorage.setItem('sb|sidebar-toggle', document.body.classList.contains('sb-sidenav-toggled'));
        });
    }
});

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
async function handleFormSubmission(event) {
    event.preventDefault();

    const formData = new FormData(event.target);
    resultFile = await postData('/form-vairo', formData, true, true);
    if (resultFile) {
        const vairoPath = document.getElementById('general-output').value;
        const encodedPath = encodeURIComponent(vairoPath);
        window.location.href = `/output?data=${encodedPath}`;
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
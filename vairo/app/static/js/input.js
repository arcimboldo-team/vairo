window.addEventListener('DOMContentLoaded', async (event) => {
    try{
        const response = await fetch('/read-yml');
        const result = await response.json();
        const inputTextarea = document.getElementById('input-text');
        if (inputTextarea && result.data !== null) {
            inputTextarea.value = JSON.parse(JSON.stringify(result.data));
        }
    } catch (error) {
        console.error('Error loading session data:', error);
    }
    const inputForm = document.getElementById('input-form');
    if (inputForm) {
        inputForm.addEventListener('submit', handleInputSubmission);
    }
});

async function handleInputSubmission(event) {
  event.preventDefault();
  const formData = new FormData(event.target);
  const resultFile = await postData('/run-vairo', formData, true, true);
  if (resultFile) {
    window.location.href = "/output";
  }
}

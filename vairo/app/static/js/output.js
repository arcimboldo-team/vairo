async function displayOutput(){

    const urlParams = new URLSearchParams(window.location.search);
    const encodedData = urlParams.get('data');

    if (encodedData) {
        const decodedData = decodeURIComponent(encodedData);
        try{
            outputInfo = await postData('/read-output', {'path': decodedData});
            const mainElement = document.getElementById('vairo-output');
            const configPathElement = document.createElement('p');
            configPathElement.textContent = `Config Path: ${outputInfo['config_path']}`;
            mainElement.appendChild(configPathElement);
        
            const configInfoElement = document.createElement('p');
            const configInfoText = outputInfo['config_info'].replace(/\n/g, '<br>');
            configInfoElement.innerHTML = `Config Info:<br> ${configInfoText}`;
            mainElement.appendChild(configInfoElement);
        
            const outputPathElement = document.createElement('p');
            outputPathElement.textContent = `Output Path: ${outputInfo['output_path']}`;
            mainElement.appendChild(outputPathElement);

        } catch (error) {
            alert('It has not been possible to read the output');
            console.error('Error:', error);
            window.location.href = `/`;
        }
    }
}


window.addEventListener('DOMContentLoaded', event => {
    displayOutput();
});
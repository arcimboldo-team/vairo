async function displayOutput(){

    const urlParams = new URLSearchParams(window.location.search);
    const encodedData = urlParams.get('data');

    if (encodedData) {
        const decodedData = decodeURIComponent(encodedData);
        try{
            outputInfo = await postData('/read-output', {'path': decodedData});
            const mainElement = document.getElementById('vairo-output');

            mainElement.innerHTML = `
                <p>Configuration file: ${outputInfo['config_path']}</p>
                <pre>${outputInfo['config_info']}</pre>
                <p>HTML output path: ${outputInfo['output_path']}</p>
            `;

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
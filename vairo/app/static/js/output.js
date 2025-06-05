function startPollingIfElementExists() {
    let consecutiveFailures = 0;
    const maxFailures = 5;

    (function poll() {
        fetch('/check-update')
            .then(response => response.json())
            .then(data => {
                consecutiveFailures = 0;
                console.log(data);
                if (data.changed && data.content) {
                    document.getElementById('vairo-output-iframe').srcdoc = data.content;
                    console.log(data.content);
                }
            })
            .catch(error => {
                consecutiveFailures++;
                console.error('Error checking for updates:', error);
                if (consecutiveFailures >= maxFailures) {
                    console.error('Server could not be located. Stopping polling.');
                    return;
                }
            })
            .finally(() => setTimeout(poll, 10000));
    })();
}

document.addEventListener('DOMContentLoaded', function() {
    startPollingIfElementExists();
});
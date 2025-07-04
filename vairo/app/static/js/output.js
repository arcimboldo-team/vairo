let pollingTimeoutId = null;
let currentController = null;

function startPollingIfElementExists() {
    if (currentController) {
        currentController.abort();
    }

    if (pollingTimeoutId) {
        clearTimeout(pollingTimeoutId);
        pollingTimeoutId = null;
    }

    currentController = new AbortController();
    const signal = currentController.signal;

    let consecutiveFailures = 0;
    const maxFailures = 5;

    (function poll() {
        fetch('/check-update', { signal })
            .then(response => response.json())
            .then(data => {
                consecutiveFailures = 0;
                if (data.error === 'Finished') {
                    console.log('Polling stopped: Finished');
                    currentController.abort();
                    return;
                }
                if (data.changed && data.content) {
                    document.getElementById('vairo-output-iframe').srcdoc = data.content;
                }
            })
            .catch(error => {
                if (error.name === 'AbortError') {
                    console.log('Polling aborted - new instance started');
                    return;
                }
                consecutiveFailures++;
                console.error('Error checking for updates:', error);
                if (consecutiveFailures >= maxFailures) {
                    console.error('Server could not be located. Stopping polling.');
                    return;
                }
            })
            .finally(() => {
                if (!signal.aborted) {
                    pollingTimeoutId = setTimeout(poll, 10000);
                }
            });
    })();
}

document.addEventListener('DOMContentLoaded', function() {
    startPollingIfElementExists();
});
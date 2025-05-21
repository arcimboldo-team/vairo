class Accordion {
    constructor(type, title) {
        this.counter = 1;
        this.type = type;
        this.namePrefix = type.charAt(0).toUpperCase() + type.slice(1);
        this.addButtonID = `${type}-add-button`;
        this.wrapper = `accordion-${this.type}`;

        $(`#${this.type}`).append(`
            <div class="card-header">
                <h2 class="card-title">${title}</h2>
            </div>
            <div class="card-body">
                <div class="accordion" id="accordion-${this.type}">
                </div>
                <div class="col-md-12 text-center" style="margin-top:15px;">
                    <button class="btn btn-success" id="${this.addButtonID}"><span class="glyphicon glyphicon-plus"></span> Add ${this.type}</button>
                </div>
            </div>
        `);

        $(`#${this.addButtonID}`).on("click", (e) => {
            e.preventDefault();
            this.addAccordionItem();
        });

        $(`#${this.wrapper}`).on("click", ".delete-table", function(e) {
            e.preventDefault();
            const id = $(this).closest('.accordion-item').attr('id');
            const parts = id.split('-');
            const idNum = parts[parts.length - 1];
            $(this).closest('.accordion-item').remove();
            deleteInDict(type, idNum)
            updatePlot();
        });
    }

    addAccordionItem() {
        const name = `${this.namePrefix} ${this.counter}`;
        $(`#${this.wrapper}`).append(`
            <div class="accordion-item" id="${this.type}-item-${this.counter}">
                <div class="accordion-heading">
                    <h2 class="accordion-header" id="${this.type}-header-${this.counter}">
                        <div class="accordion-button d-flex delete-icon-margin" type="button" data-bs-toggle="collapse" data-bs-target="#${this.type}-collapse-${this.counter}" aria-expanded="false" aria-controls="${this.type}-collapse-${this.counter}">
                            ${name}
                            <div class="ms-auto">
                                <span class="fa fa-trash-alt delete-table delete-icon-format" data-app-id="${this.type}-item-${this.counter}" data-content-id="${this.type}-item-${this.counter}" aria-hidden="true"></span>
                            </div>
                        </div>
                    </h2>
                    <div id="${this.type}-collapse-${this.counter}" class="accordion-collapse collapse show" aria-labelledby="${this.type}-header-${this.counter}" data-bs-parent="#${this.type}">
                        <div class="accordion-body">
                            <${this.type}-component ${this.type}-id=${this.counter}></${this.type}-component>
                        </div>
                    </div>
                </div>
            </div>
        `);
        this.counter++;
    }
}

// SECTION: Form State Management
const formState = {
    currentTab: 0,
    tabs: ['geometry', 'diffusivity', 'advanced', 'review'],
    completedTabs: new Set()
};

var advancedEnabled = true;
// !SECTION


// SECTION: Tab Navigation
const tabButtons = document.querySelectorAll('.tab-button');
const tabContents = document.querySelectorAll('.tab-content');
const nextBtn = document.getElementById('nextBtn');
const prevBtn = document.getElementById('prevBtn');
const submitBtn = document.getElementById('submitBtn');

function switchTab(tabIndex) {
    if (tabIndex < 0 || tabIndex >= formState.tabs.length) return;

    // Validate current tab before switching
    if (tabIndex > formState.currentTab && !validateTab(formState.currentTab)) {
        showValidationMessage('Please fill in all required fields in this tab', 'error');
        return;
    }

    // Remove active class from all buttons and contents
    tabButtons.forEach(btn => btn.classList.remove('active'));
    tabContents.forEach(content => content.classList.remove('active'));

    // Add active class to new tab
    tabButtons[tabIndex].classList.add('active');
    tabContents[tabIndex].classList.add('active');

    formState.currentTab = tabIndex;
    updateNavigationButtons();
    updateReviewTab();
}

function buttonNextStepSize() {
    if (!advancedEnabled && formState.currentTab === 1) return 2;
    return 1; 
}

function buttonPrevStepSize() {
    if (!advancedEnabled && formState.currentTab === formState.tabs.length - 1) return 2;
    return 1; 
}


tabButtons.forEach((button, index) => {
    button.addEventListener('click', (e) => {
        e.preventDefault();
        switchTab(index);
    });
});

nextBtn.addEventListener('click', () => {
    switchTab(formState.currentTab + buttonNextStepSize());
});

prevBtn.addEventListener('click', () => {
    switchTab(formState.currentTab - buttonPrevStepSize());
});

function updateNavigationButtons() {
    prevBtn.style.display = formState.currentTab > 0 ? 'block' : 'none';
    nextBtn.style.display = formState.currentTab < formState.tabs.length - 1 ? 'block' : 'none';
    submitBtn.style.display = formState.currentTab === formState.tabs.length - 1 ? 'block' : 'none';
}
// !SECTION

// SECTION: Dynamic Label Update for Geometry
const geometryRadioButtons = document.querySelectorAll('input[name="geometryType"]');
const characteristicLengthLabel = document.querySelector('label[for="charlength"]');
const labelMap = {
    slab: 'Length', 
    cylinder: 'Radius',
    sphere: 'Radius',
    custom: 'Volume/Surface Area Ratio'
};

function updateCharacteristicLengthLabel() {
    const selectedGeometry = document.querySelector('input[name="geometryType"]:checked').value;
    characteristicLengthLabel.textContent = labelMap[selectedGeometry];
}

geometryRadioButtons.forEach(radio => {
    radio.addEventListener('change', updateCharacteristicLengthLabel);
});

// Run it at least once on load
updateCharacteristicLengthLabel();

// !SECTION

// SECTION: Dynamic form requirement update for reaction type
// User does not need to specify all concentrations, it depends on reaction type.
// type 2, type 4, type 6 => disable B
// type 3, type 5, type 6 => disable D
const rxTypesRadioButtons = document.querySelectorAll('input[name="rxType"]');
const concBInput = document.getElementById('concB');
const constBLabel = document.querySelector('label[for="concB"]');
const concDInput = document.getElementById('concD');
const constDLabel = document.querySelector('label[for="concD"]');

function updateConcentrationRequirements() {
    const selectedRxType = document.querySelector('input[name="rxType"]:checked').value;

    if (['type2', 'type4', 'type6'].includes(selectedRxType)) {
        concBInput.removeAttribute('required');
        concBInput.setAttribute('disabled', true);
        concBInput.value = '';
        constBLabel.classList.remove('required');
    } else {
        concBInput.setAttribute('required', true);
        concBInput.removeAttribute('disabled');
        constBLabel.classList.add('required');
    }

    if (['type3', 'type5', 'type6'].includes(selectedRxType)) {
        concDInput.removeAttribute('required');
        concDInput.setAttribute('disabled', true);
        concDInput.value = '';
        constDLabel.classList.remove('required');
    } else {
        concDInput.setAttribute('required', true);
        concDInput.removeAttribute('disabled');
        constDLabel.classList.add('required');
    }
}

rxTypesRadioButtons.forEach(radio => {
    radio.addEventListener('change', updateConcentrationRequirements);
});

//!SECTION

// SECTION: Dynamic form requirement update for diffusivity type
// Depends on:
// - Selected diffusivity computation mode (pre-computed)
// - Previously selected reaction type
const diffusivityTypeRadioButtons = document.querySelectorAll('input[name="diffusivityType"]');
const diffusivityInputs = document.querySelectorAll('#diffusivityA, #diffusivityB, #diffusivityC, #diffusivityD');
const diffusivityLabels = document.querySelectorAll('label[for="diffusivityA"], label[for="diffusivityB"], label[for="diffusivityC"], label[for="diffusivityD"]');
const precomputedForm = document.getElementById('precomputed-props');
const wilkeChangForm = document.getElementById('wilke-chang-props');

// Remove "required" for fields in precomputed-props inputs
function removePrecomputedRequirements() {
    diffusivityInputs.forEach(input => {
        input.removeAttribute('required');
        input.setAttribute('disabled', true);
        input.value = '';
    });
    diffusivityLabels.forEach(label => {
        label.classList.remove('required');
    });
}

function addPrecomputedRequirements() {
    // Add everything back 
    diffusivityInputs.forEach(input => {
        input.setAttribute('required', true);
        input.removeAttribute('disabled');
    });
    diffusivityLabels.forEach(label => {
        label.classList.add('required');
    });

    // Then remove / disable according to reaction type
    const selectedRxType = document.querySelector('input[name="rxType"]:checked').value;

    if (['type2', 'type4', 'type6'].includes(selectedRxType)) {
        const inputB = document.getElementById('diffusivityB');
        const labelB = document.querySelector('label[for="diffusivityB"]');
        inputB.removeAttribute('required');
        inputB.setAttribute('disabled', true);
        labelB.classList.remove('required');
    }

    if (['type3', 'type5', 'type6'].includes(selectedRxType)) {
        const inputD = document.getElementById('diffusivityD');
        const labelD = document.querySelector('label[for="diffusivityD"]');
        inputD.removeAttribute('required');
        inputD.setAttribute('disabled', true);
        labelD.classList.remove('required');
    }
}


function updateDiffusivityRequirements() {
    const selected = document.querySelector('input[name="diffusivityType"]:checked').value;

    // Make each section disappear or appear according to radio button selection
    if (selected === "precomputed") {
        precomputedForm.style.display = 'block';
        wilkeChangForm.style.display = 'none';
        addPrecomputedRequirements();
    } else if (selected === "ideal" || selected === "nonideal") {
        precomputedForm.style.display = 'none';
        wilkeChangForm.style.display = 'block';
        removePrecomputedRequirements();
    } else {
        precomputedForm.style.display = 'none';
        wilkeChangForm.style.display = 'none';
        removePrecomputedRequirements();
    }
}

const advancedTabButton = document.querySelector('.tab-button[data-tab="advanced"]');
// Utility function to make the "advanced" tab appear/disappear according to diffusivity
function toggleAdvancedTab() {
    if (advancedEnabled) {
        advancedTabButton.style.display = 'none';
    } else {
        advancedTabButton.style.display = 'block';
    }
    advancedEnabled = !advancedEnabled;
}

// Ties the toggle to diffusivity selection
// advancedEnabled starts as "false" since toggle runs at least once on load
// once "nonideal" is selected, advancedEnabled becomes true and tab appears
function updateAdvancedTabVisibility() {
    const selected = document.querySelector('input[name="diffusivityType"]:checked').value;
    if (!advancedEnabled && selected === "nonideal") {
        toggleAdvancedTab();
    } else if (advancedEnabled && selected !== "nonideal") {
        toggleAdvancedTab();
    }
}


// Event listeners for diffusivity radio
// Form requirements/form section display
diffusivityTypeRadioButtons.forEach(radio => {
    radio.addEventListener('change', updateDiffusivityRequirements);
});

// Advanced tab/next and prev button step size control
diffusivityTypeRadioButtons.forEach(radio => {
    radio.addEventListener('change', updateAdvancedTabVisibility);
});

// Also need to add an event listener for the reaction selection 
rxTypesRadioButtons.forEach(radio => {
    radio.addEventListener('change', updateDiffusivityRequirements);
});

//!SECTION



// SECTION: Collapsible Sections
const collapsibleHeaders = document.querySelectorAll('.collapsible-header');

function setupCollapsibleHeader(header) {
    header.addEventListener('click', () => {
        const content = header.nextElementSibling;
        header.classList.toggle('active');
        content.classList.toggle('active');
    });
}

collapsibleHeaders.forEach(setupCollapsibleHeader);

// !SECTION

// SECTION: Wilke-Chang parameters table
// Hide rows according to selected reaction type 
const wilkeChangTableBody = document.getElementById('wilkeChangTable');

function updateWilkeChangRowVisibility() {
    const selectedRxType = document.querySelector('input[name="rxType"]:checked').value;
    const rows = wilkeChangTableBody.querySelectorAll('tr');

    rows.forEach(row => {
        const componentName = row.cells[0].textContent.trim();
        if ((componentName === 'B' && ['type2', 'type4', 'type6'].includes(selectedRxType)) ||
            (componentName === 'D' && ['type3', 'type5', 'type6'].includes(selectedRxType))) {
            row.style.display = 'none';
        } else {
            row.style.display = '';
        }
    });
}

rxTypesRadioButtons.forEach(radio => {
    radio.addEventListener('change', updateWilkeChangRowVisibility);
});

//!SECTION

//SECTION: Group contribution table with dynamic rows
const groupMatrixBody = document.getElementById('groupMatrix');
const groupParametersBody = document.getElementById('groupParameters');
const energyMatrixHead = document.getElementById('energyMatrixHead');
const energyMatrixBody = document.getElementById('energyMatrix');
const addRowBtn = document.getElementById('addRowBtn');
let groupIdCounter = 1; // Unique ID generator for groups

function addParameterRow() {
    const groupId = groupIdCounter++;
    const rowIndex = groupMatrixBody.children.length;

    const newMatrixRow = document.createElement('tr');
    newMatrixRow.dataset.groupId = groupId;
    newMatrixRow.innerHTML = `
        <td><input type="text" placeholder="Group ${rowIndex + 1}" class="group-name"></td>
        <td><input type="number" placeholder="0" step="1" min="0" name="group${rowIndex + 1}-A"></td>
        <td><input type="number" placeholder="0" step="1" min="0" name="group${rowIndex + 1}-B"></td>
        <td><input type="number" placeholder="0" step="1" min="0" name="group${rowIndex + 1}-C"></td>
        <td><input type="number" placeholder="0" step="1" min="0" name="group${rowIndex + 1}-D"></td>
        <td class="row-action"><button type="button" class="btn-danger btn-small delete-row">✕</button></td>
    `;

    
    const newParameterRow = document.createElement('tr');
    newParameterRow.dataset.groupId = groupId;
    newParameterRow.innerHTML = `
        <td><input type="text" placeholder="Group ${rowIndex + 1}" class="group-name"></td>
        <td><input type="number" placeholder="0.0" step="0.1" name="group${rowIndex + 1}-Rk"></td>
        <td><input type="number" placeholder="0.0" step="0.1" name="group${rowIndex + 1}-Qk"></td>
        `;
    

    // Energy matrix
    // Header column
    const newEnergyColumn = document.createElement('th');
    newEnergyColumn.dataset.groupId = groupId;
    newEnergyColumn.textContent = `Group ${rowIndex + 1}`;
    // New basic row  
    const newEnergyRow = document.createElement('tr');
    newEnergyRow.dataset.groupId = groupId;
    let rowCells = `<td><input type="text" placeholder="Group ${rowIndex + 1}" class="group-name"></td>` + generateEnergyMatrixCells(0, rowIndex + 1);
    newEnergyRow.innerHTML = rowCells;

    const deleteBtn = newMatrixRow.querySelector('.delete-row');
    deleteBtn.addEventListener('click', (e) => {
        e.preventDefault();
        deletePairedRows(groupId);
    });
    
    groupMatrixBody.appendChild(newMatrixRow);
    groupParametersBody.appendChild(newParameterRow);
    energyMatrixHead.appendChild(newEnergyColumn);
    energyMatrixBody.appendChild(newEnergyRow);
    updateEnergyMatrixRows(rowIndex + 1, groupId);
    updateDeleteButtons();
}

function updateEnergyMatrixRows(numGroups, newGroupId) {
    // Loop over old rows and add new cells for the new group column
    energyMatrixBody.querySelectorAll('tr').forEach((row, rowIndex) => {
        // Do not clear the contents of any pre-existing cells
        const rowGroupId = row.dataset.groupId; 
        if (rowGroupId !== newGroupId && rowIndex !== numGroups - 1) {
            let newCells = `<td><input type="number" placeholder="0.0" step="0.1" name="energy-${numGroups}-${numGroups}"></td>`;
            row.insertAdjacentHTML('beforeend', newCells);
        }
    });
}

function generateEnergyMatrixCells(start, numGroups) {
    let newCells = ``;
    for (let j = start; j < numGroups; j++) {
        newCells += `<td><input type="number" placeholder="0.0" step="0.1" name="energy-${numGroups}-${j + 1}"></td>`;
    }
    return newCells;
}

function deletePairedRows(groupId) {
    const matrixRows = groupMatrixBody.querySelector(`tr[data-group-id="${groupId}"]`);
    const parameterRows = groupParametersBody.querySelector(`tr[data-group-id="${groupId}"]`);
    const energyRow = energyMatrixBody.querySelector(`tr[data-group-id="${groupId}"]`);
    const energyColumn = energyMatrixHead.querySelector(`th[data-group-id="${groupId}"]`);
    const allEnergyColumns = Array.from(energyMatrixHead.querySelectorAll('th'));
    const deletedColumnIndex = allEnergyColumns.findIndex(th => th.dataset.groupId === String(groupId));

    // Remove last added energy column and row
    if (energyColumn) energyColumn.remove();
    if (energyRow) energyRow.remove();
    // Remove cell at the column index from all remaining rows
    if (deletedColumnIndex !== -1) {
        energyMatrixBody.querySelectorAll('tr').forEach(row => {
            const cells = Array.from(row.querySelectorAll('td'));
            if (cells[deletedColumnIndex]) {
                cells[deletedColumnIndex].remove();
            }
        });
    }

    // Finally remove the other simpler paired rows
    if (matrixRows) matrixRows.remove();
    if (parameterRows) parameterRows.remove();

    updateDeleteButtons();
}

function updateDeleteButtons() {
    const rows = groupMatrixBody.querySelectorAll('tr');
    rows.forEach(row => {
        const deleteBtn = row.querySelector('.delete-row');
        deleteBtn.style.display = rows.length > 1 ? 'block' : 'none';
    });
}

addRowBtn.addEventListener('click', (e) => {
    e.preventDefault();
    addParameterRow();
});

// Delete row buttons for initial rows
groupMatrixBody.querySelectorAll('tr').forEach(row => {
    const deleteBtn = row.querySelector('.delete-row');
    deleteBtn.addEventListener('click', (e) => {
        e.preventDefault();
        const groupId = row.dataset.groupId;
        deletePairedRows(groupId);
    });
});

updateDeleteButtons();


const groupMatrixHead = document.querySelector('#groupMatrixHead');

function updateGroupMatrixColumnVisibility() {
    const selectedRxType = document.querySelector('input[name="rxType"]:checked').value;
    const hideB = ['type2', 'type4', 'type6'].includes(selectedRxType);
    const hideD = ['type3', 'type5', 'type6'].includes(selectedRxType);
    
    // Update header visibility
    groupMatrixHead.querySelectorAll('th').forEach((th, index) => {
        if (index === 2) { // Component B
            th.style.display = hideB ? 'none' : '';
        } else if (index === 4) { // Component D
            th.style.display = hideD ? 'none' : '';
        }
    });

    // Update body visibility of certain columns 
    groupMatrixBody.querySelectorAll('tr').forEach(row => {
        row.querySelectorAll('td').forEach((td, index) => {
            if (index === 2) { // Component B
                td.style.display = hideB ? 'none' : '';
            } else if (index === 4) { // Component D
                td.style.display = hideD ? 'none' : '';
            }
        });
    });
}

updateGroupMatrixColumnVisibility();

rxTypesRadioButtons.forEach(radio => {
    radio.addEventListener('change', updateGroupMatrixColumnVisibility);
});

//!SECTION

// SECTION: Form Validation
function validateTab(tabIndex) {
    const tabId = formState.tabs[tabIndex];
    const tabContent = document.getElementById(tabId);
    const requiredFields = tabContent.querySelectorAll('[required]');
    
    let isValid = true;
    requiredFields.forEach(field => {
        if (!field.value.trim()) {
            field.classList.add('error');
            isValid = false;
        } else {
            field.classList.remove('error');
        }
    });

    return isValid;
}

// Validation Messages
function showValidationMessage(message, type = 'info') {
    const container = document.getElementById('validationContainer');
    const messageEl = document.createElement('div');
    messageEl.className = `validation-message ${type} fade-in`;
    messageEl.innerHTML = `
        <span>${message}</span>
    `;
    container.innerHTML = '';
    container.appendChild(messageEl);

    setTimeout(() => {
        messageEl.remove();
    }, 5000);
}
// !SECTION

// SECTION: Update Review Tab
function updateReviewTab() {
    if (formState.currentTab === formState.tabs.length - 1) {
        const reviewContent = document.getElementById('reviewContent');
        const formData = new FormData(document.getElementById('calculatorForm'));
        
        let reviewHTML = '<div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: var(--spacing-lg);">';
        
        for (let [key, value] of formData.entries()) {
            if (value) {
                reviewHTML += `
                    <div style="padding: var(--spacing-md); background: var(--secondary-color); border-radius: var(--radius-md);">
                        <strong style="color: var(--primary-color);">${key}</strong><br>
                        <span style="color: var(--text-secondary);">${value}</span>
                    </div>
                `;
            }
        }
        
        reviewHTML += '</div>';
        reviewContent.innerHTML = reviewHTML;
    }
}
//!SECTION


// SECTION: Form Submission to Backend
const calculatorForm = document.getElementById('calculatorForm');

calculatorForm.addEventListener('submit', (e) => {
    e.preventDefault();
    
    // Validate all tabs
    let allValid = true;
    for (let i = 0; i < formState.tabs.length - 1; i++) {
        if (!validateTab(i)) {
            allValid = false;
            break;
        }
    }

    if (!allValid) {
        showValidationMessage('Please fill in all required fields', 'error');
        return;
    }

    const formData = new FormData(calculatorForm);
    const data = Object.fromEntries(formData);
    
    console.log('Form submitted with data:', data);
    showValidationMessage('Form submitted successfully! Check console for data.', 'success');
    
    // Send the data to backend
    fetch('/calculate', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
        },
        body: JSON.stringify(data)
    });
});

//!SECTION

// Initialize
updateNavigationButtons();
toggleAdvancedTab();
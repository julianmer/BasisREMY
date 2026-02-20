<div align="center">
  <img src="assets/imgs/basisremy_logo_round.png" alt="BasisREMY Logo" width="120" style="margin-bottom: -10px;"/>
  <h1 style="margin-top: 5px; margin-bottom: 5px;">BasisREMY</h1>
  <p style="margin-top: 0px;"><em>Study-Specific Basis Set Generation for MR Spectroscopy</em></p>
  
  [![Python](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
  [![Tests](https://img.shields.io/badge/tests-passing-brightgreen.svg)]()
  [![ISMRM 2026](https://img.shields.io/badge/ISMRM-Abstract%20%2301716-lightgrey.svg)](https://submissions.mirasmart.com/ISMRM2026/)
</div>

---

A tool for generating study-specific basis sets directly from raw MRS data, integrating real pulse shapes and acquisition parameters. This project is in its early development stages, and contributions, testing, and feedback are highly welcomed!

<div align="center">
  <img src="assets/imgs/basisremy_workflow.png" alt="BasisREMY Workflow" width="750"/>
</div>

---

## Prerequisites

Before installing BasisREMY, ensure that the following are installed on your system:

- **Python**: Version 3.10 or higher. You can download the latest version from the [official Python website](https://www.python.org/downloads/).

- **Octave Runtime** (required for simulation backends using MATLAB/Octave code): Version 4.0 or higher.

  **You have two options:**
  
  1. **Docker** (Recommended) - Automatic setup, works everywhere
  2. **Local Octave** - Traditional installation
  
  **📖 See the [Octave Setup Guide](assets/OCTAVE_SETUP.md) for detailed installation instructions.**
  
  > **Note**: BasisREMY automatically detects and uses Docker if available, otherwise falls back to local Octave.

---

## Setting Up the Python Environment

It is recommended to use a virtual environment to isolate project dependencies. Start by cloning the repository.
```bash
git clone --recurse-submodules https://github.com/julianmer/BasisREMY.git
cd BasisREMY
```

### Create and Activate the Virtual Environment

**Windows:**
```bash
python -m venv --prompt basisREMY .venv
.venv\Scripts\activate
```

**macOS/Linux:**
```bash
python3 -m venv --prompt basisREMY .venv
source .venv/bin/activate
```

### Install Required Python Packages
Upgrade pip and install the required packages by running:

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

---

## Running BasisREMY
With all dependencies installed and your basisREMY virtual environment activated, and run the application:
```bash
python main.py
```
This will launch the BasisREMY GUI.

### Usage Overview
1. Data Selection (Tab 1):
   * Drag and drop your MRS data file or click to select a file. 
   * Click Process File to extract REMY data or Skip to proceed without processing.
2. Parameter Configuration (Tab 2):
   * Adjust parameters and select metabolites using the provided checkboxes. 
   * Click Simulate Basis Set when all required parameters are provided.
3. Basis Simulation (Tab 3):
   * A progress bar will display the simulation status.
   * Once complete, the basis set is created and an interactive plot is shown.

### Examples (No GUI)

Want to use BasisREMY programmatically? Check out the **[examples/](examples/)** folder!

**Quick start:**
```bash
python examples/basic_usage.py
```

The example shows how to:
- Load MRS data and extract parameters automatically
- Configure and run simulations without the GUI
- Customize metabolite lists and output settings

---

## Related References
The project will build upon the methodologies used in existing tools. Some references include:
- [REMY](https://github.com/agudmundson/mrs_in_mrs) and related literature ([nbm.70039](https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/nbm.70039))
- [FID-A](https://github.com/CIC-methods/FID-A) and related literature ([mrm.26091](https://doi.org/10.1002/mrm.26091))
- [FSL-MRS](https://github.com/wtclarke/fsl_mrs) and related literature ([mrm.28630](https://doi.org/10.1002/mrm.28630))
- [MRSCloud](https://github.com/shui5/MRSCloud) and related literature ([mrm.29370](https://doi.org/10.1002/mrm.29370))
- [BasisSetSimulation](https://github.com/arcj-hub/BasisSetSimulation/tree/main)

---

<div align="center">
  <sub>Built with ❤️ for the MRS community</sub>
</div>
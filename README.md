# BasisREMY

A tool for generating study-specific basis sets directly from raw MRS data, integrating real pulse shapes and acquisition parameters. This project is in its early development stages, and contributions, testing, and feedback are highly welcomed!


## Prerequisites

Before installing BasisREMY, ensure that the following are installed on your system:

- **Python**: Version 3.6 or higher. You can download the latest version from the [official Python website](https://www.python.org/downloads/).

- **GNU Octave**: Version 4.0 or higher. Octave is required for BasisREMY to function properly.

  - **Windows**: Download the latest MinGW version of Octave from the [GNU Octave website](https://www.gnu.org/software/octave/download.html). During installation, ensure that Octave is added to your system's PATH. You can do this by selecting the option during installation or manually adding the path to the Octave `bin` directory to your system's PATH environment variable.

  - **macOS**: Install Octave using Homebrew:

    ```bash
    brew install octave
    ```

  - **Linux**: Install Octave using your distribution's package manager. For example, on Ubuntu:

    ```bash
    sudo apt-get update
    sudo apt-get install octave
    ```

  After installation, verify that Octave is accessible from the command line by running:

  ```bash
  octave --version
    ```


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
pip install -r requirements.txt
```

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


## Related References
The project will build upon the methodologies used in existing tools. Some references include:
- [FID-A](https://github.com/CIC-methods/FID-A) and related literature ([mrm.26091](https://doi.org/10.1002/mrm.26091))
- [REMY](https://github.com/agudmundson/mrs_in_mrs) and related literature ([arXiv:2403.19594](https://arxiv.org/abs/2403.19594))
- [BasisSetSimulation](https://github.com/arcj-hub/BasisSetSimulation/tree/main)

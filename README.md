<div align="center">
  <img src="https://raw.githubusercontent.com/julianmer/BasisREMY/main/basisremy/assets/imgs/basisremy_logo_round.png" alt="BasisREMY Logo" width="120" style="margin-bottom: -10px;"/>
  <h1 style="margin-top: 5px; margin-bottom: 5px;">BasisREMY</h1>
  <p style="margin-top: 0px;"><em>A Unified Framework for Study-Specific Basis Set Generation in MR Spectroscopy</em></p>
  
  [![PyPI](https://img.shields.io/pypi/v/basisremy.svg)](https://pypi.org/project/basisremy/)
  [![Python](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/)
  [![Tests](https://github.com/julianmer/BasisREMY/actions/workflows/tests.yml/badge.svg)](https://github.com/julianmer/BasisREMY/actions/workflows/tests.yml)
  [![ISMRM 2026](https://img.shields.io/badge/ISMRM-Abstract%20%2301716-lightgrey.svg)](https://archive.ismrm.org)
</div>

---

A tool for generating study-specific basis sets directly from raw MRS data, integrating real pulse shapes and acquisition parameters. This project is in its early development stages, and contributions, testing, and feedback are highly welcomed!

<div align="center">
  <img src="https://raw.githubusercontent.com/julianmer/BasisREMY/main/basisremy/assets/imgs/basisremy_workflow_v2.png" alt="BasisREMY Workflow" width="750"/>
</div>

---

## Prerequisites

- **[`uv`](https://docs.astral.sh/uv/)** — installs Python and all dependencies into an
  isolated environment (one-line installer below); nothing is installed globally.
- **[`git`](https://git-scm.com/)** — fetches the simulation toolboxes the first time you
  simulate (and clones the repository if you develop).
- **Docker or Octave** *(only for simulating)* — [Docker Desktop](https://www.docker.com/products/docker-desktop/)
  (recommended) or a local [Octave](https://octave.org/) 4.0+; see
  [Docker / Octave Requirements](#docker--octave-requirements).

---

## Quick Start (for users)

The fastest way to run BasisREMY is with [`uv`](https://docs.astral.sh/uv/) — a single
tool that installs Python for you (no manual Python, `pip`, or conda setup required).

**1. Install `uv`** (one time):

```bash
# macOS / Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# Windows (PowerShell)
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

**2. Launch BasisREMY** — the released version, straight from PyPI:

```bash
uvx basisremy
```

Or, for repeated use, install the command once and then run it by name:

```bash
uv tool install basisremy
basisremy
```

To run the latest development state instead of the release:

```bash
uvx --from git+https://github.com/julianmer/BasisREMY basisremy
```

> **Using BasisREMY from Python instead of the GUI?** The package is on PyPI —
> `pip install basisremy` gives you the programmatic API
> (see [Examples](#examples-no-gui) below) without any GUI setup.

`uv` automatically downloads a compatible Python and all dependencies into an
isolated environment — nothing is installed globally and your data never leaves
your machine.

> **First run:** the first simulation with an engine on a new machine sets that engine up
> (its toolbox at a pinned version into `~/.basisremy` or `$BASISREMY_HOME`, and its Docker
> image if needed) and can take a few minutes; after that it starts in seconds.
> **Ran BasisREMY before?** `uvx basisremy@latest` gets the newest release.

---

## Usage Overview

BasisREMY guides you through a simple **three-step workflow**, shown as a numbered
stepper at the top of the window (**Data → Parameters → Simulate**). You can click
any completed step to go back and make changes. A light/dark theme toggle sits in
the top-right corner and follows your system setting by default.

1. **Data** — Drop your MRS data file onto the upload area (or click it to browse).
   Click **Continue** to automatically read the file and pre-fill the acquisition
   parameters, or **Skip** to configure everything manually.
2. **Parameters** — One layout for every engine: pick the **Engine** and the
   **Sequence**, optionally a **Sequence file** (any pulse or whole-sequence file),
   then check the **timings**, the **acquisition** values and the **metabolites**
   to include. Each value is marked by where it comes from: *from the data file*,
   *recommended / default* (for design values no header holds, such as the echo
   split; hover for the source), *set by you*, or *missing* (Simulate stays
   blocked until you fill it). When the required fields are set, click
   **Simulate basis set**.
3. **Simulate** — A progress bar shows the simulation status. When it finishes,
   an interactive spectrum plot appears and you can click **Export basis…** to
   save the basis set in LCModel (.basis or .RAW), jMRUI, FSL-MRS, Osprey, FID-A,
   INSPECTOR, ProFit, MARSS, MRSCloud, or SpinWizard/JET format. Format
   writing is handled by the bundled
   [MRS Basis Set Conversion Toolbox](https://github.com/igweckay/MRS-Basis-Set-Conversion-Toolbox).

> Every input has a small **(?)** help icon — hover it for a short explanation of
> that parameter.

**Sequence designer** *(new, in testing)* — the wand next to the sequence file
opens the designer: choose the sequence, its timings and a pulse per role (ideal,
a standard shape, or a pulse file in any format; the scan's own pulse where the data
header stores it), see the timeline, and save the design as a Pulseq `.seq`
(in `~/BasisREMY/sequences/`). The saved file is selected right away; the engine
chips show which engines run the design as is, with differences, or not at all.


## Examples (No GUI)

Want to use BasisREMY programmatically? Check out the **[examples/](examples/)** folder!

**Quick start:**
```bash
python examples/basic_usage.py
```

The example shows how to:
- Load MRS data and extract parameters automatically
- Configure and run simulations without the GUI
- Customize metabolite lists and output settings

In your own scripts, call `basisremy.prepare_runtime()` first and give file paths as absolute
paths: the simulation engines run from BasisREMY's runtime folder, as in the GUI.

---

## Developer Setup

Clone the repository (with submodules) and let `uv` create the environment:

```bash
git clone --recurse-submodules https://github.com/julianmer/BasisREMY.git
cd BasisREMY
uv sync            # creates .venv and installs all dependencies
uv run basisremy   # launches the GUI
```

`uv sync` reads [`pyproject.toml`](pyproject.toml), provisions a matching Python,
and installs everything into `.venv`. You do **not** need to activate the
environment — `uv run` handles that automatically. To include development tools
(pytest, coverage), run `uv sync --extra dev`.

Prefer pip? `pip install basisremy` for the Python API, or `pip install -e .` in a clone.

---

## Docker / Octave Requirements

Simulation backends execute MATLAB/Octave code and therefore need an Octave
runtime. BasisREMY supports **two interchangeable options** and selects one
automatically:

1. **Docker** (recommended) — install [Docker Desktop](https://www.docker.com/products/docker-desktop/)
   (macOS/Windows) or Docker Engine (Linux), and make sure it is **running**.
   BasisREMY pulls/builds an Octave container on demand.
2. **Local Octave** — install [Octave](https://octave.org/) 4.0+ directly; it is
   used as a fallback when Docker is unavailable.

Whether an Octave runtime is needed is **backend-dependent**: data extraction and
parameter configuration work without it; only the simulation step requires it.
See the [Octave Setup Guide](basisremy/assets/OCTAVE_SETUP.md) for detailed instructions.

The **Spant** backend needs **R** in the same way: Docker (BasisREMY builds an R
image with the `spant` package on first use) or an R installation of your own,
into which the package is installed on first use. The Vespa backend needs
neither — it runs in a Python side-environment or Docker.

Three environment variables tune the runtimes:

| Variable | Effect |
| --- | --- |
| `BASISREMY_OCTAVE_WORKERS` | Octave processes used to simulate metabolites in parallel (default: half the CPUs, at most 4). |
| `BASISREMY_NO_DOCKER` | Set to `1` to ignore Docker entirely and use only a local Octave. |
| `BASISREMY_SPANT_RUNTIME` | `docker` or `local` to force the Spant runtime (default: Docker, then your own R). |

---

## Troubleshooting

- **Window stays blank or pywebview errors** — without `pywebview` the UI falls
  back to opening in your browser; install it (`pip install pywebview`) for the
  native window. On Linux a system WebKitGTK package may be needed.
- **Simulation fails / "Octave not found"** — start Docker Desktop/Engine, or
  install local Octave. See the [Octave Setup Guide](basisremy/assets/OCTAVE_SETUP.md).
- **`uv: command not found`** — re-open your terminal after installing `uv`, or
  add its install location to your `PATH`.
- **Backends can't find `externals/...`** — the simulation toolboxes are fetched
  on first use into `~/.basisremy` (or `$BASISREMY_HOME`); make sure `git` is
  installed and you have network access the first time you run a simulation. In a
  cloned repository the existing submodules under `externals/` are used as-is.

---

## Related References
The project builds upon the methodologies used in existing tools. Some references include:
- [REMY](https://github.com/agudmundson/mrs_in_mrs) and related literature ([nbm.70039](https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/nbm.70039))
- [FID-A](https://github.com/CIC-methods/FID-A) and related literature ([mrm.26091](https://doi.org/10.1002/mrm.26091))
- [FSL-MRS](https://github.com/wtclarke/fsl_mrs) and related literature ([mrm.28630](https://doi.org/10.1002/mrm.28630))
- [MRSCloud](https://github.com/shui5/MRSCloud) and related literature ([mrm.29370](https://doi.org/10.1002/mrm.29370))
- [Vespa](https://github.com/vespa-mrs/vespa) and related literature ([mrm.29686](https://doi.org/10.1002/mrm.29686))
- [Spant](https://github.com/martin3141/spant) and related literature ([mrm.28385](https://doi.org/10.1002/mrm.28385))
- [Spinach](https://github.com/IlyaKuprov/Spinach) and related literature ([jmr.2010.11.008](https://doi.org/10.1016/j.jmr.2010.11.008))
- [MRS Basis Set Conversion Toolbox](https://github.com/igweckay/MRS-Basis-Set-Conversion-Toolbox)
- [Custom Basis Set Simulation](https://github.com/arcj-hub/BasisSetSimulation/tree/main)

---

<div align="center">
  <sub>Built with ❤️ for the MRS community</sub>
</div>
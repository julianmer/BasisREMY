####################################################################################################
#                                       parameter_registry.py                                      #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Purpose: Centralized registry mapping every parameter name used in any backend to a              #
#          human-readable description, units, typical range, and aliases. Used by:                 #
#            - gui/help_widget.py to populate "?" tooltips next to every parameter                 #
#            - core/exporters.py to write reproducibility sidecars                                 #
#                                                                                                  #
# TODO: Check and improve the descriptions for all parameters.                                     #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any


# Sentinel marking parameters that still need a human-written description
TODO = "TODO_PLACEHOLDER"


@dataclass
class ParamInfo:
    label: str                      # short display name
    description: str                # tooltip body
    units: str = ""                 # e.g. "Hz", "ms", "T", "ppm"
    typical: str = ""               # e.g. "1024 - 8192", "30 - 288 ms"
    aliases: tuple = field(default_factory=tuple)
    widget_hint: str = "entry"      # one of {"entry", "combobox", "checkbox", "file", "directory", "metabolites"}


# ---- Registry ---------------------------------------------------------------
# Keyed by canonical parameter name (matches keys in backends' mandatory_params /
# optional_params / dropdown / file_selection). Pre-filled with descriptions
# for parameters whose meaning is unambiguous from the MRS literature; ambiguous
# / backend-specific ones use TODO_PLACEHOLDER and should be filled by the user.

REGISTRY: dict[str, ParamInfo] = {
    # --- Acquisition ---------------------------------------------------------
    "Bandwidth": ParamInfo(
        label="Bandwidth",
        description=(
            "Spectral width of the acquired FID, i.e. the receiver sampling rate. "
            "Determines the frequency range covered by the spectrum (Nyquist limits "
            "to ±Bandwidth/2 around the carrier). For 1H MRS at 3T, 2000-4000 Hz is "
            "typical; 7T scans use 4000-6000 Hz."
        ),
        units="Hz",
        typical="2000 - 6000 Hz",
        aliases=("SpectralWidth", "sw", "Rx_SW"),
    ),
    "Samples": ParamInfo(
        label="Samples",
        description=(
            "Number of complex points in the acquired FID. Together with the "
            "bandwidth, sets the spectral resolution: Δf = Bandwidth / Samples. "
            "More points yield finer line shape but longer acquisition / simulation."
        ),
        units="points",
        typical="1024 - 8192",
        aliases=("NumberOfDatapoints", "n", "Rx_Points", "Npts"),
    ),
    "Bfield": ParamInfo(
        label="B0 field strength",
        description=(
            "Static magnetic field of the scanner. Determines the proton Larmor "
            "frequency (ν₀ ≈ 42.577 · B0 MHz) and indirectly chemical-shift dispersion."
        ),
        units="T",
        typical="1.5 / 3 / 7 T",
        aliases=("B0",),
    ),
    "Center Freq": ParamInfo(
        label="Centre frequency",
        description=(
            "Spectrometer / synthesizer carrier frequency. For 1H this is "
            "≈ 42.577 MHz × B0 (≈ 127.7 MHz at 3 T, ≈ 297.2 MHz at 7 T). "
            "Used for ppm ↔ Hz conversion."
        ),
        units="MHz",
        typical="63 / 127 / 297 MHz",
        aliases=("centreFreq", "SpectrometerFrequency", "lFrequency", "MRFrequency"),
    ),
    "Nucleus": ParamInfo(
        label="Nucleus",
        description="Observed NMR-active nucleus. Most clinical MRS uses 1H; 31P, 13C, 19F also supported.",
        typical="1H",
    ),
    "Linewidth": ParamInfo(
        label="Linewidth",
        description=(
            "Lorentzian full-width at half-maximum (FWHM) applied to the simulated "
            "FIDs as exponential apodization. Sets the broadening of basis spectra."
        ),
        units="Hz",
        typical="1 - 5 Hz",
        aliases=("Rx_LW", "lb", "fwhm"),
    ),

    # --- Sequence timing -----------------------------------------------------
    "TE": ParamInfo(
        label="Echo time (TE)",
        description="Total echo time of the sequence (excitation to acquisition centre).",
        units="ms",
        typical="20 - 288 ms",
    ),
    "TE2": ParamInfo(
        label="Second sub-TE",
        description=(
            "Second TE component for multi-echo sequences (e.g. PRESS uses "
            "TE = TE1 + TE2). Set to 0 for spin-echo / STEAM."
        ),
        units="ms",
        typical="0 - TE/2",
    ),
    "TM": ParamInfo(
        label="Mixing time (TM)",
        description="STEAM mixing time between the 2nd and 3rd 90° pulses.",
        units="ms",
        typical="10 - 50 ms",
    ),
    "TR": ParamInfo(
        label="Repetition time (TR)",
        description="Time between successive excitations. Affects T1-weighting; not used for basis simulation directly.",
        units="ms",
        typical="1500 - 4000 ms",
    ),
    "Tau 1": ParamInfo(
        label="Tau 1",
        description=(
            TODO + " — preliminary: first sLASER inter-pulse delay, i.e. the "
            "time between the excitation pulse and the centre of the first AFP "
            "refocusing pair. Total echo time TE = 2·(Tau1 + Tau2). Please "
            "confirm the exact convention used by the active backend."
        ),
        units="ms",
        typical="6 - 15 ms",
    ),
    "Tau 2": ParamInfo(
        label="Tau 2",
        description=(
            TODO + " — preliminary: second sLASER inter-pulse delay between the "
            "centres of the first and second AFP refocusing pairs. Total echo "
            "time TE = 2·(Tau1 + Tau2). Please confirm the exact convention "
            "used by the active backend."
        ),
        units="ms",
        typical="6 - 15 ms",
    ),

    # --- Sequence / mode -----------------------------------------------------
    "Sequence": ParamInfo(
        label="Sequence",
        description=(
            "Localization sequence used for acquisition. Each backend exposes a "
            "different sub-set (PRESS, STEAM, sLASER, LASER, MEGA-PRESS, HERMES, "
            "HERCULES, …)."
        ),
        widget_hint="combobox",
    ),
    "System": ParamInfo(
        label="Vendor / system",
        description="Scanner manufacturer (Siemens, Philips, GE, Bruker). Selects vendor-specific pulse files in some backends.",
        widget_hint="combobox",
    ),
    "Edit Frequency": ParamInfo(
        label="Editing frequency",
        description=(
            "Chemical-shift offset (ppm) of the editing pulse for J-difference "
            "edited MRS (e.g. 1.9 ppm targets GABA, 4.56 ppm targets GSH)."
        ),
        units="ppm",
        typical="1.9 (GABA), 4.56 (GSH)",
    ),

    # --- Pulse / RF ---------------------------------------------------------
    "B1max": ParamInfo(
        label="B1 max",
        description=(
            TODO + " — preliminary: peak RF amplitude of the refocusing (AFP) "
            "pulse, used by FID-A / CustomSLaser to scale the simulated waveform. "
            "FID-A typically expects Hz (γ·B1/2π). Please confirm the unit "
            "convention of the backend you are using; convert via "
            "B1[Hz] = 42.577 · B1[µT] for 1H."
        ),
        units="Hz (FID-A) / µT",
        typical="800 - 1500 Hz  ·  20 - 25 µT",
    ),
    "Flip Angle": ParamInfo(
        label="Flip angle",
        description="Nominal flip angle of the refocusing pulse (sLASER backend uses 180° AFP).",
        units="degrees",
        typical="180°",
        aliases=("ExcitationFlipAngle",),
    ),
    "RefTp": ParamInfo(
        label="Refocusing pulse duration",
        description="Duration of the refocusing RF pulse (used for slice-selection bandwidth in sLASER).",
        units="ms",
        typical="3 - 6 ms",
    ),
    "Path to Pulse": ParamInfo(
        label="Path to pulse waveform",
        description="File containing the refocusing RF pulse shape (vendor-specific .pta / .RF / .pulse / .json file).",
        widget_hint="file",
    ),
    "Vendor Pulse File": ParamInfo(
        label="Vendor pulse file",
        description=(
            "Vendor-confidential RF waveform that MRSCloud needs but does not "
            "ship in its public repo. Only required when a non-Universal_* "
            "vendor is selected."
        ),
        widget_hint="file",
    ),

    # --- Voxel geometry -----------------------------------------------------
    "thkX": ParamInfo(
        label="Slice thickness X",
        description="Voxel thickness along the X (read) direction for the sLASER spatial simulation grid.",
        units="cm",
        typical="1.5 - 3 cm",
    ),
    "thkY": ParamInfo(
        label="Slice thickness Y",
        description="Voxel thickness along the Y (phase) direction for the sLASER spatial simulation grid.",
        units="cm",
        typical="1.5 - 3 cm",
    ),
    "fovX": ParamInfo(
        label="FOV X",
        description="Spatial field-of-view along X for the simulation grid (defaults to thkX + 1 cm if unknown).",
        units="cm",
        typical="thkX + 1 cm",
    ),
    "fovY": ParamInfo(
        label="FOV Y",
        description="Spatial field-of-view along Y for the simulation grid (defaults to thkY + 1 cm if unknown).",
        units="cm",
        typical="thkY + 1 cm",
    ),
    "nX": ParamInfo(
        label="Spatial points X",
        description="Number of spatial sample points along X used by the spatially-resolved sLASER simulation. Higher = more accurate, slower.",
        typical="32 - 64",
    ),
    "nY": ParamInfo(
        label="Spatial points Y",
        description="Number of spatial sample points along Y used by the spatially-resolved sLASER simulation. Higher = more accurate, slower.",
        typical="32 - 64",
    ),

    # --- Output / I-O -------------------------------------------------------
    "Output Path": ParamInfo(
        label="Output directory",
        description="Folder where basis files, sequence JSON, and the reproducibility sidecar will be written.",
        widget_hint="directory",
    ),
    "Output Format": ParamInfo(
        label="Output format",
        description=(
            "File format for the generated basis. "
            "LCModel `.basis` and `.RAW`, jMRUI `.txt`, FSL-MRS `.json`/`.basis` "
            "directory, and Osprey `.mat` are supported via the unified Exporter."
        ),
        widget_hint="combobox",
    ),
    "Basis Name": ParamInfo(
        label="Basis name",
        description="Filename (no extension) for the output basis set.",
    ),
    "Add Ref.": ParamInfo(
        label="Add reference peak",
        description=(
            "If 'Yes', a 0-ppm reference singlet (e.g. TMS / DSS surrogate) is "
            "added to the basis. Useful for LCModel referencing."
        ),
        widget_hint="combobox",
    ),
    "Add Reference": ParamInfo(
        label="Add reference peak",
        description="Same as 'Add Ref.' — adds a 0-ppm reference singlet to the basis (FSL-MRS backend naming).",
        widget_hint="checkbox",
    ),
    "Make .raw": ParamInfo(
        label="Write LCModel .RAW per metabolite",
        description=(
            "Write one LCModel-compatible `.RAW` file per metabolite into the output "
            "directory (in addition to whatever the selected Output Format produces). "
            "Currently required by some FID-A scripts to run."
        ),
        widget_hint="combobox",
    ),
    "Template File": ParamInfo(
        label="Template sequence file",
        description=(
            "Pre-bundled FSL-MRS sequence JSON containing real RF pulse shapes "
            "for a specific scanner / TE. Accurate only if your acquisition "
            "parameters match the template's B0 and TE."
        ),
        widget_hint="combobox",
    ),
    "Custom Sequence": ParamInfo(
        label="Custom sequence JSON",
        description="Path to a user-supplied FSL-MRS sequence JSON (real pulses, custom timings).",
        widget_hint="file",
    ),
    "Auto Phase": ParamInfo(
        label="Auto phase",
        description=(
            TODO + " — preliminary: when enabled, applies a zero-order phase "
            "correction to each simulated metabolite FID so the on-resonance "
            "singlet is purely absorptive (real-positive) at t = 0. Please "
            "verify the exact behaviour for the FSL-MRS backend before relying "
            "on it for quantification."
        ),
        widget_hint="checkbox",
    ),
    "Parallel": ParamInfo(
        label="Parallel processing",
        description="Use multi-process parallelism when simulating metabolites (Python-only backends).",
        widget_hint="checkbox",
    ),
    "Metabolites": ParamInfo(
        label="Metabolites",
        description="Subset of spin-system files to include in the generated basis. At least one must be selected.",
        widget_hint="metabolites",
    ),

    # --- MRSCloud-specific --------------------------------------------------
    "Localization": ParamInfo(
        label="Localization",
        description=(
            "Voxel localization scheme used by MRSCloud. PRESS uses two refocusing "
            "pulses, sLASER uses adiabatic AFP refocusing, STEAM_7T is the 7-T "
            "STEAM variant bundled with MRSCloud."
        ),
        typical="PRESS / sLASER / STEAM_7T",
        widget_hint="combobox",
    ),
    "Field Strength": ParamInfo(
        label="Field strength preset",
        description=(
            "Coarse field-strength preset used by MRSCloud's load_parameters "
            "to pick the right pulse waveforms (1.5 T, 3 T, or 7 T). For the "
            "MRSCloud backend this preset is also the canonical B0 source — "
            "it is converted to a Larmor frequency internally as 42.577 × B0 MHz "
            "(no separate Bfield parameter is exposed)."
        ),
        typical="3T",
        widget_hint="combobox",
    ),
    "Spatial Points": ParamInfo(
        label="Spatial points",
        description=(
            "Number of spatial sampling points used by MRSCloud's 1-D projection "
            "method per spatial direction. 41 is acceptable, 101 is the ideal "
            "(but considerably slower) setting recommended by the MRSCloud authors."
        ),
        typical="41 (fast)  /  101 (ideal)",
        widget_hint="entry",
    ),
    "Edit Target": ParamInfo(
        label="Editing target",
        description=(
            "Metabolite targeted by the editing scheme of an MRSCloud edited "
            "sequence (MEGA / HERMES / HERCULES). Leave empty for un-edited "
            "acquisitions."
        ),
        typical="GABA / GSH / Lac / PE",
        widget_hint="combobox",
    ),
    "Edit On": ParamInfo(
        label="Edit-ON frequency",
        description=(
            "Chemical-shift offset (ppm) of the editing pulse in the ON sub-spectrum "
            "for MEGA-style J-difference editing. Defaults: 1.9 ppm (GABA), 4.56 ppm (GSH)."
        ),
        units="ppm",
        typical="1.9 - 4.56 ppm",
        widget_hint="entry",
    ),
    "Edit Off": ParamInfo(
        label="Edit-OFF frequency",
        description=(
            "Chemical-shift offset (ppm) of the editing pulse in the OFF sub-spectrum "
            "(typically far off-resonance, e.g. 7.5 ppm)."
        ),
        units="ppm",
        typical="7.5 ppm",
        widget_hint="entry",
    ),
    "Edit Tp": ParamInfo(
        label="Editing pulse duration",
        description=(
            "Duration of the editing RF pulse — 14 ms for MEGA-PRESS / HERMES, "
            "20 ms for HERCULES (per MRSCloud documentation)."
        ),
        units="ms",
        typical="14 ms (MEGA/HERMES)  ·  20 ms (HERCULES)",
        widget_hint="entry",
    ),

    # --- Misc REMY-extracted ------------------------------------------------
    "Model": ParamInfo(label="Scanner model", description="Scanner model string (REMY)."),
    "SoftwareVersion": ParamInfo(label="Software version", description="Scanner software version (REMY)."),
    "BodyPart": ParamInfo(label="Body part", description="Anatomy scanned (REMY)."),
    "VOI": ParamInfo(label="VOI", description="Volume-of-interest descriptor (REMY)."),
    "AnteriorPosteriorSize": ParamInfo(label="A-P voxel size", description="Voxel size in the anterior-posterior direction.", units="mm"),
    "LeftRightSize": ParamInfo(label="L-R voxel size", description="Voxel size in the left-right direction.", units="mm"),
    "CranioCaudalSize": ParamInfo(label="C-C voxel size", description="Voxel size in the cranio-caudal direction.", units="mm"),
    "NumberOfAverages": ParamInfo(label="Averages", description="Number of acquisitions averaged."),
    "WaterSuppression": ParamInfo(label="Water suppression", description="Water-suppression scheme reported by the scanner (REMY)."),
}


# ----------- helpers ---------------------------------------------------------

def get(param: str) -> ParamInfo:
    """Return the registry entry for `param`, or a TODO placeholder if missing."""
    if param in REGISTRY:
        return REGISTRY[param]
    return ParamInfo(
        label=param,
        description=TODO + f" (no registry entry for '{param}' — please add one in core/parameter_registry.py).",
    )


def tooltip_text(param: str) -> str:
    """Multi-line tooltip body suitable for rendering in a Tk Label."""
    info = get(param)
    parts = [info.description]
    meta = []
    if info.units:
        meta.append(f"Units: {info.units}")
    if info.typical:
        meta.append(f"Typical: {info.typical}")
    if meta:
        parts.append("")
        parts.append("  •  ".join(meta))
    return "\n".join(parts)


def missing_descriptions() -> list[str]:
    """Parameters whose description is still a TODO placeholder."""
    return [k for k, v in REGISTRY.items() if TODO in v.description]


def to_dict() -> dict[str, dict[str, Any]]:
    """Serialisable form for JSON sidecars."""
    return {k: asdict(v) for k, v in REGISTRY.items()}


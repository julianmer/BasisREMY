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
    "Sim Centre (ppm)": ParamInfo(
        label="Simulation centre (ppm)",
        description=(
            "Rotating-frame centre of the FID-A shaped-pulse simulation, in "
            "ppm — spins evolve relative to this chemical shift. 4.65 ppm "
            "(water) is standard for 1H; editing simulations sometimes use "
            "the editing target instead."
        ),
        units="ppm",
        typical="4.65 (water)",
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
        label="Second echo time (PRESS)",
        description=(
            "Second spin-echo time of the ideal PRESS (FID-A sim_press: "
            "TE = TE1 + TE2). Leave 0 / blank for a symmetric PRESS with TE/2 "
            "per echo; enter the scanner's TE2 to model an asymmetric product "
            "PRESS. Ignored by Spin Echo / LASER; STEAM uses TM instead."
        ),
        units="ms",
        typical="0 (symmetric) or the vendor's TE2",
    ),
    "TM": ParamInfo(
        label="Mixing time (TM)",
        description="STEAM mixing time between the 2nd and 3rd 90° pulses.",
        units="ms",
        typical="10 - 50 ms",
    ),
    "Edit Bandwidth (ppm)": ParamInfo(
        label="Edit bandwidth",
        description=(
            "Inversion band of the ideal editing pulse: spins within "
            "± bandwidth/2 of 'Edit On' are fully inverted, spins outside are "
            "untouched. Approximates the finite bandwidth of a real editing "
            "pulse."
        ),
        units="ppm",
        typical="0.5 - 1.5 ppm",
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
            "Echo time of the first PRESS spin echo (FID-A sim_press "
            "convention: excitation → first refocusing → first echo); "
            "TE = Tau 1 + Tau 2. PRESS shaped: leave blank for TE/2. "
            "CustomSLaser: only used for the short reference PRESS run that "
            "sets the ppm range (jbss); the sLASER timings themselves follow "
            "from TE and RefTp."
        ),
        units="ms",
        typical="TE/2 (e.g. 17.5 ms at TE 35 ms)",
    ),
    "Tau 2": ParamInfo(
        label="Tau 2",
        description=(
            "Echo time of the second PRESS spin echo (FID-A sim_press "
            "convention); TE = Tau 1 + Tau 2. PRESS shaped: leave blank for "
            "TE/2. CustomSLaser: only used for the reference PRESS run that "
            "sets the ppm range (jbss), not for the sLASER timings."
        ),
        units="ms",
        typical="TE/2 (e.g. 17.5 ms at TE 35 ms)",
    ),

    "Tau": ParamInfo(
        label="Echo time per echo",
        description=(
            "Echo time of each spin echo in the multi-echo train (FID-A "
            "sim_spinecho_xN: 'tau = echo time in ms'); the train is Nechoes "
            "refocusing pulses long."
        ),
        units="ms",
        typical="10 - 40 ms",
    ),
    "Nechoes": ParamInfo(
        label="Number of echoes",
        description="Number of spin echoes (refocusing pulses) in the multi-echo train (FID-A sim_spinecho_xN).",
        typical="2 - 10",
    ),
    "Delay": ParamInfo(
        label="ADC onset delay",
        description=(
            "Delay between the excitation pulse and the start of acquisition "
            "(FID-A sim_onepulse_delay); produces the first-order phase a real "
            "acquisition delay would."
        ),
        units="ms",
        typical="0 - 1 ms",
    ),
    "Pulse Phase": ParamInfo(
        label="Excitation pulse phase",
        description="Phase of the ideal excitation pulse (FID-A sim_onepulse_arbPh); 0 = x, 90 = y.",
        units="deg",
        typical="0 - 360",
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
    "Edit Bandwidth (Hz)": ParamInfo(
        label="Editing pulse bandwidth",
        description=(
            "Bandwidth of the Gaussian editing pulse in spant's MEGA-PRESS "
            "(seq_mega_press_ideal, default 110 Hz); the pulse is applied at "
            "'Edit On' for the ON and at 'Edit Off' for the OFF sub-spectrum."
        ),
        units="Hz",
        typical="80 - 140 Hz",
    ),
    "STEAM Variant": ParamInfo(
        label="STEAM simulation variant",
        description=(
            "How spant handles the STEAM mixing period: 'Standard' "
            "(seq_steam_ideal), 'Coherence filter' (seq_steam_ideal_cof, keeps "
            "only the zero-order coherences during TM, like FID-A and BasisREMY's "
            "Vespa STEAM) or 'z-rotation (Young)' (seq_steam_ideal_young, "
            "gradient simulation by z-rotations)."
        ),
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
            "Maximum B1 amplitude available at the scanner, in µT. "
            "CustomSLaser scales the loaded AFP refocusing waveform to it "
            "(jbss io_loadRFwaveform: w1max = γ·B1max)."
        ),
        units="µT",
        typical="22 µT (jbss default)",
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
    "Edit Pulse Path": ParamInfo(
        label="Editing pulse waveform",
        description=(
            "File with the frequency-selective editing RF waveform (.pta / "
            ".RF / .txt). Loaded as an inversion pulse and frequency-shifted "
            "to 'Edit On' for the ON and 'Edit Off' for the OFF sub-spectrum."
        ),
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
        description=(
            "Spatial field-of-view along X for the simulation grid. FID-A places "
            "the nX points at linspace(-fov/2, fov/2, nX), so the outermost "
            "points sit at the FOV edge: keep fovX ≥ thkX for a proper voxel "
            "average, and note that a 2-point grid samples only the FOV edges."
        ),
        units="cm",
        typical="thkX + 1 cm",
    ),
    "fovY": ParamInfo(
        label="FOV Y",
        description=(
            "Spatial field-of-view along Y for the simulation grid. FID-A places "
            "the nY points at linspace(-fov/2, fov/2, nY), so the outermost "
            "points sit at the FOV edge: keep fovY ≥ thkY for a proper voxel "
            "average, and note that a 2-point grid samples only the FOV edges."
        ),
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
            "Duration of the editing RF pulse. MRSCloud: MEGA only — HERMES / "
            "HERCULES use MRSCloud's fixed 20 ms pulses. FID-A shaped MEGA "
            "kinds: the duration the loaded waveform is played out over."
        ),
        units="ms",
        typical="14 ms (MEGA at TE 68)  ·  20 ms (HERCULES)",
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


# ----------- metabolite display names ---------------------------------------
# Full names for the abbreviations the backends use. Display-only — the keys
# passed to the simulation engines stay each library's own convention.
METABOLITE_NAMES: dict[str, str] = {
    # TODO (needs a human): meanings of 'Oac', 'AcO', 'Hist', 'EA' are
    # unverified in the upstream sources — deliberately left without
    # tooltips rather than guessed.
    "AcAc": "Acetoacetate",
    "Ace": "Acetate",
    "Acn": "Acetone",
    "Ala": "Alanine",
    "Asc": "Ascorbate (vitamin C)",
    "Asp": "Aspartate",
    "Bet": "Betaine",
    "Ch": "Choline",
    "Cit": "Citrate",
    "Cr": "Creatine",
    "Cystat": "Cystathionine",
    "EtOH": "Ethanol",
    "GABA": "γ-Aminobutyric acid",
    "GABA_gov": "γ-Aminobutyric acid (Govindaraju values)",
    "GABA_govind": "γ-Aminobutyric acid (Govindaraju values)",
    "GPC": "Glycerophosphocholine",
    "GSH": "Glutathione",
    "GSH_v2": "Glutathione (alternative spin system)",
    "Glc": "Glucose",
    "Gln": "Glutamine",
    "Glu": "Glutamate",
    "Gly": "Glycine",
    "Gua": "Guanidinoacetate",
    "H2O": "Water",
    "HCar": "Homocarnosine",
    "Hypotau": "Hypotaurine",
    "Ins": "myo-Inositol",
    "Lac": "Lactate",
    "Lip": "Lipids",
    "Lip09": "Lipid (0.9 ppm)",
    "Lip13a": "Lipid (1.3 ppm, a)",
    "Lip13b": "Lipid (1.3 ppm, b)",
    "Lip20": "Lipid (2.0 ppm)",
    "Lys": "Lysine",
    "MSM": "Methylsulfonylmethane",
    "NAA": "N-Acetylaspartate",
    "NAAG": "N-Acetylaspartylglutamate",
    "PCh": "Phosphocholine",
    "PCr": "Phosphocreatine",
    "PE": "Phosphorylethanolamine",
    "Phenyl": "Phenylalanine",
    "Pyr": "Pyruvate",
    "Pyruv": "Pyruvate",
    "Ref0ppm": "Reference singlet at 0 ppm",
    "Scyllo": "scyllo-Inositol",
    "Ser": "Serine",
    "Suc": "Succinate",
    "Tau": "Taurine",
    "Tau_govind": "Taurine (Govindaraju values)",
    "Thr": "Threonine",
    "Tryp": "Tryptophan",
    "Tyros": "Tyrosine",
    "Val": "Valine",
    "bHB": "β-Hydroxybutyrate",
    "bHG": "2-Hydroxyglutarate",
    "iLe": "Isoleucine",
    "mI": "myo-Inositol",
    "sI": "scyllo-Inositol",
    # spant names (get_mol_names) that differ from the other libraries
    "sins": "scyllo-Inositol",
    "peth": "Phosphorylethanolamine",
    "2hg": "2-Hydroxyglutarate",
    "a_glc": "α-Glucose",
    "b_glc": "β-Glucose",
    "mm09": "Macromolecule (0.9 ppm)",
    "mm12": "Macromolecule (1.2 ppm)",
    "mm14": "Macromolecule (1.4 ppm)",
    "mm17": "Macromolecule (1.7 ppm)",
    "mm20": "Macromolecule (2.0 ppm)",
}


_METABOLITE_NAMES_CI = {k.lower(): v for k, v in METABOLITE_NAMES.items()}


def metabolite_full_name(abbr: str) -> "str | None":
    """Full display name for a metabolite abbreviation, or None if unknown
    (case-insensitive: spant uses lower-case keys such as 'naa')."""
    name = METABOLITE_NAMES.get(str(abbr))
    if name is None:
        name = _METABOLITE_NAMES_CI.get(str(abbr).lower())
    return name


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


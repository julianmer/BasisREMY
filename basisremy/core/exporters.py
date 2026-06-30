####################################################################################################
#                                            exporters.py                                          #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Purpose: Unified basis-set exporter. Backends only need to return a                               #
#          { metabolite_name: complex_FID_array } dictionary; this module turns that               #
#          dict into any of the formats commonly used in the MRS community:                        #
#                                                                                                  #
#            - LCModel  .basis              (single combined file)                                  #
#            - LCModel  .RAW (per metab)    (folder of individual FIDs)                            #
#            - jMRUI    .txt (per metab)                                                           #
#            - FSL-MRS  .json (per metab)   (FSL-MRS basis directory)                              #
#            - Osprey   .mat                (BASIS struct compatible with Osprey)                  #
#                                                                                                  #
#          A reproducibility sidecar JSON is always written alongside the export,                  #
#          containing the BasisREMY version, git SHA (if available), backend +                     #
#          mode, and every parameter from core.parameter_registry.                                 #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import datetime as _dt
import hashlib as _hashlib
import json as _json
import os
import subprocess
from typing import Any, Iterable

import numpy as np


# ----------------------------- public surface --------------------------------

SUPPORTED_FORMATS = (
    "lcmodel_basis",   # single .basis file
    "lcmodel_raw",     # folder of <metab>.RAW files
    "jmrui_txt",       # folder of <metab>.txt files
    "fsl_json",        # folder of <metab>.json files (FSL-MRS basis dir)
    "osprey_mat",      # single .mat with Osprey-style BASIS struct
)

# Human-friendly labels for the GUI dropdown
FORMAT_LABELS = {
    "lcmodel_basis": "LCModel .basis (single file)",
    "lcmodel_raw":   "LCModel .RAW (one file per metabolite)",
    "jmrui_txt":     "jMRUI .txt (one file per metabolite)",
    "fsl_json":      "FSL-MRS .json basis directory",
    "osprey_mat":    "Osprey .mat (BASIS struct)",
}

FORMAT_EXTENSIONS = {
    "lcmodel_basis": ".basis",
    "lcmodel_raw":   "",            # directory
    "jmrui_txt":     "",            # directory
    "fsl_json":      "",            # directory
    "osprey_mat":    ".mat",
}


def export(basis: dict[str, np.ndarray],
           path: str,
           fmt: str,
           params: dict[str, Any] | None = None,
           *,
           extra_metadata: dict | None = None) -> str:
    """Export a basis dict to the requested format.

    Args:
        basis: { metabolite_name -> 1-D complex FID }.
        path: Output file (for single-file formats) or directory.
        fmt: One of SUPPORTED_FORMATS.
        params: Simulation parameters (must include 'Bandwidth', 'Center Freq',
                'Samples'; 'TE', 'Sequence', 'Nucleus' optional but recommended).
        extra_metadata: Free-form dict written verbatim into the sidecar JSON.

    Returns:
        Absolute path of the primary output (file or directory).
    """
    if fmt not in SUPPORTED_FORMATS:
        raise ValueError(f"Unknown export format '{fmt}'. Supported: {SUPPORTED_FORMATS}")
    params = dict(params or {})

    basis = _normalize_basis(basis)
    hdr = _make_header(params, basis)

    out = os.path.abspath(path)

    if fmt == "lcmodel_basis":
        _write_lcmodel_basis(basis, out, hdr, params)
    elif fmt == "lcmodel_raw":
        _ensure_dir(out)
        _write_lcmodel_raw_folder(basis, out, hdr, params)
    elif fmt == "jmrui_txt":
        _ensure_dir(out)
        _write_jmrui_folder(basis, out, hdr, params)
    elif fmt == "fsl_json":
        _ensure_dir(out)
        _write_fsl_json_folder(basis, out, hdr, params)
    elif fmt == "osprey_mat":
        _write_osprey_mat(basis, out, hdr, params)

    # Always emit a reproducibility sidecar
    sidecar_dir = out if os.path.isdir(out) else os.path.dirname(out)
    _ensure_dir(sidecar_dir)
    sidecar_path = os.path.join(sidecar_dir, "basis_sidecar.json")
    _write_sidecar(sidecar_path, fmt, out, basis, params, extra_metadata)

    return out


# ----------------------------- internals -------------------------------------

def _normalize_basis(basis: dict[str, Any]) -> dict[str, np.ndarray]:
    """Coerce every entry to a 1-D complex64/128 numpy array."""
    norm = {}
    for name, fid in basis.items():
        arr = np.asarray(fid)
        if arr.dtype.kind != "c":
            arr = arr.astype(np.complex128)
        if arr.ndim > 1:
            arr = arr.flatten()
        norm[str(name)] = arr
    return norm


def _b0_from_params(params: dict[str, Any]) -> float:
    """Resolve B0 [T] from any of the parameter spellings backends use.

    Priority: explicit `Bfield`/`B0` → `Field Strength` string ('3T') →
    derive from `Center Freq` (MHz) via 1H Larmor → 3.0 T fallback.
    Used by jMRUI / FSL-MRS / Osprey writers so MRSCloud-style exports
    (which only carry `Field Strength`) get a sensible B0 instead of 0.
    """
    for key in ("Bfield", "B0"):
        v = params.get(key)
        if v not in (None, "", "missing input", "Select option"):
            try:
                return float(v)
            except (TypeError, ValueError):
                pass
    fs = params.get("Field Strength")
    if fs not in (None, "", "missing input", "Select option"):
        try:
            return float(str(fs).replace("T", "").strip())
        except (TypeError, ValueError):
            pass
    cf = params.get("Center Freq") or params.get("centralFrequency")
    if cf not in (None, "", "missing input", "Select option"):
        try:
            cf_f = float(cf)
            cf_mhz = cf_f / 1e6 if cf_f > 1000 else cf_f
            return cf_mhz / 42.577
        except (TypeError, ValueError):
            pass
    return 3.0


def _make_header(params: dict[str, Any], basis: dict[str, np.ndarray]) -> dict[str, Any]:
    """Build a canonical header used by every writer."""
    if not basis:
        raise ValueError("Cannot build export header from an empty basis dict")

    bw_raw = params.get("Bandwidth") or params.get("SpectralWidth") or 2000.0
    try:
        bw = float(bw_raw)
    except (TypeError, ValueError):
        bw = 2000.0

    # Central / Larmor frequency in MHz.  Priority:
    #   1. Explicit 'Center Freq' (MHz or Hz — heuristic by magnitude)
    #   2. 'Bfield' / 'B0' in Tesla
    #   3. 'Field Strength' string ('1.5T' / '3T' / '7T') — used by MRSCloud
    cf = None
    cf_raw = params.get("Center Freq") or params.get("centralFrequency")
    if cf_raw not in (None, "", "missing input", "Select option"):
        try:
            cf_f = float(cf_raw)
            if cf_f != 0.0:
                cf = cf_f / 1e6 if cf_f > 1000 else cf_f
        except (TypeError, ValueError):
            cf = None
    if cf is None:
        cf = 42.577 * _b0_from_params(params)   # γ·B₀ [MHz]

    npts_raw = params.get("Samples")
    try:
        npts = int(npts_raw) if npts_raw not in (None, "") else next(iter(basis.values())).size
    except (TypeError, ValueError):
        npts = next(iter(basis.values())).size

    nucleus = str(params.get("Nucleus") or "1H")
    te = params.get("TE")
    try:
        te = float(te) if te not in (None, "") else None
    except (TypeError, ValueError):
        te = None

    try:
        fwhm = float(params.get("Linewidth", 1.0))
    except (TypeError, ValueError):
        fwhm = 1.0

    return {
        "bandwidth": bw,                # Hz
        "dwelltime": 1.0 / bw if bw else None,
        "centralFrequency": cf,         # MHz
        "points": npts,
        "nucleus": nucleus,
        "echotime": te,
        "fwhm": fwhm,
    }


def _ensure_dir(p: str) -> None:
    if p and not os.path.isdir(p):
        os.makedirs(p, exist_ok=True)


# ---- kbsct conversion toolbox bridge ---------------------------------------
# All on-disk format writing is delegated to the MRS Basis Set Conversion
# Toolbox (vendored as the ``kbsct`` git submodule under ``externals/``).
# BasisREMY only adapts its ``{name: complex FID}`` representation into the
# toolbox's per-metabolite "core" struct and calls the toolbox's tested writers,
# so the exported files match that community-validated implementation.

_KBSCT_WRITER_FILES = {
    "lcmodel": "write_lcmodel.py",
    "jmrui":   "write_jmrui.py",
    "fslmrs":  "write_fsLmrs.py",
    "osprey":  "write_osprey.py",
}

_KBSCT_MODULE_CACHE: dict[str, Any] = {}


def _kbsct_writers_dir() -> str:
    """Locate the toolbox ``writers/`` directory, fetching it if needed."""
    from basisremy.core.externals import ensure
    root = ensure("kbsct")
    preferred = os.path.join(root, "basis_converter", "writers")
    if os.path.isfile(os.path.join(preferred, "write_lcmodel.py")):
        return preferred
    # Layout-robust fallback: locate the writers package anywhere in the
    # checkout so an upstream reorganisation doesn't break the export.
    for dirpath, _dirs, files in os.walk(root):
        if os.path.basename(dirpath) == "writers" and "write_lcmodel.py" in files:
            return dirpath
    raise RuntimeError(f"Could not find the kbsct 'writers' directory under {root}.")


def _kbsct(module_key: str):
    """Import a kbsct writer module by file path (cached)."""
    mod = _KBSCT_MODULE_CACHE.get(module_key)
    if mod is not None:
        return mod
    import importlib.util
    fpath = os.path.join(_kbsct_writers_dir(), _KBSCT_WRITER_FILES[module_key])
    spec = importlib.util.spec_from_file_location(f"kbsct_{module_key}", fpath)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load kbsct writer module at {fpath}.")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    _KBSCT_MODULE_CACHE[module_key] = mod
    return mod


def _to_core_list(basis: dict[str, np.ndarray], hdr: dict[str, Any]) -> list[dict]:
    """Adapt a BasisREMY basis dict into kbsct 'core' structs.

    Each core carries the keys the toolbox writers expect:
    ``fid`` (complex FID), ``sw`` (spectral width, Hz), ``sf`` (Larmor
    frequency, MHz), ``n`` (points) and ``name``.
    """
    sw = float(hdr["bandwidth"])
    sf = float(hdr["centralFrequency"])   # Larmor frequency in MHz
    n = int(hdr["points"])
    return [
        {"fid": np.asarray(fid).ravel(), "sw": sw, "sf": sf, "n": n, "name": str(name)}
        for name, fid in basis.items()
    ]


# ---- LCModel .RAW per metabolite -------------------------------------------

def _write_lcmodel_raw_folder(basis, out_dir, hdr, params):
    _kbsct("lcmodel").write_lcmodel_raw_folder(_to_core_list(basis, hdr), out_dir)


# ---- LCModel .basis (single file) ------------------------------------------

def _write_lcmodel_basis(basis, out_file, hdr, params):
    if not out_file.lower().endswith(".basis"):
        out_file = out_file + ".basis"
    seq = str(params.get("Sequence") or "PRESS")
    _kbsct("lcmodel").write_lcmodel_basis(
        _to_core_list(basis, hdr), out_file,
        te=hdr["echotime"], seq=seq,
        description=f"BasisREMY {seq} basis",
    )


# ---- jMRUI text format -----------------------------------------------------

def _write_jmrui_folder(basis, out_dir, hdr, params):
    _kbsct("jmrui").write_jmrui_folder(_to_core_list(basis, hdr), out_dir)


# ---- FSL-MRS JSON basis directory ------------------------------------------

def _write_fsl_json_folder(basis, out_dir, hdr, params):
    _kbsct("fslmrs").write_fsLmrs_folder(_to_core_list(basis, hdr), out_dir)


# ---- Osprey .mat ------------------------------------------------------------
# Note: the toolbox's Osprey writer is not a 1:1 dump — it applies metabolite
# name-mapping, DC-correction, proton-equivalent rescaling and adds a synthetic
# H2O peak to produce an Osprey-compatible BASIS struct. ``add_mm=False`` keeps
# it from also injecting parametric MM/Lip basis functions.

def _write_osprey_mat(basis, out_file, hdr, params):
    if not out_file.lower().endswith(".mat"):
        out_file = out_file + ".mat"
    _ensure_dir(os.path.dirname(out_file))
    _kbsct("osprey").write_osprey(
        _to_core_list(basis, hdr), out_file,
        target_n=0, add_mm=False,
        te=float(hdr["echotime"] or 0.0),
        sequence=str(params.get("Sequence") or "unedited"),
    )


# ---- Reproducibility sidecar -----------------------------------------------

def _git_sha() -> str | None:
    try:
        sha = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
            stderr=subprocess.DEVNULL,
        ).decode().strip()
        return sha or None
    except Exception:
        return None


def _hash_file(path: str) -> str | None:
    try:
        h = _hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(65536), b""):
                h.update(chunk)
        return h.hexdigest()
    except Exception:
        return None


def _write_sidecar(path: str, fmt: str, output_path: str,
                   basis: dict[str, np.ndarray], params: dict, extra: dict | None) -> None:
    from basisremy.core import parameter_registry as _pr
    sidecar = {
        "tool": "BasisREMY",
        "tool_version": _basisremy_version(),
        "git_sha": _git_sha(),
        "timestamp_utc": _dt.datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "export_format": fmt,
        "output_path": output_path,
        "metabolites": list(basis.keys()),
        "n_points": int(next(iter(basis.values())).size) if basis else 0,
        "parameters": {k: _json_safe(v) for k, v in params.items()},
        "parameter_descriptions": {
            k: _pr.get(k).description for k in params.keys()
        },
        "pulse_file_sha256": (
            _hash_file(params["Path to Pulse"])
            if params.get("Path to Pulse") and os.path.isfile(params["Path to Pulse"]) else None
        ),
    }
    if extra:
        sidecar["extra"] = _json_safe(extra)
    with open(path, "w") as f:
        _json.dump(sidecar, f, indent=2, default=str)


def _basisremy_version() -> str:
    # No proper packaging metadata yet — use a placeholder.
    return "0.1.0-dev"


def _json_safe(v):
    if isinstance(v, (str, int, float, bool)) or v is None:
        return v
    if isinstance(v, dict):
        return {str(k): _json_safe(x) for k, x in v.items()}
    if isinstance(v, (list, tuple, set)):
        return [_json_safe(x) for x in v]
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, (np.integer, np.floating)):
        return v.item()
    return str(v)


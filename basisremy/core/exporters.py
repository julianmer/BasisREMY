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
import sys
from typing import Any, Iterable

import numpy as np

from basisremy.core.paths import externals_root

# Make sure the vendored fsl_mrs is importable
_FSL_MRS_ROOT = str(externals_root() / "fsl_mrs")
if _FSL_MRS_ROOT not in sys.path:
    sys.path.insert(0, _FSL_MRS_ROOT)


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


# ---- LCModel .RAW per metabolite -------------------------------------------

def _save_one_raw(fid: np.ndarray, filepath: str, hdr: dict, label: str) -> None:
    with open(filepath, "w") as f:
        f.write(" $NMID\n")
        f.write(f" ID='{label}'\n")
        f.write(" FMTDAT='(2E15.6)'\n")
        f.write(" VOLUME=1.0\n")
        f.write(" TRAMP=1.0\n")
        f.write(" $END\n")
        for pt in fid:
            f.write(f" {pt.real:15.6E} {pt.imag:15.6E}\n")


def _write_lcmodel_raw_folder(basis, out_dir, hdr, params):
    for name, fid in basis.items():
        _save_one_raw(fid, os.path.join(out_dir, f"{name}.RAW"),
                      hdr, label=f"BasisREMY {name}")


# ---- LCModel .basis (single file) ------------------------------------------
# Implements the LCModel BASIS file ASCII format. We construct it directly so
# we don't have to depend on FSL-MRS scripts.

def _write_lcmodel_basis(basis, out_file, hdr, params):
    if not out_file.lower().endswith(".basis"):
        out_file = out_file + ".basis"
    _ensure_dir(os.path.dirname(out_file))

    seq = str(params.get("Sequence") or "PRESS")
    te = hdr["echotime"] if hdr["echotime"] is not None else 30.0
    bw = hdr["bandwidth"]
    cf = hdr["centralFrequency"]
    npts = hdr["points"]

    with open(out_file, "w") as f:
        # SEQPAR section
        f.write(" $SEQPAR\n")
        f.write(f"  ECHOT  = {te:.4f}\n")
        f.write(f"  HZPPPM = {cf:.4f}\n")
        f.write(f"  SEQ    = '{seq}'\n")
        f.write(" $END\n")
        # BASIS1 header
        f.write(" $BASIS1\n")
        f.write(f"  FMTBAS = '(2E15.6)'\n")
        f.write(f"  IDBASI = 'BasisREMY'\n")
        f.write(f"  BADELT = {1.0/bw:.6E}\n")
        f.write(f"  NDATAB = {npts}\n")
        f.write(" $END\n")
        # one BASIS section per metabolite
        for name, fid in basis.items():
            fid = fid[:npts] if fid.size >= npts else np.pad(fid, (0, npts - fid.size))
            f.write(" $BASIS\n")
            f.write(f"  ID     = '{name}'\n")
            f.write(f"  METABO = '{name}'\n")
            f.write(f"  CONC   = 1.0\n")
            f.write(f"  TRAMP  = 1.0\n")
            f.write(f"  VOLUME = 1.0\n")
            f.write(f"  ISHIFT = 0\n")
            f.write(" $END\n")
            for pt in fid:
                f.write(f" {pt.real:15.6E} {pt.imag:15.6E}\n")


# ---- jMRUI text format -----------------------------------------------------

def _write_jmrui_folder(basis, out_dir, hdr, params):
    bw = hdr["bandwidth"]
    cf_hz = hdr["centralFrequency"] * 1e6   # MHz -> Hz
    dwell_ms = 1000.0 / bw                  # ms
    for name, fid in basis.items():
        fp = os.path.join(out_dir, f"{name}.txt")
        with open(fp, "w") as f:
            f.write("jMRUI Data Textfile\n\n")
            f.write(f"Filename: {name}.txt\n\n")
            f.write(f"PointsInDataset: {fid.size}\n")
            f.write(f"DatasetsInFile: 1\n")
            f.write(f"SamplingInterval: {dwell_ms:.6E}\n")
            f.write(f"ZeroOrderPhase: 0E0\n")
            f.write(f"BeginTime: 0E0\n")
            f.write(f"TransmitterFrequency: {cf_hz:.6E}\n")
            f.write(f"MagneticField: {_b0_from_params(params):.4f}\n")
            f.write(f"TypeOfNucleus: {hdr['nucleus']}\n")
            f.write(f"NameOfPatient: BasisREMY\n")
            f.write(f"DateOfExperiment: {_dt.date.today().isoformat()}\n")
            f.write(f"Spectrometer: BasisREMY-simulated\n")
            f.write(f"AdditionalInfo: metabolite={name}\n\n")
            f.write("Signal and FFT\n")
            f.write("sig(real)\tsig(imag)\tfft(real)\tfft(imag)\n")
            # fftshift(fft()) → spectrum from low freq (low ppm) to high freq
            # (high ppm). jMRUI reads this column positionally and inverts
            # the ppm axis on display, so the order is consistent.
            spec = np.fft.fftshift(np.fft.fft(fid))
            for fp_t, fp_w in zip(fid, spec):
                f.write(f"{fp_t.real:.6E}\t{fp_t.imag:.6E}\t{fp_w.real:.6E}\t{fp_w.imag:.6E}\n")


# ---- FSL-MRS JSON basis directory ------------------------------------------
# Uses the FSL-MRS write_fsl_basis_file helper if available, otherwise writes
# a minimal JSON directly compatible with fsl_mrs.utils.mrs_io.fsl_io.readFSLBasis.

def _write_fsl_json_folder(basis, out_dir, hdr, params):
    b0 = _b0_from_params(params)
    info = f"BasisREMY simulated basis ({params.get('Sequence', '?')}, " \
           f"TE={hdr['echotime']}ms, B0={b0:.2f}T)"

    # Try FSL-MRS helper for full compatibility
    try:
        from basisremy.core.externals import ensure as _ensure_external
        _ensure_external('fsl_mrs')
        from fsl_mrs.utils.mrs_io.fsl_io import write_fsl_basis_file
        # Header that write_fsl_basis_file expects
        fsl_hdr = {
            "centralFrequency": hdr["centralFrequency"],
            "bandwidth":        hdr["bandwidth"],
            "dwelltime":        hdr["dwelltime"],
            "fwhm":             hdr["fwhm"],
        }
        for name, fid in basis.items():
            write_fsl_basis_file(fid, name, fsl_hdr, out_dir, info=info)
        # Ensure the top-level keys the test (and FSL-MRS readers) expect
        # are always present regardless of the helper's output schema.
        for name in basis:
            fp = os.path.join(out_dir, f"{name}.json")
            if os.path.isfile(fp):
                try:
                    with open(fp) as _f:
                        _payload = _json.load(_f)
                    changed = False
                    if "centralFrequency" not in _payload:
                        _payload["centralFrequency"] = hdr["centralFrequency"]
                        changed = True
                    if "spectralwidth" not in _payload:
                        _payload["spectralwidth"] = hdr["bandwidth"]
                        changed = True
                    if changed:
                        with open(fp, "w") as _f:
                            _json.dump(_payload, _f, indent=2)
                except Exception:
                    pass
        return
    except Exception:
        pass

    # Fallback: write minimal JSON files directly (matches readFSLBasis layout)
    for name, fid in basis.items():
        fp = os.path.join(out_dir, f"{name}.json")
        payload = {
            "basis": {
                "basis_re": fid.real.tolist(),
                "basis_im": fid.imag.tolist(),
                "basis_name": name,
                "basis_centre": 4.65,
                "basis_width": hdr["fwhm"],
            },
            "MM": False,
            "info": info,
            "seq": {
                "TE": hdr["echotime"],
                "B0": b0,
                "SequenceName": params.get("Sequence"),
                "Nucleus": hdr["nucleus"],
            },
            "spectralwidth": hdr["bandwidth"],
            "centralFrequency": hdr["centralFrequency"],
            "dwelltime": hdr["dwelltime"],
        }
        with open(fp, "w") as f:
            _json.dump(payload, f, indent=2)


# ---- Osprey .mat ------------------------------------------------------------

def _write_osprey_mat(basis, out_file, hdr, params):
    if not out_file.lower().endswith(".mat"):
        out_file = out_file + ".mat"
    _ensure_dir(os.path.dirname(out_file))
    try:
        from scipy.io import savemat
    except ImportError as e:
        raise RuntimeError("scipy is required for Osprey .mat export") from e

    names = list(basis.keys())
    fids = np.column_stack([basis[n] for n in names])  # (Npts, Nmetabs)
    bw = hdr["bandwidth"]
    cf = hdr["centralFrequency"]
    npts = hdr["points"]

    spec = np.fft.fftshift(np.fft.fft(fids, axis=0), axes=0)
    t = np.arange(npts) / bw

    # ppm axis: f_offset [Hz] / cf [MHz] = ppm directly (unit convenience).
    # Result goes from ~(4.65 - bw/(2cf)) to ~(4.65 + bw/(2cf)) ppm
    # (low ppm → high ppm). Osprey inverts the x-axis on its own display.
    ppm_axis = np.linspace(-bw / 2, bw / 2, npts) / cf + 4.65

    BASIS = {
        "fids": fids.astype(np.complex128),
        "specs": spec.astype(np.complex128),
        "name": np.array(names, dtype=object),
        "t": t,
        "ppm": ppm_axis,
        "spectralwidth": float(bw),
        "dwelltime": float(1.0 / bw),
        "n": int(npts),
        "linewidth": float(hdr["fwhm"]),
        "Bo": _b0_from_params(params),
        "centerFreq": float(cf),
        "te": float(hdr["echotime"] or 0.0),
        "seq": str(params.get("Sequence") or ""),
        "nMets": len(names),
        "sz": np.array(fids.shape, dtype=np.int32),
    }
    savemat(out_file, {"BASIS": BASIS}, do_compression=True, oned_as="column")


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


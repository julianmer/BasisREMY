####################################################################################################
#                                          rf_pulses.py                                            #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 28/08/26                                                                                #
#                                                                                                  #
# Purpose: Load RF pulse waveform files and scale them to a target flip angle for the Python-side  #
#          simulators (Vespa/PyGAMMA). This is a port of FID-A's io_loadRFwaveform (plus the       #
#          headless w1max search of adapters/backends/io_loadRFwaveform.m), so a pulse file gives  #
#          the same B1 scaling here as in the Octave backends.                                     #
#                                                                                                  #
#          Formats: Siemens .pta, Varian/Agilent .RF, FID-A basic .txt (amp phase [timestep]).     #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import math
import os

import numpy as np


def read_waveform(path: str) -> np.ndarray:
    """Return an (N, 3) array of [phase_deg, amplitude, timestep] like FID-A's rf matrix."""
    ext = os.path.splitext(path)[1].lower()
    if ext == '.pta':
        return _read_pta(path)
    if ext == '.rf':
        return _read_rf(path)
    if ext == '.txt':
        return _read_txt(path)
    raise ValueError(f"Unrecognised RF pulse file '{os.path.basename(path)}' "
                     f"(supported: .pta, .RF, .txt)")


def _read_pta(path):
    """Siemens .pta: header 'KEY: value' lines, then 'amp phase ; (i)' lines (phase in rad)."""
    rows = []
    with open(path, 'r', errors='replace') as f:
        for line in f:
            if ';' not in line:
                continue
            data = line.split(';')[0].split()
            if len(data) < 2:
                continue
            try:
                amp, phase = float(data[0]), float(data[1])
            except ValueError:
                continue
            rows.append((math.degrees(phase), amp, 1.0))
    if not rows:
        raise ValueError(f"No waveform data found in {path}")
    return np.array(rows, dtype=float)


def _read_rf(path):
    """Varian/Agilent .RF: comment lines start with '#'; data 'phase amp gate'."""
    rows = []
    with open(path, 'r', errors='replace') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            parts = s.split()
            try:
                vals = [float(v) for v in parts[:3]]
            except ValueError:
                continue
            if len(vals) == 2:
                vals.append(1.0)
            rows.append(vals)
    if not rows:
        raise ValueError(f"No waveform data found in {path}")
    return np.array(rows, dtype=float)


def _read_txt(path):
    """FID-A basic .txt: columns amplitude, phase [deg], optional timestep (and gradient)."""
    data = np.loadtxt(path, comments=('#', '%'))
    if data.ndim == 1:
        data = data[None, :]
    if data.shape[1] < 2:
        raise ValueError(f"{path}: need at least amplitude and phase columns")
    amp, phase = data[:, 0], data[:, 1]
    step = data[:, 2] if data.shape[1] >= 3 else np.ones(len(amp))
    return np.column_stack([phase, amp, step])


def _bloch_mz_after_pulse(phase_deg, amp_norm, timestep, tp_s, w1_hz):
    """Mz after the pulse for each w1 (Hz), on resonance, M0 = +z, no relaxation.

    Vectorised over w1: each step is a rotation about the in-plane axis at the
    step's phase by 2*pi*w1*amp*dt.
    """
    dt = tp_s * timestep / timestep.sum()
    w1 = np.asarray(w1_hz, dtype=float)
    mx = np.zeros_like(w1); my = np.zeros_like(w1); mz = np.ones_like(w1)
    for phi_d, a, d in zip(phase_deg, amp_norm, dt):
        theta = 2.0 * math.pi * w1 * a * d
        phi = math.radians(phi_d)
        ux, uy = math.cos(phi), math.sin(phi)
        c, s = np.cos(theta), np.sin(theta)
        # Rodrigues rotation about unit axis (ux, uy, 0)
        dot = ux * mx + uy * my
        nmx = mx * c + (uy * mz) * s + ux * dot * (1 - c)
        nmy = my * c + (-ux * mz) * s + uy * dot * (1 - c)
        nmz = mz * c + (ux * my - uy * mx) * s
        mx, my, mz = nmx, nmy, nmz
    return mz


def scale_waveform(rf: np.ndarray, tp_s: float, flip='ref'):
    """FID-A io_loadRFwaveform scaling: returns (phase_deg, amp_hz, dt_s).

    flip: 'exc' (90), 'ref'/'inv' (180) or a numeric flip angle in degrees.
    Amplitude-modulated pulses are scaled by their integral; phase-modulated
    (adiabatic / GOIA) pulses by a Bloch sweep of w1max from 0 to 5 kHz that
    takes the lowest w1max reaching the target Mz (FID-A's headless search).
    """
    rf = np.array(rf, dtype=float)
    phase = rf[:, 0].copy()
    amp = rf[:, 1].copy()
    step = rf[:, 2].copy() if rf.shape[1] >= 3 else np.ones(len(amp))

    # remove 360-degree wraps (io_loadRFwaveform)
    jumps = np.diff(phase)
    for idx in np.where((np.abs(jumps) > 355) & (np.abs(jumps) < 365))[0]:
        phase[idx + 1:] -= 360.0 * np.sign(jumps[idx])
    is_phase_modulated = bool(np.any((np.round(phase) != 180) & (np.round(phase) != 0)))

    amp = amp / np.max(np.abs(amp))

    if isinstance(flip, str):
        flip_cyc = {'exc': 0.25, 'ref': 0.5, 'inv': 0.5}[flip]
        target = {'exc': 0.0, 'ref': -1.0, 'inv': -1.0}[flip]
    else:
        flip_cyc = float(flip) / 360.0
        target = math.cos(math.radians(float(flip)))

    if not is_phase_modulated:
        sign = np.where(phase > 179, -1.0, 1.0)
        int_rf = float(np.sum(amp * sign) / len(amp))
        w1max = flip_cyc / (int_rf * tp_s) if int_rf else 0.0
    else:
        sweep = np.linspace(0.0, 5000.0, 40000)
        mz = _bloch_mz_after_pulse(phase, amp, step, tp_s, sweep)
        hit = np.where(mz <= target + 0.02)[0]
        idx = hit[0] if len(hit) else int(np.argmin(np.abs(mz - target)))
        w1max = float(sweep[idx])

    dt = tp_s * step / step.sum()
    return phase, amp * w1max, dt


def load_pulse(path: str, tp_ms: float, flip='ref') -> dict:
    """File -> worker-ready pulse dict (phase in deg, amplitude in Hz, dwell in s)."""
    rf = read_waveform(path)
    phase, amp_hz, dt = scale_waveform(rf, float(tp_ms) / 1000.0, flip)
    return {'phase_deg': phase.tolist(), 'amp_hz': amp_hz.tolist(),
            'dt_s': dt.tolist(), 'name': os.path.basename(path)}

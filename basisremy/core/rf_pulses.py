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
#          Formats: Siemens .pta, Varian/Agilent .RF, FID-A basic .txt (amp phase [timestep]),     #
#          Bruker JCAMP-DX shapes (.exc / .rfc / .inv) and a MATLAB .mat holding an FID-A RF       #
#          struct (its 'waveform' field, as MRSCloud's GOIA files).                                #
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
    if ext in BRUKER_EXTS:
        return _read_bruker(path)
    if ext == '.mat':
        return _read_mat(path)
    raise ValueError(f"Unrecognised RF pulse file '{os.path.basename(path)}' "
                     f"(supported: .pta, .RF, .txt, Bruker .exc/.rfc/.inv, FID-A .mat)")


BRUKER_EXTS = ('.exc', '.rfc', '.inv')


def _read_bruker(path):
    """Bruker JCAMP-DX shape: '##XYPOINTS=(XY..XY)' then 'amplitude, phase' lines (amplitude in
    percent, phase in degrees) up to '##END' - the layout FID-A's io_readRFBruk reads."""
    rows, data = [], False
    with open(path, 'r', errors='replace') as f:
        for line in f:
            s = line.strip()
            if s.startswith('##XYPOINTS'):
                data = True
                continue
            if not data or not s:
                continue
            if s.startswith('##'):
                break
            try:
                amp, phase = (float(v) for v in s.split(',')[:2])
            except ValueError:
                continue
            rows.append((phase, amp, 1.0))
    if not rows:
        raise ValueError(f"No ##XYPOINTS waveform found in {path}")
    return np.array(rows, dtype=float)


def _read_mat(path):
    """MATLAB .mat with an FID-A RF struct: its 'waveform' [phase deg, amp, step(, gradient)]."""
    from scipy.io import loadmat
    d = loadmat(path, squeeze_me=True, struct_as_record=False)
    for name, v in d.items():
        if name.startswith('__'):
            continue
        wave = getattr(v, 'waveform', None)
        if wave is not None:
            wave = np.atleast_2d(np.asarray(wave, dtype=float))
            if wave.shape[1] >= 3:
                return wave
    raise ValueError(f"{os.path.basename(path)}: no FID-A RF struct (a variable with a "
                     f"'waveform' field) in the file")


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
    cols = [phase, amp, step]
    if data.shape[1] >= 4:                       # gradient column (GOIA), kept for the engines
        cols.append(data[:, 3])
    return np.column_stack(cols)


def _bloch_mz_after_pulse(phase_deg, amp_norm, timestep, tp_s, w1_hz, bz_hz=None):
    """Mz after the pulse for each w1 (Hz), M0 = +z, no relaxation.

    Vectorised over w1. bz_hz is the longitudinal field per step (N,) or per
    step and w1 column (N, M) in Hz: frequency offset plus gradient x position;
    None means on resonance.
    """
    dt = tp_s * timestep / timestep.sum()
    w1 = np.asarray(w1_hz, dtype=float)
    mx = np.zeros_like(w1); my = np.zeros_like(w1); mz = np.ones_like(w1)
    bz_all = None if bz_hz is None else np.asarray(bz_hz, dtype=float)
    for i, (phi_d, a, d) in enumerate(zip(phase_deg, amp_norm, dt)):
        phi = math.radians(phi_d)
        bx, by = w1 * a * math.cos(phi), w1 * a * math.sin(phi)
        bz = np.zeros_like(w1) if bz_all is None else np.broadcast_to(bz_all[i], w1.shape)
        bmag = np.sqrt(bx ** 2 + by ** 2 + bz ** 2)
        theta = 2.0 * math.pi * bmag * d
        safe = np.where(bmag > 0, bmag, 1.0)
        ux, uy, uz = bx / safe, by / safe, bz / safe
        c, s = np.cos(theta), np.sin(theta)
        # Rodrigues rotation about the unit axis (ux, uy, uz)
        dot = ux * mx + uy * my + uz * mz
        nmx = mx * c + (uy * mz - uz * my) * s + ux * dot * (1 - c)
        nmy = my * c + (uz * mx - ux * mz) * s + uy * dot * (1 - c)
        nmz = mz * c + (ux * my - uy * mx) * s + uz * dot * (1 - c)
        mx, my, mz = nmx, nmy, nmz
    return mz


GYRO_HZ_PER_G = 4257.7          # 1H, as in FID-A's bes / rf_goia


def _adiabatic_w1max(phase, amp, step, tp_s, target, grad=None):
    """w1max (Hz) of a phase-modulated pulse, as the headless io_loadRFwaveform.m.

    First the lowest w1max whose on-resonance Mz reaches the target (FID-A's
    plot-and-type step made deterministic). For an inversion (target -1) that
    is only the adiabatic threshold of the centre of the slab, so the search
    then looks at the whole slab: the Mz < 0 extent at that w1max (position
    with the pulse's gradient column, else frequency) gives the slab, and the
    lowest w1max whose mean Mz over nine points across its central 80 % is
    within 0.01 of the best achievable is taken, never below the threshold.
    """
    sweep = np.linspace(0.0, 5000.0, 40000)
    mz = _bloch_mz_after_pulse(phase, amp, step, tp_s, sweep)
    hit = np.where(mz <= target + 0.02)[0]
    idx = hit[0] if len(hit) else int(np.argmin(np.abs(mz - target)))
    w1_0 = float(sweep[idx])
    if target != -1.0 or w1_0 == 0.0:
        return w1_0
    gm = grad is not None and np.any(grad)
    axis = np.linspace(-5.0, 5.0, 2001) if gm else np.linspace(-5000.0, 5000.0, 2001)   # cm | Hz, bes window
    bz = (GYRO_HZ_PER_G * np.asarray(grad, dtype=float)[:, None] * axis[None, :] if gm
          else np.repeat(axis[None, :], len(amp), axis=0))
    prof = _bloch_mz_after_pulse(phase, amp, step, tp_s, np.full(axis.size, w1_0), bz)
    neg = axis[prof < 0]
    if neg.size < 2:
        return w1_0
    half = (neg.max() - neg.min()) / 2.0
    probes = np.linspace(-0.8, 0.8, 9) * half
    sweep = np.arange(25.0, 5000.0 + 25.0, 25.0)
    w1 = np.repeat(sweep, probes.size)
    pos = np.tile(probes, sweep.size)
    bz = (GYRO_HZ_PER_G * np.asarray(grad, dtype=float)[:, None] * pos[None, :] if gm
          else np.repeat(pos[None, :], len(amp), axis=0))
    mean = _bloch_mz_after_pulse(phase, amp, step, tp_s, w1, bz).reshape(sweep.size, probes.size).mean(axis=1)
    best = float(sweep[np.where(mean <= mean.min() + 0.01)[0][0]])
    return max(w1_0, best)


def scale_waveform(rf: np.ndarray, tp_s: float, flip='ref'):
    """FID-A io_loadRFwaveform scaling: returns (phase_deg, amp_hz, dt_s).

    flip: 'exc' (90), 'ref'/'inv' (180) or a numeric flip angle in degrees.
    Amplitude-modulated pulses are scaled by their integral; phase-modulated
    (adiabatic / GOIA) pulses by a Bloch sweep of w1max (_adiabatic_w1max): the
    on-resonance threshold, and for an inversion the slab-wide search over the
    pulse's gradient column (rf[:, 3]) or its frequency band.
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
        grad = rf[:, 3] if rf.shape[1] >= 4 else None
        w1max = _adiabatic_w1max(phase, amp, step, tp_s, target, grad)

    dt = tp_s * step / step.sum()
    return phase, amp * w1max, dt


def load_pulse(path: str, tp_ms: float, flip='ref') -> dict:
    """File -> worker-ready pulse dict (phase in deg, amplitude in Hz, dwell in s)."""
    rf = read_waveform(path)
    phase, amp_hz, dt = scale_waveform(rf, float(tp_ms) / 1000.0, flip)
    return {'phase_deg': phase.tolist(), 'amp_hz': amp_hz.tolist(),
            'dt_s': dt.tolist(), 'name': os.path.basename(path)}

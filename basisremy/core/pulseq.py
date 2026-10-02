####################################################################################################
#                                            pulseq.py                                             #
####################################################################################################
#                                                                                                  #
# Purpose: Read a Pulseq (.seq) sequence file into the parameter sheet. A .seq file carries the    #
#          exact RF waveforms, flip angles, frequency offsets, gradients and block timing of the   #
#          acquisition, so the echo times, the refocusing pulse, its duration and the slab         #
#          thicknesses come from the sequence itself instead of being typed or pre-filled.         #
#                                                                                                  #
#          The file holds no field strength and no sequence name: B0 comes from the data header,   #
#          the sequence type from the RF events (1 excitation + 2 refocusing = PRESS, + 4 = semi-  #
#          LASER, 3 excitations = STEAM, PRESS + 2 non-selective inversions = MEGA-PRESS). Echo    #
#          times follow from the RF centres (the ADC may start before the echo). Editing           #
#          frequencies take the carrier on water (4.65 ppm) and need B0 from the data header.      #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass

import numpy as np


@dataclass
class RFEvent:
    """One RF pulse of the sequence (times in ms from the start of the file)."""
    centre_ms: float
    dur_ms: float
    flip_deg: float
    phase_rad: float
    freq_offset_hz: float
    signal: np.ndarray            # complex B1 samples [Hz]
    slab_axis: str | None         # gradient axis under the pulse (x / y / z)
    slab_grad_hz_m: float | None  # slice-select gradient at the pulse centre [Hz/m]
    role: str                     # 'exc', 'ref' (slice-selective refocusing / inversion), 'edit'
                                  # (non-selective inversion, MEGA) or 'sat' (water suppression)
    adiabatic: bool               # frequency-swept (phase-modulated) waveform
    freq_ppm: float = 0.0         # offset in ppm of the carrier (Pulseq 1.5), besides freq_offset_hz
    grad_hz_m: np.ndarray | None = None  # slice gradient per RF sample when it varies (GOIA) [Hz/m]


@dataclass
class SeqInfo:
    """What the sheet needs from one acquisition of a .seq file (the first ADC)."""
    rf: list[RFEvent]
    adc_start_ms: float
    samples: int
    dwell_s: float
    name: str | None


def read_seq(path: str, acquisition: int = 0) -> SeqInfo:
    """RF events of one acquisition (those after the previous ADC) and its ADC (PyPulseq)."""
    import pypulseq as pp
    seq = pp.Sequence()
    seq.read(path)
    rf, t, n_adc = [], 0.0, 0
    for bid in sorted(seq.block_durations):
        b = seq.get_block(bid)
        if b.rf is not None:
            rf.append(_rf_event(b, t, pp))
        if b.adc is not None:
            if n_adc == acquisition:
                return SeqInfo(rf=rf, adc_start_ms=(t + b.adc.delay) * 1e3,
                               samples=int(b.adc.num_samples), dwell_s=float(b.adc.dwell),
                               name=seq.definitions.get('Name'))
            n_adc, rf = n_adc + 1, []
        t += seq.block_durations[bid]
    if n_adc == 0:
        raise ValueError(f"{os.path.basename(path)}: no ADC event in the sequence file.")
    raise ValueError(f"{os.path.basename(path)}: {n_adc} acquisition(s), not {acquisition + 1}.")


def _grad_at(g, times):
    """Gradient amplitude [Hz/m] of a block's gradient event at times within the block [s]."""
    if g.type == 'trap':
        tt = g.delay + np.array([0.0, g.rise_time, g.rise_time + g.flat_time,
                                 g.rise_time + g.flat_time + g.fall_time])
        return np.interp(times, tt, [0.0, g.amplitude, g.amplitude, 0.0], left=0.0, right=0.0)
    return np.interp(times, g.delay + np.asarray(g.tt), np.asarray(g.waveform), left=0.0, right=0.0)


def _rf_event(b, t, pp) -> RFEvent:
    """The RF pulse of block b, which starts t seconds into the file."""
    p = b.rf
    # sample spacing; the duration is samples x spacing (shape_dur ends half a sample short)
    dt = float(np.mean(np.diff(p.t))) if len(p.t) > 1 else float(p.shape_dur)
    tc, _ = pp.calc_rf_center(p)
    t_rel = p.delay + tc                               # centre within the block [s]
    axis, g_c, g_ev = None, None, None
    for ax in ('gx', 'gy', 'gz'):
        g = getattr(b, ax)
        if g is None:
            continue
        val = float(_grad_at(g, t_rel))
        if abs(val) > abs(g_c or 0.0):
            axis, g_c, g_ev = ax[1], val, g
    signal = np.asarray(p.signal, dtype=complex)
    flip = float(np.degrees(2 * np.pi * abs(np.sum(signal) * dt)))
    swept = _phase_modulated(signal)
    grad = None
    if g_ev is not None:                               # gradient-modulated (GOIA) pulse?
        gs = _grad_at(g_ev, p.delay + np.asarray(p.t))
        on = np.abs(signal) > 0.01 * np.abs(signal).max()
        if np.ptp(gs[on]) > 0.05 * np.abs(gs[on]).max():
            grad = gs
    use = str(getattr(p, 'use', 'undefined'))
    selective = axis is not None
    if use == 'excitation':
        role = 'exc'
    elif use in ('saturation', 'preparation'):
        role = 'sat'
    elif use in ('refocusing', 'inversion') or swept or flip >= 120:
        # older files: an adiabatic sweep refocuses; its signed area says nothing
        role = 'ref' if selective else 'edit'
    else:
        role = 'exc' if selective else 'sat'
    return RFEvent(
        centre_ms=(t + t_rel) * 1e3, dur_ms=len(signal) * dt * 1e3, flip_deg=flip,
        phase_rad=float(p.phase_offset), freq_offset_hz=float(p.freq_offset),
        signal=signal, slab_axis=axis, slab_grad_hz_m=abs(g_c) if g_c else None,
        role=role, adiabatic=swept, freq_ppm=float(getattr(p, 'freq_ppm', 0.0) or 0.0),
        grad_hz_m=grad)


def _phase_modulated(signal: np.ndarray) -> bool:
    """True when the phase is not just 0 / 180 degrees (a frequency sweep)."""
    on = np.abs(signal) > 0.01 * np.abs(signal).max()
    return bool(np.mean(np.abs(np.sin(np.angle(signal[on]))) > 0.1) > 0.2)


# sequence type -> FID-A shaped backend that simulates it
BACKENDS = {'PRESS': 'FidaPressShaped', 'sLASER': 'FidaSemiLaserShaped',
            'STEAM': 'FidaSteamShaped', 'MEGA-PRESS': 'FidaMegaPressShaped'}
WATER_PPM = 4.65                  # the carrier: scanners centre on water (FID-A's water shift)


def sequence_type(info: SeqInfo) -> str:
    """PRESS / semi-LASER / STEAM / MEGA-PRESS from the RF events before the readout."""
    n = tuple(sum(r.role == k for r in info.rf) for k in ('exc', 'ref', 'edit'))
    kinds = {(1, 2, 0): 'PRESS', (1, 4, 0): 'sLASER', (3, 0, 0): 'STEAM', (1, 2, 2): 'MEGA-PRESS'}
    if n in kinds:
        return kinds[n]
    raise ValueError(f"Sequence with {n[0]} excitation, {n[1]} refocusing and {n[2]} editing pulses "
                     "is not supported from a .seq file yet (PRESS, semi-LASER, STEAM and "
                     "MEGA-PRESS are).")


def echo_ms(info: SeqInfo) -> float:
    """Echo time from the RF centres: each refocusing pulse mirrors the previous echo; STEAM's
    stimulated echo follows the third 90 by the first-to-second interval (TE = 2 x that)."""
    exc = [r for r in info.rf if r.role == 'exc']
    if len(exc) == 3:
        return 2 * (exc[1].centre_ms - exc[0].centre_ms)
    echo = exc[0].centre_ms
    for r in (r for r in info.rf if r.role == 'ref'):
        echo = 2 * r.centre_ms - echo
    return echo - exc[0].centre_ms


def edit_ppm(event: RFEvent, bfield: float) -> float:
    """Frequency of an editing pulse in ppm (carrier on water)."""
    return WATER_PPM + event.freq_ppm + event.freq_offset_hz / (float(bfield) * 42.577)


def waveform_pulse(event: RFEvent, kind: str = 'ref'):
    """The RF event as a pulse_library.Pulse (FID-A waveform convention; a gradient-modulated
    pulse carries its gradient in G/cm as the fourth column)."""
    from basisremy.core.pulse_library import GYRO_HZ_PER_G, Pulse
    amp = np.abs(event.signal)
    cols = [np.degrees(np.angle(event.signal)) % 360.0, amp / amp.max(), np.ones(len(amp))]
    if event.grad_hz_m is not None:
        cols.append(np.abs(event.grad_hz_m) / GYRO_HZ_PER_G / 100.0)
    return Pulse('seq', np.column_stack(cols), kind, event.dur_ms)


def slab_cm(event: RFEvent) -> float | None:
    """Slab thickness: the pulse bandwidth (FID-A's measure) over the slice gradient; for a
    gradient-modulated pulse the Mz < 0 extent of a Bloch simulation over position with the
    file's own B1 and gradient waveforms."""
    if not event.slab_grad_hz_m:
        return None
    if event.grad_hz_m is not None:
        return _slab_gm_cm(event)
    from basisremy.core.pulse_library import bandwidth_hz
    pulse = waveform_pulse(event, 'exc' if event.role == 'exc' else 'ref')
    bw = bandwidth_hz(pulse, event.dur_ms, points=8001)
    return bw / event.slab_grad_hz_m * 100.0


def _slab_gm_cm(event: RFEvent, half_m: float = 0.1, points: int = 8001) -> float:
    x = np.linspace(-half_m, half_m, points)
    m = np.zeros((3, points)); m[2] = 1.0
    dt = event.dur_ms / 1e3 / len(event.signal)
    for b1, g in zip(event.signal, event.grad_hz_m):  # Hz, Hz/m
        bx, by, bz = b1.real, b1.imag, g * x
        bmag = np.sqrt(bx ** 2 + by ** 2 + bz ** 2)
        safe = np.where(bmag > 0, bmag, 1.0)
        ux, uy, uz = bx / safe, by / safe, bz / safe
        c, s = np.cos(2 * np.pi * bmag * dt), np.sin(2 * np.pi * bmag * dt)
        dot = ux * m[0] + uy * m[1] + uz * m[2]
        m = np.array([m[0] * c + (uy * m[2] - uz * m[1]) * s + ux * dot * (1 - c),
                      m[1] * c + (uz * m[0] - ux * m[2]) * s + uy * dot * (1 - c),
                      m[2] * c + (ux * m[1] - uy * m[0]) * s + uz * dot * (1 - c)])
    idx = np.where(m[2] < 0)[0]
    return float(x[idx[-1]] - x[idx[0]]) * 100.0 if idx.size else 0.0


def sheet_params(path: str, workdir: str, bfield: float | None = None) -> tuple[str, dict]:
    """(backend, registry-named parameters) for a FID-A shaped simulation of the .seq.

    Writes the shaped waveform (refocusing; STEAM: the second 90; MEGA: also the editing
    pulse) into workdir as FID-A .txt files. TE, the PRESS echo split, STEAM's TM, duration,
    flip angle, slabs and the readout come from the file; MEGA's editing frequencies also
    need B0 (bfield, from the data header)."""
    info = read_seq(path)
    kind = sequence_type(info)
    exc = [r for r in info.rf if r.role == 'exc']
    ref = [r for r in info.rf if r.role == 'ref']
    te = echo_ms(info)
    os.makedirs(workdir, exist_ok=True)
    # the pulses FID-A simulates as shaped (its excitation is instantaneous)
    shaped, name = (exc[1:], 'seq_excitation.txt') if kind == 'STEAM' else (ref, 'seq_refocusing.txt')
    pulse = waveform_pulse(shaped[0], 'exc' if kind == 'STEAM' else 'ref')
    params = {
        'TE': round(te, 4),
        'Samples': info.samples,
        'Bandwidth': round(1.0 / info.dwell_s, 4),
        'Path to Pulse': pulse.write(os.path.join(workdir, name)),
        'RefTp': round(shaped[0].dur_ms, 4),
        'Flip Angle': 180.0 if shaped[0].adiabatic else round(shaped[0].flip_deg, 2),
    }
    for axis, ev in zip(('X', 'Y'), (shaped[0], shaped[-1])):
        thk = slab_cm(ev)
        if thk:
            params[f'thk{axis}'] = round(thk, 3)
            params[f'fov{axis}'] = round(2.0 * thk, 3)      # grid covers the transition bands
    t0 = exc[0].centre_ms
    if kind == 'PRESS':
        te1 = 2 * (ref[0].centre_ms - t0)
        params.update({'Tau 1': round(te1, 4), 'Tau 2': round(te - te1, 4)})
    elif kind == 'STEAM':
        params['TM'] = round(exc[2].centre_ms - exc[1].centre_ms, 4)
    elif kind == 'MEGA-PRESS':
        del params['Flip Angle']                          # FID-A's MEGA-PRESS takes none
        params.update(_edit_params(path, info, workdir, bfield, te))
    else:
        # FID-A's semi-LASER places the refocusing pulses at TE/8, 3TE/8, 5TE/8, 7TE/8
        off = max(abs(r.centre_ms - t0 - k * te / 8) for r, k in zip(ref, (1, 3, 5, 7)))
        if off > 0.05:
            warnings.warn(f"{os.path.basename(path)}: refocusing pulses are up to {off:.2f} ms away "
                          "from the TE/8, 3TE/8, 5TE/8, 7TE/8 spacing FID-A's semi-LASER simulates.")
    return BACKENDS[kind], params


def _edit_params(path, info, workdir, bfield, te):
    """MEGA editing pulse, its duration and the ON / OFF frequencies (first two acquisitions)."""
    exc = next(r for r in info.rf if r.role == 'exc')
    edit = [r for r in info.rf if r.role == 'edit']
    ref = [r for r in info.rf if r.role == 'ref']
    params = {'Edit Pulse Path': waveform_pulse(edit[0], 'inv').write(
                  os.path.join(workdir, 'seq_editing.txt')),
              'Edit Tp': round(edit[0].dur_ms, 4)}
    # FID-A simulates its own timing: 90 - 5 - 180 - 17 - edit - 17 - 180 - 17 - edit (at TE 68)
    s = te / 68.0
    want = [5 * s, 22 * s, 39 * s, 56 * s]
    have = [r.centre_ms - exc.centre_ms for r in (ref[0], edit[0], ref[1], edit[1])]
    off = max(abs(a - b) for a, b in zip(have, want))
    if off > 0.05:
        warnings.warn(f"{os.path.basename(path)}: pulses are up to {off:.2f} ms away from the timing "
                      "FID-A's MEGA-PRESS simulates (its TE 68 ms layout scaled to this TE).")
    name = os.path.basename(path)
    if bfield in (None, ''):
        warnings.warn(f"{name}: the editing frequencies need B0, which a .seq file does not carry; "
                      "read the data file first. Edit On / Off keep the sheet's values.")
        return params
    try:
        other = [r for r in read_seq(path, 1).rf if r.role == 'edit']
    except ValueError:
        other = []
    ppm = [round(edit_ppm(edit[0], bfield), 3)] + ([round(edit_ppm(other[0], bfield), 3)] if other else [])
    if len(ppm) < 2 or abs(ppm[0] - ppm[1]) < 0.01:
        warnings.warn(f"{name}: the second acquisition does not edit at another frequency; "
                      "Edit Off keeps the sheet's value.")
        params['Edit On'] = ppm[0]
        return params
    upfield = [p < WATER_PPM for p in ppm]
    if upfield[0] == upfield[1]:
        warnings.warn(f"{name}: both editing frequencies ({ppm[0]:g}, {ppm[1]:g} ppm) lie on the same "
                      "side of water; the first acquisition is taken as edit ON.")
        on, off_ppm = ppm
    else:
        on, off_ppm = ppm if upfield[0] else ppm[::-1]    # ON edits upfield of water
    params.update({'Edit On': on, 'Edit Off': off_ppm})
    return params

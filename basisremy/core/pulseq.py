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
#          LASER). Echo times follow from the RF centres (the ADC may start before the echo).      #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import os
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


@dataclass
class SeqInfo:
    """What the sheet needs from one acquisition of a .seq file (the first ADC)."""
    rf: list[RFEvent]
    adc_start_ms: float
    samples: int
    dwell_s: float
    name: str | None


def read_seq(path: str) -> SeqInfo:
    """RF events up to the first ADC, and that ADC, from a Pulseq file (PyPulseq)."""
    import pypulseq as pp
    seq = pp.Sequence()
    seq.read(path)
    rf, t = [], 0.0
    for bid in sorted(seq.block_durations):
        b = seq.get_block(bid)
        if b.rf is not None:
            p = b.rf
            dt = float(np.mean(np.diff(p.t))) if len(p.t) > 1 else float(p.shape_dur)
            tc, _ = pp.calc_rf_center(p)
            t_rel = p.delay + tc                       # centre within the block [s]
            axis, g_c = None, None
            for ax in ('gx', 'gy', 'gz'):
                g = getattr(b, ax)
                if g is None:
                    continue
                if g.type == 'trap':
                    on = g.delay + g.rise_time <= t_rel <= g.delay + g.rise_time + g.flat_time
                    val = g.amplitude if on else 0.0
                else:
                    val = float(np.interp(t_rel, g.delay + np.asarray(g.tt), np.asarray(g.waveform),
                                          left=0.0, right=0.0))
                if abs(val) > abs(g_c or 0.0):
                    axis, g_c = ax[1], val
            rf.append(RFEvent(
                centre_ms=(t + t_rel) * 1e3, dur_ms=float(p.shape_dur) * 1e3,
                flip_deg=float(np.degrees(2 * np.pi * abs(np.sum(p.signal) * dt))),
                phase_rad=float(p.phase_offset), freq_offset_hz=float(p.freq_offset),
                signal=np.asarray(p.signal, dtype=complex), slab_axis=axis,
                slab_grad_hz_m=abs(g_c) if g_c else None))
        if b.adc is not None:
            return SeqInfo(rf=rf, adc_start_ms=(t + b.adc.delay) * 1e3, samples=int(b.adc.num_samples),
                           dwell_s=float(b.adc.dwell), name=seq.definitions.get('Name'))
        t += seq.block_durations[bid]
    raise ValueError(f"{os.path.basename(path)}: no ADC event in the sequence file.")


def sequence_type(info: SeqInfo) -> str:
    """PRESS / semi-LASER from the RF events before the readout."""
    exc = [r for r in info.rf if r.flip_deg < 120]
    ref = [r for r in info.rf if r.flip_deg >= 120]
    if len(exc) == 1 and len(ref) == 2:
        return 'PRESS'
    if len(exc) == 1 and len(ref) == 4:
        return 'sLASER'
    raise ValueError(f"Sequence with {len(exc)} excitation and {len(ref)} refocusing pulses "
                     "is not supported from a .seq file yet (PRESS and semi-LASER are).")


def echo_ms(info: SeqInfo) -> float:
    """Echo time from the RF centres: each refocusing pulse mirrors the previous echo."""
    exc = next(r for r in info.rf if r.flip_deg < 120)
    echo = exc.centre_ms
    for r in (r for r in info.rf if r.flip_deg >= 120):
        echo = 2 * r.centre_ms - echo
    return echo - exc.centre_ms


def waveform_pulse(event: RFEvent, kind: str = 'ref'):
    """The RF event as a pulse_library.Pulse (FID-A waveform convention)."""
    from basisremy.core.pulse_library import Pulse
    amp = np.abs(event.signal)
    wf = np.column_stack([np.degrees(np.angle(event.signal)) % 360.0, amp / amp.max(),
                          np.ones(len(amp))])
    return Pulse('seq', wf, kind, event.dur_ms)


def slab_cm(event: RFEvent) -> float | None:
    """Slab thickness: the pulse bandwidth (FID-A's measure) over the slice gradient."""
    if not event.slab_grad_hz_m:
        return None
    from basisremy.core.pulse_library import bandwidth_hz
    bw = bandwidth_hz(waveform_pulse(event), event.dur_ms, points=8001)
    return bw / event.slab_grad_hz_m * 100.0


def sheet_params(path: str, workdir: str) -> tuple[str, dict]:
    """(backend, registry-named parameters) for a FID-A shaped simulation of the .seq.

    Writes the refocusing waveform into workdir as an FID-A .txt file. TE, the PRESS
    echo split, duration, flip angle, slabs and the readout come from the file."""
    info = read_seq(path)
    kind = sequence_type(info)
    exc = next(r for r in info.rf if r.flip_deg < 120)
    ref = [r for r in info.rf if r.flip_deg >= 120]
    te = echo_ms(info)
    os.makedirs(workdir, exist_ok=True)
    pulse_file = waveform_pulse(ref[0]).write(os.path.join(workdir, 'seq_refocusing.txt'))
    slabs = [slab_cm(r) for r in (ref[0], ref[-1])]
    params = {
        'TE': round(te, 4),
        'Samples': info.samples,
        'Bandwidth': round(1.0 / info.dwell_s, 4),
        'Path to Pulse': pulse_file,
        'RefTp': round(ref[0].dur_ms, 4),
        'Flip Angle': round(ref[0].flip_deg, 2),
    }
    for axis, thk in zip(('X', 'Y'), slabs):
        if thk:
            params[f'thk{axis}'] = round(thk, 3)
            params[f'fov{axis}'] = round(2.0 * thk, 3)      # grid covers the transition bands
    if kind == 'PRESS':
        te1 = 2 * (ref[0].centre_ms - exc.centre_ms)
        params.update({'Tau 1': round(te1, 4), 'Tau 2': round(te - te1, 4)})
        return 'FidaPressShaped', params
    return 'FidaSemiLaserShaped', params

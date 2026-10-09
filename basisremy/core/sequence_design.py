####################################################################################################
#                                       sequence_design.py                                         #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 07/10/26                                                                                #
#                                                                                                  #
# Purpose: One sequence design for every engine: the sequence, its timings, a pulse per role       #
#          (ideal, a standard shape or a pulse file in any format, with its duration), the         #
#          editing targets and the voxel. The designer behind the wand edits one; it is saved      #
#          as a Pulseq .seq (RF `use` roles, [DEFINITIONS]) and read back from any .seq or         #
#          FSL-MRS / WIN sequence JSON. Per engine, plan() tells whether and how the engine runs   #
#          the design (and why not), apply() fills its sheet, and fsl_sequences() translates a     #
#          .seq into FSL-MRS sequence descriptions (one per sub-experiment).                       #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field

import numpy as np

from basisremy.core import sequence_setup as ss

DESIGNABLE = ['PRESS', 'sLASER', 'STEAM', 'Spin Echo', 'LASER', 'MEGA-PRESS', 'MEGA-sLASER',
              'MEGA-SPECIAL', 'HERMES', 'HERCULES', 'HERMES (sLASER)', 'HERCULES (sLASER)']
# Hadamard-edited sequences: the MEGA layout, four sub-experiments A-D
_LAYOUT = {'HERMES': 'MEGA-PRESS', 'HERCULES': 'MEGA-PRESS',
           'HERMES (sLASER)': 'MEGA-sLASER', 'HERCULES (sLASER)': 'MEGA-sLASER'}
# sub-experiment -> editing targets [ppm] (two = dual-lobe pulse), in MRSCloud's order A-D, as
# the FSL-MRS backend simulates them
SCHEMES = {'HERMES':   {'A': (4.56,), 'B': (1.90,), 'C': (4.56, 1.90), 'D': (7.50,)},
           'HERCULES': {'A': (4.58,), 'B': (4.18,), 'C': (4.58, 1.90), 'D': (4.18, 1.90)}}
_SCHEME_SRC = {'HERMES': "HERMES GABA / GSH (Chan et al. 2016), as MRSCloud orders A-D",
               'HERCULES': "HERCULES (Oeltzschner et al. 2019), as MRSCloud orders A-D"}
_AXES = 'xyz'                      # voxel: x = excitation (L-R), y / z = refocusing (A-P / C-C)
WATER_PPM = 4.65                   # carrier on water
GAMMA_HZ_T = 42.577e6
IDEAL_MS = 0.02                    # an ideal pulse in a .seq: 20 us, constant, non-selective
_RASTER = 1e-5                     # Pulseq block / gradient raster [s]
_KIND = {'exc': 'exc', 'ref': 'ref', 'edit': 'inv'}       # role -> how B1 is scaled
_USE = {'exc': 'excitation', 'ref': 'refocusing', 'edit': 'inversion'}
_FLIP = {'exc': 90.0, 'ref': 180.0, 'edit': 180.0}
GRADIENTS = {'max_mT_m': 80.0, 'rise_ms': 0.4}        # clinical system (200 T/m/s) without header values
_SALEH = "Saleh et al. 2019, multi-vendor universal MEGA-PRESS"
# standard pulses offered per role
STANDARD_FOR = {'exc': ['sinc-exc'],
                'ref': ['sinc-ref', 'hs4-ref', 'hs1-inv', 'goia-wurst', 'goia-hs', 'foci'],
                'edit': ['gauss-edit', 'gauss-edit-20ms']}
# timing values per sequence (ms); TE itself is the scan's
_TIMING_KEYS = {'PRESS': ('TE1', 'TE2'), 'MEGA-PRESS': ('TE1', 'TE2'),
                'sLASER': ('TE1', 'TE2', 'TE3'), 'MEGA-sLASER': ('TE1', 'TE2', 'TE3'),
                'STEAM': ('TM',), 'Spin Echo': (), 'LASER': (), 'MEGA-SPECIAL': ()}
TIMING_KEYS = {k: _TIMING_KEYS[_LAYOUT.get(k, k)] for k in DESIGNABLE}
TIMING_LABEL = {'TE1': 'TE1 (first echo)', 'TE2': 'TE2 (second echo)', 'TE3': 'TE3 (last pair)',
                'TM': 'TM (mixing time)'}


def roles(kind: str | None) -> list[str]:
    return ss.roles(kind)


def _edited(kind) -> bool:
    return 'edit' in roles(kind)


def layout(kind: str) -> str:
    """The pulse layout of ``kind``: HERMES / HERCULES use the MEGA ones."""
    return _LAYOUT.get(kind, kind)


def _family(kind: str) -> str | None:
    return next((f for f in SCHEMES if kind.startswith(f)), None)


#**************************************************************************************************#
#                                             design                                               #
#**************************************************************************************************#
@dataclass
class Design:
    """kind, TE [ms], timing {TE1/TE2/TE3/TM: ms}, pulses {role: {'source', 'dur'}} where source
    is 'ideal', 'standard:<name>' or a pulse file (any format, '#role' / '#n' selects in a whole
    sequence file) and dur its duration in ms; edit = (ON, OFF) ppm (MEGA); scheme = {sub-
    experiment: editing targets in ppm} (HERMES / HERCULES); voxel = (x, y, z) cm the selective
    pulses select (excitation x, first refocusing y, second z); bfield [T] for dual-lobe editing;
    gradients {'max_mT_m', 'rise_ms'} the scanner's gradient system where the data header gives it
    (else GRADIENTS); rec {key: source} the values that are recommendations."""
    kind: str
    te: float
    timing: dict = field(default_factory=dict)
    pulses: dict = field(default_factory=dict)
    edit: tuple = (1.9, 7.5)
    scheme: dict = field(default_factory=dict)
    voxel: tuple = (2.0, 2.0, 2.0)
    bfield: float | None = None
    samples: int | None = None
    bandwidth: float | None = None
    gradients: dict | None = None
    rec: dict = field(default_factory=dict)


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def default_duration(source: str, role: str) -> tuple[float, str]:
    """(ms, where it comes from) for a pulse source: its own design duration where it has one."""
    if source == 'ideal':
        return 0.0, ''
    try:
        from basisremy.core.sequence_view import read_pulse
        tp = read_pulse(source, _KIND[role]).tp_ms
    except Exception:                                      # noqa: BLE001 - unreadable: typical value
        tp = float('nan')
    if tp == tp and tp > 0:
        return round(tp, 4), "the pulse's own duration"
    typical = {'exc': 2.0, 'ref': 5.0, 'edit': 15.0}[role]
    return typical, "typical duration: the file does not hold one"


_VOXEL_KEYS = ('LeftRightSize', 'AnteriorPosteriorSize', 'CranioCaudalSize')   # mm, from REMY


def _header_pulses(d: Design, kind: str, header: dict | None):
    """The scan's own pulse per role where the data header stores it (design_fields.PULSE_PARAMS,
    Bruker method): its waveform, written as a Bruker shape file, at its duration; without a
    waveform the pulse stays ideal and the note gives the header's values."""
    from basisremy.remy.design_fields import PULSE_PARAMS
    found = {}
    for name, p in ((header or {}).get(PULSE_PARAMS) or {}).items():
        if p.get('role') in roles(kind) and p['role'] not in found:
            found[p['role']] = (name, p)
    for role, (name, p) in found.items():
        info = f"{p['dur_ms']:g} ms, {p['bw_hz']:g} Hz, {p['flip']:g} deg, shape '{p['shape']}'"
        if p.get('waveform'):
            d.pulses[role] = {'source': _shape_file(name, p['waveform']), 'dur': p['dur_ms']}
            d.rec.pop(f'pulse:{role}', None)
        else:
            d.rec[f'pulse:{role}'] = (f"Default: ideal. The data header gives {name}: {info}, but not its "
                                      "waveform. Most precise: the scanner's own pulse file")


def _shape_file(name: str, waveform) -> str:
    """A header waveform [(amplitude %, phase deg)] as a Bruker JCAMP shape in designs_dir()/pulses
    (named by its content, so the same pulse is written once)."""
    import hashlib
    text = ''.join(f'{a:.6e}, {ph:.6e}\n' for a, ph in waveform)
    folder = os.path.join(designs_dir(), 'pulses')
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, f"header_{name}_{hashlib.sha1(text.encode()).hexdigest()[:8]}.exc")
    if not os.path.exists(path):
        with open(path, 'w') as f:
            f.write(f'##TITLE= {name} from the data header\n##JCAMP-DX= 5.00 Bruker JCAMP library\n'
                    f'##DATA TYPE= Shape Data\n##NPOINTS= {len(waveform)}\n##XYPOINTS= (XY..XY)\n'
                    f'{text}##END=\n')
    return path


def recommend(kind: str, te: float, sheet: dict | None = None, header: dict | None = None) -> Design:
    """The recommended design for ``kind`` at ``te``: ideal excitation and refocusing, a shaped
    editing pulse (Saleh 2019), symmetric timing; sheet values (TE1/TE2, sLASER TE1-3, TM) that
    the user set or the file gave are kept, and the voxel comes from the data ``header`` (REMY's
    L-R / A-P / C-C sizes) when it has one."""
    sheet = sheet or {}
    d = Design(kind=kind, te=float(te))
    lay = layout(kind)
    sym = "Symmetric: the scan's own split is not in the data header"
    keys = TIMING_KEYS[kind]
    given = {'TE1': _num(sheet.get('Tau 1')), 'TE2': _num(sheet.get('Tau 2') or sheet.get('TE2'))}
    if lay in ('sLASER', 'MEGA-sLASER'):
        given = {f'TE{i}': _num(sheet.get(f'sLASER TE{i}')) for i in (1, 2, 3)}
    given['TM'] = _num(sheet.get('TM'))
    if kind == 'PRESS':
        d.timing, d.rec = {'TE1': te / 2, 'TE2': te / 2}, {'TE1': sym, 'TE2': sym}
    elif lay == 'MEGA-PRESS':
        te1 = round(13.1 * te / 68.0, 3)
        src = f"{_SALEH}: TE1 13.1 ms at TE 68, scaled to TE"
        d.timing, d.rec = {'TE1': te1, 'TE2': te - te1}, {'TE1': src, 'TE2': src}
    elif lay in ('sLASER', 'MEGA-sLASER'):
        src = "Symmetric spacing (TE/4, TE/2, TE/4): the scan's own is not in the data header"
        d.timing = {'TE1': te / 4, 'TE2': te / 2, 'TE3': te / 4}
        d.rec = dict.fromkeys(d.timing, src)
    elif kind == 'STEAM':
        d.timing, d.rec = {'TM': 10.0}, {'TM': "Short mixing time commonly used; this data header holds no TM"}
    if all(given.get(k) is not None for k in keys) and keys:
        if kind == 'STEAM' or abs(sum(given[k] for k in keys) - te) < 1e-3:
            d.timing = {k: given[k] for k in keys}
            for k in keys:
                d.rec.pop(k, None)
    for role in roles(kind):
        d.pulses[role] = {'source': 'ideal', 'dur': 0.0}
        d.rec[f'pulse:{role}'] = "Default: ideal. Most precise: the scanner's own pulse file"
    _header_pulses(d, kind, header)
    from basisremy.remy.design_fields import GRADIENTS as HEADER_GRADIENTS
    d.gradients = (header or {}).get(HEADER_GRADIENTS)
    if _edited(kind):
        family = _family(kind)
        if family:
            tp, src = 20.0, f"20 ms editing pulses, as in HERMES ({_SALEH})"
        else:
            tp, src = 15.0, f"15 ms editing pulse ({_SALEH})"
        fit = te / 4 - 2 * IDEAL_MS                   # between the (ideal) pulses of a pair
        if lay == 'MEGA-sLASER' and fit < tp:
            tp, src = round(math.floor(fit * 100) / 100, 2), src + ", shortened to fit TE"
        d.pulses['edit'] = {'source': 'standard:gauss-edit-20ms' if family else 'standard:gauss-edit',
                            'dur': tp}
        d.rec['pulse:edit'] = (f"Gaussian editing pulse, open stand-in for the sinc-Gaussian "
                               f"of {_SALEH}; an editing pulse is never ideal")
        d.rec['dur:edit'] = src
        if family:
            d.scheme = dict(SCHEMES[family])
            d.rec['edit'] = _SCHEME_SRC[family]
        else:
            d.rec['edit'] = f"GABA editing: ON 1.9 ppm, OFF 7.5 ppm ({_SALEH})"
    mm = [_num((header or {}).get(k)) for k in _VOXEL_KEYS]
    if all(mm):
        d.voxel = tuple(round(v / 10.0, 4) for v in mm)
    else:
        d.rec['voxel'] = "2 x 2 x 2 cm: the data header holds no voxel size"
    d.bfield = _num(sheet.get('Bfield'))
    d.samples = int(_num(sheet.get('Samples'))) if _num(sheet.get('Samples')) else None
    d.bandwidth = _num(sheet.get('Bandwidth'))
    return d


#**************************************************************************************************#
#                                            timeline                                              #
#**************************************************************************************************#
def events(d: Design) -> list[dict]:
    """[{'role', 'centre', 'dur', 'source', 'axis'}] in ms from the excitation's centre."""
    te, t, kind = d.te, d.timing, layout(d.kind)
    ev = [('exc', 0.0, 'x')]
    if kind in ('PRESS', 'MEGA-PRESS'):
        te1, te2 = t['TE1'], t['TE2']
        ev += [('ref', te1 / 2, 'y'), ('ref', te1 + te2 / 2, 'z')]
        if kind == 'MEGA-PRESS':        # symmetric about the second refocusing pulse
            ev += [('edit', te1 + te2 / 4, None), ('edit', te1 + 3 * te2 / 4, None)]
    elif kind in ('Spin Echo', 'MEGA-SPECIAL'):
        ev += [('ref', te / 2, 'y')]
        # MEGA-SPECIAL: the spin echo the ISIS add / subtract leaves (the inversion itself is not
        # part of the design), editing pulses TE/4 apart as in FID-A's run_simMegaSpecialShaped
        if kind == 'MEGA-SPECIAL':
            ev += [('edit', te / 4, None), ('edit', 3 * te / 4, None)]
    elif kind in ('sLASER', 'MEGA-sLASER'):
        te1, te2, te3 = t['TE1'], t['TE2'], t['TE3']
        c = [te1 / 2, te1 + te2 / 4, te1 + 3 * te2 / 4, te1 + te2 + te3 / 2]
        ev += [('ref', x, a) for x, a in zip(c, 'yyzz')]
        if kind == 'MEGA-sLASER':       # between the pulses of each pair
            ev += [('edit', (c[0] + c[1]) / 2, None), ('edit', (c[2] + c[3]) / 2, None)]
    elif kind == 'LASER':
        ev += [('ref', te * (2 * i + 1) / 12, a) for i, a in enumerate('xxyyzz')]
    elif kind == 'STEAM':
        ev = [('exc', 0.0, 'x'), ('exc', te / 2, 'y'), ('exc', te / 2 + t['TM'], 'z')]
    out = []
    for role, c, axis in sorted(ev, key=lambda e: e[1]):
        p = d.pulses[role]
        ideal = p['source'] == 'ideal'
        out.append({'role': role, 'centre': c, 'dur': 0.0 if ideal else float(p['dur']),
                    'source': None if ideal else p['source'], 'axis': None if ideal else axis})
    return out


def acquisitions(d: Design) -> list[tuple]:
    """[(label, editing targets in ppm)] per acquisition: MEGA ON / OFF, HERMES / HERCULES
    A-D, else one unedited acquisition."""
    if d.scheme:
        return [(k, tuple(v)) for k, v in d.scheme.items()]
    if _edited(d.kind):
        return [('ON', (d.edit[0],)), ('OFF', (d.edit[1],))]
    return [(None, ())]


def echo(d: Design) -> float:
    return d.te + (d.timing.get('TM', 0.0) if d.kind == 'STEAM' else 0.0)


def problems(d: Design) -> list[str]:
    """What keeps the design from being saved: timings that do not add up, overlapping pulses."""
    out = []
    keys = TIMING_KEYS[d.kind]
    if any(_num(d.timing.get(k)) is None for k in keys):
        return ["Set every timing value."]
    if d.kind != 'STEAM' and keys and abs(sum(d.timing[k] for k in keys) - d.te) > 1e-3:
        out.append(f"{' + '.join(keys)} = {sum(d.timing[k] for k in keys):g} ms must equal TE = {d.te:g} ms.")
    for role, p in d.pulses.items():
        if p['source'] != 'ideal' and not (_num(p['dur']) or 0) > 0:
            out.append(f"Set the {ss.ROLE_NAME[role].lower()} pulse's duration.")
    if any(len(t) > 1 for t in d.scheme.values()) and not d.bfield:
        out.append("A dual-lobe editing pulse needs the field strength: read the data file or "
                   "type it in the sheet.")
    if out:
        return out
    ev = events(d)
    for a, b in zip(ev, ev[1:]):
        gap = (b['centre'] - b['dur'] / 2) - (a['centre'] + a['dur'] / 2)
        if gap < 0:
            out.append(f"The {ss.ROLE_NAME[a['role']].lower()} and {ss.ROLE_NAME[b['role']].lower()} "
                       f"pulses overlap by {-gap:.2f} ms: shorten them or lengthen TE.")
    last = ev[-1]
    if last['centre'] + last['dur'] / 2 > echo(d) + 1e-6:
        out.append("The last pulse runs past the echo: shorten it or lengthen TE.")
    return out


#**************************************************************************************************#
#                                          Pulseq writer                                           #
#**************************************************************************************************#
def _ceil(t: float) -> float:
    return math.ceil(round(t / _RASTER, 6)) * _RASTER


def _waveform(source: str, role: str, dur_ms: float):
    """(complex B1 [Hz] on the 1 us RF raster, gradient shape [G/cm] or None, bandwidth [Hz])."""
    from basisremy.core.pulse_library import bandwidth_hz
    from basisremy.core.rf_pulses import scale_waveform
    from basisremy.core.sequence_view import read_pulse
    pulse = read_pulse(source, _KIND[role])
    phase, amp, dt = scale_waveform(pulse.waveform, dur_ms / 1e3, _KIND[role])
    edges = np.concatenate([[0.0], np.cumsum(dt)])
    mid = (edges[:-1] + edges[1:]) / 2
    n = int(round(dur_ms * 1e3))
    t = (np.arange(n) + 0.5) * 1e-6 * (edges[-1] / (n * 1e-6))
    b1 = amp * np.exp(1j * np.deg2rad(phase))          # complex: sign flips pass through zero
    signal = np.interp(t, mid, b1.real) + 1j * np.interp(t, mid, b1.imag)
    grad = None
    if pulse.is_gradient_modulated:
        grad = np.interp(t, mid, pulse.waveform[:, 3]) * (pulse.thk_cm or 1.0)
    bw = None if grad is not None or role == 'edit' else bandwidth_hz(pulse, dur_ms, points=8001)
    return signal, grad, bw


def write_seq(d: Design, path: str) -> str:
    """Write the design as a Pulseq .seq: one acquisition (MEGA: ON, then OFF) of the RF pulses
    with their slice gradients, the excitation rephasers, and the readout starting at the echo."""
    import pypulseq as pp
    errs = problems(d)
    if errs:
        raise ValueError(" ".join(errs))
    g = d.gradients or GRADIENTS
    system = pp.Opts(max_grad=g['max_mT_m'], grad_unit='mT/m', rise_time=g['rise_ms'] / 1e3,
                     adc_raster_time=1e-9)
    seq = pp.Sequence(system)
    ev = events(d)
    ideal = [r for r, p in d.pulses.items() if p['source'] == 'ideal']
    built = {}

    def rf_block(e, targets, shift):
        """RF event and gradient of one pulse, the RF delayed by ``shift`` (< one block raster)
        so its centre lands on the design time; the centre lies ``lead`` s into the block. An
        editing pulse sits at its one target (freq_ppm) or carries several as lobes."""
        role = e['role']
        phase = math.pi / 2 if role == 'ref' else 0.0
        freq_ppm = (targets[0] - WATER_PPM) if role == 'edit' and len(targets) == 1 else 0.0
        if e['source'] is None:
            n = int(round(IDEAL_MS * 1e3))
            b1 = _FLIP[role] / 360.0 / (IDEAL_MS / 1e3)       # Hz: flip = 2 pi B1 t
            rf = pp.make_arbitrary_rf(np.full(n, b1, dtype=complex), math.radians(_FLIP[role]),
                                      no_signal_scaling=True, use=_USE[role], phase_offset=phase,
                                      freq_ppm=freq_ppm, delay=shift, system=system)
            return [rf], rf.delay + IDEAL_MS / 2e3, None
        key = (e['source'], role, e['dur'])
        if key not in built:
            built[key] = _waveform(e['source'], role, e['dur'])
        signal, gshape, bw = built[key]
        dur = len(signal) * 1e-6
        if role == 'edit' and len(targets) > 1:               # dual-lobe: one lobe per target
            t = (np.arange(len(signal)) + 0.5) * 1e-6 - dur / 2
            hz = [(p - WATER_PPM) * float(d.bfield) * 42.577 for p in targets]
            signal = signal * sum(np.exp(2j * np.pi * f * t) for f in hz)
        if e['axis'] is None:                                 # non-selective (editing)
            rf = pp.make_arbitrary_rf(signal, math.radians(_FLIP[role]), no_signal_scaling=True,
                                      use=_USE[role], phase_offset=phase, freq_ppm=freq_ppm,
                                      delay=shift, system=system)
            return [rf], rf.delay + dur / 2, None
        slab_m = d.voxel[_AXES.index(e['axis'])] / 100.0
        if gshape is not None:                               # gradient-modulated (GOIA / FOCI)
            ramp = 0.2e-3
            g = gshape / slab_m * 1e-2 * 4257.7 * 100.0       # G/cm at the design slab -> Hz/m
            n_ramp = int(round(ramp / _RASTER))
            gs = g[::10][:int(round(dur / _RASTER))]
            wave = np.concatenate([np.linspace(0, gs[0], n_ramp + 1)[:-1], gs,
                                   np.linspace(gs[-1], 0, n_ramp + 1)[1:]])
            shift = 0.0               # RF and gradient stay in step (centre within one raster)
            grad = pp.make_arbitrary_grad(e['axis'], wave, system=system, first=0, last=0)
        else:
            grad = pp.make_trapezoid(e['axis'], amplitude=bw / slab_m, flat_time=_ceil(dur + shift),
                                     system=system)
            ramp = grad.rise_time
        rf = pp.make_arbitrary_rf(signal, math.radians(_FLIP[role]), no_signal_scaling=True,
                                  delay=ramp + shift, use=_USE[role], phase_offset=phase,
                                  freq_ppm=freq_ppm, system=system)
        rephase = None
        if role == 'exc' and gshape is None:                  # area after the centre, undone
            after = grad.flat_time - shift - dur / 2
            rephase = pp.make_trapezoid(e['axis'], area=-(grad.amplitude * after
                                                          + grad.amplitude * grad.fall_time / 2),
                                        system=system)
        return [rf, grad], ramp + shift + dur / 2, rephase

    def add(*events_):
        seq.add_block(*events_, pp.make_delay(_ceil(pp.calc_duration(*events_))))

    def acquisition(ppm):
        """One acquisition with the editing targets ``ppm``."""
        t = 0.0                                                # time in this acquisition [s]
        t0 = None
        for i, e in enumerate(ev):
            _, lead, _ = rf_block(e, ppm, 0.0)
            start = (e['centre'] / 1e3) - lead
            if t0 is None:
                t0 = -start                                    # first pulse starts the acquisition
                t = 0.0
                shift = 0.0
            else:
                ahead = start + t0 - t
                gap = math.floor(round(ahead / _RASTER, 6)) * _RASTER
                shift = round((ahead - gap) * 1e6) * 1e-6          # on the 1 us RF raster
                if gap < -1e-9:
                    raise ValueError(f"Pulse {i + 1} does not fit the timing (overlaps by "
                                     f"{-gap * 1e3:.2f} ms with the previous block).")
                if gap > 0:
                    add(pp.make_delay(gap))
                    t += gap
            blocks, _, rephase = rf_block(e, ppm, shift)
            add(*blocks)
            t += _ceil(pp.calc_duration(*blocks))
            if rephase is not None:
                add(rephase)
                t += _ceil(pp.calc_duration(rephase))
        gap = round((echo(d) / 1e3 + t0 - t) / _RASTER) * _RASTER
        if gap < -1e-9:
            raise ValueError("The last pulse runs past the echo.")
        if gap > 0:
            add(pp.make_delay(gap))
        samples = int(d.samples or 2048)
        dwell = 1.0 / float(d.bandwidth or 4000.0)
        add(pp.make_adc(samples, dwell=dwell, system=system))
        add(pp.make_delay(0.01))                               # end of the acquisition

    for _label, targets in acquisitions(d):
        acquisition(targets)
    seq.set_definition('Name', d.kind)
    seq.set_definition('Creator', 'BasisREMY')
    seq.set_definition('TE', d.te)
    for k, v in d.timing.items():
        seq.set_definition(k, round(float(v), 6))
    if ideal:
        seq.set_definition('IdealPulses', ' '.join(ideal))
    if _edited(d.kind):
        seq.set_definition('EditScheme', ';'.join(f"{k}={','.join(f'{p:g}' for p in v)}"
                                                  for k, v in acquisitions(d)))
    if d.bfield:
        seq.set_definition('B0', d.bfield)
    seq.set_definition('VoxelCm', ' '.join(f'{v:g}' for v in d.voxel))
    if d.gradients:
        seq.set_definition('GradientSystem', f"{d.gradients['max_mT_m']:g} {d.gradients['rise_ms']:g}")
    seq.set_definition('ADCFromScan', int(bool(d.samples and d.bandwidth)))
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    seq.write(path)
    return path


#**************************************************************************************************#
#                                     reading a sequence file                                      #
#**************************************************************************************************#
_COUNTS = {(1, 2, 0): 'PRESS', (1, 4, 0): 'sLASER', (3, 0, 0): 'STEAM', (1, 1, 0): 'Spin Echo',
           (1, 6, 0): 'LASER', (1, 2, 2): 'MEGA-PRESS', (1, 4, 2): 'MEGA-sLASER',
           (1, 1, 2): 'MEGA-SPECIAL'}


def kind_of(roles_: list[str]) -> str:
    n = tuple(sum(r == k for r in roles_) for k in ('exc', 'ref', 'edit'))
    if n in _COUNTS:
        return _COUNTS[n]
    raise ValueError(f"A sequence with {n[0]} excitation, {n[1]} refocusing and {n[2]} editing "
                     f"pulses is not one BasisREMY simulates (PRESS, sLASER, STEAM, Spin Echo, "
                     f"LASER, MEGA-PRESS, MEGA-sLASER, MEGA-SPECIAL).")


def _timing_from(kind, centres, te):
    """TE1/TE2/TE3/TM from the pulse centres (ms from the excitation)."""
    kind = layout(kind)
    exc = [c for r, c in centres if r == 'exc']
    ref = [c for r, c in centres if r == 'ref']
    if kind == 'STEAM':
        return {'TM': exc[2] - exc[1]}
    if kind in ('PRESS', 'MEGA-PRESS'):
        return {'TE1': 2 * ref[0], 'TE2': te - 2 * ref[0]}
    if kind in ('sLASER', 'MEGA-sLASER'):
        e = [0.0]
        for c in ref:
            e.append(2 * c - e[-1])
        return {'TE1': e[1], 'TE2': e[3] - e[1], 'TE3': te - e[3]}
    return {}


def read_design(path: str, bfield=None) -> Design:
    """The design a whole-sequence file holds (Pulseq .seq, or FSL-MRS / WIN sequence JSON)."""
    if path.lower().endswith('.seq'):
        return _read_seq_design(path, bfield)
    return _read_json_design(path)


def _read_seq_design(path, bfield):
    import pypulseq as pp
    from basisremy.core import pulseq
    info = pulseq.read_seq(path)
    raw = pp.Sequence()
    raw.read(path)
    defs = raw.definitions
    rf = [r for r in info.rf if r.role in ('exc', 'ref', 'edit')]
    kind = kind_of([r.role for r in rf])
    name = str(defs.get('Name', ''))
    if name in DESIGNABLE and layout(name) == kind:          # HERMES / HERCULES by name
        kind = name
    t0 = next(r.centre_ms for r in rf if r.role == 'exc')
    te = pulseq.echo_ms(info)
    centres = [(r.role, r.centre_ms - t0) for r in rf]
    ideal = set(str(defs.get('IdealPulses', '')).split())
    bfield = _num(defs.get('B0')) or _num(bfield)
    d = Design(kind=kind, te=round(float(te), 4),
               timing={k: round(float(v), 4) for k, v in _timing_from(kind, centres, te).items()})
    for role in roles(kind):
        first = next(r for r in rf if r.role == role)
        if role in ideal or first.dur_ms <= 0.1:
            d.pulses[role] = {'source': 'ideal', 'dur': 0.0}
        else:
            d.pulses[role] = {'source': f'{path}#{role}', 'dur': round(first.dur_ms, 4)}
    scheme = {}
    for part in str(defs.get('EditScheme', '')).split(';'):
        label, _, ppm = part.partition('=')
        if ppm:
            scheme[label.strip()] = tuple(float(p) for p in ppm.split(','))
    if _family(kind):
        d.scheme = scheme or dict(SCHEMES[_family(kind)])
    elif _edited(kind) and {'ON', 'OFF'} <= set(scheme):
        d.edit = (scheme['ON'][0], scheme['OFF'][0])
    elif _edited(kind):
        ppm = []
        for acq in (0, 1):
            try:
                ev = next(r for r in pulseq.read_seq(path, acq).rf if r.role == 'edit')
            except (ValueError, StopIteration):
                break
            if ev.freq_offset_hz and not bfield:
                break
            ppm.append(round(WATER_PPM + ev.freq_ppm + (ev.freq_offset_hz / (float(bfield) * 42.577)
                                                       if ev.freq_offset_hz else 0.0), 3))
        if len(ppm) == 2:
            d.edit = tuple(ppm) if ppm[0] < WATER_PPM else (ppm[1], ppm[0])
    voxel = [_num(v) for v in np.ravel(defs.get('VoxelCm', []))]       # pypulseq: an array
    if len(voxel) == 3 and all(voxel):
        d.voxel = tuple(voxel)
    else:              # from the slice gradients: excitation, first and last refocusing slab
        exc = [r for r in rf if r.role == 'exc']
        ref = [r for r in rf if r.role == 'ref']
        picks = [exc[0], ref[0] if ref else None, ref[-1] if len(ref) > 1 else None]
        if kind == 'STEAM':
            picks = exc[:3]
        slabs = [pulseq.slab_cm(r) if r is not None and r.slab_grad_hz_m else None for r in picks]
        d.voxel = tuple(round(v, 3) if v else 2.0 for v in slabs)
    d.bfield = bfield
    grad = [_num(v) for v in np.ravel(defs.get('GradientSystem', []))]   # pypulseq: an array
    if len(grad) == 2 and all(grad):
        d.gradients = {'max_mT_m': grad[0], 'rise_ms': grad[1]}
    if int(_num(defs.get('ADCFromScan')) if defs.get('ADCFromScan') is not None else 1):
        d.samples, d.bandwidth = info.samples, round(1.0 / info.dwell_s, 4)
    return d


def _read_json_design(path):
    from basisremy.core import sequence_view as sv
    seq = sv._fsl_sequence(path)
    roles_ = sv._fsl_roles(seq)
    tl = sv._seq_timeline(seq, path)
    kind = kind_of(roles_)
    centres = [(e['role'], e['centre_ms']) for e in tl['events']]
    te = tl['echo_ms'] - (centres[2][1] - centres[1][1] if kind == 'STEAM' else 0.0)
    d = Design(kind=kind, te=round(float(te), 4),
               timing={k: round(float(v), 4) for k, v in _timing_from(kind, centres, te).items()})
    for role in roles(kind):
        i = roles_.index(role)
        dur = float(seq['RF'][i]['time']) * 1e3
        d.pulses[role] = ({'source': 'ideal', 'dur': 0.0} if dur <= 0.1
                          else {'source': f'{path}#{role}', 'dur': round(dur, 4)})
    d.samples, d.bandwidth = seq.get('Rx_Points'), seq.get('Rx_SW')
    d.bfield = _num(seq.get('B0'))
    return d


#**************************************************************************************************#
#                                        engines: plan, apply                                      #
#**************************************************************************************************#
@dataclass
class Plan:
    """How one engine runs a design: 'ok' (as designed), 'approx' (runs, with the differences in
    notes) or 'no' (cannot; notes say why)."""
    status: str
    notes: list = field(default_factory=list)
    route: object = None
    values: dict = field(default_factory=dict)       # what apply() wrote into the sheet


def _symmetric(d: Design) -> bool:
    t, te = d.timing, d.te
    if d.kind in ('PRESS',):
        return abs(t['TE1'] - te / 2) < 0.01
    if layout(d.kind) in ('sLASER', 'MEGA-sLASER'):
        return all(abs(t[k] - te * f) < 0.01 for k, f in (('TE1', .25), ('TE2', .5), ('TE3', .25)))
    return True


def _timing_keys(backend, route) -> set:
    """Timing fields ``backend`` shows once ``route`` selects the sequence (sheet restored)."""
    saved = (dict(backend.mandatory_params), dict(backend.optional_params),
             getattr(backend, 'current_mode', None))
    try:
        if route.modes and backend.current_mode not in route.modes:
            backend.set_mode(route.modes[0])
        for k, v in route.sheet.items():
            (backend.mandatory_params if k in backend.mandatory_params
             or k not in backend.optional_params else backend.optional_params)[k] = v
        shown = backend.get_params_for_mode()
    finally:
        backend.mandatory_params.clear()
        backend.mandatory_params.update(saved[0])
        backend.optional_params.clear()
        backend.optional_params.update(saved[1])
        if saved[2] is not None and backend.current_mode != saved[2]:
            backend.set_mode(saved[2])
    return {k for k in ('Tau 1', 'Tau 2', 'TE2', 'sLASER TE1', 'sLASER TE2', 'sLASER TE3', 'TM')
            if k in shown}


def plan(d: Design, category: str, br=None, whole_file: str | None = None) -> Plan:
    """Whether ``category`` runs the design. FSL-MRS runs a whole sequence file as it is; the
    other engines take the pulses one by one through their own routes (core/sequence_setup)."""
    eng = ss.ENGINE_LABEL.get(category, category)
    if category == 'MRSCloud':
        if d.kind not in ('PRESS', 'MEGA-PRESS'):
            return Plan('no', [f"MRSCloud runs a design for PRESS and MEGA-PRESS; its {d.kind} keeps "
                               f"its own pulse set."])
        ref = d.pulses['ref']['source']
        if ref != 'ideal':
            from basisremy.core.sequence_view import read_pulse
            if read_pulse(ref, 'ref').is_gradient_modulated:
                return Plan('no', ["MRSCloud's PRESS takes no gradient-modulated refocusing pulse."])
        status, notes = 'ok', []
        if d.pulses['exc']['source'] != 'ideal':
            status = 'approx'
            notes.append("MRSCloud's excitation is ideal: the design's excitation pulse is not used.")
        if d.kind == 'MEGA-PRESS':
            status = 'approx'
            notes.append("MRSCloud places the editing pulses its own way (TE1 as designed).")
        return Plan(status, notes, ss.ROUTES['MRSCloud'][d.kind][0])
    if category == 'FSL-MRS':
        if d.kind not in DESIGNABLE:
            return Plan('no', [f"FSL-MRS: {d.kind} cannot be described here."])
        return Plan('ok', ["FSL-MRS runs the whole sequence: every pulse, the timing and the slabs."])
    routes = ss.ROUTES.get(category, {}).get(d.kind)
    if not routes:
        return Plan('no', [ss.why_not(category, d.kind)])
    want = {r: ss.IDEAL if p['source'] == 'ideal' else ss.SHAPED for r, p in d.pulses.items()}
    route = ss.choose(category, d.kind, want)
    status, notes = 'ok', []
    for role, kind in want.items():
        got = route.pulses.get(role)
        if got == kind:
            continue
        name = ss.ROLE_NAME[role].lower()
        if got == ss.OWN:
            status = 'approx'
            notes.append(f"{eng} uses its own {name} pulse (Gaussian), not the design's.")
        elif role == 'exc' and kind == ss.SHAPED:
            status = 'approx'
            notes.append(f"{eng} has no waveform excitation for {d.kind}: the excitation is "
                         f"simulated ideal.")
        else:
            return Plan('no', [ss.pulse_note(category, d.kind, role, kind)], route)
    # timing the engine cannot follow
    fixed = None
    if br is not None and br.backends.get(route.backend) is not None:
        b = br.backends[route.backend]
        keys = _timing_keys(b, route)
        need = {'PRESS': {'Tau 1', 'TE2'}, 'MEGA-PRESS': {'Tau 1'}, 'sLASER': {'sLASER TE1'},
                'MEGA-sLASER': {'sLASER TE1'}, 'STEAM': {'TM'}}.get(layout(d.kind), set())
        if need and not (need & keys):
            p = {**b.optional_params, **b.mandatory_params, **route.sheet}
            fixed = (ss.FIXED_TIMING.get((b.name, p.get('Sequence')))
                     or ss.FIXED_TIMING.get((b.name, None)))
            if _edited(d.kind) or not _symmetric(d) or d.kind == 'STEAM':
                status = 'approx' if status == 'ok' else status
                notes.append(fixed or f"{eng} places the pulses its own way.")
    return Plan(status, notes, route)


def apply(br, d: Design, path: str, category: str | None = None) -> Plan:
    """Fill the sheet of ``category`` (default: the active engine) from the design saved at
    ``path`` (the file the pulses are read from) and switch to it. Returns the plan; nothing
    changes when the engine cannot run it."""
    cat = category or br.backend.category
    pl = plan(d, cat, br, path)
    if pl.status == 'no':
        return pl
    if cat == 'MRSCloud':
        if br.backend.name != 'MRSCloud':
            br.set_backend('MRSCloud')
        b = br.backend
        b.mandatory_params.update(pl.route.sheet)            # Sequence, Localization
        b.optional_params['Sequence File'] = path
        vals = {'TE': d.te}
        if d.kind == 'MEGA-PRESS':
            vals.update({'Edit On': d.edit[0], 'Edit Off': d.edit[1], 'Edit Tp': d.pulses['edit']['dur']})
    elif cat == 'FSL-MRS':
        if br.backend.name != 'FSL-MRS':
            br.set_backend('FSL-MRS')
        b = br.backend
        b.set_mode('Custom')
        b.optional_params['Custom Sequence'] = path
        vals = {'TE': d.te}
        if d.kind in b.dropdown.get('Sequence', []):
            vals['Sequence'] = d.kind
    else:
        want = {r: ss.IDEAL if p['source'] == 'ideal' else ss.SHAPED for r, p in d.pulses.items()}
        ss.apply(br, cat, d.kind, want)
        b = br.backend
        vals = _sheet_values(b, d, path)
    p = {**b.optional_params, **b.mandatory_params}
    for k, v in (('Samples', d.samples), ('Bandwidth', d.bandwidth)):
        if v and k in p and p[k] in (None, '', 'missing input'):   # the scan's, when it has none
            vals[k] = v
    for k, v in vals.items():
        (b.mandatory_params if k in b.mandatory_params or k not in b.optional_params
         else b.optional_params)[k] = v
    br.from_file.setdefault(b.name, {}).update(vals)
    pl.values = vals
    return pl


def _sheet_values(b, d: Design, path: str) -> dict:
    """The sheet values of backend ``b`` that carry the design."""
    from basisremy.core import sequence_view as sv
    from basisremy.core.pulse_library import bandwidth_hz
    shown = b.get_params_for_mode()
    t = d.timing
    vals = {'TE': d.te}
    lay = layout(d.kind)
    if lay in ('PRESS', 'MEGA-PRESS'):
        vals.update({'Tau 1': t['TE1'], 'Tau 2': t['TE2'], 'TE2': t['TE2']})
    elif lay in ('sLASER', 'MEGA-sLASER'):
        vals.update({f'sLASER TE{i}': t[f'TE{i}'] for i in (1, 2, 3)})
    elif d.kind == 'STEAM':
        vals['TM'] = t['TM']
    for role, p in d.pulses.items():
        key = sv.pulse_key(b, role)
        if p['source'] == 'ideal' or not key:
            continue
        vals[key] = f'{path}#{role}'
        vals['Edit Tp' if role == 'edit' else 'RefTp'] = p['dur']
        if role != 'edit':
            vals['Flip Angle'] = _FLIP[role]
        if role == 'edit':
            try:
                pulse = sv.read_pulse(f'{path}#edit', 'inv')
                vals['Edit Bandwidth (Hz)'] = round(bandwidth_hz(pulse, p['dur']), 1)
            except Exception:                              # noqa: BLE001 - keeps the sheet's
                pass
    if _edited(d.kind) and not d.scheme:
        vals['Edit On'], vals['Edit Off'] = d.edit
    if any(p['source'] != 'ideal' for r, p in d.pulses.items() if r != 'edit'):
        # the engines' two simulated slabs: the first and the second refocusing (STEAM: the
        # second and third 90) pulse, i.e. the voxel's y and z
        for ax, size in (('X', d.voxel[1]), ('Y', d.voxel[2])):
            vals[f'thk{ax}'] = size
            vals[f'fov{ax}'] = round(2 * size, 3)          # grid covers the transition bands
    return {k: v for k, v in vals.items() if k in shown}


#**************************************************************************************************#
#                                      FSL-MRS translation                                         #
#**************************************************************************************************#
_FSL_POINTS = 500


def fsl_sequences(path: str, params: dict) -> dict:
    """{label: FSL-MRS sequence description} of a .seq, one per acquisition (MEGA: 'ON', 'OFF' by
    the editing frequency; otherwise a single None). B0, the readout and the linewidth come from
    ``params``; RF amplitudes in Hz, gradients in mT/m, rephasing of the excitation slab."""
    import pypulseq as pp
    from basisremy.core import pulseq
    bfield = float(params['Bfield'])
    hz_per_ppm = bfield * 42.577
    raw = pp.Sequence()
    raw.read(path)
    labels = [part.partition('=')[0].strip() for part in
              str(raw.definitions.get('EditScheme', '')).split(';') if '=' in part]
    out, acq = {}, 0
    while True:
        try:
            info = pulseq.read_seq(path, acq)
        except ValueError:
            break
        rf = [r for r in info.rf if r.role in ('exc', 'ref', 'edit')]
        kind = kind_of([r.role for r in rf])
        te = pulseq.echo_ms(info)
        blocks, cfilter, reph = [], [], []
        steam = kind == 'STEAM'
        # coherence orders counted back from the readout (-1): an odd number of refocusing
        # pulses (Spin Echo, MEGA-SPECIAL) starts the excitation at +1
        left = sum(r.role == 'ref' for r in rf)
        for i, r in enumerate(rf):
            # at most ~_FSL_POINTS samples (block means): denmatsim's cost grows with every sample
            step = max(1, math.ceil(len(r.signal) / _FSL_POINTS))
            n = len(r.signal) // step
            signal = r.signal[:n * step].reshape(n, step).mean(axis=1)
            t = (np.arange(n) + 0.5) / n * r.dur_ms / 1e3
            off = r.freq_offset_hz + r.freq_ppm * hz_per_ppm
            wave = signal * np.exp(1j * (r.phase_rad + 2 * np.pi * off * (t - r.dur_ms / 2e3)))
            if r.grad_hz_m is not None:
                grad = np.asarray(r.grad_hz_m)[:n * step].reshape(n, step).mean(axis=1)
                g = [[0.0] * n for _ in range(3)]
                g['xyz'.index(r.slab_axis)] = (grad / GAMMA_HZ_T * 1e3).tolist()
            else:
                g = [0.0, 0.0, 0.0]
                if r.slab_axis:
                    g['xyz'.index(r.slab_axis)] = r.slab_grad_hz_m / GAMMA_HZ_T * 1e3
            blocks.append({'time': r.dur_ms / 1e3, 'frequencyOffset': 0, 'phaseOffset': 0,
                           'amp': np.abs(wave).tolist(), 'phase': np.angle(wave).tolist(),
                           'grad': g})
            if r.role == 'edit':
                cfilter.append(None)
            elif r.role == 'exc':
                cfilter.append([1, 0, -1][i] if steam else (1 if left % 2 else -1))
            else:
                cfilter.append(-1 if left % 2 else 1)
                left -= 1
            reph.append([0.0, 0.0, 0.0])
        # excitation slabs: undo the gradient area after the centre (STEAM's second pulse: the
        # area before its centre, applied after the first pulse, as in FSL-MRS's example)
        for i, r in enumerate(rf):
            if r.role != 'exc' or not r.slab_axis or r.grad_hz_m is not None:
                continue
            area = r.slab_grad_hz_m / GAMMA_HZ_T * 1e3 * r.dur_ms / 2e3
            ax = 'xyz'.index(r.slab_axis)
            if steam and i == 1:
                reph[0][ax] -= area
            else:
                reph[i][ax] -= area
        delays = []
        for a, b in zip(rf, rf[1:]):
            delays.append(((b.centre_ms - b.dur_ms / 2) - (a.centre_ms + a.dur_ms / 2)) / 1e3)
        t0 = rf[0].centre_ms
        echo_at = t0 + te + (rf[2].centre_ms - rf[1].centre_ms if steam else 0.0)
        delays.append((echo_at - (rf[-1].centre_ms + rf[-1].dur_ms / 2)) / 1e3)
        delays = [0.0 if -1e-8 < x < 0 else x for x in delays]     # rounding, not an overlap
        if min(delays) < 0:
            raise ValueError(f"{os.path.basename(path)}: pulses overlap; FSL-MRS cannot run it.")
        # grid per axis: the slab selected along it (the whole slab plus its transition bands)
        half_mm = {}
        for r in rf:
            if r.slab_axis:
                half_mm[r.slab_axis] = max(half_mm.get(r.slab_axis, 0.0), (pulseq.slab_cm(r) or 0) * 10.0)
        res = int(_num(params.get('Spatial Points')) or 10)
        seq = {
            'sequenceName': os.path.splitext(os.path.basename(path))[0],
            'description': f'{kind} from {os.path.basename(path)} (BasisREMY .seq translation)',
            'B0': bfield, 'centralShift': WATER_PPM,
            'Rx_Points': int(float(params['Samples'])), 'Rx_SW': float(params['Bandwidth']),
            'Rx_LW': float(params.get('Linewidth') or 1.0), 'Rx_Phase': -1.5708,
            'x': [-15, 15], 'y': [-15, 15], 'z': [-15, 15], 'resolution': [8, 8, 8],
            'RFUnits': 'Hz', 'GradUnits': 'mT', 'spaceUnits': 'mm',
            'RF': blocks, 'delays': delays, 'rephaseAreas': reph, 'CoherenceFilter': cfilter,
        }
        if half_mm:
            seq['resolution'] = [res if a in half_mm else 1 for a in _AXES]
            for a, h in half_mm.items():
                seq[a] = [-h, h]
        edit = next((r for r in rf if r.role == 'edit'), None)
        label = None
        if acq < len(labels):                     # the design's own sub-experiment names
            label = labels[acq]
        elif edit is not None:
            ppm = WATER_PPM + edit.freq_ppm + edit.freq_offset_hz / hz_per_ppm
            label = 'ON' if ppm < WATER_PPM else 'OFF'
            if label in out:
                label = f'{label}{acq}'
        out[label] = seq
        acq += 1
    if not out:
        raise ValueError(f"{os.path.basename(path)}: no acquisition in the sequence file.")
    return out


#**************************************************************************************************#
#                                        saved designs                                             #
#**************************************************************************************************#
def designs_dir() -> str:
    return os.environ.get('BASISREMY_SEQUENCES_DIR') or os.path.join(
        os.path.expanduser('~'), 'BasisREMY', 'sequences')


def saved_designs() -> list[str]:
    folder = designs_dir()
    if not os.path.isdir(folder):
        return []
    files = [os.path.join(folder, f) for f in os.listdir(folder) if f.lower().endswith('.seq')]
    return sorted(files, key=os.path.getmtime, reverse=True)


def default_name(d: Design) -> str:
    shaped = [os.path.basename(str(p['source']).partition('#')[0]).replace('standard:', '')
              for r, p in d.pulses.items() if p['source'] != 'ideal' and r != 'edit']
    tag = '_'.join(dict.fromkeys(shaped)) if shaped else 'ideal'
    tag = os.path.splitext(tag)[0]
    return f"{d.kind.replace(' ', '')}_TE{d.te:g}_{tag}".replace('/', '-')

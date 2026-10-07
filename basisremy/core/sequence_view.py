####################################################################################################
#                                        sequence_view.py                                          #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 07/10/26                                                                                #
#                                                                                                  #
# Purpose: What the Sequence panel shows, independent of the GUI:                                  #
#            - read_pulse(): one reader for every pulse source BasisREMY takes (standard pulse,    #
#              .pta / .RF / .txt / Bruker / FID-A .mat waveform, Pulseq .seq, FSL-MRS or WIN       #
#              sequence / basis JSON) into a pulse_library.Pulse;                                  #
#            - timeline(): the RF events and the echo of a sequence, exact from a .seq or FSL-MRS  #
#              description, else from the sheet (TE, echo split, TM, pulse durations);             #
#            - RUNS_ON: which engine simulates which sequence with ideal pulses, with a shaped     #
#              pulse, or from a whole sequence file (read from the backends, 7 Oct 2026).          #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import json
import os

import numpy as np

from basisremy.core.pulse_library import Pulse, is_standard, make_standard, standard_name

# sequence -> {'ideal' | 'shaped' | 'file': engines}. 'shaped' = a pulse waveform per role
# (MRSCloud: its vendor / universal set, the open GOIA stand-in); 'file' = the whole sequence
# from a file (Pulseq .seq -> FID-A shaped, FSL-MRS / WIN .json -> FSL-MRS Custom).
RUNS_ON = {
    'PRESS':        {'ideal': ['FID-A', 'FSL-MRS', 'Vespa', 'spant', 'Spinach'],
                     'shaped': ['FID-A', 'Spinach', 'Vespa', 'spant', 'MRSCloud'],
                     'file': ['FID-A (.seq)', 'FSL-MRS (.json)']},
    'sLASER':       {'ideal': ['FSL-MRS', 'spant'],
                     'shaped': ['FID-A', 'Spinach', 'MRSCloud', 'jbss'],
                     'file': ['FID-A (.seq)', 'FSL-MRS (.json)']},
    'LASER':        {'ideal': ['FID-A', 'FSL-MRS', 'Spinach'], 'shaped': [],
                     'file': ['FSL-MRS (.json)']},
    'STEAM':        {'ideal': ['FID-A', 'FSL-MRS', 'Vespa', 'spant', 'Spinach'], 'shaped': ['FID-A'],
                     'file': ['FID-A (.seq)', 'FSL-MRS (.json)']},
    'Spin Echo':    {'ideal': ['FID-A', 'Vespa', 'spant', 'Spinach'], 'shaped': ['FID-A'],
                     'file': ['FSL-MRS (.json)']},
    'MEGA-PRESS':   {'ideal': ['FID-A', 'FSL-MRS', 'spant'], 'shaped': ['FID-A', 'MRSCloud'],
                     'file': ['FID-A (.seq)', 'FSL-MRS (.json)']},
    'MEGA-sLASER':  {'ideal': ['FSL-MRS'], 'shaped': ['MRSCloud'], 'file': ['FSL-MRS (.json)']},
    'MEGA-SPECIAL': {'ideal': [], 'shaped': ['FID-A'], 'file': ['FSL-MRS (.json)']},
    'HERMES':       {'ideal': ['FSL-MRS'], 'shaped': ['MRSCloud'], 'file': ['FSL-MRS (.json)']},
    'HERCULES':     {'ideal': ['FSL-MRS'], 'shaped': ['MRSCloud'], 'file': ['FSL-MRS (.json)']},
}

_EDITED = {'MEGA-PRESS', 'MEGA-sLASER', 'MEGA-SPECIAL', 'HERMES', 'HERCULES'}


def roles(kind: str | None) -> list[str]:
    """The pulse roles of a sequence (STEAM: three excitations, no refocusing)."""
    if kind is None:
        return []
    if kind == 'STEAM':
        return ['exc']
    return ['exc', 'ref'] + (['edit'] if kind in _EDITED else [])


# FID-A backends carry their sequence in the backend, not in a 'Sequence' value
_FIDA_KIND = {'FidaPressShaped': 'PRESS', 'FidaSemiLaserShaped': 'sLASER', 'FidaSteamShaped': 'STEAM',
              'FidaSpinEchoShaped': 'Spin Echo', 'FidaMegaPressShaped': 'MEGA-PRESS',
              'FidaMegaSpecialShaped': 'MEGA-SPECIAL', 'FidaLaser': 'LASER',
              'FidaMegaPressIdeal': 'MEGA-PRESS', 'FidaSpinEchoXN': 'Spin Echo',
              'SpinachPressShaped': 'PRESS', 'SpinachSemiLaserShaped': 'sLASER', 'CustomSLaser': 'sLASER'}


#**************************************************************************************************#
#                                     ideal <-> waveform pulse                                     #
#**************************************************************************************************#
# backends whose 'Path to Pulse' is the excitation (STEAM's 90s, the FID), not the refocusing
_EXC_PULSE = {'FidaSteamShaped', 'FidaOnePulse'}
# the same sequence with the waveform pulse made ideal (and back): shaped backend -> ideal one
_TO_IDEAL = {'FidaPressShaped': ('FidaIdeal', 'PRESS'), 'FidaSteamShaped': ('FidaIdeal', 'STEAM'),
             'FidaSpinEchoShaped': ('FidaIdeal', 'Spin Echo'), 'SpinachPressShaped': ('Spinach', 'PRESS')}
_TO_SHAPED = {v: k for k, v in _TO_IDEAL.items()}
# FID-A MEGA-PRESS: which mode makes one role ideal / shaped
_MEGA = 'FidaMegaPressShaped'
_FULL, _EDIT_ONLY, _REFOC_ONLY = ('Full shaped (refoc + edit)', 'Edit-only shaped (ideal refoc)',
                                  'Refoc-only shaped (ideal edit)')


def pulse_role(backend) -> str:
    """The role the backend's 'Path to Pulse' plays."""
    return 'exc' if backend.name in _EXC_PULSE else 'ref'


def pulse_key(backend, role: str) -> str | None:
    """The sheet key holding the waveform of ``role``, None when that pulse is ideal here."""
    shown = backend.get_params_for_mode()
    if role == 'edit':
        return 'Edit Pulse Path' if 'Edit Pulse Path' in shown else None
    return 'Path to Pulse' if 'Path to Pulse' in shown and pulse_role(backend) == role else None


def switch_target(backend, role: str, to_ideal: bool) -> dict | None:
    """How to make ``role`` ideal (to_ideal) or a waveform pulse in the same sequence:
    {'backend': name, 'Sequence': value, 'mode': mode} (only the keys that change), or None when
    no engine of this software has that variant."""
    p = {**backend.optional_params, **backend.mandatory_params}
    seq, name = p.get('Sequence'), backend.name
    if name in ('Vespa', 'Spant') and role == 'ref':
        if to_ideal and seq == 'PRESS shaped':
            return {'Sequence': 'PRESS'}
        if not to_ideal and seq == 'PRESS':
            return {'Sequence': 'PRESS shaped'}
        return None
    if name == _MEGA:
        mode = backend.current_mode
        table = {('ref', True): {_FULL: _EDIT_ONLY, _REFOC_ONLY: None},
                 ('edit', True): {_FULL: _REFOC_ONLY, _EDIT_ONLY: None},
                 ('ref', False): {_EDIT_ONLY: _FULL}, ('edit', False): {_REFOC_ONLY: _FULL}}
        if mode not in table[(role, to_ideal)]:
            return None
        target = table[(role, to_ideal)][mode]
        return {'mode': target} if target else {'backend': 'FidaMegaPressIdeal'}
    if name == 'FidaMegaPressIdeal' and not to_ideal and role in ('ref', 'edit'):
        return {'backend': _MEGA, 'mode': _REFOC_ONLY if role == 'ref' else _EDIT_ONLY}
    if to_ideal and name in _TO_IDEAL and role == pulse_role(backend):
        target, sequence = _TO_IDEAL[name]
        return {'backend': target, 'Sequence': sequence}
    if not to_ideal and (name, seq) in _TO_SHAPED:
        target = _TO_SHAPED[(name, seq)]
        return {'backend': target} if role == ('exc' if target in _EXC_PULSE else 'ref') else None
    return None


#**************************************************************************************************#
#                                             pulses                                               #
#**************************************************************************************************#
_KIND = {'exc': 'exc', 'ref': 'ref', 'edit': 'inv'}   # role -> how the engines scale B1


def read_pulse(spec: str, kind: str = 'ref', index: int | None = None) -> Pulse:
    """A pulse from any source BasisREMY reads. ``kind`` picks the role in a whole-sequence file
    (first excitation / refocusing / editing pulse) unless ``index`` names the RF event."""
    if is_standard(spec):
        return make_standard(standard_name(spec))
    spec, _, sel = str(spec).partition('#')       # 'file.json#edit' / 'file.seq#3' (1-based)
    if sel:
        if sel.isdigit():
            index = int(sel) - 1
        elif sel in _KIND:
            kind = sel
        else:
            raise ValueError(f"Unknown pulse selector '#{sel}' (use #exc, #ref, #edit or #<n>)")
    ext = os.path.splitext(spec)[1].lower()
    name = os.path.basename(spec)
    if ext in ('.pta', '.rf', '.txt', '.exc', '.rfc', '.inv', '.mat'):
        from basisremy.core.rf_pulses import read_waveform
        return Pulse(name, read_waveform(spec), _KIND[kind], float('nan'))
    if ext == '.seq':
        from basisremy.core import pulseq
        rf = pulseq.read_seq(spec).rf
        event = rf[index] if index is not None else next((r for r in rf if r.role == kind), None)
        if event is None:
            raise ValueError(f"{name}: no '{kind}' pulse in the sequence")
        return pulseq.waveform_pulse(event, _KIND.get(event.role, 'ref'))
    if ext == '.json':
        seq = _fsl_sequence(spec)
        roles = _fsl_roles(seq)
        i = index if index is not None else next((j for j, r in enumerate(roles) if r == kind), None)
        if i is None:
            raise ValueError(f"{name}: no '{kind}' pulse in the sequence")
        return block_pulse(seq['RF'][i], roles[i], f'{name} RF {i + 1}')
    raise ValueError(f"Unrecognised pulse source '{name}' (standard:<name>, .pta, .RF, .txt, Bruker "
                     ".exc/.rfc/.inv, FID-A .mat, Pulseq .seq, FSL-MRS / WIN .json)")


def block_pulse(block: dict, role: str, name: str = 'FSL-MRS RF') -> Pulse:
    """One FSL-MRS RF block (amp Hz, phase rad, time s) as a Pulse."""
    amp = np.asarray(block['amp'], dtype=float)
    phase = np.degrees(np.asarray(block['phase'], dtype=float)) % 360.0
    wave = np.column_stack([phase, amp / amp.max(), np.ones(len(amp))])
    return Pulse(name, wave, _KIND[role], float(block['time']) * 1e3)


def _fsl_sequence(path: str) -> dict:
    """The FSL-MRS sequence description in a sequence JSON or a basis JSON (its 'seq')."""
    with open(path) as fh:
        d = json.load(fh)
    seq = d.get('seq', d)
    if 'RF' not in seq or 'delays' not in seq:
        raise ValueError(f"{os.path.basename(path)} is not an FSL-MRS sequence description")
    return seq


def _fsl_roles(seq: dict) -> list[str]:
    """FSL-MRS RF blocks carry no role; the coherence filter tells: None after an editing pulse
    (selective, so no filter), a 0 only in STEAM (all three pulses excite), and a pulse off the
    excitation's frequency by more than 100 Hz edits too. The first pulse excites, the rest refocus."""
    cf = list(seq.get('CoherenceFilter') or [])
    cf += [1] * (len(seq['RF']) - len(cf))
    f0 = float(seq['RF'][0].get('frequencyOffset', 0))
    steam = 0 in cf
    roles = []
    for i, b in enumerate(seq['RF']):
        if cf[i] is None or abs(float(b.get('frequencyOffset', 0)) - f0) > 100:
            roles.append('edit')
        else:
            roles.append('exc' if i == 0 or steam else 'ref')
    return roles


#**************************************************************************************************#
#                                            timeline                                              #
#**************************************************************************************************#
def sequence_kind(backend) -> str | None:
    """The sequence of the current sheet in RUNS_ON's names."""
    if backend.name in _FIDA_KIND:
        return _FIDA_KIND[backend.name]
    p = {**backend.optional_params, **backend.mandatory_params}
    seq = str(p.get('Sequence') or '').replace(' shaped', '')
    if backend.name == 'MRSCloud':
        loc = 'sLASER' if 'laser' in str(p.get('Localization', '')).lower() else 'PRESS'
        return {'UnEdited': loc, 'MEGA': f'MEGA-{loc}'}.get(seq, seq or None)
    return seq if seq in RUNS_ON else None


def _num(v, default=None):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def timeline(backend, seq_file: str | None = None) -> dict:
    """{'events': [{'role', 'centre_ms', 'dur_ms', 'pulse'}], 'echo_ms', 'exact', 'kind'}; 'pulse' is
    a sheet value, (file, RF index) or an FSL-MRS RF block, None when ideal.

    exact: from a Pulseq or FSL-MRS file (the FSL-MRS backend's Custom Sequence, or seq_file);
    otherwise placed from the sheet: TE, the echo split (Tau 1 / Tau 2, sLASER TE1-3), TM and
    the pulse durations it has (RefTp, Edit Tp; an ideal pulse has none)."""
    p = {**backend.optional_params, **backend.mandatory_params}
    custom = p.get('Custom Sequence') if getattr(backend, 'current_mode', None) == 'Custom' else None
    path = seq_file or custom
    if path and os.path.exists(path):
        return _file_timeline(path)
    kind = sequence_kind(backend)
    te = _num(p.get('TE'))
    if kind is None or te is None:
        return {'events': [], 'echo_ms': None, 'exact': False, 'kind': kind}
    if backend.name == 'FSL-MRS':
        # the ideal sequence FSL-MRS will run (its own MEGA timing), not a schematic one
        try:
            edit = (_num(p.get('Edit On'), 1.9),) if kind in backend._edited_sequences else None
            seq = backend._generate_sequence_json(backend._coerce_params(dict(p)), edit)
            return {**_seq_timeline(seq, None), 'kind': kind}
        except Exception:                                   # noqa: BLE001 - incomplete sheet
            pass
    shown = backend.get_params_for_mode() if hasattr(backend, 'get_params_for_mode') else p
    ref_tp = _num(p.get('RefTp'), 0.0) if shown.get('Path to Pulse') else 0.0
    edit_tp = _num(p.get('Edit Tp'), 0.0) if shown.get('Edit Pulse Path') or 'Edit Tp' in shown else 0.0
    ev = [('exc', 0.0, 0.0)]
    taus = getattr(backend, '_TE68_TAUS', None)     # FID-A MEGA kinds: their own timing, scaled to TE
    if taus and kind in ('MEGA-PRESS', 'MEGA-SPECIAL'):
        t = 0.0
        order = ['ref', 'edit', 'ref', 'edit'] if kind == 'MEGA-PRESS' else ['edit', 'ref', 'edit']
        for role, tau in zip(order, taus):
            t += tau * te / 68.0
            ev.append((role, t, ref_tp if role == 'ref' else edit_tp))
    elif kind in ('PRESS', 'MEGA-PRESS', 'MEGA-SPECIAL', 'HERMES', 'HERCULES'):
        t1 = _num(p.get('Tau 1'), te / 2)
        t2 = _num(p.get('Tau 2'), te - t1)
        ev += [('ref', t1 / 2, ref_tp), ('ref', t1 + t2 / 2, ref_tp)]
        if kind != 'PRESS':   # editing pulses in the two halves of the second echo
            ev += [('edit', t1 + t2 / 4, edit_tp), ('edit', t1 + 3 * t2 / 4, edit_tp)]
    elif kind == 'Spin Echo':
        ev += [('ref', te / 2, ref_tp)]
    elif kind in ('sLASER', 'MEGA-sLASER'):
        split = [_num(p.get(f'sLASER TE{i}')) for i in (1, 2, 3)]
        te1, te2, te3 = split if None not in split else (te / 4, te / 2, te / 4)
        t = 0.0
        for gap in (te1 / 2, te1 / 2 + te2 / 4, te2 / 2, te2 / 4 + te3 / 2):
            t += gap
            ev.append(('ref', t, ref_tp))
        if kind == 'MEGA-sLASER':
            ev += [('edit', te / 4, edit_tp), ('edit', 3 * te / 4, edit_tp)]
    elif kind == 'LASER':
        ev += [('ref', te * (2 * i + 1) / 12, ref_tp) for i in range(6)]
    elif kind == 'STEAM':
        tm = _num(p.get('TM'), 0.0)
        exc_tp = _num(p.get('RefTp'), 0.0) if pulse_key(backend, 'exc') else 0.0
        ev = [('exc', 0.0, exc_tp), ('exc', te / 2, exc_tp), ('exc', te / 2 + tm, exc_tp)]
        te = te + tm
    pulse = {r: shown.get(pulse_key(backend, r)) if pulse_key(backend, r) else None
             for r in ('exc', 'ref', 'edit')}
    events = [{'role': r, 'centre_ms': c, 'dur_ms': d, 'pulse': pulse[r]}
              for r, c, d in sorted(ev, key=lambda e: e[1])]
    return {'events': events, 'echo_ms': te, 'exact': False, 'kind': kind}


def _file_timeline(path: str) -> dict:
    if path.lower().endswith('.seq'):
        from basisremy.core import pulseq
        info = pulseq.read_seq(path)
        rf = [r for r in info.rf if r.role in ('exc', 'ref', 'edit')]
        t0 = next(r.centre_ms for r in rf if r.role == 'exc')
        events = [{'role': r.role, 'centre_ms': r.centre_ms - t0, 'dur_ms': r.dur_ms,
                   'pulse': (path, i)} for i, r in enumerate(info.rf) if r in rf]
        return {'events': events, 'echo_ms': pulseq.echo_ms(info), 'exact': True,
                'kind': pulseq.sequence_type(info)}
    return _seq_timeline(_fsl_sequence(path), path)


def _seq_timeline(seq: dict, path: str | None) -> dict:
    """Timeline of an FSL-MRS sequence description: RF blocks and the delays between them."""
    roles = _fsl_roles(seq)
    events, t = [], 0.0
    for i, block in enumerate(seq['RF']):
        dur = float(block['time']) * 1e3
        events.append({'role': roles[i], 'centre_ms': t + dur / 2, 'dur_ms': dur,
                       'pulse': (path, i) if path else block})
        t += dur + (float(seq['delays'][i]) * 1e3 if i < len(seq['delays']) else 0.0)
    t0 = events[0]['centre_ms']
    for e in events:
        e['centre_ms'] -= t0
    # FSL-MRS acquires right after the last delay; the echo sits there
    return {'events': events, 'echo_ms': t - t0, 'exact': True, 'kind': None}

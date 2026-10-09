####################################################################################################
#                                       sequence_setup.py                                          #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 07/10/26                                                                                #
#                                                                                                  #
# Purpose: One description of what to simulate, the same for every engine: an engine, a sequence   #
#          and a pulse per role (excitation / refocusing / editing: ideal, a waveform, or the      #
#          engine's own). ROUTES maps that onto the backends as they are (backend, mode, sheet     #
#          values), so the physics and every backend's own defaults stay untouched. Also: why an   #
#          engine cannot run a sequence or a pulse choice, the timings an engine fixes itself,     #
#          and the recommended values for what a data header never holds (echo split, sLASER       #
#          spacing, mixing time, editing pulse), each with its source.                             #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

from dataclasses import dataclass, field

IDEAL, SHAPED, OWN = 'ideal', 'shaped', 'own'     # OWN: the engine's own pulse, not selectable
ROLES = ('exc', 'ref', 'edit')
ROLE_NAME = {'exc': 'Excitation', 'ref': 'Refocusing', 'edit': 'Editing'}

# display order of the sequence picker
SEQUENCES = ['PRESS', 'sLASER', 'STEAM', 'Spin Echo', 'LASER', 'MEGA-PRESS', 'MEGA-sLASER',
             'MEGA-SPECIAL', 'HERMES', 'HERCULES', 'HERMES (sLASER)', 'HERCULES (sLASER)',
             'Pulse-acquire', 'Spin Echo train']

# BasisREMY category -> engine name shown in the picker
ENGINE_LABEL = {'MRSCloud': 'MRSCloud', 'FID-A': 'FID-A', 'FSL-MRS': 'FSL-MRS', 'Vespa': 'Vespa',
                'Spant': 'spant', 'Spinach': 'Spinach', 'Custom': 'Custom'}


@dataclass(frozen=True)
class Route:
    """One way an engine runs a sequence: the backend, its mode (the first of ``modes`` is the
    default; the others are engine options that keep the same pulses) and the sheet values that
    select the sequence inside the backend."""
    backend: str
    pulses: dict
    modes: tuple = ()
    sheet: dict = field(default_factory=dict)


def _r(backend, pulses, modes=(), **sheet):
    return Route(backend, pulses, tuple(modes), sheet)


_I2 = {'exc': IDEAL, 'ref': IDEAL}
_S2 = {'exc': IDEAL, 'ref': SHAPED}
_EDIT_OWN = {'exc': IDEAL, 'ref': IDEAL, 'edit': OWN}
_MRSC = {'exc': OWN, 'ref': OWN}
_MRSC_EDIT = {'exc': OWN, 'ref': OWN, 'edit': OWN}
_MRSC_MODES = ('Universal', 'Non-Universal')

ROUTES: dict[str, dict[str, list[Route]]] = {
    'MRSCloud': {
        'PRESS':  [_r('MRSCloud', _MRSC, _MRSC_MODES, Sequence='UnEdited', Localization='PRESS')],
        'sLASER': [_r('MRSCloud', _MRSC, _MRSC_MODES, Sequence='UnEdited', Localization='sLASER')],
        'STEAM':  [_r('MRSCloud', {'exc': OWN}, _MRSC_MODES, Sequence='UnEdited', Localization='STEAM_7T')],
        'MEGA-PRESS':  [_r('MRSCloud', _MRSC_EDIT, _MRSC_MODES, Sequence='MEGA', Localization='PRESS')],
        'MEGA-sLASER': [_r('MRSCloud', _MRSC_EDIT, _MRSC_MODES, Sequence='MEGA', Localization='sLASER')],
        'HERMES':      [_r('MRSCloud', _MRSC_EDIT, _MRSC_MODES, Sequence='HERMES', Localization='PRESS')],
        'HERCULES':    [_r('MRSCloud', _MRSC_EDIT, _MRSC_MODES, Sequence='HERCULES', Localization='PRESS')],
        'HERMES (sLASER)':   [_r('MRSCloud', _MRSC_EDIT, _MRSC_MODES, Sequence='HERMES',
                                 Localization='sLASER')],
        'HERCULES (sLASER)': [_r('MRSCloud', _MRSC_EDIT, _MRSC_MODES, Sequence='HERCULES',
                                 Localization='sLASER')],
    },
    'FID-A': {
        'PRESS':     [_r('FidaIdeal', _I2, Sequence='PRESS'), _r('FidaPressShaped', _S2)],
        'sLASER':    [_r('FidaSemiLaserShaped', _S2, ('Standard', 'Phase cycled'))],
        'STEAM':     [_r('FidaIdeal', {'exc': IDEAL}, Sequence='STEAM'),
                      _r('FidaSteamShaped', {'exc': SHAPED})],
        'Spin Echo': [_r('FidaIdeal', _I2, Sequence='Spin Echo'), _r('FidaSpinEchoShaped', _S2)],
        'LASER':     [_r('FidaIdeal', _I2, Sequence='LASER')],
        'MEGA-PRESS': [_r('FidaMegaPressShaped', {'exc': IDEAL, 'ref': IDEAL, 'edit': SHAPED},
                          ('Edit-only shaped (ideal refoc)',)),
                       _r('FidaMegaPressShaped', {'exc': IDEAL, 'ref': SHAPED, 'edit': SHAPED},
                          ('Full shaped (refoc + edit)',))],
        'MEGA-SPECIAL': [_r('FidaMegaSpecialShaped', {'exc': IDEAL, 'ref': SHAPED, 'edit': SHAPED})],
        'Pulse-acquire': [_r('FidaOnePulse', {'exc': IDEAL}, ('Ideal', 'Delay', 'Arbitrary phase')),
                          _r('FidaOnePulse', {'exc': SHAPED}, ('Shaped',))],
        'Spin Echo train': [_r('FidaSpinEchoXN', _I2)],
    },
    'FSL-MRS': {
        'PRESS':  [_r('FSL-MRS', _I2, ('Simple',), Sequence='PRESS')],
        'sLASER': [_r('FSL-MRS', _I2, ('Simple',), Sequence='sLASER')],
        'STEAM':  [_r('FSL-MRS', {'exc': IDEAL}, ('Simple',), Sequence='STEAM')],
        'LASER':  [_r('FSL-MRS', _I2, ('Simple',), Sequence='LASER')],
        'MEGA-PRESS':  [_r('FSL-MRS', _EDIT_OWN, ('Simple',), Sequence='MEGA-PRESS')],
        'MEGA-sLASER': [_r('FSL-MRS', _EDIT_OWN, ('Simple',), Sequence='MEGA-sLASER')],
        'HERMES':      [_r('FSL-MRS', _EDIT_OWN, ('Simple',), Sequence='HERMES')],
        'HERCULES':    [_r('FSL-MRS', _EDIT_OWN, ('Simple',), Sequence='HERCULES')],
    },
    'Vespa': {
        'PRESS':     [_r('Vespa', _I2, Sequence='PRESS'), _r('Vespa', _S2, Sequence='PRESS shaped')],
        'STEAM':     [_r('Vespa', {'exc': IDEAL}, Sequence='STEAM')],
        'Spin Echo': [_r('Vespa', _I2, Sequence='Spin Echo')],
    },
    'Spant': {
        'PRESS':     [_r('Spant', _I2, Sequence='PRESS'), _r('Spant', _S2, Sequence='PRESS shaped')],
        'sLASER':    [_r('Spant', _I2, Sequence='sLASER')],
        'STEAM':     [_r('Spant', {'exc': IDEAL}, Sequence='STEAM')],
        'Spin Echo': [_r('Spant', _I2, Sequence='Spin Echo')],
        'MEGA-PRESS':    [_r('Spant', _EDIT_OWN, Sequence='MEGA-PRESS')],
        'Pulse-acquire': [_r('Spant', {'exc': IDEAL}, Sequence='Pulse-acquire')],
    },
    'Spinach': {
        'PRESS':     [_r('Spinach', _I2, Sequence='PRESS'), _r('SpinachPressShaped', _S2)],
        'sLASER':    [_r('SpinachSemiLaserShaped', _S2)],
        'STEAM':     [_r('Spinach', {'exc': IDEAL}, Sequence='STEAM')],
        'Spin Echo': [_r('Spinach', _I2, Sequence='Spin Echo')],
        'LASER':     [_r('Spinach', _I2, Sequence='LASER')],
    },
    'Custom': {
        'sLASER': [_r('CustomSLaser', _S2)],
    },
}

# What each engine runs in BasisREMY, and why the rest is missing: "not wired yet" where the
# engine itself could, "not in <engine>" where it cannot.
ENGINE_SCOPE = {
    'MRSCloud': "MRSCloud runs only its own sequences (PRESS, sLASER, 7 T STEAM, MEGA, HERMES, "
                "HERCULES) with its own pulse sets; no other sequence or pulse can be given to it.",
    'FID-A': "FID-A has no function for MEGA-sLASER, HERMES or HERCULES, and its semi-LASER and "
             "MEGA-SPECIAL always use a refocusing waveform; they could be built from FID-A's parts "
             "(not wired yet).",
    'FSL-MRS': "FSL-MRS runs its sequences with ideal pulses here; waveforms reach it through a whole "
               "sequence file. A waveform per pulse and Spin Echo are not wired yet (FSL-MRS itself "
               "could).",
    'Vespa': "Vespa runs PRESS (ideal, or with a refocusing waveform), STEAM and Spin Echo here; "
             "GAMMA could simulate the others (not wired yet).",
    'Spant': "spant has no LASER, MEGA-sLASER, HERMES or HERCULES function; a waveform is wired for "
             "PRESS only.",
    'Spinach': "Spinach runs ideal Spin Echo / PRESS / STEAM / LASER and shaped PRESS / sLASER here; "
               "the others are not wired yet (Spinach itself could).",
    'Custom': "Custom is a dedicated semi-LASER simulation with a refocusing waveform.",
}

# Timings an engine sets itself (no field for them in the sheet)
FIXED_TIMING = {
    ('FidaSemiLaserShaped', None): "Pulse spacing: FID-A's own semi-LASER timing for this TE.",
    ('FidaMegaPressShaped', None): "Pulse spacing: FID-A's MEGA-PRESS timing, scaled to TE.",
    ('FidaMegaSpecialShaped', None): "Pulse spacing: FID-A's MEGA-SPECIAL timing, scaled to TE.",
    ('FSL-MRS', 'PRESS'): "Echo split: symmetric (TE/2 each) in FSL-MRS.",
    ('FSL-MRS', 'MEGA-PRESS'): "Pulse spacing: FID-A's Siemens MEGA-PRESS timing, scaled to TE.",
    ('Vespa', 'PRESS'): "Echo split: symmetric (TE/2 each) in Vespa.",
    ('Vespa', 'PRESS shaped'): "Echo split: symmetric (TE/2 each) in Vespa.",
    ('Spant', 'MEGA-PRESS'): "Echo split: spant's MEGA-PRESS timing (15 / 53 ms at TE 68), scaled to TE.",
    ('MRSCloud', None): "Timing: the vendor's, from MRSCloud's parameter set.",
}


#**************************************************************************************************#
#                                        routes and state                                          #
#**************************************************************************************************#
def roles(sequence: str | None) -> list[str]:
    """The pulse roles of a sequence (STEAM and pulse-acquire only excite)."""
    if sequence is None:
        return []
    if sequence in ('STEAM', 'Pulse-acquire'):
        return ['exc']
    edited = sequence.startswith(('MEGA', 'HERMES', 'HERCULES'))
    return ['exc', 'ref'] + (['edit'] if edited else [])


def engines() -> list[str]:
    return list(ROUTES)


def sequences(category: str) -> list[str]:
    """The sequences ``category`` runs, in picker order."""
    return [s for s in SEQUENCES if s in ROUTES.get(category, {})]


def runs_on(sequence: str) -> list[str]:
    """Engines (display names) that run ``sequence``."""
    return [ENGINE_LABEL[c] for c in ROUTES if sequence in ROUTES[c]]


def levels(sequence: str) -> dict:
    """Engines running ``sequence`` with ideal pulses, a waveform per pulse, their own pulse set,
    or a whole sequence file (Pulseq .seq on FID-A, a sequence description on FSL-MRS)."""
    out = {'ideal': [], 'shaped': [], 'own': [], 'file': []}
    for cat, seqs in ROUTES.items():
        kinds = [set(r.pulses.values()) for r in seqs.get(sequence, [])]
        name = ENGINE_LABEL[cat]
        if any(SHAPED not in k and IDEAL in k for k in kinds):
            out['ideal'].append(name)
        if any(SHAPED in k for k in kinds):
            out['shaped'].append(name)
        if any(k == {OWN} for k in kinds):
            out['own'].append(name)
    if sequence in ('PRESS', 'sLASER', 'STEAM', 'MEGA-PRESS'):
        out['file'].append('FID-A (.seq)')
    out['file'].append('FSL-MRS (.json)')
    return out


def why_not(category: str, sequence: str) -> str:
    """Why ``category`` cannot run ``sequence``, and which engines can."""
    others = runs_on(sequence)
    where = f" Runs on: {', '.join(others)}." if others else ""
    return f"{ENGINE_SCOPE.get(category, '')}{where}"


def _matches(route: Route, backend) -> bool:
    if route.backend != backend.name:
        return False
    if route.modes and backend.current_mode not in route.modes:
        return False
    p = {**backend.optional_params, **backend.mandatory_params}
    return all(p.get(k) == v for k, v in route.sheet.items())


def current(br) -> tuple[str | None, Route | None]:
    """(sequence, route) of the active backend, (None, None) when its sheet names no sequence."""
    b = br.backend
    for seq, routes in ROUTES.get(b.category, {}).items():
        for route in routes:
            if _matches(route, b):
                return seq, route
    return None, None


def pulse_choices(category: str, sequence: str, role: str) -> set[str]:
    """IDEAL / SHAPED / OWN choices ``category`` offers for ``role`` in ``sequence``."""
    return {r.pulses.get(role) for r in ROUTES.get(category, {}).get(sequence, [])} - {None}


def choose(category: str, sequence: str, want: dict | None = None) -> Route | None:
    """The route closest to the wanted pulses (role -> IDEAL / SHAPED); without a wish, the one with
    the fewest waveforms, i.e. ideal wherever the engine allows."""
    routes = ROUTES.get(category, {}).get(sequence, [])
    if not routes:
        return None
    want = want or {}

    def cost(r):
        miss = sum(1 for role, kind in want.items() if r.pulses.get(role) not in (kind, OWN, None))
        return miss, sum(1 for k in r.pulses.values() if k == SHAPED)
    return min(routes, key=cost)


def apply(br, category: str, sequence: str, want: dict | None = None) -> Route | None:
    """Switch BasisREMY to ``sequence`` on ``category`` with the pulses closest to ``want``.
    Sheet values carry over between backends (BasisREMY.set_backend). None: not runnable here."""
    route = choose(category, sequence, want)
    if route is None:
        return None
    if br.backend.name != route.backend:
        br.set_backend(route.backend)
    b = br.backend
    if route.modes and b.current_mode not in route.modes:
        b.set_mode(route.modes[0])
    for k, v in route.sheet.items():
        (b.mandatory_params if k in b.mandatory_params or k not in b.optional_params
         else b.optional_params)[k] = v
    return route


def pulse_note(category: str, sequence: str, role: str, kind: str) -> str | None:
    """Why ``role`` cannot be ``kind`` (IDEAL / SHAPED) on ``category``; None when it can."""
    if kind in pulse_choices(category, sequence, role):
        return None
    name = ROLE_NAME[role].lower()
    eng = ENGINE_LABEL.get(category, category)
    if OWN in pulse_choices(category, sequence, role):
        return f"{eng} uses its own {name} pulse."
    if role == 'edit' and kind == IDEAL:
        return ("An editing pulse is never ideal: its frequency selectivity decides what is "
                "edited.")
    other = [ENGINE_LABEL[c] for c in ROUTES
             if c != category and kind in pulse_choices(c, sequence, role)]
    what = 'ideal' if kind == IDEAL else 'waveform'
    if other:
        where = f"; {', '.join(other)} {'has' if len(other) == 1 else 'have'} one"
    elif kind == SHAPED:
        where = "; no engine here has one (a whole sequence file can carry it)"
    else:
        where = ""
    return f"{eng} has no {what} {name} pulse for {sequence}{where}."


def fixed_timing(backend) -> str | None:
    p = {**backend.optional_params, **backend.mandatory_params}
    return (FIXED_TIMING.get((backend.name, p.get('Sequence')))
            or FIXED_TIMING.get((backend.name, None)))


#**************************************************************************************************#
#                                       recommended values                                         #
#**************************************************************************************************#
# Values a data header never holds. A recommended value is shown in its own colour with its
# source; the scan's own values (TE, field, points, bandwidth) are never recommended.
_SALEH = "Saleh et al. 2019, multi-vendor universal MEGA-PRESS"
_UNSET = (None, '', 'missing input', 'Select option')


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def recommended(backend, sequence: str | None) -> dict:
    """{key: (value, source)} for the design values the backend shows; empty without a TE where
    the value follows from it."""
    shown = backend.get_params_for_mode()
    p = {**backend.optional_params, **backend.mandatory_params}
    te = _num(p.get('TE'))
    out = {}

    def put(key, value, source):
        if key in shown and value is not None:
            out[key] = (value, source)

    if te is not None and sequence == 'PRESS':
        sym = "Symmetric: the scan's own split is not in the data header"
        put('Tau 1', te / 2, sym)
        put('Tau 2', te / 2, sym)
        put('TE2', te / 2, sym)
    if te is not None:
        for k, frac in (('sLASER TE1', 0.25), ('sLASER TE2', 0.5), ('sLASER TE3', 0.25)):
            put(k, te * frac, "Symmetric spacing (TE/4, TE/2, TE/4): the scan's own is not in "
                              "the data header")
    put('TM', 10.0, "Short mixing time commonly used; this data header holds no TM")
    if sequence and sequence.startswith(('MEGA', 'HERMES', 'HERCULES')):
        put('Edit On', 1.9, f"GABA editing at 1.9 ppm ({_SALEH})")
        put('Edit Off', 7.5, f"OFF at 7.5 ppm ({_SALEH})")
        hermes = sequence.startswith(('HERMES', 'HERCULES'))
        tp = 20.0 if hermes else 15.0
        src = (f"20 ms editing pulses, as in HERMES ({_SALEH})" if hermes
               else f"15 ms editing pulse ({_SALEH})")
        if backend.name == 'MRSCloud':      # its own editing pulses, their own duration
            return out
        if sequence == 'MEGA-sLASER' and te is not None and 'Edit Tp' in shown:
            fit = 2 * (te / 8 - 0.01)          # FSL-MRS: two editing pulses inside TE/4 gaps
            if fit < tp:
                tp, src = round(fit, 2), src + ", shortened to fit TE"
        put('Edit Tp', tp, src)
        put('Edit Bandwidth (Hz)', round(2.274 / (tp / 1e3), 1),
            f"FWHM of a {tp:g} ms Gaussian ({_SALEH})")
        if 'Edit Pulse Path' in shown:
            put('Edit Pulse Path', 'standard:gauss-edit',
                f"Gaussian editing pulse, open stand-in for the sinc-Gaussian of {_SALEH}")
    return out


def apply_recommended(backend, sequence: str | None, keep: set) -> dict:
    """Write the recommended values into the sheet, except for the keys in ``keep`` (set by the
    user or read from the file). Returns {key: source} of what is now recommended."""
    rec = recommended(backend, sequence)
    for key, (value, _) in rec.items():
        if key in keep:
            continue
        target = backend.mandatory_params if key in backend.mandatory_params else backend.optional_params
        target[key] = value
    return {k: s for k, (_, s) in rec.items() if k not in keep}


# the scan's own values: from the file or typed, never recommended
SCAN_KEYS = {'TE', 'Bfield', 'Samples', 'Bandwidth', 'Center Freq', 'Sequence', 'Localization',
             'System'}


def value_state(backend, key, value, user: set, from_file: dict, rec: dict) -> str:
    """'user' / 'file' / 'recommended' / 'default' / 'missing' for one sheet value."""
    if value in _UNSET:
        return 'missing'
    if key in user:
        return 'user'
    if key in from_file and from_file[key] == value:
        return 'file'
    if key in rec:
        return 'recommended'
    if key in SCAN_KEYS:
        return 'file' if key in from_file else 'user'
    return 'default'

####################################################################################################
#                                      fslmrs_backend.py                                           #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 18/02/26                                                                                #
#                                                                                                  #
# Purpose: FSL-MRS backend for basis set simulation using density matrix methods.                  #
#          Uses FSL-MRS Python package (no Octave required).                                       #
#                                                                                                  #
####################################################################################################


#*************#
#   imports   #
#*************#
import os
import sys
import json
import numpy as np
from basisremy.backends.base import Backend
from basisremy.core.paths import externals_root

# denmatsim lives at <externals>/fsl_mrs/fsl_mrs/denmatsim/
# add its parent to sys.path so 'from denmatsim import ...' works
_denmatsim_parent = str(externals_root() / 'fsl_mrs' / 'fsl_mrs')
if _denmatsim_parent not in sys.path:
    sys.path.insert(0, _denmatsim_parent)


#**************************************************************************************************#
#                                          FSL-MRS Backend                                         #
#**************************************************************************************************#
#                                                                                                  #
# Implements basis set simulation using FSL-MRS density matrix simulations.                        #
# Advantages:                                                                                      #
#   - Pure Python (no Octave/MATLAB required)                                                      #
#   - Quantum-mechanically accurate simulations                                                    #
#   - Supports custom pulse sequences via JSON                                                     #
#   - Edited sequences use a Gaussian editing pulse and return ON / OFF / DIFF sub-spectra          #
#                                                                                                  #
#**************************************************************************************************#
class FSLMRSBackend(Backend):
    def __init__(self):
        super().__init__()

        self.name = 'FSL-MRS'
        self.display_name = 'FSL-MRS'
        self.category = 'FSL-MRS'
        self.requires_octave = False  # Pure Python!

        # Mode support (overrides base class)
        self.modes = ['Simple', 'Template', 'Custom']
        self.current_mode = 'Simple'

        # Supported sequences
        # - PRESS, STEAM: Have both ideal and template options
        # - Others: Ideal pulses only (no templates yet)
        # - Custom: User provides complete file
        self.supported_sequences = [
            'PRESS', 'STEAM', 'LASER', 'sLASER',
            'MEGA-PRESS', 'HERMES', 'HERCULES', 'MEGA-sLASER',
            'Custom'
        ]

        # Predefined sequence files available in externals/fsl_mrs
        # These contain REAL pulse shapes for specific parameters!
        # Using them with different parameters may be inaccurate
        self.predefined_sequences = {
            'PRESS_7T': {
                'path': 'examplePRESS.json',
                'description': 'PRESS at 7T with real pulse shapes',
                'B0': 7.0,
                'TE': None,
                'Rx_Points': 4096,
                'Rx_SW': 6000,
                'notes': 'Real Siemens 7T PRESS pulse shapes'
            },
            'STEAM_7T_11ms': {
                'path': 'example.json',
                'description': '11ms STEAM at 7T with real pulse shapes',
                'B0': 6.98,
                'TE': 11,
                'Rx_Points': 8192,
                'Rx_SW': 6000,
                'notes': 'FMRIB 7T STEAM sequence with real pulse shapes'
            },
        }

        # Metabolites denmatsim can actually simulate (the 'sys*' systems in
        # denmatsim/spinSystems.json) — anything else would always fail with
        # "no spin system". TODO: Add macromolecules sim parts.
        self.default_metabolites = [
            'Ace', 'Ala', 'Asc', 'Asp', 'Cit', 'Cr', 'EtOH', 'GABA',
            'GABA_gov', 'GPC', 'GSH', 'GSH_v2', 'Glc', 'Gln', 'Glu', 'Gly',
            'H2O', 'Ins', 'Lac', 'NAA', 'NAAG', 'PCh', 'PCr', 'PE', 'Phenyl',
            'Ref0ppm', 'Scyllo', 'Ser', 'Tau', 'Tyros', 'bHB', 'bHG',
        ]

        # Metabolite selection for UI — the standard 1H-MRS set starts on.
        _default_on = {
            'Ala', 'Asc', 'Asp', 'Cr', 'GABA', 'GPC', 'GSH', 'Glc', 'Gln',
            'Glu', 'Gly', 'Ins', 'Lac', 'NAA', 'NAAG', 'PCh', 'PCr', 'PE',
            'Scyllo', 'Tau',
        }
        self.metabs = {name: name in _default_on
                       for name in self.default_metabolites}

        # Dropdown options (export-related options live in the Export dialog)
        self.dropdown = {
            'Sequence': [
                'PRESS', 'STEAM', 'LASER', 'sLASER',
                'MEGA-PRESS', 'HERMES', 'HERCULES', 'MEGA-sLASER'
            ],
            'Template File': [info['description'] for info in self.predefined_sequences.values()],
        }

        # File selection fields
        self.file_selection = ['Custom Sequence']

        # Mandatory parameters (no Output Path / Output Format — those belong
        # to the post-simulation Export dialog)
        # Scan-physics values have NO defaults — they must come from REMY or
        # the user (a prefilled 3 T / TE 35 would masquerade as file metadata).
        self.mandatory_params = {
            'Sequence': None,
            'Samples': None,
            'Bandwidth': None,
            'Bfield': None,
            'TE': None,
            'Nucleus': '1H',
            'Center Freq': None,
            'Metabolites': [],
        }

        # Optional parameters
        self.optional_params = {
            'TM': 10,
            'Template File': None,
            'Edit On': 1.9,          # ppm (GABA; 4.56 for GSH) — MEGA kinds only
            'Edit Off': 7.5,
            'Edit Tp': 14.0,         # editing pulse duration [ms], all edited kinds
            'Linewidth': 1.0,
            'Custom Sequence': None,
        }

        # Edited sequences: MEGA kinds take the editing frequencies from the
        # sheet; HERMES / HERCULES use their fixed schemes (see _SCHEMES).
        self._edited_sequences = {'MEGA-PRESS', 'HERMES', 'HERCULES', 'MEGA-sLASER'}
        self._mega_sequences = {'MEGA-PRESS', 'MEGA-sLASER'}

        # Sequences that use a mixing time (STEAM family). For PRESS,
        # LASER, sLASER, MEGA-* etc. TM is meaningless.
        self._tm_sequences = {'STEAM'}

        # Sequence selection changes which fields are relevant — the GUI
        # rebuilds the parameter panel whenever a key listed here changes.
        self.schema_affecting_keys = {'Sequence'}

    def map_sequence_in(self, seq: str) -> 'str | None':
        """Translate an arbitrary sequence name into FSL-MRS's vocabulary."""
        if not seq:
            return None
        s = seq.strip().lower()
        # Exact match
        for opt in self.dropdown.get('Sequence', []):
            if opt.lower() == s:
                return opt
        # Cross-backend synonyms
        if 'steam' in s:
            return 'STEAM'
        if 'hercules' in s:
            return 'HERCULES'
        if 'hermes' in s:
            return 'HERMES'
        if 'mega' in s and ('slaser' in s or 'semi' in s):
            return 'MEGA-sLASER'
        if 'mega' in s:
            return 'MEGA-PRESS'
        if 'slaser' in s or 'semi' in s:
            return 'sLASER'
        if 'laser' in s:
            return 'LASER'
        if 'press' in s or 'spin echo' in s or 'spinecho' in s or s == 'se':
            return 'PRESS'
        return None

    def get_params_for_mode(self, mode=None):
        """
        Return parameters to display in GUI for the given mode.
        Overrides base class.
        """
        if mode is None:
            mode = self.current_mode

        # Common parameters shown in every mode
        common = {
            'Metabolites': self.mandatory_params['Metabolites'],
        }

        if mode == 'Simple':
            seq = self.mandatory_params.get('Sequence')
            params = {
                'Sequence': seq,
                'Bfield':   self.mandatory_params['Bfield'],
                'TE':       self.mandatory_params['TE'],
                'Samples':  self.mandatory_params['Samples'],
                'Bandwidth': self.mandatory_params['Bandwidth'],
                'Nucleus': self.mandatory_params['Nucleus'],
                'Center Freq': self.mandatory_params['Center Freq'],
            }
            # Conditional fields — only show when relevant to the chosen
            # sequence so the parameter sheet stays uncluttered.
            if seq in self._tm_sequences:
                params['TM'] = self.optional_params['TM']
            if seq in self._mega_sequences:
                params['Edit On'] = self.optional_params['Edit On']
                params['Edit Off'] = self.optional_params['Edit Off']
            if seq in self._edited_sequences:
                params['Edit Tp'] = self.optional_params['Edit Tp']
            params['Linewidth'] = self.optional_params['Linewidth']
            params.update(common)
            return params

        elif mode == 'Template':
            return {
                'Template File': self.optional_params['Template File'],
                'Samples': self.mandatory_params['Samples'],
                'Bandwidth': self.mandatory_params['Bandwidth'],
                'Nucleus': self.mandatory_params['Nucleus'],
                'Center Freq': self.mandatory_params['Center Freq'],
                'Linewidth': self.optional_params['Linewidth'],
                **common,
            }

        elif mode == 'Custom':
            return {
                'Custom Sequence': self.optional_params['Custom Sequence'],
                'Samples': self.mandatory_params['Samples'],
                'Bandwidth': self.mandatory_params['Bandwidth'],
                'Nucleus': self.mandatory_params['Nucleus'],
                'Center Freq': self.mandatory_params['Center Freq'],
                **common,
            }

        return dict(self.mandatory_params)

    @staticmethod
    def _is_missing(value) -> bool:
        """Return True for GUI/REMY placeholder values."""
        if value is None:
            return True
        if isinstance(value, str):
            return value.strip().lower() in {"", "missing input", "select option"}
        return False

    @staticmethod
    def _first_raw(source: dict, *keys, default=None):
        """Return the first key present in *source*, preserving blank values."""
        for key in keys:
            if key in source:
                return source[key]
        return default

    @classmethod
    def _coerce_number_or_blank(cls, value, *, default=None, as_int=False):
        """Coerce numeric REMY values, preserving blanks that need user input."""
        if cls._is_missing(value):
            return "" if default is None else default
        try:
            number = float(value)
        except (TypeError, ValueError):
            return value
        return int(number) if as_int else number

    def parseREMY(self, MRSinMRS):
        """
        Parse REMY output to FSL-MRS backend parameters

        Args:
            MRSinMRS: Dictionary of parameters extracted by REMY

        Returns:
            tuple: (mandatory_params_dict, optional_params_dict)
        """
        params = {}
        opt = {}

        # Required parameters
        params['Samples'] = self._coerce_number_or_blank(
            self._first_raw(
                MRSinMRS,
                'NumberOfDatapoints',
                'Samples',
                default=self.mandatory_params['Samples'],
            ),
            as_int=True,
        )
        params['Bandwidth'] = self._coerce_number_or_blank(
            self._first_raw(
                MRSinMRS,
                'SpectralWidth',
                'Bandwidth',
                default=self.mandatory_params['Bandwidth'],
            )
        )
        params['Bfield'] = self._coerce_number_or_blank(
            self._first_raw(
                MRSinMRS,
                'B0',
                'Bfield',
                default=self.mandatory_params['Bfield'],
            )
        )
        params['TE'] = self._coerce_number_or_blank(
            self._first_raw(MRSinMRS, 'TE', default=self.mandatory_params['TE'])
        )
        params['Nucleus'] = self._first_raw(
            MRSinMRS, 'Nucleus', default=self.mandatory_params['Nucleus'],
        )

        # Calculate center frequency from field strength if not provided
        center_freq = self._first_raw(MRSinMRS, 'Center Freq')
        if not self._is_missing(center_freq):
            params['Center Freq'] = self._coerce_number_or_blank(center_freq)
        elif not self._is_missing(params['Bfield']):
            # Calculate for 1H at given field strength
            # gamma_1H = 42.577 MHz/T
            params['Center Freq'] = 42.577 * float(params['Bfield'])
        else:
            params['Center Freq'] = ""

        # Sequence detection
        sequence = self.parseProtocol(MRSinMRS.get('Protocol', ''))
        if sequence and sequence in self.supported_sequences:
            params['Sequence'] = sequence
        else:
            params['Sequence'] = None  # Will need manual selection

        params['Metabolites'] = []

        # Optional parameters
        opt['Linewidth'] = 1.0
        opt['Custom Sequence'] = None

        return params, opt

    def show_predefined_sequences(self):
        """
        Display information about available predefined sequence files
        """
        print("\n" + "="*80)
        print("Available Predefined Sequences (with REAL pulse shapes)")
        print("="*80)

        for key, info in self.predefined_sequences.items():
            print(f"\n{key}:")
            print(f"  Description: {info['description']}")
            print(f"  Field Strength: {info['B0']} T")
            print(f"  TE: {info['TE']} ms" if info['TE'] else "  TE: Depends on delays")
            print(f"  Acquisition: {info['Rx_Points']} points @ {info['Rx_SW']} Hz")
            print(f"  Notes: {info['notes']}")
            print(f"  File: {info['path']}")

        print("\n" + "="*80)
        print("⚠️  These files contain real RF pulse waveforms!")
        print("   Using them with significantly different parameters may be inaccurate.")
        print("   For other field strengths/TEs, provide a custom sequence JSON file.")
        print("="*80 + "\n")

    def parseProtocol(self, protocol):
        """
        Parse sequence name from protocol string

        Args:
            protocol: Protocol string from scanner

        Returns:
            str: Standardized sequence name or None
        """
        return self.map_sequence_in(protocol)

    def _coerce_params(self, params: dict) -> dict:
        """Return a copy of *params* with all numeric fields cast to their
        proper Python types.

        GUI entries arrive as strings; arithmetic inside _generate_sequence_json
        and run_simulation will crash with ``TypeError: unsupported operand type
        'str'`` without this guard. We coerce conservatively — only keys we
        know are numeric — and leave everything else untouched.
        """
        p = dict(params)
        float_keys = ('TE', 'Bfield', 'Bandwidth', 'TM', 'Edit On', 'Edit Off',
                      'Edit Tp', 'Linewidth', 'Center Freq')
        int_keys   = ('Samples',)
        for k in float_keys:
            if k in p and not self._is_missing(p[k]):
                try:
                    p[k] = float(p[k])
                except (TypeError, ValueError):
                    pass
        for k in int_keys:
            if k in p and not self._is_missing(p[k]):
                try:
                    p[k] = int(float(p[k]))
                except (TypeError, ValueError):
                    pass
        return p

    def _load_predefined(self, predefined_info, params):
        """Load a predefined sequence JSON and apply the user's acquisition
        overrides (Rx points / bandwidth — never the pulse shapes)."""
        seq_rel_path = predefined_info['path']
        # JSON example files live inside denmatsim itself
        denmatsim_path = str(externals_root() / 'fsl_mrs' / 'fsl_mrs' / 'denmatsim')
        seq_file_path = os.path.join(denmatsim_path, seq_rel_path)
        if not os.path.exists(seq_file_path):
            raise RuntimeError(
                f"Predefined sequence file not found: {seq_file_path}\n"
                f"Make sure FSL-MRS submodule is initialized:\n"
                f"  git submodule update --init --recursive"
            )
        print(f"  File: {seq_rel_path}")
        print(f"  Parameters: B0={predefined_info['B0']}T, TE={predefined_info['TE']}ms, "
              f"Points={predefined_info['Rx_Points']}, BW={predefined_info['Rx_SW']}Hz")
        print("  Note: This file contains REAL pulse shapes - highly accurate!")
        with open(seq_file_path, 'r') as f:
            seq_params = json.load(f)

        # The template's own echo time (pulse centre to ADC) must match the
        # sheet: its pulse shapes and delays are fixed, so a different TE
        # would silently simulate the template's TE.
        te_file = predefined_info['TE']
        if te_file is None:
            te_file = self._template_te_ms(seq_params)
        te_req = params.get('TE')
        if not self._is_missing(te_req) and abs(float(te_req) - te_file) > 1.0:
            raise ValueError(
                f"FSL-MRS: the template '{predefined_info['description']}' has a "
                f"fixed TE of {te_file:.1f} ms, but TE is {float(te_req):g} ms. "
                f"Set TE to {te_file:.1f}, or use Simple mode for an ideal-pulse "
                f"basis at {float(te_req):g} ms.")

        # Update ONLY acquisition parameters (not pulse shapes!)
        # User can override Rx_Points and Rx_SW for their specific needs
        if params['Samples'] != predefined_info['Rx_Points']:
            print(f"  ⚠️  Updating Rx_Points from {predefined_info['Rx_Points']} to {params['Samples']}")
            seq_params['Rx_Points'] = params['Samples']

        if params['Bandwidth'] != predefined_info['Rx_SW']:
            print(f"  ⚠️  Updating Rx_SW from {predefined_info['Rx_SW']} to {params['Bandwidth']}")
            seq_params['Rx_SW'] = params['Bandwidth']

        return seq_params

    @staticmethod
    def _template_te_ms(seq_params):
        """Echo time of a spin-echo type template: excitation pulse centre to
        the ADC (delays run from pulse end to next pulse start)."""
        times = [float(r['time']) for r in seq_params['RF']]
        return (sum(float(d) for d in seq_params['delays']) + sum(times)
                - times[0] / 2.0) * 1e3

    # Editing schemes: sub-experiment -> editing target(s) in ppm, applied to
    # both editing pulses (two targets = dual-lobe pulse). HERMES / HERCULES
    # mirror MRSCloud's order A-D; the differences are combined as the
    # MRSCloud backend does (DIFF1 = GABA, DIFF2 = GSH).
    _SCHEMES = {
        'HERMES':   {'A': (4.56,), 'B': (1.90,), 'C': (4.56, 1.90), 'D': (7.50,)},
        'HERCULES': {'A': (4.58,), 'B': (4.18,), 'C': (4.58, 1.90), 'D': (4.18, 1.90)},
    }
    # MEGA-PRESS pulse-centre spacings as fractions of TE (FID-A's Siemens
    # MEGA-PRESS timing: excite - 180 - edit - 180 - edit - ADC, editing pulses
    # at TE/4 and 3TE/4). FID-A lists the set for "TE = 68" but it sums to
    # 69.0 ms, so the fractions - not the ms values - are what is kept.
    _MEGA_TAU_FRACTIONS = tuple(t / 69.0001 for t in (4.545, 12.7025, 21.7975, 12.7025, 17.2526))
    _EDIT_PULSE_POINTS = 400

    @staticmethod
    def _editing_pulse(targets_ppm, tp_ms, bfield, central_shift, n=400):
        """Gaussian 180-degree editing pulse (single- or dual-lobe) as a
        denmatsim RF entry. The frequency modulation is written into the
        waveform's phase, so denmatsim sees an on-carrier pulse: its
        ``frequencyOffset`` shortcut shifts the rotating frame for the pulse
        duration without unwinding the phase afterwards, and a 10 us
        "ideal" pulse at an offset is a hard pulse on every spin."""
        tp = float(tp_ms) / 1000.0
        t = (np.arange(n) + 0.5) / n * tp
        env = np.exp(-4.0 * np.log(100.0) * ((t - tp / 2.0) / tp) ** 2)   # 1 % at the edges
        env = env * (0.5 / (env.sum() * tp / n))                            # integral 0.5 cycles = 180 deg
        wave = np.zeros(n, dtype=complex)
        for ppm in targets_ppm:
            f_hz = (float(ppm) - central_shift) * bfield * 42.577
            wave += env * np.exp(2j * np.pi * f_hz * (t - tp / 2.0))
        return {'time': tp, 'frequencyOffset': 0, 'phaseOffset': 0,
                'amp': np.abs(wave).tolist(), 'phase': np.angle(wave).tolist(),
                'grad': [0, 0, 0]}

    @staticmethod
    def _delays_from_taus(taus_s, rf):
        """denmatsim delays (pulse end to next pulse start) from pulse-centre
        spacings; the last spacing runs from the last pulse centre to the ADC."""
        durations = [float(r['time']) for r in rf]
        delays = []
        for i, tau in enumerate(taus_s):
            after = durations[i + 1] / 2.0 if i + 1 < len(durations) else 0.0
            d = tau - durations[i] / 2.0 - after
            if d < 0:
                raise ValueError(
                    "FSL-MRS: the editing pulse does not fit the timing - shorten "
                    "'Edit Tp' or lengthen TE.")
            delays.append(d)
        return delays

    def _generate_sequence_json(self, params, edit_ppm=None):
        """
        Generate FSL-MRS sequence JSON with IDEAL PULSES (FID-A style)

        ``edit_ppm`` (edited sequences only) is the tuple of editing targets
        in ppm for this sub-experiment; both editing pulses are Gaussian
        180s of ``Edit Tp`` ms at those frequencies.

        Uses instantaneous rotation operators (~10 μs) for all pulses.
        These are mathematically rigorous and standard in NMR simulation.

        For accurate simulations with real pulse shapes:
          - Use Template Mode (predefined files)
          - Use Custom Mode (provide your own JSON)

        Args:
            params: Dictionary of simulation parameters

        Returns:
            dict: Sequence definition in FSL-MRS JSON format with ideal pulses
        """
        sequence = params['Sequence']
        te = float(params['TE'])
        bandwidth = float(params['Bandwidth'])
        samples = int(float(params['Samples']))
        bfield = float(params['Bfield'])

        print(f"Generating IDEAL pulse sequence: {sequence}")
        print("  Using instantaneous rotation operators (FID-A style)")
        print("  Perfect flip angles, no realistic pulse effects")
        print("  For accurate simulations, use Template or Custom mode")

        # Rotating-frame carrier: spins evolve (and RF frequencyOffsets are
        # measured) relative to this chemical shift.
        central_shift = 4.65  # ppm

        # Base sequence structure for denmatsim
        seq_def = {
            'sequenceName': f'{sequence}_ideal',
            'description': f'Ideal {sequence} with instantaneous pulses',
            'B0': bfield,
            'centralShift': central_shift,  # ppm - typical for 1H MRS
            'Rx_Points': samples,
            'Rx_SW': bandwidth,
            'Rx_LW': 1.0 if self._is_missing(params.get('Linewidth'))
                     else float(params['Linewidth']),
            # denmatsim's FID starts at -90 deg for an on-resonance spin;
            # this receiver phase puts singlets on the real axis, so no
            # per-metabolite phasing is needed (which would flip inverted
            # multiplets such as lactate at TE 144 upright).
            'Rx_Phase': -1.5708,
            'x': [-15, 15],
            'y': [-15, 15],
            'z': [-15, 15],
            'resolution': [8, 8, 8],
            'RFUnits': 'Hz',
            'GradUnits': 'mT',
            'spaceUnits': 'mm',
        }

        # Ideal pulse duration (essentially instantaneous)
        ideal_pulse_duration = 0.00001  # 10 microseconds

        # For ideal pulses: flip_angle = 2π * B1_Hz * duration
        # So B1_Hz = flip_angle_rad / (2π * duration)
        import math
        amp_90 = (math.pi / 2) / (2 * math.pi * ideal_pulse_duration)   # ~25000 Hz
        amp_180 = math.pi / (2 * math.pi * ideal_pulse_duration)         # ~50000 Hz

        # Define sequences
        if sequence == 'PRESS':
            # PRESS: 90° - TE/4 - 180° - TE/2 - 180° - TE/4 - ACQ (between pulse
            # centres; denmatsim delays run from pulse end to next pulse start,
            # so the pulse length is taken off). With TE/2 - TE/2 - TE/2 the
            # echo formed at 1.5 x TE (verified against FID-A on Lac / Glu).
            # 3 RF pulses = 3 delays, 3 rephaseAreas, 3 CoherenceFilter
            te_s = te / 1000.0
            d_quarter = te_s / 4.0 - ideal_pulse_duration
            d_half = te_s / 2.0 - ideal_pulse_duration
            d_last = te_s / 4.0 - ideal_pulse_duration / 2.0
            seq_def.update({
                'RF': [
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_90], 'phase': [0], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                ],
                'delays': [d_quarter, d_half, d_last],
                'rephaseAreas': [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
                'CoherenceFilter': [-1, 1, -1],
            })

        elif sequence == 'STEAM':
            # STEAM: 90° - TE/2 - 90° - TM - 90° - TE/2 - ACQ
            # 3 RF pulses = 3 delays, 3 rephaseAreas, 3 CoherenceFilter
            tm_raw = params.get('TM')
            tm = 10.0 if self._is_missing(tm_raw) else float(tm_raw)  # keeps an explicit 0
            seq_def.update({
                'RF': [
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_90], 'phase': [0], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_90], 'phase': [3.14159], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_90], 'phase': [0], 'grad': [0, 0, 0]},
                ],
                'delays': [te/2000 - ideal_pulse_duration, tm/1000 - ideal_pulse_duration,
                           te/2000 - ideal_pulse_duration/2],
                'rephaseAreas': [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
                'CoherenceFilter': [1, 0, -1],
            })

        elif sequence == 'LASER':
            # LASER: 90° + 3 pairs of 180° AFP = 7 RF pulses
            # 7 RF = 7 delays, 7 rephaseAreas, 7 CoherenceFilter
            # 90 - tau/2 - 180 - tau - 180 - ... - 180 - tau/2 - ACQ with
            # tau = TE/6 (six refocusing pulses, echo at the ADC). Seven equal
            # TE/6 delays put the echo TE/12 late and made TE 7/6 too long.
            tau = te / 6000.0  # TE/6 in seconds
            d = tau - ideal_pulse_duration
            d_start = tau / 2.0 - ideal_pulse_duration
            d_end = tau / 2.0 - ideal_pulse_duration / 2.0
            seq_def.update({
                'RF': [
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_90], 'phase': [0], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                ],
                'delays': [d_start, d, d, d, d, d, d_end],
                'rephaseAreas': [[0, 0, 0]] * 7,
                'CoherenceFilter': [-1, 1, -1, 1, -1, 1, -1],
            })

        elif sequence == 'sLASER':
            # sLASER: 90° + 2 pairs of 180° AFP = 5 RF pulses
            # 5 RF = 5 delays, 5 rephaseAreas, 5 CoherenceFilter
            # 90 - TE/8 - 180 - TE/4 - 180 - TE/4 - 180 - TE/4 - 180 - TE/8 - ACQ
            # (five equal TE/4 delays made the effective TE 5/4 too long).
            te_s = te / 1000.0
            d = te_s / 4.0 - ideal_pulse_duration
            d_start = te_s / 8.0 - ideal_pulse_duration
            d_end = te_s / 8.0 - ideal_pulse_duration / 2.0
            seq_def.update({
                'RF': [
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_90], 'phase': [0], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]},
                ],
                'delays': [d_start, d, d, d, d_end],
                'rephaseAreas': [[0, 0, 0]] * 5,
                'CoherenceFilter': [-1, 1, -1, 1, -1],
            })

        elif sequence in self._edited_sequences:
            if not edit_ppm:
                raise ValueError(f"FSL-MRS: {sequence} needs editing targets (edit_ppm).")
            tp = params.get('Edit Tp')
            tp = self.optional_params['Edit Tp'] if self._is_missing(tp) else float(tp)
            edit = self._editing_pulse(edit_ppm, tp, bfield, central_shift,
                                       self._EDIT_PULSE_POINTS)
            exc = {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                   'amp': [amp_90], 'phase': [0], 'grad': [0, 0, 0]}
            ref = {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                   'amp': [amp_180], 'phase': [1.5708], 'grad': [0, 0, 0]}
            te_s = te / 1000.0
            if sequence == 'MEGA-sLASER':
                # 90 - TE/8 - AFP - TE/8 - edit - TE/8 - AFP - TE/4 - AFP - TE/8 -
                # edit - TE/8 - AFP - TE/8 - ADC: sLASER refocusing at odd
                # eighths of TE, editing pulses at TE/4 and 3TE/4 (a TE/2 apart).
                rf = [exc, ref, edit, ref, ref, edit, ref]
                taus = [te_s / 8, te_s / 8, te_s / 8, te_s / 4, te_s / 8, te_s / 8, te_s / 8]
                cfilter = [-1, 1, None, -1, 1, None, -1]
            else:
                # PRESS localisation with the MEGA-PRESS timing (MEGA-PRESS,
                # HERMES, HERCULES): 90 - 180 - edit - 180 - edit - ADC.
                rf = [exc, ref, edit, ref, edit]
                taus = [f * te_s for f in self._MEGA_TAU_FRACTIONS]
                cfilter = [-1, 1, None, -1, None]
            # No coherence filter after an editing pulse: it is selective, so
            # the coherence order of the untouched spins must survive.
            seq_def.update({
                'RF': rf,
                'delays': self._delays_from_taus(taus, rf),
                'rephaseAreas': [[0, 0, 0]] * len(rf),
                'CoherenceFilter': cfilter,
            })
            print(f"  {sequence}: Gaussian editing pulse {tp:g} ms at "
                  f"{', '.join(f'{p:g}' for p in edit_ppm)} ppm")

        else:
            # Generic single pulse excitation for unknown sequences
            seq_def.update({
                'RF': [
                    {'time': ideal_pulse_duration, 'frequencyOffset': 0, 'phaseOffset': 0,
                     'amp': [amp_90], 'phase': [0], 'grad': [0, 0, 0]},
                ],
                'delays': [te / 1000.0],
                'rephaseAreas': [[0, 0, 0]],
                'CoherenceFilter': [-1],
            })

        return seq_def

    def run_simulation(self, params, progress_callback=None, stop_event=None):
        """
        Run FSL-MRS basis set simulation using denmatsim

        Args:
            params: Dictionary of simulation parameters
            progress_callback: Optional callback function(current, total)

        Returns:
            dict: {metabolite_name: FID_array}
        """
        # Import denmatsim (path already set up at module level)
        # Fetch FSL-MRS on first use (no-op in a source checkout).
        from basisremy.core.externals import ensure
        ensure('fsl_mrs')
        try:
            from denmatsim import simseq, utils as simutils
        except ImportError:
            # An import attempted during/before the one-time fetch leaves a
            # stale namespace module in sys.modules ("unknown location") that
            # importlib.invalidate_caches() cannot clear — purge and retry.
            import importlib
            for _name in [m for m in sys.modules
                          if m == 'denmatsim' or m.startswith('denmatsim.')]:
                del sys.modules[_name]
            importlib.invalidate_caches()
            try:
                from denmatsim import simseq, utils as simutils
            except ImportError as e:
                raise RuntimeError(
                    f"denmatsim not found at {_denmatsim_parent}/denmatsim/\n"
                    f"Error: {e}\n\n"
                    f"If the one-time FSL-MRS download just finished, restart "
                    f"BasisREMY. In a source checkout, run:\n"
                    f"  git submodule update --init --recursive"
                )

        # Coerce all numeric fields from string (GUI entries) to float/int
        # before any arithmetic — prevents "unsupported operand type 'str'" crashes.
        params = self._coerce_params(params)

        # Fail with a clear message instead of a TypeError deep inside the
        # sequence matching when a required numeric never arrived (headless
        # use). Template/Custom sequence files carry their own B0/TE, so only
        # the acquisition grid is required on those paths.
        required = ['Samples', 'Bandwidth']
        using_template = (self.current_mode == 'Template'
                          and not self._is_missing(params.get('Template File')))
        using_custom = (self.current_mode == 'Custom'
                        and params.get('Custom Sequence')
                        and os.path.exists(params['Custom Sequence']))
        if not (using_template or using_custom):
            required += ['Bfield', 'TE']
        missing = [k for k in required
                   if self._is_missing(params.get(k))
                   or not isinstance(params.get(k), (int, float))]
        if missing:
            raise ValueError(
                f"FSL-MRS: missing or non-numeric parameter(s): {', '.join(missing)}")

        print("✓ denmatsim imported successfully")

        print(f"\n{'='*80}")
        print("Running FSL-MRS simulation (denmatsim)")
        print(f"{'='*80}")
        print(f"Sequence: {params['Sequence']}")
        print(f"TE: {params['TE']} ms")
        print(f"Field strength: {params['Bfield']} T")
        print(f"Metabolites: {len(params['Metabolites'])}")

        # Allocate an internal scratch directory for the sequence JSON dump
        # (final exports go through core.exporters via the GUI Export dialog)
        output_path = self.ensure_workdir()
        os.makedirs(output_path, exist_ok=True)

        # Get or generate sequence JSON. The Custom branch is gated on the
        # mode: a stale 'Custom Sequence' pick must not override the sequence
        # configured in Simple or Template mode.
        if (self.current_mode == 'Custom'
                and params.get('Custom Sequence')
                and os.path.exists(params['Custom Sequence'])):
            # User provided custom sequence file
            print(f"Using custom sequence: {params['Custom Sequence']}")
            try:
                with open(params['Custom Sequence'], 'r') as f:
                    seq_params = json.load(f)
            except (json.JSONDecodeError, UnicodeDecodeError) as exc:
                raise ValueError(
                    f"FSL-MRS: 'Custom Sequence' must be an FSL-MRS sequence "
                    f"description (JSON, like fsl_mrs/denmatsim/examplePRESS.json); "
                    f"'{os.path.basename(params['Custom Sequence'])}' is not JSON "
                    f"({exc}). RF waveform files (.pta / .RF) belong to the FID-A "
                    f"shaped backends, not here.") from exc
            for key in ('RF', 'delays', 'CoherenceFilter'):
                if key not in seq_params:
                    raise ValueError(
                        f"FSL-MRS: 'Custom Sequence' JSON has no '{key}' entry — "
                        f"it is not an FSL-MRS sequence description.")

        elif (self.current_mode == 'Template'
              and not self._is_missing(params.get('Template File'))):
            # User explicitly chose a template — honor it (no Sequence-based
            # guessing, no B0/TE second-guessing of an explicit choice).
            choice = params['Template File']
            predefined_info = next(
                (info for info in self.predefined_sequences.values()
                 if info['description'] == choice or info['path'] == choice),
                None)
            if predefined_info is None:
                raise RuntimeError(
                    f"Unknown template '{choice}'. Available: "
                    f"{[i['description'] for i in self.predefined_sequences.values()]}")
            print(f"✓ Using selected template: {predefined_info['description']}")
            try:
                bf = float(params.get('Bfield'))
            except (TypeError, ValueError):
                bf = None
            if bf is not None and abs(bf - predefined_info['B0']) > 0.5:
                raise ValueError(
                    f"FSL-MRS: the template '{predefined_info['description']}' is a "
                    f"{predefined_info['B0']} T sequence with real pulse shapes, but "
                    f"Bfield is {bf} T. Its basis would be computed at "
                    f"{predefined_info['B0']} T and plotted/exported on a {bf} T "
                    f"axis. Set Bfield (and Center Freq) to {predefined_info['B0']} T, "
                    f"or use Simple mode for an ideal-pulse basis at {bf} T.")
            seq_params = self._load_predefined(predefined_info, params)
            lw = params.get('Linewidth')
            if not self._is_missing(lw):
                try:
                    seq_params['Rx_LW'] = float(lw)
                except (TypeError, ValueError):
                    print(f"⚠️  Ignoring non-numeric Linewidth {lw!r}; keeping the "
                          f"template's {seq_params.get('Rx_LW')} Hz.")

        else:
            # Simple mode: ideal pulses at the requested TE. (Earlier versions
            # silently swapped in the 7 T real-pulse templates for PRESS /
            # STEAM near 7 T, ignoring the requested TE.)
            seq_params = None

        # Sub-experiments: one sequence per editing condition, or a single
        # unlabelled run for everything else.
        sequence = params['Sequence']
        if seq_params is not None:
            variants = {None: seq_params}
        elif sequence in self._mega_sequences:
            on, off = params.get('Edit On'), params.get('Edit Off')
            on = self.optional_params['Edit On'] if self._is_missing(on) else float(on)
            off = self.optional_params['Edit Off'] if self._is_missing(off) else float(off)
            variants = {'ON': self._generate_sequence_json(params, (on,)),
                        'OFF': self._generate_sequence_json(params, (off,))}
        elif sequence in self._SCHEMES:
            variants = {label: self._generate_sequence_json(params, targets)
                        for label, targets in self._SCHEMES[sequence].items()}
        else:
            variants = {None: self._generate_sequence_json(params)}

        # Save the sequence file(s) for reference
        for label, seq_params in variants.items():
            suffix = f'_{label}' if label else ''
            seq_file = os.path.join(output_path, f'{sequence}{suffix}_sequence.json')
            with open(seq_file, 'w') as f:
                json.dump(seq_params, f, indent=2)
            print(f"Saved sequence file: {seq_file}")

        # Load spin systems for metabolites
        try:
            spinSystems = simutils.readBuiltInSpins()
            print(f"✓ Loaded {len(spinSystems)} built-in spin systems")
        except Exception as e:
            # A zero-filled placeholder would render as a "successful" basis
            # of flat traces — fail loudly instead.
            raise RuntimeError(f"Could not load FSL-MRS spin systems: {e}")

        # Run simulation for each metabolite (and editing condition)
        basis_set = {}
        self.last_failures = {}   # metab -> reason, surfaced by the GUI
        total_metabs = len(params['Metabolites'])

        def simulate(spin_system, seq_params):
            # denmatsim spin systems are lists of sub-spin-systems
            # (e.g. NAA has acetyl + aspartate groups): sum their FIDs
            subs = spin_system if isinstance(spin_system, list) else [spin_system]
            fid = None
            for sub_sys in subs:
                sub_fid, _ax, _pmat = simseq.simseq(sub_sys, seq_params, verbose=False)
                sub_fid = sub_fid * sub_sys.get('scaleFactor', 1.0)
                fid = sub_fid if fid is None else fid + sub_fid
            return fid

        for idx, metab in enumerate(params['Metabolites'], 1):
            if stop_event and stop_event.is_set():
                print(f"  Stopped before simulating {metab} (user cancelled).")
                break
            if progress_callback:
                progress_callback(idx, total_metabs)

            print(f"\n[{idx}/{total_metabs}] Simulating {metab}...")

            sys_name = f'sys{metab}'
            if sys_name not in spinSystems:
                print(f"  Spin system '{sys_name}' not found, skipping")
                self.last_failures[metab] = f"no spin system '{sys_name}'"
                continue

            try:
                fids = {label: simulate(spinSystems[sys_name], seq)
                        for label, seq in variants.items()}
            except Exception as e:
                print(f"  Simulation failed: {e}")
                import traceback
                traceback.print_exc()
                self.last_failures[metab] = str(e)
                continue

            if None in fids:
                basis_set[metab] = fids[None]
            else:
                for label, fid in fids.items():
                    basis_set[f'{metab} ({label})'] = fid
                if sequence in self._mega_sequences:
                    basis_set[f'{metab} (DIFF)'] = fids['ON'] - fids['OFF']
                else:
                    a, b, c, d = (fids[k] for k in 'ABCD')
                    basis_set[f'{metab} (SUM)'] = a + b + c + d
                    basis_set[f'{metab} (DIFF1)'] = (b + c) - (a + d)
                    basis_set[f'{metab} (DIFF2)'] = (a + c) - (b + d)
            print(f"  Generated {len(fids)} FID(s) with {len(next(iter(fids.values())))} points")

        print(f"\n{'='*80}")
        print("Simulation complete!")
        print(f"Generated {len(basis_set)}/{total_metabs} metabolite spectra")
        print(f"Output location: {output_path}")
        print(f"{'='*80}\n")

        return basis_set

    def _save_lcmodel_raw(self, fid, filepath, params):
        """Save FID in LCModel RAW format"""
        with open(filepath, 'w') as f:
            # LCModel RAW header
            f.write(" $NMID\n")
            f.write(" ID='BasisREMY FSL-MRS simulation'\n")
            f.write(" FMTDAT='(2E15.6)'\n")
            f.write(" VOLUME=1.0\n")
            f.write(" TRAMP=1.0\n")
            f.write(" $END\n")

            # Write FID data (real, imag pairs)
            for point in fid:
                f.write(f" {point.real:15.6E} {point.imag:15.6E}\n")

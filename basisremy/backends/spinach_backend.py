####################################################################################################
#                                        spinach_backend.py                                        #
####################################################################################################
#                                                                                                  #
# Purpose: Spinach backends — ideal-pulse Spin Echo / PRESS / STEAM / LASER, shaped PRESS and      #
#          shaped semi-LASER on a spatial grid, simulated by Spinach (Kuprov; MIT) through the     #
#          Octave runtime FID-A and MRSCloud use, on FID-A's spin-system definitions, so results   #
#          compare point by point with the FID-A backends. Spinach is written for MATLAB R2024b;   #
#          it runs under Octave with the source patch adapters/backends/spinach_octave.patch and   #
#          the shims in adapters/backends/spinach_shims/ (see core/externals.py for the sparse     #
#          fetch and the patch step). One Octave adapter serves every entry: spinach_run.m.        #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import os

import numpy as np

from basisremy.backends.base import Backend
from basisremy.backends.fida_backends import _DEFAULT_FIDA_METABS, _shaped_params


class _SpinachRuntime(Backend):
    """Shared plumbing of the Spinach backends: Octave paths, argument checks,
    the per-metabolite run loop. Subclasses set the schema and `_build_args`."""

    def __init__(self):
        super().__init__()
        self.category = 'Spinach'
        self.requires_octave = True
        self.metabs = dict(_DEFAULT_FIDA_METABS)
        self.optional_params = {'Nucleus': '1H', 'TR': None}

    def _refresh_metab_list(self):
        self.mandatory_params['Metabolites'] = [k for k, v in self.metabs.items() if v]

    def map_sequence_in(self, seq: str) -> 'str | None':
        return self.parseProtocol(seq)

    def parseProtocol(self, protocol):
        return None

    def parseREMY(self, MRSinMRS):
        mandatory = {
            'Samples':   MRSinMRS.get('NumberOfDatapoints', None),
            'Bandwidth': MRSinMRS.get('SpectralWidth', None),
            'Bfield':    MRSinMRS.get('B0', None),
            'TE':        MRSinMRS.get('TE', None),
        }
        if 'Sequence' in self.mandatory_params:
            mandatory['Sequence'] = self.parseProtocol(MRSinMRS.get('Protocol', None))
        optional = {
            'Nucleus':         MRSinMRS.get('Nucleus', None),
            'TR':              MRSinMRS.get('TR', None),
            'Model':           MRSinMRS.get('Model', None),
            'SoftwareVersion': MRSinMRS.get('SoftwareVersion', None),
            'BodyPart':        MRSinMRS.get('BodyPart', None),
        }
        return mandatory, optional

    # -------------------------------------------------- Octave
    def setup_octave_paths(self, octave=None):
        octave = octave or self.octave
        if octave is None:
            raise RuntimeError("Octave not initialized.")
        from basisremy.core.externals import ensure
        from basisremy.core.paths import octave_adapters_base
        ensure('spinach')                     # sparse clone + Octave patch (one-time)
        ensure('fidA')                        # the spin-system definitions, rf_scaleGrad
        adapters_base = octave_adapters_base(octave)
        octave.eval("warning('off', 'all');")
        octave.addpath(adapters_base + '/backends/')
        self._paths = ('./externals/spinach',
                       adapters_base + '/backends/spinach_shims',
                       './externals/fidA')

    # -------------------------------------------------- helpers
    @staticmethod
    def _num(params, key):
        v = params.get(key)
        if v is None or (isinstance(v, str) and not v.strip()):
            raise ValueError(f"Spinach: '{key}' is required but empty.")
        return float(v)

    @staticmethod
    def _make_relative(path):
        """Absolute -> relative path (the Docker container sees the project tree)."""
        if path and isinstance(path, str) and os.path.isabs(path):
            try:
                return os.path.relpath(path)
            except ValueError:
                return path.replace('\\', '/')
        return path

    def _stage_pulse(self, params):
        pulse_src = params.get('Path to Pulse')
        if not pulse_src:
            raise ValueError(f"{self.name}: 'Path to Pulse' is required (refocusing waveform).")
        return self._make_relative(self._stage_into_workdir(pulse_src))

    def _build_args(self, params):
        """[kind, *numbers] handed to spinach_run.m after the three paths."""
        raise NotImplementedError

    # -------------------------------------------------- driver
    def run_simulation(self, params, progress_callback=None, stop_event=None):
        if self.octave is None:
            print("Initializing Octave runtime...")
            self.initialize_octave(prefer_docker=True)
        self.setup_octave_paths()

        args = self._build_args(params)

        def work(octave, metab):
            fid_re, fid_im, _npts, _sw, _cf = octave.feval(
                'spinach_run', metab, args[0], *self._paths, *args[1:], nout=5,
            )
            return (np.asarray(fid_re, dtype=float).flatten()
                    + 1j * np.asarray(fid_im, dtype=float).flatten())

        # metabolites run concurrently over several Octave processes
        return self.simulate_in_parallel(params.get('Metabolites') or [], work,
                                         progress_callback, stop_event)


class SpinachBackend(_SpinachRuntime):
    """Spinach through Octave: ideal-pulse sequences on FID-A's spin systems."""

    _KINDS = {'Spin Echo': 'spinecho', 'PRESS': 'press', 'STEAM': 'steam', 'LASER': 'laser'}

    def __init__(self):
        super().__init__()
        self.name = 'Spinach'
        self.display_name = 'Ideal (SE / PRESS / STEAM / LASER)'
        self.dropdown = {'Sequence': list(self._KINDS)}
        self.mandatory_params = {
            'Sequence':    None,
            'Samples':     None,
            'Bandwidth':   None,
            'Bfield':      None,
            'Linewidth':   1,
            'TE':          None,
            'TM':          10,
            'Metabolites': [],
        }
        # TM is only shown for STEAM — rebuild the panel on Sequence changes.
        self.schema_affecting_keys = {'Sequence'}
        self._refresh_metab_list()

    def get_params_for_mode(self, mode=None):
        params = dict(self.mandatory_params)
        if params.get('Sequence') != 'STEAM':
            params.pop('TM', None)
        return params

    def parseProtocol(self, protocol):
        p = str(protocol or '').lower()
        if not p:
            return None
        if 'press' in p or 'unedited' in p:
            return 'PRESS'
        if 'steam' in p:
            return 'STEAM'
        if 'laser' in p:
            return 'LASER'
        if 'spin' in p or 'se' in p:
            return 'Spin Echo'
        return None

    def _build_args(self, params):
        seq = params.get('Sequence')
        if seq not in self._KINDS:
            raise ValueError(
                f"Spinach: unrecognised Sequence '{seq}'. "
                f"Valid: {', '.join(self._KINDS)}.")
        tm = self._num(params, 'TM') if seq == 'STEAM' else 0.0
        return [
            self._KINDS[seq],
            self._num(params, 'Samples'),
            self._num(params, 'Bandwidth'),
            self._num(params, 'Bfield'),
            float(params.get('Linewidth') or 1),
            self._num(params, 'TE'),
            tm,
        ]


class SpinachPressShaped(_SpinachRuntime):
    """PRESS with a shaped refocusing pulse on a spatial grid — the fast method
    (Zhang 2017) with Spinach's shaped_pulse_xy for the pulses."""

    def __init__(self):
        super().__init__()
        self.name = 'SpinachPressShaped'
        self.display_name = 'PRESS shaped'
        self.file_selection = ['Path to Pulse']
        self.mandatory_params = _shaped_params({'Tau 1': None, 'Tau 2': None})
        self._refresh_metab_list()

    def parseProtocol(self, protocol):
        return 'PRESS' if 'press' in str(protocol or '').lower() else None

    def _build_args(self, params):
        te = params.get('TE')
        te = float(te) if te not in (None, '') else None
        tau1 = params.get('Tau 1')
        tau2 = params.get('Tau 2')
        if tau1 in (None, ''):
            tau1 = te / 2.0 if te is not None else self._num(params, 'Tau 1')
        if tau2 in (None, ''):
            tau2 = te / 2.0 if te is not None else self._num(params, 'Tau 2')
        return ['press_shaped',
                self._num(params, 'Samples'), self._num(params, 'Bandwidth'),
                self._num(params, 'Bfield'), float(params.get('Linewidth') or 1.0),
                float(tau1), float(tau2), self._stage_pulse(params),
                float(params.get('RefTp') or 5.0),
                float(params.get('thkX') or 2.0), float(params.get('thkY') or 2.0),
                float(params.get('fovX') or 3.0), float(params.get('fovY') or 3.0),
                int(float(params.get('nX') or 8)), int(float(params.get('nY') or 8)),
                float(params.get('Flip Angle') or 180.0),
                float(params.get('Sim Centre (ppm)') or 4.65)]


class SpinachSemiLaserShaped(_SpinachRuntime):
    """semi-LASER with one shaped AFP waveform for both refocusing pairs on a
    spatial grid — the fast method (Zhang 2017), Spinach pulses."""

    def __init__(self):
        super().__init__()
        self.name = 'SpinachSemiLaserShaped'
        self.display_name = 'semi-LASER shaped'
        self.file_selection = ['Path to Pulse']
        self.mandatory_params = _shaped_params()
        self._refresh_metab_list()

    def parseProtocol(self, protocol):
        return 'sLASER' if 'laser' in str(protocol or '').lower() else None

    def _build_args(self, params):
        return ['semilaser_shaped',
                self._num(params, 'Samples'), self._num(params, 'Bandwidth'),
                self._num(params, 'Bfield'), float(params.get('Linewidth') or 1.0),
                self._num(params, 'TE'), self._stage_pulse(params),
                float(params.get('RefTp') or 5.0),
                float(params.get('thkX') or 2.0), float(params.get('thkY') or 2.0),
                float(params.get('fovX') or 3.0), float(params.get('fovY') or 3.0),
                int(float(params.get('nX') or 8)), int(float(params.get('nY') or 8)),
                float(params.get('Flip Angle') or 180.0),
                float(params.get('Sim Centre (ppm)') or 4.65)]


SPINACH_BACKENDS = [SpinachBackend, SpinachPressShaped, SpinachSemiLaserShaped]

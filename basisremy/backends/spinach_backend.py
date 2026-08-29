####################################################################################################
#                                        spinach_backend.py                                        #
####################################################################################################
#                                                                                                  #
# Purpose: Spinach backend — ideal-pulse Spin Echo / PRESS / STEAM / LASER simulated by Spinach    #
#          (Kuprov; MIT) through the same Octave runtime FID-A and MRSCloud use. The spin systems  #
#          are FID-A's .mat definitions, so results compare point by point with FidaIdeal.         #
#          Spinach is written for MATLAB R2024b; it runs under Octave with the source patch        #
#          adapters/backends/spinach_octave.patch and the shims in adapters/backends/              #
#          spinach_shims/ (see core/externals.py for the sparse fetch and the patch step).         #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import numpy as np

from basisremy.backends.base import Backend
from basisremy.backends.fida_backends import _DEFAULT_FIDA_METABS


class SpinachBackend(Backend):
    """Spinach through Octave: ideal-pulse sequences on FID-A's spin systems."""

    _KINDS = {'Spin Echo': 'spinecho', 'PRESS': 'press', 'STEAM': 'steam', 'LASER': 'laser'}

    def __init__(self):
        super().__init__()
        self.name = 'Spinach'
        self.display_name = 'Ideal (SE / PRESS / STEAM / LASER)'
        self.category = 'Spinach'
        self.requires_octave = True

        self.metabs = dict(_DEFAULT_FIDA_METABS)
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
        self.optional_params = {'Nucleus': '1H', 'TR': None}
        # TM is only shown for STEAM — rebuild the panel on Sequence changes.
        self.schema_affecting_keys = {'Sequence'}
        self._refresh_metab_list()

    def _refresh_metab_list(self):
        self.mandatory_params['Metabolites'] = [k for k, v in self.metabs.items() if v]

    def get_params_for_mode(self, mode=None):
        params = dict(self.mandatory_params)
        if params.get('Sequence') != 'STEAM':
            params.pop('TM', None)
        return params

    # -------------------------------------------------- sequence mapping
    def map_sequence_in(self, seq: str) -> 'str | None':
        return self.parseProtocol(seq)

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

    def parseREMY(self, MRSinMRS):
        mandatory = {
            'Samples':   MRSinMRS.get('NumberOfDatapoints', None),
            'Bandwidth': MRSinMRS.get('SpectralWidth', None),
            'Bfield':    MRSinMRS.get('B0', None),
            'TE':        MRSinMRS.get('TE', None),
            'Sequence':  self.parseProtocol(MRSinMRS.get('Protocol', None)),
        }
        optional = {
            'Nucleus':         MRSinMRS.get('Nucleus', None),
            'TR':              MRSinMRS.get('TR', None),
            'Model':           MRSinMRS.get('Model', None),
            'SoftwareVersion': MRSinMRS.get('SoftwareVersion', None),
            'BodyPart':        MRSinMRS.get('BodyPart', None),
        }
        return mandatory, optional

    # -------------------------------------------------- Octave
    def setup_octave_paths(self):
        if self.octave is None:
            raise RuntimeError("Octave not initialized.")
        from basisremy.core.externals import ensure
        from basisremy.core.paths import octave_adapters_base
        ensure('spinach')                     # sparse clone + Octave patch (one-time)
        ensure('fidA')                        # the spin-system definitions
        adapters_base = octave_adapters_base(self.octave)
        self.octave.eval("warning('off', 'all');")
        self.octave.addpath(adapters_base + '/backends/')
        self._paths = ('./externals/spinach',
                       adapters_base + '/backends/spinach_shims',
                       './externals/fidA')

    # -------------------------------------------------- arguments
    @staticmethod
    def _num(params, key):
        v = params.get(key)
        if v is None or (isinstance(v, str) and not v.strip()):
            raise ValueError(f"Spinach: '{key}' is required but empty.")
        return float(v)

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

    # -------------------------------------------------- driver
    def run_simulation(self, params, progress_callback=None, stop_event=None):
        if self.octave is None:
            print("Initializing Octave runtime...")
            self.initialize_octave(prefer_docker=True)
        self.setup_octave_paths()

        args = self._build_args(params)
        metabs = params.get('Metabolites') or []
        basis = {}
        for i, metab in enumerate(metabs):
            if stop_event and stop_event.is_set():
                print(f"  ⏹  Stopped before simulating {metab}.")
                break
            fid_re, fid_im, _npts, _sw, _cf = self.octave.feval(
                'spinach_run', metab, *args, *self._paths, nout=5,
            )
            basis[metab] = (np.asarray(fid_re, dtype=float).flatten()
                            + 1j * np.asarray(fid_im, dtype=float).flatten())
            if progress_callback:
                progress_callback(i + 1, len(metabs))
        return basis

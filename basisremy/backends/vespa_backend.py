####################################################################################################
#                                        vespa_backend.py                                          #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 27/08/26                                                                                #
#                                                                                                  #
# Purpose: Vespa (PyGAMMA) backend for basis-set simulation. The parameter interface, REMY         #
#          parsing, and sequence mapping are live; the density-matrix simulation itself requires   #
#          PyGAMMA, which only publishes wheels for Python <= 3.9 (x86_64). The planned runtime    #
#          is a Docker image (analogous to the Octave one) so it works on any Python/platform.     #
#                                                                                                  #
####################################################################################################

import importlib.util

from basisremy.backends.base import Backend


# Metabolites the future PyGAMMA simulation will support (names follow the
# denmatsim spin-system library that ships with the fsl_mrs external, which
# will provide the chemical shifts / J-couplings).
_DEFAULT_METABS = [
    'Ala', 'Asc', 'Asp', 'Cr', 'GABA', 'GPC', 'GSH', 'Glc', 'Gln', 'Glu',
    'Gly', 'Ins', 'Lac', 'NAA', 'NAAG', 'PCh', 'PCr', 'PE', 'Scyllo', 'Tau',
]


class VespaBackend(Backend):
    """Vespa-Simulation (PyGAMMA density-matrix) backend — experimental."""

    def __init__(self):
        super().__init__()
        self.name = 'Vespa'
        self.display_name = 'Vespa (PyGAMMA)'
        self.category = 'Vespa'
        self.requires_octave = False

        self.metabs = {m: True for m in _DEFAULT_METABS}

        self.dropdown = {
            'Sequence': ['PRESS', 'STEAM', 'Spin Echo'],
        }
        # Scan-physics values have NO defaults — they must come from REMY or
        # the user, never masquerade as file metadata.
        self.mandatory_params = {
            'Sequence':    None,
            'Samples':     None,
            'Bandwidth':   None,
            'Bfield':      None,
            'TE':          None,
            'TM':          10,          # STEAM only — hidden otherwise
            'Nucleus':     '1H',
            'Center Freq': None,        # MHz
            'Metabolites': [m for m, v in self.metabs.items() if v],
        }
        self.optional_params = {
            'TR': None,
        }
        # TM is only shown for STEAM — rebuild the panel on Sequence changes.
        self.schema_affecting_keys = {'Sequence'}

    def get_params_for_mode(self, mode=None):
        params = dict(self.mandatory_params)
        if params.get('Sequence') != 'STEAM':
            params.pop('TM', None)
        return params

    # -------------------------------------------------- sequence mapping
    def map_sequence_in(self, seq: str) -> 'str | None':
        if not seq:
            return None
        s = seq.strip().lower()
        for opt in self.dropdown['Sequence']:
            if opt.lower() == s:
                return opt
        if 'mega' in s or 'hermes' in s or 'hercules' in s:
            return None   # edited sequences not supported (yet)
        if 'steam' in s:
            return 'STEAM'
        # 'UnEdited' is the MRSCloud/BigGABA name for a plain acquisition
        if 'press' in s or 'unedited' in s:
            return 'PRESS'
        if 'laser' in s:
            return None   # LASER/sLASER not supported (yet); 'laser' contains 'se'
        if 'spin' in s or 'se' in s:
            return 'Spin Echo'
        return None

    def parseProtocol(self, protocol):
        return self.map_sequence_in(str(protocol)) if protocol else None

    # -------------------------------------------------- REMY
    def parseREMY(self, MRSinMRS):
        def first(*keys, default=None):
            for k in keys:
                if k in MRSinMRS and MRSinMRS[k] not in (None, ''):
                    return MRSinMRS[k]
            return default

        def num(v):
            try:
                return float(v)
            except (TypeError, ValueError):
                return None

        mandatory = {
            'Samples':   num(first('NumberOfDatapoints', 'Samples')),
            'Bandwidth': num(first('SpectralWidth', 'Bandwidth')),
            'Bfield':    num(first('B0', 'Bfield')),
            'TE':        num(first('TE')),
            'Nucleus':   first('Nucleus'),
            'Sequence':  self.parseProtocol(first('Protocol')),
        }
        cf = num(first('Center Freq'))  # MHz everywhere in BasisREMY
        if cf is None and mandatory['Bfield'] is not None:
            cf = 42.577 * mandatory['Bfield']
        mandatory['Center Freq'] = cf

        optional = {
            'TR': num(first('TR')),
        }
        return mandatory, optional

    # -------------------------------------------------- simulation
    @staticmethod
    def pygamma_available() -> bool:
        return importlib.util.find_spec('pygamma') is not None

    def run_simulation(self, params, progress_callback=None, stop_event=None):
        if not self.pygamma_available():
            raise RuntimeError(
                "The Vespa backend needs PyGAMMA, which is not installed.\n"
                "PyGAMMA only publishes wheels for Python <= 3.9 (x86_64):\n"
                "  • a Python 3.9 (x86_64) environment: pip install pygamma\n"
                "A Docker-based PyGAMMA runtime (like the Octave one) is in\n"
                "the works so this will run out of the box on any setup."
            )
        raise NotImplementedError(
            "Vespa (PyGAMMA) simulation is under development — the parameter "
            "interface is ready; the density-matrix simulation lands next."
        )

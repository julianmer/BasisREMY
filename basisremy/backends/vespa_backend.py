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

        # 'PRESS shaped' replaces both ideal 180s with a real refocusing
        # waveform (Vespa's 'PRESS with real 180 pulses'), on-resonance,
        # single voxel position - no spatial grid.
        self.dropdown = {
            'Sequence': ['PRESS', 'PRESS shaped', 'STEAM', 'Spin Echo'],
        }
        self.file_selection = ['Path to Pulse']
        # Scan-physics values have NO defaults — they must come from REMY or
        # the user, never masquerade as file metadata.
        self.mandatory_params = {
            'Sequence':    None,
            'Samples':     None,
            'Bandwidth':   None,
            'Bfield':      None,
            'TE':          None,
            'TM':          10,          # STEAM only — hidden otherwise
            'Path to Pulse': None,      # PRESS shaped only — hidden otherwise
            'RefTp':       5.0,         # PRESS shaped only — pulse duration [ms]
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
        if params.get('Sequence') != 'PRESS shaped':
            params.pop('Path to Pulse', None)
            params.pop('RefTp', None)
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
        from basisremy.core import pygamma_manager
        return pygamma_manager.is_available()

    @staticmethod
    def _load_spin_systems():
        """denmatsim's spin-system library: {name: [sub-system, ...]}."""
        import json
        from basisremy.core.externals import ensure
        from basisremy.core.paths import externals_root
        ensure('fsl_mrs')
        spins_path = (externals_root() / 'fsl_mrs' / 'fsl_mrs'
                      / 'denmatsim' / 'spinSystems.json')
        with open(spins_path, 'r') as f:
            raw = json.load(f)
        def _as_list(v):
            return list(v) if isinstance(v, (list, tuple)) else [v]

        def _as_matrix(v):
            if not isinstance(v, (list, tuple)):
                return [[float(v)]]
            if v and not isinstance(v[0], (list, tuple)):
                return [list(v)]
            return [list(row) for row in v]

        library = {}
        for key, entry in raw.items():
            name = key[3:] if key.startswith('sys') else key
            subs = entry if isinstance(entry, list) else [entry]
            library[name] = [
                {'shifts_ppm': [float(s) for s in _as_list(sub['shifts'])],
                 'j_hz': _as_matrix(sub['j']),
                 'scale': float(sub.get('scaleFactor', 1.0))}
                for sub in subs
            ]
        return library

    def run_simulation(self, params, progress_callback=None, stop_event=None):
        from basisremy.core import pygamma_manager
        import numpy as np

        sequence = params.get('Sequence')
        if sequence not in self.dropdown['Sequence']:
            raise ValueError(
                f"Vespa: unsupported Sequence {sequence!r} — choose one of "
                f"{self.dropdown['Sequence']}.")
        if sequence == 'STEAM' and params.get('TM') in (None, ''):
            raise ValueError("Vespa: STEAM needs a mixing time 'TM' (ms).")
        pulse = None
        if sequence == 'PRESS shaped':
            path = params.get('Path to Pulse')
            if not path:
                raise ValueError("Vespa: 'PRESS shaped' needs a refocusing "
                                 "waveform ('Path to Pulse').")
            if params.get('RefTp') in (None, ''):
                raise ValueError("Vespa: 'PRESS shaped' needs the pulse "
                                 "duration 'RefTp' (ms).")
            from basisremy.core.rf_pulses import load_pulse
            pulse = load_pulse(str(path), float(params['RefTp']), 'ref')
            print(f"  RF pulse '{pulse['name']}': {len(pulse['amp_hz'])} "
                  f"steps, w1max {max(pulse['amp_hz']) / 42.577:.2f} µT")

        runtime = pygamma_manager.preferred_runtime()
        if runtime == 'docker':
            pygamma_manager.ensure_docker_image()
        else:
            pygamma_manager.ensure_env(create=True)
        print(f"Vespa (PyGAMMA) runtime: {runtime}")
        library = self._load_spin_systems()

        cf = float(params['Center Freq'])
        base_job = {
            'sequence': sequence,
            'te_ms': float(params['TE']),
            'samples': int(float(params['Samples'])),
            'bandwidth': float(params['Bandwidth']),
            'cf_mhz': cf,
            'centre_ppm': 4.65,
            'linewidth': float(params.get('Linewidth') or 1.0),
        }
        if sequence == 'STEAM':
            base_job['tm_ms'] = float(params['TM'])
        if pulse is not None:
            base_job['pulse'] = pulse

        metabs = params.get('Metabolites') or []
        self.last_failures = {}   # metab -> reason, surfaced by the GUI
        basis = {}
        for i, metab in enumerate(metabs):
            if stop_event and stop_event.is_set():
                print(f"  ⏹  Stopped before simulating {metab}.")
                break
            if metab not in library:
                print(f"  ⚠️  No denmatsim spin system for '{metab}', skipping")
                self.last_failures[metab] = 'no spin system'
                continue
            job = dict(base_job)
            job['metabolites'] = {metab: library[metab]}
            try:
                try:
                    result = pygamma_manager.run_worker(job, runtime=runtime)
                except RuntimeError as exc:
                    if runtime == 'env' and 'timed out' in str(exc) \
                            and pygamma_manager.docker_available():
                        # wedged side-env (e.g. corrupt Rosetta translation):
                        # switch to Docker now and remember for future runs
                        print("  ⚠️  side-env timed out — switching to the "
                              "Docker runtime")
                        pygamma_manager.prefer_docker(str(exc))
                        runtime = 'docker'
                        pygamma_manager.ensure_docker_image()
                        result = pygamma_manager.run_worker(job, runtime=runtime)
                    else:
                        raise
                entry = result[metab]
                basis[metab] = (np.asarray(entry['re'], dtype=float)
                                + 1j * np.asarray(entry['im'], dtype=float))
                print(f"  ✓ {metab} simulated (PyGAMMA)")
            except Exception as exc:  # noqa: BLE001
                print(f"  ✗ {metab}: {exc}")
                self.last_failures[metab] = str(exc)
            if progress_callback:
                progress_callback(i + 1, len(metabs))
        return basis

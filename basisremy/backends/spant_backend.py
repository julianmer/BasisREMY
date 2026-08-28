####################################################################################################
#                                         spant_backend.py                                         #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 28/08/26                                                                                #
#                                                                                                  #
# Purpose: spant backend — density-matrix simulation with the spant R package (Wilson, MRM 2021)   #
#          through core/spant_manager (Docker image or the user's own R, one Rscript worker).      #
#          Sequences: spant's ideal PRESS / STEAM (three variants) / sLASER / spin echo /          #
#          pulse-acquire, MEGA-PRESS with a Gaussian editing pulse (ON / OFF / DIFF entries)       #
#          and PRESS with a shaped refocusing pulse (frequency-offset sub-simulations, no spatial  #
#          grid). Metabolites come from spant's own spin-system library under spant's names.       #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import os

import numpy as np

from basisremy.backends.base import Backend


class SpantBackend(Backend):

    # spant molecule names (get_mol_names(), 1H entries). Standard brain set on.
    _METABS_DEFAULT = {
        'naa': True, 'naag': True, 'cr': True, 'pcr': True, 'gpc': True,
        'pch': True, 'ins': True, 'sins': True, 'glu': True, 'gln': True,
        'gaba': True, 'gsh': True, 'lac': True, 'tau': True, 'asp': True,
        'ala': True, 'asc': True, 'gly': True, 'glc': True, 'peth': True,
        # available, off by default
        'cho': False, 'h2o': False, 'a_glc': False, 'b_glc': False,
        'bhb': False, 'glyc': False, 'lys': False, '2hg': False, 'cit': False,
        'val': False, 'ace': False, 'pyr': False, 'suc': False, 'thr': False,
        'msm': False,
        # spant's parameterised lipid / macromolecule components
        'lip09': False, 'lip13a': False, 'lip13b': False, 'lip20': False,
        'mm09': False, 'mm12': False, 'mm14': False, 'mm17': False,
        'mm20': False, 'mm_3t': False,
        # alternative parameterisations
        'naa2': False, 'naa_rt': False, 'gaba_rt': False, 'gaba_jn': False,
        'glu_rt': False, 'lac_rt': False, 'ins_rt': False,
        'cr_ch2_rt': False, 'cr_ch3_rt': False, 'naag_ch3': False,
        'm_cr_ch2': False,
    }

    _SEQ_KEY = {
        'PRESS':         'press',
        'STEAM':         'steam',
        'sLASER':        'slaser',
        'Spin Echo':     'spin_echo',
        'MEGA-PRESS':    'mega_press',
        'PRESS shaped':  'press_shaped',
        'Pulse-acquire': 'pulse_acquire',
    }
    _STEAM_VARIANTS = {
        'Standard':           'ideal',
        'Coherence filter':   'cof',    # seq_steam_ideal_cof
        'z-rotation (Young)': 'young',  # seq_steam_ideal_young
    }
    # spant's default timing splits, scaled to the requested TE
    _SLASER_SPLIT_MS = (8.0, 11.0, 9.0)    # seq_slaser_ideal TE1 / TE2 / TE3
    _MEGA_SPLIT_MS = (15.0, 53.0)          # seq_mega_press_ideal TE1 / TE2 (TE 68)

    def __init__(self):
        super().__init__()
        self.name = 'spant'
        self.display_name = 'spant'
        self.category = 'spant'
        self.requires_octave = False
        self.metabs = dict(self._METABS_DEFAULT)
        self.dropdown = {
            'Sequence':      list(self._SEQ_KEY),
            'STEAM Variant': list(self._STEAM_VARIANTS),
        }
        self.file_selection = ['Path to Pulse']
        # Scan-physics values have NO defaults — they come from REMY or the user.
        self.mandatory_params = {
            'Sequence':      None,
            'Samples':       None,
            'Bandwidth':     None,
            'Bfield':        None,
            'Center Freq':   None,        # MHz (γ·B0 when blank)
            'TE':            None,
            'Tau 1':         None,        # PRESS: TE1, blank → TE/2
            'Tau 2':         None,        # PRESS: TE2, blank → TE/2
            'TM':            10.0,        # STEAM only
            'STEAM Variant': 'Standard',  # STEAM only
            'Edit On':       1.9,         # MEGA-PRESS only (ppm)
            'Edit Off':      7.5,
            'Edit Bandwidth (Hz)': 110.0,
            'Path to Pulse': None,        # PRESS shaped only
            'RefTp':         5.0,         # PRESS shaped only [ms]
            'Flip Angle':    180.0,       # PRESS shaped only
            'Linewidth':     1.0,
            'Metabolites':   [m for m, v in self.metabs.items() if v],
        }
        self.optional_params = {'Nucleus': '1H', 'TR': None}
        self.schema_affecting_keys = {'Sequence'}

    # ------------------------------------------------------------ schema
    _PRESS_KEYS = ('Tau 1', 'Tau 2')
    _STEAM_KEYS = ('TM', 'STEAM Variant')
    _MEGA_KEYS = ('Edit On', 'Edit Off', 'Edit Bandwidth (Hz)')
    _SHAPED_KEYS = ('Path to Pulse', 'RefTp', 'Flip Angle')

    def get_params_for_mode(self, mode=None):
        params = dict(self.mandatory_params)
        seq = params.get('Sequence')
        show = set()
        if seq in ('PRESS', 'PRESS shaped'):
            show |= set(self._PRESS_KEYS)
        if seq == 'PRESS shaped':
            show |= set(self._SHAPED_KEYS)
        if seq == 'STEAM':
            show |= set(self._STEAM_KEYS)
        if seq == 'MEGA-PRESS':
            show |= set(self._MEGA_KEYS)
        for k in (self._PRESS_KEYS + self._STEAM_KEYS + self._MEGA_KEYS
                  + self._SHAPED_KEYS):
            if k not in show:
                params.pop(k, None)
        if seq == 'Pulse-acquire':
            params.pop('TE', None)
        return params

    # ------------------------------------------------------------ sequence mapping
    def map_sequence_in(self, seq: str) -> 'str | None':
        if not seq:
            return None
        s = seq.strip().lower()
        for opt in self.dropdown['Sequence']:
            if opt.lower() == s:
                return opt
        if 'hermes' in s or 'hercules' in s:
            return None          # spant has no HERMES / HERCULES sequence
        if 'mega' in s:
            return 'MEGA-PRESS'
        if 'steam' in s:
            return 'STEAM'
        if 'slaser' in s or 'semi' in s:
            return 'sLASER'
        if 'laser' in s:
            return None          # no full LASER in spant; 'laser' contains 'se'
        # 'UnEdited' is the MRSCloud / BigGABA name for a plain acquisition
        if 'press' in s or 'unedited' in s:
            return 'PRESS'
        if 'spin' in s or 'se' in s:
            return 'Spin Echo'
        return None

    def parseProtocol(self, protocol):
        return self.map_sequence_in(str(protocol)) if protocol else None

    # ------------------------------------------------------------ REMY
    def parseREMY(self, MRSinMRS):
        def first(*keys):
            for k in keys:
                if k in MRSinMRS and MRSinMRS[k] not in (None, ''):
                    return MRSinMRS[k]
            return None

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
            'Sequence':  self.parseProtocol(first('Protocol')),
        }
        cf = num(first('Center Freq'))   # MHz everywhere in BasisREMY
        if cf is None and mandatory['Bfield'] is not None:
            cf = 42.577 * mandatory['Bfield']
        mandatory['Center Freq'] = cf
        optional = {'TR': num(first('TR')), 'Nucleus': first('Nucleus')}
        return mandatory, optional

    # ------------------------------------------------------------ job
    @staticmethod
    def _num(params, key, what=None):
        try:
            return float(params.get(key))
        except (TypeError, ValueError):
            raise ValueError(
                f"spant: '{key}' ({what or key}) must be set before simulating.") from None

    @staticmethod
    def _blank(v):
        return v is None or (isinstance(v, str) and not v.strip())

    def _build_job(self, params, metab):
        """The spant_worker.R job for one metabolite (validates the inputs)."""
        seq = params.get('Sequence')
        key = self._SEQ_KEY.get(seq)
        if key is None:
            raise ValueError(
                f"spant: Sequence must be one of {list(self._SEQ_KEY)} (got {seq!r}).")
        samples = int(self._num(params, 'Samples', 'points'))
        bandwidth = self._num(params, 'Bandwidth', 'spectral width in Hz')
        if not self._blank(params.get('Center Freq')):
            ft_hz = self._num(params, 'Center Freq', 'MHz') * 1e6
        else:
            ft_hz = 42.577e6 * self._num(params, 'Bfield', 'field strength in T')
        job = {
            'sequence': key,
            'ft_hz': ft_hz,
            'fs_hz': bandwidth,
            'n': samples,
            'ref_ppm': 4.65,
            'linewidth_hz': float(params.get('Linewidth') or 1.0),
            'metabolites': [metab],
        }
        if key == 'pulse_acquire':
            return job
        te = self._num(params, 'TE', 'echo time in ms')
        if key in ('press', 'press_shaped'):
            tau1 = params.get('Tau 1')
            tau2 = params.get('Tau 2')
            tau1 = te / 2.0 if self._blank(tau1) else float(tau1)
            tau2 = te / 2.0 if self._blank(tau2) else float(tau2)
            job.update(te1_s=tau1 / 1e3, te2_s=tau2 / 1e3)
        if key == 'press_shaped':
            pulse = params.get('Path to Pulse')
            if self._blank(pulse) or not os.path.isfile(str(pulse)):
                raise ValueError(
                    "spant: 'Path to Pulse' (refocusing waveform: .pta, Bruker "
                    "or two-column ASCII) is required for PRESS shaped.")
            ext = os.path.splitext(str(pulse))[1].lower()
            fmt = 'pta' if ext == '.pta' else 'bruker' if ext in ('.exc', '.inv', '.rfc') else 'ascii'
            job.update(pulse_file=os.path.abspath(str(pulse)), pulse_format=fmt,
                       pulse_dur_s=self._num(params, 'RefTp', 'pulse duration in ms') / 1e3,
                       refoc_flip_deg=float(params.get('Flip Angle') or 180.0))
        elif key == 'steam':
            job.update(te_s=te / 1e3,
                       tm_s=self._num(params, 'TM', 'mixing time in ms') / 1e3,
                       steam_variant=self._STEAM_VARIANTS.get(
                           params.get('STEAM Variant') or 'Standard', 'ideal'))
        elif key == 'slaser':
            total = sum(self._SLASER_SPLIT_MS)
            te1, te2, te3 = (te * s / total for s in self._SLASER_SPLIT_MS)
            job.update(te1_s=te1 / 1e3, te2_s=te2 / 1e3, te3_s=te3 / 1e3)
        elif key == 'spin_echo':
            job.update(te_s=te / 1e3)
        elif key == 'mega_press':
            total = sum(self._MEGA_SPLIT_MS)
            te1, te2 = (te * s / total for s in self._MEGA_SPLIT_MS)
            job.update(te1_s=te1 / 1e3, te2_s=te2 / 1e3,
                       edit_on_ppm=self._num(params, 'Edit On', 'ppm'),
                       edit_off_ppm=self._num(params, 'Edit Off', 'ppm'),
                       edit_bw_hz=self._num(params, 'Edit Bandwidth (Hz)', 'Hz'))
        return job

    # ------------------------------------------------------------ simulation
    @staticmethod
    def _fid(entry):
        # spant's signals rotate in the opposite sense to FID-A's (and to the
        # fft + "4.65 + f/f0" axis the GUI and exporters assume): taken as-is,
        # NAA showed up mirrored at 7.3 ppm. The complex conjugate puts every
        # resonance at its own shift.
        return np.conj(np.asarray(entry['re'], dtype=float)
                       + 1j * np.asarray(entry['im'], dtype=float))

    def run_simulation(self, params, progress_callback=None, stop_event=None):
        from basisremy.core import spant_manager
        metabs = list(params.get('Metabolites') or [])
        if not metabs:
            raise ValueError("spant: no metabolites selected.")
        self._build_job(params, metabs[0])   # validate the inputs before any work

        runtime = spant_manager.preferred_runtime()
        where = spant_manager.ensure_runtime(runtime)
        print(f"spant runtime: {where}")

        self.last_failures = {}   # metab -> reason, surfaced by the GUI
        basis = {}
        for i, metab in enumerate(metabs):
            if stop_event and stop_event.is_set():
                print(f"  ⏹  Stopped before simulating {metab}.")
                break
            try:
                job = self._build_job(params, metab)
                result = spant_manager.run_worker(job, runtime=runtime)
                entry = result['basis'][metab]
                if job['sequence'] == 'mega_press':
                    on, off = self._fid(entry['on']), self._fid(entry['off'])
                    basis[f'{metab} (ON)'] = on
                    basis[f'{metab} (OFF)'] = off
                    basis[f'{metab} (DIFF)'] = on - off
                else:
                    basis[metab] = self._fid(entry)
                print(f"  ✓ {metab} simulated (spant {result.get('spant_version', '')})")
            except Exception as exc:  # noqa: BLE001
                print(f"  ✗ {metab}: {exc}")
                self.last_failures[metab] = str(exc)
            if progress_callback:
                progress_callback(i + 1, len(metabs))
        return basis

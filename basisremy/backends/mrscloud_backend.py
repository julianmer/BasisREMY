####################################################################################################
#                                        mrscloud_backend.py                                       #
####################################################################################################
#                                                                                                  #
# Authors: G. Simegn (gsimegn1@jh.edu)                                                             #
#          J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 18/02/26                                                                                #
# Rewritten: 24/04/26                                                                              #
#                                                                                                  #
# Purpose: Per-metabolite execution goes through `adapters/backends/mrscloud_run_metab.m`, a thin  #
#          Octave adapter that mirrors `externals/mrscloud/run/run_simulations_cloud.m` and        #
#          returns the FID as plain numeric arrays so oct2py can ferry them across.                #
#                                                                                                  #
# Notes                                                                                            #
#   - Output Path / Make .raw / Add Ref. are no longer user-facing; final export is handled by the #
#     post-simulation Export dialog (core/exporters).                                              #
#   - The metabolite list matches the official MRSCloud README. Some entries (Cystat, HCar, iLe,   #
#     Lys, Glc) simulate slowly — flagged with TODO.                                               #
#   - MEGA returns '<metab> (ON)' / '(OFF)' / '(DIFF)' entries (DIFF = ON − OFF, as the FID-A      #
#     backends do). HERMES / HERCULES return the four sub-experiments '(A)' … '(D)' (editing       #
#     pulses at the scheme's offsets, in MRSCloud's order) and their '(SUM)'; the GABA / GSH       #
#     difference combinations are a TODO (see run_simulation).                                     #
#                                                                                                  #
####################################################################################################


#*************#
#   imports   #
#*************#
from __future__ import annotations

import os

import numpy as np

# own
from basisremy.backends.base import Backend


#**************************************************************************************************#
#                                         MRSCloudBackend                                          #
#**************************************************************************************************#
class MRSCloudBackend(Backend):

    # Canonical MRSCloud metabolite list (externals/mrscloud/README.md).
    # The bool indicates whether the metabolite is enabled by default in the GUI.
    _METABS_DEFAULT = {
        # Common healthy-brain metabolites (always available)
        'Asc':  True,  'Asp': True,  'Cr':  True,  'EA':   False,
        'GABA': True,  'GPC': True,  'GSH': True,  'Gln':  True,
        'Glu':  True,  'Gly': True,  'H2O': False, 'Lac':  True,
        'mI':   True,  'NAA': True,  'NAAG':True,  'PCh':  True,
        'PCr':  True,  'PE':  True,  'Ser': False, 'sI':   True,
        'Tau':  True,
        # Specific-interest (slower / niche)
        'Ala':  False, 'Ace':  False, 'AcO':  False, 'AcAc': False,
        'Cit':  False, 'Cystat': False,  # TODO slow spin system
        'HCar': False,                    # TODO slow spin system
        'Lys':  False,                    # TODO slow spin system
        'Thr':  False, 'bHG':  False, 'Tyros': False, 'Val': False,
        'Phenyl': False, 'bHB': False, 'Gua': False,
        'iLe':  False,                    # TODO slow spin system
        'Pyr':  False, 'Suc':  False, 'Tryp': False,
        # Exogenous compounds
        'EtOH': False, 'MSM':  False,
    }

    _SEQUENCES     = ['UnEdited', 'MEGA', 'HERMES', 'HERCULES']
    _LOCALIZATIONS = ['PRESS', 'sLASER', 'STEAM_7T']
    # Plain scanner vendors shown in the System dropdown. Universal-vs-vendor-
    # specific is NOT a vendor — it is the Mode (see self.modes). The stored
    # System value is always one of these; the Universal_* translation happens
    # behind the scenes at run time (see _mrscloud_vendor / run_simulation).
    _SYSTEMS       = ['Philips', 'Siemens', 'GE']
    _VENDORS       = ['Philips', 'Universal_Philips', 'Siemens', 'Universal_Siemens', 'GE']

    # ---- Pulse-file requirements --------------------------------------
    # The waveform names `externals/mrscloud/functions/load_parameters.m`
    # requests per vendor / localization / editing scheme. MRSCloud's README
    # (Remark 3): "Product sequence and rf waveform are not shared in the
    # GitHub repo" — only the pulses_universal/ set ships. Every other name
    # is vendor-confidential: the GUI asks for it with a file picker, or the
    # user drops it into externals/mrscloud/pulses/.
    #
    # Excitation: load_parameters hard-codes Philips_spredrex.pta for every
    # vendor; BasisREMY substitutes the bundled universal excitation for it
    # (see _stage_pulse_shims), so it is never listed as required.
    _REFOC_PULSES = {          # per vendor: (PRESS, sLASER) refocusing waveform
        'Philips':           ('gtst1203_sp.pta',          'Philips_GOIA_WURST_100pts.mat'),
        'Siemens':           ('orig_refoc_mao_100_4.pta', 'Philips_GOIA_WURST_100pts.mat'),
        'GE':                ('GE_rfa_3.9ms.pta',         'GE_GOIA_WURST_100pts.mat'),
        'Universal_Philips': ('univ_eddenrefo.pta',       'Philips_GOIA_WURST_100pts.mat'),
        'Universal_Siemens': ('univ_eddenrefo.pta',       'Philips_GOIA_WURST_100pts.mat'),
    }
    _EDIT_PULSES = {           # per vendor and editing scheme
        'Philips': {
            'MEGA':     ['sg100_100_0_14ms_88hz.pta'],
            'HERMES':   ['sg100_100_0_14ms_88hz.pta', 'dl_Philips_4_56_1_90.pta'],
            'HERCULES': ['sg100_100_0_14ms_88hz.pta', 'dl_Philips_4_58_1_90.pta',
                         'dl_Philips_4_18_1_90.pta'],
        },
        'Siemens': {
            'MEGA':     ['Siemens_filtered_editing.pta'],
            'HERMES':   ['Siemens_filtered_editing.pta', 'dl_Siemens_4_56_1_90.pta'],
            'HERCULES': ['Siemens_filtered_editing.pta', 'dl_Siemens_4_58_1_90.pta',
                         'dl_Siemens_4_18_1_90.pta'],
        },
        'GE': {                # upstream reuses the Philips editing waveforms
            'MEGA':     ['sg100_100_0_14ms_88hz.pta'],
            'HERMES':   ['sg100_100_0_14ms_88hz.pta', 'dl_Philips_4_56_1_90.pta'],
            'HERCULES': ['sg100_100_0_14ms_88hz.pta', 'dl_Philips_4_58_1_90.pta',
                         'dl_Philips_4_18_1_90.pta'],
        },
        'Universal_Philips': {
            'MEGA':     ['sl_univ_pulse.pta'],
            'HERMES':   ['sl_univ_pulse.pta', 'dl_Philips_univ_4_56_1_90.pta'],
            'HERCULES': ['sl_univ_pulse.pta', 'dl_Philips_4_58_1_90.pta',
                         'dl_Philips_4_18_1_90.pta'],
        },
        'Universal_Siemens': {
            'MEGA':     ['sl_univ_pulse.pta'],
            'HERMES':   ['sl_univ_pulse.pta', 'dl_Siemens_4_56_1_90.pta'],
            'HERCULES': ['sl_univ_pulse.pta', 'dl_Siemens_4_58_1_90.pta',
                         'dl_Siemens_4_18_1_90.pta'],
        },
    }
    # Shipped in externals/mrscloud/pulses_universal/.
    _BUNDLED_PULSES = {
        'sl_univ_pulse.pta', 'univ_eddenrefo.pta', 'univ_spreddenrex.pta',
        'dl_Univ_4_68_1_9_20ms.pta',
        'dl_Philips_univ_3_67_1_9_20ms.pta', 'dl_Philips_univ_3_67_4_56_20ms.pta',
        'dl_Philips_univ_4_56_1_9_20ms.pta', 'dl_Philips_univ_4_68_1_9_20ms.pta',
        'dl_Siemens_univ_3_67_1_9_20ms.pta', 'dl_Siemens_univ_3_67_4_56_20ms.pta',
        'dl_Siemens_univ_4_56_1_9_20ms.pta', 'dl_Siemens_univ_4_68_1_9_20ms.pta',
    }
    # Names load_parameters.m requests that ship under a later name: the
    # universal 4.56 / 1.9 ppm dual-lobe pulse was committed upstream on
    # 2023-03-30 as *_1_9_20ms.pta while the HERMES branches still ask for
    # *_1_90.pta, a file the public repo never had. Staged under the
    # requested name by _stage_pulse_shims.
    _PULSE_ALIASES = {
        'dl_Philips_univ_4_56_1_90.pta': 'dl_Philips_univ_4_56_1_9_20ms.pta',
        'dl_Siemens_4_56_1_90.pta':      'dl_Siemens_univ_4_56_1_9_20ms.pta',
    }

    @classmethod
    def required_pulse_files(cls, vendor: str, sequence: str, localization: str) -> list[str]:
        """Vendor-confidential pulse files MRSCloud will demand for this
        (vendor, seq, loc) combo, as paths relative to externals/mrscloud/.
        Bundled and aliased universal waveforms are not listed."""
        if vendor not in cls._REFOC_PULSES:
            return []
        names = []
        if localization == 'PRESS':
            names.append(cls._REFOC_PULSES[vendor][0])
        elif localization == 'sLASER':
            names.append(cls._REFOC_PULSES[vendor][1])
        names += cls._EDIT_PULSES[vendor].get(sequence, [])
        out = []
        for name in names:
            rel = f'pulses/{name}'
            if (name in cls._BUNDLED_PULSES or name in cls._PULSE_ALIASES
                    or rel in out):
                continue
            out.append(rel)
        return out

    @classmethod
    def missing_pulse_files(cls, vendor: str, sequence: str, localization: str,
                            mrscloud_root: str = './externals/mrscloud') -> list[str]:
        """Subset of required pulse files that are NOT present on disk."""
        return [rel for rel in cls.required_pulse_files(vendor, sequence, localization)
                if not os.path.exists(os.path.join(mrscloud_root, rel))]

    def __init__(self):
        super().__init__()
        self.name = 'MRSCloud'
        self.display_name = 'MRSCloud'
        self.category = 'MRSCloud'
        self.requires_octave = True

        # Universal (bundled, vendor-agnostic) vs vendor-specific System
        # options. Exposed in the GUI as a "Mode" selector; the choice filters
        # which vendors appear in self.dropdown['System'] (see
        # get_params_for_mode / set_mode below).
        self.modes = ['Universal', 'Non-Universal']
        self.current_mode = 'Universal'

        # Metabolite library
        self.metabs = dict(self._METABS_DEFAULT)

        # Dropdowns shown in the GUI. NOTE: the GUI rebuilds the parameter
        # panel whenever a key in `schema_affecting_keys` changes, so the
        # editing fields appear/disappear depending on `Sequence`.
        self.dropdown = {
            'System':         list(self._SYSTEMS),
            'Sequence':       list(self._SEQUENCES),
            'Localization':   list(self._LOCALIZATIONS),
        }

        # GUI tells us when these change so we can rebuild the visible
        # parameter list (e.g. show editing fields only for MEGA/HERMES/HERCULES,
        # show the pulse-file picker only when the pulse is missing).
        self.schema_affecting_keys = {'Sequence', 'Localization', 'System'}

        # Pulse-file pickers are populated dynamically by get_params_for_mode().
        self.file_selection: list[str] = []

        # Mandatory parameters — only the ones MRSCloud actually consumes.
        # Not exposed:
        #   * Center Freq   — computed from Bfield × γ inside FID-A
        #   * Linewidth     — MRSCloud hard-codes lw = 1 Hz in load_parameters
        #   * TE2           — MRSCloud uses a single TE; TE1 is set per-vendor
        # Bfield is the acquisition's own field strength (REMY or the user).
        # load_parameters.m hard-codes B0 per vendor (2.89 / 3 T; 7 T only for
        # STEAM_7T), so the adapter rebuilds the field-dependent parts at this
        # value; below 2.25 T MRSCloud's 1.5 T parameter set is used. Editing
        # fields (Edit On / Off / Tp) live in `_edit_params` and are spliced
        # in by get_params_for_mode() for MEGA only.
        self.mandatory_params = {
            'System':         None,        # vendor (must be selected)
            'Sequence':       None,        # UnEdited / MEGA / HERMES / HERCULES
            'Localization':   None,        # PRESS / sLASER / STEAM_7T
            'Bfield':         None,        # T — from REMY (B0) or the user
            'Samples':        None,
            'Bandwidth':      None,
            'TE':             None,
            'Spatial Points': 41,          # 41 acceptable, 101 ideal (slow)
            'Metabolites':    [k for k, v in self.metabs.items() if v],
        }

        # Editing parameters, MEGA only (added to mandatory_params on the
        # fly). HERMES / HERCULES are fixed schemes in MRSCloud — the offsets
        # (adapter) and the 20 ms editing pulses (load_parameters.m) cannot be
        # changed — so they expose no editing fields.
        self._edit_params = {
            'Edit On':        1.9,         # ppm
            'Edit Off':       7.5,         # ppm
            'Edit Tp':        14,          # ms
        }

        # Vendor pulse file the user must supply when not bundled.
        # Key is the GUI label; value is the canonical filename MRSCloud
        # expects (so we can rename-on-copy into the workdir).
        self._pulse_param_label = 'Vendor Pulse File'

        # Optional / REMY-extracted only — NOT shown in the GUI panel.
        self.optional_params = {
            'Nucleus': None,
            'TR':      None,
        }

    # --------------------------------------------------------------- mode/schema
    # Vendor equivalences across modes. A vendor with no counterpart in the
    # target mode (e.g. GE has no Universal waveform) maps to None and clears.
    _VENDOR_TO_UNIVERSAL = {
        'Philips':           'Universal_Philips',
        'Siemens':           'Universal_Siemens',
        'Universal_Philips': 'Universal_Philips',
        'Universal_Siemens': 'Universal_Siemens',
        'GE':                None,
    }
    _VENDOR_TO_SPECIFIC = {
        'Universal_Philips': 'Philips',
        'Universal_Siemens': 'Siemens',
        'Philips':           'Philips',
        'Siemens':           'Siemens',
        'GE':                'GE',
    }

    def _mrscloud_vendor(self, vendor=None, mode=None):
        """Translate the stored plain vendor + Mode into the vendor label that
        MRSCloud's load_parameters.m expects ('Philips' / 'Universal_Philips' /
        ...). Universal mode maps Philips/Siemens onto their Universal_* bundled
        waveforms; Non-Universal (and GE, which has no Universal set) keeps the
        plain vendor. This is the only place the Universal_* name is produced."""
        vendor = vendor if vendor is not None else self.mandatory_params.get('System')
        mode = mode or self.current_mode
        if not vendor:
            return vendor
        # Normalize any accidental Universal_* value back to the plain vendor.
        plain = self._VENDOR_TO_SPECIFIC.get(vendor, vendor)
        if mode == 'Universal':
            return self._VENDOR_TO_UNIVERSAL.get(plain) or plain
        return plain

    def set_mode(self, mode):
        """Switch between 'Universal' and 'Non-Universal'. The stored System
        (plain scanner vendor) is left untouched — the Universal_* translation
        happens behind the scenes at run time. All other fields are preserved."""
        if mode not in self.modes:
            raise ValueError(f"Unknown mode '{mode}'. Available: {self.modes}")
        self.current_mode = mode
        return self.get_params_for_mode(mode)

    def get_params_for_mode(self, mode=None):
        """Return only the parameters relevant to the current Sequence /
        Localization / System combo so the GUI never shows fields MRSCloud
        will silently ignore."""
        params = dict(self.mandatory_params)
        seq    = (self.mandatory_params.get('Sequence')       or '').strip()
        loc    = (self.mandatory_params.get('Localization')   or '').strip()
        vendor = (self.mandatory_params.get('System')         or '').strip()

        # ---- System (plain vendor) + Mode (Universal / Non-Universal) -------
        # System always stores the plain scanner vendor; the Mode decides
        # whether the bundled Universal_* waveforms are used. Normalize any
        # legacy Universal_* value back to the plain vendor (Mode now carries
        # universality) so switching modes never loses the stored System.
        if vendor in ('Universal_Philips', 'Universal_Siemens'):
            vendor = self._VENDOR_TO_SPECIFIC[vendor]
            self.mandatory_params['System'] = vendor
            self.current_mode = 'Universal'
        # GE has no Universal waveform set — force Non-Universal when selected.
        if vendor == 'GE':
            self.current_mode = 'Non-Universal'
        # The System dropdown always lists the plain scanner vendors.
        self.dropdown['System'] = list(self._SYSTEMS)

        # ---- edited sequences restrict Localization to PRESS / sLASER -------
        if seq in ('MEGA', 'HERMES', 'HERCULES'):
            restricted = ['PRESS', 'sLASER']
            self.dropdown['Localization'] = restricted
            # If a STEAM_7T localization carried over, reset to PRESS
            if loc not in restricted:
                self.mandatory_params['Localization'] = 'PRESS'
                params['Localization'] = 'PRESS'
        else:
            # UnEdited — restore all localization options
            self.dropdown['Localization'] = list(self._LOCALIZATIONS)

        # ---- splice editing fields only when relevant -----------------------
        if seq == 'MEGA':
            for k, v in self._edit_params.items():
                params.setdefault(k, v)
                self.mandatory_params.setdefault(k, v)
        else:
            # UnEdited, and HERMES / HERCULES whose editing scheme MRSCloud
            # fixes → hide everything edit-related
            for k in self._edit_params:
                params.pop(k, None)
                self.mandatory_params.pop(k, None)

        # ---- ask for the missing vendor pulse file when needed ---------------
        # Universal mode covers PRESS with bundled waveforms, but sLASER
        # (GOIA-WURST) and HERCULES (vendor dual-lobe pulses) need
        # vendor-confidential files in either mode.
        cur_loc = params.get('Localization') or loc
        self.file_selection = []
        missing = []
        if seq and vendor and cur_loc:
            missing = self.missing_pulse_files(
                self._mrscloud_vendor(vendor), seq, cur_loc)
        if missing:
            label = self._pulse_param_label
            self.file_selection.append(label)
            params[label] = self.mandatory_params.get(label)
            self.mandatory_params.setdefault(label, None)
        else:
            # everything MRSCloud will load is bundled or on disk — no picker
            params.pop(self._pulse_param_label, None)
            self.mandatory_params.pop(self._pulse_param_label, None)

        # Re-order Metabolites to the bottom for the GUI grid
        if 'Metabolites' in params:
            mets = params.pop('Metabolites')
            params['Metabolites'] = mets
        return params

    # ----------------------------------------------------------------- REMY parsing
    def parseREMY(self, MRSinMRS):
        """Map REMY-extracted metadata onto MRSCloud parameters.

        Note: `Center Freq` is not returned — MRSCloud derives it from Bfield.
        """
        bfield = MRSinMRS.get('B0', None)
        protocol = MRSinMRS.get('Protocol', None)

        # Vendor is parsed plainly (Philips / Siemens / GE). Universal-vs-
        # vendor-specific is a Mode, not a vendor: default to Universal when the
        # vendor has bundled universal waveforms; GE has none, so it stays
        # Non-Universal. The stored System keeps the plain vendor either way.
        vendor = self.parseSystem(MRSinMRS.get('Manufacturer', None))
        plain  = self._VENDOR_TO_SPECIFIC.get(vendor, vendor)
        if plain in ('Philips', 'Siemens'):
            self.current_mode = 'Universal'
        elif plain == 'GE':
            self.current_mode = 'Non-Universal'

        mandatory = {
            'Sequence':       self.parseProtocol(protocol),
            'Localization':   self.parseLocalization(protocol),
            'System':         plain,
            'Bfield':         bfield,
            'Samples':        MRSinMRS.get('NumberOfDatapoints', None),
            'Bandwidth':      MRSinMRS.get('SpectralWidth', None),
            'TE':             MRSinMRS.get('TE', None),
        }
        optional = {
            'Nucleus':         MRSinMRS.get('Nucleus', None),
            'TR':              MRSinMRS.get('TR', None),
            'Model':           MRSinMRS.get('Model', None),
            'SoftwareVersion': MRSinMRS.get('SoftwareVersion', None),
            'BodyPart':        MRSinMRS.get('BodyPart', None),
            'VOI':             MRSinMRS.get('VOI', None),
            'AnteriorPosteriorSize': MRSinMRS.get('AnteriorPosteriorSize', None),
            'LeftRightSize':         MRSinMRS.get('LeftRightSize', None),
            'CranioCaudalSize':      MRSinMRS.get('CranioCaudalSize', None),
            'NumberOfAverages':      MRSinMRS.get('NumberOfAverages', None),
            'WaterSuppression':      MRSinMRS.get('WaterSuppression', None),
        }
        return mandatory, optional

    @staticmethod
    def _parameter_set(b0: float) -> str:
        """MRSCloud parameter set (pulses, timings) for a field strength:
        load_parameters_1_5T.m below 2.25 T, load_parameters.m otherwise.
        The field itself is passed separately and applied by the adapter."""
        return '1.5T' if float(b0) < 2.25 else '3T'

    def parseProtocol(self, protocol):
        """Return the MRSCloud editing-scheme label from a raw protocol string.

        Returns one of: 'UnEdited', 'MEGA', 'HERMES', 'HERCULES', or None.
        Use parseLocalization() separately to get the localisation.
        """
        if protocol is None:
            return None
        p = str(protocol).lower()
        if 'hercules' in p:
            return 'HERCULES'
        if 'hermes' in p:
            return 'HERMES'
        if 'mega' in p:
            return 'MEGA'
        return 'UnEdited'

    def parseLocalization(self, protocol):
        """Return the MRSCloud localisation label from a raw protocol string.

        Returns one of: 'PRESS', 'sLASER', 'STEAM_7T', or None.
        """
        if protocol is None:
            return None
        p = str(protocol).lower()
        if 'steam' in p:
            return 'STEAM_7T'
        if ('slaser' in p or 'semi_laser' in p
                or 'semi-laser' in p or 'semilaser' in p):
            return 'sLASER'
        return 'PRESS'

    def parseSystem(self, system):
        """Map a vendor string onto the MRSCloud vendor list.

        Returns the canonical vendor name as stored in self.dropdown['System'].
        GE, Universal_Philips, and Universal_Siemens are passed through directly;
        Philips and Siemens are returned as-is so the GUI accurately reflects the
        scanner that acquired the data (the user can switch to a Universal_*
        variant manually if the vendor pulse files are not available).
        """
        if system is None:
            return None
        s = str(system).lower()
        # Explicit universal tags take priority
        if 'philips_universal' in s or 'universal_philips' in s:
            return 'Universal_Philips'
        if 'siemens_universal' in s or 'universal_siemens' in s:
            return 'Universal_Siemens'
        if 'philips' in s:
            return 'Philips'
        if 'siemens' in s:
            return 'Siemens'
        if 'ge' in s:
            return 'GE'
        print(f"Warning: MRSCloud unsupported vendor '{system}'.")
        return None

    def map_sequence_in(self, seq: str) -> 'str | None':
        """Translate any sequence name into MRSCloud's Sequence vocabulary."""
        if not seq:
            return None
        s = seq.strip().lower()
        # Exact match first (handles same-label backends)
        for opt in self._SEQUENCES:
            if opt.lower() == s:
                return opt
        # Cross-backend synonyms
        if 'hercules' in s:
            return 'HERCULES'
        if 'hermes' in s:
            return 'HERMES'
        if 'mega' in s:
            return 'MEGA'
        return 'UnEdited'

    # --------------------------------------------------------------- Octave paths
    def setup_octave_paths(self, octave=None):
        octave = octave or self.octave
        if octave is None:
            raise RuntimeError("Octave not initialized. Call initialize_octave() first.")
        # Fetch MRSCloud on first use (no-op in a source checkout).
        from basisremy.core.externals import ensure
        from basisremy.core.paths import octave_adapters_base
        ensure('mrscloud')
        octave.eval("warning('off', 'all');")
        # adapters/backends contains our mrscloud_run_metab.m wrapper
        octave.addpath(octave_adapters_base(octave) + '/backends/')
        # Pull in MRSCloud (functions + bundled FID-A) recursively
        octave.addpath(octave.genpath('./externals/mrscloud/functions/'))
        octave.addpath(octave.genpath('./externals/mrscloud/pulses_universal/'))
        # vendor-confidential waveforms the user placed next to the bundled ones
        if os.path.isdir('./externals/mrscloud/pulses'):
            octave.addpath('./externals/mrscloud/pulses/')
        # the run's pulse shims (extra worker sessions are set up after staging)
        if getattr(self, '_shim_path', None):
            octave.addpath(self._shim_path)

    # ----------------------------------------------------- pulse-file shimming
    def _stage_user_pulse(self, workdir: str, vendor: str, sequence: str,
                          localization: str, user_path: str | None) -> None:
        """Copy a user-supplied vendor pulse file into the workdir.

        Triggered only when `Vendor Pulse File` is exposed in the GUI
        (i.e. for non-Universal_* vendors that need a confidential pulse
        which isn't shipped with the public MRSCloud repo). The file is
        renamed to the *first* missing canonical filename so MRSCloud can
        find it via `io_loadRFwaveform`.
        """
        import os, shutil
        if not user_path:
            return
        missing = self.missing_pulse_files(vendor, sequence, localization)
        if not missing:
            return
        wanted = [os.path.basename(m) for m in missing]
        picked = os.path.basename(user_path)
        # a file named like one of the missing waveforms keeps its name;
        # anything else is taken as the first missing one
        target_name = picked if picked in wanted else wanted[0]
        dst = os.path.join(workdir, target_name)
        try:
            if os.path.abspath(user_path) != os.path.abspath(dst):
                shutil.copyfile(user_path, dst)
            print(f"  ✓ Staged user pulse '{os.path.basename(user_path)}' "
                  f"as '{target_name}' in workdir")
        except Exception as e:
            print(f"  ⚠️  Could not stage user pulse {user_path}: {e}")

    def _stage_pulse_shims(self, workdir: str) -> None:
        """Stage the bundled universal waveforms under the names MRSCloud asks for.

        Excitation: work around the hard-coded `Philips_spredrex.pta` pulse.

        `externals/mrscloud/functions/load_parameters.m` (line ~485) calls

            excWaveform = 'Philips_spredrex.pta';
            io_loadRFwaveform(excWaveform, 'exc', 0)

        for *every* vendor — the "universal" alternative `univ_spreddenrex.pta`
        is commented out. The Philips waveform is vendor-confidential and is
        NOT shipped with the public MRSCloud repo (see its README, Remark 3).

        Workaround: copy the bundled universal excitation waveform into the
        per-run workdir under the expected filename and prepend the workdir
        to the Octave search path. This way `io_loadRFwaveform` finds it
        without modifying the third-party submodule. We do this for ALL
        vendors (the user may have already staged their own copy in
        `_stage_user_pulse`, in which case we don't overwrite it).

        TODO upstream: make `load_parameters.m` branch the excWaveform on vendor.
        """
        import os, shutil
        src = os.path.abspath(
            './externals/mrscloud/pulses_universal/univ_spreddenrex.pta'
        )
        dst = os.path.join(workdir, 'Philips_spredrex.pta')
        if not os.path.exists(dst):
            if not os.path.exists(src):
                print(f"  ⚠️  Universal excitation waveform not found at {src} — "
                      f"MRSCloud will likely fail with 'File not found'.")
                return
            shutil.copyfile(src, dst)
        # Universal dual-lobe pulses that load_parameters.m requests under a
        # name the public repo never had (see _PULSE_ALIASES); a user-supplied
        # file of the requested name takes precedence.
        univ_dir = os.path.abspath('./externals/mrscloud/pulses_universal')
        vendor_dir = os.path.abspath('./externals/mrscloud/pulses')
        for wanted, shipped in self._PULSE_ALIASES.items():
            alias_dst = os.path.join(workdir, wanted)
            alias_src = os.path.join(univ_dir, shipped)
            if (not os.path.exists(os.path.join(vendor_dir, wanted))
                    and not os.path.exists(alias_dst) and os.path.exists(alias_src)):
                shutil.copyfile(alias_src, alias_dst)
        # Prepend the workdir so the shim is found *first*. The Octave runtime
        # may be running inside a Docker container that mounts the project root
        # at /workspace, so we MUST pass a path relative to the project root —
        # an absolute host path won't resolve inside the container.
        try:
            rel = os.path.relpath(workdir, start=os.path.abspath('.'))
            shim_path = './' + rel.replace('\\', '/')
            self._shim_path = shim_path
            self.octave.addpath(shim_path)
        except Exception as e:
            print(f"  ⚠️  Could not addpath({workdir}): {e}")

    @staticmethod
    def _column(re, im, i, npts):
        """Column i of the adapter's sub-spectrum matrices (N×k; a single
        column arrives squeezed to a vector)."""
        re = np.asarray(re, dtype=np.float64)
        im = np.asarray(im, dtype=np.float64)
        if re.ndim == 1:
            re, im = re[:, None], im[:, None]
        if re.shape[0] != npts or i >= re.shape[1]:
            raise RuntimeError("sub-spectrum missing from the adapter output")
        return re[:, i] + 1j * im[:, i]

    # --------------------------------------------------------------- main entry
    def run_simulation(self, params, progress_callback=None, stop_event=None):
        """Run MRSCloud per-metabolite and return { metab : 1-D complex FID }."""
        # Lazy Octave init
        if self.octave is None:
            print("Initializing Octave runtime...")
            self.initialize_octave(prefer_docker=True)
        self.setup_octave_paths()

        # Internal scratch (MRSCloud writes intermediate .mat files here)
        workdir = self.ensure_workdir()
        save_dir = workdir
        if os.path.isabs(save_dir):
            try:
                save_dir = os.path.relpath(save_dir)
            except ValueError:
                save_dir = save_dir.replace('\\', '/')

        # Pull params into local strongly-typed variables. The stored System is
        # the plain scanner vendor; translate it (with the current Mode) into
        # the label MRSCloud expects ('Philips' / 'Universal_Philips' / ...).
        vendor       = self._mrscloud_vendor(params.get('System'))
        if not vendor:
            raise RuntimeError(
                "MRSCloud: no System (scanner vendor) selected — choose "
                "Philips, Siemens, or GE before simulating.")
        sequence     = str(params.get('Sequence') or '').strip()
        localization = str(params.get('Localization') or '').strip()
        if not sequence or not localization:
            raise ValueError(
                "MRSCloud: Sequence and Localization must be selected before "
                "simulating.")
        if localization == 'STEAM_7T':
            raise ValueError(
                "MRSCloud: STEAM_7T cannot be simulated from the public MRSCloud "
                "repository — sim_signals_STEAM.m is unfinished upstream (no "
                "mixing time is ever set, and the value is used as the flip "
                "angle of the third pulse). Use FID-A's STEAM (ideal or shaped) "
                "or Vespa's STEAM instead.")

        # Stage the user-supplied vendor pulse (if any) FIRST, then drop in
        # the bundled universal excitation waveform under the name MRSCloud
        # hard-codes (Philips_spredrex.pta) so io_loadRFwaveform can find it
        # for every vendor.
        self._stage_user_pulse(workdir, vendor, sequence, localization,
                               params.get(self._pulse_param_label))
        self._stage_pulse_shims(workdir)
        try:
            bfield = float(params.get('Bfield'))
        except (TypeError, ValueError):
            raise ValueError(
                "MRSCloud: Bfield (T) must be set before simulating — it is "
                "read from the data header or entered by hand.") from None
        field_str    = self._parameter_set(bfield)
        edit_target  = ''   # adapter signature only; MRSCloud has no such input
        edit_on      = float(params.get('Edit On',  1.9))
        edit_off     = float(params.get('Edit Off', 7.5))
        edit_tp      = float(params.get('Edit Tp',  14))
        spatial      = int(float(params.get('Spatial Points', 41)))
        try:
            te = float(params.get('TE'))
        except (TypeError, ValueError):
            raise ValueError(
                "MRSCloud: TE (ms) must be set before simulating.") from None
        samples      = int(float(params.get('Samples') or 0))
        bandwidth    = float(params.get('Bandwidth') or 0)
        if samples <= 0 or bandwidth <= 0:
            raise ValueError(
                f"MRSCloud: Samples ({samples}) and Bandwidth ({bandwidth}) "
                "must be set before simulating.")
        metabs       = list(params.get('Metabolites') or [])
        self.last_failures = {}   # metab -> reason, surfaced by the GUI

        # MRSCloud's own runner fixes TE = 80 ms for HERMES / HERCULES; the
        # adapter passes the requested TE through to load_parameters instead.
        if sequence in ('HERMES', 'HERCULES') and abs(te - 80) > 1:
            print(f"  Note: MRSCloud's runner uses TE = 80 ms for {sequence}; "
                  f"simulating at your TE = {te:g} ms.")

        if not metabs:
            raise ValueError("MRSCloud: no metabolites selected.")

        multi = sequence in ('MEGA', 'HERMES', 'HERCULES')
        total = len(metabs)
        counter = {'n': 0}

        def work(octave, metab):
            """One metabolite -> its basis entries (None when it failed: the
            reason is recorded in last_failures and surfaced by the GUI)."""
            counter['n'] += 1
            print(f"[MRSCloud] {counter['n']}/{total}  simulating {metab} "
                  f"({sequence}/{localization} on {vendor}, TE={te:g} ms, "
                  f"B0={bfield:g} T, {field_str} parameter set)")
            try:
                # edited schemes: the adapter also hands back the further
                # sub-spectra (MEGA: edit-OFF; HERMES / HERCULES: B, C, D)
                out = octave.feval(
                    'mrscloud_run_metab',
                    metab, vendor, sequence, localization,
                    te, field_str, edit_target,
                    edit_on, edit_off, edit_tp, float(spatial), save_dir,
                    float(samples), float(bandwidth), float(bfield),
                    nout=7 if multi else 5,
                )
                fid = np.asarray(out[0], dtype=np.float64).flatten() \
                    + 1j * np.asarray(out[1], dtype=np.float64).flatten()
                if fid.size == 0:
                    raise RuntimeError("empty FID returned")
                entries = {}
                if sequence == 'MEGA':
                    off = self._column(out[5], out[6], 0, fid.size)
                    entries[f'{metab} (ON)'] = fid
                    entries[f'{metab} (OFF)'] = off
                    entries[f'{metab} (DIFF)'] = fid - off
                elif multi:
                    # Hadamard-encoded schemes: the four sub-experiments A–D
                    # (editing pulses at the scheme's offsets, in the order of
                    # the adapter's editON), their SUM and the two difference
                    # spectra named as Osprey does (DIFF1 = GABA, DIFF2 = GSH).
                    subs = {'A': fid}
                    for i, tag in enumerate('BCD'):
                        subs[tag] = self._column(out[5], out[6], i, fid.size)
                    for tag, arr in subs.items():
                        entries[f'{metab} ({tag})'] = arr
                    entries[f'{metab} (SUM)'] = sum(subs.values())
                    a, b, c, d = (subs[t] for t in 'ABCD')
                    if sequence == 'HERMES':
                        # A 4.56 (GSH on), B 1.90 (GABA on), C dual (both), D 7.5
                        # (neither) — measured: (B + C) − (A + D) isolates the
                        # 3 ppm GABA signal (0.69 of SUM; upstream's own
                        # B + D − A − C gives 0.02 with this pulse order)
                        entries[f'{metab} (DIFF1)'] = (b + c) - (a + d)
                        entries[f'{metab} (DIFF2)'] = (a + c) - (b + d)
                    else:
                        # HERCULES: A 4.58 (GSH on), B 4.18, C dual 4.58 + 1.9,
                        # D dual 4.18 + 1.9 — GABA is edited in C and D, GSH in
                        # A and C (from the pulse assignment; not measured)
                        entries[f'{metab} (DIFF1)'] = (c + d) - (a + b)
                        entries[f'{metab} (DIFF2)'] = (a + c) - (b + d)
                else:
                    entries[metab] = fid
                return entries
            except Exception as e:
                # Don't kill the whole run — record the failure (surfaced by
                # the GUI) and skip; a zero-filled FID would look like success.
                print(f"  ✗ {metab}: {e}")
                self.last_failures[metab] = str(e)
                return {}

        # metabolites run concurrently over several Octave processes
        per_metab = self.simulate_in_parallel(metabs, work, progress_callback, stop_event)
        basis_set: dict[str, np.ndarray] = {}
        for entries in per_metab.values():
            basis_set.update(entries)
        return basis_set

















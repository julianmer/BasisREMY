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
# Purpose: BasisREMY backend that drives the **real** MRSCloud workflow                            #
#          (externals/mrscloud) instead of the wrong FID-A `sim_lcmrawbasis`                       #
#          script that was used previously. Per-metabolite execution goes                          #
#          through `adapters/backends/mrscloud_run_metab.m`, a thin Octave                         #
#          adapter that mirrors the body of                                                        #
#          `externals/mrscloud/run/run_simulations_cloud.m` and returns                            #
#          the FID as plain numeric arrays so oct2py can ferry them across.                        #
#                                                                                                  #
# Notes                                                                                            #
#   - Output Path / Make .raw / Add Ref. are no longer user-facing; final                          #
#     export is handled by the post-simulation Export dialog (core/exporters).                     #
#   - The metabolite list matches the official MRSCloud README. Some entries                       #
#     (Cystat, HCar, iLe, Lys, Glc) simulate slowly — flagged with TODO.                           #
#   - HERMES / HERCULES sub-spectrum recombination is not yet exposed; only                        #
#     the first sub-spectrum (off / 'a') is returned. TODO marked below.                           #
#                                                                                                  #
####################################################################################################


#*************#
#   imports   #
#*************#
from __future__ import annotations

import os

import numpy as np

# own
from backends.base import Backend


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

    _LOCALIZATIONS = ['PRESS', 'sLASER', 'STEAM_7T']
    _VENDORS       = ['Philips', 'Universal_Philips', 'Siemens', 'Universal_Siemens', 'GE']

    # Combined Sequence dropdown: each label encodes both the MRSCloud editing
    # scheme *and* the localization in one user-facing string.  The GUI only
    # ever shows this one field; the backend decodes it into the two internal
    # values MRSCloud's Octave layer expects.
    _COMBINED_SEQ_OPTIONS = [
        'PRESS',
        'sLASER',
        'STEAM (7T only)',
        'MEGA-PRESS',
        'MEGA-sLASER',
        'HERMES-PRESS',
        'HERMES-sLASER',
        'HERCULES-PRESS',
        'HERCULES-sLASER',
    ]
    # label → (mrscloud_sequence, mrscloud_localization)
    _COMBINED_SEQ_MAP = {
        'PRESS':            ('UnEdited',  'PRESS'),
        'sLASER':           ('UnEdited',  'sLASER'),
        'STEAM (7T only)':  ('UnEdited',  'STEAM_7T'),
        'MEGA-PRESS':       ('MEGA',      'PRESS'),
        'MEGA-sLASER':      ('MEGA',      'sLASER'),
        'HERMES-PRESS':     ('HERMES',    'PRESS'),
        'HERMES-sLASER':    ('HERMES',    'sLASER'),
        'HERCULES-PRESS':   ('HERCULES',  'PRESS'),
        'HERCULES-sLASER':  ('HERCULES',  'sLASER'),
    }
    # Reverse map for when we need to go from internals → display label
    _COMBINED_SEQ_REVERSE = {v: k for k, v in _COMBINED_SEQ_MAP.items()}

    # ---- Pulse-file availability map ---------------------------------
    # MRSCloud's README warns: "Product sequence and rf waveform are not
    # shared in the GitHub repo." That means *most* vendor pulse files are
    # NOT bundled — only the Universal_* `univ_*.pta` waveforms ship.
    #
    # The map below is the EXACT set of `.pta` waveforms that
    # `externals/mrscloud/functions/load_parameters.m` requests via
    # `io_loadRFwaveform` for a given (vendor, sequence, localization)
    # triplet. Anything not present on disk is surfaced to the user as a
    # "Browse" file picker so they can supply the vendor-confidential file
    # themselves (and we copy it into the per-run workdir under the right
    # name before MRSCloud is invoked).
    #
    # Path is relative to externals/mrscloud/. Universal_* combos are fully
    # bundled — `missing_pulse_files()` returns [] for them.
    _PULSE_FILES_BY_VENDOR_LOC_SEQ = {
        # ------------- Philips -------------
        ('Philips',          'PRESS',  'UnEdited'): ['pulses/Philips_spredrex.pta',
                                                     'pulses/gtst1203_sp.pta'],
        ('Philips',          'PRESS',  'MEGA'):     ['pulses/Philips_spredrex.pta',
                                                     'pulses/gtst1203_sp.pta'],
        ('Philips',          'sLASER', 'UnEdited'): ['pulses/Philips_spredrex.pta',
                                                     'pulses/Philips_GOIA_WURST_100pts.mat'],
        # ------------- Siemens -------------
        ('Siemens',          'PRESS',  'UnEdited'): ['pulses/Philips_spredrex.pta',  # excitation hard-coded
                                                     'pulses/orig_refoc_mao_100_4.pta'],
        ('Siemens',          'PRESS',  'MEGA'):     ['pulses/Philips_spredrex.pta',
                                                     'pulses/orig_refoc_mao_100_4.pta'],
        # ------------- GE -------------
        ('GE',               'PRESS',  'UnEdited'): ['pulses/Philips_spredrex.pta',
                                                     'pulses/GE_rfa_3.9ms.pta'],
        # ------------- Universal_* (bundled) -------------
        ('Universal_Philips','PRESS',  'UnEdited'): [],   # univ_* shipped + shimmed
        ('Universal_Philips','PRESS',  'MEGA'):     [],
        ('Universal_Philips','PRESS',  'HERMES'):   [],
        ('Universal_Philips','PRESS',  'HERCULES'): [],
        ('Universal_Siemens','PRESS',  'UnEdited'): [],
        ('Universal_Siemens','PRESS',  'MEGA'):     [],
        # STEAM_7T uses inline pulses
        ('Universal_Philips','STEAM_7T','UnEdited'): [],
        ('Universal_Siemens','STEAM_7T','UnEdited'): [],
    }

    @classmethod
    def _decode_sequence(cls, combined: str) -> tuple[str, str]:
        """Return (mrscloud_sequence, localization) for a combined label.

        Falls back to ('UnEdited', 'PRESS') for unknown / legacy values so
        that stale params from other backends don't crash on startup.
        """
        return cls._COMBINED_SEQ_MAP.get(combined, ('UnEdited', 'PRESS'))

    @classmethod
    def required_pulse_files(cls, vendor: str, sequence: str, localization: str) -> list[str]:
        """Pulse files MRSCloud will demand for this (vendor, seq, loc) combo."""
        return list(cls._PULSE_FILES_BY_VENDOR_LOC_SEQ.get(
            (vendor, localization, sequence), []))

    @classmethod
    def missing_pulse_files(cls, vendor: str, sequence: str, localization: str,
                            mrscloud_root: str = './externals/mrscloud') -> list[str]:
        """Subset of required pulse files that are NOT present on disk."""
        import os
        out = []
        for rel in cls.required_pulse_files(vendor, sequence, localization):
            if not os.path.exists(os.path.join(mrscloud_root, rel)):
                out.append(rel)
        return out

    def __init__(self):
        super().__init__()
        self.name = 'MRSCloud'
        self.display_name = 'MRSCloud'
        self.category = 'MRSCloud'
        self.requires_octave = True

        # Metabolite library
        self.metabs = dict(self._METABS_DEFAULT)

        # Dropdowns shown in the GUI. NOTE: the GUI rebuilds the parameter
        # panel whenever a key in `schema_affecting_keys` changes, so the
        # editing fields appear/disappear depending on `Sequence`.
        self.dropdown = {
            'System':         list(self._VENDORS),
            # One combined field: localization is encoded in the option name.
            # get_params_for_mode() filters out 'STEAM (7T only)' at non-7T.
            'Sequence':       list(self._COMBINED_SEQ_OPTIONS),
            'Field Strength': ['1.5T', '3T', '7T'],
            'Edit Target':    ['', 'GABA', 'GSH', 'Lac', 'PE'],
        }

        # GUI tells us when these change so we can rebuild the visible
        # parameter list (e.g. show editing fields only for MEGA/HERMES/HERCULES,
        # show the pulse-file picker only when the pulse is missing).
        self.schema_affecting_keys = {'Sequence', 'System', 'Field Strength'}

        # Pulse-file pickers are populated dynamically by get_params_for_mode().
        self.file_selection: list[str] = []

        # Mandatory parameters — only the ones MRSCloud actually consumes.
        # 'Sequence' now encodes localization too (e.g. 'PRESS', 'MEGA-PRESS').
        # Removed:
        #   * Localization  — encoded inside Sequence combined label
        #   * Bfield        — derived from Field Strength + vendor inside
        #                     externals/mrscloud/functions/load_parameters.m
        #   * Center Freq   — same (computed from Bfield × γ inside FID-A)
        #   * Linewidth     — MRSCloud hard-codes lw = 1 Hz in load_parameters
        #   * TE2           — MRSCloud uses a single TE; TE1 is set per-vendor
        # Editing-only fields (Edit Target / On / Off / Tp) live in
        # `_edit_params` and are spliced in by get_params_for_mode().
        self.mandatory_params = {
            'System':         None,        # vendor (must be selected)
            'Sequence':       None,        # combined seq+loc (must be selected)
            'Field Strength': '3T',
            'Samples':        None,
            'Bandwidth':      None,
            'TE':             None,
            'Spatial Points': 41,          # 41 acceptable, 101 ideal (slow)
            'Metabolites':    [k for k, v in self.metabs.items() if v],
        }

        # Editing-sequence-only parameters (added to mandatory_params on the
        # fly). MRSCloud overrides editON internally for HERMES/HERCULES, so
        # those sequences only expose Edit Target + Edit Tp.
        self._edit_params = {
            'Edit Target':    '',          # GABA / GSH / Lac / PE
            'Edit On':        1.9,         # ppm (MEGA only)
            'Edit Off':       7.5,         # ppm (MEGA only)
            'Edit Tp':        14,          # ms — 14 for MEGA/HERMES, 20 for HERCULES
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
    def get_params_for_mode(self, mode=None):
        """Return only the parameters relevant to the current Sequence /
        System / Field Strength combo so the GUI never shows fields MRSCloud
        will silently ignore."""
        params = dict(self.mandatory_params)
        combined = (self.mandatory_params.get('Sequence') or '').strip()
        vendor   = (self.mandatory_params.get('System')        or '').strip()
        field    = (self.mandatory_params.get('Field Strength') or '3T').strip()

        # Decode combined label → internal (sequence, localization)
        seq_internal, loc_internal = self._decode_sequence(combined)

        # ---- filter out STEAM (7T only) when not at 7T ----------------------
        if field == '7T':
            all_combined = list(self._COMBINED_SEQ_OPTIONS)
        else:
            all_combined = [o for o in self._COMBINED_SEQ_OPTIONS
                            if o != 'STEAM (7T only)']
            if combined == 'STEAM (7T only)':
                # Reset to PRESS (safe fallback)
                self.mandatory_params['Sequence'] = 'PRESS'
                params['Sequence'] = 'PRESS'
                combined = 'PRESS'
                seq_internal, loc_internal = 'UnEdited', 'PRESS'

        self.dropdown['Sequence'] = all_combined

        # ---- splice editing fields only when relevant -----------------------
        if seq_internal == 'MEGA':
            for k, v in self._edit_params.items():
                params.setdefault(k, v)
                self.mandatory_params.setdefault(k, v)
        elif seq_internal in ('HERMES', 'HERCULES'):
            for k in ('Edit Target', 'Edit Tp'):
                params.setdefault(k, self._edit_params[k])
                self.mandatory_params.setdefault(k, self._edit_params[k])
            for k in ('Edit On', 'Edit Off'):
                params.pop(k, None)
                self.mandatory_params.pop(k, None)
        else:
            # UnEdited → hide everything edit-related
            for k in self._edit_params:
                params.pop(k, None)
                self.mandatory_params.pop(k, None)

        # ---- ask for the missing vendor pulse file when needed ---------------
        self.file_selection = []
        if combined and vendor and vendor not in ('Universal_Philips', 'Universal_Siemens'):
            missing = self.missing_pulse_files(vendor, seq_internal, loc_internal)
            if missing:
                label = self._pulse_param_label
                self.file_selection.append(label)
                params[label] = self.mandatory_params.get(label)
                self.mandatory_params.setdefault(label, None)
        else:
            self.mandatory_params.pop(self._pulse_param_label, None)

        # Re-order Metabolites to the bottom for the GUI grid
        if 'Metabolites' in params:
            mets = params.pop('Metabolites')
            params['Metabolites'] = mets
        return params

    # ----------------------------------------------------------------- REMY parsing
    def parseREMY(self, MRSinMRS):
        """Map REMY-extracted metadata onto MRSCloud parameters.

        Note: `Bfield` and `Center Freq` are intentionally NOT returned —
        MRSCloud derives them internally from `Field Strength` + vendor
        (see externals/mrscloud/functions/load_parameters.m, lines 20-36).
        """
        bfield = MRSinMRS.get('B0', None)

        mandatory = {
            'Sequence':       self.parseProtocol(MRSinMRS.get('Protocol', None)),
            'System':         self.parseSystem(MRSinMRS.get('Manufacturer', None)),
            'Field Strength': self._field_str_from_b0(bfield),
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
    def _field_str_from_b0(b0):
        if b0 is None:
            return '3T'
        try:
            b0 = float(b0)
        except (TypeError, ValueError):
            return '3T'
        if abs(b0 - 1.5) < 0.3:
            return '1.5T'
        if abs(b0 - 7.0) < 0.5:
            return '7T'
        return '3T'

    def parseProtocol(self, protocol):
        """Return a combined Sequence label (e.g. 'PRESS', 'MEGA-PRESS') from
        the raw protocol string.  This single field encodes both the MRSCloud
        editing scheme *and* the localization so the GUI only needs one
        dropdown."""
        if protocol is None:
            return None
        p = str(protocol).lower()
        is_slaser = ('slaser' in p or 'semi_laser' in p
                     or 'semi-laser' in p or 'semilaser' in p)
        # Editing scheme (highest priority), with localization suffix
        if 'hercules' in p:
            return 'HERCULES-sLASER' if is_slaser else 'HERCULES-PRESS'
        if 'hermes' in p:
            return 'HERMES-sLASER' if is_slaser else 'HERMES-PRESS'
        if 'mega' in p:
            return 'MEGA-sLASER' if is_slaser else 'MEGA-PRESS'
        # Localization-only (→ UnEdited)
        if is_slaser:
            return 'sLASER'
        if 'steam' in p:
            return 'STEAM (7T only)'
        # PRESS, UnEdited (BigGABA/MRSCloud convention), LASER all default to PRESS
        if any(tag in p for tag in ('press', 'unedited', 'laser')):
            return 'PRESS'
        print(f"Warning: MRSCloud could not infer sequence from '{protocol}'.")
        return None

    def parseSystem(self, system):
        """Map a vendor string onto the MRSCloud vendor list."""
        if system is None:
            return None
        s = str(system).lower()
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

    # --------------------------------------------------------------- Octave paths
    def setup_octave_paths(self):
        if self.octave is None:
            raise RuntimeError("Octave not initialized. Call initialize_octave() first.")
        self.octave.eval("warning('off', 'all');")
        # adapters/backends contains our mrscloud_run_metab.m wrapper
        self.octave.addpath('./adapters/backends/')
        # Pull in MRSCloud (functions + bundled FID-A) recursively
        self.octave.addpath(self.octave.genpath('./externals/mrscloud/functions/'))
        self.octave.addpath(self.octave.genpath('./externals/mrscloud/pulses_universal/'))

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
        # Drop the excitation Philips_spredrex.pta — that one is shimmed
        # automatically by `_stage_universal_excite_shim` below.
        missing = [m for m in missing if not m.endswith('Philips_spredrex.pta')]
        if not missing:
            return
        target_name = os.path.basename(missing[0])
        dst = os.path.join(workdir, target_name)
        try:
            if os.path.abspath(user_path) != os.path.abspath(dst):
                shutil.copyfile(user_path, dst)
            print(f"  ✓ Staged user pulse '{os.path.basename(user_path)}' "
                  f"as '{target_name}' in workdir")
        except Exception as e:
            print(f"  ⚠️  Could not stage user pulse {user_path}: {e}")

    def _stage_universal_excite_shim(self, workdir: str) -> None:
        """Work around the hard-coded `Philips_spredrex.pta` excitation pulse.

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
        # Prepend the workdir so the shim is found *first*. The Octave runtime
        # may be running inside a Docker container that mounts the project root
        # at /workspace, so we MUST pass a path relative to the project root —
        # an absolute host path won't resolve inside the container.
        try:
            rel = os.path.relpath(workdir, start=os.path.abspath('.'))
            shim_path = './' + rel.replace('\\', '/')
            self.octave.addpath(shim_path)
        except Exception as e:
            print(f"  ⚠️  Could not addpath({workdir}): {e}")

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

        # Pull params into local strongly-typed variables.
        # 'Sequence' is the combined label (e.g. 'PRESS', 'MEGA-PRESS').
        # Decode it into the two internal values MRSCloud's Octave layer expects.
        vendor   = str(params.get('System')        or 'Philips')
        combined = str(params.get('Sequence')      or 'PRESS')
        sequence, localization = self._decode_sequence(combined)

        # Stage the user-supplied vendor pulse (if any) FIRST, then drop in
        # the bundled universal excitation waveform under the name MRSCloud
        # hard-codes (Philips_spredrex.pta) so io_loadRFwaveform can find it
        # for every vendor.
        self._stage_user_pulse(workdir, vendor, sequence, localization,
                               params.get(self._pulse_param_label))
        self._stage_universal_excite_shim(workdir)
        field_str    = str(params.get('Field Strength')or '3T')
        edit_target  = str(params.get('Edit Target')   or '')
        edit_on      = float(params.get('Edit On',  1.9))
        edit_off     = float(params.get('Edit Off', 7.5))
        edit_tp      = float(params.get('Edit Tp',  14))
        spatial      = int(float(params.get('Spatial Points', 41)))
        te           = float(params.get('TE') or 35)
        samples      = int(float(params.get('Samples') or 0))
        bandwidth    = float(params.get('Bandwidth') or 0)
        metabs       = list(params.get('Metabolites') or [])

        # MRSCloud overrides TE for HERMES (68) / HERCULES (80) internally;
        # we keep the user value for documentation but warn if it's unusual.
        if sequence == 'HERMES' and abs(te - 68) > 1:
            print(f"  Note: MRSCloud will internally use TE=80 ms for HERMES "
                  f"(your TE={te} ms is informational only).")
        if sequence == 'HERCULES' and abs(te - 80) > 1:
            print(f"  Note: MRSCloud will internally use TE=80 ms for HERCULES "
                  f"(your TE={te} ms is informational only).")

        if not metabs:
            raise ValueError("MRSCloud: no metabolites selected.")

        basis_set: dict[str, np.ndarray] = {}
        total = len(metabs)
        for i, metab in enumerate(metabs):
            if stop_event and stop_event.is_set():
                print(f"  ⏹  Stopped before simulating {metab} (user cancelled).")
                break
            print(f"[MRSCloud] {i+1}/{total}  simulating {metab} "
                  f"({sequence}/{localization} on {vendor}, TE={te} ms, B0={field_str})")
            try:
                fid_re, fid_im, npts, _sw, _cf = self.octave.feval(
                    'mrscloud_run_metab',
                    metab, vendor, sequence, localization,
                    te, field_str, edit_target,
                    edit_on, edit_off, edit_tp, float(spatial), save_dir,
                    float(samples), float(bandwidth),
                    nout=5,
                )
                fid = np.asarray(fid_re, dtype=np.float64).flatten() \
                    + 1j * np.asarray(fid_im, dtype=np.float64).flatten()
                if fid.size == 0:
                    raise RuntimeError("empty FID returned")
                basis_set[metab] = fid
            except Exception as e:
                # Don't kill the whole run — log, store an empty FID, continue.
                # TODO surface this in the GUI summary instead of just printing.
                print(f"  ✗ {metab}: {e}")
                basis_set[metab] = np.zeros(int(params.get('Samples') or 2048),
                                            dtype=np.complex128)

            if progress_callback:
                progress_callback(i + 1, total)

        return basis_set

















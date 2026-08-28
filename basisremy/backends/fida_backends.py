####################################################################################################
#                                         fida_backends.py                                         #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 25/04/26                                                                                #
#                                                                                                  #
# Purpose: One-stop module for the entire "FID-A" backend family. Each FID-A simulation entry      #
#          point under externals/fidA/simulationTools/ is exposed as a small subclass of           #
#          FidaBackend. They all share:                                                            #
#            * the spinSystems.mat metabolite library                                              #
#            * Octave path setup                                                                   #
#            * REMY → param parsing                                                                #
#            * the per-metabolite driver loop                                                      #
#                                                                                                  #
#          Only the parameter schema and the `kind` dispatched into the shared Octave adapter      #
#          (adapters/backends/fida_run.m) differ between subclasses. FidaIdeal is the              #
#          ex-"LCModel" backend, now living natively under FID-A.                                  #
#                                                                                                  #
####################################################################################################


from __future__ import annotations

import os

import numpy as np

from basisremy.backends.base import Backend


# --------------------------------------------------------------------------- defaults
# Default metabolite map (matches the entries available in
# externals/fidA/simulationTools/metabolites/spinSystems.mat).
_DEFAULT_FIDA_METABS = {
    'Ala': False, 'Asc': True,  'Asp': False, 'Ch':  False, 'Cit': False,
    'Cr':  True,  'EtOH': False,'GABA': True, 'GPC': True,  'GSH': True,
    'Glc': True,  'Gln': True,  'Glu': True,  'Gly': True,  'H2O': False,
    'Ins': True,  'Lac': True,  'Lip': False, 'NAA': True,  'NAAG': True,
    'PCh': True,  'PCr': True,  'PE':  True,  'Phenyl': False,
    'Ref0ppm': False, 'Scyllo': True, 'Ser': False, 'Tau': True, 'Tyros': False,
}


def _shaped_params(extra: dict | None = None) -> dict:
    """Common parameter sheet for shaped 2-D-localised FID-A sims."""
    base = {
        'Samples':       None,
        'Bandwidth':     None,
        'Bfield':        None,
        'Linewidth':     1.0,
        'TE':            None,
        'Flip Angle':    180.0,
        'RefTp':         5.0,
        'thkX':          2.0,
        'thkY':          2.0,
        'fovX':          3.0,
        'fovY':          3.0,
        'nX':            8,
        'nY':            8,
        'Sim Centre (ppm)': 4.65,
        'Path to Pulse': None,
        'Metabolites':   [],
    }
    if extra:
        base.update(extra)
    return base


# =================================================================== base class
class FidaBackend(Backend):
    """Common scaffolding for FID-A simulation wrappers.

    Subclasses set:
      * ``self.name`` / ``self.display_name``
      * ``self._kind``: string passed to ``fida_run.m`` to pick the simulator
      * ``self.mandatory_params`` (and optionally ``self.dropdown`` /
        ``self.file_selection`` / ``self.modes``)
      * ``self._build_args(params, metab)``: positional args for ``fida_run``
        AFTER the metabolite name and the ``kind`` argument

    A subclass with ``self._kind == ''`` is considered a stub — the GUI shows
    it, but ``run_simulation`` raises NotImplementedError until the matching
    branch is added inside ``adapters/backends/fida_run.m``.
    """

    _kind: str = ''      # dispatch key for fida_run.m
    _is_stub: bool = False

    def __init__(self):
        super().__init__()
        self.category = 'FID-A'
        self.requires_octave = True
        self.metabs = dict(_DEFAULT_FIDA_METABS)
        self.optional_params = {'Nucleus': None, 'TR': None}

    # -------------------------------------------------- sequence mapping
    def map_sequence_in(self, seq: str) -> 'str | None':
        """Translate any sequence name into this FID-A backend's vocabulary.
        Subclasses (FidaIdeal, FidaPressShaped …) override as needed."""
        options = self.dropdown.get('Sequence', [])
        if not options or not seq:
            return None
        s = seq.strip().lower()
        for opt in options:
            if opt.lower() == s:
                return opt
        return None

    # -------------------------------------------------- helpers
    def _refresh_metab_list(self):
        active = [k for k, v in self.metabs.items() if v]
        if 'Metabolites' in self.mandatory_params:
            self.mandatory_params['Metabolites'] = active

    @staticmethod
    def _make_relative(path):
        """Convert absolute paths → relative for Docker-Octave compat."""
        if path and isinstance(path, str) and os.path.isabs(path):
            try:
                return os.path.relpath(path)
            except ValueError:
                return path.replace('\\', '/')
        return path

    # -------------------------------------------------- Octave
    def setup_octave_paths(self):
        if self.octave is None:
            raise RuntimeError("Octave not initialized.")
        # Fetch FID-A on first use (no-op in a source checkout).
        from basisremy.core.externals import ensure
        from basisremy.core.paths import octave_adapters_base
        ensure('fidA')
        adapters_base = octave_adapters_base(self.octave)
        self.octave.eval("warning('off', 'all');")
        # First add the FID-A tree recursively so nested helpers (e.g.
        # rfPulseTools/mklassenTools/bes.m, used by io_loadRFwaveform for
        # phase-modulated waveforms like GOIA) are resolvable. Without this,
        # shaped-pulse sims fail with "error: 'bes' undefined".
        self.octave.eval("addpath(genpath('./externals/fidA/'));")
        # THEN add our adapter dirs — addpath() prepends, so these now win
        # over the upstream FID-A files. We use this to:
        #   * ship a patched sim_lcmrawbasis.m
        #   * ship a non-interactive io_loadRFwaveform.m (the upstream one
        #     calls plot()/input() for phase-modulated pulses, which fails
        #     in headless Docker Octave with "ft_text_renderer: invalid
        #     bounding box, cannot render, unable to create graphics handle").
        self.octave.addpath(adapters_base + '/backends/')

    # -------------------------------------------------- REMY
    def parseREMY(self, MRSinMRS):
        mandatory = {
            'Samples':     MRSinMRS.get('NumberOfDatapoints', None),
            'Bandwidth':   MRSinMRS.get('SpectralWidth', None),
            'Bfield':      MRSinMRS.get('B0', None),
            'TE':          MRSinMRS.get('TE', None),
            'Center Freq': MRSinMRS.get('Center Freq', None),
        }
        # Only keys this backend actually exposes.
        mandatory = {k: v for k, v in mandatory.items()
                     if k in self.mandatory_params}
        optional = {
            'Nucleus':         MRSinMRS.get('Nucleus', None),
            'TR':              MRSinMRS.get('TR', None),
            'Model':           MRSinMRS.get('Model', None),
            'SoftwareVersion': MRSinMRS.get('SoftwareVersion', None),
            'BodyPart':        MRSinMRS.get('BodyPart', None),
        }
        return mandatory, optional

    def parseProtocol(self, protocol):
        return protocol

    # -------------------------------------------------- edited sub-spectra
    def _run_on_off_subspectra(self, params, base_args,
                               progress_callback=None, stop_event=None):
        """Run each metabolite twice (edit ON / OFF) and return the flat
        '<metab> (ON/OFF/DIFF)' basis convention used by edited backends."""
        metabs = params.get('Metabolites') or []
        basis = {}
        for i, metab in enumerate(metabs):
            if stop_event and stop_event.is_set():
                print(f"  ⏹  Stopped before simulating {metab}.")
                break
            fids = {}
            for label, flag in (('ON', 1), ('OFF', 0)):
                results = self.octave.feval(
                    'fida_run', metab, self.active_kind(), *base_args, flag, nout=5,
                )
                fid_re, fid_im, _npts, _sw, _cf = results
                fids[label] = (np.asarray(fid_re, dtype=float).flatten()
                               + 1j * np.asarray(fid_im, dtype=float).flatten())
            basis[f'{metab} (ON)'] = fids['ON']
            basis[f'{metab} (OFF)'] = fids['OFF']
            basis[f'{metab} (DIFF)'] = fids['ON'] - fids['OFF']
            if progress_callback:
                progress_callback(i + 1, len(metabs))
        return basis

    # -------------------------------------------------- per-subclass hook
    def _build_args(self, params, metab):
        """Positional args (AFTER metab + kind) for ``fida_run.m``."""
        raise NotImplementedError

    # -------------------------------------------------- driver
    def active_kind(self):
        """`fida_run.m` kind for the current mode (subclasses with modes
        override this)."""
        return self._kind

    def run_simulation(self, params, progress_callback=None, stop_event=None):
        if self._is_stub or not self._kind:
            raise NotImplementedError(
                f"{self.name}: this FID-A wrapper is a stub. The schema is "
                "complete but the matching branch in adapters/backends/"
                "fida_run.m has not been implemented yet. See FidaIdeal / "
                "FidaPressShaped (kinds 'ideal' / 'press_shaped') for the "
                "canonical reference."
            )

        if self.octave is None:
            print("Initializing Octave runtime...")
            self.initialize_octave(prefer_docker=True)
        self.setup_octave_paths()
        self.ensure_workdir()

        metabs = params.get('Metabolites') or []
        basis = {}
        for i, metab in enumerate(metabs):
            if stop_event and stop_event.is_set():
                print(f"  ⏹  Stopped before simulating {metab}.")
                break
            extra_args = self._build_args(params, metab)
            results = self.octave.feval(
                'fida_run', metab, self.active_kind(), *extra_args, nout=5,
            )
            fid_re, fid_im, _npts, _sw, _cf = results
            fid = (np.asarray(fid_re, dtype=float).flatten()
                   + 1j * np.asarray(fid_im, dtype=float).flatten())
            # FID-A's sim_readout stores `out.specs = fftshift(ifft(out.fids))`
            # with a ppm axis `ppm = -freq/larmor + 4.65`.  fida_run.m returns
            # out.fids directly (no conjugation applied), so the FID oscillates
            # at -(δ - centreFreq)*larmor Hz for a metabolite at δ ppm.  Our
            # GUI computes `fftshift(fft(fid))` and uses a ppm axis
            # `+freq/larmor + 4.65`.  fft of a −f0 signal peaks at −f0 →
            # maps to (−f0/larmor + 4.65) ppm — which correctly equals δ ppm
            # when centreFreq = 4.65.  No conjugation needed here.
            basis[metab] = fid
            if progress_callback:
                progress_callback(i + 1, len(metabs))
        return basis


# =================================================================== Ideal (ex-LCModel)
class FidaIdeal(FidaBackend):
    """sim_lcmrawbasis: ideal-pulse Spin Echo / PRESS / STEAM / LASER.

    Canonical "Ideal" entry under the FID-A category. This is the renamed
    successor of the historical ``LCModelBackend``: same simulator (FID-A's
    ``sim_lcmrawbasis``), different home.
    """

    _kind = 'ideal'

    def __init__(self):
        super().__init__()
        self.name = 'FidaIdeal'
        self.display_name = 'Ideal (SE / PRESS / STEAM / LASER)'

        self.dropdown = {
            'Sequence': ['Spin Echo', 'PRESS', 'STEAM', 'LASER'],
        }
        self.mandatory_params = {
            'Sequence':    None,
            'Samples':     None,
            'Bandwidth':   None,
            'Bfield':      None,
            'Linewidth':   1,
            'TE':          None,
            'TE2':         0,
            'TM':          10,
            'Metabolites': [],
            'Center Freq': None,
        }
        # TM is only shown for STEAM — rebuild the panel on Sequence changes.
        self.schema_affecting_keys = {'Sequence'}
        self._refresh_metab_list()

    def get_params_for_mode(self, mode=None):
        # TM applies to STEAM only; hide it for the other sequences.
        params = dict(self.mandatory_params)
        if params.get('Sequence') != 'STEAM':
            params.pop('TM', None)
        return params

    # ---- sequence mapping ----------------------------------------------
    def map_sequence_in(self, seq: str) -> 'str | None':
        """Translate an arbitrary sequence name into FidaIdeal's vocabulary."""
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
        if 'press' in s or 'spin echo' in s or 'spinecho' in s or s == 'se':
            return 'PRESS'
        if 'laser' in s and 'mega' not in s and 'slaser' not in s and 'semi' not in s:
            return 'LASER'
        # sLASER, MEGA-*, HERMES, HERCULES — no equivalent in FidaIdeal
        return None

    # ---- REMY ---------------------------------------------------------
    def parseREMY(self, MRSinMRS):
        mandatory, optional = super().parseREMY(MRSinMRS)
        mandatory['Sequence'] = self.parseProtocol(MRSinMRS.get('Protocol', None))
        # extra optional fields used by the export dialog
        for k in ('Manufacturer', 'NumberOfAverages', 'WaterSuppression',
                  'BodyPart', 'VOI', 'AnteriorPosteriorSize', 'LeftRightSize',
                  'CranioCaudalSize'):
            if k in MRSinMRS:
                optional[k.replace('Manufacturer', 'System')] = MRSinMRS[k]
        return mandatory, optional

    def parseProtocol(self, protocol):
        if protocol is None:
            return None
        p = str(protocol).lower()
        if 'mega' in p:
            print("Warning: FidaIdeal does not support MEGA editing — ignoring.")
        if 'slaser' in p:
            print("Warning: FidaIdeal does not support sLASER. Switch backend.")
            return None
        if 'press' in p:      return 'PRESS'
        if 'steam' in p:      return 'STEAM'
        # 'laser' must be tested before the bare 'se' substring — "laser"
        # contains "se", which used to shadow this branch entirely.
        if 'laser' in p:      return 'LASER'
        if 'spin' in p or 'se' in p: return 'Spin Echo'
        # 'UnEdited' is MRSCloud / BigGABA convention for a plain (non-edited)
        # acquisition — default to PRESS, which is by far the most common.
        if 'unedited' in p:   return 'PRESS'
        return None

    # ---- args ---------------------------------------------------------
    _VALID_SEQUENCES = {'Spin Echo': 'se', 'PRESS': 'p', 'STEAM': 'st', 'LASER': 'l'}

    @staticmethod
    def _seq_to_fida(seq):
        mapped = FidaIdeal._VALID_SEQUENCES.get(seq)
        if mapped is None:
            valid = list(FidaIdeal._VALID_SEQUENCES.keys())
            raise ValueError(
                f"FidaIdeal: unrecognised Sequence '{seq}'. "
                f"Please choose one of: {valid}. "
                f"If REMY set this from the file's Protocol (e.g. 'UnEdited'), "
                f"the protocol wasn't recognised — select the sequence manually."
            )
        return mapped

    def _build_args(self, params, metab):
        # workdir for the FID-A-side .RAW writes (kept for parity with the
        # original sim_lcmrawbasis flow; the adapter ignores it but it
        # keeps the path structure consistent across runs).
        out = self._make_relative(self.ensure_workdir()) + os.sep
        seq = params['Sequence']
        te = float(params['TE'])
        if seq == 'STEAM':
            tau1, tau2 = te, float(params.get('TM') or 10.0)   # sim_steam: (TE, TM)
        elif seq == 'PRESS':
            # sim_press evolves tau1/2 - 180 - (tau1+tau2)/2 - 180 - tau2/2 with
            # TE = tau1 + tau2. 'TE2' is the second echo; blank or 0 means a
            # symmetric PRESS (TE/2 each). Passing (TE, 0) — the old behaviour —
            # put the second 180 right before the readout, i.e. a spin echo:
            # singlets were unaffected, strongly coupled spins (Glu) were not.
            te2 = float(params.get('TE2') or 0)
            tau2 = te2 if te2 > 0 else te / 2.0
            tau1 = te - tau2
        else:
            tau1, tau2 = te, 0.0   # spin echo / LASER use tau1 only
        return [
            float(params['Samples']),
            float(params['Bandwidth']),
            float(params['Bfield']),
            float(params.get('Linewidth') or 1),
            tau1,
            tau2,
            self._seq_to_fida(seq),
            out,
        ]


# =================================================================== PRESS shaped
class FidaPressShaped(FidaBackend):
    """PRESS with a shaped refocusing pulse on a spatial grid — FID-A's fast
    method (Zhang 2017): first pulse over the x grid, second over the y grid,
    nX + nY runs with coherence-order filtering instead of nX × nY."""

    _kind = 'press_shaped'

    def __init__(self):
        super().__init__()
        self.name = 'FidaPressShaped'
        self.display_name = 'PRESS shaped'
        self.file_selection = ['Path to Pulse']
        self.mandatory_params = _shaped_params({
            'Tau 1': None,   # ms; defaults to TE/2 if blank
            'Tau 2': None,
        })
        # Move Tau 1/Tau 2 right after TE for a nicer GUI ordering.
        self._refresh_metab_list()

    def parseProtocol(self, protocol):
        if protocol is None:
            return None
        return 'PRESS' if 'press' in str(protocol).lower() else None

    def _build_args(self, params, metab):
        te   = params.get('TE')
        te   = float(te) if te is not None else None
        tau1 = params.get('Tau 1')
        tau2 = params.get('Tau 2')
        if tau1 in (None, ''): tau1 = (te / 2.0) if te is not None else 15.0
        if tau2 in (None, ''): tau2 = (te / 2.0) if te is not None else 15.0

        pulse_src = params.get('Path to Pulse')
        if not pulse_src:
            raise ValueError(
                f"{self.name}: 'Path to Pulse' is required (refocusing waveform).")
        # Stage the waveform inside the (mounted) workdir — the picked file may
        # live outside the tree the Docker container can see.
        pulse_path = self._make_relative(self._stage_into_workdir(pulse_src))

        return [
            float(params['Samples']),
            float(params['Bandwidth']),
            float(params['Bfield']),
            float(params.get('Linewidth') or 1.0),
            float(tau1), float(tau2),
            pulse_path,
            float(params.get('RefTp') or 5.0),
            float(params.get('thkX') or 2.0),
            float(params.get('thkY') or 2.0),
            float(params.get('fovX') or 3.0),
            float(params.get('fovY') or 3.0),
            int(float(params.get('nX') or 8)),
            int(float(params.get('nY') or 8)),
            float(params.get('Flip Angle') or 180.0),
            # rotating-frame centre of the shaped-pulse simulation, in ppm
            float(params.get('Sim Centre (ppm)') or 4.65),
        ]


# =================================================================== Shaped
# (A future not-yet-implemented backend can set `_is_stub = True` on its
# class to get the GUI's "(in development)" treatment and the guarded
# NotImplementedError in run_simulation.)

class FidaSemiLaserShaped(FidaBackend):
    """semi-LASER (Öz 2018) with one shaped AFP waveform for both refocusing
    pairs on a spatial grid. 'Standard' uses FID-A's fast method (Zhang 2017:
    the X pair over the x grid, the Y pair over the y grid, nX + nY runs,
    coherence-order filtered — r = 0.99 against the phase cycle at a fifth of
    the cost); 'Phase cycled' runs FID-A's explicit 4-step cycle on the full
    nX × nY grid. FID-A's own `sim_semiLASER_shaped` is not used: it ignores
    its centreFreq / flipAngle arguments (`nargin < 18` in a 14-argument
    function)."""

    _kind = 'semilaser_shaped'

    def __init__(self):
        super().__init__()
        self.name, self.display_name = 'FidaSemiLaserShaped', 'semi-LASER shaped'
        self.modes = ['Standard', 'Phase cycled']
        self.current_mode = 'Standard'
        self.file_selection = ['Path to Pulse']
        self.mandatory_params = _shaped_params()
        self._refresh_metab_list()

    def parseProtocol(self, protocol):
        if protocol is None:
            return None
        p = str(protocol).lower()
        return 'sLASER' if ('slaser' in p or 'semi' in p) else None

    def _stage_pulse(self, params):
        pulse_src = params.get('Path to Pulse')
        if not pulse_src:
            raise ValueError(
                f"{self.name}: 'Path to Pulse' is required (AFP waveform, "
                f"e.g. a GOIA pulse).")
        return self._make_relative(self._stage_into_workdir(pulse_src))

    def active_kind(self):
        # 'Phase cycled': FID-A's 4-step refocusing phase cycle
        # (run_simSemiLASERShaped_phCyc.m) - 4x the simulations.
        return ('semilaser_shaped_phcyc' if self.current_mode == 'Phase cycled'
                else 'semilaser_shaped')

    def _build_args(self, params, metab):
        return [
            float(params['Samples']),
            float(params['Bandwidth']),
            float(params['Bfield']),
            float(params.get('Linewidth') or 1.0),
            float(params['TE']),
            self._stage_pulse(params),
            float(params.get('RefTp') or 5.0),
            float(params.get('thkX') or 2.0),
            float(params.get('thkY') or 2.0),
            float(params.get('fovX') or 3.0),
            float(params.get('fovY') or 3.0),
            int(float(params.get('nX') or 8)),
            int(float(params.get('nY') or 8)),
            float(params.get('Flip Angle') or 180.0),
            float(params.get('Sim Centre (ppm)') or 4.65),
        ]


class FidaSteamShaped(FidaBackend):
    """sim_steam_shaped: STEAM with shaped 90° pulses + spatial grid."""

    _kind = 'steam_shaped'

    def __init__(self):
        super().__init__()
        self.name, self.display_name = 'FidaSteamShaped', 'STEAM shaped'
        self.file_selection = ['Path to Pulse']
        self.mandatory_params = _shaped_params({'TM': 10.0})
        # STEAM pulses are 90° excitations, not 180° refocusers
        self.mandatory_params['Flip Angle'] = 90.0
        self._refresh_metab_list()

    def parseProtocol(self, protocol):
        if protocol is None:
            return None
        return 'STEAM' if 'steam' in str(protocol).lower() else None

    def _build_args(self, params, metab):
        pulse_src = params.get('Path to Pulse')
        if not pulse_src:
            raise ValueError(
                f"{self.name}: 'Path to Pulse' is required (excitation waveform).")
        pulse_path = self._make_relative(self._stage_into_workdir(pulse_src))
        return [
            float(params['Samples']),
            float(params['Bandwidth']),
            float(params['Bfield']),
            float(params.get('Linewidth') or 1.0),
            float(params['TE']),
            float(params.get('TM') or 10.0),
            pulse_path,
            float(params.get('RefTp') or 5.0),
            float(params.get('thkX') or 2.0),
            float(params.get('thkY') or 2.0),
            float(params.get('fovX') or 3.0),
            float(params.get('fovY') or 3.0),
            int(float(params.get('nX') or 8)),
            int(float(params.get('nY') or 8)),
            float(params.get('Flip Angle') or 90.0),
            float(params.get('Sim Centre (ppm)') or 4.65),
        ]


class FidaSpinEchoShaped(FidaBackend):
    """sim_spinecho_shaped: 1-D shaped refocusing with a subtractive
    [0°, 90°] phase cycle (per FID-A's run_simSpinEchoShaped)."""

    _kind = 'spinecho_shaped'

    def __init__(self):
        super().__init__()
        self.name, self.display_name = 'FidaSpinEchoShaped', 'Spin Echo shaped'
        self.file_selection = ['Path to Pulse']
        self.mandatory_params = _shaped_params()
        for k in ('thkY', 'fovY', 'nY', 'Flip Angle'):
            self.mandatory_params.pop(k, None)
        self._refresh_metab_list()

    def _build_args(self, params, metab):
        pulse_src = params.get('Path to Pulse')
        if not pulse_src:
            raise ValueError(
                f"{self.name}: 'Path to Pulse' is required (refocusing waveform).")
        pulse_path = self._make_relative(self._stage_into_workdir(pulse_src))
        return [
            float(params['Samples']),
            float(params['Bandwidth']),
            float(params['Bfield']),
            float(params.get('Linewidth') or 1.0),
            float(params['TE']),
            pulse_path,
            float(params.get('RefTp') or 5.0),
            float(params.get('thkX') or 2.0),
            float(params.get('fovX') or 3.0),
            int(float(params.get('nX') or 8)),
        ]


class FidaMegaPressShaped(FidaBackend):
    """MEGA-PRESS with real editing and/or refocusing waveforms; each
    metabolite yields '(ON)', '(OFF)', '(DIFF)' entries. Modes: edit-only
    shaped (sim_megapress_shapedEdit, ideal refocusing), fully shaped
    (sim_megapress_shaped, both waveforms on a spatial grid) and
    refocusing-only shaped (sim_megapress_shapedRefoc, ideal editing).
    """

    _kind = 'megapress_shapededit'

    # run_simMegaPressShapedEdit.m timing at TE = 68 ms:
    # excite – t1 – 180 – t2 – edit – t3 – 180 – t4 – edit – t5 – ADC
    _TE68_TAUS = (5.0, 17.0, 17.0, 17.0, 12.0)

    def __init__(self):
        super().__init__()
        self.name, self.display_name = 'FidaMegaPressShaped', 'MEGA-PRESS shaped'
        self.modes = ['Edit-only shaped (ideal refoc)',
                      'Full shaped (refoc + edit)',
                      'Refoc-only shaped (ideal edit)']
        self.current_mode = 'Edit-only shaped (ideal refoc)'
        self.file_selection = ['Edit Pulse Path', 'Path to Pulse']
        self.mandatory_params = {
            'Samples':   None, 'Bandwidth': None, 'Bfield': None,
            'Linewidth': 1.0,  'TE':        None,   # 68 ms is the standard
            'Edit Pulse Path': None,
            'Edit Tp':         20.0,
            'Edit On':         1.9,                 # ppm (GABA); 4.56 for GSH
            'Edit Off':        7.5,
            'Edit Bandwidth (ppm)': 1.0,            # ideal-edit modes only
            'Sim Centre (ppm)': 4.65,
            # shaped-refocusing modes only (spatial grid kept small: the
            # fully shaped mode runs 32 simulations per grid point)
            'Path to Pulse': None,
            'RefTp':         5.0,
            'thkX': 2.0, 'thkY': 2.0, 'fovX': 3.0, 'fovY': 3.0,
            'nX': 4, 'nY': 4,
            'Metabolites': [],
        }
        self._refresh_metab_list()

    _MODE_KINDS = {
        'Edit-only shaped (ideal refoc)': 'megapress_shapededit',
        'Full shaped (refoc + edit)':     'megapress_shaped',
        'Refoc-only shaped (ideal edit)': 'megapress_shapedrefoc',
    }
    _GRID_KEYS = ('Path to Pulse', 'RefTp', 'thkX', 'thkY', 'fovX', 'fovY',
                  'nX', 'nY')

    def active_kind(self):
        return self._MODE_KINDS[self.current_mode]

    def get_params_for_mode(self, mode=None):
        mode = mode or self.current_mode
        params = dict(self.mandatory_params)
        if mode == 'Edit-only shaped (ideal refoc)':
            for k in self._GRID_KEYS + ('Edit Bandwidth (ppm)',):
                params.pop(k, None)
        elif mode == 'Full shaped (refoc + edit)':
            params.pop('Edit Bandwidth (ppm)', None)
            # sim_megapress_shaped simulates in a fixed 3 ppm frame
            # (fida_run.m re-references the result) — no centre to choose
            params.pop('Sim Centre (ppm)', None)
        else:  # ideal editing: no editing waveform
            for k in ('Edit Pulse Path', 'Edit Tp', 'Edit Off', 'Sim Centre (ppm)'):
                params.pop(k, None)
        return params

    def parseProtocol(self, protocol):
        if protocol is None:
            return None
        return 'MEGA-PRESS' if 'mega' in str(protocol).lower() else None

    def _staged(self, params, key, what):
        src = params.get(key)
        if not src:
            raise ValueError(f"{self.name}: '{key}' is required ({what}).")
        return self._make_relative(self._stage_into_workdir(src))

    def _grid_args(self, params):
        return [
            self._staged(params, 'Path to Pulse', 'refocusing waveform'),
            float(params.get('RefTp') or 5.0),
            float(params.get('thkX') or 2.0),
            float(params.get('thkY') or 2.0),
            float(params.get('fovX') or 3.0),
            float(params.get('fovY') or 3.0),
            int(float(params.get('nX') or 4)),
            int(float(params.get('nY') or 4)),
        ]

    def run_simulation(self, params, progress_callback=None, stop_event=None):
        if self.octave is None:
            print("Initializing Octave runtime...")
            self.initialize_octave(prefer_docker=True)
        self.setup_octave_paths()
        self.ensure_workdir()

        te = float(params['TE'])
        scale = te / 68.0
        taus = [t * scale for t in self._TE68_TAUS]
        head = [
            float(params['Samples']),
            float(params['Bandwidth']),
            float(params['Bfield']),
            float(params.get('Linewidth') or 1.0),
            *taus,
        ]
        mode = self.current_mode
        if mode == 'Edit-only shaped (ideal refoc)':
            base_args = head + [
                self._staged(params, 'Edit Pulse Path', 'editing waveform'),
                float(params.get('Edit Tp') or 20.0),
                float(params.get('Edit On') or 1.9),
                float(params.get('Edit Off') or 7.5),
                float(params.get('Sim Centre (ppm)') or 4.65),
            ]
        elif mode == 'Full shaped (refoc + edit)':
            base_args = head + [
                self._staged(params, 'Edit Pulse Path', 'editing waveform'),
                float(params.get('Edit Tp') or 20.0),
                float(params.get('Edit On') or 1.9),
                float(params.get('Edit Off') or 7.5),
                *self._grid_args(params),
                4.65,   # ignored by fida_run.m: sim_megapress_shaped uses a 3 ppm frame
            ]
        else:  # 'Refoc-only shaped (ideal edit)'
            base_args = head + [
                float(params.get('Edit On') or 1.9),
                float(params.get('Edit Bandwidth (ppm)') or 1.0),
                *self._grid_args(params),
            ]
        return self._run_on_off_subspectra(params, base_args,
                                           progress_callback, stop_event)


class FidaMegaSpecialShaped(FidaBackend):
    """sim_megaspecial_shaped: 1-D-localised MEGA-SPECIAL with a shaped,
    frequency-shifted editing pulse and one shaped refocusing pulse.
    Each metabolite yields '(ON)', '(OFF)', '(DIFF)' entries."""

    _kind = 'megaspecial_shaped'

    # run_simMegaSpecialShaped.m timing at TE = 68 ms:
    # excite – t1 – edit – t2 – 180 – t3 – edit – t4 – ADC
    _TE68_TAUS = (17.0, 17.0, 17.0, 17.0)

    def __init__(self):
        super().__init__()
        self.name, self.display_name = 'FidaMegaSpecialShaped', 'MEGA-SPECIAL shaped'
        self.file_selection = ['Path to Pulse', 'Edit Pulse Path']
        self.mandatory_params = {
            'Samples':   None, 'Bandwidth': None, 'Bfield': None,
            'Linewidth': 1.0,  'TE':        None,   # 68 ms is the standard
            'Path to Pulse':   None,                # refocusing waveform
            'RefTp':           5.0,
            'Edit Pulse Path': None,
            'Edit Tp':         14.0,
            'Edit On':         1.9,                 # ppm (GABA); 4.56 for GSH
            'Edit Off':        7.5,
            'thkX': 2.0, 'fovX': 3.0, 'nX': 8,
            'Metabolites': [],
        }
        self._refresh_metab_list()

    def parseProtocol(self, protocol):
        if protocol is None:
            return None
        p = str(protocol).lower()
        return 'MEGA-SPECIAL' if ('mega' in p and 'special' in p) else None

    def run_simulation(self, params, progress_callback=None, stop_event=None):
        if self.octave is None:
            print("Initializing Octave runtime...")
            self.initialize_octave(prefer_docker=True)
        self.setup_octave_paths()
        self.ensure_workdir()

        refoc_src = params.get('Path to Pulse')
        edit_src = params.get('Edit Pulse Path')
        if not refoc_src or not edit_src:
            raise ValueError(
                f"{self.name}: 'Path to Pulse' (refocusing) and "
                f"'Edit Pulse Path' (editing) waveforms are both required.")
        refoc_path = self._make_relative(self._stage_into_workdir(refoc_src))
        edit_path = self._make_relative(self._stage_into_workdir(edit_src))

        te = float(params['TE'])
        scale = te / 68.0
        taus = [t * scale for t in self._TE68_TAUS]
        base_args = [
            float(params['Samples']),
            float(params['Bandwidth']),
            float(params['Bfield']),
            float(params.get('Linewidth') or 1.0),
            *taus,
            edit_path,
            float(params.get('Edit Tp') or 14.0),
            float(params.get('Edit On') or 1.9),
            float(params.get('Edit Off') or 7.5),
            refoc_path,
            float(params.get('RefTp') or 5.0),
            float(params.get('thkX') or 2.0),
            float(params.get('fovX') or 3.0),
            int(float(params.get('nX') or 8)),
            4.65,   # ignored by fida_run.m: sim_megaspecial_shaped uses a 3 ppm frame
        ]
        return self._run_on_off_subspectra(params, base_args,
                                           progress_callback, stop_event)


class FidaLaser(FidaBackend):
    """sim_laser: ideal-AFP LASER with six equally spaced echoes."""

    _kind = 'laser'

    def __init__(self):
        super().__init__()
        self.name, self.display_name = 'FidaLaser', 'LASER (ideal AFP)'
        self.mandatory_params = {
            'Samples':   None, 'Bandwidth': None, 'Bfield': None,
            'Linewidth': 1.0,  'TE':        None,
            'Metabolites': [],
        }
        self._refresh_metab_list()

    def parseProtocol(self, protocol):
        if protocol is None:
            return None
        return 'LASER' if 'laser' in str(protocol).lower() else None

    def _build_args(self, params, metab):
        return [
            float(params['Samples']),
            float(params['Bandwidth']),
            float(params['Bfield']),
            float(params.get('Linewidth') or 1.0),
            float(params['TE']),
        ]


class FidaMegaPressIdeal(FidaBackend):
    """sim_megapress: MEGA-PRESS with instantaneous localization and editing.

    The ideal editing pulse inverts every spin within the editing band around
    'Edit On'; the OFF sub-spectrum applies no editing. Each metabolite yields
    three basis entries: '<metab> (ON)', '<metab> (OFF)', '<metab> (DIFF)'
    (DIFF = ON − OFF, the edited difference basis).
    """

    _kind = 'megapress_ideal'

    # Siemens MEGA-PRESS timing at TE = 68 ms:
    # 90 – t1 – 180 – t2 – edit – t3 – 180 – t4 – edit – t5 – ADC.
    # Other echo times scale this scheme proportionally (TE/68).
    _TE68_TAUS = (4.545, 12.7025, 21.7975, 12.7025, 17.2526)

    def __init__(self):
        super().__init__()
        self.name, self.display_name = 'FidaMegaPressIdeal', 'MEGA-PRESS ideal'
        self.mandatory_params = {
            'Samples':   None, 'Bandwidth': None, 'Bfield': None,
            'Linewidth': 1.0,  'TE':        None,   # 68 ms is the standard
            'Edit On':   1.9,                       # ppm (GABA); 4.56 for GSH
            'Edit Bandwidth (ppm)': 1.0,
            'Metabolites': [],
        }
        self._refresh_metab_list()

    def parseProtocol(self, protocol):
        if protocol is None:
            return None
        return 'MEGA-PRESS' if 'mega' in str(protocol).lower() else None

    def run_simulation(self, params, progress_callback=None, stop_event=None):
        if self.octave is None:
            print("Initializing Octave runtime...")
            self.initialize_octave(prefer_docker=True)
        self.setup_octave_paths()
        self.ensure_workdir()

        te = float(params['TE'])
        scale = te / 68.0
        taus = [t * scale for t in self._TE68_TAUS]
        base_args = [
            float(params['Samples']),
            float(params['Bandwidth']),
            float(params['Bfield']),
            float(params.get('Linewidth') or 1.0),
            *taus,
            float(params.get('Edit On') or 1.9),
            float(params.get('Edit Bandwidth (ppm)') or 1.0),
        ]
        return self._run_on_off_subspectra(params, base_args,
                                           progress_callback, stop_event)


class FidaSpinEchoXN(FidaBackend):
    """sim_spinecho_xN: multi-echo spin-echo train (CPMG-style)."""

    _kind = 'spinecho_xn'

    def __init__(self):
        super().__init__()
        self.name, self.display_name = 'FidaSpinEchoXN', 'Spin Echo (multi-echo)'
        self.mandatory_params = {
            'Samples':   None, 'Bandwidth': None, 'Bfield': None,
            'Linewidth': 1.0,  'Tau':       15.0, 'Nechoes': 2,
            'Metabolites': [],
        }
        self._refresh_metab_list()

    def _build_args(self, params, metab):
        return [
            float(params['Samples']),
            float(params['Bandwidth']),
            float(params['Bfield']),
            float(params.get('Linewidth') or 1.0),
            float(params.get('Tau') or 15.0),
            int(float(params.get('Nechoes') or 2)),
        ]


class FidaOnePulse(FidaBackend):
    """Pulse-acquire FID. Modes: Ideal (sim_onepulse), Shaped
    (sim_onepulse_shaped, frequency-selective excitation waveform), Delay
    (sim_onepulse_delay, ADC onset delay) and Arbitrary phase
    (sim_onepulse_arbPh)."""

    _kind = 'onepulse'

    def __init__(self):
        super().__init__()
        self.name, self.display_name = 'FidaOnePulse', 'One pulse (FID only)'
        self.modes = ['Ideal', 'Shaped', 'Delay', 'Arbitrary phase']
        self.current_mode = 'Ideal'
        self.file_selection = ['Path to Pulse']
        self.mandatory_params = {
            'Samples':   None, 'Bandwidth': None, 'Bfield': None,
            'Linewidth': 1.0,
            'Flip Angle': 90.0,       # Shaped only
            'Path to Pulse': None,    # Shaped only
            'RefTp': 5.0,             # Shaped only: pulse duration [ms]
            'Delay': 0.5,             # Delay only: ADC onset delay [ms]
            'Pulse Phase': 0.0,       # Arbitrary phase only [deg]
            'Metabolites': [],
        }
        self._refresh_metab_list()

    _MODE_KINDS = {
        'Ideal':           'onepulse',
        'Shaped':          'onepulse_shaped',
        'Delay':           'onepulse_delay',
        'Arbitrary phase': 'onepulse_arbph',
    }
    _MODE_KEYS = {
        'Ideal':           (),
        'Shaped':          ('Flip Angle', 'Path to Pulse', 'RefTp'),
        'Delay':           ('Delay',),
        'Arbitrary phase': ('Pulse Phase',),
    }

    def active_kind(self):
        return self._MODE_KINDS[self.current_mode]

    def get_params_for_mode(self, mode=None):
        mode = mode or self.current_mode
        keep = set(self._MODE_KEYS[mode])
        all_mode_keys = {k for keys in self._MODE_KEYS.values() for k in keys}
        return {k: v for k, v in self.mandatory_params.items()
                if k not in all_mode_keys or k in keep}

    def _build_args(self, params, metab):
        head = [
            float(params['Samples']),
            float(params['Bandwidth']),
            float(params['Bfield']),
            float(params.get('Linewidth') or 1.0),
        ]
        mode = self.current_mode
        if mode == 'Ideal':
            return head
        if mode == 'Delay':
            return head + [float(params.get('Delay') or 0.0)]
        if mode == 'Arbitrary phase':
            return head + [float(params.get('Pulse Phase') or 0.0)]
        pulse_src = params.get('Path to Pulse')
        if not pulse_src:
            raise ValueError(
                f"{self.name}: 'Path to Pulse' is required (excitation waveform).")
        return head + [
            self._make_relative(self._stage_into_workdir(pulse_src)),
            float(params.get('RefTp') or 5.0),
            float(params.get('Flip Angle') or 90.0),
        ]


# =================================================================== registry
# GUI dropdown order: Ideal first (most-used), then shaped variants, then
# the niche / debug entries.
FIDA_BACKENDS = [
    FidaIdeal,
    FidaPressShaped,
    FidaSemiLaserShaped,
    FidaSteamShaped,
    FidaSpinEchoShaped,
    FidaMegaPressShaped,
    FidaMegaSpecialShaped,
    FidaLaser,
    FidaMegaPressIdeal,
    FidaSpinEchoXN,
    FidaOnePulse,
]


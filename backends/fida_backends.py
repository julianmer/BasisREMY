####################################################################################################
#                                  backends/fida/fida_backends.py                                   #
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
#          (adapters/backends/fida/fida_run.m) differ between subclasses. FidaIdeal is the         #
#          ex-"LCModel" backend, now living natively under FID-A; backends/lcmodel_backend.py      #
#          stays as a thin alias for backwards compat.                                             #
#                                                                                                  #
####################################################################################################


from __future__ import annotations

import os

import numpy as np

from backends.base import Backend


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
        'Center Freq':   4.65,
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
    branch is added inside ``adapters/backends/fida/fida_run.m``.
    """

    _kind: str = ''      # dispatch key for fida_run.m
    _is_stub: bool = False

    def __init__(self):
        super().__init__()
        self.category = 'FID-A'
        self.requires_octave = True
        self.metabs = dict(_DEFAULT_FIDA_METABS)
        self.optional_params = {'Nucleus': None, 'TR': None}

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
        self.octave.addpath('./adapters/backends/')
        self.octave.addpath('./adapters/backends/fida/')

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

    # -------------------------------------------------- per-subclass hook
    def _build_args(self, params, metab):
        """Positional args (AFTER metab + kind) for ``fida_run.m``."""
        raise NotImplementedError

    # -------------------------------------------------- driver
    def run_simulation(self, params, progress_callback=None):
        if self._is_stub or not self._kind:
            raise NotImplementedError(
                f"{self.name}: this FID-A wrapper is a stub. The schema is "
                "complete but the matching branch in adapters/backends/fida/"
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
            extra_args = self._build_args(params, metab)
            results = self.octave.feval(
                'fida_run', metab, self._kind, *extra_args, nout=5,
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
            'Metabolites': [],
            'Center Freq': None,
        }
        self._refresh_metab_list()

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
        if 'press' in p:    return 'PRESS'
        if 'steam' in p:    return 'STEAM'
        if 'spin' in p or 'se' in p: return 'Spin Echo'
        if 'laser' in p:    return 'LASER'
        return None

    # ---- args ---------------------------------------------------------
    @staticmethod
    def _seq_to_fida(seq):
        return {'Spin Echo': 'se', 'PRESS': 'p',
                'STEAM': 'st', 'LASER': 'l'}.get(seq, seq)

    def _build_args(self, params, metab):
        # workdir for the FID-A-side .RAW writes (kept for parity with the
        # original sim_lcmrawbasis flow; the adapter ignores it but it
        # keeps the path structure consistent across runs).
        out = self._make_relative(self.ensure_workdir()) + os.sep
        return [
            float(params['Samples']),
            float(params['Bandwidth']),
            float(params['Bfield']),
            float(params.get('Linewidth') or 1),
            float(params['TE']),
            float(params.get('TE2') or 0),
            self._seq_to_fida(params['Sequence']),
            out,
        ]


# =================================================================== PRESS shaped
class FidaPressShaped(FidaBackend):
    """sim_press_shaped: PRESS with shaped refocusing pulses + spatial grid."""

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

        pulse_path = self._make_relative(params.get('Path to Pulse'))
        if not pulse_path:
            raise ValueError(
                f"{self.name}: 'Path to Pulse' is required (refocusing waveform).")

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
            # centreFreq for sim_press_shaped is in ppm (the rotating-frame
            # centre). REMY fills 'Center Freq' from scanner metadata in Hz
            # (e.g. 127736713). When that's the case (>1000), use the standard
            # water reference 4.65 ppm so FID-A's spin shifts and the GUI's
            # ppm axis (which adds +4.65 offset) are consistent.
            (float(params.get('Center Freq') or 4.65)
             if (float(params.get('Center Freq') or 4.65) <= 1000)
             else 4.65),
        ]


# =================================================================== Stubs
# Each stub carries a complete parameter schema for the GUI but has no Octave
# dispatch yet. They will become real backends as the matching branches are
# added to fida_run.m. `_is_stub = True` keeps the schema-only tests honest.

class _Stub(FidaBackend):
    _is_stub = True


class FidaSemiLaserShaped(_Stub):
    """sim_semiLASER_shaped / _phCyc."""
    def __init__(self):
        super().__init__()
        self.name, self.display_name = 'FidaSemiLaserShaped', 'semi-LASER shaped'
        self.modes = ['Standard', 'Phase cycled']
        self.current_mode = 'Standard'
        self.file_selection = ['Path to Pulse']
        self.mandatory_params = _shaped_params()
        self._refresh_metab_list()


class FidaSteamShaped(_Stub):
    """sim_steam_shaped."""
    def __init__(self):
        super().__init__()
        self.name, self.display_name = 'FidaSteamShaped', 'STEAM shaped'
        self.file_selection = ['Path to Pulse']
        self.mandatory_params = _shaped_params({'TM': 10.0})
        self._refresh_metab_list()


class FidaSpinEchoShaped(_Stub):
    """sim_spinecho_shaped (1-D refocusing only)."""
    def __init__(self):
        super().__init__()
        self.name, self.display_name = 'FidaSpinEchoShaped', 'Spin Echo shaped'
        self.file_selection = ['Path to Pulse']
        self.mandatory_params = _shaped_params()
        for k in ('thkY', 'fovY', 'nY'):
            self.mandatory_params.pop(k, None)
        self._refresh_metab_list()


class FidaMegaPressShaped(_Stub):
    """sim_megapress_shaped / _shapedEdit / _shapedRefoc."""
    def __init__(self):
        super().__init__()
        self.name, self.display_name = 'FidaMegaPressShaped', 'MEGA-PRESS shaped'
        self.modes = ['Full shaped (refoc + edit)',
                      'Edit-only shaped (ideal refoc)',
                      'Refoc-only shaped (ideal edit)']
        self.current_mode = 'Full shaped (refoc + edit)'
        self.file_selection = ['Path to Pulse', 'Edit Pulse Path']
        self.mandatory_params = _shaped_params({
            'Edit Pulse Path': None,
            'Edit Tp':         20.0,
            'Edit On':         1.9,
            'Edit Off':        7.5,
            'Edit Target':     'GABA',
        })
        self.dropdown = {'Edit Target': ['GABA', 'GSH', 'Lac', 'PE']}
        self._refresh_metab_list()


class FidaMegaSpecialShaped(_Stub):
    """sim_megaspecial_shaped (1-D-localised MEGA)."""
    def __init__(self):
        super().__init__()
        self.name, self.display_name = 'FidaMegaSpecialShaped', 'MEGA-SPECIAL shaped'
        self.file_selection = ['Path to Pulse', 'Edit Pulse Path']
        self.mandatory_params = _shaped_params({
            'Edit Pulse Path': None,
            'Edit Tp':         20.0,
            'Edit On':         1.9,
            'Edit Off':        7.5,
            'Edit Target':     'GABA',
        })
        for k in ('thkY', 'fovY', 'nY'):
            self.mandatory_params.pop(k, None)
        self.dropdown = {'Edit Target': ['GABA', 'GSH', 'Lac', 'PE']}
        self._refresh_metab_list()


class FidaLaser(_Stub):
    """sim_laser (ideal AFP)."""
    def __init__(self):
        super().__init__()
        self.name, self.display_name = 'FidaLaser', 'LASER (ideal AFP)'
        self.mandatory_params = {
            'Samples':   None, 'Bandwidth': None, 'Bfield': None,
            'Linewidth': 1.0,  'TE':        None,
            'Center Freq': 2.3,
            'Metabolites': [],
        }
        self._refresh_metab_list()


class FidaMegaPressIdeal(_Stub):
    """sim_megapress (ideal pulses)."""
    def __init__(self):
        super().__init__()
        self.name, self.display_name = 'FidaMegaPressIdeal', 'MEGA-PRESS ideal'
        self.mandatory_params = {
            'Samples':   None, 'Bandwidth': None, 'Bfield': None,
            'Linewidth': 1.0,  'TE':        None,
            'Edit On':   1.9,  'Edit Off':  7.5,
            'Edit Target': 'GABA',
            'Center Freq': 2.3,
            'Metabolites': [],
        }
        self.dropdown = {'Edit Target': ['GABA', 'GSH', 'Lac', 'PE']}
        self._refresh_metab_list()


class FidaSpinEchoXN(_Stub):
    """sim_spinecho_xN (multi-echo train)."""
    def __init__(self):
        super().__init__()
        self.name, self.display_name = 'FidaSpinEchoXN', 'Spin Echo (multi-echo)'
        self.mandatory_params = {
            'Samples':   None, 'Bandwidth': None, 'Bfield': None,
            'Linewidth': 1.0,  'Tau':       15.0, 'Nechoes': 2,
            'Center Freq': 2.3,
            'Metabolites': [],
        }
        self._refresh_metab_list()


class FidaOnePulse(_Stub):
    """sim_onepulse / _shaped / _delay / _arbPh (FID only)."""
    def __init__(self):
        super().__init__()
        self.name, self.display_name = 'FidaOnePulse', 'One pulse (FID only)'
        self.modes = ['Ideal', 'Shaped', 'Delay', 'Arbitrary phase']
        self.current_mode = 'Ideal'
        self.file_selection = ['Path to Pulse']
        self.mandatory_params = {
            'Samples':   None, 'Bandwidth': None, 'Bfield': None,
            'Linewidth': 1.0,
            'Flip Angle': 90.0,
            'Path to Pulse': None,
            'Center Freq': 2.3,
            'Metabolites': [],
        }
        self._refresh_metab_list()


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


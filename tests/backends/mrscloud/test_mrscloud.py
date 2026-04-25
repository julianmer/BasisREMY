####################################################################################################
#                                      test_mrscloud.py                                            #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 24/04/26                                                                                #
#                                                                                                  #
# Purpose: Test MRSCloudBackend (rewritten 24/04/26 to drive the real MRSCloud workflow            #
#          through `adapters/backends/mrscloud_run_metab.m` instead of FID-A's                     #
#          `sim_lcmrawbasis`).                                                                     #
#                                                                                                  #
#          Tests are split into:                                                                   #
#            - Pure-Python checks (always run): metabolite list, REMY parsing, params              #
#              structure, no-Output-Path leakage.                                                  #
#            - Live Octave checks (skipped if no Octave / Docker available):                       #
#              one parametrize per (sequence × localization) combo. These require                  #
#              a working MRSCloud + FID-A install.                                                 #
#                                                                                                  #
#          TODO: get vendor-specific testing data for HERMES / HERCULES / MEGA-sLASER and          #
#          extend the integration suite to cover them with real REMY-extracted parameters.         #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

from core.basisremy import BasisREMY
from backends.mrscloud_backend import MRSCloudBackend


# ----------------------------------------------------------------- helpers
def _octave_available() -> bool:
    """True if either a local oct2py-compatible Octave or the Docker fallback works."""
    try:
        from core.octave_manager import OctaveManager
        m = OctaveManager()
        return bool(m.check_local_octave_availability() or m.check_docker_availability())
    except Exception:
        return False


_OCTAVE = pytest.mark.skipif(
    not _octave_available(),
    reason="Octave (local or Docker) not available — skipping live MRSCloud simulations.",
)


# =================================================================== pure
class TestMRSCloudParameters:
    """Pure-Python checks that never require Octave."""

    @pytest.fixture
    def backend(self):
        return MRSCloudBackend()

    def test_initialization(self, backend):
        assert backend.name == 'MRSCloud'
        assert backend.requires_octave is True
        assert backend.octave is None  # lazy

    def test_in_basisremy_registry(self):
        br = BasisREMY()
        assert 'MRSCloud' in br.available_backends

    # ---- export-related leakage --------------------------------------
    def test_no_output_path_in_params(self, backend):
        """Output Path / Format / Make .raw / Add Ref. must NOT be exposed any more."""
        for forbidden in ('Output Path', 'Output Format', 'Make .raw', 'Add Ref.'):
            assert forbidden not in backend.mandatory_params, \
                f"{forbidden!r} leaked back into mandatory_params"
            assert forbidden not in backend.dropdown, \
                f"{forbidden!r} leaked back into dropdown"
            assert forbidden not in backend.optional_params, \
                f"{forbidden!r} leaked back into optional_params"

    # ---- mandatory params --------------------------------------------
    # NOTE: Bfield, Center Freq, Linewidth, TE2 were intentionally removed
    # in the 24/04/26 refactor — MRSCloud derives them internally from
    # `Field Strength` + vendor (see externals/mrscloud/functions/load_parameters.m).
    # Editing fields (Edit Target / On / Off / Tp) are spliced in by
    # get_params_for_mode() only when an editing sequence is selected, so
    # they are NOT in the default mandatory_params either.
    @pytest.mark.parametrize("key", [
        'System', 'Sequence', 'Localization', 'Field Strength',
        'Samples', 'Bandwidth', 'TE',
        'Spatial Points', 'Metabolites',
    ])
    def test_mandatory_keys_present(self, backend, key):
        assert key in backend.mandatory_params

    @pytest.mark.parametrize("forbidden_key", [
        'Bfield', 'Center Freq', 'Linewidth', 'TE2',
    ])
    def test_removed_keys_absent(self, backend, forbidden_key):
        """These were dropped in the 24/04/26 refactor — make sure they don't sneak back."""
        assert forbidden_key not in backend.mandatory_params, \
            f"{forbidden_key!r} is no longer used by MRSCloud and must not be in mandatory_params"

    # ---- dropdowns ----------------------------------------------------
    def test_dropdown_options(self, backend):
        assert 'Philips' in backend.dropdown['System']
        assert 'GE' in backend.dropdown['System']
        assert set(backend.dropdown['Sequence']) >= {'UnEdited', 'MEGA', 'HERMES', 'HERCULES'}
        assert set(backend.dropdown['Localization']) >= {'PRESS', 'sLASER', 'STEAM_7T'}
        assert '3T' in backend.dropdown['Field Strength']

    # ---- metabolite library ------------------------------------------
    def test_metabolites_match_mrscloud_readme(self, backend):
        """Default list must include every common metabolite from MRSCloud's README."""
        common_healthy = {
            'Asc','Asp','Cr','GABA','GPC','GSH','Gln','Glu','Gly',
            'Lac','mI','NAA','NAAG','PCh','PCr','PE','sI','Tau',
        }
        for m in common_healthy:
            assert m in backend.metabs, f"missing common-brain metabolite {m!r}"
        # Specific-interest list must also be available (off by default is OK)
        for m in ['Ala','Cit','Tyros','Phenyl','bHB','bHG','EtOH']:
            assert m in backend.metabs, f"missing specific-interest metabolite {m!r}"

    def test_default_selection_nonempty(self, backend):
        defaults = [k for k, v in backend.metabs.items() if v]
        assert len(defaults) > 10
        for must_have in ('NAA', 'Cr', 'Glu', 'GABA'):
            assert must_have in defaults

    # ---- REMY parsing -------------------------------------------------
    def test_parse_remy_basic(self, backend):
        remy = {
            'NumberOfDatapoints': 4096,
            'SpectralWidth':      4000,
            'B0':                 3.0,
            'TE':                 35,
            'Protocol':           'PRESS_35',
            'Manufacturer':       'Philips',
        }
        m, _o = backend.parseREMY(remy)
        assert m['Samples']        == 4096
        assert m['Bandwidth']      == 4000
        assert m['TE']             == 35
        assert m['Sequence']       == 'UnEdited'
        assert m['Localization']   == 'PRESS'
        assert m['System']         == 'Philips'
        assert m['Field Strength'] == '3T'
        # Bfield / Center Freq are intentionally NOT emitted any more — MRSCloud
        # derives them internally from Field Strength + vendor.
        assert 'Bfield'      not in m
        assert 'Center Freq' not in m

    @pytest.mark.parametrize("protocol,seq,loc", [
        ('PRESS_35',          'UnEdited', 'PRESS'),
        ('STEAM_7T_11ms',     'UnEdited', 'STEAM_7T'),
        ('semi_LASER_TE30',   'UnEdited', 'sLASER'),
        ('MEGA_PRESS_GABA',   'MEGA',     'PRESS'),
        ('HERMES_GABA_GSH',   'HERMES',   'PRESS'),
        ('HERCULES_TE80',     'HERCULES', 'PRESS'),
        ('mega_sLASER_GABA',  'MEGA',     'sLASER'),
    ])
    def test_parse_protocol_table(self, backend, protocol, seq, loc):
        assert backend.parseProtocol(protocol)     == seq
        assert backend.parseLocalization(protocol) == loc

    @pytest.mark.parametrize("vendor_in,expected", [
        ('Philips Achieva',     'Philips'),
        ('SIEMENS',             'Siemens'),
        ('Universal_Philips',   'Universal_Philips'),
        ('Universal_Siemens',   'Universal_Siemens'),
        ('GE Healthcare',       'GE'),
        ('Bruker',              None),   # unsupported
    ])
    def test_parse_system(self, backend, vendor_in, expected):
        assert backend.parseSystem(vendor_in) == expected

    @pytest.mark.parametrize("b0,expected", [
        (1.5,  '1.5T'),
        (3.0,  '3T'),
        (7.0,  '7T'),
        (None, '3T'),
        ('',   '3T'),
    ])
    def test_field_str_from_b0(self, b0, expected):
        assert MRSCloudBackend._field_str_from_b0(b0) == expected

    # ---- workdir behaviour -------------------------------------------
    def test_workdir_lazy(self, backend):
        assert backend._workdir is None
        wd = backend.ensure_workdir()
        assert os.path.isdir(wd)
        assert backend.ensure_workdir() == wd
        backend.cleanup_workdir()
        assert backend._workdir is None


# ============================================================ live Octave
@_OCTAVE
class TestMRSCloudLive:
    """Per-(sequence, localization) live runs against real MRSCloud Octave code.

    These need:
      - Octave (Docker or local) installed.
      - externals/mrscloud submodule populated.
      - externals/mrscloud/functions/FID-A populated.

    NB: MRSCloud's README states "Product sequence and rf waveform are not
    shared in the GitHub repo." That means the *only* combo with fully-bundled
    pulses is `Universal_Philips + PRESS + UnEdited` (uses
    `pulses_universal/univ_eddenrefo.pta`). All other vendor / localization
    combinations are marked `xfail` with a TODO until the user installs the
    matching pulse waveform locally (see `MRSCloudBackend.required_pulse_files`).
    """

    @pytest.fixture(scope="class")
    def backend(self):
        return MRSCloudBackend()

    _COMMON = {
        'System':         'Universal_Philips',   # only fully-bundled vendor
        'Field Strength': '3T',
        'Samples':        2048,
        'Bandwidth':      2000,
        'TE':             35,
        'Spatial Points': 41,        # use the cheaper grid in tests
        'Metabolites':    ['NAA'],
    }

    @pytest.mark.parametrize("sequence,localization,edit_target", [
        # Fully bundled — should pass
        ('UnEdited', 'PRESS',    ''),
        # Need vendor product pulses not shipped by MRSCloud — TODO
        pytest.param('UnEdited', 'sLASER',   '',
                     marks=pytest.mark.xfail(strict=False,
                         reason="TODO get Philips_GOIA_WURST_100pts.mat (not in MRSCloud public repo)")),
        pytest.param('UnEdited', 'STEAM_7T', '',
                     marks=pytest.mark.xfail(strict=False,
                         reason="TODO get STEAM_7T pulse waveforms (not in MRSCloud public repo)")),
        pytest.param('MEGA',     'PRESS',  'GABA',
                     marks=pytest.mark.xfail(strict=False,
                         reason="TODO get testing data + edit-pulse waveforms for MEGA-PRESS")),
        pytest.param('HERMES',   'PRESS',  'GABA',
                     marks=pytest.mark.xfail(strict=False,
                         reason="TODO get testing data + edit-pulse waveforms for HERMES")),
        pytest.param('HERCULES', 'PRESS',  'GABA',
                     marks=pytest.mark.xfail(strict=False,
                         reason="TODO get testing data + edit-pulse waveforms for HERCULES")),
        pytest.param('MEGA',     'sLASER', 'GABA',
                     marks=pytest.mark.xfail(strict=False,
                         reason="TODO get testing data + GOIA pulse for MEGA-sLASER")),
    ])
    def test_run_one_metab(self, backend, sequence, localization, edit_target):
        # Pre-flight: skip cleanly if MRSCloud needs a pulse file that isn't bundled.
        missing = backend.missing_pulse_files(self._COMMON['System'], sequence, localization)
        if missing:
            pytest.skip(f"MRSCloud pulse files not available locally: {missing}")

        params = {**self._COMMON,
                  'Sequence': sequence,
                  'Localization': localization,
                  'Edit Target': edit_target}
        result = backend.run_simulation(params)
        assert isinstance(result, dict) and 'NAA' in result
        fid = result['NAA']
        assert isinstance(fid, np.ndarray)
        assert fid.dtype.kind == 'c'
        assert fid.size > 0
        assert np.max(np.abs(fid)) > 0, "FID is identically zero"

    def test_run_multiple_metabs(self, backend):
        # Only run the fully-bundled combo here.
        missing = backend.missing_pulse_files('Universal_Philips', 'UnEdited', 'PRESS')
        if missing:
            pytest.skip(f"MRSCloud pulse files not available locally: {missing}")

        params = {**self._COMMON,
                  'Sequence': 'UnEdited',
                  'Localization': 'PRESS',
                  'Metabolites': ['NAA', 'Cr', 'Glu']}
        result = backend.run_simulation(params)
        assert set(result.keys()) >= {'NAA', 'Cr', 'Glu'}
        for m in ('NAA', 'Cr', 'Glu'):
            assert np.max(np.abs(result[m])) > 0


# ============================================================ integration
@_OCTAVE
@pytest.mark.skipif(
    not os.path.exists('./example_data/BigGABA_P1P_S01/S01_PRESS_35_act.SPAR'),
    reason="Example data not available",
)
class TestMRSCloudIntegration:
    """End-to-end: REMY-extracted params → MRSCloud simulation → exporter."""

    def test_full_press(self, tmp_path):
        br = BasisREMY()
        params = br.runREMY(
            import_fpath='./example_data/BigGABA_P1P_S01/S01_PRESS_35_act.SPAR'
        )
        br.set_backend('MRSCloud')
        m, _o = br.backend.parseREMY(params)

        # Force the only vendor preset whose pulses ship with the MRSCloud repo
        # (per README: "Product sequence and rf waveform are not shared").
        sim_params = {
            **br.backend.mandatory_params,
            **{k: v for k, v in m.items() if v is not None},
            'System':         'Universal_Philips',
            'Localization':   'PRESS',
            'Sequence':       'UnEdited',
            'Spatial Points': 41,
            'Metabolites':    ['NAA', 'Cr'],
        }

        missing = br.backend.missing_pulse_files(
            sim_params['System'], sim_params['Sequence'], sim_params['Localization'])
        if missing:
            pytest.skip(f"MRSCloud pulse files not available locally: {missing}")

        result = br.backend.run_simulation(sim_params)
        assert {'NAA', 'Cr'}.issubset(result)
        for k in ('NAA', 'Cr'):
            assert np.max(np.abs(result[k])) > 0

        # Round-trip through the unified exporter
        from core.exporters import export
        out = export(result, str(tmp_path / 'mrscloud.basis'),
                     'lcmodel_basis', sim_params)
        assert os.path.exists(out)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])



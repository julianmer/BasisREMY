####################################################################################################
#                                    test_param_units.py                                           #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 27/08/26                                                                                #
#                                                                                                  #
# Purpose: Test the unified parameter units ('Center Freq' in MHz, 'Sim Centre (ppm)' in ppm),     #
#          the FSL-MRS editing-pulse offsets, STEAM mixing time handling, and protocol parsing.    #
#                                                                                                  #
####################################################################################################


#*************#
#   imports   #
#*************#
import sys, os
import pytest
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from basisremy.core.basisremy import BasisREMY
from basisremy.backends.fida_backends import FidaIdeal, FidaPressShaped


class TestCenterFreqUnits:
    """'Center Freq' must be MHz everywhere; REMY headers get normalized."""

    @pytest.fixture
    def br(self):
        return BasisREMY()

    def test_philips_hz_to_mhz(self, br):
        info = br.extract_more({'synthesizer_frequency': 127766000}, 'Philips', 'spar')
        assert info['Center Freq'] == pytest.approx(127.766)

    def test_siemens_hz_to_mhz(self, br):
        info = br.extract_more({'lFrequency': 123255780}, 'Siemens', 'rda')
        assert info['Center Freq'] == pytest.approx(123.25578)

    def test_siemens_already_mhz_kept(self, br):
        info = br.extract_more({'MRFrequency': 123.255780}, 'Siemens', 'rda')
        assert info['Center Freq'] == pytest.approx(123.25578)

    def test_ge_tenths_of_hz_to_mhz(self, br):
        # Pfiles store the frequency in 0.1 Hz units and use the GE key
        info = br.extract_more({'rhr_rh_ps_mps_freq': 1277360000}, 'GE', '7')
        assert info['Center Freq'] == pytest.approx(127.736)

    def test_nifti_list_and_missing_flip_angle(self, br):
        # SpectrometerFrequency may be a list (MHz per spec);
        # ExcitationFlipAngle is optional and must not crash the import
        info = br.extract_more({'SpectrometerFrequency': [127.732, 127.732]},
                               'NIfTI', 'json')
        assert info['Center Freq'] == pytest.approx(127.732)
        assert 'ExcitationFlipAngle' not in info


class TestFslmrsSequenceJson:
    """Editing pulses sit relative to the 4.65 ppm carrier; STEAM honors TM."""

    @pytest.fixture
    def backend(self):
        br = BasisREMY()
        br.set_backend('FSL-MRS')
        return br.backend

    def _params(self, **over):
        base = {'Sequence': 'PRESS', 'Samples': 2048, 'Bandwidth': 2000,
                'Bfield': 3.0, 'TE': 68}
        base.update(over)
        return base

    def test_mega_press_edit_offset_relative_to_carrier(self, backend):
        seq = backend._generate_sequence_json(
            self._params(Sequence='MEGA-PRESS', **{'Edit Frequency': 1.9}))
        offsets = [rf['frequencyOffset'] for rf in seq['RF']
                   if rf['frequencyOffset'] != 0]
        expected = (1.9 - seq['centralShift']) * 3.0 * 42.577  # ≈ -351.3 Hz
        assert offsets, "MEGA-PRESS must contain editing pulses"
        for off in offsets:
            assert off == pytest.approx(expected)
            assert off < 0  # 1.9 ppm lies below the 4.65 ppm carrier

    def test_hermes_edit_offsets(self, backend):
        seq = backend._generate_sequence_json(self._params(Sequence='HERMES'))
        offsets = sorted(rf['frequencyOffset'] for rf in seq['RF']
                         if rf['frequencyOffset'] != 0)
        cs = seq['centralShift']
        assert offsets[0] == pytest.approx((1.9 - cs) * 3.0 * 42.577)   # GABA
        assert offsets[-1] == pytest.approx((4.56 - cs) * 3.0 * 42.577)  # GSH

    def test_steam_uses_tm(self, backend):
        seq = backend._generate_sequence_json(
            self._params(Sequence='STEAM', TE=20, TM=50))
        # TM sits in the middle delay; the 10 us ideal pulse length is taken off
        # the pulse-end-to-pulse-start delay, so compare with that tolerance.
        assert seq['delays'][1] == pytest.approx(0.05, abs=2e-5)


class TestFidaIdeal:
    """LASER parsing and the dedicated STEAM mixing-time parameter."""

    def test_laser_protocol_not_shadowed_by_se(self):
        b = FidaIdeal()
        assert b.parseProtocol('LASER') == 'LASER'
        assert b.parseProtocol('laser_te30') == 'LASER'
        assert b.parseProtocol('svs_se_30') == 'Spin Echo'
        assert b.parseProtocol('PRESS') == 'PRESS'

    def test_steam_build_args_use_tm(self):
        b = FidaIdeal()
        params = {'Sequence': 'STEAM', 'Samples': 2048, 'Bandwidth': 2000,
                  'Bfield': 3.0, 'Linewidth': 1, 'TE': 20, 'TE2': 0, 'TM': 50}
        args = b._build_args(params, 'NAA')
        assert args[5] == pytest.approx(50.0)   # tau2 == TM for STEAM

    def test_press_build_args_use_te2(self):
        b = FidaIdeal()
        params = {'Sequence': 'PRESS', 'Samples': 2048, 'Bandwidth': 2000,
                  'Bfield': 3.0, 'Linewidth': 1, 'TE': 35, 'TE2': 17, 'TM': 50}
        args = b._build_args(params, 'NAA')
        assert args[5] == pytest.approx(17.0)   # tau2 == TE2 otherwise

    def test_tm_hidden_unless_steam(self):
        b = FidaIdeal()
        b.mandatory_params['Sequence'] = 'PRESS'
        assert 'TM' not in b.get_params_for_mode()
        b.mandatory_params['Sequence'] = 'STEAM'
        assert 'TM' in b.get_params_for_mode()


class TestFidaShapedPpmCentre:
    """Shaped backends carry the ppm centre under 'Sim Centre (ppm)'."""

    def test_press_shaped_passes_ppm_centre(self, tmp_path):
        b = FidaPressShaped()
        pulse = tmp_path / "dummy.pta"
        pulse.write_text("dummy")   # staged into the workdir by _build_args
        params = dict(b.mandatory_params)
        params.update({'Samples': 2048, 'Bandwidth': 2000, 'Bfield': 3.0,
                       'TE': 30, 'Path to Pulse': str(pulse)})
        args = b._build_args(params, 'NAA')
        assert args[-1] == pytest.approx(4.65)

    def test_no_center_freq_key_in_shaped_schema(self):
        b = FidaPressShaped()
        assert 'Center Freq' not in b.mandatory_params
        assert b.mandatory_params['Sim Centre (ppm)'] == pytest.approx(4.65)


class TestParamReset:
    """A new file import starts from defaults but keeps metabolite curation."""

    def test_reset_restores_defaults_and_keeps_metabs(self):
        br = BasisREMY()
        br.set_backend('FSL-MRS')
        default_te = br.backend.mandatory_params['TE']
        br.backend.mandatory_params['TE'] = 99.9
        br.backend.mandatory_params['Metabolites'] = ['NAA']
        br.reset_backend_params()
        assert br.backend.mandatory_params['TE'] == default_te
        assert br.backend.mandatory_params['Metabolites'] == ['NAA']

    def test_reset_does_not_share_state_between_calls(self):
        br = BasisREMY()
        br.set_backend('FSL-MRS')
        br.reset_backend_params()
        br.backend.mandatory_params['TE'] = 55
        br.reset_backend_params()
        assert br.backend.mandatory_params['TE'] != 55

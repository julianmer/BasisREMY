####################################################################################################
#                                      test_fslmrs.py                                              #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 18/02/26                                                                                #
#                                                                                                  #
# Purpose: Test FSL-MRS backend for basis set simulation.                                          #
#                                                                                                  #
####################################################################################################


#*************#
#   imports   #
#*************#
import sys, os
import pytest
import numpy as np
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

from core.basisremy import BasisREMY


class TestFSLMRSBackend:
    """Test suite for FSL-MRS backend"""

    @pytest.fixture
    def basisremy(self):
        """Initialize BasisREMY with FSL-MRS backend"""
        br = BasisREMY()
        br.set_backend('FSL-MRS')
        return br

    # ------------------------------------------------------------------
    #  Initialization
    # ------------------------------------------------------------------
    def test_backend_initialization(self, basisremy):
        """Test that FSL-MRS backend initializes correctly"""
        assert basisremy.backend.name == 'FSL-MRS'
        assert basisremy.backend.requires_octave == False
        assert 'PRESS' in basisremy.backend.supported_sequences
        assert 'STEAM' in basisremy.backend.supported_sequences
        assert 'sLASER' in basisremy.backend.supported_sequences
        assert 'Custom' in basisremy.backend.supported_sequences

    def test_backend_in_available_list(self):
        """Test that FSL-MRS backend is available"""
        br = BasisREMY()
        assert 'FSL-MRS' in br.available_backends

    # ------------------------------------------------------------------
    #  Mode system
    # ------------------------------------------------------------------
    def test_modes(self, basisremy):
        """Test mode system"""
        b = basisremy.backend
        assert b.modes == ['Simple', 'Template', 'Custom']
        assert b.current_mode == 'Simple'

    def test_mode_switching(self, basisremy):
        """Test switching modes changes displayed parameters"""
        b = basisremy.backend

        b.set_mode('Simple')
        p = b.get_params_for_mode()
        assert 'Sequence' in p
        assert 'TE' in p
        assert 'Bfield' in p
        assert 'Custom Sequence' not in p

        b.set_mode('Template')
        p = b.get_params_for_mode()
        assert 'Template File' in p
        assert 'Sequence' not in p
        assert 'TE' not in p

        b.set_mode('Custom')
        p = b.get_params_for_mode()
        assert 'Custom Sequence' in p
        assert 'Sequence' not in p

    def test_invalid_mode(self, basisremy):
        """Test that invalid modes are rejected"""
        with pytest.raises(ValueError):
            basisremy.backend.set_mode('InvalidMode')

    # ------------------------------------------------------------------
    #  REMY parsing
    # ------------------------------------------------------------------
    def test_parse_remy_basic(self, basisremy):
        """Test basic REMY parameter parsing"""
        remy_params = {
            'NumberOfDatapoints': 2048,
            'SpectralWidth': 2000,
            'B0': 3.0,
            'TE': 35,
            'Nucleus': '1H',
            'Protocol': 'PRESS',
        }
        parsed, opt = basisremy.backend.parseREMY(remy_params)
        assert parsed['Samples'] == 2048
        assert parsed['Bandwidth'] == 2000
        assert parsed['Bfield'] == 3.0
        assert parsed['TE'] == 35
        assert parsed['Sequence'] == 'PRESS'

    def test_parse_remy_calculates_center_freq(self, basisremy):
        """Test center frequency calculation from B0"""
        parsed, _ = basisremy.backend.parseREMY({'B0': 3.0, 'TE': 35})
        assert abs(parsed['Center Freq'] - 42.577 * 3.0) < 0.01

    def test_parse_protocol(self, basisremy):
        """Test protocol string parsing"""
        b = basisremy.backend
        assert b.parseProtocol('PRESS_35ms') == 'PRESS'
        assert b.parseProtocol('steam') == 'STEAM'
        assert b.parseProtocol('sLASER_test') == 'sLASER'
        assert b.parseProtocol('MEGA_PRESS') == 'MEGA-PRESS'
        assert b.parseProtocol('HERMES') == 'HERMES'
        assert b.parseProtocol('unknown') is None

    # ------------------------------------------------------------------
    #  denmatsim import
    # ------------------------------------------------------------------
    def test_denmatsim_importable(self):
        """Test that denmatsim can be imported from the submodule"""
        # The import happens at module level in fslmrs_backend.py
        from denmatsim import simseq, utils as simutils
        spins = simutils.readBuiltInSpins()
        assert len(spins) > 20  # should have 30+ metabolites

    # ------------------------------------------------------------------
    #  Sequence JSON generation
    # ------------------------------------------------------------------
    @pytest.mark.parametrize("seq_name", [
        'PRESS', 'STEAM', 'sLASER', 'LASER',
        'MEGA-PRESS', 'HERMES', 'HERCULES', 'MEGA-sLASER',
    ])
    def test_generate_sequence_json(self, basisremy, seq_name):
        """Test that sequence JSONs are structurally valid for denmatsim"""
        params = {
            'Sequence': seq_name, 'TE': 35, 'Bandwidth': 2000,
            'Samples': 2048, 'Bfield': 3.0, 'TM': 10,
        }
        seq = basisremy.backend._generate_sequence_json(params)

        assert seq['B0'] == 3.0
        assert seq['Rx_SW'] == 2000
        assert seq['Rx_Points'] == 2048
        assert 'RF' in seq

        n_rf = len(seq['RF'])
        assert len(seq['delays']) == n_rf, \
            f"delays has {len(seq['delays'])} entries but there are {n_rf} RF pulses"
        assert len(seq['rephaseAreas']) == n_rf, \
            f"rephaseAreas has {len(seq['rephaseAreas'])} entries but there are {n_rf} RF pulses"
        assert len(seq['CoherenceFilter']) == n_rf, \
            f"CoherenceFilter has {len(seq['CoherenceFilter'])} entries but there are {n_rf} RF pulses"

    # ------------------------------------------------------------------
    #  Actual simulations (per sequence)
    # ------------------------------------------------------------------
    @pytest.mark.parametrize("seq_name", [
        'PRESS', 'STEAM', 'sLASER', 'LASER', 'MEGA-PRESS',
    ])
    def test_simulation_produces_signal(self, basisremy, seq_name):
        """Test that simulation produces non-zero FID for each sequence"""
        params = {
            'Sequence': seq_name, 'TE': 35, 'Bandwidth': 2000, 'Samples': 2048,
            'Bfield': 3.0, 'Nucleus': '1H', 'Center Freq': 127.7, 'TM': 10,
            'Metabolites': ['NAA'],
        }
        result = basisremy.backend.run_simulation(params)

        assert isinstance(result, dict)
        assert 'NAA' in result
        assert result['NAA'].shape == (2048,)
        assert np.max(np.abs(result['NAA'])) > 0, "FID should be non-zero"

    def test_simulation_multiple_metabolites(self, basisremy):
        """Test simulation with multiple metabolites"""
        params = {
            'Sequence': 'PRESS', 'TE': 35, 'Bandwidth': 2000, 'Samples': 2048,
            'Bfield': 3.0, 'Nucleus': '1H', 'Center Freq': 127.7,
            'Metabolites': ['NAA', 'Cr', 'Glu'],
        }
        result = basisremy.backend.run_simulation(params)

        assert len(result) == 3
        for metab in ['NAA', 'Cr', 'Glu']:
            assert metab in result
            assert result[metab].shape == (2048,)
            assert np.max(np.abs(result[metab])) > 0

    def test_simulation_writes_intermediate_artefacts(self, basisremy):
        """Backend should allocate an internal workdir; user export is separate."""
        params = {
            'Sequence': 'PRESS', 'TE': 35, 'Bandwidth': 2000, 'Samples': 2048,
            'Bfield': 3.0, 'Nucleus': '1H', 'Center Freq': 127.7,
            'Metabolites': ['NAA'],
        }
        basisremy.backend.run_simulation(params)
        wd = basisremy.backend._workdir
        assert wd is not None and os.path.isdir(wd)
        # The PRESS sequence JSON is written into the internal workdir
        assert os.path.exists(os.path.join(wd, 'PRESS_sequence.json'))

    # ------------------------------------------------------------------
    #  Parameter structure
    # ------------------------------------------------------------------
    def test_mandatory_params(self, basisremy):
        """Test mandatory parameters are defined (export-related ones are NOT)"""
        required = ['Sequence', 'Samples', 'Bandwidth', 'Bfield', 'TE',
                     'Nucleus', 'Center Freq', 'Metabolites']
        for p in required:
            assert p in basisremy.backend.mandatory_params
        for forbidden in ('Output Path', 'Output Format', 'Make .raw', 'Add Ref.'):
            assert forbidden not in basisremy.backend.mandatory_params, \
                f"{forbidden!r} leaked back into mandatory_params"

    def test_metabolite_library(self, basisremy):
        """Test metabolite library is extensive"""
        common = ['NAA', 'Cr', 'Glu', 'Gln', 'Ins', 'GABA', 'GSH', 'Lac']
        for m in common:
            assert m in basisremy.backend.default_metabolites
        assert len(basisremy.backend.default_metabolites) > 30


@pytest.mark.skipif(
    not os.path.exists('./example_data/BigGABA_P1P_S01/S01_PRESS_35_act.SPAR'),
    reason="Example data not available"
)
class TestFSLMRSIntegration:
    """Integration tests with real data"""

    def test_full_workflow(self):
        """Test complete workflow: data → REMY → FSL-MRS simulation"""
        br = BasisREMY()
        params = br.runREMY(
            import_fpath='./example_data/BigGABA_P1P_S01/S01_PRESS_35_act.SPAR'
        )
        br.set_backend('FSL-MRS')
        parsed, opt = br.backend.parseREMY(params)

        sim_params = {
            'Sequence': parsed.get('Sequence') or 'PRESS',
            'Samples': parsed.get('Samples', 2048),
            'Bandwidth': parsed.get('Bandwidth', 2000),
            'Bfield': parsed.get('Bfield', 3.0),
            'TE': parsed.get('TE', 35),
            'Nucleus': parsed.get('Nucleus', '1H'),
            'Center Freq': parsed.get('Center Freq', 127.7),
            'Metabolites': ['NAA', 'Cr'],
        }
        result = br.backend.run_simulation(sim_params)

        assert len(result) == 2
        for metab in ['NAA', 'Cr']:
            assert metab in result
            assert np.max(np.abs(result[metab])) > 0


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

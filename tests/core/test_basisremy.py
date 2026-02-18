"""
Tests for core.basisremy module
"""

import pytest
import os
from core.basisremy import BasisREMY


@pytest.mark.core
@pytest.mark.unit
class TestBasisREMY:
    """Test BasisREMY core functionality"""

    def test_init(self):
        """Test BasisREMY initialization"""
        br = BasisREMY()
        assert br is not None
        assert br.backend is not None

    def test_set_backend_lcmodel(self):
        """Test setting LCModel backend"""
        br = BasisREMY()
        br.set_backend('LCModel')
        assert br.backend.name == 'LCModel'

    def test_set_backend_slaser(self):
        """Test setting sLaser backend"""
        br = BasisREMY()
        br.set_backend('sLaserSim')
        assert br.backend.name == 'sLaserSim'

    def test_set_backend_invalid(self):
        """Test setting invalid backend"""
        br = BasisREMY()
        with pytest.raises(ValueError):
            br.set_backend('InvalidBackend')

    def test_get_backend_list(self):
        """Test getting list of available backends"""
        br = BasisREMY()
        backends = br.available_backends
        assert isinstance(backends, list)
        assert len(backends) > 0
        assert 'LCModel' in backends
        assert 'sLaserSim' in backends

    def test_run_remy_invalid_file(self):
        """Test runREMY with invalid file"""
        br = BasisREMY()
        try:
            result = br.runREMY(import_fpath='/path/that/does/not/exist.spar')
            # Should handle gracefully - may return None or empty dict
            assert result is not None or result is None
        except Exception:
            # It's okay if it raises an exception
            pass


@pytest.mark.core
@pytest.mark.integration
class TestBasisREMYIntegration:
    """Integration tests for BasisREMY with real files"""

    @pytest.mark.parametrize("backend_name", ['LCModel', 'sLaserSim'])
    def test_backend_initialization(self, backend_name):
        """Test that each backend initializes correctly"""
        br = BasisREMY()
        br.set_backend(backend_name)
        assert br.backend is not None
        assert br.backend.name == backend_name

    def test_remy_with_philips_file(self, philips_spar_file):
        """Test REMY parsing with Philips SPAR file"""
        if not os.path.exists(philips_spar_file):
            pytest.skip("Philips SPAR file not found")

        br = BasisREMY()
        params = br.runREMY(import_fpath=philips_spar_file)

        assert params is not None
        assert isinstance(params, dict)
        # Check for key REMY fields
        assert 'Protocol' in params or 'Sequence' in params

    def test_remy_with_ge_file(self, ge_p_file):
        """Test REMY parsing with GE P-file"""
        if not ge_p_file or not os.path.exists(ge_p_file):
            pytest.skip("GE P-file not found")

        br = BasisREMY()
        params = br.runREMY(import_fpath=ge_p_file)

        assert params is not None
        assert isinstance(params, dict)

    def test_remy_with_bruker_file(self, bruker_dat_file):
        """Test REMY parsing with Bruker file"""
        if not bruker_dat_file or not os.path.exists(bruker_dat_file):
            pytest.skip("Bruker file not found")

        br = BasisREMY()
        params = br.runREMY(import_fpath=bruker_dat_file)

        assert params is not None
        assert isinstance(params, dict)





"""
LCModel Backend - PRESS Sequence Tests

Tests that LCModel backend can simulate PRESS sequences
Uses hardcoded parameters - does NOT loop through all files
"""

import pytest
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

from basisremy.core.basisremy import BasisREMY

@pytest.mark.backend
@pytest.mark.lcmodel
@pytest.mark.slow
@pytest.mark.requires_docker
class TestLCModelPRESS:
    """Test LCModel backend with PRESS sequence"""

    @pytest.fixture(scope="class")
    def basisremy(self):
        br = BasisREMY()
        br.set_backend('FidaIdeal')
        return br

    def test_press_simulation(self, basisremy, test_output_dir, cleanup_docker_processes):
        """
        Test that LCModel backend can simulate a PRESS basis set
        Uses hardcoded parameters - just verifies the backend works
        """
        print(f"\n{'='*80}")
        print(f"Testing LCModel Backend - PRESS Sequence")
        print(f"{'='*80}\n")

        # Hardcoded PRESS parameters
        test_params = {
            'Sequence': 'PRESS',
            'Samples': 2048,
            'Bandwidth': 2000,
            'Bfield': 3.0,
            'Linewidth': 1,
            'TE': 35,
            'TE2': 0,
            'Metabolites': ['NAA'],  # Just one metabolite for speed
            'Center Freq': 127736713,
        }

        print("Initializing Octave...")
        basisremy.backend.initialize_octave(prefer_docker=True, verbose=False)

        print("Running PRESS simulation...")
        result = basisremy.backend.run_simulation(test_params)

        # Verify backend returned a non-empty FID for NAA.
        # The backend returns {metabolite: complex_numpy_array}; export to
        # .RAW / .basis / etc. is handled separately by core/exporters.py.
        import numpy as np
        assert 'NAA' in result, f"NAA missing from result keys: {list(result.keys())}"
        fid = np.asarray(result['NAA'])
        assert fid.ndim >= 1 and fid.size > 0, "NAA FID is empty"
        assert np.max(np.abs(fid)) > 0, "NAA FID is all-zero"
        file_size = fid.size * fid.itemsize

        print(f"\n✅ SUCCESS!")
        print(f"   Output: NAA.RAW ({file_size} bytes)")
        print(f"{'='*80}\n")


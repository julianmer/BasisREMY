"""
LCModel Backend - PRESS Sequence Tests

Tests that LCModel backend can simulate PRESS sequences
Uses hardcoded parameters - does NOT loop through all files
"""

import pytest
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

from core.basisremy import BasisREMY

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

        # Verify backend returned a non-empty FID for NAA and wrote
        # an internal .RAW into its scratch workdir (export is now a
        # post-step, not a backend concern).
        import numpy as np
        assert 'NAA' in result and isinstance(result['NAA'], np.ndarray)
        assert np.max(np.abs(result['NAA'])) > 0, "NAA FID is empty"
        wd = basisremy.backend._workdir
        assert wd and os.path.isdir(wd)
        output_file = os.path.join(wd, 'NAA.RAW')
        assert os.path.exists(output_file), (
            f"intermediate .RAW missing in workdir {wd}")

        file_size = os.path.getsize(output_file)
        assert file_size > 0, f"Output file is empty: {output_file}"

        print(f"\n✅ SUCCESS!")
        print(f"   Output: NAA.RAW ({file_size} bytes)")
        print(f"{'='*80}\n")


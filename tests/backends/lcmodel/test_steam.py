"""
LCModel Backend - STEAM Sequence Tests

Tests that LCModel backend can simulate STEAM sequences
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
class TestLCModelSTEAM:
    """Test LCModel backend with STEAM sequence"""

    @pytest.fixture(scope="class")
    def basisremy(self):
        br = BasisREMY()
        br.set_backend('LCModel')
        return br

    def test_steam_simulation(self, basisremy, test_output_dir, cleanup_docker_processes):
        """
        Test that LCModel backend can simulate a STEAM basis set
        Uses hardcoded parameters - just verifies the backend works
        """
        print(f"\n{'='*80}")
        print(f"Testing LCModel Backend - STEAM Sequence")
        print(f"{'='*80}\n")

        # Hardcoded STEAM parameters
        test_params = {
            'Sequence': 'STEAM',
            'Samples': 2048,
            'Bandwidth': 2000,
            'Bfield': 3.0,
            'Linewidth': 1,
            'TE': 20,
            'TE2': 0,
            'Add Ref.': 'No',
            'Output Path': os.path.join(test_output_dir, 'lcmodel_steam_test'),
            'Metabolites': ['NAA'],
            'Center Freq': 127736713,
            'Make .raw': 'Yes',
        }

        os.makedirs(test_params['Output Path'], exist_ok=True)

        print("Initializing Octave...")
        basisremy.backend.initialize_octave(prefer_docker=True, verbose=False)

        print("Running STEAM simulation...")
        result = basisremy.backend.run_simulation(test_params)

        # Verify output file was created
        output_file = os.path.join(test_params['Output Path'], 'NAA.RAW')

        assert os.path.exists(output_file), (
            f"❌ Output file not created: {output_file}\n"
            f"Simulation ran but produced no output."
        )

        file_size = os.path.getsize(output_file)
        assert file_size > 0, f"Output file is empty: {output_file}"

        print(f"\n✅ SUCCESS!")
        print(f"   Output: NAA.RAW ({file_size} bytes)")
        print(f"   Location: {test_params['Output Path']}")
        print(f"{'='*80}\n")



"""
sLaser Backend - sLASER Sequence Tests

Tests that sLaser backend can simulate sLASER sequences
Uses hardcoded parameters - does NOT loop through all files
"""

import pytest
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

from core.basisremy import BasisREMY


@pytest.mark.backend
@pytest.mark.slaser
@pytest.mark.slow
@pytest.mark.requires_docker
class TestSLaserSLASER:
    """Test sLaser backend with sLASER sequence"""

    @pytest.fixture(scope="class")
    def basisremy(self):
        br = BasisREMY()
        br.set_backend('sLaserSim')
        return br

    @pytest.fixture(scope="class")
    def pulse_file(self, project_root_dir):
        return os.path.join(project_root_dir, 'externals', 'jbss', 'my_pulse', 'standardized_goia.txt')

    def test_slaser_simulation(self, basisremy, pulse_file, test_output_dir, cleanup_docker_processes):
        """
        Test that sLaser backend can simulate an sLASER basis set
        Uses hardcoded parameters - just verifies the backend works
        """
        if not os.path.exists(pulse_file):
            pytest.skip(f"Pulse file not found: {pulse_file}")

        print(f"\n{'='*80}")
        print(f"Testing sLaser Backend - sLASER Sequence")
        print(f"{'='*80}\n")

        # Hardcoded sLASER parameters
        test_params = {
            'System': 'Philips',
            'Sequence': 'sLASER',
            'Basis Name': 'test_slaser.basis',
            'B1max': 22.0,
            'Flip Angle': 180.0,
            'RefTp': 4.5008,
            'Samples': 2048,
            'Bandwidth': 2000,
            'Linewidth': 1,
            'Bfield': 3.0,
            'thkX': 2.0,
            'thkY': 2.0,
            'fovX': 3.0,
            'fovY': 3.0,
            'nX': 16.0,  # Small grid for speed
            'nY': 16.0,
            'TE': 35,
            'Center Freq': 127736713,
            'Metabolites': ['NAA'],
            'Tau 1': 15.0,
            'Tau 2': 13.0,
            'Path to Pulse': pulse_file,
            'Output Path': os.path.join(test_output_dir, 'slaser_test'),
            'Make .raw': 'No',
        }

        os.makedirs(test_params['Output Path'], exist_ok=True)

        print("Initializing Octave...")
        basisremy.backend.initialize_octave(prefer_docker=True, verbose=False)

        print("Running sLASER simulation...")
        result = basisremy.backend.run_simulation(test_params)

        # Verify output file was created
        output_file = os.path.join(test_params['Output Path'], 'test_slaser.basis')

        assert os.path.exists(output_file), (
            f"❌ Output file not created: {output_file}\n"
            f"Simulation ran but produced no output."
        )

        file_size = os.path.getsize(output_file)
        assert file_size > 0, f"Output file is empty: {output_file}"

        print(f"\n✅ SUCCESS!")
        print(f"   Output: test_slaser.basis ({file_size} bytes)")
        print(f"   Location: {test_params['Output Path']}")
        print(f"{'='*80}\n")

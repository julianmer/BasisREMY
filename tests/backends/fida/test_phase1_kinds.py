"""
FID-A Phase-1 kinds — live simulation tests
Covers the newly implemented LASER (ideal AFP), Spin Echo xN, and One pulse
backends with hardcoded parameters (single metabolite for speed).
"""
import pytest
import os
import sys
import numpy as np
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))
from basisremy.core.basisremy import BasisREMY


def _assert_valid_fid(result, metab='NAA'):
    assert metab in result, f"{metab} missing from result keys: {list(result.keys())}"
    fid = np.asarray(result[metab])
    assert fid.ndim >= 1 and fid.size > 0, f"{metab} FID is empty"
    assert np.max(np.abs(fid)) > 0, f"{metab} FID is all-zero"


@pytest.mark.backend
@pytest.mark.slow
@pytest.mark.requires_docker
class TestFidaPhase1Kinds:

    def test_laser_ideal_afp(self, cleanup_docker_processes):
        br = BasisREMY()
        br.set_backend('FidaLaser')
        br.backend.initialize_octave(prefer_docker=True)
        result = br.backend.run_simulation({
            'Samples': 2048, 'Bandwidth': 2000, 'Bfield': 3.0,
            'Linewidth': 1.0, 'TE': 30, 'Metabolites': ['NAA'],
        })
        _assert_valid_fid(result)

    def test_spinecho_xn(self, cleanup_docker_processes):
        br = BasisREMY()
        br.set_backend('FidaSpinEchoXN')
        br.backend.initialize_octave(prefer_docker=True)
        result = br.backend.run_simulation({
            'Samples': 2048, 'Bandwidth': 2000, 'Bfield': 3.0,
            'Linewidth': 1.0, 'Tau': 15.0, 'Nechoes': 2,
            'Metabolites': ['NAA'],
        })
        _assert_valid_fid(result)

    def test_onepulse_ideal(self, cleanup_docker_processes):
        br = BasisREMY()
        br.set_backend('FidaOnePulse')
        br.backend.initialize_octave(prefer_docker=True)
        result = br.backend.run_simulation({
            'Samples': 2048, 'Bandwidth': 2000, 'Bfield': 3.0,
            'Linewidth': 1.0, 'Metabolites': ['NAA'],
        })
        _assert_valid_fid(result)

    def test_onepulse_non_ideal_mode_raises(self):
        br = BasisREMY()
        br.set_backend('FidaOnePulse')
        br.backend.set_mode('Shaped')
        with pytest.raises(NotImplementedError, match='Shaped'):
            br.backend._build_args({'Samples': 2048, 'Bandwidth': 2000,
                                    'Bfield': 3.0}, 'NAA')

    def test_megapress_ideal_subspectra(self, cleanup_docker_processes):
        br = BasisREMY()
        br.set_backend('FidaMegaPressIdeal')
        br.backend.initialize_octave(prefer_docker=True)
        result = br.backend.run_simulation({
            'Samples': 2048, 'Bandwidth': 2000, 'Bfield': 3.0,
            'Linewidth': 1.0, 'TE': 68,
            'Edit On': 1.9, 'Edit Bandwidth (ppm)': 1.0,
            'Metabolites': ['GABA'],
        })
        for sub in ('ON', 'OFF', 'DIFF'):
            _assert_valid_fid(result, f'GABA ({sub})')
        on = np.asarray(result['GABA (ON)'])
        off = np.asarray(result['GABA (OFF)'])
        diff = np.asarray(result['GABA (DIFF)'])
        assert np.allclose(diff, on - off)
        # editing must actually change the GABA signal
        assert np.max(np.abs(diff)) > 1e-6 * np.max(np.abs(off))

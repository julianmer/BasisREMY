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

    # ---- Phase 2: shaped kinds (small grids for speed) ----
    # FID-A grids are linspace(-fov/2, fov/2, n): with n = 2 the two points
    # sit at the FOV edge, so the FOV (1.5 cm) is kept inside the 2 cm slab —
    # with fov 3 both points were outside it and the signals were ~0.

    def test_spinecho_shaped(self, cleanup_docker_processes, project_root_dir):
        pulse = os.path.join(project_root_dir, 'externals', 'fidA',
                             'rfPulseTools', 'rfPulses', 'sampleRefocPulse.pta')
        if not os.path.exists(pulse):
            pytest.skip("sampleRefocPulse.pta not available")
        br = BasisREMY()
        br.set_backend('FidaSpinEchoShaped')
        br.backend.initialize_octave(prefer_docker=True)
        result = br.backend.run_simulation({
            'Samples': 1024, 'Bandwidth': 2000, 'Bfield': 3.0,
            'Linewidth': 1.0, 'TE': 30, 'RefTp': 5.0,
            'thkX': 2.0, 'fovX': 1.5, 'nX': 2,
            'Path to Pulse': pulse, 'Metabolites': ['NAA'],
        })
        _assert_valid_fid(result)

    def test_steam_shaped(self, cleanup_docker_processes, project_root_dir):
        pulse = os.path.join(project_root_dir, 'externals', 'fidA',
                             'rfPulseTools', 'rfPulses', 'sampleExcPulse.pta')
        if not os.path.exists(pulse):
            pytest.skip("sampleExcPulse.pta not available")
        br = BasisREMY()
        br.set_backend('FidaSteamShaped')
        br.backend.initialize_octave(prefer_docker=True)
        result = br.backend.run_simulation({
            'Samples': 1024, 'Bandwidth': 2000, 'Bfield': 3.0,
            'Linewidth': 1.0, 'TE': 20, 'TM': 10, 'RefTp': 5.0,
            'thkX': 2.0, 'thkY': 2.0, 'fovX': 1.5, 'fovY': 1.5,
            'nX': 2, 'nY': 2, 'Flip Angle': 90.0,
            'Sim Centre (ppm)': 4.65,
            'Path to Pulse': pulse, 'Metabolites': ['NAA'],
        })
        _assert_valid_fid(result)

    def test_semilaser_shaped_goia(self, cleanup_docker_processes, project_root_dir):
        pulse = os.path.join(project_root_dir, 'externals', 'jbss',
                             'my_pulse', 'standardized_goia.txt')
        if not os.path.exists(pulse):
            pytest.skip("GOIA pulse not available")
        br = BasisREMY()
        br.set_backend('FidaSemiLaserShaped')
        br.backend.initialize_octave(prefer_docker=True)
        result = br.backend.run_simulation({
            'Samples': 1024, 'Bandwidth': 2000, 'Bfield': 3.0,
            'Linewidth': 1.0, 'TE': 35, 'RefTp': 4.5008,
            'thkX': 2.0, 'thkY': 2.0, 'fovX': 1.5, 'fovY': 1.5,
            'nX': 2, 'nY': 2, 'Flip Angle': 180.0,
            'Sim Centre (ppm)': 4.65,
            'Path to Pulse': pulse, 'Metabolites': ['NAA'],
        })
        _assert_valid_fid(result)

    def test_megapress_shaped_edit(self, cleanup_docker_processes, project_root_dir):
        pulse = os.path.join(project_root_dir, 'externals', 'fidA',
                             'rfPulseTools', 'rfPulses', 'sampleEditPulse.pta')
        if not os.path.exists(pulse):
            pytest.skip("sampleEditPulse.pta not available")
        br = BasisREMY()
        br.set_backend('FidaMegaPressShaped')
        br.backend.initialize_octave(prefer_docker=True)
        result = br.backend.run_simulation({
            'Samples': 1024, 'Bandwidth': 2000, 'Bfield': 3.0,
            'Linewidth': 1.0, 'TE': 68,
            'Edit Pulse Path': pulse, 'Edit Tp': 20.0,
            'Edit On': 1.9, 'Edit Off': 7.5,
            'Sim Centre (ppm)': 4.65,
            'Metabolites': ['GABA'],
        })
        for sub in ('ON', 'OFF', 'DIFF'):
            _assert_valid_fid(result, f'GABA ({sub})')
        diff = np.asarray(result['GABA (DIFF)'])
        off = np.asarray(result['GABA (OFF)'])
        # the shaped editing pulse must actually edit GABA
        assert np.max(np.abs(diff)) > 1e-6 * np.max(np.abs(off))

    def test_megaspecial_shaped(self, cleanup_docker_processes, project_root_dir):
        rfdir = os.path.join(project_root_dir, 'externals', 'fidA',
                             'rfPulseTools', 'rfPulses')
        refoc = os.path.join(rfdir, 'sampleRefocPulse.pta')
        edit = os.path.join(rfdir, 'sampleEditPulse.pta')
        if not (os.path.exists(refoc) and os.path.exists(edit)):
            pytest.skip("sample pulses not available")
        br = BasisREMY()
        br.set_backend('FidaMegaSpecialShaped')
        br.backend.initialize_octave(prefer_docker=True)
        result = br.backend.run_simulation({
            'Samples': 1024, 'Bandwidth': 2000, 'Bfield': 3.0,
            'Linewidth': 1.0, 'TE': 68,
            'Path to Pulse': refoc, 'RefTp': 5.0,
            'Edit Pulse Path': edit, 'Edit Tp': 14.0,
            'Edit On': 1.9, 'Edit Off': 7.5,
            'thkX': 2.0, 'fovX': 1.5, 'nX': 2,
            'Sim Centre (ppm)': 4.65,
            'Metabolites': ['GABA'],
        })
        for sub in ('ON', 'OFF', 'DIFF'):
            _assert_valid_fid(result, f'GABA ({sub})')


# ------------------------------------------------------------------ modes
# The previously gated modes (roadmap step 2): every mode maps to a real
# fida_run.m kind and shows only its own parameters.

class TestFidaModes:
    def test_onepulse_modes_map_to_kinds(self):
        from basisremy.backends.fida_backends import FidaOnePulse
        b = FidaOnePulse()
        for mode, kind in (('Ideal', 'onepulse'), ('Shaped', 'onepulse_shaped'),
                           ('Delay', 'onepulse_delay'),
                           ('Arbitrary phase', 'onepulse_arbph')):
            b.current_mode = mode
            assert b.active_kind() == kind
        assert 'Delay' not in b.get_params_for_mode('Ideal')
        assert set(b.get_params_for_mode('Delay')) & {'Delay', 'Pulse Phase',
                                                       'Path to Pulse'} == {'Delay'}
        assert {'Flip Angle', 'Path to Pulse', 'RefTp'} <= set(b.get_params_for_mode('Shaped'))
        assert 'Pulse Phase' in b.get_params_for_mode('Arbitrary phase')

    def test_onepulse_shaped_requires_pulse(self):
        from basisremy.backends.fida_backends import FidaOnePulse
        b = FidaOnePulse()
        b.current_mode = 'Shaped'
        with pytest.raises(ValueError, match='Path to Pulse'):
            b._build_args({'Samples': 2048, 'Bandwidth': 2000, 'Bfield': 3.0},
                          'NAA')

    def test_ideal_press_echo_split(self):
        from basisremy.backends.fida_backends import FidaIdeal
        b = FidaIdeal()
        base = {'Sequence': 'PRESS', 'Samples': 2048, 'Bandwidth': 2000,
                'Bfield': 3.0, 'TE': 35, 'TE2': 0}
        assert b._build_args(base, 'NAA')[4:6] == [17.5, 17.5]   # symmetric by default
        assert b._build_args({**base, 'TE2': 24}, 'NAA')[4:6] == [11.0, 24.0]
        assert b._build_args({**base, 'Sequence': 'STEAM', 'TM': 10}, 'NAA')[4:6] == [35.0, 10.0]
        assert b._build_args({**base, 'Sequence': 'Spin Echo'}, 'NAA')[4:6] == [35.0, 0.0]

    def test_semilaser_modes_map_to_kinds(self):
        from basisremy.backends.fida_backends import FidaSemiLaserShaped
        b = FidaSemiLaserShaped()
        assert b.active_kind() == 'semilaser_shaped'
        b.current_mode = 'Phase cycled'
        assert b.active_kind() == 'semilaser_shaped_phcyc'

    def test_megapress_shaped_modes(self):
        from basisremy.backends.fida_backends import FidaMegaPressShaped
        b = FidaMegaPressShaped()
        assert b.active_kind() == 'megapress_shapededit'
        p = b.get_params_for_mode('Edit-only shaped (ideal refoc)')
        assert 'Path to Pulse' not in p and 'nX' not in p
        b.current_mode = 'Full shaped (refoc + edit)'
        assert b.active_kind() == 'megapress_shaped'
        p = b.get_params_for_mode()
        assert {'Edit Pulse Path', 'Path to Pulse', 'RefTp', 'nX'} <= set(p)
        b.current_mode = 'Refoc-only shaped (ideal edit)'
        assert b.active_kind() == 'megapress_shapedrefoc'
        p = b.get_params_for_mode()
        assert 'Edit Pulse Path' not in p
        assert {'Edit On', 'Edit Bandwidth (ppm)', 'Path to Pulse'} <= set(p)
        assert set(b.file_selection) == {'Edit Pulse Path', 'Path to Pulse'}


@pytest.mark.backend
@pytest.mark.slow
class TestFidaModesLive:
    """Small-grid live runs of the newly opened modes (Docker Octave)."""

    def _rf(self, project_root_dir, name):
        p = os.path.join(project_root_dir, 'externals', 'fidA',
                         'rfPulseTools', 'rfPulses', name)
        if not os.path.exists(p):
            pytest.skip(f"{name} not available")
        return p

    def _base(self):
        return {'Samples': 1024, 'Bandwidth': 2000, 'Bfield': 3.0,
                'Linewidth': 1.0, 'Metabolites': ['Cr']}

    def test_onepulse_delay_and_phase(self, cleanup_docker_processes):
        import numpy as np
        br = BasisREMY()
        br.set_backend('FidaOnePulse')
        br.backend.initialize_octave(prefer_docker=True)
        br.backend.current_mode = 'Ideal'
        ideal = br.backend.run_simulation(self._base())
        br.backend.current_mode = 'Delay'
        delayed = br.backend.run_simulation({**self._base(), 'Delay': 0.5})
        br.backend.current_mode = 'Arbitrary phase'
        phased = br.backend.run_simulation({**self._base(), 'Pulse Phase': 90.0})
        for r in (ideal, delayed, phased):
            _assert_valid_fid(r, metab='Cr')
        # same magnitude spectrum, different phase
        s = lambda r: np.abs(np.fft.fft(r['Cr']))
        assert np.corrcoef(s(ideal), s(phased))[0, 1] > 0.999
        assert not np.allclose(ideal['Cr'], phased['Cr'])
        assert not np.allclose(ideal['Cr'], delayed['Cr'])

    def test_onepulse_shaped(self, cleanup_docker_processes, project_root_dir):
        br = BasisREMY()
        br.set_backend('FidaOnePulse')
        br.backend.initialize_octave(prefer_docker=True)
        br.backend.current_mode = 'Shaped'
        result = br.backend.run_simulation({
            **self._base(), 'Flip Angle': 90.0, 'RefTp': 5.0,
            'Path to Pulse': self._rf(project_root_dir, 'sampleExcPulse.pta'),
        })
        _assert_valid_fid(result, metab='Cr')

    def test_semilaser_phase_cycled(self, cleanup_docker_processes, project_root_dir):
        br = BasisREMY()
        br.set_backend('FidaSemiLaserShaped')
        br.backend.initialize_octave(prefer_docker=True)
        br.backend.current_mode = 'Phase cycled'
        result = br.backend.run_simulation({
            **self._base(), 'TE': 30, 'RefTp': 4.5,
            'thkX': 2.0, 'thkY': 2.0, 'fovX': 1.5, 'fovY': 1.5,
            'nX': 2, 'nY': 2, 'Flip Angle': 180.0, 'Sim Centre (ppm)': 4.65,
            'Path to Pulse': self._rf(project_root_dir, 'GOIA_tthk0.01_R120.txt'),
        })
        _assert_valid_fid(result, metab='Cr')

    def test_megapress_refoc_only_shaped(self, cleanup_docker_processes, project_root_dir):
        br = BasisREMY()
        br.set_backend('FidaMegaPressShaped')
        br.backend.initialize_octave(prefer_docker=True)
        br.backend.current_mode = 'Refoc-only shaped (ideal edit)'
        result = br.backend.run_simulation({
            **self._base(), 'TE': 68, 'Edit On': 1.9, 'Edit Bandwidth (ppm)': 1.0,
            'Path to Pulse': self._rf(project_root_dir, 'sampleRefocPulse.pta'),
            'RefTp': 5.0, 'thkX': 2.0, 'thkY': 2.0, 'fovX': 1.5, 'fovY': 1.5,
            'nX': 2, 'nY': 2,
        })
        assert {'Cr (ON)', 'Cr (OFF)', 'Cr (DIFF)'} <= set(result)
        _assert_valid_fid({'Cr': result['Cr (OFF)']}, metab='Cr')

    @pytest.mark.heavy   # 256 simulations for one metabolite — not for CI
    def test_megapress_full_shaped(self, cleanup_docker_processes, project_root_dir):
        br = BasisREMY()
        br.set_backend('FidaMegaPressShaped')
        br.backend.initialize_octave(prefer_docker=True)
        br.backend.current_mode = 'Full shaped (refoc + edit)'
        result = br.backend.run_simulation({
            **self._base(), 'TE': 68,
            'Edit Pulse Path': self._rf(project_root_dir, 'sampleEditPulse.pta'),
            'Edit Tp': 20.0, 'Edit On': 1.9, 'Edit Off': 7.5,
            'Path to Pulse': self._rf(project_root_dir, 'sampleRefocPulse.pta'),
            'RefTp': 5.0, 'thkX': 2.0, 'thkY': 2.0, 'fovX': 1.5, 'fovY': 1.5,
            'nX': 2, 'nY': 2,
        })
        assert {'Cr (ON)', 'Cr (OFF)', 'Cr (DIFF)'} <= set(result)
        _assert_valid_fid({'Cr': result['Cr (OFF)']}, metab='Cr')

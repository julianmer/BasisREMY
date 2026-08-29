####################################################################################################
#                                          test_spant.py                                           #
####################################################################################################
#                                                                                                  #
# Purpose: spant backend — parameter schema, REMY parsing, the job handed to spant_worker.R,       #
#          runtime discovery, and live simulations when a local R with spant is available.         #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

from basisremy.backends.spant_backend import SpantBackend
from basisremy.core import spant_manager


def _rscript_with_spant():
    r = spant_manager.find_rscript()
    return r if r and spant_manager.spant_version(r) else None


_BASE = {'Samples': 2048, 'Bandwidth': 2000, 'Bfield': 3.0, 'TE': 35,
         'Linewidth': 1.0, 'Metabolites': ['naa']}


# ============================================================ schema
class TestSchema:

    @pytest.fixture
    def backend(self):
        return SpantBackend()

    def test_registered(self):
        from basisremy.core.basisremy import BasisREMY
        br = BasisREMY()
        assert 'Spant' in br.backends
        assert br.backends['Spant'].category == 'Spant'
        assert 'Spant' in br.categories and br.categories['Spant'] == ['Spant']

    def test_physics_blank_by_default(self, backend):
        for k in ('Samples', 'Bandwidth', 'Bfield', 'TE', 'Sequence'):
            assert backend.mandatory_params[k] is None

    @pytest.mark.parametrize("seq,shown,hidden", [
        ('PRESS', {'Tau 1', 'Tau 2'}, {'TM', 'Edit On', 'Path to Pulse'}),
        ('PRESS shaped', {'Tau 1', 'Path to Pulse', 'RefTp', 'Flip Angle'}, {'TM', 'Edit On'}),
        ('STEAM', {'TM', 'STEAM Variant'}, {'Tau 1', 'Edit On', 'Path to Pulse'}),
        ('sLASER', set(), {'Tau 1', 'TM', 'Edit On', 'Path to Pulse'}),
        ('MEGA-PRESS', {'Edit On', 'Edit Off', 'Edit Bandwidth (Hz)'}, {'Tau 1', 'TM'}),
        ('Pulse-acquire', set(), {'TE', 'Tau 1', 'TM', 'Edit On'}),
    ])
    def test_params_per_sequence(self, backend, seq, shown, hidden):
        backend.mandatory_params['Sequence'] = seq
        params = backend.get_params_for_mode()
        assert shown <= set(params)
        assert not (hidden & set(params))

    @pytest.mark.parametrize("protocol,expected", [
        ('PRESS_35', 'PRESS'), ('svs_se_30', 'Spin Echo'), ('STEAM_TE20', 'STEAM'),
        ('semi_LASER_TE30', 'sLASER'), ('MEGA_PRESS_GABA', 'MEGA-PRESS'),
        ('UnEdited', 'PRESS'), ('HERMES', None), ('LASER_TE30', None),
    ])
    def test_map_sequence_in(self, backend, protocol, expected):
        assert backend.map_sequence_in(protocol) == expected

    def test_parse_remy(self, backend):
        m, o = backend.parseREMY({'NumberOfDatapoints': 4096, 'SpectralWidth': 4000,
                                  'B0': 2.89, 'TE': 35, 'Protocol': 'PRESS_35'})
        assert m['Samples'] == 4096 and m['Bandwidth'] == 4000 and m['TE'] == 35
        assert m['Bfield'] == 2.89 and m['Sequence'] == 'PRESS'
        assert m['Center Freq'] == pytest.approx(42.577 * 2.89)
        m, _ = backend.parseREMY({'TE': 35})
        assert m['Bfield'] is None and m['Center Freq'] is None

    def test_spant_names_get_tooltips(self):
        from basisremy.core.parameter_registry import metabolite_full_name
        assert metabolite_full_name('naa') == 'N-Acetylaspartate'
        assert metabolite_full_name('sins') == 'scyllo-Inositol'


# ============================================================ job building
class TestJob:

    @pytest.fixture
    def backend(self):
        return SpantBackend()

    def test_press_defaults_to_symmetric_echoes(self, backend):
        job = backend._build_job({**_BASE, 'Sequence': 'PRESS'}, 'naa')
        assert job['sequence'] == 'press'
        assert job['te1_s'] == pytest.approx(0.0175) and job['te2_s'] == pytest.approx(0.0175)
        assert job['ft_hz'] == pytest.approx(42.577e6 * 3.0)
        assert job['fs_hz'] == 2000 and job['n'] == 2048 and job['metabolites'] == ['naa']

    def test_press_explicit_taus_and_center_freq(self, backend):
        job = backend._build_job({**_BASE, 'Sequence': 'PRESS', 'Tau 1': 10, 'Tau 2': 25,
                                  'Center Freq': 123.2}, 'cr')
        assert job['te1_s'] == pytest.approx(0.010) and job['te2_s'] == pytest.approx(0.025)
        assert job['ft_hz'] == pytest.approx(123.2e6)

    def test_slaser_scales_spant_split(self, backend):
        job = backend._build_job({**_BASE, 'Sequence': 'sLASER', 'TE': 56}, 'naa')
        assert (job['te1_s'], job['te2_s'], job['te3_s']) == \
            pytest.approx((0.016, 0.022, 0.018))

    def test_mega_press(self, backend):
        job = backend._build_job({**_BASE, 'Sequence': 'MEGA-PRESS', 'TE': 68,
                                  'Edit On': 1.9, 'Edit Off': 7.5,
                                  'Edit Bandwidth (Hz)': 110}, 'gaba')
        assert job['sequence'] == 'mega_press'
        assert job['te1_s'] == pytest.approx(0.015) and job['te2_s'] == pytest.approx(0.053)
        assert job['edit_on_ppm'] == 1.9 and job['edit_off_ppm'] == 7.5
        assert job['edit_bw_hz'] == 110

    @pytest.mark.parametrize("variant,key", [('Standard', 'ideal'),
                                             ('Coherence filter', 'cof'),
                                             ('z-rotation (Young)', 'young')])
    def test_steam_variants(self, backend, variant, key):
        job = backend._build_job({**_BASE, 'Sequence': 'STEAM', 'TE': 20, 'TM': 10,
                                  'STEAM Variant': variant}, 'naa')
        assert job['te_s'] == pytest.approx(0.020) and job['tm_s'] == pytest.approx(0.010)
        assert job['steam_variant'] == key

    def test_press_shaped_needs_pulse(self, backend, tmp_path):
        with pytest.raises(ValueError, match='Path to Pulse'):
            backend._build_job({**_BASE, 'Sequence': 'PRESS shaped'}, 'naa')
        pulse = tmp_path / 'refoc.pta'
        pulse.write_text('##\n')
        job = backend._build_job({**_BASE, 'Sequence': 'PRESS shaped',
                                  'Path to Pulse': str(pulse), 'RefTp': 5.0,
                                  'Flip Angle': 180}, 'naa')
        assert job['pulse_format'] == 'pta' and job['pulse_dur_s'] == pytest.approx(0.005)

    def test_pulse_acquire_needs_no_te(self, backend):
        job = backend._build_job({**_BASE, 'Sequence': 'Pulse-acquire', 'TE': None}, 'naa')
        assert job['sequence'] == 'pulse_acquire' and 'te_s' not in job

    @pytest.mark.parametrize("missing", ['Samples', 'Bandwidth', 'Bfield', 'TE'])
    def test_missing_inputs_raise(self, backend, missing):
        with pytest.raises(ValueError, match=missing):
            backend._build_job({**_BASE, 'Sequence': 'PRESS', missing: None}, 'naa')

    def test_unknown_sequence_raises(self, backend):
        with pytest.raises(ValueError, match='Sequence'):
            backend._build_job({**_BASE, 'Sequence': 'HERMES'}, 'naa')


# ============================================================ runtime discovery
class TestManager:

    def test_env_override_for_rscript(self, monkeypatch, tmp_path):
        fake = tmp_path / 'Rscript'
        fake.write_text('#!/bin/sh\n')
        monkeypatch.setenv('BASISREMY_RSCRIPT', str(fake))
        assert spant_manager.find_rscript() == str(fake)

    def test_env_override_for_runtime(self, monkeypatch):
        monkeypatch.setenv('BASISREMY_SPANT_RUNTIME', 'local')
        assert spant_manager.preferred_runtime() == 'local'
        monkeypatch.setenv('BASISREMY_SPANT_RUNTIME', 'docker')
        assert spant_manager.preferred_runtime() == 'docker'

    def test_docker_preferred_then_local(self, monkeypatch):
        monkeypatch.delenv('BASISREMY_SPANT_RUNTIME', raising=False)
        monkeypatch.setattr(spant_manager, 'docker_available', lambda: True)
        monkeypatch.setattr(spant_manager, 'find_rscript', lambda: '/usr/bin/Rscript')
        assert spant_manager.preferred_runtime() == 'docker'
        monkeypatch.setattr(spant_manager, 'docker_available', lambda: False)
        assert spant_manager.preferred_runtime() == 'local'
        monkeypatch.setattr(spant_manager, 'find_rscript', lambda: None)
        with pytest.raises(spant_manager.SpantUnavailable, match='Docker'):
            spant_manager.preferred_runtime()

    def test_worker_script_ships(self):
        assert spant_manager._WORKER.is_file()


# ============================================================ live (local R + spant)
@pytest.mark.backend
@pytest.mark.slow
@pytest.mark.skipif(_rscript_with_spant() is None,
                    reason="no local R with spant (set BASISREMY_RSCRIPT to use one)")
class TestSpantLive:

    @pytest.fixture(autouse=True)
    def local_runtime(self, monkeypatch):
        monkeypatch.setenv('BASISREMY_SPANT_RUNTIME', 'local')

    @staticmethod
    def _peak_ppm(fid, bw, bfield):
        spec = np.abs(np.fft.fftshift(np.fft.fft(fid)))
        ppm = np.linspace(-bw / 2, bw / 2, fid.size) / (42.577 * bfield) + 4.65
        return ppm[np.argmax(spec)]

    def test_press_naa_cr(self):
        b = SpantBackend()
        result = b.run_simulation({**_BASE, 'Sequence': 'PRESS',
                                   'Metabolites': ['naa', 'cr']})
        assert set(result) == {'naa', 'cr'} and not b.last_failures
        for m in ('naa', 'cr'):
            fid = result[m]
            assert fid.dtype.kind == 'c' and fid.size == 2048
            assert np.max(np.abs(fid)) > 0
        # NAA singlet at 2.01 ppm, Cr methyl at 3.03 ppm on the GUI's axis
        assert abs(self._peak_ppm(result['naa'], 2000, 3.0) - 2.01) < 0.05
        assert abs(self._peak_ppm(result['cr'], 2000, 3.0) - 3.03) < 0.05

    def test_mega_press_subspectra(self):
        b = SpantBackend()
        result = b.run_simulation({**_BASE, 'Sequence': 'MEGA-PRESS', 'TE': 68,
                                   'Edit On': 1.9, 'Edit Off': 7.5,
                                   'Edit Bandwidth (Hz)': 110,
                                   'Metabolites': ['gaba']})
        assert {'gaba (ON)', 'gaba (OFF)', 'gaba (DIFF)'} <= set(result)
        assert np.allclose(result['gaba (DIFF)'], result['gaba (ON)'] - result['gaba (OFF)'])
        # editing must change the GABA signal
        assert np.max(np.abs(result['gaba (DIFF)'])) > 1e-3 * np.max(np.abs(result['gaba (OFF)']))

    def test_matches_fida_ideal_press(self):
        """Same acquisition through FID-A's ideal PRESS (Docker or local
        Octave): the singlet-dominated metabolites must agree closely."""
        from basisremy.core.basisremy import BasisREMY
        from basisremy.core.octave_manager import OctaveManager
        om = OctaveManager()
        if not (om.check_docker_availability() or om.check_local_octave_availability()):
            pytest.skip("no Octave runtime for the FID-A reference")
        spant = SpantBackend().run_simulation(
            {**_BASE, 'Sequence': 'PRESS', 'Metabolites': ['naa', 'cr', 'lac', 'glu']})
        br = BasisREMY()
        br.set_backend('FidaIdeal')
        br.backend.initialize_octave(prefer_docker=True)
        fida = br.backend.run_simulation(
            {'Sequence': 'PRESS', 'Samples': 2048, 'Bandwidth': 2000, 'Bfield': 3.0,
             'Linewidth': 1.0, 'TE': 35, 'TE2': 0, 'Metabolites': ['NAA', 'Cr', 'Lac', 'Glu']})
        ppm = np.linspace(-1000, 1000, 2048) / (42.577 * 3.0) + 4.65
        win = (ppm > 0.5) & (ppm < 4.2)
        # Glu (strongly coupled) only agrees because FidaIdeal's PRESS is a
        # symmetric TE/2 + TE/2 like spant's — it guards that echo split
        for a, b in (('naa', 'NAA'), ('cr', 'Cr'), ('lac', 'Lac'), ('glu', 'Glu')):
            s1 = np.abs(np.fft.fftshift(np.fft.fft(spant[a])))
            s2 = np.abs(np.fft.fftshift(np.fft.fft(fida[b])))
            r = np.corrcoef(s1[win], s2[win])[0, 1]
            assert r > 0.99, f"{a}: r = {r:.3f} against FID-A"

    def test_docker_runtime_matches_local(self, monkeypatch):
        """The Docker image and the local R give the same numbers (image
        built on first use elsewhere; here only if it already exists)."""
        import subprocess
        from basisremy.core.octave_manager import docker_disabled
        if docker_disabled():
            pytest.skip('Docker disabled by BASISREMY_NO_DOCKER')
        try:
            have = subprocess.run(['docker', 'image', 'inspect', spant_manager._DOCKER_IMAGE],
                                  capture_output=True)
        except FileNotFoundError:
            pytest.skip('no docker binary')
        have = subprocess.run(['docker', 'image', 'inspect', spant_manager._DOCKER_IMAGE],
                              capture_output=True)
        if have.returncode != 0:
            pytest.skip("spant Docker image not built on this machine")
        local = SpantBackend().run_simulation({**_BASE, 'Sequence': 'PRESS'})
        monkeypatch.setenv('BASISREMY_SPANT_RUNTIME', 'docker')
        docker = SpantBackend().run_simulation({**_BASE, 'Sequence': 'PRESS'})
        assert np.max(np.abs(docker['naa'] - local['naa'])) < 1e-9

    @pytest.mark.parametrize("seq,extra", [
        ('STEAM', {'TE': 20, 'TM': 10, 'STEAM Variant': 'Coherence filter'}),
        ('sLASER', {'TE': 30}),
        ('Spin Echo', {'TE': 30}),
        ('Pulse-acquire', {}),
    ])
    def test_other_sequences_run(self, seq, extra):
        b = SpantBackend()
        result = b.run_simulation({**_BASE, 'Sequence': seq, **extra})
        assert 'naa' in result and np.max(np.abs(result['naa'])) > 0
        assert abs(self._peak_ppm(result['naa'], 2000, 3.0) - 2.01) < 0.05

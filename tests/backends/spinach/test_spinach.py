####################################################################################################
#                                          test_spinach.py                                         #
####################################################################################################
#                                                                                                  #
# Purpose: Tests for the Spinach backend. The unit part (registration, schema, argument            #
#          building, REMY parsing, the Octave patch file) runs everywhere; the live part needs an  #
#          Octave runtime (Docker or local) and network for the one-time sparse fetch of           #
#          Spinach, and checks the physics against FidaIdeal on the same FID-A spin systems.       #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import os
import shutil
import subprocess
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

from basisremy.backends.spinach_backend import (SpinachBackend, SpinachPressShaped,
                                                SpinachSemiLaserShaped)
from basisremy.core import externals
from basisremy.core.paths import ADAPTERS_DIR


def _octave_runtime():
    """True if Docker (any engine) or a local Octave can run the backend."""
    from basisremy.core.octave_manager import docker_disabled
    if shutil.which('octave-cli') or shutil.which('octave'):
        return True
    if docker_disabled():
        return False
    try:
        import docker
        docker.from_env().ping()
        return True
    except Exception:
        pass
    sock = os.path.expanduser('~/.orbstack/run/docker.sock')
    return os.path.exists(sock)


def _spinach_fetchable():
    """Skip the live tests only when Spinach cannot be downloaded (no network);
    a patch that does not fit the checkout is a real failure and must show."""
    try:
        externals.ensure('spinach')
        return True
    except externals.ExternalFetchError as e:
        return 'patch' in str(e)


# ------------------------------------------------------------------ unit
class TestSpinachSchema:

    def test_registered(self):
        from basisremy.core.basisremy import BasisREMY
        br = BasisREMY()
        assert 'Spinach' in br.backends
        assert br.backends['Spinach'].category == 'Spinach'
        assert br.categories['Spinach'] == ['Spinach', 'SpinachPressShaped', 'SpinachSemiLaserShaped']

    def test_schema(self):
        b = SpinachBackend()
        assert b.requires_octave
        assert b.dropdown['Sequence'] == ['Spin Echo', 'PRESS', 'STEAM', 'LASER']
        for k in ('Sequence', 'Samples', 'Bandwidth', 'Bfield', 'Linewidth', 'TE', 'TM', 'Metabolites'):
            assert k in b.mandatory_params
        # FID-A's spin systems, FID-A's default selection
        assert 'NAA' in b.metabs and 'Ins' in b.metabs and b.metabs['NAA']
        assert 'NAA' in b.mandatory_params['Metabolites']

    def test_tm_only_for_steam(self):
        b = SpinachBackend()
        b.mandatory_params['Sequence'] = 'PRESS'
        assert 'TM' not in b.get_params_for_mode()
        b.mandatory_params['Sequence'] = 'STEAM'
        assert 'TM' in b.get_params_for_mode()

    def test_build_args(self):
        b = SpinachBackend()
        base = {'Samples': 2048, 'Bandwidth': 2000, 'Bfield': 3, 'Linewidth': 2, 'TE': 35, 'TM': 12}
        assert b._build_args({**base, 'Sequence': 'PRESS'}) == ['press', 2048.0, 2000.0, 3.0, 2.0, 35.0, 0.0]
        assert b._build_args({**base, 'Sequence': 'STEAM'})[-1] == 12.0
        assert b._build_args({**base, 'Sequence': 'Spin Echo'})[0] == 'spinecho'
        assert b._build_args({**base, 'Sequence': 'LASER'})[0] == 'laser'

    def test_missing_inputs_raise(self):
        b = SpinachBackend()
        with pytest.raises(ValueError, match="'TE'"):
            b._build_args({'Sequence': 'PRESS', 'Samples': 2048, 'Bandwidth': 2000, 'Bfield': 3, 'TE': ''})
        with pytest.raises(ValueError, match='Sequence'):
            b._build_args({'Sequence': 'MEGA-PRESS', 'Samples': 2048, 'Bandwidth': 2000, 'Bfield': 3, 'TE': 35})

    def test_parse_remy_and_protocol(self):
        b = SpinachBackend()
        mand, opt = b.parseREMY({'NumberOfDatapoints': 4096, 'SpectralWidth': 4000, 'B0': 2.89,
                                 'TE': 30, 'Protocol': 'svs_se_30', 'Nucleus': '1H'})
        assert mand == {'Samples': 4096, 'Bandwidth': 4000, 'Bfield': 2.89, 'TE': 30, 'Sequence': 'Spin Echo'}
        assert opt['Nucleus'] == '1H'
        assert b.parseProtocol('PRESS_35') == 'PRESS'
        assert b.parseProtocol('steam_te20') == 'STEAM'
        assert b.parseProtocol('sLASER') == 'LASER'
        assert b.parseProtocol('') is None
        assert b.map_sequence_in('MEGA-PRESS') == 'PRESS'

    def test_patch_and_shims_ship_with_the_package(self):
        patch = ADAPTERS_DIR / 'backends' / 'spinach_octave.patch'
        shims = ADAPTERS_DIR / 'backends' / 'spinach_shims'
        assert patch.is_file() and patch.stat().st_size > 1000
        for f in ('isworkernode.m', 'pad.m', 'contains.m', 'allfinite.m', 'md5_hash.m', 'gather.m'):
            assert (shims / f).is_file(), f
        assert externals.PATCHES['spinach'].endswith('spinach_octave.patch')
        assert externals.SPARSE['spinach'] == ['kernel']


class TestSpinachShapedSchema:

    def test_schemas(self):
        p = SpinachPressShaped()
        assert p.category == 'Spinach' and p.file_selection == ['Path to Pulse']
        for k in ('TE', 'Tau 1', 'Tau 2', 'RefTp', 'thkX', 'thkY', 'fovX', 'fovY', 'nX', 'nY',
                  'Flip Angle', 'Sim Centre (ppm)', 'Metabolites'):
            assert k in p.mandatory_params, k
        s = SpinachSemiLaserShaped()
        assert 'TE' in s.mandatory_params and 'Tau 1' not in s.mandatory_params
        assert s.parseProtocol('svs_slaser') == 'sLASER' and p.parseProtocol('PRESS_35') == 'PRESS'

    def test_build_args(self, tmp_path):
        pulse = tmp_path / 'p.pta'
        pulse.write_text('x')
        base = {'Samples': 1024, 'Bandwidth': 2000, 'Bfield': 3, 'Linewidth': 1, 'TE': 30,
                'RefTp': 5, 'thkX': 2, 'thkY': 2, 'fovX': 1.5, 'fovY': 1.5, 'nX': 2, 'nY': 2,
                'Flip Angle': 180, 'Sim Centre (ppm)': 4.65, 'Path to Pulse': str(pulse)}
        a = SpinachPressShaped()._build_args(base)
        assert a[0] == 'press_shaped' and a[5:7] == [15.0, 15.0] and a[7].endswith('p.pta')
        assert a[8:] == [5.0, 2.0, 2.0, 1.5, 1.5, 2, 2, 180.0, 4.65]
        a = SpinachPressShaped()._build_args({**base, 'Tau 1': 12, 'Tau 2': 18})
        assert a[5:7] == [12.0, 18.0]
        b = SpinachSemiLaserShaped()._build_args(base)
        assert b[0] == 'semilaser_shaped' and b[5] == 30.0 and b[6].endswith('p.pta')

    def test_missing_pulse_raises(self):
        with pytest.raises(ValueError, match='Path to Pulse'):
            SpinachSemiLaserShaped()._build_args({'Samples': 1024, 'Bandwidth': 2000,
                                                  'Bfield': 3, 'TE': 30})


# ------------------------------------------------------------------ live
@pytest.mark.skipif(not _octave_runtime(), reason='needs Docker or a local Octave')
@pytest.mark.skipif(not _spinach_fetchable(), reason='Spinach could not be fetched (network)')
class TestSpinachLive:
    _BASE = {'Samples': 2048, 'Bandwidth': 2000, 'Bfield': 3.0, 'Linewidth': 1, 'TE': 35, 'TM': 10}

    @staticmethod
    def _spec(fid):
        return np.abs(np.fft.fftshift(np.fft.fft(fid)))

    @staticmethod
    def _ppm(n=2048, sw=2000, b0=3.0):
        return np.fft.fftshift(np.fft.fftfreq(n, 1 / sw)) / (b0 * 42.577) + 4.65

    def test_patch_is_applied_to_the_checkout(self):
        root = externals.ensure('spinach')
        assert os.path.isfile(os.path.join(root, 'kernel', 'overloads', '@polyadic', 'absorb_prefix.m'))
        assert not os.path.isfile(os.path.join(root, 'kernel', 'overloads', '@polyadic', 'prefix.m'))
        r = subprocess.run(['git', '-C', root, 'apply', '--check', '--reverse',
                            str(ADAPTERS_DIR / 'backends' / 'spinach_octave.patch')],
                           capture_output=True)
        assert r.returncode == 0, r.stderr.decode()[:500]

    @pytest.mark.parametrize('seq', ['PRESS', 'STEAM', 'Spin Echo', 'LASER'])
    def test_matches_fida_ideal(self, seq):
        from basisremy.backends.fida_backends import FidaIdeal
        metabs = ['NAA', 'Lac', 'Glu']
        params = {**self._BASE, 'Sequence': seq, 'Metabolites': metabs}
        sp = SpinachBackend().run_simulation(dict(params))
        fa = FidaIdeal().run_simulation({**params, 'TE2': 0, 'Center Freq': 127.7})
        ppm = self._ppm()
        for m in metabs:
            a, b = fa[m], sp[m]
            assert a.shape == b.shape == (2048,)
            # same amplitude scale and same first point (FID-A's 2^(2-nspins) normalisation).
            # The residual difference is the proton gyromagnetic ratio: FID-A rounds it to
            # 42.576 MHz/T, Spinach uses the CODATA 42.5775 — a 3.5e-5 relative frequency
            # offset, i.e. ~5 degrees of phase at the end of a 1 s FID for a resonance 3 ppm
            # off water, and ~1e-4 on zero-quantum evolution during STEAM's mixing time.
            assert abs(a[0]) == pytest.approx(abs(b[0]), rel=1e-3)
            # same complex spectrum, not just the magnitude (phase convention included)
            r = np.corrcoef(np.real(np.fft.fft(a)), np.real(np.fft.fft(b)))[0, 1]
            assert r > 0.999, f'{seq} {m}: r={r:.4f}'
            # point by point over the first eighth of the FID (128 ms: gamma phase drift
            # below 0.7 degrees, i.e. ~1% of the amplitude for a resonance 3 ppm off water)
            q = len(a) // 8
            assert np.max(np.abs(a[:q] - b[:q])) < 2e-2 * np.max(np.abs(a)), f'{seq} {m}'
            assert ppm[np.argmax(self._spec(a))] == pytest.approx(ppm[np.argmax(self._spec(b))], abs=1e-9)

    def test_naa_singlet_at_2_01_ppm(self):
        sp = SpinachBackend().run_simulation({**self._BASE, 'Sequence': 'PRESS', 'Metabolites': ['NAA']})
        ppm = self._ppm()
        assert ppm[np.argmax(self._spec(sp['NAA']))] == pytest.approx(2.008, abs=0.02)

    def test_h2o_scale_factor(self):
        # sysH2O: one spin with scaleFactor 2 -> amplitude 2 (per-proton amplitude 1)
        sp = SpinachBackend().run_simulation({**self._BASE, 'Sequence': 'Spin Echo', 'Metabolites': ['H2O']})
        assert abs(sp['H2O'][0]) == pytest.approx(2.0, rel=1e-6)

    @pytest.mark.parametrize('spinach_cls,fida_name,pulse,extra', [
        (SpinachPressShaped, 'FidaPressShaped',
         'externals/fidA/rfPulseTools/rfPulses/sampleRefocPulse.pta', {'TE': 30, 'RefTp': 5.0}),
        (SpinachSemiLaserShaped, 'FidaSemiLaserShaped',
         'externals/jbss/my_pulse/standardized_goia.txt', {'TE': 35, 'RefTp': 4.5008}),
    ], ids=['press_shaped', 'semilaser_goia'])
    def test_shaped_matches_fida(self, spinach_cls, fida_name, pulse, extra):
        import basisremy.backends.fida_backends as fb
        root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..'))
        pulse = os.path.join(root, pulse)
        if not os.path.exists(pulse):
            pytest.skip(f'{pulse} not available')
        metabs = ['NAA', 'Glu']
        params = {'Samples': 1024, 'Bandwidth': 2000, 'Bfield': 3.0, 'Linewidth': 1.0,
                  'thkX': 2.0, 'thkY': 2.0, 'fovX': 1.5, 'fovY': 1.5, 'nX': 2, 'nY': 2,
                  'Flip Angle': 180.0, 'Sim Centre (ppm)': 4.65, 'Path to Pulse': pulse,
                  'Metabolites': metabs, **extra}
        sp = spinach_cls().run_simulation(dict(params))
        fa = getattr(fb, fida_name)().run_simulation(dict(params))
        ppm = self._ppm(n=1024)
        for m in metabs:
            a, b = fa[m], sp[m]
            assert abs(a[0]) == pytest.approx(abs(b[0]), rel=1e-3)
            r = np.corrcoef(np.real(np.fft.fft(a)), np.real(np.fft.fft(b)))[0, 1]
            assert r > 0.999, f'{fida_name} {m}: r={r:.4f}'
            assert ppm[np.argmax(self._spec(a))] == pytest.approx(ppm[np.argmax(self._spec(b))], abs=1e-9)

    def test_every_fida_metabolite_simulates(self):
        b = SpinachBackend()
        metabs = list(b.metabs)                     # all 29 FID-A spin systems, incl. the
        out = b.run_simulation({**self._BASE, 'Sequence': 'PRESS', 'Metabolites': metabs})
        assert list(out) == metabs                  # v7.3-only ones (GSH, EtOH, Ref0ppm)
        for m in metabs:
            assert np.isfinite(out[m]).all() and np.abs(out[m]).max() > 0, m

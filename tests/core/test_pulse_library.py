"""core/pulse_library.py: generated standard pulses, their profiles, writers and resolver."""

import os

import numpy as np
import pytest

from basisremy.core import pulse_library as pl
from basisremy.core.rf_pulses import read_waveform, scale_waveform


@pytest.mark.parametrize('name', list(pl.STANDARD))
def test_standard_pulse_flips_on_resonance(name):
    p = pl.make_standard(name)
    assert p.waveform.shape[0] >= 256 and np.max(p.waveform[:, 1]) == pytest.approx(1.0)
    m = pl.bloch_profile(p, [0.0])
    if p.kind == 'exc':
        assert abs(m[0, 0] + 1j * m[1, 0]) > 0.99
    else:
        assert m[2, 0] < -0.95


def test_sinc_and_hs_bandwidth_follow_the_design():
    assert pl.bandwidth_hz(pl.make_standard('sinc-exc'), span_hz=8000) == pytest.approx(3900, rel=0.1)
    assert pl.bandwidth_hz(pl.make_standard('sinc-ref'), span_hz=5000) == pytest.approx(1120, rel=0.1)
    hs = pl.make_standard('hs4-ref')
    assert pl.bandwidth_hz(hs, span_hz=20000) == pytest.approx(hs.bw_hz, rel=0.05)


def test_gaussian_is_calibrated_to_its_fwhm():
    g = pl.gauss(14.0, 88.0)
    assert pl.bandwidth_hz(g, span_hz=400, points=4001) == pytest.approx(88.0, rel=0.03)
    m = pl.bloch_profile(g, [0.0, 200.0])             # on resonance; nothing 200 Hz away
    assert m[2, 0] < -0.95 and m[2, 1] > 0.9


def test_goia_has_gradient_column_and_edges_full_gradient():
    g = pl.make_standard('goia-wurst')
    assert g.is_gradient_modulated and g.waveform.shape[1] == 4
    grad = g.waveform[:, 3]
    assert grad[0] == pytest.approx(grad.max()) and grad[g.n // 2] == pytest.approx(0.1 * grad.max(), rel=0.05)
    with pytest.raises(ValueError, match='gmax'):
        pl.goia(4.5, 45.0, 0.1)


def test_writers_round_trip_through_the_loader(tmp_path):
    for name, (factory, ext) in pl.STANDARD.items():
        p = factory()
        path = p.write(str(tmp_path / f'{name}{ext}'))
        rf = read_waveform(path)
        assert rf.shape[0] == p.n
        assert np.allclose(rf[:, 1] / rf[:, 1].max(), p.waveform[:, 1], atol=1e-6)
        assert np.allclose(np.mod(rf[:, 0], 360), np.mod(p.waveform[:, 0], 360), atol=0.01)
        if p.is_gradient_modulated:
            assert np.allclose(rf[:, 3], p.waveform[:, 3], atol=1e-6)
    with pytest.raises(ValueError, match='gradient'):
        pl.make_standard('goia-wurst').to_pta(str(tmp_path / 'x.pta'))
    rf_file = pl.make_standard('sinc-ref').to_rf(str(tmp_path / 'sinc.RF'))
    assert read_waveform(rf_file).shape == (256, 3)


def test_fslmrs_and_array_exports():
    p = pl.make_standard('sinc-ref')
    block = p.to_fslmrs()
    assert block['time'] == pytest.approx(0.005) and len(block['amp']) == p.n
    assert max(block['amp']) == pytest.approx(scale_waveform(p.waveform, 0.005, 'ref')[1].max())
    assert block['grad'] == [0.0, 0.0, 0.0]
    g = pl.make_standard('goia-wurst').to_fslmrs(grad_mt_per_m=12.0)
    assert len(g['grad']) == 3 and len(g['grad'][0]) == 512 and max(g['grad'][0]) == pytest.approx(12.0)
    b1, dt = p.to_array()
    assert b1.dtype == complex and dt.sum() == pytest.approx(0.005)


def test_resolver(tmp_path):
    path = pl.resolve_pulse('standard:goia-wurst', str(tmp_path))
    assert path.endswith('std_goia-wurst.txt') and os.path.exists(path)
    assert pl.resolve_pulse('standard:goia-wurst', str(tmp_path)) == path      # cached
    assert pl.resolve_pulse('/some/file.pta', str(tmp_path)) == '/some/file.pta'
    with pytest.raises(ValueError, match='Unknown standard pulse'):
        pl.resolve_pulse('standard:nope', str(tmp_path))


@pytest.mark.parametrize('cls_path', [
    'basisremy.backends.fida_backends.FidaPressShaped',
    'basisremy.backends.fida_backends.FidaSemiLaserShaped',
    'basisremy.backends.fida_backends.FidaSteamShaped',
    'basisremy.backends.fida_backends.FidaSpinEchoShaped',
    'basisremy.backends.fida_backends.FidaMegaPressShaped',
    'basisremy.backends.fida_backends.FidaMegaSpecialShaped',
    'basisremy.backends.fida_backends.FidaOnePulse',
    'basisremy.backends.spinach_backend.SpinachPressShaped',
    'basisremy.backends.spinach_backend.SpinachSemiLaserShaped',
    'basisremy.backends.vespa_backend.VespaBackend',
    'basisremy.backends.spant_backend.SpantBackend',
    'basisremy.backends.custom_backends.CustomSLaser',
])
def test_real_pulse_backends_start_with_open_pulse(cls_path):
    """Every pulse field starts with a catalogue pulse (replaced by a vendor file)."""
    import importlib
    mod, cls = cls_path.rsplit('.', 1)
    b = getattr(importlib.import_module(mod), cls)()
    keys = [k for k in b.file_selection if k in ('Path to Pulse', 'Edit Pulse Path')]
    assert keys
    for k in keys:
        v = b.mandatory_params[k]
        assert pl.is_standard(v) and pl.standard_name(v) in pl.STANDARD, (cls, k, v)

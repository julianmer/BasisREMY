"""RF pulse loading / scaling (port of FID-A's io_loadRFwaveform)."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from basisremy.core.rf_pulses import (  # noqa: E402
    load_pulse, read_waveform, scale_waveform, _bloch_mz_after_pulse)

PULSES = os.path.join(os.path.dirname(__file__), '..', '..', 'externals',
                      'fidA', 'rfPulseTools', 'rfPulses')


def _pulse(name):
    p = os.path.join(PULSES, name)
    if not os.path.exists(p):
        pytest.skip(f"{name} not available (fidA external not fetched)")
    return p


class TestScaling:
    def test_hard_pulse_180(self):
        rf = np.column_stack([np.zeros(100), np.ones(100), np.ones(100)])
        phase, amp_hz, dt = scale_waveform(rf, 1e-3, 'ref')
        # flat pulse: w1max * Tp = 0.5 cycles for a 180
        assert amp_hz.max() * 1e-3 == pytest.approx(0.5)
        assert dt.sum() == pytest.approx(1e-3)
        # Bloch check: inverts
        mz = _bloch_mz_after_pulse(phase, amp_hz / amp_hz.max(), np.ones(100),
                                   1e-3, np.array([amp_hz.max()]))
        assert mz[0] == pytest.approx(-1.0, abs=1e-6)

    def test_hard_pulse_90(self):
        rf = np.column_stack([np.zeros(50), np.ones(50)])
        _, amp_hz, _ = scale_waveform(rf, 2e-3, 'exc')
        assert amp_hz.max() * 2e-3 == pytest.approx(0.25)

    def test_numeric_flip(self):
        rf = np.column_stack([np.zeros(50), np.ones(50)])
        _, amp_hz, _ = scale_waveform(rf, 1e-3, 45.0)
        assert amp_hz.max() * 1e-3 == pytest.approx(0.125)


class TestFiles:
    def test_pta_mao_refocusing(self):
        rf = read_waveform(_pulse('sampleRefocPulse.pta'))
        assert rf.shape == (400, 3)
        assert set(np.round(rf[:, 0]).astype(int)) <= {0, 180}   # AM pulse
        p = load_pulse(_pulse('sampleRefocPulse.pta'), 5.0, 'ref')
        amp = np.array(p['amp_hz'])
        assert len(amp) == 400 and amp.max() > 0
        # scaled AM refocusing pulse inverts on resonance
        mz = _bloch_mz_after_pulse(np.array(p['phase_deg']), amp / amp.max(),
                                   np.ones(400), 5e-3, np.array([amp.max()]))
        assert mz[0] == pytest.approx(-1.0, abs=0.05)

    def test_rf_adiabatic(self):
        p = load_pulse(_pulse('sampleAFPpulse_HS2_R15.RF'), 5.0, 'inv')
        amp = np.array(p['amp_hz'])
        assert amp.max() > 0
        mz = _bloch_mz_after_pulse(np.array(p['phase_deg']), amp / amp.max(),
                                   np.ones(len(amp)), 5e-3, np.array([amp.max()]))
        assert mz[0] < -0.95   # reached the adiabatic plateau

    def test_txt_goia(self):
        rf = read_waveform(_pulse('GOIA_tthk0.01_R120.txt'))
        assert rf.shape[1] == 3 and len(rf) > 10

    def test_unknown_extension(self, tmp_path):
        f = tmp_path / 'p.xyz'
        f.write_text('1 2\n')
        with pytest.raises(ValueError, match='Unrecognised'):
            read_waveform(str(f))

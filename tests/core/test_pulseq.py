####################################################################################################
#                                          test_pulseq.py                                          #
####################################################################################################
#                                                                                                  #
# Purpose: Pulseq (.seq) input. A PRESS sequence with known timing is written with PyPulseq and    #
#          read back: sequence type, TE and the echo split from the RF centres, readout, pulse     #
#          duration and flip angle, slabs from the slice gradients, and the sheet BasisREMY fills. #
#                                                                                                  #
####################################################################################################

import math
import os

import numpy as np
import pytest

pp = pytest.importorskip('pypulseq')

from basisremy.core import pulseq                 # noqa: E402


def _press(path, te1_ms=20.0, te_ms=35.0, samples=2048, dwell=250e-6):
    """PRESS: 90 (slab x 2 cm) - 180 (y, 3 cm) - 180 (z, 4 cm), ADC at the echo."""
    sys = pp.Opts(max_grad=30, grad_unit='mT/m', max_slew=120, slew_unit='T/m/s',
                  rf_ringdown_time=20e-6, rf_dead_time=100e-6, adc_dead_time=10e-6)
    seq = pp.Sequence(sys)
    exc, gx, _ = pp.make_sinc_pulse(flip_angle=math.pi / 2, duration=3e-3, slice_thickness=0.02,
                                    time_bw_product=8, system=sys, return_gz=True, use='excitation')
    gx.channel = 'x'
    blocks = [(exc, gx)]
    for thk, ch in ((0.03, 'y'), (0.04, 'z')):
        ref, g, _ = pp.make_sinc_pulse(flip_angle=math.pi, duration=3e-3, slice_thickness=thk,
                                       time_bw_product=8, phase_offset=math.pi / 2, system=sys,
                                       return_gz=True, use='refocusing')
        g.channel = ch
        blocks.append((ref, g))
    centre = lambda rf: rf.delay + pp.calc_rf_center(rf)[0]           # noqa: E731
    raster = sys.block_duration_raster
    t0 = centre(exc)
    targets = [t0 + te1_ms / 2e3, t0 + te1_ms / 2e3 + te_ms / 2e3]
    seq.add_block(*blocks[0])
    t = pp.calc_duration(*blocks[0])
    for (rf, g), target in zip(blocks[1:], targets):
        wait = round((target - centre(rf) - t) / raster) * raster
        seq.add_block(pp.make_delay(wait))
        seq.add_block(rf, g)
        t += wait + pp.calc_duration(rf, g)
    adc = pp.make_adc(num_samples=samples, dwell=dwell, system=sys)
    seq.add_block(pp.make_delay(round((t0 + te_ms / 1e3 - t) / raster) * raster))
    seq.add_block(adc)
    seq.write(str(path))
    return path


def test_read_press_timing_and_pulses(tmp_path):
    info = pulseq.read_seq(_press(tmp_path / 'press.seq'))
    assert pulseq.sequence_type(info) == 'PRESS'
    assert [round(r.flip_deg) for r in info.rf] == [90, 180, 180]
    assert [r.slab_axis for r in info.rf] == ['x', 'y', 'z']
    assert pulseq.echo_ms(info) == pytest.approx(35.0, abs=0.02)
    assert (info.samples, info.dwell_s) == (2048, 250e-6)
    slabs = [pulseq.slab_cm(r) for r in info.rf[1:]]
    assert slabs[1] / slabs[0] == pytest.approx(4 / 3, rel=0.01)     # same pulse, 3 vs 4 cm slab


def test_sheet_from_seq(tmp_path):
    backend, p = pulseq.sheet_params(str(_press(tmp_path / 'press.seq')), str(tmp_path / 'work'))
    assert backend == 'FidaPressShaped'
    assert p['TE'] == pytest.approx(35.0, abs=0.02)
    assert p['Tau 1'] == pytest.approx(20.0, abs=0.02) and p['Tau 2'] == pytest.approx(15.0, abs=0.02)
    assert (p['Samples'], p['Bandwidth'], p['RefTp']) == (2048, 4000.0, 3.0)
    assert p['Flip Angle'] == pytest.approx(180.0, abs=0.5)
    assert p['fovX'] > p['thkX'] > 0 and p['fovY'] > p['thkY'] > 0
    rf = np.loadtxt(p['Path to Pulse'], comments='#')
    assert rf.shape[1] == 3 and rf[:, 0].max() == pytest.approx(1.0)      # amplitude, phase, step


def test_basisremy_load_sequence(tmp_path):
    from basisremy.core.basisremy import BasisREMY
    br = BasisREMY('FidaIdeal')
    br.load_sequence(str(_press(tmp_path / 'press.seq')))
    assert br.backend.name == 'FidaPressShaped'
    assert br.backend.mandatory_params['Tau 1'] == pytest.approx(20.0, abs=0.02)
    assert os.path.exists(br.backend.mandatory_params['Path to Pulse'])


def test_unsupported_sequence_raises(tmp_path):
    sys = pp.Opts()
    seq = pp.Sequence(sys)
    seq.add_block(pp.make_block_pulse(flip_angle=math.pi / 2, duration=1e-3, system=sys))
    seq.add_block(pp.make_adc(num_samples=64, dwell=1e-4, system=sys))
    seq.write(str(tmp_path / 'fid.seq'))
    with pytest.raises(ValueError, match='not supported'):
        pulseq.sequence_type(pulseq.read_seq(str(tmp_path / 'fid.seq')))

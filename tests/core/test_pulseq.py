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


def _slaser(path, te_ms=40.0, shift_ms=0.0, samples=1024, dwell=250e-6):
    """semi-LASER: sinc 90 (x) + four HS1 adiabatic 180s (y, y, z, z) at TE/8, 3TE/8, 5TE/8,
    7TE/8, played as arbitrary waveforms without a 'use' tag; shift_ms moves the 2nd pulse."""
    from basisremy.core.pulse_library import hsn
    sys = pp.Opts(max_grad=30, grad_unit='mT/m', max_slew=120, slew_unit='T/m/s',
                  rf_ringdown_time=20e-6, rf_dead_time=100e-6, adc_dead_time=10e-6)
    seq = pp.Sequence(sys)
    exc, gx, _ = pp.make_sinc_pulse(flip_angle=math.pi / 2, duration=2.6e-3, slice_thickness=0.02,
                                    time_bw_product=8, system=sys, return_gz=True, use='excitation')
    gx.channel = 'x'
    afp = hsn(4.5, 15.0, order=1, n=450)
    signal = 700.0 * afp.waveform[:, 1] * np.exp(1j * np.deg2rad(afp.waveform[:, 0]))   # Hz
    raster = sys.block_duration_raster
    blocks = [(exc, gx)]
    for thk, ch in ((0.03, 'y'), (0.03, 'y'), (0.04, 'z'), (0.04, 'z')):
        g = pp.make_trapezoid(channel=ch, amplitude=15.0 / 4.5e-3 / thk,
                              flat_time=math.ceil(4.5e-3 / raster) * raster, system=sys)
        rf = pp.make_arbitrary_rf(signal, flip_angle=math.pi, dwell=1e-5, no_signal_scaling=True,
                                  delay=g.rise_time, system=sys)
        blocks.append((rf, g))
    centre = lambda rf: rf.delay + pp.calc_rf_center(rf)[0]           # noqa: E731
    t0 = centre(exc)
    targets = [t0 + k * te_ms / 8e3 for k in (1, 3, 5, 7)]
    targets[1] += shift_ms / 1e3
    seq.add_block(*blocks[0])
    t = pp.calc_duration(*blocks[0])
    for (rf, g), target in zip(blocks[1:], targets):
        wait = round((target - centre(rf) - t) / raster) * raster
        seq.add_block(pp.make_delay(wait))
        seq.add_block(rf, g)
        t += wait + pp.calc_duration(rf, g)
    seq.add_block(pp.make_delay(round((t0 + te_ms / 1e3 - t) / raster) * raster))
    seq.add_block(pp.make_adc(num_samples=samples, dwell=dwell, system=sys))
    seq.write(str(path))
    return path


def test_semilaser_from_seq(tmp_path):
    info = pulseq.read_seq(str(_slaser(tmp_path / 'slaser.seq')))
    assert pulseq.sequence_type(info) == 'sLASER'
    assert [r.role for r in info.rf] == ['exc', 'ref', 'ref', 'ref', 'ref']
    assert [r.adiabatic for r in info.rf] == [False, True, True, True, True]
    assert [r.slab_axis for r in info.rf] == ['x', 'y', 'y', 'z', 'z']
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('error', UserWarning)                     # FID-A spacing: no warning
        backend, p = pulseq.sheet_params(str(tmp_path / 'slaser.seq'), str(tmp_path / 'work'))
    assert backend == 'FidaSemiLaserShaped'
    assert p['TE'] == pytest.approx(40.0, abs=0.02)
    assert (p['RefTp'], p['Flip Angle'], p['Samples'], p['Bandwidth']) == (4.5, 180.0, 1024, 4000.0)
    assert p['thkY'] / p['thkX'] == pytest.approx(4 / 3, rel=0.01)
    assert 'Tau 1' not in p


def test_semilaser_off_spacing_warns(tmp_path):
    path = str(_slaser(tmp_path / 'shifted.seq', shift_ms=1.0))
    with pytest.warns(UserWarning, match='spacing'):
        pulseq.sheet_params(path, str(tmp_path / 'work'))


def _tutorial(example_data_dir, name):
    """Pulseq tutorial file (Git LFS); skip when absent or still an LFS pointer."""
    path = os.path.join(example_data_dir, 'pulseq_tutorial', name)
    if not os.path.exists(path):
        pytest.skip(f'{name} not found')
    with open(path, 'rb') as f:
        if f.read(40).startswith(b'version https://git-lfs'):
            pytest.skip(f'{name} is a Git LFS pointer (run git lfs pull)')
    return path


def test_tutorial_press_echo_matches_the_measurement(example_data_dir):
    """Pulseq tutorial PRESS (CC0): TE and the echo split from the .seq file, and the echo
    position they predict in the readout against the echo in the TWIX of the same scan."""
    info = pulseq.read_seq(_tutorial(example_data_dir, '06_PRESS_center.seq'))
    assert pulseq.sequence_type(info) == 'PRESS'
    assert pulseq.echo_ms(info) == pytest.approx(119.98, abs=0.01)
    exc = info.rf[0]
    predicted = (exc.centre_ms + pulseq.echo_ms(info) - info.adc_start_ms) / (info.dwell_s * 1e3)
    mapvbvd = pytest.importorskip('mapvbvd')
    tw = mapvbvd.mapVBVD(_tutorial(example_data_dir, '06_PRESS_center.dat'), quiet=True)
    tw = tw[-1] if isinstance(tw, list) else tw
    tw.image.squeeze = True
    tw.image.flagRemoveOS = False                                 # keep the 4096 ADC samples
    data = np.asarray(tw.image[''])
    assert data.shape[0] == info.samples
    measured = int(np.argmax(np.abs(data.reshape(info.samples, -1)).sum(axis=1)))
    assert abs(measured - predicted) < 4                          # 800.8 vs 803 (fat/water phantom)


def test_unsupported_sequence_raises(tmp_path):
    sys = pp.Opts()
    seq = pp.Sequence(sys)
    seq.add_block(pp.make_block_pulse(flip_angle=math.pi / 2, duration=1e-3, system=sys))
    seq.add_block(pp.make_adc(num_samples=64, dwell=1e-4, system=sys))
    seq.write(str(tmp_path / 'fid.seq'))
    with pytest.raises(ValueError, match='not supported'):
        pulseq.sequence_type(pulseq.read_seq(str(tmp_path / 'fid.seq')))

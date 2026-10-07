####################################################################################################
#                                     test_sequence_view.py                                        #
####################################################################################################
#                                                                                                  #
# Purpose: Tests for core/sequence_view.py: one pulse reader for every source (standard, .pta,     #
#          FSL-MRS sequence / basis JSON by role), roles from the coherence filter (MEGA editing,  #
#          STEAM), the exact timeline of an FSL-MRS description (echo after the last delay) and    #
#          the sheet timelines (PRESS echo split, sLASER TE1/TE2/TE3, FSL-MRS's own MEGA timing).  #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import json
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from basisremy.core import sequence_view as sv
from basisremy.core.pulse_library import make_standard


def _block(ms, f=0.0, n=8):
    return {'time': ms / 1e3, 'frequencyOffset': f, 'phaseOffset': 0,
            'amp': list(np.hanning(n) * 1000 + 1), 'phase': [0.0] * n, 'grad': [0, 0, 0]}


@pytest.fixture
def press_json(tmp_path):
    # WIN-like PRESS: exc 2.6, ref 5.2, ref 5.2 ms; delays put TE1 = 11, TE2 = 24 ms
    seq = {'RF': [_block(2.6), _block(5.2), _block(5.2)], 'delays': [0.0016, 0.0123, 0.0094],
           'CoherenceFilter': [-1, 1, -1], 'B0': 2.89}
    path = tmp_path / 'NAA.json'
    path.write_text(json.dumps({'seq': seq, 'basis': {}}))   # a basis JSON carries it under 'seq'
    return str(path)


def test_fsl_timeline_exact(press_json):
    t = sv.timeline(_Backend('FSL-MRS', {}), seq_file=press_json)
    assert t['exact'] and [e['role'] for e in t['events']] == ['exc', 'ref', 'ref']
    assert t['echo_ms'] == pytest.approx(35.0)
    assert t['events'][1]['centre_ms'] == pytest.approx(5.5, abs=0.05)     # TE1 / 2
    assert t['events'][2]['centre_ms'] == pytest.approx(23.0, abs=0.05)    # TE1 + TE2 / 2


def test_read_pulse_sources(press_json, tmp_path):
    p = sv.read_pulse(press_json, 'ref')
    assert p.kind == 'ref' and p.tp_ms == pytest.approx(5.2) and p.n == 8
    assert sv.read_pulse(press_json, 'exc').tp_ms == pytest.approx(2.6)
    pta = make_standard('sinc-ref').write(str(tmp_path / 'ref.pta'))
    assert sv.read_pulse(pta).n == make_standard('sinc-ref').n
    assert sv.read_pulse('standard:hs4-ref').name.startswith('hs')
    with pytest.raises(ValueError, match='edit'):
        sv.read_pulse(press_json, 'edit')
    with pytest.raises(ValueError, match='Unrecognised'):
        sv.read_pulse(str(tmp_path / 'pulse.xyz'))


def test_roles_from_coherence_filter():
    mega = {'RF': [_block(2.5), _block(11.8), _block(8.5)], 'delays': [0, 0, 0],
            'CoherenceFilter': [-1, None, 1]}
    steam = {'RF': [_block(2.0)] * 3, 'delays': [0, 0, 0], 'CoherenceFilter': [1, 0, -1]}
    offset = {'RF': [_block(2.5, -505), _block(11.8, -817), _block(8.5, -505)], 'delays': [0, 0, 0]}
    assert sv._fsl_roles(mega) == ['exc', 'edit', 'ref']
    assert sv._fsl_roles(steam) == ['exc', 'exc', 'exc']
    assert sv._fsl_roles(offset) == ['exc', 'edit', 'ref']


class _Backend:
    def __init__(self, name, params, mode='Simple'):
        self.name, self.mandatory_params, self.optional_params = name, params, {}
        self.current_mode = mode

    def get_params_for_mode(self, mode=None):
        return dict(self.mandatory_params)


def test_sheet_timelines():
    press = sv.timeline(_Backend('Spant', {'Sequence': 'PRESS', 'TE': 30, 'Tau 1': 10, 'Tau 2': 20}))
    assert [e['centre_ms'] for e in press['events']] == [0.0, 5.0, 20.0] and press['echo_ms'] == 30
    slaser = sv.timeline(_Backend('Spant', {'Sequence': 'sLASER', 'TE': 30, 'sLASER TE1': 8,
                                            'sLASER TE2': 12, 'sLASER TE3': 10}))
    assert [e['centre_ms'] for e in slaser['events']] == [0.0, 4.0, 11.0, 17.0, 25.0]
    steam = sv.timeline(_Backend('Vespa', {'Sequence': 'STEAM', 'TE': 20, 'TM': 10}))
    assert [e['role'] for e in steam['events']] == ['exc'] * 3 and steam['echo_ms'] == 30
    assert sv.timeline(_Backend('Vespa', {'Sequence': 'PRESS', 'TE': None}))['events'] == []
    assert sv.sequence_kind(_Backend('FidaSemiLaserShaped', {})) == 'sLASER'
    assert sv.sequence_kind(_Backend('MRSCloud', {'Sequence': 'MEGA', 'Localization': 'sLASER'})) == 'MEGA-sLASER'


def test_fslmrs_timeline_is_its_own_sequence():
    from basisremy.core.basisremy import BasisREMY
    b = BasisREMY('FSL-MRS').backend
    b.mandatory_params.update({'Sequence': 'MEGA-PRESS', 'TE': 68, 'Bfield': 3.0, 'Samples': 2048,
                               'Bandwidth': 2000, 'Center Freq': 127.7})
    t = sv.timeline(b)
    assert [e['role'] for e in t['events']] == ['exc', 'ref', 'edit', 'ref', 'edit']
    assert t['echo_ms'] == pytest.approx(68.0, abs=0.01)


def test_switch_targets():
    from basisremy.core.basisremy import BasisREMY
    br = BasisREMY('Vespa')
    br.backend.mandatory_params['Sequence'] = 'PRESS shaped'
    assert sv.switch_target(br.backend, 'ref', True) == {'Sequence': 'PRESS'}
    br.backend.mandatory_params['Sequence'] = 'PRESS'
    assert sv.switch_target(br.backend, 'ref', False) == {'Sequence': 'PRESS shaped'}
    br.set_backend('FidaSteamShaped')
    assert sv.pulse_role(br.backend) == 'exc' and sv.pulse_key(br.backend, 'exc') == 'Path to Pulse'
    assert sv.switch_target(br.backend, 'exc', True) == {'backend': 'FidaIdeal', 'Sequence': 'STEAM'}
    br.set_backend('FidaIdeal')
    br.backend.mandatory_params['Sequence'] = 'PRESS'
    assert sv.switch_target(br.backend, 'ref', False) == {'backend': 'FidaPressShaped'}
    assert sv.switch_target(br.backend, 'exc', False) is None
    br.set_backend('FidaMegaPressShaped')
    br.backend.set_mode('Full shaped (refoc + edit)')
    assert sv.switch_target(br.backend, 'edit', True) == {'mode': 'Refoc-only shaped (ideal edit)'}
    br.backend.set_mode('Edit-only shaped (ideal refoc)')
    assert sv.switch_target(br.backend, 'edit', True) == {'backend': 'FidaMegaPressIdeal'}
    assert sv.switch_target(br.backend, 'ref', False) == {'mode': 'Full shaped (refoc + edit)'}
    br.set_backend('FidaSemiLaserShaped')
    assert sv.switch_target(br.backend, 'ref', True) is None    # no ideal sLASER in FID-A


def _bruker(path, n=64):
    amp = 100 * np.hanning(n)
    lines = ['##TITLE= test shape', '##JCAMP-DX= 5.00 Bruker JCAMP library', '##$SHAPE_EXMODE= Refocusing',
             '##$SHAPE_TOTROT= 1.800000e+02', f'##NPOINTS= {n}', '##XYPOINTS= (XY..XY)']
    lines += [f'{a:.6e}, {180.0 if i % 9 == 0 else 0.0:.6e}' for i, a in enumerate(amp)] + ['##END=']
    path.write_text('\n'.join(lines) + '\n')
    return str(path), amp


def test_bruker_and_mat_readers(tmp_path):
    from scipy.io import savemat
    from basisremy.core.rf_pulses import read_waveform
    path, amp = _bruker(tmp_path / 'shape.rfc')
    w = read_waveform(path)
    assert w.shape == (64, 3) and np.allclose(w[:, 1], amp, atol=1e-4) and w[0, 0] == 180.0
    wave = np.column_stack([np.zeros(10), np.linspace(0, 1, 10), np.ones(10), np.ones(10)])
    savemat(str(tmp_path / 'goia.mat'), {'Sweep2': {'waveform': wave, 'tbw': 45.0}})
    assert np.allclose(read_waveform(str(tmp_path / 'goia.mat')), wave)
    savemat(str(tmp_path / 'other.mat'), {'x': np.ones(3)})
    with pytest.raises(ValueError, match='waveform'):
        read_waveform(str(tmp_path / 'other.mat'))


def test_resolve_pulse_converts_every_format(press_json, tmp_path):
    from basisremy.core.pulse_library import resolve_pulse
    from basisremy.core.rf_pulses import read_waveform
    work = tmp_path / 'work'
    work.mkdir()
    native = make_standard('sinc-ref').write(str(tmp_path / 'ref.pta'))
    assert resolve_pulse(native, str(work)) == native                 # read as it is
    out = resolve_pulse(press_json + '#exc', str(work))               # one pulse of a sequence
    assert out.endswith('NAA_exc.pta') and read_waveform(out).shape[0] == 8
    assert resolve_pulse(press_json, str(work)).endswith('NAA.pta')   # default: refocusing
    bruker, amp = _bruker(tmp_path / 'shape.exc')
    out = resolve_pulse(bruker, str(work))
    assert out.endswith('.pta') and np.allclose(read_waveform(out)[:, 1], amp / amp.max(), atol=1e-6)
    assert sv.read_pulse(press_json + '#2').tp_ms == pytest.approx(5.2)
    with pytest.raises(ValueError, match='selector'):
        sv.read_pulse(press_json + '#nope')

####################################################################################################
#                                     test_sequence_design.py                                      #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 07/10/26                                                                                #
#                                                                                                  #
# Purpose: The sequence designer's designs: written as Pulseq .seq and read back unchanged (type,  #
#          TE, timings, pulse roles and durations, editing ppm, pulse bandwidth), the FSL-MRS      #
#          translation of an ideal PRESS equals FSL-MRS's own ideal PRESS, and each engine is      #
#          told whether it runs a design (and its sheet filled with the design's values).          #
#                                                                                                  #
####################################################################################################

import pytest

pytest.importorskip('pypulseq')

from basisremy.core import pulseq                           # noqa: E402
from basisremy.core import sequence_design as sd            # noqa: E402
from basisremy.core.basisremy import BasisREMY              # noqa: E402
from basisremy.core.pulse_library import bandwidth_hz, make_standard   # noqa: E402

_SHEET = {'Samples': 2048, 'Bandwidth': 4000, 'Bfield': 2.89}


@pytest.fixture(scope='module')
def br():
    return BasisREMY()


def _shaped(kind, te):
    d = sd.recommend(kind, te, _SHEET)
    d.pulses['exc'] = {'source': 'standard:sinc-exc', 'dur': 2.0}
    if 'ref' in d.pulses:
        laser = 'LASER' in kind
        d.pulses['ref'] = {'source': 'standard:hs4-ref' if laser else 'standard:sinc-ref',
                           'dur': 3.5 if laser else 5.0}
    return d


@pytest.mark.parametrize('kind', sd.DESIGNABLE)
@pytest.mark.parametrize('shaped', [False, True])
def test_round_trip(kind, shaped, tmp_path):
    edited = 'edit' in sd.roles(kind)
    te = (100.0 if shaped and edited and 'LASER' in kind else 80.0 if shaped and 'LASER' in kind
          else 80.0 if edited else 35.0)
    d = _shaped(kind, te) if shaped else sd.recommend(kind, te, _SHEET)
    assert sd.problems(d) == []
    back = sd.read_design(sd.write_seq(d, str(tmp_path / 'd.seq')))
    assert back.kind == kind
    assert back.te == pytest.approx(te, abs=1e-3)
    for k, v in d.timing.items():
        assert back.timing[k] == pytest.approx(v, abs=1e-3)
    for role, p in d.pulses.items():
        if p['source'] == 'ideal':
            assert back.pulses[role]['source'] == 'ideal'
        else:
            assert back.pulses[role]['source'].endswith(f'#{role}')
            assert back.pulses[role]['dur'] == pytest.approx(p['dur'], abs=1e-3)
    if kind.startswith('MEGA'):
        assert back.edit == (1.9, 7.5)
    assert back.scheme == d.scheme


def test_waveform_and_slab_survive(tmp_path):
    d = _shaped('PRESS', 35.0)
    path = sd.write_seq(d, str(tmp_path / 'p.seq'))
    ref = next(r for r in pulseq.read_seq(path).rf if r.role == 'ref')
    written = pulseq.waveform_pulse(ref, 'ref')
    assert bandwidth_hz(written, 5.0, points=8001) == pytest.approx(
        bandwidth_hz(make_standard('sinc-ref'), 5.0, points=8001), abs=10)
    assert pulseq.slab_cm(ref) == pytest.approx(2.0, abs=0.05)


def test_overlap_and_sum_are_reported():
    d = sd.recommend('PRESS', 35.0, _SHEET)
    d.timing = {'TE1': 10.0, 'TE2': 20.0}
    assert any('must equal TE' in p for p in sd.problems(d))
    d = _shaped('MEGA-sLASER', 35.0)
    assert any('overlap' in p for p in sd.problems(d))


def test_fsl_translation_matches_fsl_ideal_press(tmp_path):
    from basisremy.backends.fslmrs_backend import FSLMRSBackend
    params = {'Sequence': 'PRESS', 'TE': 35, 'Bandwidth': 4000, 'Samples': 2048, 'Bfield': 2.89,
              'Linewidth': 1}
    b = FSLMRSBackend()
    own = b._generate_sequence_json(b._coerce_params(dict(params)))
    path = sd.write_seq(sd.recommend('PRESS', 35.0, _SHEET), str(tmp_path / 'i.seq'))
    tr = sd.fsl_sequences(path, params)[None]
    assert tr['CoherenceFilter'] == own['CoherenceFilter']
    # FSL-MRS's ideal pulses last 10 us, the .seq's 20 us: the gaps differ by that much
    assert tr['delays'] == pytest.approx(own['delays'], abs=2e-5)
    for a, c in zip(own['RF'], tr['RF']):
        flip = lambda r: sum(r['amp']) * r['time'] / len(r['amp'])   # noqa: E731
        assert flip(c) == pytest.approx(flip(a), rel=1e-6)
        assert c['phase'][0] == pytest.approx(a['phase'][0], abs=1e-4)


def test_fsl_translation_mega_on_off(tmp_path):
    params = {'TE': 68, 'Bandwidth': 4000, 'Samples': 2048, 'Bfield': 2.89}
    path = sd.write_seq(sd.recommend('MEGA-PRESS', 68.0, _SHEET), str(tmp_path / 'm.seq'))
    out = sd.fsl_sequences(path, params)
    assert set(out) == {'ON', 'OFF'}
    assert out['ON']['CoherenceFilter'] == [-1, 1, None, -1, None]


def test_plans(br):
    ideal = sd.recommend('PRESS', 35.0, _SHEET)
    shaped_ref = sd.recommend('PRESS', 35.0, _SHEET)
    shaped_ref.pulses['ref'] = {'source': 'standard:sinc-ref', 'dur': 5.0}
    assert sd.plan(ideal, 'MRSCloud', br).status == 'no'
    assert sd.plan(ideal, 'FSL-MRS', br).status == 'ok'
    assert sd.plan(ideal, 'Vespa', br).status == 'ok'
    assert sd.plan(shaped_ref, 'Vespa', br).status == 'ok'
    assert sd.plan(shaped_ref, 'FID-A', br).status == 'ok'
    # a waveform excitation: Vespa runs it ideal and says so
    shaped_ref.pulses['exc'] = {'source': 'standard:sinc-exc', 'dur': 2.0}
    pl = sd.plan(shaped_ref, 'Vespa', br)
    assert pl.status == 'approx' and 'excitation' in pl.notes[0]
    assert sd.plan(sd.recommend('sLASER', 35.0, _SHEET), 'Vespa', br).status == 'no'
    # FID-A's MEGA-PRESS places the pulses its own way
    assert sd.plan(sd.recommend('MEGA-PRESS', 68.0, _SHEET), 'FID-A', br).status == 'approx'


def test_apply_fills_the_engine_sheet(br, tmp_path):
    d = sd.recommend('PRESS', 35.0, _SHEET)
    d.timing = {'TE1': 11.0, 'TE2': 24.0}
    d.pulses['ref'] = {'source': 'standard:sinc-ref', 'dur': 5.0}
    path = sd.write_seq(d, str(tmp_path / 'v.seq'))
    back = sd.read_design(path)
    pl = sd.apply(br, back, path, 'Spant')
    b = br.backend
    assert pl.status == 'ok' and b.name == 'Spant'
    assert b.mandatory_params['Sequence'] == 'PRESS shaped'
    p = {**b.optional_params, **b.mandatory_params}
    assert p['Path to Pulse'] == f'{path}#ref'
    assert float(p['Tau 1']) == pytest.approx(11.0) and float(p['Tau 2']) == pytest.approx(24.0)
    pl = sd.apply(br, back, path, 'FSL-MRS')
    assert br.backend.current_mode == 'Custom'
    assert br.backend.optional_params['Custom Sequence'] == path


_HEADER = {'LeftRightSize': 30.0, 'AnteriorPosteriorSize': 25.0, 'CranioCaudalSize': 20.0}


def test_voxel_from_the_header_reaches_the_slabs(br, tmp_path):
    d = sd.recommend('PRESS', 35.0, _SHEET, _HEADER)
    assert d.voxel == (3.0, 2.5, 2.0) and 'voxel' not in d.rec
    assert sd.recommend('PRESS', 35.0, _SHEET).rec['voxel']          # no header: recommended
    d.pulses['exc'] = {'source': 'standard:sinc-exc', 'dur': 2.0}
    d.pulses['ref'] = {'source': 'standard:sinc-ref', 'dur': 5.0}
    path = sd.write_seq(d, str(tmp_path / 'v.seq'))
    rf = [r for r in pulseq.read_seq(path).rf]
    assert [round(pulseq.slab_cm(r), 1) for r in rf] == [3.0, 2.5, 2.0]
    back = sd.read_design(path)
    assert back.voxel == (3.0, 2.5, 2.0)
    sd.apply(br, back, path, 'FID-A')                 # FID-A's two slabs: the refocusing pulses
    p = {**br.backend.optional_params, **br.backend.mandatory_params}
    assert (p['thkX'], p['thkY']) == (2.5, 2.0)


@pytest.mark.parametrize('kind', ['HERMES', 'HERCULES', 'HERMES (sLASER)', 'HERCULES (sLASER)'])
def test_hadamard_designs(kind, br, tmp_path):
    sheet = {**_SHEET, 'Bfield': 2.89}
    d = sd.recommend(kind, 80.0, sheet)
    family = 'HERMES' if kind.startswith('HERMES') else 'HERCULES'
    assert d.scheme == sd.SCHEMES[family]
    # 20 ms editing pulses; between the pulses of an sLASER pair at TE 80 they get 19.96 ms
    assert d.pulses['edit']['dur'] == (19.96 if 'sLASER' in kind else 20.0)
    assert sd.problems(d) == []
    path = sd.write_seq(d, str(tmp_path / 'h.seq'))
    back = sd.read_design(path)
    assert back.kind == kind and back.scheme == sd.SCHEMES[family]
    assert back.timing == pytest.approx(d.timing, abs=1e-3)
    out = sd.fsl_sequences(path, {'TE': 80, 'Bandwidth': 4000, 'Samples': 2048, 'Bfield': 2.89})
    assert list(out) == ['A', 'B', 'C', 'D']
    assert sd.plan(back, 'FSL-MRS', br).status == 'ok'
    assert sd.plan(back, 'FID-A', br).status == 'no'
    assert sd.plan(back, 'MRSCloud', br).status == 'no'


def test_dual_lobe_edits_both_targets(tmp_path):
    import numpy as np
    d = sd.recommend('HERMES', 80.0, {**_SHEET, 'Bfield': 2.89})
    path = sd.write_seq(d, str(tmp_path / 'h.seq'))
    out = sd.fsl_sequences(path, {'TE': 80, 'Bandwidth': 4000, 'Samples': 2048, 'Bfield': 2.89})
    edit = out['C']['RF'][2]                         # 90 - 180 - edit - 180 - edit
    assert out['C']['CoherenceFilter'][2] is None
    wave = np.asarray(edit['amp']) * np.exp(1j * np.asarray(edit['phase']))
    n = len(wave)
    spec = np.abs(np.fft.fftshift(np.fft.fft(wave, 64 * n)))
    f = np.fft.fftshift(np.fft.fftfreq(64 * n, edit['time'] / n))
    peaks = sorted(f[np.argsort(spec)[-1:]].tolist() + [f[np.argmax(spec * (np.abs(f - f[np.argmax(spec)]) > 100))]])
    want = sorted((p - 4.65) * 2.89 * 42.577 for p in (4.56, 1.90))
    assert peaks == pytest.approx(want, abs=10)
    no_b0 = {k: v for k, v in _SHEET.items() if k != 'Bfield'}
    assert any('field strength' in p for p in sd.problems(sd.recommend('HERMES', 80.0, no_b0)))


def test_fsl_custom_seq_shows_field_and_linewidth(br, tmp_path):
    path = sd.write_seq(sd.recommend('PRESS', 35.0, _SHEET), str(tmp_path / 'f.seq'))
    sd.apply(br, sd.read_design(path), path, 'FSL-MRS')
    shown = br.backend.get_params_for_mode()
    assert {'TE', 'Bfield', 'Linewidth'} <= set(shown)

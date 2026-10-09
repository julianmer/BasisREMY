"""Sequence-design values from data headers (basisremy/remy/design_fields.py): the parsers on
small header dicts, and the values REMY reads from the local example files."""

import contextlib
import io
import math
import os

import numpy as np

import pytest

from basisremy.remy import design_fields as df


def test_ascconv_twix_notation_and_dicom_quotes():
    text = ('### ASCCONV BEGIN object=MrProtDataImpl ###\n'
            'tSequenceFileName = ""%CustomerSeq%\\\\svs_edit_mgs_univ""\n'
            'sWipMemBlock.alFree[7] = 3\n'
            'sWipMemBlock.adFree[8] = 4.58\n'
            'sTXSPEC.aRFPULSE[0].tName = ""ss_rf_exc""\n'
            'sPrepPulses.ucWaterSat = 0x40\n'
            '### ASCCONV END ###')
    p = df.ascconv(text)
    assert p['sWipMemBlock.alFree[7]'] == 3 and p['sWipMemBlock.adFree[8]'] == 4.58
    assert p['sTXSPEC.aRFPULSE[0].tName'] == 'ss_rf_exc'
    assert p['sPrepPulses.ucWaterSat'] == 0x40


def test_flat_twix_matches_ascconv_notation():
    phx = {('sWipMemBlock', 'alFree', '7'): 3.0, ('sTXSPEC', 'aRFPULSE', '0', 'tName'): '"ss_rf_exc"',
           ('tSequenceFileName',): '"%CustomerSeq%\\smm_svs_herc"'}
    flat = df.flat_twix(phx)
    assert set(flat) == {'sWipMemBlock.alFree[7]', 'sTXSPEC.aRFPULSE[0].tName', 'tSequenceFileName'}


def test_siemens_mega_mode_and_other_sequences():
    mega = {'tSequenceFileName': 'svs_edit_mgs_univ', 'sWipMemBlock.adFree[8]': 1.9,
            'sWipMemBlock.adFree[11]': 7.5, 'sWipMemBlock.alFree[12]': 15000}
    assert df.siemens(mega) == {'Edit On': 1.9, 'Edit Off': 7.5, 'Edit Tp': 15.0}
    # svs_edit_859G keeps its editing in other, unconfirmed slots: nothing is read
    assert df.siemens({**mega, 'tSequenceFileName': 'svs_edit_859G'}) == {}


def test_ge_gaba_only():
    hdr = {'rhi_psdname': b'gaba', 'rhi_user20': -356.0, 'rhi_user21': 356.0, 'rhi_user22': 15000.0,
           'rhr_rh_ps_mps_freq': 1277581390}
    assert df.ge(hdr) == {'Edit On': 1.89, 'Edit Off': 7.47, 'Edit Tp': 15.0}
    # oslaser uses the same user CVs for other values
    assert df.ge({**hdr, 'rhi_psdname': b'oslaser', 'rhi_user20': 10000.0, 'rhi_user21': 15000.0}) == {}


def test_bruker_line_wrapped_pulse_struct():
    method = {'$StTM': '10; ', '$VoxPul1Enum': '<Calculated>; ',
              '$VoxPul1': '(0.5, 8400, 90, Yes, 3, 4200, 0.236, 0.200, ; 0, 50, 2.73, <$VoxPul1Shape>); '}
    out = df.bruker(method)
    assert out['TM'] == 10.0
    assert out[df.RF_PULSES] == ['VoxPul1: 0.5 ms, 8400 Hz, 90 deg, Calculated']


def test_nifti_mixing_time_and_edit_pulse():
    hdr = {'MixingTime': 12.0,
           'EditPulse': {'ON': {'PulseOffset': 1.9, 'PulseDuration': 0.015},
                         'OFF': {'PulseOffset': 7.5, 'PulseDuration': 0.015}}}
    assert df.nifti(hdr) == {'TM': 12.0, 'Edit On': 1.9, 'Edit Off': 7.5, 'Edit Tp': 15.0}
    herc = {'EditPulse': {'A': {'PulseOffset': [4.58, 1.9]}, 'C': {'PulseOffset': 4.58}}}
    assert df.nifti(herc) == {df.EDIT_PPM: [4.58, 1.9]}


def test_zero_mixing_time_is_not_a_value():
    assert df.rda({'TM': '0.000000'}) == {}


FILES = {
    'example_data/spec2nii_tests/ge/pFiles/big_gaba/S01_GABA_68.7':
        {'Edit On': 1.89, 'Edit Off': 7.47, 'Edit Tp': 15.0},
    'example_data/BasisREMY_testDatasets/Dataset_00_Bruker_14T_STEAM_08/Dataset_00_Bruker_14T_STEAM_08_method':
        {'TM': 10.0},
    'example_data/spec2nii_tests/siemens/HERCULES/TIEMO/Siemens_TIEMO_HERC.dat':
        {df.EDIT_PPM: [4.58, 1.9, 3.67, 4.18]},
}


@pytest.mark.remy
@pytest.mark.parametrize('path', list(FILES))
def test_remy_reads_design_values_from_example_files(path):
    if not os.path.exists(path):
        pytest.skip('example file not present')
    from basisremy.core.basisremy import BasisREMY
    with contextlib.redirect_stdout(io.StringIO()):
        m = BasisREMY('FidaIdeal').runREMY(import_fpath=path)
    assert {k: m.get(k) for k in FILES[path]} == FILES[path]


def test_bruker_inline_waveform_and_role():
    method = {'$VoxPul1': '(0.5, 8400, 90, No, 3, 4200, 0.2, 0.2, ; 0, 50, 2.7, <$VoxPul1Shape>); ',
              '$VoxPul1Shape': '( 6 ); 0 0 100 0 ; 50 180; '}
    p = df.bruker(method)[df.PULSE_PARAMS]['VoxPul1']
    assert p['role'] == 'exc' and p['dur_ms'] == 0.5 and p['bw_hz'] == 8400
    assert p['waveform'] == [(0.0, 0.0), (100.0, 0.0), (50.0, 180.0)]


def test_designer_takes_the_header_waveform(tmp_path, monkeypatch):
    from basisremy.core import sequence_design as sd
    from basisremy.core.rf_pulses import read_waveform
    monkeypatch.setenv('BASISREMY_SEQUENCES_DIR', str(tmp_path))
    wf = [(100.0 * math.sin(math.pi * (i + 0.5) / 64), 180.0 * (i % 2)) for i in range(64)]
    header = {df.PULSE_PARAMS: {'VoxPul1': {'role': 'exc', 'dur_ms': 0.5, 'bw_hz': 8400, 'flip': 90,
                                            'shape': 'Calculated', 'waveform': wf}}}
    d = sd.recommend('STEAM', 20.0, {'TM': 10.0}, header)
    assert d.pulses['exc']['dur'] == 0.5 and 'pulse:exc' not in d.rec
    w = read_waveform(d.pulses['exc']['source'])
    assert np.allclose(w[:, 1], [a for a, _ in wf], atol=1e-4)
    assert np.allclose(w[:, 0], [ph for _, ph in wf])
    # parameters without a waveform: the pulse stays ideal, the note gives the header's values
    header[df.PULSE_PARAMS]['VoxPul1']['waveform'] = None
    d = sd.recommend('STEAM', 20.0, {'TM': 10.0}, header)
    assert d.pulses['exc']['source'] == 'ideal' and '0.5 ms, 8400 Hz' in d.rec['pulse:exc']


@pytest.mark.remy
def test_designer_uses_the_header_gradient_system(tmp_path, monkeypatch):
    """Dataset_00 (Bruker 14 T STEAM, TE 3 ms) fits only with its own gradients (992 mT/m, 0.27 ms);
    a clinical 80 mT/m / 0.4 ms system misses by 0.06 ms. The .seq keeps them for a read-back."""
    from basisremy.core import sequence_design as sd
    path = 'example_data/BasisREMY_testDatasets/Dataset_00_Bruker_14T_STEAM_08/Dataset_00_Bruker_14T_STEAM_08_method'
    if not os.path.exists(path):
        pytest.skip('example file not present')
    monkeypatch.setenv('BASISREMY_SEQUENCES_DIR', str(tmp_path))
    from basisremy.core.basisremy import BasisREMY
    with contextlib.redirect_stdout(io.StringIO()):
        header = BasisREMY('FidaIdeal').runREMY(import_fpath=path)
    g = header[df.GRADIENTS]
    assert round(g['max_mT_m'], 1) == 992.4 and g['rise_ms'] == 0.27
    d = sd.recommend('STEAM', header['TE'], {'TM': header['TM']}, header)
    seq = sd.write_seq(d, str(tmp_path / 'ds00.seq'))
    assert sd.read_design(seq).gradients == {'max_mT_m': 992.371, 'rise_ms': 0.27}
    d.gradients = None
    with pytest.raises(ValueError, match='0.06 ms'):
        sd.write_seq(d, str(tmp_path / 'clinical.seq'))

"""Sequence-design values from data headers (basisremy/remy/design_fields.py): the parsers on
small header dicts, and the values REMY reads from the local example files."""

import contextlib
import io
import os

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
    'example_data/REMY_tests/Dataset_00_Bruker_14T_STEAM_08/Dataset_00_Bruker_14T_STEAM_08_method':
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

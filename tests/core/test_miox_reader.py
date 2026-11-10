####################################################################################################
#                                       test_miox_reader.py                                        #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 10/11/26                                                                                #
#                                                                                                  #
# Purpose: MIOX opens the data file and BasisREMY's table comes out as with REMY's own readers     #
#          (twix, SPAR, GE P-file, Siemens DICOM VB / XA, RDA); a format REMY could not read       #
#          (Philips DICOM) fills the sheet; without pymiox the REMY readers run.                   #
#                                                                                                  #
####################################################################################################

import glob
import math
import os
import sys

import pytest

from basisremy.core import miox_reader
from basisremy.core.basisremy import BasisREMY

pytest.importorskip('miox')

_FILES = [('BasisREMY_testDatasets/Dataset_28_Siemens_PRESS/*.dat', 'twix'),
          ('BasisREMY_testDatasets/Dataset_13_Philips_SPAR_3T_PRESS_45_Act/*.SPAR', 'spar'),
          ('BasisREMY_testDatasets/Dataset_11_GE_P32256/*.7', 'pfile'),
          ('BasisREMY_testDatasets/Dataset_22_Siemens_Dicom_7T/*.IMA', 'ima'),
          ('spec2nii_tests/siemens/XAData/XA20/DICOM/26516628.dcm', 'xa'),
          ('BasisREMY_testDatasets/Dataset_32_Siemens_RDA_PreOn/*.rda', 'rda')]


def _file(example_data_dir, pattern):
    files = glob.glob(os.path.join(example_data_dir, pattern))
    if not files:
        pytest.skip(f'{pattern} not found')
    return files[0]


def _same(a, b):
    if isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b):
        return True
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return math.isclose(a, b, rel_tol=1e-9)
    return a == b


@pytest.mark.parametrize('pattern', [f for f, _ in _FILES], ids=[i for _, i in _FILES])
def test_miox_table_equals_remy_table(example_data_dir, pattern, monkeypatch):
    path = _file(example_data_dir, pattern)
    br = BasisREMY()
    with_miox = dict(br.runREMY(path))
    monkeypatch.setattr(miox_reader, 'read', lambda *_: None)
    with_remy = dict(br.runREMY(path))
    assert set(with_miox) == set(with_remy)
    off = [k for k in with_remy if not _same(with_miox[k], with_remy[k])]
    assert not off, {k: (with_remy[k], with_miox[k]) for k in off}
    assert with_miox['TE'] and with_miox['Manufacturer']


def test_philips_dicom_fills_the_sheet(example_data_dir):
    path = _file(example_data_dir, 'spec2nii_tests/philips/DICOM/SV_phantom_center/IM-0018-0002-0001.dcm')
    unif = BasisREMY().runREMY(path)
    assert float(unif['TE']) == pytest.approx(40.0)
    assert 'Philips' in unif['Manufacturer'] and unif['Nucleus'] == '1H'
    assert int(unif['NumberOfDatapoints']) == 1024


def test_without_pymiox(monkeypatch, example_data_dir):
    path = _file(example_data_dir, 'BasisREMY_testDatasets/Dataset_28_Siemens_PRESS/*.dat')
    monkeypatch.setitem(sys.modules, 'miox', None)
    assert miox_reader.read(path, '.dat') is None
    assert miox_reader.read(path, 'method') is None

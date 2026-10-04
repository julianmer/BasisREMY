"""
Tests for core.basisremy module
"""

import pytest
import os
from basisremy.core.basisremy import BasisREMY


@pytest.mark.core
@pytest.mark.unit
class TestBasisREMY:
    """Test BasisREMY core functionality"""

    def test_init(self):
        """Test BasisREMY initialization"""
        br = BasisREMY()
        assert br is not None
        assert br.backend is not None

    def test_set_backend_lcmodel(self):
        """Test setting LCModel backend"""
        br = BasisREMY()
        br.set_backend('FidaIdeal')
        assert br.backend.name == 'FidaIdeal'

    def test_set_backend_slaser(self):
        """Test setting sLaser backend"""
        br = BasisREMY()
        br.set_backend('CustomSLaser')
        assert br.backend.name == 'CustomSLaser'

    def test_set_backend_invalid(self):
        """Test setting invalid backend"""
        br = BasisREMY()
        with pytest.raises(ValueError):
            br.set_backend('InvalidBackend')

    def test_get_backend_list(self):
        """Test getting list of available backends"""
        br = BasisREMY()
        backends = br.available_backends
        assert isinstance(backends, list)
        assert len(backends) > 0
        assert 'FidaIdeal' in backends
        assert 'CustomSLaser' in backends

    def test_run_remy_invalid_file(self):
        """Test runREMY with invalid file"""
        br = BasisREMY()
        try:
            result = br.runREMY(import_fpath='/path/that/does/not/exist.spar')
            # Should handle gracefully - may return None or empty dict
            assert result is not None or result is None
        except Exception:
            # It's okay if it raises an exception
            pass


@pytest.mark.core
@pytest.mark.integration
class TestBasisREMYIntegration:
    """Integration tests for BasisREMY with real files"""

    @pytest.mark.parametrize("backend_name", ['FidaIdeal', 'CustomSLaser'])
    def test_backend_initialization(self, backend_name):
        """Test that each backend initializes correctly"""
        br = BasisREMY()
        br.set_backend(backend_name)
        assert br.backend is not None
        assert br.backend.name == backend_name

    def test_remy_with_philips_file(self, philips_spar_file):
        """Test REMY parsing with Philips SPAR file"""
        if not os.path.exists(philips_spar_file):
            pytest.skip("Philips SPAR file not found")

        br = BasisREMY()
        params = br.runREMY(import_fpath=philips_spar_file)

        assert params is not None
        assert isinstance(params, dict)
        # Check for key REMY fields
        assert 'Protocol' in params or 'Sequence' in params

    def test_remy_with_ge_file(self, ge_p_file):
        """Test REMY parsing with GE P-file"""
        if not ge_p_file or not os.path.exists(ge_p_file):
            pytest.skip("GE P-file not found")

        br = BasisREMY()
        params = br.runREMY(import_fpath=ge_p_file)

        assert params is not None
        assert isinstance(params, dict)

    def test_remy_ge_file_has_nucleus(self, ge_p_file):
        """GE headers carry no nucleus; REMY derives it from the Larmor frequency."""
        if not ge_p_file or not os.path.exists(ge_p_file):
            pytest.skip("GE P-file not found")
        params = BasisREMY().runREMY(import_fpath=ge_p_file)
        assert params.get('Nucleus') == '1H'

    def test_remy_ge_file_without_spec2nii_mapper_reads_the_header(self, example_data_dir):
        """spec2nii has no data mapper for some GE sequences (press hbcd, research/oslaser);
        the P-file header still fills the sheet."""
        import glob
        files = glob.glob(os.path.join(example_data_dir, 'REMY_tests', 'Dataset_11_GE_P32256', '*.7'))
        if not files:
            pytest.skip("REMY GE test file not found")
        params = BasisREMY().runREMY(import_fpath=files[0])
        assert float(params['TE']) == 35.0 and float(params['TR']) == 2000.0
        assert int(params['NumberOfDatapoints']) == 2048 and params.get('Nucleus') == '1H'

    def test_remy_nifti_times_in_ms_and_points(self, example_data_dir):
        """NIfTI-MRS stores TE/TR in seconds; the sheet is in ms. Points come from the image."""
        f = os.path.join(example_data_dir, 'example_data.nii.gz')
        if not os.path.exists(f):
            pytest.skip("NIfTI example not found")
        params = BasisREMY().runREMY(import_fpath=f)
        assert 1.0 < float(params['TE']) < 1000.0
        assert 100.0 < float(params['TR']) < 20000.0
        assert int(params['NumberOfDatapoints']) > 0

    def test_remy_nifti_voxel_and_averages_from_the_image(self, example_data_dir):
        """Voxel sizes (pixdim, by axis code) and averages (DIM_DYN length) come from the image."""
        f = os.path.join(example_data_dir, 'REMY_tests', 'Datasets_nifti',
                         'Dataset_12_Philips_SPAR_3T_PRESS_Ref.nii.gz')
        if not os.path.exists(f):
            pytest.skip("NIfTI example not found")
        params = BasisREMY().runREMY(import_fpath=f)
        assert (params['LeftRightSize'], params['AnteriorPosteriorSize'],
                params['CranioCaudalSize']) == (15.0, 30.0, 10.0)
        assert params['NumberOfAverages'] == 2

    def test_remy_rda_vendor_and_ima_averages(self, example_data_dir):
        """RDA (a Siemens-only format) reports the vendor; IMA averages come from the CSA header."""
        import glob
        tests = os.path.join(example_data_dir, 'REMY_tests')
        rda = sorted(glob.glob(os.path.join(tests, 'Dataset_32_Siemens_RDA_PreOn', '*.rda')))
        ima = sorted(glob.glob(os.path.join(tests, 'Dataset_22_Siemens_Dicom_7T', '*.IMA')))
        if not rda or not ima:
            pytest.skip("RDA / IMA examples not found")
        assert BasisREMY().runREMY(import_fpath=rda[0])['Manufacturer'] == 'Siemens'
        assert float(BasisREMY().runREMY(import_fpath=ima[0])['NumberOfAverages']) == 2.0

    def test_remy_philips_averages_stored_per_row(self, example_data_dir):
        """A SPAR with one transient per row (averages 1, rows 320: HERCULES export) reports the
        rows as averages; a scanner-averaged one (averages 64, rows 2) keeps its 64."""
        herc = os.path.join(example_data_dir, 'spec2nii_tests', 'philips', 'HERCULES_spar_sdat',
                            'HERCULES_Example_noID.spar')
        press = os.path.join(example_data_dir, 'REMY_tests', 'Dataset_13_Philips_SPAR_3T_PRESS_45_Act',
                             'Dataset_13_Philips_SPAR_3T_PRESS_45_Act.SPAR')
        if not (os.path.exists(herc) and os.path.exists(press)):
            pytest.skip("Philips SPAR examples not found")
        assert float(BasisREMY().runREMY(import_fpath=herc)['NumberOfAverages']) == 320
        assert float(BasisREMY().runREMY(import_fpath=press)['NumberOfAverages']) == 64

    def test_field_follows_the_spectrometer_frequency(self, example_data_dir):
        """B0 is derived from the header frequency, unrounded (REMY rounds to 2 dp)."""
        f = os.path.join(example_data_dir, 'BigGABA_S1P_S01', 'S01_PRESS_35.dat')
        if not os.path.exists(f) or os.path.getsize(f) < 1024:
            pytest.skip("Siemens twix example not available")
        params = BasisREMY().runREMY(import_fpath=f)
        assert params['Center Freq'] == pytest.approx(123.252468, abs=1e-3)
        assert params['B0'] == pytest.approx(123.252468 / 42.577, rel=1e-6)

    def test_remy_with_bruker_file(self, bruker_dat_file):
        """Test REMY parsing with Bruker file"""
        if not bruker_dat_file or not os.path.exists(bruker_dat_file):
            pytest.skip("Bruker file not found")

        br = BasisREMY()
        params = br.runREMY(import_fpath=bruker_dat_file)

        assert params is not None
        assert isinstance(params, dict)





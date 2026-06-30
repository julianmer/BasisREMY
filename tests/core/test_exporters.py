####################################################################################################
#                                       test_exporters.py                                          #
####################################################################################################
#                                                                                                  #
# Purpose: Unit tests for core/exporters.py — the unified basis-set writer.                        #
#          Each export format gets a smoke test that verifies:                                     #
#            • the file/directory is actually created,                                             #
#            • format-specific structural markers are present,                                     #
#            • the reproducibility sidecar JSON is written next to it,                             #
#            • MRSCloud-style params (no Bfield / Center Freq, only Field Strength)                #
#              produce a sensible header (γ·B0 fallback path).                                     #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import json
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from basisremy.core.exporters import (
    export,
    SUPPORTED_FORMATS,
    _make_header,
    _b0_from_params,
)


# ----------------------------------------------------------------- fixtures
def _synthetic_basis(npts: int = 64) -> dict[str, np.ndarray]:
    """Two metabolites: NAA at 2 ppm, Cr at 3 ppm (using a fake bandwidth/cf)."""
    t = np.arange(npts) / 2000.0  # dwell from BW=2000 Hz
    # arbitrary distinct decaying complex sinusoids
    naa = (np.exp(1j * 2 * np.pi * 50 * t) * np.exp(-t * 5)).astype(np.complex128)
    cr  = (np.exp(1j * 2 * np.pi * 100 * t) * np.exp(-t * 5)).astype(np.complex128)
    return {'NAA': naa, 'Cr': cr}


_PARAMS_MRSCLOUD = {
    # MRSCloud-style: no Bfield, no Center Freq — only Field Strength.
    'Field Strength': '3T',
    'Bandwidth':      2000,
    'Samples':        64,
    'TE':             35,
    'Sequence':       'UnEdited',
    'Nucleus':        '1H',
}

_PARAMS_LEGACY = {
    # Legacy backends still pass Bfield + Center Freq explicitly.
    'Bfield':       3.0,
    'Center Freq':  127.7,
    'Bandwidth':    2000,
    'Samples':      64,
    'TE':           35,
    'Sequence':     'PRESS',
    'Nucleus':      '1H',
    'Linewidth':    1.0,
}


# ----------------------------------------------------------------- header
class TestMakeHeader:

    def test_field_strength_only(self):
        """MRSCloud-style: only `Field Strength` → cf derived via γ·B0."""
        hdr = _make_header(
            {'Field Strength': '7T', 'Bandwidth': 4000, 'Samples': 64},
            _synthetic_basis(),
        )
        assert hdr['centralFrequency'] == pytest.approx(42.577 * 7.0, rel=1e-3)
        assert hdr['bandwidth'] == 4000.0
        assert hdr['points'] == 64
        assert hdr['nucleus'] == '1H'

    def test_bfield_only(self):
        hdr = _make_header({'Bfield': 1.5, 'Bandwidth': 2000}, _synthetic_basis())
        assert hdr['centralFrequency'] == pytest.approx(42.577 * 1.5, rel=1e-3)

    def test_explicit_center_freq_mhz(self):
        hdr = _make_header({'Center Freq': 297.2, 'Bandwidth': 4000}, _synthetic_basis())
        assert hdr['centralFrequency'] == pytest.approx(297.2)

    def test_explicit_center_freq_hz(self):
        """Heuristic: a value > 1000 is interpreted as Hz."""
        hdr = _make_header({'Center Freq': 297.2e6, 'Bandwidth': 4000}, _synthetic_basis())
        assert hdr['centralFrequency'] == pytest.approx(297.2)

    def test_defaults_when_params_empty(self):
        """Empty params shouldn't blow up — fallback to 3T defaults."""
        hdr = _make_header({}, _synthetic_basis())
        assert hdr['centralFrequency'] == pytest.approx(42.577 * 3.0, rel=1e-3)
        assert hdr['bandwidth'] == 2000.0   # default
        assert hdr['points'] == 64           # from basis fallback

    def test_empty_basis_raises(self):
        with pytest.raises(ValueError):
            _make_header({}, {})

    def test_garbage_te(self):
        """Non-numeric TE shouldn't crash the header."""
        hdr = _make_header({'TE': 'oops', 'Bandwidth': 2000}, _synthetic_basis())
        assert hdr['echotime'] is None


# ----------------------------------------------------------------- B0 helper
class TestB0FromParams:

    def test_bfield_wins(self):
        assert _b0_from_params({'Bfield': 7.0, 'Field Strength': '3T'}) == 7.0

    def test_field_strength_string(self):
        assert _b0_from_params({'Field Strength': '1.5T'}) == 1.5

    def test_field_strength_no_t_suffix(self):
        assert _b0_from_params({'Field Strength': '7'}) == 7.0

    def test_center_freq_back_derives(self):
        b0 = _b0_from_params({'Center Freq': 297.2})
        assert b0 == pytest.approx(7.0, rel=1e-2)

    def test_default_3t(self):
        assert _b0_from_params({}) == 3.0

    def test_ignores_placeholders(self):
        assert _b0_from_params({'Bfield': 'Select option',
                                'Field Strength': 'missing input'}) == 3.0


# ----------------------------------------------------------------- formats
@pytest.fixture(params=[_PARAMS_MRSCLOUD, _PARAMS_LEGACY],
                ids=['mrscloud_style', 'legacy_style'])
def params(request):
    return request.param


class TestExportFormats:

    def test_lcmodel_basis(self, tmp_path, params):
        out = export(_synthetic_basis(), str(tmp_path / 'basis.basis'),
                     'lcmodel_basis', params)
        assert os.path.exists(out)
        text = open(out).read()
        assert ' $SEQPAR' in text
        assert ' $BASIS1' in text
        # kbsct writes a $NMUSED fitting-defaults block + a $BASIS block per metab
        assert ' $NMUSED' in text
        assert text.count(' $BASIS\n') == 2
        assert "ID='NAA'" in text
        assert "ID='Cr'" in text
        # sidecar
        sidecar = tmp_path / 'basis_sidecar.json'
        assert sidecar.exists()
        sc = json.loads(sidecar.read_text())
        assert sc['export_format'] == 'lcmodel_basis'
        assert set(sc['metabolites']) == {'NAA', 'Cr'}

    def test_lcmodel_raw_folder(self, tmp_path, params):
        out_dir = tmp_path / 'raw'
        out = export(_synthetic_basis(), str(out_dir), 'lcmodel_raw', params)
        assert os.path.isdir(out)
        for name in ('NAA', 'Cr'):
            fp = out_dir / f'{name}.RAW'
            assert fp.exists(), f'{name}.RAW missing'
            assert ' $NMID' in fp.read_text()
        assert (out_dir / 'basis_sidecar.json').exists()

    def test_jmrui_txt_folder(self, tmp_path, params):
        out_dir = tmp_path / 'jmrui'
        export(_synthetic_basis(), str(out_dir), 'jmrui_txt', params)
        for name in ('NAA', 'Cr'):
            fp = out_dir / f'{name}.txt'
            assert fp.exists()
            text = fp.read_text()
            assert 'PointsInDataset: 64' in text
            assert 'SamplingInterval' in text
            # B0 must be > 0 even for MRSCloud-style params (was the bug!)
            for line in text.splitlines():
                if line.startswith('MagneticField:'):
                    val = float(line.split(':', 1)[1].strip())
                    assert val > 0, "MagneticField must be derived from Field Strength fallback"

    def test_fsl_json_folder(self, tmp_path, params):
        out_dir = tmp_path / 'fsl'
        export(_synthetic_basis(), str(out_dir), 'fsl_json', params)
        # kbsct FSL-MRS writer produces one <metab>.json per metabolite, each
        # holding a 'basis' block (real/imag FID) and a 'meta' block.
        json_files = [f for f in out_dir.glob('*.json')
                      if f.name != 'basis_sidecar.json']
        assert len(json_files) >= 2, \
            f"Expected per-metab files, got {[p.name for p in out_dir.iterdir()]}"
        for f in json_files:
            payload = json.loads(f.read_text())
            assert 'basis' in payload
            b = payload['basis']
            assert 'basis_re' in b and 'basis_im' in b
            assert b['basis_name'] in {'NAA', 'Cr'}

    def test_osprey_mat(self, tmp_path, params):
        scipy_io = pytest.importorskip('scipy.io')
        out = export(_synthetic_basis(), str(tmp_path / 'basis.mat'),
                     'osprey_mat', params)
        assert os.path.exists(out)
        loaded = scipy_io.loadmat(out, simplify_cells=True)
        assert 'BASIS' in loaded
        b = loaded['BASIS']
        # struct keys
        for required in ('fids', 'specs', 'name', 'ppm', 'Bo', 'centerFreq', 'n'):
            assert required in b, f"Osprey BASIS struct missing key {required!r}"
        # fids shape: (npts, nmetabs). The kbsct Osprey writer appends a
        # synthetic H2O peak, so there is at least one column per input metab.
        assert b['fids'].shape[0] == 64
        assert b['fids'].shape[1] >= 2
        # names may be char-array padded on the MATLAB round-trip -> strip.
        names = [str(x).strip() for x in np.atleast_1d(b['name'])]
        assert 'NAA' in names and 'Cr' in names
        # Bo should be sensible (>0) regardless of which params style we used
        assert float(b['Bo']) > 0

    def test_unknown_format_raises(self, tmp_path):
        with pytest.raises(ValueError):
            export(_synthetic_basis(), str(tmp_path / 'x'), 'not_a_format', {})


# ----------------------------------------------------------------- sidecar
class TestSidecar:

    def test_sidecar_contains_metadata(self, tmp_path):
        export(_synthetic_basis(), str(tmp_path / 'basis.basis'),
               'lcmodel_basis', _PARAMS_LEGACY,
               extra_metadata={'note': 'unit test'})
        sc = json.loads((tmp_path / 'basis_sidecar.json').read_text())
        assert sc['tool'] == 'BasisREMY'
        assert 'tool_version' in sc
        assert 'timestamp_utc' in sc
        assert sc['n_points'] == 64
        assert sc['parameters']['Sequence'] == 'PRESS'
        assert sc['extra']['note'] == 'unit test'

    def test_supported_formats_constant(self):
        assert set(SUPPORTED_FORMATS) == {
            'lcmodel_basis', 'lcmodel_raw', 'jmrui_txt', 'fsl_json', 'osprey_mat',
        }


if __name__ == '__main__':
    pytest.main([__file__, '-v'])


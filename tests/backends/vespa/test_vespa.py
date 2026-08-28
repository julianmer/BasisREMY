####################################################################################################
#                                        test_vespa.py                                             #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 27/08/26                                                                                #
#                                                                                                  #
# Purpose: Test the Vespa (PyGAMMA) backend scaffold: registry integration, parameter schema,      #
#          REMY parsing, sequence mapping, and the availability guidance at simulation time.       #
#                                                                                                  #
####################################################################################################


#*************#
#   imports   #
#*************#
import sys, os
import pytest
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

from basisremy.core.basisremy import BasisREMY
from basisremy.backends.vespa_backend import VespaBackend


class TestVespaBackend:

    def test_registered_with_own_category(self):
        br = BasisREMY()
        assert 'Vespa' in br.backends
        assert 'Vespa' in br.categories
        assert br.categories['Vespa'] == ['Vespa']

    def test_switch_to_vespa(self):
        br = BasisREMY()
        br.set_backend('Vespa')
        assert br.backend.name == 'Vespa'
        assert br.get_current_category() == 'Vespa'

    def test_no_scan_physics_defaults(self):
        b = VespaBackend()
        # Scan physics must come from REMY or the user — never fake defaults
        for key in ('Sequence', 'Samples', 'Bandwidth', 'Bfield', 'TE',
                    'Center Freq'):
            assert b.mandatory_params[key] is None, key
        assert b.mandatory_params['Nucleus'] == '1H'
        assert b.mandatory_params['Metabolites']

    def test_tm_hidden_unless_steam(self):
        b = VespaBackend()
        b.mandatory_params['Sequence'] = 'PRESS'
        assert 'TM' not in b.get_params_for_mode()
        b.mandatory_params['Sequence'] = 'STEAM'
        assert 'TM' in b.get_params_for_mode()

    def test_parse_remy_mapping(self):
        b = VespaBackend()
        mandatory, optional = b.parseREMY({
            'NumberOfDatapoints': 4096, 'SpectralWidth': 4000,
            'B0': 2.89, 'TE': 30, 'Protocol': 'svs_se_30',
            'Center Freq': 123.25, 'TR': 2000,
        })
        assert mandatory['Samples'] == 4096
        assert mandatory['Bandwidth'] == 4000
        assert mandatory['Bfield'] == pytest.approx(2.89)
        assert mandatory['Sequence'] == 'Spin Echo'
        assert mandatory['Center Freq'] == pytest.approx(123.25)
        assert optional['TR'] == 2000

    def test_center_freq_derived_from_bfield(self):
        b = VespaBackend()
        mandatory, _ = b.parseREMY({'B0': 3.0})
        assert mandatory['Center Freq'] == pytest.approx(42.577 * 3.0)

    def test_sequence_synonyms(self):
        b = VespaBackend()
        assert b.map_sequence_in('UnEdited') == 'PRESS'
        assert b.map_sequence_in('steam_te11') is None  # not supported (yet)
        assert b.map_sequence_in('MEGA-PRESS') is None  # not supported (yet)

    def test_unsupported_sequence_raises(self):
        b = VespaBackend()
        params = dict(b.mandatory_params)
        params['Sequence'] = 'STEAM'
        with pytest.raises(ValueError, match='unsupported Sequence'):
            b.run_simulation(params)


class TestPyGammaManager:
    """Worker plumbing, tested with a fake worker (no PyGAMMA needed)."""

    def _fake_worker(self, tmp_path, body):
        w = tmp_path / 'fake_worker.py'
        w.write_text(body)
        return w

    def test_run_worker_roundtrip(self, tmp_path):
        import sys as _sys
        from basisremy.core.pygamma_manager import run_worker
        worker = self._fake_worker(tmp_path, (
            "import json, sys\n"
            "job = json.load(open(sys.argv[1]))\n"
            "out = {'ok': True, 'basis': {m: {'re': [1.0], 'im': [0.0]}\n"
            "       for m in job['metabolites']}}\n"
            "json.dump(out, open(sys.argv[2], 'w'))\n"
        ))
        basis = run_worker({'metabolites': {'NAA': []}},
                           python=_sys.executable, worker=worker)
        assert basis['NAA']['re'] == [1.0]

    def test_run_worker_error_surfaces(self, tmp_path):
        import sys as _sys
        from basisremy.core.pygamma_manager import run_worker
        worker = self._fake_worker(tmp_path, (
            "import json, sys\n"
            "json.dump({'ok': False, 'error': 'boom'}, open(sys.argv[2], 'w'))\n"
        ))
        with pytest.raises(RuntimeError, match='boom'):
            run_worker({'metabolites': {}}, python=_sys.executable, worker=worker)

    def test_run_worker_timeout(self, tmp_path):
        import sys as _sys
        from basisremy.core.pygamma_manager import run_worker
        worker = self._fake_worker(tmp_path, "import time\ntime.sleep(30)\n")
        with pytest.raises(RuntimeError, match='timed out'):
            run_worker({'metabolites': {}}, timeout=1.0,
                       python=_sys.executable, worker=worker)


@pytest.mark.backend
@pytest.mark.slow
@pytest.mark.requires_docker
class TestVespaLive:
    """Live PyGAMMA simulation through the real backend (Docker runtime)."""

    def test_press_naa_cr(self):
        import numpy as np
        from basisremy.core import pygamma_manager
        if pygamma_manager.preferred_runtime() == 'env' \
                and not pygamma_manager.is_available():
            pytest.skip("no PyGAMMA runtime available")
        br = BasisREMY()
        br.set_backend('Vespa')
        basis = br.backend.run_simulation({
            'Sequence': 'PRESS', 'Samples': 2048, 'Bandwidth': 2000,
            'Bfield': 3.0, 'TE': 35, 'Nucleus': '1H',
            'Center Freq': 127.732, 'Metabolites': ['NAA', 'Cr'],
        })
        assert not br.backend.last_failures
        ppm = np.linspace(-1000, 1000, 2048) / 127.732 + 4.65
        for metab, expected in (('NAA', 2.01), ('Cr', 3.03)):
            spec = np.abs(np.fft.fftshift(np.fft.fft(basis[metab])))
            peak = ppm[np.argmax(spec)]
            assert abs(peak - expected) < 0.05, \
                f"{metab} peak at {peak:.2f}, expected {expected}"

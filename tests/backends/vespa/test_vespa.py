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
        assert b.map_sequence_in('steam_te11') == 'STEAM'
        assert b.map_sequence_in('MEGA-PRESS') is None  # not supported (yet)

    @pytest.mark.skipif(VespaBackend.pygamma_available(),
                        reason="PyGAMMA installed — guidance path not reachable")
    def test_run_without_pygamma_gives_guidance(self):
        b = VespaBackend()
        with pytest.raises(RuntimeError, match="PyGAMMA"):
            b.run_simulation(dict(b.mandatory_params))

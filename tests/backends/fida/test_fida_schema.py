####################################################################################################
#                                       test_fida_schema.py                                         #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 25/04/26                                                                                #
#                                                                                                  #
# Purpose: Schema-only tests for the FID-A backend family. No Octave required.                     #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

from basisremy.backends.fida_backends import (
    FIDA_BACKENDS,
    FidaBackend,
    FidaIdeal,
    FidaPressShaped,
    FidaSemiLaserShaped,
    FidaSteamShaped,
    FidaMegaPressShaped,
    FidaLaser,
    FidaOnePulse,
)


# ============================================================ taxonomy
@pytest.mark.parametrize("cls", FIDA_BACKENDS)
def test_every_fida_backend_is_categorized(cls):
    inst = cls()
    assert inst.category == 'FID-A'
    assert inst.requires_octave is True
    assert inst.name and inst.display_name


def test_fida_ideal_identity():
    """FidaIdeal is the canonical 'Ideal' wrapper around sim_lcmrawbasis."""
    inst = FidaIdeal()
    assert inst.name == 'FidaIdeal'
    assert inst.category == 'FID-A'
    assert 'Ideal' in inst.display_name
    assert inst._kind == 'ideal'


def test_registry_has_ideal_first():
    assert FIDA_BACKENDS[0] is FidaIdeal


# ============================================================ FidaIdeal
class TestFidaIdeal:
    def test_sequences_dropdown(self):
        b = FidaIdeal()
        assert b.dropdown['Sequence'] == ['Spin Echo', 'PRESS', 'STEAM', 'LASER']

    def test_seq_to_fida_mapping(self):
        assert FidaIdeal._seq_to_fida('Spin Echo') == 'se'
        assert FidaIdeal._seq_to_fida('PRESS')     == 'p'
        assert FidaIdeal._seq_to_fida('STEAM')     == 'st'
        assert FidaIdeal._seq_to_fida('LASER')     == 'l'

    def test_seq_to_fida_invalid_raises(self):
        with pytest.raises(ValueError, match="unrecognised Sequence"):
            FidaIdeal._seq_to_fida('UnEdited')

    def test_parse_protocol(self):
        b = FidaIdeal()
        assert b.parseProtocol('PRESS_TE30')  == 'PRESS'
        assert b.parseProtocol('svs_steam')   == 'STEAM'
        assert b.parseProtocol('UnEdited')    == 'PRESS'  # BigGABA / MRSCloud convention
        assert b.parseProtocol('svs_slaserVOI') is None  # unsupported


# ============================================================ FidaPressShaped
class TestPressShaped:
    def test_required_params(self):
        b = FidaPressShaped()
        for k in ('Samples', 'Bandwidth', 'Bfield', 'Linewidth',
                  'TE', 'RefTp', 'Path to Pulse',
                  'thkX', 'thkY', 'fovX', 'fovY', 'nX', 'nY'):
            assert k in b.mandatory_params, f"missing {k}"
        assert 'Path to Pulse' in b.file_selection
        assert b._kind == 'press_shaped'

    def test_build_args_requires_pulse_path(self):
        b = FidaPressShaped()
        b.mandatory_params.update({
            'Samples': 2048, 'Bandwidth': 4000, 'Bfield': 3.0,
            'TE': 30, 'Path to Pulse': None,
        })
        with pytest.raises(ValueError, match='Path to Pulse'):
            b._build_args(b.mandatory_params, 'NAA')

    def test_build_args_returns_positional_list(self, tmp_path):
        b = FidaPressShaped()
        pulse = tmp_path / 'sample.pta'; pulse.write_text('# fake')
        b.mandatory_params.update({
            'Samples': 2048, 'Bandwidth': 4000, 'Bfield': 3.0,
            'Linewidth': 1.0, 'TE': 30, 'Path to Pulse': str(pulse),
        })
        args = b._build_args(b.mandatory_params, 'NAA')
        # FidaBackend.run_simulation prepends (metab, kind); _build_args
        # returns just the trailing positional args. Length = 16 for press_shaped.
        assert len(args) == 16
        assert args[0] == 2048.0          # n
        assert args[1] == 4000.0          # sw
        assert args[2] == 3.0             # Bfield


# ============================================================ stubs
@pytest.mark.parametrize("cls", [
    FidaSemiLaserShaped, FidaSteamShaped, FidaMegaPressShaped,
    FidaLaser, FidaOnePulse,
])
def test_stub_run_raises_not_implemented(cls):
    b = cls()
    assert b._is_stub is True
    with pytest.raises(NotImplementedError):
        b.run_simulation({'Metabolites': ['NAA']})


# ============================================================ parseREMY
def test_parse_remy_only_returns_known_keys():
    b = FidaLaser()
    mandatory, _ = b.parseREMY({
        'NumberOfDatapoints': 2048, 'SpectralWidth': 4000,
        'B0': 3.0, 'TE': 30,
    })
    assert mandatory['Samples'] == 2048
    # FidaLaser doesn't expose Tau 1/Tau 2 — parseREMY must not inject them
    assert 'Tau 1' not in mandatory


if __name__ == '__main__':
    pytest.main([__file__, '-v'])


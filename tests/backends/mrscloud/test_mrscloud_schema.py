####################################################################################################
#                                  test_mrscloud_schema.py                                         #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 24/04/26                                                                                #
#                                                                                                  #
# Purpose: Tests for MRSCloudBackend.get_params_for_mode() — the schema-aware                      #
#          method that hides irrelevant fields and dynamically adds the                            #
#          vendor-pulse picker. Also exercises required_pulse_files /                              #
#          missing_pulse_files. No Octave needed.                                                  #
#                                                                                                  #
####################################################################################################

from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

from basisremy.backends.mrscloud_backend import MRSCloudBackend


# ============================================================ schema-aware
class TestSchemaAware:

    @pytest.fixture
    def backend(self):
        return MRSCloudBackend()

    # ---- Edit-field visibility ---------------------------------------
    def test_unedited_hides_edit_fields(self, backend):
        backend.mandatory_params['Sequence'] = 'UnEdited'
        params = backend.get_params_for_mode()
        for k in ('Edit Target', 'Edit On', 'Edit Off', 'Edit Tp'):
            assert k not in params, f"{k!r} must be hidden for UnEdited"

    def test_mega_shows_all_edit_fields(self, backend):
        backend.mandatory_params['Sequence'] = 'MEGA'
        # Need a vendor + localization compatible with editing
        backend.mandatory_params['System'] = 'Philips'
        backend.mandatory_params['Localization'] = 'PRESS'
        params = backend.get_params_for_mode()
        for k in ('Edit Target', 'Edit On', 'Edit Off', 'Edit Tp'):
            assert k in params, f"{k!r} must be shown for MEGA"

    @pytest.mark.parametrize("seq", ['HERMES', 'HERCULES'])
    def test_hermes_hercules_show_only_target_and_tp(self, backend, seq):
        backend.mandatory_params['Sequence'] = seq
        backend.mandatory_params['System'] = 'Philips'
        backend.mandatory_params['Localization'] = 'PRESS'
        params = backend.get_params_for_mode()
        assert 'Edit Target' in params
        assert 'Edit Tp'     in params
        # MRSCloud overrides offsets internally → don't show
        assert 'Edit On'     not in params
        assert 'Edit Off'    not in params

    # ---- Localization restriction ------------------------------------
    @pytest.mark.parametrize("seq", ['MEGA', 'HERMES', 'HERCULES'])
    def test_edited_seq_restricts_localization(self, backend, seq):
        backend.mandatory_params['Sequence'] = seq
        backend.mandatory_params['System'] = 'Philips'
        backend.get_params_for_mode()
        assert backend.dropdown['Localization'] == ['PRESS', 'sLASER'], \
            "Edited sequences must hide STEAM (7T only) from Localization choices"

    def test_edited_seq_replaces_steam_carryover(self, backend):
        """Switching from UnEdited+STEAM_7T → MEGA must not leave a stale value."""
        backend.mandatory_params['Sequence'] = 'UnEdited'
        backend.mandatory_params['Localization'] = 'STEAM_7T'
        backend.mandatory_params['Field Strength'] = '7T'
        backend.get_params_for_mode()  # establishes UnEdited state
        backend.mandatory_params['Sequence'] = 'MEGA'
        backend.mandatory_params['System'] = 'Philips'
        params = backend.get_params_for_mode()
        assert params['Localization'] == 'PRESS'
        assert backend.mandatory_params['Localization'] == 'PRESS'

    def test_unedited_restores_full_localization(self, backend):
        backend.mandatory_params['Sequence'] = 'MEGA'
        backend.mandatory_params['Field Strength'] = '7T'
        backend.get_params_for_mode()  # narrow
        backend.mandatory_params['Sequence'] = 'UnEdited'
        backend.get_params_for_mode()  # restore
        assert set(backend.dropdown['Localization']) == \
            {'PRESS', 'sLASER', 'STEAM_7T'}

    # ---- Vendor-pulse-file picker -----------------------------------
    def test_universal_vendor_no_picker(self, backend):
        backend.mandatory_params['Sequence'] = 'UnEdited'
        backend.mandatory_params['Localization'] = 'PRESS'
        backend.mandatory_params['System'] = 'Philips'
        backend.current_mode = 'Universal'
        backend.get_params_for_mode()
        assert backend.file_selection == [], \
            "Universal mode ships its own pulses — no picker should appear"

    def test_philips_vendor_shows_picker_when_pulses_missing(self, backend, tmp_path):
        """When the vendor-confidential pulse files aren't on disk, expose the picker."""
        backend.mandatory_params['Sequence'] = 'UnEdited'
        backend.mandatory_params['Localization'] = 'PRESS'
        backend.mandatory_params['System'] = 'Philips'
        backend.current_mode = 'Non-Universal'
        # Force the missing-file check to look at an empty mrscloud root → all missing
        # Monkey-patch via the classmethod default arg
        original = MRSCloudBackend.missing_pulse_files
        try:
            MRSCloudBackend.missing_pulse_files = classmethod(  # type: ignore
                lambda cls, v, s, l, mrscloud_root=str(tmp_path):
                    [f for f in cls.required_pulse_files(v, s, l)]
            )
            params = backend.get_params_for_mode()
            assert 'Vendor Pulse File' in backend.file_selection
            assert 'Vendor Pulse File' in params
        finally:
            MRSCloudBackend.missing_pulse_files = original  # type: ignore

    # ---- schema_affecting_keys ---------------------------------------
    def test_schema_affecting_keys_declared(self, backend):
        assert 'Sequence' in backend.schema_affecting_keys
        assert 'Localization' in backend.schema_affecting_keys
        assert 'System' in backend.schema_affecting_keys


# ============================================================ pulse map
class TestPulseFiles:

    def test_universal_combos_empty(self):
        for vendor in ('Universal_Philips', 'Universal_Siemens'):
            for seq in ('UnEdited', 'MEGA'):
                assert MRSCloudBackend.required_pulse_files(vendor, seq, 'PRESS') == [], \
                    f"{vendor}/{seq}/PRESS should not require extra pulses"

    def test_philips_press_unedited_lists_pulses(self):
        files = MRSCloudBackend.required_pulse_files('Philips', 'UnEdited', 'PRESS')
        assert any(f.endswith('.pta') for f in files)
        # Philips_spredrex.pta is hard-coded as excitation for ALL vendors
        assert any('Philips_spredrex' in f for f in files)
        assert any('gtst1203_sp' in f for f in files)

    def test_unknown_combo_returns_empty(self):
        assert MRSCloudBackend.required_pulse_files('Bruker', 'MEGA', 'sLASER') == []

    def test_missing_pulse_files_with_bogus_root(self, tmp_path):
        """All required pulses missing → all reported as missing."""
        missing = MRSCloudBackend.missing_pulse_files(
            'Philips', 'UnEdited', 'PRESS', mrscloud_root=str(tmp_path))
        assert len(missing) >= 2

    def test_missing_pulse_files_universal_empty(self, tmp_path):
        """Universal_* combos don't require any pulses → empty even with bogus root."""
        missing = MRSCloudBackend.missing_pulse_files(
            'Universal_Philips', 'UnEdited', 'PRESS', mrscloud_root=str(tmp_path))
        assert missing == []


if __name__ == '__main__':
    pytest.main([__file__, '-v'])


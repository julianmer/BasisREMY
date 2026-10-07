####################################################################################################
#                                     test_sequence_setup.py                                       #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 07/10/26                                                                                #
#                                                                                                  #
# Purpose: The engine / sequence / pulse routes land on the backends exactly as before (same       #
#          backend, mode and sheet values), ideal is preferred, impossible choices say why, and    #
#          the recommended values follow TE and never overwrite what the user or the file set.     #
#                                                                                                  #
####################################################################################################

import pytest

from basisremy.core import sequence_setup as ss
from basisremy.core.basisremy import BasisREMY


@pytest.fixture(scope='module')
def br():
    return BasisREMY()


def test_every_route_round_trips(br):
    for cat, seqs in ss.ROUTES.items():
        for seq, routes in seqs.items():
            for route in routes:
                want = {k: v for k, v in route.pulses.items() if v != ss.OWN}
                assert ss.apply(br, cat, seq, want) == route
                assert ss.current(br) == (seq, route), (cat, seq, route.backend)


def test_routes_select_the_backends_as_before(br):
    ss.apply(br, 'FID-A', 'PRESS')
    assert br.backend.name == 'FidaIdeal' and br.backend.mandatory_params['Sequence'] == 'PRESS'
    ss.apply(br, 'FID-A', 'PRESS', {'ref': ss.SHAPED})
    assert br.backend.name == 'FidaPressShaped'
    ss.apply(br, 'Vespa', 'PRESS', {'ref': ss.SHAPED})
    assert br.backend.mandatory_params['Sequence'] == 'PRESS shaped'
    ss.apply(br, 'FID-A', 'MEGA-PRESS')
    assert br.backend.current_mode == 'Edit-only shaped (ideal refoc)'
    ss.apply(br, 'FID-A', 'MEGA-PRESS', {'ref': ss.SHAPED})
    assert br.backend.current_mode == 'Full shaped (refoc + edit)'
    ss.apply(br, 'MRSCloud', 'HERMES (sLASER)')
    p = br.backend.mandatory_params
    assert (p['Sequence'], p['Localization']) == ('HERMES', 'sLASER')
    # an ideal editing pulse is never offered
    assert all(r.backend != 'FidaMegaPressIdeal' for s in ss.ROUTES.values() for rs in s.values()
               for r in rs)


def test_ideal_preferred_and_closest_kept():
    assert ss.choose('FID-A', 'STEAM').backend == 'FidaIdeal'
    assert ss.choose('FID-A', 'sLASER').backend == 'FidaSemiLaserShaped'     # no ideal sLASER
    # a wish the engine cannot meet falls back to what it has
    assert ss.choose('FSL-MRS', 'PRESS', {'ref': ss.SHAPED}).pulses['ref'] == ss.IDEAL
    assert ss.choose('Vespa', 'sLASER') is None


def test_reasons():
    assert 'never ideal' in ss.pulse_note('FID-A', 'MEGA-PRESS', 'edit', ss.IDEAL)
    note = ss.pulse_note('FSL-MRS', 'PRESS', 'ref', ss.SHAPED)
    assert 'FID-A' in note and 'Vespa' in note
    assert ss.pulse_note('FID-A', 'PRESS', 'ref', ss.SHAPED) is None
    assert 'own' in ss.pulse_note('MRSCloud', 'PRESS', 'ref', ss.IDEAL)
    assert 'FSL-MRS' in ss.why_not('Vespa', 'sLASER')
    lv = ss.levels('PRESS')
    assert 'FSL-MRS' in lv['ideal'] and 'Vespa' in lv['shaped'] and lv['own'] == ['MRSCloud']


def test_recommended_values(br):
    ss.apply(br, 'FID-A', 'PRESS', {'ref': ss.SHAPED})
    b = br.backend
    b.mandatory_params['TE'] = 30
    rec = ss.apply_recommended(b, 'PRESS', keep=set())
    assert b.mandatory_params['Tau 1'] == 15 and b.mandatory_params['Tau 2'] == 15
    assert 'Tau 1' in rec
    # the user's value stays
    b.mandatory_params['Tau 1'] = 11
    ss.apply_recommended(b, 'PRESS', keep={'Tau 1'})
    assert b.mandatory_params['Tau 1'] == 11

    ss.apply(br, 'Spant', 'sLASER')
    br.backend.mandatory_params['TE'] = 32
    ss.apply_recommended(br.backend, 'sLASER', keep=set())
    p = br.backend.mandatory_params
    assert (p['sLASER TE1'], p['sLASER TE2'], p['sLASER TE3']) == (8, 16, 8)

    ss.apply(br, 'FSL-MRS', 'MEGA-PRESS')
    rec = ss.recommended(br.backend, 'MEGA-PRESS')
    assert rec['Edit Tp'][0] == 15 and rec['Edit On'][0] == 1.9 and 'Saleh' in rec['Edit Tp'][1]
    ss.apply(br, 'FSL-MRS', 'HERMES')
    assert ss.recommended(br.backend, 'HERMES')['Edit Tp'][0] == 20
    ss.apply(br, 'FSL-MRS', 'MEGA-sLASER')
    br.backend.mandatory_params['TE'] = 80
    assert ss.recommended(br.backend, 'MEGA-sLASER')['Edit Tp'][0] == 15     # fits
    br.backend.mandatory_params['TE'] = 40      # TE/4 gaps too short for 15 ms
    assert ss.recommended(br.backend, 'MEGA-sLASER')['Edit Tp'][0] == pytest.approx(9.98)
    ss.apply(br, 'MRSCloud', 'MEGA-PRESS')
    br.backend.get_params_for_mode()
    assert 'Edit Tp' not in ss.recommended(br.backend, 'MEGA-PRESS')   # MRSCloud's own pulses
    ss.apply(br, 'Spant', 'MEGA-PRESS')
    assert ss.recommended(br.backend, 'MEGA-PRESS')['Edit Bandwidth (Hz)'][0] == pytest.approx(151.6)


def test_value_state():
    rec = {'Tau 1': 'x'}
    assert ss.value_state(None, 'TE', None, set(), {}, rec) == 'missing'
    assert ss.value_state(None, 'TE', 30, set(), {'TE': 30}, rec) == 'file'
    assert ss.value_state(None, 'TE', 35, {'TE'}, {'TE': 30}, rec) == 'user'
    assert ss.value_state(None, 'Tau 1', 15, set(), {}, rec) == 'recommended'
    assert ss.value_state(None, 'Linewidth', 1.0, set(), {}, rec) == 'default'


def test_apply_remy_records_the_file_values(br):
    br.set_backend('Vespa')
    br.apply_remy({'TE': 30, 'B0': 2.89, 'NumberOfDatapoints': 2048, 'SpectralWidth': 4000,
                   'Protocol': 'PRESS'})
    got = br.from_file['Vespa']
    assert got['TE'] == 30 and got['Sequence'] == 'PRESS' and 'TM' not in got
    br._last_mrsinmrs = None

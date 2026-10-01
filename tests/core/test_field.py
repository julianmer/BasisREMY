"""core/field.py: frequency primary, field derived; nucleus-aware."""

import pytest

from basisremy.core.field import reconcile_field, gamma_mhz_per_t


def test_frequency_wins_over_rounded_b0():
    out = reconcile_field({'Center Freq': 123.25, 'B0': 2.89, 'Nucleus': '1H'})
    assert out['Center Freq'] == 123.25
    assert out['B0'] == pytest.approx(123.25 / 42.577)


def test_b0_gives_frequency_when_no_frequency():
    out = reconcile_field({'B0': 3.0})
    assert out == {'Center Freq': pytest.approx(3.0 * 42.577), 'B0': 3.0}


def test_nucleus_changes_gamma():
    out = reconcile_field({'Center Freq': 51.7, 'Nucleus': '31P'})
    assert out['B0'] == pytest.approx(51.7 / 17.235)
    assert gamma_mhz_per_t('unknown') == 42.577


@pytest.mark.parametrize('bad', [{}, {'Center Freq': ''}, {'B0': 'missing input'}])
def test_blank_stays_blank(bad):
    assert reconcile_field(bad) == {}

"""Cross-engine agreement on ideal PRESS: FID-A is the reference, FSL-MRS
(pure Python) and Spinach (same Octave runtime as FID-A) must agree."""

import pytest

from basisremy.core.basisremy import BasisREMY
from tests.validation.cross_engine import PANEL, run_engine, agreement_table

pytestmark = pytest.mark.requires_octave_runtime

B0, BW, N, TE = 3.0, 2000.0, 2048, 35.0


@pytest.fixture(scope='module')
def results():
    br = BasisREMY('FidaIdeal')
    out = {}
    for engine in ('FidaIdeal', 'FSL-MRS', 'Spinach'):
        out[engine], _ = run_engine(br, engine, 'PRESS', TE, PANEL, B0, N, BW)
    return out


@pytest.mark.parametrize('engine', ['FSL-MRS', 'Spinach'])
@pytest.mark.parametrize('metab', PANEL)
def test_ideal_press_agrees_with_fida(results, engine, metab):
    table = agreement_table(results, 'FidaIdeal', BW, B0)
    c = table[engine][metab]
    assert c['r'] > 0.99, c
    assert abs(c['peak_a'] - c['peak_b']) < 0.02, c
